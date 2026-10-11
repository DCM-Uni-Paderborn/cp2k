!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Native contracted Gaussian fields and their consistent Cartesian/Bloch responses.
!>        Uses the stored sphi coefficients and native spherical AO order without renormalization.
! **************************************************************************************************
MODULE gaussian_basis_fields
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: twopi
   USE orbital_pointers,                ONLY: indco,&
                                              ncoset

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: gaussian_basis_fields_evaluate, gaussian_basis_fields_periodic, gaussian_basis_workspace_type

   TYPE gaussian_basis_workspace_type
      PRIVATE
      REAL(dp), ALLOCATABLE :: displacement(:, :), field(:, :), ddisplacement(:, :), &
                               gradient(:, :, :), hessian(:, :, :, :)
      REAL(dp), ALLOCATABLE :: coordinate_power(:, :, :), cartesian_value(:, :), &
                               cartesian_gradient(:, :, :), cartesian_hessian(:, :, :, :)
   END TYPE gaussian_basis_workspace_type

CONTAINS

! **************************************************************************************************
!> \brief Evaluate a native contracted basis at explicit point-minus-center displacements.
!>        No minimum-image operation is applied; an explicit cutoff can skip excluded points.
!> \param basis Initialized native basis, including its original normalization in sphi.
!> \param displacements Cartesian displacements, indexed (axis,point).
!> \param value New values, indexed (point,AO).
!> \param gradient Optional Cartesian derivatives, indexed (point,AO,axis).
!> \param hessian Optional Cartesian Hessians, indexed (point,AO,axis,axis).
!> \param cutoff_squared Optional squared cutoff radius; fields beyond it are zero.
!> \param workspace Optional caller-owned scratch arrays, with independent storage per OpenMP thread.
! **************************************************************************************************
   SUBROUTINE gaussian_basis_fields_evaluate(basis, displacements, value, gradient, hessian, cutoff_squared, workspace)
      TYPE(gto_basis_set_type), INTENT(IN)               :: basis
      REAL(dp), INTENT(IN)                               :: displacements(:, :)
      REAL(dp), INTENT(OUT)                              :: value(:, :)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: gradient(:, :, :), hessian(:, :, :, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: cutoff_squared
      TYPE(gaussian_basis_workspace_type), &
         INTENT(INOUT), OPTIONAL, TARGET                 :: workspace

      INTEGER                                            :: axis, axis2, first, ico, ipgf, ipoint, &
                                                            irow, iset, isgf, ishell, l, last, &
                                                            maxl, ncart, powers(3), reduced(3), &
                                                            thread, threads
      REAL(dp) :: alpha, exponential, polynomial, polynomial_gradient(3), &
         polynomial_hessian(3, 3), primitive_gradient(3), primitive_hessian(3, 3), r(3), &
         r_squared, weight
      REAL(dp), POINTER                                  :: cartesian_gradient(:, :), &
                                                            cartesian_hessian(:, :, :), &
                                                            cartesian_value(:), &
                                                            coordinate_power(:, :)
      TYPE(gaussian_basis_workspace_type), POINTER       :: scratch
      TYPE(gaussian_basis_workspace_type), TARGET        :: local_workspace

      CPASSERT(basis%nset > 0 .AND. basis%nsgf > 0)
      CPASSERT(ASSOCIATED(basis%sphi) .AND. ASSOCIATED(basis%zet))
      CPASSERT(UBOUND(ncoset, 1) >= MAXVAL(basis%lmax))
      CPASSERT(SIZE(displacements, 1) == 3)
      CPASSERT(ALL(SHAPE(value) == [SIZE(displacements, 2), basis%nsgf]))
      IF (PRESENT(gradient)) THEN
         CPASSERT(ALL(SHAPE(gradient) == [SIZE(value, 1), SIZE(value, 2), 3]))
         gradient = 0.0_dp
      END IF
      IF (PRESENT(hessian)) THEN
         CPASSERT(ALL(SHAPE(hessian) == [SIZE(value, 1), SIZE(value, 2), 3, 3]))
         hessian = 0.0_dp
      END IF
      value = 0.0_dp
      maxl = MAXVAL(basis%lmax)
      threads = 1
!$    threads = omp_get_max_threads()
      scratch => local_workspace
      IF (PRESENT(workspace)) scratch => workspace
      CALL prepare_polynomial_workspace(scratch, maxl, threads, PRESENT(gradient) .OR. PRESENT(hessian), PRESENT(hessian))
      !$OMP PARALLEL DEFAULT(NONE) &
      !$OMP SHARED(basis, displacements, value, gradient, hessian, cutoff_squared, ncoset, indco, scratch) &
      !$OMP PRIVATE(ipoint, r, iset, ncart, ishell, l, first, last, ipgf, alpha, exponential, &
      !$OMP         isgf, ico, irow, weight, powers, reduced, polynomial, polynomial_gradient, &
      !$OMP         polynomial_hessian, primitive_gradient, primitive_hessian, cartesian_value, &
      !$OMP         cartesian_gradient, cartesian_hessian, coordinate_power, r_squared, axis, axis2, thread)
      thread = 1
!$    thread = omp_get_thread_num() + 1
      coordinate_power(1:, 0:) => scratch%coordinate_power(:, :, thread)
      cartesian_value => scratch%cartesian_value(:, thread)
      NULLIFY (cartesian_gradient, cartesian_hessian)
      IF (PRESENT(gradient) .OR. PRESENT(hessian)) cartesian_gradient => scratch%cartesian_gradient(:, :, thread)
      IF (PRESENT(hessian)) cartesian_hessian => scratch%cartesian_hessian(:, :, :, thread)
      !$OMP DO SCHEDULE(STATIC)
      DO ipoint = 1, SIZE(displacements, 2)
         r = displacements(:, ipoint)
         r_squared = SUM(r**2)
         IF (PRESENT(cutoff_squared)) THEN
            IF (r_squared > cutoff_squared) CYCLE
         END IF
         coordinate_power(:, 0) = 1.0_dp
         DO l = 1, UBOUND(coordinate_power, 2)
            coordinate_power(:, l) = coordinate_power(:, l - 1)*r
         END DO
         ! Cartesian monomials and their derivatives do not depend on the primitive or contraction.
         DO ico = 1, SIZE(cartesian_value)
            powers = indco(:, ico)
            cartesian_value(ico) = coordinate_power(1, powers(1))*coordinate_power(2, powers(2))* &
                                   coordinate_power(3, powers(3))
            IF (.NOT. PRESENT(gradient) .AND. .NOT. PRESENT(hessian)) CYCLE
            cartesian_gradient(:, ico) = 0.0_dp
            IF (PRESENT(hessian)) cartesian_hessian(:, :, ico) = 0.0_dp
            DO axis = 1, 3
               IF (powers(axis) == 0) CYCLE
               reduced = powers
               reduced(axis) = reduced(axis) - 1
               cartesian_gradient(axis, ico) = REAL(powers(axis), dp)*coordinate_power(1, reduced(1))* &
                                               coordinate_power(2, reduced(2))*coordinate_power(3, reduced(3))
               IF (.NOT. PRESENT(hessian)) CYCLE
               DO axis2 = 1, 3
                  IF (reduced(axis2) == 0) CYCLE
                  reduced(axis2) = reduced(axis2) - 1
                  cartesian_hessian(axis, axis2, ico) = REAL(powers(axis)*(reduced(axis2) + 1), dp)* &
                                                        coordinate_power(1, reduced(1))*coordinate_power(2, reduced(2))* &
                                                        coordinate_power(3, reduced(3))
                  reduced(axis2) = reduced(axis2) + 1
               END DO
            END DO
         END DO
         DO iset = 1, basis%nset
            ncart = ncoset(basis%lmax(iset))
            DO ipgf = 1, basis%npgf(iset)
               alpha = basis%zet(ipgf, iset)
               exponential = EXP(-alpha*r_squared)
               DO ishell = 1, basis%nshell(iset)
                  l = basis%l(ishell, iset)
                  first = basis%first_sgf(ishell, iset)
                  last = basis%last_sgf(ishell, iset)
                  ! Primitive Cartesian derivatives are shared by every spherical contraction.
                  DO ico = ncoset(l - 1) + 1, ncoset(l)
                     irow = (ipgf - 1)*ncart + ico
                     polynomial = cartesian_value(ico)
                     IF (PRESENT(gradient) .OR. PRESENT(hessian)) THEN
                        polynomial_gradient = cartesian_gradient(:, ico)
                        IF (PRESENT(hessian)) polynomial_hessian = cartesian_hessian(:, :, ico)
                        DO axis = 1, 3
                           primitive_gradient(axis) = polynomial_gradient(axis) - 2.0_dp*alpha*r(axis)*polynomial
                           IF (.NOT. PRESENT(hessian)) CYCLE
                           DO axis2 = 1, 3
                              primitive_hessian(axis, axis2) = polynomial_hessian(axis, axis2) - 2.0_dp*alpha*( &
                                                        r(axis)*polynomial_gradient(axis2) + r(axis2)*polynomial_gradient(axis)) + &
                                                               4.0_dp*alpha**2*r(axis)*r(axis2)*polynomial
                              IF (axis == axis2) THEN
                                 primitive_hessian(axis, axis2) = primitive_hessian(axis, axis2) - 2.0_dp*alpha*polynomial
                              END IF
                           END DO
                        END DO
                     END IF
                     DO isgf = first, last
                        weight = basis%sphi(irow, isgf)*exponential
                        value(ipoint, isgf) = value(ipoint, isgf) + weight*polynomial
                        IF (.NOT. PRESENT(gradient) .AND. .NOT. PRESENT(hessian)) CYCLE
                        DO axis = 1, 3
                           IF (PRESENT(gradient)) THEN
                              gradient(ipoint, isgf, axis) = gradient(ipoint, isgf, axis) + weight*primitive_gradient(axis)
                           END IF
                           IF (.NOT. PRESENT(hessian)) CYCLE
                           DO axis2 = 1, 3
                              hessian(ipoint, isgf, axis, axis2) = hessian(ipoint, isgf, axis, axis2) + &
                                                                   weight*primitive_hessian(axis, axis2)
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
      !$OMP END DO
      !$OMP END PARALLEL
   END SUBROUTINE gaussian_basis_fields_evaluate

! **************************************************************************************************
!> \brief Sum explicit Gaussian images with physical Bloch gradients and joint geometry response.
!>        Caller owns the finite image set and its convergence; the set is held fixed in a response.
!> \param basis ...
!> \param points Cartesian points, indexed (axis,point).
!> \param center Physical basis center, without implicit folding.
!> \param cell Cartesian lattice vectors, one vector per column.
!> \param images Explicit integer lattice images, indexed (axis,image).
!> \param reduced_k Reduced Bloch vector; phase is exp(i*2*pi*k.image).
!> \param value New Bloch values, indexed (point,AO).
!> \param gradient Optional physical Cartesian Bloch gradients, axis last.
!> \param response Optional joint first response of the values.
!> \param gradient_response Optional response of the physical spatial gradients.
!> \param point_response Optional Cartesian point response, default zero.
!> \param center_response Optional Cartesian center response, default zero.
!> \param cell_response Optional lattice response, default zero.
!> \param reduced_k_response Optional reduced Bloch-vector response, default zero.
!> \param cell_gradient Optional lattice derivative at fixed points, centers and reduced k, axes last.
!> \param workspace Optional caller-owned scratch arrays, refilled at every invocation.
! **************************************************************************************************
   SUBROUTINE gaussian_basis_fields_periodic(basis, points, center, cell, images, reduced_k, value, gradient, &
                                             response, gradient_response, point_response, center_response, &
                                             cell_response, reduced_k_response, cell_gradient, workspace)
      TYPE(gto_basis_set_type), INTENT(IN)               :: basis
      REAL(dp), INTENT(IN)                               :: points(:, :), center(3), cell(3, 3)
      INTEGER, INTENT(IN)                                :: images(:, :)
      REAL(dp), INTENT(IN)                               :: reduced_k(3)
      COMPLEX(dp), INTENT(OUT)                           :: value(:, :)
      COMPLEX(dp), INTENT(OUT), OPTIONAL                 :: gradient(:, :, :), response(:, :), &
                                                            gradient_response(:, :, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: point_response(:, :), &
                                                            center_response(3), &
                                                            cell_response(3, 3), &
                                                            reduced_k_response(3)
      COMPLEX(dp), INTENT(OUT), OPTIONAL                 :: cell_gradient(:, :, :, :)
      TYPE(gaussian_basis_workspace_type), &
         INTENT(INOUT), OPTIONAL, TARGET                 :: workspace

      COMPLEX(dp)                                        :: phase, phase_response
      INTEGER                                            :: ao, axis, axis2, image, npoint, point
      LOGICAL                                            :: need_geometry, need_gradient, &
                                                            need_hessian
      REAL(dp)                                           :: dc(3, 3), dcenter(3), dk(3), dshift(3), &
                                                            lattice(3), shift(3)
      REAL(dp), POINTER                                  :: ddisplacement(:, :), displacement(:, :), &
                                                            field(:, :), g(:, :, :), h(:, :, :, :)
      TYPE(gaussian_basis_workspace_type), POINTER       :: scratch
      TYPE(gaussian_basis_workspace_type), TARGET        :: local_workspace

      npoint = SIZE(points, 2)
      CPASSERT(SIZE(points, 1) == 3 .AND. SIZE(images, 1) == 3)
      CPASSERT(ALL(SHAPE(value) == [npoint, basis%nsgf]))
      need_gradient = PRESENT(gradient) .OR. PRESENT(response) .OR. PRESENT(gradient_response) .OR. PRESENT(cell_gradient)
      need_hessian = PRESENT(gradient_response)
      need_geometry = PRESENT(response) .OR. PRESENT(gradient_response)
      IF (PRESENT(gradient)) THEN
         CPASSERT(ALL(SHAPE(gradient) == [npoint, basis%nsgf, 3]))
         gradient = CMPLX(0.0_dp, 0.0_dp, dp)
      END IF
      IF (PRESENT(response)) THEN
         CPASSERT(ALL(SHAPE(response) == SHAPE(value)))
         response = CMPLX(0.0_dp, 0.0_dp, dp)
      END IF
      IF (PRESENT(gradient_response)) THEN
         CPASSERT(ALL(SHAPE(gradient_response) == [npoint, basis%nsgf, 3]))
         gradient_response = CMPLX(0.0_dp, 0.0_dp, dp)
      END IF
      IF (PRESENT(cell_gradient)) THEN
         CPASSERT(ALL(SHAPE(cell_gradient) == [npoint, basis%nsgf, 3, 3]))
         cell_gradient = CMPLX(0.0_dp, 0.0_dp, dp)
      END IF
      dc = 0.0_dp
      dcenter = 0.0_dp
      dk = 0.0_dp
      IF (PRESENT(cell_response)) dc = cell_response
      IF (PRESENT(center_response)) dcenter = center_response
      IF (PRESENT(reduced_k_response)) dk = reduced_k_response
      scratch => local_workspace
      IF (PRESENT(workspace)) scratch => workspace
      CALL prepare_workspace(scratch, npoint, basis%nsgf, need_geometry, need_gradient, need_hessian)
      displacement => scratch%displacement
      field => scratch%field
      NULLIFY (ddisplacement, g, h)
      IF (need_geometry) ddisplacement => scratch%ddisplacement
      IF (need_gradient) g => scratch%gradient
      IF (need_hessian) h => scratch%hessian
      IF (PRESENT(point_response)) THEN
         CPASSERT(ALL(SHAPE(point_response) == SHAPE(points)))
      END IF
      value = CMPLX(0.0_dp, 0.0_dp, dp)
      DO image = 1, SIZE(images, 2)
         lattice = REAL(images(:, image), dp)
         shift = center + MATMUL(cell, lattice)
         DO point = 1, npoint
            displacement(:, point) = points(:, point) - shift
         END DO
         IF (need_geometry) THEN
            dshift = -(dcenter + MATMUL(dc, lattice))
            DO point = 1, npoint
               ddisplacement(:, point) = dshift
            END DO
            IF (PRESENT(point_response)) ddisplacement(:, :) = ddisplacement + point_response
         END IF
         phase = EXP(CMPLX(0.0_dp, twopi*DOT_PRODUCT(reduced_k, lattice), dp))
         phase_response = CMPLX(0.0_dp, twopi*DOT_PRODUCT(dk, lattice), dp)*phase
         IF (need_hessian) THEN
            CALL gaussian_basis_fields_evaluate(basis, displacement, field, g, h, workspace=scratch)
         ELSE IF (need_gradient) THEN
            CALL gaussian_basis_fields_evaluate(basis, displacement, field, g, workspace=scratch)
         ELSE
            CALL gaussian_basis_fields_evaluate(basis, displacement, field, workspace=scratch)
         END IF
         value = value + phase*field
         IF (PRESENT(gradient)) gradient = gradient + phase*g
         IF (PRESENT(cell_gradient)) THEN
            DO axis2 = 1, 3
               DO axis = 1, 3
                  cell_gradient(:, :, axis, axis2) = cell_gradient(:, :, axis, axis2) - phase*g(:, :, axis)*lattice(axis2)
               END DO
            END DO
         END IF
         IF (PRESENT(response)) THEN
            response = response + phase_response*field
            DO axis = 1, 3
               DO ao = 1, basis%nsgf
                  response(:, ao) = response(:, ao) + phase*ddisplacement(axis, :)*g(:, ao, axis)
               END DO
            END DO
         END IF
         IF (PRESENT(gradient_response)) THEN
            gradient_response = gradient_response + phase_response*g
            DO axis = 1, 3
               DO axis2 = 1, 3
                  DO ao = 1, basis%nsgf
                     gradient_response(:, ao, axis) = gradient_response(:, ao, axis) + &
                                                      phase*ddisplacement(axis2, :)*h(:, ao, axis, axis2)
                  END DO
               END DO
            END DO
         END IF
      END DO
   END SUBROUTINE gaussian_basis_fields_periodic

! **************************************************************************************************
!> \brief Retain only scratch storage. No geometry or field values are reused.
!> \param workspace ...
!> \param npoint ...
!> \param nao ...
!> \param need_geometry ...
!> \param need_gradient ...
!> \param need_hessian ...
! **************************************************************************************************
   SUBROUTINE prepare_workspace(workspace, npoint, nao, need_geometry, need_gradient, need_hessian)
      TYPE(gaussian_basis_workspace_type), INTENT(INOUT) :: workspace
      INTEGER, INTENT(IN)                                :: npoint, nao
      LOGICAL, INTENT(IN)                                :: need_geometry, need_gradient, &
                                                            need_hessian

      IF (ALLOCATED(workspace%field)) THEN
         IF (ANY(SHAPE(workspace%field) /= [npoint, nao])) THEN
            DEALLOCATE (workspace%displacement, workspace%field)
            IF (ALLOCATED(workspace%ddisplacement)) DEALLOCATE (workspace%ddisplacement)
            IF (ALLOCATED(workspace%gradient)) DEALLOCATE (workspace%gradient)
            IF (ALLOCATED(workspace%hessian)) DEALLOCATE (workspace%hessian)
         END IF
      END IF
      IF (.NOT. ALLOCATED(workspace%field)) THEN
         ALLOCATE (workspace%displacement(3, npoint), workspace%field(npoint, nao))
      END IF
      IF (need_geometry .AND. .NOT. ALLOCATED(workspace%ddisplacement)) THEN
         ALLOCATE (workspace%ddisplacement(3, npoint))
      END IF
      IF (need_gradient .AND. .NOT. ALLOCATED(workspace%gradient)) THEN
         ALLOCATE (workspace%gradient(npoint, nao, 3))
      END IF
      IF (need_hessian .AND. .NOT. ALLOCATED(workspace%hessian)) THEN
         ALLOCATE (workspace%hessian(npoint, nao, 3, 3))
      END IF
   END SUBROUTINE prepare_workspace

! **************************************************************************************************
!> \brief Reuse Cartesian polynomial scratch storage independently on every OpenMP thread.
!>        All entries used by the requested derivatives are refilled for each point.
!> \param workspace ...
!> \param maxl ...
!> \param threads ...
!> \param need_gradient ...
!> \param need_hessian ...
! **************************************************************************************************
   SUBROUTINE prepare_polynomial_workspace(workspace, maxl, threads, need_gradient, need_hessian)
      TYPE(gaussian_basis_workspace_type), INTENT(INOUT) :: workspace
      INTEGER, INTENT(IN)                                :: maxl, threads
      LOGICAL, INTENT(IN)                                :: need_gradient, need_hessian

      IF (ALLOCATED(workspace%coordinate_power)) THEN
         IF (ANY(SHAPE(workspace%coordinate_power) /= [3, maxl + 1, threads])) THEN
            DEALLOCATE (workspace%coordinate_power, workspace%cartesian_value)
            IF (ALLOCATED(workspace%cartesian_gradient)) DEALLOCATE (workspace%cartesian_gradient)
            IF (ALLOCATED(workspace%cartesian_hessian)) DEALLOCATE (workspace%cartesian_hessian)
         END IF
      END IF
      IF (.NOT. ALLOCATED(workspace%coordinate_power)) THEN
         ALLOCATE (workspace%coordinate_power(3, 0:maxl, threads), workspace%cartesian_value(ncoset(maxl), threads))
      END IF
      IF (need_gradient .AND. .NOT. ALLOCATED(workspace%cartesian_gradient)) THEN
         ALLOCATE (workspace%cartesian_gradient(3, ncoset(maxl), threads))
      END IF
      IF (need_hessian .AND. .NOT. ALLOCATED(workspace%cartesian_hessian)) THEN
         ALLOCATE (workspace%cartesian_hessian(3, 3, ncoset(maxl), threads))
      END IF
   END SUBROUTINE prepare_polynomial_workspace

END MODULE gaussian_basis_fields
