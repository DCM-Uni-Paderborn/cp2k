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
   USE kinds,                           ONLY: dp,&
                                              int_8
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
               IF (exponential == 0.0_dp) CYCLE
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
                        weight = basis%sphi(irow, isgf)
                        IF (weight == 0.0_dp) CYCLE
                        weight = weight*exponential
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
      LOGICAL                                            :: factored, need_geometry, need_gradient, &
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
      IF (PRESENT(point_response)) THEN
         CPASSERT(ALL(SHAPE(point_response) == SHAPE(points)))
      END IF
      CALL periodic_cartesian_product(basis, points, center, cell, images, reduced_k, &
                                      value, gradient, cell_gradient, factored, &
                                      response, point_response, dcenter, dc, dk, gradient_response)
      IF (factored) RETURN
      scratch => local_workspace
      IF (PRESENT(workspace)) scratch => workspace
      CALL prepare_workspace(scratch, npoint, basis%nsgf, need_geometry, need_gradient, need_hessian)
      displacement => scratch%displacement
      field => scratch%field
      NULLIFY (ddisplacement, g, h)
      IF (need_geometry) ddisplacement => scratch%ddisplacement
      IF (need_gradient) g => scratch%gradient
      IF (need_hessian) h => scratch%hessian
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
!> \brief Factor a complete Cartesian translation box in an orthorhombic cell.
!>        Products of one-dimensional Gaussian moments retain every supplied image and phase.
!>        Image-weighted moments give all nine cell derivatives, including shear at this cell.
!>        Other cells and image layouts return without changing the caller's fields.
!> \param basis Original native contraction and Cartesian monomial ordering.
!> \param points Evaluation coordinates, indexed (axis,point).
!> \param center Physical Gaussian center.
!> \param cell Lattice vectors in columns.
!> \param images Exact image box in the native x-fast ordering.
!> \param reduced_k Fixed reduced Bloch vector.
!> \param value New fields if factored, indexed (point,AO).
!> \param gradient Optional spatial gradients, indexed (point,AO,axis).
!> \param cell_gradient Optional cell derivatives, indexed (point,AO,axis,column).
!> \param factored Whether the supplied cell and image list permit this factorization.
!> \param response Optional joint first response of the fields.
!> \param point_response Optional motion of the evaluation points.
!> \param center_response Motion of the Gaussian center.
!> \param cell_response Lattice derivative at fixed reduced Bloch coordinates.
!> \param reduced_k_response Reduced Bloch-vector derivative.
!> \param gradient_response Optional joint response of the spatial gradients.
! **************************************************************************************************
   SUBROUTINE periodic_cartesian_product(basis, points, center, cell, images, reduced_k, &
                                         value, gradient, cell_gradient, factored, &
                                         response, point_response, center_response, cell_response, reduced_k_response, &
                                         gradient_response)
      TYPE(gto_basis_set_type), INTENT(IN)               :: basis
      REAL(dp), INTENT(IN)                               :: points(:, :), center(3), cell(3, 3)
      INTEGER, INTENT(IN)                                :: images(:, :)
      REAL(dp), INTENT(IN)                               :: reduced_k(3)
      COMPLEX(dp), INTENT(INOUT)                         :: value(:, :)
      COMPLEX(dp), INTENT(INOUT), OPTIONAL               :: gradient(:, :, :), &
                                                            cell_gradient(:, :, :, :)
      LOGICAL, INTENT(OUT)                               :: factored
      COMPLEX(dp), INTENT(INOUT), OPTIONAL               :: response(:, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: point_response(:, :)
      REAL(dp), INTENT(IN)                               :: center_response(3), cell_response(3, 3), &
                                                            reduced_k_response(3)
      COMPLEX(dp), INTENT(INOUT), OPTIONAL               :: gradient_response(:, :, :)

      COMPLEX(dp)                                        :: phase
      COMPLEX(dp), ALLOCATABLE :: derivative(:, :), hessian_component(:), image_moments(:, :, :), &
         moments(:, :, :), phase_moment(:), primitive_cell(:, :, :), primitive_gradient(:, :), &
         primitive_gradient_response(:, :), primitive_response(:), primitive_value(:)
      INTEGER :: axis, column, degree, expected(3), first, ico, image, ipgf, irow, iset, isgf, &
         ishell, l, last, lower(3), maxl, moment_degree, ncart, npoint, other, powers(3), &
         translation, upper(3), width(3)
      LOGICAL                                            :: cell_derivatives, derivatives, &
                                                            image_derivatives
      REAL(dp)                                           :: alpha, weight
      REAL(dp), ALLOCATABLE                              :: displacement(:), exponential(:), power(:)

      factored = .FALSE.
      IF (SIZE(images, 2) < 8) RETURN
      DO column = 1, 3
         DO axis = 1, 3
            IF (axis /= column .AND. cell(axis, column) /= 0.0_dp) RETURN
         END DO
      END DO
      lower = MINVAL(images, DIM=2)
      upper = MAXVAL(images, DIM=2)
      width = upper - lower + 1
      IF (ANY(width > SIZE(images, 2))) RETURN
      IF (INT(width(1), int_8)*INT(width(2), int_8)*INT(width(3), int_8) /= SIZE(images, 2, KIND=int_8)) RETURN
      DO image = 1, SIZE(images, 2)
         expected = lower + [MODULO(image - 1, width(1)), &
                             MODULO((image - 1)/width(1), width(2)), (image - 1)/(width(1)*width(2))]
         IF (ANY(images(:, image) /= expected)) RETURN
      END DO
      CPASSERT(basis%nset > 0 .AND. basis%nsgf > 0)
      CPASSERT(ASSOCIATED(basis%sphi) .AND. ASSOCIATED(basis%zet))
      maxl = MAXVAL(basis%lmax)
      CPASSERT(UBOUND(ncoset, 1) >= maxl)
      npoint = SIZE(points, 2)
      cell_derivatives = PRESENT(cell_gradient)
      IF (PRESENT(response)) cell_derivatives = cell_derivatives .OR. ANY(cell_response /= 0.0_dp)
      ! Gradient responses use image-weighted Hessians and Bloch derivatives of the spatial gradient.
      IF (PRESENT(gradient_response)) cell_derivatives = .TRUE.
      image_derivatives = cell_derivatives
      IF (PRESENT(response)) image_derivatives = image_derivatives .OR. ANY(reduced_k_response /= 0.0_dp)
      derivatives = PRESENT(gradient) .OR. cell_derivatives .OR. PRESENT(response)
      moment_degree = maxl + 1
      IF (PRESENT(gradient_response)) moment_degree = maxl + 2
      ALLOCATE (moments(npoint, 0:moment_degree, 3), image_moments(npoint, 0:moment_degree, 3), &
                derivative(npoint, 3), primitive_value(npoint), primitive_gradient(npoint, 3), &
                primitive_cell(npoint, 3, 3), displacement(npoint), exponential(npoint), power(npoint))
      value = CMPLX(0.0_dp, 0.0_dp, dp)
      IF (PRESENT(gradient)) gradient = CMPLX(0.0_dp, 0.0_dp, dp)
      IF (PRESENT(cell_gradient)) cell_gradient = CMPLX(0.0_dp, 0.0_dp, dp)
      IF (PRESENT(response)) THEN
         ALLOCATE (primitive_response(npoint), phase_moment(npoint))
         response = CMPLX(0.0_dp, 0.0_dp, dp)
      END IF
      IF (PRESENT(gradient_response)) THEN
         ALLOCATE (primitive_gradient_response(npoint, 3), hessian_component(npoint))
         gradient_response = CMPLX(0.0_dp, 0.0_dp, dp)
      END IF
      DO iset = 1, basis%nset
         ncart = ncoset(basis%lmax(iset))
         DO ipgf = 1, basis%npgf(iset)
            alpha = basis%zet(ipgf, iset)
            moments = CMPLX(0.0_dp, 0.0_dp, dp)
            image_moments = CMPLX(0.0_dp, 0.0_dp, dp)
            DO axis = 1, 3
               DO translation = lower(axis), upper(axis)
                  phase = EXP(CMPLX(0.0_dp, twopi*reduced_k(axis)*REAL(translation, dp), dp))
                  displacement(:) = points(axis, :) - center(axis) - cell(axis, axis)*REAL(translation, dp)
                  exponential(:) = EXP(-alpha*displacement**2)
                  power = 1.0_dp
                  WHERE (exponential == 0.0_dp) power = 0.0_dp
                  DO degree = 0, moment_degree
                     moments(:, degree, axis) = moments(:, degree, axis) + phase*power*exponential
                     IF (image_derivatives) THEN
                        image_moments(:, degree, axis) = image_moments(:, degree, axis) + &
                                                         REAL(translation, dp)*phase*power*exponential
                     END IF
                     power(:) = power*displacement
                  END DO
               END DO
            END DO
            DO ishell = 1, basis%nshell(iset)
               l = basis%l(ishell, iset)
               first = basis%first_sgf(ishell, iset)
               last = basis%last_sgf(ishell, iset)
               DO ico = ncoset(l - 1) + 1, ncoset(l)
                  irow = (ipgf - 1)*ncart + ico
                  powers = indco(:, ico)
                  primitive_value(:) = moments(:, powers(1), 1)*moments(:, powers(2), 2)*moments(:, powers(3), 3)
                  IF (derivatives) THEN
                     DO axis = 1, 3
                        derivative(:, axis) = -2.0_dp*alpha*moments(:, powers(axis) + 1, axis)
                        IF (powers(axis) > 0) THEN
                           derivative(:, axis) = derivative(:, axis) + powers(axis)*moments(:, powers(axis) - 1, axis)
                        END IF
                        primitive_gradient(:, axis) = derivative(:, axis)
                        DO other = 1, 3
                           IF (other == axis) CYCLE
                           primitive_gradient(:, axis) = primitive_gradient(:, axis)*moments(:, powers(other), other)
                        END DO
                     END DO
                  END IF
                  IF (cell_derivatives) THEN
                     DO column = 1, 3
                        DO axis = 1, 3
                           IF (axis == column) THEN
                              primitive_cell(:, axis, column) = 2.0_dp*alpha*image_moments(:, powers(axis) + 1, axis)
                              IF (powers(axis) > 0) THEN
                                 primitive_cell(:, axis, column) = primitive_cell(:, axis, column) - &
                                                                   powers(axis)*image_moments(:, powers(axis) - 1, axis)
                              END IF
                           ELSE
                              primitive_cell(:, axis, column) = -derivative(:, axis)*image_moments(:, powers(column), column)
                           END IF
                           DO other = 1, 3
                              IF (other == axis .OR. other == column) CYCLE
                              primitive_cell(:, axis, column) = primitive_cell(:, axis, column)*moments(:, powers(other), other)
                           END DO
                        END DO
                     END DO
                  END IF
                  IF (PRESENT(response)) THEN
                     primitive_response = CMPLX(0.0_dp, 0.0_dp, dp)
                     DO axis = 1, 3
                        primitive_response(:) = primitive_response - center_response(axis)*primitive_gradient(:, axis)
                        IF (PRESENT(point_response)) THEN
                           primitive_response(:) = primitive_response + point_response(axis, :)*primitive_gradient(:, axis)
                        END IF
                        IF (cell_derivatives) THEN
                           DO column = 1, 3
                            primitive_response(:) = primitive_response + cell_response(axis, column)*primitive_cell(:, axis, column)
                           END DO
                        END IF
                        IF (reduced_k_response(axis) /= 0.0_dp) THEN
                           phase_moment(:) = image_moments(:, powers(axis), axis)
                           DO other = 1, 3
                              IF (other == axis) CYCLE
                              phase_moment(:) = phase_moment*moments(:, powers(other), other)
                           END DO
                           primitive_response(:) = primitive_response + &
                                                   CMPLX(0.0_dp, twopi*reduced_k_response(axis), dp)*phase_moment
                        END IF
                     END DO
                  END IF
                  IF (PRESENT(gradient_response)) THEN
                     primitive_gradient_response = CMPLX(0.0_dp, 0.0_dp, dp)
                     DO axis = 1, 3
                        DO other = 1, 3
                           hessian_component(:) = factored_hessian(moments, image_moments, powers, alpha, axis, other, 0)
                           primitive_gradient_response(:, axis) = primitive_gradient_response(:, axis) - &
                                                                  center_response(other)*hessian_component
                           IF (PRESENT(point_response)) THEN
                              primitive_gradient_response(:, axis) = primitive_gradient_response(:, axis) + &
                                                                     point_response(other, :)*hessian_component
                           END IF
                           DO column = 1, 3
                              IF (cell_response(other, column) == 0.0_dp) CYCLE
                              hessian_component(:) = factored_hessian(moments, image_moments, powers, alpha, axis, other, column)
                              primitive_gradient_response(:, axis) = primitive_gradient_response(:, axis) - &
                                                                     cell_response(other, column)*hessian_component
                           END DO
                        END DO
                        DO column = 1, 3
                           IF (reduced_k_response(column) == 0.0_dp) CYCLE
                           primitive_gradient_response(:, axis) = primitive_gradient_response(:, axis) - &
                                                 CMPLX(0.0_dp, twopi*reduced_k_response(column), dp)*primitive_cell(:, axis, column)
                        END DO
                     END DO
                  END IF
                  DO isgf = first, last
                     weight = basis%sphi(irow, isgf)
                     IF (weight == 0.0_dp) CYCLE
                     value(:, isgf) = value(:, isgf) + weight*primitive_value
                     IF (PRESENT(gradient)) gradient(:, isgf, :) = gradient(:, isgf, :) + weight*primitive_gradient
                     IF (PRESENT(cell_gradient)) cell_gradient(:, isgf, :, :) = cell_gradient(:, isgf, :, :) + &
                                                                                weight*primitive_cell
                     IF (PRESENT(response)) response(:, isgf) = response(:, isgf) + weight*primitive_response
                     IF (PRESENT(gradient_response)) gradient_response(:, isgf, :) = gradient_response(:, isgf, :) + &
                                                                                     weight*primitive_gradient_response
                  END DO
               END DO
            END DO
         END DO
      END DO
      factored = .TRUE.
   END SUBROUTINE periodic_cartesian_product

! **************************************************************************************************
!> \brief Cartesian second derivative of a separable periodic primitive, optionally weighted by an image index.
!> \param moments One-dimensional periodic Gaussian polynomial moments.
!> \param image_moments The same moments multiplied by their integer image coordinate.
!> \param powers Cartesian polynomial degrees of the primitive.
!> \param alpha Primitive Gaussian exponent.
!> \param first First Cartesian differentiation axis.
!> \param second Second Cartesian differentiation axis.
!> \param image_axis Image coordinate for weighting, or zero for an unweighted Hessian.
!> \return Sum over the original image product, including its Bloch phase.
! **************************************************************************************************
   FUNCTION factored_hessian(moments, image_moments, powers, alpha, first, second, image_axis) RESULT(value)
      COMPLEX(dp), INTENT(IN), TARGET                    :: moments(:, 0:, :), &
                                                            image_moments(:, 0:, :)
      INTEGER, INTENT(IN)                                :: powers(3)
      REAL(dp), INTENT(IN)                               :: alpha
      INTEGER, INTENT(IN)                                :: first, second, image_axis
      COMPLEX(dp)                                        :: value(SIZE(moments, 1))

      COMPLEX(dp)                                        :: factor(SIZE(moments, 1))
      COMPLEX(dp), POINTER                               :: source(:, :)
      INTEGER                                            :: axis, degree, order

      value = CMPLX(1.0_dp, 0.0_dp, dp)
      DO axis = 1, 3
         source(1:, 0:) => moments(:, :, axis)
         IF (axis == image_axis) source(1:, 0:) => image_moments(:, :, axis)
         degree = powers(axis)
         order = MERGE(1, 0, axis == first) + MERGE(1, 0, axis == second)
         SELECT CASE (order)
         CASE (0)
            factor = source(:, degree)
         CASE (1)
            factor = -2.0_dp*alpha*source(:, degree + 1)
            IF (degree > 0) factor = factor + degree*source(:, degree - 1)
         CASE (2)
            factor = 4.0_dp*alpha**2*source(:, degree + 2) - 2.0_dp*alpha*(2*degree + 1)*source(:, degree)
            IF (degree > 1) factor = factor + degree*(degree - 1)*source(:, degree - 2)
         END SELECT
         value = value*factor
      END DO
   END FUNCTION factored_hessian

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
