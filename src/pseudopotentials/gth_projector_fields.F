!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Original GTH Cartesian projectors through the common analytic Gaussian field evaluator.
!>        Keep the native cprj polynomial coefficients and vprj_ppnl coupling without fitting,
!>        renormalization or additional spatial truncation. The stored radius is image-bound metadata.
! **************************************************************************************************
MODULE gth_projector_fields
   USE ai_overlap,                      ONLY: overlap_ab
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE external_potential_types,        ONLY: gth_potential_type
   USE gaussian_basis_fields,           ONLY: gaussian_basis_fields_evaluate,&
                                              gaussian_basis_fields_periodic
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: twopi
   USE orbital_pointers,                ONLY: init_orbital_pointers,&
                                              nco,&
                                              ncoset
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: gth_projector_fields_type, gth_projector_fields_init, gth_projector_basis
   PUBLIC :: gth_projector_fields_evaluate, gth_projector_fields_periodic, gth_projector_overlap

   TYPE gth_projector_fields_type
      INTEGER :: nfunction = 0
      REAL(dp) :: radius = 0.0_dp
      REAL(dp), ALLOCATABLE :: coupling(:, :)
      INTEGER, PRIVATE :: nset = 0
      INTEGER, ALLOCATABLE, PRIVATE :: lmax(:), npgf(:), nshell(:), l(:, :), first(:, :), last(:, :)
      REAL(dp), ALLOCATABLE, PRIVATE :: zet(:, :), sphi(:, :)
   END TYPE gth_projector_fields_type
CONTAINS
! **************************************************************************************************
!> \brief Own the exact native projector polynomials and their Cartesian coupling matrix.
!> \param potential Initialized scalar native GTH potential.
!> \param fields Owned immutable projector representation.
! **************************************************************************************************
   SUBROUTINE gth_projector_fields_init(potential, fields)
      TYPE(gth_potential_type), INTENT(IN)               :: potential
      TYPE(gth_projector_fields_type), INTENT(OUT)       :: fields

      INTEGER                                            :: angular, degree, first, iset, maximum, &
                                                            n, radial, total

      CPASSERT(.NOT. potential%soc)
      fields%nfunction = potential%nppnl
      fields%radius = potential%ppnl_radius
      ALLOCATE (fields%coupling(fields%nfunction, fields%nfunction))
      IF (fields%nfunction == 0) RETURN
      CPASSERT(ASSOCIATED(potential%nprj_ppnl))
      CPASSERT(ASSOCIATED(potential%alpha_ppnl))
      CPASSERT(ASSOCIATED(potential%cprj))
      CPASSERT(ASSOCIATED(potential%vprj_ppnl))
      fields%nset = COUNT(potential%nprj_ppnl(0:potential%lppnl) > 0)
      CPASSERT(fields%nset > 0)
      CALL init_orbital_pointers(MAX(0, potential%lprj_ppnl_max))
      maximum = 0
      total = 0
      DO angular = 0, potential%lppnl
         n = potential%nprj_ppnl(angular)
         IF (n == 0) CYCLE
         CPASSERT(n > 0)
         CPASSERT(potential%alpha_ppnl(angular) > 0.0_dp)
         maximum = MAX(maximum, angular + 2*(n - 1))
         total = total + n*nco(angular)
      END DO
      CPASSERT(total == fields%nfunction)
      CALL init_orbital_pointers(maximum)
      CPASSERT(SIZE(potential%cprj, 1) >= ncoset(maximum))
      CPASSERT(SIZE(potential%cprj, 2) == fields%nfunction)
      n = fields%nfunction
      CPASSERT(SIZE(potential%vprj_ppnl, 1) == n)
      CPASSERT(SIZE(potential%vprj_ppnl, 2) == n)
      n = MAXVAL(potential%nprj_ppnl)
      ALLOCATE (fields%lmax(fields%nset), fields%npgf(fields%nset), fields%nshell(fields%nset), &
                fields%zet(1, fields%nset), fields%sphi(ncoset(maximum), fields%nfunction), &
                fields%coupling(fields%nfunction, fields%nfunction), &
                fields%l(n, fields%nset), fields%first(n, fields%nset), fields%last(n, fields%nset))
      fields%coupling(:, :) = potential%vprj_ppnl
      fields%sphi(:, :) = potential%cprj(:ncoset(maximum), :)
      fields%npgf = 1
      fields%l = 0
      fields%first = 0
      fields%last = 0
      iset = 0
      first = 1
      DO angular = 0, potential%lppnl
         n = potential%nprj_ppnl(angular)
         IF (n == 0) CYCLE
         iset = iset + 1
         degree = angular + 2*(n - 1)
         fields%lmax(iset) = degree
         fields%nshell(iset) = n
         fields%zet(1, iset) = potential%alpha_ppnl(angular)
         DO radial = 1, n
            fields%l(radial, iset) = angular + 2*(radial - 1)
            fields%first(radial, iset) = first
            first = first + nco(angular)
            fields%last(radial, iset) = first - 1
         END DO
      END DO
      CPASSERT(first == fields%nfunction + 1)
   END SUBROUTINE gth_projector_fields_init
! **************************************************************************************************
!> \brief Borrow an evaluator view for this call only. The field object owns every array.
!> \param fields Initialized projector polynomials.
!> \param basis Borrowed Gaussian contraction view, never released as an owning basis.
! **************************************************************************************************
   SUBROUTINE gth_projector_basis(fields, basis)
      TYPE(gth_projector_fields_type), INTENT(IN), &
         TARGET                                          :: fields
      TYPE(gto_basis_set_type), INTENT(OUT)              :: basis

      basis%nset = fields%nset
      basis%nsgf = fields%nfunction
      basis%kind_radius = fields%radius
      basis%lmax => fields%lmax
      basis%npgf => fields%npgf
      basis%nshell => fields%nshell
      basis%l => fields%l
      basis%first_sgf => fields%first
      basis%last_sgf => fields%last
      basis%zet => fields%zet
      basis%sphi => fields%sphi
   END SUBROUTINE gth_projector_basis
! **************************************************************************************************
!> \brief Original projector values, Cartesian gradients and Hessians at point-minus-center vectors.
!> \param fields Initialized original GTH projector polynomials.
!> \param displacements Point-minus-center vectors, indexed (axis,point).
!> \param value New point-major native Cartesian projector fields.
!> \param gradient Optional point gradients, Cartesian axis last.
!> \param hessian Optional point Hessians, Cartesian axes last.
! **************************************************************************************************
   SUBROUTINE gth_projector_fields_evaluate(fields, displacements, value, gradient, hessian)
      TYPE(gth_projector_fields_type), INTENT(IN), &
         TARGET                                          :: fields
      REAL(dp), INTENT(IN)                               :: displacements(:, :)
      REAL(dp), INTENT(OUT)                              :: value(:, :)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: gradient(:, :, :), hessian(:, :, :, :)

      TYPE(gto_basis_set_type)                           :: basis

      CPASSERT(SIZE(displacements, 1) == 3)
      CPASSERT(SIZE(value, 1) == SIZE(displacements, 2))
      CPASSERT(SIZE(value, 2) == fields%nfunction)
      value = 0.0_dp
      IF (PRESENT(gradient)) gradient = 0.0_dp
      IF (PRESENT(hessian)) hessian = 0.0_dp
      IF (fields%nfunction == 0) RETURN
      CALL gth_projector_basis(fields, basis)
      CALL gaussian_basis_fields_evaluate(basis, displacements, value, gradient, hessian)
   END SUBROUTINE gth_projector_fields_evaluate
! **************************************************************************************************
!> \brief Sum every supplied physical image and retain the analytic geometry and Bloch responses.
!> \param fields Initialized native projector polynomials.
!> \param points Evaluation points in Cartesian coordinates.
!> \param center Physical projector center.
!> \param cell Lattice vectors in columns.
!> \param images Complete fixed integer image list.
!> \param reduced_k Reduced Bloch vector.
!> \param value Complex point-major native Cartesian projector values.
!> \param gradient Optional spatial gradients.
!> \param response Optional joint geometry and Bloch response.
!> \param gradient_response Optional joint response of the spatial gradients.
!> \param point_response Optional evaluation-point motion.
!> \param center_response Optional projector-center motion.
!> \param cell_response Optional cell-vector response.
!> \param reduced_k_response Optional reduced Bloch-vector response.
!> \param cell_gradient Optional nine cell derivatives at fixed Cartesian points and center.
! **************************************************************************************************
   SUBROUTINE gth_projector_fields_periodic(fields, points, center, cell, images, reduced_k, value, &
                                            gradient, response, gradient_response, point_response, center_response, &
                                            cell_response, reduced_k_response, cell_gradient)
      TYPE(gth_projector_fields_type), INTENT(IN), &
         TARGET                                          :: fields
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

      TYPE(gto_basis_set_type)                           :: basis

      CPASSERT(SIZE(points, 1) == 3)
      CPASSERT(SIZE(images, 1) == 3)
      CPASSERT(SIZE(value, 1) == SIZE(points, 2))
      CPASSERT(SIZE(value, 2) == fields%nfunction)
      value = CMPLX(0.0_dp, 0.0_dp, dp)
      IF (PRESENT(gradient)) gradient = CMPLX(0.0_dp, 0.0_dp, dp)
      IF (PRESENT(response)) response = CMPLX(0.0_dp, 0.0_dp, dp)
      IF (PRESENT(gradient_response)) gradient_response = CMPLX(0.0_dp, 0.0_dp, dp)
      IF (PRESENT(cell_gradient)) cell_gradient = CMPLX(0.0_dp, 0.0_dp, dp)
      IF (fields%nfunction == 0) RETURN
      CALL gth_projector_basis(fields, basis)
      CALL gaussian_basis_fields_periodic(basis, points, center, cell, images, reduced_k, value, gradient, &
                                          response, gradient_response, point_response, center_response, &
                                          cell_response, reduced_k_response, cell_gradient)
   END SUBROUTINE gth_projector_fields_periodic
! **************************************************************************************************
!> \brief Analytic native Gaussian overlaps against the original GTH Cartesian projector polynomials.
!>        Retain every supplied image and the native contractions. No radial quadrature or AO
!>        projection of a reconstruction field enters this operator. Derivatives move both centers,
!>        lattice images and the reduced Bloch vector in this same finite representation.
!> \param fields Original GTH projector representation.
!> \param basis Native contracted orbital basis.
!> \param center Projector center.
!> \param ao_center Orbital-basis center.
!> \param cell Lattice vectors in columns.
!> \param images Complete fixed physical image list.
!> \param reduced_k Reduced Bloch vector.
!> \param overlaps Complex projector bras against Gaussian columns.
!> \param response Optional joint geometry and Bloch response.
!> \param center_response Optional projector-center motion.
!> \param ao_response Optional orbital-center motion.
!> \param cell_response Optional cell-vector motion.
!> \param reduced_k_response Optional reduced Bloch-vector motion.
! **************************************************************************************************
   SUBROUTINE gth_projector_overlap(fields, basis, center, ao_center, cell, images, reduced_k, overlaps, &
                                    response, center_response, ao_response, cell_response, reduced_k_response)
      TYPE(gth_projector_fields_type), INTENT(IN), &
         TARGET                                          :: fields
      TYPE(gto_basis_set_type), INTENT(IN)               :: basis
      REAL(dp), INTENT(IN)                               :: center(3), ao_center(3), cell(3, 3)
      INTEGER, INTENT(IN)                                :: images(:, :)
      REAL(dp), INTENT(IN)                               :: reduced_k(3)
      COMPLEX(dp), INTENT(OUT)                           :: overlaps(:, :)
      COMPLEX(dp), INTENT(OUT), OPTIONAL                 :: response(:, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: center_response(3), ao_response(3), &
                                                            cell_response(3, 3), &
                                                            reduced_k_response(3)

      COMPLEX(dp)                                        :: phase, phase_response
      INTEGER                                            :: a, axis, b, first, image, last, na, nb, &
                                                            nprim, projector_first, projector_last
      REAL(dp)                                           :: displacement(3), dk(3), motion(3)
      REAL(dp), ALLOCATABLE                              :: contracted(:, :), derivative(:, :, :), &
                                                            native(:, :), radii(:), &
                                                            response_block(:, :)
      TYPE(gto_basis_set_type)                           :: projectors

      CPASSERT(SIZE(images, 1) == 3)
      CPASSERT(ALL(SHAPE(overlaps) == [fields%nfunction, basis%nsgf]))
      overlaps = CMPLX(0.0_dp, 0.0_dp, dp)
      IF (PRESENT(response)) THEN
         CPASSERT(ALL(SHAPE(response) == SHAPE(overlaps)))
         response = CMPLX(0.0_dp, 0.0_dp, dp)
      END IF
      IF (fields%nfunction == 0) RETURN
      CALL gth_projector_basis(fields, projectors)
      dk = 0.0_dp
      IF (PRESENT(reduced_k_response)) dk = reduced_k_response
      DO image = 1, SIZE(images, 2)
         ! The native Cartesian overlap recurrence uses center B minus center A.
         displacement = ao_center + MATMUL(cell, REAL(images(:, image), dp)) - center
         motion = 0.0_dp
         IF (PRESENT(center_response)) motion = motion + center_response
         IF (PRESENT(ao_response)) motion = motion - ao_response
         IF (PRESENT(cell_response)) motion = motion - MATMUL(cell_response, REAL(images(:, image), dp))
         phase = EXP(CMPLX(0.0_dp, twopi*DOT_PRODUCT(reduced_k, REAL(images(:, image), dp)), dp))
         phase_response = CMPLX(0.0_dp, twopi*DOT_PRODUCT(dk, REAL(images(:, image), dp)), dp)*phase
         DO a = 1, projectors%nset
            na = ncoset(projectors%lmax(a))
            projector_first = MINVAL(projectors%first_sgf(:projectors%nshell(a), a))
            projector_last = MAXVAL(projectors%last_sgf(:projectors%nshell(a), a))
            DO b = 1, basis%nset
               nprim = basis%npgf(b)
               nb = nprim*ncoset(basis%lmax(b))
               first = MINVAL(basis%first_sgf(:basis%nshell(b), b))
               last = MAXVAL(basis%last_sgf(:basis%nshell(b), b))
               ALLOCATE (native(na, nb), contracted(projector_last - projector_first + 1, last - first + 1), radii(MAX(1, nprim)))
               radii = HUGE(0.0_dp)/4.0_dp
               IF (PRESENT(response)) THEN
                  ALLOCATE (derivative(na, nb, 3), response_block(projector_last - projector_first + 1, last - first + 1))
                  CALL overlap_ab(projectors%lmax(a), 0, 1, radii(:1), projectors%zet(:, a), &
                                  basis%lmax(b), 0, nprim, radii, basis%zet(:, b), displacement, &
                                  sab=native, dab=derivative)
               ELSE
                  CALL overlap_ab(projectors%lmax(a), 0, 1, radii(:1), projectors%zet(:, a), &
                                  basis%lmax(b), 0, nprim, radii, basis%zet(:, b), displacement, sab=native)
               END IF
               contracted(:, :) = MATMUL(TRANSPOSE(projectors%sphi(:na, projector_first:projector_last)), &
                                         MATMUL(native, basis%sphi(:nb, first:last)))
               overlaps(projector_first:projector_last, first:last) = &
                  overlaps(projector_first:projector_last, first:last) + phase*contracted
               IF (PRESENT(response)) THEN
                  response_block = 0.0_dp
                  DO axis = 1, 3
      response_block(:, :) = response_block + motion(axis)*MATMUL(TRANSPOSE(projectors%sphi(:na, projector_first:projector_last)), &
                                                                                 MATMUL(derivative(:, :, axis), &
                                                                                        basis%sphi(:nb, first:last)))
                  END DO
                  response(projector_first:projector_last, first:last) = &
                     response(projector_first:projector_last, first:last) + phase_response*contracted + phase*response_block
                  DEALLOCATE (derivative, response_block)
               END IF
               DEALLOCATE (native, contracted, radii)
            END DO
         END DO
      END DO
   END SUBROUTINE gth_projector_overlap
END MODULE gth_projector_fields
