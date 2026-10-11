!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Original inverse defining matrices resolved by physical lattice translation.
!>        C_Q=sum_T <beta_T|Delta_(T+Q)> and B_Q=sum_T <beta_T|chi_(T+Q)> on the
!>        actual distributed point quadrature. Their Fourier sums recover the original
!>        finite-image Bloch defining matrices without assigning Bloch W to individual images.
! **************************************************************************************************
MODULE qs_upf_inverse_lattice
   USE gaussian_basis_fields,           ONLY: gaussian_basis_fields_periodic
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: dp
   USE qs_upf_inverse_types,            ONLY: qs_upf_inverse_lattice_type,&
                                              qs_upf_inverse_type
   USE upf_inverse_lattice,             ONLY: upf_inverse_lattice_adjoint,&
                                              upf_inverse_lattice_error_bounds,&
                                              upf_inverse_lattice_init,&
                                              upf_inverse_lattice_project
   USE upf_inverse_quadrature,          ONLY: upf_inverse_pair,&
                                              upf_inverse_pair_response
   USE upf_radial_fields,               ONLY: upf_radial_fields_periodic
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: qs_upf_inverse_lattice_kernel
   PUBLIC :: qs_upf_inverse_lattice_kernel_adjoint
   PUBLIC :: qs_upf_inverse_lattice_fields, qs_upf_inverse_lattice_field_adjoint
   PUBLIC :: qs_upf_inverse_lattice_geometry_adjoint
   PUBLIC :: qs_upf_inverse_lattice_prepare, qs_upf_inverse_lattice_adjoint

CONTAINS

! **************************************************************************************************
!> \brief Own the original defining integrals and a controlled physical inverse coefficient sequence.
!>        The original point communicator computes C and B once. Complete B support remains
!>        separate from the solution domain. Increasing physical shells are accepted only when
!>        both the sequence and uniform Bloch error bounds satisfy the requested tolerance.
!>        These bounds concern the original finite defining quadrature and its supplied images.
!> \param frame Original common point distribution, unphased fields and native geometry.
!> \param state New owned convolution, original B and retained W sequence.
!> \param tolerance Absolute Frobenius bounds for the W sequence and every Bloch sum.
!> \param info Zero on success, one for LU failure, two for solve failure, three for an
!>        inconclusive uniform inverse bound, four for insufficient physical image support.
!> \param coefficient_images Optional fixed physical solution domain instead of adaptive shells.
!> \param max_shells Maximum adaptive shell extent in periodic directions, twelve by default.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_lattice_prepare(frame, state, tolerance, info, coefficient_images, max_shells)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      TYPE(qs_upf_inverse_lattice_type), INTENT(OUT)     :: state
      REAL(dp), INTENT(IN)                               :: tolerance
      INTEGER, INTENT(OUT)                               :: info
      INTEGER, INTENT(IN), OPTIONAL                      :: coefficient_images(:, :), max_shells

      COMPLEX(dp), ALLOCATABLE                           :: coupling(:, :, :), rhs(:, :, :)
      INTEGER                                            :: extent(3), i, image, ix, iy, iz, limit, &
                                                            nao, np, q, shell, status
      INTEGER, ALLOCATABLE                               :: images(:, :)

      CPASSERT(.NOT. frame%atom_local .AND. .NOT. frame%direct)
      CPASSERT(ieee_is_finite(tolerance) .AND. tolerance > 0.0_dp)
      state%tolerance = tolerance
      state%solve_tolerance = MIN(0.01_dp*tolerance, 1.0e-12_dp)
      CALL qs_upf_inverse_lattice_kernel(frame, state%overlap_images, coupling, state%overlaps)
      np = SIZE(coupling, 1)
      nao = SIZE(state%overlaps, 2)
      limit = 12
      IF (PRESENT(max_shells)) limit = max_shells
      CPASSERT(limit >= 0)
      IF (PRESENT(coefficient_images)) THEN
         CALL check_lattice_images(frame, coefficient_images)
         CPASSERT(SIZE(coefficient_images, 2) > 0)
         limit = 0
      END IF
      DO shell = 0, limit
         IF (PRESENT(coefficient_images)) THEN
            ALLOCATE (images, SOURCE=coefficient_images)
         ELSE
            extent = shell*MIN(frame%periodicity, MAXVAL(ABS(state%overlap_images), DIM=2))
            ALLOCATE (images(3, PRODUCT(2*extent + 1)))
            image = 0
            DO iz = -extent(3), extent(3)
            DO iy = -extent(2), extent(2)
            DO ix = -extent(1), extent(1)
               image = image + 1
               images(:, image) = [ix, iy, iz]
            END DO
            END DO
            END DO
         END IF
         CALL upf_inverse_lattice_init(state%overlap_images, coupling, images, state%inverse, status)
         info = 1
         IF (status /= 0) RETURN
         ALLOCATE (rhs(np, nao, SIZE(images, 2)), state%coefficients(np, nao, SIZE(images, 2)))
         rhs = CMPLX(0.0_dp, 0.0_dp, dp)
         DO i = 1, SIZE(images, 2)
            DO q = 1, SIZE(state%overlap_images, 2)
               IF (ANY(images(:, i) /= state%overlap_images(:, q))) CYCLE
               rhs(:, :, i) = state%overlaps(:, :, q)
               EXIT
            END DO
         END DO
         CALL upf_inverse_lattice_project(state%inverse, rhs, state%coefficients, state%solve_tolerance, status)
         info = 2
         IF (status /= 0) RETURN
         CALL upf_inverse_lattice_error_bounds(state%inverse, state%coefficients, state%overlap_images, &
                                               state%overlaps, state%coefficient_error, state%bloch_error, &
                                               state%lower_bound, status)
         info = 3
         IF (status /= 0) RETURN
         IF (MAX(state%coefficient_error, state%bloch_error) <= tolerance) THEN
            state%ready = .TRUE.
            info = 0
            RETURN
         END IF
         info = 4
         IF (shell == limit) RETURN
         DEALLOCATE (images, rhs, state%coefficients)
      END DO
   END SUBROUTINE qs_upf_inverse_lattice_prepare

! **************************************************************************************************
!> \brief Return the original defining-quadrature geometry transpose of the owned W sequence.
!>        The caller supplies the complete coefficient energy covector after its RI owner sum.
!>        Outputs belong to the original point owners and are not reduced here. Complete C_Q
!>        variations and the retained B_R rows are mapped back to the original defining support.
!> \param frame Unchanged original point distribution and native radial/Gaussian fields.
!> \param state Controlled physical convolution used for the forward RI integrals.
!> \param covector Complete W derivative, indexed (projector,AO,physical image).
!> \param point_geometry New owner-local Cartesian point derivative.
!> \param atom_geometry New owner-local Cartesian derivative in native atom order.
!> \param cell_geometry New owner-local cell-vector derivative at fixed Cartesian points/atoms.
!> \param weights_covector New owner-local original measure derivative, including zero weights.
!> \param info Zero on success, otherwise the coupled transpose solve status.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_lattice_adjoint(frame, state, covector, point_geometry, atom_geometry, &
                                             cell_geometry, weights_covector, info)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      TYPE(qs_upf_inverse_lattice_type), INTENT(IN)      :: state
      COMPLEX(dp), INTENT(IN)                            :: covector(:, :, :)
      REAL(dp), INTENT(OUT)                              :: point_geometry(:, :), &
                                                            atom_geometry(:, :), &
                                                            cell_geometry(3, 3), &
                                                            weights_covector(:)
      INTEGER, INTENT(OUT)                               :: info

      COMPLEX(dp), ALLOCATABLE                           :: coupling_covector(:, :, :), &
                                                            overlaps_covector(:, :, :), &
                                                            retained_covector(:, :, :)
      INTEGER                                            :: i, q

      CPASSERT(state%ready)
      CPASSERT(.NOT. frame%atom_local .AND. .NOT. frame%direct)
      CPASSERT(ALL(SHAPE(covector) == SHAPE(state%coefficients)))
      point_geometry = 0.0_dp
      atom_geometry = 0.0_dp
      cell_geometry = 0.0_dp
      weights_covector = 0.0_dp
      ALLOCATE (retained_covector, MOLD=state%coefficients)
      ALLOCATE (coupling_covector, MOLD=state%inverse%coupling)
      ALLOCATE (overlaps_covector, MOLD=state%overlaps)
      CALL upf_inverse_lattice_adjoint(state%inverse, state%coefficients, covector, retained_covector, &
                                       coupling_covector, state%solve_tolerance, info)
      IF (info /= 0) RETURN
      overlaps_covector = CMPLX(0.0_dp, 0.0_dp, dp)
      DO q = 1, SIZE(state%overlap_images, 2)
         DO i = 1, SIZE(state%inverse%images, 2)
            IF (ANY(state%overlap_images(:, q) /= state%inverse%images(:, i))) CYCLE
            overlaps_covector(:, :, q) = retained_covector(:, :, i)
            EXIT
         END DO
      END DO
      CALL qs_upf_inverse_lattice_kernel_adjoint(frame, state%inverse%kernel_images, coupling_covector, &
                                                 overlaps_covector, point_geometry, atom_geometry, &
                                                 cell_geometry, weights_covector)
   END SUBROUTINE qs_upf_inverse_lattice_adjoint

! **************************************************************************************************
!> \brief Evaluate a localized inverse AO image from its coupled physical-image coefficients.
!>        For AO image R, phi_R=chi_R-sum_Q Delta_(R-Q)*W_Q. W_Q retains the full
!>        projector matrix and the requested AO columns. No Bloch phase, minimum image,
!>        one-center projection, or further image truncation is applied. Scratch fields
!>        contain one physical image at a time. The caller owns the coefficient tail error.
!> \param frame Original basis, radial difference fields and physical geometry.
!> \param points Cartesian evaluation points, indexed (axis,point).
!> \param ao_image Physical integer AO translation R.
!> \param coefficient_images Distinct physical translations Q in coefficient order.
!> \param coefficients W_Q, indexed (projector,AO in ao_range,image).
!> \param ao_range Inclusive native AO column range.
!> \param value New localized soft values, indexed (point,AO in ao_range).
!> \param gradient Optional Cartesian gradients, axis last.
!> \param response Optional complete geometry response of the localized values.
!> \param gradient_response Optional response of the Cartesian gradients.
!> \param coefficients_response Coupled dW_Q, required for a geometry response.
!> \param point_response Physical point displacements, required for a geometry response.
!> \param atom_response Physical displacements of all original centers.
!> \param cell_response Cell-vector displacement at fixed integer images.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_lattice_fields(frame, points, ao_image, coefficient_images, coefficients, ao_range, &
                                            value, gradient, response, gradient_response, coefficients_response, &
                                            point_response, atom_response, cell_response)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      REAL(dp), INTENT(IN)                               :: points(:, :)
      INTEGER, INTENT(IN)                                :: ao_image(3), coefficient_images(:, :)
      COMPLEX(dp), INTENT(IN)                            :: coefficients(:, :, :)
      INTEGER, INTENT(IN)                                :: ao_range(2)
      COMPLEX(dp), INTENT(OUT)                           :: value(:, :)
      COMPLEX(dp), INTENT(OUT), OPTIONAL                 :: gradient(:, :, :), response(:, :), &
                                                            gradient_response(:, :, :)
      COMPLEX(dp), INTENT(IN), OPTIONAL                  :: coefficients_response(:, :, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: point_response(:, :), &
                                                            atom_response(:, :), &
                                                            cell_response(3, 3)

      COMPLEX(dp), ALLOCATABLE :: delta(:, :), delta_gradient(:, :, :), &
         delta_gradient_response(:, :, :), delta_response(:, :), gaussian(:, :), &
         gaussian_gradient(:, :, :), gaussian_gradient_response(:, :, :), gaussian_response(:, :)
      INTEGER                                            :: atom, axis, first, last, n, nao, np, q, &
                                                            single_image(3, 1), source_first, &
                                                            source_last
      LOGICAL                                            :: derivative, gradients

      CALL check_lattice_images(frame, coefficient_images)
      CPASSERT(ALL(ao_image*frame%periodicity == ao_image))
      CPASSERT(SIZE(points, 1) == 3)
      np = SIZE(frame%coefficients, 1)
      n = SIZE(points, 2)
      nao = ao_range(2) - ao_range(1) + 1
      CPASSERT(nao > 0 .AND. ao_range(1) > 0)
      CPASSERT(ao_range(2) <= SIZE(frame%gaussian, 2))
      CPASSERT(SIZE(coefficients, 1) == np)
      CPASSERT(SIZE(coefficients, 2) == nao)
      CPASSERT(SIZE(coefficients, 3) == SIZE(coefficient_images, 2))
      CPASSERT(ALL(SHAPE(value) == [n, nao]))
      derivative = PRESENT(response)
      gradients = PRESENT(gradient)
      CPASSERT(.NOT. PRESENT(gradient_response) .OR. (derivative .AND. gradients))
      IF (gradients) THEN
         CPASSERT(ALL(SHAPE(gradient) == [n, nao, 3]))
      END IF
      IF (derivative) THEN
         CPASSERT(PRESENT(coefficients_response))
         CPASSERT(PRESENT(point_response))
         CPASSERT(PRESENT(atom_response))
         CPASSERT(PRESENT(cell_response))
         CPASSERT(ALL(SHAPE(response) == SHAPE(value)))
         CPASSERT(ALL(SHAPE(coefficients_response) == SHAPE(coefficients)))
         CPASSERT(ALL(SHAPE(point_response) == SHAPE(points)))
         CPASSERT(ALL(SHAPE(atom_response) == SHAPE(frame%positions)))
         response = CMPLX(0.0_dp, 0.0_dp, dp)
      END IF
      IF (PRESENT(gradient_response)) THEN
         CPASSERT(ALL(SHAPE(gradient_response) == [n, nao, 3]))
         gradient_response = CMPLX(0.0_dp, 0.0_dp, dp)
      END IF
      value = CMPLX(0.0_dp, 0.0_dp, dp)
      IF (gradients) gradient = CMPLX(0.0_dp, 0.0_dp, dp)
      single_image(:, 1) = ao_image
      DO atom = 1, SIZE(frame%bases)
         first = MAX(frame%ao_first(atom), ao_range(1))
         last = MIN(frame%ao_last(atom), ao_range(2))
         IF (last < first) CYCLE
         ALLOCATE (gaussian(n, frame%bases(atom)%gto_basis_set%nsgf))
         IF (gradients) ALLOCATE (gaussian_gradient(n, SIZE(gaussian, 2), 3))
         IF (derivative) ALLOCATE (gaussian_response(n, SIZE(gaussian, 2)))
         IF (PRESENT(gradient_response)) ALLOCATE (gaussian_gradient_response(n, SIZE(gaussian, 2), 3))
         IF (derivative) THEN
            CALL gaussian_basis_fields_periodic(frame%bases(atom)%gto_basis_set, points, frame%positions(:, atom), &
                                                frame%cell, single_image, [0.0_dp, 0.0_dp, 0.0_dp], gaussian, &
                                                gaussian_gradient, gaussian_response, gaussian_gradient_response, &
                                                point_response, atom_response(:, atom), cell_response)
         ELSE
            CALL gaussian_basis_fields_periodic(frame%bases(atom)%gto_basis_set, points, frame%positions(:, atom), &
                                                frame%cell, single_image, [0.0_dp, 0.0_dp, 0.0_dp], gaussian, gaussian_gradient)
         END IF
         source_first = first - frame%ao_first(atom) + 1
         source_last = last - frame%ao_first(atom) + 1
         first = first - ao_range(1) + 1
         last = last - ao_range(1) + 1
         value(:, first:last) = gaussian(:, source_first:source_last)
         IF (gradients) gradient(:, first:last, :) = gaussian_gradient(:, source_first:source_last, :)
         IF (derivative) response(:, first:last) = gaussian_response(:, source_first:source_last)
         IF (PRESENT(gradient_response)) THEN
            gradient_response(:, first:last, :) = gaussian_gradient_response(:, source_first:source_last, :)
            DEALLOCATE (gaussian_gradient_response)
         END IF
         IF (derivative) DEALLOCATE (gaussian_response)
         IF (gradients) DEALLOCATE (gaussian_gradient)
         DEALLOCATE (gaussian)
      END DO
      ALLOCATE (delta(n, np))
      IF (gradients) ALLOCATE (delta_gradient(n, np, 3))
      IF (derivative) ALLOCATE (delta_response(n, np))
      IF (PRESENT(gradient_response)) ALLOCATE (delta_gradient_response(n, np, 3))
      DO q = 1, SIZE(coefficient_images, 2)
         CALL lattice_delta_fields(frame, points, ao_image - coefficient_images(:, q), delta, &
                                   delta_gradient, delta_response, delta_gradient_response, &
                                   point_response, atom_response, cell_response)
         CALL subtract_product(delta, coefficients(:, :, q), value)
         IF (derivative) THEN
            CALL subtract_product(delta_response, coefficients(:, :, q), response)
            CALL subtract_product(delta, coefficients_response(:, :, q), response)
         END IF
         DO axis = 1, 3
            IF (gradients) CALL subtract_product(delta_gradient(:, :, axis), coefficients(:, :, q), gradient(:, :, axis))
            IF (PRESENT(gradient_response)) THEN
               CALL subtract_product(delta_gradient_response(:, :, axis), coefficients(:, :, q), &
                                     gradient_response(:, :, axis))
               CALL subtract_product(delta_gradient(:, :, axis), coefficients_response(:, :, q), &
                                     gradient_response(:, :, axis))
            END IF
         END DO
      END DO
   END SUBROUTINE qs_upf_inverse_lattice_fields

! **************************************************************************************************
!> \brief Accumulate a field product, retaining independent value and response contributions.
!>        An exact zero left field requires no matrix product or temporary result.
! **************************************************************************************************
   SUBROUTINE subtract_product(left, right, result)
      COMPLEX(dp), INTENT(IN)                            :: left(:, :), right(:, :)
      COMPLEX(dp), INTENT(INOUT)                         :: result(:, :)

      INTEGER                                            :: k, m, n

      EXTERNAL zgemm

      m = SIZE(left, 1)
      k = SIZE(left, 2)
      n = SIZE(right, 2)
      CPASSERT(SIZE(right, 1) == k .AND. ALL(SHAPE(result) == [m, n]))
      IF (m == 0 .OR. k == 0 .OR. n == 0) RETURN
      IF (ALL(left == CMPLX(0.0_dp, 0.0_dp, dp))) RETURN
      CALL zgemm('N', 'N', m, n, k, CMPLX(-1.0_dp, 0.0_dp, dp), left, m, right, k, &
                 CMPLX(1.0_dp, 0.0_dp, dp), result, m)
   END SUBROUTINE subtract_product

! **************************************************************************************************
!> \brief Transpose the localized value and gradient maps to their coupled W_Q coefficients.
!>        Covectors use Re sum(conjg(G)*dfield). Point weights and any owner reductions
!>        belong to the caller. No phase, occupation or fitting weight is introduced here.
!> \param frame Original radial difference fields and geometry.
!> \param points ...
!> \param ao_image ...
!> \param coefficient_images ...
!> \param covector Value covector, indexed (point,requested AO).
!> \param coefficients_covector New local derivative, indexed (projector,requested AO,image).
!> \param gradient_covector Optional Cartesian gradient covector, axis last.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_lattice_field_adjoint(frame, points, ao_image, coefficient_images, &
                                                   covector, coefficients_covector, gradient_covector)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      REAL(dp), INTENT(IN)                               :: points(:, :)
      INTEGER, INTENT(IN)                                :: ao_image(3), coefficient_images(:, :)
      COMPLEX(dp), INTENT(IN)                            :: covector(:, :)
      COMPLEX(dp), INTENT(OUT)                           :: coefficients_covector(:, :, :)
      COMPLEX(dp), INTENT(IN), OPTIONAL                  :: gradient_covector(:, :, :)

      COMPLEX(dp), ALLOCATABLE                           :: delta(:, :), gradient(:, :, :)
      INTEGER                                            :: axis, n, np, q

      CALL check_lattice_images(frame, coefficient_images)
      CPASSERT(ALL(ao_image*frame%periodicity == ao_image))
      CPASSERT(SIZE(points, 1) == 3 .AND. SIZE(covector, 1) == SIZE(points, 2))
      n = SIZE(points, 2)
      np = SIZE(frame%coefficients, 1)
      CPASSERT(SIZE(coefficients_covector, 1) == np)
      CPASSERT(SIZE(coefficients_covector, 2) == SIZE(covector, 2))
      CPASSERT(SIZE(coefficients_covector, 3) == SIZE(coefficient_images, 2))
      ALLOCATE (delta(n, np))
      IF (PRESENT(gradient_covector)) THEN
         CPASSERT(ALL(SHAPE(gradient_covector) == [n, SIZE(covector, 2), 3]))
         ALLOCATE (gradient(n, np, 3))
      END IF
      DO q = 1, SIZE(coefficient_images, 2)
         CALL lattice_delta_fields(frame, points, ao_image - coefficient_images(:, q), delta, gradient)
         coefficients_covector(:, :, q) = -MATMUL(CONJG(TRANSPOSE(delta)), covector)
         DO axis = 1, 3
            IF (.NOT. PRESENT(gradient_covector)) CYCLE
            coefficients_covector(:, :, q) = coefficients_covector(:, :, q) - &
                                             MATMUL(CONJG(TRANSPOSE(gradient(:, :, axis))), gradient_covector(:, :, axis))
         END DO
      END DO
   END SUBROUTINE qs_upf_inverse_lattice_field_adjoint

! **************************************************************************************************
!> \brief Return the explicit geometry covectors of a localized field at fixed coupled W.
!>        The field is gaussian_factor*chi_R-sum_Q Delta_(R-Q)*W_Q, with default factor one.
!>        The value covector contains the consumer's measure. Point, native-center and lattice
!>        outputs are owner-local real derivatives. Cell vectors vary at fixed Cartesian points
!>        and native centers. All physical translations retain their actual integer indices.
!>        The caller adds moving-point geometry and the defining inverse-coefficient response.
!> \param frame ...
!> \param points ...
!> \param ao_image ...
!> \param coefficient_images ...
!> \param coefficients Complete projector rows and selected AO columns of W_Q.
!> \param ao_range ...
!> \param covector Weighted complex value derivative, in (point,selected AO) order.
!> \param point_geometry New Cartesian point covector, in (axis,point) order.
!> \param atom_geometry New native Cartesian center covector, in (axis,atom) order.
!> \param cell_geometry New lattice covector, in (axis,cell vector) order.
!> \param gaussian_factor Optional fixed Gaussian multiplier, two for soft-plus-Gaussian fields.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_lattice_geometry_adjoint(frame, points, ao_image, coefficient_images, coefficients, &
                                                      ao_range, covector, point_geometry, atom_geometry, cell_geometry, &
                                                      gaussian_factor)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      REAL(dp), INTENT(IN)                               :: points(:, :)
      INTEGER, INTENT(IN)                                :: ao_image(3), coefficient_images(:, :)
      COMPLEX(dp), INTENT(IN)                            :: coefficients(:, :, :)
      INTEGER, INTENT(IN)                                :: ao_range(2)
      COMPLEX(dp), INTENT(IN)                            :: covector(:, :)
      REAL(dp), INTENT(OUT)                              :: point_geometry(:, :), &
                                                            atom_geometry(:, :), &
                                                            cell_geometry(3, 3)
      REAL(dp), INTENT(IN), OPTIONAL                     :: gaussian_factor

      COMPLEX(dp), ALLOCATABLE                           :: delta(:, :), delta_covector(:, :), &
                                                            delta_gradient(:, :, :), &
                                                            gaussian(:, :), &
                                                            gaussian_gradient(:, :, :)
      INTEGER                                            :: atom, axis, center, first, image(3), j, &
                                                            last, n, nao, np, q, &
                                                            single_image(3, 1), source_first, &
                                                            source_last
      REAL(dp)                                           :: factor, total
      REAL(dp), ALLOCATABLE                              :: point_component(:)

      CALL check_lattice_images(frame, coefficient_images)
      CPASSERT(ALL(ao_image*frame%periodicity == ao_image))
      CPASSERT(SIZE(points, 1) == 3)
      n = SIZE(points, 2)
      np = SIZE(frame%coefficients, 1)
      nao = ao_range(2) - ao_range(1) + 1
      CPASSERT(nao > 0 .AND. ao_range(1) > 0)
      CPASSERT(ao_range(2) <= SIZE(frame%gaussian, 2))
      CPASSERT(ALL(SHAPE(covector) == [n, nao]))
      IF (ANY(SHAPE(coefficients) /= [np, nao, SIZE(coefficient_images, 2)])) THEN
         CPABORT("Invalid lattice geometry coefficient dimensions")
      END IF
      CPASSERT(ALL(SHAPE(point_geometry) == SHAPE(points)))
      CPASSERT(ALL(SHAPE(atom_geometry) == SHAPE(frame%positions)))
      point_geometry = 0.0_dp
      atom_geometry = 0.0_dp
      cell_geometry = 0.0_dp
      factor = 1.0_dp
      IF (PRESENT(gaussian_factor)) factor = gaussian_factor
      ALLOCATE (point_component(n))
      single_image(:, 1) = ao_image
      DO atom = 1, SIZE(frame%bases)
         first = MAX(frame%ao_first(atom), ao_range(1))
         last = MIN(frame%ao_last(atom), ao_range(2))
         IF (last < first) CYCLE
         ALLOCATE (gaussian(n, frame%bases(atom)%gto_basis_set%nsgf), &
                   gaussian_gradient(n, frame%bases(atom)%gto_basis_set%nsgf, 3))
         CALL gaussian_basis_fields_periodic(frame%bases(atom)%gto_basis_set, points, frame%positions(:, atom), &
                                             frame%cell, single_image, [0.0_dp, 0.0_dp, 0.0_dp], gaussian, gaussian_gradient)
         source_first = first - frame%ao_first(atom) + 1
         source_last = last - frame%ao_first(atom) + 1
         first = first - ao_range(1) + 1
         last = last - ao_range(1) + 1
         DO axis = 1, 3
            point_component(:) = factor*REAL(SUM(CONJG(covector(:, first:last))* &
                                                 gaussian_gradient(:, source_first:source_last, axis), DIM=2), dp)
            point_geometry(axis, :) = point_geometry(axis, :) + point_component
            total = SUM(point_component)
            atom_geometry(axis, atom) = atom_geometry(axis, atom) - total
            cell_geometry(axis, :) = cell_geometry(axis, :) - total*REAL(ao_image, dp)
         END DO
         DEALLOCATE (gaussian, gaussian_gradient)
      END DO
      ALLOCATE (delta(n, np), delta_gradient(n, np, 3), delta_covector(n, np))
      DO q = 1, SIZE(coefficient_images, 2)
         image = ao_image - coefficient_images(:, q)
         CALL lattice_delta_fields(frame, points, image, delta, delta_gradient)
         delta_covector(:, :) = -MATMUL(covector, CONJG(TRANSPOSE(coefficients(:, :, q))))
         DO center = 1, SIZE(frame%center_atoms)
            atom = frame%center_atoms(center)
            first = frame%projector_first(center)
            last = frame%projector_last(center)
            DO axis = 1, 3
               point_component(:) = REAL(SUM(CONJG(delta_covector(:, first:last))* &
                                             delta_gradient(:, first:last, axis), DIM=2), dp)
               point_geometry(axis, :) = point_geometry(axis, :) + point_component
               total = SUM(point_component)
               atom_geometry(axis, atom) = atom_geometry(axis, atom) - total
               DO j = 1, 3
                  cell_geometry(axis, j) = cell_geometry(axis, j) - total*REAL(image(j), dp)
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE qs_upf_inverse_lattice_geometry_adjoint

! **************************************************************************************************
!> \brief Evaluate all original Delta columns at one unphased physical image.
!> \param frame ...
!> \param points ...
!> \param image ...
!> \param value ...
!> \param gradient ...
!> \param response ...
!> \param gradient_response ...
!> \param point_response ...
!> \param atom_response ...
!> \param cell_response ...
! **************************************************************************************************
   SUBROUTINE lattice_delta_fields(frame, points, image, value, gradient, response, gradient_response, &
                                   point_response, atom_response, cell_response)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      REAL(dp), INTENT(IN)                               :: points(:, :)
      INTEGER, INTENT(IN)                                :: image(3)
      COMPLEX(dp), INTENT(OUT)                           :: value(:, :)
      COMPLEX(dp), INTENT(OUT), OPTIONAL                 :: gradient(:, :, :), response(:, :), &
                                                            gradient_response(:, :, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: point_response(:, :), &
                                                            atom_response(:, :), &
                                                            cell_response(3, 3)

      INTEGER                                            :: atom, center, first, kind, last, point, &
                                                            single_image(3, 1)
      LOGICAL                                            :: intersects
      REAL(dp)                                           :: position(3), radius_squared

      single_image(:, 1) = image
      value = CMPLX(0.0_dp, 0.0_dp, dp)
      IF (PRESENT(gradient)) gradient = CMPLX(0.0_dp, 0.0_dp, dp)
      IF (PRESENT(response)) response = CMPLX(0.0_dp, 0.0_dp, dp)
      IF (PRESENT(gradient_response)) gradient_response = CMPLX(0.0_dp, 0.0_dp, dp)
      DO center = 1, SIZE(frame%center_atoms)
         atom = frame%center_atoms(center)
         kind = frame%center_kinds(center)
         first = frame%projector_first(center)
         last = frame%projector_last(center)
         ! Use the original radial evaluator's exact support test before allocating its work arrays.
         ! Values and geometry responses outside that same support are identically zero.
         position = frame%positions(:, atom) + MATMUL(frame%cell, REAL(image, dp))
         radius_squared = frame%delta_cache(kind)%radius**2
         intersects = .FALSE.
         DO point = 1, SIZE(points, 2)
            IF (.NOT. (SUM((points(:, point) - position)**2) >= radius_squared)) THEN
               intersects = .TRUE.
               EXIT
            END IF
         END DO
         IF (.NOT. intersects) CYCLE
         IF (PRESENT(response)) THEN
            IF (PRESENT(gradient_response)) THEN
               CALL upf_radial_fields_periodic(frame%delta_cache(kind), points, frame%positions(:, atom), frame%cell, &
                                               single_image, [0.0_dp, 0.0_dp, 0.0_dp], value(:, first:last), &
                                               gradient(:, first:last, :), response(:, first:last), &
                                               gradient_response(:, first:last, :), point_response, &
                                               atom_response(:, atom), cell_response)
            ELSE IF (PRESENT(gradient)) THEN
               CALL upf_radial_fields_periodic(frame%delta_cache(kind), points, frame%positions(:, atom), frame%cell, &
                                               single_image, [0.0_dp, 0.0_dp, 0.0_dp], value(:, first:last), &
                                               gradient=gradient(:, first:last, :), response=response(:, first:last), &
                                               points_response=point_response, center_response=atom_response(:, atom), &
                                               cell_response=cell_response)
            ELSE
               CALL upf_radial_fields_periodic(frame%delta_cache(kind), points, frame%positions(:, atom), frame%cell, &
                                               single_image, [0.0_dp, 0.0_dp, 0.0_dp], value(:, first:last), &
                                               response=response(:, first:last), points_response=point_response, &
                                               center_response=atom_response(:, atom), cell_response=cell_response)
            END IF
         ELSE
            IF (PRESENT(gradient)) THEN
               CALL upf_radial_fields_periodic(frame%delta_cache(kind), points, frame%positions(:, atom), frame%cell, &
                                               single_image, [0.0_dp, 0.0_dp, 0.0_dp], value(:, first:last), &
                                               gradient=gradient(:, first:last, :))
            ELSE
               CALL upf_radial_fields_periodic(frame%delta_cache(kind), points, frame%positions(:, atom), frame%cell, &
                                               single_image, [0.0_dp, 0.0_dp, 0.0_dp], value(:, first:last))
            END IF
         END IF
      END DO
   END SUBROUTINE lattice_delta_fields

! **************************************************************************************************
!> \brief Require distinct physical coefficient images in the original periodic subspace.
!> \param frame ...
!> \param images ...
! **************************************************************************************************
   SUBROUTINE check_lattice_images(frame, images)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      INTEGER, INTENT(IN)                                :: images(:, :)

      INTEGER                                            :: i, j

      CPASSERT(.NOT. frame%atom_local)
      CPASSERT(SIZE(images, 1) == 3 .AND. SIZE(images, 2) > 0)
      DO i = 1, SIZE(images, 2)
         CPASSERT(ALL(images(:, i)*frame%periodicity == images(:, i)))
         DO j = 1, i - 1
            CPASSERT(ANY(images(:, i) /= images(:, j)))
         END DO
      END DO
   END SUBROUTINE check_lattice_images

! **************************************************************************************************
!> \brief Resolve the finite-image defining polynomial and its joint geometry response.
!>        Point fields are streamed. Full original projector matrices and all AO columns are
!>        retained. Matrices are reduced over the original point communicator after accumulation.
!> \param frame Original coupled inverse frame and its unchanged point measure and field images.
!> \param kernel_images New distinct image differences, shared by C_Q and B_Q.
!> \param coupling New original C_Q, without the identity.
!> \param overlaps New original B_Q.
!> \param coupling_response Optional new dC_Q.
!> \param overlaps_response Optional new dB_Q, required together with dC_Q and all geometry inputs.
!> \param point_response Optional physical point displacement.
!> \param atom_response Optional physical atom displacement.
!> \param cell_response Optional cell displacement at fixed reduced wave vector.
!> \param weights_response Optional response of the original point measure.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_lattice_kernel(frame, kernel_images, coupling, overlaps, &
                                            coupling_response, overlaps_response, point_response, &
                                            atom_response, cell_response, weights_response)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: kernel_images(:, :)
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: coupling(:, :, :), overlaps(:, :, :)
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT), OPTIONAL    :: coupling_response(:, :, :), &
                                                            overlaps_response(:, :, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: point_response(:, :), atom_response(:, :), &
                                                            cell_response(3, 3), weights_response(:)

      COMPLEX(dp), ALLOCATABLE :: beta(:, :), delta(:, :), gaussian(:, :), left(:, :), pair_c(:, :), pair_b(:, :), &
                                  dbeta(:, :), ddelta(:, :), dgaussian(:, :), dleft(:, :)
      INTEGER :: i, image(3), image_count, j, n, nao, np, nq, q
      INTEGER, ALLOCATABLE :: differences(:, :), pair_index(:, :)
      LOGICAL :: derivative

      CPASSERT(.NOT. frame%atom_local)
      derivative = PRESENT(coupling_response)
      CPASSERT(derivative .EQV. PRESENT(overlaps_response))
      IF (derivative) THEN
         CPASSERT(PRESENT(point_response))
         CPASSERT(PRESENT(atom_response))
         CPASSERT(PRESENT(cell_response))
         CPASSERT(PRESENT(weights_response))
         CPASSERT(ALL(SHAPE(point_response) == SHAPE(frame%points)))
         CPASSERT(ALL(SHAPE(atom_response) == SHAPE(frame%positions)))
         CPASSERT(SIZE(weights_response) == SIZE(frame%weights))
      END IF
      np = SIZE(frame%beta, 2)
      nao = SIZE(frame%gaussian, 2)
      n = SIZE(frame%weights)
      image_count = SIZE(frame%images, 2)
      ALLOCATE (differences(3, image_count**2), pair_index(image_count, image_count))
      nq = 0
      DO i = 1, image_count
         DO j = 1, image_count
            image = frame%images(:, j) - frame%images(:, i)
            q = 1
            DO WHILE (q <= nq)
               IF (ALL(differences(:, q) == image)) EXIT
               q = q + 1
            END DO
            IF (q > nq) THEN
               nq = nq + 1
               differences(:, nq) = image
            END IF
            pair_index(j, i) = q
         END DO
      END DO
      ALLOCATE (kernel_images, SOURCE=differences(:, :nq))
      ALLOCATE (coupling(np, np, nq), overlaps(np, nao, nq), SOURCE=CMPLX(0.0_dp, 0.0_dp, dp))
      ALLOCATE (beta(n, np), left(n, np), delta(n, np), gaussian(n, nao), pair_c(np, np), pair_b(np, nao))
      IF (derivative) THEN
         ALLOCATE (coupling_response(np, np, nq), overlaps_response(np, nao, nq), SOURCE=CMPLX(0.0_dp, 0.0_dp, dp))
         ALLOCATE (dbeta(n, np), dleft(n, np), ddelta(n, np), dgaussian(n, nao))
      END IF
      DO i = 1, image_count
         CALL original_fields(frame%images(:, i), projector_only=.TRUE.)
         left(:, :) = beta
         IF (derivative) dleft(:, :) = dbeta
         ! Compact original projectors outside their support contribute exactly zero.
         ! A zero value with a nonzero supplied response must still be retained.
         IF (derivative) THEN
            IF (.NOT. ANY(left /= 0.0_dp) .AND. .NOT. ANY(dleft /= 0.0_dp)) CYCLE
         ELSE
            IF (.NOT. ANY(left /= 0.0_dp)) CYCLE
         END IF
         DO j = 1, image_count
            CALL original_fields(frame%images(:, j), projector_only=.FALSE.)
            q = pair_index(j, i)
            CALL upf_inverse_pair(left, delta, frame%weights, pair_c)
            CALL upf_inverse_pair(left, gaussian, frame%weights, pair_b)
            coupling(:, :, q) = coupling(:, :, q) + pair_c
            overlaps(:, :, q) = overlaps(:, :, q) + pair_b
            IF (derivative) THEN
               CALL upf_inverse_pair_response(left, delta, frame%weights, dleft, ddelta, weights_response, pair_c)
               CALL upf_inverse_pair_response(left, gaussian, frame%weights, dleft, dgaussian, weights_response, pair_b)
               coupling_response(:, :, q) = coupling_response(:, :, q) + pair_c
               overlaps_response(:, :, q) = overlaps_response(:, :, q) + pair_b
            END IF
         END DO
      END DO
      DO q = 1, nq
         CALL frame%point_comm%sum(coupling(:, :, q))
         CALL frame%point_comm%sum(overlaps(:, :, q))
         IF (derivative) THEN
            CALL frame%point_comm%sum(coupling_response(:, :, q))
            CALL frame%point_comm%sum(overlaps_response(:, :, q))
         END IF
      END DO

   CONTAINS

! **************************************************************************************************
!> \brief Evaluate original unphased fields at one physical image, using native column order.
!> \param physical_image ...
!> \param projector_only Evaluate the left projectors or the right difference/Gaussian fields.
! **************************************************************************************************
      SUBROUTINE original_fields(physical_image, projector_only)
      INTEGER, INTENT(IN)                                :: physical_image(3)
      LOGICAL, INTENT(IN)                                :: projector_only

      INTEGER                                            :: atom, center, first, kind, last, &
                                                            single_image(3, 1)

         single_image(:, 1) = physical_image
         DO atom = 1, SIZE(frame%bases)
            IF (projector_only) EXIT
            first = frame%ao_first(atom)
            last = frame%ao_last(atom)
            IF (derivative) THEN
               CALL gaussian_basis_fields_periodic(frame%bases(atom)%gto_basis_set, frame%points, &
                            frame%positions(:, atom), frame%cell, single_image, [0.0_dp, 0.0_dp, 0.0_dp], gaussian(:, first:last), &
                                                   response=dgaussian(:, first:last), point_response=point_response, &
                                                   center_response=atom_response(:, atom), cell_response=cell_response)
            ELSE
               CALL gaussian_basis_fields_periodic(frame%bases(atom)%gto_basis_set, frame%points, &
                              frame%positions(:, atom), frame%cell, single_image, [0.0_dp, 0.0_dp, 0.0_dp], gaussian(:, first:last))
            END IF
         END DO
         DO center = 1, SIZE(frame%center_atoms)
            atom = frame%center_atoms(center)
            kind = frame%center_kinds(center)
            first = frame%projector_first(center)
            last = frame%projector_last(center)
            IF (projector_only .AND. derivative) THEN
               CALL upf_radial_fields_periodic(frame%beta_cache(kind), frame%points, frame%positions(:, atom), &
                                               frame%cell, single_image, [0.0_dp, 0.0_dp, 0.0_dp], beta(:, first:last), &
                                               response=dbeta(:, first:last), points_response=point_response, &
                                               center_response=atom_response(:, atom), cell_response=cell_response)
            ELSE IF (derivative) THEN
               CALL upf_radial_fields_periodic(frame%delta_cache(kind), frame%points, frame%positions(:, atom), &
                                               frame%cell, single_image, [0.0_dp, 0.0_dp, 0.0_dp], delta(:, first:last), &
                                               response=ddelta(:, first:last), points_response=point_response, &
                                               center_response=atom_response(:, atom), cell_response=cell_response)
            ELSE IF (projector_only) THEN
               CALL upf_radial_fields_periodic(frame%beta_cache(kind), frame%points, frame%positions(:, atom), &
                                               frame%cell, single_image, [0.0_dp, 0.0_dp, 0.0_dp], beta(:, first:last))
            ELSE
               CALL upf_radial_fields_periodic(frame%delta_cache(kind), frame%points, frame%positions(:, atom), &
                                               frame%cell, single_image, [0.0_dp, 0.0_dp, 0.0_dp], delta(:, first:last))
            END IF
         END DO
      END SUBROUTINE original_fields
   END SUBROUTINE qs_upf_inverse_lattice_kernel

! **************************************************************************************************
!> \brief Transpose the original physical-image defining integrals into geometry and measure.
!>        Uses Re sum(conjg(G_C)*dC+conjg(G_B)*dB), with complete projector and AO matrices.
!>        Both original fields and their actual image indices are differentiated. Point and
!>        weight outputs belong to local quadrature rows. All center and cell outputs are
!>        point-owner local too. The caller reduces their complete energy contraction once.
!>        Cell derivatives hold Cartesian points and native centers fixed. Moving-point and
!>        integration-measure responses are returned separately, without division by weights.
!> \param frame Original fields, distributed quadrature and physical image domain.
!> \param kernel_images Distinct physical differences including every original image pair.
!> \param coupling_covector Complete original C_Q energy derivative in kernel-image order.
!> \param overlaps_covector Complete original B_Q energy derivative in the same order.
!> \param point_geometry New owner-local Cartesian point derivative.
!> \param atom_geometry New owner-local native Cartesian center derivative.
!> \param cell_geometry New owner-local lattice derivative at fixed Cartesian coordinates.
!> \param weights_covector New owner-local derivative of each original quadrature weight.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_lattice_kernel_adjoint(frame, kernel_images, coupling_covector, overlaps_covector, &
                                                    point_geometry, atom_geometry, cell_geometry, weights_covector)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      INTEGER, INTENT(IN)                                :: kernel_images(:, :)
      COMPLEX(dp), INTENT(IN)                            :: coupling_covector(:, :, :), &
                                                            overlaps_covector(:, :, :)
      REAL(dp), INTENT(OUT)                              :: point_geometry(:, :), atom_geometry(:, :), &
                                                            cell_geometry(3, 3), weights_covector(:)

      COMPLEX(dp), ALLOCATABLE :: beta(:, :), beta_gradient(:, :, :), delta(:, :), delta_gradient(:, :, :), &
                                  gaussian(:, :), gaussian_gradient(:, :, :), g_beta(:, :), g_delta(:, :), g_gaussian(:, :)
      INTEGER                                            :: atom, center, first, i, image(3), j, &
                                                            kind, last, n, nao, np, point, q, &
                                                            single_image(3, 1)

      CALL check_lattice_images(frame, kernel_images)
      np = SIZE(frame%beta, 2)
      nao = SIZE(frame%gaussian, 2)
      n = SIZE(frame%weights)
      IF (ANY(SHAPE(coupling_covector) /= [np, np, SIZE(kernel_images, 2)]) .OR. &
          ANY(SHAPE(overlaps_covector) /= [np, nao, SIZE(kernel_images, 2)])) THEN
         CPABORT("Invalid defining lattice covector dimensions")
      END IF
      CPASSERT(ALL(SHAPE(point_geometry) == SHAPE(frame%points)))
      CPASSERT(ALL(SHAPE(atom_geometry) == SHAPE(frame%positions)))
      CPASSERT(SIZE(weights_covector) == n)
      point_geometry = 0.0_dp
      atom_geometry = 0.0_dp
      cell_geometry = 0.0_dp
      weights_covector = 0.0_dp
      ALLOCATE (beta(n, np), beta_gradient(n, np, 3), delta(n, np), delta_gradient(n, np, 3), &
                gaussian(n, nao), gaussian_gradient(n, nao, 3), &
                g_beta(n, np), g_delta(n, np), g_gaussian(n, nao))
      DO i = 1, SIZE(frame%images, 2)
         single_image(:, 1) = frame%images(:, i)
         DO center = 1, SIZE(frame%center_atoms)
            atom = frame%center_atoms(center)
            kind = frame%center_kinds(center)
            first = frame%projector_first(center)
            last = frame%projector_last(center)
            CALL upf_radial_fields_periodic(frame%beta_cache(kind), frame%points, frame%positions(:, atom), &
                                            frame%cell, single_image, [0.0_dp, 0.0_dp, 0.0_dp], beta(:, first:last), &
                                            gradient=beta_gradient(:, first:last, :))
         END DO
         ! Outside every compact projector support, both the defining integrals and
         ! all geometry/measure derivatives vanish. Values alone do not suffice
         ! here, since a zero field can still have a nonzero spatial derivative.
         IF (ALL(beta == CMPLX(0.0_dp, 0.0_dp, dp)) .AND. &
             ALL(beta_gradient == CMPLX(0.0_dp, 0.0_dp, dp))) CYCLE
         DO j = 1, SIZE(frame%images, 2)
            image = frame%images(:, j) - frame%images(:, i)
            q = 1
            DO WHILE (q <= SIZE(kernel_images, 2))
               IF (ALL(kernel_images(:, q) == image)) EXIT
               q = q + 1
            END DO
            IF (q > SIZE(kernel_images, 2)) THEN
               CPABORT("Missing original defining lattice image difference")
            END IF
            CALL lattice_delta_fields(frame, frame%points, frame%images(:, j), delta, delta_gradient)
            single_image(:, 1) = frame%images(:, j)
            DO atom = 1, SIZE(frame%bases)
               first = frame%ao_first(atom)
               last = frame%ao_last(atom)
               CALL gaussian_basis_fields_periodic(frame%bases(atom)%gto_basis_set, frame%points, &
                                                   frame%positions(:, atom), frame%cell, single_image, [0.0_dp, 0.0_dp, 0.0_dp], &
                                                   gaussian(:, first:last), gaussian_gradient(:, first:last, :))
            END DO
            g_delta(:, :) = MATMUL(beta, coupling_covector(:, :, q))
            g_gaussian(:, :) = MATMUL(beta, overlaps_covector(:, :, q))
            weights_covector(:) = weights_covector + REAL(SUM(CONJG(g_delta)*delta, DIM=2) + &
                                                          SUM(CONJG(g_gaussian)*gaussian, DIM=2), dp)
            g_beta(:, :) = MATMUL(delta, CONJG(TRANSPOSE(coupling_covector(:, :, q)))) + &
                           MATMUL(gaussian, CONJG(TRANSPOSE(overlaps_covector(:, :, q))))
            DO point = 1, n
               g_beta(point, :) = frame%weights(point)*g_beta(point, :)
               g_delta(point, :) = frame%weights(point)*g_delta(point, :)
               g_gaussian(point, :) = frame%weights(point)*g_gaussian(point, :)
            END DO
            DO center = 1, SIZE(frame%center_atoms)
               atom = frame%center_atoms(center)
               first = frame%projector_first(center)
               last = frame%projector_last(center)
               CALL field_geometry(g_beta(:, first:last), beta_gradient(:, first:last, :), atom, frame%images(:, i))
               CALL field_geometry(g_delta(:, first:last), delta_gradient(:, first:last, :), atom, frame%images(:, j))
            END DO
            DO atom = 1, SIZE(frame%bases)
               first = frame%ao_first(atom)
               last = frame%ao_last(atom)
               CALL field_geometry(g_gaussian(:, first:last), gaussian_gradient(:, first:last, :), atom, frame%images(:, j))
            END DO
         END DO
      END DO

   CONTAINS

! **************************************************************************************************
!> \brief Contract one weighted field covector with its original Cartesian gradients.
!> \param covector ...
!> \param gradient ...
!> \param atom ...
!> \param image ...
! **************************************************************************************************
      SUBROUTINE field_geometry(covector, gradient, atom, image)
      COMPLEX(dp), INTENT(IN)                            :: covector(:, :), gradient(:, :, :)
      INTEGER, INTENT(IN)                                :: atom, image(3)

      INTEGER                                            :: axis
      REAL(dp)                                           :: component(n), total

         DO axis = 1, 3
            component(:) = REAL(SUM(CONJG(covector)*gradient(:, :, axis), DIM=2), dp)
            point_geometry(axis, :) = point_geometry(axis, :) + component
            total = SUM(component)
            atom_geometry(axis, atom) = atom_geometry(axis, atom) - total
            cell_geometry(axis, :) = cell_geometry(axis, :) - total*REAL(image, dp)
         END DO
      END SUBROUTINE field_geometry
   END SUBROUTINE qs_upf_inverse_lattice_kernel_adjoint
END MODULE qs_upf_inverse_lattice
