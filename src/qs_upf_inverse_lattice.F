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
   USE kinds,                           ONLY: dp
   USE qs_upf_inverse_types,            ONLY: qs_upf_inverse_type
   USE upf_inverse_quadrature,          ONLY: upf_inverse_pair,&
                                              upf_inverse_pair_response
   USE upf_radial_fields,               ONLY: upf_radial_fields_periodic
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: qs_upf_inverse_lattice_kernel
   PUBLIC :: qs_upf_inverse_lattice_fields, qs_upf_inverse_lattice_field_adjoint

CONTAINS

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
         value(:, :) = value - MATMUL(delta, coefficients(:, :, q))
         IF (derivative) THEN
            response(:, :) = response - MATMUL(delta_response, coefficients(:, :, q)) - &
                             MATMUL(delta, coefficients_response(:, :, q))
         END IF
         DO axis = 1, 3
            IF (gradients) gradient(:, :, axis) = gradient(:, :, axis) - &
                                                  MATMUL(delta_gradient(:, :, axis), coefficients(:, :, q))
            IF (PRESENT(gradient_response)) THEN
               gradient_response(:, :, axis) = gradient_response(:, :, axis) - &
                                               MATMUL(delta_gradient_response(:, :, axis), coefficients(:, :, q)) - &
                                               MATMUL(delta_gradient(:, :, axis), coefficients_response(:, :, q))
            END IF
         END DO
      END DO
   END SUBROUTINE qs_upf_inverse_lattice_fields

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

      INTEGER                                            :: atom, center, first, kind, last, &
                                                            single_image(3, 1)

      single_image(:, 1) = image
      DO center = 1, SIZE(frame%center_atoms)
         atom = frame%center_atoms(center)
         kind = frame%center_kinds(center)
         first = frame%projector_first(center)
         last = frame%projector_last(center)
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
         CALL original_fields(frame%images(:, i))
         left(:, :) = beta
         IF (derivative) dleft(:, :) = dbeta
         DO j = 1, image_count
            CALL original_fields(frame%images(:, j))
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
! **************************************************************************************************
      SUBROUTINE original_fields(physical_image)
      INTEGER, INTENT(IN)                                :: physical_image(3)

      INTEGER                                            :: atom, center, first, kind, last, &
                                                            single_image(3, 1)

         single_image(:, 1) = physical_image
         DO atom = 1, SIZE(frame%bases)
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
            IF (derivative) THEN
               CALL upf_radial_fields_periodic(frame%beta_cache(kind), frame%points, frame%positions(:, atom), &
                                               frame%cell, single_image, [0.0_dp, 0.0_dp, 0.0_dp], beta(:, first:last), &
                                               response=dbeta(:, first:last), points_response=point_response, &
                                               center_response=atom_response(:, atom), cell_response=cell_response)
               CALL upf_radial_fields_periodic(frame%delta_cache(kind), frame%points, frame%positions(:, atom), &
                                               frame%cell, single_image, [0.0_dp, 0.0_dp, 0.0_dp], delta(:, first:last), &
                                               response=ddelta(:, first:last), points_response=point_response, &
                                               center_response=atom_response(:, atom), cell_response=cell_response)
            ELSE
               CALL upf_radial_fields_periodic(frame%beta_cache(kind), frame%points, frame%positions(:, atom), &
                                               frame%cell, single_image, [0.0_dp, 0.0_dp, 0.0_dp], beta(:, first:last))
               CALL upf_radial_fields_periodic(frame%delta_cache(kind), frame%points, frame%positions(:, atom), &
                                               frame%cell, single_image, [0.0_dp, 0.0_dp, 0.0_dp], delta(:, first:last))
            END IF
         END DO
      END SUBROUTINE original_fields
   END SUBROUTINE qs_upf_inverse_lattice_kernel
END MODULE qs_upf_inverse_lattice
