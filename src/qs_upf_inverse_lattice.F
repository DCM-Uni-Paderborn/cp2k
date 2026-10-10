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

CONTAINS

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
