!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Native GTH nonlocal operators acting on the complete inverse Gaussian basis.
!>        GTH contributes no reconstruction difference. Its projector overlaps nevertheless change
!>        through the reconstruction at other atoms. Keep these radial overlaps outside the AO span.
! **************************************************************************************************
MODULE qs_upf_inverse_gth
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE gth_projector_fields,            ONLY: gth_projector_basis,&
                                              gth_projector_fields_init,&
                                              gth_projector_fields_type,&
                                              gth_projector_overlap
   USE kinds,                           ONLY: dp
   USE particle_types,                  ONLY: particle_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_upf_inverse,                  ONLY: qs_upf_inverse_atomic_overlap
   USE qs_upf_inverse_types,            ONLY: qs_upf_inverse_type
   USE upf_inverse,                     ONLY: upf_inverse_delta_init
   USE upf_projector_integrals,         ONLY: upf_projector_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: qs_upf_inverse_gth_nonlocal

CONTAINS

! **************************************************************************************************
!> \brief Add G^dagger D_G G with G=<p_G|chi>-<p_G|Delta>W and W=K^-1 B.
!>        Gaussian/GTH overlaps are analytic. GTH/Delta overlaps use the original UPF radial
!>        quadrature with analytic Gaussian angular integrals. Differentiate both overlaps and W.
!>        Supplied PAW projections add <p_G|Delta_PAW><beta_PAW|U chi> before forming the operator.
!> \param frame Defining geometry, original AO ordering and reconstruction coefficients.
!> \param kinds Original native GTH and UPF datasets with borrowed live orbital bases.
!> \param particles Original atom-to-kind mapping, in the frame's atom order.
!> \param operator Optional existing dense operator, incremented in place.
!>        At least one of operator or response must be supplied.
!> \param response Optional existing operator response, incremented in place.
!> \param coefficients_response Required with response, the defining response of W.
!> \param atom_response Required with response, motion of all Gaussian and projector centers.
!> \param cell_response Required with response, motion of the lattice vectors.
!> \param reduced_k_response Required with response, motion of the reduced Bloch vector.
!> \param paw_projections Optional original PAW bras against the same soft orbitals, in frame row order.
!>        Adds the full original PAW valence reconstruction before the GTH projectors act.
!>        Other rows are ignored. Frozen-core expectation values remain a separate contribution.
!> \param paw_projections_response Required with response and paw_projections, their complete joint response.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_gth_nonlocal(frame, kinds, particles, OPERATOR, response, coefficients_response, &
                                          atom_response, cell_response, reduced_k_response, &
                                          paw_projections, paw_projections_response)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      TYPE(qs_kind_type), INTENT(IN)                     :: kinds(:)
      TYPE(particle_type), INTENT(IN)                    :: particles(:)
      COMPLEX(dp), INTENT(INOUT), OPTIONAL               :: operator(:, :), response(:, :)
      COMPLEX(dp), INTENT(IN), OPTIONAL                  :: coefficients_response(:, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: atom_response(:, :), &
                                                            cell_response(3, 3), &
                                                            reduced_k_response(3)
      COMPLEX(dp), INTENT(IN), OPTIONAL                  :: paw_projections(:, :), &
                                                            paw_projections_response(:, :)

      COMPLEX(dp), ALLOCATABLE :: coupling(:, :), dcoupling(:, :), doverlap(:, :), overlap(:, :), &
         projection_change(:, :), projection_response(:, :), weighted(:, :)
      INTEGER                                            :: active, atom, basis_atom, degree, first, &
                                                            ikind, kind, last, n, nf, p, q
      LOGICAL                                            :: physical_paw, reconstruct, skip
      TYPE(gth_projector_fields_type), ALLOCATABLE, &
         TARGET                                          :: projectors(:)
      TYPE(gto_basis_set_type)                           :: native_basis
      LOGICAL                                            :: active_kind(SIZE(kinds))
      TYPE(upf_projector_type), ALLOCATABLE              :: difference(:)

      CPASSERT(SIZE(frame%positions, 2) == SIZE(particles))
      n = SIZE(frame%coefficients, 2)
      p = SIZE(frame%coefficients, 1)
      CPASSERT(PRESENT(OPERATOR) .OR. PRESENT(response))
      IF (PRESENT(OPERATOR)) THEN
         CPASSERT(ALL(SHAPE(OPERATOR) == [n, n]))
      END IF
      IF (PRESENT(response)) THEN
         CPASSERT(ALL(SHAPE(response) == [n, n]))
         CPASSERT(PRESENT(coefficients_response))
         CPASSERT(SIZE(coefficients_response, 1) == p)
         CPASSERT(SIZE(coefficients_response, 2) == n)
         CPASSERT(PRESENT(atom_response))
         CPASSERT(SIZE(atom_response, 1) == 3)
         CPASSERT(SIZE(atom_response, 2) == SIZE(particles))
         CPASSERT(PRESENT(cell_response))
         CPASSERT(PRESENT(reduced_k_response))
      END IF
      IF (PRESENT(paw_projections)) THEN
         CPASSERT(ALL(SHAPE(paw_projections) == [p, n]))
         IF (PRESENT(response)) THEN
            CPASSERT(PRESENT(paw_projections_response))
            CPASSERT(ALL(SHAPE(paw_projections_response) == [p, n]))
         END IF
      ELSE
         CPASSERT(.NOT. PRESENT(paw_projections_response))
      END IF
      active_kind = .FALSE.
      degree = -1
      DO kind = 1, SIZE(kinds)
         IF (kinds(kind)%ghost) CYCLE
         IF (.NOT. ASSOCIATED(kinds(kind)%gth_potential)) CYCLE
         ASSOCIATE (potential => kinds(kind)%gth_potential)
            CPASSERT(.NOT. potential%soc)
            IF (potential%nppnl == 0) CYCLE
            active_kind(kind) = .TRUE.
            degree = MAX(degree, potential%lprj_ppnl_max + 1)
         END ASSOCIATE
      END DO
      IF (.NOT. ANY(active_kind)) RETURN
      ALLOCATE (projectors(SIZE(kinds)), difference(SIZE(kinds)))
      DO kind = 1, SIZE(kinds)
         IF (active_kind(kind)) CALL gth_projector_fields_init(kinds(kind)%gth_potential, projectors(kind))
      END DO
      DO active = 1, SIZE(frame%center_atoms)
         kind = frame%center_kinds(active)
         IF (difference(kind)%nproj > 0) CYCLE
         CPASSERT(ASSOCIATED(kinds(kind)%sgp_potential))
         CPASSERT(ALLOCATED(kinds(kind)%sgp_potential%upf_data))
         reconstruct = .NOT. frame%direct
         IF (PRESENT(paw_projections)) THEN
            IF (kinds(kind)%sgp_potential%upf_data%is_paw) reconstruct = .TRUE.
         END IF
         IF (.NOT. reconstruct) CYCLE
         IF (.NOT. frame%direct) THEN
            IF (frame%delta_cache(kind)%radius == 0.0_dp) CYCLE
         END IF
         CALL upf_inverse_delta_init(kinds(kind)%sgp_potential%upf_data, degree, difference(kind))
      END DO
      DO atom = 1, SIZE(particles)
         CPASSERT(ASSOCIATED(particles(atom)%atomic_kind))
         ikind = particles(atom)%atomic_kind%kind_number
         CPASSERT(ikind >= 1 .AND. ikind <= SIZE(kinds))
         IF (.NOT. active_kind(ikind)) CYCLE
         nf = projectors(ikind)%nfunction
         ALLOCATE (overlap(nf, n), weighted(nf, n))
         IF (PRESENT(response)) ALLOCATE (doverlap(nf, n))
         DO basis_atom = 1, SIZE(frame%bases)
            first = frame%ao_first(basis_atom)
            last = frame%ao_last(basis_atom)
            IF (PRESENT(response)) THEN
               CALL gth_projector_overlap(projectors(ikind), frame%bases(basis_atom)%gto_basis_set, &
                                          frame%positions(:, atom), frame%positions(:, basis_atom), frame%cell, &
                                          frame%images, frame%reduced_k, overlap(:, first:last), doverlap(:, first:last), &
                                          atom_response(:, atom), atom_response(:, basis_atom), cell_response, reduced_k_response)
            ELSE
               CALL gth_projector_overlap(projectors(ikind), frame%bases(basis_atom)%gto_basis_set, &
                                          frame%positions(:, atom), frame%positions(:, basis_atom), frame%cell, &
                                          frame%images, frame%reduced_k, overlap(:, first:last))
            END IF
         END DO
         CALL gth_projector_basis(projectors(ikind), native_basis)
         DO active = 1, SIZE(frame%center_atoms)
            kind = frame%center_kinds(active)
            IF (difference(kind)%nproj == 0) CYCLE
            basis_atom = frame%center_atoms(active)
            p = frame%projector_first(active)
            q = frame%projector_last(active)
            physical_paw = .FALSE.
            IF (PRESENT(paw_projections)) physical_paw = kinds(kind)%sgp_potential%upf_data%is_paw
            IF (physical_paw) THEN
               ALLOCATE (projection_change, SOURCE=paw_projections(p:q, :))
               IF (.NOT. frame%direct) projection_change = projection_change - frame%coefficients(p:q, :)
               skip = ALL(projection_change == CMPLX(0.0_dp, 0.0_dp, dp))
               IF (PRESENT(response)) THEN
                  ALLOCATE (projection_response, SOURCE=paw_projections_response(p:q, :))
                  IF (.NOT. frame%direct) projection_response = projection_response - coefficients_response(p:q, :)
                  skip = skip .AND. ALL(projection_response == CMPLX(0.0_dp, 0.0_dp, dp))
               END IF
               IF (skip) THEN
                  DEALLOCATE (projection_change)
                  IF (ALLOCATED(projection_response)) DEALLOCATE (projection_response)
                  CYCLE
               END IF
            END IF
            ! Reverse relative images before taking the adjoint, including for asymmetric lists.
            ALLOCATE (coupling(q - p + 1, nf))
            IF (PRESENT(response)) THEN
               ALLOCATE (dcoupling(q - p + 1, nf))
               CALL qs_upf_inverse_atomic_overlap(native_basis, difference(kind), frame%positions(:, basis_atom), &
                                                  frame%positions(:, atom), frame%cell, -frame%images, frame%reduced_k, &
                                                  coupling, dcoupling, atom_response(:, basis_atom), atom_response(:, atom), &
                                                  cell_response, reduced_k_response)
               IF (physical_paw) THEN
                  doverlap = doverlap + MATMUL(CONJG(TRANSPOSE(dcoupling)), projection_change) + &
                             MATMUL(CONJG(TRANSPOSE(coupling)), projection_response)
               ELSE
                  doverlap = doverlap - MATMUL(CONJG(TRANSPOSE(dcoupling)), frame%coefficients(p:q, :)) - &
                             MATMUL(CONJG(TRANSPOSE(coupling)), coefficients_response(p:q, :))
               END IF
               DEALLOCATE (dcoupling)
            ELSE
               CALL qs_upf_inverse_atomic_overlap(native_basis, difference(kind), frame%positions(:, basis_atom), &
                                                  frame%positions(:, atom), frame%cell, -frame%images, frame%reduced_k, coupling)
            END IF
            IF (physical_paw) THEN
               overlap = overlap + MATMUL(CONJG(TRANSPOSE(coupling)), projection_change)
               DEALLOCATE (projection_change)
               IF (ALLOCATED(projection_response)) DEALLOCATE (projection_response)
            ELSE
               overlap = overlap - MATMUL(CONJG(TRANSPOSE(coupling)), frame%coefficients(p:q, :))
            END IF
            DEALLOCATE (coupling)
         END DO
         weighted = MATMUL(projectors(ikind)%coupling, overlap)
         IF (PRESENT(OPERATOR)) OPERATOR = OPERATOR + MATMUL(CONJG(TRANSPOSE(overlap)), weighted)
         IF (PRESENT(response)) THEN
            response = response + MATMUL(CONJG(TRANSPOSE(doverlap)), weighted) + &
                       MATMUL(CONJG(TRANSPOSE(overlap)), MATMUL(projectors(ikind)%coupling, doverlap))
            DEALLOCATE (doverlap)
         END IF
         DEALLOCATE (overlap, weighted)
      END DO
   END SUBROUTINE qs_upf_inverse_gth_nonlocal
END MODULE qs_upf_inverse_gth
