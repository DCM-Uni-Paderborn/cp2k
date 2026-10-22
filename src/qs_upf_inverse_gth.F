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
   USE libgrpp,                         ONLY: libgrpp_init
   USE particle_types,                  ONLY: particle_type
   USE paw_core_kinetic,                ONLY: paw_core_kinetic_build_ae
   USE paw_types,                       ONLY: qs_upf_inverse_gth_core_type => paw_core_projector_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_upf_inverse,                  ONLY: qs_upf_inverse_atomic_overlap
   USE qs_upf_inverse_types,            ONLY: qs_upf_inverse_type
   USE upf_inverse,                     ONLY: upf_inverse_delta_init
   USE upf_projector_integrals,         ONLY: upf_projector_angular_init,&
                                              upf_projector_type
   USE upf_types,                       ONLY: atom_upfpot_type,&
                                              upf_radial_weights
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: qs_upf_inverse_gth_nonlocal
   PUBLIC :: qs_upf_inverse_gth_core_type, qs_upf_inverse_gth_core_init, qs_upf_inverse_gth_core_nonlocal

CONTAINS

! **************************************************************************************************
!> \brief Add G^dagger D_G G with G=<p_G|chi>-<p_G|Delta>W and W=K^-1 B.
!>        Gaussian/GTH overlaps are analytic. GTH/Delta overlaps use the original UPF radial
!>        quadrature with analytic Gaussian angular integrals. Differentiate both overlaps and W.
!>        Supplied PAW projections add <p_G|Delta_PAW><beta_PAW|U chi> before forming the operator.
!> \param frame Defining geometry, original AO ordering and reconstruction coefficients.
!> \param kinds Original native GTH and UPF datasets with borrowed live orbital bases.
!> \param particles Original atom-to-kind mapping, in the frame's atom order.
!> \param OPERATOR Optional existing dense operator, incremented in place.
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
! **************************************************************************************************
!> \brief Keep signed original AE core orbitals and their spin-summed magnetic-state occupations.
!>        Reuse the original-data closed-core assignment and charge checks. The angular overlap
!>        service integrates these orbitals on their full original radial mesh without a fit,
!>        normalization change, or additional tail truncation. No model-core kinetic field is needed.
!> \param pot Original scalar PAW data with complete compatible signed occupied core waves.
!> \param degree Largest Cartesian Gaussian degree, including requested geometry derivatives.
!> \param core Owned original-data cache, including an explicit unavailable-data reason.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_gth_core_init(pot, degree, core)
      TYPE(atom_upfpot_type), INTENT(IN)                 :: pot
      INTEGER, INTENT(IN)                                :: degree
      TYPE(qs_upf_inverse_gth_core_type), INTENT(OUT)    :: core

      INTEGER                                            :: i, j, nr, nshell, p
      INTEGER, ALLOCATABLE                               :: angular(:), selected(:), shells(:)
      REAL(dp), ALLOCATABLE                              :: tau_ae(:), weights(:)

      core%initialized = .TRUE.
      CALL paw_core_kinetic_build_ae(pot, tau_ae, core%available, core%reason, selected)
      IF (.NOT. core%available) RETURN
      DEALLOCATE (tau_ae)
      nshell = COUNT(selected > 0)
      IF (nshell == 0) RETURN
      nr = SIZE(pot%r)
      ALLOCATE (angular(nshell), shells(nshell))
      j = 0
      DO i = 1, SIZE(selected)
         IF (selected(i) == 0) CYCLE
         j = j + 1
         shells(j) = i
         angular(j) = pot%paw%gipaw_core_orbitals%angular_momenta(i)
      END DO
      CALL upf_projector_angular_init(angular, degree, core%orbitals)
      CALL libgrpp_init()
      ALLOCATE (core%orbitals%r, SOURCE=pot%r)
      core%orbitals%radius = MAXVAL(pot%r)
      ALLOCATE (core%orbitals%weighted_beta(nr, nshell), core%occupations(core%orbitals%nproj))
      weights = upf_radial_weights(pot)*pot%r
      DO j = 1, nshell
         i = shells(j)
         ! The original file stores u(r). The radial overlap measure for R=u/r is w*r*u.
         core%orbitals%weighted_beta(:, j) = weights*pot%paw%gipaw_core_orbitals%wfc(:, i)
      END DO
      DO p = 1, core%orbitals%nproj
         j = core%orbitals%radial_index(p)
         core%occupations(p) = REAL(selected(shells(j)), dp)/REAL(2*angular(j) + 1, dp)
      END DO
   END SUBROUTINE qs_upf_inverse_gth_core_init

! **************************************************************************************************
!> \brief Add the original frozen PAW core expectation of the GTH nonlocal one-body operator.
!>        Each physical projector image contributes its own quadratic expectation. Images are not
!>        combined coherently, and no valence occupation, spin, k-point weight, or inverse map enters.
!>        This constant has no valence Hamiltonian contribution. Its geometry response moves both
!>        atomic centers and the image lattice, with the original radial measure held fixed.
!> \param frame Original positions, cell, and complete fixed physical projector image list.
!> \param kinds Original scalar GTH and PAW datasets.
!> \param particles Original atom-to-kind mapping in the frame's atom order.
!> \param cores Prepared original PAW core caches, indexed by kind.
!> \param energy Existing energy, incremented once independently of valence spin and k-point loops.
!> \param response Optional existing joint geometry response, incremented in place.
!> \param atom_response Required with response, motion of every physical atom.
!> \param cell_response Required with response, motion of the image lattice vectors.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_gth_core_nonlocal(frame, kinds, particles, cores, energy, response, &
                                               atom_response, cell_response)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      TYPE(qs_kind_type), INTENT(IN)                     :: kinds(:)
      TYPE(particle_type), INTENT(IN)                    :: particles(:)
      TYPE(qs_upf_inverse_gth_core_type), INTENT(IN)     :: cores(:)
      REAL(dp), INTENT(INOUT)                            :: energy
      REAL(dp), INTENT(INOUT), OPTIONAL                  :: response
      REAL(dp), INTENT(IN), OPTIONAL                     :: atom_response(:, :), cell_response(3, 3)

      COMPLEX(dp), ALLOCATABLE                           :: doverlap(:, :), overlap(:, :), &
                                                            weighted(:, :)
      INTEGER                                            :: atom, core_atom, core_kind, image, kind, &
                                                            ncore, nf
      TYPE(gth_projector_fields_type), ALLOCATABLE, &
         TARGET                                          :: projectors(:)
      TYPE(gto_basis_set_type)                           :: native_basis

      CPASSERT(SIZE(cores) == SIZE(kinds))
      CPASSERT(ALL(SHAPE(frame%positions) == [3, SIZE(particles)]))
      CPASSERT(SIZE(frame%images, 1) == 3)
      IF (PRESENT(response)) THEN
         CPASSERT(PRESENT(atom_response))
         CPASSERT(ALL(SHAPE(atom_response) == [3, SIZE(particles)]))
         CPASSERT(PRESENT(cell_response))
      END IF
      ALLOCATE (projectors(SIZE(kinds)))
      DO kind = 1, SIZE(kinds)
         IF (kinds(kind)%ghost) CYCLE
         IF (.NOT. ASSOCIATED(kinds(kind)%gth_potential)) CYCLE
         IF (kinds(kind)%gth_potential%nppnl == 0) CYCLE
         CALL gth_projector_fields_init(kinds(kind)%gth_potential, projectors(kind))
      END DO
      DO atom = 1, SIZE(particles)
         CPASSERT(ASSOCIATED(particles(atom)%atomic_kind))
         kind = particles(atom)%atomic_kind%kind_number
         CPASSERT(kind >= 1 .AND. kind <= SIZE(kinds))
         nf = projectors(kind)%nfunction
         IF (nf == 0) CYCLE
         CALL gth_projector_basis(projectors(kind), native_basis)
         DO core_atom = 1, SIZE(particles)
            CPASSERT(ASSOCIATED(particles(core_atom)%atomic_kind))
            core_kind = particles(core_atom)%atomic_kind%kind_number
            CPASSERT(core_kind >= 1 .AND. core_kind <= SIZE(kinds))
            IF (kinds(core_kind)%ghost) CYCLE
            IF (.NOT. ASSOCIATED(kinds(core_kind)%sgp_potential)) CYCLE
            IF (.NOT. ALLOCATED(kinds(core_kind)%sgp_potential%upf_data)) CYCLE
            IF (.NOT. kinds(core_kind)%sgp_potential%upf_data%is_paw) CYCLE
            CPASSERT(cores(core_kind)%initialized)
            IF (.NOT. cores(core_kind)%available) THEN
               CPABORT('Original PAW core orbitals are unavailable: '//TRIM(cores(core_kind)%reason))
            END IF
            ncore = cores(core_kind)%orbitals%nproj
            IF (ncore == 0) CYCLE
            ALLOCATE (overlap(ncore, nf), weighted(ncore, nf))
            IF (PRESENT(response)) ALLOCATE (doverlap(ncore, nf))
            DO image = 1, SIZE(frame%images, 2)
               IF (PRESENT(response)) THEN
                  CALL qs_upf_inverse_atomic_overlap(native_basis, cores(core_kind)%orbitals, &
                                                     frame%positions(:, core_atom), frame%positions(:, atom), frame%cell, &
                                                     frame%images(:, image:image), [0.0_dp, 0.0_dp, 0.0_dp], overlap, doverlap, &
                                                     atom_response(:, core_atom), atom_response(:, atom), cell_response)
               ELSE
                  CALL qs_upf_inverse_atomic_overlap(native_basis, cores(core_kind)%orbitals, &
                                                     frame%positions(:, core_atom), frame%positions(:, atom), frame%cell, &
                                                     frame%images(:, image:image), [0.0_dp, 0.0_dp, 0.0_dp], overlap)
               END IF
               weighted = MATMUL(overlap, projectors(kind)%coupling)
               energy = energy + SUM(cores(core_kind)%occupations*REAL(SUM(CONJG(overlap)*weighted, DIM=2), dp))
               IF (PRESENT(response)) THEN
                  response = response + 2.0_dp*SUM(cores(core_kind)%occupations* &
                                                   REAL(SUM(CONJG(doverlap)*weighted, DIM=2), dp))
               END IF
            END DO
            DEALLOCATE (overlap, weighted)
            IF (ALLOCATED(doverlap)) DEALLOCATE (doverlap)
         END DO
      END DO
   END SUBROUTINE qs_upf_inverse_gth_core_nonlocal
END MODULE qs_upf_inverse_gth
