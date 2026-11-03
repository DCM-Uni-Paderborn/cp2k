!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
!> \brief Original-coordinate auxiliary operators and physical projection/charge adjoints.
MODULE qs_upf_inverse_aux_operator
   USE admm_types,                      ONLY: admm_type,&
                                              admm_upf_ace_search_type,&
                                              get_admm_env
   USE admm_utils,                      ONLY: admm_charge_scaled_operator,&
                                              admm_purified_operator
   USE bibliography,                    ONLY: Lin2016ACE,&
                                              cite_reference
   USE cell_types,                      ONLY: cell_type
   USE cp_cfm_types,                    ONLY: cp_cfm_get_submatrix
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_copy,&
                                              dbcsr_p_type,&
                                              dbcsr_set
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_submatrix,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE hfx_ace_affine,                  ONLY: hfx_ace_affine_build,&
                                              hfx_ace_affine_evaluate
   USE hfx_ri,                          ONLY: hfx_ri_update_forces,&
                                              hfx_ri_update_ks
   USE hfx_ri_kp,                       ONLY: hfx_ri_update_forces_kp,&
                                              hfx_ri_update_ks_kp
   USE hfx_types,                       ONLY: hfx_ace_release_state,&
                                              hfx_type
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE input_constants,                 ONLY: do_admm_exch_scaling_merlot,&
                                              do_admm_exch_scaling_none,&
                                              do_admm_purify_mo_diag,&
                                              do_admm_purify_none
   USE input_cp2k_hfx,                  ONLY: ri_pmat
   USE input_section_types,             ONLY: section_vals_type
   USE integral_library_types,          ONLY: compare_potential_types
   USE kinds,                           ONLY: dp
   USE kpoint_methods,                  ONLY: kpoint_star_density,&
                                              kpoint_star_operator,&
                                              kpoint_star_operator_defect,&
                                              kpoint_star_orbitals
   USE kpoint_transitional,             ONLY: get_1d_pointer
   USE kpoint_types,                    ONLY: kpoint_env_type,&
                                              kpoint_type
   USE mathlib,                         ONLY: invert_matrix,&
                                              invmat
   USE message_passing,                 ONLY: mp_para_env_type
   USE particle_types,                  ONLY: particle_type
   USE paw_types,                       ONLY: paw_potential_type
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type
   USE qs_mo_types,                     ONLY: mo_set_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_upf_inverse_aux_xc,           ONLY: qs_upf_inverse_aux_xc_rescale
   USE qs_upf_inverse_background,       ONLY: qs_upf_inverse_background_collocate,&
                                              qs_upf_inverse_background_geometry
   USE qs_upf_inverse_charge_metric,    ONLY: qs_upf_inverse_charge_metric_contract
   USE qs_upf_inverse_coulomb,          ONLY: qs_upf_inverse_coulomb_metric_response
   USE qs_upf_inverse_grid,             ONLY: qs_upf_inverse_grid_metric_dense,&
                                              qs_upf_inverse_matrix_gather
   USE qs_upf_inverse_kpoints,          ONLY: qs_upf_inverse_kpoints_density,&
                                              qs_upf_inverse_kpoints_original_operator
   USE qs_upf_inverse_paw,              ONLY: qs_upf_inverse_paw_core_prepare,&
                                              qs_upf_inverse_paw_evaluate,&
                                              qs_upf_inverse_paw_response
   USE qs_upf_inverse_types,            ONLY: qs_upf_inverse_grid_type,&
                                              qs_upf_inverse_scf_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: qs_upf_inverse_aux_projection, qs_upf_inverse_aux_operator_evaluate
   PUBLIC :: qs_upf_inverse_aux_metric_response
   PUBLIC :: qs_upf_inverse_aux_core_exchange, qs_upf_inverse_aux_response
   PUBLIC :: qs_upf_inverse_aux_evaluate
   PUBLIC :: qs_upf_inverse_aux_ace_build, qs_upf_inverse_aux_ace_evaluate
   PUBLIC :: qs_upf_inverse_aux_ace_search_begin, qs_upf_inverse_aux_ace_search_end

CONTAINS

! **************************************************************************************************
!> \brief Begin one immutable auxiliary ACE inner functional. Its first evaluation builds the reference.
!>        Outside this interval OT always evaluates complete original RI exchange.
!> \param qs_env Native environment whose auxiliary density and exchange owners are borrowed.
!> \param active Whether this search uses auxiliary ACE. Other methods remain unchanged.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_aux_ace_search_begin(qs_env, active)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(OUT)                               :: active

      INTEGER                                            :: contribution
      TYPE(admm_type), POINTER                           :: admm
      TYPE(dft_control_type), POINTER                    :: control
      TYPE(hfx_type), POINTER                            :: exchange(:, :)

      active = .FALSE.
      CALL get_qs_env(qs_env, dft_control=control, admm_env=admm, x_data=exchange)
      IF (.NOT. control%do_admm .OR. .NOT. ASSOCIATED(admm) .OR. .NOT. ASSOCIATED(exchange)) RETURN
      IF (.NOT. qs_env%scf_control%use_ot) RETURN
      IF (SIZE(exchange, 1) == 0 .OR. SIZE(exchange, 2) == 0) RETURN
      IF (.NOT. ANY(exchange(:, 1)%use_ace)) RETURN
      CPASSERT(.NOT. admm%upf_ace_search%active)
      DO contribution = 1, SIZE(exchange, 1)
         IF (.NOT. exchange(contribution, 1)%use_ace) CYCLE
         ! Honor the native final full-exchange phase across orbital-chart rebuilds.
         IF (exchange(contribution, 1)%ace%convergence_rebuild) RETURN
      END DO
      admm%upf_ace_search = admm_upf_ace_search_type()
      admm%upf_ace_search%active = .TRUE.
      active = .TRUE.
   END SUBROUTINE qs_upf_inverse_aux_ace_search_begin

! **************************************************************************************************
!> \brief End the inner functional and release its independently owned projectors.
!>        The caller must reevaluate the complete original functional before physical acceptance.
!> \param qs_env Same native environment used to begin the search.
!> \param finish_scf Enter the native full-exchange finishing phase after inner stationarity.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_aux_ace_search_end(qs_env, finish_scf)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN), OPTIONAL                      :: finish_scf

      INTEGER                                            :: contribution
      TYPE(admm_type), POINTER                           :: admm
      TYPE(hfx_type), POINTER                            :: exchange(:, :)

      CALL get_qs_env(qs_env, admm_env=admm, x_data=exchange)
      CPASSERT(ASSOCIATED(admm))
      CPASSERT(admm%upf_ace_search%active)
      admm%upf_ace_search = admm_upf_ace_search_type()
      IF (PRESENT(finish_scf)) THEN
         IF (finish_scf) THEN
            CPASSERT(ASSOCIATED(exchange))
            DO contribution = 1, SIZE(exchange, 1)
               IF (exchange(contribution, 1)%use_ace) exchange(contribution, 1)%ace%convergence_rebuild = .TRUE.
            END DO
         END IF
      END IF
   END SUBROUTINE qs_upf_inverse_aux_ace_search_end

! **************************************************************************************************
!> \brief Evaluate one fixed inner exchange functional for all energy, gradient and curvature trials.
!>        Copy only owned local affine projectors, never distributed FM pointer descriptors.
!>        If compression fails at the reference, retain complete RI for this entire search.
!>        Both current ADMMP charge factors and their affine offsets remain part of the energy.
!> \param qs_env Current native auxiliary density, projection and physical mesh weights.
!> \param exchange Original HF repetitions, unchanged throughout the search.
!> \param geometry_changed Prepare original RI geometry on the first evaluation only.
!> \param energy Weighted inner valence exchange energy.
!> \param operators Unweighted physical-k/spin auxiliary derivative.
!> \param spin_energies Weighted inner spin energies, including complete affine offsets.
! **************************************************************************************************
   SUBROUTINE evaluate_auxiliary_search(qs_env, exchange, geometry_changed, energy, operators, spin_energies)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(hfx_type), INTENT(INOUT)                      :: exchange(:, :)
      LOGICAL, INTENT(IN)                                :: geometry_changed
      REAL(dp), INTENT(OUT)                              :: energy
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: operators(:, :, :, :)
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: spin_energies(:)

      COMPLEX(dp), ALLOCATABLE                           :: matrix(:, :)
      INTEGER                                            :: contribution, ik, nk, ns, spin
      LOGICAL                                            :: compressible, ready
      REAL(dp)                                           :: fraction, part, scale, &
                                                            unused_strain(3, 3)
      REAL(dp), ALLOCATABLE                              :: unused_gradient(:, :)
      TYPE(admm_type), POINTER                           :: admm

      CALL get_qs_env(qs_env, admm_env=admm)
      CPASSERT(ASSOCIATED(admm))
      CPASSERT(admm%upf_density%valid)
      ns = SIZE(admm%upf_density%matrices, 4)
      nk = SIZE(admm%upf_density%matrices, 3)
      ASSOCIATE (search => admm%upf_ace_search)
         CPASSERT(search%active)
         IF (.NOT. search%initialized) THEN
            ALLOCATE (search%fractions(SIZE(exchange, 1)), search%potentials(SIZE(exchange, 1)))
            compressible = .TRUE.
            DO contribution = 1, SIZE(exchange, 1)
               fraction = exchange(contribution, 1)%general_parameter%fraction
               search%fractions(contribution) = fraction
               search%potentials(contribution) = exchange(contribution, 1)%potential_parameter
               IF (.NOT. MERGE(0.5_dp*fraction, fraction, ns == 1) > EPSILON(0.0_dp)) CYCLE
               compressible = compressible .AND. exchange(contribution, 1)%use_ace
            END DO
            ready = .FALSE.
            IF (compressible) THEN
               CALL qs_upf_inverse_aux_ace_build(qs_env, exchange, geometry_changed, energy, operators, spin_energies, ready)
            ELSE
               CALL auxiliary_ri(qs_env, exchange, geometry_changed, .FALSE., .FALSE., energy, operators, &
                                 unused_gradient, unused_strain, spin_energies, allow_ace=.TRUE.)
            END IF
            search%initialized = .TRUE.
            search%compressed = ready
            IF (.NOT. ready) RETURN
            ALLOCATE (search%references(nk, ns, SIZE(exchange, 1)))
            DO contribution = 1, SIZE(exchange, 1)
               IF (.NOT. ALLOCATED(exchange(contribution, 1)%ace%auxiliary)) CYCLE
               search%references(:, :, contribution) = exchange(contribution, 1)%ace%auxiliary
            END DO
         ELSE
            CPASSERT(.NOT. geometry_changed)
            CPASSERT(SIZE(search%fractions) == SIZE(exchange, 1))
            DO contribution = 1, SIZE(exchange, 1)
               ready = search%fractions(contribution) == exchange(contribution, 1)%general_parameter%fraction
               CPASSERT(ready)
               ready = compare_potential_types(search%potentials(contribution), exchange(contribution, 1)%potential_parameter)
               CPASSERT(ready)
            END DO
            IF (.NOT. search%compressed) THEN
               CALL auxiliary_ri(qs_env, exchange, .FALSE., .FALSE., .FALSE., energy, operators, &
                                 unused_gradient, unused_strain, spin_energies, allow_ace=.TRUE.)
               RETURN
            END IF
         END IF
         ! The very first inner operator is compressed too. Mixing a full reference operator
         ! with later frozen trial operators would give inconsistent orbital-chart curvature.
         IF (ALLOCATED(operators)) DEALLOCATE (operators)
         IF (ALLOCATED(spin_energies)) DEALLOCATE (spin_energies)
         ALLOCATE (operators, MOLD=admm%upf_density%matrices)
         ALLOCATE (spin_energies(ns), SOURCE=0.0_dp)
         operators = CMPLX(0.0_dp, 0.0_dp, dp)
         DO contribution = 1, SIZE(exchange, 1)
            fraction = search%fractions(contribution)
            IF (.NOT. MERGE(0.5_dp*fraction, fraction, ns == 1) > EPSILON(0.0_dp)) CYCLE
            DO spin = 1, ns
               scale = 1.0_dp
               IF (admm%do_admmp) THEN
                  ready = ieee_is_finite(admm%gsi(spin)) .AND. admm%gsi(spin) > 0.0_dp
                  CPASSERT(ready)
                  scale = admm%gsi(spin)**2
               END IF
               DO ik = 1, nk
                  CALL hfx_ace_affine_evaluate(search%references(ik, spin, contribution), &
                                               admm%upf_density%matrices(:, :, ik, spin), part, matrix)
                  operators(:, :, ik, spin) = operators(:, :, ik, spin) + scale*matrix
                  spin_energies(spin) = spin_energies(spin) + admm%upf_aux_metric%weights(ik)*scale*part
               END DO
            END DO
         END DO
      END ASSOCIATE
      energy = SUM(spin_energies)
   END SUBROUTINE evaluate_auxiliary_search

! **************************************************************************************************
!> \brief Build owned auxiliary ACE projectors from each complete native RI exchange contribution.
!>        Return the original RI energy and operator even if compression is unavailable. The
!>        occupied interpolation factors do not replace the actual mixed RI build density.
!>        Each kernel, physical k point and spin keeps its own affine energy and projector set.
!>        Compression must preserve physical-star and little-group operator covariance.
!>        A nonzero retained density without occupied factors requires full exchange.
!> \param qs_env Native primary MOs, independent auxiliary basis and installed auxiliary density.
!> \param exchange Actual HF repetitions. Projectors belong to the first thread of each repetition.
!> \param geometry_changed Recompute the original native RI integrals when the geometry changes.
!> \param energy Original valence exchange energy, including mesh weights and ADMMP rescaling.
!> \param operators Original unweighted physical-k/spin auxiliary valence exchange derivative.
!> \param spin_energies Original spin-resolved exchange energies with the same weights and scaling.
!> \param ready All active ACE contributions have usable projectors. False retains full RI outputs.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_aux_ace_build(qs_env, exchange, geometry_changed, energy, operators, spin_energies, ready)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(hfx_type), INTENT(INOUT)                      :: exchange(:, :)
      LOGICAL, INTENT(IN)                                :: geometry_changed
      REAL(dp), INTENT(OUT)                              :: energy
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: operators(:, :, :, :)
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: spin_energies(:)
      LOGICAL, INTENT(OUT)                               :: ready

      COMPLEX(dp), ALLOCATABLE                           :: compressed(:, :, :), &
                                                            factors(:, :, :, :), matrix(:, :), &
                                                            projection(:, :, :), &
                                                            term_operator(:, :, :, :)
      INTEGER                                            :: contribution, ik, info, nk, ns, rank, &
                                                            spin
      INTEGER, ALLOCATABLE                               :: ranks(:, :)
      LOGICAL                                            :: available, selected, term_ready
      REAL(dp) :: anchor_shift, covariance_defect, covariance_scale, fraction, reference_energy, &
         scale, term_energy, unused_strain(3, 3), weight
      REAL(dp), ALLOCATABLE                              :: term_spin_energy(:), &
                                                            unused_gradient(:, :)
      TYPE(admm_type), POINTER                           :: admm

      CALL get_qs_env(qs_env, admm_env=admm)
      CPASSERT(ASSOCIATED(admm))
      CPASSERT(admm%upf_density%valid)
      CPASSERT(SIZE(exchange, 2) > 0)
      ns = SIZE(admm%upf_density%matrices, 4)
      nk = SIZE(admm%upf_density%matrices, 3)
      CALL qs_upf_inverse_aux_projection(qs_env, projection, interpolation_factors=factors, &
                                         interpolation_ranks=ranks, interpolation_ready=available)
      ALLOCATE (operators, MOLD=admm%upf_density%matrices)
      ALLOCATE (spin_energies(ns))
      operators = CMPLX(0.0_dp, 0.0_dp, dp)
      spin_energies = 0.0_dp
      energy = 0.0_dp
      ready = available
      selected = .FALSE.
      DO contribution = 1, SIZE(exchange, 1)
         CALL hfx_ace_release_state(exchange(contribution, 1)%ace, keep_motion=.TRUE.)
         fraction = exchange(contribution, 1)%general_parameter%fraction
         CPASSERT(ieee_is_finite(fraction))
         IF (fraction < 0.0_dp) CPABORT('Auxiliary RI needs nonnegative HF fractions')
         IF (.NOT. MERGE(0.5_dp*fraction, fraction, ns == 1) > EPSILON(0.0_dp)) CYCLE
         CALL auxiliary_ri(qs_env, exchange(contribution:contribution, :), geometry_changed, .FALSE., .FALSE., &
                           term_energy, term_operator, unused_gradient, unused_strain, term_spin_energy, &
                           ace_reference=.TRUE.)
         IF (exchange(contribution, 1)%use_ace) THEN
            selected = .TRUE.
            term_ready = available
            IF (available) THEN
               ALLOCATE (exchange(contribution, 1)%ace%auxiliary(nk, ns))
               DO spin = 1, ns
                  DO ik = 1, nk
                     rank = ranks(ik, spin)
                     IF (rank == 0 .AND. ANY(admm%upf_density%matrices(:, :, ik, spin) /= 0.0_dp)) THEN
                        term_ready = .FALSE.
                        CYCLE
                     END IF
                     reference_energy = 0.5_dp*REAL(SUM(CONJG(term_operator(:, :, ik, spin))* &
                                                        admm%upf_density%matrices(:, :, ik, spin)), dp)
                     CALL hfx_ace_affine_build(term_operator(:, :, ik, spin), factors(:, :rank, ik, spin), &
                                               reference_energy, exchange(contribution, 1)%ace%auxiliary(ik, spin), &
                                               info, reference_density=admm%upf_density%matrices(:, :, ik, spin))
                     IF (info /= 0) term_ready = .FALSE.
                  END DO
               END DO
            END IF
            IF (term_ready .AND. ASSOCIATED(admm%upf_inverse%bloch%kpoints)) THEN
               ALLOCATE (compressed(SIZE(term_operator, 1), SIZE(term_operator, 2), nk))
               DO spin = 1, ns
                  DO ik = 1, nk
                     CALL hfx_ace_affine_evaluate(exchange(contribution, 1)%ace%auxiliary(ik, spin), &
                                                  admm%upf_density%matrices(:, :, ik, spin), reference_energy, matrix)
                     compressed(:, :, ik) = matrix
                  END DO
                  CALL kpoint_star_operator_defect(admm%upf_inverse%bloch%star, compressed, &
                                                   covariance_defect, covariance_scale, little_group=.TRUE.)
                  IF (.NOT. ieee_is_finite(covariance_defect) .OR. .NOT. ieee_is_finite(covariance_scale)) &
                     term_ready = .FALSE.
                  IF (covariance_defect > 1.e-10_dp*MAX(1.0_dp, covariance_scale)) term_ready = .FALSE.
               END DO
               DEALLOCATE (compressed)
            END IF
            ready = ready .AND. term_ready
         END IF
         IF (.NOT. admm%do_admmp) THEN
            ! Native RI and its half-trace can differ by reduction roundoff. Anchor the
            ! affine energy to the actual RI value without modifying its AO derivative.
            anchor_shift = term_energy - SUM(term_spin_energy)
            term_spin_energy(1) = term_spin_energy(1) + anchor_shift
            IF (exchange(contribution, 1)%use_ace .AND. available) THEN
               ik = MAXLOC(admm%upf_aux_metric%weights, DIM=1)
               weight = admm%upf_aux_metric%weights(ik)
               CPASSERT(weight > 0.0_dp)
               exchange(contribution, 1)%ace%auxiliary(ik, 1)%offset = &
                  exchange(contribution, 1)%ace%auxiliary(ik, 1)%offset + anchor_shift/weight
            END IF
         END IF
         DO spin = 1, ns
            scale = 1.0_dp
            IF (admm%do_admmp) THEN
               CPASSERT(ieee_is_finite(admm%gsi(spin)))
               CPASSERT(admm%gsi(spin) > 0.0_dp)
               scale = admm%gsi(spin)**2
            END IF
            operators(:, :, :, spin) = operators(:, :, :, spin) + scale*term_operator(:, :, :, spin)
            spin_energies(spin) = spin_energies(spin) + scale*term_spin_energy(spin)
         END DO
         energy = energy + term_energy
      END DO
      IF (admm%do_admmp) energy = SUM(spin_energies)
      ready = ready .AND. selected
      DO contribution = 1, SIZE(exchange, 1)
         IF (.NOT. ready) THEN
            CALL hfx_ace_release_state(exchange(contribution, 1)%ace, keep_motion=.TRUE.)
         ELSE IF (ALLOCATED(exchange(contribution, 1)%ace%auxiliary)) THEN
            exchange(contribution, 1)%ace%built = .TRUE.
            exchange(contribution, 1)%ace%step_counter = 1
            exchange(contribution, 1)%ace%auxiliary_fraction = exchange(contribution, 1)%general_parameter%fraction
            exchange(contribution, 1)%ace%auxiliary_potential = exchange(contribution, 1)%potential_parameter
         END IF
      END DO
   END SUBROUTINE qs_upf_inverse_aux_ace_build

! **************************************************************************************************
!> \brief Evaluate frozen auxiliary ACE matrices and affine energies at the actual current RI density.
!>        Caches store raw valence exchange. Apply the current ADMMP spin charge ratio squared
!>        to both the matrix and complete affine energy. Core exchange and XC remain separate.
!>        Geometry changes, final full-exchange iterations, incompatible owners and newly populated
!>        empty blocks require full exchange. Failure returns no partial operator or spin energy.
!> \param qs_env Current auxiliary density, physical mesh weights and spin charge ratios.
!> \param exchange Fixed native HF owners. This evaluation changes neither caches nor counters.
!> \param geometry_changed Whether the AO basis, lattice or physical k-point geometry changed.
!> \param energy Weighted affine valence exchange energy, zero when the cache is unavailable.
!> \param operators Complete unweighted physical-k/spin AO derivative, unallocated on failure.
!> \param spin_energies Weighted affine spin energies, unallocated on failure.
!> \param ready All positive exchange contributions have compatible usable caches.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_aux_ace_evaluate(qs_env, exchange, geometry_changed, energy, operators, spin_energies, ready)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(hfx_type), INTENT(IN)                         :: exchange(:, :)
      LOGICAL, INTENT(IN)                                :: geometry_changed
      REAL(dp), INTENT(OUT)                              :: energy
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: operators(:, :, :, :)
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: spin_energies(:)
      LOGICAL, INTENT(OUT)                               :: ready

      COMPLEX(dp), ALLOCATABLE                           :: matrix(:, :)
      INTEGER                                            :: contribution, ik, n, nk, ns, spin
      LOGICAL                                            :: selected
      REAL(dp)                                           :: fraction, part, scale
      TYPE(admm_type), POINTER                           :: admm

      energy = 0.0_dp
      ready = .FALSE.
      IF (geometry_changed) RETURN
      CALL get_qs_env(qs_env, admm_env=admm)
      CPASSERT(ASSOCIATED(admm))
      IF (.NOT. admm%upf_density%valid) RETURN
      CPASSERT(admm%upf_aux_metric%ready)
      CPASSERT(SIZE(exchange, 2) > 0)
      n = SIZE(admm%upf_density%matrices, 1)
      nk = SIZE(admm%upf_density%matrices, 3)
      ns = SIZE(admm%upf_density%matrices, 4)
      IF (admm%do_admmp) THEN
         IF (.NOT. ALL(ieee_is_finite(admm%gsi(1:ns)))) RETURN
         IF (ANY(admm%gsi(1:ns) <= 0.0_dp)) RETURN
      END IF
      selected = .FALSE.
      DO contribution = 1, SIZE(exchange, 1)
         ASSOCIATE (hf => exchange(contribution, 1), ace => exchange(contribution, 1)%ace)
            fraction = hf%general_parameter%fraction
            CPASSERT(ieee_is_finite(fraction))
            IF (fraction < 0.0_dp) CPABORT('Auxiliary RI needs nonnegative HF fractions')
            IF (.NOT. MERGE(0.5_dp*fraction, fraction, ns == 1) > EPSILON(0.0_dp)) CYCLE
            selected = .TRUE.
            IF (.NOT. hf%use_ace .OR. .NOT. hf%do_hfx_ri) RETURN
            IF (.NOT. ASSOCIATED(hf%ri_data)) RETURN
            IF (hf%ri_data%input_flavor /= ri_pmat .OR. hf%ri_data%orb_basis_type /= 'AUX_FIT') RETURN
            IF (.NOT. ace%built .OR. ace%convergence_rebuild) RETURN
            IF (ace%auxiliary_fraction /= fraction) RETURN
            IF (.NOT. compare_potential_types(ace%auxiliary_potential, hf%potential_parameter)) RETURN
            IF (.NOT. ALLOCATED(ace%auxiliary)) RETURN
            IF (ANY(SHAPE(ace%auxiliary) /= [nk, ns])) RETURN
            DO spin = 1, ns
               DO ik = 1, nk
                  IF (.NOT. ace%auxiliary(ik, spin)%ready) RETURN
                  IF (.NOT. ALLOCATED(ace%auxiliary(ik, spin)%projectors)) RETURN
                  IF (SIZE(ace%auxiliary(ik, spin)%projectors, 1) /= n) RETURN
                  IF (SIZE(ace%auxiliary(ik, spin)%projectors, 2) == 0 .AND. &
                      ANY(admm%upf_density%matrices(:, :, ik, spin) /= 0.0_dp)) RETURN
               END DO
            END DO
         END ASSOCIATE
      END DO
      IF (.NOT. selected) RETURN
      ALLOCATE (operators, MOLD=admm%upf_density%matrices)
      ALLOCATE (spin_energies(ns))
      operators = CMPLX(0.0_dp, 0.0_dp, dp)
      spin_energies = 0.0_dp
      DO contribution = 1, SIZE(exchange, 1)
         fraction = exchange(contribution, 1)%general_parameter%fraction
         IF (.NOT. MERGE(0.5_dp*fraction, fraction, ns == 1) > EPSILON(0.0_dp)) CYCLE
         DO spin = 1, ns
            scale = 1.0_dp
            IF (admm%do_admmp) scale = admm%gsi(spin)**2
            DO ik = 1, nk
               CALL hfx_ace_affine_evaluate(exchange(contribution, 1)%ace%auxiliary(ik, spin), &
                                            admm%upf_density%matrices(:, :, ik, spin), part, matrix)
               operators(:, :, ik, spin) = operators(:, :, ik, spin) + scale*matrix
               spin_energies(spin) = spin_energies(spin) + admm%upf_aux_metric%weights(ik)*scale*part
            END DO
         END DO
      END DO
      energy = SUM(spin_energies)
      ready = .TRUE.
   END SUBROUTINE qs_upf_inverse_aux_ace_evaluate

! **************************************************************************************************
!> \brief Evaluate the complete auxiliary exchange correction using the actual native RI owners.
!>        Each HF repetition contributes once, irrespective of its thread replicas. RI uses
!>        independent temporary output matrices and retains its normal density/integral history.
!>        Core HF, selected XC, projection and physical charge normalization use that same density.
!> \param qs_env Prepared native projection, auxiliary occupations and RI basis owners.
!> \param xc_section Selected signed auxiliary XC correction.
!> \param exchange Actual native HF repetitions and their initialized AUX_FIT RI containers.
!> \param geometry_changed Rebuild each RI geometry when required by the caller.
!> \param components RI valence HF, PAW core HF, grid XC and atomic XC energies, in that order.
!> \param operators Complete primary AO derivative in original physical k/spin order, without weights.
!> \param atom_response Optional total atomic displacement, including affine cell motion.
!> \param cell_response Optional lattice derivative at fixed reduced k.
!> \param response Optional complete directional derivative. Supply the three directional arguments together.
!> \param atom_gradient Optional complete Cartesian energy gradients, evaluated with one native RI force call.
!> \param strain_gradient Optional complete affine dE/dstrain at fixed reduced k.
!> \param just_energy Complete original exchange outside an explicit immutable inner search.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_aux_evaluate(qs_env, xc_section, exchange, geometry_changed, components, operators, &
                                          atom_response, cell_response, response, atom_gradient, strain_gradient, just_energy)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: xc_section
      TYPE(hfx_type), INTENT(INOUT)                      :: exchange(:, :)
      LOGICAL, INTENT(IN)                                :: geometry_changed
      REAL(dp), INTENT(OUT)                              :: components(4)
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: operators(:, :, :, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: atom_response(:, :), cell_response(3, 3)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: response, atom_gradient(:, :), &
                                                            strain_gradient(3, 3)
      LOGICAL, INTENT(IN), OPTIONAL                      :: just_energy

      COMPLEX(dp), ALLOCATABLE                           :: ri_operator(:, :, :, :)
      INTEGER                                            :: contribution, natom
      LOGICAL                                            :: active, frozen_search, full_exchange, &
                                                            geometry, ready, rebuild, use_virial
      REAL(dp)                                           :: energy, ri_energy, ri_strain(3, 3)
      REAL(dp), ALLOCATABLE                              :: ri_gradient(:, :), ri_spin_energy(:)
      TYPE(admm_type), POINTER                           :: admm

      CPASSERT(PRESENT(response) .EQV. PRESENT(atom_response))
      CPASSERT(PRESENT(response) .EQV. PRESENT(cell_response))
      geometry = PRESENT(atom_gradient) .OR. PRESENT(strain_gradient)
      CPASSERT(.NOT. (geometry .AND. PRESENT(response)))
      use_virial = PRESENT(strain_gradient)
      geometry = geometry .OR. PRESENT(response)
      IF (PRESENT(response)) use_virial = ANY(cell_response /= 0.0_dp)
      CPASSERT(SIZE(exchange, 2) > 0)
      active = ANY(exchange(:, 1)%use_ace)
      IF (active) CALL cite_reference(Lin2016ACE)
      CALL get_qs_env(qs_env, admm_env=admm)
      CPASSERT(ASSOCIATED(admm))
      frozen_search = admm%upf_ace_search%active
      IF (frozen_search) THEN
         CPASSERT(active)
      END IF
      ! Ordinary OT evaluations, canonical acceptance and geometry use the original functional.
      ! Only an explicit inner search may reuse one immutable affine reference.
      full_exchange = geometry .OR. (active .AND. qs_env%scf_control%use_ot)
      IF (PRESENT(just_energy)) full_exchange = full_exchange .OR. just_energy
      rebuild = .FALSE.
      DO contribution = 1, SIZE(exchange, 1)
         IF (.NOT. exchange(contribution, 1)%use_ace) CYCLE
         ASSOCIATE (ace => exchange(contribution, 1)%ace)
            IF (geometry_changed .AND. .NOT. frozen_search) CALL hfx_ace_release_state(ace, keep_motion=.TRUE.)
            full_exchange = full_exchange .OR. ace%convergence_rebuild .OR. &
                            (ace%dynamic_mode .AND. ace%geometry_step == 0)
            rebuild = rebuild .OR. .NOT. ace%built .OR. &
                      MOD(ace%step_counter, MAX(1, exchange(contribution, 1)%ace_rebuild_freq)) == 0
         END ASSOCIATE
      END DO
      IF (frozen_search) THEN
         CPASSERT(.NOT. geometry)
         CALL evaluate_auxiliary_search(qs_env, exchange, geometry_changed, ri_energy, ri_operator, ri_spin_energy)
         CALL get_qs_env(qs_env, natom=natom)
         ALLOCATE (ri_gradient(3, natom), SOURCE=0.0_dp)
         ri_strain = 0.0_dp
      ELSE IF (active .AND. .NOT. full_exchange) THEN
         ready = .FALSE.
         IF (.NOT. rebuild) THEN
            CALL qs_upf_inverse_aux_ace_evaluate(qs_env, exchange, geometry_changed, &
                                                 ri_energy, ri_operator, ri_spin_energy, ready)
         END IF
         IF (.NOT. ready) THEN
            CALL qs_upf_inverse_aux_ace_build(qs_env, exchange, geometry_changed, &
                                              ri_energy, ri_operator, ri_spin_energy, ready)
         ELSE
            DO contribution = 1, SIZE(exchange, 1)
               IF (exchange(contribution, 1)%use_ace .AND. exchange(contribution, 1)%ace%built) &
                  exchange(contribution, 1)%ace%step_counter = exchange(contribution, 1)%ace%step_counter + 1
            END DO
         END IF
         CALL get_qs_env(qs_env, natom=natom)
         ALLOCATE (ri_gradient(3, natom), SOURCE=0.0_dp)
         ri_strain = 0.0_dp
      ELSE
         CALL auxiliary_ri(qs_env, exchange, geometry_changed, geometry, use_virial, &
                           ri_energy, ri_operator, ri_gradient, ri_strain, ri_spin_energy, allow_ace=active)
         IF (active .AND. .NOT. geometry) THEN
            DO contribution = 1, SIZE(exchange, 1)
               IF (.NOT. exchange(contribution, 1)%use_ace) CYCLE
               ! Full RI outside an OT inner search supersedes its compressed reference.
               IF (.NOT. qs_env%scf_control%use_ot .AND. &
                   .NOT. exchange(contribution, 1)%ace%convergence_rebuild) CYCLE
               exchange(contribution, 1)%ace%built = .FALSE.
               exchange(contribution, 1)%ace%step_counter = 0
            END DO
         END IF
      END IF
      IF (PRESENT(response)) THEN
         CALL qs_upf_inverse_aux_response(qs_env, xc_section, exchange, ri_energy, ri_operator, &
                                          ri_gradient, ri_strain, atom_response, cell_response, &
                                          energy, operators, response, components, valence_spin_energy=ri_spin_energy)
      ELSE
         components(1) = ri_energy
         CALL qs_upf_inverse_aux_operator_evaluate(qs_env, xc_section, ri_operator, &
                                                   components(3), components(4), operators, &
                                                   exchange=exchange, core_energy=components(2), &
                                                   valence_spin_energy=ri_spin_energy)
         IF (geometry) THEN
            CALL auxiliary_gradients(qs_env, xc_section, exchange, ri_energy, ri_operator, ri_gradient, ri_strain, &
                                     atom_gradient, strain_gradient, ri_spin_energy)
         END IF
      END IF
      IF (active .AND. geometry) THEN
         DO contribution = 1, SIZE(exchange, 1)
            IF (exchange(contribution, 1)%use_ace) &
               CALL hfx_ace_release_state(exchange(contribution, 1)%ace, keep_motion=.TRUE.)
         END DO
      END IF
   END SUBROUTINE qs_upf_inverse_aux_evaluate

! **************************************************************************************************
!> \brief Contract complete auxiliary geometry while reusing the native RI energy, operator and gradients.
!>        The independent XC and metric responses retain the same frozen primary occupations.
!>        Cartesian variations keep the cell fixed. Affine strain moves atoms and cell together.
!> \param qs_env ...
!> \param xc_section ...
!> \param exchange ...
!> \param ri_energy ...
!> \param ri_operator ...
!> \param ri_gradient ...
!> \param ri_strain ...
!> \param atom_gradient ...
!> \param strain_gradient ...
!> \param ri_spin_energy Physical spin-resolved valence energies, including ADMMP eta**2 and k weights.
! **************************************************************************************************
   SUBROUTINE auxiliary_gradients(qs_env, xc_section, exchange, ri_energy, ri_operator, ri_gradient, ri_strain, &
                                  atom_gradient, strain_gradient, ri_spin_energy)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: xc_section
      TYPE(hfx_type), INTENT(IN)                         :: exchange(:, :)
      REAL(dp), INTENT(IN)                               :: ri_energy
      COMPLEX(dp), INTENT(IN)                            :: ri_operator(:, :, :, :)
      REAL(dp), INTENT(IN)                               :: ri_gradient(:, :), ri_strain(3, 3)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: atom_gradient(:, :), &
                                                            strain_gradient(3, 3)
      REAL(dp), INTENT(IN)                               :: ri_spin_energy(:)

      COMPLEX(dp), ALLOCATABLE                           :: operators(:, :, :, :), &
                                                            xc_operator(:, :, :, :)
      INTEGER                                            :: atom, axis, other
      REAL(dp)                                           :: atom_energy, energy, grid_energy, &
                                                            response, strain(3, 3)
      REAL(dp), ALLOCATABLE                              :: eta_derivative(:), motion(:, :), &
                                                            positions(:, :), xc_gradient(:, :), &
                                                            xc_strain(:, :)
      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_type), POINTER                       :: atoms(:)

      CALL get_qs_env(qs_env, cell=cell, particle_set=atoms)
      ALLOCATE (motion(3, SIZE(atoms)), positions(3, SIZE(atoms)))
      DO atom = 1, SIZE(atoms)
         positions(:, atom) = atoms(atom)%r
      END DO
      ALLOCATE (eta_derivative(SIZE(ri_operator, 4)))
      IF (PRESENT(atom_gradient)) ALLOCATE (xc_gradient, MOLD=atom_gradient)
      IF (PRESENT(strain_gradient)) ALLOCATE (xc_strain(3, 3))
      CALL qs_upf_inverse_aux_xc_rescale(qs_env, xc_section, grid_energy, atom_energy, xc_operator, eta_derivative, &
                                         atom_gradient=xc_gradient, strain_gradient=xc_strain)
      IF (PRESENT(atom_gradient)) THEN
         CPASSERT(ALL(SHAPE(atom_gradient) == SHAPE(motion)))
         DO atom = 1, SIZE(atoms)
            DO axis = 1, 3
               motion = 0.0_dp
               motion(axis, atom) = 1.0_dp
               CALL compose_auxiliary_response(qs_env, exchange, ri_energy, ri_operator, ri_gradient, ri_strain, &
                                               motion, 0.0_dp*cell%hmat, grid_energy, atom_energy, xc_operator, &
                                               eta_derivative, xc_gradient(axis, atom), energy, operators, response, &
                                               valence_spin_energy=ri_spin_energy)
               atom_gradient(axis, atom) = response
            END DO
         END DO
      END IF
      IF (PRESENT(strain_gradient)) THEN
         DO other = 1, 3
            DO axis = 1, 3
               strain = 0.0_dp
               strain(axis, other) = 1.0_dp
               motion(:, :) = MATMUL(strain, positions)
               CALL compose_auxiliary_response(qs_env, exchange, ri_energy, ri_operator, ri_gradient, ri_strain, &
                                               motion, MATMUL(strain, cell%hmat), grid_energy, atom_energy, xc_operator, &
                                               eta_derivative, xc_strain(axis, other), energy, operators, response, &
                                               valence_spin_energy=ri_spin_energy)
               strain_gradient(axis, other) = response
            END DO
         END DO
      END IF
   END SUBROUTINE auxiliary_gradients

! **************************************************************************************************
!> \brief Sum native RI repetitions into owned original-coordinate outputs.
!>        Clear the temporary KS matrices separately for each repetition. Gamma obtains energy
!>        from those matrices, while k-point RI also owns a physical Bloch operator accumulator.
!>        Preserve the caller's auxiliary accumulator, its native KS matrices and primary storage.
!> \param qs_env ...
!> \param exchange ...
!> \param geometry_changed ...
!> \param geometry ...
!> \param use_virial ...
!> \param energy Globally reduced valence-valence RI HF energy.
!> \param operators Unweighted physical-k/spin auxiliary AO derivative.
!> \param atom_gradient Positive globally reduced Cartesian energy gradients, zero without geometry.
!> \param strain_gradient Positive globally reduced affine dE/dstrain, zero without requested strain.
!> \param spin_energies Spin-resolved valence energies with the same k weights and ADMMP scaling.
!> \param ace_reference Optional full RI reference for explicit ACE construction, before ADMMP scaling.
!> \param allow_ace Permit complete original RI while ACE is selected, with normal scaling and geometry.
! **************************************************************************************************
   SUBROUTINE auxiliary_ri(qs_env, exchange, geometry_changed, geometry, use_virial, &
                           energy, operators, atom_gradient, strain_gradient, spin_energies, ace_reference, allow_ace)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(hfx_type), INTENT(IN)                         :: exchange(:, :)
      LOGICAL, INTENT(IN)                                :: geometry_changed, geometry, use_virial
      REAL(dp), INTENT(OUT)                              :: energy
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: operators(:, :, :, :)
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: atom_gradient(:, :)
      REAL(dp), INTENT(OUT)                              :: strain_gradient(3, 3)
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: spin_energies(:)
      LOGICAL, INTENT(IN), OPTIONAL                      :: ace_reference, allow_ace

      COMPLEX(dp), ALLOCATABLE                           :: gamma(:, :), saved_exchange(:, :, :, :)
      INTEGER                                            :: contribution, ik, image, natom, ns, spin
      LOGICAL                                            :: bloch, permitted, reference, rescale, &
                                                            valid
      REAL(dp)                                           :: fraction, native_energy, part, &
                                                            raw_energy, term_energy, &
                                                            term_strain(3, 3)
      REAL(dp), ALLOCATABLE                              :: spin_weights(:), term_gradient(:, :), &
                                                            term_spin_energy(:)
      TYPE(admm_type), POINTER                           :: admm
      TYPE(dbcsr_p_type), POINTER                        :: density(:, :), native_ks(:, :), &
                                                            term_ks(:, :)
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_upf_inverse_scf_type), POINTER             :: state

      NULLIFY (term_ks)
      CALL get_qs_env(qs_env, admm_env=admm, natom=natom)
      CPASSERT(ASSOCIATED(admm))
      CPASSERT(admm%upf_density%valid)
      state => admm%upf_inverse
      CPASSERT(ASSOCIATED(state))
      CPASSERT(state%basis_type == 'AUX_FIT')
      CPASSERT(SIZE(exchange, 2) > 0)
      ns = SIZE(admm%upf_density%matrices, 4)
      ALLOCATE (spin_energies(ns), term_spin_energy(ns))
      spin_energies = 0.0_dp
      reference = .FALSE.
      IF (PRESENT(ace_reference)) reference = ace_reference
      permitted = reference
      IF (PRESENT(allow_ace)) permitted = permitted .OR. allow_ace
      CPASSERT(.NOT. (reference .AND. geometry))
      rescale = admm%do_admmp .AND. .NOT. reference
      IF (rescale) THEN
         ALLOCATE (spin_weights(ns))
         CPASSERT(ALL(admm%gsi(1:ns) > 0.0_dp))
         spin_weights(:) = admm%gsi(1:ns)**2
      END IF
      ! An unallocated optional spin_weights argument is absent for the other variants.
      bloch = ASSOCIATED(state%bloch%kpoints)
      ! RI image densities include every physical star contribution. The returned
      ! image adjoint has no mesh weights; transfer_operator applies the matching star adjoint.
      DO contribution = 1, SIZE(exchange, 1)
         fraction = exchange(contribution, 1)%general_parameter%fraction
         CPASSERT(ieee_is_finite(fraction))
         IF (fraction < 0.0_dp) THEN
            CPABORT('Auxiliary RI needs nonnegative HF fractions')
         END IF
         IF (.NOT. MERGE(0.5_dp*fraction, fraction, ns == 1) > EPSILON(0.0_dp)) CYCLE
         IF (.NOT. exchange(contribution, 1)%do_hfx_ri .OR. (exchange(contribution, 1)%use_ace .AND. .NOT. permitted)) THEN
            CPABORT('Original-UPF auxiliary exchange requires uncompressed RI')
         END IF
         CPASSERT(ASSOCIATED(exchange(contribution, 1)%ri_data))
         ASSOCIATE (ri => exchange(contribution, 1)%ri_data)
            IF (ri%input_flavor /= ri_pmat .OR. ri%orb_basis_type /= 'AUX_FIT') THEN
               CPABORT('Original-UPF auxiliary exchange requires RI-RHO in AUX_FIT')
            END IF
            valid = compare_potential_types(ri%hfx_pot, exchange(contribution, 1)%potential_parameter)
            CPASSERT(valid)
         END ASSOCIATE
      END DO
      CALL get_admm_env(admm, matrix_ks_aux_fit_kp=native_ks, rho_aux_fit=rho)
      CALL qs_rho_get(rho, rho_ao_kp=density)
      CPASSERT(ASSOCIATED(native_ks) .AND. ASSOCIATED(density))
      valid = ALL(SHAPE(native_ks) == SHAPE(density)) .AND. SIZE(native_ks, 1) == ns
      CPASSERT(valid)
      ALLOCATE (operators, MOLD=admm%upf_density%matrices)
      operators = CMPLX(0.0_dp, 0.0_dp, dp)
      ALLOCATE (atom_gradient(3, natom), term_gradient(3, natom))
      atom_gradient = 0.0_dp
      strain_gradient = 0.0_dp
      energy = 0.0_dp
      CALL dbcsr_allocate_matrix_set(term_ks, ns, SIZE(native_ks, 2))
      DO image = 1, SIZE(native_ks, 2)
         DO spin = 1, ns
            ALLOCATE (term_ks(spin, image)%matrix)
            CALL dbcsr_copy(term_ks(spin, image)%matrix, native_ks(spin, image)%matrix)
         END DO
      END DO
      IF (bloch) CALL MOVE_ALLOC(state%exchange_k, saved_exchange)
      DO contribution = 1, SIZE(exchange, 1)
         fraction = exchange(contribution, 1)%general_parameter%fraction
         IF (.NOT. MERGE(0.5_dp*fraction, fraction, ns == 1) > EPSILON(0.0_dp)) CYCLE
         term_spin_energy = 0.0_dp
         DO image = 1, SIZE(term_ks, 2)
            DO spin = 1, ns
               CALL dbcsr_set(term_ks(spin, image)%matrix, 0.0_dp)
            END DO
         END DO
         IF (bloch) THEN
            IF (ALLOCATED(state%exchange_k)) state%exchange_k = CMPLX(0.0_dp, 0.0_dp, dp)
            CALL hfx_ri_update_ks_kp(qs_env, exchange(contribution, 1)%ri_data, term_ks, term_energy, &
                                     density, geometry_changed, ns, fraction)
            native_energy = term_energy
            CPASSERT(ALLOCATED(state%exchange_k))
            valid = ALL(SHAPE(state%exchange_k) == SHAPE(operators))
            CPASSERT(valid)
            DO spin = 1, ns
               DO ik = 1, SIZE(operators, 3)
                  ! Screened native contractions can retain a small anti-Hermitian numerical residue.
                  ! Hermitian occupations couple only to the Hermitian energy covector.
                  state%exchange_k(:, :, ik, spin) = 0.5_dp*(state%exchange_k(:, :, ik, spin) + &
                                                             CONJG(TRANSPOSE(state%exchange_k(:, :, ik, spin))))
                  term_spin_energy(spin) = term_spin_energy(spin) + 0.5_dp*admm%upf_aux_metric%weights(ik)* &
                                    REAL(SUM(CONJG(state%exchange_k(:, :, ik, spin))*admm%upf_density%matrices(:, :, ik, spin)), dp)
               END DO
            END DO
            raw_energy = SUM(term_spin_energy)
            valid = ABS(raw_energy - native_energy) <= 1.e-9_dp*MAX(1.0_dp, ABS(native_energy))
            IF (.NOT. valid) THEN
               CPABORT('Auxiliary RI energy/operator mismatch')
            END IF
            IF (rescale) THEN
               term_spin_energy(:) = spin_weights*term_spin_energy
               term_energy = SUM(term_spin_energy)
               DO spin = 1, ns
                  operators(:, :, :, spin) = operators(:, :, :, spin) + spin_weights(spin)*state%exchange_k(:, :, :, spin)
               END DO
            END IF
            IF (.NOT. rescale) operators(:, :, :, :) = operators + state%exchange_k
            IF (geometry) THEN
               CALL hfx_ri_update_forces_kp(qs_env, exchange(contribution, 1)%ri_data, ns, fraction, density, &
                                            use_virial=use_virial, atom_gradient=term_gradient, strain_gradient=term_strain, &
                                            spin_weights=spin_weights)
            END IF
         ELSE
            CALL hfx_ri_update_ks(qs_env, exchange(contribution, 1)%ri_data, term_ks, term_energy, &
                                  rho_ao=density, geometry_did_change=geometry_changed, nspins=ns, hf_fraction=fraction)
            raw_energy = term_energy
            IF (rescale) term_energy = 0.0_dp
            DO spin = 1, ns
               CALL qs_upf_inverse_matrix_gather(state%frame%ao, term_ks(spin, 1)%matrix, gamma)
               gamma(:, :) = 0.5_dp*(gamma + CONJG(TRANSPOSE(gamma)))
               part = 0.5_dp*REAL(SUM(CONJG(gamma)*admm%upf_density%matrices(:, :, 1, spin)), dp)
               term_spin_energy(spin) = part
               IF (rescale) THEN
                  term_energy = term_energy + spin_weights(spin)*part
                  operators(:, :, 1, spin) = operators(:, :, 1, spin) + spin_weights(spin)*gamma
               ELSE
                  operators(:, :, 1, spin) = operators(:, :, 1, spin) + gamma
               END IF
            END DO
            valid = ABS(SUM(term_spin_energy) - raw_energy) <= 1.e-9_dp*MAX(1.0_dp, ABS(raw_energy))
            IF (.NOT. valid) THEN
               CPABORT('Auxiliary RI energy/operator mismatch')
            END IF
            IF (rescale) term_spin_energy(:) = spin_weights*term_spin_energy
            IF (geometry) THEN
               CALL hfx_ri_update_forces(qs_env, exchange(contribution, 1)%ri_data, ns, fraction, rho_ao=density, &
                                         use_virial=use_virial, atom_gradient=term_gradient, strain_gradient=term_strain, &
                                         spin_weights=spin_weights)
            END IF
         END IF
         energy = energy + term_energy
         spin_energies(:) = spin_energies + term_spin_energy
         IF (geometry) THEN
            atom_gradient(:, :) = atom_gradient + term_gradient
            strain_gradient = strain_gradient + term_strain
         END IF
      END DO
      IF (bloch) THEN
         IF (ALLOCATED(state%exchange_k)) DEALLOCATE (state%exchange_k)
         CALL MOVE_ALLOC(saved_exchange, state%exchange_k)
      END IF
      CALL dbcsr_deallocate_matrix_set(term_ks)
   END SUBROUTINE auxiliary_ri

! **************************************************************************************************
!> \brief Gather the actual native projection in original coordinates for every physical k.
!>        The two bases keep their own folding phases and symmetry actions. No projection is
!>        refactorized here. Native FM gathering is collective within its owning k-point group.
!> \param qs_env Prepared native primary and auxiliary projection owners.
!> \param projection New A(k), with auxiliary rows, orbital columns and no mesh weights.
!> \param reduced_projection Optional actual projection before physical star expansion.
!> \param metric_inverse Optional inverse of its retained native metric in representative order.
!> \param interpolation_factors Optional occupied projected MOs, including sqrt(f) and density scaling.
!>        These specify an interpolation subspace, not a factorization of a mixed or purified density.
!> \param interpolation_ranks Number of populated columns at each physical k point and spin.
!> \param interpolation_ready Whether all native MO factors were available and finite.
!>        Supply all three interpolation arguments together. The RI density and MOs stay unchanged.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_aux_projection(qs_env, projection, reduced_projection, metric_inverse, &
                                            interpolation_factors, interpolation_ranks, interpolation_ready)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: projection(:, :, :)
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT), OPTIONAL :: reduced_projection(:, :, :), &
         metric_inverse(:, :, :), interpolation_factors(:, :, :, :)
      INTEGER, ALLOCATABLE, INTENT(OUT), OPTIONAL        :: interpolation_ranks(:, :)
      LOGICAL, INTENT(OUT), OPTIONAL                     :: interpolation_ready

      COMPLEX(dp), ALLOCATABLE                           :: inverse(:, :), native(:, :, :), &
                                                            native_metric(:, :, :), &
                                                            original(:, :, :)
      INTEGER                                            :: ik, local, na, no
      REAL(dp)                                           :: inverse_error
      REAL(dp), ALLOCATABLE                              :: real_metric(:, :), real_part(:, :)
      TYPE(admm_type), POINTER                           :: admm
      TYPE(kpoint_env_type), POINTER                     :: kp
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_upf_inverse_scf_type), POINTER             :: auxiliary, primary

      CPASSERT(PRESENT(interpolation_factors) .EQV. PRESENT(interpolation_ranks))
      CPASSERT(PRESENT(interpolation_factors) .EQV. PRESENT(interpolation_ready))
      CALL get_qs_env(qs_env, ks_env=ks_env, admm_env=admm)
      CPASSERT(ASSOCIATED(admm))
      CALL get_ks_env(ks_env, upf_inverse=primary)
      auxiliary => admm%upf_inverse
      CPASSERT(ASSOCIATED(primary) .AND. ASSOCIATED(auxiliary))
      CPASSERT(primary%enabled)
      CPASSERT(admm%upf_aux_metric%ready .AND. admm%upf_mixed_metric%ready)
      na = admm%nao_aux_fit
      no = admm%nao_orb
      ALLOCATE (real_part(na, no))
      IF (ASSOCIATED(auxiliary%bloch%kpoints)) THEN
         kpoints => auxiliary%bloch%kpoints
         CPASSERT(ASSOCIATED(primary%bloch%kpoints, kpoints))
         ALLOCATE (native(na, no, kpoints%nkp))
         native = CMPLX(0.0_dp, 0.0_dp, dp)
         IF (PRESENT(metric_inverse)) THEN
            ALLOCATE (native_metric(na, na, kpoints%nkp), inverse(na, na), real_metric(na, na))
            native_metric = CMPLX(0.0_dp, 0.0_dp, dp)
         END IF
         DO local = 1, SIZE(kpoints%kp_aux_env)
            kp => kpoints%kp_aux_env(local)%kpoint_env
            ik = kp%nkpoint
            IF (kpoints%use_real_wfn) THEN
               CALL cp_fm_get_submatrix(kp%amat, real_part)
               IF (kpoints%para_env_kp%is_source()) native(:, :, ik) = CMPLX(real_part, 0.0_dp, dp)
            ELSE
               CALL cp_cfm_get_submatrix(kp%camat, native(:, :, ik))
               IF (.NOT. kpoints%para_env_kp%is_source()) native(:, :, ik) = CMPLX(0.0_dp, 0.0_dp, dp)
            END IF
            IF (PRESENT(metric_inverse)) THEN
               IF (kpoints%use_real_wfn) THEN
                  CPASSERT(ASSOCIATED(kp%smat%matrix_struct))
                  CALL cp_fm_get_submatrix(kp%smat, real_metric)
                  IF (kpoints%para_env_kp%is_source()) native_metric(:, :, ik) = CMPLX(real_metric, 0.0_dp, dp)
               ELSE
                  CPASSERT(ASSOCIATED(kp%csmat%matrix_struct))
                  CALL cp_cfm_get_submatrix(kp%csmat, native_metric(:, :, ik))
                  IF (.NOT. kpoints%para_env_kp%is_source()) native_metric(:, :, ik) = CMPLX(0.0_dp, 0.0_dp, dp)
               END IF
            END IF
         END DO
         DO ik = 1, kpoints%nkp
            CALL kpoints%para_env%sum(native(:, :, ik))
         END DO
         CALL qs_upf_inverse_kpoints_original_operator(auxiliary%bloch, native, original, right_frame=primary%bloch)
         CALL kpoint_star_density(auxiliary%bloch%star, original, projection, right_map=primary%bloch%star)
         IF (PRESENT(reduced_projection)) ALLOCATE (reduced_projection, SOURCE=original)
         IF (PRESENT(metric_inverse)) THEN
            DO ik = 1, kpoints%nkp
               CALL kpoints%para_env%sum(native_metric(:, :, ik))
               CALL invert_matrix(native_metric(:, :, ik), inverse, inverse_error)
               CPASSERT(inverse_error >= 0.0_dp)
               native_metric(:, :, ik) = inverse
            END DO
            CALL qs_upf_inverse_kpoints_original_operator(auxiliary%bloch, native_metric, metric_inverse)
         END IF
      ELSE
         CPASSERT(.NOT. ASSOCIATED(primary%bloch%kpoints))
         CPASSERT(ASSOCIATED(admm%A))
         CALL cp_fm_get_submatrix(admm%A, real_part)
         ALLOCATE (projection(na, no, 1))
         projection(:, :, 1) = CMPLX(real_part, 0.0_dp, dp)
         IF (PRESENT(reduced_projection)) ALLOCATE (reduced_projection, SOURCE=projection)
         IF (PRESENT(metric_inverse)) THEN
            ALLOCATE (real_metric(na, na), metric_inverse(na, na, 1))
            CALL cp_fm_get_submatrix(admm%S_inv, real_metric)
            metric_inverse(:, :, 1) = CMPLX(real_metric, 0.0_dp, dp)
         END IF
      END IF
      IF (PRESENT(interpolation_factors)) THEN
         CALL auxiliary_orbital_factors(qs_env, admm, primary, projection, &
                                        interpolation_factors, interpolation_ranks, interpolation_ready)
      END IF
   END SUBROUTINE qs_upf_inverse_aux_projection

! **************************************************************************************************
!> \brief Project the current native occupied MOs without replacing the actual RI build density.
!>        Primary coefficients and A use original atom coordinates. Apply the same physical star
!>        to both before multiplication. No mesh weight enters a per-k interpolation factor.
!> \param qs_env Native MO and communicator owners.
!> \param admm Current independent auxiliary basis and charge convention.
!> \param primary Original-coordinate primary fields and physical symmetry map.
!> \param projection Complete physical A(k), including independent row/column folding phases.
!> \param factors New auxiliary factors, padded to the maximum occupied rank.
!> \param ranks Number of occupied columns, indexed by physical k and spin.
!> \param ready All native factors are available. Failure returns no usable interpolation space.
! **************************************************************************************************
   SUBROUTINE auxiliary_orbital_factors(qs_env, admm, primary, projection, factors, ranks, ready)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(admm_type), INTENT(IN)                        :: admm
      TYPE(qs_upf_inverse_scf_type), INTENT(IN)          :: primary
      COMPLEX(dp), INTENT(IN)                            :: projection(:, :, :)
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: factors(:, :, :, :)
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: ranks(:, :)
      LOGICAL, INTENT(OUT)                               :: ready

      COMPLEX(dp), ALLOCATABLE                           :: block(:, :), coefficients(:, :, :, :), &
                                                            expanded(:, :, :)
      INTEGER                                            :: count, failed, ik, local, max_rank, na, &
                                                            nk, nl, no, nr, ns, spin
      INTEGER, ALLOCATABLE                               :: native_ranks(:, :)
      LOGICAL                                            :: available, bloch, complex_mos, owner
      REAL(dp)                                           :: scale
      TYPE(dft_control_type), POINTER                    :: control
      TYPE(kpoint_env_type), POINTER                     :: kp
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mo_set_type), POINTER                         :: mos(:), native
      TYPE(mp_para_env_type), POINTER                    :: para

      CALL get_qs_env(qs_env, dft_control=control, mos=mos, para_env=para)
      na = SIZE(projection, 1)
      no = SIZE(projection, 2)
      nk = SIZE(projection, 3)
      ns = control%nspins
      ALLOCATE (ranks(nk, ns))
      ranks = 0
      ready = .FALSE.
      IF (.NOT. admm%upf_density%valid) THEN
         ALLOCATE (factors(na, 0, nk, ns))
         RETURN
      END IF
      bloch = ASSOCIATED(primary%bloch%kpoints)
      nr = 1
      nl = 1
      complex_mos = .FALSE.
      IF (bloch) THEN
         kpoints => primary%bloch%kpoints
         nr = kpoints%nkp
         nl = SIZE(kpoints%kp_env)
         complex_mos = .NOT. kpoints%use_real_wfn
         CPASSERT(nk == SIZE(primary%bloch%weights))
      END IF
      ALLOCATE (native_ranks(nr, ns))
      native_ranks = 0
      failed = 0
      IF (.NOT. bloch) THEN
         IF (.NOT. ASSOCIATED(mos)) THEN
            failed = 1
         ELSE IF (SIZE(mos) < ns) THEN
            failed = 1
         END IF
      END IF
      CALL para%max(failed)
      IF (failed /= 0) THEN
         ALLOCATE (factors(na, 0, nk, ns))
         RETURN
      END IF
      DO spin = 1, ns
         DO local = 1, nl
            ik = 1
            owner = para%is_source()
            IF (bloch) THEN
               kp => kpoints%kp_env(local)%kpoint_env
               ik = kp%nkpoint
               IF (.NOT. ASSOCIATED(kp%mos)) THEN
                  failed = 1
                  CYCLE
               END IF
               IF (SIZE(kp%mos) < ns) THEN
                  failed = 1
                  CYCLE
               END IF
               native => kp%mos(spin)
               owner = kpoints%para_env_kp%is_source()
            ELSE
               native => mos(spin)
            END IF
            available = occupied_mos_available(native, no, complex_mos, count)
            IF (.NOT. available) failed = 1
            IF (available .AND. owner) native_ranks(ik, spin) = count
         END DO
      END DO
      CALL para%max(failed)
      CALL para%sum(native_ranks)
      IF (failed /= 0) THEN
         ALLOCATE (factors(na, 0, nk, ns))
         RETURN
      END IF
      max_rank = MAXVAL(native_ranks)
      ALLOCATE (factors(na, max_rank, nk, ns), coefficients(no, max_rank, nr, ns))
      factors = CMPLX(0.0_dp, 0.0_dp, dp)
      coefficients = CMPLX(0.0_dp, 0.0_dp, dp)
      IF (max_rank == 0) THEN
         ready = .TRUE.
         RETURN
      END IF
      DO spin = 1, ns
         DO local = 1, nl
            ik = 1
            owner = para%is_source()
            IF (bloch) THEN
               kp => kpoints%kp_env(local)%kpoint_env
               ik = kp%nkpoint
               native => kp%mos(spin)
               owner = kpoints%para_env_kp%is_source()
            ELSE
               native => mos(spin)
            END IF
            count = native_ranks(ik, spin)
            CALL gather_occupied_mos(native, complex_mos, block, available)
            IF (.NOT. available) failed = 1
            IF (available .AND. owner) coefficients(:, :count, ik, spin) = block
         END DO
      END DO
      CALL para%max(failed)
      IF (failed /= 0) RETURN
      CALL para%sum(coefficients)
      DO spin = 1, ns
         scale = 1.0_dp
         IF (admm%do_admmq .OR. admm%do_admms) THEN
            IF (.NOT. ieee_is_finite(admm%gsi(spin))) RETURN
            IF (admm%gsi(spin) <= 0.0_dp) RETURN
            scale = SQRT(admm%gsi(spin))
         END IF
         IF (bloch) THEN
            CALL kpoint_star_orbitals(primary%bloch%star, coefficients(:, :, :, spin), expanded)
            ranks(:, spin) = native_ranks(primary%bloch%star%representative, spin)
         ELSE
            ALLOCATE (expanded, SOURCE=coefficients(:, :, :, spin))
            ranks(1, spin) = native_ranks(1, spin)
         END IF
         DO ik = 1, nk
            count = ranks(ik, spin)
            factors(:, :count, ik, spin) = scale*MATMUL(projection(:, :, ik), expanded(:, :count, ik))
         END DO
         DEALLOCATE (expanded)
      END DO
      ready = .TRUE.
   END SUBROUTINE auxiliary_orbital_factors

! **************************************************************************************************
!> \brief Inspect native occupation and coefficient availability before collective MO gathering.
!> \param native Actual primary MO owner.
!> \param nao Expected primary AO dimension.
!> \param complex_mos Whether the native complex FM owner is used.
!> \param occupied Number of strictly positive occupations, without a truncation threshold.
!> \return Whether all occupied columns can be gathered without changing native storage.
! **************************************************************************************************
   LOGICAL FUNCTION occupied_mos_available(native, nao, complex_mos, occupied) RESULT(available)
      TYPE(mo_set_type), INTENT(IN)                      :: native
      INTEGER, INTENT(IN)                                :: nao
      LOGICAL, INTENT(IN)                                :: complex_mos
      INTEGER, INTENT(OUT)                               :: occupied

      available = .FALSE.
      occupied = 0
      IF (native%nao /= nao .OR. native%nmo < 0 .OR. native%nmo > nao) RETURN
      IF (native%nmo == 0) THEN
         available = .TRUE.
         RETURN
      END IF
      IF (.NOT. ASSOCIATED(native%occupation_numbers)) RETURN
      IF (SIZE(native%occupation_numbers) < native%nmo) RETURN
      IF (.NOT. ALL(ieee_is_finite(native%occupation_numbers(:native%nmo)))) RETURN
      IF (ANY(native%occupation_numbers(:native%nmo) < 0.0_dp)) RETURN
      occupied = COUNT(native%occupation_numbers(:native%nmo) > 0.0_dp)
      IF (occupied > 0) THEN
         IF (complex_mos) THEN
            IF (.NOT. ASSOCIATED(native%cmo_coeff)) RETURN
            IF (.NOT. ASSOCIATED(native%cmo_coeff%matrix_struct)) RETURN
         ELSE
            IF (.NOT. ASSOCIATED(native%mo_coeff)) RETURN
            IF (.NOT. ASSOCIATED(native%mo_coeff%matrix_struct)) RETURN
            IF (native%use_mo_coeff_b) THEN
               IF (.NOT. ASSOCIATED(native%mo_coeff_b)) RETURN
            END IF
         END IF
      END IF
      available = .TRUE.
   END FUNCTION occupied_mos_available

! **************************************************************************************************
!> \brief Gather C sqrt(f), preserving both native FM and DBCSR MO storage.
!> \param native Prepared primary MO owner with nonnegative occupations.
!> \param complex_mos Whether the complex FM owner is used.
!> \param factors New occupied interpolation columns, without k or auxiliary charge scaling.
!> \param available Whether every occupied coefficient is finite.
! **************************************************************************************************
   SUBROUTINE gather_occupied_mos(native, complex_mos, factors, available)
      TYPE(mo_set_type), INTENT(IN)                      :: native
      LOGICAL, INTENT(IN)                                :: complex_mos
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: factors(:, :)
      LOGICAL, INTENT(OUT)                               :: available

      COMPLEX(dp), ALLOCATABLE                           :: coefficients(:, :)
      INTEGER                                            :: column, first, last, occupied, orbital
      LOGICAL                                            :: temporary_used
      REAL(dp), ALLOCATABLE                              :: real_part(:, :)
      TYPE(cp_fm_type)                                   :: temporary

      available = occupied_mos_available(native, native%nao, complex_mos, occupied)
      CPASSERT(available)
      ALLOCATE (factors(native%nao, occupied))
      factors = CMPLX(0.0_dp, 0.0_dp, dp)
      IF (occupied == 0) RETURN
      temporary_used = .FALSE.
      IF (.NOT. complex_mos .AND. native%use_mo_coeff_b) THEN
         CALL cp_fm_create(temporary, native%mo_coeff%matrix_struct, name='Auxiliary occupied interpolation')
         CALL copy_dbcsr_to_fm(native%mo_coeff_b, temporary)
         temporary_used = .TRUE.
      END IF
      column = 0
      orbital = 1
      DO WHILE (orbital <= native%nmo)
         IF (native%occupation_numbers(orbital) == 0.0_dp) THEN
            orbital = orbital + 1
            CYCLE
         END IF
         first = orbital
         last = orbital
         DO WHILE (last < native%nmo)
            IF (native%occupation_numbers(last + 1) == 0.0_dp) EXIT
            last = last + 1
         END DO
         ! Unoccupied native columns may be unset. Gather only populated contiguous spans.
         ALLOCATE (coefficients(native%nao, last - first + 1))
         IF (complex_mos) THEN
            CALL cp_cfm_get_submatrix(native%cmo_coeff, coefficients, start_col=first, n_cols=last - first + 1)
         ELSE
            ALLOCATE (real_part(native%nao, last - first + 1))
            IF (temporary_used) THEN
               CALL cp_fm_get_submatrix(temporary, real_part, start_col=first, n_cols=last - first + 1)
            ELSE
               CALL cp_fm_get_submatrix(native%mo_coeff, real_part, start_col=first, n_cols=last - first + 1)
            END IF
            coefficients(:, :) = CMPLX(real_part, 0.0_dp, dp)
            DEALLOCATE (real_part)
         END IF
         available = ALL(ieee_is_finite(REAL(coefficients, dp))) .AND. ALL(ieee_is_finite(AIMAG(coefficients)))
         IF (.NOT. available) EXIT
         DO orbital = first, last
            column = column + 1
            factors(:, column) = SQRT(native%occupation_numbers(orbital))*coefficients(:, orbital - first + 1)
         END DO
         DEALLOCATE (coefficients)
         orbital = last + 1
      END DO
      IF (temporary_used) CALL cp_fm_release(temporary)
   END SUBROUTINE gather_occupied_mos

! **************************************************************************************************
!> \brief Compose original auxiliary Fock/XC and transfer their complete electronic derivative.
!>        Use the actual native projection A and its complete original charge metric. For charge
!>        normalization, Dbar_A=eta A D_O A^dagger and eta=N_O/N_A. The global spin response
!>        lambda includes the density-normalization response and explicit exchange-ratio derivatives. Geometry outputs
!>        additionally differentiate A=G_AA^-1 G_AO. Their scalar pairing applies mesh weights once.
!>        Reduce derivatives before differentiating the native projection. Its metric is the
!>        native star average, not a separately factorized metric for every physical partner.
!>        This component does not activate full ADMM SCF or supply explicit Fock/XC geometry terms.
!> \param qs_env Same native projection and auxiliary density used by the Fock producer.
!> \param xc_section Selected auxiliary XC functional, with its actual exchange scaling.
!> \param fock Original-coordinate auxiliary Fock operators at the current auxiliary density,
!>        including HF fractions and ADMMP eta**2 scaling, but no k weights or PAW core exchange.
!> \param grid_energy Auxiliary uniform-grid XC energy.
!> \param atom_energy Auxiliary original atomic residual and selected PAW XC energy.
!> \param operators New primary operators, in original physical k and spin order, without weights.
!> \param mixed_covector Optional original-coordinate covectors of G_AO, required with both below.
!> \param auxiliary_covector Combined projection and normalization covectors of G_AA.
!> \param primary_covector Normalization covectors of the actual primary solver metric G_OO.
!> \param exchange Optional actual native HF interactions. The supplied Fock operator excludes PAW core exchange.
!> \param core_energy Optional PAW core-core and core-valence HF energy, required with exchange.
!>        These terms use the auxiliary density and contribute once to its transferred operator.
!> \param valence_spin_energy Optional physical valence energies, including eta**2 and k weights.
!>        The explicit ADMMP ratio derivative then retains constant affine energy terms.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_aux_operator_evaluate(qs_env, xc_section, fock, grid_energy, atom_energy, operators, &
                                                   mixed_covector, auxiliary_covector, primary_covector, exchange, core_energy, &
                                                   valence_spin_energy)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: xc_section
      COMPLEX(dp), INTENT(IN)                            :: fock(:, :, :, :)
      REAL(dp), INTENT(OUT)                              :: grid_energy, atom_energy
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: operators(:, :, :, :)
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT), OPTIONAL    :: mixed_covector(:, :, :, :), &
                                                            auxiliary_covector(:, :, :, :), &
                                                            primary_covector(:, :, :, :)
      TYPE(hfx_type), INTENT(IN), OPTIONAL               :: exchange(:, :)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: core_energy
      REAL(dp), INTENT(IN), OPTIONAL                     :: valence_spin_energy(:)

      COMPLEX(dp), ALLOCATABLE                           :: auxiliary_operator(:, :, :, :), &
                                                            core_operator(:, :, :, :)
      LOGICAL                                            :: valid
      REAL(dp), ALLOCATABLE                              :: core_derivative(:), eta_derivative(:)

      CPASSERT(PRESENT(exchange) .EQV. PRESENT(core_energy))
      ALLOCATE (eta_derivative(SIZE(fock, 4)), core_derivative(SIZE(fock, 4)))
      CALL qs_upf_inverse_aux_xc_rescale(qs_env, xc_section, grid_energy, atom_energy, auxiliary_operator, eta_derivative)
      valid = ALL(SHAPE(fock) == SHAPE(auxiliary_operator))
      CPASSERT(valid)
      IF (PRESENT(exchange)) THEN
         CALL qs_upf_inverse_aux_core_exchange(qs_env, exchange, core_energy, core_operator, eta_derivative=core_derivative)
         eta_derivative(:) = eta_derivative + core_derivative
         auxiliary_operator(:, :, :, :) = auxiliary_operator + core_operator
      END IF
      auxiliary_operator(:, :, :, :) = auxiliary_operator + fock
      CALL hf_ratio_response(qs_env, fock, eta_derivative, valence_spin_energy)
      CALL transfer_operator(qs_env, auxiliary_operator, operators, mixed_covector, auxiliary_covector, primary_covector, &
                             eta_derivative)
   END SUBROUTINE qs_upf_inverse_aux_operator_evaluate

! **************************************************************************************************
!> \brief Add the explicit ADMMP valence-HF ratio derivative without differentiating core or XC twice.
!> \param qs_env Actual native auxiliary density and full physical-k weights.
!> \param fock Valence RI operator, including eta**2 for ADMMP and no mesh weights.
!> \param derivative Accumulated partial_eta F at fixed native auxiliary density.
!> \param valence_spin_energy Optional eta**2-scaled physical valence energies with k weights.
!>        Without these values the supplied operator is assumed to define a quadratic functional.
! **************************************************************************************************
   SUBROUTINE hf_ratio_response(qs_env, fock, derivative, valence_spin_energy)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      COMPLEX(dp), INTENT(IN)                            :: fock(:, :, :, :)
      REAL(dp), INTENT(INOUT)                            :: derivative(:)
      REAL(dp), INTENT(IN), OPTIONAL                     :: valence_spin_energy(:)

      INTEGER                                            :: ik, spin
      LOGICAL                                            :: valid
      TYPE(admm_type), POINTER                           :: admm

      CALL get_qs_env(qs_env, admm_env=admm)
      CPASSERT(ASSOCIATED(admm))
      IF (.NOT. admm%do_admmp) RETURN
      valid = ALL(SHAPE(fock) == SHAPE(admm%upf_density%matrices))
      CPASSERT(valid)
      CPASSERT(SIZE(derivative) == SIZE(fock, 4))
      IF (PRESENT(valence_spin_energy)) THEN
         valid = SIZE(valence_spin_energy) == SIZE(derivative) .AND. ALL(ieee_is_finite(valence_spin_energy))
         CPASSERT(valid)
      END IF
      DO spin = 1, SIZE(fock, 4)
         CPASSERT(admm%gsi(spin) > 0.0_dp)
         IF (PRESENT(valence_spin_energy)) THEN
            derivative(spin) = derivative(spin) + 2.0_dp*valence_spin_energy(spin)/admm%gsi(spin)
         ELSE
            DO ik = 1, SIZE(fock, 3)
               derivative(spin) = derivative(spin) + admm%upf_aux_metric%weights(ik)* &
                                  REAL(SUM(CONJG(fock(:, :, ik, spin))* &
                                           admm%upf_density%matrices(:, :, ik, spin)), dp)/admm%gsi(spin)
            END DO
         END IF
      END DO
   END SUBROUTINE hf_ratio_response

! **************************************************************************************************
!> \brief Transfer the derivative of the same auxiliary energy through the actual A and eta.
!> \param qs_env ...
!> \param auxiliary_operator Complete unweighted original-coordinate derivative at fixed auxiliary density.
!> \param operators ...
!> \param mixed_covector ...
!> \param auxiliary_covector ...
!> \param primary_covector ...
!> \param eta_derivative Optional explicit partial_eta F at fixed native density, required for P/S scaling.
! **************************************************************************************************
   SUBROUTINE transfer_operator(qs_env, auxiliary_operator, operators, &
                                mixed_covector, auxiliary_covector, primary_covector, eta_derivative)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      COMPLEX(dp), INTENT(IN)                            :: auxiliary_operator(:, :, :, :)
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: operators(:, :, :, :)
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT), OPTIONAL    :: mixed_covector(:, :, :, :), &
                                                            auxiliary_covector(:, :, :, :), &
                                                            primary_covector(:, :, :, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: eta_derivative(:)

      COMPLEX(dp), ALLOCATABLE :: auxiliary_reduced(:, :, :), corrected(:, :), &
         corrected_k(:, :, :), density(:, :, :, :), expanded(:, :, :), gamma(:, :), &
         matrix(:, :, :), metric_inverse(:, :, :), mixed_reduced(:, :, :), &
         primary_metric(:, :, :), projection(:, :, :), reduced_corrected(:, :, :), &
         reduced_density(:, :, :), reduced_projection(:, :, :), unscaled(:, :), work(:, :), z(:, :)
      INTEGER                                            :: ik, na, nk, no, ns, spin
      LOGICAL                                            :: geometry, normalized, purified, &
                                                            ratio_dependent, scaled, valid
      REAL(dp)                                           :: density_scale, eta, explicit_derivative, &
                                                            lambda
      TYPE(admm_type), POINTER                           :: admm
      TYPE(dbcsr_p_type), POINTER                        :: density_ao(:), metric(:)
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_upf_inverse_scf_type), POINTER             :: primary

      CALL get_qs_env(qs_env, ks_env=ks_env, admm_env=admm, rho=rho)
      CPASSERT(ASSOCIATED(admm))
      CPASSERT(admm%upf_density%valid)
      purified = admm%purification_method == do_admm_purify_mo_diag
      IF ((admm%purification_method /= do_admm_purify_none .AND. .NOT. purified) .OR. &
          admm%block_fit .OR. admm%block_dm) THEN
         CPABORT('This UPF ADMM purification or blocked fit is not implemented yet')
      END IF
      scaled = admm%do_admmp .OR. admm%do_admms
      IF (admm%scaling_model /= do_admm_exch_scaling_none) THEN
         CPASSERT(scaled .AND. admm%scaling_model == do_admm_exch_scaling_merlot)
      END IF
      IF (scaled) THEN
         CPASSERT(PRESENT(eta_derivative))
         IF (purified) THEN
            CPABORT('Exchange rescaling requires no purification')
         END IF
      END IF
      geometry = PRESENT(mixed_covector) .OR. PRESENT(auxiliary_covector) .OR. PRESENT(primary_covector)
      IF (geometry) THEN
         valid = PRESENT(mixed_covector) .AND. PRESENT(auxiliary_covector) .AND. PRESENT(primary_covector)
         CPASSERT(valid)
      END IF
      normalized = admm%do_admmq .OR. admm%do_admms
      ratio_dependent = normalized .OR. admm%do_admmp
      IF (normalized) THEN
         CPASSERT(admm%charge_constrain)
      END IF
      CALL get_ks_env(ks_env, upf_inverse=primary)
      CPASSERT(ASSOCIATED(primary))
      CPASSERT(primary%enabled)
      IF (purified .AND. ASSOCIATED(primary%bloch%kpoints)) THEN
         CPABORT('MO_DIAG purification requires real Gamma orbitals')
      END IF
      IF (geometry) THEN
         CALL qs_upf_inverse_aux_projection(qs_env, projection, reduced_projection, metric_inverse)
      ELSE
         CALL qs_upf_inverse_aux_projection(qs_env, projection)
      END IF
      na = SIZE(projection, 1)
      no = SIZE(projection, 2)
      nk = SIZE(projection, 3)
      ns = SIZE(admm%upf_density%matrices, 4)
      valid = ALL(SHAPE(auxiliary_operator) == [na, na, nk, ns])
      CPASSERT(valid)
      CPASSERT(SIZE(admm%upf_aux_metric%weights) == nk)
      ALLOCATE (density(no, no, nk, ns), primary_metric(no, no, nk))
      IF (ASSOCIATED(primary%bloch%kpoints)) THEN
         valid = ALL(SHAPE(primary%metric_k) == [no, no, nk])
         CPASSERT(valid)
         primary_metric(:, :, :) = primary%metric_k
         DO spin = 1, ns
            CALL qs_upf_inverse_kpoints_density(primary%bloch, spin, matrix)
            density(:, :, :, spin) = matrix
         END DO
      ELSE
         CPASSERT(nk == 1)
         metric => get_1d_pointer(primary%metric)
         CALL qs_upf_inverse_matrix_gather(primary%frame%ao, metric(1)%matrix, gamma)
         primary_metric(:, :, 1) = gamma
         CALL qs_rho_get(rho, rho_ao=density_ao)
         DO spin = 1, ns
            CALL qs_upf_inverse_matrix_gather(primary%frame%ao, density_ao(spin)%matrix, gamma)
            density(:, :, 1, spin) = gamma
         END DO
      END IF
      IF (PRESENT(eta_derivative)) THEN
         CPASSERT(SIZE(eta_derivative) == ns)
      END IF
      ALLOCATE (operators(no, no, nk, ns), unscaled(na, na), corrected(na, na), work(na, no), z(na, no), &
                corrected_k(na, na, nk))
      IF (geometry) THEN
         ALLOCATE (mixed_covector(na, no, nk, ns), auxiliary_covector(na, na, nk, ns), &
                   primary_covector(no, no, nk, ns), &
                   mixed_reduced(na, no, SIZE(reduced_projection, 3)), &
                   auxiliary_reduced(na, na, SIZE(reduced_projection, 3)))
      END IF
      DO spin = 1, ns
         eta = 1.0_dp
         lambda = 0.0_dp
         density_scale = 1.0_dp
         IF (ratio_dependent) THEN
            eta = admm%gsi(spin)
            IF (normalized) density_scale = eta
            explicit_derivative = 0.0_dp
            IF (PRESENT(eta_derivative)) explicit_derivative = eta_derivative(spin)
            CALL admm_charge_scaled_operator(auxiliary_operator(:, :, :, spin), admm%upf_density%matrices(:, :, :, spin), &
                                  admm%upf_aux_metric%values, density(:, :, :, spin), primary_metric, admm%upf_aux_metric%weights, &
                                             eta, normalized, explicit_derivative, corrected_k, lambda)
         ELSE
            corrected_k(:, :, :) = auxiliary_operator(:, :, :, spin)
         END IF
         DO ik = 1, nk
            unscaled(:, :) = admm%upf_density%matrices(:, :, ik, spin)/density_scale
            corrected(:, :) = corrected_k(:, :, ik)
            IF (purified) THEN
               IF (geometry) THEN
                  CALL transfer_purified(qs_env, spin, projection(:, :, ik), primary_metric(:, :, ik), &
                                         density(:, :, ik, spin), unscaled, corrected, operators(:, :, ik, spin), &
                                         metric_inverse(:, :, ik), mixed_covector(:, :, ik, spin), &
                                         auxiliary_covector(:, :, ik, spin), primary_covector(:, :, ik, spin))
               ELSE
                  CALL transfer_purified(qs_env, spin, projection(:, :, ik), primary_metric(:, :, ik), &
                                         density(:, :, ik, spin), unscaled, corrected, operators(:, :, ik, spin))
               END IF
               operators(:, :, ik, spin) = operators(:, :, ik, spin) + lambda*primary_metric(:, :, ik)
            ELSE
               work(:, :) = MATMUL(corrected, projection(:, :, ik))
               operators(:, :, ik, spin) = MATMUL(CONJG(TRANSPOSE(projection(:, :, ik))), work) + &
                                           lambda*primary_metric(:, :, ik)
               IF (geometry) THEN
                  auxiliary_covector(:, :, ik, spin) = CMPLX(0.0_dp, 0.0_dp, dp)
                  primary_covector(:, :, ik, spin) = CMPLX(0.0_dp, 0.0_dp, dp)
               END IF
            END IF
            ! Publish the Hermitian covector after the congruence and charge correction.
            operators(:, :, ik, spin) = 0.5_dp*(operators(:, :, ik, spin) + &
                                                CONJG(TRANSPOSE(operators(:, :, ik, spin))))
            IF (geometry) THEN
               corrected_k(:, :, ik) = corrected
               auxiliary_covector(:, :, ik, spin) = auxiliary_covector(:, :, ik, spin) - eta*lambda*unscaled
               primary_covector(:, :, ik, spin) = primary_covector(:, :, ik, spin) + lambda*density(:, :, ik, spin)
            END IF
         END DO
         IF (geometry .AND. .NOT. purified) THEN
            IF (ASSOCIATED(primary%bloch%kpoints)) THEN
               CALL kpoint_star_operator(admm%upf_inverse%bloch%star, corrected_k, reduced_corrected)
               CALL kpoint_star_operator(primary%bloch%star, density(:, :, :, spin), reduced_density)
            ELSE
               IF (ALLOCATED(reduced_corrected)) DEALLOCATE (reduced_corrected, reduced_density)
               ALLOCATE (reduced_corrected, SOURCE=corrected_k)
               ALLOCATE (reduced_density, SOURCE=density(:, :, :, spin))
            END IF
            DO ik = 1, SIZE(reduced_projection, 3)
               work(:, :) = MATMUL(reduced_corrected(:, :, ik), reduced_projection(:, :, ik))
               z(:, :) = MATMUL(metric_inverse(:, :, ik), MATMUL(work, reduced_density(:, :, ik)))
               mixed_reduced(:, :, ik) = 2.0_dp*z
               auxiliary_reduced(:, :, ik) = -MATMUL(z, CONJG(TRANSPOSE(reduced_projection(:, :, ik)))) - &
                                             MATMUL(reduced_projection(:, :, ik), CONJG(TRANSPOSE(z)))
            END DO
            IF (ASSOCIATED(primary%bloch%kpoints)) THEN
               CALL kpoint_star_density(admm%upf_inverse%bloch%star, mixed_reduced, expanded, &
                                        right_map=primary%bloch%star)
               mixed_covector(:, :, :, spin) = expanded
               CALL kpoint_star_density(admm%upf_inverse%bloch%star, auxiliary_reduced, expanded)
               auxiliary_covector(:, :, :, spin) = auxiliary_covector(:, :, :, spin) + expanded
            ELSE
               mixed_covector(:, :, :, spin) = mixed_reduced
               auxiliary_covector(:, :, :, spin) = auxiliary_covector(:, :, :, spin) + auxiliary_reduced
            END IF
         END IF
      END DO
   END SUBROUTINE transfer_operator

! **************************************************************************************************
!> \brief Use the actual native MO_DIAG factors for the original-coordinate AO derivative.
!>        This orbital purification requires a primary density from the same orthonormal MOs.
!>        Density mixing cannot silently replace those MOs by an unrelated occupation matrix.
!> \param qs_env ...
!> \param spin ...
!> \param projection ...
!> \param primary_metric ...
!> \param primary_density ...
!> \param auxiliary_density Unscaled purified density, before charge normalization.
!> \param auxiliary_operator Complete derivative including the charge-normalization multiplier.
!> \param operator ...
!> \param metric_inverse ...
!> \param mixed_covector ...
!> \param auxiliary_covector ...
!> \param primary_covector ...
! **************************************************************************************************
   SUBROUTINE transfer_purified(qs_env, spin, projection, primary_metric, primary_density, auxiliary_density, &
                                auxiliary_operator, OPERATOR, metric_inverse, mixed_covector, &
                                auxiliary_covector, primary_covector)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: spin
      COMPLEX(dp), INTENT(IN) :: projection(:, :), primary_metric(:, :), primary_density(:, :), &
         auxiliary_density(:, :), auxiliary_operator(:, :)
      COMPLEX(dp), INTENT(OUT)                           :: operator(:, :)
      COMPLEX(dp), INTENT(IN), OPTIONAL                  :: metric_inverse(:, :)
      COMPLEX(dp), INTENT(OUT), OPTIONAL                 :: mixed_covector(:, :), &
                                                            auxiliary_covector(:, :), &
                                                            primary_covector(:, :)

      COMPLEX(dp), ALLOCATABLE                           :: check(:, :), coefficients(:, :), &
                                                            gram_vectors(:, :), &
                                                            inverse_root(:, :), weighted(:, :)
      INTEGER                                            :: i, m, no
      LOGICAL                                            :: valid
      REAL(dp), ALLOCATABLE                              :: occupations(:), real_coefficients(:, :), &
                                                            real_gram(:, :)
      TYPE(admm_type), POINTER                           :: admm
      TYPE(mo_set_type), POINTER                         :: mos(:)

      CALL get_qs_env(qs_env, admm_env=admm, mos=mos)
      m = admm%nmo(spin)
      no = SIZE(primary_metric, 1)
      ALLOCATE (real_coefficients(no, m), real_gram(m, m), occupations(m))
      IF (m > 0) THEN
         CALL cp_fm_get_submatrix(mos(spin)%mo_coeff, real_coefficients, n_rows=no, n_cols=m)
         occupations(:) = mos(spin)%occupation_numbers(1:m)
      END IF
      ALLOCATE (coefficients, SOURCE=CMPLX(real_coefficients, 0.0_dp, dp))
      IF (m > 0) CALL cp_fm_get_submatrix(admm%R(spin), real_gram)
      ALLOCATE (gram_vectors, SOURCE=CMPLX(real_gram, 0.0_dp, dp))
      IF (m > 0) CALL cp_fm_get_submatrix(admm%lambda_inv_sqrt(spin), real_gram)
      ALLOCATE (inverse_root, SOURCE=CMPLX(real_gram, 0.0_dp, dp))
      ALLOCATE (check, SOURCE=MATMUL(CONJG(TRANSPOSE(coefficients)), MATMUL(primary_metric, coefficients)))
      DO i = 1, m
         check(i, i) = check(i, i) - 1.0_dp
      END DO
      IF (m > 0) THEN
         valid = MAXVAL(ABS(check)) < 1.0E-8_dp
         CPASSERT(valid)
      END IF
      ALLOCATE (weighted, SOURCE=coefficients)
      DO i = 1, m
         weighted(:, i) = weighted(:, i)*occupations(i)
      END DO
      DEALLOCATE (check)
      ALLOCATE (check, SOURCE=MATMUL(weighted, CONJG(TRANSPOSE(coefficients))))
      valid = MAXVAL(ABS(check - primary_density)) < 1.0E-8_dp*MAX(1.0_dp, MAXVAL(ABS(primary_density)))
      CPASSERT(valid)
      DEALLOCATE (weighted, check)
      ALLOCATE (weighted, SOURCE=MATMUL(projection, MATMUL(coefficients, inverse_root)))
      ALLOCATE (check, SOURCE=weighted)
      DO i = 1, m
         check(:, i) = check(:, i)*occupations(i)
      END DO
      valid = MAXVAL(ABS(MATMUL(check, CONJG(TRANSPOSE(weighted))) - auxiliary_density)) < &
              1.0E-8_dp*MAX(1.0_dp, MAXVAL(ABS(auxiliary_density)))
      CPASSERT(valid)
      CALL admm_purified_operator(projection, primary_metric, admm%upf_aux_metric%values(:, :, 1), &
                                  coefficients, occupations, gram_vectors, admm%eigvals_lambda(spin)%eigvals%data, &
                                  inverse_root, auxiliary_operator, OPERATOR, metric_inverse, mixed_covector, &
                                  auxiliary_covector, primary_covector)
   END SUBROUTINE transfer_purified

! **************************************************************************************************
!> \brief Original PAW closed-core HF energy and its auxiliary AO and geometry derivatives.
!>        Core-core exchange is constant. Core-valence exchange is linear in the original
!>        projector occupations, without a further factor one half. Sum native HF repetitions
!>        once and ignore thread replicas. Valence-valence exchange remains in the RI producer.
!>        No PAW Hartree, semilocal XC or fixed dataset reference energy is added here.
!> \param qs_env Valid current auxiliary occupations and their independent original-UPF representation.
!> \param exchange The same native HF interactions, fractions and physical range parameters as RI.
!> \param energy Selected frozen-core HF energy, including each actual PAW atom once.
!> \param operators Original physical-k/spin auxiliary AO derivative, without mesh weights.
!> \param atom_response Optional total atomic motion, including affine cell motion.
!> \param cell_response Optional physical lattice derivative at fixed reduced k coordinates.
!> \param response Optional fixed-auxiliary-occupation directional derivative. Supply all three geometry arguments.
!> \param eta_derivative Optional explicit ADMMP core-valence ratio derivative at fixed native density.
!>        Core-core energy remains unscaled. Other variants return zero explicit ratio response.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_aux_core_exchange(qs_env, exchange, energy, operators, &
                                               atom_response, cell_response, response, eta_derivative)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(hfx_type), INTENT(IN)                         :: exchange(:, :)
      REAL(dp), INTENT(OUT)                              :: energy
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: operators(:, :, :, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: atom_response(:, :), cell_response(3, 3)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: response, eta_derivative(:)

      COMPLEX(dp), ALLOCATABLE, TARGET                   :: scaled_matrices(:, :, :, :)
      COMPLEX(dp), POINTER                               :: matrices(:, :, :, :)
      INTEGER                                            :: ik, kind, nk, ns, spin
      LOGICAL                                            :: geometry
      REAL(dp)                                           :: energies(3)
      REAL(dp), ALLOCATABLE                              :: core_energy(:), reduced_k_response(:, :)
      TYPE(admm_type), POINTER                           :: admm
      TYPE(paw_potential_type), ALLOCATABLE              :: core_potential(:), potential(:)
      TYPE(qs_kind_type), POINTER                        :: kinds(:)
      TYPE(qs_upf_inverse_scf_type), POINTER             :: state
      TYPE(section_vals_type), POINTER                   :: no_functional

      CALL get_qs_env(qs_env, admm_env=admm, qs_kind_set=kinds)
      CPASSERT(ASSOCIATED(admm))
      CPASSERT(admm%upf_density%valid)
      state => admm%upf_inverse
      CPASSERT(ASSOCIATED(state))
      CPASSERT(state%basis_type == 'AUX_FIT')
      CPASSERT(SIZE(exchange, 2) > 0)
      geometry = PRESENT(response)
      CPASSERT(geometry .EQV. PRESENT(atom_response))
      CPASSERT(geometry .EQV. PRESENT(cell_response))
      energy = 0.0_dp
      ALLOCATE (operators, MOLD=admm%upf_density%matrices)
      operators = CMPLX(0.0_dp, 0.0_dp, dp)
      IF (geometry) response = 0.0_dp
      ns = SIZE(admm%upf_density%matrices, 4)
      IF (PRESENT(eta_derivative)) THEN
         CPASSERT(SIZE(eta_derivative) == ns)
         eta_derivative = 0.0_dp
      END IF
      IF (.NOT. ANY(state%paw%active)) RETURN
      nk = SIZE(admm%upf_density%matrices, 3)
      ALLOCATE (core_energy(SIZE(kinds)), SOURCE=0.0_dp)
      ALLOCATE (core_potential(SIZE(kinds)))
      DO kind = 1, SIZE(kinds)
         IF (.NOT. state%paw(kind)%active) CYCLE
         ASSOCIATE (pot => kinds(kind)%sgp_potential%upf_data, projectors => kinds(kind)%sgp_potential%upf_projectors)
            ALLOCATE (core_potential(kind)%matrix(projectors%nproj, projectors%nproj, 1))
            CALL qs_upf_inverse_paw_core_prepare(pot, projectors, state%paw(kind), exchange, ns, &
                                                 core_energy(kind), core_potential(kind)%matrix(:, :, 1))
         END ASSOCIATE
      END DO
      matrices => admm%upf_density%matrices
      IF (admm%do_admmp) THEN
         ALLOCATE (scaled_matrices, SOURCE=admm%upf_density%matrices)
         matrices => scaled_matrices
         DO spin = 1, ns
            CPASSERT(admm%gsi(spin) > 0.0_dp)
            matrices(:, :, :, spin) = admm%gsi(spin)*admm%upf_density%matrices(:, :, :, spin)
         END DO
      END IF
      NULLIFY (no_functional)
      IF (ASSOCIATED(state%bloch%kpoints)) THEN
         CALL qs_upf_inverse_paw_evaluate(state%paw, state%bloch%frames, matrices, state%bloch%weights, &
                                          no_functional, 0.0_dp, 0.0_dp, energies, operators, potential, &
                                      semilocal_xc=.FALSE., core_energy=core_energy, core_potential=core_potential, hartree=.FALSE.)
      ELSE
         CPASSERT(nk == 1)
         CALL qs_upf_inverse_paw_evaluate(state%paw, [state%frame], matrices, [1.0_dp], &
                                          no_functional, 0.0_dp, 0.0_dp, energies, operators, potential, &
                                      semilocal_xc=.FALSE., core_energy=core_energy, core_potential=core_potential, hartree=.FALSE.)
      END IF
      energy = energies(2)
      IF (geometry) THEN
         ALLOCATE (reduced_k_response(3, nk), SOURCE=0.0_dp)
         IF (ASSOCIATED(state%bloch%kpoints)) THEN
            CALL qs_upf_inverse_paw_response(state%bloch%frames, matrices, state%bloch%weights, potential, &
                                             atom_response, cell_response, reduced_k_response, response)
         ELSE
            CALL qs_upf_inverse_paw_response([state%frame], matrices, [1.0_dp], potential, &
                                             atom_response, cell_response, reduced_k_response, response)
         END IF
      END IF
      IF (admm%do_admmp) THEN
         DO spin = 1, ns
            IF (PRESENT(eta_derivative)) THEN
               DO ik = 1, nk
                  eta_derivative(spin) = eta_derivative(spin) + admm%upf_aux_metric%weights(ik)* &
                                         REAL(SUM(CONJG(operators(:, :, ik, spin))*admm%upf_density%matrices(:, :, ik, spin)), dp)
               END DO
            END IF
            operators(:, :, :, spin) = admm%gsi(spin)*operators(:, :, :, spin)
         END DO
      END IF
   END SUBROUTINE qs_upf_inverse_aux_core_exchange

! **************************************************************************************************
!> \brief Compose auxiliary RI, original PAW core, XC and projection/normalization geometry.
!>        Native RI supplies energy, its auxiliary AO derivative and globally reduced Cartesian
!>        and affine-strain derivatives at the current fixed auxiliary occupation matrices.
!>        All RI quantities include their HF fractions. Energy and geometry include physical
!>        mesh weights once, while the supplied AO derivative has no mesh weights.
!>        The energy and primary AO operator use exactly the same auxiliary functional. Its
!>        explicit geometry is combined with the implicit response of the native A and eta.
!>        This component does not activate the complete native ADMM SCF calculation.
!> \param qs_env The actual native primary and auxiliary density/projection owners.
!> \param xc_section The selected auxiliary XC correction with its actual exchange scaling.
!> \param exchange The same native HF interactions as the supplied RI contribution.
!> \param ri_energy Native valence-valence RI exchange energy at the current auxiliary density.
!> \param ri_operator Its unweighted original-coordinate auxiliary AO derivative (AO,AO,k,spin).
!> \param ri_atom_gradient Globally reduced positive Cartesian RI energy gradients (3,atom).
!> \param ri_strain_gradient Globally reduced positive dE_RI/dstrain, including affine atomic motion.
!> \param atom_response Total physical atomic motion, including affine cell motion.
!> \param cell_response Physical lattice derivative, keeping reduced k coordinates fixed.
!> \param energy Complete auxiliary RI plus core-HF and XC correction energy.
!> \param operators Complete primary AO derivative through the actual projection and charge constraint.
!> \param response Complete directional derivative at fixed original-coordinate primary occupations.
!> \param components Optional RI valence HF, PAW core HF, grid XC and atomic XC energies.
!> \param valence_spin_energy Optional physical valence energies with eta**2 and k weights, as supplied by RI.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_aux_response(qs_env, xc_section, exchange, ri_energy, ri_operator, &
                                          ri_atom_gradient, ri_strain_gradient, atom_response, cell_response, &
                                          energy, operators, response, components, valence_spin_energy)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: xc_section
      TYPE(hfx_type), INTENT(IN)                         :: exchange(:, :)
      REAL(dp), INTENT(IN)                               :: ri_energy
      COMPLEX(dp), INTENT(IN)                            :: ri_operator(:, :, :, :)
      REAL(dp), INTENT(IN)                               :: ri_atom_gradient(:, :), &
                                                            ri_strain_gradient(3, 3), &
                                                            atom_response(:, :), &
                                                            cell_response(3, 3)
      REAL(dp), INTENT(OUT)                              :: energy
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: operators(:, :, :, :)
      REAL(dp), INTENT(OUT)                              :: response
      REAL(dp), INTENT(OUT), OPTIONAL                    :: components(4)
      REAL(dp), INTENT(IN), OPTIONAL                     :: valence_spin_energy(:)

      COMPLEX(dp), ALLOCATABLE                           :: xc_operator(:, :, :, :)
      REAL(dp)                                           :: atom_energy, grid_energy, xc_response
      REAL(dp), ALLOCATABLE                              :: eta_derivative(:)

      ALLOCATE (eta_derivative(SIZE(ri_operator, 4)))
      CALL qs_upf_inverse_aux_xc_rescale(qs_env, xc_section, grid_energy, atom_energy, xc_operator, eta_derivative, &
                                         atom_response=atom_response, cell_response=cell_response, response=xc_response)
      CALL compose_auxiliary_response(qs_env, exchange, ri_energy, ri_operator, ri_atom_gradient, ri_strain_gradient, &
                                      atom_response, cell_response, grid_energy, atom_energy, xc_operator, eta_derivative, &
                                      xc_response, energy, operators, response, components, valence_spin_energy)
   END SUBROUTINE qs_upf_inverse_aux_response

! **************************************************************************************************
!> \brief Compose one analytic geometry direction from already evaluated auxiliary XC quantities.
!>        The XC energy, operator and ratio derivative are independent of this direction.
!>        RI, original core, projection and normalization responses keep their original meanings.
! **************************************************************************************************
   SUBROUTINE compose_auxiliary_response(qs_env, exchange, ri_energy, ri_operator, &
                                         ri_atom_gradient, ri_strain_gradient, atom_response, cell_response, &
                                         grid_energy, atom_energy, xc_operator, xc_eta, xc_response, &
                                         energy, operators, response, components, valence_spin_energy)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(hfx_type), INTENT(IN)                         :: exchange(:, :)
      REAL(dp), INTENT(IN)                               :: ri_energy
      COMPLEX(dp), INTENT(IN)                            :: ri_operator(:, :, :, :)
      REAL(dp), INTENT(IN) :: ri_atom_gradient(:, :), ri_strain_gradient(3, 3), &
         atom_response(:, :), cell_response(3, 3), grid_energy, atom_energy
      COMPLEX(dp), INTENT(IN)                            :: xc_operator(:, :, :, :)
      REAL(dp), INTENT(IN)                               :: xc_eta(:), xc_response
      REAL(dp), INTENT(OUT)                              :: energy
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: operators(:, :, :, :)
      REAL(dp), INTENT(OUT)                              :: response
      REAL(dp), INTENT(OUT), OPTIONAL                    :: components(4)
      REAL(dp), INTENT(IN), OPTIONAL                     :: valence_spin_energy(:)

      COMPLEX(dp), ALLOCATABLE :: auxiliary_covector(:, :, :, :), auxiliary_operator(:, :, :, :), &
         core_operator(:, :, :, :), mixed_covector(:, :, :, :), primary_covector(:, :, :, :)
      INTEGER                                            :: atom
      LOGICAL                                            :: valid
      REAL(dp)                                           :: core_energy, core_response, &
                                                            projection_response, ri_response, &
                                                            strain(3, 3)
      REAL(dp), ALLOCATABLE                              :: core_derivative(:), eta_derivative(:)
      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_type), POINTER                       :: atoms(:)

      CALL get_qs_env(qs_env, cell=cell, particle_set=atoms)
      valid = ALL(SHAPE(atom_response) == [3, SIZE(atoms)]) .AND. &
              ALL(SHAPE(ri_atom_gradient) == SHAPE(atom_response))
      CPASSERT(valid)
      valid = ieee_is_finite(ri_energy) .AND. ALL(ieee_is_finite(ri_atom_gradient)) .AND. &
              ALL(ieee_is_finite(ri_strain_gradient))
      CPASSERT(valid)
      ALLOCATE (eta_derivative(SIZE(ri_operator, 4)), core_derivative(SIZE(ri_operator, 4)))
      valid = ALL(SHAPE(ri_operator) == SHAPE(xc_operator))
      CPASSERT(valid .AND. SIZE(xc_eta) == SIZE(eta_derivative))
      ALLOCATE (auxiliary_operator, SOURCE=xc_operator)
      eta_derivative(:) = xc_eta
      CALL qs_upf_inverse_aux_core_exchange(qs_env, exchange, core_energy, core_operator, &
                                            atom_response, cell_response, core_response, eta_derivative=core_derivative)
      eta_derivative(:) = eta_derivative + core_derivative
      CALL hf_ratio_response(qs_env, ri_operator, eta_derivative, valence_spin_energy)
      auxiliary_operator(:, :, :, :) = auxiliary_operator + core_operator
      auxiliary_operator(:, :, :, :) = auxiliary_operator + ri_operator
      CALL transfer_operator(qs_env, auxiliary_operator, operators, mixed_covector, auxiliary_covector, primary_covector, &
                             eta_derivative)
      CALL qs_upf_inverse_aux_metric_response(qs_env, mixed_covector, auxiliary_covector, primary_covector, &
                                              atom_response, cell_response, projection_response)
      strain = MATMUL(cell_response, cell%h_inv)
      ri_response = SUM(ri_strain_gradient*strain)
      DO atom = 1, SIZE(atoms)
         ri_response = ri_response + SUM(ri_atom_gradient(:, atom)* &
                                         (atom_response(:, atom) - MATMUL(strain, atoms(atom)%r)))
      END DO
      energy = ri_energy + core_energy + grid_energy + atom_energy
      IF (PRESENT(components)) components = [ri_energy, core_energy, grid_energy, atom_energy]
      response = ri_response + core_response + xc_response + projection_response
   END SUBROUTINE compose_auxiliary_response

! **************************************************************************************************
!> \brief Contract projection and charge-normalization covectors with the actual moving metrics.
!>        Each native metric service applies its defining quadrature and k weights once. The
!>        primary normalization uses the installed solver metric, including the complete GAPW
!>        density metric where present. Frozen core charges never enter these valence metrics.
!>        This is the implicit projection/normalization response. Explicit auxiliary Fock and
!>        XC geometry derivatives at fixed auxiliary occupation matrices remain separate terms.
!> \param qs_env Same native primary and auxiliary owners used to obtain the covectors.
!> \param mixed_covector Original-coordinate G_AO covectors, without mesh weights.
!> \param auxiliary_covector Combined projection and normalization G_AA covectors.
!> \param primary_covector Normalization G_OO covectors in the actual solver basis.
!> \param atom_response Total atomic motion, including any affine cell motion.
!> \param cell_response Lattice-vector derivative at fixed reduced k coordinates.
!> \param response Real directional energy derivative, without a force or virial minus sign.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_aux_metric_response(qs_env, mixed_covector, auxiliary_covector, primary_covector, &
                                                 atom_response, cell_response, response)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      COMPLEX(dp), INTENT(IN)                            :: mixed_covector(:, :, :, :), &
                                                            auxiliary_covector(:, :, :, :), &
                                                            primary_covector(:, :, :, :)
      REAL(dp), INTENT(IN)                               :: atom_response(:, :), cell_response(3, 3)
      REAL(dp), INTENT(OUT)                              :: response

      COMPLEX(dp), ALLOCATABLE                           :: auxiliary(:, :, :), mixed(:, :, :), &
                                                            orbital(:, :, :)
      INTEGER                                            :: na, nk, no, ns
      LOGICAL                                            :: valid
      REAL(dp)                                           :: auxiliary_response, mixed_response, &
                                                            orbital_response
      TYPE(admm_type), POINTER                           :: admm
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_upf_inverse_scf_type), POINTER             :: primary

      CALL get_qs_env(qs_env, ks_env=ks_env, admm_env=admm)
      CPASSERT(ASSOCIATED(admm))
      CPASSERT(ASSOCIATED(admm%upf_inverse))
      CPASSERT(admm%upf_aux_metric%ready .AND. admm%upf_mixed_metric%ready)
      CALL get_ks_env(ks_env, upf_inverse=primary)
      CPASSERT(ASSOCIATED(primary))
      na = SIZE(admm%upf_aux_metric%values, 1)
      no = SIZE(admm%upf_mixed_metric%values, 2)
      nk = SIZE(admm%upf_aux_metric%weights)
      ns = SIZE(mixed_covector, 4)
      valid = ALL(SHAPE(mixed_covector) == [na, no, nk, ns]) .AND. &
              ALL(SHAPE(auxiliary_covector) == [na, na, nk, ns]) .AND. &
              ALL(SHAPE(primary_covector) == [no, no, nk, ns])
      CPASSERT(valid .AND. ns > 0)
      ALLOCATE (auxiliary(na, na, nk), mixed(na, no, nk), orbital(no, no, nk))
      auxiliary(:, :, :) = SUM(auxiliary_covector, DIM=4)
      mixed(:, :, :) = SUM(mixed_covector, DIM=4)
      orbital(:, :, :) = SUM(primary_covector, DIM=4)
      CALL qs_upf_inverse_charge_metric_contract(qs_env, admm%upf_inverse, primary, admm%upf_mixed_metric, &
                                                 mixed, atom_response, cell_response, mixed_response)
      CALL qs_upf_inverse_charge_metric_contract(qs_env, admm%upf_inverse, admm%upf_inverse, admm%upf_aux_metric, &
                                                 auxiliary, atom_response, cell_response, auxiliary_response)
      orbital_response = 0.0_dp
      IF (ANY(orbital /= CMPLX(0.0_dp, 0.0_dp, dp))) THEN
         IF (ASSOCIATED(primary%coulomb)) THEN
            CALL qs_upf_inverse_coulomb_metric_response(qs_env, primary, primary%coulomb, orbital, &
                                                        atom_response, cell_response, orbital_response)
         ELSE
            CALL primary_grid_metric_response(qs_env, primary, orbital, atom_response, cell_response, orbital_response)
         END IF
      END IF
      response = mixed_response + auxiliary_response + orbital_response
   END SUBROUTINE qs_upf_inverse_aux_metric_response

! **************************************************************************************************
!> \brief Contract the primary solver's original grid metric before a complete metric is installed.
!>        Match its native Gaussian reference and inverse/Q remainder. Grid and Gaussian source
!>        measures enter once. A later complete density metric uses its own service above.
!> \param qs_env Native geometry and primary Gaussian task owner.
!> \param state Primary inverse representation without an installed complete Coulomb metric.
!> \param matrix Hermitian valence metric covector, without mesh weights, summed over spins.
!> \param atom_response Total atomic motion.
!> \param cell_response Lattice response at fixed reduced k.
!> \param response Directional metric contraction, without a Pulay minus sign.
! **************************************************************************************************
   SUBROUTINE primary_grid_metric_response(qs_env, state, matrix, atom_response, cell_response, response)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_upf_inverse_scf_type), INTENT(INOUT), &
         TARGET                                          :: state
      COMPLEX(dp), INTENT(IN)                            :: matrix(:, :, :)
      REAL(dp), INTENT(IN)                               :: atom_response(:, :), cell_response(3, 3)
      REAL(dp), INTENT(OUT)                              :: response

      COMPLEX(dp), ALLOCATABLE                           :: derivative(:, :), metric(:, :)
      INTEGER                                            :: info, k, nk
      REAL(dp)                                           :: charge, inverse_cell(3, 3), &
                                                            strain(3, 3), trace_strain, &
                                                            virial(3, 3), weight
      REAL(dp), ALLOCATABLE                              :: forces(:, :), point_response(:, :)
      TYPE(pw_c1d_gs_type)                               :: density_g
      TYPE(pw_r3d_rs_type)                               :: density, one
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_upf_inverse_grid_type), POINTER            :: frame

      nk = 1
      frame => state%frame
      IF (ASSOCIATED(state%bloch%kpoints)) THEN
         nk = SIZE(state%bloch%weights)
         frame => state%bloch%frames(1)
      END IF
      CPASSERT(SIZE(matrix, 3) == nk)
      CPASSERT(ALL(SHAPE(atom_response) == SHAPE(frame%ao%positions)))
      inverse_cell = frame%ao%cell
      CALL invmat(inverse_cell, info)
      CPASSERT(info == 0)
      strain = MATMUL(cell_response, inverse_cell)
      trace_strain = strain(1, 1) + strain(2, 2) + strain(3, 3)
      response = 0.0_dp
      IF (state%native_gaussian) THEN
         CALL get_qs_env(qs_env, ks_env=ks_env)
         CALL one%create(state%source_grid)
         CALL density%create(state%source_grid)
         CALL density_g%create(state%source_grid)
         one%array = 1.0_dp
         CALL qs_upf_inverse_background_collocate(state%gaussian, ks_env, matrix, density, density_g, total_rho=charge)
         ALLOCATE (forces, MOLD=atom_response)
         CALL qs_upf_inverse_background_geometry(state%gaussian, qs_env, matrix, one, forces, virial)
         ! Native total_rho is the negative electronic density integral.
         response = -trace_strain*charge - SUM(virial*strain) - &
                    SUM(forces*(atom_response - MATMUL(strain, frame%ao%positions)))
         CALL one%release()
         CALL density%release()
         CALL density_g%release()
      END IF
      DO k = 1, nk
         weight = 1.0_dp
         IF (ASSOCIATED(state%bloch%kpoints)) THEN
            frame => state%bloch%frames(k)
            weight = state%bloch%weights(k)
         END IF
         ALLOCATE (point_response, SOURCE=MATMUL(strain, frame%ao%points))
         CALL qs_upf_inverse_grid_metric_dense(frame, frame%ao%radial_q, metric, &
                                               orbital_fields=.NOT. state%native_gaussian, response=derivative, &
                                               atom_response=atom_response, cell_response=cell_response, &
                                               reduced_k_response=[0.0_dp, 0.0_dp, 0.0_dp], point_response=point_response, &
                                               weights_response=trace_strain*frame%ao%weights)
         response = response + weight*REAL(SUM(CONJG(matrix(:, :, k))*derivative), dp)
         DEALLOCATE (point_response)
      END DO
   END SUBROUTINE primary_grid_metric_response
END MODULE qs_upf_inverse_aux_operator
