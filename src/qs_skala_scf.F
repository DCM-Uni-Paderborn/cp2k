!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
!> \brief Evaluate fixed-occupation native orbital trials and physical Skala roots.
MODULE qs_skala_scf
   USE cell_types,                      ONLY: cell_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_copy,&
                                              dbcsr_p_type,&
                                              dbcsr_set
   USE cp_dbcsr_operations,             ONLY: dbcsr_allocate_matrix_set
   USE cp_external_control,             ONLY: external_control
   USE cp_fm_types,                     ONLY: cp_fm_set_all
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE input_constants,                 ONLY: ot_precond_none,&
                                              smear_fermi_dirac
   USE input_section_types,             ONLY: section_get_ival,&
                                              section_get_rval,&
                                              section_vals_get_subs_vals,&
                                              section_vals_get_subs_vals2,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE kpoint_methods,                  ONLY: kpoint_star_density,&
                                              kpoint_star_operator,&
                                              kpoint_star_orbitals
   USE kpoint_types,                    ONLY: kpoint_env_type,&
                                              kpoint_type
   USE message_passing,                 ONLY: mp_comm_type,&
                                              mp_para_env_type
   USE particle_types,                  ONLY: particle_type
   USE paw_xc,                          ONLY: paw_xc_point_adjoint
   USE pw_grid_types,                   ONLY: pw_grid_type
   USE pw_types,                        ONLY: pw_r3d_rs_type
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_ks_methods,                   ONLY: qs_ks_update_qs_env
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_did_change,&
                                              qs_ks_env_type,&
                                              set_ks_env
   USE qs_mo_types,                     ONLY: mo_set_type
   USE qs_native_grid_interpolation,    ONLY: native_grid_batch_adjoint
   USE qs_rho_methods,                  ONLY: qs_rho_update_rho
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_set,&
                                              qs_rho_type
   USE qs_skala_orbitals,               ONLY: qs_skala_orbitals_canonicalize,&
                                              qs_skala_orbitals_create,&
                                              qs_skala_orbitals_evaluate,&
                                              qs_skala_orbitals_install_gamma,&
                                              qs_skala_orbitals_install_kpoint,&
                                              qs_skala_orbitals_release,&
                                              qs_skala_orbitals_type
   USE qs_upf_inverse,                  ONLY: qs_upf_inverse_evaluate
   USE qs_upf_inverse_augmentation,     ONLY: qs_upf_inverse_augmentation_integrate
   USE qs_upf_inverse_background,       ONLY: qs_upf_inverse_background_integrate
   USE qs_upf_inverse_density,          ONLY: qs_upf_inverse_kinetic_operator
   USE qs_upf_inverse_grid,             ONLY: qs_upf_inverse_grid_integrate_dense,&
                                              qs_upf_inverse_matrix_gather
   USE qs_upf_inverse_paw,              ONLY: qs_upf_inverse_paw_operator
   USE qs_upf_inverse_scf_response,     ONLY: qs_upf_inverse_scf_density_tangent_available,&
                                              qs_upf_inverse_scf_electronic_response
   USE qs_upf_inverse_types,            ONLY: qs_upf_inverse_augmentation_type,&
                                              qs_upf_inverse_scf_type,&
                                              qs_upf_inverse_type
   USE scf_control_types,               ONLY: scf_control_type
   USE skala_density_constraints,       ONLY: &
        skala_background_domain, skala_density_constraints_balance, &
        skala_density_constraints_clear_selection, skala_density_constraints_count, &
        skala_density_constraints_get, skala_density_constraints_matches, &
        skala_density_constraints_reset, skala_density_constraints_set, &
        skala_density_constraints_type, skala_density_constraints_update_weights, &
        skala_density_crossings, skala_density_descriptor, skala_density_jumps_begin, &
        skala_density_jumps_finish, skala_density_root, skala_density_samples_begin, &
        skala_density_samples_end, skala_density_samples_generation, skala_density_samples_get, &
        skala_fft_domain, skala_kinetic_descriptor, skala_last_sample_domain, skala_paw_ae_domain, &
        skala_paw_ps_domain
   USE skala_gpw_functional,            ONLY: native_skala_spin_floors
   USE skala_orbital_constraints,       ONLY: &
        skala_constrained_stationarity, skala_constraint_chart_create, &
        skala_constraint_chart_restore, skala_constraint_chart_type, skala_cusp_admissibility, &
        skala_orbital_chart_density_response, skala_orbital_chart_derivative, &
        skala_orbital_chart_operator_response, skala_orbital_chart_precondition, &
        skala_orbital_chart_response_inverse, skala_orbital_chart_scales, skala_orbital_chart_size
   USE smearing_utils,                  ONLY: Smearkp
   USE upf_inverse_quadrature,          ONLY: upf_inverse_semilocal_pair
   USE xc_derivatives,                  ONLY: xc_functionals_floor_jumps_available,&
                                              xc_functionals_get_needs,&
                                              xc_functionals_spin_floors
   USE xc_input_constants,              ONLY: xc_rho_no_smooth
   USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_type
   USE xc_util,                         ONLY: xc_pw_smooth
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   ! Atomic energy unit for temperature-independent frozen-Hamiltonian occupation curvature.
   REAL(dp), PARAMETER :: thermal_energy_unit = 1.0_dp
   PUBLIC :: qs_skala_scf_type, qs_skala_scf_create, qs_skala_scf_release, qs_skala_scf_size
   PUBLIC :: qs_skala_scf_evaluate, qs_skala_scf_inspect, qs_skala_scf_value, qs_skala_scf_values, qs_skala_scf_root
   PUBLIC :: qs_skala_scf_point_derivative
   PUBLIC :: qs_skala_scf_crossings, qs_skala_scf_crossing_roots
   PUBLIC :: qs_skala_scf_choose_root, qs_skala_scf_boundary_root
   PUBLIC :: qs_skala_scf_constraints, qs_skala_scf_stationarity
   PUBLIC :: qs_skala_scf_crossing_type
   PUBLIC :: qs_skala_scf_select, qs_skala_scf_selected, qs_skala_scf_limiting, qs_skala_scf_balance
   PUBLIC :: qs_skala_scf_admissibility, qs_skala_scf_classify
   PUBLIC :: qs_skala_scf_release_cusp
   PUBLIC :: qs_skala_scf_accept
   PUBLIC :: qs_skala_scf_coordinate_scales, qs_skala_scf_fermi_target
   PUBLIC :: qs_skala_scf_density_response, qs_skala_scf_operator_response
   PUBLIC :: qs_skala_scf_curvature, qs_skala_scf_lagrangian_curvature
   INTEGER, PARAMETER, PUBLIC :: qs_skala_positive_density_sum = -2

   INTEGER(KIND=int_8), PARAMETER :: max_point_operator_bytes = 64_int_8*1024*1024
   INTEGER, PARAMETER :: initial_point_operator_entries = 32

   TYPE physical_point_operator
      INTEGER :: domain = 0, point(3) = 0, descriptor = 0
      COMPLEX(dp), ALLOCATABLE :: matrix(:, :, :)
   END TYPE physical_point_operator

   TYPE qs_skala_scf_type
      PRIVATE
      TYPE(qs_environment_type), POINTER :: qs_env => NULL()
      TYPE(qs_ks_env_type), POINTER :: ks_env => NULL()
      TYPE(qs_rho_type), POINTER :: rho => NULL()
      TYPE(qs_energy_type), POINTER :: energy => NULL()
      TYPE(qs_upf_inverse_scf_type), POINTER :: inverse => NULL()
      TYPE(skala_density_constraints_type), POINTER :: constraints => NULL()
      TYPE(kpoint_type), POINTER :: kpoints => NULL()
      TYPE(mp_para_env_type), POINTER :: para_env => NULL()
      TYPE(qs_skala_orbitals_type), ALLOCATABLE :: orbitals(:, :)
      TYPE(physical_point_operator), ALLOCATABLE :: point_operators(:)
      INTEGER(KIND=int_8) :: point_operator_bytes = 0_int_8
      INTEGER, ALLOCATABLE :: dimensions(:, :), offsets(:, :)
      INTEGER, ALLOCATABLE :: band_dimensions(:, :), band_offsets(:, :), band_groups(:), band_spins(:)
      REAL(dp), ALLOCATABLE :: band_reference(:), band_energies(:), band_weights(:), band_occupations(:), &
                               band_response(:), band_scales(:), thermal_targets(:), thermal_mu(:), thermal_spin_entropy(:)
      INTEGER :: orbital_size = 0
      REAL(dp) :: temperature = 0.0_dp, maxocc = 0.0_dp, entropy = 0.0_dp
      LOGICAL :: thermal = .FALSE.
      REAL(dp), ALLOCATABLE :: coordinate(:), positions(:, :)
      REAL(dp) :: hmat(3, 3) = 0.0_dp
      INTEGER :: sample_generation = -1, sample_density_evaluation = -1
      LOGICAL :: bloch = .FALSE., evaluated = .FALSE., operator_current = .FALSE.
   END TYPE qs_skala_scf_type

   TYPE qs_skala_scf_crossing_type
      TYPE(skala_constraint_chart_type) :: curve
      INTEGER, ALLOCATABLE :: points(:, :), domains(:), spins(:), descriptors(:)
      REAL(dp), ALLOCATABLE :: first(:), last(:), values(:, :), targets(:), root(:)
      LOGICAL, ALLOCATABLE :: selected(:)
      LOGICAL :: encountered = .FALSE., resolved = .FALSE.
   END TYPE qs_skala_scf_crossing_type

   TYPE physical_root_context
      TYPE(qs_skala_scf_type), POINTER :: state => NULL()
      REAL(dp), ALLOCATABLE :: first(:), direction(:)
      INTEGER :: domain, point(3), spin, descriptor
   END TYPE physical_root_context

   TYPE physical_domain_samples
      INTEGER, ALLOCATABLE :: points(:, :)
      REAL(dp), ALLOCATABLE :: density(:, :), kinetic(:, :)
      LOGICAL :: available = .FALSE.
   END TYPE physical_domain_samples

   TYPE physical_constraints_context
      TYPE(qs_skala_scf_type), POINTER :: state => NULL()
      TYPE(qs_skala_scf_crossing_type), POINTER :: crossings => NULL()
      INTEGER, ALLOCATABLE :: points(:, :), domains(:), spins(:), descriptors(:)
      REAL(dp), ALLOCATABLE :: targets(:), density_floors(:)
      INTEGER :: iteration_unit = -1, iteration_offset = 0
      REAL(dp) :: gradient_tolerance = 0.0_dp
      REAL(dp) :: max_coordinate_norm = 0.0_dp
      LOGICAL :: rechart = .FALSE.
      LOGICAL :: check_external = .FALSE., admit_crossings = .FALSE., external_stopped = .FALSE.
      LOGICAL :: kinetic_boundaries = .TRUE.
   END TYPE physical_constraints_context

CONTAINS

! **************************************************************************************************
!> \brief Pull one unweighted contribution operator back to its native representative.
!> \param state Borrowed native orbital state.
!> \param operators AO operators indexed by the complete contribution map.
!> \param ik Native irreducible k-point index.
!> \param matrix Unweighted native AO operator.
! **************************************************************************************************
   SUBROUTINE native_bloch_operator(state, operators, ik, matrix)
      TYPE(qs_skala_scf_type), INTENT(IN)                :: state
      COMPLEX(dp), INTENT(IN)                            :: operators(:, :, :)
      INTEGER, INTENT(IN)                                :: ik
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: matrix(:, :)

      COMPLEX(dp), ALLOCATABLE                           :: reduced(:, :, :)

      CPASSERT(state%bloch)
      CALL kpoint_star_operator(state%inverse%bloch%star, operators, reduced, representative=ik)
      matrix = reduced(:, :, 1)
   END SUBROUTINE native_bloch_operator

! **************************************************************************************************
!> \brief Return spin-resolved contribution covectors to the native orbital coordinates.
!> \param state Borrowed native orbital state.
!> \param operators Owned unweighted (AO,AO,contribution,spin) matrices, replaced by native representatives.
! **************************************************************************************************
   SUBROUTINE reduce_bloch_operators(state, operators)
      TYPE(qs_skala_scf_type), INTENT(IN)                :: state
      COMPLEX(dp), ALLOCATABLE, INTENT(INOUT)            :: operators(:, :, :, :)

      COMPLEX(dp), ALLOCATABLE                           :: block(:, :, :), reduced(:, :, :, :)
      INTEGER                                            :: spin

      CPASSERT(state%bloch)
      IF (ALL(state%inverse%bloch%star%operation == 0)) RETURN
      ALLOCATE (reduced(SIZE(operators, 1), SIZE(operators, 2), state%kpoints%nkp, SIZE(operators, 4)))
      DO spin = 1, SIZE(operators, 4)
         CALL kpoint_star_operator(state%inverse%bloch%star, operators(:, :, :, spin), block)
         reduced(:, :, :, spin) = block
      END DO
      CALL MOVE_ALLOC(reduced, operators)
   END SUBROUTINE reduce_bloch_operators

! **************************************************************************************************
!> \brief Borrow the native QS state and complete each occupied block in its original charge metric.
!>        Collective over QS and each owning k group. Zero-weight and empty-spin blocks have no
!>        coordinates. Geometry, occupations and k weights stay fixed throughout this state.
!> \param state Owned orbital charts, borrowing a QS environment that must outlive this state.
!> \param qs_env Native USPP or PAW environment with initialized MOs and normalization operators.
!> \param valid Whether every active block has a positive complete charge metric.
!> \param thermal Include the original finite-temperature occupation variables when enabled.
! **************************************************************************************************
   SUBROUTINE qs_skala_scf_create(state, qs_env, valid, thermal)
      TYPE(qs_skala_scf_type), INTENT(INOUT)             :: state
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(OUT)                               :: valid
      LOGICAL, INTENT(IN), OPTIONAL                      :: thermal

      COMPLEX(dp), ALLOCATABLE                           :: hamiltonian(:, :), metric(:, :)
      INTEGER                                            :: failed, i, ik, local, nk, nl, ns, &
                                                            offset, spin
      LOGICAL                                            :: block_valid, complex_orbitals, owner
      REAL(dp)                                           :: weight
      REAL(dp), ALLOCATABLE                              :: diagonal(:), zero(:)
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_p_type), POINTER                        :: matrices(:, :), metrics(:)
      TYPE(kpoint_env_type), POINTER                     :: kp
      TYPE(mo_set_type), POINTER                         :: mos(:), native
      TYPE(particle_type), POINTER                       :: particles(:)
      TYPE(scf_control_type), POINTER                    :: control

      CALL qs_skala_scf_release(state)
      IF (PRESENT(thermal)) state%thermal = thermal
      CPASSERT(ASSOCIATED(qs_env))
      state%qs_env => qs_env
      CALL get_qs_env(qs_env, ks_env=state%ks_env, rho=state%rho, energy=state%energy, &
                      kpoints=state%kpoints, mos=mos, para_env=state%para_env, cell=cell, particle_set=particles, &
                      scf_control=control)
      CALL get_ks_env(state%ks_env, upf_inverse=state%inverse, skala_constraints=state%constraints)
      valid = ASSOCIATED(state%inverse) .AND. ASSOCIATED(state%constraints)
      IF (valid) valid = state%inverse%enabled
      IF (.NOT. valid) THEN
         CALL qs_skala_scf_release(state)
         RETURN
      END IF
      state%bloch = ASSOCIATED(state%inverse%bloch%kpoints)
      nk = 1
      nl = 1
      ns = SIZE(mos)
      owner = state%para_env%mepos == 0
      complex_orbitals = .FALSE.
      IF (state%bloch) THEN
         nk = SIZE(state%kpoints%wkp)
         nl = SIZE(state%kpoints%kp_env)
         ns = SIZE(state%kpoints%kp_env(1)%kpoint_env%mos)
         owner = state%kpoints%para_env_kp%mepos == 0
         complex_orbitals = .NOT. state%kpoints%use_real_wfn
      ELSE
         CALL get_ks_env(state%ks_env, matrix_metric=metrics)
         IF (state%thermal) CALL get_ks_env(state%ks_env, matrix_ks_kp=matrices)
      END IF
      IF (state%thermal) THEN
         CPASSERT(control%smear%electronic_temperature > 0.0_dp)
      END IF
      ALLOCATE (state%orbitals(nl, ns), state%dimensions(nk, ns), state%offsets(nk, ns), &
                state%positions(3, SIZE(particles)))
      state%hmat(:, :) = cell%hmat
      DO i = 1, SIZE(particles)
         state%positions(:, i) = particles(i)%r
      END DO
      state%dimensions = 0
      failed = 0
      DO spin = 1, ns
         DO local = 1, nl
            ik = 1
            weight = 1.0_dp
            native => mos(spin)
            IF (state%bloch) THEN
               kp => state%kpoints%kp_env(local)%kpoint_env
               ik = kp%nkpoint
               native => kp%mos(spin)
               weight = state%kpoints%wkp(ik)
            END IF
            CPASSERT(ieee_is_finite(weight) .AND. weight >= 0.0_dp)
            IF (weight == 0.0_dp .OR. native%nmo == 0) CYCLE
            ! Retained virtuals in an empty spin still require matching canonical eigenvalues.
            IF (state%bloch) THEN
               CALL native_bloch_operator(state, state%inverse%metric_k, ik, metric)
            ELSE
               CALL qs_upf_inverse_matrix_gather(state%inverse%frame%ao, metrics(1)%matrix, metric)
            END IF
            CALL qs_skala_orbitals_create(state%orbitals(local, spin), native, metric, &
                                          complex_orbitals, weight, block_valid, variable_occupations=state%thermal)
            IF (.NOT. block_valid) failed = 1
            IF (block_valid .AND. state%thermal) THEN
               IF (state%bloch) THEN
                  CALL native_bloch_operator(state, state%inverse%ks_k(:, :, :, spin), ik, hamiltonian)
               ELSE
                  CALL qs_upf_inverse_matrix_gather(state%inverse%frame%ao, matrices(spin, 1)%matrix, hamiltonian)
               END IF
               diagonal = REAL(SUM(CONJG(state%orbitals(local, spin)%coefficients)* &
                                   MATMUL(hamiltonian, state%orbitals(local, spin)%coefficients), DIM=1), dp)
               CALL skala_orbital_chart_precondition(state%orbitals(local, spin)%chart, diagonal, &
                                                     control%smear%electronic_temperature)
               ALLOCATE (zero(skala_orbital_chart_size(state%orbitals(local, spin)%chart)), SOURCE=0.0_dp)
               CALL qs_skala_orbitals_evaluate(state%orbitals(local, spin), zero)
               DEALLOCATE (zero)
            END IF
            IF (block_valid .AND. owner) THEN
               state%dimensions(ik, spin) = skala_orbital_chart_size(state%orbitals(local, spin)%chart)
            END IF
         END DO
      END DO
      CALL state%para_env%max(failed)
      valid = failed == 0
      IF (.NOT. valid) THEN
         CALL qs_skala_scf_release(state)
         RETURN
      END IF
      CALL state%para_env%sum(state%dimensions)
      offset = 0
      DO spin = 1, ns
         DO ik = 1, nk
            state%offsets(ik, spin) = offset
            offset = offset + state%dimensions(ik, spin)
         END DO
      END DO
      state%orbital_size = offset
      IF (state%thermal) CALL initialize_thermal_state(state, offset)
      ALLOCATE (state%coordinate(offset))
      state%coordinate = 0.0_dp
      IF (state%thermal) THEN
         CALL update_thermal_state(state, state%coordinate(state%orbital_size + 1:))
         CALL initialize_thermal_scales(state)
      END IF
   END SUBROUTINE qs_skala_scf_create

! **************************************************************************************************
!> \brief Release owned charts and buffers without releasing or resetting the borrowed native QS state.
!> \param state ...
! **************************************************************************************************
   SUBROUTINE qs_skala_scf_release(state)
      TYPE(qs_skala_scf_type), INTENT(INOUT)             :: state

      INTEGER                                            :: local, spin

      IF (ALLOCATED(state%orbitals)) THEN
         DO spin = 1, SIZE(state%orbitals, 2)
            DO local = 1, SIZE(state%orbitals, 1)
               CALL qs_skala_orbitals_release(state%orbitals(local, spin))
            END DO
         END DO
         DEALLOCATE (state%orbitals, state%dimensions, state%offsets, state%positions)
      END IF
      IF (ALLOCATED(state%band_dimensions)) THEN
         DEALLOCATE (state%band_dimensions, state%band_offsets, state%band_groups, state%band_spins, &
                     state%band_reference, state%band_energies, state%band_weights, state%band_occupations, &
                     state%band_response, state%band_scales, state%thermal_targets, state%thermal_mu, state%thermal_spin_entropy)
      END IF
      state%thermal = .FALSE.
      state%orbital_size = 0
      state%temperature = 0.0_dp
      state%maxocc = 0.0_dp
      state%entropy = 0.0_dp
      IF (ALLOCATED(state%coordinate)) DEALLOCATE (state%coordinate)
      IF (ALLOCATED(state%point_operators)) DEALLOCATE (state%point_operators)
      state%point_operator_bytes = 0_int_8
      NULLIFY (state%qs_env, state%ks_env, state%rho, state%energy, state%inverse, state%constraints, &
               state%kpoints, state%para_env)
      state%hmat = 0.0_dp
      state%bloch = .FALSE.
      state%evaluated = .FALSE.
      state%operator_current = .FALSE.
      state%sample_generation = -1
      state%sample_density_evaluation = -1
   END SUBROUTINE qs_skala_scf_release

! **************************************************************************************************
!> \brief Number of global independent fixed-occupation orbital coordinates.
!> \param state ...
!> \return ...
! **************************************************************************************************
   FUNCTION qs_skala_scf_size(state) RESULT(n)
      TYPE(qs_skala_scf_type), INTENT(IN)                :: state
      INTEGER                                            :: n

      CPASSERT(ALLOCATED(state%coordinate))
      n = SIZE(state%coordinate)
   END FUNCTION qs_skala_scf_size

! **************************************************************************************************
!> \brief Read the constant map to unscaled Cayley and dimensionless auxiliary-energy coordinates.
!>        This map is frozen at construction, including the positive thermal response factors.
!> \param state Native chart with its original particle-number groups.
!> \param scales Per-coordinate multiplier to the unscaled reference variables.
! **************************************************************************************************
   SUBROUTINE qs_skala_scf_coordinate_scales(state, scales)
      TYPE(qs_skala_scf_type), INTENT(IN)                :: state
      REAL(dp), INTENT(OUT)                              :: scales(:)

      INTEGER                                            :: first, ik, last, local, spin
      LOGICAL                                            :: owner

      CPASSERT(SIZE(scales) == qs_skala_scf_size(state))
      scales = 0.0_dp
      owner = state%para_env%mepos == 0
      IF (state%bloch) owner = state%kpoints%para_env_kp%mepos == 0
      DO spin = 1, SIZE(state%orbitals, 2)
         DO local = 1, SIZE(state%orbitals, 1)
            IF (.NOT. ASSOCIATED(state%orbitals(local, spin)%native) .OR. .NOT. owner) CYCLE
            ik = 1
            IF (state%bloch) ik = state%kpoints%kp_env(local)%kpoint_env%nkpoint
            first = state%offsets(ik, spin) + 1
            last = first + state%dimensions(ik, spin) - 1
            CALL skala_orbital_chart_scales(state%orbitals(local, spin)%chart, scales(first:last))
         END DO
      END DO
      CALL state%para_env%sum(scales)
      IF (state%thermal) scales(state%orbital_size + 1:) = state%band_scales/state%temperature
      CPASSERT(ALL(ieee_is_finite(scales)) .AND. ALL(scales > 0.0_dp))
   END SUBROUTINE qs_skala_scf_coordinate_scales

! **************************************************************************************************
!> \brief Install all owning orbital blocks, rebuild rho/tau and evaluate the complete native energy.
!>        Samples are the actual local model inputs of this trial, including atomic/background maps.
!>        Gradients are globally ordered and include each k weight once. Reevaluate the accepted
!>        coordinate after rejecting a trial; native MOs, P and fields then describe the same state.
!> \param state ...
!> \param coordinate Replicated global coordinates, in spin then global k-point order.
!> \param energy Complete native energy at fixed occupations.
!> \param gradient Optional derivative in the supplied Cayley coordinates.
!> \param physical_gradient Optional derivative in the orthonormal physical orbital tangent.
!> \param sample_inputs Retain physical point fields, default true. Disabled fields are cleared.
! **************************************************************************************************
   SUBROUTINE qs_skala_scf_evaluate(state, coordinate, energy, gradient, physical_gradient, sample_inputs)
      TYPE(qs_skala_scf_type), INTENT(INOUT)             :: state
      REAL(dp), INTENT(IN)                               :: coordinate(:)
      REAL(dp), INTENT(OUT)                              :: energy
      REAL(dp), INTENT(OUT), OPTIONAL                    :: gradient(:), physical_gradient(:)
      LOGICAL, INTENT(IN), OPTIONAL                      :: sample_inputs

      LOGICAL                                            :: need_gradient, record_inputs

      CALL install_native_trial(state, coordinate)
      record_inputs = .TRUE.
      IF (PRESENT(sample_inputs)) record_inputs = sample_inputs
      CALL skala_density_samples_begin(state%constraints)
      IF (.NOT. record_inputs) CALL skala_density_samples_end(state%constraints)
      CALL qs_ks_did_change(state%ks_env, rho_changed=.TRUE.)
      need_gradient = PRESENT(gradient) .OR. PRESENT(physical_gradient)
      CALL qs_ks_update_qs_env(state%qs_env, calculate_forces=.FALSE., &
                               just_energy=.NOT. need_gradient, print_active=.FALSE.)
      CALL skala_density_samples_end(state%constraints)
      IF (record_inputs) THEN
         CALL retain_native_samples(state)
      ELSE
         state%sample_generation = -1
         state%sample_density_evaluation = -1
      END IF
      energy = state%energy%total
      state%operator_current = need_gradient
      IF (need_gradient) CALL native_orbital_covectors(state, gradient, physical_gradient)
      state%coordinate(:) = coordinate
      state%evaluated = .TRUE.
   END SUBROUTINE qs_skala_scf_evaluate

! **************************************************************************************************
!> \brief Install the original normalized orbital trial and rebuild its native rho/tau sources.
!> \param state ...
!> \param coordinate ...
! **************************************************************************************************
   SUBROUTINE install_native_trial(state, coordinate)
      TYPE(qs_skala_scf_type), INTENT(INOUT)             :: state
      REAL(dp), INTENT(IN)                               :: coordinate(:)

      INTEGER                                            :: bfirst, blast, first, group, ik, last, &
                                                            local, spin
      LOGICAL                                            :: valid
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_p_type), POINTER                        :: densities(:, :)
      TYPE(particle_type), POINTER                       :: particles(:)

      state%sample_generation = -1
      state%sample_density_evaluation = -1
      state%operator_current = .FALSE.
      CPASSERT(ASSOCIATED(state%qs_env))
      state%inverse%orbital_factors%density_evaluation = -1
      CPASSERT(SIZE(coordinate) == qs_skala_scf_size(state))
      CPASSERT(ALL(ieee_is_finite(coordinate)))
      CALL get_qs_env(state%qs_env, cell=cell, particle_set=particles)
      CPASSERT(ALL(cell%hmat == state%hmat))
      CPASSERT(SIZE(particles) == SIZE(state%positions, 2))
      DO local = 1, SIZE(particles)
         CPASSERT(ALL(particles(local)%r == state%positions(:, local)))
      END DO
      CALL qs_rho_get(state%rho, rho_ao_kp=densities)
      IF (state%thermal) THEN
         CALL update_thermal_state(state, coordinate(state%orbital_size + 1:))
         ! CP2K stores the negative entropic correction and includes it in energy%total.
         state%energy%kTS = state%entropy
      END IF
      DO spin = 1, SIZE(state%orbitals, 2)
         DO local = 1, SIZE(state%orbitals, 1)
            IF (.NOT. ASSOCIATED(state%orbitals(local, spin)%native)) CYCLE
            valid = ALL(state%orbitals(local, spin)%native%occupation_numbers(: &
                                                                              state%orbitals(local, spin)%native%nmo) == &
                        state%orbitals(local, spin)%occupations(:state%orbitals(local, spin)%native%nmo))
            CPASSERT(valid)
            ik = 1
            IF (state%bloch) THEN
               ik = state%kpoints%kp_env(local)%kpoint_env%nkpoint
               CPASSERT(state%kpoints%wkp(ik) == state%orbitals(local, spin)%weight)
            END IF
            first = state%offsets(ik, spin) + 1
            last = first + state%dimensions(ik, spin) - 1
            IF (state%thermal) THEN
               bfirst = state%band_offsets(ik, spin) + 1
               blast = bfirst + state%band_dimensions(ik, spin) - 1
               group = state%band_groups(bfirst)
               CALL qs_skala_orbitals_evaluate(state%orbitals(local, spin), coordinate(first:last), &
                                               state%band_occupations(bfirst:blast), state%band_energies(bfirst:blast))
               state%orbitals(local, spin)%native%mu = state%thermal_mu(group)
               state%orbitals(local, spin)%native%kTS = state%thermal_spin_entropy(spin)
            ELSE
               CALL qs_skala_orbitals_evaluate(state%orbitals(local, spin), coordinate(first:last))
            END IF
            IF (state%bloch) THEN
               CALL qs_skala_orbitals_install_kpoint(state%orbitals(local, spin), state%kpoints, local, spin)
            ELSE
               CALL qs_skala_orbitals_install_gamma(state%orbitals(local, spin), &
                                                    state%inverse%frame%ao, densities(spin, 1)%matrix)
            END IF
         END DO
      END DO
      CALL qs_rho_set(state%rho, rho_r_valid=.FALSE., rho_g_valid=.FALSE., tau_r_valid=.FALSE., tau_g_valid=.FALSE.)
      CALL qs_rho_update_rho(state%rho, state%qs_env)
      CALL retain_native_orbital_factors(state)
   END SUBROUTINE install_native_trial

! **************************************************************************************************
!> \brief Publish MO factors belonging to the newly installed native density.
!>        Use the same publication after trial installation and canonical acceptance, so the
!>        stabilized point-density evaluation retains the factors of the actual current MOs.
!> \param state Native orbital blocks with their density already rebuilt.
! **************************************************************************************************
   SUBROUTINE retain_native_orbital_factors(state)
      TYPE(qs_skala_scf_type), INTENT(INOUT)             :: state

      COMPLEX(dp), ALLOCATABLE                           :: coefficients(:, :, :)
      INTEGER                                            :: ik, local, nao, nk, nmo, ns, spin
      LOGICAL                                            :: owner, valid
      REAL(dp), ALLOCATABLE                              :: occupations(:, :)

      ! Retain the original factors only after their matching density is installed.
      ! Mixed or externally replaced densities have a different evaluation stamp.
      nao = SIZE(state%orbitals(1, 1)%coefficients, 1)
      nk = SIZE(state%dimensions, 1)
      IF (state%bloch) nk = SIZE(state%inverse%bloch%weights)
      ns = SIZE(state%orbitals, 2)
      IF (.NOT. ALLOCATED(state%inverse%orbital_factors%coefficients)) THEN
         ALLOCATE (state%inverse%orbital_factors%coefficients(nao, nao, nk, ns), &
                   state%inverse%orbital_factors%occupations(nao, nk, ns))
      END IF
      valid = ALL(SHAPE(state%inverse%orbital_factors%coefficients) == [nao, nao, nk, ns])
      CPASSERT(valid)
      state%inverse%orbital_factors%coefficients = CMPLX(0.0_dp, 0.0_dp, dp)
      state%inverse%orbital_factors%occupations = 0.0_dp
      owner = state%para_env%mepos == 0
      IF (state%bloch) owner = state%kpoints%para_env_kp%mepos == 0
      DO spin = 1, ns
         DO local = 1, SIZE(state%orbitals, 1)
            IF (.NOT. ASSOCIATED(state%orbitals(local, spin)%native)) CYCLE
            IF (.NOT. owner) CYCLE
            ik = 1
            IF (state%bloch) ik = state%kpoints%kp_env(local)%kpoint_env%nkpoint
            nmo = state%orbitals(local, spin)%native%nmo
            state%inverse%orbital_factors%coefficients(:, :nmo, ik, spin) = &
               state%orbitals(local, spin)%coefficients(:, :nmo)
            state%inverse%orbital_factors%occupations(:nmo, ik, spin) = &
               state%orbitals(local, spin)%occupations(:nmo)
         END DO
      END DO
      CALL state%para_env%sum(state%inverse%orbital_factors%coefficients)
      CALL state%para_env%sum(state%inverse%orbital_factors%occupations)
      IF (state%bloch) THEN
         IF (ANY(state%inverse%bloch%star%operation /= 0)) THEN
            DO spin = 1, ns
               CALL kpoint_star_orbitals(state%inverse%bloch%star, &
                                         state%inverse%orbital_factors%coefficients(:, :, :SIZE(state%dimensions, 1), spin), &
                                         coefficients)
               state%inverse%orbital_factors%coefficients(:, :, :, spin) = coefficients
               occupations = state%inverse%orbital_factors%occupations(:, :SIZE(state%dimensions, 1), spin)
               state%inverse%orbital_factors%occupations(:, :, spin) = &
                  occupations(:, state%inverse%bloch%star%representative)
            END DO
         END IF
      END IF
      state%inverse%orbital_factors%density_evaluation = state%inverse%density_evaluations
   END SUBROUTINE retain_native_orbital_factors

! **************************************************************************************************
!> \brief Check whether the requested coordinates are already installed in the native state.
!> \param state Native orbital state.
!> \param coordinate Requested normalized orbital coordinates.
!> \return Whether the evaluated native state represents these exact coordinates.
! **************************************************************************************************
   FUNCTION at_native_coordinate(state, coordinate) RESULT(current)
      TYPE(qs_skala_scf_type), INTENT(IN)                :: state
      REAL(dp), INTENT(IN)                               :: coordinate(:)
      LOGICAL                                            :: current

      current = .FALSE.
      IF (.NOT. state%evaluated) RETURN
      IF (SIZE(state%coordinate) /= SIZE(coordinate)) RETURN
      current = ALL(state%coordinate == coordinate)
   END FUNCTION at_native_coordinate

! **************************************************************************************************
!> \brief Sample the original physical fields without evaluating the learned XC model.
!>        The partial energy and KS operator are invalid. A later energy callback rebuilds both.
!> \param state ...
!> \param coordinate ...
! **************************************************************************************************
   SUBROUTINE sample_native_constraints(state, coordinate)
      TYPE(qs_skala_scf_type), INTENT(INOUT)             :: state
      REAL(dp), INTENT(IN)                               :: coordinate(:)

      ! Constraint restoration can query the same installed fields more than once.
      ! Reuse only this complete collection at the identical coordinate and density stamp.
      IF (at_native_coordinate(state, coordinate) .AND. state%sample_generation >= 0) THEN
         IF (state%sample_generation == skala_density_samples_generation(state%constraints) .AND. &
             state%sample_density_evaluation == state%inverse%density_evaluations) THEN
            ! Preserve the field-only callback contract even when no traversal is needed.
            ! Selected derivative weights can change without changing these raw fields.
            state%operator_current = .FALSE.
            CALL qs_ks_did_change(state%ks_env, rho_changed=.TRUE.)
            RETURN
         END IF
      END IF
      IF (.NOT. at_native_coordinate(state, coordinate)) THEN
         CALL install_native_trial(state, coordinate)
      END IF
      CALL skala_density_samples_begin(state%constraints, descriptors_only=.TRUE.)
      CALL qs_ks_did_change(state%ks_env, rho_changed=.TRUE.)
      CALL qs_ks_update_qs_env(state%qs_env, calculate_forces=.FALSE., just_energy=.TRUE., print_active=.FALSE.)
      CALL skala_density_samples_end(state%constraints)
      CALL retain_native_samples(state)
      state%coordinate(:) = coordinate
      state%evaluated = .TRUE.
      state%operator_current = .FALSE.
   END SUBROUTINE sample_native_constraints

! **************************************************************************************************
!> \brief Associate a completed descriptor collection with its installed native density.
!> \param state Native orbital and physical-field state.
! **************************************************************************************************
   SUBROUTINE retain_native_samples(state)
      TYPE(qs_skala_scf_type), INTENT(INOUT)             :: state

      state%sample_generation = skala_density_samples_generation(state%constraints)
      state%sample_density_evaluation = state%inverse%density_evaluations
   END SUBROUTINE retain_native_samples

! **************************************************************************************************
!> \brief Inspect a fresh native orbital chart without changing MOs, density, fields or KS.
!>        The caller supplies an already evaluated native state. Chart construction supplies
!>        only its normalized tangent basis. No orbital installation or XC reevaluation is
!>        performed, since renormalization can change a physical boundary at rounding scale.
!> \param state Fresh chart from qs_skala_scf_create, borrowing the current native state.
!> \param energy Current native total energy.
!> \param physical_gradient Complete current KS covector in the physical orbital tangent.
! **************************************************************************************************
   SUBROUTINE qs_skala_scf_inspect(state, energy, physical_gradient)
      TYPE(qs_skala_scf_type), INTENT(INOUT)             :: state
      REAL(dp), INTENT(OUT)                              :: energy, physical_gradient(:)

      CPASSERT(ASSOCIATED(state%qs_env) .AND. .NOT. state%evaluated)
      CPASSERT(ALL(state%coordinate == 0.0_dp))
      state%operator_current = .TRUE.
      CALL native_orbital_covectors(state, physical_gradient=physical_gradient)
      energy = state%energy%total
      state%evaluated = .TRUE.
   END SUBROUTINE qs_skala_scf_inspect

! **************************************************************************************************
!> \brief Contract the current native KS operator with orbital-coordinate or physical tangents.
!> \param state Current charts and native KS matrices, unchanged by this contraction.
!> \param gradient Optional derivative in the current chart coordinates.
!> \param physical_gradient Optional covector in the orthonormal physical orbital tangent.
! **************************************************************************************************
   SUBROUTINE native_orbital_covectors(state, gradient, physical_gradient)
      TYPE(qs_skala_scf_type), INTENT(INOUT)             :: state
      REAL(dp), INTENT(OUT), OPTIONAL                    :: gradient(:), physical_gradient(:)

      COMPLEX(dp), ALLOCATABLE                           :: hamiltonian(:, :)
      INTEGER                                            :: first, ik, last, local, spin
      LOGICAL                                            :: owner
      REAL(dp), ALLOCATABLE                              :: auxiliary(:), block_gradient(:), &
                                                            expectations(:)
      TYPE(dbcsr_p_type), POINTER                        :: matrices(:, :)

      CPASSERT(state%operator_current)
      IF (state%thermal) THEN
         ALLOCATE (expectations(SIZE(state%band_reference)), SOURCE=0.0_dp)
         ALLOCATE (auxiliary(SIZE(state%band_reference)))
      END IF
      IF (PRESENT(gradient)) THEN
         CPASSERT(SIZE(gradient) == qs_skala_scf_size(state))
         gradient = 0.0_dp
      END IF
      IF (PRESENT(physical_gradient)) THEN
         CPASSERT(SIZE(physical_gradient) == qs_skala_scf_size(state))
         physical_gradient = 0.0_dp
      END IF
      owner = state%para_env%mepos == 0
      IF (state%bloch) THEN
         owner = state%kpoints%para_env_kp%mepos == 0
      ELSE
         CALL get_ks_env(state%ks_env, matrix_ks_kp=matrices)
      END IF
      DO spin = 1, SIZE(state%orbitals, 2)
         DO local = 1, SIZE(state%orbitals, 1)
            IF (.NOT. ASSOCIATED(state%orbitals(local, spin)%native)) CYCLE
            ik = 1
            IF (state%bloch) ik = state%kpoints%kp_env(local)%kpoint_env%nkpoint
            first = state%offsets(ik, spin) + 1
            last = first + state%dimensions(ik, spin) - 1
            IF (state%bloch) THEN
               CALL native_bloch_operator(state, state%inverse%ks_k(:, :, :, spin), ik, hamiltonian)
            ELSE
               CALL qs_upf_inverse_matrix_gather(state%inverse%frame%ao, matrices(spin, 1)%matrix, hamiltonian)
            END IF
            IF (.NOT. owner) CYCLE
            ALLOCATE (block_gradient(last - first + 1))
            IF (state%thermal) CALL thermal_expectations(state, hamiltonian, local, spin, expectations)
            hamiltonian = state%orbitals(local, spin)%weight*hamiltonian
            IF (PRESENT(gradient)) THEN
               CALL skala_orbital_chart_derivative(state%orbitals(local, spin)%chart, hamiltonian, block_gradient)
               gradient(first:last) = block_gradient
            END IF
            IF (PRESENT(physical_gradient)) THEN
               CALL skala_orbital_chart_derivative(state%orbitals(local, spin)%chart, hamiltonian, &
                                                   block_gradient, physical_tangent=.TRUE.)
               physical_gradient(first:last) = block_gradient
            END IF
            DEALLOCATE (block_gradient)
         END DO
      END DO
      IF (PRESENT(gradient)) CALL state%para_env%sum(gradient)
      IF (PRESENT(physical_gradient)) CALL state%para_env%sum(physical_gradient)
      IF (state%thermal) THEN
         CALL state%para_env%sum(expectations)
         IF (PRESENT(gradient)) THEN
            CALL thermal_covector(state, expectations, .TRUE., auxiliary)
            gradient(state%orbital_size + 1:) = auxiliary
         END IF
         IF (PRESENT(physical_gradient)) THEN
            CALL thermal_covector(state, expectations, .TRUE., auxiliary, physical_gradient(state%orbital_size + 1:))
         END IF
      END IF
   END SUBROUTINE native_orbital_covectors

! **************************************************************************************************
!> \brief Differentiate all native AO density blocks in fixed-geometry orbital and Fermi coordinates.
!>        The returned matrices are unweighted and globally replicated, including empty k owners.
!>        Physical occupation changes preserve the original charge separately in each Fermi group.
!> \param state Evaluated original-metric orbital state.
!> \param direction Complete real orbital and auxiliary-band coordinate direction.
!> \param density_direction New Hermitian AO density directions indexed (AO,AO,k,spin).
!>        Bloch outputs use the complete contribution map consumed by the functional.
! **************************************************************************************************
   SUBROUTINE qs_skala_scf_density_response(state, direction, density_direction)
      TYPE(qs_skala_scf_type), INTENT(INOUT)             :: state
      REAL(dp), INTENT(IN)                               :: direction(:)
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: density_direction(:, :, :, :)

      COMPLEX(dp), ALLOCATABLE                           :: block(:, :, :), expanded(:, :, :, :)
      INTEGER                                            :: bfirst, blast, first, ik, last, local, &
                                                            n, nk, ns, spin
      LOGICAL                                            :: owner, valid
      REAL(dp), ALLOCATABLE                              :: block_occupations(:), &
                                                            occupations_direction(:)

      CPASSERT(state%evaluated)
      CPASSERT(SIZE(direction) == qs_skala_scf_size(state))
      CPASSERT(ALL(ieee_is_finite(direction)))
      nk = SIZE(state%dimensions, 1)
      ns = SIZE(state%dimensions, 2)
      IF (state%bloch) THEN
         n = SIZE(state%inverse%metric_k, 1)
      ELSE
         n = SIZE(state%inverse%frame%ao%coefficients, 2)
      END IF
      ALLOCATE (density_direction(n, n, nk, ns), SOURCE=CMPLX(0.0_dp, 0.0_dp, dp))
      IF (state%thermal) CALL thermal_band_direction(state, direction(state%orbital_size + 1:), occupations_direction)
      owner = state%para_env%mepos == 0
      IF (state%bloch) owner = state%kpoints%para_env_kp%mepos == 0
      DO spin = 1, ns
         DO local = 1, SIZE(state%orbitals, 1)
            IF (.NOT. ASSOCIATED(state%orbitals(local, spin)%native)) CYCLE
            IF (.NOT. owner) CYCLE
            valid = .NOT. state%orbitals(local, spin)%canonicalized
            CPASSERT(valid)
            ik = 1
            IF (state%bloch) ik = state%kpoints%kp_env(local)%kpoint_env%nkpoint
            first = state%offsets(ik, spin) + 1
            last = first + state%dimensions(ik, spin) - 1
            IF (state%thermal) THEN
               bfirst = state%band_offsets(ik, spin) + 1
               blast = bfirst + state%band_dimensions(ik, spin) - 1
               ALLOCATE (block_occupations(SIZE(state%orbitals(local, spin)%occupations)), SOURCE=0.0_dp)
               block_occupations(:blast - bfirst + 1) = occupations_direction(bfirst:blast)
               CALL skala_orbital_chart_density_response(state%orbitals(local, spin)%chart, direction(first:last), &
                                                         density_direction(:, :, ik, spin), block_occupations)
               DEALLOCATE (block_occupations)
            ELSE
               CALL skala_orbital_chart_density_response(state%orbitals(local, spin)%chart, direction(first:last), &
                                                         density_direction(:, :, ik, spin))
            END IF
         END DO
      END DO
      CALL state%para_env%sum(density_direction)
      IF (state%bloch) THEN
         IF (ANY(state%inverse%bloch%star%operation /= 0)) THEN
            ALLOCATE (expanded(n, n, SIZE(state%inverse%bloch%weights), ns))
            DO spin = 1, ns
               CALL kpoint_star_density(state%inverse%bloch%star, density_direction(:, :, :, spin), block)
               expanded(:, :, :, spin) = block
            END DO
            CALL MOVE_ALLOC(expanded, density_direction)
         END IF
      END IF
   END SUBROUTINE qs_skala_scf_density_response

! **************************************************************************************************
!> \brief Differentiate a complete Hermitian AO covector in native orbital and Fermi coordinates.
!>        Operators and their responses are unweighted, globally replicated (AO,AO,k,spin) blocks.
!>        The original k weights enter once. Include both Cayley factors, physical occupation
!>        changes, the Fermi chemical-potential response and optional original entropy curvature.
!>        The caller supplies the electronic operator response. Geometry and coordinate scales
!>        remain fixed. Selected constraint-normal responses are separate.
!> \param state Evaluated original-metric orbital state.
!> \param direction Complete real orbital and auxiliary-band coordinate direction.
!> \param operators Unweighted Hermitian AO operators at the current state.
!> \param operators_direction Their unweighted electronic response.
!> \param response Complete coordinate-covector response.
!> \param entropy Include the original finite-temperature entropy response when true.
! **************************************************************************************************
   SUBROUTINE qs_skala_scf_operator_response(state, direction, operators, operators_direction, response, entropy)
      TYPE(qs_skala_scf_type), INTENT(INOUT)             :: state
      REAL(dp), INTENT(IN)                               :: direction(:)
      COMPLEX(dp), INTENT(IN) :: operators(:, :, :, :), operators_direction(:, :, :, :)
      REAL(dp), INTENT(OUT)                              :: response(:)
      LOGICAL, INTENT(IN), OPTIONAL                      :: entropy

      COMPLEX(dp), ALLOCATABLE                           :: density_direction(:, :), &
                                                            orbital_direction(:, :), work(:, :), &
                                                            work_direction(:, :)
      INTEGER                                            :: band, bfirst, blast, first, group, i, &
                                                            ik, last, local, n, nk, ns, spin
      INTEGER, ALLOCATABLE                               :: indices(:)
      LOGICAL                                            :: include_entropy, owner, valid
      REAL(dp)                                           :: mean, mean_direction, weight
      REAL(dp), ALLOCATABLE :: block_occupations(:), energies_direction(:), expectations(:), &
         expectations_direction(:), occupations_direction(:), residual(:), residual_direction(:), &
         weighted_direction(:)

      CPASSERT(state%evaluated)
      CPASSERT(SIZE(direction) == qs_skala_scf_size(state))
      CPASSERT(SIZE(response) == SIZE(direction))
      CPASSERT(ALL(ieee_is_finite(direction)))
      nk = SIZE(state%dimensions, 1)
      ns = SIZE(state%dimensions, 2)
      IF (state%bloch) THEN
         n = SIZE(state%inverse%metric_k, 1)
      ELSE
         n = SIZE(state%inverse%frame%ao%coefficients, 2)
      END IF
      valid = ALL(SHAPE(operators) == [n, n, nk, ns])
      CPASSERT(valid)
      valid = ALL(SHAPE(operators_direction) == SHAPE(operators))
      CPASSERT(valid)
      response = 0.0_dp
      include_entropy = .FALSE.
      IF (PRESENT(entropy)) include_entropy = entropy
      IF (state%thermal) THEN
         CALL thermal_band_direction(state, direction(state%orbital_size + 1:), occupations_direction, &
                                     weighted_direction, energies_direction)
         ALLOCATE (expectations(SIZE(state%band_reference)), expectations_direction(SIZE(state%band_reference)))
         expectations = 0.0_dp
         expectations_direction = 0.0_dp
      END IF
      owner = state%para_env%mepos == 0
      IF (state%bloch) owner = state%kpoints%para_env_kp%mepos == 0
      DO spin = 1, ns
         DO local = 1, SIZE(state%orbitals, 1)
            IF (.NOT. ASSOCIATED(state%orbitals(local, spin)%native)) CYCLE
            IF (.NOT. owner) CYCLE
            valid = .NOT. state%orbitals(local, spin)%canonicalized
            CPASSERT(valid)
            ik = 1
            IF (state%bloch) ik = state%kpoints%kp_env(local)%kpoint_env%nkpoint
            first = state%offsets(ik, spin) + 1
            last = first + state%dimensions(ik, spin) - 1
            weight = state%orbitals(local, spin)%weight
            IF (state%thermal) THEN
               bfirst = state%band_offsets(ik, spin) + 1
               blast = bfirst + state%band_dimensions(ik, spin) - 1
               ALLOCATE (block_occupations(SIZE(state%orbitals(local, spin)%occupations)), SOURCE=0.0_dp)
               block_occupations(:blast - bfirst + 1) = occupations_direction(bfirst:blast)
               ALLOCATE (density_direction(n, n), orbital_direction(n, SIZE(block_occupations)))
               CALL skala_orbital_chart_density_response(state%orbitals(local, spin)%chart, direction(first:last), &
                                                         density_direction, block_occupations, orbital_direction)
               CALL skala_orbital_chart_operator_response(state%orbitals(local, spin)%chart, direction(first:last), &
                                                          weight*operators(:, :, ik, spin), &
                                                          weight*operators_direction(:, :, ik, spin), &
                                                          response(first:last), block_occupations)
               work = MATMUL(operators(:, :, ik, spin), state%orbitals(local, spin)%coefficients)
               work_direction = MATMUL(operators_direction(:, :, ik, spin), state%orbitals(local, spin)%coefficients)
               DO band = bfirst, blast
                  i = band - bfirst + 1
                  expectations(band) = REAL(SUM(CONJG(state%orbitals(local, spin)%coefficients(:, i))*work(:, i)), dp)
                  expectations_direction(band) = 2.0_dp*REAL(SUM(CONJG(orbital_direction(:, i))*work(:, i)), dp) + &
                                           REAL(SUM(CONJG(state%orbitals(local, spin)%coefficients(:, i))*work_direction(:, i)), dp)
               END DO
               DEALLOCATE (block_occupations, density_direction, orbital_direction, work, work_direction)
            ELSE
               CALL skala_orbital_chart_operator_response(state%orbitals(local, spin)%chart, direction(first:last), &
                                                          weight*operators(:, :, ik, spin), &
                                                          weight*operators_direction(:, :, ik, spin), response(first:last))
            END IF
         END DO
      END DO
      CALL state%para_env%sum(response)
      IF (state%thermal) THEN
         CALL state%para_env%sum(expectations)
         CALL state%para_env%sum(expectations_direction)
         residual = expectations
         residual_direction = expectations_direction
         IF (include_entropy) THEN
            residual = residual - state%band_energies
            residual_direction = residual_direction - energies_direction
         END IF
         DO group = 1, SIZE(state%thermal_targets)
            indices = PACK([(i, i=1, SIZE(state%band_reference))], state%band_groups == group)
            weight = SUM(state%band_response(indices))
            IF (weight <= TINY(weight)) CYCLE
            mean = SUM(state%band_response(indices)*residual(indices))/weight
            mean_direction = (SUM(weighted_direction(indices)*(residual(indices) - mean)) + &
                              SUM(state%band_response(indices)*residual_direction(indices)))/weight
            response(state%orbital_size + indices) = -state%band_scales(indices)/state%temperature* &
                                                     (weighted_direction(indices)*(residual(indices) - mean) + &
                                                      state%band_response(indices)*(residual_direction(indices) - mean_direction))
         END DO
      END IF
   END SUBROUTINE qs_skala_scf_operator_response

! **************************************************************************************************
!> \brief Complete native electronic and entropy curvature in the original orbital/Fermi chart.
!>        Use the actual accepted KS operator for orbital-map curvature. Its density response
!>        includes the full native Skala, Hartree and PAW one-center terms. Geometry, charge
!>        metric and coordinate scales stay fixed. Selected surfaces need separate responses.
!> \param state Evaluated native state with a current complete electronic operator.
!> \param direction Joint orbital and auxiliary-band coordinate direction.
!> \param response Complete native free-energy covector response, overwritten.
!> \param available Whether every physical electronic contribution has a supplied response.
!> \param constraint_operator Optional unweighted sum of multiplier-weighted density point operators.
!>        Subtract it from the orbital part for a native density-floor Lagrangian response.
!>        Such a product must subsequently be projected onto the same feasible tangent.
! **************************************************************************************************
   SUBROUTINE qs_skala_scf_curvature(state, direction, response, available, constraint_operator)
      TYPE(qs_skala_scf_type), INTENT(INOUT)             :: state
      REAL(dp), INTENT(IN)                               :: direction(:)
      REAL(dp), INTENT(OUT)                              :: response(:)
      LOGICAL, INTENT(OUT)                               :: available
      COMPLEX(dp), INTENT(IN), OPTIONAL :: constraint_operator(:, :, :, :)

      COMPLEX(dp), ALLOCATABLE :: density_direction(:, :, :, :), gamma(:, :), &
         operators(:, :, :, :), operators_direction(:, :, :, :)
      INTEGER                                            :: blocked, spin
      TYPE(dbcsr_p_type), POINTER                        :: matrices(:, :)

      available = .FALSE.
      response = 0.0_dp
      CPASSERT(state%evaluated .AND. state%operator_current)
      CPASSERT(SIZE(direction) == qs_skala_scf_size(state))
      CPASSERT(SIZE(response) == SIZE(direction))
      blocked = 0
      IF (ANY(state%orbitals%canonicalized)) blocked = 1
      CALL state%para_env%max(blocked)
      IF (blocked /= 0) RETURN
      CALL qs_skala_scf_density_response(state, direction, density_direction)
      CALL qs_upf_inverse_scf_electronic_response(state%qs_env, density_direction, operators_direction, available, &
                                                  density_tangent=PRESENT(constraint_operator))
      IF (.NOT. available) RETURN
      ALLOCATE (operators, MOLD=operators_direction)
      IF (state%bloch) THEN
         operators = state%inverse%ks_k
         CALL reduce_bloch_operators(state, operators)
         CALL reduce_bloch_operators(state, operators_direction)
      ELSE
         CALL get_ks_env(state%ks_env, matrix_ks_kp=matrices)
         DO spin = 1, SIZE(operators, 4)
            CALL qs_upf_inverse_matrix_gather(state%inverse%frame%ao, matrices(spin, 1)%matrix, gamma)
            operators(:, :, 1, spin) = gamma
         END DO
      END IF
      IF (PRESENT(constraint_operator)) THEN
         CPASSERT(ALL(SHAPE(constraint_operator) == SHAPE(operators)))
         operators = operators - constraint_operator
      END IF
      CALL qs_skala_scf_operator_response(state, direction, operators, operators_direction, response, entropy=.TRUE.)
   END SUBROUTINE qs_skala_scf_curvature

! **************************************************************************************************
!> \brief Apply the native density-floor Lagrangian Hessian to a feasible orbital tangent.
!>        Point densities are linear in the original AO matrices at fixed geometry. Their
!>        coordinate curvature therefore follows from the same orbital operator response.
!>        The selected KS operator is retained. The default electronic kernel differs only
!>        by density-normal terms removed by the caller's final tangent projection.
!> \param state Evaluated original state with feasible registered density-floor constraints.
!> \param direction Full feasible coordinate direction, including dependent motion.
!> \param multipliers Original-unit coefficients in registered constraint order.
!> \param response Full Lagrangian Hessian product, before the final tangent projection.
!> \param available Whether the native selected surfaces and complete response are supported.
! **************************************************************************************************
   SUBROUTINE qs_skala_scf_lagrangian_curvature(state, direction, multipliers, response, available)
      TYPE(qs_skala_scf_type), INTENT(INOUT)             :: state
      REAL(dp), INTENT(IN)                               :: direction(:), multipliers(:)
      REAL(dp), INTENT(OUT)                              :: response(:)
      LOGICAL, INTENT(OUT)                               :: available

      COMPLEX(dp), ALLOCATABLE                           :: normal(:, :, :, :), operators(:, :, :)
      INTEGER                                            :: i, ns
      INTEGER, ALLOCATABLE                               :: descriptors(:), domains(:), &
                                                            points(:, :), spins(:)
      LOGICAL                                            :: valid
      REAL(dp), ALLOCATABLE                              :: targets(:), tolerances(:)

      available = .FALSE.
      response = 0.0_dp
      IF (.NOT. qs_upf_inverse_scf_density_tangent_available(state%qs_env)) RETURN
      CALL qs_skala_scf_selected(state, points, domains, spins, descriptors, targets, tolerances, valid)
      IF (.NOT. valid) RETURN
      CPASSERT(SIZE(multipliers) == SIZE(spins))
      CPASSERT(ALL(ieee_is_finite(multipliers)))
      CALL point_operators(state, domains(1), points(:, 1), descriptors(1), operators)
      ns = SIZE(state%orbitals, 2)
      ALLOCATE (normal(SIZE(operators, 1), SIZE(operators, 2), SIZE(operators, 3), ns), &
                SOURCE=CMPLX(0.0_dp, 0.0_dp, dp))
      DO i = 1, SIZE(spins)
         IF (multipliers(i) == 0.0_dp) CYCLE
         CALL point_operators(state, domains(i), points(:, i), descriptors(i), operators)
         IF (ns == 1) THEN
            CPASSERT(spins(i) == 0)
            ! A native RKS sample is one half of its spin-summed AO density.
            normal(:, :, :, 1) = normal(:, :, :, 1) + 0.5_dp*multipliers(i)*operators
         ELSE
            CPASSERT(spins(i) >= 1 .AND. spins(i) <= ns)
            normal(:, :, :, spins(i)) = normal(:, :, :, spins(i)) + multipliers(i)*operators
         END IF
      END DO
      CALL qs_skala_scf_curvature(state, direction, response, available, constraint_operator=normal)
   END SUBROUTINE qs_skala_scf_lagrangian_curvature

! **************************************************************************************************
!> \brief Physical Fermi occupation directions with the original group charge held fixed.
!> \param state ...
!> \param direction Auxiliary-band coordinate direction.
!> \param occupations_direction Physical, unweighted occupation direction.
!> \param weighted_direction Optional response of the original weighted Fermi occupation variance.
!> \param energies_direction Optional auxiliary-band energy direction.
! **************************************************************************************************
   SUBROUTINE thermal_band_direction(state, direction, occupations_direction, weighted_direction, energies_direction)
      TYPE(qs_skala_scf_type), INTENT(IN)                :: state
      REAL(dp), INTENT(IN)                               :: direction(:)
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: occupations_direction(:)
      REAL(dp), ALLOCATABLE, INTENT(OUT), OPTIONAL       :: weighted_direction(:), &
                                                            energies_direction(:)

      INTEGER                                            :: group, i
      INTEGER, ALLOCATABLE                               :: indices(:)
      REAL(dp)                                           :: mean, weight
      REAL(dp), ALLOCATABLE                              :: energy_direction(:)

      CPASSERT(state%thermal)
      CPASSERT(SIZE(direction) == SIZE(state%band_reference))
      energy_direction = state%band_scales*direction
      ALLOCATE (occupations_direction(SIZE(direction)), SOURCE=0.0_dp)
      DO group = 1, SIZE(state%thermal_targets)
         indices = PACK([(i, i=1, SIZE(direction))], state%band_groups == group)
         weight = SUM(state%band_response(indices))
         IF (weight <= TINY(weight)) CYCLE
         mean = SUM(state%band_response(indices)*energy_direction(indices))/weight
         occupations_direction(indices) = -state%band_occupations(indices)* &
                        (1.0_dp - state%band_occupations(indices)/state%maxocc)/state%temperature*(energy_direction(indices) - mean)
      END DO
      IF (PRESENT(weighted_direction)) THEN
         weighted_direction = state%band_weights*(1.0_dp - 2.0_dp*state%band_occupations/state%maxocc)* &
                              occupations_direction
      END IF
      IF (PRESENT(energies_direction)) energies_direction = energy_direction
   END SUBROUTINE thermal_band_direction

! **************************************************************************************************
!> \brief Read one actual physical model input, requiring exactly one distributed owner.
!> \param state ...
!> \param domain FFT, exact atomic, interpolated background, PAW AE or compensated PAW PS.
!> \param point FFT index, (global atom, atomic row, 1), or PAW (global atom, radial node, angular node).
!> \param spin Minus two for the positive-part density sum, minus one for the spin sum,
!>        zero for paired RKS, one or two for UKS.
!> \param descriptor Density or kinetic density.
!> \param value Physical model input, including native core and density processing.
! **************************************************************************************************
   SUBROUTINE qs_skala_scf_value(state, domain, point, spin, descriptor, value)
      TYPE(qs_skala_scf_type), INTENT(IN)                :: state
      INTEGER, INTENT(IN)                                :: domain, point(3), spin, descriptor
      REAL(dp), INTENT(OUT)                              :: value

      REAL(dp)                                           :: field(1)

      CALL qs_skala_scf_values(state, [domain], RESHAPE(point, [3, 1]), [spin], [descriptor], field)
      value = field(1)
   END SUBROUTINE qs_skala_scf_value

! **************************************************************************************************
!> \brief Read physical samples with one domain copy and one collective per output array.
!>        Every sample retains exactly one distributed owner. No AO fields are cached.
!> \param state Evaluated native orbital and physical field state.
!> \param domains Physical evaluation domain per sample.
!> \param points Global FFT or original atomic sample keys, one column per sample.
!> \param spins Physical spin selectors, with the same sum conventions as qs_skala_scf_value.
!> \param descriptors Density or kinetic density per sample.
!> \param values Actual physical inputs, including native core and density processing.
!> \param rows Optional local row cache for the same physical keys and geometry, initialized to -1.
!>        Positive rows are checked against their physical keys. Zero denotes no local owner.
! **************************************************************************************************
   SUBROUTINE qs_skala_scf_values(state, domains, points, spins, descriptors, values, rows)
      TYPE(qs_skala_scf_type), INTENT(IN)                :: state
      INTEGER, INTENT(IN)                                :: domains(:), points(:, :), spins(:), &
                                                            descriptors(:)
      REAL(dp), INTENT(OUT)                              :: values(:)
      INTEGER, INTENT(INOUT), OPTIONAL                   :: rows(:)

      INTEGER                                            :: c, column, domain, first_row, i, &
                                                            last_row, m
      INTEGER, ALLOCATABLE                               :: keys(:, :), owners(:)
      LOGICAL                                            :: available, valid
      REAL(dp), ALLOCATABLE                              :: density(:, :), kinetic(:, :)

      CPASSERT(state%evaluated)
      m = SIZE(values)
      CPASSERT(SIZE(domains) == m .AND. SIZE(spins) == m .AND. SIZE(descriptors) == m)
      CPASSERT(ALL(SHAPE(points) == [3, m]))
      CPASSERT(ALL(domains >= skala_fft_domain .AND. domains <= skala_last_sample_domain))
      valid = ALL(descriptors >= skala_density_descriptor .AND. descriptors <= skala_kinetic_descriptor)
      CPASSERT(valid)
      IF (PRESENT(rows)) THEN
         CPASSERT(SIZE(rows) == m .AND. ALL(rows >= -1))
      END IF
      values = 0.0_dp
      IF (m == 0) RETURN
      ALLOCATE (owners(m), SOURCE=0)
      DO domain = skala_fft_domain, skala_last_sample_domain
         IF (.NOT. ANY(domains == domain)) CYCLE
         CALL skala_density_samples_get(state%constraints, domain, keys, density, kinetic, available)
         IF (.NOT. available) CYCLE
         DO c = 1, m
            IF (domains(c) /= domain) CYCLE
            IF (spins(c) == qs_skala_positive_density_sum) THEN
               valid = descriptors(c) == skala_density_descriptor
               CPASSERT(valid)
            END IF
            column = spins(c)
            IF (SIZE(density, 2) == 1) THEN
               valid = spins(c) == 0 .OR. spins(c) == -1 .OR. spins(c) == qs_skala_positive_density_sum
               CPASSERT(valid)
               column = 1
            END IF
            valid = spins(c) == -1 .OR. spins(c) == qs_skala_positive_density_sum .OR. &
                    (column >= 1 .AND. column <= SIZE(density, 2))
            CPASSERT(valid)
            first_row = 1
            last_row = SIZE(keys, 2)
            IF (PRESENT(rows)) THEN
               IF (rows(c) >= 0) THEN
                  first_row = rows(c)
                  last_row = rows(c)
                  IF (first_row == 0) CYCLE
                  CPASSERT(first_row <= SIZE(keys, 2))
               ELSE
                  rows(c) = 0
               END IF
            END IF
            DO i = first_row, last_row
               IF (.NOT. ALL(keys(:, i) == points(:, c))) CYCLE
               owners(c) = owners(c) + 1
               IF (PRESENT(rows)) rows(c) = i
               IF (descriptors(c) == skala_density_descriptor) THEN
                  IF (spins(c) == qs_skala_positive_density_sum) THEN
                     values(c) = SUM(MAX(density(i, :), 0.0_dp))*2.0_dp/REAL(SIZE(density, 2), dp)
                  ELSE IF (spins(c) == -1) THEN
                     values(c) = SUM(density(i, :))*2.0_dp/REAL(SIZE(density, 2), dp)
                  ELSE
                     values(c) = density(i, column)
                  END IF
               ELSE
                  CPASSERT(SIZE(kinetic, 2) == SIZE(density, 2))
                  IF (spins(c) == -1) THEN
                     values(c) = SUM(kinetic(i, :))*2.0_dp/REAL(SIZE(kinetic, 2), dp)
                  ELSE
                     values(c) = kinetic(i, column)
                  END IF
               END IF
            END DO
         END DO
      END DO
      CALL state%para_env%sum(owners)
      CPASSERT(ALL(owners == 1))
      CALL state%para_env%sum(values)
   END SUBROUTINE qs_skala_scf_values

! **************************************************************************************************
!> \brief Differentiate one actual physical rho/tau input without changing the native KS operator.
!>        Core terms are constant at fixed geometry. The paired-spin factor and each k weight are
!>        applied once. Atomic background inputs use the transpose of their actual interpolation,
!>        without the negative sign belonging to the background energy contribution.
!> \param state Evaluated native orbital state, with current physical model samples.
!> \param domain FFT, exact atomic, interpolated background, PAW AE or compensated PAW PS.
!> \param point FFT index, (global atom, atomic row, 1), or PAW (global atom, radial node, angular node).
!> \param spin Minus two for the positive-part density sum, minus one for the spin sum,
!>        zero for paired RKS, one or two for UKS. At a zero spin density, the
!>        positive-part sum uses the below-threshold derivative.
!> \param descriptor Density or kinetic density.
!> \param derivative Global derivative in the current Cayley coordinates.
!> \param physical_derivative Optional covector in the orthonormal physical orbital tangent.
! **************************************************************************************************
   SUBROUTINE qs_skala_scf_point_derivative(state, domain, point, spin, descriptor, derivative, physical_derivative)
      TYPE(qs_skala_scf_type), INTENT(INOUT)             :: state
      INTEGER, INTENT(IN)                                :: domain, point(3), spin, descriptor
      REAL(dp), INTENT(OUT)                              :: derivative(:)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: physical_derivative(:)

      REAL(dp)                                           :: values(1)
      REAL(dp), ALLOCATABLE                              :: jacobian(:, :), physical(:, :)

      CPASSERT(SIZE(derivative) == qs_skala_scf_size(state))
      ALLOCATE (jacobian(1, SIZE(derivative)))
      IF (PRESENT(physical_derivative)) THEN
         CPASSERT(SIZE(physical_derivative) == SIZE(derivative))
         ALLOCATE (physical(1, SIZE(derivative)))
         CALL sampled_point_derivatives(state, [domain], RESHAPE(point, [3, 1]), [spin], [descriptor], &
                                        values, jacobian, physical)
         physical_derivative(:) = physical(1, :)
      ELSE
         CALL sampled_point_derivatives(state, [domain], RESHAPE(point, [3, 1]), [spin], [descriptor], &
                                        values, jacobian)
      END IF
      derivative(:) = jacobian(1, :)
   END SUBROUTINE qs_skala_scf_point_derivative

! **************************************************************************************************
!> \brief Differentiate sampled physical inputs with one collective per output matrix.
!>        Each orbital coordinate has one owning k-point/spin block. Read each domain once for
!>        the values and once per spin needed by positive-part sums. Preserve the scalar
!>        sample ownership checks, physical interpolation adjoints and original charge metric.
!> \param state Evaluated native orbital and physical field state.
!> \param domains Physical evaluation domains.
!> \param points Global physical sample keys.
!> \param spins Physical spin selectors, including total and positive-part density sums.
!> \param descriptors Density or kinetic density per sample.
!> \param values Actual sampled inputs, including fixed cores.
!> \param jacobian Derivatives in the supplied Cayley coordinates, one sample per row.
!> \param physical_jacobian Optional covectors in the orthonormal physical orbital tangent.
! **************************************************************************************************
   SUBROUTINE sampled_point_derivatives(state, domains, points, spins, descriptors, values, jacobian, physical_jacobian)
      TYPE(qs_skala_scf_type), INTENT(INOUT)             :: state
      INTEGER, INTENT(IN)                                :: domains(:), points(:, :), spins(:), &
                                                            descriptors(:)
      REAL(dp), INTENT(OUT)                              :: values(:), jacobian(:, :)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: physical_jacobian(:, :)

      COMPLEX(dp), ALLOCATABLE                           :: matrix(:, :), operators(:, :, :)
      INTEGER                                            :: c, first, i, ik, last, local, m, n, ns, &
                                                            spin
      INTEGER, ALLOCATABLE                               :: positive(:), sample_spins(:)
      LOGICAL                                            :: owner
      REAL(dp), ALLOCATABLE                              :: auxiliary(:), density(:), &
                                                            expectations(:, :), scales(:, :)

      m = SIZE(values)
      n = qs_skala_scf_size(state)
      CPASSERT(ALL(SHAPE(jacobian) == [m, n]))
      IF (PRESENT(physical_jacobian)) THEN
         CPASSERT(ALL(SHAPE(physical_jacobian) == [m, n]))
         physical_jacobian = 0.0_dp
      END IF
      CALL qs_skala_scf_values(state, domains, points, spins, descriptors, values)
      jacobian = 0.0_dp
      IF (m == 0) RETURN
      IF (state%thermal) THEN
         ALLOCATE (expectations(m, SIZE(state%band_reference)), SOURCE=0.0_dp)
         ALLOCATE (auxiliary(SIZE(state%band_reference)))
      END IF
      ns = SIZE(state%orbitals, 2)
      ALLOCATE (scales(ns, m), SOURCE=0.0_dp)
      DO c = 1, m
         IF (spins(c) == -1 .OR. spins(c) == qs_skala_positive_density_sum) THEN
            scales(:, c) = 1.0_dp
         ELSE IF (ns == 1) THEN
            scales(1, c) = 0.5_dp
         ELSE
            scales(spins(c), c) = 1.0_dp
         END IF
      END DO
      positive = PACK([(i, i=1, m)], spins == qs_skala_positive_density_sum)
      IF (SIZE(positive) > 0) THEN
         ALLOCATE (density(SIZE(positive)), sample_spins(SIZE(positive)))
         DO spin = 1, ns
            sample_spins(:) = MERGE(0, spin, ns == 1)
            CALL qs_skala_scf_values(state, domains(positive), points(:, positive), sample_spins, &
                                     descriptors(positive), density)
            WHERE (density <= 0.0_dp) scales(spin, positive) = 0.0_dp
         END DO
      END IF
      owner = state%para_env%mepos == 0
      IF (state%bloch) owner = state%kpoints%para_env_kp%mepos == 0
      DO c = 1, m
         CALL point_operators(state, domains(c), points(:, c), descriptors(c), operators)
         IF (.NOT. ALLOCATED(matrix)) ALLOCATE (matrix(SIZE(operators, 1), SIZE(operators, 2)))
         DO spin = 1, ns
            IF (scales(spin, c) == 0.0_dp) CYCLE
            DO local = 1, SIZE(state%orbitals, 1)
               IF (.NOT. ASSOCIATED(state%orbitals(local, spin)%native)) CYCLE
               ik = 1
               IF (state%bloch) ik = state%kpoints%kp_env(local)%kpoint_env%nkpoint
               first = state%offsets(ik, spin) + 1
               last = first + state%dimensions(ik, spin) - 1
               IF (.NOT. owner) CYCLE
               IF (state%thermal) CALL thermal_expectations(state, scales(spin, c)*operators(:, :, ik), &
                                                            local, spin, expectations(c, :))
               matrix(:, :) = scales(spin, c)*state%orbitals(local, spin)%weight*operators(:, :, ik)
               CALL skala_orbital_chart_derivative(state%orbitals(local, spin)%chart, matrix, jacobian(c, first:last))
               IF (PRESENT(physical_jacobian)) THEN
                  CALL skala_orbital_chart_derivative(state%orbitals(local, spin)%chart, matrix, &
                                                      physical_jacobian(c, first:last), physical_tangent=.TRUE.)
               END IF
            END DO
         END DO
      END DO
      CALL state%para_env%sum(jacobian)
      IF (PRESENT(physical_jacobian)) CALL state%para_env%sum(physical_jacobian)
      IF (state%thermal) THEN
         CALL state%para_env%sum(expectations)
         DO c = 1, m
            CALL thermal_covector(state, expectations(c, :), .FALSE., auxiliary)
            jacobian(c, state%orbital_size + 1:) = auxiliary
            IF (PRESENT(physical_jacobian)) THEN
               CALL thermal_covector(state, expectations(c, :), .FALSE., auxiliary, &
                                     physical_jacobian(c, state%orbital_size + 1:))
            END IF
         END DO
      END IF
   END SUBROUTINE sampled_point_derivatives

! **************************************************************************************************
!> \brief Reuse fixed-geometry descriptor matrices, retaining fresh orbital and spin derivatives.
!>        Every rank keeps the same keys and memory budget, including empty entry metadata.
!>        Grow the entry array within that budget. A full or oversized cache retains
!>        the original collective construction. Release together with the owning orbital charts.
!> \param state Native orbital state at its defining geometry.
!> \param domain Actual physical evaluation domain.
!> \param point Original FFT or atomic sample key.
!> \param descriptor Density or kinetic-density selector.
!> \param operators Globally reduced matrices in the original AO and k-point ordering.
! **************************************************************************************************
   SUBROUTINE point_operators(state, domain, point, descriptor, operators)
      TYPE(qs_skala_scf_type), INTENT(INOUT)             :: state
      INTEGER, INTENT(IN)                                :: domain, point(3), descriptor
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: operators(:, :, :)

      COMPLEX(dp), ALLOCATABLE                           :: reduced(:, :, :)
      INTEGER                                            :: capacity, ENTRY, free_entry, growth
      INTEGER(KIND=int_8)                                :: bytes, entry_bytes
      TYPE(physical_point_operator)                      :: layout
      TYPE(physical_point_operator), ALLOCATABLE         :: expanded(:)

      free_entry = 0
      IF (ALLOCATED(state%point_operators)) THEN
         DO ENTRY = 1, SIZE(state%point_operators)
            ASSOCIATE (cached => state%point_operators(ENTRY))
               IF (.NOT. ALLOCATED(cached%matrix)) THEN
                  IF (free_entry == 0) free_entry = ENTRY
                  CYCLE
               END IF
               IF (cached%domain /= domain .OR. cached%descriptor /= descriptor) CYCLE
               IF (.NOT. ALL(cached%point == point)) CYCLE
               operators = cached%matrix
               RETURN
            END ASSOCIATE
         END DO
      END IF
      IF (domain == skala_fft_domain .OR. domain == skala_background_domain) THEN
         CALL grid_point_operators(state, domain, point, descriptor, operators)
      ELSE IF (domain == skala_paw_ae_domain .OR. domain == skala_paw_ps_domain) THEN
         CALL paw_point_operators(state, domain, point, descriptor, operators)
      ELSE
         CALL atom_point_operators(state, point, descriptor, operators)
      END IF
      IF (state%bloch) THEN
         CALL kpoint_star_operator(state%inverse%bloch%star, operators, reduced)
         CALL MOVE_ALLOC(reduced, operators)
      END IF
      bytes = 16_int_8*SIZE(operators, KIND=int_8)
      IF (bytes > max_point_operator_bytes - state%point_operator_bytes) RETURN
      IF (free_entry == 0) THEN
         capacity = 0
         IF (ALLOCATED(state%point_operators)) capacity = SIZE(state%point_operators)
         entry_bytes = STORAGE_SIZE(layout, KIND=int_8)/8_int_8
         CPASSERT(entry_bytes > 0_int_8)
         growth = INT(MIN(INT(MAX(initial_point_operator_entries, capacity), int_8), &
                          (max_point_operator_bytes - state%point_operator_bytes - bytes)/entry_bytes))
         IF (growth <= 0) RETURN
         ALLOCATE (expanded(capacity + growth))
         DO ENTRY = 1, capacity
            expanded(ENTRY)%domain = state%point_operators(ENTRY)%domain
            expanded(ENTRY)%point = state%point_operators(ENTRY)%point
            expanded(ENTRY)%descriptor = state%point_operators(ENTRY)%descriptor
            CALL MOVE_ALLOC(state%point_operators(ENTRY)%matrix, expanded(ENTRY)%matrix)
         END DO
         CALL MOVE_ALLOC(expanded, state%point_operators)
         state%point_operator_bytes = state%point_operator_bytes + INT(growth, int_8)*entry_bytes
         free_entry = capacity + 1
      END IF
      ASSOCIATE (cached => state%point_operators(free_entry))
         cached%domain = domain
         cached%point = point
         cached%descriptor = descriptor
         cached%matrix = operators
      END ASSOCIATE
      state%point_operator_bytes = state%point_operator_bytes + bytes
   END SUBROUTINE point_operators

! **************************************************************************************************
!> \brief Construct descriptor operators for one sample on the original FFT grid.
!> \param state Native orbital and field state.
!> \param domain Physical evaluation domain.
!> \param point Global FFT indices or original atomic sample indices.
!> \param descriptor Density or kinetic density selector.
!> \param operators Descriptor matrices in the original Gaussian orbital space.
! **************************************************************************************************
   SUBROUTINE grid_point_operators(state, domain, point, descriptor, operators)
      TYPE(qs_skala_scf_type), INTENT(INOUT)             :: state
      INTEGER, INTENT(IN)                                :: domain, point(3), descriptor
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: operators(:, :, :)

      COMPLEX(dp), ALLOCATABLE                           :: matrix(:, :)
      INTEGER                                            :: axis, bounds(2, 3), i, ik, nao, nk, &
                                                            smooth
      REAL(dp), ALLOCATABLE                              :: fields(:, :, :, :), values(:, :)
      TYPE(cell_type), POINTER                           :: cell
      TYPE(pw_grid_type), POINTER                        :: grid
      TYPE(pw_r3d_rs_type)                               :: density, filtered, gradient(3), kinetic
      TYPE(pw_r3d_rs_type), POINTER                      :: rho_r(:)
      TYPE(qs_upf_inverse_augmentation_type), POINTER    :: augmentation

      CALL qs_rho_get(state%rho, rho_r=rho_r)
      grid => rho_r(1)%pw_grid
      IF (domain == skala_background_domain) grid => state%inverse%xc_grid
      CPASSERT(ASSOCIATED(grid))
      CALL density%create(grid)
      CALL kinetic%create(grid)
      density%array = 0.0_dp
      kinetic%array = 0.0_dp
      DO axis = 1, 3
         CALL gradient(axis)%create(grid)
         gradient(axis)%array = 0.0_dp
      END DO
      IF (domain == skala_fft_domain) THEN
         IF (ALL(point >= grid%bounds_local(1, :)) .AND. ALL(point <= grid%bounds_local(2, :))) THEN
            IF (descriptor == skala_density_descriptor) THEN
               density%array(point(1), point(2), point(3)) = 1.0_dp/grid%dvol
            ELSE
               kinetic%array(point(1), point(2), point(3)) = 1.0_dp/grid%dvol
            END IF
         END IF
         smooth = section_get_ival(state%qs_env%input, "DFT%XC%XC_GRID%XC_SMOOTH_RHO")
         IF (descriptor == skala_density_descriptor .AND. smooth /= xc_rho_no_smooth) THEN
            CALL filtered%create(grid)
            ! The native smoothing stencils are symmetric on the auxiliary PW grid.
            CALL xc_pw_smooth(density, filtered, smooth)
            density%array(:, :, :) = filtered%array
            CALL filtered%release()
         END IF
      ELSE
         ALLOCATE (fields(grid%npts(1), grid%npts(2), grid%npts(3), 1), &
                   values(SIZE(state%inverse%atom_grid%atoms), 1))
         fields = 0.0_dp
         values = 0.0_dp
         DO i = 1, SIZE(values, 1)
            IF (state%inverse%atom_grid%atoms(i) /= point(1)) CYCLE
            IF (state%inverse%atom_grid%atom_rows(i) == point(2)) values(i, 1) = 1.0_dp
         END DO
         CALL get_qs_env(state%qs_env, cell=cell)
         CALL native_grid_batch_adjoint(values, grid, cell, state%inverse%atom_grid%points, fields, &
                                        wrap_auxiliary_cell=.TRUE., continuous=.TRUE.)
         CALL state%para_env%sum(fields)
         bounds = grid%bounds_local - SPREAD(grid%bounds(1, :) - 1, 1, 2)
         IF (descriptor == skala_density_descriptor) THEN
            density%array(:, :, :) = fields(bounds(1, 1):bounds(2, 1), bounds(1, 2):bounds(2, 2), &
                                            bounds(1, 3):bounds(2, 3), 1)/grid%dvol
         ELSE
            kinetic%array(:, :, :) = fields(bounds(1, 1):bounds(2, 1), bounds(1, 2):bounds(2, 2), &
                                            bounds(1, 3):bounds(2, 3), 1)/grid%dvol
         END IF
      END IF
      nk = SIZE(state%dimensions, 1)
      IF (state%bloch) THEN
         nk = SIZE(state%inverse%bloch%weights)
         nao = SIZE(state%inverse%bloch%frames(1)%ao%coefficients, 2)
      ELSE
         nao = SIZE(state%inverse%frame%ao%coefficients, 2)
      END IF
      ALLOCATE (operators(nao, nao, nk))
      IF (state%inverse%gaussian_background) THEN
         IF (descriptor == skala_density_descriptor) THEN
            CALL qs_upf_inverse_background_integrate(state%inverse%coulomb%background, state%qs_env, density, operators)
         ELSE
            CALL qs_upf_inverse_background_integrate(state%inverse%coulomb%background, state%qs_env, kinetic, operators, &
                                                     kinetic=.TRUE.)
         END IF
      ELSE
         augmentation => state%inverse%augmentation
         IF (state%bloch) augmentation => state%inverse%bloch%augmentation
         IF (state%inverse%native_gaussian) THEN
            IF (descriptor == skala_density_descriptor) THEN
               CALL qs_upf_inverse_background_integrate(state%inverse%gaussian, state%qs_env, density, operators)
            ELSE
               CALL qs_upf_inverse_background_integrate(state%inverse%gaussian, state%qs_env, kinetic, operators, &
                                                        kinetic=.TRUE.)
            END IF
         END IF
         DO ik = 1, nk
            IF (state%bloch) THEN
               CALL qs_upf_inverse_grid_integrate_dense(state%inverse%bloch%frames(ik), density, gradient, kinetic, &
                                                        matrix, augmentation=augmentation, fourier_charge=.TRUE., &
                                                        kinetic=state%inverse%kinetic_cache, &
                                                        orbital_fields=.NOT. state%inverse%native_gaussian)
            ELSE
               CALL qs_upf_inverse_grid_integrate_dense(state%inverse%frame, density, gradient, kinetic, matrix, &
                                                        augmentation=augmentation, fourier_charge=.TRUE., &
                                                        kinetic=state%inverse%kinetic_cache, &
                                                        orbital_fields=.NOT. state%inverse%native_gaussian)
            END IF
            IF (state%inverse%native_gaussian) THEN
               operators(:, :, ik) = operators(:, :, ik) + matrix
            ELSE
               operators(:, :, ik) = matrix
            END IF
         END DO
      END IF
      CALL density%release()
      CALL kinetic%release()
      DO axis = 1, 3
         CALL gradient(axis)%release()
      END DO
   END SUBROUTINE grid_point_operators

! **************************************************************************************************
!> \brief Construct descriptor operators for one inverse reconstruction sample.
!> \param state Native orbital and field state.
!> \param point Global FFT indices or original atomic sample indices.
!> \param descriptor Density or kinetic density selector.
!> \param operators Descriptor matrices in the original Gaussian orbital space.
! **************************************************************************************************
   SUBROUTINE atom_point_operators(state, point, descriptor, operators)
      TYPE(qs_skala_scf_type), INTENT(INOUT)             :: state
      INTEGER, INTENT(IN)                                :: point(3), descriptor
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: operators(:, :, :)

      COMPLEX(dp), ALLOCATABLE                           :: correction(:, :), &
                                                            orbital_gradient(:, :, :), q(:, :), &
                                                            soft(:, :)
      INTEGER                                            :: i, ik, nao, nk, np, row
      REAL(dp), ALLOCATABLE                              :: density(:), gradient(:, :), kinetic(:), &
                                                            points(:, :), weights(:)
      TYPE(mp_comm_type)                                 :: comm
      TYPE(qs_upf_inverse_augmentation_type), POINTER    :: augmentation
      TYPE(qs_upf_inverse_type), POINTER                 :: frame

      CPASSERT(point(3) == 1)
      np = COUNT(state%inverse%atom_grid%atoms == point(1) .AND. &
                 state%inverse%atom_grid%atom_rows == point(2))
      nk = SIZE(state%dimensions, 1)
      IF (state%bloch) nk = SIZE(state%inverse%bloch%weights)
      frame => state%inverse%frame%ao
      IF (state%bloch) frame => state%inverse%bloch%frames(1)%ao
      nao = SIZE(frame%coefficients, 2)
      ALLOCATE (density(np), kinetic(np), gradient(np, 3), weights(np), points(3, np), &
                soft(np, nao), orbital_gradient(np, nao, 3), operators(nao, nao, nk), correction(nao, nao))
      density = 0.0_dp
      kinetic = 0.0_dp
      gradient = 0.0_dp
      weights = 1.0_dp
      IF (descriptor == skala_density_descriptor) THEN
         density = 1.0_dp
      ELSE
         kinetic = 1.0_dp
      END IF
      row = 0
      DO i = 1, SIZE(state%inverse%atom_grid%atoms)
         IF (state%inverse%atom_grid%atoms(i) /= point(1)) CYCLE
         IF (state%inverse%atom_grid%atom_rows(i) /= point(2)) CYCLE
         row = row + 1
         points(:, row) = state%inverse%atom_grid%points(:, i)
      END DO
      CALL comm%set_handle(state%para_env%get_handle())
      augmentation => state%inverse%augmentation
      IF (state%bloch) augmentation => state%inverse%bloch%augmentation
      DO ik = 1, nk
         frame => state%inverse%frame%ao
         IF (state%bloch) frame => state%inverse%bloch%frames(ik)%ao
         ! Reconstruct only the selected sample on its owner, including an empty selection elsewhere.
         CALL qs_upf_inverse_evaluate(frame, points, soft, orbital_gradient)
         CALL upf_inverse_semilocal_pair(soft, soft, orbital_gradient, orbital_gradient, &
                                         weights, density, gradient, kinetic, operators(:, :, ik), comm)
         IF (descriptor == skala_density_descriptor) THEN
            ALLOCATE (q(SIZE(frame%coefficients, 1), SIZE(frame%coefficients, 1)))
            CALL qs_upf_inverse_augmentation_integrate(frame, augmentation, density, gradient, q, &
                                                       points, weights, comm)
            operators(:, :, ik) = operators(:, :, ik) + &
                                  MATMUL(CONJG(TRANSPOSE(frame%coefficients)), MATMUL(q, frame%coefficients))
            DEALLOCATE (q)
            IF (ALLOCATED(state%inverse%paw_reconstruction)) THEN
               CALL qs_upf_inverse_paw_operator(frame, state%inverse%paw_reconstruction, augmentation, &
                                                points, weights, density, gradient, correction, comm)
               operators(:, :, ik) = operators(:, :, ik) + correction
            END IF
         ELSE
            CPASSERT(ALLOCATED(state%inverse%kinetic_cache))
            CALL qs_upf_inverse_kinetic_operator(frame, state%inverse%kinetic_cache, points, &
                                                 weights, kinetic, correction, comm)
            operators(:, :, ik) = operators(:, :, ik) + correction
            IF (ALLOCATED(state%inverse%paw_reconstruction)) THEN
               CALL qs_upf_inverse_kinetic_operator(frame, state%inverse%paw_reconstruction, points, &
                                                    weights, kinetic, correction, comm)
               operators(:, :, ik) = operators(:, :, ik) + correction
            END IF
         END IF
      END DO
   END SUBROUTINE atom_point_operators

! **************************************************************************************************
!> \brief Construct descriptor operators for one PAW partial-wave sample.
!> \param state Native orbital and field state.
!> \param domain Physical evaluation domain.
!> \param point Global FFT indices or original atomic sample indices.
!> \param descriptor Density or kinetic density selector.
!> \param operators Descriptor matrices in the original Gaussian orbital space.
! **************************************************************************************************
   SUBROUTINE paw_point_operators(state, domain, point, descriptor, operators)
      TYPE(qs_skala_scf_type), INTENT(IN)                :: state
      INTEGER, INTENT(IN)                                :: domain, point(3), descriptor
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: operators(:, :, :)

      INTEGER                                            :: center, first, i, ik, ikind, last, nao, &
                                                            nk, side
      REAL(dp), ALLOCATABLE                              :: potential(:, :)
      TYPE(qs_upf_inverse_type), POINTER                 :: frame

      nk = SIZE(state%dimensions, 1)
      IF (state%bloch) nk = SIZE(state%inverse%bloch%weights)
      frame => state%inverse%frame%ao
      IF (state%bloch) frame => state%inverse%bloch%frames(1)%ao
      center = 0
      DO i = 1, SIZE(frame%center_atoms)
         IF (frame%center_atoms(i) /= point(1)) CYCLE
         CPASSERT(center == 0)
         center = i
      END DO
      CPASSERT(center > 0)
      ikind = frame%center_kinds(center)
      CPASSERT(state%inverse%paw(ikind)%active)
      first = frame%projector_first(center)
      last = frame%projector_last(center)
      nao = SIZE(frame%coefficients, 2)
      ALLOCATE (potential(last - first + 1, last - first + 1), operators(nao, nao, nk))
      side = 1
      IF (domain == skala_paw_ps_domain) side = 2
      CALL paw_xc_point_adjoint(state%inverse%paw(ikind)%xcgrid, state%inverse%paw(ikind)%onecenter, &
                                point(2), point(3), side, descriptor == skala_kinetic_descriptor, potential)
      DO ik = 1, nk
         frame => state%inverse%frame%ao
         IF (state%bloch) frame => state%inverse%bloch%frames(ik)%ao
         operators(:, :, ik) = MATMUL(CONJG(TRANSPOSE(frame%coefficients(first:last, :))), &
                                      MATMUL(CMPLX(potential, KIND=dp), frame%coefficients(first:last, :)))
      END DO
   END SUBROUTINE paw_point_operators

! **************************************************************************************************
!> \brief Find a physical root along a normalized orbital curve from the original native fields.
!>        Density matrices and physical fields are never interpolated. On failure restore the
!>        previously evaluated coordinate. On success all native caches describe the returned root.
!> \param state ...
!> \param first Initial global orbital coordinate.
!> \param last Final global orbital coordinate.
!> \param domain ...
!> \param point ...
!> \param spin ...
!> \param descriptor ...
!> \param tolerance Absolute physical root tolerance in atomic units.
!> \param coordinate Root coordinate on success, previous evaluated coordinate on failure.
!> \param value Actual physical input at the best point found by the root search.
!> \param converged Whether the actual physical input meets the requested tolerance.
!> \param TARGET Physical density or kinetic-density threshold, zero by default.
! **************************************************************************************************
   SUBROUTINE qs_skala_scf_root(state, first, last, domain, point, spin, descriptor, tolerance, &
                                coordinate, value, converged, TARGET)
      REAL(dp), INTENT(IN), OPTIONAL                     :: target
      LOGICAL, INTENT(OUT)                               :: converged
      REAL(dp), INTENT(OUT)                              :: value, coordinate(:)
      REAL(dp), INTENT(IN)                               :: tolerance
      INTEGER, INTENT(IN)                                :: descriptor, spin, point(3), domain
      REAL(dp), INTENT(IN)                               :: last(:), first(:)
      TYPE(qs_skala_scf_type), INTENT(INOUT), TARGET     :: state

      REAL(dp)                                           :: energy, step
      REAL(dp), ALLOCATABLE                              :: gradient(:)
      TYPE(physical_root_context)                        :: context

      CPASSERT(state%evaluated)
      CPASSERT(SIZE(first) == qs_skala_scf_size(state) .AND. SIZE(last) == SIZE(first))
      CPASSERT(SIZE(coordinate) == SIZE(first))
      ALLOCATE (context%first(SIZE(first)), context%direction(SIZE(first)))
      context%state => state
      context%first(:) = first
      context%direction(:) = last - first
      context%domain = domain
      context%point = point
      context%spin = spin
      context%descriptor = descriptor
      coordinate(:) = state%coordinate
      CALL skala_density_root(evaluate_physical_input, context, 0.0_dp, 1.0_dp, tolerance, step, value, converged, &
                              TARGET=TARGET)
      IF (converged) coordinate(:) = first + step*context%direction
      ALLOCATE (gradient(SIZE(coordinate)))
      CALL qs_skala_scf_evaluate(state, coordinate, energy, gradient)
      IF (converged) CALL qs_skala_scf_value(state, domain, point, spin, descriptor, value)
   END SUBROUTINE qs_skala_scf_root

! **************************************************************************************************
!> \brief Collect bracketed physical rho/tau targets from two actual native orbital states.
!>        Samples remain on their native owners; only candidate keys and endpoint values are
!>        gathered. Fields equal to the target at both ends supply no candidate. This endpoint test does not
!>        locate roots or exclude an even number of crossings along the nonlinear orbital curve.
!>        Restore the first state, including its complete native KS operator, before returning.
!> \param state ...
!> \param first Initial global orbital coordinate.
!> \param last Final global orbital coordinate.
!> \param points Global FFT indices or (global atom, row within that atom, 1).
!> \param domains FFT, exact atomic or interpolated atomic background for each candidate.
!> \param spins Minus one for the spin sum, zero for paired RKS, one or two for UKS.
!> \param descriptors Density or kinetic density for each candidate.
!> \param endpoint_values Actual first and last physical inputs, one candidate per column.
!> \param targets Density and kinetic-density thresholds, both zero by default.
!> \param total_targets Optional thresholds for the sum of both physical spins.
!> \param positive_density_target Optional cutoff for the sum of the positive spin densities.
!> \param density_only Scan density inputs only, excluding kinetic-density samples, false by default.
!> \param spin_density_targets Distinct spin-density thresholds replacing the density entry of targets.
!> \param candidate_targets Actual threshold of each returned physical candidate.
!> \param restore_operator Restore the complete first energy and KS operator, default true.
!>        Curvature scans can restore only its fields before the next complete evaluation.
! **************************************************************************************************
   SUBROUTINE qs_skala_scf_crossings(state, first, last, points, domains, spins, descriptors, &
                                     endpoint_values, targets, total_targets, positive_density_target, density_only, &
                                     spin_density_targets, candidate_targets, restore_operator)
      TYPE(qs_skala_scf_type), INTENT(INOUT)             :: state
      REAL(dp), INTENT(IN)                               :: first(:), last(:)
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: points(:, :), domains(:), spins(:), &
                                                            descriptors(:)
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: endpoint_values(:, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: targets(2), total_targets(2), &
                                                            positive_density_target
      LOGICAL, INTENT(IN), OPTIONAL                      :: density_only
      REAL(dp), INTENT(IN), OPTIONAL                     :: spin_density_targets(:)
      REAL(dp), ALLOCATABLE, INTENT(OUT), OPTIONAL       :: candidate_targets(:)
      LOGICAL, INTENT(IN), OPTIONAL                      :: restore_operator

      INTEGER                                            :: descriptor, domain, i, level_index, &
                                                            local_count, nlevels, offset, row, &
                                                            selection, spin, total
      INTEGER, ALLOCATABLE                               :: counts(:), keys(:, :)
      LOGICAL                                            :: only_density, valid
      LOGICAL, ALLOCATABLE                               :: mask(:, :)
      REAL(dp)                                           :: energy, level, levels(2), &
                                                            positive_level, total_levels(2)
      REAL(dp), ALLOCATABLE                              :: density_levels(:), first_field(:, :), &
                                                            gradient(:), last_field(:, :), &
                                                            local_targets(:), values(:, :)
      TYPE(physical_domain_samples) :: final(skala_fft_domain:skala_last_sample_domain), &
         initial(skala_fft_domain:skala_last_sample_domain)

      valid = SIZE(first) == qs_skala_scf_size(state) .AND. SIZE(last) == SIZE(first)
      CPASSERT(valid)
      only_density = .FALSE.
      IF (PRESENT(density_only)) only_density = density_only
      levels = 0.0_dp
      IF (PRESENT(targets)) levels = targets
      CPASSERT(ALL(ieee_is_finite(levels)))
      IF (PRESENT(spin_density_targets)) THEN
         CPASSERT(ALL(ieee_is_finite(spin_density_targets)))
         ALLOCATE (density_levels(SIZE(spin_density_targets)))
         density_levels(:) = spin_density_targets
         DO i = 2, SIZE(density_levels)
            CPASSERT(ALL(density_levels(:i - 1) /= density_levels(i)))
         END DO
      ELSE
         ALLOCATE (density_levels(1), SOURCE=levels(skala_density_descriptor))
      END IF
      total_levels = 0.0_dp
      IF (PRESENT(total_targets)) total_levels = total_targets
      CPASSERT(ALL(ieee_is_finite(total_levels)))
      positive_level = 0.0_dp
      IF (PRESENT(positive_density_target)) positive_level = positive_density_target
      CPASSERT(ieee_is_finite(positive_level) .AND. positive_level >= 0.0_dp)
      CALL sample_native_constraints(state, first)
      DO domain = skala_fft_domain, skala_last_sample_domain
         CALL skala_density_samples_get(state%constraints, domain, initial(domain)%points, &
                                        initial(domain)%density, initial(domain)%kinetic, initial(domain)%available)
      END DO
      CALL sample_native_constraints(state, last)
      local_count = 0
      DO domain = skala_fft_domain, skala_last_sample_domain
         CALL skala_density_samples_get(state%constraints, domain, final(domain)%points, &
                                        final(domain)%density, final(domain)%kinetic, final(domain)%available)
         valid = initial(domain)%available .EQV. final(domain)%available
         CPASSERT(valid)
         IF (.NOT. initial(domain)%available) CYCLE
         valid = ALL(SHAPE(initial(domain)%points) == SHAPE(final(domain)%points))
         CPASSERT(valid)
         valid = ALL(initial(domain)%points == final(domain)%points)
         CPASSERT(valid)
         valid = ALL(SHAPE(initial(domain)%density) == SHAPE(final(domain)%density))
         CPASSERT(valid)
         valid = ALL(SHAPE(initial(domain)%kinetic) == SHAPE(final(domain)%kinetic))
         CPASSERT(valid)
         DO descriptor = skala_density_descriptor, skala_kinetic_descriptor
            IF (only_density .AND. descriptor /= skala_density_descriptor) CYCLE
            IF (descriptor == skala_kinetic_descriptor .AND. SIZE(initial(domain)%kinetic, 2) == 0) CYCLE
            DO selection = 0, 2
               IF (selection == 1 .AND. .NOT. PRESENT(total_targets)) CYCLE
               IF (selection == 2) THEN
                  IF (descriptor /= skala_density_descriptor .OR. .NOT. PRESENT(positive_density_target)) CYCLE
               END IF
               CALL sample_descriptor(initial(domain), descriptor, selection, first_field)
               CALL sample_descriptor(final(domain), descriptor, selection, last_field)
               ALLOCATE (mask(SIZE(first_field, 1), SIZE(first_field, 2)))
               nlevels = 1
               IF (selection == 0 .AND. descriptor == skala_density_descriptor) nlevels = SIZE(density_levels)
               DO level_index = 1, nlevels
                  level = levels(descriptor)
                  IF (selection == 0 .AND. descriptor == skala_density_descriptor) level = density_levels(level_index)
                  IF (selection == 1) level = total_levels(descriptor)
                  IF (selection == 2) level = positive_level
                  CALL skala_density_crossings(first_field, last_field, mask, level)
                  local_count = local_count + COUNT(mask)
               END DO
               DEALLOCATE (mask, first_field, last_field)
            END DO
         END DO
      END DO
      ALLOCATE (keys(6, local_count), values(2, local_count), local_targets(local_count))
      row = 0
      DO domain = skala_fft_domain, skala_last_sample_domain
         IF (.NOT. initial(domain)%available) CYCLE
         DO descriptor = skala_density_descriptor, skala_kinetic_descriptor
            IF (only_density .AND. descriptor /= skala_density_descriptor) CYCLE
            IF (descriptor == skala_kinetic_descriptor .AND. SIZE(initial(domain)%kinetic, 2) == 0) CYCLE
            DO selection = 0, 2
               IF (selection == 1 .AND. .NOT. PRESENT(total_targets)) CYCLE
               IF (selection == 2) THEN
                  IF (descriptor /= skala_density_descriptor .OR. .NOT. PRESENT(positive_density_target)) CYCLE
               END IF
               CALL sample_descriptor(initial(domain), descriptor, selection, first_field)
               CALL sample_descriptor(final(domain), descriptor, selection, last_field)
               ALLOCATE (mask(SIZE(first_field, 1), SIZE(first_field, 2)))
               nlevels = 1
               IF (selection == 0 .AND. descriptor == skala_density_descriptor) nlevels = SIZE(density_levels)
               DO level_index = 1, nlevels
                  level = levels(descriptor)
                  IF (selection == 0 .AND. descriptor == skala_density_descriptor) level = density_levels(level_index)
                  IF (selection == 1) level = total_levels(descriptor)
                  IF (selection == 2) level = positive_level
                  CALL skala_density_crossings(first_field, last_field, mask, level)
                  DO spin = 1, SIZE(mask, 2)
                     DO i = 1, SIZE(mask, 1)
                        IF (.NOT. mask(i, spin)) CYCLE
                        row = row + 1
                        keys(1, row) = domain
                        keys(2:4, row) = initial(domain)%points(:, i)
                        keys(5, row) = spin
                        IF (SIZE(mask, 2) == 1) keys(5, row) = 0
                        IF (selection == 1) keys(5, row) = -1
                        IF (selection == 2) keys(5, row) = qs_skala_positive_density_sum
                        keys(6, row) = descriptor
                        values(:, row) = [first_field(i, spin), last_field(i, spin)]
                        local_targets(row) = level
                     END DO
                  END DO
               END DO
               DEALLOCATE (mask, first_field, last_field)
            END DO
         END DO
      END DO
      CPASSERT(row == local_count)
      ALLOCATE (counts(state%para_env%num_pe))
      counts = 0
      counts(state%para_env%mepos + 1) = local_count
      CALL state%para_env%sum(counts)
      total = SUM(counts)
      offset = SUM(counts(:state%para_env%mepos))
      ALLOCATE (points(3, total), domains(total), spins(total), descriptors(total), endpoint_values(2, total))
      points = 0
      domains = 0
      spins = 0
      descriptors = 0
      endpoint_values = 0.0_dp
      points(:, offset + 1:offset + local_count) = keys(2:4, :)
      domains(offset + 1:offset + local_count) = keys(1, :)
      spins(offset + 1:offset + local_count) = keys(5, :)
      descriptors(offset + 1:offset + local_count) = keys(6, :)
      endpoint_values(:, offset + 1:offset + local_count) = values
      IF (PRESENT(candidate_targets)) THEN
         ALLOCATE (candidate_targets(total), SOURCE=0.0_dp)
         candidate_targets(offset + 1:offset + local_count) = local_targets
      END IF
      IF (total > 0) THEN
         CALL state%para_env%sum(points)
         CALL state%para_env%sum(domains)
         CALL state%para_env%sum(spins)
         CALL state%para_env%sum(descriptors)
         CALL state%para_env%sum(endpoint_values)
         IF (PRESENT(candidate_targets)) CALL state%para_env%sum(candidate_targets)
      END IF
      IF (PRESENT(restore_operator)) THEN
         IF (.NOT. restore_operator) THEN
            CALL sample_native_constraints(state, first)
            RETURN
         END IF
      END IF
      ALLOCATE (gradient(SIZE(first)))
      CALL qs_skala_scf_evaluate(state, first, energy, gradient)
   END SUBROUTINE qs_skala_scf_crossings

! **************************************************************************************************
!> \brief Locate native input roots for every crossing detected between two orbital states.
!>        Retain every physical sample, including coincident roots. Endpoint signs do not exclude
!>        additional crossings within the nonlinear orbital segment. Selection, admissibility and
!>        stationarity remain separate operations. Restore the initial native state before returning.
!> \param state Native orbital and physical input state.
!> \param first Initial normalized orbital coordinates.
!> \param last Final normalized orbital coordinates.
!> \param tolerances Absolute density and kinetic-density root tolerances in atomic units.
!> \param points Global FFT indices or native atomic sample indices, one candidate per column.
!> \param domains Physical domain of each candidate.
!> \param spins Physical spin selector of each candidate.
!> \param descriptors Density or kinetic density for each candidate.
!> \param coordinates Located root coordinates, or first for an unresolved candidate.
!> \param values Actual input at the best evaluated root trial.
!> \param converged Whether each candidate satisfies its original physical root tolerance.
!> \param targets Optional density and kinetic thresholds, zero by default.
!> \param total_targets Optional thresholds for the physical spin sums.
!> \param positive_density_target Optional threshold for the sum of positive spin densities.
!> \param density_only Scan density inputs only, false by default.
! **************************************************************************************************
   SUBROUTINE qs_skala_scf_crossing_roots(state, first, last, tolerances, points, domains, spins, &
                                          descriptors, coordinates, values, converged, targets, &
                                          total_targets, positive_density_target, density_only)
      TYPE(qs_skala_scf_type), INTENT(INOUT)             :: state
      REAL(dp), INTENT(IN)                               :: first(:), last(:), tolerances(2)
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: points(:, :), domains(:), spins(:), &
                                                            descriptors(:)
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: coordinates(:, :), values(:)
      LOGICAL, ALLOCATABLE, INTENT(OUT)                  :: converged(:)
      REAL(dp), INTENT(IN), OPTIONAL                     :: targets(2), total_targets(2), &
                                                            positive_density_target
      LOGICAL, INTENT(IN), OPTIONAL                      :: density_only

      INTEGER                                            :: i, m
      REAL(dp)                                           :: energy, level, levels(2), &
                                                            positive_level, total_levels(2)
      REAL(dp), ALLOCATABLE                              :: endpoints(:, :), gradient(:)

      CPASSERT(ALL(ieee_is_finite(tolerances)) .AND. ALL(tolerances > 0.0_dp))
      levels = 0.0_dp
      total_levels = 0.0_dp
      positive_level = 0.0_dp
      IF (PRESENT(targets)) levels = targets
      IF (PRESENT(total_targets)) total_levels = total_targets
      IF (PRESENT(positive_density_target)) positive_level = positive_density_target
      CALL qs_skala_scf_crossings(state, first, last, points, domains, spins, descriptors, endpoints, &
                                  targets, total_targets, positive_density_target, density_only)
      m = SIZE(spins)
      ALLOCATE (coordinates(SIZE(first), m), values(m), converged(m), gradient(SIZE(first)))
      DO i = 1, m
         level = levels(descriptors(i))
         IF (spins(i) == -1) level = total_levels(descriptors(i))
         IF (spins(i) == qs_skala_positive_density_sum) level = positive_level
         CALL qs_skala_scf_root(state, first, last, domains(i), points(:, i), spins(i), descriptors(i), &
                                tolerances(descriptors(i)), coordinates(:, i), values(i), converged(i), TARGET=level)
         IF (.NOT. converged(i)) coordinates(:, i) = first
      END DO
      CALL qs_skala_scf_evaluate(state, first, energy, gradient)
   END SUBROUTINE qs_skala_scf_crossing_roots

! **************************************************************************************************
!> \brief Choose the earliest supplied converged root on a native orbital trial segment.
!>        At that root retain every candidate whose actual physical value meets its tolerance,
!>        including coincident or dependent constraints. Incomplete root searches return failure
!>        and restore the first state. Unrecorded interior roots, normal admissibility and a
!>        minimum are not established by this selection.
!> \param state Native orbital and physical input state.
!> \param first Initial normalized orbital coordinates.
!> \param last Final normalized orbital coordinates.
!> \param coordinates Supplied root coordinates, one column per candidate.
!> \param points Physical sample keys in the same candidate order.
!> \param domains Physical evaluation domains.
!> \param spins Physical spin selectors.
!> \param descriptors Density or kinetic density for each candidate.
!> \param converged Whether each supplied root was located successfully.
!> \param tolerances Absolute density and kinetic-density tolerances in atomic units.
!> \param root Selected coordinates, or first on failure.
!> \param selected Candidates meeting their original tolerance at the selected root.
!> \param found Whether a complete supplied search yields a verified physical target.
!> \param targets Fixed physical target of each candidate, zero by default.
! **************************************************************************************************
   SUBROUTINE qs_skala_scf_choose_root(state, first, last, coordinates, points, domains, spins, &
                                       descriptors, converged, tolerances, root, selected, found, targets)
      TYPE(qs_skala_scf_type), INTENT(INOUT)             :: state
      REAL(dp), INTENT(IN)                               :: first(:), last(:), coordinates(:, :)
      INTEGER, INTENT(IN)                                :: points(:, :), domains(:), spins(:), &
                                                            descriptors(:)
      LOGICAL, INTENT(IN)                                :: converged(:)
      REAL(dp), INTENT(IN)                               :: tolerances(2)
      REAL(dp), INTENT(OUT)                              :: root(:)
      LOGICAL, ALLOCATABLE, INTENT(OUT)                  :: selected(:)
      LOGICAL, INTENT(OUT)                               :: found
      REAL(dp), INTENT(IN), OPTIONAL                     :: targets(:)

      INTEGER                                            :: earliest, i, m, n
      REAL(dp)                                           :: energy, length2, level, roundoff, value
      REAL(dp), ALLOCATABLE                              :: direction(:), gradient(:), steps(:)

      n = qs_skala_scf_size(state)
      m = SIZE(spins)
      CPASSERT(SIZE(first) == n .AND. SIZE(last) == n .AND. SIZE(root) == n)
      CPASSERT(ALL(SHAPE(coordinates) == [n, m]) .AND. ALL(SHAPE(points) == [3, m]))
      CPASSERT(SIZE(domains) == m .AND. SIZE(descriptors) == m .AND. SIZE(converged) == m)
      CPASSERT(ALL(descriptors >= 1 .AND. descriptors <= 2))
      CPASSERT(ALL(ieee_is_finite(coordinates)))
      CPASSERT(ALL(ieee_is_finite(tolerances)) .AND. ALL(tolerances > 0.0_dp))
      IF (PRESENT(targets)) THEN
         CPASSERT(SIZE(targets) == m)
         CPASSERT(ALL(ieee_is_finite(targets)))
      END IF
      ALLOCATE (selected(m), direction(n), gradient(n), steps(m))
      selected = .FALSE.
      found = .FALSE.
      root(:) = first
      CALL qs_skala_scf_evaluate(state, first, energy, gradient)
      IF (m == 0 .OR. .NOT. ALL(converged)) RETURN
      direction(:) = last - first
      length2 = DOT_PRODUCT(direction, direction)
      IF (length2 <= TINY(0.0_dp)) RETURN
      roundoff = 64.0_dp*EPSILON(0.0_dp)*MAX(NORM2(first), NORM2(last), TINY(0.0_dp))
      DO i = 1, m
         steps(i) = DOT_PRODUCT(coordinates(:, i) - first, direction)/length2
         IF (steps(i) < -64.0_dp*EPSILON(0.0_dp) .OR. steps(i) > 1.0_dp + 64.0_dp*EPSILON(0.0_dp)) RETURN
         IF (NORM2(coordinates(:, i) - first - steps(i)*direction) > roundoff) RETURN
      END DO
      earliest = MINLOC(steps, DIM=1)
      root(:) = coordinates(:, earliest)
      CALL qs_skala_scf_evaluate(state, root, energy, gradient)
      DO i = 1, m
         CALL qs_skala_scf_value(state, domains(i), points(:, i), spins(i), descriptors(i), value)
         level = 0.0_dp
         IF (PRESENT(targets)) level = targets(i)
         selected(i) = ABS(value - level) <= tolerances(descriptors(i))
      END DO
      found = selected(earliest)
      IF (.NOT. found) THEN
         selected = .FALSE.
         root(:) = first
         CALL qs_skala_scf_evaluate(state, first, energy, gradient)
      END IF
   END SUBROUTINE qs_skala_scf_choose_root

! **************************************************************************************************
!> \brief Resolve a shared physical endpoint bracket with all supplied candidate fields.
!>        Each trial evaluates all candidates together. No unselected candidate can lie
!>        strictly on its opposite side at the returned boundary. Endpoint brackets do not
!>        exclude additional interior crossings. Energy acceptance and minimum conditions
!>        remain separate from this physical root search. Failure restores the first state.
!> \param state Native normalized-orbital state.
!> \param first Initial orbital coordinates.
!> \param last Final orbital coordinates.
!> \param points Original global sample keys, one column per candidate.
!> \param domains Physical evaluation domains.
!> \param spins Physical spin selectors.
!> \param descriptors Density or kinetic density per candidate.
!> \param targets Fixed physical target per candidate.
!> \param tolerances Absolute physical feasibility window per candidate.
!> \param root Returned coordinates, or first on failure.
!> \param values Actual physical inputs at the returned state.
!> \param selected All candidates satisfying their physical window at the root.
!> \param found Whether the returned boundary satisfies the original physical criteria.
!> \param max_iterations Maximum shared bracket refinements, defaults to 100.
!> \param evaluations Optional number of bracket field requests and complete root checks, excluding curve restoration.
!> \param check_external Poll the normal CP2K external stop conditions, false by default.
!> \param stopped Whether an external stop was requested.
!> \param curve Optional original minimization chart, restoring active physical constraints at every root trial.
! **************************************************************************************************
   SUBROUTINE qs_skala_scf_boundary_root(state, first, last, points, domains, spins, descriptors, &
                                         targets, tolerances, root, values, selected, found, &
                                         max_iterations, evaluations, check_external, stopped, curve)
      TYPE(qs_skala_scf_type), INTENT(INOUT), TARGET     :: state
      REAL(dp), INTENT(IN)                               :: first(:), last(:)
      INTEGER, INTENT(IN)                                :: points(:, :), domains(:), spins(:), &
                                                            descriptors(:)
      REAL(dp), INTENT(IN)                               :: targets(:), tolerances(:)
      REAL(dp), INTENT(OUT)                              :: root(:), values(:)
      LOGICAL, INTENT(OUT)                               :: selected(:), found
      INTEGER, INTENT(IN), OPTIONAL                      :: max_iterations
      INTEGER, INTENT(OUT), OPTIONAL                     :: evaluations
      LOGICAL, INTENT(IN), OPTIONAL                      :: check_external
      LOGICAL, INTENT(OUT), OPTIONAL                     :: stopped
      TYPE(skala_constraint_chart_type), INTENT(IN), &
         OPTIONAL                                        :: curve

      INTEGER                                            :: count, i, iteration, limit, m, n
      INTEGER, ALLOCATABLE                               :: rows(:)
      LOGICAL                                            :: bracketed, EXTERNAL, requested, valid
      REAL(dp)                                           :: energy, fraction, left, previous_width, &
                                                            right, scale, trial_fraction
      REAL(dp), ALLOCATABLE                              :: direction(:), edge_values(:, :), &
                                                            gradient(:), initial(:), &
                                                            jacobian(:, :), residuals(:), sides(:)
      TYPE(physical_constraints_context)                 :: context

      n = qs_skala_scf_size(state)
      m = SIZE(targets)
      CPASSERT(SIZE(first) == n .AND. SIZE(last) == n .AND. SIZE(root) == n)
      CPASSERT(SIZE(values) == m .AND. SIZE(selected) == m .AND. SIZE(tolerances) == m)
      CPASSERT(ALL(ieee_is_finite(targets)) .AND. ALL(ieee_is_finite(tolerances)))
      CPASSERT(ALL(tolerances > 0.0_dp))
      limit = 100
      IF (PRESENT(max_iterations)) limit = max_iterations
      CPASSERT(limit > 0)
      EXTERNAL = .FALSE.
      IF (PRESENT(check_external)) EXTERNAL = check_external
      requested = .FALSE.
      IF (PRESENT(stopped)) stopped = .FALSE.
      ALLOCATE (initial(m), sides(m), edge_values(2, m), direction(n), gradient(n))
      ALLOCATE (rows(m), SOURCE=-1)
      IF (PRESENT(curve)) THEN
         context%state => state
         CALL qs_skala_scf_selected(state, context%points, context%domains, context%spins, context%descriptors, &
                                    context%targets, residuals, valid, require_feasible=.FALSE.)
         CPASSERT(valid)
         ALLOCATE (jacobian(SIZE(residuals), n))
      END IF
      found = .FALSE.
      selected = .FALSE.
      root(:) = first
      count = 1
      IF (m == 0) THEN
         CALL qs_skala_scf_evaluate(state, first, energy, gradient)
         IF (PRESENT(evaluations)) evaluations = count
         RETURN
      END IF
      ! Initialize the bracket from physical fields alone. Every accepted root,
      ! including an initial endpoint, still receives a complete functional check.
      CALL sample_native_constraints(state, first)
      CALL qs_skala_scf_values(state, domains, points, spins, descriptors, initial, rows)
      values(:) = initial
      selected(:) = ABS(values - targets) <= tolerances
      IF (ANY(selected)) THEN
         ! Keep the fields that passed the physical window while rebuilding full KS.
         CALL evaluate_limiting_state(state, first, energy, gradient, .TRUE.)
         count = count + 1
         CALL qs_skala_scf_values(state, domains, points, spins, descriptors, initial, rows)
         values(:) = initial
         selected(:) = ABS(values - targets) <= tolerances
         found = ANY(selected)
      END IF
      IF (.NOT. found) THEN
         sides(:) = SIGN(1.0_dp, initial - targets)
         direction(:) = last - first
         CALL sample_native_constraints(state, last)
         count = count + 1
         CALL qs_skala_scf_values(state, domains, points, spins, descriptors, values, rows)
         IF (ALL(sides*(values - targets) <= tolerances)) THEN
            left = 0.0_dp
            right = 1.0_dp
            previous_width = 2.0_dp*(right - left)
            edge_values(1, :) = sides*(initial - targets)
            edge_values(2, :) = sides*(values - targets)
            DO iteration = 1, limit
               IF (EXTERNAL) THEN
                  CALL external_control(requested, "SCF", target_time=state%qs_env%target_time, &
                                        start_time=state%qs_env%start_time)
                  IF (requested) EXIT
               END IF
               ! Predict the earliest bracketed crossing from actual endpoint fields.
               ! Safeguarding retains interval contraction. Only native fields accept a root.
               trial_fraction = 1.0_dp
               bracketed = .FALSE.
               DO i = 1, m
                  IF (edge_values(2, i) > 0.0_dp) CYCLE
                  bracketed = .TRUE.
                  scale = MAX(ABS(edge_values(1, i)), ABS(edge_values(2, i)))
                  trial_fraction = MIN(trial_fraction, (edge_values(1, i)/scale)/ &
                                       ((edge_values(1, i)/scale) - (edge_values(2, i)/scale)))
               END DO
               ! Keep an accurate interpolation near either endpoint. If the last
               ! trial did not halve the bracket, bisect this interval instead.
               IF (.NOT. bracketed .OR. right - left > 0.5_dp*previous_width) trial_fraction = 0.5_dp
               previous_width = right - left
               fraction = (1.0_dp - trial_fraction)*left + trial_fraction*right
               IF (fraction <= left .OR. fraction >= right) fraction = left + 0.5_dp*(right - left)
               IF (fraction <= left .OR. fraction >= right) EXIT
               root(:) = first + fraction*direction
               IF (PRESENT(curve)) THEN
                  CALL skala_constraint_chart_restore(curve, evaluate_native_constraints, context, root, &
                                                      residuals, jacobian, valid)
                  IF (.NOT. valid) EXIT
               END IF
               CALL sample_native_constraints(state, root)
               count = count + 1
               CALL qs_skala_scf_values(state, domains, points, spins, descriptors, values, rows)
               selected(:) = ABS(values - targets) <= tolerances
               IF (ANY(selected) .AND. ALL(sides*(values - targets) >= -tolerances)) THEN
                  CALL evaluate_limiting_state(state, root, energy, gradient, .TRUE.)
                  count = count + 1
                  CALL qs_skala_scf_values(state, domains, points, spins, descriptors, values, rows)
                  selected(:) = ABS(values - targets) <= tolerances
                  found = ANY(selected) .AND. ALL(sides*(values - targets) >= -tolerances)
                  IF (found) EXIT
               END IF
               IF (ANY(sides*(values - targets) < -tolerances)) THEN
                  right = fraction
                  edge_values(2, :) = sides*(values - targets)
               ELSE
                  left = fraction
                  edge_values(1, :) = sides*(values - targets)
               END IF
            END DO
         END IF
      END IF
      IF (.NOT. found) THEN
         root(:) = first
         selected = .FALSE.
         CALL qs_skala_scf_evaluate(state, first, energy, gradient)
         count = count + 1
         CALL qs_skala_scf_values(state, domains, points, spins, descriptors, values, rows)
      END IF
      IF (PRESENT(evaluations)) evaluations = count
      IF (PRESENT(stopped)) stopped = requested
   END SUBROUTINE qs_skala_scf_boundary_root

! **************************************************************************************************
!> \brief Select the requested density or kinetic-density values from recorded samples.
!> \param sample Recorded physical fields.
!> \param descriptor Density or kinetic density selector.
!> \param selection Selected sample rows.
!> \param field Selected descriptor values.
! **************************************************************************************************
   SUBROUTINE sample_descriptor(sample, descriptor, selection, field)
      TYPE(physical_domain_samples), INTENT(IN)          :: sample
      INTEGER, INTENT(IN)                                :: descriptor, selection
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: field(:, :)

      INTEGER                                            :: np, ns

      np = SIZE(sample%density, 1)
      ns = SIZE(sample%density, 2)
      IF (selection > 0) THEN
         ALLOCATE (field(np, 1))
         IF (descriptor == skala_density_descriptor) THEN
            IF (selection == 2) THEN
               field(:, 1) = SUM(MAX(sample%density, 0.0_dp), DIM=2)*2.0_dp/REAL(ns, dp)
            ELSE
               field(:, 1) = SUM(sample%density, DIM=2)*2.0_dp/REAL(ns, dp)
            END IF
         ELSE
            field(:, 1) = SUM(sample%kinetic, DIM=2)*2.0_dp/REAL(ns, dp)
         END IF
      ELSE
         ALLOCATE (field(np, ns))
         IF (descriptor == skala_density_descriptor) THEN
            field(:, :) = sample%density
         ELSE
            field(:, :) = sample%kinetic
         END IF
      END IF
   END SUBROUTINE sample_descriptor

! **************************************************************************************************
!> \brief Evaluate several actual physical constraints and their joint native orbital Jacobian.
!>        Reuse a complete physical-field collection only at the same coordinates and native
!>        density stamp. The energy and KS require reevaluation afterwards, including when
!>        raw fields are reused, because selected derivative weights can change independently.
!>        Density and kinetic constraints retain distinct physical units. Point Jacobians are
!>        evaluated at the current state, without interpolation or reuse from another trial.
!> \param state ...
!> \param coordinate Global normalized-orbital coordinates.
!> \param points Global FFT indices or (global atom, atomic row, 1), one column per constraint.
!> \param domains Physical evaluation domain of each constraint.
!> \param spins Zero for paired RKS, one or two for UKS.
!> \param descriptors Density or kinetic density of each constraint.
!> \param values Current physical inputs minus their targets, in the original physical units.
!> \param jacobian One constraint per row, one global orbital coordinate per column.
!> \param physical_jacobian Optional covectors in the orthonormal physical orbital tangent.
!> \param targets Optional fixed physical target per constraint, zero by default.
! **************************************************************************************************
   SUBROUTINE qs_skala_scf_constraints(state, coordinate, points, domains, spins, descriptors, &
                                       values, jacobian, physical_jacobian, targets)
      TYPE(qs_skala_scf_type), INTENT(INOUT)             :: state
      REAL(dp), INTENT(IN)                               :: coordinate(:)
      INTEGER, INTENT(IN)                                :: points(:, :), domains(:), spins(:), &
                                                            descriptors(:)
      REAL(dp), INTENT(OUT)                              :: values(:), jacobian(:, :)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: physical_jacobian(:, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: targets(:)

      INTEGER                                            :: m, n
      LOGICAL                                            :: valid

      m = SIZE(values)
      n = qs_skala_scf_size(state)
      valid = SIZE(coordinate) == n .AND. ALL(SHAPE(points) == [3, m])
      CPASSERT(valid)
      valid = SIZE(domains) == m .AND. SIZE(spins) == m .AND. SIZE(descriptors) == m
      CPASSERT(valid)
      CPASSERT(ALL(SHAPE(jacobian) == [m, n]))
      IF (PRESENT(targets)) THEN
         CPASSERT(SIZE(targets) == m)
         CPASSERT(ALL(ieee_is_finite(targets)))
      END IF
      IF (PRESENT(physical_jacobian)) THEN
         valid = ALL(SHAPE(physical_jacobian) == [m, n])
         CPASSERT(valid)
      END IF
      CALL sample_native_constraints(state, coordinate)
      CALL sampled_point_derivatives(state, domains, points, spins, descriptors, values, jacobian, physical_jacobian)
      IF (PRESENT(targets)) values(:) = values - targets
   END SUBROUTINE qs_skala_scf_constraints

! **************************************************************************************************
!> \brief Restore simultaneous physical targets and solve tangential native orbital stationarity.
!>        Every constraint and energy trial reevaluates the original native fields and functional.
!>        The caller's selected limiting derivatives remain in effect. Convergence does not
!>        establish normal cusp admissibility, tangential stability or automatic active-set SCF.
!>        On failure retain the supplied coordinates. A stop or raw crossing returns the last
!>        accepted feasible coordinates. Automatic admission returns the verified downhill root.
!>        An explicit retain_progress request also keeps the last
!>        accepted state on other failures. In each case restore their native MOs, P, fields and
!>        KS operator. On success these quantities describe the returned stationary coordinates.
!> \param state ...
!> \param coordinate Initial coordinates, replaced on convergence, a stop or a crossing event.
!> \param points ...
!> \param domains ...
!> \param spins ...
!> \param descriptors ...
!> \param physical_tolerances Absolute physical tolerance per constraint, in its own units.
!> \param gradient_tolerance Absolute norm tolerance for the implicit free-coordinate derivative.
!> \param energy Complete native energy at the returned coordinates.
!> \param converged Whether simultaneous physical targets and tangential stationarity are established.
!> \param max_iterations Optional iteration limit.
!> \param max_step Optional bound on the free-coordinate step norm.
!> \param targets Optional fixed physical target per constraint, zero by default.
!> \param iterations Optional number of attempted orbital minimization iterations.
!> \param check_external Poll the normal CP2K stop conditions, false by default.
!> \param stopped Whether an external stop was requested, independent of convergence.
!> \param crossings Optional physical threshold scan along the current constrained orbital curve.
!>        Return all endpoint candidates and their trial segment before accepting the step.
!>        A detected crossing returns the last accepted coordinate with converged false.
!>        Root location, selection and minimum conditions must be established by the caller.
!> \param max_coordinate_norm Optional coordinate-norm bound for an unconstrained minimization.
!>        A trial exceeding this bound returns the accepted state for rebuilding its Cayley chart.
!> \param rechart Whether the coordinate bound or a stalled nonminimum active set requires caller handling.
!> \param retain_progress Preserve the last accepted state on failure, default false.
!> \param crossing_density_floor Optional positive spin-density target for an unconstrained crossing scan.
!>        Without an override the native XC dispatcher supplies all active spin-density floors.
!>        Raw kinetic thresholds are selected separately from the functional's projection.
!> \param admit_crossings Resolve continuous-floor candidates after failed line searches and admit downhill normal minima.
!>        Default false preserves the raw endpoint scan. Admission does not prove tangential stationarity.
!> \param output_unit Optional native SCF output unit for stored iteration scalars.
!> \param iteration_offset Previously consumed native SCF steps for the displayed counter.
! **************************************************************************************************
   SUBROUTINE qs_skala_scf_stationarity(state, coordinate, points, domains, spins, descriptors, &
                                        physical_tolerances, gradient_tolerance, energy, converged, &
                                        max_iterations, max_step, targets, iterations, check_external, stopped, crossings, &
                                        max_coordinate_norm, rechart, retain_progress, crossing_density_floor, admit_crossings, &
                                        output_unit, iteration_offset)
      TYPE(qs_skala_scf_type), INTENT(INOUT), TARGET     :: state
      REAL(dp), INTENT(INOUT)                            :: coordinate(:)
      INTEGER, INTENT(IN)                                :: points(:, :), domains(:), spins(:), &
                                                            descriptors(:)
      REAL(dp), INTENT(IN)                               :: physical_tolerances(:), &
                                                            gradient_tolerance
      REAL(dp), INTENT(OUT)                              :: energy
      LOGICAL, INTENT(OUT)                               :: converged
      INTEGER, INTENT(IN), OPTIONAL                      :: max_iterations
      REAL(dp), INTENT(IN), OPTIONAL                     :: max_step, targets(:)
      INTEGER, INTENT(OUT), OPTIONAL                     :: iterations
      LOGICAL, INTENT(IN), OPTIONAL                      :: check_external
      LOGICAL, INTENT(OUT), OPTIONAL                     :: stopped
      TYPE(qs_skala_scf_crossing_type), INTENT(OUT), &
         OPTIONAL, TARGET                                :: crossings
      REAL(dp), INTENT(IN), OPTIONAL                     :: max_coordinate_norm
      LOGICAL, INTENT(OUT), OPTIONAL                     :: rechart
      LOGICAL, INTENT(IN), OPTIONAL                      :: retain_progress
      REAL(dp), INTENT(IN), OPTIONAL                     :: crossing_density_floor
      LOGICAL, INTENT(IN), OPTIONAL                      :: admit_crossings
      INTEGER, INTENT(IN), OPTIONAL                      :: output_unit, iteration_offset

      INTEGER                                            :: i, m
      LOGICAL                                            :: coupled_curvature, valid
      REAL(dp), ALLOCATABLE                              :: coordinate_scales(:), density_floors(:), &
                                                            gradient(:), origin(:), response(:), &
                                                            zero(:), zero_multipliers(:)
      TYPE(mp_comm_type)                                 :: comm
      TYPE(physical_constraints_context)                 :: context
      TYPE(section_vals_type), POINTER                   :: functional_section, xc_section
      TYPE(xc_rho_cflags_type)                           :: needs

      m = SIZE(physical_tolerances)
      valid = SIZE(coordinate) == qs_skala_scf_size(state) .AND. ALL(SHAPE(points) == [3, m])
      CPASSERT(valid)
      valid = SIZE(domains) == m .AND. SIZE(spins) == m .AND. SIZE(descriptors) == m
      CPASSERT(valid)
      context%state => state
      context%gradient_tolerance = gradient_tolerance
      IF (PRESENT(output_unit)) context%iteration_unit = output_unit
      IF (PRESENT(iteration_offset)) context%iteration_offset = iteration_offset
      IF (PRESENT(admit_crossings)) context%admit_crossings = admit_crossings
      IF (PRESENT(check_external)) context%check_external = check_external
      IF (PRESENT(crossings)) THEN
         context%crossings => crossings
         xc_section => section_vals_get_subs_vals(state%qs_env%input, "DFT%XC")
         functional_section => section_vals_get_subs_vals(xc_section, "XC_FUNCTIONAL")
         needs = xc_functionals_get_needs(functional_section, SIZE(state%orbitals, 2) == 2, &
                                          kinetic_density_boundary=context%kinetic_boundaries)
         CPASSERT(.NOT. context%kinetic_boundaries .OR. needs%tau .OR. needs%tau_spin)
         ALLOCATE (context%density_floors, &
                   SOURCE=xc_functionals_spin_floors(functional_section, SIZE(state%orbitals, 2) == 2, &
                                                     section_get_rval(xc_section, "DENSITY_CUTOFF")))
         CALL comm%set_handle(state%para_env%get_handle())
         density_floors = native_skala_spin_floors(xc_section, comm)
         DO i = 1, SIZE(density_floors)
            IF (ALL(context%density_floors /= density_floors(i))) THEN
               context%density_floors = [context%density_floors, density_floors(i)]
            END IF
         END DO
         IF (context%kinetic_boundaries .AND. ALL(context%density_floors /= 0.0_dp)) THEN
            ! Native tau components retain their sign boundaries alongside spin-density floors.
            CALL MOVE_ALLOC(context%density_floors, density_floors)
            ALLOCATE (context%density_floors(SIZE(density_floors) + 1))
            context%density_floors(1) = 0.0_dp
            context%density_floors(2:) = density_floors
         END IF
      END IF
      IF (PRESENT(crossing_density_floor)) THEN
         CPASSERT(PRESENT(crossings))
         CPASSERT(SIZE(state%orbitals, 2) == 2)
         CPASSERT(ieee_is_finite(crossing_density_floor) .AND. crossing_density_floor > 0.0_dp)
         DEALLOCATE (context%density_floors)
         ALLOCATE (context%density_floors(1), SOURCE=crossing_density_floor)
      END IF
      IF (PRESENT(max_coordinate_norm)) THEN
         CPASSERT(m == 0)
         CPASSERT(ieee_is_finite(max_coordinate_norm) .AND. max_coordinate_norm > 0.0_dp)
         context%max_coordinate_norm = max_coordinate_norm
      END IF
      ALLOCATE (context%points(3, m), context%domains(m), context%spins(m), context%descriptors(m), context%targets(m))
      context%points(:, :) = points
      context%domains(:) = domains
      context%spins(:) = spins
      context%descriptors(:) = descriptors
      context%targets = 0.0_dp
      IF (PRESENT(targets)) THEN
         CPASSERT(SIZE(targets) == m)
         CPASSERT(ALL(ieee_is_finite(targets)))
         context%targets(:) = targets
      END IF
      origin = coordinate
      coupled_curvature = state%thermal
      IF (.NOT. coupled_curvature) THEN
         ! Probe the appropriate complete response before the native search stalls.
         ALLOCATE (zero(SIZE(coordinate)), response(SIZE(coordinate)))
         zero = 0.0_dp
         IF (m == 0) THEN
            CALL evaluate_native_hessian(coordinate, zero, context, response, coupled_curvature)
         ELSE
            ALLOCATE (zero_multipliers(m), SOURCE=0.0_dp)
            CALL evaluate_native_lagrangian_hessian(coordinate, zero, zero_multipliers, context, response, coupled_curvature)
         END IF
         DEALLOCATE (zero, response)
      END IF
      ALLOCATE (coordinate_scales(SIZE(coordinate)))
      CALL qs_skala_scf_coordinate_scales(state, coordinate_scales)
      CALL skala_constrained_stationarity(evaluate_native_constraints, evaluate_native_gradient, context, &
                                          coordinate, physical_tolerances, gradient_tolerance, converged, &
                                          max_iterations=max_iterations, max_step=max_step, &
                                          verify_stationarity=evaluate_native_stationarity, iterations=iterations, &
                                          check_stop=check_native_stop, stopped=stopped, &
                                          check_segment=check_native_step, retain_progress=retain_progress, &
                                          check_stalled_segment=check_native_stall, coordinate_scales=coordinate_scales, &
                                          curvature_steps=coupled_curvature, check_curvature_segment=check_native_curvature, &
                                          precondition=evaluate_native_preconditioner, apply_hessian=evaluate_native_hessian, &
                                          apply_lagrangian_hessian=evaluate_native_lagrangian_hessian, &
                                          report_iteration=report_native_iteration)
      IF (context%admit_crossings .AND. ASSOCIATED(context%crossings)) THEN
         IF (context%crossings%resolved) THEN
            ! Classification has already verified this downhill root in the native fields.
            coordinate(:) = context%crossings%root
            CPASSERT(at_native_coordinate(state, coordinate))
            CPASSERT(state%operator_current)
         END IF
      END IF
      IF (at_native_coordinate(state, coordinate) .AND. state%operator_current) THEN
         energy = state%energy%total
      ELSE
         ALLOCATE (gradient(SIZE(coordinate)))
         CALL evaluate_limiting_state(state, coordinate, energy, gradient, at_native_coordinate(state, coordinate))
      END IF
      IF (PRESENT(rechart)) THEN
         rechart = context%rechart
         IF (state%thermal .AND. ANY(coordinate /= origin)) rechart = .TRUE.
      END IF
      IF (PRESENT(stopped)) stopped = stopped .OR. context%external_stopped
   END SUBROUTINE qs_skala_scf_stationarity

! **************************************************************************************************
!> \brief Register physical cusp samples using the current native grid and atomic layout.
!>        Collective over QS. Geometry and sample row counts come from the actual native state.
!>        Selection changes derivatives only. Reevaluate before using the native KS operator.
!> \param state Evaluated native orbital state.
!> \param points Global FFT indices or (global atom, atomic row, 1).
!> \param domains Physical evaluation domains, in constraint order.
!> \param spins Zero for paired RKS, one or two for UKS.
!> \param descriptors Physical density or kinetic density.
!> \param negative_weights Negative-side convex weights.
!> \param density_probe Positive physical density displacement in bohr^-3.
!> \param density_tolerance Absolute physical density window in bohr^-3.
!> \param kinetic_probe Positive physical kinetic-density displacement in bohr^-5.
!> \param kinetic_tolerance Absolute physical kinetic-density window in bohr^-5.
!> \param targets Fixed physical descriptor values, zero by default.
! **************************************************************************************************
   SUBROUTINE qs_skala_scf_select(state, points, domains, spins, descriptors, negative_weights, &
                                  density_probe, density_tolerance, kinetic_probe, kinetic_tolerance, targets)
      TYPE(qs_skala_scf_type), INTENT(INOUT)             :: state
      INTEGER, INTENT(IN)                                :: points(:, :), domains(:), spins(:), &
                                                            descriptors(:)
      REAL(dp), INTENT(IN)                               :: negative_weights(:), density_probe, &
                                                            density_tolerance, kinetic_probe, &
                                                            kinetic_tolerance
      REAL(dp), INTENT(IN), OPTIONAL                     :: targets(:)

      INTEGER                                            :: atom, center, domain, i, ikind, m
      INTEGER, ALLOCATABLE                               :: atomic_sizes(:), keys(:, :), &
                                                            paw_sizes(:, :)
      LOGICAL                                            :: available, valid
      REAL(dp), ALLOCATABLE                              :: density(:, :), kinetic(:, :)
      TYPE(pw_grid_type), POINTER                        :: grid
      TYPE(pw_r3d_rs_type), POINTER                      :: rho_r(:)
      TYPE(qs_upf_inverse_type), POINTER                 :: frame

      CPASSERT(state%evaluated)
      m = SIZE(negative_weights)
      valid = ALL(SHAPE(points) == [3, m]) .AND. SIZE(domains) == m
      CPASSERT(valid)
      CPASSERT(SIZE(spins) == m .AND. SIZE(descriptors) == m)
      CALL qs_rho_get(state%rho, rho_r=rho_r)
      grid => rho_r(1)%pw_grid
      CPASSERT(ASSOCIATED(grid))
      ALLOCATE (atomic_sizes(SIZE(state%positions, 2)))
      atomic_sizes = 0
      DO domain = skala_fft_domain + 1, skala_background_domain
         CALL skala_density_samples_get(state%constraints, domain, keys, density, kinetic, available)
         IF (.NOT. available) CYCLE
         DO i = 1, SIZE(keys, 2)
            CPASSERT(keys(1, i) >= 1 .AND. keys(1, i) <= SIZE(atomic_sizes))
            atomic_sizes(keys(1, i)) = MAX(atomic_sizes(keys(1, i)), keys(2, i))
         END DO
      END DO
      CALL state%para_env%max(atomic_sizes)
      ALLOCATE (paw_sizes(2, SIZE(state%positions, 2)), SOURCE=0)
      frame => state%inverse%frame%ao
      IF (state%bloch) frame => state%inverse%bloch%frames(1)%ao
      DO center = 1, SIZE(frame%center_atoms)
         ikind = frame%center_kinds(center)
         IF (.NOT. state%inverse%paw(ikind)%active) CYCLE
         atom = frame%center_atoms(center)
         paw_sizes(:, atom) = [state%inverse%paw(ikind)%xcgrid%nr, state%inverse%paw(ikind)%xcgrid%na]
      END DO
      CALL skala_density_constraints_set(state%constraints, points, spins, negative_weights, &
                                         density_probe, density_tolerance, SIZE(state%orbitals, 2), &
                                         grid%npts, grid%bounds, state%hmat, state%positions, &
                                         descriptors=descriptors, kinetic_probe=kinetic_probe, &
                                         kinetic_tolerance=kinetic_tolerance, domains=domains, &
                                         atomic_sizes=atomic_sizes, targets=targets, paw_sizes=paw_sizes)
   END SUBROUTINE qs_skala_scf_select

! **************************************************************************************************
!> \brief Read owned constraints and verify their current physical feasibility.
!>        An empty selection is valid. Do not reuse a selection from another geometry or FFT layout.
!> \param state Evaluated native orbital state with its original physical sample keys.
!> \param points Global physical sample keys.
!> \param domains Evaluation domain per selected sample.
!> \param spins Physical spin selectors.
!> \param descriptors Density or kinetic density per sample.
!> \param targets Original prescribed physical values.
!> \param tolerances Original absolute physical feasibility windows.
!> \param valid Whether every selected physical value remains within its original window.
!> \param require_feasible Check current values, default true. False permits restoring the same
!>        owned constraints before minimization, but still requires the original geometry and layout.
! **************************************************************************************************
   SUBROUTINE qs_skala_scf_selected(state, points, domains, spins, descriptors, targets, tolerances, valid, require_feasible)
      TYPE(qs_skala_scf_type), INTENT(IN)                :: state
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: points(:, :), domains(:), spins(:), &
                                                            descriptors(:)
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: targets(:), tolerances(:)
      LOGICAL, INTENT(OUT)                               :: valid
      LOGICAL, INTENT(IN), OPTIONAL                      :: require_feasible

      INTEGER                                            :: m
      REAL(dp)                                           :: kinetic_probe, kinetic_tolerance, probe, &
                                                            tolerance
      REAL(dp), ALLOCATABLE                              :: values(:), weights(:)
      TYPE(pw_grid_type), POINTER                        :: grid
      TYPE(pw_r3d_rs_type), POINTER                      :: rho_r(:)

      CPASSERT(state%evaluated)
      m = skala_density_constraints_count(state%constraints)
      valid = .TRUE.
      IF (m == 0) THEN
         ALLOCATE (points(3, 0), domains(0), spins(0), descriptors(0), targets(0), tolerances(0))
         RETURN
      END IF
      CALL skala_density_constraints_get(state%constraints, points, spins, weights, probe, tolerance, &
                                         descriptors=descriptors, domains=domains, targets=targets, &
                                         kinetic_probe=kinetic_probe, kinetic_tolerance=kinetic_tolerance)
      ALLOCATE (values(m), tolerances(m))
      tolerances(:) = MERGE(tolerance, kinetic_tolerance, descriptors == skala_density_descriptor)
      CALL qs_rho_get(state%rho, rho_r=rho_r)
      grid => rho_r(1)%pw_grid
      valid = skala_density_constraints_matches(state%constraints, SIZE(state%orbitals, 2), &
                                                grid%npts, grid%bounds, state%hmat, state%positions)
      IF (.NOT. valid) RETURN
      IF (PRESENT(require_feasible)) THEN
         IF (.NOT. require_feasible) RETURN
      END IF
      CALL qs_skala_scf_values(state, domains, points, spins, descriptors, values)
      valid = ALL(ABS(values - targets) <= tolerances)
   END SUBROUTINE qs_skala_scf_selected

! **************************************************************************************************
!> \brief Evaluate complete native limiting covectors at the registered physical cusp samples.
!>        Columns switch one selected descriptor from its positive to its negative limit. All
!>        other native energy derivatives remain present. Covectors use the orthonormal physical
!>        orbital tangent, including occupations and k weights. Retain the original selection
!>        weights and, by default, restore the complete native KS evaluation at the supplied coordinates.
!> \param state Native state with the physical samples already registered.
!> \param coordinate Global normalized-orbital coordinates.
!> \param positive_gradient Complete all-positive physical limiting covector.
!> \param negative_changes Negative-minus-positive changes, one column per registered sample.
!> \param physical_jacobian Optional physical sample covectors, one constraint per row.
!> \param hold_orbitals Keep the evaluated native MOs and density fixed while selecting KS limits.
!> \param restore_operator Rebuild KS with the original weights, default true. Otherwise mark it
!>        invalid until the caller evaluates the restored selection.
!> \param assemble_jumps Compose qualified jumps from primitive covectors, default true.
!>        False retains the complete KS reference evaluation for every negative limit.
! **************************************************************************************************
   SUBROUTINE qs_skala_scf_limiting(state, coordinate, positive_gradient, negative_changes, physical_jacobian, &
                                    hold_orbitals, restore_operator, assemble_jumps)
      TYPE(qs_skala_scf_type), INTENT(INOUT)             :: state
      REAL(dp), INTENT(IN)                               :: coordinate(:)
      REAL(dp), INTENT(OUT)                              :: positive_gradient(:), &
                                                            negative_changes(:, :)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: physical_jacobian(:, :)
      LOGICAL, INTENT(IN), OPTIONAL                      :: hold_orbitals, restore_operator, &
                                                            assemble_jumps

      INTEGER                                            :: i, m, n
      INTEGER, ALLOCATABLE                               :: descriptors(:), domains(:), &
                                                            jump_owners(:), points(:, :), spins(:)
      LOGICAL                                            :: capture_jumps, fixed_orbitals, restore, &
                                                            valid
      LOGICAL, ALLOCATABLE                               :: assembled(:)
      REAL(dp)                                           :: energy, kinetic_probe, &
                                                            kinetic_tolerance, probe, tolerance, &
                                                            value
      REAL(dp), ALLOCATABLE                              :: gradient(:), jump_values(:), &
                                                            original_weights(:), physical(:), &
                                                            targets(:), weights(:)
      TYPE(section_vals_type), POINTER                   :: functionals

      n = qs_skala_scf_size(state)
      fixed_orbitals = .FALSE.
      IF (state%evaluated) THEN
         CPASSERT(SIZE(coordinate) == SIZE(state%coordinate))
         ! Selecting derivatives must preserve an already checked physical boundary.
         fixed_orbitals = ALL(state%coordinate == coordinate)
      END IF
      IF (PRESENT(hold_orbitals)) fixed_orbitals = hold_orbitals
      restore = .TRUE.
      IF (PRESENT(restore_operator)) restore = restore_operator
      m = skala_density_constraints_count(state%constraints)
      CPASSERT(SIZE(coordinate) == n .AND. SIZE(positive_gradient) == n)
      CPASSERT(ALL(SHAPE(negative_changes) == [n, m]))
      IF (PRESENT(physical_jacobian)) THEN
         valid = ALL(SHAPE(physical_jacobian) == [m, n])
         CPASSERT(valid)
      END IF
      negative_changes = 0.0_dp
      IF (m == 0) THEN
         CALL evaluate_limiting_state(state, coordinate, energy, positive_gradient, fixed_orbitals)
         RETURN
      END IF
      CALL skala_density_constraints_get(state%constraints, points, spins, original_weights, &
                                         probe, tolerance, descriptors=descriptors, domains=domains, &
                                         kinetic_probe=kinetic_probe, kinetic_tolerance=kinetic_tolerance, targets=targets)
      ALLOCATE (weights(m), gradient(n), physical(n), assembled(m))
      weights = 0.0_dp
      assembled = .FALSE.
      CALL skala_density_constraints_update_weights(state%constraints, weights)
      ! Primitive readouts must cover the complete active functional combination.
      ! Unsupported combinations and nonunique owners retain complete KS evaluations.
      capture_jumps = ALL(targets == 0.0_dp) .AND. skala_only_functional(state)
      IF (ALL(targets > 0.0_dp) .AND. ALL(descriptors == skala_density_descriptor)) THEN
         functionals => section_vals_get_subs_vals(state%qs_env%input, "DFT%XC%XC_FUNCTIONAL")
         capture_jumps = xc_functionals_floor_jumps_available(functionals)
      END IF
      IF (PRESENT(assemble_jumps)) capture_jumps = capture_jumps .AND. assemble_jumps
      IF (capture_jumps) CALL skala_density_jumps_begin(state%constraints)
      CALL evaluate_limiting_state(state, coordinate, energy, positive_gradient, fixed_orbitals)
      IF (capture_jumps) THEN
         CALL skala_density_jumps_finish(state%constraints, jump_values, jump_owners)
         CPASSERT(SIZE(jump_values) == m)
         CALL state%para_env%sum(jump_values)
         CALL state%para_env%sum(jump_owners)
         assembled = jump_owners == 1
      END IF
      ! Only the derivative selection changes from here. Keep the same installed MOs and density.
      fixed_orbitals = .TRUE.
      DO i = 1, m
         CALL qs_skala_scf_value(state, domains(i), points(:, i), spins(i), descriptors(i), value)
         IF (descriptors(i) == skala_density_descriptor) THEN
            CPASSERT(ABS(value - targets(i)) <= tolerance)
         ELSE
            CPASSERT(descriptors(i) == skala_kinetic_descriptor)
            CPASSERT(ABS(value - targets(i)) <= kinetic_tolerance)
         END IF
         IF (PRESENT(physical_jacobian) .OR. assembled(i)) THEN
            CALL qs_skala_scf_point_derivative(state, domains(i), points(:, i), spins(i), descriptors(i), &
                                               gradient, physical)
            IF (PRESENT(physical_jacobian)) physical_jacobian(i, :) = physical
            IF (assembled(i)) negative_changes(:, i) = jump_values(i)*physical
         END IF
      END DO
      DO i = 1, m
         IF (assembled(i)) CYCLE
         weights(i) = 1.0_dp
         CALL skala_density_constraints_update_weights(state%constraints, weights)
         CALL evaluate_limiting_state(state, coordinate, energy, physical, fixed_orbitals)
         negative_changes(:, i) = physical - positive_gradient
         weights(i) = 0.0_dp
      END DO
      CALL skala_density_constraints_update_weights(state%constraints, original_weights)
      IF (restore) THEN
         CALL evaluate_limiting_state(state, coordinate, energy, physical, fixed_orbitals)
      ELSE
         state%operator_current = .FALSE.
         CALL qs_ks_did_change(state%ks_env, rho_changed=.TRUE.)
      END IF
   END SUBROUTINE qs_skala_scf_limiting

! **************************************************************************************************
!> \brief Require a single active native model before composing complete jumps from its covectors.
!> \param state Native orbital state.
!> \return Whether GAUXC is the only active XC functional. Missing model covectors still fall back.
! **************************************************************************************************
   FUNCTION skala_only_functional(state) RESULT(skala_only)
      TYPE(qs_skala_scf_type), INTENT(IN)                :: state
      LOGICAL                                            :: skala_only

      INTEGER                                            :: i, nactive
      LOGICAL                                            :: active
      TYPE(section_vals_type), POINTER                   :: functional, functionals

      skala_only = .FALSE.
      functionals => section_vals_get_subs_vals(state%qs_env%input, "DFT%XC%XC_FUNCTIONAL")
      IF (.NOT. ASSOCIATED(functionals)) RETURN
      nactive = 0
      i = 0
      DO
         i = i + 1
         functional => section_vals_get_subs_vals2(functionals, i_section=i)
         IF (.NOT. ASSOCIATED(functional)) EXIT
         CALL section_vals_val_get(functional, "_SECTION_PARAMETERS_", l_val=active)
         IF (.NOT. active) CYCLE
         IF (functional%section%name /= "GAUXC") RETURN
         nactive = nactive + 1
      END DO
      skala_only = nactive == 1
   END FUNCTION skala_only_functional

! **************************************************************************************************
!> \brief Evaluate a selected KS limit, optionally preserving the installed native orbital state.
!>        A canonical state must be tested without installing newly normalized chart orbitals.
!>        Only the XC derivative selection changes in the held case. Rebuild KS and physical
!>        samples from the same native density, including all atomic corrections.
!> \param state Evaluated native orbital state.
!> \param coordinate Current orbital coordinates, unchanged when holding native orbitals.
!> \param energy Complete native energy.
!> \param physical_gradient Complete native KS covector in the normalized physical tangent.
!> \param hold_orbitals Reevaluate KS without changing MOs or rebuilding the density.
! **************************************************************************************************
   SUBROUTINE evaluate_limiting_state(state, coordinate, energy, physical_gradient, hold_orbitals)
      TYPE(qs_skala_scf_type), INTENT(INOUT)             :: state
      REAL(dp), INTENT(IN)                               :: coordinate(:)
      REAL(dp), INTENT(OUT)                              :: energy, physical_gradient(:)
      LOGICAL, INTENT(IN)                                :: hold_orbitals

      IF (.NOT. hold_orbitals) THEN
         CALL qs_skala_scf_evaluate(state, coordinate, energy, physical_gradient=physical_gradient)
         RETURN
      END IF
      CPASSERT(state%evaluated)
      CPASSERT(ALL(state%coordinate == coordinate))
      CALL skala_density_samples_begin(state%constraints)
      CALL qs_ks_did_change(state%ks_env, rho_changed=.TRUE.)
      CALL qs_ks_update_qs_env(state%qs_env, calculate_forces=.FALSE., just_energy=.FALSE., print_active=.FALSE.)
      CALL skala_density_samples_end(state%constraints)
      CALL retain_native_samples(state)
      energy = state%energy%total
      state%operator_current = .TRUE.
      CALL native_orbital_covectors(state, physical_gradient=physical_gradient)
   END SUBROUTINE evaluate_limiting_state

! **************************************************************************************************
!> \brief Determine convex weights from complete native physical limiting covectors.
!>        Leave native MOs, P, fields, energy and KS consistent with the resulting weights.
!>        Box optimality alone does not establish zero orbital gradient, normal admissibility,
!>        tangential stability or automatic SCF convergence.
!> \param state Registered physical cusp samples and their current convex weights.
!> \param coordinate Global normalized-orbital coordinates.
!> \param tolerance Absolute box optimality tolerance in physical tangent energy units.
!> \param balanced_gradient Complete native physical covector at the resulting weights.
!> \param balanced Whether the convex box optimality conditions meet tolerance.
! **************************************************************************************************
   SUBROUTINE qs_skala_scf_balance(state, coordinate, tolerance, balanced_gradient, balanced)
      TYPE(qs_skala_scf_type), INTENT(INOUT)             :: state
      REAL(dp), INTENT(IN)                               :: coordinate(:), tolerance
      REAL(dp), INTENT(OUT)                              :: balanced_gradient(:)
      LOGICAL, INTENT(OUT)                               :: balanced

      INTEGER                                            :: m, n
      REAL(dp)                                           :: energy
      REAL(dp), ALLOCATABLE                              :: changes(:, :), positive(:), predicted(:)

      n = qs_skala_scf_size(state)
      m = skala_density_constraints_count(state%constraints)
      CPASSERT(SIZE(balanced_gradient) == n)
      ALLOCATE (positive(n), changes(n, m), predicted(n))
      CALL qs_skala_scf_limiting(state, coordinate, positive, changes, restore_operator=.FALSE.)
      CALL skala_density_constraints_balance(state%constraints, positive, changes, tolerance, predicted, balanced)
      CALL evaluate_limiting_state(state, coordinate, energy, balanced_gradient, .TRUE.)
   END SUBROUTINE qs_skala_scf_balance

! **************************************************************************************************
!> \brief Check native limiting derivatives along attainable physical constraint normals.
!>        Retain the original weights and the complete native MOs, P, fields, energy and KS.
!>        Proportional samples switch together. Other dependent groups return invalid.
!>        Normal admissibility alone does not
!>        establish tangential stationarity, stability or automatic SCF convergence.
!> \param state Native state with registered physical cusp samples.
!> \param coordinate Global normalized-orbital coordinates.
!> \param admissible Whether the one-sided normal derivatives have the signs required for a minimum.
!> \param release_side Minus one for negative-side descent, plus one for positive-side descent.
!> \param valid Whether the selected normals have independent proportional groups.
!> \param hold_orbitals Keep native MOs and density fixed while evaluating the normal limits.
! **************************************************************************************************
   SUBROUTINE qs_skala_scf_admissibility(state, coordinate, admissible, release_side, valid, hold_orbitals)
      TYPE(qs_skala_scf_type), INTENT(INOUT)             :: state
      REAL(dp), INTENT(IN)                               :: coordinate(:)
      LOGICAL, INTENT(OUT)                               :: admissible(:)
      INTEGER, INTENT(OUT)                               :: release_side(:)
      LOGICAL, INTENT(OUT)                               :: valid
      LOGICAL, INTENT(IN), OPTIONAL                      :: hold_orbitals

      INTEGER                                            :: m, n
      REAL(dp), ALLOCATABLE                              :: changes(:, :), jacobian(:, :), &
                                                            positive(:)

      n = qs_skala_scf_size(state)
      m = skala_density_constraints_count(state%constraints)
      CPASSERT(SIZE(coordinate) == n)
      CPASSERT(SIZE(admissible) == m .AND. SIZE(release_side) == m)
      admissible = .FALSE.
      release_side = 0
      valid = .FALSE.
      ALLOCATE (positive(n), changes(n, m), jacobian(m, n))
      CALL qs_skala_scf_limiting(state, coordinate, positive, changes, jacobian, hold_orbitals)
      CALL skala_cusp_admissibility(positive, changes, jacobian, admissible, release_side, valid)
   END SUBROUTINE qs_skala_scf_admissibility

! **************************************************************************************************
!> \brief Release one nonminimum group of active physical samples from a feasible orbital state.
!>        Proportional samples leave together in the direction selected by their coupled limiting
!>        covectors. Restore every remaining physical constraint on its original implicit curve.
!>        Accept only a finite step outside the released windows that lowers the complete native
!>        energy beyond its rounding resolution. On failure restore the full original selection.
!> \param state Native normalized orbital chart and active physical samples.
!> \param coordinate Current tangentially stationary coordinates, replaced only on release.
!> \param max_step Original maximum orbital-coordinate step, or zero for a query without release.
!> \param released Whether a nonminimum group was removed and a descent step accepted.
!> \param valid Whether the normal test and any required finite release succeeded.
!> \param check_external Monitor original CP2K stop and runtime controls, default false.
!> \param stopped Whether an external stop was requested.
! **************************************************************************************************
   SUBROUTINE qs_skala_scf_release_cusp(state, coordinate, max_step, released, valid, check_external, stopped)
      TYPE(qs_skala_scf_type), INTENT(INOUT), TARGET     :: state
      REAL(dp), INTENT(INOUT)                            :: coordinate(:)
      REAL(dp), INTENT(IN)                               :: max_step
      LOGICAL, INTENT(OUT)                               :: released, valid
      LOGICAL, INTENT(IN), OPTIONAL                      :: check_external
      LOGICAL, INTENT(OUT), OPTIONAL                     :: stopped

      INTEGER                                            :: i, info, lwork, m, n, rank, shortening
      INTEGER, ALLOCATABLE                               :: descriptors(:), domains(:), groups(:), &
                                                            points(:, :), rows(:), sides(:), &
                                                            spins(:)
      LOGICAL                                            :: feasible, requested
      LOGICAL, ALLOCATABLE                               :: admissible(:), drop(:)
      REAL(dp)                                           :: energy, fraction, kinetic_probe, &
                                                            kinetic_tolerance, probe, query(1), &
                                                            roundoff, tolerance, trial_energy
      REAL(dp), ALLOCATABLE :: changes(:, :), direction(:), jacobian(:, :), matrix(:, :), &
         physical(:, :), positive(:), residual(:), rhs(:, :), scales(:), singular(:), targets(:), &
         tolerances(:), trial(:), values(:), weights(:), work(:)
      TYPE(physical_constraints_context)                 :: context
      TYPE(skala_constraint_chart_type)                  :: curve

      released = .FALSE.
      valid = .FALSE.
      IF (PRESENT(stopped)) stopped = .FALSE.
      n = qs_skala_scf_size(state)
      CPASSERT(SIZE(coordinate) == n .AND. ieee_is_finite(max_step) .AND. max_step >= 0.0_dp)
      CALL qs_skala_scf_evaluate(state, coordinate, energy)
      CALL qs_skala_scf_selected(state, points, domains, spins, descriptors, targets, tolerances, feasible)
      IF (.NOT. feasible) RETURN
      m = SIZE(spins)
      valid = .TRUE.
      IF (m == 0) RETURN
      ALLOCATE (positive(n), changes(n, m), physical(m, n), admissible(m), sides(m), groups(m), drop(m))
      CALL qs_skala_scf_limiting(state, coordinate, positive, changes, physical)
      CALL skala_cusp_admissibility(positive, changes, physical, admissible, sides, valid, group_ids=groups)
      IF (.NOT. valid .OR. ALL(admissible)) RETURN
      valid = .FALSE.
      IF (max_step == 0.0_dp) RETURN
      i = FINDLOC(admissible, .FALSE., DIM=1)
      IF (sides(i) == 0) RETURN
      drop(:) = groups == groups(i)
      ALLOCATE (values(m), jacobian(m, n), matrix(m, n), scales(m), &
                rhs(MAX(m, n), 1), singular(MIN(m, n)), direction(n), trial(n))
      CALL qs_skala_scf_constraints(state, coordinate, points, domains, spins, descriptors, values, jacobian, &
                                    targets=targets)
      scales(:) = SQRT(SUM(jacobian**2, DIM=2))
      IF (ANY(scales <= TINY(0.0_dp))) RETURN
      matrix(:, :) = jacobian/SPREAD(scales, 2, n)
      rhs = 0.0_dp
      WHERE (drop) rhs(:m, 1) = REAL(sides, dp)
      CALL dgelss(m, n, 1, matrix, m, rhs, MAX(m, n), singular, &
                  64.0_dp*EPSILON(0.0_dp), rank, query, -1, info)
      IF (info /= 0) RETURN
      lwork = MAX(1, CEILING(query(1)))
      ALLOCATE (work(lwork))
      CALL dgelss(m, n, 1, matrix, m, rhs, MAX(m, n), singular, &
                  64.0_dp*EPSILON(0.0_dp), rank, work, lwork, info)
      IF (info /= 0 .OR. rank /= MAXVAL(groups)) RETURN
      direction(:) = rhs(:n, 1)
      residual = MATMUL(jacobian/SPREAD(scales, 2, n), direction)
      WHERE (drop) residual = residual - REAL(sides, dp)
      IF (MAXVAL(ABS(residual)) > 128.0_dp*EPSILON(0.0_dp)) RETURN
      IF (NORM2(direction) <= TINY(0.0_dp)) RETURN
      direction(:) = direction/NORM2(direction)
      CALL skala_density_constraints_get(state%constraints, points, spins, weights, probe, tolerance, &
                                         descriptors=descriptors, domains=domains, targets=targets, &
                                         kinetic_probe=kinetic_probe, kinetic_tolerance=kinetic_tolerance)
      rows = PACK([(i, i=1, m)],.NOT. drop)
      IF (SIZE(rows) > 0) THEN
         CALL qs_skala_scf_select(state, points(:, rows), domains(rows), spins(rows), descriptors(rows), weights(rows), &
                                  probe, tolerance, kinetic_probe, kinetic_tolerance, targets=targets(rows))
      ELSE
         CALL skala_density_constraints_reset(state%constraints)
      END IF
      context%state => state
      IF (PRESENT(check_external)) context%check_external = check_external
      context%points = points(:, rows)
      context%domains = domains(rows)
      context%spins = spins(rows)
      context%descriptors = descriptors(rows)
      context%targets = targets(rows)
      DEALLOCATE (physical, residual)
      ALLOCATE (physical(SIZE(rows), n), residual(SIZE(rows)))
      CALL skala_constraint_chart_create(evaluate_native_constraints, context, coordinate, tolerances(rows), curve, feasible)
      fraction = max_step
      DO shortening = 1, DIGITS(1.0_dp)
         IF (.NOT. feasible) EXIT
         CALL check_native_stop(context, requested)
         IF (requested) THEN
            IF (PRESENT(stopped)) stopped = .TRUE.
            EXIT
         END IF
         trial(:) = coordinate + fraction*direction
         IF (ALL(trial == coordinate)) EXIT
         CALL skala_constraint_chart_restore(curve, evaluate_native_constraints, context, trial, residual, physical, feasible)
         IF (feasible) THEN
            ! An empty remaining curve does not evaluate a native trial. Refresh the
            ! released fields after clearing the final selection, before testing its side.
            IF (SIZE(rows) == 0) CALL sample_native_constraints(state, trial)
            CALL qs_skala_scf_values(state, domains, points, spins, descriptors, values)
            IF (ALL(REAL(sides, dp)*(values - targets) > tolerances .OR. .NOT. drop)) THEN
               CALL qs_skala_scf_evaluate(state, trial, trial_energy)
               roundoff = 64.0_dp*EPSILON(0.0_dp)*MAX(ABS(energy), ABS(trial_energy), TINY(0.0_dp))
               IF (trial_energy < energy - roundoff) THEN
                  coordinate(:) = trial
                  released = .TRUE.
                  valid = .TRUE.
                  RETURN
               END IF
            END IF
         END IF
         fraction = 0.5_dp*fraction
         feasible = .TRUE.
      END DO
      CALL qs_skala_scf_select(state, points, domains, spins, descriptors, weights, &
                               probe, tolerance, kinetic_probe, kinetic_tolerance, targets=targets)
      CALL qs_skala_scf_evaluate(state, coordinate, energy)
   END SUBROUTINE qs_skala_scf_release_cusp

! **************************************************************************************************
!> \brief Classify candidate continuous density boundaries with the existing active physical samples.
!>        Evaluate complete native limiting covectors and coupled normal minimum conditions.
!>        A zero covector change does not add an orbital constraint. Restore original selected
!>        weights and complete KS at coordinate without reinstalling an already checked native
!>        density. Optionally retain admissible new samples.
!>        These normal conditions do not establish energy continuity or tangential stationarity.
!> \param state Evaluated native state with its original selected samples.
!> \param coordinate Feasible boundary coordinates in the original normalized orbital chart.
!> \param points Candidate global physical sample keys.
!> \param domains Original evaluation domain of each candidate.
!> \param spins Original physical spin selector.
!> \param descriptors Density or kinetic-density selector.
!> \param targets Original physical target of each candidate.
!> \param admissible Whether a nonzero limiting change satisfies the coupled normal conditions.
!> \param release_side Negative or positive physical side permitting descent for a rejected candidate.
!> \param valid Whether all physical windows and the coupled normal grouping are valid.
!> \param retain Retain the candidates only if every one is admissible, default false.
!> \param restore_operator Restore complete KS after classification, default true.
!>        False leaves KS invalid until a subsequent complete native evaluation.
! **************************************************************************************************
   SUBROUTINE qs_skala_scf_classify(state, coordinate, points, domains, spins, descriptors, targets, &
                                    admissible, release_side, valid, retain, restore_operator)
      TYPE(qs_skala_scf_type), INTENT(INOUT)             :: state
      REAL(dp), INTENT(IN)                               :: coordinate(:)
      INTEGER, INTENT(IN)                                :: points(:, :), domains(:), spins(:), &
                                                            descriptors(:)
      REAL(dp), INTENT(IN)                               :: targets(:)
      LOGICAL, INTENT(OUT)                               :: admissible(:)
      INTEGER, INTENT(OUT)                               :: release_side(:)
      LOGICAL, INTENT(OUT)                               :: valid
      LOGICAL, INTENT(IN), OPTIONAL                      :: retain, restore_operator

      INTEGER                                            :: i, m, n, old, total
      INTEGER, ALLOCATABLE :: all_descriptors(:), all_domains(:), all_points(:, :), all_spins(:), &
         old_descriptors(:), old_domains(:), old_points(:, :), old_spins(:), rows(:), sides(:)
      LOGICAL                                            :: feasible, hold_orbitals, keep, restore
      LOGICAL, ALLOCATABLE                               :: normal_admissible(:), selected(:)
      REAL(dp)                                           :: energy, kinetic_probe, &
                                                            kinetic_tolerance, probe, tolerance
      REAL(dp), ALLOCATABLE :: all_targets(:), changes(:, :), jacobian(:, :), old_targets(:), &
         old_weights(:), positive(:), values(:), weights(:)

      n = qs_skala_scf_size(state)
      m = SIZE(targets)
      valid = SIZE(coordinate) == n .AND. ALL(SHAPE(points) == [3, m])
      CPASSERT(valid)
      CPASSERT(SIZE(domains) == m)
      CPASSERT(SIZE(spins) == m)
      CPASSERT(SIZE(descriptors) == m)
      CPASSERT(SIZE(admissible) == m)
      CPASSERT(SIZE(release_side) == m)
      admissible = .FALSE.
      release_side = 0
      valid = .FALSE.
      keep = .FALSE.
      IF (PRESENT(retain)) keep = retain
      restore = .TRUE.
      IF (PRESENT(restore_operator)) restore = restore_operator
      ALLOCATE (positive(n))
      hold_orbitals = at_native_coordinate(state, coordinate)
      IF (hold_orbitals .AND. state%operator_current) THEN
         energy = state%energy%total
      ELSE
         CALL evaluate_limiting_state(state, coordinate, energy, positive, hold_orbitals)
      END IF
      CALL qs_skala_scf_selected(state, old_points, old_domains, old_spins, old_descriptors, &
                                 old_targets, values, feasible)
      IF (.NOT. feasible) RETURN
      old = SIZE(old_spins)
      ! These are the qualified absolute native density windows, not relative energy tolerances.
      probe = 1.0e-14_dp
      tolerance = 1.0e-17_dp
      kinetic_probe = 1.0e-14_dp
      kinetic_tolerance = 1.0e-17_dp
      IF (old > 0) THEN
         CALL skala_density_constraints_get(state%constraints, old_points, old_spins, old_weights, probe, tolerance, &
                                            descriptors=old_descriptors, domains=old_domains, targets=old_targets, &
                                            kinetic_probe=kinetic_probe, kinetic_tolerance=kinetic_tolerance)
      ELSE
         ALLOCATE (old_weights(0))
      END IF
      DEALLOCATE (values)
      ALLOCATE (values(m))
      CALL qs_skala_scf_values(state, domains, points, spins, descriptors, values)
      feasible = ALL(ABS(values - targets) <= MERGE(tolerance, kinetic_tolerance, descriptors == skala_density_descriptor))
      IF (.NOT. feasible) RETURN
      total = old + m
      ALLOCATE (all_points(3, total), all_domains(total), all_spins(total), all_descriptors(total), &
                all_targets(total), weights(total), normal_admissible(total), sides(total), &
                changes(n, total), jacobian(total, n), selected(total))
      all_points(:, :old) = old_points
      all_points(:, old + 1:) = points
      all_domains(:old) = old_domains
      all_domains(old + 1:) = domains
      all_spins(:old) = old_spins
      all_spins(old + 1:) = spins
      all_descriptors(:old) = old_descriptors
      all_descriptors(old + 1:) = descriptors
      all_targets(:old) = old_targets
      all_targets(old + 1:) = targets
      weights = 0.0_dp
      weights(:old) = old_weights
      CALL qs_skala_scf_select(state, all_points, all_domains, all_spins, all_descriptors, weights, &
                               probe, tolerance, kinetic_probe, kinetic_tolerance, targets=all_targets)
      CALL qs_skala_scf_limiting(state, coordinate, positive, changes, jacobian, restore_operator=.FALSE.)
      CALL skala_cusp_admissibility(positive, changes, jacobian, normal_admissible, sides, valid)
      IF (valid) THEN
         admissible(:) = normal_admissible(old + 1:)
         release_side(:) = sides(old + 1:)
         DO i = 1, m
            IF (ALL(changes(:, old + i) == 0.0_dp)) admissible(i) = .FALSE.
         END DO
      END IF
      selected = .FALSE.
      selected(:old) = .TRUE.
      IF (valid .AND. keep .AND. ALL(admissible)) selected(old + 1:) = .TRUE.
      ALLOCATE (rows, SOURCE=PACK([(i, i=1, total)], selected))
      IF (SIZE(rows) > 0) THEN
         CALL qs_skala_scf_select(state, all_points(:, rows), all_domains(rows), all_spins(rows), &
                                  all_descriptors(rows), weights(rows), probe, tolerance, kinetic_probe, kinetic_tolerance, &
                                  targets=all_targets(rows))
      ELSE
         CALL skala_density_constraints_clear_selection(state%constraints)
      END IF
      IF (restore) THEN
         CALL evaluate_limiting_state(state, coordinate, energy, positive, .TRUE.)
      ELSE
         state%operator_current = .FALSE.
         CALL qs_ks_did_change(state%ks_env, rho_changed=.TRUE.)
      END IF
   END SUBROUTINE qs_skala_scf_classify

! **************************************************************************************************
!> \brief Accept a stationary native trial with matching canonical MOs, P, W, fields, energy and KS.
!>        Check complete physical stationarity and the independent cusp-normal minimum conditions
!>        before changing native eigenvalues or W. On rejection the evaluated trial remains live.
!>        Recheck the actual canonical state in a fresh normalized chart. A rejected canonical
!>        state remains in the native environment so a caller can rechart and refine it.
!>        Success consumes the orbital charts because canonicalization changes their reference.
!>        This does not establish tangential stability or perform an automatic SCF search.
!> \param state Native charts, released only on successful acceptance.
!> \param coordinate Replicated coordinates of the trial to accept.
!> \param tolerance Bound on the complete physical orbital-gradient norm.
!> \param accepted Whether all blocks were canonicalized and the native state was rebuilt.
!> \param omitted_charge Optional thermal electron count outside retained MOs, zero until checked.
! **************************************************************************************************
   SUBROUTINE qs_skala_scf_accept(state, coordinate, tolerance, accepted, omitted_charge)
      TYPE(qs_skala_scf_type), INTENT(INOUT)             :: state
      REAL(dp), INTENT(IN)                               :: coordinate(:), tolerance
      LOGICAL, INTENT(OUT)                               :: accepted
      REAL(dp), INTENT(OUT), OPTIONAL                    :: omitted_charge

      COMPLEX(dp), ALLOCATABLE                           :: hamiltonian(:, :)
      INTEGER                                            :: component, failed, ik, local, m, spin
      INTEGER, ALLOCATABLE                               :: release_side(:), &
                                                            selected_descriptors(:), &
                                                            selected_domains(:), &
                                                            selected_points(:, :), &
                                                            selected_spins(:)
      LOGICAL                                            :: valid
      LOGICAL, ALLOCATABLE                               :: admissible(:)
      REAL(dp)                                           :: energy, fermi_error, omitted_number
      REAL(dp), ALLOCATABLE                              :: canonical_coordinate(:), gradient(:), &
                                                            physical_mu(:), selected_targets(:), &
                                                            selected_tolerances(:)
      TYPE(dbcsr_p_type), POINTER                        :: densities(:, :), matrices(:, :), &
                                                            weighted(:)
      TYPE(kpoint_env_type), POINTER                     :: kp
      TYPE(qs_skala_scf_type)                            :: canonical

      accepted = .FALSE.
      IF (PRESENT(omitted_charge)) omitted_charge = 0.0_dp
      CPASSERT(ieee_is_finite(tolerance) .AND. tolerance > 0.0_dp)
      ALLOCATE (gradient(qs_skala_scf_size(state)))
      CALL qs_skala_scf_evaluate(state, coordinate, energy, physical_gradient=gradient)
      IF (NORM2(gradient) > tolerance) RETURN
      m = skala_density_constraints_count(state%constraints)
      IF (m > 0) THEN
         ALLOCATE (admissible(m), release_side(m))
         CALL qs_skala_scf_admissibility(state, coordinate, admissible, release_side, valid)
         IF (.NOT. valid) RETURN
         IF (.NOT. ALL(admissible)) RETURN
      END IF
      IF (.NOT. state%bloch) CALL get_ks_env(state%ks_env, matrix_ks_kp=matrices)
      failed = 0
      DO spin = 1, SIZE(state%orbitals, 2)
         DO local = 1, SIZE(state%orbitals, 1)
            IF (.NOT. ASSOCIATED(state%orbitals(local, spin)%native)) CYCLE
            IF (state%bloch) THEN
               ik = state%kpoints%kp_env(local)%kpoint_env%nkpoint
               CALL native_bloch_operator(state, state%inverse%ks_k(:, :, :, spin), ik, hamiltonian)
            ELSE
               CALL qs_upf_inverse_matrix_gather(state%inverse%frame%ao, matrices(spin, 1)%matrix, hamiltonian)
            END IF
            CALL qs_skala_orbitals_canonicalize(state%orbitals(local, spin), hamiltonian, tolerance, valid)
            IF (.NOT. valid) failed = 1
         END DO
      END DO
      CALL state%para_env%max(failed)
      IF (failed /= 0) THEN
         CALL qs_skala_scf_evaluate(state, coordinate, energy, physical_gradient=gradient)
         RETURN
      END IF
      state%inverse%orbital_factors%density_evaluation = -1
      IF (.NOT. state%bloch) THEN
         CALL qs_rho_get(state%rho, rho_ao_kp=densities)
         CALL get_ks_env(state%ks_env, matrix_w=weighted)
         IF (.NOT. ASSOCIATED(weighted)) THEN
            CALL dbcsr_allocate_matrix_set(weighted, SIZE(state%orbitals, 2))
            DO spin = 1, SIZE(weighted)
               ALLOCATE (weighted(spin)%matrix)
               CALL dbcsr_copy(weighted(spin)%matrix, densities(spin, 1)%matrix, name="W MATRIX")
            END DO
            CALL set_ks_env(state%ks_env, matrix_w=weighted)
         END IF
         DO spin = 1, SIZE(weighted)
            CALL dbcsr_set(weighted(spin)%matrix, 0.0_dp)
         END DO
      ELSE
         DO local = 1, SIZE(state%kpoints%kp_env)
            kp => state%kpoints%kp_env(local)%kpoint_env
            DO spin = 1, SIZE(kp%wmat, 2)
               DO component = 1, SIZE(kp%wmat, 1)
                  CALL cp_fm_set_all(kp%wmat(component, spin), 0.0_dp)
               END DO
            END DO
         END DO
      END IF
      DO spin = 1, SIZE(state%orbitals, 2)
         DO local = 1, SIZE(state%orbitals, 1)
            IF (.NOT. ASSOCIATED(state%orbitals(local, spin)%native)) CYCLE
            IF (state%bloch) THEN
               CALL qs_skala_orbitals_install_kpoint(state%orbitals(local, spin), state%kpoints, local, spin, &
                                                     energy_weighted=.TRUE.)
            ELSE
               CALL qs_skala_orbitals_install_gamma(state%orbitals(local, spin), state%inverse%frame%ao, &
                                                    densities(spin, 1)%matrix, weighted(spin)%matrix)
            END IF
         END DO
      END DO
      CALL qs_rho_set(state%rho, rho_r_valid=.FALSE., rho_g_valid=.FALSE., tau_r_valid=.FALSE., tau_g_valid=.FALSE.)
      CALL qs_rho_update_rho(state%rho, state%qs_env)
      CALL retain_native_orbital_factors(state)
      CALL skala_density_samples_begin(state%constraints)
      CALL qs_ks_did_change(state%ks_env, rho_changed=.TRUE.)
      CALL qs_ks_update_qs_env(state%qs_env, calculate_forces=.FALSE., just_energy=.FALSE., print_active=.FALSE.)
      CALL skala_density_samples_end(state%constraints)
      CALL qs_skala_scf_create(canonical, state%qs_env, valid, thermal=state%thermal)
      IF (valid .AND. state%thermal) THEN
         ! Keep auxiliary energies and occupations while native eigenvalues and W use H.
         CPASSERT(ALL(SHAPE(canonical%band_dimensions) == SHAPE(state%band_dimensions)))
         CPASSERT(ALL(canonical%band_dimensions == state%band_dimensions))
         canonical%band_reference = state%band_energies
         CALL update_thermal_state(canonical, canonical%coordinate(canonical%orbital_size + 1:))
      END IF
      IF (valid) THEN
         ALLOCATE (canonical_coordinate(qs_skala_scf_size(canonical)), SOURCE=0.0_dp)
         CALL qs_skala_scf_inspect(canonical, energy, gradient)
         CALL qs_skala_scf_selected(canonical, selected_points, selected_domains, selected_spins, &
                                    selected_descriptors, selected_targets, selected_tolerances, valid)
         valid = valid .AND. NORM2(gradient) <= tolerance
         IF (valid .AND. m > 0) THEN
            CALL qs_skala_scf_admissibility(canonical, canonical_coordinate, admissible, release_side, valid, &
                                            hold_orbitals=.TRUE.)
            valid = valid .AND. ALL(admissible)
            CALL native_orbital_covectors(canonical, physical_gradient=gradient)
            valid = valid .AND. NORM2(gradient) <= tolerance
         END IF
      END IF
      IF (valid .AND. state%thermal) THEN
         CALL thermal_density_stationarity(canonical, tolerance, fermi_error, omitted_number, physical_mu, valid)
         IF (PRESENT(omitted_charge)) omitted_charge = omitted_number
         IF (valid) THEN
            DO spin = 1, SIZE(canonical%orbitals, 2)
               DO local = 1, SIZE(canonical%orbitals, 1)
                  IF (.NOT. ASSOCIATED(canonical%orbitals(local, spin)%native)) CYCLE
                  ik = 1
                  IF (canonical%bloch) ik = canonical%kpoints%kp_env(local)%kpoint_env%nkpoint
                  component = canonical%band_groups(canonical%band_offsets(ik, spin) + 1)
                  canonical%orbitals(local, spin)%native%mu = physical_mu(component)
               END DO
            END DO
         END IF
      END IF
      CALL qs_skala_scf_release(canonical)
      IF (.NOT. valid) THEN
         state%evaluated = .FALSE.
         RETURN
      END IF
      CALL qs_skala_scf_release(state)
      accepted = .TRUE.
   END SUBROUTINE qs_skala_scf_accept

! **************************************************************************************************
!> \brief Evaluate the original physical residuals and their coordinate derivatives.
!> \param coordinate Complete normalized orbital coordinates.
!> \param values Physical residuals relative to the supplied targets.
!> \param jacobian Derivative of the physical residuals in orbital coordinates.
!> \param context Borrowed native orbital and constraint state.
! **************************************************************************************************
   SUBROUTINE evaluate_native_constraints(coordinate, values, jacobian, context)
      REAL(dp), INTENT(IN)                               :: coordinate(:)
      REAL(dp), INTENT(OUT)                              :: values(:), jacobian(:, :)

      CLASS(*), INTENT(INOUT)                            :: context

      ! With no active constraints the gradient callback supplies the complete native trial.
      IF (SIZE(values) == 0) RETURN
      SELECT TYPE (context)
      TYPE IS (physical_constraints_context)
         CALL qs_skala_scf_constraints(context%state, coordinate, context%points, context%domains, &
                                       context%spins, context%descriptors, values, jacobian, targets=context%targets)
      CLASS DEFAULT
         CPABORT("Invalid native Skala constraint context.")
      END SELECT
   END SUBROUTINE evaluate_native_constraints

! **************************************************************************************************
!> \brief Evaluate the complete native energy and its coordinate derivative.
!> \param coordinate Complete normalized orbital coordinates.
!> \param energy Complete native energy of the trial.
!> \param gradient Energy derivative in the same orbital coordinates.
!> \param context Borrowed native orbital and constraint state.
! **************************************************************************************************
   SUBROUTINE evaluate_native_gradient(coordinate, energy, gradient, context)
      REAL(dp), INTENT(IN)                               :: coordinate(:)
      REAL(dp), INTENT(OUT)                              :: energy, gradient(:)

      CLASS(*), INTENT(INOUT)                            :: context

      SELECT TYPE (context)
      TYPE IS (physical_constraints_context)
         IF (SIZE(context%points, 2) > 0 .AND. at_native_coordinate(context%state, coordinate)) THEN
            CALL evaluate_limiting_state(context%state, coordinate, energy, gradient, .TRUE.)
            CALL native_orbital_covectors(context%state, gradient=gradient)
         ELSE
            CALL qs_skala_scf_evaluate(context%state, coordinate, energy, gradient, &
                                       sample_inputs=SIZE(context%points, 2) > 0)
         END IF
      CLASS DEFAULT
         CPABORT("Invalid native Skala gradient context.")
      END SELECT
   END SUBROUTINE evaluate_native_gradient

! **************************************************************************************************
!> \brief Verify physical orbital stationarity, retaining only tangent motion for constrained trials.
!>        Project in the orthonormal physical orbital metric, using fresh sample covectors.
!>        Project onto the span of every original normal, including dependent physical samples.
!> \param coordinate Complete normalized orbital coordinates.
!> \param context Native orbital and constraint state.
!> \param converged Whether the physical gradient meets the original tolerance.
! **************************************************************************************************
   SUBROUTINE evaluate_native_stationarity(coordinate, context, converged)
      REAL(dp), INTENT(IN)                               :: coordinate(:)

      CLASS(*), INTENT(INOUT) :: context
      LOGICAL, INTENT(OUT) :: converged

      INTEGER :: info, lwork, m, n, rank
      REAL(dp) :: energy, query(1)
      REAL(dp), ALLOCATABLE :: coefficients(:, :), gradient(:), jacobian(:, :), matrix(:, :), &
                               normals(:, :), physical(:, :), scales(:), singular(:), values(:), work(:)

      converged = .FALSE.
      SELECT TYPE (context)
      TYPE IS (physical_constraints_context)
         n = SIZE(coordinate)
         m = SIZE(context%points, 2)
         ALLOCATE (gradient(n))
         IF (m > 0) THEN
            ALLOCATE (values(m), jacobian(m, n), physical(m, n), scales(m), &
                      matrix(n, m), normals(n, m), coefficients(MAX(n, m), 1), singular(MIN(n, m)))
            CALL qs_skala_scf_constraints(context%state, coordinate, context%points, context%domains, &
                                          context%spins, context%descriptors, values, jacobian, physical, context%targets)
            scales(:) = MAX(MAXVAL(ABS(physical), DIM=2), TINY(0.0_dp))
            normals(:, :) = TRANSPOSE(physical/SPREAD(scales, 2, n))
         END IF
         IF (m > 0) THEN
            CALL evaluate_limiting_state(context%state, coordinate, energy, gradient, .TRUE.)
         ELSE
            CALL qs_skala_scf_evaluate(context%state, coordinate, energy, physical_gradient=gradient, sample_inputs=.FALSE.)
         END IF
         IF (m > 0) THEN
            ! The fitted normal component is removed from the physical covector, not from
            ! its Cayley-coordinate derivative. Scaling each constraint preserves its units.
            matrix(:, :) = normals
            coefficients = 0.0_dp
            coefficients(:n, 1) = gradient
            CALL dgelss(n, m, 1, matrix, n, coefficients, MAX(n, m), singular, &
                        64.0_dp*EPSILON(0.0_dp), rank, query, -1, info)
            IF (info /= 0) RETURN
            lwork = MAX(1, CEILING(query(1)))
            ALLOCATE (work(lwork))
            CALL dgelss(n, m, 1, matrix, n, coefficients, MAX(n, m), singular, &
                        64.0_dp*EPSILON(0.0_dp), rank, work, lwork, info)
            IF (info /= 0) RETURN
            gradient(:) = gradient - MATMUL(normals, coefficients(:m, 1))
         END IF
         converged = NORM2(gradient) <= context%gradient_tolerance
      CLASS DEFAULT
         CPABORT("Invalid native orbital convergence context.")
      END SELECT
   END SUBROUTINE evaluate_native_stationarity

! **************************************************************************************************
!> \brief Print accepted-state scalars without changing fields, orbitals or operators.
!>        The reduced covector norm is distinct from the final physical stationarity test.
!> \param iteration Attempted step in the current orbital chart, zero at its reference.
!> \param energy Last accepted energy.
!> \param gradient_norm Original-unit independent-coordinate covector norm.
!> \param moved Whether a new feasible step was accepted.
!> \param context Native state, original step offset and existing SCF print unit.
! **************************************************************************************************
   SUBROUTINE report_native_iteration(iteration, energy, gradient_norm, moved, context)
      INTEGER, INTENT(IN)                                :: iteration
      REAL(dp), INTENT(IN)                               :: energy, gradient_norm
      LOGICAL, INTENT(IN)                                :: moved

      CLASS(*), INTENT(INOUT) :: context
      CHARACTER(LEN=6) :: action

      SELECT TYPE (context)
      TYPE IS (physical_constraints_context)
         IF (context%iteration_unit < 0) RETURN
         IF (context%state%para_env%mepos /= 0) RETURN
         action = 'reject'
         IF (moved) action = 'accept'
         IF (iteration == 0) THEN
            action = 'start'
            WRITE (context%iteration_unit, '(A)') &
               ' INVERSE_SCF|   Step State             Energy [Ha]     Reduced gradient [Ha] Constraints'
         END IF
         WRITE (context%iteration_unit, '(A,I6,1X,A6,2ES24.14,1X,I8)') ' INVERSE_SCF| ', &
            context%iteration_offset + iteration, action, energy, gradient_norm, SIZE(context%points, 2)
         FLUSH (context%iteration_unit)
      END SELECT
   END SUBROUTINE report_native_iteration

! **************************************************************************************************
!> \brief Poll normal CP2K stop files, external communication and the target runtime.
!> \param context Native orbital state and optional stop monitoring.
!> \param requested Whether an external stop was requested.
! **************************************************************************************************
   SUBROUTINE check_native_stop(context, requested)
      CLASS(*), INTENT(INOUT) :: context
      LOGICAL, INTENT(OUT) :: requested

      requested = .FALSE.
      SELECT TYPE (context)
      TYPE IS (physical_constraints_context)
         IF (context%check_external) THEN
            CALL external_control(requested, "SCF", target_time=context%state%qs_env%target_time, &
                                  start_time=context%state%qs_env%start_time)
         END IF
      CLASS DEFAULT
         CPABORT("Invalid native orbital stop context.")
      END SELECT
   END SUBROUTINE check_native_stop

! **************************************************************************************************
!> \brief Check finite chart bounds before accepting an ordinary native orbital step.
!>        Raw crossing requests still inspect each trial. Automatic admission waits for a
!>        failed line search or an unresolved energy change above the gradient tolerance.
!> \param first Last accepted feasible coordinates.
!> \param last Proposed feasible coordinates.
!> \param curve Original independent/dependent coordinate partition.
!> \param context Native state and caller event output.
!> \param interrupted Whether the caller must rebuild the chart or handle a raw event.
! **************************************************************************************************
   SUBROUTINE check_native_step(first, last, curve, context, interrupted)
      REAL(dp), INTENT(IN)                               :: first(:), last(:)
      TYPE(skala_constraint_chart_type), INTENT(IN)      :: curve

      CLASS(*), INTENT(INOUT) :: context
      LOGICAL, INTENT(OUT) :: interrupted

      interrupted = .FALSE.
      SELECT TYPE (context)
      TYPE IS (physical_constraints_context)
         IF (.NOT. context%admit_crossings) THEN
            CALL check_native_segment(first, last, curve, context, interrupted)
         ELSE IF (context%max_coordinate_norm > 0.0_dp) THEN
            IF (NORM2(last) > context%max_coordinate_norm) THEN
               context%rechart = .TRUE.
               interrupted = .TRUE.
            END IF
         END IF
      CLASS DEFAULT
         CPABORT("Invalid native orbital step context.")
      END SELECT
   END SUBROUTINE check_native_step

! **************************************************************************************************
!> \brief Inspect the original feasible trial when descent fails or its energy change is unresolved.
!> \param first Last accepted feasible coordinates.
!> \param last First feasible trial before backtracking reduces its interval.
!> \param curve Original independent/dependent coordinate partition.
!> \param context Native state and physical boundary admission settings.
!> \param interrupted Whether a physical boundary or external stop requires the caller.
! **************************************************************************************************
   SUBROUTINE check_native_stall(first, last, curve, context, interrupted)
      REAL(dp), INTENT(IN)                               :: first(:), last(:)
      TYPE(skala_constraint_chart_type), INTENT(IN)      :: curve

      CLASS(*), INTENT(INOUT) :: context
      LOGICAL, INTENT(OUT) :: interrupted

      INTEGER, ALLOCATABLE :: release_side(:)
      LOGICAL, ALLOCATABLE :: admissible(:)
      LOGICAL :: valid

      interrupted = .FALSE.
      SELECT TYPE (context)
      TYPE IS (physical_constraints_context)
         IF (context%admit_crossings) THEN
            CALL check_native_segment(first, last, curve, context, interrupted)
            IF (.NOT. interrupted .AND. SIZE(context%spins) > 0) THEN
               ! A selected normal can cease to bracket a minimum before the tangent converges.
               ! Return the accepted state so the caller can test a finite downhill release.
               ALLOCATE (admissible(SIZE(context%spins)), release_side(SIZE(context%spins)))
               CALL qs_skala_scf_admissibility(context%state, first, admissible, release_side, valid)
               IF (valid .AND. ANY(.NOT. admissible)) THEN
                  context%rechart = .TRUE.
                  interrupted = .TRUE.
               END IF
            END IF
         END IF
      CLASS DEFAULT
         CPABORT("Invalid native orbital stall context.")
      END SELECT
   END SUBROUTINE check_native_stall

! **************************************************************************************************
!> \brief Apply the complete physical free-energy response at the accepted unconstrained native state.
!> \param coordinate Complete orbital and Fermi coordinates.
!> \param direction Direction in those same coordinates.
!> \param context Borrowed native state and original physical constraints.
!> \param PRODUCT Complete free-energy Hessian product.
!> \param available Whether all original electronic contributions have a response.
! **************************************************************************************************
   SUBROUTINE evaluate_native_hessian(coordinate, direction, context, PRODUCT, available)
      REAL(dp), INTENT(IN)                               :: coordinate(:), direction(:)

      CLASS(*), INTENT(INOUT)                            :: context
      REAL(dp), INTENT(OUT)                              :: PRODUCT(:)
      LOGICAL, INTENT(OUT)                               :: available

      REAL(dp)                                           :: energy
      REAL(dp), ALLOCATABLE                              :: gradient(:)

      available = .FALSE.
      PRODUCT = 0.0_dp
      SELECT TYPE (context)
      TYPE IS (physical_constraints_context)
         IF (SIZE(context%points, 2) /= 0) RETURN
         IF (.NOT. at_native_coordinate(context%state, coordinate) .OR. .NOT. context%state%operator_current) THEN
            ALLOCATE (gradient(SIZE(coordinate)))
            CALL qs_skala_scf_evaluate(context%state, coordinate, energy, gradient)
         END IF
         CALL qs_skala_scf_curvature(context%state, direction, PRODUCT, available)
      CLASS DEFAULT
         CPABORT("Invalid native Skala Hessian context.")
      END SELECT
   END SUBROUTINE evaluate_native_hessian

! **************************************************************************************************
!> \brief Supply the native original-unit density-floor Lagrangian response.
!> \param coordinate Accepted coordinates on the physical constraint curve.
!> \param direction Complete feasible coordinate direction.
!> \param multipliers Constraint coefficients in original physical units.
!> \param context Borrowed state and the same registered physical constraint ordering.
!> \param PRODUCT Full response before the caller's final tangent projection.
!> \param available Whether a complete native tangent response is supplied.
! **************************************************************************************************
   SUBROUTINE evaluate_native_lagrangian_hessian(coordinate, direction, multipliers, context, PRODUCT, available)
      REAL(dp), INTENT(IN)                               :: coordinate(:), direction(:), &
                                                            multipliers(:)

      CLASS(*), INTENT(INOUT)                            :: context
      REAL(dp), INTENT(OUT)                              :: PRODUCT(:)
      LOGICAL, INTENT(OUT)                               :: available

      INTEGER, ALLOCATABLE                               :: descriptors(:), domains(:), points(:, :), &
                                                            spins(:)
      LOGICAL                                            :: valid
      REAL(dp)                                           :: energy
      REAL(dp), ALLOCATABLE                              :: gradient(:), targets(:), tolerances(:)

      available = .FALSE.
      PRODUCT = 0.0_dp
      SELECT TYPE (context)
      TYPE IS (physical_constraints_context)
         IF (.NOT. qs_upf_inverse_scf_density_tangent_available(context%state%qs_env)) RETURN
         IF (.NOT. at_native_coordinate(context%state, coordinate) .OR. .NOT. context%state%operator_current) THEN
            ALLOCATE (gradient(SIZE(coordinate)))
            CALL qs_skala_scf_evaluate(context%state, coordinate, energy, gradient)
         END IF
         CALL qs_skala_scf_selected(context%state, points, domains, spins, descriptors, targets, tolerances, valid)
         IF (.NOT. valid) RETURN
         CPASSERT(ALL(SHAPE(points) == SHAPE(context%points)))
         CPASSERT(ALL(points == context%points) .AND. ALL(domains == context%domains))
         CPASSERT(ALL(spins == context%spins) .AND. ALL(descriptors == context%descriptors))
         CPASSERT(ALL(targets == context%targets))
         CALL qs_skala_scf_lagrangian_curvature(context%state, direction, multipliers, PRODUCT, available)
      CLASS DEFAULT
         CPABORT("Invalid native Lagrangian Hessian context.")
      END SELECT
   END SUBROUTINE evaluate_native_lagrangian_hessian

! **************************************************************************************************
!> \brief Invert current frozen-Hamiltonian occupation curvature in the fixed search coordinates.
!>        The original number projection remains in the Hessian product. Only the positive
!>        diagonal inverse is used here. Zero-response variables retain unit inverse.
!>        The caller freezes this preconditioner during each coupled Newton solve.
!> \param coordinate Accepted search coordinates.
!> \param context Borrowed physical state and native OT settings.
!> \param inverse_diagonal Positive initial inverse for BFGS-preconditioned directions.
!> \param available Whether native occupation preconditioning is enabled.
! **************************************************************************************************
   SUBROUTINE evaluate_native_preconditioner(coordinate, context, inverse_diagonal, available)
      REAL(dp), INTENT(IN)                               :: coordinate(:)

      CLASS(*), INTENT(INOUT) :: context
      REAL(dp), INTENT(OUT) :: inverse_diagonal(:)
      LOGICAL, INTENT(OUT) :: available

      COMPLEX(dp), ALLOCATABLE :: hamiltonian(:, :)
      REAL(dp), ALLOCATABLE :: diagonal(:), block_inverse(:), gradient(:)
      TYPE(dbcsr_p_type), POINTER :: matrices(:, :)
      INTEGER :: i, offset, local, spin, ik, first, last
      LOGICAL :: owner
      REAL(dp) :: energy, log_inverse, curvature_floor

      available = .FALSE.
      inverse_diagonal = 1.0_dp
      SELECT TYPE (context)
      TYPE IS (physical_constraints_context)
         IF (.NOT. ASSOCIATED(context%state%qs_env%scf_env)) RETURN
         IF (.NOT. ASSOCIATED(context%state%qs_env%scf_env%qs_ot_env)) RETURN
         IF (SIZE(context%state%qs_env%scf_env%qs_ot_env) == 0) RETURN
         IF (context%state%thermal) THEN
            IF (.NOT. context%state%qs_env%scf_env%qs_ot_env(1)%settings%occupation_preconditioner) RETURN
            curvature_floor = context%state%temperature
         ELSE
            IF (context%state%qs_env%scf_env%qs_ot_env(1)%settings%preconditioner_type == ot_precond_none) RETURN
            ! A positive gap regularizes only the search metric, not the zero-temperature energy.
            curvature_floor = context%state%qs_env%scf_env%qs_ot_env(1)%settings%energy_gap
            IF (curvature_floor <= 0.0_dp) RETURN
         END IF
         IF (.NOT. context%state%evaluated .OR. .NOT. context%state%operator_current .OR. &
             ANY(context%state%coordinate /= coordinate)) THEN
            ALLOCATE (gradient(SIZE(coordinate)))
            CALL qs_skala_scf_evaluate(context%state, coordinate, energy, gradient)
         END IF
         owner = context%state%para_env%mepos == 0
         IF (context%state%bloch) THEN
            owner = context%state%kpoints%para_env_kp%mepos == 0
         ELSE
            CALL get_ks_env(context%state%ks_env, matrix_ks_kp=matrices)
         END IF
         inverse_diagonal(:context%state%orbital_size) = 0.0_dp
         DO spin = 1, SIZE(context%state%orbitals, 2)
            DO local = 1, SIZE(context%state%orbitals, 1)
               IF (.NOT. ASSOCIATED(context%state%orbitals(local, spin)%native)) CYCLE
               ik = 1
               IF (context%state%bloch) ik = context%state%kpoints%kp_env(local)%kpoint_env%nkpoint
               IF (context%state%bloch) THEN
                  CALL native_bloch_operator(context%state, context%state%inverse%ks_k(:, :, :, spin), ik, hamiltonian)
               ELSE
                  CALL qs_upf_inverse_matrix_gather(context%state%inverse%frame%ao, &
                                                    matrices(spin, 1)%matrix, hamiltonian)
               END IF
               IF (.NOT. owner) CYCLE
               diagonal = REAL(SUM(CONJG(context%state%orbitals(local, spin)%coefficients)* &
                                   MATMUL(hamiltonian, context%state%orbitals(local, spin)%coefficients), DIM=1), dp)
               first = context%state%offsets(ik, spin) + 1
               last = first + context%state%dimensions(ik, spin) - 1
               ALLOCATE (block_inverse(last - first + 1))
               CALL skala_orbital_chart_response_inverse(context%state%orbitals(local, spin)%chart, &
                                                         diagonal, curvature_floor, block_inverse)
               inverse_diagonal(first:last) = block_inverse
               DEALLOCATE (block_inverse)
            END DO
         END DO
         CALL context%state%para_env%sum(inverse_diagonal(:context%state%orbital_size))
         IF (.NOT. context%state%thermal) THEN
            available = .TRUE.
            RETURN
         END IF
         offset = context%state%orbital_size
         DO i = 1, SIZE(context%state%band_response)
            IF (context%state%band_response(i) <= 0.0_dp) CYCLE
            log_inverse = LOG(context%state%temperature) - 2.0_dp*LOG(context%state%band_scales(i)) - &
                          LOG(context%state%band_response(i))
            inverse_diagonal(offset + i) = EXP(MAX(LOG(TINY(1.0_dp)), &
                                                   MIN(LOG(HUGE(1.0_dp)) - LOG(2.0_dp), log_inverse)))
         END DO
         available = .TRUE.
      END SELECT
   END SUBROUTINE evaluate_native_preconditioner

! **************************************************************************************************
!> \brief Keep numerical curvature probes within one physical branch.
!>        A threshold crossing requires a smaller probe, not a minimum classification.
!>        The ordinary minimization retains its complete boundary admission checks.
!> \param first Accepted feasible coordinates.
!> \param last Restored curvature probe.
!> \param curve Current independent/dependent coordinate partition.
!> \param context Native physical state and admission settings.
!> \param smooth Whether no new threshold is crossed.
!> \param interrupted Whether a boundary or external stop requires caller handling.
! **************************************************************************************************
   SUBROUTINE check_native_curvature(first, last, curve, context, smooth, interrupted)
      REAL(dp), INTENT(IN)                               :: first(:), last(:)
      TYPE(skala_constraint_chart_type), INTENT(IN)      :: curve

      CLASS(*), INTENT(INOUT) :: context
      LOGICAL, INTENT(OUT) :: smooth, interrupted

      LOGICAL :: crossed

      CALL check_native_segment(first, last, curve, context, interrupted, crossed, restore_operator=.FALSE., &
                                classify_candidates=.FALSE.)
      smooth = .NOT. crossed
      IF (crossed) THEN
         ! Reject this Hessian probe without turning it into an accepted orbital event.
         SELECT TYPE (context)
         TYPE IS (physical_constraints_context)
            context%crossings%encountered = .FALSE.
         END SELECT
         interrupted = .FALSE.
      END IF
   END SUBROUTINE check_native_curvature

! **************************************************************************************************
!> \brief Detect physical descriptor crossings on one native orbital segment.
!> \param first Accepted feasible coordinates.
!> \param last Proposed segment endpoint.
!> \param curve Current independent/dependent coordinate partition.
!> \param context Native physical state and admission settings.
!> \param interrupted Whether a crossing, rechart or external stop requires caller handling.
!> \param crossed Whether any physical threshold is crossed.
!> \param restore_operator Restore the accepted native operator after checking the segment.
!> \param classify_candidates Classify and admit detected roots when true.
! **************************************************************************************************
   SUBROUTINE check_native_segment(first, last, curve, context, interrupted, crossed, restore_operator, classify_candidates)
      REAL(dp), INTENT(IN)                               :: first(:), last(:)
      TYPE(skala_constraint_chart_type), INTENT(IN)      :: curve

      INTEGER                                            :: i, j
      INTEGER, ALLOCATABLE                               :: retained(:)
      LOGICAL, ALLOCATABLE                               :: keep(:)

      CLASS(*), INTENT(INOUT)                            :: context
      LOGICAL, INTENT(OUT)                               :: interrupted
      LOGICAL, INTENT(OUT), OPTIONAL                     :: crossed
      LOGICAL, INTENT(IN), OPTIONAL                      :: restore_operator, classify_candidates

      TYPE(qs_skala_scf_crossing_type), POINTER         :: crossings

      interrupted = .FALSE.
      IF (PRESENT(crossed)) crossed = .FALSE.
      SELECT TYPE (context)
      TYPE IS (physical_constraints_context)
         IF (context%max_coordinate_norm > 0.0_dp) THEN
            IF (NORM2(last) > context%max_coordinate_norm) THEN
               context%rechart = .TRUE.
               interrupted = .TRUE.
               RETURN
            END IF
         END IF
         IF (.NOT. ASSOCIATED(context%crossings)) RETURN
         crossings => context%crossings
         ! A classified candidate can be released without interrupting the minimizer.
         ! Its next trial must start with a fresh event rather than reallocate saved endpoints.
         crossings%encountered = .FALSE.
         crossings%resolved = .FALSE.
         IF (ALLOCATED(crossings%first)) DEALLOCATE (crossings%first)
         IF (ALLOCATED(crossings%last)) DEALLOCATE (crossings%last)
         IF (ALLOCATED(crossings%root)) DEALLOCATE (crossings%root)
         IF (ALLOCATED(crossings%selected)) DEALLOCATE (crossings%selected)
         IF (SIZE(context%density_floors) > 0) THEN
            CALL qs_skala_scf_crossings(context%state, first, last, crossings%points, crossings%domains, &
                                        crossings%spins, crossings%descriptors, crossings%values, &
                                        density_only=.NOT. context%kinetic_boundaries, &
                                        spin_density_targets=context%density_floors, candidate_targets=crossings%targets, &
                                        restore_operator=restore_operator)
         ELSE
            CALL qs_skala_scf_crossings(context%state, first, last, crossings%points, crossings%domains, &
                                        crossings%spins, crossings%descriptors, crossings%values, &
                                        density_only=.NOT. context%kinetic_boundaries, candidate_targets=crossings%targets, &
                                        restore_operator=restore_operator)
         END IF
         ! Existing active samples stay on their original constraint curve. Rounding-scale
         ! endpoint sign changes within that same window are not new physical events.
         ALLOCATE (keep(SIZE(crossings%spins)), SOURCE=.TRUE.)
         DO i = 1, SIZE(crossings%spins)
            DO j = 1, SIZE(context%spins)
               IF (crossings%domains(i) /= context%domains(j)) CYCLE
               IF (crossings%spins(i) /= context%spins(j)) CYCLE
               IF (crossings%descriptors(i) /= context%descriptors(j)) CYCLE
               IF (crossings%targets(i) /= context%targets(j)) CYCLE
               IF (ALL(crossings%points(:, i) == context%points(:, j))) keep(i) = .FALSE.
            END DO
         END DO
         retained = PACK([(i, i=1, SIZE(keep))], keep)
         crossings%points = crossings%points(:, retained)
         crossings%domains = crossings%domains(retained)
         crossings%spins = crossings%spins(retained)
         crossings%descriptors = crossings%descriptors(retained)
         crossings%targets = crossings%targets(retained)
         crossings%values = crossings%values(:, retained)
         crossings%curve = curve
         interrupted = SIZE(crossings%spins) > 0
         IF (.NOT. interrupted) RETURN
         IF (PRESENT(crossed)) crossed = .TRUE.
         ALLOCATE (crossings%first(SIZE(first)), crossings%last(SIZE(last)))
         crossings%first(:) = first
         crossings%last(:) = last
         crossings%encountered = .TRUE.
         IF (context%admit_crossings) THEN
            IF (PRESENT(classify_candidates)) THEN
               IF (.NOT. classify_candidates) RETURN
            END IF
            CALL classify_native_segment(first, last, context)
            interrupted = crossings%encountered .OR. context%external_stopped
         END IF
      CLASS DEFAULT
         CPABORT("Invalid native orbital crossing context.")
      END SELECT
   END SUBROUTINE check_native_segment

! **************************************************************************************************
!> \brief Resolve and classify further continuous-floor candidates before interrupting minimization.
!>        Keep the original selected samples and physical windows. Reject an uphill boundary or
!>        one permitting normal descent, then inspect the remaining endpoint candidates. An
!>        unresolved root leaves the original line search to test its complete energy and derivative.
!>        A failed normal classification remains an event for the caller.
!> \param first Last accepted feasible coordinates.
!> \param last Proposed feasible coordinates.
!> \param context Native state, original constraint curve and optional crossing output.
! **************************************************************************************************
   SUBROUTINE classify_native_segment(first, last, context)
      REAL(dp), INTENT(IN)                               :: first(:), last(:)
      TYPE(physical_constraints_context), INTENT(INOUT)  :: context

      INTEGER                                            :: i, old
      INTEGER, ALLOCATABLE                               :: old_descriptors(:), old_domains(:), &
                                                            old_points(:, :), old_spins(:), &
                                                            pending(:), reached(:), sides(:)
      LOGICAL                                            :: found, stopped, valid
      LOGICAL, ALLOCATABLE                               :: admissible(:), selected(:)
      REAL(dp)                                           :: energy, first_energy, kinetic_probe, &
                                                            kinetic_tolerance, probe, roundoff, &
                                                            tolerance
      REAL(dp), ALLOCATABLE                              :: gradient(:), root(:), targets(:), &
                                                            values(:), weights(:), windows(:)
      TYPE(qs_skala_scf_crossing_type), POINTER          :: event

      event => context%crossings
      old = skala_density_constraints_count(context%state%constraints)
      tolerance = 1.0e-17_dp
      kinetic_tolerance = 1.0e-17_dp
      IF (old > 0) THEN
         CALL skala_density_constraints_get(context%state%constraints, old_points, old_spins, weights, probe, tolerance, &
                                            descriptors=old_descriptors, domains=old_domains, &
                                            targets=targets, kinetic_probe=kinetic_probe, kinetic_tolerance=kinetic_tolerance)
      END IF
      ALLOCATE (pending(SIZE(event%spins)), root(SIZE(first)), gradient(SIZE(first)))
      pending(:) = [(i, i=1, SIZE(pending))]
      CALL qs_skala_scf_evaluate(context%state, first, first_energy)
      DO WHILE (SIZE(pending) > 0)
         ALLOCATE (selected(SIZE(pending)), values(SIZE(pending)), windows(SIZE(pending)))
         windows(:) = MERGE(tolerance, kinetic_tolerance, event%descriptors(pending) == skala_density_descriptor)
         CALL qs_skala_scf_boundary_root(context%state, first, last, event%points(:, pending), &
                                         event%domains(pending), event%spins(pending), event%descriptors(pending), &
                                         event%targets(pending), windows, root, values, selected, found, &
                                         check_external=context%check_external, stopped=stopped, curve=event%curve)
         IF (stopped) THEN
            context%external_stopped = .TRUE.
            EXIT
         END IF
         IF (.NOT. found) THEN
            ! An unresolved target supplies no verified constraint. Continue the original
            ! line search without adding it or relaxing physical convergence. Curvature
            ! probes still detect the crossing and must remain on a single branch.
            event%encountered = .FALSE.
            EXIT
         END IF
         ! An uphill boundary cannot be accepted, regardless of its limiting derivatives.
         ! Reject it before rebuilding every normal limit of the complete XC operator.
         CPASSERT(context%state%evaluated .AND. context%state%operator_current)
         CPASSERT(ALL(context%state%coordinate == root))
         energy = context%state%energy%total
         roundoff = 64.0_dp*EPSILON(0.0_dp)*MAX(ABS(first_energy), ABS(energy), TINY(0.0_dp))
         IF (energy > first_energy + roundoff) THEN
            pending = PACK(pending,.NOT. selected)
            DEALLOCATE (selected, values, windows)
            CYCLE
         END IF
         reached = PACK(pending, selected)
         ALLOCATE (admissible(SIZE(reached)), sides(SIZE(reached)))
         CALL qs_skala_scf_classify(context%state, root, event%points(:, reached), event%domains(reached), &
                                    event%spins(reached), event%descriptors(reached), event%targets(reached), &
                                    admissible, sides, valid, restore_operator=.FALSE.)
         IF (.NOT. valid) EXIT
         IF (ANY(admissible)) THEN
            ! Restore the original operator only when an admissible root can be returned.
            CALL evaluate_limiting_state(context%state, root, energy, gradient, .TRUE.)
            roundoff = 64.0_dp*EPSILON(0.0_dp)*MAX(ABS(first_energy), ABS(energy), TINY(0.0_dp))
            IF (energy <= first_energy + roundoff) THEN
               ALLOCATE (event%selected(SIZE(event%spins)), SOURCE=.FALSE.)
               event%selected(PACK(reached, admissible)) = .TRUE.
               ALLOCATE (event%root, SOURCE=root)
               event%resolved = .TRUE.
               EXIT
            END IF
         END IF
         pending = PACK(pending,.NOT. selected)
         DEALLOCATE (selected, values, windows, admissible, sides)
      END DO
      IF (SIZE(pending) == 0) event%encountered = .FALSE.
      IF (.NOT. event%resolved) THEN
         CALL qs_skala_scf_evaluate(context%state, first, energy)
      END IF
   END SUBROUTINE classify_native_segment

! **************************************************************************************************
!> \brief Evaluate one physical descriptor along an orbital trial segment.
!> \param step Position along the original orbital trial segment.
!> \param value Physical descriptor value at the trial.
!> \param context Borrowed native orbital and constraint state.
! **************************************************************************************************
   SUBROUTINE evaluate_physical_input(step, value, context)
      REAL(dp), INTENT(IN)                               :: step
      REAL(dp), INTENT(OUT)                              :: value

      CLASS(*), INTENT(INOUT)                           :: context

      REAL(dp), ALLOCATABLE                            :: coordinate(:)

      SELECT TYPE (context)
      TYPE IS (physical_root_context)
         ALLOCATE (coordinate(SIZE(context%first)))
         coordinate(:) = context%first + step*context%direction
         CALL sample_native_constraints(context%state, coordinate)
         CALL qs_skala_scf_value(context%state, context%domain, context%point, context%spin, context%descriptor, value)
      CLASS DEFAULT
         CPABORT("Invalid native Skala root context.")
      END SELECT
   END SUBROUTINE evaluate_physical_input

! **************************************************************************************************
!> \brief Append dimensionless auxiliary-energy coordinates to the retained native orbital charts.
!>        The original finite-temperature setup determines maximum occupation, particle counts
!>        and whether the two spin channels share a chemical potential. Global band ordering
!>        permits ragged retained spaces and includes each positive k-point weight once.
!> \param state Native orbital state receiving the thermal particle-number groups.
!> \param offset Running coordinate count, increased by the thermal band count.
! **************************************************************************************************
   SUBROUTINE initialize_thermal_state(state, offset)
      TYPE(qs_skala_scf_type), INTENT(INOUT)             :: state
      INTEGER, INTENT(INOUT)                             :: offset

      INTEGER                                            :: count, first, groups, ik, last, local, &
                                                            nk, ns, spin
      LOGICAL                                            :: owner, separate
      REAL(dp), ALLOCATABLE                              :: targets(:)
      TYPE(dft_control_type), POINTER                    :: dft
      TYPE(mo_set_type), POINTER                         :: mos(:)
      TYPE(scf_control_type), POINTER                    :: control

      CALL get_qs_env(state%qs_env, scf_control=control, dft_control=dft, mos=mos)
      CPASSERT(control%smear%do_smear .AND. control%smear%method == smear_fermi_dirac)
      CPASSERT(control%smear%electronic_temperature > 0.0_dp)
      CPASSERT(.NOT. dft%roks .AND. .NOT. dft%qs_control%cdft)
      CPASSERT(.NOT. control%gce%do_gce)
      ns = SIZE(state%orbitals, 2)
      nk = SIZE(state%dimensions, 1)
      state%temperature = control%smear%electronic_temperature
      state%maxocc = 2.0_dp/REAL(ns, dp)
      separate = ns == 1 .OR. control%smear%fixed_mag_mom > 0.0_dp .OR. dft%restricted
      groups = MERGE(ns, 1, separate)
      ALLOCATE (state%band_dimensions(nk, ns), state%band_offsets(nk, ns), targets(ns), &
                state%thermal_targets(groups), state%thermal_mu(groups), state%thermal_spin_entropy(ns))
      state%band_dimensions = 0
      owner = state%para_env%mepos == 0
      IF (state%bloch) owner = state%kpoints%para_env_kp%mepos == 0
      DO spin = 1, ns
         targets(spin) = REAL(mos(spin)%nelectron, dp)
         IF (state%bloch) targets(spin) = REAL(state%kpoints%kp_env(1)%kpoint_env%mos(spin)%nelectron, dp)
         DO local = 1, SIZE(state%orbitals, 1)
            IF (.NOT. ASSOCIATED(state%orbitals(local, spin)%native)) CYCLE
            ik = 1
            IF (state%bloch) ik = state%kpoints%kp_env(local)%kpoint_env%nkpoint
            IF (owner) state%band_dimensions(ik, spin) = state%orbitals(local, spin)%native%nmo
         END DO
      END DO
      CPASSERT(ALL(targets >= 0.0_dp))
      IF (separate) THEN
         state%thermal_targets = targets
      ELSE
         state%thermal_targets(1) = SUM(targets)
      END IF
      CALL state%para_env%sum(state%band_dimensions)
      count = 0
      DO spin = 1, ns
         DO ik = 1, nk
            state%band_offsets(ik, spin) = count
            count = count + state%band_dimensions(ik, spin)
         END DO
      END DO
      ALLOCATE (state%band_reference(count), state%band_energies(count), state%band_weights(count), &
                state%band_occupations(count), state%band_response(count), state%band_groups(count), &
                state%band_spins(count), state%band_scales(count))
      state%band_scales = SQRT(state%temperature*thermal_energy_unit)
      state%band_reference = 0.0_dp
      state%band_weights = 0.0_dp
      DO spin = 1, ns
         DO ik = 1, nk
            first = state%band_offsets(ik, spin) + 1
            last = first + state%band_dimensions(ik, spin) - 1
            state%band_groups(first:last) = MERGE(spin, 1, separate)
            state%band_spins(first:last) = spin
         END DO
         DO local = 1, SIZE(state%orbitals, 1)
            IF (.NOT. ASSOCIATED(state%orbitals(local, spin)%native)) CYCLE
            ik = 1
            IF (state%bloch) ik = state%kpoints%kp_env(local)%kpoint_env%nkpoint
            first = state%band_offsets(ik, spin) + 1
            last = first + state%band_dimensions(ik, spin) - 1
            IF (owner) THEN
               state%band_reference(first:last) = state%orbitals(local, spin)%native%eigenvalues(:last - first + 1)
               state%band_weights(first:last) = state%orbitals(local, spin)%weight
            END IF
         END DO
      END DO
      CALL state%para_env%sum(state%band_reference)
      CALL state%para_env%sum(state%band_weights)
      CPASSERT(ALL(ieee_is_finite(state%band_reference)))
      CPASSERT(ALL(state%band_weights > 0.0_dp))
      offset = offset + count
   END SUBROUTINE initialize_thermal_state

! **************************************************************************************************
!> \brief Solve the original weighted Fermi constraint in every particle-number group.
!>        Frozen positive scale factors precondition the auxiliary-energy coordinates.
!>        The physical Fermi problem and Fisher covectors are independent of these factors.
!>        Smearkp accepts one band per weighted column, preserving ragged block counts.
!> \param state Native state with frozen references and coordinate scales.
!> \param coordinate Dimensionless auxiliary-energy coordinates.
! **************************************************************************************************
   SUBROUTINE update_thermal_state(state, coordinate)
      TYPE(qs_skala_scf_type), INTENT(INOUT)             :: state
      REAL(dp), INTENT(IN)                               :: coordinate(:)

      INTEGER                                            :: group, i, spin
      INTEGER, ALLOCATABLE                               :: indices(:)
      REAL(dp)                                           :: entropy, f, term
      REAL(dp), ALLOCATABLE                              :: energies(:, :), occupations(:, :)

      CPASSERT(SIZE(coordinate) == SIZE(state%band_reference))
      state%band_energies = state%band_reference + state%band_scales*coordinate
      state%entropy = 0.0_dp
      DO group = 1, SIZE(state%thermal_targets)
         indices = PACK([(i, i=1, SIZE(coordinate))], state%band_groups == group)
         CPASSERT(SIZE(indices) > 0)
         ALLOCATE (energies(1, SIZE(indices)), occupations(1, SIZE(indices)))
         energies(1, :) = state%band_energies(indices)
         CALL Smearkp(occupations, state%thermal_mu(group), entropy, energies, state%thermal_targets(group), &
                      state%band_weights(indices), state%temperature, state%maxocc, smear_fermi_dirac)
         state%band_occupations(indices) = occupations(1, :)
         state%entropy = state%entropy + entropy
         DEALLOCATE (energies, occupations)
      END DO
      state%band_response = state%band_weights*state%band_occupations*(1.0_dp - state%band_occupations/state%maxocc)
      state%thermal_spin_entropy = 0.0_dp
      DO i = 1, SIZE(coordinate)
         spin = state%band_spins(i)
         f = state%band_occupations(i)/state%maxocc
         term = 0.0_dp
         IF (f > 0.0_dp) term = term + f*LOG(f)
         IF (f < 1.0_dp) term = term + (1.0_dp - f)*LOG(1.0_dp - f)
         state%thermal_spin_entropy(spin) = state%thermal_spin_entropy(spin) + &
                                            state%temperature*state%band_weights(i)*state%maxocc*term
      END DO
   END SUBROUTINE update_thermal_state

! **************************************************************************************************
!> \brief Freeze the reference Fermi-response normalization of the auxiliary-energy variables.
!>        With p_i=w_i f_i(1-f_i/maxocc), s_i=sqrt(kBT times one hartree divided by p_i)
!>        makes the unregularized frozen-Hamiltonian Hessian a group projector. A group-mean
!>        response floor limits scale ratios to sqrt(group size) and keeps empty bands active.
!> \param state Reference thermal occupations and fixed particle-number groups.
! **************************************************************************************************
   SUBROUTINE initialize_thermal_scales(state)
      TYPE(qs_skala_scf_type), INTENT(INOUT)             :: state

      INTEGER                                            :: group, i
      INTEGER, ALLOCATABLE                               :: indices(:)
      REAL(dp)                                           :: floor, peak

      DO group = 1, SIZE(state%thermal_targets)
         indices = PACK([(i, i=1, SIZE(state%band_response))], state%band_groups == group)
         peak = MAXVAL(state%band_response(indices))
         IF (peak <= 0.0_dp) CYCLE
         ! A group-mean floor limits reference conditioning without freezing empty bands.
         floor = MAX(TINY(1.0_dp), SUM(state%band_response(indices))/REAL(SIZE(indices), dp))
         state%band_scales(indices) = SQRT(state%temperature*thermal_energy_unit/ &
                                           MAX(state%band_response(indices), floor))
      END DO
   END SUBROUTINE initialize_thermal_scales

! **************************************************************************************************
!> \brief Store unweighted diagonal expectations for one owned retained orbital block.
!> \param state Current normalized native orbitals and global band layout.
!> \param OPERATOR AO operator acting on the retained orbitals.
!> \param local Owned k-point block.
!> \param spin Spin component of the retained block.
!> \param values Global unweighted band expectations, updated for this owned block.
! **************************************************************************************************
   SUBROUTINE thermal_expectations(state, OPERATOR, local, spin, values)
      TYPE(qs_skala_scf_type), INTENT(IN)                :: state
      COMPLEX(dp), INTENT(IN)                            :: operator(:, :)
      INTEGER, INTENT(IN)                                :: local, spin
      REAL(dp), INTENT(INOUT)                            :: values(:)

      COMPLEX(dp), ALLOCATABLE                           :: product(:, :)
      INTEGER                                            :: first, i, ik, last

      ik = 1
      IF (state%bloch) ik = state%kpoints%kp_env(local)%kpoint_env%nkpoint
      first = state%band_offsets(ik, spin) + 1
      last = first + state%band_dimensions(ik, spin) - 1
      product = MATMUL(OPERATOR, state%orbitals(local, spin)%coefficients(:, :last - first + 1))
      DO i = 1, last - first + 1
         values(first + i - 1) = REAL(SUM(CONJG(state%orbitals(local, spin)%coefficients(:, i))*PRODUCT(:, i)), dp)
      END DO
   END SUBROUTINE thermal_expectations

! **************************************************************************************************
!> \brief Include the complete fixed-N response of thermal occupations in a covector.
!>        Subtract auxiliary energies only for the free-energy gradient, which includes entropy.
!>        Descriptor normals have no entropy contribution. The same global group projection
!>        enforces the original particle/spin counts in both cases.
!> \param state Current thermal occupations, response factors and particle-number groups.
!> \param expectations Unweighted band expectations of the differentiated operator.
!> \param entropy Include the thermal entropy contribution in the free-energy covector.
!> \param derivative Covector in the dimensionless auxiliary-energy coordinates.
!> \param physical Optional covector normalized by the physical Fermi response.
! **************************************************************************************************
   SUBROUTINE thermal_covector(state, expectations, entropy, derivative, physical)
      TYPE(qs_skala_scf_type), INTENT(IN)                :: state
      REAL(dp), INTENT(IN)                               :: expectations(:)
      LOGICAL, INTENT(IN)                                :: entropy
      REAL(dp), INTENT(OUT)                              :: derivative(:)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: physical(:)

      INTEGER                                            :: group, i
      INTEGER, ALLOCATABLE                               :: indices(:)
      REAL(dp)                                           :: mean, weight
      REAL(dp), ALLOCATABLE                              :: residual(:)

      CPASSERT(SIZE(expectations) == SIZE(state%band_reference))
      CPASSERT(SIZE(derivative) == SIZE(expectations))
      residual = expectations
      IF (entropy) residual = residual - state%band_energies
      derivative = 0.0_dp
      IF (PRESENT(physical)) physical = 0.0_dp
      DO group = 1, SIZE(state%thermal_targets)
         indices = PACK([(i, i=1, SIZE(expectations))], state%band_groups == group)
         weight = SUM(state%band_response(indices))
         IF (weight <= TINY(weight)) CYCLE
         mean = SUM(state%band_response(indices)*residual(indices))/weight
         derivative(indices) = -state%band_response(indices)*(residual(indices) - mean)
         IF (PRESENT(physical)) physical(indices) = -SQRT(state%band_response(indices))*(residual(indices) - mean)
      END DO
      derivative = state%band_scales/state%temperature*derivative
   END SUBROUTINE thermal_covector

! **************************************************************************************************
!> \brief Check the complete physical Fermi density in the original charge metric.
!>        The full normalized Gaussian space supplies the missing thermal complement. This
!>        operation never fills or grows the retained MOs. Significant omitted charge rejects
!>        acceptance independently of the optimizer gradient.
!> \param state Canonical native orbitals and matching physical Hamiltonian.
!> \param tolerance Original physical gradient bound in energy units.
!> \param fermi_error Weighted Frobenius distance to the complete Fermi density.
!> \param omitted_number Thermal charge outside the retained native MO spaces.
!> \param chemical_potentials Physical Fermi levels of the original particle-number groups.
!> \param valid Whether both physical density and omitted-space bounds are met.
! **************************************************************************************************
   SUBROUTINE thermal_density_stationarity(state, tolerance, fermi_error, omitted_number, chemical_potentials, valid)
      TYPE(qs_skala_scf_type), INTENT(INOUT)             :: state
      REAL(dp), INTENT(IN)                               :: tolerance
      REAL(dp), INTENT(OUT)                              :: fermi_error, omitted_number
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: chemical_potentials(:)
      LOGICAL, INTENT(OUT)                               :: valid

      INTEGER                                            :: n, nk, ns
      REAL(dp)                                           :: error_bound
      TYPE(qs_upf_inverse_type), POINTER                 :: frame
      TYPE(scf_control_type), POINTER                    :: control

      CALL qs_skala_scf_fermi_target(state, fermi_error, omitted_number, chemical_potentials, valid)
      IF (.NOT. valid) RETURN
      frame => state%inverse%frame%ao
      IF (state%bloch) frame => state%inverse%bloch%frames(1)%ao
      n = SIZE(frame%coefficients, 2)
      nk = SIZE(state%dimensions, 1)
      ns = SIZE(state%orbitals, 2)
      CALL get_qs_env(state%qs_env, scf_control=control)
      error_bound = MAX(tolerance/state%temperature, 64.0_dp*EPSILON(1.0_dp)*SQRT(REAL(n*nk*ns, dp)))
      valid = fermi_error <= error_bound .AND. omitted_number <= control%smear%eps_fermi_dirac
   END SUBROUTINE thermal_density_stationarity

! **************************************************************************************************
!> \brief Assemble the complete frozen-Hamiltonian Fermi density in the original charge metric.
!>        Use the same spin/particle-number groups and positive k-point weights as native OT.
!>        The optional density contains unweighted AO matrices. No MOs, fields, energies or
!>        convergence criteria are changed. A finite nonlinear trial must still be checked.
!> \param state Evaluated thermal state with the current native Hamiltonian.
!> \param fermi_error Weighted Frobenius distance in the complete normalized Gaussian space.
!> \param omitted_number Fermi charge outside the retained MO spaces.
!> \param chemical_potentials Fermi levels of the original particle-number groups.
!> \param valid Whether every full-space diagonalization and filling succeeded.
!> \param density Optional complete Fermi density, indexed by AO, AO, k point and spin.
!> \param model_delta Optional frozen-H free-energy change, including the original entropy.
! **************************************************************************************************
   SUBROUTINE qs_skala_scf_fermi_target(state, fermi_error, omitted_number, chemical_potentials, valid, &
                                        density, model_delta)
      TYPE(qs_skala_scf_type), INTENT(INOUT)             :: state
      REAL(dp), INTENT(OUT)                              :: fermi_error, omitted_number
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: chemical_potentials(:)
      LOGICAL, INTENT(OUT)                               :: valid
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT), OPTIONAL    :: density(:, :, :, :)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: model_delta

      COMPLEX(dp), ALLOCATABLE                           :: hamiltonian(:, :), projected(:, :), &
                                                            residual(:, :), vectors(:, :, :, :), &
                                                            weighted(:, :), work(:)
      INTEGER                                            :: count, failed, first, group, i, ik, &
                                                            info, last, local, n, nk, nl, ns, &
                                                            retained, spin
      INTEGER, ALLOCATABLE                               :: dimensions(:, :), groups(:), indices(:), &
                                                            offsets(:, :)
      LOGICAL                                            :: owner
      REAL(dp)                                           :: correction, reference_energy, &
                                                            target_entropy, weight
      REAL(dp), ALLOCATABLE                              :: energies(:, :), filled(:, :), &
                                                            mesh_weights(:), occupations(:), &
                                                            real_work(:), spectrum(:, :, :)
      TYPE(dbcsr_p_type), POINTER                        :: matrices(:, :)

      CPASSERT(state%thermal .AND. state%evaluated .AND. state%operator_current)
      IF (PRESENT(model_delta)) model_delta = HUGE(0.0_dp)
      valid = .FALSE.
      fermi_error = 0.0_dp
      omitted_number = 0.0_dp
      nk = SIZE(state%dimensions, 1)
      ns = SIZE(state%orbitals, 2)
      nl = SIZE(state%orbitals, 1)
      n = 0
      DO spin = 1, ns
         DO local = 1, nl
            IF (ASSOCIATED(state%orbitals(local, spin)%native)) n = state%orbitals(local, spin)%native%nao
         END DO
      END DO
      CALL state%para_env%max(n)
      CPASSERT(n > 0)
      ALLOCATE (spectrum(n, nk, ns), vectors(n, n, nl, ns), projected(n, n), &
                work(MAX(1, 2*n - 1)), real_work(MAX(1, 3*n - 2)), dimensions(nk, ns), offsets(nk, ns))
      spectrum = 0.0_dp
      dimensions = 0
      owner = state%para_env%mepos == 0
      IF (state%bloch) THEN
         owner = state%kpoints%para_env_kp%mepos == 0
      ELSE
         CALL get_ks_env(state%ks_env, matrix_ks_kp=matrices)
      END IF
      failed = 0
      reference_energy = 0.0_dp
      DO spin = 1, ns
         DO local = 1, nl
            IF (.NOT. ASSOCIATED(state%orbitals(local, spin)%native)) CYCLE
            ik = 1
            IF (state%bloch) THEN
               ik = state%kpoints%kp_env(local)%kpoint_env%nkpoint
               CALL native_bloch_operator(state, state%inverse%ks_k(:, :, :, spin), ik, hamiltonian)
            ELSE
               CALL qs_upf_inverse_matrix_gather(state%inverse%frame%ao, matrices(spin, 1)%matrix, hamiltonian)
            END IF
            IF (.NOT. owner) CYCLE
            projected = MATMUL(CONJG(TRANSPOSE(state%orbitals(local, spin)%coefficients)), &
                               MATMUL(hamiltonian, state%orbitals(local, spin)%coefficients))
            IF (PRESENT(model_delta)) THEN
               DO i = 1, n
                  reference_energy = reference_energy + state%orbitals(local, spin)%weight* &
                                     state%orbitals(local, spin)%occupations(i)*REAL(projected(i, i), dp)
               END DO
            END IF
            vectors(:, :, local, spin) = 0.5_dp*(projected + CONJG(TRANSPOSE(projected)))
            CALL zheev('V', 'U', n, vectors(:, :, local, spin), n, spectrum(:, ik, spin), &
                       work, SIZE(work), real_work, info)
            IF (info /= 0) failed = 1
            dimensions(ik, spin) = n
         END DO
      END DO
      CALL state%para_env%max(failed)
      IF (failed /= 0) RETURN
      IF (PRESENT(model_delta)) CALL state%para_env%sum(reference_energy)
      CALL state%para_env%sum(spectrum)
      CALL state%para_env%sum(dimensions)
      count = 0
      DO spin = 1, ns
         DO ik = 1, nk
            offsets(ik, spin) = count
            count = count + dimensions(ik, spin)
         END DO
      END DO
      ALLOCATE (energies(1, count), occupations(count), mesh_weights(count), groups(count), &
                chemical_potentials(SIZE(state%thermal_targets)))
      DO spin = 1, ns
         DO ik = 1, nk
            IF (dimensions(ik, spin) == 0) CYCLE
            first = offsets(ik, spin) + 1
            last = first + n - 1
            energies(1, first:last) = spectrum(:, ik, spin)
            group = state%band_groups(state%band_offsets(ik, spin) + 1)
            weight = state%band_weights(state%band_offsets(ik, spin) + 1)
            groups(first:last) = group
            mesh_weights(first:last) = weight
         END DO
      END DO
      target_entropy = 0.0_dp
      DO group = 1, SIZE(state%thermal_targets)
         indices = PACK([(i, i=1, count)], groups == group)
         ALLOCATE (filled(1, SIZE(indices)))
         CALL Smearkp(filled, chemical_potentials(group), correction, energies(:, indices), state%thermal_targets(group), &
                      mesh_weights(indices), state%temperature, state%maxocc, smear_fermi_dirac)
         occupations(indices) = filled(1, :)
         target_entropy = target_entropy + correction
         DEALLOCATE (filled)
      END DO
      ALLOCATE (weighted(n, n), residual(n, n))
      IF (PRESENT(density)) THEN
         ALLOCATE (density(n, n, nk, ns), SOURCE=CMPLX(0.0_dp, 0.0_dp, dp))
      END IF
      DO spin = 1, ns
         DO local = 1, nl
            IF (.NOT. ASSOCIATED(state%orbitals(local, spin)%native)) CYCLE
            IF (.NOT. owner) CYCLE
            ik = 1
            IF (state%bloch) ik = state%kpoints%kp_env(local)%kpoint_env%nkpoint
            first = offsets(ik, spin) + 1
            DO i = 1, n
               weighted(:, i) = occupations(first + i - 1)*vectors(:, i, local, spin)
            END DO
            residual = MATMUL(weighted, CONJG(TRANSPOSE(vectors(:, :, local, spin))))
            IF (PRESENT(density)) THEN
               density(:, :, ik, spin) = MATMUL(state%orbitals(local, spin)%coefficients, &
                                                MATMUL(residual, CONJG(TRANSPOSE( &
                                                                       state%orbitals(local, spin)%coefficients))))
            END IF
            weight = state%orbitals(local, spin)%weight
            retained = state%orbitals(local, spin)%native%nmo
            DO i = retained + 1, n
               omitted_number = omitted_number + weight*MAX(0.0_dp, REAL(residual(i, i), dp))
            END DO
            DO i = 1, n
               residual(i, i) = residual(i, i) - state%orbitals(local, spin)%occupations(i)
            END DO
            fermi_error = fermi_error + weight*SUM(ABS(residual)**2)
         END DO
      END DO
      CALL state%para_env%sum(fermi_error)
      CALL state%para_env%sum(omitted_number)
      fermi_error = SQRT(MAX(0.0_dp, fermi_error))
      IF (PRESENT(density)) CALL state%para_env%sum(density)
      IF (PRESENT(model_delta)) THEN
         ! Frozen-H free-energy change. The nonlinear XC energy is evaluated by the caller.
         model_delta = SUM(mesh_weights*occupations*energies(1, :)) + target_entropy - &
                       reference_energy - state%entropy
      END IF
      valid = .TRUE.
   END SUBROUTINE qs_skala_scf_fermi_target

END MODULE qs_skala_scf
