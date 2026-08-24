!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief orbital transformations
!> \par History
!>      None
!> \author Joost VandeVondele (09.2002)
! **************************************************************************************************
MODULE qs_ot_minimizer

   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_add, dbcsr_copy, dbcsr_create, dbcsr_get_info, dbcsr_multiply, dbcsr_p_type, &
        dbcsr_release, dbcsr_scale, dbcsr_set, dbcsr_type
   USE cp_dbcsr_contrib,                ONLY: dbcsr_dot,&
                                              dbcsr_get_diag,&
                                              dbcsr_init_random,&
                                              dbcsr_scale_by_vector
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: debug_print_level
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE mathlib,                         ONLY: diamat_all
   USE message_passing,                 ONLY: mp_para_env_type
   USE preconditioner,                  ONLY: apply_preconditioner
   USE qs_ot,                           ONLY: qs_ot_get_derivative,&
                                              qs_ot_get_derivative_ref,&
                                              qs_ot_get_derivative_ref_complex,&
                                              qs_ot_transport_ref_rotation_pair
   USE qs_ot_types,                     ONLY: qs_ot_occupation_preconditioner_active,&
                                              qs_ot_reset_accepted_model_calibration,&
                                              qs_ot_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   PUBLIC  :: ot_kpoint_cg_residual_growth_guard, &
              ot_kpoint_cg_metric_history_reset, &
              ot_kpoint_cg_tail_metric_tracked, &
              ot_kpoint_diis_accepted_response_step, &
              ot_kpoint_diis_rejected_step, &
              ot_kpoint_diis_response_endpoint_ratio, &
              ot_kpoint_diis_response_residual_blend, &
              ot_kpoint_diis_rotation_response_probe_ready, &
              ot_kpoint_diis_rotation_response_ready, &
              ot_kpoint_diis_mermin_step_rejected, &
              ot_kpoint_diis_residual_step_rejected, &
              ot_kpoint_fermi_window_energy_coordinate, &
              ot_kpoint_fermi_window_direction_scale, &
              ot_kpoint_mermin_closure_response_quality, &
              ot_kpoint_mermin_full_direction_model, &
              ot_kpoint_mermin_current_guard_delta, &
              ot_kpoint_mermin_derivative_refinement, &
              ot_kpoint_mermin_energy_curvature_resolved, &
              ot_kpoint_mermin_endpoint_status, &
              ot_kpoint_mermin_secant_consistent, &
              ot_kpoint_multiplicity_direction_scale, &
              ot_kpoint_physical_residual_split, &
              ot_kpoint_mermin_candidate_restart_step, &
              ot_kpoint_mermin_recovery_active_tail, &
              ot_kpoint_mermin_recovery_deferred, &
              ot_kpoint_mermin_projected_reference_ready, &
              ot_kpoint_mermin_reference_guard_delta, &
              ot_kpoint_mermin_search_stagnant, &
              ot_kpoint_mermin_restart_step, &
              ot_kpoint_mermin_shared_spin_group, &
              ot_kpoint_mermin_tail_energy_stagnant, &
              ot_kpoint_physical_exit_delta, &
              ot_kpoint_residual_aware_rayleigh_model, &
              ot_kpoint_rayleigh_probe_position, &
              ot_kpoint_rayleigh_model_factor, &
              ot_kpoint_rayleigh_response_scale, &
              ot_kpoint_secant_geometry, &
              sync_kpoint_minimizer_control, &
              ot_mini

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_ot_minimizer'
   INTEGER, PARAMETER, PRIVATE :: ot_cg_kp_physical_control_min_samples = 8, &
                                  ot_cg_kp_coupled_closure_flavor = 42, &
                                  ot_cg_kp_energy_closure_flavor = 45
   INTEGER, PARAMETER, PRIVATE :: ot_cg_kp_mermin_closure_max_backtracks = 8, &
                                  ot_diis_kp_response_probe_min_accepted = 4
   REAL(KIND=dp), PARAMETER, PRIVATE :: ot_cg_kp_history_spike_ratio = 1.0E4_dp, &
                                        ot_cg_kp_history_jump_ratio = 1.0E2_dp, &
                                        ot_cg_kp_history_jump_beta = 1.0E1_dp, &
                                        ot_cg_kp_projector_tail_beta = 1.0_dp, &
                                        ot_cg_kp_projector_tail_norm = 5.0E-9_dp, &
                                        ot_cg_kp_projector_tail_ratio = 1.2_dp, &
                                        ot_cg_kp_history_spike_old_norm = 1.0E-8_dp, &
                                        ot_cg_kp_history_spike_new_norm = 1.0E-6_dp, &
                                        ot_cg_kp_history_jump_new_norm = 1.0E-3_dp, &
                                        ot_cg_kp_projector_jump_ratio = 8.0E1_dp, &
                                        ot_cg_kp_projector_jump_new_norm = 1.0E-4_dp, &
                                        ot_cg_kp_residual_tail_jump_ratio = 5.0E1_dp, &
                                        ot_cg_kp_residual_tail_jump_beta = 1.0E1_dp, &
                                        ot_cg_kp_residual_tail_jump_norm = 1.0E-10_dp, &
                                        ot_cg_kp_residual_tail_min_delta = 3.0E-6_dp, &
                                        ot_cg_kp_residual_tail_max_delta = 5.0E-6_dp, &
                                        ot_cg_kp_residual_tail_min_eps_factor = 3.0_dp, &
                                        ot_cg_kp_mermin_residual_growth_ratio = 1.5_dp, &
                                        ot_cg_kp_mermin_catastrophic_growth_ratio = 4.0_dp, &
                                        ot_cg_kp_mermin_residual_growth_max_delta = 1.2E-5_dp, &
                                        ot_cg_kp_mermin_recovery_tail_min_delta = 1.5E-6_dp, &
                                        ot_cg_kp_mermin_recovery_tail_max_delta = 5.0E-6_dp, &
                                        ot_cg_kp_mermin_recovery_active_tail_max_delta = 1.0E-4_dp, &
                                        ot_cg_kp_mermin_recovery_tail_min_eps_factor = 1.5_dp, &
                                        ot_cg_kp_mermin_recovery_tail_max_eps_factor = 20.0_dp, &
                                        ot_cg_kp_mermin_recovery_step_floor = 1.0E-9_dp, &
                                        ot_cg_kp_mermin_recovery_active_step_floor = 1.0E-5_dp, &
                                        ot_cg_kp_physical_raw_delta_floor_ratio = 1.0E-2_dp, &
                                        ot_cg_kp_physical_raw_delta_mermin_floor = 1.0E-8_dp, &
                                        ot_cg_kp_physical_raw_delta_tail_quality = 9.9E-1_dp, &
                                        ot_cg_kp_physical_raw_delta_tail_gap = 1.0E-2_dp, &
                                        ot_cg_kp_physical_control_eps_factor = 3.0_dp, &
                                        ot_cg_kp_physical_control_quality = 8.0E-1_dp, &
                                        ot_cg_kp_physical_control_gap = 5.0E-2_dp, &
                                        ot_cg_kp_fermi_line_search_span = 4.0_dp, &
                                        ot_cg_kp_rayleigh_response_min_alignment = 9.0E-1_dp, &
                                        ot_cg_kp_rayleigh_response_min_energy_fraction = 9.95E-1_dp, &
                                        ot_cg_kp_rayleigh_response_min_scale = 2.5E-1_dp, &
                                        ot_cg_kp_rayleigh_response_max_scale = 4.0_dp, &
                                        ot_cg_kp_residual_trust_min_scale = 1.0E-1_dp, &
                                        ot_cg_kp_residual_trust_quality_floor = 1.0E-4_dp, &
                                        ot_cg_kp_residual_history_damp_quality = 2.0E-1_dp, &
                                        ot_cg_kp_residual_spike_quality_ceiling = 0.0_dp, &
                                        ot_cg_kp_closure_response_ema_alpha = 2.5E-1_dp, &
                                        ot_cg_kp_mermin_closure_backtrack = 5.0E-1_dp, &
                                        ot_cg_kp_mermin_secant_consistency = 2.0E-1_dp, &
                                        ot_diis_kp_rotation_response_min_work_ratio = 2.5E-1_dp, &
                                        ot_diis_kp_accepted_response_recovery = 2.0_dp

   TYPE, PRIVATE :: ot_cg_history_snapshot_type
      TYPE(dbcsr_type), POINTER                         :: matrix_gx_old => NULL(), &
                                                           matrix_gx_old_im => NULL(), &
                                                           rot_mat_gx_old => NULL(), &
                                                           rot_mat_gx_old_im => NULL()
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE          :: ener_gx_old
   END TYPE ot_cg_history_snapshot_type

CONTAINS

! **************************************************************************************************
!> \brief combine minimizer-specific and complete fixed-N physical residuals for SCF exit
!>
!>        CG computes a projected physical residual while constructing its search direction,
!>        whereas DIIS and Broyden do not traverse that code path.  The complete direct audit is
!>        minimizer independent and must therefore remain visible to every endpoint decision.
!> \param minimizer_delta physical residual supplied by the active minimizer
!> \param complete_delta complete fixed-N physical-gradient RMS, or a negative sentinel
!> \return conservative physical residual for an SCF endpoint decision
! **************************************************************************************************
   PURE FUNCTION ot_kpoint_physical_exit_delta(minimizer_delta, complete_delta) RESULT(physical_delta)
      REAL(KIND=dp), INTENT(IN)                          :: minimizer_delta, complete_delta
      REAL(KIND=dp)                                      :: physical_delta

      physical_delta = MAX(0.0_dp, minimizer_delta)
      IF (complete_delta >= 0.0_dp) physical_delta = MAX(physical_delta, complete_delta)

   END FUNCTION ot_kpoint_physical_exit_delta

! **************************************************************************************************
!> \brief reject a finite DIIS step that raises the fixed-N Mermin objective
!> \param current_energy Mermin free energy at the evaluated candidate
!> \param reference_energy Mermin free energy at the preceding accepted endpoint
!> \return whether the candidate lies outside the accepted-energy trust window
! **************************************************************************************************
   PURE FUNCTION ot_kpoint_diis_mermin_step_rejected(current_energy, reference_energy) RESULT(rejected)
      REAL(KIND=dp), INTENT(IN)                          :: current_energy, reference_energy
      LOGICAL                                            :: rejected

      REAL(KIND=dp)                                      :: tolerance

      tolerance = 1.0E-11_dp*MAX(1.0_dp, ABS(current_energy), ABS(reference_energy))
      rejected = current_energy > reference_energy + tolerance

   END FUNCTION ot_kpoint_diis_mermin_step_rejected

! **************************************************************************************************
!> \brief reject an energy-unresolved DIIS step that leaves the fixed-N residual trust region
!> \param current_energy Mermin free energy at the evaluated candidate
!> \param reference_energy Mermin free energy at the preceding accepted endpoint
!> \param current_delta complete fixed-N physical residual at the candidate
!> \param reference_delta complete fixed-N physical residual at the accepted endpoint
!> \param eps_scf requested SCF threshold
!> \return whether the candidate has unresolved energy and excessive tail-residual growth
! **************************************************************************************************
   PURE FUNCTION ot_kpoint_diis_residual_step_rejected( &
      current_energy, reference_energy, current_delta, reference_delta, eps_scf) RESULT(rejected)
      REAL(KIND=dp), INTENT(IN)                          :: current_energy, reference_energy, &
                                                            current_delta, reference_delta, eps_scf
      LOGICAL                                            :: rejected

      REAL(KIND=dp)                                      :: energy_tolerance, tail_max, tail_min

      rejected = .FALSE.
      IF (current_delta <= 0.0_dp .OR. reference_delta <= 0.0_dp .OR. eps_scf <= 0.0_dp) RETURN
      energy_tolerance = 1.0E-11_dp*MAX(1.0_dp, ABS(current_energy), ABS(reference_energy))
      IF (ABS(current_energy - reference_energy) > energy_tolerance) RETURN

      tail_min = kpoint_eps_scaled_tail_min( &
                 ot_cg_kp_mermin_recovery_tail_min_delta, &
                 ot_cg_kp_mermin_recovery_tail_min_eps_factor, eps_scf)
      tail_max = kpoint_eps_scaled_tail_max( &
                 ot_cg_kp_mermin_recovery_tail_max_delta, &
                 ot_cg_kp_mermin_recovery_tail_max_eps_factor, eps_scf, tail_min)
      rejected = reference_delta <= tail_max .AND. current_delta >= tail_min .AND. &
                 current_delta/reference_delta > ot_cg_kp_mermin_residual_growth_ratio

   END FUNCTION ot_kpoint_diis_residual_step_rejected

! **************************************************************************************************
!> \brief update the base DIIS/SD trust step after rejecting an evaluated candidate
!>
!>        Pulay extrapolation and the preconditioned SD direction have separate trust semantics.
!>        A failed extrapolation discards its Pulay history but leaves the validated base step
!>        unchanged.  A failed base-direction candidate contracts that step before retrying.
!> \param previous_step current signed base step
!> \param initial_step user-requested signed base step used to define the numerical floor
!> \param eps_taylor Taylor threshold used by the minimizer
!> \param base_candidate true when the rejected candidate was generated from the SD base direction
!> \return retained or contracted signed base step
! **************************************************************************************************
   PURE FUNCTION ot_kpoint_diis_rejected_step( &
      previous_step, initial_step, eps_taylor, base_candidate) RESULT(next_step)
      REAL(KIND=dp), INTENT(IN)                          :: previous_step, initial_step, eps_taylor
      LOGICAL, INTENT(IN)                                :: base_candidate
      REAL(KIND=dp)                                      :: next_step

      REAL(KIND=dp)                                      :: step_floor

      next_step = previous_step
      IF (.NOT. base_candidate) RETURN

      step_floor = MAX(1.0E4_dp*EPSILON(1.0_dp), eps_taylor)* &
                   MAX(1.0_dp, ABS(initial_step))
      next_step = SIGN(MAX(step_floor, 0.25_dp*ABS(previous_step)), previous_step)

   END FUNCTION ot_kpoint_diis_rejected_step

! **************************************************************************************************
!> \brief recover a contracted DIIS base step after a validated finite-response candidate
!>
!>        A rejected base candidate contracts the DIIS trust step, but that transient state may
!>        otherwise survive indefinitely, including across exact minimizer restarts.  Recovery is
!>        permitted only when an accepted finite complex response realizes a useful fraction of
!>        its predicted Mermin decrease, remains inside the complete fixed-N residual trust region,
!>        and does not worsen the rotation residual.  One recovery undoes half of a factor-four
!>        rejection contraction; callers must restart Pulay history when the scale changes.
!> \param previous_step current signed base step
!> \param initial_step user-requested signed base step and recovery ceiling
!> \param accepted_drop measured Mermin decrease of the accepted candidate
!> \param model_work first-order Mermin work predicted for that candidate
!> \param physical_ratio complete fixed-N residual ratio at the accepted endpoint
!> \param rotation_ratio complex rotation-residual ratio at the accepted endpoint
!> \param response_candidate true when the accepted candidate used the finite response
!> \return retained or recovered signed base step
! **************************************************************************************************
   PURE FUNCTION ot_kpoint_diis_accepted_response_step( &
      previous_step, initial_step, accepted_drop, model_work, physical_ratio, rotation_ratio, &
      response_candidate) RESULT(next_step)
      REAL(KIND=dp), INTENT(IN)                          :: previous_step, initial_step, &
                                                            accepted_drop, model_work, &
                                                            physical_ratio, rotation_ratio
      LOGICAL, INTENT(IN)                                :: response_candidate
      REAL(KIND=dp)                                      :: next_step

      LOGICAL                                            :: model_valid, residual_valid, step_valid

      next_step = previous_step
      step_valid = response_candidate .AND. previous_step*initial_step > 0.0_dp .AND. &
                   ABS(previous_step) < HUGE(previous_step) .AND. &
                   ABS(initial_step) < HUGE(initial_step) .AND. &
                   ABS(previous_step) < ABS(initial_step)*(1.0_dp - 10.0_dp*EPSILON(1.0_dp))
      IF (.NOT. step_valid) RETURN

      model_valid = accepted_drop > 0.0_dp .AND. model_work > EPSILON(model_work) .AND. &
                    accepted_drop < HUGE(accepted_drop) .AND. model_work < HUGE(model_work) .AND. &
                    accepted_drop >= ot_diis_kp_rotation_response_min_work_ratio*model_work
      residual_valid = physical_ratio > 0.0_dp .AND. rotation_ratio > 0.0_dp .AND. &
                       physical_ratio < HUGE(physical_ratio) .AND. &
                       rotation_ratio < HUGE(rotation_ratio) .AND. &
                       physical_ratio <= SQRT(ot_cg_kp_mermin_residual_growth_ratio) .AND. &
                       rotation_ratio <= 1.0_dp + 10.0_dp*EPSILON(rotation_ratio)
      IF (.NOT. model_valid .OR. .NOT. residual_valid) RETURN

      next_step = SIGN(MIN(ABS(initial_step), &
                           ot_diis_kp_accepted_response_recovery*ABS(previous_step)), previous_step)

   END FUNCTION ot_kpoint_diis_accepted_response_step

! **************************************************************************************************
!> \brief decide whether a current Schur rotation response is a valid DIIS descent correction
!>
!>        The Pulay metric retains the raw complex rotation residual.  A finite Schur response is
!>        applied only as a current-coordinate correction and only when its complete distributed
!>        contraction with the raw gradient is a finite descent direction that retains a useful
!>        fraction of the raw first-order work.  The correction remains disabled while an external
!>        Fermi-frontier projection owns the trial state.
!> \param direction_ready every distributed k-point channel supplied the current Schur response
!> \param response_slope weighted gradient-response contraction for every fixed-N group
!> \param raw_slope weighted squared raw rotation-gradient norm for every fixed-N group
!> \param external_projection_active an external frontier projection still controls trial states
!> \return true when the response can replace the current raw rotation shift
! **************************************************************************************************
   PURE FUNCTION ot_kpoint_diis_rotation_response_ready( &
      direction_ready, response_slope, raw_slope, external_projection_active) RESULT(ready)
      LOGICAL, INTENT(IN)                                :: direction_ready
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: response_slope, raw_slope
      LOGICAL, INTENT(IN)                                :: external_projection_active
      LOGICAL                                            :: ready

      INTEGER                                            :: group

      ready = direction_ready .AND. .NOT. external_projection_active .AND. &
              SIZE(response_slope) > 0 .AND. SIZE(response_slope) == SIZE(raw_slope)
      IF (.NOT. ready) RETURN
      DO group = 1, SIZE(response_slope)
         ready = ABS(response_slope(group)) < HUGE(response_slope(group)) .AND. &
                 ABS(raw_slope(group)) < HUGE(raw_slope(group)) .AND. &
                 raw_slope(group) > EPSILON(1.0_dp) .AND. &
                 response_slope(group) >= &
                 ot_diis_kp_rotation_response_min_work_ratio*raw_slope(group) .AND. &
                 response_slope(group) > &
                 EPSILON(1.0_dp)*MAX(raw_slope(group), ABS(response_slope(group)))
         IF (.NOT. ready) RETURN
      END DO

   END FUNCTION ot_kpoint_diis_rotation_response_ready

! **************************************************************************************************
!> \brief decide whether a mature multi-group DIIS state should probe the finite rotation response
!>
!>        The previous accepted Mermin decrease is compared with the complete first-order work of
!>        the current group-resolved response.  A probe is useful only when the raw/Pulay path
!>        underdelivers that local response model while the physical residual remains unresolved.
!> \param direction_ready every fixed-N group supplied a valid finite response
!> \param response_slope weighted response work for every fixed-N group
!> \param raw_slope weighted raw rotation work for every fixed-N group
!> \param accepted_samples accepted endpoints accumulated in the current AUTO subspace
!> \param accepted_drop Mermin decrease of the preceding accepted candidate
!> \param step current base trust step
!> \param physical_delta complete fixed-N residual at the accepted candidate
!> \param eps_scf requested SCF threshold
!> \return true when the next candidate should audit a finite-response correction
! **************************************************************************************************
   PURE FUNCTION ot_kpoint_diis_rotation_response_probe_ready( &
      direction_ready, response_slope, raw_slope, accepted_samples, accepted_drop, step, &
      physical_delta, eps_scf) RESULT(ready)
      LOGICAL, INTENT(IN)                                :: direction_ready
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: response_slope, raw_slope
      INTEGER, INTENT(IN)                                :: accepted_samples
      REAL(KIND=dp), INTENT(IN)                          :: accepted_drop, step, physical_delta, &
                                                            eps_scf
      LOGICAL                                            :: ready

      REAL(KIND=dp)                                      :: response_work, response_work_threshold

      ready = ot_kpoint_diis_rotation_response_ready( &
              direction_ready, response_slope, raw_slope, .FALSE.)
      ready = ready .AND. SIZE(response_slope) > 1 .AND. &
              accepted_samples >= ot_diis_kp_response_probe_min_accepted .AND. &
              accepted_drop >= 0.0_dp .AND. eps_scf > 0.0_dp .AND. &
              physical_delta > eps_scf*(1.0_dp + 10.0_dp*EPSILON(1.0_dp))
      IF (.NOT. ready) RETURN

      response_work = ABS(step)*SUM(response_slope)
      response_work_threshold = ot_diis_kp_rotation_response_min_work_ratio*response_work
      ready = response_work > EPSILON(1.0_dp) .AND. &
              accepted_drop < response_work_threshold - &
              10.0_dp*EPSILON(1.0_dp)*MAX(1.0_dp, response_work_threshold)

   END FUNCTION ot_kpoint_diis_rotation_response_probe_ready

! **************************************************************************************************
!> \brief blend raw and finite-response DIIS rotations under an accepted-rotation model
!>
!>        The affine rotation-residual model uses the last accepted raw and finite-response
!>        candidates.  It retains the largest response fraction whose interpolated ratio does
!>        not exceed one.  Missing calibration deliberately permits a full response probe.
!> \param raw_ratio accepted rotation-residual ratio for the raw DIIS direction
!> \param response_ratio accepted rotation-residual ratio for the finite response direction
!> \return response fraction in the closed interval [0,1]
! **************************************************************************************************
   PURE FUNCTION ot_kpoint_diis_response_residual_blend(raw_ratio, response_ratio) RESULT(blend)
      REAL(KIND=dp), INTENT(IN)                          :: raw_ratio, response_ratio
      REAL(KIND=dp)                                      :: blend

      REAL(KIND=dp)                                      :: denominator

      blend = 1.0_dp
      IF (raw_ratio <= 0.0_dp .OR. response_ratio <= 0.0_dp) RETURN
      IF (ABS(raw_ratio) >= HUGE(raw_ratio) .OR. &
          ABS(response_ratio) >= HUGE(response_ratio)) RETURN
      IF (response_ratio <= 1.0_dp) RETURN
      IF (raw_ratio >= 1.0_dp) THEN
         blend = 0.0_dp
         RETURN
      END IF

      denominator = response_ratio - raw_ratio
      IF (denominator <= EPSILON(denominator)*MAX(1.0_dp, ABS(response_ratio), ABS(raw_ratio))) THEN
         blend = 0.0_dp
      ELSE
         blend = MAX(0.0_dp, MIN(1.0_dp, (1.0_dp - raw_ratio)/denominator))
      END IF

   END FUNCTION ot_kpoint_diis_response_residual_blend

! **************************************************************************************************
!> \brief recover the full-response residual ratio from an accepted blended candidate
!> \param raw_ratio accepted residual ratio for the raw direction
!> \param mixed_ratio measured residual ratio for the blended candidate
!> \param response_blend response fraction used by the blended candidate
!> \return affine estimate of the corresponding full-response residual ratio
! **************************************************************************************************
   PURE FUNCTION ot_kpoint_diis_response_endpoint_ratio( &
      raw_ratio, mixed_ratio, response_blend) RESULT(response_ratio)
      REAL(KIND=dp), INTENT(IN)                          :: raw_ratio, mixed_ratio, response_blend
      REAL(KIND=dp)                                      :: response_ratio

      response_ratio = mixed_ratio
      IF (raw_ratio <= 0.0_dp .OR. mixed_ratio <= 0.0_dp .OR. &
          response_blend <= SQRT(EPSILON(response_blend))) RETURN
      IF (ABS(raw_ratio) >= HUGE(raw_ratio) .OR. ABS(mixed_ratio) >= HUGE(mixed_ratio)) RETURN

      response_ratio = raw_ratio + (mixed_ratio - raw_ratio)/response_blend

   END FUNCTION ot_kpoint_diis_response_endpoint_ratio

! **************************************************************************************************
!> \brief Evaluate convergence of an accepted complex fixed-N Mermin endpoint.
!> \param search_delta residual retained by the minimizer
!> \param physical_delta complete fixed-N rotation/energy residual
!> \param physical_ready whether the physical residual covers every channel and rotation mode
!> \param projected_delta projected fixed-N residual
!> \param density_delta canonical-density stationarity residual
!> \param density_ready whether the canonical-density audit is complete
!> \param eps_scf requested SCF threshold
!> \param exit_delta largest complete endpoint residual
!> \param converged whether every required residual is complete and below threshold
! **************************************************************************************************
   PURE SUBROUTINE ot_kpoint_mermin_endpoint_status( &
      search_delta, physical_delta, physical_ready, projected_delta, density_delta, density_ready, &
      eps_scf, exit_delta, converged)
      REAL(KIND=dp), INTENT(IN)                          :: search_delta, physical_delta
      LOGICAL, INTENT(IN)                                :: physical_ready
      REAL(KIND=dp), INTENT(IN)                          :: projected_delta, density_delta
      LOGICAL, INTENT(IN)                                :: density_ready
      REAL(KIND=dp), INTENT(IN)                          :: eps_scf
      REAL(KIND=dp), INTENT(OUT)                         :: exit_delta
      LOGICAL, INTENT(OUT)                               :: converged

      exit_delta = MAX(0.0_dp, search_delta)
      converged = .FALSE.
      IF (eps_scf <= 0.0_dp) RETURN

      ! An incomplete reduced Schur audit is not evidence of convergence.  In
      ! particular, a restart may carry a small density/model residual while
      ! most complex rotation modes have not yet entered the active subspace.
      IF (.NOT. physical_ready .OR. .NOT. density_ready) RETURN
      exit_delta = MAX(exit_delta, MAX(0.0_dp, physical_delta), &
                       MAX(0.0_dp, projected_delta), MAX(0.0_dp, density_delta))
      converged = exit_delta < eps_scf

   END SUBROUTINE ot_kpoint_mermin_endpoint_status

! **************************************************************************************************
!> \brief Complete a rotation/energy Schur model for the actual nonlinear-CG direction.
!>
!>        The Schur solve models the fresh rotation/energy preconditioned direction.  The line
!>        search also carries the independently preconditioned orbital part and, for nonlinear CG,
!>        beta times the previous direction.  Use the exact final directional slope, add the
!>        missing fresh-direction metric, and retain the measured previous-direction curvature.
!>        The latter is the conjugate-direction approximation to the history quadratic; an unknown
!>        Hessian cross term is deliberately not invented.
!> \param partial_slope Schur rotation/energy slope before CG mixing
!> \param partial_quadratic Schur rotation/energy curvature before CG mixing
!> \param fresh_slope complete fresh preconditioned-gradient slope before CG mixing
!> \param final_slope exact descent slope of the final CG direction
!> \param beta nonlinear-CG history coefficient
!> \param previous_curvature measured curvature of the previous accepted complete direction
!> \param previous_valid whether previous_curvature is available and consistent
!> \param full_slope exact completed model slope
!> \param full_quadratic completed conjugate-direction curvature
!> \param full_drop completed unit-step quadratic decrease
!> \param optimum_position unconstrained positive quadratic minimum
! **************************************************************************************************
   PURE SUBROUTINE ot_kpoint_mermin_full_direction_model( &
      partial_slope, partial_quadratic, fresh_slope, final_slope, beta, previous_curvature, &
      previous_valid, full_slope, full_quadratic, full_drop, optimum_position)
      REAL(KIND=dp), INTENT(IN)                          :: partial_slope, partial_quadratic, &
                                                            fresh_slope, final_slope, beta, &
                                                            previous_curvature
      LOGICAL, INTENT(IN)                                :: previous_valid
      REAL(KIND=dp), INTENT(OUT)                         :: full_slope, full_quadratic, full_drop, &
                                                            optimum_position

      REAL(KIND=dp)                                      :: missing_fresh_quadratic

      full_slope = MAX(0.0_dp, final_slope)
      missing_fresh_quadratic = MAX(0.0_dp, fresh_slope - MAX(0.0_dp, partial_slope))
      full_quadratic = MAX(0.0_dp, partial_quadratic) + missing_fresh_quadratic
      IF (previous_valid .AND. beta > 0.0_dp .AND. previous_curvature > 0.0_dp) THEN
         full_quadratic = full_quadratic + beta*beta*previous_curvature
      END IF
      full_drop = full_slope - 0.5_dp*full_quadratic
      optimum_position = 0.0_dp
      IF (full_slope > EPSILON(MAX(1.0_dp, full_slope)) .AND. &
          full_quadratic > EPSILON(MAX(1.0_dp, full_quadratic))) THEN
         optimum_position = full_slope/full_quadratic
      END IF

   END SUBROUTINE ot_kpoint_mermin_full_direction_model

! **************************************************************************************************
!> \brief save nonlinear-CG gradient history before a residual-only line-search audit
!> \param qs_ot_env local OT channels
!> \param snapshot saved gradient history
! **************************************************************************************************
   SUBROUTINE capture_ot_cg_history(qs_ot_env, snapshot)
      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: qs_ot_env
      TYPE(ot_cg_history_snapshot_type), ALLOCATABLE, &
         DIMENSION(:), INTENT(OUT)                       :: snapshot

      INTEGER                                            :: ispin
      LOGICAL                                            :: do_ener, do_ks, do_rotation

      do_ks = qs_ot_env(1)%settings%ks
      do_ener = qs_ot_env(1)%settings%do_ener
      do_rotation = qs_ot_env(1)%settings%do_rotation
      ALLOCATE (snapshot(SIZE(qs_ot_env)))
      DO ispin = 1, SIZE(qs_ot_env)
         IF (do_ks) THEN
            ALLOCATE (snapshot(ispin)%matrix_gx_old)
            CALL dbcsr_create(snapshot(ispin)%matrix_gx_old, &
                              template=qs_ot_env(ispin)%matrix_gx_old)
            CALL dbcsr_copy(snapshot(ispin)%matrix_gx_old, &
                            qs_ot_env(ispin)%matrix_gx_old)
            IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
               ALLOCATE (snapshot(ispin)%matrix_gx_old_im)
               CALL dbcsr_create(snapshot(ispin)%matrix_gx_old_im, &
                                 template=qs_ot_env(ispin)%matrix_gx_old_im)
               CALL dbcsr_copy(snapshot(ispin)%matrix_gx_old_im, &
                               qs_ot_env(ispin)%matrix_gx_old_im)
            END IF
            IF (do_rotation) THEN
               ALLOCATE (snapshot(ispin)%rot_mat_gx_old)
               CALL dbcsr_create(snapshot(ispin)%rot_mat_gx_old, &
                                 template=qs_ot_env(ispin)%rot_mat_gx_old)
               CALL dbcsr_copy(snapshot(ispin)%rot_mat_gx_old, &
                               qs_ot_env(ispin)%rot_mat_gx_old)
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  ALLOCATE (snapshot(ispin)%rot_mat_gx_old_im)
                  CALL dbcsr_create(snapshot(ispin)%rot_mat_gx_old_im, &
                                    template=qs_ot_env(ispin)%rot_mat_gx_old_im)
                  CALL dbcsr_copy(snapshot(ispin)%rot_mat_gx_old_im, &
                                  qs_ot_env(ispin)%rot_mat_gx_old_im)
               END IF
            END IF
         END IF
         IF (do_ener) THEN
            ALLOCATE (snapshot(ispin)%ener_gx_old, SOURCE=qs_ot_env(ispin)%ener_gx_old)
         END IF
      END DO

   END SUBROUTINE capture_ot_cg_history

! **************************************************************************************************
!> \brief restore nonlinear-CG gradient history after a residual-only line-search audit
!> \param qs_ot_env local OT channels
!> \param snapshot saved gradient history
! **************************************************************************************************
   SUBROUTINE restore_ot_cg_history(qs_ot_env, snapshot)
      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: qs_ot_env
      TYPE(ot_cg_history_snapshot_type), ALLOCATABLE, &
         DIMENSION(:), INTENT(INOUT)                     :: snapshot

      INTEGER                                            :: ispin

      DO ispin = 1, SIZE(qs_ot_env)
         IF (ASSOCIATED(snapshot(ispin)%matrix_gx_old)) THEN
            CALL dbcsr_copy(qs_ot_env(ispin)%matrix_gx_old, &
                            snapshot(ispin)%matrix_gx_old)
            CALL dbcsr_release(snapshot(ispin)%matrix_gx_old)
            DEALLOCATE (snapshot(ispin)%matrix_gx_old)
         END IF
         IF (ASSOCIATED(snapshot(ispin)%matrix_gx_old_im)) THEN
            CALL dbcsr_copy(qs_ot_env(ispin)%matrix_gx_old_im, &
                            snapshot(ispin)%matrix_gx_old_im)
            CALL dbcsr_release(snapshot(ispin)%matrix_gx_old_im)
            DEALLOCATE (snapshot(ispin)%matrix_gx_old_im)
         END IF
         IF (ASSOCIATED(snapshot(ispin)%rot_mat_gx_old)) THEN
            CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_gx_old, &
                            snapshot(ispin)%rot_mat_gx_old)
            CALL dbcsr_release(snapshot(ispin)%rot_mat_gx_old)
            DEALLOCATE (snapshot(ispin)%rot_mat_gx_old)
         END IF
         IF (ASSOCIATED(snapshot(ispin)%rot_mat_gx_old_im)) THEN
            CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_gx_old_im, &
                            snapshot(ispin)%rot_mat_gx_old_im)
            CALL dbcsr_release(snapshot(ispin)%rot_mat_gx_old_im)
            DEALLOCATE (snapshot(ispin)%rot_mat_gx_old_im)
         END IF
         IF (ALLOCATED(snapshot(ispin)%ener_gx_old)) THEN
            qs_ot_env(ispin)%ener_gx_old(:) = snapshot(ispin)%ener_gx_old
         END IF
      END DO
      DEALLOCATE (snapshot)

   END SUBROUTINE restore_ot_cg_history

! **************************************************************************************************
!> \brief Decide whether a Schur-metric change invalidates nonlinear-CG history in the residual tail.
!> \param history_change incompatible response-metric transition was detected
!> \param beta proposed nonlinear-CG history coefficient
!> \param physical_tail_delta projected physical residual
!> \param residual_tail_min lower tail bound scaled to the requested EPS_SCF
!> \param calibration_change accepted-response candidate transition requiring an isolated probe
!> \return true when the previous conjugate direction must be discarded
! **************************************************************************************************
   PURE FUNCTION ot_kpoint_cg_metric_history_reset(history_change, beta, physical_tail_delta, &
                                                   residual_tail_min, calibration_change) &
      RESULT(reset_history)
      LOGICAL, INTENT(IN)                                :: history_change
      REAL(KIND=dp), INTENT(IN)                          :: beta, physical_tail_delta, &
                                                            residual_tail_min
      LOGICAL, INTENT(IN), OPTIONAL                      :: calibration_change
      LOGICAL                                            :: reset_history

      LOGICAL                                            :: isolate_calibration

      isolate_calibration = .FALSE.
      IF (PRESENT(calibration_change)) isolate_calibration = calibration_change
      reset_history = history_change .AND. beta > 0.0_dp .AND. &
                      (isolate_calibration .OR. &
                       (physical_tail_delta > residual_tail_min .AND. &
                        physical_tail_delta < ot_cg_kp_residual_tail_max_delta))

   END FUNCTION ot_kpoint_cg_metric_history_reset

! **************************************************************************************************
!> \brief Classify consecutive dense/fallback Schur metrics guarded by the residual-tail CG restart.
!> \param attempt current Schur attempt
!> \param last_attempt previous Schur attempt
!> \param flavor current dense candidate flavor, negative for the dense fallback
!> \param last_flavor previous dense candidate flavor
!> \param physical_tail_delta projected physical residual
!> \param residual_tail_min lower tail bound scaled to the requested EPS_SCF
!> \return true when both directions belong to one guarded tail metric family
! **************************************************************************************************
   PURE FUNCTION ot_kpoint_cg_tail_metric_tracked(attempt, last_attempt, flavor, last_flavor, &
                                                  physical_tail_delta, residual_tail_min) &
      RESULT(tracked)
      INTEGER, INTENT(IN)                                :: attempt, last_attempt, flavor, &
                                                            last_flavor
      REAL(KIND=dp), INTENT(IN)                          :: physical_tail_delta, residual_tail_min
      LOGICAL                                            :: tracked

      tracked = .FALSE.
      IF (attempt < 3 .OR. last_attempt < 3) RETURN
      tracked = flavor >= 30 .AND. last_flavor >= 30
      IF (flavor < 0 .AND. last_flavor < 0) THEN
         tracked = physical_tail_delta > residual_tail_min .AND. &
                   physical_tail_delta < ot_cg_kp_residual_tail_max_delta
      END IF

   END FUNCTION ot_kpoint_cg_tail_metric_tracked

! **************************************************************************************************
!> \brief separate a calibrated physical model residual from its conservative SCF exit floor
!> \param raw_delta preconditioned residual
!> \param physical_delta measured fixed-N physical residual
!> \param physical_ready whether all physical response groups are converged
!> \param eps_scf requested SCF threshold
!> \param quality_samples number of accepted Mermin calibration samples
!> \param accepted_quality accepted/predicted Mermin-drop quality
!> \param accepted_gap accepted-to-best line-search gap
!> \param model_delta residual used by the minimizer and response model
!> \param exit_floor_delta additional residual floor used only for the outer SCF exit
! **************************************************************************************************
   PURE SUBROUTINE ot_kpoint_physical_residual_split(raw_delta, physical_delta, physical_ready, &
                                                     eps_scf, quality_samples, accepted_quality, &
                                                     accepted_gap, model_delta, exit_floor_delta)
      REAL(KIND=dp), INTENT(IN)                          :: raw_delta, physical_delta
      LOGICAL, INTENT(IN)                                :: physical_ready
      REAL(KIND=dp), INTENT(IN)                          :: eps_scf
      INTEGER, INTENT(IN)                                :: quality_samples
      REAL(KIND=dp), INTENT(IN)                          :: accepted_quality, accepted_gap
      REAL(KIND=dp), INTENT(OUT)                         :: model_delta, exit_floor_delta

      REAL(KIND=dp)                                      :: physical_control_delta

      model_delta = MAX(0.0_dp, raw_delta)
      exit_floor_delta = 0.0_dp
      IF (physical_delta <= 0.0_dp .OR. physical_delta >= model_delta) RETURN

      IF (physical_ready) THEN
         model_delta = physical_delta
         RETURN
      END IF
      IF (eps_scf <= 0.0_dp .OR. quality_samples < ot_cg_kp_physical_control_min_samples .OR. &
          accepted_quality < ot_cg_kp_physical_control_quality .OR. &
          accepted_gap > ot_cg_kp_physical_control_gap) RETURN

      physical_control_delta = MAX(physical_delta, ot_cg_kp_physical_control_eps_factor*eps_scf)
      IF (physical_control_delta >= model_delta) RETURN
      model_delta = physical_delta
      exit_floor_delta = physical_control_delta

   END SUBROUTINE ot_kpoint_physical_residual_split

! **************************************************************************************************
!> \brief advance the auxiliary-energy coordinate inside an independent symmetric trust interval
!> \param position current auxiliary-energy line-search position
!> \param delta requested change of the shared orbital line-search position
!> \param bound nonnegative auxiliary-energy position bound, or a negative inactive marker
!> \param next_position bounded auxiliary-energy position
!> \param actual_delta change applied to the auxiliary energies
!> \param tangent finite-step fraction retained from the requested coordinate change
! **************************************************************************************************
   PURE SUBROUTINE ot_kpoint_fermi_window_energy_coordinate(position, delta, bound, &
                                                            next_position, actual_delta, tangent)
      REAL(KIND=dp), INTENT(IN)                          :: position, delta, bound
      REAL(KIND=dp), INTENT(OUT)                         :: next_position, actual_delta, tangent

      REAL(KIND=dp)                                      :: proposed, tolerance

      proposed = position + delta
      next_position = proposed
      tolerance = 64.0_dp*EPSILON(1.0_dp)*MAX(1.0_dp, ABS(position), ABS(delta), ABS(bound))
      IF (bound >= 0.0_dp) next_position = SIGN(MIN(ABS(proposed), bound), proposed)
      actual_delta = next_position - position
      IF (ABS(delta) > tolerance) THEN
         tangent = MAX(0.0_dp, MIN(1.0_dp, actual_delta/delta))
      ELSE
         tangent = 1.0_dp
      END IF

   END SUBROUTINE ot_kpoint_fermi_window_energy_coordinate

! **************************************************************************************************
!> \brief largest affine step fraction that keeps a current Fermi-frontier state in its window
!> \param distance current auxiliary energy relative to the fixed-N chemical potential
!> \param step full proposed relative-energy step after removal of the chemical-potential gauge
!> \param anchor maximum allowed absolute frontier distance
!> \return scale in [0,1], or -1 when the current state is outside the guarded frontier
! **************************************************************************************************
   PURE FUNCTION ot_kpoint_fermi_window_direction_scale(distance, step, anchor) RESULT(scale)
      REAL(KIND=dp), INTENT(IN)                          :: distance, step, anchor
      REAL(KIND=dp)                                      :: scale

      REAL(KIND=dp)                                      :: boundary, tolerance

      scale = -1.0_dp
      IF (anchor <= 0.0_dp) RETURN
      tolerance = 64.0_dp*EPSILON(1.0_dp)* &
                  MAX(anchor, ABS(distance), ABS(step), TINY(1.0_dp))
      IF (ABS(distance) > anchor + tolerance) RETURN

      boundary = MAX(0.0_dp, anchor - tolerance)
      IF (ABS(distance + step) <= boundary) THEN
         scale = 1.0_dp
      ELSE IF (step > 0.0_dp) THEN
         scale = (boundary - distance)/step
      ELSE IF (step < 0.0_dp) THEN
         scale = (-boundary - distance)/step
      ELSE
         scale = 1.0_dp
      END IF
      scale = MAX(0.0_dp, MIN(1.0_dp, scale))

   END FUNCTION ot_kpoint_fermi_window_direction_scale

! **************************************************************************************************
!> \brief classify accepted Mermin steps that invalidate the nonlinear-CG history
!> \param previous_delta coupled fixed-N residual before the accepted step
!> \param current_delta coupled fixed-N residual after the accepted step
!> \param accepted_drop measured Mermin free-energy decrease
!> \param eps_scf SCF convergence threshold
!> \param reset whether the old conjugate direction must be discarded
!> \param quality residual-consistency quality used to rescale the next trust radius
!> \param catastrophic whether the reset is outside the convergence-tail guard
!> \param strict_growth apply the ordinary trust ratio outside the convergence tail
! **************************************************************************************************
   PURE SUBROUTINE ot_kpoint_cg_residual_growth_guard(previous_delta, current_delta, accepted_drop, &
                                                      eps_scf, reset, quality, catastrophic, strict_growth)
      REAL(KIND=dp), INTENT(IN)                          :: previous_delta, current_delta, &
                                                            accepted_drop, eps_scf
      LOGICAL, INTENT(OUT)                               :: reset
      REAL(KIND=dp), INTENT(OUT)                         :: quality
      LOGICAL, INTENT(OUT)                               :: catastrophic
      LOGICAL, INTENT(IN), OPTIONAL                      :: strict_growth

      LOGICAL                                            :: strict, tail_growth
      REAL(KIND=dp)                                      :: growth_ratio, tail_max, tail_min

      reset = .FALSE.
      catastrophic = .FALSE.
      quality = 1.0_dp
      strict = .FALSE.
      IF (PRESENT(strict_growth)) strict = strict_growth
      IF (accepted_drop <= 0.0_dp .OR. &
          previous_delta <= EPSILON(1.0_dp) .OR. current_delta <= previous_delta) RETURN

      growth_ratio = current_delta/previous_delta
      tail_min = kpoint_eps_scaled_tail_min( &
                 ot_cg_kp_mermin_recovery_tail_min_delta, &
                 ot_cg_kp_mermin_recovery_tail_min_eps_factor, eps_scf)
      tail_max = kpoint_eps_scaled_tail_max( &
                 ot_cg_kp_mermin_recovery_tail_max_delta, &
                 ot_cg_kp_mermin_recovery_tail_max_eps_factor, eps_scf, tail_min)
      tail_growth = previous_delta <= tail_max .AND. current_delta >= tail_min .AND. &
                    current_delta <= ot_cg_kp_mermin_catastrophic_growth_ratio*tail_max .AND. &
                    growth_ratio > ot_cg_kp_mermin_residual_growth_ratio
      catastrophic = previous_delta > ot_cg_kp_mermin_residual_growth_max_delta .AND. &
                     growth_ratio > ot_cg_kp_mermin_catastrophic_growth_ratio
      reset = tail_growth .OR. catastrophic .OR. &
              (strict .AND. growth_ratio > ot_cg_kp_mermin_residual_growth_ratio)
      IF (reset) quality = MIN(1.0_dp, (previous_delta/current_delta)**2)
   END SUBROUTINE ot_kpoint_cg_residual_growth_guard

! **************************************************************************************************
!> \brief decode the immutable residual reference of one Mermin line search
!> \param model_delta search-model residual captured before the direction was applied
!> \param reference_delta signed stored reference; positive is projected, negative is search
!> \return positive residual used by accepted-endpoint closure
! **************************************************************************************************
   PURE FUNCTION ot_kpoint_mermin_reference_guard_delta(model_delta, reference_delta) RESULT(delta)
      REAL(KIND=dp), INTENT(IN)                          :: model_delta, reference_delta
      REAL(KIND=dp)                                      :: delta

      delta = model_delta
      IF (reference_delta /= 0.0_dp) delta = ABS(reference_delta)

   END FUNCTION ot_kpoint_mermin_reference_guard_delta

! **************************************************************************************************
!> \brief select the current residual in the same norm as the immutable line-search reference
!> \param search_delta effective minimizer/search residual
!> \param projected_delta fixed-N projected physical residual, or non-positive if unavailable
!> \param reference_delta signed stored reference; positive is projected, negative is search
!> \return norm-consistent current residual
! **************************************************************************************************
   PURE FUNCTION ot_kpoint_mermin_current_guard_delta(search_delta, projected_delta, reference_delta) &
      RESULT(delta)
      REAL(KIND=dp), INTENT(IN)                          :: search_delta, projected_delta, &
                                                            reference_delta
      REAL(KIND=dp)                                      :: delta

      delta = search_delta
      IF (reference_delta > 0.0_dp .AND. projected_delta > 0.0_dp) delta = projected_delta

   END FUNCTION ot_kpoint_mermin_current_guard_delta

! **************************************************************************************************
!> \brief detect whether distinct spin channels share one fixed-N occupation group
!> \param spin_index spin label of every local OT channel
!> \param group_index fixed-N occupation-group label of every local OT channel
!> \return whether a common chemical potential couples distinct spins
! **************************************************************************************************
   PURE FUNCTION ot_kpoint_mermin_shared_spin_group(spin_index, group_index) RESULT(shared)
      INTEGER, INTENT(IN)                                :: spin_index(:), group_index(:)
      LOGICAL                                            :: shared

      INTEGER                                            :: i, j

      shared = .FALSE.
      IF (SIZE(spin_index) /= SIZE(group_index)) RETURN
      DO i = 1, SIZE(spin_index)
         DO j = i + 1, SIZE(spin_index)
            IF (spin_index(i) /= spin_index(j) .AND. group_index(i) == group_index(j)) THEN
               shared = .TRUE.
               RETURN
            END IF
         END DO
      END DO

   END FUNCTION ot_kpoint_mermin_shared_spin_group

! **************************************************************************************************
!> \brief whether a projected physical residual is a valid immutable line-search reference
!> \param physical_ready complete physical audit already satisfies its stationarity gate
!> \param physical_groups number of complete fixed-N groups in the current endpoint audit
!> \param shared_spin_group whether distinct spin channels share one fixed-N group
!> \return whether the projected residual has the required physical group coverage
! **************************************************************************************************
   PURE FUNCTION ot_kpoint_mermin_projected_reference_ready(physical_ready, physical_groups, &
                                                            shared_spin_group) RESULT(ready)
      LOGICAL, INTENT(IN)                                :: physical_ready
      INTEGER, INTENT(IN)                                :: physical_groups
      LOGICAL, INTENT(IN)                                :: shared_spin_group
      LOGICAL                                            :: ready

      ! Common-mu UKS has one fixed-N group containing both spin channels.  Group count is endpoint
      ! coverage, not line-search history, so one rebuilt group is sufficient for that coupling.
      ready = physical_ready .OR. (physical_groups > 0 .AND. shared_spin_group)

   END FUNCTION ot_kpoint_mermin_projected_reference_ready

! **************************************************************************************************
!> \brief identify a tail response whose energy completion shares the exact fixed-N closure model
!> \param fresh_response the Schur response was rebuilt for the current accepted orbital state
!> \param closure_attempted the fixed-N closure gate has been reached in this subspace epoch
!> \param model_flavor selected coupled or energy-only closure flavor
!> \return true when the Rayleigh audit must use the contracted residual-aware search interval
! **************************************************************************************************
   PURE FUNCTION ot_kpoint_residual_aware_rayleigh_model( &
      fresh_response, closure_attempted, model_flavor) RESULT(active)
      LOGICAL, INTENT(IN)                                :: fresh_response, closure_attempted
      INTEGER, INTENT(IN)                                :: model_flavor
      LOGICAL                                            :: active

      active = fresh_response .AND. closure_attempted .AND. &
               (model_flavor == ot_cg_kp_coupled_closure_flavor .OR. &
                model_flavor == ot_cg_kp_energy_closure_flavor)

   END FUNCTION ot_kpoint_residual_aware_rayleigh_model

! **************************************************************************************************
!> \brief line-search span used to calibrate the finite fixed-orbital Rayleigh response
!> \param current_step current residual-aware line-search step
!> \param configured_step configured uncontracted line-search step
!> \param residual_aware use the current contracted step when it is available
!> \return positive affine position spanning the complete 2PNT search
! **************************************************************************************************
   PURE FUNCTION ot_kpoint_rayleigh_probe_position(current_step, configured_step, residual_aware) &
      RESULT(position)
      REAL(KIND=dp), INTENT(IN)                          :: current_step, configured_step
      LOGICAL, INTENT(IN)                                :: residual_aware
      REAL(KIND=dp)                                      :: position

      IF (residual_aware .AND. ABS(current_step) > EPSILON(current_step)) THEN
         position = ot_cg_kp_fermi_line_search_span*ABS(current_step)
      ELSE
         position = ot_cg_kp_fermi_line_search_span* &
                    MAX(ABS(current_step), ABS(configured_step))
      END IF
      position = MAX(position, SQRT(EPSILON(position)))

   END FUNCTION ot_kpoint_rayleigh_probe_position

! **************************************************************************************************
!> \brief calibrate a fixed-N Schur energy response against the exact fixed-orbital Rayleigh target
!> \param target_norm2 weighted squared norm of the Rayleigh energy correction
!> \param response_norm2 weighted squared norm of the finite Schur response over the LS span
!> \param target_response_dot weighted scalar product of target and response
!> \param energy_fraction fraction of the directional Mermin slope in the energy channel
!> \param minimum_target_norm2 dimensionally scaled floor for a thermally inactive target
!> \param maximum_scale optional upper bound for line-search reparameterization
!> \return bounded least-squares response scale, or one for a non-collinear/invalid response
! **************************************************************************************************
   PURE FUNCTION ot_kpoint_rayleigh_response_scale( &
      target_norm2, response_norm2, target_response_dot, energy_fraction, &
      minimum_target_norm2, maximum_scale) RESULT(scale)
      REAL(KIND=dp), INTENT(IN)                          :: target_norm2, response_norm2, &
                                                            target_response_dot, energy_fraction, &
                                                            minimum_target_norm2
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: maximum_scale
      REAL(KIND=dp)                                      :: scale

      REAL(KIND=dp)                                      :: alignment, norm_product, scale_upper

      scale = 1.0_dp
      IF (energy_fraction < ot_cg_kp_rayleigh_response_min_energy_fraction) RETURN
      IF (target_norm2 <= MAX(0.0_dp, minimum_target_norm2) .OR. response_norm2 <= 0.0_dp) RETURN
      IF (target_response_dot <= 0.0_dp) RETURN
      ! Form the product after taking both square roots.  The weighted residual
      ! norms are dimensionful and can legitimately be far below machine
      ! epsilon in a tight SCF tail; multiplying them first can also underflow.
      norm_product = SQRT(target_norm2)*SQRT(response_norm2)
      IF (norm_product <= 0.0_dp) RETURN
      alignment = target_response_dot/norm_product
      IF (alignment < ot_cg_kp_rayleigh_response_min_alignment) RETURN

      scale = target_response_dot/response_norm2
      scale_upper = ot_cg_kp_rayleigh_response_max_scale
      IF (PRESENT(maximum_scale)) scale_upper = MAX(scale_upper, maximum_scale)
      scale = MIN(scale_upper, &
                  MAX(ot_cg_kp_rayleigh_response_min_scale, scale))

   END FUNCTION ot_kpoint_rayleigh_response_scale

! **************************************************************************************************
!> \brief transform one directional model component after groupwise Rayleigh response scaling
!> \param energy_slope nonnegative directional slopes in the auxiliary-energy channel
!> \param orbital_slope nonnegative directional slopes in all unscaled orbital channels
!> \param response_scale groupwise auxiliary-energy response scales
!> \param power one for the linear model and two for its asymptotic quadratic model
!> \return weighted model factor, or one for inconsistent/inactive input
! **************************************************************************************************
   PURE FUNCTION ot_kpoint_rayleigh_model_factor(energy_slope, orbital_slope, response_scale, power) &
      RESULT(factor)
      REAL(KIND=dp), INTENT(IN)                          :: energy_slope(:), orbital_slope(:), &
                                                            response_scale(:)
      INTEGER, INTENT(IN)                                :: power
      REAL(KIND=dp)                                      :: factor

      REAL(KIND=dp)                                      :: denominator

      factor = 1.0_dp
      IF (SIZE(energy_slope) /= SIZE(orbital_slope) .OR. &
          SIZE(energy_slope) /= SIZE(response_scale)) RETURN
      IF (power < 1 .OR. SIZE(energy_slope) == 0) RETURN
      IF (ANY(energy_slope < 0.0_dp) .OR. ANY(orbital_slope < 0.0_dp) .OR. &
          ANY(response_scale <= 0.0_dp)) RETURN

      denominator = SUM(energy_slope + orbital_slope)
      IF (denominator <= 0.0_dp) RETURN
      factor = SUM(orbital_slope + response_scale**power*energy_slope)/denominator

   END FUNCTION ot_kpoint_rayleigh_model_factor

! **************************************************************************************************
!> \brief retain a smaller trust radius after a line search returned to its reference point
!> \param default_step configured line-search step
!> \param positions sampled and next line-search positions
!> \param energies sampled Mermin free energies
!> \param count number of evaluated samples
!> \param step_floor numerical lower bound
!> \param force_contract contract after a residual-rejected closure even if sampled energies descended
!> \return ...
! **************************************************************************************************
   PURE FUNCTION ot_kpoint_mermin_restart_step(default_step, positions, energies, count, &
                                               step_floor, force_contract) RESULT(restart_step)
      REAL(KIND=dp), INTENT(IN)                          :: default_step, positions(:), energies(:)
      INTEGER, INTENT(IN)                                :: count
      REAL(KIND=dp), INTENT(IN)                          :: step_floor
      LOGICAL, INTENT(IN), OPTIONAL                      :: force_contract
      REAL(KIND=dp)                                      :: restart_step

      INTEGER                                            :: i, nsample
      LOGICAL                                            :: failed_reference, forced
      REAL(KIND=dp)                                      :: min_positive, tolerance

      restart_step = default_step
      nsample = MIN(count, SIZE(energies), SIZE(positions) - 1)
      IF (nsample < 3) RETURN
      tolerance = 1.0E-12_dp*MAX(1.0_dp, ABS(energies(1)))
      forced = .FALSE.
      IF (PRESENT(force_contract)) forced = force_contract
      failed_reference = (forced .OR. ALL(energies(2:nsample) >= energies(1) - tolerance)) .AND. &
                         ABS(positions(nsample + 1)) <= step_floor
      IF (.NOT. failed_reference) RETURN

      min_positive = HUGE(1.0_dp)
      DO i = 2, nsample
         IF (ABS(positions(i)) > step_floor) THEN
            min_positive = MIN(min_positive, ABS(positions(i)))
         END IF
      END DO
      IF (min_positive < HUGE(1.0_dp)) THEN
         restart_step = MAX(step_floor, ot_cg_kp_mermin_closure_backtrack*min_positive)
      END IF

   END FUNCTION ot_kpoint_mermin_restart_step

! **************************************************************************************************
!> \brief choose the initial line-search scale after changing the response model
!> \param default_step configured line-search step
!> \param contracted_step step inherited from the completed line search
!> \param previous_flavor response-model flavor of the completed search
!> \param current_flavor response-model flavor of the new direction
!> \return initial step for the new search
! **************************************************************************************************
   PURE FUNCTION ot_kpoint_mermin_candidate_restart_step(default_step, contracted_step, &
                                                         previous_flavor, current_flavor) &
      RESULT(restart_step)
      REAL(KIND=dp), INTENT(IN)                          :: default_step, contracted_step
      INTEGER, INTENT(IN)                                :: previous_flavor, current_flavor
      REAL(KIND=dp)                                      :: restart_step

      restart_step = contracted_step
      ! The finite-response Direct-density direction has its own accepted-step
      ! calibration.  A contracted step from another response model is not a
      ! meaningful trust radius for its first line search.
      IF (current_flavor == 32 .AND. previous_flavor /= 32) restart_step = ABS(default_step)

   END FUNCTION ot_kpoint_mermin_candidate_restart_step

! **************************************************************************************************
!> \brief detect a Mermin line search whose next point and recent energies no longer change
!> \param positions sampled and next line-search positions
!> \param energies sampled Mermin free energies
!> \param count number of evaluated samples
!> \return whether the sampled search is numerically stagnant
! **************************************************************************************************
   PURE FUNCTION ot_kpoint_mermin_search_stagnant(positions, energies, count) RESULT(stagnant)
      REAL(KIND=dp), INTENT(IN)                          :: positions(:), energies(:)
      INTEGER, INTENT(IN)                                :: count
      LOGICAL                                            :: stagnant

      INTEGER                                            :: first, nsample
      REAL(KIND=dp)                                      :: energy_scale, energy_span, position_scale

      stagnant = .FALSE.
      nsample = MIN(count, SIZE(energies), SIZE(positions) - 1)
      IF (nsample < 3) RETURN
      first = MAX(1, nsample - 2)
      energy_scale = MAX(1.0_dp, MAXVAL(ABS(energies(first:nsample))))
      energy_span = MAXVAL(energies(first:nsample)) - MINVAL(energies(first:nsample))
      position_scale = MAX(1.0_dp, ABS(positions(nsample)), ABS(positions(nsample + 1)))
      stagnant = energy_span <= 1.0E-12_dp*energy_scale .AND. &
                 ABS(positions(nsample + 1) - positions(nsample)) <= 1.0E-12_dp*position_scale

   END FUNCTION ot_kpoint_mermin_search_stagnant

! **************************************************************************************************
!> \brief detect an active physical tail whose recent Mermin energies are at line-search resolution
!> \param energies sampled Mermin free energies
!> \param count number of evaluated samples
!> \return whether five recent energies are numerically indistinguishable
! **************************************************************************************************
   PURE FUNCTION ot_kpoint_mermin_tail_energy_stagnant(energies, count) RESULT(stagnant)
      REAL(KIND=dp), INTENT(IN)                          :: energies(:)
      INTEGER, INTENT(IN)                                :: count
      LOGICAL                                            :: stagnant

      INTEGER                                            :: first, nsample
      REAL(KIND=dp)                                      :: energy_scale, energy_span

      stagnant = .FALSE.
      nsample = MIN(count, SIZE(energies))
      IF (nsample < 5) RETURN
      first = nsample - 4
      energy_scale = MAX(1.0_dp, MAXVAL(ABS(energies(first:nsample))))
      energy_span = MAXVAL(energies(first:nsample)) - MINVAL(energies(first:nsample))
      stagnant = energy_span <= 1.0E-12_dp*energy_scale

   END FUNCTION ot_kpoint_mermin_tail_energy_stagnant

! **************************************************************************************************
!> \brief refine a bracketed Mermin line-search endpoint from directional derivatives
!>
!>        The sampled Mermin energies can become less discriminating than the analytic gradient
!>        in the fixed-N residual tail.  A line search which already brackets zero slope therefore
!>        gets one safeguarded secant refinement when the selected endpoint still violates the
!>        strong-Wolfe curvature condition.  The returned point stays strictly inside the accepted
!>        reference-to-endpoint interval.
!> \param reference_slope directional derivative at position zero
!> \param current_slope directional derivative at the selected endpoint
!> \param current_position selected endpoint position
!> \param position_floor minimum resolvable displacement from either interval boundary
!> \param refined_position safeguarded secant position, or current_position when not applicable
!> \param refine whether a derivative refinement should be evaluated
! **************************************************************************************************
   PURE SUBROUTINE ot_kpoint_mermin_derivative_refinement( &
      reference_slope, current_slope, current_position, position_floor, refined_position, refine)
      REAL(KIND=dp), INTENT(IN)                          :: reference_slope, current_slope, &
                                                            current_position, position_floor
      REAL(KIND=dp), INTENT(OUT)                         :: refined_position
      LOGICAL, INTENT(OUT)                               :: refine

      REAL(KIND=dp), PARAMETER :: strong_wolfe_curvature = 0.5_dp

      REAL(KIND=dp)                                      :: denominator, floor, secant_position

      refine = .FALSE.
      refined_position = current_position
      floor = MAX(0.0_dp, position_floor)
      IF (current_position <= 2.0_dp*floor) RETURN
      IF (reference_slope >= 0.0_dp .OR. current_slope <= 0.0_dp) RETURN
      IF (ABS(current_slope) <= strong_wolfe_curvature*ABS(reference_slope)) RETURN

      denominator = current_slope - reference_slope
      IF (denominator <= EPSILON(denominator)* &
          MAX(1.0_dp, ABS(reference_slope), ABS(current_slope))) RETURN
      secant_position = -reference_slope*current_position/denominator
      IF (secant_position <= floor .OR. secant_position >= current_position - floor) RETURN

      refined_position = secant_position
      refine = .TRUE.

   END SUBROUTINE ot_kpoint_mermin_derivative_refinement

! **************************************************************************************************
!> \brief keep an open Mermin search until its selected endpoint has a physical residual
!> \param sample_pending a selected line-search sample still awaits closure
!> \param line_search_might_be_done the energy interpolation selected an endpoint
!> \param energy_only no gradient or physical residual was evaluated at the current endpoint
!> \param history_exhausted no additional sample fits in the line-search arrays
!> \param stagnant_history recent energies and positions are numerically unchanged
!> \return whether collapsed-search recovery must be deferred
! **************************************************************************************************
   PURE FUNCTION ot_kpoint_mermin_recovery_deferred(sample_pending, line_search_might_be_done, &
                                                    energy_only, history_exhausted, stagnant_history) &
      RESULT(deferred)
      LOGICAL, INTENT(IN)                                :: sample_pending, &
                                                            line_search_might_be_done, &
                                                            energy_only, history_exhausted, &
                                                            stagnant_history
      LOGICAL                                            :: deferred

      deferred = sample_pending .AND. line_search_might_be_done .AND. &
                 .NOT. history_exhausted .AND. &
                 (.NOT. energy_only .OR. .NOT. stagnant_history)

   END FUNCTION ot_kpoint_mermin_recovery_deferred
!
! the minimizer interface
! should present all possible modes of minimization
! these include CG SD DIIS
!
!
! IN the case of nspin != 1 we have a gradient that is distributed over different qs_ot_env.
! still things remain basically the same, since there are no constraints between the different qs_ot_env
! we only should take care that the various scalar products are taken over the full vectors.
! all the information needed and collected can be stored in the fist qs_ot_env only
! (indicating that the data type for the gradient/position and minization should be separated)
!
! **************************************************************************************************
!> \brief ...
!> \param qs_ot_env ...
!> \param matrix_hc ...
!> \param matrix_hc_im ...
!> \param para_env_inter_kp ...
! **************************************************************************************************
   SUBROUTINE ot_mini(qs_ot_env, matrix_hc, matrix_hc_im, para_env_inter_kp)
      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: qs_ot_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_hc
      TYPE(dbcsr_p_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: matrix_hc_im
      TYPE(mp_para_env_type), OPTIONAL, POINTER          :: para_env_inter_kp

      CHARACTER(len=*), PARAMETER                        :: routineN = 'ot_mini'

      INTEGER                                            :: handle, ispin, nspin, &
                                                            previous_search_flavor
      LOGICAL :: candidate_search_restart, contracted_search_restart, current_state_converged, &
         do_ener, do_ks, residual_closure_rejected
      REAL(KIND=dp) :: current_exit_delta, direction_energy_max, direction_matrix_im2, &
         direction_matrix_re2, direction_rotation_im2, direction_rotation_re2, next_search_step, &
         step_floor, tmp
      TYPE(cp_logger_type), POINTER                      :: logger

      CALL timeset(routineN, handle)

      nspin = SIZE(qs_ot_env)

      do_ks = qs_ot_env(1)%settings%ks
      do_ener = qs_ot_env(1)%settings%do_ener
      residual_closure_rejected = .FALSE.
      contracted_search_restart = .FALSE.

      CALL sync_kpoint_minimizer_control(qs_ot_env, para_env_inter_kp)
      logger => cp_get_default_logger()
      next_search_step = qs_ot_env(1)%settings%ds_min
      previous_search_flavor = qs_ot_env(1)%schur_cg_flavor
      CALL recover_collapsed_mermin_search(qs_ot_env)

      IF (qs_ot_env(1)%has_complex_kpoint_state) THEN
         logger => cp_get_default_logger()
         IF (logger%iter_info%print_level >= debug_print_level) THEN
            direction_matrix_re2 = 0.0_dp
            direction_matrix_im2 = 0.0_dp
            direction_rotation_re2 = 0.0_dp
            direction_rotation_im2 = 0.0_dp
            direction_energy_max = 0.0_dp
            DO ispin = 1, nspin
               IF (ASSOCIATED(qs_ot_env(ispin)%matrix_dx)) THEN
                  CALL dbcsr_dot(qs_ot_env(ispin)%matrix_dx, qs_ot_env(ispin)%matrix_dx, tmp)
                  direction_matrix_re2 = direction_matrix_re2 + tmp
               END IF
               IF (ASSOCIATED(qs_ot_env(ispin)%matrix_dx_im)) THEN
                  CALL dbcsr_dot(qs_ot_env(ispin)%matrix_dx_im, qs_ot_env(ispin)%matrix_dx_im, tmp)
                  direction_matrix_im2 = direction_matrix_im2 + tmp
               END IF
               IF (ASSOCIATED(qs_ot_env(ispin)%rot_mat_dx)) THEN
                  CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_dx, qs_ot_env(ispin)%rot_mat_dx, tmp)
                  direction_rotation_re2 = direction_rotation_re2 + tmp
               END IF
               IF (ASSOCIATED(qs_ot_env(ispin)%rot_mat_dx_im)) THEN
                  CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_dx_im, qs_ot_env(ispin)%rot_mat_dx_im, tmp)
                  direction_rotation_im2 = direction_rotation_im2 + tmp
               END IF
               IF (ASSOCIATED(qs_ot_env(ispin)%ener_dx)) THEN
                  direction_energy_max = MAX(direction_energy_max, &
                                             MAXVAL(ABS(qs_ot_env(ispin)%ener_dx)))
               END IF
            END DO
            IF (logger%para_env%is_source()) THEN
               WRITE (cp_logger_get_default_unit_nr(logger), '(A,1X,I5,5(1X,L1),6(1X,ES16.8))') &
                  " OT line-search direction state count use_dx energy_only done linear complex "// &
                  "mre mim rre rim emax ds:", &
                  qs_ot_env(1)%line_search_count, qs_ot_env(1)%use_dx, qs_ot_env(1)%energy_only, &
                  qs_ot_env(1)%line_search_might_be_done, qs_ot_env(1)%ener_direction_linear, &
                  qs_ot_env(1)%has_complex_kpoint_state, &
                  SQRT(MAX(0.0_dp, direction_matrix_re2)), &
                  SQRT(MAX(0.0_dp, direction_matrix_im2)), &
                  SQRT(MAX(0.0_dp, direction_rotation_re2)), &
                  SQRT(MAX(0.0_dp, direction_rotation_im2)), &
                  direction_energy_max, qs_ot_env(1)%ds_min
            END IF
         END IF
      END IF

      qs_ot_env(1)%OT_METHOD_FULL = ""

      ! compute the gradient for the variables x
      IF (.NOT. qs_ot_env(1)%energy_only) THEN
         qs_ot_env(1)%gradient = 0.0_dp
         DO ispin = 1, nspin
            IF (do_ks) THEN
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  IF (.NOT. PRESENT(matrix_hc_im)) THEN
                     CPABORT("Complex k-point OT derivative requires imaginary H(k)*C(k).")
                  END IF
                  SELECT CASE (qs_ot_env(1)%settings%ot_algorithm)
                  CASE ("REF")
                     CALL qs_ot_get_derivative_ref_complex(matrix_hc(ispin)%matrix, &
                                                           matrix_hc_im(ispin)%matrix, &
                                                           qs_ot_env(ispin))
                  CASE DEFAULT
                     CPABORT("Complex k-point OT derivative is implemented for REF only")
                  END SELECT
               ELSE
                  SELECT CASE (qs_ot_env(1)%settings%ot_algorithm)
                  CASE ("TOD")
                     CALL qs_ot_get_derivative(matrix_hc(ispin)%matrix, qs_ot_env(ispin)%matrix_x, &
                                               qs_ot_env(ispin)%matrix_sx, &
                                               qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin))
                  CASE ("REF")
                     CALL qs_ot_get_derivative_ref(matrix_hc(ispin)%matrix, &
                                                   qs_ot_env(ispin)%matrix_x, qs_ot_env(ispin)%matrix_sx, &
                                                   qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin))
                  CASE DEFAULT
                     CPABORT("ALGORITHM NYI")
                  END SELECT
               END IF
            END IF
            ! and also the gradient along the direction
            IF (qs_ot_env(1)%use_dx) THEN
               IF (do_ks) THEN
                  CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin)%matrix_dx, tmp)
                  qs_ot_env(1)%gradient = qs_ot_env(1)%gradient + tmp
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx_im, qs_ot_env(ispin)%matrix_dx_im, tmp)
                     qs_ot_env(1)%gradient = qs_ot_env(1)%gradient + tmp
                  END IF
                  IF (qs_ot_env(1)%settings%do_rotation) THEN
                     CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx, qs_ot_env(ispin)%rot_mat_dx, tmp)
                     qs_ot_env(1)%gradient = qs_ot_env(1)%gradient + 0.5_dp*tmp
                     IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                        CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx_im, qs_ot_env(ispin)%rot_mat_dx_im, tmp)
                        qs_ot_env(1)%gradient = qs_ot_env(1)%gradient + 0.5_dp*tmp
                     END IF
                  END IF
               END IF
               IF (do_ener) THEN
                  tmp = DOT_PRODUCT(qs_ot_env(ispin)%ener_gx, qs_ot_env(ispin)%ener_dx)
                  qs_ot_env(1)%gradient = qs_ot_env(1)%gradient + &
                                          qs_ot_env(ispin)%ener_x_line_search_tangent*tmp
               END IF
            ELSE
               IF (do_ks) THEN
                  CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin)%matrix_gx, tmp)
                  qs_ot_env(1)%gradient = qs_ot_env(1)%gradient - tmp
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx_im, qs_ot_env(ispin)%matrix_gx_im, tmp)
                     qs_ot_env(1)%gradient = qs_ot_env(1)%gradient - tmp
                  END IF
                  IF (qs_ot_env(1)%settings%do_rotation) THEN
                     CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx, qs_ot_env(ispin)%rot_mat_gx, tmp)
                     qs_ot_env(1)%gradient = qs_ot_env(1)%gradient - 0.5_dp*tmp
                     IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                        CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx_im, qs_ot_env(ispin)%rot_mat_gx_im, tmp)
                        qs_ot_env(1)%gradient = qs_ot_env(1)%gradient - 0.5_dp*tmp
                     END IF
                  END IF
               END IF
               IF (do_ener) THEN
                  tmp = DOT_PRODUCT(qs_ot_env(ispin)%ener_gx, qs_ot_env(ispin)%ener_gx)
                  qs_ot_env(1)%gradient = qs_ot_env(1)%gradient - tmp
               END IF
            END IF
         END DO
         CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, qs_ot_env(1)%gradient)
         CALL audit_kpoint_mermin_accepted_secant(qs_ot_env, para_env_inter_kp)
      END IF

      IF (do_ener) THEN
         DO ispin = 1, nspin
            qs_ot_env(ispin)%ener_x_trusted_step_max = 0.0_dp
            qs_ot_env(ispin)%ener_x_trust_scale = 1.0_dp
         END DO
      END IF

      IF (qs_ot_env(1)%has_complex_kpoint_state .AND. &
          qs_ot_env(1)%mermin_ls_sample_pending .AND. &
          qs_ot_env(1)%line_search_count >= 3 .AND. &
          .NOT. qs_ot_env(1)%energy_only) THEN
         CALL ot_new_cg_direction(qs_ot_env, para_env_inter_kp, residual_closure_rejected, &
                                  residual_preflight_only=.TRUE.)
      END IF
      CALL recover_stagnant_mermin_preflight(qs_ot_env, residual_closure_rejected)

      CALL update_kpoint_mermin_closure_response(qs_ot_env, para_env_inter_kp, &
                                                 residual_closure_rejected)

      current_state_converged = .FALSE.
      current_exit_delta = HUGE(1.0_dp)
      IF (qs_ot_env(1)%has_complex_kpoint_state .AND. &
          qs_ot_env(1)%kpoint_eps_scf > 0.0_dp .AND. &
          current_point_is_fine(qs_ot_env)) THEN
         CALL ot_kpoint_mermin_endpoint_status( &
            qs_ot_env(1)%delta, ot_kpoint_physical_exit_delta( &
            qs_ot_env(1)%kpoint_physical_delta, qs_ot_env(1)%kpoint_complete_physical_delta), &
            qs_ot_env(1)%kpoint_physical_delta_ready, qs_ot_env(1)%kpoint_projected_delta, &
            qs_ot_env(1)%kpoint_density_stationarity_delta, &
            qs_ot_env(1)%kpoint_density_stationarity_ready, qs_ot_env(1)%kpoint_eps_scf, &
            current_exit_delta, current_state_converged)
      END IF
      IF (current_state_converged) THEN
         ! Do not launch another Mermin line search from an already converged,
         ! explicitly evaluated endpoint. A trial taken here would expose its
         ! stale pre-step residual to the outer SCF convergence check.
         CALL reset_kpoint_mermin_search(qs_ot_env, qs_ot_env(1)%settings%ds_min)
         qs_ot_env(1)%OT_METHOD_FULL = "OT FIN"
         IF (logger%iter_info%print_level >= debug_print_level .AND. &
             logger%para_env%is_source()) THEN
            WRITE (cp_logger_get_default_unit_nr(logger), '(A,2(1X,ES16.8))') &
               " OT Mermin accepted endpoint converged delta eps:", &
               current_exit_delta, qs_ot_env(1)%kpoint_eps_scf
         END IF
         CALL timestop(handle)
         RETURN
      END IF

      IF (qs_ot_env(1)%has_complex_kpoint_state .AND. &
          qs_ot_env(1)%line_search_might_be_done .AND. &
          qs_ot_env(1)%line_search_count >= 3) THEN
         step_floor = MAX(1.0E4_dp*EPSILON(1.0_dp), &
                          qs_ot_env(1)%settings%eps_taylor)* &
                      MAX(1.0_dp, ABS(qs_ot_env(1)%settings%ds_min))
         next_search_step = ot_kpoint_mermin_restart_step(qs_ot_env(1)%settings%ds_min, &
                                                          qs_ot_env(1)%ot_pos, &
                                                          qs_ot_env(1)%ot_energy, &
                                                          qs_ot_env(1)%line_search_count, step_floor, &
                                                          force_contract= &
                                                          qs_ot_env(1)%mermin_ls_closure_rejected)
         contracted_search_restart = next_search_step < &
                                     ABS(qs_ot_env(1)%settings%ds_min)
      END IF

      SELECT CASE (qs_ot_env(1)%settings%OT_METHOD)
      CASE ("CG")
         IF (current_point_is_fine(qs_ot_env)) THEN
            qs_ot_env(1)%OT_METHOD_FULL = "OT CG"
            CALL ot_new_cg_direction(qs_ot_env, para_env_inter_kp, residual_closure_rejected)
            IF (residual_closure_rejected) THEN
               qs_ot_env(1)%OT_METHOD_FULL = "OT LS"
            ELSE
               CALL store_kpoint_mermin_direction_slopes(qs_ot_env, para_env_inter_kp)
               next_search_step = ot_kpoint_mermin_candidate_restart_step( &
                                  qs_ot_env(1)%settings%ds_min, next_search_step, &
                                  previous_search_flavor, qs_ot_env(1)%schur_model_flavor)
               candidate_search_restart = qs_ot_env(1)%schur_model_flavor == 32 .AND. &
                                          previous_search_flavor /= 32
               IF (contracted_search_restart .OR. candidate_search_restart) THEN
                  qs_ot_env(1)%ds_min = next_search_step
               END IF
               qs_ot_env(1)%line_search_count = 0
            END IF
         ELSE
            qs_ot_env(1)%OT_METHOD_FULL = "OT LS"
         END IF
         CALL do_line_search(qs_ot_env)
      CASE ("SD")
         IF (current_point_is_fine(qs_ot_env)) THEN
            qs_ot_env(1)%OT_METHOD_FULL = "OT SD"
            CALL ot_new_sd_direction(qs_ot_env, para_env_inter_kp)
            CALL store_kpoint_mermin_direction_slopes(qs_ot_env, para_env_inter_kp)
            qs_ot_env(1)%line_search_count = 0
         ELSE
            qs_ot_env(1)%OT_METHOD_FULL = "OT LS"
         END IF
         CALL do_line_search(qs_ot_env)
      CASE ("DIIS")
         qs_ot_env(1)%OT_METHOD_FULL = "OT DIIS"
         CALL ot_diis_step(qs_ot_env, para_env_inter_kp)
      CASE ("BROY")
         qs_ot_env(1)%OT_METHOD_FULL = "OT BROY"
         CALL ot_broyden_step(qs_ot_env)
      CASE ("LBFG")
         IF (current_point_is_fine(qs_ot_env)) THEN
            qs_ot_env(1)%OT_METHOD_FULL = "OT LBFGS"
            CALL ot_new_lbfgs_direction(qs_ot_env, para_env_inter_kp)
            qs_ot_env(1)%line_search_count = 0
         ELSE
            qs_ot_env(1)%OT_METHOD_FULL = "OT LS"
         END IF
         CALL do_line_search(qs_ot_env)
      CASE DEFAULT
         CPABORT("OT_METHOD NYI")
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE ot_mini

! **************************************************************************************************
!> \brief Return the bounded signed residual reduction between a model and its closed endpoint.
!> \param model_delta residual retained with the accepted line-search model
!> \param current_delta residual after density/Hamiltonian closure
!> \return signed reduction in [-1, 1]
! **************************************************************************************************
   FUNCTION ot_kpoint_mermin_closure_response_quality(model_delta, current_delta) RESULT(quality)
      REAL(KIND=dp), INTENT(IN)                          :: model_delta, current_delta
      REAL(KIND=dp)                                      :: quality

      quality = 0.0_dp
      IF (model_delta <= EPSILON(model_delta) .OR. current_delta < 0.0_dp) RETURN
      quality = (model_delta - current_delta)/model_delta
      quality = MAX(-1.0_dp, MIN(1.0_dp, quality))

   END FUNCTION ot_kpoint_mermin_closure_response_quality

! **************************************************************************************************
!> \brief Audit response quality only after the accepted endpoint has completed residual closure.
!> \param qs_ot_env local OT channels
!> \param para_env_inter_kp communicator connecting k-point groups
!> \param residual_closure_rejected whether the endpoint failed the residual preflight
! **************************************************************************************************
   SUBROUTINE update_kpoint_mermin_closure_response(qs_ot_env, para_env_inter_kp, &
                                                    residual_closure_rejected)
      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: qs_ot_env
      TYPE(mp_para_env_type), OPTIONAL, POINTER          :: para_env_inter_kp
      LOGICAL, INTENT(IN)                                :: residual_closure_rejected

      INTEGER                                            :: samples
      REAL(KIND=dp)                                      :: current_delta, density_current_delta, &
                                                            density_quality, quality
      TYPE(cp_logger_type), POINTER                      :: logger

      IF (.NOT. qs_ot_env(1)%schur_closure_response_pending) RETURN
      IF (residual_closure_rejected .OR. qs_ot_env(1)%mermin_ls_closure_rejected) THEN
         qs_ot_env(1)%schur_closure_response_pending = .FALSE.
         RETURN
      END IF
      IF (.NOT. qs_ot_env(1)%mermin_ls_sample_pending) RETURN
      IF (.NOT. qs_ot_env(1)%line_search_might_be_done) RETURN
      IF (qs_ot_env(1)%line_search_count < 3) RETURN
      IF (.NOT. current_point_is_fine(qs_ot_env)) RETURN

      current_delta = MAX(0.0_dp, qs_ot_env(1)%delta)
      CALL ot_mini_max_kpoint_scalar(qs_ot_env, para_env_inter_kp, current_delta)
      quality = ot_kpoint_mermin_closure_response_quality( &
                qs_ot_env(1)%schur_closure_response_model_delta, current_delta)

      density_current_delta = -1.0_dp
      density_quality = 0.0_dp
      IF (qs_ot_env(1)%kpoint_density_stationarity_ready .AND. &
          qs_ot_env(1)%schur_closure_response_model_density_delta > EPSILON(1.0_dp)) THEN
         density_current_delta = MAX(0.0_dp, qs_ot_env(1)%kpoint_density_stationarity_delta)
         CALL ot_mini_max_kpoint_scalar(qs_ot_env, para_env_inter_kp, density_current_delta)
         density_quality = ot_kpoint_mermin_closure_response_quality( &
                           qs_ot_env(1)%schur_closure_response_model_density_delta, &
                           density_current_delta)
      END IF

      samples = MIN(qs_ot_env(1)%schur_closure_response_samples + 1, HUGE(samples) - 1)
      qs_ot_env(1)%schur_closure_response_samples = samples
      qs_ot_env(1)%schur_closure_response_quality = quality
      qs_ot_env(1)%schur_closure_response_density_quality = density_quality
      IF (samples == 1) THEN
         qs_ot_env(1)%schur_closure_response_quality_ema = quality
         qs_ot_env(1)%schur_closure_response_density_quality_ema = density_quality
      ELSE
         qs_ot_env(1)%schur_closure_response_quality_ema = &
            (1.0_dp - ot_cg_kp_closure_response_ema_alpha)* &
            qs_ot_env(1)%schur_closure_response_quality_ema + &
            ot_cg_kp_closure_response_ema_alpha*quality
         qs_ot_env(1)%schur_closure_response_density_quality_ema = &
            (1.0_dp - ot_cg_kp_closure_response_ema_alpha)* &
            qs_ot_env(1)%schur_closure_response_density_quality_ema + &
            ot_cg_kp_closure_response_ema_alpha*density_quality
      END IF

      logger => cp_get_default_logger()
      IF (logger%iter_info%print_level >= debug_print_level .AND. logger%para_env%is_source()) THEN
         WRITE (cp_logger_get_default_unit_nr(logger), '(A,1X,I6,8(1X,ES16.8))') &
            " OT Mermin closed-response audit n model current q qema dmodel dcurrent dq dqema:", &
            samples, qs_ot_env(1)%schur_closure_response_model_delta, current_delta, &
            quality, qs_ot_env(1)%schur_closure_response_quality_ema, &
            qs_ot_env(1)%schur_closure_response_model_density_delta, density_current_delta, &
            density_quality, qs_ot_env(1)%schur_closure_response_density_quality_ema
      END IF
      qs_ot_env(1)%schur_closure_response_pending = .FALSE.

   END SUBROUTINE update_kpoint_mermin_closure_response

! **************************************************************************************************
!> \brief synchronize the scalar minimizer state shared by distributed k-point channels
!> \param qs_ot_env local OT channels
!> \param para_env_inter_kp communicator connecting k-point groups
! **************************************************************************************************
   SUBROUTINE sync_kpoint_minimizer_control(qs_ot_env, para_env_inter_kp)
      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: qs_ot_env
      TYPE(mp_para_env_type), OPTIONAL, POINTER          :: para_env_inter_kp

      IF (ASSOCIATED(qs_ot_env(1)%para_env)) THEN
         CALL bcast_kpoint_minimizer_control(qs_ot_env, qs_ot_env(1)%para_env)
      END IF
      IF (PRESENT(para_env_inter_kp)) THEN
         IF (ASSOCIATED(para_env_inter_kp)) THEN
            CALL bcast_kpoint_minimizer_control(qs_ot_env, para_env_inter_kp)
         END IF
      END IF

   END SUBROUTINE sync_kpoint_minimizer_control

! **************************************************************************************************
!> \brief broadcast shared scalar minimizer state on one k-point communicator
!> \param qs_ot_env local OT channels
!> \param para_env communicator used for this synchronization pass
! **************************************************************************************************
   SUBROUTINE bcast_kpoint_minimizer_control(qs_ot_env, para_env)
      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: qs_ot_env
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL para_env%bcast(qs_ot_env(1)%ot_pos)
      CALL para_env%bcast(qs_ot_env(1)%ot_energy)
      CALL para_env%bcast(qs_ot_env(1)%ot_grad)
      CALL para_env%bcast(qs_ot_env(1)%line_search_left)
      CALL para_env%bcast(qs_ot_env(1)%line_search_right)
      CALL para_env%bcast(qs_ot_env(1)%line_search_mid)
      CALL para_env%bcast(qs_ot_env(1)%line_search_count)
      CALL para_env%bcast(qs_ot_env(1)%line_search_might_be_done)
      CALL para_env%bcast(qs_ot_env(1)%mermin_ls_sample_pending)
      CALL para_env%bcast(qs_ot_env(1)%mermin_ls_closure_rejected)
      CALL para_env%bcast(qs_ot_env(1)%mermin_ls_best_sample_index)
      CALL para_env%bcast(qs_ot_env(1)%mermin_ls_reference_energy)
      CALL para_env%bcast(qs_ot_env(1)%mermin_ls_reference_delta)
      CALL para_env%bcast(qs_ot_env(1)%mermin_ls_best_sample_energy)
      CALL para_env%bcast(qs_ot_env(1)%mermin_ls_best_sample_drop)
      CALL para_env%bcast(qs_ot_env(1)%mermin_ls_model_drop)
      CALL para_env%bcast(qs_ot_env(1)%mermin_ls_model_quality)
      CALL para_env%bcast(qs_ot_env(1)%mermin_ls_best_sample_pos)
      CALL para_env%bcast(qs_ot_env(1)%mermin_ls_next_pos)
      CALL para_env%bcast(qs_ot_env(1)%mermin_ls_accepted_quality_samples)
      CALL para_env%bcast(qs_ot_env(1)%kpoint_subspace_accepted_samples)
      CALL para_env%bcast(qs_ot_env(1)%diis_last_candidate_response)
      CALL para_env%bcast(qs_ot_env(1)%diis_last_response_blend)
      CALL para_env%bcast(qs_ot_env(1)%diis_last_response_model_work)
      CALL para_env%bcast(qs_ot_env(1)%diis_reference_rotation_gnorm)
      CALL para_env%bcast(qs_ot_env(1)%diis_raw_residual_ratio)
      CALL para_env%bcast(qs_ot_env(1)%diis_response_residual_ratio)
      CALL para_env%bcast(qs_ot_env(1)%mermin_ls_accepted_quality)
      CALL para_env%bcast(qs_ot_env(1)%mermin_ls_accepted_quality_ema)
      CALL para_env%bcast(qs_ot_env(1)%mermin_ls_accepted_gap)
      CALL para_env%bcast(qs_ot_env(1)%mermin_ls_accepted_gap_ema)
      ! Mermin component slopes are channel-local endpoint data.  Broadcasting
      ! channel one here replaces the other irreducible-k-point contributions
      ! before the accepted-secant audit sums them.
      CALL para_env%bcast(qs_ot_env(1)%delta)
      CALL para_env%bcast(qs_ot_env(1)%gnorm)
      CALL para_env%bcast(qs_ot_env(1)%gnorm_old)
      CALL para_env%bcast(qs_ot_env(1)%etotal)
      CALL para_env%bcast(qs_ot_env(1)%gradient)
      CALL para_env%bcast(qs_ot_env(1)%energy_only)
      CALL para_env%bcast(qs_ot_env(1)%use_dx)
      CALL para_env%bcast(qs_ot_env(1)%use_gx_old)
      CALL para_env%bcast(qs_ot_env(1)%OT_count)
      CALL para_env%bcast(qs_ot_env(1)%ds_min)
      CALL para_env%bcast(qs_ot_env(1)%schur_model_attempt)
      CALL para_env%bcast(qs_ot_env(1)%schur_model_flavor)
      CALL para_env%bcast(qs_ot_env(1)%schur_cg_attempt)
      CALL para_env%bcast(qs_ot_env(1)%schur_cg_flavor)
      CALL para_env%bcast(qs_ot_env(1)%schur_model_delta)
      CALL para_env%bcast(qs_ot_env(1)%schur_model_control_delta)
      CALL para_env%bcast(qs_ot_env(1)%schur_model_density_stationarity_delta)
      CALL para_env%bcast(qs_ot_env(1)%schur_density_stationarity_quality_samples)
      CALL para_env%bcast(qs_ot_env(1)%schur_density_stationarity_quality)
      CALL para_env%bcast(qs_ot_env(1)%schur_density_stationarity_quality_ema)
      CALL para_env%bcast(qs_ot_env(1)%schur_closure_response_pending)
      CALL para_env%bcast(qs_ot_env(1)%schur_closure_response_samples)
      CALL para_env%bcast(qs_ot_env(1)%schur_closure_response_model_delta)
      CALL para_env%bcast(qs_ot_env(1)%schur_closure_response_model_density_delta)
      CALL para_env%bcast(qs_ot_env(1)%schur_closure_response_quality)
      CALL para_env%bcast(qs_ot_env(1)%schur_closure_response_quality_ema)
      CALL para_env%bcast(qs_ot_env(1)%schur_closure_response_density_quality)
      CALL para_env%bcast(qs_ot_env(1)%schur_closure_response_density_quality_ema)
      CALL para_env%bcast(qs_ot_env(1)%schur_closure_shadow_pending)
      CALL para_env%bcast(qs_ot_env(1)%schur_closure_shadow_quality_samples)
      CALL para_env%bcast(qs_ot_env(1)%schur_closure_shadow_baseline_quadratic)
      CALL para_env%bcast(qs_ot_env(1)%schur_closure_shadow_quadratic)
      CALL para_env%bcast(qs_ot_env(1)%schur_closure_shadow_quality)
      CALL para_env%bcast(qs_ot_env(1)%schur_closure_shadow_quality_ema)
      CALL para_env%bcast(qs_ot_env(1)%schur_closure_shadow_advantage)
      CALL para_env%bcast(qs_ot_env(1)%schur_closure_shadow_advantage_ema)
      CALL para_env%bcast(qs_ot_env(1)%schur_closure_shadow_residual_quality)
      CALL para_env%bcast(qs_ot_env(1)%schur_closure_shadow_residual_quality_ema)
      CALL para_env%bcast(qs_ot_env(1)%schur_closure_shadow_drop_sum)
      CALL para_env%bcast(qs_ot_env(1)%schur_closure_shadow_measured_sum)
      CALL para_env%bcast(qs_ot_env(1)%schur_projector_density_full_selected_stationarity_samples)
      CALL para_env%bcast(qs_ot_env(1)%schur_projector_density_full_selected_stationarity_quality)
      CALL para_env%bcast(qs_ot_env(1)%schur_projector_density_full_selected_stationarity_quality_ema)
      CALL para_env%bcast(qs_ot_env(1)%schur_residual_history_reset_pending)
      CALL para_env%bcast(qs_ot_env(1)%schur_residual_history_quality)
      CALL para_env%bcast(qs_ot_env(1)%kpoint_physical_delta_ready)
      CALL para_env%bcast(qs_ot_env(1)%kpoint_density_stationarity_ready)
      CALL para_env%bcast(qs_ot_env(1)%kpoint_physical_delta_groups)
      CALL para_env%bcast(qs_ot_env(1)%kpoint_raw_delta)
      CALL para_env%bcast(qs_ot_env(1)%kpoint_physical_delta)
      CALL para_env%bcast(qs_ot_env(1)%kpoint_complete_physical_delta)
      CALL para_env%bcast(qs_ot_env(1)%kpoint_projected_delta)
      CALL para_env%bcast(qs_ot_env(1)%kpoint_density_stationarity_delta)
      CALL para_env%bcast(qs_ot_env(1)%kpoint_density_stationarity_max)
      CALL para_env%bcast(qs_ot_env(1)%kpoint_eps_scf)
      CALL para_env%bcast(qs_ot_env(1)%ener_direction_linear)

   END SUBROUTINE bcast_kpoint_minimizer_control

! **************************************************************************************************
!> \brief retain component slopes for an accepted-step Mermin secant
!> \param qs_ot_env ...
!> \param para_env_inter_kp ...
! **************************************************************************************************
   SUBROUTINE store_kpoint_mermin_direction_slopes(qs_ot_env, para_env_inter_kp)
      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: qs_ot_env
      TYPE(mp_para_env_type), OPTIONAL, POINTER          :: para_env_inter_kp

      INTEGER                                            :: ichannel
      LOGICAL                                            :: frozen_valid, vector_valid
      REAL(KIND=dp)                                      :: energy_slope, matrix_slope, &
                                                            rotation_slope, slope_sum, tmp
      TYPE(cp_logger_type), POINTER                      :: logger

      DO ichannel = 1, SIZE(qs_ot_env)
         qs_ot_env(ichannel)%mermin_ls_component_slope_valid = .FALSE.
         qs_ot_env(ichannel)%mermin_ls_vector_start_valid = .FALSE.
         qs_ot_env(ichannel)%mermin_frozen_h_valid = .FALSE.
         qs_ot_env(ichannel)%mermin_frozen_endpoint_valid = .FALSE.
         qs_ot_env(ichannel)%mermin_ls_matrix_slope = 0.0_dp
         qs_ot_env(ichannel)%mermin_ls_rotation_slope = 0.0_dp
         qs_ot_env(ichannel)%mermin_ls_energy_slope = 0.0_dp
      END DO
      logger => cp_get_default_logger()
      IF (.NOT. use_three_point_mermin_search(qs_ot_env)) RETURN
      IF (.NOT. qs_ot_env(1)%use_dx) RETURN

      matrix_slope = 0.0_dp
      rotation_slope = 0.0_dp
      energy_slope = 0.0_dp
      DO ichannel = 1, SIZE(qs_ot_env)
         IF (.NOT. qs_ot_env(ichannel)%has_complex_kpoint_state) CYCLE
         IF (qs_ot_env(1)%settings%ks) THEN
            IF (.NOT. ASSOCIATED(qs_ot_env(ichannel)%matrix_gx) .OR. &
                .NOT. ASSOCIATED(qs_ot_env(ichannel)%matrix_dx) .OR. &
                .NOT. ASSOCIATED(qs_ot_env(ichannel)%matrix_gx_im) .OR. &
                .NOT. ASSOCIATED(qs_ot_env(ichannel)%matrix_dx_im)) CYCLE
            CALL dbcsr_dot(qs_ot_env(ichannel)%matrix_gx, qs_ot_env(ichannel)%matrix_dx, tmp)
            qs_ot_env(ichannel)%mermin_ls_matrix_slope = tmp
            CALL dbcsr_dot(qs_ot_env(ichannel)%matrix_gx_im, qs_ot_env(ichannel)%matrix_dx_im, tmp)
            qs_ot_env(ichannel)%mermin_ls_matrix_slope = &
               qs_ot_env(ichannel)%mermin_ls_matrix_slope + tmp
         END IF
         IF (qs_ot_env(1)%settings%do_rotation) THEN
            IF (.NOT. ASSOCIATED(qs_ot_env(ichannel)%rot_mat_gx) .OR. &
                .NOT. ASSOCIATED(qs_ot_env(ichannel)%rot_mat_dx) .OR. &
                .NOT. ASSOCIATED(qs_ot_env(ichannel)%rot_mat_gx_im) .OR. &
                .NOT. ASSOCIATED(qs_ot_env(ichannel)%rot_mat_dx_im)) CYCLE
            CALL dbcsr_dot(qs_ot_env(ichannel)%rot_mat_gx, qs_ot_env(ichannel)%rot_mat_dx, tmp)
            qs_ot_env(ichannel)%mermin_ls_rotation_slope = 0.5_dp*tmp
            CALL dbcsr_dot(qs_ot_env(ichannel)%rot_mat_gx_im, qs_ot_env(ichannel)%rot_mat_dx_im, tmp)
            qs_ot_env(ichannel)%mermin_ls_rotation_slope = &
               qs_ot_env(ichannel)%mermin_ls_rotation_slope + 0.5_dp*tmp
         END IF
         IF (qs_ot_env(1)%settings%do_ener) THEN
            IF (.NOT. ASSOCIATED(qs_ot_env(ichannel)%ener_gx) .OR. &
                .NOT. ASSOCIATED(qs_ot_env(ichannel)%ener_dx)) CYCLE
            qs_ot_env(ichannel)%mermin_ls_energy_slope = &
               DOT_PRODUCT(qs_ot_env(ichannel)%ener_gx, qs_ot_env(ichannel)%ener_dx)
         END IF
         vector_valid = qs_ot_env(1)%settings%do_rotation .OR. qs_ot_env(1)%settings%do_ener
         IF (qs_ot_env(1)%settings%do_rotation) THEN
            vector_valid = vector_valid .AND. &
                           ASSOCIATED(qs_ot_env(ichannel)%rot_mat_mermin_g0) .AND. &
                           ASSOCIATED(qs_ot_env(ichannel)%rot_mat_mermin_g0_im)
            IF (vector_valid) THEN
               CALL dbcsr_copy(qs_ot_env(ichannel)%rot_mat_mermin_g0, &
                               qs_ot_env(ichannel)%rot_mat_gx)
               CALL dbcsr_copy(qs_ot_env(ichannel)%rot_mat_mermin_g0_im, &
                               qs_ot_env(ichannel)%rot_mat_gx_im)
               CALL qs_ot_transport_ref_rotation_pair( &
                  qs_ot_env(ichannel), qs_ot_env(ichannel)%rot_mat_mermin_g0, &
                  qs_ot_env(ichannel)%rot_mat_mermin_g0_im, gradient=.TRUE.)
            END IF
         END IF
         IF (qs_ot_env(1)%settings%do_ener) THEN
            vector_valid = vector_valid .AND. &
                           ASSOCIATED(qs_ot_env(ichannel)%ener_mermin_g0) .AND. &
                           ASSOCIATED(qs_ot_env(ichannel)%ener_mermin_x0)
            IF (vector_valid) THEN
               qs_ot_env(ichannel)%ener_mermin_g0(:) = qs_ot_env(ichannel)%ener_gx(:)
               qs_ot_env(ichannel)%ener_mermin_x0(:) = qs_ot_env(ichannel)%ener_x(:)
            END IF
         END IF
         frozen_valid = .FALSE.
         IF (vector_valid .AND. qs_ot_env(1)%settings%do_rotation) THEN
            CALL store_kpoint_mermin_frozen_hamiltonian(qs_ot_env(ichannel), frozen_valid)
         END IF
         qs_ot_env(ichannel)%mermin_frozen_h_valid = frozen_valid
         qs_ot_env(ichannel)%mermin_ls_vector_start_valid = vector_valid
         qs_ot_env(ichannel)%mermin_ls_component_slope_valid = .TRUE.
         matrix_slope = matrix_slope + qs_ot_env(ichannel)%mermin_ls_matrix_slope
         rotation_slope = rotation_slope + qs_ot_env(ichannel)%mermin_ls_rotation_slope
         energy_slope = energy_slope + qs_ot_env(ichannel)%mermin_ls_energy_slope
      END DO

      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, matrix_slope)
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, rotation_slope)
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, energy_slope)
      slope_sum = matrix_slope + rotation_slope + energy_slope
      IF (logger%iter_info%print_level >= debug_print_level .AND. logger%para_env%is_source()) THEN
         WRITE (cp_logger_get_default_unit_nr(logger), '(A,4(1X,ES16.8))') &
            " OT Mermin direction component slopes matrix rotation energy sum:", &
            matrix_slope, rotation_slope, energy_slope, slope_sum
      END IF

   END SUBROUTINE store_kpoint_mermin_direction_slopes

! **************************************************************************************************
!> \brief Store H at a line-search start in the fixed REF basis.
!> \param qs_ot_env local complex k-point OT channel
!> \param valid whether all required runtime matrices were available
! **************************************************************************************************
   SUBROUTINE store_kpoint_mermin_frozen_hamiltonian(qs_ot_env, valid)
      TYPE(qs_ot_type)                                   :: qs_ot_env
      LOGICAL, INTENT(OUT)                               :: valid

      TYPE(dbcsr_type)                                   :: product_im, product_re, tmp

      valid = ASSOCIATED(qs_ot_env%rot_mat_u) .AND. &
              ASSOCIATED(qs_ot_env%rot_mat_u_im) .AND. &
              ASSOCIATED(qs_ot_env%rot_mat_chc) .AND. &
              ASSOCIATED(qs_ot_env%rot_mat_chc_im) .AND. &
              ASSOCIATED(qs_ot_env%rot_mat_mermin_h0_ref) .AND. &
              ASSOCIATED(qs_ot_env%rot_mat_mermin_h0_ref_im)
      IF (.NOT. valid) RETURN

      CALL dbcsr_copy(product_re, qs_ot_env%rot_mat_chc, name="mermin_h0_product_re")
      CALL dbcsr_copy(product_im, qs_ot_env%rot_mat_chc, name="mermin_h0_product_im")
      CALL dbcsr_copy(tmp, qs_ot_env%rot_mat_chc, name="mermin_h0_tmp")

      ! product = U * (U^H H U).
      CALL dbcsr_multiply('N', 'N', 1.0_dp, qs_ot_env%rot_mat_u, &
                          qs_ot_env%rot_mat_chc, 0.0_dp, product_re)
      CALL dbcsr_multiply('N', 'N', -1.0_dp, qs_ot_env%rot_mat_u_im, &
                          qs_ot_env%rot_mat_chc_im, 1.0_dp, product_re)
      CALL dbcsr_multiply('N', 'N', 1.0_dp, qs_ot_env%rot_mat_u, &
                          qs_ot_env%rot_mat_chc_im, 0.0_dp, product_im)
      CALL dbcsr_multiply('N', 'N', 1.0_dp, qs_ot_env%rot_mat_u_im, &
                          qs_ot_env%rot_mat_chc, 1.0_dp, product_im)

      ! H_REF = product * U^H.
      CALL dbcsr_multiply('N', 'T', 1.0_dp, product_re, qs_ot_env%rot_mat_u, &
                          0.0_dp, qs_ot_env%rot_mat_mermin_h0_ref)
      CALL dbcsr_multiply('N', 'T', 1.0_dp, product_im, qs_ot_env%rot_mat_u_im, &
                          0.0_dp, tmp)
      CALL dbcsr_add(qs_ot_env%rot_mat_mermin_h0_ref, tmp, &
                     alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
      CALL dbcsr_multiply('N', 'T', 1.0_dp, product_im, qs_ot_env%rot_mat_u, &
                          0.0_dp, qs_ot_env%rot_mat_mermin_h0_ref_im)
      CALL dbcsr_multiply('N', 'T', 1.0_dp, product_re, qs_ot_env%rot_mat_u_im, &
                          0.0_dp, tmp)
      CALL dbcsr_add(qs_ot_env%rot_mat_mermin_h0_ref_im, tmp, &
                     alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)

      CALL dbcsr_release(product_im)
      CALL dbcsr_release(product_re)
      CALL dbcsr_release(tmp)

   END SUBROUTINE store_kpoint_mermin_frozen_hamiltonian

! **************************************************************************************************
!> \brief compare accepted finite Mermin curvature with endpoint component derivatives
!> \param qs_ot_env ...
!> \param para_env_inter_kp ...
! **************************************************************************************************
   SUBROUTINE audit_kpoint_mermin_accepted_secant(qs_ot_env, para_env_inter_kp)
      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: qs_ot_env
      TYPE(mp_para_env_type), OPTIONAL, POINTER          :: para_env_inter_kp

      INTEGER :: anchor_valid_i, closure_valid_i, ichannel, previous_closure_valid_i, &
         previous_valid_i, valid_i, vector_valid_i
      LOGICAL :: anchor_replaced, anchor_valid, block_secant_consistent, closure_valid, &
         energy_curvature_resolved, energy_secant_consistent, previous_closure_valid, &
         previous_valid, valid, vector_valid
      REAL(KIND=dp) :: alpha, auxiliary_curvature, block_curvature, current_energy, &
         current_matrix, current_rotation, current_total, energy_curvature, frozen_energy_sty, &
         frozen_energy_y2, frozen_rotation_sty, frozen_rotation_y2, initial_energy, &
         initial_matrix, initial_rotation, initial_total, matrix_curvature, matrix_fraction, &
         rotation_curvature, secant_alignment, secant_closure_energy_sty, &
         secant_closure_energy_y2, secant_closure_rotation_sty, secant_closure_rotation_y2, &
         secant_consistency, secant_energy_s2, secant_energy_sty, secant_energy_y2, &
         secant_orthogonal_fraction, secant_rotation_s2
      REAL(KIND=dp) :: secant_rotation_sty, secant_rotation_y2, secant_s2, secant_snorm, &
         secant_sty, secant_y2, secant_ynorm, tmp, total_curvature, vector_secant_curvature
      TYPE(cp_logger_type), POINTER                      :: logger

      IF (.NOT. qs_ot_env(1)%mermin_ls_sample_pending) RETURN
      logger => cp_get_default_logger()
      IF (.NOT. use_three_point_mermin_search(qs_ot_env)) RETURN

      alpha = qs_ot_env(1)%mermin_ls_next_pos
      valid = .TRUE.
      vector_valid = ABS(alpha) > 1.0E4_dp*EPSILON(1.0_dp)
      closure_valid = vector_valid
      initial_matrix = 0.0_dp
      initial_rotation = 0.0_dp
      initial_energy = 0.0_dp
      current_matrix = 0.0_dp
      current_rotation = 0.0_dp
      current_energy = 0.0_dp
      secant_rotation_s2 = 0.0_dp
      secant_rotation_y2 = 0.0_dp
      secant_rotation_sty = 0.0_dp
      secant_energy_s2 = 0.0_dp
      secant_energy_y2 = 0.0_dp
      secant_energy_sty = 0.0_dp
      frozen_rotation_y2 = 0.0_dp
      frozen_rotation_sty = 0.0_dp
      frozen_energy_y2 = 0.0_dp
      frozen_energy_sty = 0.0_dp
      secant_closure_rotation_y2 = 0.0_dp
      secant_closure_rotation_sty = 0.0_dp
      secant_closure_energy_y2 = 0.0_dp
      secant_closure_energy_sty = 0.0_dp

      anchor_valid = .TRUE.
      DO ichannel = 1, SIZE(qs_ot_env)
         IF (.NOT. qs_ot_env(ichannel)%has_complex_kpoint_state) CYCLE
         anchor_valid = anchor_valid .AND. &
                        qs_ot_env(ichannel)%mermin_secant_previous_ready .AND. &
                        qs_ot_env(ichannel)%mermin_closure_secant_previous_ready
         IF (qs_ot_env(1)%settings%do_rotation) THEN
            anchor_valid = anchor_valid .AND. &
                           ASSOCIATED(qs_ot_env(ichannel)%rot_mat_mermin_s_previous) .AND. &
                           ASSOCIATED(qs_ot_env(ichannel)%rot_mat_mermin_s_previous_im) .AND. &
                           ASSOCIATED(qs_ot_env(ichannel)%rot_mat_mermin_closure_y_previous) .AND. &
                           ASSOCIATED(qs_ot_env(ichannel)%rot_mat_mermin_closure_y_previous_im) .AND. &
                           ASSOCIATED(qs_ot_env(ichannel)%rot_mat_mermin_s_anchor) .AND. &
                           ASSOCIATED(qs_ot_env(ichannel)%rot_mat_mermin_s_anchor_im) .AND. &
                           ASSOCIATED(qs_ot_env(ichannel)%rot_mat_mermin_y_previous) .AND. &
                           ASSOCIATED(qs_ot_env(ichannel)%rot_mat_mermin_y_previous_im) .AND. &
                           ASSOCIATED(qs_ot_env(ichannel)%rot_mat_mermin_y_anchor) .AND. &
                           ASSOCIATED(qs_ot_env(ichannel)%rot_mat_mermin_y_anchor_im) .AND. &
                           ASSOCIATED(qs_ot_env(ichannel)%rot_mat_mermin_closure_y_anchor) .AND. &
                           ASSOCIATED(qs_ot_env(ichannel)%rot_mat_mermin_closure_y_anchor_im)
         END IF
         IF (qs_ot_env(1)%settings%do_ener) THEN
            anchor_valid = anchor_valid .AND. &
                           ASSOCIATED(qs_ot_env(ichannel)%ener_mermin_s_previous) .AND. &
                           ASSOCIATED(qs_ot_env(ichannel)%ener_mermin_y_previous) .AND. &
                           ASSOCIATED(qs_ot_env(ichannel)%ener_mermin_closure_y_previous) .AND. &
                           ASSOCIATED(qs_ot_env(ichannel)%ener_mermin_s_anchor) .AND. &
                           ASSOCIATED(qs_ot_env(ichannel)%ener_mermin_y_anchor) .AND. &
                           ASSOCIATED(qs_ot_env(ichannel)%ener_mermin_closure_y_anchor)
         END IF
      END DO
      anchor_valid_i = MERGE(1, 0, anchor_valid)
      CALL ot_mini_min_kpoint_integer(qs_ot_env, para_env_inter_kp, anchor_valid_i)
      anchor_valid = anchor_valid_i > 0
      ! Current and previous follow the two newest accepted observations.  Keep rolling the older
      ! anchor until the causal response audit validates it; after that, preserve the demonstrated
      ! response direction across collapsed tail line-search steps until structural invalidation.
      anchor_replaced = anchor_valid .AND. &
                        .NOT. qs_ot_env(1)%mermin_closure_secant_anchor_validated
      IF (anchor_replaced) THEN
         DO ichannel = 1, SIZE(qs_ot_env)
            IF (.NOT. qs_ot_env(ichannel)%has_complex_kpoint_state) CYCLE
            IF (qs_ot_env(1)%settings%do_rotation) THEN
               CALL dbcsr_copy(qs_ot_env(ichannel)%rot_mat_mermin_s_anchor, &
                               qs_ot_env(ichannel)%rot_mat_mermin_s_previous)
               CALL dbcsr_copy(qs_ot_env(ichannel)%rot_mat_mermin_s_anchor_im, &
                               qs_ot_env(ichannel)%rot_mat_mermin_s_previous_im)
               CALL dbcsr_copy(qs_ot_env(ichannel)%rot_mat_mermin_y_anchor, &
                               qs_ot_env(ichannel)%rot_mat_mermin_y_previous)
               CALL dbcsr_copy(qs_ot_env(ichannel)%rot_mat_mermin_y_anchor_im, &
                               qs_ot_env(ichannel)%rot_mat_mermin_y_previous_im)
               CALL dbcsr_copy(qs_ot_env(ichannel)%rot_mat_mermin_closure_y_anchor, &
                               qs_ot_env(ichannel)%rot_mat_mermin_closure_y_previous)
               CALL dbcsr_copy(qs_ot_env(ichannel)%rot_mat_mermin_closure_y_anchor_im, &
                               qs_ot_env(ichannel)%rot_mat_mermin_closure_y_previous_im)
            END IF
            IF (qs_ot_env(1)%settings%do_ener) THEN
               qs_ot_env(ichannel)%ener_mermin_s_anchor(:) = &
                  qs_ot_env(ichannel)%ener_mermin_s_previous(:)
               qs_ot_env(ichannel)%ener_mermin_y_anchor(:) = &
                  qs_ot_env(ichannel)%ener_mermin_y_previous(:)
               qs_ot_env(ichannel)%ener_mermin_closure_y_anchor(:) = &
                  qs_ot_env(ichannel)%ener_mermin_closure_y_previous(:)
            END IF
            qs_ot_env(ichannel)%mermin_secant_work_ratio_anchor = &
               qs_ot_env(ichannel)%mermin_secant_work_ratio_previous
         END DO
      END IF

      previous_valid = .TRUE.
      DO ichannel = 1, SIZE(qs_ot_env)
         IF (.NOT. qs_ot_env(ichannel)%has_complex_kpoint_state) CYCLE
         previous_valid = previous_valid .AND. qs_ot_env(ichannel)%mermin_secant_ready
         IF (qs_ot_env(1)%settings%do_rotation) THEN
            previous_valid = previous_valid .AND. &
                             ASSOCIATED(qs_ot_env(ichannel)%rot_mat_mermin_s) .AND. &
                             ASSOCIATED(qs_ot_env(ichannel)%rot_mat_mermin_s_im) .AND. &
                             ASSOCIATED(qs_ot_env(ichannel)%rot_mat_mermin_y) .AND. &
                             ASSOCIATED(qs_ot_env(ichannel)%rot_mat_mermin_y_im) .AND. &
                             ASSOCIATED(qs_ot_env(ichannel)%rot_mat_mermin_s_previous) .AND. &
                             ASSOCIATED(qs_ot_env(ichannel)%rot_mat_mermin_s_previous_im) .AND. &
                             ASSOCIATED(qs_ot_env(ichannel)%rot_mat_mermin_y_previous) .AND. &
                             ASSOCIATED(qs_ot_env(ichannel)%rot_mat_mermin_y_previous_im)
         END IF
         IF (qs_ot_env(1)%settings%do_ener) THEN
            previous_valid = previous_valid .AND. &
                             ASSOCIATED(qs_ot_env(ichannel)%ener_mermin_s) .AND. &
                             ASSOCIATED(qs_ot_env(ichannel)%ener_mermin_y) .AND. &
                             ASSOCIATED(qs_ot_env(ichannel)%ener_mermin_s_previous) .AND. &
                             ASSOCIATED(qs_ot_env(ichannel)%ener_mermin_y_previous)
         END IF
      END DO
      previous_valid_i = MERGE(1, 0, previous_valid)
      CALL ot_mini_min_kpoint_integer(qs_ot_env, para_env_inter_kp, previous_valid_i)
      previous_valid = previous_valid_i > 0
      IF (previous_valid) THEN
         DO ichannel = 1, SIZE(qs_ot_env)
            IF (.NOT. qs_ot_env(ichannel)%has_complex_kpoint_state) CYCLE
            IF (qs_ot_env(1)%settings%do_rotation) THEN
               CALL dbcsr_copy(qs_ot_env(ichannel)%rot_mat_mermin_s_previous, &
                               qs_ot_env(ichannel)%rot_mat_mermin_s)
               CALL dbcsr_copy(qs_ot_env(ichannel)%rot_mat_mermin_s_previous_im, &
                               qs_ot_env(ichannel)%rot_mat_mermin_s_im)
               CALL dbcsr_copy(qs_ot_env(ichannel)%rot_mat_mermin_y_previous, &
                               qs_ot_env(ichannel)%rot_mat_mermin_y)
               CALL dbcsr_copy(qs_ot_env(ichannel)%rot_mat_mermin_y_previous_im, &
                               qs_ot_env(ichannel)%rot_mat_mermin_y_im)
            END IF
            IF (qs_ot_env(1)%settings%do_ener) THEN
               qs_ot_env(ichannel)%ener_mermin_s_previous(:) = &
                  qs_ot_env(ichannel)%ener_mermin_s(:)
               qs_ot_env(ichannel)%ener_mermin_y_previous(:) = &
                  qs_ot_env(ichannel)%ener_mermin_y(:)
            END IF
            qs_ot_env(ichannel)%mermin_secant_work_ratio_previous = &
               qs_ot_env(ichannel)%mermin_secant_work_ratio
         END DO
      END IF

      previous_closure_valid = previous_valid
      DO ichannel = 1, SIZE(qs_ot_env)
         IF (.NOT. qs_ot_env(ichannel)%has_complex_kpoint_state) CYCLE
         previous_closure_valid = previous_closure_valid .AND. &
                                  qs_ot_env(ichannel)%mermin_closure_secant_ready
         IF (qs_ot_env(1)%settings%do_rotation) THEN
            previous_closure_valid = previous_closure_valid .AND. &
                                     ASSOCIATED(qs_ot_env(ichannel)%rot_mat_mermin_closure_y) .AND. &
                                     ASSOCIATED(qs_ot_env(ichannel)%rot_mat_mermin_closure_y_im) .AND. &
                                     ASSOCIATED(qs_ot_env(ichannel)%rot_mat_mermin_closure_y_previous) .AND. &
                                     ASSOCIATED(qs_ot_env(ichannel)%rot_mat_mermin_closure_y_previous_im)
         END IF
         IF (qs_ot_env(1)%settings%do_ener) THEN
            previous_closure_valid = previous_closure_valid .AND. &
                                     ASSOCIATED(qs_ot_env(ichannel)%ener_mermin_closure_y) .AND. &
                                     ASSOCIATED(qs_ot_env(ichannel)%ener_mermin_closure_y_previous)
         END IF
      END DO
      previous_closure_valid_i = MERGE(1, 0, previous_closure_valid)
      CALL ot_mini_min_kpoint_integer(qs_ot_env, para_env_inter_kp, previous_closure_valid_i)
      previous_closure_valid = previous_closure_valid_i > 0
      IF (previous_closure_valid) THEN
         DO ichannel = 1, SIZE(qs_ot_env)
            IF (.NOT. qs_ot_env(ichannel)%has_complex_kpoint_state) CYCLE
            IF (qs_ot_env(1)%settings%do_rotation) THEN
               CALL dbcsr_copy(qs_ot_env(ichannel)%rot_mat_mermin_closure_y_previous, &
                               qs_ot_env(ichannel)%rot_mat_mermin_closure_y)
               CALL dbcsr_copy(qs_ot_env(ichannel)%rot_mat_mermin_closure_y_previous_im, &
                               qs_ot_env(ichannel)%rot_mat_mermin_closure_y_im)
            END IF
            IF (qs_ot_env(1)%settings%do_ener) THEN
               qs_ot_env(ichannel)%ener_mermin_closure_y_previous(:) = &
                  qs_ot_env(ichannel)%ener_mermin_closure_y(:)
            END IF
         END DO
      END IF

      DO ichannel = 1, SIZE(qs_ot_env)
         IF (.NOT. qs_ot_env(ichannel)%has_complex_kpoint_state) CYCLE
         valid = valid .AND. qs_ot_env(ichannel)%mermin_ls_component_slope_valid
         vector_valid = vector_valid .AND. qs_ot_env(ichannel)%mermin_ls_vector_start_valid
         closure_valid = closure_valid .AND. qs_ot_env(ichannel)%mermin_frozen_h_valid .AND. &
                         qs_ot_env(ichannel)%mermin_frozen_endpoint_valid
         initial_matrix = initial_matrix + qs_ot_env(ichannel)%mermin_ls_matrix_slope
         initial_rotation = initial_rotation + qs_ot_env(ichannel)%mermin_ls_rotation_slope
         initial_energy = initial_energy + qs_ot_env(ichannel)%mermin_ls_energy_slope
         IF (.NOT. qs_ot_env(ichannel)%mermin_ls_component_slope_valid) CYCLE
         IF (qs_ot_env(1)%settings%ks) THEN
            IF (.NOT. ASSOCIATED(qs_ot_env(ichannel)%matrix_gx) .OR. &
                .NOT. ASSOCIATED(qs_ot_env(ichannel)%matrix_dx) .OR. &
                .NOT. ASSOCIATED(qs_ot_env(ichannel)%matrix_gx_im) .OR. &
                .NOT. ASSOCIATED(qs_ot_env(ichannel)%matrix_dx_im)) THEN
               valid = .FALSE.
               CYCLE
            END IF
            CALL dbcsr_dot(qs_ot_env(ichannel)%matrix_gx, qs_ot_env(ichannel)%matrix_dx, tmp)
            current_matrix = current_matrix + tmp
            CALL dbcsr_dot(qs_ot_env(ichannel)%matrix_gx_im, qs_ot_env(ichannel)%matrix_dx_im, tmp)
            current_matrix = current_matrix + tmp
         END IF
         IF (qs_ot_env(1)%settings%do_rotation) THEN
            IF (.NOT. ASSOCIATED(qs_ot_env(ichannel)%rot_mat_gx) .OR. &
                .NOT. ASSOCIATED(qs_ot_env(ichannel)%rot_mat_dx) .OR. &
                .NOT. ASSOCIATED(qs_ot_env(ichannel)%rot_mat_gx_im) .OR. &
                .NOT. ASSOCIATED(qs_ot_env(ichannel)%rot_mat_dx_im)) THEN
               valid = .FALSE.
               CYCLE
            END IF
            CALL dbcsr_dot(qs_ot_env(ichannel)%rot_mat_gx, qs_ot_env(ichannel)%rot_mat_dx, tmp)
            current_rotation = current_rotation + 0.5_dp*tmp
            CALL dbcsr_dot(qs_ot_env(ichannel)%rot_mat_gx_im, qs_ot_env(ichannel)%rot_mat_dx_im, tmp)
            current_rotation = current_rotation + 0.5_dp*tmp
         END IF
         IF (qs_ot_env(1)%settings%do_ener) THEN
            IF (.NOT. ASSOCIATED(qs_ot_env(ichannel)%ener_gx) .OR. &
                .NOT. ASSOCIATED(qs_ot_env(ichannel)%ener_dx)) THEN
               valid = .FALSE.
               CYCLE
            END IF
            current_energy = current_energy + &
                             qs_ot_env(ichannel)%ener_x_line_search_tangent* &
                             DOT_PRODUCT(qs_ot_env(ichannel)%ener_gx, qs_ot_env(ichannel)%ener_dx)
         END IF
         IF (.NOT. vector_valid) CYCLE
         IF (qs_ot_env(1)%settings%do_rotation) THEN
            vector_valid = ASSOCIATED(qs_ot_env(ichannel)%rot_mat_mermin_s) .AND. &
                           ASSOCIATED(qs_ot_env(ichannel)%rot_mat_mermin_s_im) .AND. &
                           ASSOCIATED(qs_ot_env(ichannel)%rot_mat_mermin_y) .AND. &
                           ASSOCIATED(qs_ot_env(ichannel)%rot_mat_mermin_y_im)
            IF (.NOT. vector_valid) CYCLE
            CALL dbcsr_copy(qs_ot_env(ichannel)%rot_mat_mermin_s, &
                            qs_ot_env(ichannel)%rot_mat_dx)
            CALL dbcsr_scale(qs_ot_env(ichannel)%rot_mat_mermin_s, alpha)
            CALL dbcsr_copy(qs_ot_env(ichannel)%rot_mat_mermin_s_im, &
                            qs_ot_env(ichannel)%rot_mat_dx_im)
            CALL dbcsr_scale(qs_ot_env(ichannel)%rot_mat_mermin_s_im, alpha)
            CALL qs_ot_transport_ref_rotation_pair( &
               qs_ot_env(ichannel), qs_ot_env(ichannel)%rot_mat_mermin_s, &
               qs_ot_env(ichannel)%rot_mat_mermin_s_im, gradient=.FALSE.)
            CALL dbcsr_copy(qs_ot_env(ichannel)%rot_mat_mermin_y, &
                            qs_ot_env(ichannel)%rot_mat_gx)
            CALL dbcsr_copy(qs_ot_env(ichannel)%rot_mat_mermin_y_im, &
                            qs_ot_env(ichannel)%rot_mat_gx_im)
            CALL qs_ot_transport_ref_rotation_pair( &
               qs_ot_env(ichannel), qs_ot_env(ichannel)%rot_mat_mermin_y, &
               qs_ot_env(ichannel)%rot_mat_mermin_y_im, gradient=.TRUE.)
            CALL dbcsr_add(qs_ot_env(ichannel)%rot_mat_mermin_y, &
                           qs_ot_env(ichannel)%rot_mat_mermin_g0, &
                           alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)
            CALL dbcsr_add(qs_ot_env(ichannel)%rot_mat_mermin_y_im, &
                           qs_ot_env(ichannel)%rot_mat_mermin_g0_im, &
                           alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)
            closure_valid = closure_valid .AND. &
                            ASSOCIATED(qs_ot_env(ichannel)%rot_mat_mermin_frozen_g1) .AND. &
                            ASSOCIATED(qs_ot_env(ichannel)%rot_mat_mermin_frozen_g1_im) .AND. &
                            ASSOCIATED(qs_ot_env(ichannel)%rot_mat_mermin_closure_y) .AND. &
                            ASSOCIATED(qs_ot_env(ichannel)%rot_mat_mermin_closure_y_im)
            IF (closure_valid) THEN
               CALL dbcsr_copy(qs_ot_env(ichannel)%rot_mat_mermin_closure_y, &
                               qs_ot_env(ichannel)%rot_mat_mermin_frozen_g1)
               CALL dbcsr_copy(qs_ot_env(ichannel)%rot_mat_mermin_closure_y_im, &
                               qs_ot_env(ichannel)%rot_mat_mermin_frozen_g1_im)
               CALL qs_ot_transport_ref_rotation_pair( &
                  qs_ot_env(ichannel), qs_ot_env(ichannel)%rot_mat_mermin_closure_y, &
                  qs_ot_env(ichannel)%rot_mat_mermin_closure_y_im, gradient=.TRUE.)
               CALL dbcsr_add(qs_ot_env(ichannel)%rot_mat_mermin_closure_y, &
                              qs_ot_env(ichannel)%rot_mat_mermin_g0, &
                              alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)
               CALL dbcsr_add(qs_ot_env(ichannel)%rot_mat_mermin_closure_y_im, &
                              qs_ot_env(ichannel)%rot_mat_mermin_g0_im, &
                              alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)
               CALL dbcsr_dot(qs_ot_env(ichannel)%rot_mat_mermin_closure_y, &
                              qs_ot_env(ichannel)%rot_mat_mermin_closure_y, tmp)
               frozen_rotation_y2 = frozen_rotation_y2 + 0.5_dp*tmp
               CALL dbcsr_dot(qs_ot_env(ichannel)%rot_mat_mermin_closure_y_im, &
                              qs_ot_env(ichannel)%rot_mat_mermin_closure_y_im, tmp)
               frozen_rotation_y2 = frozen_rotation_y2 + 0.5_dp*tmp
               CALL dbcsr_dot(qs_ot_env(ichannel)%rot_mat_mermin_s, &
                              qs_ot_env(ichannel)%rot_mat_mermin_closure_y, tmp)
               frozen_rotation_sty = frozen_rotation_sty + 0.5_dp*tmp
               CALL dbcsr_dot(qs_ot_env(ichannel)%rot_mat_mermin_s_im, &
                              qs_ot_env(ichannel)%rot_mat_mermin_closure_y_im, tmp)
               frozen_rotation_sty = frozen_rotation_sty + 0.5_dp*tmp
               CALL dbcsr_scale(qs_ot_env(ichannel)%rot_mat_mermin_closure_y, -1.0_dp)
               CALL dbcsr_scale(qs_ot_env(ichannel)%rot_mat_mermin_closure_y_im, -1.0_dp)
               CALL dbcsr_add(qs_ot_env(ichannel)%rot_mat_mermin_closure_y, &
                              qs_ot_env(ichannel)%rot_mat_mermin_y, &
                              alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
               CALL dbcsr_add(qs_ot_env(ichannel)%rot_mat_mermin_closure_y_im, &
                              qs_ot_env(ichannel)%rot_mat_mermin_y_im, &
                              alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
               CALL dbcsr_dot(qs_ot_env(ichannel)%rot_mat_mermin_closure_y, &
                              qs_ot_env(ichannel)%rot_mat_mermin_closure_y, tmp)
               secant_closure_rotation_y2 = secant_closure_rotation_y2 + 0.5_dp*tmp
               CALL dbcsr_dot(qs_ot_env(ichannel)%rot_mat_mermin_closure_y_im, &
                              qs_ot_env(ichannel)%rot_mat_mermin_closure_y_im, tmp)
               secant_closure_rotation_y2 = secant_closure_rotation_y2 + 0.5_dp*tmp
               CALL dbcsr_dot(qs_ot_env(ichannel)%rot_mat_mermin_s, &
                              qs_ot_env(ichannel)%rot_mat_mermin_closure_y, tmp)
               secant_closure_rotation_sty = secant_closure_rotation_sty + 0.5_dp*tmp
               CALL dbcsr_dot(qs_ot_env(ichannel)%rot_mat_mermin_s_im, &
                              qs_ot_env(ichannel)%rot_mat_mermin_closure_y_im, tmp)
               secant_closure_rotation_sty = secant_closure_rotation_sty + 0.5_dp*tmp
            END IF
            CALL dbcsr_dot(qs_ot_env(ichannel)%rot_mat_mermin_s, &
                           qs_ot_env(ichannel)%rot_mat_mermin_s, tmp)
            secant_rotation_s2 = secant_rotation_s2 + 0.5_dp*tmp
            CALL dbcsr_dot(qs_ot_env(ichannel)%rot_mat_mermin_s_im, &
                           qs_ot_env(ichannel)%rot_mat_mermin_s_im, tmp)
            secant_rotation_s2 = secant_rotation_s2 + 0.5_dp*tmp
            CALL dbcsr_dot(qs_ot_env(ichannel)%rot_mat_mermin_y, &
                           qs_ot_env(ichannel)%rot_mat_mermin_y, tmp)
            secant_rotation_y2 = secant_rotation_y2 + 0.5_dp*tmp
            CALL dbcsr_dot(qs_ot_env(ichannel)%rot_mat_mermin_y_im, &
                           qs_ot_env(ichannel)%rot_mat_mermin_y_im, tmp)
            secant_rotation_y2 = secant_rotation_y2 + 0.5_dp*tmp
            CALL dbcsr_dot(qs_ot_env(ichannel)%rot_mat_mermin_s, &
                           qs_ot_env(ichannel)%rot_mat_mermin_y, tmp)
            secant_rotation_sty = secant_rotation_sty + 0.5_dp*tmp
            CALL dbcsr_dot(qs_ot_env(ichannel)%rot_mat_mermin_s_im, &
                           qs_ot_env(ichannel)%rot_mat_mermin_y_im, tmp)
            secant_rotation_sty = secant_rotation_sty + 0.5_dp*tmp
         END IF
         IF (qs_ot_env(1)%settings%do_ener) THEN
            vector_valid = ASSOCIATED(qs_ot_env(ichannel)%ener_mermin_x0) .AND. &
                           ASSOCIATED(qs_ot_env(ichannel)%ener_mermin_s) .AND. &
                           ASSOCIATED(qs_ot_env(ichannel)%ener_mermin_y)
            IF (.NOT. vector_valid) CYCLE
            qs_ot_env(ichannel)%ener_mermin_s(:) = qs_ot_env(ichannel)%ener_x(:) - &
                                                   qs_ot_env(ichannel)%ener_mermin_x0(:)
            qs_ot_env(ichannel)%ener_mermin_y(:) = qs_ot_env(ichannel)%ener_gx(:) - &
                                                   qs_ot_env(ichannel)%ener_mermin_g0(:)
            closure_valid = closure_valid .AND. &
                            ASSOCIATED(qs_ot_env(ichannel)%ener_mermin_frozen_g1) .AND. &
                            ASSOCIATED(qs_ot_env(ichannel)%ener_mermin_closure_y)
            IF (closure_valid) THEN
               qs_ot_env(ichannel)%ener_mermin_closure_y(:) = &
                  qs_ot_env(ichannel)%ener_mermin_frozen_g1(:) - &
                  qs_ot_env(ichannel)%ener_mermin_g0(:)
               frozen_energy_y2 = frozen_energy_y2 + &
                                  DOT_PRODUCT(qs_ot_env(ichannel)%ener_mermin_closure_y, &
                                              qs_ot_env(ichannel)%ener_mermin_closure_y)
               frozen_energy_sty = frozen_energy_sty + &
                                   DOT_PRODUCT(qs_ot_env(ichannel)%ener_mermin_s, &
                                               qs_ot_env(ichannel)%ener_mermin_closure_y)
               qs_ot_env(ichannel)%ener_mermin_closure_y(:) = &
                  qs_ot_env(ichannel)%ener_mermin_y(:) - &
                  qs_ot_env(ichannel)%ener_mermin_closure_y(:)
               secant_closure_energy_y2 = secant_closure_energy_y2 + &
                                          DOT_PRODUCT(qs_ot_env(ichannel)%ener_mermin_closure_y, &
                                                      qs_ot_env(ichannel)%ener_mermin_closure_y)
               secant_closure_energy_sty = secant_closure_energy_sty + &
                                           DOT_PRODUCT(qs_ot_env(ichannel)%ener_mermin_s, &
                                                       qs_ot_env(ichannel)%ener_mermin_closure_y)
            END IF
            secant_energy_s2 = secant_energy_s2 + &
                               DOT_PRODUCT(qs_ot_env(ichannel)%ener_mermin_s, &
                                           qs_ot_env(ichannel)%ener_mermin_s)
            secant_energy_y2 = secant_energy_y2 + &
                               DOT_PRODUCT(qs_ot_env(ichannel)%ener_mermin_y, &
                                           qs_ot_env(ichannel)%ener_mermin_y)
            secant_energy_sty = secant_energy_sty + &
                                DOT_PRODUCT(qs_ot_env(ichannel)%ener_mermin_s, &
                                            qs_ot_env(ichannel)%ener_mermin_y)
         END IF
      END DO
      valid_i = MERGE(1, 0, valid)
      CALL ot_mini_min_kpoint_integer(qs_ot_env, para_env_inter_kp, valid_i)
      valid = valid_i > 0
      vector_valid_i = MERGE(1, 0, vector_valid)
      CALL ot_mini_min_kpoint_integer(qs_ot_env, para_env_inter_kp, vector_valid_i)
      vector_valid = vector_valid_i > 0
      closure_valid_i = MERGE(1, 0, closure_valid)
      CALL ot_mini_min_kpoint_integer(qs_ot_env, para_env_inter_kp, closure_valid_i)
      closure_valid = closure_valid_i > 0
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, initial_matrix)
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, initial_rotation)
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, initial_energy)
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, current_matrix)
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, current_rotation)
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, current_energy)
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, secant_rotation_s2)
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, secant_rotation_y2)
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, secant_rotation_sty)
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, secant_energy_s2)
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, secant_energy_y2)
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, secant_energy_sty)
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, frozen_rotation_y2)
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, frozen_rotation_sty)
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, frozen_energy_y2)
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, frozen_energy_sty)
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, secant_closure_rotation_y2)
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, secant_closure_rotation_sty)
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, secant_closure_energy_y2)
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, secant_closure_energy_sty)

      valid = valid .AND. ABS(alpha) > 1.0E4_dp*EPSILON(1.0_dp)
      initial_total = initial_matrix + initial_rotation + initial_energy
      current_total = current_matrix + current_rotation + current_energy
      matrix_curvature = 0.0_dp
      rotation_curvature = 0.0_dp
      auxiliary_curvature = 0.0_dp
      total_curvature = 0.0_dp
      energy_curvature = 0.0_dp
      matrix_fraction = 0.0_dp
      secant_consistency = 0.0_dp
      block_curvature = 0.0_dp
      vector_secant_curvature = 0.0_dp
      block_secant_consistent = .FALSE.
      energy_curvature_resolved = .FALSE.
      energy_secant_consistent = .TRUE.
      secant_s2 = secant_rotation_s2 + secant_energy_s2
      secant_y2 = secant_rotation_y2 + secant_energy_y2
      secant_sty = secant_rotation_sty + secant_energy_sty
      IF (.NOT. vector_valid) THEN
         secant_s2 = 0.0_dp
         secant_y2 = 0.0_dp
         secant_sty = 0.0_dp
      END IF
      CALL ot_kpoint_secant_geometry(secant_s2, secant_y2, secant_sty, secant_snorm, &
                                     secant_ynorm, secant_alignment, secant_orthogonal_fraction)
      IF (valid) THEN
         matrix_curvature = (current_matrix - initial_matrix)/alpha
         rotation_curvature = (current_rotation - initial_rotation)/alpha
         auxiliary_curvature = (current_energy - initial_energy)/alpha
         total_curvature = (current_total - initial_total)/alpha
         energy_curvature = 2.0_dp* &
                            (qs_ot_env(1)%etotal - qs_ot_env(1)%mermin_ls_reference_energy - &
                             alpha*initial_total)/(alpha*alpha)
         IF (total_curvature > EPSILON(1.0_dp)) THEN
            matrix_fraction = matrix_curvature/total_curvature
         END IF
         secant_consistency = ABS(total_curvature - energy_curvature)/ &
                              MAX(ABS(total_curvature), ABS(energy_curvature), EPSILON(1.0_dp))
         block_curvature = rotation_curvature + auxiliary_curvature
         vector_secant_curvature = secant_sty/(alpha*alpha)
         block_secant_consistent = &
            ot_kpoint_mermin_secant_consistent(block_curvature, vector_secant_curvature)
         energy_curvature_resolved = ot_kpoint_mermin_energy_curvature_resolved( &
                                     alpha, qs_ot_env(1)%mermin_ls_reference_energy, &
                                     qs_ot_env(1)%etotal, initial_total, energy_curvature, &
                                     total_curvature)
         IF (energy_curvature_resolved) THEN
            energy_secant_consistent = &
               ot_kpoint_mermin_secant_consistent(total_curvature, energy_curvature)
         END IF
      END IF

      ! The stored vector secant contains only the complex rotation and auxiliary-energy
      ! coordinates.  Validate it against the endpoint derivative curvature of that same block;
      ! the matrix contribution can dominate the total line-search curvature without being part
      ! of the response fit.  Retain the finite Mermin-energy veto only while its O(alpha**2)
      ! remainder is numerically resolved.  Accepted-drop calibration still validates every
      ! resulting response candidate against the actual Mermin line search before promotion.
      vector_valid = vector_valid .AND. valid .AND. block_secant_consistent .AND. &
                     energy_secant_consistent
      closure_valid = closure_valid .AND. vector_valid
      DO ichannel = 1, SIZE(qs_ot_env)
         IF (.NOT. qs_ot_env(ichannel)%has_complex_kpoint_state) CYCLE
         ! The current/previous/anchor buffers form a history of valid accepted secants, not
         ! of every line-search endpoint.  A zero-position or unresolved endpoint may overwrite
         ! the current work buffers, but it must not erase an older valid response sample.
         IF (anchor_replaced) THEN
            qs_ot_env(ichannel)%mermin_closure_secant_anchor_flavor = &
               qs_ot_env(ichannel)%mermin_closure_secant_previous_flavor
         END IF
         IF (previous_closure_valid) THEN
            qs_ot_env(ichannel)%mermin_closure_secant_previous_flavor = &
               qs_ot_env(ichannel)%mermin_closure_secant_flavor
         END IF
         qs_ot_env(ichannel)%mermin_closure_secant_flavor = &
            MERGE(qs_ot_env(1)%schur_model_flavor, 0, closure_valid)
         qs_ot_env(ichannel)%mermin_secant_previous_ready = &
            qs_ot_env(ichannel)%mermin_secant_previous_ready .OR. previous_valid
         qs_ot_env(ichannel)%mermin_secant_ready = vector_valid
         qs_ot_env(ichannel)%mermin_closure_secant_previous_ready = &
            qs_ot_env(ichannel)%mermin_closure_secant_previous_ready .OR. previous_closure_valid
         qs_ot_env(ichannel)%mermin_closure_secant_anchor_ready = &
            qs_ot_env(ichannel)%mermin_closure_secant_anchor_ready .OR. anchor_valid
         qs_ot_env(ichannel)%mermin_closure_secant_ready = closure_valid
         qs_ot_env(ichannel)%mermin_frozen_h_valid = .FALSE.
         qs_ot_env(ichannel)%mermin_frozen_endpoint_valid = .FALSE.
         qs_ot_env(ichannel)%mermin_direction_curvature_valid = vector_valid .AND. &
                                                                total_curvature > EPSILON(MAX(1.0_dp, ABS(total_curvature)))
         qs_ot_env(ichannel)%mermin_direction_curvature = &
            MERGE(total_curvature, 0.0_dp, &
                  qs_ot_env(ichannel)%mermin_direction_curvature_valid)
         IF (.NOT. vector_valid) CYCLE
         qs_ot_env(ichannel)%mermin_secant_sty = secant_sty
         qs_ot_env(ichannel)%mermin_secant_work_ratio = energy_curvature/total_curvature
      END DO

      IF (logger%iter_info%print_level >= debug_print_level .AND. logger%para_env%is_source()) THEN
         WRITE (cp_logger_get_default_unit_nr(logger), '(A,1X,L1,8(1X,ES16.8))') &
            " OT Mermin accepted secant slopes valid alpha imat irot iener cmat crot cener current:", &
            valid, alpha, initial_matrix, initial_rotation, initial_energy, &
            current_matrix, current_rotation, current_energy, qs_ot_env(1)%gradient
         WRITE (cp_logger_get_default_unit_nr(logger), '(A,7(1X,ES16.8))') &
            " OT Mermin accepted secant curvature matrix rotation energy gradient Mermin fraction consistency:", &
            matrix_curvature, rotation_curvature, auxiliary_curvature, total_curvature, &
            energy_curvature, matrix_fraction, secant_consistency
         WRITE (cp_logger_get_default_unit_nr(logger), '(A,1X,L1,10(1X,ES16.8))') &
            " OT Mermin accepted vector secant valid srot sener yrot yener styrot styener snorm ynorm align orth:", &
            vector_valid, secant_rotation_s2, secant_energy_s2, secant_rotation_y2, &
            secant_energy_y2, secant_rotation_sty, secant_energy_sty, secant_snorm, &
            secant_ynorm, secant_alignment, secant_orthogonal_fraction
         WRITE (cp_logger_get_default_unit_nr(logger), '(A,4(1X,L1),3(1X,ES16.8))') &
            " OT Mermin accepted secant validation block energy_resolved energy_ok vector "// &
            "block_curvature vector_curvature finite_curvature:", &
            block_secant_consistent, energy_curvature_resolved, energy_secant_consistent, &
            vector_valid, block_curvature, vector_secant_curvature, energy_curvature
         WRITE (cp_logger_get_default_unit_nr(logger), '(A,1X,L1,8(1X,ES16.8))') &
            " OT Mermin frozen-H closure valid yfrot yfener ycrot ycener sfrot sfener scrot scener:", &
            closure_valid, frozen_rotation_y2, frozen_energy_y2, &
            secant_closure_rotation_y2, secant_closure_energy_y2, &
            frozen_rotation_sty, frozen_energy_sty, &
            secant_closure_rotation_sty, secant_closure_energy_sty
      END IF

   END SUBROUTINE audit_kpoint_mermin_accepted_secant

! **************************************************************************************************
!> \brief test whether an accepted derivative secant reproduces its finite Mermin curvature
!> \param gradient_curvature curvature from endpoint directional derivatives
!> \param energy_curvature curvature from the accepted finite Mermin-energy change
!> \return true when the sample is resolved and locally quadratic
! **************************************************************************************************
   PURE LOGICAL FUNCTION ot_kpoint_mermin_secant_consistent(gradient_curvature, energy_curvature) &
      RESULT(consistent)
      REAL(KIND=dp), INTENT(IN)                          :: gradient_curvature, energy_curvature

      REAL(KIND=dp)                                      :: curvature_scale, relative_mismatch

      curvature_scale = MAX(ABS(gradient_curvature), ABS(energy_curvature))
      relative_mismatch = ABS(gradient_curvature - energy_curvature)/ &
                          MAX(curvature_scale, EPSILON(1.0_dp))
      consistent = curvature_scale > 1.0E4_dp*EPSILON(1.0_dp) .AND. &
                   relative_mismatch <= ot_cg_kp_mermin_secant_consistency

   END FUNCTION ot_kpoint_mermin_secant_consistent

! **************************************************************************************************
!> \brief test whether a finite Mermin-energy curvature is above subtraction noise
!> \param alpha accepted line-search position
!> \param reference_energy Mermin energy at the line-search origin
!> \param endpoint_energy Mermin energy at the accepted endpoint
!> \param reference_slope directional derivative at the line-search origin
!> \param energy_curvature curvature reconstructed from the finite energy difference
!> \param gradient_curvature curvature reconstructed from endpoint derivatives
!> \return true when the O(alpha**2) energy remainder is numerically resolved
! **************************************************************************************************
   PURE LOGICAL FUNCTION ot_kpoint_mermin_energy_curvature_resolved( &
      alpha, reference_energy, endpoint_energy, reference_slope, energy_curvature, &
      gradient_curvature) RESULT(resolved)
      REAL(KIND=dp), INTENT(IN)                          :: alpha, reference_energy, &
                                                            endpoint_energy, reference_slope, &
                                                            energy_curvature, gradient_curvature

      REAL(KIND=dp)                                      :: derivative_remainder, energy_scale, &
                                                            quadratic_remainder, resolution

      energy_scale = MAX(1.0_dp, ABS(reference_energy), ABS(endpoint_energy), &
                         ABS(alpha*reference_slope))
      quadratic_remainder = 0.5_dp*ABS(alpha*alpha*energy_curvature)
      derivative_remainder = 0.5_dp*ABS(alpha*alpha*gradient_curvature)
      resolution = 64.0_dp*EPSILON(1.0_dp)*energy_scale
      ! Subtraction noise can look like a large curvature after division by alpha**2.
      ! Only use the finite-energy veto when the independent endpoint derivatives also
      ! predict a comfortably resolved quadratic remainder.
      resolved = quadratic_remainder > resolution .AND. &
                 derivative_remainder > 4.0_dp*resolution .AND. &
                 quadratic_remainder == quadratic_remainder .AND. &
                 derivative_remainder == derivative_remainder

   END FUNCTION ot_kpoint_mermin_energy_curvature_resolved

! **************************************************************************************************
!> \brief Euclidean geometry of an accepted vector secant
!> \param s2 squared step norm
!> \param y2 squared gradient-difference norm
!> \param sty secant curvature s^T y
!> \param snorm step norm
!> \param ynorm gradient-difference norm
!> \param alignment cosine between s and y
!> \param orthogonal_fraction fraction of y orthogonal to s
! **************************************************************************************************
   PURE SUBROUTINE ot_kpoint_secant_geometry(s2, y2, sty, snorm, ynorm, alignment, orthogonal_fraction)
      REAL(KIND=dp), INTENT(IN)                          :: s2, y2, sty
      REAL(KIND=dp), INTENT(OUT)                         :: snorm, ynorm, alignment, &
                                                            orthogonal_fraction

      snorm = SQRT(MAX(0.0_dp, s2))
      ynorm = SQRT(MAX(0.0_dp, y2))
      alignment = 0.0_dp
      IF (snorm > SQRT(TINY(1.0_dp)) .AND. ynorm > SQRT(TINY(1.0_dp))) THEN
         alignment = MAX(-1.0_dp, MIN(1.0_dp, sty/(snorm*ynorm)))
         IF (ABS(1.0_dp - ABS(alignment)) <= 16.0_dp*EPSILON(1.0_dp)) THEN
            alignment = SIGN(1.0_dp, alignment)
         END IF
      END IF
      orthogonal_fraction = SQRT(MAX(0.0_dp, 1.0_dp - alignment*alignment))

   END SUBROUTINE ot_kpoint_secant_geometry

! **************************************************************************************************
!> \brief reduce a k-point integer with MIN when a distributed communicator is present
!> \param qs_ot_env local OT channels
!> \param para_env_inter_kp ...
!> \param value ...
! **************************************************************************************************
   SUBROUTINE ot_mini_min_kpoint_integer(qs_ot_env, para_env_inter_kp, value)
      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: qs_ot_env
      TYPE(mp_para_env_type), OPTIONAL, POINTER          :: para_env_inter_kp
      INTEGER, INTENT(INOUT)                             :: value

      IF (ASSOCIATED(qs_ot_env(1)%para_env)) CALL qs_ot_env(1)%para_env%min(value)
      IF (PRESENT(para_env_inter_kp)) THEN
         IF (ASSOCIATED(para_env_inter_kp)) CALL para_env_inter_kp%min(value)
      END IF

   END SUBROUTINE ot_mini_min_kpoint_integer

! **************************************************************************************************
!> \brief reduce a k-point integer with MAX when a distributed communicator is present
!> \param qs_ot_env local OT channels
!> \param para_env_inter_kp ...
!> \param value ...
! **************************************************************************************************
   SUBROUTINE ot_mini_max_kpoint_integer(qs_ot_env, para_env_inter_kp, value)
      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: qs_ot_env
      TYPE(mp_para_env_type), OPTIONAL, POINTER          :: para_env_inter_kp
      INTEGER, INTENT(INOUT)                             :: value

      IF (ASSOCIATED(qs_ot_env(1)%para_env)) CALL qs_ot_env(1)%para_env%max(value)
      IF (PRESENT(para_env_inter_kp)) THEN
         IF (ASSOCIATED(para_env_inter_kp)) CALL para_env_inter_kp%max(value)
      END IF

   END SUBROUTINE ot_mini_max_kpoint_integer

! **************************************************************************************************
!> \brief reduce a k-point scalar with MAX when a distributed communicator is present
!> \param qs_ot_env local OT channels
!> \param para_env_inter_kp ...
!> \param value ...
! **************************************************************************************************
   SUBROUTINE ot_mini_max_kpoint_scalar(qs_ot_env, para_env_inter_kp, value)
      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: qs_ot_env
      TYPE(mp_para_env_type), OPTIONAL, POINTER          :: para_env_inter_kp
      REAL(KIND=dp), INTENT(INOUT)                       :: value

      IF (ASSOCIATED(qs_ot_env(1)%para_env)) CALL qs_ot_env(1)%para_env%max(value)
      IF (PRESENT(para_env_inter_kp)) THEN
         IF (ASSOCIATED(para_env_inter_kp)) CALL para_env_inter_kp%max(value)
      END IF

   END SUBROUTINE ot_mini_max_kpoint_scalar

! **************************************************************************************************
!> \brief reduce a k-point scalar with MIN when a distributed communicator is present
!> \param qs_ot_env local OT channels
!> \param para_env_inter_kp ...
!> \param value ...
! **************************************************************************************************
   SUBROUTINE ot_mini_min_kpoint_scalar(qs_ot_env, para_env_inter_kp, value)
      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: qs_ot_env
      TYPE(mp_para_env_type), OPTIONAL, POINTER          :: para_env_inter_kp
      REAL(KIND=dp), INTENT(INOUT)                       :: value

      IF (ASSOCIATED(qs_ot_env(1)%para_env)) CALL qs_ot_env(1)%para_env%min(value)
      IF (PRESENT(para_env_inter_kp)) THEN
         IF (ASSOCIATED(para_env_inter_kp)) CALL para_env_inter_kp%min(value)
      END IF

   END SUBROUTINE ot_mini_min_kpoint_scalar

! **************************************************************************************************
!> \brief sum a minimizer scalar over distributed k-point groups when requested
!> \param para_env_inter_kp ...
!> \param value ...
! **************************************************************************************************
   SUBROUTINE ot_mini_sum_kpoint_scalar(para_env_inter_kp, value)
      TYPE(mp_para_env_type), OPTIONAL, POINTER          :: para_env_inter_kp
      REAL(KIND=dp), INTENT(INOUT)                       :: value

      IF (PRESENT(para_env_inter_kp)) THEN
         IF (ASSOCIATED(para_env_inter_kp)) CALL para_env_inter_kp%sum(value)
      END IF

   END SUBROUTINE ot_mini_sum_kpoint_scalar

! **************************************************************************************************
!> \brief normalize a representative k-point direction to one full-grid member
!> \param multiplicity size of the irreducible k-point star
!> \return direction scale
! **************************************************************************************************
   PURE FUNCTION ot_kpoint_multiplicity_direction_scale(multiplicity) RESULT(weight)
      INTEGER, INTENT(IN)                                :: multiplicity
      REAL(KIND=dp)                                      :: weight

      weight = 1.0_dp/REAL(MAX(1, multiplicity), KIND=dp)

   END FUNCTION ot_kpoint_multiplicity_direction_scale

! **************************************************************************************************
!> \brief RMS metric weight for an irreducible k-point channel
!> \param qs_ot_env ...
!> \return ...
! **************************************************************************************************
   FUNCTION ot_kpoint_channel_norm_weight(qs_ot_env) RESULT(weight)
      TYPE(qs_ot_type), INTENT(IN)                       :: qs_ot_env
      REAL(KIND=dp)                                      :: weight

      weight = ot_kpoint_multiplicity_direction_scale(qs_ot_env%kpoint_multiplicity)

   END FUNCTION ot_kpoint_channel_norm_weight

! **************************************************************************************************
!> \brief full-grid-equivalent variable count for a minimizer channel set
!> \param qs_ot_env ...
!> \param do_ks ...
!> \param do_ener ...
!> \return ...
! **************************************************************************************************
   FUNCTION ot_kpoint_channel_nvar(qs_ot_env, do_ks, do_ener) RESULT(nvar)
      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: qs_ot_env
      LOGICAL, INTENT(IN)                                :: do_ks, do_ener
      REAL(KIND=dp)                                      :: nvar

      INTEGER                                            :: ichannel, k, n, nener
      REAL(KIND=dp)                                      :: multiplicity

      nvar = 0.0_dp
      DO ichannel = 1, SIZE(qs_ot_env)
         k = 0
         n = 0
         nener = 0
         IF (do_ks) THEN
            CALL dbcsr_get_info(qs_ot_env(ichannel)%matrix_x, nfullrows_total=n, nfullcols_total=k)
            IF (qs_ot_env(ichannel)%has_complex_kpoint_state) k = k + k
         END IF
         IF (do_ener) nener = SIZE(qs_ot_env(ichannel)%ener_x)
         multiplicity = REAL(MAX(1, qs_ot_env(ichannel)%kpoint_multiplicity), KIND=dp)
         nvar = nvar + multiplicity* &
                REAL(INT(n, KIND=int_8)*INT(k, KIND=int_8) + nener, KIND=dp)
      END DO

   END FUNCTION ot_kpoint_channel_nvar

!
! checks if the current point is a good point for finding a new direction
! or if we should improve the line_search, if it is used
!
! **************************************************************************************************
!> \brief ...
!> \param qs_ot_env ...
!> \return ...
! **************************************************************************************************
   FUNCTION current_point_is_fine(qs_ot_env) RESULT(res)
      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: qs_ot_env
      LOGICAL                                            :: res

      INTEGER                                            :: current_position_index
      LOGICAL                                            :: derivative_refinement, residual_rejected
      REAL(KIND=dp)                                      :: acceptance_tol, position_floor, &
                                                            refined_position

      res = .FALSE.

      ! only if we have a gradient it can be fine
      IF (.NOT. qs_ot_env(1)%energy_only) THEN

         ! we have not yet started with the line search
         IF (qs_ot_env(1)%line_search_count == 0) THEN
            res = .TRUE.
            RETURN
         END IF

         IF (qs_ot_env(1)%line_search_might_be_done) THEN
            acceptance_tol = 1.0E-12_dp* &
                             MAX(1.0_dp, ABS(qs_ot_env(1)%mermin_ls_reference_energy))
            IF (use_three_point_mermin_search(qs_ot_env) .AND. &
                qs_ot_env(1)%mermin_ls_sample_pending .AND. &
                qs_ot_env(1)%etotal > qs_ot_env(1)%mermin_ls_reference_energy + acceptance_tol) THEN
               RETURN
            END IF
            residual_rejected = kpoint_mermin_residual_closure_rejected(qs_ot_env)
            IF (residual_rejected) RETURN
            IF (use_three_point_mermin_search(qs_ot_env) .AND. &
                qs_ot_env(1)%mermin_ls_sample_pending .AND. &
                qs_ot_env(1)%line_search_count == 3) THEN
               current_position_index = qs_ot_env(1)%line_search_count + 1
               CPASSERT(current_position_index <= SIZE(qs_ot_env(1)%ot_pos))
               position_floor = MAX(1.0E4_dp*EPSILON(1.0_dp), &
                                    qs_ot_env(1)%settings%eps_taylor)* &
                                MAX(1.0_dp, ABS(qs_ot_env(1)%settings%ds_min))
               CALL ot_kpoint_mermin_derivative_refinement( &
                  qs_ot_env(1)%ot_grad(1), qs_ot_env(1)%gradient, &
                  qs_ot_env(1)%ot_pos(current_position_index), position_floor, &
                  refined_position, derivative_refinement)
               IF (derivative_refinement) RETURN
            END IF
            res = .TRUE.
            RETURN
         END IF

      END IF

   END FUNCTION current_point_is_fine

! **************************************************************************************************
!> \brief reject an accepted Mermin endpoint whose coupled residual left the local trust region
!> \param qs_ot_env ...
!> \param quality optional residual-consistency quality
!> \param growth_ratio optional effective residual growth ratio
!> \return ...
! **************************************************************************************************
   FUNCTION kpoint_mermin_residual_closure_rejected(qs_ot_env, quality, growth_ratio) RESULT(rejected)
      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: qs_ot_env
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: quality, growth_ratio
      LOGICAL                                            :: rejected

      LOGICAL                                            :: catastrophic, reset
      REAL(KIND=dp)                                      :: accepted_drop, current_delta, &
                                                            local_quality, local_ratio, &
                                                            previous_delta

      rejected = .FALSE.
      local_quality = 1.0_dp
      local_ratio = 0.0_dp
      IF (.NOT. qs_ot_env(1)%has_complex_kpoint_state .OR. &
          .NOT. qs_ot_occupation_preconditioner_active(qs_ot_env(1)) .OR. &
          .NOT. qs_ot_env(1)%settings%do_rotation .OR. &
          .NOT. qs_ot_env(1)%settings%do_ener .OR. &
          .NOT. qs_ot_env(1)%mermin_ls_sample_pending) GOTO 100

      previous_delta = ot_kpoint_mermin_reference_guard_delta( &
                       qs_ot_env(1)%schur_model_delta, qs_ot_env(1)%mermin_ls_reference_delta)
      current_delta = ot_kpoint_mermin_current_guard_delta( &
                      qs_ot_env(1)%delta, kpoint_projected_tail_delta(qs_ot_env), &
                      qs_ot_env(1)%mermin_ls_reference_delta)
      accepted_drop = qs_ot_env(1)%mermin_ls_reference_energy - qs_ot_env(1)%etotal
      IF (accepted_drop <= 0.0_dp) GOTO 100
      CALL ot_kpoint_cg_residual_growth_guard(previous_delta, current_delta, accepted_drop, &
                                              qs_ot_env(1)%kpoint_eps_scf, reset, &
                                              local_quality, catastrophic, &
                                              strict_growth= &
                                              ot_kpoint_mermin_shared_spin_group( &
                                              qs_ot_env%spin_index, &
                                              qs_ot_env%occupation_group_index) .AND. &
                                              qs_ot_env(1)%mermin_ls_reference_delta > 0.0_dp)
      IF (previous_delta > EPSILON(1.0_dp)) local_ratio = current_delta/previous_delta
      ! In the EPS-scaled convergence tail, residual growth means that the
      ! finite-response endpoint left the local trust region even when its
      ! sampled Mermin energy is lower.  Backtrack that endpoint; outside the
      ! local tail retain the existing catastrophic-growth protection.
      rejected = reset

100   CONTINUE
      IF (PRESENT(quality)) quality = local_quality
      IF (PRESENT(growth_ratio)) growth_ratio = local_ratio

   END FUNCTION kpoint_mermin_residual_closure_rejected

!
! performs various kinds of line searches
!
! **************************************************************************************************
!> \brief ...
!> \param qs_ot_env ...
! **************************************************************************************************
   SUBROUTINE do_line_search(qs_ot_env)
      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: qs_ot_env

      SELECT CASE (qs_ot_env(1)%settings%line_search_method)
      CASE ("GOLD")
         CALL do_line_search_gold(qs_ot_env)
      CASE ("3PNT")
         CALL do_line_search_3pnt(qs_ot_env)
      CASE ("2PNT")
         IF (use_three_point_mermin_search(qs_ot_env)) THEN
            CALL do_line_search_3pnt(qs_ot_env)
         ELSE
            CALL do_line_search_2pnt(qs_ot_env)
         END IF
      CASE ("ADPT")
         CALL do_line_search_adapt(qs_ot_env)
      CASE ("NONE")
         CALL do_line_search_none(qs_ot_env)
      CASE DEFAULT
         CPABORT("NYI")
      END SELECT
   END SUBROUTINE do_line_search

! **************************************************************************************************
!> \brief guard the strong complex Mermin preconditioner with an extra energy point
!> \param qs_ot_env ...
!> \return ...
! **************************************************************************************************
   FUNCTION use_three_point_mermin_search(qs_ot_env) RESULT(res)
      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: qs_ot_env
      LOGICAL                                            :: res

      INTEGER                                            :: ispin

      res = .FALSE.
      ! The occupation-preconditioner metric is deliberately bootstrapped
      ! after a subspace refresh.  That transient readiness state must not
      ! switch a running coupled Mermin search between 2PNT and 3PNT or expose
      ! the large conventional 2PNT interpolation step in the metallic tail.
      IF (.NOT. qs_ot_env(1)%settings%occupation_preconditioner .OR. &
          .NOT. qs_ot_env(1)%settings%do_rotation .OR. &
          .NOT. qs_ot_env(1)%settings%do_ener) RETURN
      DO ispin = 1, SIZE(qs_ot_env)
         IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
            res = .TRUE.
            RETURN
         END IF
      END DO

   END FUNCTION use_three_point_mermin_search

! **************************************************************************************************
!> \brief scale a k-point tail lower bound with the active SCF target when available
!> \param default_delta ...
!> \param eps_factor ...
!> \param eps_scf ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION kpoint_eps_scaled_tail_min(default_delta, eps_factor, eps_scf) RESULT(delta_min)
      REAL(KIND=dp), INTENT(IN)                          :: default_delta, eps_factor, eps_scf
      REAL(KIND=dp)                                      :: delta_min

      delta_min = default_delta
      IF (eps_scf > 0.0_dp) delta_min = MIN(delta_min, eps_factor*eps_scf)
      delta_min = MAX(delta_min, 1.0E4_dp*EPSILON(1.0_dp))

   END FUNCTION kpoint_eps_scaled_tail_min

! **************************************************************************************************
!> \brief scale a k-point tail upper bound with the active SCF target when available
!> \param default_delta ...
!> \param eps_factor ...
!> \param eps_scf ...
!> \param delta_min corresponding lower tail bound
!> \return ...
! **************************************************************************************************
   PURE FUNCTION kpoint_eps_scaled_tail_max(default_delta, eps_factor, eps_scf, delta_min) &
      RESULT(delta_max)
      REAL(KIND=dp), INTENT(IN)                          :: default_delta, eps_factor, eps_scf, &
                                                            delta_min
      REAL(KIND=dp)                                      :: delta_max

      delta_max = default_delta
      IF (eps_scf > 0.0_dp) delta_max = MIN(delta_max, eps_factor*eps_scf)
      delta_max = MAX(delta_max, delta_min)

   END FUNCTION kpoint_eps_scaled_tail_max

! **************************************************************************************************
!> \brief identify the nonconverged fixed-N tail from projected or effective residuals
!> \param projected_delta fixed-N projected residual
!> \param effective_delta residual used by the minimizer
!> \param eps_scf requested SCF convergence threshold
!> \return whether either residual lies in the recoverable physical tail
! **************************************************************************************************
   PURE FUNCTION ot_kpoint_mermin_recovery_active_tail(projected_delta, effective_delta, eps_scf) &
      RESULT(active)
      REAL(KIND=dp), INTENT(IN)                          :: projected_delta, effective_delta, eps_scf
      LOGICAL                                            :: active

      REAL(KIND=dp)                                      :: tail_delta, tail_min

      tail_min = kpoint_eps_scaled_tail_min( &
                 ot_cg_kp_mermin_recovery_tail_min_delta, &
                 ot_cg_kp_mermin_recovery_tail_min_eps_factor, eps_scf)
      tail_delta = projected_delta
      IF (tail_delta <= 0.0_dp) tail_delta = effective_delta
      active = ((tail_delta > tail_min .AND. &
                 tail_delta < ot_cg_kp_mermin_recovery_active_tail_max_delta) .OR. &
                (effective_delta > tail_min .AND. &
                 effective_delta < ot_cg_kp_mermin_recovery_active_tail_max_delta))

   END FUNCTION ot_kpoint_mermin_recovery_active_tail

! **************************************************************************************************
!> \brief return the fixed-N projected residual for physical k-point tail decisions
!> \param qs_ot_env ...
!> \return ...
! **************************************************************************************************
   FUNCTION kpoint_projected_tail_delta(qs_ot_env) RESULT(delta)
      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: qs_ot_env
      REAL(KIND=dp)                                      :: delta

      delta = qs_ot_env(1)%kpoint_physical_delta
      IF (qs_ot_env(1)%kpoint_projected_delta > 0.0_dp) THEN
         delta = qs_ot_env(1)%kpoint_projected_delta
      END IF

   END FUNCTION kpoint_projected_tail_delta

! **************************************************************************************************
!> \brief clear one completed Mermin search while retaining the current physical coordinates
!> \param qs_ot_env ...
!> \param restart_step first step for the replacement search direction
! **************************************************************************************************
   SUBROUTINE reset_kpoint_mermin_search(qs_ot_env, restart_step)
      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: qs_ot_env
      REAL(KIND=dp), INTENT(IN)                          :: restart_step

      qs_ot_env(1)%ot_pos(:) = 0.0_dp
      qs_ot_env(1)%ot_energy(:) = 0.0_dp
      qs_ot_env(1)%ot_grad(:) = 0.0_dp
      qs_ot_env(1)%line_search_left = -1
      qs_ot_env(1)%line_search_right = -1
      qs_ot_env(1)%line_search_mid = -1
      qs_ot_env(1)%line_search_count = 0
      qs_ot_env(1)%line_search_might_be_done = .FALSE.
      qs_ot_env(1)%energy_only = .FALSE.
      qs_ot_env(1)%mermin_ls_sample_pending = .FALSE.
      qs_ot_env(1)%mermin_ls_closure_rejected = .FALSE.
      qs_ot_env(1)%schur_closure_response_pending = .FALSE.
      qs_ot_env(1)%schur_closure_response_model_delta = 0.0_dp
      qs_ot_env(1)%schur_closure_response_model_density_delta = -1.0_dp
      qs_ot_env(1)%ds_min = restart_step

   END SUBROUTINE reset_kpoint_mermin_search

! **************************************************************************************************
!> \brief replace a repeatedly rejected, energy-flat endpoint by a fresh physical search direction
!> \param qs_ot_env ...
!> \param residual_closure_rejected result of the current endpoint residual preflight
! **************************************************************************************************
   SUBROUTINE recover_stagnant_mermin_preflight(qs_ot_env, residual_closure_rejected)
      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: qs_ot_env
      LOGICAL, INTENT(INOUT)                             :: residual_closure_rejected

      INTEGER                                            :: count
      LOGICAL                                            :: active_tail
      REAL(KIND=dp)                                      :: mermin_effective_delta, &
                                                            mermin_tail_delta, mermin_tail_min, &
                                                            restart_step, step_floor
      TYPE(cp_logger_type), POINTER                      :: logger

      IF (.NOT. residual_closure_rejected) RETURN
      IF (.NOT. use_three_point_mermin_search(qs_ot_env)) RETURN
      count = qs_ot_env(1)%line_search_count
      IF (.NOT. ot_kpoint_mermin_tail_energy_stagnant(qs_ot_env(1)%ot_energy, count)) RETURN

      mermin_tail_min = kpoint_eps_scaled_tail_min( &
                        ot_cg_kp_mermin_recovery_tail_min_delta, &
                        ot_cg_kp_mermin_recovery_tail_min_eps_factor, &
                        qs_ot_env(1)%kpoint_eps_scf)
      mermin_tail_delta = kpoint_projected_tail_delta(qs_ot_env)
      mermin_effective_delta = qs_ot_env(1)%delta
      IF (mermin_tail_delta <= 0.0_dp) mermin_tail_delta = mermin_effective_delta
      active_tail = ot_kpoint_mermin_recovery_active_tail( &
                    mermin_tail_delta, mermin_effective_delta, qs_ot_env(1)%kpoint_eps_scf)
      IF (.NOT. active_tail) RETURN
      IF (mermin_effective_delta > mermin_tail_min .AND. &
          mermin_effective_delta < ot_cg_kp_mermin_recovery_active_tail_max_delta) THEN
         mermin_tail_delta = mermin_effective_delta
      END IF

      step_floor = MAX(1.0E4_dp*EPSILON(1.0_dp), &
                       qs_ot_env(1)%settings%eps_taylor)* &
                   MAX(1.0_dp, ABS(qs_ot_env(1)%settings%ds_min))
      step_floor = MAX(step_floor, ot_cg_kp_mermin_recovery_step_floor* &
                       MAX(1.0_dp, ABS(qs_ot_env(1)%settings%ds_min)))
      restart_step = ot_kpoint_mermin_restart_step( &
                     qs_ot_env(1)%settings%ds_min, qs_ot_env(1)%ot_pos, &
                     qs_ot_env(1)%ot_energy, count, step_floor, force_contract=.TRUE.)
      restart_step = MAX(step_floor, restart_step)

      logger => cp_get_default_logger()
      IF (logger%iter_info%print_level >= debug_print_level .AND. logger%para_env%is_source()) THEN
         WRITE (cp_logger_get_default_unit_nr(logger), '(A,1X,I5,4(1X,ES16.8))') &
            " OT Mermin LS flat residual-preflight restart count delta old_step new_step floor:", &
            count, mermin_tail_delta, qs_ot_env(1)%ds_min, restart_step, step_floor
      END IF

      qs_ot_env(1)%schur_residual_history_reset_pending = .TRUE.
      qs_ot_env(1)%schur_residual_history_quality = 0.0_dp
      CALL reset_kpoint_mermin_search(qs_ot_env, restart_step)
      residual_closure_rejected = .FALSE.

   END SUBROUTINE recover_stagnant_mermin_preflight

! **************************************************************************************************
!> \brief restart a k-point Mermin line search that has collapsed to a numerically zero step
!> \param qs_ot_env ...
! **************************************************************************************************
   SUBROUTINE recover_collapsed_mermin_search(qs_ot_env)
      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: qs_ot_env

      INTEGER                                            :: count, npos
      LOGICAL                                            :: active_tail, collapsed_sample, &
                                                            collapsed_step, history_exhausted, &
                                                            stagnant_history
      REAL(KIND=dp)                                      :: mermin_effective_delta, &
                                                            mermin_tail_delta, mermin_tail_min, &
                                                            pos_floor, pos_max, restart_step, &
                                                            step_floor
      TYPE(cp_logger_type), POINTER                      :: logger

      IF (.NOT. use_three_point_mermin_search(qs_ot_env)) RETURN
      IF (qs_ot_env(1)%line_search_count <= 0) RETURN
      count = qs_ot_env(1)%line_search_count
      history_exhausted = count + 1 >= SIZE(qs_ot_env(1)%ot_pos)
      stagnant_history = ot_kpoint_mermin_search_stagnant(qs_ot_env(1)%ot_pos, &
                                                          qs_ot_env(1)%ot_energy, count)
      IF (ot_kpoint_mermin_recovery_deferred( &
          qs_ot_env(1)%mermin_ls_sample_pending, &
          qs_ot_env(1)%line_search_might_be_done, &
          qs_ot_env(1)%energy_only, history_exhausted, stagnant_history)) RETURN

      npos = MIN(count + 1, SIZE(qs_ot_env(1)%ot_pos))
      pos_max = 0.0_dp
      IF (npos > 0) pos_max = MAXVAL(ABS(qs_ot_env(1)%ot_pos(1:npos)))

      step_floor = MAX(1.0E4_dp*EPSILON(1.0_dp), &
                       qs_ot_env(1)%settings%eps_taylor)* &
                   MAX(1.0_dp, ABS(qs_ot_env(1)%settings%ds_min))
      pos_floor = step_floor
      mermin_tail_min = kpoint_eps_scaled_tail_min( &
                        ot_cg_kp_mermin_recovery_tail_min_delta, &
                        ot_cg_kp_mermin_recovery_tail_min_eps_factor, &
                        qs_ot_env(1)%kpoint_eps_scf)
      mermin_tail_delta = kpoint_projected_tail_delta(qs_ot_env)
      mermin_effective_delta = qs_ot_env(1)%delta
      IF (mermin_tail_delta <= 0.0_dp) mermin_tail_delta = mermin_effective_delta
      active_tail = ot_kpoint_mermin_recovery_active_tail( &
                    mermin_tail_delta, mermin_effective_delta, qs_ot_env(1)%kpoint_eps_scf)
      IF (active_tail .AND. mermin_effective_delta > mermin_tail_min .AND. &
          mermin_effective_delta < ot_cg_kp_mermin_recovery_active_tail_max_delta) THEN
         mermin_tail_delta = mermin_effective_delta
      END IF
      IF (active_tail) THEN
         step_floor = MAX(step_floor, ot_cg_kp_mermin_recovery_step_floor* &
                          MAX(1.0_dp, ABS(qs_ot_env(1)%settings%ds_min)))
         pos_floor = MAX(pos_floor, ot_cg_kp_mermin_recovery_active_step_floor* &
                         MAX(1.0_dp, ABS(qs_ot_env(1)%settings%ds_min)))
      END IF
      collapsed_step = ABS(qs_ot_env(1)%ds_min) <= step_floor
      collapsed_sample = active_tail .AND. count >= 3 .AND. pos_max <= pos_floor
      IF (.NOT. (collapsed_step .OR. collapsed_sample .OR. &
                 history_exhausted .OR. stagnant_history)) RETURN

      restart_step = ot_kpoint_mermin_restart_step(qs_ot_env(1)%settings%ds_min, &
                                                   qs_ot_env(1)%ot_pos, &
                                                   qs_ot_env(1)%ot_energy, count, step_floor, &
                                                   force_contract= &
                                                   qs_ot_env(1)%mermin_ls_closure_rejected)
      IF (history_exhausted .OR. stagnant_history) THEN
         restart_step = MAX(step_floor, ot_cg_kp_mermin_closure_backtrack* &
                            ABS(qs_ot_env(1)%settings%ds_min))
      END IF

      logger => cp_get_default_logger()
      IF (logger%iter_info%print_level >= debug_print_level .AND. logger%para_env%is_source()) THEN
         WRITE (cp_logger_get_default_unit_nr(logger), '(A,5(1X,L1),2(1X,ES16.8),1X,I5,4(1X,ES16.8))') &
            " OT Mermin LS recovery step_collapsed sample_collapsed active exhausted stagnant "// &
            "ds floor count pos pmax delta restart:", &
            collapsed_step, collapsed_sample, active_tail, history_exhausted, stagnant_history, &
            qs_ot_env(1)%ds_min, step_floor, count, pos_floor, pos_max, &
            mermin_tail_delta, restart_step
      END IF

      CALL reset_kpoint_mermin_search(qs_ot_env, restart_step)

   END SUBROUTINE recover_collapsed_mermin_search

! **************************************************************************************************
!> \brief fix the k-point auxiliary-energy gauge and trust scale once per search direction
!> \param qs_ot_env ...
!> \param para_env_inter_kp ...
!> \param fresh_response replace CG energy history by the newly coupled Schur response
!> \param model_slope_factor factor for the executed directional model slope
!> \param model_quadratic_factor factor for the executed directional model curvature
! **************************************************************************************************
   SUBROUTINE prepare_linear_kpoint_energy_direction(qs_ot_env, para_env_inter_kp, fresh_response, &
                                                     model_slope_factor, model_quadratic_factor)
      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: qs_ot_env
      TYPE(mp_para_env_type), OPTIONAL, POINTER          :: para_env_inter_kp
      LOGICAL, INTENT(IN)                                :: fresh_response
      REAL(KIND=dp), INTENT(OUT)                         :: model_slope_factor, &
                                                            model_quadratic_factor

      INTEGER                                            :: active_all, frontier_spin, group, i, &
                                                            ispin, ngroups, nspins
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: frontier_group_found
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: frontier_found
      LOGICAL                                            :: active, has_rayleigh_target, &
                                                            residual_aware_probe
      REAL(KIND=dp) :: alignment, combined_probe_max, coupled_probe, coupled_rms, &
         frontier_candidate, maximum_direction, minimum_target_norm2, norm_product, probe_pos, &
         rayleigh_correction, scale, shift, slope, target_rms, weight
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: coupled_norm2, direction_max, direction_scale, &
         energy_slope, final_energy_slope, fixed_orbital_fraction, frontier_group_scale, &
         frontier_position_max, orbital_slope, radius, rayleigh_diag, rayleigh_response_scale, &
         rayleigh_shift_sum, rayleigh_target_max, response_sum, target_coupled_dot, target_norm2, &
         weighted_direction
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: frontier_scale
      TYPE(cp_logger_type), POINTER                      :: logger

      model_slope_factor = 1.0_dp
      model_quadratic_factor = 1.0_dp
      active = .FALSE.
      DO ispin = 1, SIZE(qs_ot_env)
         qs_ot_env(ispin)%ener_direction_linear = .FALSE.
         qs_ot_env(ispin)%ener_x_fermi_guard_position_max = -1.0_dp
         qs_ot_env(ispin)%ener_x_line_search_position = 0.0_dp
         qs_ot_env(ispin)%ener_x_line_search_tangent = 1.0_dp
         active = active .OR. (qs_ot_env(ispin)%has_complex_kpoint_state .AND. &
                               qs_ot_occupation_preconditioner_active(qs_ot_env(ispin)) .AND. &
                               ASSOCIATED(qs_ot_env(ispin)%ener_dx))
      END DO
      active_all = MERGE(1, 0, active)
      CALL ot_mini_min_kpoint_integer(qs_ot_env, para_env_inter_kp, active_all)
      active = active_all == 1
      IF (.NOT. active) RETURN
      residual_aware_probe = ot_kpoint_residual_aware_rayleigh_model( &
                             fresh_response, qs_ot_env(1)%kpoint_energy_closure_attempted, &
                             qs_ot_env(1)%schur_model_flavor)

      ngroups = MAX(1, MAXVAL(qs_ot_env(:)%occupation_group_index))
      nspins = MAX(1, MAXVAL(qs_ot_env(:)%spin_index))
      IF (PRESENT(para_env_inter_kp)) THEN
         IF (ASSOCIATED(para_env_inter_kp)) THEN
            CALL para_env_inter_kp%max(ngroups)
            CALL para_env_inter_kp%max(nspins)
         END IF
      END IF
      ALLOCATE (coupled_norm2(ngroups), direction_max(ngroups), direction_scale(ngroups), &
                energy_slope(ngroups), final_energy_slope(ngroups), &
                fixed_orbital_fraction(ngroups), orbital_slope(ngroups), &
                frontier_group_found(ngroups), frontier_group_scale(ngroups), &
                frontier_position_max(ngroups), frontier_found(ngroups, nspins), &
                frontier_scale(ngroups, nspins), &
                radius(ngroups), rayleigh_shift_sum(ngroups), rayleigh_target_max(ngroups), &
                rayleigh_response_scale(ngroups), response_sum(ngroups), &
                target_coupled_dot(ngroups), target_norm2(ngroups), weighted_direction(ngroups))
      coupled_norm2(:) = 0.0_dp
      direction_max(:) = 0.0_dp
      direction_scale(:) = 1.0_dp
      energy_slope(:) = 0.0_dp
      final_energy_slope(:) = 0.0_dp
      fixed_orbital_fraction(:) = 0.0_dp
      orbital_slope(:) = 0.0_dp
      radius(:) = HUGE(1.0_dp)
      frontier_group_found(:) = 0
      frontier_group_scale(:) = 1.0_dp
      frontier_position_max(:) = -1.0_dp
      frontier_found(:, :) = 0
      frontier_scale(:, :) = 0.0_dp
      rayleigh_shift_sum(:) = 0.0_dp
      rayleigh_response_scale(:) = 1.0_dp
      rayleigh_target_max(:) = 0.0_dp
      response_sum(:) = 0.0_dp
      target_coupled_dot(:) = 0.0_dp
      target_norm2(:) = 0.0_dp
      weighted_direction(:) = 0.0_dp

      DO ispin = 1, SIZE(qs_ot_env)
         IF (.NOT. qs_ot_env(ispin)%has_complex_kpoint_state) CYCLE
         group = MAX(1, MIN(ngroups, qs_ot_env(ispin)%occupation_group_index))
         ! The exact finite energy correction is rebuilt for every accepted
         ! orbital state.  Do not carry its old CG component into a new
         ! rotation direction; ener_gx_old is the current coupled Schur
         ! response and remains the history-safe base direction.
         IF (fresh_response .AND. ASSOCIATED(qs_ot_env(ispin)%ener_gx_old)) THEN
            qs_ot_env(ispin)%ener_dx(:) = -qs_ot_env(ispin)%ener_gx_old(:)
         END IF
         IF (qs_ot_env(ispin)%ener_x_trust_radius > 0.0_dp) THEN
            radius(group) = MIN(radius(group), qs_ot_env(ispin)%ener_x_trust_radius)
         END IF
         IF (.NOT. ASSOCIATED(qs_ot_env(ispin)%ener_x_trust_weight)) CYCLE
         IF (SIZE(qs_ot_env(ispin)%ener_x_trust_weight) /= SIZE(qs_ot_env(ispin)%ener_dx)) CYCLE
         has_rayleigh_target = ASSOCIATED(qs_ot_env(ispin)%rot_mat_chc) .AND. &
                               SIZE(qs_ot_env(ispin)%ener_x) == SIZE(qs_ot_env(ispin)%ener_dx)
         IF (has_rayleigh_target) THEN
            ALLOCATE (rayleigh_diag(SIZE(qs_ot_env(ispin)%ener_dx)))
            CALL dbcsr_get_diag(qs_ot_env(ispin)%rot_mat_chc, rayleigh_diag)
            IF (ASSOCIATED(qs_ot_env(ispin)%para_env)) THEN
               CALL qs_ot_env(ispin)%para_env%sum(rayleigh_diag)
            END IF
         END IF
         DO i = 1, SIZE(qs_ot_env(ispin)%ener_dx)
            weight = MAX(0.0_dp, qs_ot_env(ispin)%ener_x_trust_weight(i))
            response_sum(group) = response_sum(group) + weight
            weighted_direction(group) = weighted_direction(group) + &
                                        weight*qs_ot_env(ispin)%ener_dx(i)
            IF (has_rayleigh_target) THEN
               rayleigh_shift_sum(group) = rayleigh_shift_sum(group) + &
                                           weight*(qs_ot_env(ispin)%ener_x(i) - rayleigh_diag(i))
            END IF
         END DO
         IF (has_rayleigh_target) DEALLOCATE (rayleigh_diag)
      END DO
      IF (PRESENT(para_env_inter_kp)) THEN
         IF (ASSOCIATED(para_env_inter_kp)) THEN
            CALL para_env_inter_kp%sum(response_sum)
            CALL para_env_inter_kp%sum(weighted_direction)
            CALL para_env_inter_kp%sum(rayleigh_shift_sum)
            CALL para_env_inter_kp%min(radius)
         END IF
      END IF

      ! A common shift is the chemical-potential gauge mode.  Remove one
      ! global mode per fixed-occupation group, not one mode per k-point.
      DO ispin = 1, SIZE(qs_ot_env)
         IF (.NOT. qs_ot_env(ispin)%has_complex_kpoint_state) CYCLE
         group = MAX(1, MIN(ngroups, qs_ot_env(ispin)%occupation_group_index))
         shift = 0.0_dp
         IF (response_sum(group) > EPSILON(response_sum(group))) THEN
            shift = weighted_direction(group)/response_sum(group)
         END IF
         qs_ot_env(ispin)%ener_dx(:) = qs_ot_env(ispin)%ener_dx(:) - shift
      END DO

      ! The Rayleigh target is exact only for fixed orbitals.  Measure that
      ! asymptotic condition with directional Mermin slopes in the common
      ! physical energy unit, including both complex rotation channels.
      DO ispin = 1, SIZE(qs_ot_env)
         IF (.NOT. qs_ot_env(ispin)%has_complex_kpoint_state) CYCLE
         group = MAX(1, MIN(ngroups, qs_ot_env(ispin)%occupation_group_index))
         IF (ASSOCIATED(qs_ot_env(ispin)%matrix_gx) .AND. &
             ASSOCIATED(qs_ot_env(ispin)%matrix_dx)) THEN
            CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin)%matrix_dx, slope)
            orbital_slope(group) = orbital_slope(group) + ABS(slope)
         END IF
         IF (ASSOCIATED(qs_ot_env(ispin)%matrix_gx_im) .AND. &
             ASSOCIATED(qs_ot_env(ispin)%matrix_dx_im)) THEN
            CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx_im, qs_ot_env(ispin)%matrix_dx_im, slope)
            orbital_slope(group) = orbital_slope(group) + ABS(slope)
         END IF
         IF (ASSOCIATED(qs_ot_env(ispin)%rot_mat_gx) .AND. &
             ASSOCIATED(qs_ot_env(ispin)%rot_mat_dx)) THEN
            CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx, qs_ot_env(ispin)%rot_mat_dx, slope)
            orbital_slope(group) = orbital_slope(group) + 0.5_dp*ABS(slope)
         END IF
         IF (ASSOCIATED(qs_ot_env(ispin)%rot_mat_gx_im) .AND. &
             ASSOCIATED(qs_ot_env(ispin)%rot_mat_dx_im)) THEN
            CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx_im, qs_ot_env(ispin)%rot_mat_dx_im, slope)
            orbital_slope(group) = orbital_slope(group) + 0.5_dp*ABS(slope)
         END IF
         energy_slope(group) = energy_slope(group) + &
                               ABS(DOT_PRODUCT(qs_ot_env(ispin)%ener_gx, &
                                               qs_ot_env(ispin)%ener_dx))
      END DO
      IF (PRESENT(para_env_inter_kp)) THEN
         IF (ASSOCIATED(para_env_inter_kp)) THEN
            CALL para_env_inter_kp%sum(energy_slope)
            CALL para_env_inter_kp%sum(orbital_slope)
         END IF
      END IF
      DO group = 1, ngroups
         IF (energy_slope(group) + orbital_slope(group) > 0.0_dp) THEN
            fixed_orbital_fraction(group) = energy_slope(group)/ &
                                            (energy_slope(group) + orbital_slope(group))
         END IF
      END DO

      ! The 2PNT interpolation may place its final sample up to four times
      ! beyond the initial one.  Project the auxiliary-energy direction for
      ! that complete affine interval, while leaving the orbital/rotation
      ! components at their ordinary line-search scale.
      probe_pos = ot_kpoint_rayleigh_probe_position( &
                  qs_ot_env(1)%ds_min, qs_ot_env(1)%settings%ds_min, residual_aware_probe)

      ! Audit how far the current affine direction could retain one thermal
      ! frontier state per spin.  A single extremal state must not scale the
      ! complete distributed auxiliary-energy direction: the guarded trial-state
      ! projection below invalidates and rebuilds minimizer history when needed.
      DO ispin = 1, SIZE(qs_ot_env)
         IF (.NOT. qs_ot_env(ispin)%has_complex_kpoint_state) CYCLE
         IF (.NOT. qs_ot_env(ispin)%ener_x_fermi_guard_active) CYCLE
         IF (qs_ot_env(ispin)%ener_x_fermi_guard_anchor <= 0.0_dp) CYCLE
         group = MAX(1, MIN(ngroups, qs_ot_env(ispin)%occupation_group_index))
         frontier_spin = MAX(1, MIN(nspins, qs_ot_env(ispin)%spin_index))
         DO i = 1, SIZE(qs_ot_env(ispin)%ener_dx)
            frontier_candidate = ot_kpoint_fermi_window_direction_scale( &
                                 qs_ot_env(ispin)%ener_x(i) - qs_ot_env(ispin)%ener_x_trust_mu, &
                                 probe_pos*qs_ot_env(ispin)%ener_dx(i), &
                                 qs_ot_env(ispin)%ener_x_fermi_guard_anchor)
            IF (frontier_candidate < 0.0_dp) CYCLE
            frontier_found(group, frontier_spin) = 1
            frontier_scale(group, frontier_spin) = &
               MAX(frontier_scale(group, frontier_spin), frontier_candidate)
         END DO
      END DO
      IF (PRESENT(para_env_inter_kp)) THEN
         IF (ASSOCIATED(para_env_inter_kp)) THEN
            DO frontier_spin = 1, nspins
               DO group = 1, ngroups
                  CALL para_env_inter_kp%max(frontier_found(group, frontier_spin))
                  CALL para_env_inter_kp%max(frontier_scale(group, frontier_spin))
               END DO
            END DO
         END IF
      END IF
      DO frontier_spin = 1, nspins
         DO group = 1, ngroups
            IF (frontier_found(group, frontier_spin) > 0) THEN
               frontier_group_found(group) = 1
            END IF
            IF (frontier_found(group, frontier_spin) > 0) THEN
               frontier_group_scale(group) = &
                  MIN(frontier_group_scale(group), frontier_scale(group, frontier_spin))
            END IF
         END DO
      END DO

      ! Audit the exact fixed-orbital Rayleigh correction against the coupled
      ! Schur response.  The common response-weighted shift is one fixed-N
      ! chemical-potential gauge mode per occupation group.  A full Rayleigh
      ! replacement is deliberately not applied here: H changes with the
      ! density, so its finite response must first be represented by the joint
      ! rotation/energy model and accepted by the same Mermin line search.
      DO ispin = 1, SIZE(qs_ot_env)
         IF (.NOT. qs_ot_env(ispin)%has_complex_kpoint_state) CYCLE
         IF (.NOT. ASSOCIATED(qs_ot_env(ispin)%rot_mat_chc)) CYCLE
         IF (.NOT. ASSOCIATED(qs_ot_env(ispin)%ener_x_trust_weight)) CYCLE
         IF (SIZE(qs_ot_env(ispin)%ener_x_trust_weight) /= SIZE(qs_ot_env(ispin)%ener_dx)) CYCLE
         group = MAX(1, MIN(ngroups, qs_ot_env(ispin)%occupation_group_index))
         IF (response_sum(group) <= EPSILON(response_sum(group))) CYCLE
         shift = rayleigh_shift_sum(group)/response_sum(group)
         ALLOCATE (rayleigh_diag(SIZE(qs_ot_env(ispin)%ener_dx)))
         CALL dbcsr_get_diag(qs_ot_env(ispin)%rot_mat_chc, rayleigh_diag)
         IF (ASSOCIATED(qs_ot_env(ispin)%para_env)) THEN
            CALL qs_ot_env(ispin)%para_env%sum(rayleigh_diag)
         END IF
         DO i = 1, SIZE(qs_ot_env(ispin)%ener_dx)
            weight = MAX(0.0_dp, qs_ot_env(ispin)%ener_x_trust_weight(i))
            IF (weight <= 0.0_dp) CYCLE
            rayleigh_correction = rayleigh_diag(i) + shift - qs_ot_env(ispin)%ener_x(i)
            coupled_probe = probe_pos*qs_ot_env(ispin)%ener_dx(i)
            rayleigh_target_max(group) = MAX(rayleigh_target_max(group), ABS(rayleigh_correction))
            target_norm2(group) = target_norm2(group) + weight*rayleigh_correction**2
            coupled_norm2(group) = coupled_norm2(group) + weight*coupled_probe**2
            target_coupled_dot(group) = target_coupled_dot(group) + &
                                        weight*rayleigh_correction*coupled_probe
         END DO
         DEALLOCATE (rayleigh_diag)
      END DO
      IF (PRESENT(para_env_inter_kp)) THEN
         IF (ASSOCIATED(para_env_inter_kp)) THEN
            CALL para_env_inter_kp%max(rayleigh_target_max)
            CALL para_env_inter_kp%sum(target_norm2)
            CALL para_env_inter_kp%sum(coupled_norm2)
            CALL para_env_inter_kp%sum(target_coupled_dot)
         END IF
      END IF

      ! The coupled response already has the correct finite-complex direction in
      ! the fixed-orbital tail. Calibrate only its amplitude so that the complete
      ! 2PNT search interval reaches the exact fixed-N Rayleigh target. This
      ! leaves the REF basis and nonlinear-CG history untouched; the ordinary
      ! Mermin line search still accepts or rejects the resulting direction.
      IF (fresh_response) THEN
         DO group = 1, ngroups
            minimum_target_norm2 = 0.0_dp
            IF (radius(group) < HUGE(1.0_dp) .AND. response_sum(group) > 0.0_dp) THEN
               minimum_target_norm2 = response_sum(group)* &
                                      (SQRT(EPSILON(1.0_dp))*radius(group))**2
            END IF
            IF (residual_aware_probe) THEN
               rayleigh_response_scale(group) = ot_kpoint_rayleigh_response_scale( &
                                                target_norm2(group), coupled_norm2(group), &
                                                target_coupled_dot(group), &
                                                fixed_orbital_fraction(group), minimum_target_norm2, &
                                                maximum_scale=1.0_dp/probe_pos)
            ELSE
               rayleigh_response_scale(group) = ot_kpoint_rayleigh_response_scale( &
                                                target_norm2(group), coupled_norm2(group), &
                                                target_coupled_dot(group), &
                                                fixed_orbital_fraction(group), minimum_target_norm2)
            END IF
         END DO
         DO ispin = 1, SIZE(qs_ot_env)
            IF (.NOT. qs_ot_env(ispin)%has_complex_kpoint_state) CYCLE
            group = MAX(1, MIN(ngroups, qs_ot_env(ispin)%occupation_group_index))
            qs_ot_env(ispin)%ener_dx(:) = rayleigh_response_scale(group)* &
                                          qs_ot_env(ispin)%ener_dx(:)
         END DO
      END IF

      DO ispin = 1, SIZE(qs_ot_env)
         IF (.NOT. qs_ot_env(ispin)%has_complex_kpoint_state) CYCLE
         group = MAX(1, MIN(ngroups, qs_ot_env(ispin)%occupation_group_index))
         direction_max(group) = MAX(direction_max(group), &
                                    probe_pos*MAXVAL(ABS(qs_ot_env(ispin)%ener_dx)))
      END DO
      IF (PRESENT(para_env_inter_kp)) THEN
         IF (ASSOCIATED(para_env_inter_kp)) CALL para_env_inter_kp%max(direction_max)
      END IF

      DO group = 1, ngroups
         scale = 1.0_dp
         IF (.NOT. residual_aware_probe .AND. radius(group) < HUGE(1.0_dp) .AND. &
             direction_max(group) > radius(group)) THEN
            scale = radius(group)/direction_max(group)
         END IF
         direction_scale(group) = scale
         IF (frontier_group_found(group) > 0 .AND. scale > EPSILON(scale)) THEN
            frontier_position_max(group) = probe_pos* &
                                           MIN(1.0_dp, frontier_group_scale(group)/scale)
         END IF
      END DO

      DO ispin = 1, SIZE(qs_ot_env)
         IF (.NOT. qs_ot_env(ispin)%has_complex_kpoint_state) CYCLE
         group = MAX(1, MIN(ngroups, qs_ot_env(ispin)%occupation_group_index))
         scale = direction_scale(group)
         IF (residual_aware_probe .AND. radius(group) < HUGE(1.0_dp)) THEN
            maximum_direction = radius(group)/probe_pos
            DO i = 1, SIZE(qs_ot_env(ispin)%ener_dx)
               IF (ABS(qs_ot_env(ispin)%ener_dx(i)) > maximum_direction) THEN
                  scale = MIN(scale, maximum_direction/ABS(qs_ot_env(ispin)%ener_dx(i)))
                  qs_ot_env(ispin)%ener_dx(i) = SIGN(maximum_direction, &
                                                     qs_ot_env(ispin)%ener_dx(i))
               END IF
            END DO
         ELSE
            qs_ot_env(ispin)%ener_dx(:) = scale*qs_ot_env(ispin)%ener_dx(:)
         END IF
         qs_ot_env(ispin)%ener_x_trust_scale = MIN(qs_ot_env(ispin)%ener_x_trust_scale, scale)
         qs_ot_env(ispin)%ener_x_fermi_guard_position_max = frontier_position_max(group)
         qs_ot_env(ispin)%ener_direction_linear = .TRUE.
      END DO
      IF (fresh_response) THEN
         final_energy_slope(:) = 0.0_dp
         DO ispin = 1, SIZE(qs_ot_env)
            IF (.NOT. qs_ot_env(ispin)%has_complex_kpoint_state) CYCLE
            group = MAX(1, MIN(ngroups, qs_ot_env(ispin)%occupation_group_index))
            final_energy_slope(group) = final_energy_slope(group) + &
                                        ABS(DOT_PRODUCT(qs_ot_env(ispin)%ener_gx, &
                                                        qs_ot_env(ispin)%ener_dx))
         END DO
         IF (PRESENT(para_env_inter_kp)) THEN
            IF (ASSOCIATED(para_env_inter_kp)) CALL para_env_inter_kp%sum(final_energy_slope)
         END IF
         DO group = 1, ngroups
            IF (energy_slope(group) > EPSILON(energy_slope(group))) THEN
               rayleigh_response_scale(group) = final_energy_slope(group)/energy_slope(group)
            ELSE
               rayleigh_response_scale(group) = 1.0_dp
            END IF
         END DO
         model_slope_factor = ot_kpoint_rayleigh_model_factor( &
                              energy_slope, orbital_slope, rayleigh_response_scale, 1)
         model_quadratic_factor = ot_kpoint_rayleigh_model_factor( &
                                  energy_slope, orbital_slope, rayleigh_response_scale, 2)
      END IF

      logger => cp_get_default_logger()
      IF (logger%iter_info%print_level >= debug_print_level .AND. logger%para_env%is_source()) THEN
         DO group = 1, ngroups
            combined_probe_max = direction_max(group)
            target_rms = SQRT(MAX(0.0_dp, target_norm2(group))/ &
                              MAX(response_sum(group), EPSILON(response_sum(group))))
            coupled_rms = SQRT(MAX(0.0_dp, coupled_norm2(group))/ &
                               MAX(response_sum(group), EPSILON(response_sum(group))))
            alignment = 0.0_dp
            norm_product = SQRT(MAX(0.0_dp, target_norm2(group)))* &
                           SQRT(MAX(0.0_dp, coupled_norm2(group)))
            IF (norm_product > 0.0_dp) THEN
               alignment = target_coupled_dot(group)/norm_product
            END IF
            WRITE (cp_logger_get_default_unit_nr(logger), '(A,1X,I3,12(1X,ES16.8))') &
               " OT fixed-N Rayleigh response audit group tmax trms srms align efrac rscale "// &
               "smax radius probe response fscale pmax:", &
               group, rayleigh_target_max(group), target_rms, coupled_rms, alignment, &
               fixed_orbital_fraction(group), rayleigh_response_scale(group), combined_probe_max, &
               radius(group), probe_pos, response_sum(group), frontier_group_scale(group), &
               frontier_position_max(group)
         END DO
      END IF

      DEALLOCATE (coupled_norm2, direction_max, direction_scale, energy_slope, final_energy_slope, &
                  fixed_orbital_fraction, frontier_found, &
                  frontier_group_found, frontier_group_scale, frontier_position_max, &
                  frontier_scale, orbital_slope, radius, rayleigh_response_scale, rayleigh_shift_sum, &
                  rayleigh_target_max, response_sum, target_coupled_dot, target_norm2, &
                  weighted_direction)

   END SUBROUTINE prepare_linear_kpoint_energy_direction

! **************************************************************************************************
!> \brief moves x adding the right amount (ds) of the gradient or search direction
!> \param ds ...
!> \param qs_ot_env ...
!> \par History
!>      08.2004 created [ Joost VandeVondele ] copied here from a larger number of subroutines
! **************************************************************************************************
   SUBROUTINE take_step(ds, qs_ot_env)
      REAL(KIND=dp)                                      :: ds
      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: qs_ot_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'take_step'

      INTEGER                                            :: handle, ispin, nspin
      LOGICAL                                            :: do_ener, do_ks

      CALL timeset(routineN, handle)

      nspin = SIZE(qs_ot_env)

      do_ks = qs_ot_env(1)%settings%ks
      do_ener = qs_ot_env(1)%settings%do_ener

      ! now update x to take into account this new step
      ! either dx or -gx is the direction to use
      IF (qs_ot_env(1)%use_dx) THEN
         IF (do_ks) THEN
            DO ispin = 1, nspin
               CALL dbcsr_add(qs_ot_env(ispin)%matrix_x, qs_ot_env(ispin)%matrix_dx, &
                              alpha_scalar=1.0_dp, beta_scalar=ds)
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  CALL dbcsr_add(qs_ot_env(ispin)%matrix_x_im, qs_ot_env(ispin)%matrix_dx_im, &
                                 alpha_scalar=1.0_dp, beta_scalar=ds)
               END IF
               IF (qs_ot_env(ispin)%settings%do_rotation) THEN
                  CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_x, qs_ot_env(ispin)%rot_mat_dx, &
                                 alpha_scalar=1.0_dp, beta_scalar=ds)
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_x_im, qs_ot_env(ispin)%rot_mat_dx_im, &
                                    alpha_scalar=1.0_dp, beta_scalar=ds)
                  END IF
               END IF
            END DO
         END IF
         IF (do_ener) THEN
            DO ispin = 1, nspin
               CALL take_energy_step(qs_ot_env(ispin), ds, qs_ot_env(ispin)%ener_dx)
            END DO
         END IF
      ELSE
         IF (do_ks) THEN
            DO ispin = 1, nspin
               CALL dbcsr_add(qs_ot_env(ispin)%matrix_x, qs_ot_env(ispin)%matrix_gx, &
                              alpha_scalar=1.0_dp, beta_scalar=-ds)
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  CALL dbcsr_add(qs_ot_env(ispin)%matrix_x_im, qs_ot_env(ispin)%matrix_gx_im, &
                                 alpha_scalar=1.0_dp, beta_scalar=-ds)
               END IF
               IF (qs_ot_env(ispin)%settings%do_rotation) THEN
                  CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_x, qs_ot_env(ispin)%rot_mat_gx, &
                                 alpha_scalar=1.0_dp, beta_scalar=-ds)
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_x_im, qs_ot_env(ispin)%rot_mat_gx_im, &
                                    alpha_scalar=1.0_dp, beta_scalar=-ds)
                  END IF
               END IF
            END DO
         END IF
         IF (do_ener) THEN
            DO ispin = 1, nspin
               CALL take_energy_step(qs_ot_env(ispin), -ds, qs_ot_env(ispin)%ener_gx)
            END DO
         END IF
      END IF
      CALL timestop(handle)
   END SUBROUTINE take_step

! **************************************************************************************************
!> \brief apply an optional trust radius to auxiliary orbital-energy steps
!> \param qs_ot_env ...
!> \param alpha ...
!> \param direction ...
! **************************************************************************************************
   SUBROUTINE take_energy_step(qs_ot_env, alpha, direction)
      TYPE(qs_ot_type), INTENT(INOUT)                    :: qs_ot_env
      REAL(KIND=dp), INTENT(IN)                          :: alpha
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: direction

      INTEGER                                            :: i, n
      LOGICAL                                            :: has_response_trust
      REAL(KIND=dp) :: actual_alpha, max_projected_step, next_position, raw_step_max, &
         response_weight, step_scale, step_shift, trusted_step_max
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: limited_step, raw_step

      n = SIZE(direction)
      IF (n == 0) RETURN

      IF (qs_ot_env%ener_direction_linear) THEN
         raw_step_max = MAXVAL(ABS(alpha*direction))
         CALL ot_kpoint_fermi_window_energy_coordinate( &
            qs_ot_env%ener_x_line_search_position, alpha, &
            qs_ot_env%ener_x_fermi_guard_position_max, next_position, actual_alpha, &
            qs_ot_env%ener_x_line_search_tangent)
         qs_ot_env%ener_x(:) = qs_ot_env%ener_x(:) + actual_alpha*direction(:)
         qs_ot_env%ener_x_line_search_position = next_position
         qs_ot_env%ener_x_trusted_step_max = MAXVAL(ABS(actual_alpha*direction))
         IF (raw_step_max > EPSILON(raw_step_max)) THEN
            qs_ot_env%ener_x_trust_scale = MIN(qs_ot_env%ener_x_trust_scale, &
                                               qs_ot_env%ener_x_trusted_step_max/raw_step_max)
         END IF
         RETURN
      END IF

      ALLOCATE (raw_step(n), limited_step(n))
      raw_step(:) = alpha*direction(:)
      limited_step(:) = raw_step(:)
      raw_step_max = MAXVAL(ABS(raw_step))
      trusted_step_max = raw_step_max
      step_scale = 1.0_dp

      IF (qs_ot_env%ener_x_trust_radius > 0.0_dp) THEN
         has_response_trust = .FALSE.
         IF (ASSOCIATED(qs_ot_env%ener_x_trust_weight)) THEN
            has_response_trust = SIZE(qs_ot_env%ener_x_trust_weight) == n
         END IF

         IF (has_response_trust) THEN
            response_weight = 0.0_dp
            step_shift = 0.0_dp
            DO i = 1, n
               IF (qs_ot_env%ener_x_trust_weight(i) > 0.0_dp) THEN
                  response_weight = response_weight + qs_ot_env%ener_x_trust_weight(i)
                  step_shift = step_shift + qs_ot_env%ener_x_trust_weight(i)*raw_step(i)
               END IF
            END DO

            IF (response_weight > EPSILON(response_weight)) THEN
               step_shift = step_shift/response_weight
               max_projected_step = 0.0_dp
               DO i = 1, n
                  IF (qs_ot_env%ener_x_trust_weight(i) > 0.0_dp) THEN
                     max_projected_step = MAX(max_projected_step, ABS(raw_step(i) - step_shift))
                  END IF
               END DO
               IF (max_projected_step > qs_ot_env%ener_x_trust_radius) THEN
                  step_scale = qs_ot_env%ener_x_trust_radius/max_projected_step
               END IF
               DO i = 1, n
                  IF (qs_ot_env%ener_x_trust_weight(i) > 0.0_dp) THEN
                     limited_step(i) = step_scale*(raw_step(i) - step_shift)
                  END IF
               END DO
            END IF
         END IF

         DO i = 1, n
            IF (ABS(limited_step(i)) > qs_ot_env%ener_x_trust_radius) THEN
               limited_step(i) = SIGN(qs_ot_env%ener_x_trust_radius, limited_step(i))
            END IF
         END DO
         trusted_step_max = MAXVAL(ABS(limited_step))
         IF (raw_step_max > EPSILON(raw_step_max)) step_scale = trusted_step_max/raw_step_max
      END IF
      qs_ot_env%ener_x = qs_ot_env%ener_x + limited_step
      qs_ot_env%ener_x_trusted_step_max = trusted_step_max
      qs_ot_env%ener_x_trust_scale = MIN(qs_ot_env%ener_x_trust_scale, step_scale)

      DEALLOCATE (raw_step, limited_step)

   END SUBROUTINE take_energy_step

! implements a golden ratio search as a robust way of minimizing
! **************************************************************************************************
!> \brief ...
!> \param qs_ot_env ...
! **************************************************************************************************
   SUBROUTINE do_line_search_gold(qs_ot_env)

      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: qs_ot_env

      CHARACTER(len=*), PARAMETER :: routineN = 'do_line_search_gold'
      REAL(KIND=dp), PARAMETER                           :: gold_sec = 0.3819_dp

      INTEGER                                            :: count, handle
      REAL(KIND=dp)                                      :: ds

      CALL timeset(routineN, handle)

      qs_ot_env(1)%line_search_count = qs_ot_env(1)%line_search_count + 1
      count = qs_ot_env(1)%line_search_count
      qs_ot_env(1)%line_search_might_be_done = .FALSE.
      qs_ot_env(1)%energy_only = .TRUE.

      IF (count + 1 > SIZE(qs_ot_env(1)%OT_pos)) THEN
         ! should not happen, we pass with a warning first
         ! you can increase the size of OT_pos and the like in qs_ot_env
         CPABORT("MAX ITER EXCEEDED : FATAL")
      END IF

      IF (qs_ot_env(1)%line_search_count == 1) THEN
         qs_ot_env(1)%line_search_left = 1
         qs_ot_env(1)%line_search_right = 0
         qs_ot_env(1)%line_search_mid = 1
         qs_ot_env(1)%ot_pos(1) = 0.0_dp
         qs_ot_env(1)%ot_energy(1) = qs_ot_env(1)%etotal
         qs_ot_env(1)%ot_pos(2) = qs_ot_env(1)%ds_min/gold_sec
      ELSE
         qs_ot_env(1)%ot_energy(count) = qs_ot_env(1)%etotal
         ! it's essentially a book keeping game.
         ! keep left on the left, keep (bring) right on the right
         ! and mid in between these two
         IF (qs_ot_env(1)%line_search_right == 0) THEN ! we do not yet have the right bracket
            IF (qs_ot_env(1)%ot_energy(count - 1) < qs_ot_env(1)%ot_energy(count)) THEN
               qs_ot_env(1)%line_search_right = count
               qs_ot_env(1)%ot_pos(count + 1) = qs_ot_env(1)%ot_pos(qs_ot_env(1)%line_search_mid) + &
                                                (qs_ot_env(1)%ot_pos(qs_ot_env(1)%line_search_right) - &
                                                 qs_ot_env(1)%ot_pos(qs_ot_env(1)%line_search_mid))*gold_sec
            ELSE
               qs_ot_env(1)%line_search_left = qs_ot_env(1)%line_search_mid
               qs_ot_env(1)%line_search_mid = count
               qs_ot_env(1)%ot_pos(count + 1) = qs_ot_env(1)%ot_pos(count)/gold_sec ! expand
            END IF
         ELSE
            ! first determine where we are and construct the new triplet
            IF (qs_ot_env(1)%ot_pos(count) < qs_ot_env(1)%ot_pos(qs_ot_env(1)%line_search_mid)) THEN
               IF (qs_ot_env(1)%ot_energy(count) < qs_ot_env(1)%ot_energy(qs_ot_env(1)%line_search_mid)) THEN
                  qs_ot_env(1)%line_search_right = qs_ot_env(1)%line_search_mid
                  qs_ot_env(1)%line_search_mid = count
               ELSE
                  qs_ot_env(1)%line_search_left = count
               END IF
            ELSE
               IF (qs_ot_env(1)%ot_energy(count) < qs_ot_env(1)%ot_energy(qs_ot_env(1)%line_search_mid)) THEN
                  qs_ot_env(1)%line_search_left = qs_ot_env(1)%line_search_mid
                  qs_ot_env(1)%line_search_mid = count
               ELSE
                  qs_ot_env(1)%line_search_right = count
               END IF
            END IF
            ! now find the new point in the largest section
            IF ((qs_ot_env(1)%ot_pos(qs_ot_env(1)%line_search_right) &
                 - qs_ot_env(1)%ot_pos(qs_ot_env(1)%line_search_mid)) > &
                (qs_ot_env(1)%ot_pos(qs_ot_env(1)%line_search_mid) &
                 - qs_ot_env(1)%ot_pos(qs_ot_env(1)%line_search_left))) THEN
               qs_ot_env(1)%ot_pos(count + 1) = &
                  qs_ot_env(1)%ot_pos(qs_ot_env(1)%line_search_mid) + &
                  gold_sec*(qs_ot_env(1)%ot_pos(qs_ot_env(1)%line_search_right) &
                            - qs_ot_env(1)%ot_pos(qs_ot_env(1)%line_search_mid))
            ELSE
               qs_ot_env(1)%ot_pos(count + 1) = &
                  qs_ot_env(1)%ot_pos(qs_ot_env(1)%line_search_left) + &
                  gold_sec*(qs_ot_env(1)%ot_pos(qs_ot_env(1)%line_search_mid) &
                            - qs_ot_env(1)%ot_pos(qs_ot_env(1)%line_search_left))
            END IF
            ! check for termination
            IF (((qs_ot_env(1)%ot_pos(qs_ot_env(1)%line_search_right) &
                  - qs_ot_env(1)%ot_pos(qs_ot_env(1)%line_search_mid)) < &
                 qs_ot_env(1)%ds_min*qs_ot_env(1)%settings%gold_target) .AND. &
                ((qs_ot_env(1)%ot_pos(qs_ot_env(1)%line_search_mid) &
                  - qs_ot_env(1)%ot_pos(qs_ot_env(1)%line_search_left)) < &
                 qs_ot_env(1)%ds_min*qs_ot_env(1)%settings%gold_target)) THEN
               qs_ot_env(1)%energy_only = .FALSE.
               qs_ot_env(1)%line_search_might_be_done = .TRUE.
            END IF
         END IF
      END IF
      ds = qs_ot_env(1)%OT_pos(count + 1) - qs_ot_env(1)%OT_pos(count)
      qs_ot_env(1)%ds_min = qs_ot_env(1)%OT_pos(count + 1)

      CALL take_step(ds, qs_ot_env)

      CALL timestop(handle)

   END SUBROUTINE do_line_search_gold

! **************************************************************************************************
!> \brief ...
!> \param qs_ot_env ...
! **************************************************************************************************
   SUBROUTINE do_line_search_3pnt(qs_ot_env)

      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: qs_ot_env

      CHARACTER(len=*), PARAMETER :: routineN = 'do_line_search_3pnt'
      REAL(KIND=dp), PARAMETER                           :: acceptance_tol = 1.0E-12_dp

      INTEGER                                            :: best_index, closure_backtracks, count, &
                                                            handle, i
      LOGICAL :: derivative_refinement, mermin_boundary_extend, mermin_closure_reset, &
         mermin_residual_closure, mermin_search, mermin_uphill_closure
      REAL(KIND=dp) :: best_drop, best_energy, closure_pos_floor, denom, &
         derivative_refinement_pos, ds, fa, fb, fc, model_drop, model_quality, nom, pos, &
         residual_growth_ratio, residual_quality, slope_mismatch, slope_model, tol, val, xa, xb, xc
      TYPE(cp_logger_type), POINTER                      :: logger

      CALL timeset(routineN, handle)

      qs_ot_env(1)%line_search_might_be_done = .FALSE.
      qs_ot_env(1)%energy_only = .TRUE.

      ! a three point interpolation based on the energy
      qs_ot_env(1)%line_search_count = qs_ot_env(1)%line_search_count + 1
      count = qs_ot_env(1)%line_search_count
      IF (count + 1 > SIZE(qs_ot_env(1)%OT_pos)) CPABORT("MAX ITER EXCEEDED : FATAL")
      qs_ot_env(1)%ot_energy(count) = qs_ot_env(1)%etotal
      SELECT CASE (count)
      CASE (1)
         qs_ot_env(1)%mermin_ls_closure_rejected = .FALSE.
         IF (kpoint_projected_tail_delta(qs_ot_env) > 0.0_dp .AND. &
             ot_kpoint_mermin_projected_reference_ready( &
             qs_ot_env(1)%kpoint_physical_delta_ready, &
             qs_ot_env(1)%kpoint_physical_delta_groups, &
             ot_kpoint_mermin_shared_spin_group( &
             qs_ot_env%spin_index, qs_ot_env%occupation_group_index))) THEN
            qs_ot_env(1)%mermin_ls_reference_delta = kpoint_projected_tail_delta(qs_ot_env)
         ELSE
            qs_ot_env(1)%mermin_ls_reference_delta = &
               -MAX(qs_ot_env(1)%delta, EPSILON(1.0_dp))
         END IF
         qs_ot_env(1)%ot_pos(count) = 0.0_dp
         qs_ot_env(1)%ot_grad(count) = qs_ot_env(1)%gradient
         qs_ot_env(1)%ot_pos(count + 1) = qs_ot_env(1)%ds_min*0.8_dp
      CASE (2)
         IF (qs_ot_env(1)%OT_energy(count) > qs_ot_env(1)%OT_energy(count - 1)) THEN
            qs_ot_env(1)%OT_pos(count + 1) = qs_ot_env(1)%ds_min*0.5_dp
         ELSE
            qs_ot_env(1)%OT_pos(count + 1) = qs_ot_env(1)%ds_min*1.4_dp
         END IF
      CASE (3)
         xa = qs_ot_env(1)%OT_pos(1)
         xb = qs_ot_env(1)%OT_pos(2)
         xc = qs_ot_env(1)%OT_pos(3)
         fa = qs_ot_env(1)%OT_energy(1)
         fb = qs_ot_env(1)%OT_energy(2)
         fc = qs_ot_env(1)%OT_energy(3)
         nom = (xb - xa)**2*(fb - fc) - (xb - xc)**2*(fb - fa)
         denom = (xb - xa)*(fb - fc) - (xb - xc)*(fb - fa)
         IF (ABS(denom) <= 1.0E-18_dp*MAX(ABS(fb - fc), ABS(fb - fa))) THEN
            pos = xb
         ELSE
            pos = xb - 0.5_dp*nom/denom ! position of the stationary point
         END IF
         val = (pos - xa)*(pos - xb)*fc/((xc - xa)*(xc - xb)) + &
               (pos - xb)*(pos - xc)*fa/((xa - xb)*(xa - xc)) + &
               (pos - xc)*(pos - xa)*fb/((xb - xc)*(xb - xa))
         slope_model = fa*(2.0_dp*xa - xb - xc)/((xa - xb)*(xa - xc)) + &
                       fb*(xa - xc)/((xb - xa)*(xb - xc)) + &
                       fc*(xa - xb)/((xc - xa)*(xc - xb))
         slope_mismatch = slope_model - qs_ot_env(1)%ot_grad(1)
         best_index = 1
         DO i = 2, count
            IF (qs_ot_env(1)%OT_energy(i) < qs_ot_env(1)%OT_energy(best_index)) best_index = i
         END DO
         best_energy = qs_ot_env(1)%OT_energy(best_index)
         model_drop = fa - val
         best_drop = fa - best_energy
         model_quality = 0.0_dp
         IF (model_drop > EPSILON(model_drop)) model_quality = best_drop/model_drop
         tol = acceptance_tol*MAX(1.0_dp, ABS(best_energy))
         mermin_search = use_three_point_mermin_search(qs_ot_env)
         mermin_boundary_extend = mermin_search .AND. best_index == count .AND. &
                                  fc < fb .AND. fb < fa .AND. val >= best_energy - tol
         IF (mermin_boundary_extend) THEN
            qs_ot_env(1)%OT_pos(count + 1) = MAXVAL(qs_ot_env(1)%OT_pos(1:3))*1.4_dp
         ELSE IF (mermin_search .AND. best_index == 1) THEN
            ! No sampled point descends from the current Mermin reference.
            ! Close the line search at the reference coordinates and verify
            ! that state once more after the density/Hamiltonian update.
            qs_ot_env(1)%OT_pos(count + 1) = qs_ot_env(1)%OT_pos(best_index)
         ELSE IF (mermin_search .AND. best_index /= 1 .AND. val >= best_energy - tol) THEN
            qs_ot_env(1)%OT_pos(count + 1) = qs_ot_env(1)%OT_pos(best_index)
         ELSE IF (val < fa .AND. val <= fb .AND. val <= fc) THEN ! OK, we go to a minimum
            ! we take a guard against too large steps
            qs_ot_env(1)%OT_pos(count + 1) = MAX(MAXVAL(qs_ot_env(1)%OT_pos(1:3))*0.01_dp, &
                                                 MIN(pos, MAXVAL(qs_ot_env(1)%OT_pos(1:3))*4.0_dp))
         ELSE ! just take an extended step
            qs_ot_env(1)%OT_pos(count + 1) = MAXVAL(qs_ot_env(1)%OT_pos(1:3))*2.0_dp
         END IF
         IF (mermin_search .AND. .NOT. mermin_boundary_extend) THEN
            qs_ot_env(1)%mermin_ls_sample_pending = .TRUE.
            qs_ot_env(1)%mermin_ls_best_sample_index = best_index
            qs_ot_env(1)%mermin_ls_reference_energy = fa
            qs_ot_env(1)%mermin_ls_best_sample_energy = best_energy
            qs_ot_env(1)%mermin_ls_best_sample_drop = best_drop
            qs_ot_env(1)%mermin_ls_model_drop = model_drop
            qs_ot_env(1)%mermin_ls_model_quality = model_quality
            qs_ot_env(1)%mermin_ls_best_sample_pos = qs_ot_env(1)%OT_pos(best_index)
            qs_ot_env(1)%mermin_ls_next_pos = qs_ot_env(1)%OT_pos(count + 1)
         END IF
         IF (mermin_search) THEN
            logger => cp_get_default_logger()
            IF (logger%iter_info%print_level >= debug_print_level .AND. logger%para_env%is_source()) THEN
               WRITE (cp_logger_get_default_unit_nr(logger), '(A,1X,I5,10(1X,ES16.8))') &
                  " OT Mermin LS model idx pos model best mdrop bdrop q slope0 fit mismatch next:", &
                  best_index, pos, val, best_energy, model_drop, best_drop, model_quality, &
                  qs_ot_env(1)%ot_grad(1), slope_model, slope_mismatch, qs_ot_env(1)%OT_pos(count + 1)
            END IF
         END IF
         qs_ot_env(1)%energy_only = mermin_boundary_extend
         qs_ot_env(1)%line_search_might_be_done = .NOT. mermin_boundary_extend
      CASE (4:)
         IF (.NOT. use_three_point_mermin_search(qs_ot_env)) CPABORT("NYI")
         tol = acceptance_tol*MAX(1.0_dp, ABS(qs_ot_env(1)%mermin_ls_reference_energy))
         mermin_uphill_closure = qs_ot_env(1)%mermin_ls_sample_pending .AND. &
                                 qs_ot_env(1)%OT_energy(count) > &
                                 qs_ot_env(1)%mermin_ls_reference_energy + tol
         mermin_residual_closure = kpoint_mermin_residual_closure_rejected( &
                                   qs_ot_env, residual_quality, residual_growth_ratio)
         mermin_closure_reset = .FALSE.
         derivative_refinement = .FALSE.
         IF (count == 4 .AND. qs_ot_env(1)%mermin_ls_sample_pending .AND. &
             .NOT. mermin_uphill_closure .AND. .NOT. mermin_residual_closure) THEN
            closure_pos_floor = MAX(1.0E4_dp*EPSILON(1.0_dp), &
                                    qs_ot_env(1)%settings%eps_taylor)* &
                                MAX(1.0_dp, ABS(qs_ot_env(1)%settings%ds_min))
            CALL ot_kpoint_mermin_derivative_refinement( &
               qs_ot_env(1)%ot_grad(1), qs_ot_env(1)%gradient, &
               qs_ot_env(1)%ot_pos(count), closure_pos_floor, &
               derivative_refinement_pos, derivative_refinement)
         END IF
         IF (derivative_refinement) THEN
            qs_ot_env(1)%OT_pos(count + 1) = derivative_refinement_pos
            qs_ot_env(1)%mermin_ls_next_pos = derivative_refinement_pos
            logger => cp_get_default_logger()
            IF (logger%iter_info%print_level >= debug_print_level .AND. &
                logger%para_env%is_source()) THEN
               WRITE (cp_logger_get_default_unit_nr(logger), '(A,5(1X,ES16.8))') &
                  " OT Mermin LS derivative refinement reference current old new floor:", &
                  qs_ot_env(1)%ot_grad(1), qs_ot_env(1)%gradient, &
                  qs_ot_env(1)%ot_pos(count), derivative_refinement_pos, closure_pos_floor
            END IF
         ELSE IF (mermin_residual_closure) THEN
            closure_backtracks = count - 3
            closure_pos_floor = MAX(1.0E4_dp*EPSILON(1.0_dp), &
                                    qs_ot_env(1)%settings%eps_taylor)* &
                                MAX(1.0_dp, ABS(qs_ot_env(1)%settings%ds_min))
            mermin_closure_reset = closure_backtracks >= ot_cg_kp_mermin_closure_max_backtracks .OR. &
                                   ABS(qs_ot_env(1)%OT_pos(count)) <= closure_pos_floor
            IF (mermin_closure_reset) THEN
               qs_ot_env(1)%OT_pos(count + 1) = 0.0_dp
               qs_ot_env(1)%mermin_ls_sample_pending = .FALSE.
               qs_ot_env(1)%mermin_ls_closure_rejected = .TRUE.
               qs_ot_env(1)%schur_residual_history_reset_pending = .TRUE.
               qs_ot_env(1)%schur_residual_history_quality = 0.0_dp
               ! The accepted candidate model has exhausted its residual trust
               ! region even after repeated affine backtracking.  Retain the
               ! physical secants and coordinates, but relearn the candidate
               ! calibration together with the fresh steepest-descent history.
               CALL qs_ot_reset_accepted_model_calibration(qs_ot_env(1))
            ELSE
               qs_ot_env(1)%OT_pos(count + 1) = &
                  ot_cg_kp_mermin_closure_backtrack*qs_ot_env(1)%OT_pos(count)
               qs_ot_env(1)%mermin_ls_next_pos = qs_ot_env(1)%OT_pos(count + 1)
               qs_ot_env(1)%schur_residual_history_reset_pending = .TRUE.
               qs_ot_env(1)%schur_residual_history_quality = &
                  MIN(qs_ot_env(1)%schur_residual_history_quality, residual_quality)
            END IF
            logger => cp_get_default_logger()
            IF (logger%iter_info%print_level >= debug_print_level .AND. logger%para_env%is_source()) THEN
               WRITE (cp_logger_get_default_unit_nr(logger), '(A,1X,I5,1X,L1,7(1X,ES16.8))') &
                  " OT Mermin LS residual backtrack count reset ref current ratio quality pos next floor:", &
                  closure_backtracks, mermin_closure_reset, &
                  qs_ot_env(1)%mermin_ls_reference_energy, qs_ot_env(1)%OT_energy(count), &
                  residual_growth_ratio, residual_quality, qs_ot_env(1)%OT_pos(count), &
                  qs_ot_env(1)%OT_pos(count + 1), closure_pos_floor
            END IF
         ELSE IF (qs_ot_env(1)%mermin_ls_sample_pending .AND. .NOT. mermin_uphill_closure) THEN
            best_index = count
            best_energy = qs_ot_env(1)%OT_energy(count)
            best_drop = qs_ot_env(1)%mermin_ls_reference_energy - best_energy
            qs_ot_env(1)%OT_pos(count + 1) = qs_ot_env(1)%OT_pos(count)
            qs_ot_env(1)%mermin_ls_best_sample_index = best_index
            qs_ot_env(1)%mermin_ls_best_sample_energy = best_energy
            qs_ot_env(1)%mermin_ls_best_sample_drop = MAX(0.0_dp, best_drop)
            IF (qs_ot_env(1)%mermin_ls_model_drop > EPSILON(1.0_dp)) THEN
               qs_ot_env(1)%mermin_ls_model_quality = &
                  MAX(0.0_dp, best_drop)/qs_ot_env(1)%mermin_ls_model_drop
            END IF
            qs_ot_env(1)%mermin_ls_best_sample_pos = qs_ot_env(1)%OT_pos(count)
            qs_ot_env(1)%mermin_ls_next_pos = qs_ot_env(1)%OT_pos(count)
         ELSE IF (mermin_uphill_closure) THEN
            closure_backtracks = count - 3
            closure_pos_floor = MAX(1.0E4_dp*EPSILON(1.0_dp), &
                                    qs_ot_env(1)%settings%eps_taylor)* &
                                MAX(1.0_dp, ABS(qs_ot_env(1)%settings%ds_min))
            mermin_closure_reset = closure_backtracks >= ot_cg_kp_mermin_closure_max_backtracks .OR. &
                                   ABS(qs_ot_env(1)%OT_pos(count)) <= closure_pos_floor
            IF (mermin_closure_reset) THEN
               qs_ot_env(1)%OT_pos(count + 1) = 0.0_dp
               qs_ot_env(1)%mermin_ls_sample_pending = .FALSE.
               qs_ot_env(1)%mermin_ls_closure_rejected = .TRUE.
               qs_ot_env(1)%schur_residual_history_reset_pending = .TRUE.
               qs_ot_env(1)%schur_residual_history_quality = 0.0_dp
            ELSE
               qs_ot_env(1)%OT_pos(count + 1) = &
                  ot_cg_kp_mermin_closure_backtrack*qs_ot_env(1)%OT_pos(count)
               qs_ot_env(1)%mermin_ls_next_pos = qs_ot_env(1)%OT_pos(count + 1)
            END IF
            logger => cp_get_default_logger()
            IF (logger%iter_info%print_level >= debug_print_level .AND. logger%para_env%is_source()) THEN
               WRITE (cp_logger_get_default_unit_nr(logger), '(A,1X,I5,1X,L1,5(1X,ES16.8))') &
                  " OT Mermin LS closure backtrack count reset ref current pos next:", &
                  closure_backtracks, mermin_closure_reset, &
                  qs_ot_env(1)%mermin_ls_reference_energy, qs_ot_env(1)%OT_energy(count), &
                  qs_ot_env(1)%OT_pos(count), qs_ot_env(1)%OT_pos(count + 1), closure_pos_floor
            END IF
         ELSE IF (count == 4 .AND. qs_ot_env(1)%OT_pos(4) > qs_ot_env(1)%OT_pos(3) .AND. &
                  qs_ot_env(1)%OT_energy(3) < qs_ot_env(1)%OT_energy(2) .AND. &
                  qs_ot_env(1)%OT_energy(2) < qs_ot_env(1)%OT_energy(1)) THEN
            best_index = 3
            IF (qs_ot_env(1)%OT_energy(4) < qs_ot_env(1)%OT_energy(3)) best_index = 4
            best_energy = qs_ot_env(1)%OT_energy(best_index)
            best_drop = qs_ot_env(1)%OT_energy(1) - best_energy
            qs_ot_env(1)%OT_pos(count + 1) = qs_ot_env(1)%OT_pos(best_index)
            IF (best_index /= 1) THEN
               qs_ot_env(1)%mermin_ls_sample_pending = .TRUE.
               qs_ot_env(1)%mermin_ls_best_sample_index = best_index
               qs_ot_env(1)%mermin_ls_reference_energy = qs_ot_env(1)%OT_energy(1)
               qs_ot_env(1)%mermin_ls_best_sample_energy = best_energy
               qs_ot_env(1)%mermin_ls_best_sample_drop = best_drop
               qs_ot_env(1)%mermin_ls_model_drop = best_drop
               qs_ot_env(1)%mermin_ls_model_quality = 1.0_dp
               qs_ot_env(1)%mermin_ls_best_sample_pos = qs_ot_env(1)%OT_pos(best_index)
               qs_ot_env(1)%mermin_ls_next_pos = qs_ot_env(1)%OT_pos(count + 1)
            END IF
         ELSE
            ! Restarted or legacy Mermin 3-point histories can already be past the
            ! endpoint.  Accept the best sampled coordinates rather than a
            ! potentially uphill current point.
            best_index = 1
            DO i = 2, count
               IF (qs_ot_env(1)%OT_energy(i) < qs_ot_env(1)%OT_energy(best_index)) THEN
                  best_index = i
               END IF
            END DO
            best_energy = qs_ot_env(1)%OT_energy(best_index)
            best_drop = qs_ot_env(1)%OT_energy(1) - best_energy
            qs_ot_env(1)%OT_pos(count + 1) = qs_ot_env(1)%OT_pos(best_index)
            IF (best_index /= 1) THEN
               qs_ot_env(1)%mermin_ls_sample_pending = .TRUE.
               qs_ot_env(1)%mermin_ls_best_sample_index = best_index
               qs_ot_env(1)%mermin_ls_reference_energy = qs_ot_env(1)%OT_energy(1)
               qs_ot_env(1)%mermin_ls_best_sample_energy = best_energy
               qs_ot_env(1)%mermin_ls_best_sample_drop = best_drop
               qs_ot_env(1)%mermin_ls_model_drop = best_drop
               qs_ot_env(1)%mermin_ls_model_quality = 1.0_dp
               qs_ot_env(1)%mermin_ls_best_sample_pos = qs_ot_env(1)%OT_pos(best_index)
               qs_ot_env(1)%mermin_ls_next_pos = qs_ot_env(1)%OT_pos(count + 1)
            ELSE
               qs_ot_env(1)%mermin_ls_sample_pending = .FALSE.
            END IF
         END IF
         IF (derivative_refinement) THEN
            qs_ot_env(1)%energy_only = .FALSE.
            qs_ot_env(1)%line_search_might_be_done = .TRUE.
         ELSE IF (mermin_residual_closure) THEN
            ! Residual backtracking needs a gradient at every candidate endpoint.
            qs_ot_env(1)%energy_only = .FALSE.
            qs_ot_env(1)%line_search_might_be_done = mermin_closure_reset
         ELSE
            qs_ot_env(1)%energy_only = mermin_uphill_closure .AND. .NOT. mermin_closure_reset
            qs_ot_env(1)%line_search_might_be_done = .NOT. mermin_uphill_closure .OR. mermin_closure_reset
         END IF
      CASE DEFAULT
         CPABORT("NYI")
      END SELECT
      ds = qs_ot_env(1)%OT_pos(count + 1) - qs_ot_env(1)%OT_pos(count)
      qs_ot_env(1)%ds_min = qs_ot_env(1)%OT_pos(count + 1)

      CALL take_step(ds, qs_ot_env)

      CALL timestop(handle)

   END SUBROUTINE do_line_search_3pnt

! **************************************************************************************************
!> \brief ...
!> \param qs_ot_env ...
! **************************************************************************************************
   SUBROUTINE do_line_search_2pnt(qs_ot_env)

      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: qs_ot_env

      CHARACTER(len=*), PARAMETER :: routineN = 'do_line_search_2pnt'

      INTEGER                                            :: count, handle
      REAL(KIND=dp)                                      :: a, b, c, ds, pos, val, x0, x1

      CALL timeset(routineN, handle)

      qs_ot_env(1)%line_search_might_be_done = .FALSE.
      qs_ot_env(1)%energy_only = .TRUE.

      ! a three point interpolation based on the energy
      qs_ot_env(1)%line_search_count = qs_ot_env(1)%line_search_count + 1
      count = qs_ot_env(1)%line_search_count
      qs_ot_env(1)%ot_energy(count) = qs_ot_env(1)%etotal
      SELECT CASE (count)
      CASE (1)
         qs_ot_env(1)%ot_pos(count) = 0.0_dp
         qs_ot_env(1)%ot_grad(count) = qs_ot_env(1)%gradient
         qs_ot_env(1)%ot_pos(count + 1) = qs_ot_env(1)%ds_min*1.0_dp
      CASE (2)
         x0 = 0.0_dp
         c = qs_ot_env(1)%ot_energy(1)
         b = qs_ot_env(1)%ot_grad(1)
         x1 = qs_ot_env(1)%ot_pos(2)
         a = (qs_ot_env(1)%ot_energy(2) - b*x1 - c)/(x1**2)
         IF (a <= 0.0_dp) a = 1.0E-15_dp
         pos = -b/(2.0_dp*a)
         val = a*pos**2 + b*pos + c
         qs_ot_env(1)%energy_only = .FALSE.
         qs_ot_env(1)%line_search_might_be_done = .TRUE.
         IF (val < qs_ot_env(1)%ot_energy(1) .AND. val <= qs_ot_env(1)%ot_energy(2)) THEN
            ! we go to a minimum, but ...
            ! we take a guard against too large steps
            qs_ot_env(1)%OT_pos(count + 1) = MAX(MAXVAL(qs_ot_env(1)%OT_pos(1:2))*0.01_dp, &
                                                 MIN(pos, MAXVAL(qs_ot_env(1)%OT_pos(1:2))*4.0_dp))
         ELSE ! just take an extended step
            qs_ot_env(1)%OT_pos(count + 1) = MAXVAL(qs_ot_env(1)%OT_pos(1:2))*2.0_dp
         END IF
      CASE DEFAULT
         CPABORT("NYI")
      END SELECT
      ds = qs_ot_env(1)%OT_pos(count + 1) - qs_ot_env(1)%OT_pos(count)
      qs_ot_env(1)%ds_min = qs_ot_env(1)%OT_pos(count + 1)

      CALL take_step(ds, qs_ot_env)

      CALL timestop(handle)

   END SUBROUTINE do_line_search_2pnt

! **************************************************************************************************
!> \brief ...
!> \param qs_ot_env ...
! **************************************************************************************************
   SUBROUTINE do_line_search_adapt(qs_ot_env)

      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: qs_ot_env

      CHARACTER(len=*), PARAMETER :: routineN = 'do_line_search_adapt'
      REAL(KIND=dp), PARAMETER                           :: grow_factor = 2.0_dp, &
                                                            shrink_factor = 0.5_dp

      INTEGER                                            :: count, handle, il, im, ir
      REAL(KIND=dp)                                      :: a, b, c, denom, ds, el, em, er, &
                                                            step_size, xl, xm, xr

      CALL timeset(routineN, handle)

      qs_ot_env(1)%line_search_count = qs_ot_env(1)%line_search_count + 1
      count = qs_ot_env(1)%line_search_count
      qs_ot_env(1)%line_search_might_be_done = .FALSE.
      qs_ot_env(1)%energy_only = .TRUE.

      IF (count + 1 > SIZE(qs_ot_env(1)%OT_pos)) THEN
         ! should not happen, we pass with a warning first
         ! you can increase the size of OT_pos and the like in qs_ot_env
         CPABORT("MAX ITER EXCEEDED : FATAL")
      END IF

      ! Perform an adaptive linesearch
      IF (qs_ot_env(1)%line_search_count == 1) THEN
         qs_ot_env(1)%line_search_left = 1
         qs_ot_env(1)%line_search_right = 0
         qs_ot_env(1)%line_search_mid = 1
         qs_ot_env(1)%ot_Pos(1) = 0.0_dp
         qs_ot_env(1)%ot_energy(1) = qs_ot_env(1)%etotal
         qs_ot_env(1)%ot_Pos(2) = qs_ot_env(1)%ds_min*grow_factor
      ELSE
         qs_ot_env(1)%ot_energy(count) = qs_ot_env(1)%etotal
         ! it's essentially a book keeping game.
         ! keep left on the left, keep (bring) right on the right
         ! and mid in between these two
         IF (qs_ot_env(1)%line_search_right == 0) THEN ! we do not yet have the right bracket
            IF (qs_ot_env(1)%ot_energy(count - 1) < qs_ot_env(1)%ot_energy(count)) THEN
               qs_ot_env(1)%line_search_right = count
               qs_ot_env(1)%ot_Pos(count + 1) = qs_ot_env(1)%ot_Pos(qs_ot_env(1)%line_search_mid) + &
                                                (qs_ot_env(1)%ot_Pos(qs_ot_env(1)%line_search_right) - &
                                                 qs_ot_env(1)%ot_Pos(qs_ot_env(1)%line_search_mid))*shrink_factor
            ELSE
               ! expand further
               qs_ot_env(1)%line_search_left = qs_ot_env(1)%line_search_mid
               qs_ot_env(1)%line_search_mid = count
               qs_ot_env(1)%ot_Pos(count + 1) = qs_ot_env(1)%ot_Pos(count)*grow_factor
            END IF
         ELSE
            ! first determine where we are and construct the new triplet
            IF (qs_ot_env(1)%ot_pos(count) < qs_ot_env(1)%ot_pos(qs_ot_env(1)%line_search_mid)) THEN
               IF (qs_ot_env(1)%ot_energy(count) < qs_ot_env(1)%ot_energy(qs_ot_env(1)%line_search_mid)) THEN
                  qs_ot_env(1)%line_search_right = qs_ot_env(1)%line_search_mid
                  qs_ot_env(1)%line_search_mid = count
               ELSE
                  qs_ot_env(1)%line_search_left = count
               END IF
            ELSE
               IF (qs_ot_env(1)%ot_energy(count) < qs_ot_env(1)%ot_energy(qs_ot_env(1)%line_search_mid)) THEN
                  qs_ot_env(1)%line_search_left = qs_ot_env(1)%line_search_mid
                  qs_ot_env(1)%line_search_mid = count
               ELSE
                  qs_ot_env(1)%line_search_right = count
               END IF
            END IF
            il = qs_ot_env(1)%line_search_left
            im = qs_ot_env(1)%line_search_mid
            ir = qs_ot_env(1)%line_search_right
            xl = qs_ot_env(1)%OT_pos(il)
            xm = qs_ot_env(1)%OT_pos(im)
            xr = qs_ot_env(1)%OT_pos(ir)
            el = qs_ot_env(1)%ot_energy(il)
            em = qs_ot_env(1)%ot_energy(im)
            er = qs_ot_env(1)%ot_energy(ir)
            IF (em < el) THEN
               IF (er < em) THEN
                  !extend search
                  qs_ot_env(1)%ot_Pos(count + 1) = qs_ot_env(1)%ot_Pos(ir)*grow_factor
               ELSE
                  ! Cramer's rule
                  denom = (xl - xm)*(xl - xr)*(xm - xr)
                  a = (xr*(em - el) + xm*(el - er) + xl*(er - em))/denom
                  b = (xr**2*(el - em) + xm**2*(er - el) + xl**2*(em - er))/denom
                  c = (xm*xr*(xm - xr)*el + xr*xl*(xr - xl)*em + xr*xm*(xr - xm)*er)/denom

                  IF (ABS(a) /= 0.0_dp) THEN
                     step_size = -b/(2.0_dp*a)
                  ELSE
                     step_size = 0.0_dp
                  END IF
                  CPASSERT(step_size >= 0.0_dp)
                  qs_ot_env(1)%ot_Pos(count + 1) = step_size
                  qs_ot_env(1)%line_search_might_be_done = .TRUE.
                  qs_ot_env(1)%energy_only = .FALSE.
               END IF
            ELSE
               ! contract search
               qs_ot_env(1)%ot_Pos(count + 1) = qs_ot_env(1)%ot_Pos(im) + &
                                                (qs_ot_env(1)%ot_Pos(ir) - qs_ot_env(1)%ot_Pos(im))*shrink_factor
            END IF

         END IF
      END IF
      ds = qs_ot_env(1)%OT_pos(count + 1) - qs_ot_env(1)%OT_pos(count)
      qs_ot_env(1)%ds_min = qs_ot_env(1)%OT_pos(count + 1)

      CALL take_step(ds, qs_ot_env)

      CALL timestop(handle)

   END SUBROUTINE do_line_search_adapt

! **************************************************************************************************
!> \brief ...
!> \param qs_ot_env ...
! **************************************************************************************************
   SUBROUTINE do_line_search_none(qs_ot_env)
      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: qs_ot_env

      CALL take_step(qs_ot_env(1)%ds_min, qs_ot_env)

   END SUBROUTINE do_line_search_none

!
! creates a new SD direction, using the preconditioner if associated
! also updates the gradient for line search
!

! **************************************************************************************************
!> \brief ...
!> \param qs_ot_env ...
!> \param para_env_inter_kp ...
! **************************************************************************************************
   SUBROUTINE ot_new_sd_direction(qs_ot_env, para_env_inter_kp)
      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: qs_ot_env
      TYPE(mp_para_env_type), OPTIONAL, POINTER          :: para_env_inter_kp

      CHARACTER(len=*), PARAMETER :: routineN = 'ot_new_sd_direction'

      INTEGER                                            :: handle, ispin, nspin
      LOGICAL                                            :: do_ener, do_ks
      REAL(KIND=dp)                                      :: gnorm_metric, norm_weight, nvar, tmp
      TYPE(cp_logger_type), POINTER                      :: logger

      CALL timeset(routineN, handle)

!***SCP

      nspin = SIZE(qs_ot_env)
      logger => cp_get_default_logger()
      do_ks = qs_ot_env(1)%settings%ks
      do_ener = qs_ot_env(1)%settings%do_ener
      gnorm_metric = 0.0_dp

      IF (ASSOCIATED(qs_ot_env(1)%preconditioner)) THEN
         IF (.NOT. qs_ot_env(1)%use_dx) CPABORT("use dx")
         qs_ot_env(1)%gnorm = 0.0_dp
         IF (do_ks) THEN
            DO ispin = 1, nspin
               norm_weight = ot_kpoint_channel_norm_weight(qs_ot_env(ispin))
               IF (qs_ot_env(ispin)%has_complex_kpoint_state .AND. &
                   ASSOCIATED(qs_ot_env(ispin)%preconditioner%complex_fm)) THEN
                  CALL apply_preconditioner(qs_ot_env(ispin)%preconditioner, &
                                            qs_ot_env(ispin)%matrix_gx, &
                                            qs_ot_env(ispin)%matrix_gx_im, &
                                            qs_ot_env(ispin)%matrix_dx, &
                                            qs_ot_env(ispin)%matrix_dx_im)
               ELSE
                  CALL apply_preconditioner(qs_ot_env(ispin)%preconditioner, &
                                            qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin)%matrix_dx)
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL apply_preconditioner(qs_ot_env(ispin)%preconditioner, &
                                               qs_ot_env(ispin)%matrix_gx_im, &
                                               qs_ot_env(ispin)%matrix_dx_im)
                  END IF
               END IF
               CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin)%matrix_dx, tmp)
               qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + tmp
               gnorm_metric = gnorm_metric + norm_weight*tmp
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx_im, qs_ot_env(ispin)%matrix_dx_im, tmp)
                  qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + tmp
                  gnorm_metric = gnorm_metric + norm_weight*tmp
               END IF
            END DO
            IF (qs_ot_env(1)%gnorm < 0.0_dp) THEN
               logger => cp_get_default_logger()
               WRITE (cp_logger_get_default_unit_nr(logger), *) "WARNING Preconditioner not positive definite !"
            END IF
            DO ispin = 1, nspin
               CALL dbcsr_scale(qs_ot_env(ispin)%matrix_dx, -1.0_dp)
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  CALL dbcsr_scale(qs_ot_env(ispin)%matrix_dx_im, -1.0_dp)
               END IF
            END DO
            IF (qs_ot_env(1)%settings%do_rotation) THEN
               DO ispin = 1, nspin
                  norm_weight = ot_kpoint_channel_norm_weight(qs_ot_env(ispin))
                  ! right now no preconditioner yet
                  CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_dx, qs_ot_env(ispin)%rot_mat_gx)
                  CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx, qs_ot_env(ispin)%rot_mat_dx, tmp)
                  ! added 0.5, because we have (antisymmetry) only half the number of variables
                  qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + 0.5_dp*tmp
                  gnorm_metric = gnorm_metric + norm_weight*0.5_dp*tmp
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_dx_im, qs_ot_env(ispin)%rot_mat_gx_im)
                     CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx_im, qs_ot_env(ispin)%rot_mat_dx_im, tmp)
                     qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + 0.5_dp*tmp
                     gnorm_metric = gnorm_metric + norm_weight*0.5_dp*tmp
                  END IF
               END DO
               DO ispin = 1, nspin
                  CALL dbcsr_scale(qs_ot_env(ispin)%rot_mat_dx, -1.0_dp)
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_scale(qs_ot_env(ispin)%rot_mat_dx_im, -1.0_dp)
                  END IF
               END DO
            END IF
         END IF
         IF (do_ener) THEN
            DO ispin = 1, nspin
               norm_weight = ot_kpoint_channel_norm_weight(qs_ot_env(ispin))
               qs_ot_env(ispin)%ener_dx = qs_ot_env(ispin)%ener_gx
               tmp = DOT_PRODUCT(qs_ot_env(ispin)%ener_dx, qs_ot_env(ispin)%ener_gx)
               qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + tmp
               gnorm_metric = gnorm_metric + norm_weight*tmp
               qs_ot_env(ispin)%ener_dx = -qs_ot_env(ispin)%ener_dx
            END DO
         END IF
      ELSE
         qs_ot_env(1)%gnorm = 0.0_dp
         IF (do_ks) THEN
            DO ispin = 1, nspin
               norm_weight = ot_kpoint_channel_norm_weight(qs_ot_env(ispin))
               CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin)%matrix_gx, tmp)
               qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + tmp
               gnorm_metric = gnorm_metric + norm_weight*tmp
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx_im, qs_ot_env(ispin)%matrix_gx_im, tmp)
                  qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + tmp
                  gnorm_metric = gnorm_metric + norm_weight*tmp
               END IF
            END DO
            IF (qs_ot_env(1)%settings%do_rotation) THEN
               DO ispin = 1, nspin
                  norm_weight = ot_kpoint_channel_norm_weight(qs_ot_env(ispin))
                  CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx, qs_ot_env(ispin)%rot_mat_gx, tmp)
                  ! added 0.5, because we have (antisymmetry) only half the number of variables
                  qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + 0.5_dp*tmp
                  gnorm_metric = gnorm_metric + norm_weight*0.5_dp*tmp
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx_im, qs_ot_env(ispin)%rot_mat_gx_im, tmp)
                     qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + 0.5_dp*tmp
                     gnorm_metric = gnorm_metric + norm_weight*0.5_dp*tmp
                  END IF
               END DO
            END IF
         END IF
         IF (do_ener) THEN
            DO ispin = 1, nspin
               norm_weight = ot_kpoint_channel_norm_weight(qs_ot_env(ispin))
               tmp = DOT_PRODUCT(qs_ot_env(ispin)%ener_gx, qs_ot_env(ispin)%ener_gx)
               qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + tmp
               gnorm_metric = gnorm_metric + norm_weight*tmp
            END DO
         END IF
      END IF

      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, qs_ot_env(1)%gnorm)
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, gnorm_metric)
      nvar = ot_kpoint_channel_nvar(qs_ot_env, do_ks, do_ener)
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, nvar)
      ! Handling the case of no free variables to optimize
      IF (nvar > 0.0_dp) THEN
         qs_ot_env(1)%delta = SQRT(ABS(gnorm_metric)/nvar)
         qs_ot_env(1)%gradient = -qs_ot_env(1)%gnorm
      ELSE
         qs_ot_env(1)%delta = 0.0_dp
         qs_ot_env(1)%gradient = 0.0_dp
      END IF

      CALL timestop(handle)

   END SUBROUTINE ot_new_sd_direction

!
! creates a new CG direction. Implements Polak-Ribierre variant
! using the preconditioner if associated
! also updates the gradient for line search
!
! **************************************************************************************************
!> \brief ...
!> \param qs_ot_env ...
!> \param para_env_inter_kp ...
!> \param residual_closure_rejected retain the current line-search direction for residual backtracking
!> \param residual_preflight_only update the current residual without replacing the search direction
! **************************************************************************************************
   SUBROUTINE ot_new_cg_direction(qs_ot_env, para_env_inter_kp, residual_closure_rejected, &
                                  residual_preflight_only)
      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: qs_ot_env
      TYPE(mp_para_env_type), OPTIONAL, POINTER          :: para_env_inter_kp
      LOGICAL, INTENT(OUT), OPTIONAL                     :: residual_closure_rejected
      LOGICAL, INTENT(IN), OPTIONAL                      :: residual_preflight_only

      CHARACTER(len=*), PARAMETER :: routineN = 'ot_new_cg_direction'

      INTEGER :: catastrophic_reset_any, direction_ready_all, growth_reset_any, handle, ispin, &
         matrix_cols, matrix_rows, nspin, schur_attempt, schur_flavor, schur_last_attempt, &
         schur_last_flavor
      LOGICAL :: complex_kpoint_history, debug_diagnostic, debug_diagnostic_source, do_ener, &
         do_ks, history_spike_reset, kpoint_preconditioned_direction, &
         legacy_preconditioned_gradient, mermin_residual_catastrophic_reset, &
         mermin_residual_growth_reset, negative_beta, physical_raw_gap_guard, preflight_only, &
         projected_raw_tail_release, reset_to_sd, residual_history_reset, &
         residual_spike_history_reset, residual_tail_spike_reset, schur_calibration_change, &
         schur_history_change, schur_history_mixed, schur_projector_jump_reset, &
         schur_projector_tail_reset
      REAL(KIND=dp) :: beta_pr, beta_raw, full_direction_drop, full_direction_optimum, &
         full_direction_quadratic, full_direction_slope, gnorm_cross, gnorm_ener, gnorm_matrix_im, &
         gnorm_matrix_re, gnorm_metric, gnorm_ratio, gnorm_rot_im, gnorm_rot_re, &
         matrix_gnorm_metric, matrix_nvar, mermin_residual_growth_quality, mermin_tail_min, &
         norm_weight, nvar, occ_scale_avg, occ_scale_max, occ_scale_min, &
         partial_direction_quadratic, partial_direction_slope, physical_control_delta, &
         physical_delta, physical_gnorm_metric, physical_tail_delta, previous_delta, &
         raw_channel_rms, raw_delta, raw_guard_delta
      REAL(KIND=dp) :: raw_matrix_gnorm_im, raw_matrix_gnorm_metric_im, &
         raw_matrix_gnorm_metric_re, raw_matrix_gnorm_re, rayleigh_model_quadratic_factor, &
         rayleigh_model_slope_factor, residual_growth_delta, residual_spike_quality, &
         residual_tail_min, residual_trust_quality, residual_trust_scale, sample_accepted_drop, &
         sample_gap, sample_guard_drop, sample_quality, sample_residual_growth_ratio, &
         search_channel_rms, search_direction_rms, test_down, tmp
      TYPE(ot_cg_history_snapshot_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: history_snapshot
      TYPE(cp_logger_type), POINTER                      :: logger
      REAL(KIND=dp), DIMENSION(SIZE(qs_ot_env))          :: channel_direction_metric, &
                                                            channel_raw_metric, &
                                                            channel_search_metric

      CALL timeset(routineN, handle)

      IF (PRESENT(residual_closure_rejected)) residual_closure_rejected = .FALSE.
      preflight_only = .FALSE.
      IF (PRESENT(residual_preflight_only)) preflight_only = residual_preflight_only

      nspin = SIZE(qs_ot_env)
      logger => cp_get_default_logger()
      debug_diagnostic = logger%iter_info%print_level >= debug_print_level
      debug_diagnostic_source = debug_diagnostic .AND. logger%para_env%is_source()

      do_ks = qs_ot_env(1)%settings%ks
      do_ener = qs_ot_env(1)%settings%do_ener
      IF (preflight_only) CALL capture_ot_cg_history(qs_ot_env, history_snapshot)
      complex_kpoint_history = .FALSE.
      kpoint_preconditioned_direction = .FALSE.
      rayleigh_model_slope_factor = 1.0_dp
      rayleigh_model_quadratic_factor = 1.0_dp
      gnorm_metric = 0.0_dp
      DO ispin = 1, nspin
         complex_kpoint_history = complex_kpoint_history .OR. qs_ot_env(ispin)%has_complex_kpoint_state
      END DO
      legacy_preconditioned_gradient = .NOT. complex_kpoint_history .AND. &
                                       .NOT. qs_ot_env(1)%settings%do_rotation .AND. .NOT. do_ener
      IF (complex_kpoint_history) THEN
         kpoint_preconditioned_direction = .TRUE.
         DO ispin = 1, nspin
            IF (qs_ot_env(ispin)%has_complex_kpoint_state .AND. &
                .NOT. qs_ot_env(ispin)%kpoint_preconditioned_direction_ready) THEN
               kpoint_preconditioned_direction = .FALSE.
            END IF
         END DO
      END IF
      direction_ready_all = MERGE(1, 0, kpoint_preconditioned_direction)
      CALL ot_mini_min_kpoint_integer(qs_ot_env, para_env_inter_kp, direction_ready_all)
      kpoint_preconditioned_direction = direction_ready_all == 1
      history_spike_reset = .FALSE.
      negative_beta = .FALSE.
      mermin_residual_catastrophic_reset = .FALSE.
      mermin_residual_growth_reset = .FALSE.
      mermin_residual_growth_quality = 1.0_dp
      residual_history_reset = .FALSE.
      residual_spike_history_reset = .FALSE.
      residual_spike_quality = 1.0_dp
      residual_tail_spike_reset = .FALSE.
      residual_trust_quality = 1.0_dp
      mermin_tail_min = kpoint_eps_scaled_tail_min( &
                        ot_cg_kp_mermin_recovery_tail_min_delta, &
                        ot_cg_kp_mermin_recovery_tail_min_eps_factor, &
                        qs_ot_env(1)%kpoint_eps_scf)
      residual_tail_min = kpoint_eps_scaled_tail_min( &
                          ot_cg_kp_residual_tail_min_delta, &
                          ot_cg_kp_residual_tail_min_eps_factor, &
                          qs_ot_env(1)%kpoint_eps_scf)
      reset_to_sd = .FALSE.
      schur_attempt = qs_ot_env(1)%schur_model_attempt
      schur_flavor = qs_ot_env(1)%schur_model_flavor
      schur_last_attempt = qs_ot_env(1)%schur_cg_attempt
      schur_last_flavor = qs_ot_env(1)%schur_cg_flavor
      ! Expose Split/Dense Schur model changes next to CG beta.  Dense
      ! sub-flavor changes are diagnosed here.  Explicitly mixed gauge windows
      ! and activation/deactivation of the accepted-secant Hessian reset CG;
      ! ordinary calibrated candidate changes keep their accumulated history.
      schur_history_change = qs_ot_occupation_preconditioner_active(qs_ot_env(1)) .AND. &
                             qs_ot_env(1)%settings%do_rotation .AND. &
                             (((schur_attempt /= schur_last_attempt) .AND. &
                               schur_attempt /= 0 .AND. schur_last_attempt /= 0 .AND. &
                               (schur_attempt <= 2 .OR. schur_last_attempt <= 2)) .OR. &
                              (schur_attempt == 3 .AND. schur_last_attempt == 3 .AND. &
                               schur_flavor /= schur_last_flavor .AND. &
                               schur_flavor /= 0 .AND. schur_last_flavor /= 0 .AND. &
                               (schur_flavor < 0 .OR. schur_last_flavor < 0 .OR. &
                                schur_flavor == 42 .OR. schur_last_flavor == 42 .OR. &
                                schur_flavor == 43 .OR. schur_last_flavor == 43 .OR. &
                                schur_flavor == 44 .OR. schur_last_flavor == 44 .OR. &
                                schur_flavor == 45 .OR. schur_last_flavor == 45 .OR. &
                                schur_flavor == 46 .OR. schur_last_flavor == 46 .OR. &
                                schur_flavor == 47 .OR. schur_last_flavor == 47 .OR. &
                                schur_flavor == 48 .OR. schur_last_flavor == 48)))
      schur_calibration_change = schur_history_change .AND. &
                                 (schur_flavor == 43 .OR. schur_last_flavor == 43 .OR. &
                                  schur_flavor == 44 .OR. schur_last_flavor == 44 .OR. &
                                  schur_flavor == 46 .OR. schur_last_flavor == 46 .OR. &
                                  schur_flavor == 47 .OR. schur_last_flavor == 47 .OR. &
                                  schur_flavor == 48 .OR. schur_last_flavor == 48)
      schur_history_mixed = .FALSE.

      gnorm_matrix_re = 0.0_dp
      gnorm_matrix_im = 0.0_dp
      raw_matrix_gnorm_re = 0.0_dp
      raw_matrix_gnorm_im = 0.0_dp
      raw_matrix_gnorm_metric_re = 0.0_dp
      raw_matrix_gnorm_metric_im = 0.0_dp
      gnorm_rot_re = 0.0_dp
      gnorm_rot_im = 0.0_dp
      gnorm_ener = 0.0_dp
      channel_direction_metric(:) = 0.0_dp
      channel_raw_metric(:) = 0.0_dp
      channel_search_metric(:) = 0.0_dp

      gnorm_cross = 0.0_dp
      IF (do_ks) THEN
         DO ispin = 1, nspin
            norm_weight = ot_kpoint_channel_norm_weight(qs_ot_env(ispin))
            CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin)%matrix_gx_old, tmp)
            gnorm_cross = gnorm_cross + tmp
            CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin)%matrix_gx, tmp)
            raw_matrix_gnorm_re = raw_matrix_gnorm_re + tmp
            raw_matrix_gnorm_metric_re = raw_matrix_gnorm_metric_re + norm_weight*tmp
            channel_raw_metric(ispin) = channel_raw_metric(ispin) + norm_weight*tmp
            IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
               CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx_im, qs_ot_env(ispin)%matrix_gx_old_im, tmp)
               gnorm_cross = gnorm_cross + tmp
               CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx_im, qs_ot_env(ispin)%matrix_gx_im, tmp)
               raw_matrix_gnorm_im = raw_matrix_gnorm_im + tmp
               raw_matrix_gnorm_metric_im = raw_matrix_gnorm_metric_im + norm_weight*tmp
               channel_raw_metric(ispin) = channel_raw_metric(ispin) + norm_weight*tmp
            END IF
         END DO
         IF (qs_ot_env(1)%settings%do_rotation) THEN
            DO ispin = 1, nspin
               CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx, qs_ot_env(ispin)%rot_mat_gx_old, tmp)
               ! added 0.5, because we have (antisymmetry) only half the number of variables
               gnorm_cross = gnorm_cross + 0.5_dp*tmp
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx_im, qs_ot_env(ispin)%rot_mat_gx_old_im, tmp)
                  gnorm_cross = gnorm_cross + 0.5_dp*tmp
               END IF
            END DO
         END IF
      END IF
      IF (do_ener) THEN
         DO ispin = 1, nspin
            tmp = DOT_PRODUCT(qs_ot_env(ispin)%ener_gx, qs_ot_env(ispin)%ener_gx_old)
            gnorm_cross = gnorm_cross + tmp
         END DO
      END IF
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, gnorm_cross)
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, raw_matrix_gnorm_re)
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, raw_matrix_gnorm_im)
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, raw_matrix_gnorm_metric_re)
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, raw_matrix_gnorm_metric_im)

      IF (ASSOCIATED(qs_ot_env(1)%preconditioner)) THEN

         DO ispin = 1, nspin
            norm_weight = ot_kpoint_channel_norm_weight(qs_ot_env(ispin))
            IF (qs_ot_env(ispin)%has_complex_kpoint_state .AND. &
                ASSOCIATED(qs_ot_env(ispin)%preconditioner%complex_fm)) THEN
               CALL apply_preconditioner(qs_ot_env(ispin)%preconditioner, &
                                         qs_ot_env(ispin)%matrix_gx, &
                                         qs_ot_env(ispin)%matrix_gx_im, &
                                         qs_ot_env(ispin)%matrix_gx_old, &
                                         qs_ot_env(ispin)%matrix_gx_old_im)
            ELSE
               CALL apply_preconditioner(qs_ot_env(ispin)%preconditioner, &
                                         qs_ot_env(ispin)%matrix_gx, &
                                         qs_ot_env(ispin)%matrix_gx_old)
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  CALL apply_preconditioner(qs_ot_env(ispin)%preconditioner, &
                                            qs_ot_env(ispin)%matrix_gx_im, &
                                            qs_ot_env(ispin)%matrix_gx_old_im)
               END IF
            END IF
            IF (qs_ot_env(ispin)%has_complex_kpoint_state .AND. &
                ASSOCIATED(qs_ot_env(ispin)%matrix_occ_precond_scale)) THEN
               CALL dbcsr_scale_by_vector(qs_ot_env(ispin)%matrix_gx_old, &
                                          alpha=qs_ot_env(ispin)%matrix_occ_precond_scale, &
                                          side='right')
            END IF
            IF (qs_ot_env(ispin)%has_kpoint_context) THEN
               CALL dbcsr_scale(qs_ot_env(ispin)%matrix_gx_old, norm_weight)
            END IF
            IF (qs_ot_env(ispin)%has_complex_kpoint_state .AND. &
                ASSOCIATED(qs_ot_env(ispin)%matrix_occ_precond_scale)) THEN
               CALL dbcsr_scale_by_vector(qs_ot_env(ispin)%matrix_gx_old_im, &
                                          alpha=qs_ot_env(ispin)%matrix_occ_precond_scale, &
                                          side='right')
            END IF
            IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
               CALL dbcsr_scale(qs_ot_env(ispin)%matrix_gx_old_im, norm_weight)
            END IF
         END DO
         qs_ot_env(1)%gnorm = 0.0_dp
         IF (do_ks) THEN
            DO ispin = 1, nspin
               norm_weight = ot_kpoint_channel_norm_weight(qs_ot_env(ispin))
               CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin)%matrix_gx_old, tmp)
               channel_search_metric(ispin) = channel_search_metric(ispin) + tmp
               qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + tmp
               gnorm_matrix_re = gnorm_matrix_re + tmp
               IF (kpoint_preconditioned_direction .OR. &
                   qs_ot_env(ispin)%has_kpoint_context) THEN
                  gnorm_metric = gnorm_metric + tmp
               ELSE
                  gnorm_metric = gnorm_metric + norm_weight*tmp
               END IF
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx_im, qs_ot_env(ispin)%matrix_gx_old_im, tmp)
                  channel_search_metric(ispin) = channel_search_metric(ispin) + tmp
                  qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + tmp
                  gnorm_matrix_im = gnorm_matrix_im + tmp
                  IF (kpoint_preconditioned_direction .OR. &
                      qs_ot_env(ispin)%has_kpoint_context) THEN
                     gnorm_metric = gnorm_metric + tmp
                  ELSE
                     gnorm_metric = gnorm_metric + norm_weight*tmp
                  END IF
               END IF
               IF (debug_diagnostic) THEN
                  CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx_old, &
                                 qs_ot_env(ispin)%matrix_gx_old, tmp)
                  channel_direction_metric(ispin) = channel_direction_metric(ispin) + tmp
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx_old_im, &
                                    qs_ot_env(ispin)%matrix_gx_old_im, tmp)
                     channel_direction_metric(ispin) = channel_direction_metric(ispin) + tmp
                  END IF
               END IF
            END DO
            IF (qs_ot_env(1)%gnorm < 0.0_dp) THEN
               WRITE (cp_logger_get_default_unit_nr(logger), *) "WARNING Preconditioner not positive definite !"
            END IF
            IF (legacy_preconditioned_gradient) THEN
               ! Preserve the established Gamma-point CG metric after computing g.P^-1.g.
               ! Coupled metallic OT keeps the raw gradient for its physical Mermin slope.
               DO ispin = 1, nspin
                  CALL dbcsr_copy(qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin)%matrix_gx_old)
               END DO
            END IF
            IF (qs_ot_env(1)%settings%do_rotation) THEN
               DO ispin = 1, nspin
                  norm_weight = ot_kpoint_channel_norm_weight(qs_ot_env(ispin))
                  IF (kpoint_preconditioned_direction) THEN
                     CPASSERT(ASSOCIATED(qs_ot_env(ispin)%rot_mat_pgx))
                     CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx, qs_ot_env(ispin)%rot_mat_pgx, tmp)
                     CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_gx_old, qs_ot_env(ispin)%rot_mat_pgx)
                  ELSE
                     ! right now no preconditioner yet
                     CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_gx_old, qs_ot_env(ispin)%rot_mat_gx)
                     IF (qs_ot_env(ispin)%has_kpoint_context) THEN
                        CALL dbcsr_scale(qs_ot_env(ispin)%rot_mat_gx_old, norm_weight)
                     END IF
                     CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx, qs_ot_env(ispin)%rot_mat_gx_old, tmp)
                  END IF
                  ! added 0.5, because we have (antisymmetry) only half the number of variables
                  qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + 0.5_dp*tmp
                  gnorm_rot_re = gnorm_rot_re + 0.5_dp*tmp
                  IF (kpoint_preconditioned_direction .OR. &
                      qs_ot_env(ispin)%has_kpoint_context) THEN
                     gnorm_metric = gnorm_metric + 0.5_dp*tmp
                  ELSE
                     gnorm_metric = gnorm_metric + norm_weight*0.5_dp*tmp
                  END IF
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     IF (kpoint_preconditioned_direction) THEN
                        CPASSERT(ASSOCIATED(qs_ot_env(ispin)%rot_mat_pgx_im))
                        CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx_im, qs_ot_env(ispin)%rot_mat_pgx_im, tmp)
                        CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_gx_old_im, qs_ot_env(ispin)%rot_mat_pgx_im)
                     ELSE
                        CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_gx_old_im, qs_ot_env(ispin)%rot_mat_gx_im)
                        CALL dbcsr_scale(qs_ot_env(ispin)%rot_mat_gx_old_im, norm_weight)
                        CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx_im, qs_ot_env(ispin)%rot_mat_gx_old_im, tmp)
                     END IF
                     qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + 0.5_dp*tmp
                     gnorm_rot_im = gnorm_rot_im + 0.5_dp*tmp
                     IF (kpoint_preconditioned_direction .OR. &
                         qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                        gnorm_metric = gnorm_metric + 0.5_dp*tmp
                     ELSE
                        gnorm_metric = gnorm_metric + norm_weight*0.5_dp*tmp
                     END IF
                  END IF
               END DO
            END IF
         END IF
         IF (do_ener) THEN
            DO ispin = 1, nspin
               norm_weight = ot_kpoint_channel_norm_weight(qs_ot_env(ispin))
               IF (kpoint_preconditioned_direction) THEN
                  CPASSERT(ASSOCIATED(qs_ot_env(ispin)%ener_pgx))
                  tmp = DOT_PRODUCT(qs_ot_env(ispin)%ener_gx, qs_ot_env(ispin)%ener_pgx)
                  qs_ot_env(ispin)%ener_gx_old = qs_ot_env(ispin)%ener_pgx
               ELSE
                  qs_ot_env(ispin)%ener_gx_old = norm_weight*qs_ot_env(ispin)%ener_gx
                  tmp = DOT_PRODUCT(qs_ot_env(ispin)%ener_gx, qs_ot_env(ispin)%ener_gx_old)
               END IF
               qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + tmp
               gnorm_ener = gnorm_ener + tmp
               IF (kpoint_preconditioned_direction .OR. &
                   qs_ot_env(ispin)%has_kpoint_context) THEN
                  gnorm_metric = gnorm_metric + tmp
               ELSE
                  gnorm_metric = gnorm_metric + norm_weight*tmp
               END IF
            END DO
         END IF
      ELSE
         IF (do_ks) THEN
            qs_ot_env(1)%gnorm = 0.0_dp
            DO ispin = 1, nspin
               norm_weight = ot_kpoint_channel_norm_weight(qs_ot_env(ispin))
               CALL dbcsr_copy(qs_ot_env(ispin)%matrix_gx_old, qs_ot_env(ispin)%matrix_gx)
               IF (qs_ot_env(ispin)%has_kpoint_context) THEN
                  CALL dbcsr_scale(qs_ot_env(ispin)%matrix_gx_old, norm_weight)
               END IF
               CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin)%matrix_gx_old, tmp)
               channel_search_metric(ispin) = channel_search_metric(ispin) + tmp
               qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + tmp
               gnorm_matrix_re = gnorm_matrix_re + tmp
               gnorm_metric = gnorm_metric + tmp
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  CALL dbcsr_copy(qs_ot_env(ispin)%matrix_gx_old_im, qs_ot_env(ispin)%matrix_gx_im)
                  CALL dbcsr_scale(qs_ot_env(ispin)%matrix_gx_old_im, norm_weight)
                  CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx_im, qs_ot_env(ispin)%matrix_gx_old_im, tmp)
                  channel_search_metric(ispin) = channel_search_metric(ispin) + tmp
                  qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + tmp
                  gnorm_matrix_im = gnorm_matrix_im + tmp
                  gnorm_metric = gnorm_metric + tmp
               END IF
               IF (debug_diagnostic) THEN
                  CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx_old, &
                                 qs_ot_env(ispin)%matrix_gx_old, tmp)
                  channel_direction_metric(ispin) = channel_direction_metric(ispin) + tmp
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx_old_im, &
                                    qs_ot_env(ispin)%matrix_gx_old_im, tmp)
                     channel_direction_metric(ispin) = channel_direction_metric(ispin) + tmp
                  END IF
               END IF
            END DO
            IF (qs_ot_env(1)%settings%do_rotation) THEN
               DO ispin = 1, nspin
                  norm_weight = ot_kpoint_channel_norm_weight(qs_ot_env(ispin))
                  CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_gx_old, qs_ot_env(ispin)%rot_mat_gx)
                  IF (qs_ot_env(ispin)%has_kpoint_context) THEN
                     CALL dbcsr_scale(qs_ot_env(ispin)%rot_mat_gx_old, norm_weight)
                  END IF
                  CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx, qs_ot_env(ispin)%rot_mat_gx_old, tmp)
                  ! added 0.5, because we have (antisymmetry) only half the number of variables
                  qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + 0.5_dp*tmp
                  gnorm_rot_re = gnorm_rot_re + 0.5_dp*tmp
                  gnorm_metric = gnorm_metric + 0.5_dp*tmp
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_gx_old_im, qs_ot_env(ispin)%rot_mat_gx_im)
                     CALL dbcsr_scale(qs_ot_env(ispin)%rot_mat_gx_old_im, norm_weight)
                     CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx_im, qs_ot_env(ispin)%rot_mat_gx_old_im, tmp)
                     qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + 0.5_dp*tmp
                     gnorm_rot_im = gnorm_rot_im + 0.5_dp*tmp
                     gnorm_metric = gnorm_metric + 0.5_dp*tmp
                  END IF
               END DO
            END IF
         END IF
         IF (do_ener) THEN
            DO ispin = 1, nspin
               norm_weight = ot_kpoint_channel_norm_weight(qs_ot_env(ispin))
               qs_ot_env(ispin)%ener_gx_old = norm_weight*qs_ot_env(ispin)%ener_gx
               tmp = DOT_PRODUCT(qs_ot_env(ispin)%ener_gx, qs_ot_env(ispin)%ener_gx_old)
               qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + tmp
               gnorm_ener = gnorm_ener + tmp
               gnorm_metric = gnorm_metric + tmp
            END DO
         END IF
      END IF

      IF (debug_diagnostic .AND. do_ks) THEN
         DO ispin = 1, nspin
            CALL dbcsr_get_info(qs_ot_env(ispin)%matrix_gx, &
                                nfullrows_total=matrix_rows, nfullcols_total=matrix_cols)
            matrix_nvar = REAL(INT(matrix_rows, KIND=int_8)*INT(matrix_cols, KIND=int_8), KIND=dp)
            IF (qs_ot_env(ispin)%has_complex_kpoint_state) matrix_nvar = 2.0_dp*matrix_nvar
            matrix_nvar = REAL(MAX(1, qs_ot_env(ispin)%kpoint_multiplicity), KIND=dp)*matrix_nvar
            raw_channel_rms = SQRT(MAX(0.0_dp, channel_raw_metric(ispin))/MAX(matrix_nvar, 1.0_dp))
            search_channel_rms = &
               SQRT(ABS(channel_search_metric(ispin))/MAX(matrix_nvar, 1.0_dp))
            search_direction_rms = &
               SQRT(MAX(0.0_dp, channel_direction_metric(ispin))/MAX(matrix_nvar, 1.0_dp))
            occ_scale_min = 1.0_dp
            occ_scale_avg = 1.0_dp
            occ_scale_max = 1.0_dp
            IF (ASSOCIATED(qs_ot_env(ispin)%matrix_occ_precond_scale)) THEN
               occ_scale_min = MINVAL(qs_ot_env(ispin)%matrix_occ_precond_scale)
               occ_scale_avg = SUM(qs_ot_env(ispin)%matrix_occ_precond_scale)/ &
                               REAL(MAX(1, SIZE(qs_ot_env(ispin)%matrix_occ_precond_scale)), KIND=dp)
               occ_scale_max = MAXVAL(qs_ot_env(ispin)%matrix_occ_precond_scale)
            END IF
            IF (debug_diagnostic_source) THEN
               WRITE (cp_logger_get_default_unit_nr(logger), '(A,3(1X,I4),1X,L1,1X,I4,6(1X,ES16.8))') &
                  " OT CG matrix channel spin k mult occ-active samples raw search direction omin oavg omax:", &
                  qs_ot_env(ispin)%spin_index, qs_ot_env(ispin)%kpoint_index, &
                  qs_ot_env(ispin)%kpoint_multiplicity, &
                  qs_ot_occupation_preconditioner_active(qs_ot_env(ispin)), &
                  qs_ot_env(ispin)%kpoint_occ_preconditioner_ready_samples, &
                  raw_channel_rms, search_channel_rms, search_direction_rms, &
                  occ_scale_min, occ_scale_avg, occ_scale_max
            END IF
         END DO
      END IF

      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, qs_ot_env(1)%gnorm)
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, gnorm_metric)
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, gnorm_matrix_re)
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, gnorm_matrix_im)
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, gnorm_rot_re)
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, gnorm_rot_im)
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, gnorm_ener)
      physical_control_delta = 0.0_dp
      nvar = ot_kpoint_channel_nvar(qs_ot_env, do_ks, do_ener)
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, nvar)
      ! Handling the case of no free variables to optimize
      IF (nvar > 0.0_dp) THEN
         raw_delta = SQRT(ABS(gnorm_metric)/nvar)
         qs_ot_env(1)%delta = raw_delta
         qs_ot_env(1)%kpoint_raw_delta = raw_delta
         qs_ot_env(1)%kpoint_physical_delta = -1.0_dp
         qs_ot_env(1)%kpoint_projected_delta = -1.0_dp
         IF (complex_kpoint_history .AND. qs_ot_env(1)%settings%occupation_preconditioner .AND. &
             qs_ot_env(1)%settings%do_rotation .AND. do_ener .AND. &
             qs_ot_env(1)%kpoint_physical_delta_groups > 0) THEN
            matrix_gnorm_metric = MAX(0.0_dp, raw_matrix_gnorm_metric_re) + &
                                  MAX(0.0_dp, raw_matrix_gnorm_metric_im)
            physical_gnorm_metric = matrix_gnorm_metric + MAX(0.0_dp, qs_ot_env(1)%kpoint_physical_gnorm)
            physical_delta = SQRT(ABS(physical_gnorm_metric)/nvar)
            qs_ot_env(1)%kpoint_physical_delta = physical_delta
            qs_ot_env(1)%kpoint_projected_delta = physical_delta
            CALL ot_kpoint_physical_residual_split( &
               raw_delta, physical_delta, qs_ot_env(1)%kpoint_physical_delta_ready, &
               qs_ot_env(1)%kpoint_eps_scf, qs_ot_env(1)%mermin_ls_accepted_quality_samples, &
               qs_ot_env(1)%mermin_ls_accepted_quality_ema, qs_ot_env(1)%mermin_ls_accepted_gap_ema, &
               qs_ot_env(1)%delta, physical_control_delta)
         END IF
         beta_pr = (qs_ot_env(1)%gnorm - gnorm_cross)/qs_ot_env(1)%gnorm_old
      ELSE
         qs_ot_env(1)%delta = 0.0_dp
         qs_ot_env(1)%kpoint_raw_delta = 0.0_dp
         qs_ot_env(1)%kpoint_physical_delta = 0.0_dp
         qs_ot_env(1)%kpoint_projected_delta = 0.0_dp
         beta_pr = 0.0_dp
      END IF
      IF (complex_kpoint_history .AND. qs_ot_env(1)%mermin_ls_sample_pending) THEN
         sample_accepted_drop = qs_ot_env(1)%mermin_ls_reference_energy - qs_ot_env(1)%etotal
         sample_gap = qs_ot_env(1)%etotal - qs_ot_env(1)%mermin_ls_best_sample_energy
         sample_quality = 0.0_dp
         IF (qs_ot_env(1)%mermin_ls_best_sample_drop > EPSILON(1.0_dp)) THEN
            sample_quality = sample_accepted_drop/qs_ot_env(1)%mermin_ls_best_sample_drop
         END IF
         sample_guard_drop = sample_accepted_drop
         IF (qs_ot_env(1)%mermin_ls_model_quality < 0.5_dp) THEN
            sample_guard_drop = MAX(sample_guard_drop, qs_ot_env(1)%mermin_ls_model_drop)
         END IF
         physical_raw_gap_guard = complex_kpoint_history .AND. &
                                  qs_ot_occupation_preconditioner_active(qs_ot_env(1)) .AND. &
                                  qs_ot_env(1)%settings%do_rotation .AND. do_ener .AND. &
                                  qs_ot_env(1)%kpoint_physical_delta_ready .AND. &
                                  qs_ot_env(1)%kpoint_physical_delta_groups > 0 .AND. &
                                  qs_ot_env(1)%kpoint_physical_delta > 0.0_dp .AND. &
                                  qs_ot_env(1)%kpoint_raw_delta > qs_ot_env(1)%delta .AND. &
                                  sample_guard_drop > ot_cg_kp_physical_raw_delta_mermin_floor
         IF (physical_raw_gap_guard) THEN
            raw_guard_delta = ot_cg_kp_physical_raw_delta_floor_ratio*qs_ot_env(1)%kpoint_raw_delta
            ! Release the raw preconditioned floor only in the calibrated tail,
            ! where the fixed-N projected Mermin residual already satisfies
            ! EPS_SCF and accepted line-search samples remain predictive.
            projected_raw_tail_release = qs_ot_env(1)%kpoint_eps_scf > 0.0_dp .AND. &
                                         qs_ot_env(1)%kpoint_projected_delta > 0.0_dp .AND. &
                                         qs_ot_env(1)%kpoint_projected_delta <= &
                                         qs_ot_env(1)%kpoint_eps_scf .AND. &
                                         qs_ot_env(1)%mermin_ls_accepted_quality_samples >= &
                                         ot_cg_kp_physical_control_min_samples .AND. &
                                         qs_ot_env(1)%mermin_ls_accepted_quality_ema >= &
                                         ot_cg_kp_physical_raw_delta_tail_quality .AND. &
                                         qs_ot_env(1)%mermin_ls_accepted_gap_ema <= &
                                         ot_cg_kp_physical_raw_delta_tail_gap
            IF (projected_raw_tail_release) THEN
               raw_guard_delta = MIN(raw_guard_delta, qs_ot_env(1)%kpoint_projected_delta)
            END IF
            IF (raw_guard_delta > qs_ot_env(1)%delta) THEN
               qs_ot_env(1)%delta = raw_guard_delta
               qs_ot_env(1)%kpoint_physical_delta = MAX(qs_ot_env(1)%kpoint_physical_delta, raw_guard_delta)
               IF (logger%iter_info%print_level >= debug_print_level .AND. logger%para_env%is_source()) THEN
                  WRITE (cp_logger_get_default_unit_nr(logger), '(A,8(1X,ES16.8))') &
                     " OT CG physical residual raw-gap guard raw floor projected adrop gdrop mdrop mq delta:", &
                     qs_ot_env(1)%kpoint_raw_delta, raw_guard_delta, &
                     qs_ot_env(1)%kpoint_projected_delta, sample_accepted_drop, &
                     sample_guard_drop, qs_ot_env(1)%mermin_ls_model_drop, &
                     qs_ot_env(1)%mermin_ls_model_quality, qs_ot_env(1)%delta
               END IF
            END IF
         END IF
         physical_tail_delta = kpoint_projected_tail_delta(qs_ot_env)
         residual_growth_delta = ot_kpoint_mermin_current_guard_delta( &
                                 qs_ot_env(1)%delta, physical_tail_delta, &
                                 qs_ot_env(1)%mermin_ls_reference_delta)
         previous_delta = ot_kpoint_mermin_reference_guard_delta( &
                          qs_ot_env(1)%schur_model_delta, &
                          qs_ot_env(1)%mermin_ls_reference_delta)
         IF (qs_ot_occupation_preconditioner_active(qs_ot_env(1)) .AND. &
             qs_ot_env(1)%settings%do_rotation .AND. do_ener) THEN
            CALL ot_kpoint_cg_residual_growth_guard(previous_delta, &
                                                    residual_growth_delta, sample_accepted_drop, &
                                                    qs_ot_env(1)%kpoint_eps_scf, &
                                                    mermin_residual_growth_reset, &
                                                    mermin_residual_growth_quality, &
                                                    mermin_residual_catastrophic_reset, &
                                                    strict_growth= &
                                                    ot_kpoint_mermin_shared_spin_group( &
                                                    qs_ot_env%spin_index, &
                                                    qs_ot_env%occupation_group_index) .AND. &
                                                    qs_ot_env(1)%mermin_ls_reference_delta > 0.0_dp)
         END IF
         ! Residual guards depend on distributed floating-point reductions.  A
         ! rank-local threshold crossing must not split the CG control flow or
         ! leave only part of a DBCSR subgroup with updated minimizer history.
         growth_reset_any = MERGE(1, 0, mermin_residual_growth_reset)
         catastrophic_reset_any = MERGE(1, 0, mermin_residual_catastrophic_reset)
         CALL ot_mini_max_kpoint_integer(qs_ot_env, para_env_inter_kp, growth_reset_any)
         CALL ot_mini_max_kpoint_integer(qs_ot_env, para_env_inter_kp, catastrophic_reset_any)
         CALL ot_mini_min_kpoint_scalar(qs_ot_env, para_env_inter_kp, &
                                        mermin_residual_growth_quality)
         mermin_residual_growth_reset = growth_reset_any == 1
         mermin_residual_catastrophic_reset = catastrophic_reset_any == 1
         sample_residual_growth_ratio = 0.0_dp
         IF (previous_delta > EPSILON(1.0_dp)) THEN
            sample_residual_growth_ratio = residual_growth_delta/previous_delta
         END IF
         IF (mermin_residual_growth_reset) THEN
            qs_ot_env(1)%schur_residual_history_reset_pending = .TRUE.
            qs_ot_env(1)%schur_residual_history_quality = &
               MIN(qs_ot_env(1)%schur_residual_history_quality, mermin_residual_growth_quality)
            IF (logger%iter_info%print_level >= debug_print_level .AND. logger%para_env%is_source()) THEN
               WRITE (cp_logger_get_default_unit_nr(logger), '(A,1X,L1,6(1X,ES16.8))') &
                  " OT Mermin LS residual growth reset catastrophic prev next ratio quality adrop bdrop:", &
                  mermin_residual_catastrophic_reset, &
                  previous_delta, residual_growth_delta, sample_residual_growth_ratio, &
                  mermin_residual_growth_quality, sample_accepted_drop, &
                  qs_ot_env(1)%mermin_ls_best_sample_drop
            END IF
         END IF
         IF (logger%iter_info%print_level >= debug_print_level .AND. logger%para_env%is_source()) THEN
            WRITE (cp_logger_get_default_unit_nr(logger), '(A,1X,I5,12(1X,ES16.8))') &
               " OT Mermin LS accepted-state audit idx ref best accepted bdrop adrop gap q mdrop mq bpos next delta:", &
               qs_ot_env(1)%mermin_ls_best_sample_index, qs_ot_env(1)%mermin_ls_reference_energy, &
               qs_ot_env(1)%mermin_ls_best_sample_energy, qs_ot_env(1)%etotal, &
               qs_ot_env(1)%mermin_ls_best_sample_drop, sample_accepted_drop, sample_gap, &
               sample_quality, qs_ot_env(1)%mermin_ls_model_drop, qs_ot_env(1)%mermin_ls_model_quality, &
               qs_ot_env(1)%mermin_ls_best_sample_pos, qs_ot_env(1)%mermin_ls_next_pos, qs_ot_env(1)%delta
         END IF
         IF (PRESENT(residual_closure_rejected)) THEN
            residual_closure_rejected = mermin_residual_catastrophic_reset .AND. &
                                        qs_ot_env(1)%line_search_count >= 3
         END IF
         IF (preflight_only) THEN
            CALL restore_ot_cg_history(qs_ot_env, history_snapshot)
            CALL timestop(handle)
            RETURN
         ELSE IF (mermin_residual_catastrophic_reset .AND. &
                  qs_ot_env(1)%line_search_might_be_done .AND. &
                  qs_ot_env(1)%line_search_count >= 3 .AND. PRESENT(residual_closure_rejected)) THEN
            CALL timestop(handle)
            RETURN
         END IF
         qs_ot_env(1)%mermin_ls_sample_pending = .FALSE.
      END IF
      IF (preflight_only) THEN
         CALL restore_ot_cg_history(qs_ot_env, history_snapshot)
         CALL timestop(handle)
         RETURN
      END IF
      IF (logger%iter_info%print_level >= debug_print_level .AND. logger%para_env%is_source()) THEN
         WRITE (cp_logger_get_default_unit_nr(logger), '(A,6(1X,ES16.8))') &
            " OT CG channel norms C_re C_im R_re R_im eps total:", &
            gnorm_matrix_re, gnorm_matrix_im, gnorm_rot_re, gnorm_rot_im, &
            gnorm_ener, qs_ot_env(1)%gnorm
         IF (qs_ot_env(1)%kpoint_physical_delta_groups > 0) THEN
            WRITE (cp_logger_get_default_unit_nr(logger), '(A,1X,L1,1X,I6,9(1X,ES16.8))') &
               " OT CG physical residual ready groups raw model exit physical projected matrix phys_rot phys_eps gauge:", &
               qs_ot_env(1)%kpoint_physical_delta_ready, qs_ot_env(1)%kpoint_physical_delta_groups, &
               qs_ot_env(1)%kpoint_raw_delta, qs_ot_env(1)%delta, physical_control_delta, &
               qs_ot_env(1)%kpoint_physical_delta, &
               qs_ot_env(1)%kpoint_projected_delta, &
               MAX(0.0_dp, raw_matrix_gnorm_re) + MAX(0.0_dp, raw_matrix_gnorm_im), &
               qs_ot_env(1)%kpoint_physical_rot_gnorm, &
               qs_ot_env(1)%kpoint_physical_energy_gnorm, qs_ot_env(1)%kpoint_gauge_gnorm
         END IF
         IF (qs_ot_env(1)%kpoint_physical_delta_groups > 0) THEN
            WRITE (cp_logger_get_default_unit_nr(logger), '(A,4(1X,ES16.8))') &
               " OT CG matrix residual raw_re raw_im precond_re precond_im:", &
               raw_matrix_gnorm_re, raw_matrix_gnorm_im, gnorm_matrix_re, gnorm_matrix_im
         END IF
      END IF
      beta_raw = beta_pr
      negative_beta = beta_pr < 0.0_dp
      beta_pr = MAX(beta_pr, 0.0_dp) ! reset to SD
      gnorm_ratio = qs_ot_env(1)%gnorm/MAX(qs_ot_env(1)%gnorm_old, EPSILON(1.0_dp))
      physical_tail_delta = kpoint_projected_tail_delta(qs_ot_env)
      ! The marked Schur flavor transitions change the metric used to map the
      ! physical gradient into a search direction.  The previous conjugate
      ! direction is therefore not a history vector in the current metric.
      ! Keeping beta across this transition can remain downhill locally while
      ! amplifying the physical rotation residual in the following SCF step.
      ! Retain the productive nonlinear-CG history outside the residual tail.
      IF (ot_kpoint_cg_metric_history_reset( &
          schur_history_change, beta_pr, physical_tail_delta, residual_tail_min, &
          schur_calibration_change)) THEN
         beta_pr = 0.0_dp
         reset_to_sd = .TRUE.
      END IF
      ! Projector-like Schur flavors can validly alternate in the metallic tail.
      ! Do not carry a huge CG beta across that local response-model change.
      schur_projector_jump_reset = complex_kpoint_history .AND. &
                                   qs_ot_occupation_preconditioner_active(qs_ot_env(1)) .AND. &
                                   qs_ot_env(1)%settings%do_rotation .AND. do_ener .AND. &
                                   beta_pr > ot_cg_kp_history_jump_beta .AND. &
                                   qs_ot_env(1)%gnorm > ot_cg_kp_projector_jump_new_norm .AND. &
                                   gnorm_ratio > ot_cg_kp_projector_jump_ratio .AND. &
                                   schur_attempt >= 3 .AND. schur_last_attempt >= 3 .AND. &
                                   schur_flavor >= 35 .AND. schur_last_flavor >= 35
      schur_projector_tail_reset = complex_kpoint_history .AND. &
                                   qs_ot_occupation_preconditioner_active(qs_ot_env(1)) .AND. &
                                   qs_ot_env(1)%settings%do_rotation .AND. do_ener .AND. &
                                   ot_kpoint_cg_tail_metric_tracked( &
                                   schur_attempt, schur_last_attempt, &
                                   schur_flavor, schur_last_flavor, &
                                   physical_tail_delta, residual_tail_min) .AND. &
                                   beta_pr > ot_cg_kp_projector_tail_beta .AND. &
                                   qs_ot_env(1)%gnorm > ot_cg_kp_projector_tail_norm .AND. &
                                   gnorm_ratio > ot_cg_kp_projector_tail_ratio
      residual_tail_spike_reset = complex_kpoint_history .AND. &
                                  qs_ot_occupation_preconditioner_active(qs_ot_env(1)) .AND. &
                                  qs_ot_env(1)%settings%do_rotation .AND. do_ener .AND. &
                                  qs_ot_env(1)%kpoint_physical_delta_ready .AND. &
                                  physical_tail_delta > residual_tail_min .AND. &
                                  physical_tail_delta < ot_cg_kp_residual_tail_max_delta .AND. &
                                  beta_pr > ot_cg_kp_residual_tail_jump_beta .AND. &
                                  qs_ot_env(1)%gnorm > ot_cg_kp_residual_tail_jump_norm .AND. &
                                  gnorm_ratio > ot_cg_kp_residual_tail_jump_ratio
      history_spike_reset = complex_kpoint_history .AND. &
                            qs_ot_occupation_preconditioner_active(qs_ot_env(1)) .AND. &
                            qs_ot_env(1)%settings%do_rotation .AND. do_ener .AND. &
                            beta_pr > 0.0_dp .AND. &
                            (schur_projector_jump_reset .OR. schur_projector_tail_reset .OR. &
                             residual_tail_spike_reset .OR. &
                             (qs_ot_env(1)%gnorm > ot_cg_kp_history_spike_new_norm .AND. &
                              ((qs_ot_env(1)%gnorm_old < ot_cg_kp_history_spike_old_norm .AND. &
                                gnorm_ratio > ot_cg_kp_history_spike_ratio) .OR. &
                               (beta_pr > ot_cg_kp_history_jump_beta .AND. &
                                qs_ot_env(1)%gnorm > ot_cg_kp_history_jump_new_norm .AND. &
                                gnorm_ratio > ot_cg_kp_history_jump_ratio))))
      IF (history_spike_reset) THEN
         beta_pr = 0.0_dp
         reset_to_sd = .TRUE.
      END IF
      residual_spike_history_reset = history_spike_reset .AND. &
                                     qs_ot_env(1)%kpoint_physical_delta_ready .AND. &
                                     physical_tail_delta > &
                                     mermin_tail_min .AND. &
                                     physical_tail_delta < &
                                     ot_cg_kp_mermin_recovery_tail_max_delta .AND. &
                                     (schur_projector_tail_reset .OR. residual_tail_spike_reset)
      IF (residual_spike_history_reset) THEN
         residual_spike_quality = ot_cg_kp_residual_spike_quality_ceiling
         IF (gnorm_ratio > EPSILON(1.0_dp)) THEN
            residual_spike_quality = MIN(residual_spike_quality, 1.0_dp/SQRT(gnorm_ratio))
         END IF
         qs_ot_env(1)%schur_residual_history_reset_pending = .TRUE.
         qs_ot_env(1)%schur_residual_history_quality = &
            MIN(qs_ot_env(1)%schur_residual_history_quality, residual_spike_quality)
         IF (logger%iter_info%print_level >= debug_print_level .AND. logger%para_env%is_source()) THEN
            WRITE (cp_logger_get_default_unit_nr(logger), '(A,4(1X,ES16.8),2(1X,I6))') &
               " OT CG residual spike trust reset q ratio delta rscale flavor last:", &
               residual_spike_quality, gnorm_ratio, physical_tail_delta, &
               MAX(ot_cg_kp_residual_trust_min_scale, &
                   MIN(1.0_dp, SQRT(MAX(residual_spike_quality, &
                                        ot_cg_kp_residual_trust_quality_floor)))), &
               schur_flavor, schur_last_flavor
         END IF
      END IF
      residual_history_reset = complex_kpoint_history .AND. &
                               qs_ot_occupation_preconditioner_active(qs_ot_env(1)) .AND. &
                               qs_ot_env(1)%settings%do_rotation .AND. do_ener .AND. &
                               qs_ot_env(1)%schur_residual_history_reset_pending
      residual_trust_scale = 1.0_dp
      IF (residual_history_reset) THEN
         residual_trust_quality = MAX(0.0_dp, MIN(1.0_dp, qs_ot_env(1)%schur_residual_history_quality))
         residual_trust_scale = SQRT(MAX(residual_trust_quality, &
                                         ot_cg_kp_residual_trust_quality_floor))
         residual_trust_scale = MAX(ot_cg_kp_residual_trust_min_scale, MIN(1.0_dp, residual_trust_scale))
         IF (qs_ot_env(1)%ds_min > 0.0_dp) qs_ot_env(1)%ds_min = qs_ot_env(1)%ds_min*residual_trust_scale
      END IF
      IF (residual_history_reset .AND. beta_pr > 0.0_dp) THEN
         IF (residual_trust_quality >= ot_cg_kp_residual_history_damp_quality) THEN
            beta_pr = beta_pr*residual_trust_scale
         ELSE
            beta_pr = 0.0_dp
            reset_to_sd = .TRUE.
         END IF
      END IF
      qs_ot_env(1)%schur_residual_history_reset_pending = .FALSE.

      test_down = 0.0_dp
      IF (do_ks) THEN
         DO ispin = 1, nspin
            CALL dbcsr_add(qs_ot_env(ispin)%matrix_dx, qs_ot_env(ispin)%matrix_gx_old, &
                           alpha_scalar=beta_pr, beta_scalar=-1.0_dp)
            CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin)%matrix_dx, tmp)
            test_down = test_down + tmp
            IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
               CALL dbcsr_add(qs_ot_env(ispin)%matrix_dx_im, qs_ot_env(ispin)%matrix_gx_old_im, &
                              alpha_scalar=beta_pr, beta_scalar=-1.0_dp)
               CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx_im, qs_ot_env(ispin)%matrix_dx_im, tmp)
               test_down = test_down + tmp
            END IF
            IF (qs_ot_env(1)%settings%do_rotation) THEN
               CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_dx, qs_ot_env(ispin)%rot_mat_gx_old, &
                              alpha_scalar=beta_pr, beta_scalar=-1.0_dp)
               CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx, qs_ot_env(ispin)%rot_mat_dx, tmp)
               test_down = test_down + 0.5_dp*tmp
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_dx_im, qs_ot_env(ispin)%rot_mat_gx_old_im, &
                                 alpha_scalar=beta_pr, beta_scalar=-1.0_dp)
                  CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx_im, qs_ot_env(ispin)%rot_mat_dx_im, tmp)
                  test_down = test_down + 0.5_dp*tmp
               END IF
            END IF
         END DO
      END IF
      IF (do_ener) THEN
         DO ispin = 1, nspin
            qs_ot_env(ispin)%ener_dx = beta_pr*qs_ot_env(ispin)%ener_dx - qs_ot_env(ispin)%ener_gx_old
         END DO
         IF (complex_kpoint_history) THEN
            CALL prepare_linear_kpoint_energy_direction(qs_ot_env, para_env_inter_kp, &
                                                        kpoint_preconditioned_direction, &
                                                        rayleigh_model_slope_factor, &
                                                        rayleigh_model_quadratic_factor)
         END IF
         DO ispin = 1, nspin
            tmp = DOT_PRODUCT(qs_ot_env(ispin)%ener_gx, qs_ot_env(ispin)%ener_dx)
            test_down = test_down + tmp
         END DO
      END IF
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, test_down)

      IF (test_down >= 0.0_dp) THEN ! reset to SD
         reset_to_sd = .TRUE.
         beta_pr = 0.0_dp
         IF (do_ks) THEN
            DO ispin = 1, nspin
               CALL dbcsr_add(qs_ot_env(ispin)%matrix_dx, qs_ot_env(ispin)%matrix_gx_old, &
                              alpha_scalar=beta_pr, beta_scalar=-1.0_dp)
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  CALL dbcsr_add(qs_ot_env(ispin)%matrix_dx_im, qs_ot_env(ispin)%matrix_gx_old_im, &
                                 alpha_scalar=beta_pr, beta_scalar=-1.0_dp)
               END IF
               IF (qs_ot_env(1)%settings%do_rotation) THEN
                  CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_dx, &
                                 qs_ot_env(ispin)%rot_mat_gx_old, alpha_scalar=beta_pr, beta_scalar=-1.0_dp)
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_dx_im, &
                                    qs_ot_env(ispin)%rot_mat_gx_old_im, alpha_scalar=beta_pr, beta_scalar=-1.0_dp)
                  END IF
               END IF
            END DO
         END IF
         IF (do_ener) THEN
            DO ispin = 1, nspin
               qs_ot_env(ispin)%ener_dx = -qs_ot_env(ispin)%ener_gx_old
            END DO
            IF (complex_kpoint_history) THEN
               CALL prepare_linear_kpoint_energy_direction(qs_ot_env, para_env_inter_kp, &
                                                           kpoint_preconditioned_direction, &
                                                           rayleigh_model_slope_factor, &
                                                           rayleigh_model_quadratic_factor)
            END IF
         END IF
         test_down = 0.0_dp
         IF (do_ks) THEN
            DO ispin = 1, nspin
               CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin)%matrix_dx, tmp)
               test_down = test_down + tmp
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx_im, qs_ot_env(ispin)%matrix_dx_im, tmp)
                  test_down = test_down + tmp
               END IF
               IF (qs_ot_env(1)%settings%do_rotation) THEN
                  CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx, qs_ot_env(ispin)%rot_mat_dx, tmp)
                  test_down = test_down + 0.5_dp*tmp
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx_im, qs_ot_env(ispin)%rot_mat_dx_im, tmp)
                     test_down = test_down + 0.5_dp*tmp
                  END IF
               END IF
            END DO
         END IF
         IF (do_ener) THEN
            DO ispin = 1, nspin
               test_down = test_down + DOT_PRODUCT(qs_ot_env(ispin)%ener_gx, &
                                                   qs_ot_env(ispin)%ener_dx)
            END DO
         END IF
         CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, test_down)
      END IF
      IF (complex_kpoint_history .AND. qs_ot_env(1)%schur_model_pending) THEN
         qs_ot_env(1)%schur_model_slope = rayleigh_model_slope_factor* &
                                          qs_ot_env(1)%schur_model_slope
         qs_ot_env(1)%schur_model_quadratic = rayleigh_model_quadratic_factor* &
                                              qs_ot_env(1)%schur_model_quadratic
         IF (qs_ot_env(1)%schur_closure_shadow_pending) THEN
            qs_ot_env(1)%schur_closure_shadow_baseline_quadratic = &
               rayleigh_model_quadratic_factor* &
               qs_ot_env(1)%schur_closure_shadow_baseline_quadratic
            qs_ot_env(1)%schur_closure_shadow_quadratic = rayleigh_model_quadratic_factor* &
                                                          qs_ot_env(1)%schur_closure_shadow_quadratic
         END IF
         partial_direction_slope = qs_ot_env(1)%schur_model_slope
         partial_direction_quadratic = qs_ot_env(1)%schur_model_quadratic
         CALL ot_kpoint_mermin_full_direction_model( &
            partial_direction_slope, partial_direction_quadratic, qs_ot_env(1)%gnorm, &
            -test_down, beta_pr, qs_ot_env(1)%mermin_direction_curvature, &
            qs_ot_env(1)%mermin_direction_curvature_valid, full_direction_slope, &
            full_direction_quadratic, full_direction_drop, full_direction_optimum)
         IF (logger%iter_info%print_level >= debug_print_level .AND. logger%para_env%is_source()) THEN
            WRITE (cp_logger_get_default_unit_nr(logger), '(A,1X,L1,8(1X,ES16.8))') &
               " OT Mermin full-direction model history base_s base_q final_s final_q drop opt beta prev_q:", &
               qs_ot_env(1)%mermin_direction_curvature_valid, partial_direction_slope, &
               partial_direction_quadratic, full_direction_slope, full_direction_quadratic, &
               full_direction_drop, full_direction_optimum, beta_pr, &
               qs_ot_env(1)%mermin_direction_curvature
         END IF
      END IF
      schur_history_mixed = schur_history_change .AND. beta_pr > 0.0_dp
      qs_ot_env(1)%schur_history_mixed = schur_history_mixed
      IF (logger%iter_info%print_level >= debug_print_level .AND. logger%para_env%is_source()) THEN
         WRITE (cp_logger_get_default_unit_nr(logger), '(A,8(1X,ES16.8),4(1X,I6),6(1X,L1))') &
            " OT CG beta raw used gnorm old cross down ratio rscale schur_att last flavor last neg "// &
            "hist_change mixed reset spike resid_reset:", &
            beta_raw, beta_pr, qs_ot_env(1)%gnorm, qs_ot_env(1)%gnorm_old, gnorm_cross, test_down, gnorm_ratio, &
            residual_trust_scale, &
            schur_attempt, schur_last_attempt, schur_flavor, schur_last_flavor, negative_beta, &
            schur_history_change, schur_history_mixed, reset_to_sd, history_spike_reset, residual_history_reset
      END IF
      ! since we change the direction we have to adjust the gradient
      IF (legacy_preconditioned_gradient) THEN
         qs_ot_env(1)%gradient = beta_pr*qs_ot_env(1)%gradient - qs_ot_env(1)%gnorm
      ELSE
         qs_ot_env(1)%gradient = test_down
      END IF
      qs_ot_env(1)%gnorm_old = qs_ot_env(1)%gnorm
      qs_ot_env(1)%schur_cg_attempt = schur_attempt
      qs_ot_env(1)%schur_cg_flavor = schur_flavor
      IF (kpoint_preconditioned_direction) THEN
         DO ispin = 1, nspin
            qs_ot_env(ispin)%kpoint_preconditioned_direction_ready = .FALSE.
         END DO
      END IF

      CALL timestop(handle)

   END SUBROUTINE ot_new_cg_direction

! **************************************************************************************************
!> \brief restart Pulay history from the current accepted endpoint after a singular solve
!> \param qs_ot_env local OT channels
!> \param current_index circular-history slot containing the current endpoint and SD direction
! **************************************************************************************************
   SUBROUTINE ot_diis_restart_current_anchor(qs_ot_env, current_index)
      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: qs_ot_env
      INTEGER, INTENT(IN)                                :: current_index

      INTEGER                                            :: ispin
      LOGICAL                                            :: do_ener, do_ks

      do_ks = qs_ot_env(1)%settings%ks
      do_ener = qs_ot_env(1)%settings%do_ener
      IF (current_index /= 1) THEN
         DO ispin = 1, SIZE(qs_ot_env)
            IF (do_ks) THEN
               CALL dbcsr_copy(qs_ot_env(ispin)%matrix_h_x(1)%matrix, &
                               qs_ot_env(ispin)%matrix_h_x(current_index)%matrix)
               CALL dbcsr_copy(qs_ot_env(ispin)%matrix_h_e(1)%matrix, &
                               qs_ot_env(ispin)%matrix_h_e(current_index)%matrix)
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  CALL dbcsr_copy(qs_ot_env(ispin)%matrix_h_x_im(1)%matrix, &
                                  qs_ot_env(ispin)%matrix_h_x_im(current_index)%matrix)
                  CALL dbcsr_copy(qs_ot_env(ispin)%matrix_h_e_im(1)%matrix, &
                                  qs_ot_env(ispin)%matrix_h_e_im(current_index)%matrix)
               END IF
               IF (qs_ot_env(ispin)%settings%do_rotation) THEN
                  CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_h_x(1)%matrix, &
                                  qs_ot_env(ispin)%rot_mat_h_x(current_index)%matrix)
                  CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_h_e(1)%matrix, &
                                  qs_ot_env(ispin)%rot_mat_h_e(current_index)%matrix)
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_h_x_im(1)%matrix, &
                                     qs_ot_env(ispin)%rot_mat_h_x_im(current_index)%matrix)
                     CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_h_e_im(1)%matrix, &
                                     qs_ot_env(ispin)%rot_mat_h_e_im(current_index)%matrix)
                  END IF
               END IF
            END IF
            IF (do_ener) THEN
               qs_ot_env(ispin)%ener_h_x(1, :) = qs_ot_env(ispin)%ener_h_x(current_index, :)
               qs_ot_env(ispin)%ener_h_e(1, :) = qs_ot_env(ispin)%ener_h_e(current_index, :)
            END IF
         END DO
      END IF
      qs_ot_env(1)%diis_iter = 1
      qs_ot_env(1)%ls_diis(:, :) = 0.0_dp
      qs_ot_env(1)%lss_diis(:, :) = 0.0_dp
      qs_ot_env(1)%c_diis(:) = 0.0_dp
      qs_ot_env(1)%energy_h(:) = 0.0_dp
      qs_ot_env(1)%energy_h(1) = qs_ot_env(1)%etotal

   END SUBROUTINE ot_diis_restart_current_anchor

! **************************************************************************************************
!> \brief Builds an L-BFGS direction in the fixed OT product chart.
!>
!> The selected OT preconditioner is used as the initial inverse Hessian. History slot m+1 stores
!> the previous accepted point and gradient; slots 1:m store the circular (s,y) secant history.
!> Complex orbital and rotation channels and fixed-N auxiliary energies use the same product
!> metric and inter-k-point reductions as the other OT minimizers. Curvature-rejected updates and
!> non-descent directions reset the history.
!> \param qs_ot_env OT environments for all spin channels
!> \param para_env_inter_kp communicator spanning distributed k-point groups
! **************************************************************************************************
   SUBROUTINE ot_new_lbfgs_direction(qs_ot_env, para_env_inter_kp)
      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: qs_ot_env
      TYPE(mp_para_env_type), OPTIONAL, POINTER          :: para_env_inter_kp

      CHARACTER(len=*), PARAMETER :: routineN = 'ot_new_lbfgs_direction'
      REAL(KIND=dp), PARAMETER                           :: curvature_tol = 1.0E-10_dp

      INTEGER :: coupled_h0_all, handle, i, ispin, itmp, j, k, m_history, n, newest, &
                 nhistory, nspin
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: history_index
      LOGICAL                                            :: complex_kpoint_history, do_ener, &
                                                            use_current_coupled_h0
      REAL(KIND=dp) :: beta, block_scale, g_dot_z, gamma, model_quadratic_factor, &
                       model_slope_factor, norm_weight, nvariables, ss, stq, sy, test_down, tmp, yy, yy_newest
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: alpha, rho

      CALL timeset(routineN, handle)

      IF (.NOT. qs_ot_env(1)%settings%ks) THEN
         CPABORT("MINIMIZER LBFGS currently requires OT orbital variables")
      END IF

      nspin = SIZE(qs_ot_env)
      m_history = qs_ot_env(1)%settings%diis_m
      do_ener = qs_ot_env(1)%settings%do_ener
      complex_kpoint_history = ANY(qs_ot_env(:)%has_complex_kpoint_state)

      ! Form a new secant pair relative to the preceding accepted line-search point.
      IF (qs_ot_env(1)%diis_iter > 0) THEN
         j = MOD(qs_ot_env(1)%diis_iter - 1, m_history) + 1
         ss = 0.0_dp
         sy = 0.0_dp
         yy = 0.0_dp
         DO ispin = 1, nspin
            CALL dbcsr_copy(qs_ot_env(ispin)%matrix_h_x(j)%matrix, qs_ot_env(ispin)%matrix_x)
            CALL dbcsr_add(qs_ot_env(ispin)%matrix_h_x(j)%matrix, &
                           qs_ot_env(ispin)%matrix_h_x(m_history + 1)%matrix, &
                           alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)
            CALL dbcsr_copy(qs_ot_env(ispin)%matrix_h_e(j)%matrix, qs_ot_env(ispin)%matrix_gx)
            CALL dbcsr_add(qs_ot_env(ispin)%matrix_h_e(j)%matrix, &
                           qs_ot_env(ispin)%matrix_h_e(m_history + 1)%matrix, &
                           alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)
            CALL dbcsr_dot(qs_ot_env(ispin)%matrix_h_x(j)%matrix, &
                           qs_ot_env(ispin)%matrix_h_x(j)%matrix, tmp)
            ss = ss + tmp
            CALL dbcsr_dot(qs_ot_env(ispin)%matrix_h_x(j)%matrix, &
                           qs_ot_env(ispin)%matrix_h_e(j)%matrix, tmp)
            sy = sy + tmp
            CALL dbcsr_dot(qs_ot_env(ispin)%matrix_h_e(j)%matrix, &
                           qs_ot_env(ispin)%matrix_h_e(j)%matrix, tmp)
            yy = yy + tmp
            IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
               CALL dbcsr_copy(qs_ot_env(ispin)%matrix_h_x_im(j)%matrix, &
                               qs_ot_env(ispin)%matrix_x_im)
               CALL dbcsr_add(qs_ot_env(ispin)%matrix_h_x_im(j)%matrix, &
                              qs_ot_env(ispin)%matrix_h_x_im(m_history + 1)%matrix, &
                              alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)
               CALL dbcsr_copy(qs_ot_env(ispin)%matrix_h_e_im(j)%matrix, &
                               qs_ot_env(ispin)%matrix_gx_im)
               CALL dbcsr_add(qs_ot_env(ispin)%matrix_h_e_im(j)%matrix, &
                              qs_ot_env(ispin)%matrix_h_e_im(m_history + 1)%matrix, &
                              alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)
               CALL dbcsr_dot(qs_ot_env(ispin)%matrix_h_x_im(j)%matrix, &
                              qs_ot_env(ispin)%matrix_h_x_im(j)%matrix, tmp)
               ss = ss + tmp
               CALL dbcsr_dot(qs_ot_env(ispin)%matrix_h_x_im(j)%matrix, &
                              qs_ot_env(ispin)%matrix_h_e_im(j)%matrix, tmp)
               sy = sy + tmp
               CALL dbcsr_dot(qs_ot_env(ispin)%matrix_h_e_im(j)%matrix, &
                              qs_ot_env(ispin)%matrix_h_e_im(j)%matrix, tmp)
               yy = yy + tmp
            END IF
            IF (qs_ot_env(1)%settings%do_rotation) THEN
               CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_h_x(j)%matrix, qs_ot_env(ispin)%rot_mat_x)
               CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_h_x(j)%matrix, &
                              qs_ot_env(ispin)%rot_mat_h_x(m_history + 1)%matrix, &
                              alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)
               CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_h_e(j)%matrix, qs_ot_env(ispin)%rot_mat_gx)
               CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_h_e(j)%matrix, &
                              qs_ot_env(ispin)%rot_mat_h_e(m_history + 1)%matrix, &
                              alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)
               CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_h_x(j)%matrix, &
                              qs_ot_env(ispin)%rot_mat_h_x(j)%matrix, tmp)
               ss = ss + 0.5_dp*tmp
               CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_h_x(j)%matrix, &
                              qs_ot_env(ispin)%rot_mat_h_e(j)%matrix, tmp)
               sy = sy + 0.5_dp*tmp
               CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_h_e(j)%matrix, &
                              qs_ot_env(ispin)%rot_mat_h_e(j)%matrix, tmp)
               yy = yy + 0.5_dp*tmp
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_h_x_im(j)%matrix, &
                                  qs_ot_env(ispin)%rot_mat_x_im)
                  CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_h_x_im(j)%matrix, &
                                 qs_ot_env(ispin)%rot_mat_h_x_im(m_history + 1)%matrix, &
                                 alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)
                  CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_h_e_im(j)%matrix, &
                                  qs_ot_env(ispin)%rot_mat_gx_im)
                  CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_h_e_im(j)%matrix, &
                                 qs_ot_env(ispin)%rot_mat_h_e_im(m_history + 1)%matrix, &
                                 alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)
                  CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_h_x_im(j)%matrix, &
                                 qs_ot_env(ispin)%rot_mat_h_x_im(j)%matrix, tmp)
                  ss = ss + 0.5_dp*tmp
                  CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_h_x_im(j)%matrix, &
                                 qs_ot_env(ispin)%rot_mat_h_e_im(j)%matrix, tmp)
                  sy = sy + 0.5_dp*tmp
                  CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_h_e_im(j)%matrix, &
                                 qs_ot_env(ispin)%rot_mat_h_e_im(j)%matrix, tmp)
                  yy = yy + 0.5_dp*tmp
               END IF
            END IF
            IF (do_ener) THEN
               qs_ot_env(ispin)%ener_h_x(j, :) = qs_ot_env(ispin)%ener_x(:) - &
                                                 qs_ot_env(ispin)%ener_h_x(m_history + 1, :)
               qs_ot_env(ispin)%ener_h_e(j, :) = qs_ot_env(ispin)%ener_gx(:) - &
                                                 qs_ot_env(ispin)%ener_h_e(m_history + 1, :)
               ss = ss + DOT_PRODUCT(qs_ot_env(ispin)%ener_h_x(j, :), &
                                     qs_ot_env(ispin)%ener_h_x(j, :))
               sy = sy + DOT_PRODUCT(qs_ot_env(ispin)%ener_h_x(j, :), &
                                     qs_ot_env(ispin)%ener_h_e(j, :))
               yy = yy + DOT_PRODUCT(qs_ot_env(ispin)%ener_h_e(j, :), &
                                     qs_ot_env(ispin)%ener_h_e(j, :))
            END IF
         END DO
         CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, ss)
         CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, sy)
         CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, yy)
         IF (ieee_is_finite(sy) .AND. ieee_is_finite(ss) .AND. ieee_is_finite(yy) .AND. &
             ss > 0.0_dp .AND. yy > 0.0_dp .AND. sy > curvature_tol*SQRT(ss*yy)) THEN
            qs_ot_env(1)%lbfgs_rho(j) = 1.0_dp/sy
            qs_ot_env(1)%lbfgs_yy(j) = yy
            qs_ot_env(1)%diis_iter = qs_ot_env(1)%diis_iter + 1
         ELSE
            ! The current point remains the reference, but all older pairs are discarded.
            qs_ot_env(1)%diis_iter = 1
         END IF
      ELSE
         qs_ot_env(1)%diis_iter = 1
      END IF

      DO ispin = 1, nspin
         CALL dbcsr_copy(qs_ot_env(ispin)%matrix_h_x(m_history + 1)%matrix, &
                         qs_ot_env(ispin)%matrix_x)
         CALL dbcsr_copy(qs_ot_env(ispin)%matrix_h_e(m_history + 1)%matrix, &
                         qs_ot_env(ispin)%matrix_gx)
         CALL dbcsr_copy(qs_ot_env(ispin)%matrix_gx_old, qs_ot_env(ispin)%matrix_gx)
         IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
            CALL dbcsr_copy(qs_ot_env(ispin)%matrix_h_x_im(m_history + 1)%matrix, &
                            qs_ot_env(ispin)%matrix_x_im)
            CALL dbcsr_copy(qs_ot_env(ispin)%matrix_h_e_im(m_history + 1)%matrix, &
                            qs_ot_env(ispin)%matrix_gx_im)
            CALL dbcsr_copy(qs_ot_env(ispin)%matrix_gx_old_im, qs_ot_env(ispin)%matrix_gx_im)
         END IF
         IF (qs_ot_env(1)%settings%do_rotation) THEN
            CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_h_x(m_history + 1)%matrix, &
                            qs_ot_env(ispin)%rot_mat_x)
            CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_h_e(m_history + 1)%matrix, &
                            qs_ot_env(ispin)%rot_mat_gx)
            CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_gx_old, qs_ot_env(ispin)%rot_mat_gx)
            IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
               CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_h_x_im(m_history + 1)%matrix, &
                               qs_ot_env(ispin)%rot_mat_x_im)
               CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_h_e_im(m_history + 1)%matrix, &
                               qs_ot_env(ispin)%rot_mat_gx_im)
               CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_gx_old_im, &
                               qs_ot_env(ispin)%rot_mat_gx_im)
            END IF
         END IF
         IF (do_ener) THEN
            qs_ot_env(ispin)%ener_h_x(m_history + 1, :) = qs_ot_env(ispin)%ener_x(:)
            qs_ot_env(ispin)%ener_h_e(m_history + 1, :) = qs_ot_env(ispin)%ener_gx(:)
            qs_ot_env(ispin)%ener_gx_old(:) = qs_ot_env(ispin)%ener_gx(:)
         END IF
      END DO

      nhistory = MIN(qs_ot_env(1)%diis_iter - 1, m_history)
      ALLOCATE (alpha(MAX(1, nhistory)), rho(MAX(1, nhistory)), &
                history_index(MAX(1, nhistory)))
      yy_newest = 0.0_dp

      ! First loop, newest to oldest: q <- V^T g.
      IF (nhistory > 0) THEN
         newest = MOD(qs_ot_env(1)%diis_iter - 2, m_history) + 1
         yy_newest = qs_ot_env(1)%lbfgs_yy(newest)
         DO i = 1, nhistory
            history_index(i) = MODULO(newest - i, m_history) + 1
            stq = 0.0_dp
            DO ispin = 1, nspin
               CALL dbcsr_dot(qs_ot_env(ispin)%matrix_h_x(history_index(i))%matrix, &
                              qs_ot_env(ispin)%matrix_gx_old, tmp)
               stq = stq + tmp
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  CALL dbcsr_dot(qs_ot_env(ispin)%matrix_h_x_im(history_index(i))%matrix, &
                                 qs_ot_env(ispin)%matrix_gx_old_im, tmp)
                  stq = stq + tmp
               END IF
               IF (qs_ot_env(1)%settings%do_rotation) THEN
                  CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_h_x(history_index(i))%matrix, &
                                 qs_ot_env(ispin)%rot_mat_gx_old, tmp)
                  stq = stq + 0.5_dp*tmp
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_h_x_im(history_index(i))%matrix, &
                                    qs_ot_env(ispin)%rot_mat_gx_old_im, tmp)
                     stq = stq + 0.5_dp*tmp
                  END IF
               END IF
               IF (do_ener) stq = stq + &
                                  DOT_PRODUCT(qs_ot_env(ispin)%ener_h_x(history_index(i), :), &
                                              qs_ot_env(ispin)%ener_gx_old(:))
            END DO
            CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, stq)
            rho(i) = qs_ot_env(1)%lbfgs_rho(history_index(i))
            alpha(i) = rho(i)*stq
            DO ispin = 1, nspin
               CALL dbcsr_add(qs_ot_env(ispin)%matrix_gx_old, &
                              qs_ot_env(ispin)%matrix_h_e(history_index(i))%matrix, &
                              alpha_scalar=1.0_dp, beta_scalar=-alpha(i))
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  CALL dbcsr_add(qs_ot_env(ispin)%matrix_gx_old_im, &
                                 qs_ot_env(ispin)%matrix_h_e_im(history_index(i))%matrix, &
                                 alpha_scalar=1.0_dp, beta_scalar=-alpha(i))
               END IF
               IF (qs_ot_env(1)%settings%do_rotation) THEN
                  CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_gx_old, &
                                 qs_ot_env(ispin)%rot_mat_h_e(history_index(i))%matrix, &
                                 alpha_scalar=1.0_dp, beta_scalar=-alpha(i))
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_gx_old_im, &
                                    qs_ot_env(ispin)%rot_mat_h_e_im(history_index(i))%matrix, &
                                    alpha_scalar=1.0_dp, beta_scalar=-alpha(i))
                  END IF
               END IF
               IF (do_ener) qs_ot_env(ispin)%ener_gx_old(:) = &
                  qs_ot_env(ispin)%ener_gx_old(:) - &
                  alpha(i)*qs_ot_env(ispin)%ener_h_e(history_index(i), :)
            END DO
         END DO
      END IF

      gamma = 1.0_dp
      IF (nhistory > 0 .AND. yy_newest > 0.0_dp) gamma = 1.0_dp/(rho(1)*yy_newest)
      use_current_coupled_h0 = complex_kpoint_history .AND. nhistory == 0
      DO ispin = 1, nspin
         IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
            use_current_coupled_h0 = use_current_coupled_h0 .AND. &
                                     qs_ot_env(ispin)%kpoint_preconditioned_direction_ready
         END IF
      END DO
      coupled_h0_all = MERGE(1, 0, use_current_coupled_h0)
      CALL ot_mini_min_kpoint_integer(qs_ot_env, para_env_inter_kp, coupled_h0_all)
      use_current_coupled_h0 = coupled_h0_all == 1
      CALL apply_lbfgs_h0(use_current_coupled_h0)

      ! Second loop, oldest to newest: z <- V H0 V^T g.
      DO i = nhistory, 1, -1
         beta = 0.0_dp
         DO ispin = 1, nspin
            CALL dbcsr_dot(qs_ot_env(ispin)%matrix_h_e(history_index(i))%matrix, &
                           qs_ot_env(ispin)%matrix_dx, tmp)
            beta = beta + tmp
            IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
               CALL dbcsr_dot(qs_ot_env(ispin)%matrix_h_e_im(history_index(i))%matrix, &
                              qs_ot_env(ispin)%matrix_dx_im, tmp)
               beta = beta + tmp
            END IF
            IF (qs_ot_env(1)%settings%do_rotation) THEN
               CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_h_e(history_index(i))%matrix, &
                              qs_ot_env(ispin)%rot_mat_dx, tmp)
               beta = beta + 0.5_dp*tmp
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_h_e_im(history_index(i))%matrix, &
                                 qs_ot_env(ispin)%rot_mat_dx_im, tmp)
                  beta = beta + 0.5_dp*tmp
               END IF
            END IF
            IF (do_ener) beta = beta + &
                                DOT_PRODUCT(qs_ot_env(ispin)%ener_h_e(history_index(i), :), &
                                            qs_ot_env(ispin)%ener_dx(:))
         END DO
         CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, beta)
         beta = rho(i)*beta
         DO ispin = 1, nspin
            CALL dbcsr_add(qs_ot_env(ispin)%matrix_dx, &
                           qs_ot_env(ispin)%matrix_h_x(history_index(i))%matrix, &
                           alpha_scalar=1.0_dp, beta_scalar=alpha(i) - beta)
            IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
               CALL dbcsr_add(qs_ot_env(ispin)%matrix_dx_im, &
                              qs_ot_env(ispin)%matrix_h_x_im(history_index(i))%matrix, &
                              alpha_scalar=1.0_dp, beta_scalar=alpha(i) - beta)
            END IF
            IF (qs_ot_env(1)%settings%do_rotation) THEN
               CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_dx, &
                              qs_ot_env(ispin)%rot_mat_h_x(history_index(i))%matrix, &
                              alpha_scalar=1.0_dp, beta_scalar=alpha(i) - beta)
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_dx_im, &
                                 qs_ot_env(ispin)%rot_mat_h_x_im(history_index(i))%matrix, &
                                 alpha_scalar=1.0_dp, beta_scalar=alpha(i) - beta)
               END IF
            END IF
            IF (do_ener) qs_ot_env(ispin)%ener_dx(:) = qs_ot_env(ispin)%ener_dx(:) + &
                                                       (alpha(i) - beta)*qs_ot_env(ispin)%ener_h_x(history_index(i), :)
         END DO
      END DO

      CALL scale_lbfgs_direction(-1.0_dp)
      IF (do_ener .AND. complex_kpoint_history) THEN
         CALL prepare_linear_kpoint_energy_direction(qs_ot_env, para_env_inter_kp, .FALSE., &
                                                     model_slope_factor, model_quadratic_factor)
      END IF
      CALL compute_lbfgs_direction_slope(test_down)

      ! Preserve a guaranteed descent direction even if the supplied initial preconditioner is
      ! indefinite or numerical noise corrupts the secant recursion.
      IF (.NOT. ieee_is_finite(test_down) .OR. test_down >= 0.0_dp) THEN
         qs_ot_env(1)%diis_iter = 1
         qs_ot_env(1)%OT_METHOD_FULL = "OT L-SD"
         CALL copy_lbfgs_current_gradient_to_q()
         CALL apply_lbfgs_h0(use_current_coupled_h0)
         CALL scale_lbfgs_direction(-1.0_dp)
         IF (do_ener .AND. complex_kpoint_history) THEN
            CALL prepare_linear_kpoint_energy_direction(qs_ot_env, para_env_inter_kp, .FALSE., &
                                                        model_slope_factor, model_quadratic_factor)
         END IF
         CALL compute_lbfgs_direction_slope(test_down)
         IF (.NOT. ieee_is_finite(test_down) .OR. test_down >= 0.0_dp) THEN
            DO ispin = 1, nspin
               CALL dbcsr_copy(qs_ot_env(ispin)%matrix_dx, qs_ot_env(ispin)%matrix_gx)
               CALL dbcsr_scale(qs_ot_env(ispin)%matrix_dx, -1.0_dp)
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  CALL dbcsr_copy(qs_ot_env(ispin)%matrix_dx_im, qs_ot_env(ispin)%matrix_gx_im)
                  CALL dbcsr_scale(qs_ot_env(ispin)%matrix_dx_im, -1.0_dp)
               END IF
               IF (qs_ot_env(1)%settings%do_rotation) THEN
                  CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_dx, qs_ot_env(ispin)%rot_mat_gx)
                  CALL dbcsr_scale(qs_ot_env(ispin)%rot_mat_dx, -1.0_dp)
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_dx_im, qs_ot_env(ispin)%rot_mat_gx_im)
                     CALL dbcsr_scale(qs_ot_env(ispin)%rot_mat_dx_im, -1.0_dp)
                  END IF
               END IF
               IF (do_ener) qs_ot_env(ispin)%ener_dx(:) = -qs_ot_env(ispin)%ener_gx(:)
            END DO
            IF (do_ener .AND. complex_kpoint_history) THEN
               CALL prepare_linear_kpoint_energy_direction(qs_ot_env, para_env_inter_kp, .FALSE., &
                                                           model_slope_factor, model_quadratic_factor)
            END IF
            CALL compute_lbfgs_direction_slope(test_down)
         END IF
      END IF

      g_dot_z = MAX(0.0_dp, -test_down)
      qs_ot_env(1)%gnorm = g_dot_z
      qs_ot_env(1)%gradient = test_down

      nvariables = 0.0_dp
      DO ispin = 1, nspin
         CALL dbcsr_get_info(qs_ot_env(ispin)%matrix_x, nfullrows_total=n, nfullcols_total=k)
         norm_weight = REAL(MAX(1, qs_ot_env(ispin)%kpoint_multiplicity), KIND=dp)
         nvariables = nvariables + norm_weight* &
                      REAL(INT(n, KIND=int_8)*INT(k, KIND=int_8), KIND=dp)
         IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
            nvariables = nvariables + norm_weight* &
                         REAL(INT(n, KIND=int_8)*INT(k, KIND=int_8), KIND=dp)
         END IF
         IF (qs_ot_env(1)%settings%do_rotation) THEN
            IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
               nvariables = nvariables + norm_weight* &
                            REAL(INT(k, KIND=int_8)*INT(k, KIND=int_8), KIND=dp)
            ELSE
               nvariables = nvariables + norm_weight* &
                            REAL(INT(k, KIND=int_8)*INT(k - 1, KIND=int_8)/2_int_8, KIND=dp)
            END IF
         END IF
         IF (do_ener) THEN
            itmp = SIZE(qs_ot_env(ispin)%ener_x)
            nvariables = nvariables + norm_weight*REAL(itmp, KIND=dp)
         END IF
      END DO
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, nvariables)
      IF (nvariables > 0.0_dp) THEN
         qs_ot_env(1)%delta = SQRT(g_dot_z/nvariables)
      ELSE
         qs_ot_env(1)%delta = 0.0_dp
         qs_ot_env(1)%gradient = 0.0_dp
      END IF

      IF (use_current_coupled_h0) qs_ot_env(:)%kpoint_preconditioned_direction_ready = .FALSE.

      DEALLOCATE (alpha, rho, history_index)
      CALL timestop(handle)

   CONTAINS

! **************************************************************************************************
!> \brief product-space slope g^T p, including real, imaginary, rotation, and energy blocks
!> \param slope globally reduced slope
! **************************************************************************************************
      SUBROUTINE compute_lbfgs_direction_slope(slope)
      REAL(KIND=dp), INTENT(OUT)                         :: slope

         slope = 0.0_dp
         DO ispin = 1, nspin
            CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin)%matrix_dx, tmp)
            slope = slope + tmp
            IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
               CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx_im, qs_ot_env(ispin)%matrix_dx_im, tmp)
               slope = slope + tmp
            END IF
            IF (qs_ot_env(1)%settings%do_rotation) THEN
               CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx, qs_ot_env(ispin)%rot_mat_dx, tmp)
               slope = slope + 0.5_dp*tmp
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx_im, qs_ot_env(ispin)%rot_mat_dx_im, tmp)
                  slope = slope + 0.5_dp*tmp
               END IF
            END IF
            IF (do_ener) slope = slope + &
                                 DOT_PRODUCT(qs_ot_env(ispin)%ener_gx, qs_ot_env(ispin)%ener_dx)
         END DO
         CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, slope)
      END SUBROUTINE compute_lbfgs_direction_slope

! **************************************************************************************************
!> \brief copy the current product gradient into the L-BFGS q workspaces
! **************************************************************************************************
      SUBROUTINE copy_lbfgs_current_gradient_to_q()
         DO ispin = 1, nspin
            CALL dbcsr_copy(qs_ot_env(ispin)%matrix_gx_old, qs_ot_env(ispin)%matrix_gx)
            IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
               CALL dbcsr_copy(qs_ot_env(ispin)%matrix_gx_old_im, qs_ot_env(ispin)%matrix_gx_im)
            END IF
            IF (qs_ot_env(1)%settings%do_rotation) THEN
               CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_gx_old, qs_ot_env(ispin)%rot_mat_gx)
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_gx_old_im, qs_ot_env(ispin)%rot_mat_gx_im)
               END IF
            END IF
            IF (do_ener) qs_ot_env(ispin)%ener_gx_old(:) = qs_ot_env(ispin)%ener_gx(:)
         END DO
      END SUBROUTINE copy_lbfgs_current_gradient_to_q

! **************************************************************************************************
!> \brief apply the fixed-chart initial inverse Hessian to the q workspaces
!> \param coupled_h0 use the already built coupled rotation/energy response for q=g
! **************************************************************************************************
      SUBROUTINE apply_lbfgs_h0(coupled_h0)
      LOGICAL, INTENT(IN)                                :: coupled_h0

         DO ispin = 1, nspin
            norm_weight = ot_kpoint_channel_norm_weight(qs_ot_env(ispin))
            IF (ASSOCIATED(qs_ot_env(ispin)%preconditioner)) THEN
               IF (qs_ot_env(ispin)%has_complex_kpoint_state .AND. &
                   ASSOCIATED(qs_ot_env(ispin)%preconditioner%complex_fm)) THEN
                  CALL apply_preconditioner(qs_ot_env(ispin)%preconditioner, &
                                            qs_ot_env(ispin)%matrix_gx_old, &
                                            qs_ot_env(ispin)%matrix_gx_old_im, &
                                            qs_ot_env(ispin)%matrix_dx, &
                                            qs_ot_env(ispin)%matrix_dx_im)
               ELSE
                  CALL apply_preconditioner(qs_ot_env(ispin)%preconditioner, &
                                            qs_ot_env(ispin)%matrix_gx_old, &
                                            qs_ot_env(ispin)%matrix_dx)
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL apply_preconditioner(qs_ot_env(ispin)%preconditioner, &
                                               qs_ot_env(ispin)%matrix_gx_old_im, &
                                               qs_ot_env(ispin)%matrix_dx_im)
                  END IF
               END IF
               IF (qs_ot_env(ispin)%has_complex_kpoint_state .AND. &
                   ASSOCIATED(qs_ot_env(ispin)%matrix_occ_precond_scale)) THEN
                  CALL dbcsr_scale_by_vector(qs_ot_env(ispin)%matrix_dx, &
                                             alpha=qs_ot_env(ispin)%matrix_occ_precond_scale, &
                                             side='right')
                  CALL dbcsr_scale_by_vector(qs_ot_env(ispin)%matrix_dx_im, &
                                             alpha=qs_ot_env(ispin)%matrix_occ_precond_scale, &
                                             side='right')
               END IF
               IF (qs_ot_env(ispin)%has_kpoint_context) THEN
                  CALL dbcsr_scale(qs_ot_env(ispin)%matrix_dx, norm_weight)
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_scale(qs_ot_env(ispin)%matrix_dx_im, norm_weight)
                  END IF
               END IF
               block_scale = 1.0_dp
            ELSE
               CALL dbcsr_copy(qs_ot_env(ispin)%matrix_dx, qs_ot_env(ispin)%matrix_gx_old)
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  CALL dbcsr_copy(qs_ot_env(ispin)%matrix_dx_im, qs_ot_env(ispin)%matrix_gx_old_im)
               END IF
               block_scale = gamma
               IF (qs_ot_env(ispin)%has_kpoint_context) block_scale = block_scale*norm_weight
               CALL dbcsr_scale(qs_ot_env(ispin)%matrix_dx, block_scale)
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  CALL dbcsr_scale(qs_ot_env(ispin)%matrix_dx_im, block_scale)
               END IF
            END IF

            IF (qs_ot_env(1)%settings%do_rotation) THEN
               IF (coupled_h0 .AND. qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_dx, qs_ot_env(ispin)%rot_mat_pgx)
                  CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_dx_im, qs_ot_env(ispin)%rot_mat_pgx_im)
               ELSE
                  CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_dx, qs_ot_env(ispin)%rot_mat_gx_old)
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_dx_im, &
                                     qs_ot_env(ispin)%rot_mat_gx_old_im)
                  END IF
                  block_scale = MERGE(1.0_dp, gamma, ASSOCIATED(qs_ot_env(ispin)%preconditioner))
                  IF (qs_ot_env(ispin)%has_kpoint_context) block_scale = block_scale*norm_weight
                  CALL dbcsr_scale(qs_ot_env(ispin)%rot_mat_dx, block_scale)
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_scale(qs_ot_env(ispin)%rot_mat_dx_im, block_scale)
                  END IF
               END IF
            END IF
            IF (do_ener) THEN
               IF (coupled_h0 .AND. qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  qs_ot_env(ispin)%ener_dx(:) = qs_ot_env(ispin)%ener_pgx(:)
               ELSE
                  block_scale = MERGE(1.0_dp, gamma, ASSOCIATED(qs_ot_env(ispin)%preconditioner))
                  IF (qs_ot_env(ispin)%has_kpoint_context) block_scale = block_scale*norm_weight
                  qs_ot_env(ispin)%ener_dx(:) = block_scale*qs_ot_env(ispin)%ener_gx_old(:)
               END IF
            END IF
         END DO
      END SUBROUTINE apply_lbfgs_h0

! **************************************************************************************************
!> \brief scale every component of the current L-BFGS product direction
!> \param factor common scale
! **************************************************************************************************
      SUBROUTINE scale_lbfgs_direction(factor)
      REAL(KIND=dp), INTENT(IN)                          :: factor

         DO ispin = 1, nspin
            CALL dbcsr_scale(qs_ot_env(ispin)%matrix_dx, factor)
            IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
               CALL dbcsr_scale(qs_ot_env(ispin)%matrix_dx_im, factor)
            END IF
            IF (qs_ot_env(1)%settings%do_rotation) THEN
               CALL dbcsr_scale(qs_ot_env(ispin)%rot_mat_dx, factor)
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  CALL dbcsr_scale(qs_ot_env(ispin)%rot_mat_dx_im, factor)
               END IF
            END IF
            IF (do_ener) qs_ot_env(ispin)%ener_dx(:) = factor*qs_ot_env(ispin)%ener_dx(:)
         END DO
      END SUBROUTINE scale_lbfgs_direction

   END SUBROUTINE ot_new_lbfgs_direction

! **************************************************************************************************
!> \brief ...
!> \param qs_ot_env ...
!> \param para_env_inter_kp communicator spanning distributed k-point groups
! **************************************************************************************************
   SUBROUTINE ot_diis_step(qs_ot_env, para_env_inter_kp)
      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: qs_ot_env
      TYPE(mp_para_env_type), OPTIONAL, POINTER          :: para_env_inter_kp

      CHARACTER(len=*), PARAMETER                        :: routineN = 'ot_diis_step'

      INTEGER                                            :: diis_bound, diis_m, direction_ready_all, &
                                                            group, handle, i, info, ispin, j, &
                                                            last_index, nfixed_groups, nspin
      LOGICAL :: current_state_converged, diis_trust_recovered, do_ener, do_ks, do_ot_sd, &
                 kpoint_preconditioned_direction, kpoint_preconditioned_energy, &
                 kpoint_preconditioned_rotation, kpoint_response_probe, mermin_energy_rejected, &
                 mermin_residual_rejected, mermin_step_rejected, previous_response_candidate
      REAL(KIND=dp) :: accepted_drop, accepted_physical_ratio, accepted_residual_ratio, &
                       contracted_step, current_exit_delta, current_physical_delta, current_rotation_gnorm, &
                       norm_weight, nvar, overlap, previous_response_model_work, previous_step, recovered_step, &
                       reference_energy, reference_physical_delta, reference_rotation_gnorm, rejected_energy, &
                       response_blend, response_model_work, rotation_raw_slope, rotation_response_slope, tmp, &
                       tr_xnew_gx, tr_xold_gx
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: rotation_group_raw_slope, &
                                                            rotation_group_response_slope
      TYPE(cp_logger_type), POINTER                      :: logger

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()

      do_ks = qs_ot_env(1)%settings%ks
      do_ener = qs_ot_env(1)%settings%do_ener
      nspin = SIZE(qs_ot_env)

      kpoint_preconditioned_direction = .FALSE.
      IF (ANY(qs_ot_env(:)%has_complex_kpoint_state)) THEN
         kpoint_preconditioned_direction = .TRUE.
         DO ispin = 1, nspin
            IF (qs_ot_env(ispin)%has_complex_kpoint_state .AND. &
                .NOT. qs_ot_env(ispin)%kpoint_preconditioned_direction_ready) THEN
               kpoint_preconditioned_direction = .FALSE.
            END IF
         END DO
      END IF
      direction_ready_all = MERGE(1, 0, kpoint_preconditioned_direction)
      CALL ot_mini_min_kpoint_integer(qs_ot_env, para_env_inter_kp, direction_ready_all)
      kpoint_preconditioned_direction = direction_ready_all == 1
      ! Keep the raw complex rotation residual in the Pulay metric.  The
      ! sign-aware Schur response is admitted separately as the current
      ! coordinate correction, so indefinite directions never enter history.
      nfixed_groups = qs_ot_env(1)%kpoint_physical_delta_groups
      ALLOCATE (rotation_group_raw_slope(MAX(1, nfixed_groups)))
      ALLOCATE (rotation_group_response_slope(MAX(1, nfixed_groups)))
      rotation_group_raw_slope(:) = 0.0_dp
      rotation_group_response_slope(:) = 0.0_dp
      rotation_raw_slope = 0.0_dp
      rotation_response_slope = 0.0_dp
      IF (kpoint_preconditioned_direction .AND. nfixed_groups > 0 .AND. &
          qs_ot_env(1)%settings%do_rotation) THEN
         DO ispin = 1, nspin
            group = qs_ot_env(ispin)%occupation_group_index
            CPASSERT(group >= 1 .AND. group <= nfixed_groups)
            norm_weight = ot_kpoint_channel_norm_weight(qs_ot_env(ispin))
            CPASSERT(ASSOCIATED(qs_ot_env(ispin)%rot_mat_pgx))
            CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx, qs_ot_env(ispin)%rot_mat_gx, tmp)
            rotation_group_raw_slope(group) = &
               rotation_group_raw_slope(group) + 0.5_dp*norm_weight*tmp
            CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx, qs_ot_env(ispin)%rot_mat_pgx, tmp)
            rotation_group_response_slope(group) = &
               rotation_group_response_slope(group) + 0.5_dp*tmp
            IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
               CPASSERT(ASSOCIATED(qs_ot_env(ispin)%rot_mat_pgx_im))
               CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx_im, &
                              qs_ot_env(ispin)%rot_mat_gx_im, tmp)
               rotation_group_raw_slope(group) = &
                  rotation_group_raw_slope(group) + 0.5_dp*norm_weight*tmp
               CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx_im, &
                              qs_ot_env(ispin)%rot_mat_pgx_im, tmp)
               rotation_group_response_slope(group) = &
                  rotation_group_response_slope(group) + 0.5_dp*tmp
            END IF
         END DO
         DO group = 1, nfixed_groups
            CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, rotation_group_raw_slope(group))
            CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, rotation_group_response_slope(group))
         END DO
         rotation_raw_slope = SUM(rotation_group_raw_slope(1:nfixed_groups))
         rotation_response_slope = SUM(rotation_group_response_slope(1:nfixed_groups))
      END IF
      ! DIIS candidates are protected by the accepted-Mermin audit.  The SCF driver therefore
      ! monitors the frontier without externally replacing these coordinates.
      kpoint_preconditioned_rotation = ot_kpoint_diis_rotation_response_ready( &
                                       kpoint_preconditioned_direction .AND. nfixed_groups > 0, &
                                       rotation_group_response_slope(1:MAX(1, nfixed_groups)), &
                                       rotation_group_raw_slope(1:MAX(1, nfixed_groups)), &
                                       nfixed_groups > 1)
      kpoint_response_probe = .FALSE.
      kpoint_preconditioned_energy = kpoint_preconditioned_direction
      IF (logger%iter_info%print_level >= debug_print_level .AND. logger%para_env%is_source()) THEN
         WRITE (cp_logger_get_default_unit_nr(logger), '(A,1X,L1,2(1X,ES16.8))') &
            " OT DIIS current rotation response ready raw_slope response_slope:", &
            kpoint_preconditioned_rotation, rotation_raw_slope, rotation_response_slope
         DO group = 1, nfixed_groups
            WRITE (cp_logger_get_default_unit_nr(logger), '(A,1X,I4,2(1X,ES16.8))') &
               " OT DIIS fixed-N rotation response group raw_slope response_slope:", &
               group, rotation_group_raw_slope(group), rotation_group_response_slope(group)
         END DO
      END IF

      diis_m = qs_ot_env(1)%settings%diis_m

      mermin_step_rejected = .FALSE.
      mermin_energy_rejected = .FALSE.
      mermin_residual_rejected = .FALSE.
      diis_trust_recovered = .FALSE.
      accepted_drop = -1.0_dp
      accepted_physical_ratio = -1.0_dp
      accepted_residual_ratio = -1.0_dp
      current_physical_delta = -1.0_dp
      current_rotation_gnorm = -1.0_dp
      reference_physical_delta = -1.0_dp
      reference_rotation_gnorm = -1.0_dp
      response_blend = 1.0_dp
      IF (qs_ot_env(1)%has_complex_kpoint_state .AND. &
          qs_ot_env(1)%settings%occupation_preconditioner .AND. &
          qs_ot_env(1)%settings%do_rotation .AND. do_ener .AND. &
          qs_ot_env(1)%diis_iter > 0) THEN
         last_index = MOD(qs_ot_env(1)%diis_iter - 1, diis_m) + 1
         reference_energy = qs_ot_env(1)%energy_h(last_index)
         rejected_energy = qs_ot_env(1)%etotal
         current_physical_delta = ot_kpoint_physical_exit_delta( &
                                  qs_ot_env(1)%kpoint_physical_delta, &
                                  qs_ot_env(1)%kpoint_complete_physical_delta)
         reference_physical_delta = qs_ot_env(1)%mermin_ls_reference_delta
         current_rotation_gnorm = qs_ot_env(1)%kpoint_physical_rot_gnorm
         reference_rotation_gnorm = qs_ot_env(1)%diis_reference_rotation_gnorm
         mermin_energy_rejected = &
            ot_kpoint_diis_mermin_step_rejected(rejected_energy, reference_energy)
         mermin_residual_rejected = ot_kpoint_diis_residual_step_rejected( &
                                    rejected_energy, reference_energy, current_physical_delta, &
                                    reference_physical_delta, qs_ot_env(1)%kpoint_eps_scf)
         mermin_step_rejected = mermin_energy_rejected .OR. mermin_residual_rejected
         accepted_drop = MAX(0.0_dp, reference_energy - rejected_energy)
         IF (.NOT. mermin_step_rejected .AND. current_physical_delta > 0.0_dp .AND. &
             reference_physical_delta > 0.0_dp) THEN
            accepted_physical_ratio = current_physical_delta/reference_physical_delta
         END IF
         IF (.NOT. mermin_step_rejected .AND. current_rotation_gnorm > 0.0_dp .AND. &
             reference_rotation_gnorm > 0.0_dp) THEN
            accepted_residual_ratio = SQRT(current_rotation_gnorm/reference_rotation_gnorm)
            IF (qs_ot_env(1)%diis_last_candidate_response) THEN
               qs_ot_env(1)%diis_response_residual_ratio = &
                  ot_kpoint_diis_response_endpoint_ratio( &
                  qs_ot_env(1)%diis_raw_residual_ratio, accepted_residual_ratio, &
                  qs_ot_env(1)%diis_last_response_blend)
            ELSE
               qs_ot_env(1)%diis_raw_residual_ratio = accepted_residual_ratio
            END IF
         END IF
      END IF
      previous_step = qs_ot_env(1)%ds_min
      previous_response_model_work = qs_ot_env(1)%diis_last_response_model_work
      previous_response_candidate = qs_ot_env(1)%diis_last_candidate_response
      recovered_step = ot_kpoint_diis_accepted_response_step( &
                       previous_step, qs_ot_env(1)%settings%ds_min, accepted_drop, &
                       previous_response_model_work, accepted_physical_ratio, &
                       accepted_residual_ratio, previous_response_candidate .AND. &
                       .NOT. mermin_step_rejected)
      diis_trust_recovered = ABS(recovered_step - previous_step) > &
                             10.0_dp*EPSILON(previous_step)* &
                             MAX(1.0_dp, ABS(previous_step), ABS(recovered_step))
      IF (diis_trust_recovered) THEN
         DO ispin = 1, nspin
            qs_ot_env(ispin)%ds_min = recovered_step
            qs_ot_env(ispin)%diis_last_candidate_sd = .TRUE.
            qs_ot_env(ispin)%diis_last_candidate_response = .FALSE.
            qs_ot_env(ispin)%diis_last_response_blend = 0.0_dp
            qs_ot_env(ispin)%diis_last_response_model_work = 0.0_dp
         END DO
         ! The stored Pulay errors carry the old trust scale.  Rebuild history
         ! from this accepted endpoint before using the recovered base step.
         qs_ot_env(1)%diis_iter = 0
         qs_ot_env(1)%ls_diis(:, :) = 0.0_dp
         qs_ot_env(1)%lss_diis(:, :) = 0.0_dp
         qs_ot_env(1)%c_diis(:) = 0.0_dp
         qs_ot_env(1)%energy_h(:) = 0.0_dp
      END IF
      IF (.NOT. mermin_step_rejected .AND. .NOT. diis_trust_recovered .AND. &
          nfixed_groups > 1) THEN
         kpoint_response_probe = ot_kpoint_diis_rotation_response_probe_ready( &
                                 kpoint_preconditioned_direction, &
                                 rotation_group_response_slope(1:nfixed_groups), &
                                 rotation_group_raw_slope(1:nfixed_groups), &
                                 qs_ot_env(1)%kpoint_subspace_accepted_samples, accepted_drop, &
                                 qs_ot_env(1)%ds_min, current_physical_delta, &
                                 qs_ot_env(1)%kpoint_eps_scf)
         IF (kpoint_response_probe) THEN
            response_blend = ot_kpoint_diis_response_residual_blend( &
                             qs_ot_env(1)%diis_raw_residual_ratio, &
                             qs_ot_env(1)%diis_response_residual_ratio)
         END IF
         kpoint_preconditioned_rotation = kpoint_response_probe .AND. response_blend > 0.0_dp
      END IF
      response_model_work = ABS(qs_ot_env(1)%ds_min)* &
                            (rotation_raw_slope + response_blend* &
                             (rotation_response_slope - rotation_raw_slope))
      IF (logger%iter_info%print_level >= debug_print_level .AND. logger%para_env%is_source()) THEN
         WRITE (cp_logger_get_default_unit_nr(logger), '(A,2(1X,L1),4(1X,ES16.8),1X,I6)') &
            " OT DIIS rotation response probe ready active adrop model_work total_metric "// &
            "rotation_metric samples:", &
            kpoint_response_probe, kpoint_preconditioned_rotation, accepted_drop, &
            response_model_work, current_physical_delta, current_rotation_gnorm, &
            qs_ot_env(1)%kpoint_subspace_accepted_samples
         WRITE (cp_logger_get_default_unit_nr(logger), '(A,5(1X,ES16.8),1X,L1)') &
            " OT DIIS response residual calibration accepted raw response blend last_blend last_response:", &
            accepted_residual_ratio, qs_ot_env(1)%diis_raw_residual_ratio, &
            qs_ot_env(1)%diis_response_residual_ratio, response_blend, &
            qs_ot_env(1)%diis_last_response_blend, &
            qs_ot_env(1)%diis_last_candidate_response
         WRITE (cp_logger_get_default_unit_nr(logger), '(A,2(1X,L1),6(1X,ES16.8))') &
            " OT DIIS trust recovery previous_response accepted old_step new_step model_work "// &
            "physical_ratio rotation_ratio:", &
            previous_response_candidate, diis_trust_recovered, accepted_drop, &
            previous_step, recovered_step, previous_response_model_work, &
            accepted_physical_ratio, accepted_residual_ratio
      END IF
      IF (mermin_step_rejected) THEN
         previous_step = qs_ot_env(1)%ds_min
         contracted_step = ot_kpoint_diis_rejected_step( &
                           previous_step, qs_ot_env(1)%settings%ds_min, &
                           qs_ot_env(1)%settings%eps_taylor, &
                           qs_ot_env(1)%diis_last_candidate_sd)
         IF (do_ks) THEN
            DO ispin = 1, nspin
               CALL dbcsr_copy(qs_ot_env(ispin)%matrix_x, &
                               qs_ot_env(ispin)%matrix_h_x(last_index)%matrix)
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  CALL dbcsr_copy(qs_ot_env(ispin)%matrix_x_im, &
                                  qs_ot_env(ispin)%matrix_h_x_im(last_index)%matrix)
               END IF
               IF (qs_ot_env(ispin)%settings%do_rotation) THEN
                  CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_x, &
                                  qs_ot_env(ispin)%rot_mat_h_x(last_index)%matrix)
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_x_im, &
                                     qs_ot_env(ispin)%rot_mat_h_x_im(last_index)%matrix)
                  END IF
               END IF
            END DO
         END IF
         IF (do_ener) THEN
            DO ispin = 1, nspin
               qs_ot_env(ispin)%ener_x(:) = qs_ot_env(ispin)%ener_h_x(last_index, :)
            END DO
         END IF
         DO ispin = 1, nspin
            qs_ot_env(ispin)%ds_min = contracted_step
            qs_ot_env(ispin)%diis_last_candidate_sd = .TRUE.
            qs_ot_env(ispin)%diis_last_candidate_response = .FALSE.
            qs_ot_env(ispin)%diis_last_response_blend = 0.0_dp
            qs_ot_env(ispin)%diis_last_response_model_work = 0.0_dp
         END DO
         qs_ot_env(1)%diis_iter = 0
         qs_ot_env(1)%ls_diis(:, :) = 0.0_dp
         qs_ot_env(1)%lss_diis(:, :) = 0.0_dp
         qs_ot_env(1)%c_diis(:) = 0.0_dp
         qs_ot_env(1)%energy_h(:) = 0.0_dp
         qs_ot_env(1)%kpoint_physical_delta_ready = .FALSE.
         qs_ot_env(1)%kpoint_density_stationarity_ready = .FALSE.
         qs_ot_env(1)%kpoint_complete_physical_delta = -1.0_dp
         IF (qs_ot_env(1)%kpoint_eps_scf > 0.0_dp) THEN
            qs_ot_env(1)%delta = MAX(qs_ot_env(1)%delta, &
                                     qs_ot_env(1)%kpoint_eps_scf* &
                                     (1.0_dp + 10.0_dp*EPSILON(1.0_dp)))
         END IF
         qs_ot_env(1)%OT_METHOD_FULL = "OT REJ"
         IF (logger%iter_info%print_level >= debug_print_level .AND. logger%para_env%is_source()) THEN
            WRITE (cp_logger_get_default_unit_nr(logger), '(A,6(1X,ES16.8),2(1X,L1))') &
               " OT DIIS Mermin candidate rejected energy reference residual reference_residual "// &
               "old_step new_step energy_reject residual_reject:", &
               rejected_energy, reference_energy, current_physical_delta, reference_physical_delta, &
               previous_step, contracted_step, mermin_energy_rejected, mermin_residual_rejected
         END IF
         CALL timestop(handle)
         RETURN
      END IF

      IF (qs_ot_env(1)%has_complex_kpoint_state .AND. &
          qs_ot_env(1)%settings%occupation_preconditioner .AND. &
          qs_ot_env(1)%settings%do_rotation .AND. do_ener .AND. &
          qs_ot_env(1)%diis_iter > 0) THEN
         qs_ot_env(1)%kpoint_subspace_accepted_samples = &
            MIN(qs_ot_env(1)%kpoint_subspace_accepted_samples + 1, &
                HUGE(qs_ot_env(1)%kpoint_subspace_accepted_samples) - 1)
      END IF

      IF (qs_ot_env(1)%diis_iter < diis_m) THEN
         diis_bound = qs_ot_env(1)%diis_iter + 1
      ELSE
         diis_bound = diis_m
      END IF

      j = MOD(qs_ot_env(1)%diis_iter, diis_m) + 1 ! index in the circular array
      qs_ot_env(1)%energy_h(j) = qs_ot_env(1)%etotal
      ! DIIS has no line search, so reuse its restart-serialized accepted-endpoint residual slot.
      qs_ot_env(1)%mermin_ls_reference_delta = ot_kpoint_physical_exit_delta( &
                                               qs_ot_env(1)%kpoint_physical_delta, &
                                               qs_ot_env(1)%kpoint_complete_physical_delta)
      qs_ot_env(1)%diis_reference_rotation_gnorm = qs_ot_env(1)%kpoint_physical_rot_gnorm

      ! copy the position and the error vector in the diis buffers

      IF (do_ks) THEN
         DO ispin = 1, nspin
            CALL dbcsr_copy(qs_ot_env(ispin)%matrix_h_x(j)%matrix, qs_ot_env(ispin)%matrix_x)
            IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
               CALL dbcsr_copy(qs_ot_env(ispin)%matrix_h_x_im(j)%matrix, qs_ot_env(ispin)%matrix_x_im)
            END IF
            IF (qs_ot_env(ispin)%settings%do_rotation) THEN
               CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_h_x(j)%matrix, qs_ot_env(ispin)%rot_mat_x)
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_h_x_im(j)%matrix, qs_ot_env(ispin)%rot_mat_x_im)
               END IF
            END IF
         END DO
      END IF
      IF (do_ener) THEN
         DO ispin = 1, nspin
            qs_ot_env(ispin)%ener_h_x(j, :) = qs_ot_env(ispin)%ener_x(:)
         END DO
      END IF
      IF (ASSOCIATED(qs_ot_env(1)%preconditioner)) THEN
         qs_ot_env(1)%gnorm = 0.0_dp
         IF (do_ks) THEN
            DO ispin = 1, nspin
               norm_weight = ot_kpoint_channel_norm_weight(qs_ot_env(ispin))
               IF (qs_ot_env(ispin)%has_complex_kpoint_state .AND. &
                   ASSOCIATED(qs_ot_env(ispin)%preconditioner%complex_fm)) THEN
                  CALL apply_preconditioner(qs_ot_env(ispin)%preconditioner, &
                                            qs_ot_env(ispin)%matrix_gx, &
                                            qs_ot_env(ispin)%matrix_gx_im, &
                                            qs_ot_env(ispin)%matrix_h_e(j)%matrix, &
                                            qs_ot_env(ispin)%matrix_h_e_im(j)%matrix)
               ELSE
                  CALL apply_preconditioner(qs_ot_env(ispin)%preconditioner, &
                                            qs_ot_env(ispin)%matrix_gx, &
                                            qs_ot_env(ispin)%matrix_h_e(j)%matrix)
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL apply_preconditioner(qs_ot_env(ispin)%preconditioner, &
                                               qs_ot_env(ispin)%matrix_gx_im, &
                                               qs_ot_env(ispin)%matrix_h_e_im(j)%matrix)
                  END IF
               END IF
               IF (qs_ot_env(ispin)%has_complex_kpoint_state .AND. &
                   ASSOCIATED(qs_ot_env(ispin)%matrix_occ_precond_scale)) THEN
                  CALL dbcsr_scale_by_vector(qs_ot_env(ispin)%matrix_h_e(j)%matrix, &
                                             alpha=qs_ot_env(ispin)%matrix_occ_precond_scale, &
                                             side='right')
                  CALL dbcsr_scale_by_vector(qs_ot_env(ispin)%matrix_h_e_im(j)%matrix, &
                                             alpha=qs_ot_env(ispin)%matrix_occ_precond_scale, &
                                             side='right')
               END IF
               IF (qs_ot_env(ispin)%has_kpoint_context) THEN
                  CALL dbcsr_scale(qs_ot_env(ispin)%matrix_h_e(j)%matrix, norm_weight)
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_scale(qs_ot_env(ispin)%matrix_h_e_im(j)%matrix, norm_weight)
                  END IF
               END IF
               CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin)%matrix_h_e(j)%matrix, &
                              tmp)
               qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + tmp
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx_im, &
                                 qs_ot_env(ispin)%matrix_h_e_im(j)%matrix, tmp)
                  qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + tmp
               END IF
            END DO
            IF (qs_ot_env(1)%gnorm < 0.0_dp) THEN
               WRITE (cp_logger_get_default_unit_nr(logger), *) "WARNING Preconditioner not positive definite !"
            END IF
            DO ispin = 1, nspin
               CALL dbcsr_scale(qs_ot_env(ispin)%matrix_h_e(j)%matrix, -qs_ot_env(1)%ds_min)
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  CALL dbcsr_scale(qs_ot_env(ispin)%matrix_h_e_im(j)%matrix, -qs_ot_env(1)%ds_min)
               END IF
            END DO
            IF (qs_ot_env(1)%settings%do_rotation) THEN
               DO ispin = 1, nspin
                  norm_weight = ot_kpoint_channel_norm_weight(qs_ot_env(ispin))
                  CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_h_e(j)%matrix, &
                                  qs_ot_env(ispin)%rot_mat_gx)
                  IF (qs_ot_env(ispin)%has_kpoint_context) THEN
                     CALL dbcsr_scale(qs_ot_env(ispin)%rot_mat_h_e(j)%matrix, norm_weight)
                  END IF
                  CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx, qs_ot_env(ispin)%rot_mat_h_e(j)%matrix, &
                                 tmp)
                  qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + 0.5_dp*tmp
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_h_e_im(j)%matrix, &
                                     qs_ot_env(ispin)%rot_mat_gx_im)
                     CALL dbcsr_scale(qs_ot_env(ispin)%rot_mat_h_e_im(j)%matrix, norm_weight)
                     CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx_im, &
                                    qs_ot_env(ispin)%rot_mat_h_e_im(j)%matrix, tmp)
                     qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + 0.5_dp*tmp
                  END IF
               END DO
               DO ispin = 1, nspin
                  CALL dbcsr_scale(qs_ot_env(ispin)%rot_mat_h_e(j)%matrix, -qs_ot_env(1)%ds_min)
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_scale(qs_ot_env(ispin)%rot_mat_h_e_im(j)%matrix, -qs_ot_env(1)%ds_min)
                  END IF
               END DO
            END IF
         END IF
         IF (do_ener) THEN
            DO ispin = 1, nspin
               norm_weight = ot_kpoint_channel_norm_weight(qs_ot_env(ispin))
               IF (kpoint_preconditioned_energy) THEN
                  CPASSERT(ASSOCIATED(qs_ot_env(ispin)%ener_pgx))
                  qs_ot_env(ispin)%ener_h_e(j, :) = qs_ot_env(ispin)%ener_pgx(:)
               ELSE
                  qs_ot_env(ispin)%ener_h_e(j, :) = qs_ot_env(ispin)%ener_gx(:)
                  IF (qs_ot_env(ispin)%has_kpoint_context) THEN
                     qs_ot_env(ispin)%ener_h_e(j, :) = norm_weight*qs_ot_env(ispin)%ener_h_e(j, :)
                  END IF
               END IF
               tmp = DOT_PRODUCT(qs_ot_env(ispin)%ener_h_e(j, :), qs_ot_env(ispin)%ener_gx(:))
               qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + tmp
               qs_ot_env(ispin)%ener_h_e(j, :) = -qs_ot_env(1)%ds_min*qs_ot_env(ispin)%ener_h_e(j, :)
            END DO
         END IF
      ELSE
         qs_ot_env(1)%gnorm = 0.0_dp
         IF (do_ks) THEN
            DO ispin = 1, nspin
               norm_weight = ot_kpoint_channel_norm_weight(qs_ot_env(ispin))
               CALL dbcsr_copy(qs_ot_env(ispin)%matrix_h_e(j)%matrix, qs_ot_env(ispin)%matrix_gx)
               IF (qs_ot_env(ispin)%has_kpoint_context) THEN
                  CALL dbcsr_scale(qs_ot_env(ispin)%matrix_h_e(j)%matrix, norm_weight)
               END IF
               CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin)%matrix_h_e(j)%matrix, tmp)
               qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + tmp
               CALL dbcsr_scale(qs_ot_env(ispin)%matrix_h_e(j)%matrix, -qs_ot_env(1)%ds_min)
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  CALL dbcsr_copy(qs_ot_env(ispin)%matrix_h_e_im(j)%matrix, &
                                  qs_ot_env(ispin)%matrix_gx_im)
                  CALL dbcsr_scale(qs_ot_env(ispin)%matrix_h_e_im(j)%matrix, norm_weight)
                  CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx_im, &
                                 qs_ot_env(ispin)%matrix_h_e_im(j)%matrix, tmp)
                  qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + tmp
                  CALL dbcsr_scale(qs_ot_env(ispin)%matrix_h_e_im(j)%matrix, -qs_ot_env(1)%ds_min)
               END IF
            END DO
            IF (qs_ot_env(1)%settings%do_rotation) THEN
               DO ispin = 1, nspin
                  norm_weight = ot_kpoint_channel_norm_weight(qs_ot_env(ispin))
                  CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_h_e(j)%matrix, qs_ot_env(ispin)%rot_mat_gx)
                  IF (qs_ot_env(ispin)%has_kpoint_context) THEN
                     CALL dbcsr_scale(qs_ot_env(ispin)%rot_mat_h_e(j)%matrix, norm_weight)
                  END IF
                  CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx, &
                                 qs_ot_env(ispin)%rot_mat_h_e(j)%matrix, tmp)
                  qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + 0.5_dp*tmp
                  CALL dbcsr_scale(qs_ot_env(ispin)%rot_mat_h_e(j)%matrix, -qs_ot_env(1)%ds_min)
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_h_e_im(j)%matrix, &
                                     qs_ot_env(ispin)%rot_mat_gx_im)
                     CALL dbcsr_scale(qs_ot_env(ispin)%rot_mat_h_e_im(j)%matrix, norm_weight)
                     CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx_im, &
                                    qs_ot_env(ispin)%rot_mat_h_e_im(j)%matrix, tmp)
                     qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + 0.5_dp*tmp
                     CALL dbcsr_scale(qs_ot_env(ispin)%rot_mat_h_e_im(j)%matrix, -qs_ot_env(1)%ds_min)
                  END IF
               END DO
            END IF
         END IF
         IF (do_ener) THEN
            DO ispin = 1, nspin
               norm_weight = ot_kpoint_channel_norm_weight(qs_ot_env(ispin))
               qs_ot_env(ispin)%ener_h_e(j, :) = qs_ot_env(ispin)%ener_gx(:)
               IF (qs_ot_env(ispin)%has_kpoint_context) THEN
                  qs_ot_env(ispin)%ener_h_e(j, :) = norm_weight*qs_ot_env(ispin)%ener_h_e(j, :)
               END IF
               tmp = DOT_PRODUCT(qs_ot_env(ispin)%ener_gx(:), qs_ot_env(ispin)%ener_h_e(j, :))
               qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + tmp
               qs_ot_env(ispin)%ener_h_e(j, :) = -qs_ot_env(1)%ds_min*qs_ot_env(ispin)%ener_h_e(j, :)
            END DO
         END IF
      END IF
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, qs_ot_env(1)%gnorm)
      nvar = ot_kpoint_channel_nvar(qs_ot_env, do_ks, do_ener)
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, nvar)
      ! Handling the case of no free variables to optimize
      IF (nvar > 0.0_dp) THEN
         qs_ot_env(1)%delta = SQRT(ABS(qs_ot_env(1)%gnorm)/nvar)
         qs_ot_env(1)%gradient = -qs_ot_env(1)%gnorm
      ELSE
         qs_ot_env(1)%delta = 0.0_dp
         qs_ot_env(1)%gradient = 0.0_dp
      END IF

      ! make the diis matrix and solve it
      DO i = 1, diis_bound
         ! I think there are two possible options, with and without preconditioner
         ! as a metric
         ! the second option seems most logical to me, and it seems marginally faster
         ! in some of the tests
         IF (.FALSE.) THEN
            qs_ot_env(1)%ls_diis(i, j) = 0.0_dp
            IF (do_ks) THEN
               DO ispin = 1, nspin
                  norm_weight = ot_kpoint_channel_norm_weight(qs_ot_env(ispin))
                  CALL dbcsr_dot(qs_ot_env(ispin)%matrix_h_e(j)%matrix, &
                                 qs_ot_env(ispin)%matrix_h_e(i)%matrix, &
                                 tmp)
                  qs_ot_env(1)%ls_diis(i, j) = qs_ot_env(1)%ls_diis(i, j) + norm_weight*tmp
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_dot(qs_ot_env(ispin)%matrix_h_e_im(j)%matrix, &
                                    qs_ot_env(ispin)%matrix_h_e_im(i)%matrix, tmp)
                     qs_ot_env(1)%ls_diis(i, j) = qs_ot_env(1)%ls_diis(i, j) + norm_weight*tmp
                  END IF
                  IF (qs_ot_env(ispin)%settings%do_rotation) THEN
                     CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_h_e(j)%matrix, &
                                    qs_ot_env(ispin)%rot_mat_h_e(i)%matrix, &
                                    tmp)
                     qs_ot_env(1)%ls_diis(i, j) = qs_ot_env(1)%ls_diis(i, j) + &
                                                  norm_weight*0.5_dp*tmp
                     IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                        CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_h_e_im(j)%matrix, &
                                       qs_ot_env(ispin)%rot_mat_h_e_im(i)%matrix, tmp)
                        qs_ot_env(1)%ls_diis(i, j) = qs_ot_env(1)%ls_diis(i, j) + &
                                                     norm_weight*0.5_dp*tmp
                     END IF
                  END IF
               END DO
            END IF
            IF (do_ener) THEN
               DO ispin = 1, nspin
                  norm_weight = ot_kpoint_channel_norm_weight(qs_ot_env(ispin))
                  tmp = DOT_PRODUCT(qs_ot_env(ispin)%ener_h_e(j, :), qs_ot_env(ispin)%ener_h_e(i, :))
                  qs_ot_env(1)%ls_diis(i, j) = qs_ot_env(1)%ls_diis(i, j) + norm_weight*tmp
               END DO
            END IF
         ELSE
            qs_ot_env(1)%ls_diis(i, j) = 0.0_dp
            IF (do_ks) THEN
               DO ispin = 1, nspin
                  CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx, &
                                 qs_ot_env(ispin)%matrix_h_e(i)%matrix, &
                                 tmp)
                  qs_ot_env(1)%ls_diis(i, j) = qs_ot_env(1)%ls_diis(i, j) - &
                                               qs_ot_env(1)%ds_min*tmp
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx_im, &
                                    qs_ot_env(ispin)%matrix_h_e_im(i)%matrix, tmp)
                     qs_ot_env(1)%ls_diis(i, j) = qs_ot_env(1)%ls_diis(i, j) - &
                                                  qs_ot_env(1)%ds_min*tmp
                  END IF
                  IF (qs_ot_env(ispin)%settings%do_rotation) THEN
                     CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx, &
                                    qs_ot_env(ispin)%rot_mat_h_e(i)%matrix, &
                                    tmp)
                     qs_ot_env(1)%ls_diis(i, j) = qs_ot_env(1)%ls_diis(i, j) - &
                                                  qs_ot_env(1)%ds_min*0.5_dp*tmp
                     IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                        CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx_im, &
                                       qs_ot_env(ispin)%rot_mat_h_e_im(i)%matrix, tmp)
                        qs_ot_env(1)%ls_diis(i, j) = qs_ot_env(1)%ls_diis(i, j) - &
                                                     qs_ot_env(1)%ds_min*0.5_dp*tmp
                     END IF
                  END IF
               END DO
            END IF
            IF (do_ener) THEN
               DO ispin = 1, nspin
                  tmp = DOT_PRODUCT(qs_ot_env(ispin)%ener_gx(:), qs_ot_env(ispin)%ener_h_e(i, :))
                  qs_ot_env(1)%ls_diis(i, j) = qs_ot_env(1)%ls_diis(i, j) - &
                                               qs_ot_env(1)%ds_min*tmp
               END DO
            END IF
         END IF
         CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, qs_ot_env(1)%ls_diis(i, j))
         qs_ot_env(1)%ls_diis(j, i) = qs_ot_env(1)%ls_diis(i, j)
         qs_ot_env(1)%ls_diis(i, diis_bound + 1) = 1.0_dp
         qs_ot_env(1)%ls_diis(diis_bound + 1, i) = 1.0_dp
         qs_ot_env(1)%c_diis(i) = 0.0_dp
      END DO
      qs_ot_env(1)%ls_diis(diis_bound + 1, diis_bound + 1) = 0.0_dp
      qs_ot_env(1)%c_diis(diis_bound + 1) = 1.0_dp

      current_state_converged = .FALSE.
      current_exit_delta = HUGE(1.0_dp)
      IF (qs_ot_env(1)%has_complex_kpoint_state .AND. qs_ot_env(1)%kpoint_eps_scf > 0.0_dp) THEN
         CALL ot_kpoint_mermin_endpoint_status( &
            qs_ot_env(1)%delta, ot_kpoint_physical_exit_delta( &
            qs_ot_env(1)%kpoint_physical_delta, qs_ot_env(1)%kpoint_complete_physical_delta), &
            qs_ot_env(1)%kpoint_physical_delta_ready, qs_ot_env(1)%kpoint_projected_delta, &
            qs_ot_env(1)%kpoint_density_stationarity_delta, &
            qs_ot_env(1)%kpoint_density_stationarity_ready, qs_ot_env(1)%kpoint_eps_scf, &
            current_exit_delta, current_state_converged)
      END IF
      IF (current_state_converged) THEN
         ! Retain the audited endpoint in DIIS history, but do not generate a
         ! post-convergence state that could later be written to a restart.
         qs_ot_env(1)%OT_METHOD_FULL = "OT FIN"
         IF (logger%iter_info%print_level >= debug_print_level .AND. logger%para_env%is_source()) THEN
            WRITE (cp_logger_get_default_unit_nr(logger), '(A,2(1X,ES16.8))') &
               " OT DIIS accepted endpoint converged delta eps:", &
               current_exit_delta, qs_ot_env(1)%kpoint_eps_scf
         END IF
         CALL timestop(handle)
         RETURN
      END IF

      ! put in buffer, dgesv destroys
      qs_ot_env(1)%lss_diis = qs_ot_env(1)%ls_diis

      CALL DGESV(diis_bound + 1, 1, qs_ot_env(1)%lss_diis, diis_m + 1, qs_ot_env(1)%ipivot, &
                 qs_ot_env(1)%c_diis, diis_m + 1, info)

      IF (info /= 0) THEN
         do_ot_sd = .TRUE.
         WRITE (cp_logger_get_default_unit_nr(logger), *) "Singular DIIS matrix"
         ! Keep one physical endpoint and its preconditioned SD direction.  Leaving
         ! diis_iter unchanged retries the same singular Pulay system forever.
         CALL ot_diis_restart_current_anchor(qs_ot_env, j)
      ELSE
         do_ot_sd = .FALSE.
         IF (do_ks) THEN
            DO ispin = 1, nspin
               ! OK, add the vectors now
               CALL dbcsr_set(qs_ot_env(ispin)%matrix_x, 0.0_dp)
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  CALL dbcsr_set(qs_ot_env(ispin)%matrix_x_im, 0.0_dp)
               END IF
               DO i = 1, diis_bound
                  CALL dbcsr_add(qs_ot_env(ispin)%matrix_x, &
                                 qs_ot_env(ispin)%matrix_h_e(i)%matrix, &
                                 alpha_scalar=1.0_dp, beta_scalar=qs_ot_env(1)%c_diis(i))
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_add(qs_ot_env(ispin)%matrix_x_im, &
                                    qs_ot_env(ispin)%matrix_h_e_im(i)%matrix, &
                                    alpha_scalar=1.0_dp, beta_scalar=qs_ot_env(1)%c_diis(i))
                  END IF
               END DO
               DO i = 1, diis_bound
                  CALL dbcsr_add(qs_ot_env(ispin)%matrix_x, &
                                 qs_ot_env(ispin)%matrix_h_x(i)%matrix, &
                                 alpha_scalar=1.0_dp, beta_scalar=qs_ot_env(1)%c_diis(i))
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_add(qs_ot_env(ispin)%matrix_x_im, &
                                    qs_ot_env(ispin)%matrix_h_x_im(i)%matrix, &
                                    alpha_scalar=1.0_dp, beta_scalar=qs_ot_env(1)%c_diis(i))
                  END IF
               END DO
               IF (qs_ot_env(ispin)%settings%do_rotation) THEN
                  CALL dbcsr_set(qs_ot_env(ispin)%rot_mat_x, 0.0_dp)
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_set(qs_ot_env(ispin)%rot_mat_x_im, 0.0_dp)
                  END IF
                  DO i = 1, diis_bound
                     CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_x, &
                                    qs_ot_env(ispin)%rot_mat_h_e(i)%matrix, &
                                    alpha_scalar=1.0_dp, beta_scalar=qs_ot_env(1)%c_diis(i))
                     IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                        CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_x_im, &
                                       qs_ot_env(ispin)%rot_mat_h_e_im(i)%matrix, &
                                       alpha_scalar=1.0_dp, beta_scalar=qs_ot_env(1)%c_diis(i))
                     END IF
                  END DO
                  DO i = 1, diis_bound
                     CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_x, &
                                    qs_ot_env(ispin)%rot_mat_h_x(i)%matrix, &
                                    alpha_scalar=1.0_dp, beta_scalar=qs_ot_env(1)%c_diis(i))
                     IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                        CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_x_im, &
                                       qs_ot_env(ispin)%rot_mat_h_x_im(i)%matrix, &
                                       alpha_scalar=1.0_dp, beta_scalar=qs_ot_env(1)%c_diis(i))
                     END IF
                  END DO
                  IF (kpoint_preconditioned_rotation) THEN
                     norm_weight = ot_kpoint_channel_norm_weight(qs_ot_env(ispin))
                     CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_x, &
                                    qs_ot_env(ispin)%rot_mat_pgx, &
                                    alpha_scalar=1.0_dp, &
                                    beta_scalar=-response_blend*qs_ot_env(1)%ds_min)
                     CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_x, &
                                    qs_ot_env(ispin)%rot_mat_gx, alpha_scalar=1.0_dp, &
                                    beta_scalar=response_blend*qs_ot_env(1)%ds_min*norm_weight)
                     IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                        CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_x_im, &
                                       qs_ot_env(ispin)%rot_mat_pgx_im, &
                                       alpha_scalar=1.0_dp, &
                                       beta_scalar=-response_blend*qs_ot_env(1)%ds_min)
                        CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_x_im, &
                                       qs_ot_env(ispin)%rot_mat_gx_im, alpha_scalar=1.0_dp, &
                                       beta_scalar=response_blend*qs_ot_env(1)%ds_min*norm_weight)
                     END IF
                  END IF
               END IF
            END DO
         END IF
         IF (do_ener) THEN
            DO ispin = 1, nspin
               qs_ot_env(ispin)%ener_x(:) = 0.0_dp
               DO i = 1, diis_bound
                  qs_ot_env(ispin)%ener_x(:) = qs_ot_env(ispin)%ener_x(:) &
                                               + qs_ot_env(1)%c_diis(i)*qs_ot_env(ispin)%ener_h_e(i, :)
               END DO
               DO i = 1, diis_bound
                  qs_ot_env(ispin)%ener_x(:) = qs_ot_env(ispin)%ener_x(:) &
                                               + qs_ot_env(1)%c_diis(i)*qs_ot_env(ispin)%ener_h_x(i, :)
               END DO
            END DO
         END IF
         qs_ot_env(1)%diis_iter = qs_ot_env(1)%diis_iter + 1
         IF (qs_ot_env(1)%settings%safer_diis) THEN
            ! now, final check, is the step in fact in the direction of the -gradient ?
            ! if not we're walking towards a sadle point, and should avoid that
            ! the direction of the step is x_new-x_old
            tr_xold_gx = 0.0_dp
            tr_xnew_gx = 0.0_dp
            IF (do_ks) THEN
               DO ispin = 1, nspin
                  CALL dbcsr_dot(qs_ot_env(ispin)%matrix_h_x(j)%matrix, &
                                 qs_ot_env(ispin)%matrix_gx, tmp)
                  tr_xold_gx = tr_xold_gx + tmp
                  CALL dbcsr_dot(qs_ot_env(ispin)%matrix_x, &
                                 qs_ot_env(ispin)%matrix_gx, tmp)
                  tr_xnew_gx = tr_xnew_gx + tmp
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_dot(qs_ot_env(ispin)%matrix_h_x_im(j)%matrix, &
                                    qs_ot_env(ispin)%matrix_gx_im, tmp)
                     tr_xold_gx = tr_xold_gx + tmp
                     CALL dbcsr_dot(qs_ot_env(ispin)%matrix_x_im, &
                                    qs_ot_env(ispin)%matrix_gx_im, tmp)
                     tr_xnew_gx = tr_xnew_gx + tmp
                  END IF
                  IF (qs_ot_env(ispin)%settings%do_rotation) THEN
                     CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_h_x(j)%matrix, &
                                    qs_ot_env(ispin)%rot_mat_gx, tmp)
                     tr_xold_gx = tr_xold_gx + 0.5_dp*tmp
                     CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_x, &
                                    qs_ot_env(ispin)%rot_mat_gx, tmp)
                     tr_xnew_gx = tr_xnew_gx + 0.5_dp*tmp
                     IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                        CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_h_x_im(j)%matrix, &
                                       qs_ot_env(ispin)%rot_mat_gx_im, tmp)
                        tr_xold_gx = tr_xold_gx + 0.5_dp*tmp
                        CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_x_im, &
                                       qs_ot_env(ispin)%rot_mat_gx_im, tmp)
                        tr_xnew_gx = tr_xnew_gx + 0.5_dp*tmp
                     END IF
                  END IF
               END DO
            END IF
            IF (do_ener) THEN
               DO ispin = 1, nspin
                  tmp = DOT_PRODUCT(qs_ot_env(ispin)%ener_h_x(j, :), qs_ot_env(ispin)%ener_gx(:))
                  tr_xold_gx = tr_xold_gx + tmp
                  tmp = DOT_PRODUCT(qs_ot_env(ispin)%ener_x(:), qs_ot_env(ispin)%ener_gx(:))
                  tr_xnew_gx = tr_xnew_gx + tmp
               END DO
            END IF
            CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, tr_xold_gx)
            CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, tr_xnew_gx)
            overlap = (tr_xnew_gx - tr_xold_gx)
            ! OK, bad luck, take a SD step along the preconditioned gradient
            IF (overlap > 0.0_dp) THEN
               do_ot_sd = .TRUE.
            END IF
         END IF
      END IF

      IF (do_ot_sd) THEN
         qs_ot_env(1)%OT_METHOD_FULL = "OT SD"
         IF (do_ks) THEN
            DO ispin = 1, nspin
               CALL dbcsr_set(qs_ot_env(ispin)%matrix_x, 0.0_dp)
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  CALL dbcsr_set(qs_ot_env(ispin)%matrix_x_im, 0.0_dp)
               END IF
               CALL dbcsr_add(qs_ot_env(ispin)%matrix_x, &
                              qs_ot_env(ispin)%matrix_h_e(j)%matrix, &
                              1.0_dp, 1.0_dp)
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  CALL dbcsr_add(qs_ot_env(ispin)%matrix_x_im, &
                                 qs_ot_env(ispin)%matrix_h_e_im(j)%matrix, &
                                 1.0_dp, 1.0_dp)
               END IF
               CALL dbcsr_add(qs_ot_env(ispin)%matrix_x, &
                              qs_ot_env(ispin)%matrix_h_x(j)%matrix, &
                              1.0_dp, 1.0_dp)
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  CALL dbcsr_add(qs_ot_env(ispin)%matrix_x_im, &
                                 qs_ot_env(ispin)%matrix_h_x_im(j)%matrix, &
                                 1.0_dp, 1.0_dp)
               END IF
               IF (qs_ot_env(ispin)%settings%do_rotation) THEN
                  CALL dbcsr_set(qs_ot_env(ispin)%rot_mat_x, 0.0_dp)
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_set(qs_ot_env(ispin)%rot_mat_x_im, 0.0_dp)
                  END IF
                  CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_x, &
                                 qs_ot_env(ispin)%rot_mat_h_e(j)%matrix, &
                                 1.0_dp, 1.0_dp)
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_x_im, &
                                    qs_ot_env(ispin)%rot_mat_h_e_im(j)%matrix, &
                                    1.0_dp, 1.0_dp)
                  END IF
                  CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_x, &
                                 qs_ot_env(ispin)%rot_mat_h_x(j)%matrix, &
                                 1.0_dp, 1.0_dp)
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_x_im, &
                                    qs_ot_env(ispin)%rot_mat_h_x_im(j)%matrix, &
                                    1.0_dp, 1.0_dp)
                  END IF
                  IF (kpoint_preconditioned_rotation) THEN
                     norm_weight = ot_kpoint_channel_norm_weight(qs_ot_env(ispin))
                     CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_x, &
                                    qs_ot_env(ispin)%rot_mat_pgx, &
                                    alpha_scalar=1.0_dp, &
                                    beta_scalar=-response_blend*qs_ot_env(1)%ds_min)
                     CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_x, &
                                    qs_ot_env(ispin)%rot_mat_gx, alpha_scalar=1.0_dp, &
                                    beta_scalar=response_blend*qs_ot_env(1)%ds_min*norm_weight)
                     IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                        CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_x_im, &
                                       qs_ot_env(ispin)%rot_mat_pgx_im, &
                                       alpha_scalar=1.0_dp, &
                                       beta_scalar=-response_blend*qs_ot_env(1)%ds_min)
                        CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_x_im, &
                                       qs_ot_env(ispin)%rot_mat_gx_im, alpha_scalar=1.0_dp, &
                                       beta_scalar=response_blend*qs_ot_env(1)%ds_min*norm_weight)
                     END IF
                  END IF
               END IF
            END DO
         END IF
         IF (do_ener) THEN
            DO ispin = 1, nspin
               qs_ot_env(ispin)%ener_x(:) = 0._dp
               qs_ot_env(ispin)%ener_x(:) = qs_ot_env(ispin)%ener_x(:) + qs_ot_env(ispin)%ener_h_e(j, :)
               qs_ot_env(ispin)%ener_x(:) = qs_ot_env(ispin)%ener_x(:) + qs_ot_env(ispin)%ener_h_x(j, :)
            END DO
         END IF
      END IF

      DO ispin = 1, nspin
         qs_ot_env(ispin)%diis_last_candidate_sd = do_ot_sd .OR. diis_bound == 1
         qs_ot_env(ispin)%diis_last_candidate_response = &
            kpoint_response_probe .AND. kpoint_preconditioned_rotation
         qs_ot_env(ispin)%diis_last_response_blend = &
            MERGE(response_blend, 0.0_dp, qs_ot_env(ispin)%diis_last_candidate_response)
         qs_ot_env(ispin)%diis_last_response_model_work = &
            MERGE(response_model_work, 0.0_dp, qs_ot_env(ispin)%diis_last_candidate_response)
      END DO

      CALL timestop(handle)

   END SUBROUTINE ot_diis_step

! **************************************************************************************************
!> \brief Energy minimizer by Broyden's method
!> \param qs_ot_env variable to control minimizer behaviour
!> \author Kurt Baarman (09.2010)
! **************************************************************************************************
   SUBROUTINE ot_broyden_step(qs_ot_env)
      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: qs_ot_env

      INTEGER                                            :: diis_bound, diis_m, i, ispin, itmp, j, &
                                                            k, n, nspin
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: circ_index
      LOGICAL                                            :: adaptive_sigma, do_ener, do_ks, &
                                                            enable_flip, forget_history
      REAL(KIND=dp)                                      :: beta, eta, gamma, omega, sigma, &
                                                            sigma_dec, sigma_min, tmp, tmp2
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: f, x
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: G, S
      TYPE(cp_logger_type), POINTER                      :: logger

      eta = qs_ot_env(1)%settings%broyden_eta
      omega = qs_ot_env(1)%settings%broyden_omega
      sigma_dec = qs_ot_env(1)%settings%broyden_sigma_decrease
      sigma_min = qs_ot_env(1)%settings%broyden_sigma_min
      forget_history = qs_ot_env(1)%settings%broyden_forget_history
      adaptive_sigma = qs_ot_env(1)%settings%broyden_adaptive_sigma
      enable_flip = qs_ot_env(1)%settings%broyden_enable_flip

      do_ks = qs_ot_env(1)%settings%ks
      do_ener = qs_ot_env(1)%settings%do_ener

      beta = qs_ot_env(1)%settings%broyden_beta
      gamma = qs_ot_env(1)%settings%broyden_gamma
      IF (adaptive_sigma) THEN
         IF (qs_ot_env(1)%broyden_adaptive_sigma < 0.0_dp) THEN
            sigma = qs_ot_env(1)%settings%broyden_sigma
         ELSE
            sigma = qs_ot_env(1)%broyden_adaptive_sigma
         END IF
      ELSE
         sigma = qs_ot_env(1)%settings%broyden_sigma
      END IF

      ! simplify our life....
      IF (do_ener .OR. .NOT. do_ks .OR. qs_ot_env(1)%settings%do_rotation .OR. &
          ANY(qs_ot_env(:)%has_complex_kpoint_state)) THEN
         CPABORT("Not yet implemented")
      END IF
      !
      nspin = SIZE(qs_ot_env)

      diis_m = qs_ot_env(1)%settings%diis_m

      IF (qs_ot_env(1)%diis_iter < diis_m) THEN
         diis_bound = qs_ot_env(1)%diis_iter + 1
      ELSE
         diis_bound = diis_m
      END IF

      ! We want x:s, f:s and one random vector
      k = 2*diis_bound + 1
      ALLOCATE (S(k, k))
      ALLOCATE (G(k, k))
      ALLOCATE (f(k))
      ALLOCATE (x(k))
      ALLOCATE (circ_index(diis_bound))
      G = 0.0
      DO i = 1, k
         G(i, i) = sigma
      END DO
      S = 0.0

      j = MOD(qs_ot_env(1)%diis_iter, diis_m) + 1 ! index in the circular array

      DO ispin = 1, nspin
         CALL dbcsr_copy(qs_ot_env(ispin)%matrix_h_x(j)%matrix, qs_ot_env(ispin)%matrix_x)
      END DO

      IF (ASSOCIATED(qs_ot_env(1)%preconditioner)) THEN
         qs_ot_env(1)%gnorm = 0.0_dp
         DO ispin = 1, nspin
            CALL apply_preconditioner(qs_ot_env(ispin)%preconditioner, &
                                      qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin)%matrix_h_e(j)%matrix)
            CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin)%matrix_h_e(j)%matrix, &
                           tmp)
            qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + tmp
         END DO
         IF (qs_ot_env(1)%gnorm < 0.0_dp) THEN
            WRITE (cp_logger_get_default_unit_nr(logger), *) "WARNING Preconditioner not positive definite !"
         END IF
         DO ispin = 1, nspin
            CALL dbcsr_scale(qs_ot_env(ispin)%matrix_h_e(j)%matrix, -1.0_dp)
         END DO
      ELSE
         qs_ot_env(1)%gnorm = 0.0_dp
         DO ispin = 1, nspin
            CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin)%matrix_gx, tmp)
            qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + tmp
            CALL dbcsr_add(qs_ot_env(ispin)%matrix_h_e(j)%matrix, &
                           qs_ot_env(ispin)%matrix_gx, alpha_scalar=0.0_dp, beta_scalar=-1.0_dp)
         END DO
      END IF

      k = 0
      n = 0
      CALL dbcsr_get_info(qs_ot_env(1)%matrix_x, nfullrows_total=n)
      DO ispin = 1, nspin
         CALL dbcsr_get_info(qs_ot_env(ispin)%matrix_x, nfullcols_total=itmp)
         k = k + itmp
         IF (qs_ot_env(ispin)%has_complex_kpoint_state) k = k + itmp
      END DO

      ! Handling the case of no free variables to optimize
      IF (INT(n, KIND=int_8)*INT(k, KIND=int_8) /= 0) THEN
         qs_ot_env(1)%delta = SQRT(ABS(qs_ot_env(1)%gnorm)/(INT(n, KIND=int_8)*INT(k, KIND=int_8)))
         qs_ot_env(1)%gradient = -qs_ot_env(1)%gnorm
      ELSE
         qs_ot_env(1)%delta = 0.0_dp
         qs_ot_env(1)%gradient = 0.0_dp
      END IF

      IF (diis_bound == diis_m) THEN
         DO i = 1, diis_bound
            circ_index(i) = MOD(j + i - 1, diis_m) + 1
         END DO
      ELSE
         DO i = 1, diis_bound
            circ_index(i) = i
         END DO
      END IF

      S = 0.0_dp
      DO ispin = 1, nspin
         CALL dbcsr_init_random(qs_ot_env(ispin)%matrix_x)
         DO i = 1, diis_bound
            CALL dbcsr_dot( &
               qs_ot_env(ispin)%matrix_h_x(circ_index(i))%matrix, &
               qs_ot_env(ispin)%matrix_h_x(circ_index(i))%matrix, tmp)
            S(i, i) = S(i, i) + tmp
            CALL dbcsr_dot( &
               qs_ot_env(ispin)%matrix_h_e(circ_index(i))%matrix, &
               qs_ot_env(ispin)%matrix_h_e(circ_index(i))%matrix, tmp)
            S(i + diis_bound, i + diis_bound) = S(i + diis_bound, i + diis_bound) + tmp
            CALL dbcsr_dot( &
               qs_ot_env(ispin)%matrix_h_x(circ_index(i))%matrix, &
               qs_ot_env(ispin)%matrix_x, tmp)
            S(i, 2*diis_bound + 1) = S(i, 2*diis_bound + 1) + tmp
            S(i, 2*diis_bound + 1) = S(2*diis_bound + 1, i)
            CALL dbcsr_dot( &
               qs_ot_env(ispin)%matrix_h_e(circ_index(i))%matrix, &
               qs_ot_env(ispin)%matrix_x, tmp)
            S(i + diis_bound, 2*diis_bound + 1) = S(i + diis_bound, 2*diis_bound + 1) + tmp
            S(i + diis_bound, 2*diis_bound + 1) = S(2*diis_bound + 1, diis_bound + i)
            DO k = (i + 1), diis_bound
               CALL dbcsr_dot( &
                  qs_ot_env(ispin)%matrix_h_x(circ_index(i))%matrix, &
                  qs_ot_env(ispin)%matrix_h_x(circ_index(k))%matrix, &
                  tmp)
               S(i, k) = S(i, k) + tmp
               S(k, i) = S(i, k)
               CALL dbcsr_dot( &
                  qs_ot_env(ispin)%matrix_h_e(circ_index(i))%matrix, &
                  qs_ot_env(ispin)%matrix_h_e(circ_index(k))%matrix, &
                  tmp)
               S(diis_bound + i, diis_bound + k) = S(diis_bound + i, diis_bound + k) + tmp
               S(diis_bound + k, diis_bound + i) = S(diis_bound + i, diis_bound + k)
            END DO
            DO k = 1, diis_bound
               CALL dbcsr_dot( &
                  qs_ot_env(ispin)%matrix_h_x(circ_index(i))%matrix, &
                  qs_ot_env(ispin)%matrix_h_e(circ_index(k))%matrix, tmp)
               S(i, k + diis_bound) = S(i, k + diis_bound) + tmp
               S(k + diis_bound, i) = S(i, k + diis_bound)
            END DO
         END DO
         CALL dbcsr_dot(qs_ot_env(ispin)%matrix_x, qs_ot_env(ispin)%matrix_x, tmp)
         S(2*diis_bound + 1, 2*diis_bound + 1) = S(2*diis_bound + 1, 2*diis_bound + 1) + tmp
      END DO

      ! normalize
      k = 2*diis_bound + 1
      tmp = SQRT(S(k, k))
      S(k, :) = S(k, :)/tmp
      S(:, k) = S(:, k)/tmp

      IF (diis_bound > 1) THEN
         tmp = 0.0_dp
         tmp2 = 0.0_dp
         i = diis_bound
         DO ispin = 1, nspin
            ! dot product of differences
            CALL dbcsr_dot( &
               qs_ot_env(ispin)%matrix_h_x(circ_index(i))%matrix, &
               qs_ot_env(ispin)%matrix_h_e(circ_index(i))%matrix, &
               tmp)
            tmp2 = tmp2 + tmp
            CALL dbcsr_dot( &
               qs_ot_env(ispin)%matrix_h_x(circ_index(i - 1))%matrix, &
               qs_ot_env(ispin)%matrix_h_e(circ_index(i))%matrix, &
               tmp)
            tmp2 = tmp2 - tmp
            CALL dbcsr_dot( &
               qs_ot_env(ispin)%matrix_h_x(circ_index(i))%matrix, &
               qs_ot_env(ispin)%matrix_h_e(circ_index(i - 1))%matrix, &
               tmp)
            tmp2 = tmp2 - tmp
            CALL dbcsr_dot( &
               qs_ot_env(ispin)%matrix_h_x(circ_index(i - 1))%matrix, &
               qs_ot_env(ispin)%matrix_h_e(circ_index(i - 1))%matrix, &
               tmp)
            tmp2 = tmp2 + tmp
         END DO
         qs_ot_env(1)%c_broy(i - 1) = tmp2
      END IF

      qs_ot_env(1)%energy_h(j) = qs_ot_env(1)%etotal

      ! If we went uphill, do backtracking line search
      i = MINLOC(qs_ot_env(1)%energy_h(1:diis_bound), dim=1)
      IF (i /= j) THEN
         sigma = sigma_dec*sigma
         qs_ot_env(1)%OT_METHOD_FULL = "OT BTRK"
         DO ispin = 1, nspin
            CALL dbcsr_set(qs_ot_env(ispin)%matrix_x, 0.0_dp)
            CALL dbcsr_add(qs_ot_env(ispin)%matrix_x, &
                           qs_ot_env(ispin)%matrix_h_x(i)%matrix, &
                           alpha_scalar=1.0_dp, beta_scalar=(1.0_dp - gamma))
            CALL dbcsr_add(qs_ot_env(ispin)%matrix_x, &
                           qs_ot_env(ispin)%matrix_h_x(circ_index(diis_bound))%matrix, &
                           alpha_scalar=1.0_dp, beta_scalar=gamma)
         END DO
      ELSE
         ! Construct G
         DO i = 2, diis_bound
            f = 0.0
            x = 0.0
            ! f is df_i
            x(i) = 1.0
            x(i - 1) = -1.0
            ! x is dx_i
            f(diis_bound + i) = 1.0
            f(diis_bound + i - 1) = -1.0
            tmp = 1.0_dp
            ! We want a pos def Hessian
            IF (enable_flip) THEN
               IF (qs_ot_env(1)%c_broy(i - 1) > 0) THEN
                  !qs_ot_env(1)%OT_METHOD_FULL="OT FLIP"
                  tmp = -1.0_dp
               END IF
            END IF

            ! get dx-Gdf
            x(:) = tmp*x - MATMUL(G, f)
            ! dfSdf
            ! we calculate matmul(S, f) twice. They're small...
            tmp = DOT_PRODUCT(f, MATMUL(S, f))
            ! NOTE THAT S IS SYMMETRIC !!!
            f(:) = MATMUL(S, f)/tmp
            ! the spread is an outer vector product
            G(:, :) = G + SPREAD(x, dim=2, ncopies=SIZE(f))*SPREAD(f, dim=1, ncopies=SIZE(x))
         END DO
         f = 0.0_dp
         f(2*diis_bound) = 1.0_dp
         x(:) = -beta*MATMUL(G, f)

         ! OK, add the vectors now, this sums up to the proposed step
         DO ispin = 1, nspin
            CALL dbcsr_set(qs_ot_env(ispin)%matrix_x, 0.0_dp)
            DO i = 1, diis_bound
               CALL dbcsr_add(qs_ot_env(ispin)%matrix_x, &
                              qs_ot_env(ispin)%matrix_h_e(circ_index(i))%matrix, &
                              alpha_scalar=1.0_dp, beta_scalar=-x(i + diis_bound))
            END DO
            DO i = 1, diis_bound
               CALL dbcsr_add(qs_ot_env(ispin)%matrix_x, &
                              qs_ot_env(ispin)%matrix_h_x(circ_index(i))%matrix, &
                              alpha_scalar=1.0_dp, beta_scalar=x(i))
            END DO
         END DO

         IF (adaptive_sigma) THEN
            tmp = new_sigma(G, S, diis_bound)
            !tmp = tmp * qs_ot_env ( 1 ) % settings % broyden_sigma
            tmp = tmp*eta
            sigma = MIN(omega*sigma, tmp)
         END IF

         ! compute the inner product of direction of the step and gradient
         tmp = 0.0_dp
         DO ispin = 1, nspin
            CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx, &
                           qs_ot_env(ispin)%matrix_x, &
                           tmp2)
            tmp = tmp + tmp2
         END DO

         DO ispin = 1, nspin
            ! if the direction of the step is not in direction of the gradient,
            ! change step sign
            IF (tmp >= 0.0_dp) THEN
               qs_ot_env(1)%OT_METHOD_FULL = "OT TURN"
               IF (forget_history) THEN
                  qs_ot_env(1)%diis_iter = 0
               END IF
               sigma = sigma*sigma_dec
               CALL dbcsr_add(qs_ot_env(ispin)%matrix_x, &
                              qs_ot_env(ispin)%matrix_h_x(circ_index(diis_bound))%matrix, &
                              alpha_scalar=-1.0_dp, beta_scalar=1.0_dp)
            ELSE
               CALL dbcsr_add(qs_ot_env(ispin)%matrix_x, &
                              qs_ot_env(ispin)%matrix_h_x(circ_index(diis_bound))%matrix, &
                              alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
            END IF
         END DO
      END IF

      ! get rid of S, G, f, x, circ_index for next round
      DEALLOCATE (S, G, f, x, circ_index)

      ! update for next round
      qs_ot_env(1)%diis_iter = qs_ot_env(1)%diis_iter + 1
      qs_ot_env(1)%broyden_adaptive_sigma = MAX(sigma, sigma_min)

   END SUBROUTINE ot_broyden_step

! **************************************************************************************************
!> \brief ...
!> \param G ...
!> \param S ...
!> \param n ...
!> \return ...
! **************************************************************************************************
   FUNCTION new_sigma(G, S, n) RESULT(sigma)
!
! Calculate new sigma from eigenvalues of full size G by Arnoldi.
!
! **************************************************************************************************

      REAL(KIND=dp), DIMENSION(:, :)                     :: G, S
      INTEGER                                            :: n
      REAL(KIND=dp)                                      :: sigma

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigv
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: H

      ALLOCATE (H(n, n))
      CALL hess_G(G, S, H, n)
      ALLOCATE (eigv(n))
      CALL diamat_all(H(1:n, 1:n), eigv)

      SELECT CASE (1)
      CASE (1)
         ! This estimator seems to work well. No theory.
         sigma = SUM(ABS(eigv**2))/SUM(ABS(eigv))
      CASE (2)
         ! Estimator based on Frobenius norm minimizer
         sigma = SUM(ABS(eigv))/MAX(1, SIZE(eigv))
      CASE (3)
         ! Estimator based on induced 2-norm
         sigma = (MAXVAL(ABS(eigv)) + MINVAL(ABS(eigv)))*0.5_dp
      END SELECT

      DEALLOCATE (H, eigv)
   END FUNCTION new_sigma

! **************************************************************************************************
!> \brief ...
!> \param G ...
!> \param S ...
!> \param H ...
!> \param n ...
! **************************************************************************************************
   SUBROUTINE hess_G(G, S, H, n)
!
! Make a hessenberg out of G into H. Cf Arnoldi.
! Inner product is weighted by S.
! Possible lucky breakdown at n.
!
! **************************************************************************************************
      REAL(KIND=dp), DIMENSION(:, :)                     :: G, S, H
      INTEGER                                            :: n

      INTEGER                                            :: i, j, k
      REAL(KIND=dp)                                      :: tmp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: v
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: Q

      i = SIZE(G, 1)
      k = SIZE(H, 1)
      ALLOCATE (Q(i, k))
      ALLOCATE (v(i))
      H = 0.0_dp
      Q = 0.0_dp

      Q(:, 1) = 1.0_dp
      tmp = SQRT(DOT_PRODUCT(Q(:, 1), MATMUL(S, Q(:, 1))))
      Q(:, :) = Q(:, :)/tmp

      DO i = 1, k
         v(:) = MATMUL(G, Q(:, i))
         DO j = 1, i
            H(j, i) = DOT_PRODUCT(Q(:, j), MATMUL(S, v))
            v(:) = v - H(j, i)*Q(:, j)
         END DO
         IF (i < k) THEN
            tmp = DOT_PRODUCT(v, MATMUL(S, v))
            IF (tmp <= 0.0_dp) THEN
               n = i
               EXIT
            END IF
            tmp = SQRT(tmp)
            ! Lucky breakdown
            IF (ABS(tmp) < 1e-9_dp) THEN
               n = i
               EXIT
            END IF
            H(i + 1, i) = tmp
            Q(:, i + 1) = v/H(i + 1, i)
         END IF
      END DO

      DEALLOCATE (Q, v)
   END SUBROUTINE hess_G

END MODULE qs_ot_minimizer
