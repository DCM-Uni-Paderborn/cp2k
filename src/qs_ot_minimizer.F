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

   USE cp_dbcsr_api,                    ONLY: dbcsr_add,&
                                              dbcsr_copy,&
                                              dbcsr_get_info,&
                                              dbcsr_p_type,&
                                              dbcsr_scale,&
                                              dbcsr_set
   USE cp_dbcsr_contrib,                ONLY: dbcsr_dot,&
                                              dbcsr_get_diag,&
                                              dbcsr_init_random,&
                                              dbcsr_scale_by_vector
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: high_print_level
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE mathlib,                         ONLY: diamat_all
   USE message_passing,                 ONLY: mp_para_env_type
   USE preconditioner,                  ONLY: apply_preconditioner
   USE qs_ot,                           ONLY: qs_ot_get_derivative,&
                                              qs_ot_get_derivative_ref,&
                                              qs_ot_get_derivative_ref_complex
   USE qs_ot_types,                     ONLY: qs_ot_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   PUBLIC  :: ot_mini

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_ot_minimizer'
   INTEGER, PARAMETER, PRIVATE :: ot_cg_kp_physical_control_min_samples = 8
   INTEGER, PARAMETER, PRIVATE :: ot_cg_kp_mermin_closure_max_backtracks = 8
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
                                        ot_cg_kp_mermin_residual_growth_max_delta = 1.2E-5_dp, &
                                        ot_cg_kp_mermin_recovery_tail_min_delta = 1.5E-6_dp, &
                                        ot_cg_kp_mermin_recovery_tail_max_delta = 5.0E-6_dp, &
                                        ot_cg_kp_mermin_recovery_active_tail_max_delta = 1.0E-4_dp, &
                                        ot_cg_kp_mermin_recovery_tail_min_eps_factor = 1.5_dp, &
                                        ot_cg_kp_mermin_recovery_step_floor = 1.0E-9_dp, &
                                        ot_cg_kp_mermin_recovery_active_step_floor = 1.0E-5_dp, &
                                        ot_cg_kp_physical_raw_delta_floor_ratio = 1.0E-2_dp, &
                                        ot_cg_kp_physical_raw_delta_mermin_floor = 1.0E-8_dp, &
                                        ot_cg_kp_physical_raw_delta_tail_quality = 9.9E-1_dp, &
                                        ot_cg_kp_physical_raw_delta_tail_gap = 1.0E-2_dp, &
                                        ot_cg_kp_physical_control_eps_factor = 3.0_dp, &
                                        ot_cg_kp_physical_control_quality = 8.0E-1_dp, &
                                        ot_cg_kp_physical_control_gap = 5.0E-2_dp, &
                                        ot_cg_kp_residual_trust_min_scale = 1.0E-1_dp, &
                                        ot_cg_kp_residual_trust_quality_floor = 1.0E-4_dp, &
                                        ot_cg_kp_residual_history_damp_quality = 2.0E-1_dp, &
                                        ot_cg_kp_residual_spike_quality_ceiling = 0.0_dp, &
                                        ot_cg_kp_mermin_closure_backtrack = 5.0E-1_dp

CONTAINS
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

      INTEGER                                            :: handle, ispin, nspin
      LOGICAL                                            :: do_ener, do_ks
      REAL(KIND=dp) :: direction_energy_max, direction_matrix_im2, direction_matrix_re2, &
         direction_rotation_im2, direction_rotation_re2, tmp
      TYPE(cp_logger_type), POINTER                      :: logger

      CALL timeset(routineN, handle)

      nspin = SIZE(qs_ot_env)

      do_ks = qs_ot_env(1)%settings%ks
      do_ener = qs_ot_env(1)%settings%do_ener

      CALL recover_collapsed_mermin_search(qs_ot_env)

      IF (qs_ot_env(1)%has_complex_kpoint_state) THEN
         logger => cp_get_default_logger()
         IF (logger%iter_info%print_level >= high_print_level .AND. logger%para_env%is_source()) THEN
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
               IF (ASSOCIATED(qs_ot_env(ispin)%ener_dx)) &
                  direction_energy_max = MAX(direction_energy_max, &
                                             MAXVAL(ABS(qs_ot_env(ispin)%ener_dx)))
            END DO
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

      qs_ot_env(1)%OT_METHOD_FULL = ""

      ! compute the gradient for the variables x
      IF (.NOT. qs_ot_env(1)%energy_only) THEN
         qs_ot_env(1)%gradient = 0.0_dp
         DO ispin = 1, nspin
            IF (do_ks) THEN
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  IF (.NOT. PRESENT(matrix_hc_im)) &
                     CPABORT("Complex k-point OT derivative requires imaginary H(k)*C(k).")
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
                  qs_ot_env(1)%gradient = qs_ot_env(1)%gradient + tmp
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
            qs_ot_env(ispin)%ener_x_step_max = 0.0_dp
            qs_ot_env(ispin)%ener_x_trusted_step_max = 0.0_dp
            qs_ot_env(ispin)%ener_x_trust_scale = 1.0_dp
         END DO
      END IF

      SELECT CASE (qs_ot_env(1)%settings%OT_METHOD)
      CASE ("CG")
         IF (current_point_is_fine(qs_ot_env)) THEN
            qs_ot_env(1)%OT_METHOD_FULL = "OT CG"
            CALL ot_new_cg_direction(qs_ot_env, para_env_inter_kp)
            CALL store_kpoint_mermin_direction_slopes(qs_ot_env, para_env_inter_kp)
            qs_ot_env(1)%line_search_count = 0
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
         CALL ot_diis_step(qs_ot_env)
      CASE ("BROY")
         qs_ot_env(1)%OT_METHOD_FULL = "OT BROY"
         CALL ot_broyden_step(qs_ot_env)
      CASE DEFAULT
         CPABORT("OT_METHOD NYI")
      END SELECT

      IF (do_ener) CALL print_energy_step_diagnostic(qs_ot_env)

      CALL timestop(handle)

   END SUBROUTINE ot_mini

! **************************************************************************************************
!> \brief retain component slopes for an accepted-step Mermin secant
!> \param qs_ot_env ...
!> \param para_env_inter_kp ...
! **************************************************************************************************
   SUBROUTINE store_kpoint_mermin_direction_slopes(qs_ot_env, para_env_inter_kp)
      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: qs_ot_env
      TYPE(mp_para_env_type), OPTIONAL, POINTER          :: para_env_inter_kp

      INTEGER                                            :: ichannel
      REAL(KIND=dp) :: energy_slope, matrix_slope, rotation_slope, slope_sum, tmp
      TYPE(cp_logger_type), POINTER                      :: logger

      DO ichannel = 1, SIZE(qs_ot_env)
         qs_ot_env(ichannel)%mermin_ls_component_slope_valid = .FALSE.
         qs_ot_env(ichannel)%mermin_ls_matrix_slope = 0.0_dp
         qs_ot_env(ichannel)%mermin_ls_rotation_slope = 0.0_dp
         qs_ot_env(ichannel)%mermin_ls_energy_slope = 0.0_dp
      END DO
      logger => cp_get_default_logger()
      IF (logger%iter_info%print_level < high_print_level) RETURN
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
         qs_ot_env(ichannel)%mermin_ls_component_slope_valid = .TRUE.
         matrix_slope = matrix_slope + qs_ot_env(ichannel)%mermin_ls_matrix_slope
         rotation_slope = rotation_slope + qs_ot_env(ichannel)%mermin_ls_rotation_slope
         energy_slope = energy_slope + qs_ot_env(ichannel)%mermin_ls_energy_slope
      END DO

      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, matrix_slope)
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, rotation_slope)
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, energy_slope)
      slope_sum = matrix_slope + rotation_slope + energy_slope
      IF (logger%iter_info%print_level >= high_print_level .AND. logger%para_env%is_source()) &
         WRITE (cp_logger_get_default_unit_nr(logger), '(A,4(1X,ES16.8))') &
            " OT Mermin direction component slopes matrix rotation energy sum:", &
            matrix_slope, rotation_slope, energy_slope, slope_sum

   END SUBROUTINE store_kpoint_mermin_direction_slopes

! **************************************************************************************************
!> \brief compare accepted finite Mermin curvature with endpoint component derivatives
!> \param qs_ot_env ...
!> \param para_env_inter_kp ...
! **************************************************************************************************
   SUBROUTINE audit_kpoint_mermin_accepted_secant(qs_ot_env, para_env_inter_kp)
      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: qs_ot_env
      TYPE(mp_para_env_type), OPTIONAL, POINTER          :: para_env_inter_kp

      INTEGER                                            :: ichannel, valid_i
      LOGICAL                                            :: valid
      REAL(KIND=dp) :: alpha, current_energy, current_matrix, current_rotation, current_total, &
         energy_curvature, initial_energy, initial_matrix, initial_rotation, initial_total, &
         matrix_curvature, matrix_fraction, rotation_curvature, auxiliary_curvature, &
         secant_consistency, total_curvature, tmp
      TYPE(cp_logger_type), POINTER                      :: logger

      IF (.NOT. qs_ot_env(1)%mermin_ls_sample_pending) RETURN
      logger => cp_get_default_logger()
      IF (logger%iter_info%print_level < high_print_level) RETURN
      IF (.NOT. use_three_point_mermin_search(qs_ot_env)) RETURN

      valid = .TRUE.
      initial_matrix = 0.0_dp
      initial_rotation = 0.0_dp
      initial_energy = 0.0_dp
      current_matrix = 0.0_dp
      current_rotation = 0.0_dp
      current_energy = 0.0_dp
      DO ichannel = 1, SIZE(qs_ot_env)
         IF (.NOT. qs_ot_env(ichannel)%has_complex_kpoint_state) CYCLE
         valid = valid .AND. qs_ot_env(ichannel)%mermin_ls_component_slope_valid
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
                             DOT_PRODUCT(qs_ot_env(ichannel)%ener_gx, qs_ot_env(ichannel)%ener_dx)
         END IF
      END DO
      valid_i = MERGE(1, 0, valid)
      CALL ot_mini_min_kpoint_integer(para_env_inter_kp, valid_i)
      valid = valid_i > 0
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, initial_matrix)
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, initial_rotation)
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, initial_energy)
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, current_matrix)
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, current_rotation)
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, current_energy)

      alpha = qs_ot_env(1)%mermin_ls_next_pos
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
      IF (valid) THEN
         matrix_curvature = (current_matrix - initial_matrix)/alpha
         rotation_curvature = (current_rotation - initial_rotation)/alpha
         auxiliary_curvature = (current_energy - initial_energy)/alpha
         total_curvature = (current_total - initial_total)/alpha
         energy_curvature = 2.0_dp* &
                            (qs_ot_env(1)%etotal - qs_ot_env(1)%mermin_ls_reference_energy - &
                             alpha*initial_total)/(alpha*alpha)
         IF (total_curvature > EPSILON(1.0_dp)) &
            matrix_fraction = matrix_curvature/total_curvature
         secant_consistency = ABS(total_curvature - energy_curvature)/ &
                              MAX(ABS(total_curvature), ABS(energy_curvature), EPSILON(1.0_dp))
      END IF

      IF (logger%iter_info%print_level >= high_print_level .AND. logger%para_env%is_source()) THEN
         WRITE (cp_logger_get_default_unit_nr(logger), '(A,1X,L1,8(1X,ES16.8))') &
            " OT Mermin accepted secant slopes valid alpha imat irot iener cmat crot cener current:", &
            valid, alpha, initial_matrix, initial_rotation, initial_energy, &
            current_matrix, current_rotation, current_energy, qs_ot_env(1)%gradient
         WRITE (cp_logger_get_default_unit_nr(logger), '(A,7(1X,ES16.8))') &
            " OT Mermin accepted secant curvature matrix rotation energy gradient Mermin fraction consistency:", &
            matrix_curvature, rotation_curvature, auxiliary_curvature, total_curvature, &
            energy_curvature, matrix_fraction, secant_consistency
      END IF

   END SUBROUTINE audit_kpoint_mermin_accepted_secant

! **************************************************************************************************
!> \brief reduce a k-point integer with MIN when a distributed communicator is present
!> \param para_env_inter_kp ...
!> \param value ...
! **************************************************************************************************
   SUBROUTINE ot_mini_min_kpoint_integer(para_env_inter_kp, value)
      TYPE(mp_para_env_type), OPTIONAL, POINTER          :: para_env_inter_kp
      INTEGER, INTENT(INOUT)                             :: value

      IF (PRESENT(para_env_inter_kp)) THEN
         IF (ASSOCIATED(para_env_inter_kp)) CALL para_env_inter_kp%min(value)
      END IF

   END SUBROUTINE ot_mini_min_kpoint_integer

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
!> \brief RMS metric weight for an irreducible k-point channel
!> \param qs_ot_env ...
!> \return ...
! **************************************************************************************************
   FUNCTION ot_kpoint_channel_norm_weight(qs_ot_env) RESULT(weight)
      TYPE(qs_ot_type), INTENT(IN)                       :: qs_ot_env
      REAL(KIND=dp)                                      :: weight

      weight = 1.0_dp/REAL(MAX(1, qs_ot_env%kpoint_multiplicity), KIND=dp)

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

      REAL(KIND=dp)                                      :: acceptance_tol

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
                qs_ot_env(1)%etotal > qs_ot_env(1)%mermin_ls_reference_energy + acceptance_tol) &
               RETURN
            res = .TRUE.
            RETURN
         END IF

      END IF

   END FUNCTION current_point_is_fine

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
      IF (.NOT. qs_ot_env(1)%settings%occupation_preconditioner) RETURN
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
   FUNCTION kpoint_eps_scaled_tail_min(default_delta, eps_factor, eps_scf) RESULT(delta_min)
      REAL(KIND=dp), INTENT(IN)                          :: default_delta, eps_factor, eps_scf
      REAL(KIND=dp)                                      :: delta_min

      delta_min = default_delta
      IF (eps_scf > 0.0_dp) delta_min = MIN(delta_min, eps_factor*eps_scf)
      delta_min = MAX(delta_min, 1.0E4_dp*EPSILON(1.0_dp))

   END FUNCTION kpoint_eps_scaled_tail_min

! **************************************************************************************************
!> \brief return the fixed-N projected residual for physical k-point tail decisions
!> \param qs_ot_env ...
!> \return ...
! **************************************************************************************************
   FUNCTION kpoint_projected_tail_delta(qs_ot_env) RESULT(delta)
      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: qs_ot_env
      REAL(KIND=dp)                                      :: delta

      delta = qs_ot_env(1)%kpoint_physical_delta
      IF (qs_ot_env(1)%kpoint_projected_delta > 0.0_dp) &
         delta = qs_ot_env(1)%kpoint_projected_delta

   END FUNCTION kpoint_projected_tail_delta

! **************************************************************************************************
!> \brief restart a k-point Mermin line search that has collapsed to a numerically zero step
!> \param qs_ot_env ...
! **************************************************************************************************
   SUBROUTINE recover_collapsed_mermin_search(qs_ot_env)
      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: qs_ot_env

      INTEGER                                            :: count, npos
      LOGICAL                                            :: active_tail, collapsed_sample, &
                                                            collapsed_step
      REAL(KIND=dp)                                      :: mermin_effective_delta, &
                                                            mermin_tail_delta, mermin_tail_min, &
                                                            pos_floor, pos_max, step_floor
      TYPE(cp_logger_type), POINTER                      :: logger

      IF (.NOT. use_three_point_mermin_search(qs_ot_env)) RETURN
      IF (qs_ot_env(1)%line_search_count <= 0) RETURN
      IF (qs_ot_env(1)%mermin_ls_sample_pending .AND. &
          qs_ot_env(1)%line_search_might_be_done) RETURN

      count = qs_ot_env(1)%line_search_count
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
      active_tail = ((mermin_tail_delta > mermin_tail_min .AND. &
                      mermin_tail_delta < ot_cg_kp_mermin_recovery_active_tail_max_delta) .OR. &
                     (mermin_effective_delta > mermin_tail_min .AND. &
                      mermin_effective_delta < ot_cg_kp_mermin_recovery_active_tail_max_delta))
      IF (active_tail .AND. mermin_effective_delta > mermin_tail_min .AND. &
          mermin_effective_delta < ot_cg_kp_mermin_recovery_active_tail_max_delta) &
         mermin_tail_delta = mermin_effective_delta
      IF (active_tail) THEN
         step_floor = MAX(step_floor, ot_cg_kp_mermin_recovery_step_floor* &
                          MAX(1.0_dp, ABS(qs_ot_env(1)%settings%ds_min)))
         pos_floor = MAX(pos_floor, ot_cg_kp_mermin_recovery_active_step_floor* &
                         MAX(1.0_dp, ABS(qs_ot_env(1)%settings%ds_min)))
      END IF
      collapsed_step = ABS(qs_ot_env(1)%ds_min) <= step_floor
      collapsed_sample = active_tail .AND. count >= 3 .AND. pos_max <= pos_floor
      IF (.NOT. (collapsed_step .OR. collapsed_sample)) RETURN

      logger => cp_get_default_logger()
      IF (logger%iter_info%print_level >= high_print_level .AND. logger%para_env%is_source()) THEN
         WRITE (cp_logger_get_default_unit_nr(logger), '(A,1X,L1,2(1X,ES16.8),1X,I5,3(1X,ES16.8))') &
            " OT Mermin LS recovery collapsed active ds floor count pos pmax delta:", &
            active_tail, qs_ot_env(1)%ds_min, step_floor, count, pos_floor, pos_max, mermin_tail_delta
      END IF

      qs_ot_env(1)%ot_pos(:) = 0.0_dp
      qs_ot_env(1)%ot_energy(:) = 0.0_dp
      qs_ot_env(1)%ot_grad(:) = 0.0_dp
      qs_ot_env(1)%line_search_left = -1
      qs_ot_env(1)%line_search_right = -1
      qs_ot_env(1)%line_search_mid = -1
      qs_ot_env(1)%line_search_count = 0
      qs_ot_env(1)%line_search_might_be_done = .FALSE.
      qs_ot_env(1)%energy_only = .FALSE.
      qs_ot_env(1)%ds_min = qs_ot_env(1)%settings%ds_min

   END SUBROUTINE recover_collapsed_mermin_search

! **************************************************************************************************
!> \brief fix the k-point auxiliary-energy gauge and trust scale once per search direction
!> \param qs_ot_env ...
!> \param para_env_inter_kp ...
! **************************************************************************************************
   SUBROUTINE prepare_linear_kpoint_energy_direction(qs_ot_env, para_env_inter_kp)
      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: qs_ot_env
      TYPE(mp_para_env_type), OPTIONAL, POINTER          :: para_env_inter_kp

      INTEGER                                            :: group, i, ispin, ngroups
      LOGICAL                                            :: active, has_rayleigh_target
      REAL(KIND=dp) :: alignment, combined_probe_max, coupled_probe, coupled_rms, probe_pos, &
         rayleigh_correction, scale, shift, target_rms, weight
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: direction_max, radius, rayleigh_diag, &
         rayleigh_shift_sum, rayleigh_target_max, response_sum, target_coupled_dot, &
         target_norm2, coupled_norm2, weighted_direction
      TYPE(cp_logger_type), POINTER                      :: logger

      active = .FALSE.
      DO ispin = 1, SIZE(qs_ot_env)
         qs_ot_env(ispin)%ener_direction_linear = .FALSE.
         active = active .OR. (qs_ot_env(ispin)%has_complex_kpoint_state .AND. &
                               qs_ot_env(ispin)%settings%occupation_preconditioner .AND. &
                               ASSOCIATED(qs_ot_env(ispin)%ener_dx))
      END DO
      IF (.NOT. active) RETURN

      ngroups = MAX(1, MAXVAL(qs_ot_env(:)%occupation_group_index))
      IF (PRESENT(para_env_inter_kp)) THEN
         IF (ASSOCIATED(para_env_inter_kp)) CALL para_env_inter_kp%max(ngroups)
      END IF
      ALLOCATE (coupled_norm2(ngroups), direction_max(ngroups), radius(ngroups), &
                rayleigh_shift_sum(ngroups), rayleigh_target_max(ngroups), &
                response_sum(ngroups), target_coupled_dot(ngroups), target_norm2(ngroups), &
                weighted_direction(ngroups))
      coupled_norm2(:) = 0.0_dp
      direction_max(:) = 0.0_dp
      radius(:) = HUGE(1.0_dp)
      rayleigh_shift_sum(:) = 0.0_dp
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
         IF (ASSOCIATED(qs_ot_env(ispin)%ener_gx_old)) &
            qs_ot_env(ispin)%ener_dx(:) = -qs_ot_env(ispin)%ener_gx_old(:)
         IF (qs_ot_env(ispin)%ener_x_trust_radius > 0.0_dp) &
            radius(group) = MIN(radius(group), qs_ot_env(ispin)%ener_x_trust_radius)
         IF (.NOT. ASSOCIATED(qs_ot_env(ispin)%ener_x_trust_weight)) CYCLE
         IF (SIZE(qs_ot_env(ispin)%ener_x_trust_weight) /= SIZE(qs_ot_env(ispin)%ener_dx)) CYCLE
         has_rayleigh_target = ASSOCIATED(qs_ot_env(ispin)%rot_mat_chc) .AND. &
                                 SIZE(qs_ot_env(ispin)%ener_x) == SIZE(qs_ot_env(ispin)%ener_dx)
         IF (has_rayleigh_target) THEN
            ALLOCATE (rayleigh_diag(SIZE(qs_ot_env(ispin)%ener_dx)))
            CALL dbcsr_get_diag(qs_ot_env(ispin)%rot_mat_chc, rayleigh_diag)
         END IF
         DO i = 1, SIZE(qs_ot_env(ispin)%ener_dx)
            weight = MAX(0.0_dp, qs_ot_env(ispin)%ener_x_trust_weight(i))
            response_sum(group) = response_sum(group) + weight
            weighted_direction(group) = weighted_direction(group) + &
                                        weight*qs_ot_env(ispin)%ener_dx(i)
            IF (has_rayleigh_target) &
               rayleigh_shift_sum(group) = rayleigh_shift_sum(group) + &
                                             weight*(qs_ot_env(ispin)%ener_x(i) - rayleigh_diag(i))
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
         IF (response_sum(group) > EPSILON(response_sum(group))) &
            shift = weighted_direction(group)/response_sum(group)
         qs_ot_env(ispin)%ener_dx(:) = qs_ot_env(ispin)%ener_dx(:) - shift
      END DO

      probe_pos = 0.8_dp*ABS(qs_ot_env(1)%ds_min)
      IF (probe_pos <= EPSILON(probe_pos)) &
         probe_pos = 0.8_dp*ABS(qs_ot_env(1)%settings%ds_min)
      probe_pos = MAX(probe_pos, SQRT(EPSILON(probe_pos)))

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

      DO ispin = 1, SIZE(qs_ot_env)
         IF (.NOT. qs_ot_env(ispin)%has_complex_kpoint_state) CYCLE
         group = MAX(1, MIN(ngroups, qs_ot_env(ispin)%occupation_group_index))
         direction_max(group) = MAX(direction_max(group), &
                                    probe_pos*MAXVAL(ABS(qs_ot_env(ispin)%ener_dx)))
      END DO
      IF (PRESENT(para_env_inter_kp)) THEN
         IF (ASSOCIATED(para_env_inter_kp)) CALL para_env_inter_kp%max(direction_max)
      END IF

      DO ispin = 1, SIZE(qs_ot_env)
         IF (.NOT. qs_ot_env(ispin)%has_complex_kpoint_state) CYCLE
         group = MAX(1, MIN(ngroups, qs_ot_env(ispin)%occupation_group_index))
         scale = 1.0_dp
         IF (radius(group) < HUGE(1.0_dp) .AND. direction_max(group) > radius(group)) &
            scale = radius(group)/direction_max(group)
         qs_ot_env(ispin)%ener_dx(:) = scale*qs_ot_env(ispin)%ener_dx(:)
         qs_ot_env(ispin)%ener_x_trust_scale = MIN(qs_ot_env(ispin)%ener_x_trust_scale, scale)
         qs_ot_env(ispin)%ener_direction_linear = .TRUE.
      END DO

      logger => cp_get_default_logger()
      IF (logger%iter_info%print_level >= high_print_level .AND. logger%para_env%is_source()) THEN
         DO group = 1, ngroups
            combined_probe_max = direction_max(group)
            target_rms = SQRT(MAX(0.0_dp, target_norm2(group))/ &
                              MAX(response_sum(group), EPSILON(response_sum(group))))
            coupled_rms = SQRT(MAX(0.0_dp, coupled_norm2(group))/ &
                               MAX(response_sum(group), EPSILON(response_sum(group))))
            alignment = 0.0_dp
            IF (target_norm2(group) > EPSILON(target_norm2(group)) .AND. &
                coupled_norm2(group) > EPSILON(coupled_norm2(group))) &
               alignment = target_coupled_dot(group)/ &
                           SQRT(target_norm2(group)*coupled_norm2(group))
            WRITE (cp_logger_get_default_unit_nr(logger), '(A,1X,I3,8(1X,ES16.8))') &
               " OT fixed-N Rayleigh response audit group tmax trms srms align smax radius probe response:", &
               group, rayleigh_target_max(group), target_rms, coupled_rms, alignment, &
               combined_probe_max, radius(group), probe_pos, response_sum(group)
         END DO
      END IF

      DEALLOCATE (coupled_norm2, direction_max, radius, rayleigh_shift_sum, rayleigh_target_max, &
                  response_sum, target_coupled_dot, target_norm2, weighted_direction)

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
      REAL(KIND=dp)                                      :: max_projected_step, raw_step_max, &
                                                            response_weight, step_scale, &
                                                            step_shift, trusted_step_max
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: limited_step, raw_step

      n = SIZE(direction)
      IF (n == 0) RETURN

      IF (qs_ot_env%ener_direction_linear) THEN
         raw_step_max = MAXVAL(ABS(alpha*direction))
         qs_ot_env%ener_x(:) = qs_ot_env%ener_x(:) + alpha*direction(:)
         qs_ot_env%ener_x_step_max = raw_step_max
         qs_ot_env%ener_x_trusted_step_max = raw_step_max
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
         IF (ASSOCIATED(qs_ot_env%ener_x_trust_weight)) &
            has_response_trust = SIZE(qs_ot_env%ener_x_trust_weight) == n

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
                  IF (qs_ot_env%ener_x_trust_weight(i) > 0.0_dp) &
                     max_projected_step = MAX(max_projected_step, ABS(raw_step(i) - step_shift))
               END DO
               IF (max_projected_step > qs_ot_env%ener_x_trust_radius) &
                  step_scale = qs_ot_env%ener_x_trust_radius/max_projected_step
               DO i = 1, n
                  IF (qs_ot_env%ener_x_trust_weight(i) > 0.0_dp) &
                     limited_step(i) = step_scale*(raw_step(i) - step_shift)
               END DO
            END IF
         END IF

         DO i = 1, n
            IF (ABS(limited_step(i)) > qs_ot_env%ener_x_trust_radius) &
               limited_step(i) = SIGN(qs_ot_env%ener_x_trust_radius, limited_step(i))
         END DO
         trusted_step_max = MAXVAL(ABS(limited_step))
         IF (raw_step_max > EPSILON(raw_step_max)) step_scale = trusted_step_max/raw_step_max
      END IF
      qs_ot_env%ener_x = qs_ot_env%ener_x + limited_step
      qs_ot_env%ener_x_step_max = raw_step_max
      qs_ot_env%ener_x_trusted_step_max = trusted_step_max
      qs_ot_env%ener_x_trust_scale = MIN(qs_ot_env%ener_x_trust_scale, step_scale)

      DEALLOCATE (raw_step, limited_step)

   END SUBROUTINE take_energy_step

! **************************************************************************************************
!> \brief high-print diagnostic for auxiliary-energy trust-region steps
!> \param qs_ot_env ...
! **************************************************************************************************
   SUBROUTINE print_energy_step_diagnostic(qs_ot_env)
      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: qs_ot_env

      INTEGER                                            :: ispin
      REAL(KIND=dp) :: entropy, frac_weight, min_abs_eps_mu, raw_step_max, rayleigh_edge_gap, &
         rayleigh_gap, scale_min, tether_energy, tether_grad_max, tether_shift_abs, &
         tether_strength, tether_weight_sum, trusted_step_max
      TYPE(cp_logger_type), POINTER                      :: logger

      logger => cp_get_default_logger()
      IF (logger%iter_info%print_level < high_print_level) RETURN

      raw_step_max = 0.0_dp
      trusted_step_max = 0.0_dp
      scale_min = 1.0_dp
      entropy = 0.0_dp
      frac_weight = 0.0_dp
      min_abs_eps_mu = HUGE(1.0_dp)
      rayleigh_gap = 0.0_dp
      rayleigh_edge_gap = 0.0_dp
      tether_energy = 0.0_dp
      tether_strength = 0.0_dp
      tether_weight_sum = 0.0_dp
      tether_shift_abs = 0.0_dp
      tether_grad_max = 0.0_dp
      DO ispin = 1, SIZE(qs_ot_env)
         IF (qs_ot_env(ispin)%ener_x_trust_radius <= 0.0_dp) CYCLE
         raw_step_max = MAX(raw_step_max, qs_ot_env(ispin)%ener_x_step_max)
         trusted_step_max = MAX(trusted_step_max, qs_ot_env(ispin)%ener_x_trusted_step_max)
         scale_min = MIN(scale_min, qs_ot_env(ispin)%ener_x_trust_scale)
         entropy = entropy + qs_ot_env(ispin)%ener_x_trust_entropy
         frac_weight = frac_weight + qs_ot_env(ispin)%ener_x_trust_frac_weight
         min_abs_eps_mu = MIN(min_abs_eps_mu, qs_ot_env(ispin)%ener_x_trust_min_abs_eps_mu)
         rayleigh_gap = MAX(rayleigh_gap, qs_ot_env(ispin)%ener_x_trust_rayleigh_gap)
         rayleigh_edge_gap = MAX(rayleigh_edge_gap, qs_ot_env(ispin)%ener_x_trust_rayleigh_edge_gap)
         tether_energy = tether_energy + qs_ot_env(ispin)%canonical_tether_energy
         tether_strength = MAX(tether_strength, qs_ot_env(ispin)%canonical_tether_strength)
         tether_weight_sum = tether_weight_sum + qs_ot_env(ispin)%canonical_tether_weight_sum
         tether_shift_abs = MAX(tether_shift_abs, ABS(qs_ot_env(ispin)%canonical_tether_shift))
         tether_grad_max = MAX(tether_grad_max, qs_ot_env(ispin)%canonical_tether_grad_max)
      END DO

      IF (raw_step_max > 0.0_dp .AND. logger%para_env%is_source()) THEN
         WRITE (cp_logger_get_default_unit_nr(logger), '(A,9(1X,ES16.8),1X,I5)') &
            " OT ener_x step raw kept scale ds kTS frac min_abs ray_gap edge_gap ls:", &
            raw_step_max, trusted_step_max, scale_min, qs_ot_env(1)%ds_min, entropy, &
            frac_weight, min_abs_eps_mu, rayleigh_gap, rayleigh_edge_gap, &
            qs_ot_env(1)%line_search_count
         IF (tether_weight_sum > 0.0_dp) THEN
            WRITE (cp_logger_get_default_unit_nr(logger), '(A,5(1X,ES16.8))') &
               " OT canonical tether energy strength weight shift grad:", &
               tether_energy, tether_strength, tether_weight_sum, tether_shift_abs, tether_grad_max
         END IF
      END IF

   END SUBROUTINE print_energy_step_diagnostic

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
      LOGICAL :: mermin_boundary_extend, mermin_closure_reset, mermin_search, &
         mermin_uphill_closure
      REAL(KIND=dp) :: best_drop, best_energy, closure_pos_floor, denom, ds, fa, fb, fc, &
         model_drop, model_quality, nom, pos, slope_mismatch, slope_model, tol, val, xa, xb, xc
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
            IF (logger%iter_info%print_level >= high_print_level .AND. logger%para_env%is_source()) THEN
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
         mermin_closure_reset = .FALSE.
         IF (qs_ot_env(1)%mermin_ls_sample_pending .AND. .NOT. mermin_uphill_closure) THEN
            best_index = count
            best_energy = qs_ot_env(1)%OT_energy(count)
            best_drop = qs_ot_env(1)%mermin_ls_reference_energy - best_energy
            qs_ot_env(1)%OT_pos(count + 1) = qs_ot_env(1)%OT_pos(count)
            qs_ot_env(1)%mermin_ls_best_sample_index = best_index
            qs_ot_env(1)%mermin_ls_best_sample_energy = best_energy
            qs_ot_env(1)%mermin_ls_best_sample_drop = MAX(0.0_dp, best_drop)
            IF (qs_ot_env(1)%mermin_ls_model_drop > EPSILON(1.0_dp)) &
               qs_ot_env(1)%mermin_ls_model_quality = &
               MAX(0.0_dp, best_drop)/qs_ot_env(1)%mermin_ls_model_drop
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
            IF (logger%iter_info%print_level >= high_print_level .AND. logger%para_env%is_source()) THEN
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
               IF (qs_ot_env(1)%OT_energy(i) < qs_ot_env(1)%OT_energy(best_index)) &
                  best_index = i
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
         qs_ot_env(1)%energy_only = mermin_uphill_closure .AND. .NOT. mermin_closure_reset
         qs_ot_env(1)%line_search_might_be_done = .NOT. mermin_uphill_closure .OR. mermin_closure_reset
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
               CALL apply_preconditioner(qs_ot_env(ispin)%preconditioner, &
                                         qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin)%matrix_dx)
               CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin)%matrix_dx, tmp)
               qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + tmp
               gnorm_metric = gnorm_metric + norm_weight*tmp
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  CALL apply_preconditioner(qs_ot_env(ispin)%preconditioner, &
                                            qs_ot_env(ispin)%matrix_gx_im, qs_ot_env(ispin)%matrix_dx_im)
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
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) &
                  CALL dbcsr_scale(qs_ot_env(ispin)%matrix_dx_im, -1.0_dp)
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
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) &
                     CALL dbcsr_scale(qs_ot_env(ispin)%rot_mat_dx_im, -1.0_dp)
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
! **************************************************************************************************
   SUBROUTINE ot_new_cg_direction(qs_ot_env, para_env_inter_kp)
      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: qs_ot_env
      TYPE(mp_para_env_type), OPTIONAL, POINTER          :: para_env_inter_kp

      CHARACTER(len=*), PARAMETER :: routineN = 'ot_new_cg_direction'

      INTEGER                                            :: handle, ispin, nspin, schur_attempt, &
                                                            schur_flavor, schur_last_attempt, &
                                                            schur_last_flavor
      LOGICAL :: complex_kpoint_history, do_ener, do_ks, history_spike_reset, &
         kpoint_preconditioned_direction, &
         mermin_residual_growth_reset, negative_beta, physical_raw_gap_guard, &
         projected_raw_tail_release, reset_to_sd, residual_history_reset, &
         residual_spike_history_reset, residual_tail_spike_reset, schur_history_change, &
         schur_history_mixed, schur_projector_jump_reset, schur_projector_tail_reset
      REAL(KIND=dp) :: beta_pr, beta_raw, gnorm_cross, gnorm_ener, gnorm_matrix_im, &
         gnorm_matrix_re, gnorm_metric, gnorm_ratio, gnorm_rot_im, gnorm_rot_re, &
         matrix_gnorm_metric, mermin_tail_min, norm_weight, nvar, physical_control_delta, &
         physical_delta, physical_gnorm_metric, physical_tail_delta, raw_delta, raw_guard_delta, &
         raw_matrix_gnorm_im, raw_matrix_gnorm_metric_im, raw_matrix_gnorm_metric_re, &
         raw_matrix_gnorm_re, residual_growth_delta, residual_spike_quality, residual_tail_min, &
         residual_trust_quality, residual_trust_scale, sample_accepted_drop, sample_gap, &
         sample_guard_drop, sample_quality
      REAL(KIND=dp) :: sample_residual_growth_ratio, test_down, tmp
      TYPE(cp_logger_type), POINTER                      :: logger

      CALL timeset(routineN, handle)

      nspin = SIZE(qs_ot_env)
      logger => cp_get_default_logger()

      do_ks = qs_ot_env(1)%settings%ks
      do_ener = qs_ot_env(1)%settings%do_ener
      complex_kpoint_history = .FALSE.
      kpoint_preconditioned_direction = .FALSE.
      gnorm_metric = 0.0_dp
      DO ispin = 1, nspin
         complex_kpoint_history = complex_kpoint_history .OR. qs_ot_env(ispin)%has_complex_kpoint_state
      END DO
      IF (complex_kpoint_history) THEN
         kpoint_preconditioned_direction = .TRUE.
         DO ispin = 1, nspin
            IF (qs_ot_env(ispin)%has_complex_kpoint_state .AND. &
                .NOT. qs_ot_env(ispin)%kpoint_preconditioned_direction_ready) &
               kpoint_preconditioned_direction = .FALSE.
         END DO
      END IF
      history_spike_reset = .FALSE.
      negative_beta = .FALSE.
      mermin_residual_growth_reset = .FALSE.
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
      ! sub-flavor changes are diagnosed here, but only explicitly mixed
      ! sub-flavor windows suppress accepted-step quality calibration.
      schur_history_change = qs_ot_env(1)%settings%occupation_preconditioner .AND. &
                             qs_ot_env(1)%settings%do_rotation .AND. &
                             (((schur_attempt /= schur_last_attempt) .AND. &
                               schur_attempt /= 0 .AND. schur_last_attempt /= 0 .AND. &
                               (schur_attempt <= 2 .OR. schur_last_attempt <= 2)) .OR. &
                              (schur_attempt == 3 .AND. schur_last_attempt == 3 .AND. &
                               schur_flavor /= schur_last_flavor .AND. &
                               schur_flavor /= 0 .AND. schur_last_flavor /= 0 .AND. &
                               (schur_flavor < 0 .OR. schur_last_flavor < 0 .OR. &
                                schur_flavor == 42 .OR. schur_last_flavor == 42)))
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

      gnorm_cross = 0.0_dp
      IF (do_ks) THEN
         DO ispin = 1, nspin
            norm_weight = ot_kpoint_channel_norm_weight(qs_ot_env(ispin))
            CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin)%matrix_gx_old, tmp)
            gnorm_cross = gnorm_cross + tmp
            CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin)%matrix_gx, tmp)
            raw_matrix_gnorm_re = raw_matrix_gnorm_re + tmp
            raw_matrix_gnorm_metric_re = raw_matrix_gnorm_metric_re + norm_weight*tmp
            IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
               CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx_im, qs_ot_env(ispin)%matrix_gx_old_im, tmp)
               gnorm_cross = gnorm_cross + tmp
               CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx_im, qs_ot_env(ispin)%matrix_gx_im, tmp)
               raw_matrix_gnorm_im = raw_matrix_gnorm_im + tmp
               raw_matrix_gnorm_metric_im = raw_matrix_gnorm_metric_im + norm_weight*tmp
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
            CALL apply_preconditioner(qs_ot_env(ispin)%preconditioner, &
                                      qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin)%matrix_gx_old)
            IF (qs_ot_env(ispin)%has_complex_kpoint_state .AND. &
                ASSOCIATED(qs_ot_env(ispin)%matrix_occ_precond_scale)) &
               CALL dbcsr_scale_by_vector(qs_ot_env(ispin)%matrix_gx_old, &
                                          alpha=qs_ot_env(ispin)%matrix_occ_precond_scale, &
                                          side='right')
            IF (kpoint_preconditioned_direction) &
               CALL dbcsr_scale(qs_ot_env(ispin)%matrix_gx_old, norm_weight)
            IF (qs_ot_env(ispin)%has_complex_kpoint_state) &
               CALL apply_preconditioner(qs_ot_env(ispin)%preconditioner, &
                                         qs_ot_env(ispin)%matrix_gx_im, qs_ot_env(ispin)%matrix_gx_old_im)
            IF (qs_ot_env(ispin)%has_complex_kpoint_state .AND. &
                ASSOCIATED(qs_ot_env(ispin)%matrix_occ_precond_scale)) &
               CALL dbcsr_scale_by_vector(qs_ot_env(ispin)%matrix_gx_old_im, &
                                          alpha=qs_ot_env(ispin)%matrix_occ_precond_scale, &
                                          side='right')
            IF (kpoint_preconditioned_direction .AND. qs_ot_env(ispin)%has_complex_kpoint_state) &
               CALL dbcsr_scale(qs_ot_env(ispin)%matrix_gx_old_im, norm_weight)
         END DO
         qs_ot_env(1)%gnorm = 0.0_dp
         IF (do_ks) THEN
            DO ispin = 1, nspin
               norm_weight = ot_kpoint_channel_norm_weight(qs_ot_env(ispin))
               CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin)%matrix_gx_old, tmp)
               qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + tmp
               gnorm_matrix_re = gnorm_matrix_re + tmp
               IF (kpoint_preconditioned_direction) THEN
                  gnorm_metric = gnorm_metric + tmp
               ELSE
                  gnorm_metric = gnorm_metric + norm_weight*tmp
               END IF
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx_im, qs_ot_env(ispin)%matrix_gx_old_im, tmp)
                  qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + tmp
                  gnorm_matrix_im = gnorm_matrix_im + tmp
                  IF (kpoint_preconditioned_direction) THEN
                     gnorm_metric = gnorm_metric + tmp
                  ELSE
                     gnorm_metric = gnorm_metric + norm_weight*tmp
                  END IF
               END IF
            END DO
            IF (qs_ot_env(1)%gnorm < 0.0_dp) THEN
               WRITE (cp_logger_get_default_unit_nr(logger), *) "WARNING Preconditioner not positive definite !"
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
                     CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx, qs_ot_env(ispin)%rot_mat_gx_old, tmp)
                  END IF
                  ! added 0.5, because we have (antisymmetry) only half the number of variables
                  qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + 0.5_dp*tmp
                  gnorm_rot_re = gnorm_rot_re + 0.5_dp*tmp
                  IF (kpoint_preconditioned_direction) THEN
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
                        CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx_im, qs_ot_env(ispin)%rot_mat_gx_old_im, tmp)
                     END IF
                     qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + 0.5_dp*tmp
                     gnorm_rot_im = gnorm_rot_im + 0.5_dp*tmp
                     IF (kpoint_preconditioned_direction) THEN
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
                  qs_ot_env(ispin)%ener_gx_old = qs_ot_env(ispin)%ener_gx
                  tmp = DOT_PRODUCT(qs_ot_env(ispin)%ener_gx, qs_ot_env(ispin)%ener_gx_old)
               END IF
               qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + tmp
               gnorm_ener = gnorm_ener + tmp
               IF (kpoint_preconditioned_direction) THEN
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
               CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin)%matrix_gx, tmp)
               qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + tmp
               gnorm_matrix_re = gnorm_matrix_re + tmp
               gnorm_metric = gnorm_metric + norm_weight*tmp
               CALL dbcsr_copy(qs_ot_env(ispin)%matrix_gx_old, qs_ot_env(ispin)%matrix_gx)
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx_im, qs_ot_env(ispin)%matrix_gx_im, tmp)
                  qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + tmp
                  gnorm_matrix_im = gnorm_matrix_im + tmp
                  gnorm_metric = gnorm_metric + norm_weight*tmp
                  CALL dbcsr_copy(qs_ot_env(ispin)%matrix_gx_old_im, qs_ot_env(ispin)%matrix_gx_im)
               END IF
            END DO
            IF (qs_ot_env(1)%settings%do_rotation) THEN
               DO ispin = 1, nspin
                  norm_weight = ot_kpoint_channel_norm_weight(qs_ot_env(ispin))
                  CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx, qs_ot_env(ispin)%rot_mat_gx, tmp)
                  ! added 0.5, because we have (antisymmetry) only half the number of variables
                  qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + 0.5_dp*tmp
                  gnorm_rot_re = gnorm_rot_re + 0.5_dp*tmp
                  gnorm_metric = gnorm_metric + norm_weight*0.5_dp*tmp
                  CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_gx_old, qs_ot_env(ispin)%rot_mat_gx)
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx_im, qs_ot_env(ispin)%rot_mat_gx_im, tmp)
                     qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + 0.5_dp*tmp
                     gnorm_rot_im = gnorm_rot_im + 0.5_dp*tmp
                     gnorm_metric = gnorm_metric + norm_weight*0.5_dp*tmp
                     CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_gx_old_im, qs_ot_env(ispin)%rot_mat_gx_im)
                  END IF
               END DO
            END IF
         END IF
         IF (do_ener) THEN
            DO ispin = 1, nspin
               norm_weight = ot_kpoint_channel_norm_weight(qs_ot_env(ispin))
               tmp = DOT_PRODUCT(qs_ot_env(ispin)%ener_gx, qs_ot_env(ispin)%ener_gx)
               qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + tmp
               gnorm_ener = gnorm_ener + tmp
               gnorm_metric = gnorm_metric + norm_weight*tmp
               qs_ot_env(ispin)%ener_gx_old = qs_ot_env(ispin)%ener_gx
            END DO
         END IF
      END IF

      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, qs_ot_env(1)%gnorm)
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, gnorm_metric)
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, gnorm_matrix_re)
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, gnorm_matrix_im)
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, gnorm_rot_re)
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, gnorm_rot_im)
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, gnorm_ener)
      nvar = ot_kpoint_channel_nvar(qs_ot_env, do_ks, do_ener)
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, nvar)
      ! Handling the case of no free variables to optimize
      IF (nvar > 0.0_dp) THEN
         qs_ot_env(1)%delta = SQRT(ABS(gnorm_metric)/nvar)
         raw_delta = qs_ot_env(1)%delta
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
            IF (qs_ot_env(1)%kpoint_physical_delta_ready) THEN
               IF (physical_delta < qs_ot_env(1)%delta) qs_ot_env(1)%delta = physical_delta
            ELSE IF (qs_ot_env(1)%kpoint_eps_scf > 0.0_dp .AND. &
                     physical_delta > 0.0_dp .AND. physical_delta < qs_ot_env(1)%delta .AND. &
                     qs_ot_env(1)%mermin_ls_accepted_quality_samples >= &
                     ot_cg_kp_physical_control_min_samples .AND. &
                     qs_ot_env(1)%mermin_ls_accepted_quality_ema >= &
                     ot_cg_kp_physical_control_quality .AND. &
                     qs_ot_env(1)%mermin_ls_accepted_gap_ema <= &
                     ot_cg_kp_physical_control_gap) THEN
               ! Keep non-ready tail steps above EPS_SCF, but steer trust/projector logic
               ! with the measured fixed-N residual instead of a preconditioned raw spike.
               physical_control_delta = MAX(physical_delta, &
                                            ot_cg_kp_physical_control_eps_factor* &
                                            qs_ot_env(1)%kpoint_eps_scf)
               IF (physical_control_delta < qs_ot_env(1)%delta) qs_ot_env(1)%delta = physical_control_delta
            END IF
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
         IF (qs_ot_env(1)%mermin_ls_best_sample_drop > EPSILON(1.0_dp)) &
            sample_quality = sample_accepted_drop/qs_ot_env(1)%mermin_ls_best_sample_drop
         sample_guard_drop = sample_accepted_drop
         IF (qs_ot_env(1)%mermin_ls_model_quality < 0.5_dp) &
            sample_guard_drop = MAX(sample_guard_drop, qs_ot_env(1)%mermin_ls_model_drop)
         physical_raw_gap_guard = complex_kpoint_history .AND. &
                                  qs_ot_env(1)%settings%occupation_preconditioner .AND. &
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
            IF (projected_raw_tail_release) &
               raw_guard_delta = MIN(raw_guard_delta, qs_ot_env(1)%kpoint_projected_delta)
            IF (raw_guard_delta > qs_ot_env(1)%delta) THEN
               qs_ot_env(1)%delta = raw_guard_delta
               qs_ot_env(1)%kpoint_physical_delta = MAX(qs_ot_env(1)%kpoint_physical_delta, raw_guard_delta)
               IF (logger%iter_info%print_level >= high_print_level .AND. logger%para_env%is_source()) THEN
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
         residual_growth_delta = qs_ot_env(1)%delta
         IF (physical_tail_delta > 0.0_dp) residual_growth_delta = physical_tail_delta
         sample_residual_growth_ratio = 0.0_dp
         IF (qs_ot_env(1)%schur_model_delta > EPSILON(1.0_dp)) &
            sample_residual_growth_ratio = residual_growth_delta/qs_ot_env(1)%schur_model_delta
         mermin_residual_growth_reset = &
            qs_ot_env(1)%settings%occupation_preconditioner .AND. &
            qs_ot_env(1)%settings%do_rotation .AND. do_ener .AND. &
            sample_accepted_drop > ot_cg_kp_mermin_recovery_step_floor .AND. &
            qs_ot_env(1)%schur_model_delta >= ot_cg_kp_mermin_recovery_tail_min_delta .AND. &
            qs_ot_env(1)%schur_model_delta <= ot_cg_kp_mermin_residual_growth_max_delta .AND. &
            residual_growth_delta > ot_cg_kp_mermin_residual_growth_ratio* &
            qs_ot_env(1)%schur_model_delta
         IF (mermin_residual_growth_reset) THEN
            qs_ot_env(1)%schur_residual_history_reset_pending = .TRUE.
            qs_ot_env(1)%schur_residual_history_quality = 0.0_dp
            IF (logger%iter_info%print_level >= high_print_level .AND. logger%para_env%is_source()) THEN
               WRITE (cp_logger_get_default_unit_nr(logger), '(A,5(1X,ES16.8))') &
                  " OT Mermin LS residual growth reset prev next ratio adrop bdrop:", &
                  qs_ot_env(1)%schur_model_delta, residual_growth_delta, sample_residual_growth_ratio, &
                  sample_accepted_drop, qs_ot_env(1)%mermin_ls_best_sample_drop
            END IF
         END IF
         IF (logger%iter_info%print_level >= high_print_level .AND. logger%para_env%is_source()) THEN
            WRITE (cp_logger_get_default_unit_nr(logger), '(A,1X,I5,12(1X,ES16.8))') &
               " OT Mermin LS accepted-state audit idx ref best accepted bdrop adrop gap q mdrop mq bpos next delta:", &
               qs_ot_env(1)%mermin_ls_best_sample_index, qs_ot_env(1)%mermin_ls_reference_energy, &
               qs_ot_env(1)%mermin_ls_best_sample_energy, qs_ot_env(1)%etotal, &
               qs_ot_env(1)%mermin_ls_best_sample_drop, sample_accepted_drop, sample_gap, &
               sample_quality, qs_ot_env(1)%mermin_ls_model_drop, qs_ot_env(1)%mermin_ls_model_quality, &
               qs_ot_env(1)%mermin_ls_best_sample_pos, qs_ot_env(1)%mermin_ls_next_pos, qs_ot_env(1)%delta
         END IF
         qs_ot_env(1)%mermin_ls_sample_pending = .FALSE.
      END IF
      IF (logger%iter_info%print_level >= high_print_level .AND. logger%para_env%is_source()) THEN
         WRITE (cp_logger_get_default_unit_nr(logger), '(A,6(1X,ES16.8))') &
            " OT CG channel norms C_re C_im R_re R_im eps total:", &
            gnorm_matrix_re, gnorm_matrix_im, gnorm_rot_re, gnorm_rot_im, &
            gnorm_ener, qs_ot_env(1)%gnorm
         IF (qs_ot_env(1)%kpoint_physical_delta_groups > 0) &
            WRITE (cp_logger_get_default_unit_nr(logger), '(A,1X,L1,1X,I6,8(1X,ES16.8))') &
            " OT CG physical residual ready groups raw effective guarded projected matrix phys_rot phys_eps gauge:", &
            qs_ot_env(1)%kpoint_physical_delta_ready, qs_ot_env(1)%kpoint_physical_delta_groups, &
            qs_ot_env(1)%kpoint_raw_delta, qs_ot_env(1)%delta, qs_ot_env(1)%kpoint_physical_delta, &
            qs_ot_env(1)%kpoint_projected_delta, &
            MAX(0.0_dp, raw_matrix_gnorm_re) + MAX(0.0_dp, raw_matrix_gnorm_im), &
            qs_ot_env(1)%kpoint_physical_rot_gnorm, &
            qs_ot_env(1)%kpoint_physical_energy_gnorm, qs_ot_env(1)%kpoint_gauge_gnorm
         IF (qs_ot_env(1)%kpoint_physical_delta_groups > 0) &
            WRITE (cp_logger_get_default_unit_nr(logger), '(A,4(1X,ES16.8))') &
            " OT CG matrix residual raw_re raw_im precond_re precond_im:", &
            raw_matrix_gnorm_re, raw_matrix_gnorm_im, gnorm_matrix_re, gnorm_matrix_im
      END IF
      beta_raw = beta_pr
      negative_beta = beta_pr < 0.0_dp
      beta_pr = MAX(beta_pr, 0.0_dp) ! reset to SD
      gnorm_ratio = qs_ot_env(1)%gnorm/MAX(qs_ot_env(1)%gnorm_old, EPSILON(1.0_dp))
      physical_tail_delta = kpoint_projected_tail_delta(qs_ot_env)
      ! Projector-like Schur flavors can validly alternate in the metallic tail.
      ! Do not carry a huge CG beta across that local response-model change.
      schur_projector_jump_reset = complex_kpoint_history .AND. &
                                   qs_ot_env(1)%settings%occupation_preconditioner .AND. &
                                   qs_ot_env(1)%settings%do_rotation .AND. do_ener .AND. &
                                   beta_pr > ot_cg_kp_history_jump_beta .AND. &
                                   qs_ot_env(1)%gnorm > ot_cg_kp_projector_jump_new_norm .AND. &
                                   gnorm_ratio > ot_cg_kp_projector_jump_ratio .AND. &
                                   schur_attempt >= 3 .AND. schur_last_attempt >= 3 .AND. &
                                   schur_flavor >= 35 .AND. schur_last_flavor >= 35
      schur_projector_tail_reset = complex_kpoint_history .AND. &
                                   qs_ot_env(1)%settings%occupation_preconditioner .AND. &
                                   qs_ot_env(1)%settings%do_rotation .AND. do_ener .AND. &
                                   schur_attempt >= 3 .AND. schur_last_attempt >= 3 .AND. &
                                   schur_flavor >= 30 .AND. schur_last_flavor >= 30 .AND. &
                                   beta_pr > ot_cg_kp_projector_tail_beta .AND. &
                                   qs_ot_env(1)%gnorm > ot_cg_kp_projector_tail_norm .AND. &
                                   gnorm_ratio > ot_cg_kp_projector_tail_ratio
      residual_tail_spike_reset = complex_kpoint_history .AND. &
                                  qs_ot_env(1)%settings%occupation_preconditioner .AND. &
                                  qs_ot_env(1)%settings%do_rotation .AND. do_ener .AND. &
                                  qs_ot_env(1)%kpoint_physical_delta_ready .AND. &
                                  physical_tail_delta > residual_tail_min .AND. &
                                  physical_tail_delta < ot_cg_kp_residual_tail_max_delta .AND. &
                                  beta_pr > ot_cg_kp_residual_tail_jump_beta .AND. &
                                  qs_ot_env(1)%gnorm > ot_cg_kp_residual_tail_jump_norm .AND. &
                                  gnorm_ratio > ot_cg_kp_residual_tail_jump_ratio
      history_spike_reset = complex_kpoint_history .AND. &
                            qs_ot_env(1)%settings%occupation_preconditioner .AND. &
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
         IF (gnorm_ratio > EPSILON(1.0_dp)) &
            residual_spike_quality = MIN(residual_spike_quality, 1.0_dp/SQRT(gnorm_ratio))
         qs_ot_env(1)%schur_residual_history_reset_pending = .TRUE.
         qs_ot_env(1)%schur_residual_history_quality = &
            MIN(qs_ot_env(1)%schur_residual_history_quality, residual_spike_quality)
         IF (logger%iter_info%print_level >= high_print_level .AND. logger%para_env%is_source()) THEN
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
                               qs_ot_env(1)%settings%occupation_preconditioner .AND. &
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
         IF (kpoint_preconditioned_direction) &
            CALL prepare_linear_kpoint_energy_direction(qs_ot_env, para_env_inter_kp)
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
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) &
                  CALL dbcsr_add(qs_ot_env(ispin)%matrix_dx_im, qs_ot_env(ispin)%matrix_gx_old_im, &
                                 alpha_scalar=beta_pr, beta_scalar=-1.0_dp)
               IF (qs_ot_env(1)%settings%do_rotation) THEN
                  CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_dx, &
                                 qs_ot_env(ispin)%rot_mat_gx_old, alpha_scalar=beta_pr, beta_scalar=-1.0_dp)
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) &
                     CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_dx_im, &
                                    qs_ot_env(ispin)%rot_mat_gx_old_im, alpha_scalar=beta_pr, beta_scalar=-1.0_dp)
               END IF
            END DO
         END IF
         IF (do_ener) THEN
            DO ispin = 1, nspin
               qs_ot_env(ispin)%ener_dx = -qs_ot_env(ispin)%ener_gx_old
            END DO
            IF (kpoint_preconditioned_direction) &
               CALL prepare_linear_kpoint_energy_direction(qs_ot_env, para_env_inter_kp)
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
      schur_history_mixed = schur_history_change .AND. beta_pr > 0.0_dp
      qs_ot_env(1)%schur_history_mixed = schur_history_mixed
      IF (logger%iter_info%print_level >= high_print_level .AND. logger%para_env%is_source()) THEN
         WRITE (cp_logger_get_default_unit_nr(logger), '(A,8(1X,ES16.8),4(1X,I6),6(1X,L1))') &
            " OT CG beta raw used gnorm old cross down ratio rscale schur_att last flavor last neg "// &
            "hist_change mixed reset spike resid_reset:", &
            beta_raw, beta_pr, qs_ot_env(1)%gnorm, qs_ot_env(1)%gnorm_old, gnorm_cross, test_down, gnorm_ratio, &
            residual_trust_scale, &
            schur_attempt, schur_last_attempt, schur_flavor, schur_last_flavor, negative_beta, &
            schur_history_change, schur_history_mixed, reset_to_sd, history_spike_reset, residual_history_reset
      END IF
      ! since we change the direction we have to adjust the gradient
      qs_ot_env(1)%gradient = test_down
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
!> \brief ...
!> \param qs_ot_env ...
! **************************************************************************************************
   SUBROUTINE ot_diis_step(qs_ot_env)
      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: qs_ot_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'ot_diis_step'

      INTEGER                                            :: diis_bound, diis_m, handle, i, info, &
                                                            ispin, itmp, j, k, n, nener, nspin
      LOGICAL                                            :: do_ener, do_ks, do_ot_sd
      REAL(KIND=dp)                                      :: overlap, tmp, tr_xnew_gx, tr_xold_gx
      TYPE(cp_logger_type), POINTER                      :: logger

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()

      do_ks = qs_ot_env(1)%settings%ks
      do_ener = qs_ot_env(1)%settings%do_ener
      nspin = SIZE(qs_ot_env)

      diis_m = qs_ot_env(1)%settings%diis_m

      IF (qs_ot_env(1)%diis_iter < diis_m) THEN
         diis_bound = qs_ot_env(1)%diis_iter + 1
      ELSE
         diis_bound = diis_m
      END IF

      j = MOD(qs_ot_env(1)%diis_iter, diis_m) + 1 ! index in the circular array

      ! copy the position and the error vector in the diis buffers

      IF (do_ks) THEN
         DO ispin = 1, nspin
            CALL dbcsr_copy(qs_ot_env(ispin)%matrix_h_x(j)%matrix, qs_ot_env(ispin)%matrix_x)
            IF (qs_ot_env(ispin)%has_complex_kpoint_state) &
               CALL dbcsr_copy(qs_ot_env(ispin)%matrix_h_x_im(j)%matrix, qs_ot_env(ispin)%matrix_x_im)
            IF (qs_ot_env(ispin)%settings%do_rotation) THEN
               CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_h_x(j)%matrix, qs_ot_env(ispin)%rot_mat_x)
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) &
                  CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_h_x_im(j)%matrix, qs_ot_env(ispin)%rot_mat_x_im)
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
               CALL apply_preconditioner(qs_ot_env(ispin)%preconditioner, &
                                         qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin)%matrix_h_e(j)%matrix)
               CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin)%matrix_h_e(j)%matrix, &
                              tmp)
               qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + tmp
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  CALL apply_preconditioner(qs_ot_env(ispin)%preconditioner, &
                                            qs_ot_env(ispin)%matrix_gx_im, &
                                            qs_ot_env(ispin)%matrix_h_e_im(j)%matrix)
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
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) &
                  CALL dbcsr_scale(qs_ot_env(ispin)%matrix_h_e_im(j)%matrix, -qs_ot_env(1)%ds_min)
            END DO
            IF (qs_ot_env(1)%settings%do_rotation) THEN
               DO ispin = 1, nspin
                  CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_h_e(j)%matrix, qs_ot_env(ispin)%rot_mat_gx)
                  CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx, qs_ot_env(ispin)%rot_mat_h_e(j)%matrix, &
                                 tmp)
                  qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + 0.5_dp*tmp
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_h_e_im(j)%matrix, &
                                     qs_ot_env(ispin)%rot_mat_gx_im)
                     CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx_im, &
                                    qs_ot_env(ispin)%rot_mat_h_e_im(j)%matrix, tmp)
                     qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + 0.5_dp*tmp
                  END IF
               END DO
               DO ispin = 1, nspin
                  CALL dbcsr_scale(qs_ot_env(ispin)%rot_mat_h_e(j)%matrix, -qs_ot_env(1)%ds_min)
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) &
                     CALL dbcsr_scale(qs_ot_env(ispin)%rot_mat_h_e_im(j)%matrix, -qs_ot_env(1)%ds_min)
               END DO
            END IF
         END IF
         IF (do_ener) THEN
            DO ispin = 1, nspin
               qs_ot_env(ispin)%ener_h_e(j, :) = qs_ot_env(ispin)%ener_gx(:)
               tmp = DOT_PRODUCT(qs_ot_env(ispin)%ener_h_e(j, :), qs_ot_env(ispin)%ener_gx(:))
               qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + tmp
               qs_ot_env(ispin)%ener_h_e(j, :) = -qs_ot_env(1)%ds_min*qs_ot_env(ispin)%ener_h_e(j, :)
            END DO
         END IF
      ELSE
         qs_ot_env(1)%gnorm = 0.0_dp
         IF (do_ks) THEN
            DO ispin = 1, nspin
               CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin)%matrix_gx, tmp)
               qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + tmp
               CALL dbcsr_add(qs_ot_env(ispin)%matrix_h_e(j)%matrix, &
                              qs_ot_env(ispin)%matrix_gx, alpha_scalar=0.0_dp, beta_scalar=-qs_ot_env(1)%ds_min)
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx_im, qs_ot_env(ispin)%matrix_gx_im, tmp)
                  qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + tmp
                  CALL dbcsr_add(qs_ot_env(ispin)%matrix_h_e_im(j)%matrix, &
                                 qs_ot_env(ispin)%matrix_gx_im, alpha_scalar=0.0_dp, &
                                 beta_scalar=-qs_ot_env(1)%ds_min)
               END IF
            END DO
            IF (qs_ot_env(1)%settings%do_rotation) THEN
               DO ispin = 1, nspin
                  CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx, qs_ot_env(ispin)%rot_mat_gx, tmp)
                  qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + 0.5_dp*tmp
                  CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_h_e(j)%matrix, &
                                 qs_ot_env(ispin)%rot_mat_gx, alpha_scalar=0.0_dp, beta_scalar=-qs_ot_env(1)%ds_min)
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx_im, qs_ot_env(ispin)%rot_mat_gx_im, tmp)
                     qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + 0.5_dp*tmp
                     CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_h_e_im(j)%matrix, &
                                    qs_ot_env(ispin)%rot_mat_gx_im, alpha_scalar=0.0_dp, &
                                    beta_scalar=-qs_ot_env(1)%ds_min)
                  END IF
               END DO
            END IF
         END IF
         IF (do_ener) THEN
            DO ispin = 1, nspin
               tmp = DOT_PRODUCT(qs_ot_env(ispin)%ener_gx(:), qs_ot_env(ispin)%ener_gx(:))
               qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + tmp
               qs_ot_env(ispin)%ener_h_e(j, :) = -qs_ot_env(1)%ds_min*qs_ot_env(ispin)%ener_gx(:)
            END DO
         END IF
      END IF
      k = 0
      n = 0
      nener = 0
      IF (do_ks) THEN
         CALL dbcsr_get_info(qs_ot_env(1)%matrix_x, nfullrows_total=n)
         DO ispin = 1, nspin
            CALL dbcsr_get_info(qs_ot_env(ispin)%matrix_x, nfullcols_total=itmp)
            k = k + itmp
            IF (qs_ot_env(ispin)%has_complex_kpoint_state) k = k + itmp
         END DO
      END IF
      IF (do_ener) THEN
         DO ispin = 1, nspin
            nener = nener + SIZE(qs_ot_env(ispin)%ener_x)
         END DO
      END IF
      ! Handling the case of no free variables to optimize
      IF (INT(n, KIND=int_8)*INT(k, KIND=int_8) + nener /= 0) THEN
         qs_ot_env(1)%delta = SQRT(ABS(qs_ot_env(1)%gnorm)/(INT(n, KIND=int_8)*INT(k, KIND=int_8) + nener))
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
                  CALL dbcsr_dot(qs_ot_env(ispin)%matrix_h_e(j)%matrix, &
                                 qs_ot_env(ispin)%matrix_h_e(i)%matrix, &
                                 tmp)
                  qs_ot_env(1)%ls_diis(i, j) = qs_ot_env(1)%ls_diis(i, j) + tmp
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_dot(qs_ot_env(ispin)%matrix_h_e_im(j)%matrix, &
                                    qs_ot_env(ispin)%matrix_h_e_im(i)%matrix, tmp)
                     qs_ot_env(1)%ls_diis(i, j) = qs_ot_env(1)%ls_diis(i, j) + tmp
                  END IF
                  IF (qs_ot_env(ispin)%settings%do_rotation) THEN
                     CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_h_e(j)%matrix, &
                                    qs_ot_env(ispin)%rot_mat_h_e(i)%matrix, &
                                    tmp)
                     qs_ot_env(1)%ls_diis(i, j) = qs_ot_env(1)%ls_diis(i, j) + 0.5_dp*tmp
                     IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                        CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_h_e_im(j)%matrix, &
                                       qs_ot_env(ispin)%rot_mat_h_e_im(i)%matrix, tmp)
                        qs_ot_env(1)%ls_diis(i, j) = qs_ot_env(1)%ls_diis(i, j) + 0.5_dp*tmp
                     END IF
                  END IF
               END DO
            END IF
            IF (do_ener) THEN
               DO ispin = 1, nspin
                  tmp = DOT_PRODUCT(qs_ot_env(ispin)%ener_h_e(j, :), qs_ot_env(ispin)%ener_h_e(i, :))
                  qs_ot_env(1)%ls_diis(i, j) = qs_ot_env(1)%ls_diis(i, j) + tmp
               END DO
            END IF
         ELSE
            qs_ot_env(1)%ls_diis(i, j) = 0.0_dp
            IF (do_ks) THEN
               DO ispin = 1, nspin
                  CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx, &
                                 qs_ot_env(ispin)%matrix_h_e(i)%matrix, &
                                 tmp)
                  qs_ot_env(1)%ls_diis(i, j) = qs_ot_env(1)%ls_diis(i, j) - qs_ot_env(1)%ds_min*tmp
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx_im, &
                                    qs_ot_env(ispin)%matrix_h_e_im(i)%matrix, tmp)
                     qs_ot_env(1)%ls_diis(i, j) = qs_ot_env(1)%ls_diis(i, j) - qs_ot_env(1)%ds_min*tmp
                  END IF
                  IF (qs_ot_env(ispin)%settings%do_rotation) THEN
                     CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx, &
                                    qs_ot_env(ispin)%rot_mat_h_e(i)%matrix, &
                                    tmp)
                     qs_ot_env(1)%ls_diis(i, j) = qs_ot_env(1)%ls_diis(i, j) - qs_ot_env(1)%ds_min*0.5_dp*tmp
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
                  qs_ot_env(1)%ls_diis(i, j) = qs_ot_env(1)%ls_diis(i, j) - qs_ot_env(1)%ds_min*tmp
               END DO
            END IF
         END IF
         qs_ot_env(1)%ls_diis(j, i) = qs_ot_env(1)%ls_diis(i, j)
         qs_ot_env(1)%ls_diis(i, diis_bound + 1) = 1.0_dp
         qs_ot_env(1)%ls_diis(diis_bound + 1, i) = 1.0_dp
         qs_ot_env(1)%c_diis(i) = 0.0_dp
      END DO
      qs_ot_env(1)%ls_diis(diis_bound + 1, diis_bound + 1) = 0.0_dp
      qs_ot_env(1)%c_diis(diis_bound + 1) = 1.0_dp
      ! put in buffer, dgesv destroys
      qs_ot_env(1)%lss_diis = qs_ot_env(1)%ls_diis

      CALL DGESV(diis_bound + 1, 1, qs_ot_env(1)%lss_diis, diis_m + 1, qs_ot_env(1)%ipivot, &
                 qs_ot_env(1)%c_diis, diis_m + 1, info)

      IF (info /= 0) THEN
         do_ot_sd = .TRUE.
         WRITE (cp_logger_get_default_unit_nr(logger), *) "Singular DIIS matrix"
      ELSE
         do_ot_sd = .FALSE.
         IF (do_ks) THEN
            DO ispin = 1, nspin
               ! OK, add the vectors now
               CALL dbcsr_set(qs_ot_env(ispin)%matrix_x, 0.0_dp)
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) &
                  CALL dbcsr_set(qs_ot_env(ispin)%matrix_x_im, 0.0_dp)
               DO i = 1, diis_bound
                  CALL dbcsr_add(qs_ot_env(ispin)%matrix_x, &
                                 qs_ot_env(ispin)%matrix_h_e(i)%matrix, &
                                 alpha_scalar=1.0_dp, beta_scalar=qs_ot_env(1)%c_diis(i))
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) &
                     CALL dbcsr_add(qs_ot_env(ispin)%matrix_x_im, &
                                    qs_ot_env(ispin)%matrix_h_e_im(i)%matrix, &
                                    alpha_scalar=1.0_dp, beta_scalar=qs_ot_env(1)%c_diis(i))
               END DO
               DO i = 1, diis_bound
                  CALL dbcsr_add(qs_ot_env(ispin)%matrix_x, &
                                 qs_ot_env(ispin)%matrix_h_x(i)%matrix, &
                                 alpha_scalar=1.0_dp, beta_scalar=qs_ot_env(1)%c_diis(i))
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) &
                     CALL dbcsr_add(qs_ot_env(ispin)%matrix_x_im, &
                                    qs_ot_env(ispin)%matrix_h_x_im(i)%matrix, &
                                    alpha_scalar=1.0_dp, beta_scalar=qs_ot_env(1)%c_diis(i))
               END DO
               IF (qs_ot_env(ispin)%settings%do_rotation) THEN
                  CALL dbcsr_set(qs_ot_env(ispin)%rot_mat_x, 0.0_dp)
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) &
                     CALL dbcsr_set(qs_ot_env(ispin)%rot_mat_x_im, 0.0_dp)
                  DO i = 1, diis_bound
                     CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_x, &
                                    qs_ot_env(ispin)%rot_mat_h_e(i)%matrix, &
                                    alpha_scalar=1.0_dp, beta_scalar=qs_ot_env(1)%c_diis(i))
                     IF (qs_ot_env(ispin)%has_complex_kpoint_state) &
                        CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_x_im, &
                                       qs_ot_env(ispin)%rot_mat_h_e_im(i)%matrix, &
                                       alpha_scalar=1.0_dp, beta_scalar=qs_ot_env(1)%c_diis(i))
                  END DO
                  DO i = 1, diis_bound
                     CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_x, &
                                    qs_ot_env(ispin)%rot_mat_h_x(i)%matrix, &
                                    alpha_scalar=1.0_dp, beta_scalar=qs_ot_env(1)%c_diis(i))
                     IF (qs_ot_env(ispin)%has_complex_kpoint_state) &
                        CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_x_im, &
                                       qs_ot_env(ispin)%rot_mat_h_x_im(i)%matrix, &
                                       alpha_scalar=1.0_dp, beta_scalar=qs_ot_env(1)%c_diis(i))
                  END DO
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
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) &
                  CALL dbcsr_set(qs_ot_env(ispin)%matrix_x_im, 0.0_dp)
               CALL dbcsr_add(qs_ot_env(ispin)%matrix_x, &
                              qs_ot_env(ispin)%matrix_h_e(j)%matrix, &
                              1.0_dp, 1.0_dp)
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) &
                  CALL dbcsr_add(qs_ot_env(ispin)%matrix_x_im, &
                                 qs_ot_env(ispin)%matrix_h_e_im(j)%matrix, &
                                 1.0_dp, 1.0_dp)
               CALL dbcsr_add(qs_ot_env(ispin)%matrix_x, &
                              qs_ot_env(ispin)%matrix_h_x(j)%matrix, &
                              1.0_dp, 1.0_dp)
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) &
                  CALL dbcsr_add(qs_ot_env(ispin)%matrix_x_im, &
                                 qs_ot_env(ispin)%matrix_h_x_im(j)%matrix, &
                                 1.0_dp, 1.0_dp)
               IF (qs_ot_env(ispin)%settings%do_rotation) THEN
                  CALL dbcsr_set(qs_ot_env(ispin)%rot_mat_x, 0.0_dp)
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) &
                     CALL dbcsr_set(qs_ot_env(ispin)%rot_mat_x_im, 0.0_dp)
                  CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_x, &
                                 qs_ot_env(ispin)%rot_mat_h_e(j)%matrix, &
                                 1.0_dp, 1.0_dp)
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) &
                     CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_x_im, &
                                    qs_ot_env(ispin)%rot_mat_h_e_im(j)%matrix, &
                                    1.0_dp, 1.0_dp)
                  CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_x, &
                                 qs_ot_env(ispin)%rot_mat_h_x(j)%matrix, &
                                 1.0_dp, 1.0_dp)
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) &
                     CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_x_im, &
                                    qs_ot_env(ispin)%rot_mat_h_x_im(j)%matrix, &
                                    1.0_dp, 1.0_dp)
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
