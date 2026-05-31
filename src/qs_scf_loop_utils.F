!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
! **************************************************************************************************
!> \brief Utility routines for qs_scf
! **************************************************************************************************
MODULE qs_scf_loop_utils
   USE cp_control_types,                ONLY: dft_control_type,&
                                              hairy_probes_type
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_add, dbcsr_copy, dbcsr_create, dbcsr_deallocate_matrix, dbcsr_get_info, &
        dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, dbcsr_iterator_start, &
        dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_multiply, dbcsr_p_type, dbcsr_release, &
        dbcsr_scale, dbcsr_set, dbcsr_type, dbcsr_type_antisymmetric, dbcsr_type_symmetric
   USE cp_dbcsr_contrib,                ONLY: dbcsr_dot,&
                                              dbcsr_get_diag,&
                                              dbcsr_scale_by_vector,&
                                              dbcsr_set_diag
   USE cp_dbcsr_diag,                   ONLY: cp_dbcsr_heevd
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              cp_dbcsr_sm_fm_multiply
   USE cp_external_control,             ONLY: external_control
   USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_release,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type,&
                                              cp_to_string
   USE cp_output_handling,              ONLY: high_print_level
   USE input_constants,                 ONLY: smear_fermi_dirac,&
                                              smear_gaussian,&
                                              smear_mp,&
                                              smear_mv
   USE input_section_types,             ONLY: section_vals_type
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE kpoint_methods,                  ONLY: kpoint_density_matrices,&
                                              kpoint_density_transform,&
                                              kpoint_set_mo_occupation,&
                                              rskp_transform
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_env_type,&
                                              kpoint_type
   USE mathconstants,                   ONLY: rootpi,&
                                              sqrt2,&
                                              sqrthalf
   USE mathlib,                         ONLY: diamat_all,&
                                              diag_complex
   USE message_passing,                 ONLY: mp_para_env_type
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE qs_density_matrices,             ONLY: calculate_density_matrix
   USE qs_density_mixing_types,         ONLY: broyden_mixing_nr,&
                                              direct_mixing_nr,&
                                              gspace_mixing_nr,&
                                              modified_broyden_mixing_nr,&
                                              multisecant_mixing_nr,&
                                              no_mixing_nr,&
                                              pulay_mixing_nr
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_fb_env_methods,               ONLY: fb_env_do_diag
   USE qs_gspace_mixing,                ONLY: gspace_mixing
   USE qs_ks_types,                     ONLY: qs_ks_did_change,&
                                              qs_ks_env_type
   USE qs_mixing_utils,                 ONLY: self_consistency_check
   USE qs_mo_occupation,                ONLY: set_mo_occupation
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
   USE qs_mom_methods,                  ONLY: do_mom_diag
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_ot,                           ONLY: qs_ot_complex_exp_frechet_kernel,&
                                              qs_ot_generate_rotation_complex,&
                                              qs_ot_get_orbitals_ref_complex,&
                                              qs_ot_rot_mat_derivative_complex
   USE qs_ot_minimizer,                 ONLY: ot_mini
   USE qs_ot_scf,                       ONLY: ot_scf_destroy,&
                                              ot_scf_mini
   USE qs_ot_types,                     ONLY: qs_ot_channel_index,&
                                              qs_ot_check_channel_context,&
                                              qs_ot_type
   USE qs_outer_scf,                    ONLY: outer_loop_gradient
   USE qs_rho_methods,                  ONLY: qs_rho_update_rho
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_scf_diagonalization,          ONLY: do_block_davidson_diag,&
                                              do_block_krylov_diag,&
                                              do_general_diag,&
                                              do_general_diag_kp,&
                                              do_ot_diag,&
                                              do_roks_diag,&
                                              do_scf_diag_subspace,&
                                              do_special_diag
   USE qs_scf_methods,                  ONLY: scf_env_density_mixing
   USE qs_scf_output,                   ONLY: qs_scf_print_summary
   USE qs_scf_types,                    ONLY: &
        block_davidson_diag_method_nr, block_krylov_diag_method_nr, filter_matrix_diag_method_nr, &
        general_diag_method_nr, ot_diag_method_nr, ot_method_nr, qs_scf_env_type, &
        smeagol_method_nr, special_diag_method_nr
   USE scf_control_types,               ONLY: scf_control_type,&
                                              smear_type
   USE smeagol_interface,               ONLY: run_smeagol_emtrans
   USE smearing_utils,                  ONLY: Smearkp,&
                                              Smearkp2
   USE tblite_interface,                ONLY: tb_native_scc_mixer_active
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_scf_loop_utils'
   REAL(KIND=dp), PARAMETER, PRIVATE :: ot_kp_canonical_density_max_weight = 2.0E-2_dp, &
                                        ot_kp_canonical_projector_delta = 5.0E-3_dp, &
                                        ot_kp_canonical_projector_max_weight = 0.10_dp, &
                                        ot_kp_fermi_window_tail_fraction = 3.22954647E-2_dp, &
                                        ot_kp_fermi_window_min_frac = 1.0E-3_dp, &
                                        ot_kp_fermi_window_outer_factor = 8.0_dp, &
                                        ot_kp_schur_response_trust_min_factor = 1.0E-1_dp
   INTEGER, PARAMETER, PRIVATE :: ot_kp_schur_quality_min_samples = 3
   REAL(KIND=dp), PARAMETER, PRIVATE :: ot_kp_schur_quality_ema_alpha = 0.30_dp, &
                                        ot_kp_schur_quality_max_ratio = 2.0_dp, &
                                        ot_kp_schur_quality_min_ratio = 1.0E-1_dp
   INTEGER, PARAMETER, PRIVATE :: ot_kp_schur_calibration_min_samples = 4
   INTEGER, PARAMETER, PRIVATE :: ot_kp_schur_constrained_min_samples = 6
   INTEGER, PARAMETER, PRIVATE :: ot_kp_schur_projector_response_min_samples = 8
   REAL(KIND=dp), PARAMETER, PRIVATE :: ot_kp_schur_calibration_floor = 2.5E-1_dp, &
                                        ot_kp_schur_calibration_max_scale = 1.2_dp, &
                                        ot_kp_schur_calibrated_quality_max_ratio = 2.5_dp, &
                                        ot_kp_schur_calibrated_quality_min_ratio = 2.5E-1_dp, &
                                        ot_kp_schur_constrained_advantage_ratio = 1.05_dp, &
                                        ot_kp_schur_constrained_drop_advantage_ratio = 1.02_dp, &
                                        ot_kp_schur_projector_response_advantage_ratio = 1.10_dp, &
                                        ot_kp_schur_projector_response_floor_delta = 3.0E-6_dp, &
                                        ot_kp_schur_projector_response_full_delta = 1.2E-5_dp, &
                                        ot_kp_schur_projector_response_quality_max_ratio = 1.25_dp, &
                                        ot_kp_schur_projector_response_quality_min_ratio = 1.0E-1_dp, &
                                        ot_kp_schur_projector_response_max_weight = 1.0E-1_dp, &
                                        ot_kp_schur_projector_response_ratio_gate = 2.0E-2_dp, &
                                        ot_kp_schur_projector_mode_trust_max_penalty = 2.5E-1_dp, &
                                        ot_kp_schur_projector_coupled_tail_max_weight = 0.10_dp, &
                                        ot_kp_schur_projector_coupled_tail_delta = 1.0E-4_dp

   PUBLIC :: qs_scf_set_loop_flags, &
             qs_scf_new_mos, qs_scf_new_mos_kp, &
             qs_scf_density_mixing, qs_scf_check_inner_exit, &
             qs_scf_check_outer_exit, qs_scf_inner_finalize, qs_scf_rho_update

CONTAINS

! **************************************************************************************************
!> \brief computes properties for a given hamiltonian using the current wfn
!> \param scf_env ...
!> \param diis_step ...
!> \param energy_only ...
!> \param just_energy ...
!> \param exit_inner_loop ...
! **************************************************************************************************
   SUBROUTINE qs_scf_set_loop_flags(scf_env, diis_step, &
                                    energy_only, just_energy, exit_inner_loop)

      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      LOGICAL                                            :: diis_step, energy_only, just_energy, &
                                                            exit_inner_loop

! Some flags needed to be set at the beginning of the loop

      diis_step = .FALSE.
      energy_only = .FALSE.
      just_energy = .FALSE.

      ! SCF loop, optimisation of the wfn coefficients
      ! qs_env%rho%rho_r and qs_env%rho%rho_g should be up to date here

      scf_env%iter_count = 0
      exit_inner_loop = .FALSE.

   END SUBROUTINE qs_scf_set_loop_flags

! **************************************************************************************************
!> \brief takes known energy and derivatives and produces new wfns
!>        and or density matrix
!> \param qs_env ...
!> \param scf_env ...
!> \param scf_control ...
!> \param scf_section ...
!> \param diis_step ...
!> \param energy_only ...
!> \param probe ...
! **************************************************************************************************
   SUBROUTINE qs_scf_new_mos(qs_env, scf_env, scf_control, scf_section, diis_step, &
                             energy_only, probe, added_mos_auto_grow)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(section_vals_type), POINTER                   :: scf_section
      LOGICAL                                            :: diis_step, energy_only
      TYPE(hairy_probes_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: probe
      LOGICAL, INTENT(OUT), OPTIONAL                     :: added_mos_auto_grow

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'qs_scf_new_mos'

      INTEGER                                            :: handle, ispin
      LOGICAL                                            :: disable_diis, has_unit_metric, &
                                                            my_added_mos_auto_grow, skip_diag_sub
      REAL(KIND=dp)                                      :: saved_eps_diis
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho

      CALL timeset(routineN, handle)

      NULLIFY (energy, ks_env, matrix_ks, matrix_s, rho, mos, dft_control)
      my_added_mos_auto_grow = .FALSE.
      IF (PRESENT(added_mos_auto_grow)) added_mos_auto_grow = .FALSE.

      CALL get_qs_env(qs_env=qs_env, &
                      matrix_s=matrix_s, energy=energy, &
                      ks_env=ks_env, &
                      matrix_ks=matrix_ks, rho=rho, mos=mos, &
                      dft_control=dft_control, &
                      has_unit_metric=has_unit_metric)
      scf_env%iter_param = 0.0_dp
      disable_diis = dft_control%qs_control%xtb_control%do_tblite .AND. &
                     tb_native_scc_mixer_active(dft_control)
      IF (disable_diis) THEN
         saved_eps_diis = scf_control%eps_diis
         scf_control%eps_diis = 0.0_dp
      END IF

      ! transfer total_zeff_corr from qs_env to scf_env only if
      ! correct_el_density_dip is switched on [SGh]
      IF (dft_control%correct_el_density_dip) THEN
         scf_env%sum_zeff_corr = qs_env%total_zeff_corr
         IF (ABS(qs_env%total_zeff_corr) > 0.0_dp) THEN
            IF (scf_env%method /= general_diag_method_nr) THEN
               CALL cp_abort(__LOCATION__, &
                             "Please use ALGORITHM STANDARD in "// &
                             "SCF%DIAGONALIZATION if "// &
                             "CORE_CORRECTION /= 0.0 and "// &
                             "SURFACE_DIPOLE_CORRECTION TRUE ")
            ELSEIF (dft_control%roks) THEN
               CALL cp_abort(__LOCATION__, &
                             "Combination of "// &
                             "CORE_CORRECTION /= 0.0 and "// &
                             "SURFACE_DIPOLE_CORRECTION TRUE "// &
                             "is not implemented with ROKS")
            ELSEIF (scf_control%diagonalization%mom) THEN
               CALL cp_abort(__LOCATION__, &
                             "Combination of "// &
                             "CORE_CORRECTION /= 0.0 and "// &
                             "SURFACE_DIPOLE_CORRECTION TRUE "// &
                             "is not implemented with SCF%MOM")
            END IF
         END IF
      END IF

      SELECT CASE (scf_env%method)
      CASE DEFAULT
         CALL cp_abort(__LOCATION__, &
                       "unknown scf method: "// &
                       cp_to_string(scf_env%method))

         ! *************************************************************************
         ! Filter matrix diagonalisation: ugly implementation at this point of time
         ! *************************************************************************
      CASE (filter_matrix_diag_method_nr)

         IF (ABS(qs_env%total_zeff_corr) > 0.0_dp) THEN
            CALL cp_abort(__LOCATION__, &
                          "CORE_CORRECTION /= 0.0 plus SURFACE_DIPOLE_CORRECTION TRUE "// &
                          "requires SCF%DIAGONALIZATION: ALGORITHM STANDARD")
         END IF
         CALL fb_env_do_diag(scf_env%filter_matrix_env, qs_env, &
                             matrix_ks, matrix_s, scf_section, diis_step)

         ! Diagonlization in non orthonormal case
      CASE (general_diag_method_nr)
         IF (dft_control%roks) THEN
            CALL do_roks_diag(scf_env, mos, matrix_ks, matrix_s, &
                              scf_control, scf_section, diis_step, &
                              has_unit_metric)
         ELSE
            IF (scf_control%diagonalization%mom) THEN
               CALL do_mom_diag(scf_env, mos, matrix_ks, &
                                matrix_s, scf_control, scf_section, &
                                diis_step)
            ELSE
               IF (dft_control%hairy_probes .EQV. .TRUE.) THEN
                  CALL do_general_diag(scf_env, mos, matrix_ks, &
                                       matrix_s, scf_control, scf_section, &
                                       diis_step, &
                                       probe)
               ELSE
                  CALL do_general_diag(scf_env, mos, matrix_ks, &
                                       matrix_s, scf_control, scf_section, &
                                       diis_step)
               END IF
            END IF
            IF (scf_control%do_diag_sub) THEN
               skip_diag_sub = (scf_env%subspace_env%eps_diag_sub > 0.0_dp) .AND. &
                               (scf_env%iter_count == 1 .OR. scf_env%iter_delta > scf_env%subspace_env%eps_diag_sub)
               IF (.NOT. skip_diag_sub) THEN
                  CALL do_scf_diag_subspace(qs_env, scf_env, scf_env%subspace_env, mos, rho, &
                                            ks_env, scf_section, scf_control)
               END IF
            END IF
         END IF
         ! Diagonlization in orthonormal case
      CASE (special_diag_method_nr)
         IF (dft_control%roks) THEN
            CALL do_roks_diag(scf_env, mos, matrix_ks, matrix_s, &
                              scf_control, scf_section, diis_step, &
                              has_unit_metric)
         ELSE
            CALL do_special_diag(scf_env, mos, matrix_ks, &
                                 scf_control, scf_section, &
                                 diis_step)
         END IF
         ! OT diagonalization
      CASE (ot_diag_method_nr)
         CALL do_ot_diag(scf_env, mos, matrix_ks, matrix_s, &
                         scf_control, scf_section, diis_step)
         ! Block Krylov diagonlization
      CASE (block_krylov_diag_method_nr)
         IF ((scf_env%krylov_space%eps_std_diag > 0.0_dp) .AND. &
             (scf_env%iter_count == 1 .OR. scf_env%iter_delta > scf_env%krylov_space%eps_std_diag)) THEN
            IF (scf_env%krylov_space%always_check_conv) THEN
               CALL do_block_krylov_diag(scf_env, mos, matrix_ks, &
                                         scf_control, scf_section, check_moconv_only=.TRUE.)
            END IF
            CALL do_general_diag(scf_env, mos, matrix_ks, &
                                 matrix_s, scf_control, scf_section, diis_step)
         ELSE
            CALL do_block_krylov_diag(scf_env, mos, matrix_ks, &
                                      scf_control, scf_section)
         END IF
         IF (scf_control%do_diag_sub) THEN
            skip_diag_sub = (scf_env%subspace_env%eps_diag_sub > 0.0_dp) .AND. &
                            (scf_env%iter_count == 1 .OR. scf_env%iter_delta > scf_env%subspace_env%eps_diag_sub)
            IF (.NOT. skip_diag_sub) THEN
               CALL do_scf_diag_subspace(qs_env, scf_env, scf_env%subspace_env, mos, rho, &
                                         ks_env, scf_section, scf_control)
            END IF
         END IF
         ! Block Davidson diagonlization
      CASE (block_davidson_diag_method_nr)
         CALL do_block_davidson_diag(qs_env, scf_env, mos, matrix_ks, matrix_s, scf_control, &
                                     scf_section, .FALSE.)
         ! OT without diagonlization. Needs special treatment for SCP runs
      CASE (ot_method_nr)
         CALL qs_scf_loop_do_ot(qs_env, scf_env, scf_control%smear, mos, rho, &
                                qs_env%mo_derivs, energy%total, &
                                matrix_s, energy_only=energy_only, has_unit_metric=has_unit_metric, &
                                added_mos_auto_grow=my_added_mos_auto_grow)
      END SELECT
      IF (disable_diis) scf_control%eps_diis = saved_eps_diis
      IF (my_added_mos_auto_grow) THEN
         IF (PRESENT(added_mos_auto_grow)) added_mos_auto_grow = .TRUE.
         CALL timestop(handle)
         RETURN
      END IF

      energy%kTS = 0.0_dp
      energy%efermi = 0.0_dp
      CALL get_qs_env(qs_env, mos=mos)
      DO ispin = 1, SIZE(mos)
         energy%kTS = energy%kTS + mos(ispin)%kTS
         energy%efermi = energy%efermi + mos(ispin)%mu
      END DO
      energy%efermi = energy%efermi/REAL(SIZE(mos), KIND=dp)

      CALL timestop(handle)

   END SUBROUTINE qs_scf_new_mos

! **************************************************************************************************
!> \brief Updates MOs and density matrix using diagonalization
!>        Kpoint code
!> \param qs_env ...
!> \param scf_env ...
!> \param scf_control ...
!> \param diis_step ...
!> \param probe ...
! **************************************************************************************************
   SUBROUTINE qs_scf_new_mos_kp(qs_env, scf_env, scf_control, diis_step, probe, added_mos_auto_grow)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(scf_control_type), POINTER                    :: scf_control
      LOGICAL                                            :: diis_step
      TYPE(hairy_probes_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: probe
      LOGICAL, INTENT(OUT), OPTIONAL                     :: added_mos_auto_grow

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'qs_scf_new_mos_kp'

      INTEGER                                            :: handle, ispin
      LOGICAL                                            :: disable_diis, has_unit_metric, &
                                                            my_added_mos_auto_grow
      REAL(dp)                                           :: diis_error, saved_eps_diis
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks, matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mo_set_type), DIMENSION(:, :), POINTER        :: mos
      TYPE(qs_energy_type), POINTER                      :: energy

      CALL timeset(routineN, handle)

      NULLIFY (dft_control, kpoints, matrix_ks, matrix_s)
      my_added_mos_auto_grow = .FALSE.
      IF (PRESENT(added_mos_auto_grow)) added_mos_auto_grow = .FALSE.

      CALL get_qs_env(qs_env=qs_env, dft_control=dft_control, kpoints=kpoints)
      scf_env%iter_param = 0.0_dp
      disable_diis = dft_control%qs_control%xtb_control%do_tblite .AND. &
                     tb_native_scc_mixer_active(dft_control)
      IF (disable_diis) THEN
         saved_eps_diis = scf_control%eps_diis
         scf_control%eps_diis = 0.0_dp
      END IF

      IF (dft_control%roks) &
         CPABORT("KP code: ROKS method not available: ")

      SELECT CASE (scf_env%method)
      CASE DEFAULT
         CALL cp_abort(__LOCATION__, &
                       "KP code: Unknown scf method: "// &
                       cp_to_string(scf_env%method))
      CASE (general_diag_method_nr)
         ! Diagonlization in non orthonormal case
         CALL get_qs_env(qs_env, matrix_ks_kp=matrix_ks, matrix_s_kp=matrix_s)
         IF (dft_control%hairy_probes .EQV. .TRUE.) THEN
            scf_control%smear%do_smear = .FALSE.
            CALL do_general_diag_kp(matrix_ks, matrix_s, kpoints, scf_env, scf_control, .TRUE., &
                                    diis_step, diis_error, qs_env, probe)
         ELSE
            CALL do_general_diag_kp(matrix_ks, matrix_s, kpoints, scf_env, scf_control, .TRUE., &
                                    diis_step, diis_error, qs_env, &
                                    added_mos_auto_grow=my_added_mos_auto_grow)
         END IF
         IF (my_added_mos_auto_grow) THEN
            IF (PRESENT(added_mos_auto_grow)) added_mos_auto_grow = .TRUE.
            IF (disable_diis) scf_control%eps_diis = saved_eps_diis
            CALL timestop(handle)
            RETURN
         END IF
         IF (diis_step) THEN
            scf_env%iter_param = diis_error
            scf_env%iter_method = "DIIS/Diag."
         ELSE
            IF (scf_env%mixing_method == 0) THEN
               scf_env%iter_method = "NoMix/Diag."
            ELSE IF (scf_env%mixing_method == 1) THEN
               scf_env%iter_param = scf_env%p_mix_alpha
               scf_env%iter_method = "P_Mix/Diag."
            ELSEIF (scf_env%mixing_method > 1) THEN
               scf_env%iter_param = scf_env%mixing_store%alpha
               scf_env%iter_method = TRIM(scf_env%mixing_store%iter_method)//"/Diag."
            END IF
         END IF
      CASE (special_diag_method_nr)
         CALL get_qs_env(qs_env=qs_env, has_unit_metric=has_unit_metric)
         CPASSERT(has_unit_metric)
         ! Diagonlization in orthonormal case
         CALL cp_abort(__LOCATION__, &
                       "KP code: Scf method not available: "// &
                       cp_to_string(scf_env%method))
      CASE (ot_diag_method_nr, &
            block_krylov_diag_method_nr, &
            block_davidson_diag_method_nr)
         CALL cp_abort(__LOCATION__, &
                       "KP code: Scf method not available: "// &
                       cp_to_string(scf_env%method))
      CASE (ot_method_nr)
         CALL get_qs_env(qs_env, matrix_ks_kp=matrix_ks, matrix_s_kp=matrix_s)
         CALL qs_scf_loop_do_ot_kp(qs_env, scf_env, scf_control, matrix_ks, matrix_s, &
                                   added_mos_auto_grow=my_added_mos_auto_grow)
         IF (my_added_mos_auto_grow) THEN
            IF (PRESENT(added_mos_auto_grow)) added_mos_auto_grow = .TRUE.
            IF (disable_diis) scf_control%eps_diis = saved_eps_diis
            CALL timestop(handle)
            RETURN
         END IF
      CASE (smeagol_method_nr)
         ! SMEAGOL interface
         diis_step = .FALSE.
         IF (scf_env%mixing_method == 0) THEN
            scf_env%iter_method = "NoMix/SMGL"
         ELSE IF (scf_env%mixing_method == 1) THEN
            scf_env%iter_param = scf_env%p_mix_alpha
            scf_env%iter_method = "P_Mix/SMGL"
         ELSE IF (scf_env%mixing_method > 1) THEN
            scf_env%iter_param = scf_env%mixing_store%alpha
            scf_env%iter_method = TRIM(scf_env%mixing_store%iter_method)//"/SMGL"
         END IF
         CALL run_smeagol_emtrans(qs_env, last=.FALSE., iter=scf_env%iter_count, rho_ao_kp=scf_env%p_mix_new)
      END SELECT
      IF (disable_diis) scf_control%eps_diis = saved_eps_diis

      CALL get_qs_env(qs_env=qs_env, energy=energy)
      energy%kTS = 0.0_dp
      energy%efermi = 0.0_dp
      mos => kpoints%kp_env(1)%kpoint_env%mos
      DO ispin = 1, SIZE(mos, 2)
         energy%kTS = energy%kTS + mos(1, ispin)%kTS
         energy%efermi = energy%efermi + mos(1, ispin)%mu
      END DO
      energy%efermi = energy%efermi/REAL(SIZE(mos, 2), KIND=dp)

      CALL timestop(handle)

   END SUBROUTINE qs_scf_new_mos_kp

! **************************************************************************************************
!> \brief K-point OT driver for local complex spin/k-point channels
!> \param qs_env ...
!> \param scf_env ...
!> \param scf_control ...
!> \param matrix_ks ...
!> \param matrix_s ...
!> \param added_mos_auto_grow ...
! **************************************************************************************************
   SUBROUTINE qs_scf_loop_do_ot_kp(qs_env, scf_env, scf_control, matrix_ks, matrix_s, &
                                   added_mos_auto_grow)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks, matrix_s
      LOGICAL, INTENT(OUT)                               :: added_mos_auto_grow

      CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_scf_loop_do_ot_kp'

      INTEGER                                            :: first_channel, handle, ikpoint, ispin, &
                                                            local_channel, local_kpoint, nkpoint, &
                                                            nlocal_channels, nmo, nspin_ot
      INTEGER, DIMENSION(2)                              :: kp_range
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: fermi_window_guarded, use_real_wfn
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: orbital_scaling_factor
      REAL(KIND=dp)                                      :: ener_nondiag, maxocc, mu
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: trial_mu
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvalues, occupation_numbers, wkp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: trial_eigenvalues, trial_occupation
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: xkp
      TYPE(cp_fm_type)                                   :: hc_im, hc_re
      TYPE(cp_fm_type), DIMENSION(:), POINTER            :: fmwork
      TYPE(cp_fm_type), POINTER                          :: mo_coeff, mo_coeff_im
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_hc, matrix_hc_im
      TYPE(dbcsr_type)                                   :: density_hc_im, density_hc_re
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: density_target, rho_ao_kp
      TYPE(dbcsr_type), POINTER                          :: matrix_k_im, matrix_k_re, matrix_s_im, &
                                                            matrix_s_re
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_env_type), POINTER                     :: kp
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mp_para_env_type), POINTER                    :: para_env_inter_kp
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: local_ot_env
      TYPE(qs_rho_type), POINTER                         :: rho

      CALL timeset(routineN, handle)

      NULLIFY (cell_to_index, density_target, dft_control, eigenvalues, energy, fmwork, kpoints, local_ot_env, matrix_hc, &
               matrix_hc_im, matrix_k_im, matrix_k_re, matrix_s_im, matrix_s_re, mo_coeff, &
               mo_coeff_im, occupation_numbers, para_env_inter_kp, rho, rho_ao_kp, sab_nl, wkp, xkp)
      added_mos_auto_grow = .FALSE.

      CPASSERT(ASSOCIATED(scf_env%qs_ot_env))
      CALL get_qs_env(qs_env=qs_env, dft_control=dft_control, energy=energy, kpoints=kpoints, rho=rho)
      CALL get_kpoint_info(kpoints, nkp=nkpoint, xkp=xkp, wkp=wkp, use_real_wfn=use_real_wfn, &
                           kp_range=kp_range, para_env_inter_kp=para_env_inter_kp, &
                           cell_to_index=cell_to_index, sab_nl=sab_nl)
      CPASSERT(ASSOCIATED(sab_nl))
      CPASSERT(ASSOCIATED(matrix_ks))
      CPASSERT(ASSOCIATED(matrix_s))

      IF (use_real_wfn) THEN
         CALL cp_abort(__LOCATION__, &
                       "K-point OT driver currently supports complex k-point wavefunctions only.")
      END IF
      IF (scf_env%qs_ot_env(1)%settings%ot_algorithm /= "REF") THEN
         CALL cp_abort(__LOCATION__, &
                       "K-point OT driver currently supports OT%ALGORITHM REF only.")
      END IF
      kp => kpoints%kp_env(1)%kpoint_env
      nspin_ot = SIZE(kp%mos, 2)
      CALL qs_ot_check_channel_context(scf_env%qs_ot_env, nspin_ot, &
                                       nkpoint=nkpoint, &
                                       restricted=ANY(scf_env%qs_ot_env(:)%restricted), &
                                       require_kpoint=.TRUE., &
                                       kp_range=kp_range, &
                                       wkp=wkp, &
                                       require_local_state=.TRUE., &
                                       require_complex_state=.TRUE.)

      first_channel = qs_ot_channel_index(1, kp_range(1), nspin_ot)
      nlocal_channels = SIZE(kpoints%kp_env)*nspin_ot
      local_ot_env => scf_env%qs_ot_env(first_channel:first_channel + nlocal_channels - 1)
      ener_nondiag = 0.0_dp

      ALLOCATE (matrix_hc(nlocal_channels))
      ALLOCATE (matrix_hc_im(nlocal_channels))
      DO local_channel = 1, nlocal_channels
         NULLIFY (matrix_hc(local_channel)%matrix)
         NULLIFY (matrix_hc_im(local_channel)%matrix)
         ALLOCATE (matrix_hc(local_channel)%matrix)
         ALLOCATE (matrix_hc_im(local_channel)%matrix)
         CALL dbcsr_copy(matrix_hc(local_channel)%matrix, local_ot_env(local_channel)%matrix_x, name="matrix_hc")
         CALL dbcsr_copy(matrix_hc_im(local_channel)%matrix, local_ot_env(local_channel)%matrix_x_im, &
                         name="matrix_hc_im")
         CALL dbcsr_set(matrix_hc(local_channel)%matrix, 0.0_dp)
         CALL dbcsr_set(matrix_hc_im(local_channel)%matrix, 0.0_dp)
      END DO

      ALLOCATE (matrix_k_re, matrix_k_im, matrix_s_re, matrix_s_im)
      CALL dbcsr_create(matrix_k_re, template=matrix_ks(1, 1)%matrix, &
                        matrix_type=dbcsr_type_symmetric)
      CALL dbcsr_create(matrix_k_im, template=matrix_ks(1, 1)%matrix, &
                        matrix_type=dbcsr_type_antisymmetric)
      CALL dbcsr_create(matrix_s_re, template=matrix_s(1, 1)%matrix, &
                        matrix_type=dbcsr_type_symmetric)
      CALL dbcsr_create(matrix_s_im, template=matrix_s(1, 1)%matrix, &
                        matrix_type=dbcsr_type_antisymmetric)
      CALL cp_dbcsr_alloc_block_from_nbl(matrix_k_re, sab_nl)
      CALL cp_dbcsr_alloc_block_from_nbl(matrix_k_im, sab_nl)
      CALL cp_dbcsr_alloc_block_from_nbl(matrix_s_re, sab_nl)
      CALL cp_dbcsr_alloc_block_from_nbl(matrix_s_im, sab_nl)

      DO local_kpoint = 1, SIZE(kpoints%kp_env)
         kp => kpoints%kp_env(local_kpoint)%kpoint_env
         ikpoint = kp%nkpoint
         CALL dbcsr_set(matrix_s_re, 0.0_dp)
         CALL dbcsr_set(matrix_s_im, 0.0_dp)
         CALL rskp_transform(rmatrix=matrix_s_re, cmatrix=matrix_s_im, rsmat=matrix_s, &
                             ispin=1, xkp=xkp(1:3, ikpoint), &
                             cell_to_index=cell_to_index, sab_nl=sab_nl)
         DO ispin = 1, nspin_ot
            local_channel = qs_ot_channel_index(ispin, ikpoint, nspin_ot) - first_channel + 1
            CALL get_mo_set(kp%mos(1, ispin), mo_coeff=mo_coeff)
            CALL get_mo_set(kp%mos(2, ispin), mo_coeff=mo_coeff_im)
            CALL qs_ot_get_orbitals_ref_complex(matrix_hc(local_channel)%matrix, &
                                                matrix_hc_im(local_channel)%matrix, &
                                                matrix_s_re, matrix_s_im, &
                                                local_ot_env(local_channel))
            CALL copy_dbcsr_to_fm(matrix_hc(local_channel)%matrix, mo_coeff)
            CALL copy_dbcsr_to_fm(matrix_hc_im(local_channel)%matrix, mo_coeff_im)
         END DO
      END DO

      DO local_kpoint = 1, SIZE(kpoints%kp_env)
         kp => kpoints%kp_env(local_kpoint)%kpoint_env
         ikpoint = kp%nkpoint
         DO ispin = 1, nspin_ot
            local_channel = qs_ot_channel_index(ispin, ikpoint, nspin_ot) - first_channel + 1
            CALL dbcsr_set(matrix_k_re, 0.0_dp)
            CALL dbcsr_set(matrix_k_im, 0.0_dp)
            CALL rskp_transform(rmatrix=matrix_k_re, cmatrix=matrix_k_im, rsmat=matrix_ks, &
                                ispin=ispin, xkp=xkp(1:3, ikpoint), &
                                cell_to_index=cell_to_index, sab_nl=sab_nl)

            CALL get_mo_set(kp%mos(1, ispin), maxocc=maxocc, mo_coeff=mo_coeff, nmo=nmo)
            CALL get_mo_set(kp%mos(2, ispin), mo_coeff=mo_coeff_im)
            CALL cp_fm_create(hc_re, mo_coeff%matrix_struct)
            CALL cp_fm_create(hc_im, mo_coeff%matrix_struct)

            CALL cp_dbcsr_sm_fm_multiply(matrix_k_re, mo_coeff, hc_re, nmo)
            CALL cp_dbcsr_sm_fm_multiply(matrix_k_im, mo_coeff_im, hc_re, nmo, &
                                         alpha=-1.0_dp, beta=1.0_dp)
            CALL cp_dbcsr_sm_fm_multiply(matrix_k_re, mo_coeff_im, hc_im, nmo)
            CALL cp_dbcsr_sm_fm_multiply(matrix_k_im, mo_coeff, hc_im, nmo, &
                                         alpha=1.0_dp, beta=1.0_dp)

            CALL copy_fm_to_dbcsr(hc_re, matrix_hc(local_channel)%matrix)
            CALL copy_fm_to_dbcsr(hc_im, matrix_hc_im(local_channel)%matrix)

            CALL qs_scf_loop_project_ot_kp_energy_subproblem(local_ot_env(local_channel), &
                                                             matrix_hc(local_channel)%matrix, &
                                                             matrix_hc_im(local_channel)%matrix, &
                                                             scf_control%smear)

            CALL cp_fm_release(hc_re)
            CALL cp_fm_release(hc_im)
         END DO
      END DO

      CALL qs_scf_loop_collect_ot_kp_trial_smearing(kpoints, local_ot_env, first_channel, nspin_ot, &
                                                    kp_range, wkp, para_env_inter_kp, scf_control%smear, &
                                                    trial_eigenvalues, trial_occupation, trial_mu)
      CALL qs_scf_loop_guard_ot_kp_trial_fermi_window(kpoints, local_ot_env, first_channel, &
                                                      nspin_ot, kp_range, wkp, para_env_inter_kp, &
                                                      scf_control%smear, trial_eigenvalues, &
                                                      trial_occupation, trial_mu, fermi_window_guarded)
      IF (fermi_window_guarded) THEN
         DEALLOCATE (trial_eigenvalues, trial_occupation, trial_mu)
         CALL qs_scf_loop_collect_ot_kp_trial_smearing(kpoints, local_ot_env, first_channel, nspin_ot, &
                                                       kp_range, wkp, para_env_inter_kp, scf_control%smear, &
                                                       trial_eigenvalues, trial_occupation, trial_mu)
      END IF

      DO local_kpoint = 1, SIZE(kpoints%kp_env)
         kp => kpoints%kp_env(local_kpoint)%kpoint_env
         ikpoint = kp%nkpoint
         DO ispin = 1, nspin_ot
            local_channel = qs_ot_channel_index(ispin, ikpoint, nspin_ot) - first_channel + 1
            CALL get_mo_set(kp%mos(1, ispin), maxocc=maxocc, nmo=nmo)

            CALL qs_scf_loop_prepare_ot_kp_channel_gradient(local_ot_env(local_channel), &
                                                            matrix_hc(local_channel)%matrix, &
                                                            matrix_hc_im(local_channel)%matrix, &
                                                            trial_occupation(1:nmo, ikpoint, ispin), &
                                                            trial_eigenvalues(1:nmo, ikpoint, ispin), &
                                                            maxocc, trial_mu(ispin), &
                                                            wkp(ikpoint), scf_control%smear, &
                                                            ener_nondiag)

            ALLOCATE (orbital_scaling_factor(nmo))
            ! Keep orbital gradients on the fixed-N Mermin surface. The occupation
            ! preconditioner below acts on energy variables, not by unweighting empty bands.
            orbital_scaling_factor(:) = 2.0_dp*wkp(ikpoint)*trial_occupation(1:nmo, ikpoint, ispin)
            CALL dbcsr_scale_by_vector(matrix_hc(local_channel)%matrix, &
                                       alpha=orbital_scaling_factor, side='right')
            CALL dbcsr_scale_by_vector(matrix_hc_im(local_channel)%matrix, &
                                       alpha=orbital_scaling_factor, side='right')
            DEALLOCATE (orbital_scaling_factor)
         END DO
      END DO

      local_ot_env(1)%etotal = energy%total + ener_nondiag

      CALL qs_scf_loop_prepare_ot_kp_energy_gradient(kpoints, local_ot_env, first_channel, nspin_ot, &
                                                     kp_range, wkp, para_env_inter_kp, scf_control%smear, &
                                                     trial_eigenvalues, trial_occupation, trial_mu)
      CALL qs_scf_loop_prepare_ot_kp_energy_trust(kpoints, local_ot_env, first_channel, nspin_ot, &
                                                  wkp, scf_control%smear, &
                                                  trial_eigenvalues, trial_occupation, trial_mu)

      CALL ot_mini(local_ot_env, matrix_hc, matrix_hc_im=matrix_hc_im)
      scf_env%iter_delta = local_ot_env(1)%delta

      DO local_kpoint = 1, SIZE(kpoints%kp_env)
         kp => kpoints%kp_env(local_kpoint)%kpoint_env
         ikpoint = kp%nkpoint
         CALL dbcsr_set(matrix_s_re, 0.0_dp)
         CALL dbcsr_set(matrix_s_im, 0.0_dp)
         CALL rskp_transform(rmatrix=matrix_s_re, cmatrix=matrix_s_im, rsmat=matrix_s, &
                             ispin=1, xkp=xkp(1:3, ikpoint), &
                             cell_to_index=cell_to_index, sab_nl=sab_nl)
         DO ispin = 1, nspin_ot
            local_channel = qs_ot_channel_index(ispin, ikpoint, nspin_ot) - first_channel + 1
            CALL get_mo_set(kp%mos(1, ispin), mo_coeff=mo_coeff, eigenvalues=eigenvalues, nmo=nmo)
            CALL get_mo_set(kp%mos(2, ispin), mo_coeff=mo_coeff_im)
            CALL qs_ot_get_orbitals_ref_complex(matrix_hc(local_channel)%matrix, &
                                                matrix_hc_im(local_channel)%matrix, &
                                                matrix_s_re, matrix_s_im, &
                                                local_ot_env(local_channel))
            CALL copy_dbcsr_to_fm(matrix_hc(local_channel)%matrix, mo_coeff)
            CALL copy_dbcsr_to_fm(matrix_hc_im(local_channel)%matrix, mo_coeff_im)
            IF (kpoint_ot_project_energy_subproblem(local_ot_env(local_channel), scf_control%smear)) THEN
               CALL dbcsr_set(matrix_k_re, 0.0_dp)
               CALL dbcsr_set(matrix_k_im, 0.0_dp)
               CALL rskp_transform(rmatrix=matrix_k_re, cmatrix=matrix_k_im, rsmat=matrix_ks, &
                                   ispin=ispin, xkp=xkp(1:3, ikpoint), &
                                   cell_to_index=cell_to_index, sab_nl=sab_nl)
               CALL cp_fm_create(hc_re, mo_coeff%matrix_struct)
               CALL cp_fm_create(hc_im, mo_coeff%matrix_struct)

               CALL cp_dbcsr_sm_fm_multiply(matrix_k_re, mo_coeff, hc_re, nmo)
               CALL cp_dbcsr_sm_fm_multiply(matrix_k_im, mo_coeff_im, hc_re, nmo, &
                                            alpha=-1.0_dp, beta=1.0_dp)
               CALL cp_dbcsr_sm_fm_multiply(matrix_k_re, mo_coeff_im, hc_im, nmo)
               CALL cp_dbcsr_sm_fm_multiply(matrix_k_im, mo_coeff, hc_im, nmo, &
                                            alpha=1.0_dp, beta=1.0_dp)

               CALL dbcsr_copy(density_hc_re, matrix_hc(local_channel)%matrix, name="ot_kp_density_hc_re")
               CALL dbcsr_copy(density_hc_im, matrix_hc_im(local_channel)%matrix, name="ot_kp_density_hc_im")
               CALL copy_fm_to_dbcsr(hc_re, density_hc_re)
               CALL copy_fm_to_dbcsr(hc_im, density_hc_im)
               CALL qs_scf_loop_project_ot_kp_final_energies(local_ot_env(local_channel), &
                                                             density_hc_re, density_hc_im, &
                                                             scf_control%smear)
               CALL dbcsr_release(density_hc_re)
               CALL dbcsr_release(density_hc_im)
               CALL cp_fm_release(hc_re)
               CALL cp_fm_release(hc_im)
            END IF
            IF (local_ot_env(local_channel)%settings%do_ener) THEN
               eigenvalues(:) = local_ot_env(local_channel)%ener_x(:)
               kp%mos(2, ispin)%eigenvalues(:) = local_ot_env(local_channel)%ener_x(:)
            END IF
         END DO
      END DO

      CALL kpoint_set_mo_occupation(kpoints, scf_control%smear, &
                                    added_mos_auto=scf_control%added_mos_auto, &
                                    added_mos_auto_grow=added_mos_auto_grow)
      IF (.NOT. added_mos_auto_grow) THEN
         CALL qs_scf_loop_guard_ot_kp_fermi_window(kpoints, local_ot_env, first_channel, &
                                                   nspin_ot, kp_range, wkp, para_env_inter_kp, &
                                                   scf_control%smear, fermi_window_guarded)
         IF (fermi_window_guarded) THEN
            CALL kpoint_set_mo_occupation(kpoints, scf_control%smear, &
                                          added_mos_auto=scf_control%added_mos_auto, &
                                          added_mos_auto_grow=added_mos_auto_grow)
         END IF
      END IF
      CALL qs_scf_loop_print_ot_kp_smearing(kpoints, local_ot_env, wkp, para_env_inter_kp)
      IF (.NOT. added_mos_auto_grow) THEN
         fmwork => scf_env%scf_work1
         IF (ASSOCIATED(scf_env%p_mix_new)) THEN
            density_target => scf_env%p_mix_new
         ELSE
            CALL qs_rho_get(rho, rho_ao_kp=rho_ao_kp)
            density_target => rho_ao_kp
         END IF
         CALL kpoint_density_matrices(kpoints)
         IF (kpoint_ot_project_energy_subproblem(local_ot_env(1), scf_control%smear)) THEN
            CALL qs_scf_loop_apply_ot_kp_projected_density(kpoints, local_ot_env, first_channel, &
                                                           nspin_ot, kp_range, scf_control%smear)
         END IF
         CALL kpoint_density_transform(kpoints, density_target, .FALSE., &
                                       matrix_s(1, 1)%matrix, sab_nl, fmwork, &
                                       qs_ot_env=scf_env%qs_ot_env)
         scf_env%iter_method = scf_env%qs_ot_env(first_channel)%OT_METHOD_FULL
         scf_env%iter_param = scf_env%qs_ot_env(first_channel)%ds_min
         qs_env%broyden_adaptive_sigma = scf_env%qs_ot_env(first_channel)%broyden_adaptive_sigma
      END IF

      IF (ALLOCATED(trial_eigenvalues)) DEALLOCATE (trial_eigenvalues)
      IF (ALLOCATED(trial_occupation)) DEALLOCATE (trial_occupation)
      IF (ALLOCATED(trial_mu)) DEALLOCATE (trial_mu)

      DO local_channel = 1, nlocal_channels
         CALL dbcsr_release(matrix_hc(local_channel)%matrix)
         CALL dbcsr_release(matrix_hc_im(local_channel)%matrix)
         DEALLOCATE (matrix_hc(local_channel)%matrix)
         DEALLOCATE (matrix_hc_im(local_channel)%matrix)
      END DO
      DEALLOCATE (matrix_hc, matrix_hc_im)

      CALL dbcsr_deallocate_matrix(matrix_k_re)
      CALL dbcsr_deallocate_matrix(matrix_k_im)
      CALL dbcsr_deallocate_matrix(matrix_s_re)
      CALL dbcsr_deallocate_matrix(matrix_s_im)

      CALL timestop(handle)

   END SUBROUTINE qs_scf_loop_do_ot_kp

! **************************************************************************************************
!> \brief collect the fixed-N trial occupations used by k-point OT
!> \param kpoints ...
!> \param local_ot_env ...
!> \param first_channel ...
!> \param nspin_ot ...
!> \param kp_range ...
!> \param wkp ...
!> \param para_env_inter_kp ...
!> \param smear ...
!> \param trial_eigenvalues ...
!> \param trial_occupation ...
!> \param trial_mu ...
! **************************************************************************************************
   SUBROUTINE qs_scf_loop_collect_ot_kp_trial_smearing(kpoints, local_ot_env, first_channel, nspin_ot, &
                                                       kp_range, wkp, para_env_inter_kp, smear, &
                                                       trial_eigenvalues, trial_occupation, trial_mu)

      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: local_ot_env
      INTEGER, INTENT(IN)                                :: first_channel, nspin_ot
      INTEGER, DIMENSION(2), INTENT(IN)                  :: kp_range
      REAL(KIND=dp), DIMENSION(:), POINTER               :: wkp
      TYPE(mp_para_env_type), POINTER                    :: para_env_inter_kp
      TYPE(smear_type), INTENT(IN)                       :: smear
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), INTENT(OUT) :: trial_eigenvalues, &
                                                            trial_occupation
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: trial_mu

      INTEGER                                            :: ikpoint, ispin, local_channel, &
                                                            local_kpoint, method, ne_a, ne_b, &
                                                            nelectron, nkpoint, nmo
      REAL(KIND=dp)                                      :: kTS, mu, nel, sigma
      REAL(KIND=dp), DIMENSION(2)                        :: mus
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvalues
      TYPE(kpoint_env_type), POINTER                     :: kp

      NULLIFY (eigenvalues)

      CALL get_kpoint_info(kpoints, nkp=nkpoint)
      kp => kpoints%kp_env(1)%kpoint_env
      CALL get_mo_set(kp%mos(1, 1), nmo=nmo, nelectron=nelectron)
      ne_a = nelectron
      ne_b = 0
      IF (nspin_ot == 2) CALL get_mo_set(kp%mos(1, 2), nelectron=ne_b)

      ALLOCATE (trial_eigenvalues(nmo, nkpoint, nspin_ot))
      ALLOCATE (trial_occupation(nmo, nkpoint, nspin_ot))
      ALLOCATE (trial_mu(nspin_ot))
      trial_eigenvalues(:, :, :) = 0.0_dp
      trial_occupation(:, :, :) = 0.0_dp
      trial_mu(:) = 0.0_dp
      mus(:) = 0.0_dp

      DO local_kpoint = 1, SIZE(kpoints%kp_env)
         kp => kpoints%kp_env(local_kpoint)%kpoint_env
         ikpoint = kp%nkpoint
         CPASSERT(ikpoint >= kp_range(1) .AND. ikpoint <= kp_range(2))
         DO ispin = 1, nspin_ot
            local_channel = qs_ot_channel_index(ispin, ikpoint, nspin_ot) - first_channel + 1
            CALL get_mo_set(kp%mos(1, ispin), eigenvalues=eigenvalues)
            IF (local_ot_env(local_channel)%settings%do_ener) THEN
               trial_eigenvalues(1:nmo, ikpoint, ispin) = local_ot_env(local_channel)%ener_x(1:nmo)
            ELSE
               trial_eigenvalues(1:nmo, ikpoint, ispin) = eigenvalues(1:nmo)
            END IF
         END DO
      END DO
      CALL para_env_inter_kp%sum(trial_eigenvalues)

      IF (smear%do_smear) THEN
         SELECT CASE (smear%method)
         CASE (smear_fermi_dirac)
            sigma = smear%electronic_temperature
            method = smear_fermi_dirac
         CASE (smear_gaussian, smear_mp, smear_mv)
            sigma = smear%smearing_width
            method = smear%method
         CASE DEFAULT
            CPABORT("kpoints: Selected smearing not (yet) supported")
         END SELECT
      ELSE
         sigma = 0.0_dp
         method = smear_gaussian
      END IF

      IF (nspin_ot == 1) THEN
         nel = REAL(nelectron, KIND=dp)
         CALL Smearkp(trial_occupation(:, :, 1), mus(1), kTS, trial_eigenvalues(:, :, 1), nel, wkp, &
                      sigma, 2.0_dp, method)
      ELSE IF (smear%fixed_mag_mom > 0.0_dp .OR. .NOT. smear%do_smear) THEN
         nel = REAL(ne_a, KIND=dp)
         CALL Smearkp(trial_occupation(:, :, 1), mus(1), kTS, trial_eigenvalues(:, :, 1), nel, wkp, &
                      sigma, 1.0_dp, method)
         nel = REAL(ne_b, KIND=dp)
         CALL Smearkp(trial_occupation(:, :, 2), mus(2), kTS, trial_eigenvalues(:, :, 2), nel, wkp, &
                      sigma, 1.0_dp, method)
      ELSE
         nel = REAL(ne_a, KIND=dp) + REAL(ne_b, KIND=dp)
         CALL Smearkp2(trial_occupation(:, :, :), mu, kTS, trial_eigenvalues(:, :, :), nel, wkp, &
                       sigma, method)
         mus(1:2) = mu
      END IF
      trial_mu(:) = mus(1:nspin_ot)

   END SUBROUTINE qs_scf_loop_collect_ot_kp_trial_smearing

! **************************************************************************************************
!> \brief keep trial k-point OT auxiliary energies inside a finite-temperature frontier window
!> \param kpoints ...
!> \param local_ot_env ...
!> \param first_channel ...
!> \param nspin_ot ...
!> \param kp_range ...
!> \param wkp ...
!> \param para_env_inter_kp ...
!> \param smear ...
!> \param trial_eigenvalues ...
!> \param trial_occupation ...
!> \param trial_mu ...
!> \param guarded ...
! **************************************************************************************************
   SUBROUTINE qs_scf_loop_guard_ot_kp_trial_fermi_window(kpoints, local_ot_env, first_channel, &
                                                         nspin_ot, kp_range, wkp, para_env_inter_kp, &
                                                         smear, trial_eigenvalues, trial_occupation, &
                                                         trial_mu, guarded)

      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: local_ot_env
      INTEGER, INTENT(IN)                                :: first_channel, nspin_ot
      INTEGER, DIMENSION(2), INTENT(IN)                  :: kp_range
      REAL(KIND=dp), DIMENSION(:), POINTER               :: wkp
      TYPE(mp_para_env_type), POINTER                    :: para_env_inter_kp
      TYPE(smear_type), INTENT(IN)                       :: smear
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: trial_eigenvalues
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: trial_occupation
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: trial_mu
      LOGICAL, INTENT(OUT)                               :: guarded

      INTEGER                                            :: adjusted, group, ikpoint, imo, ispin, &
                                                            local_channel, local_kpoint, ngroups, &
                                                            nkpoint, nmo
      LOGICAL                                            :: guard_group, lower_frontier, upper_frontier
      REAL(KIND=dp)                                      :: anchor_window, distance_tol, eps, maxocc, &
                                                            occ_tail, outer_window, sigma
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: frac_weight, lower_dist, max_partial, &
                                                            min_abs_eps_mu, upper_dist
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(kpoint_env_type), POINTER                     :: kp

      guarded = .FALSE.
      IF (.NOT. local_ot_env(1)%settings%do_ener) RETURN
      IF (.NOT. local_ot_env(1)%settings%occupation_preconditioner) RETURN
      IF (.NOT. smear%do_smear .OR. smear%method /= smear_fermi_dirac) RETURN

      sigma = smear%electronic_temperature
      IF (sigma <= 1.0E-14_dp) RETURN
      anchor_window = kpoint_ot_fermi_window_anchor(sigma)
      outer_window = ot_kp_fermi_window_outer_factor*sigma

      nmo = SIZE(trial_eigenvalues, 1)
      nkpoint = SIZE(trial_eigenvalues, 2)
      ngroups = 1
      IF (nspin_ot == 2 .AND. smear%fixed_mag_mom > 0.0_dp) ngroups = 2
      ALLOCATE (frac_weight(ngroups), lower_dist(ngroups), max_partial(ngroups), &
                min_abs_eps_mu(ngroups), upper_dist(ngroups))
      frac_weight(:) = 0.0_dp
      lower_dist(:) = HUGE(1.0_dp)
      max_partial(:) = 0.0_dp
      min_abs_eps_mu(:) = HUGE(1.0_dp)
      upper_dist(:) = HUGE(1.0_dp)

      kp => kpoints%kp_env(1)%kpoint_env
      DO ikpoint = 1, nkpoint
         DO ispin = 1, nspin_ot
            group = MERGE(ispin, 1, ngroups == nspin_ot)
            CALL get_mo_set(kp%mos(1, ispin), maxocc=maxocc)
            DO imo = 1, nmo
               occ_tail = MAX(0.0_dp, MIN(trial_occupation(imo, ikpoint, ispin), &
                                          maxocc - trial_occupation(imo, ikpoint, ispin)))
               frac_weight(group) = frac_weight(group) + wkp(ikpoint)*occ_tail
               max_partial(group) = MAX(max_partial(group), occ_tail)
               min_abs_eps_mu(group) = MIN(min_abs_eps_mu(group), &
                                           ABS(trial_eigenvalues(imo, ikpoint, ispin) - trial_mu(ispin)))
               IF (trial_eigenvalues(imo, ikpoint, ispin) <= trial_mu(ispin)) &
                  lower_dist(group) = MIN(lower_dist(group), &
                                          trial_mu(ispin) - trial_eigenvalues(imo, ikpoint, ispin))
               IF (trial_eigenvalues(imo, ikpoint, ispin) >= trial_mu(ispin)) &
                  upper_dist(group) = MIN(upper_dist(group), &
                                          trial_eigenvalues(imo, ikpoint, ispin) - trial_mu(ispin))
            END DO
         END DO
      END DO

      adjusted = 0
      DO local_kpoint = 1, SIZE(kpoints%kp_env)
         kp => kpoints%kp_env(local_kpoint)%kpoint_env
         ikpoint = kp%nkpoint
         CPASSERT(ikpoint >= kp_range(1) .AND. ikpoint <= kp_range(2))
         DO ispin = 1, nspin_ot
            group = MERGE(ispin, 1, ngroups == nspin_ot)
            guard_group = lower_dist(group) < HUGE(1.0_dp) .AND. &
                          upper_dist(group) < HUGE(1.0_dp) .AND. &
                          min_abs_eps_mu(group) > anchor_window .AND. &
                          (min_abs_eps_mu(group) < outer_window .OR. &
                           frac_weight(group) < ot_kp_fermi_window_min_frac)
            IF (.NOT. guard_group) CYCLE

            local_channel = qs_ot_channel_index(ispin, ikpoint, nspin_ot) - first_channel + 1
            distance_tol = MAX(1.0E-10_dp, 1.0E-8_dp*MAX(lower_dist(group), upper_dist(group)))
            DO imo = 1, nmo
               lower_frontier = trial_eigenvalues(imo, ikpoint, ispin) <= trial_mu(ispin) .AND. &
                                ABS(trial_mu(ispin) - trial_eigenvalues(imo, ikpoint, ispin) - &
                                    lower_dist(group)) <= distance_tol
               upper_frontier = trial_eigenvalues(imo, ikpoint, ispin) >= trial_mu(ispin) .AND. &
                                ABS(trial_eigenvalues(imo, ikpoint, ispin) - trial_mu(ispin) - &
                                    upper_dist(group)) <= distance_tol
               IF (lower_frontier .AND. lower_dist(group) > anchor_window) THEN
                  eps = trial_mu(ispin) - anchor_window
               ELSEIF (upper_frontier .AND. upper_dist(group) > anchor_window) THEN
                  eps = trial_mu(ispin) + anchor_window
               ELSE
                  CYCLE
               END IF
               trial_eigenvalues(imo, ikpoint, ispin) = eps
               local_ot_env(local_channel)%ener_x(imo) = eps
               adjusted = adjusted + 1
            END DO
         END DO
      END DO

      CALL para_env_inter_kp%sum(adjusted)
      guarded = adjusted > 0

      logger => cp_get_default_logger()
      IF (guarded .AND. logger%iter_info%print_level >= high_print_level .AND. &
          logger%para_env%is_source()) THEN
         WRITE (cp_logger_get_default_unit_nr(logger), '(A,4(1X,ES16.8),1X,I5)') &
            " OT k-point trial Fermi-window guard sigma anchor frac min_abs adjusted:", &
            sigma, anchor_window, MINVAL(frac_weight), MINVAL(min_abs_eps_mu), adjusted
      END IF

      DEALLOCATE (frac_weight, lower_dist, max_partial, min_abs_eps_mu, upper_dist)

   END SUBROUTINE qs_scf_loop_guard_ot_kp_trial_fermi_window

! **************************************************************************************************
!> \brief whether k-point OT should eliminate auxiliary energies by exact fixed-N Mermin projection
!> \param qs_ot_env ...
!> \param smear ...
! **************************************************************************************************
   FUNCTION kpoint_ot_project_energy_subproblem(qs_ot_env, smear) RESULT(do_projection)
      TYPE(qs_ot_type), INTENT(IN)                       :: qs_ot_env
      TYPE(smear_type), INTENT(IN)                       :: smear
      LOGICAL                                            :: do_projection

      do_projection = qs_ot_env%settings%do_ener .AND. &
                      qs_ot_env%settings%do_rotation .AND. &
                      smear%do_smear .AND. smear%method == smear_fermi_dirac

   END FUNCTION kpoint_ot_project_energy_subproblem

! **************************************************************************************************
!> \brief Fermi-window guard anchor from a target FD frontier tail occupation
!> \param sigma ...
! **************************************************************************************************
   FUNCTION kpoint_ot_fermi_window_anchor(sigma) RESULT(anchor)
      REAL(KIND=dp), INTENT(IN)                          :: sigma
      REAL(KIND=dp)                                      :: anchor

      REAL(KIND=dp)                                      :: tail_fraction

      tail_fraction = MIN(0.5_dp - EPSILON(1.0_dp), &
                          MAX(EPSILON(1.0_dp), ot_kp_fermi_window_tail_fraction))
      anchor = sigma*LOG((1.0_dp - tail_fraction)/tail_fraction)

   END FUNCTION kpoint_ot_fermi_window_anchor

! **************************************************************************************************
!> \brief residual-gated canonical projector blend weight for k-point OT
!> \param delta ...
!> \param max_weight ...
! **************************************************************************************************
   FUNCTION kpoint_ot_canonical_projector_weight(delta, max_weight) RESULT(projector_weight)
      REAL(KIND=dp), INTENT(IN)                          :: delta, max_weight
      REAL(KIND=dp)                                      :: projector_weight

      REAL(KIND=dp)                                      :: delta_scale

      projector_weight = 0.0_dp
      IF (delta > 0.0_dp .AND. max_weight > 0.0_dp) THEN
         delta_scale = ot_kp_canonical_projector_delta/MAX(delta, ot_kp_canonical_projector_delta)
         projector_weight = max_weight*delta_scale**2
      END IF

   END FUNCTION kpoint_ot_canonical_projector_weight

! **************************************************************************************************
!> \brief project k-point OT auxiliary energies onto the exact fixed-N Mermin subproblem
!> \param qs_ot_env ...
!> \param matrix_hc ...
!> \param matrix_hc_im ...
!> \param smear ...
! **************************************************************************************************
   SUBROUTINE qs_scf_loop_project_ot_kp_energy_subproblem(qs_ot_env, matrix_hc, matrix_hc_im, smear)

      TYPE(qs_ot_type)                                   :: qs_ot_env
      TYPE(dbcsr_type), POINTER                          :: matrix_hc, matrix_hc_im
      TYPE(smear_type), INTENT(IN)                       :: smear

      INTEGER                                            :: nmo
      REAL(KIND=dp)                                      :: gauge_shift
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: canonical_energies
      TYPE(dbcsr_type), POINTER                          :: b_im, b_re, c_im, c_re, chc_im, chc_re, &
                                                            tmp
      TYPE(dbcsr_type)                                   :: evec_im, evec_re

      IF (.NOT. kpoint_ot_project_energy_subproblem(qs_ot_env, smear)) RETURN

      CPASSERT(qs_ot_env%has_complex_kpoint_state)
      CPASSERT(ASSOCIATED(matrix_hc))
      CPASSERT(ASSOCIATED(matrix_hc_im))
      CPASSERT(ASSOCIATED(qs_ot_env%rot_mat_chc))
      CPASSERT(ASSOCIATED(qs_ot_env%rot_mat_chc_im))
      CPASSERT(ASSOCIATED(qs_ot_env%rot_mat_u_im))

      c_re => qs_ot_env%matrix_x
      c_im => qs_ot_env%matrix_x_im
      chc_re => qs_ot_env%rot_mat_chc
      chc_im => qs_ot_env%rot_mat_chc_im
      b_re => qs_ot_env%matrix_buf1
      b_im => qs_ot_env%buf2_k_k_nosym
      tmp => qs_ot_env%buf3_k_k_nosym

      CALL qs_ot_generate_rotation_complex(qs_ot_env)

      ! K = C_base(k)^H * H(k)C_current(k).
      CALL dbcsr_multiply('T', 'N', 1.0_dp, c_re, matrix_hc, 0.0_dp, chc_re)
      CALL dbcsr_multiply('T', 'N', 1.0_dp, c_im, matrix_hc_im, 0.0_dp, tmp)
      CALL dbcsr_add(chc_re, tmp, alpha_scalar=1.0_dp, beta_scalar=1.0_dp)

      CALL dbcsr_multiply('T', 'N', 1.0_dp, c_re, matrix_hc_im, 0.0_dp, chc_im)
      CALL dbcsr_multiply('T', 'N', 1.0_dp, c_im, matrix_hc, 0.0_dp, tmp)
      CALL dbcsr_add(chc_im, tmp, alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)

      ! CHC = C_current(k)^H * H(k)C_current(k) = U(k)^H * K.
      CALL dbcsr_multiply('T', 'N', 1.0_dp, qs_ot_env%rot_mat_u, chc_re, 0.0_dp, b_re)
      CALL dbcsr_multiply('T', 'N', 1.0_dp, qs_ot_env%rot_mat_u_im, chc_im, 0.0_dp, tmp)
      CALL dbcsr_add(b_re, tmp, alpha_scalar=1.0_dp, beta_scalar=1.0_dp)

      CALL dbcsr_multiply('T', 'N', 1.0_dp, qs_ot_env%rot_mat_u, chc_im, 0.0_dp, b_im)
      CALL dbcsr_multiply('T', 'N', 1.0_dp, qs_ot_env%rot_mat_u_im, chc_re, 0.0_dp, tmp)
      CALL dbcsr_add(b_im, tmp, alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)

      CALL dbcsr_copy(chc_re, b_re)
      CALL dbcsr_copy(chc_im, b_im)

      nmo = SIZE(qs_ot_env%ener_x)
      ALLOCATE (canonical_energies(nmo))
      CALL dbcsr_copy(evec_re, chc_re, name="ot_kp_chc_evec_re")
      CALL dbcsr_copy(evec_im, chc_im, name="ot_kp_chc_evec_im")
      CALL cp_dbcsr_heevd(matrix_re=chc_re, &
                          matrix_im=chc_im, &
                          eigenvectors_re=evec_re, &
                          eigenvectors_im=evec_im, &
                          eigenvalues=canonical_energies, &
                          para_env=qs_ot_env%para_env, &
                          blacs_env=qs_ot_env%blacs_env)

      gauge_shift = SUM(qs_ot_env%ener_x(:) - canonical_energies(:))/REAL(nmo, KIND=dp)
      qs_ot_env%ener_x(:) = canonical_energies(:) + gauge_shift
      CALL dbcsr_release(evec_re)
      CALL dbcsr_release(evec_im)
      DEALLOCATE (canonical_energies)

   END SUBROUTINE qs_scf_loop_project_ot_kp_energy_subproblem

! **************************************************************************************************
!> \brief project final k-point OT energies after the accepted rotation
!> \param qs_ot_env ...
!> \param matrix_hc ...
!> \param matrix_hc_im ...
!> \param smear ...
! **************************************************************************************************
   SUBROUTINE qs_scf_loop_project_ot_kp_final_energies(qs_ot_env, matrix_hc, matrix_hc_im, smear)

      TYPE(qs_ot_type)                                   :: qs_ot_env
      TYPE(dbcsr_type)                                   :: matrix_hc, matrix_hc_im
      TYPE(smear_type), INTENT(IN)                       :: smear

      INTEGER                                            :: nmo
      REAL(KIND=dp)                                      :: gauge_shift
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: canonical_energies
      TYPE(dbcsr_type), POINTER                          :: b_im, b_re, c_im, c_re, chc_im, chc_re, &
                                                            tmp
      TYPE(dbcsr_type)                                   :: evec_im, evec_re

      IF (.NOT. kpoint_ot_project_energy_subproblem(qs_ot_env, smear)) RETURN

      CPASSERT(qs_ot_env%has_complex_kpoint_state)
      CPASSERT(ASSOCIATED(qs_ot_env%rot_mat_chc))
      CPASSERT(ASSOCIATED(qs_ot_env%rot_mat_chc_im))
      CPASSERT(ASSOCIATED(qs_ot_env%rot_mat_u_im))

      c_re => qs_ot_env%matrix_x
      c_im => qs_ot_env%matrix_x_im
      chc_re => qs_ot_env%rot_mat_chc
      chc_im => qs_ot_env%rot_mat_chc_im
      b_re => qs_ot_env%matrix_buf1
      b_im => qs_ot_env%buf2_k_k_nosym
      tmp => qs_ot_env%buf3_k_k_nosym

      CALL qs_ot_generate_rotation_complex(qs_ot_env)

      ! K = C_base(k)^H * H(k)C_current(k).
      CALL dbcsr_multiply('T', 'N', 1.0_dp, c_re, matrix_hc, 0.0_dp, chc_re)
      CALL dbcsr_multiply('T', 'N', 1.0_dp, c_im, matrix_hc_im, 0.0_dp, tmp)
      CALL dbcsr_add(chc_re, tmp, alpha_scalar=1.0_dp, beta_scalar=1.0_dp)

      CALL dbcsr_multiply('T', 'N', 1.0_dp, c_re, matrix_hc_im, 0.0_dp, chc_im)
      CALL dbcsr_multiply('T', 'N', 1.0_dp, c_im, matrix_hc, 0.0_dp, tmp)
      CALL dbcsr_add(chc_im, tmp, alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)

      ! CHC = C_current(k)^H * H(k)C_current(k) = U(k)^H * K.
      CALL dbcsr_multiply('T', 'N', 1.0_dp, qs_ot_env%rot_mat_u, chc_re, 0.0_dp, b_re)
      CALL dbcsr_multiply('T', 'N', 1.0_dp, qs_ot_env%rot_mat_u_im, chc_im, 0.0_dp, tmp)
      CALL dbcsr_add(b_re, tmp, alpha_scalar=1.0_dp, beta_scalar=1.0_dp)

      CALL dbcsr_multiply('T', 'N', 1.0_dp, qs_ot_env%rot_mat_u, chc_im, 0.0_dp, b_im)
      CALL dbcsr_multiply('T', 'N', 1.0_dp, qs_ot_env%rot_mat_u_im, chc_re, 0.0_dp, tmp)
      CALL dbcsr_add(b_im, tmp, alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)

      CALL dbcsr_copy(chc_re, b_re)
      CALL dbcsr_copy(chc_im, b_im)

      nmo = SIZE(qs_ot_env%ener_x)
      ALLOCATE (canonical_energies(nmo))
      CALL dbcsr_copy(evec_re, chc_re, name="ot_kp_density_chc_evec_re")
      CALL dbcsr_copy(evec_im, chc_im, name="ot_kp_density_chc_evec_im")
      CALL cp_dbcsr_heevd(matrix_re=chc_re, &
                          matrix_im=chc_im, &
                          eigenvectors_re=evec_re, &
                          eigenvectors_im=evec_im, &
                          eigenvalues=canonical_energies, &
                          para_env=qs_ot_env%para_env, &
                          blacs_env=qs_ot_env%blacs_env)

      gauge_shift = SUM(qs_ot_env%ener_x(:) - canonical_energies(:))/REAL(nmo, KIND=dp)
      qs_ot_env%ener_x(:) = canonical_energies(:) + gauge_shift

      CALL dbcsr_release(evec_re)
      CALL dbcsr_release(evec_im)
      DEALLOCATE (canonical_energies)

   END SUBROUTINE qs_scf_loop_project_ot_kp_final_energies

! **************************************************************************************************
!> \brief damp the k-point OT density toward the canonical fixed-N Mermin projector
!> \param kpoints ...
!> \param local_ot_env ...
!> \param first_channel ...
!> \param nspin_ot ...
!> \param kp_range ...
!> \param smear ...
! **************************************************************************************************
   SUBROUTINE qs_scf_loop_apply_ot_kp_projected_density(kpoints, local_ot_env, first_channel, &
                                                        nspin_ot, kp_range, smear)

      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: local_ot_env
      INTEGER, INTENT(IN)                                :: first_channel, nspin_ot
      INTEGER, DIMENSION(2), INTENT(IN)                  :: kp_range
      TYPE(smear_type), INTENT(IN)                       :: smear

      INTEGER                                            :: ikpoint, ispin, local_channel, &
                                                            local_kpoint, nao, nblend, nmo
      REAL(KIND=dp)                                      :: density_weight, max_weight, sum_weight
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: canonical_energies
      REAL(KIND=dp), DIMENSION(:), POINTER               :: occupation_numbers
      TYPE(cp_fm_type)                                   :: can_im, can_re, evec_fm_im, &
                                                            evec_fm_re, work
      TYPE(cp_fm_type), POINTER                          :: cpmat, mo_coeff, mo_coeff_im, rpmat
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_type)                                   :: chc_im, chc_re, evec_im, evec_re
      TYPE(kpoint_env_type), POINTER                     :: kp

      logger => cp_get_default_logger()
      NULLIFY (cpmat, mo_coeff, mo_coeff_im, occupation_numbers, rpmat)
      max_weight = 0.0_dp
      sum_weight = 0.0_dp
      nblend = 0

      DO local_kpoint = 1, SIZE(kpoints%kp_env)
         kp => kpoints%kp_env(local_kpoint)%kpoint_env
         ikpoint = kp%nkpoint
         CPASSERT(ikpoint >= kp_range(1) .AND. ikpoint <= kp_range(2))
         DO ispin = 1, nspin_ot
            local_channel = qs_ot_channel_index(ispin, ikpoint, nspin_ot) - first_channel + 1
            density_weight = kpoint_ot_canonical_projector_weight(local_ot_env(local_channel)%delta, &
                                                                  ot_kp_canonical_density_max_weight)
            IF (density_weight <= 0.0_dp) CYCLE
            IF (.NOT. kpoint_ot_project_energy_subproblem(local_ot_env(local_channel), smear)) THEN
               CALL cp_abort(__LOCATION__, &
                             "Inconsistent k-point OT canonical density projector request.")
            END IF

            CALL get_mo_set(kp%mos(1, ispin), mo_coeff=mo_coeff, nmo=nmo, nao=nao, &
                            occupation_numbers=occupation_numbers)
            CALL get_mo_set(kp%mos(2, ispin), mo_coeff=mo_coeff_im)

            CALL dbcsr_copy(chc_re, local_ot_env(local_channel)%rot_mat_chc, &
                            name="ot_kp_density_chc_re")
            CALL dbcsr_copy(chc_im, local_ot_env(local_channel)%rot_mat_chc_im, &
                            name="ot_kp_density_chc_im")
            CALL dbcsr_copy(evec_re, chc_re, name="ot_kp_density_evec_re")
            CALL dbcsr_copy(evec_im, chc_im, name="ot_kp_density_evec_im")
            ALLOCATE (canonical_energies(nmo))
            CALL cp_dbcsr_heevd(matrix_re=chc_re, &
                                matrix_im=chc_im, &
                                eigenvectors_re=evec_re, &
                                eigenvectors_im=evec_im, &
                                eigenvalues=canonical_energies, &
                                para_env=local_ot_env(local_channel)%para_env, &
                                blacs_env=local_ot_env(local_channel)%blacs_env)

            CALL cp_fm_create(evec_fm_re, mo_coeff%matrix_struct, nrow=nmo, ncol=nmo)
            CALL cp_fm_create(evec_fm_im, mo_coeff%matrix_struct, nrow=nmo, ncol=nmo)
            CALL cp_fm_create(can_re, mo_coeff%matrix_struct)
            CALL cp_fm_create(can_im, mo_coeff%matrix_struct)
            CALL cp_fm_create(work, mo_coeff%matrix_struct)
            CALL copy_dbcsr_to_fm(evec_re, evec_fm_re)
            CALL copy_dbcsr_to_fm(evec_im, evec_fm_im)

            ! The canonical projector enters only as a damped density blend.  The
            ! MO/OT variables stay in the REF gauge used by the minimizer history.
            CALL parallel_gemm("N", "N", nao, nmo, nmo, 1.0_dp, mo_coeff, evec_fm_re, 0.0_dp, can_re)
            CALL parallel_gemm("N", "N", nao, nmo, nmo, -1.0_dp, mo_coeff_im, evec_fm_im, 1.0_dp, can_re)
            CALL parallel_gemm("N", "N", nao, nmo, nmo, 1.0_dp, mo_coeff, evec_fm_im, 0.0_dp, can_im)
            CALL parallel_gemm("N", "N", nao, nmo, nmo, 1.0_dp, mo_coeff_im, evec_fm_re, 1.0_dp, can_im)

            rpmat => kp%pmat(1, ispin)
            cpmat => kp%pmat(2, ispin)
            CALL cp_fm_to_fm(can_re, work)
            CALL cp_fm_column_scale(work, occupation_numbers)
            CALL parallel_gemm("N", "T", nao, nao, nmo, density_weight, can_re, work, &
                               1.0_dp - density_weight, rpmat)
            CALL parallel_gemm("N", "T", nao, nao, nmo, density_weight, can_im, work, &
                               1.0_dp - density_weight, cpmat)
            CALL parallel_gemm("N", "T", nao, nao, nmo, -density_weight, work, can_im, 1.0_dp, cpmat)

            CALL cp_fm_to_fm(can_im, work)
            CALL cp_fm_column_scale(work, occupation_numbers)
            CALL parallel_gemm("N", "T", nao, nao, nmo, density_weight, can_im, work, 1.0_dp, rpmat)

            max_weight = MAX(max_weight, density_weight)
            sum_weight = sum_weight + density_weight
            nblend = nblend + 1

            CALL cp_fm_release(evec_fm_re)
            CALL cp_fm_release(evec_fm_im)
            CALL cp_fm_release(can_re)
            CALL cp_fm_release(can_im)
            CALL cp_fm_release(work)
            CALL dbcsr_release(chc_re)
            CALL dbcsr_release(chc_im)
            CALL dbcsr_release(evec_re)
            CALL dbcsr_release(evec_im)
            DEALLOCATE (canonical_energies)
         END DO
      END DO

      IF (nblend > 0 .AND. logger%iter_info%print_level >= high_print_level .AND. &
          logger%para_env%is_source()) THEN
         WRITE (cp_logger_get_default_unit_nr(logger), '(A,2(1X,ES16.8),1X,I5,1X,I5)') &
            " OT canonical density projector max avg channels ls:", &
            max_weight, sum_weight/REAL(nblend, KIND=dp), nblend, &
            local_ot_env(1)%line_search_count
      END IF

   END SUBROUTINE qs_scf_loop_apply_ot_kp_projected_density

! **************************************************************************************************
!> \brief keep k-point OT auxiliary energies inside a finite-temperature frontier window
!> \param kpoints ...
!> \param local_ot_env ...
!> \param first_channel ...
!> \param nspin_ot ...
!> \param kp_range ...
!> \param wkp ...
!> \param para_env_inter_kp ...
!> \param smear ...
!> \param guarded ...
! **************************************************************************************************
   SUBROUTINE qs_scf_loop_guard_ot_kp_fermi_window(kpoints, local_ot_env, first_channel, &
                                                   nspin_ot, kp_range, wkp, para_env_inter_kp, &
                                                   smear, guarded)

      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: local_ot_env
      INTEGER, INTENT(IN)                                :: first_channel, nspin_ot
      INTEGER, DIMENSION(2), INTENT(IN)                  :: kp_range
      REAL(KIND=dp), DIMENSION(:), POINTER               :: wkp
      TYPE(mp_para_env_type), POINTER                    :: para_env_inter_kp
      TYPE(smear_type), INTENT(IN)                       :: smear
      LOGICAL, INTENT(OUT)                               :: guarded

      INTEGER                                            :: adjusted, group, ikpoint, imo, ispin, &
                                                            local_channel, local_kpoint, ngroups, nmo
      LOGICAL                                            :: guard_group, lower_frontier, upper_frontier
      REAL(KIND=dp)                                      :: anchor_window, distance_tol, eps, maxocc, &
                                                            occ_tail, outer_window, sigma
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: frac_weight, lower_dist, max_partial, &
                                                            min_abs_eps_mu, mu_count, mu_group, &
                                                            mu_sum, upper_dist
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvalues, occupation_numbers
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(kpoint_env_type), POINTER                     :: kp

      guarded = .FALSE.
      IF (.NOT. local_ot_env(1)%settings%do_ener) RETURN
      IF (.NOT. local_ot_env(1)%settings%occupation_preconditioner) RETURN
      IF (.NOT. smear%do_smear .OR. smear%method /= smear_fermi_dirac) RETURN

      sigma = smear%electronic_temperature
      IF (sigma <= 1.0E-14_dp) RETURN
      anchor_window = kpoint_ot_fermi_window_anchor(sigma)
      outer_window = ot_kp_fermi_window_outer_factor*sigma

      ngroups = 1
      IF (nspin_ot == 2 .AND. smear%fixed_mag_mom > 0.0_dp) ngroups = 2
      ALLOCATE (frac_weight(ngroups), lower_dist(ngroups), max_partial(ngroups), &
                min_abs_eps_mu(ngroups), mu_count(ngroups), mu_group(ngroups), &
                mu_sum(ngroups), upper_dist(ngroups))
      frac_weight(:) = 0.0_dp
      lower_dist(:) = HUGE(1.0_dp)
      max_partial(:) = 0.0_dp
      min_abs_eps_mu(:) = HUGE(1.0_dp)
      mu_count(:) = 0.0_dp
      mu_group(:) = 0.0_dp
      mu_sum(:) = 0.0_dp
      upper_dist(:) = HUGE(1.0_dp)

      DO local_kpoint = 1, SIZE(kpoints%kp_env)
         kp => kpoints%kp_env(local_kpoint)%kpoint_env
         ikpoint = kp%nkpoint
         CPASSERT(ikpoint >= kp_range(1) .AND. ikpoint <= kp_range(2))
         DO ispin = 1, nspin_ot
            group = MERGE(ispin, 1, ngroups == nspin_ot)
            mu_sum(group) = mu_sum(group) + kp%mos(1, ispin)%mu
            mu_count(group) = mu_count(group) + 1.0_dp
         END DO
      END DO

      CALL para_env_inter_kp%sum(mu_sum)
      CALL para_env_inter_kp%sum(mu_count)
      DO group = 1, ngroups
         IF (mu_count(group) > 0.0_dp) mu_group(group) = mu_sum(group)/mu_count(group)
      END DO

      DO local_kpoint = 1, SIZE(kpoints%kp_env)
         kp => kpoints%kp_env(local_kpoint)%kpoint_env
         ikpoint = kp%nkpoint
         DO ispin = 1, nspin_ot
            group = MERGE(ispin, 1, ngroups == nspin_ot)
            CALL get_mo_set(kp%mos(1, ispin), maxocc=maxocc, nmo=nmo, eigenvalues=eigenvalues, &
                            occupation_numbers=occupation_numbers)
            DO imo = 1, nmo
               occ_tail = MAX(0.0_dp, MIN(occupation_numbers(imo), maxocc - occupation_numbers(imo)))
               frac_weight(group) = frac_weight(group) + wkp(ikpoint)*occ_tail
               max_partial(group) = MAX(max_partial(group), occ_tail)
               min_abs_eps_mu(group) = MIN(min_abs_eps_mu(group), ABS(eigenvalues(imo) - mu_group(group)))
               IF (eigenvalues(imo) <= mu_group(group)) &
                  lower_dist(group) = MIN(lower_dist(group), mu_group(group) - eigenvalues(imo))
               IF (eigenvalues(imo) >= mu_group(group)) &
                  upper_dist(group) = MIN(upper_dist(group), eigenvalues(imo) - mu_group(group))
            END DO
         END DO
      END DO

      DO group = 1, ngroups
         CALL para_env_inter_kp%sum(frac_weight(group))
         CALL para_env_inter_kp%max(max_partial(group))
         CALL para_env_inter_kp%min(min_abs_eps_mu(group))
         CALL para_env_inter_kp%min(lower_dist(group))
         CALL para_env_inter_kp%min(upper_dist(group))
      END DO

      adjusted = 0
      DO local_kpoint = 1, SIZE(kpoints%kp_env)
         kp => kpoints%kp_env(local_kpoint)%kpoint_env
         ikpoint = kp%nkpoint
         DO ispin = 1, nspin_ot
            group = MERGE(ispin, 1, ngroups == nspin_ot)
            guard_group = lower_dist(group) < HUGE(1.0_dp) .AND. &
                          upper_dist(group) < HUGE(1.0_dp) .AND. &
                          min_abs_eps_mu(group) > anchor_window .AND. &
                          (min_abs_eps_mu(group) < outer_window .OR. &
                           frac_weight(group) < ot_kp_fermi_window_min_frac)
            IF (.NOT. guard_group) CYCLE

            local_channel = qs_ot_channel_index(ispin, ikpoint, nspin_ot) - first_channel + 1
            CALL get_mo_set(kp%mos(1, ispin), nmo=nmo, eigenvalues=eigenvalues)
            distance_tol = MAX(1.0E-10_dp, 1.0E-8_dp*MAX(lower_dist(group), upper_dist(group)))
            DO imo = 1, nmo
               lower_frontier = eigenvalues(imo) <= mu_group(group) .AND. &
                                ABS(mu_group(group) - eigenvalues(imo) - lower_dist(group)) <= distance_tol
               upper_frontier = eigenvalues(imo) >= mu_group(group) .AND. &
                                ABS(eigenvalues(imo) - mu_group(group) - upper_dist(group)) <= distance_tol
               IF (lower_frontier .AND. lower_dist(group) > anchor_window) THEN
                  eps = mu_group(group) - anchor_window
               ELSEIF (upper_frontier .AND. upper_dist(group) > anchor_window) THEN
                  eps = mu_group(group) + anchor_window
               ELSE
                  CYCLE
               END IF
               eigenvalues(imo) = eps
               kp%mos(2, ispin)%eigenvalues(imo) = eps
               local_ot_env(local_channel)%ener_x(imo) = eps
               adjusted = adjusted + 1
            END DO
         END DO
      END DO

      CALL para_env_inter_kp%sum(adjusted)
      guarded = adjusted > 0

      logger => cp_get_default_logger()
      IF (guarded .AND. logger%iter_info%print_level >= high_print_level .AND. &
          logger%para_env%is_source()) THEN
         WRITE (cp_logger_get_default_unit_nr(logger), '(A,4(1X,ES16.8),1X,I5)') &
            " OT k-point Fermi-window guard sigma anchor frac min_abs adjusted:", &
            sigma, anchor_window, MINVAL(frac_weight), MINVAL(min_abs_eps_mu), adjusted
      END IF

      DEALLOCATE (frac_weight, lower_dist, max_partial, min_abs_eps_mu, mu_count, mu_group, &
                  mu_sum, upper_dist)

   END SUBROUTINE qs_scf_loop_guard_ot_kp_fermi_window

! **************************************************************************************************
!> \brief high-print diagnostic for k-point OT smearing quality
!> \param kpoints ...
!> \param wkp ...
!> \param para_env_inter_kp ...
! **************************************************************************************************
   SUBROUTINE qs_scf_loop_print_ot_kp_smearing(kpoints, local_ot_env, wkp, para_env_inter_kp)

      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: local_ot_env
      REAL(KIND=dp), DIMENSION(:), POINTER               :: wkp
      TYPE(mp_para_env_type), POINTER                    :: para_env_inter_kp

      INTEGER                                            :: ikpoint, imo, ispin, local_channel, &
                                                            local_kpoint, nmo, nspin
      REAL(KIND=dp)                                      :: density_projector_fraction, &
                                                            density_projector_weight, fractional_weight, &
                                                            kTS, max_partial, maxocc, min_fermi_dist, &
                                                            min_trust_scale, mu, occ_tail, residual_delta, &
                                                            tether_energy, tether_grad_max, trusted_step_max, &
                                                            weighted_electrons
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvalues, occupation_numbers
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(kpoint_env_type), POINTER                     :: kp

      logger => cp_get_default_logger()
      IF (logger%iter_info%print_level < high_print_level) RETURN

      fractional_weight = 0.0_dp
      weighted_electrons = 0.0_dp
      max_partial = 0.0_dp
      min_fermi_dist = HUGE(1.0_dp)
      kTS = 0.0_dp
      mu = 0.0_dp
      residual_delta = 0.0_dp
      tether_energy = 0.0_dp
      tether_grad_max = 0.0_dp
      trusted_step_max = 0.0_dp
      min_trust_scale = 1.0_dp

      kp => kpoints%kp_env(1)%kpoint_env
      nspin = SIZE(kp%mos, 2)
      DO ispin = 1, nspin
         kTS = kTS + kp%mos(1, ispin)%kTS
         mu = mu + kp%mos(1, ispin)%mu
      END DO
      mu = mu/REAL(nspin, KIND=dp)

      DO local_kpoint = 1, SIZE(kpoints%kp_env)
         kp => kpoints%kp_env(local_kpoint)%kpoint_env
         ikpoint = kp%nkpoint
         DO ispin = 1, nspin
            CALL get_mo_set(kp%mos(1, ispin), maxocc=maxocc, nmo=nmo, eigenvalues=eigenvalues, &
                            occupation_numbers=occupation_numbers)
            DO imo = 1, nmo
               occ_tail = MAX(0.0_dp, MIN(occupation_numbers(imo), maxocc - occupation_numbers(imo)))
               max_partial = MAX(max_partial, occ_tail)
               fractional_weight = fractional_weight + wkp(ikpoint)*occ_tail
               weighted_electrons = weighted_electrons + wkp(ikpoint)*occupation_numbers(imo)
               min_fermi_dist = MIN(min_fermi_dist, ABS(eigenvalues(imo) - kp%mos(1, ispin)%mu))
            END DO
         END DO
      END DO

      CALL para_env_inter_kp%sum(fractional_weight)
      CALL para_env_inter_kp%sum(weighted_electrons)
      CALL para_env_inter_kp%max(max_partial)
      CALL para_env_inter_kp%min(min_fermi_dist)

      DO local_channel = 1, SIZE(local_ot_env)
         residual_delta = MAX(residual_delta, local_ot_env(local_channel)%delta)
         tether_energy = tether_energy + local_ot_env(local_channel)%canonical_tether_energy
         tether_grad_max = MAX(tether_grad_max, local_ot_env(local_channel)%canonical_tether_grad_max)
         trusted_step_max = MAX(trusted_step_max, local_ot_env(local_channel)%ener_x_trusted_step_max)
         IF (local_ot_env(local_channel)%ener_x_trust_radius > 0.0_dp) &
            min_trust_scale = MIN(min_trust_scale, local_ot_env(local_channel)%ener_x_trust_scale)
      END DO
      CALL para_env_inter_kp%max(residual_delta)
      CALL para_env_inter_kp%sum(tether_energy)
      CALL para_env_inter_kp%max(tether_grad_max)
      CALL para_env_inter_kp%max(trusted_step_max)
      CALL para_env_inter_kp%min(min_trust_scale)

      density_projector_weight = kpoint_ot_canonical_projector_weight(residual_delta, &
                                                                      ot_kp_canonical_density_max_weight)
      density_projector_fraction = 0.0_dp
      IF (ot_kp_canonical_density_max_weight > 0.0_dp) &
         density_projector_fraction = density_projector_weight/ot_kp_canonical_density_max_weight

      IF (logger%para_env%is_source()) THEN
         WRITE (cp_logger_get_default_unit_nr(logger), '(A,6(1X,ES16.8))') &
            " OT k-point smearing kTS mu weighted_N frac_weight max_partial min_abs_eps_mu:", &
            kTS, mu, weighted_electrons, fractional_weight, max_partial, min_fermi_dist
         WRITE (cp_logger_get_default_unit_nr(logger), '(A,7(1X,ES16.8),1X,I5)') &
            " OT k-point residual model delta dproj frac trusted scale tether_E tether_g ls:", &
            residual_delta, density_projector_weight, density_projector_fraction, trusted_step_max, &
            min_trust_scale, tether_energy, tether_grad_max, local_ot_env(1)%line_search_count
      END IF

   END SUBROUTINE qs_scf_loop_print_ot_kp_smearing

! **************************************************************************************************
!> \brief prepare per-channel Fermi-window trust data for k-point OT energy variables
!> \param kpoints ...
!> \param local_ot_env ...
!> \param first_channel ...
!> \param nspin_ot ...
!> \param wkp ...
!> \param smear ...
!> \param trial_eigenvalues ...
!> \param trial_occupation ...
!> \param trial_mu ...
! **************************************************************************************************
   SUBROUTINE qs_scf_loop_prepare_ot_kp_energy_trust(kpoints, local_ot_env, first_channel, nspin_ot, &
                                                     wkp, smear, trial_eigenvalues, trial_occupation, &
                                                     trial_mu)

      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: local_ot_env
      INTEGER, INTENT(IN)                                :: first_channel, nspin_ot
      REAL(KIND=dp), DIMENSION(:), POINTER               :: wkp
      TYPE(smear_type), INTENT(IN)                       :: smear
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: trial_eigenvalues, trial_occupation
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: trial_mu

      REAL(KIND=dp), PARAMETER                           :: edge_window_factor = 8.0_dp, &
                                                            partial_tol = 1.0E-8_dp, &
                                                            trust_radius_factor = 4.0_dp

      INTEGER                                            :: ikpoint, imo, ispin, local_channel, &
                                                            local_kpoint, nmo
      LOGICAL                                            :: active, has_trust_weight
      REAL(KIND=dp)                                      :: edge_window, fi_norm, frac_weight, gvalue, &
                                                            maxocc, min_abs_eps_mu, mu, occ_tail, &
                                                            rayleigh_edge_gap, rayleigh_gap, sigma, &
                                                            trust_radius
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: rayleigh_diag
      TYPE(kpoint_env_type), POINTER                     :: kp

      DO local_channel = 1, SIZE(local_ot_env)
         local_ot_env(local_channel)%ener_x_trust_radius = -1.0_dp
         local_ot_env(local_channel)%ener_x_trust_edge_window = -1.0_dp
         local_ot_env(local_channel)%ener_x_trust_mu = 0.0_dp
         local_ot_env(local_channel)%ener_x_trust_sigma = 0.0_dp
         local_ot_env(local_channel)%ener_x_trust_frac_weight = 0.0_dp
         local_ot_env(local_channel)%ener_x_trust_entropy = 0.0_dp
         local_ot_env(local_channel)%ener_x_trust_min_abs_eps_mu = HUGE(1.0_dp)
         local_ot_env(local_channel)%ener_x_trust_rayleigh_gap = 0.0_dp
         local_ot_env(local_channel)%ener_x_trust_rayleigh_edge_gap = 0.0_dp
         local_ot_env(local_channel)%ener_x_step_max = 0.0_dp
         local_ot_env(local_channel)%ener_x_trusted_step_max = 0.0_dp
         local_ot_env(local_channel)%ener_x_trust_scale = 1.0_dp
         IF (ASSOCIATED(local_ot_env(local_channel)%ener_x_trust_weight)) &
            local_ot_env(local_channel)%ener_x_trust_weight(:) = 0.0_dp
      END DO

      IF (.NOT. local_ot_env(1)%settings%do_ener) RETURN
      IF (.NOT. local_ot_env(1)%settings%occupation_preconditioner) RETURN
      IF (.NOT. smear%do_smear .OR. smear%method /= smear_fermi_dirac) RETURN

      sigma = smear%electronic_temperature
      IF (sigma <= 1.0E-14_dp) RETURN
      trust_radius = trust_radius_factor*sigma
      edge_window = edge_window_factor*sigma

      DO local_kpoint = 1, SIZE(kpoints%kp_env)
         kp => kpoints%kp_env(local_kpoint)%kpoint_env
         ikpoint = kp%nkpoint
         DO ispin = 1, nspin_ot
            local_channel = qs_ot_channel_index(ispin, ikpoint, nspin_ot) - first_channel + 1
            CALL get_mo_set(kp%mos(1, ispin), maxocc=maxocc, nmo=nmo)
            mu = trial_mu(ispin)
            has_trust_weight = .FALSE.
            IF (ASSOCIATED(local_ot_env(local_channel)%ener_x_trust_weight)) &
               has_trust_weight = SIZE(local_ot_env(local_channel)%ener_x_trust_weight) >= nmo
            local_ot_env(local_channel)%ener_x_trust_radius = trust_radius
            local_ot_env(local_channel)%ener_x_trust_edge_window = edge_window
            local_ot_env(local_channel)%ener_x_trust_mu = mu
            local_ot_env(local_channel)%ener_x_trust_sigma = sigma
            local_ot_env(local_channel)%ener_x_trust_entropy = kp%mos(1, ispin)%kTS
            frac_weight = 0.0_dp
            min_abs_eps_mu = HUGE(1.0_dp)
            rayleigh_gap = 0.0_dp
            rayleigh_edge_gap = 0.0_dp
            ALLOCATE (rayleigh_diag(nmo))
            CALL dbcsr_get_diag(local_ot_env(local_channel)%rot_mat_chc, rayleigh_diag)

            DO imo = 1, nmo
               occ_tail = MAX(0.0_dp, MIN(trial_occupation(imo, ikpoint, ispin), &
                                          maxocc - trial_occupation(imo, ikpoint, ispin)))
               frac_weight = frac_weight + wkp(ikpoint)*occ_tail
               min_abs_eps_mu = MIN(min_abs_eps_mu, ABS(trial_eigenvalues(imo, ikpoint, ispin) - mu))
               rayleigh_gap = MAX(rayleigh_gap, ABS(rayleigh_diag(imo) - local_ot_env(local_channel)%ener_x(imo)))
               active = occ_tail > partial_tol*MAX(maxocc, 1.0_dp) .OR. &
                        ABS(trial_eigenvalues(imo, ikpoint, ispin) - mu) < edge_window
               IF (active) rayleigh_edge_gap = MAX(rayleigh_edge_gap, &
                                                   ABS(rayleigh_diag(imo) - local_ot_env(local_channel)%ener_x(imo)))
               IF (active .AND. has_trust_weight) THEN
                  fi_norm = trial_occupation(imo, ikpoint, ispin)/MAX(maxocc, EPSILON(maxocc))
                  gvalue = maxocc*fi_norm*(1.0_dp - fi_norm)/sigma
                  local_ot_env(local_channel)%ener_x_trust_weight(imo) = wkp(ikpoint)*MAX(gvalue, 0.0_dp)
               END IF
            END DO
            DEALLOCATE (rayleigh_diag)
            local_ot_env(local_channel)%ener_x_trust_frac_weight = frac_weight
            local_ot_env(local_channel)%ener_x_trust_min_abs_eps_mu = min_abs_eps_mu
            local_ot_env(local_channel)%ener_x_trust_rayleigh_gap = rayleigh_gap
            local_ot_env(local_channel)%ener_x_trust_rayleigh_edge_gap = rayleigh_edge_gap
         END DO
      END DO

   END SUBROUTINE qs_scf_loop_prepare_ot_kp_energy_trust

! **************************************************************************************************
!> \brief active Fermi-window response for the k-point OT canonical tether
!> \param occupation ...
!> \param eigenvalue ...
!> \param mu ...
!> \param sigma ...
!> \param maxocc ...
! **************************************************************************************************
   FUNCTION kpoint_ot_canonical_tether_response(occupation, eigenvalue, mu, sigma, maxocc) RESULT(response)
      REAL(KIND=dp), INTENT(IN)                          :: occupation, eigenvalue, mu, sigma, &
                                                            maxocc
      REAL(KIND=dp)                                      :: response

      REAL(KIND=dp), PARAMETER                           :: edge_window_factor = 8.0_dp, &
                                                            partial_tol = 1.0E-8_dp, &
                                                            response_floor_fraction = 1.0E-2_dp

      LOGICAL                                            :: active
      REAL(KIND=dp)                                      :: edge_window, fi_norm, gvalue, occ_tail, &
                                                            response_floor

      response = 0.0_dp
      IF (sigma <= 1.0E-14_dp .OR. maxocc <= EPSILON(maxocc)) RETURN

      edge_window = edge_window_factor*sigma
      occ_tail = MAX(0.0_dp, MIN(occupation, maxocc - occupation))
      active = occ_tail > partial_tol*MAX(maxocc, 1.0_dp) .OR. ABS(eigenvalue - mu) < edge_window
      IF (.NOT. active) RETURN

      fi_norm = occupation/MAX(maxocc, EPSILON(maxocc))
      gvalue = maxocc*fi_norm*(1.0_dp - fi_norm)/sigma
      response_floor = response_floor_fraction*maxocc/sigma
      response = MAX(gvalue, response_floor)

   END FUNCTION kpoint_ot_canonical_tether_response

! **************************************************************************************************
!> \brief adaptive prefactor for the k-point OT canonical tether
!> \param qs_ot_env ...
!> \param edge_gap ...
!> \param sigma ...
! **************************************************************************************************
   FUNCTION kpoint_ot_canonical_tether_strength(qs_ot_env, edge_gap, sigma) RESULT(strength)
      TYPE(qs_ot_type), INTENT(IN)                       :: qs_ot_env
      REAL(KIND=dp), INTENT(IN)                          :: edge_gap, sigma
      REAL(KIND=dp)                                      :: strength

      REAL(KIND=dp)                                      :: base_strength, gap_scale, max_strength

      strength = 0.0_dp
      IF (.NOT. qs_ot_env%settings%canonical_tether) RETURN

      base_strength = qs_ot_env%settings%canonical_tether_strength
      max_strength = qs_ot_env%settings%canonical_tether_max_strength
      IF (base_strength <= 0.0_dp .OR. max_strength <= 0.0_dp) RETURN

      strength = MIN(base_strength, max_strength)
      IF (.NOT. qs_ot_env%settings%canonical_tether_adaptive) RETURN

      ! Keep line-search trial points on the same model; adapt only when a fresh
      ! accepted-gradient point is building the next search direction.
      IF (qs_ot_env%line_search_count > 0) THEN
         IF (qs_ot_env%canonical_tether_strength > 0.0_dp) &
            strength = MIN(max_strength, qs_ot_env%canonical_tether_strength)
         RETURN
      END IF
      IF (sigma <= 1.0E-14_dp .OR. edge_gap <= 0.0_dp) RETURN

      gap_scale = SQRT(MAX(1.0_dp, edge_gap/MAX(16.0_dp*sigma, EPSILON(sigma))))
      strength = MIN(max_strength, base_strength*gap_scale)

   END FUNCTION kpoint_ot_canonical_tether_strength

! **************************************************************************************************
!> \brief response-weighted Fermi-window tether weight for k-point OT auxiliary energies
!> \param occupation ...
!> \param eigenvalue ...
!> \param mu ...
!> \param sigma ...
!> \param maxocc ...
!> \param kpoint_weight ...
!> \param tether_strength ...
! **************************************************************************************************
   FUNCTION kpoint_ot_canonical_tether_weight(occupation, eigenvalue, mu, sigma, maxocc, &
                                              kpoint_weight, tether_strength) RESULT(weight)
      REAL(KIND=dp), INTENT(IN)                          :: occupation, eigenvalue, mu, sigma, &
                                                            maxocc, kpoint_weight, tether_strength
      REAL(KIND=dp)                                      :: weight

      weight = tether_strength*kpoint_weight* &
               kpoint_ot_canonical_tether_response(occupation, eigenvalue, mu, sigma, maxocc)

   END FUNCTION kpoint_ot_canonical_tether_weight

! **************************************************************************************************
!> \brief prepare complex rotation and band-energy gradients for one local k-point OT channel
!> \param qs_ot_env ...
!> \param matrix_hc ...
!> \param matrix_hc_im ...
!> \param occupation_numbers ...
!> \param kpoint_weight ...
! **************************************************************************************************
   SUBROUTINE qs_scf_loop_prepare_ot_kp_channel_gradient(qs_ot_env, matrix_hc, matrix_hc_im, &
                                                         occupation_numbers, eigenvalues, maxocc, mu, &
                                                         kpoint_weight, smear, nondiag_energy)

      TYPE(qs_ot_type)                                   :: qs_ot_env
      TYPE(dbcsr_type), POINTER                          :: matrix_hc, matrix_hc_im
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: occupation_numbers
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: eigenvalues
      REAL(KIND=dp), INTENT(IN)                          :: maxocc, mu, kpoint_weight
      TYPE(smear_type), INTENT(IN)                       :: smear
      REAL(KIND=dp), INTENT(INOUT)                       :: nondiag_energy

      INTEGER                                            :: imo, nmo
      LOGICAL                                            :: project_energy_subproblem
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: canonical_energies, diagonal, &
                                                            orbital_scaling_factor, projector_diag, &
                                                            tether_diag
      REAL(KIND=dp)                                      :: canonical_projector_weight, canonical_tether_energy, &
                                                            delta, nondiag_trace_im, nondiag_trace_re, &
                                                            nondiag_weight, rayleigh_edge_gap, sigma, tether_response, &
                                                            tether_shift, tether_strength, tether_weight, &
                                                            tether_weight_sum
      TYPE(dbcsr_type)                                   :: canonical_evec_im, canonical_evec_re, &
                                                            diag_im, diag_re, occ_im, occ_re, &
                                                            weighted_evec_im, weighted_evec_re
      TYPE(dbcsr_type), POINTER                          :: b_im, b_re, c_im, c_re, chc_im, chc_re, tmp, &
                                                            tmp2

      IF (.NOT. qs_ot_env%settings%do_rotation) RETURN

      CPASSERT(qs_ot_env%has_complex_kpoint_state)
      CPASSERT(ASSOCIATED(matrix_hc))
      CPASSERT(ASSOCIATED(matrix_hc_im))
      CPASSERT(ASSOCIATED(qs_ot_env%rot_mat_chc))
      CPASSERT(ASSOCIATED(qs_ot_env%rot_mat_chc_im))
      CPASSERT(ASSOCIATED(qs_ot_env%rot_mat_u_im))

      c_re => qs_ot_env%matrix_x
      c_im => qs_ot_env%matrix_x_im
      chc_re => qs_ot_env%rot_mat_chc
      chc_im => qs_ot_env%rot_mat_chc_im
      b_re => qs_ot_env%matrix_buf1
      b_im => qs_ot_env%buf2_k_k_nosym
      tmp => qs_ot_env%buf3_k_k_nosym
      tmp2 => qs_ot_env%buf4_k_k_nosym
      project_energy_subproblem = kpoint_ot_project_energy_subproblem(qs_ot_env, smear)

      CALL qs_ot_generate_rotation_complex(qs_ot_env)

      ! K = C_base(k)^H * H(k)C_current(k).
      CALL dbcsr_multiply('T', 'N', 1.0_dp, c_re, matrix_hc, 0.0_dp, chc_re)
      CALL dbcsr_multiply('T', 'N', 1.0_dp, c_im, matrix_hc_im, 0.0_dp, tmp)
      CALL dbcsr_add(chc_re, tmp, alpha_scalar=1.0_dp, beta_scalar=1.0_dp)

      CALL dbcsr_multiply('T', 'N', 1.0_dp, c_re, matrix_hc_im, 0.0_dp, chc_im)
      CALL dbcsr_multiply('T', 'N', 1.0_dp, c_im, matrix_hc, 0.0_dp, tmp)
      CALL dbcsr_add(chc_im, tmp, alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)

      ! CHC = C_current(k)^H * H(k)C_current(k) = U(k)^H * K.
      CALL dbcsr_multiply('T', 'N', 1.0_dp, qs_ot_env%rot_mat_u, chc_re, 0.0_dp, b_re)
      CALL dbcsr_multiply('T', 'N', 1.0_dp, qs_ot_env%rot_mat_u_im, chc_im, 0.0_dp, tmp)
      CALL dbcsr_add(b_re, tmp, alpha_scalar=1.0_dp, beta_scalar=1.0_dp)

      CALL dbcsr_multiply('T', 'N', 1.0_dp, qs_ot_env%rot_mat_u, chc_im, 0.0_dp, b_im)
      CALL dbcsr_multiply('T', 'N', 1.0_dp, qs_ot_env%rot_mat_u_im, chc_re, 0.0_dp, tmp)
      CALL dbcsr_add(b_im, tmp, alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)

      CALL dbcsr_copy(chc_re, b_re)
      CALL dbcsr_copy(chc_im, b_im)

      IF (qs_ot_env%settings%do_ener) THEN
         CALL dbcsr_get_diag(chc_re, qs_ot_env%ener_gx)
         IF (project_energy_subproblem) &
            qs_ot_env%ener_gx(:) = eigenvalues(:)
      END IF

      IF (project_energy_subproblem) THEN
         nmo = SIZE(occupation_numbers)
         ALLOCATE (canonical_energies(nmo), orbital_scaling_factor(nmo))
         orbital_scaling_factor(:) = 2.0_dp*kpoint_weight*occupation_numbers(:)

         CALL dbcsr_copy(diag_re, chc_re, name="ot_kp_chc_diag_re")
         CALL dbcsr_copy(diag_im, chc_im, name="ot_kp_chc_diag_im")
         CALL dbcsr_copy(canonical_evec_re, chc_re, name="ot_kp_chc_canon_re")
         CALL dbcsr_copy(canonical_evec_im, chc_im, name="ot_kp_chc_canon_im")
         CALL cp_dbcsr_heevd(matrix_re=diag_re, &
                             matrix_im=diag_im, &
                             eigenvectors_re=canonical_evec_re, &
                             eigenvectors_im=canonical_evec_im, &
                             eigenvalues=canonical_energies, &
                             para_env=qs_ot_env%para_env, &
                             blacs_env=qs_ot_env%blacs_env)

         CALL dbcsr_copy(weighted_evec_re, canonical_evec_re, name="ot_kp_occ_evec_re")
         CALL dbcsr_scale_by_vector(weighted_evec_re, alpha=orbital_scaling_factor, side='right')
         CALL dbcsr_copy(weighted_evec_im, canonical_evec_im, name="ot_kp_occ_evec_im")
         CALL dbcsr_scale_by_vector(weighted_evec_im, alpha=orbital_scaling_factor, side='right')
         CALL dbcsr_copy(occ_re, chc_re, name="ot_kp_occ_re")
         CALL dbcsr_set(occ_re, 0.0_dp)
         CALL dbcsr_copy(occ_im, chc_im, name="ot_kp_occ_im")
         CALL dbcsr_set(occ_im, 0.0_dp)

         ! P_occ = V*diag(2*w_k*f)*V^H.  Blend only a small off-diagonal part into
         ! the rotation gradient, so the REF coordinates and minimizer history keep
         ! the stable diagonal occupation gauge.
         CALL dbcsr_multiply('N', 'T', 1.0_dp, weighted_evec_re, canonical_evec_re, 0.0_dp, occ_re)
         CALL dbcsr_multiply('N', 'T', 1.0_dp, weighted_evec_im, canonical_evec_im, 1.0_dp, occ_re)
         CALL dbcsr_multiply('N', 'T', 1.0_dp, weighted_evec_im, canonical_evec_re, 0.0_dp, occ_im)
         CALL dbcsr_multiply('N', 'T', -1.0_dp, weighted_evec_re, canonical_evec_im, 1.0_dp, occ_im)
         canonical_projector_weight = kpoint_ot_canonical_projector_weight(qs_ot_env%delta, &
                                                                           ot_kp_canonical_projector_max_weight)
         CALL dbcsr_scale(occ_re, canonical_projector_weight)
         CALL dbcsr_scale(occ_im, canonical_projector_weight)
         ALLOCATE (projector_diag(nmo))
         CALL dbcsr_get_diag(occ_re, projector_diag)
         projector_diag(:) = projector_diag(:) + &
                             (1.0_dp - canonical_projector_weight)*orbital_scaling_factor(:)
         CALL dbcsr_set_diag(occ_re, projector_diag)
         projector_diag(:) = 0.0_dp
         CALL dbcsr_set_diag(occ_im, projector_diag)
         DEALLOCATE (projector_diag)

         CALL dbcsr_multiply('N', 'N', 1.0_dp, chc_re, occ_re, 0.0_dp, b_re)
         CALL dbcsr_multiply('N', 'N', -1.0_dp, chc_im, occ_im, 1.0_dp, b_re)
         CALL dbcsr_multiply('N', 'N', 1.0_dp, chc_re, occ_im, 0.0_dp, b_im)
         CALL dbcsr_multiply('N', 'N', 1.0_dp, chc_im, occ_re, 1.0_dp, b_im)

         CALL dbcsr_release(canonical_evec_re)
         CALL dbcsr_release(canonical_evec_im)
         CALL dbcsr_release(diag_re)
         CALL dbcsr_release(diag_im)
         CALL dbcsr_release(occ_re)
         CALL dbcsr_release(occ_im)
         CALL dbcsr_release(weighted_evec_re)
         CALL dbcsr_release(weighted_evec_im)
         DEALLOCATE (canonical_energies, orbital_scaling_factor)
      ELSE
         CALL dbcsr_copy(b_re, chc_re)
         CALL dbcsr_copy(b_im, chc_im)

         ALLOCATE (orbital_scaling_factor(SIZE(occupation_numbers)))
         orbital_scaling_factor(:) = 2.0_dp*kpoint_weight*occupation_numbers(:)
         CALL dbcsr_scale_by_vector(b_re, alpha=orbital_scaling_factor, side='right')
         CALL dbcsr_scale_by_vector(b_im, alpha=orbital_scaling_factor, side='right')
         DEALLOCATE (orbital_scaling_factor)
      END IF

      ! dF/dU for the complex occupied-space rotation U(k).  The finite anti-Hermitian
      ! Frechet derivative below consumes the complete real/imaginary payload.
      CALL dbcsr_multiply('N', 'N', 1.0_dp, qs_ot_env%rot_mat_u, b_re, 0.0_dp, qs_ot_env%rot_mat_dedu)
      IF (ASSOCIATED(qs_ot_env%rot_mat_u_im)) THEN
         CALL dbcsr_multiply('N', 'N', 1.0_dp, qs_ot_env%rot_mat_u_im, b_im, 0.0_dp, tmp)
         CALL dbcsr_add(qs_ot_env%rot_mat_dedu, tmp, alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)

         CALL dbcsr_multiply('N', 'N', 1.0_dp, qs_ot_env%rot_mat_u, b_im, 0.0_dp, qs_ot_env%rot_mat_dedu_im)
         CALL dbcsr_multiply('N', 'N', 1.0_dp, qs_ot_env%rot_mat_u_im, b_re, 0.0_dp, tmp2)
         CALL dbcsr_add(qs_ot_env%rot_mat_dedu_im, tmp2, alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
      END IF

      qs_ot_env%canonical_tether_energy = 0.0_dp
      qs_ot_env%canonical_tether_weight_sum = 0.0_dp
      qs_ot_env%canonical_tether_shift = 0.0_dp
      qs_ot_env%canonical_tether_grad_max = 0.0_dp

      IF (qs_ot_env%settings%do_ener .AND. qs_ot_env%settings%occupation_preconditioner .AND. &
          qs_ot_env%settings%canonical_tether .AND. &
          .NOT. project_energy_subproblem .AND. &
          smear%do_smear .AND. smear%method == smear_fermi_dirac) THEN
         sigma = smear%electronic_temperature
         nmo = SIZE(qs_ot_env%ener_x)
         IF (sigma > 1.0E-14_dp .AND. SIZE(eigenvalues) >= nmo .AND. SIZE(occupation_numbers) >= nmo) THEN
            ALLOCATE (diagonal(nmo), tether_diag(nmo))
            CALL dbcsr_get_diag(chc_re, diagonal)
            tether_diag(:) = 0.0_dp
            tether_shift = 0.0_dp
            rayleigh_edge_gap = 0.0_dp
            tether_weight_sum = 0.0_dp
            DO imo = 1, nmo
               tether_response = kpoint_ot_canonical_tether_response(occupation_numbers(imo), eigenvalues(imo), &
                                                                      mu, sigma, maxocc)
               IF (tether_response > 0.0_dp) &
                  rayleigh_edge_gap = MAX(rayleigh_edge_gap, ABS(diagonal(imo) - qs_ot_env%ener_x(imo)))
            END DO
            tether_strength = kpoint_ot_canonical_tether_strength(qs_ot_env, rayleigh_edge_gap, sigma)
            DO imo = 1, nmo
               tether_weight = kpoint_ot_canonical_tether_weight(occupation_numbers(imo), eigenvalues(imo), &
                                                                  mu, sigma, maxocc, kpoint_weight, &
                                                                  tether_strength)
               IF (tether_weight <= 0.0_dp) CYCLE
               delta = diagonal(imo) - qs_ot_env%ener_x(imo)
               tether_shift = tether_shift + tether_weight*delta
               tether_weight_sum = tether_weight_sum + tether_weight
            END DO
            IF (tether_weight_sum > EPSILON(tether_weight_sum)) THEN
               tether_shift = tether_shift/tether_weight_sum
               canonical_tether_energy = 0.0_dp
               DO imo = 1, nmo
                  tether_weight = kpoint_ot_canonical_tether_weight(occupation_numbers(imo), eigenvalues(imo), &
                                                                     mu, sigma, maxocc, kpoint_weight, &
                                                                     tether_strength)
                  IF (tether_weight <= 0.0_dp) CYCLE
                  delta = diagonal(imo) - qs_ot_env%ener_x(imo) - tether_shift
                  tether_diag(imo) = tether_weight*delta
                  canonical_tether_energy = canonical_tether_energy + 0.5_dp*tether_weight*delta**2
               END DO
               nondiag_energy = nondiag_energy + canonical_tether_energy
               qs_ot_env%canonical_tether_energy = canonical_tether_energy
               qs_ot_env%canonical_tether_strength = tether_strength
               qs_ot_env%canonical_tether_weight_sum = tether_weight_sum
               qs_ot_env%canonical_tether_shift = tether_shift
               qs_ot_env%canonical_tether_grad_max = MAXVAL(ABS(tether_diag))

               CALL dbcsr_copy(b_re, chc_re)
               CALL dbcsr_set(b_re, 0.0_dp)
               CALL dbcsr_set_diag(b_re, tether_diag)
               CALL dbcsr_multiply('N', 'N', 1.0_dp, qs_ot_env%rot_mat_u, chc_re, 0.0_dp, tmp)
               CALL dbcsr_multiply('N', 'N', -1.0_dp, qs_ot_env%rot_mat_u_im, chc_im, 1.0_dp, tmp)
               CALL dbcsr_multiply('N', 'N', 1.0_dp, qs_ot_env%rot_mat_u, chc_im, 0.0_dp, tmp2)
               CALL dbcsr_multiply('N', 'N', 1.0_dp, qs_ot_env%rot_mat_u_im, chc_re, 1.0_dp, tmp2)
               CALL dbcsr_multiply('N', 'N', 2.0_dp, tmp, b_re, 1.0_dp, qs_ot_env%rot_mat_dedu)
               CALL dbcsr_multiply('N', 'N', 2.0_dp, tmp2, b_re, 1.0_dp, qs_ot_env%rot_mat_dedu_im)
            END IF
            DEALLOCATE (diagonal, tether_diag)
         END IF
      END IF

      IF (qs_ot_env%settings%add_nondiag_energy) THEN
         CPASSERT(qs_ot_env%settings%do_ener)
         nondiag_weight = qs_ot_env%settings%nondiag_energy_strength*kpoint_weight

         CALL dbcsr_copy(b_re, chc_re)
         CALL dbcsr_copy(b_im, chc_im)
         ALLOCATE (diagonal(SIZE(qs_ot_env%ener_x)))
         CALL dbcsr_get_diag(b_re, diagonal)
         diagonal(:) = diagonal(:) - qs_ot_env%ener_x(:)
         CALL dbcsr_set_diag(b_re, diagonal)
         DEALLOCATE (diagonal)

         CALL dbcsr_dot(b_re, b_re, nondiag_trace_re)
         CALL dbcsr_dot(b_im, b_im, nondiag_trace_im)
         nondiag_energy = nondiag_energy + &
                          0.5_dp*nondiag_weight*(nondiag_trace_re + nondiag_trace_im)

         ! Complex version of the Gamma-point NONDIAG_ENERGY tether:
         ! d/dU 0.5*k*||U^H A U - eps||^2 = 2*k*A*U*(U^H A U - eps).
         ! Since CHC=U^H A U, A*U is recovered as U*CHC for unitary U.
         CALL dbcsr_multiply('N', 'N', 1.0_dp, qs_ot_env%rot_mat_u, chc_re, 0.0_dp, tmp)
         CALL dbcsr_multiply('N', 'N', -1.0_dp, qs_ot_env%rot_mat_u_im, chc_im, 1.0_dp, tmp)
         CALL dbcsr_multiply('N', 'N', 1.0_dp, qs_ot_env%rot_mat_u, chc_im, 0.0_dp, tmp2)
         CALL dbcsr_multiply('N', 'N', 1.0_dp, qs_ot_env%rot_mat_u_im, chc_re, 1.0_dp, tmp2)

         CALL dbcsr_multiply('N', 'N', 2.0_dp*nondiag_weight, tmp, b_re, 1.0_dp, &
                             qs_ot_env%rot_mat_dedu)
         CALL dbcsr_multiply('N', 'N', -2.0_dp*nondiag_weight, tmp2, b_im, 1.0_dp, &
                             qs_ot_env%rot_mat_dedu)
         CALL dbcsr_multiply('N', 'N', 2.0_dp*nondiag_weight, tmp, b_im, 1.0_dp, &
                             qs_ot_env%rot_mat_dedu_im)
         CALL dbcsr_multiply('N', 'N', 2.0_dp*nondiag_weight, tmp2, b_re, 1.0_dp, &
                             qs_ot_env%rot_mat_dedu_im)
      END IF

      CALL qs_ot_rot_mat_derivative_complex(qs_ot_env)
      ! The k-point gradient is assembled for the C_base*U channel convention.
      ! Convert the complex rotation adjoint to the minimizer variable convention.
      CALL dbcsr_scale(qs_ot_env%rot_mat_gx, -1.0_dp)
      CALL dbcsr_scale(qs_ot_env%rot_mat_gx_im, -1.0_dp)

   END SUBROUTINE qs_scf_loop_prepare_ot_kp_channel_gradient

! **************************************************************************************************
!> \brief convert k-point band-energy derivatives to fixed-N Mermin free-energy gradients
!> \param kpoints ...
!> \param local_ot_env ...
!> \param first_channel ...
!> \param nspin_ot ...
!> \param kp_range ...
!> \param wkp ...
!> \param para_env_inter_kp ...
!> \param smear ...
! **************************************************************************************************
   SUBROUTINE qs_scf_loop_prepare_ot_kp_energy_gradient(kpoints, local_ot_env, first_channel, nspin_ot, &
                                                        kp_range, wkp, para_env_inter_kp, smear, &
                                                        trial_eigenvalues, trial_occupation, trial_mu)

      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: local_ot_env
      INTEGER, INTENT(IN)                                :: first_channel, nspin_ot
      INTEGER, DIMENSION(2), INTENT(IN)                  :: kp_range
      REAL(KIND=dp), DIMENSION(:), POINTER               :: wkp
      TYPE(mp_para_env_type), POINTER                    :: para_env_inter_kp
      TYPE(smear_type), INTENT(IN)                       :: smear
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: trial_eigenvalues, trial_occupation
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: trial_mu

      INTEGER                                            :: group, ikpoint, imo, ispin, local_channel, &
                                                            local_kpoint, ngroups, nmo
      REAL(KIND=dp), PARAMETER                           :: curvature_floor = 1.0_dp

      REAL(KIND=dp)                                      :: delta, gvalue, maxocc, mu, sigma, &
                                                            rayleigh_edge_gap, tether_response, &
                                                            tether_shift, tether_strength, &
                                                            tether_weight, tether_weight_sum, v
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: canonical_diag, gdotv, gsum, &
                                                            nondiag_diag, precond_count, &
                                                            precond_shift, precond_sum
      LOGICAL                                            :: active_preconditioner
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: schur_group_done
      TYPE(kpoint_env_type), POINTER                     :: kp

      IF (.NOT. local_ot_env(1)%settings%do_ener) RETURN

      ngroups = 1
      IF (nspin_ot == 2 .AND. smear%fixed_mag_mom > 0.0_dp) ngroups = 2
      ALLOCATE (gsum(ngroups), gdotv(ngroups), precond_count(ngroups), &
                precond_shift(ngroups), precond_sum(ngroups), schur_group_done(ngroups))
      gsum(:) = 0.0_dp
      gdotv(:) = 0.0_dp
      precond_count(:) = 0.0_dp
      precond_shift(:) = 0.0_dp
      precond_sum(:) = 0.0_dp
      schur_group_done(:) = .FALSE.
      active_preconditioner = local_ot_env(1)%settings%occupation_preconditioner .AND. &
                              smear%do_smear .AND. smear%method == smear_fermi_dirac

      IF (smear%do_smear) THEN
         SELECT CASE (smear%method)
         CASE (smear_fermi_dirac)
            sigma = smear%electronic_temperature
         CASE (smear_gaussian, smear_mp, smear_mv)
            sigma = smear%smearing_width
         CASE DEFAULT
            sigma = 0.0_dp
         END SELECT
      ELSE
         sigma = 0.0_dp
      END IF

      ! For FD smearing, g_i = -df_i/deps_i. Fixed electron number removes the
      ! chemical-potential mode, giving the projected response
      ! w_i*g_i*(<v>_g - v_i) with <v>_g = sum(w_i*g_i*v_i)/sum(w_i*g_i).
      DO local_kpoint = 1, SIZE(kpoints%kp_env)
         kp => kpoints%kp_env(local_kpoint)%kpoint_env
         ikpoint = kp%nkpoint
         DO ispin = 1, nspin_ot
            group = MERGE(ispin, 1, ngroups == nspin_ot)
            local_channel = qs_ot_channel_index(ispin, ikpoint, nspin_ot) - first_channel + 1
            CALL get_mo_set(kp%mos(1, ispin), maxocc=maxocc, nmo=nmo)
            mu = trial_mu(ispin)
            DO imo = 1, nmo
               gvalue = kpoint_ot_smearing_g(trial_occupation(imo, ikpoint, ispin), &
                                             trial_eigenvalues(imo, ikpoint, ispin), mu, &
                                             sigma, maxocc, smear%method)
               v = local_ot_env(local_channel)%ener_gx(imo) - trial_eigenvalues(imo, ikpoint, ispin) + mu
               gsum(group) = gsum(group) + wkp(ikpoint)*gvalue
               gdotv(group) = gdotv(group) + wkp(ikpoint)*gvalue*v
            END DO
         END DO
      END DO

      CALL para_env_inter_kp%sum(gsum)
      CALL para_env_inter_kp%sum(gdotv)

      DO local_kpoint = 1, SIZE(kpoints%kp_env)
         kp => kpoints%kp_env(local_kpoint)%kpoint_env
         ikpoint = kp%nkpoint
         DO ispin = 1, nspin_ot
            group = MERGE(ispin, 1, ngroups == nspin_ot)
            local_channel = qs_ot_channel_index(ispin, ikpoint, nspin_ot) - first_channel + 1
            CALL get_mo_set(kp%mos(1, ispin), maxocc=maxocc, nmo=nmo)
            mu = trial_mu(ispin)
            DO imo = 1, nmo
               gvalue = kpoint_ot_smearing_g(trial_occupation(imo, ikpoint, ispin), &
                                             trial_eigenvalues(imo, ikpoint, ispin), mu, &
                                             sigma, maxocc, smear%method)
               IF (ABS(gsum(group)) > EPSILON(gsum(group))) THEN
                  v = local_ot_env(local_channel)%ener_gx(imo) - trial_eigenvalues(imo, ikpoint, ispin) + mu
                  local_ot_env(local_channel)%ener_gx(imo) = &
                     wkp(ikpoint)*gvalue*(gdotv(group)/gsum(group) - v)
               ELSE
                  local_ot_env(local_channel)%ener_gx(imo) = 0.0_dp
               END IF
            END DO
         END DO
      END DO

      IF (local_ot_env(1)%settings%add_nondiag_energy) THEN
         DO local_kpoint = 1, SIZE(kpoints%kp_env)
            kp => kpoints%kp_env(local_kpoint)%kpoint_env
            ikpoint = kp%nkpoint
            DO ispin = 1, nspin_ot
               local_channel = qs_ot_channel_index(ispin, ikpoint, nspin_ot) - first_channel + 1
               CALL get_mo_set(kp%mos(1, ispin), nmo=nmo)
               ALLOCATE (nondiag_diag(nmo))
               CALL dbcsr_get_diag(local_ot_env(local_channel)%rot_mat_chc, nondiag_diag)
               nondiag_diag(:) = nondiag_diag(:) - local_ot_env(local_channel)%ener_x(:)
               local_ot_env(local_channel)%ener_gx(:) = local_ot_env(local_channel)%ener_gx(:) - &
                  local_ot_env(local_channel)%settings%nondiag_energy_strength*wkp(ikpoint)*nondiag_diag(:)
               DEALLOCATE (nondiag_diag)
            END DO
         END DO
      END IF

      IF (active_preconditioner .AND. local_ot_env(1)%settings%canonical_tether) THEN
         DO local_kpoint = 1, SIZE(kpoints%kp_env)
            kp => kpoints%kp_env(local_kpoint)%kpoint_env
            ikpoint = kp%nkpoint
            DO ispin = 1, nspin_ot
               local_channel = qs_ot_channel_index(ispin, ikpoint, nspin_ot) - first_channel + 1
               CALL get_mo_set(kp%mos(1, ispin), maxocc=maxocc, nmo=nmo)
               mu = trial_mu(ispin)
               ALLOCATE (canonical_diag(nmo))
               CALL dbcsr_get_diag(local_ot_env(local_channel)%rot_mat_chc, canonical_diag)
               tether_shift = 0.0_dp
               rayleigh_edge_gap = 0.0_dp
               tether_weight_sum = 0.0_dp
               DO imo = 1, nmo
                  tether_response = kpoint_ot_canonical_tether_response(trial_occupation(imo, ikpoint, ispin), &
                                                                        trial_eigenvalues(imo, ikpoint, ispin), mu, &
                                                                        sigma, maxocc)
                  IF (tether_response > 0.0_dp) &
                     rayleigh_edge_gap = MAX(rayleigh_edge_gap, &
                                             ABS(canonical_diag(imo) - local_ot_env(local_channel)%ener_x(imo)))
               END DO
               tether_strength = kpoint_ot_canonical_tether_strength(local_ot_env(local_channel), &
                                                                     rayleigh_edge_gap, sigma)
               DO imo = 1, nmo
                  tether_weight = kpoint_ot_canonical_tether_weight(trial_occupation(imo, ikpoint, ispin), &
                                                                     trial_eigenvalues(imo, ikpoint, ispin), &
                                                                     mu, sigma, maxocc, wkp(ikpoint), &
                                                                     tether_strength)
                  IF (tether_weight <= 0.0_dp) CYCLE
                  delta = canonical_diag(imo) - local_ot_env(local_channel)%ener_x(imo)
                  tether_shift = tether_shift + tether_weight*delta
                  tether_weight_sum = tether_weight_sum + tether_weight
               END DO
               IF (tether_weight_sum > EPSILON(tether_weight_sum)) THEN
                  tether_shift = tether_shift/tether_weight_sum
                  DO imo = 1, nmo
                     tether_weight = kpoint_ot_canonical_tether_weight(trial_occupation(imo, ikpoint, ispin), &
                                                                        trial_eigenvalues(imo, ikpoint, ispin), &
                                                                        mu, sigma, maxocc, wkp(ikpoint), &
                                                                        tether_strength)
                     IF (tether_weight <= 0.0_dp) CYCLE
                     delta = canonical_diag(imo) - local_ot_env(local_channel)%ener_x(imo) - tether_shift
                     local_ot_env(local_channel)%ener_gx(imo) = &
                        local_ot_env(local_channel)%ener_gx(imo) - tether_weight*delta
                  END DO
               END IF
               DEALLOCATE (canonical_diag)
            END DO
         END DO
      END IF

      IF (active_preconditioner) THEN
         CALL precondition_kpoint_ot_schur(schur_group_done)
         CALL precondition_kpoint_ot_energy_diagonal(schur_group_done)
      END IF

      IF (kpoint_ot_project_energy_subproblem(local_ot_env(1), smear)) THEN
         DO local_channel = 1, SIZE(local_ot_env)
            IF (ASSOCIATED(local_ot_env(local_channel)%ener_gx)) &
               local_ot_env(local_channel)%ener_gx(:) = 0.0_dp
            IF (ASSOCIATED(local_ot_env(local_channel)%ener_dx)) &
               local_ot_env(local_channel)%ener_dx(:) = 0.0_dp
            IF (ASSOCIATED(local_ot_env(local_channel)%ener_gx_old)) &
               local_ot_env(local_channel)%ener_gx_old(:) = 0.0_dp
         END DO
      END IF

      DEALLOCATE (gsum, gdotv, precond_count, precond_shift, precond_sum, schur_group_done)

      CPASSERT(kp_range(2) - kp_range(1) + 1 == SIZE(kpoints%kp_env))

   CONTAINS

      SUBROUTINE update_kpoint_ot_schur_line_search_quality()

      INTEGER                                           :: local_channel, samples
      LOGICAL                                           :: mermin_search
      REAL(KIND=dp)                                     :: calibrated_predicted_drop, calibrated_quality_ratio, &
                                                            constrained_calibrated_predicted_drop, &
                                                            constrained_calibrated_quality_ratio, &
                                                            constrained_predicted_drop, constrained_quality_ratio, &
                                                            measured_drop, predicted_drop, quality_ratio, &
                                                            occgap_shadow_predicted_drop, &
                                                            occgap_shadow_quality_ratio, &
                                                            projector_mode_predicted_drop, &
                                                            projector_mode_quality_ratio, &
                                                            projector_mode_response_drop_ratio, &
                                                               projector_mode_response_negative, &
                                                               projector_mode_step_negative, &
                                                               projector_mode_trust, &
                                                               projector_physical_predicted_drop, &
                                                               projector_physical_quality_ratio, &
                                                               projector_response_predicted_drop, &
                                                               projector_response_quality_ratio, &
                                                            residual_current_delta, residual_drop, &
                                                            residual_quality_ratio, &
                                                            split_shadow_calibrated_predicted_drop, &
                                                            split_shadow_calibrated_quality_ratio, &
                                                            split_shadow_predicted_drop, split_shadow_quality_ratio
      TYPE(cp_logger_type), POINTER                     :: logger

         IF (.NOT. local_ot_env(1)%schur_model_pending) RETURN
         IF (local_ot_env(1)%line_search_count <= 0) RETURN
         IF (.NOT. local_ot_env(1)%line_search_might_be_done) RETURN

         mermin_search = local_ot_env(1)%settings%occupation_preconditioner
         IF (mermin_search) THEN
            mermin_search = .FALSE.
            DO local_channel = 1, SIZE(local_ot_env)
               IF (local_ot_env(local_channel)%has_complex_kpoint_state) THEN
                  mermin_search = .TRUE.
                  EXIT
               END IF
            END DO
         END IF

         predicted_drop = local_ot_env(1)%schur_model_drop
         calibrated_predicted_drop = local_ot_env(1)%schur_model_calibrated_drop
         IF (calibrated_predicted_drop <= EPSILON(calibrated_predicted_drop)) &
            calibrated_predicted_drop = predicted_drop
         measured_drop = local_ot_env(1)%OT_energy(1) - local_ot_env(1)%etotal
         local_ot_env(1)%schur_model_measured_drop = measured_drop
         residual_current_delta = MAX(local_ot_env(1)%delta, 0.0_dp)
         residual_drop = MAX(0.0_dp, local_ot_env(1)%schur_model_delta - residual_current_delta)
         local_ot_env(1)%schur_model_residual_drop = residual_drop
         residual_quality_ratio = 0.0_dp
         IF (mermin_search .AND. local_ot_env(1)%schur_model_delta > EPSILON(residual_current_delta)) &
            residual_quality_ratio = residual_drop/local_ot_env(1)%schur_model_delta
         IF (local_ot_env(1)%schur_history_mixed) THEN
            logger => cp_get_default_logger()
            IF (logger%iter_info%print_level >= high_print_level .AND. logger%para_env%is_source()) THEN
               WRITE (cp_logger_get_default_unit_nr(logger), '(A,1X,I6,3(1X,ES16.8))') &
                  " OT Schur LS quality skipped mixed CG history attempt pred cal meas:", &
                  local_ot_env(1)%schur_model_attempt, predicted_drop, &
                  calibrated_predicted_drop, measured_drop
            END IF
            local_ot_env(1)%schur_model_pending = .FALSE.
            local_ot_env(1)%schur_constrained_pending = .FALSE.
            local_ot_env(1)%schur_split_shadow_pending = .FALSE.
               local_ot_env(1)%schur_occgap_shadow_pending = .FALSE.
               local_ot_env(1)%schur_projector_mode_pending = .FALSE.
               local_ot_env(1)%schur_projector_response_pending = .FALSE.
               local_ot_env(1)%schur_projector_physical_pending = .FALSE.
               local_ot_env(1)%schur_history_mixed = .FALSE.
            RETURN
         END IF
         quality_ratio = 0.0_dp
         calibrated_quality_ratio = 0.0_dp
         IF (mermin_search .AND. predicted_drop > EPSILON(predicted_drop)) &
            quality_ratio = measured_drop/predicted_drop
         IF (mermin_search .AND. calibrated_predicted_drop > EPSILON(calibrated_predicted_drop)) &
            calibrated_quality_ratio = measured_drop/calibrated_predicted_drop

         IF (mermin_search .AND. predicted_drop > EPSILON(predicted_drop) .AND. measured_drop > 0.0_dp) THEN
            local_ot_env(1)%schur_model_quality = quality_ratio
            local_ot_env(1)%schur_model_calibrated_quality = calibrated_quality_ratio
            samples = MIN(local_ot_env(1)%schur_model_quality_samples + 1, HUGE(samples) - 1)
            local_ot_env(1)%schur_model_quality_samples = samples
            IF (samples == 1) THEN
               local_ot_env(1)%schur_model_quality_ema = quality_ratio
               local_ot_env(1)%schur_model_calibrated_quality_ema = calibrated_quality_ratio
            ELSE
               local_ot_env(1)%schur_model_quality_ema = &
                  (1.0_dp - ot_kp_schur_quality_ema_alpha)*local_ot_env(1)%schur_model_quality_ema + &
                  ot_kp_schur_quality_ema_alpha*quality_ratio
               local_ot_env(1)%schur_model_calibrated_quality_ema = &
                  (1.0_dp - ot_kp_schur_quality_ema_alpha)* &
                  local_ot_env(1)%schur_model_calibrated_quality_ema + &
                  ot_kp_schur_quality_ema_alpha*calibrated_quality_ratio
            END IF
            local_ot_env(1)%schur_model_quality_ready = &
               samples >= ot_kp_schur_quality_min_samples .AND. &
               local_ot_env(1)%schur_model_quality_ema >= ot_kp_schur_quality_min_ratio .AND. &
               local_ot_env(1)%schur_model_quality_ema <= ot_kp_schur_quality_max_ratio .AND. &
               quality_ratio >= 0.5_dp*ot_kp_schur_quality_min_ratio .AND. &
               quality_ratio <= 2.0_dp*ot_kp_schur_quality_max_ratio
         ELSE
            local_ot_env(1)%schur_model_quality = 0.0_dp
            local_ot_env(1)%schur_model_quality_ema = 0.0_dp
            local_ot_env(1)%schur_model_calibrated_quality = 0.0_dp
            local_ot_env(1)%schur_model_calibrated_quality_ema = 0.0_dp
            local_ot_env(1)%schur_model_quality_samples = 0
            local_ot_env(1)%schur_model_quality_ready = .FALSE.
         END IF
         constrained_predicted_drop = local_ot_env(1)%schur_constrained_drop
         constrained_calibrated_predicted_drop = local_ot_env(1)%schur_constrained_calibrated_drop
         IF (constrained_calibrated_predicted_drop <= EPSILON(constrained_calibrated_predicted_drop)) &
            constrained_calibrated_predicted_drop = constrained_predicted_drop
         constrained_quality_ratio = 0.0_dp
         constrained_calibrated_quality_ratio = 0.0_dp
         IF (mermin_search .AND. local_ot_env(1)%schur_constrained_pending .AND. &
             constrained_predicted_drop > EPSILON(constrained_predicted_drop) .AND. measured_drop > 0.0_dp) THEN
            constrained_quality_ratio = measured_drop/constrained_predicted_drop
            IF (constrained_calibrated_predicted_drop > EPSILON(constrained_calibrated_predicted_drop)) &
               constrained_calibrated_quality_ratio = measured_drop/constrained_calibrated_predicted_drop
            local_ot_env(1)%schur_constrained_quality = constrained_quality_ratio
            local_ot_env(1)%schur_constrained_calibrated_quality = constrained_calibrated_quality_ratio
            samples = MIN(local_ot_env(1)%schur_constrained_quality_samples + 1, HUGE(samples) - 1)
            local_ot_env(1)%schur_constrained_quality_samples = samples
            IF (samples == 1) THEN
               local_ot_env(1)%schur_constrained_quality_ema = constrained_quality_ratio
               local_ot_env(1)%schur_constrained_calibrated_quality_ema = constrained_calibrated_quality_ratio
            ELSE
               local_ot_env(1)%schur_constrained_quality_ema = &
                  (1.0_dp - ot_kp_schur_quality_ema_alpha)*local_ot_env(1)%schur_constrained_quality_ema + &
                  ot_kp_schur_quality_ema_alpha*constrained_quality_ratio
               local_ot_env(1)%schur_constrained_calibrated_quality_ema = &
                  (1.0_dp - ot_kp_schur_quality_ema_alpha)* &
                  local_ot_env(1)%schur_constrained_calibrated_quality_ema + &
                  ot_kp_schur_quality_ema_alpha*constrained_calibrated_quality_ratio
            END IF
            local_ot_env(1)%schur_constrained_drop_sum = &
               local_ot_env(1)%schur_constrained_drop_sum + constrained_predicted_drop
            local_ot_env(1)%schur_constrained_calibrated_drop_sum = &
               local_ot_env(1)%schur_constrained_calibrated_drop_sum + constrained_calibrated_predicted_drop
            local_ot_env(1)%schur_constrained_measured_sum = &
               local_ot_env(1)%schur_constrained_measured_sum + measured_drop
         ELSE IF (local_ot_env(1)%schur_constrained_pending) THEN
            local_ot_env(1)%schur_constrained_quality = 0.0_dp
            local_ot_env(1)%schur_constrained_quality_ema = 0.0_dp
            local_ot_env(1)%schur_constrained_calibrated_quality = 0.0_dp
            local_ot_env(1)%schur_constrained_calibrated_quality_ema = 0.0_dp
            local_ot_env(1)%schur_constrained_quality_samples = 0
         END IF
         split_shadow_predicted_drop = local_ot_env(1)%schur_split_shadow_drop
         split_shadow_calibrated_predicted_drop = local_ot_env(1)%schur_split_shadow_calibrated_drop
         IF (split_shadow_calibrated_predicted_drop <= EPSILON(split_shadow_calibrated_predicted_drop)) &
            split_shadow_calibrated_predicted_drop = split_shadow_predicted_drop
         split_shadow_quality_ratio = 0.0_dp
         split_shadow_calibrated_quality_ratio = 0.0_dp
         IF (mermin_search .AND. local_ot_env(1)%schur_split_shadow_pending .AND. &
             split_shadow_predicted_drop > EPSILON(split_shadow_predicted_drop) .AND. measured_drop > 0.0_dp) THEN
            split_shadow_quality_ratio = measured_drop/split_shadow_predicted_drop
            IF (split_shadow_calibrated_predicted_drop > EPSILON(split_shadow_calibrated_predicted_drop)) &
               split_shadow_calibrated_quality_ratio = measured_drop/split_shadow_calibrated_predicted_drop
            local_ot_env(1)%schur_split_shadow_quality = split_shadow_quality_ratio
            local_ot_env(1)%schur_split_shadow_calibrated_quality = split_shadow_calibrated_quality_ratio
            samples = MIN(local_ot_env(1)%schur_split_shadow_quality_samples + 1, HUGE(samples) - 1)
            local_ot_env(1)%schur_split_shadow_quality_samples = samples
            IF (samples == 1) THEN
               local_ot_env(1)%schur_split_shadow_quality_ema = split_shadow_quality_ratio
               local_ot_env(1)%schur_split_shadow_calibrated_quality_ema = split_shadow_calibrated_quality_ratio
            ELSE
               local_ot_env(1)%schur_split_shadow_quality_ema = &
                  (1.0_dp - ot_kp_schur_quality_ema_alpha)*local_ot_env(1)%schur_split_shadow_quality_ema + &
                  ot_kp_schur_quality_ema_alpha*split_shadow_quality_ratio
               local_ot_env(1)%schur_split_shadow_calibrated_quality_ema = &
                  (1.0_dp - ot_kp_schur_quality_ema_alpha)* &
                  local_ot_env(1)%schur_split_shadow_calibrated_quality_ema + &
                  ot_kp_schur_quality_ema_alpha*split_shadow_calibrated_quality_ratio
            END IF
            local_ot_env(1)%schur_split_shadow_drop_sum = &
               local_ot_env(1)%schur_split_shadow_drop_sum + split_shadow_predicted_drop
            local_ot_env(1)%schur_split_shadow_calibrated_drop_sum = &
               local_ot_env(1)%schur_split_shadow_calibrated_drop_sum + split_shadow_calibrated_predicted_drop
            local_ot_env(1)%schur_split_shadow_measured_sum = &
               local_ot_env(1)%schur_split_shadow_measured_sum + measured_drop
         END IF
         occgap_shadow_predicted_drop = local_ot_env(1)%schur_occgap_shadow_drop
         occgap_shadow_quality_ratio = 0.0_dp
         IF (mermin_search .AND. local_ot_env(1)%schur_occgap_shadow_pending .AND. &
             occgap_shadow_predicted_drop > EPSILON(occgap_shadow_predicted_drop) .AND. measured_drop > 0.0_dp) THEN
            occgap_shadow_quality_ratio = measured_drop/occgap_shadow_predicted_drop
            local_ot_env(1)%schur_occgap_shadow_quality = occgap_shadow_quality_ratio
            samples = MIN(local_ot_env(1)%schur_occgap_shadow_quality_samples + 1, HUGE(samples) - 1)
            local_ot_env(1)%schur_occgap_shadow_quality_samples = samples
            IF (samples == 1) THEN
               local_ot_env(1)%schur_occgap_shadow_quality_ema = occgap_shadow_quality_ratio
            ELSE
               local_ot_env(1)%schur_occgap_shadow_quality_ema = &
                  (1.0_dp - ot_kp_schur_quality_ema_alpha)* &
                  local_ot_env(1)%schur_occgap_shadow_quality_ema + &
                  ot_kp_schur_quality_ema_alpha*occgap_shadow_quality_ratio
            END IF
            local_ot_env(1)%schur_occgap_shadow_drop_sum = &
               local_ot_env(1)%schur_occgap_shadow_drop_sum + occgap_shadow_predicted_drop
            local_ot_env(1)%schur_occgap_shadow_measured_sum = &
               local_ot_env(1)%schur_occgap_shadow_measured_sum + measured_drop
         END IF
         IF (mermin_search .AND. predicted_drop > EPSILON(predicted_drop)) THEN
            SELECT CASE (local_ot_env(1)%schur_model_attempt)
            CASE (1, 2)
               local_ot_env(1)%schur_split_quality = quality_ratio
               local_ot_env(1)%schur_split_calibrated_quality = calibrated_quality_ratio
               samples = MIN(local_ot_env(1)%schur_split_quality_samples + 1, HUGE(samples) - 1)
               local_ot_env(1)%schur_split_quality_samples = samples
               IF (samples == 1) THEN
                  local_ot_env(1)%schur_split_quality_ema = quality_ratio
                  local_ot_env(1)%schur_split_calibrated_quality_ema = calibrated_quality_ratio
               ELSE
                  local_ot_env(1)%schur_split_quality_ema = &
                     (1.0_dp - ot_kp_schur_quality_ema_alpha)*local_ot_env(1)%schur_split_quality_ema + &
                     ot_kp_schur_quality_ema_alpha*quality_ratio
                  local_ot_env(1)%schur_split_calibrated_quality_ema = &
                     (1.0_dp - ot_kp_schur_quality_ema_alpha)* &
                     local_ot_env(1)%schur_split_calibrated_quality_ema + &
                     ot_kp_schur_quality_ema_alpha*calibrated_quality_ratio
               END IF
               local_ot_env(1)%schur_split_drop_sum = local_ot_env(1)%schur_split_drop_sum + predicted_drop
               local_ot_env(1)%schur_split_calibrated_drop_sum = &
                  local_ot_env(1)%schur_split_calibrated_drop_sum + calibrated_predicted_drop
               local_ot_env(1)%schur_split_measured_sum = local_ot_env(1)%schur_split_measured_sum + measured_drop
               local_ot_env(1)%schur_split_residual_quality = residual_quality_ratio
               IF (samples == 1) THEN
                  local_ot_env(1)%schur_split_residual_quality_ema = residual_quality_ratio
               ELSE
                  local_ot_env(1)%schur_split_residual_quality_ema = &
                     (1.0_dp - ot_kp_schur_quality_ema_alpha)* &
                     local_ot_env(1)%schur_split_residual_quality_ema + &
                     ot_kp_schur_quality_ema_alpha*residual_quality_ratio
               END IF
            CASE (3)
               local_ot_env(1)%schur_coupled_quality = quality_ratio
               samples = MIN(local_ot_env(1)%schur_coupled_quality_samples + 1, HUGE(samples) - 1)
               local_ot_env(1)%schur_coupled_quality_samples = samples
               IF (samples == 1) THEN
                  local_ot_env(1)%schur_coupled_quality_ema = quality_ratio
               ELSE
                  local_ot_env(1)%schur_coupled_quality_ema = &
                     (1.0_dp - ot_kp_schur_quality_ema_alpha)*local_ot_env(1)%schur_coupled_quality_ema + &
                     ot_kp_schur_quality_ema_alpha*quality_ratio
               END IF
               local_ot_env(1)%schur_coupled_drop_sum = local_ot_env(1)%schur_coupled_drop_sum + predicted_drop
               local_ot_env(1)%schur_coupled_measured_sum = local_ot_env(1)%schur_coupled_measured_sum + measured_drop
               local_ot_env(1)%schur_coupled_residual_quality = residual_quality_ratio
               IF (samples == 1) THEN
                  local_ot_env(1)%schur_coupled_residual_quality_ema = residual_quality_ratio
               ELSE
                  local_ot_env(1)%schur_coupled_residual_quality_ema = &
                     (1.0_dp - ot_kp_schur_quality_ema_alpha)* &
                     local_ot_env(1)%schur_coupled_residual_quality_ema + &
                     ot_kp_schur_quality_ema_alpha*residual_quality_ratio
               END IF
            END SELECT
         END IF

         projector_mode_predicted_drop = local_ot_env(1)%schur_projector_mode_drop
         projector_mode_quality_ratio = 0.0_dp
         projector_mode_response_drop_ratio = 0.0_dp
         projector_mode_step_negative = 0.0_dp
         projector_mode_response_negative = 0.0_dp
         projector_mode_trust = 1.0_dp
         IF (projector_mode_predicted_drop > EPSILON(projector_mode_predicted_drop)) THEN
            projector_mode_step_negative = &
               local_ot_env(1)%schur_projector_mode_step_negative/projector_mode_predicted_drop
            projector_mode_trust = &
               local_ot_env(1)%schur_projector_mode_trust/projector_mode_predicted_drop
            IF (local_ot_env(1)%schur_projector_mode_response_drop > &
                EPSILON(local_ot_env(1)%schur_projector_mode_response_drop)) &
               projector_mode_response_negative = &
                  local_ot_env(1)%schur_projector_mode_response_negative/ &
                  local_ot_env(1)%schur_projector_mode_response_drop
            projector_mode_response_drop_ratio = &
               local_ot_env(1)%schur_projector_mode_response_drop/projector_mode_predicted_drop
         END IF
         IF (mermin_search .AND. local_ot_env(1)%schur_projector_mode_pending .AND. &
             projector_mode_predicted_drop > EPSILON(projector_mode_predicted_drop) .AND. &
             measured_drop > 0.0_dp) THEN
            projector_mode_quality_ratio = measured_drop/projector_mode_predicted_drop
            local_ot_env(1)%schur_projector_mode_quality = projector_mode_quality_ratio
            samples = MIN(local_ot_env(1)%schur_projector_mode_quality_samples + 1, HUGE(samples) - 1)
            local_ot_env(1)%schur_projector_mode_quality_samples = samples
            IF (samples == 1) THEN
               local_ot_env(1)%schur_projector_mode_quality_ema = projector_mode_quality_ratio
               local_ot_env(1)%schur_projector_mode_step_negative_ema = projector_mode_step_negative
               local_ot_env(1)%schur_projector_mode_response_negative_ema = projector_mode_response_negative
               local_ot_env(1)%schur_projector_mode_trust_ema = projector_mode_trust
            ELSE
               local_ot_env(1)%schur_projector_mode_quality_ema = &
                  (1.0_dp - ot_kp_schur_quality_ema_alpha)* &
                  local_ot_env(1)%schur_projector_mode_quality_ema + &
                  ot_kp_schur_quality_ema_alpha*projector_mode_quality_ratio
               local_ot_env(1)%schur_projector_mode_step_negative_ema = &
                  (1.0_dp - ot_kp_schur_quality_ema_alpha)* &
                  local_ot_env(1)%schur_projector_mode_step_negative_ema + &
                  ot_kp_schur_quality_ema_alpha*projector_mode_step_negative
               local_ot_env(1)%schur_projector_mode_response_negative_ema = &
                  (1.0_dp - ot_kp_schur_quality_ema_alpha)* &
                  local_ot_env(1)%schur_projector_mode_response_negative_ema + &
                  ot_kp_schur_quality_ema_alpha*projector_mode_response_negative
               local_ot_env(1)%schur_projector_mode_trust_ema = &
                  (1.0_dp - ot_kp_schur_quality_ema_alpha)* &
                  local_ot_env(1)%schur_projector_mode_trust_ema + &
                  ot_kp_schur_quality_ema_alpha*projector_mode_trust
            END IF
            local_ot_env(1)%schur_projector_mode_drop_sum = &
               local_ot_env(1)%schur_projector_mode_drop_sum + projector_mode_predicted_drop
            local_ot_env(1)%schur_projector_mode_measured_sum = &
               local_ot_env(1)%schur_projector_mode_measured_sum + measured_drop
            local_ot_env(1)%schur_projector_mode_response_drop_sum = &
               local_ot_env(1)%schur_projector_mode_response_drop_sum + &
               local_ot_env(1)%schur_projector_mode_response_drop
         END IF

         projector_response_predicted_drop = local_ot_env(1)%schur_projector_response_drop
         projector_response_quality_ratio = 0.0_dp
         IF (mermin_search .AND. local_ot_env(1)%schur_projector_response_pending .AND. &
             projector_response_predicted_drop > EPSILON(projector_response_predicted_drop) .AND. &
             measured_drop > 0.0_dp) THEN
            projector_response_quality_ratio = measured_drop/projector_response_predicted_drop
            local_ot_env(1)%schur_projector_response_quality = projector_response_quality_ratio
            samples = MIN(local_ot_env(1)%schur_projector_response_quality_samples + 1, HUGE(samples) - 1)
            local_ot_env(1)%schur_projector_response_quality_samples = samples
            IF (samples == 1) THEN
               local_ot_env(1)%schur_projector_response_quality_ema = projector_response_quality_ratio
               local_ot_env(1)%schur_projector_response_weight_ema = &
                  local_ot_env(1)%schur_projector_response_weight
               local_ot_env(1)%schur_projector_response_negative_ema = &
                  local_ot_env(1)%schur_projector_response_negative
            ELSE
               local_ot_env(1)%schur_projector_response_quality_ema = &
                  (1.0_dp - ot_kp_schur_quality_ema_alpha)* &
                  local_ot_env(1)%schur_projector_response_quality_ema + &
                  ot_kp_schur_quality_ema_alpha*projector_response_quality_ratio
               local_ot_env(1)%schur_projector_response_weight_ema = &
                  (1.0_dp - ot_kp_schur_quality_ema_alpha)* &
                  local_ot_env(1)%schur_projector_response_weight_ema + &
                  ot_kp_schur_quality_ema_alpha*local_ot_env(1)%schur_projector_response_weight
               local_ot_env(1)%schur_projector_response_negative_ema = &
                  (1.0_dp - ot_kp_schur_quality_ema_alpha)* &
                  local_ot_env(1)%schur_projector_response_negative_ema + &
                  ot_kp_schur_quality_ema_alpha*local_ot_env(1)%schur_projector_response_negative
            END IF
            local_ot_env(1)%schur_projector_response_drop_sum = &
               local_ot_env(1)%schur_projector_response_drop_sum + projector_response_predicted_drop
            local_ot_env(1)%schur_projector_response_measured_sum = &
               local_ot_env(1)%schur_projector_response_measured_sum + measured_drop
         END IF
            local_ot_env(1)%schur_projector_response_candidate_enabled = &
               local_ot_env(1)%schur_projector_response_quality_samples >= &
               ot_kp_schur_projector_response_min_samples .AND. &
               local_ot_env(1)%schur_projector_response_quality_ema >= &
            ot_kp_schur_projector_response_quality_min_ratio .AND. &
            local_ot_env(1)%schur_projector_response_quality_ema <= &
            ot_kp_schur_projector_response_quality_max_ratio .AND. &
            local_ot_env(1)%schur_projector_response_quality >= &
            0.5_dp*ot_kp_schur_projector_response_quality_min_ratio .AND. &
               local_ot_env(1)%schur_projector_response_quality <= &
               2.0_dp*ot_kp_schur_projector_response_quality_max_ratio

            projector_physical_predicted_drop = local_ot_env(1)%schur_projector_physical_drop
            projector_physical_quality_ratio = 0.0_dp
            IF (mermin_search .AND. local_ot_env(1)%schur_projector_physical_pending .AND. &
                projector_physical_predicted_drop > EPSILON(projector_physical_predicted_drop) .AND. &
                measured_drop > 0.0_dp) THEN
               projector_physical_quality_ratio = measured_drop/projector_physical_predicted_drop
               local_ot_env(1)%schur_projector_physical_quality = projector_physical_quality_ratio
               samples = MIN(local_ot_env(1)%schur_projector_physical_quality_samples + 1, HUGE(samples) - 1)
               local_ot_env(1)%schur_projector_physical_quality_samples = samples
               IF (samples == 1) THEN
                  local_ot_env(1)%schur_projector_physical_quality_ema = projector_physical_quality_ratio
                  local_ot_env(1)%schur_projector_physical_weight_ema = &
                     local_ot_env(1)%schur_projector_physical_weight
                  local_ot_env(1)%schur_projector_physical_negative_ema = &
                     local_ot_env(1)%schur_projector_physical_negative
               ELSE
                  local_ot_env(1)%schur_projector_physical_quality_ema = &
                     (1.0_dp - ot_kp_schur_quality_ema_alpha)* &
                     local_ot_env(1)%schur_projector_physical_quality_ema + &
                     ot_kp_schur_quality_ema_alpha*projector_physical_quality_ratio
                  local_ot_env(1)%schur_projector_physical_weight_ema = &
                     (1.0_dp - ot_kp_schur_quality_ema_alpha)* &
                     local_ot_env(1)%schur_projector_physical_weight_ema + &
                     ot_kp_schur_quality_ema_alpha*local_ot_env(1)%schur_projector_physical_weight
                  local_ot_env(1)%schur_projector_physical_negative_ema = &
                     (1.0_dp - ot_kp_schur_quality_ema_alpha)* &
                     local_ot_env(1)%schur_projector_physical_negative_ema + &
                     ot_kp_schur_quality_ema_alpha*local_ot_env(1)%schur_projector_physical_negative
               END IF
               local_ot_env(1)%schur_projector_physical_drop_sum = &
                  local_ot_env(1)%schur_projector_physical_drop_sum + projector_physical_predicted_drop
               local_ot_env(1)%schur_projector_physical_measured_sum = &
                  local_ot_env(1)%schur_projector_physical_measured_sum + measured_drop
            END IF
            ! A too-large physical quality ratio means the shadow model is
            ! conservative, not unstable; selection still uses the clamped
            ! expected-drop quality below.
            local_ot_env(1)%schur_projector_physical_candidate_enabled = &
               local_ot_env(1)%schur_projector_physical_quality_samples >= &
               ot_kp_schur_projector_response_min_samples .AND. &
               local_ot_env(1)%schur_projector_physical_quality_ema >= &
               ot_kp_schur_projector_response_quality_min_ratio .AND. &
               local_ot_env(1)%schur_projector_physical_quality >= &
               0.5_dp*ot_kp_schur_projector_response_quality_min_ratio

            local_ot_env(1)%schur_split_candidate_enabled = local_ot_env(1)%schur_model_quality_ready
         IF (local_ot_env(1)%schur_split_quality_samples >= ot_kp_schur_calibration_min_samples) &
            local_ot_env(1)%schur_split_candidate_enabled = local_ot_env(1)%schur_split_candidate_enabled .AND. &
               local_ot_env(1)%schur_split_calibrated_quality_ema >= &
               ot_kp_schur_calibrated_quality_min_ratio .AND. &
               local_ot_env(1)%schur_split_calibrated_quality_ema <= &
               ot_kp_schur_calibrated_quality_max_ratio .AND. &
               local_ot_env(1)%schur_split_calibrated_quality >= &
               0.5_dp*ot_kp_schur_calibrated_quality_min_ratio .AND. &
               local_ot_env(1)%schur_split_calibrated_quality <= &
               2.0_dp*ot_kp_schur_calibrated_quality_max_ratio
         local_ot_env(1)%schur_constrained_candidate_enabled = &
            local_ot_env(1)%schur_split_candidate_enabled .AND. &
            local_ot_env(1)%schur_constrained_quality_samples >= ot_kp_schur_constrained_min_samples .AND. &
            local_ot_env(1)%schur_constrained_calibrated_quality_ema >= &
            ot_kp_schur_constrained_advantage_ratio* &
            MAX(ot_kp_schur_calibrated_quality_min_ratio, local_ot_env(1)%schur_split_calibrated_quality_ema) .AND. &
            local_ot_env(1)%schur_constrained_calibrated_quality_ema <= &
            ot_kp_schur_calibrated_quality_max_ratio .AND. &
            local_ot_env(1)%schur_constrained_calibrated_quality >= &
            0.5_dp*ot_kp_schur_calibrated_quality_min_ratio .AND. &
            local_ot_env(1)%schur_constrained_calibrated_quality <= &
            2.0_dp*ot_kp_schur_calibrated_quality_max_ratio
         local_ot_env(1)%schur_model_pending = .FALSE.
         local_ot_env(1)%schur_constrained_pending = .FALSE.
         local_ot_env(1)%schur_split_shadow_pending = .FALSE.
            local_ot_env(1)%schur_occgap_shadow_pending = .FALSE.
            local_ot_env(1)%schur_projector_mode_pending = .FALSE.
            local_ot_env(1)%schur_projector_response_pending = .FALSE.
            local_ot_env(1)%schur_projector_physical_pending = .FALSE.

         logger => cp_get_default_logger()
         IF (logger%iter_info%print_level >= high_print_level .AND. logger%para_env%is_source()) THEN
            WRITE (cp_logger_get_default_unit_nr(logger), '(A,1X,I6,9(1X,ES16.8),1X,I6,2(1X,L1))') &
               " OT Schur LS quality attempt raw_pred cal_pred meas raw_ratio raw_ema cal_ratio cal_ema "// &
               "gap_ratio delta samples ready enabled:", &
               local_ot_env(1)%schur_model_attempt, &
               predicted_drop, calibrated_predicted_drop, measured_drop, quality_ratio, &
               local_ot_env(1)%schur_model_quality_ema, calibrated_quality_ratio, &
               local_ot_env(1)%schur_model_calibrated_quality_ema, &
               local_ot_env(1)%schur_model_gap_ratio, local_ot_env(1)%delta, &
               local_ot_env(1)%schur_model_quality_samples, &
               local_ot_env(1)%schur_model_quality_ready, local_ot_env(1)%schur_split_candidate_enabled
            WRITE (cp_logger_get_default_unit_nr(logger), '(A,2(1X,I6),14(1X,ES16.8))') &
               " OT Schur LS history split_n coupled_n split_ratio split_ema split_cal split_cal_ema "// &
               "split_resid split_resid_ema split_pred split_meas "// &
               "coupled_ratio coupled_ema coupled_resid coupled_resid_ema coupled_pred coupled_meas:", &
               local_ot_env(1)%schur_split_quality_samples, local_ot_env(1)%schur_coupled_quality_samples, &
               local_ot_env(1)%schur_split_quality, local_ot_env(1)%schur_split_quality_ema, &
               local_ot_env(1)%schur_split_calibrated_quality, &
               local_ot_env(1)%schur_split_calibrated_quality_ema, &
               local_ot_env(1)%schur_split_residual_quality, &
               local_ot_env(1)%schur_split_residual_quality_ema, &
               local_ot_env(1)%schur_split_drop_sum, local_ot_env(1)%schur_split_measured_sum, &
               local_ot_env(1)%schur_coupled_quality, local_ot_env(1)%schur_coupled_quality_ema, &
               local_ot_env(1)%schur_coupled_residual_quality, &
               local_ot_env(1)%schur_coupled_residual_quality_ema, &
               local_ot_env(1)%schur_coupled_drop_sum, local_ot_env(1)%schur_coupled_measured_sum
            WRITE (cp_logger_get_default_unit_nr(logger), '(A,1X,I6,9(1X,ES16.8))') &
               " OT Schur LS split shadow n ratio ema cal cal_ema pred cpred meas latest_cal:", &
               local_ot_env(1)%schur_split_shadow_quality_samples, &
               local_ot_env(1)%schur_split_shadow_quality, local_ot_env(1)%schur_split_shadow_quality_ema, &
               local_ot_env(1)%schur_split_shadow_calibrated_quality, &
               local_ot_env(1)%schur_split_shadow_calibrated_quality_ema, &
               local_ot_env(1)%schur_split_shadow_drop_sum, &
               local_ot_env(1)%schur_split_shadow_calibrated_drop_sum, &
               local_ot_env(1)%schur_split_shadow_measured_sum, &
               split_shadow_calibrated_quality_ratio
            WRITE (cp_logger_get_default_unit_nr(logger), '(A,1X,I6,6(1X,ES16.8))') &
               " OT Schur LS occ-gap shadow n ratio ema pred meas gap latest:", &
               local_ot_env(1)%schur_occgap_shadow_quality_samples, &
               local_ot_env(1)%schur_occgap_shadow_quality, &
               local_ot_env(1)%schur_occgap_shadow_quality_ema, &
               local_ot_env(1)%schur_occgap_shadow_drop_sum, &
               local_ot_env(1)%schur_occgap_shadow_measured_sum, &
               local_ot_env(1)%schur_occgap_shadow_gap_ratio, &
               occgap_shadow_predicted_drop
            WRITE (cp_logger_get_default_unit_nr(logger), '(A,1X,I6,10(1X,ES16.8),1X,L1)') &
               " OT Schur LS constrained n ratio ema cal cal_ema pred cpred meas gap split_cal_ema split_cpred enabled:", &
               local_ot_env(1)%schur_constrained_quality_samples, &
               local_ot_env(1)%schur_constrained_quality, local_ot_env(1)%schur_constrained_quality_ema, &
               local_ot_env(1)%schur_constrained_calibrated_quality, &
               local_ot_env(1)%schur_constrained_calibrated_quality_ema, &
               local_ot_env(1)%schur_constrained_drop_sum, &
               local_ot_env(1)%schur_constrained_calibrated_drop_sum, &
               local_ot_env(1)%schur_constrained_measured_sum, local_ot_env(1)%schur_constrained_gap_ratio, &
               local_ot_env(1)%schur_split_calibrated_quality_ema, &
               local_ot_env(1)%schur_split_calibrated_drop_sum, &
               local_ot_env(1)%schur_constrained_candidate_enabled
            WRITE (cp_logger_get_default_unit_nr(logger), '(A,1X,I6,11(1X,ES16.8))') &
               " OT Schur LS projector modes n ratio ema pred meas rdrop rratio step rneg trust tema latest:", &
               local_ot_env(1)%schur_projector_mode_quality_samples, &
               local_ot_env(1)%schur_projector_mode_quality, &
               local_ot_env(1)%schur_projector_mode_quality_ema, &
               local_ot_env(1)%schur_projector_mode_drop_sum, &
               local_ot_env(1)%schur_projector_mode_measured_sum, &
               local_ot_env(1)%schur_projector_mode_response_drop_sum, &
               projector_mode_response_drop_ratio, projector_mode_step_negative, &
               projector_mode_response_negative, projector_mode_trust, &
               local_ot_env(1)%schur_projector_mode_trust_ema, &
               local_ot_env(1)%schur_projector_mode_response_drop
            WRITE (cp_logger_get_default_unit_nr(logger), '(A,1X,I6,7(1X,ES16.8),1X,L1)') &
               " OT Schur LS projector response n ratio ema pred meas weight wema neg enabled:", &
               local_ot_env(1)%schur_projector_response_quality_samples, &
               local_ot_env(1)%schur_projector_response_quality, &
               local_ot_env(1)%schur_projector_response_quality_ema, &
               local_ot_env(1)%schur_projector_response_drop_sum, &
               local_ot_env(1)%schur_projector_response_measured_sum, &
               local_ot_env(1)%schur_projector_response_weight, &
                  local_ot_env(1)%schur_projector_response_weight_ema, &
                  local_ot_env(1)%schur_projector_response_negative_ema, &
                  local_ot_env(1)%schur_projector_response_candidate_enabled
               WRITE (cp_logger_get_default_unit_nr(logger), '(A,1X,I6,7(1X,ES16.8),1X,L1)') &
                  " OT Schur LS projector physical n ratio ema pred meas weight wema neg enabled:", &
                  local_ot_env(1)%schur_projector_physical_quality_samples, &
                  local_ot_env(1)%schur_projector_physical_quality, &
                  local_ot_env(1)%schur_projector_physical_quality_ema, &
                  local_ot_env(1)%schur_projector_physical_drop_sum, &
                  local_ot_env(1)%schur_projector_physical_measured_sum, &
                  local_ot_env(1)%schur_projector_physical_weight, &
                  local_ot_env(1)%schur_projector_physical_weight_ema, &
                  local_ot_env(1)%schur_projector_physical_negative_ema, &
                  local_ot_env(1)%schur_projector_physical_candidate_enabled
            END IF

      END SUBROUTINE update_kpoint_ot_schur_line_search_quality

      SUBROUTINE precondition_kpoint_ot_energy_diagonal(skip_group)
      LOGICAL, DIMENSION(:), INTENT(IN)                  :: skip_group

      INTEGER                                            :: group, ikpoint, imo, ispin, &
                                                            local_channel, local_kpoint, nmo
      REAL(KIND=dp)                                      :: curvature, gvalue, maxocc, mu, &
                                                            rayleigh_edge_gap, tether_response, &
                                                            tether_strength, tether_weight
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: channel_diag
      TYPE(kpoint_env_type), POINTER                     :: kp

         precond_count(:) = 0.0_dp
         precond_shift(:) = 0.0_dp
         precond_sum(:) = 0.0_dp

         DO local_kpoint = 1, SIZE(kpoints%kp_env)
            kp => kpoints%kp_env(local_kpoint)%kpoint_env
            ikpoint = kp%nkpoint
            DO ispin = 1, nspin_ot
               group = MERGE(ispin, 1, ngroups == nspin_ot)
               IF (skip_group(group)) CYCLE
               local_channel = qs_ot_channel_index(ispin, ikpoint, nspin_ot) - first_channel + 1
               CALL get_mo_set(kp%mos(1, ispin), maxocc=maxocc, nmo=nmo)
               mu = trial_mu(ispin)
               tether_strength = 0.0_dp
               IF (local_ot_env(local_channel)%settings%canonical_tether) THEN
                  ALLOCATE (channel_diag(nmo))
                  CALL dbcsr_get_diag(local_ot_env(local_channel)%rot_mat_chc, channel_diag)
                  rayleigh_edge_gap = 0.0_dp
                  DO imo = 1, nmo
                     tether_response = kpoint_ot_canonical_tether_response(trial_occupation(imo, ikpoint, ispin), &
                                                                           trial_eigenvalues(imo, ikpoint, ispin), mu, &
                                                                           sigma, maxocc)
                     IF (tether_response > 0.0_dp) &
                        rayleigh_edge_gap = MAX(rayleigh_edge_gap, &
                                                ABS(channel_diag(imo) - local_ot_env(local_channel)%ener_x(imo)))
                  END DO
                  tether_strength = kpoint_ot_canonical_tether_strength(local_ot_env(local_channel), &
                                                                        rayleigh_edge_gap, sigma)
                  DEALLOCATE (channel_diag)
               END IF
               DO imo = 1, nmo
                  gvalue = kpoint_ot_smearing_g(trial_occupation(imo, ikpoint, ispin), &
                                                trial_eigenvalues(imo, ikpoint, ispin), mu, &
                                                sigma, maxocc, smear%method)
                  tether_weight = kpoint_ot_canonical_tether_weight(trial_occupation(imo, ikpoint, ispin), &
                                                                     trial_eigenvalues(imo, ikpoint, ispin), &
                                                                     mu, sigma, &
                                                                     maxocc, wkp(ikpoint), tether_strength)
                  curvature = wkp(ikpoint)*gvalue + tether_weight
                  IF (curvature > EPSILON(curvature)) THEN
                     local_ot_env(local_channel)%ener_gx(imo) = &
                        local_ot_env(local_channel)%ener_gx(imo)/MAX(curvature, curvature_floor)
                     precond_sum(group) = precond_sum(group) + local_ot_env(local_channel)%ener_gx(imo)
                     precond_count(group) = precond_count(group) + 1.0_dp
                  ELSE
                     local_ot_env(local_channel)%ener_gx(imo) = 0.0_dp
                  END IF
               END DO
            END DO
         END DO

         CALL para_env_inter_kp%sum(precond_sum)
         CALL para_env_inter_kp%sum(precond_count)
         DO group = 1, ngroups
            IF (skip_group(group)) CYCLE
            IF (precond_count(group) > 0.0_dp) precond_shift(group) = precond_sum(group)/precond_count(group)
         END DO

         DO local_kpoint = 1, SIZE(kpoints%kp_env)
            kp => kpoints%kp_env(local_kpoint)%kpoint_env
            ikpoint = kp%nkpoint
            DO ispin = 1, nspin_ot
               group = MERGE(ispin, 1, ngroups == nspin_ot)
               IF (skip_group(group)) CYCLE
               local_channel = qs_ot_channel_index(ispin, ikpoint, nspin_ot) - first_channel + 1
               CALL get_mo_set(kp%mos(1, ispin), maxocc=maxocc, nmo=nmo)
               mu = trial_mu(ispin)
               tether_strength = 0.0_dp
               IF (local_ot_env(local_channel)%settings%canonical_tether) THEN
                  ALLOCATE (channel_diag(nmo))
                  CALL dbcsr_get_diag(local_ot_env(local_channel)%rot_mat_chc, channel_diag)
                  rayleigh_edge_gap = 0.0_dp
                  DO imo = 1, nmo
                     tether_response = kpoint_ot_canonical_tether_response(trial_occupation(imo, ikpoint, ispin), &
                                                                           trial_eigenvalues(imo, ikpoint, ispin), mu, &
                                                                           sigma, maxocc)
                     IF (tether_response > 0.0_dp) &
                        rayleigh_edge_gap = MAX(rayleigh_edge_gap, &
                                                ABS(channel_diag(imo) - local_ot_env(local_channel)%ener_x(imo)))
                  END DO
                  tether_strength = kpoint_ot_canonical_tether_strength(local_ot_env(local_channel), &
                                                                        rayleigh_edge_gap, sigma)
                  DEALLOCATE (channel_diag)
               END IF
               DO imo = 1, nmo
                  gvalue = kpoint_ot_smearing_g(trial_occupation(imo, ikpoint, ispin), &
                                                trial_eigenvalues(imo, ikpoint, ispin), mu, &
                                                sigma, maxocc, smear%method)
                  tether_weight = kpoint_ot_canonical_tether_weight(trial_occupation(imo, ikpoint, ispin), &
                                                                     trial_eigenvalues(imo, ikpoint, ispin), &
                                                                     mu, sigma, &
                                                                     maxocc, wkp(ikpoint), tether_strength)
                  curvature = wkp(ikpoint)*gvalue + tether_weight
                  IF (curvature > EPSILON(curvature)) &
                     local_ot_env(local_channel)%ener_gx(imo) = &
                     local_ot_env(local_channel)%ener_gx(imo) - precond_shift(group)
               END DO
            END DO
         END DO
      END SUBROUTINE precondition_kpoint_ot_energy_diagonal

      SUBROUTINE precondition_kpoint_ot_schur(group_done)
      LOGICAL, DIMENSION(:), INTENT(OUT)                 :: group_done

      INTEGER                                            :: group

         group_done(:) = .FALSE.
         IF (.NOT. local_ot_env(1)%settings%do_rotation) RETURN
         CALL update_kpoint_ot_schur_line_search_quality()
         IF (.NOT. local_ot_env(1)%energy_only .AND. &
             (local_ot_env(1)%line_search_count == 0 .OR. local_ot_env(1)%line_search_might_be_done)) THEN
            local_ot_env(1)%schur_model_pending = .FALSE.
            local_ot_env(1)%schur_constrained_pending = .FALSE.
            local_ot_env(1)%schur_split_shadow_pending = .FALSE.
               local_ot_env(1)%schur_occgap_shadow_pending = .FALSE.
               local_ot_env(1)%schur_projector_mode_pending = .FALSE.
               local_ot_env(1)%schur_projector_response_pending = .FALSE.
               local_ot_env(1)%schur_projector_physical_pending = .FALSE.
               local_ot_env(1)%schur_history_mixed = .FALSE.
            local_ot_env(1)%schur_model_attempt = 0
            local_ot_env(1)%schur_model_drop = 0.0_dp
            local_ot_env(1)%schur_model_calibrated_drop = 0.0_dp
            local_ot_env(1)%schur_model_gap_ratio = 0.0_dp
            local_ot_env(1)%schur_model_delta = 0.0_dp
            local_ot_env(1)%schur_model_residual_drop = 0.0_dp
            local_ot_env(1)%schur_constrained_drop = 0.0_dp
            local_ot_env(1)%schur_constrained_calibrated_drop = 0.0_dp
            local_ot_env(1)%schur_constrained_gap_ratio = 0.0_dp
            local_ot_env(1)%schur_split_shadow_drop = 0.0_dp
            local_ot_env(1)%schur_split_shadow_calibrated_drop = 0.0_dp
            local_ot_env(1)%schur_occgap_shadow_drop = 0.0_dp
            local_ot_env(1)%schur_occgap_shadow_gap_ratio = 0.0_dp
            local_ot_env(1)%schur_projector_mode_drop = 0.0_dp
            local_ot_env(1)%schur_projector_mode_response_drop = 0.0_dp
            local_ot_env(1)%schur_projector_mode_step_negative = 0.0_dp
            local_ot_env(1)%schur_projector_mode_response_negative = 0.0_dp
            local_ot_env(1)%schur_projector_mode_trust = 0.0_dp
            local_ot_env(1)%schur_projector_response_drop = 0.0_dp
               local_ot_env(1)%schur_projector_response_weight = 0.0_dp
               local_ot_env(1)%schur_projector_response_correction_drop = 0.0_dp
               local_ot_env(1)%schur_projector_response_negative = 0.0_dp
               local_ot_env(1)%schur_projector_physical_drop = 0.0_dp
               local_ot_env(1)%schur_projector_physical_weight = 0.0_dp
               local_ot_env(1)%schur_projector_physical_correction_drop = 0.0_dp
               local_ot_env(1)%schur_projector_physical_negative = 0.0_dp
            END IF
         DO group = 1, ngroups
            CALL precondition_kpoint_ot_schur_group(group, group_done(group))
         END DO
      END SUBROUTINE precondition_kpoint_ot_schur

      SUBROUTINE precondition_kpoint_ot_schur_group(target_group, group_done)
      INTEGER, INTENT(IN)                                :: target_group
      LOGICAL, INTENT(OUT)                               :: group_done

      INTEGER, PARAMETER                                 :: max_schur_rot = 800
      INTEGER, PARAMETER                                 :: reject_dgesv = 1, reject_model_mismatch = 8, &
                                                            reject_nondescent = 2, reject_norm = 3, &
                                                            reject_response = 4, reject_scale = 5, &
                                                            reject_model_gap = 6, reject_quality_gate = 7
      REAL(KIND=dp), PARAMETER                           :: max_norm_ratio = 1.0E6_dp, &
                                                            schur_frontier_curvature_floor = 5.0E-1_dp, &
                                                            schur_min_coupled_scale = 5.0E-1_dp, &
                                                            schur_model_gap_factor = 5.0E-1_dp, &
                                                            schur_response_penalty_scale = 1.0_dp, &
                                                            schur_response_trust_factor = 2.5E-1_dp

      INTEGER                                            :: channel_var_start, group, ie, ikpoint, &
                                                            imo, info, ispin, ivar, jmo, &
                                                            local_channel, local_kpoint, &
                                                            nchannels, ne, nmo, nrot, r, s, schur_attempt, &
                                                            selected_attempt, split_selected_attempt, &
                                                            schur_used_attempt, constrained_reject_reason, &
                                                            coupled_reject_reason, dense_reject_reason, &
                                                            reject_reason
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: channel_nmo, e_offset, ipiv, ipiv_energy, &
                                                            var_channel, var_ei, var_ej, var_i, &
                                                            var_j, var_kind
      LOGICAL                                            :: attempt_failed, exact_attempt_response, failed, &
                                                            new_search_point, response_model_available, solved, &
                                                            constrained_candidate_enabled, &
                                                            constrained_candidate_ready, constrained_candidate_selected, &
                                                               exact_energy_model_ready, exact_physical_response_ready, &
                                                               projector_physical_candidate_ready, &
                                                               projector_physical_candidate_selected, &
                                                               projector_physical_direct_candidate_ready, &
                                                               projector_physical_direct_candidate_selected, &
                                                               projector_physical_direct_candidate_used, &
                                                               projector_physical_shadow_enabled, &
                                                               projector_response_candidate_ready, &
                                                               projector_response_candidate_selected, &
                                                            split_candidate_ready, split_candidate_selected, &
                                                            split_candidate_diagnostic, split_candidate_enabled, &
                                                            trust_rejected
      REAL(KIND=dp) :: accepted_model_calibrated_drop, accepted_model_drop, accepted_model_gap_ratio, &
         accepted_occgap_shadow_drop, accepted_occgap_shadow_gap_ratio, &
         constrained_diag_shift_avg, constrained_direction_scale, constrained_model_calibrated_drop, &
         constrained_model_drop, constrained_model_gap, constrained_model_gap_drop, constrained_model_gap_ratio, &
         constrained_expected_drop, constrained_expected_drop_ratio, constrained_gdotp, &
         constrained_occgap_shadow_drop, constrained_occgap_shadow_gap_ratio, &
         constrained_model_gap_trust, constrained_normp, constrained_norm_ratio, constrained_occ_response_diag_avg, &
         constrained_response_abs_step, constrained_response_abs_trust_ratio, constrained_response_penalty_diag_avg, &
         constrained_response_step, constrained_response_trust_ratio, constrained_response_weighted_step, &
         constrained_response_weighted_trust_ratio, constrained_rot_diag_avg, constrained_tether_response_diag_avg, &
         constrained_trust_quality, coupled_diag_shift_avg, coupled_direction_scale, coupled_model_calibrated_drop, &
         coupled_model_drop, &
         coupled_model_gap, coupled_model_gap_drop, coupled_occ_response_diag_avg, coupled_rot_diag_avg, &
         coupled_tether_response_diag_avg, &
         coupled_model_gap_ratio, coupled_model_gap_trust, coupled_norm_ratio, &
         coupled_response_abs_step, coupled_response_abs_trust_ratio, coupled_response_step, &
         coupled_response_trust_ratio, coupled_response_weighted_step, coupled_response_weighted_trust_ratio, &
         coupled_gdot_energy, coupled_gdot_rotation, coupled_norm_energy, coupled_norm_rotation, &
         coupled_scale_quality, curv_projection_sum, dense_direction_scale, dense_expected_drop, dense_term, &
         dense_gdot_energy, dense_gdot_rotation, dense_model_gap, dense_model_gap_drop, dense_model_gap_ratio, &
         dense_model_gap_trust, dense_norm_energy, dense_norm_rotation, dense_quality_for_selection, &
         dense_residual_quality, dense_scale_quality, dense_trust_quality, diag_shift, diag_shift_avg, direction_scale, &
         energy_gauge_weight, exact_response_diag_avg, gdotp, gdot_energy, &
         gdot_rotation, gnorm, &
         accepted_projector_coupled_weight, density_projector_fraction, dense_projector_coupled_weight, diag_shift_max, &
         diag_shift_sum, gvalue, maxocc, model_calibrated_drop, model_drop, model_gap, model_gap_drop, &
         model_gap_trust, &
         model_gap_ratio, mu, norm_energy, norm_ratio, norm_rotation, normp, occ_curv_projection_sum, &
         occ_projected_term, occgap_shadow_abs_step, occgap_shadow_drop, occgap_shadow_gap_ratio, &
         occgap_shadow_quality, occgap_shadow_step, &
         occ_response_diag_sum, projected_term, projector_coupled_weight, rayleigh_edge_gap, residual_delta, &
         physical_response_gap_avg, physical_response_gap_max, physical_response_gap_norm, &
         projector_hessian_abs_avg, projector_hessian_delta_diag_avg, projector_hessian_delta_max, &
         projector_hessian_diag_avg, projector_physical_abs_avg, projector_physical_delta_diag_avg, &
         projector_physical_delta_max, projector_physical_diag_avg, projector_gauge_abs_avg, &
         projector_gauge_delta_max, projector_diag_compensation, projector_mode_damped_abs_avg, &
         projector_mode_eval_max, projector_mode_eval_min, projector_mode_negative_fraction, projector_mode_quality, &
         accepted_projector_mode_drop, accepted_projector_mode_response_drop, &
         accepted_projector_mode_response_negative_fraction, accepted_projector_mode_step_negative_fraction, &
         accepted_projector_mode_trust_quality, projector_mode_history_quality, &
            projector_mode_history_response_negative, projector_mode_history_response_ratio, &
            accepted_projector_response_drop, accepted_projector_response_negative, &
            accepted_projector_response_weight, &
            projector_physical_direction_scale, projector_physical_expected_drop, &
            projector_physical_gdot_energy, projector_physical_gdot_rotation, projector_physical_gdotp, &
            projector_physical_model_calibrated_drop, projector_physical_model_drop, &
            projector_physical_model_gap, projector_physical_model_gap_drop, projector_physical_model_gap_ratio, &
            projector_physical_model_gap_trust, projector_physical_norm_energy, projector_physical_norm_rotation, &
            projector_physical_normp, projector_physical_norm_ratio, projector_physical_quality_for_selection, &
            projector_physical_direct_diag_compensation, projector_physical_direct_diag_shift_avg, &
            projector_physical_direct_diag_shift_sum, projector_physical_direct_direction_scale, &
            projector_physical_direct_expected_drop, &
            projector_physical_direct_gdot_energy, projector_physical_direct_gdot_rotation, &
            projector_physical_direct_gdotp, projector_physical_direct_model_calibrated_drop, &
            projector_physical_direct_model_drop, projector_physical_direct_model_gap, &
            projector_physical_direct_model_gap_drop, projector_physical_direct_model_gap_ratio, &
            projector_physical_direct_model_gap_trust, projector_physical_direct_norm_energy, &
            projector_physical_direct_norm_rotation, projector_physical_direct_normp, &
            projector_physical_direct_norm_ratio, projector_physical_direct_response_abs_step, &
            projector_physical_direct_response_abs_trust_ratio, &
            projector_physical_direct_response_weighted_step, &
            projector_physical_direct_response_weighted_trust_ratio, projector_physical_direct_scale_quality, &
            projector_physical_direct_trust_quality, &
            projector_physical_response_abs_step, projector_physical_response_abs_trust_ratio, &
            projector_physical_response_weighted_step, projector_physical_response_weighted_trust_ratio, &
            projector_physical_scale_quality, projector_physical_trust_quality, &
            projector_physical_response_correction_drop, projector_physical_response_correction_negative, &
            projector_physical_response_negative_coeff, projector_physical_response_positive_coeff, &
            projector_physical_response_cross_drop, projector_physical_response_negative_drop, &
            projector_physical_response_positive_drop, projector_physical_response_raw_negative_fraction, &
            projector_physical_response_weight, projector_response_direction_scale, &
            projector_response_expected_drop, projector_response_gdot_energy, projector_response_gdot_rotation, &
         projector_response_gdotp, projector_response_model_calibrated_drop, projector_response_model_drop, &
         projector_response_model_gap, projector_response_model_gap_drop, projector_response_model_gap_ratio, &
         projector_response_model_gap_trust, projector_response_norm_energy, projector_response_norm_rotation, &
         projector_response_normp, projector_response_norm_ratio, projector_response_quality_for_selection, &
         projector_response_response_abs_step, projector_response_response_abs_trust_ratio, &
         projector_response_response_weighted_step, projector_response_response_weighted_trust_ratio, &
         projector_response_scale_quality, projector_response_trust_quality, &
         projector_mode_response_correction_drop, projector_mode_response_correction_negative, &
         projector_mode_response_negative_coeff, projector_mode_response_positive_coeff, &
         projector_mode_response_cross_drop, projector_mode_response_negative_drop, &
         projector_mode_response_positive_drop, &
         projector_mode_response_drop_ratio, projector_mode_response_weight, &
         projector_mode_response_drop, projector_mode_response_negative_fraction, &
         projector_mode_response_raw_negative_fraction, &
         projector_mode_step_negative_fraction, projector_mode_trust_quality, &
         response_penalty_diag_sum, &
         response_penalty_term, response_shift, rot_diag_sum, projector_response_gap_avg, &
         projector_response_gap_max, projector_response_gap_norm, &
         response_mixed_diag_avg, &
         response_abs_step, response_abs_trust, response_abs_trust_ratio, response_step_max, response_trust, &
         response_trust_base, response_trust_factor, response_trust_ratio, response_weighted_step, &
         response_weighted_trust_ratio, row_abs, row_floor, rot_floor_min, rx_shift, rx_weight, scale_quality, &
         schur_diag_min, split_density_projector_fraction, split_diag_shift_max, split_diag_shift_sum, &
         split_direction_scale, &
         split_gdotp, split_gdot_energy, split_gdot_rotation, split_model_calibrated_drop, split_model_drop, &
         split_model_gap_ratio, split_norm_energy, split_norm_rotation, split_normp, &
         split_norm_ratio, split_occgap_shadow_drop, split_occgap_shadow_gap_ratio, &
         split_residual_delta, split_response_abs_step, split_response_abs_trust, &
         split_response_abs_trust_ratio, split_response_step, split_response_trust, split_response_trust_factor, &
         split_response_weighted_step, split_response_weighted_trust_ratio, split_rot_floor_min, split_schur_diag_min, &
         split_scale_quality, split_curvature_scale, split_expected_drop, split_quality_for_calibration, &
         split_quality_for_selection, split_dense_expected_ratio, split_residual_quality, &
         residual_reference_quality, constrained_quality_for_selection, &
         constrained_gdot_energy, constrained_gdot_rotation, constrained_norm_energy, constrained_norm_rotation, &
         constrained_scale_quality, split_trust_quality, tether_projected_term, tether_response, tether_response_diag_sum, &
         tether_strength, tether_weight, z_shift, z_weight
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: curv, grad_e, grad_rot, occ_curv, rhs, &
                                                            constrained_sol, constrained_z_e, &
                                                            projector_physical_direct_model_e, &
                                                            projector_physical_direct_rx_e, &
                                                            projector_physical_direct_sol, &
                                                            projector_physical_direct_z_e, &
                                                            projector_physical_model_e, projector_physical_rx_e, &
                                                            projector_physical_sol, projector_physical_z_e, &
                                                            projector_response_rx_e, projector_response_sol, &
                                                            projector_response_z_e, rot_floor, rtd, rx_e, rx_model, &
                                                            sol, split_sol, split_z_e, tether_curv, z_e, &
                                                            z_exact_raw_trial, z_exact_physical_trial, &
                                                            z_exact_trial, z_trial
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: grad_im, grad_re, h_im, h_re, response, &
                                                            energy_cross, energy_cross_physical, energy_cross_raw, &
                                                            energy_hessian, response_coupled, response_exact, &
                                                            response_exact_physical, response_exact_raw, response_occ, &
                                                            response_occ_physical, response_split, rot_hessian, &
                                                            rot_hessian_coupled, rot_hessian_projector, &
                                                            rot_hessian_projector_physical, schur
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(kpoint_env_type), POINTER                     :: kp

         group_done = .FALSE.
         logger => cp_get_default_logger()
         IF (ABS(gsum(target_group)) <= EPSILON(gsum(target_group))) RETURN

            new_search_point = .NOT. local_ot_env(1)%energy_only .AND. &
                               (local_ot_env(1)%line_search_count == 0 .OR. &
                                local_ot_env(1)%line_search_might_be_done)
            projector_physical_shadow_enabled = logger%iter_info%print_level >= high_print_level .OR. &
                                                local_ot_env(1)%schur_projector_physical_candidate_enabled
            IF (new_search_point) &
               projector_physical_shadow_enabled = projector_physical_shadow_enabled .OR. &
                                                   local_ot_env(1)%schur_projector_physical_quality_samples > 0 .OR. &
                                                   local_ot_env(1)%schur_projector_mode_quality_samples >= &
                                                   ot_kp_schur_calibration_min_samples

            nchannels = SIZE(local_ot_env)
         ALLOCATE (channel_nmo(nchannels), e_offset(nchannels))
         channel_nmo(:) = 0
         e_offset(:) = 0

         ne = 0
         nrot = 0
         DO local_kpoint = 1, SIZE(kpoints%kp_env)
            kp => kpoints%kp_env(local_kpoint)%kpoint_env
            ikpoint = kp%nkpoint
            DO ispin = 1, nspin_ot
               group = MERGE(ispin, 1, ngroups == nspin_ot)
               IF (group /= target_group) CYCLE
               local_channel = qs_ot_channel_index(ispin, ikpoint, nspin_ot) - first_channel + 1
               CALL get_mo_set(kp%mos(1, ispin), nmo=nmo)
               e_offset(local_channel) = ne + 1
               channel_nmo(local_channel) = nmo
               ne = ne + nmo
               nrot = nrot + nmo*(nmo - 1)
            END DO
         END DO

         IF (ne == 0 .OR. nrot == 0 .OR. nrot > max_schur_rot) THEN
            DEALLOCATE (channel_nmo, e_offset)
            RETURN
         END IF

         ALLOCATE (curv(ne), grad_e(ne), occ_curv(ne), tether_curv(ne), z_e(ne), z_exact_raw_trial(ne), &
                      z_exact_physical_trial(ne), z_exact_trial(ne), z_trial(ne), rx_e(ne), rx_model(ne), &
                      constrained_z_e(ne), projector_physical_direct_model_e(ne), &
                      projector_physical_direct_rx_e(ne), projector_physical_direct_z_e(ne), &
                      projector_physical_model_e(ne), projector_physical_rx_e(ne), &
                      projector_physical_z_e(ne), projector_response_rx_e(ne), &
                      projector_response_z_e(ne), split_z_e(ne))
            ALLOCATE (grad_rot(nrot), rhs(nrot), rtd(nrot), sol(nrot), constrained_sol(nrot), &
                      projector_physical_direct_sol(nrot), projector_physical_sol(nrot), &
                      projector_response_sol(nrot), split_sol(nrot))
         ALLOCATE (var_channel(nrot), var_ei(nrot), var_ej(nrot), var_i(nrot), var_j(nrot), var_kind(nrot))
         ALLOCATE (rot_floor(nrot))
         ALLOCATE (energy_cross(ne, nrot), energy_cross_physical(ne, nrot), energy_cross_raw(ne, nrot), &
                   energy_hessian(ne, ne), response(ne, nrot), response_coupled(ne, nrot), &
                   response_exact(ne, nrot), response_exact_physical(ne, nrot), &
                   response_exact_raw(ne, nrot), response_occ(ne, nrot), response_occ_physical(ne, nrot), &
                   response_split(ne, nrot), &
                   rot_hessian(nrot, nrot), rot_hessian_coupled(nrot, nrot), &
                   rot_hessian_projector(nrot, nrot), rot_hessian_projector_physical(nrot, nrot), &
                   schur(nrot, nrot), ipiv(nrot))
         ALLOCATE (ipiv_energy(ne))

         curv(:) = 0.0_dp
         grad_e(:) = 0.0_dp
         occ_curv(:) = 0.0_dp
         tether_curv(:) = 0.0_dp
         z_e(:) = 0.0_dp
         z_exact_physical_trial(:) = 0.0_dp
         z_exact_raw_trial(:) = 0.0_dp
         z_exact_trial(:) = 0.0_dp
            rx_e(:) = 0.0_dp
            rx_model(:) = 0.0_dp
            projector_physical_direct_model_e(:) = 0.0_dp
            projector_physical_direct_rx_e(:) = 0.0_dp
            projector_physical_direct_z_e(:) = 0.0_dp
            projector_physical_model_e(:) = 0.0_dp
            projector_physical_rx_e(:) = 0.0_dp
            projector_physical_z_e(:) = 0.0_dp
            projector_response_rx_e(:) = 0.0_dp
            projector_response_z_e(:) = 0.0_dp
         grad_rot(:) = 0.0_dp
         rhs(:) = 0.0_dp
         rtd(:) = 0.0_dp
            sol(:) = 0.0_dp
            projector_physical_direct_sol(:) = 0.0_dp
            projector_physical_sol(:) = 0.0_dp
            projector_response_sol(:) = 0.0_dp
         energy_cross(:, :) = 0.0_dp
         energy_cross_physical(:, :) = 0.0_dp
         energy_cross_raw(:, :) = 0.0_dp
         energy_hessian(:, :) = 0.0_dp
         response(:, :) = 0.0_dp
         response_coupled(:, :) = 0.0_dp
         response_exact(:, :) = 0.0_dp
         response_exact_physical(:, :) = 0.0_dp
         response_exact_raw(:, :) = 0.0_dp
         response_occ(:, :) = 0.0_dp
         response_occ_physical(:, :) = 0.0_dp
         response_split(:, :) = 0.0_dp
         rot_hessian(:, :) = 0.0_dp
         rot_hessian_coupled(:, :) = 0.0_dp
         rot_hessian_projector(:, :) = 0.0_dp
         rot_hessian_projector_physical(:, :) = 0.0_dp
         rot_floor(:) = curvature_floor
         diag_shift_max = 0.0_dp
         diag_shift_sum = 0.0_dp
         response_abs_step = 0.0_dp
         response_abs_trust = HUGE(1.0_dp)
         response_step_max = 0.0_dp
         response_trust = HUGE(1.0_dp)
         response_trust_base = HUGE(1.0_dp)
         response_trust_factor = 1.0_dp
         residual_delta = MAX(local_ot_env(1)%delta, 0.0_dp)
         residual_reference_quality = 0.0_dp
         density_projector_fraction = 0.0_dp
         physical_response_gap_avg = 0.0_dp
         physical_response_gap_max = 0.0_dp
         rot_floor_min = curvature_floor
         schur_diag_min = HUGE(1.0_dp)
         constrained_diag_shift_avg = 0.0_dp
         constrained_direction_scale = 0.0_dp
         constrained_expected_drop = 0.0_dp
         constrained_expected_drop_ratio = 0.0_dp
         constrained_gdotp = 0.0_dp
         constrained_gdot_energy = 0.0_dp
         constrained_gdot_rotation = 0.0_dp
         constrained_norm_energy = 0.0_dp
         constrained_norm_rotation = 0.0_dp
         constrained_model_drop = 0.0_dp
         constrained_normp = 0.0_dp
         constrained_norm_ratio = 0.0_dp
         constrained_scale_quality = 0.0_dp
         constrained_reject_reason = 0
         constrained_model_gap = 0.0_dp
         constrained_model_gap_drop = 0.0_dp
         constrained_model_gap_ratio = 0.0_dp
         constrained_occgap_shadow_drop = 0.0_dp
         constrained_occgap_shadow_gap_ratio = 0.0_dp
         constrained_model_gap_trust = HUGE(1.0_dp)
         constrained_model_calibrated_drop = 0.0_dp
         constrained_occ_response_diag_avg = 0.0_dp
         constrained_response_penalty_diag_avg = 0.0_dp
         constrained_rot_diag_avg = 0.0_dp
         constrained_tether_response_diag_avg = 0.0_dp
         constrained_response_abs_step = 0.0_dp
         constrained_response_abs_trust_ratio = 0.0_dp
         constrained_response_step = 0.0_dp
         constrained_response_trust_ratio = 0.0_dp
         constrained_response_weighted_step = 0.0_dp
         constrained_response_weighted_trust_ratio = 0.0_dp
         constrained_trust_quality = 1.0_dp
         coupled_diag_shift_avg = 0.0_dp
         accepted_model_drop = 0.0_dp
         accepted_model_calibrated_drop = 0.0_dp
         accepted_model_gap_ratio = 0.0_dp
         accepted_occgap_shadow_drop = 0.0_dp
         accepted_occgap_shadow_gap_ratio = 0.0_dp
         accepted_projector_coupled_weight = 0.0_dp
         accepted_projector_mode_drop = 0.0_dp
         accepted_projector_mode_response_drop = 0.0_dp
         accepted_projector_mode_response_negative_fraction = 0.0_dp
         accepted_projector_mode_step_negative_fraction = 0.0_dp
         accepted_projector_mode_trust_quality = 1.0_dp
         accepted_projector_response_drop = 0.0_dp
         accepted_projector_response_negative = 0.0_dp
         accepted_projector_response_weight = 0.0_dp
         coupled_direction_scale = 0.0_dp
         coupled_gdot_energy = 0.0_dp
         coupled_gdot_rotation = 0.0_dp
         coupled_model_drop = 0.0_dp
         coupled_norm_energy = 0.0_dp
         coupled_norm_rotation = 0.0_dp
         coupled_norm_ratio = 0.0_dp
         coupled_scale_quality = 0.0_dp
         coupled_reject_reason = 0
         coupled_model_gap = 0.0_dp
         coupled_model_gap_drop = 0.0_dp
         coupled_model_gap_ratio = 0.0_dp
         coupled_model_gap_trust = HUGE(1.0_dp)
         coupled_model_calibrated_drop = 0.0_dp
         coupled_occ_response_diag_avg = 0.0_dp
         coupled_rot_diag_avg = 0.0_dp
         coupled_tether_response_diag_avg = 0.0_dp
         coupled_response_abs_step = 0.0_dp
         coupled_response_abs_trust_ratio = 0.0_dp
         coupled_response_step = 0.0_dp
         coupled_response_trust_ratio = 0.0_dp
         coupled_response_weighted_step = 0.0_dp
         coupled_response_weighted_trust_ratio = 0.0_dp
         dense_direction_scale = 0.0_dp
         dense_expected_drop = 0.0_dp
         dense_gdot_energy = 0.0_dp
         dense_gdot_rotation = 0.0_dp
         dense_model_gap = 0.0_dp
         dense_model_gap_drop = 0.0_dp
         dense_model_gap_ratio = 0.0_dp
         dense_model_gap_trust = 0.0_dp
         dense_norm_energy = 0.0_dp
         dense_norm_rotation = 0.0_dp
         dense_projector_coupled_weight = 0.0_dp
         dense_quality_for_selection = 1.0_dp
         dense_reject_reason = 0
         dense_residual_quality = 1.0_dp
         dense_scale_quality = 0.0_dp
         dense_trust_quality = 1.0_dp
         energy_gauge_weight = 0.0_dp
         exact_energy_model_ready = .FALSE.
         exact_response_diag_avg = 0.0_dp
         response_mixed_diag_avg = 0.0_dp
         projector_hessian_abs_avg = 0.0_dp
         projector_hessian_delta_diag_avg = 0.0_dp
         projector_hessian_delta_max = 0.0_dp
         projector_hessian_diag_avg = 0.0_dp
         projector_physical_abs_avg = 0.0_dp
         projector_physical_delta_diag_avg = 0.0_dp
         projector_physical_delta_max = 0.0_dp
         projector_physical_diag_avg = 0.0_dp
         projector_gauge_abs_avg = 0.0_dp
         projector_gauge_delta_max = 0.0_dp
         projector_diag_compensation = 0.0_dp
         projector_mode_damped_abs_avg = 0.0_dp
         projector_mode_eval_max = 0.0_dp
         projector_mode_eval_min = 0.0_dp
            projector_mode_negative_fraction = 0.0_dp
            projector_mode_quality = 0.0_dp
            projector_mode_history_quality = 1.0_dp
            projector_mode_history_response_negative = 0.0_dp
            projector_mode_history_response_ratio = 0.0_dp
            projector_physical_direct_diag_compensation = 0.0_dp
            projector_physical_direct_diag_shift_avg = 0.0_dp
            projector_physical_direct_diag_shift_sum = 0.0_dp
            projector_physical_direct_direction_scale = 0.0_dp
            projector_physical_direct_expected_drop = 0.0_dp
            projector_physical_direct_gdot_energy = 0.0_dp
            projector_physical_direct_gdot_rotation = 0.0_dp
            projector_physical_direct_gdotp = 0.0_dp
            projector_physical_direct_model_calibrated_drop = 0.0_dp
            projector_physical_direct_model_drop = 0.0_dp
            projector_physical_direct_model_gap = 0.0_dp
            projector_physical_direct_model_gap_drop = 0.0_dp
            projector_physical_direct_model_gap_ratio = 0.0_dp
            projector_physical_direct_model_gap_trust = HUGE(1.0_dp)
            projector_physical_direct_norm_energy = 0.0_dp
            projector_physical_direct_norm_rotation = 0.0_dp
            projector_physical_direct_normp = 0.0_dp
            projector_physical_direct_norm_ratio = 0.0_dp
            projector_physical_direct_response_abs_step = 0.0_dp
            projector_physical_direct_response_abs_trust_ratio = 0.0_dp
            projector_physical_direct_response_weighted_step = 0.0_dp
            projector_physical_direct_response_weighted_trust_ratio = 0.0_dp
            projector_physical_direct_scale_quality = 0.0_dp
            projector_physical_direct_trust_quality = 1.0_dp
            projector_physical_direction_scale = 0.0_dp
            projector_physical_expected_drop = 0.0_dp
            projector_physical_gdot_energy = 0.0_dp
            projector_physical_gdot_rotation = 0.0_dp
            projector_physical_gdotp = 0.0_dp
            projector_physical_model_calibrated_drop = 0.0_dp
            projector_physical_model_drop = 0.0_dp
            projector_physical_model_gap = 0.0_dp
            projector_physical_model_gap_drop = 0.0_dp
            projector_physical_model_gap_ratio = 0.0_dp
            projector_physical_model_gap_trust = HUGE(1.0_dp)
            projector_physical_norm_energy = 0.0_dp
            projector_physical_norm_rotation = 0.0_dp
            projector_physical_normp = 0.0_dp
            projector_physical_norm_ratio = 0.0_dp
            projector_physical_quality_for_selection = 0.0_dp
            projector_physical_response_abs_step = 0.0_dp
            projector_physical_response_abs_trust_ratio = 0.0_dp
            projector_physical_response_weighted_step = 0.0_dp
            projector_physical_response_weighted_trust_ratio = 0.0_dp
            projector_physical_scale_quality = 0.0_dp
            projector_physical_trust_quality = 1.0_dp
            projector_physical_response_correction_drop = 0.0_dp
            projector_physical_response_correction_negative = 0.0_dp
            projector_physical_response_negative_coeff = 0.0_dp
            projector_physical_response_positive_coeff = 0.0_dp
            projector_physical_response_cross_drop = 0.0_dp
            projector_physical_response_negative_drop = 0.0_dp
            projector_physical_response_positive_drop = 0.0_dp
            projector_physical_response_raw_negative_fraction = 0.0_dp
            projector_physical_response_weight = 0.0_dp
            projector_response_direction_scale = 0.0_dp
         projector_response_expected_drop = 0.0_dp
         projector_response_gdot_energy = 0.0_dp
         projector_response_gdot_rotation = 0.0_dp
         projector_response_gdotp = 0.0_dp
         projector_response_model_calibrated_drop = 0.0_dp
         projector_response_model_drop = 0.0_dp
         projector_response_model_gap = 0.0_dp
         projector_response_model_gap_drop = 0.0_dp
         projector_response_model_gap_ratio = 0.0_dp
         projector_response_model_gap_trust = HUGE(1.0_dp)
         projector_response_norm_energy = 0.0_dp
         projector_response_norm_rotation = 0.0_dp
         projector_response_normp = 0.0_dp
         projector_response_norm_ratio = 0.0_dp
         projector_response_quality_for_selection = 0.0_dp
         projector_response_response_abs_step = 0.0_dp
         projector_response_response_abs_trust_ratio = 0.0_dp
         projector_response_response_weighted_step = 0.0_dp
         projector_response_response_weighted_trust_ratio = 0.0_dp
         projector_response_scale_quality = 0.0_dp
         projector_response_trust_quality = 1.0_dp
         projector_mode_response_correction_drop = 0.0_dp
         projector_mode_response_correction_negative = 0.0_dp
         projector_mode_response_negative_coeff = 0.0_dp
         projector_mode_response_positive_coeff = 0.0_dp
         projector_mode_response_cross_drop = 0.0_dp
         projector_mode_response_negative_drop = 0.0_dp
         projector_mode_response_positive_drop = 0.0_dp
         projector_mode_response_drop = 0.0_dp
         projector_mode_response_drop_ratio = 0.0_dp
         projector_mode_response_negative_fraction = 0.0_dp
         projector_mode_response_raw_negative_fraction = 0.0_dp
         projector_mode_step_negative_fraction = 0.0_dp
         projector_mode_trust_quality = 1.0_dp
         projector_mode_response_weight = 0.0_dp
         model_gap = 0.0_dp
         model_gap_drop = 0.0_dp
         model_gap_ratio = 0.0_dp
         model_gap_trust = HUGE(1.0_dp)
         model_calibrated_drop = 0.0_dp
         occgap_shadow_abs_step = 0.0_dp
         occgap_shadow_drop = 0.0_dp
         occgap_shadow_gap_ratio = 0.0_dp
         occgap_shadow_quality = 1.0_dp
         occgap_shadow_step = 0.0_dp
         split_curvature_scale = 1.0_dp
         split_direction_scale = 0.0_dp
         split_dense_expected_ratio = 0.0_dp
         split_expected_drop = 0.0_dp
         split_gdotp = 0.0_dp
         split_gdot_energy = 0.0_dp
         split_gdot_rotation = 0.0_dp
         split_model_calibrated_drop = 0.0_dp
         split_model_drop = 0.0_dp
         split_model_gap_ratio = 0.0_dp
         split_occgap_shadow_drop = 0.0_dp
         split_occgap_shadow_gap_ratio = 0.0_dp
         split_norm_energy = 0.0_dp
         split_norm_rotation = 0.0_dp
         split_normp = 0.0_dp
         split_norm_ratio = 0.0_dp
         split_scale_quality = 0.0_dp
         split_response_abs_step = 0.0_dp
         split_response_abs_trust_ratio = 0.0_dp
         split_response_step = 0.0_dp
         split_response_weighted_step = 0.0_dp
         split_response_weighted_trust_ratio = 0.0_dp
         split_quality_for_calibration = 1.0_dp
         split_quality_for_selection = ot_kp_schur_calibrated_quality_min_ratio
         split_residual_quality = 1.0_dp
         constrained_quality_for_selection = ot_kp_schur_calibrated_quality_min_ratio
         split_trust_quality = 1.0_dp
         split_density_projector_fraction = 0.0_dp
         split_diag_shift_max = 0.0_dp
         split_diag_shift_sum = 0.0_dp
         split_residual_delta = 0.0_dp
         split_response_abs_trust = HUGE(1.0_dp)
         split_response_trust = HUGE(1.0_dp)
         split_response_trust_factor = 1.0_dp
         split_rot_floor_min = curvature_floor
         split_schur_diag_min = HUGE(1.0_dp)
            exact_physical_response_ready = .FALSE.
            failed = .FALSE.
            constrained_candidate_ready = .FALSE.
            constrained_candidate_selected = .FALSE.
            projector_physical_candidate_ready = .FALSE.
            projector_physical_candidate_selected = .FALSE.
            projector_physical_direct_candidate_ready = .FALSE.
            projector_physical_direct_candidate_selected = .FALSE.
            projector_physical_direct_candidate_used = .FALSE.
            projector_response_candidate_ready = .FALSE.
            projector_response_candidate_selected = .FALSE.
         split_candidate_ready = .FALSE.
            split_candidate_selected = .FALSE.
            split_selected_attempt = 0
            ivar = 0
            constrained_sol(:) = 0.0_dp
            constrained_z_e(:) = 0.0_dp
            projector_physical_direct_sol(:) = 0.0_dp
            projector_physical_direct_z_e(:) = 0.0_dp
            projector_physical_sol(:) = 0.0_dp
            projector_physical_z_e(:) = 0.0_dp
            split_sol(:) = 0.0_dp
            split_z_e(:) = 0.0_dp

         DO local_kpoint = 1, SIZE(kpoints%kp_env)
            kp => kpoints%kp_env(local_kpoint)%kpoint_env
            ikpoint = kp%nkpoint
            DO ispin = 1, nspin_ot
               group = MERGE(ispin, 1, ngroups == nspin_ot)
               IF (group /= target_group) CYCLE
               local_channel = qs_ot_channel_index(ispin, ikpoint, nspin_ot) - first_channel + 1
               CALL get_mo_set(kp%mos(1, ispin), maxocc=maxocc, nmo=nmo)
               mu = trial_mu(ispin)
               ALLOCATE (h_re(nmo, nmo), h_im(nmo, nmo), grad_re(nmo, nmo), grad_im(nmo, nmo))
               CALL dbcsr_square_to_dense(local_ot_env(local_channel)%rot_mat_chc, h_re)
               CALL dbcsr_square_to_dense(local_ot_env(local_channel)%rot_mat_chc_im, h_im)
               CALL dbcsr_square_to_dense(local_ot_env(local_channel)%rot_mat_gx, grad_re)
               CALL dbcsr_square_to_dense(local_ot_env(local_channel)%rot_mat_gx_im, grad_im)
               CALL local_ot_env(local_channel)%para_env%sum(h_re)
               CALL local_ot_env(local_channel)%para_env%sum(h_im)
               CALL local_ot_env(local_channel)%para_env%sum(grad_re)
               CALL local_ot_env(local_channel)%para_env%sum(grad_im)

               tether_strength = 0.0_dp
               IF (local_ot_env(local_channel)%settings%canonical_tether) THEN
                  rayleigh_edge_gap = 0.0_dp
                  DO imo = 1, nmo
                     tether_response = kpoint_ot_canonical_tether_response(trial_occupation(imo, ikpoint, ispin), &
                                                                           trial_eigenvalues(imo, ikpoint, ispin), mu, &
                                                                           sigma, maxocc)
                     IF (tether_response > 0.0_dp) &
                        rayleigh_edge_gap = MAX(rayleigh_edge_gap, &
                                                ABS(h_re(imo, imo) - local_ot_env(local_channel)%ener_x(imo)))
                  END DO
                  tether_strength = kpoint_ot_canonical_tether_strength(local_ot_env(local_channel), &
                                                                        rayleigh_edge_gap, sigma)
               END IF

               DO imo = 1, nmo
                  ie = e_offset(local_channel) + imo - 1
                  gvalue = kpoint_ot_smearing_g(trial_occupation(imo, ikpoint, ispin), &
                                                trial_eigenvalues(imo, ikpoint, ispin), mu, &
                                                sigma, maxocc, smear%method)
                  tether_weight = kpoint_ot_canonical_tether_weight(trial_occupation(imo, ikpoint, ispin), &
                                                                     trial_eigenvalues(imo, ikpoint, ispin), &
                                                                     mu, sigma, &
                                                                     maxocc, wkp(ikpoint), tether_strength)
                  occ_curv(ie) = wkp(ikpoint)*gvalue
                  tether_curv(ie) = tether_weight
                  curv(ie) = occ_curv(ie) + tether_curv(ie)
                  grad_e(ie) = local_ot_env(local_channel)%ener_gx(imo)
               END DO

               channel_var_start = ivar + 1
               DO imo = 1, nmo - 1
                  DO jmo = imo + 1, nmo
                     ivar = ivar + 1
                     var_channel(ivar) = local_channel
                     var_i(ivar) = imo
                     var_j(ivar) = jmo
                     var_kind(ivar) = 1
                     var_ei(ivar) = e_offset(local_channel) + imo - 1
                     var_ej(ivar) = e_offset(local_channel) + jmo - 1
                     grad_rot(ivar) = 2.0_dp*grad_re(imo, jmo)
                     IF (MAX(curv(var_ei(ivar)), curv(var_ej(ivar))) > EPSILON(1.0_dp)) &
                        rot_floor(ivar) = schur_frontier_curvature_floor

                     ivar = ivar + 1
                     var_channel(ivar) = local_channel
                     var_i(ivar) = imo
                     var_j(ivar) = jmo
                     var_kind(ivar) = 2
                     var_ei(ivar) = e_offset(local_channel) + imo - 1
                     var_ej(ivar) = e_offset(local_channel) + jmo - 1
                     grad_rot(ivar) = 2.0_dp*grad_im(imo, jmo)
                     IF (MAX(curv(var_ei(ivar)), curv(var_ej(ivar))) > EPSILON(1.0_dp)) &
                        rot_floor(ivar) = schur_frontier_curvature_floor
                  END DO
               END DO
                  CALL build_finite_ref_rotation_model(local_channel, nmo, channel_var_start, ivar, &
                                                       e_offset(local_channel), h_re, h_im, &
                                                       trial_occupation(1:nmo, ikpoint, ispin), wkp(ikpoint), &
                                                       occ_curv(e_offset(local_channel):e_offset(local_channel) + nmo - 1), &
                                                       var_i, var_j, var_kind, &
                                                       projector_physical_shadow_enabled, &
                                                       rot_hessian, rot_hessian_coupled, rot_hessian_projector, &
                                                       rot_hessian_projector_physical, response, response_occ, &
                                                       response_occ_physical, failed)
               DEALLOCATE (h_re, h_im, grad_re, grad_im)
            END DO
         END DO

            CPASSERT(ivar == nrot)
            IF (failed) THEN
               DEALLOCATE (curv, grad_e, occ_curv, tether_curv, z_e, z_exact_physical_trial, z_exact_raw_trial, &
                           z_exact_trial, z_trial, rx_e, rx_model, constrained_z_e, &
                           projector_physical_direct_model_e, projector_physical_direct_rx_e, &
                           projector_physical_direct_z_e, projector_physical_model_e, &
                           projector_physical_rx_e, projector_physical_z_e, &
                           projector_response_rx_e, projector_response_z_e, split_z_e)
               DEALLOCATE (grad_rot, rhs, rtd, sol, constrained_sol, projector_physical_direct_sol, &
                           projector_physical_sol, projector_response_sol, split_sol)
            DEALLOCATE (var_channel, var_ei, var_ej, var_i, var_j, var_kind)
            DEALLOCATE (rot_floor, energy_cross, energy_cross_physical, energy_cross_raw, energy_hessian, &
                        response, response_coupled, response_exact, response_exact_physical, &
                        response_exact_raw, response_occ, response_occ_physical, response_split, rot_hessian, &
                        rot_hessian_coupled, rot_hessian_projector, rot_hessian_projector_physical, schur, &
                        ipiv, ipiv_energy, channel_nmo, e_offset)
            RETURN
         END IF

         z_shift = 0.0_dp
         z_weight = 0.0_dp
         DO ie = 1, ne
            IF (curv(ie) > EPSILON(curv(ie))) THEN
               z_e(ie) = grad_e(ie)/MAX(curv(ie), curvature_floor)
               z_shift = z_shift + curv(ie)*z_e(ie)
               z_weight = z_weight + curv(ie)
            END IF
         END DO
         CALL para_env_inter_kp%sum(z_shift)
         CALL para_env_inter_kp%sum(z_weight)
         IF (z_weight > EPSILON(z_weight)) z_e(:) = z_e(:) - z_shift/z_weight

         z_trial(:) = z_e(:)
         rot_hessian(:, :) = 0.5_dp*(rot_hessian(:, :) + TRANSPOSE(rot_hessian(:, :)))
         rot_hessian_coupled(:, :) = 0.5_dp*(rot_hessian_coupled(:, :) + &
                                             TRANSPOSE(rot_hessian_coupled(:, :)))
         rot_hessian_projector(:, :) = 0.5_dp*(rot_hessian_projector(:, :) + &
                                               TRANSPOSE(rot_hessian_projector(:, :)))
         rot_hessian_projector_physical(:, :) = 0.5_dp*(rot_hessian_projector_physical(:, :) + &
                                                        TRANSPOSE(rot_hessian_projector_physical(:, :)))
         DO r = 1, nrot
            projector_hessian_diag_avg = projector_hessian_diag_avg + rot_hessian_projector(r, r)
            projector_hessian_delta_diag_avg = projector_hessian_delta_diag_avg + &
                                               rot_hessian_projector(r, r) - rot_hessian(r, r)
            projector_physical_diag_avg = projector_physical_diag_avg + rot_hessian_projector_physical(r, r)
            projector_physical_delta_diag_avg = projector_physical_delta_diag_avg + &
                                                rot_hessian_projector_physical(r, r) - rot_hessian(r, r)
            DO s = 1, nrot
               projector_hessian_abs_avg = projector_hessian_abs_avg + &
                                           ABS(rot_hessian_projector(r, s) - rot_hessian(r, s))
               projector_hessian_delta_max = MAX(projector_hessian_delta_max, &
                                                 ABS(rot_hessian_projector(r, s) - rot_hessian(r, s)))
               projector_physical_abs_avg = projector_physical_abs_avg + &
                                            ABS(rot_hessian_projector_physical(r, s) - rot_hessian(r, s))
               projector_physical_delta_max = MAX(projector_physical_delta_max, &
                                                  ABS(rot_hessian_projector_physical(r, s) - rot_hessian(r, s)))
               projector_gauge_abs_avg = projector_gauge_abs_avg + &
                                         ABS(rot_hessian_projector(r, s) - rot_hessian_projector_physical(r, s))
               projector_gauge_delta_max = MAX(projector_gauge_delta_max, &
                                               ABS(rot_hessian_projector(r, s) - rot_hessian_projector_physical(r, s)))
            END DO
         END DO
         projector_hessian_diag_avg = projector_hessian_diag_avg/REAL(MAX(1, nrot), KIND=dp)
         projector_hessian_delta_diag_avg = projector_hessian_delta_diag_avg/REAL(MAX(1, nrot), KIND=dp)
         projector_hessian_abs_avg = projector_hessian_abs_avg/REAL(MAX(1, nrot*nrot), KIND=dp)
         projector_physical_diag_avg = projector_physical_diag_avg/REAL(MAX(1, nrot), KIND=dp)
         projector_physical_delta_diag_avg = projector_physical_delta_diag_avg/REAL(MAX(1, nrot), KIND=dp)
         projector_physical_abs_avg = projector_physical_abs_avg/REAL(MAX(1, nrot*nrot), KIND=dp)
         projector_gauge_abs_avg = projector_gauge_abs_avg/REAL(MAX(1, nrot*nrot), KIND=dp)
         curv_projection_sum = SUM(curv)
         CALL para_env_inter_kp%sum(curv_projection_sum)
         occ_curv_projection_sum = SUM(occ_curv)
         CALL para_env_inter_kp%sum(occ_curv_projection_sum)
         ! The split candidate follows the same diagonal/spectral projector blend as
         ! the finite complex rotation gradient: occupation curvature sees the
         ! finite projector response, while the channel-local tether still sees the
         ! full Rayleigh diagonal it anchors.
         response_coupled(:, :) = response(:, :)
         response_split(:, :) = 0.0_dp
         DO r = 1, nrot
            DO ie = 1, ne
               IF (curv(ie) > EPSILON(curv(ie))) &
                  response_split(ie, r) = (occ_curv(ie)*response_occ(ie, r) + &
                                           tether_curv(ie)*response(ie, r))/curv(ie)
            END DO
         END DO
         DO r = 1, nrot
            response_shift = DOT_PRODUCT(curv, response_coupled(:, r))
            CALL para_env_inter_kp%sum(response_shift)
            IF (curv_projection_sum > EPSILON(curv_projection_sum)) &
               response_coupled(:, r) = response_coupled(:, r) - response_shift/curv_projection_sum
            DO ie = 1, ne
               IF (curv(ie) <= EPSILON(curv(ie))) response_coupled(ie, r) = 0.0_dp
            END DO
            CALL project_kpoint_ot_schur_response(response_split(:, r), occ_curv, tether_curv, &
                                                  occ_curv_projection_sum, channel_nmo, e_offset)
         END DO
         projector_response_gap_avg = 0.0_dp
         projector_response_gap_max = 0.0_dp
         physical_response_gap_avg = 0.0_dp
         physical_response_gap_max = 0.0_dp
         DO r = 1, nrot
            projector_response_gap_norm = DOT_PRODUCT(occ_curv(:)*(response_occ(:, r) - response(:, r)), &
                                                      response_occ(:, r) - response(:, r))
            CALL para_env_inter_kp%sum(projector_response_gap_norm)
            IF (occ_curv_projection_sum > EPSILON(occ_curv_projection_sum)) &
               projector_response_gap_avg = projector_response_gap_avg + &
                  SQRT(MAX(0.0_dp, projector_response_gap_norm/occ_curv_projection_sum))
            projector_response_gap_max = MAX(projector_response_gap_max, &
                                             MAXVAL(ABS(response_occ(:, r) - response(:, r))))
            physical_response_gap_norm = DOT_PRODUCT(occ_curv(:)*(response_occ_physical(:, r) - response_occ(:, r)), &
                                                     response_occ_physical(:, r) - response_occ(:, r))
            CALL para_env_inter_kp%sum(physical_response_gap_norm)
            IF (occ_curv_projection_sum > EPSILON(occ_curv_projection_sum)) &
               physical_response_gap_avg = physical_response_gap_avg + &
                  SQRT(MAX(0.0_dp, physical_response_gap_norm/occ_curv_projection_sum))
            physical_response_gap_max = MAX(physical_response_gap_max, &
                                            MAXVAL(ABS(response_occ_physical(:, r) - response_occ(:, r))))
         END DO
         projector_response_gap_avg = projector_response_gap_avg/REAL(MAX(1, nrot), KIND=dp)
         physical_response_gap_avg = physical_response_gap_avg/REAL(MAX(1, nrot), KIND=dp)
         CALL para_env_inter_kp%max(projector_response_gap_max)
         CALL para_env_inter_kp%max(physical_response_gap_max)
         CALL build_kpoint_ot_exact_energy_response(occ_curv, tether_curv, response_occ, response, grad_e, &
                                                    channel_nmo, e_offset, energy_hessian, energy_cross, &
                                                    response_exact, z_exact_trial, ipiv_energy, &
                                                    exact_energy_model_ready, energy_gauge_weight)
         exact_physical_response_ready = .FALSE.
         IF (exact_energy_model_ready .AND. projector_physical_shadow_enabled) &
            CALL build_kpoint_ot_exact_energy_response(occ_curv, tether_curv, response_occ_physical, response, &
                                                       grad_e, channel_nmo, e_offset, energy_hessian, &
                                                       energy_cross_physical, response_exact_physical, &
                                                       z_exact_physical_trial, ipiv_energy, &
                                                       exact_physical_response_ready, energy_gauge_weight)
         IF (exact_energy_model_ready) &
            CALL build_kpoint_ot_exact_energy_response(occ_curv, tether_curv, response, response, grad_e, &
                                                       channel_nmo, e_offset, energy_hessian, energy_cross_raw, &
                                                       response_exact_raw, z_exact_raw_trial, ipiv_energy, &
                                                       exact_energy_model_ready, energy_gauge_weight)
         IF (exact_energy_model_ready) THEN
            z_shift = DOT_PRODUCT(curv, z_exact_trial)
            CALL para_env_inter_kp%sum(z_shift)
            IF (curv_projection_sum > EPSILON(curv_projection_sum)) &
               z_exact_trial(:) = z_exact_trial(:) - z_shift/curv_projection_sum
            DO r = 1, nrot
               exact_response_diag_avg = exact_response_diag_avg + &
                                         DOT_PRODUCT(energy_cross(:, r), response_exact(:, r))
               response_mixed_diag_avg = response_mixed_diag_avg + &
                                         0.5_dp*(DOT_PRODUCT(energy_cross(:, r), response_exact_raw(:, r)) + &
                                                 DOT_PRODUCT(energy_cross_raw(:, r), response_exact(:, r)))
            END DO
            exact_response_diag_avg = exact_response_diag_avg/REAL(MAX(1, nrot), KIND=dp)
            response_mixed_diag_avg = response_mixed_diag_avg/REAL(MAX(1, nrot), KIND=dp)
         END IF

         IF (sigma > EPSILON(sigma)) THEN
            response_trust_base = schur_response_trust_factor*kpoint_ot_fermi_window_anchor(sigma)
            response_abs_trust = ot_kp_fermi_window_outer_factor*sigma
            IF (residual_delta > 0.0_dp .AND. ot_kp_canonical_density_max_weight > 0.0_dp) THEN
               density_projector_fraction = &
                  kpoint_ot_canonical_projector_weight(residual_delta, ot_kp_canonical_density_max_weight)/ &
                  ot_kp_canonical_density_max_weight
               density_projector_fraction = MIN(1.0_dp, MAX(0.0_dp, density_projector_fraction))
               response_trust_factor = MAX(ot_kp_schur_response_trust_min_factor, &
                                           1.0_dp - density_projector_fraction)
            END IF
            response_trust = response_trust_base*response_trust_factor
         END IF
         solved = .FALSE.
         trust_rejected = .FALSE.
         schur_used_attempt = 0
         ! Evaluate the split fixed-N/tether Schur model under HIGH audit, or
         ! after accepted Mermin line-search drops validate the model quality.
         split_candidate_enabled = new_search_point .AND. local_ot_env(1)%schur_split_candidate_enabled
         constrained_candidate_enabled = new_search_point .AND. &
                                         local_ot_env(1)%schur_constrained_candidate_enabled
         split_candidate_diagnostic = logger%iter_info%print_level >= high_print_level .OR. &
                                      split_candidate_enabled .OR. constrained_candidate_enabled
         IF (local_ot_env(1)%schur_split_quality_samples >= ot_kp_schur_calibration_min_samples) THEN
            split_quality_for_calibration = &
               kpoint_ot_schur_history_quality(local_ot_env(1)%schur_split_drop_sum, &
                                               local_ot_env(1)%schur_split_measured_sum, &
                                               local_ot_env(1)%schur_split_quality_ema, &
                                               local_ot_env(1)%schur_split_quality_samples, &
                                               ot_kp_schur_calibration_min_samples, &
                                               local_ot_env(1)%schur_split_quality_ema, &
                                               ot_kp_schur_calibration_floor, 1.0_dp)
            split_quality_for_calibration = MAX(split_quality_for_calibration, &
                                                MAX(ot_kp_schur_calibration_floor, &
                                                    MIN(1.0_dp, local_ot_env(1)%schur_split_quality_ema)))
         ELSE IF (local_ot_env(1)%schur_model_quality_samples >= ot_kp_schur_calibration_min_samples) THEN
            ! Bootstrap Split curvature from the accepted dense/coupled Mermin
            ! model quality.  Otherwise Split candidates that are already known
            ! to be over-optimistic can be rejected by the OT z.z=g.z scale gate
            ! before they ever collect their own calibration samples.
            split_quality_for_calibration = &
               kpoint_ot_schur_history_quality(local_ot_env(1)%schur_coupled_drop_sum, &
                                               local_ot_env(1)%schur_coupled_measured_sum, &
                                               local_ot_env(1)%schur_model_quality_ema, &
                                               local_ot_env(1)%schur_model_quality_samples, &
                                               ot_kp_schur_calibration_min_samples, &
                                               local_ot_env(1)%schur_model_quality_ema, &
                                               ot_kp_schur_calibration_floor, 1.0_dp)
            split_quality_for_calibration = MAX(split_quality_for_calibration, &
                                                MAX(ot_kp_schur_calibration_floor, &
                                                    MIN(1.0_dp, local_ot_env(1)%schur_model_quality_ema)))
         END IF
         IF (local_ot_env(1)%schur_coupled_quality_samples >= ot_kp_schur_calibration_min_samples) THEN
            dense_quality_for_selection = &
               kpoint_ot_schur_history_quality(local_ot_env(1)%schur_coupled_drop_sum, &
                                               local_ot_env(1)%schur_coupled_measured_sum, &
                                               local_ot_env(1)%schur_coupled_quality_ema, &
                                               local_ot_env(1)%schur_coupled_quality_samples, &
                                               ot_kp_schur_calibration_min_samples, 1.0_dp, &
                                               ot_kp_schur_calibrated_quality_min_ratio, &
                                               ot_kp_schur_calibrated_quality_max_ratio)
            dense_quality_for_selection = MAX(dense_quality_for_selection, &
                                              MAX(ot_kp_schur_calibrated_quality_min_ratio, &
                                                  MIN(ot_kp_schur_calibrated_quality_max_ratio, &
                                                      local_ot_env(1)%schur_coupled_quality_ema)))
         END IF
         projector_mode_history_quality = dense_quality_for_selection
         IF (local_ot_env(1)%schur_projector_mode_quality_samples >= ot_kp_schur_calibration_min_samples) THEN
            projector_mode_history_quality = &
               kpoint_ot_schur_history_quality(local_ot_env(1)%schur_projector_mode_drop_sum, &
                                               local_ot_env(1)%schur_projector_mode_measured_sum, &
                                               local_ot_env(1)%schur_projector_mode_quality_ema, &
                                               local_ot_env(1)%schur_projector_mode_quality_samples, &
                                               ot_kp_schur_calibration_min_samples, &
                                               dense_quality_for_selection, &
                                               ot_kp_schur_calibrated_quality_min_ratio, &
                                               ot_kp_schur_calibrated_quality_max_ratio)
         END IF
         IF (local_ot_env(1)%schur_projector_mode_drop_sum > &
             EPSILON(local_ot_env(1)%schur_projector_mode_drop_sum)) &
            projector_mode_history_response_ratio = &
               local_ot_env(1)%schur_projector_mode_response_drop_sum/ &
               local_ot_env(1)%schur_projector_mode_drop_sum
         projector_mode_history_response_negative = &
            local_ot_env(1)%schur_projector_mode_response_negative_ema
            projector_response_quality_for_selection = &
               kpoint_ot_schur_history_quality(local_ot_env(1)%schur_projector_response_drop_sum, &
                                               local_ot_env(1)%schur_projector_response_measured_sum, &
                                               local_ot_env(1)%schur_projector_response_quality_ema, &
                                            local_ot_env(1)%schur_projector_response_quality_samples, &
                                            ot_kp_schur_projector_response_min_samples, &
                                               0.0_dp, &
                                               ot_kp_schur_projector_response_quality_min_ratio, &
                                               ot_kp_schur_projector_response_quality_max_ratio)
            projector_physical_quality_for_selection = &
               kpoint_ot_schur_history_quality(local_ot_env(1)%schur_projector_physical_drop_sum, &
                                               local_ot_env(1)%schur_projector_physical_measured_sum, &
                                               local_ot_env(1)%schur_projector_physical_quality_ema, &
                                               local_ot_env(1)%schur_projector_physical_quality_samples, &
                                               ot_kp_schur_projector_response_min_samples, &
                                               0.0_dp, &
                                               ot_kp_schur_projector_response_quality_min_ratio, &
                                               ot_kp_schur_projector_response_quality_max_ratio)
            split_quality_for_selection = &
            kpoint_ot_schur_history_quality(local_ot_env(1)%schur_split_calibrated_drop_sum, &
                                            local_ot_env(1)%schur_split_measured_sum, &
                                            local_ot_env(1)%schur_split_calibrated_quality_ema, &
                                            local_ot_env(1)%schur_split_quality_samples, &
                                            ot_kp_schur_calibration_min_samples, &
                                            local_ot_env(1)%schur_split_calibrated_quality_ema, &
                                            ot_kp_schur_calibrated_quality_min_ratio, &
                                            ot_kp_schur_calibrated_quality_max_ratio)
         constrained_quality_for_selection = &
            kpoint_ot_schur_history_quality(local_ot_env(1)%schur_constrained_calibrated_drop_sum, &
                                            local_ot_env(1)%schur_constrained_measured_sum, &
                                            local_ot_env(1)%schur_constrained_calibrated_quality_ema, &
                                            local_ot_env(1)%schur_constrained_quality_samples, &
                                            ot_kp_schur_constrained_min_samples, &
                                            local_ot_env(1)%schur_constrained_calibrated_quality_ema, &
                                            ot_kp_schur_calibrated_quality_min_ratio, &
                                            ot_kp_schur_calibrated_quality_max_ratio)
         IF (local_ot_env(1)%schur_split_quality_samples >= ot_kp_schur_calibration_min_samples .OR. &
             local_ot_env(1)%schur_coupled_quality_samples >= ot_kp_schur_calibration_min_samples) &
            residual_reference_quality = ot_kp_schur_quality_min_ratio
         IF (local_ot_env(1)%schur_split_quality_samples >= ot_kp_schur_calibration_min_samples) &
            residual_reference_quality = MAX(residual_reference_quality, &
                                             local_ot_env(1)%schur_split_residual_quality_ema)
         IF (local_ot_env(1)%schur_coupled_quality_samples >= ot_kp_schur_calibration_min_samples) &
            residual_reference_quality = MAX(residual_reference_quality, &
                                             local_ot_env(1)%schur_coupled_residual_quality_ema)
         IF (residual_reference_quality > EPSILON(residual_reference_quality)) THEN
            IF (local_ot_env(1)%schur_split_quality_samples >= ot_kp_schur_calibration_min_samples) &
               split_residual_quality = &
                  MAX(ot_kp_schur_calibrated_quality_min_ratio, &
                      MIN(1.0_dp, local_ot_env(1)%schur_split_residual_quality_ema/residual_reference_quality))
            IF (local_ot_env(1)%schur_coupled_quality_samples >= ot_kp_schur_calibration_min_samples) &
               dense_residual_quality = &
                  MAX(ot_kp_schur_calibrated_quality_min_ratio, &
                      MIN(1.0_dp, local_ot_env(1)%schur_coupled_residual_quality_ema/residual_reference_quality))
            split_quality_for_selection = split_quality_for_selection*split_residual_quality
            dense_quality_for_selection = dense_quality_for_selection*dense_residual_quality
         END IF
         IF (local_ot_env(1)%schur_split_quality_samples >= ot_kp_schur_calibration_min_samples .OR. &
             local_ot_env(1)%schur_model_quality_samples >= ot_kp_schur_calibration_min_samples) THEN
            split_curvature_scale = MIN(ot_kp_schur_calibration_max_scale, &
                                        1.0_dp/SQRT(split_quality_for_calibration))
         END IF
         ! Attempt order: response-constrained split shadow, raw split,
         ! dense fixed-N coupled response, and the conservative raw fallback.
         DO schur_attempt = 1, 4
            attempt_failed = .FALSE.
            reject_reason = 0
            direction_scale = 0.0_dp
            scale_quality = 0.0_dp
            gdotp = 0.0_dp
            gdot_energy = 0.0_dp
            gdot_rotation = 0.0_dp
            gnorm = 0.0_dp
            norm_energy = 0.0_dp
            norm_rotation = 0.0_dp
            norm_ratio = 0.0_dp
            normp = 0.0_dp
            exact_attempt_response = exact_energy_model_ready .AND. schur_attempt <= 3
            response_model_available = schur_attempt <= 2 .OR. &
                                       (schur_attempt == 3 .AND. exact_attempt_response)
            projector_coupled_weight = 0.0_dp
            projector_diag_compensation = 0.0_dp
            projector_mode_damped_abs_avg = 0.0_dp
            projector_mode_eval_max = 0.0_dp
            projector_mode_eval_min = 0.0_dp
            projector_mode_negative_fraction = 0.0_dp
            projector_mode_quality = 0.0_dp
            projector_mode_response_correction_drop = 0.0_dp
            projector_mode_response_correction_negative = 0.0_dp
            projector_mode_response_negative_coeff = 0.0_dp
            projector_mode_response_positive_coeff = 0.0_dp
            projector_mode_response_cross_drop = 0.0_dp
            projector_mode_response_negative_drop = 0.0_dp
            projector_mode_response_positive_drop = 0.0_dp
            projector_mode_response_drop = 0.0_dp
            projector_mode_response_drop_ratio = 0.0_dp
            projector_mode_response_negative_fraction = 0.0_dp
            projector_mode_response_raw_negative_fraction = 0.0_dp
            projector_mode_step_negative_fraction = 0.0_dp
               projector_mode_trust_quality = 1.0_dp
               projector_mode_response_weight = 0.0_dp
               projector_physical_candidate_ready = .FALSE.
               projector_physical_candidate_selected = .FALSE.
               projector_physical_direct_candidate_ready = .FALSE.
               projector_physical_direct_candidate_selected = .FALSE.
               projector_physical_direct_candidate_used = .FALSE.
               projector_physical_direct_diag_compensation = 0.0_dp
               projector_physical_direct_diag_shift_avg = 0.0_dp
               projector_physical_direct_diag_shift_sum = 0.0_dp
               projector_physical_direct_direction_scale = 0.0_dp
               projector_physical_direct_expected_drop = 0.0_dp
               projector_physical_direct_model_calibrated_drop = 0.0_dp
               projector_physical_direct_model_drop = 0.0_dp
               projector_physical_direct_model_gap = 0.0_dp
               projector_physical_direct_model_gap_drop = 0.0_dp
               projector_physical_direct_model_gap_ratio = 0.0_dp
               projector_physical_direct_model_gap_trust = HUGE(1.0_dp)
               projector_physical_direct_response_abs_step = 0.0_dp
               projector_physical_direct_response_abs_trust_ratio = 0.0_dp
               projector_physical_direct_response_weighted_step = 0.0_dp
               projector_physical_direct_response_weighted_trust_ratio = 0.0_dp
               projector_physical_expected_drop = 0.0_dp
               projector_physical_model_calibrated_drop = 0.0_dp
               projector_physical_model_drop = 0.0_dp
               projector_physical_model_gap = 0.0_dp
               projector_physical_model_gap_drop = 0.0_dp
               projector_physical_model_gap_ratio = 0.0_dp
               projector_physical_model_gap_trust = HUGE(1.0_dp)
               projector_response_candidate_ready = .FALSE.
            projector_response_candidate_selected = .FALSE.
            projector_response_expected_drop = 0.0_dp
            projector_response_model_calibrated_drop = 0.0_dp
            projector_response_model_drop = 0.0_dp
            projector_response_model_gap_ratio = 0.0_dp
            IF (exact_attempt_response) THEN
               z_e(:) = z_exact_trial(:)
            ELSE
               z_e(:) = z_trial(:)
            END IF
            IF (schur_attempt == 1 .AND. logger%iter_info%print_level < high_print_level .AND. &
                .NOT. split_candidate_enabled .AND. .NOT. constrained_candidate_enabled) THEN
               constrained_reject_reason = reject_quality_gate
               CYCLE
            END IF
            IF (schur_attempt <= 2 .AND. .NOT. split_candidate_diagnostic) THEN
               IF (schur_attempt == 1) THEN
                  constrained_reject_reason = reject_quality_gate
               ELSE
                  coupled_reject_reason = reject_quality_gate
               END IF
               CYCLE
            END IF
            IF (schur_attempt <= 2) THEN
               schur(:, :) = rot_hessian_coupled(:, :) + rot_hessian_projector(:, :) - rot_hessian(:, :)
            ELSE IF (schur_attempt == 3) THEN
               IF (exact_attempt_response) &
                  projector_coupled_weight = kpoint_ot_projector_coupled_weight(residual_delta, &
                     dense_quality_for_selection, local_ot_env(1)%schur_coupled_quality_samples, &
                     density_projector_fraction, projector_hessian_diag_avg, projector_hessian_delta_diag_avg)
               schur(:, :) = rot_hessian_coupled(:, :) + &
                  projector_coupled_weight*(rot_hessian_projector(:, :) - rot_hessian(:, :))
               IF (residual_delta <= ot_kp_schur_projector_coupled_tail_delta) &
                  projector_diag_compensation = &
                     MAX(0.0_dp, -projector_coupled_weight*projector_hessian_delta_diag_avg)
            ELSE
               schur(:, :) = rot_hessian(:, :)
            END IF
            DO r = 1, nrot
               IF (schur_attempt <= 2) THEN
                  IF (exact_energy_model_ready) THEN
                     rhs(r) = grad_rot(r) + DOT_PRODUCT(response_exact(:, r), grad_e)
                  ELSE
                     rhs(r) = grad_rot(r) + DOT_PRODUCT(response_split(:, r), grad_e)
                  END IF
               ELSE IF (schur_attempt == 3) THEN
                  IF (exact_attempt_response) THEN
                     rhs(r) = grad_rot(r) + DOT_PRODUCT(response_exact(:, r), grad_e)
                  ELSE
                     rtd(r) = DOT_PRODUCT(curv, response_coupled(:, r))
                     rhs(r) = grad_rot(r) + DOT_PRODUCT(response_coupled(:, r), grad_e)
                  END IF
               ELSE
                  rtd(r) = DOT_PRODUCT(curv, response(:, r))
                  rhs(r) = grad_rot(r) + DOT_PRODUCT(response(:, r), grad_e)
               END IF
            END DO
            diag_shift_sum = 0.0_dp
            diag_shift_max = 0.0_dp
            response_abs_step = 0.0_dp
            response_abs_trust_ratio = 0.0_dp
            response_step_max = 0.0_dp
            response_weighted_step = 0.0_dp
            response_weighted_trust_ratio = 0.0_dp
            model_drop = 0.0_dp
            model_calibrated_drop = 0.0_dp
            model_gap = 0.0_dp
            model_gap_drop = 0.0_dp
            model_gap_ratio = 0.0_dp
            model_gap_trust = HUGE(1.0_dp)
            occgap_shadow_abs_step = 0.0_dp
            occgap_shadow_drop = 0.0_dp
            occgap_shadow_gap_ratio = 0.0_dp
            occgap_shadow_quality = 1.0_dp
            occgap_shadow_step = 0.0_dp
            occ_response_diag_sum = 0.0_dp
            response_penalty_diag_sum = 0.0_dp
            rot_diag_sum = 0.0_dp
            tether_response_diag_sum = 0.0_dp
            rot_floor_min = curvature_floor
            schur_diag_min = HUGE(1.0_dp)

            DO s = 1, nrot
               DO r = 1, nrot
                  response_penalty_term = 0.0_dp
                  IF (schur_attempt <= 2) THEN
                     CALL kpoint_ot_mixed_projected_cross_term(response_occ(:, r), response_occ(:, s), &
                                                               response(:, r), response(:, s), &
                                                               occ_curv, tether_curv, occ_curv_projection_sum, &
                                                               channel_nmo, e_offset, projected_term, &
                                                               occ_projected_term, tether_projected_term)
                     IF (exact_energy_model_ready) &
                        projected_term = 0.5_dp*(DOT_PRODUCT(energy_cross(:, r), response_exact_raw(:, s)) + &
                                                DOT_PRODUCT(energy_cross_raw(:, r), response_exact(:, s)))
                     IF (schur_attempt == 1) response_penalty_term = &
                        schur_response_penalty_scale*projected_term
                     IF (r == s) THEN
                        rot_diag_sum = rot_diag_sum + schur(r, r)
                        occ_response_diag_sum = occ_response_diag_sum + occ_projected_term
                        tether_response_diag_sum = tether_response_diag_sum + tether_projected_term
                        response_penalty_diag_sum = response_penalty_diag_sum + response_penalty_term
                     END IF
                  ELSE IF (schur_attempt == 3) THEN
                     IF (exact_attempt_response) THEN
                        projected_term = 0.5_dp*(DOT_PRODUCT(energy_cross(:, r), response_exact_raw(:, s)) + &
                                                DOT_PRODUCT(energy_cross_raw(:, r), response_exact(:, s)))
                     ELSE
                        dense_term = DOT_PRODUCT(curv(:)*response_coupled(:, r), response_coupled(:, s))
                        projected_term = dense_term
                        IF (curv_projection_sum > EPSILON(curv_projection_sum)) &
                           projected_term = projected_term - rtd(r)*rtd(s)/curv_projection_sum
                     END IF
                  ELSE
                     dense_term = DOT_PRODUCT(curv(:)*response(:, r), response(:, s))
                     projected_term = dense_term
                     IF (curv_projection_sum > EPSILON(curv_projection_sum)) &
                        projected_term = projected_term - rtd(r)*rtd(s)/curv_projection_sum
                  END IF
                  schur(r, s) = schur(r, s) - projected_term + response_penalty_term
               END DO
            END DO
            IF (schur_attempt <= 2) schur(:, :) = split_curvature_scale*schur(:, :)
            DO r = 1, nrot
               IF (projector_diag_compensation > 0.0_dp) schur(r, r) = schur(r, r) + projector_diag_compensation
               row_floor = curvature_floor
               IF (schur_attempt <= 3) THEN
                  row_floor = rot_floor(r)
                  IF (schur_attempt == 3 .AND. exact_attempt_response .AND. &
                      exact_physical_response_ready .AND. &
                      residual_delta <= ot_kp_schur_projector_coupled_tail_delta .AND. &
                      rot_floor(r) <= schur_frontier_curvature_floor + EPSILON(schur_frontier_curvature_floor)) &
                     row_floor = kpoint_ot_physical_rotation_floor(sigma, schur(r, r), &
                                                                   rot_hessian_projector_physical(r, r), &
                                                                   rot_floor(r))
               END IF
               rot_floor_min = MIN(rot_floor_min, row_floor)
               IF (schur_attempt <= 3) THEN
                  ! First try the coupled finite-rotation Schur block itself.  The
                  ! line-search safeguards below reject non-descent or overlarge
                  ! response directions; the Gershgorin model remains the fallback.
                  diag_shift = MAX(0.0_dp, row_floor - schur(r, r))
               ELSE
                  row_abs = SUM(ABS(schur(r, :))) - ABS(schur(r, r))
                  diag_shift = MAX(0.0_dp, row_abs + row_floor - schur(r, r))
               END IF
               IF (diag_shift > 0.0_dp) THEN
                  schur(r, r) = schur(r, r) + diag_shift
                  diag_shift_sum = diag_shift_sum + diag_shift
                  diag_shift_max = MAX(diag_shift_max, diag_shift)
               END IF
               schur_diag_min = MIN(schur_diag_min, schur(r, r))
            END DO

            sol(:) = rhs(:)
            CALL dgesv(nrot, 1, schur, nrot, ipiv, sol, nrot, info)
            IF (info /= 0) THEN
               attempt_failed = .TRUE.
               reject_reason = reject_dgesv
            END IF

            IF (.NOT. attempt_failed) THEN
               IF (schur_attempt <= 2) THEN
                  IF (exact_energy_model_ready) THEN
                     rx_e(:) = MATMUL(response_exact_raw, sol)
                     rx_model(:) = MATMUL(response_exact, sol)
                  ELSE
                     rx_e(:) = MATMUL(response_split, sol)
                     rx_model(:) = MATMUL(response_coupled, sol)
                  END IF
               ELSE IF (schur_attempt == 3) THEN
                  IF (exact_attempt_response) THEN
                     rx_e(:) = MATMUL(response_exact_raw, sol)
                     rx_model(:) = MATMUL(response_exact, sol)
                  ELSE
                     rx_e(:) = MATMUL(response_coupled, sol)
                     rx_model(:) = rx_e(:)
                  END IF
               ELSE
                  rx_e(:) = MATMUL(response, sol)
                  rx_model(:) = rx_e(:)
               END IF
               rx_shift = 0.0_dp
               rx_weight = 0.0_dp
               DO ie = 1, ne
                  rx_shift = rx_shift + curv(ie)*rx_e(ie)
                  rx_weight = rx_weight + curv(ie)
               END DO
               CALL para_env_inter_kp%sum(rx_shift)
               CALL para_env_inter_kp%sum(rx_weight)
               IF (rx_weight > EPSILON(rx_weight)) rx_e(:) = rx_e(:) - rx_shift/rx_weight
               IF (response_model_available) THEN
                  rx_shift = 0.0_dp
                  rx_weight = 0.0_dp
                  DO ie = 1, ne
                     rx_shift = rx_shift + curv(ie)*rx_model(ie)
                     rx_weight = rx_weight + curv(ie)
                  END DO
                  CALL para_env_inter_kp%sum(rx_shift)
                  CALL para_env_inter_kp%sum(rx_weight)
                  IF (rx_weight > EPSILON(rx_weight)) rx_model(:) = rx_model(:) - rx_shift/rx_weight
               END IF
               gnorm = DOT_PRODUCT(grad_e, grad_e) + DOT_PRODUCT(grad_rot, grad_rot)
               IF (schur_attempt == 3 .AND. exact_attempt_response .AND. &
                   projector_coupled_weight > EPSILON(projector_coupled_weight) .AND. &
                   local_ot_env(1)%schur_projector_mode_quality_samples >= &
                   ot_kp_schur_calibration_min_samples) THEN
                  CALL audit_kpoint_ot_projector_response_modes(rot_hessian, rot_hessian_projector, sol, &
                                                                response_exact_raw, response_exact, curv, &
                                                                projector_coupled_weight, residual_delta, &
                                                                projector_mode_history_quality, &
                                                                local_ot_env(1)%schur_projector_mode_trust_ema, &
                                                                  projector_mode_history_response_ratio, grad_e, &
                                                                  z_e + rx_e, rx_e, &
                                                                projector_response_rx_e, projector_mode_response_weight, &
                                                                projector_mode_response_correction_drop, &
                                                                projector_mode_response_correction_negative, &
                                                                projector_mode_response_raw_negative_fraction, &
                                                                projector_mode_response_positive_drop, &
                                                                projector_mode_response_negative_drop, &
                                                                projector_mode_response_cross_drop, &
                                                                projector_mode_response_positive_coeff, &
                                                                projector_mode_response_negative_coeff)
                  IF (projector_mode_response_weight > EPSILON(projector_mode_response_weight)) THEN
                     projector_response_z_e(:) = z_e(:) + projector_response_rx_e(:)
                     projector_response_gdot_energy = DOT_PRODUCT(grad_e, projector_response_z_e)
                     projector_response_gdot_rotation = DOT_PRODUCT(grad_rot, sol)
                     projector_response_norm_energy = DOT_PRODUCT(projector_response_z_e, projector_response_z_e)
                     projector_response_norm_rotation = DOT_PRODUCT(sol, sol)
                     projector_response_gdotp = projector_response_gdot_energy + projector_response_gdot_rotation
                     projector_response_normp = projector_response_norm_energy + projector_response_norm_rotation
                     projector_response_norm_ratio = projector_response_normp/MAX(gnorm, EPSILON(1.0_dp))
                     IF (projector_response_gdotp > 0.0_dp .AND. &
                         projector_response_normp > EPSILON(projector_response_normp) .AND. &
                         projector_response_norm_ratio <= max_norm_ratio) THEN
                        projector_response_direction_scale = &
                           projector_response_gdotp/MAX(projector_response_normp, EPSILON(projector_response_normp))
                        projector_response_scale_quality = &
                           kpoint_ot_schur_scale_quality(projector_response_direction_scale, schur_min_coupled_scale)
                        projector_response_z_e(:) = projector_response_direction_scale*projector_response_z_e(:)
                        projector_response_sol(:) = projector_response_direction_scale*sol(:)
                        projector_response_model_drop = 0.5_dp*projector_response_direction_scale* &
                                                        projector_response_gdotp
                        projector_response_response_abs_step = &
                           projector_response_direction_scale*MAXVAL(ABS(projector_response_rx_e))
                        CALL para_env_inter_kp%max(projector_response_response_abs_step)
                        projector_response_response_weighted_step = projector_response_direction_scale* &
                           kpoint_ot_weighted_response_step(projector_response_rx_e, curv, curv_projection_sum)
                        projector_response_model_gap = &
                           projector_response_direction_scale*MAXVAL(ABS(rx_model - projector_response_rx_e))
                        projector_response_model_gap_drop = &
                           kpoint_ot_response_mismatch_drop(projector_response_rx_e, rx_model, curv, &
                                                            projector_response_direction_scale)
                        projector_response_model_calibrated_drop = &
                           MAX(0.0_dp, projector_response_model_drop - projector_response_model_gap_drop)
                        projector_response_model_gap_trust = schur_model_gap_factor*response_trust
                        projector_response_model_gap_ratio = 0.0_dp
                        IF (projector_response_model_gap_trust > EPSILON(projector_response_model_gap_trust)) &
                           projector_response_model_gap_ratio = projector_response_model_gap/ &
                                                                projector_response_model_gap_trust
                        projector_response_response_weighted_trust_ratio = 0.0_dp
                        IF (response_trust > EPSILON(response_trust)) &
                           projector_response_response_weighted_trust_ratio = &
                              projector_response_response_weighted_step/response_trust
                        projector_response_response_abs_trust_ratio = 0.0_dp
                        IF (response_abs_trust > EPSILON(response_abs_trust)) &
                           projector_response_response_abs_trust_ratio = &
                              projector_response_response_abs_step/response_abs_trust
                        projector_response_trust_quality = &
                           kpoint_ot_schur_trust_quality(projector_response_response_weighted_trust_ratio, &
                                                         projector_response_response_abs_trust_ratio, &
                                                         projector_response_model_gap_ratio)
                        IF (projector_response_direction_scale < schur_min_coupled_scale) &
                           projector_response_model_calibrated_drop = &
                              projector_response_scale_quality*projector_response_model_calibrated_drop
                        projector_response_expected_drop = &
                           kpoint_ot_schur_expected_drop(projector_response_model_calibrated_drop, &
                                                         projector_response_quality_for_selection, &
                                                         projector_response_trust_quality)
                        projector_response_candidate_ready = &
                           projector_response_model_calibrated_drop > EPSILON(projector_response_model_drop) .AND. &
                           projector_response_response_weighted_step <= response_trust .AND. &
                           projector_response_response_abs_step <= response_abs_trust
                     END IF
                     END IF
                  END IF
                  IF (schur_attempt == 3 .AND. exact_attempt_response .AND. exact_physical_response_ready .AND. &
                      projector_physical_shadow_enabled .AND. &
                      projector_coupled_weight > EPSILON(projector_coupled_weight) .AND. &
                      local_ot_env(1)%schur_projector_mode_quality_samples >= &
                      ot_kp_schur_calibration_min_samples) THEN
                     projector_physical_model_e(:) = MATMUL(response_exact_physical, sol)
                     rx_shift = DOT_PRODUCT(curv, projector_physical_model_e)
                     CALL para_env_inter_kp%sum(rx_shift)
                     IF (rx_weight > EPSILON(rx_weight)) &
                        projector_physical_model_e(:) = projector_physical_model_e(:) - rx_shift/rx_weight
                     CALL audit_kpoint_ot_projector_response_modes(rot_hessian_projector, &
                                                                   rot_hessian_projector_physical, sol, &
                                                                   response_exact_raw, response_exact_physical, curv, &
                                                                   projector_coupled_weight, residual_delta, &
                                                                   projector_mode_history_quality, &
                                                                   local_ot_env(1)%schur_projector_mode_trust_ema, &
                                                                  projector_mode_history_response_ratio, grad_e, &
                                                                  z_e + rx_e, rx_e, &
                                                                   projector_physical_rx_e, projector_physical_response_weight, &
                                                                   projector_physical_response_correction_drop, &
                                                                   projector_physical_response_correction_negative, &
                                                                   projector_physical_response_raw_negative_fraction, &
                                                                   projector_physical_response_positive_drop, &
                                                                   projector_physical_response_negative_drop, &
                                                                   projector_physical_response_cross_drop, &
                                                                   projector_physical_response_positive_coeff, &
                                                                   projector_physical_response_negative_coeff)
                     IF (projector_physical_response_weight > EPSILON(projector_physical_response_weight)) THEN
                        projector_physical_z_e(:) = z_e(:) + projector_physical_rx_e(:)
                        projector_physical_gdot_energy = DOT_PRODUCT(grad_e, projector_physical_z_e)
                        projector_physical_gdot_rotation = DOT_PRODUCT(grad_rot, sol)
                        projector_physical_norm_energy = DOT_PRODUCT(projector_physical_z_e, projector_physical_z_e)
                        projector_physical_norm_rotation = DOT_PRODUCT(sol, sol)
                        projector_physical_gdotp = projector_physical_gdot_energy + projector_physical_gdot_rotation
                        projector_physical_normp = projector_physical_norm_energy + projector_physical_norm_rotation
                        projector_physical_norm_ratio = projector_physical_normp/MAX(gnorm, EPSILON(1.0_dp))
                        IF (projector_physical_gdotp > 0.0_dp .AND. &
                            projector_physical_normp > EPSILON(projector_physical_normp) .AND. &
                            projector_physical_norm_ratio <= max_norm_ratio) THEN
                           projector_physical_direction_scale = &
                              projector_physical_gdotp/MAX(projector_physical_normp, EPSILON(projector_physical_normp))
                           projector_physical_scale_quality = &
                              kpoint_ot_schur_scale_quality(projector_physical_direction_scale, schur_min_coupled_scale)
                           projector_physical_z_e(:) = projector_physical_direction_scale*projector_physical_z_e(:)
                           projector_physical_sol(:) = projector_physical_direction_scale*sol(:)
                           projector_physical_model_drop = 0.5_dp*projector_physical_direction_scale* &
                                                           projector_physical_gdotp
                           projector_physical_response_abs_step = &
                              projector_physical_direction_scale*MAXVAL(ABS(projector_physical_rx_e))
                           CALL para_env_inter_kp%max(projector_physical_response_abs_step)
                           projector_physical_response_weighted_step = projector_physical_direction_scale* &
                              kpoint_ot_weighted_response_step(projector_physical_rx_e, curv, curv_projection_sum)
                           projector_physical_model_gap = &
                              projector_physical_direction_scale* &
                              MAXVAL(ABS(projector_physical_model_e - projector_physical_rx_e))
                           projector_physical_model_gap_drop = &
                              kpoint_ot_response_mismatch_drop(projector_physical_rx_e, projector_physical_model_e, curv, &
                                                               projector_physical_direction_scale)
                           projector_physical_model_calibrated_drop = &
                              MAX(0.0_dp, projector_physical_model_drop - projector_physical_model_gap_drop)
                           projector_physical_model_gap_trust = schur_model_gap_factor*response_trust
                           projector_physical_model_gap_ratio = 0.0_dp
                           IF (projector_physical_model_gap_trust > EPSILON(projector_physical_model_gap_trust)) &
                              projector_physical_model_gap_ratio = projector_physical_model_gap/ &
                                                                   projector_physical_model_gap_trust
                           projector_physical_response_weighted_trust_ratio = 0.0_dp
                           IF (response_trust > EPSILON(response_trust)) &
                              projector_physical_response_weighted_trust_ratio = &
                                 projector_physical_response_weighted_step/response_trust
                           projector_physical_response_abs_trust_ratio = 0.0_dp
                           IF (response_abs_trust > EPSILON(response_abs_trust)) &
                              projector_physical_response_abs_trust_ratio = &
                                 projector_physical_response_abs_step/response_abs_trust
                           projector_physical_trust_quality = &
                              kpoint_ot_schur_trust_quality(projector_physical_response_weighted_trust_ratio, &
                                                            projector_physical_response_abs_trust_ratio, &
                                                            projector_physical_model_gap_ratio)
                           IF (projector_physical_direction_scale < schur_min_coupled_scale) &
                              projector_physical_model_calibrated_drop = &
                                 projector_physical_scale_quality*projector_physical_model_calibrated_drop
                           projector_physical_expected_drop = &
                              kpoint_ot_schur_expected_drop(projector_physical_model_calibrated_drop, &
                                                            projector_physical_quality_for_selection, &
                                                            projector_physical_trust_quality)
                                             projector_physical_candidate_ready = &
                                                projector_physical_model_calibrated_drop > &
                                                EPSILON(projector_physical_model_drop) .AND. &
                                                projector_physical_response_weighted_step <= response_trust .AND. &
                                                projector_physical_response_abs_step <= response_abs_trust
                                          END IF
                                       END IF
                                    END IF
                                    IF (schur_attempt == 3 .AND. exact_attempt_response .AND. &
                                        exact_physical_response_ready .AND. &
                                        projector_physical_shadow_enabled .AND. &
                                        projector_coupled_weight > EPSILON(projector_coupled_weight) .AND. &
                                        local_ot_env(1)%schur_projector_mode_quality_samples >= &
                                        ot_kp_schur_calibration_min_samples) THEN
                                       schur(:, :) = rot_hessian_coupled(:, :) + &
                                                     projector_coupled_weight*(rot_hessian_projector_physical(:, :) - &
                                                                               rot_hessian(:, :))
                                       DO r = 1, nrot
                                          projector_physical_direct_sol(r) = grad_rot(r) + &
                                                                             DOT_PRODUCT(response_exact_physical(:, r), grad_e)
                                       END DO
                                       DO s = 1, nrot
                                          DO r = 1, nrot
                                             projected_term = 0.5_dp* &
                                                (DOT_PRODUCT(energy_cross_physical(:, r), response_exact_raw(:, s)) + &
                                                 DOT_PRODUCT(energy_cross_raw(:, r), response_exact_physical(:, s)))
                                             schur(r, s) = schur(r, s) - projected_term
                                          END DO
                                       END DO
                                       projector_physical_direct_diag_compensation = 0.0_dp
                                       IF (residual_delta <= ot_kp_schur_projector_coupled_tail_delta) &
                                          projector_physical_direct_diag_compensation = &
                                             MAX(0.0_dp, -projector_coupled_weight*projector_physical_delta_diag_avg)
                                       projector_physical_direct_diag_shift_sum = 0.0_dp
                                       DO r = 1, nrot
                                          IF (projector_physical_direct_diag_compensation > 0.0_dp) &
                                             schur(r, r) = schur(r, r) + projector_physical_direct_diag_compensation
                                          ! The exact fixed-N Mermin Schur complement has genuinely soft
                                          ! frontier rotation modes.  Do not re-harden the physical direct
                                          ! candidate with the generic frontier floor once its own diagonal
                                          ! curvature is already positive.
                                          row_floor = kpoint_ot_physical_rotation_floor(sigma, schur(r, r), &
                                                                                       rot_hessian_projector_physical(r, r), &
                                                                                       rot_floor(r))
                                          diag_shift = MAX(0.0_dp, row_floor - schur(r, r))
                                          IF (diag_shift > 0.0_dp) THEN
                                             schur(r, r) = schur(r, r) + diag_shift
                                             projector_physical_direct_diag_shift_sum = &
                                                projector_physical_direct_diag_shift_sum + diag_shift
                                          END IF
                                       END DO
                                       projector_physical_direct_diag_shift_avg = &
                                          projector_physical_direct_diag_shift_sum/REAL(MAX(1, nrot), KIND=dp)
                                       CALL dgesv(nrot, 1, schur, nrot, ipiv, projector_physical_direct_sol, nrot, info)
                                       IF (info == 0) THEN
                                          projector_physical_direct_rx_e(:) = &
                                             MATMUL(response_exact_raw, projector_physical_direct_sol)
                                          projector_physical_direct_model_e(:) = MATMUL(response_exact_physical, &
                                                                                         projector_physical_direct_sol)
                                          rx_shift = DOT_PRODUCT(curv, projector_physical_direct_rx_e)
                                          CALL para_env_inter_kp%sum(rx_shift)
                                          IF (curv_projection_sum > EPSILON(curv_projection_sum)) &
                                             projector_physical_direct_rx_e(:) = projector_physical_direct_rx_e(:) - &
                                                                                 rx_shift/curv_projection_sum
                                          rx_shift = DOT_PRODUCT(curv, projector_physical_direct_model_e)
                                          CALL para_env_inter_kp%sum(rx_shift)
                                          IF (curv_projection_sum > EPSILON(curv_projection_sum)) &
                                             projector_physical_direct_model_e(:) = projector_physical_direct_model_e(:) - &
                                                                                    rx_shift/curv_projection_sum
                                          projector_physical_direct_z_e(:) = z_e(:) + projector_physical_direct_rx_e(:)
                                          projector_physical_direct_gdot_energy = &
                                             DOT_PRODUCT(grad_e, projector_physical_direct_z_e)
                                          projector_physical_direct_gdot_rotation = &
                                             DOT_PRODUCT(grad_rot, projector_physical_direct_sol)
                                          projector_physical_direct_norm_energy = DOT_PRODUCT(projector_physical_direct_z_e, &
                                                                                              projector_physical_direct_z_e)
                                          projector_physical_direct_norm_rotation = DOT_PRODUCT(projector_physical_direct_sol, &
                                                                                                projector_physical_direct_sol)
                                          projector_physical_direct_gdotp = projector_physical_direct_gdot_energy + &
                                                                            projector_physical_direct_gdot_rotation
                                          projector_physical_direct_normp = projector_physical_direct_norm_energy + &
                                                                            projector_physical_direct_norm_rotation
                                          projector_physical_direct_norm_ratio = projector_physical_direct_normp/ &
                                                                                 MAX(gnorm, EPSILON(1.0_dp))
                                          IF (projector_physical_direct_gdotp > 0.0_dp .AND. &
                                              projector_physical_direct_normp > EPSILON(projector_physical_direct_normp) .AND. &
                                              projector_physical_direct_norm_ratio <= max_norm_ratio) THEN
                                             projector_physical_direct_direction_scale = &
                                                projector_physical_direct_gdotp/ &
                                                MAX(projector_physical_direct_normp, EPSILON(projector_physical_direct_normp))
                                             projector_physical_direct_scale_quality = &
                                                kpoint_ot_schur_scale_quality(projector_physical_direct_direction_scale, &
                                                                              schur_min_coupled_scale)
                                             projector_physical_direct_z_e(:) = projector_physical_direct_direction_scale* &
                                                                                 projector_physical_direct_z_e(:)
                                             projector_physical_direct_sol(:) = projector_physical_direct_direction_scale* &
                                                                                 projector_physical_direct_sol(:)
                                             projector_physical_direct_model_drop = 0.5_dp* &
                                                projector_physical_direct_direction_scale*projector_physical_direct_gdotp
                                             projector_physical_direct_response_abs_step = &
                                                projector_physical_direct_direction_scale* &
                                                MAXVAL(ABS(projector_physical_direct_rx_e))
                                             CALL para_env_inter_kp%max(projector_physical_direct_response_abs_step)
                                             projector_physical_direct_response_weighted_step = &
                                                projector_physical_direct_direction_scale* &
                                                kpoint_ot_weighted_response_step(projector_physical_direct_rx_e, curv, &
                                                                                 curv_projection_sum)
                                             projector_physical_direct_model_gap = &
                                                projector_physical_direct_direction_scale* &
                                                MAXVAL(ABS(projector_physical_direct_model_e - projector_physical_direct_rx_e))
                                             projector_physical_direct_model_gap_drop = &
                                                kpoint_ot_response_mismatch_drop(projector_physical_direct_rx_e, &
                                                                                 projector_physical_direct_model_e, curv, &
                                                                                 projector_physical_direct_direction_scale)
                                             projector_physical_direct_model_calibrated_drop = &
                                                MAX(0.0_dp, projector_physical_direct_model_drop - &
                                                            projector_physical_direct_model_gap_drop)
                                             projector_physical_direct_model_gap_trust = schur_model_gap_factor*response_trust
                                             IF (projector_physical_direct_model_gap_trust > &
                                                 EPSILON(projector_physical_direct_model_gap_trust)) &
                                                projector_physical_direct_model_gap_ratio = &
                                                   projector_physical_direct_model_gap/ &
                                                   projector_physical_direct_model_gap_trust
                                             IF (response_trust > EPSILON(response_trust)) &
                                                projector_physical_direct_response_weighted_trust_ratio = &
                                                   projector_physical_direct_response_weighted_step/response_trust
                                             IF (response_abs_trust > EPSILON(response_abs_trust)) &
                                                projector_physical_direct_response_abs_trust_ratio = &
                                                   projector_physical_direct_response_abs_step/response_abs_trust
                                             projector_physical_direct_trust_quality = &
                                                kpoint_ot_schur_trust_quality( &
                                                   projector_physical_direct_response_weighted_trust_ratio, &
                                                   projector_physical_direct_response_abs_trust_ratio, &
                                                   projector_physical_direct_model_gap_ratio)
                                             ! The physical-direct candidate is selected only after
                                             ! accepted-drop history proves it conservative.  Keep the
                                             ! z.z=g.z scale diagnostic, but do not double-discount
                                             ! its model drop here.
                                             projector_physical_direct_expected_drop = &
                                                kpoint_ot_schur_expected_drop(projector_physical_direct_model_calibrated_drop, &
                                                                              projector_physical_quality_for_selection, &
                                                                              projector_physical_direct_trust_quality)
                                             projector_physical_direct_candidate_ready = &
                                                projector_physical_direct_model_calibrated_drop > &
                                                EPSILON(projector_physical_direct_model_drop) .AND. &
                                                projector_physical_direct_response_weighted_step <= response_trust .AND. &
                                                projector_physical_direct_response_abs_step <= response_abs_trust
                                          END IF
                                       END IF
                                       IF (projector_physical_direct_candidate_ready) THEN
                                          projector_physical_candidate_ready = .TRUE.
                                          projector_physical_direct_candidate_used = .TRUE.
                                          projector_physical_model_e(:) = projector_physical_direct_model_e(:)
                                          projector_physical_rx_e(:) = projector_physical_direct_rx_e(:)
                                          projector_physical_z_e(:) = projector_physical_direct_z_e(:)
                                          projector_physical_sol(:) = projector_physical_direct_sol(:)
                                          projector_physical_direction_scale = projector_physical_direct_direction_scale
                                          projector_physical_expected_drop = projector_physical_direct_expected_drop
                                          projector_physical_gdot_energy = projector_physical_direct_gdot_energy
                                          projector_physical_gdot_rotation = projector_physical_direct_gdot_rotation
                                          projector_physical_gdotp = projector_physical_direct_gdotp
                                          projector_physical_model_calibrated_drop = &
                                             projector_physical_direct_model_calibrated_drop
                                          projector_physical_model_drop = projector_physical_direct_model_drop
                                          projector_physical_model_gap = projector_physical_direct_model_gap
                                          projector_physical_model_gap_drop = projector_physical_direct_model_gap_drop
                                          projector_physical_model_gap_ratio = projector_physical_direct_model_gap_ratio
                                          projector_physical_model_gap_trust = projector_physical_direct_model_gap_trust
                                          projector_physical_norm_energy = projector_physical_direct_norm_energy
                                          projector_physical_norm_rotation = projector_physical_direct_norm_rotation
                                          projector_physical_normp = projector_physical_direct_normp
                                          projector_physical_norm_ratio = projector_physical_direct_norm_ratio
                                          projector_physical_response_abs_step = projector_physical_direct_response_abs_step
                                          projector_physical_response_abs_trust_ratio = &
                                             projector_physical_direct_response_abs_trust_ratio
                                          projector_physical_response_weighted_step = &
                                             projector_physical_direct_response_weighted_step
                                          projector_physical_response_weighted_trust_ratio = &
                                             projector_physical_direct_response_weighted_trust_ratio
                                          projector_physical_scale_quality = projector_physical_direct_scale_quality
                                          projector_physical_trust_quality = projector_physical_direct_trust_quality
                                          IF (projector_physical_response_weight <= EPSILON(projector_physical_response_weight)) &
                                             projector_physical_response_weight = projector_coupled_weight
                                       END IF
                                    END IF
                                    z_e(:) = z_e(:) + rx_e(:)

               gdot_energy = DOT_PRODUCT(grad_e, z_e)
               gdot_rotation = DOT_PRODUCT(grad_rot, sol)
               norm_energy = DOT_PRODUCT(z_e, z_e)
               norm_rotation = DOT_PRODUCT(sol, sol)
               gdotp = gdot_energy + gdot_rotation
               gnorm = DOT_PRODUCT(grad_e, grad_e) + DOT_PRODUCT(grad_rot, grad_rot)
               normp = norm_energy + norm_rotation
               norm_ratio = normp/MAX(gnorm, EPSILON(1.0_dp))
               IF (gdotp <= 0.0_dp .OR. gdotp /= gdotp .OR. normp /= normp) THEN
                  attempt_failed = .TRUE.
                  reject_reason = reject_nondescent
               END IF
               IF (.NOT. attempt_failed .AND. norm_ratio > max_norm_ratio) THEN
                  attempt_failed = .TRUE.
                  reject_reason = reject_norm
               END IF
               IF (.NOT. attempt_failed) THEN
                  ! The generic OT line searches assume that the stored search vector z has
                  ! z.z = g.z.  Rescale the external Schur direction to preserve that slope.
                  direction_scale = gdotp/MAX(normp, EPSILON(normp))
                  IF (direction_scale <= 0.0_dp .OR. direction_scale /= direction_scale) THEN
                     attempt_failed = .TRUE.
                     reject_reason = reject_nondescent
                  ELSE
                     scale_quality = kpoint_ot_schur_scale_quality(direction_scale, schur_min_coupled_scale)
                     z_e(:) = direction_scale*z_e(:)
                     sol(:) = direction_scale*sol(:)
                     model_drop = 0.5_dp*direction_scale*gdotp
                     model_calibrated_drop = model_drop
                     IF (schur_attempt <= 2) model_calibrated_drop = &
                        split_quality_for_calibration*model_drop
                     response_abs_step = direction_scale*MAXVAL(ABS(rx_e))
                     CALL para_env_inter_kp%max(response_abs_step)
                     response_step_max = response_abs_step
                     response_weighted_step = direction_scale* &
                                              kpoint_ot_weighted_response_step(rx_e, curv, curv_projection_sum)
                     IF (response_trust > EPSILON(response_trust)) &
                        response_weighted_trust_ratio = response_weighted_step/response_trust
                     IF (response_abs_trust > EPSILON(response_abs_trust)) &
                        response_abs_trust_ratio = response_abs_step/response_abs_trust
                     IF (response_model_available) THEN
                        model_gap = direction_scale*MAXVAL(ABS(rx_model - rx_e))
                        model_gap_drop = kpoint_ot_response_mismatch_drop(rx_e, rx_model, curv, direction_scale)
                        IF (schur_attempt <= 2) THEN
                           model_calibrated_drop = &
                              MAX(0.0_dp, split_quality_for_calibration*model_drop - model_gap_drop)
                        ELSE
                           model_calibrated_drop = MAX(0.0_dp, model_drop - model_gap_drop)
                        END IF
                        model_gap_trust = schur_model_gap_factor*response_trust
                        IF (model_gap_trust > EPSILON(model_gap_trust)) model_gap_ratio = model_gap/model_gap_trust
                        IF (schur_attempt == 3 .AND. exact_attempt_response .AND. &
                            projector_coupled_weight > 0.0_dp .AND. &
                            (local_ot_env(1)%schur_projector_mode_quality_samples < &
                             ot_kp_schur_calibration_min_samples .OR. split_candidate_ready .OR. &
                             logger%iter_info%print_level >= high_print_level)) THEN
                           CALL audit_kpoint_ot_projector_modes(rot_hessian, rot_hessian_projector, sol, &
                                                                response_exact_raw, response_exact, curv, &
                                                                model_calibrated_drop, projector_coupled_weight, &
                                                                projector_mode_history_quality, dense_residual_quality, &
                                                                projector_mode_quality, projector_mode_negative_fraction, &
                                                                projector_mode_damped_abs_avg, projector_mode_eval_min, &
                                                                projector_mode_eval_max, &
                                                                projector_mode_step_negative_fraction, &
                                                                projector_mode_response_negative_fraction, &
                                                                projector_mode_response_drop, &
                                                                projector_mode_response_drop_ratio, &
                                                                projector_mode_trust_quality)
                           model_calibrated_drop = projector_mode_trust_quality*model_calibrated_drop
                        END IF
                     END IF
                     IF (schur_attempt <= 2 .AND. .NOT. attempt_failed .AND. &
                         model_calibrated_drop <= EPSILON(model_drop)) THEN
                        attempt_failed = .TRUE.
                        trust_rejected = .TRUE.
                        reject_reason = reject_model_mismatch
                     END IF
                     IF (schur_attempt <= 3 .AND. direction_scale < schur_min_coupled_scale) THEN
                        IF (schur_attempt == 3 .AND. exact_attempt_response .AND. &
                            projector_coupled_weight > 0.0_dp) THEN
                           model_calibrated_drop = scale_quality*model_calibrated_drop
                        ELSE
                           attempt_failed = .TRUE.
                           reject_reason = reject_scale
                        END IF
                     END IF
                     IF (schur_attempt <= 3 .AND. .NOT. attempt_failed) THEN
                        rx_model(:) = MATMUL(response_occ - response, sol)
                        occgap_shadow_abs_step = MAXVAL(ABS(rx_model))
                        CALL para_env_inter_kp%max(occgap_shadow_abs_step)
                        occgap_shadow_step = kpoint_ot_weighted_response_step(rx_model, occ_curv, &
                                                                              occ_curv_projection_sum)
                        IF (response_trust > EPSILON(response_trust)) &
                           occgap_shadow_gap_ratio = MAX(occgap_shadow_gap_ratio, &
                                                         occgap_shadow_step/response_trust)
                        IF (response_abs_trust > EPSILON(response_abs_trust)) &
                           occgap_shadow_gap_ratio = MAX(occgap_shadow_gap_ratio, &
                                                         occgap_shadow_abs_step/response_abs_trust)
                        occgap_shadow_quality = 1.0_dp/MAX(1.0_dp, occgap_shadow_gap_ratio)
                        occgap_shadow_drop = kpoint_ot_schur_expected_drop(model_calibrated_drop, &
                                                                           1.0_dp, occgap_shadow_quality)
                     END IF
                     ! Trust the Schur response in the fixed-N Mermin metric.  A
                     ! single weakly-curved band may have a large absolute energy
                     ! response; the broad outer guard below catches real Fermi-window
                     ! escapes without rejecting curvature-weighted rotation steps.
                     IF (schur_attempt <= 3 .AND. .NOT. attempt_failed .AND. &
                         response_weighted_step > response_trust) THEN
                        attempt_failed = .TRUE.
                        trust_rejected = .TRUE.
                        reject_reason = reject_response
                     END IF
                     IF (schur_attempt <= 3 .AND. .NOT. attempt_failed .AND. &
                         response_abs_step > response_abs_trust) THEN
                        attempt_failed = .TRUE.
                        trust_rejected = .TRUE.
                        reject_reason = reject_response
                     END IF
                     IF (schur_attempt <= 2 .AND. .NOT. attempt_failed .AND. &
                         model_gap > model_gap_trust) THEN
                        attempt_failed = .TRUE.
                        trust_rejected = .TRUE.
                        reject_reason = reject_model_gap
                     END IF
                     IF (schur_attempt <= 2 .AND. .NOT. attempt_failed .AND. &
                         .NOT. split_candidate_enabled) THEN
                        attempt_failed = .TRUE.
                        reject_reason = reject_quality_gate
                     END IF
                  END IF
               END IF
            END IF

            IF (schur_attempt == 1 .AND. .NOT. attempt_failed) THEN
               constrained_candidate_ready = .TRUE.
               constrained_sol(:) = sol(:)
               constrained_z_e(:) = z_e(:)
               constrained_gdotp = gdotp
               constrained_gdot_energy = gdot_energy
               constrained_gdot_rotation = gdot_rotation
               constrained_norm_energy = norm_energy
               constrained_norm_rotation = norm_rotation
               constrained_normp = normp
               constrained_scale_quality = scale_quality
               constrained_occgap_shadow_drop = occgap_shadow_drop
               constrained_occgap_shadow_gap_ratio = occgap_shadow_gap_ratio
               ! Delay the constrained candidate until the raw Split candidate for
               ! the same search point gives an absolute Mermin-drop comparator.
               attempt_failed = .TRUE.
               reject_reason = reject_quality_gate
            END IF

            IF (schur_attempt == 2) THEN
               coupled_diag_shift_avg = diag_shift_sum/REAL(MAX(1, nrot), KIND=dp)
               coupled_direction_scale = direction_scale
               coupled_gdot_energy = gdot_energy
               coupled_gdot_rotation = gdot_rotation
               coupled_model_drop = model_drop
               coupled_model_calibrated_drop = model_calibrated_drop
               coupled_model_gap = model_gap
               coupled_model_gap_drop = model_gap_drop
               coupled_model_gap_ratio = model_gap_ratio
               coupled_model_gap_trust = model_gap_trust
               coupled_norm_energy = norm_energy
               coupled_norm_rotation = norm_rotation
               coupled_norm_ratio = norm_ratio
               coupled_scale_quality = scale_quality
               coupled_occ_response_diag_avg = occ_response_diag_sum/REAL(MAX(1, nrot), KIND=dp)
               coupled_rot_diag_avg = rot_diag_sum/REAL(MAX(1, nrot), KIND=dp)
               coupled_tether_response_diag_avg = tether_response_diag_sum/REAL(MAX(1, nrot), KIND=dp)
               coupled_reject_reason = reject_reason
               coupled_response_abs_step = response_abs_step
               coupled_response_abs_trust_ratio = response_abs_trust_ratio
               coupled_response_step = response_step_max
               coupled_response_trust_ratio = 0.0_dp
               IF (response_trust > EPSILON(response_trust)) &
                  coupled_response_trust_ratio = response_step_max/response_trust
               coupled_response_weighted_step = response_weighted_step
               coupled_response_weighted_trust_ratio = response_weighted_trust_ratio
            ELSE IF (schur_attempt == 1) THEN
               constrained_diag_shift_avg = diag_shift_sum/REAL(MAX(1, nrot), KIND=dp)
               constrained_direction_scale = direction_scale
               constrained_gdot_energy = gdot_energy
               constrained_gdot_rotation = gdot_rotation
               constrained_model_drop = model_drop
               constrained_model_calibrated_drop = model_calibrated_drop
               constrained_model_gap = model_gap
               constrained_model_gap_drop = model_gap_drop
               constrained_model_gap_ratio = model_gap_ratio
               constrained_model_gap_trust = model_gap_trust
               constrained_norm_energy = norm_energy
               constrained_norm_rotation = norm_rotation
               constrained_norm_ratio = norm_ratio
               constrained_scale_quality = scale_quality
               constrained_occ_response_diag_avg = occ_response_diag_sum/REAL(MAX(1, nrot), KIND=dp)
               constrained_response_penalty_diag_avg = response_penalty_diag_sum/REAL(MAX(1, nrot), KIND=dp)
               constrained_rot_diag_avg = rot_diag_sum/REAL(MAX(1, nrot), KIND=dp)
               constrained_tether_response_diag_avg = tether_response_diag_sum/REAL(MAX(1, nrot), KIND=dp)
               constrained_reject_reason = reject_reason
               constrained_response_abs_step = response_abs_step
               constrained_response_abs_trust_ratio = response_abs_trust_ratio
               constrained_response_step = response_step_max
               constrained_response_trust_ratio = 0.0_dp
               IF (response_trust > EPSILON(response_trust)) &
                  constrained_response_trust_ratio = response_step_max/response_trust
               constrained_response_weighted_step = response_weighted_step
               constrained_response_weighted_trust_ratio = response_weighted_trust_ratio
               IF (new_search_point .AND. split_candidate_enabled .AND. &
                   model_calibrated_drop > EPSILON(model_drop) .AND. &
                   (reject_reason == 0 .OR. reject_reason == reject_quality_gate)) THEN
                  local_ot_env(1)%schur_constrained_pending = .TRUE.
                  local_ot_env(1)%schur_constrained_drop = &
                     local_ot_env(1)%schur_constrained_drop + model_drop
                  local_ot_env(1)%schur_constrained_calibrated_drop = &
                     local_ot_env(1)%schur_constrained_calibrated_drop + model_calibrated_drop
                  local_ot_env(1)%schur_constrained_gap_ratio = &
                     MAX(local_ot_env(1)%schur_constrained_gap_ratio, model_gap_ratio)
               END IF
            END IF

            IF (schur_attempt == 2 .AND. constrained_candidate_ready .AND. &
                model_calibrated_drop > EPSILON(model_calibrated_drop)) THEN
               split_trust_quality = kpoint_ot_schur_trust_quality(response_weighted_trust_ratio, &
                                                                   response_abs_trust_ratio, &
                                                                   model_gap_ratio)
               constrained_trust_quality = &
                  kpoint_ot_schur_trust_quality(constrained_response_weighted_trust_ratio, &
                                                constrained_response_abs_trust_ratio, &
                                                constrained_model_gap_ratio)
               split_expected_drop = kpoint_ot_schur_expected_drop(model_calibrated_drop, &
                                                                   split_quality_for_selection, &
                                                                   split_trust_quality)
               constrained_expected_drop = &
                  kpoint_ot_schur_expected_drop(constrained_model_calibrated_drop, &
                                               constrained_quality_for_selection, &
                                               constrained_trust_quality)
               IF (split_expected_drop > EPSILON(split_expected_drop)) &
                  constrained_expected_drop_ratio = constrained_expected_drop/split_expected_drop
               IF (.NOT. attempt_failed .AND. constrained_candidate_enabled .AND. &
                   constrained_expected_drop_ratio >= ot_kp_schur_constrained_drop_advantage_ratio) THEN
                  constrained_candidate_selected = .TRUE.
                  z_e(:) = constrained_z_e(:)
                  sol(:) = constrained_sol(:)
                  direction_scale = constrained_direction_scale
                  gdotp = constrained_gdotp
                  gdot_energy = constrained_gdot_energy
                  gdot_rotation = constrained_gdot_rotation
                  norm_energy = constrained_norm_energy
                  norm_rotation = constrained_norm_rotation
                  normp = constrained_normp
                  norm_ratio = constrained_norm_ratio
                  scale_quality = constrained_scale_quality
                  model_drop = constrained_model_drop
                  model_calibrated_drop = constrained_model_calibrated_drop
                  model_gap_ratio = constrained_model_gap_ratio
                  occgap_shadow_drop = constrained_occgap_shadow_drop
                  occgap_shadow_gap_ratio = constrained_occgap_shadow_gap_ratio
                  response_step_max = constrained_response_step
                  response_abs_step = constrained_response_abs_step
                  response_abs_trust_ratio = constrained_response_abs_trust_ratio
                  response_weighted_step = constrained_response_weighted_step
                  response_weighted_trust_ratio = constrained_response_weighted_trust_ratio
                  constrained_reject_reason = 0
               END IF
            END IF
            IF (schur_attempt == 2 .AND. .NOT. attempt_failed) THEN
               IF (.NOT. constrained_candidate_selected) THEN
                  split_trust_quality = kpoint_ot_schur_trust_quality(response_weighted_trust_ratio, &
                                                                      response_abs_trust_ratio, &
                                                                      model_gap_ratio)
                  split_expected_drop = &
                     kpoint_ot_schur_expected_drop(model_calibrated_drop, split_quality_for_selection, &
                                                  split_trust_quality)
               END IF
               split_candidate_ready = .TRUE.
               split_selected_attempt = MERGE(1, 2, constrained_candidate_selected)
               split_z_e(:) = z_e(:)
               split_sol(:) = sol(:)
               split_direction_scale = direction_scale
               split_gdotp = gdotp
               split_gdot_energy = gdot_energy
               split_gdot_rotation = gdot_rotation
               split_norm_energy = norm_energy
               split_norm_rotation = norm_rotation
               split_normp = normp
               split_norm_ratio = norm_ratio
               split_scale_quality = scale_quality
               split_model_drop = model_drop
               split_model_calibrated_drop = model_calibrated_drop
               split_model_gap_ratio = model_gap_ratio
               split_occgap_shadow_drop = occgap_shadow_drop
               split_occgap_shadow_gap_ratio = occgap_shadow_gap_ratio
               split_density_projector_fraction = density_projector_fraction
               split_diag_shift_max = diag_shift_max
               split_diag_shift_sum = diag_shift_sum
               split_residual_delta = residual_delta
               split_response_step = response_step_max
               split_response_abs_step = response_abs_step
               split_response_abs_trust = response_abs_trust
               split_response_abs_trust_ratio = response_abs_trust_ratio
               split_response_trust = response_trust
               split_response_trust_factor = response_trust_factor
               split_response_weighted_step = response_weighted_step
               split_response_weighted_trust_ratio = response_weighted_trust_ratio
               split_rot_floor_min = rot_floor_min
               split_schur_diag_min = schur_diag_min
               attempt_failed = .TRUE.
               reject_reason = reject_quality_gate
            ELSE IF (schur_attempt == 3 .AND. split_candidate_ready) THEN
               dense_direction_scale = direction_scale
               dense_gdot_energy = gdot_energy
               dense_gdot_rotation = gdot_rotation
               dense_model_gap = model_gap
               dense_model_gap_drop = model_gap_drop
               dense_model_gap_ratio = model_gap_ratio
               dense_model_gap_trust = model_gap_trust
               dense_norm_energy = norm_energy
               dense_norm_rotation = norm_rotation
               dense_projector_coupled_weight = projector_coupled_weight
               dense_reject_reason = reject_reason
               dense_scale_quality = scale_quality
               dense_trust_quality = kpoint_ot_schur_trust_quality(response_weighted_trust_ratio, &
                                                                   response_abs_trust_ratio, &
                                                                   model_gap_ratio)
               IF (.NOT. attempt_failed) &
                  dense_expected_drop = kpoint_ot_schur_expected_drop(model_calibrated_drop, &
                                                                     dense_quality_for_selection, &
                                                                     dense_trust_quality)
               IF (.NOT. attempt_failed .AND. local_ot_env(1)%schur_projector_response_candidate_enabled .AND. &
                   projector_response_candidate_ready .AND. &
                   projector_response_expected_drop > EPSILON(projector_response_expected_drop) .AND. &
                   projector_response_expected_drop >= &
                   ot_kp_schur_projector_response_advantage_ratio*MAX(dense_expected_drop, EPSILON(dense_expected_drop))) THEN
                  projector_response_candidate_selected = .TRUE.
                  z_e(:) = projector_response_z_e(:)
                  sol(:) = projector_response_sol(:)
                  direction_scale = projector_response_direction_scale
                  gdotp = projector_response_gdotp
                  gdot_energy = projector_response_gdot_energy
                  gdot_rotation = projector_response_gdot_rotation
                  norm_energy = projector_response_norm_energy
                  norm_rotation = projector_response_norm_rotation
                  normp = projector_response_normp
                  norm_ratio = projector_response_norm_ratio
                  scale_quality = projector_response_scale_quality
                  model_drop = projector_response_model_drop
                  model_calibrated_drop = projector_response_model_calibrated_drop
                  model_gap = projector_response_model_gap
                  model_gap_drop = projector_response_model_gap_drop
                  model_gap_ratio = projector_response_model_gap_ratio
                  model_gap_trust = projector_response_model_gap_trust
                  response_abs_step = projector_response_response_abs_step
                  response_abs_trust_ratio = projector_response_response_abs_trust_ratio
                  response_weighted_step = projector_response_response_weighted_step
                     response_weighted_trust_ratio = projector_response_response_weighted_trust_ratio
                     dense_expected_drop = projector_response_expected_drop
                  END IF
                  IF (.NOT. attempt_failed .AND. local_ot_env(1)%schur_projector_physical_candidate_enabled .AND. &
                      projector_physical_candidate_ready .AND. &
                      projector_physical_expected_drop > EPSILON(projector_physical_expected_drop) .AND. &
                      projector_physical_expected_drop >= &
                      ot_kp_schur_projector_response_advantage_ratio*MAX(dense_expected_drop, EPSILON(dense_expected_drop))) THEN
                     projector_physical_candidate_selected = .TRUE.
                     projector_physical_direct_candidate_selected = projector_physical_direct_candidate_used
                     projector_response_candidate_selected = .FALSE.
                     z_e(:) = projector_physical_z_e(:)
                     sol(:) = projector_physical_sol(:)
                     direction_scale = projector_physical_direction_scale
                     gdotp = projector_physical_gdotp
                     gdot_energy = projector_physical_gdot_energy
                     gdot_rotation = projector_physical_gdot_rotation
                     norm_energy = projector_physical_norm_energy
                     norm_rotation = projector_physical_norm_rotation
                     normp = projector_physical_normp
                     norm_ratio = projector_physical_norm_ratio
                     scale_quality = projector_physical_scale_quality
                     model_drop = projector_physical_model_drop
                     model_calibrated_drop = projector_physical_model_calibrated_drop
                     model_gap = projector_physical_model_gap
                     model_gap_drop = projector_physical_model_gap_drop
                     model_gap_ratio = projector_physical_model_gap_ratio
                     model_gap_trust = projector_physical_model_gap_trust
                     response_abs_step = projector_physical_response_abs_step
                     response_abs_trust_ratio = projector_physical_response_abs_trust_ratio
                     response_weighted_step = projector_physical_response_weighted_step
                     response_weighted_trust_ratio = projector_physical_response_weighted_trust_ratio
                     dense_expected_drop = projector_physical_expected_drop
                  END IF
                  split_dense_expected_ratio = 0.0_dp
               IF (dense_expected_drop > EPSILON(dense_expected_drop)) &
                  split_dense_expected_ratio = split_expected_drop/dense_expected_drop
               split_candidate_selected = attempt_failed .OR. split_expected_drop > dense_expected_drop
               IF (new_search_point .AND. split_candidate_ready .AND. .NOT. split_candidate_selected .AND. &
                   split_model_drop > EPSILON(split_model_drop)) THEN
                  local_ot_env(1)%schur_split_shadow_pending = .TRUE.
                  local_ot_env(1)%schur_split_shadow_drop = &
                     local_ot_env(1)%schur_split_shadow_drop + split_model_drop
                  local_ot_env(1)%schur_split_shadow_calibrated_drop = &
                     local_ot_env(1)%schur_split_shadow_calibrated_drop + split_model_calibrated_drop
                  END IF
                  IF (split_candidate_selected) THEN
                     attempt_failed = .FALSE.
                     projector_physical_candidate_selected = .FALSE.
                     projector_physical_direct_candidate_selected = .FALSE.
                     projector_response_candidate_selected = .FALSE.
                     z_e(:) = split_z_e(:)
                  sol(:) = split_sol(:)
                  direction_scale = split_direction_scale
                  gdotp = split_gdotp
                  gdot_energy = split_gdot_energy
                  gdot_rotation = split_gdot_rotation
                  norm_energy = split_norm_energy
                  norm_rotation = split_norm_rotation
                  normp = split_normp
                  norm_ratio = split_norm_ratio
                  scale_quality = split_scale_quality
                  model_drop = split_model_drop
                  model_calibrated_drop = split_model_calibrated_drop
                  model_gap_ratio = split_model_gap_ratio
                  occgap_shadow_drop = split_occgap_shadow_drop
                  occgap_shadow_gap_ratio = split_occgap_shadow_gap_ratio
                  density_projector_fraction = split_density_projector_fraction
                  diag_shift_max = split_diag_shift_max
                  diag_shift_sum = split_diag_shift_sum
                  residual_delta = split_residual_delta
                  response_step_max = split_response_step
                  response_abs_step = split_response_abs_step
                  response_abs_trust = split_response_abs_trust
                  response_abs_trust_ratio = split_response_abs_trust_ratio
                  response_trust = split_response_trust
                  response_trust_factor = split_response_trust_factor
                  response_weighted_step = split_response_weighted_step
                  response_weighted_trust_ratio = split_response_weighted_trust_ratio
                  rot_floor_min = split_rot_floor_min
                  schur_diag_min = split_schur_diag_min
               END IF
            END IF

            IF (schur_attempt == 3 .AND. .NOT. split_candidate_ready) THEN
               dense_direction_scale = direction_scale
               dense_gdot_energy = gdot_energy
               dense_gdot_rotation = gdot_rotation
               dense_model_gap = model_gap
               dense_model_gap_drop = model_gap_drop
               dense_model_gap_ratio = model_gap_ratio
               dense_model_gap_trust = model_gap_trust
               dense_norm_energy = norm_energy
               dense_norm_rotation = norm_rotation
               dense_projector_coupled_weight = projector_coupled_weight
               dense_reject_reason = reject_reason
               dense_scale_quality = scale_quality
               dense_trust_quality = kpoint_ot_schur_trust_quality(response_weighted_trust_ratio, &
                                                                  response_abs_trust_ratio, &
                                                                  model_gap_ratio)
               IF (.NOT. attempt_failed) &
                  dense_expected_drop = kpoint_ot_schur_expected_drop(model_calibrated_drop, &
                                                                     dense_quality_for_selection, &
                                                                     dense_trust_quality)
               IF (.NOT. attempt_failed .AND. local_ot_env(1)%schur_projector_response_candidate_enabled .AND. &
                   projector_response_candidate_ready .AND. &
                   projector_response_expected_drop > EPSILON(projector_response_expected_drop) .AND. &
                   projector_response_expected_drop >= &
                   ot_kp_schur_projector_response_advantage_ratio*MAX(dense_expected_drop, EPSILON(dense_expected_drop))) THEN
                  projector_response_candidate_selected = .TRUE.
                  z_e(:) = projector_response_z_e(:)
                  sol(:) = projector_response_sol(:)
                  direction_scale = projector_response_direction_scale
                  gdotp = projector_response_gdotp
                  gdot_energy = projector_response_gdot_energy
                  gdot_rotation = projector_response_gdot_rotation
                  norm_energy = projector_response_norm_energy
                  norm_rotation = projector_response_norm_rotation
                  normp = projector_response_normp
                  norm_ratio = projector_response_norm_ratio
                  scale_quality = projector_response_scale_quality
                  model_drop = projector_response_model_drop
                  model_calibrated_drop = projector_response_model_calibrated_drop
                  model_gap = projector_response_model_gap
                  model_gap_drop = projector_response_model_gap_drop
                  model_gap_ratio = projector_response_model_gap_ratio
                  model_gap_trust = projector_response_model_gap_trust
                  response_abs_step = projector_response_response_abs_step
                  response_abs_trust_ratio = projector_response_response_abs_trust_ratio
                  response_weighted_step = projector_response_response_weighted_step
                     response_weighted_trust_ratio = projector_response_response_weighted_trust_ratio
                     dense_expected_drop = projector_response_expected_drop
                  END IF
                  IF (.NOT. attempt_failed .AND. local_ot_env(1)%schur_projector_physical_candidate_enabled .AND. &
                      projector_physical_candidate_ready .AND. &
                      projector_physical_expected_drop > EPSILON(projector_physical_expected_drop) .AND. &
                      projector_physical_expected_drop >= &
                      ot_kp_schur_projector_response_advantage_ratio*MAX(dense_expected_drop, EPSILON(dense_expected_drop))) THEN
                     projector_physical_candidate_selected = .TRUE.
                     projector_physical_direct_candidate_selected = projector_physical_direct_candidate_used
                     projector_response_candidate_selected = .FALSE.
                     z_e(:) = projector_physical_z_e(:)
                     sol(:) = projector_physical_sol(:)
                     direction_scale = projector_physical_direction_scale
                     gdotp = projector_physical_gdotp
                     gdot_energy = projector_physical_gdot_energy
                     gdot_rotation = projector_physical_gdot_rotation
                     norm_energy = projector_physical_norm_energy
                     norm_rotation = projector_physical_norm_rotation
                     normp = projector_physical_normp
                     norm_ratio = projector_physical_norm_ratio
                     scale_quality = projector_physical_scale_quality
                     model_drop = projector_physical_model_drop
                     model_calibrated_drop = projector_physical_model_calibrated_drop
                     model_gap = projector_physical_model_gap
                     model_gap_drop = projector_physical_model_gap_drop
                     model_gap_ratio = projector_physical_model_gap_ratio
                     model_gap_trust = projector_physical_model_gap_trust
                     response_abs_step = projector_physical_response_abs_step
                     response_abs_trust_ratio = projector_physical_response_abs_trust_ratio
                     response_weighted_step = projector_physical_response_weighted_step
                     response_weighted_trust_ratio = projector_physical_response_weighted_trust_ratio
                     dense_expected_drop = projector_physical_expected_drop
                  END IF
               END IF

            IF (.NOT. attempt_failed) THEN
               solved = .TRUE.
               selected_attempt = schur_attempt
               IF (split_candidate_selected) selected_attempt = split_selected_attempt
               schur_used_attempt = selected_attempt
               accepted_model_drop = model_drop
               accepted_model_calibrated_drop = model_calibrated_drop
               accepted_model_gap_ratio = model_gap_ratio
               accepted_occgap_shadow_drop = occgap_shadow_drop
               accepted_occgap_shadow_gap_ratio = occgap_shadow_gap_ratio
               accepted_projector_coupled_weight = projector_coupled_weight
               IF (split_candidate_selected) accepted_projector_coupled_weight = 1.0_dp
               IF (selected_attempt == 3 .AND. projector_coupled_weight > 0.0_dp .AND. &
                   (projector_mode_response_drop > 0.0_dp .OR. &
                    projector_mode_response_drop_ratio > 0.0_dp)) THEN
                  accepted_projector_mode_drop = model_calibrated_drop
                  accepted_projector_mode_response_drop = projector_mode_response_drop
                  accepted_projector_mode_response_negative_fraction = projector_mode_response_negative_fraction
                  accepted_projector_mode_step_negative_fraction = projector_mode_step_negative_fraction
                  accepted_projector_mode_trust_quality = projector_mode_trust_quality
               END IF
                  IF (selected_attempt == 3 .AND. projector_response_candidate_ready .AND. &
                      .NOT. projector_physical_candidate_selected .AND. &
                      projector_response_model_calibrated_drop > EPSILON(projector_response_model_calibrated_drop)) THEN
                     accepted_projector_response_drop = projector_response_model_calibrated_drop
                     accepted_projector_response_weight = projector_mode_response_weight
                     accepted_projector_response_negative = projector_mode_response_correction_negative
                  END IF
                  EXIT
            END IF
         END DO
         failed = .NOT. solved

         IF (new_search_point .AND. accepted_model_drop > EPSILON(accepted_model_drop)) THEN
            IF (.NOT. local_ot_env(1)%schur_model_pending) THEN
               local_ot_env(1)%schur_model_attempt = schur_used_attempt
            ELSE IF (local_ot_env(1)%schur_model_attempt /= schur_used_attempt) THEN
               local_ot_env(1)%schur_model_attempt = -1
            END IF
            local_ot_env(1)%schur_model_pending = .TRUE.
            local_ot_env(1)%schur_model_drop = local_ot_env(1)%schur_model_drop + accepted_model_drop
            local_ot_env(1)%schur_model_calibrated_drop = &
               local_ot_env(1)%schur_model_calibrated_drop + accepted_model_calibrated_drop
            local_ot_env(1)%schur_model_gap_ratio = MAX(local_ot_env(1)%schur_model_gap_ratio, &
                                                        accepted_model_gap_ratio)
            local_ot_env(1)%schur_model_delta = MAX(local_ot_env(1)%schur_model_delta, residual_delta)
            IF (accepted_occgap_shadow_drop > EPSILON(accepted_occgap_shadow_drop)) THEN
               local_ot_env(1)%schur_occgap_shadow_pending = .TRUE.
               local_ot_env(1)%schur_occgap_shadow_drop = &
                  local_ot_env(1)%schur_occgap_shadow_drop + accepted_occgap_shadow_drop
               local_ot_env(1)%schur_occgap_shadow_gap_ratio = &
                  MAX(local_ot_env(1)%schur_occgap_shadow_gap_ratio, accepted_occgap_shadow_gap_ratio)
            END IF
            IF (accepted_projector_mode_drop > EPSILON(accepted_projector_mode_drop)) THEN
               local_ot_env(1)%schur_projector_mode_pending = .TRUE.
               local_ot_env(1)%schur_projector_mode_drop = &
                  local_ot_env(1)%schur_projector_mode_drop + accepted_projector_mode_drop
               local_ot_env(1)%schur_projector_mode_response_drop = &
                  local_ot_env(1)%schur_projector_mode_response_drop + accepted_projector_mode_response_drop
               local_ot_env(1)%schur_projector_mode_step_negative = &
                  local_ot_env(1)%schur_projector_mode_step_negative + &
                  accepted_projector_mode_step_negative_fraction*accepted_projector_mode_drop
               local_ot_env(1)%schur_projector_mode_response_negative = &
                  local_ot_env(1)%schur_projector_mode_response_negative + &
                  accepted_projector_mode_response_negative_fraction*accepted_projector_mode_response_drop
               local_ot_env(1)%schur_projector_mode_trust = &
                  local_ot_env(1)%schur_projector_mode_trust + &
                  accepted_projector_mode_trust_quality*accepted_projector_mode_drop
            END IF
            IF (accepted_projector_response_drop > EPSILON(accepted_projector_response_drop)) THEN
               local_ot_env(1)%schur_projector_response_pending = .TRUE.
               local_ot_env(1)%schur_projector_response_drop = &
                  local_ot_env(1)%schur_projector_response_drop + accepted_projector_response_drop
               local_ot_env(1)%schur_projector_response_weight = &
                  MAX(local_ot_env(1)%schur_projector_response_weight, accepted_projector_response_weight)
               local_ot_env(1)%schur_projector_response_correction_drop = &
                  local_ot_env(1)%schur_projector_response_correction_drop + &
                  projector_mode_response_correction_drop
                  local_ot_env(1)%schur_projector_response_negative = &
                     MAX(local_ot_env(1)%schur_projector_response_negative, accepted_projector_response_negative)
               END IF
               IF (projector_physical_candidate_ready .AND. &
                   projector_physical_model_calibrated_drop > EPSILON(projector_physical_model_calibrated_drop)) THEN
                  local_ot_env(1)%schur_projector_physical_pending = .TRUE.
                  local_ot_env(1)%schur_projector_physical_drop = &
                     local_ot_env(1)%schur_projector_physical_drop + projector_physical_model_calibrated_drop
                  local_ot_env(1)%schur_projector_physical_weight = &
                     MAX(local_ot_env(1)%schur_projector_physical_weight, projector_physical_response_weight)
                  local_ot_env(1)%schur_projector_physical_correction_drop = &
                     local_ot_env(1)%schur_projector_physical_correction_drop + &
                     projector_physical_response_correction_drop
                  local_ot_env(1)%schur_projector_physical_negative = &
                     MAX(local_ot_env(1)%schur_projector_physical_negative, &
                         projector_physical_response_correction_negative)
               END IF
            END IF

         IF (.NOT. failed) THEN
            DO local_kpoint = 1, SIZE(kpoints%kp_env)
               kp => kpoints%kp_env(local_kpoint)%kpoint_env
               ikpoint = kp%nkpoint
               DO ispin = 1, nspin_ot
                  group = MERGE(ispin, 1, ngroups == nspin_ot)
                  IF (group /= target_group) CYCLE
                  local_channel = qs_ot_channel_index(ispin, ikpoint, nspin_ot) - first_channel + 1
                  nmo = channel_nmo(local_channel)
                  DO imo = 1, nmo
                     ie = e_offset(local_channel) + imo - 1
                     local_ot_env(local_channel)%ener_gx(imo) = z_e(ie)
                  END DO
               END DO
            END DO

            DO local_channel = 1, nchannels
               nmo = channel_nmo(local_channel)
               IF (nmo <= 0) CYCLE
               ALLOCATE (grad_re(nmo, nmo), grad_im(nmo, nmo))
               CALL dbcsr_square_to_dense(local_ot_env(local_channel)%rot_mat_gx, grad_re)
               CALL dbcsr_square_to_dense(local_ot_env(local_channel)%rot_mat_gx_im, grad_im)
               CALL local_ot_env(local_channel)%para_env%sum(grad_re)
               CALL local_ot_env(local_channel)%para_env%sum(grad_im)
               DO r = 1, nrot
                  IF (var_channel(r) /= local_channel) CYCLE
                  imo = var_i(r)
                  jmo = var_j(r)
                  SELECT CASE (var_kind(r))
                  CASE (1)
                     grad_re(imo, jmo) = sol(r)
                     grad_re(jmo, imo) = -sol(r)
                  CASE (2)
                     grad_im(imo, jmo) = sol(r)
                     grad_im(jmo, imo) = sol(r)
                  END SELECT
               END DO
               CALL dense_to_dbcsr_square(local_ot_env(local_channel)%rot_mat_gx, grad_re)
               CALL dense_to_dbcsr_square(local_ot_env(local_channel)%rot_mat_gx_im, grad_im)
               DEALLOCATE (grad_re, grad_im)
            END DO
            group_done = .TRUE.
            IF (logger%iter_info%print_level >= high_print_level .AND. logger%para_env%is_source()) THEN
               diag_shift_avg = diag_shift_sum/REAL(MAX(1, nrot), KIND=dp)
               response_abs_trust_ratio = 0.0_dp
               IF (response_abs_trust > EPSILON(response_abs_trust)) &
                  response_abs_trust_ratio = response_abs_step/response_abs_trust
               response_trust_ratio = 0.0_dp
               IF (response_trust > EPSILON(response_trust)) &
                  response_trust_ratio = response_step_max/response_trust
               response_weighted_trust_ratio = 0.0_dp
               IF (response_trust > EPSILON(response_trust)) &
                  response_weighted_trust_ratio = response_weighted_step/response_trust
               WRITE (cp_logger_get_default_unit_nr(logger), '(A,4(1X,I6),19(1X,ES16.8),1(1X,L1))') &
                  " OT Schur ne nrot att scale gdotp normp drop resp trust rtrust wresp wratio atrust aratio fac "// &
                  "delta proj dsum davg dmax dmin floor reject:", &
                  target_group, ne, nrot, schur_used_attempt, direction_scale, gdotp, normp, &
                  model_drop, response_step_max, response_trust, response_trust_ratio, response_weighted_step, &
                  response_weighted_trust_ratio, response_abs_trust, response_abs_trust_ratio, &
                  response_trust_factor, residual_delta, density_projector_fraction, diag_shift_sum, diag_shift_avg, &
                  diag_shift_max, schur_diag_min, rot_floor_min, trust_rejected
               IF (split_candidate_diagnostic) THEN
                  WRITE (cp_logger_get_default_unit_nr(logger), '(A,3(1X,I6),3(1X,L1),27(1X,ES16.8))') &
                  " OT Schur split audit group code samples gate new energy scale drop cdrop pred meas "// &
                  "qratio qema cqual cema calib qcal resp rtrust rratio wresp wratio atrust aratio mgap mdrop "// &
                  "mtrust mgratio norm "// &
                  "davg rdiag odiag tdiag:", &
                  target_group, coupled_reject_reason, local_ot_env(1)%schur_model_quality_samples, &
                  split_candidate_enabled, new_search_point, local_ot_env(1)%energy_only, &
                  coupled_direction_scale, coupled_model_drop, coupled_model_calibrated_drop, &
                  local_ot_env(1)%schur_model_drop, &
                  local_ot_env(1)%schur_model_measured_drop, local_ot_env(1)%schur_model_quality, &
                  local_ot_env(1)%schur_model_quality_ema, local_ot_env(1)%schur_split_calibrated_quality, &
                  local_ot_env(1)%schur_split_calibrated_quality_ema, &
                  split_curvature_scale, split_quality_for_calibration, &
                  coupled_response_step, response_trust, &
                  coupled_response_trust_ratio, coupled_response_weighted_step, &
                  coupled_response_weighted_trust_ratio, response_abs_trust, coupled_response_abs_trust_ratio, &
                  coupled_model_gap, coupled_model_gap_drop, coupled_model_gap_trust, &
                  coupled_model_gap_ratio, coupled_norm_ratio, &
                  coupled_diag_shift_avg, coupled_rot_diag_avg, coupled_occ_response_diag_avg, &
                  coupled_tether_response_diag_avg
                  IF (constrained_reject_reason /= 0 .OR. &
                      constrained_model_drop > EPSILON(constrained_model_drop)) &
                     WRITE (cp_logger_get_default_unit_nr(logger), '(A,3(1X,I6),3(1X,L1),25(1X,ES16.8))') &
                     " OT Schur constrained audit group code samples gate new energy scale drop cdrop "// &
                     "resp rtrust rratio wresp wratio atrust aratio mgap mdrop mtrust mgratio norm davg rdiag "// &
                     "odiag tdiag pdiag "// &
                     "calib qcal split_exp constr_exp exp_ratio:", &
                     target_group, constrained_reject_reason, local_ot_env(1)%schur_model_quality_samples, &
                     split_candidate_enabled, new_search_point, local_ot_env(1)%energy_only, &
                     constrained_direction_scale, constrained_model_drop, constrained_model_calibrated_drop, &
                     constrained_response_step, response_trust, constrained_response_trust_ratio, &
                     constrained_response_weighted_step, constrained_response_weighted_trust_ratio, response_abs_trust, &
                     constrained_response_abs_trust_ratio, constrained_model_gap, constrained_model_gap_drop, &
                     constrained_model_gap_trust, constrained_model_gap_ratio, constrained_norm_ratio, &
                     constrained_diag_shift_avg, &
                     constrained_rot_diag_avg, constrained_occ_response_diag_avg, constrained_tether_response_diag_avg, &
                     constrained_response_penalty_diag_avg, split_curvature_scale, split_quality_for_calibration, &
                     split_expected_drop, constrained_expected_drop, constrained_expected_drop_ratio
                  WRITE (cp_logger_get_default_unit_nr(logger), '(A,2(1X,L1),6(1X,ES16.8))') &
                     " OT Schur exact energy response ready physical gauge diag occ_gap occ_gap_max phys_gap phys_gap_max:", &
                     exact_energy_model_ready, exact_physical_response_ready, energy_gauge_weight, &
                     exact_response_diag_avg, projector_response_gap_avg, projector_response_gap_max, &
                     physical_response_gap_avg, physical_response_gap_max
                  WRITE (cp_logger_get_default_unit_nr(logger), '(A,1X,ES16.8)') &
                     " OT Schur mixed response diag:", response_mixed_diag_avg
                  WRITE (cp_logger_get_default_unit_nr(logger), '(A,4(1X,ES16.8))') &
                     " OT Schur projector Hessian diag delta_diag abs_avg delta_max:", &
                     projector_hessian_diag_avg, projector_hessian_delta_diag_avg, &
                     projector_hessian_abs_avg, projector_hessian_delta_max
                  WRITE (cp_logger_get_default_unit_nr(logger), '(A,6(1X,ES16.8))') &
                     " OT Schur projector physical Hessian diag delta_diag abs_avg delta_max gauge_avg gauge_max:", &
                     projector_physical_diag_avg, projector_physical_delta_diag_avg, &
                     projector_physical_abs_avg, projector_physical_delta_max, &
                     projector_gauge_abs_avg, projector_gauge_delta_max
                  WRITE (cp_logger_get_default_unit_nr(logger), '(A,4(1X,ES16.8))') &
                     " OT Schur occ-gap shadow cdrop odrop ratio gap:", &
                     accepted_model_calibrated_drop, accepted_occgap_shadow_drop, &
                     1.0_dp/MAX(1.0_dp, accepted_occgap_shadow_gap_ratio), &
                     accepted_occgap_shadow_gap_ratio
                  WRITE (cp_logger_get_default_unit_nr(logger), '(A,5(1X,ES16.8))') &
                     " OT Schur projector curvature residual tail eligible dense accepted:", &
                     residual_delta, ot_kp_schur_projector_coupled_tail_delta, &
                     MERGE(1.0_dp, 0.0_dp, residual_delta > ot_kp_schur_projector_coupled_tail_delta), &
                     dense_projector_coupled_weight, accepted_projector_coupled_weight
                  WRITE (cp_logger_get_default_unit_nr(logger), '(A,10(1X,ES16.8))') &
                     " OT Schur projector modes q neg damp emin emax step rneg rdrop rratio trust:", &
                     projector_mode_quality, projector_mode_negative_fraction, projector_mode_damped_abs_avg, &
                     projector_mode_eval_min, projector_mode_eval_max, projector_mode_step_negative_fraction, &
                     projector_mode_response_negative_fraction, projector_mode_response_drop, &
                     projector_mode_response_drop_ratio, &
                     projector_mode_trust_quality
                  WRITE (cp_logger_get_default_unit_nr(logger), '(A,11(1X,ES16.8))') &
                     " OT Schur projector response shadow w drop sneg rneg pdrop ndrop xdrop cpos cneg hratio hrneg:", &
                     projector_mode_response_weight, projector_mode_response_correction_drop, &
                     projector_mode_response_correction_negative, projector_mode_response_raw_negative_fraction, &
                     projector_mode_response_positive_drop, projector_mode_response_negative_drop, &
                     projector_mode_response_cross_drop, projector_mode_response_positive_coeff, &
                     projector_mode_response_negative_coeff, &
                     projector_mode_history_response_ratio, projector_mode_history_response_negative
                     WRITE (cp_logger_get_default_unit_nr(logger), '(A,2(1X,L1),7(1X,ES16.8))') &
                        " OT Schur projector response candidate ready selected exp qsel cdrop mgap rratio aratio weight:", &
                        projector_response_candidate_ready, projector_response_candidate_selected, &
                        projector_response_expected_drop, projector_response_quality_for_selection, &
                        projector_response_model_calibrated_drop, projector_response_model_gap, &
                        projector_response_response_weighted_trust_ratio, &
                        projector_response_response_abs_trust_ratio, projector_mode_response_weight
                     WRITE (cp_logger_get_default_unit_nr(logger), '(A,11(1X,ES16.8))') &
                        " OT Schur projector physical shadow w drop sneg rneg pdrop ndrop xdrop cpos cneg qsel negema:", &
                        projector_physical_response_weight, projector_physical_response_correction_drop, &
                        projector_physical_response_correction_negative, projector_physical_response_raw_negative_fraction, &
                        projector_physical_response_positive_drop, projector_physical_response_negative_drop, &
                        projector_physical_response_cross_drop, projector_physical_response_positive_coeff, &
                        projector_physical_response_negative_coeff, projector_physical_quality_for_selection, &
                        local_ot_env(1)%schur_projector_physical_negative_ema
                     WRITE (cp_logger_get_default_unit_nr(logger), '(A,2(1X,L1),8(1X,ES16.8))') &
                        " OT Schur projector phys direct ready selected exp qsel cdrop mgap rratio aratio comp dshift:", &
                        projector_physical_direct_candidate_ready, projector_physical_direct_candidate_selected, &
                        projector_physical_direct_expected_drop, projector_physical_quality_for_selection, &
                        projector_physical_direct_model_calibrated_drop, projector_physical_direct_model_gap, &
                        projector_physical_direct_response_weighted_trust_ratio, &
                        projector_physical_direct_response_abs_trust_ratio, &
                        projector_physical_direct_diag_compensation, &
                        projector_physical_direct_diag_shift_avg
                     WRITE (cp_logger_get_default_unit_nr(logger), '(A,8(1X,ES16.8))') &
                        " OT Schur projector phys direct model scale drop gapdrop gdotp normp nratio tq sq:", &
                        projector_physical_direct_direction_scale, projector_physical_direct_model_drop, &
                        projector_physical_direct_model_gap_drop, projector_physical_direct_gdotp, &
                        projector_physical_direct_normp, projector_physical_direct_norm_ratio, &
                        projector_physical_direct_trust_quality, projector_physical_direct_scale_quality
                     WRITE (cp_logger_get_default_unit_nr(logger), '(A,2(1X,L1),7(1X,ES16.8))') &
                        " OT Schur projector physical candidate ready selected exp qsel cdrop mgap rratio aratio weight:", &
                        projector_physical_candidate_ready, projector_physical_candidate_selected, &
                        projector_physical_expected_drop, projector_physical_quality_for_selection, &
                        projector_physical_model_calibrated_drop, projector_physical_model_gap, &
                        projector_physical_response_weighted_trust_ratio, &
                        projector_physical_response_abs_trust_ratio, projector_physical_response_weight
                     WRITE (cp_logger_get_default_unit_nr(logger), '(A,3(1X,I6),3(1X,L1),17(1X,ES16.8))') &
                     " OT Schur candidate selection used split_used dense_code split_ready constr_selected "// &
                     "split_selected split_exp dense_exp dense_q drop cdrop split_qsel constr_qsel "// &
                     "split_tq dense_tq constr_tq exp_ratio dense_mgap dense_mdrop dense_mtrust dense_mratio "// &
                     "split_rq dense_rq:", &
                     schur_used_attempt, split_selected_attempt, dense_reject_reason, &
                     split_candidate_ready, constrained_candidate_selected, split_candidate_selected, &
                     split_expected_drop, dense_expected_drop, &
                     dense_quality_for_selection, accepted_model_drop, accepted_model_calibrated_drop, &
                     split_quality_for_selection, constrained_quality_for_selection, split_trust_quality, &
                     dense_trust_quality, constrained_trust_quality, split_dense_expected_ratio, &
                     dense_model_gap, dense_model_gap_drop, dense_model_gap_trust, dense_model_gap_ratio, &
                     split_residual_quality, dense_residual_quality
                  WRITE (cp_logger_get_default_unit_nr(logger), '(A,4(1X,I6),16(1X,ES16.8))') &
                     " OT Schur scale audit group split_code constr_code used floor split_scale split_sq "// &
                     "split_ge split_gr split_ne split_nr constr_scale constr_sq constr_ge constr_gr constr_ne "// &
                     "constr_nr dense_scale dense_sq dense_q:", &
                     target_group, coupled_reject_reason, constrained_reject_reason, schur_used_attempt, &
                     schur_min_coupled_scale, coupled_direction_scale, coupled_scale_quality, &
                     coupled_gdot_energy, coupled_gdot_rotation, coupled_norm_energy, coupled_norm_rotation, &
                     constrained_direction_scale, constrained_scale_quality, constrained_gdot_energy, &
                     constrained_gdot_rotation, constrained_norm_energy, constrained_norm_rotation, &
                     dense_direction_scale, dense_scale_quality, dense_quality_for_selection
               END IF
            END IF
         END IF

            DEALLOCATE (curv, grad_e, occ_curv, tether_curv, z_e, z_exact_physical_trial, z_exact_raw_trial, &
                        z_exact_trial, z_trial, rx_e, rx_model, constrained_z_e, &
                        projector_physical_direct_model_e, projector_physical_direct_rx_e, &
                        projector_physical_direct_z_e, projector_physical_model_e, projector_physical_rx_e, &
                        projector_physical_z_e, &
                        projector_response_rx_e, projector_response_z_e, split_z_e)
            DEALLOCATE (grad_rot, rhs, rtd, sol, constrained_sol, projector_physical_direct_sol, &
                        projector_physical_sol, projector_response_sol, split_sol)
         DEALLOCATE (var_channel, var_ei, var_ej, var_i, var_j, var_kind)
         DEALLOCATE (rot_floor, energy_cross, energy_cross_physical, energy_cross_raw, energy_hessian, &
                     response, response_coupled, response_exact, response_exact_physical, response_exact_raw, &
                     response_occ, response_occ_physical, response_split, rot_hessian, &
                     rot_hessian_coupled, rot_hessian_projector, rot_hessian_projector_physical, schur, &
                     ipiv, ipiv_energy, channel_nmo, e_offset)
      END SUBROUTINE precondition_kpoint_ot_schur_group

      SUBROUTINE project_kpoint_ot_schur_response(response_vec, occ_curv, tether_curv, &
                                                  occ_curv_sum, channel_nmo, e_offset)
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)        :: response_vec
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)           :: occ_curv, tether_curv
      REAL(KIND=dp), INTENT(IN)                         :: occ_curv_sum
      INTEGER, DIMENSION(:), INTENT(IN)                 :: channel_nmo, e_offset

      INTEGER                                           :: first_energy, ie, local_channel, nmo, pass
      REAL(KIND=dp)                                     :: response_shift, weight_sum

         ! The fixed-N Mermin projector is global over the spin/k group, while
         ! the canonical tether gauge is local to each OT channel.
         DO pass = 1, 2
            IF (occ_curv_sum > EPSILON(occ_curv_sum)) THEN
               response_shift = DOT_PRODUCT(occ_curv, response_vec)
               CALL para_env_inter_kp%sum(response_shift)
               response_vec(:) = response_vec(:) - response_shift/occ_curv_sum
            END IF
            DO local_channel = 1, SIZE(channel_nmo)
               nmo = channel_nmo(local_channel)
               IF (nmo <= 0) CYCLE
               first_energy = e_offset(local_channel)
               weight_sum = SUM(tether_curv(first_energy:first_energy + nmo - 1))
               IF (weight_sum <= EPSILON(weight_sum)) CYCLE
               response_shift = DOT_PRODUCT(tether_curv(first_energy:first_energy + nmo - 1), &
                                            response_vec(first_energy:first_energy + nmo - 1))
               response_vec(first_energy:first_energy + nmo - 1) = &
                  response_vec(first_energy:first_energy + nmo - 1) - response_shift/weight_sum
            END DO
         END DO

         DO ie = 1, SIZE(response_vec)
            IF (occ_curv(ie) + tether_curv(ie) <= EPSILON(occ_curv(ie) + tether_curv(ie))) &
               response_vec(ie) = 0.0_dp
         END DO
      END SUBROUTINE project_kpoint_ot_schur_response

      SUBROUTINE kpoint_ot_mixed_projected_term(response_left, response_right, occ_curv, tether_curv, &
                                                occ_curv_sum, channel_nmo, e_offset, projected_term, &
                                                occ_projected_term, tether_projected_term)
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)           :: response_left, response_right, occ_curv, tether_curv
      REAL(KIND=dp), INTENT(IN)                         :: occ_curv_sum
      INTEGER, DIMENSION(:), INTENT(IN)                 :: channel_nmo, e_offset
      REAL(KIND=dp), INTENT(OUT)                        :: projected_term
      REAL(KIND=dp), INTENT(OUT), OPTIONAL              :: occ_projected_term, tether_projected_term

      INTEGER                                           :: first_energy, local_channel, nmo
      REAL(KIND=dp)                                     :: left_shift, occ_term, right_shift, tether_term, &
                                                            weight_sum

         occ_term = DOT_PRODUCT(occ_curv(:)*response_left(:), response_right(:))
         IF (occ_curv_sum > EPSILON(occ_curv_sum)) THEN
            left_shift = DOT_PRODUCT(occ_curv, response_left)
            right_shift = DOT_PRODUCT(occ_curv, response_right)
            occ_term = occ_term - left_shift*right_shift/occ_curv_sum
         END IF
         tether_term = DOT_PRODUCT(tether_curv(:)*response_left(:), response_right(:))
         DO local_channel = 1, SIZE(channel_nmo)
            nmo = channel_nmo(local_channel)
            IF (nmo <= 0) CYCLE
            first_energy = e_offset(local_channel)
            weight_sum = SUM(tether_curv(first_energy:first_energy + nmo - 1))
            IF (weight_sum <= EPSILON(weight_sum)) CYCLE
            left_shift = DOT_PRODUCT(tether_curv(first_energy:first_energy + nmo - 1), &
                                     response_left(first_energy:first_energy + nmo - 1))
            right_shift = DOT_PRODUCT(tether_curv(first_energy:first_energy + nmo - 1), &
                                      response_right(first_energy:first_energy + nmo - 1))
            tether_term = tether_term - left_shift*right_shift/weight_sum
         END DO
         projected_term = occ_term + tether_term
         IF (PRESENT(occ_projected_term)) occ_projected_term = occ_term
         IF (PRESENT(tether_projected_term)) tether_projected_term = tether_term
      END SUBROUTINE kpoint_ot_mixed_projected_term

      SUBROUTINE kpoint_ot_mixed_projected_cross_term(response_occ_left, response_occ_right, &
                                                      response_tether_left, response_tether_right, &
                                                      occ_curv, tether_curv, occ_curv_sum, &
                                                      channel_nmo, e_offset, projected_term, &
                                                      occ_projected_term, tether_projected_term)
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)           :: response_occ_left, response_occ_right, &
                                                            response_tether_left, response_tether_right, &
                                                            occ_curv, tether_curv
      REAL(KIND=dp), INTENT(IN)                         :: occ_curv_sum
      INTEGER, DIMENSION(:), INTENT(IN)                 :: channel_nmo, e_offset
      REAL(KIND=dp), INTENT(OUT)                        :: projected_term
      REAL(KIND=dp), INTENT(OUT), OPTIONAL              :: occ_projected_term, tether_projected_term

      INTEGER                                           :: first_energy, local_channel, nmo
      REAL(KIND=dp)                                     :: left_shift, occ_term, right_shift, tether_term, &
                                                            weight_sum

         occ_term = DOT_PRODUCT(occ_curv(:)*response_occ_left(:), response_occ_right(:))
         IF (occ_curv_sum > EPSILON(occ_curv_sum)) THEN
            left_shift = DOT_PRODUCT(occ_curv, response_occ_left)
            right_shift = DOT_PRODUCT(occ_curv, response_occ_right)
            occ_term = occ_term - left_shift*right_shift/occ_curv_sum
         END IF
         tether_term = DOT_PRODUCT(tether_curv(:)*response_tether_left(:), response_tether_right(:))
         DO local_channel = 1, SIZE(channel_nmo)
            nmo = channel_nmo(local_channel)
            IF (nmo <= 0) CYCLE
            first_energy = e_offset(local_channel)
            weight_sum = SUM(tether_curv(first_energy:first_energy + nmo - 1))
            IF (weight_sum <= EPSILON(weight_sum)) CYCLE
            left_shift = DOT_PRODUCT(tether_curv(first_energy:first_energy + nmo - 1), &
                                     response_tether_left(first_energy:first_energy + nmo - 1))
            right_shift = DOT_PRODUCT(tether_curv(first_energy:first_energy + nmo - 1), &
                                      response_tether_right(first_energy:first_energy + nmo - 1))
            tether_term = tether_term - left_shift*right_shift/weight_sum
         END DO
         projected_term = occ_term + tether_term
         IF (PRESENT(occ_projected_term)) occ_projected_term = occ_term
         IF (PRESENT(tether_projected_term)) tether_projected_term = tether_term
      END SUBROUTINE kpoint_ot_mixed_projected_cross_term

      SUBROUTINE build_kpoint_ot_exact_energy_response(occ_curv, tether_curv, response_occ, &
                                                       response_tether, grad_e, channel_nmo, e_offset, &
                                                       energy_hessian, energy_cross, response_exact, z_exact, &
                                                       ipiv_energy, ready, gauge_weight)
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)           :: occ_curv, tether_curv, grad_e
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)        :: response_occ, response_tether
      INTEGER, DIMENSION(:), INTENT(IN)                 :: channel_nmo, e_offset
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)     :: energy_hessian, energy_cross, response_exact
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)          :: z_exact
      INTEGER, DIMENSION(:), INTENT(INOUT)              :: ipiv_energy
      LOGICAL, INTENT(OUT)                              :: ready
      REAL(KIND=dp), INTENT(OUT)                        :: gauge_weight

      INTEGER                                           :: first_energy, ie, info, je, local_channel, &
                                                           ne_local, nmo, nrhs, nrot_local, r
      REAL(KIND=dp)                                     :: channel_occ_sum, curv_max, occ_sum, &
                                                           response_shift, weight_sum
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)       :: rhs

         ready = .FALSE.
         gauge_weight = 0.0_dp
         ne_local = SIZE(grad_e)
         nrot_local = SIZE(response_occ, 2)
         IF (ne_local <= 0 .OR. nrot_local <= 0) RETURN

         energy_hessian(:, :) = 0.0_dp
         energy_cross(:, :) = 0.0_dp
         response_exact(:, :) = 0.0_dp
         z_exact(:) = 0.0_dp

         DO ie = 1, ne_local
            energy_hessian(ie, ie) = occ_curv(ie) + tether_curv(ie)
         END DO

         occ_sum = SUM(occ_curv)
         CALL para_env_inter_kp%sum(occ_sum)
         IF (occ_sum > EPSILON(occ_sum)) THEN
            DO je = 1, ne_local
               DO ie = 1, ne_local
                  energy_hessian(ie, je) = energy_hessian(ie, je) - &
                                           occ_curv(ie)*occ_curv(je)/occ_sum
               END DO
            END DO
            DO r = 1, nrot_local
               response_shift = DOT_PRODUCT(occ_curv, response_occ(:, r))
               CALL para_env_inter_kp%sum(response_shift)
               energy_cross(:, r) = energy_cross(:, r) + &
                                    occ_curv(:)*(response_occ(:, r) - response_shift/occ_sum)
            END DO
         END IF

         DO local_channel = 1, SIZE(channel_nmo)
            nmo = channel_nmo(local_channel)
            IF (nmo <= 0) CYCLE
            first_energy = e_offset(local_channel)
            weight_sum = SUM(tether_curv(first_energy:first_energy + nmo - 1))
            IF (weight_sum <= EPSILON(weight_sum)) CYCLE
            DO je = first_energy, first_energy + nmo - 1
               DO ie = first_energy, first_energy + nmo - 1
                  energy_hessian(ie, je) = energy_hessian(ie, je) - &
                                           tether_curv(ie)*tether_curv(je)/weight_sum
               END DO
            END DO
            DO r = 1, nrot_local
               response_shift = DOT_PRODUCT(tether_curv(first_energy:first_energy + nmo - 1), &
                                            response_tether(first_energy:first_energy + nmo - 1, r))
               energy_cross(first_energy:first_energy + nmo - 1, r) = &
                  energy_cross(first_energy:first_energy + nmo - 1, r) + &
                  tether_curv(first_energy:first_energy + nmo - 1)* &
                  (response_tether(first_energy:first_energy + nmo - 1, r) - response_shift/weight_sum)
            END DO
         END DO

         curv_max = MAXVAL(occ_curv(:) + tether_curv(:))
         CALL para_env_inter_kp%max(curv_max)
         IF (curv_max <= 0.0_dp) RETURN
         gauge_weight = curv_max/REAL(MAX(1, ne_local), KIND=dp)
         IF (occ_sum > EPSILON(occ_sum)) THEN
            ! Lift the global fixed-N gauge.  The response RHS is fixed-N
            ! projected, so this rank-one term only pins the null direction.
            energy_hessian(:, :) = energy_hessian(:, :) + gauge_weight
         END IF
         DO local_channel = 1, SIZE(channel_nmo)
            nmo = channel_nmo(local_channel)
            IF (nmo <= 0) CYCLE
            first_energy = e_offset(local_channel)
            weight_sum = SUM(tether_curv(first_energy:first_energy + nmo - 1))
            channel_occ_sum = SUM(occ_curv(first_energy:first_energy + nmo - 1))
            IF (weight_sum <= EPSILON(weight_sum) .OR. channel_occ_sum > EPSILON(channel_occ_sum)) CYCLE
            ! If the smearing curvature is absent for a channel, the canonical
            ! tether projector leaves a channel-local constant-energy gauge.
            DO je = first_energy, first_energy + nmo - 1
               DO ie = first_energy, first_energy + nmo - 1
                  energy_hessian(ie, je) = energy_hessian(ie, je) + gauge_weight
               END DO
            END DO
         END DO
         DO ie = 1, ne_local
            IF (occ_curv(ie) + tether_curv(ie) > EPSILON(curv_max)) CYCLE
            energy_hessian(ie, :) = 0.0_dp
            energy_hessian(:, ie) = 0.0_dp
            energy_hessian(ie, ie) = 1.0_dp
            energy_cross(ie, :) = 0.0_dp
         END DO

         nrhs = nrot_local + 1
         ALLOCATE (rhs(ne_local, nrhs))
         rhs(:, 1:nrot_local) = energy_cross(:, 1:nrot_local)
         rhs(:, nrhs) = grad_e(:)
         DO ie = 1, ne_local
            IF (occ_curv(ie) + tether_curv(ie) <= EPSILON(curv_max)) rhs(ie, nrhs) = 0.0_dp
         END DO
         CALL dgesv(ne_local, nrhs, energy_hessian, ne_local, ipiv_energy, rhs, ne_local, info)
         IF (info == 0) THEN
            response_exact(:, 1:nrot_local) = rhs(:, 1:nrot_local)
            z_exact(:) = rhs(:, nrhs)
            ready = .TRUE.
         END IF
         DEALLOCATE (rhs)

      END SUBROUTINE build_kpoint_ot_exact_energy_response

      FUNCTION kpoint_ot_weighted_response_step(response_vec, curv, curv_sum) RESULT(step)
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)           :: response_vec, curv
      REAL(KIND=dp), INTENT(IN)                         :: curv_sum
      REAL(KIND=dp)                                     :: step

      REAL(KIND=dp)                                     :: weighted_norm

         weighted_norm = DOT_PRODUCT(curv(:)*response_vec(:), response_vec(:))
         CALL para_env_inter_kp%sum(weighted_norm)
         IF (curv_sum > EPSILON(curv_sum)) THEN
            step = SQRT(MAX(0.0_dp, weighted_norm/curv_sum))
         ELSE
            step = MAXVAL(ABS(response_vec))
            CALL para_env_inter_kp%max(step)
         END IF

      END FUNCTION kpoint_ot_weighted_response_step

      FUNCTION kpoint_ot_physical_rotation_floor(sigma, row_diag, physical_diag, frontier_floor) RESULT(row_floor)
      REAL(KIND=dp), INTENT(IN)                         :: sigma, row_diag, physical_diag, frontier_floor
      REAL(KIND=dp)                                     :: row_floor

         row_floor = SQRT(EPSILON(1.0_dp))* &
                     MAX(1.0_dp, ABS(row_diag), ABS(physical_diag))
         IF (sigma > EPSILON(sigma)) &
            row_floor = MAX(row_floor, &
                            ot_kp_schur_response_trust_min_factor*kpoint_ot_fermi_window_anchor(sigma))
         row_floor = MIN(frontier_floor, row_floor)

      END FUNCTION kpoint_ot_physical_rotation_floor

      FUNCTION kpoint_ot_schur_scale_quality(direction_scale, scale_floor) RESULT(quality)
      REAL(KIND=dp), INTENT(IN)                         :: direction_scale, scale_floor
      REAL(KIND=dp)                                     :: quality

         IF (scale_floor > EPSILON(scale_floor)) THEN
            quality = MIN(1.0_dp, MAX(0.0_dp, direction_scale)/scale_floor)
         ELSE
            quality = 1.0_dp
         END IF

      END FUNCTION kpoint_ot_schur_scale_quality

      FUNCTION kpoint_ot_schur_history_quality(predicted_sum, measured_sum, ema_quality, &
                                               samples, min_samples, default_quality, &
                                               min_quality, max_quality) RESULT(quality)
      REAL(KIND=dp), INTENT(IN)                         :: predicted_sum, measured_sum, ema_quality, &
                                                            default_quality, min_quality, max_quality
      INTEGER, INTENT(IN)                               :: samples, min_samples
      REAL(KIND=dp)                                     :: quality

         quality = default_quality
         IF (samples >= min_samples) THEN
            IF (predicted_sum > EPSILON(predicted_sum) .AND. measured_sum > 0.0_dp) THEN
               quality = measured_sum/predicted_sum
            ELSE IF (ema_quality > 0.0_dp) THEN
               quality = ema_quality
            END IF
         END IF
         quality = MAX(min_quality, MIN(max_quality, quality))

      END FUNCTION kpoint_ot_schur_history_quality

      FUNCTION kpoint_ot_schur_trust_quality(weighted_trust_ratio, abs_trust_ratio, &
                                             gap_trust_ratio) RESULT(quality)
      REAL(KIND=dp), INTENT(IN)                         :: weighted_trust_ratio, abs_trust_ratio, gap_trust_ratio
      REAL(KIND=dp)                                     :: quality, trust_ratio

         trust_ratio = MAX(1.0_dp, weighted_trust_ratio, abs_trust_ratio, gap_trust_ratio)
         quality = 1.0_dp/trust_ratio

      END FUNCTION kpoint_ot_schur_trust_quality

      FUNCTION kpoint_ot_schur_expected_drop(model_drop, history_quality, trust_quality) RESULT(expected_drop)
      REAL(KIND=dp), INTENT(IN)                         :: model_drop, history_quality, trust_quality
      REAL(KIND=dp)                                     :: expected_drop

         expected_drop = MAX(0.0_dp, model_drop)*MAX(0.0_dp, history_quality)* &
                         MAX(0.0_dp, MIN(1.0_dp, trust_quality))

      END FUNCTION kpoint_ot_schur_expected_drop

!> \brief accepted-drop calibrated finite Projector-Hessian weight for the low residual tail
!> \param residual_delta current OT residual
!> \param history_quality measured/predicted dense-coupled Mermin drop quality
!> \param samples number of accepted dense-coupled quality samples
!> \param density_projector_fraction active canonical density-projector fraction
!> \param projector_diag_avg average finite Projector-Hessian diagonal
!> \param projector_delta_diag_avg average finite Projector-Hessian diagonal shift
!> \return full finite Projector weight above the validated tail gate, otherwise
!>         a small Mermin-mismatch-calibrated correction
      FUNCTION kpoint_ot_projector_coupled_weight(residual_delta, history_quality, samples, &
                                                  density_projector_fraction, projector_diag_avg, &
                                                  projector_delta_diag_avg) RESULT(projector_weight)
      REAL(KIND=dp), INTENT(IN)                         :: residual_delta, history_quality, &
                                                            density_projector_fraction, projector_diag_avg, &
                                                            projector_delta_diag_avg
      INTEGER, INTENT(IN)                               :: samples
      REAL(KIND=dp)                                     :: diag_ratio, projector_weight, quality_weight

         IF (residual_delta > ot_kp_schur_projector_coupled_tail_delta) THEN
            projector_weight = 1.0_dp
            RETURN
         END IF

         projector_weight = 0.0_dp
         IF (samples < ot_kp_schur_calibration_min_samples) RETURN

         quality_weight = (1.0_dp - MIN(1.0_dp, MAX(0.0_dp, history_quality)))/ &
                          MAX(EPSILON(1.0_dp), 1.0_dp - ot_kp_schur_calibrated_quality_min_ratio)
         quality_weight = MAX(0.0_dp, MIN(1.0_dp, quality_weight))
         diag_ratio = ABS(projector_delta_diag_avg)/MAX(ABS(projector_diag_avg), EPSILON(1.0_dp))
         projector_weight = ot_kp_schur_projector_coupled_tail_max_weight* &
                            quality_weight*MIN(1.0_dp, MAX(0.0_dp, density_projector_fraction))/ &
                            MAX(1.0_dp, 1.0_dp + diag_ratio)

      END FUNCTION kpoint_ot_projector_coupled_weight

!> \brief audit an accepted-history calibrated finite-response correction in Projector modes
!> \param rot_hessian diagonal/canonical rotation Hessian
!> \param rot_hessian_projector finite complex Projector-Hessian
!> \param step_rot Schur rotation step before line-search normalization
!> \param response_raw conservative response used by the accepted update
!> \param response_exact finite fixed-N/tether response model
!> \param curv energy curvature weights
!> \param projector_weight scalar finite Projector-Hessian tail weight
!> \param residual_delta current OT residual
!> \param history_quality accepted Mermin-drop quality of Projector-mode predictions
!> \param history_trust accepted Projector-mode trust EMA
!> \param history_response_ratio accepted response-mismatch drop over predicted drop
!> \param grad_e Mermin energy-variable gradient
!> \param energy_base energy search vector before finite-response correction
!> \param response_base conservative response vector used by the accepted update
!> \param response_candidate finite-response candidate vector
!> \param applied_weight step-weighted average shadow correction weight
!> \param correction_drop curvature-weighted norm of the applied correction
!> \param correction_negative_fraction fraction of sign-aware correction in negative Projector modes
!> \param raw_negative_fraction raw negative finite Projector-response mode fraction
!> \param positive_drop positive-mode contribution to the signed two-mode model
!> \param negative_drop negative-mode contribution to the signed two-mode model
!> \param cross_drop positive/negative cross contribution to the signed two-mode model
!> \param positive_coeff bounded positive-subspace response coefficient
!> \param negative_coeff bounded negative-subspace response coefficient
      SUBROUTINE audit_kpoint_ot_projector_response_modes(rot_hessian, rot_hessian_projector, step_rot, &
                                                          response_raw, response_exact, curv, projector_weight, &
                                                          residual_delta, history_quality, history_trust, &
                                                          history_response_ratio, grad_e, energy_base, response_base, &
                                                          response_candidate, applied_weight, &
                                                          correction_drop, correction_negative_fraction, &
                                                          raw_negative_fraction, positive_drop, negative_drop, cross_drop, &
                                                          positive_coeff, negative_coeff)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: rot_hessian, rot_hessian_projector, &
                                                            response_raw, response_exact
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: curv, energy_base, grad_e, step_rot
      REAL(KIND=dp), INTENT(IN)                          :: projector_weight, residual_delta, history_quality, &
                                                            history_trust, history_response_ratio
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: response_base
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: response_candidate
      REAL(KIND=dp), INTENT(OUT)                         :: applied_weight, correction_drop, &
                                                            correction_negative_fraction, cross_drop, negative_drop, &
                                                            negative_coeff, positive_coeff, positive_drop, &
                                                            raw_negative_fraction

      INTEGER                                            :: imode, nrot_local
      REAL(KIND=dp)                                      :: base_weight, coeff_det, coeff_negative, coeff_positive, &
                                                            eval_abs_avg, h12, h11, h22, g1, g2, mode_norm, mode_sign, &
                                                            mode_weight, negative_balance, negative_norm, &
                                                            negative_weight_sum, positive_sum, positive_weight_sum, &
                                                            raw_negative_sum, response_total_sum, residual_weight, &
                                                            spectral_weight, step_weight, step_weight_sum
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: mode_coeff, mode_eval, mode_norms, mode_rx, &
                                                            response_corr, response_neg, response_pos
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: mode_vec

         applied_weight = 0.0_dp
         correction_drop = 0.0_dp
         correction_negative_fraction = 0.0_dp
         cross_drop = 0.0_dp
         negative_coeff = 0.0_dp
         negative_drop = 0.0_dp
         positive_coeff = 0.0_dp
         positive_drop = 0.0_dp
         raw_negative_fraction = 0.0_dp
         response_candidate(:) = response_base(:)
         IF (projector_weight <= EPSILON(projector_weight)) RETURN
         IF (history_quality < ot_kp_schur_calibrated_quality_min_ratio) RETURN
         IF (history_trust <= EPSILON(history_trust)) RETURN

         nrot_local = SIZE(rot_hessian, 1)
         IF (nrot_local <= 0) RETURN
         IF (SIZE(response_base) /= SIZE(curv)) RETURN
         IF (SIZE(energy_base) /= SIZE(curv)) RETURN
         IF (SIZE(grad_e) /= SIZE(curv)) RETURN

         residual_weight = (residual_delta - ot_kp_schur_projector_response_floor_delta)/ &
                           MAX(EPSILON(1.0_dp), &
                               ot_kp_schur_projector_response_full_delta - &
                               ot_kp_schur_projector_response_floor_delta)
         residual_weight = MAX(0.0_dp, MIN(1.0_dp, residual_weight))
         base_weight = ot_kp_schur_projector_response_max_weight*residual_weight* &
                       MIN(1.0_dp, MAX(0.0_dp, projector_weight))* &
                       MIN(1.0_dp, MAX(0.0_dp, history_quality))* &
                       MIN(1.0_dp, MAX(0.0_dp, history_trust))/ &
                       MAX(1.0_dp, history_response_ratio/ &
                           ot_kp_schur_projector_response_ratio_gate)
         IF (base_weight <= EPSILON(base_weight)) RETURN

         ALLOCATE (mode_coeff(nrot_local), mode_eval(nrot_local), mode_norms(nrot_local), &
                   mode_rx(SIZE(curv)), response_corr(SIZE(curv)), response_neg(SIZE(curv)), &
                   response_pos(SIZE(curv)), mode_vec(nrot_local, nrot_local))
         response_corr(:) = 0.0_dp
         response_neg(:) = 0.0_dp
         response_pos(:) = 0.0_dp
         mode_vec(:, :) = 0.5_dp*((rot_hessian_projector(:, :) - rot_hessian(:, :)) + &
                                  TRANSPOSE(rot_hessian_projector(:, :) - rot_hessian(:, :)))
         CALL diamat_all(mode_vec, mode_eval)
         mode_coeff(:) = MATMUL(TRANSPOSE(mode_vec), step_rot)
         eval_abs_avg = SUM(ABS(mode_eval))/REAL(MAX(1, nrot_local), KIND=dp)
         response_total_sum = 0.0_dp
         raw_negative_sum = 0.0_dp
         mode_norms(:) = 0.0_dp
         DO imode = 1, nrot_local
            IF (ABS(mode_coeff(imode)) <= EPSILON(mode_coeff(imode))) CYCLE
            mode_rx(:) = mode_coeff(imode)*MATMUL(response_exact - response_raw, mode_vec(:, imode))
            mode_norm = DOT_PRODUCT(curv(:)*mode_rx(:), mode_rx(:))
            CALL para_env_inter_kp%sum(mode_norm)
            mode_norms(imode) = MAX(0.0_dp, mode_norm)
            response_total_sum = response_total_sum + mode_norms(imode)
            IF (mode_eval(imode) < 0.0_dp) raw_negative_sum = raw_negative_sum + mode_norms(imode)
         END DO
         IF (response_total_sum > EPSILON(response_total_sum)) &
            raw_negative_fraction = raw_negative_sum/response_total_sum
         positive_sum = MAX(0.0_dp, response_total_sum - raw_negative_sum)
         negative_balance = 0.0_dp
         IF (raw_negative_sum > EPSILON(raw_negative_sum)) &
            negative_balance = SQRT(positive_sum/raw_negative_sum)
         negative_balance = MIN(1.0_dp, MAX(0.0_dp, negative_balance))

         step_weight_sum = DOT_PRODUCT(mode_coeff, mode_coeff)
         negative_weight_sum = 0.0_dp
         positive_weight_sum = 0.0_dp
         DO imode = 1, nrot_local
            IF (ABS(mode_coeff(imode)) <= EPSILON(mode_coeff(imode))) CYCLE
            spectral_weight = 1.0_dp
            mode_sign = 1.0_dp
            IF (mode_eval(imode) < 0.0_dp) THEN
               mode_sign = -1.0_dp
               spectral_weight = 1.0_dp/MAX(1.0_dp, ABS(mode_eval(imode))/ &
                                            MAX(eval_abs_avg, EPSILON(eval_abs_avg)))
               spectral_weight = spectral_weight*negative_balance
            END IF
            mode_weight = base_weight*MAX(0.0_dp, MIN(1.0_dp, spectral_weight))
            IF (mode_weight <= EPSILON(mode_weight)) CYCLE
            mode_rx(:) = mode_coeff(imode)*MATMUL(response_exact - response_raw, mode_vec(:, imode))
            step_weight = mode_coeff(imode)*mode_coeff(imode)
            IF (mode_sign > 0.0_dp) THEN
               response_pos(:) = response_pos(:) + mode_weight*mode_rx(:)
               positive_weight_sum = positive_weight_sum + mode_weight*step_weight
            ELSE
               response_neg(:) = response_neg(:) + mode_weight*mode_rx(:)
               negative_weight_sum = negative_weight_sum + mode_weight*step_weight
            END IF
         END DO
         h11 = DOT_PRODUCT(curv(:)*response_pos(:), response_pos(:))
         h12 = -DOT_PRODUCT(curv(:)*response_pos(:), response_neg(:))
         h22 = DOT_PRODUCT(curv(:)*response_neg(:), response_neg(:))
         g1 = DOT_PRODUCT(grad_e(:) - curv(:)*energy_base(:), response_pos(:))
         g2 = -DOT_PRODUCT(grad_e(:) - curv(:)*energy_base(:), response_neg(:))
         CALL para_env_inter_kp%sum(h11)
         CALL para_env_inter_kp%sum(h12)
         CALL para_env_inter_kp%sum(h22)
         CALL para_env_inter_kp%sum(g1)
         CALL para_env_inter_kp%sum(g2)
         coeff_positive = 1.0_dp
         coeff_negative = 1.0_dp
         coeff_det = h11*h22 - h12*h12
         IF (coeff_det > EPSILON(coeff_det)) THEN
            coeff_positive = (g1*h22 - g2*h12)/coeff_det
            coeff_negative = (h11*g2 - h12*g1)/coeff_det
            coeff_positive = MIN(1.0_dp, MAX(0.0_dp, coeff_positive))
            coeff_negative = MIN(1.0_dp, MAX(0.0_dp, coeff_negative))
         END IF
         positive_coeff = coeff_positive
         negative_coeff = coeff_negative
         response_corr(:) = coeff_positive*response_pos(:) - coeff_negative*response_neg(:)
         response_candidate(:) = response_base(:) + response_corr(:)
         positive_drop = 0.5_dp*coeff_positive*coeff_positive*MAX(0.0_dp, h11)
         negative_drop = 0.5_dp*coeff_negative*coeff_negative*MAX(0.0_dp, h22)
         cross_drop = coeff_positive*coeff_negative*h12
         correction_drop = DOT_PRODUCT(curv(:)*response_corr(:), response_corr(:))
         negative_norm = DOT_PRODUCT(curv(:)*(coeff_negative*response_neg(:)), coeff_negative*response_neg(:))
         CALL para_env_inter_kp%sum(correction_drop)
         CALL para_env_inter_kp%sum(negative_norm)
         correction_drop = 0.5_dp*MAX(0.0_dp, correction_drop)
         IF (step_weight_sum > EPSILON(step_weight_sum)) &
            applied_weight = (coeff_positive*positive_weight_sum + coeff_negative*negative_weight_sum)/step_weight_sum
         IF (correction_drop > EPSILON(correction_drop)) &
            correction_negative_fraction = negative_norm/(2.0_dp*correction_drop)
         DEALLOCATE (mode_coeff, mode_eval, mode_norms, mode_rx, response_corr, response_neg, response_pos, mode_vec)

      END SUBROUTINE audit_kpoint_ot_projector_response_modes

!> \brief audit a Schur step in the finite Projector-Hessian eigenbasis
!> \param rot_hessian diagonal/canonical rotation Hessian
!> \param rot_hessian_projector finite complex Projector-Hessian
!> \param step_rot scaled Schur rotation step
!> \param response_actual actual energy response used by the candidate
!> \param response_model calibrated response model used for mismatch
!> \param curv energy curvature weights
!> \param candidate_drop calibrated candidate Mermin drop before mode trust
!> \param projector_weight scalar accepted-drop calibrated projector weight
!> \param history_quality accepted Mermin-drop quality of recent coupled steps
!> \param residual_quality residual-quality calibration of recent coupled steps
!> \param mode_quality accepted/residual quality used for Projector-Hessian diagnostics
!> \param mode_negative_fraction fraction of Projector-Hessian delta eigenmodes with negative curvature
!> \param mode_damped_abs_avg average absolute spectral curvature correction
!> \param mode_eval_min minimum Projector-Hessian delta eigenvalue
!> \param mode_eval_max maximum Projector-Hessian delta eigenvalue
!> \param mode_step_negative_fraction fraction of the rotation step in negative eigenmodes
!> \param mode_response_negative_fraction fraction of response mismatch in negative eigenmodes
!> \param mode_response_drop total mode-decomposed response-mismatch drop
!> \param mode_response_drop_ratio response-mismatch drop relative to candidate drop
!> \param mode_trust_quality conservative quality factor for candidate expected drop
      SUBROUTINE audit_kpoint_ot_projector_modes(rot_hessian, rot_hessian_projector, step_rot, &
                                                 response_actual, response_model, curv, candidate_drop, projector_weight, &
                                                 history_quality, residual_quality, mode_quality, &
                                                 mode_negative_fraction, mode_damped_abs_avg, &
                                                 mode_eval_min, mode_eval_max, mode_step_negative_fraction, &
                                                 mode_response_negative_fraction, mode_response_drop, &
                                                 mode_response_drop_ratio, mode_trust_quality)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: rot_hessian, rot_hessian_projector, &
                                                            response_actual, response_model
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: curv, step_rot
      REAL(KIND=dp), INTENT(IN)                          :: candidate_drop, projector_weight, history_quality, &
                                                            residual_quality
      REAL(KIND=dp), INTENT(OUT)                         :: mode_quality, mode_negative_fraction, &
                                                            mode_damped_abs_avg, mode_eval_min, mode_eval_max, &
                                                            mode_step_negative_fraction, &
                                                            mode_response_negative_fraction, mode_response_drop, &
                                                            mode_response_drop_ratio, mode_trust_quality

      INTEGER                                            :: imode, nrot_local
      REAL(KIND=dp)                                      :: damped_abs_sum, filtered_eval, mode_response_norm, &
                                                            mode_risk, mode_soft_shift, negative_count, &
                                                            response_negative_sum, response_total_sum, &
                                                            step_negative_sum, step_total_sum
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: mode_coeff, mode_eval, mode_rx
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: mode_vec

         mode_quality = 0.0_dp
         mode_negative_fraction = 0.0_dp
         mode_damped_abs_avg = 0.0_dp
         mode_eval_min = 0.0_dp
         mode_eval_max = 0.0_dp
         mode_step_negative_fraction = 0.0_dp
         mode_response_negative_fraction = 0.0_dp
         mode_response_drop = 0.0_dp
         mode_response_drop_ratio = 0.0_dp
         mode_trust_quality = 1.0_dp
         IF (projector_weight <= EPSILON(projector_weight)) RETURN

         nrot_local = SIZE(rot_hessian, 1)
         IF (nrot_local <= 0) RETURN

         ALLOCATE (mode_coeff(nrot_local), mode_eval(nrot_local), mode_rx(SIZE(curv)), &
                   mode_vec(nrot_local, nrot_local))
         mode_vec(:, :) = 0.5_dp*((rot_hessian_projector(:, :) - rot_hessian(:, :)) + &
                                  TRANSPOSE(rot_hessian_projector(:, :) - rot_hessian(:, :)))
         CALL diamat_all(mode_vec, mode_eval)
         mode_eval_min = MINVAL(mode_eval)
         mode_eval_max = MAXVAL(mode_eval)

         mode_quality = MIN(1.0_dp, MAX(0.0_dp, history_quality))* &
                        MIN(1.0_dp, MAX(0.0_dp, residual_quality))
         mode_soft_shift = MAX(0.0_dp, -SUM(mode_eval)/REAL(nrot_local, KIND=dp))
         negative_count = 0.0_dp
         damped_abs_sum = 0.0_dp
         response_negative_sum = 0.0_dp
         response_total_sum = 0.0_dp
         step_negative_sum = 0.0_dp
         mode_coeff(:) = MATMUL(TRANSPOSE(mode_vec), step_rot)
         step_total_sum = DOT_PRODUCT(mode_coeff, mode_coeff)
         DO imode = 1, nrot_local
            filtered_eval = mode_eval(imode) + mode_soft_shift
            IF (filtered_eval < 0.0_dp) THEN
               negative_count = negative_count + 1.0_dp
               step_negative_sum = step_negative_sum + mode_coeff(imode)*mode_coeff(imode)
            END IF
            damped_abs_sum = damped_abs_sum + ABS(mode_eval(imode) - filtered_eval)
            mode_rx(:) = mode_coeff(imode)*MATMUL(response_model - response_actual, mode_vec(:, imode))
            mode_response_norm = DOT_PRODUCT(curv(:)*mode_rx(:), mode_rx(:))
            CALL para_env_inter_kp%sum(mode_response_norm)
            response_total_sum = response_total_sum + mode_response_norm
            IF (filtered_eval < 0.0_dp) response_negative_sum = response_negative_sum + mode_response_norm
         END DO

         mode_negative_fraction = negative_count/REAL(MAX(1, nrot_local), KIND=dp)
         mode_damped_abs_avg = projector_weight*damped_abs_sum/REAL(MAX(1, nrot_local), KIND=dp)
         IF (step_total_sum > EPSILON(step_total_sum)) &
            mode_step_negative_fraction = step_negative_sum/step_total_sum
         IF (response_total_sum > EPSILON(response_total_sum)) &
            mode_response_negative_fraction = response_negative_sum/response_total_sum
         mode_response_drop = 0.5_dp*MAX(0.0_dp, response_total_sum)
         IF (candidate_drop > EPSILON(candidate_drop)) &
            mode_response_drop_ratio = mode_response_drop/candidate_drop
         mode_risk = MAX(mode_response_negative_fraction, &
                         mode_step_negative_fraction*MIN(1.0_dp, mode_response_drop_ratio)) - mode_quality
         mode_trust_quality = 1.0_dp - ot_kp_schur_projector_mode_trust_max_penalty* &
                              projector_weight*MAX(0.0_dp, mode_risk)
         mode_trust_quality = MAX(1.0_dp - ot_kp_schur_projector_mode_trust_max_penalty*projector_weight, &
                                  MIN(1.0_dp, mode_trust_quality))
         DEALLOCATE (mode_coeff, mode_eval, mode_rx, mode_vec)

      END SUBROUTINE audit_kpoint_ot_projector_modes

      FUNCTION kpoint_ot_response_mismatch_drop(response_actual, response_model, curv, scale) RESULT(drop)
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)           :: response_actual, response_model, curv
      REAL(KIND=dp), INTENT(IN)                         :: scale
      REAL(KIND=dp)                                     :: drop

         drop = SUM(curv(:)*(response_model(:) - response_actual(:))**2)
         CALL para_env_inter_kp%sum(drop)
         drop = 0.5_dp*scale*scale*MAX(0.0_dp, drop)

      END FUNCTION kpoint_ot_response_mismatch_drop

      SUBROUTINE build_finite_ref_rotation_model(local_channel, nmo, first_var, last_var, first_energy, &
                                                 chc_re, chc_im, occupation_numbers, kpoint_weight, &
                                                 occ_curv_channel, var_i, var_j, var_kind, compute_projector_shadow, &
                                                 rot_hessian, rot_hessian_coupled, rot_hessian_projector, &
                                                 rot_hessian_projector_physical, response, response_occ, &
                                                 response_occ_physical, failed)
      INTEGER, INTENT(IN)                                :: local_channel, nmo, first_var, last_var, &
                                                            first_energy
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: chc_re, chc_im
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: occupation_numbers, occ_curv_channel
      REAL(KIND=dp), INTENT(IN)                          :: kpoint_weight
      INTEGER, DIMENSION(:), INTENT(IN)                  :: var_i, var_j, var_kind
      LOGICAL, INTENT(IN)                                :: compute_projector_shadow
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: rot_hessian, rot_hessian_coupled, &
                                                            rot_hessian_projector, rot_hessian_projector_physical, &
                                                            response, response_occ, response_occ_physical
      LOGICAL, INTENT(INOUT)                             :: failed

      REAL(KIND=dp), PARAMETER                           :: fd_step = 1.0E-4_dp

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: base, chc, chc_m, chc_p, dchc, docc_mat, &
                                                            occ_basis, occ_mat, occ_mat_coupled, occ_mat_m, occ_mat_p, &
                                                            projector_evec, u
      INTEGER                                            :: ie, imo, jmo, r, s
      REAL(KIND=dp)                                      :: projector_response, projector_weight, raw_response
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eps_m, eps_p, projector_evals, projector_evals_m, &
                                                            projector_evals_p
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: grad_m_im, grad_m_re, grad_p_im, &
                                                            grad_p_re, u_im, u_re, x_im, x_m_im, &
                                                            x_m_re, x_p_im, x_p_re, x_re

         IF (failed) RETURN
         IF (first_var > last_var) RETURN

         ALLOCATE (x_re(nmo, nmo), x_im(nmo, nmo), u_re(nmo, nmo), u_im(nmo, nmo))
         CALL dbcsr_square_to_dense(local_ot_env(local_channel)%rot_mat_x, x_re)
         CALL dbcsr_square_to_dense(local_ot_env(local_channel)%rot_mat_x_im, x_im)
         CALL dbcsr_square_to_dense(local_ot_env(local_channel)%rot_mat_u, u_re)
         CALL dbcsr_square_to_dense(local_ot_env(local_channel)%rot_mat_u_im, u_im)
         CALL local_ot_env(local_channel)%para_env%sum(x_re)
         CALL local_ot_env(local_channel)%para_env%sum(x_im)
         CALL local_ot_env(local_channel)%para_env%sum(u_re)
         CALL local_ot_env(local_channel)%para_env%sum(u_im)

         ALLOCATE (chc(nmo, nmo), u(nmo, nmo), base(nmo, nmo), &
                   chc_p(nmo, nmo), chc_m(nmo, nmo), dchc(nmo, nmo), docc_mat(nmo, nmo), &
                   occ_mat(nmo, nmo), occ_mat_p(nmo, nmo), occ_mat_m(nmo, nmo), &
                   occ_mat_coupled(nmo, nmo), occ_basis(nmo, nmo), projector_evec(nmo, nmo))
         chc = CMPLX(chc_re, chc_im, KIND=dp)
         u = CMPLX(u_re, u_im, KIND=dp)
         base = MATMUL(u, MATMUL(chc, CONJG(TRANSPOSE(u))))
         CALL dense_history_safe_occ_projector(chc, occupation_numbers, kpoint_weight, &
                                               local_ot_env(local_channel)%delta, occ_mat)
         occ_mat_coupled(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
         DO imo = 1, nmo
            occ_mat_coupled(imo, imo) = CMPLX(2.0_dp*kpoint_weight*occupation_numbers(imo), 0.0_dp, KIND=dp)
         END DO
         projector_weight = kpoint_ot_canonical_projector_weight(local_ot_env(local_channel)%delta, &
                                                                 ot_kp_canonical_projector_max_weight)
         IF (projector_weight > 0.0_dp) THEN
            ALLOCATE (projector_evals(nmo))
            CALL diag_complex(chc, projector_evec, projector_evals)
            DEALLOCATE (projector_evals)
         END IF

         ALLOCATE (x_p_re(nmo, nmo), x_p_im(nmo, nmo), x_m_re(nmo, nmo), x_m_im(nmo, nmo))
         ALLOCATE (grad_p_re(nmo, nmo), grad_p_im(nmo, nmo), grad_m_re(nmo, nmo), grad_m_im(nmo, nmo))
         ALLOCATE (eps_p(nmo), eps_m(nmo), projector_evals_p(nmo), projector_evals_m(nmo))

         DO ie = 1, nmo
            CALL dense_history_safe_occ_projector_basis(projector_evec, projector_weight, ie, occ_basis)
            CALL dense_ref_rotation_gradient(x_re, x_im, base, occ_basis, grad_p_re, grad_p_im, eps_p)
            DO r = first_var, last_var
               imo = var_i(r)
               jmo = var_j(r)
               SELECT CASE (var_kind(r))
               CASE (1)
                  response_occ(first_energy + ie - 1, r) = 2.0_dp*grad_p_re(imo, jmo)
               CASE (2)
                  response_occ(first_energy + ie - 1, r) = 2.0_dp*grad_p_im(imo, jmo)
               END SELECT
            END DO
         END DO

         DO s = first_var, last_var
            x_p_re(:, :) = x_re(:, :)
            x_p_im(:, :) = x_im(:, :)
            x_m_re(:, :) = x_re(:, :)
            x_m_im(:, :) = x_im(:, :)
            imo = var_i(s)
            jmo = var_j(s)
            SELECT CASE (var_kind(s))
            CASE (1)
               x_p_re(imo, jmo) = x_p_re(imo, jmo) + fd_step
               x_p_re(jmo, imo) = x_p_re(jmo, imo) - fd_step
               x_m_re(imo, jmo) = x_m_re(imo, jmo) - fd_step
               x_m_re(jmo, imo) = x_m_re(jmo, imo) + fd_step
            CASE (2)
               x_p_im(imo, jmo) = x_p_im(imo, jmo) + fd_step
               x_p_im(jmo, imo) = x_p_im(jmo, imo) + fd_step
               x_m_im(imo, jmo) = x_m_im(imo, jmo) - fd_step
               x_m_im(jmo, imo) = x_m_im(jmo, imo) - fd_step
            END SELECT

            CALL dense_ref_rotation_gradient(x_p_re, x_p_im, base, occ_mat, grad_p_re, grad_p_im, eps_p)
            CALL dense_ref_rotation_gradient(x_m_re, x_m_im, base, occ_mat, grad_m_re, grad_m_im, eps_m)

            DO r = first_var, last_var
               imo = var_i(r)
               jmo = var_j(r)
               SELECT CASE (var_kind(r))
               CASE (1)
                  rot_hessian(r, s) = 2.0_dp*(grad_p_re(imo, jmo) - grad_m_re(imo, jmo))/(2.0_dp*fd_step)
               CASE (2)
                  rot_hessian(r, s) = 2.0_dp*(grad_p_im(imo, jmo) - grad_m_im(imo, jmo))/(2.0_dp*fd_step)
               END SELECT
            END DO
            CALL dense_rotated_chc(x_p_re, x_p_im, base, chc_p)
            CALL dense_rotated_chc(x_m_re, x_m_im, base, chc_m)
            IF (compute_projector_shadow .AND. projector_weight > 0.0_dp) THEN
               CALL diag_complex(chc_p, projector_evec, projector_evals_p)
               CALL diag_complex(chc_m, projector_evec, projector_evals_m)
            END IF
            DO ie = 1, nmo
               raw_response = (REAL(chc_p(ie, ie), KIND=dp) - REAL(chc_m(ie, ie), KIND=dp))/(2.0_dp*fd_step)
               projector_response = raw_response
               IF (compute_projector_shadow .AND. projector_weight > 0.0_dp) &
                  projector_response = (projector_evals_p(ie) - projector_evals_m(ie))/(2.0_dp*fd_step)
               response_occ_physical(first_energy + ie - 1, s) = &
                  (1.0_dp - projector_weight)*raw_response + projector_weight*projector_response
            END DO
            CALL dense_history_safe_occ_projector(chc_p, occupation_numbers, kpoint_weight, &
                                                  local_ot_env(local_channel)%delta, occ_mat_p)
            CALL dense_history_safe_occ_projector(chc_m, occupation_numbers, kpoint_weight, &
                                                  local_ot_env(local_channel)%delta, occ_mat_m)
            CALL dense_ref_rotation_gradient(x_p_re, x_p_im, base, occ_mat_p, grad_p_re, grad_p_im, eps_p)
            CALL dense_ref_rotation_gradient(x_m_re, x_m_im, base, occ_mat_m, grad_m_re, grad_m_im, eps_m)
            DO r = first_var, last_var
               imo = var_i(r)
               jmo = var_j(r)
               SELECT CASE (var_kind(r))
               CASE (1)
                  rot_hessian_projector(r, s) = &
                     2.0_dp*(grad_p_re(imo, jmo) - grad_m_re(imo, jmo))/(2.0_dp*fd_step)
               CASE (2)
                  rot_hessian_projector(r, s) = &
                     2.0_dp*(grad_p_im(imo, jmo) - grad_m_im(imo, jmo))/(2.0_dp*fd_step)
               END SELECT
            END DO
            IF (compute_projector_shadow) THEN
               dchc(:, :) = (chc_p(:, :) - chc_m(:, :))/(2.0_dp*fd_step)
               CALL dense_history_safe_occ_projector_derivative(chc, dchc, occupation_numbers, kpoint_weight, &
                                                                local_ot_env(local_channel)%delta, docc_mat)
               DO ie = 1, nmo
                  IF (ABS(occ_curv_channel(ie)) > curvature_floor) THEN
                     response_occ_physical(first_energy + ie - 1, s) = &
                        response_occ_physical(first_energy + ie - 1, s) - &
                        REAL(docc_mat(ie, ie), KIND=dp)/occ_curv_channel(ie)
                     ! The diagonal projector response is represented by the
                     ! energy-response channel above.  Keeping it in the
                     ! rotation Hessian as well double-counts the same density
                     ! response and creates an artificial physical row floor.
                     docc_mat(ie, ie) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
                  END IF
               END DO
               CALL dense_ref_rotation_gradient(x_re, x_im, base, docc_mat, grad_p_re, grad_p_im, eps_p)
               DO r = first_var, last_var
                  imo = var_i(r)
                  jmo = var_j(r)
                  SELECT CASE (var_kind(r))
                  CASE (1)
                     rot_hessian_projector_physical(r, s) = rot_hessian(r, s) + 2.0_dp*grad_p_re(imo, jmo)
                  CASE (2)
                     rot_hessian_projector_physical(r, s) = rot_hessian(r, s) + 2.0_dp*grad_p_im(imo, jmo)
                  END SELECT
               END DO
            ELSE
               DO r = first_var, last_var
                  rot_hessian_projector_physical(r, s) = rot_hessian_projector(r, s)
               END DO
            END IF
            CALL dense_ref_rotation_gradient(x_p_re, x_p_im, base, occ_mat_coupled, grad_p_re, &
                                             grad_p_im, eps_p)
            CALL dense_ref_rotation_gradient(x_m_re, x_m_im, base, occ_mat_coupled, grad_m_re, &
                                             grad_m_im, eps_m)
            DO r = first_var, last_var
               imo = var_i(r)
               jmo = var_j(r)
               SELECT CASE (var_kind(r))
               CASE (1)
                  rot_hessian_coupled(r, s) = &
                     2.0_dp*(grad_p_re(imo, jmo) - grad_m_re(imo, jmo))/(2.0_dp*fd_step)
               CASE (2)
                  rot_hessian_coupled(r, s) = &
                     2.0_dp*(grad_p_im(imo, jmo) - grad_m_im(imo, jmo))/(2.0_dp*fd_step)
               END SELECT
            END DO
            DO ie = 1, nmo
               raw_response = (eps_p(ie) - eps_m(ie))/(2.0_dp*fd_step)
               response(first_energy + ie - 1, s) = raw_response
            END DO
         END DO

         DEALLOCATE (x_re, x_im, u_re, u_im, chc, u, base, chc_p, chc_m, dchc, docc_mat, &
                     occ_basis, occ_mat, occ_mat_p, occ_mat_m, occ_mat_coupled, projector_evec)
         DEALLOCATE (x_p_re, x_p_im, x_m_re, x_m_im, grad_p_re, grad_p_im, grad_m_re, grad_m_im)
         DEALLOCATE (eps_p, eps_m, projector_evals_p, projector_evals_m)
      END SUBROUTINE build_finite_ref_rotation_model

      SUBROUTINE dense_history_safe_occ_projector_basis(evec, projector_weight, iband, occ_mat)
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: evec
      REAL(KIND=dp), INTENT(IN)                          :: projector_weight
      INTEGER, INTENT(IN)                                :: iband
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(OUT)     :: occ_mat

      INTEGER                                            :: i, j, n

         n = SIZE(evec, 1)
         occ_mat(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)

         IF (projector_weight > 0.0_dp) THEN
            DO j = 1, n
               DO i = 1, n
                  occ_mat(i, j) = projector_weight*evec(i, iband)*CONJG(evec(j, iband))
               END DO
            END DO
         END IF

         occ_mat(iband, iband) = occ_mat(iband, iband) + &
                                  CMPLX(1.0_dp - projector_weight, 0.0_dp, KIND=dp)

      END SUBROUTINE dense_history_safe_occ_projector_basis

      SUBROUTINE dense_rotated_chc(x_re, x_im, base, chc)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: x_re, x_im
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: base
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(OUT)     :: chc

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: evec, generator, u
      INTEGER                                            :: j, n
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: evals

         n = SIZE(x_re, 1)
         ALLOCATE (generator(n, n), evec(n, n), u(n, n), evals(n))

         generator(:, :) = CMPLX(-x_im(:, :), x_re(:, :), KIND=dp)
         CALL diag_complex(generator, evec, evals)

         u(:, :) = evec(:, :)
         DO j = 1, n
            u(:, j) = EXP(CMPLX(0.0_dp, -evals(j), KIND=dp))*u(:, j)
         END DO
         u(:, :) = MATMUL(u, CONJG(TRANSPOSE(evec)))

         chc(:, :) = MATMUL(CONJG(TRANSPOSE(u)), MATMUL(base, u))

         DEALLOCATE (generator, evec, u, evals)
      END SUBROUTINE dense_rotated_chc

      SUBROUTINE dense_history_safe_occ_projector(chc, occupation_numbers, kpoint_weight, delta, occ_mat)
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: chc
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: occupation_numbers
      REAL(KIND=dp), INTENT(IN)                          :: kpoint_weight, delta
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(OUT)     :: occ_mat

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: evec, weighted_evec
      INTEGER                                            :: i, n
      REAL(KIND=dp)                                      :: projector_weight
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: evals, occ_scale

         n = SIZE(chc, 1)
         ALLOCATE (evec(n, n), weighted_evec(n, n), evals(n), occ_scale(n))

         occ_scale(:) = 2.0_dp*kpoint_weight*occupation_numbers(1:n)
         CALL diag_complex(chc, evec, evals)
         weighted_evec(:, :) = evec(:, :)
         DO i = 1, n
            weighted_evec(:, i) = occ_scale(i)*weighted_evec(:, i)
         END DO
         occ_mat(:, :) = MATMUL(weighted_evec, CONJG(TRANSPOSE(evec)))

         projector_weight = kpoint_ot_canonical_projector_weight(delta, ot_kp_canonical_projector_max_weight)
         occ_mat(:, :) = projector_weight*occ_mat(:, :)
         DO i = 1, n
            occ_mat(i, i) = CMPLX(REAL(occ_mat(i, i), KIND=dp) + &
                                  (1.0_dp - projector_weight)*occ_scale(i), 0.0_dp, KIND=dp)
         END DO

         DEALLOCATE (evec, weighted_evec, evals, occ_scale)
      END SUBROUTINE dense_history_safe_occ_projector

      SUBROUTINE dense_history_safe_occ_projector_derivative(chc, dchc, occupation_numbers, kpoint_weight, &
                                                             delta, docc_mat)
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: chc, dchc
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: occupation_numbers
      REAL(KIND=dp), INTENT(IN)                          :: kpoint_weight, delta
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(OUT)     :: docc_mat

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: evec, kernel, work
      INTEGER                                            :: i, j, n
      REAL(KIND=dp)                                      :: coeff, denom, gap_tol, occ_tol, projector_weight, &
                                                            scale
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: evals, occ_scale

         n = SIZE(chc, 1)
         docc_mat(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)

         projector_weight = kpoint_ot_canonical_projector_weight(delta, ot_kp_canonical_projector_max_weight)
         IF (projector_weight <= 0.0_dp) RETURN

         ALLOCATE (evec(n, n), kernel(n, n), work(n, n), evals(n), occ_scale(n))

         occ_scale(:) = 2.0_dp*kpoint_weight*occupation_numbers(1:n)
         CALL diag_complex(chc, evec, evals)
         work(:, :) = MATMUL(CONJG(TRANSPOSE(evec)), MATMUL(dchc, evec))

         scale = MAX(1.0_dp, MAXVAL(ABS(evals(:))))
         gap_tol = SQRT(EPSILON(1.0_dp))*scale
         occ_tol = SQRT(EPSILON(1.0_dp))*MAX(1.0_dp, MAXVAL(ABS(occ_scale(:))))
         kernel(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
         DO j = 1, n
            DO i = 1, n
               IF (i == j) CYCLE
               denom = evals(i) - evals(j)
               IF (ABS(denom) > gap_tol) THEN
                  coeff = (occ_scale(i) - occ_scale(j))/denom
               ELSE IF (ABS(occ_scale(i) - occ_scale(j)) <= occ_tol) THEN
                  coeff = 0.0_dp
               ELSE
                  coeff = (occ_scale(i) - occ_scale(j))/SIGN(gap_tol, denom)
               END IF
               kernel(i, j) = coeff*work(i, j)
            END DO
         END DO

         docc_mat(:, :) = projector_weight*MATMUL(evec, MATMUL(kernel, CONJG(TRANSPOSE(evec))))
         docc_mat(:, :) = 0.5_dp*(docc_mat(:, :) + CONJG(TRANSPOSE(docc_mat(:, :))))

         DEALLOCATE (evec, kernel, work, evals, occ_scale)

      END SUBROUTINE dense_history_safe_occ_projector_derivative

      SUBROUTINE dense_ref_rotation_gradient(x_re, x_im, base, occ_mat, grad_re, grad_im, eps_diag)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: x_re, x_im
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: base
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: occ_mat
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: grad_re, grad_im
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: eps_diag

      COMPLEX(KIND=dp)                                   :: kernel
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: dedu, evec, frechet, generator, inner, &
                                                            outer, u, work
      INTEGER                                            :: i, j, n
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: evals

         n = SIZE(x_re, 1)
         ALLOCATE (generator(n, n), evec(n, n), u(n, n), work(n, n), dedu(n, n))
         ALLOCATE (inner(n, n), outer(n, n), frechet(n, n), evals(n))

         generator(:, :) = CMPLX(-x_im(:, :), x_re(:, :), KIND=dp)
         CALL diag_complex(generator, evec, evals)

         u(:, :) = evec(:, :)
         DO j = 1, n
            u(:, j) = EXP(CMPLX(0.0_dp, -evals(j), KIND=dp))*u(:, j)
         END DO
         u(:, :) = MATMUL(u, CONJG(TRANSPOSE(evec)))

         dedu(:, :) = MATMUL(base, MATMUL(u, occ_mat))

         work(:, :) = MATMUL(CONJG(TRANSPOSE(u)), MATMUL(base, u))
         DO i = 1, n
            eps_diag(i) = REAL(work(i, i), KIND=dp)
         END DO

         inner(:, :) = MATMUL(CONJG(TRANSPOSE(evec)), MATMUL(dedu, evec))
         DO j = 1, n
            DO i = 1, n
               kernel = CONJG(qs_ot_complex_exp_frechet_kernel(evals(i), evals(j)))
               outer(i, j) = inner(i, j)*kernel
            END DO
         END DO
         frechet(:, :) = MATMUL(evec, MATMUL(outer, CONJG(TRANSPOSE(evec))))

         grad_re(:, :) = REAL(frechet, KIND=dp) - TRANSPOSE(REAL(frechet, KIND=dp))
         grad_im(:, :) = AIMAG(frechet) + TRANSPOSE(AIMAG(frechet))

         DEALLOCATE (generator, evec, u, work, dedu, inner, outer, frechet, evals)
      END SUBROUTINE dense_ref_rotation_gradient

      SUBROUTINE dbcsr_square_to_dense(matrix, dense)
      TYPE(dbcsr_type), POINTER                          :: matrix
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: dense

      INTEGER                                            :: bcol, brow, i, ig, j, jg
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_size, row_blk_size
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block
      TYPE(dbcsr_iterator_type)                          :: iter

         dense(:, :) = 0.0_dp
         NULLIFY (row_blk_size, col_blk_size)
         CALL dbcsr_get_info(matrix, row_blk_size=row_blk_size, col_blk_size=col_blk_size)
         CALL dbcsr_iterator_start(iter, matrix)
         DO WHILE (dbcsr_iterator_blocks_left(iter))
            CALL dbcsr_iterator_next_block(iter, brow, bcol, block)
            DO j = 1, SIZE(block, 2)
               jg = dbcsr_block_offset(col_blk_size, bcol) + j - 1
               IF (jg < 1 .OR. jg > SIZE(dense, 2)) CYCLE
               DO i = 1, SIZE(block, 1)
                  ig = dbcsr_block_offset(row_blk_size, brow) + i - 1
                  IF (ig >= 1 .AND. ig <= SIZE(dense, 1)) dense(ig, jg) = block(i, j)
               END DO
            END DO
         END DO
         CALL dbcsr_iterator_stop(iter)
      END SUBROUTINE dbcsr_square_to_dense

      SUBROUTINE dense_to_dbcsr_square(matrix, dense)
      TYPE(dbcsr_type), POINTER                          :: matrix
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: dense

      INTEGER                                            :: bcol, brow, i, ig, j, jg
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_size, row_blk_size
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block
      TYPE(dbcsr_iterator_type)                          :: iter

         NULLIFY (row_blk_size, col_blk_size)
         CALL dbcsr_get_info(matrix, row_blk_size=row_blk_size, col_blk_size=col_blk_size)
         CALL dbcsr_iterator_start(iter, matrix)
         DO WHILE (dbcsr_iterator_blocks_left(iter))
            CALL dbcsr_iterator_next_block(iter, brow, bcol, block)
            DO j = 1, SIZE(block, 2)
               jg = dbcsr_block_offset(col_blk_size, bcol) + j - 1
               DO i = 1, SIZE(block, 1)
                  ig = dbcsr_block_offset(row_blk_size, brow) + i - 1
                  IF (ig >= 1 .AND. ig <= SIZE(dense, 1) .AND. jg >= 1 .AND. jg <= SIZE(dense, 2)) &
                     block(i, j) = dense(ig, jg)
               END DO
            END DO
         END DO
         CALL dbcsr_iterator_stop(iter)
      END SUBROUTINE dense_to_dbcsr_square

      FUNCTION dbcsr_block_offset(block_size, block_index) RESULT(offset)
      INTEGER, DIMENSION(:), POINTER                     :: block_size
      INTEGER, INTENT(IN)                                :: block_index
      INTEGER                                            :: offset

         offset = 1
         IF (block_index > 1) offset = offset + SUM(block_size(1:block_index - 1))
      END FUNCTION dbcsr_block_offset

      FUNCTION kpoint_ot_smearing_g(occupation, eigenvalue, mu, sigma, maxocc, method) RESULT(gvalue)
      REAL(KIND=dp), INTENT(IN)                          :: occupation, eigenvalue, mu, sigma, maxocc
      INTEGER, INTENT(IN)                                :: method
      REAL(KIND=dp)                                      :: gvalue

      REAL(KIND=dp)                                      :: expu2, expx2, fi_norm, u, x

         IF (sigma <= 1.0E-14_dp .OR. maxocc <= EPSILON(maxocc)) THEN
            gvalue = 0.0_dp
            RETURN
         END IF

         x = (eigenvalue - mu)/sigma
         SELECT CASE (method)
         CASE (smear_fermi_dirac)
            fi_norm = occupation/maxocc
            gvalue = maxocc*fi_norm*(1.0_dp - fi_norm)/sigma
         CASE (smear_gaussian)
            expx2 = EXP(-x*x)
            gvalue = maxocc/(sigma*rootpi)*expx2
         CASE (smear_mp)
            expx2 = EXP(-x*x)
            gvalue = maxocc*(3.0_dp - 2.0_dp*x*x)/(2.0_dp*sigma*rootpi)*expx2
         CASE (smear_mv)
            u = x + sqrthalf
            expu2 = EXP(-u*u)
            gvalue = maxocc*(2.0_dp + sqrt2*x)/(sigma*rootpi)*expu2
         CASE DEFAULT
            gvalue = 0.0_dp
         END SELECT
      END FUNCTION kpoint_ot_smearing_g

   END SUBROUTINE qs_scf_loop_prepare_ot_kp_energy_gradient

! **************************************************************************************************
!> \brief the inner loop of scf, specific to using to the orbital transformation method
!>       basically, in goes the ks matrix out goes a new p matrix
!> \param qs_env ...
!> \param scf_env ...
!> \param smear ...
!> \param mos ...
!> \param rho ...
!> \param mo_derivs ...
!> \param total_energy ...
!> \param matrix_s ...
!> \param energy_only ...
!> \param has_unit_metric ...
!> \par History
!>      03.2006 created [Joost VandeVondele]
!>      2013    moved from qs_scf [Florian Schiffmann]
! **************************************************************************************************
   SUBROUTINE qs_scf_loop_do_ot(qs_env, scf_env, smear, mos, rho, mo_derivs, total_energy, &
                                matrix_s, energy_only, has_unit_metric, added_mos_auto_grow)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(smear_type), POINTER                          :: smear
      TYPE(mo_set_type), DIMENSION(:), INTENT(INOUT)     :: mos
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mo_derivs
      REAL(KIND=dp), INTENT(IN)                          :: total_energy
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      LOGICAL, INTENT(INOUT)                             :: energy_only
      LOGICAL, INTENT(IN)                                :: has_unit_metric
      LOGICAL, INTENT(OUT), OPTIONAL                     :: added_mos_auto_grow

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'qs_scf_loop_do_ot'

      INTEGER                                            :: handle, ispin
      LOGICAL                                            :: my_added_mos_auto_grow
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: rho_ao
      TYPE(dbcsr_type), POINTER                          :: orthogonality_metric
      TYPE(scf_control_type), POINTER                    :: scf_control

      CALL timeset(routineN, handle)
      NULLIFY (rho_ao)
      NULLIFY (scf_control)
      my_added_mos_auto_grow = .FALSE.
      IF (PRESENT(added_mos_auto_grow)) added_mos_auto_grow = .FALSE.

      CALL qs_rho_get(rho, rho_ao=rho_ao)
      CALL get_qs_env(qs_env=qs_env, scf_control=scf_control)

      IF (has_unit_metric) THEN
         NULLIFY (orthogonality_metric)
      ELSE
         orthogonality_metric => matrix_s(1)%matrix
      END IF

      ! in case of LSD the first spin qs_ot_env will drive the minimization
      ! in the case of a restricted calculation, it will make sure the spin orbitals are equal

      CALL ot_scf_mini(mos, mo_derivs, smear, orthogonality_metric, &
                       total_energy, energy_only, scf_env%iter_delta, &
                       scf_env%qs_ot_env)

      DO ispin = 1, SIZE(mos)
         CALL set_mo_occupation(mo_set=mos(ispin), smear=smear, &
                                emit_warnings=(.NOT. energy_only .AND. &
                                               .NOT. scf_env%qs_ot_env(1)%settings%do_ener))
      END DO
      IF (.NOT. energy_only) &
         CALL check_added_mos_auto_occupation(mos, scf_control, smear, my_added_mos_auto_grow)
      IF (my_added_mos_auto_grow) THEN
         IF (PRESENT(added_mos_auto_grow)) added_mos_auto_grow = .TRUE.
         CALL timestop(handle)
         RETURN
      END IF

      DO ispin = 1, SIZE(mos)
         CALL calculate_density_matrix(mos(ispin), &
                                       rho_ao(ispin)%matrix, &
                                       use_dbcsr=.TRUE.)
      END DO

      scf_env%iter_method = scf_env%qs_ot_env(1)%OT_METHOD_FULL
      scf_env%iter_param = scf_env%qs_ot_env(1)%ds_min
      qs_env%broyden_adaptive_sigma = scf_env%qs_ot_env(1)%broyden_adaptive_sigma

      CALL timestop(handle)

   END SUBROUTINE qs_scf_loop_do_ot

! **************************************************************************************************
!> \brief checks whether an automatic smearing virtual-space buffer is still sufficient
!> \param mos ...
!> \param scf_control ...
!> \param smear ...
! **************************************************************************************************
   SUBROUTINE check_added_mos_auto_occupation(mos, scf_control, smear, added_mos_auto_grow)

      TYPE(mo_set_type), DIMENSION(:), INTENT(IN)        :: mos
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(smear_type), POINTER                          :: smear
      LOGICAL, INTENT(OUT)                               :: added_mos_auto_grow

      INTEGER                                            :: ispin, nmo
      REAL(KIND=dp)                                      :: last_occ, threshold

      added_mos_auto_grow = .FALSE.
      IF (.NOT. ASSOCIATED(scf_control)) RETURN
      IF (.NOT. ANY(scf_control%added_mos_auto)) RETURN
      IF (.NOT. ASSOCIATED(smear)) RETURN
      IF (.NOT. smear%do_smear) RETURN

      threshold = MAX(smear%eps_fermi_dirac, 10.0_dp*EPSILON(1.0_dp))
      DO ispin = 1, MIN(SIZE(mos), SIZE(scf_control%added_mos_auto))
         IF (.NOT. scf_control%added_mos_auto(ispin)) CYCLE
         IF (.NOT. ASSOCIATED(mos(ispin)%occupation_numbers)) CYCLE

         nmo = mos(ispin)%nmo
         IF (nmo < 1) CYCLE
         last_occ = ABS(mos(ispin)%occupation_numbers(nmo))
         IF (last_occ > threshold) THEN
            IF (nmo >= mos(ispin)%nao) THEN
               CALL cp_abort(__LOCATION__, &
                             "ADDED_MOS AUTO exhausted the AO basis but the highest MO is still occupied. "// &
                             "Use a larger basis or a lower electronic temperature.")
            ELSE
               added_mos_auto_grow = .TRUE.
               RETURN
            END IF
         END IF
      END DO

   END SUBROUTINE check_added_mos_auto_occupation

! **************************************************************************************************
!> \brief Performs the requested density mixing if any needed
!> \param scf_env   Holds SCF environment information
!> \param rho       All data for the electron density
!> \param para_env  Parallel environment
!> \param diis_step Did we do a DIIS step?
! **************************************************************************************************
   SUBROUTINE qs_scf_density_mixing(scf_env, rho, para_env, diis_step)
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(mp_para_env_type), POINTER                    :: para_env
      LOGICAL                                            :: diis_step

      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rho_ao_kp

      NULLIFY (rho_ao_kp)

      CALL qs_rho_get(rho, rho_ao_kp=rho_ao_kp)

      SELECT CASE (scf_env%mixing_method)
      CASE (direct_mixing_nr)
         CALL scf_env_density_mixing(scf_env%p_mix_new, &
                                     scf_env%mixing_store, rho_ao_kp, para_env, scf_env%iter_delta, scf_env%iter_count, &
                                     diis=diis_step)
      CASE (gspace_mixing_nr, pulay_mixing_nr, broyden_mixing_nr, modified_broyden_mixing_nr, &
            multisecant_mixing_nr)
         ! Compute the difference p_out-p_in
         CALL self_consistency_check(rho_ao_kp, scf_env%p_delta, para_env, scf_env%p_mix_new, &
                                     delta=scf_env%iter_delta)
      CASE (no_mixing_nr)
      CASE DEFAULT
         CALL cp_abort(__LOCATION__, &
                       "unknown scf mixing method: "// &
                       cp_to_string(scf_env%mixing_method))
      END SELECT

   END SUBROUTINE qs_scf_density_mixing

! **************************************************************************************************
!> \brief checks whether exit conditions for outer loop are satisfied
!> \param qs_env ...
!> \param scf_env ...
!> \param scf_control ...
!> \param should_stop ...
!> \param outer_loop_converged ...
!> \param exit_outer_loop ...
! **************************************************************************************************
   SUBROUTINE qs_scf_check_outer_exit(qs_env, scf_env, scf_control, should_stop, &
                                      outer_loop_converged, exit_outer_loop)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(scf_control_type), POINTER                    :: scf_control
      LOGICAL                                            :: should_stop, outer_loop_converged, &
                                                            exit_outer_loop

      REAL(KIND=dp)                                      :: outer_loop_eps

      outer_loop_converged = .TRUE.
      IF (scf_control%outer_scf%have_scf) THEN
         ! We have an outer SCF loop...
         scf_env%outer_scf%iter_count = scf_env%outer_scf%iter_count + 1
         outer_loop_converged = .FALSE.

         CALL outer_loop_gradient(qs_env, scf_env)
         ! Multiple constraints: get largest deviation
         outer_loop_eps = SQRT(MAXVAL(scf_env%outer_scf%gradient(:, scf_env%outer_scf%iter_count)**2))

         IF (outer_loop_eps < scf_control%outer_scf%eps_scf) outer_loop_converged = .TRUE.
      END IF

      exit_outer_loop = should_stop .OR. outer_loop_converged .OR. &
                        scf_env%outer_scf%iter_count > scf_control%outer_scf%max_scf

   END SUBROUTINE qs_scf_check_outer_exit

! **************************************************************************************************
!> \brief checks whether exit conditions for inner loop are satisfied
!> \param qs_env ...
!> \param scf_env ...
!> \param scf_control ...
!> \param should_stop ...
!> \param just_energy ...
!> \param exit_inner_loop ...
!> \param inner_loop_converged ...
!> \param output_unit ...
! **************************************************************************************************
   SUBROUTINE qs_scf_check_inner_exit(qs_env, scf_env, scf_control, should_stop, just_energy, &
                                      exit_inner_loop, inner_loop_converged, output_unit)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(scf_control_type), POINTER                    :: scf_control
      LOGICAL                                            :: should_stop, just_energy, &
                                                            exit_inner_loop, inner_loop_converged
      INTEGER                                            :: output_unit

      inner_loop_converged = .FALSE.
      exit_inner_loop = .FALSE.

      CALL external_control(should_stop, "SCF", target_time=qs_env%target_time, &
                            start_time=qs_env%start_time)
      IF (scf_env%iter_delta < scf_control%eps_scf) THEN
         IF (output_unit > 0) THEN
            WRITE (UNIT=output_unit, FMT="(/,T3,A,I5,A/)") &
               "*** SCF run converged in ", scf_env%iter_count, " steps ***"
         END IF
         inner_loop_converged = .TRUE.
         exit_inner_loop = .TRUE.
      ELSE IF (should_stop .OR. scf_env%iter_count >= scf_control%max_scf) THEN
         inner_loop_converged = .FALSE.
         IF (just_energy) THEN
            exit_inner_loop = .FALSE.
         ELSE
            exit_inner_loop = .TRUE.
            IF (output_unit > 0) THEN
               WRITE (UNIT=output_unit, FMT="(/,T3,A,I5,A/)") &
                  "Leaving inner SCF loop after reaching ", scf_env%iter_count, " steps."
            END IF
         END IF
      END IF

   END SUBROUTINE qs_scf_check_inner_exit

! **************************************************************************************************
!> \brief undoing density mixing. Important upon convergence
!> \param scf_env ...
!> \param rho ...
!> \param dft_control ...
!> \param para_env ...
!> \param diis_step ...
! **************************************************************************************************
   SUBROUTINE qs_scf_undo_mixing(scf_env, rho, dft_control, para_env, diis_step)
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mp_para_env_type), POINTER                    :: para_env
      LOGICAL                                            :: diis_step

      CHARACTER(len=default_string_length)               :: name
      INTEGER                                            :: ic, ispin, nc
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rho_ao_kp

      NULLIFY (rho_ao_kp)

      IF (scf_env%mixing_method > 0) THEN
         CALL qs_rho_get(rho, rho_ao_kp=rho_ao_kp)
         nc = SIZE(scf_env%p_mix_new, 2)
         SELECT CASE (scf_env%mixing_method)
         CASE (direct_mixing_nr)
            CALL scf_env_density_mixing(scf_env%p_mix_new, scf_env%mixing_store, &
                                        rho_ao_kp, para_env, scf_env%iter_delta, &
                                        scf_env%iter_count, diis=diis_step, &
                                        invert=.TRUE.)
            DO ic = 1, nc
               DO ispin = 1, dft_control%nspins
                  CALL dbcsr_get_info(rho_ao_kp(ispin, ic)%matrix, name=name) ! keep the name
                  CALL dbcsr_copy(rho_ao_kp(ispin, ic)%matrix, scf_env%p_mix_new(ispin, ic)%matrix, name=name)
               END DO
            END DO
         CASE (gspace_mixing_nr, pulay_mixing_nr, broyden_mixing_nr, modified_broyden_mixing_nr, &
               multisecant_mixing_nr)
            DO ic = 1, nc
               DO ispin = 1, dft_control%nspins
                  CALL dbcsr_get_info(rho_ao_kp(ispin, ic)%matrix, name=name) ! keep the name
                  CALL dbcsr_copy(rho_ao_kp(ispin, ic)%matrix, scf_env%p_mix_new(ispin, ic)%matrix, name=name)
               END DO
            END DO
         END SELECT
      END IF
   END SUBROUTINE qs_scf_undo_mixing

! **************************************************************************************************
!> \brief Performs the updates rho (takes care of mixing as well)
!> \param rho ...
!> \param qs_env ...
!> \param scf_env ...
!> \param ks_env ...
!> \param mix_rho ...
! **************************************************************************************************
   SUBROUTINE qs_scf_rho_update(rho, qs_env, scf_env, ks_env, mix_rho)
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      LOGICAL, INTENT(IN)                                :: mix_rho

      TYPE(mp_para_env_type), POINTER                    :: para_env

      NULLIFY (para_env)
      CALL get_qs_env(qs_env, para_env=para_env)
      ! ** update qs_env%rho
      CALL qs_rho_update_rho(rho, qs_env=qs_env)
      ! ** Density mixing through density matrix or on the reciprocal space grid (exclusive)
      IF (mix_rho) THEN
         CALL gspace_mixing(qs_env, scf_env%mixing_method, scf_env%mixing_store, rho, &
                            para_env, scf_env%iter_count)

      END IF
      CALL qs_ks_did_change(ks_env, rho_changed=.TRUE.)

   END SUBROUTINE qs_scf_rho_update

! **************************************************************************************************
!> \brief Performs the necessary steps before leaving innner scf loop
!> \param scf_env ...
!> \param qs_env ...
!> \param diis_step ...
!> \param output_unit ...
! **************************************************************************************************
   SUBROUTINE qs_scf_inner_finalize(scf_env, qs_env, diis_step, output_unit)
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL                                            :: diis_step
      INTEGER, INTENT(IN)                                :: output_unit

      LOGICAL                                            :: do_kpoints
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho

      NULLIFY (energy, rho, dft_control, ks_env)

      CALL get_qs_env(qs_env=qs_env, energy=energy, ks_env=ks_env, &
                      rho=rho, dft_control=dft_control, para_env=para_env, &
                      do_kpoints=do_kpoints)

      CALL cleanup_scf_loop(scf_env)

      ! now, print out energies and charges corresponding to the obtained wfn
      ! (this actually is not 100% consistent at this point)!
      CALL qs_scf_print_summary(output_unit, qs_env)

      CALL qs_scf_undo_mixing(scf_env, rho, dft_control, para_env, diis_step)

      !   *** update rspace rho since the mo changed
      !   *** this might not always be needed (i.e. no post calculation / no forces )
      !   *** but guarantees that rho and wfn are consistent at this point
      CALL qs_scf_rho_update(rho, qs_env, scf_env, ks_env, mix_rho=.FALSE.)

   END SUBROUTINE qs_scf_inner_finalize

! **************************************************************************************************
!> \brief perform cleanup operations at the end of an scf loop
!> \param scf_env ...
!> \par History
!>      03.2006 created [Joost VandeVondele]
! **************************************************************************************************
   SUBROUTINE cleanup_scf_loop(scf_env)
      TYPE(qs_scf_env_type), INTENT(INOUT)               :: scf_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'cleanup_scf_loop'

      INTEGER                                            :: handle, ispin
      LOGICAL                                            :: keep_kpoint_ot_for_restart

      CALL timeset(routineN, handle)

      SELECT CASE (scf_env%method)
      CASE (ot_method_nr)
         keep_kpoint_ot_for_restart = ANY(scf_env%qs_ot_env(:)%kpoint_index > 0)
         IF (.NOT. keep_kpoint_ot_for_restart) THEN
            DO ispin = 1, SIZE(scf_env%qs_ot_env)
               CALL ot_scf_destroy(scf_env%qs_ot_env(ispin))
            END DO
            DEALLOCATE (scf_env%qs_ot_env)
         END IF
      CASE (ot_diag_method_nr)
         !
      CASE (general_diag_method_nr)
         !
      CASE (special_diag_method_nr)
         !
      CASE (block_krylov_diag_method_nr, block_davidson_diag_method_nr)
         !
      CASE (filter_matrix_diag_method_nr)
         !
      CASE (smeagol_method_nr)
         !
      CASE DEFAULT
         CALL cp_abort(__LOCATION__, &
                       "unknown scf method method:"// &
                       cp_to_string(scf_env%method))
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE cleanup_scf_loop

END MODULE qs_scf_loop_utils
