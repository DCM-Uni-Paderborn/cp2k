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
        dbcsr_multiply, dbcsr_p_type, dbcsr_release, dbcsr_set, dbcsr_transposed, dbcsr_type, &
        dbcsr_type_antisymmetric, dbcsr_type_symmetric
   USE cp_dbcsr_contrib,                ONLY: dbcsr_get_diag,&
                                              dbcsr_scale_by_vector
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              cp_dbcsr_sm_fm_multiply
   USE cp_external_control,             ONLY: external_control
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_to_string
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
   USE message_passing,                 ONLY: mp_para_env_type
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
   USE qs_ot,                           ONLY: qs_ot_get_orbitals_ref_complex
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
   USE tblite_interface,                ONLY: tb_native_scc_mixer_active
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_scf_loop_utils'

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
      LOGICAL                                            :: use_real_wfn
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: orbital_scaling_factor
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvalues, occupation_numbers, wkp
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: xkp
      TYPE(cp_fm_type)                                   :: hc_im, hc_re
      TYPE(cp_fm_type), DIMENSION(:), POINTER            :: fmwork
      TYPE(cp_fm_type), POINTER                          :: mo_coeff, mo_coeff_im
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_hc, matrix_hc_im
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

      CALL timeset(routineN, handle)

      NULLIFY (cell_to_index, dft_control, energy, fmwork, kpoints, local_ot_env, matrix_hc, &
               matrix_hc_im, matrix_k_im, matrix_k_re, matrix_s_im, matrix_s_re, mo_coeff, &
               mo_coeff_im, occupation_numbers, para_env_inter_kp, sab_nl, wkp, xkp)
      added_mos_auto_grow = .FALSE.

      CPASSERT(ASSOCIATED(scf_env%qs_ot_env))
      CALL get_qs_env(qs_env=qs_env, dft_control=dft_control, energy=energy, kpoints=kpoints)
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
      local_ot_env(1)%etotal = energy%total

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
         DO ispin = 1, nspin_ot
            local_channel = qs_ot_channel_index(ispin, ikpoint, nspin_ot) - first_channel + 1
            CALL dbcsr_set(matrix_k_re, 0.0_dp)
            CALL dbcsr_set(matrix_k_im, 0.0_dp)
            CALL rskp_transform(rmatrix=matrix_k_re, cmatrix=matrix_k_im, rsmat=matrix_ks, &
                                ispin=ispin, xkp=xkp(1:3, ikpoint), &
                                cell_to_index=cell_to_index, sab_nl=sab_nl)

            CALL get_mo_set(kp%mos(1, ispin), mo_coeff=mo_coeff, nmo=nmo, &
                            occupation_numbers=occupation_numbers)
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

            CALL qs_scf_loop_prepare_ot_kp_channel_gradient(local_ot_env(local_channel), &
                                                            matrix_hc(local_channel)%matrix, &
                                                            matrix_hc_im(local_channel)%matrix, &
                                                            occupation_numbers, wkp(ikpoint))

            ALLOCATE (orbital_scaling_factor(SIZE(occupation_numbers)))
            IF (local_ot_env(local_channel)%settings%occupation_preconditioner) THEN
               orbital_scaling_factor(:) = 2.0_dp*wkp(ikpoint)
            ELSE
               orbital_scaling_factor(:) = 2.0_dp*wkp(ikpoint)*occupation_numbers(:)
            END IF
            CALL dbcsr_scale_by_vector(matrix_hc(local_channel)%matrix, &
                                       alpha=orbital_scaling_factor, side='right')
            CALL dbcsr_scale_by_vector(matrix_hc_im(local_channel)%matrix, &
                                       alpha=orbital_scaling_factor, side='right')
            DEALLOCATE (orbital_scaling_factor)

            CALL cp_fm_release(hc_re)
            CALL cp_fm_release(hc_im)
         END DO
      END DO

      CALL qs_scf_loop_prepare_ot_kp_energy_gradient(kpoints, local_ot_env, first_channel, nspin_ot, &
                                                     kp_range, wkp, para_env_inter_kp, scf_control%smear)

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
            CALL get_mo_set(kp%mos(1, ispin), mo_coeff=mo_coeff, eigenvalues=eigenvalues)
            CALL get_mo_set(kp%mos(2, ispin), mo_coeff=mo_coeff_im)
            CALL qs_ot_get_orbitals_ref_complex(matrix_hc(local_channel)%matrix, &
                                                matrix_hc_im(local_channel)%matrix, &
                                                matrix_s_re, matrix_s_im, &
                                                local_ot_env(local_channel))
            CALL copy_dbcsr_to_fm(matrix_hc(local_channel)%matrix, mo_coeff)
            CALL copy_dbcsr_to_fm(matrix_hc_im(local_channel)%matrix, mo_coeff_im)
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
         fmwork => scf_env%scf_work1
         CALL kpoint_density_matrices(kpoints)
         CALL kpoint_density_transform(kpoints, scf_env%p_mix_new, .FALSE., &
                                       matrix_s(1, 1)%matrix, sab_nl, fmwork, &
                                       qs_ot_env=scf_env%qs_ot_env)
         scf_env%iter_method = scf_env%qs_ot_env(first_channel)%OT_METHOD_FULL
         scf_env%iter_param = scf_env%qs_ot_env(first_channel)%ds_min
         qs_env%broyden_adaptive_sigma = scf_env%qs_ot_env(first_channel)%broyden_adaptive_sigma
      END IF

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
!> \brief prepare complex rotation and band-energy gradients for one local k-point OT channel
!> \param qs_ot_env ...
!> \param matrix_hc ...
!> \param matrix_hc_im ...
!> \param occupation_numbers ...
!> \param kpoint_weight ...
! **************************************************************************************************
   SUBROUTINE qs_scf_loop_prepare_ot_kp_channel_gradient(qs_ot_env, matrix_hc, matrix_hc_im, &
                                                         occupation_numbers, kpoint_weight)

      TYPE(qs_ot_type)                                   :: qs_ot_env
      TYPE(dbcsr_type), POINTER                          :: matrix_hc, matrix_hc_im
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: occupation_numbers
      REAL(KIND=dp), INTENT(IN)                          :: kpoint_weight

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: orbital_scaling_factor
      TYPE(dbcsr_type), POINTER                          :: b_im, b_re, c_im, c_re, chc_im, chc_re, &
                                                            tmp, tmp2

      IF (.NOT. qs_ot_env%settings%do_rotation) RETURN

      CPASSERT(qs_ot_env%has_complex_kpoint_state)
      CPASSERT(ASSOCIATED(matrix_hc))
      CPASSERT(ASSOCIATED(matrix_hc_im))
      CPASSERT(ASSOCIATED(qs_ot_env%rot_mat_chc))
      CPASSERT(ASSOCIATED(qs_ot_env%rot_mat_chc_im))

      c_re => qs_ot_env%matrix_x
      c_im => qs_ot_env%matrix_x_im
      chc_re => qs_ot_env%rot_mat_chc
      chc_im => qs_ot_env%rot_mat_chc_im
      tmp => qs_ot_env%buf3_k_k_nosym
      tmp2 => qs_ot_env%buf4_k_k_nosym

      ! CHC = C(k)^H * H(k)C(k), kept unoccupied for the band-energy gradient.
      CALL dbcsr_multiply('T', 'N', 1.0_dp, c_re, matrix_hc, 0.0_dp, chc_re)
      CALL dbcsr_multiply('T', 'N', 1.0_dp, c_im, matrix_hc_im, 0.0_dp, tmp)
      CALL dbcsr_add(chc_re, tmp, alpha_scalar=1.0_dp, beta_scalar=1.0_dp)

      CALL dbcsr_multiply('T', 'N', 1.0_dp, c_re, matrix_hc_im, 0.0_dp, chc_im)
      CALL dbcsr_multiply('T', 'N', 1.0_dp, c_im, matrix_hc, 0.0_dp, tmp)
      CALL dbcsr_add(chc_im, tmp, alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)

      IF (qs_ot_env%settings%do_ener) THEN
         CALL dbcsr_get_diag(chc_re, qs_ot_env%ener_gx)
      END IF

      b_re => qs_ot_env%matrix_buf1
      b_im => qs_ot_env%buf2_k_k_nosym
      CALL dbcsr_copy(b_re, chc_re)
      CALL dbcsr_copy(b_im, chc_im)

      ALLOCATE (orbital_scaling_factor(SIZE(occupation_numbers)))
      orbital_scaling_factor(:) = 2.0_dp*kpoint_weight*occupation_numbers(:)
      CALL dbcsr_scale_by_vector(b_re, alpha=orbital_scaling_factor, side='right')
      CALL dbcsr_scale_by_vector(b_im, alpha=orbital_scaling_factor, side='right')
      DEALLOCATE (orbital_scaling_factor)

      ! dF/dU for the complex occupied-space rotation U(k).  This mirrors the Gamma path, but keeps
      ! the real and imaginary channels separate so a later full anti-Hermitian Frechet derivative can
      ! use the complete payload.
      CALL dbcsr_multiply('N', 'N', 1.0_dp, qs_ot_env%rot_mat_u, b_re, 0.0_dp, qs_ot_env%rot_mat_dedu)
      IF (ASSOCIATED(qs_ot_env%rot_mat_u_im)) THEN
         CALL dbcsr_multiply('N', 'N', 1.0_dp, qs_ot_env%rot_mat_u_im, b_im, 0.0_dp, tmp)
         CALL dbcsr_add(qs_ot_env%rot_mat_dedu, tmp, alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)

         CALL dbcsr_multiply('N', 'N', 1.0_dp, qs_ot_env%rot_mat_u, b_im, 0.0_dp, qs_ot_env%rot_mat_dedu_im)
         CALL dbcsr_multiply('N', 'N', 1.0_dp, qs_ot_env%rot_mat_u_im, b_re, 0.0_dp, tmp2)
         CALL dbcsr_add(qs_ot_env%rot_mat_dedu_im, tmp2, alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
      END IF

      CALL dbcsr_transposed(tmp, qs_ot_env%rot_mat_dedu, shallow_data_copy=.FALSE., &
                            transpose_distribution=.FALSE.)
      CALL dbcsr_copy(qs_ot_env%rot_mat_gx, tmp)
      CALL dbcsr_add(qs_ot_env%rot_mat_gx, qs_ot_env%rot_mat_dedu, alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)

      IF (ASSOCIATED(qs_ot_env%rot_mat_gx_im) .AND. ASSOCIATED(qs_ot_env%rot_mat_dedu_im)) THEN
         CALL dbcsr_transposed(tmp2, qs_ot_env%rot_mat_dedu_im, shallow_data_copy=.FALSE., &
                               transpose_distribution=.FALSE.)
         CALL dbcsr_copy(qs_ot_env%rot_mat_gx_im, tmp2)
         CALL dbcsr_add(qs_ot_env%rot_mat_gx_im, qs_ot_env%rot_mat_dedu_im, &
                        alpha_scalar=-1.0_dp, beta_scalar=-1.0_dp)
      END IF

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
                                                        kp_range, wkp, para_env_inter_kp, smear)

      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: local_ot_env
      INTEGER, INTENT(IN)                                :: first_channel, nspin_ot
      INTEGER, DIMENSION(2), INTENT(IN)                  :: kp_range
      REAL(KIND=dp), DIMENSION(:), POINTER               :: wkp
      TYPE(mp_para_env_type), POINTER                    :: para_env_inter_kp
      TYPE(smear_type), INTENT(IN)                       :: smear

      INTEGER                                            :: group, ikpoint, imo, ispin, local_channel, &
                                                            local_kpoint, ngroups, nmo
      REAL(KIND=dp)                                      :: curvature, gvalue, maxocc, mu, sigma, v
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: gdotv, gsum
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvalues, occupation_numbers
      TYPE(kpoint_env_type), POINTER                     :: kp

      IF (.NOT. local_ot_env(1)%settings%do_ener) RETURN

      ngroups = 1
      IF (nspin_ot == 2 .AND. smear%fixed_mag_mom > 0.0_dp) ngroups = 2
      ALLOCATE (gsum(ngroups), gdotv(ngroups))
      gsum(:) = 0.0_dp
      gdotv(:) = 0.0_dp

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

      DO local_kpoint = 1, SIZE(kpoints%kp_env)
         kp => kpoints%kp_env(local_kpoint)%kpoint_env
         ikpoint = kp%nkpoint
         DO ispin = 1, nspin_ot
            group = MERGE(ispin, 1, ngroups == nspin_ot)
            local_channel = qs_ot_channel_index(ispin, ikpoint, nspin_ot) - first_channel + 1
            CALL get_mo_set(kp%mos(1, ispin), maxocc=maxocc, nmo=nmo, eigenvalues=eigenvalues, &
                            occupation_numbers=occupation_numbers, mu=mu)
            DO imo = 1, nmo
               gvalue = kpoint_ot_smearing_g(occupation_numbers(imo), eigenvalues(imo), mu, &
                                             sigma, maxocc, smear%method)
               v = local_ot_env(local_channel)%ener_gx(imo) - eigenvalues(imo) + mu
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
            CALL get_mo_set(kp%mos(1, ispin), maxocc=maxocc, nmo=nmo, eigenvalues=eigenvalues, &
                            occupation_numbers=occupation_numbers, mu=mu)
            DO imo = 1, nmo
               gvalue = kpoint_ot_smearing_g(occupation_numbers(imo), eigenvalues(imo), mu, &
                                             sigma, maxocc, smear%method)
               IF (ABS(gsum(group)) > EPSILON(gsum(group))) THEN
                  v = local_ot_env(local_channel)%ener_gx(imo) - eigenvalues(imo) + mu
                  local_ot_env(local_channel)%ener_gx(imo) = &
                     wkp(ikpoint)*gvalue*(gdotv(group)/gsum(group) - v)
               ELSE
                  local_ot_env(local_channel)%ener_gx(imo) = 0.0_dp
               END IF

               IF (local_ot_env(local_channel)%settings%occupation_preconditioner .AND. &
                   smear%do_smear .AND. smear%method == smear_fermi_dirac) THEN
                  curvature = wkp(ikpoint)*gvalue
                  IF (curvature > 1.0_dp) &
                     local_ot_env(local_channel)%ener_gx(imo) = &
                     local_ot_env(local_channel)%ener_gx(imo)/curvature
               END IF
            END DO
         END DO
      END DO

      DEALLOCATE (gsum, gdotv)

      CPASSERT(kp_range(2) - kp_range(1) + 1 == SIZE(kpoints%kp_env))

   CONTAINS

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

      CALL timeset(routineN, handle)

      SELECT CASE (scf_env%method)
      CASE (ot_method_nr)
         DO ispin = 1, SIZE(scf_env%qs_ot_env)
            CALL ot_scf_destroy(scf_env%qs_ot_env(ispin))
         END DO
         DEALLOCATE (scf_env%qs_ot_env)
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
