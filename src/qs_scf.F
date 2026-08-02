!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines for the Quickstep SCF run.
!> \par History
!>      - Joost VandeVondele (02.2002)
!>           added code for: incremental (pab and gvg) update
!>                            initialisation (init_cube, l_info)
!>      - Joost VandeVondele (02.2002)
!>           called the poisson code of the classical part
!>           this takes into account the spherical cutoff and allows for
!>           isolated systems
!>      - Joost VandeVondele (02.2002)
!>           added multiple grid feature
!>           changed to spherical cutoff consistently (?)
!>           therefore removed the gradient correct functionals
!>      - updated with the new QS data structures (10.04.02,MK)
!>      - copy_matrix replaced by transfer_matrix (11.04.02,MK)
!>      - nrebuild_rho and nrebuild_gvg unified (12.04.02,MK)
!>      - set_mo_occupation for smearing of the MO occupation numbers
!>        (17.04.02,MK)
!>      - MO level shifting added (22.04.02,MK)
!>      - Usage of TYPE mo_set_p_type
!>      - Joost VandeVondele (05.2002)
!>            added cholesky based diagonalisation
!>      - 05.2002 added pao method [fawzi]
!>      - parallel FFT (JGH 22.05.2002)
!>      - 06.2002 moved KS matrix construction to qs_build_KS_matrix.F [fawzi]
!>      - started to include more LSD (01.2003,Joost VandeVondele)
!>      - 02.2003 scf_env [fawzi]
!>      - got rid of nrebuild (01.2004, Joost VandeVondele)
!>      - 10.2004 removed pao [fawzi]
!>      - 03.2006 large cleaning action [Joost VandeVondele]
!>      - High-spin ROKS added (05.04.06,MK)
!>      - Mandes (10.2013)
!>        intermediate energy communication with external communicator added
!>      - kpoints (08.2014, JGH)
!>      - unified k-point and gamma-point code (2014.11) [Ole Schuett]
!>      - added extra SCF loop for CDFT constraints (12.2015) [Nico Holmberg]
!> \author Matthias Krack (30.04.2001)
! **************************************************************************************************
MODULE qs_scf
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE cell_types,                      ONLY: cell_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_add, dbcsr_copy, dbcsr_create, dbcsr_deallocate_matrix, dbcsr_get_info, &
        dbcsr_init_p, dbcsr_multiply, dbcsr_p_type, dbcsr_release_p, dbcsr_set, dbcsr_type, &
        dbcsr_type_no_symmetry
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              copy_fm_to_dbcsr_bc,&
                                              cp_dbcsr_m_by_n_from_row_template,&
                                              cp_dbcsr_m_by_n_from_template,&
                                              cp_dbcsr_sm_fm_multiply,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_files,                        ONLY: close_file
   USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_init_random,&
                                              cp_fm_release,&
                                              cp_fm_set_submatrix,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_add_default_logger,&
                                              cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_release,&
                                              cp_logger_type,&
                                              cp_rm_default_logger,&
                                              cp_to_string
   USE cp_output_handling,              ONLY: cp_add_iter_level,&
                                              cp_iterate,&
                                              cp_p_file,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr,&
                                              cp_rm_iter_level
   USE cp_result_methods,               ONLY: get_results,&
                                              test_for_result
   USE cp_result_types,                 ONLY: cp_result_type
   USE ec_env_types,                    ONLY: energy_correction_type
   USE input_constants,                 ONLY: &
        broyden_type_1, broyden_type_1_explicit, broyden_type_1_explicit_ls, broyden_type_1_ls, &
        broyden_type_2, broyden_type_2_explicit, broyden_type_2_explicit_ls, broyden_type_2_ls, &
        cdft2ot, history_guess, ot2cdft, ot_precond_full_all, ot_precond_full_single, &
        ot_precond_full_single_inverse, ot_precond_none, ot_precond_s_inverse, &
        outer_scf_becke_constraint, outer_scf_hirshfeld_constraint, outer_scf_optimizer_broyden, &
        outer_scf_optimizer_newton_ls, tblite_scc_mixer_tblite
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp
   USE kpoint_io,                       ONLY: write_kpoints_restart
   USE kpoint_methods,                  ONLY: kpoint_initialize_mo_set,&
                                              kpoint_initialize_mos
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_env_type,&
                                              kpoint_type
   USE machine,                         ONLY: m_flush,&
                                              m_walltime
   USE mathlib,                         ONLY: invert_matrix
   USE message_passing,                 ONLY: mp_comm_type,&
                                              mp_para_env_type
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: evolt
   USE preconditioner,                  ONLY: make_preconditioner,&
                                              make_preconditioner_complex_full_all,&
                                              prepare_preconditioner,&
                                              restart_preconditioner
   USE preconditioner_types,            ONLY: destroy_preconditioner,&
                                              init_preconditioner
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE qs_block_davidson_types,         ONLY: block_davidson_deallocate
   USE qs_cdft_scf_utils,               ONLY: build_diagonal_jacobian,&
                                              create_tmp_logger,&
                                              initialize_inverse_jacobian,&
                                              prepare_jacobian_stencil,&
                                              print_inverse_jacobian,&
                                              restart_inverse_jacobian
   USE qs_cdft_types,                   ONLY: cdft_control_type
   USE qs_charge_mixing,                ONLY: charge_mixing_scc_error
   USE qs_charges_types,                ONLY: qs_charges_type
   USE qs_density_matrices,             ONLY: calculate_density_matrix
   USE qs_density_mixing_types,         ONLY: gspace_mixing_nr
   USE qs_diis,                         ONLY: qs_diis_b_clear,&
                                              qs_diis_b_clear_kp,&
                                              qs_diis_b_create,&
                                              qs_diis_b_create_kp
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type,&
                                              set_qs_env
   USE qs_integrate_potential,          ONLY: integrate_v_rspace
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_ks_methods,                   ONLY: evaluate_core_matrix_traces,&
                                              qs_ks_update_qs_env
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_did_change,&
                                              qs_ks_env_type
   USE qs_matrix_pools,                 ONLY: mpools_get,&
                                              mpools_rebuild_fm_pools,&
                                              mpools_release
   USE qs_mo_io,                        ONLY: write_mo_set_to_restart
   USE qs_mo_methods,                   ONLY: make_basis_simple,&
                                              make_basis_sm
   USE qs_mo_occupation,                ONLY: set_mo_occupation
   USE qs_mo_types,                     ONLY: allocate_mo_set,&
                                              deallocate_mo_set,&
                                              duplicate_mo_set,&
                                              get_mo_set,&
                                              init_mo_set,&
                                              mo_set_restrict,&
                                              mo_set_type,&
                                              reassign_allocated_mos
   USE qs_ot,                           ONLY: qs_ot_get_orbitals_ref_complex,&
                                              qs_ot_get_p,&
                                              qs_ot_new_preconditioner
   USE qs_ot_scf,                       ONLY: init_ot_smearing_energies,&
                                              init_ot_smearing_frontier_energies,&
                                              ot_scf_destroy,&
                                              ot_scf_init,&
                                              ot_scf_read_input
   USE qs_ot_types,                     ONLY: qs_ot_allocate,&
                                              qs_ot_channel_index,&
                                              qs_ot_check_channel_context,&
                                              qs_ot_init,&
                                              qs_ot_number_of_channels,&
                                              qs_ot_reset_accepted_model_calibration,&
                                              qs_ot_set_context,&
                                              qs_ot_type
   USE qs_outer_scf,                    ONLY: outer_loop_gradient,&
                                              outer_loop_optimize,&
                                              outer_loop_purge_history,&
                                              outer_loop_switch,&
                                              outer_loop_update_qs_env
   USE qs_rho_methods,                  ONLY: qs_rho_update_rho
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_scf_diagonalization,          ONLY: do_general_diag_kp
   USE qs_scf_initialization,           ONLY: qs_scf_env_initialize
   USE qs_scf_loop_utils,               ONLY: &
        qs_scf_check_inner_exit, qs_scf_check_outer_exit, qs_scf_density_mixing, &
        qs_scf_inner_finalize, qs_scf_kp_added_mos_auto_post_resize_steps, &
        qs_scf_kp_rebased_calibration_samples, qs_scf_new_mos, qs_scf_new_mos_kp, &
        qs_scf_rho_update, qs_scf_set_loop_flags
   USE qs_scf_output,                   ONLY: qs_scf_cdft_info,&
                                              qs_scf_cdft_initial_info,&
                                              qs_scf_gce_info,&
                                              qs_scf_loop_info,&
                                              qs_scf_loop_print,&
                                              qs_scf_outer_loop_info,&
                                              qs_scf_write_mos
   USE qs_scf_post_scf,                 ONLY: qs_scf_compute_properties
   USE qs_scf_types,                    ONLY: &
        block_davidson_diag_method_nr, block_krylov_diag_method_nr, filter_matrix_diag_method_nr, &
        general_diag_method_nr, ot_diag_method_nr, ot_method_nr, qs_ot_restart_dense_matrix_type, &
        qs_scf_env_type, qs_scf_ot_restart_clear, smeagol_method_nr, special_diag_method_nr
   USE qs_wf_history_methods,           ONLY: wfi_purge_history,&
                                              wfi_update
   USE scf_control_types,               ONLY: scf_control_type
   USE smeagol_interface,               ONLY: run_smeagol_bulktrans,&
                                              run_smeagol_emtrans
   USE tblite_interface,                ONLY: tb_get_energy,&
                                              tb_native_scc_mixer_active,&
                                              tb_scf_mixer_error,&
                                              tb_update_charges
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_scf'
   LOGICAL, PRIVATE                     :: reuse_precond = .FALSE.
   LOGICAL, PRIVATE                     :: used_history = .FALSE.

   PUBLIC :: scf, scf_env_cleanup, scf_env_do_scf, cdft_scf, init_scf_loop

CONTAINS

! **************************************************************************************************
!> \brief perform an scf procedure in the given qs_env
!> \param qs_env the qs_environment where to perform the scf procedure
!> \param has_converged ...
!> \param total_scf_steps ...
!> \par History
!>      02.2003 introduced scf_env, moved real work to scf_env_do_scf [fawzi]
!> \author fawzi
!> \note
! **************************************************************************************************
   SUBROUTINE scf(qs_env, has_converged, total_scf_steps)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(OUT), OPTIONAL                     :: has_converged
      INTEGER, INTENT(OUT), OPTIONAL                     :: total_scf_steps

      INTEGER                                            :: ihistory, max_scf_tmp, tsteps
      LOGICAL                                            :: converged, outer_scf_loop, should_stop
      LOGICAL, SAVE                                      :: first_step_flag = .TRUE.
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: gradient_history, variable_history
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(section_vals_type), POINTER                   :: dft_section, input, scf_section

      NULLIFY (scf_env)
      logger => cp_get_default_logger()
      CPASSERT(ASSOCIATED(qs_env))
      IF (PRESENT(has_converged)) THEN
         has_converged = .FALSE.
      END IF
      IF (PRESENT(total_scf_steps)) THEN
         total_scf_steps = 0
      END IF
      CALL get_qs_env(qs_env, scf_env=scf_env, input=input, &
                      dft_control=dft_control, scf_control=scf_control)
      IF (scf_control%max_scf > 0) THEN

         dft_section => section_vals_get_subs_vals(input, "DFT")
         scf_section => section_vals_get_subs_vals(dft_section, "SCF")

         IF (.NOT. ASSOCIATED(scf_env)) THEN
            CALL qs_scf_env_initialize(qs_env, scf_env)
            ! Moved here from qs_scf_env_initialize to be able to have more scf_env
            CALL set_qs_env(qs_env, scf_env=scf_env)
         ELSE
            CALL qs_scf_env_initialize(qs_env, scf_env)
         END IF

         IF ((scf_control%density_guess == history_guess) .AND. (first_step_flag)) THEN
            max_scf_tmp = scf_control%max_scf
            scf_control%max_scf = 1
            outer_scf_loop = scf_control%outer_scf%have_scf
            scf_control%outer_scf%have_scf = .FALSE.
         END IF

         IF (.NOT. dft_control%qs_control%cdft) THEN
            CALL scf_env_do_scf(scf_env=scf_env, scf_control=scf_control, qs_env=qs_env, &
                                converged=converged, should_stop=should_stop, total_scf_steps=tsteps)
         ELSE
            ! Third SCF loop needed for CDFT with OT to properly restart OT inner loop
            CALL cdft_scf(qs_env=qs_env, should_stop=should_stop)
         END IF

         ! If SCF has not converged, then we should not start MP2
         IF (ASSOCIATED(qs_env%mp2_env)) qs_env%mp2_env%hf_fail = .NOT. converged

         ! Add the converged outer_scf SCF gradient(s)/variable(s) to history
         IF (scf_control%outer_scf%have_scf) THEN
            ihistory = scf_env%outer_scf%iter_count
            CALL get_qs_env(qs_env, gradient_history=gradient_history, &
                            variable_history=variable_history)
            ! We only store the latest two values
            gradient_history(:, 1) = gradient_history(:, 2)
            gradient_history(:, 2) = scf_env%outer_scf%gradient(:, ihistory)
            variable_history(:, 1) = variable_history(:, 2)
            variable_history(:, 2) = scf_env%outer_scf%variables(:, ihistory)
            ! Reset flag
            IF (used_history) used_history = .FALSE.
            ! Update a counter and check if the Jacobian should be deallocated
            IF (ASSOCIATED(scf_env%outer_scf%inv_jacobian)) THEN
               scf_control%outer_scf%cdft_opt_control%ijacobian(2) = scf_control%outer_scf%cdft_opt_control%ijacobian(2) + 1
               IF (scf_control%outer_scf%cdft_opt_control%ijacobian(2) >= &
                   scf_control%outer_scf%cdft_opt_control%jacobian_freq(2) .AND. &
                   scf_control%outer_scf%cdft_opt_control%jacobian_freq(2) > 0) THEN
                  scf_env%outer_scf%deallocate_jacobian = .TRUE.
               END IF
            END IF
         END IF
         !   *** add the converged wavefunction to the wavefunction history
         IF ((ASSOCIATED(qs_env%wf_history)) .AND. &
             ((scf_control%density_guess /= history_guess) .OR. &
              (.NOT. first_step_flag))) THEN
            IF (.NOT. dft_control%qs_control%cdft) THEN
               CALL wfi_update(qs_env%wf_history, qs_env=qs_env, dt=1.0_dp)
            ELSE
               IF (dft_control%qs_control%cdft_control%should_purge) THEN
                  CALL wfi_purge_history(qs_env)
                  CALL outer_loop_purge_history(qs_env)
                  dft_control%qs_control%cdft_control%should_purge = .FALSE.
               ELSE
                  CALL wfi_update(qs_env%wf_history, qs_env=qs_env, dt=1.0_dp)
               END IF
            END IF
         ELSE IF ((scf_control%density_guess == history_guess) .AND. &
                  (first_step_flag)) THEN
            scf_control%max_scf = max_scf_tmp
            scf_control%outer_scf%have_scf = outer_scf_loop
            first_step_flag = .FALSE.
         END IF

         ! *** compute properties that depend on the converged wavefunction
         IF (.NOT. (should_stop)) CALL qs_scf_compute_properties(qs_env)

         ! *** SMEAGOL interface ***
         IF (.NOT. (should_stop)) THEN
            ! compute properties that depend on the converged wavefunction ..
            CALL run_smeagol_emtrans(qs_env, last=.TRUE., iter=0)
            ! .. or save matrices related to bulk leads
            CALL run_smeagol_bulktrans(qs_env)
         END IF

         ! *** cleanup
         CALL scf_env_cleanup(scf_env)
         IF (dft_control%qs_control%cdft) THEN
            CALL cdft_control_cleanup(dft_control%qs_control%cdft_control)
         END IF

         IF (PRESENT(has_converged)) THEN
            has_converged = converged
         END IF
         IF (PRESENT(total_scf_steps)) THEN
            total_scf_steps = tsteps
         END IF

      END IF

   END SUBROUTINE scf

! **************************************************************************************************
!> \brief perform an scf loop
!> \param scf_env the scf_env where to perform the scf procedure
!> \param scf_control ...
!> \param qs_env the qs_env, the scf_env lives in
!> \param converged will be true / false if converged is reached
!> \param should_stop ...
!> \param total_scf_steps ...
!> \par History
!>      long history, see cvs and qs_scf module history
!>      02.2003 introduced scf_env [fawzi]
!>      09.2005 Frozen density approximation [TdK]
!>      06.2007 Check for SCF iteration count early [jgh]
!>      10.2019 switch_surf_dip [SGh]
!> \author Matthias Krack
!> \note
! **************************************************************************************************
   SUBROUTINE scf_env_do_scf(scf_env, scf_control, qs_env, converged, should_stop, total_scf_steps)

      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(OUT)                               :: converged, should_stop
      INTEGER, INTENT(OUT)                               :: total_scf_steps

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'scf_env_do_scf'
      INTEGER, PARAMETER :: max_ot_kp_subspace_refreshes = 3, max_ot_kp_total_refreshes = 4

      CHARACTER(LEN=default_string_length)               :: description, name
      INTEGER :: ext_master_id, handle, handle2, i_tmp, ic, ispin, iter_count, &
         ot_kp_subspace_refresh_count, ot_kp_subspace_refresh_iter_count, output_unit, &
         post_added_mos_auto_resize_steps, post_ot_kp_subspace_refresh_steps, &
         scf_energy_message_tag, total_steps
      LOGICAL :: added_mos_auto_grow, density_full_step, diis_step, do_kpoints, energy_only, &
         exit_inner_loop, exit_outer_loop, inner_loop_converged, &
         internal_tblite_density_full_step, internal_tblite_mixer, just_energy, &
         ot_kp_subspace_refresh, outer_loop_converged, settle_added_mos_frontier, &
         tblite_native_mixer
      REAL(KIND=dp)                                      :: t1, t2
      REAL(KIND=dp), DIMENSION(3)                        :: res_val_3
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_result_type), POINTER                      :: results
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rho_ao_kp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(energy_correction_type), POINTER              :: ec_env
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos, mos_last_converged
      TYPE(mp_comm_type)                                 :: external_comm
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(qs_charges_type), POINTER                     :: qs_charges
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(section_vals_type), POINTER                   :: dft_section, input, scf_section

! One rebuild seeds k-point OT from the first physical Hamiltonian; retain
! two later adaptive rebuilds for collapsed metallic subspaces and reserve
! one final canonical rebase for a calibrated orbital-dominated tail.

      CALL timeset(routineN, handle)

      NULLIFY (cell, dft_control, rho, energy, &
               logger, qs_charges, ks_env, mos, atomic_kind_set, qs_kind_set, &
               particle_set, dft_section, input, &
               scf_section, para_env, results, kpoints, pw_env, rho_ao_kp, mos_last_converged)

      CPASSERT(ASSOCIATED(scf_env))
      CPASSERT(ASSOCIATED(qs_env))

      logger => cp_get_default_logger()
      t1 = m_walltime()

      CALL get_qs_env(qs_env=qs_env, &
                      cell=cell, &
                      energy=energy, &
                      particle_set=particle_set, &
                      qs_charges=qs_charges, &
                      ks_env=ks_env, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      rho=rho, &
                      mos=mos, &
                      input=input, &
                      dft_control=dft_control, &
                      do_kpoints=do_kpoints, &
                      kpoints=kpoints, &
                      results=results, &
                      pw_env=pw_env, &
                      para_env=para_env)
      tblite_native_mixer = dft_control%qs_control%xtb_control%do_tblite .AND. &
                            scf_env%method /= ot_method_nr .AND. &
                            tb_native_scc_mixer_active(dft_control)
      internal_tblite_mixer = (dft_control%qs_control%dftb .AND. &
                               dft_control%qs_control%dftb_control%tblite_scc_mixer == tblite_scc_mixer_tblite) .OR. &
                              (dft_control%qs_control%xtb .AND. &
                               .NOT. dft_control%qs_control%xtb_control%do_tblite .AND. &
                               dft_control%qs_control%xtb_control%tblite_scc_mixer == tblite_scc_mixer_tblite)
      internal_tblite_density_full_step = dft_control%qs_control%xtb .AND. &
                                          .NOT. dft_control%qs_control%xtb_control%do_tblite .AND. &
                                          dft_control%qs_control%xtb_control%tblite_scc_mixer == tblite_scc_mixer_tblite

      CALL qs_rho_get(rho, rho_ao_kp=rho_ao_kp)

      dft_section => section_vals_get_subs_vals(input, "DFT")
      scf_section => section_vals_get_subs_vals(dft_section, "SCF")

      output_unit = cp_print_key_unit_nr(logger, scf_section, "PRINT%PROGRAM_RUN_INFO", &
                                         extension=".scfLog")

      IF (scf_control%gce%do_gce .AND. output_unit > 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T2,78('-'))")
         WRITE (UNIT=output_unit, FMT="(T31,A)") &
            "GRAND-CANONICAL SCF"
         WRITE (UNIT=output_unit, FMT="(T20,A,F12.6,A)") &
            "Target work function (TWF):", &
            evolt*scf_control%gce%target_workfunction, " eV"
         WRITE (UNIT=output_unit, FMT="(T2,78('-'))")
      END IF

      IF (output_unit > 0) WRITE (UNIT=output_unit, FMT="(/,/,T2,A)") &
         "SCF WAVEFUNCTION OPTIMIZATION"

      ! when switch_surf_dip is switched on, indicate storing mos from the last converged step
      IF (dft_control%switch_surf_dip) THEN
         CALL get_qs_env(qs_env, mos_last_converged=mos_last_converged)
         DO ispin = 1, dft_control%nspins
            CALL reassign_allocated_mos(mos(ispin), mos_last_converged(ispin))
         END DO
         IF (output_unit > 0) WRITE (UNIT=output_unit, FMT="(/,/,T2,A)") &
            "COPIED mos_last_converged ---> mos"
      END IF

      IF ((output_unit > 0) .AND. (.NOT. scf_control%use_ot)) THEN
         WRITE (UNIT=output_unit, &
                FMT="(/,T3,A,T12,A,T31,A,T39,A,T59,A,T75,A,/,T3,A)") &
            "Step", "Update method", "Time", "Convergence", "Total energy", "Change", &
            REPEAT("-", 78)
      END IF
      CALL cp_add_iter_level(logger%iter_info, "QS_SCF")

      ! check for external communicator and if the intermediate energy should be sent
      res_val_3(:) = -1.0_dp
      description = "[EXT_SCF_ENER_COMM]"
      IF (test_for_result(results, description=description)) THEN
         CALL get_results(results, description=description, &
                          values=res_val_3, n_entries=i_tmp)
         CPASSERT(i_tmp == 3)
         IF (ALL(res_val_3(:) <= 0.0)) THEN
            CALL cp_abort(__LOCATION__, &
                          " Trying to access result ("//TRIM(description)// &
                          ") which is not correctly stored.")
         END IF
         CALL external_comm%set_handle(NINT(res_val_3(1)))
      END IF
      ext_master_id = NINT(res_val_3(2))
      scf_energy_message_tag = NINT(res_val_3(3))

      ! *** outer loop of the scf, can treat other variables,
      ! *** such as lagrangian multipliers
      scf_env%outer_scf%iter_count = 0
      iter_count = 0
      post_added_mos_auto_resize_steps = 0
      ot_kp_subspace_refresh_count = 0
      ot_kp_subspace_refresh_iter_count = 0
      post_ot_kp_subspace_refresh_steps = 0
      settle_added_mos_frontier = .FALSE.
      total_steps = 0
      energy%tot_old = 0.0_dp

      scf_outer_loop: DO

         CALL init_scf_loop(scf_env=scf_env, qs_env=qs_env, &
                            scf_section=scf_section)

         CALL qs_scf_set_loop_flags(scf_env, diis_step, &
                                    energy_only, just_energy, exit_inner_loop)
         IF (post_ot_kp_subspace_refresh_steps > 0) THEN
            scf_env%iter_count = ot_kp_subspace_refresh_iter_count
         END IF

         ! decide whether to switch off dipole correction for convergence purposes
         dft_control%surf_dip_correct_switch = dft_control%correct_surf_dip
         IF ((dft_control%correct_surf_dip) .AND. (scf_control%outer_scf%have_scf) .AND. &
             (scf_env%outer_scf%iter_count > FLOOR(scf_control%outer_scf%max_scf/2.0_dp))) THEN
            IF (dft_control%switch_surf_dip) THEN
               dft_control%surf_dip_correct_switch = .FALSE.
               IF (output_unit > 0) WRITE (UNIT=output_unit, FMT="(/,/,T2,A)") &
                  "SURFACE DIPOLE CORRECTION switched off"
            END IF
         END IF

         scf_loop: DO

            CALL timeset(routineN//"_inner_loop", handle2)

            IF (.NOT. just_energy) scf_env%iter_count = scf_env%iter_count + 1
            iter_count = iter_count + 1
            CALL cp_iterate(logger%iter_info, last=.FALSE., iter_nr=iter_count)

            IF (output_unit > 0) CALL m_flush(output_unit)

            total_steps = total_steps + 1
            just_energy = energy_only

            CALL qs_ks_update_qs_env(qs_env, just_energy=just_energy, &
                                     calculate_forces=.FALSE.)

            ! print 'heavy weight' or relatively expensive quantities
            CALL qs_scf_loop_print(qs_env, scf_env, para_env)

            added_mos_auto_grow = .FALSE.
            ot_kp_subspace_refresh = .FALSE.
            IF (do_kpoints) THEN
               ! kpoints
               IF (dft_control%hairy_probes .EQV. .TRUE.) THEN
                  scf_control%smear%do_smear = .FALSE.
                  CALL qs_scf_new_mos_kp(qs_env, scf_env, scf_control, diis_step, dft_control%probe)
               ELSE
                  CALL qs_scf_new_mos_kp(qs_env, scf_env, scf_control, diis_step, &
                                         added_mos_auto_grow=added_mos_auto_grow, &
                                         ot_kp_subspace_refresh=ot_kp_subspace_refresh, &
                                         allow_ot_kp_subspace_refresh= &
                                         ot_kp_subspace_refresh_count < max_ot_kp_subspace_refreshes, &
                                         allow_ot_kp_orbital_tail_rebase= &
                                         ot_kp_subspace_refresh_count < max_ot_kp_total_refreshes, &
                                         settle_ot_kp_subspace= &
                                         post_ot_kp_subspace_refresh_steps > 0 .OR. &
                                         (settle_added_mos_frontier .AND. &
                                          post_added_mos_auto_resize_steps == &
                                          qs_scf_kp_added_mos_auto_post_resize_steps()))
               END IF
            ELSE
               ! Gamma points only
               IF (dft_control%hairy_probes .EQV. .TRUE.) THEN
                  scf_control%smear%do_smear = .FALSE.
                  CALL qs_scf_new_mos(qs_env, scf_env, scf_control, scf_section, diis_step, energy_only, &
                                      dft_control%probe, added_mos_auto_grow=added_mos_auto_grow)
               ELSE
                  CALL qs_scf_new_mos(qs_env, scf_env, scf_control, scf_section, diis_step, energy_only, &
                                      added_mos_auto_grow=added_mos_auto_grow)
               END IF
            END IF

            IF (ot_kp_subspace_refresh) THEN
               CALL qs_scf_refresh_ot_kp_subspace(qs_env, scf_env, output_unit)
               ot_kp_subspace_refresh_count = ot_kp_subspace_refresh_count + 1
               post_ot_kp_subspace_refresh_steps = 1
               total_steps = MAX(0, total_steps - 1)
               iter_count = MAX(0, iter_count - 1)
               IF (.NOT. just_energy) scf_env%iter_count = MAX(0, scf_env%iter_count - 1)
               ot_kp_subspace_refresh_iter_count = scf_env%iter_count
               CALL timestop(handle2)
               CYCLE scf_outer_loop
            END IF

            IF (post_ot_kp_subspace_refresh_steps > 0) THEN
               post_ot_kp_subspace_refresh_steps = post_ot_kp_subspace_refresh_steps - 1
            END IF

            IF (added_mos_auto_grow) THEN
               IF (do_kpoints) THEN
                  CALL qs_scf_grow_added_mos_auto_kp( &
                     qs_env, scf_env, scf_control, output_unit, &
                     reseed_global_basis=ot_kp_subspace_refresh_count > 0, &
                     settle_frontier=settle_added_mos_frontier)
               ELSE
                  CALL qs_scf_grow_added_mos_auto(qs_env, scf_env, scf_control, output_unit)
               END IF
               total_steps = MAX(0, total_steps - 1)
               iter_count = MAX(0, iter_count - 1)
               IF (.NOT. just_energy) scf_env%iter_count = MAX(0, scf_env%iter_count - 1)
               IF (do_kpoints) THEN
                  post_added_mos_auto_resize_steps = qs_scf_kp_added_mos_auto_post_resize_steps()
               ELSE
                  post_added_mos_auto_resize_steps = 1
               END IF
               CALL timestop(handle2)
               CYCLE scf_outer_loop
            END IF

            ! Print requested MO information (can be computationally expensive with OT)
            CALL qs_scf_write_mos(qs_env, scf_env, final_mos=.FALSE.)

            IF (dft_control%qs_control%xtb_control%do_tblite) THEN
               IF (scf_env%method == ot_method_nr) THEN
                  CALL tb_update_charges(qs_env, dft_control, qs_env%tb_tblite, .FALSE., .TRUE.)
                  CALL evaluate_core_matrix_traces(qs_env)
               ELSE
                  CPASSERT(scf_env%mixing_method > 0)
                  CALL tb_update_charges(qs_env, dft_control, qs_env%tb_tblite, .FALSE., .FALSE.)
                  CALL evaluate_core_matrix_traces(qs_env, rho_ao_ext=scf_env%p_mix_new)
               END IF
               CALL tb_get_energy(qs_env, qs_env%tb_tblite, energy)
            END IF

            density_full_step = diis_step .OR. tblite_native_mixer .OR. internal_tblite_density_full_step
            CALL qs_scf_density_mixing(scf_env, rho, para_env, density_full_step)
            IF (dft_control%qs_control%xtb_control%do_tblite .AND. &
                .NOT. (dft_control%qs_control%do_ls_scf .OR. scf_control%use_ot)) THEN
               scf_env%iter_delta = MAX(scf_env%iter_delta, &
                                        tb_scf_mixer_error(dft_control, qs_env%tb_tblite, &
                                                           scf_control%eps_scf))
            END IF
            IF (dft_control%qs_control%dftb .OR. &
                (dft_control%qs_control%xtb .AND. .NOT. dft_control%qs_control%xtb_control%do_tblite)) THEN
               scf_env%iter_delta = MAX(scf_env%iter_delta, &
                                        charge_mixing_scc_error(scf_env%mixing_store, scf_control%eps_scf))
            END IF
            ! K-point ADDED_MOS AUTO resizes and diagonalizes the virtual
            ! space against the current density.  With OT/no density mixing
            ! the orbital residual can be below EPS_SCF before the density
            ! has reacted in the enlarged MO space, so keep the convergence
            ! gate closed for the short AUTO probe window.
            IF (post_added_mos_auto_resize_steps > 0 .AND. .NOT. just_energy) THEN
               scf_env%iter_delta = MAX(scf_env%iter_delta, &
                                        scf_control%eps_scf*(1.0_dp + 10.0_dp*EPSILON(1.0_dp)))
               post_added_mos_auto_resize_steps = post_added_mos_auto_resize_steps - 1
            END IF
            IF (tblite_native_mixer) THEN
               scf_env%iter_param = dft_control%qs_control%xtb_control%tblite_mixer_damping
               scf_env%iter_method = "TBLite/Diag"
            ELSE IF (internal_tblite_mixer) THEN
               scf_env%iter_method = "TBLite/Diag"
               IF (dft_control%qs_control%dftb) THEN
                  scf_env%iter_param = dft_control%qs_control%dftb_control%tblite_mixer_damping
               ELSE
                  scf_env%iter_param = dft_control%qs_control%xtb_control%tblite_mixer_damping
               END IF
            END IF

            t2 = m_walltime()

            CALL qs_scf_loop_info(scf_env, output_unit, just_energy, t1, t2, energy)

            IF (scf_control%gce%do_gce) THEN
               CALL qs_scf_gce_info(output_unit, qs_env, just_energy)
            END IF

            IF (.NOT. just_energy) energy%tot_old = energy%total

            ! check for external communicator and if the intermediate energy should be sent
            IF (scf_energy_message_tag > 0) THEN
               CALL external_comm%send(energy%total, ext_master_id, scf_energy_message_tag)
            END IF

            CALL qs_scf_check_inner_exit(qs_env, scf_env, scf_control, should_stop, just_energy, &
                                         exit_inner_loop, inner_loop_converged, output_unit)

            ! In case we decide to exit we perform few more check to see if this one
            ! is really the last SCF step
            IF (exit_inner_loop) THEN

               CALL qs_scf_inner_finalize(scf_env, qs_env, density_full_step, output_unit)

               CALL qs_scf_check_outer_exit(qs_env, scf_env, scf_control, should_stop, &
                                            outer_loop_converged, exit_outer_loop)

               ! Let's tag the last SCF cycle so we can print informations only of the last step
               IF (exit_outer_loop) CALL cp_iterate(logger%iter_info, last=.TRUE., iter_nr=iter_count)

            END IF

            IF (do_kpoints) THEN
               CALL write_kpoints_restart(rho_ao_kp, kpoints, scf_env, dft_section, &
                                          particle_set, qs_kind_set, cell)
               IF (exit_inner_loop .AND. ASSOCIATED(scf_env%qs_ot_env)) THEN
                  IF (ANY(scf_env%qs_ot_env(:)%kpoint_index > 0)) THEN
                     DO ispin = 1, SIZE(scf_env%qs_ot_env)
                        CALL ot_scf_destroy(scf_env%qs_ot_env(ispin))
                     END DO
                     DEALLOCATE (scf_env%qs_ot_env)
                     NULLIFY (scf_env%qs_ot_env)
                  END IF
               END IF
            ELSE
               ! Write wavefunction restart file
               IF (scf_env%method == ot_method_nr) THEN
                  ! With OT: provide the Kohn-Sham matrix for the calculation of the MO eigenvalues
                  CALL get_ks_env(ks_env=ks_env, matrix_ks=matrix_ks)
                  CALL write_mo_set_to_restart(mos, particle_set, dft_section=dft_section, qs_kind_set=qs_kind_set, &
                                               matrix_ks=matrix_ks)
               ELSE
                  CALL write_mo_set_to_restart(mos, particle_set, dft_section=dft_section, qs_kind_set=qs_kind_set)
               END IF

            END IF

            ! Exit if we have finished with the SCF inner loop
            IF (exit_inner_loop) THEN
               CALL timestop(handle2)
               EXIT scf_loop
            END IF

            IF (.NOT. BTEST(cp_print_key_should_output(logger%iter_info, &
                                                       scf_section, "PRINT%ITERATION_INFO/TIME_CUMUL"), cp_p_file)) THEN
               t1 = m_walltime()
            END IF

            ! mixing methods have the new density matrix in p_mix_new
            IF (scf_env%mixing_method > 0) THEN
               DO ic = 1, SIZE(rho_ao_kp, 2)
                  DO ispin = 1, dft_control%nspins
                     CALL dbcsr_get_info(rho_ao_kp(ispin, ic)%matrix, name=name) ! keep the name
                     CALL dbcsr_copy(rho_ao_kp(ispin, ic)%matrix, scf_env%p_mix_new(ispin, ic)%matrix, name=name)
                  END DO
               END DO
            END IF

            CALL qs_scf_rho_update(rho, qs_env, scf_env, ks_env, &
                                   mix_rho=scf_env%mixing_method >= gspace_mixing_nr)

            CALL timestop(handle2)

         END DO scf_loop

         IF (.NOT. scf_control%outer_scf%have_scf) EXIT scf_outer_loop

         ! In case we use the OUTER SCF loop let's print some info..
         CALL qs_scf_outer_loop_info(output_unit, scf_control, scf_env, &
                                     energy, total_steps, should_stop, outer_loop_converged)

         ! Save MOs to converged MOs if outer_loop_converged and surf_dip_correct_switch is true
         IF (exit_outer_loop) THEN
            IF ((dft_control%switch_surf_dip) .AND. (outer_loop_converged) .AND. &
                (dft_control%surf_dip_correct_switch)) THEN
               DO ispin = 1, dft_control%nspins
                  CALL reassign_allocated_mos(mos_last_converged(ispin), mos(ispin))
               END DO
               IF (output_unit > 0) WRITE (UNIT=output_unit, FMT="(/,/,T2,A)") &
                  "COPIED mos ---> mos_last_converged"
            END IF
         END IF

         IF (exit_outer_loop) EXIT scf_outer_loop

         !
         CALL outer_loop_optimize(scf_env, scf_control)
         CALL outer_loop_update_qs_env(qs_env, scf_env)
         CALL qs_ks_did_change(ks_env, potential_changed=.TRUE.)

      END DO scf_outer_loop

      converged = inner_loop_converged .AND. outer_loop_converged
      total_scf_steps = total_steps

      IF (dft_control%qs_control%cdft) THEN
         dft_control%qs_control%cdft_control%total_steps = &
            dft_control%qs_control%cdft_control%total_steps + total_steps
      END IF

      IF (.NOT. converged) THEN
         IF (scf_control%ignore_convergence_failure .OR. should_stop) THEN
            CALL cp_warn(__LOCATION__, "SCF run NOT converged")
         ELSE
            CALL cp_abort(__LOCATION__, &
                          "SCF run NOT converged. To continue the calculation "// &
                          "regardless, please set the keyword IGNORE_CONVERGENCE_FAILURE.")
         END IF
      END IF

      ! Skip Harris functional calculation if ground-state is NOT converged
      IF (qs_env%energy_correction) THEN
         CALL get_qs_env(qs_env, ec_env=ec_env)
         ec_env%do_skip = .FALSE.
         IF (ec_env%skip_ec .AND. .NOT. converged) ec_env%do_skip = .TRUE.
      END IF

      ! if needed copy mo_coeff dbcsr->fm for later use in post_scf!fm->dbcsr
      DO ispin = 1, SIZE(mos) !fm -> dbcsr
         IF (mos(ispin)%use_mo_coeff_b) THEN !fm->dbcsr
            IF (.NOT. ASSOCIATED(mos(ispin)%mo_coeff_b)) THEN
               !fm->dbcsr
               CPABORT("mo_coeff_b is not allocated")
            END IF !fm->dbcsr
            CALL copy_dbcsr_to_fm(mos(ispin)%mo_coeff_b, & !fm->dbcsr
                                  mos(ispin)%mo_coeff) !fm -> dbcsr
         END IF !fm->dbcsr
      END DO !fm -> dbcsr

      CALL cp_rm_iter_level(logger%iter_info, level_name="QS_SCF")
      CALL timestop(handle)

   END SUBROUTINE scf_env_do_scf

! **************************************************************************************************
!> \brief grow the automatic smearing virtual space and reset OT objects
!> \param qs_env ...
!> \param scf_env ...
!> \param scf_control ...
!> \param output_unit ...
! **************************************************************************************************
   SUBROUTINE qs_scf_grow_added_mos_auto(qs_env, scf_env, scf_control, output_unit)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(scf_control_type), POINTER                    :: scf_control
      INTEGER, INTENT(IN)                                :: output_unit

      CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_scf_grow_added_mos_auto'

      INTEGER                                            :: base_nmo, current_added, grow_by, &
                                                            handle, ispin, j, nmo_mat, nspins, &
                                                            old_nmo
      INTEGER, DIMENSION(2)                              :: new_added, new_nmo
      LOGICAL                                            :: has_unit_metric, need_resize
      REAL(KIND=dp)                                      :: energy_step, flexible_electron_count, &
                                                            last_occ, maxocc, n_el_f, threshold
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvalues, occupation_numbers, &
                                                            old_eigenvalues, old_occupation_numbers
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_pool_p_type), DIMENSION(:), POINTER     :: ao_mo_fm_pools
      TYPE(cp_fm_type), POINTER                          :: mo_coeff, old_mo_coeff
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_type), ALLOCATABLE, DIMENSION(:)       :: old_mos
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      NULLIFY (ao_mo_fm_pools, blacs_env, dft_control, eigenvalues, matrix_s, mo_coeff, mos, &
               occupation_numbers, old_eigenvalues, old_mo_coeff, old_occupation_numbers, para_env)

      CALL get_qs_env(qs_env=qs_env, blacs_env=blacs_env, dft_control=dft_control, &
                      has_unit_metric=has_unit_metric, matrix_s=matrix_s, mos=mos, para_env=para_env)

      CPASSERT(ASSOCIATED(mos))
      CPASSERT(ASSOCIATED(dft_control))
      CPASSERT(ASSOCIATED(scf_control))
      CPASSERT(ASSOCIATED(matrix_s))

      nspins = dft_control%nspins
      new_added = scf_control%added_mos
      new_nmo(:) = 0
      need_resize = .FALSE.
      threshold = MAX(scf_control%smear%eps_fermi_dirac, 10.0_dp*EPSILON(1.0_dp))

      DO ispin = 1, nspins
         new_nmo(ispin) = mos(ispin)%nmo
         IF (.NOT. scf_control%added_mos_auto(ispin)) CYCLE
         IF (.NOT. ASSOCIATED(mos(ispin)%occupation_numbers)) CYCLE
         IF (mos(ispin)%nmo < 1) CYCLE

         last_occ = ABS(mos(ispin)%occupation_numbers(mos(ispin)%nmo))
         IF (last_occ <= threshold) CYCLE
         IF (mos(ispin)%nmo >= mos(ispin)%nao) THEN
            CALL cp_abort(__LOCATION__, &
                          "ADDED_MOS AUTO exhausted the AO basis but the highest MO is still occupied. "// &
                          "Use a larger basis or a lower electronic temperature.")
         END IF

         current_added = MAX(0, scf_control%added_mos(ispin))
         base_nmo = MAX(0, mos(ispin)%nmo - current_added)
         grow_by = MAX(4, current_added)
         new_added(ispin) = MIN(mos(ispin)%nao - base_nmo, current_added + grow_by)
         new_nmo(ispin) = base_nmo + new_added(ispin)
         IF (new_nmo(ispin) <= mos(ispin)%nmo) THEN
            CALL cp_abort(__LOCATION__, &
                          "ADDED_MOS AUTO cannot grow further although the highest MO is still occupied.")
         END IF
         need_resize = .TRUE.
      END DO

      IF (nspins == 2) THEN
         IF (new_nmo(2) > new_nmo(1)) THEN
            IF (.NOT. scf_control%added_mos_auto(1)) THEN
               CALL cp_abort(__LOCATION__, &
                             "ADDED_MOS AUTO needs to grow beta beyond alpha, but alpha uses explicit ADDED_MOS.")
            END IF
            base_nmo = MAX(0, mos(1)%nmo - MAX(0, scf_control%added_mos(1)))
            IF (new_nmo(2) > mos(1)%nao) THEN
               CALL cp_abort(__LOCATION__, &
                             "ADDED_MOS AUTO exhausted the alpha AO basis while matching the beta buffer.")
            END IF
            new_nmo(1) = new_nmo(2)
            new_added(1) = new_nmo(1) - base_nmo
            need_resize = .TRUE.
         END IF
      END IF

      IF (.NOT. need_resize) THEN
         CALL timestop(handle)
         RETURN
      END IF

      scf_control%added_mos(1:nspins) = new_added(1:nspins)
      IF (output_unit > 0) THEN
         IF (nspins == 2) THEN
            WRITE (UNIT=output_unit, FMT="(T2,A,2I5)") &
               "ADDED_MOS AUTO: growing virtual-space buffer to:", scf_control%added_mos(1:nspins)
         ELSE
            WRITE (UNIT=output_unit, FMT="(T2,A,I0)") &
               "ADDED_MOS AUTO: growing virtual-space buffer to: ", scf_control%added_mos(1)
         END IF
      END IF

      ALLOCATE (old_mos(nspins))
      DO ispin = 1, nspins
         CALL duplicate_mo_set(old_mos(ispin), mos(ispin))
      END DO

      DO ispin = 1, nspins
         CALL get_mo_set(old_mos(ispin), maxocc=maxocc, n_el_f=n_el_f, &
                         flexible_electron_count=flexible_electron_count)
         CALL deallocate_mo_set(mos(ispin))
         CALL allocate_mo_set(mo_set=mos(ispin), nao=old_mos(ispin)%nao, nmo=new_nmo(ispin), &
                              nelectron=old_mos(ispin)%nelectron, n_el_f=n_el_f, maxocc=maxocc, &
                              flexible_electron_count=flexible_electron_count)
         mos(ispin)%use_mo_coeff_b = old_mos(ispin)%use_mo_coeff_b
      END DO

      CALL mpools_rebuild_fm_pools(qs_env%mpools, mos=mos, blacs_env=blacs_env, para_env=para_env)
      CALL mpools_get(qs_env%mpools, ao_mo_fm_pools=ao_mo_fm_pools)

      DO ispin = 1, nspins
         CALL init_mo_set(mos(ispin), fm_pool=ao_mo_fm_pools(ispin)%pool, &
                          name="qs_env%mo"//TRIM(ADJUSTL(cp_to_string(ispin))))
         CALL get_mo_set(mos(ispin), mo_coeff=mo_coeff, eigenvalues=eigenvalues, &
                         occupation_numbers=occupation_numbers)
         CALL get_mo_set(old_mos(ispin), mo_coeff=old_mo_coeff, eigenvalues=old_eigenvalues, &
                         occupation_numbers=old_occupation_numbers)

         old_nmo = old_mos(ispin)%nmo
         CALL cp_fm_init_random(mo_coeff, mos(ispin)%nmo)
         CALL cp_fm_to_fm(old_mo_coeff, mo_coeff, old_nmo)
         IF (has_unit_metric) THEN
            CALL make_basis_simple(mo_coeff, mos(ispin)%nmo)
         ELSE
            CALL make_basis_sm(mo_coeff, mos(ispin)%nmo, matrix_s(1)%matrix)
         END IF

         eigenvalues(1:old_nmo) = old_eigenvalues(1:old_nmo)
         occupation_numbers(:) = 0.0_dp
         occupation_numbers(1:old_nmo) = old_occupation_numbers(1:old_nmo)
         IF (mos(ispin)%nmo > old_nmo) THEN
            energy_step = MAX(scf_control%smear%electronic_temperature, 1.0E-3_dp)
            DO j = old_nmo + 1, mos(ispin)%nmo
               eigenvalues(j) = eigenvalues(old_nmo) + energy_step*REAL(j - old_nmo, KIND=dp)
            END DO
         END IF
         mos(ispin)%homo = old_mos(ispin)%homo
         mos(ispin)%lfomo = old_mos(ispin)%lfomo
         mos(ispin)%kTS = old_mos(ispin)%kTS
         mos(ispin)%mu = old_mos(ispin)%mu
         mos(ispin)%uniform_occupation = old_mos(ispin)%uniform_occupation
      END DO

      IF (dft_control%restricted) CALL mo_set_restrict(mos)

      DO ispin = 1, nspins
         IF (.NOT. mos(ispin)%use_mo_coeff_b) CYCLE
         CALL get_mo_set(mos(ispin), mo_coeff=mo_coeff)
         CALL dbcsr_init_p(mos(ispin)%mo_coeff_b)
         CALL cp_dbcsr_m_by_n_from_row_template(mos(ispin)%mo_coeff_b, template=matrix_s(1)%matrix, &
                                                n=mos(ispin)%nmo, sym=dbcsr_type_no_symmetry)
         CALL copy_fm_to_dbcsr(mo_coeff, mos(ispin)%mo_coeff_b)
      END DO

      IF (ASSOCIATED(qs_env%mo_derivs)) THEN
         DO ispin = 1, SIZE(qs_env%mo_derivs)
            IF (ASSOCIATED(qs_env%mo_derivs(ispin)%matrix)) THEN
               CALL dbcsr_release_p(qs_env%mo_derivs(ispin)%matrix)
            END IF
         END DO
         DEALLOCATE (qs_env%mo_derivs)
         NULLIFY (qs_env%mo_derivs)
      END IF
      IF (qs_env%requires_mo_derivs) THEN
         nmo_mat = nspins
         IF (dft_control%restricted) nmo_mat = 1
         ALLOCATE (qs_env%mo_derivs(nmo_mat))
         DO ispin = 1, nmo_mat
            NULLIFY (qs_env%mo_derivs(ispin)%matrix)
            CALL dbcsr_init_p(qs_env%mo_derivs(ispin)%matrix)
            CALL dbcsr_create(qs_env%mo_derivs(ispin)%matrix, template=mos(ispin)%mo_coeff_b, &
                              name="mo_derivs", matrix_type=dbcsr_type_no_symmetry)
         END DO
      END IF

      IF (ASSOCIATED(scf_env%qs_ot_env)) THEN
         DO ispin = 1, SIZE(scf_env%qs_ot_env)
            CALL ot_scf_destroy(scf_env%qs_ot_env(ispin))
         END DO
         DEALLOCATE (scf_env%qs_ot_env)
         NULLIFY (scf_env%qs_ot_env)
      END IF
      reuse_precond = .FALSE.

      DO ispin = 1, nspins
         CALL deallocate_mo_set(old_mos(ispin))
      END DO
      DEALLOCATE (old_mos)

      CALL timestop(handle)

   END SUBROUTINE qs_scf_grow_added_mos_auto

! **************************************************************************************************
!> \brief discard a collapsed k-point OT subspace before rebuilding it by full diagonalization
!> \param qs_env ...
!> \param scf_env ...
!> \param output_unit ...
! **************************************************************************************************
   SUBROUTINE qs_scf_refresh_ot_kp_subspace(qs_env, scf_env, output_unit)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      INTEGER, INTENT(IN)                                :: output_unit

      INTEGER                                            :: ot_channel
      TYPE(kpoint_type), POINTER                         :: kpoints

      CPASSERT(ASSOCIATED(qs_env))
      CPASSERT(ASSOCIATED(scf_env))
      NULLIFY (kpoints)
      CALL get_qs_env(qs_env=qs_env, kpoints=kpoints)

      IF (output_unit > 0) THEN
         WRITE (UNIT=output_unit, FMT="(T2,A)") &
            "K-point OT: rebuilding physical virtual subspace by full KS diagonalization."
      END IF

      ! Capture the bounded scalar prior while the old channels still exist.
      ! Everything that depends on their orbital basis is destroyed below.
      CALL qs_scf_save_kpoint_ot_density_calibration(scf_env, subspace_refresh_prior=.TRUE.)

      ! The next outer-loop initialization performs do_general_diag_kp on the
      ! current Hamiltonian before allocating new channel-local OT state.  Old
      ! line-search/history and restart data refer to the collapsed subspace and
      ! must not be applied to the refreshed eigenvectors.
      IF (ASSOCIATED(scf_env%qs_ot_env)) THEN
         DO ot_channel = 1, SIZE(scf_env%qs_ot_env)
            CALL ot_scf_destroy(scf_env%qs_ot_env(ot_channel))
         END DO
         DEALLOCATE (scf_env%qs_ot_env)
         NULLIFY (scf_env%qs_ot_env)
      END IF
      IF (ASSOCIATED(kpoints)) THEN
         IF (ASSOCIATED(kpoints%scf_diis_buffer)) THEN
            CALL qs_diis_b_clear_kp(kpoints%scf_diis_buffer)
         END IF
      END IF
      CALL qs_scf_ot_restart_clear(scf_env)
      reuse_precond = .FALSE.

   END SUBROUTINE qs_scf_refresh_ot_kp_subspace

! **************************************************************************************************
!> \brief restore the k-point ADDED_MOS AUTO size serialized with an OT restart
!> \param qs_env ...
!> \param scf_env ...
!> \param scf_control ...
!> \param output_unit ...
! **************************************************************************************************
   SUBROUTINE qs_scf_restore_added_mos_auto_kp(qs_env, scf_env, scf_control, output_unit)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(scf_control_type), POINTER                    :: scf_control
      INTEGER, INTENT(IN)                                :: output_unit

      INTEGER                                            :: ichannel, ispin, nspins, restart_nmo
      INTEGER, DIMENSION(2)                              :: current_nmo, target_nmo
      LOGICAL                                            :: incompatible, need_resize
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos

      NULLIFY (dft_control, mos)
      IF (.NOT. scf_env%qs_ot_restart%valid) RETURN
      IF (.NOT. ASSOCIATED(scf_env%qs_ot_restart%channel)) RETURN

      CALL get_qs_env(qs_env=qs_env, dft_control=dft_control, mos=mos)
      CPASSERT(ASSOCIATED(dft_control))
      CPASSERT(ASSOCIATED(mos))
      nspins = dft_control%nspins
      CPASSERT(nspins >= 1 .AND. nspins <= SIZE(target_nmo))

      current_nmo(:) = 0
      target_nmo(:) = 0
      DO ispin = 1, nspins
         current_nmo(ispin) = mos(ispin)%nmo
      END DO

      incompatible = .FALSE.
      DO ichannel = 1, SIZE(scf_env%qs_ot_restart%channel)
         ASSOCIATE (channel => scf_env%qs_ot_restart%channel(ichannel))
            IF (.NOT. channel%valid) CYCLE
            ispin = channel%spin_index
            IF (ispin < 1 .OR. ispin > nspins .OR. channel%nmo <= 0) THEN
               incompatible = .TRUE.
               CYCLE
            END IF
            restart_nmo = channel%nmo
            IF (target_nmo(ispin) == 0) THEN
               target_nmo(ispin) = restart_nmo
            ELSE IF (target_nmo(ispin) /= restart_nmo) THEN
               incompatible = .TRUE.
            END IF
         END ASSOCIATE
      END DO
      IF (incompatible .OR. ANY(target_nmo(1:nspins) <= 0)) RETURN

      DO
         need_resize = .FALSE.
         DO ispin = 1, nspins
            current_nmo(ispin) = mos(ispin)%nmo
            IF (target_nmo(ispin) <= current_nmo(ispin)) CYCLE
            IF (.NOT. scf_control%added_mos_auto(ispin)) RETURN
            IF (target_nmo(ispin) > mos(ispin)%nao) THEN
               CALL cp_abort(__LOCATION__, &
                             "K-point OT restart requests more MOs than available in the AO basis.")
            END IF
            need_resize = .TRUE.
         END DO
         IF (.NOT. need_resize) EXIT

         CALL qs_scf_grow_added_mos_auto_kp(qs_env, scf_env, scf_control, output_unit, &
                                            restart_target_nmo=target_nmo(1:nspins))
      END DO

   END SUBROUTINE qs_scf_restore_added_mos_auto_kp

! **************************************************************************************************
!> \brief grow the automatic k-point smearing virtual space and reset k-point objects
!> \param qs_env ...
!> \param scf_env ...
!> \param scf_control ...
!> \param output_unit ...
!> \param restart_target_nmo exact MO dimensions requested by an AUTO restart
!> \param reseed_global_basis discard stale global placeholder orbitals after a k-point refresh
!> \param settle_frontier whether the rebuilt density must settle after restoring a thermal frontier
! **************************************************************************************************
   SUBROUTINE qs_scf_grow_added_mos_auto_kp(qs_env, scf_env, scf_control, output_unit, &
                                            restart_target_nmo, reseed_global_basis, settle_frontier)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(scf_control_type), POINTER                    :: scf_control
      INTEGER, INTENT(IN)                                :: output_unit
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL        :: restart_target_nmo
      LOGICAL, INTENT(IN), OPTIONAL                      :: reseed_global_basis
      LOGICAL, INTENT(OUT), OPTIONAL                     :: settle_frontier

      CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_scf_grow_added_mos_auto_kp'

      INTEGER                                            :: base_nmo, current_added, grow_by, &
                                                            handle, ic, ik, ispin, nmo_mat, &
                                                            nspins, old_nmo, target_nmo
      INTEGER, DIMENSION(2)                              :: new_added, new_nmo
      LOGICAL                                            :: diag_step, global_frontier_lost, &
                                                            has_unit_metric, need_resize, &
                                                            reseed_basis
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: restore_ot_ener_x
      REAL(KIND=dp)                                      :: energy_step, flexible_electron_count, &
                                                            maxocc, n_el_f
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: old_kpoint_kTS
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvalues, occupation_numbers, &
                                                            old_eigenvalues, old_occupation_numbers
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_pool_p_type), DIMENSION(:), POINTER     :: ao_mo_fm_pools
      TYPE(cp_fm_type), POINTER                          :: mo_coeff, old_mo_coeff
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mo_derivs
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks, matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_env_type), POINTER                     :: kp
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mo_set_type), ALLOCATABLE, DIMENSION(:)       :: old_mos
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)
      IF (PRESENT(settle_frontier)) settle_frontier = .FALSE.

      NULLIFY (ao_mo_fm_pools, blacs_env, dft_control, eigenvalues, kpoints, matrix_ks, matrix_s, &
               mo_coeff, mo_derivs, mos, occupation_numbers, old_eigenvalues, old_mo_coeff, &
               old_occupation_numbers, para_env)

      CALL get_qs_env(qs_env=qs_env, blacs_env=blacs_env, dft_control=dft_control, &
                      kpoints=kpoints, matrix_ks_kp=matrix_ks, matrix_s_kp=matrix_s, &
                      mo_derivs=mo_derivs, mos=mos, para_env=para_env, &
                      has_unit_metric=has_unit_metric)

      CPASSERT(ASSOCIATED(mos))
      CPASSERT(ASSOCIATED(kpoints))
      CPASSERT(ASSOCIATED(dft_control))
      CPASSERT(ASSOCIATED(scf_control))
      CPASSERT(ASSOCIATED(matrix_ks))
      CPASSERT(ASSOCIATED(matrix_s))

      nspins = dft_control%nspins
      reseed_basis = .FALSE.
      IF (PRESENT(reseed_global_basis)) reseed_basis = reseed_global_basis
      new_added = scf_control%added_mos
      new_nmo(:) = 0
      need_resize = .FALSE.

      DO ispin = 1, nspins
         old_nmo = mos(ispin)%nmo
         new_nmo(ispin) = old_nmo
         IF (.NOT. scf_control%added_mos_auto(ispin)) CYCLE

         current_added = MAX(0, scf_control%added_mos(ispin))
         base_nmo = MAX(0, old_nmo - current_added)
         IF (PRESENT(restart_target_nmo)) THEN
            CPASSERT(SIZE(restart_target_nmo) >= nspins)
            IF (restart_target_nmo(ispin) > mos(ispin)%nao) THEN
               CALL cp_abort(__LOCATION__, &
                             "K-point OT restart requests more MOs than available in the AO basis.")
            END IF
            grow_by = MAX(4, current_added)
            new_added(ispin) = MAX(current_added, &
                                   MIN(restart_target_nmo(ispin) - base_nmo, &
                                       current_added + grow_by))
         ELSE
            grow_by = MAX(4, current_added)
            new_added(ispin) = MIN(mos(ispin)%nao - base_nmo, current_added + grow_by)
         END IF
         new_nmo(ispin) = base_nmo + new_added(ispin)
         IF (new_nmo(ispin) > old_nmo) need_resize = .TRUE.
      END DO

      IF (nspins == 2) THEN
         target_nmo = MAX(new_nmo(1), new_nmo(2))
         DO ispin = 1, nspins
            IF (target_nmo > mos(ispin)%nao) THEN
               CALL cp_abort(__LOCATION__, &
                             "K-point ADDED_MOS AUTO exhausted the AO basis while matching spin bands.")
            END IF
            IF (target_nmo > mos(ispin)%nmo .AND. .NOT. scf_control%added_mos_auto(ispin)) THEN
               CALL cp_abort(__LOCATION__, &
                             "K-point ADDED_MOS AUTO needs to grow a spin channel with explicit ADDED_MOS.")
            END IF
            base_nmo = MAX(0, mos(ispin)%nmo - MAX(0, scf_control%added_mos(ispin)))
            new_nmo(ispin) = target_nmo
            new_added(ispin) = target_nmo - base_nmo
            IF (new_nmo(ispin) > mos(ispin)%nmo) need_resize = .TRUE.
         END DO
      END IF

      IF (.NOT. need_resize) THEN
         IF (PRESENT(restart_target_nmo)) THEN
            CALL timestop(handle)
            RETURN
         END IF
         CALL cp_abort(__LOCATION__, &
                       "K-point ADDED_MOS AUTO cannot grow further although the highest band is still occupied.")
      END IF

      scf_control%added_mos(1:nspins) = new_added(1:nspins)
      IF (output_unit > 0) THEN
         IF (PRESENT(restart_target_nmo) .AND. nspins == 2) THEN
            WRITE (UNIT=output_unit, FMT="(T2,A,2I5)") &
               "K-point OT restart: restoring ADDED_MOS AUTO virtual-space buffer to:", &
               scf_control%added_mos(1:nspins)
         ELSE IF (PRESENT(restart_target_nmo)) THEN
            WRITE (UNIT=output_unit, FMT="(T2,A,I0)") &
               "K-point OT restart: restoring ADDED_MOS AUTO virtual-space buffer to: ", &
               scf_control%added_mos(1)
         ELSE IF (nspins == 2) THEN
            WRITE (UNIT=output_unit, FMT="(T2,A,2I5)") &
               "K-point ADDED_MOS AUTO: growing virtual-space buffer to:", scf_control%added_mos(1:nspins)
         ELSE
            WRITE (UNIT=output_unit, FMT="(T2,A,I0)") &
               "K-point ADDED_MOS AUTO: growing virtual-space buffer to: ", scf_control%added_mos(1)
         END IF
      END IF

      ALLOCATE (old_mos(nspins))
      DO ispin = 1, nspins
         CALL duplicate_mo_set(old_mos(ispin), mos(ispin))
      END DO
      ALLOCATE (old_kpoint_kTS(nspins, SIZE(kpoints%kp_env)))
      DO ik = 1, SIZE(kpoints%kp_env)
         kp => kpoints%kp_env(ik)%kpoint_env
         CPASSERT(ASSOCIATED(kp))
         CPASSERT(ASSOCIATED(kp%mos))
         DO ispin = 1, nspins
            old_kpoint_kTS(ispin, ik) = kp%mos(1, ispin)%kTS
         END DO
      END DO
      DO ispin = 1, nspins
         CALL get_mo_set(old_mos(ispin), maxocc=maxocc, n_el_f=n_el_f, &
                         flexible_electron_count=flexible_electron_count)
         CALL deallocate_mo_set(mos(ispin))
         CALL allocate_mo_set(mo_set=mos(ispin), nao=old_mos(ispin)%nao, nmo=new_nmo(ispin), &
                              nelectron=old_mos(ispin)%nelectron, n_el_f=n_el_f, maxocc=maxocc, &
                              flexible_electron_count=flexible_electron_count)
         mos(ispin)%use_mo_coeff_b = old_mos(ispin)%use_mo_coeff_b
      END DO

      CALL mpools_rebuild_fm_pools(qs_env%mpools, mos=mos, blacs_env=blacs_env, para_env=para_env)
      CALL mpools_get(qs_env%mpools, ao_mo_fm_pools=ao_mo_fm_pools)

      DO ispin = 1, nspins
         CALL init_mo_set(mos(ispin), fm_pool=ao_mo_fm_pools(ispin)%pool, &
                          name="qs_env%mo"//TRIM(ADJUSTL(cp_to_string(ispin))))
         CALL get_mo_set(mos(ispin), mo_coeff=mo_coeff, eigenvalues=eigenvalues, &
                         occupation_numbers=occupation_numbers)
         CALL get_mo_set(old_mos(ispin), mo_coeff=old_mo_coeff, eigenvalues=old_eigenvalues, &
                         occupation_numbers=old_occupation_numbers)

         old_nmo = old_mos(ispin)%nmo
         CALL cp_fm_init_random(mo_coeff, mos(ispin)%nmo)
         ! Restart sizing and post-refresh AUTO growth can see global placeholder
         ! MOs that were not reconstructed from the current k-point orbitals.
         ! A full diagonalization follows immediately, so seed a complete
         ! independent basis instead of preserving rank-deficient columns.
         IF (.NOT. PRESENT(restart_target_nmo) .AND. .NOT. reseed_basis) THEN
            CALL cp_fm_to_fm(old_mo_coeff, mo_coeff, old_nmo)
         END IF
         IF (has_unit_metric) THEN
            CALL make_basis_simple(mo_coeff, mos(ispin)%nmo)
         ELSE
            CALL make_basis_sm(mo_coeff, mos(ispin)%nmo, matrix_s(1, 1)%matrix)
         END IF

         eigenvalues(1:old_nmo) = old_eigenvalues(1:old_nmo)
         occupation_numbers(:) = 0.0_dp
         occupation_numbers(1:old_nmo) = old_occupation_numbers(1:old_nmo)
         IF (mos(ispin)%nmo > old_nmo) THEN
            energy_step = MAX(scf_control%smear%electronic_temperature, 1.0E-3_dp)
            DO ic = old_nmo + 1, mos(ispin)%nmo
               eigenvalues(ic) = eigenvalues(old_nmo) + energy_step*REAL(ic - old_nmo, KIND=dp)
            END DO
         END IF
         mos(ispin)%homo = old_mos(ispin)%homo
         mos(ispin)%lfomo = old_mos(ispin)%lfomo
         mos(ispin)%kTS = old_mos(ispin)%kTS
         mos(ispin)%mu = old_mos(ispin)%mu
         mos(ispin)%uniform_occupation = old_mos(ispin)%uniform_occupation
      END DO

      IF (dft_control%restricted) CALL mo_set_restrict(mos)

      DO ispin = 1, nspins
         IF (.NOT. mos(ispin)%use_mo_coeff_b) CYCLE
         CALL get_mo_set(mos(ispin), mo_coeff=mo_coeff)
         CALL dbcsr_init_p(mos(ispin)%mo_coeff_b)
         CALL cp_dbcsr_m_by_n_from_row_template(mos(ispin)%mo_coeff_b, template=matrix_s(1, 1)%matrix, &
                                                n=mos(ispin)%nmo, sym=dbcsr_type_no_symmetry)
         CALL copy_fm_to_dbcsr(mo_coeff, mos(ispin)%mo_coeff_b)
      END DO

      IF (ASSOCIATED(mo_derivs)) THEN
         DO ispin = 1, SIZE(mo_derivs)
            IF (ASSOCIATED(mo_derivs(ispin)%matrix)) THEN
               CALL dbcsr_release_p(mo_derivs(ispin)%matrix)
            END IF
         END DO
         DEALLOCATE (mo_derivs)
         NULLIFY (mo_derivs)
         CALL set_qs_env(qs_env, mo_derivs=mo_derivs)
      END IF
      IF (qs_env%requires_mo_derivs) THEN
         nmo_mat = nspins
         IF (dft_control%restricted) nmo_mat = 1
         ALLOCATE (mo_derivs(nmo_mat))
         DO ispin = 1, nmo_mat
            NULLIFY (mo_derivs(ispin)%matrix)
            CALL dbcsr_init_p(mo_derivs(ispin)%matrix)
            CALL dbcsr_create(mo_derivs(ispin)%matrix, template=mos(ispin)%mo_coeff_b, &
                              name="mo_derivs", matrix_type=dbcsr_type_no_symmetry)
         END DO
         CALL set_qs_env(qs_env, mo_derivs=mo_derivs)
      END IF

      DO ik = 1, SIZE(kpoints%kp_env)
         kp => kpoints%kp_env(ik)%kpoint_env
         IF (ASSOCIATED(kp%mos)) THEN
            DO ispin = 1, SIZE(kp%mos, 2)
               DO ic = 1, SIZE(kp%mos, 1)
                  CALL deallocate_mo_set(kp%mos(ic, ispin))
               END DO
            END DO
            DEALLOCATE (kp%mos)
            NULLIFY (kp%mos)
         END IF
         CALL cp_fm_release(kp%pmat)
         CALL cp_fm_release(kp%wmat)
         CALL cp_fm_release(kp%ot_hmat)
         CALL cp_fm_release(kp%ot_smat)
      END DO
      CALL mpools_release(kpoints%mpools)
      CALL kpoint_initialize_mos(kpoints, mos)
      CALL kpoint_initialize_mo_set(kpoints)

      diag_step = .FALSE.
      CALL do_general_diag_kp(matrix_ks, matrix_s, kpoints, scf_env, scf_control, .FALSE., diag_step, &
                              store_ot_matrices=.TRUE.)
      ALLOCATE (restore_ot_ener_x(SIZE(scf_env%qs_ot_env)))
      CALL qs_scf_detect_kpoint_ot_frontier_loss( &
         kpoints, scf_env%qs_ot_env, old_kpoint_kTS, restore_ot_ener_x, global_frontier_lost)
      IF (PRESENT(settle_frontier)) settle_frontier = global_frontier_lost

      IF (ASSOCIATED(kpoints%scf_diis_buffer)) THEN
         CALL qs_diis_b_clear_kp(kpoints%scf_diis_buffer)
      END IF

      CALL qs_scf_save_kpoint_ot_density_calibration(scf_env, restore_ot_ener_x)
      DEALLOCATE (restore_ot_ener_x)

      IF (ASSOCIATED(scf_env%qs_ot_env)) THEN
         DO ispin = 1, SIZE(scf_env%qs_ot_env)
            CALL ot_scf_destroy(scf_env%qs_ot_env(ispin))
         END DO
         DEALLOCATE (scf_env%qs_ot_env)
         NULLIFY (scf_env%qs_ot_env)
      END IF
      reuse_precond = .FALSE.

      DO ispin = 1, nspins
         CALL deallocate_mo_set(old_mos(ispin))
      END DO
      DEALLOCATE (old_mos)
      DEALLOCATE (old_kpoint_kTS)

      CALL timestop(handle)

   END SUBROUTINE qs_scf_grow_added_mos_auto_kp

! **************************************************************************************************
!> \brief detect loss of a protected metallic frontier during a k-point AUTO resize
!>
!>        The resize diagonalization supplies the enlarged canonical orbital basis and the
!>        energies of newly added bands.  If it also removes established finite-temperature
!>        support, the accepted auxiliary-energy prefix must survive the second initialization
!>        diagonalization in the following outer SCF cycle.  Seed protected local channels before
!>        that diagonalization to retain band ordering and gauge, and persist all channel prefixes
!>        once a weighted thermal frontier has been established.  Request an additional
!>        density-settling step only when the MPI-reduced entropy also collapses globally.
!> \param kpoints resized k-point environment
!> \param qs_ot_env accepted pre-resize OT channels
!> \param old_kpoint_kTS entropy correction before the resize diagonalization
!> \param restore_ener_x channels whose accepted auxiliary-energy prefix must be restored
!> \param global_frontier_lost whether the weighted total entropy collapsed
! **************************************************************************************************
   SUBROUTINE qs_scf_detect_kpoint_ot_frontier_loss( &
      kpoints, qs_ot_env, old_kpoint_kTS, restore_ener_x, global_frontier_lost)
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: qs_ot_env
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: old_kpoint_kTS
      LOGICAL, DIMENSION(:), INTENT(OUT)                 :: restore_ener_x
      LOGICAL, INTENT(OUT), OPTIONAL                     :: global_frontier_lost

      INTEGER                                            :: ic, ikpoint_local, ispin, ncopy, nspin, &
                                                            ot_channel
      LOGICAL                                            :: global_collapse, lost_thermal_frontier, &
                                                            thermal_frontier_established
      REAL(KIND=dp)                                      :: new_weighted_kTS, old_weighted_kTS
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvalues
      TYPE(kpoint_env_type), POINTER                     :: kp

      NULLIFY (eigenvalues, kp)
      restore_ener_x(:) = .FALSE.
      IF (PRESENT(global_frontier_lost)) global_frontier_lost = .FALSE.
      IF (.NOT. ASSOCIATED(qs_ot_env)) RETURN
      CPASSERT(SIZE(restore_ener_x) == SIZE(qs_ot_env))
      CPASSERT(SIZE(old_kpoint_kTS, 2) == SIZE(kpoints%kp_env))

      old_weighted_kTS = 0.0_dp
      new_weighted_kTS = 0.0_dp
      DO ikpoint_local = 1, SIZE(kpoints%kp_env)
         kp => kpoints%kp_env(ikpoint_local)%kpoint_env
         CPASSERT(ASSOCIATED(kp))
         CPASSERT(ASSOCIATED(kp%mos))
         nspin = MIN(SIZE(old_kpoint_kTS, 1), SIZE(kp%mos, 2))
         DO ispin = 1, nspin
            old_weighted_kTS = old_weighted_kTS + &
                               kp%wkp*ABS(old_kpoint_kTS(ispin, ikpoint_local))
            new_weighted_kTS = new_weighted_kTS + &
                               kp%wkp*ABS(kp%mos(1, ispin)%kTS)
         END DO
      END DO
      CPASSERT(ASSOCIATED(kpoints%para_env_inter_kp))
      CALL kpoints%para_env_inter_kp%sum(old_weighted_kTS)
      CALL kpoints%para_env_inter_kp%sum(new_weighted_kTS)
      global_collapse = old_weighted_kTS > 1.0E-12_dp .AND. &
                        new_weighted_kTS < 0.25_dp*old_weighted_kTS
      thermal_frontier_established = old_weighted_kTS > 1.0E-12_dp
      IF (PRESENT(global_frontier_lost)) global_frontier_lost = global_collapse

      DO ikpoint_local = 1, SIZE(kpoints%kp_env)
         kp => kpoints%kp_env(ikpoint_local)%kpoint_env
         CPASSERT(ASSOCIATED(kp))
         CPASSERT(ASSOCIATED(kp%mos))
         DO ot_channel = 1, SIZE(qs_ot_env)
            IF (.NOT. qs_ot_env(ot_channel)%state_allocated) CYCLE
            IF (.NOT. ASSOCIATED(qs_ot_env(ot_channel)%ener_x)) CYCLE
            IF (qs_ot_env(ot_channel)%kpoint_index /= kp%nkpoint) CYCLE
            IF (qs_ot_env(ot_channel)%spin_index < 1 .OR. &
                qs_ot_env(ot_channel)%spin_index > SIZE(kp%mos, 2)) CYCLE
            lost_thermal_frontier = &
               ABS(old_kpoint_kTS(qs_ot_env(ot_channel)%spin_index, ikpoint_local)) > 1.0E-12_dp .AND. &
               ABS(kp%mos(1, qs_ot_env(ot_channel)%spin_index)%kTS) < &
               0.25_dp*ABS(old_kpoint_kTS(qs_ot_env(ot_channel)%spin_index, ikpoint_local))
            IF (lost_thermal_frontier .OR. thermal_frontier_established) THEN
               qs_ot_env(ot_channel)%ener_x_fermi_frontier_armed = .TRUE.
            END IF
            IF (.NOT. qs_ot_env(ot_channel)%ener_x_fermi_frontier_armed) CYCLE

            ! Seed the second diagonalization with the accepted ordering and gauge.
            ! A globally thermal reduced mesh is free to replace this seed with its
            ! new canonical spectrum; a global collapse is restored afterwards.
            DO ic = 1, SIZE(kp%mos, 1)
               CALL get_mo_set(kp%mos(ic, qs_ot_env(ot_channel)%spin_index), &
                               eigenvalues=eigenvalues)
               ncopy = MIN(SIZE(eigenvalues), SIZE(qs_ot_env(ot_channel)%ener_x))
               IF (ncopy > 0) THEN
                  eigenvalues(1:ncopy) = qs_ot_env(ot_channel)%ener_x(1:ncopy)
               END IF
            END DO
            restore_ener_x(ot_channel) = thermal_frontier_established
         END DO
      END DO

   END SUBROUTINE qs_scf_detect_kpoint_ot_frontier_loss

! **************************************************************************************************
!> \brief cache accepted response calibration and energy coordinates across k-point OT rebuilds
!> \param scf_env ...
!> \param restore_ener_x channels whose accepted auxiliary-energy prefix must be restored
!> \param subspace_refresh_prior retain only a bounded scalar prior after an orbital rebase
! **************************************************************************************************
   SUBROUTINE qs_scf_save_kpoint_ot_density_calibration( &
      scf_env, restore_ener_x, subspace_refresh_prior)
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      LOGICAL, DIMENSION(:), INTENT(IN), OPTIONAL        :: restore_ener_x
      LOGICAL, INTENT(IN), OPTIONAL                      :: subspace_refresh_prior

      INTEGER                                            :: ot_channel
      LOGICAL                                            :: refresh_prior

      refresh_prior = .FALSE.
      IF (PRESENT(subspace_refresh_prior)) refresh_prior = subspace_refresh_prior

      IF (ASSOCIATED(scf_env%qs_ot_schur_density_calibration%channel)) THEN
         DEALLOCATE (scf_env%qs_ot_schur_density_calibration%channel)
         NULLIFY (scf_env%qs_ot_schur_density_calibration%channel)
      END IF
      scf_env%qs_ot_schur_density_calibration%valid = .FALSE.
      scf_env%qs_ot_schur_density_calibration%nchannel = 0

      IF (.NOT. ASSOCIATED(scf_env%qs_ot_env)) RETURN
      IF (PRESENT(restore_ener_x)) THEN
         CPASSERT(SIZE(restore_ener_x) == SIZE(scf_env%qs_ot_env))
      END IF

      ALLOCATE (scf_env%qs_ot_schur_density_calibration%channel(SIZE(scf_env%qs_ot_env)))
      scf_env%qs_ot_schur_density_calibration%valid = .TRUE.
      scf_env%qs_ot_schur_density_calibration%nchannel = SIZE(scf_env%qs_ot_env)

      DO ot_channel = 1, SIZE(scf_env%qs_ot_env)
         ASSOCIATE (cache => scf_env%qs_ot_schur_density_calibration%channel(ot_channel), &
                    ot_env => scf_env%qs_ot_env(ot_channel))
            cache%subspace_refresh_prior = refresh_prior
            IF (refresh_prior) THEN
               cache%valid = ot_env%has_kpoint_context .AND. &
                             (ot_env%schur_coupled_quality_samples > 0 .OR. &
                              ot_env%mermin_ls_accepted_quality_samples > 0)
            ELSE
               cache%valid = ot_env%has_kpoint_context .AND. &
                             (ot_env%ener_x_fermi_frontier_armed .OR. &
                              ot_env%schur_projector_density_quality_samples > 0 .OR. &
                              ot_env%schur_projector_density_audit_samples > 0 .OR. &
                              ot_env%schur_projector_density_descent_quality_samples > 0 .OR. &
                              ot_env%schur_projector_physical_quality_samples > 0 .OR. &
                              ot_env%schur_projector_physical_selected_quality_samples > 0 .OR. &
                              ot_env%schur_projector_physical_direct_quality_samples > 0 .OR. &
                              ot_env%schur_projector_physical_direct_selected_quality_samples > 0 .OR. &
                              ot_env%schur_projector_phys_direct_density_sel_quality_samples > 0 .OR. &
                              ot_env%schur_projector_density_full_selected_quality_samples > 0 .OR. &
                              ot_env%schur_projector_density_full_direct_quality_samples > 0 .OR. &
                              ot_env%schur_projector_density_fdir_sel_quality_samples > 0 .OR. &
                              ot_env%schur_projector_density_full_global_quality_samples > 0 .OR. &
                              ot_env%schur_projector_density_full_global_sel_quality_samples > 0 .OR. &
                              ot_env%schur_density_stationarity_quality_samples > 0 .OR. &
                              ot_env%schur_projector_density_full_selected_stationarity_samples > 0 .OR. &
                              ot_env%mermin_ls_accepted_quality_samples > 0 .OR. &
                              ot_env%schur_projector_mixed_direct_quality_samples > 0)
            END IF
            IF (.NOT. cache%valid) CYCLE
            cache%spin_index = ot_env%spin_index
            cache%kpoint_index = ot_env%kpoint_index
            IF (refresh_prior) THEN
               cache%coupled_quality_samples = qs_scf_kp_rebased_calibration_samples( &
                                               ot_env%schur_coupled_quality_samples)
               cache%coupled_quality = ot_env%schur_coupled_quality
               cache%coupled_quality_ema = ot_env%schur_coupled_quality_ema
               cache%coupled_calibrated_quality = ot_env%schur_coupled_calibrated_quality
               cache%coupled_calibrated_quality_ema = &
                  ot_env%schur_coupled_calibrated_quality_ema
               cache%coupled_residual_quality = ot_env%schur_coupled_residual_quality
               cache%coupled_residual_quality_ema = ot_env%schur_coupled_residual_quality_ema
               cache%coupled_drop_sum = 0.0_dp
               cache%coupled_measured_sum = 0.0_dp
               cache%coupled_calibrated_drop_sum = 0.0_dp
               cache%mermin_ls_accepted_quality_samples = &
                  qs_scf_kp_rebased_calibration_samples( &
                  ot_env%mermin_ls_accepted_quality_samples)
               cache%mermin_ls_accepted_quality = ot_env%mermin_ls_accepted_quality
               cache%mermin_ls_accepted_quality_ema = ot_env%mermin_ls_accepted_quality_ema
               cache%mermin_ls_accepted_gap = ot_env%mermin_ls_accepted_gap
               cache%mermin_ls_accepted_gap_ema = ot_env%mermin_ls_accepted_gap_ema
               CYCLE
            END IF
            cache%ener_x_fermi_frontier_armed = ot_env%ener_x_fermi_frontier_armed
            cache%occupation_preconditioner_active = &
               ot_env%kpoint_occ_preconditioner_active
            cache%occupation_preconditioner_ready_samples = &
               ot_env%kpoint_occ_preconditioner_ready_samples
            IF (PRESENT(restore_ener_x)) THEN
               IF (restore_ener_x(ot_channel) .AND. ASSOCIATED(ot_env%ener_x)) THEN
                  ALLOCATE (cache%ener_x(SIZE(ot_env%ener_x)))
                  cache%ener_x(:) = ot_env%ener_x(:)
               END IF
            END IF
            cache%density_quality_samples = ot_env%schur_projector_density_quality_samples
            cache%density_audit_samples = ot_env%schur_projector_density_audit_samples
            cache%density_descent_quality_samples = &
               ot_env%schur_projector_density_descent_quality_samples
            cache%physical_quality_samples = ot_env%schur_projector_physical_quality_samples
            cache%physical_selected_quality_samples = &
               ot_env%schur_projector_physical_selected_quality_samples
            cache%physical_direct_quality_samples = &
               ot_env%schur_projector_physical_direct_quality_samples
            cache%physical_direct_selected_quality_samples = &
               ot_env%schur_projector_physical_direct_selected_quality_samples
            cache%physical_direct_density_sel_quality_samples = &
               ot_env%schur_projector_phys_direct_density_sel_quality_samples
            cache%mixed_direct_quality_samples = &
               ot_env%schur_projector_mixed_direct_quality_samples
            cache%density_full_selected_quality_samples = &
               ot_env%schur_projector_density_full_selected_quality_samples
            cache%density_stationarity_quality_samples = &
               ot_env%schur_density_stationarity_quality_samples
            cache%density_full_selected_stationarity_samples = &
               ot_env%schur_projector_density_full_selected_stationarity_samples
            cache%density_full_direct_quality_samples = &
               ot_env%schur_projector_density_full_direct_quality_samples
            cache%density_fdir_sel_quality_samples = &
               ot_env%schur_projector_density_fdir_sel_quality_samples
            cache%density_full_global_quality_samples = &
               ot_env%schur_projector_density_full_global_quality_samples
            cache%density_full_global_sel_quality_samples = &
               ot_env%schur_projector_density_full_global_sel_quality_samples
            cache%mermin_ls_accepted_quality_samples = &
               ot_env%mermin_ls_accepted_quality_samples
            cache%density_quality = ot_env%schur_projector_density_quality
            cache%density_quality_ema = ot_env%schur_projector_density_quality_ema
            cache%density_drop_sum = ot_env%schur_projector_density_drop_sum
            cache%density_measured_sum = ot_env%schur_projector_density_measured_sum
            cache%density_lowcurv_ema = ot_env%schur_projector_density_lowcurv_ema
            cache%density_lowcurv_quality_ema = &
               ot_env%schur_projector_density_lowcurv_quality_ema
            cache%density_audit_amplification_ema = &
               ot_env%schur_projector_density_audit_amplification_ema
            cache%density_audit_amplification_quality_ema = &
               ot_env%schur_projector_density_audit_amplification_quality_ema
            cache%density_audit_lowcurv_ema = &
               ot_env%schur_projector_density_audit_lowcurv_ema
            cache%density_audit_lowcurv_quality_ema = &
               ot_env%schur_projector_density_audit_lowcurv_quality_ema
            cache%density_descent_quality = ot_env%schur_projector_density_descent_quality
            cache%density_descent_quality_ema = &
               ot_env%schur_projector_density_descent_quality_ema
            cache%density_descent_drop_sum = &
               ot_env%schur_projector_density_descent_drop_sum
            cache%density_descent_measured_sum = &
               ot_env%schur_projector_density_descent_measured_sum
            cache%density_descent_coeff_ema = &
               ot_env%schur_projector_density_descent_coeff_ema
            cache%density_descent_support = ot_env%schur_projector_density_descent_support
            cache%density_descent_support_ema = &
               ot_env%schur_projector_density_descent_support_ema
            cache%density_descent_support_ready = &
               ot_env%schur_projector_density_descent_support_ready
            cache%density_full_selected_quality = &
               ot_env%schur_projector_density_full_selected_quality
            cache%density_full_selected_quality_ema = &
               ot_env%schur_projector_density_full_selected_quality_ema
            cache%density_full_selected_residual_quality = &
               ot_env%schur_projector_density_full_selected_residual_quality
            cache%density_full_selected_residual_quality_ema = &
               ot_env%schur_projector_density_full_selected_residual_quality_ema
            cache%density_stationarity_quality = &
               ot_env%schur_density_stationarity_quality
            cache%density_stationarity_quality_ema = &
               ot_env%schur_density_stationarity_quality_ema
            cache%density_full_selected_stationarity_quality = &
               ot_env%schur_projector_density_full_selected_stationarity_quality
            cache%density_full_selected_stationarity_quality_ema = &
               ot_env%schur_projector_density_full_selected_stationarity_quality_ema
            cache%density_full_selected_drop_sum = &
               ot_env%schur_projector_density_full_selected_drop_sum
            cache%density_full_selected_measured_sum = &
               ot_env%schur_projector_density_full_selected_measured_sum
            cache%density_full_direct_quality = &
               ot_env%schur_projector_density_full_direct_quality
            cache%density_full_direct_quality_ema = &
               ot_env%schur_projector_density_full_direct_quality_ema
            cache%density_full_direct_drop_sum = &
               ot_env%schur_projector_density_full_direct_drop_sum
            cache%density_full_direct_measured_sum = &
               ot_env%schur_projector_density_full_direct_measured_sum
            cache%density_full_direct_trust_ema = &
               ot_env%schur_projector_density_full_direct_trust_ema
            cache%density_full_direct_weighted_trust_ema = &
               ot_env%schur_projector_density_full_direct_weighted_trust_ema
            cache%density_full_direct_abs_trust_ema = &
               ot_env%schur_projector_density_full_direct_abs_trust_ema
            cache%density_fdir_sel_quality = &
               ot_env%schur_projector_density_fdir_sel_quality
            cache%density_fdir_sel_quality_ema = &
               ot_env%schur_projector_density_fdir_sel_quality_ema
            cache%density_fdir_sel_resid_quality = &
               ot_env%schur_projector_density_fdir_sel_resid_quality
            cache%density_fdir_sel_resid_quality_ema = &
               ot_env%schur_projector_density_fdir_sel_resid_quality_ema
            cache%density_fdir_sel_drop_sum = &
               ot_env%schur_projector_density_fdir_sel_drop_sum
            cache%density_fdir_sel_measured_sum = &
               ot_env%schur_projector_density_fdir_sel_measured_sum
            cache%density_full_global_quality = &
               ot_env%schur_projector_density_full_global_quality
            cache%density_full_global_quality_ema = &
               ot_env%schur_projector_density_full_global_quality_ema
            cache%density_full_global_drop_sum = &
               ot_env%schur_projector_density_full_global_drop_sum
            cache%density_full_global_measured_sum = &
               ot_env%schur_projector_density_full_global_measured_sum
            cache%density_full_global_trust_ema = &
               ot_env%schur_projector_density_full_global_trust_ema
            cache%density_full_global_weighted_trust_ema = &
               ot_env%schur_projector_density_full_global_weighted_trust_ema
            cache%density_full_global_abs_trust_ema = &
               ot_env%schur_projector_density_full_global_abs_trust_ema
            cache%density_full_global_selected = &
               ot_env%schur_projector_density_full_global_selected
            cache%density_full_global_sel_quality = &
               ot_env%schur_projector_density_full_global_sel_quality
            cache%density_full_global_sel_quality_ema = &
               ot_env%schur_projector_density_full_global_sel_quality_ema
            cache%density_full_global_sel_residual_quality = &
               ot_env%schur_projector_density_full_global_sel_resid_quality
            cache%density_full_global_sel_residual_quality_ema = &
               ot_env%schur_projector_density_full_global_sel_resid_quality_ema
            cache%density_full_global_sel_drop_sum = &
               ot_env%schur_projector_density_full_global_sel_drop_sum
            cache%density_full_global_sel_measured_sum = &
               ot_env%schur_projector_density_full_global_sel_measured_sum
            cache%physical_quality = ot_env%schur_projector_physical_quality
            cache%physical_quality_ema = ot_env%schur_projector_physical_quality_ema
            cache%physical_drop_sum = ot_env%schur_projector_physical_drop_sum
            cache%physical_measured_sum = ot_env%schur_projector_physical_measured_sum
            cache%physical_selected_quality = ot_env%schur_projector_physical_selected_quality
            cache%physical_selected_quality_ema = &
               ot_env%schur_projector_physical_selected_quality_line_ema
            cache%physical_selected_residual_quality = &
               ot_env%schur_projector_physical_selected_residual_quality
            cache%physical_selected_residual_quality_ema = &
               ot_env%schur_projector_physical_selected_residual_quality_ema
            cache%physical_selected_drop_sum = &
               ot_env%schur_projector_physical_selected_drop_sum
            cache%physical_selected_measured_sum = &
               ot_env%schur_projector_physical_selected_measured_sum
            cache%physical_trace_actual_ema = ot_env%schur_projector_physical_trace_actual_ema
            cache%physical_trace_raw_ema = ot_env%schur_projector_physical_trace_raw_ema
            cache%physical_trace_active_ema = ot_env%schur_projector_physical_trace_active_ema
            cache%physical_trace_actual_diag_ema = &
               ot_env%schur_projector_physical_trace_actual_diag_ema
            cache%physical_trace_actual_free_ema = &
               ot_env%schur_projector_physical_trace_actual_free_ema
            cache%physical_trace_raw_diag_ema = &
               ot_env%schur_projector_physical_trace_raw_diag_ema
            cache%physical_trace_raw_free_ema = &
               ot_env%schur_projector_physical_trace_raw_free_ema
            cache%physical_trace_active_diag_ema = &
               ot_env%schur_projector_physical_trace_active_diag_ema
            cache%physical_trace_active_free_ema = &
               ot_env%schur_projector_physical_trace_active_free_ema
            cache%physical_trace_invariant_ema = &
               ot_env%schur_projector_physical_trace_invariant_ema
            cache%physical_trace_invariant_diag_ema = &
               ot_env%schur_projector_physical_trace_invariant_diag_ema
            cache%physical_trace_invariant_free_ema = &
               ot_env%schur_projector_physical_trace_invariant_free_ema
            cache%physical_trace_actual_quality_ema = &
               ot_env%schur_projector_physical_trace_actual_quality_ema
            cache%physical_trace_raw_quality_ema = &
               ot_env%schur_projector_physical_trace_raw_quality_ema
            cache%physical_trace_active_quality_ema = &
               ot_env%schur_projector_physical_trace_active_quality_ema
            cache%physical_trace_actual_diag_quality_ema = &
               ot_env%schur_projector_physical_trace_actual_diag_quality_ema
            cache%physical_trace_actual_free_quality_ema = &
               ot_env%schur_projector_physical_trace_actual_free_quality_ema
            cache%physical_trace_raw_diag_quality_ema = &
               ot_env%schur_projector_physical_trace_raw_diag_quality_ema
            cache%physical_trace_raw_free_quality_ema = &
               ot_env%schur_projector_physical_trace_raw_free_quality_ema
            cache%physical_trace_active_diag_quality_ema = &
               ot_env%schur_projector_physical_trace_active_diag_quality_ema
            cache%physical_trace_active_free_quality_ema = &
               ot_env%schur_projector_physical_trace_active_free_quality_ema
            cache%physical_trace_invariant_quality_ema = &
               ot_env%schur_projector_physical_trace_invariant_quality_ema
            cache%physical_trace_invariant_diag_quality_ema = &
               ot_env%schur_projector_physical_trace_invariant_diag_quality_ema
            cache%physical_trace_invariant_free_quality_ema = &
               ot_env%schur_projector_physical_trace_invariant_free_quality_ema
            cache%physical_direct_quality = ot_env%schur_projector_physical_direct_quality
            cache%physical_direct_quality_shadow_ema = &
               ot_env%schur_projector_physical_direct_quality_shadow_ema
            cache%physical_direct_drop_sum = ot_env%schur_projector_physical_direct_drop_sum
            cache%physical_direct_measured_sum = &
               ot_env%schur_projector_physical_direct_measured_sum
            cache%physical_direct_selected_quality = &
               ot_env%schur_projector_physical_direct_selected_quality
            cache%physical_direct_selected_quality_ema = &
               ot_env%schur_projector_physical_direct_selected_quality_ema
            cache%physical_direct_selected_residual_quality_ema = &
               ot_env%schur_projector_physical_direct_selected_residual_quality_ema
            cache%physical_direct_selected_drop_sum = &
               ot_env%schur_projector_physical_direct_selected_drop_sum
            cache%physical_direct_selected_measured_sum = &
               ot_env%schur_projector_physical_direct_selected_measured_sum
            cache%physical_direct_density_sel_quality = &
               ot_env%schur_projector_phys_direct_density_sel_quality
            cache%physical_direct_density_sel_quality_ema = &
               ot_env%schur_projector_phys_direct_density_sel_quality_ema
            cache%physical_direct_density_sel_resid_quality = &
               ot_env%schur_projector_phys_direct_density_sel_resid_quality
            cache%physical_direct_density_sel_resid_quality_ema = &
               ot_env%schur_projector_phys_direct_density_sel_resid_quality_ema
            cache%physical_direct_density_sel_drop_sum = &
               ot_env%schur_projector_phys_direct_density_sel_drop_sum
            cache%physical_direct_density_sel_measured_sum = &
               ot_env%schur_projector_phys_direct_density_sel_measured_sum
            cache%mixed_direct_quality = ot_env%schur_projector_mixed_direct_quality
            cache%mixed_direct_quality_ema = ot_env%schur_projector_mixed_direct_quality_ema
            cache%mixed_direct_drop_sum = ot_env%schur_projector_mixed_direct_drop_sum
            cache%mixed_direct_measured_sum = ot_env%schur_projector_mixed_direct_measured_sum
            cache%mixed_direct_weight_ema = ot_env%schur_projector_mixed_direct_weight_ema
            cache%mermin_ls_accepted_quality = ot_env%mermin_ls_accepted_quality
            cache%mermin_ls_accepted_quality_ema = ot_env%mermin_ls_accepted_quality_ema
            cache%mermin_ls_accepted_gap = ot_env%mermin_ls_accepted_gap
            cache%mermin_ls_accepted_gap_ema = ot_env%mermin_ls_accepted_gap_ema
         END ASSOCIATE
      END DO

   END SUBROUTINE qs_scf_save_kpoint_ot_density_calibration

! **************************************************************************************************
!> \brief restore scalar density-response calibration after a k-point ADDED_MOS AUTO rebuild
!> \param qs_ot_env ...
!> \param scf_env ...
! **************************************************************************************************
   SUBROUTINE qs_scf_restore_kpoint_ot_density_calibration(qs_ot_env, scf_env)
      TYPE(qs_ot_type)                                   :: qs_ot_env
      TYPE(qs_scf_env_type), POINTER                     :: scf_env

      INTEGER                                            :: ot_channel

      IF (.NOT. scf_env%qs_ot_schur_density_calibration%valid) RETURN
      IF (.NOT. ASSOCIATED(scf_env%qs_ot_schur_density_calibration%channel)) RETURN
      IF (.NOT. qs_ot_env%has_kpoint_context) RETURN

      DO ot_channel = 1, SIZE(scf_env%qs_ot_schur_density_calibration%channel)
         ASSOCIATE (cache => scf_env%qs_ot_schur_density_calibration%channel(ot_channel))
            IF (.NOT. cache%valid) CYCLE
            IF (cache%spin_index /= qs_ot_env%spin_index) CYCLE
            IF (cache%kpoint_index /= qs_ot_env%kpoint_index) CYCLE

            IF (cache%subspace_refresh_prior) THEN
               qs_ot_env%schur_coupled_quality_samples = cache%coupled_quality_samples
               qs_ot_env%schur_coupled_quality = cache%coupled_quality
               qs_ot_env%schur_coupled_quality_ema = cache%coupled_quality_ema
               qs_ot_env%schur_coupled_calibrated_quality = cache%coupled_calibrated_quality
               qs_ot_env%schur_coupled_calibrated_quality_ema = &
                  cache%coupled_calibrated_quality_ema
               qs_ot_env%schur_coupled_residual_quality = cache%coupled_residual_quality
               qs_ot_env%schur_coupled_residual_quality_ema = &
                  cache%coupled_residual_quality_ema
               qs_ot_env%schur_coupled_drop_sum = 0.0_dp
               qs_ot_env%schur_coupled_measured_sum = 0.0_dp
               qs_ot_env%schur_coupled_calibrated_drop_sum = 0.0_dp
               qs_ot_env%mermin_ls_accepted_quality_samples = &
                  cache%mermin_ls_accepted_quality_samples
               qs_ot_env%mermin_ls_accepted_quality = cache%mermin_ls_accepted_quality
               qs_ot_env%mermin_ls_accepted_quality_ema = cache%mermin_ls_accepted_quality_ema
               qs_ot_env%mermin_ls_accepted_gap = cache%mermin_ls_accepted_gap
               qs_ot_env%mermin_ls_accepted_gap_ema = cache%mermin_ls_accepted_gap_ema
               RETURN
            END IF

            qs_ot_env%schur_projector_density_pending = .FALSE.
            qs_ot_env%schur_projector_density_audit_pending = .FALSE.
            qs_ot_env%schur_projector_density_descent_pending = .FALSE.
            qs_ot_env%schur_projector_density_full_pending = .FALSE.
            qs_ot_env%schur_projector_density_full_direct_pending = .FALSE.
            qs_ot_env%schur_projector_density_full_global_pending = .FALSE.
            qs_ot_env%schur_projector_physical_pending = .FALSE.
            qs_ot_env%schur_projector_density_drop = 0.0_dp
            qs_ot_env%schur_projector_density_lowcurv = 0.0_dp
            qs_ot_env%schur_projector_density_audit_amplification = 0.0_dp
            qs_ot_env%schur_projector_density_audit_lowcurv = 0.0_dp
            qs_ot_env%schur_projector_density_descent_drop = 0.0_dp
            qs_ot_env%schur_projector_density_descent_coeff = 0.0_dp
            qs_ot_env%schur_projector_density_descent_support = cache%density_descent_support
            qs_ot_env%schur_projector_density_full_drop = 0.0_dp
            qs_ot_env%schur_projector_density_full_selected = 0.0_dp
            qs_ot_env%schur_projector_density_full_trust = 1.0_dp
            qs_ot_env%schur_projector_density_full_weighted_trust = 0.0_dp
            qs_ot_env%schur_projector_density_full_abs_trust = 0.0_dp
            qs_ot_env%schur_projector_density_full_direct_drop = 0.0_dp
            qs_ot_env%schur_projector_density_fdir_selected = 0.0_dp
            qs_ot_env%schur_projector_density_full_direct_trust = 1.0_dp
            qs_ot_env%schur_projector_density_full_direct_weighted_trust = 0.0_dp
            qs_ot_env%schur_projector_density_full_direct_abs_trust = 0.0_dp
            qs_ot_env%schur_projector_density_full_global_drop = 0.0_dp
            qs_ot_env%schur_projector_density_full_global_trust = 1.0_dp
            qs_ot_env%schur_projector_density_full_global_weighted_trust = 0.0_dp
            qs_ot_env%schur_projector_density_full_global_abs_trust = 0.0_dp
            qs_ot_env%schur_projector_physical_drop = 0.0_dp
            qs_ot_env%schur_projector_physical_metric = 0.0_dp
            qs_ot_env%schur_projector_physical_weight = 0.0_dp
            qs_ot_env%schur_projector_physical_correction_drop = 0.0_dp
            qs_ot_env%schur_projector_physical_compensation = 0.0_dp
            qs_ot_env%schur_projector_physical_negative = 0.0_dp
            qs_ot_env%schur_projector_physical_direct = 0.0_dp
            qs_ot_env%schur_projector_physical_selected = 0.0_dp
            qs_ot_env%schur_projector_physical_shift = 0.0_dp
            qs_ot_env%schur_projector_physical_trace_actual = 0.0_dp
            qs_ot_env%schur_projector_physical_trace_raw = 0.0_dp
            qs_ot_env%schur_projector_physical_trace_active = 0.0_dp
            qs_ot_env%schur_projector_physical_trace_actual_diag = 0.0_dp
            qs_ot_env%schur_projector_physical_trace_actual_free = 0.0_dp
            qs_ot_env%schur_projector_physical_trace_raw_diag = 0.0_dp
            qs_ot_env%schur_projector_physical_trace_raw_free = 0.0_dp
            qs_ot_env%schur_projector_physical_trace_active_diag = 0.0_dp
            qs_ot_env%schur_projector_physical_trace_active_free = 0.0_dp
            qs_ot_env%schur_projector_physical_trace_invariant = 0.0_dp
            qs_ot_env%schur_projector_physical_trace_invariant_diag = 0.0_dp
            qs_ot_env%schur_projector_physical_trace_invariant_free = 0.0_dp
            qs_ot_env%schur_projector_mixed_direct_drop = 0.0_dp
            qs_ot_env%schur_projector_mixed_direct_weight = 0.0_dp

            qs_ot_env%schur_projector_density_quality_samples = cache%density_quality_samples
            qs_ot_env%schur_projector_density_audit_samples = cache%density_audit_samples
            qs_ot_env%schur_projector_density_descent_quality_samples = &
               cache%density_descent_quality_samples
            qs_ot_env%schur_projector_physical_quality_samples = cache%physical_quality_samples
            qs_ot_env%schur_projector_physical_selected_quality_samples = &
               cache%physical_selected_quality_samples
            qs_ot_env%schur_projector_physical_direct_quality_samples = &
               cache%physical_direct_quality_samples
            qs_ot_env%schur_projector_physical_direct_selected_quality_samples = &
               cache%physical_direct_selected_quality_samples
            qs_ot_env%schur_projector_phys_direct_density_sel_quality_samples = &
               cache%physical_direct_density_sel_quality_samples
            qs_ot_env%schur_projector_mixed_direct_quality_samples = &
               cache%mixed_direct_quality_samples
            qs_ot_env%schur_projector_density_full_selected_quality_samples = &
               cache%density_full_selected_quality_samples
            qs_ot_env%schur_density_stationarity_quality_samples = &
               cache%density_stationarity_quality_samples
            qs_ot_env%schur_projector_density_full_selected_stationarity_samples = &
               cache%density_full_selected_stationarity_samples
            qs_ot_env%schur_projector_density_full_direct_quality_samples = &
               cache%density_full_direct_quality_samples
            qs_ot_env%schur_projector_density_fdir_sel_quality_samples = &
               cache%density_fdir_sel_quality_samples
            qs_ot_env%schur_projector_density_full_global_quality_samples = &
               cache%density_full_global_quality_samples
            qs_ot_env%schur_projector_density_full_global_sel_quality_samples = &
               cache%density_full_global_sel_quality_samples
            qs_ot_env%mermin_ls_accepted_quality_samples = &
               cache%mermin_ls_accepted_quality_samples
            qs_ot_env%schur_projector_density_quality = cache%density_quality
            qs_ot_env%schur_projector_density_quality_ema = cache%density_quality_ema
            qs_ot_env%schur_projector_density_drop_sum = cache%density_drop_sum
            qs_ot_env%schur_projector_density_measured_sum = cache%density_measured_sum
            qs_ot_env%schur_projector_density_lowcurv_ema = cache%density_lowcurv_ema
            qs_ot_env%schur_projector_density_lowcurv_quality_ema = &
               cache%density_lowcurv_quality_ema
            qs_ot_env%schur_projector_density_audit_amplification_ema = &
               cache%density_audit_amplification_ema
            qs_ot_env%schur_projector_density_audit_amplification_quality_ema = &
               cache%density_audit_amplification_quality_ema
            qs_ot_env%schur_projector_density_audit_lowcurv_ema = &
               cache%density_audit_lowcurv_ema
            qs_ot_env%schur_projector_density_audit_lowcurv_quality_ema = &
               cache%density_audit_lowcurv_quality_ema
            qs_ot_env%schur_projector_density_descent_quality = cache%density_descent_quality
            qs_ot_env%schur_projector_density_descent_quality_ema = &
               cache%density_descent_quality_ema
            qs_ot_env%schur_projector_density_descent_drop_sum = cache%density_descent_drop_sum
            qs_ot_env%schur_projector_density_descent_measured_sum = &
               cache%density_descent_measured_sum
            qs_ot_env%schur_projector_density_descent_coeff_ema = &
               cache%density_descent_coeff_ema
            qs_ot_env%schur_projector_density_descent_support_ema = &
               cache%density_descent_support_ema
            qs_ot_env%schur_projector_density_descent_support_ready = &
               cache%density_descent_support_ready
            qs_ot_env%schur_projector_density_full_selected_quality = &
               cache%density_full_selected_quality
            qs_ot_env%schur_projector_density_full_selected_quality_ema = &
               cache%density_full_selected_quality_ema
            qs_ot_env%schur_projector_density_full_selected_residual_quality = &
               cache%density_full_selected_residual_quality
            qs_ot_env%schur_projector_density_full_selected_residual_quality_ema = &
               cache%density_full_selected_residual_quality_ema
            qs_ot_env%schur_density_stationarity_quality = &
               cache%density_stationarity_quality
            qs_ot_env%schur_density_stationarity_quality_ema = &
               cache%density_stationarity_quality_ema
            qs_ot_env%schur_projector_density_full_selected_stationarity_quality = &
               cache%density_full_selected_stationarity_quality
            qs_ot_env%schur_projector_density_full_selected_stationarity_quality_ema = &
               cache%density_full_selected_stationarity_quality_ema
            qs_ot_env%schur_projector_density_full_selected_drop_sum = &
               cache%density_full_selected_drop_sum
            qs_ot_env%schur_projector_density_full_selected_measured_sum = &
               cache%density_full_selected_measured_sum
            qs_ot_env%schur_projector_density_full_direct_quality = &
               cache%density_full_direct_quality
            qs_ot_env%schur_projector_density_full_direct_quality_ema = &
               cache%density_full_direct_quality_ema
            qs_ot_env%schur_projector_density_full_direct_drop_sum = &
               cache%density_full_direct_drop_sum
            qs_ot_env%schur_projector_density_full_direct_measured_sum = &
               cache%density_full_direct_measured_sum
            qs_ot_env%schur_projector_density_full_direct_trust_ema = &
               cache%density_full_direct_trust_ema
            qs_ot_env%schur_projector_density_full_direct_weighted_trust_ema = &
               cache%density_full_direct_weighted_trust_ema
            qs_ot_env%schur_projector_density_full_direct_abs_trust_ema = &
               cache%density_full_direct_abs_trust_ema
            qs_ot_env%schur_projector_density_fdir_sel_quality = &
               cache%density_fdir_sel_quality
            qs_ot_env%schur_projector_density_fdir_sel_quality_ema = &
               cache%density_fdir_sel_quality_ema
            qs_ot_env%schur_projector_density_fdir_sel_resid_quality = &
               cache%density_fdir_sel_resid_quality
            qs_ot_env%schur_projector_density_fdir_sel_resid_quality_ema = &
               cache%density_fdir_sel_resid_quality_ema
            qs_ot_env%schur_projector_density_fdir_sel_drop_sum = &
               cache%density_fdir_sel_drop_sum
            qs_ot_env%schur_projector_density_fdir_sel_measured_sum = &
               cache%density_fdir_sel_measured_sum
            qs_ot_env%schur_projector_density_full_global_quality = &
               cache%density_full_global_quality
            qs_ot_env%schur_projector_density_full_global_quality_ema = &
               cache%density_full_global_quality_ema
            qs_ot_env%schur_projector_density_full_global_drop_sum = &
               cache%density_full_global_drop_sum
            qs_ot_env%schur_projector_density_full_global_measured_sum = &
               cache%density_full_global_measured_sum
            qs_ot_env%schur_projector_density_full_global_trust_ema = &
               cache%density_full_global_trust_ema
            qs_ot_env%schur_projector_density_full_global_weighted_trust_ema = &
               cache%density_full_global_weighted_trust_ema
            qs_ot_env%schur_projector_density_full_global_abs_trust_ema = &
               cache%density_full_global_abs_trust_ema
            qs_ot_env%schur_projector_density_full_global_selected = &
               cache%density_full_global_selected
            qs_ot_env%schur_projector_density_full_global_sel_quality = &
               cache%density_full_global_sel_quality
            qs_ot_env%schur_projector_density_full_global_sel_quality_ema = &
               cache%density_full_global_sel_quality_ema
            qs_ot_env%schur_projector_density_full_global_sel_resid_quality = &
               cache%density_full_global_sel_residual_quality
            qs_ot_env%schur_projector_density_full_global_sel_resid_quality_ema = &
               cache%density_full_global_sel_residual_quality_ema
            qs_ot_env%schur_projector_density_full_global_sel_drop_sum = &
               cache%density_full_global_sel_drop_sum
            qs_ot_env%schur_projector_density_full_global_sel_measured_sum = &
               cache%density_full_global_sel_measured_sum
            qs_ot_env%schur_projector_physical_quality = cache%physical_quality
            qs_ot_env%schur_projector_physical_quality_ema = cache%physical_quality_ema
            qs_ot_env%schur_projector_physical_drop_sum = cache%physical_drop_sum
            qs_ot_env%schur_projector_physical_measured_sum = cache%physical_measured_sum
            qs_ot_env%schur_projector_physical_selected_quality = cache%physical_selected_quality
            qs_ot_env%schur_projector_physical_selected_quality_line_ema = &
               cache%physical_selected_quality_ema
            qs_ot_env%schur_projector_physical_selected_residual_quality = &
               cache%physical_selected_residual_quality
            qs_ot_env%schur_projector_physical_selected_residual_quality_ema = &
               cache%physical_selected_residual_quality_ema
            qs_ot_env%schur_projector_physical_selected_drop_sum = &
               cache%physical_selected_drop_sum
            qs_ot_env%schur_projector_physical_selected_measured_sum = &
               cache%physical_selected_measured_sum
            qs_ot_env%schur_projector_physical_trace_actual_ema = cache%physical_trace_actual_ema
            qs_ot_env%schur_projector_physical_trace_raw_ema = cache%physical_trace_raw_ema
            qs_ot_env%schur_projector_physical_trace_active_ema = cache%physical_trace_active_ema
            qs_ot_env%schur_projector_physical_trace_actual_diag_ema = &
               cache%physical_trace_actual_diag_ema
            qs_ot_env%schur_projector_physical_trace_actual_free_ema = &
               cache%physical_trace_actual_free_ema
            qs_ot_env%schur_projector_physical_trace_raw_diag_ema = &
               cache%physical_trace_raw_diag_ema
            qs_ot_env%schur_projector_physical_trace_raw_free_ema = &
               cache%physical_trace_raw_free_ema
            qs_ot_env%schur_projector_physical_trace_active_diag_ema = &
               cache%physical_trace_active_diag_ema
            qs_ot_env%schur_projector_physical_trace_active_free_ema = &
               cache%physical_trace_active_free_ema
            qs_ot_env%schur_projector_physical_trace_invariant_ema = &
               cache%physical_trace_invariant_ema
            qs_ot_env%schur_projector_physical_trace_invariant_diag_ema = &
               cache%physical_trace_invariant_diag_ema
            qs_ot_env%schur_projector_physical_trace_invariant_free_ema = &
               cache%physical_trace_invariant_free_ema
            qs_ot_env%schur_projector_physical_trace_actual_quality_ema = &
               cache%physical_trace_actual_quality_ema
            qs_ot_env%schur_projector_physical_trace_raw_quality_ema = &
               cache%physical_trace_raw_quality_ema
            qs_ot_env%schur_projector_physical_trace_active_quality_ema = &
               cache%physical_trace_active_quality_ema
            qs_ot_env%schur_projector_physical_trace_actual_diag_quality_ema = &
               cache%physical_trace_actual_diag_quality_ema
            qs_ot_env%schur_projector_physical_trace_actual_free_quality_ema = &
               cache%physical_trace_actual_free_quality_ema
            qs_ot_env%schur_projector_physical_trace_raw_diag_quality_ema = &
               cache%physical_trace_raw_diag_quality_ema
            qs_ot_env%schur_projector_physical_trace_raw_free_quality_ema = &
               cache%physical_trace_raw_free_quality_ema
            qs_ot_env%schur_projector_physical_trace_active_diag_quality_ema = &
               cache%physical_trace_active_diag_quality_ema
            qs_ot_env%schur_projector_physical_trace_active_free_quality_ema = &
               cache%physical_trace_active_free_quality_ema
            qs_ot_env%schur_projector_physical_trace_invariant_quality_ema = &
               cache%physical_trace_invariant_quality_ema
            qs_ot_env%schur_projector_physical_trace_invariant_diag_quality_ema = &
               cache%physical_trace_invariant_diag_quality_ema
            qs_ot_env%schur_projector_physical_trace_invariant_free_quality_ema = &
               cache%physical_trace_invariant_free_quality_ema
            qs_ot_env%schur_projector_physical_direct_quality = cache%physical_direct_quality
            qs_ot_env%schur_projector_physical_direct_quality_shadow_ema = &
               cache%physical_direct_quality_shadow_ema
            qs_ot_env%schur_projector_physical_direct_drop_sum = cache%physical_direct_drop_sum
            qs_ot_env%schur_projector_physical_direct_measured_sum = &
               cache%physical_direct_measured_sum
            qs_ot_env%schur_projector_physical_direct_selected_quality = &
               cache%physical_direct_selected_quality
            qs_ot_env%schur_projector_physical_direct_selected_quality_ema = &
               cache%physical_direct_selected_quality_ema
            qs_ot_env%schur_projector_physical_direct_selected_residual_quality_ema = &
               cache%physical_direct_selected_residual_quality_ema
            qs_ot_env%schur_projector_physical_direct_selected_drop_sum = &
               cache%physical_direct_selected_drop_sum
            qs_ot_env%schur_projector_physical_direct_selected_measured_sum = &
               cache%physical_direct_selected_measured_sum
            qs_ot_env%schur_projector_phys_direct_density_sel_quality = &
               cache%physical_direct_density_sel_quality
            qs_ot_env%schur_projector_phys_direct_density_sel_quality_ema = &
               cache%physical_direct_density_sel_quality_ema
            qs_ot_env%schur_projector_phys_direct_density_sel_resid_quality = &
               cache%physical_direct_density_sel_resid_quality
            qs_ot_env%schur_projector_phys_direct_density_sel_resid_quality_ema = &
               cache%physical_direct_density_sel_resid_quality_ema
            qs_ot_env%schur_projector_phys_direct_density_sel_drop_sum = &
               cache%physical_direct_density_sel_drop_sum
            qs_ot_env%schur_projector_phys_direct_density_sel_measured_sum = &
               cache%physical_direct_density_sel_measured_sum
            qs_ot_env%schur_projector_mixed_direct_quality = cache%mixed_direct_quality
            qs_ot_env%schur_projector_mixed_direct_quality_ema = cache%mixed_direct_quality_ema
            qs_ot_env%schur_projector_mixed_direct_drop_sum = cache%mixed_direct_drop_sum
            qs_ot_env%schur_projector_mixed_direct_measured_sum = cache%mixed_direct_measured_sum
            qs_ot_env%schur_projector_mixed_direct_weight_ema = cache%mixed_direct_weight_ema
            qs_ot_env%mermin_ls_accepted_quality = cache%mermin_ls_accepted_quality
            qs_ot_env%mermin_ls_accepted_quality_ema = cache%mermin_ls_accepted_quality_ema
            qs_ot_env%mermin_ls_accepted_gap = cache%mermin_ls_accepted_gap
            qs_ot_env%mermin_ls_accepted_gap_ema = cache%mermin_ls_accepted_gap_ema
            qs_ot_env%ener_x_fermi_frontier_armed = cache%ener_x_fermi_frontier_armed
            IF (cache%ener_x_fermi_frontier_armed) THEN
               qs_ot_env%kpoint_occ_preconditioner_active = &
                  cache%occupation_preconditioner_active
               qs_ot_env%kpoint_occ_preconditioner_ready_samples = &
                  cache%occupation_preconditioner_ready_samples
            END IF
            EXIT
         END ASSOCIATE
      END DO

   END SUBROUTINE qs_scf_restore_kpoint_ot_density_calibration

! **************************************************************************************************
!> \brief restore a protected auxiliary-energy prefix in the enlarged canonical energy gauge
!> \param qs_ot_env allocated k-point OT channel
!> \param scf_env SCF environment holding the transient AUTO-resize cache
! **************************************************************************************************
   SUBROUTINE qs_scf_restore_kpoint_ot_frontier_energies(qs_ot_env, scf_env)
      TYPE(qs_ot_type)                                   :: qs_ot_env
      TYPE(qs_scf_env_type), POINTER                     :: scf_env

      INTEGER                                            :: ncopy, ot_channel
      REAL(KIND=dp)                                      :: prefix_shift

      IF (.NOT. scf_env%qs_ot_schur_density_calibration%valid) RETURN
      IF (.NOT. ASSOCIATED(scf_env%qs_ot_schur_density_calibration%channel)) RETURN
      IF (.NOT. ASSOCIATED(qs_ot_env%ener_x)) RETURN
      IF (.NOT. qs_ot_env%has_kpoint_context) RETURN

      DO ot_channel = 1, SIZE(scf_env%qs_ot_schur_density_calibration%channel)
         ASSOCIATE (cache => scf_env%qs_ot_schur_density_calibration%channel(ot_channel))
            IF (.NOT. cache%valid) CYCLE
            IF (cache%spin_index /= qs_ot_env%spin_index) CYCLE
            IF (cache%kpoint_index /= qs_ot_env%kpoint_index) CYCLE
            IF (.NOT. ALLOCATED(cache%ener_x)) RETURN
            ncopy = MIN(SIZE(cache%ener_x), SIZE(qs_ot_env%ener_x))
            IF (ncopy > 0) THEN
               ! A uniform shift of all auxiliary energies is a fixed-N gauge.
               ! Preserve the accepted prefix spacings, but anchor that prefix
               ! to the freshly diagonalized canonical tail so the reported
               ! chemical potential retains the physical KS energy reference.
               prefix_shift = qs_ot_env%ener_x(ncopy) - cache%ener_x(ncopy)
               qs_ot_env%ener_x(1:ncopy) = cache%ener_x(1:ncopy) + prefix_shift
            END IF
            RETURN
         END ASSOCIATE
      END DO

   END SUBROUTINE qs_scf_restore_kpoint_ot_frontier_energies

! **************************************************************************************************
!> \brief clear cached scalar density-response calibration
!> \param scf_env ...
! **************************************************************************************************
   SUBROUTINE qs_scf_clear_kpoint_ot_density_calibration(scf_env)
      TYPE(qs_scf_env_type), POINTER                     :: scf_env

      IF (ASSOCIATED(scf_env%qs_ot_schur_density_calibration%channel)) THEN
         DEALLOCATE (scf_env%qs_ot_schur_density_calibration%channel)
         NULLIFY (scf_env%qs_ot_schur_density_calibration%channel)
      END IF
      scf_env%qs_ot_schur_density_calibration%valid = .FALSE.
      scf_env%qs_ot_schur_density_calibration%nchannel = 0

   END SUBROUTINE qs_scf_clear_kpoint_ot_density_calibration

! **************************************************************************************************
!> \brief inits those objects needed if you want to restart the scf with, say
!>        only a new initial guess, or different density functional or ...
!>        this will happen just before the scf loop starts
!> \param scf_env ...
!> \param qs_env ...
!> \param scf_section ...
!> \par History
!>      03.2006 created [Joost VandeVondele]
! **************************************************************************************************
   SUBROUTINE init_scf_loop(scf_env, qs_env, scf_section)

      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: scf_section

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'init_scf_loop'

      INTEGER                                            :: handle, ispin, nkpoint, nmo, nspin_ot
      INTEGER, DIMENSION(2)                              :: kp_range
      LOGICAL                                            :: do_kpoints, do_rotation, &
                                                            has_unit_metric, init_diag_step
      REAL(KIND=dp), DIMENSION(:), POINTER               :: wkp
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_s
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks_kp, matrix_s_kp
      TYPE(dbcsr_type), POINTER                          :: orthogonality_metric
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(scf_control_type), POINTER                    :: scf_control

      CALL timeset(routineN, handle)

      NULLIFY (scf_control, matrix_s, matrix_ks_kp, matrix_s_kp, matrix_ks, dft_control, mos, mo_coeff, &
               kpoints, wkp)

      CPASSERT(ASSOCIATED(scf_env))
      CPASSERT(ASSOCIATED(qs_env))

      CALL get_qs_env(qs_env=qs_env, &
                      scf_control=scf_control, &
                      dft_control=dft_control, &
                      do_kpoints=do_kpoints, &
                      kpoints=kpoints, &
                      mos=mos)

      nkpoint = 1
      kp_range = 0
      nspin_ot = dft_control%nspins
      IF (dft_control%restricted) nspin_ot = 1
      IF (do_kpoints) CALL get_kpoint_info(kpoints, nkp=nkpoint, kp_range=kp_range, wkp=wkp)

      ! if using mo_coeff_b then copy to fm
      DO ispin = 1, SIZE(mos) !fm->dbcsr
         IF (mos(1)%use_mo_coeff_b) THEN !fm->dbcsr
            CALL copy_dbcsr_to_fm(mos(ispin)%mo_coeff_b, mos(ispin)%mo_coeff) !fm->dbcsr
         END IF !fm->dbcsr
      END DO !fm->dbcsr

      ! this just guarantees that all mo_occupations match the eigenvalues, if smear
      DO ispin = 1, dft_control%nspins
         ! do not reset mo_occupations if the maximum overlap method is in use
         IF (.NOT. scf_control%diagonalization%mom) THEN
            !if the hair probes section is present, this sends hairy_probes to set_mo_occupation subroutine
            !and switches off the standard smearing
            IF (dft_control%hairy_probes .EQV. .TRUE.) THEN
               IF (scf_env%outer_scf%iter_count > 0) THEN
                  scf_control%smear%do_smear = .FALSE.
                  CALL set_mo_occupation(mo_set=mos(ispin), &
                                         smear=scf_control%smear, &
                                         probe=dft_control%probe)
               END IF
            ELSE
               IF (.NOT. scf_control%gce%do_gce) THEN
                  CALL set_mo_occupation(mo_set=mos(ispin), &
                                         smear=scf_control%smear, &
                                         emit_warnings=.NOT. do_kpoints)
               ELSE
                  CALL set_mo_occupation(mo_set=mos(ispin), &
                                         smear=scf_control%smear, &
                                         gce=scf_control%gce, &
                                         emit_warnings=.NOT. do_kpoints)
               END IF
            END IF
         END IF
      END DO

      SELECT CASE (scf_env%method)
      CASE DEFAULT

         CPABORT("Unknown SCF method <"//TRIM(cp_to_string(scf_env%method))//"> found. Check the code!")

      CASE (filter_matrix_diag_method_nr)

         IF (.NOT. scf_env%skip_diis) THEN
            IF (.NOT. ASSOCIATED(scf_env%scf_diis_buffer)) THEN
               ALLOCATE (scf_env%scf_diis_buffer)
               CALL qs_diis_b_create(scf_env%scf_diis_buffer, nbuffer=scf_control%max_diis)
            END IF
            CALL qs_diis_b_clear(scf_env%scf_diis_buffer)
         END IF

      CASE (general_diag_method_nr, special_diag_method_nr, block_krylov_diag_method_nr, smeagol_method_nr)
         IF (.NOT. scf_env%skip_diis) THEN
            IF (do_kpoints) THEN
               IF (.NOT. ASSOCIATED(kpoints%scf_diis_buffer)) THEN
                  ALLOCATE (kpoints%scf_diis_buffer)
                  CALL qs_diis_b_create_kp(kpoints%scf_diis_buffer, nbuffer=scf_control%max_diis)
               END IF
               CALL qs_diis_b_clear_kp(kpoints%scf_diis_buffer)
            ELSE
               IF (.NOT. ASSOCIATED(scf_env%scf_diis_buffer)) THEN
                  ALLOCATE (scf_env%scf_diis_buffer)
                  CALL qs_diis_b_create(scf_env%scf_diis_buffer, nbuffer=scf_control%max_diis)
               END IF
               CALL qs_diis_b_clear(scf_env%scf_diis_buffer)
            END IF
         END IF

      CASE (ot_diag_method_nr)
         CALL get_qs_env(qs_env, matrix_ks=matrix_ks, matrix_s=matrix_s)

         IF (.NOT. scf_env%skip_diis) THEN
            IF (.NOT. ASSOCIATED(scf_env%scf_diis_buffer)) THEN
               ALLOCATE (scf_env%scf_diis_buffer)
               CALL qs_diis_b_create(scf_env%scf_diis_buffer, nbuffer=scf_control%max_diis)
            END IF
            CALL qs_diis_b_clear(scf_env%scf_diis_buffer)
         END IF

         ! disable DFTB and SE for now
         IF (dft_control%qs_control%dftb .OR. &
             dft_control%qs_control%xtb .OR. &
             dft_control%qs_control%semi_empirical) THEN
            CPABORT("DFTB and SE not available with OT/DIAG")
         END IF

         ! if an old preconditioner is still around (i.e. outer SCF is active),
         ! remove it if this could be worthwhile
         CALL restart_preconditioner(qs_env, scf_env%ot_preconditioner, &
                                     scf_control%diagonalization%ot_settings%preconditioner_type, &
                                     dft_control%nspins)

         CALL prepare_preconditioner(qs_env, mos, matrix_ks, matrix_s, scf_env%ot_preconditioner, &
                                     scf_control%diagonalization%ot_settings%preconditioner_type, &
                                     scf_control%diagonalization%ot_settings%precond_solver_type, &
                                     scf_control%diagonalization%ot_settings%energy_gap, dft_control%nspins)

      CASE (block_davidson_diag_method_nr)
         ! Preconditioner initialized within the loop, when required
      CASE (ot_method_nr)
         IF (do_kpoints) THEN
            IF (.NOT. ASSOCIATED(scf_env%qs_ot_env)) THEN
               CALL get_qs_env(qs_env, matrix_ks_kp=matrix_ks_kp, matrix_s_kp=matrix_s_kp)
               CPASSERT(ASSOCIATED(matrix_ks_kp))
               CPASSERT(ASSOCIATED(matrix_s_kp))
               init_diag_step = .FALSE.
               CALL do_general_diag_kp(matrix_ks_kp, matrix_s_kp, kpoints, scf_env, scf_control, &
                                       .FALSE., init_diag_step, store_ot_matrices=.TRUE.)
               CALL qs_scf_restore_added_mos_auto_kp(qs_env, scf_env, scf_control, &
                                                     cp_logger_get_default_io_unit())
               CALL allocate_qs_ot_envs(scf_env=scf_env, &
                                        scf_control=scf_control, &
                                        dft_control=dft_control, &
                                        scf_section=scf_section, &
                                        do_kpoints=do_kpoints, &
                                        nkpoint=nkpoint, &
                                        nspin_ot=nspin_ot, &
                                        kp_range=kp_range, &
                                        wkp=wkp, &
                                        kpoints=kpoints)
               CALL allocate_qs_ot_kpoint_state(qs_env=qs_env, &
                                                scf_env=scf_env, &
                                                scf_control=scf_control, &
                                                dft_control=dft_control, &
                                                kpoints=kpoints, &
                                                matrix_s_kp=matrix_s_kp, &
                                                nspin_ot=nspin_ot)
               CALL prepare_qs_ot_kpoint_preconditioners(qs_env=qs_env, &
                                                         scf_env=scf_env, &
                                                         kpoints=kpoints, &
                                                         matrix_ks_kp=matrix_ks_kp, &
                                                         matrix_s_kp=matrix_s_kp, &
                                                         nspin_ot=nspin_ot)
            END IF

            CALL timestop(handle)
            RETURN
         END IF

         CALL get_qs_env(qs_env, &
                         has_unit_metric=has_unit_metric, &
                         matrix_s=matrix_s, &
                         matrix_ks=matrix_ks)

         ! reortho the wavefunctions if we are having an outer scf and
         ! this is not the first iteration
         ! this is useful to avoid the build-up of numerical noise
         ! however, we can not play this trick if restricted (don't mix non-equivalent orbs)
         IF (scf_control%do_outer_scf_reortho) THEN
            IF (scf_control%outer_scf%have_scf .AND. .NOT. dft_control%restricted) THEN
               IF (scf_env%outer_scf%iter_count > 0) THEN
                  DO ispin = 1, dft_control%nspins
                     CALL get_mo_set(mo_set=mos(ispin), mo_coeff=mo_coeff, nmo=nmo)
                     IF (has_unit_metric) THEN
                        CALL make_basis_simple(mo_coeff, nmo)
                     ELSE
                        CALL make_basis_sm(mo_coeff, nmo, matrix_s(1)%matrix)
                     END IF
                  END DO
               END IF
            END IF
         ELSE
            ! dont need any dirty trick for the numerically stable irac algorithm.
         END IF

         IF (.NOT. ASSOCIATED(scf_env%qs_ot_env)) THEN

            CALL allocate_qs_ot_envs(scf_env=scf_env, &
                                     scf_control=scf_control, &
                                     dft_control=dft_control, &
                                     scf_section=scf_section, &
                                     do_kpoints=do_kpoints, &
                                     nkpoint=nkpoint, &
                                     nspin_ot=nspin_ot, &
                                     kp_range=kp_range, &
                                     wkp=wkp)

            ! might need the KS matrix to init properly
            CALL qs_ks_update_qs_env(qs_env, just_energy=.FALSE., &
                                     calculate_forces=.FALSE.)

            ! if an old preconditioner is still around (i.e. outer SCF is active),
            ! remove it if this could be worthwhile
            IF (.NOT. reuse_precond) THEN
               CALL restart_preconditioner(qs_env, scf_env%ot_preconditioner, &
                                           scf_env%qs_ot_env(1)%settings%preconditioner_type, &
                                           dft_control%nspins)
            END IF

            !
            ! preconditioning still needs to be done correctly with has_unit_metric
            ! notice that a big part of the preconditioning (S^-1) is fine anyhow
            !
            IF (has_unit_metric) THEN
               NULLIFY (orthogonality_metric)
            ELSE
               orthogonality_metric => matrix_s(1)%matrix
            END IF

            IF (.NOT. reuse_precond) THEN
               CALL prepare_preconditioner(qs_env, mos, matrix_ks, matrix_s, scf_env%ot_preconditioner, &
                                           scf_env%qs_ot_env(1)%settings%preconditioner_type, &
                                           scf_env%qs_ot_env(1)%settings%precond_solver_type, &
                                           scf_env%qs_ot_env(1)%settings%energy_gap, dft_control%nspins, &
                                           has_unit_metric=has_unit_metric, &
                                           chol_type=scf_env%qs_ot_env(1)%settings%cholesky_type)
            END IF
            IF (reuse_precond) reuse_precond = .FALSE.

            CALL ot_scf_init(mo_array=mos, matrix_s=orthogonality_metric, &
                             broyden_adaptive_sigma=qs_env%broyden_adaptive_sigma, &
                             qs_ot_env=scf_env%qs_ot_env, matrix_ks=matrix_ks(1)%matrix, &
                             smear=scf_control%smear)

            SELECT CASE (scf_env%qs_ot_env(1)%settings%preconditioner_type)
            CASE (ot_precond_none)
            CASE (ot_precond_full_all, ot_precond_full_single_inverse)
               DO ispin = 1, SIZE(scf_env%qs_ot_env)
                  CALL qs_ot_new_preconditioner(scf_env%qs_ot_env(ispin), &
                                                scf_env%ot_preconditioner(ispin)%preconditioner)
               END DO
            CASE (ot_precond_s_inverse, ot_precond_full_single)
               DO ispin = 1, SIZE(scf_env%qs_ot_env)
                  CALL qs_ot_new_preconditioner(scf_env%qs_ot_env(ispin), &
                                                scf_env%ot_preconditioner(1)%preconditioner)
               END DO
            CASE DEFAULT
               DO ispin = 1, SIZE(scf_env%qs_ot_env)
                  CALL qs_ot_new_preconditioner(scf_env%qs_ot_env(ispin), &
                                                scf_env%ot_preconditioner(1)%preconditioner)
               END DO
            END SELECT
         END IF

         ! if we have non-uniform occupations we should be using rotation
         do_rotation = scf_env%qs_ot_env(1)%settings%do_rotation
         DO ispin = 1, SIZE(mos)
            IF (.NOT. mos(ispin)%uniform_occupation) THEN
               CPASSERT(do_rotation)
            END IF
         END DO
      END SELECT

      ! another safety check
      IF (dft_control%low_spin_roks) THEN
         CPASSERT(scf_env%method == ot_method_nr)
         do_rotation = scf_env%qs_ot_env(1)%settings%do_rotation
         CPASSERT(do_rotation)
      END IF

      CALL timestop(handle)

   END SUBROUTINE init_scf_loop

! **************************************************************************************************
!> \brief allocate OT environments and label their spin/k-point channel identity
!> \param scf_env ...
!> \param scf_control ...
!> \param dft_control ...
!> \param scf_section ...
!> \param do_kpoints ...
!> \param nkpoint ...
!> \param nspin_ot ...
!> \param kp_range ...
!> \param wkp ...
!> \param kpoints ...
! **************************************************************************************************
   SUBROUTINE allocate_qs_ot_envs(scf_env, scf_control, dft_control, scf_section, &
                                  do_kpoints, nkpoint, nspin_ot, kp_range, wkp, kpoints)
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(section_vals_type), POINTER                   :: scf_section
      LOGICAL, INTENT(IN)                                :: do_kpoints
      INTEGER, INTENT(IN)                                :: nkpoint, nspin_ot
      INTEGER, DIMENSION(2), INTENT(IN)                  :: kp_range
      REAL(KIND=dp), DIMENSION(:), POINTER               :: wkp
      TYPE(kpoint_type), OPTIONAL, POINTER               :: kpoints

      INTEGER                                            :: ikpoint, ispin, kpoint_multiplicity, &
                                                            local_kpoint, number_of_OT_envs, &
                                                            ot_channel
      LOGICAL                                            :: do_rotation, is_full_all
      REAL(KIND=dp)                                      :: kpoint_weight

      CPASSERT(.NOT. ASSOCIATED(scf_env%qs_ot_env))

      ! Restricted calculations require just one set of OT orbitals per k-point.
      number_of_OT_envs = qs_ot_number_of_channels(dft_control%nspins, &
                                                   nkpoint=nkpoint, &
                                                   restricted=dft_control%restricted)

      ALLOCATE (scf_env%qs_ot_env(number_of_OT_envs))

      ! XXX Joost XXX should disentangle reading input from this part
      IF (scf_env%outer_scf%iter_count > 0) THEN
         IF (scf_env%iter_delta < scf_control%eps_diis) THEN
            scf_env%qs_ot_env(1)%settings%ot_state = 1
         END IF
      END IF

      CALL ot_scf_read_input(scf_env%qs_ot_env, scf_section)

      IF (scf_env%outer_scf%iter_count > 0) THEN
         IF (scf_env%qs_ot_env(1)%settings%ot_state == 1) THEN
            scf_control%max_scf = MAX(scf_env%qs_ot_env(1)%settings%max_scf_diis, &
                                      scf_control%max_scf)
         END IF
      END IF

      ! Keep a note that we are restricted.
      IF (dft_control%restricted) THEN
         scf_env%qs_ot_env(:)%restricted = .TRUE.
         ! requires rotation
         IF (.NOT. scf_env%qs_ot_env(1)%settings%do_rotation) THEN
            CALL cp_abort(__LOCATION__, &
                          "Restricted calculation with OT requires orbital rotation. Please "// &
                          "activate the OT%ROTATION keyword!")
         END IF
      ELSE
         scf_env%qs_ot_env(:)%restricted = .FALSE.
      END IF

      ! This will rotate the MOs to be eigen states, which is not compatible with rotation.
      ! e.g. mo_derivs here do not yet include potentially different occupation numbers.
      do_rotation = scf_env%qs_ot_env(1)%settings%do_rotation
      ! Only full all needs rotation.
      is_full_all = scf_env%qs_ot_env(1)%settings%preconditioner_type == ot_precond_full_all
      IF (do_rotation .AND. is_full_all .AND. .NOT. do_kpoints) THEN
         CPABORT('PRECONDITIONER FULL_ALL is not compatible with ROTATION.')
      END IF

      DO ikpoint = 1, nkpoint
         local_kpoint = 0
         IF (do_kpoints) THEN
            IF (ikpoint >= kp_range(1) .AND. ikpoint <= kp_range(2)) THEN
               local_kpoint = ikpoint - kp_range(1) + 1
            END IF
         END IF
         DO ispin = 1, nspin_ot
            ot_channel = qs_ot_channel_index(ispin, ikpoint, nspin_ot)
            IF (do_kpoints) THEN
               kpoint_weight = 1.0_dp
               IF (ASSOCIATED(wkp)) kpoint_weight = wkp(ikpoint)
               kpoint_multiplicity = 1
               IF (PRESENT(kpoints)) THEN
                  IF (ASSOCIATED(kpoints)) THEN
                     IF (ASSOCIATED(kpoints%kp_sym)) THEN
                        IF (ikpoint <= SIZE(kpoints%kp_sym)) THEN
                           IF (ASSOCIATED(kpoints%kp_sym(ikpoint)%kpoint_sym)) THEN
                              kpoint_multiplicity = MAX(1, kpoints%kp_sym(ikpoint)%kpoint_sym%nwght)
                           END IF
                        END IF
                     END IF
                  END IF
               END IF
               CALL qs_ot_set_context(scf_env%qs_ot_env(ot_channel), &
                                      spin_index=ispin, &
                                      kpoint_index=ikpoint, &
                                      local_kpoint_index=local_kpoint, &
                                      kpoint_weight=kpoint_weight, &
                                      kpoint_multiplicity=kpoint_multiplicity)
            ELSE
               CALL qs_ot_set_context(scf_env%qs_ot_env(ot_channel), spin_index=ispin)
            END IF
         END DO
      END DO

      CALL qs_ot_check_channel_context(scf_env%qs_ot_env, dft_control%nspins, &
                                       nkpoint=nkpoint, &
                                       restricted=dft_control%restricted, &
                                       require_kpoint=do_kpoints, &
                                       kp_range=kp_range, &
                                       wkp=wkp)

   END SUBROUTINE allocate_qs_ot_envs

! **************************************************************************************************
!> \brief allocate local k-point OT minimizer state for labelled spin/k-point channels
!> \param qs_env ...
!> \param scf_env ...
!> \param scf_control ...
!> \param dft_control ...
!> \param kpoints ...
!> \param matrix_s_kp ...
!> \param nspin_ot ...
! **************************************************************************************************
   SUBROUTINE allocate_qs_ot_kpoint_state(qs_env, scf_env, scf_control, dft_control, &
                                          kpoints, matrix_s_kp, nspin_ot)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s_kp
      INTEGER, INTENT(IN)                                :: nspin_ot

      CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_qs_ot_kpoint_state'

      INTEGER                                            :: handle, ikpoint, ikpoint_local, ispin, &
                                                            nao, nelectron, nkpoint, nmo, &
                                                            ot_channel, restart_compatible_int
      LOGICAL                                            :: restart_compatible, use_real_wfn
      REAL(KIND=dp)                                      :: maxocc
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvalues
      TYPE(cp_fm_type)                                   :: sc_coeff, sc_coeff_im
      TYPE(cp_fm_type), POINTER                          :: mo_coeff, mo_coeff_im
      TYPE(dbcsr_type), POINTER                          :: matrix_sk_im, matrix_sk_re
      TYPE(kpoint_env_type), POINTER                     :: kp
      TYPE(mo_set_type), POINTER                         :: mo_set, mo_set_im

      CALL timeset(routineN, handle)

      NULLIFY (eigenvalues, kp, matrix_sk_im, matrix_sk_re, mo_coeff, mo_coeff_im, &
               mo_set, mo_set_im)

      CPASSERT(ASSOCIATED(qs_env))
      CPASSERT(ASSOCIATED(scf_env%qs_ot_env))
      CPASSERT(ASSOCIATED(kpoints))
      CPASSERT(ASSOCIATED(kpoints%kp_env))
      CPASSERT(ASSOCIATED(matrix_s_kp))
      CPASSERT(ASSOCIATED(matrix_s_kp(1, 1)%matrix))
      CPASSERT(.NOT. dft_control%restricted .OR. nspin_ot == 1)
      IF (scf_env%qs_ot_env(1)%settings%ot_algorithm /= "REF") THEN
         CPABORT("K-point OT currently requires OT%ALGORITHM IRAC.")
      END IF

      CALL get_kpoint_info(kpoints, nkp=nkpoint, use_real_wfn=use_real_wfn)

      ! A distributed minimizer restart is meaningful only if every channel can
      ! restore the same accepted-step state.  Reject the payload globally
      ! before touching any local channel if the current MO layout differs.
      restart_compatible = scf_env%qs_ot_restart%valid .AND. &
                           ASSOCIATED(scf_env%qs_ot_restart%channel)
      IF (restart_compatible) THEN
         restart_compatible = SIZE(scf_env%qs_ot_restart%channel) == SIZE(scf_env%qs_ot_env)
      END IF
      IF (restart_compatible) THEN
         DO ikpoint_local = 1, SIZE(kpoints%kp_env)
            kp => kpoints%kp_env(ikpoint_local)%kpoint_env
            CPASSERT(ASSOCIATED(kp))
            ikpoint = kp%nkpoint
            DO ispin = 1, nspin_ot
               ot_channel = qs_ot_channel_index(ispin, ikpoint, nspin_ot)
               mo_set => kp%mos(1, ispin)
               CALL get_mo_set(mo_set=mo_set, nmo=nmo)
               restart_compatible = restart_compatible .AND. &
                                    scf_env%qs_ot_restart%channel(ot_channel)%valid .AND. &
                                    scf_env%qs_ot_restart%channel(ot_channel)%spin_index == ispin .AND. &
                                    scf_env%qs_ot_restart%channel(ot_channel)%kpoint_index == ikpoint .AND. &
                                    scf_env%qs_ot_restart%channel(ot_channel)%nmo == nmo
            END DO
         END DO
      END IF
      restart_compatible_int = MERGE(1, 0, restart_compatible)
      CPASSERT(ASSOCIATED(kpoints%para_env))
      CALL kpoints%para_env%min(restart_compatible_int)
      IF (scf_env%qs_ot_restart%valid .AND. restart_compatible_int == 0) THEN
         CALL cp_warn(__LOCATION__, &
                      "Ignoring k-point OT restart payload incompatible with the current channel layout")
      END IF
      scf_env%qs_ot_restart%valid = restart_compatible_int /= 0

      DO ikpoint_local = 1, SIZE(kpoints%kp_env)
         kp => kpoints%kp_env(ikpoint_local)%kpoint_env
         CPASSERT(ASSOCIATED(kp))
         CPASSERT(ASSOCIATED(kp%mos))
         CPASSERT(ASSOCIATED(kp%ot_smat))
         CPASSERT(nspin_ot <= SIZE(kp%mos, 2))
         ikpoint = kp%nkpoint
         CPASSERT(ikpoint >= 1)
         CALL dbcsr_init_p(matrix_sk_re)
         CALL copy_fm_to_dbcsr_bc(kp%ot_smat(1), matrix_sk_re)
         IF (.NOT. use_real_wfn) THEN
            CPASSERT(SIZE(kp%ot_smat) >= 2)
            CALL dbcsr_init_p(matrix_sk_im)
            CALL copy_fm_to_dbcsr_bc(kp%ot_smat(2), matrix_sk_im)
         END IF
         DO ispin = 1, nspin_ot
            ot_channel = qs_ot_channel_index(ispin, ikpoint, nspin_ot)
            CPASSERT(.NOT. scf_env%qs_ot_env(ot_channel)%state_allocated)

            NULLIFY (eigenvalues, mo_coeff, mo_coeff_im, mo_set_im)
            mo_set => kp%mos(1, ispin)
            CALL get_mo_set(mo_set=mo_set, mo_coeff=mo_coeff, eigenvalues=eigenvalues, &
                            maxocc=maxocc, nao=nao, nelectron=nelectron, nmo=nmo)
            CPASSERT(ASSOCIATED(mo_coeff))
            CPASSERT(nao > 0)
            CPASSERT(nmo > 0)
            CALL qs_ot_allocate(scf_env%qs_ot_env(ot_channel), &
                                matrix_sk_re, &
                                mo_coeff%matrix_struct)

            CALL cp_fm_create(sc_coeff, mo_coeff%matrix_struct)
            CALL cp_dbcsr_sm_fm_multiply(matrix_sk_re, mo_coeff, sc_coeff, nmo)

            CALL dbcsr_set(scf_env%qs_ot_env(ot_channel)%matrix_c0, 0.0_dp)
            CALL dbcsr_set(scf_env%qs_ot_env(ot_channel)%matrix_sc0, 0.0_dp)
            CALL copy_fm_to_dbcsr(mo_coeff, scf_env%qs_ot_env(ot_channel)%matrix_c0)
            CALL copy_fm_to_dbcsr(sc_coeff, scf_env%qs_ot_env(ot_channel)%matrix_sc0)

            IF (.NOT. use_real_wfn) THEN
               CPASSERT(SIZE(kp%mos, 1) >= 2)
               mo_set_im => kp%mos(2, ispin)
               CALL get_mo_set(mo_set=mo_set_im, mo_coeff=mo_coeff_im)
               CPASSERT(ASSOCIATED(mo_coeff_im))
               CALL allocate_qs_ot_kpoint_complex_state(scf_env%qs_ot_env(ot_channel), &
                                                        matrix_sk_re)
               CALL cp_fm_create(sc_coeff_im, mo_coeff%matrix_struct)
               CALL cp_dbcsr_sm_fm_multiply(matrix_sk_im, mo_coeff_im, sc_coeff, nmo, &
                                            alpha=-1.0_dp, beta=1.0_dp)
               CALL cp_dbcsr_sm_fm_multiply(matrix_sk_re, mo_coeff_im, sc_coeff_im, nmo)
               CALL cp_dbcsr_sm_fm_multiply(matrix_sk_im, mo_coeff, sc_coeff_im, nmo, &
                                            alpha=1.0_dp, beta=1.0_dp)
               CALL dbcsr_set(scf_env%qs_ot_env(ot_channel)%matrix_c0_im, 0.0_dp)
               CALL dbcsr_set(scf_env%qs_ot_env(ot_channel)%matrix_sc0_im, 0.0_dp)
               CALL copy_fm_to_dbcsr(mo_coeff_im, scf_env%qs_ot_env(ot_channel)%matrix_c0_im)
               CALL copy_fm_to_dbcsr(sc_coeff_im, scf_env%qs_ot_env(ot_channel)%matrix_sc0_im)
            END IF

            CALL qs_ot_init(scf_env%qs_ot_env(ot_channel))
            CALL qs_scf_restore_kpoint_ot_density_calibration(scf_env%qs_ot_env(ot_channel), scf_env)
            scf_env%qs_ot_env(ot_channel)%broyden_adaptive_sigma = qs_env%broyden_adaptive_sigma

            CALL dbcsr_set(scf_env%qs_ot_env(ot_channel)%matrix_x, 0.0_dp)
            CALL dbcsr_set(scf_env%qs_ot_env(ot_channel)%matrix_sx, 0.0_dp)
            IF (scf_env%qs_ot_env(ot_channel)%has_complex_kpoint_state) THEN
               CALL dbcsr_set(scf_env%qs_ot_env(ot_channel)%matrix_x_im, 0.0_dp)
               CALL dbcsr_set(scf_env%qs_ot_env(ot_channel)%matrix_sx_im, 0.0_dp)
            END IF

            IF (scf_env%qs_ot_env(ot_channel)%settings%do_rotation) THEN
               CALL dbcsr_set(scf_env%qs_ot_env(ot_channel)%rot_mat_x, 0.0_dp)
               IF (scf_env%qs_ot_env(ot_channel)%has_complex_kpoint_state) THEN
                  CALL dbcsr_set(scf_env%qs_ot_env(ot_channel)%rot_mat_u_im, 0.0_dp)
                  CALL dbcsr_set(scf_env%qs_ot_env(ot_channel)%rot_mat_x_im, 0.0_dp)
               END IF
            END IF

            IF (scf_env%qs_ot_env(ot_channel)%settings%do_ener) THEN
               CALL init_ot_smearing_energies(kp%mos(1, ispin), scf_control%smear)
               IF (ASSOCIATED(eigenvalues)) THEN
                  CPASSERT(SIZE(scf_env%qs_ot_env(ot_channel)%ener_x) == SIZE(eigenvalues))
                  scf_env%qs_ot_env(ot_channel)%ener_x = eigenvalues
                  IF (nkpoint == 1 .AND. scf_env%iter_count <= 1 .AND. &
                      .NOT. scf_env%qs_ot_restart%valid) THEN
                     CALL init_ot_smearing_frontier_energies( &
                        scf_env%qs_ot_env(ot_channel)%ener_x, nelectron, maxocc, scf_control%smear)
                  END IF
               ELSE
                  scf_env%qs_ot_env(ot_channel)%ener_x = 0.0_dp
               END IF
            END IF
            CALL qs_scf_restore_kpoint_ot_frontier_energies( &
               scf_env%qs_ot_env(ot_channel), scf_env)
            SELECT CASE (scf_env%qs_ot_env(ot_channel)%settings%ot_algorithm)
            CASE ("TOD")
               CALL qs_ot_get_p(scf_env%qs_ot_env(ot_channel)%matrix_x, &
                                scf_env%qs_ot_env(ot_channel)%matrix_sx, &
                                scf_env%qs_ot_env(ot_channel))
            CASE ("REF")
               CALL dbcsr_copy(scf_env%qs_ot_env(ot_channel)%matrix_x, &
                               scf_env%qs_ot_env(ot_channel)%matrix_c0)
               CALL dbcsr_copy(scf_env%qs_ot_env(ot_channel)%matrix_sx, &
                               scf_env%qs_ot_env(ot_channel)%matrix_sc0)
               IF (scf_env%qs_ot_env(ot_channel)%has_complex_kpoint_state) THEN
                  CALL dbcsr_copy(scf_env%qs_ot_env(ot_channel)%matrix_x_im, &
                                  scf_env%qs_ot_env(ot_channel)%matrix_c0_im)
                  CALL dbcsr_copy(scf_env%qs_ot_env(ot_channel)%matrix_sx_im, &
                                  scf_env%qs_ot_env(ot_channel)%matrix_sc0_im)
                  CALL qs_ot_get_orbitals_ref_complex(scf_env%qs_ot_env(ot_channel)%matrix_c0, &
                                                      scf_env%qs_ot_env(ot_channel)%matrix_c0_im, &
                                                      matrix_sk_re, matrix_sk_im, &
                                                      scf_env%qs_ot_env(ot_channel))
               END IF
            CASE DEFAULT
               CPABORT("Algorithm not yet implemented")
            END SELECT
            CALL apply_qs_ot_kpoint_restart(scf_env%qs_ot_env(ot_channel), &
                                            scf_env, ot_channel, nmo)
            CALL refresh_qs_ot_kpoint_metric_image(scf_env%qs_ot_env(ot_channel), &
                                                   matrix_sk_re, matrix_sk_im)

            CPASSERT(scf_env%qs_ot_env(ot_channel)%state_allocated)
            IF (.NOT. use_real_wfn) THEN
               CPASSERT(scf_env%qs_ot_env(ot_channel)%has_complex_kpoint_state)
            END IF
            CALL cp_fm_release(sc_coeff)
            IF (.NOT. use_real_wfn) CALL cp_fm_release(sc_coeff_im)
         END DO
         CALL dbcsr_release_p(matrix_sk_re)
         IF (ASSOCIATED(matrix_sk_im)) CALL dbcsr_release_p(matrix_sk_im)
      END DO

      CALL qs_ot_check_channel_context(scf_env%qs_ot_env, dft_control%nspins, &
                                       nkpoint=kpoints%nkp, &
                                       restricted=ANY(scf_env%qs_ot_env(:)%restricted), &
                                       require_kpoint=.TRUE., &
                                       kp_range=kpoints%kp_range, &
                                       require_local_state=.TRUE., &
                                       require_complex_state=.NOT. use_real_wfn)

      CALL qs_scf_clear_kpoint_ot_density_calibration(scf_env)

      CALL timestop(handle)

   END SUBROUTINE allocate_qs_ot_kpoint_state

! **************************************************************************************************
!> \brief rebuild S(k)X after restoring the independent k-point OT coordinates
!> \param qs_ot_env allocated OT channel
!> \param matrix_sk_re real part of the k-point overlap matrix
!> \param matrix_sk_im imaginary part of the k-point overlap matrix
! **************************************************************************************************
   SUBROUTINE refresh_qs_ot_kpoint_metric_image(qs_ot_env, matrix_sk_re, matrix_sk_im)
      TYPE(qs_ot_type)                                   :: qs_ot_env
      TYPE(dbcsr_type), POINTER                          :: matrix_sk_re, matrix_sk_im

      TYPE(dbcsr_type), POINTER                          :: tmp

      CPASSERT(ASSOCIATED(matrix_sk_re))
      CPASSERT(ASSOCIATED(qs_ot_env%matrix_x))
      CPASSERT(ASSOCIATED(qs_ot_env%matrix_sx))

      CALL dbcsr_multiply('N', 'N', 1.0_dp, matrix_sk_re, qs_ot_env%matrix_x, &
                          0.0_dp, qs_ot_env%matrix_sx)
      IF (qs_ot_env%has_complex_kpoint_state) THEN
         CPASSERT(ASSOCIATED(matrix_sk_im))
         CPASSERT(ASSOCIATED(qs_ot_env%matrix_x_im))
         CPASSERT(ASSOCIATED(qs_ot_env%matrix_sx_im))
         CPASSERT(ASSOCIATED(qs_ot_env%buf1_n_k))
         tmp => qs_ot_env%buf1_n_k

         CALL dbcsr_multiply('N', 'N', 1.0_dp, matrix_sk_im, qs_ot_env%matrix_x_im, &
                             0.0_dp, tmp)
         CALL dbcsr_add(qs_ot_env%matrix_sx, tmp, alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)

         CALL dbcsr_multiply('N', 'N', 1.0_dp, matrix_sk_re, qs_ot_env%matrix_x_im, &
                             0.0_dp, qs_ot_env%matrix_sx_im)
         CALL dbcsr_multiply('N', 'N', 1.0_dp, matrix_sk_im, qs_ot_env%matrix_x, &
                             0.0_dp, tmp)
         CALL dbcsr_add(qs_ot_env%matrix_sx_im, tmp, alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
      END IF

   END SUBROUTINE refresh_qs_ot_kpoint_metric_image

! **************************************************************************************************
!> \brief prepare and attach orbital preconditioners for complex k-point OT channels
!> \param qs_env ...
!> \param scf_env ...
!> \param kpoints ...
!> \param matrix_ks_kp ...
!> \param matrix_s_kp ...
!> \param nspin_ot ...
! **************************************************************************************************
   SUBROUTINE prepare_qs_ot_kpoint_preconditioners(qs_env, scf_env, kpoints, &
                                                   matrix_ks_kp, matrix_s_kp, nspin_ot)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks_kp, matrix_s_kp
      INTEGER, INTENT(IN)                                :: nspin_ot

      CHARACTER(LEN=*), PARAMETER :: routineN = 'prepare_qs_ot_kpoint_preconditioners'

      INTEGER                                            :: handle, ikpoint, ispin, local_kpoint, &
                                                            n_ot_channels, ot_channel, prec_type, &
                                                            ref_channel
      LOGICAL                                            :: use_real_wfn
      TYPE(dbcsr_type), POINTER                          :: matrix_k_im, matrix_k_re, matrix_s_im, &
                                                            matrix_s_re
      TYPE(kpoint_env_type), POINTER                     :: kp

      CALL timeset(routineN, handle)

      NULLIFY (kp, matrix_k_im, matrix_k_re, matrix_s_im, matrix_s_re)

      CPASSERT(ASSOCIATED(qs_env))
      CPASSERT(ASSOCIATED(scf_env%qs_ot_env))
      CPASSERT(ASSOCIATED(kpoints))
      CPASSERT(ASSOCIATED(matrix_ks_kp))
      CPASSERT(ASSOCIATED(matrix_s_kp))
      CPASSERT(ASSOCIATED(matrix_ks_kp(1, 1)%matrix))
      CPASSERT(ASSOCIATED(matrix_s_kp(1, 1)%matrix))

      prec_type = scf_env%qs_ot_env(1)%settings%preconditioner_type
      SELECT CASE (prec_type)
      CASE (ot_precond_none)
         CALL timestop(handle)
         RETURN
      CASE (ot_precond_full_single_inverse)
         ! Complex k-point channels are promoted below to the exact Hermitian
         ! FULL_ALL construction.  Applying a real-part inverse separately to
         ! the real and imaginary gradients is not rank invariant.
      CASE (ot_precond_full_all)
         ! Built directly from complex H(k), S(k), and C(k), without rotating
         ! the OT reference gauge.
      CASE DEFAULT
         ! Other preconditioners keep the historical unpreconditioned k-point OT path.
         CALL timestop(handle)
         RETURN
      END SELECT

      n_ot_channels = SIZE(scf_env%qs_ot_env)
      IF (ASSOCIATED(scf_env%ot_preconditioner)) THEN
         DO ot_channel = 1, SIZE(scf_env%ot_preconditioner)
            CALL destroy_preconditioner(scf_env%ot_preconditioner(ot_channel)%preconditioner)
            DEALLOCATE (scf_env%ot_preconditioner(ot_channel)%preconditioner)
         END DO
         DEALLOCATE (scf_env%ot_preconditioner)
         NULLIFY (scf_env%ot_preconditioner)
      END IF

      ref_channel = 0
      DO ot_channel = 1, n_ot_channels
         IF (scf_env%qs_ot_env(ot_channel)%state_allocated) THEN
            ref_channel = ot_channel
            EXIT
         END IF
      END DO
      CPASSERT(ref_channel > 0)
      ALLOCATE (scf_env%ot_preconditioner(n_ot_channels))
      DO ot_channel = 1, n_ot_channels
         ALLOCATE (scf_env%ot_preconditioner(ot_channel)%preconditioner)
         IF (scf_env%qs_ot_env(ot_channel)%state_allocated) THEN
            CALL init_preconditioner(scf_env%ot_preconditioner(ot_channel)%preconditioner, &
                                     para_env=scf_env%qs_ot_env(ot_channel)%para_env, &
                                     blacs_env=scf_env%qs_ot_env(ot_channel)%blacs_env)
         ELSE
            CALL init_preconditioner(scf_env%ot_preconditioner(ot_channel)%preconditioner, &
                                     para_env=scf_env%qs_ot_env(ref_channel)%para_env, &
                                     blacs_env=scf_env%qs_ot_env(ref_channel)%blacs_env)
         END IF
      END DO
      CPASSERT(ASSOCIATED(scf_env%ot_preconditioner))
      CPASSERT(SIZE(scf_env%ot_preconditioner) >= n_ot_channels)

      CALL get_kpoint_info(kpoints, use_real_wfn=use_real_wfn)

      DO local_kpoint = 1, SIZE(kpoints%kp_env)
         kp => kpoints%kp_env(local_kpoint)%kpoint_env
         CPASSERT(ASSOCIATED(kp))
         CPASSERT(ASSOCIATED(kp%mos))
         CPASSERT(ASSOCIATED(kp%ot_hmat))
         CPASSERT(ASSOCIATED(kp%ot_smat))
         ikpoint = kp%nkpoint
         CALL dbcsr_init_p(matrix_s_re)
         CALL copy_fm_to_dbcsr_bc(kp%ot_smat(1), matrix_s_re)
         IF (.NOT. use_real_wfn) THEN
            CPASSERT(SIZE(kp%ot_smat) >= 2)
            CALL dbcsr_init_p(matrix_s_im)
            CALL copy_fm_to_dbcsr_bc(kp%ot_smat(2), matrix_s_im)
         END IF

         DO ispin = 1, nspin_ot
            ot_channel = qs_ot_channel_index(ispin, ikpoint, nspin_ot)
            CPASSERT(scf_env%qs_ot_env(ot_channel)%state_allocated)
            CALL dbcsr_init_p(matrix_k_re)
            CALL copy_fm_to_dbcsr_bc(kp%ot_hmat(1, ispin), matrix_k_re)
            IF (.NOT. use_real_wfn) THEN
               CPASSERT(SIZE(kp%ot_hmat, 1) >= 2)
               CALL dbcsr_init_p(matrix_k_im)
               CALL copy_fm_to_dbcsr_bc(kp%ot_hmat(2, ispin), matrix_k_im)
            END IF
            IF ((prec_type == ot_precond_full_all .OR. &
                 prec_type == ot_precond_full_single_inverse) .AND. .NOT. use_real_wfn) THEN
               CALL make_preconditioner_complex_full_all( &
                  scf_env%ot_preconditioner(ot_channel)%preconditioner, &
                  matrix_k_re, matrix_k_im, matrix_s_re, matrix_s_im, &
                  kp%mos(1, ispin), kp%mos(2, ispin), &
                  scf_env%qs_ot_env(1)%settings%energy_gap)
            ELSE
               CALL make_preconditioner(scf_env%ot_preconditioner(ot_channel)%preconditioner, &
                                        prec_type, &
                                        scf_env%qs_ot_env(1)%settings%precond_solver_type, &
                                        matrix_h=matrix_k_re, &
                                        matrix_s=matrix_s_re, &
                                        mo_set=kp%mos(1, ispin), &
                                        energy_gap=scf_env%qs_ot_env(1)%settings%energy_gap, &
                                        convert_precond_to_dbcsr=.TRUE., &
                                        chol_type=scf_env%qs_ot_env(1)%settings%cholesky_type)
            END IF
            CALL qs_ot_new_preconditioner(scf_env%qs_ot_env(ot_channel), &
                                          scf_env%ot_preconditioner(ot_channel)%preconditioner)
            CALL dbcsr_release_p(matrix_k_re)
            IF (ASSOCIATED(matrix_k_im)) CALL dbcsr_release_p(matrix_k_im)
         END DO
         CALL dbcsr_release_p(matrix_s_re)
         IF (ASSOCIATED(matrix_s_im)) CALL dbcsr_release_p(matrix_s_im)
      END DO

      CALL timestop(handle)

   END SUBROUTINE prepare_qs_ot_kpoint_preconditioners

! **************************************************************************************************
!> \brief apply cached k-point OT restart data to an allocated channel
!> \param qs_ot_env allocated OT channel
!> \param scf_env SCF environment carrying the restart cache
!> \param ot_channel channel index
!> \param nmo number of MOs in the allocated channel
! **************************************************************************************************
   SUBROUTINE apply_qs_ot_kpoint_restart(qs_ot_env, scf_env, ot_channel, nmo)
      TYPE(qs_ot_type)                                   :: qs_ot_env
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      INTEGER, INTENT(IN)                                :: ot_channel, nmo

      INTEGER                                            :: ihist, ndiis, shared_channel
      LOGICAL                                            :: restore_schur_calibration, restored

      qs_ot_env%mermin_ls_component_slope_valid = .FALSE.
      IF (.NOT. scf_env%qs_ot_restart%valid) RETURN
      IF (.NOT. ASSOCIATED(scf_env%qs_ot_restart%channel)) RETURN
      IF (ot_channel < 1 .OR. ot_channel > SIZE(scf_env%qs_ot_restart%channel)) RETURN
      IF (.NOT. scf_env%qs_ot_restart%channel(ot_channel)%valid) RETURN
      IF (scf_env%qs_ot_restart%channel(ot_channel)%spin_index /= qs_ot_env%spin_index) RETURN
      IF (scf_env%qs_ot_restart%channel(ot_channel)%kpoint_index /= qs_ot_env%kpoint_index) RETURN
      IF (scf_env%qs_ot_restart%channel(ot_channel)%nmo /= nmo) RETURN

      ! Matrix and vector history is channel resolved, but OT line-search and
      ! Schur calibration scalars belong to the global minimizer.  Serial
      ! writers historically only advanced those scalars in the first OT
      ! channel, so every distributed k-point group must restore them from the
      ! same canonical restart channel.
      shared_channel = 0
      DO ihist = 1, SIZE(scf_env%qs_ot_restart%channel)
         IF (scf_env%qs_ot_restart%channel(ihist)%valid) THEN
            shared_channel = ihist
            EXIT
         END IF
      END DO
      IF (shared_channel == 0) RETURN

      ndiis = scf_env%qs_ot_restart%channel(ot_channel)%diis_m
      restored = .FALSE.
      ! Version 20 keeps the accepted-step and CG semantics invariant across a
      ! restart, so its serialized response calibration belongs to the restored
      ! Mermin path.  Older payloads are deliberately relearned.
      restore_schur_calibration = scf_env%qs_ot_restart%format_version >= 20
      IF (scf_env%qs_ot_restart%channel(ot_channel)%has_matrix_x .AND. &
          ASSOCIATED(qs_ot_env%matrix_x)) THEN
         CALL copy_qs_ot_restart_dense_matrix( &
            scf_env%qs_ot_restart%channel(ot_channel)%matrix_x_dense, qs_ot_env%matrix_x, qs_ot_env)
         restored = .TRUE.
      END IF
      IF (qs_ot_env%has_complex_kpoint_state .AND. &
          scf_env%qs_ot_restart%channel(ot_channel)%has_matrix_x_im .AND. &
          ASSOCIATED(qs_ot_env%matrix_x_im)) THEN
         CALL copy_qs_ot_restart_dense_matrix( &
            scf_env%qs_ot_restart%channel(ot_channel)%matrix_x_im_dense, qs_ot_env%matrix_x_im, qs_ot_env)
         restored = .TRUE.
      END IF
      IF (qs_ot_env%settings%do_rotation .AND. &
          scf_env%qs_ot_restart%channel(ot_channel)%has_rot_mat_x .AND. &
          ASSOCIATED(qs_ot_env%rot_mat_x)) THEN
         CALL copy_qs_ot_restart_dense_matrix( &
            scf_env%qs_ot_restart%channel(ot_channel)%rot_mat_x_dense, qs_ot_env%rot_mat_x, qs_ot_env)
         restored = .TRUE.
      END IF
      IF (qs_ot_env%settings%do_rotation .AND. qs_ot_env%has_complex_kpoint_state .AND. &
          scf_env%qs_ot_restart%channel(ot_channel)%has_rot_mat_x_im .AND. &
          ASSOCIATED(qs_ot_env%rot_mat_x_im)) THEN
         CALL copy_qs_ot_restart_dense_matrix( &
            scf_env%qs_ot_restart%channel(ot_channel)%rot_mat_x_im_dense, qs_ot_env%rot_mat_x_im, qs_ot_env)
         restored = .TRUE.
      END IF
      IF (qs_ot_env%settings%do_ener .AND. &
          scf_env%qs_ot_restart%channel(ot_channel)%has_ener_x .AND. &
          ASSOCIATED(qs_ot_env%ener_x)) THEN
         IF (SIZE(qs_ot_env%ener_x) == SIZE(scf_env%qs_ot_restart%channel(ot_channel)%ener_x)) THEN
            qs_ot_env%ener_x = scf_env%qs_ot_restart%channel(ot_channel)%ener_x
            restored = .TRUE.
         END IF
      END IF
      IF (scf_env%qs_ot_restart%channel(ot_channel)%has_matrix_dx .AND. &
          ASSOCIATED(qs_ot_env%matrix_dx)) THEN
         CALL copy_qs_ot_restart_dense_matrix( &
            scf_env%qs_ot_restart%channel(ot_channel)%matrix_dx_dense, qs_ot_env%matrix_dx, qs_ot_env)
         restored = .TRUE.
      END IF
      IF (scf_env%qs_ot_restart%channel(ot_channel)%has_matrix_dx_im .AND. &
          ASSOCIATED(qs_ot_env%matrix_dx_im)) THEN
         CALL copy_qs_ot_restart_dense_matrix( &
            scf_env%qs_ot_restart%channel(ot_channel)%matrix_dx_im_dense, qs_ot_env%matrix_dx_im, qs_ot_env)
         restored = .TRUE.
      END IF
      IF (scf_env%qs_ot_restart%channel(ot_channel)%has_matrix_gx_old .AND. &
          ASSOCIATED(qs_ot_env%matrix_gx_old)) THEN
         CALL copy_qs_ot_restart_dense_matrix( &
            scf_env%qs_ot_restart%channel(ot_channel)%matrix_gx_old_dense, qs_ot_env%matrix_gx_old, qs_ot_env)
         restored = .TRUE.
      END IF
      IF (scf_env%qs_ot_restart%channel(ot_channel)%has_matrix_gx_old_im .AND. &
          ASSOCIATED(qs_ot_env%matrix_gx_old_im)) THEN
         CALL copy_qs_ot_restart_dense_matrix( &
            scf_env%qs_ot_restart%channel(ot_channel)%matrix_gx_old_im_dense, qs_ot_env%matrix_gx_old_im, qs_ot_env)
         restored = .TRUE.
      END IF
      IF (qs_ot_env%settings%do_rotation .AND. &
          scf_env%qs_ot_restart%channel(ot_channel)%has_rot_mat_dx .AND. &
          ASSOCIATED(qs_ot_env%rot_mat_dx)) THEN
         CALL copy_qs_ot_restart_dense_matrix( &
            scf_env%qs_ot_restart%channel(ot_channel)%rot_mat_dx_dense, qs_ot_env%rot_mat_dx, qs_ot_env)
         restored = .TRUE.
      END IF
      IF (qs_ot_env%settings%do_rotation .AND. qs_ot_env%has_complex_kpoint_state .AND. &
          scf_env%qs_ot_restart%channel(ot_channel)%has_rot_mat_dx_im .AND. &
          ASSOCIATED(qs_ot_env%rot_mat_dx_im)) THEN
         CALL copy_qs_ot_restart_dense_matrix( &
            scf_env%qs_ot_restart%channel(ot_channel)%rot_mat_dx_im_dense, qs_ot_env%rot_mat_dx_im, qs_ot_env)
         restored = .TRUE.
      END IF
      IF (qs_ot_env%settings%do_rotation .AND. &
          scf_env%qs_ot_restart%channel(ot_channel)%has_rot_mat_gx_old .AND. &
          ASSOCIATED(qs_ot_env%rot_mat_gx_old)) THEN
         CALL copy_qs_ot_restart_dense_matrix( &
            scf_env%qs_ot_restart%channel(ot_channel)%rot_mat_gx_old_dense, qs_ot_env%rot_mat_gx_old, qs_ot_env)
         restored = .TRUE.
      END IF
      IF (qs_ot_env%settings%do_rotation .AND. qs_ot_env%has_complex_kpoint_state .AND. &
          scf_env%qs_ot_restart%channel(ot_channel)%has_rot_mat_gx_old_im .AND. &
          ASSOCIATED(qs_ot_env%rot_mat_gx_old_im)) THEN
         CALL copy_qs_ot_restart_dense_matrix( &
            scf_env%qs_ot_restart%channel(ot_channel)%rot_mat_gx_old_im_dense, &
            qs_ot_env%rot_mat_gx_old_im, qs_ot_env)
         restored = .TRUE.
      END IF
      IF (qs_ot_env%settings%do_ener .AND. &
          scf_env%qs_ot_restart%channel(ot_channel)%has_ener_dx .AND. &
          ASSOCIATED(qs_ot_env%ener_dx)) THEN
         IF (SIZE(qs_ot_env%ener_dx) == SIZE(scf_env%qs_ot_restart%channel(ot_channel)%ener_dx)) THEN
            qs_ot_env%ener_dx = scf_env%qs_ot_restart%channel(ot_channel)%ener_dx
            restored = .TRUE.
         END IF
      END IF
      IF (qs_ot_env%settings%do_ener .AND. &
          scf_env%qs_ot_restart%channel(ot_channel)%has_ener_gx_old .AND. &
          ASSOCIATED(qs_ot_env%ener_gx_old)) THEN
         IF (SIZE(qs_ot_env%ener_gx_old) == SIZE(scf_env%qs_ot_restart%channel(ot_channel)%ener_gx_old)) THEN
            qs_ot_env%ener_gx_old = scf_env%qs_ot_restart%channel(ot_channel)%ener_gx_old
            restored = .TRUE.
         END IF
      END IF
      IF (scf_env%qs_ot_restart%channel(ot_channel)%has_diis_scalars .AND. &
          ASSOCIATED(qs_ot_env%ls_diis) .AND. ndiis == SIZE(qs_ot_env%energy_h)) THEN
         qs_ot_env%ls_diis = scf_env%qs_ot_restart%channel(ot_channel)%ls_diis
         qs_ot_env%lss_diis = scf_env%qs_ot_restart%channel(ot_channel)%lss_diis
         qs_ot_env%c_diis = scf_env%qs_ot_restart%channel(ot_channel)%c_diis
         qs_ot_env%c_broy = scf_env%qs_ot_restart%channel(ot_channel)%c_broy
         qs_ot_env%energy_h = scf_env%qs_ot_restart%channel(ot_channel)%energy_h
         qs_ot_env%ipivot = scf_env%qs_ot_restart%channel(ot_channel)%ipivot
         restored = .TRUE.
      END IF
      IF (scf_env%qs_ot_restart%channel(ot_channel)%has_matrix_history .AND. &
          ASSOCIATED(qs_ot_env%matrix_h_e) .AND. ndiis == SIZE(qs_ot_env%matrix_h_e)) THEN
         DO ihist = 1, ndiis
            CALL copy_qs_ot_restart_dense_matrix( &
               scf_env%qs_ot_restart%channel(ot_channel)%matrix_h_e_dense(ihist), &
               qs_ot_env%matrix_h_e(ihist)%matrix, qs_ot_env)
            CALL copy_qs_ot_restart_dense_matrix( &
               scf_env%qs_ot_restart%channel(ot_channel)%matrix_h_x_dense(ihist), &
               qs_ot_env%matrix_h_x(ihist)%matrix, qs_ot_env)
         END DO
         restored = .TRUE.
      END IF
      IF (scf_env%qs_ot_restart%channel(ot_channel)%has_matrix_history_im .AND. &
          ASSOCIATED(qs_ot_env%matrix_h_e_im) .AND. ndiis == SIZE(qs_ot_env%matrix_h_e_im)) THEN
         DO ihist = 1, ndiis
            CALL copy_qs_ot_restart_dense_matrix( &
               scf_env%qs_ot_restart%channel(ot_channel)%matrix_h_e_im_dense(ihist), &
               qs_ot_env%matrix_h_e_im(ihist)%matrix, qs_ot_env)
            CALL copy_qs_ot_restart_dense_matrix( &
               scf_env%qs_ot_restart%channel(ot_channel)%matrix_h_x_im_dense(ihist), &
               qs_ot_env%matrix_h_x_im(ihist)%matrix, qs_ot_env)
         END DO
         restored = .TRUE.
      END IF
      IF (qs_ot_env%settings%do_rotation .AND. &
          scf_env%qs_ot_restart%channel(ot_channel)%has_rot_history .AND. &
          ASSOCIATED(qs_ot_env%rot_mat_h_e) .AND. ndiis == SIZE(qs_ot_env%rot_mat_h_e)) THEN
         DO ihist = 1, ndiis
            CALL copy_qs_ot_restart_dense_matrix( &
               scf_env%qs_ot_restart%channel(ot_channel)%rot_mat_h_e_dense(ihist), &
               qs_ot_env%rot_mat_h_e(ihist)%matrix, qs_ot_env)
            CALL copy_qs_ot_restart_dense_matrix( &
               scf_env%qs_ot_restart%channel(ot_channel)%rot_mat_h_x_dense(ihist), &
               qs_ot_env%rot_mat_h_x(ihist)%matrix, qs_ot_env)
         END DO
         restored = .TRUE.
      END IF
      IF (qs_ot_env%settings%do_rotation .AND. qs_ot_env%has_complex_kpoint_state .AND. &
          scf_env%qs_ot_restart%channel(ot_channel)%has_rot_history_im .AND. &
          ASSOCIATED(qs_ot_env%rot_mat_h_e_im) .AND. ndiis == SIZE(qs_ot_env%rot_mat_h_e_im)) THEN
         DO ihist = 1, ndiis
            CALL copy_qs_ot_restart_dense_matrix( &
               scf_env%qs_ot_restart%channel(ot_channel)%rot_mat_h_e_im_dense(ihist), &
               qs_ot_env%rot_mat_h_e_im(ihist)%matrix, qs_ot_env)
            CALL copy_qs_ot_restart_dense_matrix( &
               scf_env%qs_ot_restart%channel(ot_channel)%rot_mat_h_x_im_dense(ihist), &
               qs_ot_env%rot_mat_h_x_im(ihist)%matrix, qs_ot_env)
         END DO
         restored = .TRUE.
      END IF
      IF (qs_ot_env%settings%do_ener .AND. &
          scf_env%qs_ot_restart%channel(ot_channel)%has_ener_history .AND. &
          ASSOCIATED(qs_ot_env%ener_h_e) .AND. ndiis == SIZE(qs_ot_env%ener_h_e, 1)) THEN
         IF (SIZE(qs_ot_env%ener_h_e, 2) == SIZE(scf_env%qs_ot_restart%channel(ot_channel)%ener_h_e, 2)) THEN
            qs_ot_env%ener_h_e = scf_env%qs_ot_restart%channel(ot_channel)%ener_h_e
            qs_ot_env%ener_h_x = scf_env%qs_ot_restart%channel(ot_channel)%ener_h_x
            restored = .TRUE.
         END IF
      END IF

      IF (restored) THEN
         qs_ot_env%ot_pos = scf_env%qs_ot_restart%channel(shared_channel)%ot_pos
         qs_ot_env%ot_energy = scf_env%qs_ot_restart%channel(shared_channel)%ot_energy
         qs_ot_env%ot_grad = scf_env%qs_ot_restart%channel(shared_channel)%ot_grad
         qs_ot_env%line_search_left = scf_env%qs_ot_restart%channel(shared_channel)%line_search_left
         qs_ot_env%line_search_right = scf_env%qs_ot_restart%channel(shared_channel)%line_search_right
         qs_ot_env%line_search_mid = scf_env%qs_ot_restart%channel(shared_channel)%line_search_mid
         qs_ot_env%line_search_count = scf_env%qs_ot_restart%channel(shared_channel)%line_search_count
         qs_ot_env%line_search_might_be_done = &
            scf_env%qs_ot_restart%channel(shared_channel)%line_search_might_be_done
         qs_ot_env%delta = scf_env%qs_ot_restart%channel(shared_channel)%delta
         IF (scf_env%qs_ot_restart%format_version >= 27) THEN
            qs_ot_env%mermin_ls_reference_delta = &
               scf_env%qs_ot_restart%channel(shared_channel)%mermin_ls_reference_delta
         ELSE
            qs_ot_env%mermin_ls_reference_delta = qs_ot_env%delta
         END IF
         qs_ot_env%gnorm = scf_env%qs_ot_restart%channel(shared_channel)%gnorm
         qs_ot_env%gnorm_old = scf_env%qs_ot_restart%channel(shared_channel)%gnorm_old
         qs_ot_env%etotal = scf_env%qs_ot_restart%channel(shared_channel)%etotal
         qs_ot_env%gradient = scf_env%qs_ot_restart%channel(shared_channel)%gradient
         qs_ot_env%energy_only = scf_env%qs_ot_restart%channel(shared_channel)%energy_only
         IF (scf_env%qs_ot_restart%format_version >= 19) THEN
            qs_ot_env%ener_direction_linear = &
               scf_env%qs_ot_restart%channel(shared_channel)%ener_direction_linear
         ELSE IF (scf_env%qs_ot_restart%format_version >= 17) THEN
            ! v17/v18 serialized the already trust-projected ener_dx and all
            ! line-search positions, but omitted this reversibility flag.
            qs_ot_env%ener_direction_linear = &
               scf_env%qs_ot_restart%channel(shared_channel)%has_ener_dx .AND. &
               scf_env%qs_ot_restart%channel(shared_channel)%line_search_count > 0
         END IF
         qs_ot_env%diis_iter = scf_env%qs_ot_restart%channel(shared_channel)%diis_iter
         qs_ot_env%ds_min = scf_env%qs_ot_restart%channel(shared_channel)%ds_min
         qs_ot_env%broyden_adaptive_sigma = &
            scf_env%qs_ot_restart%channel(shared_channel)%broyden_adaptive_sigma
         IF (scf_env%qs_ot_restart%format_version >= 18) THEN
            qs_ot_env%mermin_ls_sample_pending = &
               scf_env%qs_ot_restart%channel(shared_channel)%mermin_ls_sample_pending
            qs_ot_env%mermin_ls_best_sample_index = &
               scf_env%qs_ot_restart%channel(shared_channel)%mermin_ls_best_sample_index
            qs_ot_env%mermin_ls_reference_energy = &
               scf_env%qs_ot_restart%channel(shared_channel)%mermin_ls_reference_energy
            qs_ot_env%mermin_ls_best_sample_energy = &
               scf_env%qs_ot_restart%channel(shared_channel)%mermin_ls_best_sample_energy
            qs_ot_env%mermin_ls_best_sample_drop = &
               scf_env%qs_ot_restart%channel(shared_channel)%mermin_ls_best_sample_drop
            qs_ot_env%mermin_ls_model_drop = &
               scf_env%qs_ot_restart%channel(shared_channel)%mermin_ls_model_drop
            qs_ot_env%mermin_ls_model_quality = &
               scf_env%qs_ot_restart%channel(shared_channel)%mermin_ls_model_quality
            qs_ot_env%mermin_ls_best_sample_pos = &
               scf_env%qs_ot_restart%channel(shared_channel)%mermin_ls_best_sample_pos
            qs_ot_env%mermin_ls_next_pos = &
               scf_env%qs_ot_restart%channel(shared_channel)%mermin_ls_next_pos
            qs_ot_env%mermin_ls_accepted_quality_samples = &
               scf_env%qs_ot_restart%channel(shared_channel)%mermin_ls_accepted_quality_samples
            qs_ot_env%mermin_ls_accepted_quality = &
               scf_env%qs_ot_restart%channel(shared_channel)%mermin_ls_accepted_quality
            qs_ot_env%mermin_ls_accepted_quality_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%mermin_ls_accepted_quality_ema
            qs_ot_env%mermin_ls_accepted_gap = &
               scf_env%qs_ot_restart%channel(shared_channel)%mermin_ls_accepted_gap
            qs_ot_env%mermin_ls_accepted_gap_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%mermin_ls_accepted_gap_ema
         END IF
         IF (restore_schur_calibration .AND. scf_env%qs_ot_restart%format_version >= 4) THEN
            qs_ot_env%schur_projector_density_pending = .FALSE.
            qs_ot_env%schur_projector_density_audit_pending = .FALSE.
            qs_ot_env%schur_projector_density_descent_pending = .FALSE.
            qs_ot_env%schur_projector_physical_pending = .FALSE.
            qs_ot_env%schur_projector_physical_direct_selected = 0.0_dp
            qs_ot_env%schur_projector_density_descent_selected = 0.0_dp
            qs_ot_env%schur_projector_response_selected = 0.0_dp
            qs_ot_env%schur_projector_density_full_selected = 0.0_dp
            qs_ot_env%schur_projector_physical_drop = 0.0_dp
            qs_ot_env%schur_projector_physical_metric = 0.0_dp
            qs_ot_env%schur_projector_physical_weight = 0.0_dp
            qs_ot_env%schur_projector_physical_correction_drop = 0.0_dp
            qs_ot_env%schur_projector_physical_compensation = 0.0_dp
            qs_ot_env%schur_projector_physical_negative = 0.0_dp
            qs_ot_env%schur_projector_physical_direct = 0.0_dp
            qs_ot_env%schur_projector_physical_selected = 0.0_dp
            qs_ot_env%schur_projector_physical_shift = 0.0_dp
            qs_ot_env%schur_projector_physical_trace_actual = 0.0_dp
            qs_ot_env%schur_projector_physical_trace_raw = 0.0_dp
            qs_ot_env%schur_projector_physical_trace_active = 0.0_dp
            qs_ot_env%schur_projector_physical_trace_actual_diag = 0.0_dp
            qs_ot_env%schur_projector_physical_trace_actual_free = 0.0_dp
            qs_ot_env%schur_projector_physical_trace_raw_diag = 0.0_dp
            qs_ot_env%schur_projector_physical_trace_raw_free = 0.0_dp
            qs_ot_env%schur_projector_physical_trace_active_diag = 0.0_dp
            qs_ot_env%schur_projector_physical_trace_active_free = 0.0_dp
            qs_ot_env%schur_projector_physical_trace_invariant = 0.0_dp
            qs_ot_env%schur_projector_physical_trace_invariant_diag = 0.0_dp
            qs_ot_env%schur_projector_physical_trace_invariant_free = 0.0_dp
            qs_ot_env%schur_projector_density_quality_samples = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_quality_samples
            qs_ot_env%schur_projector_density_audit_samples = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_audit_samples
            qs_ot_env%schur_projector_density_descent_quality_samples = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_descent_quality_samples
            qs_ot_env%schur_projector_density_descent_selected_quality_samples = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_descent_selected_quality_samples
            qs_ot_env%schur_projector_response_selected_quality_samples = &
               scf_env%qs_ot_restart%channel(shared_channel)%response_selected_quality_samples
            qs_ot_env%schur_projector_density_full_selected_quality_samples = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_full_selected_quality_samples
            qs_ot_env%schur_projector_physical_quality_samples = &
               scf_env%qs_ot_restart%channel(shared_channel)%physical_quality_samples
            qs_ot_env%schur_projector_physical_selected_quality_samples = &
               scf_env%qs_ot_restart%channel(shared_channel)%physical_selected_quality_samples
            qs_ot_env%schur_projector_physical_direct_quality_samples = &
               scf_env%qs_ot_restart%channel(shared_channel)%physical_direct_quality_samples
            qs_ot_env%schur_projector_physical_direct_selected_quality_samples = &
               scf_env%qs_ot_restart%channel(shared_channel)%physical_direct_selected_quality_samples
            qs_ot_env%schur_projector_phys_direct_density_sel_quality_samples = &
               scf_env%qs_ot_restart%channel(shared_channel)%phys_direct_density_sel_quality_samples
            qs_ot_env%schur_projector_mixed_direct_quality_samples = &
               scf_env%qs_ot_restart%channel(shared_channel)%mixed_direct_quality_samples
            qs_ot_env%schur_projector_spectral_selected_quality_samples = &
               scf_env%qs_ot_restart%channel(shared_channel)%spectral_selected_quality_samples
            qs_ot_env%schur_projector_spectral_selected_failure_samples = &
               scf_env%qs_ot_restart%channel(shared_channel)%spectral_selected_failure_samples
            qs_ot_env%schur_projector_density_quality = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_quality
            qs_ot_env%schur_projector_density_quality_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_quality_ema
            qs_ot_env%schur_projector_density_drop_sum = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_drop_sum
            qs_ot_env%schur_projector_density_measured_sum = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_measured_sum
            qs_ot_env%schur_projector_density_lowcurv_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_lowcurv_ema
            qs_ot_env%schur_projector_density_lowcurv_quality_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_lowcurv_quality_ema
            qs_ot_env%schur_projector_density_audit_amplification_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_audit_amplification_ema
            qs_ot_env%schur_projector_density_audit_amplification_quality_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_audit_amplification_quality_ema
            qs_ot_env%schur_projector_density_audit_lowcurv_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_audit_lowcurv_ema
            qs_ot_env%schur_projector_density_audit_lowcurv_quality_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_audit_lowcurv_quality_ema
            qs_ot_env%schur_projector_density_descent_quality = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_descent_quality
            qs_ot_env%schur_projector_density_descent_quality_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_descent_quality_ema
            qs_ot_env%schur_projector_density_descent_drop_sum = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_descent_drop_sum
            qs_ot_env%schur_projector_density_descent_measured_sum = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_descent_measured_sum
            qs_ot_env%schur_projector_density_descent_coeff_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_descent_coeff_ema
            qs_ot_env%schur_projector_density_descent_support = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_descent_support
            qs_ot_env%schur_projector_density_descent_support_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_descent_support_ema
            qs_ot_env%schur_projector_density_descent_support_ready = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_descent_support_ready
            qs_ot_env%schur_projector_physical_quality = &
               scf_env%qs_ot_restart%channel(shared_channel)%physical_quality
            qs_ot_env%schur_projector_physical_quality_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%physical_quality_ema
            qs_ot_env%schur_projector_physical_drop_sum = &
               scf_env%qs_ot_restart%channel(shared_channel)%physical_drop_sum
            qs_ot_env%schur_projector_physical_measured_sum = &
               scf_env%qs_ot_restart%channel(shared_channel)%physical_measured_sum
            qs_ot_env%schur_projector_physical_selected_quality = &
               scf_env%qs_ot_restart%channel(shared_channel)%physical_selected_quality
            qs_ot_env%schur_projector_physical_selected_quality_line_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%physical_selected_quality_ema
            qs_ot_env%schur_projector_physical_selected_residual_quality = &
               scf_env%qs_ot_restart%channel(shared_channel)%physical_selected_residual_quality
            qs_ot_env%schur_projector_physical_selected_residual_quality_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%physical_selected_residual_quality_ema
            qs_ot_env%schur_projector_physical_selected_drop_sum = &
               scf_env%qs_ot_restart%channel(shared_channel)%physical_selected_drop_sum
            qs_ot_env%schur_projector_physical_selected_measured_sum = &
               scf_env%qs_ot_restart%channel(shared_channel)%physical_selected_measured_sum
            qs_ot_env%schur_projector_physical_trace_actual_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%physical_trace_actual_ema
            qs_ot_env%schur_projector_physical_trace_raw_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%physical_trace_raw_ema
            qs_ot_env%schur_projector_physical_trace_active_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%physical_trace_active_ema
            qs_ot_env%schur_projector_physical_trace_actual_diag_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%physical_trace_actual_diag_ema
            qs_ot_env%schur_projector_physical_trace_actual_free_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%physical_trace_actual_free_ema
            qs_ot_env%schur_projector_physical_trace_raw_diag_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%physical_trace_raw_diag_ema
            qs_ot_env%schur_projector_physical_trace_raw_free_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%physical_trace_raw_free_ema
            qs_ot_env%schur_projector_physical_trace_active_diag_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%physical_trace_active_diag_ema
            qs_ot_env%schur_projector_physical_trace_active_free_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%physical_trace_active_free_ema
            qs_ot_env%schur_projector_physical_trace_invariant_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%physical_trace_invariant_ema
            qs_ot_env%schur_projector_physical_trace_invariant_diag_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%physical_trace_invariant_diag_ema
            qs_ot_env%schur_projector_physical_trace_invariant_free_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%physical_trace_invariant_free_ema
            qs_ot_env%schur_projector_physical_trace_actual_quality_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%physical_trace_actual_quality_ema
            qs_ot_env%schur_projector_physical_trace_raw_quality_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%physical_trace_raw_quality_ema
            qs_ot_env%schur_projector_physical_trace_active_quality_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%physical_trace_active_quality_ema
            qs_ot_env%schur_projector_physical_trace_actual_diag_quality_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%physical_trace_actual_diag_quality_ema
            qs_ot_env%schur_projector_physical_trace_actual_free_quality_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%physical_trace_actual_free_quality_ema
            qs_ot_env%schur_projector_physical_trace_raw_diag_quality_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%physical_trace_raw_diag_quality_ema
            qs_ot_env%schur_projector_physical_trace_raw_free_quality_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%physical_trace_raw_free_quality_ema
            qs_ot_env%schur_projector_physical_trace_active_diag_quality_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%physical_trace_active_diag_quality_ema
            qs_ot_env%schur_projector_physical_trace_active_free_quality_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%physical_trace_active_free_quality_ema
            qs_ot_env%schur_projector_physical_trace_invariant_quality_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%physical_trace_invariant_quality_ema
            qs_ot_env%schur_projector_physical_trace_invariant_diag_quality_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%physical_trace_invariant_diag_quality_ema
            qs_ot_env%schur_projector_physical_trace_invariant_free_quality_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%physical_trace_invariant_free_quality_ema
            qs_ot_env%schur_projector_density_descent_selected_quality = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_descent_selected_quality
            qs_ot_env%schur_projector_density_descent_selected_quality_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_descent_selected_quality_ema
            qs_ot_env%schur_projector_density_descent_selected_residual_quality = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_descent_selected_residual_quality
            qs_ot_env%schur_projector_density_descent_selected_residual_quality_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_descent_selected_residual_quality_ema
            qs_ot_env%schur_projector_density_descent_selected_drop_sum = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_descent_selected_drop_sum
            qs_ot_env%schur_projector_density_descent_selected_measured_sum = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_descent_selected_measured_sum
            qs_ot_env%schur_projector_response_selected_quality = &
               scf_env%qs_ot_restart%channel(shared_channel)%response_selected_quality
            qs_ot_env%schur_projector_response_selected_quality_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%response_selected_quality_ema
            qs_ot_env%schur_projector_response_selected_residual_quality = &
               scf_env%qs_ot_restart%channel(shared_channel)%response_selected_residual_quality
            qs_ot_env%schur_projector_response_selected_residual_quality_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%response_selected_residual_quality_ema
            qs_ot_env%schur_projector_response_selected_drop_sum = &
               scf_env%qs_ot_restart%channel(shared_channel)%response_selected_drop_sum
            qs_ot_env%schur_projector_response_selected_measured_sum = &
               scf_env%qs_ot_restart%channel(shared_channel)%response_selected_measured_sum
            qs_ot_env%schur_projector_density_full_selected_quality = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_full_selected_quality
            qs_ot_env%schur_projector_density_full_selected_quality_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_full_selected_quality_ema
            qs_ot_env%schur_projector_density_full_selected_residual_quality = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_full_selected_residual_quality
            qs_ot_env%schur_projector_density_full_selected_residual_quality_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_full_selected_residual_quality_ema
            qs_ot_env%schur_projector_density_full_selected_drop_sum = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_full_selected_drop_sum
            qs_ot_env%schur_projector_density_full_selected_measured_sum = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_full_selected_measured_sum
            qs_ot_env%schur_projector_physical_direct_quality = &
               scf_env%qs_ot_restart%channel(shared_channel)%physical_direct_quality
            qs_ot_env%schur_projector_physical_direct_quality_shadow_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%physical_direct_quality_shadow_ema
            qs_ot_env%schur_projector_physical_direct_drop_sum = &
               scf_env%qs_ot_restart%channel(shared_channel)%physical_direct_drop_sum
            qs_ot_env%schur_projector_physical_direct_measured_sum = &
               scf_env%qs_ot_restart%channel(shared_channel)%physical_direct_measured_sum
            qs_ot_env%schur_projector_physical_direct_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%physical_direct_ema
            qs_ot_env%schur_projector_physical_direct_selected_quality = &
               scf_env%qs_ot_restart%channel(shared_channel)%physical_direct_selected_quality
            qs_ot_env%schur_projector_physical_direct_selected_quality_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%physical_direct_selected_quality_ema
            qs_ot_env%schur_projector_physical_direct_selected_residual_quality = &
               scf_env%qs_ot_restart%channel(shared_channel)%physical_direct_selected_residual_quality
            qs_ot_env%schur_projector_physical_direct_selected_residual_quality_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%physical_direct_selected_residual_quality_ema
            qs_ot_env%schur_projector_physical_direct_selected_drop_sum = &
               scf_env%qs_ot_restart%channel(shared_channel)%physical_direct_selected_drop_sum
            qs_ot_env%schur_projector_physical_direct_selected_measured_sum = &
               scf_env%qs_ot_restart%channel(shared_channel)%physical_direct_selected_measured_sum
            qs_ot_env%schur_projector_phys_direct_density_sel_quality = &
               scf_env%qs_ot_restart%channel(shared_channel)%phys_direct_density_sel_quality
            qs_ot_env%schur_projector_phys_direct_density_sel_quality_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%phys_direct_density_sel_quality_ema
            qs_ot_env%schur_projector_phys_direct_density_sel_resid_quality = &
               scf_env%qs_ot_restart%channel(shared_channel)%phys_direct_density_sel_resid_quality
            qs_ot_env%schur_projector_phys_direct_density_sel_resid_quality_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%phys_direct_density_sel_resid_quality_ema
            qs_ot_env%schur_projector_phys_direct_density_sel_drop_sum = &
               scf_env%qs_ot_restart%channel(shared_channel)%phys_direct_density_sel_drop_sum
            qs_ot_env%schur_projector_phys_direct_density_sel_measured_sum = &
               scf_env%qs_ot_restart%channel(shared_channel)%phys_direct_density_sel_measured_sum
            qs_ot_env%schur_projector_mixed_direct_quality = &
               scf_env%qs_ot_restart%channel(shared_channel)%mixed_direct_quality
            qs_ot_env%schur_projector_mixed_direct_quality_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%mixed_direct_quality_ema
            qs_ot_env%schur_projector_mixed_direct_drop_sum = &
               scf_env%qs_ot_restart%channel(shared_channel)%mixed_direct_drop_sum
            qs_ot_env%schur_projector_mixed_direct_measured_sum = &
               scf_env%qs_ot_restart%channel(shared_channel)%mixed_direct_measured_sum
            qs_ot_env%schur_projector_mixed_direct_weight = &
               scf_env%qs_ot_restart%channel(shared_channel)%mixed_direct_weight
            qs_ot_env%schur_projector_mixed_direct_weight_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%mixed_direct_weight_ema
            qs_ot_env%schur_projector_spectral_selected_quality = &
               scf_env%qs_ot_restart%channel(shared_channel)%spectral_selected_quality
            qs_ot_env%schur_projector_spectral_selected_quality_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%spectral_selected_quality_ema
            qs_ot_env%schur_projector_spectral_selected_residual_quality = &
               scf_env%qs_ot_restart%channel(shared_channel)%spectral_selected_residual_quality
            qs_ot_env%schur_projector_spectral_selected_residual_quality_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%spectral_selected_residual_quality_ema
            qs_ot_env%schur_projector_spectral_selected_drop_sum = &
               scf_env%qs_ot_restart%channel(shared_channel)%spectral_selected_drop_sum
            qs_ot_env%schur_projector_spectral_selected_measured_sum = &
               scf_env%qs_ot_restart%channel(shared_channel)%spectral_selected_measured_sum
         END IF
         IF (scf_env%qs_ot_restart%format_version >= 18) THEN
            qs_ot_env%schur_projector_density_full_pending = .FALSE.
            qs_ot_env%schur_projector_density_full_direct_pending = .FALSE.
            qs_ot_env%schur_projector_density_full_global_pending = .FALSE.
            qs_ot_env%schur_projector_density_full_drop = 0.0_dp
            qs_ot_env%schur_projector_density_full_selected = 0.0_dp
            qs_ot_env%schur_projector_density_full_direct_drop = 0.0_dp
            qs_ot_env%schur_projector_density_fdir_selected = 0.0_dp
            qs_ot_env%schur_projector_density_full_global_drop = 0.0_dp
            qs_ot_env%schur_projector_density_full_selected_quality_samples = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_full_selected_quality_samples
            qs_ot_env%schur_projector_density_full_direct_quality_samples = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_full_direct_quality_samples
            qs_ot_env%schur_projector_density_fdir_sel_quality_samples = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_fdir_sel_quality_samples
            qs_ot_env%schur_projector_density_full_global_quality_samples = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_full_global_quality_samples
            qs_ot_env%schur_projector_density_full_global_sel_quality_samples = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_full_global_sel_quality_samples
            qs_ot_env%schur_projector_density_full_selected_quality = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_full_selected_quality
            qs_ot_env%schur_projector_density_full_selected_quality_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_full_selected_quality_ema
            qs_ot_env%schur_projector_density_full_selected_residual_quality = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_full_selected_residual_quality
            qs_ot_env%schur_projector_density_full_selected_residual_quality_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_full_selected_residual_quality_ema
            qs_ot_env%schur_projector_density_full_selected_drop_sum = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_full_selected_drop_sum
            qs_ot_env%schur_projector_density_full_selected_measured_sum = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_full_selected_measured_sum
            qs_ot_env%schur_projector_density_full_direct_quality = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_full_direct_quality
            qs_ot_env%schur_projector_density_full_direct_quality_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_full_direct_quality_ema
            qs_ot_env%schur_projector_density_full_direct_drop_sum = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_full_direct_drop_sum
            qs_ot_env%schur_projector_density_full_direct_measured_sum = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_full_direct_measured_sum
            qs_ot_env%schur_projector_density_full_direct_trust = 1.0_dp
            qs_ot_env%schur_projector_density_full_direct_trust_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_full_direct_trust_ema
            qs_ot_env%schur_projector_density_full_direct_weighted_trust = 0.0_dp
            qs_ot_env%schur_projector_density_full_direct_weighted_trust_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_full_direct_weighted_trust_ema
            qs_ot_env%schur_projector_density_full_direct_abs_trust = 0.0_dp
            qs_ot_env%schur_projector_density_full_direct_abs_trust_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_full_direct_abs_trust_ema
            qs_ot_env%schur_projector_density_fdir_sel_quality = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_fdir_sel_quality
            qs_ot_env%schur_projector_density_fdir_sel_quality_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_fdir_sel_quality_ema
            qs_ot_env%schur_projector_density_fdir_sel_resid_quality = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_fdir_sel_resid_quality
            qs_ot_env%schur_projector_density_fdir_sel_resid_quality_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_fdir_sel_resid_quality_ema
            qs_ot_env%schur_projector_density_fdir_sel_drop_sum = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_fdir_sel_drop_sum
            qs_ot_env%schur_projector_density_fdir_sel_measured_sum = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_fdir_sel_measured_sum
            qs_ot_env%schur_projector_density_full_global_quality = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_full_global_quality
            qs_ot_env%schur_projector_density_full_global_quality_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_full_global_quality_ema
            qs_ot_env%schur_projector_density_full_global_drop_sum = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_full_global_drop_sum
            qs_ot_env%schur_projector_density_full_global_measured_sum = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_full_global_measured_sum
            qs_ot_env%schur_projector_density_full_global_trust = 1.0_dp
            qs_ot_env%schur_projector_density_full_global_trust_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_full_global_trust_ema
            qs_ot_env%schur_projector_density_full_global_weighted_trust = 0.0_dp
            qs_ot_env%schur_projector_density_full_global_weighted_trust_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_full_global_weighted_trust_ema
            qs_ot_env%schur_projector_density_full_global_abs_trust = 0.0_dp
            qs_ot_env%schur_projector_density_full_global_abs_trust_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_full_global_abs_trust_ema
            qs_ot_env%schur_projector_density_full_global_selected = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_full_global_selected
            qs_ot_env%schur_projector_density_full_global_sel_quality = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_full_global_sel_quality
            qs_ot_env%schur_projector_density_full_global_sel_quality_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_full_global_sel_quality_ema
            qs_ot_env%schur_projector_density_full_global_sel_resid_quality = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_full_global_sel_residual_quality
            qs_ot_env%schur_projector_density_full_global_sel_resid_quality_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_full_global_sel_residual_quality_ema
            qs_ot_env%schur_projector_density_full_global_sel_drop_sum = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_full_global_sel_drop_sum
            qs_ot_env%schur_projector_density_full_global_sel_measured_sum = &
               scf_env%qs_ot_restart%channel(shared_channel)%density_full_global_sel_measured_sum
         END IF
         IF (scf_env%qs_ot_restart%format_version >= 20) THEN
            qs_ot_env%schur_model_attempt = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_int(1)
            qs_ot_env%schur_cg_attempt = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_int(2)
            qs_ot_env%schur_model_flavor = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_int(3)
            qs_ot_env%schur_cg_flavor = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_int(4)
            qs_ot_env%schur_model_quality_samples = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_int(5)
            qs_ot_env%schur_constrained_quality_samples = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_int(6)
            qs_ot_env%schur_split_quality_samples = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_int(7)
            qs_ot_env%schur_coupled_quality_samples = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_int(8)
            qs_ot_env%schur_split_shadow_quality_samples = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_int(9)
            qs_ot_env%schur_occgap_shadow_quality_samples = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_int(10)
            qs_ot_env%schur_dense_metric_shadow_quality_samples = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_int(11)
            qs_ot_env%schur_step_gauge_quality_samples = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_int(12)
            qs_ot_env%kpoint_physical_delta_groups = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_int(13)

            qs_ot_env%schur_model_pending = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_log(1)
            qs_ot_env%schur_model_quality_ready = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_log(2)
            qs_ot_env%schur_constrained_pending = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_log(3)
            qs_ot_env%schur_constrained_candidate_enabled = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_log(4)
            qs_ot_env%schur_split_shadow_pending = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_log(5)
            qs_ot_env%schur_occgap_shadow_pending = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_log(6)
            qs_ot_env%schur_dense_metric_shadow_pending = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_log(7)
            qs_ot_env%schur_history_mixed = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_log(8)
            qs_ot_env%schur_split_candidate_enabled = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_log(9)
            qs_ot_env%schur_residual_history_reset_pending = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_log(10)
            qs_ot_env%kpoint_physical_delta_ready = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_log(11)

            qs_ot_env%schur_model_drop = scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(1)
            qs_ot_env%schur_model_calibrated_drop = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(2)
            qs_ot_env%schur_model_measured_drop = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(3)
            qs_ot_env%schur_model_quality = scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(4)
            qs_ot_env%schur_model_quality_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(5)
            qs_ot_env%schur_model_calibrated_quality = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(6)
            qs_ot_env%schur_model_calibrated_quality_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(7)
            qs_ot_env%schur_model_gap_ratio = scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(8)
            qs_ot_env%schur_model_delta = scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(9)
            qs_ot_env%schur_model_control_delta = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(10)
            qs_ot_env%schur_model_residual_drop = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(11)
            qs_ot_env%schur_residual_history_quality = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(12)
            qs_ot_env%schur_step_gap_norm = scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(13)
            qs_ot_env%schur_step_gauge_norm = scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(14)
            qs_ot_env%schur_step_fermi_norm = scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(15)
            qs_ot_env%schur_step_gauge_ratio = scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(16)
            qs_ot_env%schur_step_gauge_ratio_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(17)
            qs_ot_env%schur_step_gauge_residual_quality = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(18)
            qs_ot_env%schur_step_gauge_residual_quality_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(19)
            qs_ot_env%kpoint_raw_delta = scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(20)
            qs_ot_env%kpoint_physical_delta = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(21)
            qs_ot_env%kpoint_projected_delta = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(22)
            qs_ot_env%kpoint_eps_scf = scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(23)
            qs_ot_env%kpoint_physical_gnorm = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(24)
            qs_ot_env%kpoint_physical_rot_gnorm = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(25)
            qs_ot_env%kpoint_physical_energy_gnorm = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(26)
            qs_ot_env%kpoint_gauge_gnorm = scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(27)
            qs_ot_env%schur_constrained_drop = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(28)
            qs_ot_env%schur_constrained_calibrated_drop = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(29)
            qs_ot_env%schur_constrained_quality = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(30)
            qs_ot_env%schur_constrained_quality_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(31)
            qs_ot_env%schur_constrained_calibrated_quality = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(32)
            qs_ot_env%schur_constrained_calibrated_quality_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(33)
            qs_ot_env%schur_constrained_drop_sum = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(34)
            qs_ot_env%schur_constrained_measured_sum = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(35)
            qs_ot_env%schur_constrained_calibrated_drop_sum = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(36)
            qs_ot_env%schur_constrained_gap_ratio = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(37)
            qs_ot_env%schur_split_quality = scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(38)
            qs_ot_env%schur_split_quality_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(39)
            qs_ot_env%schur_split_calibrated_quality = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(40)
            qs_ot_env%schur_split_calibrated_quality_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(41)
            qs_ot_env%schur_split_residual_quality = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(42)
            qs_ot_env%schur_split_residual_quality_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(43)
            qs_ot_env%schur_split_drop_sum = scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(44)
            qs_ot_env%schur_split_measured_sum = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(45)
            qs_ot_env%schur_split_calibrated_drop_sum = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(46)
            qs_ot_env%schur_split_shadow_drop = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(47)
            qs_ot_env%schur_split_shadow_calibrated_drop = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(48)
            qs_ot_env%schur_split_shadow_quality = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(49)
            qs_ot_env%schur_split_shadow_quality_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(50)
            qs_ot_env%schur_split_shadow_calibrated_quality = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(51)
            qs_ot_env%schur_split_shadow_calibrated_quality_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(52)
            qs_ot_env%schur_split_shadow_drop_sum = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(53)
            qs_ot_env%schur_split_shadow_measured_sum = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(54)
            qs_ot_env%schur_split_shadow_calibrated_drop_sum = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(55)
            qs_ot_env%schur_occgap_shadow_drop = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(56)
            qs_ot_env%schur_occgap_shadow_quality = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(57)
            qs_ot_env%schur_occgap_shadow_quality_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(58)
            qs_ot_env%schur_occgap_shadow_drop_sum = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(59)
            qs_ot_env%schur_occgap_shadow_measured_sum = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(60)
            qs_ot_env%schur_occgap_shadow_gap_ratio = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(61)
            qs_ot_env%schur_coupled_quality = scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(62)
            qs_ot_env%schur_coupled_quality_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(63)
            qs_ot_env%schur_coupled_calibrated_quality = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(64)
            qs_ot_env%schur_coupled_calibrated_quality_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(65)
            qs_ot_env%schur_coupled_residual_quality = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(66)
            qs_ot_env%schur_coupled_residual_quality_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(67)
            qs_ot_env%schur_coupled_drop_sum = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(68)
            qs_ot_env%schur_coupled_measured_sum = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(69)
            qs_ot_env%schur_coupled_calibrated_drop_sum = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(70)
            qs_ot_env%schur_dense_metric_shadow_drop = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(71)
            qs_ot_env%schur_dense_metric_shadow_quality = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(72)
            qs_ot_env%schur_dense_metric_shadow_quality_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(73)
            qs_ot_env%schur_dense_metric_shadow_drop_sum = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(74)
            qs_ot_env%schur_dense_metric_shadow_measured_sum = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(75)
         END IF
         IF (scf_env%qs_ot_restart%format_version >= 21) THEN
            qs_ot_env%schur_projector_mode_quality_samples = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_int(14)
            qs_ot_env%schur_projector_mode_pending = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_log(12)
            qs_ot_env%schur_projector_mode_drop = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(76)
            qs_ot_env%schur_projector_mode_response_drop = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(77)
            qs_ot_env%schur_projector_mode_step_negative = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(78)
            qs_ot_env%schur_projector_mode_response_negative = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(79)
            qs_ot_env%schur_projector_mode_trust = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(80)
            qs_ot_env%schur_projector_mode_quality = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(81)
            qs_ot_env%schur_projector_mode_quality_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(82)
            qs_ot_env%schur_projector_mode_drop_sum = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(83)
            qs_ot_env%schur_projector_mode_measured_sum = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(84)
            qs_ot_env%schur_projector_mode_response_drop_sum = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(85)
            qs_ot_env%schur_projector_mode_step_negative_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(86)
            qs_ot_env%schur_projector_mode_response_negative_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(87)
            qs_ot_env%schur_projector_mode_trust_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(88)
         END IF
         IF (scf_env%qs_ot_restart%format_version >= 22) THEN
            qs_ot_env%kpoint_occ_preconditioner_ready_samples = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_int(15)
            qs_ot_env%kpoint_occ_preconditioner_active = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_log(13)
            qs_ot_env%kpoint_occ_metric_residual_quality_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(89)
         ELSE IF (qs_ot_env%has_complex_kpoint_state .AND. &
                  qs_ot_env%settings%occupation_preconditioner .AND. &
                  scf_env%qs_ot_restart%channel(1)%mermin_ls_accepted_quality_samples > 0 .AND. &
                  scf_env%qs_ot_restart%channel(1)%delta > 0.0_dp .AND. &
                  scf_env%qs_ot_restart%channel(1)%delta <= 7.5E-2_dp) THEN
            ! Before v22 the complex k-point occupation preconditioner was
            ! implicitly active.  Preserve a locally calibrated old restart,
            ! but relearn the new residual-aware orbital metric.
            qs_ot_env%kpoint_occ_preconditioner_active = .TRUE.
            qs_ot_env%kpoint_occ_preconditioner_ready_samples = 2
            qs_ot_env%kpoint_occ_metric_residual_quality_ema = 0.0_dp
         END IF
         IF (scf_env%qs_ot_restart%format_version >= 24) THEN
            qs_ot_env%schur_model_density_stationarity_delta = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(90)
            qs_ot_env%schur_density_stationarity_quality = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(91)
            qs_ot_env%schur_density_stationarity_quality_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(92)
            qs_ot_env%schur_density_stationarity_quality_samples = &
               MAX(0, NINT(scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(93)))
            qs_ot_env%schur_projector_density_full_selected_stationarity_quality = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(94)
            qs_ot_env%schur_projector_density_full_selected_stationarity_quality_ema = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(95)
            qs_ot_env%schur_projector_density_full_selected_stationarity_samples = &
               MAX(0, NINT(scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(96)))
         END IF
         IF (scf_env%qs_ot_restart%format_version >= 25) THEN
            qs_ot_env%schur_model_slope = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(97)
            qs_ot_env%schur_model_quadratic = &
               scf_env%qs_ot_restart%channel(shared_channel)%adaptive_real(98)
         ELSE
            ! Version 24 and older calibrated accepted directions against the
            ! unit-step quadratic model even when the line search accepted at
            ! another coordinate.  Keep orbitals and minimizer history, but
            ! relearn only the incompatible accepted-direction statistics.
            CALL qs_ot_reset_accepted_model_calibration(qs_ot_env)
         END IF
         IF (scf_env%qs_ot_restart%format_version >= 20 .AND. &
             scf_env%qs_ot_restart%format_version < 31 .AND. &
             qs_ot_env%has_complex_kpoint_state .AND. &
             qs_ot_env%settings%occupation_preconditioner .AND. &
             qs_ot_env%settings%do_rotation .AND. qs_ot_env%settings%do_ener .AND. &
             qs_ot_env%line_search_count > 0) THEN
            ! Version 27 adds the immutable projected residual at the start of
            ! a Mermin line search. Version 28 reconnects the regular
            ! canonical-density Newton tangent. Version 29 replaces, rather
            ! than duplicates, the density part of that response. Version 30
            ! projects the rotation floor only in the separate density
            ! candidate. Version 31 distinguishes a complete projected
            ! reference from the bootstrap search residual. An older open
            ! search therefore belongs to an incompatible candidate model.
            ! Accept its current orbitals as a fresh point and discard only
            ! the open direction.
            qs_ot_env%line_search_left = -1
            qs_ot_env%line_search_right = -1
            qs_ot_env%line_search_mid = -1
            qs_ot_env%line_search_count = 0
            qs_ot_env%line_search_might_be_done = .FALSE.
            qs_ot_env%energy_only = .FALSE.
            qs_ot_env%mermin_ls_sample_pending = .FALSE.
            qs_ot_env%mermin_ls_closure_rejected = .FALSE.
            qs_ot_env%mermin_ls_best_sample_index = 0
            qs_ot_env%mermin_ls_reference_delta = 0.0_dp
            qs_ot_env%ot_pos(:) = -1.0_dp
            qs_ot_env%ot_energy(:) = -1.0_dp
            qs_ot_env%ot_grad(:) = -1.0_dp
            qs_ot_env%gnorm = 0.0_dp
            qs_ot_env%gnorm_old = 1.0_dp
            qs_ot_env%gradient = 0.0_dp
            qs_ot_env%ds_min = qs_ot_env%settings%ds_min
            IF (ASSOCIATED(qs_ot_env%matrix_dx)) CALL dbcsr_set(qs_ot_env%matrix_dx, 0.0_dp)
            IF (ASSOCIATED(qs_ot_env%matrix_dx_im)) CALL dbcsr_set(qs_ot_env%matrix_dx_im, 0.0_dp)
            IF (ASSOCIATED(qs_ot_env%matrix_gx_old)) CALL dbcsr_set(qs_ot_env%matrix_gx_old, 0.0_dp)
            IF (ASSOCIATED(qs_ot_env%matrix_gx_old_im)) CALL dbcsr_set(qs_ot_env%matrix_gx_old_im, 0.0_dp)
            IF (ASSOCIATED(qs_ot_env%rot_mat_dx)) CALL dbcsr_set(qs_ot_env%rot_mat_dx, 0.0_dp)
            IF (ASSOCIATED(qs_ot_env%rot_mat_dx_im)) CALL dbcsr_set(qs_ot_env%rot_mat_dx_im, 0.0_dp)
            IF (ASSOCIATED(qs_ot_env%rot_mat_gx_old)) CALL dbcsr_set(qs_ot_env%rot_mat_gx_old, 0.0_dp)
            IF (ASSOCIATED(qs_ot_env%rot_mat_gx_old_im)) CALL dbcsr_set(qs_ot_env%rot_mat_gx_old_im, 0.0_dp)
            IF (ASSOCIATED(qs_ot_env%ener_dx)) qs_ot_env%ener_dx(:) = 0.0_dp
            IF (ASSOCIATED(qs_ot_env%ener_gx_old)) qs_ot_env%ener_gx_old(:) = 0.0_dp
            qs_ot_env%ener_direction_linear = .FALSE.
         END IF
         IF (scf_env%qs_ot_restart%format_version < 30) THEN
            ! The v30 exact-density candidate uses a modally projected floor
            ! and calibrates residual quality against canonical density
            ! stationarity. Retain physical secants and minimizer coordinates,
            ! but relearn the incompatible accepted-candidate statistics.
            CALL qs_ot_reset_accepted_model_calibration(qs_ot_env)
         END IF
         IF (scf_env%qs_ot_restart%format_version < 20) THEN
            ! Version 20 keeps the structural CG direction flags invariant.
            ! Older k-point runs could switch USE_DX off during a canonical
            ! bootstrap or overlap guard, so their open line-search directions
            ! cannot be continued with a freshly allocated CG environment.
            qs_ot_env%line_search_left = -1
            qs_ot_env%line_search_right = -1
            qs_ot_env%line_search_mid = -1
            qs_ot_env%line_search_count = 0
            qs_ot_env%line_search_might_be_done = .FALSE.
            qs_ot_env%energy_only = .FALSE.
            qs_ot_env%mermin_ls_sample_pending = .FALSE.
            qs_ot_env%mermin_ls_best_sample_index = 0
            qs_ot_env%ot_pos(:) = -1.0_dp
            qs_ot_env%ot_energy(:) = -1.0_dp
            qs_ot_env%ot_grad(:) = -1.0_dp
            qs_ot_env%gnorm = 0.0_dp
            qs_ot_env%gnorm_old = 1.0_dp
            qs_ot_env%gradient = 0.0_dp
            qs_ot_env%diis_iter = 0
            qs_ot_env%ds_min = qs_ot_env%settings%ds_min
            IF (ASSOCIATED(qs_ot_env%matrix_dx)) CALL dbcsr_set(qs_ot_env%matrix_dx, 0.0_dp)
            IF (ASSOCIATED(qs_ot_env%matrix_dx_im)) CALL dbcsr_set(qs_ot_env%matrix_dx_im, 0.0_dp)
            IF (ASSOCIATED(qs_ot_env%matrix_gx_old)) CALL dbcsr_set(qs_ot_env%matrix_gx_old, 0.0_dp)
            IF (ASSOCIATED(qs_ot_env%matrix_gx_old_im)) CALL dbcsr_set(qs_ot_env%matrix_gx_old_im, 0.0_dp)
            IF (ASSOCIATED(qs_ot_env%rot_mat_dx)) CALL dbcsr_set(qs_ot_env%rot_mat_dx, 0.0_dp)
            IF (ASSOCIATED(qs_ot_env%rot_mat_dx_im)) CALL dbcsr_set(qs_ot_env%rot_mat_dx_im, 0.0_dp)
            IF (ASSOCIATED(qs_ot_env%rot_mat_gx_old)) CALL dbcsr_set(qs_ot_env%rot_mat_gx_old, 0.0_dp)
            IF (ASSOCIATED(qs_ot_env%rot_mat_gx_old_im)) CALL dbcsr_set(qs_ot_env%rot_mat_gx_old_im, 0.0_dp)
            IF (ASSOCIATED(qs_ot_env%ener_dx)) qs_ot_env%ener_dx(:) = 0.0_dp
            IF (ASSOCIATED(qs_ot_env%ener_gx_old)) qs_ot_env%ener_gx_old(:) = 0.0_dp
            qs_ot_env%ener_direction_linear = .FALSE.
         END IF
      END IF

   END SUBROUTINE apply_qs_ot_kpoint_restart

! **************************************************************************************************
!> \brief copy replicated restart data into the owning k-point communicator
!> \param source replicated dense restart matrix
!> \param matrix local OT DBCSR matrix
!> \param qs_ot_env owning OT channel
! **************************************************************************************************
   SUBROUTINE copy_qs_ot_restart_dense_matrix(source, matrix, qs_ot_env)
      TYPE(qs_ot_restart_dense_matrix_type), INTENT(IN)  :: source
      TYPE(dbcsr_type), POINTER                          :: matrix
      TYPE(qs_ot_type), INTENT(IN)                       :: qs_ot_env

      TYPE(cp_fm_struct_type), POINTER                   :: local_fm_struct
      TYPE(cp_fm_type)                                   :: local_fm

      CPASSERT(ALLOCATED(source%data))
      CPASSERT(ASSOCIATED(matrix))
      CPASSERT(ASSOCIATED(qs_ot_env%blacs_env))

      NULLIFY (local_fm_struct)
      CALL cp_fm_struct_create(local_fm_struct, context=qs_ot_env%blacs_env, &
                               nrow_global=SIZE(source%data, 1), &
                               ncol_global=SIZE(source%data, 2))
      CALL cp_fm_create(local_fm, local_fm_struct)
      CALL cp_fm_set_submatrix(local_fm, source%data)
      CALL copy_fm_to_dbcsr(local_fm, matrix)
      CALL cp_fm_release(local_fm)
      CALL cp_fm_struct_release(local_fm_struct)

   END SUBROUTINE copy_qs_ot_restart_dense_matrix

! **************************************************************************************************
!> \brief allocate imaginary orbital images for complex k-point OT channel state
!> \param qs_ot_env ...
!> \param matrix_s ...
! **************************************************************************************************
   SUBROUTINE allocate_qs_ot_kpoint_complex_state(qs_ot_env, matrix_s)
      TYPE(qs_ot_type)                                   :: qs_ot_env
      TYPE(dbcsr_type), POINTER                          :: matrix_s

      INTEGER                                            :: i, m_diis, nmo, nmo_ortho

      CPASSERT(qs_ot_env%state_allocated)
      CPASSERT(.NOT. qs_ot_env%has_complex_kpoint_state)
      CPASSERT(ASSOCIATED(matrix_s))

      m_diis = qs_ot_env%settings%diis_m

      CALL dbcsr_get_info(qs_ot_env%matrix_c0, nfullcols_total=nmo)
      CALL dbcsr_get_info(qs_ot_env%matrix_sc0, nfullcols_total=nmo_ortho)

      CALL dbcsr_init_p(qs_ot_env%matrix_c0_im)
      CALL cp_dbcsr_m_by_n_from_row_template(qs_ot_env%matrix_c0_im, template=matrix_s, n=nmo, &
                                             sym=dbcsr_type_no_symmetry)
      CALL dbcsr_init_p(qs_ot_env%matrix_sc0_im)
      CALL cp_dbcsr_m_by_n_from_row_template(qs_ot_env%matrix_sc0_im, template=matrix_s, n=nmo_ortho, &
                                             sym=dbcsr_type_no_symmetry)
      CALL dbcsr_init_p(qs_ot_env%matrix_ref_orth)
      CALL cp_dbcsr_m_by_n_from_row_template(qs_ot_env%matrix_ref_orth, template=matrix_s, n=nmo, &
                                             sym=dbcsr_type_no_symmetry)
      CALL dbcsr_init_p(qs_ot_env%matrix_ref_orth_im)
      CALL cp_dbcsr_m_by_n_from_row_template(qs_ot_env%matrix_ref_orth_im, template=matrix_s, n=nmo, &
                                             sym=dbcsr_type_no_symmetry)
      CALL dbcsr_init_p(qs_ot_env%matrix_x_im)
      CALL cp_dbcsr_m_by_n_from_row_template(qs_ot_env%matrix_x_im, template=matrix_s, n=nmo, &
                                             sym=dbcsr_type_no_symmetry)
      CALL dbcsr_init_p(qs_ot_env%matrix_sx_im)
      CALL cp_dbcsr_m_by_n_from_row_template(qs_ot_env%matrix_sx_im, template=matrix_s, n=nmo, &
                                             sym=dbcsr_type_no_symmetry)
      CALL dbcsr_init_p(qs_ot_env%matrix_gx_im)
      CALL cp_dbcsr_m_by_n_from_row_template(qs_ot_env%matrix_gx_im, template=matrix_s, n=nmo, &
                                             sym=dbcsr_type_no_symmetry)
      IF (qs_ot_env%use_dx) THEN
         CALL dbcsr_init_p(qs_ot_env%matrix_dx_im)
         CALL cp_dbcsr_m_by_n_from_row_template(qs_ot_env%matrix_dx_im, template=matrix_s, n=nmo, &
                                                sym=dbcsr_type_no_symmetry)
      END IF
      IF (qs_ot_env%use_gx_old) THEN
         CALL dbcsr_init_p(qs_ot_env%matrix_gx_old_im)
         CALL cp_dbcsr_m_by_n_from_row_template(qs_ot_env%matrix_gx_old_im, template=matrix_s, n=nmo, &
                                                sym=dbcsr_type_no_symmetry)
      END IF

      IF (qs_ot_env%settings%ot_method == "DIIS" .OR. &
          qs_ot_env%settings%ot_method == "BROY") THEN
         CALL dbcsr_allocate_matrix_set(qs_ot_env%matrix_h_e_im, m_diis)
         CALL dbcsr_allocate_matrix_set(qs_ot_env%matrix_h_x_im, m_diis)
         DO i = 1, m_diis
            CALL dbcsr_init_p(qs_ot_env%matrix_h_x_im(i)%matrix)
            CALL cp_dbcsr_m_by_n_from_row_template(qs_ot_env%matrix_h_x_im(i)%matrix, template=matrix_s, &
                                                   n=nmo, sym=dbcsr_type_no_symmetry)
            CALL dbcsr_init_p(qs_ot_env%matrix_h_e_im(i)%matrix)
            CALL cp_dbcsr_m_by_n_from_row_template(qs_ot_env%matrix_h_e_im(i)%matrix, template=matrix_s, &
                                                   n=nmo, sym=dbcsr_type_no_symmetry)
         END DO
      END IF

      IF (qs_ot_env%settings%do_rotation) THEN
         CALL dbcsr_init_p(qs_ot_env%rot_mat_pgx)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%rot_mat_pgx, template=matrix_s, m=nmo, n=nmo, &
                                            sym=dbcsr_type_no_symmetry)
         CALL dbcsr_set(qs_ot_env%rot_mat_pgx, 0.0_dp)
         CALL dbcsr_init_p(qs_ot_env%rot_mat_u_im)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%rot_mat_u_im, template=matrix_s, m=nmo, n=nmo, &
                                            sym=dbcsr_type_no_symmetry)
         CALL dbcsr_init_p(qs_ot_env%rot_mat_x_im)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%rot_mat_x_im, template=matrix_s, m=nmo, n=nmo, &
                                            sym=dbcsr_type_no_symmetry)
         CALL dbcsr_init_p(qs_ot_env%rot_mat_dedu_im)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%rot_mat_dedu_im, template=matrix_s, m=nmo, n=nmo, &
                                            sym=dbcsr_type_no_symmetry)
         CALL dbcsr_init_p(qs_ot_env%rot_mat_chc_im)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%rot_mat_chc_im, template=matrix_s, m=nmo, n=nmo, &
                                            sym=dbcsr_type_no_symmetry)
         CALL dbcsr_init_p(qs_ot_env%rot_mat_gx_im)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%rot_mat_gx_im, template=matrix_s, m=nmo, n=nmo, &
                                            sym=dbcsr_type_no_symmetry)
         CALL dbcsr_init_p(qs_ot_env%rot_mat_pgx_im)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%rot_mat_pgx_im, template=matrix_s, m=nmo, n=nmo, &
                                            sym=dbcsr_type_no_symmetry)
         CALL dbcsr_set(qs_ot_env%rot_mat_pgx_im, 0.0_dp)
         IF (qs_ot_env%use_dx) THEN
            CALL dbcsr_init_p(qs_ot_env%rot_mat_dx_im)
            CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%rot_mat_dx_im, template=matrix_s, m=nmo, n=nmo, &
                                               sym=dbcsr_type_no_symmetry)
         END IF
         IF (qs_ot_env%use_gx_old) THEN
            CALL dbcsr_init_p(qs_ot_env%rot_mat_gx_old_im)
            CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%rot_mat_gx_old_im, template=matrix_s, m=nmo, n=nmo, &
                                               sym=dbcsr_type_no_symmetry)
         END IF
         CALL dbcsr_init_p(qs_ot_env%rot_mat_mermin_g0)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%rot_mat_mermin_g0, template=matrix_s, m=nmo, n=nmo, &
                                            sym=dbcsr_type_no_symmetry)
         CALL dbcsr_init_p(qs_ot_env%rot_mat_mermin_g0_im)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%rot_mat_mermin_g0_im, template=matrix_s, m=nmo, n=nmo, &
                                            sym=dbcsr_type_no_symmetry)
         CALL dbcsr_init_p(qs_ot_env%rot_mat_mermin_s)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%rot_mat_mermin_s, template=matrix_s, m=nmo, n=nmo, &
                                            sym=dbcsr_type_no_symmetry)
         CALL dbcsr_init_p(qs_ot_env%rot_mat_mermin_s_im)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%rot_mat_mermin_s_im, template=matrix_s, m=nmo, n=nmo, &
                                            sym=dbcsr_type_no_symmetry)
         CALL dbcsr_init_p(qs_ot_env%rot_mat_mermin_y)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%rot_mat_mermin_y, template=matrix_s, m=nmo, n=nmo, &
                                            sym=dbcsr_type_no_symmetry)
         CALL dbcsr_init_p(qs_ot_env%rot_mat_mermin_y_im)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%rot_mat_mermin_y_im, template=matrix_s, m=nmo, n=nmo, &
                                            sym=dbcsr_type_no_symmetry)
         CALL dbcsr_init_p(qs_ot_env%rot_mat_mermin_s_previous)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%rot_mat_mermin_s_previous, template=matrix_s, &
                                            m=nmo, n=nmo, sym=dbcsr_type_no_symmetry)
         CALL dbcsr_init_p(qs_ot_env%rot_mat_mermin_s_previous_im)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%rot_mat_mermin_s_previous_im, template=matrix_s, &
                                            m=nmo, n=nmo, sym=dbcsr_type_no_symmetry)
         CALL dbcsr_init_p(qs_ot_env%rot_mat_mermin_y_previous)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%rot_mat_mermin_y_previous, template=matrix_s, &
                                            m=nmo, n=nmo, sym=dbcsr_type_no_symmetry)
         CALL dbcsr_init_p(qs_ot_env%rot_mat_mermin_y_previous_im)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%rot_mat_mermin_y_previous_im, template=matrix_s, &
                                            m=nmo, n=nmo, sym=dbcsr_type_no_symmetry)
         CALL dbcsr_init_p(qs_ot_env%rot_mat_mermin_h0_ref)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%rot_mat_mermin_h0_ref, template=matrix_s, &
                                            m=nmo, n=nmo, sym=dbcsr_type_no_symmetry)
         CALL dbcsr_init_p(qs_ot_env%rot_mat_mermin_h0_ref_im)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%rot_mat_mermin_h0_ref_im, template=matrix_s, &
                                            m=nmo, n=nmo, sym=dbcsr_type_no_symmetry)
         CALL dbcsr_init_p(qs_ot_env%rot_mat_mermin_frozen_g1)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%rot_mat_mermin_frozen_g1, template=matrix_s, &
                                            m=nmo, n=nmo, sym=dbcsr_type_no_symmetry)
         CALL dbcsr_init_p(qs_ot_env%rot_mat_mermin_frozen_g1_im)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%rot_mat_mermin_frozen_g1_im, template=matrix_s, &
                                            m=nmo, n=nmo, sym=dbcsr_type_no_symmetry)
         CALL dbcsr_init_p(qs_ot_env%rot_mat_mermin_frozen_reduced_y)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%rot_mat_mermin_frozen_reduced_y, template=matrix_s, &
                                            m=nmo, n=nmo, sym=dbcsr_type_no_symmetry)
         CALL dbcsr_set(qs_ot_env%rot_mat_mermin_frozen_reduced_y, 0.0_dp)
         CALL dbcsr_init_p(qs_ot_env%rot_mat_mermin_frozen_reduced_y_im)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%rot_mat_mermin_frozen_reduced_y_im, template=matrix_s, &
                                            m=nmo, n=nmo, sym=dbcsr_type_no_symmetry)
         CALL dbcsr_set(qs_ot_env%rot_mat_mermin_frozen_reduced_y_im, 0.0_dp)
         CALL dbcsr_init_p(qs_ot_env%rot_mat_mermin_closure_y)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%rot_mat_mermin_closure_y, template=matrix_s, &
                                            m=nmo, n=nmo, sym=dbcsr_type_no_symmetry)
         CALL dbcsr_init_p(qs_ot_env%rot_mat_mermin_closure_y_im)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%rot_mat_mermin_closure_y_im, template=matrix_s, &
                                            m=nmo, n=nmo, sym=dbcsr_type_no_symmetry)
         CALL dbcsr_init_p(qs_ot_env%rot_mat_mermin_closure_y_previous)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%rot_mat_mermin_closure_y_previous, template=matrix_s, &
                                            m=nmo, n=nmo, sym=dbcsr_type_no_symmetry)
         CALL dbcsr_init_p(qs_ot_env%rot_mat_mermin_closure_y_previous_im)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%rot_mat_mermin_closure_y_previous_im, template=matrix_s, &
                                            m=nmo, n=nmo, sym=dbcsr_type_no_symmetry)
         CALL dbcsr_init_p(qs_ot_env%rot_mat_mermin_s_anchor)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%rot_mat_mermin_s_anchor, template=matrix_s, &
                                            m=nmo, n=nmo, sym=dbcsr_type_no_symmetry)
         CALL dbcsr_init_p(qs_ot_env%rot_mat_mermin_s_anchor_im)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%rot_mat_mermin_s_anchor_im, template=matrix_s, &
                                            m=nmo, n=nmo, sym=dbcsr_type_no_symmetry)
         CALL dbcsr_init_p(qs_ot_env%rot_mat_mermin_y_anchor)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%rot_mat_mermin_y_anchor, template=matrix_s, &
                                            m=nmo, n=nmo, sym=dbcsr_type_no_symmetry)
         CALL dbcsr_init_p(qs_ot_env%rot_mat_mermin_y_anchor_im)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%rot_mat_mermin_y_anchor_im, template=matrix_s, &
                                            m=nmo, n=nmo, sym=dbcsr_type_no_symmetry)
         CALL dbcsr_init_p(qs_ot_env%rot_mat_mermin_closure_y_anchor)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%rot_mat_mermin_closure_y_anchor, template=matrix_s, &
                                            m=nmo, n=nmo, sym=dbcsr_type_no_symmetry)
         CALL dbcsr_init_p(qs_ot_env%rot_mat_mermin_closure_y_anchor_im)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%rot_mat_mermin_closure_y_anchor_im, template=matrix_s, &
                                            m=nmo, n=nmo, sym=dbcsr_type_no_symmetry)
         IF (qs_ot_env%settings%ot_method == "DIIS") THEN
            CALL dbcsr_allocate_matrix_set(qs_ot_env%rot_mat_h_e_im, m_diis)
            CALL dbcsr_allocate_matrix_set(qs_ot_env%rot_mat_h_x_im, m_diis)
            DO i = 1, m_diis
               CALL dbcsr_init_p(qs_ot_env%rot_mat_h_x_im(i)%matrix)
               CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%rot_mat_h_x_im(i)%matrix, template=matrix_s, &
                                                  m=nmo, n=nmo, sym=dbcsr_type_no_symmetry)
               CALL dbcsr_init_p(qs_ot_env%rot_mat_h_e_im(i)%matrix)
               CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%rot_mat_h_e_im(i)%matrix, template=matrix_s, &
                                                  m=nmo, n=nmo, sym=dbcsr_type_no_symmetry)
            END DO
         END IF
      END IF

      IF (qs_ot_env%settings%do_ener) THEN
         ALLOCATE (qs_ot_env%ener_pgx(nmo))
         qs_ot_env%ener_pgx(:) = 0.0_dp
         ALLOCATE (qs_ot_env%ener_mermin_g0(nmo), qs_ot_env%ener_mermin_x0(nmo), &
                   qs_ot_env%ener_mermin_s(nmo), &
                   qs_ot_env%ener_mermin_y(nmo), qs_ot_env%ener_mermin_s_previous(nmo), &
                   qs_ot_env%ener_mermin_y_previous(nmo), qs_ot_env%ener_mermin_frozen_diag(nmo), &
                   qs_ot_env%ener_mermin_frozen_g1(nmo), qs_ot_env%ener_mermin_closure_y(nmo), &
                   qs_ot_env%ener_mermin_closure_y_previous(nmo), &
                   qs_ot_env%ener_mermin_s_anchor(nmo), &
                   qs_ot_env%ener_mermin_y_anchor(nmo), &
                   qs_ot_env%ener_mermin_closure_y_anchor(nmo))
         qs_ot_env%ener_mermin_g0(:) = 0.0_dp
         qs_ot_env%ener_mermin_x0(:) = 0.0_dp
         qs_ot_env%ener_mermin_s(:) = 0.0_dp
         qs_ot_env%ener_mermin_y(:) = 0.0_dp
         qs_ot_env%ener_mermin_s_previous(:) = 0.0_dp
         qs_ot_env%ener_mermin_y_previous(:) = 0.0_dp
         qs_ot_env%ener_mermin_frozen_diag(:) = 0.0_dp
         qs_ot_env%ener_mermin_frozen_g1(:) = 0.0_dp
         qs_ot_env%ener_mermin_closure_y(:) = 0.0_dp
         qs_ot_env%ener_mermin_closure_y_previous(:) = 0.0_dp
         qs_ot_env%ener_mermin_s_anchor(:) = 0.0_dp
         qs_ot_env%ener_mermin_y_anchor(:) = 0.0_dp
         qs_ot_env%ener_mermin_closure_y_anchor(:) = 0.0_dp
      END IF

      ALLOCATE (qs_ot_env%kpoint_density_pair_residual(nmo, nmo))
      qs_ot_env%kpoint_density_pair_residual(:, :) = 0.0_dp
      ALLOCATE (qs_ot_env%kpoint_canonical_occupation(nmo))
      qs_ot_env%kpoint_canonical_occupation(:) = 0.0_dp
      ALLOCATE (qs_ot_env%kpoint_canonical_eigenvalue(nmo))
      qs_ot_env%kpoint_canonical_eigenvalue(:) = 0.0_dp
      qs_ot_env%kpoint_canonical_occupation_ready = .FALSE.
      qs_ot_env%kpoint_canonical_mu = 0.0_dp

      qs_ot_env%has_complex_kpoint_state = .TRUE.
      qs_ot_env%kpoint_preconditioned_direction_ready = .FALSE.
      qs_ot_env%ener_direction_linear = .FALSE.

   END SUBROUTINE allocate_qs_ot_kpoint_complex_state

! **************************************************************************************************
!> \brief perform cleanup operations (like releasing temporary storage)
!>      at the end of the scf
!> \param scf_env ...
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE scf_env_cleanup(scf_env)
      TYPE(qs_scf_env_type), INTENT(INOUT)               :: scf_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'scf_env_cleanup'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      ! Release SCF work storage
      CALL cp_fm_release(scf_env%scf_work1)

      IF (ASSOCIATED(scf_env%scf_work1_red)) THEN
         CALL cp_fm_release(scf_env%scf_work1_red)
      END IF
      IF (ASSOCIATED(scf_env%scf_work2)) THEN
         CALL cp_fm_release(scf_env%scf_work2)
         DEALLOCATE (scf_env%scf_work2)
         NULLIFY (scf_env%scf_work2)
      END IF
      IF (ASSOCIATED(scf_env%scf_work2_red)) THEN
         CALL cp_fm_release(scf_env%scf_work2_red)
         DEALLOCATE (scf_env%scf_work2_red)
         NULLIFY (scf_env%scf_work2_red)
      END IF
      IF (ASSOCIATED(scf_env%ortho)) THEN
         CALL cp_fm_release(scf_env%ortho)
         DEALLOCATE (scf_env%ortho)
         NULLIFY (scf_env%ortho)
      END IF
      IF (ASSOCIATED(scf_env%ortho_red)) THEN
         CALL cp_fm_release(scf_env%ortho_red)
         DEALLOCATE (scf_env%ortho_red)
         NULLIFY (scf_env%ortho_red)
      END IF
      IF (ASSOCIATED(scf_env%ortho_m1)) THEN
         CALL cp_fm_release(scf_env%ortho_m1)
         DEALLOCATE (scf_env%ortho_m1)
         NULLIFY (scf_env%ortho_m1)
      END IF
      IF (ASSOCIATED(scf_env%ortho_m1_red)) THEN
         CALL cp_fm_release(scf_env%ortho_m1_red)
         DEALLOCATE (scf_env%ortho_m1_red)
         NULLIFY (scf_env%ortho_m1_red)
      END IF

      IF (ASSOCIATED(scf_env%ortho_dbcsr)) THEN
         CALL dbcsr_deallocate_matrix(scf_env%ortho_dbcsr)
      END IF
      IF (ASSOCIATED(scf_env%buf1_dbcsr)) THEN
         CALL dbcsr_deallocate_matrix(scf_env%buf1_dbcsr)
      END IF
      IF (ASSOCIATED(scf_env%buf2_dbcsr)) THEN
         CALL dbcsr_deallocate_matrix(scf_env%buf2_dbcsr)
      END IF

      IF (ASSOCIATED(scf_env%p_mix_new)) THEN
         CALL dbcsr_deallocate_matrix_set(scf_env%p_mix_new)
      END IF

      IF (ASSOCIATED(scf_env%p_delta)) THEN
         CALL dbcsr_deallocate_matrix_set(scf_env%p_delta)
      END IF

      ! Method dependent cleanup
      SELECT CASE (scf_env%method)
      CASE (ot_method_nr)
         !
      CASE (ot_diag_method_nr)
         !
      CASE (general_diag_method_nr)
         !
      CASE (special_diag_method_nr)
         !
      CASE (block_krylov_diag_method_nr)
      CASE (block_davidson_diag_method_nr)
         CALL block_davidson_deallocate(scf_env%block_davidson_env)
      CASE (filter_matrix_diag_method_nr)
         !
      CASE (smeagol_method_nr)
         !
      CASE DEFAULT
         CPABORT("unknown scf method method:"//cp_to_string(scf_env%method))
      END SELECT

      IF (ASSOCIATED(scf_env%outer_scf%variables)) THEN
         DEALLOCATE (scf_env%outer_scf%variables)
      END IF
      IF (ASSOCIATED(scf_env%outer_scf%count)) THEN
         DEALLOCATE (scf_env%outer_scf%count)
      END IF
      IF (ASSOCIATED(scf_env%outer_scf%gradient)) THEN
         DEALLOCATE (scf_env%outer_scf%gradient)
      END IF
      IF (ASSOCIATED(scf_env%outer_scf%energy)) THEN
         DEALLOCATE (scf_env%outer_scf%energy)
      END IF
      IF (ASSOCIATED(scf_env%outer_scf%inv_jacobian) .AND. &
          scf_env%outer_scf%deallocate_jacobian) THEN
         DEALLOCATE (scf_env%outer_scf%inv_jacobian)
      END IF

      CALL timestop(handle)

   END SUBROUTINE scf_env_cleanup

! **************************************************************************************************
!> \brief perform a CDFT scf procedure in the given qs_env
!> \param qs_env the qs_environment where to perform the scf procedure
!> \param should_stop flag determining if calculation should stop
!> \par History
!>      12.2015 Created
!> \author Nico Holmberg
! **************************************************************************************************
   SUBROUTINE cdft_scf(qs_env, should_stop)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(OUT)                               :: should_stop

      CHARACTER(len=*), PARAMETER                        :: routineN = 'cdft_scf'

      INTEGER                                            :: handle, iatom, ispin, ivar, nmo, nvar, &
                                                            output_unit, tsteps
      LOGICAL                                            :: cdft_loop_converged, converged, &
                                                            exit_cdft_loop, first_iteration, &
                                                            my_uocc, uniform_occupation
      REAL(KIND=dp), DIMENSION(:), POINTER               :: mo_occupations
      TYPE(cdft_control_type), POINTER                   :: cdft_control
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s, rho_ao
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(section_vals_type), POINTER                   :: dft_section, input, scf_section

      NULLIFY (scf_env, ks_env, energy, rho, matrix_s, rho_ao, cdft_control, logger, &
               dft_control, pw_env, auxbas_pw_pool, energy, ks_env, scf_env, dft_section, &
               input, scf_section, scf_control, mos, mo_occupations)
      logger => cp_get_default_logger()

      CPASSERT(ASSOCIATED(qs_env))
      CALL get_qs_env(qs_env, scf_env=scf_env, energy=energy, &
                      dft_control=dft_control, scf_control=scf_control, &
                      ks_env=ks_env, input=input)

      CALL timeset(routineN//"_loop", handle)
      dft_section => section_vals_get_subs_vals(input, "DFT")
      scf_section => section_vals_get_subs_vals(dft_section, "SCF")
      output_unit = cp_print_key_unit_nr(logger, scf_section, "PRINT%PROGRAM_RUN_INFO", &
                                         extension=".scfLog")
      first_iteration = .TRUE.

      cdft_control => dft_control%qs_control%cdft_control

      scf_env%outer_scf%iter_count = 0
      cdft_control%total_steps = 0

      ! Write some info about the CDFT calculation
      IF (output_unit > 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,/,T2,A)") &
            "CDFT EXTERNAL SCF WAVEFUNCTION OPTIMIZATION"
         CALL qs_scf_cdft_initial_info(output_unit, cdft_control)
      END IF
      IF (cdft_control%reuse_precond) THEN
         reuse_precond = .FALSE.
         cdft_control%nreused = 0
      END IF
      cdft_outer_loop: DO
         ! Change outer_scf settings to OT settings
         CALL outer_loop_switch(scf_env, scf_control, cdft_control, cdft2ot)
         ! Solve electronic structure with fixed value of constraint
         CALL scf_env_do_scf(scf_env=scf_env, scf_control=scf_control, qs_env=qs_env, &
                             converged=converged, should_stop=should_stop, total_scf_steps=tsteps)
         ! Decide whether to reuse the preconditioner on the next iteration
         IF (cdft_control%reuse_precond) THEN
            ! For convergence in exactly one step, the preconditioner is always reused (assuming max_reuse > 0)
            ! usually this means that the electronic structure has already converged to the correct state
            ! but the constraint optimizer keeps jumping over the optimal solution
            IF (scf_env%outer_scf%iter_count == 1 .AND. scf_env%iter_count == 1 &
                .AND. cdft_control%total_steps /= 1) THEN
               cdft_control%nreused = cdft_control%nreused - 1
            END IF
            ! SCF converged in less than precond_freq steps
            IF (scf_env%outer_scf%iter_count == 1 .AND. scf_env%iter_count <= cdft_control%precond_freq .AND. &
                cdft_control%total_steps /= 1 .AND. cdft_control%nreused < cdft_control%max_reuse) THEN
               reuse_precond = .TRUE.
               cdft_control%nreused = cdft_control%nreused + 1
            ELSE
               reuse_precond = .FALSE.
               cdft_control%nreused = 0
            END IF
         END IF
         ! Update history purging counters
         IF (first_iteration .AND. cdft_control%purge_history) THEN
            cdft_control%istep = cdft_control%istep + 1
            IF (scf_env%outer_scf%iter_count > 1) THEN
               cdft_control%nbad_conv = cdft_control%nbad_conv + 1
               IF (cdft_control%nbad_conv >= cdft_control%purge_freq .AND. &
                   cdft_control%istep >= cdft_control%purge_offset) THEN
                  cdft_control%nbad_conv = 0
                  cdft_control%istep = 0
                  cdft_control%should_purge = .TRUE.
               END IF
            END IF
         END IF
         first_iteration = .FALSE.
         ! Change outer_scf settings to CDFT settings
         CALL outer_loop_switch(scf_env, scf_control, cdft_control, ot2cdft)
         CALL qs_scf_check_outer_exit(qs_env, scf_env, scf_control, should_stop, &
                                      cdft_loop_converged, exit_cdft_loop)
         CALL qs_scf_cdft_info(output_unit, scf_control, scf_env, cdft_control, &
                               energy, cdft_control%total_steps, &
                               should_stop, cdft_loop_converged, cdft_loop=.TRUE.)
         IF (exit_cdft_loop) EXIT cdft_outer_loop
         ! Check if the inverse Jacobian needs to be calculated
         CALL qs_calculate_inverse_jacobian(qs_env)
         ! Check if a line search should be performed to find an optimal step size for the optimizer
         CALL qs_cdft_line_search(qs_env)
         ! Optimize constraint
         CALL outer_loop_optimize(scf_env, scf_control)
         CALL outer_loop_update_qs_env(qs_env, scf_env)
         CALL qs_ks_did_change(ks_env, potential_changed=.TRUE.)
      END DO cdft_outer_loop

      cdft_control%ienergy = cdft_control%ienergy + 1

      ! Store needed arrays for ET coupling calculation
      IF (cdft_control%do_et) THEN
         CALL get_qs_env(qs_env=qs_env, matrix_s=matrix_s, mos=mos)
         nvar = SIZE(cdft_control%target)
         ! Matrix representation of weight function
         ALLOCATE (cdft_control%wmat(nvar))
         DO ivar = 1, nvar
            CALL dbcsr_init_p(cdft_control%wmat(ivar)%matrix)
            CALL dbcsr_copy(cdft_control%wmat(ivar)%matrix, matrix_s(1)%matrix, &
                            name="ET_RESTRAINT_MATRIX")
            CALL dbcsr_set(cdft_control%wmat(ivar)%matrix, 0.0_dp)
            CALL integrate_v_rspace(cdft_control%group(ivar)%weight, &
                                    hmat=cdft_control%wmat(ivar), qs_env=qs_env, &
                                    calculate_forces=.FALSE., &
                                    gapw=dft_control%qs_control%gapw)
         END DO
         ! Overlap matrix
         CALL dbcsr_init_p(cdft_control%matrix_s%matrix)
         CALL dbcsr_copy(cdft_control%matrix_s%matrix, matrix_s(1)%matrix, &
                         name="OVERLAP")
         ! Molecular orbital coefficients
         NULLIFY (cdft_control%mo_coeff)
         ALLOCATE (cdft_control%mo_coeff(dft_control%nspins))
         DO ispin = 1, dft_control%nspins
            CALL cp_fm_create(matrix=cdft_control%mo_coeff(ispin), &
                              matrix_struct=qs_env%mos(ispin)%mo_coeff%matrix_struct, &
                              name="MO_COEFF_A"//TRIM(ADJUSTL(cp_to_string(ispin)))//"MATRIX")
            CALL cp_fm_to_fm(qs_env%mos(ispin)%mo_coeff, &
                             cdft_control%mo_coeff(ispin))
         END DO
         ! Density matrix
         IF (cdft_control%calculate_metric) THEN
            CALL get_qs_env(qs_env, rho=rho)
            CALL qs_rho_get(rho, rho_ao=rho_ao)
            ALLOCATE (cdft_control%matrix_p(dft_control%nspins))
            DO ispin = 1, dft_control%nspins
               NULLIFY (cdft_control%matrix_p(ispin)%matrix)
               CALL dbcsr_init_p(cdft_control%matrix_p(ispin)%matrix)
               CALL dbcsr_copy(cdft_control%matrix_p(ispin)%matrix, rho_ao(ispin)%matrix, &
                               name="DENSITY MATRIX")
            END DO
         END IF
         ! Copy occupation numbers if non-uniform occupation
         uniform_occupation = .TRUE.
         DO ispin = 1, dft_control%nspins
            CALL get_mo_set(mo_set=mos(ispin), uniform_occupation=my_uocc)
            uniform_occupation = uniform_occupation .AND. my_uocc
         END DO
         IF (.NOT. uniform_occupation) THEN
            ALLOCATE (cdft_control%occupations(dft_control%nspins))
            DO ispin = 1, dft_control%nspins
               CALL get_mo_set(mo_set=mos(ispin), &
                               nmo=nmo, &
                               occupation_numbers=mo_occupations)
               ALLOCATE (cdft_control%occupations(ispin)%array(nmo))
               cdft_control%occupations(ispin)%array(1:nmo) = mo_occupations(1:nmo)
            END DO
         END IF
      END IF

      ! Deallocate constraint storage if forces are not needed
      ! In case of a simulation with multiple force_evals,
      ! deallocate only if weight function should not be copied to different force_evals
      IF (.NOT. (cdft_control%save_pot .OR. cdft_control%transfer_pot)) THEN
         CALL get_qs_env(qs_env, pw_env=pw_env)
         CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)
         DO iatom = 1, SIZE(cdft_control%group)
            CALL auxbas_pw_pool%give_back_pw(cdft_control%group(iatom)%weight)
            DEALLOCATE (cdft_control%group(iatom)%weight)
         END DO
         IF (cdft_control%atomic_charges) THEN
            DO iatom = 1, cdft_control%natoms
               CALL auxbas_pw_pool%give_back_pw(cdft_control%charge(iatom))
            END DO
            DEALLOCATE (cdft_control%charge)
         END IF
         IF (cdft_control%type == outer_scf_becke_constraint .AND. &
             cdft_control%becke_control%cavity_confine) THEN
            IF (.NOT. ASSOCIATED(cdft_control%becke_control%cavity_mat)) THEN
               CALL auxbas_pw_pool%give_back_pw(cdft_control%becke_control%cavity)
            ELSE
               DEALLOCATE (cdft_control%becke_control%cavity_mat)
            END IF
         ELSE IF (cdft_control%type == outer_scf_hirshfeld_constraint) THEN
            IF (ASSOCIATED(cdft_control%hirshfeld_control%hirshfeld_env%fnorm)) THEN
               CALL auxbas_pw_pool%give_back_pw(cdft_control%hirshfeld_control%hirshfeld_env%fnorm)
            END IF
         END IF
         IF (ASSOCIATED(cdft_control%charges_fragment)) DEALLOCATE (cdft_control%charges_fragment)
         cdft_control%need_pot = .TRUE.
         cdft_control%external_control = .FALSE.
      END IF

      CALL timestop(handle)

   END SUBROUTINE cdft_scf

! **************************************************************************************************
!> \brief perform cleanup operations for cdft_control
!> \param cdft_control container for the external CDFT SCF loop variables
!> \par History
!>      12.2015 created [Nico Holmberg]
!> \author Nico Holmberg
! **************************************************************************************************
   SUBROUTINE cdft_control_cleanup(cdft_control)
      TYPE(cdft_control_type), POINTER                   :: cdft_control

      IF (ASSOCIATED(cdft_control%constraint%variables)) THEN
         DEALLOCATE (cdft_control%constraint%variables)
      END IF
      IF (ASSOCIATED(cdft_control%constraint%count)) THEN
         DEALLOCATE (cdft_control%constraint%count)
      END IF
      IF (ASSOCIATED(cdft_control%constraint%gradient)) THEN
         DEALLOCATE (cdft_control%constraint%gradient)
      END IF
      IF (ASSOCIATED(cdft_control%constraint%energy)) THEN
         DEALLOCATE (cdft_control%constraint%energy)
      END IF
      IF (ASSOCIATED(cdft_control%constraint%inv_jacobian) .AND. &
          cdft_control%constraint%deallocate_jacobian) THEN
         DEALLOCATE (cdft_control%constraint%inv_jacobian)
      END IF

   END SUBROUTINE cdft_control_cleanup

! **************************************************************************************************
!> \brief Calculates the finite difference inverse Jacobian
!> \param qs_env the qs_environment_type where to compute the Jacobian
!> \par History
!>      01.2017 created [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE qs_calculate_inverse_jacobian(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_calculate_inverse_jacobian'

      CHARACTER(len=default_path_length)                 :: project_name
      INTEGER                                            :: counter, handle, i, ispin, iter_count, &
                                                            iwork, j, max_scf, nspins, nsteps, &
                                                            nvar, nwork, output_unit, pwork, &
                                                            tsteps, twork
      LOGICAL                                            :: converged, explicit_jacobian, &
                                                            should_build, should_stop, &
                                                            use_md_history
      REAL(KIND=dp)                                      :: inv_error, step_size
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: coeff, dh, step_multiplier
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: jacobian
      REAL(KIND=dp), DIMENSION(:), POINTER               :: energy
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: gradient, inv_jacobian
      TYPE(cdft_control_type), POINTER                   :: cdft_control
      TYPE(cp_logger_type), POINTER                      :: logger, tmp_logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: p_rmpv
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rho_ao_kp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_type), ALLOCATABLE, DIMENSION(:)       :: mos_stashed
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(qs_energy_type), POINTER                      :: energy_qs
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(scf_control_type), POINTER                    :: scf_control

      NULLIFY (energy, gradient, p_rmpv, rho_ao_kp, mos, rho, &
               ks_env, scf_env, scf_control, dft_control, cdft_control, &
               inv_jacobian, para_env, tmp_logger, energy_qs)
      logger => cp_get_default_logger()

      CPASSERT(ASSOCIATED(qs_env))
      CALL get_qs_env(qs_env, scf_env=scf_env, ks_env=ks_env, &
                      scf_control=scf_control, mos=mos, rho=rho, &
                      dft_control=dft_control, &
                      para_env=para_env, energy=energy_qs)
      explicit_jacobian = .FALSE.
      should_build = .FALSE.
      use_md_history = .FALSE.
      iter_count = scf_env%outer_scf%iter_count
      ! Quick exit if optimizer does not require Jacobian
      IF (.NOT. ASSOCIATED(scf_control%outer_scf%cdft_opt_control)) RETURN
      ! Check if Jacobian should be calculated and initialize
      CALL timeset(routineN, handle)
      CALL initialize_inverse_jacobian(scf_control, scf_env, explicit_jacobian, should_build, used_history)
      IF (scf_control%outer_scf%cdft_opt_control%jacobian_restart) THEN
         ! Restart from previously calculated inverse Jacobian
         should_build = .FALSE.
         CALL restart_inverse_jacobian(qs_env)
      END IF
      IF (should_build) THEN
         scf_env%outer_scf%deallocate_jacobian = .FALSE.
         ! Actually need to (re)build the Jacobian
         IF (explicit_jacobian) THEN
            ! Build Jacobian with finite differences
            cdft_control => dft_control%qs_control%cdft_control
            IF (.NOT. ASSOCIATED(cdft_control)) THEN
               CALL cp_abort(__LOCATION__, &
                             "Optimizers that need the explicit Jacobian can"// &
                             " only be used together with a valid CDFT constraint.")
            END IF
            ! Redirect output from Jacobian calculation to a new file by creating a temporary logger
            project_name = logger%iter_info%project_name
            CALL create_tmp_logger(para_env, project_name, "-JacobianInfo.out", output_unit, tmp_logger)
            ! Save last converged state so we can roll back to it (mo_coeff and some outer_loop variables)
            nspins = dft_control%nspins
            ALLOCATE (mos_stashed(nspins))
            DO ispin = 1, nspins
               CALL duplicate_mo_set(mos_stashed(ispin), mos(ispin))
            END DO
            CALL qs_rho_get(rho, rho_ao_kp=rho_ao_kp)
            p_rmpv => rho_ao_kp(:, 1)
            ! Allocate work
            nvar = SIZE(scf_env%outer_scf%variables, 1)
            max_scf = scf_control%outer_scf%max_scf + 1
            ALLOCATE (gradient(nvar, max_scf))
            gradient = scf_env%outer_scf%gradient
            ALLOCATE (energy(max_scf))
            energy = scf_env%outer_scf%energy
            ALLOCATE (jacobian(nvar, nvar))
            jacobian = 0.0_dp
            nsteps = cdft_control%total_steps
            ! Setup finite difference scheme
            CALL prepare_jacobian_stencil(qs_env, output_unit, nwork, pwork, coeff, step_multiplier, dh)
            twork = pwork - nwork
            DO i = 1, nvar
               jacobian(i, :) = coeff(0)*scf_env%outer_scf%gradient(i, iter_count)
            END DO
            ! Calculate the Jacobian by perturbing each Lagrangian and recalculating the energy self-consistently
            CALL cp_add_default_logger(tmp_logger)
            DO i = 1, nvar
               IF (output_unit > 0) THEN
                  WRITE (output_unit, FMT="(A)") " "
                  WRITE (output_unit, FMT="(A)") " #####################################"
                  WRITE (output_unit, '(A,I3,A,I3,A)') &
                     " ###  Constraint        ", i, " of ", nvar, " ###"
                  WRITE (output_unit, FMT="(A)") " #####################################"
               END IF
               counter = 0
               DO iwork = nwork, pwork
                  IF (iwork == 0) CYCLE
                  counter = counter + 1
                  IF (output_unit > 0) THEN
                     WRITE (output_unit, FMT="(A)") " #####################################"
                     WRITE (output_unit, '(A,I3,A,I3,A)') &
                        " ###  Energy evaluation ", counter, " of ", twork, " ###"
                     WRITE (output_unit, FMT="(A)") " #####################################"
                  END IF
                  IF (SIZE(scf_control%outer_scf%cdft_opt_control%jacobian_step) == 1) THEN
                     step_size = scf_control%outer_scf%cdft_opt_control%jacobian_step(1)
                  ELSE
                     step_size = scf_control%outer_scf%cdft_opt_control%jacobian_step(i)
                  END IF
                  scf_env%outer_scf%variables(:, iter_count + 1) = scf_env%outer_scf%variables(:, iter_count)
                  scf_env%outer_scf%variables(i, iter_count + 1) = scf_env%outer_scf%variables(i, iter_count) + &
                                                                   step_multiplier(iwork)*step_size
                  CALL outer_loop_update_qs_env(qs_env, scf_env)
                  CALL qs_ks_did_change(ks_env, potential_changed=.TRUE.)
                  CALL outer_loop_switch(scf_env, scf_control, cdft_control, cdft2ot)
                  CALL scf_env_do_scf(scf_env=scf_env, scf_control=scf_control, qs_env=qs_env, &
                                      converged=converged, should_stop=should_stop, total_scf_steps=tsteps)
                  CALL outer_loop_switch(scf_env, scf_control, cdft_control, ot2cdft)
                  ! Update (iter_count + 1) element of gradient and print constraint info
                  scf_env%outer_scf%iter_count = scf_env%outer_scf%iter_count + 1
                  CALL outer_loop_gradient(qs_env, scf_env)
                  CALL qs_scf_cdft_info(output_unit, scf_control, scf_env, cdft_control, &
                                        energy_qs, cdft_control%total_steps, &
                                        should_stop=.FALSE., outer_loop_converged=.FALSE., cdft_loop=.FALSE.)
                  scf_env%outer_scf%iter_count = scf_env%outer_scf%iter_count - 1
                  ! Update Jacobian
                  DO j = 1, nvar
                     jacobian(j, i) = jacobian(j, i) + coeff(iwork)*scf_env%outer_scf%gradient(j, iter_count + 1)
                  END DO
                  ! Reset everything to last converged state
                  scf_env%outer_scf%variables(:, iter_count + 1) = 0.0_dp
                  scf_env%outer_scf%gradient = gradient
                  scf_env%outer_scf%energy = energy
                  cdft_control%total_steps = nsteps
                  DO ispin = 1, nspins
                     CALL deallocate_mo_set(mos(ispin))
                     CALL duplicate_mo_set(mos(ispin), mos_stashed(ispin))
                     CALL calculate_density_matrix(mos(ispin), &
                                                   p_rmpv(ispin)%matrix)
                  END DO
                  CALL qs_rho_update_rho(rho, qs_env=qs_env)
                  CALL qs_ks_did_change(qs_env%ks_env, rho_changed=.TRUE.)
               END DO
            END DO
            CALL cp_rm_default_logger()
            CALL cp_logger_release(tmp_logger)
            ! Finalize and invert Jacobian
            DO j = 1, nvar
               DO i = 1, nvar
                  jacobian(i, j) = jacobian(i, j)/dh(j)
               END DO
            END DO
            IF (.NOT. ASSOCIATED(scf_env%outer_scf%inv_jacobian)) THEN
               ALLOCATE (scf_env%outer_scf%inv_jacobian(nvar, nvar))
            END IF
            inv_jacobian => scf_env%outer_scf%inv_jacobian
            CALL invert_matrix(jacobian, inv_jacobian, inv_error)
            scf_control%outer_scf%cdft_opt_control%broyden_update = .FALSE.
            ! Release temporary storage
            DO ispin = 1, nspins
               CALL deallocate_mo_set(mos_stashed(ispin))
            END DO
            DEALLOCATE (mos_stashed, jacobian, gradient, energy, coeff, step_multiplier, dh)
            IF (output_unit > 0) THEN
               WRITE (output_unit, FMT="(/,A)") &
                  " ================================== JACOBIAN CALCULATED =================================="
               CALL close_file(unit_number=output_unit)
            END IF
         ELSE
            ! Build a strictly diagonal Jacobian from history and invert it
            CALL build_diagonal_jacobian(qs_env, used_history)
         END IF
      END IF
      IF (ASSOCIATED(scf_env%outer_scf%inv_jacobian) .AND. para_env%is_source()) THEN
         ! Write restart file for inverse Jacobian
         CALL print_inverse_jacobian(logger, scf_env%outer_scf%inv_jacobian, iter_count)
      END IF
      ! Update counter
      scf_control%outer_scf%cdft_opt_control%ijacobian(1) = scf_control%outer_scf%cdft_opt_control%ijacobian(1) + 1
      CALL timestop(handle)

   END SUBROUTINE qs_calculate_inverse_jacobian

! **************************************************************************************************
!> \brief Perform backtracking line search to find the optimal step size for the CDFT constraint
!>        optimizer. Assumes that the CDFT gradient function is a smooth function of the constraint
!>        variables.
!> \param qs_env the qs_environment_type where to perform the line search
!> \par History
!>      02.2017 created [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE qs_cdft_line_search(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_cdft_line_search'

      CHARACTER(len=default_path_length)                 :: project_name
      INTEGER                                            :: handle, i, ispin, iter_count, &
                                                            max_linesearch, max_scf, nspins, &
                                                            nsteps, nvar, output_unit, tsteps
      LOGICAL :: continue_ls, continue_ls_exit, converged, do_linesearch, found_solution, &
         reached_maxls, should_exit, should_stop, sign_changed
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: positive_sign
      REAL(KIND=dp)                                      :: alpha, alpha_ls, factor, norm_ls
      REAL(KIND=dp), DIMENSION(:), POINTER               :: energy
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: gradient, inv_jacobian
      REAL(KIND=dp), EXTERNAL                            :: dnrm2
      TYPE(cdft_control_type), POINTER                   :: cdft_control
      TYPE(cp_logger_type), POINTER                      :: logger, tmp_logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: p_rmpv
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rho_ao_kp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(qs_energy_type), POINTER                      :: energy_qs
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(scf_control_type), POINTER                    :: scf_control

      CALL timeset(routineN, handle)

      NULLIFY (energy, gradient, p_rmpv, rho_ao_kp, mos, rho, &
               ks_env, scf_env, scf_control, dft_control, &
               cdft_control, inv_jacobian, para_env, &
               tmp_logger, energy_qs)
      logger => cp_get_default_logger()

      CPASSERT(ASSOCIATED(qs_env))
      CALL get_qs_env(qs_env, scf_env=scf_env, ks_env=ks_env, &
                      scf_control=scf_control, mos=mos, rho=rho, &
                      dft_control=dft_control, &
                      para_env=para_env, energy=energy_qs)
      do_linesearch = .FALSE.
      SELECT CASE (scf_control%outer_scf%optimizer)
      CASE DEFAULT
         do_linesearch = .FALSE.
      CASE (outer_scf_optimizer_newton_ls)
         do_linesearch = .TRUE.
      CASE (outer_scf_optimizer_broyden)
         SELECT CASE (scf_control%outer_scf%cdft_opt_control%broyden_type)
         CASE (broyden_type_1, broyden_type_2, broyden_type_1_explicit, broyden_type_2_explicit)
            do_linesearch = .FALSE.
         CASE (broyden_type_1_ls, broyden_type_1_explicit_ls, broyden_type_2_ls, broyden_type_2_explicit_ls)
            cdft_control => dft_control%qs_control%cdft_control
            IF (.NOT. ASSOCIATED(cdft_control)) THEN
               CALL cp_abort(__LOCATION__, &
                             "Optimizers that perform a line search can"// &
                             " only be used together with a valid CDFT constraint")
            END IF
            IF (ASSOCIATED(scf_env%outer_scf%inv_jacobian)) THEN
               do_linesearch = .TRUE.
            END IF
         END SELECT
      END SELECT
      IF (do_linesearch) THEN
         BLOCK
            TYPE(mo_set_type), DIMENSION(:), ALLOCATABLE :: mos_ls, mos_stashed
            cdft_control => dft_control%qs_control%cdft_control
            IF (.NOT. ASSOCIATED(cdft_control)) THEN
               CALL cp_abort(__LOCATION__, &
                             "Optimizers that perform a line search can"// &
                             " only be used together with a valid CDFT constraint")
            END IF
            CPASSERT(ASSOCIATED(scf_env%outer_scf%inv_jacobian))
            CPASSERT(ASSOCIATED(scf_control%outer_scf%cdft_opt_control))
            alpha = scf_control%outer_scf%cdft_opt_control%newton_step_save
            iter_count = scf_env%outer_scf%iter_count
            ! Redirect output from line search procedure to a new file by creating a temporary logger
            project_name = logger%iter_info%project_name
            CALL create_tmp_logger(para_env, project_name, "-LineSearch.out", output_unit, tmp_logger)
            ! Save last converged state so we can roll back to it (mo_coeff and some outer_loop variables)
            nspins = dft_control%nspins
            ALLOCATE (mos_stashed(nspins))
            DO ispin = 1, nspins
               CALL duplicate_mo_set(mos_stashed(ispin), mos(ispin))
            END DO
            CALL qs_rho_get(rho, rho_ao_kp=rho_ao_kp)
            p_rmpv => rho_ao_kp(:, 1)
            nsteps = cdft_control%total_steps
            ! Allocate work
            nvar = SIZE(scf_env%outer_scf%variables, 1)
            max_scf = scf_control%outer_scf%max_scf + 1
            max_linesearch = scf_control%outer_scf%cdft_opt_control%max_ls
            continue_ls = scf_control%outer_scf%cdft_opt_control%continue_ls
            factor = scf_control%outer_scf%cdft_opt_control%factor_ls
            continue_ls_exit = .FALSE.
            found_solution = .FALSE.
            ALLOCATE (gradient(nvar, max_scf))
            gradient = scf_env%outer_scf%gradient
            ALLOCATE (energy(max_scf))
            energy = scf_env%outer_scf%energy
            reached_maxls = .FALSE.
            ! Broyden optimizers: perform update of inv_jacobian if necessary
            IF (scf_control%outer_scf%cdft_opt_control%broyden_update) THEN
               CALL outer_loop_optimize(scf_env, scf_control)
               ! Reset the variables and prevent a reupdate of inv_jacobian
               scf_env%outer_scf%variables(:, iter_count + 1) = 0
               scf_control%outer_scf%cdft_opt_control%broyden_update = .FALSE.
            END IF
            ! Print some info
            IF (output_unit > 0) THEN
               WRITE (output_unit, FMT="(/,A)") &
                  " ================================== LINE SEARCH STARTED  =================================="
               WRITE (output_unit, FMT="(A,I5,A)") &
                  " Evaluating optimal step size for optimizer using a maximum of", max_linesearch, " steps"
               IF (continue_ls) THEN
                  WRITE (output_unit, FMT="(A)") &
                     " Line search continues until best step size is found or max steps are reached"
               END IF
               WRITE (output_unit, '(/,A,F5.3)') &
                  " Initial step size: ", alpha
               WRITE (output_unit, '(/,A,F5.3)') &
                  " Step size update factor: ", factor
               WRITE (output_unit, '(/,A,I10,A,I10)') &
                  " Energy evaluation: ", cdft_control%ienergy, ", CDFT SCF iteration: ", iter_count
            END IF
            ! Perform backtracking line search
            CALL cp_add_default_logger(tmp_logger)
            DO i = 1, max_linesearch
               IF (output_unit > 0) THEN
                  WRITE (output_unit, FMT="(A)") " "
                  WRITE (output_unit, FMT="(A)") " #####################################"
                  WRITE (output_unit, '(A,I10,A)') &
                     " ###  Line search step: ", i, " ###"
                  WRITE (output_unit, FMT="(A)") " #####################################"
               END IF
               inv_jacobian => scf_env%outer_scf%inv_jacobian
               ! Newton update of CDFT variables with a step size of alpha
               scf_env%outer_scf%variables(:, iter_count + 1) = scf_env%outer_scf%variables(:, iter_count) - alpha* &
                                                                MATMUL(inv_jacobian, scf_env%outer_scf%gradient(:, iter_count))
               ! With updated CDFT variables, perform SCF
               CALL outer_loop_update_qs_env(qs_env, scf_env)
               CALL qs_ks_did_change(ks_env, potential_changed=.TRUE.)
               CALL outer_loop_switch(scf_env, scf_control, cdft_control, cdft2ot)
               CALL scf_env_do_scf(scf_env=scf_env, scf_control=scf_control, qs_env=qs_env, &
                                   converged=converged, should_stop=should_stop, total_scf_steps=tsteps)
               CALL outer_loop_switch(scf_env, scf_control, cdft_control, ot2cdft)
               ! Update (iter_count + 1) element of gradient and print constraint info
               scf_env%outer_scf%iter_count = scf_env%outer_scf%iter_count + 1
               CALL outer_loop_gradient(qs_env, scf_env)
               CALL qs_scf_cdft_info(output_unit, scf_control, scf_env, cdft_control, &
                                     energy_qs, cdft_control%total_steps, &
                                     should_stop=.FALSE., outer_loop_converged=.FALSE., cdft_loop=.FALSE.)
               scf_env%outer_scf%iter_count = scf_env%outer_scf%iter_count - 1
               ! Store sign of initial gradient for each variable for continue_ls
               IF (continue_ls .AND. .NOT. ALLOCATED(positive_sign)) THEN
                  ALLOCATE (positive_sign(nvar))
                  DO ispin = 1, nvar
                     positive_sign(ispin) = scf_env%outer_scf%gradient(ispin, iter_count + 1) >= 0.0_dp
                  END DO
               END IF
               ! Check if the L2 norm of the gradient decreased
               inv_jacobian => scf_env%outer_scf%inv_jacobian
               IF (dnrm2(nvar, scf_env%outer_scf%gradient(:, iter_count + 1), 1) < &
                   dnrm2(nvar, scf_env%outer_scf%gradient(:, iter_count), 1)) THEN
                  ! Optimal step size found
                  IF (.NOT. continue_ls) THEN
                     should_exit = .TRUE.
                  ELSE
                     ! But line search continues for at least one more iteration in an attempt to find a better solution
                     ! if max number of steps is not exceeded
                     IF (found_solution) THEN
                        ! Check if the norm also decreased w.r.t. to previously found solution
                        IF (dnrm2(nvar, scf_env%outer_scf%gradient(:, iter_count + 1), 1) > norm_ls) THEN
                           ! Norm increased => accept previous solution and exit
                           continue_ls_exit = .TRUE.
                        END IF
                     END IF
                     ! Store current state and the value of alpha
                     IF (.NOT. continue_ls_exit) THEN
                        should_exit = .FALSE.
                        alpha_ls = alpha
                        found_solution = .TRUE.
                        norm_ls = dnrm2(nvar, scf_env%outer_scf%gradient(:, iter_count + 1), 1)
                        ! Check if the sign of the gradient has changed for all variables (w.r.t initial gradient)
                        ! In this case we should exit because further line search steps will just increase the norm
                        sign_changed = .TRUE.
                        DO ispin = 1, nvar
                           sign_changed = sign_changed .AND. (positive_sign(ispin) .NEQV. &
                                                              scf_env%outer_scf%gradient(ispin, iter_count + 1) >= 0.0_dp)
                        END DO
                        IF (.NOT. ALLOCATED(mos_ls)) THEN
                           ALLOCATE (mos_ls(nspins))
                        ELSE
                           DO ispin = 1, nspins
                              CALL deallocate_mo_set(mos_ls(ispin))
                           END DO
                        END IF
                        DO ispin = 1, nspins
                           CALL duplicate_mo_set(mos_ls(ispin), mos(ispin))
                        END DO
                        alpha = alpha*factor
                        ! Exit on last iteration
                        IF (i == max_linesearch) continue_ls_exit = .TRUE.
                        ! Exit if constraint target is satisfied to requested tolerance
                        IF (SQRT(MAXVAL(scf_env%outer_scf%gradient(:, scf_env%outer_scf%iter_count + 1)**2)) < &
                            scf_control%outer_scf%eps_scf) THEN
                           continue_ls_exit = .TRUE.
                        END IF
                        ! Exit if line search jumped over the optimal step length
                        IF (sign_changed) continue_ls_exit = .TRUE.
                     END IF
                  END IF
               ELSE
                  ! Gradient increased => alpha is too large (if the gradient function is smooth)
                  should_exit = .FALSE.
                  ! Update alpha using Armijo's scheme
                  alpha = alpha*factor
               END IF
               IF (continue_ls_exit) THEN
                  ! Continuation of line search did not yield a better alpha, use previously located solution and exit
                  alpha = alpha_ls
                  DO ispin = 1, nspins
                     CALL deallocate_mo_set(mos(ispin))
                     CALL duplicate_mo_set(mos(ispin), mos_ls(ispin))
                     CALL calculate_density_matrix(mos(ispin), &
                                                   p_rmpv(ispin)%matrix)
                     CALL deallocate_mo_set(mos_ls(ispin))
                  END DO
                  CALL qs_rho_update_rho(rho, qs_env=qs_env)
                  CALL qs_ks_did_change(qs_env%ks_env, rho_changed=.TRUE.)
                  DEALLOCATE (mos_ls)
                  should_exit = .TRUE.
               END IF
               ! Reached max steps and SCF converged: continue with last iterated step size
               IF (.NOT. should_exit .AND. &
                   (i == max_linesearch .AND. converged .AND. .NOT. found_solution)) THEN
                  should_exit = .TRUE.
                  reached_maxls = .TRUE.
                  alpha = alpha*(1.0_dp/factor)
               END IF
               ! Reset outer SCF environment to last converged state
               scf_env%outer_scf%variables(:, iter_count + 1) = 0.0_dp
               scf_env%outer_scf%gradient = gradient
               scf_env%outer_scf%energy = energy
               ! Exit line search if a suitable step size was found
               IF (should_exit) EXIT
               ! Reset the electronic structure
               cdft_control%total_steps = nsteps
               DO ispin = 1, nspins
                  CALL deallocate_mo_set(mos(ispin))
                  CALL duplicate_mo_set(mos(ispin), mos_stashed(ispin))
                  CALL calculate_density_matrix(mos(ispin), &
                                                p_rmpv(ispin)%matrix)
               END DO
               CALL qs_rho_update_rho(rho, qs_env=qs_env)
               CALL qs_ks_did_change(qs_env%ks_env, rho_changed=.TRUE.)
            END DO
            scf_control%outer_scf%cdft_opt_control%newton_step = alpha
            IF (.NOT. should_exit) THEN
               CALL cp_warn(__LOCATION__, &
                            "Line search did not converge. CDFT SCF proceeds with fixed step size.")
               scf_control%outer_scf%cdft_opt_control%newton_step = scf_control%outer_scf%cdft_opt_control%newton_step_save
            END IF
            IF (reached_maxls) THEN
               CALL cp_warn(__LOCATION__, &
                            "Line search did not converge. CDFT SCF proceeds with lasted iterated step size.")
            END IF
            CALL cp_rm_default_logger()
            CALL cp_logger_release(tmp_logger)
            ! Release temporary storage
            DO ispin = 1, nspins
               CALL deallocate_mo_set(mos_stashed(ispin))
            END DO
            DEALLOCATE (mos_stashed, gradient, energy)
            IF (ALLOCATED(positive_sign)) DEALLOCATE (positive_sign)
            IF (output_unit > 0) THEN
               WRITE (output_unit, FMT="(/,A)") &
                  " ================================== LINE SEARCH COMPLETE =================================="
               CALL close_file(unit_number=output_unit)
            END IF
         END BLOCK
      END IF

      CALL timestop(handle)

   END SUBROUTINE qs_cdft_line_search

END MODULE qs_scf
