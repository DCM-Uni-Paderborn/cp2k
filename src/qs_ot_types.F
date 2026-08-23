!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief orbital transformations
!> \par History
!>      Added Taylor expansion based computation of the matrix functions (01.2004)
!>      added additional rotation variables for non-equivalent occupied orbs (08.2004)
!> \author Joost VandeVondele (06.2002)
! **************************************************************************************************
MODULE qs_ot_types
   USE bibliography,                    ONLY: VandeVondele2003,&
                                              Weber2008,&
                                              cite_reference
   USE cp_blacs_env,                    ONLY: cp_blacs_env_release,&
                                              cp_blacs_env_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_init_p,&
                                              dbcsr_p_type,&
                                              dbcsr_release_p,&
                                              dbcsr_set,&
                                              dbcsr_type,&
                                              dbcsr_type_no_symmetry
   USE cp_dbcsr_contrib,                ONLY: dbcsr_add_on_diag
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_m_by_n_from_row_template,&
                                              cp_dbcsr_m_by_n_from_template,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_fm_struct,                    ONLY: cp_fm_struct_get,&
                                              cp_fm_struct_type
   USE input_constants,                 ONLY: &
        cholesky_reduce, ls_2pnt, ls_3pnt, ls_adapt, ls_gold, ls_none, ot_algo_irac, &
        ot_algo_taylor_or_diag, ot_chol_irac, ot_lwdn_irac, ot_mini_broyden, ot_mini_cg, &
        ot_mini_diis, ot_mini_sd, ot_poly_irac, ot_precond_fermi_low_rank, ot_precond_full_all, &
        ot_precond_full_kinetic, ot_precond_full_single, ot_precond_full_single_inverse, &
        ot_precond_none, ot_precond_s_inverse, ot_precond_solver_chebyshev, &
        ot_precond_solver_default, ot_precond_solver_direct, ot_precond_solver_inv_chol, &
        ot_precond_solver_update
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_release,&
                                              mp_para_env_type
   USE preconditioner_types,            ONLY: preconditioner_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   PUBLIC  :: qs_ot_type
   PUBLIC  :: qs_ot_settings_type
   PUBLIC  :: qs_ot_destroy
   PUBLIC  :: qs_ot_allocate
   PUBLIC  :: qs_ot_channel_index
   PUBLIC  :: qs_ot_check_channel_context
   PUBLIC  :: qs_ot_init
   PUBLIC  :: qs_ot_reset_accepted_model_calibration
   PUBLIC  :: qs_ot_number_of_channels
   PUBLIC  :: qs_ot_occupation_preconditioner_active
   PUBLIC  :: qs_ot_kpoint_preconditioner_supported
   PUBLIC  :: qs_ot_kpoint_preconditioner_solver_supported
   PUBLIC  :: ot_metallic_variable_defaults
   PUBLIC  :: qs_ot_restart_signature_nchar
   PUBLIC  :: qs_ot_restart_signature_nint
   PUBLIC  :: qs_ot_restart_signature_nlog
   PUBLIC  :: qs_ot_restart_signature_nreal
   PUBLIC  :: qs_ot_settings_restart_compatible
   PUBLIC  :: qs_ot_settings_restart_signature
   PUBLIC  :: qs_ot_settings_init
   PUBLIC  :: qs_ot_set_context
   PUBLIC  :: ot_readwrite_input

! **************************************************************************************************
!> \brief notice, this variable needs to be copyable, needed for spins as e.g. in qs_ot_scf
! **************************************************************************************************
   TYPE qs_ot_settings_type
      LOGICAL           :: do_rotation = .FALSE., do_ener = .FALSE.
      LOGICAL           :: ks = .FALSE.
      CHARACTER(LEN=4)  :: ot_method = ""
      CHARACTER(LEN=3)  :: ot_algorithm = ""
      CHARACTER(LEN=4)  :: line_search_method = ""
      CHARACTER(LEN=20) :: preconditioner_name = ""
      INTEGER           :: preconditioner_type = -1
      INTEGER           :: cholesky_type = -1
      INTEGER           :: ot_state = 0
      CHARACTER(LEN=20) :: precond_solver_name = ""
      INTEGER           :: precond_solver_type = -1
      INTEGER           :: chebyshev_degree = 8
      LOGICAL           :: safer_diis = .FALSE.
      REAL(KIND=dp)     :: ds_min = -1.0_dp
      REAL(KIND=dp)     :: energy_gap = -1.0_dp
      INTEGER           :: diis_m = -1
      INTEGER           :: max_scf_diis = 0
      REAL(KIND=dp)     :: gold_target = -1.0_dp
      REAL(KIND=dp)     :: eps_taylor = -1.0_dp ! minimum accuracy of Taylor expansion
      INTEGER           :: max_taylor = -1 ! maximum order of Taylor expansion before switching to diagonalization
      INTEGER           :: irac_degree = -1 ! used to control the refinement polynomial degree
      INTEGER           :: max_irac = -1 ! maximum number of iteration for refinement
      REAL(KIND=dp)     :: eps_irac = -1.0_dp ! target accuracy for refinement
      REAL(KIND=dp)     :: eps_irac_quick_exit = -1.0_dp
      REAL(KIND=dp)     :: eps_irac_filter_matrix = -1.0_dp
      REAL(KIND=dp)     :: eps_irac_switch = -1.0_dp
      LOGICAL           :: on_the_fly_loc = .FALSE.
      CHARACTER(LEN=4)  :: ortho_irac = ""
      LOGICAL           :: occupation_preconditioner = .FALSE., add_nondiag_energy = .FALSE.
      LOGICAL           :: canonical_tether = .TRUE., canonical_tether_adaptive = .TRUE.
      REAL(KIND=dp)     :: nondiag_energy_strength = -1.0_dp
      REAL(KIND=dp)     :: canonical_tether_strength = 1.0E-3_dp
      REAL(KIND=dp)     :: canonical_tether_max_strength = 1.5E-3_dp
      REAL(KIND=dp)     :: broyden_beta = -1.0_dp, broyden_gamma = -1.0_dp, broyden_sigma = -1.0_dp
      REAL(KIND=dp)     :: broyden_eta = -1.0_dp, broyden_omega = -1.0_dp, broyden_sigma_decrease = -1.0_dp
      REAL(KIND=dp)     :: broyden_sigma_min = -1.0_dp
      LOGICAL           :: broyden_forget_history = .FALSE., broyden_adaptive_sigma = .FALSE.
      LOGICAL           :: broyden_enable_flip = .FALSE.
   END TYPE qs_ot_settings_type

! **************************************************************************************************
   TYPE qs_ot_type
      ! this sets the method to be used
      TYPE(qs_ot_settings_type) :: settings = qs_ot_settings_type()
      LOGICAL                   :: restricted = .FALSE.
      INTEGER                   :: spin_index = 0, kpoint_index = 0, local_kpoint_index = 0
      INTEGER                   :: occupation_group_index = 1
      INTEGER                   :: kpoint_multiplicity = 1
      REAL(KIND=dp)             :: kpoint_weight = 1.0_dp
      LOGICAL                   :: has_kpoint_context = .FALSE.
      LOGICAL                   :: state_allocated = .FALSE.
      LOGICAL                   :: has_complex_kpoint_state = .FALSE.

      ! first part of the variables, for occupied subspace invariant optimisation

      ! add a preconditioner matrix. should be symmetric and positive definite
      ! the type of this matrix might change in the future
      TYPE(preconditioner_type), POINTER :: preconditioner => NULL()

      ! these will/might change during iterations

      ! OT / TOD
      TYPE(dbcsr_type), POINTER :: matrix_p => NULL()
      TYPE(dbcsr_type), POINTER :: matrix_r => NULL()
      TYPE(dbcsr_type), POINTER :: matrix_sinp => NULL()
      TYPE(dbcsr_type), POINTER :: matrix_cosp => NULL()
      TYPE(dbcsr_type), POINTER :: matrix_sinp_b => NULL()
      TYPE(dbcsr_type), POINTER :: matrix_cosp_b => NULL()
      TYPE(dbcsr_type), POINTER :: matrix_buf1 => NULL()
      TYPE(dbcsr_type), POINTER :: matrix_buf2 => NULL()
      TYPE(dbcsr_type), POINTER :: matrix_buf3 => NULL()
      TYPE(dbcsr_type), POINTER :: matrix_buf4 => NULL()
      TYPE(dbcsr_type), POINTER :: matrix_os => NULL()
      TYPE(dbcsr_type), POINTER :: matrix_buf1_ortho => NULL()
      TYPE(dbcsr_type), POINTER :: matrix_buf2_ortho => NULL()

      REAL(KIND=dp), DIMENSION(:), POINTER :: evals => NULL()
      REAL(KIND=dp), DIMENSION(:), POINTER :: dum => NULL()

      ! matrix os valid
      LOGICAL :: os_valid = .FALSE.

      ! for efficient/parallel writing to the blacs_matrix
      TYPE(mp_para_env_type), POINTER :: para_env => NULL()
      TYPE(cp_blacs_env_type), POINTER :: blacs_env => NULL()

      ! mo-like vectors
      TYPE(dbcsr_type), POINTER :: matrix_c0 => NULL(), matrix_sc0 => NULL(), matrix_psc0 => NULL()
      TYPE(dbcsr_type), POINTER :: matrix_c0_im => NULL(), matrix_sc0_im => NULL()

      ! OT / IR
      TYPE(dbcsr_type), POINTER :: buf1_k_k_nosym => NULL(), buf2_k_k_nosym => NULL(), &
                                   buf3_k_k_nosym => NULL(), buf4_k_k_nosym => NULL(), &
                                   buf1_k_k_sym => NULL(), buf2_k_k_sym => NULL(), &
                                   buf3_k_k_sym => NULL(), buf4_k_k_sym => NULL(), &
                                   p_k_k_sym => NULL(), buf1_n_k => NULL(), buf1_n_k_dp => NULL()

      ! only here for the ease of programming. These will have to be supplied
      ! explicitly at all times
      TYPE(dbcsr_type), POINTER :: matrix_x => NULL(), matrix_sx => NULL(), matrix_gx => NULL()
      TYPE(dbcsr_type), POINTER :: matrix_x_im => NULL(), matrix_sx_im => NULL(), matrix_gx_im => NULL()
      TYPE(dbcsr_type), POINTER :: matrix_ref_orth => NULL(), matrix_ref_orth_im => NULL()
      TYPE(dbcsr_type), POINTER :: matrix_dx => NULL(), matrix_gx_old => NULL()
      TYPE(dbcsr_type), POINTER :: matrix_dx_im => NULL(), matrix_gx_old_im => NULL()

      LOGICAL :: use_gx_old = .FALSE., use_dx = .FALSE.

      TYPE(dbcsr_p_type), DIMENSION(:), POINTER :: matrix_h_e => NULL(), matrix_h_x => NULL()
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER :: matrix_h_e_im => NULL(), matrix_h_x_im => NULL()

      REAL(KIND=dp), DIMENSION(:, :), POINTER  :: ls_diis => NULL()
      REAL(KIND=dp), DIMENSION(:, :), POINTER  :: lss_diis => NULL()
      REAL(KIND=dp), DIMENSION(:), POINTER  :: c_diis => NULL()
      REAL(KIND=dp), DIMENSION(:), POINTER  :: c_broy => NULL()
      REAL(KIND=dp), DIMENSION(:), POINTER  :: energy_h => NULL()
      INTEGER, DIMENSION(:), POINTER         :: ipivot => NULL()

      REAL(KIND=dp)    :: ot_pos(53) = -1.0_dp, ot_energy(53) = -1.0_dp, ot_grad(53) = -1.0_dp ! HARD LIMIT FOR THE LS
      INTEGER          :: line_search_left = -1, line_search_right = -1, line_search_mid = -1
      INTEGER          :: line_search_count = -1
      LOGICAL          :: line_search_might_be_done = .FALSE.
      LOGICAL          :: mermin_ls_sample_pending = .FALSE.
      LOGICAL          :: mermin_ls_closure_rejected = .FALSE.
      INTEGER          :: mermin_ls_best_sample_index = 0
      REAL(KIND=dp)    :: mermin_ls_reference_energy = 0.0_dp, mermin_ls_best_sample_energy = 0.0_dp
      REAL(KIND=dp)    :: mermin_ls_reference_delta = 0.0_dp
      REAL(KIND=dp)    :: mermin_ls_best_sample_drop = 0.0_dp, mermin_ls_model_drop = 0.0_dp
      REAL(KIND=dp)    :: mermin_ls_model_quality = 0.0_dp, mermin_ls_best_sample_pos = 0.0_dp
      REAL(KIND=dp)    :: mermin_ls_next_pos = 0.0_dp
      INTEGER          :: mermin_ls_accepted_quality_samples = 0
      ! Accepted Mermin endpoints in the current k-point orbital subspace.
      ! Unlike response calibration, this runtime-only counter is deliberately
      ! not restored across an ADDED_MOS AUTO rebuild or external restart.
      INTEGER          :: kpoint_subspace_accepted_samples = 0
      REAL(KIND=dp)    :: mermin_ls_accepted_quality = 1.0_dp
      REAL(KIND=dp)    :: mermin_ls_accepted_quality_ema = 1.0_dp
      REAL(KIND=dp)    :: mermin_ls_accepted_gap = 0.0_dp
      REAL(KIND=dp)    :: mermin_ls_accepted_gap_ema = 0.0_dp
      ! Runtime directional derivatives at the start of a Mermin line search.
      ! They intentionally are not restart payload: a restarted in-flight
      ! search must collect a fresh direction before auditing its curvature.
      LOGICAL          :: mermin_ls_component_slope_valid = .FALSE.
      LOGICAL          :: mermin_ls_vector_start_valid = .FALSE.
      LOGICAL          :: mermin_secant_ready = .FALSE.
      ! One additional accepted secant is runtime-only and is relearned after restart/resize.
      LOGICAL          :: mermin_secant_previous_ready = .FALSE.
      LOGICAL          :: mermin_frozen_h_valid = .FALSE.
      LOGICAL          :: mermin_frozen_endpoint_valid = .FALSE.
      LOGICAL          :: mermin_closure_secant_ready = .FALSE.
      LOGICAL          :: mermin_closure_secant_previous_ready = .FALSE.
      ! One older accepted closure secant is runtime-only and provides an
      ! independent density-response direction when the latest pair is collinear.
      LOGICAL          :: mermin_closure_secant_anchor_ready = .FALSE.
      ! A causally cross-validated anchor is retained until structural history invalidation.
      LOGICAL          :: mermin_closure_secant_anchor_validated = .FALSE.
      ! Candidate identities follow the runtime-only current/previous/anchor
      ! closure secants so density-specific response samples survive intervening CG steps.
      INTEGER          :: mermin_closure_secant_flavor = 0
      INTEGER          :: mermin_closure_secant_previous_flavor = 0
      INTEGER          :: mermin_closure_secant_anchor_flavor = 0
      ! Endpoint-local reduced frozen-H response.  This runtime history is
      ! relearned after restart or a change of the channel MO dimension.
      LOGICAL          :: mermin_frozen_reduced_target_ready = .FALSE.
      REAL(KIND=dp)    :: mermin_ls_matrix_slope = 0.0_dp
      REAL(KIND=dp)    :: mermin_ls_rotation_slope = 0.0_dp
      REAL(KIND=dp)    :: mermin_ls_energy_slope = 0.0_dp
      REAL(KIND=dp)    :: mermin_secant_sty = 0.0_dp
      ! The accepted finite-energy curvature divided by the endpoint-gradient
      ! curvature follows the runtime-only current/previous/anchor secants.
      REAL(KIND=dp)    :: mermin_secant_work_ratio = 1.0_dp
      REAL(KIND=dp)    :: mermin_secant_work_ratio_previous = 1.0_dp
      REAL(KIND=dp)    :: mermin_secant_work_ratio_anchor = 1.0_dp
      ! Curvature of the complete accepted matrix/rotation/energy direction.
      ! This is runtime-only because it is meaningful only with the current CG history.
      LOGICAL          :: mermin_direction_curvature_valid = .FALSE.
      REAL(KIND=dp)    :: mermin_direction_curvature = 0.0_dp
      REAL(KIND=dp)    :: delta = -1.0_dp, gnorm = -1.0_dp, gnorm_old = -1.0_dp, etotal = -1.0_dp, gradient = -1.0_dp
      LOGICAL          :: energy_only = .FALSE.
      INTEGER          :: diis_iter = -1
      ! A rejected Pulay extrapolation invalidates its history, but must not
      ! shrink the underlying SD trust step.  Only a rejected candidate built
      ! from the base preconditioned direction contracts that step.
      LOGICAL          :: diis_last_candidate_sd = .TRUE.
      CHARACTER(LEN=8) :: OT_METHOD_FULL = ""
      INTEGER          :: OT_count = -1
      REAL(KIND=dp)    :: ds_min = -1.0_dp
      REAL(KIND=dp)    :: broyden_adaptive_sigma = -1.0_dp

      ! Complex k-point OT first approaches the local fixed-N Mermin basin
      ! with the unpreconditioned rotation/energy model.  This runtime flag
      ! distinguishes that bootstrap from the user request in settings.
      LOGICAL          :: kpoint_occ_preconditioner_active = .FALSE.
      INTEGER          :: kpoint_occ_preconditioner_ready_samples = 0
      REAL(KIND=dp)    :: kpoint_occ_metric_residual_quality_ema = 0.0_dp
      LOGICAL          :: orbital_schur_response_ready = .FALSE.

      LOGICAL          :: do_taylor = .FALSE.
      INTEGER          :: taylor_order = -1
      REAL(KIND=dp)    :: largest_eval_upper_bound = -1.0_dp

      ! second part of the variables, if an explicit rotation is required as well
      TYPE(dbcsr_type), POINTER :: rot_mat_u => NULL() ! rotation matrix
      TYPE(dbcsr_type), POINTER :: rot_mat_u_im => NULL() ! imaginary part of complex k-point rotation matrix
      TYPE(dbcsr_type), POINTER :: rot_mat_x => NULL() ! antisymmetric matrix that parametrises rot_matrix_u
      TYPE(dbcsr_type), POINTER :: rot_mat_x_im => NULL() ! symmetric imaginary part of anti-Hermitian k-point step
      TYPE(dbcsr_type), POINTER :: rot_mat_dedu => NULL() ! derivative of the total energy wrt to u
      TYPE(dbcsr_type), POINTER :: rot_mat_dedu_im => NULL()
      TYPE(dbcsr_type), POINTER :: rot_mat_chc => NULL() ! for convencience, the matrix c^T H c
      TYPE(dbcsr_type), POINTER :: rot_mat_chc_im => NULL()

      TYPE(dbcsr_p_type), DIMENSION(:), POINTER :: rot_mat_h_e => NULL()
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER :: rot_mat_h_x => NULL()
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER :: rot_mat_h_e_im => NULL()
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER :: rot_mat_h_x_im => NULL()
      TYPE(dbcsr_type), POINTER :: rot_mat_gx => NULL()
      TYPE(dbcsr_type), POINTER :: rot_mat_gx_im => NULL()
      TYPE(dbcsr_type), POINTER :: rot_mat_pgx => NULL()
      TYPE(dbcsr_type), POINTER :: rot_mat_pgx_im => NULL()
      TYPE(dbcsr_type), POINTER :: rot_mat_gx_old => NULL()
      TYPE(dbcsr_type), POINTER :: rot_mat_gx_old_im => NULL()
      TYPE(dbcsr_type), POINTER :: rot_mat_dx => NULL()
      TYPE(dbcsr_type), POINTER :: rot_mat_dx_im => NULL()
      TYPE(dbcsr_type), POINTER :: rot_mat_mermin_g0 => NULL()
      TYPE(dbcsr_type), POINTER :: rot_mat_mermin_g0_im => NULL()
      TYPE(dbcsr_type), POINTER :: rot_mat_mermin_s => NULL()
      TYPE(dbcsr_type), POINTER :: rot_mat_mermin_s_im => NULL()
      TYPE(dbcsr_type), POINTER :: rot_mat_mermin_y => NULL()
      TYPE(dbcsr_type), POINTER :: rot_mat_mermin_y_im => NULL()
      TYPE(dbcsr_type), POINTER :: rot_mat_mermin_s_previous => NULL()
      TYPE(dbcsr_type), POINTER :: rot_mat_mermin_s_previous_im => NULL()
      TYPE(dbcsr_type), POINTER :: rot_mat_mermin_y_previous => NULL()
      TYPE(dbcsr_type), POINTER :: rot_mat_mermin_y_previous_im => NULL()
      TYPE(dbcsr_type), POINTER :: rot_mat_mermin_h0_ref => NULL()
      TYPE(dbcsr_type), POINTER :: rot_mat_mermin_h0_ref_im => NULL()
      TYPE(dbcsr_type), POINTER :: rot_mat_mermin_frozen_g1 => NULL()
      TYPE(dbcsr_type), POINTER :: rot_mat_mermin_frozen_g1_im => NULL()
      TYPE(dbcsr_type), POINTER :: rot_mat_mermin_frozen_reduced_y => NULL()
      TYPE(dbcsr_type), POINTER :: rot_mat_mermin_frozen_reduced_y_im => NULL()
      TYPE(dbcsr_type), POINTER :: rot_mat_mermin_closure_y => NULL()
      TYPE(dbcsr_type), POINTER :: rot_mat_mermin_closure_y_im => NULL()
      TYPE(dbcsr_type), POINTER :: rot_mat_mermin_closure_y_previous => NULL()
      TYPE(dbcsr_type), POINTER :: rot_mat_mermin_closure_y_previous_im => NULL()
      TYPE(dbcsr_type), POINTER :: rot_mat_mermin_s_anchor => NULL()
      TYPE(dbcsr_type), POINTER :: rot_mat_mermin_s_anchor_im => NULL()
      TYPE(dbcsr_type), POINTER :: rot_mat_mermin_y_anchor => NULL()
      TYPE(dbcsr_type), POINTER :: rot_mat_mermin_y_anchor_im => NULL()
      TYPE(dbcsr_type), POINTER :: rot_mat_mermin_closure_y_anchor => NULL()
      TYPE(dbcsr_type), POINTER :: rot_mat_mermin_closure_y_anchor_im => NULL()

      REAL(KIND=dp), DIMENSION(:), POINTER :: rot_mat_evals => NULL()
      TYPE(dbcsr_type), POINTER :: rot_mat_evec_re => NULL()
      TYPE(dbcsr_type), POINTER :: rot_mat_evec_im => NULL()

      ! third part of the variables, if we need to optimize orbital energies
      REAL(KIND=dp) :: ener_x_trust_radius = -1.0_dp, ener_x_trust_mu = 0.0_dp
      REAL(KIND=dp) :: ener_x_trusted_step_max = 0.0_dp
      REAL(KIND=dp) :: ener_x_trust_scale = 1.0_dp
      LOGICAL       :: ener_x_fermi_guard_active = .FALSE.
      LOGICAL       :: ener_x_fermi_frontier_armed = .FALSE.
      REAL(KIND=dp) :: ener_x_fermi_guard_anchor = -1.0_dp
      REAL(KIND=dp) :: ener_x_fermi_guard_position_max = -1.0_dp
      REAL(KIND=dp) :: ener_x_line_search_position = 0.0_dp
      REAL(KIND=dp) :: ener_x_line_search_tangent = 1.0_dp
      REAL(KIND=dp) :: canonical_tether_strength = 0.0_dp
      REAL(KIND=dp) :: added_mos_auto_basis_signal = 0.0_dp
      ! Runtime accepted-step calibration only.  Restart files persist the OT
      ! orbitals/history, while this Schur quality model is relearned after read.
      LOGICAL :: schur_model_pending = .FALSE., schur_model_quality_ready = .FALSE.
      LOGICAL :: schur_constrained_pending = .FALSE., schur_constrained_candidate_enabled = .FALSE.
      LOGICAL :: schur_split_shadow_pending = .FALSE.
      LOGICAL :: schur_occgap_shadow_pending = .FALSE.
      LOGICAL :: schur_vector_secant_pending = .FALSE.
      LOGICAL :: schur_vector_secant_candidate_enabled = .FALSE.
      LOGICAL :: schur_multisecant_pending = .FALSE.
      LOGICAL :: schur_multisecant_candidate_enabled = .FALSE.
      LOGICAL :: schur_prediction_pending = .FALSE.
      LOGICAL :: schur_prediction_candidate_enabled = .FALSE.
      LOGICAL :: schur_reduced_secant_pending = .FALSE.
      LOGICAL :: schur_reduced_secant_shadow_validated = .FALSE.
      LOGICAL :: schur_closure_candidate_pending = .FALSE.
      LOGICAL :: schur_closure_candidate_validated = .FALSE.
      LOGICAL :: schur_projector_mode_pending = .FALSE.
      LOGICAL :: schur_projector_offdiag_pending = .FALSE.
      LOGICAL :: schur_projector_inactive_offdiag_pending = .FALSE.
      LOGICAL :: schur_projector_response_pending = .FALSE.
      LOGICAL :: schur_projector_response_candidate_enabled = .FALSE.
      LOGICAL :: schur_projector_plateau_gap_pending = .FALSE.
      LOGICAL :: schur_projector_physical_pending = .FALSE.
      LOGICAL :: schur_projector_mixed_direct_pending = .FALSE.
      LOGICAL :: schur_projector_spectral_pending = .FALSE.
      LOGICAL :: schur_projector_physical_candidate_enabled = .FALSE.
      LOGICAL :: schur_projector_spectral_candidate_enabled = .FALSE.
      LOGICAL :: schur_projector_diagref_pending = .FALSE.
      LOGICAL :: schur_projector_density_pending = .FALSE.
      LOGICAL :: schur_projector_density_audit_pending = .FALSE.
      LOGICAL :: schur_projector_density_descent_pending = .FALSE.
      LOGICAL :: schur_projector_density_descent_support_ready = .FALSE.
      LOGICAL :: schur_projector_density_full_pending = .FALSE.
      LOGICAL :: schur_projector_density_full_candidate_enabled = .FALSE.
      LOGICAL :: schur_projector_density_full_direct_pending = .FALSE.
      LOGICAL :: schur_projector_density_full_direct_candidate_enabled = .FALSE.
      LOGICAL :: schur_projector_density_full_global_pending = .FALSE.
      LOGICAL :: schur_dense_metric_shadow_pending = .FALSE.
      LOGICAL :: schur_history_mixed = .FALSE., schur_split_candidate_enabled = .FALSE.
      LOGICAL :: schur_density_probe_tail_armed = .FALSE.
      INTEGER :: schur_model_attempt = 0, schur_cg_attempt = 0, schur_model_flavor = 0, &
                 schur_last_accepted_flavor = 0, schur_cg_flavor = 0, &
                 schur_model_quality_samples = 0
      INTEGER :: schur_constrained_quality_samples = 0
      INTEGER :: schur_split_quality_samples = 0, schur_coupled_quality_samples = 0
      INTEGER :: schur_split_shadow_quality_samples = 0
      INTEGER :: schur_vector_secant_quality_samples = 0
      INTEGER :: schur_multisecant_quality_samples = 0
      INTEGER :: schur_prediction_quality_samples = 0
      INTEGER :: schur_reduced_secant_quality_samples = 0
      INTEGER :: schur_reduced_secant_model_groups = 0
      INTEGER :: schur_reduced_secant_selected_groups = 0
      INTEGER :: schur_closure_candidate_quality_samples = 0
      INTEGER :: schur_closure_candidate_model_groups = 0
      INTEGER :: schur_closure_candidate_selected_groups = 0
      INTEGER :: schur_projector_mode_quality_samples = 0
      INTEGER :: schur_projector_offdiag_quality_samples = 0
      INTEGER :: schur_projector_offdiag_selected_quality_samples = 0
      INTEGER :: schur_projector_inactive_offdiag_quality_samples = 0
      INTEGER :: schur_projector_inactive_offdiag_selected_quality_samples = 0
      INTEGER :: schur_projector_response_quality_samples = 0
      INTEGER :: schur_projector_response_selected_quality_samples = 0
      INTEGER :: schur_projector_plateau_gap_quality_samples = 0
      INTEGER :: schur_projector_physical_quality_samples = 0
      INTEGER :: schur_projector_physical_selected_quality_samples = 0
      INTEGER :: schur_projector_physical_direct_quality_samples = 0
      INTEGER :: schur_projector_physical_direct_selected_quality_samples = 0
      INTEGER :: schur_projector_phys_direct_density_sel_quality_samples = 0
      INTEGER :: schur_projector_phys_direct_density_fshadow_quality_samples = 0
      INTEGER :: schur_projector_phys_direct_density_cshadow_quality_samples = 0
      INTEGER :: schur_projector_mixed_direct_quality_samples = 0
      INTEGER :: schur_projector_spectral_quality_samples = 0
      INTEGER :: schur_projector_spectral_selected_quality_samples = 0
      INTEGER :: schur_projector_spectral_selected_failure_samples = 0
      INTEGER :: schur_projector_diagref_quality_samples = 0
      INTEGER :: schur_projector_density_quality_samples = 0
      INTEGER :: schur_projector_density_audit_samples = 0
      INTEGER :: schur_projector_density_descent_quality_samples = 0
      INTEGER :: schur_projector_density_descent_incremental_quality_samples = 0
      INTEGER :: schur_projector_density_descent_selected_quality_samples = 0
      INTEGER :: schur_projector_density_full_quality_samples = 0
      INTEGER :: schur_projector_density_full_selected_quality_samples = 0
      INTEGER :: schur_projector_density_full_direct_quality_samples = 0
      INTEGER :: schur_projector_density_fdir_sel_quality_samples = 0
      INTEGER :: schur_projector_density_full_global_quality_samples = 0
      INTEGER :: schur_projector_density_full_global_sel_quality_samples = 0
      INTEGER :: schur_dense_metric_shadow_quality_samples = 0
      REAL(KIND=dp) :: schur_model_drop = 0.0_dp, schur_model_calibrated_drop = 0.0_dp
      REAL(KIND=dp) :: schur_model_slope = 0.0_dp, schur_model_quadratic = 0.0_dp
      REAL(KIND=dp) :: schur_model_measured_drop = 0.0_dp
      REAL(KIND=dp) :: schur_model_quality = 0.0_dp, schur_model_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_model_calibrated_quality = 0.0_dp, schur_model_calibrated_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_model_gap_ratio = 0.0_dp
      REAL(KIND=dp) :: schur_model_delta = 0.0_dp, schur_model_control_delta = 0.0_dp
      REAL(KIND=dp) :: schur_model_density_stationarity_delta = -1.0_dp
      REAL(KIND=dp) :: schur_model_physical_delta = -1.0_dp
      REAL(KIND=dp) :: schur_model_residual_drop = 0.0_dp
      INTEGER :: schur_density_stationarity_quality_samples = 0
      REAL(KIND=dp) :: schur_density_stationarity_quality = 0.0_dp
      REAL(KIND=dp) :: schur_density_stationarity_quality_ema = 0.0_dp
      LOGICAL :: schur_closure_response_pending = .FALSE.
      INTEGER :: schur_closure_response_samples = 0
      REAL(KIND=dp) :: schur_closure_response_model_delta = 0.0_dp
      REAL(KIND=dp) :: schur_closure_response_model_density_delta = -1.0_dp
      REAL(KIND=dp) :: schur_closure_response_quality = 0.0_dp
      REAL(KIND=dp) :: schur_closure_response_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_closure_response_density_quality = 0.0_dp
      REAL(KIND=dp) :: schur_closure_response_density_quality_ema = 0.0_dp
      LOGICAL :: schur_closure_shadow_pending = .FALSE.
      INTEGER :: schur_closure_shadow_quality_samples = 0
      REAL(KIND=dp) :: schur_closure_shadow_baseline_quadratic = 0.0_dp
      REAL(KIND=dp) :: schur_closure_shadow_quadratic = 0.0_dp
      REAL(KIND=dp) :: schur_closure_shadow_quality = 0.0_dp
      REAL(KIND=dp) :: schur_closure_shadow_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_closure_shadow_advantage = 0.0_dp
      REAL(KIND=dp) :: schur_closure_shadow_advantage_ema = 0.0_dp
      REAL(KIND=dp) :: schur_closure_shadow_residual_quality = 0.0_dp
      REAL(KIND=dp) :: schur_closure_shadow_residual_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_closure_shadow_drop_sum = 0.0_dp
      REAL(KIND=dp) :: schur_closure_shadow_measured_sum = 0.0_dp
      LOGICAL :: schur_closure_full_shadow_pending = .FALSE.
      INTEGER :: schur_closure_full_shadow_quality_samples = 0
      REAL(KIND=dp) :: schur_closure_full_shadow_baseline_quadratic = 0.0_dp
      REAL(KIND=dp) :: schur_closure_full_shadow_quadratic = 0.0_dp
      REAL(KIND=dp) :: schur_closure_full_shadow_quality = 0.0_dp
      REAL(KIND=dp) :: schur_closure_full_shadow_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_closure_full_shadow_advantage = 0.0_dp
      REAL(KIND=dp) :: schur_closure_full_shadow_advantage_ema = 0.0_dp
      ! Runtime-only accepted-step audit of the finite closure-response predictor.  Promotion
      ! eligibility stays separate so rejected models can be calibrated without becoming active.
      LOGICAL :: schur_closure_vector_shadow_pending = .FALSE.
      LOGICAL :: schur_closure_vector_shadow_consistent = .FALSE.
      LOGICAL :: schur_closure_vector_shadow_work_consistent = .FALSE.
      INTEGER :: schur_closure_vector_shadow_quality_samples = 0
      REAL(KIND=dp) :: schur_closure_vector_shadow_baseline_quadratic = 0.0_dp
      REAL(KIND=dp) :: schur_closure_vector_shadow_quadratic = 0.0_dp
      REAL(KIND=dp) :: schur_closure_vector_shadow_quality = 0.0_dp
      REAL(KIND=dp) :: schur_closure_vector_shadow_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_closure_vector_shadow_advantage = 0.0_dp
      REAL(KIND=dp) :: schur_closure_vector_shadow_advantage_ema = 0.0_dp
      INTEGER :: schur_closure_work_shadow_quality_samples = 0
      REAL(KIND=dp) :: schur_closure_work_shadow_quality = 0.0_dp
      REAL(KIND=dp) :: schur_closure_work_shadow_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_closure_work_shadow_advantage = 0.0_dp
      REAL(KIND=dp) :: schur_closure_work_shadow_advantage_ema = 0.0_dp
      LOGICAL :: schur_residual_history_reset_pending = .FALSE.
      REAL(KIND=dp) :: schur_residual_history_quality = 1.0_dp
      INTEGER :: schur_step_gauge_quality_samples = 0
      REAL(KIND=dp) :: schur_step_gap_norm = 0.0_dp, schur_step_gauge_norm = 0.0_dp
      REAL(KIND=dp) :: schur_step_fermi_norm = 0.0_dp, schur_step_gauge_ratio = 0.0_dp
      REAL(KIND=dp) :: schur_step_gauge_ratio_ema = 0.0_dp
      REAL(KIND=dp) :: schur_step_gauge_residual_quality = 1.0_dp
      REAL(KIND=dp) :: schur_step_gauge_residual_quality_ema = 1.0_dp
      LOGICAL :: kpoint_physical_delta_ready = .FALSE.
      LOGICAL :: kpoint_density_stationarity_ready = .FALSE.
      INTEGER :: kpoint_physical_delta_groups = 0
      REAL(KIND=dp) :: kpoint_raw_delta = -1.0_dp, kpoint_physical_delta = -1.0_dp
      REAL(KIND=dp) :: kpoint_complete_physical_delta = -1.0_dp
      REAL(KIND=dp) :: kpoint_projected_delta = -1.0_dp
      REAL(KIND=dp) :: kpoint_density_stationarity_delta = -1.0_dp
      REAL(KIND=dp) :: kpoint_density_stationarity_diag = -1.0_dp
      REAL(KIND=dp) :: kpoint_density_stationarity_offdiag = -1.0_dp
      REAL(KIND=dp) :: kpoint_density_stationarity_max = -1.0_dp
      LOGICAL :: kpoint_energy_closure_attempted = .FALSE.
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :) :: kpoint_density_pair_residual
      LOGICAL :: kpoint_canonical_occupation_ready = .FALSE.
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: kpoint_canonical_occupation
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: kpoint_canonical_eigenvalue
      REAL(KIND=dp) :: kpoint_canonical_mu = 0.0_dp
      REAL(KIND=dp) :: kpoint_eps_scf = -1.0_dp
      REAL(KIND=dp) :: kpoint_physical_gnorm = 0.0_dp, kpoint_physical_rot_gnorm = 0.0_dp
      REAL(KIND=dp) :: kpoint_full_physical_rot_gnorm = 0.0_dp
      REAL(KIND=dp) :: kpoint_physical_energy_gnorm = 0.0_dp, kpoint_gauge_gnorm = 0.0_dp
      REAL(KIND=dp) :: schur_constrained_drop = 0.0_dp, schur_constrained_calibrated_drop = 0.0_dp
      REAL(KIND=dp) :: schur_constrained_quality = 0.0_dp, schur_constrained_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_constrained_calibrated_quality = 0.0_dp
      REAL(KIND=dp) :: schur_constrained_calibrated_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_constrained_drop_sum = 0.0_dp, schur_constrained_measured_sum = 0.0_dp
      REAL(KIND=dp) :: schur_constrained_calibrated_drop_sum = 0.0_dp
      REAL(KIND=dp) :: schur_constrained_gap_ratio = 0.0_dp
      REAL(KIND=dp) :: schur_split_quality = 0.0_dp, schur_split_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_split_calibrated_quality = 0.0_dp
      REAL(KIND=dp) :: schur_split_calibrated_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_split_residual_quality = 0.0_dp
      REAL(KIND=dp) :: schur_split_residual_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_split_drop_sum = 0.0_dp, schur_split_measured_sum = 0.0_dp
      REAL(KIND=dp) :: schur_split_calibrated_drop_sum = 0.0_dp
      REAL(KIND=dp) :: schur_split_shadow_drop = 0.0_dp, schur_split_shadow_calibrated_drop = 0.0_dp
      REAL(KIND=dp) :: schur_split_shadow_quality = 0.0_dp, schur_split_shadow_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_split_shadow_calibrated_quality = 0.0_dp
      REAL(KIND=dp) :: schur_split_shadow_calibrated_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_split_shadow_drop_sum = 0.0_dp, schur_split_shadow_measured_sum = 0.0_dp
      REAL(KIND=dp) :: schur_split_shadow_calibrated_drop_sum = 0.0_dp
      REAL(KIND=dp) :: schur_occgap_shadow_drop = 0.0_dp, schur_occgap_shadow_quality = 0.0_dp
      REAL(KIND=dp) :: schur_occgap_shadow_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_occgap_shadow_drop_sum = 0.0_dp
      REAL(KIND=dp) :: schur_occgap_shadow_measured_sum = 0.0_dp
      REAL(KIND=dp) :: schur_occgap_shadow_gap_ratio = 0.0_dp
      INTEGER :: schur_occgap_shadow_quality_samples = 0
      REAL(KIND=dp) :: schur_coupled_quality = 0.0_dp, schur_coupled_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_coupled_calibrated_quality = 0.0_dp
      REAL(KIND=dp) :: schur_coupled_calibrated_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_coupled_residual_quality = 0.0_dp
      REAL(KIND=dp) :: schur_coupled_residual_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_coupled_drop_sum = 0.0_dp, schur_coupled_measured_sum = 0.0_dp
      REAL(KIND=dp) :: schur_coupled_calibrated_drop_sum = 0.0_dp
      REAL(KIND=dp) :: schur_vector_secant_drop = 0.0_dp
      REAL(KIND=dp) :: schur_vector_secant_quality = 0.0_dp
      REAL(KIND=dp) :: schur_vector_secant_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_vector_secant_residual_quality = 0.0_dp
      REAL(KIND=dp) :: schur_vector_secant_residual_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_vector_secant_drop_sum = 0.0_dp
      REAL(KIND=dp) :: schur_vector_secant_measured_sum = 0.0_dp
      REAL(KIND=dp) :: schur_multisecant_drop = 0.0_dp
      REAL(KIND=dp) :: schur_multisecant_quality = 0.0_dp
      REAL(KIND=dp) :: schur_multisecant_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_multisecant_residual_quality = 0.0_dp
      REAL(KIND=dp) :: schur_multisecant_residual_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_multisecant_drop_sum = 0.0_dp
      REAL(KIND=dp) :: schur_multisecant_measured_sum = 0.0_dp
      REAL(KIND=dp) :: schur_prediction_drop = 0.0_dp
      REAL(KIND=dp) :: schur_prediction_quality = 0.0_dp
      REAL(KIND=dp) :: schur_prediction_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_prediction_residual_quality = 0.0_dp
      REAL(KIND=dp) :: schur_prediction_residual_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_prediction_drop_sum = 0.0_dp
      REAL(KIND=dp) :: schur_prediction_measured_sum = 0.0_dp
      REAL(KIND=dp) :: schur_reduced_secant_drop = 0.0_dp
      REAL(KIND=dp) :: schur_reduced_secant_quality = 0.0_dp
      REAL(KIND=dp) :: schur_reduced_secant_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_reduced_secant_residual_quality = 0.0_dp
      REAL(KIND=dp) :: schur_reduced_secant_residual_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_reduced_secant_drop_sum = 0.0_dp
      REAL(KIND=dp) :: schur_reduced_secant_measured_sum = 0.0_dp
      REAL(KIND=dp) :: schur_reduced_secant_rot_residual_start_sq = 0.0_dp
      REAL(KIND=dp) :: schur_reduced_secant_rot_residual_delta_dot = 0.0_dp
      REAL(KIND=dp) :: schur_reduced_secant_rot_residual_delta_sq = 0.0_dp
      REAL(KIND=dp) :: schur_closure_candidate_drop = 0.0_dp
      REAL(KIND=dp) :: schur_closure_candidate_quality = 0.0_dp
      REAL(KIND=dp) :: schur_closure_candidate_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_closure_candidate_residual_quality = 0.0_dp
      REAL(KIND=dp) :: schur_closure_candidate_residual_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_closure_candidate_residual_gain = 0.0_dp
      REAL(KIND=dp) :: schur_closure_candidate_residual_gain_ema = 0.0_dp
      REAL(KIND=dp) :: schur_closure_candidate_residual_coverage = 0.0_dp
      REAL(KIND=dp) :: schur_closure_candidate_residual_coverage_ema = 0.0_dp
      REAL(KIND=dp) :: schur_closure_candidate_drop_sum = 0.0_dp
      REAL(KIND=dp) :: schur_closure_candidate_measured_sum = 0.0_dp
      REAL(KIND=dp) :: schur_closure_candidate_rot_residual_start_sq = 0.0_dp
      REAL(KIND=dp) :: schur_closure_candidate_rot_residual_delta_dot = 0.0_dp
      REAL(KIND=dp) :: schur_closure_candidate_rot_residual_delta_sq = 0.0_dp
      REAL(KIND=dp) :: schur_dense_metric_shadow_drop = 0.0_dp
      REAL(KIND=dp) :: schur_dense_metric_shadow_quality = 0.0_dp
      REAL(KIND=dp) :: schur_dense_metric_shadow_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_dense_metric_shadow_drop_sum = 0.0_dp
      REAL(KIND=dp) :: schur_dense_metric_shadow_measured_sum = 0.0_dp
      REAL(KIND=dp) :: schur_projector_mode_drop = 0.0_dp
      REAL(KIND=dp) :: schur_projector_mode_response_drop = 0.0_dp
      REAL(KIND=dp) :: schur_projector_mode_step_negative = 0.0_dp
      REAL(KIND=dp) :: schur_projector_mode_response_negative = 0.0_dp
      REAL(KIND=dp) :: schur_projector_mode_trust = 0.0_dp
      REAL(KIND=dp) :: schur_projector_mode_quality = 0.0_dp
      REAL(KIND=dp) :: schur_projector_mode_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_mode_drop_sum = 0.0_dp
      REAL(KIND=dp) :: schur_projector_mode_measured_sum = 0.0_dp
      REAL(KIND=dp) :: schur_projector_mode_response_drop_sum = 0.0_dp
      REAL(KIND=dp) :: schur_projector_mode_step_negative_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_mode_response_negative_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_mode_trust_ema = 1.0_dp
      REAL(KIND=dp) :: schur_projector_offdiag_drop = 0.0_dp
      REAL(KIND=dp) :: schur_projector_offdiag_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_offdiag_drop_sum = 0.0_dp
      REAL(KIND=dp) :: schur_projector_offdiag_measured_sum = 0.0_dp
      REAL(KIND=dp) :: schur_projector_offdiag_signed_ratio = 0.0_dp
      REAL(KIND=dp) :: schur_projector_offdiag_signed_ratio_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_offdiag_abs_ratio = 0.0_dp
      REAL(KIND=dp) :: schur_projector_offdiag_abs_ratio_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_offdiag_negative = 0.0_dp
      REAL(KIND=dp) :: schur_projector_offdiag_negative_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_offdiag_step_negative = 0.0_dp
      REAL(KIND=dp) :: schur_projector_offdiag_step_negative_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_offdiag_selected = 0.0_dp
      REAL(KIND=dp) :: schur_projector_offdiag_selected_quality = 0.0_dp
      REAL(KIND=dp) :: schur_projector_offdiag_selected_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_offdiag_selected_residual_quality = 0.0_dp
      REAL(KIND=dp) :: schur_projector_offdiag_selected_residual_quality_ema = 1.0_dp
      REAL(KIND=dp) :: schur_projector_offdiag_selected_drop_sum = 0.0_dp
      REAL(KIND=dp) :: schur_projector_offdiag_selected_measured_sum = 0.0_dp
      REAL(KIND=dp) :: schur_projector_inactive_offdiag_drop = 0.0_dp
      REAL(KIND=dp) :: schur_projector_inactive_offdiag_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_inactive_offdiag_drop_sum = 0.0_dp
      REAL(KIND=dp) :: schur_projector_inactive_offdiag_measured_sum = 0.0_dp
      REAL(KIND=dp) :: schur_projector_inactive_offdiag_signed_ratio = 0.0_dp
      REAL(KIND=dp) :: schur_projector_inactive_offdiag_signed_ratio_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_inactive_offdiag_abs_ratio = 0.0_dp
      REAL(KIND=dp) :: schur_projector_inactive_offdiag_abs_ratio_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_inactive_offdiag_negative = 0.0_dp
      REAL(KIND=dp) :: schur_projector_inactive_offdiag_negative_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_inactive_offdiag_step_negative = 0.0_dp
      REAL(KIND=dp) :: schur_projector_inactive_offdiag_step_negative_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_inactive_offdiag_selected = 0.0_dp
      REAL(KIND=dp) :: schur_projector_inactive_offdiag_selected_quality = 0.0_dp
      REAL(KIND=dp) :: schur_projector_inactive_offdiag_selected_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_inactive_offdiag_selected_residual_quality = 0.0_dp
      REAL(KIND=dp) :: schur_projector_inactive_offdiag_selected_residual_quality_ema = 1.0_dp
      REAL(KIND=dp) :: schur_projector_inactive_offdiag_selected_drop_sum = 0.0_dp
      REAL(KIND=dp) :: schur_projector_inactive_offdiag_selected_measured_sum = 0.0_dp
      REAL(KIND=dp) :: schur_projector_response_drop = 0.0_dp
      REAL(KIND=dp) :: schur_projector_response_quality = 0.0_dp
      REAL(KIND=dp) :: schur_projector_response_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_response_drop_sum = 0.0_dp
      REAL(KIND=dp) :: schur_projector_response_measured_sum = 0.0_dp
      REAL(KIND=dp) :: schur_projector_response_weight = 0.0_dp
      REAL(KIND=dp) :: schur_projector_response_weight_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_response_correction_drop = 0.0_dp
      REAL(KIND=dp) :: schur_projector_response_negative = 0.0_dp
      REAL(KIND=dp) :: schur_projector_response_negative_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_response_selected = 0.0_dp
      REAL(KIND=dp) :: schur_projector_response_selected_quality = 0.0_dp
      REAL(KIND=dp) :: schur_projector_response_selected_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_response_selected_residual_quality = 0.0_dp
      REAL(KIND=dp) :: schur_projector_response_selected_residual_quality_ema = 1.0_dp
      REAL(KIND=dp) :: schur_projector_response_selected_drop_sum = 0.0_dp
      REAL(KIND=dp) :: schur_projector_response_selected_measured_sum = 0.0_dp
      REAL(KIND=dp) :: schur_projector_plateau_gap_quadratic = 0.0_dp
      REAL(KIND=dp) :: schur_projector_plateau_gap_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_plateau_gap_advantage_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_plateau_gap_drop_sum = 0.0_dp
      REAL(KIND=dp) :: schur_projector_plateau_gap_measured_sum = 0.0_dp
      REAL(KIND=dp) :: schur_projector_plateau_gap_weight = 0.0_dp
      REAL(KIND=dp) :: schur_projector_plateau_gap_weight_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_plateau_gap_negative = 0.0_dp
      REAL(KIND=dp) :: schur_projector_plateau_gap_negative_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_plateau_gap_step = 0.0_dp
      REAL(KIND=dp) :: schur_projector_plateau_gap_step_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_diagref_drop = 0.0_dp
      REAL(KIND=dp) :: schur_projector_diagref_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_diagref_drop_sum = 0.0_dp
      REAL(KIND=dp) :: schur_projector_diagref_measured_sum = 0.0_dp
      REAL(KIND=dp) :: schur_projector_diagref_weight = 0.0_dp
      REAL(KIND=dp) :: schur_projector_diagref_weight_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_diagref_correction_drop = 0.0_dp
      REAL(KIND=dp) :: schur_projector_diagref_correction_ratio_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_diagref_negative = 0.0_dp
      REAL(KIND=dp) :: schur_projector_diagref_negative_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_drop = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_quality = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_drop_sum = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_measured_sum = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_lowcurv = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_lowcurv_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_lowcurv_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_audit_amplification = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_audit_amplification_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_audit_amplification_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_audit_lowcurv = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_audit_lowcurv_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_audit_lowcurv_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_descent_drop = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_descent_quality = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_descent_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_descent_drop_sum = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_descent_measured_sum = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_descent_base = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_descent_incremental = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_descent_incremental_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_descent_incremental_drop_sum = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_descent_incremental_measured_sum = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_descent_coeff = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_descent_coeff_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_descent_support = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_descent_support_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_descent_selected = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_descent_selected_quality = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_descent_selected_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_descent_selected_residual_quality = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_descent_selected_residual_quality_ema = 1.0_dp
      REAL(KIND=dp) :: schur_projector_density_descent_selected_drop_sum = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_descent_selected_measured_sum = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_full_drop = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_full_quality = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_full_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_full_drop_sum = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_full_measured_sum = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_full_trust = 1.0_dp
      REAL(KIND=dp) :: schur_projector_density_full_trust_ema = 1.0_dp
      REAL(KIND=dp) :: schur_projector_density_full_weighted_trust = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_full_weighted_trust_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_full_abs_trust = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_full_abs_trust_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_full_selected = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_full_selected_quality = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_full_selected_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_full_selected_residual_quality = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_full_selected_residual_quality_ema = 1.0_dp
      INTEGER :: schur_projector_density_full_selected_stationarity_samples = 0
      REAL(KIND=dp) :: schur_projector_density_full_selected_stationarity_quality = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_full_selected_stationarity_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_full_selected_drop_sum = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_full_selected_measured_sum = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_full_direct_drop = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_full_direct_quality = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_full_direct_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_full_direct_drop_sum = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_full_direct_measured_sum = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_full_direct_trust = 1.0_dp
      REAL(KIND=dp) :: schur_projector_density_full_direct_trust_ema = 1.0_dp
      REAL(KIND=dp) :: schur_projector_density_full_direct_weighted_trust = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_full_direct_weighted_trust_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_full_direct_abs_trust = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_full_direct_abs_trust_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_fdir_selected = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_fdir_sel_quality = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_fdir_sel_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_fdir_sel_resid_quality = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_fdir_sel_resid_quality_ema = 1.0_dp
      REAL(KIND=dp) :: schur_projector_density_fdir_sel_drop_sum = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_fdir_sel_measured_sum = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_full_global_drop = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_full_global_quality = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_full_global_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_full_global_drop_sum = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_full_global_measured_sum = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_full_global_trust = 1.0_dp
      REAL(KIND=dp) :: schur_projector_density_full_global_trust_ema = 1.0_dp
      REAL(KIND=dp) :: schur_projector_density_full_global_weighted_trust = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_full_global_weighted_trust_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_full_global_abs_trust = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_full_global_abs_trust_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_full_global_selected = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_full_global_sel_quality = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_full_global_sel_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_full_global_sel_resid_quality = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_full_global_sel_resid_quality_ema = 1.0_dp
      REAL(KIND=dp) :: schur_projector_density_full_global_sel_drop_sum = 0.0_dp
      REAL(KIND=dp) :: schur_projector_density_full_global_sel_measured_sum = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_drop = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_quality = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_drop_sum = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_measured_sum = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_metric = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_metric_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_weight = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_weight_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_correction_drop = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_correction_ratio_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_compensation = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_compensation_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_negative = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_negative_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_direct = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_direct_quality = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_direct_quality_shadow_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_direct_drop_sum = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_direct_measured_sum = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_direct_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_direct_selected = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_direct_selected_quality = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_direct_selected_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_direct_selected_residual_quality = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_direct_selected_residual_quality_ema = 1.0_dp
      REAL(KIND=dp) :: schur_projector_physical_direct_selected_drop_sum = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_direct_selected_measured_sum = 0.0_dp
      REAL(KIND=dp) :: schur_projector_phys_direct_density_sel = 0.0_dp
      REAL(KIND=dp) :: schur_projector_phys_direct_density_sel_quality = 0.0_dp
      REAL(KIND=dp) :: schur_projector_phys_direct_density_sel_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_phys_direct_density_sel_resid_quality = 0.0_dp
      REAL(KIND=dp) :: schur_projector_phys_direct_density_sel_resid_quality_ema = 1.0_dp
      REAL(KIND=dp) :: schur_projector_phys_direct_density_sel_drop_sum = 0.0_dp
      REAL(KIND=dp) :: schur_projector_phys_direct_density_sel_measured_sum = 0.0_dp
      REAL(KIND=dp) :: schur_projector_phys_direct_density_fshadow = 0.0_dp
      REAL(KIND=dp) :: schur_projector_phys_direct_density_fshadow_quality = 0.0_dp
      REAL(KIND=dp) :: schur_projector_phys_direct_density_fshadow_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_phys_direct_density_fshadow_resid_quality = 0.0_dp
      REAL(KIND=dp) :: schur_projector_phys_direct_density_fshadow_resid_quality_ema = 1.0_dp
      REAL(KIND=dp) :: schur_projector_phys_direct_density_fshadow_coverage = 0.0_dp
      REAL(KIND=dp) :: schur_projector_phys_direct_density_fshadow_coverage_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_phys_direct_density_fshadow_drop_sum = 0.0_dp
      REAL(KIND=dp) :: schur_projector_phys_direct_density_fshadow_measured_sum = 0.0_dp
      REAL(KIND=dp) :: schur_projector_phys_direct_density_cshadow = 0.0_dp
      REAL(KIND=dp) :: schur_projector_phys_direct_density_cshadow_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_phys_direct_density_cshadow_resid_quality = 0.0_dp
      REAL(KIND=dp) :: schur_projector_phys_direct_density_cshadow_resid_quality_ema = 1.0_dp
      REAL(KIND=dp) :: schur_projector_phys_direct_density_cshadow_coverage = 0.0_dp
      REAL(KIND=dp) :: schur_projector_phys_direct_density_cshadow_coverage_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_phys_direct_density_cshadow_drop_sum = 0.0_dp
      REAL(KIND=dp) :: schur_projector_phys_direct_density_cshadow_measured_sum = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_selected = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_selected_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_selected_quality = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_selected_quality_line_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_selected_residual_quality = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_selected_residual_quality_ema = 1.0_dp
      REAL(KIND=dp) :: schur_projector_physical_selected_drop_sum = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_selected_measured_sum = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_shift = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_shift_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_trace_actual = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_trace_actual_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_trace_raw = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_trace_raw_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_trace_active = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_trace_active_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_trace_actual_diag = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_trace_actual_diag_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_trace_actual_free = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_trace_actual_free_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_trace_raw_diag = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_trace_raw_diag_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_trace_raw_free = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_trace_raw_free_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_trace_active_diag = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_trace_active_diag_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_trace_active_free = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_trace_active_free_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_trace_invariant = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_trace_invariant_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_trace_invariant_diag = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_trace_invariant_diag_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_trace_invariant_free = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_trace_invariant_free_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_weight_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_metric_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_correction_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_compensation_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_negative_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_direct_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_selected_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_shift_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_trace_actual_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_trace_raw_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_trace_active_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_trace_actual_diag_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_trace_actual_free_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_trace_raw_diag_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_trace_raw_free_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_trace_active_diag_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_trace_active_free_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_trace_invariant_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_trace_invariant_diag_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_physical_trace_invariant_free_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_mixed_direct_drop = 0.0_dp
      REAL(KIND=dp) :: schur_projector_mixed_direct_quality = 0.0_dp
      REAL(KIND=dp) :: schur_projector_mixed_direct_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_mixed_direct_drop_sum = 0.0_dp
      REAL(KIND=dp) :: schur_projector_mixed_direct_measured_sum = 0.0_dp
      REAL(KIND=dp) :: schur_projector_mixed_direct_weight = 0.0_dp
      REAL(KIND=dp) :: schur_projector_mixed_direct_weight_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_spectral_drop = 0.0_dp
      REAL(KIND=dp) :: schur_projector_spectral_quality = 0.0_dp
      REAL(KIND=dp) :: schur_projector_spectral_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_spectral_drop_sum = 0.0_dp
      REAL(KIND=dp) :: schur_projector_spectral_measured_sum = 0.0_dp
      REAL(KIND=dp) :: schur_projector_spectral_selected_quality = 0.0_dp
      REAL(KIND=dp) :: schur_projector_spectral_selected_quality_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_spectral_selected_residual_quality = 0.0_dp
      REAL(KIND=dp) :: schur_projector_spectral_selected_residual_quality_ema = 1.0_dp
      REAL(KIND=dp) :: schur_projector_spectral_selected_drop_sum = 0.0_dp
      REAL(KIND=dp) :: schur_projector_spectral_selected_measured_sum = 0.0_dp
      REAL(KIND=dp) :: schur_projector_spectral_dot = 0.0_dp
      REAL(KIND=dp) :: schur_projector_spectral_dot_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_spectral_norm = 0.0_dp
      REAL(KIND=dp) :: schur_projector_spectral_norm_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_spectral_floor = 0.0_dp
      REAL(KIND=dp) :: schur_projector_spectral_floor_ema = 0.0_dp
      REAL(KIND=dp) :: schur_projector_spectral_emin = 0.0_dp
      REAL(KIND=dp) :: schur_projector_spectral_emin_ema = 0.0_dp
      REAL(KIND=dp), POINTER, DIMENSION(:) :: ener_x => NULL()
      REAL(KIND=dp), POINTER, DIMENSION(:) :: ener_dx => NULL()
      REAL(KIND=dp), POINTER, DIMENSION(:) :: ener_gx => NULL()
      REAL(KIND=dp), POINTER, DIMENSION(:) :: ener_pgx => NULL()
      REAL(KIND=dp), POINTER, DIMENSION(:) :: ener_gx_old => NULL()
      REAL(KIND=dp), POINTER, DIMENSION(:) :: ener_mermin_g0 => NULL()
      REAL(KIND=dp), POINTER, DIMENSION(:) :: ener_mermin_x0 => NULL()
      REAL(KIND=dp), POINTER, DIMENSION(:) :: ener_mermin_s => NULL()
      REAL(KIND=dp), POINTER, DIMENSION(:) :: ener_mermin_y => NULL()
      REAL(KIND=dp), POINTER, DIMENSION(:) :: ener_mermin_s_previous => NULL()
      REAL(KIND=dp), POINTER, DIMENSION(:) :: ener_mermin_y_previous => NULL()
      REAL(KIND=dp), POINTER, DIMENSION(:) :: ener_mermin_frozen_diag => NULL()
      REAL(KIND=dp), POINTER, DIMENSION(:) :: ener_mermin_frozen_g1 => NULL()
      REAL(KIND=dp), POINTER, DIMENSION(:) :: ener_mermin_closure_y => NULL()
      REAL(KIND=dp), POINTER, DIMENSION(:) :: ener_mermin_closure_y_previous => NULL()
      REAL(KIND=dp), POINTER, DIMENSION(:) :: ener_mermin_s_anchor => NULL()
      REAL(KIND=dp), POINTER, DIMENSION(:) :: ener_mermin_y_anchor => NULL()
      REAL(KIND=dp), POINTER, DIMENSION(:) :: ener_mermin_closure_y_anchor => NULL()
      REAL(KIND=dp), POINTER, DIMENSION(:) :: ener_x_trust_weight => NULL()
      REAL(KIND=dp), POINTER, DIMENSION(:) :: matrix_occ_precond_scale => NULL()
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :) :: orbital_schur_response_metric
      REAL(KIND=dp), POINTER, DIMENSION(:, :) :: ener_h_e => NULL()
      REAL(KIND=dp), POINTER, DIMENSION(:, :) :: ener_h_x => NULL()
      LOGICAL :: kpoint_preconditioned_direction_ready = .FALSE.
      LOGICAL :: ener_direction_linear = .FALSE.
   END TYPE qs_ot_type

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_ot_types'
   INTEGER, PARAMETER :: qs_ot_restart_signature_nchar = 3, &
                         qs_ot_restart_signature_nint = 4, &
                         qs_ot_restart_signature_nlog = 6, &
                         qs_ot_restart_signature_nreal = 4

CONTAINS

! **************************************************************************************************
!> \brief resolve coupled ROTATION/ENERGIES defaults for metallic OT
!> \param do_rotation requested or defaulted rotation state
!> \param do_ener requested or defaulted auxiliary-energy state
!> \param rotation_explicit whether ROTATION was supplied explicitly
!> \param energies_explicit whether ENERGIES was supplied explicitly
!> \param selected_rotation_default whether ROTATION was enabled by this routine
!> \param selected_energies_default whether ENERGIES was enabled by this routine
! **************************************************************************************************
   PURE SUBROUTINE ot_metallic_variable_defaults( &
      do_rotation, do_ener, rotation_explicit, energies_explicit, &
      selected_rotation_default, selected_energies_default)
      LOGICAL, INTENT(INOUT)                             :: do_rotation, do_ener
      LOGICAL, INTENT(IN)                                :: rotation_explicit, energies_explicit
      LOGICAL, INTENT(OUT)                               :: selected_rotation_default, &
                                                            selected_energies_default

      selected_rotation_default = .FALSE.
      selected_energies_default = .FALSE.
      IF (.NOT. rotation_explicit .AND. .NOT. energies_explicit) THEN
         do_rotation = .TRUE.
         do_ener = .TRUE.
         selected_rotation_default = .TRUE.
         selected_energies_default = .TRUE.
      ELSE IF (.NOT. rotation_explicit) THEN
         do_rotation = do_ener
         selected_rotation_default = do_rotation
      ELSE IF (.NOT. energies_explicit) THEN
         do_ener = do_rotation
         selected_energies_default = do_ener
      END IF

   END SUBROUTINE ot_metallic_variable_defaults

! **************************************************************************************************
!> \brief pack settings that define the meaning of serialized OT minimizer history
!> \param settings current OT settings
!> \param char_signature character-valued signature fields
!> \param int_signature integer-valued signature fields
!> \param log_signature logical-valued signature fields
!> \param real_signature real-valued signature fields
! **************************************************************************************************
   PURE SUBROUTINE qs_ot_settings_restart_signature(settings, char_signature, int_signature, &
                                                    log_signature, real_signature)
      TYPE(qs_ot_settings_type), INTENT(IN)              :: settings
      CHARACTER(LEN=20), &
         DIMENSION(qs_ot_restart_signature_nchar), &
         INTENT(OUT)                                     :: char_signature
      INTEGER, DIMENSION(qs_ot_restart_signature_nint), &
         INTENT(OUT)                                     :: int_signature
      LOGICAL, DIMENSION(qs_ot_restart_signature_nlog), &
         INTENT(OUT)                                     :: log_signature
      REAL(KIND=dp), &
         DIMENSION(qs_ot_restart_signature_nreal), &
         INTENT(OUT)                                     :: real_signature

      char_signature = [CHARACTER(LEN=20) :: settings%ot_method, settings%ot_algorithm, &
                        settings%line_search_method]
      int_signature = [settings%preconditioner_type, settings%cholesky_type, &
                       settings%precond_solver_type, settings%diis_m]
      log_signature = [settings%do_rotation, settings%do_ener, &
                       settings%occupation_preconditioner, settings%add_nondiag_energy, &
                       settings%canonical_tether, settings%canonical_tether_adaptive]
      real_signature = [settings%energy_gap, settings%nondiag_energy_strength, &
                        settings%canonical_tether_strength, settings%canonical_tether_max_strength]

   END SUBROUTINE qs_ot_settings_restart_signature

! **************************************************************************************************
!> \brief test whether serialized minimizer history has the same OT coordinate model
!> \param settings current OT settings
!> \param restart_char serialized character-valued signature fields
!> \param restart_int serialized integer-valued signature fields
!> \param restart_log serialized logical-valued signature fields
!> \param restart_real serialized real-valued signature fields
!> \return true when the history can be continued without changing its mathematical meaning
! **************************************************************************************************
   PURE FUNCTION qs_ot_settings_restart_compatible(settings, restart_char, restart_int, &
                                                   restart_log, restart_real) RESULT(compatible)
      TYPE(qs_ot_settings_type), INTENT(IN)              :: settings
      CHARACTER(LEN=20), &
         DIMENSION(qs_ot_restart_signature_nchar), &
         INTENT(IN)                                      :: restart_char
      INTEGER, DIMENSION(qs_ot_restart_signature_nint), &
         INTENT(IN)                                      :: restart_int
      LOGICAL, DIMENSION(qs_ot_restart_signature_nlog), &
         INTENT(IN)                                      :: restart_log
      REAL(KIND=dp), &
         DIMENSION(qs_ot_restart_signature_nreal), &
         INTENT(IN)                                      :: restart_real
      LOGICAL                                            :: compatible

      CHARACTER(LEN=20), &
         DIMENSION(qs_ot_restart_signature_nchar)        :: current_char
      INTEGER                                            :: i
      INTEGER, DIMENSION(qs_ot_restart_signature_nint)   :: current_int
      LOGICAL, DIMENSION(qs_ot_restart_signature_nlog)   :: current_log
      REAL(KIND=dp)                                      :: scale
      REAL(KIND=dp), &
         DIMENSION(qs_ot_restart_signature_nreal)        :: current_real

      CALL qs_ot_settings_restart_signature(settings, current_char, current_int, &
                                            current_log, current_real)
      compatible = ALL(current_char == restart_char) .AND. &
                   ALL(current_int == restart_int) .AND. &
                   ALL(current_log .EQV. restart_log)
      IF (.NOT. compatible) RETURN

      DO i = 1, qs_ot_restart_signature_nreal
         scale = MAX(1.0_dp, ABS(current_real(i)), ABS(restart_real(i)))
         compatible = compatible .AND. &
                      ABS(current_real(i) - restart_real(i)) <= 1.0E-12_dp*scale
      END DO

   END FUNCTION qs_ot_settings_restart_compatible

! **************************************************************************************************
!> \brief Return whether the occupation preconditioner belongs to the current runtime model.
!> \param qs_ot_env ...
!> \return Requested state for real/Gamma OT, staged state for complex k-point OT.
! **************************************************************************************************
   PURE FUNCTION qs_ot_occupation_preconditioner_active(qs_ot_env) RESULT(active)
      TYPE(qs_ot_type), INTENT(IN)                       :: qs_ot_env
      LOGICAL                                            :: active

      active = qs_ot_env%settings%occupation_preconditioner
      IF (active .AND. qs_ot_env%has_complex_kpoint_state) THEN
         active = qs_ot_env%kpoint_occ_preconditioner_active
      END IF

   END FUNCTION qs_ot_occupation_preconditioner_active

! **************************************************************************************************
!> \brief Return whether an orbital preconditioner has an implementation for a k-point wavefunction.
!> \param preconditioner_type selected OT preconditioner
!> \param use_real_wfn whether the k-point channel stores a real wavefunction
!> \return true for the preconditioners implemented by the complex k-point OT driver
! **************************************************************************************************
   PURE ELEMENTAL FUNCTION qs_ot_kpoint_preconditioner_supported(preconditioner_type, use_real_wfn) &
      RESULT(supported)
      INTEGER, INTENT(IN)                                :: preconditioner_type
      LOGICAL, INTENT(IN)                                :: use_real_wfn
      LOGICAL                                            :: supported

      IF (use_real_wfn) THEN
         supported = .FALSE.
      ELSE
         SELECT CASE (preconditioner_type)
         CASE (ot_precond_none, ot_precond_full_single_inverse, ot_precond_full_all, &
               ot_precond_fermi_low_rank)
            supported = .TRUE.
         CASE DEFAULT
            supported = .FALSE.
         END SELECT
      END IF

   END FUNCTION qs_ot_kpoint_preconditioner_supported

! **************************************************************************************************
!> \brief Return whether a solver is implemented for a k-point OT preconditioner.
!> \param preconditioner_type selected OT preconditioner
!> \param solver_type selected preconditioner solver
!> \param use_real_wfn whether the k-point channel stores a real wavefunction
!> \return true when the complex builder implements the complete combination
! **************************************************************************************************
   PURE ELEMENTAL FUNCTION qs_ot_kpoint_preconditioner_solver_supported( &
      preconditioner_type, solver_type, use_real_wfn) RESULT(supported)
      INTEGER, INTENT(IN)                                :: preconditioner_type, solver_type
      LOGICAL, INTENT(IN)                                :: use_real_wfn
      LOGICAL                                            :: supported

      supported = qs_ot_kpoint_preconditioner_supported(preconditioner_type, use_real_wfn)
      IF (.NOT. supported) RETURN

      SELECT CASE (preconditioner_type)
      CASE (ot_precond_none)
         supported = .TRUE.
      CASE (ot_precond_full_single_inverse)
         supported = solver_type == ot_precond_solver_default .OR. &
                     solver_type == ot_precond_solver_inv_chol
      CASE (ot_precond_full_all, ot_precond_fermi_low_rank)
         supported = solver_type == ot_precond_solver_default
      CASE DEFAULT
         supported = .FALSE.
      END SELECT

   END FUNCTION qs_ot_kpoint_preconditioner_solver_supported

! **************************************************************************************************
!> \brief sets default values for the settings type
!> \param settings ...
!> \par History
!>      10.2004 created [Joost VandeVondele]
! **************************************************************************************************
   SUBROUTINE qs_ot_settings_init(settings)
      TYPE(qs_ot_settings_type)                          :: settings

      settings%ot_method = "CG"
      settings%ot_algorithm = "TOD"
      settings%diis_m = 7
      settings%preconditioner_name = "FULL_KINETIC"
      settings%preconditioner_type = ot_precond_full_kinetic
      settings%cholesky_type = cholesky_reduce
      settings%precond_solver_name = "CHOLESKY_INVERSE"
      settings%precond_solver_type = ot_precond_solver_inv_chol
      settings%chebyshev_degree = 8
      settings%line_search_method = "2PNT"
      settings%ds_min = 0.15_dp
      settings%safer_diis = .TRUE.
      settings%energy_gap = 0.2_dp
      settings%eps_taylor = 1.0E-16_dp
      settings%max_taylor = 4
      settings%gold_target = 0.01_dp
      settings%do_rotation = .FALSE.
      settings%do_ener = .FALSE.
      settings%irac_degree = 4
      settings%max_irac = 50
      settings%max_scf_diis = 0
      settings%eps_irac = 1.0E-10_dp
      settings%eps_irac_quick_exit = 1.0E-5_dp
      settings%eps_irac_switch = 1.0E-2
      settings%eps_irac_filter_matrix = 0.0_dp
      settings%on_the_fly_loc = .FALSE.
      settings%ortho_irac = "CHOL"
      settings%ks = .TRUE.
      settings%occupation_preconditioner = .FALSE.
      settings%add_nondiag_energy = .FALSE.
      settings%nondiag_energy_strength = 0.0_dp
      settings%canonical_tether = .TRUE.
      settings%canonical_tether_adaptive = .TRUE.
      settings%canonical_tether_strength = 1.0E-3_dp
      settings%canonical_tether_max_strength = 1.5E-3_dp

   END SUBROUTINE qs_ot_settings_init

! **************************************************************************************************
!> \brief label an OT environment by spin and optional irreducible k-point context
!> \param qs_ot_env ...
!> \param spin_index ...
!> \param kpoint_index ...
!> \param local_kpoint_index ...
!> \param kpoint_weight ...
!> \param kpoint_multiplicity ...
! **************************************************************************************************
   SUBROUTINE qs_ot_set_context(qs_ot_env, spin_index, kpoint_index, local_kpoint_index, kpoint_weight, &
                                kpoint_multiplicity)
      TYPE(qs_ot_type)                                   :: qs_ot_env
      INTEGER, INTENT(IN), OPTIONAL                      :: spin_index, kpoint_index, &
                                                            local_kpoint_index
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: kpoint_weight
      INTEGER, INTENT(IN), OPTIONAL                      :: kpoint_multiplicity

      IF (PRESENT(spin_index)) qs_ot_env%spin_index = spin_index
      IF (PRESENT(kpoint_index)) qs_ot_env%kpoint_index = kpoint_index
      IF (PRESENT(local_kpoint_index)) qs_ot_env%local_kpoint_index = local_kpoint_index
      IF (PRESENT(kpoint_weight)) qs_ot_env%kpoint_weight = kpoint_weight
      IF (PRESENT(kpoint_multiplicity)) THEN
         qs_ot_env%kpoint_multiplicity = MAX(1, kpoint_multiplicity)
      END IF
      qs_ot_env%has_kpoint_context = PRESENT(kpoint_index) .OR. PRESENT(local_kpoint_index)

   END SUBROUTINE qs_ot_set_context

! **************************************************************************************************
!> \brief number of OT optimization channels for spin and k-point resolved state
!> \param nspin ...
!> \param nkpoint ...
!> \param restricted ...
!> \return ...
! **************************************************************************************************
   INTEGER FUNCTION qs_ot_number_of_channels(nspin, nkpoint, restricted)
      INTEGER, INTENT(IN)                                :: nspin
      INTEGER, INTENT(IN), OPTIONAL                      :: nkpoint
      LOGICAL, INTENT(IN), OPTIONAL                      :: restricted

      INTEGER                                            :: nkpoint_eff, nspin_eff

      nspin_eff = nspin
      IF (PRESENT(restricted)) THEN
         IF (restricted) nspin_eff = 1
      END IF
      nkpoint_eff = 1
      IF (PRESENT(nkpoint)) nkpoint_eff = nkpoint

      CPASSERT(nspin_eff >= 1)
      CPASSERT(nkpoint_eff >= 1)
      qs_ot_number_of_channels = nspin_eff*nkpoint_eff

   END FUNCTION qs_ot_number_of_channels

! **************************************************************************************************
!> \brief flat OT channel index for a spin/k-point pair
!> \param ispin ...
!> \param ikpoint ...
!> \param nspin ...
!> \return ...
! **************************************************************************************************
   INTEGER FUNCTION qs_ot_channel_index(ispin, ikpoint, nspin)
      INTEGER, INTENT(IN)                                :: ispin, ikpoint, nspin

      CPASSERT(ispin >= 1)
      CPASSERT(ikpoint >= 1)
      CPASSERT(nspin >= 1)
      CPASSERT(ispin <= nspin)
      qs_ot_channel_index = (ikpoint - 1)*nspin + ispin

   END FUNCTION qs_ot_channel_index

! **************************************************************************************************
!> \brief validate the flat OT channel identity for spin and optional irreducible k-points
!> \param qs_ot_env ...
!> \param nspin ...
!> \param nkpoint ...
!> \param restricted ...
!> \param require_kpoint ...
!> \param kp_range ...
!> \param wkp ...
!> \param require_local_state ...
!> \param require_complex_state ...
! **************************************************************************************************
   SUBROUTINE qs_ot_check_channel_context(qs_ot_env, nspin, nkpoint, restricted, require_kpoint, &
                                          kp_range, wkp, require_local_state, require_complex_state)
      TYPE(qs_ot_type), DIMENSION(:), INTENT(IN)         :: qs_ot_env
      INTEGER, INTENT(IN)                                :: nspin
      INTEGER, INTENT(IN), OPTIONAL                      :: nkpoint
      LOGICAL, INTENT(IN), OPTIONAL                      :: restricted, require_kpoint
      INTEGER, DIMENSION(2), INTENT(IN), OPTIONAL        :: kp_range
      REAL(KIND=dp), DIMENSION(:), OPTIONAL, POINTER     :: wkp
      LOGICAL, INTENT(IN), OPTIONAL                      :: require_local_state, &
                                                            require_complex_state

      INTEGER                                            :: expected_local_kpoint, ikpoint, ispin, &
                                                            nkpoint_eff, nspin_eff, ot_channel
      LOGICAL                                            :: check_complex_state, &
                                                            check_kpoint_context, &
                                                            check_local_state, check_weights, &
                                                            context_ok, restricted_eff
      REAL(KIND=dp)                                      :: weight_tol

      nkpoint_eff = 1
      IF (PRESENT(nkpoint)) nkpoint_eff = nkpoint
      nspin_eff = nspin
      restricted_eff = .FALSE.
      IF (PRESENT(restricted)) THEN
         restricted_eff = restricted
         IF (restricted_eff) nspin_eff = 1
      END IF

      check_kpoint_context = (nkpoint_eff > 1)
      IF (PRESENT(require_kpoint)) check_kpoint_context = require_kpoint
      check_weights = .FALSE.
      IF (PRESENT(wkp)) check_weights = ASSOCIATED(wkp)
      check_local_state = .FALSE.
      IF (PRESENT(require_local_state)) check_local_state = require_local_state
      check_complex_state = .FALSE.
      IF (PRESENT(require_complex_state)) check_complex_state = require_complex_state

      context_ok = (SIZE(qs_ot_env) == qs_ot_number_of_channels(nspin, &
                                                                nkpoint=nkpoint_eff, &
                                                                restricted=restricted_eff))
      CPASSERT(context_ok)
      DO ikpoint = 1, nkpoint_eff
         expected_local_kpoint = 0
         IF (PRESENT(kp_range)) THEN
            IF (ikpoint >= kp_range(1) .AND. ikpoint <= kp_range(2)) THEN
               expected_local_kpoint = ikpoint - kp_range(1) + 1
            END IF
         END IF
         DO ispin = 1, nspin_eff
            ot_channel = qs_ot_channel_index(ispin, ikpoint, nspin_eff)
            CPASSERT(qs_ot_env(ot_channel)%spin_index == ispin)
            IF (check_kpoint_context) THEN
               CPASSERT(qs_ot_env(ot_channel)%has_kpoint_context)
               CPASSERT(qs_ot_env(ot_channel)%kpoint_index == ikpoint)
               IF (PRESENT(kp_range)) THEN
                  context_ok = (qs_ot_env(ot_channel)%local_kpoint_index == expected_local_kpoint)
                  CPASSERT(context_ok)
               END IF
               IF (check_weights) THEN
                  weight_tol = 1000.0_dp*EPSILON(1.0_dp)*MAX(1.0_dp, ABS(wkp(ikpoint)))
                  context_ok = (ABS(qs_ot_env(ot_channel)%kpoint_weight - wkp(ikpoint)) <= weight_tol)
                  CPASSERT(context_ok)
               END IF
            END IF
            IF (check_local_state) THEN
               IF (expected_local_kpoint > 0 .OR. .NOT. PRESENT(kp_range)) THEN
                  CPASSERT(qs_ot_env(ot_channel)%state_allocated)
                  IF (check_complex_state) THEN
                     CPASSERT(qs_ot_env(ot_channel)%has_complex_kpoint_state)
                  END IF
               ELSE
                  CPASSERT(.NOT. qs_ot_env(ot_channel)%state_allocated)
               END IF
            END IF
         END DO
      END DO

   END SUBROUTINE qs_ot_check_channel_context

! **************************************************************************************************
!> \brief init matrices, needs c0 and sc0 so that c0*sc0=1
!> \param qs_ot_env ...
! **************************************************************************************************
   SUBROUTINE qs_ot_init(qs_ot_env)
      TYPE(qs_ot_type)                                   :: qs_ot_env

      qs_ot_env%OT_energy(:) = 0.0_dp
      qs_ot_env%OT_pos(:) = 0.0_dp
      qs_ot_env%OT_grad(:) = 0.0_dp
      qs_ot_env%line_search_count = 0
      qs_ot_env%mermin_ls_sample_pending = .FALSE.
      qs_ot_env%mermin_ls_closure_rejected = .FALSE.
      qs_ot_env%mermin_ls_best_sample_index = 0
      qs_ot_env%mermin_ls_reference_energy = 0.0_dp
      qs_ot_env%mermin_ls_reference_delta = 0.0_dp
      qs_ot_env%mermin_ls_best_sample_energy = 0.0_dp
      qs_ot_env%mermin_ls_best_sample_drop = 0.0_dp
      qs_ot_env%mermin_ls_model_drop = 0.0_dp
      qs_ot_env%mermin_ls_model_quality = 0.0_dp
      qs_ot_env%mermin_ls_best_sample_pos = 0.0_dp
      qs_ot_env%mermin_ls_next_pos = 0.0_dp
      qs_ot_env%mermin_ls_accepted_quality_samples = 0
      qs_ot_env%kpoint_subspace_accepted_samples = 0
      qs_ot_env%mermin_ls_accepted_quality = 1.0_dp
      qs_ot_env%mermin_ls_accepted_quality_ema = 1.0_dp
      qs_ot_env%mermin_ls_accepted_gap = 0.0_dp
      qs_ot_env%mermin_ls_accepted_gap_ema = 0.0_dp
      qs_ot_env%mermin_ls_vector_start_valid = .FALSE.
      qs_ot_env%mermin_secant_ready = .FALSE.
      qs_ot_env%mermin_secant_previous_ready = .FALSE.
      qs_ot_env%mermin_frozen_h_valid = .FALSE.
      qs_ot_env%mermin_frozen_endpoint_valid = .FALSE.
      qs_ot_env%mermin_closure_secant_ready = .FALSE.
      qs_ot_env%mermin_closure_secant_previous_ready = .FALSE.
      qs_ot_env%mermin_closure_secant_anchor_ready = .FALSE.
      qs_ot_env%mermin_closure_secant_anchor_validated = .FALSE.
      qs_ot_env%mermin_closure_secant_flavor = 0
      qs_ot_env%mermin_closure_secant_previous_flavor = 0
      qs_ot_env%mermin_closure_secant_anchor_flavor = 0
      qs_ot_env%mermin_frozen_reduced_target_ready = .FALSE.
      qs_ot_env%mermin_secant_sty = 0.0_dp
      qs_ot_env%mermin_secant_work_ratio = 1.0_dp
      qs_ot_env%mermin_secant_work_ratio_previous = 1.0_dp
      qs_ot_env%mermin_secant_work_ratio_anchor = 1.0_dp
      qs_ot_env%mermin_direction_curvature_valid = .FALSE.
      qs_ot_env%mermin_direction_curvature = 0.0_dp

      qs_ot_env%energy_only = .FALSE.
      qs_ot_env%gnorm_old = 1.0_dp
      qs_ot_env%diis_iter = 0
      qs_ot_env%diis_last_candidate_sd = .TRUE.
      qs_ot_env%ds_min = qs_ot_env%settings%ds_min
      qs_ot_env%os_valid = .FALSE.
      qs_ot_env%kpoint_occ_preconditioner_active = .FALSE.
      qs_ot_env%kpoint_occ_preconditioner_ready_samples = 0
      qs_ot_env%kpoint_occ_metric_residual_quality_ema = 0.0_dp
      qs_ot_env%orbital_schur_response_ready = .FALSE.
      qs_ot_env%schur_model_pending = .FALSE.
      qs_ot_env%schur_model_quality_ready = .FALSE.
      qs_ot_env%schur_constrained_pending = .FALSE.
      qs_ot_env%schur_constrained_candidate_enabled = .FALSE.
      qs_ot_env%schur_split_shadow_pending = .FALSE.
      qs_ot_env%schur_occgap_shadow_pending = .FALSE.
      qs_ot_env%schur_vector_secant_pending = .FALSE.
      qs_ot_env%schur_vector_secant_candidate_enabled = .FALSE.
      qs_ot_env%schur_multisecant_pending = .FALSE.
      qs_ot_env%schur_multisecant_candidate_enabled = .FALSE.
      qs_ot_env%schur_prediction_pending = .FALSE.
      qs_ot_env%schur_prediction_candidate_enabled = .FALSE.
      qs_ot_env%schur_reduced_secant_pending = .FALSE.
      qs_ot_env%schur_reduced_secant_shadow_validated = .FALSE.
      qs_ot_env%schur_closure_candidate_pending = .FALSE.
      qs_ot_env%schur_closure_candidate_validated = .FALSE.
      qs_ot_env%schur_projector_mode_pending = .FALSE.
      qs_ot_env%schur_projector_offdiag_pending = .FALSE.
      qs_ot_env%schur_projector_inactive_offdiag_pending = .FALSE.
      qs_ot_env%schur_projector_response_pending = .FALSE.
      qs_ot_env%schur_projector_response_candidate_enabled = .FALSE.
      qs_ot_env%schur_projector_plateau_gap_pending = .FALSE.
      qs_ot_env%schur_projector_physical_pending = .FALSE.
      qs_ot_env%schur_projector_mixed_direct_pending = .FALSE.
      qs_ot_env%schur_projector_spectral_pending = .FALSE.
      qs_ot_env%schur_projector_physical_candidate_enabled = .FALSE.
      qs_ot_env%schur_projector_spectral_candidate_enabled = .FALSE.
      qs_ot_env%schur_projector_diagref_pending = .FALSE.
      qs_ot_env%schur_projector_density_pending = .FALSE.
      qs_ot_env%schur_projector_density_audit_pending = .FALSE.
      qs_ot_env%schur_projector_density_descent_pending = .FALSE.
      qs_ot_env%schur_projector_density_descent_support_ready = .FALSE.
      qs_ot_env%schur_projector_density_full_pending = .FALSE.
      qs_ot_env%schur_projector_density_full_candidate_enabled = .FALSE.
      qs_ot_env%schur_projector_density_full_direct_pending = .FALSE.
      qs_ot_env%schur_projector_density_full_direct_candidate_enabled = .FALSE.
      qs_ot_env%schur_projector_density_full_global_pending = .FALSE.
      qs_ot_env%schur_dense_metric_shadow_pending = .FALSE.
      qs_ot_env%schur_history_mixed = .FALSE.
      qs_ot_env%schur_split_candidate_enabled = .FALSE.
      qs_ot_env%schur_density_probe_tail_armed = .FALSE.
      qs_ot_env%schur_model_attempt = 0
      qs_ot_env%schur_cg_attempt = 0
      qs_ot_env%schur_model_flavor = 0
      qs_ot_env%schur_last_accepted_flavor = 0
      qs_ot_env%schur_cg_flavor = 0
      qs_ot_env%schur_model_quality_samples = 0
      qs_ot_env%schur_constrained_quality_samples = 0
      qs_ot_env%schur_split_quality_samples = 0
      qs_ot_env%schur_coupled_quality_samples = 0
      qs_ot_env%schur_split_shadow_quality_samples = 0
      qs_ot_env%schur_vector_secant_quality_samples = 0
      qs_ot_env%schur_multisecant_quality_samples = 0
      qs_ot_env%schur_prediction_quality_samples = 0
      qs_ot_env%schur_reduced_secant_quality_samples = 0
      qs_ot_env%schur_reduced_secant_model_groups = 0
      qs_ot_env%schur_reduced_secant_selected_groups = 0
      qs_ot_env%schur_closure_candidate_quality_samples = 0
      qs_ot_env%schur_closure_candidate_model_groups = 0
      qs_ot_env%schur_closure_candidate_selected_groups = 0
      qs_ot_env%schur_occgap_shadow_quality_samples = 0
      qs_ot_env%schur_projector_mode_quality_samples = 0
      qs_ot_env%schur_projector_offdiag_quality_samples = 0
      qs_ot_env%schur_projector_offdiag_selected_quality_samples = 0
      qs_ot_env%schur_projector_inactive_offdiag_quality_samples = 0
      qs_ot_env%schur_projector_inactive_offdiag_selected_quality_samples = 0
      qs_ot_env%schur_projector_response_quality_samples = 0
      qs_ot_env%schur_projector_response_selected_quality_samples = 0
      qs_ot_env%schur_projector_plateau_gap_quality_samples = 0
      qs_ot_env%schur_projector_physical_quality_samples = 0
      qs_ot_env%schur_projector_physical_selected_quality_samples = 0
      qs_ot_env%schur_projector_physical_direct_quality_samples = 0
      qs_ot_env%schur_projector_physical_direct_selected_quality_samples = 0
      qs_ot_env%schur_projector_phys_direct_density_sel_quality_samples = 0
      qs_ot_env%schur_projector_phys_direct_density_fshadow_quality_samples = 0
      qs_ot_env%schur_projector_phys_direct_density_cshadow_quality_samples = 0
      qs_ot_env%schur_projector_mixed_direct_quality_samples = 0
      qs_ot_env%schur_projector_spectral_quality_samples = 0
      qs_ot_env%schur_projector_spectral_selected_quality_samples = 0
      qs_ot_env%schur_projector_spectral_selected_failure_samples = 0
      qs_ot_env%schur_projector_diagref_quality_samples = 0
      qs_ot_env%schur_projector_density_quality_samples = 0
      qs_ot_env%schur_projector_density_audit_samples = 0
      qs_ot_env%schur_projector_density_descent_quality_samples = 0
      qs_ot_env%schur_projector_density_descent_incremental_quality_samples = 0
      qs_ot_env%schur_projector_density_descent_selected_quality_samples = 0
      qs_ot_env%schur_projector_density_full_quality_samples = 0
      qs_ot_env%schur_projector_density_full_selected_quality_samples = 0
      qs_ot_env%schur_projector_density_full_direct_quality_samples = 0
      qs_ot_env%schur_projector_density_fdir_sel_quality_samples = 0
      qs_ot_env%schur_projector_density_full_global_quality_samples = 0
      qs_ot_env%schur_projector_density_full_global_sel_quality_samples = 0
      qs_ot_env%schur_dense_metric_shadow_quality_samples = 0
      qs_ot_env%schur_model_drop = 0.0_dp
      qs_ot_env%schur_model_calibrated_drop = 0.0_dp
      qs_ot_env%schur_model_slope = 0.0_dp
      qs_ot_env%schur_model_quadratic = 0.0_dp
      qs_ot_env%schur_model_measured_drop = 0.0_dp
      qs_ot_env%schur_model_quality = 0.0_dp
      qs_ot_env%schur_model_quality_ema = 0.0_dp
      qs_ot_env%schur_model_calibrated_quality = 0.0_dp
      qs_ot_env%schur_model_calibrated_quality_ema = 0.0_dp
      qs_ot_env%schur_model_gap_ratio = 0.0_dp
      qs_ot_env%schur_model_delta = 0.0_dp
      qs_ot_env%schur_model_control_delta = 0.0_dp
      qs_ot_env%schur_model_density_stationarity_delta = -1.0_dp
      qs_ot_env%schur_model_physical_delta = -1.0_dp
      qs_ot_env%schur_model_residual_drop = 0.0_dp
      qs_ot_env%schur_density_stationarity_quality_samples = 0
      qs_ot_env%schur_density_stationarity_quality = 0.0_dp
      qs_ot_env%schur_density_stationarity_quality_ema = 0.0_dp
      qs_ot_env%schur_closure_response_pending = .FALSE.
      qs_ot_env%schur_closure_response_samples = 0
      qs_ot_env%schur_closure_response_model_delta = 0.0_dp
      qs_ot_env%schur_closure_response_model_density_delta = -1.0_dp
      qs_ot_env%schur_closure_response_quality = 0.0_dp
      qs_ot_env%schur_closure_response_quality_ema = 0.0_dp
      qs_ot_env%schur_closure_response_density_quality = 0.0_dp
      qs_ot_env%schur_closure_response_density_quality_ema = 0.0_dp
      qs_ot_env%schur_closure_shadow_pending = .FALSE.
      qs_ot_env%schur_closure_shadow_quality_samples = 0
      qs_ot_env%schur_closure_shadow_baseline_quadratic = 0.0_dp
      qs_ot_env%schur_closure_shadow_quadratic = 0.0_dp
      qs_ot_env%schur_closure_shadow_quality = 0.0_dp
      qs_ot_env%schur_closure_shadow_quality_ema = 0.0_dp
      qs_ot_env%schur_closure_shadow_advantage = 0.0_dp
      qs_ot_env%schur_closure_shadow_advantage_ema = 0.0_dp
      qs_ot_env%schur_closure_shadow_residual_quality = 0.0_dp
      qs_ot_env%schur_closure_shadow_residual_quality_ema = 0.0_dp
      qs_ot_env%schur_closure_shadow_drop_sum = 0.0_dp
      qs_ot_env%schur_closure_shadow_measured_sum = 0.0_dp
      qs_ot_env%schur_closure_full_shadow_pending = .FALSE.
      qs_ot_env%schur_closure_full_shadow_quality_samples = 0
      qs_ot_env%schur_closure_full_shadow_baseline_quadratic = 0.0_dp
      qs_ot_env%schur_closure_full_shadow_quadratic = 0.0_dp
      qs_ot_env%schur_closure_full_shadow_quality = 0.0_dp
      qs_ot_env%schur_closure_full_shadow_quality_ema = 0.0_dp
      qs_ot_env%schur_closure_full_shadow_advantage = 0.0_dp
      qs_ot_env%schur_closure_full_shadow_advantage_ema = 0.0_dp
      qs_ot_env%schur_closure_vector_shadow_pending = .FALSE.
      qs_ot_env%schur_closure_vector_shadow_consistent = .FALSE.
      qs_ot_env%schur_closure_vector_shadow_work_consistent = .FALSE.
      qs_ot_env%schur_closure_vector_shadow_quality_samples = 0
      qs_ot_env%schur_closure_vector_shadow_baseline_quadratic = 0.0_dp
      qs_ot_env%schur_closure_vector_shadow_quadratic = 0.0_dp
      qs_ot_env%schur_closure_vector_shadow_quality = 0.0_dp
      qs_ot_env%schur_closure_vector_shadow_quality_ema = 0.0_dp
      qs_ot_env%schur_closure_vector_shadow_advantage = 0.0_dp
      qs_ot_env%schur_closure_vector_shadow_advantage_ema = 0.0_dp
      qs_ot_env%schur_closure_work_shadow_quality_samples = 0
      qs_ot_env%schur_closure_work_shadow_quality = 0.0_dp
      qs_ot_env%schur_closure_work_shadow_quality_ema = 0.0_dp
      qs_ot_env%schur_closure_work_shadow_advantage = 0.0_dp
      qs_ot_env%schur_closure_work_shadow_advantage_ema = 0.0_dp
      qs_ot_env%schur_residual_history_reset_pending = .FALSE.
      qs_ot_env%schur_residual_history_quality = 1.0_dp
      qs_ot_env%schur_step_gauge_quality_samples = 0
      qs_ot_env%schur_step_gap_norm = 0.0_dp
      qs_ot_env%schur_step_gauge_norm = 0.0_dp
      qs_ot_env%schur_step_fermi_norm = 0.0_dp
      qs_ot_env%schur_step_gauge_ratio = 0.0_dp
      qs_ot_env%schur_step_gauge_ratio_ema = 0.0_dp
      qs_ot_env%schur_step_gauge_residual_quality = 1.0_dp
      qs_ot_env%schur_step_gauge_residual_quality_ema = 1.0_dp
      qs_ot_env%kpoint_physical_delta_ready = .FALSE.
      qs_ot_env%kpoint_density_stationarity_ready = .FALSE.
      qs_ot_env%kpoint_physical_delta_groups = 0
      qs_ot_env%kpoint_raw_delta = -1.0_dp
      qs_ot_env%kpoint_physical_delta = -1.0_dp
      qs_ot_env%kpoint_complete_physical_delta = -1.0_dp
      qs_ot_env%kpoint_projected_delta = -1.0_dp
      qs_ot_env%kpoint_density_stationarity_delta = -1.0_dp
      qs_ot_env%kpoint_density_stationarity_diag = -1.0_dp
      qs_ot_env%kpoint_density_stationarity_offdiag = -1.0_dp
      qs_ot_env%kpoint_density_stationarity_max = -1.0_dp
      qs_ot_env%kpoint_energy_closure_attempted = .FALSE.
      IF (ALLOCATED(qs_ot_env%kpoint_density_pair_residual)) THEN
         qs_ot_env%kpoint_density_pair_residual(:, :) = 0.0_dp
      END IF
      qs_ot_env%kpoint_canonical_occupation_ready = .FALSE.
      qs_ot_env%kpoint_canonical_mu = 0.0_dp
      IF (ALLOCATED(qs_ot_env%kpoint_canonical_occupation)) THEN
         qs_ot_env%kpoint_canonical_occupation(:) = 0.0_dp
      END IF
      IF (ALLOCATED(qs_ot_env%kpoint_canonical_eigenvalue)) THEN
         qs_ot_env%kpoint_canonical_eigenvalue(:) = 0.0_dp
      END IF
      qs_ot_env%kpoint_eps_scf = -1.0_dp
      qs_ot_env%kpoint_physical_gnorm = 0.0_dp
      qs_ot_env%kpoint_physical_rot_gnorm = 0.0_dp
      qs_ot_env%kpoint_full_physical_rot_gnorm = 0.0_dp
      qs_ot_env%kpoint_physical_energy_gnorm = 0.0_dp
      qs_ot_env%kpoint_gauge_gnorm = 0.0_dp
      qs_ot_env%schur_constrained_drop = 0.0_dp
      qs_ot_env%schur_constrained_calibrated_drop = 0.0_dp
      qs_ot_env%schur_constrained_quality = 0.0_dp
      qs_ot_env%schur_constrained_quality_ema = 0.0_dp
      qs_ot_env%schur_constrained_calibrated_quality = 0.0_dp
      qs_ot_env%schur_constrained_calibrated_quality_ema = 0.0_dp
      qs_ot_env%schur_constrained_drop_sum = 0.0_dp
      qs_ot_env%schur_constrained_measured_sum = 0.0_dp
      qs_ot_env%schur_constrained_calibrated_drop_sum = 0.0_dp
      qs_ot_env%schur_constrained_gap_ratio = 0.0_dp
      qs_ot_env%schur_split_quality = 0.0_dp
      qs_ot_env%schur_split_quality_ema = 0.0_dp
      qs_ot_env%schur_split_calibrated_quality = 0.0_dp
      qs_ot_env%schur_split_calibrated_quality_ema = 0.0_dp
      qs_ot_env%schur_split_residual_quality = 0.0_dp
      qs_ot_env%schur_split_residual_quality_ema = 0.0_dp
      qs_ot_env%schur_split_drop_sum = 0.0_dp
      qs_ot_env%schur_split_measured_sum = 0.0_dp
      qs_ot_env%schur_split_calibrated_drop_sum = 0.0_dp
      qs_ot_env%schur_split_shadow_drop = 0.0_dp
      qs_ot_env%schur_split_shadow_calibrated_drop = 0.0_dp
      qs_ot_env%schur_split_shadow_quality = 0.0_dp
      qs_ot_env%schur_split_shadow_quality_ema = 0.0_dp
      qs_ot_env%schur_split_shadow_calibrated_quality = 0.0_dp
      qs_ot_env%schur_split_shadow_calibrated_quality_ema = 0.0_dp
      qs_ot_env%schur_split_shadow_drop_sum = 0.0_dp
      qs_ot_env%schur_split_shadow_measured_sum = 0.0_dp
      qs_ot_env%schur_split_shadow_calibrated_drop_sum = 0.0_dp
      qs_ot_env%schur_occgap_shadow_drop = 0.0_dp
      qs_ot_env%schur_occgap_shadow_quality = 0.0_dp
      qs_ot_env%schur_occgap_shadow_quality_ema = 0.0_dp
      qs_ot_env%schur_occgap_shadow_drop_sum = 0.0_dp
      qs_ot_env%schur_occgap_shadow_measured_sum = 0.0_dp
      qs_ot_env%schur_occgap_shadow_gap_ratio = 0.0_dp
      qs_ot_env%schur_coupled_quality = 0.0_dp
      qs_ot_env%schur_coupled_quality_ema = 0.0_dp
      qs_ot_env%schur_coupled_calibrated_quality = 0.0_dp
      qs_ot_env%schur_coupled_calibrated_quality_ema = 0.0_dp
      qs_ot_env%schur_coupled_residual_quality = 0.0_dp
      qs_ot_env%schur_coupled_residual_quality_ema = 0.0_dp
      qs_ot_env%schur_coupled_drop_sum = 0.0_dp
      qs_ot_env%schur_coupled_measured_sum = 0.0_dp
      qs_ot_env%schur_coupled_calibrated_drop_sum = 0.0_dp
      qs_ot_env%schur_vector_secant_drop = 0.0_dp
      qs_ot_env%schur_vector_secant_quality = 0.0_dp
      qs_ot_env%schur_vector_secant_quality_ema = 0.0_dp
      qs_ot_env%schur_vector_secant_residual_quality = 0.0_dp
      qs_ot_env%schur_vector_secant_residual_quality_ema = 0.0_dp
      qs_ot_env%schur_vector_secant_drop_sum = 0.0_dp
      qs_ot_env%schur_vector_secant_measured_sum = 0.0_dp
      qs_ot_env%schur_multisecant_drop = 0.0_dp
      qs_ot_env%schur_multisecant_quality = 0.0_dp
      qs_ot_env%schur_multisecant_quality_ema = 0.0_dp
      qs_ot_env%schur_multisecant_residual_quality = 0.0_dp
      qs_ot_env%schur_multisecant_residual_quality_ema = 0.0_dp
      qs_ot_env%schur_multisecant_drop_sum = 0.0_dp
      qs_ot_env%schur_multisecant_measured_sum = 0.0_dp
      qs_ot_env%schur_prediction_drop = 0.0_dp
      qs_ot_env%schur_prediction_quality = 0.0_dp
      qs_ot_env%schur_prediction_quality_ema = 0.0_dp
      qs_ot_env%schur_prediction_residual_quality = 0.0_dp
      qs_ot_env%schur_prediction_residual_quality_ema = 0.0_dp
      qs_ot_env%schur_prediction_drop_sum = 0.0_dp
      qs_ot_env%schur_prediction_measured_sum = 0.0_dp
      qs_ot_env%schur_reduced_secant_drop = 0.0_dp
      qs_ot_env%schur_reduced_secant_quality = 0.0_dp
      qs_ot_env%schur_reduced_secant_quality_ema = 0.0_dp
      qs_ot_env%schur_reduced_secant_residual_quality = 0.0_dp
      qs_ot_env%schur_reduced_secant_residual_quality_ema = 0.0_dp
      qs_ot_env%schur_reduced_secant_drop_sum = 0.0_dp
      qs_ot_env%schur_reduced_secant_measured_sum = 0.0_dp
      qs_ot_env%schur_reduced_secant_rot_residual_start_sq = 0.0_dp
      qs_ot_env%schur_reduced_secant_rot_residual_delta_dot = 0.0_dp
      qs_ot_env%schur_reduced_secant_rot_residual_delta_sq = 0.0_dp
      qs_ot_env%schur_closure_candidate_drop = 0.0_dp
      qs_ot_env%schur_closure_candidate_quality = 0.0_dp
      qs_ot_env%schur_closure_candidate_quality_ema = 0.0_dp
      qs_ot_env%schur_closure_candidate_residual_quality = 0.0_dp
      qs_ot_env%schur_closure_candidate_residual_quality_ema = 0.0_dp
      qs_ot_env%schur_closure_candidate_residual_gain = 0.0_dp
      qs_ot_env%schur_closure_candidate_residual_gain_ema = 0.0_dp
      qs_ot_env%schur_closure_candidate_residual_coverage = 0.0_dp
      qs_ot_env%schur_closure_candidate_residual_coverage_ema = 0.0_dp
      qs_ot_env%schur_closure_candidate_drop_sum = 0.0_dp
      qs_ot_env%schur_closure_candidate_measured_sum = 0.0_dp
      qs_ot_env%schur_closure_candidate_rot_residual_start_sq = 0.0_dp
      qs_ot_env%schur_closure_candidate_rot_residual_delta_dot = 0.0_dp
      qs_ot_env%schur_closure_candidate_rot_residual_delta_sq = 0.0_dp
      qs_ot_env%schur_dense_metric_shadow_drop = 0.0_dp
      qs_ot_env%schur_dense_metric_shadow_quality = 0.0_dp
      qs_ot_env%schur_dense_metric_shadow_quality_ema = 0.0_dp
      qs_ot_env%schur_dense_metric_shadow_drop_sum = 0.0_dp
      qs_ot_env%schur_dense_metric_shadow_measured_sum = 0.0_dp
      qs_ot_env%schur_projector_mode_drop = 0.0_dp
      qs_ot_env%schur_projector_mode_response_drop = 0.0_dp
      qs_ot_env%schur_projector_mode_step_negative = 0.0_dp
      qs_ot_env%schur_projector_mode_response_negative = 0.0_dp
      qs_ot_env%schur_projector_mode_trust = 0.0_dp
      qs_ot_env%schur_projector_mode_quality = 0.0_dp
      qs_ot_env%schur_projector_mode_quality_ema = 0.0_dp
      qs_ot_env%schur_projector_mode_drop_sum = 0.0_dp
      qs_ot_env%schur_projector_mode_measured_sum = 0.0_dp
      qs_ot_env%schur_projector_mode_response_drop_sum = 0.0_dp
      qs_ot_env%schur_projector_mode_step_negative_ema = 0.0_dp
      qs_ot_env%schur_projector_mode_response_negative_ema = 0.0_dp
      qs_ot_env%schur_projector_mode_trust_ema = 1.0_dp
      qs_ot_env%schur_projector_offdiag_drop = 0.0_dp
      qs_ot_env%schur_projector_offdiag_quality_ema = 0.0_dp
      qs_ot_env%schur_projector_offdiag_drop_sum = 0.0_dp
      qs_ot_env%schur_projector_offdiag_measured_sum = 0.0_dp
      qs_ot_env%schur_projector_offdiag_signed_ratio = 0.0_dp
      qs_ot_env%schur_projector_offdiag_signed_ratio_ema = 0.0_dp
      qs_ot_env%schur_projector_offdiag_abs_ratio = 0.0_dp
      qs_ot_env%schur_projector_offdiag_abs_ratio_ema = 0.0_dp
      qs_ot_env%schur_projector_offdiag_negative = 0.0_dp
      qs_ot_env%schur_projector_offdiag_negative_ema = 0.0_dp
      qs_ot_env%schur_projector_offdiag_step_negative = 0.0_dp
      qs_ot_env%schur_projector_offdiag_step_negative_ema = 0.0_dp
      qs_ot_env%schur_projector_offdiag_selected = 0.0_dp
      qs_ot_env%schur_projector_offdiag_selected_quality = 0.0_dp
      qs_ot_env%schur_projector_offdiag_selected_quality_ema = 0.0_dp
      qs_ot_env%schur_projector_offdiag_selected_residual_quality = 0.0_dp
      qs_ot_env%schur_projector_offdiag_selected_residual_quality_ema = 1.0_dp
      qs_ot_env%schur_projector_offdiag_selected_drop_sum = 0.0_dp
      qs_ot_env%schur_projector_offdiag_selected_measured_sum = 0.0_dp
      qs_ot_env%schur_projector_inactive_offdiag_drop = 0.0_dp
      qs_ot_env%schur_projector_inactive_offdiag_quality_ema = 0.0_dp
      qs_ot_env%schur_projector_inactive_offdiag_drop_sum = 0.0_dp
      qs_ot_env%schur_projector_inactive_offdiag_measured_sum = 0.0_dp
      qs_ot_env%schur_projector_inactive_offdiag_signed_ratio = 0.0_dp
      qs_ot_env%schur_projector_inactive_offdiag_signed_ratio_ema = 0.0_dp
      qs_ot_env%schur_projector_inactive_offdiag_abs_ratio = 0.0_dp
      qs_ot_env%schur_projector_inactive_offdiag_abs_ratio_ema = 0.0_dp
      qs_ot_env%schur_projector_inactive_offdiag_negative = 0.0_dp
      qs_ot_env%schur_projector_inactive_offdiag_negative_ema = 0.0_dp
      qs_ot_env%schur_projector_inactive_offdiag_step_negative = 0.0_dp
      qs_ot_env%schur_projector_inactive_offdiag_step_negative_ema = 0.0_dp
      qs_ot_env%schur_projector_inactive_offdiag_selected = 0.0_dp
      qs_ot_env%schur_projector_inactive_offdiag_selected_quality = 0.0_dp
      qs_ot_env%schur_projector_inactive_offdiag_selected_quality_ema = 0.0_dp
      qs_ot_env%schur_projector_inactive_offdiag_selected_residual_quality = 0.0_dp
      qs_ot_env%schur_projector_inactive_offdiag_selected_residual_quality_ema = 1.0_dp
      qs_ot_env%schur_projector_inactive_offdiag_selected_drop_sum = 0.0_dp
      qs_ot_env%schur_projector_inactive_offdiag_selected_measured_sum = 0.0_dp
      qs_ot_env%schur_projector_response_drop = 0.0_dp
      qs_ot_env%schur_projector_response_quality = 0.0_dp
      qs_ot_env%schur_projector_response_quality_ema = 0.0_dp
      qs_ot_env%schur_projector_response_drop_sum = 0.0_dp
      qs_ot_env%schur_projector_response_measured_sum = 0.0_dp
      qs_ot_env%schur_projector_response_weight = 0.0_dp
      qs_ot_env%schur_projector_response_weight_ema = 0.0_dp
      qs_ot_env%schur_projector_response_correction_drop = 0.0_dp
      qs_ot_env%schur_projector_response_negative = 0.0_dp
      qs_ot_env%schur_projector_response_negative_ema = 0.0_dp
      qs_ot_env%schur_projector_response_selected = 0.0_dp
      qs_ot_env%schur_projector_response_selected_quality = 0.0_dp
      qs_ot_env%schur_projector_response_selected_quality_ema = 0.0_dp
      qs_ot_env%schur_projector_response_selected_residual_quality = 0.0_dp
      qs_ot_env%schur_projector_response_selected_residual_quality_ema = 1.0_dp
      qs_ot_env%schur_projector_response_selected_drop_sum = 0.0_dp
      qs_ot_env%schur_projector_response_selected_measured_sum = 0.0_dp
      qs_ot_env%schur_projector_plateau_gap_quadratic = 0.0_dp
      qs_ot_env%schur_projector_plateau_gap_quality_ema = 0.0_dp
      qs_ot_env%schur_projector_plateau_gap_advantage_ema = 0.0_dp
      qs_ot_env%schur_projector_plateau_gap_drop_sum = 0.0_dp
      qs_ot_env%schur_projector_plateau_gap_measured_sum = 0.0_dp
      qs_ot_env%schur_projector_plateau_gap_weight = 0.0_dp
      qs_ot_env%schur_projector_plateau_gap_weight_ema = 0.0_dp
      qs_ot_env%schur_projector_plateau_gap_negative = 0.0_dp
      qs_ot_env%schur_projector_plateau_gap_negative_ema = 0.0_dp
      qs_ot_env%schur_projector_plateau_gap_step = 0.0_dp
      qs_ot_env%schur_projector_plateau_gap_step_ema = 0.0_dp
      qs_ot_env%schur_projector_diagref_drop = 0.0_dp
      qs_ot_env%schur_projector_diagref_quality_ema = 0.0_dp
      qs_ot_env%schur_projector_diagref_drop_sum = 0.0_dp
      qs_ot_env%schur_projector_diagref_measured_sum = 0.0_dp
      qs_ot_env%schur_projector_diagref_weight = 0.0_dp
      qs_ot_env%schur_projector_diagref_weight_ema = 0.0_dp
      qs_ot_env%schur_projector_diagref_correction_drop = 0.0_dp
      qs_ot_env%schur_projector_diagref_correction_ratio_ema = 0.0_dp
      qs_ot_env%schur_projector_diagref_negative = 0.0_dp
      qs_ot_env%schur_projector_diagref_negative_ema = 0.0_dp
      qs_ot_env%schur_projector_density_drop = 0.0_dp
      qs_ot_env%schur_projector_density_quality = 0.0_dp
      qs_ot_env%schur_projector_density_quality_ema = 0.0_dp
      qs_ot_env%schur_projector_density_drop_sum = 0.0_dp
      qs_ot_env%schur_projector_density_measured_sum = 0.0_dp
      qs_ot_env%schur_projector_density_lowcurv = 0.0_dp
      qs_ot_env%schur_projector_density_lowcurv_ema = 0.0_dp
      qs_ot_env%schur_projector_density_lowcurv_quality_ema = 0.0_dp
      qs_ot_env%schur_projector_density_audit_amplification = 0.0_dp
      qs_ot_env%schur_projector_density_audit_amplification_ema = 0.0_dp
      qs_ot_env%schur_projector_density_audit_amplification_quality_ema = 0.0_dp
      qs_ot_env%schur_projector_density_audit_lowcurv = 0.0_dp
      qs_ot_env%schur_projector_density_audit_lowcurv_ema = 0.0_dp
      qs_ot_env%schur_projector_density_audit_lowcurv_quality_ema = 0.0_dp
      qs_ot_env%schur_projector_density_descent_drop = 0.0_dp
      qs_ot_env%schur_projector_density_descent_quality = 0.0_dp
      qs_ot_env%schur_projector_density_descent_quality_ema = 0.0_dp
      qs_ot_env%schur_projector_density_descent_drop_sum = 0.0_dp
      qs_ot_env%schur_projector_density_descent_measured_sum = 0.0_dp
      qs_ot_env%schur_projector_density_descent_base = 0.0_dp
      qs_ot_env%schur_projector_density_descent_incremental = 0.0_dp
      qs_ot_env%schur_projector_density_descent_incremental_quality_ema = 0.0_dp
      qs_ot_env%schur_projector_density_descent_incremental_drop_sum = 0.0_dp
      qs_ot_env%schur_projector_density_descent_incremental_measured_sum = 0.0_dp
      qs_ot_env%schur_projector_density_descent_coeff = 0.0_dp
      qs_ot_env%schur_projector_density_descent_coeff_ema = 0.0_dp
      qs_ot_env%schur_projector_density_descent_support = 0.0_dp
      qs_ot_env%schur_projector_density_descent_support_ema = 0.0_dp
      qs_ot_env%schur_projector_density_descent_selected = 0.0_dp
      qs_ot_env%schur_projector_density_descent_selected_quality = 0.0_dp
      qs_ot_env%schur_projector_density_descent_selected_quality_ema = 0.0_dp
      qs_ot_env%schur_projector_density_descent_selected_residual_quality = 0.0_dp
      qs_ot_env%schur_projector_density_descent_selected_residual_quality_ema = 1.0_dp
      qs_ot_env%schur_projector_density_descent_selected_drop_sum = 0.0_dp
      qs_ot_env%schur_projector_density_descent_selected_measured_sum = 0.0_dp
      qs_ot_env%schur_projector_density_full_drop = 0.0_dp
      qs_ot_env%schur_projector_density_full_quality = 0.0_dp
      qs_ot_env%schur_projector_density_full_quality_ema = 0.0_dp
      qs_ot_env%schur_projector_density_full_drop_sum = 0.0_dp
      qs_ot_env%schur_projector_density_full_measured_sum = 0.0_dp
      qs_ot_env%schur_projector_density_full_trust = 1.0_dp
      qs_ot_env%schur_projector_density_full_trust_ema = 1.0_dp
      qs_ot_env%schur_projector_density_full_weighted_trust = 0.0_dp
      qs_ot_env%schur_projector_density_full_weighted_trust_ema = 0.0_dp
      qs_ot_env%schur_projector_density_full_abs_trust = 0.0_dp
      qs_ot_env%schur_projector_density_full_abs_trust_ema = 0.0_dp
      qs_ot_env%schur_projector_density_full_selected = 0.0_dp
      qs_ot_env%schur_projector_density_full_selected_quality = 0.0_dp
      qs_ot_env%schur_projector_density_full_selected_quality_ema = 0.0_dp
      qs_ot_env%schur_projector_density_full_selected_residual_quality = 0.0_dp
      qs_ot_env%schur_projector_density_full_selected_residual_quality_ema = 1.0_dp
      qs_ot_env%schur_projector_density_full_selected_stationarity_samples = 0
      qs_ot_env%schur_projector_density_full_selected_stationarity_quality = 0.0_dp
      qs_ot_env%schur_projector_density_full_selected_stationarity_quality_ema = 0.0_dp
      qs_ot_env%schur_projector_density_full_selected_drop_sum = 0.0_dp
      qs_ot_env%schur_projector_density_full_selected_measured_sum = 0.0_dp
      qs_ot_env%schur_projector_density_full_direct_drop = 0.0_dp
      qs_ot_env%schur_projector_density_full_direct_quality = 0.0_dp
      qs_ot_env%schur_projector_density_full_direct_quality_ema = 0.0_dp
      qs_ot_env%schur_projector_density_full_direct_drop_sum = 0.0_dp
      qs_ot_env%schur_projector_density_full_direct_measured_sum = 0.0_dp
      qs_ot_env%schur_projector_density_full_direct_trust = 1.0_dp
      qs_ot_env%schur_projector_density_full_direct_trust_ema = 1.0_dp
      qs_ot_env%schur_projector_density_full_direct_weighted_trust = 0.0_dp
      qs_ot_env%schur_projector_density_full_direct_weighted_trust_ema = 0.0_dp
      qs_ot_env%schur_projector_density_full_direct_abs_trust = 0.0_dp
      qs_ot_env%schur_projector_density_full_direct_abs_trust_ema = 0.0_dp
      qs_ot_env%schur_projector_density_fdir_selected = 0.0_dp
      qs_ot_env%schur_projector_density_fdir_sel_quality = 0.0_dp
      qs_ot_env%schur_projector_density_fdir_sel_quality_ema = 0.0_dp
      qs_ot_env%schur_projector_density_fdir_sel_resid_quality = 0.0_dp
      qs_ot_env%schur_projector_density_fdir_sel_resid_quality_ema = 1.0_dp
      qs_ot_env%schur_projector_density_fdir_sel_drop_sum = 0.0_dp
      qs_ot_env%schur_projector_density_fdir_sel_measured_sum = 0.0_dp
      qs_ot_env%schur_projector_density_full_global_drop = 0.0_dp
      qs_ot_env%schur_projector_density_full_global_quality = 0.0_dp
      qs_ot_env%schur_projector_density_full_global_quality_ema = 0.0_dp
      qs_ot_env%schur_projector_density_full_global_drop_sum = 0.0_dp
      qs_ot_env%schur_projector_density_full_global_measured_sum = 0.0_dp
      qs_ot_env%schur_projector_density_full_global_trust = 1.0_dp
      qs_ot_env%schur_projector_density_full_global_trust_ema = 1.0_dp
      qs_ot_env%schur_projector_density_full_global_weighted_trust = 0.0_dp
      qs_ot_env%schur_projector_density_full_global_weighted_trust_ema = 0.0_dp
      qs_ot_env%schur_projector_density_full_global_abs_trust = 0.0_dp
      qs_ot_env%schur_projector_density_full_global_abs_trust_ema = 0.0_dp
      qs_ot_env%schur_projector_density_full_global_selected = 0.0_dp
      qs_ot_env%schur_projector_density_full_global_sel_quality = 0.0_dp
      qs_ot_env%schur_projector_density_full_global_sel_quality_ema = 0.0_dp
      qs_ot_env%schur_projector_density_full_global_sel_resid_quality = 0.0_dp
      qs_ot_env%schur_projector_density_full_global_sel_resid_quality_ema = 1.0_dp
      qs_ot_env%schur_projector_density_full_global_sel_drop_sum = 0.0_dp
      qs_ot_env%schur_projector_density_full_global_sel_measured_sum = 0.0_dp
      qs_ot_env%schur_projector_physical_drop = 0.0_dp
      qs_ot_env%schur_projector_physical_quality = 0.0_dp
      qs_ot_env%schur_projector_physical_quality_ema = 0.0_dp
      qs_ot_env%schur_projector_physical_drop_sum = 0.0_dp
      qs_ot_env%schur_projector_physical_measured_sum = 0.0_dp
      qs_ot_env%schur_projector_physical_metric = 0.0_dp
      qs_ot_env%schur_projector_physical_metric_ema = 0.0_dp
      qs_ot_env%schur_projector_physical_weight = 0.0_dp
      qs_ot_env%schur_projector_physical_weight_ema = 0.0_dp
      qs_ot_env%schur_projector_physical_correction_drop = 0.0_dp
      qs_ot_env%schur_projector_physical_correction_ratio_ema = 0.0_dp
      qs_ot_env%schur_projector_physical_compensation = 0.0_dp
      qs_ot_env%schur_projector_physical_compensation_ema = 0.0_dp
      qs_ot_env%schur_projector_physical_negative = 0.0_dp
      qs_ot_env%schur_projector_physical_negative_ema = 0.0_dp
      qs_ot_env%schur_projector_physical_direct = 0.0_dp
      qs_ot_env%schur_projector_physical_direct_quality = 0.0_dp
      qs_ot_env%schur_projector_physical_direct_quality_shadow_ema = 0.0_dp
      qs_ot_env%schur_projector_physical_direct_drop_sum = 0.0_dp
      qs_ot_env%schur_projector_physical_direct_measured_sum = 0.0_dp
      qs_ot_env%schur_projector_physical_direct_ema = 0.0_dp
      qs_ot_env%schur_projector_physical_direct_selected = 0.0_dp
      qs_ot_env%schur_projector_physical_direct_selected_quality = 0.0_dp
      qs_ot_env%schur_projector_physical_direct_selected_quality_ema = 0.0_dp
      qs_ot_env%schur_projector_physical_direct_selected_residual_quality = 0.0_dp
      qs_ot_env%schur_projector_physical_direct_selected_residual_quality_ema = 1.0_dp
      qs_ot_env%schur_projector_physical_direct_selected_drop_sum = 0.0_dp
      qs_ot_env%schur_projector_physical_direct_selected_measured_sum = 0.0_dp
      qs_ot_env%schur_projector_phys_direct_density_sel = 0.0_dp
      qs_ot_env%schur_projector_phys_direct_density_sel_quality = 0.0_dp
      qs_ot_env%schur_projector_phys_direct_density_sel_quality_ema = 0.0_dp
      qs_ot_env%schur_projector_phys_direct_density_sel_resid_quality = 0.0_dp
      qs_ot_env%schur_projector_phys_direct_density_sel_resid_quality_ema = 1.0_dp
      qs_ot_env%schur_projector_phys_direct_density_sel_drop_sum = 0.0_dp
      qs_ot_env%schur_projector_phys_direct_density_sel_measured_sum = 0.0_dp
      qs_ot_env%schur_projector_phys_direct_density_fshadow = 0.0_dp
      qs_ot_env%schur_projector_phys_direct_density_fshadow_quality = 0.0_dp
      qs_ot_env%schur_projector_phys_direct_density_fshadow_quality_ema = 0.0_dp
      qs_ot_env%schur_projector_phys_direct_density_fshadow_resid_quality = 0.0_dp
      qs_ot_env%schur_projector_phys_direct_density_fshadow_resid_quality_ema = 1.0_dp
      qs_ot_env%schur_projector_phys_direct_density_fshadow_coverage = 0.0_dp
      qs_ot_env%schur_projector_phys_direct_density_fshadow_coverage_ema = 0.0_dp
      qs_ot_env%schur_projector_phys_direct_density_fshadow_drop_sum = 0.0_dp
      qs_ot_env%schur_projector_phys_direct_density_fshadow_measured_sum = 0.0_dp
      qs_ot_env%schur_projector_phys_direct_density_cshadow = 0.0_dp
      qs_ot_env%schur_projector_phys_direct_density_cshadow_quality_ema = 0.0_dp
      qs_ot_env%schur_projector_phys_direct_density_cshadow_resid_quality = 0.0_dp
      qs_ot_env%schur_projector_phys_direct_density_cshadow_resid_quality_ema = 1.0_dp
      qs_ot_env%schur_projector_phys_direct_density_cshadow_coverage = 0.0_dp
      qs_ot_env%schur_projector_phys_direct_density_cshadow_coverage_ema = 0.0_dp
      qs_ot_env%schur_projector_phys_direct_density_cshadow_drop_sum = 0.0_dp
      qs_ot_env%schur_projector_phys_direct_density_cshadow_measured_sum = 0.0_dp
      qs_ot_env%schur_projector_physical_selected = 0.0_dp
      qs_ot_env%schur_projector_physical_selected_ema = 0.0_dp
      qs_ot_env%schur_projector_physical_selected_quality = 0.0_dp
      qs_ot_env%schur_projector_physical_selected_quality_line_ema = 0.0_dp
      qs_ot_env%schur_projector_physical_selected_residual_quality = 0.0_dp
      qs_ot_env%schur_projector_physical_selected_residual_quality_ema = 1.0_dp
      qs_ot_env%schur_projector_physical_selected_drop_sum = 0.0_dp
      qs_ot_env%schur_projector_physical_selected_measured_sum = 0.0_dp
      qs_ot_env%schur_projector_physical_shift = 0.0_dp
      qs_ot_env%schur_projector_physical_shift_ema = 0.0_dp
      qs_ot_env%schur_projector_physical_trace_actual = 0.0_dp
      qs_ot_env%schur_projector_physical_trace_actual_ema = 0.0_dp
      qs_ot_env%schur_projector_physical_trace_raw = 0.0_dp
      qs_ot_env%schur_projector_physical_trace_raw_ema = 0.0_dp
      qs_ot_env%schur_projector_physical_trace_active = 0.0_dp
      qs_ot_env%schur_projector_physical_trace_active_ema = 0.0_dp
      qs_ot_env%schur_projector_physical_trace_actual_diag = 0.0_dp
      qs_ot_env%schur_projector_physical_trace_actual_diag_ema = 0.0_dp
      qs_ot_env%schur_projector_physical_trace_actual_free = 0.0_dp
      qs_ot_env%schur_projector_physical_trace_actual_free_ema = 0.0_dp
      qs_ot_env%schur_projector_physical_trace_raw_diag = 0.0_dp
      qs_ot_env%schur_projector_physical_trace_raw_diag_ema = 0.0_dp
      qs_ot_env%schur_projector_physical_trace_raw_free = 0.0_dp
      qs_ot_env%schur_projector_physical_trace_raw_free_ema = 0.0_dp
      qs_ot_env%schur_projector_physical_trace_active_diag = 0.0_dp
      qs_ot_env%schur_projector_physical_trace_active_diag_ema = 0.0_dp
      qs_ot_env%schur_projector_physical_trace_active_free = 0.0_dp
      qs_ot_env%schur_projector_physical_trace_active_free_ema = 0.0_dp
      qs_ot_env%schur_projector_physical_trace_invariant = 0.0_dp
      qs_ot_env%schur_projector_physical_trace_invariant_ema = 0.0_dp
      qs_ot_env%schur_projector_physical_trace_invariant_diag = 0.0_dp
      qs_ot_env%schur_projector_physical_trace_invariant_diag_ema = 0.0_dp
      qs_ot_env%schur_projector_physical_trace_invariant_free = 0.0_dp
      qs_ot_env%schur_projector_physical_trace_invariant_free_ema = 0.0_dp
      qs_ot_env%schur_projector_physical_weight_quality_ema = 0.0_dp
      qs_ot_env%schur_projector_physical_metric_quality_ema = 0.0_dp
      qs_ot_env%schur_projector_physical_correction_quality_ema = 0.0_dp
      qs_ot_env%schur_projector_physical_compensation_quality_ema = 0.0_dp
      qs_ot_env%schur_projector_physical_negative_quality_ema = 0.0_dp
      qs_ot_env%schur_projector_physical_direct_quality_ema = 0.0_dp
      qs_ot_env%schur_projector_physical_selected_quality_ema = 0.0_dp
      qs_ot_env%schur_projector_physical_shift_quality_ema = 0.0_dp
      qs_ot_env%schur_projector_physical_trace_actual_quality_ema = 0.0_dp
      qs_ot_env%schur_projector_physical_trace_raw_quality_ema = 0.0_dp
      qs_ot_env%schur_projector_physical_trace_active_quality_ema = 0.0_dp
      qs_ot_env%schur_projector_physical_trace_actual_diag_quality_ema = 0.0_dp
      qs_ot_env%schur_projector_physical_trace_actual_free_quality_ema = 0.0_dp
      qs_ot_env%schur_projector_physical_trace_raw_diag_quality_ema = 0.0_dp
      qs_ot_env%schur_projector_physical_trace_raw_free_quality_ema = 0.0_dp
      qs_ot_env%schur_projector_physical_trace_active_diag_quality_ema = 0.0_dp
      qs_ot_env%schur_projector_physical_trace_active_free_quality_ema = 0.0_dp
      qs_ot_env%schur_projector_physical_trace_invariant_quality_ema = 0.0_dp
      qs_ot_env%schur_projector_physical_trace_invariant_diag_quality_ema = 0.0_dp
      qs_ot_env%schur_projector_physical_trace_invariant_free_quality_ema = 0.0_dp
      qs_ot_env%schur_projector_mixed_direct_drop = 0.0_dp
      qs_ot_env%schur_projector_mixed_direct_quality = 0.0_dp
      qs_ot_env%schur_projector_mixed_direct_quality_ema = 0.0_dp
      qs_ot_env%schur_projector_mixed_direct_drop_sum = 0.0_dp
      qs_ot_env%schur_projector_mixed_direct_measured_sum = 0.0_dp
      qs_ot_env%schur_projector_mixed_direct_weight = 0.0_dp
      qs_ot_env%schur_projector_mixed_direct_weight_ema = 0.0_dp
      qs_ot_env%schur_projector_spectral_drop = 0.0_dp
      qs_ot_env%schur_projector_spectral_quality = 0.0_dp
      qs_ot_env%schur_projector_spectral_quality_ema = 0.0_dp
      qs_ot_env%schur_projector_spectral_drop_sum = 0.0_dp
      qs_ot_env%schur_projector_spectral_measured_sum = 0.0_dp
      qs_ot_env%schur_projector_spectral_selected_quality = 0.0_dp
      qs_ot_env%schur_projector_spectral_selected_quality_ema = 0.0_dp
      qs_ot_env%schur_projector_spectral_selected_residual_quality = 0.0_dp
      qs_ot_env%schur_projector_spectral_selected_residual_quality_ema = 1.0_dp
      qs_ot_env%schur_projector_spectral_selected_drop_sum = 0.0_dp
      qs_ot_env%schur_projector_spectral_selected_measured_sum = 0.0_dp
      qs_ot_env%schur_projector_spectral_dot = 0.0_dp
      qs_ot_env%schur_projector_spectral_dot_ema = 0.0_dp
      qs_ot_env%schur_projector_spectral_norm = 0.0_dp
      qs_ot_env%schur_projector_spectral_norm_ema = 0.0_dp
      qs_ot_env%schur_projector_spectral_floor = 0.0_dp
      qs_ot_env%schur_projector_spectral_floor_ema = 0.0_dp
      qs_ot_env%schur_projector_spectral_emin = 0.0_dp
      qs_ot_env%schur_projector_spectral_emin_ema = 0.0_dp

      CALL dbcsr_set(qs_ot_env%matrix_gx, 0.0_dp)
      IF (ASSOCIATED(qs_ot_env%matrix_gx_im)) THEN
         CALL dbcsr_set(qs_ot_env%matrix_gx_im, 0.0_dp)
      END IF

      IF (qs_ot_env%use_dx) THEN
         CALL dbcsr_set(qs_ot_env%matrix_dx, 0.0_dp)
         IF (ASSOCIATED(qs_ot_env%matrix_dx_im)) THEN
            CALL dbcsr_set(qs_ot_env%matrix_dx_im, 0.0_dp)
         END IF
      END IF

      IF (qs_ot_env%use_gx_old) THEN
         CALL dbcsr_set(qs_ot_env%matrix_gx_old, 0.0_dp)
         IF (ASSOCIATED(qs_ot_env%matrix_gx_old_im)) THEN
            CALL dbcsr_set(qs_ot_env%matrix_gx_old_im, 0.0_dp)
         END IF
      END IF

      IF (qs_ot_env%settings%do_rotation) THEN
         CALL dbcsr_set(qs_ot_env%rot_mat_u, 0.0_dp)
         CALL dbcsr_add_on_diag(qs_ot_env%rot_mat_u, 1.0_dp)
         CALL dbcsr_set(qs_ot_env%rot_mat_x, 0.0_dp)
         CALL dbcsr_set(qs_ot_env%rot_mat_gx, 0.0_dp)
         IF (ASSOCIATED(qs_ot_env%rot_mat_u_im)) THEN
            CALL dbcsr_set(qs_ot_env%rot_mat_u_im, 0.0_dp)
         END IF
         IF (ASSOCIATED(qs_ot_env%rot_mat_x_im)) THEN
            CALL dbcsr_set(qs_ot_env%rot_mat_x_im, 0.0_dp)
         END IF
         IF (ASSOCIATED(qs_ot_env%rot_mat_gx_im)) THEN
            CALL dbcsr_set(qs_ot_env%rot_mat_gx_im, 0.0_dp)
         END IF
         IF (ASSOCIATED(qs_ot_env%rot_mat_pgx)) THEN
            CALL dbcsr_set(qs_ot_env%rot_mat_pgx, 0.0_dp)
         END IF
         IF (ASSOCIATED(qs_ot_env%rot_mat_pgx_im)) THEN
            CALL dbcsr_set(qs_ot_env%rot_mat_pgx_im, 0.0_dp)
         END IF
         IF (qs_ot_env%use_dx) THEN
            CALL dbcsr_set(qs_ot_env%rot_mat_dx, 0.0_dp)
            IF (ASSOCIATED(qs_ot_env%rot_mat_dx_im)) THEN
               CALL dbcsr_set(qs_ot_env%rot_mat_dx_im, 0.0_dp)
            END IF
         END IF
         IF (qs_ot_env%use_gx_old) THEN
            CALL dbcsr_set(qs_ot_env%rot_mat_gx_old, 0.0_dp)
            IF (ASSOCIATED(qs_ot_env%rot_mat_gx_old_im)) THEN
               CALL dbcsr_set(qs_ot_env%rot_mat_gx_old_im, 0.0_dp)
            END IF
         END IF
      END IF
      IF (qs_ot_env%settings%do_ener) THEN
         qs_ot_env%ener_gx(:) = 0.0_dp
         IF (ASSOCIATED(qs_ot_env%ener_pgx)) qs_ot_env%ener_pgx(:) = 0.0_dp
         qs_ot_env%ener_x_trust_radius = -1.0_dp
         qs_ot_env%ener_x_trust_mu = 0.0_dp
         qs_ot_env%ener_x_trusted_step_max = 0.0_dp
         qs_ot_env%ener_x_trust_scale = 1.0_dp
         qs_ot_env%ener_x_fermi_guard_active = .FALSE.
         qs_ot_env%ener_x_fermi_frontier_armed = .FALSE.
         qs_ot_env%ener_x_fermi_guard_anchor = -1.0_dp
         qs_ot_env%ener_x_fermi_guard_position_max = -1.0_dp
         qs_ot_env%ener_x_line_search_position = 0.0_dp
         qs_ot_env%ener_x_line_search_tangent = 1.0_dp
         qs_ot_env%canonical_tether_strength = 0.0_dp
         qs_ot_env%added_mos_auto_basis_signal = 0.0_dp
         IF (ASSOCIATED(qs_ot_env%ener_x_trust_weight)) qs_ot_env%ener_x_trust_weight(:) = 0.0_dp
         IF (ASSOCIATED(qs_ot_env%matrix_occ_precond_scale)) THEN
            qs_ot_env%matrix_occ_precond_scale(:) = 1.0_dp
         END IF
         IF (qs_ot_env%use_dx) THEN
            qs_ot_env%ener_dx(:) = 0.0_dp
         END IF
         IF (qs_ot_env%use_gx_old) THEN
            qs_ot_env%ener_gx_old(:) = 0.0_dp
         END IF
      END IF
      qs_ot_env%kpoint_preconditioned_direction_ready = .FALSE.
      qs_ot_env%ener_direction_linear = .FALSE.

   END SUBROUTINE qs_ot_init

! **************************************************************************************************
!> \brief Reset accepted-direction model calibration while retaining physical and minimizer history.
!> \param qs_ot_env ...
! **************************************************************************************************
   SUBROUTINE qs_ot_reset_accepted_model_calibration(qs_ot_env)
      TYPE(qs_ot_type)                                   :: qs_ot_env

      qs_ot_env%schur_model_pending = .FALSE.
      qs_ot_env%schur_model_quality_ready = .FALSE.
      qs_ot_env%schur_model_attempt = 0
      qs_ot_env%schur_model_flavor = 0
      qs_ot_env%schur_last_accepted_flavor = 0
      qs_ot_env%schur_model_drop = 0.0_dp
      qs_ot_env%schur_model_calibrated_drop = 0.0_dp
      qs_ot_env%schur_model_slope = 0.0_dp
      qs_ot_env%schur_model_quadratic = 0.0_dp
      qs_ot_env%schur_model_measured_drop = 0.0_dp
      qs_ot_env%schur_model_quality = 0.0_dp
      qs_ot_env%schur_model_quality_ema = 0.0_dp
      qs_ot_env%schur_model_calibrated_quality = 0.0_dp
      qs_ot_env%schur_model_calibrated_quality_ema = 0.0_dp
      qs_ot_env%schur_model_quality_samples = 0
      qs_ot_env%schur_prediction_pending = .FALSE.
      qs_ot_env%schur_prediction_drop = 0.0_dp
      qs_ot_env%schur_reduced_secant_pending = .FALSE.
      qs_ot_env%schur_reduced_secant_drop = 0.0_dp
      qs_ot_env%schur_reduced_secant_model_groups = 0
      qs_ot_env%schur_reduced_secant_selected_groups = 0
      qs_ot_env%schur_reduced_secant_rot_residual_start_sq = 0.0_dp
      qs_ot_env%schur_reduced_secant_rot_residual_delta_dot = 0.0_dp
      qs_ot_env%schur_reduced_secant_rot_residual_delta_sq = 0.0_dp
      qs_ot_env%schur_closure_candidate_pending = .FALSE.
      qs_ot_env%schur_closure_candidate_drop = 0.0_dp
      qs_ot_env%schur_closure_candidate_model_groups = 0
      qs_ot_env%schur_closure_candidate_selected_groups = 0
      qs_ot_env%schur_closure_candidate_rot_residual_start_sq = 0.0_dp
      qs_ot_env%schur_closure_candidate_rot_residual_delta_dot = 0.0_dp
      qs_ot_env%schur_closure_candidate_rot_residual_delta_sq = 0.0_dp
      qs_ot_env%schur_closure_response_pending = .FALSE.
      qs_ot_env%schur_closure_response_samples = 0
      qs_ot_env%schur_closure_response_model_delta = 0.0_dp
      qs_ot_env%schur_closure_response_model_density_delta = -1.0_dp
      qs_ot_env%schur_closure_response_quality = 0.0_dp
      qs_ot_env%schur_closure_response_quality_ema = 0.0_dp
      qs_ot_env%schur_closure_response_density_quality = 0.0_dp
      qs_ot_env%schur_closure_response_density_quality_ema = 0.0_dp
      qs_ot_env%schur_closure_shadow_pending = .FALSE.
      qs_ot_env%schur_closure_shadow_quality_samples = 0
      qs_ot_env%schur_closure_shadow_baseline_quadratic = 0.0_dp
      qs_ot_env%schur_closure_shadow_quadratic = 0.0_dp
      qs_ot_env%schur_closure_shadow_quality = 0.0_dp
      qs_ot_env%schur_closure_shadow_quality_ema = 0.0_dp
      qs_ot_env%schur_closure_shadow_advantage = 0.0_dp
      qs_ot_env%schur_closure_shadow_advantage_ema = 0.0_dp
      qs_ot_env%schur_closure_shadow_residual_quality = 0.0_dp
      qs_ot_env%schur_closure_shadow_residual_quality_ema = 0.0_dp
      qs_ot_env%schur_closure_shadow_drop_sum = 0.0_dp
      qs_ot_env%schur_closure_shadow_measured_sum = 0.0_dp
      qs_ot_env%schur_closure_full_shadow_pending = .FALSE.
      qs_ot_env%schur_closure_full_shadow_quality_samples = 0
      qs_ot_env%schur_closure_full_shadow_baseline_quadratic = 0.0_dp
      qs_ot_env%schur_closure_full_shadow_quadratic = 0.0_dp
      qs_ot_env%schur_closure_full_shadow_quality = 0.0_dp
      qs_ot_env%schur_closure_full_shadow_quality_ema = 0.0_dp
      qs_ot_env%schur_closure_full_shadow_advantage = 0.0_dp
      qs_ot_env%schur_closure_full_shadow_advantage_ema = 0.0_dp
      qs_ot_env%schur_closure_vector_shadow_pending = .FALSE.
      qs_ot_env%schur_closure_vector_shadow_consistent = .FALSE.
      qs_ot_env%schur_closure_vector_shadow_work_consistent = .FALSE.
      qs_ot_env%schur_closure_vector_shadow_quality_samples = 0
      qs_ot_env%schur_closure_vector_shadow_baseline_quadratic = 0.0_dp
      qs_ot_env%schur_closure_vector_shadow_quadratic = 0.0_dp
      qs_ot_env%schur_closure_vector_shadow_quality = 0.0_dp
      qs_ot_env%schur_closure_vector_shadow_quality_ema = 0.0_dp
      qs_ot_env%schur_closure_vector_shadow_advantage = 0.0_dp
      qs_ot_env%schur_closure_vector_shadow_advantage_ema = 0.0_dp
      qs_ot_env%schur_closure_work_shadow_quality_samples = 0
      qs_ot_env%schur_closure_work_shadow_quality = 0.0_dp
      qs_ot_env%schur_closure_work_shadow_quality_ema = 0.0_dp
      qs_ot_env%schur_closure_work_shadow_advantage = 0.0_dp
      qs_ot_env%schur_closure_work_shadow_advantage_ema = 0.0_dp

      qs_ot_env%schur_constrained_candidate_enabled = .FALSE.
      qs_ot_env%schur_constrained_quality = 0.0_dp
      qs_ot_env%schur_constrained_quality_ema = 0.0_dp
      qs_ot_env%schur_constrained_calibrated_quality = 0.0_dp
      qs_ot_env%schur_constrained_calibrated_quality_ema = 0.0_dp
      qs_ot_env%schur_constrained_drop_sum = 0.0_dp
      qs_ot_env%schur_constrained_measured_sum = 0.0_dp
      qs_ot_env%schur_constrained_calibrated_drop_sum = 0.0_dp
      qs_ot_env%schur_constrained_quality_samples = 0
      qs_ot_env%schur_split_candidate_enabled = .FALSE.
      qs_ot_env%schur_split_quality = 0.0_dp
      qs_ot_env%schur_split_quality_ema = 0.0_dp
      qs_ot_env%schur_split_calibrated_quality = 0.0_dp
      qs_ot_env%schur_split_calibrated_quality_ema = 0.0_dp
      qs_ot_env%schur_split_residual_quality = 0.0_dp
      qs_ot_env%schur_split_residual_quality_ema = 0.0_dp
      qs_ot_env%schur_split_drop_sum = 0.0_dp
      qs_ot_env%schur_split_measured_sum = 0.0_dp
      qs_ot_env%schur_split_calibrated_drop_sum = 0.0_dp
      qs_ot_env%schur_split_quality_samples = 0
      qs_ot_env%schur_coupled_quality = 0.0_dp
      qs_ot_env%schur_coupled_quality_ema = 0.0_dp
      qs_ot_env%schur_coupled_calibrated_quality = 0.0_dp
      qs_ot_env%schur_coupled_calibrated_quality_ema = 0.0_dp
      qs_ot_env%schur_coupled_residual_quality = 0.0_dp
      qs_ot_env%schur_coupled_residual_quality_ema = 0.0_dp
      qs_ot_env%schur_coupled_drop_sum = 0.0_dp
      qs_ot_env%schur_coupled_measured_sum = 0.0_dp
      qs_ot_env%schur_coupled_calibrated_drop_sum = 0.0_dp
      qs_ot_env%schur_coupled_quality_samples = 0

      qs_ot_env%schur_vector_secant_candidate_enabled = .FALSE.
      qs_ot_env%schur_vector_secant_quality_samples = 0
      qs_ot_env%schur_vector_secant_quality = 0.0_dp
      qs_ot_env%schur_vector_secant_quality_ema = 0.0_dp
      qs_ot_env%schur_vector_secant_residual_quality = 0.0_dp
      qs_ot_env%schur_vector_secant_residual_quality_ema = 0.0_dp
      qs_ot_env%schur_vector_secant_drop_sum = 0.0_dp
      qs_ot_env%schur_vector_secant_measured_sum = 0.0_dp
      qs_ot_env%schur_multisecant_candidate_enabled = .FALSE.
      qs_ot_env%schur_multisecant_quality_samples = 0
      qs_ot_env%schur_multisecant_quality = 0.0_dp
      qs_ot_env%schur_multisecant_quality_ema = 0.0_dp
      qs_ot_env%schur_multisecant_residual_quality = 0.0_dp
      qs_ot_env%schur_multisecant_residual_quality_ema = 0.0_dp
      qs_ot_env%schur_multisecant_drop_sum = 0.0_dp
      qs_ot_env%schur_multisecant_measured_sum = 0.0_dp
      qs_ot_env%schur_prediction_candidate_enabled = .FALSE.
      qs_ot_env%schur_prediction_quality_samples = 0
      qs_ot_env%schur_prediction_quality = 0.0_dp
      qs_ot_env%schur_prediction_quality_ema = 0.0_dp
      qs_ot_env%schur_prediction_residual_quality = 0.0_dp
      qs_ot_env%schur_prediction_residual_quality_ema = 0.0_dp
      qs_ot_env%schur_prediction_drop_sum = 0.0_dp
      qs_ot_env%schur_prediction_measured_sum = 0.0_dp
      qs_ot_env%schur_reduced_secant_shadow_validated = .FALSE.
      qs_ot_env%schur_reduced_secant_quality_samples = 0
      qs_ot_env%schur_reduced_secant_model_groups = 0
      qs_ot_env%schur_reduced_secant_selected_groups = 0
      qs_ot_env%schur_reduced_secant_quality = 0.0_dp
      qs_ot_env%schur_reduced_secant_quality_ema = 0.0_dp
      qs_ot_env%schur_reduced_secant_residual_quality = 0.0_dp
      qs_ot_env%schur_reduced_secant_residual_quality_ema = 0.0_dp
      qs_ot_env%schur_reduced_secant_drop_sum = 0.0_dp
      qs_ot_env%schur_reduced_secant_measured_sum = 0.0_dp
      qs_ot_env%schur_reduced_secant_rot_residual_start_sq = 0.0_dp
      qs_ot_env%schur_reduced_secant_rot_residual_delta_dot = 0.0_dp
      qs_ot_env%schur_reduced_secant_rot_residual_delta_sq = 0.0_dp
      qs_ot_env%schur_closure_candidate_validated = .FALSE.
      qs_ot_env%schur_closure_candidate_quality_samples = 0
      qs_ot_env%schur_closure_candidate_model_groups = 0
      qs_ot_env%schur_closure_candidate_selected_groups = 0
      qs_ot_env%schur_closure_candidate_quality = 0.0_dp
      qs_ot_env%schur_closure_candidate_quality_ema = 0.0_dp
      qs_ot_env%schur_closure_candidate_residual_quality = 0.0_dp
      qs_ot_env%schur_closure_candidate_residual_quality_ema = 0.0_dp
      qs_ot_env%schur_closure_candidate_residual_gain = 0.0_dp
      qs_ot_env%schur_closure_candidate_residual_gain_ema = 0.0_dp
      qs_ot_env%schur_closure_candidate_residual_coverage = 0.0_dp
      qs_ot_env%schur_closure_candidate_residual_coverage_ema = 0.0_dp
      qs_ot_env%schur_closure_candidate_drop_sum = 0.0_dp
      qs_ot_env%schur_closure_candidate_measured_sum = 0.0_dp
      qs_ot_env%schur_closure_candidate_rot_residual_start_sq = 0.0_dp
      qs_ot_env%schur_closure_candidate_rot_residual_delta_dot = 0.0_dp
      qs_ot_env%schur_closure_candidate_rot_residual_delta_sq = 0.0_dp

      CALL reset_selected_calibration( &
         qs_ot_env%schur_projector_offdiag_selected_quality_samples, &
         qs_ot_env%schur_projector_offdiag_selected_quality, &
         qs_ot_env%schur_projector_offdiag_selected_quality_ema, &
         qs_ot_env%schur_projector_offdiag_selected_residual_quality, &
         qs_ot_env%schur_projector_offdiag_selected_residual_quality_ema, &
         qs_ot_env%schur_projector_offdiag_selected_drop_sum, &
         qs_ot_env%schur_projector_offdiag_selected_measured_sum)
      CALL reset_selected_calibration( &
         qs_ot_env%schur_projector_inactive_offdiag_selected_quality_samples, &
         qs_ot_env%schur_projector_inactive_offdiag_selected_quality, &
         qs_ot_env%schur_projector_inactive_offdiag_selected_quality_ema, &
         qs_ot_env%schur_projector_inactive_offdiag_selected_residual_quality, &
         qs_ot_env%schur_projector_inactive_offdiag_selected_residual_quality_ema, &
         qs_ot_env%schur_projector_inactive_offdiag_selected_drop_sum, &
         qs_ot_env%schur_projector_inactive_offdiag_selected_measured_sum)
      CALL reset_selected_calibration( &
         qs_ot_env%schur_projector_response_selected_quality_samples, &
         qs_ot_env%schur_projector_response_selected_quality, &
         qs_ot_env%schur_projector_response_selected_quality_ema, &
         qs_ot_env%schur_projector_response_selected_residual_quality, &
         qs_ot_env%schur_projector_response_selected_residual_quality_ema, &
         qs_ot_env%schur_projector_response_selected_drop_sum, &
         qs_ot_env%schur_projector_response_selected_measured_sum)
      CALL reset_selected_calibration( &
         qs_ot_env%schur_projector_physical_selected_quality_samples, &
         qs_ot_env%schur_projector_physical_selected_quality, &
         qs_ot_env%schur_projector_physical_selected_quality_line_ema, &
         qs_ot_env%schur_projector_physical_selected_residual_quality, &
         qs_ot_env%schur_projector_physical_selected_residual_quality_ema, &
         qs_ot_env%schur_projector_physical_selected_drop_sum, &
         qs_ot_env%schur_projector_physical_selected_measured_sum)
      qs_ot_env%schur_projector_physical_selected_quality_ema = 0.0_dp
      CALL reset_selected_calibration( &
         qs_ot_env%schur_projector_physical_direct_selected_quality_samples, &
         qs_ot_env%schur_projector_physical_direct_selected_quality, &
         qs_ot_env%schur_projector_physical_direct_selected_quality_ema, &
         qs_ot_env%schur_projector_physical_direct_selected_residual_quality, &
         qs_ot_env%schur_projector_physical_direct_selected_residual_quality_ema, &
         qs_ot_env%schur_projector_physical_direct_selected_drop_sum, &
         qs_ot_env%schur_projector_physical_direct_selected_measured_sum)
      CALL reset_selected_calibration( &
         qs_ot_env%schur_projector_phys_direct_density_sel_quality_samples, &
         qs_ot_env%schur_projector_phys_direct_density_sel_quality, &
         qs_ot_env%schur_projector_phys_direct_density_sel_quality_ema, &
         qs_ot_env%schur_projector_phys_direct_density_sel_resid_quality, &
         qs_ot_env%schur_projector_phys_direct_density_sel_resid_quality_ema, &
         qs_ot_env%schur_projector_phys_direct_density_sel_drop_sum, &
         qs_ot_env%schur_projector_phys_direct_density_sel_measured_sum)
      CALL reset_selected_calibration( &
         qs_ot_env%schur_projector_density_descent_selected_quality_samples, &
         qs_ot_env%schur_projector_density_descent_selected_quality, &
         qs_ot_env%schur_projector_density_descent_selected_quality_ema, &
         qs_ot_env%schur_projector_density_descent_selected_residual_quality, &
         qs_ot_env%schur_projector_density_descent_selected_residual_quality_ema, &
         qs_ot_env%schur_projector_density_descent_selected_drop_sum, &
         qs_ot_env%schur_projector_density_descent_selected_measured_sum)
      CALL reset_selected_calibration( &
         qs_ot_env%schur_projector_density_full_selected_quality_samples, &
         qs_ot_env%schur_projector_density_full_selected_quality, &
         qs_ot_env%schur_projector_density_full_selected_quality_ema, &
         qs_ot_env%schur_projector_density_full_selected_residual_quality, &
         qs_ot_env%schur_projector_density_full_selected_residual_quality_ema, &
         qs_ot_env%schur_projector_density_full_selected_drop_sum, &
         qs_ot_env%schur_projector_density_full_selected_measured_sum)
      CALL reset_selected_calibration( &
         qs_ot_env%schur_projector_density_fdir_sel_quality_samples, &
         qs_ot_env%schur_projector_density_fdir_sel_quality, &
         qs_ot_env%schur_projector_density_fdir_sel_quality_ema, &
         qs_ot_env%schur_projector_density_fdir_sel_resid_quality, &
         qs_ot_env%schur_projector_density_fdir_sel_resid_quality_ema, &
         qs_ot_env%schur_projector_density_fdir_sel_drop_sum, &
         qs_ot_env%schur_projector_density_fdir_sel_measured_sum)
      CALL reset_selected_calibration( &
         qs_ot_env%schur_projector_density_full_global_sel_quality_samples, &
         qs_ot_env%schur_projector_density_full_global_sel_quality, &
         qs_ot_env%schur_projector_density_full_global_sel_quality_ema, &
         qs_ot_env%schur_projector_density_full_global_sel_resid_quality, &
         qs_ot_env%schur_projector_density_full_global_sel_resid_quality_ema, &
         qs_ot_env%schur_projector_density_full_global_sel_drop_sum, &
         qs_ot_env%schur_projector_density_full_global_sel_measured_sum)
      qs_ot_env%schur_projector_spectral_selected_quality_samples = 0
      qs_ot_env%schur_projector_spectral_selected_failure_samples = 0
      qs_ot_env%schur_projector_spectral_selected_quality = 0.0_dp
      qs_ot_env%schur_projector_spectral_selected_quality_ema = 0.0_dp
      qs_ot_env%schur_projector_spectral_selected_residual_quality = 0.0_dp
      qs_ot_env%schur_projector_spectral_selected_residual_quality_ema = 1.0_dp
      qs_ot_env%schur_projector_spectral_selected_drop_sum = 0.0_dp
      qs_ot_env%schur_projector_spectral_selected_measured_sum = 0.0_dp

   CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param samples ...
!> \param quality ...
!> \param quality_ema ...
!> \param residual_quality ...
!> \param residual_quality_ema ...
!> \param drop_sum ...
!> \param measured_sum ...
! **************************************************************************************************
      SUBROUTINE reset_selected_calibration(samples, quality, quality_ema, residual_quality, &
                                            residual_quality_ema, drop_sum, measured_sum)
      INTEGER, INTENT(OUT)                               :: samples
      REAL(KIND=dp), INTENT(OUT)                         :: quality, quality_ema, residual_quality, &
                                                            residual_quality_ema, drop_sum, &
                                                            measured_sum

         samples = 0
         quality = 0.0_dp
         quality_ema = 0.0_dp
         residual_quality = 0.0_dp
         residual_quality_ema = 1.0_dp
         drop_sum = 0.0_dp
         measured_sum = 0.0_dp

      END SUBROUTINE reset_selected_calibration

   END SUBROUTINE qs_ot_reset_accepted_model_calibration

! **************************************************************************************************
!> \brief allocates the data in qs_ot_env, for a calculation with fm_struct_ref
!>        ortho_k allows for specifying an additional orthogonal subspace (i.e. c will
!>        be kept orthogonal provided c0 was, used in qs_ot_eigensolver)
!> \param qs_ot_env ...
!> \param matrix_s ...
!> \param fm_struct_ref ...
!> \param ortho_k ...
! **************************************************************************************************
   SUBROUTINE qs_ot_allocate(qs_ot_env, matrix_s, fm_struct_ref, ortho_k)
      TYPE(qs_ot_type)                                   :: qs_ot_env
      TYPE(dbcsr_type), POINTER                          :: matrix_s
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_ref
      INTEGER, OPTIONAL                                  :: ortho_k

      INTEGER                                            :: i, k, m_diis, my_ortho_k, n, ncoef
      TYPE(cp_blacs_env_type), POINTER                   :: context
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL cite_reference(VandeVondele2003)

      CPASSERT(.NOT. qs_ot_env%state_allocated)
      qs_ot_env%has_complex_kpoint_state = .FALSE.
      NULLIFY (qs_ot_env%preconditioner)
      NULLIFY (qs_ot_env%matrix_psc0)
      NULLIFY (qs_ot_env%para_env)
      NULLIFY (qs_ot_env%blacs_env)

      CALL cp_fm_struct_get(fm_struct_ref, nrow_global=n, ncol_global=k, &
                            para_env=para_env, context=context)

      qs_ot_env%para_env => para_env
      qs_ot_env%blacs_env => context
      CALL para_env%retain()
      CALL context%retain()

      IF (PRESENT(ortho_k)) THEN
         my_ortho_k = ortho_k
      ELSE
         my_ortho_k = k
      END IF

      m_diis = qs_ot_env%settings%diis_m

      qs_ot_env%use_gx_old = .FALSE.
      qs_ot_env%use_dx = .FALSE.

      SELECT CASE (qs_ot_env%settings%ot_method)
      CASE ("SD")
         ! nothing
      CASE ("CG")
         qs_ot_env%use_gx_old = .TRUE.
         qs_ot_env%use_dx = .TRUE.
      CASE ("DIIS", "BROY")
         IF (m_diis < 1) CPABORT("m_diis less than one")
      CASE DEFAULT
         CPABORT("Unknown option")
      END SELECT

      IF (qs_ot_env%settings%ot_method == "DIIS" .OR. &
          qs_ot_env%settings%ot_method == "BROY") THEN
         ALLOCATE (qs_ot_env%ls_diis(m_diis + 1, m_diis + 1))
         qs_ot_env%ls_diis = 0.0_dp
         ALLOCATE (qs_ot_env%lss_diis(m_diis + 1, m_diis + 1))
         ALLOCATE (qs_ot_env%c_diis(m_diis + 1))
         ALLOCATE (qs_ot_env%c_broy(m_diis))
         ALLOCATE (qs_ot_env%energy_h(m_diis))
         ALLOCATE (qs_ot_env%ipivot(m_diis + 1))
      END IF

      ALLOCATE (qs_ot_env%evals(k))
      ALLOCATE (qs_ot_env%dum(k))

      NULLIFY (qs_ot_env%matrix_os)
      NULLIFY (qs_ot_env%matrix_buf1_ortho)
      NULLIFY (qs_ot_env%matrix_buf2_ortho)
      NULLIFY (qs_ot_env%matrix_p)
      NULLIFY (qs_ot_env%matrix_r)
      NULLIFY (qs_ot_env%matrix_sinp)
      NULLIFY (qs_ot_env%matrix_cosp)
      NULLIFY (qs_ot_env%matrix_sinp_b)
      NULLIFY (qs_ot_env%matrix_cosp_b)
      NULLIFY (qs_ot_env%matrix_buf1)
      NULLIFY (qs_ot_env%matrix_buf2)
      NULLIFY (qs_ot_env%matrix_buf3)
      NULLIFY (qs_ot_env%matrix_buf4)
      NULLIFY (qs_ot_env%matrix_c0)
      NULLIFY (qs_ot_env%matrix_sc0)
      NULLIFY (qs_ot_env%matrix_c0_im)
      NULLIFY (qs_ot_env%matrix_sc0_im)
      NULLIFY (qs_ot_env%matrix_x)
      NULLIFY (qs_ot_env%matrix_sx)
      NULLIFY (qs_ot_env%matrix_x_im)
      NULLIFY (qs_ot_env%matrix_sx_im)
      NULLIFY (qs_ot_env%matrix_ref_orth)
      NULLIFY (qs_ot_env%matrix_ref_orth_im)
      NULLIFY (qs_ot_env%matrix_gx)
      NULLIFY (qs_ot_env%matrix_gx_im)
      NULLIFY (qs_ot_env%matrix_gx_old)
      NULLIFY (qs_ot_env%matrix_gx_old_im)
      NULLIFY (qs_ot_env%matrix_dx)
      NULLIFY (qs_ot_env%matrix_dx_im)
      NULLIFY (qs_ot_env%buf1_k_k_nosym)
      NULLIFY (qs_ot_env%buf2_k_k_nosym)
      NULLIFY (qs_ot_env%buf3_k_k_nosym)
      NULLIFY (qs_ot_env%buf4_k_k_nosym)
      NULLIFY (qs_ot_env%buf1_k_k_sym)
      NULLIFY (qs_ot_env%buf2_k_k_sym)
      NULLIFY (qs_ot_env%buf3_k_k_sym)
      NULLIFY (qs_ot_env%buf4_k_k_sym)
      NULLIFY (qs_ot_env%buf1_n_k)
      NULLIFY (qs_ot_env%buf1_n_k_dp)
      NULLIFY (qs_ot_env%p_k_k_sym)

      ! COMMON MATRICES
      CALL dbcsr_init_p(qs_ot_env%matrix_c0)
      CALL cp_dbcsr_m_by_n_from_row_template(qs_ot_env%matrix_c0, template=matrix_s, n=k, &
                                             sym=dbcsr_type_no_symmetry)

      CALL dbcsr_init_p(qs_ot_env%matrix_sc0)
      CALL cp_dbcsr_m_by_n_from_row_template(qs_ot_env%matrix_sc0, template=matrix_s, n=my_ortho_k, &
                                             sym=dbcsr_type_no_symmetry)

      CALL dbcsr_init_p(qs_ot_env%matrix_x)
      CALL cp_dbcsr_m_by_n_from_row_template(qs_ot_env%matrix_x, template=matrix_s, n=k, &
                                             sym=dbcsr_type_no_symmetry)

      CALL dbcsr_init_p(qs_ot_env%matrix_sx)
      CALL cp_dbcsr_m_by_n_from_row_template(qs_ot_env%matrix_sx, template=matrix_s, n=k, &
                                             sym=dbcsr_type_no_symmetry)

      CALL dbcsr_init_p(qs_ot_env%matrix_gx)
      CALL cp_dbcsr_m_by_n_from_row_template(qs_ot_env%matrix_gx, template=matrix_s, n=k, &
                                             sym=dbcsr_type_no_symmetry)

      IF (qs_ot_env%use_dx) THEN
         CALL dbcsr_init_p(qs_ot_env%matrix_dx)
         CALL cp_dbcsr_m_by_n_from_row_template(qs_ot_env%matrix_dx, template=matrix_s, n=k, &
                                                sym=dbcsr_type_no_symmetry)
      END IF

      IF (qs_ot_env%use_gx_old) THEN
         CALL dbcsr_init_p(qs_ot_env%matrix_gx_old)
         CALL cp_dbcsr_m_by_n_from_row_template(qs_ot_env%matrix_gx_old, template=matrix_s, n=k, &
                                                sym=dbcsr_type_no_symmetry)
      END IF

      SELECT CASE (qs_ot_env%settings%ot_algorithm)
      CASE ("TOD")
         CALL dbcsr_init_p(qs_ot_env%matrix_p)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%matrix_p, template=matrix_s, m=k, n=k, &
                                            sym=dbcsr_type_no_symmetry)

         CALL dbcsr_init_p(qs_ot_env%matrix_r)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%matrix_r, template=matrix_s, m=k, n=k, &
                                            sym=dbcsr_type_no_symmetry)

         CALL dbcsr_init_p(qs_ot_env%matrix_sinp)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%matrix_sinp, template=matrix_s, m=k, n=k, &
                                            sym=dbcsr_type_no_symmetry)

         CALL dbcsr_init_p(qs_ot_env%matrix_cosp)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%matrix_cosp, template=matrix_s, m=k, n=k, &
                                            sym=dbcsr_type_no_symmetry)

         CALL dbcsr_init_p(qs_ot_env%matrix_sinp_b)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%matrix_sinp_b, template=matrix_s, m=k, n=k, &
                                            sym=dbcsr_type_no_symmetry)

         CALL dbcsr_init_p(qs_ot_env%matrix_cosp_b)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%matrix_cosp_b, template=matrix_s, m=k, n=k, &
                                            sym=dbcsr_type_no_symmetry)

         CALL dbcsr_init_p(qs_ot_env%matrix_buf1)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%matrix_buf1, template=matrix_s, m=k, n=k, &
                                            sym=dbcsr_type_no_symmetry)

         CALL dbcsr_init_p(qs_ot_env%matrix_buf2)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%matrix_buf2, template=matrix_s, m=k, n=k, &
                                            sym=dbcsr_type_no_symmetry)

         CALL dbcsr_init_p(qs_ot_env%matrix_buf3)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%matrix_buf3, template=matrix_s, m=k, n=k, &
                                            sym=dbcsr_type_no_symmetry)

         CALL dbcsr_init_p(qs_ot_env%matrix_buf4)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%matrix_buf4, template=matrix_s, m=k, n=k, &
                                            sym=dbcsr_type_no_symmetry)

         CALL dbcsr_init_p(qs_ot_env%matrix_os)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%matrix_os, template=matrix_s, m=my_ortho_k, n=my_ortho_k, &
                                            sym=dbcsr_type_no_symmetry)

         CALL dbcsr_init_p(qs_ot_env%matrix_buf1_ortho)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%matrix_buf1_ortho, template=matrix_s, m=my_ortho_k, n=k, &
                                            sym=dbcsr_type_no_symmetry)

         CALL dbcsr_init_p(qs_ot_env%matrix_buf2_ortho)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%matrix_buf2_ortho, template=matrix_s, m=my_ortho_k, n=k, &
                                            sym=dbcsr_type_no_symmetry)

      CASE ("REF")
         CALL dbcsr_init_p(qs_ot_env%buf1_k_k_nosym)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%buf1_k_k_nosym, template=matrix_s, m=k, n=k, &
                                            sym=dbcsr_type_no_symmetry)

         CALL dbcsr_init_p(qs_ot_env%buf2_k_k_nosym)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%buf2_k_k_nosym, template=matrix_s, m=k, n=k, &
                                            sym=dbcsr_type_no_symmetry)

         CALL dbcsr_init_p(qs_ot_env%buf3_k_k_nosym)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%buf3_k_k_nosym, template=matrix_s, m=k, n=k, &
                                            sym=dbcsr_type_no_symmetry)

         CALL dbcsr_init_p(qs_ot_env%buf4_k_k_nosym)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%buf4_k_k_nosym, template=matrix_s, m=k, n=k, &
                                            sym=dbcsr_type_no_symmetry)

         ! It claims to be symmetric but to avoid dbcsr confusion nonsymmetric is kept
         CALL dbcsr_init_p(qs_ot_env%buf1_k_k_sym)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%buf1_k_k_sym, template=matrix_s, m=k, n=k, &
                                            sym=dbcsr_type_no_symmetry)

         CALL dbcsr_init_p(qs_ot_env%buf2_k_k_sym)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%buf2_k_k_sym, template=matrix_s, m=k, n=k, &
                                            sym=dbcsr_type_no_symmetry)

         CALL dbcsr_init_p(qs_ot_env%buf3_k_k_sym)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%buf3_k_k_sym, template=matrix_s, m=k, n=k, &
                                            sym=dbcsr_type_no_symmetry)
         !
         CALL dbcsr_init_p(qs_ot_env%buf4_k_k_sym)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%buf4_k_k_sym, template=matrix_s, m=k, n=k, &
                                            sym=dbcsr_type_no_symmetry)
         !
         CALL dbcsr_init_p(qs_ot_env%p_k_k_sym)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%p_k_k_sym, template=matrix_s, m=k, n=k, &
                                            sym=dbcsr_type_no_symmetry)
         !
         CALL dbcsr_init_p(qs_ot_env%buf1_n_k)
         CALL cp_dbcsr_m_by_n_from_row_template(qs_ot_env%buf1_n_k, template=matrix_s, n=k, &
                                                sym=dbcsr_type_no_symmetry)
         !
         CALL dbcsr_init_p(qs_ot_env%matrix_buf1)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%matrix_buf1, template=matrix_s, m=k, n=k, &
                                            sym=dbcsr_type_no_symmetry)

      END SELECT

      IF (qs_ot_env%settings%ot_method == "DIIS" .OR. &
          qs_ot_env%settings%ot_method == "BROY") THEN
         NULLIFY (qs_ot_env%matrix_h_e)
         NULLIFY (qs_ot_env%matrix_h_x)
         CALL dbcsr_allocate_matrix_set(qs_ot_env%matrix_h_e, m_diis)
         CALL dbcsr_allocate_matrix_set(qs_ot_env%matrix_h_x, m_diis)
         DO i = 1, m_diis
            CALL dbcsr_init_p(qs_ot_env%matrix_h_x(i)%matrix)
            CALL cp_dbcsr_m_by_n_from_row_template(qs_ot_env%matrix_h_x(i)%matrix, template=matrix_s, n=k, &
                                                   sym=dbcsr_type_no_symmetry)

            CALL dbcsr_init_p(qs_ot_env%matrix_h_e(i)%matrix)
            CALL cp_dbcsr_m_by_n_from_row_template(qs_ot_env%matrix_h_e(i)%matrix, template=matrix_s, n=k, &
                                                   sym=dbcsr_type_no_symmetry)
         END DO
      END IF

      NULLIFY (qs_ot_env%rot_mat_u, qs_ot_env%rot_mat_u_im, qs_ot_env%rot_mat_x, qs_ot_env%rot_mat_x_im, &
               qs_ot_env%rot_mat_h_e, qs_ot_env%rot_mat_h_x, qs_ot_env%rot_mat_h_e_im, &
               qs_ot_env%rot_mat_h_x_im, qs_ot_env%rot_mat_gx, qs_ot_env%rot_mat_gx_im, &
               qs_ot_env%rot_mat_pgx, qs_ot_env%rot_mat_pgx_im, &
               qs_ot_env%rot_mat_gx_old, qs_ot_env%rot_mat_gx_old_im, qs_ot_env%rot_mat_dx, &
               qs_ot_env%rot_mat_dx_im, qs_ot_env%rot_mat_mermin_g0, qs_ot_env%rot_mat_mermin_g0_im, &
               qs_ot_env%rot_mat_mermin_s, qs_ot_env%rot_mat_mermin_s_im, &
               qs_ot_env%rot_mat_mermin_y, qs_ot_env%rot_mat_mermin_y_im, &
               qs_ot_env%rot_mat_mermin_s_previous, qs_ot_env%rot_mat_mermin_s_previous_im, &
               qs_ot_env%rot_mat_mermin_y_previous, qs_ot_env%rot_mat_mermin_y_previous_im, &
               qs_ot_env%rot_mat_mermin_h0_ref, qs_ot_env%rot_mat_mermin_h0_ref_im, &
               qs_ot_env%rot_mat_mermin_frozen_g1, qs_ot_env%rot_mat_mermin_frozen_g1_im, &
               qs_ot_env%rot_mat_mermin_frozen_reduced_y, &
               qs_ot_env%rot_mat_mermin_frozen_reduced_y_im, &
               qs_ot_env%rot_mat_mermin_closure_y, qs_ot_env%rot_mat_mermin_closure_y_im, &
               qs_ot_env%rot_mat_mermin_closure_y_previous, &
               qs_ot_env%rot_mat_mermin_closure_y_previous_im, &
               qs_ot_env%rot_mat_mermin_s_anchor, qs_ot_env%rot_mat_mermin_s_anchor_im, &
               qs_ot_env%rot_mat_mermin_y_anchor, qs_ot_env%rot_mat_mermin_y_anchor_im, &
               qs_ot_env%rot_mat_mermin_closure_y_anchor, &
               qs_ot_env%rot_mat_mermin_closure_y_anchor_im, &
               qs_ot_env%rot_mat_evals, qs_ot_env%rot_mat_dedu, &
               qs_ot_env%rot_mat_dedu_im, qs_ot_env%rot_mat_chc, qs_ot_env%rot_mat_chc_im, &
               qs_ot_env%rot_mat_evec_re, qs_ot_env%rot_mat_evec_im)

      IF (qs_ot_env%settings%do_rotation) THEN
         CALL dbcsr_init_p(qs_ot_env%rot_mat_u)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%rot_mat_u, template=matrix_s, m=k, n=k, &
                                            sym=dbcsr_type_no_symmetry)

         CALL dbcsr_init_p(qs_ot_env%rot_mat_x)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%rot_mat_x, template=matrix_s, m=k, n=k, &
                                            sym=dbcsr_type_no_symmetry)

         CALL dbcsr_init_p(qs_ot_env%rot_mat_dedu)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%rot_mat_dedu, template=matrix_s, m=k, n=k, &
                                            sym=dbcsr_type_no_symmetry)

         CALL dbcsr_init_p(qs_ot_env%rot_mat_chc)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%rot_mat_chc, template=matrix_s, m=k, n=k, &
                                            sym=dbcsr_type_no_symmetry)

         IF (qs_ot_env%settings%ot_method == "DIIS") THEN
            CALL dbcsr_allocate_matrix_set(qs_ot_env%rot_mat_h_e, m_diis)
            CALL dbcsr_allocate_matrix_set(qs_ot_env%rot_mat_h_x, m_diis)
            DO i = 1, m_diis
               CALL dbcsr_init_p(qs_ot_env%rot_mat_h_e(i)%matrix)
               CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%rot_mat_h_e(i)%matrix, template=matrix_s, m=k, n=k, &
                                                  sym=dbcsr_type_no_symmetry)

               CALL dbcsr_init_p(qs_ot_env%rot_mat_h_x(i)%matrix)
               CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%rot_mat_h_x(i)%matrix, template=matrix_s, m=k, n=k, &
                                                  sym=dbcsr_type_no_symmetry)
            END DO
         END IF

         ALLOCATE (qs_ot_env%rot_mat_evals(k))
         CALL dbcsr_init_p(qs_ot_env%rot_mat_evec_re)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%rot_mat_evec_re, template=matrix_s, m=k, n=k, &
                                            sym=dbcsr_type_no_symmetry)
         CALL dbcsr_init_p(qs_ot_env%rot_mat_evec_im)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%rot_mat_evec_im, template=matrix_s, m=k, n=k, &
                                            sym=dbcsr_type_no_symmetry)

         CALL dbcsr_init_p(qs_ot_env%rot_mat_gx)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%rot_mat_gx, template=matrix_s, m=k, n=k, &
                                            sym=dbcsr_type_no_symmetry)

         IF (qs_ot_env%use_gx_old) THEN
            CALL dbcsr_init_p(qs_ot_env%rot_mat_gx_old)
            CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%rot_mat_gx_old, template=matrix_s, m=k, n=k, &
                                               sym=dbcsr_type_no_symmetry)
         END IF

         IF (qs_ot_env%use_dx) THEN
            CALL dbcsr_init_p(qs_ot_env%rot_mat_dx)
            CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%rot_mat_dx, template=matrix_s, m=k, n=k, &
                                               sym=dbcsr_type_no_symmetry)
         END IF

      END IF

      IF (qs_ot_env%settings%do_ener) THEN
         ncoef = k
         ALLOCATE (qs_ot_env%ener_x(ncoef))
         ALLOCATE (qs_ot_env%ener_x_trust_weight(ncoef))
         ALLOCATE (qs_ot_env%matrix_occ_precond_scale(ncoef))
         qs_ot_env%ener_x_trust_weight(:) = 0.0_dp
         qs_ot_env%matrix_occ_precond_scale(:) = 1.0_dp

         IF (qs_ot_env%settings%ot_method == "DIIS") THEN
            ALLOCATE (qs_ot_env%ener_h_e(m_diis, ncoef))
            ALLOCATE (qs_ot_env%ener_h_x(m_diis, ncoef))
         END IF

         ALLOCATE (qs_ot_env%ener_gx(ncoef))

         IF (qs_ot_env%use_gx_old) THEN
            ALLOCATE (qs_ot_env%ener_gx_old(ncoef))
         END IF

         IF (qs_ot_env%use_dx) THEN
            ALLOCATE (qs_ot_env%ener_dx(ncoef))
            qs_ot_env%ener_dx = 0.0_dp
         END IF
      END IF

      qs_ot_env%state_allocated = .TRUE.

   END SUBROUTINE qs_ot_allocate

! **************************************************************************************************
!> \brief deallocates data
!> \param qs_ot_env ...
! **************************************************************************************************
   SUBROUTINE qs_ot_destroy(qs_ot_env)
      TYPE(qs_ot_type)                                   :: qs_ot_env

      IF (.NOT. qs_ot_env%state_allocated) RETURN

      CALL mp_para_env_release(qs_ot_env%para_env)
      CALL cp_blacs_env_release(qs_ot_env%blacs_env)

      IF (ASSOCIATED(qs_ot_env%evals)) DEALLOCATE (qs_ot_env%evals)
      IF (ASSOCIATED(qs_ot_env%dum)) DEALLOCATE (qs_ot_env%dum)

      IF (ASSOCIATED(qs_ot_env%matrix_os)) CALL dbcsr_release_p(qs_ot_env%matrix_os)
      IF (ASSOCIATED(qs_ot_env%matrix_p)) CALL dbcsr_release_p(qs_ot_env%matrix_p)
      IF (ASSOCIATED(qs_ot_env%matrix_cosp)) CALL dbcsr_release_p(qs_ot_env%matrix_cosp)
      IF (ASSOCIATED(qs_ot_env%matrix_sinp)) CALL dbcsr_release_p(qs_ot_env%matrix_sinp)
      IF (ASSOCIATED(qs_ot_env%matrix_r)) CALL dbcsr_release_p(qs_ot_env%matrix_r)
      IF (ASSOCIATED(qs_ot_env%matrix_cosp_b)) CALL dbcsr_release_p(qs_ot_env%matrix_cosp_b)
      IF (ASSOCIATED(qs_ot_env%matrix_sinp_b)) CALL dbcsr_release_p(qs_ot_env%matrix_sinp_b)
      IF (ASSOCIATED(qs_ot_env%matrix_buf1)) CALL dbcsr_release_p(qs_ot_env%matrix_buf1)
      IF (ASSOCIATED(qs_ot_env%matrix_buf2)) CALL dbcsr_release_p(qs_ot_env%matrix_buf2)
      IF (ASSOCIATED(qs_ot_env%matrix_buf3)) CALL dbcsr_release_p(qs_ot_env%matrix_buf3)
      IF (ASSOCIATED(qs_ot_env%matrix_buf4)) CALL dbcsr_release_p(qs_ot_env%matrix_buf4)
      IF (ASSOCIATED(qs_ot_env%matrix_buf1_ortho)) CALL dbcsr_release_p(qs_ot_env%matrix_buf1_ortho)
      IF (ASSOCIATED(qs_ot_env%matrix_buf2_ortho)) CALL dbcsr_release_p(qs_ot_env%matrix_buf2_ortho)
      IF (ASSOCIATED(qs_ot_env%matrix_c0)) CALL dbcsr_release_p(qs_ot_env%matrix_c0)
      IF (ASSOCIATED(qs_ot_env%matrix_sc0)) CALL dbcsr_release_p(qs_ot_env%matrix_sc0)
      IF (ASSOCIATED(qs_ot_env%matrix_psc0)) CALL dbcsr_release_p(qs_ot_env%matrix_psc0)
      IF (ASSOCIATED(qs_ot_env%matrix_c0_im)) CALL dbcsr_release_p(qs_ot_env%matrix_c0_im)
      IF (ASSOCIATED(qs_ot_env%matrix_sc0_im)) CALL dbcsr_release_p(qs_ot_env%matrix_sc0_im)
      IF (ASSOCIATED(qs_ot_env%matrix_x)) CALL dbcsr_release_p(qs_ot_env%matrix_x)
      IF (ASSOCIATED(qs_ot_env%matrix_sx)) CALL dbcsr_release_p(qs_ot_env%matrix_sx)
      IF (ASSOCIATED(qs_ot_env%matrix_x_im)) CALL dbcsr_release_p(qs_ot_env%matrix_x_im)
      IF (ASSOCIATED(qs_ot_env%matrix_sx_im)) CALL dbcsr_release_p(qs_ot_env%matrix_sx_im)
      IF (ASSOCIATED(qs_ot_env%matrix_ref_orth)) CALL dbcsr_release_p(qs_ot_env%matrix_ref_orth)
      IF (ASSOCIATED(qs_ot_env%matrix_ref_orth_im)) CALL dbcsr_release_p(qs_ot_env%matrix_ref_orth_im)
      IF (ASSOCIATED(qs_ot_env%matrix_gx)) CALL dbcsr_release_p(qs_ot_env%matrix_gx)
      IF (ASSOCIATED(qs_ot_env%matrix_gx_im)) CALL dbcsr_release_p(qs_ot_env%matrix_gx_im)
      IF (ASSOCIATED(qs_ot_env%matrix_dx)) CALL dbcsr_release_p(qs_ot_env%matrix_dx)
      IF (ASSOCIATED(qs_ot_env%matrix_dx_im)) CALL dbcsr_release_p(qs_ot_env%matrix_dx_im)
      IF (ASSOCIATED(qs_ot_env%matrix_gx_old)) CALL dbcsr_release_p(qs_ot_env%matrix_gx_old)
      IF (ASSOCIATED(qs_ot_env%matrix_gx_old_im)) CALL dbcsr_release_p(qs_ot_env%matrix_gx_old_im)
      IF (ASSOCIATED(qs_ot_env%buf1_k_k_nosym)) CALL dbcsr_release_p(qs_ot_env%buf1_k_k_nosym)
      IF (ASSOCIATED(qs_ot_env%buf2_k_k_nosym)) CALL dbcsr_release_p(qs_ot_env%buf2_k_k_nosym)
      IF (ASSOCIATED(qs_ot_env%buf3_k_k_nosym)) CALL dbcsr_release_p(qs_ot_env%buf3_k_k_nosym)
      IF (ASSOCIATED(qs_ot_env%buf4_k_k_nosym)) CALL dbcsr_release_p(qs_ot_env%buf4_k_k_nosym)
      IF (ASSOCIATED(qs_ot_env%p_k_k_sym)) CALL dbcsr_release_p(qs_ot_env%p_k_k_sym)
      IF (ASSOCIATED(qs_ot_env%buf1_k_k_sym)) CALL dbcsr_release_p(qs_ot_env%buf1_k_k_sym)
      IF (ASSOCIATED(qs_ot_env%buf2_k_k_sym)) CALL dbcsr_release_p(qs_ot_env%buf2_k_k_sym)
      IF (ASSOCIATED(qs_ot_env%buf3_k_k_sym)) CALL dbcsr_release_p(qs_ot_env%buf3_k_k_sym)
      IF (ASSOCIATED(qs_ot_env%buf4_k_k_sym)) CALL dbcsr_release_p(qs_ot_env%buf4_k_k_sym)
      IF (ASSOCIATED(qs_ot_env%buf1_n_k)) CALL dbcsr_release_p(qs_ot_env%buf1_n_k)
      IF (ASSOCIATED(qs_ot_env%buf1_n_k_dp)) CALL dbcsr_release_p(qs_ot_env%buf1_n_k_dp)

      IF (qs_ot_env%settings%ot_method == "DIIS" .OR. &
          qs_ot_env%settings%ot_method == "BROY") THEN
         CALL dbcsr_deallocate_matrix_set(qs_ot_env%matrix_h_x)
         CALL dbcsr_deallocate_matrix_set(qs_ot_env%matrix_h_e)
         IF (ASSOCIATED(qs_ot_env%matrix_h_x_im)) CALL dbcsr_deallocate_matrix_set(qs_ot_env%matrix_h_x_im)
         IF (ASSOCIATED(qs_ot_env%matrix_h_e_im)) CALL dbcsr_deallocate_matrix_set(qs_ot_env%matrix_h_e_im)
         DEALLOCATE (qs_ot_env%ls_diis)
         DEALLOCATE (qs_ot_env%lss_diis)
         DEALLOCATE (qs_ot_env%c_diis)
         DEALLOCATE (qs_ot_env%c_broy)
         DEALLOCATE (qs_ot_env%energy_h)
         DEALLOCATE (qs_ot_env%ipivot)
      END IF

      IF (qs_ot_env%settings%do_rotation) THEN

         IF (ASSOCIATED(qs_ot_env%rot_mat_u)) CALL dbcsr_release_p(qs_ot_env%rot_mat_u)
         IF (ASSOCIATED(qs_ot_env%rot_mat_u_im)) CALL dbcsr_release_p(qs_ot_env%rot_mat_u_im)
         IF (ASSOCIATED(qs_ot_env%rot_mat_x)) CALL dbcsr_release_p(qs_ot_env%rot_mat_x)
         IF (ASSOCIATED(qs_ot_env%rot_mat_x_im)) CALL dbcsr_release_p(qs_ot_env%rot_mat_x_im)
         IF (ASSOCIATED(qs_ot_env%rot_mat_dedu)) CALL dbcsr_release_p(qs_ot_env%rot_mat_dedu)
         IF (ASSOCIATED(qs_ot_env%rot_mat_dedu_im)) CALL dbcsr_release_p(qs_ot_env%rot_mat_dedu_im)
         IF (ASSOCIATED(qs_ot_env%rot_mat_chc)) CALL dbcsr_release_p(qs_ot_env%rot_mat_chc)
         IF (ASSOCIATED(qs_ot_env%rot_mat_chc_im)) CALL dbcsr_release_p(qs_ot_env%rot_mat_chc_im)

         IF (qs_ot_env%settings%ot_method == "DIIS") THEN
            CALL dbcsr_deallocate_matrix_set(qs_ot_env%rot_mat_h_x)
            CALL dbcsr_deallocate_matrix_set(qs_ot_env%rot_mat_h_e)
            IF (ASSOCIATED(qs_ot_env%rot_mat_h_x_im)) CALL dbcsr_deallocate_matrix_set(qs_ot_env%rot_mat_h_x_im)
            IF (ASSOCIATED(qs_ot_env%rot_mat_h_e_im)) CALL dbcsr_deallocate_matrix_set(qs_ot_env%rot_mat_h_e_im)
         END IF

         IF (ASSOCIATED(qs_ot_env%rot_mat_evals)) DEALLOCATE (qs_ot_env%rot_mat_evals)
         IF (ASSOCIATED(qs_ot_env%rot_mat_evec_re)) CALL dbcsr_release_p(qs_ot_env%rot_mat_evec_re)
         IF (ASSOCIATED(qs_ot_env%rot_mat_evec_im)) CALL dbcsr_release_p(qs_ot_env%rot_mat_evec_im)
         IF (ASSOCIATED(qs_ot_env%rot_mat_gx)) CALL dbcsr_release_p(qs_ot_env%rot_mat_gx)
         IF (ASSOCIATED(qs_ot_env%rot_mat_gx_im)) CALL dbcsr_release_p(qs_ot_env%rot_mat_gx_im)
         IF (ASSOCIATED(qs_ot_env%rot_mat_pgx)) CALL dbcsr_release_p(qs_ot_env%rot_mat_pgx)
         IF (ASSOCIATED(qs_ot_env%rot_mat_pgx_im)) CALL dbcsr_release_p(qs_ot_env%rot_mat_pgx_im)
         IF (ASSOCIATED(qs_ot_env%rot_mat_gx_old)) CALL dbcsr_release_p(qs_ot_env%rot_mat_gx_old)
         IF (ASSOCIATED(qs_ot_env%rot_mat_gx_old_im)) CALL dbcsr_release_p(qs_ot_env%rot_mat_gx_old_im)
         IF (ASSOCIATED(qs_ot_env%rot_mat_dx)) CALL dbcsr_release_p(qs_ot_env%rot_mat_dx)
         IF (ASSOCIATED(qs_ot_env%rot_mat_dx_im)) CALL dbcsr_release_p(qs_ot_env%rot_mat_dx_im)
         IF (ASSOCIATED(qs_ot_env%rot_mat_mermin_g0)) CALL dbcsr_release_p(qs_ot_env%rot_mat_mermin_g0)
         IF (ASSOCIATED(qs_ot_env%rot_mat_mermin_g0_im)) CALL dbcsr_release_p(qs_ot_env%rot_mat_mermin_g0_im)
         IF (ASSOCIATED(qs_ot_env%rot_mat_mermin_s)) CALL dbcsr_release_p(qs_ot_env%rot_mat_mermin_s)
         IF (ASSOCIATED(qs_ot_env%rot_mat_mermin_s_im)) CALL dbcsr_release_p(qs_ot_env%rot_mat_mermin_s_im)
         IF (ASSOCIATED(qs_ot_env%rot_mat_mermin_y)) CALL dbcsr_release_p(qs_ot_env%rot_mat_mermin_y)
         IF (ASSOCIATED(qs_ot_env%rot_mat_mermin_y_im)) CALL dbcsr_release_p(qs_ot_env%rot_mat_mermin_y_im)
         IF (ASSOCIATED(qs_ot_env%rot_mat_mermin_s_previous)) THEN
            CALL dbcsr_release_p(qs_ot_env%rot_mat_mermin_s_previous)
         END IF
         IF (ASSOCIATED(qs_ot_env%rot_mat_mermin_s_previous_im)) THEN
            CALL dbcsr_release_p(qs_ot_env%rot_mat_mermin_s_previous_im)
         END IF
         IF (ASSOCIATED(qs_ot_env%rot_mat_mermin_y_previous)) THEN
            CALL dbcsr_release_p(qs_ot_env%rot_mat_mermin_y_previous)
         END IF
         IF (ASSOCIATED(qs_ot_env%rot_mat_mermin_y_previous_im)) THEN
            CALL dbcsr_release_p(qs_ot_env%rot_mat_mermin_y_previous_im)
         END IF
         IF (ASSOCIATED(qs_ot_env%rot_mat_mermin_h0_ref)) CALL dbcsr_release_p(qs_ot_env%rot_mat_mermin_h0_ref)
         IF (ASSOCIATED(qs_ot_env%rot_mat_mermin_h0_ref_im)) CALL dbcsr_release_p(qs_ot_env%rot_mat_mermin_h0_ref_im)
         IF (ASSOCIATED(qs_ot_env%rot_mat_mermin_frozen_g1)) CALL dbcsr_release_p(qs_ot_env%rot_mat_mermin_frozen_g1)
         IF (ASSOCIATED(qs_ot_env%rot_mat_mermin_frozen_g1_im)) THEN
            CALL dbcsr_release_p(qs_ot_env%rot_mat_mermin_frozen_g1_im)
         END IF
         IF (ASSOCIATED(qs_ot_env%rot_mat_mermin_frozen_reduced_y)) THEN
            CALL dbcsr_release_p(qs_ot_env%rot_mat_mermin_frozen_reduced_y)
         END IF
         IF (ASSOCIATED(qs_ot_env%rot_mat_mermin_frozen_reduced_y_im)) THEN
            CALL dbcsr_release_p(qs_ot_env%rot_mat_mermin_frozen_reduced_y_im)
         END IF
         IF (ASSOCIATED(qs_ot_env%rot_mat_mermin_closure_y)) CALL dbcsr_release_p(qs_ot_env%rot_mat_mermin_closure_y)
         IF (ASSOCIATED(qs_ot_env%rot_mat_mermin_closure_y_im)) THEN
            CALL dbcsr_release_p(qs_ot_env%rot_mat_mermin_closure_y_im)
         END IF
         IF (ASSOCIATED(qs_ot_env%rot_mat_mermin_closure_y_previous)) THEN
            CALL dbcsr_release_p(qs_ot_env%rot_mat_mermin_closure_y_previous)
         END IF
         IF (ASSOCIATED(qs_ot_env%rot_mat_mermin_closure_y_previous_im)) THEN
            CALL dbcsr_release_p(qs_ot_env%rot_mat_mermin_closure_y_previous_im)
         END IF
         IF (ASSOCIATED(qs_ot_env%rot_mat_mermin_s_anchor)) THEN
            CALL dbcsr_release_p(qs_ot_env%rot_mat_mermin_s_anchor)
         END IF
         IF (ASSOCIATED(qs_ot_env%rot_mat_mermin_s_anchor_im)) THEN
            CALL dbcsr_release_p(qs_ot_env%rot_mat_mermin_s_anchor_im)
         END IF
         IF (ASSOCIATED(qs_ot_env%rot_mat_mermin_y_anchor)) THEN
            CALL dbcsr_release_p(qs_ot_env%rot_mat_mermin_y_anchor)
         END IF
         IF (ASSOCIATED(qs_ot_env%rot_mat_mermin_y_anchor_im)) THEN
            CALL dbcsr_release_p(qs_ot_env%rot_mat_mermin_y_anchor_im)
         END IF
         IF (ASSOCIATED(qs_ot_env%rot_mat_mermin_closure_y_anchor)) THEN
            CALL dbcsr_release_p(qs_ot_env%rot_mat_mermin_closure_y_anchor)
         END IF
         IF (ASSOCIATED(qs_ot_env%rot_mat_mermin_closure_y_anchor_im)) THEN
            CALL dbcsr_release_p(qs_ot_env%rot_mat_mermin_closure_y_anchor_im)
         END IF
      END IF

      IF (qs_ot_env%settings%do_ener) THEN
         IF (ASSOCIATED(qs_ot_env%ener_x)) DEALLOCATE (qs_ot_env%ener_x)
         IF (ASSOCIATED(qs_ot_env%ener_x_trust_weight)) DEALLOCATE (qs_ot_env%ener_x_trust_weight)
         IF (ASSOCIATED(qs_ot_env%matrix_occ_precond_scale)) THEN
            DEALLOCATE (qs_ot_env%matrix_occ_precond_scale)
         END IF
         IF (ASSOCIATED(qs_ot_env%ener_gx)) DEALLOCATE (qs_ot_env%ener_gx)
         IF (ASSOCIATED(qs_ot_env%ener_pgx)) DEALLOCATE (qs_ot_env%ener_pgx)
         IF (ASSOCIATED(qs_ot_env%ener_mermin_g0)) DEALLOCATE (qs_ot_env%ener_mermin_g0)
         IF (ASSOCIATED(qs_ot_env%ener_mermin_x0)) DEALLOCATE (qs_ot_env%ener_mermin_x0)
         IF (ASSOCIATED(qs_ot_env%ener_mermin_s)) DEALLOCATE (qs_ot_env%ener_mermin_s)
         IF (ASSOCIATED(qs_ot_env%ener_mermin_y)) DEALLOCATE (qs_ot_env%ener_mermin_y)
         IF (ASSOCIATED(qs_ot_env%ener_mermin_s_previous)) DEALLOCATE (qs_ot_env%ener_mermin_s_previous)
         IF (ASSOCIATED(qs_ot_env%ener_mermin_y_previous)) DEALLOCATE (qs_ot_env%ener_mermin_y_previous)
         IF (ASSOCIATED(qs_ot_env%ener_mermin_frozen_diag)) DEALLOCATE (qs_ot_env%ener_mermin_frozen_diag)
         IF (ASSOCIATED(qs_ot_env%ener_mermin_frozen_g1)) DEALLOCATE (qs_ot_env%ener_mermin_frozen_g1)
         IF (ASSOCIATED(qs_ot_env%ener_mermin_closure_y)) DEALLOCATE (qs_ot_env%ener_mermin_closure_y)
         IF (ASSOCIATED(qs_ot_env%ener_mermin_closure_y_previous)) THEN
            DEALLOCATE (qs_ot_env%ener_mermin_closure_y_previous)
         END IF
         IF (ASSOCIATED(qs_ot_env%ener_mermin_s_anchor)) DEALLOCATE (qs_ot_env%ener_mermin_s_anchor)
         IF (ASSOCIATED(qs_ot_env%ener_mermin_y_anchor)) DEALLOCATE (qs_ot_env%ener_mermin_y_anchor)
         IF (ASSOCIATED(qs_ot_env%ener_mermin_closure_y_anchor)) THEN
            DEALLOCATE (qs_ot_env%ener_mermin_closure_y_anchor)
         END IF
         IF (qs_ot_env%settings%ot_method == "DIIS") THEN
            IF (ASSOCIATED(qs_ot_env%ener_h_x)) DEALLOCATE (qs_ot_env%ener_h_x)
            IF (ASSOCIATED(qs_ot_env%ener_h_e)) DEALLOCATE (qs_ot_env%ener_h_e)
         END IF
         IF (qs_ot_env%use_dx) THEN
            IF (ASSOCIATED(qs_ot_env%ener_dx)) DEALLOCATE (qs_ot_env%ener_dx)
         END IF
         IF (qs_ot_env%use_gx_old) THEN
            IF (ASSOCIATED(qs_ot_env%ener_gx_old)) DEALLOCATE (qs_ot_env%ener_gx_old)
         END IF
      END IF

      IF (ALLOCATED(qs_ot_env%orbital_schur_response_metric)) THEN
         DEALLOCATE (qs_ot_env%orbital_schur_response_metric)
      END IF
      IF (ALLOCATED(qs_ot_env%kpoint_density_pair_residual)) THEN
         DEALLOCATE (qs_ot_env%kpoint_density_pair_residual)
      END IF
      IF (ALLOCATED(qs_ot_env%kpoint_canonical_occupation)) THEN
         DEALLOCATE (qs_ot_env%kpoint_canonical_occupation)
      END IF
      IF (ALLOCATED(qs_ot_env%kpoint_canonical_eigenvalue)) THEN
         DEALLOCATE (qs_ot_env%kpoint_canonical_eigenvalue)
      END IF

      qs_ot_env%state_allocated = .FALSE.
      qs_ot_env%has_complex_kpoint_state = .FALSE.
      qs_ot_env%kpoint_occ_preconditioner_active = .FALSE.
      qs_ot_env%kpoint_occ_preconditioner_ready_samples = 0
      qs_ot_env%kpoint_subspace_accepted_samples = 0
      qs_ot_env%kpoint_occ_metric_residual_quality_ema = 0.0_dp
      qs_ot_env%orbital_schur_response_ready = .FALSE.
      qs_ot_env%kpoint_preconditioned_direction_ready = .FALSE.
      qs_ot_env%ener_direction_linear = .FALSE.

   END SUBROUTINE qs_ot_destroy

! **************************************************************************************************
!> \brief ...
!> \param settings ...
!> \param ot_section ...
!> \param output_unit ...
!> \param kpoint_default_irac Select IRAC when ALGORITHM was not explicitly given.
!> \param smearing_default_metallic_variables Select ROTATION and ENERGIES as a coupled default.
!> \param smearing_default_occupation_preconditioner Select the fixed-N preconditioner when omitted.
!> \param smearing_default_preconditioner Select the orbital preconditioner when omitted.
!> \param smearing_default_kpoint_full_all Select the complex k-point orbital preconditioner.
! **************************************************************************************************
   SUBROUTINE ot_readwrite_input(settings, ot_section, output_unit, kpoint_default_irac, &
                                 smearing_default_metallic_variables, &
                                 smearing_default_occupation_preconditioner, smearing_default_preconditioner, &
                                 smearing_default_kpoint_full_all)
      TYPE(qs_ot_settings_type)                          :: settings
      TYPE(section_vals_type), POINTER                   :: ot_section
      INTEGER, INTENT(IN)                                :: output_unit
      LOGICAL, INTENT(IN), OPTIONAL :: kpoint_default_irac, smearing_default_metallic_variables, &
         smearing_default_occupation_preconditioner, smearing_default_preconditioner, &
         smearing_default_kpoint_full_all

      CHARACTER(len=*), PARAMETER :: routineN = 'ot_readwrite_input'

      INTEGER                                            :: handle, ls_method, ot_algorithm, &
                                                            ot_method, ot_ortho_irac
      LOGICAL :: algorithm_explicit, energies_explicit, kpoint_smearing_context, &
         occupation_preconditioner_explicit, preconditioner_explicit, rotation_explicit, &
         selected_energies_default, selected_kpoint_default, &
         selected_occupation_preconditioner_default, selected_preconditioner_default, &
         selected_rotation_default, smearing_context

      CALL timeset(routineN, handle)

      ! choose algorithm
      selected_kpoint_default = .FALSE.
      selected_rotation_default = .FALSE.
      selected_energies_default = .FALSE.
      selected_occupation_preconditioner_default = .FALSE.
      selected_preconditioner_default = .FALSE.
      kpoint_smearing_context = .FALSE.
      smearing_context = .FALSE.
      IF (PRESENT(kpoint_default_irac)) kpoint_smearing_context = kpoint_default_irac
      IF (PRESENT(smearing_default_preconditioner)) THEN
         smearing_context = smearing_default_preconditioner
      END IF
      CALL section_vals_val_get(ot_section, "ALGORITHM", i_val=ot_algorithm, explicit=algorithm_explicit)
      IF (PRESENT(kpoint_default_irac)) THEN
         IF (kpoint_default_irac .AND. .NOT. algorithm_explicit) THEN
            ot_algorithm = ot_algo_irac
            selected_kpoint_default = .TRUE.
         END IF
      END IF
      SELECT CASE (ot_algorithm)
      CASE (ot_algo_taylor_or_diag)
         settings%ot_algorithm = "TOD"
      CASE (ot_algo_irac)
         CALL cite_reference(Weber2008)
         settings%ot_algorithm = "REF"
      CASE DEFAULT
         CPABORT("Value unknown")
      END SELECT

      ! irac input
      CALL section_vals_val_get(ot_section, "IRAC_DEGREE", i_val=settings%irac_degree)
      IF (settings%irac_degree < 2 .OR. settings%irac_degree > 4) THEN
         CPABORT("READ OT IRAC_DEGREE: Value unknown")
      END IF
      CALL section_vals_val_get(ot_section, "MAX_IRAC", i_val=settings%max_irac)
      IF (settings%max_irac < 1) THEN
         CPABORT("READ OT MAX_IRAC: VALUE MUST BE GREATER THAN ZERO")
      END IF
      CALL section_vals_val_get(ot_section, "EPS_IRAC_FILTER_MATRIX", r_val=settings%eps_irac_filter_matrix)
      CALL section_vals_val_get(ot_section, "EPS_IRAC", r_val=settings%eps_irac)
      IF (settings%eps_irac < 0.0_dp) THEN
         CPABORT("READ OT EPS_IRAC: VALUE MUST BE GREATER THAN ZERO")
      END IF
      CALL section_vals_val_get(ot_section, "EPS_IRAC_QUICK_EXIT", r_val=settings%eps_irac_quick_exit)
      IF (settings%eps_irac_quick_exit < 0.0_dp) THEN
         CPABORT("READ OT EPS_IRAC_QUICK_EXIT: VALUE MUST BE GREATER THAN ZERO")
      END IF

      CALL section_vals_val_get(ot_section, "EPS_IRAC_SWITCH", r_val=settings%eps_irac_switch)
      IF (settings%eps_irac_switch < 0.0_dp) THEN
         CPABORT("READ OT EPS_IRAC_SWITCH: VALUE MUST BE GREATER THAN ZERO")
      END IF

      CALL section_vals_val_get(ot_section, "ORTHO_IRAC", i_val=ot_ortho_irac)
      SELECT CASE (ot_ortho_irac)
      CASE (ot_chol_irac)
         settings%ortho_irac = "CHOL"
      CASE (ot_poly_irac)
         settings%ortho_irac = "POLY"
      CASE (ot_lwdn_irac)
         settings%ortho_irac = "LWDN"
      CASE DEFAULT
         CPABORT("READ OT ORTHO_IRAC: Value unknown")
      END SELECT

      CALL section_vals_val_get(ot_section, "ON_THE_FLY_LOC", l_val=settings%on_the_fly_loc)

      CALL section_vals_val_get(ot_section, "MINIMIZER", i_val=ot_method)
      ! overwrite input if ot_state is set
      IF (settings%ot_state == 1) THEN
         ot_method = ot_mini_diis
      END IF
      ! compatibility
      SELECT CASE (ot_method)
      CASE (ot_mini_sd)
         settings%ot_method = "SD"
      CASE (ot_mini_cg)
         settings%ot_method = "CG"
      CASE (ot_mini_diis)
         settings%ot_method = "DIIS"
         CALL section_vals_val_get(ot_section, "N_HISTORY_VEC", i_val=settings%diis_m)
      CASE (ot_mini_broyden)
         CALL section_vals_val_get(ot_section, "N_HISTORY_VEC", i_val=settings%diis_m)
         CALL section_vals_val_get(ot_section, "BROYDEN_BETA", r_val=settings%broyden_beta)
         CALL section_vals_val_get(ot_section, "BROYDEN_GAMMA", r_val=settings%broyden_gamma)
         CALL section_vals_val_get(ot_section, "BROYDEN_SIGMA", r_val=settings%broyden_sigma)
         CALL section_vals_val_get(ot_section, "BROYDEN_ETA", r_val=settings%broyden_eta)
         CALL section_vals_val_get(ot_section, "BROYDEN_OMEGA", r_val=settings%broyden_omega)
         CALL section_vals_val_get(ot_section, "BROYDEN_SIGMA_DECREASE", r_val=settings%broyden_sigma_decrease)
         CALL section_vals_val_get(ot_section, "BROYDEN_SIGMA_MIN", r_val=settings%broyden_sigma_min)
         CALL section_vals_val_get(ot_section, "BROYDEN_FORGET_HISTORY", l_val=settings%broyden_forget_history)
         CALL section_vals_val_get(ot_section, "BROYDEN_ADAPTIVE_SIGMA", l_val=settings%broyden_adaptive_sigma)
         CALL section_vals_val_get(ot_section, "BROYDEN_ENABLE_FLIP", l_val=settings%broyden_enable_flip)
         settings%ot_method = "BROY"
      CASE DEFAULT
         CPABORT("READ OTSCF MINIMIZER: Value unknown")
      END SELECT
      CALL section_vals_val_get(ot_section, "SAFER_DIIS", l_val=settings%safer_diis)
      CALL section_vals_val_get(ot_section, "LINESEARCH", i_val=ls_method)
      SELECT CASE (ls_method)
      CASE (ls_none)
         settings%line_search_method = "NONE"
      CASE (ls_2pnt)
         settings%line_search_method = "2PNT"
      CASE (ls_3pnt)
         settings%line_search_method = "3PNT"
      CASE (ls_adapt)
         settings%line_search_method = "ADPT"
      CASE (ls_gold)
         settings%line_search_method = "GOLD"
         CALL section_vals_val_get(ot_section, "GOLD_TARGET", r_val=settings%gold_target)
      CASE DEFAULT
         CPABORT("READ OTSCF LS: Value unknown")
      END SELECT

      CALL section_vals_val_get(ot_section, "PRECOND_SOLVER", i_val=settings%precond_solver_type)
      SELECT CASE (settings%precond_solver_type)
      CASE (ot_precond_solver_default)
         settings%precond_solver_name = "DEFAULT"
      CASE (ot_precond_solver_inv_chol)
         settings%precond_solver_name = "INVERSE_CHOLESKY"
      CASE (ot_precond_solver_direct)
         settings%precond_solver_name = "DIRECT"
      CASE (ot_precond_solver_update)
         settings%precond_solver_name = "INVERSE_UPDATE"
      CASE (ot_precond_solver_chebyshev)
         settings%precond_solver_name = "CHEBYSHEV"
      CASE DEFAULT
         CPABORT("READ OTSCF SOLVER: Value unknown")
      END SELECT
      CALL section_vals_val_get(ot_section, "CHEBYSHEV_DEGREE", i_val=settings%chebyshev_degree)
      IF (settings%chebyshev_degree < 1 .OR. settings%chebyshev_degree > 100) THEN
         CPABORT("READ OT CHEBYSHEV_DEGREE: Value must be between 1 and 100")
      END IF

      CALL section_vals_val_get(ot_section, "MAX_SCF_DIIS", i_val=settings%max_scf_diis)

      !If these values are negative we will set them "optimal" for a given precondtioner below
      CALL section_vals_val_get(ot_section, "STEPSIZE", r_val=settings%ds_min)
      CALL section_vals_val_get(ot_section, "ENERGY_GAP", r_val=settings%energy_gap)

      CALL section_vals_val_get(ot_section, "PRECONDITIONER", i_val=settings%preconditioner_type, &
                                explicit=preconditioner_explicit)
      IF (PRESENT(smearing_default_preconditioner)) THEN
         IF (smearing_default_preconditioner .AND. .NOT. preconditioner_explicit) THEN
            settings%preconditioner_type = ot_precond_full_single_inverse
            IF (PRESENT(smearing_default_kpoint_full_all)) THEN
               IF (smearing_default_kpoint_full_all) settings%preconditioner_type = ot_precond_full_all
            END IF
            selected_preconditioner_default = .TRUE.
         END IF
      END IF
      SELECT CASE (settings%preconditioner_type)
      CASE (ot_precond_none)
         settings%preconditioner_name = "NONE"
         IF (settings%ds_min < 0.0_dp) settings%ds_min = 0.15_dp
         IF (settings%energy_gap < 0.0_dp) settings%energy_gap = 0.2_dp
      CASE (ot_precond_full_single)
         settings%preconditioner_name = "FULL_SINGLE"
         IF (settings%ds_min < 0.0_dp) settings%ds_min = 0.15_dp
         IF (settings%energy_gap < 0.0_dp) settings%energy_gap = 0.2_dp
      CASE (ot_precond_full_single_inverse)
         settings%preconditioner_name = "FULL_SINGLE_INVERSE"
         IF (settings%ds_min < 0.0_dp) settings%ds_min = 0.08_dp
         IF (settings%energy_gap < 0.0_dp) settings%energy_gap = 0.08_dp
      CASE (ot_precond_full_all)
         settings%preconditioner_name = "FULL_ALL"
         IF (settings%ds_min < 0.0_dp) THEN
            settings%ds_min = 0.15_dp
            IF (kpoint_smearing_context .AND. smearing_context) settings%ds_min = 0.08_dp
         END IF
         IF (settings%energy_gap < 0.0_dp) settings%energy_gap = 0.08_dp
      CASE (ot_precond_fermi_low_rank)
         settings%preconditioner_name = "FERMI_LOW_RANK"
         IF (settings%ds_min < 0.0_dp) settings%ds_min = 0.12_dp
         IF (settings%energy_gap < 0.0_dp) settings%energy_gap = 0.08_dp
      CASE (ot_precond_full_kinetic)
         settings%preconditioner_name = "FULL_KINETIC"
         IF (settings%ds_min < 0.0_dp) settings%ds_min = 0.15_dp
         IF (settings%energy_gap < 0.0_dp) settings%energy_gap = 0.2_dp
      CASE (ot_precond_s_inverse)
         settings%preconditioner_name = "FULL_S_INVERSE"
         IF (settings%ds_min < 0.0_dp) settings%ds_min = 0.15_dp
         IF (settings%energy_gap < 0.0_dp) settings%energy_gap = 0.2_dp
      CASE DEFAULT
         CPABORT("READ OTSCF PRECONDITIONER: Value unknown")
      END SELECT
      CALL section_vals_val_get(ot_section, "CHOLESKY", i_val=settings%cholesky_type)
      CALL section_vals_val_get(ot_section, "EPS_TAYLOR", r_val=settings%eps_taylor)
      CALL section_vals_val_get(ot_section, "MAX_TAYLOR", i_val=settings%max_taylor)
      CALL section_vals_val_get(ot_section, "ROTATION", l_val=settings%do_rotation, &
                                explicit=rotation_explicit)
      CALL section_vals_val_get(ot_section, "ENERGIES", l_val=settings%do_ener, &
                                explicit=energies_explicit)
      IF (PRESENT(smearing_default_metallic_variables)) THEN
         IF (smearing_default_metallic_variables) THEN
            CALL ot_metallic_variable_defaults( &
               settings%do_rotation, settings%do_ener, rotation_explicit, energies_explicit, &
               selected_rotation_default, selected_energies_default)
         END IF
      END IF
      CALL section_vals_val_get(ot_section, "OCCUPATION_PRECONDITIONER", &
                                l_val=settings%occupation_preconditioner, &
                                explicit=occupation_preconditioner_explicit)
      IF (PRESENT(smearing_default_occupation_preconditioner)) THEN
         IF (smearing_default_occupation_preconditioner .AND. &
             settings%do_rotation .AND. settings%do_ener .AND. &
             .NOT. occupation_preconditioner_explicit) THEN
            settings%occupation_preconditioner = .TRUE.
            selected_occupation_preconditioner_default = .TRUE.
         END IF
      END IF
      CALL section_vals_val_get(ot_section, "NONDIAG_ENERGY", l_val=settings%add_nondiag_energy)
      CALL section_vals_val_get(ot_section, "NONDIAG_ENERGY_STRENGTH", &
                                r_val=settings%nondiag_energy_strength)
      CALL section_vals_val_get(ot_section, "CANONICAL_TETHER", l_val=settings%canonical_tether)
      CALL section_vals_val_get(ot_section, "CANONICAL_TETHER_ADAPTIVE", &
                                l_val=settings%canonical_tether_adaptive)
      CALL section_vals_val_get(ot_section, "CANONICAL_TETHER_STRENGTH", &
                                r_val=settings%canonical_tether_strength)
      CALL section_vals_val_get(ot_section, "CANONICAL_TETHER_MAX_STRENGTH", &
                                r_val=settings%canonical_tether_max_strength)
      IF (settings%do_ener .AND. .NOT. settings%do_rotation) THEN
         CPABORT("OT%ENERGIES requires OT%ROTATION.")
      END IF
      IF (settings%occupation_preconditioner .AND. .NOT. settings%do_ener) THEN
         CPABORT("OT%OCCUPATION_PRECONDITIONER requires OT%ENERGIES.")
      END IF
      IF (settings%add_nondiag_energy .AND. .NOT. settings%do_ener) THEN
         CPABORT("OT%NONDIAG_ENERGY requires OT%ENERGIES.")
      END IF
      IF (settings%canonical_tether .AND. .NOT. settings%occupation_preconditioner) THEN
         settings%canonical_tether = .FALSE.
      END IF
      IF (settings%canonical_tether_strength < 0.0_dp .OR. &
          settings%canonical_tether_max_strength < 0.0_dp) THEN
         CPABORT("OT%CANONICAL_TETHER strengths must be non-negative.")
      END IF
      IF (settings%canonical_tether_max_strength < settings%canonical_tether_strength) THEN
         settings%canonical_tether_max_strength = settings%canonical_tether_strength
      END IF

      ! write OT output

      IF (output_unit > 0) THEN
         WRITE (output_unit, '(/,A)') "  ----------------------------------- OT ---------------------------------------"
         IF (selected_kpoint_default) THEN
            WRITE (output_unit, '(A)') "  K-point OT default: ALGORITHM IRAC"
         END IF
         IF (selected_rotation_default) THEN
            IF (kpoint_smearing_context) THEN
               WRITE (output_unit, '(A)') "  K-point OT smearing default: ROTATION ON"
            ELSE
               WRITE (output_unit, '(A)') "  OT smearing default: ROTATION ON"
            END IF
         END IF
         IF (selected_energies_default) THEN
            IF (kpoint_smearing_context) THEN
               WRITE (output_unit, '(A)') "  K-point OT smearing default: ENERGIES ON"
            ELSE
               WRITE (output_unit, '(A)') "  OT smearing default: ENERGIES ON"
            END IF
         END IF
         IF (selected_occupation_preconditioner_default) THEN
            IF (kpoint_smearing_context) THEN
               WRITE (output_unit, '(A)') &
                  "  K-point OT smearing default: OCCUPATION_PRECONDITIONER ON"
            ELSE
               WRITE (output_unit, '(A)') &
                  "  OT smearing default: OCCUPATION_PRECONDITIONER ON"
            END IF
         END IF
         IF (selected_preconditioner_default) THEN
            IF (kpoint_smearing_context) THEN
               WRITE (output_unit, '(A)') &
                  "  K-point OT smearing default: PRECONDITIONER "//TRIM(settings%preconditioner_name)
            ELSE
               WRITE (output_unit, '(A)') &
                  "  OT smearing default: PRECONDITIONER "//TRIM(settings%preconditioner_name)
            END IF
         END IF
         IF (settings%do_rotation) THEN
            WRITE (output_unit, '(A)') "  Allowing for rotations "
         END IF
         IF (settings%do_ener) THEN
            WRITE (output_unit, '(A)') "  Optimizing orbital energies "
         END IF
         IF (settings%occupation_preconditioner) THEN
            WRITE (output_unit, '(A)') "  Using occupation preconditioner "
         END IF
         IF (settings%add_nondiag_energy) THEN
            WRITE (output_unit, '(A,E14.5)') &
               "  Non-diagonal energy strength:", settings%nondiag_energy_strength
         END IF
         IF (settings%canonical_tether) THEN
            WRITE (output_unit, '(A,E14.5,A,L1,A,E14.5)') &
               "  Canonical energy tether strength:", settings%canonical_tether_strength, &
               " adaptive:", settings%canonical_tether_adaptive, &
               " max:", settings%canonical_tether_max_strength
         END IF
         SELECT CASE (settings%OT_METHOD)
         CASE ("SD")
            WRITE (output_unit, '(A)') "  Minimizer      : SD                  : steepest descent"
         CASE ("CG")
            WRITE (output_unit, '(A)') "  Minimizer      : CG                  : conjugate gradient"
         CASE ("DIIS")
            WRITE (output_unit, '(A)') "  Minimizer      : DIIS                : direct inversion"
            WRITE (output_unit, '(A)') "                                         in the iterative subspace"
            WRITE (output_unit, '(A,I3,A)') "                                         using ", settings%diis_m, " DIIS vectors"
            IF (settings%safer_diis) THEN
               WRITE (output_unit, '(A,I3,A)') "                                         safer DIIS on"
            ELSE
               WRITE (output_unit, '(A,I3,A)') "                                         safer DIIS off"
            END IF
         CASE ("BROY")
            WRITE (output_unit, '(A)') "  Minimizer      : BROYDEN             : Broyden "
            WRITE (output_unit, '(A,F16.8)') "                   BETA                : ", settings%broyden_beta
            WRITE (output_unit, '(A,F16.8)') "                   GAMMA               : ", settings%broyden_gamma
            WRITE (output_unit, '(A,F16.8)') "                   SIGMA               : ", settings%broyden_sigma
            WRITE (output_unit, '(A,I3,A)') "                            using      : - ", &
               settings%diis_m, " BROYDEN vectors"
         CASE DEFAULT
            WRITE (output_unit, '(3A)') "  Minimizer      :      ", settings%OT_METHOD, " : UNKNOWN"
         END SELECT
         SELECT CASE (settings%preconditioner_name)
         CASE ("FULL_SINGLE")
            WRITE (output_unit, '(A)') "  Preconditioner : FULL_SINGLE         : diagonalization based"
         CASE ("FULL_SINGLE_INVERSE")
            WRITE (output_unit, '(A,/,A)') "  Preconditioner : FULL_SINGLE_INVERSE : inversion of ", &
               "                                         H + eS - 2*(Sc)(c^T*H*c+const)(Sc)^T"
         CASE ("FULL_ALL")
            WRITE (output_unit, '(A)') "  Preconditioner : FULL_ALL            : diagonalization, state selective"
         CASE ("FERMI_LOW_RANK")
            WRITE (output_unit, '(A)') &
               "  Preconditioner : FERMI_LOW_RANK       : bounded Fermi-window spectral correction"
         CASE ("FULL_KINETIC")
            WRITE (output_unit, '(A)') "  Preconditioner : FULL_KINETIC        : inversion of T + eS"
         CASE ("FULL_S_INVERSE")
            WRITE (output_unit, '(A)') "  Preconditioner : FULL_S_INVERSE      : cholesky inversion of S"
         CASE ("SPARSE_DIAG")
            WRITE (output_unit, '(A)') &
               "  Preconditioner : SPARSE_DIAG    : diagonal atomic block diagonalization"
         CASE ("SPARSE_KINETIC")
            WRITE (output_unit, '(A)') "  Preconditioner : SPARSE_KINETIC      : sparse linear solver for T + eS"
         CASE ("NONE")
            WRITE (output_unit, '(A)') "  Preconditioner : NONE"
         CASE DEFAULT
            WRITE (output_unit, '(3A)') "  Preconditioner : ", settings%preconditioner_name, " : UNKNOWN"
         END SELECT

         WRITE (output_unit, '(A)') "  Precond_solver : "//TRIM(settings%precond_solver_name)
         IF (settings%precond_solver_type == ot_precond_solver_chebyshev) THEN
            WRITE (output_unit, '(A,I14)') "  Chebyshev degree:", settings%chebyshev_degree
         END IF

         IF (settings%OT_METHOD == "SD" .OR. settings%OT_METHOD == "CG") THEN
            SELECT CASE (settings%line_search_method)
            CASE ("2PNT")
               WRITE (output_unit, '(A)') "  Line search    : 2PNT                : 2 energies, one gradient"
            CASE ("3PNT")
               WRITE (output_unit, '(A)') "  Line search    : 3PNT                : 3 energies"
            CASE ("GOLD")
               WRITE (output_unit, '(A)') "  Line search    : GOLD                : bracketing and golden section search"
               WRITE (output_unit, '(A,F14.8)') "                   target rel accuracy : ", settings%gold_target
            CASE ("NONE")
               WRITE (output_unit, '(A)') "  Line search    : NONE"
            CASE DEFAULT
               WRITE (output_unit, '(3A)') "  Line search : ", settings%line_search_method, " : UNKNOWN"
            END SELECT
         END IF
         WRITE (output_unit, '(A,F14.8,T49,A,F14.8)') "  stepsize       :", settings%ds_min, &
            "  energy_gap     :", settings%energy_gap
         IF (settings%ot_algorithm == 'TOD') THEN
            WRITE (output_unit, '(A,E14.5,T49,A,I14)') "  eps_taylor     :", settings%eps_taylor, &
               "  max_taylor     :", settings%max_taylor
         END IF
         IF (settings%ot_algorithm == 'REF') THEN
            WRITE (output_unit, '(A,1X,A,T49,A,I14)') "  ortho_irac     :", settings%ortho_irac, &
               "  irac_degree    :", settings%irac_degree
            WRITE (output_unit, '(A,I14,T49,A,E14.5)') "  max_irac       :", settings%max_irac, &
               "  eps_irac       :", settings%eps_irac
            WRITE (output_unit, '(A,E14.5,T49,A,E10.3)') "  eps_irac_switch:", settings%eps_irac_switch, &
               "  eps_irac_quick_exit:", settings%eps_irac_quick_exit
            WRITE (output_unit, '(A,L2)') "  on_the_fly_loc :", settings%on_the_fly_loc
         END IF
         WRITE (output_unit, '(A)') "  ----------------------------------- OT ---------------------------------------"
         WRITE (UNIT=output_unit, &
                FMT="(/,T3,A,T12,A,T31,A,T39,A,T59,A,T75,A,/,T3,A)") &
            "Step", "Update method", "Time", "Convergence", "Total energy", "Change", &
            REPEAT("-", 78)
      END IF

      CALL timestop(handle)

   END SUBROUTINE ot_readwrite_input

END MODULE qs_ot_types
