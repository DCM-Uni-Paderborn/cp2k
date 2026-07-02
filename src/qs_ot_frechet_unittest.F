!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief finite-difference checks for the complex OT rotation Frechet kernel
! **************************************************************************************************
PROGRAM qs_ot_frechet_unittest
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_create, dbcsr_distribution_new, dbcsr_distribution_release, dbcsr_distribution_type, &
        dbcsr_finalize, dbcsr_finalize_lib, dbcsr_get_block_p, dbcsr_init_lib, dbcsr_put_block, &
        dbcsr_release, dbcsr_reserve_blocks, dbcsr_type, dbcsr_type_no_symmetry
   USE kinds,                           ONLY: dp
   USE machine,                         ONLY: default_output_unit
   USE mathlib,                         ONLY: diag_complex,&
                                              diamat_all
   USE message_passing,                 ONLY: mp_comm_type,&
                                              mp_world_finalize,&
                                              mp_world_init
   USE qs_ot,                           ONLY: qs_ot_apply_complex_frechet_dbcsr,&
                                              qs_ot_complex_exp_frechet_kernel

   IMPLICIT NONE

   INTEGER                                            :: io_unit, mynode, nfail
   REAL(KIND=dp), PARAMETER                           :: projector_response_negative_balance_floor = 1.5E-1_dp, &
                                                         projector_response_negative_balance_max = 5.0E-1_dp
   TYPE(mp_comm_type)                                 :: mp_comm

   nfail = 0
   CALL mp_world_init(mp_comm)
   mynode = mp_comm%mepos
   io_unit = -1
   IF (mynode == 0) io_unit = default_output_unit
   CALL dbcsr_init_lib(mp_comm%get_handle(), io_unit)

   CALL check_offdiag("separated offdiag", 1.3_dp, -0.4_dp, &
                      CMPLX(0.7_dp, -0.2_dp, KIND=dp), 5.0E-12_dp, nfail)
   CALL check_offdiag("near-degenerate offdiag", 0.2_dp, 0.2000003_dp, &
                      CMPLX(-0.3_dp, 0.5_dp, KIND=dp), 5.0E-12_dp, nfail)
   CALL check_offdiag("degenerate offdiag", -0.6_dp, -0.6_dp, &
                      CMPLX(0.4_dp, 0.1_dp, KIND=dp), 5.0E-12_dp, nfail)
   CALL check_diag("diagonal", 0.8_dp, CMPLX(-0.2_dp, 0.35_dp, KIND=dp), &
                   5.0E-8_dp, nfail)
   CALL check_sparse_dbcsr_pattern("sparse DBCSR-style block pattern", 5.0E-14_dp, nfail)
   CALL check_dbcsr_sparse_pattern("DBCSR sparse block pattern", 5.0E-14_dp, nfail, mp_comm)
   CALL check_weighted_kpoint_mermin_gradient("weighted k-point Mermin gradient", 1.0E-7_dp, nfail)
   CALL check_kpoint_mermin_energy_projection("fixed-N Mermin energy projection", 5.0E-11_dp, nfail)
   CALL check_fixed_n_preconditioner_gauge("fixed-N preconditioner gauge", 5.0E-15_dp, nfail)
   CALL check_fixed_n_tether_preconditioner_gauge("fixed-N/tether preconditioner gauge", &
                                                  5.0E-15_dp, nfail)
   CALL check_weighted_kpoint_tether_preconditioner_gauge("weighted k-point/tether preconditioner gauge", &
                                                          5.0E-15_dp, nfail)
   CALL check_complex_kpoint_mermin_coupled_gradient("complex k-point Mermin coupled gradient", &
                                                     2.0E-7_dp, nfail)
   CALL check_complex_fixed_n_schur_hessian("complex fixed-N Schur Hessian", &
                                            2.0E-5_dp, nfail)
   CALL check_complex_fixed_n_response_channel("complex fixed-N response channel", &
                                               2.0E-7_dp, nfail)
   CALL check_complex_projector_blended_response_channel("complex projector-blended response channel", &
                                                         5.0E-9_dp, nfail)
   CALL check_complex_projector_physical_response_channel("complex projector physical response channel", &
                                                          5.0E-9_dp, nfail)
   CALL check_complex_projector_physical_hessian_split("complex projector physical Hessian split", &
                                                       2.0E-6_dp, nfail)
   CALL check_mixed_fixed_n_tether_projector("mixed fixed-N/tether Schur projector", &
                                             5.0E-14_dp, nfail)
   CALL check_mixed_fixed_n_tether_cross_projector("mixed fixed-N/tether Schur cross projector", &
                                                   5.0E-14_dp, nfail)
   CALL check_mixed_fixed_n_asymmetric_response_schur("mixed fixed-N asymmetric response Schur", &
                                                      5.0E-14_dp, nfail)
   CALL check_mixed_fixed_n_tether_exact_response("mixed fixed-N/tether exact response", &
                                                  5.0E-12_dp, nfail)
   CALL check_mixed_fixed_n_tether_exact_projection("mixed fixed-N/tether exact projection", &
                                                    5.0E-12_dp, nfail)
   CALL check_mixed_fixed_n_tether_exact_schur("mixed fixed-N/tether exact Schur", &
                                               5.0E-11_dp, nfail)
   CALL check_indefinite_schur_spectral_fallback("indefinite Schur spectral fallback", &
                                                 5.0E-12_dp, nfail)
   CALL check_indefinite_schur_absolute_metric("indefinite Schur absolute metric", &
                                               5.0E-12_dp, nfail)
   CALL check_complex_projector_rotation_gradient("complex projector rotation gradient", &
                                                  5.0E-8_dp, nfail)
   CALL check_complex_occ_projector_hessian("complex occupation-projector Hessian", &
                                            2.0E-6_dp, nfail)
   CALL check_complex_projector_covariant_curvature("complex projector covariant curvature", &
                                                    5.0E-7_dp, nfail)
   CALL check_complex_projector_covariant_gauge_split("complex projector covariant gauge split", &
                                                      5.0E-7_dp, nfail)
   CALL check_complex_projector_density_hessian_split("complex projector density Hessian split", &
                                                      5.0E-7_dp, nfail)
   CALL check_complex_projector_physical_invariant_base( &
      "complex projector physical invariant base", 5.0E-7_dp, nfail)
   CALL check_complex_occ_projector_gauge_separation("complex occupation-projector gauge separation", &
                                                     1.0E-7_dp, nfail)
   CALL check_complex_occ_projector_curvature_cap("complex occupation-projector curvature cap", &
                                                  5.0E-10_dp, nfail)
   CALL check_occ_projector_active_inert_pair_cap( &
      "complex occupation-projector active/inert pair cap", 5.0E-10_dp, nfail)
   CALL check_occ_projector_far_active_inert_pair_max_cap( &
      "complex occupation-projector far active/inert pair max cap", 5.0E-10_dp, nfail)
   CALL check_complex_occ_projector_gauge_plateau("complex occupation-projector gauge plateau", &
                                                  5.0E-10_dp, nfail)
   CALL check_occ_projector_plateau_caps_far_slope( &
      "complex occupation-projector plateau caps non-Fermi slope", 5.0E-10_dp, nfail)
   CALL check_complex_occ_projector_gauge_covariance("complex occupation-projector gauge covariance", &
                                                     5.0E-10_dp, nfail)
   CALL check_complex_projector_density_response_projection("complex projector density-response projection", &
                                                            2.0E-8_dp, nfail)
   CALL check_complex_projector_density_gauge_plateau_projection( &
      "complex projector density gauge-plateau projection", 5.0E-10_dp, nfail)
   CALL check_complex_projector_active_density_response_filter( &
      "complex projector active density-response filter", 5.0E-8_dp, nfail)
   CALL check_complex_projector_active_density_safe_inverse( &
      "complex projector active density safe inverse", 5.0E-14_dp, nfail)
   CALL check_complex_projector_bounded_density_mix( &
      "complex projector bounded density-response mix", 5.0E-14_dp, nfail)
   CALL check_projector_bounded_density_response_helper( &
      "complex projector bounded density-response helper", 5.0E-14_dp, nfail)
   CALL check_projector_physical_density_response_helper( &
      "complex projector physical density-response helper", 5.0E-14_dp, nfail)
   CALL check_projector_residual_gauge_density_response( &
      "complex projector residual gauge density-response", 5.0E-14_dp, nfail)
   CALL check_projector_physical_density_centered_sign( &
      "complex projector physical density centered-sign blend", 5.0E-14_dp, nfail)
   CALL check_projector_physical_density_centered_magnitude( &
      "complex projector physical density centered-magnitude blend", 5.0E-14_dp, nfail)
   CALL check_projector_physical_density_global_alignment( &
      "complex projector physical density global alignment", 5.0E-14_dp, nfail)
   CALL check_projector_physical_density_global_low_curv( &
      "complex projector physical density global low-curv response", 5.0E-14_dp, nfail)
   CALL check_projector_physical_density_response_convexity( &
      "complex projector physical density-response convexity", 5.0E-14_dp, nfail)
   CALL check_projector_density_floor_audit_response_helper( &
      "complex projector density floor-audit response helper", 5.0E-14_dp, nfail)
   CALL check_projector_physical_density_safe_source( &
      "complex projector physical density safe source", 5.0E-14_dp, nfail)
   CALL check_projector_active_density_total_invariant( &
      "complex projector active density total invariant", 5.0E-8_dp, nfail)
   CALL check_complex_projector_density_raw_audit_split( &
      "complex projector density raw-audit split", 5.0E-8_dp, nfail)
   CALL check_complex_projector_density_eigen_response_split( &
      "complex projector density eigen-response split", 2.0E-8_dp, nfail)
   CALL check_complex_projector_density_cluster_split_response( &
      "complex projector density cluster-split response", 5.0E-12_dp, nfail)
   CALL check_complex_projector_density_near_cluster_response( &
      "complex projector density near-cluster response", 5.0E-11_dp, nfail)
   CALL check_projector_density_local_mean_dP_rejection( &
      "complex projector density local-mean dP rejection", 5.0E-14_dp, nfail)
   CALL check_projector_density_plateau_trust_centering( &
      "complex projector density plateau-trust centering", 5.0E-14_dp, nfail)
   CALL check_complex_projector_density_eigen_floor_filter( &
      "complex projector density eigen floor filter", 5.0E-14_dp, nfail)
   CALL check_projector_active_rotation_pair_weight( &
      "complex projector active rotation pair weight", 5.0E-14_dp, nfail)
   CALL check_projector_physical_pair_projector_weight( &
      "complex projector physical pair weight", 5.0E-14_dp, nfail)
   CALL check_projector_density_descent_gap_projection( &
      "complex projector density descent-gap projection", 5.0E-14_dp, nfail)
   CALL check_projector_density_active_filter_descent_recovery( &
      "complex projector density active-filter descent recovery", 5.0E-14_dp, nfail)
   CALL check_projector_density_descent_history_optimum_cap( &
      "complex projector density-descent history optimum cap", 5.0E-14_dp, nfail)
   CALL check_projector_density_metric_descent_recovery( &
      "complex projector density metric-descent recovery", 5.0E-13_dp, nfail)
   CALL check_projector_density_mixed_metric_descent_recovery( &
      "complex projector density mixed-metric descent recovery", 5.0E-13_dp, nfail)
   CALL check_projector_density_energy_descent_fallback( &
      "complex projector density energy-descent fallback", 5.0E-13_dp, nfail)
   CALL check_projector_density_mixed_metric_coupled_schur( &
      "complex projector density mixed-metric coupled Schur", 5.0E-13_dp, nfail)
   CALL check_projector_response_mixed_metric_scale( &
      "complex projector response mixed-metric scale", 5.0E-13_dp, nfail)
   CALL check_projector_response_projected_linear_term( &
      "complex projector response projected linear term", 5.0E-13_dp, nfail)
   CALL check_projector_density_candidate_metric_model( &
      "complex projector density candidate metric model", 5.0E-13_dp, nfail)
   CALL check_projector_density_metric_coupled_schur( &
      "complex projector density metric-coupled Schur", 5.0E-13_dp, nfail)
   CALL check_complex_projector_offdiag_gauge_plateau_hessian( &
      "complex projector offdiag gauge-plateau Hessian", 5.0E-10_dp, nfail)
   CALL check_complex_projector_density_curvature_addback("complex projector density curvature addback", &
                                                          5.0E-12_dp, nfail)
   CALL check_complex_projector_density_coupled_schur("complex projector density coupled Schur", &
                                                      5.0E-12_dp, nfail)
   CALL check_projector_physical_direct_density_shadow( &
      "complex projector physical-direct density shadow", 5.0E-12_dp, nfail)
   CALL check_projector_density_candidate_total_response( &
      "complex projector density candidate total response", 5.0E-12_dp, nfail)
   CALL check_projector_density_direct_history_mismatch( &
      "complex projector density direct history mismatch", 5.0E-14_dp, nfail)
   CALL check_projector_density_descent_selected_rotation( &
      "complex projector density-descent selected rotation", 5.0E-14_dp, nfail)
   CALL check_projector_metric_response_projection("complex projector metric response projection", &
                                                   5.0E-13_dp, nfail)
   CALL check_projector_density_global_fixed_n_projection( &
      "complex projector density global fixed-N projection", 5.0E-14_dp, nfail)
   CALL check_projector_density_exact_source_projection( &
      "complex projector density exact source projection", 5.0E-12_dp, nfail)
   CALL check_projector_density_floor_global_projection( &
      "complex projector density floor global projection", 5.0E-14_dp, nfail)
   CALL check_projector_density_direct_floor_response( &
      "complex projector density direct floor response", 5.0E-14_dp, nfail)
   CALL check_projector_density_floor_bootstrap_weight( &
      "complex projector density floor bootstrap weight", 5.0E-14_dp, nfail)
   CALL check_projector_density_global_metric_trust( &
      "complex projector density global metric trust", 5.0E-14_dp, nfail)
   CALL check_projector_density_audit_trust_quality( &
      "complex projector density audit trust quality", 5.0E-14_dp, nfail)
   CALL check_projector_direct_trace_compensation("complex projector direct trace compensation", &
                                                  5.0E-14_dp, nfail)
   CALL check_projector_trace_split_history_quality("complex projector trace-split Direct quality", &
                                                    5.0E-14_dp, nfail)
   CALL check_projector_trace_invariant_history_quality( &
      "complex projector invariant-trace Direct trust", 5.0E-14_dp, nfail)
   CALL check_projector_direct_selected_history_quality( &
      "complex projector Direct selected-history quality", 5.0E-14_dp, nfail)
   CALL check_projector_inactive_offdiag_selected_tail_quality( &
      "complex projector inactive-offdiag selected tail quality", 5.0E-14_dp, nfail)
   CALL check_projector_offdiag_tail_residual_quality( &
      "complex projector offdiag tail residual quality", 5.0E-14_dp, nfail)
   CALL check_projector_response_underprediction_quality( &
      "complex projector response underprediction quality", 5.0E-14_dp, nfail)
   CALL check_projector_response_adverse_fraction( &
      "complex projector response adverse negative fraction", 5.0E-14_dp, nfail)
   CALL check_projector_response_residual_weight( &
      "complex projector response raw/effective residual weight", 5.0E-14_dp, nfail)
   CALL check_projector_direct_tail_selection_quality( &
      "complex projector Direct tail selection quality", 5.0E-14_dp, nfail)
   CALL check_projector_direct_tail_probe_bootstrap( &
      "complex projector Direct tail-probe bootstrap", 5.0E-14_dp, nfail)
   CALL check_projector_density_full_direct_tail_probe( &
      "complex projector density full-Direct tail probe", 5.0E-14_dp, nfail)
   CALL check_projector_density_descent_tail_probe( &
      "complex projector density-descent tail probe", 5.0E-14_dp, nfail)
   CALL check_projector_density_full_direct_probe_drop( &
      "complex projector density full-Direct probe drop", 5.0E-14_dp, nfail)
   CALL check_projector_density_selected_tail_mix( &
      "complex projector density selected tail mix", 5.0E-14_dp, nfail)
   CALL check_projector_dense_metric_tail_bridge( &
      "complex projector dense metric tail bridge", 5.0E-14_dp, nfail)
   CALL check_projector_residual_gauge_floor( &
      "complex projector residual gauge floor", 5.0E-14_dp, nfail)
   CALL check_projector_trace_active_mix_weight("complex projector trace-active Direct mix", &
                                                5.0E-14_dp, nfail)
   CALL check_projector_trace_active_basis_weight("complex projector trace-active Direct basis", &
                                                  5.0E-14_dp, nfail)
   CALL check_projector_trace_diag_component_scale("complex projector trace-diag Direct scale", &
                                                   5.0E-14_dp, nfail)
   CALL check_projector_trace_local_compensation("complex projector trace local compensation", &
                                                 5.0E-14_dp, nfail)
   CALL check_projector_direct_history_step_scale("complex projector Direct history step scale", &
                                                  5.0E-14_dp, nfail)
   CALL check_projector_direct_scaled_metric_norm("complex projector Direct scaled metric norm", &
                                                  5.0E-14_dp, nfail)
   CALL check_projector_spectral_metric_norm_scale("complex projector spectral metric norm scale", &
                                                   5.0E-14_dp, nfail)
   CALL check_projector_spectral_selected_history_promotion( &
      "complex projector spectral selected-history promotion", 5.0E-14_dp, nfail)
   CALL check_projector_direct_history_consistency_scale( &
      "complex projector Direct history consistency scale", 5.0E-14_dp, nfail)
   CALL check_projector_mixed_direct_promotion("complex projector mixed Direct promotion", &
                                               5.0E-14_dp, nfail)
   CALL check_projector_coupled_tail_weight("complex projector coupled tail weight", &
                                            5.0E-14_dp, nfail)
   CALL check_projector_tail_active_delta_weight("complex projector active-delta tail weight", &
                                                 5.0E-14_dp, nfail)
   CALL check_projector_direct_delta_eigen_copy("complex projector Direct delta eigen copy", &
                                                5.0E-13_dp, nfail)
   CALL check_projector_active_residual_filter("complex projector active residual filter", &
                                               5.0E-14_dp, nfail)
   CALL check_projector_active_direct_delta_shadow("complex projector active Direct delta shadow", &
                                                   5.0E-14_dp, nfail)
   CALL check_projector_direct_trace_model_gap("complex projector direct trace model gap", &
                                               5.0E-14_dp, nfail)
   CALL check_projector_mixed_direct_response_blend("complex projector mixed Direct response blend", &
                                                    5.0E-14_dp, nfail)
   CALL check_projector_direct_floor_density_schur( &
      "complex projector Direct floor-density Schur blend", 5.0E-13_dp, nfail)
   CALL check_projector_direct_density_curvature_mix( &
      "complex projector Direct density curvature mix", 5.0E-14_dp, nfail)
   CALL check_projector_modal_response_sign_model("complex projector modal response sign model", &
                                                  5.0E-12_dp, nfail)
   CALL check_projector_modal_response_negative_floor( &
      "complex projector modal response negative-only floor", 5.0E-12_dp, nfail)
   CALL check_projector_negative_step_balance("complex projector negative-step balance", &
                                              5.0E-14_dp, nfail)
   CALL check_projector_tail_gauge_quality("complex projector tail gauge quality", &
                                           5.0E-14_dp, nfail)
   CALL check_projector_offdiag_tail_abs_trust( &
      "complex projector offdiag tail abs trust", 5.0E-14_dp, nfail)
   CALL check_projector_raw_tail_rotation_floor( &
      "complex projector raw-tail rotation floor", 5.0E-14_dp, nfail)
   CALL check_monotone_occupation_slope_projection("complex occupation monotone slope projection", &
                                                   5.0E-15_dp, nfail)
   CALL check_history_safe_occ_projector_derivative("history-safe occupation-projector derivative", &
                                                    5.0E-12_dp, nfail)
   CALL check_history_safe_occ_projector_gauge_physical_split( &
      "history-safe occupation-projector gauge/physical split", 5.0E-8_dp, nfail)
   CALL check_history_safe_occ_projector_physical_fd( &
      "history-safe occupation-projector physical finite difference", 5.0E-8_dp, nfail)
   CALL check_history_safe_occ_projector_slope_cap_fd( &
      "history-safe occupation-projector slope-cap finite difference", 5.0E-8_dp, nfail)
   CALL check_history_safe_occ_projector_adaptive_fd( &
      "history-safe occupation-projector adaptive finite difference", 5.0E-8_dp, nfail)
   CALL check_history_safe_occ_projector_diagonal_response( &
      "history-safe occupation-projector diagonal response", 5.0E-14_dp, nfail)
   CALL check_history_safe_occ_projector_fd_diagonal_response( &
      "history-safe occupation-projector fd diagonal response", 5.0E-8_dp, nfail)
   CALL check_history_safe_occ_projector_fd_mixed_response( &
      "history-safe occupation-projector fd mixed response", 5.0E-8_dp, nfail)
   CALL check_history_safe_occ_projector_fd_fixed_n_trace( &
      "history-safe occupation-projector fd fixed-N trace", 5.0E-8_dp, nfail)
   CALL check_history_safe_occ_projector_thermal_cluster_gap( &
      "history-safe occupation-projector thermal cluster gap", 5.0E-8_dp, nfail)
   CALL check_history_safe_occ_projector_thermal_value_window( &
      "history-safe occupation-projector thermal value window", 5.0E-8_dp, nfail)
   CALL check_history_safe_occ_projector_thermal_no_fd_offblock( &
      "history-safe occupation-projector thermal no-fd offblock", 5.0E-8_dp, nfail)
   CALL check_history_safe_occ_projector_diagonal_cluster_fd( &
      "history-safe occupation-projector diagonal cluster finite difference", 5.0E-8_dp, nfail)
   CALL check_history_safe_occ_projector_degenerate_physical_split( &
      "history-safe occupation-projector degenerate physical split", 5.0E-8_dp, nfail)
   CALL check_history_safe_occ_projector_modal_curvature_split( &
      "history-safe occupation-projector modal curvature split", 5.0E-8_dp, nfail)
   CALL check_history_safe_occ_projector_modal_mu_closure( &
      "history-safe occupation-projector modal mu closure", 5.0E-8_dp, nfail)
   CALL check_history_safe_occ_projector_modal_curvature_gauge_bound( &
      "history-safe occupation-projector modal curvature gauge bound", 5.0E-8_dp, nfail)
   CALL check_history_safe_occ_projector_modal_curvature_active_tail( &
      "history-safe occupation-projector modal curvature active tail", 5.0E-8_dp, nfail)
   CALL check_history_safe_occ_projector_cluster_slope_split( &
      "history-safe occupation-projector cluster slope split", 5.0E-8_dp, nfail)
   CALL check_history_safe_occ_projector_component_filter( &
      "history-safe occupation-projector component filter", 5.0E-8_dp, nfail)
   CALL check_history_safe_occ_projector_fixed_n_cluster_shape( &
      "history-safe occupation-projector fixed-N cluster shape", 5.0E-8_dp, nfail)
   CALL check_history_safe_occ_projector_covariant_fd( &
      "history-safe occupation-projector covariant finite difference", 5.0E-8_dp, nfail)
   CALL check_history_safe_occ_projector_cluster_offblock_response( &
      "history-safe occupation-projector cluster offblock response", 5.0E-12_dp, nfail)
   CALL check_history_safe_occ_projector_cluster_offblock_modal_cap( &
      "history-safe occupation-projector cluster offblock modal cap", 5.0E-12_dp, nfail)
   CALL check_history_safe_occ_projector_cluster_pair_response( &
      "history-safe occupation-projector cluster-pair response", 5.0E-12_dp, nfail)
   CALL check_history_safe_occ_projector_cluster_modal_occ_shape( &
      "history-safe occupation-projector cluster modal occupation", 5.0E-12_dp, nfail)
   CALL check_history_safe_occ_projector_cluster_kernel_invariants( &
      "history-safe occupation-projector cluster kernel invariants", 5.0E-12_dp, nfail)
   CALL check_complex_nondiag_tether_gradient("complex non-diagonal tether gradient", &
                                              5.0E-9_dp, nfail)
   CALL check_complex_canonical_tether_gradient("complex canonical tether gradient", &
                                                5.0E-9_dp, nfail)

   IF (nfail == 0) PRINT *, "qs_ot_frechet_unittest: OK"

   CALL dbcsr_finalize_lib()
   CALL mp_world_finalize()

   IF (nfail > 0) ERROR STOP "qs_ot_frechet_unittest failed"

CONTAINS

! **************************************************************************************************
!> \brief Check L_exp(diag(-i*e1,-i*e2), E_12*z)
!> \param label ...
!> \param e1 ...
!> \param e2 ...
!> \param z ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_offdiag(label, e1, e2, z, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: e1, e2
      COMPLEX(KIND=dp), INTENT(IN)                       :: z
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      REAL(KIND=dp), PARAMETER                           :: eps = 1.0E-7_dp

      COMPLEX(KIND=dp)                                   :: fd, ref
      COMPLEX(KIND=dp), DIMENSION(2, 2)                  :: m0, meps, u0, ueps
      REAL(KIND=dp)                                      :: err

      m0(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      m0(1, 1) = CMPLX(0.0_dp, -e1, KIND=dp)
      m0(2, 2) = CMPLX(0.0_dp, -e2, KIND=dp)

      meps = m0
      meps(1, 2) = eps*z

      CALL dense_expm2(m0, u0)
      CALL dense_expm2(meps, ueps)

      fd = (ueps(1, 2) - u0(1, 2))/eps
      ref = qs_ot_complex_exp_frechet_kernel(e1, e2)*z
      err = ABS(fd - ref)

      IF (err > tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,ES12.5)') TRIM(label), "failed with error", err
      ELSE
         WRITE (*, '(A,1X,A,1X,ES12.5)') TRIM(label), "OK with error", err
      END IF

   END SUBROUTINE check_offdiag

! **************************************************************************************************
!> \brief Check L_exp(diag(-i*e,-i*e), E_11*z)
!> \param label ...
!> \param e1 ...
!> \param z ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_diag(label, e1, z, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: e1
      COMPLEX(KIND=dp), INTENT(IN)                       :: z
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      REAL(KIND=dp), PARAMETER                           :: eps = 1.0E-7_dp

      COMPLEX(KIND=dp)                                   :: fd, ref
      COMPLEX(KIND=dp), DIMENSION(2, 2)                  :: m0, meps, u0, ueps
      REAL(KIND=dp)                                      :: err

      m0(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      m0(1, 1) = CMPLX(0.0_dp, -e1, KIND=dp)
      m0(2, 2) = CMPLX(0.0_dp, -e1, KIND=dp)

      meps = m0
      meps(1, 1) = meps(1, 1) + eps*z

      CALL dense_expm2(m0, u0)
      CALL dense_expm2(meps, ueps)

      fd = (ueps(1, 1) - u0(1, 1))/eps
      ref = qs_ot_complex_exp_frechet_kernel(e1, e1)*z
      err = ABS(fd - ref)

      IF (err > tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,ES12.5)') TRIM(label), "failed with error", err
      ELSE
         WRITE (*, '(A,1X,A,1X,ES12.5)') TRIM(label), "OK with error", err
      END IF

   END SUBROUTINE check_diag

! **************************************************************************************************
!> \brief Check finite-generator Frechet application with mismatched Re/Im block patterns.
!>        This mirrors the DBCSR iterator path where real-only, imaginary-only,
!>        mixed, and absent sparse blocks must be handled consistently.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_sparse_dbcsr_pattern(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 4, nb = 2

      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: dense_ref, sparse_out, z
      LOGICAL, DIMENSION(nb, nb)                         :: mask_im, mask_re
      REAL(KIND=dp)                                      :: err
      REAL(KIND=dp), DIMENSION(n)                        :: evals

      evals = [0.7_dp, -0.2_dp, 1.4_dp, -1.1_dp]
      z(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      z(1, 1) = CMPLX(0.3_dp, -0.4_dp, KIND=dp)
      z(1, 2) = CMPLX(-0.2_dp, 0.1_dp, KIND=dp)
      z(2, 1) = CMPLX(0.5_dp, 0.8_dp, KIND=dp)
      z(2, 2) = CMPLX(-0.7_dp, 0.2_dp, KIND=dp)
      z(1, 3) = CMPLX(0.0_dp, 0.6_dp, KIND=dp)
      z(1, 4) = CMPLX(0.0_dp, -0.3_dp, KIND=dp)
      z(2, 3) = CMPLX(0.0_dp, 0.4_dp, KIND=dp)
      z(2, 4) = CMPLX(0.0_dp, 0.2_dp, KIND=dp)
      z(3, 1) = CMPLX(-0.9_dp, 0.0_dp, KIND=dp)
      z(3, 2) = CMPLX(0.4_dp, 0.0_dp, KIND=dp)
      z(4, 1) = CMPLX(0.2_dp, 0.0_dp, KIND=dp)
      z(4, 2) = CMPLX(-0.1_dp, 0.0_dp, KIND=dp)

      mask_re(:, :) = .FALSE.
      mask_im(:, :) = .FALSE.
      mask_re(1, 1) = .TRUE.  ! real and imaginary block both present
      mask_im(1, 1) = .TRUE.
      mask_im(1, 2) = .TRUE.  ! imaginary-only block
      mask_re(2, 1) = .TRUE.  ! real-only block
      ! block (2,2) is absent in both patterns

      CALL apply_dense_frechet(evals, z, dense_ref)
      CALL apply_sparse_block_frechet(evals, z, mask_re, mask_im, sparse_out)
      err = MAXVAL(ABS(sparse_out - dense_ref))

      IF (err > tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,ES12.5)') TRIM(label), "failed with error", err
      ELSE
         WRITE (*, '(A,1X,A,1X,ES12.5)') TRIM(label), "OK with error", err
      END IF

   END SUBROUTINE check_sparse_dbcsr_pattern

! **************************************************************************************************
!> \brief Check the production DBCSR Frechet path with real/imaginary-only blocks.
!> \param label ...
!> \param tol ...
!> \param nfail ...
!> \param mp_comm ...
! **************************************************************************************************
   SUBROUTINE check_dbcsr_sparse_pattern(label, tol, nfail, mp_comm)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail
      TYPE(mp_comm_type), INTENT(IN)                     :: mp_comm

      INTEGER, PARAMETER                                 :: n = 4

      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: dbcsr_out, dense_ref, z
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_size, col_dist, row_blk_size, &
                                                            row_dist
      INTEGER, DIMENSION(:, :), POINTER                  :: pgrid
      LOGICAL                                            :: found
      REAL(KIND=dp)                                      :: err
      REAL(KIND=dp), DIMENSION(2, 2)                     :: block
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block_p
      REAL(KIND=dp), DIMENSION(n)                        :: evals
      TYPE(dbcsr_distribution_type)                      :: dist
      TYPE(dbcsr_type)                                   :: inner_im, inner_re, outer_im, outer_re

      evals = [0.7_dp, -0.2_dp, 1.4_dp, -1.1_dp]
      CALL setup_sparse_pattern_data(z)
      CALL apply_dense_frechet(evals, z, dense_ref)

      ALLOCATE (pgrid(0:0, 0:0), row_dist(2), col_dist(2), row_blk_size(2), col_blk_size(2))
      pgrid(:, :) = 0
      row_dist(:) = 0
      col_dist(:) = 0
      row_blk_size(:) = 2
      col_blk_size(:) = 2
      CALL dbcsr_distribution_new(dist, group=mp_comm%get_handle(), pgrid=pgrid, &
                                  row_dist=row_dist, col_dist=col_dist)

      CALL dbcsr_create(inner_re, "inner_re", dist, dbcsr_type_no_symmetry, row_blk_size, col_blk_size)
      CALL dbcsr_create(inner_im, "inner_im", dist, dbcsr_type_no_symmetry, row_blk_size, col_blk_size)
      CALL dbcsr_reserve_blocks(inner_re, rows=[1, 2], cols=[1, 1])
      CALL dbcsr_reserve_blocks(inner_im, rows=[1, 1], cols=[1, 2])
      CALL dbcsr_finalize(inner_re)
      CALL dbcsr_finalize(inner_im)

      block = REAL(z(1:2, 1:2), KIND=dp)
      CALL dbcsr_put_block(inner_re, 1, 1, block)
      block = AIMAG(z(1:2, 1:2))
      CALL dbcsr_put_block(inner_im, 1, 1, block)
      block = REAL(z(3:4, 1:2), KIND=dp)
      CALL dbcsr_put_block(inner_re, 2, 1, block)
      block = AIMAG(z(1:2, 3:4))
      CALL dbcsr_put_block(inner_im, 1, 2, block)

      CALL qs_ot_apply_complex_frechet_dbcsr(evals, inner_re, inner_im, outer_re, outer_im)

      dbcsr_out(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      CALL dbcsr_get_block_p(outer_re, 1, 1, block_p, found)
      IF (found) dbcsr_out(1:2, 1:2) = dbcsr_out(1:2, 1:2) + CMPLX(block_p, 0.0_dp, KIND=dp)
      CALL dbcsr_get_block_p(outer_im, 1, 1, block_p, found)
      IF (found) dbcsr_out(1:2, 1:2) = dbcsr_out(1:2, 1:2) + CMPLX(0.0_dp, block_p, KIND=dp)
      CALL dbcsr_get_block_p(outer_re, 1, 2, block_p, found)
      IF (found) dbcsr_out(1:2, 3:4) = dbcsr_out(1:2, 3:4) + CMPLX(block_p, 0.0_dp, KIND=dp)
      CALL dbcsr_get_block_p(outer_im, 1, 2, block_p, found)
      IF (found) dbcsr_out(1:2, 3:4) = dbcsr_out(1:2, 3:4) + CMPLX(0.0_dp, block_p, KIND=dp)
      CALL dbcsr_get_block_p(outer_re, 2, 1, block_p, found)
      IF (found) dbcsr_out(3:4, 1:2) = dbcsr_out(3:4, 1:2) + CMPLX(block_p, 0.0_dp, KIND=dp)
      CALL dbcsr_get_block_p(outer_im, 2, 1, block_p, found)
      IF (found) dbcsr_out(3:4, 1:2) = dbcsr_out(3:4, 1:2) + CMPLX(0.0_dp, block_p, KIND=dp)
      CALL dbcsr_get_block_p(outer_re, 2, 2, block_p, found)
      IF (found) dbcsr_out(3:4, 3:4) = dbcsr_out(3:4, 3:4) + CMPLX(block_p, 0.0_dp, KIND=dp)
      CALL dbcsr_get_block_p(outer_im, 2, 2, block_p, found)
      IF (found) dbcsr_out(3:4, 3:4) = dbcsr_out(3:4, 3:4) + CMPLX(0.0_dp, block_p, KIND=dp)

      err = MAXVAL(ABS(dbcsr_out - dense_ref))
      CALL dbcsr_release(outer_re)
      CALL dbcsr_release(outer_im)

      CALL qs_ot_apply_complex_frechet_dbcsr(evals, inner_re, inner_im, outer_re, outer_im, &
                                             adjoint=.TRUE.)
      CALL apply_dense_frechet(evals, z, dense_ref, adjoint=.TRUE.)

      dbcsr_out(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      CALL dbcsr_get_block_p(outer_re, 1, 1, block_p, found)
      IF (found) dbcsr_out(1:2, 1:2) = dbcsr_out(1:2, 1:2) + CMPLX(block_p, 0.0_dp, KIND=dp)
      CALL dbcsr_get_block_p(outer_im, 1, 1, block_p, found)
      IF (found) dbcsr_out(1:2, 1:2) = dbcsr_out(1:2, 1:2) + CMPLX(0.0_dp, block_p, KIND=dp)
      CALL dbcsr_get_block_p(outer_re, 1, 2, block_p, found)
      IF (found) dbcsr_out(1:2, 3:4) = dbcsr_out(1:2, 3:4) + CMPLX(block_p, 0.0_dp, KIND=dp)
      CALL dbcsr_get_block_p(outer_im, 1, 2, block_p, found)
      IF (found) dbcsr_out(1:2, 3:4) = dbcsr_out(1:2, 3:4) + CMPLX(0.0_dp, block_p, KIND=dp)
      CALL dbcsr_get_block_p(outer_re, 2, 1, block_p, found)
      IF (found) dbcsr_out(3:4, 1:2) = dbcsr_out(3:4, 1:2) + CMPLX(block_p, 0.0_dp, KIND=dp)
      CALL dbcsr_get_block_p(outer_im, 2, 1, block_p, found)
      IF (found) dbcsr_out(3:4, 1:2) = dbcsr_out(3:4, 1:2) + CMPLX(0.0_dp, block_p, KIND=dp)
      CALL dbcsr_get_block_p(outer_re, 2, 2, block_p, found)
      IF (found) dbcsr_out(3:4, 3:4) = dbcsr_out(3:4, 3:4) + CMPLX(block_p, 0.0_dp, KIND=dp)
      CALL dbcsr_get_block_p(outer_im, 2, 2, block_p, found)
      IF (found) dbcsr_out(3:4, 3:4) = dbcsr_out(3:4, 3:4) + CMPLX(0.0_dp, block_p, KIND=dp)

      err = MAX(err, MAXVAL(ABS(dbcsr_out - dense_ref)))
      IF (err > tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,ES12.5)') TRIM(label), "failed with error", err
      ELSE
         WRITE (*, '(A,1X,A,1X,ES12.5)') TRIM(label), "OK with error", err
      END IF

      CALL dbcsr_release(inner_re)
      CALL dbcsr_release(inner_im)
      CALL dbcsr_release(outer_re)
      CALL dbcsr_release(outer_im)
      CALL dbcsr_distribution_release(dist)
      DEALLOCATE (pgrid, row_dist, col_dist, row_blk_size, col_blk_size)

   END SUBROUTINE check_dbcsr_sparse_pattern

! **************************************************************************************************
!> \brief Fill the dense companion to the sparse Frechet pattern tests.
!> \param z ...
! **************************************************************************************************
   SUBROUTINE setup_sparse_pattern_data(z)
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(OUT)     :: z

      z(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      z(1, 1) = CMPLX(0.3_dp, -0.4_dp, KIND=dp)
      z(1, 2) = CMPLX(-0.2_dp, 0.1_dp, KIND=dp)
      z(2, 1) = CMPLX(0.5_dp, 0.8_dp, KIND=dp)
      z(2, 2) = CMPLX(-0.7_dp, 0.2_dp, KIND=dp)
      z(1, 3) = CMPLX(0.0_dp, 0.6_dp, KIND=dp)
      z(1, 4) = CMPLX(0.0_dp, -0.3_dp, KIND=dp)
      z(2, 3) = CMPLX(0.0_dp, 0.4_dp, KIND=dp)
      z(2, 4) = CMPLX(0.0_dp, 0.2_dp, KIND=dp)
      z(3, 1) = CMPLX(-0.9_dp, 0.0_dp, KIND=dp)
      z(3, 2) = CMPLX(0.4_dp, 0.0_dp, KIND=dp)
      z(4, 1) = CMPLX(0.2_dp, 0.0_dp, KIND=dp)
      z(4, 2) = CMPLX(-0.1_dp, 0.0_dp, KIND=dp)

   END SUBROUTINE setup_sparse_pattern_data

! **************************************************************************************************
!> \brief Dense reference application of the complex exponential Frechet kernel.
!> \param evals ...
!> \param z ...
!> \param out ...
!> \param adjoint ...
! **************************************************************************************************
   SUBROUTINE apply_dense_frechet(evals, z, out, adjoint)
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: evals
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: z
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(OUT)     :: out
      LOGICAL, INTENT(IN), OPTIONAL                      :: adjoint

      COMPLEX(KIND=dp)                                   :: kernel
      INTEGER                                            :: i, j
      LOGICAL                                            :: use_adjoint

      use_adjoint = .FALSE.
      IF (PRESENT(adjoint)) use_adjoint = adjoint

      DO i = 1, SIZE(z, 1)
      DO j = 1, SIZE(z, 2)
         kernel = qs_ot_complex_exp_frechet_kernel(evals(i), evals(j))
         IF (use_adjoint) kernel = CONJG(kernel)
         out(i, j) = z(i, j)*kernel
      END DO
      END DO

   END SUBROUTINE apply_dense_frechet

! **************************************************************************************************
!> \brief Sparse block-pattern application equivalent to the DBCSR block iterator logic.
!> \param evals ...
!> \param z ...
!> \param mask_re ...
!> \param mask_im ...
!> \param out ...
! **************************************************************************************************
   SUBROUTINE apply_sparse_block_frechet(evals, z, mask_re, mask_im, out)
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: evals
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: z
      LOGICAL, DIMENSION(:, :), INTENT(IN)               :: mask_re, mask_im
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(OUT)     :: out

      COMPLEX(KIND=dp)                                   :: cval
      INTEGER                                            :: bc, br, i, i1, i2, j, j1, j2

      out(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      DO br = 1, SIZE(mask_re, 1)
      DO bc = 1, SIZE(mask_re, 2)
         IF (.NOT. (mask_re(br, bc) .OR. mask_im(br, bc))) CYCLE
         i1 = (br - 1)*2 + 1
         i2 = br*2
         j1 = (bc - 1)*2 + 1
         j2 = bc*2
         DO i = i1, i2
         DO j = j1, j2
            cval = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
            IF (mask_re(br, bc)) cval = cval + CMPLX(REAL(z(i, j), KIND=dp), 0.0_dp, KIND=dp)
            IF (mask_im(br, bc)) cval = cval + CMPLX(0.0_dp, AIMAG(z(i, j)), KIND=dp)
            out(i, j) = cval*qs_ot_complex_exp_frechet_kernel(evals(i), evals(j))
         END DO
         END DO
      END DO
      END DO

   END SUBROUTINE apply_sparse_block_frechet

! **************************************************************************************************
!> \brief Check the weighted fixed-N Mermin gradient used for k-point OT energy variables.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_weighted_kpoint_mermin_gradient(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: nkp = 2, nmo = 4
      REAL(KIND=dp), PARAMETER                           :: fd_step = 1.0E-6_dp, maxocc = 2.0_dp, &
                                                            nel = 5.0_dp, sigma = 0.025_dp

      INTEGER                                            :: ik, imo
      REAL(KIND=dp)                                      :: err, fd, ref
      REAL(KIND=dp), DIMENSION(nkp)                      :: wk
      REAL(KIND=dp), DIMENSION(nmo, nkp)                 :: eps, eps_m, eps_p, grad, h

      wk = [0.25_dp, 0.75_dp]
      eps(:, 1) = [-0.20_dp, -0.060_dp, -0.005_dp, 0.080_dp]
      eps(:, 2) = [-0.16_dp, -0.040_dp, 0.015_dp, 0.110_dp]
      h(:, 1) = [-0.19_dp, -0.045_dp, 0.020_dp, 0.100_dp]
      h(:, 2) = [-0.15_dp, -0.055_dp, 0.035_dp, 0.090_dp]

      CALL dense_kpoint_mermin_gradient(eps, h, wk, nel, sigma, maxocc, grad)

      err = 0.0_dp
      DO ik = 1, nkp
         DO imo = 1, nmo
            eps_p(:, :) = eps(:, :)
            eps_m(:, :) = eps(:, :)
            eps_p(imo, ik) = eps_p(imo, ik) + fd_step
            eps_m(imo, ik) = eps_m(imo, ik) - fd_step
            fd = (dense_kpoint_mermin_value(eps_p, h, wk, nel, sigma, maxocc) - &
                  dense_kpoint_mermin_value(eps_m, h, wk, nel, sigma, maxocc))/(2.0_dp*fd_step)
            ref = grad(imo, ik)
            err = MAX(err, ABS(fd - ref))
         END DO
      END DO

      IF (err > tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,ES12.5)') TRIM(label), "failed with error", err
      ELSE
         WRITE (*, '(A,1X,A,1X,ES12.5)') TRIM(label), "OK with error", err
      END IF

   END SUBROUTINE check_weighted_kpoint_mermin_gradient

! **************************************************************************************************
!> \brief Check the exact auxiliary-energy projection for the fixed-N Mermin subproblem.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_kpoint_mermin_energy_projection(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: nkp = 2, nmo = 4
      REAL(KIND=dp), PARAMETER                           :: maxocc = 2.0_dp, nel = 4.35_dp, &
                                                            sigma = 0.030_dp

      REAL(KIND=dp)                                      :: err, gauge_shift
      REAL(KIND=dp), DIMENSION(nkp)                      :: wk
      REAL(KIND=dp), DIMENSION(nmo, nkp)                 :: eps, grad, h, projected

      wk = [0.35_dp, 0.65_dp]
      eps(:, 1) = [-0.34_dp, -0.07_dp, 0.025_dp, 0.18_dp]
      eps(:, 2) = [-0.29_dp, -0.02_dp, 0.075_dp, 0.22_dp]
      h(:, 1) = [-0.31_dp, -0.035_dp, 0.045_dp, 0.16_dp]
      h(:, 2) = [-0.27_dp, -0.055_dp, 0.055_dp, 0.20_dp]

      gauge_shift = SUM(eps(:, :) - h(:, :))/REAL(nkp*nmo, KIND=dp)
      projected(:, :) = h(:, :) + gauge_shift
      CALL dense_kpoint_mermin_gradient(projected, h, wk, nel, sigma, maxocc, grad)
      err = MAXVAL(ABS(grad))

      IF (err > tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,ES12.5)') TRIM(label), "failed with error", err
      ELSE
         WRITE (*, '(A,1X,A,1X,ES12.5)') TRIM(label), "OK with error", err
      END IF

   END SUBROUTINE check_kpoint_mermin_energy_projection

! **************************************************************************************************
!> \brief Check that the diagonal energy preconditioner removes the fixed-N gauge
!>        in the occupation-curvature metric, not by an arithmetic band average.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_fixed_n_preconditioner_gauge(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 5
      REAL(KIND=dp), PARAMETER                           :: curvature_floor = 1.0_dp

      REAL(KIND=dp)                                      :: err, shift, wrong_err
      REAL(KIND=dp), DIMENSION(n)                        :: curv, grad, step, wrong_step

      curv = [0.15_dp, 0.8_dp, 2.4_dp, 6.0_dp, 11.0_dp]
      grad = [0.32_dp, -0.71_dp, 0.14_dp, 0.55_dp, -0.28_dp]

      step(:) = grad(:)/MAX(curv(:), curvature_floor)
      shift = SUM(curv(:)*step(:))/SUM(curv(:))
      step(:) = step(:) - shift
      err = ABS(SUM(curv(:)*step(:)))

      wrong_step(:) = grad(:)/MAX(curv(:), curvature_floor)
      wrong_step(:) = wrong_step(:) - SUM(wrong_step(:))/REAL(n, KIND=dp)
      wrong_err = ABS(SUM(curv(:)*wrong_step(:)))

      IF (err > tol .OR. wrong_err <= 1.0E-3_dp) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,2(1X,ES12.5))') TRIM(label), "failed with errors", err, wrong_err
      ELSE
         WRITE (*, '(A,1X,A,2(1X,ES12.5))') TRIM(label), "OK with errors", err, wrong_err
      END IF

   END SUBROUTINE check_fixed_n_preconditioner_gauge

! **************************************************************************************************
!> \brief Check that the diagonal energy preconditioner removes the fixed-N gauge
!>        in the combined occupation-plus-tether curvature metric.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_fixed_n_tether_preconditioner_gauge(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 6
      REAL(KIND=dp), PARAMETER                           :: curvature_floor = 1.0_dp

      REAL(KIND=dp)                                      :: err, inactive_err, shift, wrong_err
      REAL(KIND=dp), DIMENSION(n)                        :: curv, grad, occ_curv, step, tether_curv, &
                                                            wrong_step

      occ_curv = [0.15_dp, 0.8_dp, 2.4_dp, 6.0_dp, 11.0_dp, 0.0_dp]
      tether_curv = [0.0_dp, 0.35_dp, 0.7_dp, 0.0_dp, 1.6_dp, 0.0_dp]
      grad = [0.32_dp, -0.71_dp, 0.14_dp, 0.55_dp, -0.28_dp, 0.62_dp]
      curv(:) = occ_curv(:) + tether_curv(:)

      step(:) = 0.0_dp
      WHERE (curv(:) > EPSILON(curv(:)))
         step(:) = grad(:)/MAX(curv(:), curvature_floor)
      END WHERE
      shift = SUM(curv(:)*step(:))/SUM(curv(:))
      WHERE (curv(:) > EPSILON(curv(:)))
         step(:) = step(:) - shift
      END WHERE
      err = ABS(SUM(curv(:)*step(:)))
      inactive_err = ABS(step(n))

      wrong_step(:) = 0.0_dp
      WHERE (curv(:) > EPSILON(curv(:)))
         wrong_step(:) = grad(:)/MAX(curv(:), curvature_floor)
      END WHERE
      shift = SUM(occ_curv(:)*wrong_step(:))/SUM(occ_curv(:))
      WHERE (curv(:) > EPSILON(curv(:)))
         wrong_step(:) = wrong_step(:) - shift
      END WHERE
      wrong_err = ABS(SUM(curv(:)*wrong_step(:)))

      IF (err > tol .OR. inactive_err > tol .OR. wrong_err <= 1.0E-3_dp) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,3(1X,ES12.5))') TRIM(label), "failed with errors", &
            err, inactive_err, wrong_err
      ELSE
         WRITE (*, '(A,1X,A,3(1X,ES12.5))') TRIM(label), "OK with errors", &
            err, inactive_err, wrong_err
      END IF

   END SUBROUTINE check_fixed_n_tether_preconditioner_gauge

! **************************************************************************************************
!> \brief Check that the k-point energy preconditioner removes the fixed-N
!>        gauge in the combined irreducible-k-point and canonical-tether
!>        curvature metric.  Dropping either the k-point weights or the tether
!>        curvature leaves a residual chemical-potential component.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_weighted_kpoint_tether_preconditioner_gauge(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 7
      REAL(KIND=dp), PARAMETER                           :: curvature_floor = 1.0_dp

      REAL(KIND=dp)                                      :: err, no_tether_err, no_weight_err, &
                                                            shift, wrong_err
      REAL(KIND=dp), DIMENSION(n)                        :: curv, grad, gvalue, step, tether_curv, &
                                                            wkp, wrong_curv, wrong_step

      wkp = [0.35_dp, 0.35_dp, 0.65_dp, 0.65_dp, 1.0_dp, 0.5_dp, 0.5_dp]
      gvalue = [2.7_dp, 0.4_dp, 6.0_dp, 0.2_dp, 1.5_dp, 0.0_dp, 4.5_dp]
      tether_curv = [0.08_dp, 0.0_dp, 0.30_dp, 0.12_dp, 0.0_dp, 0.0_dp, 0.20_dp]
      grad = [0.42_dp, -0.31_dp, 0.26_dp, 0.75_dp, -0.18_dp, 0.66_dp, 0.11_dp]
      curv(:) = wkp(:)*gvalue(:) + tether_curv(:)

      step(:) = 0.0_dp
      WHERE (curv(:) > EPSILON(curv(:)))
         step(:) = grad(:)/MAX(curv(:), curvature_floor)
      END WHERE
      shift = SUM(curv(:)*step(:))/SUM(curv(:))
      WHERE (curv(:) > EPSILON(curv(:)))
         step(:) = step(:) - shift
      END WHERE
      err = ABS(SUM(curv(:)*step(:)))

      wrong_curv(:) = gvalue(:) + tether_curv(:)
      wrong_step(:) = 0.0_dp
      WHERE (wrong_curv(:) > EPSILON(wrong_curv(:)))
         wrong_step(:) = grad(:)/MAX(wrong_curv(:), curvature_floor)
      END WHERE
      shift = SUM(wrong_curv(:)*wrong_step(:))/SUM(wrong_curv(:))
      WHERE (wrong_curv(:) > EPSILON(wrong_curv(:)))
         wrong_step(:) = wrong_step(:) - shift
      END WHERE
      no_weight_err = ABS(SUM(curv(:)*wrong_step(:)))

      wrong_curv(:) = wkp(:)*gvalue(:)
      wrong_step(:) = 0.0_dp
      WHERE (wrong_curv(:) > EPSILON(wrong_curv(:)))
         wrong_step(:) = grad(:)/MAX(wrong_curv(:), curvature_floor)
      END WHERE
      shift = SUM(wrong_curv(:)*wrong_step(:))/SUM(wrong_curv(:))
      WHERE (wrong_curv(:) > EPSILON(wrong_curv(:)))
         wrong_step(:) = wrong_step(:) - shift
      END WHERE
      no_tether_err = ABS(SUM(curv(:)*wrong_step(:)))

      wrong_step(:) = 0.0_dp
      WHERE (curv(:) > EPSILON(curv(:)))
         wrong_step(:) = grad(:)/MAX(curv(:), curvature_floor)
      END WHERE
      shift = SUM(wrong_step(:))/REAL(n, KIND=dp)
      WHERE (curv(:) > EPSILON(curv(:)))
         wrong_step(:) = wrong_step(:) - shift
      END WHERE
      wrong_err = ABS(SUM(curv(:)*wrong_step(:)))

      IF (err > tol .OR. ABS(step(6)) > tol .OR. no_weight_err <= 1.0E-3_dp .OR. &
          no_tether_err <= 1.0E-3_dp .OR. wrong_err <= 1.0E-3_dp) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,5(1X,ES12.5))') TRIM(label), "failed with errors", &
            err, no_weight_err, no_tether_err, wrong_err, ABS(step(6))
      ELSE
         WRITE (*, '(A,1X,A,5(1X,ES12.5))') TRIM(label), "OK with errors", &
            err, no_weight_err, no_tether_err, wrong_err, ABS(step(6))
      END IF

   END SUBROUTINE check_weighted_kpoint_tether_preconditioner_gauge

! **************************************************************************************************
!> \brief Dense weighted k-point Mermin free-energy model.
!> \param eps ...
!> \param h ...
!> \param wk ...
!> \param nel ...
!> \param sigma ...
!> \param maxocc ...
!> \return ...
! **************************************************************************************************
   FUNCTION dense_kpoint_mermin_value(eps, h, wk, nel, sigma, maxocc) RESULT(value)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: eps, h
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: wk
      REAL(KIND=dp), INTENT(IN)                          :: nel, sigma, maxocc
      REAL(KIND=dp)                                      :: value

      INTEGER                                            :: ik
      REAL(KIND=dp)                                      :: kTS, mu
      REAL(KIND=dp), &
         DIMENSION(SIZE(eps, 1), SIZE(eps, 2))           :: occ

      CALL dense_fermi_kpoint(occ, mu, kTS, eps, nel, wk, sigma, maxocc)
      value = kTS
      DO ik = 1, SIZE(eps, 2)
         value = value + wk(ik)*SUM(occ(:, ik)*h(:, ik))
      END DO

   END FUNCTION dense_kpoint_mermin_value

! **************************************************************************************************
!> \brief Dense weighted fixed-N Mermin gradient.
!> \param eps ...
!> \param h ...
!> \param wk ...
!> \param nel ...
!> \param sigma ...
!> \param maxocc ...
!> \param grad ...
! **************************************************************************************************
   SUBROUTINE dense_kpoint_mermin_gradient(eps, h, wk, nel, sigma, maxocc, grad)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: eps, h
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: wk
      REAL(KIND=dp), INTENT(IN)                          :: nel, sigma, maxocc
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: grad

      INTEGER                                            :: ik, imo
      REAL(KIND=dp)                                      :: gdotv, gsum, gvalue, kTS, mu, v
      REAL(KIND=dp), &
         DIMENSION(SIZE(eps, 1), SIZE(eps, 2))           :: occ

      CALL dense_fermi_kpoint(occ, mu, kTS, eps, nel, wk, sigma, maxocc)

      gsum = 0.0_dp
      gdotv = 0.0_dp
      DO ik = 1, SIZE(eps, 2)
         DO imo = 1, SIZE(eps, 1)
            gvalue = maxocc*(occ(imo, ik)/maxocc)*(1.0_dp - occ(imo, ik)/maxocc)/sigma
            v = h(imo, ik) - eps(imo, ik) + mu
            gsum = gsum + wk(ik)*gvalue
            gdotv = gdotv + wk(ik)*gvalue*v
         END DO
      END DO

      grad(:, :) = 0.0_dp
      IF (ABS(gsum) > EPSILON(gsum)) THEN
         DO ik = 1, SIZE(eps, 2)
            DO imo = 1, SIZE(eps, 1)
               gvalue = maxocc*(occ(imo, ik)/maxocc)*(1.0_dp - occ(imo, ik)/maxocc)/sigma
               v = h(imo, ik) - eps(imo, ik) + mu
               grad(imo, ik) = wk(ik)*gvalue*(gdotv/gsum - v)
            END DO
         END DO
      END IF

   END SUBROUTINE dense_kpoint_mermin_gradient

! **************************************************************************************************
!> \brief Check the coupled complex rotation and auxiliary-energy Mermin gradient.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_complex_kpoint_mermin_coupled_gradient(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 3, nkp = 2
      REAL(KIND=dp), PARAMETER                           :: fd_step = 1.0E-6_dp, maxocc = 2.0_dp, &
                                                            nel = 3.20_dp, sigma = 0.035_dp

      COMPLEX(KIND=dp), DIMENSION(n, n, nkp)             :: base
      INTEGER                                            :: i, ik, j
      REAL(KIND=dp)                                      :: err, fd, value_m, value_p
      REAL(KIND=dp), DIMENSION(n, n, nkp)                :: grad_im, grad_re, x_im, x_m_im, x_m_re, &
                                                            x_p_im, x_p_re, x_re
      REAL(KIND=dp), DIMENSION(n, nkp)                   :: eps, eps_m, eps_p, grad_eps
      REAL(KIND=dp), DIMENSION(nkp)                      :: wk

      wk = [0.40_dp, 0.60_dp]
      eps(:, 1) = [-0.23_dp, -0.035_dp, 0.090_dp]
      eps(:, 2) = [-0.19_dp, -0.010_dp, 0.125_dp]

      base(:, :, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      base(1, 1, 1) = CMPLX(-0.26_dp, 0.0_dp, KIND=dp)
      base(2, 2, 1) = CMPLX(-0.015_dp, 0.0_dp, KIND=dp)
      base(3, 3, 1) = CMPLX(0.14_dp, 0.0_dp, KIND=dp)
      base(1, 2, 1) = CMPLX(0.045_dp, -0.025_dp, KIND=dp)
      base(1, 3, 1) = CMPLX(-0.030_dp, 0.035_dp, KIND=dp)
      base(2, 3, 1) = CMPLX(0.020_dp, 0.030_dp, KIND=dp)
      base(2, 1, 1) = CONJG(base(1, 2, 1))
      base(3, 1, 1) = CONJG(base(1, 3, 1))
      base(3, 2, 1) = CONJG(base(2, 3, 1))

      base(1, 1, 2) = CMPLX(-0.22_dp, 0.0_dp, KIND=dp)
      base(2, 2, 2) = CMPLX(0.005_dp, 0.0_dp, KIND=dp)
      base(3, 3, 2) = CMPLX(0.18_dp, 0.0_dp, KIND=dp)
      base(1, 2, 2) = CMPLX(-0.035_dp, 0.015_dp, KIND=dp)
      base(1, 3, 2) = CMPLX(0.025_dp, 0.020_dp, KIND=dp)
      base(2, 3, 2) = CMPLX(0.030_dp, -0.025_dp, KIND=dp)
      base(2, 1, 2) = CONJG(base(1, 2, 2))
      base(3, 1, 2) = CONJG(base(1, 3, 2))
      base(3, 2, 2) = CONJG(base(2, 3, 2))

      x_re(:, :, :) = 0.0_dp
      x_re(1, 2, 1) = 0.055_dp
      x_re(1, 3, 1) = -0.025_dp
      x_re(2, 3, 1) = 0.035_dp
      x_re(1, 2, 2) = -0.040_dp
      x_re(1, 3, 2) = 0.030_dp
      x_re(2, 3, 2) = 0.020_dp
      DO ik = 1, nkp
         DO j = 2, n
            DO i = 1, j - 1
               x_re(j, i, ik) = -x_re(i, j, ik)
            END DO
         END DO
      END DO

      x_im(:, :, :) = 0.0_dp
      x_im(1, 1, 1) = 0.020_dp
      x_im(2, 2, 1) = -0.035_dp
      x_im(3, 3, 1) = 0.015_dp
      x_im(1, 2, 1) = -0.018_dp
      x_im(1, 3, 1) = 0.012_dp
      x_im(2, 3, 1) = 0.026_dp
      x_im(1, 1, 2) = -0.015_dp
      x_im(2, 2, 2) = 0.030_dp
      x_im(3, 3, 2) = -0.010_dp
      x_im(1, 2, 2) = 0.014_dp
      x_im(1, 3, 2) = -0.022_dp
      x_im(2, 3, 2) = 0.018_dp
      DO ik = 1, nkp
         DO j = 2, n
            DO i = 1, j - 1
               x_im(j, i, ik) = x_im(i, j, ik)
            END DO
         END DO
      END DO

      CALL dense_complex_kpoint_mermin_gradient(x_re, x_im, base, eps, wk, nel, sigma, maxocc, &
                                                grad_re, grad_im, grad_eps, value_p)

      err = 0.0_dp
      DO ik = 1, nkp
         DO i = 1, n
            eps_p(:, :) = eps(:, :)
            eps_m(:, :) = eps(:, :)
            eps_p(i, ik) = eps_p(i, ik) + fd_step
            eps_m(i, ik) = eps_m(i, ik) - fd_step
            fd = (dense_complex_kpoint_mermin_value(x_re, x_im, base, eps_p, wk, nel, sigma, maxocc) - &
                  dense_complex_kpoint_mermin_value(x_re, x_im, base, eps_m, wk, nel, sigma, maxocc))/ &
                 (2.0_dp*fd_step)
            err = MAX(err, ABS(fd - grad_eps(i, ik)))
         END DO
      END DO

      DO ik = 1, nkp
         DO j = 2, n
            DO i = 1, j - 1
               x_p_re(:, :, :) = x_re(:, :, :)
               x_m_re(:, :, :) = x_re(:, :, :)
               x_p_im(:, :, :) = x_im(:, :, :)
               x_m_im(:, :, :) = x_im(:, :, :)
               x_p_re(i, j, ik) = x_p_re(i, j, ik) + fd_step
               x_p_re(j, i, ik) = x_p_re(j, i, ik) - fd_step
               x_m_re(i, j, ik) = x_m_re(i, j, ik) - fd_step
               x_m_re(j, i, ik) = x_m_re(j, i, ik) + fd_step
               value_p = dense_complex_kpoint_mermin_value(x_p_re, x_p_im, base, eps, wk, nel, &
                                                           sigma, maxocc)
               value_m = dense_complex_kpoint_mermin_value(x_m_re, x_m_im, base, eps, wk, nel, &
                                                           sigma, maxocc)
               fd = (value_p - value_m)/(2.0_dp*fd_step)
               err = MAX(err, ABS(fd - grad_re(i, j, ik)))

               x_p_re(:, :, :) = x_re(:, :, :)
               x_m_re(:, :, :) = x_re(:, :, :)
               x_p_im(:, :, :) = x_im(:, :, :)
               x_m_im(:, :, :) = x_im(:, :, :)
               x_p_im(i, j, ik) = x_p_im(i, j, ik) + fd_step
               x_p_im(j, i, ik) = x_p_im(j, i, ik) + fd_step
               x_m_im(i, j, ik) = x_m_im(i, j, ik) - fd_step
               x_m_im(j, i, ik) = x_m_im(j, i, ik) - fd_step
               value_p = dense_complex_kpoint_mermin_value(x_p_re, x_p_im, base, eps, wk, nel, &
                                                           sigma, maxocc)
               value_m = dense_complex_kpoint_mermin_value(x_m_re, x_m_im, base, eps, wk, nel, &
                                                           sigma, maxocc)
               fd = (value_p - value_m)/(2.0_dp*fd_step)
               err = MAX(err, ABS(fd - grad_im(i, j, ik)))
            END DO
         END DO

         DO i = 1, n
            x_p_re(:, :, :) = x_re(:, :, :)
            x_m_re(:, :, :) = x_re(:, :, :)
            x_p_im(:, :, :) = x_im(:, :, :)
            x_m_im(:, :, :) = x_im(:, :, :)
            x_p_im(i, i, ik) = x_p_im(i, i, ik) + fd_step
            x_m_im(i, i, ik) = x_m_im(i, i, ik) - fd_step
            value_p = dense_complex_kpoint_mermin_value(x_p_re, x_p_im, base, eps, wk, nel, &
                                                        sigma, maxocc)
            value_m = dense_complex_kpoint_mermin_value(x_m_re, x_m_im, base, eps, wk, nel, &
                                                        sigma, maxocc)
            fd = (value_p - value_m)/(2.0_dp*fd_step)
            err = MAX(err, ABS(fd - 0.5_dp*grad_im(i, i, ik)))
         END DO
      END DO

      IF (err > tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,ES12.5)') TRIM(label), "failed with error", err
      ELSE
         WRITE (*, '(A,1X,A,1X,ES12.5)') TRIM(label), "OK with error", err
      END IF

   END SUBROUTINE check_complex_kpoint_mermin_coupled_gradient

! **************************************************************************************************
!> \brief Check the reduced fixed-N Schur Hessian against finite differences.
!>        The reduced problem uses eps=h(x), so the dense formula must include
!>        the fixed-occupation rotation Hessian and the fixed-N energy response
!>        projector, H_red = H_xx - R^T P_N R.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_complex_fixed_n_schur_hessian(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 3, nkp = 1, nrot = n*(n - 1)
      REAL(KIND=dp), PARAMETER                           :: fd_step = 1.0E-5_dp, maxocc = 2.0_dp, &
                                                            nel = 3.15_dp, sigma = 0.045_dp

      COMPLEX(KIND=dp), DIMENSION(n, n, nkp)             :: base
      INTEGER                                            :: i, ivar, j, r, s
      INTEGER, DIMENSION(nrot)                           :: var_i, var_j, var_kind
      REAL(KIND=dp)                                      :: dense_term, err, kTS, mu, &
                                                            projected_term, rtd_r, rtd_s, value
      REAL(KIND=dp), DIMENSION(n)                        :: curv, response_r, response_s
      REAL(KIND=dp), DIMENSION(n, n, nkp)                :: grad_im, grad_re, x_im, x_m_im, x_m_re, &
                                                            x_p_im, x_p_re, x_re
      REAL(KIND=dp), DIMENSION(n, nkp)                   :: eps_base, eps_m, eps_p, grad_eps, &
                                                            h_base, h_m, h_p, occ
      REAL(KIND=dp), DIMENSION(nkp)                      :: wk
      REAL(KIND=dp), DIMENSION(nrot)                     :: grad_m_vec, grad_p_vec
      REAL(KIND=dp), DIMENSION(nrot, nrot)               :: hfixed, hred, schur

      wk = [1.0_dp]
      base(:, :, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      base(1, 1, 1) = CMPLX(-0.24_dp, 0.0_dp, KIND=dp)
      base(2, 2, 1) = CMPLX(0.01_dp, 0.0_dp, KIND=dp)
      base(3, 3, 1) = CMPLX(0.19_dp, 0.0_dp, KIND=dp)
      base(1, 2, 1) = CMPLX(0.045_dp, -0.025_dp, KIND=dp)
      base(1, 3, 1) = CMPLX(-0.025_dp, 0.030_dp, KIND=dp)
      base(2, 3, 1) = CMPLX(0.035_dp, 0.020_dp, KIND=dp)
      base(2, 1, 1) = CONJG(base(1, 2, 1))
      base(3, 1, 1) = CONJG(base(1, 3, 1))
      base(3, 2, 1) = CONJG(base(2, 3, 1))

      x_re(:, :, :) = 0.0_dp
      x_re(1, 2, 1) = 0.050_dp
      x_re(1, 3, 1) = -0.030_dp
      x_re(2, 3, 1) = 0.025_dp
      DO j = 2, n
         DO i = 1, j - 1
            x_re(j, i, 1) = -x_re(i, j, 1)
         END DO
      END DO

      x_im(:, :, :) = 0.0_dp
      x_im(1, 1, 1) = 0.015_dp
      x_im(2, 2, 1) = -0.025_dp
      x_im(3, 3, 1) = 0.010_dp
      x_im(1, 2, 1) = -0.020_dp
      x_im(1, 3, 1) = 0.015_dp
      x_im(2, 3, 1) = 0.018_dp
      DO j = 2, n
         DO i = 1, j - 1
            x_im(j, i, 1) = x_im(i, j, 1)
         END DO
      END DO

      ivar = 0
      DO j = 2, n
         DO i = 1, j - 1
            ivar = ivar + 1
            var_i(ivar) = i
            var_j(ivar) = j
            var_kind(ivar) = 1
            ivar = ivar + 1
            var_i(ivar) = i
            var_j(ivar) = j
            var_kind(ivar) = 2
         END DO
      END DO

      CALL dense_complex_kpoint_hdiag(x_re, x_im, base, h_base)
      eps_base(:, :) = h_base(:, :)
      CALL dense_fermi_kpoint(occ, mu, kTS, eps_base, nel, wk, sigma, maxocc)
      DO i = 1, n
         curv(i) = wk(1)*maxocc*(occ(i, 1)/maxocc)*(1.0_dp - occ(i, 1)/maxocc)/sigma
      END DO

      DO s = 1, nrot
         x_p_re(:, :, :) = x_re(:, :, :)
         x_m_re(:, :, :) = x_re(:, :, :)
         x_p_im(:, :, :) = x_im(:, :, :)
         x_m_im(:, :, :) = x_im(:, :, :)
         CALL perturb_complex_rotation_var(x_p_re, x_p_im, var_i(s), var_j(s), var_kind(s), fd_step)
         CALL perturb_complex_rotation_var(x_m_re, x_m_im, var_i(s), var_j(s), var_kind(s), -fd_step)

         CALL dense_complex_kpoint_mermin_gradient(x_p_re, x_p_im, base, eps_base, wk, nel, &
                                                   sigma, maxocc, grad_re, grad_im, grad_eps, value)
         CALL extract_rotation_gradient_vector(grad_re, grad_im, var_i, var_j, var_kind, grad_p_vec)
         CALL dense_complex_kpoint_mermin_gradient(x_m_re, x_m_im, base, eps_base, wk, nel, &
                                                   sigma, maxocc, grad_re, grad_im, grad_eps, value)
         CALL extract_rotation_gradient_vector(grad_re, grad_im, var_i, var_j, var_kind, grad_m_vec)
         hfixed(:, s) = (grad_p_vec(:) - grad_m_vec(:))/(2.0_dp*fd_step)

         CALL dense_complex_kpoint_hdiag(x_p_re, x_p_im, base, h_p)
         CALL dense_complex_kpoint_hdiag(x_m_re, x_m_im, base, h_m)
         eps_p(:, :) = h_p(:, :)
         eps_m(:, :) = h_m(:, :)
         response_s(:) = (h_p(:, 1) - h_m(:, 1))/(2.0_dp*fd_step)

         CALL dense_complex_kpoint_mermin_gradient(x_p_re, x_p_im, base, eps_p, wk, nel, &
                                                   sigma, maxocc, grad_re, grad_im, grad_eps, value)
         CALL extract_rotation_gradient_vector(grad_re, grad_im, var_i, var_j, var_kind, grad_p_vec)
         CALL dense_complex_kpoint_mermin_gradient(x_m_re, x_m_im, base, eps_m, wk, nel, &
                                                   sigma, maxocc, grad_re, grad_im, grad_eps, value)
         CALL extract_rotation_gradient_vector(grad_re, grad_im, var_i, var_j, var_kind, grad_m_vec)
         hred(:, s) = (grad_p_vec(:) - grad_m_vec(:))/(2.0_dp*fd_step)

         DO r = 1, nrot
            response_r(:) = 0.0_dp
            x_p_re(:, :, :) = x_re(:, :, :)
            x_m_re(:, :, :) = x_re(:, :, :)
            x_p_im(:, :, :) = x_im(:, :, :)
            x_m_im(:, :, :) = x_im(:, :, :)
            CALL perturb_complex_rotation_var(x_p_re, x_p_im, var_i(r), var_j(r), var_kind(r), fd_step)
            CALL perturb_complex_rotation_var(x_m_re, x_m_im, var_i(r), var_j(r), var_kind(r), -fd_step)
            CALL dense_complex_kpoint_hdiag(x_p_re, x_p_im, base, h_p)
            CALL dense_complex_kpoint_hdiag(x_m_re, x_m_im, base, h_m)
            response_r(:) = (h_p(:, 1) - h_m(:, 1))/(2.0_dp*fd_step)
            dense_term = DOT_PRODUCT(curv(:)*response_r(:), response_s(:))
            rtd_r = DOT_PRODUCT(curv(:), response_r(:))
            rtd_s = DOT_PRODUCT(curv(:), response_s(:))
            projected_term = dense_term - rtd_r*rtd_s/SUM(curv(:))
            schur(r, s) = hfixed(r, s) - projected_term
         END DO
      END DO

      err = MAXVAL(ABS(schur(:, :) - hred(:, :)))

      IF (err > tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,ES12.5)') TRIM(label), "failed with error", err
      ELSE
         WRITE (*, '(A,1X,A,1X,ES12.5)') TRIM(label), "OK with error", err
      END IF

   END SUBROUTINE check_complex_fixed_n_schur_hessian

! **************************************************************************************************
!> \brief Check the Schur response channel against finite fixed-N Mermin occupation response.
!>        The cross block is d grad_eps / d theta = -P_N R.  The production
!>        preconditioner uses the projected pseudo-inverse R = -P_N^+ d grad_eps/d theta,
!>        so null occupation-curvature modes cannot drive the rotation RHS.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_complex_fixed_n_response_channel(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 3, nkp = 1, nrot = n*(n - 1)
      REAL(KIND=dp), PARAMETER                           :: fd_step = 1.0E-6_dp, maxocc = 2.0_dp, &
                                                            nel = 3.05_dp, sigma = 0.080_dp

      COMPLEX(KIND=dp), DIMENSION(n, n, nkp)             :: base
      INTEGER                                            :: i, ivar, j, s
      INTEGER, DIMENSION(nrot)                           :: var_i, var_j, var_kind
      REAL(KIND=dp)                                      :: curv_sum, err, kTS, mu, response_shift
      REAL(KIND=dp), DIMENSION(n)                        :: cross_fd, curv, dh, projected_response, &
                                                            ref_cross, response_from_cross
      REAL(KIND=dp), DIMENSION(n, n, nkp)                :: x_im, x_m_im, x_m_re, x_p_im, x_p_re, &
                                                            x_re
      REAL(KIND=dp), DIMENSION(n, nkp)                   :: eps_base, grad_eps_m, grad_eps_p, &
                                                            h_base, h_m, h_p, occ
      REAL(KIND=dp), DIMENSION(nkp)                      :: wk

      wk = [1.0_dp]
      base(:, :, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      base(1, 1, 1) = CMPLX(-0.20_dp, 0.0_dp, KIND=dp)
      base(2, 2, 1) = CMPLX(0.02_dp, 0.0_dp, KIND=dp)
      base(3, 3, 1) = CMPLX(0.16_dp, 0.0_dp, KIND=dp)
      base(1, 2, 1) = CMPLX(0.050_dp, -0.030_dp, KIND=dp)
      base(1, 3, 1) = CMPLX(-0.020_dp, 0.025_dp, KIND=dp)
      base(2, 3, 1) = CMPLX(0.030_dp, 0.035_dp, KIND=dp)
      base(2, 1, 1) = CONJG(base(1, 2, 1))
      base(3, 1, 1) = CONJG(base(1, 3, 1))
      base(3, 2, 1) = CONJG(base(2, 3, 1))

      x_re(:, :, :) = 0.0_dp
      x_re(1, 2, 1) = 0.045_dp
      x_re(1, 3, 1) = -0.020_dp
      x_re(2, 3, 1) = 0.030_dp
      DO j = 2, n
         DO i = 1, j - 1
            x_re(j, i, 1) = -x_re(i, j, 1)
         END DO
      END DO

      x_im(:, :, :) = 0.0_dp
      x_im(1, 1, 1) = 0.010_dp
      x_im(2, 2, 1) = -0.020_dp
      x_im(3, 3, 1) = 0.018_dp
      x_im(1, 2, 1) = -0.015_dp
      x_im(1, 3, 1) = 0.020_dp
      x_im(2, 3, 1) = 0.012_dp
      DO j = 2, n
         DO i = 1, j - 1
            x_im(j, i, 1) = x_im(i, j, 1)
         END DO
      END DO

      ivar = 0
      DO j = 2, n
         DO i = 1, j - 1
            ivar = ivar + 1
            var_i(ivar) = i
            var_j(ivar) = j
            var_kind(ivar) = 1
            ivar = ivar + 1
            var_i(ivar) = i
            var_j(ivar) = j
            var_kind(ivar) = 2
         END DO
      END DO

      CALL dense_complex_kpoint_hdiag(x_re, x_im, base, h_base)
      eps_base(:, 1) = h_base(:, 1) + [0.018_dp, -0.012_dp, 0.006_dp]
      CALL dense_fermi_kpoint(occ, mu, kTS, eps_base, nel, wk, sigma, maxocc)
      DO i = 1, n
         curv(i) = wk(1)*maxocc*(occ(i, 1)/maxocc)*(1.0_dp - occ(i, 1)/maxocc)/sigma
      END DO
      curv_sum = SUM(curv(:))

      err = 0.0_dp
      DO s = 1, nrot
         x_p_re(:, :, :) = x_re(:, :, :)
         x_m_re(:, :, :) = x_re(:, :, :)
         x_p_im(:, :, :) = x_im(:, :, :)
         x_m_im(:, :, :) = x_im(:, :, :)
         CALL perturb_complex_rotation_var(x_p_re, x_p_im, var_i(s), var_j(s), var_kind(s), fd_step)
         CALL perturb_complex_rotation_var(x_m_re, x_m_im, var_i(s), var_j(s), var_kind(s), -fd_step)

         CALL dense_complex_kpoint_hdiag(x_p_re, x_p_im, base, h_p)
         CALL dense_complex_kpoint_hdiag(x_m_re, x_m_im, base, h_m)
         CALL dense_kpoint_mermin_gradient(eps_base, h_p, wk, nel, sigma, maxocc, grad_eps_p)
         CALL dense_kpoint_mermin_gradient(eps_base, h_m, wk, nel, sigma, maxocc, grad_eps_m)

         dh(:) = (h_p(:, 1) - h_m(:, 1))/(2.0_dp*fd_step)
         response_shift = SUM(curv(:)*dh(:))/curv_sum
         projected_response(:) = dh(:) - response_shift
         ref_cross(:) = curv(:)*(response_shift - dh(:))
         cross_fd(:) = (grad_eps_p(:, 1) - grad_eps_m(:, 1))/(2.0_dp*fd_step)

         response_from_cross(:) = 0.0_dp
         DO i = 1, n
            IF (curv(i) > EPSILON(curv(i))) response_from_cross(i) = -cross_fd(i)/curv(i)
         END DO
         response_shift = SUM(curv(:)*response_from_cross(:))/curv_sum
         response_from_cross(:) = response_from_cross(:) - response_shift

         err = MAX(err, MAXVAL(ABS(cross_fd(:) - ref_cross(:))))
         err = MAX(err, MAXVAL(ABS(response_from_cross(:) - projected_response(:))))
      END DO

      IF (err > tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,ES12.5)') TRIM(label), "failed with error", err
      ELSE
         WRITE (*, '(A,1X,A,1X,ES12.5)') TRIM(label), "OK with error", err
      END IF

   END SUBROUTINE check_complex_fixed_n_response_channel

! **************************************************************************************************
!> \brief Check the response blend used by the finite complex split Schur model.
!>        The occupation-projector response is the derivative of the rotation
!>        gradient with respect to the band occupation weights.  For a purely
!>        diagonal projector this reduces to d[U^H H U]_ii/dtheta; with a
!>        spectral projector blend it also contains the finite complex projector
!>        rotation response.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_complex_projector_blended_response_channel(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 3, nkp = 1, nrot = n*(n - 1)
      REAL(KIND=dp), PARAMETER                           :: active_tol = 1.0E-4_dp, &
                                                            fd_occ_step = 1.0E-6_dp, &
                                                            fd_step = 1.0E-6_dp, &
                                                            projector_weight = 0.35_dp

      COMPLEX(KIND=dp), DIMENSION(n, n, nkp)             :: base
      INTEGER                                            :: i, ivar, j, s
      INTEGER, DIMENSION(nrot)                           :: var_i, var_j, var_kind
      REAL(KIND=dp)                                      :: active_err, err
      REAL(KIND=dp), DIMENSION(n)                        :: occ_scale, occ_scale_m, occ_scale_p
      REAL(KIND=dp), DIMENSION(n, n)                     :: grad_im_2d, grad_re_2d
      REAL(KIND=dp), DIMENSION(n, n, nkp)                :: x_im, x_m_im, x_m_re, x_p_im, x_p_re, &
                                                            x_re
      REAL(KIND=dp), DIMENSION(n, nkp)                   :: h_m, h_p
      REAL(KIND=dp), DIMENSION(n, nrot)                  :: diag_response, projector_response, &
                                                            response_fd, zero_response
      REAL(KIND=dp), DIMENSION(nrot)                     :: grad_m_vec, grad_p_vec

      base(:, :, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      base(1, 1, 1) = CMPLX(-0.28_dp, 0.0_dp, KIND=dp)
      base(2, 2, 1) = CMPLX(0.03_dp, 0.0_dp, KIND=dp)
      base(3, 3, 1) = CMPLX(0.22_dp, 0.0_dp, KIND=dp)
      base(1, 2, 1) = CMPLX(0.060_dp, -0.020_dp, KIND=dp)
      base(1, 3, 1) = CMPLX(-0.025_dp, 0.035_dp, KIND=dp)
      base(2, 3, 1) = CMPLX(0.040_dp, 0.025_dp, KIND=dp)
      base(2, 1, 1) = CONJG(base(1, 2, 1))
      base(3, 1, 1) = CONJG(base(1, 3, 1))
      base(3, 2, 1) = CONJG(base(2, 3, 1))

      x_re(:, :, :) = 0.0_dp
      x_re(1, 2, 1) = 0.040_dp
      x_re(1, 3, 1) = -0.035_dp
      x_re(2, 3, 1) = 0.020_dp
      DO j = 2, n
         DO i = 1, j - 1
            x_re(j, i, 1) = -x_re(i, j, 1)
         END DO
      END DO

      x_im(:, :, :) = 0.0_dp
      x_im(1, 1, 1) = 0.018_dp
      x_im(2, 2, 1) = -0.022_dp
      x_im(3, 3, 1) = 0.014_dp
      x_im(1, 2, 1) = -0.018_dp
      x_im(1, 3, 1) = 0.016_dp
      x_im(2, 3, 1) = 0.026_dp
      DO j = 2, n
         DO i = 1, j - 1
            x_im(j, i, 1) = x_im(i, j, 1)
         END DO
      END DO

      ivar = 0
      DO j = 2, n
         DO i = 1, j - 1
            ivar = ivar + 1
            var_i(ivar) = i
            var_j(ivar) = j
            var_kind(ivar) = 1
            ivar = ivar + 1
            var_i(ivar) = i
            var_j(ivar) = j
            var_kind(ivar) = 2
         END DO
      END DO

      occ_scale(:) = [1.42_dp, 0.77_dp, 0.31_dp]

      CALL dense_blended_projector_response(x_re(:, :, 1), x_im(:, :, 1), base(:, :, 1), &
                                            0.0_dp, var_i, var_j, var_kind, zero_response)
      CALL dense_blended_projector_response(x_re(:, :, 1), x_im(:, :, 1), base(:, :, 1), &
                                            projector_weight, var_i, var_j, var_kind, projector_response)
      response_fd(:, :) = 0.0_dp
      err = 0.0_dp
      active_err = 0.0_dp
      DO s = 1, nrot
         x_p_re(:, :, :) = x_re(:, :, :)
         x_m_re(:, :, :) = x_re(:, :, :)
         x_p_im(:, :, :) = x_im(:, :, :)
         x_m_im(:, :, :) = x_im(:, :, :)
         CALL perturb_complex_rotation_var(x_p_re, x_p_im, var_i(s), var_j(s), var_kind(s), fd_step)
         CALL perturb_complex_rotation_var(x_m_re, x_m_im, var_i(s), var_j(s), var_kind(s), -fd_step)

         CALL dense_complex_kpoint_hdiag(x_p_re, x_p_im, base, h_p)
         CALL dense_complex_kpoint_hdiag(x_m_re, x_m_im, base, h_m)
         diag_response(:, s) = (h_p(:, 1) - h_m(:, 1))/(2.0_dp*fd_step)
         err = MAX(err, MAXVAL(ABS(zero_response(:, s) - diag_response(:, s))))
         active_err = MAX(active_err, MAXVAL(ABS(projector_response(:, s) - diag_response(:, s))))
      END DO

      DO i = 1, n
         occ_scale_p(:) = occ_scale(:)
         occ_scale_m(:) = occ_scale(:)
         occ_scale_p(i) = occ_scale_p(i) + fd_occ_step
         occ_scale_m(i) = occ_scale_m(i) - fd_occ_step
         CALL dense_blended_projector_rotation_gradient(x_re(:, :, 1), x_im(:, :, 1), base(:, :, 1), &
                                                        occ_scale_p, projector_weight, grad_re_2d, &
                                                        grad_im_2d, h_p(:, 1))
         CALL extract_rotation_gradient_vector_2d(grad_re_2d, grad_im_2d, var_i, var_j, &
                                                  var_kind, grad_p_vec)
         CALL dense_blended_projector_rotation_gradient(x_re(:, :, 1), x_im(:, :, 1), base(:, :, 1), &
                                                        occ_scale_m, projector_weight, grad_re_2d, &
                                                        grad_im_2d, h_m(:, 1))
         CALL extract_rotation_gradient_vector_2d(grad_re_2d, grad_im_2d, var_i, var_j, &
                                                  var_kind, grad_m_vec)
         response_fd(i, :) = 2.0_dp*(grad_p_vec(:) - grad_m_vec(:))/(2.0_dp*fd_occ_step)
      END DO

      err = MAX(err, MAXVAL(ABS(projector_response(:, :) - response_fd(:, :))))

      IF (err > tol .OR. active_err < active_tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,2(1X,ES12.5))') TRIM(label), "failed with error/active", err, active_err
      ELSE
         WRITE (*, '(A,1X,A,2(1X,ES12.5))') TRIM(label), "OK with error/active", err, active_err
      END IF

   END SUBROUTINE check_complex_projector_blended_response_channel

! **************************************************************************************************
!> \brief Check the finite-complex physical occupation response against eigenvalue differences.
!>        This is the gauge-stable companion to the band-projector basis response:
!>        the diagonal part remains the Rayleigh response, while the spectral part
!>        follows finite eigenvalue shifts rather than the instantaneous projector basis.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_complex_projector_physical_response_channel(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 3, nkp = 1, nrot = n*(n - 1)
      REAL(KIND=dp), PARAMETER                           :: active_tol = 1.0E-5_dp, &
                                                            fd_step = 1.0E-6_dp, &
                                                            projector_weight = 0.70_dp

      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: chc_m, chc_p, occ_mat_m, occ_mat_p, u
      COMPLEX(KIND=dp), DIMENSION(n, n, nkp)             :: base
      INTEGER                                            :: i, ivar, j, s
      INTEGER, DIMENSION(nrot)                           :: var_i, var_j, var_kind
      REAL(KIND=dp)                                      :: active_err, active_projector_weight, &
                                                            curv_max, err
      REAL(KIND=dp), DIMENSION(n)                        :: occ_curv, occ_scale
      REAL(KIND=dp), DIMENSION(n, n, nkp)                :: x_im, x_m_im, x_m_re, x_p_im, x_p_re, &
                                                            x_re
      REAL(KIND=dp), DIMENSION(n, nkp)                   :: diag_m, diag_p, eig_m, eig_p
      REAL(KIND=dp), DIMENSION(n, nrot)                  :: basis_response, physical_response, &
                                                            reference

      base(:, :, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      base(1, 1, 1) = CMPLX(-0.015_dp, 0.0_dp, KIND=dp)
      base(2, 2, 1) = CMPLX(0.011_dp, 0.0_dp, KIND=dp)
      base(3, 3, 1) = CMPLX(0.38_dp, 0.0_dp, KIND=dp)
      base(1, 2, 1) = CMPLX(0.020_dp, -0.035_dp, KIND=dp)
      base(1, 3, 1) = CMPLX(-0.045_dp, 0.015_dp, KIND=dp)
      base(2, 3, 1) = CMPLX(0.030_dp, 0.040_dp, KIND=dp)
      base(2, 1, 1) = CONJG(base(1, 2, 1))
      base(3, 1, 1) = CONJG(base(1, 3, 1))
      base(3, 2, 1) = CONJG(base(2, 3, 1))

      x_re(:, :, :) = 0.0_dp
      x_re(1, 2, 1) = 0.055_dp
      x_re(1, 3, 1) = -0.025_dp
      x_re(2, 3, 1) = 0.030_dp
      DO j = 2, n
         DO i = 1, j - 1
            x_re(j, i, 1) = -x_re(i, j, 1)
         END DO
      END DO

      x_im(:, :, :) = 0.0_dp
      x_im(1, 1, 1) = 0.012_dp
      x_im(2, 2, 1) = -0.018_dp
      x_im(3, 3, 1) = 0.021_dp
      x_im(1, 2, 1) = -0.028_dp
      x_im(1, 3, 1) = 0.019_dp
      x_im(2, 3, 1) = 0.024_dp
      DO j = 2, n
         DO i = 1, j - 1
            x_im(j, i, 1) = x_im(i, j, 1)
         END DO
      END DO

      ivar = 0
      DO j = 2, n
         DO i = 1, j - 1
            ivar = ivar + 1
            var_i(ivar) = i
            var_j(ivar) = j
            var_kind(ivar) = 1
            ivar = ivar + 1
            var_i(ivar) = i
            var_j(ivar) = j
            var_kind(ivar) = 2
         END DO
      END DO

      occ_curv(:) = [0.31_dp, 0.72_dp, 0.19_dp]
      occ_scale(:) = [1.42_dp, 0.77_dp, 0.31_dp]
      curv_max = MAXVAL(occ_curv)
      CALL dense_blended_projector_response(x_re(:, :, 1), x_im(:, :, 1), base(:, :, 1), &
                                            projector_weight, var_i, var_j, var_kind, basis_response)
      CALL dense_blended_projector_physical_response(x_re(:, :, 1), x_im(:, :, 1), base(:, :, 1), &
                                                     projector_weight, occ_scale, occ_curv, var_i, var_j, var_kind, &
                                                     physical_response)
      DO s = 1, nrot
         x_p_re(:, :, :) = x_re(:, :, :)
         x_m_re(:, :, :) = x_re(:, :, :)
         x_p_im(:, :, :) = x_im(:, :, :)
         x_m_im(:, :, :) = x_im(:, :, :)
         CALL perturb_complex_rotation_var(x_p_re, x_p_im, var_i(s), var_j(s), var_kind(s), fd_step)
         CALL perturb_complex_rotation_var(x_m_re, x_m_im, var_i(s), var_j(s), var_kind(s), -fd_step)
         CALL dense_complex_kpoint_hdiag(x_p_re, x_p_im, base, diag_p)
         CALL dense_complex_kpoint_hdiag(x_m_re, x_m_im, base, diag_m)
         CALL dense_complex_kpoint_eigvals(x_p_re, x_p_im, base, eig_p)
         CALL dense_complex_kpoint_eigvals(x_m_re, x_m_im, base, eig_m)
         CALL dense_rotation_u(x_p_re(:, :, 1), x_p_im(:, :, 1), u)
         chc_p(:, :) = MATMUL(CONJG(TRANSPOSE(u)), MATMUL(base(:, :, 1), u))
         CALL dense_rotation_u(x_m_re(:, :, 1), x_m_im(:, :, 1), u)
         chc_m(:, :) = MATMUL(CONJG(TRANSPOSE(u)), MATMUL(base(:, :, 1), u))
         CALL dense_spectral_occ_projector_from_chc(chc_p, occ_scale, occ_mat_p)
         CALL dense_spectral_occ_projector_from_chc(chc_m, occ_scale, occ_mat_m)
         active_projector_weight = dense_pair_projector_weight(projector_weight, &
                                                               occ_curv(var_i(s)), &
                                                               occ_curv(var_j(s)), curv_max)
         reference(:, s) = (1.0_dp - projector_weight)*(diag_p(:, 1) - diag_m(:, 1))/(2.0_dp*fd_step) + &
                           projector_weight*(eig_p(:, 1) - eig_m(:, 1))/(2.0_dp*fd_step)
         DO i = 1, n
            reference(i, s) = reference(i, s) - active_projector_weight* &
                              0.5_dp*REAL(occ_mat_p(i, i) - occ_mat_m(i, i), KIND=dp)/ &
                              (2.0_dp*fd_step*occ_curv(i))
         END DO
      END DO

      err = MAXVAL(ABS(physical_response(:, :) - reference(:, :)))
      active_err = MAXVAL(ABS(physical_response(:, :) - basis_response(:, :)))

      IF (err > tol .OR. active_err < active_tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,2(1X,ES12.5))') TRIM(label), "failed with error/active", err, active_err
      ELSE
         WRITE (*, '(A,1X,A,2(1X,ES12.5))') TRIM(label), "OK with error/active", err, active_err
      END IF

   END SUBROUTINE check_complex_projector_physical_response_channel

! **************************************************************************************************
!> \brief Check that the physical Projector-Hessian keeps off-diagonal occupation response.
!>        The diagonal density response is held fixed here because the production
!>        Schur model consumes it in the auxiliary-energy response channel.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_complex_projector_physical_hessian_split(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 3, nrot = n*(n - 1)
      REAL(KIND=dp), PARAMETER                           :: active_tol = 1.0E-5_dp, &
                                                            fd_step = 1.0E-6_dp, &
                                                            projector_weight = 0.63_dp

      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: base, chc, chc_m, chc_p, dchc, docc_mat, &
                                                            occ0, occ_m, occ_p, u
      INTEGER                                            :: i, ivar, j, r, s
      INTEGER, DIMENSION(nrot)                           :: var_i, var_j, var_kind
      REAL(KIND=dp)                                      :: active, err
      REAL(KIND=dp), DIMENSION(n)                        :: eps_diag, occ_scale
      REAL(KIND=dp), DIMENSION(n, n)                     :: grad_docc_im, grad_docc_re, grad_m_im, &
                                                            grad_m_re, grad_p_im, grad_p_re, x_im, &
                                                            x_m_im, x_m_re, x_p_im, x_p_re, x_re
      REAL(KIND=dp), DIMENSION(nrot, nrot)               :: fixed_hessian, model_hessian, &
                                                            reference_hessian

      base(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      base(1, 1) = CMPLX(-0.021_dp, 0.0_dp, KIND=dp)
      base(2, 2) = CMPLX(0.014_dp, 0.0_dp, KIND=dp)
      base(3, 3) = CMPLX(0.41_dp, 0.0_dp, KIND=dp)
      base(1, 2) = CMPLX(0.026_dp, -0.031_dp, KIND=dp)
      base(1, 3) = CMPLX(-0.052_dp, 0.017_dp, KIND=dp)
      base(2, 3) = CMPLX(0.028_dp, 0.037_dp, KIND=dp)
      base(2, 1) = CONJG(base(1, 2))
      base(3, 1) = CONJG(base(1, 3))
      base(3, 2) = CONJG(base(2, 3))

      x_re(:, :) = 0.0_dp
      x_re(1, 2) = 0.051_dp
      x_re(1, 3) = -0.029_dp
      x_re(2, 3) = 0.034_dp
      DO j = 2, n
         DO i = 1, j - 1
            x_re(j, i) = -x_re(i, j)
         END DO
      END DO

      x_im(:, :) = 0.0_dp
      x_im(1, 1) = 0.016_dp
      x_im(2, 2) = -0.015_dp
      x_im(3, 3) = 0.023_dp
      x_im(1, 2) = -0.024_dp
      x_im(1, 3) = 0.021_dp
      x_im(2, 3) = 0.019_dp
      DO j = 2, n
         DO i = 1, j - 1
            x_im(j, i) = x_im(i, j)
         END DO
      END DO

      ivar = 0
      DO j = 2, n
         DO i = 1, j - 1
            ivar = ivar + 1
            var_i(ivar) = i
            var_j(ivar) = j
            var_kind(ivar) = 1
            ivar = ivar + 1
            var_i(ivar) = i
            var_j(ivar) = j
            var_kind(ivar) = 2
         END DO
      END DO

      occ_scale(:) = [1.46_dp, 0.81_dp, 0.28_dp]
      CALL dense_rotation_u(x_re, x_im, u)
      chc(:, :) = MATMUL(CONJG(TRANSPOSE(u)), MATMUL(base, u))
      CALL dense_spectral_occ_projector_from_chc(chc, occ_scale, occ0)
      occ0(:, :) = projector_weight*occ0(:, :)
      DO i = 1, n
         occ0(i, i) = occ0(i, i) + CMPLX((1.0_dp - projector_weight)*occ_scale(i), 0.0_dp, KIND=dp)
      END DO

      DO s = 1, nrot
         x_p_re(:, :) = x_re(:, :)
         x_m_re(:, :) = x_re(:, :)
         x_p_im(:, :) = x_im(:, :)
         x_m_im(:, :) = x_im(:, :)
         CALL perturb_complex_rotation_var_2d(x_p_re, x_p_im, var_i(s), var_j(s), var_kind(s), fd_step)
         CALL perturb_complex_rotation_var_2d(x_m_re, x_m_im, var_i(s), var_j(s), var_kind(s), -fd_step)

         CALL dense_projector_rotation_gradient(x_p_re, x_p_im, base, occ0, grad_p_re, grad_p_im, eps_diag)
         CALL dense_projector_rotation_gradient(x_m_re, x_m_im, base, occ0, grad_m_re, grad_m_im, eps_diag)
         DO r = 1, nrot
            SELECT CASE (var_kind(r))
            CASE (1)
               fixed_hessian(r, s) = &
                  2.0_dp*(grad_p_re(var_i(r), var_j(r)) - grad_m_re(var_i(r), var_j(r)))/(2.0_dp*fd_step)
            CASE (2)
               fixed_hessian(r, s) = &
                  2.0_dp*(grad_p_im(var_i(r), var_j(r)) - grad_m_im(var_i(r), var_j(r)))/(2.0_dp*fd_step)
            END SELECT
         END DO

         CALL dense_rotation_u(x_p_re, x_p_im, u)
         chc_p(:, :) = MATMUL(CONJG(TRANSPOSE(u)), MATMUL(base, u))
         CALL dense_rotation_u(x_m_re, x_m_im, u)
         chc_m(:, :) = MATMUL(CONJG(TRANSPOSE(u)), MATMUL(base, u))

         CALL dense_spectral_occ_projector_from_chc(chc_p, occ_scale, occ_p)
         CALL dense_spectral_occ_projector_from_chc(chc_m, occ_scale, occ_m)
         occ_p(:, :) = projector_weight*occ_p(:, :)
         occ_m(:, :) = projector_weight*occ_m(:, :)
         DO i = 1, n
            occ_p(i, i) = occ0(i, i)
            occ_m(i, i) = occ0(i, i)
         END DO
         CALL dense_projector_rotation_gradient(x_p_re, x_p_im, base, occ_p, grad_p_re, grad_p_im, eps_diag)
         CALL dense_projector_rotation_gradient(x_m_re, x_m_im, base, occ_m, grad_m_re, grad_m_im, eps_diag)

         dchc(:, :) = (chc_p(:, :) - chc_m(:, :))/(2.0_dp*fd_step)
         CALL dense_spectral_occ_projector_derivative_chc(chc, dchc, occ_scale, docc_mat)
         docc_mat(:, :) = projector_weight*docc_mat(:, :)
         DO i = 1, n
            docc_mat(i, i) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
         END DO
         CALL dense_projector_rotation_gradient(x_re, x_im, base, docc_mat, &
                                                grad_docc_re, grad_docc_im, eps_diag)

         DO r = 1, nrot
            SELECT CASE (var_kind(r))
            CASE (1)
               reference_hessian(r, s) = &
                  2.0_dp*(grad_p_re(var_i(r), var_j(r)) - grad_m_re(var_i(r), var_j(r)))/(2.0_dp*fd_step)
               model_hessian(r, s) = fixed_hessian(r, s) + &
                                     2.0_dp*grad_docc_re(var_i(r), var_j(r))
            CASE (2)
               reference_hessian(r, s) = &
                  2.0_dp*(grad_p_im(var_i(r), var_j(r)) - grad_m_im(var_i(r), var_j(r)))/(2.0_dp*fd_step)
               model_hessian(r, s) = fixed_hessian(r, s) + &
                                     2.0_dp*grad_docc_im(var_i(r), var_j(r))
            END SELECT
         END DO
      END DO

      err = MAXVAL(ABS(model_hessian(:, :) - reference_hessian(:, :)))
      active = MAXVAL(ABS(model_hessian(:, :) - fixed_hessian(:, :)))
      IF (err > tol .OR. active < active_tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,2(1X,ES12.5))') TRIM(label), "failed with error/active", err, active
      ELSE
         WRITE (*, '(A,1X,A,2(1X,ES12.5))') TRIM(label), "OK with error/active", err, active
      END IF

   END SUBROUTINE check_complex_projector_physical_hessian_split

! **************************************************************************************************
!> \brief Dense derivative of the blended projector rotation gradient with respect to occupations.
!> \param x_re ...
!> \param x_im ...
!> \param base ...
!> \param projector_weight ...
!> \param var_i ...
!> \param var_j ...
!> \param var_kind ...
!> \param response ...
! **************************************************************************************************
   SUBROUTINE dense_blended_projector_response(x_re, x_im, base, projector_weight, &
                                               var_i, var_j, var_kind, response)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: x_re, x_im
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: base
      REAL(KIND=dp), INTENT(IN)                          :: projector_weight
      INTEGER, DIMENSION(:), INTENT(IN)                  :: var_i, var_j, var_kind
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: response

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: chc, evec, occ_mat, u
      INTEGER                                            :: band, i, j, n, r
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eps_diag, evals
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: grad_im, grad_re

      n = SIZE(x_re, 1)
      ALLOCATE (chc(n, n), evec(n, n), occ_mat(n, n), u(n, n), evals(n))
      ALLOCATE (grad_re(n, n), grad_im(n, n), eps_diag(n))

      CALL dense_rotation_u(x_re, x_im, u)
      chc(:, :) = MATMUL(CONJG(TRANSPOSE(u)), MATMUL(base, u))
      IF (projector_weight > 0.0_dp) CALL diag_complex(chc, evec, evals)

      DO band = 1, n
         occ_mat(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
         IF (projector_weight > 0.0_dp) THEN
            DO j = 1, n
               DO i = 1, n
                  occ_mat(i, j) = projector_weight*evec(i, band)*CONJG(evec(j, band))
               END DO
            END DO
         END IF
         occ_mat(band, band) = occ_mat(band, band) + &
                               CMPLX(1.0_dp - projector_weight, 0.0_dp, KIND=dp)
         CALL dense_projector_rotation_gradient(x_re, x_im, base, occ_mat, grad_re, grad_im, eps_diag)
         DO r = 1, SIZE(var_i)
            SELECT CASE (var_kind(r))
            CASE (1)
               response(band, r) = 2.0_dp*grad_re(var_i(r), var_j(r))
            CASE (2)
               response(band, r) = 2.0_dp*grad_im(var_i(r), var_j(r))
            END SELECT
         END DO
      END DO

      DEALLOCATE (chc, evec, occ_mat, u, evals, grad_re, grad_im, eps_diag)

   END SUBROUTINE dense_blended_projector_response

! **************************************************************************************************
!> \brief Dense finite-complex physical occupation response from Rayleigh/eigenvalue shifts.
!> \param x_re ...
!> \param x_im ...
!> \param base ...
!> \param projector_weight ...
!> \param occ_scale ...
!> \param occ_curv ...
!> \param var_i ...
!> \param var_j ...
!> \param var_kind ...
!> \param response ...
! **************************************************************************************************
   SUBROUTINE dense_blended_projector_physical_response(x_re, x_im, base, projector_weight, occ_scale, &
                                                        occ_curv, var_i, var_j, var_kind, response)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: x_re, x_im
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: base
      REAL(KIND=dp), INTENT(IN)                          :: projector_weight
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: occ_scale, occ_curv
      INTEGER, DIMENSION(:), INTENT(IN)                  :: var_i, var_j, var_kind
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: response

      REAL(KIND=dp), PARAMETER                           :: fd_step = 1.0E-6_dp

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: chc_m, chc_p, evec, occ_mat_m, &
                                                            occ_mat_p, u
      INTEGER                                            :: band, n, s
      REAL(KIND=dp)                                      :: active_projector_weight, curv_max, &
                                                            projector_response, raw_response
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: evals_m, evals_p
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: x_m_im, x_m_re, x_p_im, x_p_re

      n = SIZE(x_re, 1)
      curv_max = MAXVAL(occ_curv)
      ALLOCATE (chc_p(n, n), chc_m(n, n), evec(n, n), occ_mat_p(n, n), occ_mat_m(n, n), u(n, n), &
                evals_p(n), evals_m(n))
      ALLOCATE (x_p_re(n, n), x_p_im(n, n), x_m_re(n, n), x_m_im(n, n))

      DO s = 1, SIZE(var_i)
         x_p_re(:, :) = x_re(:, :)
         x_p_im(:, :) = x_im(:, :)
         x_m_re(:, :) = x_re(:, :)
         x_m_im(:, :) = x_im(:, :)
         CALL perturb_complex_rotation_var_2d(x_p_re, x_p_im, var_i(s), var_j(s), var_kind(s), fd_step)
         CALL perturb_complex_rotation_var_2d(x_m_re, x_m_im, var_i(s), var_j(s), var_kind(s), -fd_step)

         CALL dense_rotation_u(x_p_re, x_p_im, u)
         chc_p(:, :) = MATMUL(CONJG(TRANSPOSE(u)), MATMUL(base, u))
         CALL dense_rotation_u(x_m_re, x_m_im, u)
         chc_m(:, :) = MATMUL(CONJG(TRANSPOSE(u)), MATMUL(base, u))
         IF (projector_weight > 0.0_dp) THEN
            CALL diag_complex(chc_p, evec, evals_p)
            CALL diag_complex(chc_m, evec, evals_m)
            CALL dense_spectral_occ_projector_from_chc(chc_p, occ_scale, occ_mat_p)
            CALL dense_spectral_occ_projector_from_chc(chc_m, occ_scale, occ_mat_m)
         END IF
         active_projector_weight = dense_pair_projector_weight(projector_weight, occ_curv(var_i(s)), &
                                                               occ_curv(var_j(s)), curv_max)

         DO band = 1, n
            raw_response = (REAL(chc_p(band, band), KIND=dp) - &
                            REAL(chc_m(band, band), KIND=dp))/(2.0_dp*fd_step)
            projector_response = raw_response
            IF (projector_weight > 0.0_dp) &
               projector_response = (evals_p(band) - evals_m(band))/(2.0_dp*fd_step)
            response(band, s) = (1.0_dp - projector_weight)*raw_response + &
                                projector_weight*projector_response
            IF (active_projector_weight > 0.0_dp .AND. ABS(occ_curv(band)) > EPSILON(occ_curv(band))) &
               response(band, s) = response(band, s) - active_projector_weight* &
                                   0.5_dp*REAL(occ_mat_p(band, band) - occ_mat_m(band, band), KIND=dp)/ &
                                   (2.0_dp*fd_step*occ_curv(band))
         END DO
      END DO

      DEALLOCATE (chc_p, chc_m, evec, occ_mat_p, occ_mat_m, u, evals_p, evals_m, &
                  x_p_re, x_p_im, x_m_re, x_m_im)

   END SUBROUTINE dense_blended_projector_physical_response

! **************************************************************************************************
!> \brief Dense gradient for the same blended projector model with arbitrary occupation weights.
!> \param x_re ...
!> \param x_im ...
!> \param base ...
!> \param occ_scale ...
!> \param projector_weight ...
!> \param grad_re ...
!> \param grad_im ...
!> \param eps_diag ...
! **************************************************************************************************
   SUBROUTINE dense_blended_projector_rotation_gradient(x_re, x_im, base, occ_scale, &
                                                        projector_weight, grad_re, grad_im, eps_diag)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: x_re, x_im
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: base
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: occ_scale
      REAL(KIND=dp), INTENT(IN)                          :: projector_weight
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: grad_re, grad_im
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: eps_diag

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: occ_mat
      INTEGER                                            :: i, n

      n = SIZE(x_re, 1)
      ALLOCATE (occ_mat(n, n))
      IF (projector_weight > 0.0_dp) THEN
         CALL dense_spectral_occ_projector(x_re, x_im, base, occ_scale, occ_mat)
         occ_mat(:, :) = projector_weight*occ_mat(:, :)
      ELSE
         occ_mat(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      END IF
      DO i = 1, n
         occ_mat(i, i) = occ_mat(i, i) + &
                         CMPLX((1.0_dp - projector_weight)*occ_scale(i), 0.0_dp, KIND=dp)
      END DO
      CALL dense_projector_rotation_gradient(x_re, x_im, base, occ_mat, grad_re, grad_im, eps_diag)
      DEALLOCATE (occ_mat)

   END SUBROUTINE dense_blended_projector_rotation_gradient

! **************************************************************************************************
!> \brief Check the split Schur projector for global fixed-N Mermin response plus channel-local tether.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_mixed_fixed_n_tether_projector(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: nchannels = 2, ne = 5, nrot = 3

      INTEGER                                            :: first_energy, i, j, local_channel, nmo, &
                                                            r, s
      INTEGER, DIMENSION(nchannels)                      :: channel_nmo, e_offset
      REAL(KIND=dp)                                      :: err, occ_sum, projected_term, ref_term, &
                                                            shift_left, shift_right, weight_sum
      REAL(KIND=dp), DIMENSION(ne)                       :: occ_curv, tether_curv
      REAL(KIND=dp), DIMENSION(ne, ne)                   :: hproj
      REAL(KIND=dp), DIMENSION(ne, nrot)                 :: response
      REAL(KIND=dp), DIMENSION(nrot, nrot)               :: projected, reference

      channel_nmo = [3, 2]
      e_offset = [1, 4]
      occ_curv = [0.31_dp, 0.74_dp, 0.16_dp, 0.47_dp, 0.22_dp]
      tether_curv = [0.05_dp, 0.02_dp, 0.08_dp, 0.07_dp, 0.03_dp]
      response(:, 1) = [0.13_dp, -0.21_dp, 0.09_dp, 0.31_dp, -0.17_dp]
      response(:, 2) = [-0.05_dp, 0.18_dp, 0.24_dp, -0.11_dp, 0.07_dp]
      response(:, 3) = [0.22_dp, 0.04_dp, -0.19_dp, 0.15_dp, -0.08_dp]

      hproj(:, :) = 0.0_dp
      DO i = 1, ne
         hproj(i, i) = occ_curv(i) + tether_curv(i)
      END DO
      occ_sum = SUM(occ_curv)
      DO j = 1, ne
         DO i = 1, ne
            hproj(i, j) = hproj(i, j) - occ_curv(i)*occ_curv(j)/occ_sum
         END DO
      END DO
      DO local_channel = 1, nchannels
         nmo = channel_nmo(local_channel)
         first_energy = e_offset(local_channel)
         weight_sum = SUM(tether_curv(first_energy:first_energy + nmo - 1))
         DO j = first_energy, first_energy + nmo - 1
            DO i = first_energy, first_energy + nmo - 1
               hproj(i, j) = hproj(i, j) - tether_curv(i)*tether_curv(j)/weight_sum
            END DO
         END DO
      END DO

      DO s = 1, nrot
         DO r = 1, nrot
            ref_term = DOT_PRODUCT(response(:, r), MATMUL(hproj, response(:, s)))
            projected_term = DOT_PRODUCT((occ_curv(:) + tether_curv(:))*response(:, r), response(:, s))

            shift_left = DOT_PRODUCT(occ_curv, response(:, r))
            shift_right = DOT_PRODUCT(occ_curv, response(:, s))
            projected_term = projected_term - shift_left*shift_right/occ_sum

            DO local_channel = 1, nchannels
               nmo = channel_nmo(local_channel)
               first_energy = e_offset(local_channel)
               weight_sum = SUM(tether_curv(first_energy:first_energy + nmo - 1))
               shift_left = DOT_PRODUCT(tether_curv(first_energy:first_energy + nmo - 1), &
                                        response(first_energy:first_energy + nmo - 1, r))
               shift_right = DOT_PRODUCT(tether_curv(first_energy:first_energy + nmo - 1), &
                                         response(first_energy:first_energy + nmo - 1, s))
               projected_term = projected_term - shift_left*shift_right/weight_sum
            END DO

            reference(r, s) = ref_term
            projected(r, s) = projected_term
         END DO
      END DO

      err = MAXVAL(ABS(projected(:, :) - reference(:, :)))
      IF (err > tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,ES12.5)') TRIM(label), "failed with error", err
      ELSE
         WRITE (*, '(A,1X,A,1X,ES12.5)') TRIM(label), "OK with error", err
      END IF

   END SUBROUTINE check_mixed_fixed_n_tether_projector

! **************************************************************************************************
!> \brief Check split Schur projector with distinct occupation and tether response channels.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_mixed_fixed_n_tether_cross_projector(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: nchannels = 2, ne = 5, nrot = 3

      INTEGER                                            :: first_energy, i, j, local_channel, nmo, &
                                                            r, s
      INTEGER, DIMENSION(nchannels)                      :: channel_nmo, e_offset
      REAL(KIND=dp)                                      :: err, occ_sum, projected_term, ref_term, &
                                                            shift_left, shift_right, weight_sum
      REAL(KIND=dp), DIMENSION(ne)                       :: occ_curv, tether_curv
      REAL(KIND=dp), DIMENSION(ne, ne)                   :: hocc, htether
      REAL(KIND=dp), DIMENSION(ne, nrot)                 :: response_occ, response_tether
      REAL(KIND=dp), DIMENSION(nrot, nrot)               :: projected, reference

      channel_nmo = [3, 2]
      e_offset = [1, 4]
      occ_curv = [0.31_dp, 0.74_dp, 0.16_dp, 0.47_dp, 0.22_dp]
      tether_curv = [0.05_dp, 0.02_dp, 0.08_dp, 0.07_dp, 0.03_dp]
      response_occ(:, 1) = [0.13_dp, -0.21_dp, 0.09_dp, 0.31_dp, -0.17_dp]
      response_occ(:, 2) = [-0.05_dp, 0.18_dp, 0.24_dp, -0.11_dp, 0.07_dp]
      response_occ(:, 3) = [0.22_dp, 0.04_dp, -0.19_dp, 0.15_dp, -0.08_dp]
      response_tether(:, 1) = [0.06_dp, -0.09_dp, 0.18_dp, 0.27_dp, -0.14_dp]
      response_tether(:, 2) = [-0.12_dp, 0.21_dp, 0.05_dp, -0.08_dp, 0.16_dp]
      response_tether(:, 3) = [0.19_dp, -0.02_dp, -0.11_dp, 0.23_dp, -0.04_dp]

      hocc(:, :) = 0.0_dp
      htether(:, :) = 0.0_dp
      DO i = 1, ne
         hocc(i, i) = occ_curv(i)
         htether(i, i) = tether_curv(i)
      END DO
      occ_sum = SUM(occ_curv)
      DO j = 1, ne
         DO i = 1, ne
            hocc(i, j) = hocc(i, j) - occ_curv(i)*occ_curv(j)/occ_sum
         END DO
      END DO
      DO local_channel = 1, nchannels
         nmo = channel_nmo(local_channel)
         first_energy = e_offset(local_channel)
         weight_sum = SUM(tether_curv(first_energy:first_energy + nmo - 1))
         DO j = first_energy, first_energy + nmo - 1
            DO i = first_energy, first_energy + nmo - 1
               htether(i, j) = htether(i, j) - tether_curv(i)*tether_curv(j)/weight_sum
            END DO
         END DO
      END DO

      DO s = 1, nrot
         DO r = 1, nrot
            ref_term = DOT_PRODUCT(response_occ(:, r), MATMUL(hocc, response_occ(:, s))) + &
                       DOT_PRODUCT(response_tether(:, r), MATMUL(htether, response_tether(:, s)))

            projected_term = DOT_PRODUCT(occ_curv(:)*response_occ(:, r), response_occ(:, s))
            shift_left = DOT_PRODUCT(occ_curv, response_occ(:, r))
            shift_right = DOT_PRODUCT(occ_curv, response_occ(:, s))
            projected_term = projected_term - shift_left*shift_right/occ_sum

            projected_term = projected_term + DOT_PRODUCT(tether_curv(:)*response_tether(:, r), &
                                                          response_tether(:, s))
            DO local_channel = 1, nchannels
               nmo = channel_nmo(local_channel)
               first_energy = e_offset(local_channel)
               weight_sum = SUM(tether_curv(first_energy:first_energy + nmo - 1))
               shift_left = DOT_PRODUCT(tether_curv(first_energy:first_energy + nmo - 1), &
                                        response_tether(first_energy:first_energy + nmo - 1, r))
               shift_right = DOT_PRODUCT(tether_curv(first_energy:first_energy + nmo - 1), &
                                         response_tether(first_energy:first_energy + nmo - 1, s))
               projected_term = projected_term - shift_left*shift_right/weight_sum
            END DO

            reference(r, s) = ref_term
            projected(r, s) = projected_term
         END DO
      END DO

      err = MAXVAL(ABS(projected(:, :) - reference(:, :)))
      IF (err > tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,ES12.5)') TRIM(label), "failed with error", err
      ELSE
         WRITE (*, '(A,1X,A,1X,ES12.5)') TRIM(label), "OK with error", err
      END IF

   END SUBROUTINE check_mixed_fixed_n_tether_cross_projector

! **************************************************************************************************
!> \brief Check the symmetrized Schur term for distinct rotation/energy responses.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_mixed_fixed_n_asymmetric_response_schur(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: ne = 5, nrot = 3

      INTEGER                                            :: i, j, r, s
      REAL(KIND=dp)                                      :: active_err, err, occ_sum, shift, &
                                                            symmetry_err
      REAL(KIND=dp), DIMENSION(ne)                       :: occ_curv
      REAL(KIND=dp), DIMENSION(ne, ne)                   :: hocc
      REAL(KIND=dp), DIMENSION(ne, nrot)                 :: b_occ, b_raw, response_occ, &
                                                            response_occ_projected, response_raw, &
                                                            response_raw_projected
      REAL(KIND=dp), DIMENSION(nrot, nrot)               :: mixed, occ_only, reference

      occ_curv = [0.31_dp, 0.74_dp, 0.16_dp, 0.47_dp, 0.22_dp]
      response_occ(:, 1) = [0.13_dp, -0.21_dp, 0.09_dp, 0.31_dp, -0.17_dp]
      response_occ(:, 2) = [-0.05_dp, 0.18_dp, 0.24_dp, -0.11_dp, 0.07_dp]
      response_occ(:, 3) = [0.22_dp, 0.04_dp, -0.19_dp, 0.15_dp, -0.08_dp]
      response_raw(:, 1) = [0.10_dp, -0.16_dp, 0.14_dp, 0.22_dp, -0.11_dp]
      response_raw(:, 2) = [-0.01_dp, 0.15_dp, 0.20_dp, -0.03_dp, 0.02_dp]
      response_raw(:, 3) = [0.17_dp, 0.08_dp, -0.13_dp, 0.10_dp, -0.03_dp]

      hocc(:, :) = 0.0_dp
      DO i = 1, ne
         hocc(i, i) = occ_curv(i)
      END DO
      occ_sum = SUM(occ_curv)
      DO j = 1, ne
         DO i = 1, ne
            hocc(i, j) = hocc(i, j) - occ_curv(i)*occ_curv(j)/occ_sum
         END DO
      END DO

      response_occ_projected(:, :) = response_occ(:, :)
      response_raw_projected(:, :) = response_raw(:, :)
      DO r = 1, nrot
         shift = DOT_PRODUCT(occ_curv, response_occ_projected(:, r))/occ_sum
         response_occ_projected(:, r) = response_occ_projected(:, r) - shift
         shift = DOT_PRODUCT(occ_curv, response_raw_projected(:, r))/occ_sum
         response_raw_projected(:, r) = response_raw_projected(:, r) - shift
      END DO

      b_occ(:, :) = MATMUL(hocc, response_occ)
      b_raw(:, :) = MATMUL(hocc, response_raw)
      DO s = 1, nrot
         DO r = 1, nrot
            mixed(r, s) = 0.5_dp*(DOT_PRODUCT(b_occ(:, r), response_raw_projected(:, s)) + &
                                  DOT_PRODUCT(b_raw(:, r), response_occ_projected(:, s)))
            reference(r, s) = 0.5_dp*(DOT_PRODUCT(response_occ(:, r), MATMUL(hocc, response_raw(:, s))) + &
                                      DOT_PRODUCT(response_raw(:, r), MATMUL(hocc, response_occ(:, s))))
            occ_only(r, s) = DOT_PRODUCT(response_occ(:, r), MATMUL(hocc, response_occ(:, s)))
         END DO
      END DO

      symmetry_err = MAXVAL(ABS(mixed(:, :) - TRANSPOSE(mixed(:, :))))
      active_err = MAXVAL(ABS(mixed(:, :) - occ_only(:, :)))
      err = MAX(MAXVAL(ABS(mixed(:, :) - reference(:, :))), symmetry_err)

      IF (err > tol .OR. active_err < 1.0E-6_dp) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,2(1X,ES12.5))') TRIM(label), "failed with error/active", err, active_err
      ELSE
         WRITE (*, '(A,1X,A,2(1X,ES12.5))') TRIM(label), "OK with error/active", err, active_err
      END IF

   END SUBROUTINE check_mixed_fixed_n_asymmetric_response_schur

! **************************************************************************************************
!> \brief Check exact dense elimination for the mixed fixed-N/tether energy block.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_mixed_fixed_n_tether_exact_response(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: nchannels = 2, ne = 6, nrot = 3

      INTEGER                                            :: first_energy, i, info, j, local_channel, &
                                                            nmo, r
      INTEGER, DIMENSION(nchannels)                      :: channel_nmo, e_offset
      INTEGER, DIMENSION(ne)                             :: ipiv
      REAL(KIND=dp)                                      :: active_err, channel_occ_sum, curv_max, &
                                                            err, gauge_weight, ref_value, shift, &
                                                            test_value, weight_sum
      REAL(KIND=dp), DIMENSION(ne)                       :: bq, occ_curv, tether_curv, x
      REAL(KIND=dp), DIMENSION(ne, ne)                   :: h, h_ref
      REAL(KIND=dp), DIMENSION(ne, nrot)                 :: b, response_occ, response_tether, rhs
      REAL(KIND=dp), DIMENSION(nrot)                     :: q
      REAL(KIND=dp), DIMENSION(nrot, nrot)               :: exact

      channel_nmo = [3, 3]
      e_offset = [1, 4]
      occ_curv = [0.31_dp, 0.74_dp, 0.16_dp, 0.0_dp, 0.0_dp, 0.0_dp]
      tether_curv = [0.05_dp, 0.02_dp, 0.08_dp, 0.07_dp, 0.03_dp, 0.0_dp]
      response_occ(:, 1) = [0.13_dp, -0.21_dp, 0.09_dp, 0.31_dp, -0.17_dp, 0.05_dp]
      response_occ(:, 2) = [-0.05_dp, 0.18_dp, 0.24_dp, -0.11_dp, 0.07_dp, -0.09_dp]
      response_occ(:, 3) = [0.22_dp, 0.04_dp, -0.19_dp, 0.15_dp, -0.08_dp, 0.12_dp]
      response_tether(:, 1) = [0.06_dp, -0.09_dp, 0.18_dp, 0.27_dp, -0.14_dp, 0.11_dp]
      response_tether(:, 2) = [-0.12_dp, 0.21_dp, 0.05_dp, -0.08_dp, 0.16_dp, -0.03_dp]
      response_tether(:, 3) = [0.19_dp, -0.02_dp, -0.11_dp, 0.23_dp, -0.04_dp, 0.08_dp]
      q = [0.21_dp, -0.13_dp, 0.17_dp]

      h(:, :) = 0.0_dp
      b(:, :) = 0.0_dp
      DO i = 1, ne
         h(i, i) = occ_curv(i) + tether_curv(i)
      END DO
      weight_sum = SUM(occ_curv)
      DO j = 1, ne
         DO i = 1, ne
            h(i, j) = h(i, j) - occ_curv(i)*occ_curv(j)/weight_sum
         END DO
      END DO
      DO r = 1, nrot
         shift = DOT_PRODUCT(occ_curv, response_occ(:, r))/weight_sum
         b(:, r) = b(:, r) + occ_curv(:)*(response_occ(:, r) - shift)
      END DO
      DO local_channel = 1, nchannels
         nmo = channel_nmo(local_channel)
         first_energy = e_offset(local_channel)
         weight_sum = SUM(tether_curv(first_energy:first_energy + nmo - 1))
         DO j = first_energy, first_energy + nmo - 1
            DO i = first_energy, first_energy + nmo - 1
               h(i, j) = h(i, j) - tether_curv(i)*tether_curv(j)/weight_sum
            END DO
         END DO
         DO r = 1, nrot
            shift = DOT_PRODUCT(tether_curv(first_energy:first_energy + nmo - 1), &
                                response_tether(first_energy:first_energy + nmo - 1, r))/weight_sum
            b(first_energy:first_energy + nmo - 1, r) = b(first_energy:first_energy + nmo - 1, r) + &
                                                        tether_curv(first_energy:first_energy + nmo - 1)* &
                                                        (response_tether(first_energy:first_energy + nmo - 1, r) - shift)
         END DO
      END DO

      curv_max = MAXVAL(occ_curv + tether_curv)
      gauge_weight = curv_max/REAL(ne, KIND=dp)
      h(:, :) = h(:, :) + gauge_weight
      DO local_channel = 1, nchannels
         nmo = channel_nmo(local_channel)
         first_energy = e_offset(local_channel)
         weight_sum = SUM(tether_curv(first_energy:first_energy + nmo - 1))
         channel_occ_sum = SUM(occ_curv(first_energy:first_energy + nmo - 1))
         IF (weight_sum <= EPSILON(weight_sum) .OR. channel_occ_sum > EPSILON(channel_occ_sum)) CYCLE
         DO j = first_energy, first_energy + nmo - 1
            DO i = first_energy, first_energy + nmo - 1
               h(i, j) = h(i, j) + gauge_weight
            END DO
         END DO
      END DO
      DO i = 1, ne
         IF (occ_curv(i) + tether_curv(i) > EPSILON(curv_max)) CYCLE
         h(i, :) = 0.0_dp
         h(:, i) = 0.0_dp
         h(i, i) = 1.0_dp
         b(i, :) = 0.0_dp
      END DO
      rhs(:, :) = b(:, :)
      h_ref(:, :) = h(:, :)
      CALL dgesv(ne, nrot, h, ne, ipiv, rhs, ne, info)
      IF (info /= 0) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,I6)') TRIM(label), "failed in dgesv", info
         RETURN
      END IF

      exact(:, :) = MATMUL(TRANSPOSE(b), rhs)
      active_err = ABS(exact(1, 1) - DOT_PRODUCT(response_occ(:, 1), occ_curv(:)*response_occ(:, 1)))
      bq(:) = MATMUL(b, q)
      x(:) = MATMUL(rhs, q)
      ref_value = -0.5_dp*DOT_PRODUCT(q, MATMUL(exact, q))
      test_value = 0.5_dp*DOT_PRODUCT(x, MATMUL(h_ref, x)) - DOT_PRODUCT(x, bq)
      err = MAX(MAXVAL(ABS(exact(:, :) - TRANSPOSE(exact(:, :)))), ABS(ref_value - test_value))

      IF (err > tol .OR. active_err < 1.0E-4_dp) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,2(1X,ES12.5))') TRIM(label), "failed with error/active", err, active_err
      ELSE
         WRITE (*, '(A,1X,A,2(1X,ES12.5))') TRIM(label), "OK with error/active", err, active_err
      END IF

   END SUBROUTINE check_mixed_fixed_n_tether_exact_response

! **************************************************************************************************
!> \brief Check that exact dense energy response is returned in the mixed response gauge.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_mixed_fixed_n_tether_exact_projection(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: nchannels = 2, ne = 6, nrot = 3

      INTEGER                                            :: info, r
      INTEGER, DIMENSION(nchannels)                      :: channel_nmo, e_offset
      INTEGER, DIMENSION(ne)                             :: ipiv
      LOGICAL                                            :: projected, projected_ok
      REAL(KIND=dp)                                      :: err, projected_occ_residual, &
                                                            projected_tether_residual, &
                                                            raw_occ_residual, raw_tether_residual
      REAL(KIND=dp), DIMENSION(ne)                       :: occ_curv, projected_response, &
                                                            raw_response, tether_curv
      REAL(KIND=dp), DIMENSION(ne, ne)                   :: h
      REAL(KIND=dp), DIMENSION(ne, nrot)                 :: b, response_occ, response_tether, rhs

      channel_nmo = [3, 3]
      e_offset = [1, 4]
      occ_curv = [0.31_dp, 0.74_dp, 0.16_dp, 0.0_dp, 0.0_dp, 0.0_dp]
      tether_curv = [0.05_dp, 0.02_dp, 0.08_dp, 0.07_dp, 0.03_dp, 0.0_dp]
      response_occ(:, 1) = [0.13_dp, -0.21_dp, 0.09_dp, 0.31_dp, -0.17_dp, 0.05_dp]
      response_occ(:, 2) = [-0.05_dp, 0.18_dp, 0.24_dp, -0.11_dp, 0.07_dp, -0.09_dp]
      response_occ(:, 3) = [0.22_dp, 0.04_dp, -0.19_dp, 0.15_dp, -0.08_dp, 0.12_dp]
      response_tether(:, 1) = [0.06_dp, -0.09_dp, 0.18_dp, 0.27_dp, -0.14_dp, 0.11_dp]
      response_tether(:, 2) = [-0.12_dp, 0.21_dp, 0.05_dp, -0.08_dp, 0.16_dp, -0.03_dp]
      response_tether(:, 3) = [0.19_dp, -0.02_dp, -0.11_dp, 0.23_dp, -0.04_dp, 0.08_dp]

      CALL build_test_mixed_energy_block(channel_nmo, e_offset, response_occ, response_tether, h, b)
      rhs(:, :) = b(:, :)
      CALL dgesv(ne, nrot, h, ne, ipiv, rhs, ne, info)
      IF (info /= 0) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,I6)') TRIM(label), "failed in dgesv", info
         RETURN
      END IF

      projected_ok = .TRUE.
      raw_occ_residual = 0.0_dp
      raw_tether_residual = 0.0_dp
      projected_occ_residual = 0.0_dp
      projected_tether_residual = 0.0_dp
      DO r = 1, nrot
         raw_response(:) = rhs(:, r)
         projected_response(:) = raw_response(:)
         CALL dense_project_response_metric_distributed(projected_response, occ_curv, tether_curv, &
                                                        channel_nmo, e_offset, projected)
         projected_ok = projected_ok .AND. projected
         raw_occ_residual = MAX(raw_occ_residual, ABS(DOT_PRODUCT(occ_curv, raw_response)))
         raw_tether_residual = MAX(raw_tether_residual, &
                                   dense_response_tether_residual(raw_response, tether_curv, channel_nmo, e_offset))
         projected_occ_residual = MAX(projected_occ_residual, &
                                      ABS(DOT_PRODUCT(occ_curv, projected_response)))
         projected_tether_residual = MAX(projected_tether_residual, &
                                         dense_response_tether_residual(projected_response, tether_curv, &
                                                                        channel_nmo, e_offset))
      END DO

      err = MAX(projected_occ_residual, projected_tether_residual)
      IF (err > tol .OR. .NOT. projected_ok .OR. raw_occ_residual < 1.0E-3_dp .OR. &
          raw_tether_residual < 1.0E-4_dp) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "failed with raw/projected occ/tether", raw_occ_residual, raw_tether_residual, &
            projected_occ_residual, projected_tether_residual
      ELSE
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "OK with raw/projected occ/tether", raw_occ_residual, raw_tether_residual, &
            projected_occ_residual, projected_tether_residual
      END IF

   END SUBROUTINE check_mixed_fixed_n_tether_exact_projection

! **************************************************************************************************
!> \brief Check exact energy-block elimination against the full augmented Schur system.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_mixed_fixed_n_tether_exact_schur(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: nchannels = 2, ne = 6, nrot = 3, &
                                                            ntot = ne + nrot

      INTEGER                                            :: info
      INTEGER, DIMENSION(MAX(ne, ntot))                  :: ipiv
      INTEGER, DIMENSION(nchannels)                      :: channel_nmo, e_offset
      REAL(KIND=dp)                                      :: err
      REAL(KIND=dp), DIMENSION(ne)                       :: grad_e, z_aug, z_red
      REAL(KIND=dp), DIMENSION(ne, ne)                   :: h
      REAL(KIND=dp), DIMENSION(ne, nrot)                 :: b, response_occ, response_tether, &
                                                            rhs_energy
      REAL(KIND=dp), DIMENSION(nrot)                     :: grad_rot, q_aug, q_red
      REAL(KIND=dp), DIMENSION(nrot, nrot)               :: a, schur
      REAL(KIND=dp), DIMENSION(ntot)                     :: rhs_aug
      REAL(KIND=dp), DIMENSION(ntot, ntot)               :: aug

      channel_nmo = [3, 3]
      e_offset = [1, 4]
      grad_e = [0.11_dp, -0.07_dp, 0.04_dp, -0.09_dp, 0.13_dp, 0.03_dp]
      grad_rot = [0.17_dp, -0.19_dp, 0.08_dp]
      response_occ(:, 1) = [0.13_dp, -0.21_dp, 0.09_dp, 0.31_dp, -0.17_dp, 0.05_dp]
      response_occ(:, 2) = [-0.05_dp, 0.18_dp, 0.24_dp, -0.11_dp, 0.07_dp, -0.09_dp]
      response_occ(:, 3) = [0.22_dp, 0.04_dp, -0.19_dp, 0.15_dp, -0.08_dp, 0.12_dp]
      response_tether(:, 1) = [0.06_dp, -0.09_dp, 0.18_dp, 0.27_dp, -0.14_dp, 0.11_dp]
      response_tether(:, 2) = [-0.12_dp, 0.21_dp, 0.05_dp, -0.08_dp, 0.16_dp, -0.03_dp]
      response_tether(:, 3) = [0.19_dp, -0.02_dp, -0.11_dp, 0.23_dp, -0.04_dp, 0.08_dp]

      CALL build_test_mixed_energy_block(channel_nmo, e_offset, response_occ, response_tether, h, b)

      a(:, :) = RESHAPE([1.37_dp, 0.04_dp, -0.03_dp, &
                         0.04_dp, 1.11_dp, 0.05_dp, &
                         -0.03_dp, 0.05_dp, 0.93_dp], SHAPE(a))
      rhs_energy(:, 1:nrot) = b(:, 1:nrot)
      CALL dgesv(ne, nrot, h, ne, ipiv(1:ne), rhs_energy, ne, info)
      IF (info /= 0) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,I6)') TRIM(label), "failed energy dgesv", info
         RETURN
      END IF

      schur(:, :) = a(:, :) - MATMUL(TRANSPOSE(b), rhs_energy)
      z_red(:) = grad_e(:)
      CALL build_test_mixed_energy_block(channel_nmo, e_offset, response_occ, response_tether, h, b)
      CALL dgesv(ne, 1, h, ne, ipiv(1:ne), z_red, ne, info)
      IF (info /= 0) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,I6)') TRIM(label), "failed gradient dgesv", info
         RETURN
      END IF
      q_red(:) = grad_rot(:) + MATMUL(TRANSPOSE(rhs_energy), grad_e)
      CALL dgesv(nrot, 1, schur, nrot, ipiv(1:nrot), q_red, nrot, info)
      IF (info /= 0) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,I6)') TRIM(label), "failed Schur dgesv", info
         RETURN
      END IF
      z_red(:) = z_red(:) + MATMUL(rhs_energy, q_red)

      CALL build_test_mixed_energy_block(channel_nmo, e_offset, response_occ, response_tether, h, b)
      aug(:, :) = 0.0_dp
      aug(1:ne, 1:ne) = h(:, :)
      aug(1:ne, ne + 1:ntot) = -b(:, :)
      aug(ne + 1:ntot, 1:ne) = -TRANSPOSE(b)
      aug(ne + 1:ntot, ne + 1:ntot) = a(:, :)
      rhs_aug(1:ne) = grad_e(:)
      rhs_aug(ne + 1:ntot) = grad_rot(:)
      CALL dgesv(ntot, 1, aug, ntot, ipiv(1:ntot), rhs_aug, ntot, info)
      IF (info /= 0) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,I6)') TRIM(label), "failed augmented dgesv", info
         RETURN
      END IF
      z_aug(:) = rhs_aug(1:ne)
      q_aug(:) = rhs_aug(ne + 1:ntot)

      err = MAX(MAXVAL(ABS(z_red - z_aug)), MAXVAL(ABS(q_red - q_aug)))

      IF (err > tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,ES12.5)') TRIM(label), "failed with error", err
      ELSE
         WRITE (*, '(A,1X,A,1X,ES12.5)') TRIM(label), "OK with error", err
      END IF

   END SUBROUTINE check_mixed_fixed_n_tether_exact_schur

! **************************************************************************************************
!> \brief Check positive-subspace rescue for an indefinite Schur rotation block.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_indefinite_schur_spectral_fallback(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 3

      INTEGER                                            :: i, info
      INTEGER, DIMENSION(n)                              :: ipiv
      REAL(KIND=dp)                                      :: direct_dot, err, fallback_dot, floor, &
                                                            mode_rhs, negative_projection, quad
      REAL(KIND=dp), DIMENSION(n)                        :: eval, rhs, sol_direct, sol_spectral
      REAL(KIND=dp), DIMENSION(n, n)                     :: schur, schur_fact, schur_modes

      schur(:, :) = 0.0_dp
      schur(1, 1) = 2.0_dp
      schur(2, 2) = -0.5_dp
      schur(3, 3) = 1.0_dp
      rhs(:) = [1.0_dp, 1.0_dp, 0.5_dp]

      schur_fact(:, :) = schur(:, :)
      sol_direct(:) = rhs(:)
      CALL dgesv(n, 1, schur_fact, n, ipiv, sol_direct, n, info)
      IF (info /= 0) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,I6)') TRIM(label), "failed direct dgesv", info
         RETURN
      END IF
      direct_dot = DOT_PRODUCT(rhs, sol_direct)

      schur_modes(:, :) = schur(:, :)
      CALL diamat_all(schur_modes, eval)
      floor = SQRT(EPSILON(1.0_dp))*MAX(1.0_dp, MAXVAL(ABS(eval)))
      sol_spectral(:) = 0.0_dp
      negative_projection = 0.0_dp
      DO i = 1, n
         mode_rhs = DOT_PRODUCT(schur_modes(:, i), rhs)
         IF (eval(i) > floor) THEN
            sol_spectral(:) = sol_spectral(:) + (mode_rhs/eval(i))*schur_modes(:, i)
         ELSE IF (eval(i) < 0.0_dp) THEN
            negative_projection = negative_projection + mode_rhs*mode_rhs
         END IF
      END DO
      fallback_dot = DOT_PRODUCT(rhs, sol_spectral)
      quad = DOT_PRODUCT(sol_spectral, MATMUL(schur, sol_spectral))
      negative_projection = 0.0_dp
      DO i = 1, n
         IF (eval(i) < 0.0_dp) &
            negative_projection = negative_projection + DOT_PRODUCT(schur_modes(:, i), sol_spectral)**2
      END DO

      err = MAX(ABS(fallback_dot - quad), negative_projection)
      IF (err > tol .OR. direct_dot >= -tol .OR. fallback_dot <= tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "failed with error/direct/fallback/neg", err, direct_dot, fallback_dot, negative_projection
      ELSE
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "OK with error/direct/fallback/neg", err, direct_dot, fallback_dot, negative_projection
      END IF

   END SUBROUTINE check_indefinite_schur_spectral_fallback

! **************************************************************************************************
!> \brief Check that indefinite Schur model norms use modal absolute curvature.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_indefinite_schur_absolute_metric(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 4

      INTEGER                                            :: i, neg_idx
      REAL(KIND=dp)                                      :: abs_metric, err, negative_abs_metric, &
                                                            negative_clamped_metric, &
                                                            negative_signed_metric, signed_metric
      REAL(KIND=dp), DIMENSION(n)                        :: coeff, eval, negative_step, step
      REAL(KIND=dp), DIMENSION(n, n)                     :: q, schur, schur_modes

      schur(:, :) = 0.0_dp
      schur(1, 1) = 2.0_dp
      schur(2, 2) = -1.8_dp
      schur(3, 3) = 0.25_dp
      schur(4, 4) = -0.22_dp
      q(:, :) = 0.0_dp
      DO i = 1, n
         q(i, i) = 1.0_dp
      END DO
      q(1, 1) = COS(0.31_dp)
      q(2, 2) = q(1, 1)
      q(1, 2) = -SIN(0.31_dp)
      q(2, 1) = -q(1, 2)
      q(3, 3) = COS(-0.27_dp)
      q(4, 4) = q(3, 3)
      q(3, 4) = -SIN(-0.27_dp)
      q(4, 3) = -q(3, 4)
      schur(:, :) = MATMUL(q, MATMUL(schur, TRANSPOSE(q)))
      step(:) = MATMUL(q, [0.70_dp, 0.72_dp, -0.31_dp, -0.33_dp])

      signed_metric = DOT_PRODUCT(step, MATMUL(schur, step))
      schur_modes(:, :) = schur(:, :)
      CALL diamat_all(schur_modes, eval)
      coeff(:) = MATMUL(TRANSPOSE(schur_modes), step)
      abs_metric = 0.0_dp
      DO i = 1, n
         abs_metric = abs_metric + ABS(eval(i))*coeff(i)*coeff(i)
      END DO
      neg_idx = MINLOC(eval, DIM=1)
      negative_step(:) = schur_modes(:, neg_idx)
      negative_signed_metric = DOT_PRODUCT(negative_step, MATMUL(schur, negative_step))
      negative_clamped_metric = MAX(0.0_dp, negative_signed_metric)
      negative_abs_metric = ABS(eval(neg_idx))

      err = ABS(abs_metric - SUM(ABS(eval(:))*coeff(:)*coeff(:)))
      IF (err > tol .OR. ABS(signed_metric) > 0.30_dp*abs_metric .OR. &
          abs_metric <= 1.0_dp .OR. negative_signed_metric >= -tol .OR. &
          negative_clamped_metric > tol .OR. negative_abs_metric <= tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,6(1X,ES12.5))') TRIM(label), &
            "failed with error/signed/abs/nsign/nclamp/nabs", err, signed_metric, abs_metric, &
            negative_signed_metric, negative_clamped_metric, negative_abs_metric
      ELSE
         WRITE (*, '(A,1X,A,6(1X,ES12.5))') TRIM(label), &
            "OK with error/signed/abs/nsign/nclamp/nabs", err, signed_metric, abs_metric, &
            negative_signed_metric, negative_clamped_metric, negative_abs_metric
      END IF

   END SUBROUTINE check_indefinite_schur_absolute_metric

! **************************************************************************************************
!> \brief Build the dense mixed fixed-N/tether energy block used by exact-response tests.
!> \param channel_nmo ...
!> \param e_offset ...
!> \param response_occ ...
!> \param response_tether ...
!> \param h ...
!> \param b ...
! **************************************************************************************************
   SUBROUTINE build_test_mixed_energy_block(channel_nmo, e_offset, response_occ, response_tether, h, b)
      INTEGER, DIMENSION(:), INTENT(IN)                  :: channel_nmo, e_offset
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: response_occ, response_tether
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: h, b

      INTEGER                                            :: first_energy, i, j, local_channel, &
                                                            nchannels, ne, nmo, nrot, r
      REAL(KIND=dp)                                      :: channel_occ_sum, curv_max, gauge_weight, &
                                                            shift, weight_sum
      REAL(KIND=dp), DIMENSION(SIZE(h, 1))               :: occ_curv, tether_curv

      ne = SIZE(h, 1)
      nrot = SIZE(b, 2)
      nchannels = SIZE(channel_nmo)
      occ_curv = [0.31_dp, 0.74_dp, 0.16_dp, 0.0_dp, 0.0_dp, 0.0_dp]
      tether_curv = [0.05_dp, 0.02_dp, 0.08_dp, 0.07_dp, 0.03_dp, 0.0_dp]
      h(:, :) = 0.0_dp
      b(:, :) = 0.0_dp
      DO i = 1, ne
         h(i, i) = occ_curv(i) + tether_curv(i)
      END DO
      weight_sum = SUM(occ_curv)
      DO j = 1, ne
         DO i = 1, ne
            h(i, j) = h(i, j) - occ_curv(i)*occ_curv(j)/weight_sum
         END DO
      END DO
      DO r = 1, nrot
         shift = DOT_PRODUCT(occ_curv, response_occ(:, r))/weight_sum
         b(:, r) = b(:, r) + occ_curv(:)*(response_occ(:, r) - shift)
      END DO
      DO local_channel = 1, nchannels
         nmo = channel_nmo(local_channel)
         first_energy = e_offset(local_channel)
         weight_sum = SUM(tether_curv(first_energy:first_energy + nmo - 1))
         DO j = first_energy, first_energy + nmo - 1
            DO i = first_energy, first_energy + nmo - 1
               h(i, j) = h(i, j) - tether_curv(i)*tether_curv(j)/weight_sum
            END DO
         END DO
         DO r = 1, nrot
            shift = DOT_PRODUCT(tether_curv(first_energy:first_energy + nmo - 1), &
                                response_tether(first_energy:first_energy + nmo - 1, r))/weight_sum
            b(first_energy:first_energy + nmo - 1, r) = b(first_energy:first_energy + nmo - 1, r) + &
                                                        tether_curv(first_energy:first_energy + nmo - 1)* &
                                                        (response_tether(first_energy:first_energy + nmo - 1, r) - shift)
         END DO
      END DO

      curv_max = MAXVAL(occ_curv + tether_curv)
      gauge_weight = curv_max/REAL(ne, KIND=dp)
      h(:, :) = h(:, :) + gauge_weight
      DO local_channel = 1, nchannels
         nmo = channel_nmo(local_channel)
         first_energy = e_offset(local_channel)
         weight_sum = SUM(tether_curv(first_energy:first_energy + nmo - 1))
         channel_occ_sum = SUM(occ_curv(first_energy:first_energy + nmo - 1))
         IF (weight_sum <= EPSILON(weight_sum) .OR. channel_occ_sum > EPSILON(channel_occ_sum)) CYCLE
         DO j = first_energy, first_energy + nmo - 1
            DO i = first_energy, first_energy + nmo - 1
               h(i, j) = h(i, j) + gauge_weight
            END DO
         END DO
      END DO
      DO i = 1, ne
         IF (occ_curv(i) + tether_curv(i) > EPSILON(curv_max)) CYCLE
         h(i, :) = 0.0_dp
         h(:, i) = 0.0_dp
         h(i, i) = 1.0_dp
         b(i, :) = 0.0_dp
      END DO

   END SUBROUTINE build_test_mixed_energy_block

! **************************************************************************************************
!> \brief Diagonal of U^H H U for each dense complex k-point channel.
!> \param x_re ...
!> \param x_im ...
!> \param base ...
!> \param hdiag ...
! **************************************************************************************************
   SUBROUTINE dense_complex_kpoint_hdiag(x_re, x_im, base, hdiag)
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: x_re, x_im
      COMPLEX(KIND=dp), DIMENSION(:, :, :), INTENT(IN)   :: base
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: hdiag

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: chc, u
      INTEGER                                            :: i, ik, n, nkp

      n = SIZE(x_re, 1)
      nkp = SIZE(x_re, 3)
      ALLOCATE (chc(n, n), u(n, n))
      DO ik = 1, nkp
         CALL dense_rotation_u(x_re(:, :, ik), x_im(:, :, ik), u)
         chc(:, :) = MATMUL(CONJG(TRANSPOSE(u)), MATMUL(base(:, :, ik), u))
         DO i = 1, n
            hdiag(i, ik) = REAL(chc(i, i), KIND=dp)
         END DO
      END DO
      DEALLOCATE (chc, u)

   END SUBROUTINE dense_complex_kpoint_hdiag

! **************************************************************************************************
!> \brief Eigenvalues of U^H H U for each dense complex k-point channel.
!> \param x_re ...
!> \param x_im ...
!> \param base ...
!> \param eigvals ...
! **************************************************************************************************
   SUBROUTINE dense_complex_kpoint_eigvals(x_re, x_im, base, eigvals)
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: x_re, x_im
      COMPLEX(KIND=dp), DIMENSION(:, :, :), INTENT(IN)   :: base
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: eigvals

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: chc, evec, u
      INTEGER                                            :: ik, n, nkp

      n = SIZE(x_re, 1)
      nkp = SIZE(x_re, 3)
      ALLOCATE (chc(n, n), evec(n, n), u(n, n))
      DO ik = 1, nkp
         CALL dense_rotation_u(x_re(:, :, ik), x_im(:, :, ik), u)
         chc(:, :) = MATMUL(CONJG(TRANSPOSE(u)), MATMUL(base(:, :, ik), u))
         CALL diag_complex(chc, evec, eigvals(:, ik))
      END DO
      DEALLOCATE (chc, evec, u)

   END SUBROUTINE dense_complex_kpoint_eigvals

! **************************************************************************************************
!> \brief Perturb one off-diagonal real or imaginary anti-Hermitian rotation variable.
!> \param x_re ...
!> \param x_im ...
!> \param i ...
!> \param j ...
!> \param kind ...
!> \param step ...
! **************************************************************************************************
   SUBROUTINE perturb_complex_rotation_var(x_re, x_im, i, j, kind, step)
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: x_re, x_im
      INTEGER, INTENT(IN)                                :: i, j, kind
      REAL(KIND=dp), INTENT(IN)                          :: step

      SELECT CASE (kind)
      CASE (1)
         x_re(i, j, 1) = x_re(i, j, 1) + step
         x_re(j, i, 1) = x_re(j, i, 1) - step
      CASE (2)
         x_im(i, j, 1) = x_im(i, j, 1) + step
         x_im(j, i, 1) = x_im(j, i, 1) + step
      END SELECT

   END SUBROUTINE perturb_complex_rotation_var

! **************************************************************************************************
!> \brief Perturb one packed complex rotation variable in a dense 2-D fixture.
!> \param x_re ...
!> \param x_im ...
!> \param i ...
!> \param j ...
!> \param kind ...
!> \param step ...
! **************************************************************************************************
   SUBROUTINE perturb_complex_rotation_var_2d(x_re, x_im, i, j, kind, step)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: x_re, x_im
      INTEGER, INTENT(IN)                                :: i, j, kind
      REAL(KIND=dp), INTENT(IN)                          :: step

      SELECT CASE (kind)
      CASE (1)
         x_re(i, j) = x_re(i, j) + step
         x_re(j, i) = x_re(j, i) - step
      CASE (2)
         x_im(i, j) = x_im(i, j) + step
         x_im(j, i) = x_im(j, i) + step
      END SELECT

   END SUBROUTINE perturb_complex_rotation_var_2d

! **************************************************************************************************
!> \brief Pack off-diagonal real/imaginary rotation gradients into a vector.
!> \param grad_re ...
!> \param grad_im ...
!> \param var_i ...
!> \param var_j ...
!> \param var_kind ...
!> \param grad_vec ...
! **************************************************************************************************
   SUBROUTINE extract_rotation_gradient_vector(grad_re, grad_im, var_i, var_j, var_kind, grad_vec)
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: grad_re, grad_im
      INTEGER, DIMENSION(:), INTENT(IN)                  :: var_i, var_j, var_kind
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: grad_vec

      INTEGER                                            :: r

      DO r = 1, SIZE(var_i)
         SELECT CASE (var_kind(r))
         CASE (1)
            grad_vec(r) = grad_re(var_i(r), var_j(r), 1)
         CASE (2)
            grad_vec(r) = grad_im(var_i(r), var_j(r), 1)
         END SELECT
      END DO

   END SUBROUTINE extract_rotation_gradient_vector

! **************************************************************************************************
!> \brief Pack off-diagonal real/imaginary rotation gradients from a dense 2-D fixture.
!> \param grad_re ...
!> \param grad_im ...
!> \param var_i ...
!> \param var_j ...
!> \param var_kind ...
!> \param grad_vec ...
! **************************************************************************************************
   SUBROUTINE extract_rotation_gradient_vector_2d(grad_re, grad_im, var_i, var_j, var_kind, grad_vec)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: grad_re, grad_im
      INTEGER, DIMENSION(:), INTENT(IN)                  :: var_i, var_j, var_kind
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: grad_vec

      INTEGER                                            :: r

      DO r = 1, SIZE(var_i)
         SELECT CASE (var_kind(r))
         CASE (1)
            grad_vec(r) = grad_re(var_i(r), var_j(r))
         CASE (2)
            grad_vec(r) = grad_im(var_i(r), var_j(r))
         END SELECT
      END DO

   END SUBROUTINE extract_rotation_gradient_vector_2d

! **************************************************************************************************
!> \brief Dense complex k-point Mermin value with occupations from auxiliary energies.
!> \param x_re ...
!> \param x_im ...
!> \param base ...
!> \param eps ...
!> \param wk ...
!> \param nel ...
!> \param sigma ...
!> \param maxocc ...
!> \return ...
! **************************************************************************************************
   FUNCTION dense_complex_kpoint_mermin_value(x_re, x_im, base, eps, wk, nel, sigma, maxocc) &
      RESULT(value)
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: x_re, x_im
      COMPLEX(KIND=dp), DIMENSION(:, :, :), INTENT(IN)   :: base
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: eps
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: wk
      REAL(KIND=dp), INTENT(IN)                          :: nel, sigma, maxocc
      REAL(KIND=dp)                                      :: value

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: chc, u
      INTEGER                                            :: i, ik, n, nkp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: h

      n = SIZE(eps, 1)
      nkp = SIZE(eps, 2)
      ALLOCATE (chc(n, n), u(n, n), h(n, nkp))
      DO ik = 1, nkp
         CALL dense_rotation_u(x_re(:, :, ik), x_im(:, :, ik), u)
         chc(:, :) = MATMUL(CONJG(TRANSPOSE(u)), MATMUL(base(:, :, ik), u))
         DO i = 1, n
            h(i, ik) = REAL(chc(i, i), KIND=dp)
         END DO
      END DO
      value = dense_kpoint_mermin_value(eps, h, wk, nel, sigma, maxocc)
      DEALLOCATE (chc, u, h)

   END FUNCTION dense_complex_kpoint_mermin_value

! **************************************************************************************************
!> \brief Dense coupled Mermin gradients for complex k-point rotations and auxiliary energies.
!> \param x_re ...
!> \param x_im ...
!> \param base ...
!> \param eps ...
!> \param wk ...
!> \param nel ...
!> \param sigma ...
!> \param maxocc ...
!> \param grad_re ...
!> \param grad_im ...
!> \param grad_eps ...
!> \param value ...
! **************************************************************************************************
   SUBROUTINE dense_complex_kpoint_mermin_gradient(x_re, x_im, base, eps, wk, nel, sigma, maxocc, &
                                                   grad_re, grad_im, grad_eps, value)
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: x_re, x_im
      COMPLEX(KIND=dp), DIMENSION(:, :, :), INTENT(IN)   :: base
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: eps
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: wk
      REAL(KIND=dp), INTENT(IN)                          :: nel, sigma, maxocc
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(OUT)     :: grad_re, grad_im
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: grad_eps
      REAL(KIND=dp), INTENT(OUT)                         :: value

      INTEGER                                            :: ik, n, nkp
      REAL(KIND=dp)                                      :: kTS, mu
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: occ_scale
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: h, occ

      n = SIZE(eps, 1)
      nkp = SIZE(eps, 2)
      ALLOCATE (h(n, nkp), occ(n, nkp), occ_scale(n))
      CALL dense_fermi_kpoint(occ, mu, kTS, eps, nel, wk, sigma, maxocc)
      DO ik = 1, nkp
         occ_scale(:) = 2.0_dp*wk(ik)*occ(:, ik)
         CALL dense_ref_rotation_gradient(x_re(:, :, ik), x_im(:, :, ik), base(:, :, ik), occ_scale, &
                                          grad_re(:, :, ik), grad_im(:, :, ik), h(:, ik))
      END DO
      CALL dense_kpoint_mermin_gradient(eps, h, wk, nel, sigma, maxocc, grad_eps)
      value = dense_kpoint_mermin_value(eps, h, wk, nel, sigma, maxocc)
      DEALLOCATE (h, occ, occ_scale)

   END SUBROUTINE dense_complex_kpoint_mermin_gradient

! **************************************************************************************************
!> \brief Check a non-diagonal occupation projector in the complex REF rotation gradient.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_complex_projector_rotation_gradient(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 3
      REAL(KIND=dp), PARAMETER                           :: fd_step = 1.0E-6_dp

      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: base, occ_mat
      INTEGER                                            :: i, j
      REAL(KIND=dp)                                      :: err, fd, value_m, value_p
      REAL(KIND=dp), DIMENSION(n)                        :: eps_diag
      REAL(KIND=dp), DIMENSION(n, n)                     :: grad_im, grad_re, x_im, x_m_im, x_m_re, &
                                                            x_p_im, x_p_re, x_re

      base(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      base(1, 1) = CMPLX(-0.31_dp, 0.0_dp, KIND=dp)
      base(2, 2) = CMPLX(-0.04_dp, 0.0_dp, KIND=dp)
      base(3, 3) = CMPLX(0.16_dp, 0.0_dp, KIND=dp)
      base(1, 2) = CMPLX(0.05_dp, -0.02_dp, KIND=dp)
      base(1, 3) = CMPLX(-0.03_dp, 0.04_dp, KIND=dp)
      base(2, 3) = CMPLX(0.02_dp, 0.03_dp, KIND=dp)
      base(2, 1) = CONJG(base(1, 2))
      base(3, 1) = CONJG(base(1, 3))
      base(3, 2) = CONJG(base(2, 3))

      occ_mat(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      occ_mat(1, 1) = CMPLX(1.60_dp, 0.0_dp, KIND=dp)
      occ_mat(2, 2) = CMPLX(0.72_dp, 0.0_dp, KIND=dp)
      occ_mat(3, 3) = CMPLX(0.18_dp, 0.0_dp, KIND=dp)
      occ_mat(1, 2) = CMPLX(0.08_dp, -0.03_dp, KIND=dp)
      occ_mat(1, 3) = CMPLX(-0.04_dp, 0.02_dp, KIND=dp)
      occ_mat(2, 3) = CMPLX(0.03_dp, 0.05_dp, KIND=dp)
      occ_mat(2, 1) = CONJG(occ_mat(1, 2))
      occ_mat(3, 1) = CONJG(occ_mat(1, 3))
      occ_mat(3, 2) = CONJG(occ_mat(2, 3))

      x_re(:, :) = 0.0_dp
      x_re(1, 2) = 0.06_dp
      x_re(1, 3) = -0.03_dp
      x_re(2, 3) = 0.04_dp
      DO j = 2, n
         DO i = 1, j - 1
            x_re(j, i) = -x_re(i, j)
         END DO
      END DO

      x_im(:, :) = 0.0_dp
      x_im(1, 1) = 0.02_dp
      x_im(2, 2) = -0.04_dp
      x_im(3, 3) = 0.01_dp
      x_im(1, 2) = -0.02_dp
      x_im(1, 3) = 0.01_dp
      x_im(2, 3) = 0.03_dp
      DO j = 2, n
         DO i = 1, j - 1
            x_im(j, i) = x_im(i, j)
         END DO
      END DO

      CALL dense_projector_rotation_gradient(x_re, x_im, base, occ_mat, grad_re, grad_im, eps_diag)

      err = 0.0_dp
      DO j = 2, n
         DO i = 1, j - 1
            x_p_re(:, :) = x_re(:, :)
            x_m_re(:, :) = x_re(:, :)
            x_p_im(:, :) = x_im(:, :)
            x_m_im(:, :) = x_im(:, :)
            x_p_re(i, j) = x_p_re(i, j) + fd_step
            x_p_re(j, i) = x_p_re(j, i) - fd_step
            x_m_re(i, j) = x_m_re(i, j) - fd_step
            x_m_re(j, i) = x_m_re(j, i) + fd_step
            value_p = dense_projector_rotation_value(x_p_re, x_p_im, base, occ_mat)
            value_m = dense_projector_rotation_value(x_m_re, x_m_im, base, occ_mat)
            fd = (value_p - value_m)/(2.0_dp*fd_step)
            err = MAX(err, ABS(fd - grad_re(i, j)))

            x_p_re(:, :) = x_re(:, :)
            x_m_re(:, :) = x_re(:, :)
            x_p_im(:, :) = x_im(:, :)
            x_m_im(:, :) = x_im(:, :)
            x_p_im(i, j) = x_p_im(i, j) + fd_step
            x_p_im(j, i) = x_p_im(j, i) + fd_step
            x_m_im(i, j) = x_m_im(i, j) - fd_step
            x_m_im(j, i) = x_m_im(j, i) - fd_step
            value_p = dense_projector_rotation_value(x_p_re, x_p_im, base, occ_mat)
            value_m = dense_projector_rotation_value(x_m_re, x_m_im, base, occ_mat)
            fd = (value_p - value_m)/(2.0_dp*fd_step)
            err = MAX(err, ABS(fd - grad_im(i, j)))
         END DO
      END DO

      DO i = 1, n
         x_p_re(:, :) = x_re(:, :)
         x_m_re(:, :) = x_re(:, :)
         x_p_im(:, :) = x_im(:, :)
         x_m_im(:, :) = x_im(:, :)
         x_p_im(i, i) = x_p_im(i, i) + fd_step
         x_m_im(i, i) = x_m_im(i, i) - fd_step
         value_p = dense_projector_rotation_value(x_p_re, x_p_im, base, occ_mat)
         value_m = dense_projector_rotation_value(x_m_re, x_m_im, base, occ_mat)
         fd = (value_p - value_m)/(2.0_dp*fd_step)
         err = MAX(err, ABS(fd - 0.5_dp*grad_im(i, i)))
      END DO

      IF (err > tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,ES12.5)') TRIM(label), "failed with error", err
      ELSE
         WRITE (*, '(A,1X,A,1X,ES12.5)') TRIM(label), "OK with error", err
      END IF

   END SUBROUTINE check_complex_projector_rotation_gradient

! **************************************************************************************************
!> \brief Check that the finite rotation Hessian includes the moving occupation projector.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_complex_occ_projector_hessian(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 3, nrot = n*(n - 1)
      REAL(KIND=dp), PARAMETER                           :: active_tol = 1.0E-7_dp, &
                                                            fd_step = 2.0E-5_dp

      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: base, chc, chc_m, chc_p, dchc, &
                                                            docc_diag, docc_mat, docc_offdiag, &
                                                            occ_mat, u
      INTEGER                                            :: i, j, r, s
      INTEGER, DIMENSION(nrot)                           :: var_i, var_j, var_kind
      REAL(KIND=dp)                                      :: active_diag_err, active_err, &
                                                            active_offdiag_err, err, fd, &
                                                            split_err, symmetry_err, value_m, &
                                                            value_p
      REAL(KIND=dp), DIMENSION(n)                        :: eps_diag, occ_scale
      REAL(KIND=dp), DIMENSION(n, n)                     :: grad_im, grad_re, x_im, x_m_im, x_m_re, &
                                                            x_p_im, x_p_re, x_re
      REAL(KIND=dp), DIMENSION(nrot)                     :: grad_m_vec, grad_p_vec
      REAL(KIND=dp), DIMENSION(nrot, nrot)               :: hdiag, hfrechet, hfrozen, hmove, &
                                                            hoffdiag, hsplit

      base(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      base(1, 1) = CMPLX(-0.37_dp, 0.0_dp, KIND=dp)
      base(2, 2) = CMPLX(-0.08_dp, 0.0_dp, KIND=dp)
      base(3, 3) = CMPLX(0.21_dp, 0.0_dp, KIND=dp)
      base(1, 2) = CMPLX(0.07_dp, -0.03_dp, KIND=dp)
      base(1, 3) = CMPLX(-0.04_dp, 0.05_dp, KIND=dp)
      base(2, 3) = CMPLX(0.03_dp, 0.04_dp, KIND=dp)
      base(2, 1) = CONJG(base(1, 2))
      base(3, 1) = CONJG(base(1, 3))
      base(3, 2) = CONJG(base(2, 3))

      occ_scale(:) = [1.54_dp, 0.83_dp, 0.25_dp]

      x_re(:, :) = 0.0_dp
      x_re(1, 2) = 0.05_dp
      x_re(1, 3) = -0.025_dp
      x_re(2, 3) = 0.035_dp
      DO j = 2, n
         DO i = 1, j - 1
            x_re(j, i) = -x_re(i, j)
         END DO
      END DO

      x_im(:, :) = 0.0_dp
      x_im(1, 1) = 0.015_dp
      x_im(2, 2) = -0.03_dp
      x_im(3, 3) = 0.02_dp
      x_im(1, 2) = -0.015_dp
      x_im(1, 3) = 0.02_dp
      x_im(2, 3) = 0.025_dp
      DO j = 2, n
         DO i = 1, j - 1
            x_im(j, i) = x_im(i, j)
         END DO
      END DO

      r = 0
      DO j = 2, n
         DO i = 1, j - 1
            r = r + 1
            var_i(r) = i
            var_j(r) = j
            var_kind(r) = 1
            r = r + 1
            var_i(r) = i
            var_j(r) = j
            var_kind(r) = 2
         END DO
      END DO

      CALL dense_spectral_projector_gradient(x_re, x_im, base, occ_scale, grad_re, grad_im, eps_diag)

      err = 0.0_dp
      hdiag(:, :) = 0.0_dp
      hoffdiag(:, :) = 0.0_dp
      DO s = 1, nrot
         x_p_re(:, :) = x_re(:, :)
         x_m_re(:, :) = x_re(:, :)
         x_p_im(:, :) = x_im(:, :)
         x_m_im(:, :) = x_im(:, :)
         CALL perturb_complex_rotation_var_2d(x_p_re, x_p_im, var_i(s), var_j(s), var_kind(s), fd_step)
         CALL perturb_complex_rotation_var_2d(x_m_re, x_m_im, var_i(s), var_j(s), var_kind(s), -fd_step)
         value_p = dense_spectral_projector_value(x_p_re, x_p_im, base, occ_scale)
         value_m = dense_spectral_projector_value(x_m_re, x_m_im, base, occ_scale)
         fd = (value_p - value_m)/(2.0_dp*fd_step)
         SELECT CASE (var_kind(s))
         CASE (1)
            err = MAX(err, ABS(fd - grad_re(var_i(s), var_j(s))))
         CASE (2)
            err = MAX(err, ABS(fd - grad_im(var_i(s), var_j(s))))
         END SELECT
      END DO

      CALL dense_rotation_u(x_re, x_im, u)
      chc(:, :) = MATMUL(CONJG(TRANSPOSE(u)), MATMUL(base, u))
      CALL dense_spectral_occ_projector(x_re, x_im, base, occ_scale, occ_mat)
      DO s = 1, nrot
         x_p_re(:, :) = x_re(:, :)
         x_m_re(:, :) = x_re(:, :)
         x_p_im(:, :) = x_im(:, :)
         x_m_im(:, :) = x_im(:, :)
         CALL perturb_complex_rotation_var_2d(x_p_re, x_p_im, var_i(s), var_j(s), var_kind(s), fd_step)
         CALL perturb_complex_rotation_var_2d(x_m_re, x_m_im, var_i(s), var_j(s), var_kind(s), -fd_step)

         CALL dense_spectral_projector_gradient(x_p_re, x_p_im, base, occ_scale, &
                                                grad_re, grad_im, eps_diag)
         CALL extract_rotation_gradient_vector_2d(grad_re, grad_im, var_i, var_j, var_kind, grad_p_vec)
         CALL dense_spectral_projector_gradient(x_m_re, x_m_im, base, occ_scale, &
                                                grad_re, grad_im, eps_diag)
         CALL extract_rotation_gradient_vector_2d(grad_re, grad_im, var_i, var_j, var_kind, grad_m_vec)
         hmove(:, s) = (grad_p_vec(:) - grad_m_vec(:))/(2.0_dp*fd_step)

         CALL dense_projector_rotation_gradient(x_p_re, x_p_im, base, occ_mat, grad_re, grad_im, eps_diag)
         CALL extract_rotation_gradient_vector_2d(grad_re, grad_im, var_i, var_j, var_kind, grad_p_vec)
         CALL dense_projector_rotation_gradient(x_m_re, x_m_im, base, occ_mat, grad_re, grad_im, eps_diag)
         CALL extract_rotation_gradient_vector_2d(grad_re, grad_im, var_i, var_j, var_kind, grad_m_vec)
         hfrozen(:, s) = (grad_p_vec(:) - grad_m_vec(:))/(2.0_dp*fd_step)

         CALL dense_rotation_u(x_p_re, x_p_im, u)
         chc_p(:, :) = MATMUL(CONJG(TRANSPOSE(u)), MATMUL(base, u))
         CALL dense_rotation_u(x_m_re, x_m_im, u)
         chc_m(:, :) = MATMUL(CONJG(TRANSPOSE(u)), MATMUL(base, u))
         dchc(:, :) = (chc_p(:, :) - chc_m(:, :))/(2.0_dp*fd_step)
         CALL dense_spectral_occ_projector_derivative_chc(chc, dchc, occ_scale, docc_mat)
         CALL dense_projector_rotation_gradient(x_re, x_im, base, docc_mat, grad_re, grad_im, eps_diag)
         CALL extract_rotation_gradient_vector_2d(grad_re, grad_im, var_i, var_j, var_kind, grad_p_vec)
         hfrechet(:, s) = hfrozen(:, s) + grad_p_vec(:)
         docc_diag(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
         docc_offdiag(:, :) = docc_mat(:, :)
         DO i = 1, n
            docc_diag(i, i) = docc_mat(i, i)
            docc_offdiag(i, i) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
         END DO
         CALL dense_projector_rotation_gradient(x_re, x_im, base, docc_diag, grad_re, grad_im, eps_diag)
         CALL extract_rotation_gradient_vector_2d(grad_re, grad_im, var_i, var_j, var_kind, grad_p_vec)
         hdiag(:, s) = grad_p_vec(:)
         CALL dense_projector_rotation_gradient(x_re, x_im, base, docc_offdiag, grad_re, grad_im, eps_diag)
         CALL extract_rotation_gradient_vector_2d(grad_re, grad_im, var_i, var_j, var_kind, grad_p_vec)
         hoffdiag(:, s) = grad_p_vec(:)
      END DO

      hsplit(:, :) = hfrozen(:, :) + hoffdiag(:, :)
      symmetry_err = MAXVAL(ABS(hmove(:, :) - TRANSPOSE(hmove(:, :))))
      active_err = MAXVAL(ABS(hmove(:, :) - hfrozen(:, :)))
      active_diag_err = MAXVAL(ABS(hdiag(:, :)))
      active_offdiag_err = MAXVAL(ABS(hoffdiag(:, :)))
      split_err = MAXVAL(ABS(hmove(:, :) - (hsplit(:, :) + hdiag(:, :))))
      err = MAX(err, symmetry_err)
      err = MAX(err, MAXVAL(ABS(hmove(:, :) - hfrechet(:, :))))
      err = MAX(err, split_err)

      IF (err > tol .OR. active_err < active_tol .OR. active_diag_err < active_tol .OR. &
          active_offdiag_err < active_tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), "failed with error/active/diag/offdiag", &
            err, active_err, active_diag_err, active_offdiag_err
      ELSE
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), "OK with error/active/diag/offdiag", &
            err, active_err, active_diag_err, active_offdiag_err
      END IF

   END SUBROUTINE check_complex_occ_projector_hessian

! **************************************************************************************************
!> \brief Check the covariant curvature removed by a spectral projector blend.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_complex_projector_covariant_curvature(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 3, nrot = n*(n - 1)
      REAL(KIND=dp), PARAMETER                           :: active_tol = 1.0E-4_dp, &
                                                            fd_step = 2.0E-5_dp, &
                                                            projector_weight = 0.30_dp

      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: base, occ_mat
      INTEGER                                            :: i, j, r, s
      INTEGER, DIMENSION(nrot)                           :: var_i, var_j, var_kind
      REAL(KIND=dp)                                      :: active, err, shift_err
      REAL(KIND=dp), DIMENSION(n)                        :: eps_diag, occ_scale
      REAL(KIND=dp), DIMENSION(n, n)                     :: grad_im, grad_re, x_im, x_m_im, x_m_re, &
                                                            x_p_im, x_p_re, x_re
      REAL(KIND=dp), DIMENSION(nrot)                     :: grad_m_vec, grad_p_vec
      REAL(KIND=dp), DIMENSION(nrot, nrot)               :: hfrozen, hmoving

      base(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      base(1, 1) = CMPLX(-0.42_dp, 0.0_dp, KIND=dp)
      base(2, 2) = CMPLX(-0.03_dp, 0.0_dp, KIND=dp)
      base(3, 3) = CMPLX(0.26_dp, 0.0_dp, KIND=dp)
      occ_scale(:) = [1.62_dp, 0.78_dp, 0.19_dp]

      x_re(:, :) = 0.0_dp
      x_im(:, :) = 0.0_dp
      occ_mat(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      DO i = 1, n
         occ_mat(i, i) = CMPLX(occ_scale(i), 0.0_dp, KIND=dp)
      END DO

      r = 0
      DO j = 2, n
         DO i = 1, j - 1
            r = r + 1
            var_i(r) = i
            var_j(r) = j
            var_kind(r) = 1
            r = r + 1
            var_i(r) = i
            var_j(r) = j
            var_kind(r) = 2
         END DO
      END DO

      DO s = 1, nrot
         x_p_re(:, :) = x_re(:, :)
         x_m_re(:, :) = x_re(:, :)
         x_p_im(:, :) = x_im(:, :)
         x_m_im(:, :) = x_im(:, :)
         CALL perturb_complex_rotation_var_2d(x_p_re, x_p_im, var_i(s), var_j(s), var_kind(s), fd_step)
         CALL perturb_complex_rotation_var_2d(x_m_re, x_m_im, var_i(s), var_j(s), var_kind(s), -fd_step)

         CALL dense_projector_rotation_gradient(x_p_re, x_p_im, base, occ_mat, grad_re, grad_im, eps_diag)
         CALL extract_rotation_gradient_vector_2d(grad_re, grad_im, var_i, var_j, var_kind, grad_p_vec)
         CALL dense_projector_rotation_gradient(x_m_re, x_m_im, base, occ_mat, grad_re, grad_im, eps_diag)
         CALL extract_rotation_gradient_vector_2d(grad_re, grad_im, var_i, var_j, var_kind, grad_m_vec)
         hfrozen(:, s) = (grad_p_vec(:) - grad_m_vec(:))/(2.0_dp*fd_step)

         CALL dense_blended_projector_rotation_gradient(x_p_re, x_p_im, base, occ_scale, &
                                                        projector_weight, grad_re, grad_im, eps_diag)
         CALL extract_rotation_gradient_vector_2d(grad_re, grad_im, var_i, var_j, var_kind, grad_p_vec)
         CALL dense_blended_projector_rotation_gradient(x_m_re, x_m_im, base, occ_scale, &
                                                        projector_weight, grad_re, grad_im, eps_diag)
         CALL extract_rotation_gradient_vector_2d(grad_re, grad_im, var_i, var_j, var_kind, grad_m_vec)
         hmoving(:, s) = (grad_p_vec(:) - grad_m_vec(:))/(2.0_dp*fd_step)
      END DO

      err = MAXVAL(ABS(hmoving(:, :) - (1.0_dp - projector_weight)*hfrozen(:, :)))
      shift_err = MAXVAL(ABS((hmoving(:, :) - hfrozen(:, :)) + projector_weight*hfrozen(:, :)))
      active = MAXVAL(ABS(projector_weight*hfrozen(:, :)))

      IF (err > tol .OR. shift_err > tol .OR. active < active_tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,3(1X,ES12.5))') TRIM(label), "failed with error/shift/active", &
            err, shift_err, active
      ELSE
         WRITE (*, '(A,1X,A,3(1X,ES12.5))') TRIM(label), "OK with error/shift/active", &
            err, shift_err, active
      END IF

   END SUBROUTINE check_complex_projector_covariant_curvature

! **************************************************************************************************
!> \brief Check that the spectral projector covariance is removed in the diagonal REF gauge.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_complex_projector_covariant_gauge_split(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 3, nrot = n*(n - 1)
      REAL(KIND=dp), PARAMETER                           :: active_tol = 1.0E-5_dp, &
                                                            fd_step = 2.0E-5_dp, &
                                                            projector_weight = 0.35_dp

      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: base, occ_blend, occ_diag, occ_spectral
      INTEGER                                            :: i, j, r, s
      INTEGER, DIMENSION(nrot)                           :: var_i, var_j, var_kind
      REAL(KIND=dp)                                      :: err, gauge_err
      REAL(KIND=dp), DIMENSION(n)                        :: eps_diag, occ_scale
      REAL(KIND=dp), DIMENSION(n, n)                     :: grad_im, grad_re, x_im, x_m_im, x_m_re, &
                                                            x_p_im, x_p_re, x_re
      REAL(KIND=dp), DIMENSION(nrot)                     :: grad_m_vec, grad_p_vec
      REAL(KIND=dp), DIMENSION(nrot, nrot)               :: hblend_frozen, hdiag_frozen, hmoving

      base(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      base(1, 1) = CMPLX(-0.37_dp, 0.0_dp, KIND=dp)
      base(2, 2) = CMPLX(-0.06_dp, 0.0_dp, KIND=dp)
      base(3, 3) = CMPLX(0.24_dp, 0.0_dp, KIND=dp)
      base(1, 2) = CMPLX(0.045_dp, -0.030_dp, KIND=dp)
      base(1, 3) = CMPLX(-0.035_dp, 0.025_dp, KIND=dp)
      base(2, 3) = CMPLX(0.025_dp, 0.040_dp, KIND=dp)
      base(2, 1) = CONJG(base(1, 2))
      base(3, 1) = CONJG(base(1, 3))
      base(3, 2) = CONJG(base(2, 3))
      occ_scale(:) = [1.58_dp, 0.81_dp, 0.22_dp]

      x_re(:, :) = 0.0_dp
      x_re(1, 2) = 0.052_dp
      x_re(1, 3) = -0.021_dp
      x_re(2, 3) = 0.034_dp
      DO j = 2, n
         DO i = 1, j - 1
            x_re(j, i) = -x_re(i, j)
         END DO
      END DO

      x_im(:, :) = 0.0_dp
      x_im(1, 1) = 0.018_dp
      x_im(2, 2) = -0.026_dp
      x_im(3, 3) = 0.013_dp
      x_im(1, 2) = -0.019_dp
      x_im(1, 3) = 0.023_dp
      x_im(2, 3) = 0.016_dp
      DO j = 2, n
         DO i = 1, j - 1
            x_im(j, i) = x_im(i, j)
         END DO
      END DO

      r = 0
      DO j = 2, n
         DO i = 1, j - 1
            r = r + 1
            var_i(r) = i
            var_j(r) = j
            var_kind(r) = 1
            r = r + 1
            var_i(r) = i
            var_j(r) = j
            var_kind(r) = 2
         END DO
      END DO

      occ_diag(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      DO i = 1, n
         occ_diag(i, i) = CMPLX(occ_scale(i), 0.0_dp, KIND=dp)
      END DO
      CALL dense_spectral_occ_projector(x_re, x_im, base, occ_scale, occ_spectral)
      occ_blend(:, :) = projector_weight*occ_spectral(:, :) + (1.0_dp - projector_weight)*occ_diag(:, :)

      DO s = 1, nrot
         x_p_re(:, :) = x_re(:, :)
         x_m_re(:, :) = x_re(:, :)
         x_p_im(:, :) = x_im(:, :)
         x_m_im(:, :) = x_im(:, :)
         CALL perturb_complex_rotation_var_2d(x_p_re, x_p_im, var_i(s), var_j(s), var_kind(s), fd_step)
         CALL perturb_complex_rotation_var_2d(x_m_re, x_m_im, var_i(s), var_j(s), var_kind(s), -fd_step)

         CALL dense_projector_rotation_gradient(x_p_re, x_p_im, base, occ_diag, grad_re, grad_im, eps_diag)
         CALL extract_rotation_gradient_vector_2d(grad_re, grad_im, var_i, var_j, var_kind, grad_p_vec)
         CALL dense_projector_rotation_gradient(x_m_re, x_m_im, base, occ_diag, grad_re, grad_im, eps_diag)
         CALL extract_rotation_gradient_vector_2d(grad_re, grad_im, var_i, var_j, var_kind, grad_m_vec)
         hdiag_frozen(:, s) = (grad_p_vec(:) - grad_m_vec(:))/(2.0_dp*fd_step)

         CALL dense_projector_rotation_gradient(x_p_re, x_p_im, base, occ_blend, grad_re, grad_im, eps_diag)
         CALL extract_rotation_gradient_vector_2d(grad_re, grad_im, var_i, var_j, var_kind, grad_p_vec)
         CALL dense_projector_rotation_gradient(x_m_re, x_m_im, base, occ_blend, grad_re, grad_im, eps_diag)
         CALL extract_rotation_gradient_vector_2d(grad_re, grad_im, var_i, var_j, var_kind, grad_m_vec)
         hblend_frozen(:, s) = (grad_p_vec(:) - grad_m_vec(:))/(2.0_dp*fd_step)

         CALL dense_blended_projector_rotation_gradient(x_p_re, x_p_im, base, occ_scale, &
                                                        projector_weight, grad_re, grad_im, eps_diag)
         CALL extract_rotation_gradient_vector_2d(grad_re, grad_im, var_i, var_j, var_kind, grad_p_vec)
         CALL dense_blended_projector_rotation_gradient(x_m_re, x_m_im, base, occ_scale, &
                                                        projector_weight, grad_re, grad_im, eps_diag)
         CALL extract_rotation_gradient_vector_2d(grad_re, grad_im, var_i, var_j, var_kind, grad_m_vec)
         hmoving(:, s) = (grad_p_vec(:) - grad_m_vec(:))/(2.0_dp*fd_step)
      END DO

      err = MAXVAL(ABS(hmoving(:, :) - (1.0_dp - projector_weight)*hdiag_frozen(:, :)))
      gauge_err = MAXVAL(ABS(hmoving(:, :) - (1.0_dp - projector_weight)*hblend_frozen(:, :)))

      IF (err > tol .OR. gauge_err < active_tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,2(1X,ES12.5))') TRIM(label), "failed with error/gauge", &
            err, gauge_err
      ELSE
         WRITE (*, '(A,1X,A,2(1X,ES12.5))') TRIM(label), "OK with error/gauge", err, gauge_err
      END IF

   END SUBROUTINE check_complex_projector_covariant_gauge_split

! **************************************************************************************************
!> \brief Check that diagonal density response is not part of the physical rotation Hessian.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_complex_projector_density_hessian_split(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 3, nrot = n*(n - 1)
      REAL(KIND=dp), PARAMETER                           :: active_tol = 1.0E-5_dp, &
                                                            fd_step = 2.0E-5_dp, &
                                                            projector_weight = 0.35_dp

      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: base, chc, chc_m, chc_p, dchc, &
                                                            docc_diag, docc_mat, occ_diag, u
      INTEGER                                            :: i, j, r, s
      INTEGER, DIMENSION(nrot)                           :: var_i, var_j, var_kind
      REAL(KIND=dp)                                      :: density_active, err, wrong_err
      REAL(KIND=dp), DIMENSION(n)                        :: eps_diag, occ_scale
      REAL(KIND=dp), DIMENSION(n, n)                     :: grad_im, grad_re, x_im, x_m_im, x_m_re, &
                                                            x_p_im, x_p_re, x_re
      REAL(KIND=dp), DIMENSION(nrot)                     :: grad_m_vec, grad_p_vec
      REAL(KIND=dp), DIMENSION(nrot, nrot)               :: hdensity, hdiag_frozen, hmoving

      base(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      base(1, 1) = CMPLX(-0.39_dp, 0.0_dp, KIND=dp)
      base(2, 2) = CMPLX(-0.05_dp, 0.0_dp, KIND=dp)
      base(3, 3) = CMPLX(0.27_dp, 0.0_dp, KIND=dp)
      base(1, 2) = CMPLX(0.050_dp, -0.032_dp, KIND=dp)
      base(1, 3) = CMPLX(-0.037_dp, 0.024_dp, KIND=dp)
      base(2, 3) = CMPLX(0.028_dp, 0.043_dp, KIND=dp)
      base(2, 1) = CONJG(base(1, 2))
      base(3, 1) = CONJG(base(1, 3))
      base(3, 2) = CONJG(base(2, 3))
      occ_scale(:) = [1.60_dp, 0.79_dp, 0.21_dp]

      x_re(:, :) = 0.0_dp
      x_re(1, 2) = 0.055_dp
      x_re(1, 3) = -0.023_dp
      x_re(2, 3) = 0.032_dp
      DO j = 2, n
         DO i = 1, j - 1
            x_re(j, i) = -x_re(i, j)
         END DO
      END DO

      x_im(:, :) = 0.0_dp
      x_im(1, 1) = 0.017_dp
      x_im(2, 2) = -0.024_dp
      x_im(3, 3) = 0.014_dp
      x_im(1, 2) = -0.021_dp
      x_im(1, 3) = 0.022_dp
      x_im(2, 3) = 0.018_dp
      DO j = 2, n
         DO i = 1, j - 1
            x_im(j, i) = x_im(i, j)
         END DO
      END DO

      r = 0
      DO j = 2, n
         DO i = 1, j - 1
            r = r + 1
            var_i(r) = i
            var_j(r) = j
            var_kind(r) = 1
            r = r + 1
            var_i(r) = i
            var_j(r) = j
            var_kind(r) = 2
         END DO
      END DO

      occ_diag(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      DO i = 1, n
         occ_diag(i, i) = CMPLX(occ_scale(i), 0.0_dp, KIND=dp)
      END DO
      CALL dense_rotation_u(x_re, x_im, u)
      chc(:, :) = MATMUL(CONJG(TRANSPOSE(u)), MATMUL(base, u))

      DO s = 1, nrot
         x_p_re(:, :) = x_re(:, :)
         x_m_re(:, :) = x_re(:, :)
         x_p_im(:, :) = x_im(:, :)
         x_m_im(:, :) = x_im(:, :)
         CALL perturb_complex_rotation_var_2d(x_p_re, x_p_im, var_i(s), var_j(s), var_kind(s), fd_step)
         CALL perturb_complex_rotation_var_2d(x_m_re, x_m_im, var_i(s), var_j(s), var_kind(s), -fd_step)

         CALL dense_projector_rotation_gradient(x_p_re, x_p_im, base, occ_diag, grad_re, grad_im, eps_diag)
         CALL extract_rotation_gradient_vector_2d(grad_re, grad_im, var_i, var_j, var_kind, grad_p_vec)
         CALL dense_projector_rotation_gradient(x_m_re, x_m_im, base, occ_diag, grad_re, grad_im, eps_diag)
         CALL extract_rotation_gradient_vector_2d(grad_re, grad_im, var_i, var_j, var_kind, grad_m_vec)
         hdiag_frozen(:, s) = (grad_p_vec(:) - grad_m_vec(:))/(2.0_dp*fd_step)

         CALL dense_blended_projector_rotation_gradient(x_p_re, x_p_im, base, occ_scale, &
                                                        projector_weight, grad_re, grad_im, eps_diag)
         CALL extract_rotation_gradient_vector_2d(grad_re, grad_im, var_i, var_j, var_kind, grad_p_vec)
         CALL dense_blended_projector_rotation_gradient(x_m_re, x_m_im, base, occ_scale, &
                                                        projector_weight, grad_re, grad_im, eps_diag)
         CALL extract_rotation_gradient_vector_2d(grad_re, grad_im, var_i, var_j, var_kind, grad_m_vec)
         hmoving(:, s) = (grad_p_vec(:) - grad_m_vec(:))/(2.0_dp*fd_step)

         CALL dense_rotation_u(x_p_re, x_p_im, u)
         chc_p(:, :) = MATMUL(CONJG(TRANSPOSE(u)), MATMUL(base, u))
         CALL dense_rotation_u(x_m_re, x_m_im, u)
         chc_m(:, :) = MATMUL(CONJG(TRANSPOSE(u)), MATMUL(base, u))
         dchc(:, :) = (chc_p(:, :) - chc_m(:, :))/(2.0_dp*fd_step)
         CALL dense_spectral_occ_projector_derivative_chc(chc, dchc, occ_scale, docc_mat)
         docc_diag(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
         DO i = 1, n
            docc_diag(i, i) = projector_weight*docc_mat(i, i)
         END DO
         CALL dense_projector_rotation_gradient(x_re, x_im, base, docc_diag, grad_re, grad_im, eps_diag)
         CALL extract_rotation_gradient_vector_2d(grad_re, grad_im, var_i, var_j, var_kind, hdensity(:, s))
      END DO

      err = MAXVAL(ABS(hmoving(:, :) - (1.0_dp - projector_weight)*hdiag_frozen(:, :)))
      density_active = MAXVAL(ABS(hdensity(:, :)))
      wrong_err = MAXVAL(ABS(hmoving(:, :) - &
                             ((1.0_dp - projector_weight)*hdiag_frozen(:, :) + hdensity(:, :))))

      IF (err > tol .OR. density_active < active_tol .OR. wrong_err < active_tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,3(1X,ES12.5))') TRIM(label), "failed with error/density/wrong", &
            err, density_active, wrong_err
      ELSE
         WRITE (*, '(A,1X,A,3(1X,ES12.5))') TRIM(label), "OK with error/density/wrong", &
            err, density_active, wrong_err
      END IF

   END SUBROUTINE check_complex_projector_density_hessian_split

! **************************************************************************************************
!> \brief Check that the physical projector base keeps the moving spectral projector invariant.
!>        The spectral projector part of Tr[P(CHC) CHC] is unitary-invariant for fixed occupations.
!>        Therefore the no-response rotation Hessian is the finite moving-projector Hessian,
!>        equivalently the diagonal fixed-occupation Hessian scaled by the non-projector weight.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_complex_projector_physical_invariant_base(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 3, nrot = n*(n - 1)
      REAL(KIND=dp), PARAMETER                           :: active_tol = 1.0E-5_dp, &
                                                            fd_step = 2.0E-5_dp, &
                                                            projector_weight = 0.58_dp

      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: base, chc, occ_blended, occ_diag, u
      INTEGER                                            :: i, j, r, s
      INTEGER, DIMENSION(nrot)                           :: var_i, var_j, var_kind
      REAL(KIND=dp)                                      :: correct_err, legacy_gap
      REAL(KIND=dp), DIMENSION(n)                        :: eps_diag, occ_scale
      REAL(KIND=dp), DIMENSION(n, n)                     :: grad_im, grad_re, x_im, x_m_im, x_m_re, &
                                                            x_p_im, x_p_re, x_re
      REAL(KIND=dp), DIMENSION(nrot)                     :: grad_m_vec, grad_p_vec
      REAL(KIND=dp), DIMENSION(nrot, nrot)               :: hdiag_frozen, hfixed_blended, hmoving

      base(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      base(1, 1) = CMPLX(-0.34_dp, 0.0_dp, KIND=dp)
      base(2, 2) = CMPLX(0.03_dp, 0.0_dp, KIND=dp)
      base(3, 3) = CMPLX(0.29_dp, 0.0_dp, KIND=dp)
      base(1, 2) = CMPLX(0.047_dp, -0.029_dp, KIND=dp)
      base(1, 3) = CMPLX(-0.041_dp, 0.026_dp, KIND=dp)
      base(2, 3) = CMPLX(0.031_dp, 0.039_dp, KIND=dp)
      base(2, 1) = CONJG(base(1, 2))
      base(3, 1) = CONJG(base(1, 3))
      base(3, 2) = CONJG(base(2, 3))
      occ_scale(:) = [1.48_dp, 0.83_dp, 0.25_dp]

      x_re(:, :) = 0.0_dp
      x_re(1, 2) = 0.047_dp
      x_re(1, 3) = -0.031_dp
      x_re(2, 3) = 0.029_dp
      DO j = 2, n
         DO i = 1, j - 1
            x_re(j, i) = -x_re(i, j)
         END DO
      END DO

      x_im(:, :) = 0.0_dp
      x_im(1, 1) = 0.019_dp
      x_im(2, 2) = -0.017_dp
      x_im(3, 3) = 0.021_dp
      x_im(1, 2) = -0.023_dp
      x_im(1, 3) = 0.018_dp
      x_im(2, 3) = 0.022_dp
      DO j = 2, n
         DO i = 1, j - 1
            x_im(j, i) = x_im(i, j)
         END DO
      END DO

      r = 0
      DO j = 2, n
         DO i = 1, j - 1
            r = r + 1
            var_i(r) = i
            var_j(r) = j
            var_kind(r) = 1
            r = r + 1
            var_i(r) = i
            var_j(r) = j
            var_kind(r) = 2
         END DO
      END DO

      occ_diag(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      DO i = 1, n
         occ_diag(i, i) = CMPLX(occ_scale(i), 0.0_dp, KIND=dp)
      END DO
      CALL dense_rotation_u(x_re, x_im, u)
      chc(:, :) = MATMUL(CONJG(TRANSPOSE(u)), MATMUL(base, u))
      CALL dense_spectral_occ_projector_from_chc(chc, occ_scale, occ_blended)
      occ_blended(:, :) = projector_weight*occ_blended(:, :)
      DO i = 1, n
         occ_blended(i, i) = occ_blended(i, i) + &
                             CMPLX((1.0_dp - projector_weight)*occ_scale(i), 0.0_dp, KIND=dp)
      END DO

      DO s = 1, nrot
         x_p_re(:, :) = x_re(:, :)
         x_m_re(:, :) = x_re(:, :)
         x_p_im(:, :) = x_im(:, :)
         x_m_im(:, :) = x_im(:, :)
         CALL perturb_complex_rotation_var_2d(x_p_re, x_p_im, var_i(s), var_j(s), var_kind(s), fd_step)
         CALL perturb_complex_rotation_var_2d(x_m_re, x_m_im, var_i(s), var_j(s), var_kind(s), -fd_step)

         CALL dense_projector_rotation_gradient(x_p_re, x_p_im, base, occ_diag, grad_re, grad_im, eps_diag)
         CALL extract_rotation_gradient_vector_2d(grad_re, grad_im, var_i, var_j, var_kind, grad_p_vec)
         CALL dense_projector_rotation_gradient(x_m_re, x_m_im, base, occ_diag, grad_re, grad_im, eps_diag)
         CALL extract_rotation_gradient_vector_2d(grad_re, grad_im, var_i, var_j, var_kind, grad_m_vec)
         hdiag_frozen(:, s) = (grad_p_vec(:) - grad_m_vec(:))/(2.0_dp*fd_step)

         CALL dense_projector_rotation_gradient(x_p_re, x_p_im, base, occ_blended, grad_re, grad_im, eps_diag)
         CALL extract_rotation_gradient_vector_2d(grad_re, grad_im, var_i, var_j, var_kind, grad_p_vec)
         CALL dense_projector_rotation_gradient(x_m_re, x_m_im, base, occ_blended, grad_re, grad_im, eps_diag)
         CALL extract_rotation_gradient_vector_2d(grad_re, grad_im, var_i, var_j, var_kind, grad_m_vec)
         hfixed_blended(:, s) = (grad_p_vec(:) - grad_m_vec(:))/(2.0_dp*fd_step)

         CALL dense_blended_projector_rotation_gradient(x_p_re, x_p_im, base, occ_scale, &
                                                        projector_weight, grad_re, grad_im, eps_diag)
         CALL extract_rotation_gradient_vector_2d(grad_re, grad_im, var_i, var_j, var_kind, grad_p_vec)
         CALL dense_blended_projector_rotation_gradient(x_m_re, x_m_im, base, occ_scale, &
                                                        projector_weight, grad_re, grad_im, eps_diag)
         CALL extract_rotation_gradient_vector_2d(grad_re, grad_im, var_i, var_j, var_kind, grad_m_vec)
         hmoving(:, s) = (grad_p_vec(:) - grad_m_vec(:))/(2.0_dp*fd_step)
      END DO

      correct_err = MAXVAL(ABS(hmoving(:, :) - (1.0_dp - projector_weight)*hdiag_frozen(:, :)))
      legacy_gap = MAXVAL(ABS(hmoving(:, :) - (1.0_dp - projector_weight)*hfixed_blended(:, :)))

      IF (correct_err > tol .OR. legacy_gap < active_tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,2(1X,ES12.5))') TRIM(label), "failed with error/legacy_gap", &
            correct_err, legacy_gap
      ELSE
         WRITE (*, '(A,1X,A,2(1X,ES12.5))') TRIM(label), "OK with error/legacy_gap", &
            correct_err, legacy_gap
      END IF

   END SUBROUTINE check_complex_projector_physical_invariant_base

! **************************************************************************************************
!> \brief Separate projector-basis gauge noise from physical occupation response.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_complex_occ_projector_gauge_separation(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 3
      REAL(KIND=dp), PARAMETER                           :: fd_step = 1.0E-8_dp, gap = 2.0E-5_dp, &
                                                            slope = -0.65_dp

      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: chc, dchc, docc_mat, fd_mat, occ_mat_m, &
                                                            occ_mat_p
      REAL(KIND=dp)                                      :: err, gauge_err
      REAL(KIND=dp), DIMENSION(n)                        :: occ_equal, occ_linear

      chc(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      chc(1, 1) = CMPLX(0.00_dp, 0.0_dp, KIND=dp)
      chc(2, 2) = CMPLX(gap, 0.0_dp, KIND=dp)
      chc(3, 3) = CMPLX(0.41_dp, 0.0_dp, KIND=dp)

      dchc(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      dchc(1, 2) = CMPLX(0.31_dp, -0.17_dp, KIND=dp)
      dchc(2, 1) = CONJG(dchc(1, 2))

      occ_equal(:) = [0.74_dp, 0.74_dp, 0.18_dp]
      CALL dense_spectral_occ_projector_derivative_chc(chc, dchc, occ_equal, docc_mat)
      CALL dense_spectral_occ_projector_from_chc(chc + fd_step*dchc, occ_equal, occ_mat_p)
      CALL dense_spectral_occ_projector_from_chc(chc - fd_step*dchc, occ_equal, occ_mat_m)
      fd_mat(:, :) = (occ_mat_p(:, :) - occ_mat_m(:, :))/(2.0_dp*fd_step)
      gauge_err = MAX(MAXVAL(ABS(docc_mat(:, :))), MAXVAL(ABS(fd_mat(:, :))))

      occ_linear(:) = occ_equal(:)
      occ_linear(2) = occ_linear(1) + slope*gap
      CALL dense_spectral_occ_projector_derivative_chc(chc, dchc, occ_linear, docc_mat)
      CALL dense_spectral_occ_projector_from_chc(chc + fd_step*dchc, occ_linear, occ_mat_p)
      CALL dense_spectral_occ_projector_from_chc(chc - fd_step*dchc, occ_linear, occ_mat_m)
      fd_mat(:, :) = (occ_mat_p(:, :) - occ_mat_m(:, :))/(2.0_dp*fd_step)
      err = MAXVAL(ABS(docc_mat(:, :) - fd_mat(:, :)))

      IF (err > tol .OR. gauge_err > tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,2(1X,ES12.5))') TRIM(label), "failed with error/gauge", err, gauge_err
      ELSE
         WRITE (*, '(A,1X,A,2(1X,ES12.5))') TRIM(label), "OK with error/gauge", err, gauge_err
      END IF

   END SUBROUTINE check_complex_occ_projector_gauge_separation

! **************************************************************************************************
!> \brief Bound near-degenerate projector-basis response by the physical occupation curvature.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_complex_occ_projector_curvature_cap(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 3
      REAL(KIND=dp), PARAMETER                           :: fd_step = 1.0E-8_dp, gap = 2.0E-5_dp

      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: chc, dchc, docc_mat, fd_mat, occ_mat_m, &
                                                            occ_mat_p
      REAL(KIND=dp)                                      :: cap, capped_coeff, err, raw_coeff
      REAL(KIND=dp), DIMENSION(n)                        :: occ_curv, occ_equal, occ_noisy, &
                                                            occ_smooth

      chc(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      chc(1, 1) = CMPLX(0.00_dp, 0.0_dp, KIND=dp)
      chc(2, 2) = CMPLX(gap, 0.0_dp, KIND=dp)
      chc(3, 3) = CMPLX(0.39_dp, 0.0_dp, KIND=dp)

      dchc(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      dchc(1, 2) = CMPLX(0.23_dp, -0.14_dp, KIND=dp)
      dchc(2, 1) = CONJG(dchc(1, 2))

      occ_curv(:) = [0.031_dp, 0.027_dp, 0.004_dp]
      cap = 2.0_dp*MAX(occ_curv(1), occ_curv(2))

      occ_equal(:) = [0.75_dp, 0.75_dp, 0.12_dp]
      CALL dense_spectral_occ_projector_derivative_chc(chc, dchc, occ_equal, docc_mat, occ_curv)
      err = MAXVAL(ABS(docc_mat(:, :)))

      occ_smooth(:) = [0.75_dp, 0.75_dp - 0.25_dp*cap*gap, 0.12_dp]
      CALL dense_spectral_occ_projector_derivative_chc(chc, dchc, occ_smooth, docc_mat, occ_curv)
      CALL dense_spectral_occ_projector_from_chc(chc + fd_step*dchc, occ_smooth, occ_mat_p)
      CALL dense_spectral_occ_projector_from_chc(chc - fd_step*dchc, occ_smooth, occ_mat_m)
      fd_mat(:, :) = (occ_mat_p(:, :) - occ_mat_m(:, :))/(2.0_dp*fd_step)
      err = MAX(err, MAXVAL(ABS(docc_mat(:, :) - fd_mat(:, :))))

      occ_noisy(:) = [0.75_dp, 0.21_dp, 0.12_dp]
      CALL dense_spectral_occ_projector_derivative_chc(chc, dchc, occ_noisy, docc_mat, occ_curv)
      raw_coeff = (occ_noisy(1) - occ_noisy(2))/(REAL(chc(1, 1), KIND=dp) - &
                                                 REAL(chc(2, 2), KIND=dp))
      capped_coeff = REAL(docc_mat(1, 2), KIND=dp)/REAL(dchc(1, 2), KIND=dp)
      err = MAX(err, ABS(capped_coeff + cap))

      IF (err > tol .OR. ABS(raw_coeff) < 10.0_dp*cap) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), "failed with error/raw/cap/capped", &
            err, raw_coeff, cap, capped_coeff
      ELSE
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), "OK with error/raw/cap/capped", &
            err, raw_coeff, cap, capped_coeff
      END IF

   END SUBROUTINE check_complex_occ_projector_curvature_cap

! **************************************************************************************************
!> \brief Cap active/inert occupation-projector slopes by shared pair curvature.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_occ_projector_active_inert_pair_cap(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 3
      REAL(KIND=dp), PARAMETER                           :: gap = 6.0E-5_dp, gauge_gap = 3.0E-5_dp

      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: chc, dchc, docc_mat
      REAL(KIND=dp)                                      :: capped_coeff, err, max_cap, raw_coeff, &
                                                            shared_cap
      REAL(KIND=dp), DIMENSION(n)                        :: occ_curv, occ_noisy

      chc(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      chc(1, 1) = CMPLX(-0.31_dp, 0.0_dp, KIND=dp)
      chc(2, 2) = CMPLX(0.20_dp, 0.0_dp, KIND=dp)
      chc(3, 3) = CMPLX(0.20_dp + gap, 0.0_dp, KIND=dp)

      dchc(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      dchc(2, 3) = CMPLX(0.19_dp, -0.11_dp, KIND=dp)
      dchc(3, 2) = CONJG(dchc(2, 3))

      occ_curv(:) = [0.042_dp, 0.040_dp, 1.0E-5_dp]
      occ_noisy(:) = [1.30_dp, 0.72_dp, 0.08_dp]
      shared_cap = dense_shared_occ_slope(occ_curv(2), occ_curv(3))
      max_cap = 2.0_dp*MAX(occ_curv(2), occ_curv(3))

      CALL dense_spectral_occ_projector_derivative_chc(chc, dchc, occ_noisy, docc_mat, occ_curv, gauge_gap)
      raw_coeff = (occ_noisy(2) - occ_noisy(3))/(REAL(chc(2, 2), KIND=dp) - &
                                                 REAL(chc(3, 3), KIND=dp))
      capped_coeff = REAL(docc_mat(2, 3), KIND=dp)/REAL(dchc(2, 3), KIND=dp)
      err = ABS(capped_coeff + shared_cap)

      IF (err > tol .OR. ABS(raw_coeff) < 10.0_dp*max_cap .OR. shared_cap >= 0.1_dp*max_cap) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,5(1X,ES12.5))') TRIM(label), &
            "failed with error/raw/shared/max/capped", err, raw_coeff, shared_cap, max_cap, capped_coeff
      ELSE
         WRITE (*, '(A,1X,A,5(1X,ES12.5))') TRIM(label), &
            "OK with error/raw/shared/max/capped", err, raw_coeff, shared_cap, max_cap, capped_coeff
      END IF

   END SUBROUTINE check_occ_projector_active_inert_pair_cap

! **************************************************************************************************
!> \brief Keep far active/inert occupation-projector slopes out of the near-gauge cap.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_occ_projector_far_active_inert_pair_max_cap(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 3
      REAL(KIND=dp), PARAMETER                           :: gap = 0.20_dp, gauge_gap = 3.0E-5_dp

      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: chc, dchc, docc_mat
      REAL(KIND=dp)                                      :: capped_coeff, err, max_cap, raw_coeff, &
                                                            shared_cap
      REAL(KIND=dp), DIMENSION(n)                        :: occ_curv, occ_noisy

      chc(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      chc(1, 1) = CMPLX(-0.08_dp, 0.0_dp, KIND=dp)
      chc(2, 2) = CMPLX(-0.08_dp + gap, 0.0_dp, KIND=dp)
      chc(3, 3) = CMPLX(0.43_dp, 0.0_dp, KIND=dp)

      dchc(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      dchc(1, 2) = CMPLX(0.19_dp, -0.13_dp, KIND=dp)
      dchc(2, 1) = CONJG(dchc(1, 2))

      occ_curv(:) = [0.040_dp, 1.0E-5_dp, 0.003_dp]
      occ_noisy(:) = [0.90_dp, 0.02_dp, 0.01_dp]
      shared_cap = dense_shared_occ_slope(occ_curv(1), occ_curv(2))
      max_cap = 2.0_dp*MAX(occ_curv(1), occ_curv(2))

      CALL dense_spectral_occ_projector_derivative_chc(chc, dchc, occ_noisy, docc_mat, occ_curv, gauge_gap)
      raw_coeff = (occ_noisy(1) - occ_noisy(2))/(REAL(chc(1, 1), KIND=dp) - &
                                                 REAL(chc(2, 2), KIND=dp))
      capped_coeff = REAL(docc_mat(1, 2), KIND=dp)/REAL(dchc(1, 2), KIND=dp)
      err = ABS(capped_coeff + max_cap)

      IF (err > tol .OR. ABS(raw_coeff) < 10.0_dp*max_cap .OR. &
          shared_cap >= 0.05_dp*max_cap .OR. ABS(gap) <= 4.0_dp*gauge_gap .OR. &
          ABS(capped_coeff + shared_cap) < 0.5_dp*max_cap) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,5(1X,ES12.5))') TRIM(label), &
            "failed with error/raw/shared/max/capped", err, raw_coeff, shared_cap, max_cap, capped_coeff
      ELSE
         WRITE (*, '(A,1X,A,5(1X,ES12.5))') TRIM(label), &
            "OK with error/raw/shared/max/capped", err, raw_coeff, shared_cap, max_cap, capped_coeff
      END IF

   END SUBROUTINE check_occ_projector_far_active_inert_pair_max_cap

! **************************************************************************************************
!> \brief Plateau only anomalous near-degenerate projector-basis occupation jumps.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_complex_occ_projector_gauge_plateau(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 3
      REAL(KIND=dp), PARAMETER                           :: fd_step = 1.0E-8_dp, gap = 2.0E-5_dp, &
                                                            gauge_gap = 3.0E-5_dp

      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: chc, dchc, docc_mat, fd_mat, occ_mat_m, &
                                                            occ_mat_p
      REAL(KIND=dp)                                      :: cap, err, plateau_coeff, raw_coeff
      REAL(KIND=dp), DIMENSION(n)                        :: occ_curv, occ_equal, occ_noisy, &
                                                            occ_smooth

      chc(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      chc(1, 1) = CMPLX(0.00_dp, 0.0_dp, KIND=dp)
      chc(2, 2) = CMPLX(gap, 0.0_dp, KIND=dp)
      chc(3, 3) = CMPLX(0.39_dp, 0.0_dp, KIND=dp)

      dchc(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      dchc(1, 2) = CMPLX(0.23_dp, -0.14_dp, KIND=dp)
      dchc(2, 1) = CONJG(dchc(1, 2))

      occ_curv(:) = [0.031_dp, 0.027_dp, 0.004_dp]
      cap = 2.0_dp*MAX(occ_curv(1), occ_curv(2))

      occ_equal(:) = [0.75_dp, 0.75_dp, 0.12_dp]
      CALL dense_spectral_occ_projector_derivative_chc(chc, dchc, occ_equal, docc_mat, occ_curv, gauge_gap)
      err = MAXVAL(ABS(docc_mat(:, :)))

      occ_smooth(:) = [0.75_dp, 0.75_dp - 0.25_dp*cap*gap, 0.12_dp]
      CALL dense_spectral_occ_projector_derivative_chc(chc, dchc, occ_smooth, docc_mat, occ_curv, gauge_gap)
      CALL dense_spectral_occ_projector_from_chc(chc + fd_step*dchc, occ_smooth, occ_mat_p)
      CALL dense_spectral_occ_projector_from_chc(chc - fd_step*dchc, occ_smooth, occ_mat_m)
      fd_mat(:, :) = (occ_mat_p(:, :) - occ_mat_m(:, :))/(2.0_dp*fd_step)
      err = MAX(err, MAXVAL(ABS(docc_mat(:, :) - fd_mat(:, :))))

      occ_noisy(:) = [0.75_dp, 0.21_dp, 0.12_dp]
      CALL dense_spectral_occ_projector_derivative_chc(chc, dchc, occ_noisy, docc_mat, occ_curv, gauge_gap)
      raw_coeff = (occ_noisy(1) - occ_noisy(2))/(REAL(chc(1, 1), KIND=dp) - &
                                                 REAL(chc(2, 2), KIND=dp))
      plateau_coeff = REAL(docc_mat(1, 2), KIND=dp)/REAL(dchc(1, 2), KIND=dp)
      err = MAX(err, ABS(plateau_coeff))

      IF (err > tol .OR. ABS(raw_coeff) < 10.0_dp*cap) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), "failed with error/raw/cap/plateau", &
            err, raw_coeff, cap, plateau_coeff
      ELSE
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), "OK with error/raw/cap/plateau", &
            err, raw_coeff, cap, plateau_coeff
      END IF

   END SUBROUTINE check_complex_occ_projector_gauge_plateau

! **************************************************************************************************
!> \brief Keep the gauge plateau local, but thermally cap non-Fermi divided-difference slopes.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_occ_projector_plateau_caps_far_slope(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 3
      REAL(KIND=dp), PARAMETER                           :: gap = 1.5E-4_dp, gauge_gap = 3.0E-5_dp

      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: chc, dchc, docc_plateau, docc_raw
      REAL(KIND=dp)                                      :: cap, err, plateau_coeff, raw_coeff
      REAL(KIND=dp), DIMENSION(n)                        :: occ_curv, occ_noisy

      chc(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      chc(1, 1) = CMPLX(0.00_dp, 0.0_dp, KIND=dp)
      chc(2, 2) = CMPLX(gap, 0.0_dp, KIND=dp)
      chc(3, 3) = CMPLX(0.39_dp, 0.0_dp, KIND=dp)

      dchc(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      dchc(1, 2) = CMPLX(0.23_dp, -0.14_dp, KIND=dp)
      dchc(2, 1) = CONJG(dchc(1, 2))

      occ_curv(:) = [0.031_dp, 0.027_dp, 0.004_dp]
      occ_noisy(:) = [0.75_dp, 0.21_dp, 0.12_dp]
      cap = dense_projector_pair_slope_cap(-gap, gauge_gap, occ_curv(1), occ_curv(2), &
                                           MAXVAL(MAX(0.0_dp, occ_curv(:))))

      CALL dense_spectral_occ_projector_derivative_chc(chc, dchc, occ_noisy, docc_raw)
      CALL dense_spectral_occ_projector_derivative_chc(chc, dchc, occ_noisy, docc_plateau, &
                                                       occ_curv, gauge_gap)
      raw_coeff = REAL(docc_raw(1, 2), KIND=dp)/REAL(dchc(1, 2), KIND=dp)
      plateau_coeff = REAL(docc_plateau(1, 2), KIND=dp)/REAL(dchc(1, 2), KIND=dp)
      err = ABS(plateau_coeff + cap)

      IF (err > tol .OR. gap <= gauge_gap .OR. ABS(raw_coeff) < 10.0_dp*cap) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,5(1X,ES12.5))') TRIM(label), "failed with error/raw/cap/capped/gap", &
            err, raw_coeff, plateau_coeff, cap, gap
      ELSE
         WRITE (*, '(A,1X,A,5(1X,ES12.5))') TRIM(label), "OK with error/raw/cap/capped/gap", &
            err, raw_coeff, plateau_coeff, cap, gap
      END IF

   END SUBROUTINE check_occ_projector_plateau_caps_far_slope

! **************************************************************************************************
!> \brief Keep the physical dP response covariant while removing near-degenerate gauge noise.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_complex_occ_projector_gauge_covariance(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 4
      REAL(KIND=dp), PARAMETER                           :: active_tol = 1.0E-4_dp, &
                                                            gauge_gap = 3.0E-5_dp

      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: chc, chc_q, dchc, dchc_q, dchc_shift, &
                                                            docc_mat, docc_q, docc_raw, docc_ref, q
      INTEGER                                            :: i
      REAL(KIND=dp)                                      :: active, cov_err, err, gauge_err, &
                                                            trace_err
      REAL(KIND=dp), DIMENSION(n)                        :: occ_curv, occ_scale
      REAL(KIND=dp), DIMENSION(n, n)                     :: x_im, x_re

      chc(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      chc(1, 1) = CMPLX(-0.010000_dp, 0.0_dp, KIND=dp)
      chc(2, 2) = CMPLX(-0.009982_dp, 0.0_dp, KIND=dp)
      chc(3, 3) = CMPLX(0.21_dp, 0.0_dp, KIND=dp)
      chc(4, 4) = CMPLX(0.43_dp, 0.0_dp, KIND=dp)

      dchc(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      dchc(1, 2) = CMPLX(0.24_dp, -0.18_dp, KIND=dp)
      dchc(1, 3) = CMPLX(-0.06_dp, 0.04_dp, KIND=dp)
      dchc(2, 4) = CMPLX(0.05_dp, 0.03_dp, KIND=dp)
      dchc(3, 4) = CMPLX(-0.07_dp, 0.02_dp, KIND=dp)
      DO i = 2, n
         dchc(i, 1:i - 1) = CONJG(dchc(1:i - 1, i))
      END DO

      occ_scale(:) = [1.52_dp, 0.83_dp, 0.31_dp, 0.08_dp]
      occ_curv(:) = [0.033_dp, 0.030_dp, 0.011_dp, 0.002_dp]

      CALL dense_spectral_occ_projector_derivative_chc(chc, dchc, occ_scale, docc_raw)
      CALL dense_spectral_occ_projector_derivative_chc(chc, dchc, occ_scale, docc_mat, &
                                                       occ_curv, gauge_gap)

      dchc_shift(:, :) = dchc(:, :)
      DO i = 1, n
         dchc_shift(i, i) = dchc_shift(i, i) + CMPLX(0.37_dp, 0.0_dp, KIND=dp)
      END DO
      CALL dense_spectral_occ_projector_derivative_chc(chc, dchc_shift, occ_scale, docc_ref, &
                                                       occ_curv, gauge_gap)

      x_re(:, :) = 0.0_dp
      x_re(1, 2) = 0.11_dp
      x_re(1, 3) = -0.07_dp
      x_re(2, 4) = 0.05_dp
      x_re(3, 4) = -0.04_dp
      DO i = 2, n
         x_re(i, 1:i - 1) = -x_re(1:i - 1, i)
      END DO

      x_im(:, :) = 0.0_dp
      x_im(1, 1) = 0.03_dp
      x_im(2, 2) = -0.02_dp
      x_im(3, 3) = 0.01_dp
      x_im(4, 4) = -0.04_dp
      x_im(1, 2) = -0.05_dp
      x_im(1, 4) = 0.04_dp
      x_im(2, 3) = 0.02_dp
      x_im(3, 4) = -0.03_dp
      DO i = 2, n
         x_im(i, 1:i - 1) = x_im(1:i - 1, i)
      END DO

      CALL dense_rotation_u(x_re, x_im, q)
      chc_q(:, :) = MATMUL(CONJG(TRANSPOSE(q)), MATMUL(chc, q))
      dchc_q(:, :) = MATMUL(CONJG(TRANSPOSE(q)), MATMUL(dchc, q))
      CALL dense_spectral_occ_projector_derivative_chc(chc_q, dchc_q, occ_scale, docc_q, &
                                                       occ_curv, gauge_gap)
      docc_ref(:, :) = MATMUL(CONJG(TRANSPOSE(q)), MATMUL(docc_mat, q))

      active = MAXVAL(ABS(docc_raw(:, :) - docc_mat(:, :)))
      gauge_err = MAXVAL(ABS(docc_ref(:, :) - docc_q(:, :)))
      cov_err = gauge_err
      trace_err = ABS(REAL(SUM([(docc_mat(i, i), i=1, n)]), KIND=dp))
      err = MAX(gauge_err, trace_err)
      CALL dense_spectral_occ_projector_derivative_chc(chc, dchc_shift, occ_scale, docc_ref, &
                                                       occ_curv, gauge_gap)
      err = MAX(err, MAXVAL(ABS(docc_mat(:, :) - docc_ref(:, :))))

      IF (err > tol .OR. active < active_tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "failed with error/active/cov/trace", err, active, cov_err, trace_err
      ELSE
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "OK with error/active/cov/trace", err, active, cov_err, trace_err
      END IF

   END SUBROUTINE check_complex_occ_projector_gauge_covariance

! **************************************************************************************************
!> \brief Check that diagonal projector-density response enters through the dense fixed-N projection.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_complex_projector_density_response_projection(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 4, nrot = 3
      REAL(KIND=dp), PARAMETER                           :: fd_step = 1.0E-7_dp

      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: chc, dchc, docc_fd, docc_mat, occ_mat_m, &
                                                            occ_mat_p
      INTEGER                                            :: i, s
      REAL(KIND=dp)                                      :: active, corr_shift, err, occ_sum, &
                                                            raw_shift
      REAL(KIND=dp), DIMENSION(n)                        :: density_diag, occ_curv, occ_scale
      REAL(KIND=dp), DIMENSION(n, nrot)                  :: b_current, b_reference, raw_response, &
                                                            response_current

      chc(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      chc(1, 1) = CMPLX(-0.31_dp, 0.0_dp, KIND=dp)
      chc(2, 2) = CMPLX(-0.04_dp, 0.0_dp, KIND=dp)
      chc(3, 3) = CMPLX(0.18_dp, 0.0_dp, KIND=dp)
      chc(4, 4) = CMPLX(0.43_dp, 0.0_dp, KIND=dp)
      chc(1, 2) = CMPLX(0.045_dp, -0.025_dp, KIND=dp)
      chc(1, 3) = CMPLX(-0.020_dp, 0.030_dp, KIND=dp)
      chc(1, 4) = CMPLX(0.025_dp, 0.010_dp, KIND=dp)
      chc(2, 3) = CMPLX(0.035_dp, 0.020_dp, KIND=dp)
      chc(2, 4) = CMPLX(-0.030_dp, 0.015_dp, KIND=dp)
      chc(3, 4) = CMPLX(0.020_dp, -0.035_dp, KIND=dp)
      DO i = 2, n
         chc(i, 1:i - 1) = CONJG(chc(1:i - 1, i))
      END DO

      occ_scale = [1.64_dp, 1.08_dp, 0.47_dp, 0.12_dp]
      occ_curv = [0.38_dp, 0.71_dp, 0.29_dp, 0.16_dp]
      raw_response(:, 1) = [0.12_dp, -0.05_dp, 0.19_dp, -0.08_dp]
      raw_response(:, 2) = [-0.17_dp, 0.09_dp, 0.04_dp, 0.13_dp]
      raw_response(:, 3) = [0.06_dp, 0.21_dp, -0.11_dp, -0.03_dp]

      err = 0.0_dp
      active = 0.0_dp
      occ_sum = SUM(occ_curv)
      DO s = 1, nrot
         dchc(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
         SELECT CASE (s)
         CASE (1)
            dchc(1, 2) = CMPLX(0.27_dp, -0.19_dp, KIND=dp)
            dchc(2, 4) = CMPLX(-0.08_dp, 0.11_dp, KIND=dp)
         CASE (2)
            dchc(1, 3) = CMPLX(-0.22_dp, 0.14_dp, KIND=dp)
            dchc(3, 4) = CMPLX(0.09_dp, 0.16_dp, KIND=dp)
         CASE (3)
            dchc(1, 4) = CMPLX(0.18_dp, 0.07_dp, KIND=dp)
            dchc(2, 3) = CMPLX(-0.13_dp, -0.10_dp, KIND=dp)
         END SELECT
         DO i = 2, n
            dchc(i, 1:i - 1) = CONJG(dchc(1:i - 1, i))
         END DO
         CALL dense_spectral_occ_projector_derivative_chc(chc, dchc, occ_scale, docc_mat)
         CALL dense_spectral_occ_projector_from_chc(chc + fd_step*dchc, occ_scale, occ_mat_p)
         CALL dense_spectral_occ_projector_from_chc(chc - fd_step*dchc, occ_scale, occ_mat_m)
         docc_fd(:, :) = (occ_mat_p(:, :) - occ_mat_m(:, :))/(2.0_dp*fd_step)
         density_diag(:) = 0.0_dp
         DO i = 1, n
            density_diag(i) = REAL(docc_mat(i, i), KIND=dp)
         END DO

         response_current(:, s) = raw_response(:, s) - 0.5_dp*density_diag(:)/occ_curv(:)
         raw_shift = DOT_PRODUCT(occ_curv, response_current(:, s))/occ_sum
         b_current(:, s) = occ_curv(:)*(response_current(:, s) - raw_shift)

         raw_shift = DOT_PRODUCT(occ_curv, raw_response(:, s))/occ_sum
         corr_shift = SUM(density_diag(:))/occ_sum
         b_reference(:, s) = occ_curv(:)*(raw_response(:, s) - raw_shift) - &
                             0.5_dp*(density_diag(:) - occ_curv(:)*corr_shift)

         err = MAX(err, MAXVAL(ABS(docc_mat(:, :) - docc_fd(:, :))))
         err = MAX(err, MAXVAL(ABS(b_current(:, s) - b_reference(:, s))))
         active = MAX(active, MAXVAL(ABS(b_current(:, s) - &
                                         occ_curv(:)*(raw_response(:, s) - raw_shift))))
      END DO

      IF (err > tol .OR. active < 1.0E-4_dp) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,2(1X,ES12.5))') TRIM(label), "failed with error/active", err, active
      ELSE
         WRITE (*, '(A,1X,A,2(1X,ES12.5))') TRIM(label), "OK with error/active", err, active
      END IF

   END SUBROUTINE check_complex_projector_density_response_projection

! **************************************************************************************************
!> \brief Check that near-degenerate projector gauge motion is removed before fixed-N density response.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_complex_projector_density_gauge_plateau_projection(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 4
      REAL(KIND=dp), PARAMETER                           :: active_tol = 1.0E-4_dp, &
                                                            gauge_gap = 5.0E-5_dp

      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: chc, chc_eig, dchc, dchc_eig, &
                                                            docc_plateau, docc_raw, q
      INTEGER                                            :: i
      REAL(KIND=dp) :: err, occ_sum, physical_active, physical_capped_norm, physical_err, &
         physical_raw_norm, plateau_norm, raw_norm, shift
      REAL(KIND=dp), DIMENSION(n)                        :: occ_curv, occ_scale, projected_response
      REAL(KIND=dp), DIMENSION(n, n)                     :: x_im, x_re

      chc_eig(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      chc_eig(1, 1) = CMPLX(-0.010000_dp, 0.0_dp, KIND=dp)
      chc_eig(2, 2) = CMPLX(-0.009982_dp, 0.0_dp, KIND=dp)
      chc_eig(3, 3) = CMPLX(0.23_dp, 0.0_dp, KIND=dp)
      chc_eig(4, 4) = CMPLX(0.51_dp, 0.0_dp, KIND=dp)

      x_re(:, :) = 0.0_dp
      x_re(1, 2) = 0.13_dp
      x_re(1, 3) = -0.06_dp
      x_re(2, 4) = 0.04_dp
      x_re(3, 4) = -0.05_dp
      DO i = 2, n
         x_re(i, 1:i - 1) = -x_re(1:i - 1, i)
      END DO
      x_im(:, :) = 0.0_dp
      x_im(1, 1) = 0.02_dp
      x_im(2, 2) = -0.03_dp
      x_im(3, 3) = 0.01_dp
      x_im(4, 4) = -0.01_dp
      x_im(1, 2) = -0.04_dp
      x_im(1, 4) = 0.03_dp
      x_im(2, 3) = 0.02_dp
      DO i = 2, n
         x_im(i, 1:i - 1) = x_im(1:i - 1, i)
      END DO
      CALL dense_rotation_u(x_re, x_im, q)
      chc(:, :) = MATMUL(q, MATMUL(chc_eig, CONJG(TRANSPOSE(q))))

      occ_scale(:) = [1.70_dp, 0.52_dp, 0.23_dp, 0.07_dp]
      occ_curv(:) = [0.038_dp, 0.035_dp, 0.018_dp, 0.006_dp]
      occ_sum = SUM(occ_curv(:))

      dchc_eig(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      dchc_eig(1, 2) = CMPLX(0.21_dp, -0.17_dp, KIND=dp)
      dchc_eig(2, 1) = CONJG(dchc_eig(1, 2))
      dchc(:, :) = MATMUL(q, MATMUL(dchc_eig, CONJG(TRANSPOSE(q))))
      CALL dense_spectral_occ_projector_derivative_chc(chc, dchc, occ_scale, docc_raw)
      CALL dense_spectral_occ_projector_derivative_chc(chc, dchc, occ_scale, docc_plateau, &
                                                       occ_curv, gauge_gap)

      DO i = 1, n
         projected_response(i) = 0.5_dp*REAL(docc_raw(i, i), KIND=dp)/occ_curv(i)
      END DO
      shift = DOT_PRODUCT(occ_curv(:), projected_response(:))/occ_sum
      raw_norm = MAXVAL(ABS(occ_curv(:)*(projected_response(:) - shift)))
      DO i = 1, n
         projected_response(i) = 0.5_dp*REAL(docc_plateau(i, i), KIND=dp)/occ_curv(i)
      END DO
      shift = DOT_PRODUCT(occ_curv(:), projected_response(:))/occ_sum
      plateau_norm = MAXVAL(ABS(occ_curv(:)*(projected_response(:) - shift)))

      dchc_eig(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      dchc_eig(3, 4) = CMPLX(-0.16_dp, 0.11_dp, KIND=dp)
      dchc_eig(4, 3) = CONJG(dchc_eig(3, 4))
      dchc(:, :) = MATMUL(q, MATMUL(dchc_eig, CONJG(TRANSPOSE(q))))
      CALL dense_spectral_occ_projector_derivative_chc(chc, dchc, occ_scale, docc_raw)
      CALL dense_spectral_occ_projector_derivative_chc(chc, dchc, occ_scale, docc_plateau, &
                                                       occ_curv, gauge_gap)
      physical_raw_norm = MAXVAL(ABS(docc_raw(:, :)))
      physical_capped_norm = MAXVAL(ABS(docc_plateau(:, :)))
      physical_err = MAX(0.0_dp, physical_capped_norm - physical_raw_norm)
      physical_active = 0.0_dp
      DO i = 1, n
         projected_response(i) = 0.5_dp*REAL(docc_plateau(i, i), KIND=dp)/occ_curv(i)
      END DO
      shift = DOT_PRODUCT(occ_curv(:), projected_response(:))/occ_sum
      physical_active = MAXVAL(ABS(occ_curv(:)*(projected_response(:) - shift)))

      err = MAX(plateau_norm, physical_err)
      IF (err > tol .OR. raw_norm < active_tol .OR. physical_active < active_tol .OR. &
          physical_capped_norm >= 0.5_dp*physical_raw_norm) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "failed with error/raw/plateau/physical", err, raw_norm, plateau_norm, physical_active
      ELSE
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "OK with error/raw/plateau/physical", err, raw_norm, plateau_norm, physical_active
      END IF

   END SUBROUTINE check_complex_projector_density_gauge_plateau_projection

! **************************************************************************************************
!> \brief Check that the production density-response filter is separate from finite dP physics.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_complex_projector_active_density_response_filter(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 4
      REAL(KIND=dp), PARAMETER                           :: active_tol = 1.0E-4_dp, &
                                                            fd_step = 1.0E-7_dp, &
                                                            gauge_gap = 5.0E-5_dp, &
                                                            projector_weight = 0.37_dp

      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: chc, chc_eig, dchc, dchc_eig, docc_fd, &
                                                            docc_safe, occ_mat_m, occ_mat_p, q
      INTEGER                                            :: i
      REAL(KIND=dp) :: active_bound_err, active_norm, active_plateau_gap, active_shift, &
         active_sum, active_weight, fixed_n_leak, full_err, full_norm, gauge_leak_err, gauge_norm, &
         gauge_shift, low_ratio_err, low_weight, occ_curv_max, occ_sum, old_fixed_n_leak, &
         plateau_fixed_n_leak, plateau_norm, plateau_shift
      REAL(KIND=dp), DIMENSION(n) :: active_response, full_response, occ_curv, occ_plateau, &
         occ_scale, plateau_response, reference_response, response_shifted, weight_response
      REAL(KIND=dp), DIMENSION(n, n)                     :: x_im, x_re

      chc_eig(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      chc_eig(1, 1) = CMPLX(-0.010000_dp, 0.0_dp, KIND=dp)
      chc_eig(2, 2) = CMPLX(-0.009982_dp, 0.0_dp, KIND=dp)
      chc_eig(3, 3) = CMPLX(0.23_dp, 0.0_dp, KIND=dp)
      chc_eig(4, 4) = CMPLX(0.51_dp, 0.0_dp, KIND=dp)

      x_re(:, :) = 0.0_dp
      x_re(1, 2) = 0.13_dp
      x_re(1, 3) = -0.06_dp
      x_re(2, 4) = 0.04_dp
      x_re(3, 4) = -0.05_dp
      DO i = 2, n
         x_re(i, 1:i - 1) = -x_re(1:i - 1, i)
      END DO
      x_im(:, :) = 0.0_dp
      x_im(1, 1) = 0.02_dp
      x_im(2, 2) = -0.03_dp
      x_im(3, 3) = 0.01_dp
      x_im(4, 4) = -0.01_dp
      x_im(1, 2) = -0.04_dp
      x_im(1, 4) = 0.03_dp
      x_im(2, 3) = 0.02_dp
      DO i = 2, n
         x_im(i, 1:i - 1) = x_im(1:i - 1, i)
      END DO
      CALL dense_rotation_u(x_re, x_im, q)
      chc(:, :) = MATMUL(q, MATMUL(chc_eig, CONJG(TRANSPOSE(q))))

      occ_scale(:) = [1.70_dp, 0.52_dp, 0.24_dp, 0.2316_dp]
      occ_plateau(:) = occ_scale(:)
      occ_plateau(1:2) = SUM(occ_scale(1:2))/2.0_dp
      occ_curv(:) = [0.038_dp, 0.035_dp, 0.018_dp, 0.0019_dp]
      occ_curv_max = MAXVAL(occ_curv(:))

      dchc_eig(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      dchc_eig(3, 4) = CMPLX(-0.16_dp, 0.11_dp, KIND=dp)
      dchc_eig(4, 3) = CONJG(dchc_eig(3, 4))
      dchc(:, :) = MATMUL(q, MATMUL(dchc_eig, CONJG(TRANSPOSE(q))))

      CALL dense_history_safe_occ_projector_derivative_chc(chc, dchc, occ_scale, projector_weight, &
                                                           docc_safe, occ_curv, gauge_gap)
      CALL dense_spectral_occ_projector_from_chc(chc + fd_step*dchc, occ_plateau, occ_mat_p)
      CALL dense_spectral_occ_projector_from_chc(chc - fd_step*dchc, occ_plateau, occ_mat_m)
      docc_fd(:, :) = projector_weight*(occ_mat_p(:, :) - occ_mat_m(:, :))/(2.0_dp*fd_step)

      active_response(:) = 0.0_dp
      full_response(:) = 0.0_dp
      reference_response(:) = 0.0_dp
      weight_response(:) = 0.0_dp
      DO i = 1, n
         full_response(i) = 0.5_dp*REAL(docc_safe(i, i), KIND=dp)/occ_curv(i)
         reference_response(i) = 0.5_dp*REAL(docc_fd(i, i), KIND=dp)/occ_curv(i)
         active_weight = 0.0_dp
         IF (occ_curv(i) > EPSILON(occ_curv(i)) .AND. occ_curv_max > EPSILON(occ_curv_max)) &
            active_weight = SQRT(MIN(1.0_dp, MAX(0.0_dp, occ_curv(i)/occ_curv_max)))
         weight_response(i) = active_weight
      END DO
      occ_sum = SUM(occ_curv(:))
      plateau_shift = DOT_PRODUCT(occ_curv(:), full_response(:))/occ_sum
      plateau_response(:) = full_response(:) - plateau_shift
      response_shifted(:) = plateau_response(:)
      old_fixed_n_leak = ABS(DOT_PRODUCT(occ_curv(:), weight_response(:)*response_shifted(:)))/occ_sum
      active_sum = DOT_PRODUCT(occ_curv(:), weight_response(:))
      active_shift = DOT_PRODUCT(occ_curv(:)*weight_response(:), full_response(:))/active_sum
      response_shifted(:) = full_response(:) - active_shift
      active_response(:) = weight_response(:)*(full_response(:) - active_shift)
      plateau_fixed_n_leak = ABS(DOT_PRODUCT(occ_curv(:), plateau_response(:)))/occ_sum
      fixed_n_leak = ABS(DOT_PRODUCT(occ_curv(:), active_response(:)))/occ_sum

      full_err = MAXVAL(ABS(full_response(:) - reference_response(:)))
      active_bound_err = MAXVAL(MAX(0.0_dp, ABS(active_response(:)) - ABS(response_shifted(:))))
      full_norm = MAXVAL(ABS(response_shifted(:)))
      plateau_norm = MAXVAL(ABS(plateau_response(:)))
      active_norm = MAXVAL(ABS(active_response(:)))
      active_plateau_gap = MAXVAL(ABS(plateau_response(:) - active_response(:)))
      low_weight = SQRT(occ_curv(n)/occ_curv_max)
      low_ratio_err = ABS(active_response(n)/response_shifted(n) - low_weight)
      gauge_shift = 0.42_dp
      gauge_norm = MAXVAL(ABS(weight_response(:)*(0.42_dp - gauge_shift)))
      gauge_leak_err = MAXVAL(ABS(weight_response(:)*0.42_dp))

      IF (full_err > tol .OR. active_bound_err > tol .OR. low_ratio_err > tol .OR. &
          fixed_n_leak > tol .OR. plateau_fixed_n_leak > tol .OR. old_fixed_n_leak < active_tol .OR. &
          gauge_norm > tol .OR. gauge_leak_err < active_tol .OR. &
          full_norm < active_tol .OR. active_norm >= full_norm .OR. &
          plateau_norm <= active_norm .OR. active_plateau_gap < active_tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,12(1X,ES12.5))') TRIM(label), &
            "failed with error/bound/low/fixed/pfixed/old/gauge/leak/full/plateau/active/gap", &
            full_err, active_bound_err, low_ratio_err, fixed_n_leak, plateau_fixed_n_leak, &
            old_fixed_n_leak, gauge_norm, gauge_leak_err, full_norm, plateau_norm, active_norm, &
            active_plateau_gap
      ELSE
         WRITE (*, '(A,1X,A,12(1X,ES12.5))') TRIM(label), &
            "OK with error/bound/low/fixed/pfixed/old/gauge/leak/full/plateau/active/gap", &
            full_err, active_bound_err, low_ratio_err, fixed_n_leak, plateau_fixed_n_leak, &
            old_fixed_n_leak, gauge_norm, gauge_leak_err, full_norm, plateau_norm, active_norm, &
            active_plateau_gap
      END IF

   END SUBROUTINE check_complex_projector_active_density_response_filter

! **************************************************************************************************
!> \brief Check that the active density response does not amplify low-curvature projector noise.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_complex_projector_active_density_safe_inverse(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 4
      REAL(KIND=dp), PARAMETER :: active_tol = 1.0E-4_dp, response_ratio_gate = 2.0E-2_dp

      INTEGER                                            :: i
      REAL(KIND=dp)                                      :: active_floor, active_shift, active_sum, &
                                                            err, fixed_n_leak, low_suppression, &
                                                            old_shift, old_sum
      REAL(KIND=dp), DIMENSION(n)                        :: active_response, active_weight, &
                                                            docc_diag, occ_curv, old_response, &
                                                            raw_column, safe_column

      occ_curv(:) = [4.0E-2_dp, 1.2E-2_dp, 1.0E-3_dp, 4.0E-8_dp]
      docc_diag(:) = [-8.0E-3_dp, 2.4E-3_dp, -2.0E-4_dp, 1.6E-5_dp]
      active_floor = response_ratio_gate*MAXVAL(occ_curv(:))

      active_weight(:) = 0.0_dp
      raw_column(:) = 0.0_dp
      safe_column(:) = 0.0_dp
      DO i = 1, n
         raw_column(i) = 0.5_dp*docc_diag(i)/occ_curv(i)
         safe_column(i) = 0.5_dp*docc_diag(i)/MAX(occ_curv(i), active_floor)
         active_weight(i) = SQRT(MIN(1.0_dp, MAX(0.0_dp, occ_curv(i)/MAXVAL(occ_curv(:)))))
      END DO

      active_sum = DOT_PRODUCT(occ_curv(:), active_weight(:))
      active_shift = DOT_PRODUCT(occ_curv(:)*active_weight(:), safe_column(:))/active_sum
      active_response(:) = active_weight(:)*(safe_column(:) - active_shift)

      old_sum = DOT_PRODUCT(occ_curv(:), active_weight(:))
      old_shift = DOT_PRODUCT(occ_curv(:)*active_weight(:), raw_column(:))/old_sum
      old_response(:) = active_weight(:)*(raw_column(:) - old_shift)

      fixed_n_leak = ABS(DOT_PRODUCT(occ_curv(:), active_response(:)))
      low_suppression = ABS(active_response(n))/MAX(ABS(old_response(n)), EPSILON(1.0_dp))
      err = fixed_n_leak

      IF (err > tol .OR. MAXVAL(ABS(active_response(1:2))) < active_tol .OR. &
          ABS(old_response(n)) < active_tol .OR. low_suppression > 2.5E-2_dp) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "failed with error/active/oldlow/suppress", err, MAXVAL(ABS(active_response(1:2))), &
            ABS(old_response(n)), low_suppression
      ELSE
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "OK with error/active/oldlow/suppress", err, MAXVAL(ABS(active_response(1:2))), &
            ABS(old_response(n)), low_suppression
      END IF

   END SUBROUTINE check_complex_projector_active_density_safe_inverse

! **************************************************************************************************
!> \brief Check that bounded physical density-response mix recovers active curvature without gauge blow-up.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_complex_projector_bounded_density_mix(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 4
      REAL(KIND=dp), PARAMETER :: active_tol = 1.0E-4_dp, response_ratio_gate = 2.0E-2_dp

      INTEGER                                            :: i
      REAL(KIND=dp) :: active_floor, active_shift, active_sum, fixed_n_leak, low_suppression, &
         mid_recovery, mix_shift, mix_sum, mix_weight, old_shift, old_sum
      REAL(KIND=dp), DIMENSION(n) :: active_response, active_weight, docc_diag, mixed_column, &
         mixed_response, occ_curv, old_response, raw_column, safe_column

      occ_curv(:) = [4.0E-2_dp, 1.2E-2_dp, 4.0E-4_dp, 4.0E-8_dp]
      docc_diag(:) = [-8.0E-3_dp, 2.4E-3_dp, -1.6E-4_dp, 1.6E-5_dp]
      active_floor = response_ratio_gate*MAXVAL(occ_curv(:))

      active_weight(:) = 0.0_dp
      raw_column(:) = 0.0_dp
      safe_column(:) = 0.0_dp
      mixed_column(:) = 0.0_dp
      DO i = 1, n
         raw_column(i) = 0.5_dp*docc_diag(i)/occ_curv(i)
         safe_column(i) = 0.5_dp*docc_diag(i)/MAX(occ_curv(i), active_floor)
         mix_weight = occ_curv(i)/MAX(occ_curv(i), active_floor)
         mix_weight = MIN(1.0_dp, MAX(0.0_dp, mix_weight))
         mixed_column(i) = safe_column(i) + mix_weight*(raw_column(i) - safe_column(i))
         active_weight(i) = SQRT(MIN(1.0_dp, MAX(0.0_dp, occ_curv(i)/MAXVAL(occ_curv(:)))))
      END DO

      active_sum = DOT_PRODUCT(occ_curv(:), active_weight(:))
      active_shift = DOT_PRODUCT(occ_curv(:)*active_weight(:), safe_column(:))/active_sum
      active_response(:) = active_weight(:)*(safe_column(:) - active_shift)

      mix_sum = DOT_PRODUCT(occ_curv(:), active_weight(:))
      mix_shift = DOT_PRODUCT(occ_curv(:)*active_weight(:), mixed_column(:))/mix_sum
      mixed_response(:) = active_weight(:)*(mixed_column(:) - mix_shift)

      old_sum = DOT_PRODUCT(occ_curv(:), active_weight(:))
      old_shift = DOT_PRODUCT(occ_curv(:)*active_weight(:), raw_column(:))/old_sum
      old_response(:) = active_weight(:)*(raw_column(:) - old_shift)

      fixed_n_leak = ABS(DOT_PRODUCT(occ_curv(:), mixed_response(:)))
      mid_recovery = ABS(mixed_response(3)) - ABS(active_response(3))
      low_suppression = ABS(mixed_response(n))/MAX(ABS(old_response(n)), EPSILON(1.0_dp))

      IF (fixed_n_leak > tol .OR. mid_recovery < active_tol .OR. &
          ABS(old_response(n)) < active_tol .OR. low_suppression > 2.5E-2_dp .OR. &
          MAXVAL(ABS(mixed_response(1:2) - old_response(1:2))) > active_tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,5(1X,ES12.5))') TRIM(label), &
            "failed with fixed/mid/oldlow/suppress/highgap", fixed_n_leak, mid_recovery, &
            ABS(old_response(n)), low_suppression, MAXVAL(ABS(mixed_response(1:2) - old_response(1:2)))
      ELSE
         WRITE (*, '(A,1X,A,5(1X,ES12.5))') TRIM(label), &
            "OK with fixed/mid/oldlow/suppress/highgap", fixed_n_leak, mid_recovery, &
            ABS(old_response(n)), low_suppression, MAXVAL(ABS(mixed_response(1:2) - old_response(1:2)))
      END IF

   END SUBROUTINE check_complex_projector_bounded_density_mix

! **************************************************************************************************
!> \brief Check the production bounded density-response helper algebra.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_bounded_density_response_helper(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 4
      REAL(KIND=dp), PARAMETER :: active_tol = 1.0E-4_dp, response_ratio_gate = 2.0E-2_dp

      INTEGER                                            :: i
      REAL(KIND=dp)                                      :: active_fraction, err, inert_suppression, &
                                                            occ_curv_max
      REAL(KIND=dp), DIMENSION(n)                        :: expected, occ_curv, raw_response, &
                                                            response

      occ_curv(:) = [4.0E-2_dp, 1.2E-2_dp, 4.0E-4_dp, 4.0E-8_dp]
      raw_response(:) = [0.10_dp, -0.08_dp, 0.06_dp, -25.0_dp]
      occ_curv_max = MAXVAL(occ_curv(:))

      DO i = 1, n
         response(i) = dense_bounded_density_response(raw_response(i), occ_curv(i), occ_curv_max, &
                                                      response_ratio_gate)
         active_fraction = occ_curv(i)/MAX(occ_curv(i), response_ratio_gate*occ_curv_max)
         active_fraction = MIN(1.0_dp, MAX(0.0_dp, active_fraction))
         expected(i) = SQRT(MIN(1.0_dp, MAX(0.0_dp, occ_curv(i)/occ_curv_max)))*raw_response(i)* &
                       (2.0_dp*active_fraction - active_fraction*active_fraction)
      END DO

      err = MAXVAL(ABS(response(:) - expected(:)))
      inert_suppression = ABS(response(n))/MAX(ABS(raw_response(n)), EPSILON(1.0_dp))

      IF (err > tol .OR. ABS(response(1) - raw_response(1)) > tol .OR. &
          ABS(response(3)) < active_tol .OR. inert_suppression > 1.0E-5_dp) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "failed with error/active/mid/inert", err, response(1), response(3), inert_suppression
      ELSE
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "OK with error/active/mid/inert", err, response(1), response(3), inert_suppression
      END IF

   END SUBROUTINE check_projector_bounded_density_response_helper

! **************************************************************************************************
!> \brief Check the gauge-safe physical density-response helper algebra.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_physical_density_response_helper(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 4
      REAL(KIND=dp), PARAMETER :: active_tol = 1.0E-4_dp, response_ratio_gate = 2.0E-2_dp

      INTEGER                                            :: i
      REAL(KIND=dp) :: active_fraction, active_recovery, bounded_response, err, inert_suppression, &
         mid_recovery, occ_curv_max, projector_abs, response_shape, safe_abs, safe_blend, &
         safe_cap, safe_component, safe_rejection, zero_active, zero_inert
      REAL(KIND=dp), DIMENSION(n)                        :: expected, occ_curv, projector_response, &
                                                            response, safe_response

      occ_curv(:) = [4.0E-2_dp, 1.2E-2_dp, 4.0E-4_dp, 4.0E-8_dp]
      projector_response(:) = [0.010_dp, -0.008_dp, 0.006_dp, -0.004_dp]
      safe_response(:) = [0.10_dp, -0.08_dp, 0.06_dp, -25.0_dp]
      occ_curv_max = MAXVAL(occ_curv(:))

      DO i = 1, n
         response(i) = dense_physical_density_response(projector_response(i), safe_response(i), &
                                                       occ_curv(i), occ_curv_max, response_ratio_gate)
         active_fraction = occ_curv(i)/MAX(occ_curv(i), response_ratio_gate*occ_curv_max)
         active_fraction = MIN(1.0_dp, MAX(0.0_dp, active_fraction))
         response_shape = 2.0_dp*active_fraction - active_fraction*active_fraction
         safe_component = projector_response(i)
         IF (projector_response(i)*safe_response(i) > 0.0_dp .OR. &
             (ABS(projector_response(i)) <= EPSILON(projector_response(i)) .AND. &
              ABS(safe_response(i)) > EPSILON(safe_response(i)))) THEN
            projector_abs = ABS(projector_response(i))
            safe_abs = ABS(safe_response(i))
            safe_cap = MAX(4.0_dp*projector_abs, response_shape*safe_abs)
            safe_component = SIGN(MIN(safe_abs, safe_cap), safe_response(i))
            safe_blend = active_fraction
            safe_component = (1.0_dp - safe_blend)*projector_response(i) + &
                             safe_blend*safe_component
         END IF
         expected(i) = response_shape*safe_component
      END DO

      err = MAXVAL(ABS(response(:) - expected(:)))
      active_recovery = ABS(response(1)) - ABS(projector_response(1))
      safe_rejection = ABS(response(1) - safe_response(1))
      bounded_response = dense_bounded_density_response(projector_response(3), occ_curv(3), &
                                                        occ_curv_max, response_ratio_gate)
      mid_recovery = ABS(response(3)) - ABS(dense_bounded_density_response(projector_response(3), &
                                                                           occ_curv(3), occ_curv_max, &
                                                                           response_ratio_gate))
      inert_suppression = ABS(response(n))/MAX(ABS(safe_response(n)), EPSILON(1.0_dp))
      zero_active = dense_physical_density_response(0.0_dp, 0.125_dp, occ_curv(1), &
                                                    occ_curv_max, response_ratio_gate)
      zero_inert = dense_physical_density_response(0.0_dp, 0.125_dp, occ_curv(n), &
                                                   occ_curv_max, response_ratio_gate)

      IF (err > tol .OR. active_recovery < active_tol .OR. safe_rejection > tol .OR. &
          mid_recovery <= active_tol .OR. ABS(bounded_response) <= 0.0_dp .OR. &
          inert_suppression > 1.0E-4_dp .OR. ABS(zero_active - 0.125_dp) > tol .OR. &
          ABS(zero_inert) > 1.0E-4_dp*0.125_dp) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,7(1X,ES12.5))') TRIM(label), &
            "failed with error/active/mid/safegap/inert/zeroa/zeroi", err, response(1), mid_recovery, &
            safe_rejection, inert_suppression, zero_active, zero_inert
      ELSE
         WRITE (*, '(A,1X,A,7(1X,ES12.5))') TRIM(label), &
            "OK with error/active/mid/safegap/inert/zeroa/zeroi", err, response(1), mid_recovery, &
            safe_rejection, inert_suppression, zero_active, zero_inert
      END IF

   END SUBROUTINE check_projector_physical_density_response_helper

! **************************************************************************************************
!> \brief Check residual-aware filtering of unsupported projector-density response.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_residual_gauge_density_response(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      REAL(KIND=dp), PARAMETER :: response_ratio_gate = 2.0E-2_dp

      REAL(KIND=dp)                                      :: aligned_response, base_response, &
                                                            damped_response, err, good_response, &
                                                            occ_curv, occ_curv_max, &
                                                            partial_response

      occ_curv_max = 4.0E-2_dp
      occ_curv = 4.0E-4_dp
      base_response = dense_physical_density_response(1.0_dp, -0.10_dp, occ_curv, occ_curv_max, &
                                                      response_ratio_gate, 1.0_dp, -0.10_dp, 0.0_dp)
      good_response = residual_gauge_density_response(base_response, 1.0_dp, -0.10_dp, &
                                                      occ_curv, occ_curv_max, response_ratio_gate, &
                                                      1.0_dp, 1.0_dp, -0.10_dp, 0.0_dp)
      damped_response = residual_gauge_density_response(base_response, 1.0_dp, -0.10_dp, &
                                                        occ_curv, occ_curv_max, response_ratio_gate, &
                                                        0.0_dp, 1.0_dp, -0.10_dp, 0.0_dp)
      partial_response = residual_gauge_density_response(base_response, 1.0_dp, 0.25_dp, &
                                                         occ_curv, occ_curv_max, response_ratio_gate, &
                                                         0.0_dp, 1.0_dp, 0.25_dp, 0.0_dp)
      aligned_response = residual_gauge_density_response(base_response, 1.0_dp, 2.0_dp, &
                                                         occ_curv, occ_curv_max, response_ratio_gate, &
                                                         0.0_dp, 1.0_dp, 2.0_dp, 0.0_dp)

      err = ABS(good_response - base_response)
      err = MAX(err, ABS(damped_response - 0.75_dp*base_response))
      err = MAX(err, ABS(partial_response - 0.8125_dp*base_response))
      err = MAX(err, ABS(aligned_response - base_response))

      IF (err > tol .OR. ABS(damped_response) >= ABS(good_response) .OR. &
          ABS(partial_response) <= ABS(damped_response) .OR. &
          ABS(partial_response) >= ABS(good_response)) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,5(1X,ES12.5))') TRIM(label), &
            "failed with error/base/good/damped/partial", err, base_response, good_response, &
            damped_response, partial_response
      ELSE
         WRITE (*, '(A,1X,A,5(1X,ES12.5))') TRIM(label), &
            "OK with error/base/good/damped/partial", err, base_response, good_response, &
            damped_response, partial_response
      END IF

   END SUBROUTINE check_projector_residual_gauge_density_response

! **************************************************************************************************
!> \brief Check that the density-response blend compares fixed-N centered signs.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_physical_density_centered_sign(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      REAL(KIND=dp), PARAMETER :: response_ratio_gate = 2.0E-2_dp

      REAL(KIND=dp)                                      :: centered_response, err, occ_curv, &
                                                            occ_curv_max, raw_response, &
                                                            safe_response, uncentered_response

      occ_curv_max = 4.0E-2_dp
      occ_curv = occ_curv_max
      raw_response = 0.04_dp
      safe_response = -0.08_dp

      uncentered_response = dense_physical_density_response(raw_response, safe_response, &
                                                            occ_curv, occ_curv_max, response_ratio_gate)
      centered_response = dense_physical_density_response(raw_response, safe_response, &
                                                          occ_curv, occ_curv_max, response_ratio_gate, &
                                                          -0.02_dp, -0.04_dp)
      err = ABS(centered_response + 0.04_dp)

      IF (err > tol .OR. ABS(uncentered_response - raw_response) > tol .OR. &
          centered_response*safe_response <= 0.0_dp) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "failed with error/uncentered/centered/safe", err, uncentered_response, &
            centered_response, safe_response
      ELSE
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "OK with error/uncentered/centered/safe", err, uncentered_response, &
            centered_response, safe_response
      END IF

   END SUBROUTINE check_projector_physical_density_centered_sign

! **************************************************************************************************
!> \brief Check that fixed-N centered columns also set physical response magnitudes.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_physical_density_centered_magnitude(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      REAL(KIND=dp), PARAMETER :: response_ratio_gate = 2.0E-2_dp

      REAL(KIND=dp)                                      :: centered_response, err, legacy_gap, &
                                                            occ_curv, occ_curv_max, raw_centered, &
                                                            raw_response, safe_centered, safe_response

      occ_curv_max = 4.0E-2_dp
      occ_curv = occ_curv_max
      raw_response = 0.52_dp
      safe_response = 0.49_dp
      raw_centered = 0.02_dp
      safe_centered = -0.01_dp

      centered_response = dense_physical_density_response(raw_response, safe_response, &
                                                          occ_curv, occ_curv_max, response_ratio_gate, &
                                                          raw_centered, safe_centered)
      err = ABS(centered_response - safe_centered)
      legacy_gap = ABS(centered_response + MIN(ABS(raw_response), ABS(safe_response)))

      IF (err > tol .OR. ABS(centered_response) > ABS(raw_centered) + tol .OR. legacy_gap < 0.25_dp) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "failed with error/response/centered/legacygap", err, centered_response, &
            safe_centered, legacy_gap
      ELSE
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "OK with error/response/centered/legacygap", err, centered_response, &
            safe_centered, legacy_gap
      END IF

   END SUBROUTINE check_projector_physical_density_centered_magnitude

! **************************************************************************************************
!> \brief Check that a globally coherent fixed-N density response can recover local safe components.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_physical_density_global_alignment(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      REAL(KIND=dp), PARAMETER :: response_ratio_gate = 2.0E-2_dp

      REAL(KIND=dp)                                      :: aligned_response, blocked_expected, &
                                                            blocked_response, expected_response, occ_curv, &
                                                            occ_curv_max, raw_response, &
                                                            safe_response

      occ_curv_max = 4.0E-2_dp
      occ_curv = 2.0E-2_dp
      raw_response = 0.004_dp
      safe_response = -0.040_dp

      blocked_response = dense_physical_density_response(raw_response, safe_response, &
                                                         occ_curv, occ_curv_max, response_ratio_gate, &
                                                         0.002_dp, -0.020_dp, 0.0_dp)
      aligned_response = dense_physical_density_response(raw_response, safe_response, &
                                                         occ_curv, occ_curv_max, response_ratio_gate, &
                                                         0.002_dp, -0.020_dp, 0.80_dp)
      blocked_expected = -0.002_dp
      expected_response = -0.020_dp

      IF (ABS(blocked_response - blocked_expected) > tol .OR. &
          ABS(aligned_response - expected_response) > tol .OR. &
          ABS(aligned_response) <= ABS(blocked_response)) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "failed with blocked/blockexp/aligned/expected", blocked_response, blocked_expected, &
            aligned_response, expected_response
      ELSE
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "OK with blocked/blockexp/aligned/expected", blocked_response, blocked_expected, &
            aligned_response, expected_response
      END IF

   END SUBROUTINE check_projector_physical_density_global_alignment

! **************************************************************************************************
!> \brief Check that globally coherent fixed-N density response is not damped twice.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_physical_density_global_low_curv(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      REAL(KIND=dp), PARAMETER :: response_ratio_gate = 2.0E-2_dp

      REAL(KIND=dp) :: active_fraction, expected_response, global_response, local_response, &
         occ_curv, occ_curv_max, raw_response, response_shape, safe_response

      occ_curv_max = 4.0E-2_dp
      occ_curv = 4.0E-5_dp
      raw_response = 1.0E-5_dp
      safe_response = 1.0E-1_dp

      active_fraction = occ_curv/MAX(occ_curv, response_ratio_gate*occ_curv_max)
      active_fraction = MIN(1.0_dp, MAX(0.0_dp, active_fraction))
      response_shape = 2.0_dp*active_fraction - active_fraction*active_fraction
      expected_response = response_shape*((1.0_dp - active_fraction)*raw_response + &
                                          active_fraction*SQRT(response_shape)*safe_response)

      local_response = dense_physical_density_response(raw_response, safe_response, &
                                                       occ_curv, occ_curv_max, response_ratio_gate)
      global_response = dense_physical_density_response(raw_response, safe_response, &
                                                        occ_curv, occ_curv_max, response_ratio_gate, &
                                                        raw_response, safe_response, 0.80_dp)

      IF (ABS(global_response - expected_response) > tol .OR. &
          global_response <= 2.0_dp*local_response .OR. &
          global_response >= response_shape*safe_response) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,5(1X,ES12.5))') TRIM(label), &
            "failed with global/expected/local/shape/safe", global_response, expected_response, &
            local_response, response_shape, safe_response
      ELSE
         WRITE (*, '(A,1X,A,5(1X,ES12.5))') TRIM(label), &
            "OK with global/expected/local/shape/safe", global_response, expected_response, &
            local_response, response_shape, safe_response
      END IF

   END SUBROUTINE check_projector_physical_density_global_low_curv

! **************************************************************************************************
!> \brief Check that conflicting density responses are blended, not extrapolated.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_physical_density_response_convexity(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      REAL(KIND=dp), PARAMETER :: response_ratio_gate = 2.0E-2_dp

      REAL(KIND=dp) :: active_fraction, bounded_response, density_response, expected_response, &
         occ_curv, occ_curv_max, response_shape, safe_rejection, safe_response

      occ_curv_max = 4.0E-2_dp
      occ_curv = 4.0E-4_dp
      bounded_response = dense_bounded_density_response(1.0_dp, occ_curv, occ_curv_max, &
                                                        response_ratio_gate)
      safe_response = -1.0_dp
      density_response = dense_physical_density_response(1.0_dp, safe_response, occ_curv, occ_curv_max, &
                                                         response_ratio_gate)

      active_fraction = occ_curv/MAX(occ_curv, response_ratio_gate*occ_curv_max)
      active_fraction = MIN(1.0_dp, MAX(0.0_dp, active_fraction))
      response_shape = 2.0_dp*active_fraction - active_fraction*active_fraction
      expected_response = response_shape
      safe_rejection = ABS(density_response - safe_response)

      IF (ABS(bounded_response) < tol .OR. &
          ABS(density_response - expected_response) > tol .OR. &
          density_response <= bounded_response .OR. safe_rejection < 0.5_dp) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,5(1X,ES12.5))') TRIM(label), &
            "failed with bounded/response/expected/reject/safe", bounded_response, density_response, &
            expected_response, safe_rejection, safe_response
      ELSE
         WRITE (*, '(A,1X,A,5(1X,ES12.5))') TRIM(label), &
            "OK with bounded/response/expected/reject/safe", bounded_response, density_response, &
            expected_response, safe_rejection, safe_response
      END IF

   END SUBROUTINE check_projector_physical_density_response_convexity

! **************************************************************************************************
!> \brief Check the shadow-only finite-floor density-response helper algebra.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_density_floor_audit_response_helper(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 4
      REAL(KIND=dp), PARAMETER :: active_tol = 1.0E-4_dp, response_ratio_gate = 2.0E-2_dp

      INTEGER                                            :: i
      REAL(KIND=dp)                                      :: active_fraction, err, finite_floor, &
                                                            occ_curv_max, physical_recovery, &
                                                            response_shape, safe_rejection
      REAL(KIND=dp), DIMENSION(n)                        :: expected, floor_response, occ_curv, &
                                                            physical_response, projector_response, &
                                                            safe_response

      occ_curv(:) = [4.0E-2_dp, 1.2E-2_dp, 4.0E-4_dp, 4.0E-8_dp]
      projector_response(:) = [0.010_dp, -0.008_dp, 0.006_dp, -0.004_dp]
      safe_response(:) = [0.10_dp, -0.08_dp, 0.06_dp, -25.0_dp]
      occ_curv_max = MAXVAL(occ_curv(:))

      DO i = 1, n
         physical_response(i) = dense_physical_density_response(projector_response(i), safe_response(i), &
                                                                occ_curv(i), occ_curv_max, response_ratio_gate)
         floor_response(i) = dense_physical_density_floor_audit_response(projector_response(i), safe_response(i), &
                                                                         occ_curv(i), occ_curv_max, &
                                                                         response_ratio_gate)
         active_fraction = occ_curv(i)/MAX(occ_curv(i), response_ratio_gate*occ_curv_max)
         active_fraction = MIN(1.0_dp, MAX(0.0_dp, active_fraction))
         response_shape = 2.0_dp*active_fraction - active_fraction*active_fraction
         expected(i) = response_shape*projector_response(i)
      END DO

      err = MAXVAL(ABS(floor_response(:) - expected(:)))
      finite_floor = ABS(floor_response(n))/MAX(ABS(projector_response(n)), EPSILON(1.0_dp))
      physical_recovery = ABS(physical_response(3)) - ABS(floor_response(3))
      safe_rejection = MAXVAL(ABS(floor_response(:) - safe_response(:)))

      IF (err > tol .OR. ABS(floor_response(1) - physical_response(1)) < active_tol .OR. &
          finite_floor < 5.0E-5_dp .OR. physical_recovery <= active_tol .OR. &
          safe_rejection < active_tol .OR. &
          MAXVAL(ABS(floor_response(:))) > MAXVAL(ABS(projector_response(:))) + tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,5(1X,ES12.5))') TRIM(label), &
            "failed with error/floor/recovery/reject/max", err, finite_floor, physical_recovery, &
            safe_rejection, MAXVAL(ABS(floor_response(:)))
      ELSE
         WRITE (*, '(A,1X,A,5(1X,ES12.5))') TRIM(label), &
            "OK with error/floor/recovery/reject/max", err, finite_floor, physical_recovery, &
            safe_rejection, MAXVAL(ABS(floor_response(:)))
      END IF

   END SUBROUTINE check_projector_density_floor_audit_response_helper

! **************************************************************************************************
!> \brief Check that production uses the gauge-safe diagonal dP source for physical density response.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_physical_density_safe_source(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 4
      REAL(KIND=dp), PARAMETER :: active_tol = 1.0E-4_dp, gauge_gap = 5.0E-5_dp, &
         projector_weight = 0.42_dp, response_ratio_gate = 2.0E-2_dp

      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: chc, dchc, docc_raw, docc_safe
      INTEGER                                            :: i
      REAL(KIND=dp)                                      :: collapsed, err, occ_curv_max, retained, &
                                                            safe_norm, spectral_response
      REAL(KIND=dp), DIMENSION(n)                        :: occ_curv, occ_scale, projector_response, &
                                                            raw_density, safe_density, shift

      occ_scale(:) = [1.58_dp, 0.64_dp, 0.31_dp, 0.11_dp]
      occ_curv(:) = [0.040_dp, 0.032_dp, 0.013_dp, 0.002_dp]
      shift(:) = [0.11_dp, -0.07_dp, 0.03_dp, -0.02_dp]
      occ_curv_max = MAXVAL(occ_curv(:))

      chc(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      dchc(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      DO i = 1, n
         chc(i, i) = CMPLX(0.08_dp*REAL(i, KIND=dp), 0.0_dp, KIND=dp)
         dchc(i, i) = CMPLX(shift(i), 0.0_dp, KIND=dp)
      END DO

      CALL dense_history_safe_occ_projector_derivative_chc(chc, dchc, occ_scale, projector_weight, docc_raw)
      CALL dense_history_safe_occ_projector_derivative_chc(chc, dchc, occ_scale, projector_weight, &
                                                           docc_safe, occ_curv, gauge_gap)

      DO i = 1, n
         projector_response(i) = -projector_weight*shift(i)
         raw_density(i) = 0.5_dp*REAL(docc_raw(i, i), KIND=dp)/occ_curv(i)
         safe_density(i) = 0.5_dp*REAL(docc_safe(i, i), KIND=dp)/occ_curv(i)
      END DO

      retained = dense_physical_density_response(projector_response(1), safe_density(1), &
                                                 occ_curv(1), occ_curv_max, response_ratio_gate)
      collapsed = dense_physical_density_response(projector_response(1), raw_density(1), &
                                                  occ_curv(1), occ_curv_max, response_ratio_gate)
      spectral_response = projector_response(1)
      safe_norm = MAXVAL(ABS(safe_density(:)))
      err = MAXVAL(ABS(raw_density(:)))

      IF (err > tol .OR. ABS(retained - safe_density(1)) > tol .OR. &
          ABS(collapsed - spectral_response) > tol .OR. safe_norm < active_tol .OR. &
          ABS(retained - spectral_response) < active_tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,5(1X,ES12.5))') TRIM(label), &
            "failed with raw/retained/collapsed/safe/spec", err, retained, collapsed, safe_norm, &
            spectral_response
      ELSE
         WRITE (*, '(A,1X,A,5(1X,ES12.5))') TRIM(label), &
            "OK with raw/retained/collapsed/safe/spec", err, retained, collapsed, safe_norm, &
            spectral_response
      END IF

   END SUBROUTINE check_projector_physical_density_safe_source

! **************************************************************************************************
!> \brief Check that the active density split preserves the finite physical response.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_active_density_total_invariant(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 4
      REAL(KIND=dp), PARAMETER                           :: active_tol = 1.0E-4_dp, &
                                                            fd_step = 1.0E-7_dp, &
                                                            gauge_gap = 5.0E-5_dp, &
                                                            projector_weight = 0.42_dp

      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: chc, chc_eig, dchc, dchc_eig, docc_safe, &
                                                            q, work
      INTEGER                                            :: i
      REAL(KIND=dp) :: active_fixed_n_leak, active_shift, active_sum, active_weight, audit_gap, &
         finite_total_err, full_fixed_n_leak, full_shift, full_sum, occ_curv_max, total_active_gap
      REAL(KIND=dp), DIMENSION(n) :: active_density_response, density_column, &
         finite_physical_response, full_density_response, occ_curv, occ_scale, occ_work, &
         physical_occ_response, projector_evals_m, projector_evals_p, total_response
      REAL(KIND=dp), DIMENSION(n, n)                     :: x_im, x_re

      chc_eig(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      chc_eig(1, 1) = CMPLX(-0.018000_dp, 0.0_dp, KIND=dp)
      chc_eig(2, 2) = CMPLX(-0.017986_dp, 0.0_dp, KIND=dp)
      chc_eig(3, 3) = CMPLX(0.19_dp, 0.0_dp, KIND=dp)
      chc_eig(4, 4) = CMPLX(0.47_dp, 0.0_dp, KIND=dp)

      x_re(:, :) = 0.0_dp
      x_re(1, 2) = 0.10_dp
      x_re(1, 3) = -0.07_dp
      x_re(2, 4) = 0.05_dp
      x_re(3, 4) = -0.04_dp
      DO i = 2, n
         x_re(i, 1:i - 1) = -x_re(1:i - 1, i)
      END DO
      x_im(:, :) = 0.0_dp
      x_im(1, 1) = 0.01_dp
      x_im(2, 2) = -0.02_dp
      x_im(3, 3) = 0.015_dp
      x_im(1, 2) = -0.03_dp
      x_im(1, 4) = 0.04_dp
      x_im(2, 3) = 0.02_dp
      DO i = 2, n
         x_im(i, 1:i - 1) = x_im(1:i - 1, i)
      END DO
      CALL dense_rotation_u(x_re, x_im, q)
      chc(:, :) = MATMUL(q, MATMUL(chc_eig, CONJG(TRANSPOSE(q))))

      occ_scale(:) = [1.68_dp, 0.48_dp, 0.25_dp, 0.17_dp]
      occ_curv(:) = [0.040_dp, 0.032_dp, 0.013_dp, 0.0017_dp]
      occ_curv_max = MAXVAL(occ_curv(:))

      dchc_eig(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      dchc_eig(1, 2) = CMPLX(0.08_dp, -0.05_dp, KIND=dp)
      dchc_eig(2, 1) = CONJG(dchc_eig(1, 2))
      dchc_eig(3, 4) = CMPLX(-0.14_dp, 0.09_dp, KIND=dp)
      dchc_eig(4, 3) = CONJG(dchc_eig(3, 4))
      dchc(:, :) = MATMUL(q, MATMUL(dchc_eig, CONJG(TRANSPOSE(q))))

      CALL dense_history_safe_occ_projector_derivative_chc(chc, dchc, occ_scale, projector_weight, &
                                                           docc_safe, occ_curv, gauge_gap)
      CALL diag_complex(chc + fd_step*dchc, work, projector_evals_p)
      CALL diag_complex(chc - fd_step*dchc, work, projector_evals_m)

      full_sum = SUM(occ_curv(:))
      full_shift = 0.0_dp
      active_sum = 0.0_dp
      active_shift = 0.0_dp
      DO i = 1, n
         finite_physical_response(i) = (1.0_dp - projector_weight)*REAL(dchc(i, i), KIND=dp) + &
                                       projector_weight*(projector_evals_p(i) - projector_evals_m(i))/(2.0_dp*fd_step)
         density_column(i) = 0.5_dp*REAL(docc_safe(i, i), KIND=dp)/occ_curv(i)
         full_shift = full_shift + occ_curv(i)*density_column(i)
         active_weight = SQRT(MIN(1.0_dp, MAX(0.0_dp, occ_curv(i)/occ_curv_max)))
         occ_work(i) = active_weight
         active_sum = active_sum + occ_curv(i)*active_weight
         active_shift = active_shift + occ_curv(i)*active_weight*density_column(i)
      END DO
      full_shift = full_shift/full_sum
      active_shift = active_shift/active_sum

      DO i = 1, n
         full_density_response(i) = density_column(i) - full_shift
         active_density_response(i) = occ_work(i)*(density_column(i) - active_shift)
         physical_occ_response(i) = finite_physical_response(i) - active_density_response(i)
      END DO
      total_response(:) = physical_occ_response(:) + active_density_response(:)

      finite_total_err = MAXVAL(ABS(total_response(:) - finite_physical_response(:)))
      active_fixed_n_leak = ABS(DOT_PRODUCT(occ_curv(:), active_density_response(:)))/full_sum
      full_fixed_n_leak = ABS(DOT_PRODUCT(occ_curv(:), full_density_response(:)))/full_sum
      audit_gap = MAXVAL(ABS(full_density_response(:) - active_density_response(:)))
      total_active_gap = MAXVAL(ABS(physical_occ_response(:) + full_density_response(:) - &
                                    finite_physical_response(:)))

      IF (finite_total_err > tol .OR. active_fixed_n_leak > tol .OR. full_fixed_n_leak > tol .OR. &
          audit_gap < active_tol .OR. total_active_gap < active_tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,5(1X,ES12.5))') TRIM(label), &
            "failed with total/afixed/ffixed/audit/wrongtotal", &
            finite_total_err, active_fixed_n_leak, full_fixed_n_leak, audit_gap, total_active_gap
      ELSE
         WRITE (*, '(A,1X,A,5(1X,ES12.5))') TRIM(label), &
            "OK with total/afixed/ffixed/audit/wrongtotal", &
            finite_total_err, active_fixed_n_leak, full_fixed_n_leak, audit_gap, total_active_gap
      END IF

   END SUBROUTINE check_projector_active_density_total_invariant

! **************************************************************************************************
!> \brief Check that raw projector-basis density motion stays in audit, not production response.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_complex_projector_density_raw_audit_split(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 4
      REAL(KIND=dp), PARAMETER                           :: fd_step = 1.0E-7_dp, &
                                                            gauge_gap = 5.0E-5_dp, &
                                                            projector_weight = 0.55_dp

      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: chc, chc_eig, dchc, dchc_eig, docc_raw, &
                                                            docc_safe, q, work
      INTEGER                                            :: i
      REAL(KIND=dp) :: active_norm, active_shift, active_sum, active_weight, audit_gap, err, &
         fixed_n_leak, occ_curv_max, occ_sum, production_norm, raw_norm, response_shift
      REAL(KIND=dp), DIMENSION(n) :: active_response, audit_response, occ_curv, occ_scale, &
         physical_column, production_full, projector_evals_m, projector_evals_p
      REAL(KIND=dp), DIMENSION(n, n)                     :: x_im, x_re

      chc_eig(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      chc_eig(1, 1) = CMPLX(-0.020000_dp, 0.0_dp, KIND=dp)
      chc_eig(2, 2) = CMPLX(-0.019987_dp, 0.0_dp, KIND=dp)
      chc_eig(3, 3) = CMPLX(0.17_dp, 0.0_dp, KIND=dp)
      chc_eig(4, 4) = CMPLX(0.44_dp, 0.0_dp, KIND=dp)

      x_re(:, :) = 0.0_dp
      x_re(1, 2) = 0.11_dp
      x_re(1, 3) = -0.05_dp
      x_re(2, 4) = 0.06_dp
      x_re(3, 4) = -0.04_dp
      DO i = 2, n
         x_re(i, 1:i - 1) = -x_re(1:i - 1, i)
      END DO
      x_im(:, :) = 0.0_dp
      x_im(1, 1) = 0.02_dp
      x_im(2, 2) = -0.01_dp
      x_im(1, 3) = -0.03_dp
      x_im(2, 4) = 0.04_dp
      DO i = 2, n
         x_im(i, 1:i - 1) = x_im(1:i - 1, i)
      END DO
      CALL dense_rotation_u(x_re, x_im, q)
      chc(:, :) = MATMUL(q, MATMUL(chc_eig, CONJG(TRANSPOSE(q))))

      occ_scale(:) = [1.76_dp, 0.36_dp, 0.19_dp, 0.08_dp]
      occ_curv(:) = [0.045_dp, 0.034_dp, 0.010_dp, 0.002_dp]
      dchc_eig(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      dchc_eig(1, 2) = CMPLX(0.19_dp, -0.07_dp, KIND=dp)
      dchc_eig(2, 1) = CONJG(dchc_eig(1, 2))
      dchc(:, :) = MATMUL(q, MATMUL(dchc_eig, CONJG(TRANSPOSE(q))))

      CALL dense_history_safe_occ_projector_derivative_chc(chc, dchc, occ_scale, projector_weight, &
                                                           docc_raw)
      CALL dense_history_safe_occ_projector_derivative_chc(chc, dchc, occ_scale, projector_weight, &
                                                           docc_safe, occ_curv, gauge_gap)
      CALL diag_complex(chc + fd_step*dchc, work, projector_evals_p)
      CALL diag_complex(chc - fd_step*dchc, work, projector_evals_m)

      DO i = 1, n
         audit_response(i) = 0.5_dp*REAL(docc_raw(i, i), KIND=dp)/occ_curv(i)
         physical_column(i) = -projector_weight* &
                              (projector_evals_p(i) - projector_evals_m(i))/(2.0_dp*fd_step)
      END DO
      occ_sum = SUM(occ_curv(:))
      response_shift = DOT_PRODUCT(occ_curv(:), physical_column(:))/occ_sum
      production_full(:) = physical_column(:) - response_shift

      occ_curv_max = MAXVAL(occ_curv(:))
      active_sum = 0.0_dp
      active_shift = 0.0_dp
      DO i = 1, n
         active_weight = SQRT(MIN(1.0_dp, MAX(0.0_dp, occ_curv(i)/occ_curv_max)))
         active_sum = active_sum + occ_curv(i)*active_weight
         active_shift = active_shift + occ_curv(i)*active_weight*physical_column(i)
      END DO
      IF (active_sum > EPSILON(active_sum)) active_shift = active_shift/active_sum
      DO i = 1, n
         active_weight = SQRT(MIN(1.0_dp, MAX(0.0_dp, occ_curv(i)/occ_curv_max)))
         active_response(i) = active_weight*(physical_column(i) - active_shift)
      END DO

      raw_norm = MAXVAL(ABS(audit_response(:)))
      production_norm = MAXVAL(ABS(production_full(:)))
      active_norm = MAXVAL(ABS(active_response(:)))
      audit_gap = MAXVAL(ABS(audit_response(:) - production_full(:)))
      fixed_n_leak = ABS(DOT_PRODUCT(occ_curv(:), production_full(:)))/occ_sum
      err = fixed_n_leak

      IF (err > tol .OR. raw_norm < 1.0E-2_dp .OR. &
          production_norm >= 0.25_dp*raw_norm .OR. active_norm > production_norm + tol .OR. &
          audit_gap < 0.5_dp*raw_norm) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,5(1X,ES12.5))') TRIM(label), &
            "failed with error/raw/prod/active/gap", err, raw_norm, production_norm, active_norm, audit_gap
      ELSE
         WRITE (*, '(A,1X,A,5(1X,ES12.5))') TRIM(label), &
            "OK with error/raw/prod/active/gap", err, raw_norm, production_norm, active_norm, audit_gap
      END IF

   END SUBROUTINE check_complex_projector_density_raw_audit_split

   ! **************************************************************************************************
   !> \brief Check that production density response follows finite eigenvalue response, not dP gauge diagonal.
   !> \param label ...
   !> \param tol ...
   !> \param nfail ...
   ! **************************************************************************************************
! **************************************************************************************************
!> \brief ...
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_complex_projector_density_eigen_response_split(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 4
      REAL(KIND=dp), PARAMETER                           :: active_tol = 1.0E-4_dp, &
                                                            fd_step = 1.0E-7_dp, &
                                                            gauge_gap = 5.0E-5_dp, &
                                                            projector_weight = 0.55_dp

      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: chc, chc_eig, dchc, dchc_eig, docc_safe, &
                                                            q, work
      INTEGER                                            :: i
      REAL(KIND=dp)                                      :: err, gauge_gap_norm, gauge_new_norm, &
                                                            gauge_old_norm, mu_shift, &
                                                            physical_err, physical_leak, &
                                                            physical_norm
      REAL(KIND=dp), DIMENSION(n) :: eig_response, expected_column, expected_full, new_column, &
         new_full, occ_curv, occ_scale, old_column, old_full, projector_evals_m, &
         projector_evals_p, shift
      REAL(KIND=dp), DIMENSION(n, n)                     :: x_im, x_re

      chc_eig(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      chc_eig(1, 1) = CMPLX(-0.020000_dp, 0.0_dp, KIND=dp)
      chc_eig(2, 2) = CMPLX(-0.019987_dp, 0.0_dp, KIND=dp)
      chc_eig(3, 3) = CMPLX(0.17_dp, 0.0_dp, KIND=dp)
      chc_eig(4, 4) = CMPLX(0.44_dp, 0.0_dp, KIND=dp)

      x_re(:, :) = 0.0_dp
      x_re(1, 2) = 0.11_dp
      x_re(1, 3) = -0.05_dp
      x_re(2, 4) = 0.06_dp
      x_re(3, 4) = -0.04_dp
      DO i = 2, n
         x_re(i, 1:i - 1) = -x_re(1:i - 1, i)
      END DO
      x_im(:, :) = 0.0_dp
      x_im(1, 1) = 0.02_dp
      x_im(2, 2) = -0.01_dp
      x_im(1, 3) = -0.03_dp
      x_im(2, 4) = 0.04_dp
      DO i = 2, n
         x_im(i, 1:i - 1) = x_im(1:i - 1, i)
      END DO
      CALL dense_rotation_u(x_re, x_im, q)
      chc(:, :) = MATMUL(q, MATMUL(chc_eig, CONJG(TRANSPOSE(q))))

      occ_scale(:) = [1.76_dp, 0.36_dp, 0.19_dp, 0.08_dp]
      occ_curv(:) = [0.045_dp, 0.034_dp, 0.010_dp, 0.002_dp]

      dchc_eig(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      dchc_eig(1, 2) = CMPLX(0.19_dp, -0.07_dp, KIND=dp)
      dchc_eig(2, 1) = CONJG(dchc_eig(1, 2))
      dchc(:, :) = MATMUL(q, MATMUL(dchc_eig, CONJG(TRANSPOSE(q))))
      CALL dense_history_safe_occ_projector_derivative_chc(chc, dchc, occ_scale, projector_weight, &
                                                           docc_safe, occ_curv, gauge_gap)
      CALL diag_complex(chc + fd_step*dchc, work, projector_evals_p)
      CALL diag_complex(chc - fd_step*dchc, work, projector_evals_m)
      eig_response(:) = (projector_evals_p(:) - projector_evals_m(:))/(2.0_dp*fd_step)
      new_column(:) = -projector_weight*eig_response(:)
      mu_shift = DOT_PRODUCT(occ_curv(:), new_column(:))/SUM(occ_curv(:))
      new_full(:) = new_column(:) - mu_shift
      DO i = 1, n
         old_column(i) = 0.5_dp*REAL(docc_safe(i, i), KIND=dp)/occ_curv(i)
      END DO
      mu_shift = DOT_PRODUCT(occ_curv(:), old_column(:))/SUM(occ_curv(:))
      old_full(:) = old_column(:) - mu_shift
      gauge_new_norm = MAXVAL(ABS(new_full(:)))
      gauge_old_norm = MAXVAL(ABS(old_full(:)))
      gauge_gap_norm = MAXVAL(ABS(old_full(:) - new_full(:)))

      shift(:) = [-0.002_dp, -0.002_dp, 0.030_dp, -0.015_dp]
      dchc_eig(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      DO i = 1, n
         dchc_eig(i, i) = CMPLX(shift(i), 0.0_dp, KIND=dp)
      END DO
      dchc(:, :) = MATMUL(q, MATMUL(dchc_eig, CONJG(TRANSPOSE(q))))
      CALL diag_complex(chc + fd_step*dchc, work, projector_evals_p)
      CALL diag_complex(chc - fd_step*dchc, work, projector_evals_m)
      eig_response(:) = (projector_evals_p(:) - projector_evals_m(:))/(2.0_dp*fd_step)
      new_column(:) = -projector_weight*eig_response(:)
      expected_column(:) = -projector_weight*shift(:)
      mu_shift = DOT_PRODUCT(occ_curv(:), new_column(:))/SUM(occ_curv(:))
      new_full(:) = new_column(:) - mu_shift
      mu_shift = DOT_PRODUCT(occ_curv(:), expected_column(:))/SUM(occ_curv(:))
      expected_full(:) = expected_column(:) - mu_shift
      physical_err = MAXVAL(ABS(new_full(:) - expected_full(:)))
      physical_leak = ABS(DOT_PRODUCT(occ_curv(:), new_full(:)))/SUM(occ_curv(:))
      physical_norm = MAXVAL(ABS(new_full(:)))
      err = MAX(gauge_new_norm, physical_err, physical_leak)

      IF (err > tol .OR. gauge_old_norm < active_tol .OR. gauge_gap_norm < active_tol .OR. &
          physical_norm < active_tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,7(1X,ES12.5))') TRIM(label), &
            "failed with error/gnew/gold/gap/perr/pleak/pnorm", err, gauge_new_norm, &
            gauge_old_norm, gauge_gap_norm, physical_err, physical_leak, physical_norm
      ELSE
         WRITE (*, '(A,1X,A,7(1X,ES12.5))') TRIM(label), &
            "OK with error/gnew/gold/gap/perr/pleak/pnorm", err, gauge_new_norm, &
            gauge_old_norm, gauge_gap_norm, physical_err, physical_leak, physical_norm
      END IF

   END SUBROUTINE check_complex_projector_density_eigen_response_split

! **************************************************************************************************
!> \brief Check that exact-cluster density response survives sorted-eigenvalue FD collapse.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_complex_projector_density_cluster_split_response(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 3
      REAL(KIND=dp), PARAMETER                           :: active_tol = 1.0E-4_dp, &
                                                            collapse_tol = 1.0E-6_dp, &
                                                            fd_step = 1.0E-7_dp, &
                                                            gauge_gap = 5.0E-5_dp, &
                                                            projector_weight = 0.55_dp, &
                                                            response_ratio_gate = 2.0E-2_dp

      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: chc, dchc, docc_safe, work
      INTEGER                                            :: i
      LOGICAL                                            :: fd_used
      REAL(KIND=dp)                                      :: alignment, collapse_ratio, err, &
                                                            filtered_err, occ_curv_max, occ_sum, &
                                                            projector_norm, safe_norm
      REAL(KIND=dp), DIMENSION(n) :: density_projector, density_projector_centered, &
         density_safe, density_safe_centered, eig_response, filtered_response, occ_curv, occ_scale, &
         physical_response, projector_evals_m, projector_evals_p

      chc(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      dchc(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      chc(1, 1) = CMPLX(-0.010_dp, 0.0_dp, KIND=dp)
      chc(2, 2) = CMPLX(-0.010_dp, 0.0_dp, KIND=dp)
      chc(3, 3) = CMPLX(0.420_dp, 0.0_dp, KIND=dp)

      dchc(1, 1) = CMPLX(0.080_dp, 0.0_dp, KIND=dp)
      dchc(2, 2) = CMPLX(-0.080_dp, 0.0_dp, KIND=dp)
      dchc(1, 2) = CMPLX(0.060_dp, -0.040_dp, KIND=dp)
      dchc(2, 1) = CONJG(dchc(1, 2))

      occ_scale(:) = [1.10_dp, 0.40_dp, 0.05_dp]
      occ_curv(:) = [0.045_dp, 0.045_dp, 0.004_dp]
      occ_curv_max = MAXVAL(occ_curv(:))

      CALL dense_history_safe_occ_projector_derivative_chc(chc, dchc, occ_scale, projector_weight, &
                                                           docc_safe, occ_curv, gauge_gap, fd_used)
      CALL diag_complex(chc + fd_step*dchc, work, projector_evals_p)
      CALL diag_complex(chc - fd_step*dchc, work, projector_evals_m)
      eig_response(:) = (projector_evals_p(:) - projector_evals_m(:))/(2.0_dp*fd_step)
      density_projector(:) = -projector_weight*eig_response(:)
      density_safe(:) = 0.0_dp
      DO i = 1, n
         IF (occ_curv(i) > EPSILON(occ_curv(i))) &
            density_safe(i) = 0.5_dp*REAL(docc_safe(i, i), KIND=dp)/occ_curv(i)
      END DO

      density_projector_centered(:) = density_projector(:)
      density_safe_centered(:) = density_safe(:)
      occ_sum = SUM(MAX(0.0_dp, occ_curv(:)))
      IF (occ_sum > EPSILON(occ_sum)) THEN
         density_projector_centered(:) = density_projector_centered(:) - &
            DOT_PRODUCT(MAX(0.0_dp, occ_curv(:)), density_projector_centered(:))/occ_sum
         density_safe_centered(:) = density_safe_centered(:) - &
            DOT_PRODUCT(MAX(0.0_dp, occ_curv(:)), density_safe_centered(:))/occ_sum
      END IF

      projector_norm = MAXVAL(ABS(density_projector_centered(:)))
      safe_norm = MAXVAL(ABS(density_safe_centered(:)))
      collapse_ratio = projector_norm/MAX(safe_norm, EPSILON(safe_norm))
      alignment = 0.0_dp
      IF (projector_norm > EPSILON(projector_norm) .AND. safe_norm > EPSILON(safe_norm)) &
         alignment = DOT_PRODUCT(MAX(0.0_dp, occ_curv(:))*density_projector_centered(:), &
                                 density_safe_centered(:))/ &
                     SQRT(DOT_PRODUCT(MAX(0.0_dp, occ_curv(:)), density_projector_centered(:)**2)* &
                          DOT_PRODUCT(MAX(0.0_dp, occ_curv(:)), density_safe_centered(:)**2))

      DO i = 1, n
         physical_response(i) = dense_physical_density_response(density_projector(i), density_safe(i), &
                                                                occ_curv(i), occ_curv_max, &
                                                                response_ratio_gate, &
                                                                density_projector_centered(i), &
                                                                density_safe_centered(i), alignment)
         filtered_response(i) = residual_gauge_density_response(physical_response(i), &
                                                                density_projector(i), density_safe(i), &
                                                                occ_curv(i), occ_curv_max, &
                                                                response_ratio_gate, 0.0_dp, &
                                                                density_projector_centered(i), &
                                                                density_safe_centered(i), alignment)
      END DO
      err = MAXVAL(ABS(physical_response(1:2) - density_safe_centered(1:2)))
      filtered_err = MAXVAL(ABS(filtered_response(1:2) - physical_response(1:2)))

      IF (.NOT. fd_used .OR. collapse_ratio > collapse_tol .OR. safe_norm < active_tol .OR. &
          err > tol .OR. filtered_err > tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,5(1X,ES12.5),1(1X,L1))') TRIM(label), &
            "failed with collapse/safe/err/ferr/align/fd", collapse_ratio, safe_norm, err, &
            filtered_err, alignment, fd_used
      ELSE
         WRITE (*, '(A,1X,A,5(1X,ES12.5),1(1X,L1))') TRIM(label), &
            "OK with collapse/safe/err/ferr/align/fd", collapse_ratio, safe_norm, err, &
            filtered_err, alignment, fd_used
      END IF

   END SUBROUTINE check_complex_projector_density_cluster_split_response

! **************************************************************************************************
!> \brief Check that near-cluster off-diagonal occupation response survives in a rotated REF basis.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_complex_projector_density_near_cluster_response(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 3
      REAL(KIND=dp), PARAMETER                           :: active_tol = 1.0E-4_dp, &
                                                            collapse_tol = 1.0E-6_dp, &
                                                            fd_step = 1.0E-7_dp, &
                                                            gauge_gap = 5.0E-5_dp, &
                                                            projector_weight = 0.55_dp, &
                                                            response_ratio_gate = 2.0E-2_dp

      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: chc, chc_eig, dchc, dchc_eig, &
                                                            docc_safe, q, work
      INTEGER                                            :: i
      LOGICAL                                            :: fd_used
      REAL(KIND=dp)                                      :: alignment, collapse_ratio, err, &
                                                            filtered_err, occ_curv_max, occ_sum, &
                                                            projector_norm, safe_norm, unit
      REAL(KIND=dp), DIMENSION(n) :: density_projector, density_projector_centered, &
         density_safe, density_safe_centered, eig_response, filtered_response, occ_curv, occ_scale, &
         physical_response, projector_evals_m, projector_evals_p

      unit = 1.0_dp/SQRT(2.0_dp)
      q(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      q(1, 1) = CMPLX(unit, 0.0_dp, KIND=dp)
      q(2, 1) = CMPLX(unit, 0.0_dp, KIND=dp)
      q(1, 2) = CMPLX(unit, 0.0_dp, KIND=dp)
      q(2, 2) = CMPLX(-unit, 0.0_dp, KIND=dp)
      q(3, 3) = CMPLX(1.0_dp, 0.0_dp, KIND=dp)

      chc_eig(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      dchc_eig(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      chc_eig(1, 1) = CMPLX(-0.010000_dp, 0.0_dp, KIND=dp)
      chc_eig(2, 2) = CMPLX(-0.009982_dp, 0.0_dp, KIND=dp)
      chc_eig(3, 3) = CMPLX(0.420_dp, 0.0_dp, KIND=dp)
      dchc_eig(1, 2) = CMPLX(0.060_dp, -0.040_dp, KIND=dp)
      dchc_eig(2, 1) = CONJG(dchc_eig(1, 2))
      chc(:, :) = MATMUL(q, MATMUL(chc_eig, CONJG(TRANSPOSE(q))))
      dchc(:, :) = MATMUL(q, MATMUL(dchc_eig, CONJG(TRANSPOSE(q))))

      occ_scale(:) = [1.10_dp, 0.40_dp, 0.05_dp]
      occ_curv(:) = [0.045_dp, 0.045_dp, 0.004_dp]
      occ_curv_max = MAXVAL(occ_curv(:))

      CALL dense_history_safe_occ_projector_derivative_chc(chc, dchc, occ_scale, projector_weight, &
                                                           docc_safe, occ_curv, gauge_gap, fd_used)
      CALL diag_complex(chc + fd_step*dchc, work, projector_evals_p)
      CALL diag_complex(chc - fd_step*dchc, work, projector_evals_m)
      eig_response(:) = (projector_evals_p(:) - projector_evals_m(:))/(2.0_dp*fd_step)
      density_projector(:) = -projector_weight*eig_response(:)
      density_safe(:) = 0.0_dp
      DO i = 1, n
         IF (occ_curv(i) > EPSILON(occ_curv(i))) &
            density_safe(i) = 0.5_dp*REAL(docc_safe(i, i), KIND=dp)/occ_curv(i)
      END DO

      density_projector_centered(:) = density_projector(:)
      density_safe_centered(:) = density_safe(:)
      occ_sum = SUM(MAX(0.0_dp, occ_curv(:)))
      IF (occ_sum > EPSILON(occ_sum)) THEN
         density_projector_centered(:) = density_projector_centered(:) - &
            DOT_PRODUCT(MAX(0.0_dp, occ_curv(:)), density_projector_centered(:))/occ_sum
         density_safe_centered(:) = density_safe_centered(:) - &
            DOT_PRODUCT(MAX(0.0_dp, occ_curv(:)), density_safe_centered(:))/occ_sum
      END IF

      projector_norm = MAXVAL(ABS(density_projector_centered(:)))
      safe_norm = MAXVAL(ABS(density_safe_centered(:)))
      collapse_ratio = projector_norm/MAX(safe_norm, EPSILON(safe_norm))
      alignment = 0.0_dp
      IF (projector_norm > EPSILON(projector_norm) .AND. safe_norm > EPSILON(safe_norm)) &
         alignment = DOT_PRODUCT(MAX(0.0_dp, occ_curv(:))*density_projector_centered(:), &
                                 density_safe_centered(:))/ &
                     SQRT(DOT_PRODUCT(MAX(0.0_dp, occ_curv(:)), density_projector_centered(:)**2)* &
                          DOT_PRODUCT(MAX(0.0_dp, occ_curv(:)), density_safe_centered(:)**2))

      DO i = 1, n
         physical_response(i) = dense_physical_density_response(density_projector(i), density_safe(i), &
                                                                occ_curv(i), occ_curv_max, &
                                                                response_ratio_gate, &
                                                                density_projector_centered(i), &
                                                                density_safe_centered(i), alignment)
         filtered_response(i) = residual_gauge_density_response(physical_response(i), &
                                                                density_projector(i), density_safe(i), &
                                                                occ_curv(i), occ_curv_max, &
                                                                response_ratio_gate, 0.0_dp, &
                                                                density_projector_centered(i), &
                                                                density_safe_centered(i), alignment)
      END DO
      err = MAXVAL(ABS(physical_response(1:2) - density_safe_centered(1:2)))
      filtered_err = MAXVAL(ABS(filtered_response(1:2) - physical_response(1:2)))

      IF (.NOT. fd_used .OR. collapse_ratio > collapse_tol .OR. safe_norm < active_tol .OR. &
          err > tol .OR. filtered_err > tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,5(1X,ES12.5),1(1X,L1))') TRIM(label), &
            "failed with collapse/safe/err/ferr/align/fd", collapse_ratio, safe_norm, err, &
            filtered_err, alignment, fd_used
      ELSE
         WRITE (*, '(A,1X,A,5(1X,ES12.5),1(1X,L1))') TRIM(label), &
            "OK with collapse/safe/err/ferr/align/fd", collapse_ratio, safe_norm, err, &
            filtered_err, alignment, fd_used
      END IF

   END SUBROUTINE check_complex_projector_density_near_cluster_response

! **************************************************************************************************
!> \brief Check that local-mean-matched dP diagonals are not a valid density-response replacement.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_density_local_mean_dP_rejection(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: nchannels = 2, ne = 6

      INTEGER                                            :: first_energy, local_channel, nmo
      INTEGER, DIMENSION(nchannels)                      :: channel_nmo, e_offset
      LOGICAL                                            :: dP_projected, eigen_projected
      REAL(KIND=dp)                                      :: charge_dP, charge_eigen, err, &
                                                            local_mean_err, metric_gap, shape_gap, &
                                                            shift, weight_sum
      REAL(KIND=dp), DIMENSION(ne) :: dP_projected_response, eigen_projected_response, &
         eigen_response, local_mean_dP, occ_curv, raw_dP, response_gap, tether_curv

      channel_nmo(:) = [3, 3]
      e_offset(:) = [1, 4]
      occ_curv(:) = [0.045_dp, 0.030_dp, 0.006_dp, 0.026_dp, 0.015_dp, 0.004_dp]
      tether_curv(:) = 0.0_dp
      eigen_response(:) = [0.10_dp, -0.05_dp, 0.02_dp, -0.12_dp, 0.04_dp, 0.03_dp]
      raw_dP(:) = [1.70_dp, -0.90_dp, 0.25_dp, -1.10_dp, 0.80_dp, -0.45_dp]

      ! This mirrors the tempting but harmful replacement: use the back-transformed
      ! dP diagonal shape and only restore each local channel's weighted mean.
      local_mean_dP(:) = raw_dP(:)
      local_mean_err = 0.0_dp
      DO local_channel = 1, nchannels
         nmo = channel_nmo(local_channel)
         first_energy = e_offset(local_channel)
         weight_sum = SUM(occ_curv(first_energy:first_energy + nmo - 1))
         shift = (DOT_PRODUCT(occ_curv(first_energy:first_energy + nmo - 1), &
                              eigen_response(first_energy:first_energy + nmo - 1)) - &
                  DOT_PRODUCT(occ_curv(first_energy:first_energy + nmo - 1), &
                              local_mean_dP(first_energy:first_energy + nmo - 1)))/weight_sum
         local_mean_dP(first_energy:first_energy + nmo - 1) = &
            local_mean_dP(first_energy:first_energy + nmo - 1) + shift
         local_mean_err = MAX(local_mean_err, &
                              ABS(DOT_PRODUCT(occ_curv(first_energy:first_energy + nmo - 1), &
                                              local_mean_dP(first_energy:first_energy + nmo - 1) - &
                                              eigen_response(first_energy:first_energy + nmo - 1))))
      END DO

      eigen_projected_response(:) = eigen_response(:)
      CALL dense_project_response_metric(eigen_projected_response, occ_curv, tether_curv, &
                                         channel_nmo, e_offset, eigen_projected)
      dP_projected_response(:) = local_mean_dP(:)
      CALL dense_project_response_metric(dP_projected_response, occ_curv, tether_curv, &
                                         channel_nmo, e_offset, dP_projected)

      response_gap(:) = dP_projected_response(:) - eigen_projected_response(:)
      CALL dense_projected_mixed_metric(response_gap, response_gap, occ_curv, tether_curv, &
                                        channel_nmo, e_offset, metric_gap)
      shape_gap = MAXVAL(ABS(response_gap(:)))
      charge_eigen = ABS(DOT_PRODUCT(occ_curv(:), eigen_projected_response(:)))
      charge_dP = ABS(DOT_PRODUCT(occ_curv(:), dP_projected_response(:)))
      err = MAX(local_mean_err, charge_eigen, charge_dP)

      IF (.NOT. eigen_projected .OR. .NOT. dP_projected .OR. err > tol .OR. &
          shape_gap < 0.50_dp .OR. metric_gap < 1.0E-2_dp) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,5(1X,ES12.5))') TRIM(label), &
            "failed with error/shape/metric/qe/qd", err, shape_gap, metric_gap, charge_eigen, charge_dP
      ELSE
         WRITE (*, '(A,1X,A,5(1X,ES12.5))') TRIM(label), &
            "OK with error/shape/metric/qe/qd", err, shape_gap, metric_gap, charge_eigen, charge_dP
      END IF

   END SUBROUTINE check_projector_density_local_mean_dP_rejection

! **************************************************************************************************
!> \brief Check that plateau trust removes local gauge shifts but full response keeps global transfer.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_density_plateau_trust_centering(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: nchannels = 2, ne = 6

      INTEGER                                            :: first_energy, local_channel, nmo
      INTEGER, DIMENSION(nchannels)                      :: channel_nmo, e_offset
      LOGICAL                                            :: full_projected
      REAL(KIND=dp)                                      :: err, full_metric, local_metric, &
                                                            raw_local_norm, shift, weight_sum
      REAL(KIND=dp), DIMENSION(ne)                       :: full_global, occ_curv, plateau_trust, &
                                                            projector_response, tether_curv

      channel_nmo(:) = [3, 3]
      e_offset(:) = [1, 4]
      occ_curv(:) = [0.045_dp, 0.030_dp, 0.006_dp, 0.026_dp, 0.015_dp, 0.004_dp]
      tether_curv(:) = 0.0_dp
      projector_response(:) = [0.31_dp, 0.31_dp, 0.31_dp, -0.18_dp, -0.18_dp, -0.18_dp]

      plateau_trust(:) = projector_response(:)
      raw_local_norm = 0.0_dp
      DO local_channel = 1, nchannels
         nmo = channel_nmo(local_channel)
         first_energy = e_offset(local_channel)
         weight_sum = SUM(occ_curv(first_energy:first_energy + nmo - 1))
         shift = DOT_PRODUCT(occ_curv(first_energy:first_energy + nmo - 1), &
                             plateau_trust(first_energy:first_energy + nmo - 1))/weight_sum
         raw_local_norm = MAX(raw_local_norm, ABS(shift))
         plateau_trust(first_energy:first_energy + nmo - 1) = &
            plateau_trust(first_energy:first_energy + nmo - 1) - shift
      END DO

      full_global(:) = projector_response(:)
      CALL dense_project_response_metric(full_global, occ_curv, tether_curv, &
                                         channel_nmo, e_offset, full_projected)
      CALL dense_projected_mixed_metric(plateau_trust, plateau_trust, occ_curv, tether_curv, &
                                        channel_nmo, e_offset, local_metric)
      CALL dense_projected_mixed_metric(full_global, full_global, occ_curv, tether_curv, &
                                        channel_nmo, e_offset, full_metric)
      err = MAXVAL(ABS(plateau_trust(:)))
      err = MAX(err, ABS(DOT_PRODUCT(occ_curv(:), full_global(:))))

      IF (.NOT. full_projected .OR. err > tol .OR. raw_local_norm < 1.0E-1_dp .OR. &
          local_metric > tol .OR. full_metric < 1.0E-3_dp) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "failed with error/local/full/raw", err, local_metric, full_metric, raw_local_norm
      ELSE
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "OK with error/local/full/raw", err, local_metric, full_metric, raw_local_norm
      END IF

   END SUBROUTINE check_projector_density_plateau_trust_centering

! **************************************************************************************************
!> \brief Check that gauge-free eigenvalue density response is floor-limited but not active-weight damped.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_complex_projector_density_eigen_floor_filter(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 4
      REAL(KIND=dp), PARAMETER :: active_tol = 1.0E-4_dp, response_ratio_gate = 2.0E-2_dp

      INTEGER                                            :: i
      REAL(KIND=dp) :: active_floor, fixed_n_leak, full_norm, full_shift, mix_weight, new_shift, &
         new_sum, old_shift, old_sum, recovery, weight
      REAL(KIND=dp), DIMENSION(n)                        :: active_column, eig_column, &
                                                            full_response, mixed_column, &
                                                            new_response, occ_curv, old_response

      occ_curv(:) = [4.0E-2_dp, 2.0E-2_dp, 4.0E-4_dp, 4.0E-5_dp]
      eig_column(:) = [3.0E-3_dp, -2.0E-3_dp, 1.5E-3_dp, -1.0E-3_dp]
      active_floor = response_ratio_gate*MAXVAL(occ_curv(:))

      active_column(:) = 0.0_dp
      mixed_column(:) = 0.0_dp
      DO i = 1, n
         active_column(i) = eig_column(i)*occ_curv(i)/MAX(occ_curv(i), active_floor)
         mix_weight = occ_curv(i)/MAX(occ_curv(i), active_floor)
         mix_weight = MIN(1.0_dp, MAX(0.0_dp, mix_weight))
         mixed_column(i) = active_column(i) + mix_weight*(eig_column(i) - active_column(i))
      END DO

      full_shift = DOT_PRODUCT(occ_curv(:), eig_column(:))/SUM(occ_curv(:))
      full_response(:) = eig_column(:) - full_shift
      full_norm = MAXVAL(ABS(full_response(:)))

      new_sum = SUM(occ_curv(:))
      new_shift = DOT_PRODUCT(occ_curv(:), mixed_column(:))/new_sum
      new_response(:) = mixed_column(:) - new_shift
      fixed_n_leak = ABS(DOT_PRODUCT(occ_curv(:), new_response(:)))/new_sum

      old_sum = 0.0_dp
      old_shift = 0.0_dp
      DO i = 1, n
         weight = SQRT(MIN(1.0_dp, MAX(0.0_dp, occ_curv(i)/MAXVAL(occ_curv(:)))))
         old_sum = old_sum + occ_curv(i)*weight
         old_shift = old_shift + occ_curv(i)*weight*mixed_column(i)
      END DO
      old_shift = old_shift/old_sum
      DO i = 1, n
         weight = SQRT(MIN(1.0_dp, MAX(0.0_dp, occ_curv(i)/MAXVAL(occ_curv(:)))))
         old_response(i) = weight*(mixed_column(i) - old_shift)
      END DO

      recovery = MAXVAL(ABS(new_response(3:4)))/MAX(MAXVAL(ABS(old_response(3:4))), EPSILON(1.0_dp))

      IF (fixed_n_leak > tol .OR. full_norm < active_tol .OR. &
          MAXVAL(ABS(new_response(:))) > full_norm + active_tol .OR. recovery < 2.0_dp) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "failed with fixed/full/new/recovery", fixed_n_leak, full_norm, &
            MAXVAL(ABS(new_response(:))), recovery
      ELSE
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "OK with fixed/full/new/recovery", fixed_n_leak, full_norm, &
            MAXVAL(ABS(new_response(:))), recovery
      END IF

   END SUBROUTINE check_complex_projector_density_eigen_floor_filter

   ! **************************************************************************************************
   !> \brief Check that active rotation weights require both bands to be in the Fermi window.
   !> \param label ...
   !> \param tol ...
   !> \param nfail ...
   ! **************************************************************************************************
! **************************************************************************************************
!> \brief ...
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_active_rotation_pair_weight(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      REAL(KIND=dp), PARAMETER                           :: curv_max = 0.040_dp, &
                                                            projector_weight = 0.60_dp

      REAL(KIND=dp)                                      :: active_inert, diag_delta_ratio, err, &
                                                            full_active, half_active, &
                                                            old_active_inert, reversed, zero_weight

      full_active = dense_pair_projector_weight(projector_weight, curv_max, curv_max, curv_max)
      half_active = dense_pair_projector_weight(projector_weight, 0.25_dp*curv_max, &
                                                0.25_dp*curv_max, curv_max)
      active_inert = dense_pair_projector_weight(projector_weight, curv_max, &
                                                 1.0E-4_dp*curv_max, curv_max)
      reversed = dense_pair_projector_weight(projector_weight, 1.0E-4_dp*curv_max, &
                                             curv_max, curv_max)
      zero_weight = dense_pair_projector_weight(projector_weight, curv_max, 0.0_dp, curv_max)
      old_active_inert = projector_weight
      diag_delta_ratio = active_inert/MAX(old_active_inert, EPSILON(old_active_inert))

      err = ABS(full_active - projector_weight)
      err = MAX(err, ABS(half_active - 0.25_dp*projector_weight))
      err = MAX(err, ABS(active_inert - reversed), ABS(zero_weight))

      IF (err > tol .OR. diag_delta_ratio >= 0.05_dp) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,6(1X,ES12.5))') TRIM(label), &
            "failed with error/full/half/active-inert/zero/diagratio", err, full_active, half_active, &
            active_inert, zero_weight, diag_delta_ratio
      ELSE
         WRITE (*, '(A,1X,A,6(1X,ES12.5))') TRIM(label), &
            "OK with error/full/half/active-inert/zero/diagratio", err, full_active, half_active, &
            active_inert, zero_weight, diag_delta_ratio
      END IF

   END SUBROUTINE check_projector_active_rotation_pair_weight

! **************************************************************************************************
!> \brief Check that trusted density history can lift physical gap pairs without opening gauge pairs.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_physical_pair_projector_weight(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      REAL(KIND=dp), PARAMETER                           :: projector_weight = 0.60_dp

      REAL(KIND=dp)                                      :: err, fermi_weight, full_gap, gauge_pair, &
                                                            half_gap, low_history_gap, no_history, &
                                                            over_history, response_lifted_pair, &
                                                            response_lifted_no_history, trusted_gap

      fermi_weight = 0.02_dp*projector_weight
      trusted_gap = dense_physical_pair_projector_weight(projector_weight, fermi_weight, &
                                                         0.80_dp, 0.50_dp)
      half_gap = dense_physical_pair_projector_weight(projector_weight, fermi_weight, &
                                                      0.50_dp, 0.50_dp)
      full_gap = dense_physical_pair_projector_weight(projector_weight, fermi_weight, &
                                                      1.00_dp, 0.50_dp)
      low_history_gap = dense_physical_pair_projector_weight(projector_weight, fermi_weight, &
                                                             1.00_dp, 0.10_dp)
      gauge_pair = dense_physical_pair_projector_weight(projector_weight, fermi_weight, &
                                                        0.00_dp, 0.50_dp)
      no_history = dense_physical_pair_projector_weight(projector_weight, fermi_weight, &
                                                        0.80_dp, 0.00_dp)
      over_history = dense_physical_pair_projector_weight(projector_weight, fermi_weight, &
                                                          1.00_dp, 2.00_dp)
      response_lifted_pair = dense_physical_pair_projector_weight(projector_weight, fermi_weight, &
                                                                  0.00_dp, 0.50_dp, 0.70_dp)
      response_lifted_no_history = dense_physical_pair_projector_weight(projector_weight, fermi_weight, &
                                                                        0.00_dp, 0.00_dp, 1.00_dp)

      err = ABS(trusted_gap - 0.80_dp*SQRT(SQRT(0.50_dp))*projector_weight)
      err = MAX(err, ABS(half_gap - 0.50_dp*SQRT(SQRT(0.50_dp))*projector_weight))
      err = MAX(err, ABS(full_gap - SQRT(SQRT(0.50_dp))*projector_weight))
      err = MAX(err, ABS(low_history_gap - SQRT(SQRT(0.10_dp))*projector_weight))
      err = MAX(err, ABS(gauge_pair - fermi_weight))
      err = MAX(err, ABS(no_history - fermi_weight))
      err = MAX(err, ABS(over_history - projector_weight))
      err = MAX(err, ABS(response_lifted_pair - 0.70_dp*SQRT(SQRT(0.50_dp))*projector_weight))
      err = MAX(err, ABS(response_lifted_no_history - fermi_weight))

      IF (err > tol .OR. trusted_gap <= fermi_weight .OR. gauge_pair > trusted_gap .OR. &
          no_history > trusted_gap .OR. low_history_gap <= 0.10_dp*projector_weight .OR. &
          response_lifted_pair <= fermi_weight .OR. response_lifted_no_history > fermi_weight) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,9(1X,ES12.5))') TRIM(label), &
            "failed with error/trusted/half/full/lowhist/gauge/nohist/resp/nresp", err, trusted_gap, &
            half_gap, full_gap, low_history_gap, gauge_pair, no_history, response_lifted_pair, &
            response_lifted_no_history
      ELSE
         WRITE (*, '(A,1X,A,9(1X,ES12.5))') TRIM(label), &
            "OK with error/trusted/half/full/lowhist/gauge/nohist/resp/nresp", err, trusted_gap, &
            half_gap, full_gap, low_history_gap, gauge_pair, no_history, response_lifted_pair, &
            response_lifted_no_history
      END IF

   END SUBROUTINE check_projector_physical_pair_projector_weight

! **************************************************************************************************
!> \brief Check that density-gap response keeps only fixed-N Mermin descent components.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_density_descent_gap_projection(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: ne = 5

      REAL(KIND=dp)                                      :: coeff, denom, err, gdot, neg_norm, &
                                                            null_norm, phys_norm, raw_norm
      REAL(KIND=dp), DIMENSION(ne)                       :: descent, gap_null, gap_phys, grad, &
                                                            projected

      grad(:) = [0.31_dp, -0.24_dp, 0.17_dp, -0.09_dp, 0.05_dp]
      descent(:) = [0.42_dp, -0.18_dp, 0.11_dp, -0.04_dp, 0.02_dp]
      gap_null(:) = [0.80_dp, -0.10_dp, -0.30_dp, 0.45_dp, -0.20_dp]
      gap_null(:) = gap_null(:) - DOT_PRODUCT(grad, gap_null)/DOT_PRODUCT(grad, grad)*grad(:)

      CALL dense_project_descent_response(gap_null, descent, grad, projected, coeff, gdot, denom)
      null_norm = MAXVAL(ABS(gap_null(:)))
      err = MAXVAL(ABS(projected(:)))
      err = MAX(err, ABS(gdot))
      err = MAX(err, ABS(coeff))

      gap_phys(:) = gap_null(:) + 0.23_dp*descent(:)
      CALL dense_project_descent_response(gap_phys, descent, grad, projected, coeff, gdot, denom)
      raw_norm = MAXVAL(ABS(gap_phys(:)))
      phys_norm = MAXVAL(ABS(projected(:)))
      err = MAX(err, MAXVAL(ABS(projected(:) - 0.23_dp*descent(:))))
      err = MAX(err, ABS(coeff - 0.23_dp))
      err = MAX(err, ABS(gdot - 0.23_dp*denom))

      gap_phys(:) = gap_null(:) - 0.17_dp*descent(:)
      CALL dense_project_descent_response(gap_phys, descent, grad, projected, coeff, gdot, denom)
      err = MAX(err, MAXVAL(ABS(projected(:))))
      err = MAX(err, ABS(coeff))
      err = MAX(err, ABS(gdot + 0.17_dp*denom))
      CALL dense_project_descent_response(gap_phys, descent, grad, projected, coeff, gdot, denom, &
                                          flip_negative=.TRUE.)
      neg_norm = MAXVAL(ABS(projected(:)))
      err = MAX(err, MAXVAL(ABS(projected(:) - 0.17_dp*descent(:))))
      err = MAX(err, ABS(coeff - 0.17_dp))
      err = MAX(err, ABS(gdot + 0.17_dp*denom))

      IF (err > tol .OR. null_norm < 1.0E-1_dp .OR. phys_norm >= raw_norm .OR. neg_norm <= 0.0_dp) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,5(1X,ES12.5))') TRIM(label), &
            "failed with error/null/raw/phys/neg", err, null_norm, raw_norm, phys_norm, neg_norm
      ELSE
         WRITE (*, '(A,1X,A,5(1X,ES12.5))') TRIM(label), &
            "OK with error/null/raw/phys/neg", err, null_norm, raw_norm, phys_norm, neg_norm
      END IF

   END SUBROUTINE check_projector_density_descent_gap_projection

! **************************************************************************************************
!> \brief Check that the active density filter can leave a recoverable Mermin-descent gap.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_density_active_filter_descent_recovery(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: ne = 5
      REAL(KIND=dp), PARAMETER                           :: expected_coeff = 0.31_dp

      INTEGER                                            :: i
      REAL(KIND=dp)                                      :: coeff, denom, err, gap_weighted, gdot, &
                                                            projected_weighted, &
                                                            recoverable_fraction
      REAL(KIND=dp), DIMENSION(ne)                       :: active_response, active_weight, descent, &
                                                            full_response, gap, gap_null, grad, &
                                                            occ_curv, projected

      occ_curv(:) = [0.040_dp, 0.020_dp, 0.010_dp, 0.005_dp, 0.002_dp]
      active_weight(:) = SQRT(occ_curv(:)/MAXVAL(occ_curv(:)))
      grad(:) = [0.0_dp, 0.40_dp, -0.20_dp, 0.10_dp, -0.05_dp]
      descent(:) = [0.0_dp, 0.20_dp, -0.10_dp, 0.04_dp, -0.02_dp]

      gap_null(:) = [0.0_dp, 0.10_dp, 0.20_dp, -0.10_dp, 0.05_dp]
      gap_null(:) = gap_null(:) - DOT_PRODUCT(grad, gap_null)/DOT_PRODUCT(grad, grad)*grad(:)
      gap(:) = gap_null(:) + expected_coeff*descent(:)

      full_response(:) = 0.0_dp
      DO i = 1, ne
         IF (1.0_dp - active_weight(i) > EPSILON(active_weight(i))) &
            full_response(i) = gap(i)/(1.0_dp - active_weight(i))
      END DO
      active_response(:) = active_weight(:)*full_response(:)
      gap(:) = full_response(:) - active_response(:)

      CALL dense_project_descent_response(gap, descent, grad, projected, coeff, gdot, denom)
      gap_weighted = SQRT(MAX(0.0_dp, DOT_PRODUCT(occ_curv(:)*gap(:), gap(:))/SUM(occ_curv(:))))
      projected_weighted = SQRT(MAX(0.0_dp, DOT_PRODUCT(occ_curv(:)*projected(:), projected(:))/SUM(occ_curv(:))))
      recoverable_fraction = 0.0_dp
      IF (gap_weighted > EPSILON(gap_weighted)) recoverable_fraction = projected_weighted/gap_weighted

      err = MAXVAL(ABS(projected(:) - expected_coeff*descent(:)))
      err = MAX(err, ABS(coeff - expected_coeff))
      err = MAX(err, ABS(gdot - expected_coeff*denom))

      IF (err > tol .OR. projected_weighted <= 0.0_dp .OR. &
          recoverable_fraction <= 0.0_dp .OR. recoverable_fraction >= 1.0_dp) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "failed with error/gap/projected/recover", err, gap_weighted, projected_weighted, &
            recoverable_fraction
      ELSE
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "OK with error/gap/projected/recover", err, gap_weighted, projected_weighted, &
            recoverable_fraction
      END IF

   END SUBROUTINE check_projector_density_active_filter_descent_recovery

! **************************************************************************************************
!> \brief Check density-descent history calibration does not overstep the quadratic optimum.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_density_descent_history_optimum_cap(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      REAL(KIND=dp), PARAMETER :: gdotp = 0.30_dp, normp = 0.40_dp, &
         requested_history_scale = 1.60_dp

      REAL(KIND=dp)                                      :: applied_history_scale, base_drop, &
                                                            base_scale, err, legacy_drop, &
                                                            scaled_drop, scaled_scale

      base_scale = gdotp/normp
      base_drop = base_scale*gdotp - 0.5_dp*base_scale*base_scale*normp
      applied_history_scale = MIN(1.0_dp, requested_history_scale)
      scaled_scale = applied_history_scale*base_scale
      scaled_drop = scaled_scale*gdotp - 0.5_dp*scaled_scale*scaled_scale*normp
      legacy_drop = requested_history_scale*base_drop

      err = ABS(scaled_drop - base_drop)

      IF (err > tol .OR. scaled_drop <= 0.0_dp .OR. legacy_drop <= scaled_drop .OR. &
          applied_history_scale > 1.0_dp) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "failed with error/base/scaled/legacy", err, base_drop, scaled_drop, legacy_drop
      ELSE
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "OK with error/base/scaled/legacy", err, base_drop, scaled_drop, legacy_drop
      END IF

   END SUBROUTINE check_projector_density_descent_history_optimum_cap

! **************************************************************************************************
!> \brief Recover fixed-N density-gap response that is hidden from the linear gradient dot product.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_density_metric_descent_recovery(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: ne = 5
      REAL(KIND=dp), PARAMETER                           :: expected_coeff = 0.27_dp

      REAL(KIND=dp) :: coeff_grad, coeff_metric, denom_grad, err, gdot_grad, gdot_metric, &
         metric_cross, metric_denom, null_dot, old_norm, recovered_norm
      REAL(KIND=dp), DIMENSION(ne)                       :: descent, grad, metric_null, occ_curv, &
                                                            projected_grad, projected_metric, &
                                                            response

      occ_curv(:) = [0.045_dp, 0.026_dp, 0.018_dp, 0.010_dp, 0.004_dp]
      descent(:) = [0.18_dp, -0.10_dp, 0.06_dp, -0.03_dp, 0.015_dp]
      grad(:) = [0.36_dp, -0.21_dp, 0.09_dp, -0.05_dp, 0.02_dp]
      metric_null(:) = [-0.08_dp, 0.17_dp, 0.04_dp, -0.11_dp, 0.05_dp]

      CALL dense_projected_metric(metric_null, descent, occ_curv, metric_cross)
      CALL dense_projected_metric(descent, descent, occ_curv, metric_denom)
      metric_null(:) = metric_null(:) - metric_cross/MAX(metric_denom, EPSILON(metric_denom))*descent(:)
      null_dot = DOT_PRODUCT(grad, metric_null)
      denom_grad = DOT_PRODUCT(grad, descent)
      response(:) = expected_coeff*descent(:) - expected_coeff*denom_grad/null_dot*metric_null(:)

      CALL dense_project_descent_response(response, descent, grad, projected_grad, coeff_grad, gdot_grad, &
                                          denom_grad)
      CALL dense_project_metric_descent_response(response, descent, grad, occ_curv, projected_metric, &
                                                 coeff_metric, gdot_metric, denom_grad)

      CALL dense_projected_metric(response, descent, occ_curv, metric_cross)
      CALL dense_projected_metric(descent, descent, occ_curv, metric_denom)
      old_norm = MAXVAL(ABS(projected_grad(:)))
      recovered_norm = MAXVAL(ABS(projected_metric(:)))
      err = MAXVAL(ABS(projected_metric(:) - expected_coeff*descent(:)))
      err = MAX(err, ABS(coeff_metric - expected_coeff))
      err = MAX(err, ABS(gdot_metric))
      err = MAX(err, ABS(metric_cross/metric_denom - expected_coeff))

      IF (err > tol .OR. old_norm > tol .OR. recovered_norm <= 0.0_dp) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,5(1X,ES12.5))') TRIM(label), &
            "failed with error/old/recovered/gdot/mcoeff", err, old_norm, recovered_norm, &
            gdot_metric, coeff_metric
      ELSE
         WRITE (*, '(A,1X,A,5(1X,ES12.5))') TRIM(label), &
            "OK with error/old/recovered/gdot/mcoeff", err, old_norm, recovered_norm, &
            gdot_metric, coeff_metric
      END IF

   END SUBROUTINE check_projector_density_metric_descent_recovery

! **************************************************************************************************
!> \brief Recover density-gap response hidden in the fixed-N/tether Mermin metric.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_density_mixed_metric_descent_recovery(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: nchannel = 2, ne = 6
      REAL(KIND=dp), PARAMETER                           :: expected_coeff = 0.29_dp

      INTEGER, DIMENSION(nchannel)                       :: channel_nmo, e_offset
      REAL(KIND=dp) :: coeff_grad, coeff_metric, denom_grad, err, gdot_grad, gdot_metric, &
         metric_cross, metric_denom, null_dot, old_norm, recovered_norm
      REAL(KIND=dp), DIMENSION(ne)                       :: descent, grad, metric_null, occ_curv, &
                                                            projected_grad, projected_metric, &
                                                            response, tether_curv

      channel_nmo(:) = [3, 3]
      e_offset(:) = [1, 4]
      occ_curv(:) = [0.045_dp, 0.026_dp, 0.018_dp, 0.010_dp, 0.004_dp, 0.003_dp]
      tether_curv(:) = [0.008_dp, 0.012_dp, 0.010_dp, 0.006_dp, 0.009_dp, 0.005_dp]
      descent(:) = [0.18_dp, -0.10_dp, 0.06_dp, -0.03_dp, 0.015_dp, -0.007_dp]
      grad(:) = [0.36_dp, -0.21_dp, 0.09_dp, -0.05_dp, 0.02_dp, -0.01_dp]
      metric_null(:) = [-0.08_dp, 0.17_dp, 0.04_dp, -0.11_dp, 0.05_dp, -0.03_dp]

      CALL dense_projected_mixed_metric(metric_null, descent, occ_curv, tether_curv, channel_nmo, &
                                        e_offset, metric_cross)
      CALL dense_projected_mixed_metric(descent, descent, occ_curv, tether_curv, channel_nmo, &
                                        e_offset, metric_denom)
      metric_null(:) = metric_null(:) - metric_cross/MAX(metric_denom, EPSILON(metric_denom))*descent(:)
      null_dot = DOT_PRODUCT(grad, metric_null)
      denom_grad = DOT_PRODUCT(grad, descent)
      response(:) = expected_coeff*descent(:) - expected_coeff*denom_grad/null_dot*metric_null(:)

      CALL dense_project_descent_response(response, descent, grad, projected_grad, coeff_grad, &
                                          gdot_grad, denom_grad)
      CALL dense_project_mixed_metric_descent_response(response, descent, grad, occ_curv, tether_curv, &
                                                       channel_nmo, e_offset, projected_metric, &
                                                       coeff_metric, gdot_metric, denom_grad)

      CALL dense_projected_mixed_metric(response, descent, occ_curv, tether_curv, channel_nmo, &
                                        e_offset, metric_cross)
      CALL dense_projected_mixed_metric(descent, descent, occ_curv, tether_curv, channel_nmo, &
                                        e_offset, metric_denom)
      old_norm = MAXVAL(ABS(projected_grad(:)))
      recovered_norm = MAXVAL(ABS(projected_metric(:)))
      err = MAXVAL(ABS(projected_metric(:) - expected_coeff*descent(:)))
      err = MAX(err, ABS(coeff_metric - expected_coeff))
      err = MAX(err, ABS(gdot_metric))
      err = MAX(err, ABS(metric_cross/metric_denom - expected_coeff))

      IF (err > tol .OR. old_norm > tol .OR. recovered_norm <= 0.0_dp) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,5(1X,ES12.5))') TRIM(label), &
            "failed with error/old/recovered/gdot/mcoeff", err, old_norm, recovered_norm, &
            gdot_metric, coeff_metric
      ELSE
         WRITE (*, '(A,1X,A,5(1X,ES12.5))') TRIM(label), &
            "OK with error/old/recovered/gdot/mcoeff", err, old_norm, recovered_norm, &
            gdot_metric, coeff_metric
      END IF

   END SUBROUTINE check_projector_density_mixed_metric_descent_recovery

! **************************************************************************************************
!> \brief Recover density-gap response when the active response basis has collapsed.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_density_energy_descent_fallback(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: nchannel = 2, ne = 6
      REAL(KIND=dp), PARAMETER                           :: expected_coeff = 0.34_dp

      INTEGER, DIMENSION(nchannel)                       :: channel_nmo, e_offset
      REAL(KIND=dp) :: coeff_active, coeff_energy, denom_active, denom_energy, err, gdot_active, &
         gdot_energy, primary_norm, recovered_norm
      REAL(KIND=dp), DIMENSION(ne)                       :: active_descent, density_gap, &
                                                            energy_descent, grad, occ_curv, &
                                                            projected_active, projected_energy, &
                                                            tether_curv

      channel_nmo(:) = [3, 3]
      e_offset(:) = [1, 4]
      occ_curv(:) = [0.043_dp, 0.027_dp, 0.016_dp, 0.009_dp, 0.004_dp, 0.002_dp]
      tether_curv(:) = [0.007_dp, 0.011_dp, 0.009_dp, 0.005_dp, 0.008_dp, 0.004_dp]
      grad(:) = [0.41_dp, -0.20_dp, 0.12_dp, -0.06_dp, 0.025_dp, -0.012_dp]
      energy_descent(:) = [0.19_dp, -0.11_dp, 0.055_dp, -0.030_dp, 0.014_dp, -0.006_dp]

      active_descent(:) = 0.0_dp
      density_gap(:) = expected_coeff*energy_descent(:)

      CALL dense_project_mixed_metric_descent_response(density_gap, active_descent, grad, occ_curv, &
                                                       tether_curv, channel_nmo, e_offset, &
                                                       projected_active, coeff_active, gdot_active, &
                                                       denom_active)
      CALL dense_project_mixed_metric_descent_response(density_gap, energy_descent, grad, occ_curv, &
                                                       tether_curv, channel_nmo, e_offset, &
                                                       projected_energy, coeff_energy, gdot_energy, &
                                                       denom_energy)

      primary_norm = MAXVAL(ABS(projected_active(:)))
      recovered_norm = MAXVAL(ABS(projected_energy(:)))
      err = MAXVAL(ABS(projected_energy(:) - density_gap(:)))
      err = MAX(err, ABS(coeff_energy - expected_coeff))
      err = MAX(err, ABS(coeff_active), ABS(denom_active))

      IF (err > tol .OR. primary_norm > tol .OR. recovered_norm <= 0.0_dp .OR. &
          gdot_energy <= 0.0_dp .OR. denom_energy <= 0.0_dp) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,5(1X,ES12.5))') TRIM(label), &
            "failed with error/primary/recovered/gdot/denom", err, primary_norm, recovered_norm, &
            gdot_energy, denom_energy
      ELSE
         WRITE (*, '(A,1X,A,5(1X,ES12.5))') TRIM(label), &
            "OK with error/primary/recovered/gdot/denom", err, primary_norm, recovered_norm, &
            gdot_energy, denom_energy
      END IF

   END SUBROUTINE check_projector_density_energy_descent_fallback

! **************************************************************************************************
!> \brief Hidden density-gap response must enter the coupled mixed Mermin Schur curvature.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_density_mixed_metric_coupled_schur(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: nchannel = 2, ne = 6, nrot = 2

      INTEGER                                            :: info, r, s
      INTEGER, DIMENSION(nchannel)                       :: channel_nmo, e_offset
      INTEGER, DIMENSION(nrot)                           :: ipiv
      REAL(KIND=dp) :: blind_coeff, blind_denom, blind_gdot, curvature_gap, drop_active, &
         drop_full, drop_gap, err, metric_coeff, metric_cross, metric_denom, metric_term, rhs_gap, &
         solution_gap
      REAL(KIND=dp), DIMENSION(ne)                       :: grad_e, occ_curv, projected, tether_curv
      REAL(KIND=dp), DIMENSION(ne, nrot)                 :: density_gap, response_active, &
                                                            response_full, seed_gap
      REAL(KIND=dp), DIMENSION(nrot)                     :: grad_rot, rhs_active, rhs_full, &
                                                            sol_active, sol_full
      REAL(KIND=dp), DIMENSION(nrot, nrot)               :: rot_base, schur_active, schur_fact, &
                                                            schur_full

      channel_nmo(:) = [3, 3]
      e_offset(:) = [1, 4]
      occ_curv(:) = [0.045_dp, 0.026_dp, 0.018_dp, 0.010_dp, 0.004_dp, 0.003_dp]
      tether_curv(:) = [0.008_dp, 0.012_dp, 0.010_dp, 0.006_dp, 0.009_dp, 0.005_dp]
      grad_e(:) = [0.33_dp, -0.23_dp, 0.13_dp, -0.07_dp, 0.03_dp, -0.015_dp]
      grad_rot(:) = [0.12_dp, -0.09_dp]
      response_active(:, 1) = [0.16_dp, -0.09_dp, 0.05_dp, -0.03_dp, 0.015_dp, -0.008_dp]
      response_active(:, 2) = [-0.07_dp, 0.12_dp, -0.04_dp, 0.08_dp, -0.02_dp, 0.01_dp]
      seed_gap(:, 1) = 3.0_dp*[0.05_dp, 0.14_dp, -0.09_dp, 0.03_dp, -0.02_dp, 0.012_dp]
      seed_gap(:, 2) = 3.0_dp*[-0.02_dp, 0.08_dp, 0.11_dp, -0.06_dp, 0.04_dp, -0.03_dp]
      rot_base(:, :) = RESHAPE([0.72_dp, 0.05_dp, &
                                0.05_dp, 0.68_dp], [nrot, nrot])

      density_gap(:, :) = seed_gap(:, :)
      DO s = 1, nrot
         density_gap(:, s) = density_gap(:, s) - &
                             DOT_PRODUCT(grad_e(:), density_gap(:, s))/DOT_PRODUCT(grad_e(:), grad_e(:))*grad_e(:)
      END DO
      response_full(:, :) = response_active(:, :) + density_gap(:, :)

      schur_active(:, :) = rot_base(:, :)
      schur_full(:, :) = rot_base(:, :)
      DO s = 1, nrot
         DO r = 1, nrot
            CALL dense_projected_mixed_metric(response_active(:, r), response_active(:, s), occ_curv, &
                                              tether_curv, channel_nmo, e_offset, metric_term)
            schur_active(r, s) = schur_active(r, s) - metric_term
            CALL dense_projected_mixed_metric(response_full(:, r), response_full(:, s), occ_curv, &
                                              tether_curv, channel_nmo, e_offset, metric_term)
            schur_full(r, s) = schur_full(r, s) - metric_term
         END DO
      END DO

      rhs_active(:) = grad_rot(:) + MATMUL(TRANSPOSE(response_active), grad_e)
      rhs_full(:) = grad_rot(:) + MATMUL(TRANSPOSE(response_full), grad_e)
      CALL dense_project_descent_response(density_gap(:, 1), response_active(:, 1), grad_e, projected, &
                                          blind_coeff, blind_gdot, blind_denom)
      CALL dense_projected_mixed_metric(density_gap(:, 1), response_active(:, 1), occ_curv, tether_curv, &
                                        channel_nmo, e_offset, metric_cross)
      CALL dense_projected_mixed_metric(response_active(:, 1), response_active(:, 1), occ_curv, tether_curv, &
                                        channel_nmo, e_offset, metric_denom)
      metric_coeff = metric_cross/MAX(metric_denom, EPSILON(metric_denom))

      sol_active(:) = rhs_active(:)
      schur_fact(:, :) = schur_active(:, :)
      CALL dgesv(nrot, 1, schur_fact, nrot, ipiv, sol_active, nrot, info)
      IF (info /= 0) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,I6)') TRIM(label), "failed in active dgesv", info
         RETURN
      END IF

      sol_full(:) = rhs_full(:)
      schur_fact(:, :) = schur_full(:, :)
      CALL dgesv(nrot, 1, schur_fact, nrot, ipiv, sol_full, nrot, info)
      IF (info /= 0) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,I6)') TRIM(label), "failed in full dgesv", info
         RETURN
      END IF

      rhs_gap = MAXVAL(ABS(rhs_full(:) - rhs_active(:)))
      curvature_gap = MAXVAL(ABS(schur_full(:, :) - schur_active(:, :)))
      solution_gap = MAXVAL(ABS(sol_full(:) - sol_active(:)))
      drop_active = 0.5_dp*DOT_PRODUCT(rhs_active(:), sol_active(:))
      drop_full = 0.5_dp*DOT_PRODUCT(rhs_full(:), sol_full(:))
      drop_gap = ABS(drop_full - drop_active)
      err = rhs_gap
      err = MAX(err, MAXVAL(ABS(projected(:))))
      err = MAX(err, ABS(blind_gdot), ABS(blind_coeff))

      IF (err > tol .OR. curvature_gap <= 1.0E-4_dp .OR. solution_gap <= 1.0E-4_dp .OR. &
          drop_gap <= 1.0E-5_dp .OR. metric_coeff <= 1.0E-2_dp) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,5(1X,ES12.5))') TRIM(label), &
            "failed with error/curv/step/drop/mcoeff", err, curvature_gap, solution_gap, drop_gap, metric_coeff
      ELSE
         WRITE (*, '(A,1X,A,5(1X,ES12.5))') TRIM(label), &
            "OK with error/curv/step/drop/mcoeff", err, curvature_gap, solution_gap, drop_gap, metric_coeff
      END IF

   END SUBROUTINE check_projector_density_mixed_metric_coupled_schur

! **************************************************************************************************
!> \brief Projector-response candidate scaling must ignore fixed-N/tether gauge energy shifts.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_response_mixed_metric_scale(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: nchannel = 2, ne = 6

      INTEGER, DIMENSION(nchannel)                       :: channel_nmo, e_offset
      REAL(KIND=dp) :: drop_euclidean, drop_mixed, err, gauge_dot, gdotp, norm_energy_euclidean, &
         norm_energy_mixed, norm_euclidean, norm_mixed, rotation_gdot, rotation_norm, &
         scale_euclidean, scale_mixed
      REAL(KIND=dp), DIMENSION(ne)                       :: gauge, grad_e, occ_curv, response, &
                                                            response_physical, tether_curv

      channel_nmo(:) = [3, 3]
      e_offset(:) = [1, 4]
      occ_curv(:) = [0.045_dp, 0.026_dp, 0.018_dp, 0.010_dp, 0.004_dp, 0.003_dp]
      tether_curv(:) = [0.008_dp, 0.012_dp, 0.010_dp, 0.006_dp, 0.009_dp, 0.005_dp]
      grad_e(:) = [0.24_dp, -0.14_dp, -0.10_dp, 0.11_dp, -0.04_dp, -0.07_dp]
      response_physical(:) = [0.10_dp, -0.06_dp, -0.04_dp, 0.03_dp, -0.01_dp, -0.02_dp]
      gauge(:) = 0.80_dp
      response(:) = response_physical(:) + gauge(:)
      rotation_gdot = 0.04_dp
      rotation_norm = 0.20_dp

      gauge_dot = DOT_PRODUCT(grad_e, gauge)
      gdotp = DOT_PRODUCT(grad_e, response) + rotation_gdot
      norm_energy_euclidean = DOT_PRODUCT(response, response)
      CALL dense_projected_mixed_metric(response, response, occ_curv, tether_curv, channel_nmo, &
                                        e_offset, norm_energy_mixed)
      norm_energy_mixed = MAX(0.0_dp, norm_energy_mixed)
      norm_euclidean = norm_energy_euclidean + rotation_norm
      norm_mixed = norm_energy_mixed + rotation_norm
      scale_euclidean = gdotp/MAX(norm_euclidean, EPSILON(norm_euclidean))
      scale_mixed = gdotp/MAX(norm_mixed, EPSILON(norm_mixed))
      drop_euclidean = 0.5_dp*scale_euclidean*gdotp
      drop_mixed = 0.5_dp*scale_mixed*gdotp

      err = ABS(gauge_dot)
      IF (err > tol .OR. gdotp <= 0.0_dp .OR. norm_energy_mixed >= 0.05_dp*norm_energy_euclidean .OR. &
          scale_mixed <= 8.0_dp*scale_euclidean .OR. drop_mixed <= 8.0_dp*drop_euclidean) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,5(1X,ES12.5))') TRIM(label), &
            "failed with error/enorm/mnorm/escale/mscale", err, norm_energy_euclidean, norm_energy_mixed, &
            scale_euclidean, scale_mixed
      ELSE
         WRITE (*, '(A,1X,A,5(1X,ES12.5))') TRIM(label), &
            "OK with error/enorm/mnorm/escale/mscale", err, norm_energy_euclidean, norm_energy_mixed, &
            scale_euclidean, scale_mixed
      END IF

   END SUBROUTINE check_projector_response_mixed_metric_scale

! **************************************************************************************************
!> \brief Projector-response modal coefficients must use the projected Mermin linear term.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_response_projected_linear_term(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: nchannel = 2, ne = 6

      INTEGER                                            :: first_energy, i, j, local_channel, nmo
      INTEGER, DIMENSION(nchannel)                       :: channel_nmo, e_offset
      REAL(KIND=dp) :: base_neg, base_pos, coeff_det, coeff_neg, coeff_pos, err, g_neg, g_pos, &
         h_neg_neg, h_pos_neg, h_pos_pos, occ_sum, raw_coeff_neg, raw_coeff_pos, raw_g_neg, &
         raw_g_pos, raw_gap, target_neg, target_pos, weight_sum
      REAL(KIND=dp), DIMENSION(ne)                       :: corr, curv, energy_base, grad_e, &
                                                            occ_curv, response_neg, response_pos, &
                                                            tether_curv
      REAL(KIND=dp), DIMENSION(ne, ne)                   :: hproj

      channel_nmo(:) = [3, 3]
      e_offset(:) = [1, 4]
      occ_curv(:) = [0.045_dp, 0.026_dp, 0.018_dp, 0.010_dp, 0.004_dp, 0.003_dp]
      tether_curv(:) = [0.008_dp, 0.012_dp, 0.010_dp, 0.006_dp, 0.009_dp, 0.005_dp]
      curv(:) = occ_curv(:) + tether_curv(:)
      response_pos(:) = [0.12_dp, -0.07_dp, -0.05_dp, 0.04_dp, -0.015_dp, -0.025_dp]
      response_neg(:) = [-0.06_dp, 0.03_dp, 0.02_dp, 0.08_dp, -0.04_dp, -0.03_dp]
      energy_base(:) = [0.78_dp, 0.67_dp, 0.65_dp, -0.43_dp, -0.46_dp, -0.46_dp]
      target_pos = 0.65_dp
      target_neg = 0.25_dp
      corr(:) = target_pos*response_pos(:) + target_neg*response_neg(:)

      hproj(:, :) = 0.0_dp
      DO i = 1, ne
         hproj(i, i) = curv(i)
      END DO
      occ_sum = SUM(occ_curv(:))
      DO j = 1, ne
         DO i = 1, ne
            hproj(i, j) = hproj(i, j) - occ_curv(i)*occ_curv(j)/occ_sum
         END DO
      END DO
      DO local_channel = 1, nchannel
         nmo = channel_nmo(local_channel)
         first_energy = e_offset(local_channel)
         weight_sum = SUM(tether_curv(first_energy:first_energy + nmo - 1))
         DO j = first_energy, first_energy + nmo - 1
            DO i = first_energy, first_energy + nmo - 1
               hproj(i, j) = hproj(i, j) - tether_curv(i)*tether_curv(j)/weight_sum
            END DO
         END DO
      END DO
      grad_e(:) = MATMUL(hproj, energy_base(:) + corr(:))

      CALL dense_projected_mixed_metric(response_pos, response_pos, occ_curv, tether_curv, channel_nmo, &
                                        e_offset, h_pos_pos)
      CALL dense_projected_mixed_metric(response_pos, response_neg, occ_curv, tether_curv, channel_nmo, &
                                        e_offset, h_pos_neg)
      CALL dense_projected_mixed_metric(response_neg, response_neg, occ_curv, tether_curv, channel_nmo, &
                                        e_offset, h_neg_neg)
      CALL dense_projected_mixed_metric(energy_base, response_pos, occ_curv, tether_curv, channel_nmo, &
                                        e_offset, base_pos)
      CALL dense_projected_mixed_metric(energy_base, response_neg, occ_curv, tether_curv, channel_nmo, &
                                        e_offset, base_neg)

      g_pos = DOT_PRODUCT(grad_e(:), response_pos(:)) - base_pos
      g_neg = DOT_PRODUCT(grad_e(:), response_neg(:)) - base_neg
      coeff_det = h_pos_pos*h_neg_neg - h_pos_neg*h_pos_neg
      coeff_pos = (g_pos*h_neg_neg - g_neg*h_pos_neg)/coeff_det
      coeff_neg = (h_pos_pos*g_neg - h_pos_neg*g_pos)/coeff_det

      raw_g_pos = DOT_PRODUCT(grad_e(:) - curv(:)*energy_base(:), response_pos(:))
      raw_g_neg = DOT_PRODUCT(grad_e(:) - curv(:)*energy_base(:), response_neg(:))
      raw_coeff_pos = (raw_g_pos*h_neg_neg - raw_g_neg*h_pos_neg)/coeff_det
      raw_coeff_neg = (h_pos_pos*raw_g_neg - h_pos_neg*raw_g_pos)/coeff_det

      err = MAX(ABS(coeff_pos - target_pos), ABS(coeff_neg - target_neg))
      err = MAX(err, ABS(g_pos - (target_pos*h_pos_pos + target_neg*h_pos_neg)))
      err = MAX(err, ABS(g_neg - (target_pos*h_pos_neg + target_neg*h_neg_neg)))
      raw_gap = MAX(ABS(raw_coeff_pos - target_pos), ABS(raw_coeff_neg - target_neg))

      IF (err > tol .OR. raw_gap < 0.1_dp .OR. coeff_det <= EPSILON(coeff_det)) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,7(1X,ES12.5))') TRIM(label), &
            "failed with error/rawgap/cpos/cneg/rpos/rneg/det", err, raw_gap, coeff_pos, coeff_neg, &
            raw_coeff_pos, raw_coeff_neg, coeff_det
      ELSE
         WRITE (*, '(A,1X,A,7(1X,ES12.5))') TRIM(label), &
            "OK with error/rawgap/cpos/cneg/rpos/rneg/det", err, raw_gap, coeff_pos, coeff_neg, &
            raw_coeff_pos, raw_coeff_neg, coeff_det
      END IF

   END SUBROUTINE check_projector_response_projected_linear_term

! **************************************************************************************************
!> \brief Hidden density-gap response must enter the coupled fixed-N Schur curvature.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_density_metric_coupled_schur(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: ne = 5, nrot = 2

      INTEGER                                            :: info, r, s
      INTEGER, DIMENSION(nrot)                           :: ipiv
      REAL(KIND=dp) :: blind_coeff, blind_denom, blind_gdot, curvature_gap, drop_active, &
         drop_full, drop_gap, err, metric_term, rhs_gap, solution_gap
      REAL(KIND=dp), DIMENSION(ne)                       :: grad_e, occ_curv, projected
      REAL(KIND=dp), DIMENSION(ne, nrot)                 :: density_gap, response_active, &
                                                            response_full, seed_gap
      REAL(KIND=dp), DIMENSION(nrot)                     :: grad_rot, rhs_active, rhs_full, &
                                                            sol_active, sol_full
      REAL(KIND=dp), DIMENSION(nrot, nrot)               :: rot_base, schur_active, schur_fact, &
                                                            schur_full

      occ_curv(:) = [0.045_dp, 0.026_dp, 0.018_dp, 0.010_dp, 0.004_dp]
      grad_e(:) = [0.33_dp, -0.23_dp, 0.13_dp, -0.07_dp, 0.03_dp]
      grad_rot(:) = [0.12_dp, -0.09_dp]
      response_active(:, 1) = [0.16_dp, -0.09_dp, 0.05_dp, -0.03_dp, 0.015_dp]
      response_active(:, 2) = [-0.07_dp, 0.12_dp, -0.04_dp, 0.08_dp, -0.02_dp]
      seed_gap(:, 1) = 3.0_dp*[0.05_dp, 0.14_dp, -0.09_dp, 0.03_dp, -0.02_dp]
      seed_gap(:, 2) = 3.0_dp*[-0.02_dp, 0.08_dp, 0.11_dp, -0.06_dp, 0.04_dp]
      rot_base(:, :) = RESHAPE([0.72_dp, 0.05_dp, &
                                0.05_dp, 0.68_dp], [nrot, nrot])

      density_gap(:, :) = seed_gap(:, :)
      DO s = 1, nrot
         density_gap(:, s) = density_gap(:, s) - &
                             DOT_PRODUCT(grad_e(:), density_gap(:, s))/DOT_PRODUCT(grad_e(:), grad_e(:))*grad_e(:)
      END DO
      response_full(:, :) = response_active(:, :) + density_gap(:, :)

      schur_active(:, :) = rot_base(:, :)
      schur_full(:, :) = rot_base(:, :)
      DO s = 1, nrot
         DO r = 1, nrot
            CALL dense_projected_metric(response_active(:, r), response_active(:, s), occ_curv, metric_term)
            schur_active(r, s) = schur_active(r, s) - metric_term
            CALL dense_projected_metric(response_full(:, r), response_full(:, s), occ_curv, metric_term)
            schur_full(r, s) = schur_full(r, s) - metric_term
         END DO
      END DO

      rhs_active(:) = grad_rot(:) + MATMUL(TRANSPOSE(response_active), grad_e)
      rhs_full(:) = grad_rot(:) + MATMUL(TRANSPOSE(response_full), grad_e)
      CALL dense_project_descent_response(density_gap(:, 1), response_active(:, 1), grad_e, projected, &
                                          blind_coeff, blind_gdot, blind_denom)

      sol_active(:) = rhs_active(:)
      schur_fact(:, :) = schur_active(:, :)
      CALL dgesv(nrot, 1, schur_fact, nrot, ipiv, sol_active, nrot, info)
      IF (info /= 0) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,I6)') TRIM(label), "failed in active dgesv", info
         RETURN
      END IF

      sol_full(:) = rhs_full(:)
      schur_fact(:, :) = schur_full(:, :)
      CALL dgesv(nrot, 1, schur_fact, nrot, ipiv, sol_full, nrot, info)
      IF (info /= 0) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,I6)') TRIM(label), "failed in full dgesv", info
         RETURN
      END IF

      rhs_gap = MAXVAL(ABS(rhs_full(:) - rhs_active(:)))
      curvature_gap = MAXVAL(ABS(schur_full(:, :) - schur_active(:, :)))
      solution_gap = MAXVAL(ABS(sol_full(:) - sol_active(:)))
      drop_active = 0.5_dp*DOT_PRODUCT(rhs_active(:), sol_active(:))
      drop_full = 0.5_dp*DOT_PRODUCT(rhs_full(:), sol_full(:))
      drop_gap = ABS(drop_full - drop_active)
      err = rhs_gap
      err = MAX(err, MAXVAL(ABS(projected(:))))
      err = MAX(err, ABS(blind_gdot), ABS(blind_coeff))

      IF (err > tol .OR. curvature_gap <= 1.0E-4_dp .OR. solution_gap <= 1.0E-4_dp .OR. &
          drop_gap <= 1.0E-5_dp) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "failed with error/curv/step/drop", err, curvature_gap, solution_gap, drop_gap
      ELSE
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "OK with error/curv/step/drop", err, curvature_gap, solution_gap, drop_gap
      END IF

   END SUBROUTINE check_projector_density_metric_coupled_schur

! **************************************************************************************************
!> \brief Remove near-degenerate gauge noise from the off-diagonal Projector-Hessian channel.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_complex_projector_offdiag_gauge_plateau_hessian(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 4, nrot = n*(n - 1)
      REAL(KIND=dp), PARAMETER                           :: active_tol = 1.0E-4_dp, &
                                                            gauge_gap = 3.0E-5_dp

      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: chc, chc_eig, dchc, dchc_eig, &
                                                            docc_plateau, docc_raw, q
      INTEGER                                            :: i, ivar, j
      INTEGER, DIMENSION(nrot)                           :: var_i, var_j, var_kind
      REAL(KIND=dp) :: err, physical_active, physical_capped_norm, physical_err, &
         physical_raw_norm, plateau_norm, raw_gap_norm, raw_norm
      REAL(KIND=dp), DIMENSION(n)                        :: eps_diag, occ_curv, occ_scale
      REAL(KIND=dp), DIMENSION(n, n)                     :: grad_im, grad_re, x_im, x_re
      REAL(KIND=dp), DIMENSION(nrot)                     :: grad_physical_vec, grad_plateau_vec, &
                                                            grad_raw_vec

      chc_eig(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      chc_eig(1, 1) = CMPLX(-0.010000_dp, 0.0_dp, KIND=dp)
      chc_eig(2, 2) = CMPLX(-0.009982_dp, 0.0_dp, KIND=dp)
      chc_eig(3, 3) = CMPLX(0.24_dp, 0.0_dp, KIND=dp)
      chc_eig(4, 4) = CMPLX(0.52_dp, 0.0_dp, KIND=dp)

      x_re(:, :) = 0.0_dp
      x_re(1, 2) = 0.12_dp
      x_re(1, 3) = -0.06_dp
      x_re(2, 4) = 0.05_dp
      x_re(3, 4) = -0.04_dp
      DO i = 2, n
         x_re(i, 1:i - 1) = -x_re(1:i - 1, i)
      END DO
      x_im(:, :) = 0.0_dp
      x_im(1, 1) = 0.02_dp
      x_im(2, 2) = -0.03_dp
      x_im(3, 3) = 0.01_dp
      x_im(4, 4) = -0.02_dp
      x_im(1, 2) = -0.04_dp
      x_im(1, 4) = 0.03_dp
      x_im(2, 3) = 0.02_dp
      DO i = 2, n
         x_im(i, 1:i - 1) = x_im(1:i - 1, i)
      END DO
      CALL dense_rotation_u(x_re, x_im, q)
      chc(:, :) = MATMUL(q, MATMUL(chc_eig, CONJG(TRANSPOSE(q))))

      ivar = 0
      DO j = 2, n
         DO i = 1, j - 1
            ivar = ivar + 1
            var_i(ivar) = i
            var_j(ivar) = j
            var_kind(ivar) = 1
            ivar = ivar + 1
            var_i(ivar) = i
            var_j(ivar) = j
            var_kind(ivar) = 2
         END DO
      END DO

      occ_scale(:) = [1.70_dp, 0.46_dp, 0.24_dp, 0.08_dp]
      occ_curv(:) = [0.038_dp, 0.035_dp, 0.018_dp, 0.006_dp]

      dchc_eig(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      dchc_eig(1, 2) = CMPLX(0.25_dp, -0.18_dp, KIND=dp)
      dchc_eig(2, 1) = CONJG(dchc_eig(1, 2))
      dchc(:, :) = MATMUL(q, MATMUL(dchc_eig, CONJG(TRANSPOSE(q))))
      CALL dense_spectral_occ_projector_derivative_chc(chc, dchc, occ_scale, docc_raw)
      CALL dense_spectral_occ_projector_derivative_chc(chc, dchc, occ_scale, docc_plateau, &
                                                       occ_curv, gauge_gap)
      DO i = 1, n
         docc_raw(i, i) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
         docc_plateau(i, i) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      END DO
      x_re(:, :) = 0.0_dp
      x_im(:, :) = 0.0_dp
      CALL dense_projector_rotation_gradient(x_re, x_im, chc, docc_raw, grad_re, grad_im, eps_diag)
      CALL extract_rotation_gradient_vector_2d(grad_re, grad_im, var_i, var_j, var_kind, grad_raw_vec)
      CALL dense_projector_rotation_gradient(x_re, x_im, chc, docc_plateau, grad_re, grad_im, eps_diag)
      CALL extract_rotation_gradient_vector_2d(grad_re, grad_im, var_i, var_j, var_kind, grad_plateau_vec)
      raw_norm = MAXVAL(ABS(grad_raw_vec(:)))
      plateau_norm = MAXVAL(ABS(grad_plateau_vec(:)))
      raw_gap_norm = MAXVAL(ABS(grad_raw_vec(:) - grad_plateau_vec(:)))

      dchc_eig(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      dchc_eig(3, 4) = CMPLX(-0.17_dp, 0.12_dp, KIND=dp)
      dchc_eig(4, 3) = CONJG(dchc_eig(3, 4))
      dchc(:, :) = MATMUL(q, MATMUL(dchc_eig, CONJG(TRANSPOSE(q))))
      CALL dense_spectral_occ_projector_derivative_chc(chc, dchc, occ_scale, docc_raw)
      CALL dense_spectral_occ_projector_derivative_chc(chc, dchc, occ_scale, docc_plateau, &
                                                       occ_curv, gauge_gap)
      DO i = 1, n
         docc_raw(i, i) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
         docc_plateau(i, i) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      END DO
      CALL dense_projector_rotation_gradient(x_re, x_im, chc, docc_raw, grad_re, grad_im, eps_diag)
      CALL extract_rotation_gradient_vector_2d(grad_re, grad_im, var_i, var_j, var_kind, grad_physical_vec)
      CALL dense_projector_rotation_gradient(x_re, x_im, chc, docc_plateau, grad_re, grad_im, eps_diag)
      CALL extract_rotation_gradient_vector_2d(grad_re, grad_im, var_i, var_j, var_kind, grad_plateau_vec)
      physical_raw_norm = MAXVAL(ABS(grad_physical_vec(:)))
      physical_capped_norm = MAXVAL(ABS(grad_plateau_vec(:)))
      physical_err = MAX(0.0_dp, physical_capped_norm - physical_raw_norm)
      physical_active = MAXVAL(ABS(grad_plateau_vec(:)))

      err = MAX(plateau_norm, physical_err)
      IF (err > tol .OR. raw_norm < active_tol .OR. raw_gap_norm < active_tol .OR. &
          physical_active < active_tol .OR. &
          physical_capped_norm >= 0.5_dp*physical_raw_norm) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,5(1X,ES12.5))') TRIM(label), &
            "failed with error/raw/gap/plateau/physical", err, raw_norm, raw_gap_norm, plateau_norm, physical_active
      ELSE
         WRITE (*, '(A,1X,A,5(1X,ES12.5))') TRIM(label), &
            "OK with error/raw/gap/plateau/physical", err, raw_norm, raw_gap_norm, plateau_norm, physical_active
      END IF

   END SUBROUTINE check_complex_projector_offdiag_gauge_plateau_hessian

! **************************************************************************************************
!> \brief Check the fixed-N low-rank curvature from diagonal projector-density response.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_complex_projector_density_curvature_addback(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: ne = 4, nrot = 3

      INTEGER                                            :: i, info, j, s
      INTEGER, DIMENSION(ne)                             :: ipiv
      REAL(KIND=dp)                                      :: err, gauge_shift, occ_sum, qdrop, &
                                                            ref_drop, shift
      REAL(KIND=dp), DIMENSION(ne)                       :: bq, occ_curv, tether_curv, xq
      REAL(KIND=dp), DIMENSION(ne, ne)                   :: h, h_ref
      REAL(KIND=dp), DIMENSION(ne, nrot)                 :: b, b_gauge, density_response, rhs
      REAL(KIND=dp), DIMENSION(nrot)                     :: qvec
      REAL(KIND=dp), DIMENSION(nrot, nrot)               :: curv, curv_gauge

      occ_curv(:) = [0.38_dp, 0.71_dp, 0.29_dp, 0.16_dp]
      tether_curv(:) = [0.05_dp, 0.02_dp, 0.03_dp, 0.01_dp]
      density_response(:, 1) = [0.09_dp, -0.04_dp, 0.16_dp, -0.08_dp]
      density_response(:, 2) = [-0.11_dp, 0.07_dp, 0.03_dp, 0.10_dp]
      density_response(:, 3) = [0.05_dp, 0.13_dp, -0.09_dp, -0.02_dp]
      qvec(:) = [0.43_dp, -0.28_dp, 0.19_dp]

      h(:, :) = 0.0_dp
      DO i = 1, ne
         h(i, i) = occ_curv(i) + tether_curv(i)
      END DO
      occ_sum = SUM(occ_curv)
      DO j = 1, ne
         DO i = 1, ne
            h(i, j) = h(i, j) - occ_curv(i)*occ_curv(j)/occ_sum
         END DO
      END DO
      h(:, :) = h(:, :) + MAXVAL(occ_curv(:) + tether_curv(:))/REAL(ne, KIND=dp)
      b(:, :) = 0.0_dp
      DO s = 1, nrot
         shift = DOT_PRODUCT(occ_curv, density_response(:, s))/occ_sum
         b(:, s) = occ_curv(:)*(density_response(:, s) - shift)
      END DO

      rhs(:, :) = b(:, :)
      h_ref(:, :) = h(:, :)
      CALL dgesv(ne, nrot, h, ne, ipiv, rhs, ne, info)
      IF (info /= 0) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,I6)') TRIM(label), "failed in dgesv", info
         RETURN
      END IF
      curv(:, :) = MATMUL(TRANSPOSE(b), rhs)

      gauge_shift = 0.37_dp
      b_gauge(:, :) = 0.0_dp
      DO s = 1, nrot
         shift = DOT_PRODUCT(occ_curv, density_response(:, s) + gauge_shift)/occ_sum
         b_gauge(:, s) = occ_curv(:)*(density_response(:, s) + gauge_shift - shift)
      END DO
      h(:, :) = h_ref(:, :)
      rhs(:, :) = b_gauge(:, :)
      CALL dgesv(ne, nrot, h, ne, ipiv, rhs, ne, info)
      IF (info /= 0) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,I6)') TRIM(label), "failed in gauge dgesv", info
         RETURN
      END IF
      curv_gauge(:, :) = MATMUL(TRANSPOSE(b_gauge), rhs)

      bq(:) = MATMUL(b, qvec)
      xq(:) = MATMUL(rhs, qvec)
      qdrop = 0.5_dp*DOT_PRODUCT(qvec, MATMUL(curv, qvec))
      ref_drop = 0.5_dp*DOT_PRODUCT(bq, xq)
      err = MAXVAL(ABS(curv(:, :) - TRANSPOSE(curv(:, :))))
      err = MAX(err, MAXVAL(ABS(curv(:, :) - curv_gauge(:, :))))
      err = MAX(err, ABS(qdrop - ref_drop))

      IF (err > tol .OR. qdrop <= 1.0E-8_dp) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,2(1X,ES12.5))') TRIM(label), "failed with error/drop", err, qdrop
      ELSE
         WRITE (*, '(A,1X,A,2(1X,ES12.5))') TRIM(label), "OK with error/drop", err, qdrop
      END IF

   END SUBROUTINE check_complex_projector_density_curvature_addback

! **************************************************************************************************
!> \brief Check the algebra of promoting diagonal density response into the Schur block.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_complex_projector_density_coupled_schur(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: ne = 4, nrot = 3, ntot = ne + nrot

      INTEGER                                            :: i, info, j, s
      INTEGER, DIMENSION(ne)                             :: ipiv
      INTEGER, DIMENSION(ntot)                           :: ipiv_aug
      REAL(KIND=dp)                                      :: aug_err, err, gauge_shift, norm_delta, &
                                                            occ_sum, shift
      REAL(KIND=dp), DIMENSION(ne)                       :: grad_e, occ_curv, tether_curv, z_aug, &
                                                            z_red
      REAL(KIND=dp), DIMENSION(ne, ne)                   :: h, h_ref
      REAL(KIND=dp), DIMENSION(ne, nrot)                 :: b_density, b_phys, b_total, &
                                                            density_response, phys_response, &
                                                            rhs_density, rhs_phys, rhs_total
      REAL(KIND=dp), DIMENSION(nrot)                     :: grad_rot, q_aug, q_red
      REAL(KIND=dp), DIMENSION(nrot, nrot)               :: a, combined_schur, curv_density, &
                                                            expected_combined, rot_base, &
                                                            schur_solve, split_schur
      REAL(KIND=dp), DIMENSION(ntot)                     :: rhs_aug
      REAL(KIND=dp), DIMENSION(ntot, ntot)               :: aug

      occ_curv(:) = [0.44_dp, 0.67_dp, 0.31_dp, 0.18_dp]
      tether_curv(:) = [0.03_dp, 0.04_dp, 0.02_dp, 0.01_dp]
      phys_response(:, 1) = [0.18_dp, -0.07_dp, 0.04_dp, -0.02_dp]
      phys_response(:, 2) = [-0.09_dp, 0.11_dp, -0.03_dp, 0.08_dp]
      phys_response(:, 3) = [0.05_dp, 0.02_dp, -0.12_dp, 0.06_dp]
      density_response(:, 1) = [0.07_dp, -0.03_dp, 0.10_dp, -0.06_dp]
      density_response(:, 2) = [-0.04_dp, 0.09_dp, 0.02_dp, 0.05_dp]
      density_response(:, 3) = [0.03_dp, 0.08_dp, -0.05_dp, -0.01_dp]
      grad_e(:) = [0.14_dp, -0.06_dp, 0.03_dp, -0.10_dp]
      grad_rot(:) = [0.18_dp, -0.13_dp, 0.09_dp]
      rot_base(:, :) = RESHAPE([0.63_dp, 0.04_dp, -0.02_dp, &
                                0.04_dp, 0.51_dp, 0.03_dp, &
                                -0.02_dp, 0.03_dp, 0.58_dp], [nrot, nrot])

      h(:, :) = 0.0_dp
      DO i = 1, ne
         h(i, i) = occ_curv(i) + tether_curv(i)
      END DO
      occ_sum = SUM(occ_curv)
      DO j = 1, ne
         DO i = 1, ne
            h(i, j) = h(i, j) - occ_curv(i)*occ_curv(j)/occ_sum
         END DO
      END DO
      h(:, :) = h(:, :) + MAXVAL(occ_curv(:) + tether_curv(:))/REAL(ne, KIND=dp)
      h_ref(:, :) = h(:, :)

      DO s = 1, nrot
         shift = DOT_PRODUCT(occ_curv, phys_response(:, s))/occ_sum
         b_phys(:, s) = occ_curv(:)*(phys_response(:, s) - shift)
         shift = DOT_PRODUCT(occ_curv, density_response(:, s))/occ_sum
         b_density(:, s) = occ_curv(:)*(density_response(:, s) - shift)
      END DO
      b_total(:, :) = b_phys(:, :) + b_density(:, :)

      rhs_phys(:, :) = b_phys(:, :)
      h(:, :) = h_ref(:, :)
      CALL dgesv(ne, nrot, h, ne, ipiv, rhs_phys, ne, info)
      IF (info /= 0) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,I6)') TRIM(label), "failed in phys dgesv", info
         RETURN
      END IF

      rhs_density(:, :) = b_density(:, :)
      h(:, :) = h_ref(:, :)
      CALL dgesv(ne, nrot, h, ne, ipiv, rhs_density, ne, info)
      IF (info /= 0) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,I6)') TRIM(label), "failed in density dgesv", info
         RETURN
      END IF

      rhs_total(:, :) = b_total(:, :)
      h(:, :) = h_ref(:, :)
      CALL dgesv(ne, nrot, h, ne, ipiv, rhs_total, ne, info)
      IF (info /= 0) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,I6)') TRIM(label), "failed in total dgesv", info
         RETURN
      END IF

      curv_density(:, :) = MATMUL(TRANSPOSE(b_density), rhs_density)
      split_schur(:, :) = rot_base(:, :) + curv_density(:, :) - MATMUL(TRANSPOSE(b_phys), rhs_phys)
      combined_schur(:, :) = rot_base(:, :) + curv_density(:, :) - MATMUL(TRANSPOSE(b_total), rhs_total)
      expected_combined(:, :) = split_schur(:, :) - MATMUL(TRANSPOSE(b_phys), rhs_density) - &
                                MATMUL(TRANSPOSE(b_density), rhs_phys) - curv_density(:, :)
      norm_delta = MAXVAL(ABS(combined_schur(:, :) - split_schur(:, :)))

      a(:, :) = rot_base(:, :) + curv_density(:, :)
      z_red(:) = grad_e(:)
      h(:, :) = h_ref(:, :)
      CALL dgesv(ne, 1, h, ne, ipiv, z_red, ne, info)
      IF (info /= 0) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,I6)') TRIM(label), "failed in gradient dgesv", info
         RETURN
      END IF
      q_red(:) = grad_rot(:) + MATMUL(TRANSPOSE(rhs_total), grad_e)
      schur_solve(:, :) = combined_schur(:, :)
      CALL dgesv(nrot, 1, schur_solve, nrot, ipiv(1:nrot), q_red, nrot, info)
      IF (info /= 0) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,I6)') TRIM(label), "failed in coupled Schur dgesv", info
         RETURN
      END IF
      z_red(:) = z_red(:) + MATMUL(rhs_total, q_red)

      aug(:, :) = 0.0_dp
      aug(1:ne, 1:ne) = h_ref(:, :)
      aug(1:ne, ne + 1:ntot) = -b_total(:, :)
      aug(ne + 1:ntot, 1:ne) = -TRANSPOSE(b_total(:, :))
      aug(ne + 1:ntot, ne + 1:ntot) = a(:, :)
      rhs_aug(1:ne) = grad_e(:)
      rhs_aug(ne + 1:ntot) = grad_rot(:)
      CALL dgesv(ntot, 1, aug, ntot, ipiv_aug, rhs_aug, ntot, info)
      IF (info /= 0) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,I6)') TRIM(label), "failed in coupled augmented dgesv", info
         RETURN
      END IF
      z_aug(:) = rhs_aug(1:ne)
      q_aug(:) = rhs_aug(ne + 1:ntot)
      aug_err = MAX(MAXVAL(ABS(z_red - z_aug)), MAXVAL(ABS(q_red - q_aug)))

      gauge_shift = 0.41_dp
      DO s = 1, nrot
         shift = DOT_PRODUCT(occ_curv, density_response(:, s) + gauge_shift)/occ_sum
         b_total(:, s) = b_phys(:, s) + occ_curv(:)*(density_response(:, s) + gauge_shift - shift)
      END DO
      rhs_total(:, :) = b_total(:, :)
      h(:, :) = h_ref(:, :)
      CALL dgesv(ne, nrot, h, ne, ipiv, rhs_total, ne, info)
      IF (info /= 0) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,I6)') TRIM(label), "failed in gauge dgesv", info
         RETURN
      END IF

      err = MAXVAL(ABS(combined_schur(:, :) - expected_combined(:, :)))
      err = MAX(err, aug_err)
      err = MAX(err, MAXVAL(ABS(combined_schur(:, :) - &
                                (rot_base(:, :) + curv_density(:, :) - &
                                 MATMUL(TRANSPOSE(b_total), rhs_total)))))

      IF (err > tol .OR. norm_delta <= 1.0E-6_dp) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,2(1X,ES12.5))') TRIM(label), "failed with error/delta", err, norm_delta
      ELSE
         WRITE (*, '(A,1X,A,2(1X,ES12.5))') TRIM(label), "OK with error/delta", err, norm_delta
      END IF

   END SUBROUTINE check_complex_projector_density_coupled_schur

! **************************************************************************************************
!> \brief Check that physical-direct shadows augment Schur, RHS and model response with density response.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_physical_direct_density_shadow(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: ne = 4, nrot = 3, ntot = ne + nrot

      INTEGER                                            :: i, info, j, s
      INTEGER, DIMENSION(ne)                             :: ipiv
      INTEGER, DIMENSION(nrot)                           :: ipiv_rot
      INTEGER, DIMENSION(ntot)                           :: ipiv_aug
      REAL(KIND=dp)                                      :: aug_err, err, legacy_gap, occ_sum, shift
      REAL(KIND=dp), DIMENSION(ne)                       :: grad_e, occ_curv, z_aug, z_full, z_legacy
      REAL(KIND=dp), DIMENSION(ne, ne)                   :: h, h_ref
      REAL(KIND=dp), DIMENSION(ne, nrot)                 :: b_density, b_phys, b_total, &
                                                            density_response, phys_response, &
                                                            rhs_density, rhs_phys, rhs_total
      REAL(KIND=dp), DIMENSION(nrot)                     :: grad_rot, q_aug, q_full, q_legacy
      REAL(KIND=dp), DIMENSION(nrot, nrot)               :: rot_delta, rot_metric, schur_full, &
                                                            schur_legacy
      REAL(KIND=dp), DIMENSION(ntot)                     :: rhs_aug
      REAL(KIND=dp), DIMENSION(ntot, ntot)               :: aug

      occ_curv(:) = [0.41_dp, 0.72_dp, 0.29_dp, 0.16_dp]
      phys_response(:, 1) = [0.16_dp, -0.08_dp, 0.05_dp, -0.01_dp]
      phys_response(:, 2) = [-0.07_dp, 0.13_dp, -0.04_dp, 0.09_dp]
      phys_response(:, 3) = [0.06_dp, 0.03_dp, -0.11_dp, 0.04_dp]
      density_response(:, 1) = [0.08_dp, -0.02_dp, 0.11_dp, -0.05_dp]
      density_response(:, 2) = [-0.03_dp, 0.10_dp, 0.03_dp, 0.04_dp]
      density_response(:, 3) = [0.04_dp, 0.07_dp, -0.06_dp, -0.02_dp]
      grad_e(:) = [0.12_dp, -0.05_dp, 0.04_dp, -0.09_dp]
      grad_rot(:) = [0.19_dp, -0.11_dp, 0.08_dp]
      rot_delta(:, :) = RESHAPE([0.12_dp, 0.03_dp, -0.01_dp, &
                                 0.03_dp, 0.08_dp, 0.02_dp, &
                                 -0.01_dp, 0.02_dp, 0.10_dp], [nrot, nrot])
      rot_metric(:, :) = RESHAPE([0.69_dp, 0.05_dp, -0.03_dp, &
                                  0.05_dp, 0.56_dp, 0.04_dp, &
                                  -0.03_dp, 0.04_dp, 0.61_dp], [nrot, nrot]) + rot_delta(:, :)

      h(:, :) = 0.0_dp
      DO i = 1, ne
         h(i, i) = occ_curv(i) + 0.02_dp*REAL(i, KIND=dp)
      END DO
      occ_sum = SUM(occ_curv)
      DO j = 1, ne
         DO i = 1, ne
            h(i, j) = h(i, j) - occ_curv(i)*occ_curv(j)/occ_sum
         END DO
      END DO
      h(:, :) = h(:, :) + MAXVAL(occ_curv(:))/REAL(ne, KIND=dp)
      h_ref(:, :) = h(:, :)

      DO s = 1, nrot
         shift = DOT_PRODUCT(occ_curv, phys_response(:, s))/occ_sum
         b_phys(:, s) = occ_curv(:)*(phys_response(:, s) - shift)
         shift = DOT_PRODUCT(occ_curv, density_response(:, s))/occ_sum
         b_density(:, s) = occ_curv(:)*(density_response(:, s) - shift)
      END DO
      b_total(:, :) = b_phys(:, :) + b_density(:, :)

      rhs_phys(:, :) = b_phys(:, :)
      h(:, :) = h_ref(:, :)
      CALL dgesv(ne, nrot, h, ne, ipiv, rhs_phys, ne, info)
      IF (info /= 0) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,I6)') TRIM(label), "failed in phys dgesv", info
         RETURN
      END IF

      rhs_density(:, :) = b_density(:, :)
      h(:, :) = h_ref(:, :)
      CALL dgesv(ne, nrot, h, ne, ipiv, rhs_density, ne, info)
      IF (info /= 0) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,I6)') TRIM(label), "failed in density dgesv", info
         RETURN
      END IF
      rhs_total(:, :) = rhs_phys(:, :) + rhs_density(:, :)

      schur_full(:, :) = rot_metric(:, :) - MATMUL(TRANSPOSE(b_total), rhs_total)
      schur_legacy(:, :) = rot_metric(:, :) - MATMUL(TRANSPOSE(b_phys), rhs_phys)

      q_full(:) = grad_rot(:) + MATMUL(TRANSPOSE(rhs_total), grad_e)
      q_legacy(:) = grad_rot(:) + MATMUL(TRANSPOSE(rhs_phys), grad_e)
      CALL dgesv(nrot, 1, schur_full, nrot, ipiv_rot, q_full, nrot, info)
      IF (info /= 0) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,I6)') TRIM(label), "failed in full shadow dgesv", info
         RETURN
      END IF
      CALL dgesv(nrot, 1, schur_legacy, nrot, ipiv_rot, q_legacy, nrot, info)
      IF (info /= 0) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,I6)') TRIM(label), "failed in legacy shadow dgesv", info
         RETURN
      END IF

      h(:, :) = h_ref(:, :)
      z_full(:) = grad_e(:)
      CALL dgesv(ne, 1, h, ne, ipiv, z_full, ne, info)
      IF (info /= 0) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,I6)') TRIM(label), "failed in z dgesv", info
         RETURN
      END IF
      z_full(:) = z_full(:) + MATMUL(rhs_total, q_full)
      z_legacy(:) = MATMUL(rhs_phys, q_legacy)

      aug(:, :) = 0.0_dp
      aug(1:ne, 1:ne) = h_ref(:, :)
      aug(1:ne, ne + 1:ntot) = -b_total(:, :)
      aug(ne + 1:ntot, 1:ne) = -TRANSPOSE(b_total(:, :))
      aug(ne + 1:ntot, ne + 1:ntot) = rot_metric(:, :)
      rhs_aug(1:ne) = grad_e(:)
      rhs_aug(ne + 1:ntot) = grad_rot(:)
      CALL dgesv(ntot, 1, aug, ntot, ipiv_aug, rhs_aug, ntot, info)
      IF (info /= 0) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,I6)') TRIM(label), "failed in augmented dgesv", info
         RETURN
      END IF
      z_aug(:) = rhs_aug(1:ne)
      q_aug(:) = rhs_aug(ne + 1:ntot)

      aug_err = MAX(MAXVAL(ABS(z_full(:) - z_aug(:))), MAXVAL(ABS(q_full(:) - q_aug(:))))
      legacy_gap = MAX(MAXVAL(ABS(q_full(:) - q_legacy(:))), &
                       MAXVAL(ABS(MATMUL(rhs_total, q_full) - z_legacy(:))))
      err = aug_err

      IF (err > tol .OR. legacy_gap <= 1.0E-5_dp) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,2(1X,ES12.5))') TRIM(label), "failed with error/legacy_gap", &
            err, legacy_gap
      ELSE
         WRITE (*, '(A,1X,A,2(1X,ES12.5))') TRIM(label), "OK with error/legacy_gap", err, legacy_gap
      END IF

   END SUBROUTINE check_projector_physical_direct_density_shadow

! **************************************************************************************************
!> \brief Check that density-response candidates keep the physical response in the accepted z step.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_density_candidate_total_response(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: ne = 4, nrot = 2, ntot = ne + nrot

      INTEGER                                            :: i, info, j, s
      INTEGER, DIMENSION(ne)                             :: ipiv
      INTEGER, DIMENSION(nrot)                           :: ipiv_rot
      INTEGER, DIMENSION(ntot)                           :: ipiv_aug
      REAL(KIND=dp)                                      :: err, filter_gap, gap, occ_sum, shift
      REAL(KIND=dp), DIMENSION(ne)                       :: grad_e, occ_curv, z_aug, z_base, &
                                                            z_density_only, z_filtered, z_total
      REAL(KIND=dp), DIMENSION(ne, ne)                   :: h, h_ref
      REAL(KIND=dp), DIMENSION(ne, nrot)                 :: b_density, b_phys, b_total, &
                                                            density_response, phys_response, &
                                                            rhs_density, rhs_phys, rhs_total
      REAL(KIND=dp), DIMENSION(nrot)                     :: filter_weight, grad_rot, q_aug, &
                                                            q_filtered, q_total
      REAL(KIND=dp), DIMENSION(nrot, nrot)               :: rot_metric, schur_total
      REAL(KIND=dp), DIMENSION(ntot)                     :: rhs_aug
      REAL(KIND=dp), DIMENSION(ntot, ntot)               :: aug

      occ_curv(:) = [0.52_dp, 0.37_dp, 0.21_dp, 0.09_dp]
      phys_response(:, 1) = [0.19_dp, -0.11_dp, 0.07_dp, -0.03_dp]
      phys_response(:, 2) = [-0.05_dp, 0.14_dp, -0.08_dp, 0.04_dp]
      density_response(:, 1) = [0.06_dp, -0.04_dp, 0.12_dp, -0.07_dp]
      density_response(:, 2) = [-0.08_dp, 0.03_dp, 0.09_dp, 0.02_dp]
      grad_e(:) = [0.15_dp, -0.08_dp, 0.05_dp, -0.04_dp]
      grad_rot(:) = [0.17_dp, -0.13_dp]
      rot_metric(:, :) = RESHAPE([0.74_dp, 0.06_dp, &
                                  0.06_dp, 0.66_dp], [nrot, nrot])

      h(:, :) = 0.0_dp
      DO i = 1, ne
         h(i, i) = occ_curv(i) + 0.03_dp*REAL(i, KIND=dp)
      END DO
      occ_sum = SUM(occ_curv)
      DO j = 1, ne
         DO i = 1, ne
            h(i, j) = h(i, j) - occ_curv(i)*occ_curv(j)/occ_sum
         END DO
      END DO
      h(:, :) = h(:, :) + MAXVAL(occ_curv(:))/REAL(ne, KIND=dp)
      h_ref(:, :) = h(:, :)

      DO s = 1, nrot
         shift = DOT_PRODUCT(occ_curv, phys_response(:, s))/occ_sum
         b_phys(:, s) = occ_curv(:)*(phys_response(:, s) - shift)
         shift = DOT_PRODUCT(occ_curv, density_response(:, s))/occ_sum
         b_density(:, s) = occ_curv(:)*(density_response(:, s) - shift)
      END DO
      b_total(:, :) = b_phys(:, :) + b_density(:, :)

      rhs_phys(:, :) = b_phys(:, :)
      h(:, :) = h_ref(:, :)
      CALL dgesv(ne, nrot, h, ne, ipiv, rhs_phys, ne, info)
      IF (info /= 0) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,I6)') TRIM(label), "failed in phys dgesv", info
         RETURN
      END IF

      rhs_density(:, :) = b_density(:, :)
      h(:, :) = h_ref(:, :)
      CALL dgesv(ne, nrot, h, ne, ipiv, rhs_density, ne, info)
      IF (info /= 0) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,I6)') TRIM(label), "failed in density dgesv", info
         RETURN
      END IF
      rhs_total(:, :) = rhs_phys(:, :) + rhs_density(:, :)

      z_base(:) = grad_e(:)
      h(:, :) = h_ref(:, :)
      CALL dgesv(ne, 1, h, ne, ipiv, z_base, ne, info)
      IF (info /= 0) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,I6)') TRIM(label), "failed in z dgesv", info
         RETURN
      END IF

      schur_total(:, :) = rot_metric(:, :) - MATMUL(TRANSPOSE(b_total), rhs_total)
      q_total(:) = grad_rot(:) + MATMUL(TRANSPOSE(rhs_total), grad_e)
      CALL dgesv(nrot, 1, schur_total, nrot, ipiv_rot, q_total, nrot, info)
      IF (info /= 0) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,I6)') TRIM(label), "failed in Schur dgesv", info
         RETURN
      END IF

      z_total(:) = z_base(:) + MATMUL(rhs_total, q_total)
      z_density_only(:) = z_base(:) + MATMUL(rhs_density, q_total)
      filter_weight(:) = [1.0_dp, 0.12_dp]
      q_filtered(:) = filter_weight(:)*q_total(:)
      z_filtered(:) = z_base(:) + MATMUL(rhs_total, q_filtered)

      aug(:, :) = 0.0_dp
      aug(1:ne, 1:ne) = h_ref(:, :)
      aug(1:ne, ne + 1:ntot) = -b_total(:, :)
      aug(ne + 1:ntot, 1:ne) = -TRANSPOSE(b_total(:, :))
      aug(ne + 1:ntot, ne + 1:ntot) = rot_metric(:, :)
      rhs_aug(1:ne) = grad_e(:)
      rhs_aug(ne + 1:ntot) = grad_rot(:)
      CALL dgesv(ntot, 1, aug, ntot, ipiv_aug, rhs_aug, ntot, info)
      IF (info /= 0) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,I6)') TRIM(label), "failed in augmented dgesv", info
         RETURN
      END IF
      z_aug(:) = rhs_aug(1:ne)
      q_aug(:) = rhs_aug(ne + 1:ntot)

      err = MAX(MAXVAL(ABS(z_total(:) - z_aug(:))), MAXVAL(ABS(q_total(:) - q_aug(:))))
      gap = MAXVAL(ABS(z_total(:) - z_density_only(:)))
      filter_gap = MAX(MAXVAL(ABS(z_total(:) - z_filtered(:))), &
                       MAXVAL(ABS(q_total(:) - q_filtered(:))))

      IF (err > tol .OR. gap <= 1.0E-5_dp .OR. filter_gap <= 1.0E-4_dp) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,3(1X,ES12.5))') TRIM(label), &
            "failed with error/omitted_gap/filter_gap", err, gap, filter_gap
      ELSE
         WRITE (*, '(A,1X,A,3(1X,ES12.5))') TRIM(label), &
            "OK with error/omitted_gap/filter_gap", err, gap, filter_gap
      END IF

   END SUBROUTINE check_projector_density_candidate_total_response

! **************************************************************************************************
!> \brief Check that Direct full-density history is calibrated against the accepted response.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_density_direct_history_mismatch(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      REAL(KIND=dp)                                      :: accepted_drop, corrected_drop, &
                                                            corrected_quality, err, mismatch_drop, &
                                                            naive_drop, naive_quality
      REAL(KIND=dp), DIMENSION(4)                        :: accepted_response, curv, model_response

      curv(:) = [0.45_dp, 0.31_dp, 0.18_dp, 0.07_dp]
      accepted_response(:) = [0.22_dp, -0.14_dp, 0.07_dp, -0.03_dp]
      model_response(:) = accepted_response(:) + [1.0E-3_dp, -2.0E-3_dp, 1.0E-3_dp, 0.0_dp]
      accepted_drop = 1.0E-3_dp
      naive_drop = 2.0E-3_dp

      mismatch_drop = 0.5_dp*SUM(curv(:)*(model_response(:) - accepted_response(:))**2)
      corrected_drop = MAX(0.0_dp, accepted_drop - mismatch_drop)
      naive_quality = accepted_drop/naive_drop
      corrected_quality = accepted_drop/corrected_drop
      err = ABS(corrected_quality - 1.0_dp)

      IF (err > 5.0E-3_dp .OR. naive_quality >= 0.6_dp .OR. &
          corrected_drop <= 0.95_dp*accepted_drop) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "failed with error/naive/corrected/drop", err, naive_quality, &
            corrected_quality, corrected_drop
      ELSE
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "OK with error/naive/corrected/drop", err, naive_quality, &
            corrected_quality, corrected_drop
      END IF

   END SUBROUTINE check_projector_density_direct_history_mismatch

! **************************************************************************************************
!> \brief Check that density-descent selection keeps its modeled rotation step.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_density_descent_selected_rotation(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: ne = 4, nrot = 2

      REAL(KIND=dp)                                      :: err, gap, model_gdotp, model_normp, &
                                                            scale, selected_gdotp, selected_normp, &
                                                            wrong_gdotp, wrong_normp
      REAL(KIND=dp), DIMENSION(ne)                       :: density_rx_e, grad_e, z_e, z_phys
      REAL(KIND=dp), DIMENSION(nrot)                     :: density_sol, grad_rot, physical_sol

      scale = 0.37_dp
      grad_e(:) = [0.21_dp, -0.17_dp, 0.09_dp, -0.06_dp]
      grad_rot(:) = [0.31_dp, -0.22_dp]
      z_phys(:) = [0.18_dp, -0.08_dp, 0.05_dp, -0.02_dp]
      density_rx_e(:) = [0.07_dp, -0.04_dp, 0.03_dp, -0.01_dp]
      physical_sol(:) = [0.26_dp, -0.15_dp]
      density_sol(:) = [-0.11_dp, 0.29_dp]
      z_e(:) = z_phys(:) + density_rx_e(:)

      model_gdotp = DOT_PRODUCT(grad_e, z_e) + DOT_PRODUCT(grad_rot, physical_sol)
      model_normp = DOT_PRODUCT(z_e, z_e) + DOT_PRODUCT(physical_sol, physical_sol)
      selected_gdotp = DOT_PRODUCT(grad_e, z_e) + DOT_PRODUCT(grad_rot, physical_sol)
      selected_normp = DOT_PRODUCT(z_e, z_e) + DOT_PRODUCT(physical_sol, physical_sol)
      wrong_gdotp = DOT_PRODUCT(grad_e, z_e) + DOT_PRODUCT(grad_rot, density_sol)
      wrong_normp = DOT_PRODUCT(z_e, z_e) + DOT_PRODUCT(density_sol, density_sol)

      err = MAX(ABS(scale*selected_gdotp - scale*model_gdotp), &
                ABS(scale*scale*selected_normp - scale*scale*model_normp))
      gap = MAX(ABS(scale*wrong_gdotp - scale*model_gdotp), &
                ABS(scale*scale*wrong_normp - scale*scale*model_normp))

      IF (err > tol .OR. gap <= 1.0E-3_dp) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,2(1X,ES12.5))') TRIM(label), "failed with error/wrong_gap", err, gap
      ELSE
         WRITE (*, '(A,1X,A,2(1X,ES12.5))') TRIM(label), "OK with error/wrong_gap", err, gap
      END IF

   END SUBROUTINE check_projector_density_descent_selected_rotation

! **************************************************************************************************
!> \brief Check the trust-limited Mermin quadratic used by density-response candidates.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_density_candidate_metric_model(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: nchannels = 1, ne = 5, nrot = 2

      INTEGER                                            :: info
      INTEGER, DIMENSION(nchannels)                      :: channel_nmo, e_offset
      INTEGER, DIMENSION(nrot)                           :: ipiv
      REAL(KIND=dp) :: curvature_scale, err, gdotp, indef_abs_norm, indef_signed_norm, &
         legacy_drop, metric_drop, norm_energy, norm_rotation, normp, response_trust, &
         soft_effective_norm, soft_euclid_norm, soft_metric_norm, solved_euclid_drop, &
         solved_euclid_norm, solved_metric_drop, solved_metric_norm, solved_metric_normp, &
         solved_scale, trust_scale, unit_weighted_step
      REAL(KIND=dp), DIMENSION(ne)                       :: grad_e, occ_curv, tether_curv, z_e
      REAL(KIND=dp), DIMENSION(nrot)                     :: grad_rot, rhs_rot, rhs_soft, rhs_solve, &
                                                            sol, sol_soft, sol_solve
      REAL(KIND=dp), DIMENSION(nrot, nrot)               :: schur_fact, schur_metric, schur_soft

      channel_nmo(:) = [ne]
      e_offset(:) = [1]
      occ_curv(:) = [0.47_dp, 0.23_dp, 0.11_dp, 0.07_dp, 0.03_dp]
      tether_curv(:) = [0.19_dp, 0.13_dp, 0.05_dp, 0.02_dp, 0.01_dp]
      z_e(:) = [0.31_dp, -0.22_dp, 0.18_dp, -0.07_dp, 0.11_dp]
      grad_e(:) = [0.42_dp, -0.26_dp, 0.19_dp, -0.13_dp, 0.09_dp]
      sol(:) = [0.34_dp, -0.21_dp]
      grad_rot(:) = [0.57_dp, -0.36_dp]
      schur_metric(:, :) = RESHAPE([1.80_dp, -0.14_dp, &
                                    -0.14_dp, 1.35_dp], [nrot, nrot])

      CALL dense_projected_mixed_metric(z_e, z_e, occ_curv, tether_curv, channel_nmo, e_offset, norm_energy)
      rhs_rot(:) = MATMUL(schur_metric, sol)
      norm_rotation = DOT_PRODUCT(rhs_rot, sol)
      normp = norm_energy + norm_rotation
      gdotp = DOT_PRODUCT(grad_e, z_e) + DOT_PRODUCT(grad_rot, sol)
      unit_weighted_step = 0.92_dp
      response_trust = 0.25_dp*unit_weighted_step
      trust_scale = MIN(1.0_dp, 0.8_dp*response_trust/unit_weighted_step)
      curvature_scale = gdotp/MAX(normp, EPSILON(normp))
      metric_drop = MIN(trust_scale, curvature_scale)*gdotp - &
                    0.5_dp*MIN(trust_scale, curvature_scale)**2*normp
      legacy_drop = 0.5_dp*MIN(trust_scale, curvature_scale)*gdotp

      err = ABS(norm_rotation - DOT_PRODUCT(sol, MATMUL(schur_metric, sol)))
      err = MAX(err, ABS(metric_drop - &
                         (MIN(trust_scale, curvature_scale)*gdotp - &
                          0.5_dp*MIN(trust_scale, curvature_scale)**2*normp)))

      rhs_solve(:) = [0.62_dp, -0.29_dp]
      schur_fact(:, :) = schur_metric(:, :)
      sol_solve(:) = rhs_solve(:)
      CALL dgesv(nrot, 1, schur_fact, nrot, ipiv, sol_solve, nrot, info)
      IF (info /= 0) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,I6)') TRIM(label), "failed in density-candidate dgesv", info
         RETURN
      END IF
      solved_metric_norm = DOT_PRODUCT(rhs_solve, sol_solve)
      solved_euclid_norm = DOT_PRODUCT(sol_solve, sol_solve)
      solved_metric_normp = norm_energy + solved_metric_norm
      gdotp = DOT_PRODUCT(grad_e, z_e) + DOT_PRODUCT(grad_rot, sol_solve)
      solved_scale = MIN(0.45_dp, gdotp/MAX(solved_metric_normp, EPSILON(solved_metric_normp)))
      solved_metric_drop = solved_scale*gdotp - 0.5_dp*solved_scale*solved_scale* &
                           solved_metric_normp
      solved_euclid_drop = solved_scale*gdotp - 0.5_dp*solved_scale*solved_scale* &
                           (norm_energy + solved_euclid_norm)
      err = MAX(err, ABS(solved_metric_norm - &
                         DOT_PRODUCT(sol_solve, MATMUL(schur_metric, sol_solve))))

      schur_soft(:, :) = RESHAPE([0.08_dp, 0.00_dp, &
                                  0.00_dp, 0.04_dp], [nrot, nrot])
      sol_soft(:) = [1.0_dp, -0.5_dp]
      rhs_soft(:) = MATMUL(schur_soft, sol_soft)
      soft_metric_norm = DOT_PRODUCT(rhs_soft, sol_soft)
      soft_euclid_norm = DOT_PRODUCT(sol_soft, sol_soft)
      soft_effective_norm = MAX(soft_metric_norm, soft_euclid_norm)
      err = MAX(err, ABS(soft_effective_norm - soft_euclid_norm))

      schur_soft(:, :) = RESHAPE([-0.08_dp, 0.00_dp, &
                                  0.00_dp, 0.04_dp], [nrot, nrot])
      sol_soft(:) = [1.0_dp, 0.5_dp]
      rhs_soft(:) = MATMUL(schur_soft, sol_soft)
      indef_signed_norm = DOT_PRODUCT(rhs_soft, sol_soft)
      indef_abs_norm = ABS(schur_soft(1, 1))*sol_soft(1)**2 + &
                       ABS(schur_soft(2, 2))*sol_soft(2)**2
      err = MAX(err, ABS(indef_abs_norm - 9.0E-2_dp))

      IF (err > tol .OR. trust_scale >= curvature_scale .OR. metric_drop <= legacy_drop .OR. &
          norm_rotation <= DOT_PRODUCT(sol, sol) .OR. &
          solved_metric_norm <= solved_euclid_norm .OR. &
          solved_euclid_drop <= solved_metric_drop .OR. &
          soft_metric_norm >= soft_euclid_norm .OR. indef_signed_norm >= 0.0_dp .OR. &
          indef_abs_norm <= ABS(indef_signed_norm)) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,11(1X,ES12.5))') TRIM(label), &
            "failed with err/scale/curv/drop/legacy/rnorm/snorm/sdrop/isgn/iabs", err, trust_scale, &
            curvature_scale, metric_drop, legacy_drop, norm_rotation, solved_metric_norm, &
            solved_metric_drop, indef_signed_norm, indef_abs_norm
      ELSE
         WRITE (*, '(A,1X,A,11(1X,ES12.5))') TRIM(label), &
            "OK with err/scale/curv/drop/legacy/rnorm/snorm/sdrop/isgn/iabs", err, trust_scale, &
            curvature_scale, metric_drop, legacy_drop, norm_rotation, solved_metric_norm, &
            solved_metric_drop, indef_signed_norm, indef_abs_norm
      END IF

   END SUBROUTINE check_projector_density_candidate_metric_model

! **************************************************************************************************
!> \brief Check exact metric projection of coupled fixed-N and local tether response gauges.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_metric_response_projection(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: nchannels = 2, ne = 6

      INTEGER, DIMENSION(nchannels)                      :: channel_nmo, e_offset
      LOGICAL                                            :: metric_projected, reference_projected, &
                                                            single_projected
      REAL(KIND=dp)                                      :: err, legacy_residual, reference_err, &
                                                            single_err, tether_residual
      REAL(KIND=dp), DIMENSION(ne)                       :: legacy, legacy_occ_only, metric, &
                                                            metric_occ_only, metric_reference, &
                                                            occ_curv, response, tether_curv, &
                                                            tether_zero

      channel_nmo(:) = [3, 3]
      e_offset(:) = [1, 4]
      occ_curv(:) = [0.90_dp, 0.20_dp, 0.05_dp, 0.70_dp, 0.12_dp, 0.03_dp]
      tether_curv(:) = [0.40_dp, 0.05_dp, 0.02_dp, 0.25_dp, 0.14_dp, 0.03_dp]
      tether_zero(:) = 0.0_dp
      response(:) = [0.40_dp, -0.20_dp, 0.70_dp, -0.35_dp, 0.18_dp, 0.60_dp]

      legacy(:) = response(:)
      metric(:) = response(:)
      metric_reference(:) = response(:)
      CALL dense_project_response_two_pass(legacy, occ_curv, tether_curv, channel_nmo, e_offset)
      CALL dense_project_response_metric_distributed(metric, occ_curv, tether_curv, channel_nmo, e_offset, &
                                                     metric_projected)
      CALL dense_project_response_metric(metric_reference, occ_curv, tether_curv, channel_nmo, e_offset, &
                                         reference_projected)
      err = ABS(DOT_PRODUCT(occ_curv(:), metric(:)))
      tether_residual = dense_response_tether_residual(metric, tether_curv, channel_nmo, e_offset)
      err = MAX(err, tether_residual)
      reference_err = MAXVAL(ABS(metric(:) - metric_reference(:)))
      legacy_residual = ABS(DOT_PRODUCT(occ_curv(:), legacy(:)))

      legacy_occ_only(:) = response(:)
      metric_occ_only(:) = response(:)
      CALL dense_project_response_two_pass(legacy_occ_only, occ_curv, tether_zero, channel_nmo, e_offset)
      CALL dense_project_response_metric_distributed(metric_occ_only, occ_curv, tether_zero, channel_nmo, e_offset, &
                                                     single_projected)
      single_err = MAXVAL(ABS(metric_occ_only(:) - legacy_occ_only(:)))
      single_err = MAX(single_err, ABS(DOT_PRODUCT(occ_curv(:), metric_occ_only(:))))

      IF (.NOT. metric_projected .OR. .NOT. reference_projected .OR. .NOT. single_projected .OR. &
          err > tol .OR. reference_err > tol .OR. single_err > tol .OR. legacy_residual < 1.0E-2_dp) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,5(1X,ES12.5))') TRIM(label), &
            "failed with error/tether/legacy/single/reference", err, tether_residual, legacy_residual, &
            single_err, reference_err
      ELSE
         WRITE (*, '(A,1X,A,5(1X,ES12.5))') TRIM(label), &
            "OK with error/tether/legacy/single/reference", err, tether_residual, legacy_residual, &
            single_err, reference_err
      END IF

   END SUBROUTINE check_projector_metric_response_projection

! **************************************************************************************************
!> \brief Check that full-density response keeps global fixed-N k-point charge transfer.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_density_global_fixed_n_projection(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: nchannels = 2, ne = 4

      INTEGER                                            :: first_energy, local_channel, nmo
      INTEGER, DIMENSION(nchannels)                      :: channel_nmo, e_offset
      LOGICAL                                            :: projected
      REAL(KIND=dp)                                      :: charge, err, local_gap, local_norm, &
                                                            shift, weight_sum
      REAL(KIND=dp), DIMENSION(ne)                       :: global_projected, local_projected, &
                                                            occ_curv, raw_response, tether_curv

      channel_nmo(:) = [2, 2]
      e_offset(:) = [1, 3]
      occ_curv(:) = [0.60_dp, 0.40_dp, 0.25_dp, 0.75_dp]
      tether_curv(:) = 0.0_dp
      raw_response(:) = [1.0_dp, 1.0_dp, -1.0_dp, -1.0_dp]

      global_projected(:) = raw_response(:)
      CALL dense_project_response_metric_distributed(global_projected, occ_curv, tether_curv, &
                                                     channel_nmo, e_offset, projected)

      local_projected(:) = raw_response(:)
      DO local_channel = 1, nchannels
         nmo = channel_nmo(local_channel)
         first_energy = e_offset(local_channel)
         weight_sum = SUM(occ_curv(first_energy:first_energy + nmo - 1))
         shift = DOT_PRODUCT(occ_curv(first_energy:first_energy + nmo - 1), &
                             local_projected(first_energy:first_energy + nmo - 1))/weight_sum
         local_projected(first_energy:first_energy + nmo - 1) = &
            local_projected(first_energy:first_energy + nmo - 1) - shift
      END DO

      err = MAXVAL(ABS(global_projected(:) - raw_response(:)))
      charge = ABS(DOT_PRODUCT(occ_curv(:), global_projected(:)))
      local_gap = MAXVAL(ABS(raw_response(:) - local_projected(:)))
      local_norm = MAXVAL(ABS(local_projected(:)))

      IF (.NOT. projected .OR. err > tol .OR. charge > tol .OR. local_gap < 0.5_dp .OR. &
          local_norm > tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "failed with error/charge/localgap/localnorm", err, charge, local_gap, local_norm
      ELSE
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "OK with error/charge/localgap/localnorm", err, charge, local_gap, local_norm
      END IF

   END SUBROUTINE check_projector_density_global_fixed_n_projection

! **************************************************************************************************
!> \brief Check that density sources entering the exact response are fixed-N/tether projected first.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_density_exact_source_projection(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: nchannels = 2, ne = 6, nrot = 2

      INTEGER                                            :: info, r
      INTEGER, DIMENSION(nchannels)                      :: channel_nmo, e_offset
      INTEGER, DIMENSION(ne)                             :: ipiv
      LOGICAL                                            :: projected
      REAL(KIND=dp)                                      :: err, exact_gap, exact_residual, raw_occ_residual, &
                                                            raw_tether_residual, source_gap, source_residual
      REAL(KIND=dp), DIMENSION(ne)                       :: occ_curv, tether_curv
      REAL(KIND=dp), DIMENSION(ne, ne)                   :: h_projected, h_raw
      REAL(KIND=dp), DIMENSION(ne, nrot)                 :: b_projected, b_raw, exact_projected, exact_raw, &
                                                            projected_density, raw_density, zero_tether

      channel_nmo(:) = [3, 3]
      e_offset(:) = [1, 4]
      occ_curv(:) = [0.31_dp, 0.74_dp, 0.16_dp, 0.0_dp, 0.0_dp, 0.0_dp]
      tether_curv(:) = [0.05_dp, 0.02_dp, 0.08_dp, 0.07_dp, 0.03_dp, 0.0_dp]
      raw_density(:, 1) = [0.31_dp, -0.17_dp, 0.22_dp, 0.48_dp, -0.36_dp, 0.18_dp]
      raw_density(:, 2) = [-0.22_dp, 0.19_dp, -0.11_dp, -0.40_dp, 0.27_dp, -0.15_dp]
      zero_tether(:, :) = 0.0_dp

      projected_density(:, :) = raw_density(:, :)
      source_residual = 0.0_dp
      raw_occ_residual = 0.0_dp
      raw_tether_residual = 0.0_dp
      DO r = 1, nrot
         raw_occ_residual = MAX(raw_occ_residual, ABS(DOT_PRODUCT(occ_curv(:), raw_density(:, r))))
         raw_tether_residual = MAX(raw_tether_residual, &
                                   dense_response_tether_residual(raw_density(:, r), tether_curv, &
                                                                  channel_nmo, e_offset))
         CALL dense_project_response_metric_distributed(projected_density(:, r), occ_curv, tether_curv, &
                                                        channel_nmo, e_offset, projected)
         IF (.NOT. projected) source_residual = HUGE(source_residual)
         source_residual = MAX(source_residual, ABS(DOT_PRODUCT(occ_curv(:), projected_density(:, r))))
         source_residual = MAX(source_residual, &
                               dense_response_tether_residual(projected_density(:, r), tether_curv, &
                                                              channel_nmo, e_offset))
      END DO
      source_gap = MAXVAL(ABS(projected_density(:, :) - raw_density(:, :)))

      CALL build_test_mixed_energy_block(channel_nmo, e_offset, raw_density, zero_tether, h_raw, b_raw)
      CALL build_test_mixed_energy_block(channel_nmo, e_offset, projected_density, zero_tether, &
                                         h_projected, b_projected)
      exact_raw(:, :) = b_raw(:, :)
      exact_projected(:, :) = b_projected(:, :)
      CALL dgesv(ne, nrot, h_raw, ne, ipiv, exact_raw, ne, info)
      IF (info /= 0) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,I6)') TRIM(label), "failed raw dgesv", info
         RETURN
      END IF
      CALL dgesv(ne, nrot, h_projected, ne, ipiv, exact_projected, ne, info)
      IF (info /= 0) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,I6)') TRIM(label), "failed projected dgesv", info
         RETURN
      END IF

      exact_residual = 0.0_dp
      DO r = 1, nrot
         CALL dense_project_response_metric_distributed(exact_raw(:, r), occ_curv, tether_curv, &
                                                        channel_nmo, e_offset, projected)
         CALL dense_project_response_metric_distributed(exact_projected(:, r), occ_curv, tether_curv, &
                                                        channel_nmo, e_offset, projected)
         exact_residual = MAX(exact_residual, ABS(DOT_PRODUCT(occ_curv(:), exact_projected(:, r))))
         exact_residual = MAX(exact_residual, &
                              dense_response_tether_residual(exact_projected(:, r), tether_curv, &
                                                             channel_nmo, e_offset))
      END DO
      exact_gap = MAXVAL(ABS(exact_projected(:, :) - exact_raw(:, :)))
      err = MAX(source_residual, exact_residual)

      IF (err > tol .OR. raw_occ_residual < 1.0E-2_dp .OR. raw_tether_residual < 1.0E-3_dp .OR. &
          source_gap < 1.0E-2_dp .OR. exact_gap < 1.0E-5_dp) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,5(1X,ES12.5))') TRIM(label), &
            "failed with err/rawocc/rawtether/sourcegap/exactgap", err, raw_occ_residual, &
            raw_tether_residual, source_gap, exact_gap
      ELSE
         WRITE (*, '(A,1X,A,5(1X,ES12.5))') TRIM(label), &
            "OK with err/rawocc/rawtether/sourcegap/exactgap", err, raw_occ_residual, &
            raw_tether_residual, source_gap, exact_gap
      END IF

   END SUBROUTINE check_projector_density_exact_source_projection

! **************************************************************************************************
!> \brief Check that floor-limited density response is suitable for global fixed-N projection.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_density_floor_global_projection(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: nchannels = 2, ne = 4
      REAL(KIND=dp), PARAMETER :: response_ratio_gate = 2.0E-2_dp

      INTEGER, DIMENSION(nchannels)                      :: channel_nmo, e_offset
      LOGICAL                                            :: floor_projected, raw_projected_ok
      REAL(KIND=dp)                                      :: charge, floor_norm, raw_norm, retained
      REAL(KIND=dp), DIMENSION(ne)                       :: floor_response, occ_curv, &
                                                            projector_response, raw_response, &
                                                            safe_response, tether_curv

      channel_nmo(:) = [2, 2]
      e_offset(:) = [1, 3]
      occ_curv(:) = [4.0E-2_dp, 1.2E-2_dp, 4.0E-4_dp, 4.0E-8_dp]
      tether_curv(:) = 0.0_dp
      projector_response(:) = [0.010_dp, -0.008_dp, 0.006_dp, -0.004_dp]
      safe_response(:) = [0.10_dp, -0.08_dp, 0.06_dp, -25.0_dp]

      floor_response(:) = 0.0_dp
      raw_response(:) = safe_response(:)
      floor_response(1) = dense_physical_density_floor_audit_response( &
                          projector_response(1), safe_response(1), occ_curv(1), MAXVAL(occ_curv(:)), response_ratio_gate)
      floor_response(2) = dense_physical_density_floor_audit_response( &
                          projector_response(2), safe_response(2), occ_curv(2), MAXVAL(occ_curv(:)), response_ratio_gate)
      floor_response(3) = dense_physical_density_floor_audit_response( &
                          projector_response(3), safe_response(3), occ_curv(3), MAXVAL(occ_curv(:)), response_ratio_gate)
      floor_response(4) = dense_physical_density_floor_audit_response( &
                          projector_response(4), safe_response(4), occ_curv(4), MAXVAL(occ_curv(:)), response_ratio_gate)

      CALL dense_project_response_metric_distributed(floor_response, occ_curv, tether_curv, &
                                                     channel_nmo, e_offset, floor_projected)
      CALL dense_project_response_metric_distributed(raw_response, occ_curv, tether_curv, &
                                                     channel_nmo, e_offset, raw_projected_ok)

      charge = ABS(DOT_PRODUCT(occ_curv(:), floor_response(:)))
      floor_norm = MAXVAL(ABS(floor_response(:)))
      raw_norm = MAXVAL(ABS(raw_response(:)))
      retained = ABS(floor_response(3))/MAX(ABS(projector_response(3)), EPSILON(1.0_dp))

      IF (.NOT. floor_projected .OR. .NOT. raw_projected_ok .OR. charge > tol .OR. &
          floor_norm >= 0.1_dp*raw_norm .OR. retained < 1.0E-3_dp) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "failed with charge/floor/raw/retained", charge, floor_norm, raw_norm, retained
      ELSE
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "OK with charge/floor/raw/retained", charge, floor_norm, raw_norm, retained
      END IF

   END SUBROUTINE check_projector_density_floor_global_projection

! **************************************************************************************************
!> \brief Check that Direct full-density calibration uses floor-limited density response.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_density_direct_floor_response(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: nchannels = 2, ne = 4, nrot = 2
      REAL(KIND=dp), PARAMETER :: response_ratio_gate = 2.0E-2_dp

      INTEGER                                            :: r
      INTEGER, DIMENSION(nchannels)                      :: channel_nmo, e_offset
      LOGICAL                                            :: column_projected, floor_projected, &
                                                            raw_projected
      REAL(KIND=dp)                                      :: charge, floor_abs, raw_abs, response_gap
      REAL(KIND=dp), DIMENSION(ne)                       :: floor_model, occ_curv, raw_model, &
                                                            tether_curv
      REAL(KIND=dp), DIMENSION(ne, nrot)                 :: floor_density, phys_response, &
                                                            projector_response, raw_density
      REAL(KIND=dp), DIMENSION(nrot)                     :: direct_sol

      channel_nmo(:) = [2, 2]
      e_offset(:) = [1, 3]
      occ_curv(:) = [3.0E-1_dp, 1.5E-1_dp, 5.0E-4_dp, 5.0E-8_dp]
      tether_curv(:) = 0.0_dp
      direct_sol(:) = [0.70_dp, -0.40_dp]
      phys_response(:, 1) = [0.04_dp, -0.02_dp, 0.01_dp, -0.01_dp]
      phys_response(:, 2) = [-0.03_dp, 0.01_dp, -0.02_dp, 0.02_dp]
      projector_response(:, 1) = [0.015_dp, -0.012_dp, 0.006_dp, -0.005_dp]
      projector_response(:, 2) = [-0.014_dp, 0.010_dp, -0.004_dp, 0.003_dp]
      raw_density(:, 1) = [0.02_dp, -0.015_dp, 25.0_dp, -20.0_dp]
      raw_density(:, 2) = [-0.03_dp, 0.010_dp, -18.0_dp, 15.0_dp]

      floor_projected = .TRUE.
      raw_projected = .TRUE.
      DO r = 1, nrot
         floor_density(1, r) = dense_physical_density_floor_audit_response( &
                               projector_response(1, r), raw_density(1, r), &
                               occ_curv(1), MAXVAL(occ_curv(:)), response_ratio_gate)
         floor_density(2, r) = dense_physical_density_floor_audit_response( &
                               projector_response(2, r), raw_density(2, r), &
                               occ_curv(2), MAXVAL(occ_curv(:)), response_ratio_gate)
         floor_density(3, r) = dense_physical_density_floor_audit_response( &
                               projector_response(3, r), raw_density(3, r), &
                               occ_curv(3), MAXVAL(occ_curv(:)), response_ratio_gate)
         floor_density(4, r) = dense_physical_density_floor_audit_response( &
                               projector_response(4, r), raw_density(4, r), &
                               occ_curv(4), MAXVAL(occ_curv(:)), response_ratio_gate)
         CALL dense_project_response_metric_distributed(floor_density(:, r), occ_curv, tether_curv, &
                                                        channel_nmo, e_offset, column_projected)
         floor_projected = floor_projected .AND. column_projected
         CALL dense_project_response_metric_distributed(raw_density(:, r), occ_curv, tether_curv, &
                                                        channel_nmo, e_offset, column_projected)
         raw_projected = raw_projected .AND. column_projected
      END DO

      floor_model(:) = MATMUL(phys_response, direct_sol) + MATMUL(floor_density, direct_sol)
      raw_model(:) = MATMUL(phys_response, direct_sol) + MATMUL(raw_density, direct_sol)
      charge = ABS(DOT_PRODUCT(occ_curv(:), MATMUL(floor_density, direct_sol)))
      floor_abs = MAXVAL(ABS(floor_model(:)))
      raw_abs = MAXVAL(ABS(raw_model(:)))
      response_gap = raw_abs/MAX(floor_abs, EPSILON(floor_abs))

      IF (.NOT. floor_projected .OR. .NOT. raw_projected .OR. charge > tol .OR. &
          floor_abs >= 0.1_dp*raw_abs .OR. response_gap < 20.0_dp) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "failed with charge/floor/raw/gap", charge, floor_abs, raw_abs, response_gap
      ELSE
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "OK with charge/floor/raw/gap", charge, floor_abs, raw_abs, response_gap
      END IF

   END SUBROUTINE check_projector_density_direct_floor_response

! **************************************************************************************************
!> \brief Check the conservative Direct floor-density bootstrap gate.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_density_floor_bootstrap_weight(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      REAL(KIND=dp) :: active_source, boot, err, expected_boot, expected_over, high_lowcurv, &
         high_residual, low_mermin, no_exact, no_floor, no_samples, over_mermin

      boot = dense_direct_floor_density_bootstrap_weight(5.0E-6_dp, 0.96_dp, 6, &
                                                         2.0E-5_dp, 1.0E-6_dp, &
                                                         2.0E-2_dp, 0.0_dp)
      over_mermin = dense_direct_floor_density_bootstrap_weight(5.0E-6_dp, 1.50_dp, 6, &
                                                                2.0E-5_dp, 1.0E-6_dp, &
                                                                2.0E-2_dp, 0.0_dp)
      high_residual = dense_direct_floor_density_bootstrap_weight(2.0E-5_dp, 0.96_dp, 6, &
                                                                  2.0E-5_dp, 1.0E-6_dp, &
                                                                  2.0E-2_dp, 0.0_dp)
      no_samples = dense_direct_floor_density_bootstrap_weight(5.0E-6_dp, 0.96_dp, 3, &
                                                               2.0E-5_dp, 1.0E-6_dp, &
                                                               2.0E-2_dp, 0.0_dp)
      low_mermin = dense_direct_floor_density_bootstrap_weight(5.0E-6_dp, 0.70_dp, 6, &
                                                               2.0E-5_dp, 1.0E-6_dp, &
                                                               2.0E-2_dp, 0.0_dp)
      high_lowcurv = dense_direct_floor_density_bootstrap_weight(5.0E-6_dp, 0.96_dp, 6, &
                                                                 2.0E-5_dp, 1.0E-6_dp, &
                                                                 3.0E-1_dp, 0.0_dp)
      active_source = dense_direct_floor_density_bootstrap_weight(5.0E-6_dp, 0.96_dp, 6, &
                                                                  2.0E-5_dp, 1.0E-6_dp, &
                                                                  2.0E-2_dp, 2.0E-5_dp)
      no_floor = dense_direct_floor_density_bootstrap_weight(5.0E-6_dp, 0.96_dp, 6, &
                                                             0.0_dp, 1.0E-6_dp, &
                                                             2.0E-2_dp, 0.0_dp)
      no_exact = dense_direct_floor_density_bootstrap_weight(5.0E-6_dp, 0.96_dp, 6, &
                                                             2.0E-5_dp, 0.0_dp, &
                                                             2.0E-2_dp, 0.0_dp)

      expected_boot = 2.5E-2_dp*0.96_dp*(1.0_dp - 2.0E-2_dp/2.5E-1_dp)
      expected_over = 2.5E-2_dp*(1.0_dp - 2.0E-2_dp/2.5E-1_dp)
      err = ABS(boot - expected_boot)
      err = MAX(err, ABS(over_mermin - expected_over))

      IF (err > tol .OR. boot <= 0.0_dp .OR. over_mermin <= boot .OR. &
          high_residual > tol .OR. no_samples > tol .OR. low_mermin > tol .OR. &
          high_lowcurv > tol .OR. active_source > tol .OR. no_floor > tol .OR. &
          no_exact > tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,10(1X,ES12.5))') TRIM(label), &
            "failed with err/boot/over/hr/ns/lm/hl/as/nf/ne", err, boot, over_mermin, &
            high_residual, no_samples, low_mermin, high_lowcurv, active_source, no_floor, no_exact
      ELSE
         WRITE (*, '(A,1X,A,10(1X,ES12.5))') TRIM(label), &
            "OK with err/boot/over/hr/ns/lm/hl/as/nf/ne", err, boot, over_mermin, &
            high_residual, no_samples, low_mermin, high_lowcurv, active_source, no_floor, no_exact
      END IF

   END SUBROUTINE check_projector_density_floor_bootstrap_weight

! **************************************************************************************************
!> \brief Check that global fixed-N density-response trust follows the Mermin metric.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_density_global_metric_trust(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      REAL(KIND=dp) :: abs_step, abs_trust, abs_trust_ratio, err, legacy_quality, metric_quality, &
         response_trust, weighted_step, weighted_trust_ratio
      REAL(KIND=dp), DIMENSION(4)                        :: curv, response

      curv(:) = [1.0E-10_dp, 0.40_dp, 0.40_dp, 0.20_dp]
      response(:) = [2.0E4_dp, 0.20_dp, -0.15_dp, 0.05_dp]
      response_trust = 0.50_dp
      abs_trust = response_trust

      weighted_step = SQRT(SUM(curv(:)*response(:)**2)/SUM(curv(:)))
      abs_step = MAXVAL(ABS(response(:)))
      weighted_trust_ratio = weighted_step/response_trust
      abs_trust_ratio = abs_step/abs_trust
      legacy_quality = 1.0_dp/MAX(1.0_dp, weighted_trust_ratio, abs_trust_ratio)
      metric_quality = 1.0_dp/MAX(1.0_dp, weighted_trust_ratio)

      err = ABS(metric_quality - 1.0_dp)
      IF (err > tol .OR. weighted_trust_ratio >= 1.0_dp .OR. &
          abs_trust_ratio <= 1.0E4_dp .OR. legacy_quality >= 1.0E-3_dp) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,5(1X,ES12.5))') TRIM(label), &
            "failed with err/wtr/atr/metric/legacy", err, weighted_trust_ratio, &
            abs_trust_ratio, metric_quality, legacy_quality
      ELSE
         WRITE (*, '(A,1X,A,5(1X,ES12.5))') TRIM(label), &
            "OK with err/wtr/atr/metric/legacy", err, weighted_trust_ratio, &
            abs_trust_ratio, metric_quality, legacy_quality
      END IF

   END SUBROUTINE check_projector_density_global_metric_trust

! **************************************************************************************************
!> \brief Check that low-curvature density-audit trust measures consistency, not amplification.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_density_audit_trust_quality(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      REAL(KIND=dp) :: consistent_quality, err, expected_huge_quality, expected_low_quality, &
         expected_over_quality, huge_quality, inactive_quality, low_amp_over_quality, low_quality, &
         over_quality, young_quality

      young_quality = dense_density_audit_trust_quality(0.95_dp, 0.95_dp, &
                                                        10.0_dp, 10.0_dp, 7)
      inactive_quality = dense_density_audit_trust_quality(0.20_dp, 0.20_dp, &
                                                           10.0_dp, 10.0_dp, 16)
      consistent_quality = dense_density_audit_trust_quality(0.80_dp, 0.80_dp, &
                                                             8.0_dp, 8.0_dp, 16)
      low_quality = dense_density_audit_trust_quality(0.80_dp, 0.16_dp, &
                                                      8.0_dp, 1.60_dp, 16)
      over_quality = dense_density_audit_trust_quality(0.999_dp, 52.947_dp, &
                                                       53.0_dp, 2809.0_dp, 32)
      low_amp_over_quality = dense_density_audit_trust_quality(0.999_dp, 52.947_dp, &
                                                               0.112_dp, 4.704_dp, 32)
      huge_quality = dense_density_audit_trust_quality(1.0_dp, 1.0E4_dp, &
                                                       1.0E6_dp, 1.0E12_dp, 32)

      expected_low_quality = 1.0_dp - ((0.80_dp - 0.25_dp)/0.75_dp)* &
                             (LOG(8.0_dp)/LOG(1.0E6_dp))*(1.0_dp - 0.20_dp)
      expected_over_quality = 1.0_dp - ((0.999_dp - 0.25_dp)/0.75_dp)* &
                              (LOG(53.0_dp)/LOG(1.0E6_dp))*(1.0_dp - 0.10_dp)
      expected_huge_quality = 0.10_dp

      err = ABS(young_quality - 1.0_dp)
      err = MAX(err, ABS(inactive_quality - 1.0_dp))
      err = MAX(err, ABS(consistent_quality - 1.0_dp))
      err = MAX(err, ABS(low_quality - expected_low_quality))
      err = MAX(err, ABS(over_quality - expected_over_quality))
      err = MAX(err, ABS(low_amp_over_quality - expected_over_quality))
      err = MAX(err, ABS(huge_quality - expected_huge_quality))

      IF (err > tol .OR. low_quality >= consistent_quality .OR. &
          over_quality >= low_quality .OR. &
          ABS(low_amp_over_quality - over_quality) > tol .OR. huge_quality >= over_quality) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,8(1X,ES12.5))') TRIM(label), &
            "failed with error/young/inactive/consistent/low/over/lowamp/huge", &
            err, young_quality, inactive_quality, consistent_quality, low_quality, &
            over_quality, low_amp_over_quality, huge_quality
      ELSE
         WRITE (*, '(A,1X,A,8(1X,ES12.5))') TRIM(label), &
            "OK with error/young/inactive/consistent/low/over/lowamp/huge", &
            err, young_quality, inactive_quality, consistent_quality, low_quality, &
            over_quality, low_amp_over_quality, huge_quality
      END IF

   END SUBROUTINE check_projector_density_audit_trust_quality

! **************************************************************************************************
!> \brief Check that Direct trace compensation follows the selected delta model.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_direct_trace_compensation(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      REAL(KIND=dp), PARAMETER                           :: active_delta = -0.027_dp, &
                                                            coupled_weight = 0.20_dp, &
                                                            full_delta = -0.181_dp, &
                                                            reference_delta = 0.0_dp

      INTEGER                                            :: i
      REAL(KIND=dp) :: active_comp, active_diag_after, active_negative_ref_comp, err, full_comp, &
         old_active_diag_after, positive_comp, scaled_active_after, scaled_active_comp, &
         trace_scale, uniform_after, uniform_comp
      REAL(KIND=dp), DIMENSION(4)                        :: row_after, row_comp, row_delta, &
                                                            row_reference

      full_comp = dense_direct_trace_compensation(coupled_weight, full_delta, reference_delta)
      active_comp = dense_direct_trace_compensation(coupled_weight, active_delta, reference_delta)
      active_diag_after = active_delta + active_comp/coupled_weight
      old_active_diag_after = active_delta + full_comp/coupled_weight
      active_negative_ref_comp = dense_direct_trace_compensation(coupled_weight, active_delta, full_delta)
      positive_comp = dense_direct_trace_compensation(coupled_weight, 0.05_dp, 0.01_dp)
      trace_scale = dense_trace_diag_component_scale(-1.0_dp, 0.75_dp, 1.0_dp, 0.8_dp, 16)
      scaled_active_comp = trace_scale*active_comp
      scaled_active_after = active_delta + scaled_active_comp/coupled_weight
      row_delta(:) = [-0.20_dp, -0.04_dp, 0.005_dp, 0.06_dp]
      row_reference(:) = [0.0_dp, -0.10_dp, 0.01_dp, 0.10_dp]
      DO i = 1, SIZE(row_delta)
         row_comp(i) = dense_direct_trace_compensation(coupled_weight, row_delta(i), row_reference(i))
         row_after(i) = row_delta(i) + row_comp(i)/coupled_weight
      END DO
      uniform_comp = dense_direct_trace_compensation(coupled_weight, &
                                                     SUM(row_delta(:))/REAL(SIZE(row_delta), KIND=dp), &
                                                     SUM(row_reference(:))/REAL(SIZE(row_reference), KIND=dp))
      uniform_after = row_delta(3) + uniform_comp/coupled_weight

      err = ABS(active_diag_after)
      err = MAX(err, ABS(full_delta + full_comp/coupled_weight))
      err = MAX(err, ABS(active_negative_ref_comp - active_comp))
      err = MAX(err, ABS(positive_comp - coupled_weight*(0.05_dp - 0.01_dp)))
      err = MAX(err, ABS(scaled_active_comp - trace_scale*active_comp))
      err = MAX(err, ABS(row_after(1)), ABS(row_after(2)), ABS(row_after(3) - row_delta(3)))
      err = MAX(err, ABS(row_after(4) - row_delta(4)))

      IF (err > tol .OR. old_active_diag_after < 0.1_dp .OR. active_comp >= full_comp .OR. &
          scaled_active_after >= active_diag_after .OR. scaled_active_after <= active_delta .OR. &
          uniform_after <= row_delta(3)) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,9(1X,ES12.5))') TRIM(label), &
            "failed with error/full/active/negref/old/positive/scaled/row/uniform", err, full_comp, &
            active_comp, active_negative_ref_comp, old_active_diag_after, positive_comp, &
            scaled_active_after, row_comp(3), uniform_after
      ELSE
         WRITE (*, '(A,1X,A,9(1X,ES12.5))') TRIM(label), &
            "OK with error/full/active/negref/old/positive/scaled/row/uniform", err, full_comp, active_comp, &
            active_negative_ref_comp, old_active_diag_after, positive_comp, scaled_active_after, row_comp(3), &
            uniform_after
      END IF

   END SUBROUTINE check_projector_direct_trace_compensation

! **************************************************************************************************
!> \brief Check that diagonal-dominated Direct traces do not masquerade as failed split histories.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_trace_split_history_quality(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      REAL(KIND=dp) :: balanced_quality, diag_tail_quality, err, fallback_quality, &
         no_samples_quality, poor_tail_quality, two_component_quality

      fallback_quality = 0.42_dp
      no_samples_quality = dense_trace_split_history_quality(1.0_dp, -0.4_dp, &
                                                             0.5_dp, 0.2_dp, 2, &
                                                             fallback_quality)
      diag_tail_quality = dense_trace_split_history_quality(1.0_dp, -1.0E-3_dp, &
                                                            0.05_dp, 0.0_dp, 16, &
                                                            fallback_quality)
      poor_tail_quality = dense_trace_split_history_quality(1.0_dp, -1.0E-3_dp, &
                                                            0.02_dp, 0.0_dp, 16, &
                                                            0.0_dp)
      balanced_quality = dense_trace_split_history_quality(0.75_dp, -0.25_dp, &
                                                           0.60_dp, 0.10_dp, 16, &
                                                           fallback_quality)
      two_component_quality = dense_trace_split_history_quality(0.75_dp, 0.25_dp, &
                                                                0.60_dp, 0.10_dp, 16, &
                                                                fallback_quality)

      err = ABS(no_samples_quality - fallback_quality)
      err = MAX(err, ABS(diag_tail_quality - fallback_quality))
      err = MAX(err, ABS(poor_tail_quality - 0.10_dp))
      err = MAX(err, ABS(balanced_quality - 0.35_dp))
      err = MAX(err, ABS(two_component_quality - 0.10_dp))

      IF (err > tol .OR. diag_tail_quality <= balanced_quality .OR. &
          two_component_quality >= balanced_quality) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,6(1X,ES12.5))') TRIM(label), &
            "failed with error/no/diag/poor/bal/two", err, no_samples_quality, &
            diag_tail_quality, poor_tail_quality, balanced_quality, two_component_quality
      ELSE
         WRITE (*, '(A,1X,A,6(1X,ES12.5))') TRIM(label), &
            "OK with error/no/diag/poor/bal/two", err, no_samples_quality, &
            diag_tail_quality, poor_tail_quality, balanced_quality, two_component_quality
      END IF

   END SUBROUTINE check_projector_trace_split_history_quality

! **************************************************************************************************
!> \brief Check invariant-trace quality caps and mismatch trust for Direct candidates.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_trace_invariant_history_quality(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      REAL(KIND=dp) :: capped_quality, err, expected_trust, fallback_quality, good_quality, &
         immature_quality, immature_trust, mature_trust, neutral_trust

      fallback_quality = 0.40_dp
      immature_quality = dense_trace_invariant_history_quality( &
                         -0.40_dp, -0.02_dp, -0.55_dp, 0.15_dp, &
                         -0.03_dp, 0.01_dp, 3, fallback_quality)
      capped_quality = dense_trace_invariant_history_quality( &
                       -0.40_dp, -0.02_dp, -0.55_dp, 0.15_dp, &
                       -0.03_dp, 0.01_dp, 16, fallback_quality)
      good_quality = dense_trace_invariant_history_quality( &
                     -0.50_dp, -0.45_dp, -0.80_dp, 0.30_dp, &
                     -0.72_dp, 0.27_dp, 16, fallback_quality)
      immature_trust = dense_direct_trace_mismatch_trust(1.0E-4_dp, 2.0E-3_dp, &
                                                         capped_quality, 3)
      neutral_trust = dense_direct_trace_mismatch_trust(1.0E-4_dp, 2.0E-3_dp, &
                                                        1.0_dp, 16)
      mature_trust = dense_direct_trace_mismatch_trust(1.0E-4_dp, 2.0E-3_dp, &
                                                       capped_quality, 16)
      expected_trust = 1.0E-4_dp/(1.0E-4_dp + (1.0_dp - capped_quality)*2.0E-3_dp)

      err = ABS(immature_quality - fallback_quality)
      err = MAX(err, ABS(capped_quality - 0.10_dp))
      err = MAX(err, ABS(good_quality - fallback_quality))
      err = MAX(err, ABS(immature_trust - 1.0_dp))
      err = MAX(err, ABS(neutral_trust - 1.0_dp))
      err = MAX(err, ABS(mature_trust - expected_trust))

      IF (err > tol .OR. mature_trust >= 0.10_dp .OR. capped_quality >= fallback_quality) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,8(1X,ES12.5))') TRIM(label), &
            "failed with error/imm/cap/good/itrust/ntrust/mtrust/expected", &
            err, immature_quality, capped_quality, good_quality, immature_trust, neutral_trust, &
            mature_trust, expected_trust
      ELSE
         WRITE (*, '(A,1X,A,8(1X,ES12.5))') TRIM(label), &
            "OK with error/imm/cap/good/itrust/ntrust/mtrust/expected", &
            err, immature_quality, capped_quality, good_quality, immature_trust, neutral_trust, &
            mature_trust, expected_trust
      END IF

   END SUBROUTINE check_projector_trace_invariant_history_quality

! **************************************************************************************************
!> \brief Check that sparse selected Direct samples do not erase non-circular trace history.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_direct_selected_history_quality(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      REAL(KIND=dp), PARAMETER                           :: max_quality = 1.25_dp, &
                                                            min_quality = 0.10_dp

      REAL(KIND=dp) :: capped_floor_quality, early_bad_quality, early_bootstrap_quality, &
         early_capped_bootstrap_quality, early_floor_quality, early_good_quality, err, &
         mature_bad_quality, mature_bootstrap_quality, mature_floor_quality, &
         no_samples_bootstrap_quality, no_samples_floor_quality, no_samples_quality, &
         probe_early_quality, probe_mature_quality, recent_bootstrap_quality, &
         stale_bootstrap_quality

      no_samples_quality = dense_direct_selected_history_quality(1.0_dp, 10.0_dp, 1.0_dp, &
                                                                 0.05_dp, 3, 4, 8, &
                                                                 min_quality, max_quality)
      early_bad_quality = dense_direct_selected_history_quality(1.0_dp, 10.0_dp, 1.0_dp, &
                                                                0.05_dp, 4, 4, 8, &
                                                                min_quality, max_quality)
      mature_bad_quality = dense_direct_selected_history_quality(1.0_dp, 10.0_dp, 1.0_dp, &
                                                                 0.05_dp, 16, 4, 8, &
                                                                 min_quality, max_quality)
      early_good_quality = dense_direct_selected_history_quality(1.0_dp, 10.0_dp, 20.0_dp, &
                                                                 1.20_dp, 4, 4, 8, &
                                                                 min_quality, max_quality)
      probe_early_quality = dense_direct_selected_probe_quality(min_quality, 0.053_dp, &
                                                                4, 4, 8, max_quality)
      probe_mature_quality = dense_direct_selected_probe_quality(min_quality, 0.053_dp, &
                                                                 16, 4, 8, max_quality)
      no_samples_floor_quality = dense_direct_floor_selected_quality(max_quality, 0.053_dp, 3, 4, 8, &
                                                                     max_quality)
      early_floor_quality = dense_direct_floor_selected_quality(max_quality, 0.053_dp, 4, 4, 8, &
                                                                max_quality)
      mature_floor_quality = dense_direct_floor_selected_quality(max_quality, 0.053_dp, 16, 4, 8, &
                                                                 max_quality)
      capped_floor_quality = dense_direct_floor_selected_quality(0.439428811_dp, 0.053_dp, 3, 4, 8, &
                                                                 max_quality)
      no_samples_bootstrap_quality = dense_floor_bootstrap_history_quality(max_quality, &
                                                                           0.03_dp, 0.04_dp, 7, 8, &
                                                                           min_quality, max_quality)
      early_bootstrap_quality = dense_floor_bootstrap_history_quality(max_quality, &
                                                                      0.154524141_dp, 1.53795753_dp, 75, 8, &
                                                                      min_quality, max_quality)
      recent_bootstrap_quality = dense_floor_bootstrap_history_quality(max_quality, &
                                                                       0.171683247_dp, 0.440422594_dp, 79, 8, &
                                                                       min_quality, max_quality)
      early_capped_bootstrap_quality = dense_floor_bootstrap_history_quality(0.745910488_dp, &
                                                                             0.151939554_dp, 0.431667199_dp, 3, &
                                                                             1, min_quality, max_quality)
      mature_bootstrap_quality = dense_floor_bootstrap_history_quality(max_quality, &
                                                                       0.0309592052_dp, 0.0805830466_dp, 85, 8, &
                                                                       min_quality, max_quality)
      stale_bootstrap_quality = dense_floor_bootstrap_history_quality(max_quality, &
                                                                      0.0_dp, 0.0805830466_dp, 85, 8, &
                                                                      min_quality, max_quality)

      err = ABS(no_samples_quality - 1.0_dp)
      err = MAX(err, ABS(early_bad_quality - (1.0_dp - (1.0_dp/13.0_dp)*(1.0_dp - min_quality))))
      err = MAX(err, ABS(mature_bad_quality - min_quality))
      err = MAX(err, ABS(early_good_quality - (1.0_dp + (1.0_dp/13.0_dp)*(max_quality - 1.0_dp))))
      err = MAX(err, ABS(probe_early_quality - ((12.0_dp*min_quality + 0.053_dp)/13.0_dp)))
      err = MAX(err, ABS(probe_mature_quality - 0.053_dp))
      err = MAX(err, ABS(no_samples_floor_quality - max_quality))
      err = MAX(err, ABS(early_floor_quality - ((12.0_dp*max_quality + 0.053_dp)/13.0_dp)))
      err = MAX(err, ABS(mature_floor_quality - 0.053_dp))
      err = MAX(err, ABS(capped_floor_quality - 0.439428811_dp))
      err = MAX(err, ABS(no_samples_bootstrap_quality - max_quality))
      err = MAX(err, ABS(early_bootstrap_quality - SQRT(0.154524141_dp*max_quality)))
      err = MAX(err, ABS(recent_bootstrap_quality - SQRT(0.171683247_dp*0.440422594_dp)))
      err = MAX(err, ABS(early_capped_bootstrap_quality - &
                         SQRT(0.151939554_dp*0.431667199_dp)))
      err = MAX(err, ABS(mature_bootstrap_quality - &
                         MAX(min_quality, SQRT(0.0309592052_dp*min_quality))))
      err = MAX(err, ABS(stale_bootstrap_quality - max_quality))

      IF (err > tol .OR. early_bad_quality <= mature_bad_quality .OR. &
          early_bad_quality <= min_quality .OR. early_good_quality <= no_samples_quality .OR. &
          probe_early_quality >= min_quality .OR. probe_mature_quality >= min_quality .OR. &
          early_floor_quality <= mature_floor_quality .OR. mature_floor_quality >= min_quality .OR. &
          capped_floor_quality >= max_quality .OR. &
          early_bootstrap_quality >= max_quality .OR. recent_bootstrap_quality >= early_bootstrap_quality .OR. &
          early_capped_bootstrap_quality >= recent_bootstrap_quality .OR. &
          mature_bootstrap_quality >= recent_bootstrap_quality .OR. &
          stale_bootstrap_quality <= mature_bootstrap_quality) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,17(1X,ES12.5))') TRIM(label), &
            "failed err/no/early/mature/good/probe/mprobe/floor0/floor1/floorm/fcap/boot0/boot1/boot2/bootecap/bootm/stale", &
            err, no_samples_quality, &
            early_bad_quality, mature_bad_quality, early_good_quality, probe_early_quality, &
            probe_mature_quality, no_samples_floor_quality, early_floor_quality, mature_floor_quality, &
            capped_floor_quality, no_samples_bootstrap_quality, early_bootstrap_quality, &
            recent_bootstrap_quality, early_capped_bootstrap_quality, mature_bootstrap_quality, &
            stale_bootstrap_quality
      ELSE
         WRITE (*, '(A,1X,A,17(1X,ES12.5))') TRIM(label), &
            "OK err/no/early/mature/good/probe/mprobe/floor0/floor1/floorm/fcap/boot0/boot1/boot2/bootecap/bootm/stale", &
            err, no_samples_quality, &
            early_bad_quality, mature_bad_quality, early_good_quality, probe_early_quality, &
            probe_mature_quality, no_samples_floor_quality, early_floor_quality, mature_floor_quality, &
            capped_floor_quality, no_samples_bootstrap_quality, early_bootstrap_quality, &
            recent_bootstrap_quality, early_capped_bootstrap_quality, mature_bootstrap_quality, &
            stale_bootstrap_quality
      END IF

   END SUBROUTINE check_projector_direct_selected_history_quality

! **************************************************************************************************
!> \brief Check mature inactive-offdiag selected history can lose the residual-tail floor.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_inactive_offdiag_selected_tail_quality(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      REAL(KIND=dp), PARAMETER                           :: cal_quality = 0.25_dp, &
                                                            max_quality = 1.25_dp, &
                                                            tail_delta = 1.0E-4_dp

      REAL(KIND=dp)                                      :: al_tail_quality, err, expected_al_tail, &
                                                            good_tail_quality, immature_quality, &
                                                            outside_tail_quality, weak_tail_quality, &
                                                            zero_tail_quality

      immature_quality = dense_inactive_offdiag_selected_tail_quality( &
         0.75_dp, 0.75_dp, 1.0_dp, 0.75_dp, 0.02_dp, 0.02_dp, 0.25_dp, 4, 4, &
         tail_delta, tail_delta, cal_quality, max_quality)
      outside_tail_quality = dense_inactive_offdiag_selected_tail_quality( &
         0.75_dp, 0.75_dp, 1.0_dp, 0.75_dp, 0.02_dp, 0.02_dp, 0.25_dp, 16, 4, &
         2.0_dp*tail_delta, tail_delta, cal_quality, max_quality)
      good_tail_quality = dense_inactive_offdiag_selected_tail_quality( &
         0.75_dp, 0.75_dp, 1.0_dp, 0.75_dp, 0.50_dp, 0.50_dp, 0.25_dp, 16, 4, &
         tail_delta, tail_delta, cal_quality, max_quality)
      al_tail_quality = dense_inactive_offdiag_selected_tail_quality( &
         0.120821737_dp, 0.120821737_dp, 2.38717371E-4_dp, 2.84589460E-5_dp, &
         0.276863603_dp, 0.169491253_dp, 0.10_dp, 50, 4, tail_delta, tail_delta, &
         cal_quality, max_quality)
      weak_tail_quality = dense_inactive_offdiag_selected_tail_quality( &
         0.10_dp, 0.064148488_dp, 2.11544262E-4_dp, 2.03786650E-5_dp, &
         0.240690437_dp, 0.274936610_dp, 0.10_dp, 21, 4, tail_delta, tail_delta, &
         cal_quality, max_quality)
      zero_tail_quality = dense_inactive_offdiag_selected_tail_quality( &
         0.120821737_dp, 0.120821737_dp, 2.38717371E-4_dp, 2.84589460E-5_dp, &
         0.0_dp, 0.0_dp, 0.10_dp, 50, 4, tail_delta, tail_delta, cal_quality, &
         max_quality)

      expected_al_tail = 0.120821737_dp*SQRT((2.84589460E-5_dp/2.38717371E-4_dp)/ &
                                             cal_quality)* &
                         SQRT(0.169491253_dp/cal_quality)
      err = ABS(immature_quality - 0.75_dp)
      err = MAX(err, ABS(outside_tail_quality - 0.75_dp))
      err = MAX(err, ABS(good_tail_quality - 0.75_dp))
      err = MAX(err, ABS(al_tail_quality - expected_al_tail))
      err = MAX(err, ABS(weak_tail_quality))
      err = MAX(err, ABS(zero_tail_quality))

      IF (err > tol .OR. al_tail_quality >= 0.10_dp .OR. weak_tail_quality > tol .OR. &
          zero_tail_quality > tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,7(1X,ES12.5))') TRIM(label), &
            "failed err/immature/outside/good/al/weak/zero", err, immature_quality, &
            outside_tail_quality, good_tail_quality, al_tail_quality, weak_tail_quality, &
            zero_tail_quality
      ELSE
         WRITE (*, '(A,1X,A,7(1X,ES12.5))') TRIM(label), &
            "OK err/immature/outside/good/al/weak/zero", err, immature_quality, &
            outside_tail_quality, good_tail_quality, al_tail_quality, weak_tail_quality, &
            zero_tail_quality
      END IF

   END SUBROUTINE check_projector_inactive_offdiag_selected_tail_quality

! **************************************************************************************************
!> \brief Check residual-tail cap for offdiag projector history quality.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_offdiag_tail_residual_quality(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      REAL(KIND=dp), PARAMETER                           :: cal_quality = 0.25_dp, &
                                                            max_quality = 1.25_dp, &
                                                            tail_delta = 1.0E-4_dp

      REAL(KIND=dp)                                      :: err, expected_moderate, expected_near, &
                                                            expected_weak, good_quality, immature_quality, &
                                                            near_gate_quality, &
                                                            moderate_quality, outside_tail_quality, &
                                                            weak_quality, zero_quality

      immature_quality = dense_offdiag_tail_residual_quality( &
         1.25_dp, 0.0_dp, 0.0_dp, 4, tail_delta, tail_delta, cal_quality, max_quality)
      outside_tail_quality = dense_offdiag_tail_residual_quality( &
         1.25_dp, 0.0_dp, 0.0_dp, 16, 2.0_dp*tail_delta, tail_delta, cal_quality, max_quality)
      good_quality = dense_offdiag_tail_residual_quality( &
         1.25_dp, 0.50_dp, 0.50_dp, 16, tail_delta, tail_delta, cal_quality, max_quality)
      moderate_quality = dense_offdiag_tail_residual_quality( &
         1.25_dp, 0.04_dp, 0.05_dp, 16, tail_delta, tail_delta, cal_quality, max_quality)
      near_gate_quality = dense_offdiag_tail_residual_quality( &
         1.25_dp, 0.018_dp, 0.020_dp, 104, tail_delta, tail_delta, cal_quality, &
         max_quality)
      weak_quality = dense_offdiag_tail_residual_quality( &
         1.25_dp, 0.00536993433_dp, 0.00883931517_dp, 104, tail_delta, tail_delta, &
         cal_quality, max_quality)
      zero_quality = dense_offdiag_tail_residual_quality( &
         1.25_dp, 0.0_dp, 0.0_dp, 104, tail_delta, tail_delta, cal_quality, max_quality)

      expected_moderate = max_quality*(0.5_dp*(0.05_dp + 0.04_dp)/cal_quality)
      expected_near = max_quality*(0.5_dp*(0.020_dp + 0.018_dp)/cal_quality)
      err = ABS(near_gate_quality - expected_near)
      expected_weak = max_quality*(0.5_dp*(0.00883931517_dp + 0.00536993433_dp)/ &
                                   cal_quality)
      err = MAX(err, ABS(immature_quality - max_quality))
      err = MAX(err, ABS(outside_tail_quality - max_quality))
      err = MAX(err, ABS(good_quality - max_quality))
      err = MAX(err, ABS(moderate_quality - expected_moderate))
      err = MAX(err, ABS(weak_quality - expected_weak))
      err = MAX(err, ABS(zero_quality))

      IF (err > tol .OR. near_gate_quality >= 0.10_dp .OR. weak_quality >= 0.10_dp .OR. &
          zero_quality > tol .OR. &
          moderate_quality <= weak_quality) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,8(1X,ES12.5))') TRIM(label), &
            "failed err/immature/outside/good/moderate/near/weak/zero", err, immature_quality, &
            outside_tail_quality, good_quality, moderate_quality, near_gate_quality, &
            weak_quality, zero_quality
      ELSE
         WRITE (*, '(A,1X,A,8(1X,ES12.5))') TRIM(label), &
            "OK err/immature/outside/good/moderate/near/weak/zero", err, immature_quality, &
            outside_tail_quality, good_quality, moderate_quality, near_gate_quality, &
            weak_quality, zero_quality
      END IF

   END SUBROUTINE check_projector_offdiag_tail_residual_quality

! **************************************************************************************************
!> \brief Check that accepted Mermin history treats large response ratios as underprediction.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_response_underprediction_quality(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      REAL(KIND=dp), PARAMETER :: max_quality = 1.25_dp, min_quality = 0.10_dp, &
         underpredict_max_quality = 20.0_dp

      LOGICAL :: negative_normal_disabled, normal_enabled, overcap_disabled, recovered_enabled, &
         recovered_negative_disabled, underpredict_enabled, unstable_disabled, weak_disabled
      REAL(KIND=dp)                                      :: err, legacy_quality, underpredict_quality

      normal_enabled = dense_projector_response_candidate_enabled(8, 1.0_dp, 1.0_dp, &
                                                                  10.0_dp, 10.0_dp, 0.02_dp, &
                                                                  min_quality, max_quality, &
                                                                  underpredict_max_quality)
      negative_normal_disabled = dense_projector_response_candidate_enabled(8, 1.0_dp, 1.0_dp, &
                                                                            10.0_dp, 10.0_dp, 0.95_dp, &
                                                                            min_quality, max_quality, &
                                                                            underpredict_max_quality)
      underpredict_enabled = dense_projector_response_candidate_enabled(12, 14.0_dp, 16.0_dp, &
                                                                        1.0_dp, 18.0_dp, 0.02_dp, &
                                                                        min_quality, max_quality, &
                                                                        underpredict_max_quality)
      recovered_enabled = dense_projector_response_candidate_enabled(10, 0.16_dp, 1.30_dp, &
                                                                     30.0_dp, 6.0_dp, 0.02_dp, &
                                                                     min_quality, max_quality, &
                                                                     underpredict_max_quality)
      recovered_negative_disabled = dense_projector_response_candidate_enabled(10, 0.16_dp, 1.30_dp, &
                                                                               30.0_dp, 6.0_dp, 0.30_dp, &
                                                                               min_quality, max_quality, &
                                                                               underpredict_max_quality)
      unstable_disabled = dense_projector_response_candidate_enabled(12, 14.0_dp, 16.0_dp, &
                                                                     1.0_dp, 18.0_dp, 0.30_dp, &
                                                                     min_quality, max_quality, &
                                                                     underpredict_max_quality)
      overcap_disabled = dense_projector_response_candidate_enabled(12, 24.0_dp, 16.0_dp, &
                                                                    1.0_dp, 25.0_dp, 0.02_dp, &
                                                                    min_quality, max_quality, &
                                                                    underpredict_max_quality)
      weak_disabled = dense_projector_response_candidate_enabled(12, 0.03_dp, 0.16_dp, &
                                                                 1.0_dp, 18.0_dp, 0.02_dp, &
                                                                 min_quality, max_quality, &
                                                                 underpredict_max_quality)

      legacy_quality = dense_schur_history_quality(1.0_dp, 18.0_dp, 16.0_dp, 12, 8, 0.0_dp, &
                                                   min_quality, max_quality)
      underpredict_quality = dense_schur_history_quality(1.0_dp, 18.0_dp, 16.0_dp, 12, 8, 0.0_dp, &
                                                         min_quality, underpredict_max_quality)

      err = ABS(legacy_quality - max_quality)
      err = MAX(err, ABS(underpredict_quality - 18.0_dp))

      IF (err > tol .OR. .NOT. normal_enabled .OR. negative_normal_disabled .OR. &
          .NOT. underpredict_enabled .OR. .NOT. recovered_enabled .OR. recovered_negative_disabled .OR. &
          unstable_disabled .OR. overcap_disabled .OR. weak_disabled .OR. &
          underpredict_quality <= legacy_quality) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,3(1X,ES12.5),8(1X,L1))') TRIM(label), &
            "failed err/legacy/under/norm/nneg/under/rec/rneg/eneg/over/weak", err, legacy_quality, &
            underpredict_quality, normal_enabled, negative_normal_disabled, underpredict_enabled, &
            recovered_enabled, recovered_negative_disabled, unstable_disabled, overcap_disabled, weak_disabled
      ELSE
         WRITE (*, '(A,1X,A,3(1X,ES12.5),8(1X,L1))') TRIM(label), &
            "OK err/legacy/under/norm/nneg/under/rec/rneg/eneg/over/weak", err, legacy_quality, &
            underpredict_quality, normal_enabled, negative_normal_disabled, underpredict_enabled, &
            recovered_enabled, recovered_negative_disabled, unstable_disabled, overcap_disabled, weak_disabled
      END IF

   END SUBROUTINE check_projector_response_underprediction_quality

! **************************************************************************************************
!> \brief Check that negative Projector-response modes are only adverse when anti-descent.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_response_adverse_fraction(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      REAL(KIND=dp) :: aligned_fraction, anti_descent_fraction, clipped_fraction, err, &
         zero_correction_fraction, zero_negative_fraction

      aligned_fraction = dense_projector_response_adverse_fraction(4.0_dp, 5.0_dp, 1.0_dp)
      anti_descent_fraction = dense_projector_response_adverse_fraction(4.0_dp, 5.0_dp, -1.0_dp)
      clipped_fraction = dense_projector_response_adverse_fraction(7.0_dp, 5.0_dp, -1.0_dp)
      zero_correction_fraction = dense_projector_response_adverse_fraction(4.0_dp, 0.0_dp, -1.0_dp)
      zero_negative_fraction = dense_projector_response_adverse_fraction(0.0_dp, 5.0_dp, -1.0_dp)

      err = ABS(aligned_fraction)
      err = MAX(err, ABS(anti_descent_fraction - 0.8_dp))
      err = MAX(err, ABS(clipped_fraction - 1.0_dp))
      err = MAX(err, ABS(zero_correction_fraction))
      err = MAX(err, ABS(zero_negative_fraction))

      IF (err > tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,6(1X,ES12.5))') TRIM(label), &
            "failed err/aligned/anti/clipped/zcorr/zneg", err, aligned_fraction, &
            anti_descent_fraction, clipped_fraction, zero_correction_fraction, zero_negative_fraction
      ELSE
         WRITE (*, '(A,1X,A,6(1X,ES12.5))') TRIM(label), &
            "OK err/aligned/anti/clipped/zcorr/zneg", err, aligned_fraction, &
            anti_descent_fraction, clipped_fraction, zero_correction_fraction, zero_negative_fraction
      END IF

   END SUBROUTINE check_projector_response_adverse_fraction

! **************************************************************************************************
!> \brief Check raw/effective residual weighting for finite Projector-response tails.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_response_residual_weight(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      REAL(KIND=dp), PARAMETER :: floor_delta = 3.0E-6_dp, full_delta = 1.2E-5_dp

      REAL(KIND=dp) :: both_tiny, err, no_raw, raw_capped, raw_smaller, raw_tail
      REAL(KIND=dp) :: expected_both_tiny, expected_raw_capped, expected_raw_smaller, &
         expected_raw_tail

      no_raw = dense_projector_response_residual_weight(9.0E-6_dp, -1.0_dp)
      both_tiny = dense_projector_response_residual_weight(1.0E-7_dp, 1.0E-7_dp)
      raw_tail = dense_projector_response_residual_weight(1.0E-6_dp, 6.4E-5_dp)
      raw_capped = dense_projector_response_residual_weight(4.0E-6_dp, 1.0E-2_dp)
      raw_smaller = dense_projector_response_residual_weight(6.0E-6_dp, 2.0E-6_dp)

      expected_both_tiny = (1.0E-7_dp/full_delta)*SQRT(1.0E-7_dp/floor_delta)
      expected_raw_tail = SQRT(1.0E-6_dp*6.4E-5_dp)/full_delta
      expected_raw_capped = 1.0_dp
      expected_raw_smaller = 6.0E-6_dp/full_delta

      err = ABS(no_raw - 0.75_dp)
      err = MAX(err, ABS(both_tiny - expected_both_tiny))
      err = MAX(err, ABS(raw_tail - expected_raw_tail))
      err = MAX(err, ABS(raw_capped - expected_raw_capped))
      err = MAX(err, ABS(raw_smaller - expected_raw_smaller))

      IF (err > tol .OR. raw_tail <= no_raw*0.8_dp .OR. raw_smaller >= raw_tail) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,6(1X,ES12.5))') TRIM(label), &
            "failed err/no/btiny/rtail/rcap/rsmall", err, no_raw, both_tiny, raw_tail, &
            raw_capped, raw_smaller
      ELSE
         WRITE (*, '(A,1X,A,6(1X,ES12.5))') TRIM(label), &
            "OK err/no/btiny/rtail/rcap/rsmall", err, no_raw, both_tiny, raw_tail, &
            raw_capped, raw_smaller
      END IF

   END SUBROUTINE check_projector_response_residual_weight

! **************************************************************************************************
!> \brief Check that a trust-safe Direct tail candidate can collect selected history.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_direct_tail_probe_bootstrap(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      LOGICAL :: bootstrap_ready, enabled_veto, promoted_veto, selected_bad, selected_good, &
         selected_low_positive, selected_low_wait, selected_stalled_probe, selected_stalled_wait, &
         trust_veto
      REAL(KIND=dp)                                      :: err

      bootstrap_ready = dense_direct_tail_probe_ready(.TRUE., .FALSE., .FALSE., &
                                                      2.0E-5_dp, 4, 0, 0.0_dp, 0.0_dp, &
                                                      0.10_dp, 0.30_dp, 3.0E-3_dp, 2.0E-2_dp, 0.0_dp)
      promoted_veto = dense_direct_tail_probe_ready(.TRUE., .TRUE., .FALSE., &
                                                    2.0E-5_dp, 4, 0, 0.0_dp, 0.0_dp, &
                                                    0.10_dp, 0.30_dp, 3.0E-3_dp, 2.0E-2_dp, 0.0_dp)
      enabled_veto = dense_direct_tail_probe_ready(.TRUE., .FALSE., .TRUE., &
                                                   2.0E-5_dp, 4, 0, 0.0_dp, 0.0_dp, &
                                                   0.10_dp, 0.30_dp, 3.0E-3_dp, 2.0E-2_dp, 0.0_dp)
      trust_veto = dense_direct_tail_probe_ready(.TRUE., .FALSE., .FALSE., &
                                                 2.0E-5_dp, 4, 0, 0.0_dp, 0.0_dp, &
                                                 0.10_dp, 0.30_dp, 8.0E-1_dp, 2.0E-2_dp, 0.0_dp)
      selected_bad = dense_direct_tail_probe_ready(.TRUE., .FALSE., .FALSE., &
                                                   2.0E-5_dp, 4, 5, 0.0_dp, 0.0_dp, &
                                                   0.10_dp, 0.30_dp, 3.0E-3_dp, 2.0E-2_dp, 0.0_dp)
      selected_low_positive = dense_direct_tail_probe_ready(.TRUE., .FALSE., .FALSE., &
                                                            2.0E-5_dp, 16, 5, 0.056_dp, 0.053_dp, &
                                                            0.10_dp, 0.30_dp, 3.0E-3_dp, 2.0E-2_dp, 0.0_dp)
      selected_low_wait = dense_direct_tail_probe_ready(.TRUE., .FALSE., .FALSE., &
                                                        2.0E-5_dp, 17, 5, 0.056_dp, 0.053_dp, &
                                                        0.10_dp, 0.30_dp, 3.0E-3_dp, 2.0E-2_dp, 0.0_dp)
      selected_stalled_probe = dense_direct_tail_probe_ready(.TRUE., .FALSE., .FALSE., &
                                                             2.0E-5_dp, 24, 5, 0.056_dp, 0.030_dp, &
                                                             0.10_dp, 0.10_dp, 3.0E-3_dp, 2.0E-2_dp, 0.0_dp)
      selected_stalled_wait = dense_direct_tail_probe_ready(.TRUE., .FALSE., .FALSE., &
                                                            2.0E-5_dp, 25, 5, 0.056_dp, 0.030_dp, &
                                                            0.10_dp, 0.10_dp, 3.0E-3_dp, 2.0E-2_dp, 0.0_dp)
      selected_good = dense_direct_tail_probe_ready(.TRUE., .FALSE., .FALSE., &
                                                    2.0E-5_dp, 16, 5, 0.14_dp, 0.08_dp, &
                                                    0.10_dp, 0.30_dp, 3.0E-3_dp, 2.0E-2_dp, 0.0_dp)
      err = 0.0_dp
      IF (.NOT. bootstrap_ready .OR. promoted_veto .OR. enabled_veto .OR. trust_veto .OR. &
          selected_bad .OR. selected_low_wait .OR. .NOT. selected_low_positive .OR. &
          .NOT. selected_stalled_probe .OR. selected_stalled_wait .OR. .NOT. selected_good) err = 1.0_dp

      IF (err > tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,10(1X,L1))') TRIM(label), &
            "failed with boot/prom/en/trust/bad/wait/low/stall/wait/good", bootstrap_ready, &
            promoted_veto, enabled_veto, trust_veto, selected_bad, selected_low_wait, selected_low_positive, &
            selected_stalled_probe, selected_stalled_wait, selected_good
      ELSE
         WRITE (*, '(A,1X,A,10(1X,L1))') TRIM(label), &
            "OK with boot/prom/en/trust/bad/wait/low/stall/wait/good", bootstrap_ready, &
            promoted_veto, enabled_veto, trust_veto, selected_bad, selected_low_wait, selected_low_positive, &
            selected_stalled_probe, selected_stalled_wait, selected_good
      END IF

   END SUBROUTINE check_projector_direct_tail_probe_bootstrap

! **************************************************************************************************
!> \brief Check that full-density Direct candidates can collect close-call tail history.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_density_full_direct_tail_probe(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      LOGICAL :: low_drop, low_quality, low_samples, ready, selected_bad, selected_floor, &
         selected_ok, selected_recent_bad, stale_interval_veto, stale_refresh, tail_veto, weak_drop
      REAL(KIND=dp)                                      :: err

      ready = dense_density_full_direct_tail_probe(5.0E-5_dp, 8.5E-8_dp, 1.0E-7_dp, &
                                                   0.20_dp, 8, 0.0_dp, 0.0_dp, 0, 0.40_dp)
      tail_veto = dense_density_full_direct_tail_probe(2.0E-4_dp, 8.5E-8_dp, 1.0E-7_dp, &
                                                       0.20_dp, 8, 0.0_dp, 0.0_dp, 0, 0.40_dp)
      low_samples = dense_density_full_direct_tail_probe(5.0E-5_dp, 8.5E-8_dp, 1.0E-7_dp, &
                                                         0.20_dp, 7, 0.0_dp, 0.0_dp, 0, 0.40_dp)
      low_quality = dense_density_full_direct_tail_probe(5.0E-5_dp, 8.5E-8_dp, 1.0E-7_dp, &
                                                         0.09_dp, 8, 0.0_dp, 0.0_dp, 0, 0.40_dp)
      low_drop = dense_density_full_direct_tail_probe(5.0E-5_dp, 1.0E-13_dp, 1.0E-7_dp, &
                                                      0.20_dp, 8, 0.0_dp, 0.0_dp, 0, 0.40_dp)
      weak_drop = dense_density_full_direct_tail_probe(5.0E-5_dp, 7.9E-8_dp, 1.0E-7_dp, &
                                                       0.20_dp, 8, 0.0_dp, 0.0_dp, 0, 0.40_dp)
      selected_ok = dense_density_full_direct_tail_probe(5.0E-5_dp, 8.5E-8_dp, 1.0E-7_dp, &
                                                         0.20_dp, 8, 0.40_dp, 0.40_dp, 4, 0.40_dp)
      selected_bad = dense_density_full_direct_tail_probe(5.0E-5_dp, 8.5E-8_dp, 1.0E-7_dp, &
                                                          0.20_dp, 8, 0.40_dp, 0.20_dp, 4, 0.70_dp)
      selected_recent_bad = dense_density_full_direct_tail_probe(5.0E-5_dp, 8.5E-8_dp, 1.0E-7_dp, &
                                                                 0.20_dp, 8, 0.05_dp, 0.40_dp, 4, 0.40_dp)
      selected_floor = dense_density_full_direct_tail_probe(5.0E-5_dp, 8.5E-8_dp, 1.0E-7_dp, &
                                                            0.20_dp, 8, 0.21_dp, 0.21_dp, 4, 0.0_dp)
      stale_refresh = dense_density_full_direct_tail_probe(5.0E-5_dp, 2.5E-8_dp, 1.0E-7_dp, &
                                                           0.25_dp, 8, 0.05_dp, 0.10_dp, 4, 0.25_dp)
      stale_interval_veto = dense_density_full_direct_tail_probe(5.0E-5_dp, 2.5E-8_dp, 1.0E-7_dp, &
                                                                 0.25_dp, 9, 0.05_dp, 0.10_dp, 4, 0.25_dp)

      err = 0.0_dp
      IF (.NOT. ready .OR. tail_veto .OR. low_samples .OR. low_quality .OR. low_drop .OR. &
          weak_drop .OR. .NOT. selected_ok .OR. selected_bad .OR. selected_recent_bad .OR. &
          .NOT. selected_floor .OR. .NOT. stale_refresh .OR. stale_interval_veto) err = 1.0_dp

      IF (err > tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,12(1X,L1))') TRIM(label), &
            "failed with ready/tail/samples/quality/drop/weak/sel/bad/recent/floor/stale/int", ready, tail_veto, &
            low_samples, low_quality, low_drop, weak_drop, selected_ok, selected_bad, selected_recent_bad, &
            selected_floor, stale_refresh, stale_interval_veto
      ELSE
         WRITE (*, '(A,1X,A,12(1X,L1))') TRIM(label), &
            "OK with ready/tail/samples/quality/drop/weak/sel/bad/recent/floor/stale/int", ready, tail_veto, &
            low_samples, low_quality, low_drop, weak_drop, selected_ok, selected_bad, selected_recent_bad, &
            selected_floor, stale_refresh, stale_interval_veto
      END IF

   END SUBROUTINE check_projector_density_full_direct_tail_probe

! **************************************************************************************************
!> \brief Check that fixed-N density-descent candidates can collect residual-tail history.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_density_descent_tail_probe(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      LOGICAL :: low_drop, low_quality, low_recovery, ready, selected_bad, selected_ok, &
         selected_recent_bad, tail_veto, trust_veto, weak_drop
      REAL(KIND=dp)                                      :: err

      ready = dense_density_descent_tail_probe(5.0E-5_dp, 8.5E-8_dp, 1.0E-7_dp, &
                                               0.06_dp, 0.80_dp, 1.0_dp, 1.0_dp, &
                                               0.0_dp, 0.0_dp, 0, 0.40_dp)
      tail_veto = dense_density_descent_tail_probe(2.0E-4_dp, 8.5E-8_dp, 1.0E-7_dp, &
                                                   0.06_dp, 0.80_dp, 1.0_dp, 1.0_dp, &
                                                   0.0_dp, 0.0_dp, 0, 0.40_dp)
      low_quality = dense_density_descent_tail_probe(5.0E-5_dp, 8.5E-8_dp, 1.0E-7_dp, &
                                                     0.04_dp, 0.80_dp, 1.0_dp, 1.0_dp, &
                                                     0.0_dp, 0.0_dp, 0, 0.40_dp)
      low_recovery = dense_density_descent_tail_probe(5.0E-5_dp, 8.5E-8_dp, 1.0E-7_dp, &
                                                      0.06_dp, 0.49_dp, 1.0_dp, 1.0_dp, &
                                                      0.0_dp, 0.0_dp, 0, 0.40_dp)
      low_drop = dense_density_descent_tail_probe(5.0E-5_dp, 1.0E-13_dp, 1.0E-7_dp, &
                                                  0.06_dp, 0.80_dp, 1.0_dp, 1.0_dp, &
                                                  0.0_dp, 0.0_dp, 0, 0.40_dp)
      weak_drop = dense_density_descent_tail_probe(5.0E-5_dp, 7.9E-8_dp, 1.0E-7_dp, &
                                                   0.06_dp, 0.80_dp, 1.0_dp, 1.0_dp, &
                                                   0.0_dp, 0.0_dp, 0, 0.40_dp)
      trust_veto = dense_density_descent_tail_probe(5.0E-5_dp, 8.5E-8_dp, 1.0E-7_dp, &
                                                    0.06_dp, 0.80_dp, 1.01_dp, 1.0_dp, &
                                                    0.0_dp, 0.0_dp, 0, 0.40_dp)
      selected_ok = dense_density_descent_tail_probe(5.0E-5_dp, 8.5E-8_dp, 1.0E-7_dp, &
                                                     0.06_dp, 0.80_dp, 1.0_dp, 1.0_dp, &
                                                     0.40_dp, 0.40_dp, 4, 0.40_dp)
      selected_bad = dense_density_descent_tail_probe(5.0E-5_dp, 8.5E-8_dp, 1.0E-7_dp, &
                                                      0.06_dp, 0.80_dp, 1.0_dp, 1.0_dp, &
                                                      0.40_dp, 0.20_dp, 4, 0.40_dp)
      selected_recent_bad = dense_density_descent_tail_probe(5.0E-5_dp, 8.5E-8_dp, 1.0E-7_dp, &
                                                             0.06_dp, 0.80_dp, 1.0_dp, 1.0_dp, &
                                                             0.05_dp, 0.40_dp, 4, 0.40_dp)
      err = 0.0_dp
      IF (.NOT. ready .OR. tail_veto .OR. low_quality .OR. low_recovery .OR. low_drop .OR. &
          weak_drop .OR. trust_veto .OR. .NOT. selected_ok .OR. selected_bad .OR. &
          selected_recent_bad) err = 1.0_dp

      IF (err > tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,10(1X,L1))') TRIM(label), &
            "failed with ready/tail/quality/recovery/drop/weak/trust/sel/bad/recent", ready, tail_veto, &
            low_quality, low_recovery, low_drop, weak_drop, trust_veto, selected_ok, selected_bad, &
            selected_recent_bad
      ELSE
         WRITE (*, '(A,1X,A,10(1X,L1))') TRIM(label), &
            "OK with ready/tail/quality/recovery/drop/weak/trust/sel/bad/recent", ready, tail_veto, &
            low_quality, low_recovery, low_drop, weak_drop, trust_veto, selected_ok, selected_bad, &
            selected_recent_bad
      END IF

   END SUBROUTINE check_projector_density_descent_tail_probe

! **************************************************************************************************
!> \brief Check that residual-tail Direct shadow drops are visible only to trust-safe probes.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_density_full_direct_probe_drop(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      REAL(KIND=dp)                                      :: abs_trust_keep, capped_keep, err, &
                                                            low_shadow, outside_tail, probe_drop, &
                                                            weighted_trust_keep

      probe_drop = dense_density_full_direct_probe_drop(5.0E-5_dp, 0.0_dp, 8.0E-8_dp, &
                                                        1.0E-2_dp, 1.0E-2_dp)
      outside_tail = dense_density_full_direct_probe_drop(2.0E-4_dp, 0.0_dp, 8.0E-8_dp, &
                                                          1.0E-2_dp, 1.0E-2_dp)
      capped_keep = dense_density_full_direct_probe_drop(5.0E-5_dp, 3.0E-8_dp, 2.0E-8_dp, &
                                                         1.0E-2_dp, 1.0E-2_dp)
      low_shadow = dense_density_full_direct_probe_drop(5.0E-5_dp, 0.0_dp, 1.0E-13_dp, &
                                                        1.0E-2_dp, 1.0E-2_dp)
      weighted_trust_keep = dense_density_full_direct_probe_drop(5.0E-5_dp, 2.0E-8_dp, &
                                                                 8.0E-8_dp, 8.0E-1_dp, 1.0E-2_dp)
      abs_trust_keep = dense_density_full_direct_probe_drop(5.0E-5_dp, 2.0E-8_dp, 8.0E-8_dp, &
                                                            1.0E-2_dp, 8.0E-1_dp)

      err = MAX(ABS(probe_drop - 8.0E-8_dp), ABS(outside_tail), ABS(capped_keep - 3.0E-8_dp), &
                ABS(low_shadow), ABS(weighted_trust_keep - 2.0E-8_dp), &
                ABS(abs_trust_keep - 2.0E-8_dp))

      IF (err > tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,6(1X,ES12.5))') TRIM(label), &
            "failed with probe/out/cap/low/wtrust/atrust", probe_drop, outside_tail, capped_keep, &
            low_shadow, weighted_trust_keep, abs_trust_keep
      ELSE
         WRITE (*, '(A,1X,A,6(1X,ES12.5))') TRIM(label), &
            "OK with probe/out/cap/low/wtrust/atrust", probe_drop, outside_tail, capped_keep, &
            low_shadow, weighted_trust_keep, abs_trust_keep
      END IF

   END SUBROUTINE check_projector_density_full_direct_probe_drop

! **************************************************************************************************
!> \brief Check that selected Direct-density history only gives a tiny non-promoting tail mix.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_density_selected_tail_mix(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      REAL(KIND=dp) :: early_expected, early_weight, err, good_mix_quality, good_strict_quality, &
         mature_weight, strict_zero_quality, tail_veto, weak_quality, young_history, zero_mix_quality

      strict_zero_quality = dense_selected_tail_residual_quality(0.0_dp, 0.1188258236_dp)
      zero_mix_quality = dense_selected_tail_mix_residual_quality(0.0_dp, 0.1188258236_dp)
      good_strict_quality = dense_selected_tail_residual_quality(0.25_dp, 0.25_dp)
      good_mix_quality = dense_selected_tail_mix_residual_quality(0.25_dp, 0.25_dp)
      tail_veto = dense_direct_density_selected_mix_weight(2.0E-4_dp, 0.0267028466_dp, &
                                                           zero_mix_quality, 4)
      young_history = dense_direct_density_selected_mix_weight(2.0E-5_dp, 0.0267028466_dp, &
                                                               zero_mix_quality, 3)
      weak_quality = dense_direct_density_selected_mix_weight(2.0E-5_dp, 0.024_dp, &
                                                              zero_mix_quality, 4)
      early_weight = dense_direct_density_selected_mix_weight(2.0E-5_dp, 0.0267028466_dp, &
                                                              zero_mix_quality, 4)
      early_expected = 0.25_dp*0.95_dp*(1.0_dp/13.0_dp)* &
                       SQRT(0.0267028466_dp)*SQRT(zero_mix_quality)
      mature_weight = dense_direct_density_selected_mix_weight(2.0E-5_dp, 1.0_dp, 1.0_dp, 20)

      err = ABS(strict_zero_quality)
      err = MAX(err, ABS(zero_mix_quality - 0.0594129118_dp))
      err = MAX(err, ABS(good_strict_quality - 0.25_dp))
      err = MAX(err, ABS(good_mix_quality - 0.5_dp))
      err = MAX(err, ABS(tail_veto))
      err = MAX(err, ABS(young_history))
      err = MAX(err, ABS(weak_quality))
      err = MAX(err, ABS(early_weight - early_expected))
      err = MAX(err, ABS(mature_weight - 0.2375_dp))

      IF (err > tol .OR. early_weight <= 0.0_dp .OR. early_weight >= 1.0E-3_dp) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,10(1X,ES12.5))') TRIM(label), &
            "failed err/strict0/mix0/strictg/mixg/tail/young/weak/early/expected/mature", &
            err, strict_zero_quality, zero_mix_quality, good_strict_quality, good_mix_quality, &
            tail_veto, young_history, weak_quality, early_weight, early_expected, mature_weight
      ELSE
         WRITE (*, '(A,1X,A,10(1X,ES12.5))') TRIM(label), &
            "OK err/strict0/mix0/strictg/mixg/tail/young/weak/early/expected/mature", &
            err, strict_zero_quality, zero_mix_quality, good_strict_quality, good_mix_quality, &
            tail_veto, young_history, weak_quality, early_weight, early_expected, mature_weight
      END IF

   END SUBROUTINE check_projector_density_selected_tail_mix

! **************************************************************************************************
!> \brief Check the residual-tail dense metric-response bridge gate.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_dense_metric_tail_bridge(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      LOGICAL                                            :: calibrated_rescue, low_ema, low_mermin, &
                                                            low_quality, low_residual, &
                                                            low_samples, raw_stale, ready, &
                                                            tail_veto
      REAL(KIND=dp)                                      :: err

      ready = dense_metric_tail_bridge_ready(5.0E-5_dp, 8.1E-2_dp, 9.7E-2_dp, 22, &
                                             2.5E-1_dp, 9.9E-1_dp, 40)
      tail_veto = dense_metric_tail_bridge_ready(2.0E-4_dp, 8.1E-2_dp, 9.7E-2_dp, 22, &
                                                 2.5E-1_dp, 9.9E-1_dp, 40)
      low_samples = dense_metric_tail_bridge_ready(5.0E-5_dp, 8.1E-2_dp, 9.7E-2_dp, 3, &
                                                   2.5E-1_dp, 9.9E-1_dp, 40)
      low_quality = dense_metric_tail_bridge_ready(5.0E-5_dp, 4.9E-2_dp, 9.7E-2_dp, 22, &
                                                   2.5E-1_dp, 9.9E-1_dp, 40)
      low_ema = dense_metric_tail_bridge_ready(5.0E-5_dp, 8.1E-2_dp, 7.4E-2_dp, 22, &
                                               2.5E-1_dp, 9.9E-1_dp, 40)
      low_residual = dense_metric_tail_bridge_ready(5.0E-5_dp, 8.1E-2_dp, 9.7E-2_dp, 22, &
                                                    2.4E-1_dp, 9.9E-1_dp, 40)
      low_mermin = dense_metric_tail_bridge_ready(5.0E-5_dp, 8.1E-2_dp, 9.7E-2_dp, 22, &
                                                  2.5E-1_dp, 7.9E-1_dp, 40)
      raw_stale = dense_metric_tail_bridge_ready(5.0E-5_dp, 3.0E-2_dp, 4.0E-2_dp, 22, &
                                                 9.0E-1_dp, 9.9E-1_dp, 40)
      calibrated_rescue = dense_metric_tail_bridge_ready(5.0E-5_dp, 8.0E-1_dp, 6.0E-1_dp, 22, &
                                                         9.0E-1_dp, 9.9E-1_dp, 40)

      err = 0.0_dp
      IF (.NOT. ready .OR. tail_veto .OR. low_samples .OR. low_quality .OR. low_ema .OR. &
          low_residual .OR. low_mermin .OR. raw_stale .OR. .NOT. calibrated_rescue) err = 1.0_dp

      IF (err > tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,9(1X,L1))') TRIM(label), &
            "failed with ready/tail/samples/quality/ema/residual/mermin/raw/cal", ready, tail_veto, &
            low_samples, low_quality, low_ema, low_residual, low_mermin, raw_stale, calibrated_rescue
      ELSE
         WRITE (*, '(A,1X,A,9(1X,L1))') TRIM(label), &
            "OK with ready/tail/samples/quality/ema/residual/mermin/raw/cal", ready, tail_veto, &
            low_samples, low_quality, low_ema, low_residual, low_mermin, raw_stale, calibrated_rescue
      END IF

   END SUBROUTINE check_projector_dense_metric_tail_bridge

! **************************************************************************************************
!> \brief Check the residual-aware gauge rotation floor.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_residual_gauge_floor(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      REAL(KIND=dp)                                      :: err, gap_pair, gauge_pair, good_history, &
                                                            mixed_pair

      good_history = residual_gauge_rotation_floor(0.5_dp, 1.5_dp, 0.0_dp, 1.0_dp)
      gap_pair = residual_gauge_rotation_floor(0.5_dp, 1.5_dp, 1.0_dp, 0.25_dp)
      gauge_pair = residual_gauge_rotation_floor(0.5_dp, 1.5_dp, 0.0_dp, 0.25_dp)
      mixed_pair = residual_gauge_rotation_floor(0.5_dp, 1.5_dp, 0.4_dp, 0.5_dp)

      err = MAX(ABS(good_history - 0.5_dp), ABS(gap_pair - 0.5_dp))
      err = MAX(err, ABS(gauge_pair - 1.25_dp))
      err = MAX(err, ABS(mixed_pair - 0.8_dp))

      IF (err > tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,5(1X,ES12.5))') TRIM(label), &
            "failed with error/good/gap/gauge/mixed", err, good_history, gap_pair, gauge_pair, mixed_pair
      ELSE
         WRITE (*, '(A,1X,A,5(1X,ES12.5))') TRIM(label), &
            "OK with error/good/gap/gauge/mixed", err, good_history, gap_pair, gauge_pair, mixed_pair
      END IF

   END SUBROUTINE check_projector_residual_gauge_floor

! **************************************************************************************************
!> \brief Check that raw trace-gauge motion promotes the active Direct Hessian blend.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_trace_active_mix_weight(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      REAL(KIND=dp), PARAMETER :: active_trace = -4.0E-3_dp, history_weight = 0.25_dp, &
         poor_quality = 0.05_dp, quality = 1.20_dp, raw_trace = -1.0_dp

      REAL(KIND=dp) :: err, expected_gauge, gauge_mix, mature_mix, no_raw_mix, poor_mix, &
         saturated_mix, selected_poor_mix, shadow_only_mix, young_mix

      young_mix = dense_trace_active_mix_weight(raw_trace, active_trace, history_weight, &
                                                quality, 2, quality, 16)
      poor_mix = dense_trace_active_mix_weight(raw_trace, active_trace, history_weight, &
                                               poor_quality, 16, quality, 16)
      selected_poor_mix = dense_trace_active_mix_weight(raw_trace, active_trace, history_weight, &
                                                        quality, 16, poor_quality, 16)
      shadow_only_mix = dense_trace_active_mix_weight(raw_trace, active_trace, history_weight, &
                                                      quality, 16, 0.0_dp, 0)
      gauge_mix = dense_trace_active_mix_weight(raw_trace, active_trace, history_weight, &
                                                quality, 16, quality, 16)
      mature_mix = dense_trace_active_mix_weight(raw_trace, raw_trace, history_weight, &
                                                 quality, 16, quality, 16)
      saturated_mix = dense_trace_active_mix_weight(raw_trace, 0.0_dp, 0.90_dp, &
                                                    quality, 16, quality, 16)
      no_raw_mix = dense_trace_active_mix_weight(0.0_dp, active_trace, history_weight, &
                                                 quality, 16, quality, 16)

      expected_gauge = 0.95_dp*(1.0_dp - ABS(active_trace/raw_trace))
      err = ABS(young_mix - history_weight)
      err = MAX(err, ABS(poor_mix - history_weight))
      err = MAX(err, ABS(selected_poor_mix - history_weight))
      err = MAX(err, ABS(shadow_only_mix - history_weight))
      err = MAX(err, ABS(gauge_mix - expected_gauge))
      err = MAX(err, ABS(mature_mix - history_weight))
      err = MAX(err, ABS(saturated_mix - 0.95_dp))
      err = MAX(err, ABS(no_raw_mix - history_weight))

      IF (err > tol .OR. gauge_mix <= history_weight .OR. saturated_mix > 0.95_dp) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,7(1X,ES12.5))') TRIM(label), &
            "failed with error/young/poor/selpoor/gauge/mature/saturated", err, young_mix, poor_mix, &
            selected_poor_mix, gauge_mix, mature_mix, saturated_mix
      ELSE
         WRITE (*, '(A,1X,A,7(1X,ES12.5))') TRIM(label), &
            "OK with error/young/poor/selpoor/gauge/mature/saturated", err, young_mix, poor_mix, &
            selected_poor_mix, gauge_mix, mature_mix, saturated_mix
      END IF

   END SUBROUTINE check_projector_trace_active_mix_weight

! **************************************************************************************************
!> \brief Check non-promoting active-basis Direct trace blending.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_trace_active_basis_weight(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      REAL(KIND=dp), PARAMETER                           :: active_trace = -5.5E-2_dp, &
                                                            history_quality = 4.7E-2_dp, &
                                                            raw_trace = -1.0_dp

      REAL(KIND=dp)                                      :: al_like_weight, err, matched_weight, &
                                                            no_raw_weight, poor_weight, &
                                                            saturated_weight, young_weight

      young_weight = dense_trace_active_basis_weight(raw_trace, active_trace, history_quality, 4)
      poor_weight = dense_trace_active_basis_weight(raw_trace, active_trace, 2.0E-2_dp, 16)
      no_raw_weight = dense_trace_active_basis_weight(0.0_dp, active_trace, history_quality, 16)
      matched_weight = dense_trace_active_basis_weight(raw_trace, raw_trace, history_quality, 16)
      al_like_weight = dense_trace_active_basis_weight(raw_trace, active_trace, history_quality, 16)
      saturated_weight = dense_trace_active_basis_weight(raw_trace, 0.0_dp, 1.0_dp, 16)

      err = ABS(young_weight)
      err = MAX(err, ABS(poor_weight))
      err = MAX(err, ABS(no_raw_weight))
      err = MAX(err, ABS(matched_weight))
      err = MAX(err, ABS(saturated_weight - 0.25_dp))

      IF (err > tol .OR. al_like_weight <= 0.18_dp .OR. al_like_weight >= 0.22_dp .OR. &
          saturated_weight > 0.25_dp) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,6(1X,ES12.5))') TRIM(label), &
            "failed with error/young/poor/matched/al/sat", err, young_weight, poor_weight, &
            matched_weight, al_like_weight, saturated_weight
      ELSE
         WRITE (*, '(A,1X,A,6(1X,ES12.5))') TRIM(label), &
            "OK with error/young/poor/matched/al/sat", err, young_weight, poor_weight, &
            matched_weight, al_like_weight, saturated_weight
      END IF

   END SUBROUTINE check_projector_trace_active_basis_weight

! **************************************************************************************************
!> \brief Check fixed-N trace-balance damping of the Direct trace diagonal component.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_trace_diag_component_scale(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      REAL(KIND=dp) :: balanced_scale, err, expected_gauge, gauge_scale, no_samples_scale, &
         poor_diag_scale, zero_free_poor_scale, zero_free_scale

      no_samples_scale = dense_trace_diag_component_scale(-1.0_dp, 0.75_dp, 1.0_dp, 1.0_dp, 2)
      zero_free_scale = dense_trace_diag_component_scale(-1.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 16)
      zero_free_poor_scale = dense_trace_diag_component_scale(-1.0_dp, 0.0_dp, 0.04_dp, 1.0_dp, 16)
      balanced_scale = dense_trace_diag_component_scale(0.75_dp, 0.25_dp, 1.0_dp, 0.1_dp, 16)
      gauge_scale = dense_trace_diag_component_scale(-1.0_dp, 0.75_dp, 1.0_dp, 0.8_dp, 16)
      poor_diag_scale = dense_trace_diag_component_scale(-1.0_dp, 0.75_dp, 0.08_dp, 0.8_dp, 16)

      expected_gauge = SQRT(0.25_dp/1.75_dp)
      expected_gauge = SQRT(MAX(0.10_dp, expected_gauge*expected_gauge))
      err = ABS(no_samples_scale - 1.0_dp)
      err = MAX(err, ABS(zero_free_scale - 1.0_dp))
      err = MAX(err, ABS(zero_free_poor_scale - SQRT(0.10_dp)))
      err = MAX(err, ABS(balanced_scale - 1.0_dp))
      err = MAX(err, ABS(gauge_scale - expected_gauge))
      err = MAX(err, ABS(poor_diag_scale - SQRT(0.10_dp)))

      IF (err > tol .OR. gauge_scale >= balanced_scale .OR. &
          poor_diag_scale > gauge_scale .OR. zero_free_poor_scale >= zero_free_scale) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,7(1X,ES12.5))') TRIM(label), &
            "failed with error/no/free/zpoor/bal/gauge/poor", err, no_samples_scale, zero_free_scale, &
            zero_free_poor_scale, balanced_scale, gauge_scale, poor_diag_scale
      ELSE
         WRITE (*, '(A,1X,A,7(1X,ES12.5))') TRIM(label), &
            "OK with error/no/free/zpoor/bal/gauge/poor", err, no_samples_scale, zero_free_scale, &
            zero_free_poor_scale, balanced_scale, gauge_scale, poor_diag_scale
      END IF

   END SUBROUTINE check_projector_trace_diag_component_scale

! **************************************************************************************************
!> \brief Check quality-gated scaling of local Direct trace compensation.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_trace_local_compensation(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: ntrace = 4
      REAL(KIND=dp), PARAMETER :: base_diag = -0.20_dp, coupled_weight = 0.50_dp, &
         local_comp = 0.08_dp, mixed_weight = 0.75_dp, quality_floor = 0.25_dp

      INTEGER                                            :: itrace
      REAL(KIND=dp) :: err, expected_active_shift, global_active_comp, global_inactive_shift, &
         high_old, high_scale, high_sep, high_with_comp, local_active_shift, local_inactive_shift, &
         low_scale, low_with_comp, no_comp, trace_scale
      REAL(KIND=dp), DIMENSION(ntrace)                   :: active_delta, local_active_comp, &
                                                            reference_delta

      low_scale = SQRT(quality_floor)
      high_scale = 0.90_dp
      high_sep = (high_scale - low_scale)/(1.0_dp - low_scale)
      no_comp = dense_trace_scaled_diag_with_local_comp(base_diag, 0.0_dp, &
                                                        high_scale, quality_floor)
      low_with_comp = dense_trace_scaled_diag_with_local_comp(base_diag, local_comp, &
                                                              low_scale, quality_floor)
      high_with_comp = dense_trace_scaled_diag_with_local_comp(base_diag, local_comp, &
                                                               high_scale, quality_floor)
      high_old = high_scale*(base_diag + local_comp)
      err = ABS(low_with_comp - low_scale*(base_diag + local_comp))
      err = MAX(err, ABS(high_with_comp - &
                         (base_diag + local_comp + &
                          (high_scale - 1.0_dp)*(base_diag + (1.0_dp - high_sep)*local_comp))))
      err = MAX(err, ABS(no_comp - high_scale*base_diag))

      trace_scale = 0.80_dp
      active_delta(:) = [-0.20_dp, -0.01_dp, 0.01_dp, 0.02_dp]
      reference_delta(:) = [0.10_dp, 0.03_dp, 0.01_dp, 0.03_dp]
      DO itrace = 1, ntrace
         local_active_comp(itrace) = trace_scale* &
                                     dense_direct_trace_compensation(coupled_weight, &
                                                                     active_delta(itrace), &
                                                                     reference_delta(itrace))
      END DO
      global_active_comp = trace_scale* &
                           dense_direct_trace_compensation(coupled_weight, &
                                                           SUM(active_delta(:))/REAL(ntrace, KIND=dp), &
                                                           SUM(reference_delta(:))/REAL(ntrace, KIND=dp))
      local_active_shift = mixed_weight*local_active_comp(1)
      local_inactive_shift = mixed_weight*local_active_comp(ntrace)
      global_inactive_shift = mixed_weight*global_active_comp
      expected_active_shift = mixed_weight*trace_scale* &
                              dense_direct_trace_compensation(coupled_weight, active_delta(1), &
                                                              reference_delta(1))
      err = MAX(err, ABS(local_active_shift - expected_active_shift))
      err = MAX(err, ABS(local_inactive_shift))

      IF (err > tol .OR. high_with_comp <= high_old .OR. &
          local_active_shift <= global_inactive_shift .OR. &
          global_inactive_shift <= tol .OR. local_inactive_shift > tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,9(1X,ES12.5))') TRIM(label), &
            "failed with error/no/low/high/old/sep/lact/linact/ginact", err, no_comp, &
            low_with_comp, high_with_comp, high_old, high_sep, local_active_shift, &
            local_inactive_shift, global_inactive_shift
      ELSE
         WRITE (*, '(A,1X,A,9(1X,ES12.5))') TRIM(label), &
            "OK with error/no/low/high/old/sep/lact/linact/ginact", err, no_comp, &
            low_with_comp, high_with_comp, high_old, high_sep, local_active_shift, &
            local_inactive_shift, global_inactive_shift
      END IF

   END SUBROUTINE check_projector_trace_local_compensation

! **************************************************************************************************
!> \brief Check accepted-drop scaling of the physical Direct candidate.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_direct_history_step_scale(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      REAL(KIND=dp), PARAMETER                           :: max_quality = 1.25_dp, &
                                                            scale_floor = 0.25_dp

      REAL(KIND=dp) :: capped_growth_scale, err, floor_scale, good_scale, growth_scale, &
         mature_floor_scale, mature_ratio_scale, over_scale, ratio_scale, recent_growth_scale, &
         recent_scale, selection_limited_growth_scale, young_growth_scale, young_scale

      young_scale = dense_direct_history_step_scale(4.0_dp, 1.0_dp, 0.25_dp, 3, 4, &
                                                    8, scale_floor, max_quality)
      ratio_scale = dense_direct_history_step_scale(10.0_dp, 4.0_dp, 0.20_dp, 4, 4, &
                                                    8, scale_floor, max_quality)
      floor_scale = dense_direct_history_step_scale(10.0_dp, 0.01_dp, 0.02_dp, 4, 4, &
                                                    8, scale_floor, max_quality)
      mature_ratio_scale = dense_direct_history_step_scale(10.0_dp, 4.0_dp, 0.20_dp, 16, 4, &
                                                           8, scale_floor, max_quality)
      mature_floor_scale = dense_direct_history_step_scale(10.0_dp, 0.01_dp, 0.02_dp, 16, 4, &
                                                           8, scale_floor, max_quality)
      good_scale = dense_direct_history_step_scale(10.0_dp, 10.0_dp, 0.70_dp, 4, 4, &
                                                   8, scale_floor, max_quality)
      over_scale = dense_direct_history_step_scale(10.0_dp, 20.0_dp, 1.20_dp, 4, 4, &
                                                   8, scale_floor, max_quality)
      recent_scale = dense_direct_history_step_scale(10.0_dp, 4.0_dp, 1.20_dp, 16, 4, &
                                                     8, scale_floor, max_quality)
      young_growth_scale = dense_direct_history_growth_scale(4.0_dp, 1.0_dp, 0.25_dp, 3, 4, &
                                                             8, max_quality)
      growth_scale = dense_direct_history_growth_scale(10.0_dp, 20.0_dp, 1.20_dp, 4, 4, &
                                                       8, max_quality)
      selection_limited_growth_scale = dense_direct_history_growth_scale(10.0_dp, 4.0_dp, &
                                                                         0.40_dp, 16, 4, &
                                                                         8, max_quality)
      recent_growth_scale = dense_direct_history_growth_scale(10.0_dp, 4.0_dp, 1.20_dp, &
                                                              16, 4, 8, max_quality)
      capped_growth_scale = dense_direct_history_growth_scale(1.0_dp, 10.0_dp, 5.0_dp, 16, 4, &
                                                              8, max_quality)

      err = ABS(young_scale - 1.0_dp)
      err = MAX(err, ABS(ratio_scale - &
                         (1.0_dp - (1.0_dp/13.0_dp)*(1.0_dp - SQRT(0.40_dp)))))
      err = MAX(err, ABS(floor_scale - &
                         (1.0_dp - (1.0_dp/13.0_dp)*(1.0_dp - SQRT(scale_floor)))))
      err = MAX(err, ABS(mature_ratio_scale - SQRT(0.40_dp)))
      err = MAX(err, ABS(mature_floor_scale - SQRT(scale_floor)))
      err = MAX(err, ABS(good_scale - 1.0_dp))
      err = MAX(err, ABS(over_scale - 1.0_dp))
      err = MAX(err, ABS(recent_scale - 1.0_dp))
      err = MAX(err, ABS(young_growth_scale - 1.0_dp))
      err = MAX(err, ABS(growth_scale - &
                         (1.0_dp + (1.0_dp/13.0_dp)*(SQRT(max_quality) - 1.0_dp))))
      err = MAX(err, ABS(selection_limited_growth_scale - 1.0_dp))
      err = MAX(err, ABS(recent_growth_scale - SQRT(1.20_dp)))
      err = MAX(err, ABS(capped_growth_scale - SQRT(max_quality)))

      IF (err > tol .OR. ratio_scale >= 1.0_dp .OR. floor_scale >= ratio_scale .OR. &
          mature_ratio_scale >= ratio_scale .OR. mature_floor_scale >= floor_scale .OR. &
          recent_scale < 1.0_dp .OR. &
          growth_scale <= selection_limited_growth_scale .OR. &
          recent_growth_scale <= selection_limited_growth_scale .OR. &
          capped_growth_scale <= growth_scale) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,14(1X,ES12.5))') TRIM(label), &
            "failed with error/young/ratio/floor/mratio/mfloor/good/over/recent/ygrow/grow/selgrow/rgrow/cgrow", &
            err, young_scale, ratio_scale, floor_scale, mature_ratio_scale, mature_floor_scale, &
            good_scale, over_scale, recent_scale, young_growth_scale, growth_scale, &
            selection_limited_growth_scale, recent_growth_scale, capped_growth_scale
      ELSE
         WRITE (*, '(A,1X,A,14(1X,ES12.5))') TRIM(label), &
            "OK with error/young/ratio/floor/mratio/mfloor/good/over/recent/ygrow/grow/selgrow/rgrow/cgrow", &
            err, young_scale, ratio_scale, floor_scale, mature_ratio_scale, mature_floor_scale, &
            good_scale, over_scale, recent_scale, young_growth_scale, growth_scale, &
            selection_limited_growth_scale, recent_growth_scale, capped_growth_scale
      END IF

   END SUBROUTINE check_projector_direct_history_step_scale

! **************************************************************************************************
!> \brief Check that residual-tail Direct selection does not double-count correlated histories.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_direct_tail_selection_quality(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      REAL(KIND=dp)                                      :: err, no_history, product_collapse, &
                                                            shadow_floor, shadow_quality, &
                                                            stale_selected_quality, trace_quality

      trace_quality = 7.55E-2_dp
      shadow_quality = 7.53E-2_dp
      stale_selected_quality = 3.6E-4_dp

      shadow_floor = dense_direct_tail_selection_quality(trace_quality, &
                                                         MAX(shadow_quality, stale_selected_quality))
      product_collapse = trace_quality*stale_selected_quality
      no_history = dense_direct_tail_selection_quality(trace_quality, 1.0_dp)

      err = ABS(shadow_floor - shadow_quality)
      err = MAX(err, ABS(no_history - trace_quality))

      IF (err > tol .OR. shadow_floor <= 100.0_dp*product_collapse .OR. &
          shadow_floor > trace_quality) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,5(1X,ES12.5))') TRIM(label), &
            "failed with error/trace/shadow/stale/product", err, trace_quality, shadow_floor, &
            stale_selected_quality, product_collapse
      ELSE
         WRITE (*, '(A,1X,A,5(1X,ES12.5))') TRIM(label), &
            "OK with error/trace/shadow/stale/product", err, trace_quality, shadow_floor, &
            stale_selected_quality, product_collapse
      END IF

   END SUBROUTINE check_projector_direct_tail_selection_quality

! **************************************************************************************************
!> \brief Check that scaled Direct steps use the Schur quadratic, not the unscaled RHS product.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_direct_scaled_metric_norm(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 3
      REAL(KIND=dp), PARAMETER                           :: damp_scale = 0.55_dp, &
                                                            grow_scale = 1.15_dp

      REAL(KIND=dp) :: damp_expected, damp_new_drop, damp_old_drop, damp_rhs_quad, &
         damp_schur_quad, err, grow_expected, grow_new_drop, grow_old_drop, grow_rhs_quad, &
         grow_schur_quad
      REAL(KIND=dp), DIMENSION(n)                        :: q_base, q_damp, q_grow, rhs
      REAL(KIND=dp), DIMENSION(n, n)                     :: h

      h(:, :) = RESHAPE([0.72_dp, 0.06_dp, -0.04_dp, &
                         0.06_dp, 0.48_dp, 0.05_dp, &
                         -0.04_dp, 0.05_dp, 0.63_dp], [n, n])
      q_base(:) = [0.18_dp, -0.11_dp, 0.07_dp]
      rhs(:) = MATMUL(h, q_base)

      q_damp(:) = damp_scale*q_base(:)
      damp_rhs_quad = DOT_PRODUCT(rhs, q_damp)
      damp_schur_quad = DOT_PRODUCT(q_damp, MATMUL(h, q_damp))
      damp_expected = damp_scale*damp_rhs_quad
      damp_old_drop = damp_rhs_quad - 0.5_dp*damp_rhs_quad
      damp_new_drop = damp_rhs_quad - 0.5_dp*damp_schur_quad

      q_grow(:) = grow_scale*q_base(:)
      grow_rhs_quad = DOT_PRODUCT(rhs, q_grow)
      grow_schur_quad = DOT_PRODUCT(q_grow, MATMUL(h, q_grow))
      grow_expected = grow_scale*grow_rhs_quad
      grow_old_drop = grow_rhs_quad - 0.5_dp*grow_rhs_quad
      grow_new_drop = grow_rhs_quad - 0.5_dp*grow_schur_quad

      err = ABS(damp_schur_quad - damp_expected)
      err = MAX(err, ABS(grow_schur_quad - grow_expected))

      IF (err > tol .OR. damp_new_drop <= damp_old_drop .OR. grow_new_drop >= grow_old_drop) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,7(1X,ES12.5))') TRIM(label), &
            "failed with error/drhs/dschur/ddrop/grhs/gschur/gdrop", err, &
            damp_rhs_quad, damp_schur_quad, damp_new_drop - damp_old_drop, &
            grow_rhs_quad, grow_schur_quad, grow_new_drop - grow_old_drop
      ELSE
         WRITE (*, '(A,1X,A,7(1X,ES12.5))') TRIM(label), &
            "OK with error/drhs/dschur/ddrop/grhs/gschur/gdrop", err, &
            damp_rhs_quad, damp_schur_quad, damp_new_drop - damp_old_drop, &
            grow_rhs_quad, grow_schur_quad, grow_new_drop - grow_old_drop
      END IF

   END SUBROUTINE check_projector_direct_scaled_metric_norm

! **************************************************************************************************
!> \brief Check that spectral candidates keep model curvature unscaled from the accepted step norm.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_spectral_metric_norm_scale(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      REAL(KIND=dp), PARAMETER                           :: gdotp = 0.46_dp, gnorm = 2.5_dp, &
                                                            metric_energy = 0.21_dp, &
                                                            metric_rotation = 0.34_dp, &
                                                            scale = 0.37_dp

      REAL(KIND=dp)                                      :: actual_normp, err, model_drop, &
                                                            model_normp, norm_ratio, wrong_drop

      model_normp = metric_energy + metric_rotation
      actual_normp = scale*scale*model_normp
      model_drop = scale*gdotp - 0.5_dp*actual_normp
      norm_ratio = actual_normp/gnorm

      ! Reusing the accepted-step norm as the unscaled model curvature would
      ! apply the line-search scale twice and overestimate the spectral drop.
      wrong_drop = scale*gdotp - 0.5_dp*scale*scale*actual_normp
      err = ABS(model_drop - (scale*gdotp - 0.5_dp*scale*scale*model_normp))
      err = MAX(err, ABS(norm_ratio - scale*scale*model_normp/gnorm))

      IF (err > tol .OR. wrong_drop <= model_drop .OR. actual_normp >= model_normp) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,5(1X,ES12.5))') TRIM(label), &
            "failed with error/model/actual/drop/wrong", err, model_normp, actual_normp, &
            model_drop, wrong_drop
      ELSE
         WRITE (*, '(A,1X,A,5(1X,ES12.5))') TRIM(label), &
            "OK with error/model/actual/drop/wrong", err, model_normp, actual_normp, &
            model_drop, wrong_drop
      END IF

   END SUBROUTINE check_projector_spectral_metric_norm_scale

! **************************************************************************************************
!> \brief Check that mature spectral shadow history can promote a good early selected step.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_spectral_selected_history_promotion(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      LOGICAL                                            :: early_promotes, failure_veto, &
                                                            no_selected_veto, recent_veto, &
                                                            shadow_veto
      REAL(KIND=dp)                                      :: err

      early_promotes = dense_spectral_selected_history_promotion( &
         1, 0, 0.67_dp, 0.076_dp, 0.076_dp, 0.076_dp, 2.0E-6_dp, 1.2E-7_dp, &
         16, 2.0E-6_dp, 1.2E-7_dp, 0.31_dp, 0.12_dp, 3.2E-4_dp, 1.0_dp, 1.0_dp)
      no_selected_veto = dense_spectral_selected_history_promotion( &
         0, 0, 0.67_dp, 0.076_dp, 0.076_dp, 0.076_dp, 2.0E-6_dp, 1.2E-7_dp, &
         16, 2.0E-6_dp, 1.2E-7_dp, 0.31_dp, 0.12_dp, 3.2E-4_dp, 1.0_dp, 1.0_dp)
      recent_veto = dense_spectral_selected_history_promotion( &
         1, 0, 0.67_dp, 0.020_dp, 0.076_dp, 0.020_dp, 2.0E-6_dp, 1.2E-7_dp, &
         16, 2.0E-6_dp, 1.2E-7_dp, 0.31_dp, 0.12_dp, 3.2E-4_dp, 1.0_dp, 1.0_dp)
      shadow_veto = dense_spectral_selected_history_promotion( &
         1, 0, 0.67_dp, 0.076_dp, 0.076_dp, 0.076_dp, 2.0E-6_dp, 1.2E-7_dp, &
         7, 2.0E-6_dp, 1.2E-7_dp, 0.31_dp, 0.12_dp, 3.2E-4_dp, 1.0_dp, 1.0_dp)
      failure_veto = dense_spectral_selected_history_promotion( &
         1, 2, 0.67_dp, 0.076_dp, 0.076_dp, 0.076_dp, 2.0E-6_dp, 1.2E-7_dp, &
         16, 2.0E-6_dp, 1.2E-7_dp, 0.31_dp, 0.12_dp, 3.2E-4_dp, 1.0_dp, 1.0_dp)

      err = 0.0_dp
      IF (.NOT. early_promotes .OR. no_selected_veto .OR. recent_veto .OR. shadow_veto .OR. &
          failure_veto) err = 1.0_dp

      IF (err > tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,5(1X,L1))') TRIM(label), &
            "failed with early/no-selected/recent/shadow/failure", early_promotes, &
            no_selected_veto, recent_veto, shadow_veto, failure_veto
      ELSE
         WRITE (*, '(A,1X,A,5(1X,L1))') TRIM(label), &
            "OK with early/no-selected/recent/shadow/failure", early_promotes, &
            no_selected_veto, recent_veto, shadow_veto, failure_veto
      END IF

   END SUBROUTINE check_projector_spectral_selected_history_promotion

! **************************************************************************************************
!> \brief Check accepted-drop consistency scaling of mixed Direct history.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_direct_history_consistency_scale(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      REAL(KIND=dp), PARAMETER                           :: max_quality = 1.25_dp

      REAL(KIND=dp) :: both_high_quality, both_high_scale, boundary_scale, ema_high_quality, &
         ema_high_scale, err, good_scale, inactive_scale, raw_high_quality, raw_high_scale

      good_scale = dense_direct_history_consistency_scale(1.0_dp, 0.8_dp, max_quality)
      boundary_scale = dense_direct_history_consistency_scale(max_quality, max_quality, &
                                                              max_quality)
      raw_high_scale = dense_direct_history_consistency_scale(5.0_dp, 1.0_dp, max_quality)
      ema_high_scale = dense_direct_history_consistency_scale(1.0_dp, 20.0_dp, max_quality)
      both_high_scale = dense_direct_history_consistency_scale(5.0_dp, 20.0_dp, max_quality)
      inactive_scale = dense_direct_history_consistency_scale(-1.0_dp, 0.0_dp, max_quality)
      raw_high_quality = raw_high_scale*MIN(1.0_dp, 5.0_dp)
      ema_high_quality = ema_high_scale*MIN(1.0_dp, 20.0_dp)
      both_high_quality = both_high_scale*MIN(1.0_dp, 20.0_dp)

      err = ABS(good_scale - 1.0_dp)
      err = MAX(err, ABS(boundary_scale - 1.0_dp))
      err = MAX(err, ABS(raw_high_scale - SQRT(max_quality/5.0_dp)))
      err = MAX(err, ABS(ema_high_scale - SQRT(max_quality/20.0_dp)))
      err = MAX(err, ABS(both_high_scale - ema_high_scale))
      err = MAX(err, ABS(inactive_scale - 1.0_dp))
      err = MAX(err, ABS(raw_high_quality - raw_high_scale))
      err = MAX(err, ABS(ema_high_quality - ema_high_scale))
      err = MAX(err, ABS(both_high_quality - both_high_scale))

      IF (err > tol .OR. raw_high_scale >= good_scale .OR. &
          ema_high_scale >= raw_high_scale .OR. both_high_scale /= ema_high_scale .OR. &
          raw_high_quality >= good_scale .OR. ema_high_quality >= raw_high_quality .OR. &
          both_high_quality /= ema_high_quality) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,10(1X,ES12.5))') TRIM(label), &
            "failed with error/good/boundary/raw/ema/both/inactive/rawq/emaq/bothq", err, &
            good_scale, boundary_scale, raw_high_scale, ema_high_scale, both_high_scale, &
            inactive_scale, raw_high_quality, ema_high_quality, both_high_quality
      ELSE
         WRITE (*, '(A,1X,A,10(1X,ES12.5))') TRIM(label), &
            "OK with error/good/boundary/raw/ema/both/inactive/rawq/emaq/bothq", err, good_scale, &
            boundary_scale, raw_high_scale, ema_high_scale, both_high_scale, inactive_scale, &
            raw_high_quality, ema_high_quality, both_high_quality
      END IF

   END SUBROUTINE check_projector_direct_history_consistency_scale

! **************************************************************************************************
!> \brief Check accepted-history promotion of mixed Direct curvature.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_mixed_direct_promotion(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: min_samples = 8
      REAL(KIND=dp), PARAMETER                           :: min_quality = 0.10_dp

      LOGICAL :: aggregate_promoted, ema_low_blocked, first_window_blocked, long_shadow_promoted, &
         low_blocked, selected_young_blocked, tail_selected_blocked, tail_shadow_blocked, &
         young_blocked
      REAL(KIND=dp)                                      :: promoted_value

      aggregate_promoted = dense_mixed_direct_promoted(42, min_samples, 0.102135336_dp, &
                                                       0.102135336_dp, 0.0949032371_dp, 0.2_dp, &
                                                       min_samples, min_quality)
      long_shadow_promoted = dense_mixed_direct_promoted(80, min_samples, 0.055_dp, 0.055_dp, &
                                                         0.055_dp, 0.0_dp, 0, min_quality)
      tail_selected_blocked = dense_mixed_direct_promoted(16, min_samples, 0.5_dp, 0.5_dp, &
                                                          0.5_dp, 0.5_dp, min_samples, min_quality)
      tail_shadow_blocked = dense_mixed_direct_promoted(16, min_samples, 0.5_dp, 0.5_dp, &
                                                        0.5_dp, 0.0_dp, 0, min_quality)
      young_blocked = dense_mixed_direct_promoted(7, min_samples, 0.5_dp, 0.5_dp, 0.5_dp, 0.5_dp, &
                                                  min_samples, min_quality)
      first_window_blocked = dense_mixed_direct_promoted(15, min_samples, 0.5_dp, 0.5_dp, &
                                                         0.5_dp, 0.5_dp, min_samples, min_quality)
      selected_young_blocked = dense_mixed_direct_promoted(42, min_samples, 0.5_dp, 0.5_dp, &
                                                           0.5_dp, 0.5_dp, min_samples - 1, min_quality)
      low_blocked = dense_mixed_direct_promoted(42, min_samples, 0.099_dp, 0.11_dp, 0.11_dp, &
                                                0.2_dp, min_samples, min_quality)
      ema_low_blocked = dense_mixed_direct_promoted(42, min_samples, 0.11_dp, 0.11_dp, &
                                                    0.049_dp, 0.2_dp, min_samples, min_quality)
      promoted_value = MERGE(1.0_dp, 0.0_dp, aggregate_promoted)

      IF (.NOT. aggregate_promoted .OR. .NOT. long_shadow_promoted .OR. tail_selected_blocked .OR. &
          tail_shadow_blocked .OR. young_blocked .OR. first_window_blocked .OR. &
          selected_young_blocked .OR. low_blocked .OR. ema_low_blocked .OR. &
          ABS(promoted_value - 1.0_dp) > tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,9(1X,L1),1(1X,ES12.5))') TRIM(label), &
            "failed with aggregate/long_shadow/tail/shot/young/window/selyoung/low/ema/value", &
            aggregate_promoted, long_shadow_promoted, tail_selected_blocked, tail_shadow_blocked, young_blocked, &
            first_window_blocked, &
            selected_young_blocked, low_blocked, ema_low_blocked, promoted_value
      ELSE
         WRITE (*, '(A,1X,A,9(1X,L1),1(1X,ES12.5))') TRIM(label), &
            "OK with aggregate/long_shadow/tail/shot/young/window/selyoung/low/ema/value", &
            aggregate_promoted, long_shadow_promoted, tail_selected_blocked, tail_shadow_blocked, young_blocked, &
            first_window_blocked, &
            selected_young_blocked, low_blocked, ema_low_blocked, promoted_value
      END IF

   END SUBROUTINE check_projector_mixed_direct_promotion

! **************************************************************************************************
!> \brief Check accepted-drop mismatch weighting of the finite Projector-Hessian tail.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_coupled_tail_weight(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      REAL(KIND=dp), PARAMETER                           :: density_fraction = 0.80_dp, &
                                                            diag_avg = 0.50_dp, &
                                                            diag_delta = 0.25_dp

      REAL(KIND=dp) :: above_tail_weight, curvature_floor_weight, diag_damped_weight, err, &
         expected_curvature_floor_weight, expected_mismatch_weight, expected_saturated_weight, &
         good_weight, high_weight, low_weight, saturated_weight, scaled_weight, young_weight

      above_tail_weight = dense_projector_coupled_tail_weight(2.0E-4_dp, 1.0_dp, 0, &
                                                              density_fraction, diag_avg, diag_delta)
      young_weight = dense_projector_coupled_tail_weight(5.0E-5_dp, 0.50_dp, 3, &
                                                         density_fraction, diag_avg, diag_delta)
      good_weight = dense_projector_coupled_tail_weight(5.0E-5_dp, 1.0_dp, 4, &
                                                        density_fraction, diag_avg, diag_delta)
      low_weight = dense_projector_coupled_tail_weight(5.0E-5_dp, 0.625_dp, 4, &
                                                       density_fraction, diag_avg, diag_delta)
      high_weight = dense_projector_coupled_tail_weight(5.0E-5_dp, 1.75_dp, 4, &
                                                        density_fraction, diag_avg, diag_delta)
      saturated_weight = dense_projector_coupled_tail_weight(5.0E-5_dp, 10.0_dp, 4, &
                                                             density_fraction, diag_avg, diag_delta)
      diag_damped_weight = dense_projector_coupled_tail_weight(5.0E-5_dp, 1.75_dp, 4, &
                                                               density_fraction, diag_avg, 1.50_dp)
      scaled_weight = dense_projector_coupled_tail_weight(5.0E-5_dp, 1.75_dp, 4, &
                                                          0.5_dp*density_fraction, diag_avg, diag_delta)
      curvature_floor_weight = dense_projector_coupled_tail_weight(5.0E-5_dp, 1.75_dp, 4, &
                                                                   1.0E-3_dp, diag_avg, diag_delta)

      expected_mismatch_weight = 0.10_dp*0.50_dp*density_fraction/1.50_dp
      expected_saturated_weight = 0.10_dp*density_fraction/1.50_dp
      expected_curvature_floor_weight = 0.10_dp*0.50_dp*(0.25_dp*(0.50_dp/1.50_dp))/1.50_dp

      err = ABS(above_tail_weight - 1.0_dp)
      err = MAX(err, ABS(young_weight))
      err = MAX(err, ABS(good_weight))
      err = MAX(err, ABS(low_weight - expected_mismatch_weight))
      err = MAX(err, ABS(high_weight - expected_mismatch_weight))
      err = MAX(err, ABS(saturated_weight - expected_saturated_weight))
      err = MAX(err, ABS(scaled_weight - 0.5_dp*expected_mismatch_weight))
      err = MAX(err, ABS(curvature_floor_weight - expected_curvature_floor_weight))

      IF (err > tol .OR. high_weight <= good_weight .OR. low_weight <= good_weight .OR. &
          saturated_weight < high_weight .OR. diag_damped_weight >= high_weight .OR. &
          curvature_floor_weight <= scaled_weight*1.0E-3_dp/density_fraction) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,9(1X,ES12.5))') TRIM(label), &
            "failed with error/above/young/good/low/high/saturated/diag/curv", err, above_tail_weight, &
            young_weight, good_weight, low_weight, high_weight, saturated_weight, diag_damped_weight, &
            curvature_floor_weight
      ELSE
         WRITE (*, '(A,1X,A,9(1X,ES12.5))') TRIM(label), &
            "OK with error/above/young/good/low/high/saturated/diag/curv", err, above_tail_weight, &
            young_weight, good_weight, low_weight, high_weight, saturated_weight, diag_damped_weight, &
            curvature_floor_weight
      END IF

   END SUBROUTINE check_projector_coupled_tail_weight

! **************************************************************************************************
!> \brief Check conversion from absolute density-history weight to Projector-delta blend fraction.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_tail_active_delta_weight(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      REAL(KIND=dp)                                      :: err, full_blend, no_density, no_projector, &
                                                            partial_blend, saturated_blend

      no_projector = dense_projector_tail_active_delta_weight(0.0_dp, 0.05_dp)
      no_density = dense_projector_tail_active_delta_weight(0.10_dp, 0.0_dp)
      partial_blend = dense_projector_tail_active_delta_weight(0.10_dp, 0.025_dp)
      full_blend = dense_projector_tail_active_delta_weight(0.10_dp, 0.10_dp)
      saturated_blend = dense_projector_tail_active_delta_weight(0.05_dp, 0.20_dp)

      err = ABS(no_projector)
      err = MAX(err, ABS(no_density))
      err = MAX(err, ABS(partial_blend - 0.25_dp))
      err = MAX(err, ABS(full_blend - 1.0_dp))
      err = MAX(err, ABS(saturated_blend - 1.0_dp))

      IF (err > tol .OR. partial_blend <= no_density .OR. full_blend < partial_blend .OR. &
          saturated_blend < full_blend) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,6(1X,ES12.5))') TRIM(label), &
            "failed with error/noproj/nodens/partial/full/saturated", err, no_projector, &
            no_density, partial_blend, full_blend, saturated_blend
      ELSE
         WRITE (*, '(A,1X,A,6(1X,ES12.5))') TRIM(label), &
            "OK with error/noproj/nodens/partial/full/saturated", err, no_projector, &
            no_density, partial_blend, full_blend, saturated_blend
      END IF

   END SUBROUTINE check_projector_tail_active_delta_weight

! **************************************************************************************************
!> \brief Check that Direct-Delta eigenspaces do not overwrite the production curvature matrix.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_direct_delta_eigen_copy(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 4

      INTEGER                                            :: i
      REAL(KIND=dp)                                      :: err, matrix_quad, mode_quad
      REAL(KIND=dp), DIMENSION(n)                        :: coeff, eval, step
      REAL(KIND=dp), DIMENSION(n, n)                     :: delta, delta_saved, eigvec

      delta(:, :) = RESHAPE([0.35_dp, 0.04_dp, -0.02_dp, 0.03_dp, &
                             0.04_dp, -0.18_dp, 0.07_dp, -0.01_dp, &
                             -0.02_dp, 0.07_dp, 0.12_dp, 0.05_dp, &
                             0.03_dp, -0.01_dp, 0.05_dp, -0.09_dp], [n, n])
      delta(:, :) = 0.5_dp*(delta(:, :) + TRANSPOSE(delta(:, :)))
      delta_saved(:, :) = delta(:, :)
      step(:) = [0.42_dp, -0.31_dp, 0.27_dp, 0.16_dp]

      eigvec(:, :) = delta(:, :)
      CALL diamat_all(eigvec, eval)
      coeff(:) = MATMUL(TRANSPOSE(eigvec), step)
      matrix_quad = DOT_PRODUCT(step, MATMUL(delta, step))
      mode_quad = 0.0_dp
      DO i = 1, n
         mode_quad = mode_quad + eval(i)*coeff(i)*coeff(i)
      END DO

      err = MAXVAL(ABS(delta(:, :) - delta_saved(:, :)))
      err = MAX(err, ABS(matrix_quad - mode_quad))
      IF (err > tol .OR. ABS(matrix_quad) < 1.0E-4_dp .OR. &
          MAXVAL(ABS(delta(:, :) - eigvec(:, :))) < 1.0E-4_dp) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,3(1X,ES12.5))') TRIM(label), &
            "failed with error/matrix_quad/mode_quad", err, matrix_quad, mode_quad
      ELSE
         WRITE (*, '(A,1X,A,3(1X,ES12.5))') TRIM(label), &
            "OK with error/matrix_quad/mode_quad", err, matrix_quad, mode_quad
      END IF

   END SUBROUTINE check_projector_direct_delta_eigen_copy

! **************************************************************************************************
!> \brief Check that active Direct curvature filters local residuals, not fixed-N cross-channel terms.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_active_residual_filter(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 3

      INTEGER                                            :: r, s
      INTEGER, DIMENSION(n)                              :: channel
      REAL(KIND=dp) :: active_cov, active_diag_avg, active_weight, cov, diag_change, err, &
         expected_diag_avg, legacy_diag_avg, local_change, residual
      REAL(KIND=dp), DIMENSION(n)                        :: active_projector_weight, projector_weight
      REAL(KIND=dp), DIMENSION(n, n)                     :: active, expected, physical, rot

      channel(:) = [1, 1, 2]
      projector_weight(:) = [0.50_dp, 0.20_dp, 0.40_dp]
      active_projector_weight(:) = [0.10_dp, 0.00_dp, 0.40_dp]
      rot(:, :) = RESHAPE([0.31_dp, 0.04_dp, -0.02_dp, &
                           0.04_dp, 0.27_dp, 0.03_dp, &
                           -0.02_dp, 0.03_dp, 0.18_dp], [n, n])
      physical(:, :) = rot(:, :)
      physical(1, 1) = physical(1, 1) + 0.35_dp
      physical(1, 2) = physical(1, 2) - 0.12_dp
      physical(2, 1) = physical(1, 2)
      physical(2, 2) = physical(2, 2) + 0.08_dp
      physical(1, 3) = physical(1, 3) + 0.21_dp
      physical(3, 1) = physical(1, 3)

      active(:, :) = physical(:, :)
      expected(:, :) = physical(:, :)
      DO s = 1, n
         DO r = 1, n
            IF (channel(r) /= channel(s)) CYCLE
            cov = -0.5_dp*(projector_weight(r) + projector_weight(s))*rot(r, s)
            active_cov = -0.5_dp*(active_projector_weight(r) + active_projector_weight(s))*rot(r, s)
            active_weight = 0.0_dp
            IF (projector_weight(r) + projector_weight(s) > EPSILON(active_weight)) &
               active_weight = (active_projector_weight(r) + active_projector_weight(s))/ &
                               MAX(projector_weight(r) + projector_weight(s), EPSILON(1.0_dp))
            active_weight = MIN(1.0_dp, MAX(0.0_dp, active_weight))
            residual = physical(r, s) - rot(r, s) - cov
            active(r, s) = rot(r, s) + active_cov + active_weight*residual
            expected(r, s) = active(r, s)
         END DO
      END DO
      active_diag_avg = 0.0_dp
      expected_diag_avg = 0.0_dp
      legacy_diag_avg = 0.0_dp
      DO r = 1, n
         cov = -projector_weight(r)*rot(r, r)
         active_cov = -active_projector_weight(r)*rot(r, r)
         active_weight = 0.0_dp
         IF (projector_weight(r) > EPSILON(active_weight)) &
            active_weight = active_projector_weight(r)/MAX(projector_weight(r), EPSILON(1.0_dp))
         active_weight = MIN(1.0_dp, MAX(0.0_dp, active_weight))
         residual = physical(r, r) - rot(r, r) - cov
         active_diag_avg = active_diag_avg + active(r, r) - rot(r, r)
         expected_diag_avg = expected_diag_avg + active_cov + active_weight*residual
         legacy_diag_avg = legacy_diag_avg + active_cov + residual
      END DO
      active_diag_avg = active_diag_avg/REAL(n, KIND=dp)
      expected_diag_avg = expected_diag_avg/REAL(n, KIND=dp)
      legacy_diag_avg = legacy_diag_avg/REAL(n, KIND=dp)

      err = MAXVAL(ABS(active(:, :) - expected(:, :)))
      err = MAX(err, ABS(active_diag_avg - expected_diag_avg))
      err = MAX(err, ABS(active(1, 3) - physical(1, 3)))
      local_change = ABS(active(1, 2) - physical(1, 2))
      diag_change = ABS(active(1, 1) - physical(1, 1))

      IF (err > tol .OR. local_change < 1.0E-3_dp .OR. diag_change < 1.0E-3_dp .OR. &
          ABS(active(1, 3) - physical(1, 3)) > tol .OR. &
          ABS(legacy_diag_avg - active_diag_avg) < 1.0E-3_dp) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,7(1X,ES12.5))') TRIM(label), &
            "failed with error/local/diag/full/cross/active_avg/legacy_avg", &
            err, local_change, diag_change, physical(1, 2), active(1, 3), &
            active_diag_avg, legacy_diag_avg
      ELSE
         WRITE (*, '(A,1X,A,7(1X,ES12.5))') TRIM(label), &
            "OK with error/local/diag/full/cross/active_avg/legacy_avg", &
            err, local_change, diag_change, physical(1, 2), active(1, 3), &
            active_diag_avg, legacy_diag_avg
      END IF

   END SUBROUTINE check_projector_active_residual_filter

! **************************************************************************************************
!> \brief Check that the active Direct shadow is distinct from the full physical delta.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_active_direct_delta_shadow(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 3

      INTEGER                                            :: r, s
      INTEGER, DIMENSION(n)                              :: channel
      REAL(KIND=dp)                                      :: active_cov, active_diag_avg, &
                                                            active_weight, cov, direct_diag_avg, &
                                                            err, full_diag_avg, residual
      REAL(KIND=dp), DIMENSION(n)                        :: active_projector_weight, projector_weight
      REAL(KIND=dp), DIMENSION(n, n)                     :: active, direct_delta, full_delta, &
                                                            physical, rot

      channel(:) = [1, 1, 2]
      projector_weight(:) = [0.50_dp, 0.20_dp, 0.40_dp]
      active_projector_weight(:) = [0.10_dp, 0.00_dp, 0.40_dp]
      rot(:, :) = RESHAPE([0.31_dp, 0.04_dp, -0.02_dp, &
                           0.04_dp, 0.27_dp, 0.03_dp, &
                           -0.02_dp, 0.03_dp, 0.18_dp], [n, n])
      physical(:, :) = rot(:, :)
      physical(1, 1) = physical(1, 1) + 0.35_dp
      physical(1, 2) = physical(1, 2) - 0.12_dp
      physical(2, 1) = physical(1, 2)
      physical(2, 2) = physical(2, 2) + 0.08_dp
      physical(1, 3) = physical(1, 3) + 0.21_dp
      physical(3, 1) = physical(1, 3)

      active(:, :) = physical(:, :)
      DO s = 1, n
         DO r = 1, n
            IF (channel(r) /= channel(s)) CYCLE
            cov = -0.5_dp*(projector_weight(r) + projector_weight(s))*rot(r, s)
            active_cov = -0.5_dp*(active_projector_weight(r) + active_projector_weight(s))*rot(r, s)
            active_weight = 0.0_dp
            IF (projector_weight(r) + projector_weight(s) > EPSILON(active_weight)) &
               active_weight = (active_projector_weight(r) + active_projector_weight(s))/ &
                               MAX(projector_weight(r) + projector_weight(s), EPSILON(1.0_dp))
            active_weight = MIN(1.0_dp, MAX(0.0_dp, active_weight))
            residual = physical(r, s) - rot(r, s) - cov
            active(r, s) = rot(r, s) + active_cov + active_weight*residual
         END DO
      END DO

      full_delta(:, :) = 0.5_dp*((physical(:, :) - rot(:, :)) + TRANSPOSE(physical(:, :) - rot(:, :)))
      direct_delta(:, :) = 0.5_dp*((active(:, :) - rot(:, :)) + TRANSPOSE(active(:, :) - rot(:, :)))
      full_diag_avg = 0.0_dp
      active_diag_avg = 0.0_dp
      direct_diag_avg = 0.0_dp
      DO r = 1, n
         full_diag_avg = full_diag_avg + full_delta(r, r)
         active_diag_avg = active_diag_avg + active(r, r) - rot(r, r)
         direct_diag_avg = direct_diag_avg + direct_delta(r, r)
      END DO
      full_diag_avg = full_diag_avg/REAL(n, KIND=dp)
      active_diag_avg = active_diag_avg/REAL(n, KIND=dp)
      direct_diag_avg = direct_diag_avg/REAL(n, KIND=dp)

      err = ABS(direct_diag_avg - active_diag_avg)
      err = MAX(err, ABS(direct_delta(1, 3) - full_delta(1, 3)))

      IF (err > tol .OR. ABS(full_diag_avg - direct_diag_avg) < 1.0E-3_dp .OR. &
          ABS(direct_delta(1, 2) - full_delta(1, 2)) < 1.0E-3_dp) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,5(1X,ES12.5))') TRIM(label), &
            "failed with error/direct/full/active/local", &
            err, direct_diag_avg, full_diag_avg, active_diag_avg, direct_delta(1, 2)
      ELSE
         WRITE (*, '(A,1X,A,5(1X,ES12.5))') TRIM(label), &
            "OK with error/direct/full/active/local", &
            err, direct_diag_avg, full_diag_avg, active_diag_avg, direct_delta(1, 2)
      END IF

   END SUBROUTINE check_projector_active_direct_delta_shadow

! **************************************************************************************************
!> \brief Check that Direct trace mismatch contributes even when projected response vectors match.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_direct_trace_model_gap(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      REAL(KIND=dp), PARAMETER                           :: scale = 0.74_dp, &
                                                            trace_active = -0.23_dp, &
                                                            trace_actual = 0.04_dp, &
                                                            trace_invariant = 0.86_dp

      REAL(KIND=dp) :: err, expected_mixed_trust, expected_trace_gap, invariant_gap, &
         invariant_trust, mixed_trust, mixed_trust_gap, model_drop, trace_gap, vector_gap
      REAL(KIND=dp), DIMENSION(4)                        :: curv, response_actual, response_model

      curv(:) = [0.55_dp, 0.23_dp, 0.12_dp, 0.07_dp]
      response_actual(:) = [0.02_dp, -0.01_dp, 0.03_dp, -0.04_dp]
      response_model(:) = response_actual(:)

      vector_gap = 0.5_dp*scale*scale*SUM(curv(:)*(response_model(:) - response_actual(:))**2)
      trace_gap = dense_direct_trace_model_gap_drop(trace_actual, trace_active, scale)
      invariant_gap = dense_direct_trace_model_gap_drop(trace_actual, trace_invariant, scale)
      mixed_trust_gap = MIN(trace_gap, invariant_gap)
      expected_trace_gap = 0.5_dp*scale*scale*ABS(trace_actual - trace_active)
      model_drop = 1.0_dp
      mixed_trust = dense_direct_trace_mismatch_trust(model_drop, mixed_trust_gap, 0.10_dp, 16)
      invariant_trust = dense_direct_trace_mismatch_trust(model_drop, invariant_gap, 0.10_dp, 16)
      expected_mixed_trust = model_drop/(model_drop + 0.90_dp*trace_gap)

      err = ABS(trace_gap - expected_trace_gap)
      err = MAX(err, ABS(mixed_trust - expected_mixed_trust))
      IF (err > tol .OR. vector_gap > tol .OR. trace_gap <= 0.0_dp .OR. &
          trace_gap/model_drop <= 0.05_dp .OR. ABS(mixed_trust_gap - trace_gap) > tol .OR. &
          mixed_trust <= invariant_trust) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,8(1X,ES12.5))') TRIM(label), &
            "failed with error/vector/trace/ratio/invgap/tgap/mtrust/itrust", &
            err, vector_gap, trace_gap, trace_gap/model_drop, invariant_gap, mixed_trust_gap, &
            mixed_trust, invariant_trust
      ELSE
         WRITE (*, '(A,1X,A,8(1X,ES12.5))') TRIM(label), &
            "OK with error/vector/trace/ratio/invgap/tgap/mtrust/itrust", &
            err, vector_gap, trace_gap, trace_gap/model_drop, invariant_gap, mixed_trust_gap, &
            mixed_trust, invariant_trust
      END IF

   END SUBROUTINE check_projector_direct_trace_model_gap

! **************************************************************************************************
!> \brief Check that mixed Direct response diagnostics use the same blend as the Hessian.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_mixed_direct_response_blend(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: ne = 4, nrot = 3
      REAL(KIND=dp), PARAMETER                           :: direct_weight = 0.40_dp, &
                                                            z_weight = 0.35_dp

      REAL(KIND=dp)                                      :: endpoint_err, err, full_err, legacy_gap, &
                                                            mixed_gap, zero_err
      REAL(KIND=dp), DIMENSION(ne)                       :: actual, curv, expected, model_full, &
                                                            model_mixed, model_zero
      REAL(KIND=dp), DIMENSION(ne, nrot)                 :: density_direct, density_floor, &
                                                            density_raw, physical, raw
      REAL(KIND=dp), DIMENSION(nrot)                     :: step

      curv(:) = [0.42_dp, 0.18_dp, 0.31_dp, 0.07_dp]
      step(:) = [0.27_dp, -0.34_dp, 0.18_dp]
      raw(:, :) = RESHAPE([ &
                          0.08_dp, -0.03_dp, 0.02_dp, -0.01_dp, &
                          -0.04_dp, 0.06_dp, -0.02_dp, 0.03_dp, &
                          0.01_dp, -0.05_dp, 0.04_dp, 0.02_dp], [ne, nrot])
      physical(:, :) = raw(:, :)
      physical(:, 1) = physical(:, 1) + [0.18_dp, -0.09_dp, 0.04_dp, -0.06_dp]
      physical(:, 2) = physical(:, 2) + [-0.07_dp, 0.13_dp, 0.08_dp, -0.02_dp]
      physical(:, 3) = physical(:, 3) + [0.05_dp, 0.02_dp, -0.11_dp, 0.09_dp]
      density_raw(:, :) = RESHAPE([ &
                                  -0.03_dp, 0.06_dp, -0.02_dp, 0.04_dp, &
                                  0.07_dp, -0.01_dp, 0.05_dp, -0.03_dp, &
                                  -0.02_dp, 0.04_dp, -0.06_dp, 0.01_dp], [ne, nrot])
      density_floor(:, :) = RESHAPE([ &
                                    0.04_dp, -0.01_dp, 0.05_dp, -0.02_dp, &
                                    -0.02_dp, 0.08_dp, -0.04_dp, 0.06_dp, &
                                    0.03_dp, -0.07_dp, 0.02_dp, -0.05_dp], [ne, nrot])
      density_direct(:, :) = (1.0_dp - direct_weight)*density_raw(:, :) + &
                             direct_weight*density_floor(:, :)

      model_zero(:) = MATMUL(raw, step)
      model_full(:) = MATMUL(physical, step) + MATMUL(density_direct, step)
      expected(:) = (1.0_dp - z_weight)*MATMUL(raw, step) + &
                    z_weight*MATMUL(physical, step) + z_weight*MATMUL(density_direct, step)
      model_mixed(:) = MATMUL(raw, step)
      endpoint_err = MAXVAL(ABS(model_mixed(:) - model_zero(:)))
      model_mixed(:) = MATMUL(physical, step) + MATMUL(density_direct, step)
      endpoint_err = MAX(endpoint_err, MAXVAL(ABS(model_mixed(:) - model_full(:))))
      model_mixed(:) = expected(:)
      actual(:) = expected(:)

      err = MAXVAL(ABS(model_mixed(:) - expected(:)))
      zero_err = MAXVAL(ABS(model_zero(:) - MATMUL(raw, step)))
      full_err = MAXVAL(ABS(model_full(:) - MATMUL(physical, step) - MATMUL(density_direct, step)))
      mixed_gap = 0.5_dp*SUM(curv(:)*(model_mixed(:) - actual(:))**2)
      legacy_gap = 0.5_dp*SUM(curv(:)* &
                              ((1.0_dp - z_weight)*MATMUL(raw, step) + &
                               z_weight*MATMUL(physical, step) + &
                               z_weight*MATMUL(density_raw, step) - actual(:))**2)

      IF (err > tol .OR. endpoint_err > tol .OR. zero_err > tol .OR. full_err > tol .OR. mixed_gap > tol .OR. &
          legacy_gap < 1.0E-6_dp) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,6(1X,ES12.5))') TRIM(label), &
            "failed with error/endpoint/zero/full/mixedgap/legacygap", &
            err, endpoint_err, zero_err, full_err, mixed_gap, legacy_gap
      ELSE
         WRITE (*, '(A,1X,A,6(1X,ES12.5))') TRIM(label), &
            "OK with error/endpoint/zero/full/mixedgap/legacygap", &
            err, endpoint_err, zero_err, full_err, mixed_gap, legacy_gap
      END IF

   END SUBROUTINE check_projector_mixed_direct_response_blend

! **************************************************************************************************
!> \brief Check that the floor-density Direct blend is Schur/RHS consistent.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_direct_floor_density_schur(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: ne = 4, nrot = 3

      INTEGER                                            :: r, s
      REAL(KIND=dp)                                      :: base_self, direct_self, err, legacy_gap, &
                                                            rhs_err, w
      REAL(KIND=dp), DIMENSION(ne)                       :: grad_e
      REAL(KIND=dp), DIMENSION(ne, nrot)                 :: b_density, b_direct, b_floor, b_phys, &
                                                            r_density, r_direct, r_floor, r_phys
      REAL(KIND=dp), DIMENSION(nrot)                     :: grad_rot, rhs_direct, rhs_expected
      REAL(KIND=dp), DIMENSION(nrot, nrot)               :: rot_metric, schur_direct, &
                                                            schur_expected, schur_legacy

      w = 0.35_dp
      grad_e(:) = [0.11_dp, -0.07_dp, 0.05_dp, -0.03_dp]
      grad_rot(:) = [0.21_dp, -0.18_dp, 0.09_dp]
      rot_metric(:, :) = RESHAPE([ &
                                 0.91_dp, 0.05_dp, -0.04_dp, &
                                 0.05_dp, 0.76_dp, 0.03_dp, &
                                 -0.04_dp, 0.03_dp, 0.69_dp], [nrot, nrot])
      b_phys(:, :) = RESHAPE([ &
                             0.13_dp, -0.04_dp, 0.02_dp, -0.01_dp, &
                             -0.08_dp, 0.12_dp, -0.03_dp, 0.05_dp, &
                             0.04_dp, -0.06_dp, 0.09_dp, -0.02_dp], [ne, nrot])
      r_phys(:, :) = RESHAPE([ &
                             0.17_dp, -0.05_dp, 0.03_dp, -0.02_dp, &
                             -0.06_dp, 0.10_dp, -0.04_dp, 0.07_dp, &
                             0.02_dp, -0.08_dp, 0.11_dp, -0.01_dp], [ne, nrot])
      b_density(:, :) = RESHAPE([ &
                                -0.03_dp, 0.07_dp, -0.01_dp, 0.04_dp, &
                                0.05_dp, -0.02_dp, 0.08_dp, -0.03_dp, &
                                -0.04_dp, 0.03_dp, -0.07_dp, 0.02_dp], [ne, nrot])
      r_density(:, :) = RESHAPE([ &
                                -0.02_dp, 0.05_dp, -0.03_dp, 0.06_dp, &
                                0.07_dp, -0.01_dp, 0.04_dp, -0.05_dp, &
                                -0.01_dp, 0.06_dp, -0.08_dp, 0.03_dp], [ne, nrot])
      b_floor(:, :) = RESHAPE([ &
                              -0.01_dp, 0.03_dp, 0.02_dp, -0.02_dp, &
                              0.04_dp, -0.04_dp, 0.05_dp, -0.01_dp, &
                              -0.03_dp, 0.02_dp, -0.04_dp, 0.04_dp], [ne, nrot])
      r_floor(:, :) = RESHAPE([ &
                              0.01_dp, 0.02_dp, -0.01_dp, -0.03_dp, &
                              0.03_dp, -0.05_dp, 0.02_dp, 0.01_dp, &
                              -0.02_dp, 0.01_dp, -0.03_dp, 0.05_dp], [ne, nrot])

      b_direct(:, :) = (1.0_dp - w)*b_density(:, :) + w*b_floor(:, :)
      r_direct(:, :) = (1.0_dp - w)*r_density(:, :) + w*r_floor(:, :)

      DO r = 1, nrot
         rhs_direct(r) = grad_rot(r) + DOT_PRODUCT(r_phys(:, r), grad_e) + &
                         DOT_PRODUCT(r_direct(:, r), grad_e)
         rhs_expected(r) = grad_rot(r) + DOT_PRODUCT(r_phys(:, r) + r_direct(:, r), grad_e)
      END DO
      schur_direct(:, :) = rot_metric(:, :)
      schur_legacy(:, :) = rot_metric(:, :)
      schur_expected(:, :) = rot_metric(:, :)
      DO s = 1, nrot
         DO r = 1, nrot
            base_self = 0.5_dp* &
                        (DOT_PRODUCT(b_density(:, r), r_density(:, s)) + &
                         DOT_PRODUCT(b_density(:, s), r_density(:, r)))
            direct_self = 0.5_dp* &
                          (DOT_PRODUCT(b_direct(:, r), r_direct(:, s)) + &
                           DOT_PRODUCT(b_direct(:, s), r_direct(:, r)))
            schur_direct(r, s) = schur_direct(r, s) + direct_self - base_self
            schur_direct(r, s) = schur_direct(r, s) - 0.5_dp* &
                                 (DOT_PRODUCT(b_phys(:, r), r_phys(:, s)) + DOT_PRODUCT(b_phys(:, s), r_phys(:, r)))
            schur_direct(r, s) = schur_direct(r, s) - 0.5_dp* &
                                 (DOT_PRODUCT(b_phys(:, r), r_direct(:, s)) + DOT_PRODUCT(b_direct(:, r), r_phys(:, s)) + &
                                  DOT_PRODUCT(b_direct(:, r), r_direct(:, s)) + DOT_PRODUCT(b_phys(:, s), r_direct(:, r)) + &
                                  DOT_PRODUCT(b_direct(:, s), r_phys(:, r)) + DOT_PRODUCT(b_direct(:, s), r_direct(:, r)))
            schur_legacy(r, s) = schur_legacy(r, s) - 0.5_dp* &
                                 (DOT_PRODUCT(b_phys(:, r), r_phys(:, s)) + DOT_PRODUCT(b_phys(:, s), r_phys(:, r)))
            schur_legacy(r, s) = schur_legacy(r, s) - 0.5_dp* &
                                 (DOT_PRODUCT(b_phys(:, r), r_direct(:, s)) + DOT_PRODUCT(b_direct(:, r), r_phys(:, s)) + &
                                  DOT_PRODUCT(b_direct(:, r), r_direct(:, s)) + DOT_PRODUCT(b_phys(:, s), r_direct(:, r)) + &
                                  DOT_PRODUCT(b_direct(:, s), r_phys(:, r)) + DOT_PRODUCT(b_direct(:, s), r_direct(:, r)))
            schur_expected(r, s) = schur_expected(r, s) + direct_self - base_self - 0.5_dp* &
                                   (DOT_PRODUCT(b_phys(:, r) + b_direct(:, r), r_phys(:, s) + r_direct(:, s)) + &
                                    DOT_PRODUCT(b_phys(:, s) + b_direct(:, s), r_phys(:, r) + r_direct(:, r)))
         END DO
      END DO

      err = MAXVAL(ABS(schur_direct(:, :) - schur_expected(:, :)))
      rhs_err = MAXVAL(ABS(rhs_direct(:) - rhs_expected(:)))
      legacy_gap = MAXVAL(ABS(schur_legacy(:, :) - schur_expected(:, :)))
      IF (err > tol .OR. rhs_err > tol .OR. legacy_gap <= 1.0E-6_dp) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,3(1X,ES12.5))') TRIM(label), &
            "failed with schur/rhs/legacy", err, rhs_err, legacy_gap
      ELSE
         WRITE (*, '(A,1X,A,3(1X,ES12.5))') TRIM(label), &
            "OK with schur/rhs/legacy", err, rhs_err, legacy_gap
      END IF

   END SUBROUTINE check_projector_direct_floor_density_schur

! **************************************************************************************************
!> \brief Check the sign-aware finite Projector-Hessian response model in mixed curvature modes.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_modal_response_sign_model(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: ne = 3, nrot = 4

      INTEGER                                            :: i, imode
      REAL(KIND=dp) :: active, applied_weight_dead, applied_weight_diag, applied_weight_ref, &
         applied_weight_rot, base_weight, c, cross_drop_dead, cross_drop_diag, cross_drop_rot, &
         dead_active, drop_dead, drop_diag, drop_ref, drop_rot, err, eval_abs_avg, mode_norm, &
         mode_weight, neg_balance, neg_coeff_dead, neg_coeff_diag, neg_coeff_rot, neg_drop_dead, &
         neg_drop_diag, neg_drop_rot, neg_fraction_dead, neg_fraction_diag, neg_fraction_ref, &
         neg_fraction_rot, neg_norm_diag, neg_weight_sum, pos_coeff_dead, pos_coeff_diag, &
         pos_coeff_rot, pos_drop_dead, pos_drop_diag, pos_drop_rot, pos_weight_sum
      REAL(KIND=dp) :: raw_neg_fraction_dead, raw_neg_fraction_diag, raw_neg_fraction_ref, &
         raw_neg_fraction_rot, raw_neg_sum, response_total_sum, s, spectral_weight, &
         step_weight_sum, target_neg, target_pos
      REAL(KIND=dp), DIMENSION(ne) :: cand_dead, cand_diag, cand_ref, cand_rot, corr_ref, curv, &
         energy_base, grad_e, mode_rx, response_base, response_neg_ref, response_pos_ref
      REAL(KIND=dp), DIMENSION(ne, nrot)                 :: response_delta_diag, response_delta_rot, &
                                                            response_exact_diag, &
                                                            response_exact_rot, response_raw
      REAL(KIND=dp), DIMENSION(nrot)                     :: eval_diag, step_diag, step_rot
      REAL(KIND=dp), DIMENSION(nrot, nrot)               :: delta_diag, delta_rot, q, rot_hessian

      base_weight = 0.42_dp
      target_pos = 0.63_dp
      target_neg = 0.37_dp
      curv(:) = [0.52_dp, 0.81_dp, 0.34_dp]
      response_base(:) = [0.08_dp, -0.04_dp, 0.03_dp]
      energy_base(:) = [0.03_dp, -0.02_dp, 0.01_dp]
      eval_diag(:) = [0.20_dp, -0.32_dp, 0.07_dp, -0.18_dp]
      step_diag(:) = [0.40_dp, -0.25_dp, 0.30_dp, 0.18_dp]

      response_delta_diag(:, 1) = [0.20_dp, -0.10_dp, 0.06_dp]
      response_delta_diag(:, 2) = [-0.12_dp, 0.04_dp, 0.15_dp]
      response_delta_diag(:, 3) = [0.05_dp, 0.11_dp, -0.07_dp]
      response_delta_diag(:, 4) = [0.09_dp, -0.08_dp, 0.04_dp]
      response_raw(:, :) = 0.0_dp
      response_exact_diag(:, :) = response_delta_diag(:, :)

      response_pos_ref(:) = 0.0_dp
      response_neg_ref(:) = 0.0_dp
      response_total_sum = 0.0_dp
      raw_neg_sum = 0.0_dp
      eval_abs_avg = SUM(ABS(eval_diag(:)))/REAL(nrot, KIND=dp)
      DO imode = 1, nrot
         mode_rx(:) = step_diag(imode)*response_delta_diag(:, imode)
         mode_norm = DOT_PRODUCT(curv(:)*mode_rx(:), mode_rx(:))
         response_total_sum = response_total_sum + mode_norm
         IF (eval_diag(imode) < 0.0_dp) raw_neg_sum = raw_neg_sum + mode_norm
      END DO
      raw_neg_fraction_ref = raw_neg_sum/response_total_sum
      neg_balance = dense_projector_negative_response_balance(response_total_sum - raw_neg_sum, raw_neg_sum)

      step_weight_sum = DOT_PRODUCT(step_diag(:), step_diag(:))
      pos_weight_sum = 0.0_dp
      neg_weight_sum = 0.0_dp
      DO imode = 1, nrot
         mode_rx(:) = step_diag(imode)*response_delta_diag(:, imode)
         spectral_weight = 1.0_dp
         IF (eval_diag(imode) < 0.0_dp) THEN
            spectral_weight = 1.0_dp/MAX(1.0_dp, ABS(eval_diag(imode))/eval_abs_avg)
            spectral_weight = spectral_weight*neg_balance
         END IF
         mode_weight = base_weight*MAX(0.0_dp, MIN(1.0_dp, spectral_weight))
         IF (eval_diag(imode) < 0.0_dp) THEN
            response_neg_ref(:) = response_neg_ref(:) + mode_weight*mode_rx(:)
            neg_weight_sum = neg_weight_sum + mode_weight*step_diag(imode)*step_diag(imode)
         ELSE
            response_pos_ref(:) = response_pos_ref(:) + mode_weight*mode_rx(:)
            pos_weight_sum = pos_weight_sum + mode_weight*step_diag(imode)*step_diag(imode)
         END IF
      END DO

      corr_ref(:) = target_pos*response_pos_ref(:) + target_neg*response_neg_ref(:)
      cand_ref(:) = response_base(:) + corr_ref(:)
      grad_e(:) = curv(:)*energy_base(:) + curv(:)*corr_ref(:)
      drop_ref = 0.5_dp*DOT_PRODUCT(curv(:)*corr_ref(:), corr_ref(:))
      neg_norm_diag = DOT_PRODUCT(curv(:)*(target_neg*response_neg_ref(:)), &
                                  target_neg*response_neg_ref(:))
      neg_fraction_ref = neg_norm_diag/(2.0_dp*drop_ref)
      applied_weight_ref = (target_pos*pos_weight_sum + target_neg*neg_weight_sum)/step_weight_sum

      rot_hessian(:, :) = 0.0_dp
      delta_diag(:, :) = 0.0_dp
      DO i = 1, nrot
         delta_diag(i, i) = eval_diag(i)
      END DO
      CALL dense_projector_modal_response_candidate(rot_hessian, delta_diag, step_diag, response_raw, &
                                                    response_exact_diag, curv, grad_e, energy_base, &
                                                    response_base, base_weight, cand_diag, applied_weight_diag, &
                                                    drop_diag, neg_fraction_diag, raw_neg_fraction_diag, &
                                                    pos_drop_diag, neg_drop_diag, cross_drop_diag, &
                                                    pos_coeff_diag, neg_coeff_diag)

      q(:, :) = 0.0_dp
      DO i = 1, nrot
         q(i, i) = 1.0_dp
      END DO
      c = COS(0.37_dp)
      s = SIN(0.37_dp)
      q(1, 1) = c
      q(2, 2) = c
      q(1, 2) = -s
      q(2, 1) = s
      c = COS(-0.28_dp)
      s = SIN(-0.28_dp)
      q(3, 3) = c
      q(4, 4) = c
      q(3, 4) = -s
      q(4, 3) = s
      delta_rot(:, :) = MATMUL(q, MATMUL(delta_diag, TRANSPOSE(q)))
      step_rot(:) = MATMUL(q, step_diag)
      response_delta_rot(:, :) = MATMUL(response_delta_diag, TRANSPOSE(q))
      response_exact_rot(:, :) = response_delta_rot(:, :)
      CALL dense_projector_modal_response_candidate(rot_hessian, delta_rot, step_rot, response_raw, &
                                                    response_exact_rot, curv, grad_e, energy_base, &
                                                    response_base, base_weight, cand_rot, applied_weight_rot, &
                                                    drop_rot, neg_fraction_rot, raw_neg_fraction_rot, &
                                                    pos_drop_rot, neg_drop_rot, cross_drop_rot, &
                                                    pos_coeff_rot, neg_coeff_rot)
      CALL dense_projector_modal_response_candidate(rot_hessian, delta_diag, step_diag, response_exact_diag, &
                                                    response_exact_diag, curv, grad_e, energy_base, &
                                                    response_base, base_weight, cand_dead, applied_weight_dead, &
                                                    drop_dead, neg_fraction_dead, raw_neg_fraction_dead, &
                                                    pos_drop_dead, neg_drop_dead, cross_drop_dead, &
                                                    pos_coeff_dead, neg_coeff_dead)
      dead_active = MAXVAL(ABS(cand_dead(:) - response_base(:)))

      err = MAXVAL(ABS(cand_diag(:) - cand_ref(:)))
      err = MAX(err, MAXVAL(ABS(cand_rot(:) - cand_ref(:))))
      err = MAX(err, dead_active, drop_dead, raw_neg_fraction_dead, &
                neg_fraction_dead, pos_drop_dead, neg_drop_dead, cross_drop_dead)
      err = MAX(err, ABS(drop_diag - drop_ref), ABS(drop_rot - drop_ref))
      err = MAX(err, ABS(drop_diag - (pos_drop_diag + neg_drop_diag + cross_drop_diag)), &
                ABS(drop_rot - (pos_drop_rot + neg_drop_rot + cross_drop_rot)))
      err = MAX(err, ABS(raw_neg_fraction_diag - raw_neg_fraction_ref), &
                ABS(raw_neg_fraction_rot - raw_neg_fraction_ref))
      err = MAX(err, ABS(neg_fraction_diag - neg_fraction_ref), &
                ABS(neg_fraction_rot - neg_fraction_ref))
      err = MAX(err, ABS(applied_weight_diag - applied_weight_ref), &
                ABS(applied_weight_rot - applied_weight_ref))
      err = MAX(err, ABS(pos_coeff_diag - target_pos), ABS(pos_coeff_rot - target_pos))
      err = MAX(err, ABS(neg_coeff_diag - target_neg), ABS(neg_coeff_rot - target_neg))

      active = MAXVAL(ABS(cand_ref(:) - response_base(:)))
      IF (err > tol .OR. active < 1.0E-4_dp .OR. raw_neg_fraction_ref < 0.1_dp .OR. &
          raw_neg_fraction_ref > 0.9_dp .OR. neg_fraction_ref <= 0.0_dp) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,5(1X,ES12.5))') TRIM(label), &
            "failed with error/active/dead/rawneg/corrneg", &
            err, active, dead_active, raw_neg_fraction_ref, neg_fraction_ref
      ELSE
         WRITE (*, '(A,1X,A,5(1X,ES12.5))') TRIM(label), &
            "OK with error/active/dead/rawneg/corrneg", err, active, dead_active, &
            raw_neg_fraction_ref, neg_fraction_ref
      END IF

   END SUBROUTINE check_projector_modal_response_sign_model

! **************************************************************************************************
!> \brief Negative-only Projector modes must not be suppressed to zero.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_modal_response_negative_floor(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: ne = 2, nrot = 3

      INTEGER                                            :: i, imode
      REAL(KIND=dp) :: active, applied_weight, applied_weight_ref, base_weight, correction_drop, &
         correction_drop_ref, cross_drop, err, eval_abs_avg, mode_weight, negative_balance_ref, &
         negative_coeff, negative_drop, negative_fraction, negative_fraction_ref, &
         negative_weight_sum, positive_coeff, positive_drop, raw_negative, raw_negative_sum, &
         spectral_weight, step_weight_sum
      REAL(KIND=dp), DIMENSION(ne)                       :: candidate, candidate_ref, curv, &
                                                            energy_base, grad_e, mode_rx, &
                                                            response_base, response_neg_ref
      REAL(KIND=dp), DIMENSION(ne, nrot)                 :: response_delta, response_raw
      REAL(KIND=dp), DIMENSION(nrot)                     :: evals, step_rot
      REAL(KIND=dp), DIMENSION(nrot, nrot)               :: delta, rot_hessian

      base_weight = 0.35_dp
      curv(:) = [0.60_dp, 0.90_dp]
      energy_base(:) = [0.02_dp, -0.01_dp]
      response_base(:) = [0.01_dp, -0.03_dp]
      evals(:) = [-0.40_dp, -0.18_dp, -0.07_dp]
      step_rot(:) = [0.32_dp, -0.24_dp, 0.19_dp]
      response_delta(:, 1) = [0.15_dp, -0.04_dp]
      response_delta(:, 2) = [-0.08_dp, 0.11_dp]
      response_delta(:, 3) = [0.05_dp, 0.07_dp]
      response_raw(:, :) = 0.0_dp

      response_neg_ref(:) = 0.0_dp
      eval_abs_avg = SUM(ABS(evals(:)))/REAL(nrot, KIND=dp)
      raw_negative_sum = 0.0_dp
      DO imode = 1, nrot
         mode_rx(:) = step_rot(imode)*response_delta(:, imode)
         raw_negative_sum = raw_negative_sum + MAX(0.0_dp, DOT_PRODUCT(curv(:)*mode_rx(:), mode_rx(:)))
      END DO
      negative_balance_ref = dense_projector_negative_response_balance(0.0_dp, raw_negative_sum)
      step_weight_sum = DOT_PRODUCT(step_rot(:), step_rot(:))
      negative_weight_sum = 0.0_dp
      DO imode = 1, nrot
         spectral_weight = 1.0_dp/MAX(1.0_dp, ABS(evals(imode))/eval_abs_avg)
         spectral_weight = spectral_weight*negative_balance_ref
         mode_weight = base_weight*MAX(0.0_dp, MIN(1.0_dp, spectral_weight))
         mode_rx(:) = step_rot(imode)*response_delta(:, imode)
         response_neg_ref(:) = response_neg_ref(:) + mode_weight*mode_rx(:)
         negative_weight_sum = negative_weight_sum + mode_weight*step_rot(imode)*step_rot(imode)
      END DO
      candidate_ref(:) = response_base(:) + response_neg_ref(:)
      grad_e(:) = curv(:)*energy_base(:) + curv(:)*response_neg_ref(:)
      correction_drop_ref = 0.5_dp*DOT_PRODUCT(curv(:)*response_neg_ref(:), response_neg_ref(:))
      negative_fraction_ref = 1.0_dp
      applied_weight_ref = negative_weight_sum/step_weight_sum

      rot_hessian(:, :) = 0.0_dp
      delta(:, :) = 0.0_dp
      DO i = 1, nrot
         delta(i, i) = evals(i)
      END DO
      CALL dense_projector_modal_response_candidate(rot_hessian, delta, step_rot, response_raw, &
                                                    response_delta, curv, grad_e, energy_base, response_base, &
                                                    base_weight, candidate, applied_weight, correction_drop, &
                                                    negative_fraction, raw_negative, positive_drop, negative_drop, &
                                                    cross_drop, positive_coeff, negative_coeff)

      active = MAXVAL(ABS(candidate(:) - response_base(:)))
      err = MAXVAL(ABS(candidate(:) - candidate_ref(:)))
      err = MAX(err, ABS(applied_weight - applied_weight_ref))
      err = MAX(err, ABS(correction_drop - correction_drop_ref))
      err = MAX(err, ABS(negative_fraction - negative_fraction_ref))
      err = MAX(err, ABS(raw_negative - 1.0_dp), ABS(positive_drop), ABS(cross_drop), &
                ABS(positive_coeff - 1.0_dp), ABS(negative_coeff - 1.0_dp))

      IF (err > tol .OR. active < 1.0E-5_dp .OR. negative_drop <= 0.0_dp) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "failed with error/active/drop/rawneg", err, active, negative_drop, raw_negative
      ELSE
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "OK with error/active/drop/rawneg", err, active, negative_drop, raw_negative
      END IF

   END SUBROUTINE check_projector_modal_response_negative_floor

! **************************************************************************************************
!> \brief Negative-step history may lift the inactive Projector shadow only coherently.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_negative_step_balance(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      REAL(KIND=dp), PARAMETER :: response_ratio_gate = 2.0E-2_dp

      REAL(KIND=dp)                                      :: err, full_negative, high_ratio, &
                                                            neutral_step, noisy_history, &
                                                            partial_negative, positive_step

      full_negative = dense_projector_negative_step_balance(1.0_dp, 1.0_dp, 0.0_dp, &
                                                            1.0_dp, response_ratio_gate)
      partial_negative = dense_projector_negative_step_balance(1.0_dp, 1.0_dp, 0.0_dp, &
                                                               0.56_dp, response_ratio_gate)
      neutral_step = dense_projector_negative_step_balance(1.0_dp, 1.0_dp, 0.0_dp, &
                                                           0.50_dp, response_ratio_gate)
      positive_step = dense_projector_negative_step_balance(1.0_dp, 1.0_dp, 0.0_dp, &
                                                            0.25_dp, response_ratio_gate)
      high_ratio = dense_projector_negative_step_balance(1.0_dp, 1.0_dp, &
                                                         16.0_dp*response_ratio_gate, &
                                                         1.0_dp, response_ratio_gate)
      noisy_history = dense_projector_negative_step_balance(0.25_dp, 0.16_dp, &
                                                            4.0_dp*response_ratio_gate, &
                                                            1.0_dp, response_ratio_gate)

      err = ABS(full_negative - projector_response_negative_balance_max)
      err = MAX(err, ABS(partial_negative - SQRT(0.12_dp)))
      err = MAX(err, ABS(neutral_step - projector_response_negative_balance_floor))
      err = MAX(err, ABS(positive_step - projector_response_negative_balance_floor))
      err = MAX(err, ABS(high_ratio - 0.25_dp))
      err = MAX(err, ABS(noisy_history - 0.04_dp))

      IF (err > tol .OR. full_negative <= neutral_step .OR. partial_negative <= neutral_step .OR. &
          high_ratio >= full_negative .OR. noisy_history >= neutral_step) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,6(1X,ES12.5))') TRIM(label), &
            "failed with error/full/partial/neutral/highratio/noisy", &
            err, full_negative, partial_negative, neutral_step, high_ratio, noisy_history
      ELSE
         WRITE (*, '(A,1X,A,6(1X,ES12.5))') TRIM(label), &
            "OK with error/full/partial/neutral/highratio/noisy", &
            err, full_negative, partial_negative, neutral_step, high_ratio, noisy_history
      END IF

   END SUBROUTINE check_projector_negative_step_balance

! **************************************************************************************************
!> \brief Check that residual-tail gauge spikes are not hidden by a low EMA.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_tail_gauge_quality(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      REAL(KIND=dp)                                      :: early, err, good_residual, neutral, &
                                                            persistent, plateau_adjusted, &
                                                            plateau_raw, spike

      early = dense_schur_step_gauge_quality(38.0_dp, 0.66_dp, 0.0_dp, 0.0_dp, 3)
      spike = dense_schur_step_gauge_quality(38.0_dp, 0.66_dp, 0.0_dp, 0.005_dp, 52)
      persistent = dense_schur_step_gauge_quality(2.0_dp, 2.2_dp, 0.0_dp, 0.0_dp, 52)
      neutral = dense_schur_step_gauge_quality(0.8_dp, 0.66_dp, 0.0_dp, 0.005_dp, 52)
      good_residual = dense_schur_step_gauge_quality(38.0_dp, 0.66_dp, 0.30_dp, 0.28_dp, 52)
      plateau_adjusted = dense_schur_step_gauge_quality(4.44_dp, 3.15_dp, 0.25_dp, 0.152_dp, 104)
      plateau_raw = dense_schur_step_gauge_quality( &
                    4.44_dp, 3.15_dp, &
                    dense_step_gauge_residual_quality(0.0_dp, 0.25_dp), &
                    dense_step_gauge_residual_quality(0.0_dp, 0.152_dp), 104)

      err = ABS(early - 1.0_dp)
      err = MAX(err, ABS(spike - 0.35_dp))
      err = MAX(err, ABS(persistent - 0.93_dp))
      err = MAX(err, ABS(neutral - 1.0_dp))
      err = MAX(err, ABS(good_residual - 1.0_dp))
      err = MAX(err, ABS(plateau_raw - 0.706_dp))

      IF (err > tol .OR. spike >= persistent .OR. spike >= neutral .OR. &
          persistent >= neutral .OR. good_residual < 1.0_dp .OR. &
          plateau_raw >= plateau_adjusted .OR. plateau_raw > 0.75_dp) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,7(1X,ES12.5))') TRIM(label), &
            "failed with early/spike/persistent/neutral/good/padj/praw", early, spike, persistent, &
            neutral, good_residual, plateau_adjusted, plateau_raw
      ELSE
         WRITE (*, '(A,1X,A,7(1X,ES12.5))') TRIM(label), &
            "OK with early/spike/persistent/neutral/good/padj/praw", early, spike, persistent, &
            neutral, good_residual, plateau_adjusted, plateau_raw
      END IF

   END SUBROUTINE check_projector_tail_gauge_quality

! **************************************************************************************************
!> \brief Check residual-tail absolute trust for offdiag projector response.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_offdiag_tail_abs_trust(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      REAL(KIND=dp), PARAMETER                           :: response_abs_trust = 8.0E-3_dp, &
                                                            sigma = 1.0E-3_dp, tail_delta = 1.0E-4_dp

      REAL(KIND=dp)                                      :: capped_by_existing, err, good_history, &
                                                            late_tail, outside_tail, weak_tail

      outside_tail = dense_offdiag_tail_abs_trust(response_abs_trust, 2.0_dp*tail_delta, sigma, &
                                                  0.70_dp, tail_delta)
      good_history = dense_offdiag_tail_abs_trust(response_abs_trust, 5.0E-5_dp, sigma, &
                                                  0.98_dp, tail_delta)
      weak_tail = dense_offdiag_tail_abs_trust(response_abs_trust, 5.0E-5_dp, sigma, &
                                               0.70_dp, tail_delta)
      late_tail = dense_offdiag_tail_abs_trust(response_abs_trust, 1.0E-6_dp, sigma, &
                                               0.70_dp, tail_delta)
      capped_by_existing = dense_offdiag_tail_abs_trust(1.0E-5_dp, 5.0E-5_dp, sigma, &
                                                        0.70_dp, tail_delta)

      err = ABS(outside_tail - response_abs_trust)
      err = MAX(err, ABS(good_history - response_abs_trust))
      err = MAX(err, ABS(weak_tail - SQRT(5.0E-8_dp)))
      err = MAX(err, ABS(late_tail - 0.5_dp*SQRT(1.0E-7_dp)))
      err = MAX(err, ABS(capped_by_existing - 1.0E-5_dp))

      IF (err > tol .OR. weak_tail >= response_abs_trust .OR. late_tail > weak_tail) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,6(1X,ES12.5))') TRIM(label), &
            "failed err/outside/good/weak/late/capped", err, outside_tail, good_history, &
            weak_tail, late_tail, capped_by_existing
      ELSE
         WRITE (*, '(A,1X,A,6(1X,ES12.5))') TRIM(label), &
            "OK err/outside/good/weak/late/capped", err, outside_tail, good_history, &
            weak_tail, late_tail, capped_by_existing
      END IF

   END SUBROUTINE check_projector_offdiag_tail_abs_trust

! **************************************************************************************************
!> \brief Raw Schur fallback must keep gauge pair floors in the metallic residual tail.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_raw_tail_rotation_floor(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      REAL(KIND=dp)                                      :: above_tail, base_tail, err, &
                                                            frontier_tail, gauge_tail

      above_tail = dense_raw_tail_rotation_floor(1.0_dp, 1.5_dp, 2.0E-4_dp, 1.0E-4_dp)
      base_tail = dense_raw_tail_rotation_floor(1.0_dp, 1.0_dp, 1.0E-6_dp, 1.0E-4_dp)
      frontier_tail = dense_raw_tail_rotation_floor(1.0_dp, 0.5_dp, 1.0E-6_dp, 1.0E-4_dp)
      gauge_tail = dense_raw_tail_rotation_floor(1.0_dp, 1.5_dp, 1.0E-6_dp, 1.0E-4_dp)

      err = ABS(above_tail - 1.0_dp)
      err = MAX(err, ABS(base_tail - 1.0_dp))
      err = MAX(err, ABS(frontier_tail - 1.0_dp))
      err = MAX(err, ABS(gauge_tail - 1.5_dp))

      IF (err > tol .OR. gauge_tail <= base_tail .OR. frontier_tail /= base_tail) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "failed with above/base/frontier/gauge", above_tail, base_tail, &
            frontier_tail, gauge_tail
      ELSE
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "OK with above/base/frontier/gauge", above_tail, base_tail, &
            frontier_tail, gauge_tail
      END IF

   END SUBROUTINE check_projector_raw_tail_rotation_floor

! **************************************************************************************************
!> \brief Accepted density-response history should harden the active-window Direct blend only in tail.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_projector_direct_density_curvature_mix(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      REAL(KIND=dp)                                      :: above_tail, capped_tail, err, relative_tail, tail

      above_tail = dense_physical_direct_density_mix(0.20_dp, 0.60_dp, 2.0E-4_dp, &
                                                     1.0E-4_dp, 0.95_dp)
      tail = dense_physical_direct_density_mix(0.20_dp, 0.60_dp, 1.0E-6_dp, &
                                               1.0E-4_dp, 0.95_dp)
      capped_tail = dense_physical_direct_density_mix(-0.10_dp, 2.0_dp, 1.0E-6_dp, &
                                                      1.0E-4_dp, 0.95_dp)
      relative_tail = dense_physical_direct_density_mix(0.02_dp, 0.025_dp, 1.0E-6_dp, &
                                                        1.0E-4_dp, 0.95_dp, projector_weight=0.05_dp)

      err = ABS(above_tail - 0.20_dp)
      err = MAX(err, ABS(tail - 0.60_dp))
      err = MAX(err, ABS(capped_tail - 0.95_dp))
      err = MAX(err, ABS(relative_tail - 0.50_dp))

      IF (err > tol .OR. tail <= above_tail .OR. relative_tail <= tail*0.5_dp) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "failed with above/tail/capped/relative", above_tail, tail, capped_tail, relative_tail
      ELSE
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "OK with above/tail/capped/relative", above_tail, tail, capped_tail, relative_tail
      END IF

   END SUBROUTINE check_projector_direct_density_curvature_mix

! **************************************************************************************************
!> \brief Dense copy of the residual-tail gauge-step quality.
!> \param gauge_ratio ...
!> \param gauge_ratio_ema ...
!> \param residual_quality ...
!> \param residual_quality_ema ...
!> \param samples ...
!> \return quality in [0,1]
! **************************************************************************************************
   FUNCTION dense_schur_step_gauge_quality(gauge_ratio, gauge_ratio_ema, residual_quality, &
                                           residual_quality_ema, samples) RESULT(quality)
      REAL(KIND=dp), INTENT(IN)                          :: gauge_ratio, gauge_ratio_ema, &
                                                            residual_quality, residual_quality_ema
      INTEGER, INTENT(IN)                                :: samples
      REAL(KIND=dp)                                      :: quality

      REAL(KIND=dp)                                      :: gauge_signal, gauge_weight, &
                                                            residual_signal, residual_weight

      quality = 1.0_dp
      IF (samples < 4) RETURN
      gauge_signal = MAX(0.0_dp, gauge_ratio, gauge_ratio_ema)
      residual_signal = MIN(MAX(0.0_dp, residual_quality), MAX(0.0_dp, residual_quality_ema))
      IF (gauge_signal <= 1.5_dp .OR. residual_signal >= 0.25_dp) RETURN

      gauge_weight = MIN(1.0_dp, MAX(0.0_dp, (gauge_signal - 1.5_dp)/6.5_dp))
      residual_weight = MIN(1.0_dp, MAX(0.0_dp, (0.25_dp - residual_signal)/0.25_dp))
      quality = 1.0_dp - 0.65_dp*gauge_weight*residual_weight
      quality = MIN(1.0_dp, MAX(0.25_dp, quality))

   END FUNCTION dense_schur_step_gauge_quality

! **************************************************************************************************
!> \brief Dense copy of absolute residual-tail trust for offdiag projector response.
!> \param response_abs_trust ...
!> \param residual_delta ...
!> \param sigma ...
!> \param step_gauge_quality ...
!> \param tail_delta ...
!> \return absolute offdiag response trust
! **************************************************************************************************
   FUNCTION dense_offdiag_tail_abs_trust(response_abs_trust, residual_delta, sigma, &
                                         step_gauge_quality, tail_delta) RESULT(abs_trust)
      REAL(KIND=dp), INTENT(IN)                          :: response_abs_trust, residual_delta, sigma, &
                                                            step_gauge_quality, tail_delta
      REAL(KIND=dp)                                      :: abs_trust

      REAL(KIND=dp)                                      :: sigma_eff, tail_cap, tail_floor

      abs_trust = response_abs_trust
      IF (response_abs_trust <= EPSILON(response_abs_trust)) RETURN
      IF (residual_delta > tail_delta) RETURN
      IF (step_gauge_quality >= 0.95_dp) RETURN

      sigma_eff = MAX(sigma, EPSILON(1.0_dp))
      tail_cap = SQRT(MAX(0.0_dp, residual_delta)*sigma_eff)
      tail_floor = 0.5_dp*SQRT(tail_delta*sigma_eff)
      tail_cap = MAX(tail_cap, tail_floor)
      abs_trust = MIN(response_abs_trust, tail_cap)

   END FUNCTION dense_offdiag_tail_abs_trust

! **************************************************************************************************
!> \brief Dense copy of the strict residual quality used for gauge-like tail rotations.
!> \param raw_residual_quality ...
!> \param residual_quality ...
!> \return quality in [0,+inf]
! **************************************************************************************************
   FUNCTION dense_step_gauge_residual_quality(raw_residual_quality, residual_quality) RESULT(gauge_quality)
      REAL(KIND=dp), INTENT(IN)                          :: raw_residual_quality, residual_quality
      REAL(KIND=dp)                                      :: gauge_quality

      gauge_quality = MIN(MAX(0.0_dp, raw_residual_quality), &
                          MAX(0.0_dp, residual_quality))

   END FUNCTION dense_step_gauge_residual_quality

! **************************************************************************************************
!> \brief Dense copy of the raw-fallback residual-tail rotation floor.
!> \param base_floor ...
!> \param pair_floor ...
!> \param residual_delta ...
!> \param tail_delta ...
!> \return raw-fallback row floor
! **************************************************************************************************
   FUNCTION dense_raw_tail_rotation_floor(base_floor, pair_floor, residual_delta, tail_delta) &
      RESULT(row_floor)
      REAL(KIND=dp), INTENT(IN)                          :: base_floor, pair_floor, residual_delta, &
                                                            tail_delta
      REAL(KIND=dp)                                      :: row_floor

      row_floor = base_floor
      IF (residual_delta <= tail_delta) row_floor = MAX(row_floor, pair_floor)

   END FUNCTION dense_raw_tail_rotation_floor

! **************************************************************************************************
!> \brief Dense mirror of the physical Direct density-history curvature blend.
!> \param active_weight ...
!> \param density_history_weight ...
!> \param residual_delta ...
!> \param tail_delta ...
!> \param max_weight ...
!> \param projector_weight ...
!> \return trace-active curvature mix weight
! **************************************************************************************************
   FUNCTION dense_physical_direct_density_mix(active_weight, density_history_weight, residual_delta, &
                                              tail_delta, max_weight, projector_weight) RESULT(mix_weight)
      REAL(KIND=dp), INTENT(IN)                          :: active_weight, density_history_weight, &
                                                            residual_delta, tail_delta, max_weight
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: projector_weight
      REAL(KIND=dp)                                      :: mix_weight

      REAL(KIND=dp)                                      :: density_weight

      mix_weight = MIN(max_weight, MAX(0.0_dp, active_weight))
      IF (residual_delta > tail_delta) RETURN

      density_weight = MAX(0.0_dp, density_history_weight)
      IF (PRESENT(projector_weight)) &
         density_weight = MAX(density_weight, dense_projector_tail_active_delta_weight(projector_weight, &
                                                                                      density_history_weight))
      mix_weight = MAX(mix_weight, MIN(max_weight, density_weight))

   END FUNCTION dense_physical_direct_density_mix

! **************************************************************************************************
!> \brief Dense copy of the sign-aware negative Projector-response mode balance.
!> \param positive_sum ...
!> \param raw_negative_sum ...
!> \return negative-mode balance in [0,1]
! **************************************************************************************************
   FUNCTION dense_projector_negative_response_balance(positive_sum, raw_negative_sum) RESULT(balance)
      REAL(KIND=dp), INTENT(IN)                          :: positive_sum, raw_negative_sum
      REAL(KIND=dp)                                      :: balance

      REAL(KIND=dp)                                      :: coherent_negative, modal_floor, &
                                                            negative_fraction, response_total, &
                                                            scalar_floor

      balance = 0.0_dp
      IF (raw_negative_sum <= EPSILON(raw_negative_sum)) RETURN

      balance = SQRT(MAX(0.0_dp, positive_sum)/raw_negative_sum)
      scalar_floor = projector_response_negative_balance_floor
      response_total = positive_sum + raw_negative_sum
      IF (response_total > EPSILON(response_total)) THEN
         negative_fraction = raw_negative_sum/response_total
         coherent_negative = MAX(0.0_dp, 2.0_dp*negative_fraction - 1.0_dp)
         modal_floor = SQRT(coherent_negative)
         scalar_floor = MAX(scalar_floor, MIN(projector_response_negative_balance_max, modal_floor))
      END IF
      balance = MAX(balance, scalar_floor)
      balance = MIN(1.0_dp, MAX(0.0_dp, balance))

   END FUNCTION dense_projector_negative_response_balance

! **************************************************************************************************
!> \brief Dense copy of the sign-aware negative Projector-step history balance.
!> \param history_quality ...
!> \param history_trust ...
!> \param history_response_ratio ...
!> \param step_negative_fraction ...
!> \param response_ratio_gate ...
!> \return negative-history balance in [0,1]
! **************************************************************************************************
   FUNCTION dense_projector_negative_step_balance(history_quality, history_trust, history_response_ratio, &
                                                  step_negative_fraction, response_ratio_gate) RESULT(balance)
      REAL(KIND=dp), INTENT(IN)                          :: history_quality, history_trust, &
                                                            history_response_ratio, &
                                                            step_negative_fraction, &
                                                            response_ratio_gate
      REAL(KIND=dp)                                      :: balance

      REAL(KIND=dp)                                      :: coherent_negative, history_scale, &
                                                            modal_floor, scalar_floor, &
                                                            step_negative

      history_scale = MIN(1.0_dp, MAX(0.0_dp, history_trust))/ &
                      SQRT(MAX(1.0_dp, history_response_ratio/response_ratio_gate))
      scalar_floor = projector_response_negative_balance_floor* &
                     MIN(1.0_dp, MAX(0.0_dp, history_quality))*history_scale
      step_negative = MIN(1.0_dp, MAX(0.0_dp, step_negative_fraction))
      coherent_negative = MAX(0.0_dp, 2.0_dp*step_negative - 1.0_dp)
      modal_floor = SQRT(coherent_negative)* &
                    SQRT(MIN(1.0_dp, MAX(0.0_dp, history_quality)))*history_scale
      balance = MAX(scalar_floor, MIN(projector_response_negative_balance_max, modal_floor))
      balance = MIN(1.0_dp, MAX(0.0_dp, balance))

   END FUNCTION dense_projector_negative_step_balance

! **************************************************************************************************
!> \brief Dense copy of the finite Projector-Hessian modal response algebra.
!> \param rot_hessian ...
!> \param rot_hessian_projector ...
!> \param step_rot ...
!> \param response_raw ...
!> \param response_exact ...
!> \param curv ...
!> \param grad_e ...
!> \param energy_base ...
!> \param response_base ...
!> \param base_weight ...
!> \param response_candidate ...
!> \param applied_weight ...
!> \param correction_drop ...
!> \param correction_negative_fraction ...
!> \param raw_negative_fraction ...
!> \param positive_drop ...
!> \param negative_drop ...
!> \param cross_drop ...
!> \param positive_coeff ...
!> \param negative_coeff ...
! **************************************************************************************************
   SUBROUTINE dense_projector_modal_response_candidate(rot_hessian, rot_hessian_projector, step_rot, &
                                                       response_raw, response_exact, curv, grad_e, energy_base, &
                                                       response_base, base_weight, response_candidate, applied_weight, &
                                                       correction_drop, correction_negative_fraction, &
                                                       raw_negative_fraction, positive_drop, negative_drop, cross_drop, &
                                                       positive_coeff, negative_coeff)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: rot_hessian, rot_hessian_projector
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: step_rot
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: response_raw, response_exact
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: curv, grad_e, energy_base, response_base
      REAL(KIND=dp), INTENT(IN)                          :: base_weight
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: response_candidate
      REAL(KIND=dp), INTENT(OUT) :: applied_weight, correction_drop, correction_negative_fraction, &
         raw_negative_fraction, positive_drop, negative_drop, cross_drop, positive_coeff, &
         negative_coeff

      INTEGER                                            :: imode, nrot_local
      REAL(KIND=dp) :: coeff_det, coeff_negative, coeff_positive, eval_abs_avg, g1, g2, h11, h12, &
         h22, mode_norm, mode_sign, mode_weight, negative_balance, negative_norm, &
         negative_weight_sum, positive_sum, positive_weight_sum, raw_negative_sum, &
         response_total_sum, spectral_weight, step_weight, step_weight_sum
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: mode_coeff, mode_eval, mode_rx, &
                                                            response_corr, response_neg, &
                                                            response_pos
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
      IF (base_weight <= EPSILON(base_weight)) RETURN

      nrot_local = SIZE(rot_hessian, 1)
      ALLOCATE (mode_coeff(nrot_local), mode_eval(nrot_local), mode_rx(SIZE(curv)), &
                response_corr(SIZE(curv)), response_neg(SIZE(curv)), response_pos(SIZE(curv)), &
                mode_vec(nrot_local, nrot_local))
      response_corr(:) = 0.0_dp
      response_neg(:) = 0.0_dp
      response_pos(:) = 0.0_dp
      mode_vec(:, :) = 0.5_dp*((rot_hessian_projector(:, :) - rot_hessian(:, :)) + &
                               TRANSPOSE(rot_hessian_projector(:, :) - rot_hessian(:, :)))
      CALL diamat_all(mode_vec, mode_eval)
      mode_coeff(:) = MATMUL(TRANSPOSE(mode_vec), step_rot)
      eval_abs_avg = SUM(ABS(mode_eval(:)))/REAL(MAX(1, nrot_local), KIND=dp)

      response_total_sum = 0.0_dp
      raw_negative_sum = 0.0_dp
      DO imode = 1, nrot_local
         IF (ABS(mode_coeff(imode)) <= EPSILON(mode_coeff(imode))) CYCLE
         mode_rx(:) = mode_coeff(imode)*MATMUL(response_exact - response_raw, mode_vec(:, imode))
         mode_norm = DOT_PRODUCT(curv(:)*mode_rx(:), mode_rx(:))
         response_total_sum = response_total_sum + MAX(0.0_dp, mode_norm)
         IF (mode_eval(imode) < 0.0_dp) raw_negative_sum = raw_negative_sum + MAX(0.0_dp, mode_norm)
      END DO
      IF (response_total_sum > EPSILON(response_total_sum)) &
         raw_negative_fraction = raw_negative_sum/response_total_sum
      positive_sum = MAX(0.0_dp, response_total_sum - raw_negative_sum)
      negative_balance = dense_projector_negative_response_balance(positive_sum, raw_negative_sum)

      step_weight_sum = DOT_PRODUCT(mode_coeff(:), mode_coeff(:))
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
      h12 = DOT_PRODUCT(curv(:)*response_pos(:), response_neg(:))
      h22 = DOT_PRODUCT(curv(:)*response_neg(:), response_neg(:))
      g1 = DOT_PRODUCT(grad_e(:) - curv(:)*energy_base(:), response_pos(:))
      g2 = DOT_PRODUCT(grad_e(:) - curv(:)*energy_base(:), response_neg(:))
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
      response_corr(:) = coeff_positive*response_pos(:) + coeff_negative*response_neg(:)
      response_candidate(:) = response_base(:) + response_corr(:)
      positive_drop = 0.5_dp*coeff_positive*coeff_positive*MAX(0.0_dp, h11)
      negative_drop = 0.5_dp*coeff_negative*coeff_negative*MAX(0.0_dp, h22)
      cross_drop = coeff_positive*coeff_negative*h12
      correction_drop = 0.5_dp*MAX(0.0_dp, DOT_PRODUCT(curv(:)*response_corr(:), response_corr(:)))
      negative_norm = DOT_PRODUCT(curv(:)*(coeff_negative*response_neg(:)), coeff_negative*response_neg(:))
      IF (step_weight_sum > EPSILON(step_weight_sum)) &
         applied_weight = (coeff_positive*positive_weight_sum + coeff_negative*negative_weight_sum)/step_weight_sum
      IF (correction_drop > EPSILON(correction_drop)) &
         correction_negative_fraction = negative_norm/(2.0_dp*correction_drop)
      DEALLOCATE (mode_coeff, mode_eval, mode_rx, response_corr, response_neg, response_pos, mode_vec)

   END SUBROUTINE dense_projector_modal_response_candidate

! **************************************************************************************************
!> \brief Dense mirror of the Direct trace compensation model.
!> \param coupled_weight ...
!> \param direct_delta_diag_avg ...
!> \param reference_delta_diag_avg ...
!> \return ...
! **************************************************************************************************
   FUNCTION dense_direct_trace_compensation(coupled_weight, direct_delta_diag_avg, reference_delta_diag_avg) &
      RESULT(compensation)
      REAL(KIND=dp), INTENT(IN)                          :: coupled_weight, direct_delta_diag_avg, &
                                                            reference_delta_diag_avg
      REAL(KIND=dp)                                      :: compensation

      REAL(KIND=dp)                                      :: positive_reference

      positive_reference = MAX(0.0_dp, reference_delta_diag_avg)
      compensation = MAX(0.0_dp, -coupled_weight*direct_delta_diag_avg)
      IF (direct_delta_diag_avg > positive_reference) &
         compensation = MAX(compensation, coupled_weight*(direct_delta_diag_avg - positive_reference))

   END FUNCTION dense_direct_trace_compensation

! **************************************************************************************************
!> \brief Dense mirror of accepted Direct trace split quality.
!> \param trace_diag_ema ...
!> \param trace_free_ema ...
!> \param trace_diag_quality_ema ...
!> \param trace_free_quality_ema ...
!> \param samples ...
!> \param fallback_quality ...
!> \return ...
! **************************************************************************************************
   FUNCTION dense_trace_split_history_quality(trace_diag_ema, trace_free_ema, trace_diag_quality_ema, &
                                              trace_free_quality_ema, samples, fallback_quality) &
      RESULT(history_quality)
      REAL(KIND=dp), INTENT(IN)                          :: trace_diag_ema, trace_free_ema, &
                                                            trace_diag_quality_ema, &
                                                            trace_free_quality_ema
      INTEGER, INTENT(IN)                                :: samples
      REAL(KIND=dp), INTENT(IN)                          :: fallback_quality
      REAL(KIND=dp)                                      :: history_quality

      REAL(KIND=dp) :: balance_quality, component_fraction, component_quality, diag_abs, &
         diag_quality, dominant_quality, free_abs, free_quality, trace_abs_sum

      history_quality = fallback_quality
      IF (samples < 8) RETURN

      diag_abs = ABS(trace_diag_ema)
      free_abs = ABS(trace_free_ema)
      trace_abs_sum = diag_abs + free_abs
      IF (trace_abs_sum <= EPSILON(trace_abs_sum)) RETURN

      diag_quality = 0.0_dp
      free_quality = 0.0_dp
      IF (diag_abs > EPSILON(diag_abs)) diag_quality = ABS(trace_diag_quality_ema)/diag_abs
      IF (free_abs > EPSILON(free_abs)) free_quality = ABS(trace_free_quality_ema)/free_abs
      component_fraction = MIN(diag_abs, free_abs)/trace_abs_sum
      IF (component_fraction < 0.10_dp) THEN
         dominant_quality = MERGE(diag_quality, free_quality, diag_abs >= free_abs)
         IF (dominant_quality > EPSILON(dominant_quality)) &
            history_quality = MAX(history_quality, MIN(2.5_dp, MAX(0.10_dp, dominant_quality)))
         RETURN
      END IF

      component_quality = (diag_abs*diag_quality + free_abs*free_quality)/trace_abs_sum
      IF (component_quality <= EPSILON(component_quality)) RETURN

      balance_quality = 1.0_dp - ABS(trace_diag_ema + trace_free_ema)/trace_abs_sum
      balance_quality = MIN(1.0_dp, MAX(0.0_dp, balance_quality))

      history_quality = MIN(2.5_dp, MAX(0.10_dp, component_quality))
      history_quality = history_quality*MAX(0.10_dp, balance_quality)
      history_quality = MIN(2.5_dp, MAX(0.10_dp, history_quality))

   END FUNCTION dense_trace_split_history_quality

! **************************************************************************************************
!> \brief Dense mirror of invariant-trace Direct history quality.
!> \param trace_ema ...
!> \param trace_quality_ema ...
!> \param trace_diag_ema ...
!> \param trace_free_ema ...
!> \param trace_diag_quality_ema ...
!> \param trace_free_quality_ema ...
!> \param samples ...
!> \param fallback_quality ...
!> \return ...
! **************************************************************************************************
   FUNCTION dense_trace_invariant_history_quality(trace_ema, trace_quality_ema, trace_diag_ema, &
                                                  trace_free_ema, trace_diag_quality_ema, &
                                                  trace_free_quality_ema, samples, fallback_quality) &
      RESULT(history_quality)
      REAL(KIND=dp), INTENT(IN) :: trace_ema, trace_quality_ema, trace_diag_ema, trace_free_ema, &
         trace_diag_quality_ema, trace_free_quality_ema
      INTEGER, INTENT(IN)                                :: samples
      REAL(KIND=dp), INTENT(IN)                          :: fallback_quality
      REAL(KIND=dp)                                      :: history_quality

      REAL(KIND=dp)                                      :: split_quality, trace_abs, trace_quality

      history_quality = fallback_quality
      IF (samples < 8) RETURN

      trace_abs = ABS(trace_ema)
      trace_quality = 0.0_dp
      IF (trace_abs > EPSILON(trace_abs)) trace_quality = ABS(trace_quality_ema)/trace_abs
      split_quality = dense_trace_split_history_quality(trace_diag_ema, trace_free_ema, &
                                                        trace_diag_quality_ema, trace_free_quality_ema, &
                                                        samples, fallback_quality)
      IF (trace_quality > EPSILON(trace_quality)) THEN
         trace_quality = MIN(1.0_dp, MAX(0.10_dp, trace_quality))
         history_quality = MIN(history_quality, trace_quality)
      END IF
      history_quality = MIN(history_quality, split_quality)
      history_quality = MIN(1.0_dp, MAX(0.10_dp, history_quality))

   END FUNCTION dense_trace_invariant_history_quality

! **************************************************************************************************
!> \brief Dense mirror of Schur accepted-history quality clipping.
!> \param predicted_sum ...
!> \param measured_sum ...
!> \param ema_quality ...
!> \param samples ...
!> \param min_samples ...
!> \param default_quality ...
!> \param min_quality ...
!> \param max_quality ...
!> \return ...
! **************************************************************************************************
   FUNCTION dense_schur_history_quality(predicted_sum, measured_sum, ema_quality, samples, min_samples, &
                                        default_quality, min_quality, max_quality) RESULT(quality)
      REAL(KIND=dp), INTENT(IN)                          :: predicted_sum, measured_sum, ema_quality
      INTEGER, INTENT(IN)                                :: samples, min_samples
      REAL(KIND=dp), INTENT(IN)                          :: default_quality, min_quality, max_quality
      REAL(KIND=dp)                                      :: quality

      quality = default_quality
      IF (samples >= min_samples) THEN
         IF (predicted_sum > EPSILON(predicted_sum) .AND. measured_sum > 0.0_dp) THEN
            quality = measured_sum/predicted_sum
         ELSE IF (ema_quality > 0.0_dp) THEN
            quality = ema_quality
         END IF
      END IF
      quality = MAX(min_quality, MIN(max_quality, quality))

   END FUNCTION dense_schur_history_quality

! **************************************************************************************************
!> \brief Dense mirror of inactive-offdiag selected residual-tail quality.
!> \param base_quality ...
!> \param selected_mermin_quality_ema ...
!> \param selected_predicted_sum ...
!> \param selected_measured_sum ...
!> \param selected_residual_quality ...
!> \param selected_residual_quality_ema ...
!> \param residual_reference_quality ...
!> \param selected_quality_samples ...
!> \param min_samples ...
!> \param residual_delta ...
!> \param tail_delta ...
!> \param calibrated_min_quality ...
!> \param max_quality ...
!> \return selected inactive-offdiag quality for candidate selection
! **************************************************************************************************
   FUNCTION dense_inactive_offdiag_selected_tail_quality(base_quality, selected_mermin_quality_ema, &
                                                        selected_predicted_sum, selected_measured_sum, &
                                                        selected_residual_quality, &
                                                        selected_residual_quality_ema, &
                                                        residual_reference_quality, &
                                                        selected_quality_samples, min_samples, &
                                                        residual_delta, tail_delta, &
                                                        calibrated_min_quality, max_quality) &
      RESULT(quality)
      REAL(KIND=dp), INTENT(IN)                          :: base_quality, selected_mermin_quality_ema, &
                                                            selected_predicted_sum, selected_measured_sum, &
                                                            selected_residual_quality, &
                                                            selected_residual_quality_ema, &
                                                            residual_reference_quality, residual_delta, &
                                                            tail_delta, calibrated_min_quality, max_quality
      INTEGER, INTENT(IN)                                :: min_samples, selected_quality_samples
      REAL(KIND=dp)                                      :: quality

      REAL(KIND=dp)                                      :: mermin_factor, mermin_quality, residual_factor, &
                                                            residual_ref, residual_quality

      quality = MAX(0.0_dp, MIN(max_quality, base_quality))
      IF (residual_delta > tail_delta) RETURN
      IF (selected_quality_samples < 2*min_samples) RETURN

      mermin_quality = quality
      IF (selected_predicted_sum > EPSILON(selected_predicted_sum) .AND. &
          selected_measured_sum > 0.0_dp) mermin_quality = selected_measured_sum/selected_predicted_sum
      IF (selected_mermin_quality_ema > 0.0_dp) &
         mermin_quality = MIN(mermin_quality, selected_mermin_quality_ema)
      mermin_quality = MAX(0.0_dp, MIN(max_quality, mermin_quality))
      IF (mermin_quality < 0.10_dp) THEN
         quality = 0.0_dp
         RETURN
      END IF

      residual_ref = MAX(residual_reference_quality, calibrated_min_quality, EPSILON(1.0_dp))
      residual_quality = MAX(0.0_dp, MIN(1.0_dp, selected_residual_quality_ema))
      IF (selected_residual_quality > EPSILON(selected_residual_quality)) &
         residual_quality = MIN(residual_quality, &
                                0.5_dp*(residual_quality + &
                                        MAX(0.0_dp, MIN(1.0_dp, selected_residual_quality))))

      mermin_factor = MIN(1.0_dp, SQRT(MAX(0.0_dp, mermin_quality)/ &
                                       MAX(calibrated_min_quality, EPSILON(1.0_dp))))
      residual_factor = MIN(1.0_dp, SQRT(residual_quality/residual_ref))

      quality = quality*mermin_factor*residual_factor
      quality = MAX(0.0_dp, MIN(max_quality, quality))

   END FUNCTION dense_inactive_offdiag_selected_tail_quality

! **************************************************************************************************
!> \brief Dense mirror of offdiag residual-tail history cap.
!> \param base_quality ...
!> \param residual_quality ...
!> \param residual_quality_ema ...
!> \param residual_samples ...
!> \param residual_delta ...
!> \param tail_delta ...
!> \param calibrated_min_quality ...
!> \param max_quality ...
!> \return offdiag candidate selection quality capped by residual consistency
! **************************************************************************************************
   FUNCTION dense_offdiag_tail_residual_quality(base_quality, residual_quality, residual_quality_ema, &
                                                residual_samples, residual_delta, tail_delta, &
                                                calibrated_min_quality, max_quality) RESULT(quality)
      REAL(KIND=dp), INTENT(IN)                          :: base_quality, residual_quality, &
                                                            residual_quality_ema, residual_delta, &
                                                            tail_delta, calibrated_min_quality, &
                                                            max_quality
      INTEGER, INTENT(IN)                                :: residual_samples
      REAL(KIND=dp)                                      :: quality

      REAL(KIND=dp)                                      :: residual_cap, residual_ref, residual_signal

      quality = MAX(0.0_dp, MIN(max_quality, base_quality))
      IF (residual_delta > tail_delta) RETURN
      IF (residual_samples < 8) RETURN

      residual_signal = MAX(0.0_dp, MIN(1.0_dp, residual_quality_ema))
      IF (residual_quality > EPSILON(residual_quality)) &
         residual_signal = MIN(residual_signal, &
                               0.5_dp*(residual_signal + &
                                       MAX(0.0_dp, MIN(1.0_dp, residual_quality))))
      IF (residual_signal <= EPSILON(residual_signal)) THEN
         quality = 0.0_dp
         RETURN
      END IF

      residual_ref = MAX(calibrated_min_quality, EPSILON(1.0_dp))
      residual_cap = max_quality*MIN(1.0_dp, residual_signal/residual_ref)
      quality = MIN(quality, residual_cap)
      quality = MAX(0.0_dp, MIN(max_quality, quality))

   END FUNCTION dense_offdiag_tail_residual_quality

! **************************************************************************************************
!> \brief Dense mirror of projector-response candidate gate.
!> \param samples ...
!> \param quality ...
!> \param quality_ema ...
!> \param predicted_sum ...
!> \param measured_sum ...
!> \param negative_ema ...
!> \param min_quality ...
!> \param max_quality ...
!> \param underpredict_max_quality ...
!> \return ...
! **************************************************************************************************
   FUNCTION dense_projector_response_candidate_enabled(samples, quality, quality_ema, predicted_sum, &
                                                       measured_sum, negative_ema, min_quality, max_quality, &
                                                       underpredict_max_quality) RESULT(enabled)
      INTEGER, INTENT(IN)                                :: samples
      REAL(KIND=dp), INTENT(IN) :: quality, quality_ema, predicted_sum, measured_sum, &
         negative_ema, min_quality, max_quality, underpredict_max_quality
      LOGICAL                                            :: enabled

      LOGICAL                                            :: aggregate_recovered_quality, &
                                                            normal_quality, underpredicted_quality
      REAL(KIND=dp)                                      :: aggregate_quality

      enabled = .FALSE.
      IF (samples < 8) RETURN
      IF (quality_ema < min_quality) RETURN
      IF (quality < 0.5_dp*min_quality) RETURN

      aggregate_quality = 0.0_dp
      IF (predicted_sum > EPSILON(predicted_sum) .AND. measured_sum > 0.0_dp) &
         aggregate_quality = measured_sum/predicted_sum

      normal_quality = quality_ema <= max_quality .AND. quality <= 2.0_dp*max_quality .AND. &
                       negative_ema <= projector_response_negative_balance_max
      aggregate_recovered_quality = aggregate_quality >= min_quality .AND. aggregate_quality <= max_quality .AND. &
                                    quality >= min_quality .AND. quality <= 2.0_dp*max_quality .AND. &
                                    negative_ema <= projector_response_negative_balance_floor
      underpredicted_quality = quality_ema > max_quality .AND. &
                               quality_ema <= underpredict_max_quality .AND. &
                               quality <= underpredict_max_quality .AND. &
                               predicted_sum > EPSILON(predicted_sum) .AND. &
                               measured_sum > predicted_sum .AND. &
                               negative_ema <= projector_response_negative_balance_floor

      enabled = normal_quality .OR. aggregate_recovered_quality .OR. underpredicted_quality

   END FUNCTION dense_projector_response_candidate_enabled

! **************************************************************************************************
!> \brief Dense mirror of finite Projector-response residual taper.
!> \param residual_delta ...
!> \param raw_residual_delta ...
!> \return ...
! **************************************************************************************************
   FUNCTION dense_projector_response_residual_weight(residual_delta, raw_residual_delta) RESULT(residual_weight)
      REAL(KIND=dp), INTENT(IN)                          :: residual_delta, raw_residual_delta
      REAL(KIND=dp)                                      :: residual_weight

      REAL(KIND=dp), PARAMETER :: coupled_tail_delta = 1.0E-4_dp, floor_delta = 3.0E-6_dp, &
         full_delta = 1.2E-5_dp

      REAL(KIND=dp)                                      :: raw_signal, response_delta

      response_delta = MAX(0.0_dp, residual_delta)
      IF (raw_residual_delta > response_delta) THEN
         raw_signal = MIN(MAX(0.0_dp, raw_residual_delta), coupled_tail_delta)
         response_delta = MAX(response_delta, SQRT(response_delta*raw_signal))
      END IF

      residual_weight = response_delta/MAX(EPSILON(1.0_dp), full_delta)
      residual_weight = MAX(0.0_dp, MIN(1.0_dp, residual_weight))
      IF (response_delta < floor_delta) &
         residual_weight = residual_weight* &
                           SQRT(MAX(0.0_dp, response_delta/MAX(EPSILON(1.0_dp), floor_delta)))

   END FUNCTION dense_projector_response_residual_weight

! **************************************************************************************************
!> \brief Dense mirror of adverse negative-mode Projector-response history.
!> \param negative_norm ...
!> \param correction_norm ...
!> \param negative_gdot ...
!> \return ...
! **************************************************************************************************
   FUNCTION dense_projector_response_adverse_fraction(negative_norm, correction_norm, &
                                                      negative_gdot) RESULT(fraction)
      REAL(KIND=dp), INTENT(IN)                          :: negative_norm, correction_norm, &
                                                            negative_gdot
      REAL(KIND=dp)                                      :: fraction

      fraction = 0.0_dp
      IF (negative_norm <= EPSILON(negative_norm)) RETURN
      IF (correction_norm <= EPSILON(correction_norm)) RETURN
      IF (negative_gdot >= 0.0_dp) RETURN

      fraction = MAX(0.0_dp, MIN(1.0_dp, negative_norm/correction_norm))

   END FUNCTION dense_projector_response_adverse_fraction

! **************************************************************************************************
!> \brief Dense mirror of selected-history Direct quality maturity.
!> \param base_quality ...
!> \param predicted_sum ...
!> \param measured_sum ...
!> \param ema_quality ...
!> \param samples ...
!> \param min_samples ...
!> \param response_min_samples ...
!> \param min_quality ...
!> \param max_quality ...
!> \return ...
! **************************************************************************************************
   FUNCTION dense_direct_selected_history_quality(base_quality, predicted_sum, measured_sum, ema_quality, &
                                                  samples, min_samples, response_min_samples, min_quality, &
                                                  max_quality) RESULT(quality)
      REAL(KIND=dp), INTENT(IN)                          :: base_quality, predicted_sum, &
                                                            measured_sum, ema_quality
      INTEGER, INTENT(IN)                                :: samples, min_samples, &
                                                            response_min_samples
      REAL(KIND=dp), INTENT(IN)                          :: min_quality, max_quality
      REAL(KIND=dp)                                      :: quality

      INTEGER                                            :: maturity_samples
      REAL(KIND=dp)                                      :: maturity, selected_quality

      quality = base_quality
      IF (samples < min_samples) RETURN

      selected_quality = base_quality
      IF (predicted_sum > EPSILON(predicted_sum) .AND. measured_sum > 0.0_dp) THEN
         selected_quality = measured_sum/predicted_sum
      ELSE IF (ema_quality > 0.0_dp) THEN
         selected_quality = ema_quality
      END IF
      selected_quality = MAX(min_quality, MIN(max_quality, selected_quality))

      maturity_samples = MAX(1, 2*response_min_samples - min_samples + 1)
      maturity = REAL(MAX(0, samples - min_samples + 1), KIND=dp)/ &
                 REAL(maturity_samples, KIND=dp)
      maturity = MIN(1.0_dp, MAX(0.0_dp, maturity))

      quality = (1.0_dp - maturity)*base_quality + maturity*selected_quality
      quality = MAX(min_quality, MIN(max_quality, quality))

   END FUNCTION dense_direct_selected_history_quality

! **************************************************************************************************
!> \brief Dense mirror of selected-history Direct probe quality.
!> \param base_quality ...
!> \param selected_quality ...
!> \param samples ...
!> \param min_samples ...
!> \param response_min_samples ...
!> \param max_quality ...
!> \return maturity-blended probe quality without a promotion-quality floor
! **************************************************************************************************
   FUNCTION dense_direct_selected_probe_quality(base_quality, selected_quality, samples, min_samples, &
                                                response_min_samples, max_quality) RESULT(quality)
      REAL(KIND=dp), INTENT(IN)                          :: base_quality, selected_quality
      INTEGER, INTENT(IN)                                :: samples, min_samples, &
                                                            response_min_samples
      REAL(KIND=dp), INTENT(IN)                          :: max_quality
      REAL(KIND=dp)                                      :: quality

      INTEGER                                            :: maturity_samples
      REAL(KIND=dp)                                      :: maturity, selected_quality_clipped

      quality = base_quality
      IF (samples < min_samples) RETURN

      selected_quality_clipped = MAX(0.0_dp, MIN(max_quality, selected_quality))
      maturity_samples = MAX(1, 2*response_min_samples - min_samples + 1)
      maturity = REAL(MAX(0, samples - min_samples + 1), KIND=dp)/ &
                 REAL(maturity_samples, KIND=dp)
      maturity = MIN(1.0_dp, MAX(0.0_dp, maturity))

      quality = (1.0_dp - maturity)*base_quality + maturity*selected_quality_clipped
      quality = MAX(0.0_dp, MIN(max_quality, quality))

   END FUNCTION dense_direct_selected_probe_quality

! **************************************************************************************************
!> \brief Dense mirror of floor-calibrated selected Direct probe quality.
!> \param floor_quality ...
!> \param selected_quality ...
!> \param samples ...
!> \param min_samples ...
!> \param response_min_samples ...
!> \param max_quality ...
!> \return floor-calibrated probe quality with selected-history consistency
! **************************************************************************************************
   FUNCTION dense_direct_floor_selected_quality(floor_quality, selected_quality, samples, min_samples, &
                                                response_min_samples, max_quality) &
      RESULT(quality)
      REAL(KIND=dp), INTENT(IN)                          :: floor_quality, selected_quality
      INTEGER, INTENT(IN)                                :: samples, min_samples, &
                                                            response_min_samples
      REAL(KIND=dp), INTENT(IN)                          :: max_quality
      REAL(KIND=dp)                                      :: quality

      quality = floor_quality

      quality = dense_direct_selected_probe_quality(quality, selected_quality, samples, &
                                                    min_samples, response_min_samples, max_quality)

   END FUNCTION dense_direct_floor_selected_quality

! **************************************************************************************************
!> \brief Dense mirror of recent accepted-drop bootstrap quality.
!> \param history_quality ...
!> \param recent_quality ...
!> \param ema_quality ...
!> \param samples ...
!> \param min_samples ...
!> \param min_quality ...
!> \param max_quality ...
!> \return aggregate quality capped by recent/EMA Mermin-drop consistency
! **************************************************************************************************
   FUNCTION dense_floor_bootstrap_history_quality(history_quality, recent_quality, ema_quality, &
                                                  samples, min_samples, min_quality, max_quality) &
      RESULT(quality)
      REAL(KIND=dp), INTENT(IN)                          :: history_quality, recent_quality, &
                                                            ema_quality
      INTEGER, INTENT(IN)                                :: samples, min_samples
      REAL(KIND=dp), INTENT(IN)                          :: min_quality, max_quality
      REAL(KIND=dp)                                      :: quality

      REAL(KIND=dp)                                      :: recent_cap, recent_clipped, recent_ref

      quality = MAX(0.0_dp, MIN(max_quality, history_quality))
      IF (samples < min_samples) RETURN

      recent_clipped = MAX(0.0_dp, MIN(max_quality, recent_quality))
      IF (recent_clipped <= EPSILON(recent_clipped)) RETURN

      recent_ref = quality
      IF (ema_quality > EPSILON(ema_quality)) recent_ref = MIN(recent_ref, &
                                                               MAX(min_quality, MIN(max_quality, ema_quality)))
      recent_cap = SQRT(recent_clipped*MAX(min_quality, recent_ref))
      recent_cap = MAX(min_quality, MIN(max_quality, recent_cap))
      quality = MIN(quality, recent_cap)

   END FUNCTION dense_floor_bootstrap_history_quality

! **************************************************************************************************
!> \brief Dense mirror of the Direct floor-density bootstrap gate.
!> \param residual_delta ...
!> \param mermin_quality_ema ...
!> \param mermin_samples ...
!> \param floor_abs_avg ...
!> \param floor_exact_abs_avg ...
!> \param floor_low_curv_fraction ...
!> \param projected_abs_avg ...
!> \return conservative Direct floor-density bootstrap weight
! **************************************************************************************************
   FUNCTION dense_direct_floor_density_bootstrap_weight(residual_delta, mermin_quality_ema, &
                                                        mermin_samples, floor_abs_avg, &
                                                        floor_exact_abs_avg, floor_low_curv_fraction, &
                                                        projected_abs_avg) RESULT(weight)
      REAL(KIND=dp), INTENT(IN)                          :: residual_delta, mermin_quality_ema, &
                                                            floor_abs_avg, floor_exact_abs_avg, &
                                                            floor_low_curv_fraction, projected_abs_avg
      INTEGER, INTENT(IN)                                :: mermin_samples
      REAL(KIND=dp)                                      :: weight

      REAL(KIND=dp), PARAMETER                           :: bootstrap_weight = 2.5E-2_dp, &
                                                            lowcurv_max = 2.5E-1_dp, &
                                                            mermin_min = 8.0E-1_dp, &
                                                            response_full_delta = 1.2E-5_dp, &
                                                            trace_active_max_weight = 9.5E-1_dp
      INTEGER, PARAMETER                                 :: calibration_min_samples = 4
      REAL(KIND=dp)                                      :: low_curv_quality, source_limit

      weight = 0.0_dp
      IF (residual_delta > response_full_delta) RETURN
      IF (mermin_samples < calibration_min_samples) RETURN
      IF (mermin_quality_ema < mermin_min) RETURN
      IF (floor_abs_avg <= EPSILON(floor_abs_avg)) RETURN
      IF (floor_exact_abs_avg <= EPSILON(floor_exact_abs_avg)) RETURN
      IF (floor_low_curv_fraction > lowcurv_max) RETURN
      source_limit = MAX(1.0E1_dp*floor_exact_abs_avg, 1.0E-1_dp*floor_abs_avg)
      IF (projected_abs_avg > source_limit) RETURN

      low_curv_quality = 1.0_dp - floor_low_curv_fraction/lowcurv_max
      low_curv_quality = MAX(0.0_dp, MIN(1.0_dp, low_curv_quality))
      weight = bootstrap_weight*MIN(1.0_dp, mermin_quality_ema)*low_curv_quality
      weight = MIN(0.5_dp*trace_active_max_weight, MAX(0.0_dp, weight))

   END FUNCTION dense_direct_floor_density_bootstrap_weight

! **************************************************************************************************
!> \brief Dense mirror of the physical Direct tail-probe bootstrap gate.
!> \param candidate_used ...
!> \param promoted ...
!> \param candidate_enabled ...
!> \param residual_delta ...
!> \param coupled_quality_samples ...
!> \param selected_quality_samples ...
!> \param selected_quality_ema ...
!> \param selected_quality ...
!> \param shadow_quality ...
!> \param coupled_residual_quality_ema ...
!> \param weighted_trust_ratio ...
!> \param abs_trust_ratio ...
!> \param model_gap_ratio ...
!> \return true if a trust-safe Direct tail probe may collect selected history
! **************************************************************************************************
   FUNCTION dense_direct_tail_probe_ready(candidate_used, promoted, candidate_enabled, residual_delta, &
                                          coupled_quality_samples, selected_quality_samples, &
                                          selected_quality_ema, selected_quality, shadow_quality, &
                                          coupled_residual_quality_ema, weighted_trust_ratio, &
                                          abs_trust_ratio, model_gap_ratio) &
      RESULT(probe_ready)
      LOGICAL, INTENT(IN)                                :: candidate_used, promoted, &
                                                            candidate_enabled
      REAL(KIND=dp), INTENT(IN)                          :: residual_delta
      INTEGER, INTENT(IN)                                :: coupled_quality_samples, &
                                                            selected_quality_samples
      REAL(KIND=dp), INTENT(IN) :: selected_quality_ema, selected_quality, shadow_quality, &
         coupled_residual_quality_ema, weighted_trust_ratio, abs_trust_ratio, model_gap_ratio
      LOGICAL                                            :: probe_ready

      LOGICAL                                            :: residual_tail_stalled, &
                                                            selected_quality_positive

      probe_ready = .FALSE.
      IF (.NOT. candidate_used) RETURN
      IF (promoted .OR. candidate_enabled) RETURN
      IF (residual_delta > 1.0E-4_dp) RETURN
      IF (coupled_quality_samples < 4) RETURN
      IF (shadow_quality < 1.0E-1_dp) RETURN
      IF (weighted_trust_ratio > 7.5E-1_dp) RETURN
      IF (abs_trust_ratio > 7.5E-1_dp) RETURN
      IF (model_gap_ratio > 7.5E-1_dp) RETURN

      IF (selected_quality_samples < 4) THEN
         probe_ready = .TRUE.
         RETURN
      END IF

      selected_quality_positive = selected_quality_ema >= 5.0E-2_dp .AND. selected_quality >= 5.0E-2_dp
      probe_ready = selected_quality_positive .AND. MOD(coupled_quality_samples, 16) == 0
      IF (.NOT. probe_ready) THEN
         residual_tail_stalled = coupled_residual_quality_ema < 0.25_dp
         selected_quality_positive = selected_quality_ema >= 5.0E-2_dp .AND. selected_quality >= 2.5E-2_dp
         probe_ready = residual_tail_stalled .AND. selected_quality_positive .AND. &
                       MOD(coupled_quality_samples, 8) == 0
      END IF

   END FUNCTION dense_direct_tail_probe_ready

! **************************************************************************************************
!> \brief Dense mirror of the spectral selected-history candidate promotion gate.
!> \param selected_samples ...
!> \param selected_failures ...
!> \param selected_quality_ema ...
!> \param selected_quality ...
!> \param selected_residual_quality_ema ...
!> \param selected_residual_quality ...
!> \param selected_drop_sum ...
!> \param selected_measured_sum ...
!> \param shadow_samples ...
!> \param shadow_drop_sum ...
!> \param shadow_measured_sum ...
!> \param shadow_quality_ema ...
!> \param shadow_quality ...
!> \param shadow_floor_ema ...
!> \param shadow_dot_ema ...
!> \param shadow_norm_ema ...
!> \return true if selected spectral history may promote the candidate
! **************************************************************************************************
   FUNCTION dense_spectral_selected_history_promotion(selected_samples, selected_failures, &
                                                      selected_quality_ema, selected_quality, &
                                                      selected_residual_quality_ema, &
                                                      selected_residual_quality, selected_drop_sum, &
                                                      selected_measured_sum, shadow_samples, &
                                                      shadow_drop_sum, shadow_measured_sum, &
                                                      shadow_quality_ema, shadow_quality, &
                                                      shadow_floor_ema, shadow_dot_ema, &
                                                      shadow_norm_ema) RESULT(enabled)
      INTEGER, INTENT(IN)                                :: selected_failures, selected_samples, &
                                                            shadow_samples
      REAL(KIND=dp), INTENT(IN)                          :: selected_drop_sum, selected_measured_sum, &
                                                            selected_quality, selected_quality_ema, &
                                                            selected_residual_quality, &
                                                            selected_residual_quality_ema, &
                                                            shadow_dot_ema, shadow_drop_sum, &
                                                            shadow_floor_ema, shadow_measured_sum, &
                                                            shadow_norm_ema, shadow_quality, &
                                                            shadow_quality_ema
      LOGICAL                                            :: enabled

      LOGICAL                                            :: ema_ok, recent_ok, selected_sums_ok, &
                                                            shadow_history_ok

      enabled = .FALSE.
      IF (selected_failures >= 2) RETURN
      IF (shadow_dot_ema <= 0.0_dp .OR. shadow_norm_ema <= 0.0_dp) RETURN

      selected_sums_ok = selected_drop_sum > EPSILON(selected_drop_sum) .AND. &
                         selected_measured_sum > 0.0_dp
      ema_ok = selected_quality_ema >= 1.0E-1_dp .OR. &
               selected_residual_quality_ema >= 1.0E-1_dp
      recent_ok = selected_quality >= 5.0E-2_dp .OR. selected_residual_quality >= 1.0E-1_dp

      IF (selected_samples >= 8 .AND. selected_sums_ok .AND. ema_ok .AND. recent_ok) THEN
         enabled = .TRUE.
         RETURN
      END IF

      IF (selected_samples >= 3 .AND. selected_sums_ok .AND. &
          selected_residual_quality_ema >= 1.0E-1_dp .AND. &
          selected_residual_quality >= 1.0E-1_dp .AND. ema_ok .AND. recent_ok) THEN
         enabled = .TRUE.
         RETURN
      END IF

      shadow_history_ok = shadow_samples >= 8 .AND. &
                          shadow_drop_sum > EPSILON(shadow_drop_sum) .AND. &
                          shadow_measured_sum > 0.0_dp .AND. &
                          shadow_floor_ema > EPSILON(shadow_floor_ema) .AND. &
                          shadow_quality_ema >= 1.0E-1_dp .AND. shadow_quality >= 5.0E-2_dp
      enabled = selected_samples > 0 .AND. selected_sums_ok .AND. ema_ok .AND. recent_ok .AND. &
                shadow_history_ok

   END FUNCTION dense_spectral_selected_history_promotion

! **************************************************************************************************
!> \brief Dense mirror of the full-density Direct residual-tail close-call gate.
!> \param residual_delta ...
!> \param candidate_expected_drop ...
!> \param reference_expected_drop ...
!> \param quality_ema ...
!> \param samples ...
!> \param selected_residual_quality ...
!> \param selected_residual_quality_ema ...
!> \param selected_samples ...
!> \param dense_residual_quality ...
!> \return true if the candidate is a calibrated residual-tail close call
! **************************************************************************************************
   FUNCTION dense_density_full_direct_tail_probe(residual_delta, candidate_expected_drop, &
                                                 reference_expected_drop, quality_ema, samples, &
                                                 selected_residual_quality, selected_residual_quality_ema, &
                                                 selected_samples, &
                                                 dense_residual_quality) RESULT(do_probe)
      REAL(KIND=dp), INTENT(IN)                          :: residual_delta, candidate_expected_drop, &
                                                            reference_expected_drop, quality_ema
      INTEGER, INTENT(IN)                                :: samples
      REAL(KIND=dp), INTENT(IN)                          :: selected_residual_quality, &
                                                            selected_residual_quality_ema
      INTEGER, INTENT(IN)                                :: selected_samples
      REAL(KIND=dp), INTENT(IN)                          :: dense_residual_quality
      LOGICAL                                            :: do_probe

      LOGICAL                                            :: stale_residual_probe
      REAL(KIND=dp)                                      :: required_drop_ratio

      do_probe = .FALSE.
      IF (residual_delta > 1.0E-4_dp) RETURN
      IF (samples < 8) RETURN
      IF (quality_ema < 1.0E-1_dp) RETURN
      IF (candidate_expected_drop <= 1.0E-12_dp) RETURN

      stale_residual_probe = .FALSE.
      IF (selected_samples >= 4) THEN
         stale_residual_probe = dense_residual_quality <= 5.0E-1_dp .AND. &
                                selected_residual_quality_ema < 8.0E-1_dp* &
                                MAX(dense_residual_quality, 2.5E-1_dp) .AND. &
                                quality_ema >= 2.0E-1_dp .AND. &
                                MOD(samples, 8) == 0
      END IF

      required_drop_ratio = 8.0E-1_dp
      IF (stale_residual_probe) required_drop_ratio = MIN(required_drop_ratio, 2.0E-1_dp)
      IF (candidate_expected_drop < required_drop_ratio*MAX(reference_expected_drop, 1.0E-12_dp)) RETURN

      IF (selected_samples >= 4) THEN
         IF (.NOT. stale_residual_probe) THEN
            IF (selected_residual_quality_ema < 8.0E-1_dp*MAX(dense_residual_quality, 2.5E-1_dp)) &
               RETURN
            IF (selected_residual_quality < 4.0E-1_dp*MAX(dense_residual_quality, 2.5E-1_dp)) &
               RETURN
         END IF
      END IF

      do_probe = .TRUE.

   END FUNCTION dense_density_full_direct_tail_probe

! **************************************************************************************************
!> \brief Dense mirror of the density-descent residual-tail close-call gate.
!> \param residual_delta ...
!> \param candidate_expected_drop ...
!> \param reference_expected_drop ...
!> \param quality_for_selection ...
!> \param gap_recoverable_fraction ...
!> \param weighted_trust_ratio ...
!> \param abs_trust_ratio ...
!> \param selected_residual_quality ...
!> \param selected_residual_quality_ema ...
!> \param selected_samples ...
!> \param dense_residual_quality ...
!> \return true if the density-descent candidate is a trust-safe residual-tail close call
! **************************************************************************************************
   FUNCTION dense_density_descent_tail_probe(residual_delta, candidate_expected_drop, &
                                             reference_expected_drop, quality_for_selection, &
                                             gap_recoverable_fraction, weighted_trust_ratio, &
                                             abs_trust_ratio, selected_residual_quality, &
                                             selected_residual_quality_ema, selected_samples, &
                                             dense_residual_quality) RESULT(do_probe)
      REAL(KIND=dp), INTENT(IN)                          :: residual_delta, candidate_expected_drop, &
                                                            reference_expected_drop, quality_for_selection, &
                                                            gap_recoverable_fraction, weighted_trust_ratio, &
                                                            abs_trust_ratio, selected_residual_quality, &
                                                            selected_residual_quality_ema
      INTEGER, INTENT(IN)                                :: selected_samples
      REAL(KIND=dp), INTENT(IN)                          :: dense_residual_quality
      LOGICAL                                            :: do_probe

      REAL(KIND=dp)                                      :: residual_reference

      do_probe = .FALSE.
      IF (residual_delta > 1.0E-4_dp) RETURN
      IF (quality_for_selection < 5.0E-2_dp) RETURN
      IF (gap_recoverable_fraction < 5.0E-1_dp) RETURN
      IF (candidate_expected_drop <= 1.0E-12_dp) RETURN
      IF (candidate_expected_drop < 8.0E-1_dp*MAX(reference_expected_drop, 1.0E-12_dp)) RETURN
      IF (weighted_trust_ratio > 1.0_dp) RETURN
      IF (abs_trust_ratio > 1.0_dp) RETURN

      IF (selected_samples >= 4) THEN
         residual_reference = MAX(dense_residual_quality, 2.5E-1_dp)
         IF (selected_residual_quality_ema < 8.0E-1_dp*residual_reference) RETURN
         IF (selected_residual_quality < 4.0E-1_dp*residual_reference) RETURN
      END IF

      do_probe = .TRUE.

   END FUNCTION dense_density_descent_tail_probe

! **************************************************************************************************
!> \brief Dense mirror of the residual-tail Direct-density shadow probe drop.
!> \param residual_delta ...
!> \param capped_expected_drop ...
!> \param shadow_expected_drop ...
!> \param weighted_trust_ratio ...
!> \param abs_trust_ratio ...
!> \return ...
! **************************************************************************************************
   FUNCTION dense_density_full_direct_probe_drop(residual_delta, capped_expected_drop, &
                                                 shadow_expected_drop, weighted_trust_ratio, &
                                                 abs_trust_ratio) RESULT(probe_drop)
      REAL(KIND=dp), INTENT(IN)                          :: residual_delta, capped_expected_drop, &
                                                            shadow_expected_drop, &
                                                            weighted_trust_ratio, abs_trust_ratio
      REAL(KIND=dp)                                      :: probe_drop

      probe_drop = MAX(0.0_dp, capped_expected_drop)
      IF (residual_delta > 1.0E-4_dp) RETURN
      IF (shadow_expected_drop <= probe_drop) RETURN
      IF (shadow_expected_drop <= 1.0E-12_dp) RETURN
      IF (weighted_trust_ratio > 7.5E-1_dp) RETURN
      IF (abs_trust_ratio > 7.5E-1_dp) RETURN

      probe_drop = shadow_expected_drop

   END FUNCTION dense_density_full_direct_probe_drop

! **************************************************************************************************
!> \brief Dense mirror of the strict selected-step residual tail cap.
!> \param selected_residual_quality ...
!> \param selected_residual_quality_ema ...
!> \return ...
! **************************************************************************************************
   FUNCTION dense_selected_tail_residual_quality(selected_residual_quality, &
                                                 selected_residual_quality_ema) RESULT(quality)
      REAL(KIND=dp), INTENT(IN)                          :: selected_residual_quality, &
                                                            selected_residual_quality_ema
      REAL(KIND=dp)                                      :: quality

      quality = MAX(0.0_dp, MIN(1.0_dp, selected_residual_quality_ema))
      quality = quality*SQRT(MAX(0.0_dp, MIN(1.0_dp, selected_residual_quality/2.5E-1_dp)))

   END FUNCTION dense_selected_tail_residual_quality

! **************************************************************************************************
!> \brief Dense mirror of the EMA-stabilized selected residual tail mix cap.
!> \param selected_residual_quality ...
!> \param selected_residual_quality_ema ...
!> \return ...
! **************************************************************************************************
   FUNCTION dense_selected_tail_mix_residual_quality(selected_residual_quality, &
                                                     selected_residual_quality_ema) RESULT(quality)
      REAL(KIND=dp), INTENT(IN)                          :: selected_residual_quality, &
                                                            selected_residual_quality_ema
      REAL(KIND=dp)                                      :: quality

      REAL(KIND=dp)                                      :: recent_quality

      quality = MAX(0.0_dp, MIN(1.0_dp, selected_residual_quality_ema))
      IF (quality <= EPSILON(quality)) RETURN

      IF (selected_residual_quality > EPSILON(selected_residual_quality)) THEN
         recent_quality = SQRT(MAX(0.0_dp, MIN(1.0_dp, selected_residual_quality/2.5E-1_dp)))
         quality = SQRT(quality*recent_quality)
      ELSE
         quality = 0.5_dp*quality
      END IF
      quality = MAX(0.0_dp, MIN(1.0_dp, quality))

   END FUNCTION dense_selected_tail_mix_residual_quality

! **************************************************************************************************
!> \brief Dense mirror of selected Direct-density active-Hessian tail mix weight.
!> \param residual_delta ...
!> \param selected_mermin_quality ...
!> \param selected_residual_quality ...
!> \param selected_samples ...
!> \return ...
! **************************************************************************************************
   FUNCTION dense_direct_density_selected_mix_weight(residual_delta, selected_mermin_quality, &
                                                     selected_residual_quality, selected_samples) RESULT(weight)
      REAL(KIND=dp), INTENT(IN)                          :: residual_delta, selected_mermin_quality, &
                                                            selected_residual_quality
      INTEGER, INTENT(IN)                                :: selected_samples
      REAL(KIND=dp)                                      :: weight

      INTEGER                                            :: maturity_samples
      REAL(KIND=dp)                                      :: maturity, mermin_weight, residual_weight

      weight = 0.0_dp
      IF (residual_delta > 1.0E-4_dp) RETURN
      IF (selected_samples < 4) RETURN
      IF (selected_mermin_quality < 2.5E-2_dp) RETURN
      IF (selected_residual_quality <= EPSILON(selected_residual_quality)) RETURN

      maturity_samples = 13
      maturity = REAL(MAX(0, selected_samples - 4 + 1), KIND=dp)/REAL(maturity_samples, KIND=dp)
      maturity = MIN(1.0_dp, MAX(0.0_dp, maturity))
      mermin_weight = SQRT(MIN(1.0_dp, MAX(0.0_dp, selected_mermin_quality)))
      residual_weight = SQRT(MIN(1.0_dp, MAX(0.0_dp, selected_residual_quality)))
      weight = 0.25_dp*0.95_dp*maturity*mermin_weight*residual_weight
      weight = MIN(0.25_dp*0.95_dp, MAX(0.0_dp, weight))

   END FUNCTION dense_direct_density_selected_mix_weight

! **************************************************************************************************
!> \brief Dense mirror of the residual-tail dense metric-response bridge gate.
!> \param residual_delta ...
!> \param coupled_quality ...
!> \param coupled_quality_ema ...
!> \param coupled_samples ...
!> \param dense_residual_quality ...
!> \param mermin_quality_ema ...
!> \param mermin_samples ...
!> \return ...
! **************************************************************************************************
   FUNCTION dense_metric_tail_bridge_ready(residual_delta, coupled_quality, coupled_quality_ema, &
                                           coupled_samples, dense_residual_quality, &
                                           mermin_quality_ema, mermin_samples) RESULT(bridge_ready)
      REAL(KIND=dp), INTENT(IN)                          :: residual_delta, coupled_quality, &
                                                            coupled_quality_ema
      INTEGER, INTENT(IN)                                :: coupled_samples
      REAL(KIND=dp), INTENT(IN)                          :: dense_residual_quality, &
                                                            mermin_quality_ema
      INTEGER, INTENT(IN)                                :: mermin_samples
      LOGICAL                                            :: bridge_ready

      bridge_ready = .FALSE.
      IF (residual_delta > 1.0E-4_dp) RETURN
      IF (coupled_samples < 4) RETURN
      IF (mermin_samples < 4) RETURN
      IF (coupled_quality < 5.0E-2_dp) RETURN
      IF (coupled_quality_ema < 7.5E-2_dp) RETURN
      IF (dense_residual_quality < 2.5E-1_dp) RETURN
      IF (mermin_quality_ema < 8.0E-1_dp) RETURN

      bridge_ready = .TRUE.

   END FUNCTION dense_metric_tail_bridge_ready

! **************************************************************************************************
!> \brief Dense mirror of the residual-tail gauge rotation floor.
!> \param base_floor ...
!> \param gauge_floor ...
!> \param occ_gap_weight ...
!> \param step_gauge_quality ...
!> \return ...
! **************************************************************************************************
   FUNCTION residual_gauge_rotation_floor(base_floor, gauge_floor, occ_gap_weight, &
                                          step_gauge_quality) RESULT(rotation_floor)
      REAL(KIND=dp), INTENT(IN)                          :: base_floor, gauge_floor, occ_gap_weight, &
                                                            step_gauge_quality
      REAL(KIND=dp)                                      :: rotation_floor

      REAL(KIND=dp)                                      :: gap_fraction, residual_weight

      rotation_floor = MAX(0.0_dp, base_floor)
      IF (gauge_floor <= rotation_floor) RETURN

      gap_fraction = 1.0_dp - MIN(1.0_dp, MAX(0.0_dp, occ_gap_weight))
      residual_weight = 1.0_dp - MIN(1.0_dp, MAX(0.0_dp, step_gauge_quality))
      IF (gap_fraction <= EPSILON(gap_fraction) .OR. &
          residual_weight <= EPSILON(residual_weight)) RETURN

      rotation_floor = rotation_floor + gap_fraction*residual_weight*(gauge_floor - rotation_floor)

   END FUNCTION residual_gauge_rotation_floor

! **************************************************************************************************
!> \brief Dense mirror of the residual-tail gauge density-response filter.
!> \param density_response ...
!> \param density_projector_response ...
!> \param density_safe_response ...
!> \param occ_curv ...
!> \param occ_curv_max ...
!> \param response_ratio_gate ...
!> \param step_gauge_quality ...
!> \param density_projector_centered_response ...
!> \param density_safe_centered_response ...
!> \param density_response_alignment ...
!> \return ...
! **************************************************************************************************
   FUNCTION residual_gauge_density_response(density_response, density_projector_response, &
                                            density_safe_response, occ_curv, occ_curv_max, &
                                            response_ratio_gate, step_gauge_quality, &
                                            density_projector_centered_response, &
                                            density_safe_centered_response, &
                                            density_response_alignment) RESULT(filtered_response)
      REAL(KIND=dp), INTENT(IN) :: density_response, density_projector_response, &
         density_safe_response, occ_curv, occ_curv_max, response_ratio_gate, step_gauge_quality
      REAL(KIND=dp), INTENT(IN), OPTIONAL :: density_projector_centered_response, &
         density_safe_centered_response, density_response_alignment
      REAL(KIND=dp)                                      :: filtered_response

      LOGICAL                                            :: globally_supported_response, &
                                                            locally_supported_response
      REAL(KIND=dp) :: active_fraction, alignment, damping, projector_abs, residual_weight, &
         response_support, safe_abs, sign_projector_response, sign_safe_response

      filtered_response = density_response
      IF (occ_curv <= EPSILON(occ_curv)) RETURN
      IF (occ_curv_max <= EPSILON(occ_curv_max)) RETURN

      residual_weight = 1.0_dp - MIN(1.0_dp, MAX(0.0_dp, step_gauge_quality))
      IF (residual_weight <= EPSILON(residual_weight)) RETURN

      active_fraction = occ_curv/MAX(occ_curv, response_ratio_gate*occ_curv_max)
      active_fraction = MIN(1.0_dp, MAX(0.0_dp, active_fraction))
      sign_projector_response = density_projector_response
      sign_safe_response = density_safe_response
      IF (PRESENT(density_projector_centered_response)) &
         sign_projector_response = density_projector_centered_response
      IF (PRESENT(density_safe_centered_response)) &
         sign_safe_response = density_safe_centered_response
      alignment = 0.0_dp
      IF (PRESENT(density_response_alignment)) &
         alignment = MIN(1.0_dp, MAX(-1.0_dp, density_response_alignment))

      projector_abs = ABS(sign_projector_response)
      safe_abs = ABS(sign_safe_response)
      locally_supported_response = sign_projector_response*sign_safe_response > 0.0_dp .OR. &
                                   projector_abs <= EPSILON(projector_abs)
      globally_supported_response = alignment > 0.25_dp .AND. &
                                    safe_abs > EPSILON(safe_abs) .AND. &
                                    projector_abs <= alignment*safe_abs
      response_support = 0.0_dp
      IF (locally_supported_response .OR. globally_supported_response) THEN
         IF (projector_abs > EPSILON(projector_abs)) THEN
            response_support = MIN(1.0_dp, safe_abs/projector_abs)
         ELSE
            response_support = 1.0_dp
         END IF
      END IF

      damping = 1.0_dp - residual_weight*(1.0_dp - response_support)* &
                (1.0_dp - active_fraction)**2
      filtered_response = MAX(0.0_dp, MIN(1.0_dp, damping))*density_response

   END FUNCTION residual_gauge_density_response

! **************************************************************************************************
!> \brief Dense mirror of the trace-active Direct Hessian blend weight.
!> \param raw_trace_ema ...
!> \param active_trace_ema ...
!> \param history_weight ...
!> \param history_quality ...
!> \param samples ...
!> \param selected_quality ...
!> \param selected_samples ...
!> \return ...
! **************************************************************************************************
   FUNCTION dense_trace_active_mix_weight(raw_trace_ema, active_trace_ema, history_weight, &
                                          history_quality, samples, selected_quality, &
                                          selected_samples) RESULT(mix_weight)
      REAL(KIND=dp), INTENT(IN)                          :: raw_trace_ema, active_trace_ema, &
                                                            history_weight, history_quality
      INTEGER, INTENT(IN)                                :: samples
      REAL(KIND=dp), INTENT(IN)                          :: selected_quality
      INTEGER, INTENT(IN)                                :: selected_samples
      REAL(KIND=dp)                                      :: mix_weight

      REAL(KIND=dp)                                      :: active_fraction, gauge_fraction, &
                                                            quality_weight, raw_abs

      mix_weight = MIN(0.75_dp, MAX(0.0_dp, history_weight))
      IF (samples < 8) RETURN
      IF (history_quality < 0.10_dp) RETURN
      IF (selected_samples < 8) RETURN
      IF (selected_quality < 0.10_dp) RETURN

      raw_abs = ABS(raw_trace_ema)
      IF (raw_abs <= EPSILON(raw_abs)) RETURN

      active_fraction = ABS(active_trace_ema)/raw_abs
      active_fraction = MIN(1.0_dp, MAX(0.0_dp, active_fraction))
      gauge_fraction = 1.0_dp - active_fraction
      quality_weight = MIN(1.0_dp, MAX(0.0_dp, history_quality))

      mix_weight = MAX(mix_weight, 0.95_dp*gauge_fraction*quality_weight)
      mix_weight = MIN(0.95_dp, MAX(0.0_dp, mix_weight))

   END FUNCTION dense_trace_active_mix_weight

! **************************************************************************************************
!> \brief Dense mirror of the non-promoting active-basis Direct trace blend.
!> \param raw_trace_ema ...
!> \param active_trace_ema ...
!> \param history_quality ...
!> \param samples ...
!> \return ...
! **************************************************************************************************
   FUNCTION dense_trace_active_basis_weight(raw_trace_ema, active_trace_ema, history_quality, &
                                            samples) RESULT(basis_weight)
      REAL(KIND=dp), INTENT(IN)                          :: raw_trace_ema, active_trace_ema, &
                                                            history_quality
      INTEGER, INTENT(IN)                                :: samples
      REAL(KIND=dp)                                      :: basis_weight

      REAL(KIND=dp)                                      :: active_fraction, gauge_fraction, &
                                                            quality_weight, raw_abs

      basis_weight = 0.0_dp
      IF (samples < 8) RETURN
      IF (history_quality < 2.5E-2_dp) RETURN

      raw_abs = ABS(raw_trace_ema)
      IF (raw_abs <= EPSILON(raw_abs)) RETURN

      active_fraction = ABS(active_trace_ema)/raw_abs
      active_fraction = MIN(1.0_dp, MAX(0.0_dp, active_fraction))
      gauge_fraction = 1.0_dp - active_fraction
      quality_weight = SQRT(MIN(1.0_dp, MAX(0.0_dp, history_quality)))

      basis_weight = MIN(0.25_dp, gauge_fraction*quality_weight)

   END FUNCTION dense_trace_active_basis_weight

! **************************************************************************************************
!> \brief Dense mirror of trace-bearing Direct curvature component scaling.
!> \param trace_diag_ema ...
!> \param trace_free_ema ...
!> \param trace_diag_quality_ema ...
!> \param trace_free_quality_ema ...
!> \param samples ...
!> \return ...
! **************************************************************************************************
   FUNCTION dense_trace_diag_component_scale(trace_diag_ema, trace_free_ema, trace_diag_quality_ema, &
                                             trace_free_quality_ema, samples) RESULT(component_scale)
      REAL(KIND=dp), INTENT(IN)                          :: trace_diag_ema, trace_free_ema, &
                                                            trace_diag_quality_ema, &
                                                            trace_free_quality_ema
      INTEGER, INTENT(IN)                                :: samples
      REAL(KIND=dp)                                      :: component_scale

      REAL(KIND=dp)                                      :: balance_scale, diag_abs, diag_quality, &
                                                            free_abs, free_quality, quality_ratio, &
                                                            trace_balance, trace_norm

      component_scale = 1.0_dp
      IF (samples < 8) RETURN

      diag_abs = ABS(trace_diag_ema)
      IF (diag_abs <= EPSILON(diag_abs)) RETURN
      diag_quality = ABS(trace_diag_quality_ema)/diag_abs
      free_abs = ABS(trace_free_ema)
      IF (free_abs <= EPSILON(free_abs)) THEN
         component_scale = MIN(component_scale, SQRT(MIN(1.0_dp, MAX(0.10_dp, diag_quality))))
         RETURN
      END IF

      free_quality = ABS(trace_free_quality_ema)/free_abs
      IF (free_quality <= EPSILON(free_quality)) THEN
         component_scale = MIN(component_scale, SQRT(MIN(1.0_dp, MAX(0.10_dp, diag_quality))))
         RETURN
      END IF
      trace_norm = diag_abs + free_abs
      trace_balance = ABS(trace_diag_ema + trace_free_ema)/MAX(trace_norm, EPSILON(trace_norm))
      trace_balance = MIN(1.0_dp, MAX(0.0_dp, trace_balance))
      balance_scale = SQRT(MAX(0.10_dp, trace_balance))
      component_scale = MIN(component_scale, balance_scale)
      IF (diag_quality >= free_quality) RETURN

      quality_ratio = MAX(0.10_dp, diag_quality/free_quality)
      component_scale = MIN(component_scale, SQRT(MIN(1.0_dp, quality_ratio)))

   END FUNCTION dense_trace_diag_component_scale

! **************************************************************************************************
!> \brief Dense mirror of Direct trace scaling with quality-gated local compensation.
!> \param base_diag_avg ...
!> \param local_comp_avg ...
!> \param trace_scale ...
!> \param quality_floor ...
!> \return scaled trace diagonal average
! **************************************************************************************************
   FUNCTION dense_trace_scaled_diag_with_local_comp(base_diag_avg, local_comp_avg, &
                                                    trace_scale, quality_floor) RESULT(scaled_diag_avg)
      REAL(KIND=dp), INTENT(IN)                          :: base_diag_avg, local_comp_avg, &
                                                            trace_scale, quality_floor
      REAL(KIND=dp)                                      :: scaled_diag_avg

      REAL(KIND=dp)                                      :: separate_weight, trace_avg

      separate_weight = MAX(0.0_dp, trace_scale - SQRT(quality_floor))/ &
                        MAX(EPSILON(1.0_dp), 1.0_dp - SQRT(quality_floor))
      separate_weight = MIN(1.0_dp, separate_weight)
      trace_avg = (1.0_dp - separate_weight)*(base_diag_avg + local_comp_avg) + &
                  separate_weight*base_diag_avg
      scaled_diag_avg = base_diag_avg + local_comp_avg + (trace_scale - 1.0_dp)*trace_avg

   END FUNCTION dense_trace_scaled_diag_with_local_comp

! **************************************************************************************************
!> \brief Dense mirror of selected-history Direct step scaling.
!> \param predicted_sum ...
!> \param measured_sum ...
!> \param ema_quality ...
!> \param samples ...
!> \param min_samples ...
!> \param response_min_samples ...
!> \param scale_floor ...
!> \param max_quality ...
!> \return ...
! **************************************************************************************************
   FUNCTION dense_direct_history_step_scale(predicted_sum, measured_sum, ema_quality, samples, &
                                            min_samples, response_min_samples, scale_floor, &
                                            max_quality) RESULT(step_scale)
      REAL(KIND=dp), INTENT(IN)                          :: predicted_sum, measured_sum, ema_quality
      INTEGER, INTENT(IN)                                :: samples, min_samples, &
                                                            response_min_samples
      REAL(KIND=dp), INTENT(IN)                          :: scale_floor, max_quality
      REAL(KIND=dp)                                      :: step_scale

      INTEGER                                            :: maturity_samples
      REAL(KIND=dp)                                      :: history_quality, maturity, target_scale

      history_quality = 1.0_dp
      IF (samples >= min_samples) THEN
         IF (predicted_sum > EPSILON(predicted_sum) .AND. measured_sum > 0.0_dp) THEN
            history_quality = measured_sum/predicted_sum
            IF (ema_quality > 0.0_dp) history_quality = MAX(history_quality, ema_quality)
         ELSE IF (ema_quality > 0.0_dp) THEN
            history_quality = ema_quality
         END IF
      END IF
      history_quality = MAX(scale_floor, MIN(max_quality, history_quality))
      target_scale = SQRT(MIN(1.0_dp, history_quality))
      maturity_samples = MAX(1, 2*response_min_samples - min_samples + 1)
      maturity = REAL(MAX(0, samples - min_samples + 1), KIND=dp)/ &
                 REAL(maturity_samples, KIND=dp)
      maturity = MIN(1.0_dp, MAX(0.0_dp, maturity))
      step_scale = 1.0_dp - maturity*(1.0_dp - target_scale)

   END FUNCTION dense_direct_history_step_scale

! **************************************************************************************************
!> \brief Dense mirror of correlated Direct tail selection quality.
!> \param base_quality ...
!> \param history_quality ...
!> \return ...
! **************************************************************************************************
   FUNCTION dense_direct_tail_selection_quality(base_quality, history_quality) RESULT(quality)
      REAL(KIND=dp), INTENT(IN)                          :: base_quality, history_quality
      REAL(KIND=dp)                                      :: quality

      REAL(KIND=dp)                                      :: clipped_history

      quality = MAX(0.0_dp, MIN(1.0_dp, base_quality))
      IF (history_quality <= 0.0_dp) RETURN

      clipped_history = MAX(0.0_dp, MIN(1.0_dp, history_quality))
      quality = MIN(quality, clipped_history)

   END FUNCTION dense_direct_tail_selection_quality

! **************************************************************************************************
!> \brief Dense mirror of accepted-drop Direct history growth.
!> \param predicted_sum ...
!> \param measured_sum ...
!> \param ema_quality ...
!> \param samples ...
!> \param min_samples ...
!> \param response_min_samples ...
!> \param max_quality ...
!> \return ...
! **************************************************************************************************
   FUNCTION dense_direct_history_growth_scale(predicted_sum, measured_sum, ema_quality, samples, &
                                              min_samples, response_min_samples, max_quality) RESULT(step_scale)
      REAL(KIND=dp), INTENT(IN)                          :: predicted_sum, measured_sum, ema_quality
      INTEGER, INTENT(IN)                                :: samples, min_samples, &
                                                            response_min_samples
      REAL(KIND=dp), INTENT(IN)                          :: max_quality
      REAL(KIND=dp)                                      :: step_scale

      INTEGER                                            :: maturity_samples
      REAL(KIND=dp)                                      :: history_quality, maturity, target_scale

      step_scale = 1.0_dp
      IF (samples < min_samples) RETURN

      history_quality = 1.0_dp
      IF (predicted_sum > EPSILON(predicted_sum) .AND. measured_sum > 0.0_dp) THEN
         history_quality = measured_sum/predicted_sum
         IF (ema_quality > 0.0_dp) history_quality = MAX(history_quality, ema_quality)
      ELSE IF (ema_quality > 0.0_dp) THEN
         history_quality = ema_quality
      END IF
      history_quality = MAX(1.0_dp, MIN(max_quality, history_quality))
      target_scale = SQRT(history_quality)
      maturity_samples = MAX(1, 2*response_min_samples - min_samples + 1)
      maturity = REAL(MAX(0, samples - min_samples + 1), KIND=dp)/ &
                 REAL(maturity_samples, KIND=dp)
      maturity = MIN(1.0_dp, MAX(0.0_dp, maturity))
      step_scale = 1.0_dp + maturity*(target_scale - 1.0_dp)

   END FUNCTION dense_direct_history_growth_scale

! **************************************************************************************************
!> \brief Dense mirror of mixed-Direct accepted-drop consistency scaling.
!> \param drop_ratio ...
!> \param ema_quality ...
!> \param max_quality ...
!> \return ...
! **************************************************************************************************
   FUNCTION dense_direct_history_consistency_scale(drop_ratio, ema_quality, max_quality) &
      RESULT(consistency_scale)
      REAL(KIND=dp), INTENT(IN)                          :: drop_ratio, ema_quality, max_quality
      REAL(KIND=dp)                                      :: consistency_scale

      consistency_scale = 1.0_dp
      IF (drop_ratio > max_quality) &
         consistency_scale = MIN(consistency_scale, SQRT(max_quality/drop_ratio))
      IF (ema_quality > max_quality) &
         consistency_scale = MIN(consistency_scale, SQRT(max_quality/ema_quality))

   END FUNCTION dense_direct_history_consistency_scale

! **************************************************************************************************
!> \brief Dense mirror of mixed active/physical Direct promotion.
!> \param samples ...
!> \param min_samples ...
!> \param drop_ratio ...
!> \param history_quality ...
!> \param ema_quality ...
!> \param selected_quality ...
!> \param selected_samples ...
!> \param min_quality ...
!> \return ...
! **************************************************************************************************
   FUNCTION dense_mixed_direct_promoted(samples, min_samples, drop_ratio, history_quality, &
                                        ema_quality, selected_quality, selected_samples, &
                                        min_quality) RESULT(promoted)
      INTEGER, INTENT(IN)                                :: samples, min_samples
      REAL(KIND=dp), INTENT(IN)                          :: drop_ratio, history_quality, &
                                                            ema_quality, selected_quality
      INTEGER, INTENT(IN)                                :: selected_samples
      REAL(KIND=dp), INTENT(IN)                          :: min_quality
      LOGICAL                                            :: promoted

      INTEGER                                            :: promotion_samples, &
                                                            shadow_promotion_samples
      LOGICAL                                            :: selected_promoted, shadow_promoted

      promotion_samples = MAX(min_samples, 5*min_samples)
      shadow_promotion_samples = MAX(promotion_samples, 10*min_samples)
      selected_promoted = samples >= promotion_samples .AND. &
                          drop_ratio >= min_quality .AND. &
                          history_quality >= min_quality .AND. &
                          ema_quality >= 0.9_dp*min_quality .AND. &
                          selected_samples >= min_samples .AND. &
                          selected_quality >= min_quality
      shadow_promoted = samples >= shadow_promotion_samples .AND. &
                        drop_ratio >= 0.5_dp*min_quality .AND. &
                        history_quality >= 0.5_dp*min_quality .AND. &
                        ema_quality >= 0.5_dp*min_quality
      promoted = selected_promoted .OR. shadow_promoted

   END FUNCTION dense_mixed_direct_promoted

! **************************************************************************************************
!> \brief Dense mirror of the low-residual finite Projector-Hessian tail weight.
!> \param residual_delta ...
!> \param history_quality ...
!> \param samples ...
!> \param density_projector_fraction ...
!> \param projector_diag_avg ...
!> \param projector_delta_diag_avg ...
!> \return ...
! **************************************************************************************************
   FUNCTION dense_projector_coupled_tail_weight(residual_delta, history_quality, samples, &
                                                density_projector_fraction, projector_diag_avg, &
                                                projector_delta_diag_avg) RESULT(projector_weight)
      REAL(KIND=dp), INTENT(IN)                          :: residual_delta, history_quality
      INTEGER, INTENT(IN)                                :: samples
      REAL(KIND=dp), INTENT(IN)                          :: density_projector_fraction, &
                                                            projector_diag_avg, &
                                                            projector_delta_diag_avg
      REAL(KIND=dp)                                      :: projector_weight

      REAL(KIND=dp)                                      :: clipped_quality, curvature_fraction, &
                                                            density_fraction, diag_ratio, &
                                                            quality_weight

      IF (residual_delta > 1.0E-4_dp) THEN
         projector_weight = 1.0_dp
         RETURN
      END IF

      projector_weight = 0.0_dp
      IF (samples < 4) RETURN

      IF (history_quality <= 1.0_dp) THEN
         clipped_quality = MIN(1.0_dp, MAX(0.0_dp, history_quality))
         quality_weight = (1.0_dp - clipped_quality)/MAX(EPSILON(1.0_dp), 1.0_dp - 0.25_dp)
      ELSE
         clipped_quality = MIN(2.5_dp, history_quality)
         quality_weight = (clipped_quality - 1.0_dp)/MAX(EPSILON(1.0_dp), 2.5_dp - 1.0_dp)
      END IF
      quality_weight = MIN(1.0_dp, MAX(0.0_dp, quality_weight))
      diag_ratio = ABS(projector_delta_diag_avg)/MAX(ABS(projector_diag_avg), EPSILON(1.0_dp))
      curvature_fraction = diag_ratio/MAX(1.0_dp, 1.0_dp + diag_ratio)
      density_fraction = MAX(MIN(1.0_dp, MAX(0.0_dp, density_projector_fraction)), &
                             0.25_dp*curvature_fraction)
      projector_weight = 0.10_dp*quality_weight*density_fraction/ &
                         MAX(1.0_dp, 1.0_dp + diag_ratio)

   END FUNCTION dense_projector_coupled_tail_weight

! **************************************************************************************************
!> \brief Dense mirror of active density-response Projector-delta blend fraction.
!> \param projector_weight ...
!> \param density_history_weight ...
!> \return ...
! **************************************************************************************************
   FUNCTION dense_projector_tail_active_delta_weight(projector_weight, density_history_weight) RESULT(delta_weight)
      REAL(KIND=dp), INTENT(IN)                          :: projector_weight, density_history_weight
      REAL(KIND=dp)                                      :: delta_weight

      delta_weight = 0.0_dp
      IF (projector_weight <= EPSILON(projector_weight)) RETURN
      IF (density_history_weight <= EPSILON(density_history_weight)) RETURN

      delta_weight = MIN(1.0_dp, MAX(0.0_dp, density_history_weight)/ &
                         MAX(projector_weight, EPSILON(projector_weight)))

   END FUNCTION dense_projector_tail_active_delta_weight

! **************************************************************************************************
!> \brief Dense mirror of accepted-drop trust for low-curvature density response.
!> \param lowcurv_ema ...
!> \param lowcurv_quality_ema ...
!> \param amplification_ema ...
!> \param amplification_quality_ema ...
!> \param samples ...
!> \return ...
! **************************************************************************************************
   FUNCTION dense_density_audit_trust_quality(lowcurv_ema, lowcurv_quality_ema, &
                                              amplification_ema, amplification_quality_ema, &
                                              samples) RESULT(trust_quality)
      REAL(KIND=dp), INTENT(IN)                          :: lowcurv_ema, lowcurv_quality_ema, &
                                                            amplification_ema, &
                                                            amplification_quality_ema
      INTEGER, INTENT(IN)                                :: samples
      REAL(KIND=dp)                                      :: trust_quality

      INTEGER, PARAMETER                                 :: min_samples = 8
      REAL(KIND=dp), PARAMETER                           :: min_quality = 0.10_dp

      REAL(KIND=dp) :: amplification_quality, amplification_weight, component_quality, &
         lowcurv_quality, lowcurv_weight, raw_amplification_quality, raw_lowcurv_quality

      trust_quality = 1.0_dp
      IF (samples < min_samples) RETURN
      IF (lowcurv_ema <= 0.25_dp) RETURN

      lowcurv_quality = 0.0_dp
      amplification_quality = 0.0_dp
      raw_lowcurv_quality = 0.0_dp
      raw_amplification_quality = 0.0_dp
      IF (ABS(lowcurv_ema) > EPSILON(lowcurv_ema)) &
         raw_lowcurv_quality = ABS(lowcurv_quality_ema)/ABS(lowcurv_ema)
      IF (ABS(amplification_ema) > EPSILON(amplification_ema)) &
         raw_amplification_quality = ABS(amplification_quality_ema)/ABS(amplification_ema)
      IF (raw_lowcurv_quality > EPSILON(1.0_dp)) &
         lowcurv_quality = MIN(raw_lowcurv_quality, 1.0_dp/raw_lowcurv_quality)
      IF (raw_amplification_quality > EPSILON(1.0_dp)) &
         amplification_quality = MIN(raw_amplification_quality, 1.0_dp/raw_amplification_quality)
      component_quality = lowcurv_quality
      IF (amplification_quality > EPSILON(amplification_quality)) &
         component_quality = MIN(component_quality, amplification_quality)
      IF (component_quality <= EPSILON(component_quality)) RETURN

      component_quality = MIN(1.0_dp, MAX(min_quality, component_quality))
      lowcurv_weight = MIN(1.0_dp, MAX(0.0_dp, (lowcurv_ema - 0.25_dp)/0.75_dp))
      amplification_weight = MIN(1.0_dp, MAX(0.0_dp, LOG(MAX(1.0_dp, amplification_ema))/ &
                                             LOG(1.0E6_dp)))
      amplification_weight = MAX(amplification_weight, &
                                 MIN(1.0_dp, MAX(0.0_dp, &
                                                 LOG(MAX(1.0_dp, raw_lowcurv_quality, &
                                                         raw_amplification_quality))/LOG(1.0E6_dp))))
      trust_quality = 1.0_dp - lowcurv_weight*amplification_weight*(1.0_dp - component_quality)
      trust_quality = MIN(1.0_dp, MAX(min_quality, trust_quality))

   END FUNCTION dense_density_audit_trust_quality

! **************************************************************************************************
!> \brief Dense mirror of the Direct trace model-gap term.
!> \param trace_actual_quad ...
!> \param trace_model_quad ...
!> \param scale ...
!> \return ...
! **************************************************************************************************
   FUNCTION dense_direct_trace_model_gap_drop(trace_actual_quad, trace_model_quad, scale) RESULT(drop)
      REAL(KIND=dp), INTENT(IN)                          :: trace_actual_quad, trace_model_quad, &
                                                            scale
      REAL(KIND=dp)                                      :: drop

      drop = 0.5_dp*scale*scale*ABS(trace_actual_quad - trace_model_quad)

   END FUNCTION dense_direct_trace_model_gap_drop

! **************************************************************************************************
!> \brief Dense mirror of Direct trace-mismatch trust.
!> \param model_drop ...
!> \param trace_gap_drop ...
!> \param history_quality ...
!> \param samples ...
!> \return ...
! **************************************************************************************************
   FUNCTION dense_direct_trace_mismatch_trust(model_drop, trace_gap_drop, history_quality, samples) &
      RESULT(trust_quality)
      REAL(KIND=dp), INTENT(IN)                          :: model_drop, trace_gap_drop, &
                                                            history_quality
      INTEGER, INTENT(IN)                                :: samples
      REAL(KIND=dp)                                      :: trust_quality

      REAL(KIND=dp)                                      :: history_weight

      trust_quality = 1.0_dp
      IF (samples < 8) RETURN
      IF (model_drop <= EPSILON(model_drop)) RETURN
      IF (trace_gap_drop <= EPSILON(trace_gap_drop)) RETURN

      history_weight = 1.0_dp - MIN(1.0_dp, MAX(0.0_dp, history_quality))
      IF (history_weight <= EPSILON(history_weight)) RETURN
      trust_quality = model_drop/(model_drop + history_weight*trace_gap_drop)
      trust_quality = MIN(1.0_dp, MAX(0.0_dp, trust_quality))

   END FUNCTION dense_direct_trace_mismatch_trust

! **************************************************************************************************
!> \brief Project a density response onto the fixed-N Mermin descent direction.
!> \param response_vec ...
!> \param descent_vec ...
!> \param grad_vec ...
!> \param projected_vec ...
!> \param coeff ...
!> \param gdot ...
!> \param denom ...
!> \param flip_negative ...
! **************************************************************************************************
   SUBROUTINE dense_project_descent_response(response_vec, descent_vec, grad_vec, projected_vec, coeff, gdot, &
                                             denom, flip_negative)
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: response_vec, descent_vec, grad_vec
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: projected_vec
      REAL(KIND=dp), INTENT(OUT)                         :: coeff, gdot, denom
      LOGICAL, INTENT(IN), OPTIONAL                      :: flip_negative

      projected_vec(:) = 0.0_dp
      coeff = 0.0_dp
      gdot = DOT_PRODUCT(grad_vec, response_vec)
      denom = DOT_PRODUCT(grad_vec, descent_vec)
      IF (denom <= EPSILON(denom)) RETURN
      IF (gdot <= 0.0_dp) THEN
         IF (.NOT. PRESENT(flip_negative)) RETURN
         IF (.NOT. flip_negative) RETURN
         IF (gdot >= -EPSILON(denom)) RETURN
      END IF

      coeff = ABS(gdot)/denom
      projected_vec(:) = coeff*descent_vec(:)

   END SUBROUTINE dense_project_descent_response

! **************************************************************************************************
!> \brief Project a density response onto the fixed-N Mermin metric descent direction.
!> \param response_vec ...
!> \param descent_vec ...
!> \param grad_vec ...
!> \param occ_curv ...
!> \param projected_vec ...
!> \param coeff ...
!> \param gdot ...
!> \param denom ...
!> \param flip_negative ...
! **************************************************************************************************
   SUBROUTINE dense_project_metric_descent_response(response_vec, descent_vec, grad_vec, occ_curv, projected_vec, &
                                                    coeff, gdot, denom, flip_negative)
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: response_vec, descent_vec, grad_vec, &
                                                            occ_curv
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: projected_vec
      REAL(KIND=dp), INTENT(OUT)                         :: coeff, gdot, denom
      LOGICAL, INTENT(IN), OPTIONAL                      :: flip_negative

      REAL(KIND=dp)                                      :: grad_coeff, metric_coeff, metric_cross, &
                                                            metric_denom

      projected_vec(:) = 0.0_dp
      coeff = 0.0_dp
      gdot = DOT_PRODUCT(grad_vec, response_vec)
      denom = DOT_PRODUCT(grad_vec, descent_vec)
      IF (denom <= EPSILON(denom)) RETURN

      grad_coeff = 0.0_dp
      IF (gdot > 0.0_dp) THEN
         grad_coeff = gdot/denom
      ELSE IF (PRESENT(flip_negative)) THEN
         IF (flip_negative .AND. gdot < -EPSILON(denom)) grad_coeff = ABS(gdot)/denom
      END IF

      metric_coeff = 0.0_dp
      CALL dense_projected_metric(response_vec, descent_vec, occ_curv, metric_cross)
      CALL dense_projected_metric(descent_vec, descent_vec, occ_curv, metric_denom)
      IF (metric_denom > EPSILON(metric_denom)) THEN
         IF (metric_cross > 0.0_dp) THEN
            metric_coeff = metric_cross/metric_denom
         ELSE IF (PRESENT(flip_negative)) THEN
            IF (flip_negative .AND. metric_cross < -EPSILON(metric_denom)) &
               metric_coeff = ABS(metric_cross)/metric_denom
         END IF
      END IF

      coeff = MAX(grad_coeff, metric_coeff)
      IF (coeff <= 0.0_dp) RETURN
      projected_vec(:) = coeff*descent_vec(:)

   END SUBROUTINE dense_project_metric_descent_response

! **************************************************************************************************
!> \brief Project a density response onto the mixed fixed-N/tether Mermin descent direction.
!> \param response_vec ...
!> \param descent_vec ...
!> \param grad_vec ...
!> \param occ_curv ...
!> \param tether_curv ...
!> \param channel_nmo ...
!> \param e_offset ...
!> \param projected_vec ...
!> \param coeff ...
!> \param gdot ...
!> \param denom ...
!> \param flip_negative ...
! **************************************************************************************************
   SUBROUTINE dense_project_mixed_metric_descent_response(response_vec, descent_vec, grad_vec, occ_curv, &
                                                          tether_curv, channel_nmo, e_offset, projected_vec, &
                                                          coeff, gdot, denom, flip_negative)
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: response_vec, descent_vec, grad_vec, &
                                                            occ_curv, tether_curv
      INTEGER, DIMENSION(:), INTENT(IN)                  :: channel_nmo, e_offset
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: projected_vec
      REAL(KIND=dp), INTENT(OUT)                         :: coeff, gdot, denom
      LOGICAL, INTENT(IN), OPTIONAL                      :: flip_negative

      REAL(KIND=dp)                                      :: grad_coeff, metric_coeff, metric_cross, &
                                                            metric_denom

      projected_vec(:) = 0.0_dp
      coeff = 0.0_dp
      gdot = DOT_PRODUCT(grad_vec, response_vec)
      denom = DOT_PRODUCT(grad_vec, descent_vec)
      IF (denom <= EPSILON(denom)) RETURN

      grad_coeff = 0.0_dp
      IF (gdot > 0.0_dp) THEN
         grad_coeff = gdot/denom
      ELSE IF (PRESENT(flip_negative)) THEN
         IF (flip_negative .AND. gdot < -EPSILON(denom)) grad_coeff = ABS(gdot)/denom
      END IF

      metric_coeff = 0.0_dp
      CALL dense_projected_mixed_metric(response_vec, descent_vec, occ_curv, tether_curv, channel_nmo, &
                                        e_offset, metric_cross)
      CALL dense_projected_mixed_metric(descent_vec, descent_vec, occ_curv, tether_curv, channel_nmo, &
                                        e_offset, metric_denom)
      IF (metric_denom > EPSILON(metric_denom)) THEN
         IF (metric_cross > 0.0_dp) THEN
            metric_coeff = metric_cross/metric_denom
         ELSE IF (PRESENT(flip_negative)) THEN
            IF (flip_negative .AND. metric_cross < -EPSILON(metric_denom)) &
               metric_coeff = ABS(metric_cross)/metric_denom
         END IF
      END IF

      coeff = MAX(grad_coeff, metric_coeff)
      IF (coeff <= 0.0_dp) RETURN
      projected_vec(:) = coeff*descent_vec(:)

   END SUBROUTINE dense_project_mixed_metric_descent_response

! **************************************************************************************************
!> \brief Fixed-N dense Mermin metric used by the density-response projection fixture.
!> \param response_left ...
!> \param response_right ...
!> \param occ_curv ...
!> \param projected_term ...
! **************************************************************************************************
   SUBROUTINE dense_projected_metric(response_left, response_right, occ_curv, projected_term)
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: response_left, response_right, occ_curv
      REAL(KIND=dp), INTENT(OUT)                         :: projected_term

      REAL(KIND=dp)                                      :: left_shift, occ_sum, right_shift

      occ_sum = SUM(occ_curv(:))
      IF (occ_sum <= EPSILON(occ_sum)) THEN
         projected_term = 0.0_dp
         RETURN
      END IF

      left_shift = DOT_PRODUCT(occ_curv(:), response_left(:))/occ_sum
      right_shift = DOT_PRODUCT(occ_curv(:), response_right(:))/occ_sum
      projected_term = DOT_PRODUCT(occ_curv(:)*(response_left(:) - left_shift), &
                                   response_right(:) - right_shift)

   END SUBROUTINE dense_projected_metric

! **************************************************************************************************
!> \brief Mixed fixed-N/tether dense Mermin metric for response projection fixtures.
!> \param response_left ...
!> \param response_right ...
!> \param occ_curv ...
!> \param tether_curv ...
!> \param channel_nmo ...
!> \param e_offset ...
!> \param projected_term ...
! **************************************************************************************************
   SUBROUTINE dense_projected_mixed_metric(response_left, response_right, occ_curv, tether_curv, &
                                           channel_nmo, e_offset, projected_term)
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: response_left, response_right, occ_curv, &
                                                            tether_curv
      INTEGER, DIMENSION(:), INTENT(IN)                  :: channel_nmo, e_offset
      REAL(KIND=dp), INTENT(OUT)                         :: projected_term

      INTEGER                                            :: first_energy, local_channel, nmo
      REAL(KIND=dp)                                      :: left_shift, occ_sum, occ_term, &
                                                            right_shift, tether_term, weight_sum

      occ_term = DOT_PRODUCT(occ_curv(:)*response_left(:), response_right(:))
      occ_sum = SUM(occ_curv(:))
      IF (occ_sum > EPSILON(occ_sum)) THEN
         left_shift = DOT_PRODUCT(occ_curv(:), response_left(:))
         right_shift = DOT_PRODUCT(occ_curv(:), response_right(:))
         occ_term = occ_term - left_shift*right_shift/occ_sum
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

   END SUBROUTINE dense_projected_mixed_metric

! **************************************************************************************************
!> \brief Legacy two-pass dense response-gauge projection.
!> \param response_vec ...
!> \param occ_curv ...
!> \param tether_curv ...
!> \param channel_nmo ...
!> \param e_offset ...
! **************************************************************************************************
   SUBROUTINE dense_project_response_two_pass(response_vec, occ_curv, tether_curv, channel_nmo, e_offset)
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: response_vec
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: occ_curv, tether_curv
      INTEGER, DIMENSION(:), INTENT(IN)                  :: channel_nmo, e_offset

      INTEGER                                            :: first_energy, local_channel, nmo, pass
      REAL(KIND=dp)                                      :: occ_sum, response_shift, weight_sum

      occ_sum = SUM(occ_curv(:))
      DO pass = 1, 2
         IF (occ_sum > EPSILON(occ_sum)) THEN
            response_shift = DOT_PRODUCT(occ_curv(:), response_vec(:))
            response_vec(:) = response_vec(:) - response_shift/occ_sum
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

   END SUBROUTINE dense_project_response_two_pass

! **************************************************************************************************
!> \brief Dense mirror of the distributed curvature-metric response-gauge projection.
!> \param response_vec ...
!> \param occ_curv ...
!> \param tether_curv ...
!> \param channel_nmo ...
!> \param e_offset ...
!> \param projected ...
! **************************************************************************************************
   SUBROUTINE dense_project_response_metric_distributed(response_vec, occ_curv, tether_curv, channel_nmo, e_offset, &
                                                        projected)
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: response_vec
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: occ_curv, tether_curv
      INTEGER, DIMENSION(:), INTENT(IN)                  :: channel_nmo, e_offset
      LOGICAL, INTENT(OUT)                               :: projected

      INTEGER                                            :: first_energy, ie, local_channel, nmo
      REAL(KIND=dp) :: channel_alpha, channel_beta, channel_rhs, curv_max, denom, &
         global_alpha2_over_beta, global_alpha_rhs_over_beta, lambda_occ, max_lambda, occ_metric, &
         occ_rhs, occ_sum, total_curv
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: channel_lambda, inv_curv

      projected = .FALSE.
      curv_max = MAXVAL(occ_curv(:) + tether_curv(:))
      IF (curv_max <= EPSILON(curv_max)) RETURN

      ALLOCATE (channel_lambda(SIZE(channel_nmo)), inv_curv(SIZE(response_vec)))
      channel_lambda(:) = 0.0_dp
      DO ie = 1, SIZE(response_vec)
         total_curv = occ_curv(ie) + tether_curv(ie)
         IF (total_curv > EPSILON(curv_max)) THEN
            inv_curv(ie) = 1.0_dp/total_curv
         ELSE
            inv_curv(ie) = 0.0_dp
         END IF
      END DO

      occ_metric = DOT_PRODUCT(occ_curv(:)*inv_curv(:), occ_curv(:))
      occ_rhs = DOT_PRODUCT(occ_curv(:), response_vec(:))
      global_alpha2_over_beta = 0.0_dp
      global_alpha_rhs_over_beta = 0.0_dp
      DO local_channel = 1, SIZE(channel_nmo)
         nmo = channel_nmo(local_channel)
         IF (nmo <= 0) CYCLE
         first_energy = e_offset(local_channel)
         channel_alpha = DOT_PRODUCT(occ_curv(first_energy:first_energy + nmo - 1)* &
                                     inv_curv(first_energy:first_energy + nmo - 1), &
                                     tether_curv(first_energy:first_energy + nmo - 1))
         channel_beta = DOT_PRODUCT(tether_curv(first_energy:first_energy + nmo - 1)* &
                                    inv_curv(first_energy:first_energy + nmo - 1), &
                                    tether_curv(first_energy:first_energy + nmo - 1))
         channel_rhs = DOT_PRODUCT(tether_curv(first_energy:first_energy + nmo - 1), &
                                   response_vec(first_energy:first_energy + nmo - 1))
         IF (channel_beta > EPSILON(curv_max)) THEN
            global_alpha2_over_beta = global_alpha2_over_beta + channel_alpha*channel_alpha/channel_beta
            global_alpha_rhs_over_beta = global_alpha_rhs_over_beta + channel_alpha*channel_rhs/channel_beta
         END IF
      END DO

      occ_sum = SUM(occ_curv(:))
      lambda_occ = 0.0_dp
      IF (occ_sum > EPSILON(occ_sum)) THEN
         denom = occ_metric - global_alpha2_over_beta
         IF (denom <= SQRT(EPSILON(1.0_dp))*MAX(1.0_dp, occ_metric)) THEN
            DEALLOCATE (channel_lambda, inv_curv)
            RETURN
         END IF
         lambda_occ = (occ_rhs - global_alpha_rhs_over_beta)/denom
      END IF

      DO local_channel = 1, SIZE(channel_nmo)
         nmo = channel_nmo(local_channel)
         IF (nmo <= 0) CYCLE
         first_energy = e_offset(local_channel)
         channel_alpha = DOT_PRODUCT(occ_curv(first_energy:first_energy + nmo - 1)* &
                                     inv_curv(first_energy:first_energy + nmo - 1), &
                                     tether_curv(first_energy:first_energy + nmo - 1))
         channel_beta = DOT_PRODUCT(tether_curv(first_energy:first_energy + nmo - 1)* &
                                    inv_curv(first_energy:first_energy + nmo - 1), &
                                    tether_curv(first_energy:first_energy + nmo - 1))
         channel_rhs = DOT_PRODUCT(tether_curv(first_energy:first_energy + nmo - 1), &
                                   response_vec(first_energy:first_energy + nmo - 1))
         IF (channel_beta > EPSILON(curv_max)) &
            channel_lambda(local_channel) = (channel_rhs - channel_alpha*lambda_occ)/channel_beta
      END DO

      max_lambda = 0.0_dp
      DO local_channel = 1, SIZE(channel_nmo)
         max_lambda = MAX(max_lambda, ABS(channel_lambda(local_channel)))
      END DO
      IF (ABS(lambda_occ) > EPSILON(lambda_occ) .OR. max_lambda > EPSILON(curv_max)) THEN
         response_vec(:) = response_vec(:) - lambda_occ*inv_curv(:)*occ_curv(:)
         DO local_channel = 1, SIZE(channel_nmo)
            nmo = channel_nmo(local_channel)
            IF (nmo <= 0) CYCLE
            first_energy = e_offset(local_channel)
            response_vec(first_energy:first_energy + nmo - 1) = &
               response_vec(first_energy:first_energy + nmo - 1) - &
               channel_lambda(local_channel)*inv_curv(first_energy:first_energy + nmo - 1)* &
               tether_curv(first_energy:first_energy + nmo - 1)
         END DO
      END IF
      projected = .TRUE.
      DEALLOCATE (channel_lambda, inv_curv)

   END SUBROUTINE dense_project_response_metric_distributed

! **************************************************************************************************
!> \brief Dense curvature-metric response-gauge projection.
!> \param response_vec ...
!> \param occ_curv ...
!> \param tether_curv ...
!> \param channel_nmo ...
!> \param e_offset ...
!> \param projected ...
! **************************************************************************************************
   SUBROUTINE dense_project_response_metric(response_vec, occ_curv, tether_curv, channel_nmo, e_offset, projected)
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: response_vec
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: occ_curv, tether_curv
      INTEGER, DIMENSION(:), INTENT(IN)                  :: channel_nmo, e_offset
      LOGICAL, INTENT(OUT)                               :: projected

      INTEGER                                            :: first_energy, ie, info, local_channel, &
                                                            nconstraints, nmo, row
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ipiv
      REAL(KIND=dp)                                      :: curv_max, occ_sum, total_curv, weight_sum
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: inv_curv, lambda, rhs
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: constraint, gram

      projected = .FALSE.
      occ_sum = SUM(occ_curv(:))
      nconstraints = 0
      IF (occ_sum > EPSILON(occ_sum)) nconstraints = nconstraints + 1
      DO local_channel = 1, SIZE(channel_nmo)
         nmo = channel_nmo(local_channel)
         IF (nmo <= 0) CYCLE
         first_energy = e_offset(local_channel)
         weight_sum = SUM(tether_curv(first_energy:first_energy + nmo - 1))
         IF (weight_sum > EPSILON(weight_sum)) nconstraints = nconstraints + 1
      END DO
      IF (nconstraints <= 0) RETURN

      ALLOCATE (constraint(nconstraints, SIZE(response_vec)), gram(nconstraints, nconstraints), &
                inv_curv(SIZE(response_vec)), ipiv(nconstraints), lambda(nconstraints), rhs(nconstraints))
      constraint(:, :) = 0.0_dp
      row = 0
      IF (occ_sum > EPSILON(occ_sum)) THEN
         row = row + 1
         constraint(row, :) = occ_curv(:)
      END IF
      DO local_channel = 1, SIZE(channel_nmo)
         nmo = channel_nmo(local_channel)
         IF (nmo <= 0) CYCLE
         first_energy = e_offset(local_channel)
         weight_sum = SUM(tether_curv(first_energy:first_energy + nmo - 1))
         IF (weight_sum <= EPSILON(weight_sum)) CYCLE
         row = row + 1
         constraint(row, first_energy:first_energy + nmo - 1) = &
            tether_curv(first_energy:first_energy + nmo - 1)
      END DO

      curv_max = MAXVAL(occ_curv(:) + tether_curv(:))
      IF (curv_max <= EPSILON(curv_max)) THEN
         DEALLOCATE (constraint, gram, inv_curv, ipiv, lambda, rhs)
         RETURN
      END IF
      DO ie = 1, SIZE(response_vec)
         total_curv = occ_curv(ie) + tether_curv(ie)
         IF (total_curv > EPSILON(curv_max)) THEN
            inv_curv(ie) = 1.0_dp/total_curv
         ELSE
            inv_curv(ie) = 0.0_dp
         END IF
      END DO

      DO row = 1, nconstraints
         rhs(row) = DOT_PRODUCT(constraint(row, :), response_vec(:))
         DO ie = 1, nconstraints
            gram(row, ie) = DOT_PRODUCT(constraint(row, :)*inv_curv(:), constraint(ie, :))
         END DO
      END DO
      lambda(:) = rhs(:)
      CALL dgesv(nconstraints, 1, gram, nconstraints, ipiv, lambda, nconstraints, info)
      IF (info == 0) THEN
         DO row = 1, nconstraints
            response_vec(:) = response_vec(:) - lambda(row)*inv_curv(:)*constraint(row, :)
         END DO
         projected = .TRUE.
      END IF
      DEALLOCATE (constraint, gram, inv_curv, ipiv, lambda, rhs)

   END SUBROUTINE dense_project_response_metric

! **************************************************************************************************
!> \brief Maximum channel-local tether residual for a dense response vector.
!> \param response_vec ...
!> \param tether_curv ...
!> \param channel_nmo ...
!> \param e_offset ...
!> \return ...
! **************************************************************************************************
   FUNCTION dense_response_tether_residual(response_vec, tether_curv, channel_nmo, e_offset) RESULT(residual)
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: response_vec, tether_curv
      INTEGER, DIMENSION(:), INTENT(IN)                  :: channel_nmo, e_offset
      REAL(KIND=dp)                                      :: residual

      INTEGER                                            :: first_energy, local_channel, nmo

      residual = 0.0_dp
      DO local_channel = 1, SIZE(channel_nmo)
         nmo = channel_nmo(local_channel)
         IF (nmo <= 0) CYCLE
         first_energy = e_offset(local_channel)
         residual = MAX(residual, ABS(DOT_PRODUCT(tether_curv(first_energy:first_energy + nmo - 1), &
                                                  response_vec(first_energy:first_energy + nmo - 1))))
      END DO

   END FUNCTION dense_response_tether_residual

! **************************************************************************************************
!> \brief Check monotone projection for noisy occupation slopes before divided differences.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_monotone_occupation_slope_projection(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 5

      INTEGER                                            :: i, j
      REAL(KIND=dp)                                      :: coeff, err, positive_slope
      REAL(KIND=dp), DIMENSION(n)                        :: evals, expected, monotone_input, &
                                                            noisy_input, projected, &
                                                            projected_monotone

      noisy_input = [1.40_dp, 0.82_dp, 0.87_dp, 0.35_dp, 0.36_dp]
      expected = [1.40_dp, 0.845_dp, 0.845_dp, 0.355_dp, 0.355_dp]
      monotone_input = [1.40_dp, 0.92_dp, 0.54_dp, 0.31_dp, 0.08_dp]
      evals = [-0.42_dp, -0.11_dp, 0.02_dp, 0.19_dp, 0.64_dp]

      CALL dense_monotone_occupation_slopes(noisy_input, SQRT(EPSILON(1.0_dp)), projected)
      CALL dense_monotone_occupation_slopes(monotone_input, SQRT(EPSILON(1.0_dp)), projected_monotone)

      err = MAXVAL(ABS(projected - expected))
      err = MAX(err, ABS(SUM(projected) - SUM(noisy_input)))
      err = MAX(err, MAXVAL(ABS(projected_monotone - monotone_input)))
      positive_slope = 0.0_dp
      DO j = 1, n
         DO i = 1, n
            IF (i == j) CYCLE
            coeff = (projected(i) - projected(j))/(evals(i) - evals(j))
            positive_slope = MAX(positive_slope, coeff)
         END DO
      END DO

      IF (err > tol .OR. positive_slope > tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,2(1X,ES12.5))') TRIM(label), "failed with error/positive", &
            err, positive_slope
      ELSE
         WRITE (*, '(A,1X,A,2(1X,ES12.5))') TRIM(label), "OK with error/positive", &
            err, positive_slope
      END IF

   END SUBROUTINE check_monotone_occupation_slope_projection

! **************************************************************************************************
!> \brief Check the production-style monotone occupation preprojection used by k-point OT.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_history_safe_occ_projector_derivative(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 4
      REAL(KIND=dp), PARAMETER                           :: projector_weight = 0.37_dp

      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: chc, dchc, docc_raw, docc_ref, docc_safe
      INTEGER                                            :: i
      REAL(KIND=dp)                                      :: active, err
      REAL(KIND=dp), DIMENSION(n)                        :: occ_mono, occ_noisy, occ_projected

      chc(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      chc(1, 1) = CMPLX(-0.42_dp, 0.0_dp, KIND=dp)
      chc(2, 2) = CMPLX(-0.11_dp, 0.0_dp, KIND=dp)
      chc(3, 3) = CMPLX(0.02_dp, 0.0_dp, KIND=dp)
      chc(4, 4) = CMPLX(0.47_dp, 0.0_dp, KIND=dp)

      dchc(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      dchc(1, 2) = CMPLX(0.17_dp, -0.05_dp, KIND=dp)
      dchc(1, 3) = CMPLX(-0.09_dp, 0.03_dp, KIND=dp)
      dchc(2, 3) = CMPLX(0.21_dp, -0.13_dp, KIND=dp)
      dchc(2, 4) = CMPLX(-0.04_dp, -0.02_dp, KIND=dp)
      dchc(3, 4) = CMPLX(0.08_dp, 0.06_dp, KIND=dp)
      DO i = 2, n
         dchc(i, 1:i - 1) = CONJG(dchc(1:i - 1, i))
      END DO

      occ_mono(:) = [1.40_dp, 0.92_dp, 0.54_dp, 0.11_dp]
      CALL dense_history_safe_occ_projector_derivative_chc(chc, dchc, occ_mono, projector_weight, docc_safe)
      CALL dense_spectral_occ_projector_derivative_chc(chc, dchc, occ_mono, docc_ref)
      err = MAXVAL(ABS(docc_safe(:, :) - projector_weight*docc_ref(:, :)))

      occ_noisy(:) = [1.40_dp, 0.82_dp, 0.87_dp, 0.35_dp]
      CALL dense_monotone_occupation_slopes(occ_noisy, SQRT(EPSILON(1.0_dp)), occ_projected)
      CALL dense_history_safe_occ_projector_derivative_chc(chc, dchc, occ_noisy, projector_weight, docc_safe)
      CALL dense_spectral_occ_projector_derivative_chc(chc, dchc, occ_projected, docc_ref)
      CALL dense_spectral_occ_projector_derivative_chc(chc, dchc, occ_noisy, docc_raw)
      err = MAX(err, MAXVAL(ABS(docc_safe(:, :) - projector_weight*docc_ref(:, :))))
      active = MAXVAL(ABS(docc_safe(:, :) - projector_weight*docc_raw(:, :)))

      IF (err > tol .OR. active < 1.0E-4_dp) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,2(1X,ES12.5))') TRIM(label), "failed with error/active", err, active
      ELSE
         WRITE (*, '(A,1X,A,2(1X,ES12.5))') TRIM(label), "OK with error/active", err, active
      END IF

   END SUBROUTINE check_history_safe_occ_projector_derivative

! **************************************************************************************************
!> \brief Check that history-safe dP removes gauge noise but keeps physical response.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_history_safe_occ_projector_gauge_physical_split(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 4
      REAL(KIND=dp), PARAMETER                           :: active_tol = 1.0E-4_dp, &
                                                            fd_step = 1.0E-7_dp, &
                                                            gauge_gap = 5.0E-5_dp, &
                                                            projector_weight = 0.37_dp

      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: chc, chc_eig, dchc, dchc_eig, docc_fd, &
                                                            docc_raw, docc_ref, docc_safe, &
                                                            occ_mat_m, occ_mat_p, q
      INTEGER                                            :: i
      REAL(KIND=dp) :: err, gauge_err, mixed_err, mixed_safe, occ_tol, physical_active, &
         physical_err, physical_near_err, raw_gauge, raw_mixed, raw_transitive, safe_gauge, &
         transitive_err, value_err, value_raw
      REAL(KIND=dp), DIMENSION(n) :: evals_ref, occ_cluster, occ_curv, occ_curv_mixed, &
         occ_filtered, occ_mixed, occ_near_physical, occ_plateau, occ_scale
      REAL(KIND=dp), DIMENSION(n, n)                     :: x_im, x_re

      chc_eig(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      chc_eig(1, 1) = CMPLX(-0.010000_dp, 0.0_dp, KIND=dp)
      chc_eig(2, 2) = CMPLX(-0.009982_dp, 0.0_dp, KIND=dp)
      chc_eig(3, 3) = CMPLX(0.23_dp, 0.0_dp, KIND=dp)
      chc_eig(4, 4) = CMPLX(0.51_dp, 0.0_dp, KIND=dp)

      x_re(:, :) = 0.0_dp
      x_re(1, 2) = 0.13_dp
      x_re(1, 3) = -0.06_dp
      x_re(2, 4) = 0.04_dp
      x_re(3, 4) = -0.05_dp
      DO i = 2, n
         x_re(i, 1:i - 1) = -x_re(1:i - 1, i)
      END DO
      x_im(:, :) = 0.0_dp
      x_im(1, 1) = 0.02_dp
      x_im(2, 2) = -0.03_dp
      x_im(3, 3) = 0.01_dp
      x_im(4, 4) = -0.01_dp
      x_im(1, 2) = -0.04_dp
      x_im(1, 4) = 0.03_dp
      x_im(2, 3) = 0.02_dp
      DO i = 2, n
         x_im(i, 1:i - 1) = x_im(1:i - 1, i)
      END DO
      CALL dense_rotation_u(x_re, x_im, q)
      chc(:, :) = MATMUL(q, MATMUL(chc_eig, CONJG(TRANSPOSE(q))))

      occ_scale(:) = [1.70_dp, 0.52_dp, 0.24_dp, 0.2316_dp]
      occ_curv(:) = [0.038_dp, 0.035_dp, 0.018_dp, 0.006_dp]
      occ_plateau(:) = occ_scale(:)
      occ_plateau(1:2) = SUM(occ_scale(1:2))/2.0_dp
      DO i = 1, n
         evals_ref(i) = REAL(chc_eig(i, i), KIND=dp)
      END DO

      dchc_eig(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      dchc_eig(1, 2) = CMPLX(0.21_dp, -0.17_dp, KIND=dp)
      dchc_eig(2, 1) = CONJG(dchc_eig(1, 2))
      dchc(:, :) = MATMUL(q, MATMUL(dchc_eig, CONJG(TRANSPOSE(q))))
      CALL dense_spectral_occ_projector_derivative_chc(chc, dchc, occ_scale, docc_raw)
      CALL dense_history_safe_occ_projector_derivative_chc(chc, dchc, occ_scale, projector_weight, &
                                                           docc_safe, occ_curv, gauge_gap)
      occ_filtered(:) = occ_scale(:)
      occ_tol = SQRT(EPSILON(1.0_dp))*MAX(1.0_dp, MAXVAL(ABS(occ_filtered(:))))
      CALL dense_near_degenerate_occupation_plateau(evals_ref, occ_curv, gauge_gap, occ_tol, occ_filtered)
      CALL dense_spectral_occ_projector_from_chc(chc + fd_step*dchc, occ_filtered, occ_mat_p)
      CALL dense_spectral_occ_projector_from_chc(chc - fd_step*dchc, occ_filtered, occ_mat_m)
      docc_fd(:, :) = projector_weight*(occ_mat_p(:, :) - occ_mat_m(:, :))/(2.0_dp*fd_step)
      gauge_err = MAXVAL(ABS(docc_safe(:, :) - docc_fd(:, :)))
      raw_gauge = MAXVAL(ABS(projector_weight*docc_raw(:, :)))
      safe_gauge = MAXVAL(ABS(docc_safe(:, :)))
      CALL dense_history_safe_occ_projector_from_chc(chc + fd_step*dchc, occ_scale, projector_weight, &
                                                     occ_mat_p, occ_curv, gauge_gap)
      CALL dense_history_safe_occ_projector_from_chc(chc - fd_step*dchc, occ_scale, projector_weight, &
                                                     occ_mat_m, occ_curv, gauge_gap)
      docc_fd(:, :) = (occ_mat_p(:, :) - occ_mat_m(:, :))/(2.0_dp*fd_step)
      value_err = MAXVAL(ABS(docc_safe(:, :) - docc_fd(:, :)))
      CALL dense_history_safe_occ_projector_from_chc(chc + fd_step*dchc, occ_scale, projector_weight, &
                                                     occ_mat_p)
      CALL dense_history_safe_occ_projector_from_chc(chc - fd_step*dchc, occ_scale, projector_weight, &
                                                     occ_mat_m)
      docc_fd(:, :) = (occ_mat_p(:, :) - occ_mat_m(:, :))/(2.0_dp*fd_step)
      value_raw = MAXVAL(ABS(docc_fd(:, :)))

      occ_curv_mixed(:) = occ_curv(:)
      occ_curv_mixed(2) = 1.0E-5_dp
      occ_mixed(:) = occ_scale(:)
      occ_mixed(1) = 1.10_dp
      occ_mixed(2) = occ_mixed(1) - 1.0E-6_dp
      CALL dense_history_safe_occ_projector_derivative_chc(chc, dchc, occ_mixed, projector_weight, &
                                                           docc_safe, occ_curv_mixed, gauge_gap)
      CALL dense_spectral_occ_projector_derivative_chc(chc, dchc, occ_mixed, docc_raw)
      occ_filtered(:) = occ_mixed(:)
      occ_tol = SQRT(EPSILON(1.0_dp))*MAX(1.0_dp, MAXVAL(ABS(occ_filtered(:))))
      CALL dense_near_degenerate_occupation_plateau(evals_ref, occ_curv_mixed, gauge_gap, occ_tol, occ_filtered)
      CALL dense_spectral_occ_projector_from_chc(chc + fd_step*dchc, occ_filtered, occ_mat_p)
      CALL dense_spectral_occ_projector_from_chc(chc - fd_step*dchc, occ_filtered, occ_mat_m)
      docc_fd(:, :) = projector_weight*(occ_mat_p(:, :) - occ_mat_m(:, :))/(2.0_dp*fd_step)
      mixed_err = MAXVAL(ABS(docc_safe(:, :) - docc_fd(:, :)))
      mixed_safe = MAXVAL(ABS(docc_safe(:, :)))
      raw_mixed = MAXVAL(ABS(projector_weight*docc_raw(:, :)))

      dchc_eig(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      dchc_eig(3, 4) = CMPLX(-0.16_dp, 0.11_dp, KIND=dp)
      dchc_eig(4, 3) = CONJG(dchc_eig(3, 4))
      dchc(:, :) = MATMUL(q, MATMUL(dchc_eig, CONJG(TRANSPOSE(q))))
      CALL dense_history_safe_occ_projector_derivative_chc(chc, dchc, occ_scale, projector_weight, &
                                                           docc_safe, occ_curv, gauge_gap)
      CALL dense_spectral_occ_projector_from_chc(chc + fd_step*dchc, occ_plateau, occ_mat_p)
      CALL dense_spectral_occ_projector_from_chc(chc - fd_step*dchc, occ_plateau, occ_mat_m)
      docc_fd(:, :) = projector_weight*(occ_mat_p(:, :) - occ_mat_m(:, :))/(2.0_dp*fd_step)
      physical_err = MAXVAL(ABS(docc_safe(:, :) - docc_fd(:, :)))
      physical_active = MAXVAL(ABS(docc_safe(:, :)))

      occ_near_physical(:) = occ_scale(:)
      occ_near_physical(1) = 1.10_dp
      occ_near_physical(2) = occ_near_physical(1) - &
                             dense_shared_occ_slope(occ_curv(1), occ_curv(2))* &
                             (REAL(chc_eig(2, 2), KIND=dp) - REAL(chc_eig(1, 1), KIND=dp))
      dchc_eig(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      dchc_eig(1, 2) = CMPLX(-0.12_dp, 0.08_dp, KIND=dp)
      dchc_eig(2, 1) = CONJG(dchc_eig(1, 2))
      dchc(:, :) = MATMUL(q, MATMUL(dchc_eig, CONJG(TRANSPOSE(q))))
      CALL dense_history_safe_occ_projector_derivative_chc(chc, dchc, occ_near_physical, projector_weight, &
                                                           docc_safe, occ_curv, gauge_gap)
      CALL dense_spectral_occ_projector_derivative_chc(chc, dchc, occ_near_physical, docc_raw)
      physical_near_err = MAXVAL(ABS(docc_safe(:, :) - projector_weight*docc_raw(:, :)))

      chc_eig(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      chc_eig(1, 1) = CMPLX(-0.010000_dp, 0.0_dp, KIND=dp)
      chc_eig(2, 2) = CMPLX(-0.009970_dp, 0.0_dp, KIND=dp)
      chc_eig(3, 3) = CMPLX(-0.009940_dp, 0.0_dp, KIND=dp)
      chc_eig(4, 4) = CMPLX(0.51_dp, 0.0_dp, KIND=dp)
      chc(:, :) = MATMUL(q, MATMUL(chc_eig, CONJG(TRANSPOSE(q))))
      DO i = 1, n
         evals_ref(i) = REAL(chc_eig(i, i), KIND=dp)
      END DO
      occ_cluster(:) = occ_scale(:)
      occ_tol = SQRT(EPSILON(1.0_dp))*MAX(1.0_dp, MAXVAL(ABS(occ_cluster(:))))
      CALL dense_near_degenerate_occupation_plateau(evals_ref, occ_curv, gauge_gap, occ_tol, occ_cluster)
      dchc_eig(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      dchc_eig(1, 3) = CMPLX(0.18_dp, -0.09_dp, KIND=dp)
      dchc_eig(3, 1) = CONJG(dchc_eig(1, 3))
      dchc(:, :) = MATMUL(q, MATMUL(dchc_eig, CONJG(TRANSPOSE(q))))
      CALL dense_history_safe_occ_projector_derivative_chc(chc, dchc, occ_scale, projector_weight, &
                                                           docc_safe, occ_curv, gauge_gap)
      CALL dense_spectral_occ_projector_derivative_chc(chc, dchc, occ_scale, docc_raw)
      raw_transitive = MAXVAL(ABS(projector_weight*docc_raw(:, :)))
      CALL dense_spectral_occ_projector_derivative_chc(chc, dchc, occ_cluster, docc_ref, &
                                                       occ_curv, gauge_gap)
      docc_ref(:, :) = projector_weight*docc_ref(:, :)
      transitive_err = MAXVAL(ABS(docc_safe(:, :) - docc_ref(:, :)))

      err = MAX(gauge_err, mixed_err, physical_err, physical_near_err, transitive_err, value_err)
      IF (raw_gauge < active_tol .OR. raw_mixed < active_tol .OR. raw_transitive < active_tol .OR. &
          safe_gauge >= 0.20_dp*raw_gauge .OR. mixed_safe >= 0.20_dp*raw_mixed .OR. &
          physical_active < active_tol .OR. value_raw < active_tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,12(1X,ES12.5))') TRIM(label), &
            "failed with error/gerr/raw/safe/merr/raw_mixed/mixed/physical/near/transitive/value/raw_value", &
            err, gauge_err, raw_gauge, safe_gauge, mixed_err, raw_mixed, mixed_safe, &
            physical_active, physical_near_err, transitive_err, value_err, value_raw
      ELSE
         WRITE (*, '(A,1X,A,12(1X,ES12.5))') TRIM(label), &
            "OK with error/gerr/raw/safe/merr/raw_mixed/mixed/physical/near/transitive/value/raw_value", &
            err, gauge_err, raw_gauge, safe_gauge, mixed_err, raw_mixed, mixed_safe, &
            physical_active, physical_near_err, transitive_err, value_err, value_raw
      END IF

   END SUBROUTINE check_history_safe_occ_projector_gauge_physical_split

! **************************************************************************************************
!> \brief Check that finite-difference history-safe dP keeps physical near-Fermi slopes.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_history_safe_occ_projector_physical_fd(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 4
      REAL(KIND=dp), PARAMETER                           :: active_tol = 1.0E-4_dp, &
                                                            fd_step = 1.0E-7_dp, &
                                                            gauge_gap = 5.0E-5_dp, &
                                                            projector_weight = 0.37_dp

      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: chc, chc_eig, dchc, dchc_eig, docc_fd, &
                                                            docc_raw, docc_safe, occ_mat_m, &
                                                            occ_mat_p, q
      INTEGER                                            :: i
      REAL(KIND=dp)                                      :: err, gap, physical_active, &
                                                            physical_raw_err, physical_slope, &
                                                            raw_gauge, safe_gauge
      REAL(KIND=dp), DIMENSION(n)                        :: occ_curv, occ_gauge, occ_physical
      REAL(KIND=dp), DIMENSION(n, n)                     :: x_im, x_re

      chc_eig(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      chc_eig(1, 1) = CMPLX(-0.010000_dp, 0.0_dp, KIND=dp)
      chc_eig(2, 2) = CMPLX(-0.009982_dp, 0.0_dp, KIND=dp)
      chc_eig(3, 3) = CMPLX(0.23_dp, 0.0_dp, KIND=dp)
      chc_eig(4, 4) = CMPLX(0.51_dp, 0.0_dp, KIND=dp)

      x_re(:, :) = 0.0_dp
      x_re(1, 2) = 0.13_dp
      x_re(1, 3) = -0.06_dp
      x_re(2, 4) = 0.04_dp
      x_re(3, 4) = -0.05_dp
      DO i = 2, n
         x_re(i, 1:i - 1) = -x_re(1:i - 1, i)
      END DO
      x_im(:, :) = 0.0_dp
      x_im(1, 1) = 0.02_dp
      x_im(2, 2) = -0.03_dp
      x_im(3, 3) = 0.01_dp
      x_im(4, 4) = -0.01_dp
      x_im(1, 2) = -0.04_dp
      x_im(1, 4) = 0.03_dp
      x_im(2, 3) = 0.02_dp
      DO i = 2, n
         x_im(i, 1:i - 1) = x_im(1:i - 1, i)
      END DO
      CALL dense_rotation_u(x_re, x_im, q)
      chc(:, :) = MATMUL(q, MATMUL(chc_eig, CONJG(TRANSPOSE(q))))

      dchc_eig(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      dchc_eig(1, 2) = CMPLX(0.21_dp, -0.17_dp, KIND=dp)
      dchc_eig(2, 1) = CONJG(dchc_eig(1, 2))
      dchc(:, :) = MATMUL(q, MATMUL(dchc_eig, CONJG(TRANSPOSE(q))))

      occ_curv(:) = [0.038_dp, 0.035_dp, 0.018_dp, 0.006_dp]
      occ_gauge(:) = [1.10_dp, 0.50_dp, 0.24_dp, 0.08_dp]

      CALL dense_history_safe_occ_projector_derivative_chc(chc, dchc, occ_gauge, projector_weight, &
                                                           docc_safe, occ_curv, gauge_gap)
      CALL dense_spectral_occ_projector_derivative_chc(chc, dchc, occ_gauge, docc_raw)
      CALL dense_history_safe_occ_projector_from_chc(chc + fd_step*dchc, occ_gauge, projector_weight, &
                                                     occ_mat_p, occ_curv, gauge_gap)
      CALL dense_history_safe_occ_projector_from_chc(chc - fd_step*dchc, occ_gauge, projector_weight, &
                                                     occ_mat_m, occ_curv, gauge_gap)
      docc_fd(:, :) = (occ_mat_p(:, :) - occ_mat_m(:, :))/(2.0_dp*fd_step)
      err = MAXVAL(ABS(docc_safe(:, :) - docc_fd(:, :)))
      raw_gauge = MAXVAL(ABS(projector_weight*docc_raw(:, :)))
      safe_gauge = MAXVAL(ABS(docc_safe(:, :)))

      gap = REAL(chc_eig(2, 2), KIND=dp) - REAL(chc_eig(1, 1), KIND=dp)
      physical_slope = -0.5_dp*dense_shared_occ_slope(occ_curv(1), occ_curv(2))
      occ_physical(:) = [1.10_dp, 0.0_dp, 0.24_dp, 0.08_dp]
      occ_physical(2) = occ_physical(1) + physical_slope*gap
      CALL dense_history_safe_occ_projector_derivative_chc(chc, dchc, occ_physical, projector_weight, &
                                                           docc_safe, occ_curv, gauge_gap)
      CALL dense_spectral_occ_projector_derivative_chc(chc, dchc, occ_physical, docc_raw)
      CALL dense_history_safe_occ_projector_from_chc(chc + fd_step*dchc, occ_physical, &
                                                     projector_weight, occ_mat_p, occ_curv, gauge_gap)
      CALL dense_history_safe_occ_projector_from_chc(chc - fd_step*dchc, occ_physical, &
                                                     projector_weight, occ_mat_m, occ_curv, gauge_gap)
      docc_fd(:, :) = (occ_mat_p(:, :) - occ_mat_m(:, :))/(2.0_dp*fd_step)
      err = MAX(err, MAXVAL(ABS(docc_safe(:, :) - docc_fd(:, :))))
      physical_raw_err = MAXVAL(ABS(docc_safe(:, :) - projector_weight*docc_raw(:, :)))
      physical_active = MAXVAL(ABS(docc_safe(:, :)))

      IF (raw_gauge < active_tol .OR. safe_gauge >= 0.20_dp*raw_gauge .OR. &
          physical_active < active_tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,5(1X,ES12.5))') TRIM(label), &
            "failed with error/rawerr/rawgauge/safegauge/physical", err, physical_raw_err, &
            raw_gauge, safe_gauge, physical_active
      ELSE
         WRITE (*, '(A,1X,A,5(1X,ES12.5))') TRIM(label), &
            "OK with error/rawerr/rawgauge/safegauge/physical", err, physical_raw_err, &
            raw_gauge, safe_gauge, physical_active
      END IF

   END SUBROUTINE check_history_safe_occ_projector_physical_fd

! **************************************************************************************************
!> \brief Check that the finite-difference projector value caps smooth gauge slopes.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_history_safe_occ_projector_slope_cap_fd(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 3
      REAL(KIND=dp), PARAMETER                           :: active_tol = 1.0E-4_dp, &
                                                            fd_step = 1.0E-7_dp, &
                                                            gauge_gap = 5.0E-5_dp, &
                                                            projector_weight = 0.37_dp

      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: chc, dchc, docc_fd, docc_raw, docc_ref, &
                                                            docc_safe, occ_mat_m, occ_mat_p
      INTEGER                                            :: i
      REAL(KIND=dp)                                      :: active, err, gap, occ_tol, raw_gap, &
                                                            slope_cap
      REAL(KIND=dp), DIMENSION(n)                        :: evals, occ_curv, occ_filtered, occ_scale

      evals(:) = [-0.010000_dp, -0.009982_dp, 0.42_dp]
      occ_curv(:) = [0.038_dp, 0.035_dp, 0.006_dp]
      gap = evals(2) - evals(1)
      slope_cap = dense_shared_occ_slope(occ_curv(1), occ_curv(2))
      occ_scale(:) = [1.20_dp, 1.20_dp - 2.5_dp*slope_cap*gap, 0.08_dp]

      chc(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      dchc(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      DO i = 1, n
         chc(i, i) = CMPLX(evals(i), 0.0_dp, KIND=dp)
      END DO
      dchc(1, 2) = CMPLX(0.24_dp, -0.18_dp, KIND=dp)
      dchc(2, 1) = CONJG(dchc(1, 2))

      occ_filtered(:) = occ_scale(:)
      occ_tol = SQRT(EPSILON(1.0_dp))*MAX(1.0_dp, MAXVAL(ABS(occ_scale(:))))
      CALL dense_near_degenerate_occupation_plateau(evals, occ_curv, gauge_gap, occ_tol, occ_filtered)

      CALL dense_history_safe_occ_projector_derivative_chc(chc, dchc, occ_scale, projector_weight, &
                                                           docc_safe, occ_curv, gauge_gap)
      CALL dense_history_safe_occ_projector_from_chc(chc + fd_step*dchc, occ_scale, &
                                                     projector_weight, occ_mat_p, occ_curv, gauge_gap)
      CALL dense_history_safe_occ_projector_from_chc(chc - fd_step*dchc, occ_scale, &
                                                     projector_weight, occ_mat_m, occ_curv, gauge_gap)
      docc_fd(:, :) = (occ_mat_p(:, :) - occ_mat_m(:, :))/(2.0_dp*fd_step)
      CALL dense_spectral_occ_projector_derivative_chc(chc, dchc, occ_filtered, docc_ref, &
                                                       occ_curv, gauge_gap)
      docc_ref(:, :) = projector_weight*docc_ref(:, :)
      CALL dense_spectral_occ_projector_derivative_chc(chc, dchc, occ_scale, docc_raw)

      err = MAX(MAXVAL(ABS(docc_safe(:, :) - docc_fd(:, :))), &
                MAXVAL(ABS(docc_safe(:, :) - docc_ref(:, :))))
      raw_gap = MAXVAL(ABS(projector_weight*docc_raw(:, :) - docc_ref(:, :)))
      active = MAXVAL(ABS(docc_ref(:, :)))

      IF (raw_gap < active_tol .OR. active < active_tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,3(1X,ES12.5))') TRIM(label), &
            "failed with error/rawgap/active", err, raw_gap, active
      ELSE
         WRITE (*, '(A,1X,A,3(1X,ES12.5))') TRIM(label), &
            "OK with error/rawgap/active", err, raw_gap, active
      END IF

   END SUBROUTINE check_history_safe_occ_projector_slope_cap_fd

! **************************************************************************************************
!> \brief Check that the gauge-cluster finite-difference probe stays inside the cluster.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_history_safe_occ_projector_adaptive_fd(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 3
      REAL(KIND=dp), PARAMETER                           :: active_tol = 1.0E-4_dp, &
                                                            gauge_gap = 5.0E-5_dp, &
                                                            old_fd_step = 1.0E-7_dp, &
                                                            projector_weight = 0.37_dp

      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: chc, dchc, docc_adaptive, docc_fd, &
                                                            docc_old, docc_ref, occ_mat_m, &
                                                            occ_mat_p
      INTEGER                                            :: i
      REAL(KIND=dp)                                      :: adaptive_step, err, occ_tol, old_gap, &
                                                            scale, step_ratio, work_norm
      REAL(KIND=dp), DIMENSION(n)                        :: evals, occ_curv, occ_filtered, occ_scale

      evals(:) = [-0.010000_dp, -0.009982_dp, 0.34_dp]
      occ_scale(:) = [1.20_dp, 0.37_dp, 0.08_dp]
      occ_curv(:) = [0.038_dp, 0.035_dp, 0.006_dp]

      chc(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      dchc(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      DO i = 1, n
         chc(i, i) = CMPLX(evals(i), 0.0_dp, KIND=dp)
      END DO
      dchc(1, 2) = CMPLX(2.0E4_dp, -1.5E4_dp, KIND=dp)
      dchc(2, 1) = CONJG(dchc(1, 2))

      work_norm = MAXVAL(ABS(dchc(:, :)))
      scale = MAX(1.0_dp, MAXVAL(ABS(evals(:))))
      adaptive_step = MIN(old_fd_step, 1.0E-2_dp*gauge_gap/work_norm)
      adaptive_step = MAX(adaptive_step, 10.0_dp*EPSILON(1.0_dp)*scale/work_norm)

      CALL dense_history_safe_occ_projector_derivative_chc(chc, dchc, occ_scale, projector_weight, &
                                                           docc_adaptive, occ_curv, gauge_gap)
      occ_filtered(:) = occ_scale(:)
      occ_tol = SQRT(EPSILON(1.0_dp))*MAX(1.0_dp, MAXVAL(ABS(occ_filtered(:))))
      CALL dense_near_degenerate_occupation_plateau(evals, occ_curv, gauge_gap, occ_tol, occ_filtered)
      CALL dense_spectral_occ_projector_derivative_chc(chc, dchc, occ_filtered, docc_ref, &
                                                       occ_curv, gauge_gap)
      docc_ref(:, :) = projector_weight*docc_ref(:, :)
      CALL dense_history_safe_occ_projector_from_chc(chc + adaptive_step*dchc, occ_scale, &
                                                     projector_weight, occ_mat_p, occ_curv, gauge_gap)
      CALL dense_history_safe_occ_projector_from_chc(chc - adaptive_step*dchc, occ_scale, &
                                                     projector_weight, occ_mat_m, occ_curv, gauge_gap)
      docc_fd(:, :) = (occ_mat_p(:, :) - occ_mat_m(:, :))/(2.0_dp*adaptive_step)

      CALL dense_history_safe_occ_projector_from_chc(chc + old_fd_step*dchc, occ_scale, projector_weight, &
                                                     occ_mat_p, occ_curv, gauge_gap)
      CALL dense_history_safe_occ_projector_from_chc(chc - old_fd_step*dchc, occ_scale, projector_weight, &
                                                     occ_mat_m, occ_curv, gauge_gap)
      docc_old(:, :) = (occ_mat_p(:, :) - occ_mat_m(:, :))/(2.0_dp*old_fd_step)

      err = MAXVAL(ABS(docc_adaptive(:, :) - docc_ref(:, :)))
      old_gap = MAXVAL(ABS(docc_old(:, :) - docc_adaptive(:, :)))
      step_ratio = adaptive_step/old_fd_step

      IF (old_gap < active_tol .OR. step_ratio >= 0.1_dp) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,3(1X,ES12.5))') TRIM(label), &
            "failed with error/oldgap/ratio", err, old_gap, step_ratio
      ELSE
         WRITE (*, '(A,1X,A,3(1X,ES12.5))') TRIM(label), &
            "OK with error/oldgap/ratio", err, old_gap, step_ratio
      END IF

   END SUBROUTINE check_history_safe_occ_projector_adaptive_fd

! **************************************************************************************************
!> \brief Check the explicit diagonal Fermi response in the history-safe dP model.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_history_safe_occ_projector_diagonal_response(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 4
      REAL(KIND=dp), PARAMETER                           :: active_tol = 1.0E-4_dp, &
                                                            gauge_gap = 5.0E-5_dp, &
                                                            projector_weight = 0.37_dp

      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: chc, dchc, docc_safe, expected
      INTEGER                                            :: i
      REAL(KIND=dp)                                      :: active, err, occ_sum, offdiag, shift_mu
      REAL(KIND=dp), DIMENSION(n)                        :: evals, occ_curv, occ_scale, shift

      evals(:) = [-0.010000_dp, -0.009982_dp, 0.22_dp, 0.47_dp]
      occ_scale(:) = [1.18_dp, 0.41_dp, 0.23_dp, 0.08_dp]
      occ_curv(:) = [0.038_dp, 0.035_dp, 0.018_dp, 0.006_dp]
      shift(:) = [0.11_dp, -0.07_dp, 0.03_dp, -0.02_dp]

      chc(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      dchc(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      expected(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      occ_sum = SUM(occ_curv(:))
      shift_mu = DOT_PRODUCT(occ_curv(:), shift(:))/occ_sum
      DO i = 1, n
         chc(i, i) = CMPLX(evals(i), 0.0_dp, KIND=dp)
         dchc(i, i) = CMPLX(shift(i), 0.0_dp, KIND=dp)
         expected(i, i) = CMPLX(-2.0_dp*projector_weight*occ_curv(i)*(shift(i) - shift_mu), &
                                0.0_dp, KIND=dp)
      END DO

      CALL dense_history_safe_occ_projector_derivative_chc(chc, dchc, occ_scale, projector_weight, &
                                                           docc_safe, occ_curv, gauge_gap)
      err = MAXVAL(ABS(docc_safe(:, :) - expected(:, :)))
      active = MAXVAL(ABS(expected(:, :)))
      offdiag = MAXVAL(ABS(docc_safe(:, :) - &
                           CMPLX(REAL(docc_safe(:, :), KIND=dp), 0.0_dp, KIND=dp)))

      IF (err > tol .OR. active < active_tol .OR. offdiag > tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,3(1X,ES12.5))') TRIM(label), &
            "failed with error/active/imag", err, active, offdiag
      ELSE
         WRITE (*, '(A,1X,A,3(1X,ES12.5))') TRIM(label), &
            "OK with error/active/imag", err, active, offdiag
      END IF

   END SUBROUTINE check_history_safe_occ_projector_diagonal_response

! **************************************************************************************************
!> \brief Check that the gauge finite-difference dP branch keeps the fixed-N diagonal response.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_history_safe_occ_projector_fd_diagonal_response(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 3
      REAL(KIND=dp), PARAMETER                           :: active_tol = 1.0E-4_dp, &
                                                            fd_step = 1.0E-7_dp, &
                                                            gauge_gap = 5.0E-5_dp, &
                                                            projector_weight = 0.37_dp

      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: chc, dchc, docc_fd, docc_safe, expected, &
                                                            expected_old, occ_mat_m, occ_mat_p
      INTEGER                                            :: i
      REAL(KIND=dp)                                      :: charge_leak, err, missing, occ_sum, &
                                                            shift_mu
      REAL(KIND=dp), DIMENSION(n)                        :: evals, occ_curv, occ_scale, shift

      evals(:) = [-0.010000_dp, -0.009982_dp, 0.22_dp]
      occ_scale(:) = [1.18_dp, 0.41_dp, 0.23_dp]
      occ_curv(:) = [0.038_dp, 0.035_dp, 0.018_dp]
      shift(:) = [0.81_dp, -0.49_dp, 0.23_dp]

      chc(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      dchc(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      expected(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      DO i = 1, n
         chc(i, i) = CMPLX(evals(i), 0.0_dp, KIND=dp)
         dchc(i, i) = CMPLX(shift(i), 0.0_dp, KIND=dp)
      END DO
      dchc(1, 2) = CMPLX(0.24_dp, -0.18_dp, KIND=dp)
      dchc(2, 1) = CONJG(dchc(1, 2))

      CALL dense_history_safe_occ_projector_derivative_chc(chc, dchc, occ_scale, projector_weight, &
                                                           docc_safe, occ_curv, gauge_gap)
      CALL dense_history_safe_occ_projector_from_chc(chc + fd_step*dchc, occ_scale, &
                                                     projector_weight, occ_mat_p, occ_curv, gauge_gap)
      CALL dense_history_safe_occ_projector_from_chc(chc - fd_step*dchc, occ_scale, &
                                                     projector_weight, occ_mat_m, occ_curv, gauge_gap)
      docc_fd(:, :) = (occ_mat_p(:, :) - occ_mat_m(:, :))/(2.0_dp*fd_step)
      expected(:, :) = docc_fd(:, :)
      expected_old(:, :) = docc_fd(:, :)

      occ_sum = SUM(occ_curv(:))
      shift_mu = DOT_PRODUCT(occ_curv(:), shift(:))/occ_sum
      DO i = 1, n
         expected(i, i) = CMPLX(-2.0_dp*projector_weight*occ_curv(i)*(shift(i) - shift_mu), &
                                0.0_dp, KIND=dp)
      END DO

      err = MAXVAL(ABS(docc_safe(:, :) - expected(:, :)))
      missing = MAXVAL(ABS(expected(:, :) - expected_old(:, :)))
      charge_leak = 0.0_dp
      DO i = 1, n
         charge_leak = charge_leak + REAL(expected(i, i), KIND=dp)
      END DO
      charge_leak = ABS(charge_leak)

      IF (missing < active_tol .OR. charge_leak > tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,3(1X,ES12.5))') TRIM(label), &
            "failed with error/missing/charge", err, missing, charge_leak
      ELSE
         WRITE (*, '(A,1X,A,3(1X,ES12.5))') TRIM(label), &
            "OK with error/missing/charge", err, missing, charge_leak
      END IF

   END SUBROUTINE check_history_safe_occ_projector_fd_diagonal_response

! **************************************************************************************************
!> \brief Check the mixed gauge/physical finite-difference projector dP branch.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_history_safe_occ_projector_fd_mixed_response(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 4
      REAL(KIND=dp), PARAMETER                           :: active_tol = 1.0E-4_dp, &
                                                            fd_step = 1.0E-7_dp, &
                                                            gauge_gap = 5.0E-5_dp, &
                                                            projector_weight = 0.37_dp

      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: chc, dchc, docc_fd, docc_raw, docc_safe, &
                                                            expected, expected_old, occ_mat_m, &
                                                            occ_mat_p
      INTEGER                                            :: i
      REAL(KIND=dp)                                      :: charge_leak, diagonal_fix, err, &
                                                            gauge_pair, occ_sum, physical_pair, &
                                                            raw_gauge, shift_mu
      REAL(KIND=dp), DIMENSION(n)                        :: evals, occ_curv, occ_scale, shift

      evals(:) = [-0.010000_dp, -0.009982_dp, 0.18_dp, 0.37_dp]
      occ_scale(:) = [1.20_dp, 0.40_dp, 0.23_dp, 0.06_dp]
      occ_curv(:) = [0.038_dp, 0.035_dp, 0.018_dp, 0.006_dp]
      shift(:) = [0.81_dp, -0.49_dp, 0.23_dp, -0.08_dp]

      chc(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      dchc(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      DO i = 1, n
         chc(i, i) = CMPLX(evals(i), 0.0_dp, KIND=dp)
         dchc(i, i) = CMPLX(shift(i), 0.0_dp, KIND=dp)
      END DO
      dchc(1, 2) = CMPLX(0.24_dp, -0.18_dp, KIND=dp)
      dchc(2, 1) = CONJG(dchc(1, 2))
      dchc(3, 4) = CMPLX(-0.19_dp, 0.11_dp, KIND=dp)
      dchc(4, 3) = CONJG(dchc(3, 4))

      CALL dense_history_safe_occ_projector_derivative_chc(chc, dchc, occ_scale, projector_weight, &
                                                           docc_safe, occ_curv, gauge_gap)
      CALL dense_history_safe_occ_projector_from_chc(chc + fd_step*dchc, occ_scale, projector_weight, &
                                                     occ_mat_p, occ_curv, gauge_gap)
      CALL dense_history_safe_occ_projector_from_chc(chc - fd_step*dchc, occ_scale, projector_weight, &
                                                     occ_mat_m, occ_curv, gauge_gap)
      docc_fd(:, :) = (occ_mat_p(:, :) - occ_mat_m(:, :))/(2.0_dp*fd_step)
      expected(:, :) = docc_fd(:, :)
      expected_old(:, :) = docc_fd(:, :)

      occ_sum = SUM(occ_curv(:))
      shift_mu = DOT_PRODUCT(occ_curv(:), shift(:))/occ_sum
      DO i = 1, n
         expected(i, i) = CMPLX(-2.0_dp*projector_weight*occ_curv(i)*(shift(i) - shift_mu), &
                                0.0_dp, KIND=dp)
      END DO

      CALL dense_spectral_occ_projector_derivative_chc(chc, dchc, occ_scale, docc_raw)
      err = MAXVAL(ABS(docc_safe(:, :) - expected(:, :)))
      diagonal_fix = MAXVAL(ABS(expected_old(:, :) - expected(:, :)))
      raw_gauge = MAX(ABS(projector_weight*docc_raw(1, 2)), ABS(projector_weight*docc_raw(2, 1)))
      gauge_pair = MAX(ABS(docc_safe(1, 2)), ABS(docc_safe(2, 1)))
      physical_pair = MAX(ABS(docc_safe(3, 4)), ABS(docc_safe(4, 3)))
      charge_leak = 0.0_dp
      DO i = 1, n
         charge_leak = charge_leak + REAL(expected(i, i), KIND=dp)
      END DO
      charge_leak = ABS(charge_leak)

      IF (diagonal_fix < active_tol .OR. raw_gauge < active_tol .OR. &
          gauge_pair >= 0.01_dp*raw_gauge .OR. physical_pair < active_tol .OR. charge_leak > tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,6(1X,ES12.5))') TRIM(label), &
            "failed with error/diagfix/rawgauge/gauge/physical/charge", err, diagonal_fix, &
            raw_gauge, gauge_pair, physical_pair, charge_leak
      ELSE
         WRITE (*, '(A,1X,A,6(1X,ES12.5))') TRIM(label), &
            "OK with error/diagfix/rawgauge/gauge/physical/charge", err, diagonal_fix, &
            raw_gauge, gauge_pair, physical_pair, charge_leak
      END IF

   END SUBROUTINE check_history_safe_occ_projector_fd_mixed_response

! **************************************************************************************************
!> \brief Check that the finite-difference projector dP branch uses the fixed-N curvature metric.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_history_safe_occ_projector_fd_fixed_n_trace(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 5
      REAL(KIND=dp), PARAMETER                           :: active_tol = 1.0E-4_dp, &
                                                            gauge_gap = 5.0E-5_dp, &
                                                            projector_weight = 0.37_dp

      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: chc, dchc, dchc_diag, docc_ref, &
                                                            docc_safe, uniform_expected
      INTEGER                                            :: i
      REAL(KIND=dp)                                      :: charge_leak, diag_err, occ_sum, &
                                                            shift_mu, uniform_gap, uniform_shift
      REAL(KIND=dp), DIMENSION(n)                        :: evals, occ_curv, occ_scale, shift

      evals(:) = [-0.012000_dp, -0.011982_dp, 0.045000_dp, 0.045030_dp, 0.29_dp]
      occ_scale(:) = [1.22_dp, 0.38_dp, 0.31_dp, 0.27_dp, 0.04_dp]
      occ_curv(:) = [0.041_dp, 0.019_dp, 0.013_dp, 0.007_dp, 0.002_dp]
      shift(:) = [0.61_dp, -0.27_dp, 0.42_dp, -0.19_dp, 0.08_dp]

      chc(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      dchc(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      dchc_diag(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      DO i = 1, n
         chc(i, i) = CMPLX(evals(i), 0.0_dp, KIND=dp)
         dchc(i, i) = CMPLX(shift(i), 0.0_dp, KIND=dp)
         dchc_diag(i, i) = dchc(i, i)
      END DO
      dchc(1, 2) = CMPLX(0.33_dp, -0.25_dp, KIND=dp)
      dchc(2, 1) = CONJG(dchc(1, 2))
      dchc(3, 4) = CMPLX(-0.29_dp, 0.17_dp, KIND=dp)
      dchc(4, 3) = CONJG(dchc(3, 4))

      CALL dense_history_safe_occ_projector_derivative_chc(chc, dchc_diag, occ_scale, projector_weight, &
                                                           docc_ref, occ_curv, gauge_gap)
      CALL dense_history_safe_occ_projector_derivative_chc(chc, dchc, occ_scale, projector_weight, &
                                                           docc_safe, occ_curv, gauge_gap)
      uniform_expected(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      occ_sum = SUM(occ_curv(:))
      shift_mu = DOT_PRODUCT(occ_curv(:), shift(:))/occ_sum
      uniform_shift = -2.0_dp*projector_weight*DOT_PRODUCT(occ_curv(:), shift(:))/REAL(n, KIND=dp)
      DO i = 1, n
         uniform_expected(i, i) = CMPLX(-2.0_dp*projector_weight*occ_curv(i)*shift(i) - &
                                        uniform_shift, 0.0_dp, KIND=dp)
      END DO

      diag_err = 0.0_dp
      charge_leak = 0.0_dp
      DO i = 1, n
         diag_err = MAX(diag_err, ABS(docc_safe(i, i) - docc_ref(i, i)))
         charge_leak = charge_leak + REAL(docc_safe(i, i), KIND=dp)
      END DO
      uniform_gap = MAXVAL(ABS(uniform_expected(:, :) - docc_ref(:, :)))
      charge_leak = ABS(charge_leak)

      IF (uniform_gap < active_tol .OR. charge_leak > tol .OR. &
          ABS(shift_mu) < active_tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "failed with diagerr/mugap/uniformgap/charge", diag_err, ABS(shift_mu), uniform_gap, charge_leak
      ELSE
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "OK with diagerr/mugap/uniformgap/charge", diag_err, ABS(shift_mu), uniform_gap, charge_leak
      END IF

   END SUBROUTINE check_history_safe_occ_projector_fd_fixed_n_trace

! **************************************************************************************************
!> \brief Check that active Fermi-response pairs can cluster beyond the fixed gauge gap.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_history_safe_occ_projector_thermal_cluster_gap(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 3
      REAL(KIND=dp), PARAMETER :: active_tol = 1.0E-4_dp, fd_step = 1.0E-7_dp, &
         gauge_gap = 5.0E-5_dp, projector_weight = 0.37_dp, projector_window_gap = 8.0E-3_dp

      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: chc, dchc, docc_expected, docc_fd, &
                                                            docc_raw, docc_safe, occ_mat_m, &
                                                            occ_mat_p
      INTEGER                                            :: i
      LOGICAL                                            :: fd_used
      REAL(KIND=dp)                                      :: err, fd_gap, gap, occ_tol, old_pair_gap, &
                                                            pair_gap, raw_noise, safe_active
      REAL(KIND=dp), DIMENSION(n)                        :: evals, occ_curv, occ_scale, occ_slope

      evals(:) = [-0.010000_dp, -0.004000_dp, 0.31_dp]
      occ_scale(:) = [1.20_dp, 0.32_dp, 0.06_dp]
      occ_curv(:) = [650.0_dp, 610.0_dp, 0.004_dp]

      chc(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      dchc(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      DO i = 1, n
         chc(i, i) = CMPLX(evals(i), 0.0_dp, KIND=dp)
      END DO
      dchc(1, 2) = CMPLX(0.19_dp, -0.12_dp, KIND=dp)
      dchc(2, 1) = CONJG(dchc(1, 2))

      CALL dense_history_safe_occ_projector_derivative_chc(chc, dchc, occ_scale, projector_weight, &
                                                           docc_safe, occ_curv, gauge_gap, fd_used, &
                                                           projector_window_gap)
      CALL dense_history_safe_occ_projector_from_chc(chc + fd_step*dchc, occ_scale, projector_weight, &
                                                     occ_mat_p, occ_curv, gauge_gap)
      CALL dense_history_safe_occ_projector_from_chc(chc - fd_step*dchc, occ_scale, projector_weight, &
                                                     occ_mat_m, occ_curv, gauge_gap)
      docc_fd(:, :) = (occ_mat_p(:, :) - occ_mat_m(:, :))/(2.0_dp*fd_step)
      CALL dense_spectral_occ_projector_derivative_chc(chc, dchc, occ_scale, docc_raw)

      occ_tol = SQRT(EPSILON(1.0_dp))*MAX(1.0_dp, MAXVAL(ABS(occ_scale(:))))
      CALL dense_monotone_occupation_slopes(occ_scale, occ_tol, occ_slope)
      CALL dense_near_degenerate_occupation_plateau(evals, occ_curv, gauge_gap, occ_tol, occ_slope, &
                                                    projector_window_gap)
      docc_expected(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      CALL dense_fixed_n_projector_cluster_kernel(evals, dchc, occ_curv, occ_slope, gauge_gap, &
                                                  projector_weight, docc_expected, projector_window_gap)

      gap = evals(2) - evals(1)
      old_pair_gap = dense_projector_occupation_pair_gap(gauge_gap, occ_curv(1), occ_curv(2), &
                                                         MAXVAL(occ_curv(:)))
      pair_gap = dense_projector_occupation_pair_gap(gauge_gap, occ_curv(1), occ_curv(2), &
                                                     MAXVAL(occ_curv(:)), projector_window_gap)
      err = MAXVAL(ABS(docc_safe(:, :) - docc_expected(:, :)))
      fd_gap = MAXVAL(ABS(docc_fd(:, :) - docc_expected(:, :)))
      raw_noise = MAXVAL(ABS(projector_weight*docc_raw(:, :) - docc_safe(:, :)))
      safe_active = MAXVAL(ABS(docc_safe(:, :)))

      IF (.NOT. fd_used .OR. gap <= old_pair_gap .OR. gap > pair_gap .OR. err > tol .OR. &
          fd_gap < active_tol .OR. raw_noise < active_tol .OR. safe_active < active_tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,7(1X,ES12.5),1(1X,L1))') TRIM(label), &
            "failed with err/fdgap/gap/oldpgap/pgap/rawnoise/active/fd", err, fd_gap, gap, &
            old_pair_gap, pair_gap, raw_noise, safe_active, fd_used
      ELSE
         WRITE (*, '(A,1X,A,7(1X,ES12.5),1(1X,L1))') TRIM(label), &
            "OK with err/fdgap/gap/oldpgap/pgap/rawnoise/active/fd", err, fd_gap, gap, &
            old_pair_gap, pair_gap, raw_noise, safe_active, fd_used
      END IF

   END SUBROUTINE check_history_safe_occ_projector_thermal_cluster_gap

! **************************************************************************************************
!> \brief Check that the finite projector value itself uses the thermal cluster window.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_history_safe_occ_projector_thermal_value_window(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 3
      REAL(KIND=dp), PARAMETER :: active_tol = 1.0E-4_dp, &
         gauge_gap = 5.0E-5_dp, projector_weight = 0.37_dp, projector_window_gap = 8.0E-3_dp

      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: chc, occ_expected, occ_raw, occ_safe
      INTEGER                                            :: i
      REAL(KIND=dp)                                      :: err, gap, occ_tol, pair_gap, raw_gap
      REAL(KIND=dp), DIMENSION(n)                        :: evals, occ_curv, occ_scale, occ_window

      evals(:) = [-0.010000_dp, -0.004000_dp, 0.31_dp]
      occ_scale(:) = [1.20_dp, 0.20_dp, 0.06_dp]
      occ_curv(:) = [10.0_dp, 9.0_dp, 0.004_dp]

      chc(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      occ_expected(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      DO i = 1, n
         chc(i, i) = CMPLX(evals(i), 0.0_dp, KIND=dp)
      END DO

      CALL dense_history_safe_occ_projector_from_chc(chc, occ_scale, projector_weight, occ_safe, &
                                                     occ_curv, gauge_gap, projector_window_gap)
      CALL dense_history_safe_occ_projector_from_chc(chc, occ_scale, projector_weight, occ_raw, &
                                                     occ_curv, gauge_gap)

      occ_tol = SQRT(EPSILON(1.0_dp))*MAX(1.0_dp, MAXVAL(ABS(occ_scale(:))))
      CALL dense_monotone_occupation_slopes(occ_scale, occ_tol, occ_window)
      CALL dense_near_degenerate_occupation_plateau(evals, occ_curv, gauge_gap, occ_tol, occ_window, &
                                                    projector_window_gap)
      DO i = 1, n
         occ_expected(i, i) = CMPLX(projector_weight*occ_window(i) + &
                                    (1.0_dp - projector_weight)*occ_scale(i), 0.0_dp, KIND=dp)
      END DO

      err = MAXVAL(ABS(occ_safe(:, :) - occ_expected(:, :)))
      raw_gap = MAXVAL(ABS(occ_safe(:, :) - occ_raw(:, :)))
      gap = evals(2) - evals(1)
      pair_gap = dense_projector_occupation_pair_gap(gauge_gap, occ_curv(1), occ_curv(2), &
                                                     MAXVAL(occ_curv(:)), projector_window_gap)

      IF (gap > pair_gap .OR. err > tol .OR. raw_gap < active_tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "failed with err/rawgap/gap/pairgap", err, raw_gap, gap, pair_gap
      ELSE
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "OK with err/rawgap/gap/pairgap", err, raw_gap, gap, pair_gap
      END IF

   END SUBROUTINE check_history_safe_occ_projector_thermal_value_window

! **************************************************************************************************
!> \brief Check that the no-FD projector derivative still uses the thermal cluster window.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_history_safe_occ_projector_thermal_no_fd_offblock(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 3
      REAL(KIND=dp), PARAMETER :: active_tol = 1.0E-4_dp, &
         gauge_gap = 5.0E-5_dp, projector_weight = 0.37_dp, projector_window_gap = 8.0E-3_dp

      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: chc, dchc, docc_expected, &
                                                            docc_raw, docc_safe
      INTEGER                                            :: i
      LOGICAL                                            :: fd_used
      REAL(KIND=dp)                                      :: active, err, gap, occ_tol, pair_gap, &
                                                            raw_gap
      REAL(KIND=dp), DIMENSION(n)                        :: evals, occ_curv, occ_raw, &
                                                            occ_scale, occ_window

      evals(:) = [-0.010000_dp, -0.004000_dp, 0.31_dp]
      occ_scale(:) = [1.20_dp, 0.20_dp, 0.06_dp]
      occ_curv(:) = [10.0_dp, 9.0_dp, 0.004_dp]

      chc(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      dchc(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      DO i = 1, n
         chc(i, i) = CMPLX(evals(i), 0.0_dp, KIND=dp)
      END DO
      dchc(1, 3) = CMPLX(0.13_dp, -0.08_dp, KIND=dp)
      dchc(2, 3) = CMPLX(-0.07_dp, 0.05_dp, KIND=dp)
      dchc(3, 1) = CONJG(dchc(1, 3))
      dchc(3, 2) = CONJG(dchc(2, 3))

      CALL dense_history_safe_occ_projector_derivative_chc(chc, dchc, occ_scale, projector_weight, &
                                                           docc_safe, occ_curv, gauge_gap, fd_used, &
                                                           projector_window_gap)

      occ_tol = SQRT(EPSILON(1.0_dp))*MAX(1.0_dp, MAXVAL(ABS(occ_scale(:))))
      CALL dense_monotone_occupation_slopes(occ_scale, occ_tol, occ_raw)
      occ_window(:) = occ_raw(:)
      CALL dense_near_degenerate_occupation_plateau(evals, occ_curv, gauge_gap, occ_tol, occ_window, &
                                                    projector_window_gap)
      CALL dense_spectral_occ_projector_derivative_chc(chc, dchc, occ_window, docc_expected, &
                                                       occ_curv, gauge_gap, projector_window_gap)
      CALL dense_spectral_occ_projector_derivative_chc(chc, dchc, occ_raw, docc_raw, occ_curv, gauge_gap)
      docc_expected(:, :) = projector_weight*docc_expected(:, :)
      docc_raw(:, :) = projector_weight*docc_raw(:, :)

      err = MAXVAL(ABS(docc_safe(:, :) - docc_expected(:, :)))
      raw_gap = MAXVAL(ABS(docc_expected(:, :) - docc_raw(:, :)))
      active = MAXVAL(ABS(docc_expected(:, :)))
      gap = evals(2) - evals(1)
      pair_gap = dense_projector_occupation_pair_gap(gauge_gap, occ_curv(1), occ_curv(2), &
                                                     MAXVAL(occ_curv(:)), projector_window_gap)

      IF (fd_used .OR. gap > pair_gap .OR. err > tol .OR. raw_gap < active_tol .OR. &
          active < active_tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,5(1X,ES12.5),1(1X,L1))') TRIM(label), &
            "failed with err/rawgap/active/gap/pairgap/fd", err, raw_gap, active, gap, pair_gap, fd_used
      ELSE
         WRITE (*, '(A,1X,A,5(1X,ES12.5),1(1X,L1))') TRIM(label), &
            "OK with err/rawgap/active/gap/pairgap/fd", err, raw_gap, active, gap, pair_gap, fd_used
      END IF

   END SUBROUTINE check_history_safe_occ_projector_thermal_no_fd_offblock

! **************************************************************************************************
!> \brief Check that diagonal splitting inside a near-degenerate cluster uses the cluster path.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_history_safe_occ_projector_diagonal_cluster_fd(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 4
      REAL(KIND=dp), PARAMETER                           :: active_tol = 1.0E-4_dp, &
                                                            gauge_gap = 5.0E-5_dp, &
                                                            projector_weight = 0.37_dp

      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: chc, dchc, docc_expected, docc_safe
      INTEGER                                            :: i, j
      LOGICAL                                            :: fd_used
      REAL(KIND=dp)                                      :: charge_leak, expected_err, offdiag_norm, &
                                                            split_gap
      REAL(KIND=dp), DIMENSION(n)                        :: evals, occ_curv, occ_scale, occ_slope

      evals(:) = [-0.010000_dp, -0.009982_dp, 0.23_dp, 0.47_dp]
      occ_scale(:) = [1.18_dp, 0.39_dp, 0.20_dp, 0.07_dp]
      occ_curv(:) = [0.038_dp, 0.018_dp, 0.018_dp, 0.006_dp]

      chc(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      dchc(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      DO i = 1, n
         chc(i, i) = CMPLX(evals(i), 0.0_dp, KIND=dp)
      END DO
      dchc(1, 1) = CMPLX(0.11_dp, 0.0_dp, KIND=dp)
      dchc(2, 2) = CMPLX(-0.05_dp, 0.0_dp, KIND=dp)
      dchc(3, 3) = CMPLX(0.02_dp, 0.0_dp, KIND=dp)
      dchc(4, 4) = CMPLX(-0.01_dp, 0.0_dp, KIND=dp)

      CALL dense_history_safe_occ_projector_derivative_chc(chc, dchc, occ_scale, projector_weight, &
                                                           docc_safe, occ_curv, gauge_gap, fd_used)

      occ_slope(:) = occ_scale(:)
      CALL dense_monotone_occupation_slopes(occ_scale, &
                                            SQRT(EPSILON(1.0_dp))*MAX(1.0_dp, MAXVAL(ABS(occ_scale(:)))), &
                                            occ_slope)
      CALL dense_near_degenerate_occupation_plateau(evals, occ_curv, gauge_gap, &
                                                    SQRT(EPSILON(1.0_dp))* &
                                                    MAX(1.0_dp, MAXVAL(ABS(occ_scale(:)))), occ_slope)
      docc_expected(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      CALL dense_fixed_n_projector_cluster_kernel(evals, dchc, occ_curv, occ_slope, gauge_gap, &
                                                  projector_weight, docc_expected)

      expected_err = MAXVAL(ABS(docc_safe(:, :) - docc_expected(:, :)))
      charge_leak = 0.0_dp
      offdiag_norm = 0.0_dp
      DO i = 1, n
         charge_leak = charge_leak + REAL(docc_safe(i, i), KIND=dp)
         DO j = 1, n
            IF (i /= j) offdiag_norm = MAX(offdiag_norm, ABS(dchc(i, j)))
         END DO
      END DO
      charge_leak = ABS(charge_leak)
      split_gap = ABS(REAL(dchc(1, 1), KIND=dp) - REAL(dchc(2, 2), KIND=dp))

      IF (.NOT. fd_used .OR. expected_err > tol .OR. charge_leak > tol .OR. &
          offdiag_norm > tol .OR. split_gap < active_tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,5(1X,ES12.5),1(1X,L1))') TRIM(label), &
            "failed with err/charge/offdiag/split/tol/fd", expected_err, charge_leak, &
            offdiag_norm, split_gap, tol, fd_used
      ELSE
         WRITE (*, '(A,1X,A,5(1X,ES12.5),1(1X,L1))') TRIM(label), &
            "OK with err/charge/offdiag/split/tol/fd", expected_err, charge_leak, &
            offdiag_norm, split_gap, tol, fd_used
      END IF

   END SUBROUTINE check_history_safe_occ_projector_diagonal_cluster_fd

! **************************************************************************************************
!> \brief Check that a near-degenerate physical cluster uses degenerate perturbation theory.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_history_safe_occ_projector_degenerate_physical_split(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 4
      REAL(KIND=dp), PARAMETER                           :: active_tol = 1.0E-4_dp, &
                                                            gauge_gap = 5.0E-5_dp, &
                                                            projector_weight = 0.37_dp

      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: chc, cluster_kernel, cluster_response, &
                                                            cluster_work, dchc, docc_safe, evec
      INTEGER                                            :: i
      REAL(KIND=dp)                                      :: block_err, charge_leak, legacy_gap, &
                                                            mu_den, mu_num, mu_shift, outside_err
      REAL(KIND=dp), DIMENSION(n)                        :: evals, modal_curv, occ_curv, occ_scale, &
                                                            split_eval

      evals(:) = [-0.010000_dp, -0.009982_dp, 0.23_dp, 0.47_dp]
      occ_scale(:) = [1.18_dp, 0.39_dp, 0.20_dp, 0.07_dp]
      occ_curv(:) = [0.038_dp, 0.018_dp, 0.018_dp, 0.006_dp]

      chc(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      dchc(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      DO i = 1, n
         chc(i, i) = CMPLX(evals(i), 0.0_dp, KIND=dp)
      END DO
      dchc(1, 1) = CMPLX(0.11_dp, 0.0_dp, KIND=dp)
      dchc(2, 2) = CMPLX(-0.05_dp, 0.0_dp, KIND=dp)
      dchc(3, 3) = CMPLX(0.02_dp, 0.0_dp, KIND=dp)
      dchc(4, 4) = CMPLX(-0.01_dp, 0.0_dp, KIND=dp)
      dchc(1, 2) = CMPLX(0.24_dp, -0.18_dp, KIND=dp)
      dchc(2, 1) = CONJG(dchc(1, 2))

      CALL dense_history_safe_occ_projector_derivative_chc(chc, dchc, occ_scale, projector_weight, &
                                                           docc_safe, occ_curv, gauge_gap)

      cluster_work(1:2, 1:2) = dchc(1:2, 1:2)
      CALL diag_complex(cluster_work(1:2, 1:2), evec(1:2, 1:2), split_eval(1:2))
      CALL dense_fixed_n_projector_modal_curvature(evec(1:2, 1:2), occ_curv(1:2), modal_curv(1:2))
      mu_num = modal_curv(1)*split_eval(1) + modal_curv(2)*split_eval(2) + &
               occ_curv(3)*REAL(dchc(3, 3), KIND=dp) + occ_curv(4)*REAL(dchc(4, 4), KIND=dp)
      mu_den = modal_curv(1) + modal_curv(2) + occ_curv(3) + occ_curv(4)
      mu_shift = mu_num/mu_den
      cluster_kernel(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      DO i = 1, 2
         cluster_kernel(i, i) = CMPLX(-2.0_dp*projector_weight*modal_curv(i)* &
                                      (split_eval(i) - mu_shift), 0.0_dp, KIND=dp)
      END DO
      cluster_response(1:2, 1:2) = MATMUL(evec(1:2, 1:2), &
                                          MATMUL(cluster_kernel(1:2, 1:2), &
                                                 CONJG(TRANSPOSE(evec(1:2, 1:2)))))

      block_err = MAXVAL(ABS(docc_safe(1:2, 1:2) - cluster_response(1:2, 1:2)))
      outside_err = MAX(ABS(docc_safe(3, 3) - &
                            CMPLX(-2.0_dp*projector_weight*occ_curv(3)* &
                                  (REAL(dchc(3, 3), KIND=dp) - mu_shift), 0.0_dp, KIND=dp)), &
                        ABS(docc_safe(4, 4) - &
                            CMPLX(-2.0_dp*projector_weight*occ_curv(4)* &
                                  (REAL(dchc(4, 4), KIND=dp) - mu_shift), 0.0_dp, KIND=dp)))
      charge_leak = ABS(REAL(docc_safe(1, 1) + docc_safe(2, 2) + docc_safe(3, 3) + &
                             docc_safe(4, 4), KIND=dp))
      legacy_gap = ABS(docc_safe(1, 1) - CMPLX(-2.0_dp*projector_weight*occ_curv(1)* &
                                               (REAL(dchc(1, 1), KIND=dp) - mu_shift), &
                                               0.0_dp, KIND=dp))

      IF (block_err > tol .OR. outside_err > tol .OR. charge_leak > tol .OR. legacy_gap < active_tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "failed with block/outside/charge/legacy", block_err, outside_err, charge_leak, legacy_gap
      ELSE
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "OK with block/outside/charge/legacy", block_err, outside_err, charge_leak, legacy_gap
      END IF

   END SUBROUTINE check_history_safe_occ_projector_degenerate_physical_split

! **************************************************************************************************
!> \brief Check that mixed active/inert clusters keep modal physical curvature.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_history_safe_occ_projector_modal_curvature_split(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 4
      REAL(KIND=dp), PARAMETER                           :: active_tol = 1.0E-4_dp, &
                                                            gauge_gap = 5.0E-5_dp, &
                                                            projector_weight = 0.37_dp

      COMPLEX(KIND=dp), DIMENSION(n, n) :: chc, cluster_kernel, cluster_min_kernel, &
         cluster_min_response, cluster_response, cluster_work, dchc, docc_safe, evec
      INTEGER                                            :: a, i
      REAL(KIND=dp) :: block_err, charge_leak, curv_floor, expected_norm, min_norm, modal_gain, &
         mu_den, mu_min_den, mu_min_num, mu_num, mu_shift, mu_shift_min
      REAL(KIND=dp), DIMENSION(n)                        :: evals, modal_curv, occ_curv, occ_scale, &
                                                            split_eval

      evals(:) = [-0.010000_dp, -0.009982_dp, 0.24_dp, 0.49_dp]
      occ_scale(:) = [1.15_dp, 0.24_dp, 0.19_dp, 0.05_dp]
      occ_curv(:) = [0.050_dp, 0.000_dp, 0.018_dp, 0.006_dp]

      chc(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      dchc(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      DO i = 1, n
         chc(i, i) = CMPLX(evals(i), 0.0_dp, KIND=dp)
      END DO
      dchc(1, 1) = CMPLX(0.04_dp, 0.0_dp, KIND=dp)
      dchc(2, 2) = CMPLX(-0.02_dp, 0.0_dp, KIND=dp)
      dchc(3, 3) = CMPLX(0.015_dp, 0.0_dp, KIND=dp)
      dchc(4, 4) = CMPLX(-0.01_dp, 0.0_dp, KIND=dp)
      dchc(1, 2) = CMPLX(0.27_dp, -0.21_dp, KIND=dp)
      dchc(2, 1) = CONJG(dchc(1, 2))

      CALL dense_history_safe_occ_projector_derivative_chc(chc, dchc, occ_scale, projector_weight, &
                                                           docc_safe, occ_curv, gauge_gap)

      cluster_work(1:2, 1:2) = dchc(1:2, 1:2)
      CALL diag_complex(cluster_work(1:2, 1:2), evec(1:2, 1:2), split_eval(1:2))
      modal_curv(1:2) = 0.0_dp
      CALL dense_fixed_n_projector_modal_curvature(evec(1:2, 1:2), occ_curv(1:2), modal_curv(1:2))

      mu_num = modal_curv(1)*split_eval(1) + modal_curv(2)*split_eval(2) + &
               occ_curv(3)*REAL(dchc(3, 3), KIND=dp) + occ_curv(4)*REAL(dchc(4, 4), KIND=dp)
      mu_den = modal_curv(1) + modal_curv(2) + occ_curv(3) + occ_curv(4)
      mu_shift = mu_num/mu_den
      cluster_kernel(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      DO i = 1, 2
         cluster_kernel(i, i) = CMPLX(-2.0_dp*projector_weight*modal_curv(i)* &
                                      (split_eval(i) - mu_shift), 0.0_dp, KIND=dp)
      END DO
      cluster_response(1:2, 1:2) = MATMUL(evec(1:2, 1:2), &
                                          MATMUL(cluster_kernel(1:2, 1:2), &
                                                 CONJG(TRANSPOSE(evec(1:2, 1:2)))))

      curv_floor = MIN(occ_curv(1), occ_curv(2))
      mu_min_num = curv_floor*(split_eval(1) + split_eval(2)) + &
                   occ_curv(3)*REAL(dchc(3, 3), KIND=dp) + occ_curv(4)*REAL(dchc(4, 4), KIND=dp)
      mu_min_den = 2.0_dp*curv_floor + occ_curv(3) + occ_curv(4)
      mu_shift_min = mu_min_num/mu_min_den
      cluster_min_kernel(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      DO i = 1, 2
         cluster_min_kernel(i, i) = CMPLX(-2.0_dp*projector_weight*curv_floor* &
                                          (split_eval(i) - mu_shift_min), 0.0_dp, KIND=dp)
      END DO
      cluster_min_response(1:2, 1:2) = MATMUL(evec(1:2, 1:2), &
                                              MATMUL(cluster_min_kernel(1:2, 1:2), &
                                                     CONJG(TRANSPOSE(evec(1:2, 1:2)))))

      block_err = MAXVAL(ABS(docc_safe(1:2, 1:2) - cluster_response(1:2, 1:2)))
      charge_leak = ABS(REAL(docc_safe(1, 1) + docc_safe(2, 2) + docc_safe(3, 3) + &
                             docc_safe(4, 4), KIND=dp))
      expected_norm = MAXVAL(ABS(cluster_response(1:2, 1:2)))
      min_norm = MAXVAL(ABS(cluster_min_response(1:2, 1:2)))
      modal_gain = expected_norm - min_norm

      IF (block_err > tol .OR. charge_leak > tol .OR. modal_gain < active_tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,5(1X,ES12.5))') TRIM(label), &
            "failed with block/charge/gain/expected/min", block_err, charge_leak, modal_gain, &
            expected_norm, min_norm
      ELSE
         WRITE (*, '(A,1X,A,5(1X,ES12.5))') TRIM(label), &
            "OK with block/charge/gain/expected/min", block_err, charge_leak, modal_gain, &
            expected_norm, min_norm
      END IF

   END SUBROUTINE check_history_safe_occ_projector_modal_curvature_split

! **************************************************************************************************
!> \brief Check that damped modal response curvature also closes the fixed-N mu projection.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_history_safe_occ_projector_modal_mu_closure(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 4
      REAL(KIND=dp), PARAMETER                           :: active_tol = 1.0E-5_dp, &
                                                            gauge_gap = 5.0E-5_dp, &
                                                            projector_weight = 0.37_dp

      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: chc, cluster_work, dchc, docc_safe, evec
      INTEGER                                            :: a, i
      REAL(KIND=dp)                                      :: charge_leak, mu_raw, mu_response, &
                                                            raw_charge_leak, response_charge_leak, &
                                                            response_gap
      REAL(KIND=dp), DIMENSION(n)                        :: evals, modal_curv, modal_raw, occ_curv, &
                                                            occ_scale, split_eval

      evals(:) = [-0.010000_dp, -0.009982_dp, 0.24_dp, 0.49_dp]
      occ_scale(:) = [1.15_dp, 0.24_dp, 0.19_dp, 0.05_dp]
      occ_curv(:) = [0.038_dp, 1.0E-5_dp, 0.018_dp, 0.006_dp]

      chc(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      dchc(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      DO i = 1, n
         chc(i, i) = CMPLX(evals(i), 0.0_dp, KIND=dp)
      END DO
      dchc(1, 1) = CMPLX(0.04_dp, 0.0_dp, KIND=dp)
      dchc(2, 2) = CMPLX(-0.02_dp, 0.0_dp, KIND=dp)
      dchc(3, 3) = CMPLX(0.015_dp, 0.0_dp, KIND=dp)
      dchc(4, 4) = CMPLX(-0.01_dp, 0.0_dp, KIND=dp)
      dchc(1, 2) = CMPLX(0.27_dp, -0.21_dp, KIND=dp)
      dchc(2, 1) = CONJG(dchc(1, 2))

      CALL dense_history_safe_occ_projector_derivative_chc(chc, dchc, occ_scale, projector_weight, &
                                                           docc_safe, occ_curv, gauge_gap)

      cluster_work(1:2, 1:2) = dchc(1:2, 1:2)
      CALL diag_complex(cluster_work(1:2, 1:2), evec(1:2, 1:2), split_eval(1:2))
      CALL dense_fixed_n_projector_modal_curvature(evec(1:2, 1:2), occ_curv(1:2), modal_curv(1:2))

      modal_raw(:) = 0.0_dp
      DO a = 1, 2
         DO i = 1, 2
            modal_raw(a) = modal_raw(a) + ABS(evec(i, a))**2*MAX(0.0_dp, occ_curv(i))
         END DO
      END DO

      mu_response = (modal_curv(1)*split_eval(1) + modal_curv(2)*split_eval(2) + &
                     occ_curv(3)*REAL(dchc(3, 3), KIND=dp) + &
                     occ_curv(4)*REAL(dchc(4, 4), KIND=dp))/ &
                    (modal_curv(1) + modal_curv(2) + occ_curv(3) + occ_curv(4))
      mu_raw = (modal_raw(1)*split_eval(1) + modal_raw(2)*split_eval(2) + &
                occ_curv(3)*REAL(dchc(3, 3), KIND=dp) + &
                occ_curv(4)*REAL(dchc(4, 4), KIND=dp))/ &
               (modal_raw(1) + modal_raw(2) + occ_curv(3) + occ_curv(4))

      response_charge_leak = 0.0_dp
      raw_charge_leak = 0.0_dp
      DO a = 1, 2
         response_charge_leak = response_charge_leak - &
                                2.0_dp*projector_weight*modal_curv(a)*(split_eval(a) - mu_response)
         raw_charge_leak = raw_charge_leak - &
                           2.0_dp*projector_weight*modal_curv(a)*(split_eval(a) - mu_raw)
      END DO
      DO i = 3, n
         response_charge_leak = response_charge_leak - &
                                2.0_dp*projector_weight*occ_curv(i)* &
                                (REAL(dchc(i, i), KIND=dp) - mu_response)
         raw_charge_leak = raw_charge_leak - &
                           2.0_dp*projector_weight*occ_curv(i)* &
                           (REAL(dchc(i, i), KIND=dp) - mu_raw)
      END DO
      charge_leak = ABS(REAL(docc_safe(1, 1) + docc_safe(2, 2) + docc_safe(3, 3) + &
                             docc_safe(4, 4), KIND=dp))
      response_charge_leak = ABS(response_charge_leak)
      raw_charge_leak = ABS(raw_charge_leak)
      response_gap = ABS(mu_response - mu_raw)

      IF (charge_leak > tol .OR. response_charge_leak > tol .OR. &
          raw_charge_leak < active_tol .OR. response_gap < active_tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "failed with charge/response/raw/gap", charge_leak, response_charge_leak, &
            raw_charge_leak, response_gap
      ELSE
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "OK with charge/response/raw/gap", charge_leak, response_charge_leak, &
            raw_charge_leak, response_gap
      END IF

   END SUBROUTINE check_history_safe_occ_projector_modal_mu_closure

! **************************************************************************************************
!> \brief Check that active/inert cluster modes stay gauge bounded.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_history_safe_occ_projector_modal_curvature_gauge_bound(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 3
      REAL(KIND=dp), PARAMETER                           :: active_tol = 1.0E-5_dp, &
                                                            ratio_limit = 0.10_dp

      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: evec
      INTEGER                                            :: a, b
      REAL(KIND=dp)                                      :: mean_curv, mean_gap, modal_max, raw_max, &
                                                            raw_modal, raw_ratio, unit
      REAL(KIND=dp), DIMENSION(n)                        :: modal_curv, occ_curv

      unit = 1.0_dp/SQRT(2.0_dp)
      evec(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      evec(1, 1) = CMPLX(unit, 0.0_dp, KIND=dp)
      evec(2, 1) = CMPLX(unit, 0.0_dp, KIND=dp)
      evec(1, 2) = CMPLX(unit, 0.0_dp, KIND=dp)
      evec(2, 2) = CMPLX(-unit, 0.0_dp, KIND=dp)
      evec(3, 3) = CMPLX(1.0_dp, 0.0_dp, KIND=dp)

      occ_curv(:) = [0.038_dp, 1.0E-5_dp, 0.004_dp]
      CALL dense_fixed_n_projector_modal_curvature(evec, occ_curv, modal_curv)

      raw_max = 0.0_dp
      DO a = 1, n
         raw_modal = 0.0_dp
         DO b = 1, n
            raw_modal = raw_modal + ABS(evec(b, a))**2*MAX(0.0_dp, occ_curv(b))
         END DO
         raw_max = MAX(raw_max, raw_modal)
      END DO
      mean_curv = SUM(MAX(0.0_dp, occ_curv(:)))/REAL(n, KIND=dp)
      modal_max = MAXVAL(modal_curv(:))
      raw_ratio = modal_max/MAX(raw_max, EPSILON(raw_max))
      mean_gap = mean_curv - modal_max

      IF (modal_max < active_tol .OR. raw_ratio >= ratio_limit .OR. mean_gap < active_tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "failed with modal/rawratio/meangap/raw", modal_max, raw_ratio, mean_gap, raw_max
      ELSE
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "OK with modal/rawratio/meangap/raw", modal_max, raw_ratio, mean_gap, raw_max
      END IF

   END SUBROUTINE check_history_safe_occ_projector_modal_curvature_gauge_bound

! **************************************************************************************************
!> \brief Check that inert admixtures do not collapse active modal curvature.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_history_safe_occ_projector_modal_curvature_active_tail(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 2, n3 = 3
      REAL(KIND=dp), PARAMETER :: active_fraction = 0.99_dp, active_pair_high_fraction = 0.50_dp, &
         active_pair_inert_fraction = 0.01_dp, active_pair_mid_fraction = 0.49_dp, &
         active_pair_tol = 0.98_dp, active_tol = 0.98_dp, dominant_fraction = 0.75_dp, &
         dominant_tol = 0.85_dp, edge_fraction = 0.55_dp, edge_tol = 0.70_dp, &
         mixed_fraction = 0.65_dp, mixed_tol = 0.80_dp, moderate_fraction = 0.90_dp, &
         moderate_tol = 0.80_dp, raw_tol = 1.0E-12_dp

      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: evec
      COMPLEX(KIND=dp), DIMENSION(n3, n3)                :: evec3
      REAL(KIND=dp) :: active_pair_ratio, active_pair_raw, active_ratio, dominant_ratio, &
         edge_ratio, inert_ratio, mixed_ratio, moderate_ratio, norm_active_pair, raw_active, &
         raw_inert, unit_active, unit_inert, v_high, v_inert, v_mid
      REAL(KIND=dp), DIMENSION(n)                        :: modal_curv, occ_curv
      REAL(KIND=dp), DIMENSION(n3)                       :: modal_curv3, occ_curv3

      unit_active = SQRT(active_fraction)
      unit_inert = SQRT(1.0_dp - active_fraction)
      evec(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      evec(1, 1) = CMPLX(unit_active, 0.0_dp, KIND=dp)
      evec(2, 1) = CMPLX(unit_inert, 0.0_dp, KIND=dp)
      evec(1, 2) = CMPLX(-unit_inert, 0.0_dp, KIND=dp)
      evec(2, 2) = CMPLX(unit_active, 0.0_dp, KIND=dp)

      occ_curv(:) = [0.038_dp, 1.0E-5_dp]
      CALL dense_fixed_n_projector_modal_curvature(evec, occ_curv, modal_curv)

      raw_active = active_fraction*occ_curv(1) + (1.0_dp - active_fraction)*occ_curv(2)
      raw_inert = (1.0_dp - active_fraction)*occ_curv(1) + active_fraction*occ_curv(2)
      active_ratio = modal_curv(1)/MAX(raw_active, EPSILON(raw_active))
      inert_ratio = modal_curv(2)/MAX(raw_inert, EPSILON(raw_inert))

      IF (active_ratio < active_tol .OR. active_ratio > 1.0_dp + raw_tol .OR. &
          inert_ratio > 0.25_dp .OR. &
          modal_curv(1) <= modal_curv(2) .OR. raw_active <= 0.0_dp .OR. &
          raw_inert <= 0.0_dp) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "failed with active/inert/mactive/minert", active_ratio, inert_ratio, &
            modal_curv(1), modal_curv(2)
      ELSE
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "OK with active/inert/mactive/minert", active_ratio, inert_ratio, &
            modal_curv(1), modal_curv(2)
      END IF

      v_high = SQRT(active_pair_high_fraction)
      v_mid = SQRT(active_pair_mid_fraction)
      v_inert = SQRT(active_pair_inert_fraction)
      norm_active_pair = SQRT(active_pair_high_fraction + active_pair_mid_fraction)
      evec3(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      evec3(1, 1) = CMPLX(v_high, 0.0_dp, KIND=dp)
      evec3(2, 1) = CMPLX(v_mid, 0.0_dp, KIND=dp)
      evec3(3, 1) = CMPLX(v_inert, 0.0_dp, KIND=dp)
      evec3(1, 2) = CMPLX(v_mid/norm_active_pair, 0.0_dp, KIND=dp)
      evec3(2, 2) = CMPLX(-v_high/norm_active_pair, 0.0_dp, KIND=dp)
      evec3(1, 3) = CMPLX(v_inert*v_high/norm_active_pair, 0.0_dp, KIND=dp)
      evec3(2, 3) = CMPLX(v_inert*v_mid/norm_active_pair, 0.0_dp, KIND=dp)
      evec3(3, 3) = CMPLX(-norm_active_pair, 0.0_dp, KIND=dp)

      occ_curv3(:) = [0.038_dp, 0.012_dp, 1.0E-5_dp]
      CALL dense_fixed_n_projector_modal_curvature(evec3, occ_curv3, modal_curv3)

      active_pair_raw = active_pair_high_fraction*occ_curv3(1) + &
                        active_pair_mid_fraction*occ_curv3(2) + &
                        active_pair_inert_fraction*occ_curv3(3)
      active_pair_ratio = modal_curv3(1)/MAX(active_pair_raw, EPSILON(active_pair_raw))

      IF (active_pair_ratio < active_pair_tol .OR. active_pair_ratio > 1.0_dp + raw_tol .OR. &
          modal_curv3(3) > 0.25_dp*modal_curv3(1) .OR. active_pair_raw <= 0.0_dp) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "failed with active-pair/inert/mactive/minert", active_pair_ratio, &
            modal_curv3(3)/MAX(modal_curv3(1), EPSILON(modal_curv3(1))), modal_curv3(1), &
            modal_curv3(3)
      ELSE
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "OK with active-pair/inert/mactive/minert", active_pair_ratio, &
            modal_curv3(3)/MAX(modal_curv3(1), EPSILON(modal_curv3(1))), modal_curv3(1), &
            modal_curv3(3)
      END IF

      unit_active = SQRT(moderate_fraction)
      unit_inert = SQRT(1.0_dp - moderate_fraction)
      evec(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      evec(1, 1) = CMPLX(unit_active, 0.0_dp, KIND=dp)
      evec(2, 1) = CMPLX(unit_inert, 0.0_dp, KIND=dp)
      evec(1, 2) = CMPLX(-unit_inert, 0.0_dp, KIND=dp)
      evec(2, 2) = CMPLX(unit_active, 0.0_dp, KIND=dp)

      CALL dense_fixed_n_projector_modal_curvature(evec, occ_curv, modal_curv)

      raw_active = moderate_fraction*occ_curv(1) + (1.0_dp - moderate_fraction)*occ_curv(2)
      raw_inert = (1.0_dp - moderate_fraction)*occ_curv(1) + moderate_fraction*occ_curv(2)
      moderate_ratio = modal_curv(1)/MAX(raw_active, EPSILON(raw_active))
      inert_ratio = modal_curv(2)/MAX(raw_inert, EPSILON(raw_inert))

      IF (moderate_ratio < moderate_tol .OR. moderate_ratio > 1.0_dp + raw_tol .OR. &
          inert_ratio > 0.25_dp .OR. &
          modal_curv(1) <= modal_curv(2) .OR. raw_active <= 0.0_dp .OR. &
          raw_inert <= 0.0_dp) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "failed with moderate/inert/mactive/minert", moderate_ratio, inert_ratio, &
            modal_curv(1), modal_curv(2)
      ELSE
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "OK with moderate/inert/mactive/minert", moderate_ratio, inert_ratio, &
            modal_curv(1), modal_curv(2)
      END IF

      unit_active = SQRT(dominant_fraction)
      unit_inert = SQRT(1.0_dp - dominant_fraction)
      evec(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      evec(1, 1) = CMPLX(unit_active, 0.0_dp, KIND=dp)
      evec(2, 1) = CMPLX(unit_inert, 0.0_dp, KIND=dp)
      evec(1, 2) = CMPLX(-unit_inert, 0.0_dp, KIND=dp)
      evec(2, 2) = CMPLX(unit_active, 0.0_dp, KIND=dp)

      CALL dense_fixed_n_projector_modal_curvature(evec, occ_curv, modal_curv)

      raw_active = dominant_fraction*occ_curv(1) + (1.0_dp - dominant_fraction)*occ_curv(2)
      raw_inert = (1.0_dp - dominant_fraction)*occ_curv(1) + dominant_fraction*occ_curv(2)
      dominant_ratio = modal_curv(1)/MAX(raw_active, EPSILON(raw_active))
      inert_ratio = modal_curv(2)/MAX(raw_inert, EPSILON(raw_inert))

      IF (dominant_ratio < dominant_tol .OR. dominant_ratio > 1.0_dp + raw_tol .OR. &
          inert_ratio > 0.25_dp .OR. &
          modal_curv(1) <= modal_curv(2) .OR. raw_active <= 0.0_dp .OR. &
          raw_inert <= 0.0_dp) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "failed with dominant/inert/mactive/minert", dominant_ratio, inert_ratio, &
            modal_curv(1), modal_curv(2)
      ELSE
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "OK with dominant/inert/mactive/minert", dominant_ratio, inert_ratio, &
            modal_curv(1), modal_curv(2)
      END IF

      unit_active = SQRT(mixed_fraction)
      unit_inert = SQRT(1.0_dp - mixed_fraction)
      evec(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      evec(1, 1) = CMPLX(unit_active, 0.0_dp, KIND=dp)
      evec(2, 1) = CMPLX(unit_inert, 0.0_dp, KIND=dp)
      evec(1, 2) = CMPLX(-unit_inert, 0.0_dp, KIND=dp)
      evec(2, 2) = CMPLX(unit_active, 0.0_dp, KIND=dp)

      CALL dense_fixed_n_projector_modal_curvature(evec, occ_curv, modal_curv)

      raw_active = mixed_fraction*occ_curv(1) + (1.0_dp - mixed_fraction)*occ_curv(2)
      raw_inert = (1.0_dp - mixed_fraction)*occ_curv(1) + mixed_fraction*occ_curv(2)
      mixed_ratio = modal_curv(1)/MAX(raw_active, EPSILON(raw_active))
      inert_ratio = modal_curv(2)/MAX(raw_inert, EPSILON(raw_inert))

      IF (mixed_ratio < mixed_tol .OR. mixed_ratio > 1.0_dp + raw_tol .OR. &
          inert_ratio > 0.25_dp .OR. &
          modal_curv(1) <= modal_curv(2) .OR. raw_active <= 0.0_dp .OR. &
          raw_inert <= 0.0_dp) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "failed with mixed/inert/mactive/minert", mixed_ratio, inert_ratio, &
            modal_curv(1), modal_curv(2)
      ELSE
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "OK with mixed/inert/mactive/minert", mixed_ratio, inert_ratio, &
            modal_curv(1), modal_curv(2)
      END IF

      unit_active = SQRT(edge_fraction)
      unit_inert = SQRT(1.0_dp - edge_fraction)
      evec(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      evec(1, 1) = CMPLX(unit_active, 0.0_dp, KIND=dp)
      evec(2, 1) = CMPLX(unit_inert, 0.0_dp, KIND=dp)
      evec(1, 2) = CMPLX(-unit_inert, 0.0_dp, KIND=dp)
      evec(2, 2) = CMPLX(unit_active, 0.0_dp, KIND=dp)

      CALL dense_fixed_n_projector_modal_curvature(evec, occ_curv, modal_curv)

      raw_active = edge_fraction*occ_curv(1) + (1.0_dp - edge_fraction)*occ_curv(2)
      raw_inert = (1.0_dp - edge_fraction)*occ_curv(1) + edge_fraction*occ_curv(2)
      edge_ratio = modal_curv(1)/MAX(raw_active, EPSILON(raw_active))
      inert_ratio = modal_curv(2)/MAX(raw_inert, EPSILON(raw_inert))

      IF (edge_ratio < edge_tol .OR. edge_ratio > 1.0_dp + raw_tol .OR. &
          inert_ratio > 0.25_dp .OR. &
          modal_curv(1) <= modal_curv(2) .OR. raw_active <= 0.0_dp .OR. &
          raw_inert <= 0.0_dp) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "failed with edge/inert/mactive/minert", edge_ratio, inert_ratio, &
            modal_curv(1), modal_curv(2)
      ELSE
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "OK with edge/inert/mactive/minert", edge_ratio, inert_ratio, &
            modal_curv(1), modal_curv(2)
      END IF

   END SUBROUTINE check_history_safe_occ_projector_modal_curvature_active_tail

! **************************************************************************************************
!> \brief Check that one gauge jump in a near-degenerate cluster does not bias physical dP slopes.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_history_safe_occ_projector_cluster_slope_split(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 4
      REAL(KIND=dp), PARAMETER                           :: active_tol = 1.0E-4_dp, &
                                                            fd_step = 1.0E-7_dp, &
                                                            gauge_gap = 5.0E-5_dp, &
                                                            projector_weight = 0.37_dp

      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: chc, dchc, docc_fd, docc_raw, docc_ref, &
                                                            docc_safe, occ_mat_m, occ_mat_p
      INTEGER                                            :: i
      REAL(KIND=dp)                                      :: cluster_avg, cluster_center, err, &
                                                            fd_err, physical_slope, raw_noise, &
                                                            retained
      REAL(KIND=dp), DIMENSION(n)                        :: evals, occ_curv, occ_expected, occ_noisy

      evals(:) = [-0.010000_dp, -0.009982_dp, -0.009964_dp, 0.42_dp]
      chc(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      DO i = 1, n
         chc(i, i) = CMPLX(evals(i), 0.0_dp, KIND=dp)
      END DO

      dchc(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      dchc(1, 2) = CMPLX(0.28_dp, -0.19_dp, KIND=dp)
      dchc(2, 3) = CMPLX(-0.16_dp, 0.13_dp, KIND=dp)
      dchc(1, 3) = CMPLX(0.10_dp, 0.07_dp, KIND=dp)
      DO i = 2, n
         dchc(i, 1:i - 1) = CONJG(dchc(1:i - 1, i))
      END DO

      occ_curv(:) = [0.038_dp, 0.035_dp, 0.032_dp, 0.006_dp]
      cluster_center = SUM(evals(1:3))/3.0_dp
      cluster_avg = 1.10_dp
      physical_slope = -0.5_dp*dense_shared_occ_slope(occ_curv(2), occ_curv(3))
      occ_expected(:) = [cluster_avg + physical_slope*(evals(1) - cluster_center), &
                         cluster_avg + physical_slope*(evals(2) - cluster_center), &
                         cluster_avg + physical_slope*(evals(3) - cluster_center), &
                         0.18_dp]
      occ_noisy(:) = occ_expected(:)
      occ_noisy(1) = occ_noisy(1) + 0.23_dp

      CALL dense_history_safe_occ_projector_derivative_chc(chc, dchc, occ_noisy, projector_weight, &
                                                           docc_safe, occ_curv, gauge_gap)
      CALL dense_spectral_occ_projector_derivative_chc(chc, dchc, occ_expected, docc_ref, &
                                                       occ_curv, gauge_gap)
      docc_ref(:, :) = projector_weight*docc_ref(:, :)
      CALL dense_spectral_occ_projector_derivative_chc(chc, dchc, occ_noisy, docc_raw)
      CALL dense_history_safe_occ_projector_from_chc(chc + fd_step*dchc, occ_noisy, projector_weight, &
                                                     occ_mat_p, occ_curv, gauge_gap)
      CALL dense_history_safe_occ_projector_from_chc(chc - fd_step*dchc, occ_noisy, projector_weight, &
                                                     occ_mat_m, occ_curv, gauge_gap)
      docc_fd(:, :) = (occ_mat_p(:, :) - occ_mat_m(:, :))/(2.0_dp*fd_step)

      err = MAXVAL(ABS(docc_safe(:, :) - docc_ref(:, :)))
      fd_err = MAXVAL(ABS(docc_safe(:, :) - docc_fd(:, :)))
      err = MAX(err, fd_err)
      raw_noise = MAXVAL(ABS(projector_weight*docc_raw(:, :) - docc_ref(:, :)))
      retained = MAXVAL(ABS(docc_safe(:, :)))

      IF (raw_noise < active_tol .OR. retained < active_tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "failed with error/fderr/rawnoise/retained", err, fd_err, raw_noise, retained
      ELSE
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "OK with error/fderr/rawnoise/retained", err, fd_err, raw_noise, retained
      END IF

   END SUBROUTINE check_history_safe_occ_projector_cluster_slope_split

! **************************************************************************************************
!> \brief Check that one noisy near-degenerate pair filters the full connected component.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_history_safe_occ_projector_component_filter(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 4
      REAL(KIND=dp), PARAMETER                           :: active_tol = 1.0E-4_dp, &
                                                            fd_step = 1.0E-7_dp, &
                                                            gauge_gap = 5.0E-5_dp, &
                                                            projector_weight = 0.37_dp

      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: chc, dchc, docc_fd, docc_legacy_fd, &
                                                            docc_safe, occ_mat_m, occ_mat_p
      INTEGER                                            :: i
      REAL(KIND=dp)                                      :: component_err, component_span, &
                                                            eval_center, legacy_err, legacy_span, &
                                                            occ_tol, pair_plateau, pair_slope
      REAL(KIND=dp), DIMENSION(n)                        :: evals, occ_component, occ_curv, &
                                                            occ_pair_local, occ_scale

      chc(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      evals(:) = [-0.010000_dp, -0.009982_dp, -0.009960_dp, 0.42_dp]
      DO i = 1, n
         chc(i, i) = CMPLX(evals(i), 0.0_dp, KIND=dp)
      END DO
      dchc(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      dchc(2, 3) = CMPLX(0.18_dp, -0.11_dp, KIND=dp)
      dchc(3, 2) = CONJG(dchc(2, 3))

      occ_scale(:) = [1.20_dp, 0.40_dp, 0.399999_dp, 0.11_dp]
      occ_curv(:) = [0.038_dp, 0.035_dp, 0.018_dp, 0.006_dp]
      occ_tol = SQRT(EPSILON(1.0_dp))*MAX(1.0_dp, MAXVAL(ABS(occ_scale(:))))

      occ_component(:) = occ_scale(:)
      CALL dense_near_degenerate_occupation_plateau(evals, occ_curv, gauge_gap, occ_tol, occ_component)

      occ_pair_local(:) = occ_scale(:)
      pair_plateau = 0.5_dp*(occ_pair_local(1) + occ_pair_local(2))
      eval_center = 0.5_dp*(evals(1) + evals(2))
      pair_slope = -dense_shared_occ_slope(occ_curv(1), occ_curv(2))
      occ_pair_local(1) = pair_plateau + pair_slope*(evals(1) - eval_center)
      occ_pair_local(2) = pair_plateau + pair_slope*(evals(2) - eval_center)

      CALL dense_history_safe_occ_projector_derivative_chc(chc, dchc, occ_scale, projector_weight, &
                                                           docc_safe, occ_curv, gauge_gap)
      CALL dense_spectral_occ_projector_from_chc(chc + fd_step*dchc, occ_component, occ_mat_p)
      CALL dense_spectral_occ_projector_from_chc(chc - fd_step*dchc, occ_component, occ_mat_m)
      docc_fd(:, :) = projector_weight*(occ_mat_p(:, :) - occ_mat_m(:, :))/(2.0_dp*fd_step)
      component_err = MAXVAL(ABS(docc_safe(:, :) - docc_fd(:, :)))

      CALL dense_spectral_occ_projector_from_chc(chc + fd_step*dchc, occ_pair_local, occ_mat_p)
      CALL dense_spectral_occ_projector_from_chc(chc - fd_step*dchc, occ_pair_local, occ_mat_m)
      docc_legacy_fd(:, :) = projector_weight*(occ_mat_p(:, :) - occ_mat_m(:, :))/(2.0_dp*fd_step)
      legacy_err = MAXVAL(ABS(docc_safe(:, :) - docc_legacy_fd(:, :)))
      component_span = MAXVAL(occ_component(1:3)) - MINVAL(occ_component(1:3))
      legacy_span = MAXVAL(occ_pair_local(1:3)) - MINVAL(occ_pair_local(1:3))

      IF (legacy_err < active_tol .OR. component_span >= 0.20_dp*legacy_span) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "failed with error/legacy/span/legacy_span", component_err, legacy_err, &
            component_span, legacy_span
      ELSE
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "OK with error/legacy/span/legacy_span", component_err, legacy_err, &
            component_span, legacy_span
      END IF

   END SUBROUTINE check_history_safe_occ_projector_component_filter

! **************************************************************************************************
!> \brief Check that a fully noisy near-degenerate cluster keeps fixed-N Mermin shape.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_history_safe_occ_projector_fixed_n_cluster_shape(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 4
      REAL(KIND=dp), PARAMETER                           :: active_tol = 1.0E-4_dp, &
                                                            fd_step = 1.0E-7_dp, &
                                                            gauge_gap = 5.0E-5_dp, &
                                                            projector_weight = 0.37_dp

      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: chc, dchc, docc_fd, docc_ref, docc_safe, &
                                                            occ_mat_m, occ_mat_p
      INTEGER                                            :: i
      REAL(KIND=dp)                                      :: active, cluster_sum, err, occ_tol, &
                                                            shape_err, shape_sum, span, total_err
      REAL(KIND=dp), DIMENSION(n)                        :: evals, expected, occ_curv, occ_filtered, &
                                                            occ_scale

      evals(:) = [-0.010000_dp, -0.009982_dp, -0.009960_dp, 0.42_dp]
      chc(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      DO i = 1, n
         chc(i, i) = CMPLX(evals(i), 0.0_dp, KIND=dp)
      END DO

      dchc(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      dchc(1, 3) = CMPLX(0.19_dp, -0.14_dp, KIND=dp)
      dchc(3, 1) = CONJG(dchc(1, 3))

      occ_scale(:) = [1.20_dp, 0.41_dp, 0.10_dp, 0.08_dp]
      occ_curv(:) = [0.038_dp, 0.035_dp, 0.018_dp, 0.006_dp]
      occ_tol = SQRT(EPSILON(1.0_dp))*MAX(1.0_dp, MAXVAL(ABS(occ_scale(:))))
      CALL dense_monotone_occupation_slopes(occ_scale, occ_tol, occ_filtered)
      cluster_sum = SUM(occ_filtered(1:3))
      CALL dense_near_degenerate_occupation_plateau(evals, occ_curv, gauge_gap, occ_tol, occ_filtered)

      expected(1) = 0.0_dp
      DO i = 2, 3
         expected(i) = expected(i - 1) - &
                       dense_shared_occ_slope(occ_curv(i - 1), occ_curv(i))*(evals(i) - evals(i - 1))
      END DO
      shape_sum = SUM(expected(1:3))
      expected(1:3) = expected(1:3) + (cluster_sum - shape_sum)/3.0_dp
      expected(4) = occ_scale(4)
      total_err = ABS(SUM(occ_filtered(1:3)) - cluster_sum)
      shape_err = MAXVAL(ABS(occ_filtered(1:3) - expected(1:3)))
      span = MAXVAL(occ_filtered(1:3)) - MINVAL(occ_filtered(1:3))

      CALL dense_history_safe_occ_projector_derivative_chc(chc, dchc, occ_scale, projector_weight, &
                                                           docc_safe, occ_curv, gauge_gap)
      CALL dense_spectral_occ_projector_derivative_chc(chc, dchc, occ_filtered, docc_ref, &
                                                       occ_curv, gauge_gap)
      docc_ref(:, :) = projector_weight*docc_ref(:, :)
      CALL dense_history_safe_occ_projector_from_chc(chc + fd_step*dchc, occ_scale, projector_weight, &
                                                     occ_mat_p, occ_curv, gauge_gap)
      CALL dense_history_safe_occ_projector_from_chc(chc - fd_step*dchc, occ_scale, projector_weight, &
                                                     occ_mat_m, occ_curv, gauge_gap)
      docc_fd(:, :) = (occ_mat_p(:, :) - occ_mat_m(:, :))/(2.0_dp*fd_step)
      err = MAX(MAXVAL(ABS(docc_safe(:, :) - docc_ref(:, :))), shape_err)
      active = MAXVAL(ABS(docc_safe(:, :)))

      IF (total_err > tol .OR. active < active_tol .OR. &
          MAXVAL(ABS(docc_fd(:, :) - docc_ref(:, :))) < active_tol .OR. span < 1.0E-7_dp) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,5(1X,ES12.5))') TRIM(label), &
            "failed with error/total/span/active/shape", err, total_err, span, active, shape_err
      ELSE
         WRITE (*, '(A,1X,A,5(1X,ES12.5))') TRIM(label), &
            "OK with error/total/span/active/shape", err, total_err, span, active, shape_err
      END IF

   END SUBROUTINE check_history_safe_occ_projector_fixed_n_cluster_shape

! **************************************************************************************************
!> \brief Check that the history-safe dP model is covariant under a complex REF-basis gauge.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_history_safe_occ_projector_covariant_fd(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 4
      REAL(KIND=dp), PARAMETER                           :: active_tol = 1.0E-4_dp, &
                                                            fd_step = 1.0E-7_dp, &
                                                            gauge_gap = 5.0E-5_dp, &
                                                            projector_weight = 0.37_dp

      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: chc, chc_eig, dchc, dchc_eig, docc_back, &
                                                            docc_fd, docc_raw, docc_ref, docc_rot, &
                                                            occ_mat_m, occ_mat_p, q
      INTEGER                                            :: i
      REAL(KIND=dp)                                      :: active, cov_err, diag_gap, fd_err, &
                                                            occ_sum, raw_noise, ref_fd_err, &
                                                            shift_mu
      REAL(KIND=dp), DIMENSION(n)                        :: evals, expected_diag, occ_curv, &
                                                            occ_scale, shift
      REAL(KIND=dp), DIMENSION(n, n)                     :: x_im, x_re

      evals(:) = [-0.010000_dp, -0.009982_dp, -0.009960_dp, 0.39_dp]
      occ_scale(:) = [1.22_dp, 0.47_dp, 0.13_dp, 0.05_dp]
      occ_curv(:) = [0.038_dp, 0.035_dp, 0.032_dp, 0.006_dp]

      chc_eig(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      dchc_eig(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      DO i = 1, n
         chc_eig(i, i) = CMPLX(evals(i), 0.0_dp, KIND=dp)
      END DO
      dchc_eig(1, 1) = CMPLX(0.19_dp, 0.0_dp, KIND=dp)
      dchc_eig(2, 2) = CMPLX(-0.07_dp, 0.0_dp, KIND=dp)
      dchc_eig(3, 3) = CMPLX(0.04_dp, 0.0_dp, KIND=dp)
      dchc_eig(4, 4) = CMPLX(0.02_dp, 0.0_dp, KIND=dp)
      dchc_eig(1, 2) = CMPLX(0.24_dp, -0.18_dp, KIND=dp)
      dchc_eig(2, 3) = CMPLX(-0.16_dp, 0.11_dp, KIND=dp)
      dchc_eig(1, 3) = CMPLX(0.09_dp, 0.08_dp, KIND=dp)
      dchc_eig(2, 4) = CMPLX(-0.04_dp, 0.03_dp, KIND=dp)
      DO i = 2, n
         dchc_eig(i, 1:i - 1) = CONJG(dchc_eig(1:i - 1, i))
      END DO
      DO i = 1, n
         shift(i) = REAL(dchc_eig(i, i), KIND=dp)
      END DO
      occ_sum = SUM(occ_curv(:))
      shift_mu = DOT_PRODUCT(occ_curv(:), shift(:))/occ_sum
      DO i = 1, n
         expected_diag(i) = -2.0_dp*projector_weight*occ_curv(i)*(shift(i) - shift_mu)
      END DO

      x_re(:, :) = 0.0_dp
      x_re(1, 2) = 0.21_dp
      x_re(1, 3) = -0.11_dp
      x_re(2, 3) = 0.07_dp
      x_re(1, 4) = 0.03_dp
      x_re(2, 4) = -0.02_dp
      DO i = 2, n
         x_re(i, 1:i - 1) = -x_re(1:i - 1, i)
      END DO
      x_im(:, :) = 0.0_dp
      x_im(1, 1) = 0.03_dp
      x_im(2, 2) = -0.02_dp
      x_im(3, 3) = 0.01_dp
      x_im(4, 4) = -0.01_dp
      x_im(1, 2) = -0.04_dp
      x_im(1, 3) = 0.02_dp
      x_im(2, 3) = 0.03_dp
      x_im(3, 4) = -0.02_dp
      DO i = 2, n
         x_im(i, 1:i - 1) = x_im(1:i - 1, i)
      END DO
      CALL dense_rotation_u(x_re, x_im, q)

      chc(:, :) = MATMUL(q, MATMUL(chc_eig, CONJG(TRANSPOSE(q))))
      dchc(:, :) = MATMUL(q, MATMUL(dchc_eig, CONJG(TRANSPOSE(q))))

      CALL dense_history_safe_occ_projector_derivative_chc(chc_eig, dchc_eig, occ_scale, &
                                                           projector_weight, docc_ref, &
                                                           occ_curv, gauge_gap)
      CALL dense_history_safe_occ_projector_derivative_chc(chc, dchc, occ_scale, &
                                                           projector_weight, docc_rot, &
                                                           occ_curv, gauge_gap)
      docc_back(:, :) = MATMUL(CONJG(TRANSPOSE(q)), MATMUL(docc_rot, q))
      cov_err = MAXVAL(ABS(docc_back(:, :) - docc_ref(:, :)))

      CALL dense_history_safe_occ_projector_from_chc(chc_eig + fd_step*dchc_eig, occ_scale, &
                                                     projector_weight, occ_mat_p, occ_curv, gauge_gap)
      CALL dense_history_safe_occ_projector_from_chc(chc_eig - fd_step*dchc_eig, occ_scale, &
                                                     projector_weight, occ_mat_m, occ_curv, gauge_gap)
      docc_fd(:, :) = (occ_mat_p(:, :) - occ_mat_m(:, :))/(2.0_dp*fd_step)
      diag_gap = 0.0_dp
      DO i = 1, n
         diag_gap = MAX(diag_gap, ABS(docc_fd(i, i) - CMPLX(expected_diag(i), 0.0_dp, KIND=dp)))
         docc_fd(i, i) = CMPLX(expected_diag(i), 0.0_dp, KIND=dp)
      END DO
      ref_fd_err = MAXVAL(ABS(docc_ref(:, :) - docc_fd(:, :)))

      CALL dense_history_safe_occ_projector_from_chc(chc + fd_step*dchc, occ_scale, projector_weight, &
                                                     occ_mat_p, occ_curv, gauge_gap)
      CALL dense_history_safe_occ_projector_from_chc(chc - fd_step*dchc, occ_scale, projector_weight, &
                                                     occ_mat_m, occ_curv, gauge_gap)
      docc_fd(:, :) = (occ_mat_p(:, :) - occ_mat_m(:, :))/(2.0_dp*fd_step)
      docc_back(:, :) = MATMUL(CONJG(TRANSPOSE(q)), MATMUL(docc_fd, q))
      DO i = 1, n
         docc_back(i, i) = CMPLX(expected_diag(i), 0.0_dp, KIND=dp)
      END DO
      fd_err = MAXVAL(ABS(docc_ref(:, :) - docc_back(:, :)))

      CALL dense_spectral_occ_projector_derivative_chc(chc_eig, dchc_eig, occ_scale, docc_raw)
      raw_noise = MAXVAL(ABS(projector_weight*docc_raw(:, :) - docc_ref(:, :)))
      active = MAXVAL(ABS(docc_ref(:, :)))

      IF (cov_err > tol .OR. diag_gap < active_tol .OR. raw_noise < active_tol .OR. &
          active < active_tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,6(1X,ES12.5))') TRIM(label), &
            "failed with cov/fd/reffd/diaggap/rawnoise/active", cov_err, fd_err, ref_fd_err, &
            diag_gap, raw_noise, active
      ELSE
         WRITE (*, '(A,1X,A,6(1X,ES12.5))') TRIM(label), &
            "OK with cov/fd/reffd/diaggap/rawnoise/active", cov_err, fd_err, ref_fd_err, &
            diag_gap, raw_noise, active
      END IF

   END SUBROUTINE check_history_safe_occ_projector_covariant_fd

! **************************************************************************************************
!> \brief Check that near-degenerate cluster-to-singleton blocks use a covariant divided difference.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_history_safe_occ_projector_cluster_offblock_response(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 4
      REAL(KIND=dp), PARAMETER                           :: gauge_gap = 5.0E-5_dp, &
                                                            projector_weight = 0.37_dp

      COMPLEX(KIND=dp), DIMENSION(2)                     :: cluster_col, expected_col, modal_col
      COMPLEX(KIND=dp), DIMENSION(2, 2)                  :: cluster_evec, cluster_work
      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: kernel, work
      INTEGER                                            :: a, b, i
      REAL(KIND=dp)                                      :: coeff, component_weight, denom, err, &
                                                            herm_err, old_coeff, old_gap, &
                                                            scalar_gap, slope_cap
      REAL(KIND=dp), DIMENSION(2)                        :: cluster_eval, modal_curv, modal_energy, &
                                                            modal_occ
      REAL(KIND=dp), DIMENSION(n)                        :: evals, occ_curv, occ_slope

      evals(:) = [-0.010000_dp, -0.009982_dp, 0.21_dp, 0.48_dp]
      occ_slope(:) = [0.82_dp, 0.78_dp, 0.17_dp, 0.04_dp]
      occ_curv(:) = [0.038_dp, 0.034_dp, 0.011_dp, 0.004_dp]

      work(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      work(1, 1) = CMPLX(0.06_dp, 0.0_dp, KIND=dp)
      work(2, 2) = CMPLX(-0.04_dp, 0.0_dp, KIND=dp)
      work(1, 2) = CMPLX(0.22_dp, -0.15_dp, KIND=dp)
      work(1, 3) = CMPLX(-0.18_dp, 0.09_dp, KIND=dp)
      work(2, 3) = CMPLX(0.07_dp, -0.05_dp, KIND=dp)
      work(3, 4) = CMPLX(0.03_dp, 0.02_dp, KIND=dp)
      DO i = 2, n
         work(i, 1:i - 1) = CONJG(work(1:i - 1, i))
      END DO

      kernel(:, :) = CMPLX(9.0_dp, -3.0_dp, KIND=dp)
      CALL dense_fixed_n_projector_cluster_kernel(evals, work, occ_curv, occ_slope, gauge_gap, &
                                                  projector_weight, kernel)

      cluster_work(:, :) = 0.5_dp*(work(1:2, 1:2) + CONJG(TRANSPOSE(work(1:2, 1:2))))
      CALL diag_complex(cluster_work, cluster_evec, cluster_eval)
      CALL dense_fixed_n_projector_modal_curvature(cluster_evec, occ_curv(1:2), modal_curv)
      modal_energy(:) = 0.0_dp
      DO a = 1, 2
         DO b = 1, 2
            component_weight = ABS(cluster_evec(b, a))**2
            modal_energy(a) = modal_energy(a) + component_weight*evals(b)
         END DO
      END DO
      CALL dense_projector_cluster_modal_occupation(evals(1:2), occ_curv(1:2), occ_slope(1:2), &
                                                    modal_energy, modal_occ, modal_curv)
      cluster_col(:) = work(1:2, 3)
      modal_col(:) = MATMUL(CONJG(TRANSPOSE(cluster_evec)), cluster_col)
      DO a = 1, 2
         denom = modal_energy(a) - evals(3)
         IF (ABS(denom) > SQRT(EPSILON(1.0_dp))*MAX(1.0_dp, MAXVAL(ABS(evals(:))))) THEN
            coeff = (modal_occ(a) - occ_slope(3))/denom
         ELSE
            coeff = 0.0_dp
         END IF
         coeff = MIN(0.0_dp, coeff)
         slope_cap = dense_projector_pair_slope_cap(denom, gauge_gap, modal_curv(a), occ_curv(3), &
                                                     MAXVAL(MAX(0.0_dp, occ_curv(:))))
         coeff = MAX(coeff, -slope_cap)
         modal_col(a) = coeff*modal_col(a)
      END DO
      expected_col(:) = projector_weight*MATMUL(cluster_evec, modal_col)

      old_coeff = ((occ_slope(1) + occ_slope(2))/2.0_dp - occ_slope(3))/ &
                  ((evals(1) + evals(2))/2.0_dp - evals(3))
      old_coeff = MIN(0.0_dp, old_coeff)
      old_coeff = MAX(old_coeff, -2.0_dp*MAX((occ_curv(1) + occ_curv(2))/2.0_dp, occ_curv(3)))

      err = MAX(ABS(kernel(1, 3) - expected_col(1)), &
                ABS(kernel(2, 3) - expected_col(2)))
      herm_err = MAX(ABS(kernel(3, 1) - CONJG(kernel(1, 3))), &
                     ABS(kernel(3, 2) - CONJG(kernel(2, 3))))
      old_gap = MAX(ABS(kernel(1, 3) - CMPLX(9.0_dp, -3.0_dp, KIND=dp)), &
                    ABS(kernel(2, 3) - CMPLX(9.0_dp, -3.0_dp, KIND=dp)))
      scalar_gap = MAX(ABS(expected_col(1) - projector_weight*old_coeff*work(1, 3)), &
                       ABS(expected_col(2) - projector_weight*old_coeff*work(2, 3)))

      IF (err > tol .OR. herm_err > tol .OR. old_gap < 1.0_dp .OR. scalar_gap < 1.0E-5_dp) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "failed with error/hermitian/oldgap/scalargap", err, herm_err, old_gap, scalar_gap
      ELSE
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "OK with error/hermitian/oldgap/scalargap", err, herm_err, old_gap, scalar_gap
      END IF

   END SUBROUTINE check_history_safe_occ_projector_cluster_offblock_response

! **************************************************************************************************
!> \brief Check that near thermal offblocks use the shared modal occupation slope cap.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_history_safe_occ_projector_cluster_offblock_modal_cap(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 4
      REAL(KIND=dp), PARAMETER                           :: active_tol = 1.0E-5_dp, &
                                                            gauge_gap = 5.0E-5_dp, &
                                                            projector_weight = 0.37_dp

      COMPLEX(KIND=dp), DIMENSION(2)                     :: cluster_col, expected_col, expected_modal_col, &
                                                            legacy_col, legacy_modal_col
      COMPLEX(KIND=dp), DIMENSION(2, 2)                  :: cluster_evec, cluster_work
      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: kernel, work
      INTEGER                                            :: a, b, i, iact
      REAL(KIND=dp)                                      :: cap_near, cap_old, cap_ratio, coeff, &
                                                            component_weight, denom, err, herm_err, raw_gap, &
                                                            slope_cap
      REAL(KIND=dp), DIMENSION(2)                        :: cluster_eval, modal_curv, modal_energy, &
                                                            modal_occ
      REAL(KIND=dp), DIMENSION(n)                        :: evals, occ_curv, occ_slope

      evals(:) = [-0.010000_dp, -0.009990_dp, -0.009910_dp, 0.42_dp]
      occ_slope(:) = [1.02_dp, 1.01_dp, 0.18_dp, 0.04_dp]
      occ_curv(:) = [0.038_dp, 1.0E-5_dp, 1.0E-5_dp, 0.004_dp]

      work(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      work(1, 1) = CMPLX(0.12_dp, 0.0_dp, KIND=dp)
      work(2, 2) = CMPLX(-0.08_dp, 0.0_dp, KIND=dp)
      work(3, 3) = CMPLX(0.03_dp, 0.0_dp, KIND=dp)
      work(1, 2) = CMPLX(2.0E-4_dp, -1.0E-4_dp, KIND=dp)
      work(1, 3) = CMPLX(0.21_dp, -0.06_dp, KIND=dp)
      work(2, 3) = CMPLX(0.03_dp, 0.01_dp, KIND=dp)
      DO i = 2, n
         work(i, 1:i - 1) = CONJG(work(1:i - 1, i))
      END DO

      kernel(:, :) = CMPLX(9.0_dp, -3.0_dp, KIND=dp)
      CALL dense_fixed_n_projector_cluster_kernel(evals, work, occ_curv, occ_slope, gauge_gap, &
                                                  projector_weight, kernel)

      cluster_work(:, :) = 0.5_dp*(work(1:2, 1:2) + CONJG(TRANSPOSE(work(1:2, 1:2))))
      CALL diag_complex(cluster_work, cluster_evec, cluster_eval)
      CALL dense_fixed_n_projector_modal_curvature(cluster_evec, occ_curv(1:2), modal_curv)
      modal_energy(:) = 0.0_dp
      DO a = 1, 2
         DO b = 1, 2
            component_weight = ABS(cluster_evec(b, a))**2
            modal_energy(a) = modal_energy(a) + component_weight*evals(b)
         END DO
      END DO
      CALL dense_projector_cluster_modal_occupation(evals(1:2), occ_curv(1:2), occ_slope(1:2), &
                                                    modal_energy, modal_occ, modal_curv)

      cluster_col(:) = work(1:2, 3)
      expected_modal_col(:) = MATMUL(CONJG(TRANSPOSE(cluster_evec)), cluster_col)
      legacy_modal_col(:) = expected_modal_col(:)
      DO a = 1, 2
         denom = modal_energy(a) - evals(3)
         IF (ABS(denom) > SQRT(EPSILON(1.0_dp))*MAX(1.0_dp, MAXVAL(ABS(evals(:))))) THEN
            coeff = (modal_occ(a) - occ_slope(3))/denom
         ELSE
            coeff = 0.0_dp
         END IF
         coeff = MIN(0.0_dp, coeff)
         slope_cap = dense_projector_pair_slope_cap(denom, gauge_gap, modal_curv(a), occ_curv(3), &
                                                     MAXVAL(MAX(0.0_dp, occ_curv(:))))
         expected_modal_col(a) = MAX(coeff, -slope_cap)*expected_modal_col(a)
         slope_cap = 2.0_dp*MAX(0.0_dp, MAX(modal_curv(a), occ_curv(3)))
         legacy_modal_col(a) = MAX(coeff, -slope_cap)*legacy_modal_col(a)
      END DO
      expected_col(:) = projector_weight*MATMUL(cluster_evec, expected_modal_col)
      legacy_col(:) = projector_weight*MATMUL(cluster_evec, legacy_modal_col)

      iact = 1
      IF (modal_curv(2) > modal_curv(1)) iact = 2
      denom = modal_energy(iact) - evals(3)
      cap_near = dense_projector_pair_slope_cap(denom, gauge_gap, modal_curv(iact), occ_curv(3), &
                                                MAXVAL(MAX(0.0_dp, occ_curv(:))))
      cap_old = 2.0_dp*MAX(0.0_dp, MAX(modal_curv(iact), occ_curv(3)))
      cap_ratio = 0.0_dp
      IF (cap_old > EPSILON(cap_old)) cap_ratio = cap_near/cap_old

      err = MAX(ABS(kernel(1, 3) - expected_col(1)), ABS(kernel(2, 3) - expected_col(2)))
      herm_err = MAX(ABS(kernel(3, 1) - CONJG(kernel(1, 3))), &
                     ABS(kernel(3, 2) - CONJG(kernel(2, 3))))
      raw_gap = MAXVAL(ABS(expected_col(:) - legacy_col(:)))

      IF (err > tol .OR. herm_err > tol .OR. raw_gap < active_tol .OR. cap_ratio > 0.10_dp) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "failed with err/herm/rawgap/capratio", err, herm_err, raw_gap, cap_ratio
      ELSE
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "OK with err/herm/rawgap/capratio", err, herm_err, raw_gap, cap_ratio
      END IF

   END SUBROUTINE check_history_safe_occ_projector_cluster_offblock_modal_cap

! **************************************************************************************************
!> \brief Check that near-degenerate cluster-to-cluster blocks use a covariant divided difference.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_history_safe_occ_projector_cluster_pair_response(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 4
      REAL(KIND=dp), PARAMETER                           :: gauge_gap = 5.0E-5_dp, &
                                                            projector_weight = 0.37_dp

      COMPLEX(KIND=dp), DIMENSION(2, 2)                  :: cluster_evec_a, cluster_evec_b, &
                                                            cluster_work_a, cluster_work_b, &
                                                            expected_block, modal_block, &
                                                            scalar_block
      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: kernel, work
      INTEGER                                            :: a, b, i
      REAL(KIND=dp) :: coeff, component_weight, denom, err, herm_err, old_gap, scalar_curv_a, &
         scalar_curv_b, scalar_energy_a, scalar_energy_b, scalar_gap, scalar_occ_a, scalar_occ_b, &
         slope_cap
      REAL(KIND=dp), DIMENSION(2) :: cluster_eval_a, cluster_eval_b, modal_curv_a, modal_curv_b, &
         modal_energy_a, modal_energy_b, modal_occ_a, modal_occ_b
      REAL(KIND=dp), DIMENSION(n)                        :: evals, occ_curv, occ_slope

      evals(:) = [-0.010000_dp, -0.009982_dp, 0.210000_dp, 0.210018_dp]
      occ_slope(:) = [0.82_dp, 0.78_dp, 0.17_dp, 0.14_dp]
      occ_curv(:) = [0.038_dp, 0.034_dp, 0.011_dp, 0.009_dp]

      work(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      work(1, 1) = CMPLX(0.06_dp, 0.0_dp, KIND=dp)
      work(2, 2) = CMPLX(-0.04_dp, 0.0_dp, KIND=dp)
      work(1, 2) = CMPLX(0.22_dp, -0.15_dp, KIND=dp)
      work(3, 3) = CMPLX(-0.03_dp, 0.0_dp, KIND=dp)
      work(4, 4) = CMPLX(0.08_dp, 0.0_dp, KIND=dp)
      work(3, 4) = CMPLX(-0.17_dp, 0.11_dp, KIND=dp)
      work(1, 3) = CMPLX(-0.18_dp, 0.09_dp, KIND=dp)
      work(1, 4) = CMPLX(0.05_dp, -0.13_dp, KIND=dp)
      work(2, 3) = CMPLX(0.07_dp, -0.05_dp, KIND=dp)
      work(2, 4) = CMPLX(-0.09_dp, 0.04_dp, KIND=dp)
      DO i = 2, n
         work(i, 1:i - 1) = CONJG(work(1:i - 1, i))
      END DO

      kernel(:, :) = CMPLX(9.0_dp, -3.0_dp, KIND=dp)
      CALL dense_fixed_n_projector_cluster_kernel(evals, work, occ_curv, occ_slope, gauge_gap, &
                                                  projector_weight, kernel)

      cluster_work_a(:, :) = 0.5_dp*(work(1:2, 1:2) + CONJG(TRANSPOSE(work(1:2, 1:2))))
      cluster_work_b(:, :) = 0.5_dp*(work(3:4, 3:4) + CONJG(TRANSPOSE(work(3:4, 3:4))))
      CALL diag_complex(cluster_work_a, cluster_evec_a, cluster_eval_a)
      CALL diag_complex(cluster_work_b, cluster_evec_b, cluster_eval_b)
      CALL dense_fixed_n_projector_modal_curvature(cluster_evec_a, occ_curv(1:2), modal_curv_a)
      CALL dense_fixed_n_projector_modal_curvature(cluster_evec_b, occ_curv(3:4), modal_curv_b)
      modal_energy_a(:) = 0.0_dp
      modal_energy_b(:) = 0.0_dp
      DO a = 1, 2
         DO b = 1, 2
            component_weight = ABS(cluster_evec_a(b, a))**2
            modal_energy_a(a) = modal_energy_a(a) + component_weight*evals(b)
            component_weight = ABS(cluster_evec_b(b, a))**2
            modal_energy_b(a) = modal_energy_b(a) + component_weight*evals(2 + b)
         END DO
      END DO
      CALL dense_projector_cluster_modal_occupation(evals(1:2), occ_curv(1:2), occ_slope(1:2), &
                                                    modal_energy_a, modal_occ_a, modal_curv_a)
      CALL dense_projector_cluster_modal_occupation(evals(3:4), occ_curv(3:4), occ_slope(3:4), &
                                                    modal_energy_b, modal_occ_b, modal_curv_b)
      modal_block(:, :) = MATMUL(CONJG(TRANSPOSE(cluster_evec_a)), MATMUL(work(1:2, 3:4), cluster_evec_b))
      DO b = 1, 2
         DO a = 1, 2
            denom = modal_energy_a(a) - modal_energy_b(b)
            IF (ABS(denom) > SQRT(EPSILON(1.0_dp))*MAX(1.0_dp, MAXVAL(ABS(evals(:))))) THEN
               coeff = (modal_occ_a(a) - modal_occ_b(b))/denom
            ELSE
               coeff = 0.0_dp
            END IF
            coeff = MIN(0.0_dp, coeff)
            slope_cap = dense_projector_pair_slope_cap(denom, gauge_gap, modal_curv_a(a), modal_curv_b(b), &
                                                        MAXVAL(MAX(0.0_dp, occ_curv(:))))
            coeff = MAX(coeff, -slope_cap)
            modal_block(a, b) = coeff*modal_block(a, b)
         END DO
      END DO
      expected_block(:, :) = projector_weight*MATMUL(cluster_evec_a, &
                                                     MATMUL(modal_block, CONJG(TRANSPOSE(cluster_evec_b))))

      scalar_energy_a = SUM(evals(1:2))/2.0_dp
      scalar_energy_b = SUM(evals(3:4))/2.0_dp
      scalar_occ_a = SUM(occ_slope(1:2))/2.0_dp
      scalar_occ_b = SUM(occ_slope(3:4))/2.0_dp
      scalar_curv_a = SUM(MAX(0.0_dp, occ_curv(1:2)))/2.0_dp
      scalar_curv_b = SUM(MAX(0.0_dp, occ_curv(3:4)))/2.0_dp
      denom = scalar_energy_a - scalar_energy_b
      IF (ABS(denom) > SQRT(EPSILON(1.0_dp))*MAX(1.0_dp, MAXVAL(ABS(evals(:))))) THEN
         coeff = (scalar_occ_a - scalar_occ_b)/denom
      ELSE
         coeff = 0.0_dp
      END IF
      coeff = MIN(0.0_dp, coeff)
      slope_cap = 2.0_dp*MAX(0.0_dp, MAX(scalar_curv_a, scalar_curv_b))
      coeff = MAX(coeff, -slope_cap)
      scalar_block(:, :) = projector_weight*coeff*work(1:2, 3:4)

      err = MAXVAL(ABS(kernel(1:2, 3:4) - expected_block(:, :)))
      herm_err = MAXVAL(ABS(kernel(3:4, 1:2) - CONJG(TRANSPOSE(kernel(1:2, 3:4)))))
      old_gap = MAXVAL(ABS(kernel(1:2, 3:4) - CMPLX(9.0_dp, -3.0_dp, KIND=dp)))
      scalar_gap = MAXVAL(ABS(expected_block(:, :) - scalar_block(:, :)))

      IF (err > tol .OR. herm_err > tol .OR. old_gap < 1.0_dp .OR. scalar_gap < 1.0E-5_dp) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "failed with error/hermitian/oldgap/scalargap", err, herm_err, old_gap, scalar_gap
      ELSE
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "OK with error/hermitian/oldgap/scalargap", err, herm_err, old_gap, scalar_gap
      END IF

   END SUBROUTINE check_history_safe_occ_projector_cluster_pair_response

! **************************************************************************************************
!> \brief Check that cluster modal occupations are physical fixed-N shapes, not basis projections.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_history_safe_occ_projector_cluster_modal_occ_shape(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 4
      REAL(KIND=dp), PARAMETER                           :: active_tol = 1.0E-5_dp, &
                                                            gauge_gap = 5.0E-5_dp, &
                                                            projector_weight = 0.37_dp

      COMPLEX(KIND=dp), DIMENSION(3)                     :: expected_new, expected_old, modal_col, &
                                                            response_col_new, response_col_old
      COMPLEX(KIND=dp), DIMENSION(3, 3)                  :: cluster_evec, cluster_work
      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: kernel, work
      INTEGER                                            :: a, b, i
      REAL(KIND=dp)                                      :: charge_err, coeff_new, coeff_old, &
                                                            component_weight, curv_sum, denom, err, &
                                                            herm_err, mu_response, raw_gap, shape_err, slope_cap
      REAL(KIND=dp), DIMENSION(3)                        :: cluster_eval, modal_curv, modal_energy, &
                                                            modal_occ_expected, modal_occ_new, modal_occ_old
      REAL(KIND=dp), DIMENSION(n)                        :: evals, occ_curv, occ_slope

      evals(:) = [-0.010000_dp, -0.009990_dp, -0.009980_dp, 1.23_dp]
      occ_slope(:) = [1.18_dp, 0.43_dp, 0.92_dp, 0.08_dp]
      occ_curv(:) = [0.42_dp, 0.36_dp, 0.31_dp, 0.06_dp]

      work(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      work(1, 1) = CMPLX(0.06_dp, 0.0_dp, KIND=dp)
      work(2, 2) = CMPLX(-0.04_dp, 0.0_dp, KIND=dp)
      work(3, 3) = CMPLX(0.025_dp, 0.0_dp, KIND=dp)
      work(1, 2) = CMPLX(0.22_dp, -0.15_dp, KIND=dp)
      work(1, 3) = CMPLX(-0.13_dp, 0.11_dp, KIND=dp)
      work(2, 3) = CMPLX(0.16_dp, 0.07_dp, KIND=dp)
      work(1, 4) = CMPLX(-0.18_dp, 0.09_dp, KIND=dp)
      work(2, 4) = CMPLX(0.07_dp, -0.05_dp, KIND=dp)
      work(3, 4) = CMPLX(0.11_dp, 0.04_dp, KIND=dp)
      DO i = 2, n
         work(i, 1:i - 1) = CONJG(work(1:i - 1, i))
      END DO

      kernel(:, :) = CMPLX(9.0_dp, -3.0_dp, KIND=dp)
      CALL dense_fixed_n_projector_cluster_kernel(evals, work, occ_curv, occ_slope, gauge_gap, &
                                                  projector_weight, kernel)

      cluster_work(:, :) = 0.5_dp*(work(1:3, 1:3) + CONJG(TRANSPOSE(work(1:3, 1:3))))
      CALL diag_complex(cluster_work, cluster_evec, cluster_eval)
      CALL dense_fixed_n_projector_modal_curvature(cluster_evec, occ_curv(1:3), modal_curv)
      modal_energy(:) = 0.0_dp
      modal_occ_old(:) = 0.0_dp
      DO a = 1, 3
         DO b = 1, 3
            component_weight = ABS(cluster_evec(b, a))**2
            modal_energy(a) = modal_energy(a) + component_weight*evals(b)
            modal_occ_old(a) = modal_occ_old(a) + component_weight*occ_slope(b)
         END DO
      END DO
      CALL dense_projector_cluster_modal_occupation(evals(1:3), occ_curv(1:3), occ_slope(1:3), &
                                                    modal_energy, modal_occ_new, modal_curv)
      curv_sum = SUM(MAX(0.0_dp, modal_curv(:)))
      mu_response = DOT_PRODUCT(MAX(0.0_dp, modal_curv(:)), modal_energy(:))/curv_sum
      DO a = 1, 3
         modal_occ_expected(a) = SUM(occ_slope(1:3))/3.0_dp - &
                                 2.0_dp*MAX(0.0_dp, modal_curv(a))* &
                                 (modal_energy(a) - mu_response)
      END DO

      modal_col(:) = MATMUL(CONJG(TRANSPOSE(cluster_evec)), work(1:3, 4))
      response_col_new(:) = modal_col(:)
      response_col_old(:) = modal_col(:)
      DO a = 1, 3
         denom = modal_energy(a) - evals(4)
         coeff_new = 0.0_dp
         coeff_old = 0.0_dp
         IF (ABS(denom) > SQRT(EPSILON(1.0_dp))*MAX(1.0_dp, MAXVAL(ABS(evals(:))))) THEN
            coeff_new = (modal_occ_new(a) - occ_slope(4))/denom
            coeff_old = (modal_occ_old(a) - occ_slope(4))/denom
         END IF
         coeff_new = MIN(0.0_dp, coeff_new)
         coeff_old = MIN(0.0_dp, coeff_old)
         slope_cap = dense_projector_pair_slope_cap(denom, gauge_gap, modal_curv(a), occ_curv(4), &
                                                     MAXVAL(MAX(0.0_dp, occ_curv(:))))
         coeff_new = MAX(coeff_new, -slope_cap)
         coeff_old = MAX(coeff_old, -slope_cap)
         response_col_new(a) = coeff_new*response_col_new(a)
         response_col_old(a) = coeff_old*response_col_old(a)
      END DO
      expected_new(:) = projector_weight*MATMUL(cluster_evec, response_col_new)
      expected_old(:) = projector_weight*MATMUL(cluster_evec, response_col_old)

      err = MAXVAL(ABS(kernel(1:3, 4) - expected_new(:)))
      herm_err = MAXVAL(ABS(kernel(4, 1:3) - CONJG(kernel(1:3, 4))))
      raw_gap = MAXVAL(ABS(expected_new(:) - expected_old(:)))
      charge_err = ABS(SUM(modal_occ_new(:)) - SUM(occ_slope(1:3)))
      shape_err = MAXVAL(ABS(modal_occ_new(:) - modal_occ_expected(:)))

      IF (err > tol .OR. herm_err > tol .OR. charge_err > tol .OR. shape_err > tol .OR. &
          raw_gap < active_tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,5(1X,ES12.5))') TRIM(label), &
            "failed with err/herm/charge/shape/rawgap", err, herm_err, charge_err, shape_err, raw_gap
      ELSE
         WRITE (*, '(A,1X,A,5(1X,ES12.5))') TRIM(label), &
            "OK with err/herm/charge/shape/rawgap", err, herm_err, charge_err, shape_err, raw_gap
      END IF

   END SUBROUTINE check_history_safe_occ_projector_cluster_modal_occ_shape

! **************************************************************************************************
!> \brief Check mixed cluster/singleton fixed-N kernel invariants independently of the formula mirror.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_history_safe_occ_projector_cluster_kernel_invariants(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 5
      REAL(KIND=dp), PARAMETER :: active_tol = 1.0E-5_dp, gauge_gap = 5.0E-5_dp, &
         projector_weight = 0.37_dp, projector_window_gap = 8.0E-3_dp

      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: kernel, kernel_back, kernel_q, q, work, &
                                                            work_q
      COMPLEX(KIND=dp)                                   :: phase1, phase2, sentinel
      INTEGER                                            :: i
      REAL(KIND=dp)                                      :: active, c1, c2, cov_err, herm_err, s1, s2, &
                                                            trace_leak
      REAL(KIND=dp), DIMENSION(n)                        :: evals, occ_curv, occ_slope

      evals(:) = [-0.010000_dp, -0.010000_dp, 0.105000_dp, 0.105000_dp, 0.42_dp]
      occ_slope(:) = [0.86_dp, 0.86_dp, 0.30_dp, 0.30_dp, 0.04_dp]
      occ_curv(:) = [0.052_dp, 0.052_dp, 0.027_dp, 0.027_dp, 0.006_dp]

      work(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      work(1, 1) = CMPLX(0.08_dp, 0.0_dp, KIND=dp)
      work(2, 2) = CMPLX(-0.03_dp, 0.0_dp, KIND=dp)
      work(3, 3) = CMPLX(0.04_dp, 0.0_dp, KIND=dp)
      work(4, 4) = CMPLX(-0.02_dp, 0.0_dp, KIND=dp)
      work(5, 5) = CMPLX(0.015_dp, 0.0_dp, KIND=dp)
      work(1, 2) = CMPLX(0.14_dp, -0.09_dp, KIND=dp)
      work(3, 4) = CMPLX(-0.11_dp, 0.07_dp, KIND=dp)
      work(1, 3) = CMPLX(0.05_dp, -0.03_dp, KIND=dp)
      work(1, 4) = CMPLX(-0.04_dp, -0.02_dp, KIND=dp)
      work(2, 3) = CMPLX(0.03_dp, 0.06_dp, KIND=dp)
      work(2, 4) = CMPLX(-0.02_dp, 0.05_dp, KIND=dp)
      work(1, 5) = CMPLX(0.04_dp, 0.01_dp, KIND=dp)
      work(3, 5) = CMPLX(-0.03_dp, 0.02_dp, KIND=dp)
      DO i = 2, n
         work(i, 1:i - 1) = CONJG(work(1:i - 1, i))
      END DO

      sentinel = CMPLX(9.0_dp, -3.0_dp, KIND=dp)
      kernel(:, :) = sentinel
      CALL dense_fixed_n_projector_cluster_kernel(evals, work, occ_curv, occ_slope, gauge_gap, &
                                                  projector_weight, kernel, projector_window_gap)

      trace_leak = 0.0_dp
      DO i = 1, n
         trace_leak = trace_leak + REAL(kernel(i, i), KIND=dp)
      END DO
      trace_leak = ABS(trace_leak)
      herm_err = MAXVAL(ABS(kernel(:, :) - CONJG(TRANSPOSE(kernel(:, :)))))
      active = MAXVAL(ABS(kernel(:, :)))

      q(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      c1 = COS(0.37_dp)
      s1 = SIN(0.37_dp)
      phase1 = CMPLX(COS(0.43_dp), SIN(0.43_dp), KIND=dp)
      q(1, 1) = CMPLX(c1, 0.0_dp, KIND=dp)
      q(1, 2) = phase1*s1
      q(2, 1) = -CONJG(phase1)*s1
      q(2, 2) = CMPLX(c1, 0.0_dp, KIND=dp)
      c2 = COS(0.29_dp)
      s2 = SIN(0.29_dp)
      phase2 = CMPLX(COS(-0.31_dp), SIN(-0.31_dp), KIND=dp)
      q(3, 3) = CMPLX(c2, 0.0_dp, KIND=dp)
      q(3, 4) = phase2*s2
      q(4, 3) = -CONJG(phase2)*s2
      q(4, 4) = CMPLX(c2, 0.0_dp, KIND=dp)
      q(5, 5) = CMPLX(1.0_dp, 0.0_dp, KIND=dp)

      work_q(:, :) = MATMUL(CONJG(TRANSPOSE(q)), MATMUL(work, q))
      kernel_q(:, :) = sentinel
      CALL dense_fixed_n_projector_cluster_kernel(evals, work_q, occ_curv, occ_slope, gauge_gap, &
                                                  projector_weight, kernel_q, projector_window_gap)
      kernel_back(:, :) = MATMUL(q, MATMUL(kernel_q, CONJG(TRANSPOSE(q))))
      cov_err = MAXVAL(ABS(kernel_back(:, :) - kernel(:, :)))

      IF (trace_leak > tol .OR. herm_err > tol .OR. cov_err > tol .OR. active < active_tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "failed with trace/herm/cov/active", trace_leak, herm_err, cov_err, active
      ELSE
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), &
            "OK with trace/herm/cov/active", trace_leak, herm_err, cov_err, active
      END IF

   END SUBROUTINE check_history_safe_occ_projector_cluster_kernel_invariants

! **************************************************************************************************
!> \brief Dense PAVA projection to non-increasing occupation slopes.
!> \param occ_scale ...
!> \param occ_tol ...
!> \param occ_mono ...
! **************************************************************************************************
   SUBROUTINE dense_monotone_occupation_slopes(occ_scale, occ_tol, occ_mono)
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: occ_scale
      REAL(KIND=dp), INTENT(IN)                          :: occ_tol
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: occ_mono

      INTEGER                                            :: i, iblock, idx, n, nblock
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: block_end, block_start
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: block_level, block_weight

      n = SIZE(occ_scale)
      IF (n <= 0) RETURN

      ALLOCATE (block_end(n), block_level(n), block_start(n), block_weight(n))
      nblock = 0
      DO i = 1, n
         nblock = nblock + 1
         block_start(nblock) = i
         block_end(nblock) = i
         block_level(nblock) = occ_scale(i)
         block_weight(nblock) = 1.0_dp
         DO WHILE (nblock > 1 .AND. &
                   block_level(nblock - 1) + occ_tol < block_level(nblock))
            block_level(nblock - 1) = &
               (block_weight(nblock - 1)*block_level(nblock - 1) + &
                block_weight(nblock)*block_level(nblock))/ &
               (block_weight(nblock - 1) + block_weight(nblock))
            block_weight(nblock - 1) = block_weight(nblock - 1) + block_weight(nblock)
            block_end(nblock - 1) = block_end(nblock)
            nblock = nblock - 1
         END DO
      END DO

      DO iblock = 1, nblock
         DO idx = block_start(iblock), block_end(iblock)
            occ_mono(idx) = block_level(iblock)
         END DO
      END DO
      DEALLOCATE (block_end, block_level, block_start, block_weight)

   END SUBROUTINE dense_monotone_occupation_slopes

! **************************************************************************************************
!> \brief Dense value 1/2 sum_i f_i lambda_i[U^H H U].
!> \param x_re ...
!> \param x_im ...
!> \param base ...
!> \param occ_scale ...
!> \return ...
! **************************************************************************************************
   FUNCTION dense_spectral_projector_value(x_re, x_im, base, occ_scale) RESULT(value)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: x_re, x_im
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: base
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: occ_scale
      REAL(KIND=dp)                                      :: value

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: chc, evec, u
      INTEGER                                            :: n
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: evals

      n = SIZE(x_re, 1)
      ALLOCATE (chc(n, n), evec(n, n), u(n, n), evals(n))
      CALL dense_rotation_u(x_re, x_im, u)
      chc(:, :) = MATMUL(CONJG(TRANSPOSE(u)), MATMUL(base, u))
      CALL diag_complex(chc, evec, evals)
      value = 0.5_dp*SUM(occ_scale(:)*evals(:))
      DEALLOCATE (chc, evec, u, evals)

   END FUNCTION dense_spectral_projector_value

! **************************************************************************************************
!> \brief Dense spectral occupation projector in the current rotated REF basis.
!> \param x_re ...
!> \param x_im ...
!> \param base ...
!> \param occ_scale ...
!> \param occ_mat ...
! **************************************************************************************************
   SUBROUTINE dense_spectral_occ_projector(x_re, x_im, base, occ_scale, occ_mat)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: x_re, x_im
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: base
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: occ_scale
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(OUT)     :: occ_mat

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: chc, u
      INTEGER                                            :: n

      n = SIZE(x_re, 1)
      ALLOCATE (chc(n, n), u(n, n))
      CALL dense_rotation_u(x_re, x_im, u)
      chc(:, :) = MATMUL(CONJG(TRANSPOSE(u)), MATMUL(base, u))
      CALL dense_spectral_occ_projector_from_chc(chc, occ_scale, occ_mat)
      DEALLOCATE (chc, u)

   END SUBROUTINE dense_spectral_occ_projector

! **************************************************************************************************
!> \brief Dense spectral occupation projector for a supplied Hermitian CHC matrix.
!> \param chc ...
!> \param occ_scale ...
!> \param occ_mat ...
! **************************************************************************************************
   SUBROUTINE dense_spectral_occ_projector_from_chc(chc, occ_scale, occ_mat)
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: chc
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: occ_scale
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(OUT)     :: occ_mat

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: evec, weighted_evec
      INTEGER                                            :: i, n
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: evals

      n = SIZE(chc, 1)
      ALLOCATE (evec(n, n), weighted_evec(n, n), evals(n))
      CALL diag_complex(chc, evec, evals)
      weighted_evec(:, :) = evec(:, :)
      DO i = 1, n
         weighted_evec(:, i) = occ_scale(i)*weighted_evec(:, i)
      END DO
      occ_mat(:, :) = MATMUL(weighted_evec, CONJG(TRANSPOSE(evec)))
      DEALLOCATE (evec, weighted_evec, evals)

   END SUBROUTINE dense_spectral_occ_projector_from_chc

! **************************************************************************************************
!> \brief Production-like history-safe occupation projector for one dense CHC fixture.
!> \param chc ...
!> \param occ_scale ...
!> \param projector_weight ...
!> \param occ_mat ...
!> \param occ_curv ...
!> \param gauge_gap ...
!> \param projector_window_gap ...
! **************************************************************************************************
   SUBROUTINE dense_history_safe_occ_projector_from_chc(chc, occ_scale, projector_weight, occ_mat, &
                                                        occ_curv, gauge_gap, projector_window_gap)
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: chc
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: occ_scale
      REAL(KIND=dp), INTENT(IN)                          :: projector_weight
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(OUT)     :: occ_mat
      REAL(KIND=dp), DIMENSION(:), INTENT(IN), OPTIONAL  :: occ_curv
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: gauge_gap
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: projector_window_gap

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: evec, weighted_evec
      INTEGER                                            :: i, n
      REAL(KIND=dp)                                      :: occ_tol
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: evals, occ_work

      n = SIZE(chc, 1)
      ALLOCATE (evec(n, n), weighted_evec(n, n), evals(n), occ_work(n))

      CALL diag_complex(chc, evec, evals)
      occ_work(:) = occ_scale(:)
      IF (PRESENT(occ_curv) .AND. PRESENT(gauge_gap)) THEN
         occ_tol = SQRT(EPSILON(1.0_dp))*MAX(1.0_dp, MAXVAL(ABS(occ_scale(:))))
         CALL dense_monotone_occupation_slopes(occ_scale, occ_tol, occ_work)
         CALL dense_near_degenerate_occupation_plateau(evals, occ_curv, gauge_gap, occ_tol, occ_work, &
                                                       projector_window_gap)
      END IF
      weighted_evec(:, :) = evec(:, :)
      DO i = 1, n
         weighted_evec(:, i) = occ_work(i)*weighted_evec(:, i)
      END DO
      occ_mat(:, :) = projector_weight*MATMUL(weighted_evec, CONJG(TRANSPOSE(evec)))
      DO i = 1, n
         occ_mat(i, i) = CMPLX(REAL(occ_mat(i, i), KIND=dp) + &
                               (1.0_dp - projector_weight)*occ_scale(i), 0.0_dp, KIND=dp)
      END DO

      DEALLOCATE (evec, weighted_evec, evals, occ_work)

   END SUBROUTINE dense_history_safe_occ_projector_from_chc

! **************************************************************************************************
!> \brief Shared pair slope for near-degenerate occupation-gauge filtering.
!> \param occ_curv_i ...
!> \param occ_curv_j ...
!> \return ...
! **************************************************************************************************
   FUNCTION dense_shared_occ_slope(occ_curv_i, occ_curv_j) RESULT(slope)
      REAL(KIND=dp), INTENT(IN)                          :: occ_curv_i, occ_curv_j
      REAL(KIND=dp)                                      :: slope

      slope = 2.0_dp*SQRT(MAX(0.0_dp, occ_curv_i)*MAX(0.0_dp, occ_curv_j))

   END FUNCTION dense_shared_occ_slope

! **************************************************************************************************
!> \brief Production-like local gap for physical occupation-response projector clusters.
!> \param gauge_gap ...
!> \param occ_curv_i ...
!> \param occ_curv_j ...
!> \param occ_curv_max ...
!> \param projector_window_gap ...
!> \return thermally widened pair gap for active Fermi-response modes
! **************************************************************************************************
   FUNCTION dense_projector_occupation_pair_gap(gauge_gap, occ_curv_i, occ_curv_j, &
                                                occ_curv_max, projector_window_gap) RESULT(pair_gap)
      REAL(KIND=dp), INTENT(IN)                          :: gauge_gap, occ_curv_i, occ_curv_j, &
                                                            occ_curv_max
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: projector_window_gap
      REAL(KIND=dp)                                      :: pair_gap

      REAL(KIND=dp), PARAMETER :: thermal_gap_fraction = 2.5E-1_dp, &
         thermal_gap_max_factor = 5.0E1_dp

      REAL(KIND=dp)                                      :: active_fraction, max_gap, pair_slope, &
                                                            physical_gap, thermal_gap

      pair_gap = gauge_gap
      IF (gauge_gap <= EPSILON(gauge_gap)) RETURN
      IF (occ_curv_max <= EPSILON(occ_curv_max)) RETURN

      physical_gap = gauge_gap
      IF (PRESENT(projector_window_gap)) physical_gap = MAX(physical_gap, projector_window_gap)
      max_gap = MAX(thermal_gap_max_factor*gauge_gap, physical_gap)
      thermal_gap = thermal_gap_fraction/occ_curv_max
      thermal_gap = MAX(thermal_gap, physical_gap)
      thermal_gap = MIN(max_gap, thermal_gap)

      pair_slope = dense_shared_occ_slope(occ_curv_i, occ_curv_j)
      active_fraction = pair_slope/occ_curv_max
      active_fraction = MIN(1.0_dp, MAX(0.0_dp, active_fraction))
      thermal_gap = MAX(gauge_gap, thermal_gap)
      pair_gap = gauge_gap + active_fraction*(thermal_gap - gauge_gap)

   END FUNCTION dense_projector_occupation_pair_gap

! **************************************************************************************************
!> \brief Production-like divided-difference slope cap for finite occupation-projector response pairs.
!> \param denom ...
!> \param gauge_gap ...
!> \param occ_curv_i ...
!> \param occ_curv_j ...
!> \param occ_curv_max ...
!> \param projector_window_gap ...
!> \return slope cap with shared-slope filtering inside the thermal response window
! **************************************************************************************************
   FUNCTION dense_projector_pair_slope_cap(denom, gauge_gap, occ_curv_i, occ_curv_j, occ_curv_max, &
                                           projector_window_gap) RESULT(slope_cap)
      REAL(KIND=dp), INTENT(IN)                          :: denom, gauge_gap, occ_curv_i, occ_curv_j, &
                                                            occ_curv_max
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: projector_window_gap
      REAL(KIND=dp)                                      :: slope_cap

      REAL(KIND=dp)                                      :: local_gap

      slope_cap = 2.0_dp*MAX(0.0_dp, MAX(occ_curv_i, occ_curv_j))
      IF (slope_cap <= EPSILON(slope_cap)) RETURN

      IF (PRESENT(projector_window_gap)) THEN
         local_gap = dense_projector_occupation_pair_gap(gauge_gap, occ_curv_i, occ_curv_j, &
                                                         occ_curv_max, projector_window_gap)
      ELSE
         local_gap = dense_projector_occupation_pair_gap(gauge_gap, occ_curv_i, occ_curv_j, &
                                                         occ_curv_max)
      END IF
      IF (ABS(denom) <= local_gap) slope_cap = MIN(slope_cap, dense_shared_occ_slope(occ_curv_i, occ_curv_j))

   END FUNCTION dense_projector_pair_slope_cap

! **************************************************************************************************
!> \brief occupation-window projector weight for one dense rotation-pair fixture
!> \param projector_weight ...
!> \param curv_i ...
!> \param curv_j ...
!> \param curv_max ...
!> \return ...
! **************************************************************************************************
   FUNCTION dense_pair_projector_weight(projector_weight, curv_i, curv_j, curv_max) RESULT(pair_weight)
      REAL(KIND=dp), INTENT(IN)                          :: projector_weight, curv_i, curv_j, &
                                                            curv_max
      REAL(KIND=dp)                                      :: pair_weight

      REAL(KIND=dp)                                      :: active_fraction

      pair_weight = 0.0_dp
      IF (projector_weight <= 0.0_dp) RETURN
      IF (curv_max <= EPSILON(curv_max)) RETURN

      active_fraction = dense_shared_occ_slope(curv_i, curv_j)/(2.0_dp*curv_max)
      active_fraction = MIN(1.0_dp, MAX(0.0_dp, active_fraction))
      pair_weight = projector_weight*active_fraction

   END FUNCTION dense_pair_projector_weight

! **************************************************************************************************
!> \brief history-qualified physical projector weight for one dense rotation variable
!> \param projector_weight ...
!> \param active_weight ...
!> \param occ_gap_weight ...
!> \param density_weight ...
!> \param density_response_weight optional fixed-N density-response character
!> \return ...
! **************************************************************************************************
   FUNCTION dense_physical_pair_projector_weight(projector_weight, active_weight, occ_gap_weight, &
                                                 density_weight, density_response_weight) RESULT(pair_weight)
      REAL(KIND=dp), INTENT(IN)                          :: projector_weight, active_weight, &
                                                            occ_gap_weight, density_weight
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: density_response_weight
      REAL(KIND=dp)                                      :: pair_weight

      REAL(KIND=dp)                                      :: gap_fraction, qualified_response_weight, &
                                                            response_character, response_weight

      pair_weight = MAX(0.0_dp, active_weight)
      IF (projector_weight <= 0.0_dp) RETURN

      gap_fraction = MIN(1.0_dp, MAX(0.0_dp, occ_gap_weight))
      response_character = gap_fraction
      IF (PRESENT(density_response_weight)) &
         response_character = MAX(response_character, &
                                  MIN(1.0_dp, MAX(0.0_dp, density_response_weight)))
      response_weight = MIN(1.0_dp, MAX(0.0_dp, density_weight))
      qualified_response_weight = SQRT(SQRT(response_weight))
      pair_weight = MAX(pair_weight, projector_weight*response_character*qualified_response_weight)
      pair_weight = MIN(projector_weight, MAX(0.0_dp, pair_weight))

   END FUNCTION dense_physical_pair_projector_weight

! **************************************************************************************************
!> \brief active Fermi-window weight for dense projector-density response fixtures.
!> \param occ_curv ...
!> \param occ_curv_max ...
!> \return ...
! **************************************************************************************************
   FUNCTION dense_active_density_response_weight(occ_curv, occ_curv_max) RESULT(response_weight)
      REAL(KIND=dp), INTENT(IN)                          :: occ_curv, occ_curv_max
      REAL(KIND=dp)                                      :: response_weight

      response_weight = 0.0_dp
      IF (occ_curv <= EPSILON(occ_curv)) RETURN
      IF (occ_curv_max <= EPSILON(occ_curv_max)) RETURN

      response_weight = SQRT(MIN(1.0_dp, MAX(0.0_dp, occ_curv/occ_curv_max)))

   END FUNCTION dense_active_density_response_weight

! **************************************************************************************************
!> \brief bounded diagonal projector-density response matching the production Schur helper.
!> \param density_projector_response ...
!> \param occ_curv ...
!> \param occ_curv_max ...
!> \param response_ratio_gate ...
!> \return ...
! **************************************************************************************************
   FUNCTION dense_bounded_density_response(density_projector_response, occ_curv, occ_curv_max, &
                                           response_ratio_gate) RESULT(density_response)
      REAL(KIND=dp), INTENT(IN)                          :: density_projector_response, occ_curv, &
                                                            occ_curv_max, response_ratio_gate
      REAL(KIND=dp)                                      :: density_response

      REAL(KIND=dp)                                      :: active_fraction

      density_response = 0.0_dp
      IF (occ_curv <= EPSILON(occ_curv)) RETURN
      IF (occ_curv_max <= EPSILON(occ_curv_max)) RETURN

      active_fraction = occ_curv/MAX(occ_curv, response_ratio_gate*occ_curv_max)
      active_fraction = MIN(1.0_dp, MAX(0.0_dp, active_fraction))
      density_response = dense_active_density_response_weight(occ_curv, occ_curv_max)* &
                         density_projector_response* &
                         (2.0_dp*active_fraction - active_fraction*active_fraction)

   END FUNCTION dense_bounded_density_response

! **************************************************************************************************
!> \brief gauge-free physical density response matching the production Schur helper.
!> \param density_projector_response ...
!> \param density_safe_response audit-only dP/curvature response ...
!> \param occ_curv ...
!> \param occ_curv_max ...
!> \param response_ratio_gate ...
!> \param density_projector_centered_response ...
!> \param density_safe_centered_response ...
!> \param density_response_alignment ...
!> \return ...
! **************************************************************************************************
   FUNCTION dense_physical_density_response(density_projector_response, density_safe_response, occ_curv, &
                                            occ_curv_max, response_ratio_gate, &
                                            density_projector_centered_response, &
                                            density_safe_centered_response, &
                                            density_response_alignment) RESULT(density_response)
      REAL(KIND=dp), INTENT(IN)                          :: density_projector_response, &
                                                            density_safe_response, occ_curv, &
                                                            occ_curv_max, response_ratio_gate
      REAL(KIND=dp), INTENT(IN), OPTIONAL :: density_projector_centered_response, &
         density_safe_centered_response, density_response_alignment
      REAL(KIND=dp)                                      :: density_response

      LOGICAL                                            :: anti_aligned_safe_response, &
                                                            globally_aligned_safe_response, &
                                                            locally_aligned_safe_response
      REAL(KIND=dp) :: active_fraction, alignment, projector_abs, projector_component, response_shape, &
         safe_abs, safe_blend, safe_cap, safe_component, sign_projector_response, sign_safe_response

      density_response = 0.0_dp
      IF (occ_curv <= EPSILON(occ_curv)) RETURN
      IF (occ_curv_max <= EPSILON(occ_curv_max)) RETURN

      active_fraction = occ_curv/MAX(occ_curv, response_ratio_gate*occ_curv_max)
      active_fraction = MIN(1.0_dp, MAX(0.0_dp, active_fraction))
      response_shape = 2.0_dp*active_fraction - active_fraction*active_fraction
      projector_component = density_projector_response
      sign_projector_response = density_projector_response
      sign_safe_response = density_safe_response
      IF (PRESENT(density_projector_centered_response)) THEN
         projector_component = density_projector_centered_response
         sign_projector_response = density_projector_centered_response
      END IF
      IF (PRESENT(density_safe_centered_response)) &
         sign_safe_response = density_safe_centered_response
      safe_component = projector_component
      projector_abs = ABS(projector_component)
      safe_abs = ABS(sign_safe_response)
      alignment = 0.0_dp
      IF (PRESENT(density_response_alignment)) &
         alignment = MIN(1.0_dp, MAX(-1.0_dp, density_response_alignment))
      locally_aligned_safe_response = sign_projector_response*sign_safe_response > 0.0_dp .OR. &
                                      (ABS(sign_projector_response) <= EPSILON(sign_projector_response) .AND. &
                                       ABS(sign_safe_response) > EPSILON(sign_safe_response))
      globally_aligned_safe_response = alignment > 0.25_dp .AND. &
                                       safe_abs > EPSILON(safe_abs) .AND. &
                                       projector_abs <= alignment*safe_abs
      anti_aligned_safe_response = PRESENT(density_projector_centered_response) .AND. &
                                   PRESENT(density_safe_centered_response) .AND. &
                                   sign_projector_response*sign_safe_response < 0.0_dp .AND. &
                                   safe_abs > EPSILON(safe_abs) .AND. &
                                   projector_abs > EPSILON(projector_abs)
      IF (locally_aligned_safe_response .OR. globally_aligned_safe_response) THEN
         IF (globally_aligned_safe_response) THEN
            safe_cap = MAX(4.0_dp*projector_abs, SQRT(response_shape)*safe_abs)
         ELSE
            safe_cap = MAX(4.0_dp*projector_abs, response_shape*safe_abs)
         END IF
         safe_component = SIGN(MIN(safe_abs, safe_cap), sign_safe_response)
         safe_blend = active_fraction
         safe_component = (1.0_dp - safe_blend)*projector_component + &
                          safe_blend*safe_component
      ELSE IF (anti_aligned_safe_response) THEN
         safe_component = SIGN(MIN(projector_abs, safe_abs), sign_safe_response)
      END IF
      density_response = response_shape*safe_component

   END FUNCTION dense_physical_density_response

! **************************************************************************************************
!> \brief shadow-only finite-floor density response matching the production Schur audit helper.
!> \param density_projector_response ...
!> \param density_safe_response audit-only dP/curvature response ...
!> \param occ_curv ...
!> \param occ_curv_max ...
!> \param response_ratio_gate ...
!> \return ...
! **************************************************************************************************
   FUNCTION dense_physical_density_floor_audit_response(density_projector_response, density_safe_response, &
                                                        occ_curv, occ_curv_max, &
                                                        response_ratio_gate) RESULT(density_response)
      REAL(KIND=dp), INTENT(IN)                          :: density_projector_response, &
                                                            density_safe_response, occ_curv, &
                                                            occ_curv_max, response_ratio_gate
      REAL(KIND=dp)                                      :: density_response

      REAL(KIND=dp)                                      :: active_fraction, response_shape

      density_response = 0.0_dp
      IF (occ_curv <= EPSILON(occ_curv)) RETURN
      IF (occ_curv_max <= EPSILON(occ_curv_max)) RETURN

      active_fraction = occ_curv/MAX(occ_curv, response_ratio_gate*occ_curv_max)
      active_fraction = MIN(1.0_dp, MAX(0.0_dp, active_fraction))
      response_shape = 2.0_dp*active_fraction - active_fraction*active_fraction
      density_response = response_shape*density_projector_response

   END FUNCTION dense_physical_density_floor_audit_response

! **************************************************************************************************
!> \brief Fréchet derivative of the dense spectral occupation projector.
!> \param chc ...
!> \param dchc ...
!> \param occ_scale ...
!> \param docc_mat ...
!> \param occ_curv ...
!> \param gauge_gap ...
!> \param projector_window_gap ...
! **************************************************************************************************
   SUBROUTINE dense_spectral_occ_projector_derivative_chc(chc, dchc, occ_scale, docc_mat, occ_curv, gauge_gap, &
                                                          projector_window_gap)
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: chc, dchc
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: occ_scale
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(OUT)     :: docc_mat
      REAL(KIND=dp), DIMENSION(:), INTENT(IN), OPTIONAL  :: occ_curv
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: gauge_gap
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: projector_window_gap

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: evec, kernel, work
      INTEGER                                            :: i, j, n
      REAL(KIND=dp) :: allowed_jump, coeff, denom, diag_response_shift, gap_tol, occ_curv_sum, &
         occ_jump, occ_tol, physical_slope, scale, slope_cap
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: evals, occ_work

      n = SIZE(chc, 1)
      ALLOCATE (evec(n, n), kernel(n, n), work(n, n), evals(n), occ_work(n))

      CALL diag_complex(chc, evec, evals)
      work(:, :) = MATMUL(CONJG(TRANSPOSE(evec)), MATMUL(dchc, evec))
      scale = MAX(1.0_dp, MAXVAL(ABS(evals(:))))
      gap_tol = SQRT(EPSILON(1.0_dp))*scale
      occ_tol = SQRT(EPSILON(1.0_dp))*MAX(1.0_dp, MAXVAL(ABS(occ_scale(:))))
      occ_work(:) = occ_scale(:)
      diag_response_shift = 0.0_dp
      IF (PRESENT(occ_curv)) THEN
         occ_curv_sum = SUM(MAX(0.0_dp, occ_curv(:)))
         IF (occ_curv_sum > EPSILON(occ_curv_sum)) THEN
            DO i = 1, n
               diag_response_shift = diag_response_shift + &
                                     MAX(0.0_dp, occ_curv(i))*REAL(work(i, i), KIND=dp)
            END DO
            diag_response_shift = diag_response_shift/occ_curv_sum
         END IF
      END IF
      kernel(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      DO j = 1, n
         DO i = 1, n
            IF (i == j) THEN
               IF (PRESENT(occ_curv)) THEN
                  coeff = -2.0_dp*MAX(0.0_dp, occ_curv(i))
                  IF (ABS(coeff) > EPSILON(coeff)) &
                     kernel(i, i) = CMPLX(coeff*(REAL(work(i, i), KIND=dp) - diag_response_shift), &
                                          0.0_dp, KIND=dp)
               END IF
               CYCLE
            END IF
            denom = evals(i) - evals(j)
            physical_slope = 0.0_dp
            IF (PRESENT(occ_curv) .AND. PRESENT(gauge_gap)) THEN
               IF (ABS(denom) <= gauge_gap) THEN
                  physical_slope = dense_shared_occ_slope(occ_curv(i), occ_curv(j))
                  allowed_jump = 4.0_dp*physical_slope*MAX(ABS(denom), EPSILON(1.0_dp)) + occ_tol
                  occ_jump = ABS(occ_work(i) - occ_work(j))
                  IF (occ_jump > allowed_jump) CYCLE
               END IF
            END IF
            IF (ABS(denom) > gap_tol) THEN
               coeff = (occ_work(i) - occ_work(j))/denom
            ELSE IF (ABS(occ_work(i) - occ_work(j)) <= occ_tol) THEN
               coeff = 0.0_dp
            ELSE
               coeff = (occ_work(i) - occ_work(j))/SIGN(gap_tol, denom)
            END IF
            coeff = MIN(0.0_dp, coeff)
            IF (PRESENT(occ_curv)) THEN
               IF (PRESENT(gauge_gap)) THEN
                  IF (PRESENT(projector_window_gap)) THEN
                     slope_cap = dense_projector_pair_slope_cap(denom, gauge_gap, occ_curv(i), occ_curv(j), &
                                                                 MAXVAL(MAX(0.0_dp, occ_curv(:))), &
                                                                 projector_window_gap)
                  ELSE
                     slope_cap = dense_projector_pair_slope_cap(denom, gauge_gap, occ_curv(i), occ_curv(j), &
                                                                 MAXVAL(MAX(0.0_dp, occ_curv(:))))
                  END IF
               ELSE
                  slope_cap = 2.0_dp*MAX(0.0_dp, MAX(occ_curv(i), occ_curv(j)))
               END IF
               IF (slope_cap > EPSILON(slope_cap)) THEN
                  coeff = MAX(coeff, -slope_cap)
               ELSE
                  coeff = 0.0_dp
               END IF
            END IF
            kernel(i, j) = coeff*work(i, j)
         END DO
      END DO

      docc_mat(:, :) = MATMUL(evec, MATMUL(kernel, CONJG(TRANSPOSE(evec))))
      docc_mat(:, :) = 0.5_dp*(docc_mat(:, :) + CONJG(TRANSPOSE(docc_mat(:, :))))
      DEALLOCATE (evec, kernel, work, evals, occ_work)

   END SUBROUTINE dense_spectral_occ_projector_derivative_chc

! **************************************************************************************************
!> \brief Production-like history-safe derivative: monotone slopes before projector dP.
!> \param chc ...
!> \param dchc ...
!> \param occ_scale ...
!> \param projector_weight ...
!> \param docc_mat ...
!> \param occ_curv ...
!> \param gauge_gap ...
!> \param fd_cluster_used ...
!> \param projector_window_gap ...
! **************************************************************************************************
   SUBROUTINE dense_history_safe_occ_projector_derivative_chc(chc, dchc, occ_scale, projector_weight, &
                                                              docc_mat, occ_curv, gauge_gap, fd_cluster_used, &
                                                              projector_window_gap)
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: chc, dchc
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: occ_scale
      REAL(KIND=dp), INTENT(IN)                          :: projector_weight
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(OUT)     :: docc_mat
      REAL(KIND=dp), DIMENSION(:), INTENT(IN), OPTIONAL  :: occ_curv
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: gauge_gap
      LOGICAL, INTENT(OUT), OPTIONAL                     :: fd_cluster_used
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: projector_window_gap

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: evec, kernel, occ_mat_m, occ_mat_p, work
      INTEGER                                            :: i, idx, j, n
      LOGICAL                                            :: same_gauge_cluster, use_projector_fd
      REAL(KIND=dp)                                      :: gap_tol, local_gap, occ_curv_max, &
                                                            occ_tol, projector_fd_step, scale, &
                                                            work_norm, work_tol
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: evals, occ_work

      docc_mat(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      IF (PRESENT(fd_cluster_used)) fd_cluster_used = .FALSE.
      IF (projector_weight <= 0.0_dp) RETURN

      ALLOCATE (occ_work(SIZE(occ_scale)))
      occ_tol = SQRT(EPSILON(1.0_dp))*MAX(1.0_dp, MAXVAL(ABS(occ_scale(:))))
      CALL dense_monotone_occupation_slopes(occ_scale, occ_tol, occ_work)
      IF (PRESENT(occ_curv) .AND. PRESENT(gauge_gap)) THEN
         n = SIZE(occ_scale)
         ALLOCATE (evec(n, n), evals(n), kernel(n, n), work(n, n))
         CALL diag_complex(chc, evec, evals)
         work(:, :) = MATMUL(CONJG(TRANSPOSE(evec)), MATMUL(dchc, evec))
         use_projector_fd = .FALSE.
         work_norm = MAXVAL(ABS(work(:, :)))
         work_tol = SQRT(EPSILON(1.0_dp))*MAX(1.0_dp, work_norm)
         scale = MAX(1.0_dp, MAXVAL(ABS(evals(:))))
         gap_tol = SQRT(EPSILON(1.0_dp))*scale
         occ_curv_max = MAXVAL(MAX(0.0_dp, occ_curv(:)))
         DO j = 1, n
            DO i = 1, n
               same_gauge_cluster = i /= j
               IF (same_gauge_cluster) THEN
                  DO idx = MIN(i, j), MAX(i, j) - 1
                     local_gap = dense_projector_occupation_pair_gap(gauge_gap, occ_curv(idx), &
                                                                     occ_curv(idx + 1), occ_curv_max, &
                                                                     projector_window_gap)
                     IF (ABS(evals(idx + 1) - evals(idx)) > local_gap) same_gauge_cluster = .FALSE.
                  END DO
               END IF
               IF (same_gauge_cluster .AND. &
                   (ABS(work(i, j)) > work_tol .OR. &
                    ABS(REAL(work(i, i), KIND=dp) - REAL(work(j, j), KIND=dp)) > work_tol)) &
                  use_projector_fd = .TRUE.
            END DO
         END DO
         CALL dense_near_degenerate_occupation_plateau(evals, occ_curv, gauge_gap, occ_tol, &
                                                       occ_work, projector_window_gap)
         IF (use_projector_fd) THEN
            IF (PRESENT(fd_cluster_used)) fd_cluster_used = .TRUE.
            projector_fd_step = 1.0E-7_dp
            IF (work_norm > EPSILON(work_norm)) THEN
               projector_fd_step = MIN(projector_fd_step, 1.0E-2_dp*gauge_gap/work_norm)
               projector_fd_step = MAX(projector_fd_step, &
                                       10.0_dp*EPSILON(1.0_dp)*MAX(1.0_dp, scale)/work_norm)
            END IF
            ALLOCATE (occ_mat_p(n, n), occ_mat_m(n, n))
            CALL dense_history_safe_occ_projector_from_chc(chc + projector_fd_step*dchc, occ_scale, &
                                                           projector_weight, occ_mat_p, occ_curv, gauge_gap, &
                                                           projector_window_gap)
            CALL dense_history_safe_occ_projector_from_chc(chc - projector_fd_step*dchc, occ_scale, &
                                                           projector_weight, occ_mat_m, occ_curv, gauge_gap, &
                                                           projector_window_gap)
            docc_mat(:, :) = (occ_mat_p(:, :) - occ_mat_m(:, :))/(2.0_dp*projector_fd_step)
            kernel(:, :) = MATMUL(CONJG(TRANSPOSE(evec)), MATMUL(docc_mat, evec))
            CALL dense_near_degenerate_occupation_plateau(evals, occ_curv, gauge_gap, occ_tol, occ_work, &
                                                          projector_window_gap)
            CALL dense_fixed_n_projector_cluster_kernel(evals, work, occ_curv, occ_work, gauge_gap, &
                                                        projector_weight, kernel, projector_window_gap)
            docc_mat(:, :) = MATMUL(evec, MATMUL(kernel, CONJG(TRANSPOSE(evec))))
            docc_mat(:, :) = 0.5_dp*(docc_mat(:, :) + CONJG(TRANSPOSE(docc_mat(:, :))))
            DEALLOCATE (evec, evals, kernel, work, occ_mat_p, occ_mat_m, occ_work)
            RETURN
         END IF
         CALL dense_near_degenerate_occupation_plateau(evals, occ_curv, gauge_gap, occ_tol, occ_work, &
                                                       projector_window_gap)
         DEALLOCATE (evec, evals, kernel, work)
      END IF
      IF (PRESENT(occ_curv) .AND. PRESENT(gauge_gap)) THEN
         CALL dense_spectral_occ_projector_derivative_chc(chc, dchc, occ_work, docc_mat, occ_curv, gauge_gap, &
                                                          projector_window_gap)
      ELSE
         CALL dense_spectral_occ_projector_derivative_chc(chc, dchc, occ_work, docc_mat)
      END IF
      docc_mat(:, :) = projector_weight*docc_mat(:, :)
      DEALLOCATE (occ_work)

   END SUBROUTINE dense_history_safe_occ_projector_derivative_chc

! **************************************************************************************************
!> \brief fixed-N physical occupations for dense finite Projector cluster modes.
!> \param evals ...
!> \param occ_curv ...
!> \param occ_slope ...
!> \param modal_energy ...
!> \param modal_occ ...
!> \param modal_curv ...
! **************************************************************************************************
   SUBROUTINE dense_projector_cluster_modal_occupation(evals, occ_curv, occ_slope, modal_energy, modal_occ, &
                                                       modal_curv)
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: evals, occ_curv, occ_slope, modal_energy
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: modal_occ
      REAL(KIND=dp), DIMENSION(:), INTENT(IN), OPTIONAL  :: modal_curv

      INTEGER                                            :: i, m
      REAL(KIND=dp)                                      :: curv, curv_sum, delta_e, eval_center, fit_den, &
                                                            fit_num, mu_shift, occ_mean, physical_slope, &
                                                            slope_cap

      m = SIZE(evals)
      modal_occ(:) = 0.0_dp
      IF (m <= 0) RETURN

      IF (PRESENT(modal_curv)) THEN
         curv_sum = 0.0_dp
         mu_shift = 0.0_dp
         DO i = 1, m
            curv = MAX(0.0_dp, modal_curv(i))
            curv_sum = curv_sum + curv
            mu_shift = mu_shift + curv*modal_energy(i)
         END DO
         IF (curv_sum > EPSILON(curv_sum)) THEN
            mu_shift = mu_shift/curv_sum
            occ_mean = SUM(occ_slope(1:m))/REAL(m, KIND=dp)
            DO i = 1, m
               modal_occ(i) = occ_mean - 2.0_dp*MAX(0.0_dp, modal_curv(i))* &
                              (modal_energy(i) - mu_shift)
            END DO
            RETURN
         END IF
      END IF

      occ_mean = SUM(occ_slope(1:m))/REAL(m, KIND=dp)
      eval_center = SUM(evals(1:m))/REAL(m, KIND=dp)
      fit_den = 0.0_dp
      fit_num = 0.0_dp
      DO i = 1, m
         delta_e = evals(i) - eval_center
         fit_den = fit_den + delta_e*delta_e
         fit_num = fit_num + delta_e*(occ_slope(i) - occ_mean)
      END DO
      physical_slope = 0.0_dp
      IF (fit_den > EPSILON(fit_den)) physical_slope = fit_num/fit_den
      physical_slope = MIN(0.0_dp, physical_slope)

      slope_cap = HUGE(1.0_dp)
      DO i = 1, m - 1
         slope_cap = MIN(slope_cap, dense_shared_occ_slope(occ_curv(i), occ_curv(i + 1)))
      END DO
      IF (slope_cap < HUGE(1.0_dp)/2.0_dp) THEN
         physical_slope = MAX(physical_slope, -MAX(0.0_dp, slope_cap))
      ELSE
         physical_slope = 0.0_dp
      END IF

      DO i = 1, SIZE(modal_occ)
         modal_occ(i) = occ_mean + physical_slope*(modal_energy(i) - eval_center)
      END DO

   END SUBROUTINE dense_projector_cluster_modal_occupation

! **************************************************************************************************
!> \brief Replace near-degenerate projector-gauge blocks by a fixed-N physical response.
!> \param evals ...
!> \param work ...
!> \param occ_curv ...
!> \param occ_slope ...
!> \param gauge_gap ...
!> \param projector_weight ...
!> \param kernel ...
!> \param projector_window_gap ...
! **************************************************************************************************
   SUBROUTINE dense_fixed_n_projector_cluster_kernel(evals, work, occ_curv, occ_slope, gauge_gap, &
                                                     projector_weight, kernel, projector_window_gap)
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: evals
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: work
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: occ_curv, occ_slope
      REAL(KIND=dp), INTENT(IN)                          :: gauge_gap, projector_weight
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(INOUT)   :: kernel
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: projector_window_gap

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:)        :: cluster_col, modal_col, response_col
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: cluster_evec, cluster_evec_b, &
                                                            cluster_kernel, cluster_work, &
                                                            cluster_work_b, modal_block, &
                                                            response_block
      INTEGER                                            :: a, b, first, first_b, j, last, last_b, &
                                                            m, m_b, n
      LOGICAL                                            :: outside_gauge_cluster
      REAL(KIND=dp)                                      :: component_weight, curv, denom, gap_tol, &
                                                            local_gap, mu_den, mu_num, mu_shift, &
                                                            occ_curv_max, offdiag_coeff, slope_cap
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: cluster_eval, cluster_eval_b, &
                                                            modal_curv, modal_curv_b, &
                                                            modal_energy, modal_energy_b, &
                                                            modal_occ, modal_occ_b

      n = SIZE(evals)
      IF (n <= 0) RETURN
      IF (projector_weight <= 0.0_dp) RETURN
      gap_tol = SQRT(EPSILON(1.0_dp))*MAX(1.0_dp, MAXVAL(ABS(evals(:))))
      occ_curv_max = MAXVAL(MAX(0.0_dp, occ_curv(:)))

      mu_num = 0.0_dp
      mu_den = 0.0_dp
      first = 1
      DO WHILE (first <= n)
         last = first
         DO WHILE (last < n)
            local_gap = dense_projector_occupation_pair_gap(gauge_gap, occ_curv(last), occ_curv(last + 1), &
                                                            occ_curv_max, projector_window_gap)
            IF (ABS(evals(last + 1) - evals(last)) > local_gap) EXIT
            last = last + 1
         END DO
         m = last - first + 1
         IF (m > 1) THEN
            ALLOCATE (cluster_work(m, m), cluster_evec(m, m), cluster_eval(m), modal_curv(m))
            cluster_work(:, :) = 0.5_dp*(work(first:last, first:last) + &
                                         CONJG(TRANSPOSE(work(first:last, first:last))))
            CALL diag_complex(cluster_work, cluster_evec, cluster_eval)
            CALL dense_fixed_n_projector_modal_curvature(cluster_evec, occ_curv(first:last), modal_curv)
            DO a = 1, m
               mu_num = mu_num + modal_curv(a)*cluster_eval(a)
               mu_den = mu_den + modal_curv(a)
            END DO
            DEALLOCATE (cluster_work, cluster_evec, cluster_eval, modal_curv)
         ELSE
            curv = MAX(0.0_dp, occ_curv(first))
            IF (curv > EPSILON(curv)) THEN
               mu_num = mu_num + curv*REAL(work(first, first), KIND=dp)
               mu_den = mu_den + curv
            END IF
         END IF
         first = last + 1
      END DO
      mu_shift = 0.0_dp
      IF (mu_den > EPSILON(mu_den)) mu_shift = mu_num/mu_den

      first = 1
      DO WHILE (first <= n)
         last = first
         DO WHILE (last < n)
            local_gap = dense_projector_occupation_pair_gap(gauge_gap, occ_curv(last), occ_curv(last + 1), &
                                                            occ_curv_max, projector_window_gap)
            IF (ABS(evals(last + 1) - evals(last)) > local_gap) EXIT
            last = last + 1
         END DO
         m = last - first + 1
         IF (m > 1) THEN
            kernel(first:last, first:last) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
            ALLOCATE (cluster_col(m), cluster_work(m, m), cluster_evec(m, m), cluster_kernel(m, m), &
                      modal_col(m), modal_energy(m), modal_occ(m), response_block(m, m), &
                      response_col(m), cluster_eval(m), modal_curv(m))
            cluster_work(:, :) = 0.5_dp*(work(first:last, first:last) + &
                                         CONJG(TRANSPOSE(work(first:last, first:last))))
            CALL diag_complex(cluster_work, cluster_evec, cluster_eval)
            CALL dense_fixed_n_projector_modal_curvature(cluster_evec, occ_curv(first:last), modal_curv)
            modal_energy(:) = 0.0_dp
            DO a = 1, m
               DO b = 1, m
                  component_weight = ABS(cluster_evec(b, a))**2
                  modal_energy(a) = modal_energy(a) + component_weight*evals(first + b - 1)
               END DO
            END DO
            CALL dense_projector_cluster_modal_occupation(evals(first:last), occ_curv(first:last), &
                                                          occ_slope(first:last), modal_energy, modal_occ, &
                                                          modal_curv)
            cluster_kernel(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
            DO a = 1, m
               cluster_kernel(a, a) = CMPLX(-2.0_dp*projector_weight*modal_curv(a)* &
                                            (cluster_eval(a) - mu_shift), 0.0_dp, KIND=dp)
            END DO
            response_block(:, :) = MATMUL(cluster_evec, &
                                          MATMUL(cluster_kernel, CONJG(TRANSPOSE(cluster_evec))))
            DO b = 1, m
               DO a = 1, m
                  kernel(first + a - 1, first + b - 1) = response_block(a, b)
               END DO
            END DO
            DO j = 1, n
               IF (j >= first .AND. j <= last) CYCLE
               outside_gauge_cluster = .FALSE.
               IF (j > 1) outside_gauge_cluster = &
                  outside_gauge_cluster .OR. ABS(evals(j) - evals(j - 1)) <= &
                  dense_projector_occupation_pair_gap(gauge_gap, occ_curv(j - 1), occ_curv(j), &
                                                      occ_curv_max, projector_window_gap)
               IF (j < n) outside_gauge_cluster = &
                  outside_gauge_cluster .OR. ABS(evals(j + 1) - evals(j)) <= &
                  dense_projector_occupation_pair_gap(gauge_gap, occ_curv(j), occ_curv(j + 1), &
                                                      occ_curv_max, projector_window_gap)
               IF (outside_gauge_cluster) CYCLE
               cluster_col(:) = work(first:last, j)
               modal_col(:) = MATMUL(CONJG(TRANSPOSE(cluster_evec)), cluster_col)
               DO a = 1, m
                  denom = modal_energy(a) - evals(j)
                  IF (ABS(denom) > gap_tol) THEN
                     offdiag_coeff = (modal_occ(a) - occ_slope(j))/denom
                  ELSE
                     offdiag_coeff = 0.0_dp
                  END IF
                  offdiag_coeff = MIN(0.0_dp, offdiag_coeff)
                  IF (PRESENT(projector_window_gap)) THEN
                     slope_cap = dense_projector_pair_slope_cap(denom, gauge_gap, modal_curv(a), occ_curv(j), &
                                                                occ_curv_max, projector_window_gap)
                  ELSE
                     slope_cap = dense_projector_pair_slope_cap(denom, gauge_gap, modal_curv(a), occ_curv(j), &
                                                                occ_curv_max)
                  END IF
                  IF (slope_cap > EPSILON(slope_cap)) THEN
                     offdiag_coeff = MAX(offdiag_coeff, -slope_cap)
                  ELSE
                     offdiag_coeff = 0.0_dp
                  END IF
                  modal_col(a) = offdiag_coeff*modal_col(a)
               END DO
               response_col(:) = projector_weight*MATMUL(cluster_evec, modal_col)
               DO a = 1, m
                  kernel(first + a - 1, j) = response_col(a)
                  kernel(j, first + a - 1) = CONJG(kernel(first + a - 1, j))
               END DO
            END DO
            DEALLOCATE (cluster_col, cluster_work, cluster_evec, cluster_kernel, modal_col, modal_energy, &
                        modal_occ, response_block, response_col, cluster_eval, modal_curv)
         ELSE
            curv = MAX(0.0_dp, occ_curv(first))
            kernel(first, first) = CMPLX(-2.0_dp*projector_weight*curv* &
                                         (REAL(work(first, first), KIND=dp) - mu_shift), &
                                         0.0_dp, KIND=dp)
         END IF
         first = last + 1
      END DO

      first = 1
      DO WHILE (first <= n)
         last = first
         DO WHILE (last < n)
            local_gap = dense_projector_occupation_pair_gap(gauge_gap, occ_curv(last), occ_curv(last + 1), &
                                                            occ_curv_max, projector_window_gap)
            IF (ABS(evals(last + 1) - evals(last)) > local_gap) EXIT
            last = last + 1
         END DO
         m = last - first + 1
         IF (m > 1) THEN
            ALLOCATE (cluster_work(m, m), cluster_evec(m, m), cluster_eval(m), modal_curv(m), &
                      modal_energy(m), modal_occ(m))
            cluster_work(:, :) = 0.5_dp*(work(first:last, first:last) + &
                                         CONJG(TRANSPOSE(work(first:last, first:last))))
            CALL diag_complex(cluster_work, cluster_evec, cluster_eval)
            CALL dense_fixed_n_projector_modal_curvature(cluster_evec, occ_curv(first:last), modal_curv)
            modal_energy(:) = 0.0_dp
            DO a = 1, m
               DO b = 1, m
                  component_weight = ABS(cluster_evec(b, a))**2
                  modal_energy(a) = modal_energy(a) + component_weight*evals(first + b - 1)
               END DO
            END DO
            CALL dense_projector_cluster_modal_occupation(evals(first:last), occ_curv(first:last), &
                                                          occ_slope(first:last), modal_energy, modal_occ, &
                                                          modal_curv)
            first_b = last + 1
            DO WHILE (first_b <= n)
               last_b = first_b
               DO WHILE (last_b < n)
                  local_gap = dense_projector_occupation_pair_gap(gauge_gap, occ_curv(last_b), &
                                                                  occ_curv(last_b + 1), occ_curv_max, &
                                                                  projector_window_gap)
                  IF (ABS(evals(last_b + 1) - evals(last_b)) > local_gap) EXIT
                  last_b = last_b + 1
               END DO
               m_b = last_b - first_b + 1
               IF (m_b > 1) THEN
                  ALLOCATE (cluster_work_b(m_b, m_b), cluster_evec_b(m_b, m_b), &
                            cluster_eval_b(m_b), modal_curv_b(m_b), modal_energy_b(m_b), &
                            modal_occ_b(m_b), modal_block(m, m_b), response_block(m, m_b))
                  cluster_work_b(:, :) = 0.5_dp*(work(first_b:last_b, first_b:last_b) + &
                                                 CONJG(TRANSPOSE(work(first_b:last_b, first_b:last_b))))
                  CALL diag_complex(cluster_work_b, cluster_evec_b, cluster_eval_b)
                  CALL dense_fixed_n_projector_modal_curvature(cluster_evec_b, occ_curv(first_b:last_b), &
                                                               modal_curv_b)
                  modal_energy_b(:) = 0.0_dp
                  DO a = 1, m_b
                     DO b = 1, m_b
                        component_weight = ABS(cluster_evec_b(b, a))**2
                        modal_energy_b(a) = modal_energy_b(a) + component_weight*evals(first_b + b - 1)
                     END DO
                  END DO
                  CALL dense_projector_cluster_modal_occupation(evals(first_b:last_b), &
                                                                occ_curv(first_b:last_b), &
                                                                occ_slope(first_b:last_b), modal_energy_b, &
                                                                modal_occ_b, modal_curv_b)
                  modal_block(:, :) = MATMUL(CONJG(TRANSPOSE(cluster_evec)), &
                                             MATMUL(work(first:last, first_b:last_b), cluster_evec_b))
                  DO b = 1, m_b
                     DO a = 1, m
                        denom = modal_energy(a) - modal_energy_b(b)
                        IF (ABS(denom) > gap_tol) THEN
                           offdiag_coeff = (modal_occ(a) - modal_occ_b(b))/denom
                        ELSE
                           offdiag_coeff = 0.0_dp
                        END IF
                        offdiag_coeff = MIN(0.0_dp, offdiag_coeff)
                        IF (PRESENT(projector_window_gap)) THEN
                           slope_cap = dense_projector_pair_slope_cap(denom, gauge_gap, modal_curv(a), &
                                                                      modal_curv_b(b), occ_curv_max, &
                                                                      projector_window_gap)
                        ELSE
                           slope_cap = dense_projector_pair_slope_cap(denom, gauge_gap, modal_curv(a), &
                                                                      modal_curv_b(b), occ_curv_max)
                        END IF
                        IF (slope_cap > EPSILON(slope_cap)) THEN
                           offdiag_coeff = MAX(offdiag_coeff, -slope_cap)
                        ELSE
                           offdiag_coeff = 0.0_dp
                        END IF
                        modal_block(a, b) = offdiag_coeff*modal_block(a, b)
                     END DO
                  END DO
                  response_block(:, :) = projector_weight*MATMUL(cluster_evec, &
                                                                 MATMUL(modal_block, CONJG(TRANSPOSE(cluster_evec_b))))
                  DO b = 1, m_b
                     DO a = 1, m
                        kernel(first + a - 1, first_b + b - 1) = response_block(a, b)
                        kernel(first_b + b - 1, first + a - 1) = CONJG(response_block(a, b))
                     END DO
                  END DO
                  DEALLOCATE (cluster_work_b, cluster_evec_b, cluster_eval_b, modal_curv_b, &
                              modal_energy_b, modal_occ_b, modal_block, response_block)
               END IF
               first_b = last_b + 1
            END DO
            DEALLOCATE (cluster_work, cluster_evec, cluster_eval, modal_curv, modal_energy, modal_occ)
         END IF
         first = last + 1
      END DO

   END SUBROUTINE dense_fixed_n_projector_cluster_kernel

! **************************************************************************************************
!> \brief Active-selective modal curvature in one near-degenerate Projector cluster.
!> \param cluster_evec ...
!> \param occ_curv ...
!> \param modal_curv ...
! **************************************************************************************************
   SUBROUTINE dense_fixed_n_projector_modal_curvature(cluster_evec, occ_curv, modal_curv)
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: cluster_evec
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: occ_curv
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: modal_curv

      REAL(KIND=dp), PARAMETER :: active_component_fraction_min = 0.60_dp, &
         active_component_majority_min = 0.50_dp, active_component_transition_min = 0.50_dp, &
         active_curvature_fraction = 0.25_dp, significant_component_fraction = 0.5_dp

      INTEGER                                            :: a, b, m
      REAL(KIND=dp) :: active_blend, active_component_fraction, component_weight, curv, &
         curv_floor, curv_max, max_component_weight, modal_damped, modal_floor, modal_raw

      m = SIZE(occ_curv)
      modal_curv(:) = 0.0_dp
      IF (m <= 0) RETURN

      curv_floor = HUGE(1.0_dp)
      curv_max = 0.0_dp
      DO b = 1, m
         curv = MAX(0.0_dp, occ_curv(b))
         IF (curv > EPSILON(curv)) curv_floor = MIN(curv_floor, curv)
         curv_max = MAX(curv_max, curv)
      END DO
      IF (curv_floor >= HUGE(1.0_dp)/2.0_dp) RETURN

      DO a = 1, m
         modal_raw = 0.0_dp
         active_component_fraction = 0.0_dp
         max_component_weight = 0.0_dp
         DO b = 1, m
            component_weight = ABS(cluster_evec(b, a))**2
            curv = MAX(0.0_dp, occ_curv(b))
            max_component_weight = MAX(max_component_weight, component_weight)
            modal_raw = modal_raw + component_weight*curv
            IF (curv >= active_curvature_fraction*curv_max) &
               active_component_fraction = active_component_fraction + component_weight
         END DO
         IF (modal_raw > EPSILON(modal_raw)) THEN
            modal_floor = HUGE(1.0_dp)
            DO b = 1, m
               component_weight = ABS(cluster_evec(b, a))**2
               curv = MAX(0.0_dp, occ_curv(b))
               IF (curv > EPSILON(curv) .AND. &
                   component_weight >= significant_component_fraction*max_component_weight) &
                  modal_floor = MIN(modal_floor, curv)
            END DO
            IF (modal_floor >= HUGE(1.0_dp)/2.0_dp) modal_floor = curv_floor
            IF (modal_raw < 0.25_dp*curv_max) modal_floor = MIN(modal_floor, curv_floor)
            modal_damped = MIN(modal_raw, SQRT(modal_floor*modal_raw))
            IF (active_component_fraction > active_component_majority_min) &
               modal_damped = MAX(modal_damped, active_component_fraction*modal_raw)
            IF (active_component_fraction >= active_component_fraction_min) THEN
               modal_curv(a) = modal_raw
            ELSE IF (active_component_fraction > active_component_transition_min .AND. &
                     active_component_fraction < active_component_fraction_min) THEN
               active_blend = (active_component_fraction - active_component_transition_min)/ &
                              (active_component_fraction_min - active_component_transition_min)
               active_blend = MIN(1.0_dp, MAX(0.0_dp, active_blend))
               modal_curv(a) = (1.0_dp - active_blend)*modal_damped + active_blend*modal_raw
            ELSE
               modal_curv(a) = modal_damped
            END IF
         END IF
      END DO

   END SUBROUTINE dense_fixed_n_projector_modal_curvature

! **************************************************************************************************
!> \brief ...
!> \param evals ...
!> \param occ_curv ...
!> \param gauge_gap ...
!> \param occ_tol ...
!> \param occ_mono ...
!> \param projector_window_gap ...
! **************************************************************************************************
   SUBROUTINE dense_near_degenerate_occupation_plateau(evals, occ_curv, gauge_gap, occ_tol, occ_mono, &
                                                       projector_window_gap)
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: evals, occ_curv
      REAL(KIND=dp), INTENT(IN)                          :: gauge_gap, occ_tol
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: occ_mono
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: projector_window_gap

      INTEGER                                            :: first, idx, last, n
      LOGICAL                                            :: needs_filter
      REAL(KIND=dp) :: allowed_jump, cluster_sum, eval_center, fit_den, fit_num, gap, local_gap, &
         occ_curv_max, occ_curv_sum, occ_jump, pair_slope, pair_weight, physical_slope, plateau, &
         shape_sum, slope_cap, slope_pair

      n = SIZE(evals)
      IF (n <= 1) RETURN
      IF (gauge_gap <= EPSILON(gauge_gap)) RETURN

      occ_curv_max = MAXVAL(MAX(0.0_dp, occ_curv(:)))
      first = 1
      DO WHILE (first < n)
         last = first
         needs_filter = .FALSE.
         DO WHILE (last < n)
            gap = ABS(evals(last + 1) - evals(last))
            local_gap = dense_projector_occupation_pair_gap(gauge_gap, occ_curv(last), occ_curv(last + 1), &
                                                            occ_curv_max, projector_window_gap)
            IF (gap > local_gap) EXIT
            slope_cap = dense_shared_occ_slope(occ_curv(last), occ_curv(last + 1))
            allowed_jump = slope_cap*MAX(gap, EPSILON(1.0_dp)) + occ_tol
            occ_jump = ABS(occ_mono(last + 1) - occ_mono(last))
            IF (occ_jump > allowed_jump) needs_filter = .TRUE.
            last = last + 1
         END DO
         IF (last > first .AND. needs_filter) THEN
            cluster_sum = SUM(occ_mono(first:last))
            plateau = cluster_sum/REAL(last - first + 1, KIND=dp)
            eval_center = SUM(evals(first:last))/REAL(last - first + 1, KIND=dp)
            slope_cap = HUGE(1.0_dp)
            fit_num = 0.0_dp
            fit_den = 0.0_dp
            DO idx = first, last - 1
               gap = ABS(evals(idx + 1) - evals(idx))
               slope_pair = dense_shared_occ_slope(occ_curv(idx), occ_curv(idx + 1))
               slope_cap = MIN(slope_cap, slope_pair)
               allowed_jump = slope_pair*MAX(gap, EPSILON(1.0_dp)) + occ_tol
               occ_jump = ABS(occ_mono(idx + 1) - occ_mono(idx))
               IF (gap <= EPSILON(gap)) CYCLE
               IF (occ_jump > allowed_jump) CYCLE
               pair_slope = (occ_mono(idx + 1) - occ_mono(idx))/(evals(idx + 1) - evals(idx))
               pair_slope = MIN(0.0_dp, pair_slope)
               pair_slope = MAX(pair_slope, -MAX(0.0_dp, slope_pair))
               pair_weight = gap*gap
               fit_num = fit_num + pair_weight*pair_slope
               fit_den = fit_den + pair_weight
            END DO
            physical_slope = 0.0_dp
            IF (fit_den > EPSILON(fit_den)) physical_slope = fit_num/fit_den
            physical_slope = MAX(physical_slope, -MAX(0.0_dp, slope_cap))
            IF (fit_den > EPSILON(fit_den) .OR. slope_cap <= EPSILON(slope_cap)) THEN
               DO idx = first, last
                  occ_mono(idx) = plateau + physical_slope*(evals(idx) - eval_center)
               END DO
            ELSE
               occ_curv_sum = 0.0_dp
               DO idx = first, last
                  occ_curv_sum = occ_curv_sum + MAX(0.0_dp, occ_curv(idx))
               END DO
               IF (occ_curv_sum > EPSILON(occ_curv_sum)) THEN
                  occ_mono(first) = 0.0_dp
                  DO idx = first + 1, last
                     slope_pair = dense_shared_occ_slope(occ_curv(idx - 1), occ_curv(idx))
                     occ_mono(idx) = occ_mono(idx - 1) - slope_pair*(evals(idx) - evals(idx - 1))
                  END DO
                  shape_sum = SUM(occ_mono(first:last))
                  occ_mono(first:last) = occ_mono(first:last) + &
                                         (cluster_sum - shape_sum)/REAL(last - first + 1, KIND=dp)
               ELSE
                  DO idx = first, last
                     occ_mono(idx) = plateau + physical_slope*(evals(idx) - eval_center)
                  END DO
               END IF
            END IF
         END IF
         first = last + 1
      END DO

   END SUBROUTINE dense_near_degenerate_occupation_plateau

! **************************************************************************************************
!> \brief Gradient of the dense spectral occupation-projector value.
!> \param x_re ...
!> \param x_im ...
!> \param base ...
!> \param occ_scale ...
!> \param grad_re ...
!> \param grad_im ...
!> \param eps_diag ...
! **************************************************************************************************
   SUBROUTINE dense_spectral_projector_gradient(x_re, x_im, base, occ_scale, grad_re, grad_im, eps_diag)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: x_re, x_im
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: base
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: occ_scale
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: grad_re, grad_im
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: eps_diag

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: occ_mat
      INTEGER                                            :: n

      n = SIZE(x_re, 1)
      ALLOCATE (occ_mat(n, n))
      CALL dense_spectral_occ_projector(x_re, x_im, base, occ_scale, occ_mat)
      CALL dense_projector_rotation_gradient(x_re, x_im, base, occ_mat, grad_re, grad_im, eps_diag)
      DEALLOCATE (occ_mat)

   END SUBROUTINE dense_spectral_projector_gradient

! **************************************************************************************************
!> \brief Dense value 1/2 Tr[P U^H H U] for a doubled Hermitian occupation projector.
!> \param x_re ...
!> \param x_im ...
!> \param base ...
!> \param occ_mat ...
!> \return ...
! **************************************************************************************************
   FUNCTION dense_projector_rotation_value(x_re, x_im, base, occ_mat) RESULT(value)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: x_re, x_im
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: base, occ_mat
      REAL(KIND=dp)                                      :: value

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: chc, u
      INTEGER                                            :: n

      n = SIZE(x_re, 1)
      ALLOCATE (chc(n, n), u(n, n))
      CALL dense_rotation_u(x_re, x_im, u)
      chc(:, :) = MATMUL(CONJG(TRANSPOSE(u)), MATMUL(base, u))
      value = 0.5_dp*REAL(SUM(CONJG(occ_mat(:, :))*chc(:, :)), KIND=dp)
      DEALLOCATE (chc, u)

   END FUNCTION dense_projector_rotation_value

! **************************************************************************************************
!> \brief Dense gradient of 1/2 Tr[P U^H H U] for a doubled Hermitian occupation projector.
!> \param x_re ...
!> \param x_im ...
!> \param base ...
!> \param occ_mat ...
!> \param grad_re ...
!> \param grad_im ...
!> \param eps_diag ...
! **************************************************************************************************
   SUBROUTINE dense_projector_rotation_gradient(x_re, x_im, base, occ_mat, grad_re, grad_im, eps_diag)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: x_re, x_im
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: base, occ_mat
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
      CALL dense_rotation_u_from_eigensystem(evec, evals, u)

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

   END SUBROUTINE dense_projector_rotation_gradient

! **************************************************************************************************
!> \brief Check the complex k-point NONDIAG_ENERGY tether gradient.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_complex_nondiag_tether_gradient(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 3
      REAL(KIND=dp), PARAMETER                           :: fd_step = 1.0E-6_dp, kweight = 0.75_dp, &
                                                            strength = 0.6_dp

      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: base
      INTEGER                                            :: i, j
      REAL(KIND=dp)                                      :: err, fd, value_m, value_p
      REAL(KIND=dp), DIMENSION(n)                        :: eps, grad_eps
      REAL(KIND=dp), DIMENSION(n, n)                     :: grad_im, grad_re, x_im, x_m_im, x_m_re, &
                                                            x_p_im, x_p_re, x_re

      base(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      base(1, 1) = CMPLX(-0.30_dp, 0.0_dp, KIND=dp)
      base(2, 2) = CMPLX(-0.02_dp, 0.0_dp, KIND=dp)
      base(3, 3) = CMPLX(0.12_dp, 0.0_dp, KIND=dp)
      base(1, 2) = CMPLX(0.04_dp, -0.03_dp, KIND=dp)
      base(2, 1) = CONJG(base(1, 2))
      base(1, 3) = CMPLX(-0.02_dp, 0.05_dp, KIND=dp)
      base(3, 1) = CONJG(base(1, 3))
      base(2, 3) = CMPLX(0.03_dp, 0.01_dp, KIND=dp)
      base(3, 2) = CONJG(base(2, 3))

      x_re(:, :) = 0.0_dp
      x_re(1, 2) = 0.08_dp
      x_re(2, 1) = -x_re(1, 2)
      x_re(1, 3) = -0.04_dp
      x_re(3, 1) = -x_re(1, 3)
      x_re(2, 3) = 0.02_dp
      x_re(3, 2) = -x_re(2, 3)

      x_im(:, :) = 0.0_dp
      x_im(1, 1) = 0.03_dp
      x_im(2, 2) = -0.05_dp
      x_im(3, 3) = 0.01_dp
      x_im(1, 2) = -0.02_dp
      x_im(2, 1) = x_im(1, 2)
      x_im(1, 3) = 0.01_dp
      x_im(3, 1) = x_im(1, 3)
      x_im(2, 3) = 0.04_dp
      x_im(3, 2) = x_im(2, 3)

      eps = [-0.24_dp, -0.01_dp, 0.09_dp]
      CALL dense_nondiag_tether_gradient(x_re, x_im, base, eps, strength, kweight, &
                                         grad_re, grad_im, grad_eps, value_p)

      err = 0.0_dp
      DO i = 1, n
         eps(i) = eps(i) + fd_step
         value_p = dense_nondiag_tether_value(x_re, x_im, base, eps, strength, kweight)
         eps(i) = eps(i) - 2.0_dp*fd_step
         value_m = dense_nondiag_tether_value(x_re, x_im, base, eps, strength, kweight)
         eps(i) = eps(i) + fd_step
         fd = (value_p - value_m)/(2.0_dp*fd_step)
         err = MAX(err, ABS(fd - grad_eps(i)))
      END DO

      DO j = 2, n
         DO i = 1, j - 1
            x_p_re(:, :) = x_re(:, :)
            x_m_re(:, :) = x_re(:, :)
            x_p_im(:, :) = x_im(:, :)
            x_m_im(:, :) = x_im(:, :)
            x_p_re(i, j) = x_p_re(i, j) + fd_step
            x_p_re(j, i) = x_p_re(j, i) - fd_step
            x_m_re(i, j) = x_m_re(i, j) - fd_step
            x_m_re(j, i) = x_m_re(j, i) + fd_step
            value_p = dense_nondiag_tether_value(x_p_re, x_p_im, base, eps, strength, kweight)
            value_m = dense_nondiag_tether_value(x_m_re, x_m_im, base, eps, strength, kweight)
            fd = (value_p - value_m)/(2.0_dp*fd_step)
            err = MAX(err, ABS(fd - grad_re(i, j)))

            x_p_re(:, :) = x_re(:, :)
            x_m_re(:, :) = x_re(:, :)
            x_p_im(:, :) = x_im(:, :)
            x_m_im(:, :) = x_im(:, :)
            x_p_im(i, j) = x_p_im(i, j) + fd_step
            x_p_im(j, i) = x_p_im(j, i) + fd_step
            x_m_im(i, j) = x_m_im(i, j) - fd_step
            x_m_im(j, i) = x_m_im(j, i) - fd_step
            value_p = dense_nondiag_tether_value(x_p_re, x_p_im, base, eps, strength, kweight)
            value_m = dense_nondiag_tether_value(x_m_re, x_m_im, base, eps, strength, kweight)
            fd = (value_p - value_m)/(2.0_dp*fd_step)
            err = MAX(err, ABS(fd - grad_im(i, j)))
         END DO
      END DO

      DO i = 1, n
         x_p_re(:, :) = x_re(:, :)
         x_m_re(:, :) = x_re(:, :)
         x_p_im(:, :) = x_im(:, :)
         x_m_im(:, :) = x_im(:, :)
         x_p_im(i, i) = x_p_im(i, i) + fd_step
         x_m_im(i, i) = x_m_im(i, i) - fd_step
         value_p = dense_nondiag_tether_value(x_p_re, x_p_im, base, eps, strength, kweight)
         value_m = dense_nondiag_tether_value(x_m_re, x_m_im, base, eps, strength, kweight)
         fd = (value_p - value_m)/(2.0_dp*fd_step)
         err = MAX(err, ABS(fd - 0.5_dp*grad_im(i, i)))
      END DO

      IF (err > tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,ES12.5)') TRIM(label), "failed with error", err
      ELSE
         WRITE (*, '(A,1X,A,1X,ES12.5)') TRIM(label), "OK with error", err
      END IF

   END SUBROUTINE check_complex_nondiag_tether_gradient

! **************************************************************************************************
!> \brief Check the gauge-projected complex canonical tether gradient.
!> \param label ...
!> \param tol ...
!> \param nfail ...
! **************************************************************************************************
   SUBROUTINE check_complex_canonical_tether_gradient(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 3
      REAL(KIND=dp), PARAMETER                           :: fd_step = 1.0E-6_dp

      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: base
      INTEGER                                            :: i, j
      REAL(KIND=dp)                                      :: err, fd, value_m, value_p
      REAL(KIND=dp), DIMENSION(n)                        :: eps, grad_eps, weights
      REAL(KIND=dp), DIMENSION(n, n)                     :: grad_im, grad_re, x_im, x_m_im, x_m_re, &
                                                            x_p_im, x_p_re, x_re

      base(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      base(1, 1) = CMPLX(-0.31_dp, 0.0_dp, KIND=dp)
      base(2, 2) = CMPLX(-0.04_dp, 0.0_dp, KIND=dp)
      base(3, 3) = CMPLX(0.17_dp, 0.0_dp, KIND=dp)
      base(1, 2) = CMPLX(0.05_dp, -0.02_dp, KIND=dp)
      base(2, 1) = CONJG(base(1, 2))
      base(1, 3) = CMPLX(-0.03_dp, 0.04_dp, KIND=dp)
      base(3, 1) = CONJG(base(1, 3))
      base(2, 3) = CMPLX(0.02_dp, 0.03_dp, KIND=dp)
      base(3, 2) = CONJG(base(2, 3))

      x_re(:, :) = 0.0_dp
      x_re(1, 2) = 0.07_dp
      x_re(2, 1) = -x_re(1, 2)
      x_re(1, 3) = -0.03_dp
      x_re(3, 1) = -x_re(1, 3)
      x_re(2, 3) = 0.04_dp
      x_re(3, 2) = -x_re(2, 3)

      x_im(:, :) = 0.0_dp
      x_im(1, 1) = 0.02_dp
      x_im(2, 2) = -0.04_dp
      x_im(3, 3) = 0.03_dp
      x_im(1, 2) = -0.01_dp
      x_im(2, 1) = x_im(1, 2)
      x_im(1, 3) = 0.03_dp
      x_im(3, 1) = x_im(1, 3)
      x_im(2, 3) = 0.02_dp
      x_im(3, 2) = x_im(2, 3)

      eps = [-0.27_dp, -0.03_dp, 0.13_dp]
      weights = [0.42_dp, 0.73_dp, 0.28_dp]
      CALL dense_canonical_tether_gradient(x_re, x_im, base, eps, weights, &
                                           grad_re, grad_im, grad_eps, value_p)

      err = 0.0_dp
      DO i = 1, n
         eps(i) = eps(i) + fd_step
         value_p = dense_canonical_tether_value(x_re, x_im, base, eps, weights)
         eps(i) = eps(i) - 2.0_dp*fd_step
         value_m = dense_canonical_tether_value(x_re, x_im, base, eps, weights)
         eps(i) = eps(i) + fd_step
         fd = (value_p - value_m)/(2.0_dp*fd_step)
         err = MAX(err, ABS(fd - grad_eps(i)))
      END DO

      DO j = 2, n
         DO i = 1, j - 1
            x_p_re(:, :) = x_re(:, :)
            x_m_re(:, :) = x_re(:, :)
            x_p_im(:, :) = x_im(:, :)
            x_m_im(:, :) = x_im(:, :)
            x_p_re(i, j) = x_p_re(i, j) + fd_step
            x_p_re(j, i) = x_p_re(j, i) - fd_step
            x_m_re(i, j) = x_m_re(i, j) - fd_step
            x_m_re(j, i) = x_m_re(j, i) + fd_step
            value_p = dense_canonical_tether_value(x_p_re, x_p_im, base, eps, weights)
            value_m = dense_canonical_tether_value(x_m_re, x_m_im, base, eps, weights)
            fd = (value_p - value_m)/(2.0_dp*fd_step)
            err = MAX(err, ABS(fd - grad_re(i, j)))

            x_p_re(:, :) = x_re(:, :)
            x_m_re(:, :) = x_re(:, :)
            x_p_im(:, :) = x_im(:, :)
            x_m_im(:, :) = x_im(:, :)
            x_p_im(i, j) = x_p_im(i, j) + fd_step
            x_p_im(j, i) = x_p_im(j, i) + fd_step
            x_m_im(i, j) = x_m_im(i, j) - fd_step
            x_m_im(j, i) = x_m_im(j, i) - fd_step
            value_p = dense_canonical_tether_value(x_p_re, x_p_im, base, eps, weights)
            value_m = dense_canonical_tether_value(x_m_re, x_m_im, base, eps, weights)
            fd = (value_p - value_m)/(2.0_dp*fd_step)
            err = MAX(err, ABS(fd - grad_im(i, j)))
         END DO
      END DO

      DO i = 1, n
         x_p_re(:, :) = x_re(:, :)
         x_m_re(:, :) = x_re(:, :)
         x_p_im(:, :) = x_im(:, :)
         x_m_im(:, :) = x_im(:, :)
         x_p_im(i, i) = x_p_im(i, i) + fd_step
         x_m_im(i, i) = x_m_im(i, i) - fd_step
         value_p = dense_canonical_tether_value(x_p_re, x_p_im, base, eps, weights)
         value_m = dense_canonical_tether_value(x_m_re, x_m_im, base, eps, weights)
         fd = (value_p - value_m)/(2.0_dp*fd_step)
         err = MAX(err, ABS(fd - 0.5_dp*grad_im(i, i)))
      END DO

      IF (err > tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,ES12.5)') TRIM(label), "failed with error", err
      ELSE
         WRITE (*, '(A,1X,A,1X,ES12.5)') TRIM(label), "OK with error", err
      END IF

   END SUBROUTINE check_complex_canonical_tether_gradient

! **************************************************************************************************
!> \brief Dense gauge-projected complex canonical tether value.
!> \param x_re ...
!> \param x_im ...
!> \param base ...
!> \param eps ...
!> \param weights ...
!> \return ...
! **************************************************************************************************
   FUNCTION dense_canonical_tether_value(x_re, x_im, base, eps, weights) RESULT(value)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: x_re, x_im
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: base
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: eps, weights
      REAL(KIND=dp)                                      :: value

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: chc, u
      INTEGER                                            :: i, n
      REAL(KIND=dp)                                      :: shift, weight_sum
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: delta

      n = SIZE(x_re, 1)
      ALLOCATE (chc(n, n), u(n, n), delta(n))
      CALL dense_rotation_u(x_re, x_im, u)
      chc(:, :) = MATMUL(CONJG(TRANSPOSE(u)), MATMUL(base, u))
      weight_sum = SUM(weights)
      shift = 0.0_dp
      DO i = 1, n
         delta(i) = REAL(chc(i, i), KIND=dp) - eps(i)
         shift = shift + weights(i)*delta(i)
      END DO
      IF (weight_sum > EPSILON(weight_sum)) shift = shift/weight_sum
      delta(:) = delta(:) - shift
      value = 0.5_dp*SUM(weights(:)*delta(:)**2)
      DEALLOCATE (chc, u, delta)

   END FUNCTION dense_canonical_tether_value

! **************************************************************************************************
!> \brief Dense gauge-projected complex canonical tether gradients.
!> \param x_re ...
!> \param x_im ...
!> \param base ...
!> \param eps ...
!> \param weights ...
!> \param grad_re ...
!> \param grad_im ...
!> \param grad_eps ...
!> \param value ...
! **************************************************************************************************
   SUBROUTINE dense_canonical_tether_gradient(x_re, x_im, base, eps, weights, &
                                              grad_re, grad_im, grad_eps, value)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: x_re, x_im
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: base
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: eps, weights
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: grad_re, grad_im
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: grad_eps
      REAL(KIND=dp), INTENT(OUT)                         :: value

      COMPLEX(KIND=dp)                                   :: kernel
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: chc, dedu, dmat, evec, frechet, &
                                                            generator, inner, outer, u
      INTEGER                                            :: i, j, n
      REAL(KIND=dp)                                      :: shift, weight_sum
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: delta, evals

      n = SIZE(x_re, 1)
      ALLOCATE (chc(n, n), dedu(n, n), dmat(n, n), evec(n, n), frechet(n, n), generator(n, n))
      ALLOCATE (inner(n, n), outer(n, n), u(n, n), delta(n), evals(n))

      generator(:, :) = CMPLX(-x_im(:, :), x_re(:, :), KIND=dp)
      CALL diag_complex(generator, evec, evals)
      CALL dense_rotation_u_from_eigensystem(evec, evals, u)

      chc(:, :) = MATMUL(CONJG(TRANSPOSE(u)), MATMUL(base, u))
      weight_sum = SUM(weights)
      shift = 0.0_dp
      DO i = 1, n
         delta(i) = REAL(chc(i, i), KIND=dp) - eps(i)
         shift = shift + weights(i)*delta(i)
      END DO
      IF (weight_sum > EPSILON(weight_sum)) shift = shift/weight_sum
      delta(:) = delta(:) - shift
      value = 0.5_dp*SUM(weights(:)*delta(:)**2)
      grad_eps(:) = -weights(:)*delta(:)

      dmat(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      DO i = 1, n
         dmat(i, i) = CMPLX(weights(i)*delta(i), 0.0_dp, KIND=dp)
      END DO
      dedu(:, :) = 2.0_dp*MATMUL(base, MATMUL(u, dmat))

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

      DEALLOCATE (chc, dedu, dmat, evec, frechet, generator, inner, outer, u, delta, evals)

   END SUBROUTINE dense_canonical_tether_gradient

! **************************************************************************************************
!> \brief Dense complex NONDIAG_ENERGY tether value.
!> \param x_re ...
!> \param x_im ...
!> \param base ...
!> \param eps ...
!> \param strength ...
!> \param kweight ...
!> \return ...
! **************************************************************************************************
   FUNCTION dense_nondiag_tether_value(x_re, x_im, base, eps, strength, kweight) RESULT(value)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: x_re, x_im
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: base
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: eps
      REAL(KIND=dp), INTENT(IN)                          :: strength, kweight
      REAL(KIND=dp)                                      :: value

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: dmat, u
      INTEGER                                            :: i, n

      n = SIZE(x_re, 1)
      ALLOCATE (dmat(n, n), u(n, n))
      CALL dense_rotation_u(x_re, x_im, u)
      dmat(:, :) = MATMUL(CONJG(TRANSPOSE(u)), MATMUL(base, u))
      DO i = 1, n
         dmat(i, i) = dmat(i, i) - eps(i)
      END DO
      value = 0.5_dp*strength*kweight*SUM(ABS(dmat(:, :))**2)
      DEALLOCATE (dmat, u)

   END FUNCTION dense_nondiag_tether_value

! **************************************************************************************************
!> \brief Dense complex NONDIAG_ENERGY tether gradients.
!> \param x_re ...
!> \param x_im ...
!> \param base ...
!> \param eps ...
!> \param strength ...
!> \param kweight ...
!> \param grad_re ...
!> \param grad_im ...
!> \param grad_eps ...
!> \param value ...
! **************************************************************************************************
   SUBROUTINE dense_nondiag_tether_gradient(x_re, x_im, base, eps, strength, kweight, &
                                            grad_re, grad_im, grad_eps, value)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: x_re, x_im
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: base
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: eps
      REAL(KIND=dp), INTENT(IN)                          :: strength, kweight
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: grad_re, grad_im
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: grad_eps
      REAL(KIND=dp), INTENT(OUT)                         :: value

      COMPLEX(KIND=dp)                                   :: kernel
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: dedu, dmat, evec, frechet, generator, &
                                                            inner, outer, u
      INTEGER                                            :: i, j, n
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: evals

      n = SIZE(x_re, 1)
      ALLOCATE (dmat(n, n), dedu(n, n), evec(n, n), frechet(n, n), generator(n, n))
      ALLOCATE (inner(n, n), outer(n, n), u(n, n), evals(n))

      generator(:, :) = CMPLX(-x_im(:, :), x_re(:, :), KIND=dp)
      CALL diag_complex(generator, evec, evals)
      CALL dense_rotation_u_from_eigensystem(evec, evals, u)

      dmat(:, :) = MATMUL(CONJG(TRANSPOSE(u)), MATMUL(base, u))
      DO i = 1, n
         dmat(i, i) = dmat(i, i) - eps(i)
      END DO
      value = 0.5_dp*strength*kweight*SUM(ABS(dmat(:, :))**2)
      DO i = 1, n
         grad_eps(i) = -strength*kweight*REAL(dmat(i, i), KIND=dp)
      END DO

      dedu(:, :) = 2.0_dp*strength*kweight*MATMUL(base, MATMUL(u, dmat))
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

      DEALLOCATE (dmat, dedu, evec, frechet, generator, inner, outer, u, evals)

   END SUBROUTINE dense_nondiag_tether_gradient

! **************************************************************************************************
!> \brief Dense gradient of sum_i occ_scale_i [U^H H U]_ii.
!> \param x_re ...
!> \param x_im ...
!> \param base ...
!> \param occ_scale ...
!> \param grad_re ...
!> \param grad_im ...
!> \param eps_diag ...
! **************************************************************************************************
   SUBROUTINE dense_ref_rotation_gradient(x_re, x_im, base, occ_scale, grad_re, grad_im, eps_diag)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: x_re, x_im
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: base
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: occ_scale
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
      CALL dense_rotation_u_from_eigensystem(evec, evals, u)

      work(:, :) = MATMUL(base, u)
      dedu(:, :) = work(:, :)
      DO j = 1, n
         dedu(:, j) = occ_scale(j)*dedu(:, j)
      END DO

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

! **************************************************************************************************
!> \brief Build U=exp(-i X) for the dense anti-Hermitian k-point generator convention.
!> \param x_re ...
!> \param x_im ...
!> \param u ...
! **************************************************************************************************
   SUBROUTINE dense_rotation_u(x_re, x_im, u)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: x_re, x_im
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(OUT)     :: u

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: evec, generator
      INTEGER                                            :: n
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: evals

      n = SIZE(x_re, 1)
      ALLOCATE (evec(n, n), generator(n, n), evals(n))
      generator(:, :) = CMPLX(-x_im(:, :), x_re(:, :), KIND=dp)
      CALL diag_complex(generator, evec, evals)
      CALL dense_rotation_u_from_eigensystem(evec, evals, u)
      DEALLOCATE (evec, generator, evals)

   END SUBROUTINE dense_rotation_u

! **************************************************************************************************
!> \brief Build U=V exp(-i evals) V^H.
!> \param evec ...
!> \param evals ...
!> \param u ...
! **************************************************************************************************
   SUBROUTINE dense_rotation_u_from_eigensystem(evec, evals, u)
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: evec
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: evals
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(OUT)     :: u

      INTEGER                                            :: j

      u(:, :) = evec(:, :)
      DO j = 1, SIZE(evals)
         u(:, j) = EXP(CMPLX(0.0_dp, -evals(j), KIND=dp))*u(:, j)
      END DO
      u(:, :) = MATMUL(u, CONJG(TRANSPOSE(evec)))

   END SUBROUTINE dense_rotation_u_from_eigensystem

! **************************************************************************************************
!> \brief Small local Fermi-Dirac k-point smearing helper for this standalone unit test.
!> \param occ ...
!> \param mu ...
!> \param kTS ...
!> \param eps ...
!> \param nel ...
!> \param wk ...
!> \param sigma ...
!> \param maxocc ...
! **************************************************************************************************
   SUBROUTINE dense_fermi_kpoint(occ, mu, kTS, eps, nel, wk, sigma, maxocc)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: occ
      REAL(KIND=dp), INTENT(OUT)                         :: mu, kTS
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: eps
      REAL(KIND=dp), INTENT(IN)                          :: nel
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: wk
      REAL(KIND=dp), INTENT(IN)                          :: sigma, maxocc

      REAL(KIND=dp), PARAMETER                           :: epsocc = 1.0E-12_dp

      INTEGER                                            :: iter
      REAL(KIND=dp)                                      :: de, mu_max, mu_min, n_now

      de = MAX(sigma*LOG((1.0_dp - epsocc)/epsocc), 0.5_dp)
      mu_min = MINVAL(eps) - de
      mu_max = MAXVAL(eps) + de
      DO iter = 1, 400
         mu = 0.5_dp*(mu_min + mu_max)
         CALL dense_fermi_kpoint_eval(occ, n_now, kTS, eps, mu, wk, sigma, maxocc)
         IF (ABS(n_now - nel) < nel*epsocc) EXIT
         IF (n_now <= nel) THEN
            mu_min = mu
         ELSE
            mu_max = mu
         END IF
      END DO
      mu = 0.5_dp*(mu_min + mu_max)
      CALL dense_fermi_kpoint_eval(occ, n_now, kTS, eps, mu, wk, sigma, maxocc)

   END SUBROUTINE dense_fermi_kpoint

! **************************************************************************************************
!> \brief Evaluate Fermi-Dirac occupations and -TS correction for a fixed chemical potential.
!> \param occ ...
!> \param nel ...
!> \param kTS ...
!> \param eps ...
!> \param mu ...
!> \param wk ...
!> \param sigma ...
!> \param maxocc ...
! **************************************************************************************************
   SUBROUTINE dense_fermi_kpoint_eval(occ, nel, kTS, eps, mu, wk, sigma, maxocc)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: occ
      REAL(KIND=dp), INTENT(OUT)                         :: nel, kTS
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: eps
      REAL(KIND=dp), INTENT(IN)                          :: mu
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: wk
      REAL(KIND=dp), INTENT(IN)                          :: sigma, maxocc

      INTEGER                                            :: ik, imo
      REAL(KIND=dp)                                      :: arg, term1, term2, tmp, tmp2, tmp3, &
                                                            tmp4, tmplog

      nel = 0.0_dp
      kTS = 0.0_dp
      DO ik = 1, SIZE(eps, 2)
         DO imo = 1, SIZE(eps, 1)
            IF (eps(imo, ik) > mu) THEN
               arg = -(eps(imo, ik) - mu)/sigma
               tmp = EXP(arg)
               tmp4 = tmp + 1.0_dp
               tmp2 = tmp/tmp4
               tmp3 = 1.0_dp/tmp4
               tmplog = -LOG(tmp4)
               term1 = tmp2*(arg + tmplog)
               term2 = tmp3*tmplog
            ELSE
               arg = (eps(imo, ik) - mu)/sigma
               tmp = EXP(arg)
               tmp4 = tmp + 1.0_dp
               tmp2 = 1.0_dp/tmp4
               tmp3 = tmp/tmp4
               tmplog = -LOG(tmp4)
               term1 = tmp2*tmplog
               term2 = tmp3*(arg + tmplog)
            END IF
            occ(imo, ik) = maxocc*tmp2
            nel = nel + wk(ik)*occ(imo, ik)
            kTS = kTS + sigma*maxocc*(term1 + term2)*wk(ik)
         END DO
      END DO

   END SUBROUTINE dense_fermi_kpoint_eval

! **************************************************************************************************
!> \brief Small dense 2x2 matrix exponential via Taylor series.
!> \param a ...
!> \param expa ...
! **************************************************************************************************
   SUBROUTINE dense_expm2(a, expa)
      COMPLEX(KIND=dp), DIMENSION(2, 2), INTENT(IN)      :: a
      COMPLEX(KIND=dp), DIMENSION(2, 2), INTENT(OUT)     :: expa

      COMPLEX(KIND=dp), DIMENSION(2, 2)                  :: term
      INTEGER                                            :: i, n

      expa(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      term(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      DO i = 1, 2
         expa(i, i) = CMPLX(1.0_dp, 0.0_dp, KIND=dp)
         term(i, i) = CMPLX(1.0_dp, 0.0_dp, KIND=dp)
      END DO

      DO n = 1, 80
         term = MATMUL(term, a)/REAL(n, KIND=dp)
         expa = expa + term
         IF (MAXVAL(ABS(term)) < 1.0E-20_dp) EXIT
      END DO

   END SUBROUTINE dense_expm2

END PROGRAM qs_ot_frechet_unittest
