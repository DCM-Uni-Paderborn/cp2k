!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Experimental CP2K-native GPW real-space-grid path for SKALA TorchScript models.
! **************************************************************************************************
MODULE skala_gpw_functional
   USE OMP_LIB,                         ONLY: omp_get_max_threads,&
                                              omp_set_num_threads
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_array_utils,                  ONLY: cp_3d_r_cp_type
   USE cp_log_handling,                 ONLY: cp_logger_get_default_io_unit
   USE input_section_types,             ONLY: section_get_rval,&
                                              section_vals_get_subs_vals,&
                                              section_vals_get_subs_vals2,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_path_length,&
                                              dp,&
                                              int_8
   USE message_passing,                 ONLY: mp_comm_type
   USE offload_api,                     ONLY: offload_set_chosen_device
   USE particle_types,                  ONLY: particle_type
   USE pw_grid_types,                   ONLY: pw_grid_type
   USE pw_methods,                      ONLY: pw_scale,&
                                              pw_zero
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE qs_grid_atom,                    ONLY: grid_atom_type
   USE skala_density_constraints,       ONLY: &
        skala_atom_domain, skala_background_domain, skala_density_constraints_count, &
        skala_density_constraints_get, skala_density_constraints_matches, &
        skala_density_constraints_type, skala_density_descriptor, skala_density_jumps_active, &
        skala_density_jumps_put, skala_density_samples_active, skala_density_samples_only, &
        skala_density_samples_put, skala_fft_domain, skala_kinetic_descriptor
   USE skala_gpw_features,              ONLY: &
        skala_gpw_atom_partition_hard, skala_gpw_atom_partition_smooth, &
        skala_gpw_atom_subchunk_count, skala_gpw_atom_subchunk_layout, skala_gpw_feature_build, &
        skala_gpw_feature_build_atom_subchunk_bounds, skala_gpw_feature_release, &
        skala_gpw_feature_type, skala_gpw_periodic_partition_response, &
        skala_gpw_periodic_weight_derivative, skala_gpw_smooth_partition_derivatives, &
        smooth_partition_atomic_weight_scale_derivative
   USE skala_torch_api,                 ONLY: skala_torch_model_density_floors,&
                                              skala_torch_model_descriptor_limits,&
                                              skala_torch_model_get_exc,&
                                              skala_torch_model_load,&
                                              skala_torch_model_release,&
                                              skala_torch_model_type
   USE string_utilities,                ONLY: uppercase
   USE torch_api,                       ONLY: &
        torch_cuda_device_count, torch_cuda_is_available, torch_dict_create, torch_dict_insert, &
        torch_dict_release, torch_dict_type, torch_tensor_backward_scalar, torch_tensor_data_ptr, &
        torch_tensor_from_array, torch_tensor_grad, torch_tensor_grad_batch3, &
        torch_tensor_hessian_product3, torch_tensor_narrow, torch_tensor_release, &
        torch_tensor_reset_from_array, torch_tensor_to_device_leaf, torch_tensor_type, &
        torch_use_cuda
   USE xc_input_constants,              ONLY: skala_gapw_atom_composite_grid,&
                                              skala_gapw_common_grid,&
                                              skala_gapw_cp2k_default,&
                                              skala_gapw_direct_valence,&
                                              skala_gapw_paw_one_center,&
                                              skala_gapw_paw_one_center_split
   USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_type
   USE xc_rho_set_types,                ONLY: xc_rho_set_create,&
                                              xc_rho_set_get,&
                                              xc_rho_set_release,&
                                              xc_rho_set_type,&
                                              xc_rho_set_update
   USE xc_util,                         ONLY: xc_pw_divergence,&
                                              xc_requires_tmp_g
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'skala_gpw_functional'
   INTEGER, PARAMETER, PRIVATE          :: atom_chunk_auto_max_rows = 400000, &
                                           atom_chunk_auto_min_rows = 100000, &
                                           atom_chunk_auto_row_quantum = 100000, &
                                           ncollapsed_grad_per_point = 5, ngrad_per_point = 10
   INTEGER, PARAMETER, PUBLIC           :: skala_gapw_density_partition_hard_minus_soft = 1, &
                                           skala_gapw_density_partition_hard_only = 2, &
                                           skala_gapw_density_partition_soft_only = 3, &
                                           skala_gapw_density_partition_none = 4

   TYPE skala_gapw_atom_cuda_tensor_cache_type
      TYPE(torch_tensor_type)                            :: atomic_coords_t, atomic_grid_weights_t, &
                                                            density_t, grad_t, grid_coords_t, &
                                                            grid_weights_t, kin_t
   END TYPE skala_gapw_atom_cuda_tensor_cache_type

   PUBLIC :: native_skala_spin_floors
   PUBLIC :: build_vxc_from_feature_grads, ensure_native_skala_grid_scope, get_gauxc_section, &
             skala_gapw_atom_composite_energy, skala_gapw_atom_vxc_of_r, &
             native_skala_gapw_atom_composite_requested, &
             native_skala_uses_atom_composite_grid, &
             native_skala_gapw_composite_direct_ao, native_skala_gapw_composite_reference, &
             native_skala_gapw_density_partition, skala_gpw_eval, skala_gpw_weight_derivative, &
             skala_gapw_representation, xc_section_uses_native_skala_evaluator, &
             xc_section_uses_native_skala_grid, xc_section_uses_gauxc_model, &
             xc_section_uses_skala_model

   TYPE(skala_torch_model_type), SAVE                  :: cached_model
   TYPE(skala_gapw_atom_cuda_tensor_cache_type), SAVE, TARGET :: cached_atom_cuda_tensors
   CHARACTER(len=default_path_length), SAVE            :: cached_model_path = ""
   LOGICAL, SAVE                                       :: cached_model_loaded = .FALSE.
   LOGICAL, SAVE                                       :: cached_model_fp64 = .FALSE.
   INTEGER, SAVE                                       :: cached_model_cuda_device = -3
   INTEGER, SAVE                                       :: logged_cuda_device = -3, &
                                                          logged_cuda_device_count = -1, &
                                                          logged_cuda_nproc = -1, &
                                                          logged_cuda_request = -3

CONTAINS

! **************************************************************************************************
!> \brief Return true if the GAUXC subsection requests the CP2K-native GPW grid path.
!> \param xc_section ...
!> \return ...
! **************************************************************************************************
   FUNCTION xc_section_uses_native_skala_grid(xc_section) RESULT(uses_native_grid)
      TYPE(section_vals_type), INTENT(IN), POINTER       :: xc_section
      LOGICAL                                            :: uses_native_grid

      TYPE(section_vals_type), POINTER                   :: gauxc_section

      uses_native_grid = .FALSE.
      gauxc_section => get_gauxc_section(xc_section)
      IF (ASSOCIATED(gauxc_section)) THEN
         CALL section_vals_val_get(gauxc_section, "NATIVE_GRID", l_val=uses_native_grid)
      END IF

   END FUNCTION xc_section_uses_native_skala_grid

! **************************************************************************************************
!> \brief Return the pseudopotential GAPW representation selected for an active model.
!> \param xc_section ...
!> \return ...
! **************************************************************************************************
   FUNCTION skala_gapw_representation(xc_section) RESULT(representation)
      TYPE(section_vals_type), INTENT(IN), POINTER       :: xc_section
      INTEGER                                            :: representation

      TYPE(section_vals_type), POINTER                   :: gauxc_section

      representation = skala_gapw_cp2k_default
      IF (.NOT. xc_section_uses_gauxc_model(xc_section)) RETURN

      gauxc_section => get_gauxc_section(xc_section)
      CPASSERT(ASSOCIATED(gauxc_section))
      CALL section_vals_val_get(gauxc_section, "PSEUDOPOTENTIAL_GAPW_REPRESENTATION", &
                                i_val=representation)

      SELECT CASE (representation)
      CASE (skala_gapw_direct_valence, skala_gapw_paw_one_center, &
            skala_gapw_cp2k_default, skala_gapw_paw_one_center_split)
         CONTINUE
      CASE DEFAULT
         CALL cp_abort(__LOCATION__, "Unknown pseudopotential GAPW representation.")
      END SELECT

   END FUNCTION skala_gapw_representation

! **************************************************************************************************
!> \brief Return true when SKALA must be evaluated by the CP2K-native grid machinery.
!> \param xc_section ...
!> \return ...
! **************************************************************************************************
   FUNCTION xc_section_uses_native_skala_evaluator(xc_section) RESULT(uses_native_evaluator)
      TYPE(section_vals_type), INTENT(IN), POINTER       :: xc_section
      LOGICAL                                            :: uses_native_evaluator

      uses_native_evaluator = xc_section_uses_native_skala_grid(xc_section) .OR. &
                              skala_gapw_representation(xc_section) == &
                              skala_gapw_paw_one_center

   END FUNCTION xc_section_uses_native_skala_evaluator

! **************************************************************************************************
!> \brief Return true if native SKALA should use the full GAPW ORB density on one common grid.
!> \param xc_section ...
!> \return ...
! **************************************************************************************************
   FUNCTION native_skala_gapw_composite_reference(xc_section) RESULT(use_composite_reference)
      TYPE(section_vals_type), INTENT(IN), POINTER       :: xc_section
      LOGICAL                                            :: use_composite_reference

      LOGICAL                                            :: native_grid
      TYPE(section_vals_type), POINTER                   :: gauxc_section

      use_composite_reference = .FALSE.
      native_grid = .FALSE.
      gauxc_section => get_gauxc_section(xc_section)
      IF (ASSOCIATED(gauxc_section)) THEN
         CALL section_vals_val_get(gauxc_section, "NATIVE_GRID", l_val=native_grid)
         CALL section_vals_val_get(gauxc_section, "NATIVE_GRID_GAPW_COMPOSITE_REFERENCE", &
                                   l_val=use_composite_reference)
      END IF
      use_composite_reference = native_grid .AND. use_composite_reference

   END FUNCTION native_skala_gapw_composite_reference

! **************************************************************************************************
!> \brief Return true when the explicit atom-centered composite reference is requested.
!> \param xc_section ...
!> \return ...
! **************************************************************************************************
   FUNCTION native_skala_gapw_atom_composite_requested(xc_section) RESULT(use_atom_composite)
      TYPE(section_vals_type), INTENT(IN), POINTER       :: xc_section
      LOGICAL                                            :: use_atom_composite

      LOGICAL                                            :: native_grid
      TYPE(section_vals_type), POINTER                   :: gauxc_section

      use_atom_composite = .FALSE.
      native_grid = .FALSE.
      gauxc_section => get_gauxc_section(xc_section)
      IF (ASSOCIATED(gauxc_section)) THEN
         CALL section_vals_val_get(gauxc_section, "NATIVE_GRID", l_val=native_grid)
         CALL section_vals_val_get(gauxc_section, &
                                   "NATIVE_GRID_GAPW_ATOM_COMPOSITE_REFERENCE", &
                                   l_val=use_atom_composite)
      END IF
      use_atom_composite = native_grid .AND. use_atom_composite

   END FUNCTION native_skala_gapw_atom_composite_requested

! **************************************************************************************************
!> \brief Return true when native Skala uses atom-centered grids.
!> \param xc_section ...
!> \return ...
! **************************************************************************************************
   FUNCTION native_skala_uses_atom_composite_grid(xc_section) RESULT(use_atom_composite)
      TYPE(section_vals_type), INTENT(IN), POINTER       :: xc_section
      LOGICAL                                            :: use_atom_composite

      INTEGER                                            :: composite_grid
      TYPE(section_vals_type), POINTER                   :: gauxc_section

      use_atom_composite = .FALSE.
      composite_grid = skala_gapw_atom_composite_grid
      gauxc_section => get_gauxc_section(xc_section)
      IF (ASSOCIATED(gauxc_section)) THEN
         CALL section_vals_val_get(gauxc_section, "NATIVE_GRID_LAYOUT", &
                                   i_val=composite_grid)
      END IF
      SELECT CASE (composite_grid)
      CASE (skala_gapw_common_grid)
         CONTINUE
      CASE (skala_gapw_atom_composite_grid)
         use_atom_composite = .TRUE.
      CASE DEFAULT
         CALL cp_abort(__LOCATION__, "Unknown native-grid layout.")
      END SELECT

   END FUNCTION native_skala_uses_atom_composite_grid

! **************************************************************************************************
!> \brief Return true if the GAPW composite reference uses direct full-ORB collocation.
!> \param xc_section ...
!> \return ...
! **************************************************************************************************
   FUNCTION native_skala_gapw_composite_direct_ao(xc_section) RESULT(use_direct_ao)
      TYPE(section_vals_type), INTENT(IN), POINTER       :: xc_section
      LOGICAL                                            :: use_direct_ao

      TYPE(section_vals_type), POINTER                   :: gauxc_section

      use_direct_ao = .FALSE.
      gauxc_section => get_gauxc_section(xc_section)
      IF (ASSOCIATED(gauxc_section)) THEN
         CALL section_vals_val_get(gauxc_section, "NATIVE_GRID_GAPW_COMPOSITE_DIRECT_AO", &
                                   l_val=use_direct_ao)
      END IF
      use_direct_ao = native_skala_gapw_composite_reference(xc_section) .AND. use_direct_ao

   END FUNCTION native_skala_gapw_composite_direct_ao

! **************************************************************************************************
!> \brief Return true if the GAUXC subsection requests a model evaluation.
!> \param xc_section ...
!> \return ...
! **************************************************************************************************
   FUNCTION xc_section_uses_gauxc_model(xc_section) RESULT(uses_gauxc_model)
      TYPE(section_vals_type), INTENT(IN), POINTER       :: xc_section
      LOGICAL                                            :: uses_gauxc_model

      CHARACTER(len=default_path_length)                 :: model_key, model_name, xc_key, xc_name
      TYPE(section_vals_type), POINTER                   :: gauxc_section

      uses_gauxc_model = .FALSE.
      gauxc_section => get_gauxc_section(xc_section)
      IF (ASSOCIATED(gauxc_section)) THEN
         CALL section_vals_val_get(gauxc_section, "MODEL", c_val=model_name)
         CALL section_vals_val_get(gauxc_section, "FUNCTIONAL", c_val=xc_name)
         model_key = ADJUSTL(model_name)
         xc_key = ADJUSTL(xc_name)
         CALL uppercase(model_key)
         CALL uppercase(xc_key)
         uses_gauxc_model = (TRIM(model_key) /= "" .AND. TRIM(model_key) /= "NONE" .AND. &
                             TRIM(model_key) /= TRIM(xc_key))
      END IF

   END FUNCTION xc_section_uses_gauxc_model

! **************************************************************************************************
!> \brief Return true if the XC section requests a SKALA model (GauXC or native/FTorch).
!> \param xc_section ...
!> \return ...
! **************************************************************************************************
   FUNCTION xc_section_uses_skala_model(xc_section) RESULT(uses_skala_model)
      TYPE(section_vals_type), INTENT(IN), POINTER       :: xc_section
      LOGICAL                                            :: uses_skala_model

      CHARACTER(len=default_path_length)                 :: model_key, model_name
      TYPE(section_vals_type), POINTER                   :: gauxc_section

      uses_skala_model = .FALSE.
      gauxc_section => get_gauxc_section(xc_section)
      IF (.NOT. ASSOCIATED(gauxc_section)) RETURN

      CALL section_vals_val_get(gauxc_section, "MODEL", c_val=model_name)
      model_key = ADJUSTL(model_name)
      CALL uppercase(model_key)
      uses_skala_model = INDEX(TRIM(model_key), "SKALA") > 0

   END FUNCTION xc_section_uses_skala_model

! **************************************************************************************************
!> \brief Return the hard/soft GAPW one-center density partition for native SKALA.
!> \param xc_section ...
!> \return ...
! **************************************************************************************************
   FUNCTION native_skala_gapw_density_partition(xc_section) RESULT(partition)
      TYPE(section_vals_type), INTENT(IN), POINTER       :: xc_section
      INTEGER                                            :: partition

      TYPE(section_vals_type), POINTER                   :: gauxc_section

      partition = skala_gapw_density_partition_hard_minus_soft
      gauxc_section => get_gauxc_section(xc_section)
      IF (ASSOCIATED(gauxc_section)) THEN
         CALL section_vals_val_get(gauxc_section, "NATIVE_GRID_GAPW_DENSITY_PARTITION", &
                                   i_val=partition)
      END IF

      SELECT CASE (partition)
      CASE (skala_gapw_density_partition_hard_minus_soft, &
            skala_gapw_density_partition_hard_only, &
            skala_gapw_density_partition_soft_only, &
            skala_gapw_density_partition_none)
         CONTINUE
      CASE DEFAULT
         CALL cp_abort(__LOCATION__, &
                       "Unknown GAUXC%NATIVE_GRID_GAPW_DENSITY_PARTITION value.")
      END SELECT

   END FUNCTION native_skala_gapw_density_partition

! **************************************************************************************************
!> \brief Enforce the currently implemented native SKALA GPW input scope.
!> \param xc_section ...
! **************************************************************************************************
   SUBROUTINE ensure_native_skala_grid_scope(xc_section)
      TYPE(section_vals_type), INTENT(IN), POINTER       :: xc_section

      CHARACTER(len=default_path_length)                 :: model_key, model_name
      INTEGER                                            :: ifun, nfun
      LOGICAL                                            :: native_grid
      TYPE(section_vals_type), POINTER                   :: functionals, gauxc_section, xc_fun

      NULLIFY (gauxc_section)
      IF (.NOT. ASSOCIATED(xc_section)) THEN
         CPABORT("Native SKALA GPW requires an XC section")
      END IF

      functionals => section_vals_get_subs_vals(xc_section, "XC_FUNCTIONAL")
      IF (.NOT. ASSOCIATED(functionals)) THEN
         CPABORT("Native SKALA GPW requires an XC_FUNCTIONAL section")
      END IF

      nfun = 0
      ifun = 0
      DO
         ifun = ifun + 1
         xc_fun => section_vals_get_subs_vals2(functionals, i_section=ifun)
         IF (.NOT. ASSOCIATED(xc_fun)) EXIT
         nfun = nfun + 1
         IF (xc_fun%section%name == "GAUXC") gauxc_section => xc_fun
      END DO

      IF (.NOT. ASSOCIATED(gauxc_section)) THEN
         CPABORT("Native SKALA GPW requires an XC_FUNCTIONAL%GAUXC section")
      END IF
      IF (nfun /= 1) THEN
         CPABORT("Native SKALA GPW requires GAUXC to be the only XC functional")
      END IF

      CALL section_vals_val_get(gauxc_section, "NATIVE_GRID", l_val=native_grid)
      IF (.NOT. native_grid .AND. &
          .NOT. xc_section_uses_native_skala_evaluator(xc_section)) RETURN

      CALL section_vals_val_get(gauxc_section, "MODEL", c_val=model_name)
      model_key = ADJUSTL(model_name)
      CALL uppercase(model_key)
      IF (TRIM(model_key) == "NONE" .OR. TRIM(model_key) == "") THEN
         CPABORT("Native SKALA GPW requires GAUXC%MODEL SKALA or a TorchScript model path")
      END IF

   END SUBROUTINE ensure_native_skala_grid_scope

! **************************************************************************************************
!> \brief Evaluate SKALA energy and first derivatives on a CP2K GPW grid.
!> \param vxc_rho ...
!> \param vxc_tau ...
!> \param exc ...
!> \param rho_r ...
!> \param rho_g ...
!> \param tau ...
!> \param xc_section ...
!> \param weights ...
!> \param pw_pool ...
!> \param particle_set ...
!> \param cell ...
!> \param compute_virial ...
!> \param virial_xc ...
!> \param just_energy ...
!> \param atom_force ...
!> \param density_constraints Selected density subgradients of this KS environment.
!> \param rho_direction Optional fixed-geometry density direction on a differentiable model branch.
!> \param rho_direction_g Matching reciprocal density for the original spectral-gradient map.
!> \param tau_direction Matching kinetic-energy-density direction.
!> \param vxc_rho_response New density-potential response, with the same weighted adjoint as vxc_rho.
!> \param vxc_tau_response New kinetic-density-potential response. Selected nonsmooth FFT limits are excluded.
! **************************************************************************************************
   SUBROUTINE skala_gpw_eval(vxc_rho, vxc_tau, exc, rho_r, rho_g, tau, xc_section, &
                             weights, pw_pool, particle_set, cell, compute_virial, virial_xc, &
                             just_energy, atom_force, density_constraints, rho_direction, rho_direction_g, &
                             tau_direction, vxc_rho_response, vxc_tau_response)
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: vxc_rho, vxc_tau
      REAL(KIND=dp), INTENT(OUT)                         :: exc
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: rho_r
      TYPE(pw_c1d_gs_type), DIMENSION(:), POINTER        :: rho_g
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: tau
      TYPE(section_vals_type), POINTER                   :: xc_section
      TYPE(pw_r3d_rs_type), POINTER                      :: weights
      TYPE(pw_pool_type), POINTER                        :: pw_pool
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(cell_type), POINTER                           :: cell
      LOGICAL, INTENT(IN)                                :: compute_virial
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(OUT)        :: virial_xc
      LOGICAL, INTENT(IN), OPTIONAL                      :: just_energy
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT), &
         OPTIONAL                                        :: atom_force
      TYPE(skala_density_constraints_type), OPTIONAL, &
         POINTER                                         :: density_constraints
      TYPE(pw_r3d_rs_type), OPTIONAL, POINTER            :: rho_direction(:)
      TYPE(pw_c1d_gs_type), OPTIONAL, POINTER            :: rho_direction_g(:)
      TYPE(pw_r3d_rs_type), OPTIONAL, POINTER            :: tau_direction(:), vxc_rho_response(:), &
                                                            vxc_tau_response(:)

      CHARACTER(len=default_path_length)                 :: model_path
      INTEGER :: i, ipt, ispin, iw, j, k, native_grid_atom_chunk_max_rows, &
         native_grid_atom_partition, native_grid_atom_subchunks, native_grid_cuda_device, nspins, &
         omp_max_threads_restore, phase_handle, selected_cuda_device, xc_deriv_method_id, &
         xc_rho_smooth_id
      INTEGER, DIMENSION(2, 3)                           :: bo
      LOGICAL :: has_atom_chunk_work, have_atom_coord_grad, lsd, my_just_energy, &
         native_grid_atom_chunk_routing, native_grid_atom_chunks, native_grid_diagnostics, &
         native_grid_use_cuda, need_response, needs_atom_force, response_collapsed_rks, &
         use_atom_subchunks, use_density_constraints
      REAL(dp), ALLOCATABLE :: density_response(:, :), grad_response(:, :, :), kin_response(:, :), &
         model_density_direction(:, :), model_grad_direction(:, :, :), model_kin_direction(:, :)
      REAL(KIND=dp)                                      :: context_radius, density_contraction, &
                                                            tau_contraction
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: density_grad, kin_grad
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: grad_grad
      REAL(KIND=dp), DIMENSION(3, 3)                     :: virial_before
      TYPE(section_vals_type), POINTER                   :: gauxc_section
      TYPE(skala_gpw_feature_type)                       :: direction_features, features
      TYPE(torch_tensor_type) :: atom_coord_grad_t, atomic_grid_weight_grad_t, exc_tensor, &
         grid_coord_grad_t, grid_weight_grad_t, response_tensors(3)
      TYPE(xc_rho_cflags_type)                           :: needs
      TYPE(xc_rho_set_type)                              :: direction_rho_set, rho_set

      virial_xc = 0.0_dp
      exc = 0.0_dp
      my_just_energy = .FALSE.
      IF (PRESENT(just_energy)) my_just_energy = just_energy
      needs_atom_force = PRESENT(atom_force)
      IF (needs_atom_force) atom_force = 0.0_dp
      need_response = PRESENT(vxc_rho_response)
      CPASSERT(PRESENT(vxc_tau_response) .EQV. need_response)
      CPASSERT(PRESENT(rho_direction) .EQV. need_response)
      CPASSERT(PRESENT(rho_direction_g) .EQV. need_response)
      CPASSERT(PRESENT(tau_direction) .EQV. need_response)
      IF (need_response) THEN
         CPASSERT(.NOT. my_just_energy)
         CPASSERT(.NOT. compute_virial)
         CPASSERT(.NOT. needs_atom_force)
         CPASSERT(ASSOCIATED(rho_r) .AND. ASSOCIATED(rho_g) .AND. ASSOCIATED(tau))
         CPASSERT(ASSOCIATED(rho_direction))
         CPASSERT(ASSOCIATED(rho_direction_g))
         CPASSERT(ASSOCIATED(tau_direction))
         CPASSERT(.NOT. ASSOCIATED(vxc_rho_response))
         CPASSERT(.NOT. ASSOCIATED(vxc_tau_response))
         CPASSERT(SIZE(rho_direction) == SIZE(rho_r))
         CPASSERT(SIZE(rho_direction_g) == SIZE(rho_r))
         CPASSERT(SIZE(tau_direction) == SIZE(rho_r))
         DO i = 1, SIZE(rho_r)
            CPASSERT(ASSOCIATED(rho_direction(i)%pw_grid, rho_r(i)%pw_grid))
            CPASSERT(ASSOCIATED(rho_direction_g(i)%pw_grid, rho_g(i)%pw_grid))
            CPASSERT(ASSOCIATED(tau_direction(i)%pw_grid, tau(i)%pw_grid))
         END DO
         IF (PRESENT(density_constraints)) THEN
            IF (ASSOCIATED(density_constraints)) THEN
               i = skala_density_constraints_count(density_constraints, skala_fft_domain)
               CPASSERT(i == 0)
            END IF
         END IF
      END IF
      have_atom_coord_grad = .FALSE.
      omp_max_threads_restore = omp_get_max_threads()

      IF (compute_virial .AND. my_just_energy) THEN
         CALL cp_abort(__LOCATION__, &
                       "Native SKALA GPW stress/virial requires feature gradients.")
      END IF
      IF (.NOT. ASSOCIATED(rho_g)) THEN
         CALL cp_abort(__LOCATION__, &
                       "Native SKALA GPW requires the reciprocal-space density to form density gradients.")
      END IF
      IF (.NOT. ASSOCIATED(tau)) THEN
         CALL cp_abort(__LOCATION__, &
                       "Native SKALA GPW requires the kinetic-energy density.")
      END IF

      nspins = SIZE(rho_r)
      lsd = (nspins /= 1)
      use_density_constraints = .FALSE.
      IF (PRESENT(density_constraints)) THEN
         IF (ASSOCIATED(density_constraints)) THEN
            use_density_constraints = density_constraints_match_geometry( &
                                      density_constraints, rho_r, particle_set, cell)
         END IF
      END IF
      CALL get_skala_model_path(xc_section, model_path)
      gauxc_section => get_gauxc_section(xc_section)
      CALL section_vals_val_get(gauxc_section, "NATIVE_GRID_USE_CUDA", l_val=native_grid_use_cuda)
      CALL section_vals_val_get(gauxc_section, "NATIVE_GRID_CUDA_DEVICE", &
                                i_val=native_grid_cuda_device)
      CALL section_vals_val_get(gauxc_section, "NATIVE_GRID_ATOM_CHUNKS", &
                                l_val=native_grid_atom_chunks)
      CALL section_vals_val_get(gauxc_section, "NATIVE_GRID_ATOM_CHUNK_ROUTING", &
                                l_val=native_grid_atom_chunk_routing)
      CALL section_vals_val_get(gauxc_section, "NATIVE_GRID_ATOM_CHUNK_MAX_ROWS", &
                                i_val=native_grid_atom_chunk_max_rows)
      CALL section_vals_val_get(gauxc_section, "NATIVE_GRID_ATOM_PARTITION", &
                                i_val=native_grid_atom_partition)
      context_radius = section_get_rval(gauxc_section, "NATIVE_GRID_CONTEXT_RADIUS")
      SELECT CASE (native_grid_atom_partition)
      CASE (1)
         native_grid_atom_partition = skala_gpw_atom_partition_hard
      CASE (2)
         native_grid_atom_partition = skala_gpw_atom_partition_smooth
      CASE DEFAULT
         CALL cp_abort(__LOCATION__, &
                       "Unknown GAUXC%NATIVE_GRID_ATOM_PARTITION value.")
      END SELECT
      native_grid_atom_chunk_routing = native_grid_atom_chunk_routing .OR. native_grid_atom_chunks
      native_grid_atom_chunks = native_grid_atom_chunks .OR. native_grid_atom_chunk_routing
      IF (native_grid_atom_chunk_max_rows < -1) THEN
         CALL cp_abort(__LOCATION__, &
                       "GAUXC%NATIVE_GRID_ATOM_CHUNK_MAX_ROWS must be -1, zero, or positive.")
      END IF
      IF (needs_atom_force .OR. compute_virial) THEN
         IF (native_grid_atom_partition == skala_gpw_atom_partition_hard) THEN
            native_grid_atom_partition = skala_gpw_atom_partition_smooth
         END IF
         native_grid_atom_chunk_routing = .FALSE.
         native_grid_atom_chunks = .FALSE.
      END IF
      IF (use_density_constraints .AND. .NOT. my_just_energy) THEN
         ! Both limiting covectors refer to the complete original feature dictionary.
         native_grid_atom_chunk_routing = .FALSE.
         native_grid_atom_chunks = .FALSE.
      END IF
      ! The portable SKALA export used by the regtests builds ragged-index tensors on CPU.
      CALL torch_use_cuda(native_grid_use_cuda)
      selected_cuda_device = configure_native_grid_cuda( &
                             native_grid_use_cuda, native_grid_cuda_device, rho_r(1)%pw_grid%para%group)
      CALL ensure_model_loaded(model_path, selected_cuda_device, xc_section)

      IF (lsd) THEN
         needs%rho_spin = .TRUE.
         needs%drho_spin = .TRUE.
         needs%tau_spin = .TRUE.
      ELSE
         needs%rho = .TRUE.
         needs%drho = .TRUE.
         needs%tau = .TRUE.
      END IF

      CALL section_vals_val_get(xc_section, "XC_GRID%XC_DERIV", i_val=xc_deriv_method_id)
      CALL section_vals_val_get(xc_section, "XC_GRID%XC_SMOOTH_RHO", i_val=xc_rho_smooth_id)

      CALL xc_rho_set_create(rho_set, &
                             rho_r(1)%pw_grid%bounds_local, &
                             rho_cutoff=section_get_rval(xc_section, "density_cutoff"), &
                             drho_cutoff=section_get_rval(xc_section, "gradient_cutoff"), &
                             tau_cutoff=section_get_rval(xc_section, "tau_cutoff"))
      CALL xc_rho_set_update(rho_set, rho_r, rho_g, tau, needs, &
                             xc_deriv_method_id, xc_rho_smooth_id, pw_pool)

      IF (PRESENT(density_constraints)) THEN
         IF (ASSOCIATED(density_constraints)) THEN
            CALL record_fft_descriptors(density_constraints, rho_set, rho_r)
            IF (skala_density_samples_only(density_constraints)) THEN
               CPASSERT(my_just_energy)
               CPASSERT(.NOT. compute_virial)
               CPASSERT(.NOT. needs_atom_force)
               CALL xc_rho_set_release(rho_set, pw_pool=pw_pool)
               CALL torch_use_cuda(.TRUE.)
               CALL omp_set_num_threads(omp_max_threads_restore)
               RETURN
            END IF
         END IF
      END IF

      IF (need_response) THEN
         ! Capture independent directional storage before the base feature build resets
         ! the shared dynamic tensors. Both maps use the same fixed geometry and layout.
         CALL xc_rho_set_create(direction_rho_set, rho_r(1)%pw_grid%bounds_local, &
                                rho_cutoff=section_get_rval(xc_section, "density_cutoff"), &
                                drho_cutoff=section_get_rval(xc_section, "gradient_cutoff"), &
                                tau_cutoff=section_get_rval(xc_section, "tau_cutoff"))
         CALL xc_rho_set_update(direction_rho_set, rho_direction, rho_direction_g, tau_direction, needs, &
                                xc_deriv_method_id, xc_rho_smooth_id, pw_pool)
         CALL skala_gpw_feature_build(direction_features, direction_rho_set, rho_direction, particle_set, cell, &
                                      requires_grad=.FALSE., weights=weights, use_atom_chunks=native_grid_atom_chunks, &
         route_atom_chunks=native_grid_atom_chunk_routing, atom_partition=native_grid_atom_partition, context_radius=context_radius)
         response_collapsed_rks = direction_features%uses_collapsed_rks_dynamic
         IF (direction_features%uses_atom_chunks) THEN
            IF (direction_features%chunk_feature_count > 0) THEN
               ALLOCATE (model_density_direction, SOURCE=direction_features%chunk_density)
               ALLOCATE (model_grad_direction, SOURCE=direction_features%chunk_grad)
               ALLOCATE (model_kin_direction, SOURCE=direction_features%chunk_kin)
            ELSE
               i = 2
               IF (response_collapsed_rks) i = 1
               ALLOCATE (model_density_direction(0, i), model_grad_direction(0, 3, i), model_kin_direction(0, i))
            END IF
         ELSE
            ALLOCATE (model_density_direction, SOURCE=direction_features%density)
            ALLOCATE (model_grad_direction, SOURCE=direction_features%grad)
            ALLOCATE (model_kin_direction, SOURCE=direction_features%kin)
         END IF
         CALL skala_gpw_feature_release(direction_features)
         CALL xc_rho_set_release(direction_rho_set, pw_pool=pw_pool)
      END IF
      CALL skala_gpw_feature_build(features, rho_set, rho_r, particle_set, cell, &
                                   requires_grad=(.NOT. my_just_energy), weights=weights, &
                                   requires_coordinate_grad=(needs_atom_force .OR. compute_virial), &
                                   requires_stress_grad=compute_virial, &
                                   use_atom_chunks=native_grid_atom_chunks, &
                                   route_atom_chunks=native_grid_atom_chunk_routing, &
                                   atom_partition=native_grid_atom_partition, context_radius=context_radius)
      IF (need_response) THEN
         CPASSERT(response_collapsed_rks .EQV. features%uses_collapsed_rks_dynamic)
      END IF
      CALL section_vals_val_get(gauxc_section, "NATIVE_GRID_DIAGNOSTICS", l_val=native_grid_diagnostics)
      IF (native_grid_diagnostics) THEN
         CALL print_native_grid_diagnostics(features, rho_r(1)%pw_grid%para%group%mepos == 0)
      END IF

      IF (features%uses_atom_chunks .AND. native_grid_atom_chunk_max_rows == -1) THEN
         native_grid_atom_chunk_max_rows = auto_atom_chunk_max_rows(features, &
                                                                    rho_r(1)%pw_grid%para%group)
      END IF
      IF (native_grid_diagnostics .AND. features%uses_atom_chunks .AND. &
          rho_r(1)%pw_grid%para%group%mepos == 0) THEN
         iw = cp_logger_get_default_io_unit()
         IF (iw > 0) THEN
            WRITE (UNIT=iw, FMT="(T2,A,1X,I0)") &
               "SKALA_GPW| Native grid atom chunk max rows", native_grid_atom_chunk_max_rows
         END IF
      END IF
      native_grid_atom_subchunks = 1
      IF (features%uses_atom_chunks .AND. native_grid_atom_chunk_max_rows > 0) THEN
         native_grid_atom_subchunks = skala_gpw_atom_subchunk_count(native_grid_atom_chunk_max_rows)
         CALL rho_r(1)%pw_grid%para%group%max(native_grid_atom_subchunks)
      END IF
      use_atom_subchunks = features%uses_atom_chunks .AND. native_grid_atom_subchunks > 1
      has_atom_chunk_work = .NOT. features%uses_atom_chunks .OR. features%chunk_feature_count > 0
      exc = 0.0_dp
      IF (use_atom_subchunks) THEN
         IF (need_response) THEN
            CALL evaluate_atom_subchunks(features, rho_r(1)%pw_grid%para%group, native_grid_atom_chunk_max_rows, &
                                         compute_grads=.TRUE., exc=exc, density_grad=density_grad, grad_grad=grad_grad, &
                                         kin_grad=kin_grad, collapse_spin_grads=(nspins == 1), &
                                         density_direction=model_density_direction, grad_direction=model_grad_direction, &
                                         kin_direction=model_kin_direction, density_response=density_response, &
                                         grad_response=grad_response, kin_response=kin_response)
         ELSE
            CALL evaluate_atom_subchunks(features, rho_r(1)%pw_grid%para%group, &
                                         native_grid_atom_chunk_max_rows, &
                                         compute_grads=(.NOT. my_just_energy), exc=exc, &
                                         density_grad=density_grad, grad_grad=grad_grad, &
                                         kin_grad=kin_grad, collapse_spin_grads=(nspins == 1))
         END IF
      ELSE IF (has_atom_chunk_work) THEN
         CALL skala_torch_model_get_exc(cached_model, features%inputs, &
                                        features%grid_weights_t, exc_tensor, exc)
      END IF
      IF (features%uses_atom_chunks) CALL rho_r(1)%pw_grid%para%group%sum(exc)

      IF (.NOT. my_just_energy) THEN
         IF (.NOT. use_atom_subchunks) THEN
            IF (need_response) THEN
               IF (has_atom_chunk_work) THEN
                  CALL feature_hessian_product(features, exc_tensor, model_density_direction, model_grad_direction, &
                                               model_kin_direction, response_tensors)
               END IF
               IF (features%uses_atom_chunks) THEN
                  CALL fetch_and_gather_atom_chunk_grads(features, rho_r(1)%pw_grid%para%group, &
                                                         density_response, grad_response, kin_response, response_tensors)
               ELSE
                  CALL fetch_local_feature_grads(features, density_response, grad_response, kin_response, response_tensors)
               END IF
               IF (has_atom_chunk_work) THEN
                  DO i = 1, 3
                     CALL torch_tensor_release(response_tensors(i))
                  END DO
               END IF
            END IF
            IF (has_atom_chunk_work) THEN
               CALL timeset("skala_gpw_backward", phase_handle)
               CALL torch_tensor_backward_scalar(exc_tensor)
               CALL timestop(phase_handle)

               IF (compute_virial) THEN
                  IF (native_grid_diagnostics) virial_before = virial_xc
                  CALL build_weight_virial(virial_xc, features, exc, grid_weight_grad_t, &
                                           atomic_grid_weight_grad_t, &
                                           rho_r(1)%pw_grid%para%group%mepos == 0, &
                                           native_grid_diagnostics)
                  IF (native_grid_diagnostics) THEN
                     CALL print_virial_delta("weight-residual", virial_xc - virial_before, &
                                             rho_r(1)%pw_grid%para%group%mepos == 0)
                  END IF
               END IF
            END IF

            CALL timeset("skala_gpw_grad_fetch", phase_handle)
            IF (features%uses_atom_chunks) THEN
               CALL fetch_and_gather_atom_chunk_grads(features, rho_r(1)%pw_grid%para%group, &
                                                      density_grad, grad_grad, kin_grad)
            ELSE
               CALL fetch_local_feature_grads(features, density_grad, grad_grad, kin_grad)
            END IF
            CALL timestop(phase_handle)
         END IF
         IF (use_density_constraints) THEN
            CALL apply_descriptor_constraints(density_constraints, features, rho_r, density_grad, kin_grad)
         END IF
         IF (needs_atom_force) THEN
            CALL add_explicit_coordinate_force(atom_force, features, atom_coord_grad_t, &
                                               rho_r(1)%pw_grid%para%group%mepos == 0)
            IF (features%atom_partition == skala_gpw_atom_partition_smooth) THEN
               CALL add_smooth_partition_force(atom_force, features, particle_set, cell, rho_r, &
                                               grid_weight_grad_t, atomic_grid_weight_grad_t)
            END IF
            have_atom_coord_grad = .TRUE.
         END IF

         CALL timeset("skala_gpw_vxc_unpack", phase_handle)
         IF (native_grid_diagnostics) THEN
            bo = rho_r(1)%pw_grid%bounds_local
            density_contraction = 0.0_dp
            tau_contraction = 0.0_dp
            ipt = 0
            DO k = bo(1, 3), bo(2, 3)
               DO j = bo(1, 2), bo(2, 2)
                  DO i = bo(1, 1), bo(2, 1)
                     ipt = ipt + 1
                     IF (nspins == 1) THEN
                        density_contraction = density_contraction + rho_r(1)%array(i, j, k)* &
                                              0.5_dp*(density_grad(ipt, 1) + density_grad(ipt, 2))
                        tau_contraction = tau_contraction + tau(1)%array(i, j, k)* &
                                          0.5_dp*(kin_grad(ipt, 1) + kin_grad(ipt, 2))
                     ELSE
                        DO ispin = 1, nspins
                           density_contraction = density_contraction + &
                                                 rho_r(ispin)%array(i, j, k)*density_grad(ipt, ispin)
                           tau_contraction = tau_contraction + &
                                             tau(ispin)%array(i, j, k)*kin_grad(ipt, ispin)
                        END DO
                     END IF
                  END DO
               END DO
            END DO
            CALL rho_r(1)%pw_grid%para%group%sum(density_contraction)
            CALL rho_r(1)%pw_grid%para%group%sum(tau_contraction)
            IF (rho_r(1)%pw_grid%para%group%mepos == 0) THEN
               iw = cp_logger_get_default_io_unit()
               IF (iw > 0) THEN
                  WRITE (iw, "(T2,A,1X,ES20.10)") &
                     "SKALA_GPW| XC density-gradient contraction", density_contraction
                  WRITE (iw, "(T2,A,1X,ES20.10)") &
                     "SKALA_GPW| XC kinetic-gradient contraction", tau_contraction
               END IF
            END IF
         END IF
         IF (compute_virial) THEN
            IF (native_grid_diagnostics) virial_before = virial_xc
            CALL build_virial_from_feature_grads(virial_xc, rho_set, rho_r, grad_grad)
            IF (native_grid_diagnostics) THEN
               CALL print_virial_delta("feature-gradient", virial_xc - virial_before, &
                                       rho_r(1)%pw_grid%para%group%mepos == 0)
               virial_before = virial_xc
            END IF
            IF (.NOT. have_atom_coord_grad) THEN
               CALL torch_tensor_grad(features%coarse_0_atomic_coords_t, atom_coord_grad_t)
               have_atom_coord_grad = .TRUE.
            END IF
            CALL build_static_coordinate_virial(virial_xc, features, atom_coord_grad_t, &
                                                grid_coord_grad_t, &
                                                rho_r(1)%pw_grid%para%group%mepos == 0, &
                                                native_grid_diagnostics)
            IF (native_grid_diagnostics) THEN
               CALL print_virial_delta("static-coordinates", virial_xc - virial_before, &
                                       rho_r(1)%pw_grid%para%group%mepos == 0)
               virial_before = virial_xc
            END IF
            IF (features%atom_partition == skala_gpw_atom_partition_smooth) THEN
               CALL build_smooth_partition_virial(virial_xc, features, particle_set, cell, rho_r, &
                                                  grid_weight_grad_t, atomic_grid_weight_grad_t)
               IF (native_grid_diagnostics) THEN
                  CALL print_virial_delta("smooth-partition", virial_xc - virial_before, &
                                          rho_r(1)%pw_grid%para%group%mepos == 0)
                  virial_before = virial_xc
               END IF
            END IF
         END IF
         CALL omp_set_num_threads(omp_max_threads_restore)
         CALL build_vxc_from_feature_grads(vxc_rho, vxc_tau, rho_r, pw_pool, &
                                           density_grad, grad_grad, kin_grad, &
                                           xc_deriv_method_id)
         IF (need_response) THEN
            CALL build_vxc_from_feature_grads(vxc_rho_response, vxc_tau_response, rho_r, pw_pool, &
                                              density_response, grad_response, kin_response, xc_deriv_method_id)
         END IF
         CALL timestop(phase_handle)

         CALL timeset("skala_gpw_grad_release", phase_handle)
         DEALLOCATE (density_grad, grad_grad, kin_grad)
         IF (have_atom_coord_grad) CALL torch_tensor_release(atom_coord_grad_t)
         CALL timestop(phase_handle)
      END IF

      CALL timeset("skala_gpw_cleanup", phase_handle)
      IF (.NOT. use_atom_subchunks .AND. has_atom_chunk_work) CALL torch_tensor_release(exc_tensor)
      CALL skala_gpw_feature_release(features)
      CALL xc_rho_set_release(rho_set, pw_pool=pw_pool)
      CALL torch_use_cuda(.TRUE.)
      CALL omp_set_num_threads(omp_max_threads_restore)
      CALL timestop(phase_handle)

   END SUBROUTINE skala_gpw_eval

! **************************************************************************************************
!> \brief Evaluate the derivative of native SKALA XC energy with respect to CP2K's external
!>        real-space integration-weight multiplier.
!> \param weight_deriv_r ...
!> \param rho_r ...
!> \param rho_g ...
!> \param tau ...
!> \param xc_section ...
!> \param weights ...
!> \param pw_pool ...
!> \param particle_set ...
!> \param cell ...
! **************************************************************************************************
   SUBROUTINE skala_gpw_weight_derivative(weight_deriv_r, rho_r, rho_g, tau, xc_section, &
                                          weights, pw_pool, particle_set, cell)
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: weight_deriv_r
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: rho_r
      TYPE(pw_c1d_gs_type), DIMENSION(:), POINTER        :: rho_g
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: tau
      TYPE(section_vals_type), POINTER                   :: xc_section
      TYPE(pw_r3d_rs_type), POINTER                      :: weights
      TYPE(pw_pool_type), POINTER                        :: pw_pool
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(cell_type), POINTER                           :: cell

      CHARACTER(len=default_path_length)                 :: model_path
      INTEGER :: feature_begin, feature_end, feature_pos, i, iatom, j, k, local_row, &
         native_grid_atom_partition, native_grid_cuda_device, natom, nspins, &
         omp_max_threads_restore, row, selected_cuda_device, xc_deriv_method_id, xc_rho_smooth_id
      INTEGER, DIMENSION(2, 3)                           :: bo
      LOGICAL                                            :: lsd, native_grid_use_cuda
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: included
      REAL(KIND=dp)                                      :: context_radius, exc, local_derivative
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: partition_weights
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: atom_coords_pbc
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: dweights_datom, dweights_dstrain
      REAL(KIND=dp), DIMENSION(3)                        :: grid_point
      REAL(KIND=dp), DIMENSION(:), POINTER               :: grid_weight_grad
      TYPE(section_vals_type), POINTER                   :: gauxc_section
      TYPE(skala_gpw_feature_type)                       :: features
      TYPE(torch_tensor_type)                            :: exc_tensor, grid_weight_grad_t
      TYPE(xc_rho_cflags_type)                           :: needs
      TYPE(xc_rho_set_type)                              :: rho_set

      CPASSERT(ASSOCIATED(rho_r))
      CPASSERT(ASSOCIATED(rho_g))
      CPASSERT(ASSOCIATED(tau))
      CALL pw_zero(weight_deriv_r)
      omp_max_threads_restore = omp_get_max_threads()

      nspins = SIZE(rho_r)
      lsd = (nspins /= 1)
      CALL get_skala_model_path(xc_section, model_path)
      gauxc_section => get_gauxc_section(xc_section)
      CALL section_vals_val_get(gauxc_section, "NATIVE_GRID_USE_CUDA", l_val=native_grid_use_cuda)
      CALL section_vals_val_get(gauxc_section, "NATIVE_GRID_CUDA_DEVICE", &
                                i_val=native_grid_cuda_device)
      CALL section_vals_val_get(gauxc_section, "NATIVE_GRID_ATOM_PARTITION", &
                                i_val=native_grid_atom_partition)
      context_radius = section_get_rval(gauxc_section, "NATIVE_GRID_CONTEXT_RADIUS")
      SELECT CASE (native_grid_atom_partition)
      CASE (1)
         native_grid_atom_partition = skala_gpw_atom_partition_hard
      CASE (2)
         native_grid_atom_partition = skala_gpw_atom_partition_smooth
      CASE DEFAULT
         CALL cp_abort(__LOCATION__, &
                       "Unknown GAUXC%NATIVE_GRID_ATOM_PARTITION value.")
      END SELECT
      IF (native_grid_atom_partition == skala_gpw_atom_partition_hard) THEN
         native_grid_atom_partition = skala_gpw_atom_partition_smooth
      END IF

      CALL torch_use_cuda(native_grid_use_cuda)
      selected_cuda_device = configure_native_grid_cuda( &
                             native_grid_use_cuda, native_grid_cuda_device, rho_r(1)%pw_grid%para%group)
      CALL ensure_model_loaded(model_path, selected_cuda_device, xc_section)

      IF (lsd) THEN
         needs%rho_spin = .TRUE.
         needs%drho_spin = .TRUE.
         needs%tau_spin = .TRUE.
      ELSE
         needs%rho = .TRUE.
         needs%drho = .TRUE.
         needs%tau = .TRUE.
      END IF

      CALL section_vals_val_get(xc_section, "XC_GRID%XC_DERIV", i_val=xc_deriv_method_id)
      CALL section_vals_val_get(xc_section, "XC_GRID%XC_SMOOTH_RHO", i_val=xc_rho_smooth_id)

      CALL xc_rho_set_create(rho_set, &
                             rho_r(1)%pw_grid%bounds_local, &
                             rho_cutoff=section_get_rval(xc_section, "density_cutoff"), &
                             drho_cutoff=section_get_rval(xc_section, "gradient_cutoff"), &
                             tau_cutoff=section_get_rval(xc_section, "tau_cutoff"))
      CALL xc_rho_set_update(rho_set, rho_r, rho_g, tau, needs, &
                             xc_deriv_method_id, xc_rho_smooth_id, pw_pool)

      CALL skala_gpw_feature_build(features, rho_set, rho_r, particle_set, cell, &
                                   requires_grad=.FALSE., weights=weights, &
                                   requires_coordinate_grad=.FALSE., &
                                   requires_stress_grad=.TRUE., &
                                   use_atom_chunks=.FALSE., route_atom_chunks=.FALSE., &
                                   atom_partition=native_grid_atom_partition, context_radius=context_radius)
      CALL skala_torch_model_get_exc(cached_model, features%inputs, features%grid_weights_t, &
                                     exc_tensor, exc)
      CALL torch_tensor_backward_scalar(exc_tensor)
      NULLIFY (grid_weight_grad)
      CALL torch_tensor_grad(features%grid_weights_t, grid_weight_grad_t)
      CALL torch_tensor_data_ptr(grid_weight_grad_t, grid_weight_grad)

      natom = SIZE(particle_set)
      IF (native_grid_atom_partition == skala_gpw_atom_partition_smooth .AND. .NOT. features%periodic_context) THEN
         ALLOCATE (atom_coords_pbc(3, natom), included(natom), partition_weights(natom), &
                   dweights_datom(3, natom, natom), dweights_dstrain(3, 3, natom))
         DO iatom = 1, natom
            atom_coords_pbc(:, iatom) = pbc(particle_set(iatom)%r, cell, positive_range=.TRUE.)
         END DO
      END IF

      IF (features%periodic_context) THEN
         CALL skala_gpw_periodic_weight_derivative(features, particle_set, cell, rho_r(1)%pw_grid, &
                                                   grid_weight_grad, weight_deriv_r)
      ELSE
         bo = rho_r(1)%pw_grid%bounds_local
         local_row = 0
         DO k = bo(1, 3), bo(2, 3)
            DO j = bo(1, 2), bo(2, 2)
               DO i = bo(1, 1), bo(2, 1)
                  local_row = local_row + 1
                  feature_begin = features%local_feature_offsets(local_row)
                  feature_end = features%local_feature_offsets(local_row + 1) - 1
                  local_derivative = 0.0_dp
                  IF (native_grid_atom_partition == skala_gpw_atom_partition_hard) THEN
                     CPASSERT(feature_end == feature_begin)
                     row = features%local_feature_rows(feature_begin)
                     local_derivative = grid_weight_grad(row)
                  ELSE
                     grid_point = native_grid_coordinate(rho_r(1)%pw_grid, [i, j, k])
                     CALL skala_gpw_smooth_partition_derivatives( &
                        grid_point, atom_coords_pbc, cell, partition_weights, included, &
                        dweights_datom, dweights_dstrain)
                     CPASSERT(feature_end - feature_begin + 1 == COUNT(included))
                     feature_pos = feature_begin
                     DO iatom = 1, natom
                        IF (.NOT. included(iatom)) CYCLE
                        row = features%local_feature_rows(feature_pos)
                        local_derivative = local_derivative + &
                                           partition_weights(iatom)*grid_weight_grad(row)
                        feature_pos = feature_pos + 1
                     END DO
                     CPASSERT(feature_pos == feature_end + 1)
                  END IF
                  weight_deriv_r%array(i, j, k) = local_derivative
               END DO
            END DO
         END DO
         CPASSERT(local_row == features%nflat_local)
      END IF

      IF (ALLOCATED(atom_coords_pbc)) THEN
         DEALLOCATE (atom_coords_pbc, dweights_datom, dweights_dstrain, included, &
                     partition_weights)
      END IF
      CALL torch_tensor_release(grid_weight_grad_t)
      CALL torch_tensor_release(exc_tensor)
      CALL skala_gpw_feature_release(features)
      CALL xc_rho_set_release(rho_set, pw_pool=pw_pool)
      CALL torch_use_cuda(.TRUE.)
      CALL omp_set_num_threads(omp_max_threads_restore)

   END SUBROUTINE skala_gpw_weight_derivative

! **************************************************************************************************
!> \brief Evaluate SKALA on a GAPW one-center atomic grid.
!> \param xc_section ...
!> \param grid_atom ...
!> \param group ...
!> \param atom_coord ...
!> \param rho ...
!> \param drho ...
!> \param tau ...
!> \param weights ...
!> \param lsd ...
!> \param nspins ...
!> \param na ...
!> \param nr ...
!> \param exc ...
!> \param vxc ...
!> \param vxg ...
!> \param vtau ...
!> \param energy_only ...
!> \param atom_force ...
!> \param atom_virial ...
! **************************************************************************************************
   SUBROUTINE skala_gapw_atom_vxc_of_r(xc_section, grid_atom, group, atom_coord, &
                                       rho, drho, tau, weights, lsd, nspins, na, nr, &
                                       exc, vxc, vxg, vtau, energy_only, atom_force, atom_virial)
      TYPE(section_vals_type), POINTER                   :: xc_section
      TYPE(grid_atom_type), POINTER                      :: grid_atom

      CLASS(mp_comm_type), INTENT(IN)                    :: group
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: atom_coord
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: rho, tau, vxc, vtau
      REAL(KIND=dp), DIMENSION(:, :, :, :), POINTER      :: drho, vxg
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: weights
      LOGICAL, INTENT(IN)                                :: lsd
      INTEGER, INTENT(IN)                                :: nspins, na, nr
      REAL(KIND=dp), INTENT(OUT)                         :: exc
      LOGICAL, INTENT(IN), OPTIONAL                      :: energy_only
      REAL(KIND=dp), DIMENSION(3), INTENT(OUT), &
         OPTIONAL                                        :: atom_force
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(OUT), &
         OPTIONAL                                        :: atom_virial

      CHARACTER(len=default_path_length)                 :: model_path
      INTEGER                                            :: ia, idir, ir, native_grid_cuda_device, &
                                                            jdir, nflat, omp_max_threads_restore, row, &
                                                            selected_cuda_device
      INTEGER(KIND=int_8), ALLOCATABLE, DIMENSION(:)     :: atomic_grid_sizes
      INTEGER(KIND=int_8), ALLOCATABLE, DIMENSION(:, :)  :: atomic_grid_size_bound_shape
      LOGICAL                                            :: need_coord_grad, my_energy_only, native_grid_use_cuda
      REAL(KIND=dp)                                      :: tmp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: atomic_grid_weights, grid_weights
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: coarse_0_atomic_coords, density, &
                                                            grid_coords, kin
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: grad
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: atom_coord_grad, density_grad, &
                                                            grid_coord_grad, kin_grad
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: grad_grad
      TYPE(section_vals_type), POINTER                   :: gauxc_section
      TYPE(torch_dict_type)                              :: inputs
      TYPE(torch_tensor_type)                            :: atomic_grid_size_bound_shape_t, &
                                                            atomic_grid_sizes_t, &
                                                            atomic_grid_weights_t, &
                                                            atom_coord_grad_t, &
                                                            coarse_0_atomic_coords_t, density_t, &
                                                            density_grad_t, exc_tensor, grad_t, &
                                                            grad_grad_t, grid_coord_grad_t, &
                                                            grid_coords_t, grid_weights_t, kin_t, &
                                                            kin_grad_t

      CPASSERT(ASSOCIATED(xc_section))
      CPASSERT(ASSOCIATED(grid_atom))
      CPASSERT(ASSOCIATED(rho))
      CPASSERT(ASSOCIATED(drho))
      CPASSERT(ASSOCIATED(tau))
      omp_max_threads_restore = omp_get_max_threads()

      my_energy_only = .FALSE.
      IF (PRESENT(energy_only)) my_energy_only = energy_only
      need_coord_grad = PRESENT(atom_force) .OR. PRESENT(atom_virial)
      exc = 0.0_dp
      IF (PRESENT(atom_force)) atom_force = 0.0_dp
      IF (PRESENT(atom_virial)) atom_virial = 0.0_dp
      IF (.NOT. my_energy_only) THEN
         vxc = 0.0_dp
         vxg = 0.0_dp
         vtau = 0.0_dp
      END IF

      CALL get_skala_model_path(xc_section, model_path)
      gauxc_section => get_gauxc_section(xc_section)
      CPASSERT(ASSOCIATED(gauxc_section))
      CALL section_vals_val_get(gauxc_section, "NATIVE_GRID_USE_CUDA", l_val=native_grid_use_cuda)
      CALL section_vals_val_get(gauxc_section, "NATIVE_GRID_CUDA_DEVICE", &
                                i_val=native_grid_cuda_device)
      CALL torch_use_cuda(native_grid_use_cuda)
      selected_cuda_device = configure_native_grid_cuda( &
                             native_grid_use_cuda, native_grid_cuda_device, group)
      CALL ensure_model_loaded(model_path, selected_cuda_device, xc_section)

      nflat = na*nr
      ALLOCATE (density(nflat, 2), grad(nflat, 3, 2), kin(nflat, 2), &
                grid_coords(3, nflat), grid_weights(nflat), &
                atomic_grid_weights(nflat), atomic_grid_sizes(1), &
                coarse_0_atomic_coords(3, 1), atomic_grid_size_bound_shape(0, nflat))
      density = 0.0_dp
      grad = 0.0_dp
      kin = 0.0_dp
      grid_coords = 0.0_dp
      grid_weights = 0.0_dp
      atomic_grid_weights = 0.0_dp
      atomic_grid_sizes(1) = INT(nflat, KIND=int_8)
      atomic_grid_size_bound_shape = 0_int_8
      coarse_0_atomic_coords(:, 1) = atom_coord

      row = 0
      DO ir = 1, nr
         DO ia = 1, na
            row = row + 1
            grid_coords(1, row) = atom_coord(1) + grid_atom%rad(ir)* &
                                  grid_atom%sin_pol(ia)*grid_atom%cos_azi(ia)
            grid_coords(2, row) = atom_coord(2) + grid_atom%rad(ir)* &
                                  grid_atom%sin_pol(ia)*grid_atom%sin_azi(ia)
            grid_coords(3, row) = atom_coord(3) + grid_atom%rad(ir)*grid_atom%cos_pol(ia)
            grid_weights(row) = weights(ia, ir)
            atomic_grid_weights(row) = weights(ia, ir)
            IF (nspins == 1) THEN
               density(row, :) = 0.5_dp*rho(ia, ir, 1)
               DO idir = 1, 3
                  grad(row, idir, :) = 0.5_dp*drho(idir, ia, ir, 1)
               END DO
               kin(row, :) = 0.5_dp*tau(ia, ir, 1)
            ELSE
               density(row, :) = rho(ia, ir, 1:2)
               DO idir = 1, 3
                  grad(row, idir, :) = drho(idir, ia, ir, 1:2)
               END DO
               kin(row, :) = tau(ia, ir, 1:2)
            END IF
         END DO
      END DO

      CALL torch_tensor_from_array(grid_coords_t, grid_coords)
      CALL torch_tensor_to_device_leaf(grid_coords_t, need_coord_grad)
      CALL torch_tensor_from_array(grid_weights_t, grid_weights)
      CALL torch_tensor_to_device_leaf(grid_weights_t, .FALSE.)
      CALL torch_tensor_from_array(atomic_grid_weights_t, atomic_grid_weights)
      CALL torch_tensor_to_device_leaf(atomic_grid_weights_t, .FALSE.)
      CALL torch_tensor_from_array(atomic_grid_sizes_t, atomic_grid_sizes)
      CALL torch_tensor_to_device_leaf(atomic_grid_sizes_t, .FALSE.)
      CALL torch_tensor_from_array(atomic_grid_size_bound_shape_t, &
                                   atomic_grid_size_bound_shape)
      CALL torch_tensor_to_device_leaf(atomic_grid_size_bound_shape_t, .FALSE.)
      CALL torch_tensor_from_array(coarse_0_atomic_coords_t, coarse_0_atomic_coords)
      CALL torch_tensor_to_device_leaf(coarse_0_atomic_coords_t, need_coord_grad)
      CALL torch_tensor_from_array(density_t, density)
      CALL torch_tensor_to_device_leaf(density_t,.NOT. my_energy_only)
      CALL torch_tensor_from_array(grad_t, grad)
      CALL torch_tensor_to_device_leaf(grad_t,.NOT. my_energy_only)
      CALL torch_tensor_from_array(kin_t, kin)
      CALL torch_tensor_to_device_leaf(kin_t,.NOT. my_energy_only)

      CALL torch_dict_create(inputs)
      CALL torch_dict_insert(inputs, "grid_coords", grid_coords_t)
      CALL torch_dict_insert(inputs, "grid_weights", grid_weights_t)
      CALL torch_dict_insert(inputs, "atomic_grid_weights", atomic_grid_weights_t)
      CALL torch_dict_insert(inputs, "atomic_grid_sizes", atomic_grid_sizes_t)
      CALL torch_dict_insert(inputs, "atomic_grid_size_bound_shape", &
                             atomic_grid_size_bound_shape_t)
      CALL torch_dict_insert(inputs, "density", density_t)
      CALL torch_dict_insert(inputs, "grad", grad_t)
      CALL torch_dict_insert(inputs, "kin", kin_t)
      CALL torch_dict_insert(inputs, "coarse_0_atomic_coords", coarse_0_atomic_coords_t)

      CALL skala_torch_model_get_exc(cached_model, inputs, grid_weights_t, exc_tensor, exc)

      IF (.NOT. my_energy_only) THEN
         NULLIFY (atom_coord_grad, density_grad, grad_grad, grid_coord_grad, kin_grad)
         CALL torch_tensor_backward_scalar(exc_tensor)
         IF (need_coord_grad) THEN
            CALL torch_tensor_grad(grid_coords_t, grid_coord_grad_t)
            CALL torch_tensor_grad(coarse_0_atomic_coords_t, atom_coord_grad_t)
            CALL torch_tensor_data_ptr(grid_coord_grad_t, grid_coord_grad)
            CALL torch_tensor_data_ptr(atom_coord_grad_t, atom_coord_grad)
            IF (PRESENT(atom_force)) THEN
               atom_force(:) = atom_coord_grad(:, 1)
               DO row = 1, nflat
                  atom_force(:) = atom_force(:) + grid_coord_grad(:, row)
               END DO
            END IF
            IF (PRESENT(atom_virial)) THEN
               DO row = 1, nflat
                  DO idir = 1, 3
                     DO jdir = 1, 3
                        tmp = grid_coord_grad(idir, row)*coarse_0_atomic_coords(jdir, 1)
                        atom_virial(idir, jdir) = atom_virial(idir, jdir) + tmp
                     END DO
                  END DO
               END DO
               DO idir = 1, 3
                  DO jdir = 1, 3
                     tmp = atom_coord_grad(idir, 1)*coarse_0_atomic_coords(jdir, 1)
                     atom_virial(idir, jdir) = atom_virial(idir, jdir) + tmp
                  END DO
               END DO
            END IF
         END IF
         CALL torch_tensor_grad_batch3(density_t, grad_t, kin_t, density_grad_t, &
                                       grad_grad_t, kin_grad_t)
         CALL torch_tensor_data_ptr(density_grad_t, density_grad)
         CALL torch_tensor_data_ptr(grad_grad_t, grad_grad)
         CALL torch_tensor_data_ptr(kin_grad_t, kin_grad)

         row = 0
         DO ir = 1, nr
            DO ia = 1, na
               row = row + 1
               IF (lsd) THEN
                  vxc(ia, ir, 1:2) = density_grad(row, 1:2)
                  DO idir = 1, 3
                     vxg(idir, ia, ir, 1:2) = grad_grad(row, idir, 1:2)
                  END DO
                  vtau(ia, ir, 1:2) = kin_grad(row, 1:2)
               ELSE
                  vxc(ia, ir, 1) = 0.5_dp*(density_grad(row, 1) + density_grad(row, 2))
                  DO idir = 1, 3
                     vxg(idir, ia, ir, 1) = &
                        0.5_dp*(grad_grad(row, idir, 1) + grad_grad(row, idir, 2))
                  END DO
                  vtau(ia, ir, 1) = 0.5_dp*(kin_grad(row, 1) + kin_grad(row, 2))
               END IF
            END DO
         END DO

         CALL torch_tensor_release(density_grad_t)
         CALL torch_tensor_release(grad_grad_t)
         CALL torch_tensor_release(kin_grad_t)
         IF (need_coord_grad) THEN
            CALL torch_tensor_release(grid_coord_grad_t)
            CALL torch_tensor_release(atom_coord_grad_t)
         END IF
      END IF

      CALL torch_tensor_release(exc_tensor)
      CALL torch_tensor_release(density_t)
      CALL torch_tensor_release(grad_t)
      CALL torch_tensor_release(kin_t)
      CALL torch_tensor_release(grid_coords_t)
      CALL torch_tensor_release(grid_weights_t)
      CALL torch_tensor_release(atomic_grid_weights_t)
      CALL torch_tensor_release(atomic_grid_sizes_t)
      CALL torch_tensor_release(atomic_grid_size_bound_shape_t)
      CALL torch_tensor_release(coarse_0_atomic_coords_t)
      CALL torch_dict_release(inputs)
      DEALLOCATE (atomic_grid_size_bound_shape, atomic_grid_sizes, atomic_grid_weights, &
                  coarse_0_atomic_coords, density, grad, grid_coords, grid_weights, kin)
      CALL torch_use_cuda(.TRUE.)
      CALL omp_set_num_threads(omp_max_threads_restore)

   END SUBROUTINE skala_gapw_atom_vxc_of_r

! **************************************************************************************************
!> \brief Evaluate a rank-local set of complete atom blocks and sum their SKALA energies.
!> \param xc_section ...
!> \param group ...
!> \param density ...
!> \param grad ...
!> \param kin ...
!> \param grid_coords ...
!> \param grid_weights ...
!> \param atomic_grid_weights ...
!> \param atomic_grid_sizes ...
!> \param atomic_coords ...
!> \param exc ...
!> \param density_grad_out ...
!> \param grad_grad_out ...
!> \param kin_grad_out ...
!> \param grid_coord_grad_out ...
!> \param grid_weight_grad_out ...
!> \param atomic_grid_weight_grad_out ...
!> \param atom_coord_grad_out ...
!> \param density_constraints Selected physical descriptor derivatives of the KS environment.
!> \param constraint_domain Exact atomic or interpolated-background evaluation identifier.
!> \param atomic_ids Global atom identifiers in the order of atomic_grid_sizes.
!> \param density_direction Joint density direction on a differentiable model branch.
!> \param grad_direction Matching density-gradient direction.
!> \param kin_direction Matching kinetic-density direction.
!> \param density_response_out Response of the weighted density covector.
!> \param grad_response_out Response of the weighted density-gradient covector.
!> \param kin_response_out Response of the weighted kinetic-density covector.
!>        Responses require first derivatives and no selected nonsmooth descriptor limits.
! **************************************************************************************************
   SUBROUTINE skala_gapw_atom_composite_energy(xc_section, group, density, grad, kin, &
                                               grid_coords, grid_weights, atomic_grid_weights, &
                                               atomic_grid_sizes, atomic_coords, exc, &
                                               density_grad_out, grad_grad_out, kin_grad_out, &
                                               grid_coord_grad_out, grid_weight_grad_out, &
                                               atomic_grid_weight_grad_out, atom_coord_grad_out, &
                                               density_constraints, constraint_domain, atomic_ids, &
                                               density_direction, grad_direction, kin_direction, &
                                               density_response_out, grad_response_out, kin_response_out)
      TYPE(section_vals_type), POINTER                   :: xc_section

      CLASS(mp_comm_type), INTENT(IN)                    :: group
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), INTENT(IN) :: density, grid_coords, &
                                                                 atomic_coords, kin
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), INTENT(IN) :: grad
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), INTENT(IN) :: grid_weights, atomic_grid_weights
      INTEGER(KIND=int_8), ALLOCATABLE, DIMENSION(:), INTENT(IN) :: atomic_grid_sizes
      REAL(KIND=dp), INTENT(OUT)                         :: exc
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), INTENT(OUT), OPTIONAL :: density_grad_out, &
                                                                            kin_grad_out, &
                                                                            grid_coord_grad_out, &
                                                                            atom_coord_grad_out
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), INTENT(OUT), OPTIONAL :: grad_grad_out
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), INTENT(OUT), OPTIONAL :: grid_weight_grad_out, &
                                                                         atomic_grid_weight_grad_out
      TYPE(skala_density_constraints_type), OPTIONAL, POINTER :: density_constraints
      INTEGER, INTENT(IN), OPTIONAL                      :: constraint_domain, atomic_ids(:)
      REAL(dp), INTENT(IN), OPTIONAL                     :: density_direction(:, :), grad_direction(:, :, :), &
                                                            kin_direction(:, :)
      REAL(dp), ALLOCATABLE, INTENT(OUT), OPTIONAL        :: density_response_out(:, :), &
                                                             grad_response_out(:, :, :), kin_response_out(:, :)

      CHARACTER(len=default_path_length)                 :: model_path
      INTEGER :: atom_begin, atom_count, chunk_max_rows, chunk_row_count, chunk_row_start, &
                 ichunk, local_natom, local_nrow, max_grid_size, native_grid_cuda_device, &
                 nchunks, omp_max_threads_restore, phase_handle, selected_cuda_device
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: chunk_atom_begin, chunk_atom_count, &
                                                            chunk_row_begin, chunk_row_counts
      INTEGER, ALLOCATABLE                               :: selected_rows(:), selected_spins(:), &
                                                            selected_descriptors(:)
      INTEGER(KIND=int_8), ALLOCATABLE, DIMENSION(:, :)  :: atomic_grid_size_bound_shape
      LOGICAL :: active_rank, capture_jumps, native_grid_atom_chunks, native_grid_use_cuda, &
                 need_coordinate_derivatives, need_derivatives, need_response, use_constraints
      LOGICAL, ALLOCATABLE                               :: selected_applied(:)
      REAL(KIND=dp)                                      :: chunk_exc, probes(2), tolerances(2)
      REAL(dp), ALLOCATABLE                              :: selected_weights(:), selected_targets(:), selected_limits(:, :), &
                                                            selected_jumps(:, :), jump_values(:)
      INTEGER                                            :: constraint, first_spin, last_spin
      REAL(KIND=dp), DIMENSION(:), POINTER               :: atomic_grid_weight_grad, grid_weight_grad
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: atom_coord_grad, density_grad, &
                                                            grid_coord_grad, kin_grad
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: grad_grad
      TYPE(section_vals_type), POINTER                   :: gauxc_section
      TYPE(torch_dict_type)                              :: inputs
      TYPE(torch_tensor_type)                            :: atomic_grid_size_bound_shape_t, &
                                                            atomic_grid_sizes_t, &
                                                            atomic_grid_sizes_chunk_t, &
                                                            atomic_grid_weight_grad_t, &
                                                            atomic_grid_weights_chunk_t, &
                                                            atom_coord_grad_t, atomic_coords_chunk_t, &
                                                            density_chunk_t, density_grad_t, &
                                                            exc_tensor, grad_chunk_t, grad_grad_t, &
                                                            grid_coord_grad_t, grid_coords_chunk_t, &
                                                            grid_weight_grad_t, grid_weights_chunk_t, &
                                                            kin_chunk_t, kin_grad_t
      TYPE(torch_tensor_type), TARGET                    :: atomic_coords_local_t, &
                                                            atomic_grid_weights_local_t, &
                                                            density_local_t, grad_local_t, &
                                                            grid_coords_local_t, grid_weights_local_t, &
                                                            kin_local_t
      TYPE(torch_tensor_type), POINTER                   :: atomic_coords_t, atomic_grid_weights_t, &
                                                            density_t, grad_t, grid_coords_t, &
                                                            grid_weights_t, kin_t
      TYPE(torch_tensor_type)                            :: response_fields(3), response_directions(3), response_tensors(3)
      REAL(dp), ALLOCATABLE                              :: density_direction_buffer(:, :), &
                                                            grad_direction_buffer(:, :, :), kin_direction_buffer(:, :)
      REAL(dp), POINTER                                  :: density_response(:, :), grad_response(:, :, :), kin_response(:, :)
      INTEGER                                            :: component

      CPASSERT(ASSOCIATED(xc_section))
      CPASSERT(SIZE(density, 1) == SIZE(grid_weights))
      CPASSERT(SIZE(density, 1) == SIZE(atomic_grid_weights))
      CPASSERT(SIZE(density, 1) == SIZE(grid_coords, 2))
      CPASSERT(SIZE(density, 1) == SIZE(grad, 1))
      CPASSERT(SIZE(density, 1) == SIZE(kin, 1))
      CPASSERT(SIZE(density, 2) == 2)
      CPASSERT(SIZE(grad, 2) == 3)
      CPASSERT(SIZE(grad, 3) == 2)
      CPASSERT(SIZE(kin, 2) == 2)
      CPASSERT(SIZE(atomic_grid_sizes) == SIZE(atomic_coords, 2))
      CPASSERT(SUM(atomic_grid_sizes) == INT(SIZE(density, 1), KIND=int_8))
      need_derivatives = PRESENT(density_grad_out)
      CPASSERT(PRESENT(grad_grad_out) .EQV. need_derivatives)
      CPASSERT(PRESENT(kin_grad_out) .EQV. need_derivatives)
      need_response = PRESENT(density_response_out)
      CPASSERT(PRESENT(grad_response_out) .EQV. need_response)
      CPASSERT(PRESENT(kin_response_out) .EQV. need_response)
      CPASSERT(PRESENT(density_direction) .EQV. need_response)
      CPASSERT(PRESENT(grad_direction) .EQV. need_response)
      CPASSERT(PRESENT(kin_direction) .EQV. need_response)
      IF (need_response) THEN
         CPASSERT(need_derivatives)
         CPASSERT(ALL(SHAPE(density_direction) == SHAPE(density)))
         CPASSERT(ALL(SHAPE(grad_direction) == SHAPE(grad)))
         CPASSERT(ALL(SHAPE(kin_direction) == SHAPE(kin)))
         ALLOCATE (density_response_out(SIZE(density, 1), 2), &
                   grad_response_out(SIZE(density, 1), 3, 2), kin_response_out(SIZE(density, 1), 2))
         density_response_out = 0.0_dp
         grad_response_out = 0.0_dp
         kin_response_out = 0.0_dp
      END IF
      need_coordinate_derivatives = PRESENT(grid_coord_grad_out)
      CPASSERT(PRESENT(grid_weight_grad_out) .EQV. need_coordinate_derivatives)
      CPASSERT(PRESENT(atomic_grid_weight_grad_out) .EQV. need_coordinate_derivatives)
      CPASSERT(PRESENT(atom_coord_grad_out) .EQV. need_coordinate_derivatives)
      CPASSERT((.NOT. need_coordinate_derivatives) .OR. need_derivatives)
      omp_max_threads_restore = omp_get_max_threads()

      local_nrow = SIZE(density, 1)
      local_natom = SIZE(atomic_grid_sizes)
      active_rank = local_natom > 0
      CPASSERT(active_rank .EQV. (local_nrow > 0))
      use_constraints = .FALSE.
      capture_jumps = .FALSE.
      IF (PRESENT(density_constraints) .AND. need_derivatives) THEN
         IF (ASSOCIATED(density_constraints)) THEN
            CPASSERT(PRESENT(constraint_domain) .AND. PRESENT(atomic_ids))
            use_constraints = skala_density_constraints_count(density_constraints, constraint_domain) > 0
            IF (use_constraints) THEN
               CALL prepare_atom_constraints(density_constraints, constraint_domain, atomic_ids, &
                                             atomic_grid_sizes, group, selected_rows, selected_spins, &
                                             selected_descriptors, selected_weights, selected_targets, probes, tolerances)
               ALLOCATE (selected_limits(SIZE(selected_rows), 2), selected_applied(SIZE(selected_rows)))
               selected_limits = 0.0_dp
               selected_applied = .FALSE.
               capture_jumps = skala_density_jumps_active(density_constraints) .AND. &
                               (constraint_domain == skala_atom_domain .OR. &
                                constraint_domain == skala_background_domain)
               IF (capture_jumps) THEN
                  ALLOCATE (selected_jumps(SIZE(selected_rows), 2), jump_values(SIZE(selected_rows)))
                  selected_jumps = 0.0_dp
                  jump_values = 0.0_dp
               END IF
            END IF
         END IF
      END IF

      IF (need_response) THEN
         CPASSERT(.NOT. use_constraints)
      END IF
      CALL get_skala_model_path(xc_section, model_path)
      gauxc_section => get_gauxc_section(xc_section)
      CPASSERT(ASSOCIATED(gauxc_section))
      CALL section_vals_val_get(gauxc_section, "NATIVE_GRID_USE_CUDA", l_val=native_grid_use_cuda)
      CALL section_vals_val_get(gauxc_section, "NATIVE_GRID_CUDA_DEVICE", &
                                i_val=native_grid_cuda_device)
      CALL section_vals_val_get(gauxc_section, "NATIVE_GRID_ATOM_CHUNKS", &
                                l_val=native_grid_atom_chunks)
      CALL section_vals_val_get(gauxc_section, "NATIVE_GRID_ATOM_CHUNK_MAX_ROWS", &
                                i_val=chunk_max_rows)
      IF (chunk_max_rows < -1) THEN
         CALL cp_abort(__LOCATION__, &
                       "GAUXC%NATIVE_GRID_ATOM_CHUNK_MAX_ROWS must be -1, zero, or positive.")
      END IF
      IF (.NOT. native_grid_atom_chunks) chunk_max_rows = 0
      IF (chunk_max_rows == -1) THEN
         chunk_max_rows = auto_atom_composite_chunk_max_rows(atomic_grid_sizes, group)
      END IF
      CALL torch_use_cuda(native_grid_use_cuda)
      selected_cuda_device = configure_native_grid_cuda( &
                             native_grid_use_cuda, native_grid_cuda_device, group)
      exc = 0.0_dp
      IF (active_rank) THEN
         CALL timeset("skala_atom_tensor_setup", phase_handle)
         CALL ensure_model_loaded(model_path, selected_cuda_device, xc_section)
         CALL atom_composite_chunk_layout(atomic_grid_sizes, chunk_max_rows, &
                                          chunk_atom_begin, chunk_atom_count, &
                                          chunk_row_begin, chunk_row_counts)
         nchunks = SIZE(chunk_atom_begin)
         CPASSERT(nchunks > 0)

         IF (native_grid_use_cuda) THEN
            grid_coords_t => cached_atom_cuda_tensors%grid_coords_t
            grid_weights_t => cached_atom_cuda_tensors%grid_weights_t
            atomic_grid_weights_t => cached_atom_cuda_tensors%atomic_grid_weights_t
            atomic_coords_t => cached_atom_cuda_tensors%atomic_coords_t
            density_t => cached_atom_cuda_tensors%density_t
            grad_t => cached_atom_cuda_tensors%grad_t
            kin_t => cached_atom_cuda_tensors%kin_t
            CALL torch_tensor_reset_from_array( &
               grid_coords_t, grid_coords, need_coordinate_derivatives)
            CALL torch_tensor_reset_from_array( &
               grid_weights_t, grid_weights, need_coordinate_derivatives)
            CALL torch_tensor_reset_from_array( &
               atomic_grid_weights_t, atomic_grid_weights, need_coordinate_derivatives)
            CALL torch_tensor_reset_from_array( &
               atomic_coords_t, atomic_coords, need_coordinate_derivatives)
            CALL torch_tensor_reset_from_array(density_t, density, need_derivatives)
            CALL torch_tensor_reset_from_array(grad_t, grad, need_derivatives)
            CALL torch_tensor_reset_from_array(kin_t, kin, need_derivatives)
         ELSE
            grid_coords_t => grid_coords_local_t
            grid_weights_t => grid_weights_local_t
            atomic_grid_weights_t => atomic_grid_weights_local_t
            atomic_coords_t => atomic_coords_local_t
            density_t => density_local_t
            grad_t => grad_local_t
            kin_t => kin_local_t
            CALL torch_tensor_from_array(grid_coords_t, grid_coords)
            CALL torch_tensor_to_device_leaf(grid_coords_t, need_coordinate_derivatives)
            CALL torch_tensor_from_array(grid_weights_t, grid_weights)
            CALL torch_tensor_to_device_leaf(grid_weights_t, need_coordinate_derivatives)
            CALL torch_tensor_from_array(atomic_grid_weights_t, atomic_grid_weights)
            CALL torch_tensor_to_device_leaf(atomic_grid_weights_t, need_coordinate_derivatives)
            CALL torch_tensor_from_array(atomic_coords_t, atomic_coords)
            CALL torch_tensor_to_device_leaf(atomic_coords_t, need_coordinate_derivatives)
            CALL torch_tensor_from_array(density_t, density)
            CALL torch_tensor_to_device_leaf(density_t, need_derivatives)
            CALL torch_tensor_from_array(grad_t, grad)
            CALL torch_tensor_to_device_leaf(grad_t, need_derivatives)
            CALL torch_tensor_from_array(kin_t, kin)
            CALL torch_tensor_to_device_leaf(kin_t, need_derivatives)
         END IF
         CALL torch_tensor_from_array(atomic_grid_sizes_t, atomic_grid_sizes)
         CALL torch_tensor_to_device_leaf(atomic_grid_sizes_t, .FALSE.)
         CALL timestop(phase_handle)

         DO ichunk = 1, nchunks
            atom_begin = chunk_atom_begin(ichunk)
            atom_count = chunk_atom_count(ichunk)
            chunk_row_start = chunk_row_begin(ichunk)
            chunk_row_count = chunk_row_counts(ichunk)
            max_grid_size = INT(MAXVAL(atomic_grid_sizes( &
                                       atom_begin:atom_begin + atom_count - 1)))
            ALLOCATE (atomic_grid_size_bound_shape(0, max_grid_size))
            atomic_grid_size_bound_shape = 0_int_8

            CALL torch_tensor_narrow(grid_coords_t, 0, chunk_row_start - 1, &
                                     chunk_row_count, grid_coords_chunk_t)
            CALL torch_tensor_narrow(grid_weights_t, 0, chunk_row_start - 1, &
                                     chunk_row_count, grid_weights_chunk_t)
            CALL torch_tensor_narrow(atomic_grid_weights_t, 0, chunk_row_start - 1, &
                                     chunk_row_count, atomic_grid_weights_chunk_t)
            CALL torch_tensor_narrow(atomic_grid_sizes_t, 0, atom_begin - 1, &
                                     atom_count, atomic_grid_sizes_chunk_t)
            CALL torch_tensor_narrow(atomic_coords_t, 0, atom_begin - 1, &
                                     atom_count, atomic_coords_chunk_t)
            CALL torch_tensor_narrow(density_t, 1, chunk_row_start - 1, &
                                     chunk_row_count, density_chunk_t)
            CALL torch_tensor_narrow(grad_t, 2, chunk_row_start - 1, &
                                     chunk_row_count, grad_chunk_t)
            CALL torch_tensor_narrow(kin_t, 1, chunk_row_start - 1, &
                                     chunk_row_count, kin_chunk_t)
            CALL torch_tensor_from_array(atomic_grid_size_bound_shape_t, &
                                         atomic_grid_size_bound_shape)
            CALL torch_tensor_to_device_leaf(atomic_grid_size_bound_shape_t, .FALSE.)

            CALL torch_dict_create(inputs)
            CALL torch_dict_insert(inputs, "grid_coords", grid_coords_chunk_t)
            CALL torch_dict_insert(inputs, "grid_weights", grid_weights_chunk_t)
            CALL torch_dict_insert(inputs, "atomic_grid_weights", atomic_grid_weights_chunk_t)
            CALL torch_dict_insert(inputs, "atomic_grid_sizes", atomic_grid_sizes_chunk_t)
            CALL torch_dict_insert(inputs, "atomic_grid_size_bound_shape", &
                                   atomic_grid_size_bound_shape_t)
            CALL torch_dict_insert(inputs, "density", density_chunk_t)
            CALL torch_dict_insert(inputs, "grad", grad_chunk_t)
            CALL torch_dict_insert(inputs, "kin", kin_chunk_t)
            CALL torch_dict_insert(inputs, "coarse_0_atomic_coords", atomic_coords_chunk_t)

            CALL timeset("skala_atom_torch_forward", phase_handle)
            CALL skala_torch_model_get_exc(cached_model, inputs, grid_weights_chunk_t, &
                                           exc_tensor, chunk_exc)
            CALL timestop(phase_handle)
            exc = exc + chunk_exc
            IF (need_response) THEN
               ! The three primitive fields enter the same nonlocal energy graph.
               ! Differentiate them jointly before the ordinary backward releases it.
               ! Keep contiguous buffers alive while Torch borrows their storage.
               ALLOCATE (density_direction_buffer, &
                         SOURCE=density_direction(chunk_row_start:chunk_row_start + chunk_row_count - 1, :))
               ALLOCATE (grad_direction_buffer, &
                         SOURCE=grad_direction(chunk_row_start:chunk_row_start + chunk_row_count - 1, :, :))
               ALLOCATE (kin_direction_buffer, &
                         SOURCE=kin_direction(chunk_row_start:chunk_row_start + chunk_row_count - 1, :))
               CALL torch_tensor_from_array(response_directions(1), density_direction_buffer)
               CALL torch_tensor_from_array(response_directions(2), grad_direction_buffer)
               CALL torch_tensor_from_array(response_directions(3), kin_direction_buffer)
               response_fields = [density_chunk_t, grad_chunk_t, kin_chunk_t]
               CALL torch_tensor_hessian_product3(exc_tensor, response_fields, response_directions, response_tensors)
               NULLIFY (density_response, grad_response, kin_response)
               CALL torch_tensor_data_ptr(response_tensors(1), density_response)
               CALL torch_tensor_data_ptr(response_tensors(2), grad_response)
               CALL torch_tensor_data_ptr(response_tensors(3), kin_response)
               density_response_out(chunk_row_start:chunk_row_start + chunk_row_count - 1, :) = density_response
               grad_response_out(chunk_row_start:chunk_row_start + chunk_row_count - 1, :, :) = grad_response
               kin_response_out(chunk_row_start:chunk_row_start + chunk_row_count - 1, :) = kin_response
               DO component = 1, 3
                  CALL torch_tensor_release(response_tensors(component))
                  CALL torch_tensor_release(response_directions(component))
               END DO
               DEALLOCATE (density_direction_buffer, grad_direction_buffer, kin_direction_buffer)
            END IF
            IF (need_derivatives) THEN
               CALL timeset("skala_atom_torch_backward", phase_handle)
               CALL torch_tensor_backward_scalar(exc_tensor)
               CALL timestop(phase_handle)
            END IF
            IF (use_constraints) THEN
               IF (capture_jumps) THEN
                  CALL atom_descriptor_limits(inputs, density(chunk_row_start:chunk_row_start + chunk_row_count - 1, :), &
                                              kin(chunk_row_start:chunk_row_start + chunk_row_count - 1, :), &
                                              chunk_row_start, selected_rows, selected_spins, selected_descriptors, &
                                              selected_weights, selected_targets, probes, tolerances, &
                                              selected_limits, selected_applied, selected_jumps)
               ELSE
                  CALL atom_descriptor_limits(inputs, density(chunk_row_start:chunk_row_start + chunk_row_count - 1, :), &
                                              kin(chunk_row_start:chunk_row_start + chunk_row_count - 1, :), &
                                              chunk_row_start, selected_rows, selected_spins, selected_descriptors, &
                                              selected_weights, selected_targets, probes, tolerances, &
                                              selected_limits, selected_applied)
               END IF
            END IF

            CALL torch_tensor_release(exc_tensor)
            CALL torch_dict_release(inputs)
            CALL torch_tensor_release(grid_coords_chunk_t)
            CALL torch_tensor_release(grid_weights_chunk_t)
            CALL torch_tensor_release(atomic_grid_weights_chunk_t)
            CALL torch_tensor_release(atomic_grid_sizes_chunk_t)
            CALL torch_tensor_release(atomic_coords_chunk_t)
            CALL torch_tensor_release(density_chunk_t)
            CALL torch_tensor_release(grad_chunk_t)
            CALL torch_tensor_release(kin_chunk_t)
            CALL torch_tensor_release(atomic_grid_size_bound_shape_t)
            DEALLOCATE (atomic_grid_size_bound_shape)
         END DO
      END IF
      CALL group%sum(exc)

      IF (need_derivatives) THEN
         ALLOCATE (density_grad_out(local_nrow, 2), grad_grad_out(local_nrow, 3, 2), &
                   kin_grad_out(local_nrow, 2))
         density_grad_out = 0.0_dp
         grad_grad_out = 0.0_dp
         kin_grad_out = 0.0_dp
         IF (active_rank) THEN
            NULLIFY (density_grad, grad_grad, kin_grad)
            CALL timeset("skala_atom_grad_fetch", phase_handle)
            CALL torch_tensor_grad_batch3(density_t, grad_t, kin_t, density_grad_t, &
                                          grad_grad_t, kin_grad_t)
            CALL torch_tensor_data_ptr(density_grad_t, density_grad)
            CALL torch_tensor_data_ptr(grad_grad_t, grad_grad)
            CALL torch_tensor_data_ptr(kin_grad_t, kin_grad)
            density_grad_out(:, :) = density_grad
            grad_grad_out(:, :, :) = grad_grad
            kin_grad_out(:, :) = kin_grad
            CALL timestop(phase_handle)
         END IF
         IF (use_constraints) THEN
            DO constraint = 1, SIZE(selected_rows)
               IF (.NOT. selected_applied(constraint)) CYCLE
               first_spin = MAX(1, selected_spins(constraint))
               last_spin = selected_spins(constraint)
               IF (last_spin == 0) last_spin = 2
               IF (capture_jumps) THEN
                  jump_values(constraint) = SUM(selected_jumps(constraint, first_spin:last_spin))
               END IF
               IF (selected_descriptors(constraint) == skala_density_descriptor) THEN
                  density_grad_out(selected_rows(constraint), first_spin:last_spin) = &
                     selected_limits(constraint, first_spin:last_spin)
               ELSE
                  kin_grad_out(selected_rows(constraint), first_spin:last_spin) = &
                     selected_limits(constraint, first_spin:last_spin)
               END IF
            END DO
            IF (capture_jumps) THEN
               ! Torch covectors already contain the quadrature weights. The interpolated
               ! background is subtracted from the full GAPW energy and operator.
               IF (constraint_domain == skala_background_domain) jump_values = -jump_values
               CALL skala_density_jumps_put(density_constraints, constraint_domain, jump_values, selected_applied)
            END IF
         END IF

         IF (need_coordinate_derivatives) THEN
            IF (active_rank) CALL timeset("skala_atom_coord_grad_fetch", phase_handle)
            ALLOCATE (grid_coord_grad_out(3, local_nrow), grid_weight_grad_out(local_nrow), &
                      atomic_grid_weight_grad_out(local_nrow), atom_coord_grad_out(3, local_natom))
            grid_coord_grad_out = 0.0_dp
            grid_weight_grad_out = 0.0_dp
            atomic_grid_weight_grad_out = 0.0_dp
            atom_coord_grad_out = 0.0_dp
            IF (active_rank) THEN
               NULLIFY (atomic_grid_weight_grad, atom_coord_grad, grid_coord_grad, grid_weight_grad)
               CALL torch_tensor_grad(grid_coords_t, grid_coord_grad_t)
               CALL torch_tensor_grad(grid_weights_t, grid_weight_grad_t)
               CALL torch_tensor_grad(atomic_grid_weights_t, atomic_grid_weight_grad_t)
               CALL torch_tensor_grad(atomic_coords_t, atom_coord_grad_t)
               CALL torch_tensor_data_ptr(grid_coord_grad_t, grid_coord_grad)
               CALL torch_tensor_data_ptr(grid_weight_grad_t, grid_weight_grad)
               CALL torch_tensor_data_ptr(atomic_grid_weight_grad_t, atomic_grid_weight_grad)
               CALL torch_tensor_data_ptr(atom_coord_grad_t, atom_coord_grad)
               grid_coord_grad_out(:, :) = grid_coord_grad
               grid_weight_grad_out(:) = grid_weight_grad
               atomic_grid_weight_grad_out(:) = atomic_grid_weight_grad
               atom_coord_grad_out(:, :) = atom_coord_grad
               CALL torch_tensor_release(grid_coord_grad_t)
               CALL torch_tensor_release(grid_weight_grad_t)
               CALL torch_tensor_release(atomic_grid_weight_grad_t)
               CALL torch_tensor_release(atom_coord_grad_t)
               CALL timestop(phase_handle)
            END IF
         END IF
         IF (active_rank) THEN
            CALL torch_tensor_release(density_grad_t)
            CALL torch_tensor_release(grad_grad_t)
            CALL torch_tensor_release(kin_grad_t)
         END IF
      END IF

      IF (active_rank) THEN
         CALL timeset("skala_atom_tensor_release", phase_handle)
         IF (.NOT. native_grid_use_cuda) THEN
            CALL torch_tensor_release(density_t)
            CALL torch_tensor_release(grad_t)
            CALL torch_tensor_release(kin_t)
            CALL torch_tensor_release(grid_coords_t)
            CALL torch_tensor_release(grid_weights_t)
            CALL torch_tensor_release(atomic_grid_weights_t)
            CALL torch_tensor_release(atomic_coords_t)
         END IF
         CALL torch_tensor_release(atomic_grid_sizes_t)
         DEALLOCATE (chunk_atom_begin, chunk_atom_count, chunk_row_begin, chunk_row_counts)
         CALL timestop(phase_handle)
      END IF
      CALL torch_use_cuda(.TRUE.)
      CALL omp_set_num_threads(omp_max_threads_restore)

   END SUBROUTINE skala_gapw_atom_composite_energy

! **************************************************************************************************
!> \brief Match selected global atomic rows to their single complete-atom owner.
!> \param state ...
!> \param domain ...
!> \param atomic_ids ...
!> \param atomic_sizes ...
!> \param group ...
!> \param rows Zero on ranks that do not own the selected atom.
!> \param spins ...
!> \param descriptors ...
!> \param weights ...
!> \param targets Registered physical target of each selected descriptor.
!> \param probes ...
!> \param tolerances ...
! **************************************************************************************************
   SUBROUTINE prepare_atom_constraints(state, domain, atomic_ids, atomic_sizes, group, &
                                       rows, spins, descriptors, weights, targets, probes, tolerances)
      TYPE(skala_density_constraints_type), INTENT(IN)   :: state
      INTEGER, INTENT(IN)                                :: domain, atomic_ids(:)
      INTEGER(int_8), INTENT(IN)                         :: atomic_sizes(:)

      CLASS(mp_comm_type), INTENT(IN)                    :: group
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: rows(:), spins(:), descriptors(:)
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: weights(:), targets(:)
      REAL(dp), INTENT(OUT)                              :: probes(2), tolerances(2)

      INTEGER                                            :: atom, c, offset
      INTEGER, ALLOCATABLE                               :: owners(:), points(:, :)

      CPASSERT(SIZE(atomic_ids) == SIZE(atomic_sizes))
      CALL skala_density_constraints_get(state, points, spins, weights, probes(1), tolerances(1), &
                                         descriptors, probes(2), tolerances(2), domain=domain, targets=targets)
      ALLOCATE (rows(SIZE(spins)), owners(SIZE(spins)))
      rows = 0
      owners = 0
      offset = 0
      DO atom = 1, SIZE(atomic_ids)
         DO c = 1, SIZE(spins)
            IF (points(1, c) /= atomic_ids(atom)) CYCLE
            CPASSERT(points(2, c) <= atomic_sizes(atom))
            rows(c) = offset + points(2, c)
            owners(c) = owners(c) + 1
         END DO
         offset = offset + INT(atomic_sizes(atom))
      END DO
      CALL group%sum(owners)
      CPASSERT(ALL(owners == 1))
   END SUBROUTINE prepare_atom_constraints

! **************************************************************************************************
!> \brief Obtain selected density and tau limits from the unchanged complete atomic model inputs.
!>        The original energy and all other covectors remain unchanged. Away from the physical
!>        target, the original smooth derivative is retained. A paired spin constraint probes both spins.
!> \param inputs ...
!> \param density ...
!> \param kinetic ...
!> \param first_row ...
!> \param rows ...
!> \param spins ...
!> \param descriptors ...
!> \param weights ...
!> \param targets Registered physical target of each selected descriptor.
!> \param probes ...
!> \param tolerances ...
!> \param limits ...
!> \param applied ...
!> \param jumps Optional already weighted negative-minus-positive primitive covectors.
! **************************************************************************************************
   SUBROUTINE atom_descriptor_limits(inputs, density, kinetic, first_row, rows, spins, descriptors, &
                                     weights, targets, probes, tolerances, limits, applied, jumps)
      TYPE(torch_dict_type), INTENT(IN)                  :: inputs
      REAL(dp), INTENT(IN)                               :: density(:, :), kinetic(:, :)
      INTEGER, INTENT(IN)                                :: first_row, rows(:), spins(:), &
                                                            descriptors(:)
      REAL(dp), INTENT(IN)                               :: weights(:), targets(:), probes(2), &
                                                            tolerances(2)
      REAL(dp), INTENT(INOUT)                            :: limits(:, :)
      LOGICAL, INTENT(INOUT)                             :: applied(:)
      REAL(dp), INTENT(INOUT), OPTIONAL                  :: jumps(:, :)

      INTEGER                                            :: c, descriptor, first_spin, last_spin, n, &
                                                            row, spin
      INTEGER, ALLOCATABLE                               :: active(:, :)
      LOGICAL                                            :: branches(2)
      LOGICAL, ALLOCATABLE                               :: selected(:)
      REAL(dp)                                           :: error
      REAL(dp), ALLOCATABLE                              :: active_targets(:), negative(:, :), &
                                                            positive(:, :)

      CPASSERT(SIZE(targets) == SIZE(rows))
      ALLOCATE (selected(SIZE(rows)))
      DO descriptor = skala_density_descriptor, skala_kinetic_descriptor
         selected = .FALSE.
         n = 0
         DO c = 1, SIZE(rows)
            IF (descriptors(c) /= descriptor .OR. rows(c) == 0) CYCLE
            row = rows(c) - first_row + 1
            IF (row < 1 .OR. row > SIZE(density, 1)) CYCLE
            first_spin = MAX(1, spins(c))
            last_spin = spins(c)
            IF (last_spin == 0) last_spin = 2
            IF (descriptor == skala_density_descriptor) THEN
               error = MAXVAL(ABS(density(row, first_spin:last_spin) - targets(c)))
            ELSE
               error = MAXVAL(ABS(kinetic(row, first_spin:last_spin) - targets(c)))
            END IF
            IF (error > tolerances(descriptor)) CYCLE
            selected(c) = .TRUE.
            n = n + last_spin - first_spin + 1
         END DO
         IF (n == 0) CYCLE
         ALLOCATE (active(2, n), active_targets(n))
         n = 0
         DO c = 1, SIZE(rows)
            IF (.NOT. selected(c)) CYCLE
            first_spin = MAX(1, spins(c))
            last_spin = spins(c)
            IF (last_spin == 0) last_spin = 2
            DO spin = first_spin, last_spin
               n = n + 1
               active(:, n) = [rows(c) - first_row + 1, spin]
               active_targets(n) = targets(c)
            END DO
         END DO
         ! Keep the full probe input, but evaluate only covectors used by the selected weights.
         branches = [ANY(selected .AND. weights /= 0.0_dp), ANY(selected .AND. weights /= 1.0_dp)]
         IF (PRESENT(jumps)) branches = .TRUE.
         IF (descriptor == skala_density_descriptor) THEN
            CALL skala_torch_model_descriptor_limits(cached_model, inputs, "density", density, &
                                                     active, probes(descriptor), negative, positive, branches=branches, &
                                                     targets=active_targets, adjacent=cached_model_fp64)
         ELSE
            CALL skala_torch_model_descriptor_limits(cached_model, inputs, "kin", kinetic, &
                                                     active, probes(descriptor), negative, positive, branches=branches, &
                                                     targets=active_targets, adjacent=cached_model_fp64)
         END IF
         DO c = 1, SIZE(rows)
            IF (.NOT. selected(c)) CYCLE
            row = rows(c) - first_row + 1
            limits(c, :) = weights(c)*negative(row, :) + (1.0_dp - weights(c))*positive(row, :)
            IF (PRESENT(jumps)) jumps(c, :) = negative(row, :) - positive(row, :)
            applied(c) = .TRUE.
         END DO
         DEALLOCATE (active, active_targets)
      END DO
   END SUBROUTINE atom_descriptor_limits

! **************************************************************************************************
!> \brief Select an automatic row cap for atom-composite Torch calls.
!> \param atomic_grid_sizes ...
!> \param group ...
!> \return ...
! **************************************************************************************************
   FUNCTION auto_atom_composite_chunk_max_rows(atomic_grid_sizes, group) RESULT(max_rows)
      INTEGER(KIND=int_8), DIMENSION(:), INTENT(IN)      :: atomic_grid_sizes

      CLASS(mp_comm_type), INTENT(IN)                    :: group

      INTEGER                                            :: local_max_atom_rows, local_natom, &
                                                            local_padded_rows, max_rows

      local_padded_rows = 0
      local_natom = SIZE(atomic_grid_sizes)
      IF (local_natom > 0) THEN
         local_max_atom_rows = INT(MAXVAL(atomic_grid_sizes))
         IF (local_natom > atom_chunk_auto_max_rows .OR. &
             local_max_atom_rows > atom_chunk_auto_max_rows/local_natom) THEN
            local_padded_rows = atom_chunk_auto_max_rows + 1
         ELSE
            local_padded_rows = local_natom*local_max_atom_rows
         END IF
      END IF
      CALL group%max(local_padded_rows)
      IF (local_padded_rows <= atom_chunk_auto_max_rows) THEN
         max_rows = 0
         RETURN
      END IF
      max_rows = atom_chunk_auto_max_rows

   END FUNCTION auto_atom_composite_chunk_max_rows

! **************************************************************************************************
!> \brief Split contiguous atom-composite rows without dividing an atomic block.
!> \param atomic_grid_sizes ...
!> \param max_rows ...
!> \param atom_begin ...
!> \param atom_count ...
!> \param row_begin ...
!> \param row_count ...
! **************************************************************************************************
   SUBROUTINE atom_composite_chunk_layout(atomic_grid_sizes, max_rows, atom_begin, atom_count, &
                                          row_begin, row_count)
      INTEGER(KIND=int_8), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: atomic_grid_sizes
      INTEGER, INTENT(IN)                                :: max_rows
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: atom_begin, atom_count, row_begin, &
                                                            row_count

      INTEGER :: atom_begin_tmp, atom_count_tmp, atom_rows, chunk_atoms, chunk_max_atom_rows, &
         iatom, insert_at, natom, nchunks, padded_rows_tmp, row_begin_tmp, row_count_tmp, rows, &
         subchunk
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: padded_rows

      natom = SIZE(atomic_grid_sizes)
      CPASSERT(natom > 0)
      CPASSERT(ALL(atomic_grid_sizes > 0_int_8))
      IF (max_rows <= 0) THEN
         ALLOCATE (atom_begin(1), atom_count(1), row_begin(1), row_count(1))
         atom_begin(1) = 1
         atom_count(1) = natom
         row_begin(1) = 1
         row_count(1) = INT(SUM(atomic_grid_sizes))
         RETURN
      END IF

      nchunks = 1
      chunk_atoms = 0
      chunk_max_atom_rows = 0
      DO iatom = 1, natom
         atom_rows = INT(atomic_grid_sizes(iatom))
         IF (chunk_atoms > 0 .AND. &
             MAX(chunk_max_atom_rows, atom_rows) > max_rows/(chunk_atoms + 1)) THEN
            nchunks = nchunks + 1
            chunk_atoms = 0
            chunk_max_atom_rows = 0
         END IF
         chunk_atoms = chunk_atoms + 1
         chunk_max_atom_rows = MAX(chunk_max_atom_rows, atom_rows)
      END DO
      ALLOCATE (atom_begin(nchunks), atom_count(nchunks), row_begin(nchunks), row_count(nchunks))

      subchunk = 1
      atom_begin(subchunk) = 1
      row_begin(subchunk) = 1
      chunk_atoms = 0
      chunk_max_atom_rows = 0
      rows = 0
      DO iatom = 1, natom
         atom_rows = INT(atomic_grid_sizes(iatom))
         IF (chunk_atoms > 0 .AND. &
             MAX(chunk_max_atom_rows, atom_rows) > max_rows/(chunk_atoms + 1)) THEN
            atom_count(subchunk) = iatom - atom_begin(subchunk)
            row_count(subchunk) = rows
            subchunk = subchunk + 1
            atom_begin(subchunk) = iatom
            row_begin(subchunk) = row_begin(subchunk - 1) + row_count(subchunk - 1)
            chunk_atoms = 0
            chunk_max_atom_rows = 0
            rows = 0
         END IF
         chunk_atoms = chunk_atoms + 1
         chunk_max_atom_rows = MAX(chunk_max_atom_rows, atom_rows)
         rows = rows + atom_rows
      END DO
      atom_count(subchunk) = natom - atom_begin(subchunk) + 1
      row_count(subchunk) = rows
      CPASSERT(subchunk == nchunks)

      ALLOCATE (padded_rows(nchunks))
      DO subchunk = 1, nchunks
         iatom = atom_begin(subchunk) + atom_count(subchunk) - 1
         padded_rows(subchunk) = atom_count(subchunk)* &
                                 INT(MAXVAL(atomic_grid_sizes(atom_begin(subchunk):iatom)))
      END DO

      ! Largest blocks run first so CUDA allocations can be reused by all following calls.
      DO subchunk = 2, nchunks
         padded_rows_tmp = padded_rows(subchunk)
         atom_begin_tmp = atom_begin(subchunk)
         atom_count_tmp = atom_count(subchunk)
         row_begin_tmp = row_begin(subchunk)
         row_count_tmp = row_count(subchunk)
         insert_at = subchunk
         DO WHILE (insert_at > 1 .AND. padded_rows(insert_at - 1) < padded_rows_tmp)
            padded_rows(insert_at) = padded_rows(insert_at - 1)
            atom_begin(insert_at) = atom_begin(insert_at - 1)
            atom_count(insert_at) = atom_count(insert_at - 1)
            row_begin(insert_at) = row_begin(insert_at - 1)
            row_count(insert_at) = row_count(insert_at - 1)
            insert_at = insert_at - 1
         END DO
         padded_rows(insert_at) = padded_rows_tmp
         atom_begin(insert_at) = atom_begin_tmp
         atom_count(insert_at) = atom_count_tmp
         row_begin(insert_at) = row_begin_tmp
         row_count(insert_at) = row_count_tmp
      END DO
      DEALLOCATE (padded_rows)

   END SUBROUTINE atom_composite_chunk_layout

! **************************************************************************************************
!> \brief Add the explicit SKALA derivative with respect to atom-center coordinates.
!> \param atom_force ...
!> \param features ...
!> \param atom_coord_grad_t ...
!> \param root_rank ...
! **************************************************************************************************
   SUBROUTINE add_explicit_coordinate_force(atom_force, features, atom_coord_grad_t, root_rank)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: atom_force
      TYPE(skala_gpw_feature_type), INTENT(IN)           :: features
      TYPE(torch_tensor_type), INTENT(INOUT)             :: atom_coord_grad_t
      LOGICAL, INTENT(IN)                                :: root_rank

      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: atom_coord_grad

      NULLIFY (atom_coord_grad)
      CALL torch_tensor_grad(features%coarse_0_atomic_coords_t, atom_coord_grad_t)
      IF (root_rank) THEN
         CALL torch_tensor_data_ptr(atom_coord_grad_t, atom_coord_grad)
         CPASSERT(SIZE(atom_force, 1) == SIZE(atom_coord_grad, 1))
         CPASSERT(SIZE(atom_coord_grad, 2) >= SIZE(atom_force, 2))
         CPASSERT(SIZE(atom_coord_grad, 2) <= SIZE(atom_force, 2) + 1)
         atom_force(:, :) = atom_force(:, :) + atom_coord_grad(:, :SIZE(atom_force, 2))
      END IF

   END SUBROUTINE add_explicit_coordinate_force

! **************************************************************************************************
!> \brief Add the force from SMOOTH native-grid atom partition weights.
!> \param atom_force ...
!> \param features ...
!> \param particle_set ...
!> \param cell ...
!> \param rho_r ...
!> \param grid_weight_grad_t ...
!> \param atomic_grid_weight_grad_t ...
! **************************************************************************************************
   SUBROUTINE add_smooth_partition_force(atom_force, features, particle_set, cell, rho_r, &
                                         grid_weight_grad_t, atomic_grid_weight_grad_t)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: atom_force
      TYPE(skala_gpw_feature_type), INTENT(IN)           :: features
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: rho_r
      TYPE(torch_tensor_type), INTENT(INOUT)             :: grid_weight_grad_t, &
                                                            atomic_grid_weight_grad_t

      INTEGER                                            :: feature_begin, feature_end, feature_pos, &
                                                            i, iatom, j, jatom, k, local_row, &
                                                            natom, row
      INTEGER, DIMENSION(2, 3)                           :: bo
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: included
      REAL(dp)                                           :: periodic_strain_response(3, 3)
      REAL(dp), ALLOCATABLE                              :: periodic_atom_response(:, :)
      REAL(KIND=dp)                                      :: grid_base_weight, weight_grad
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: weights
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: atom_coords_pbc
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: dweights_datom, dweights_dstrain
      REAL(KIND=dp), DIMENSION(3)                        :: grid_point
      REAL(KIND=dp), DIMENSION(:), POINTER               :: atomic_grid_weight_grad, grid_weight_grad

      NULLIFY (atomic_grid_weight_grad, grid_weight_grad)
      CALL torch_tensor_grad(features%grid_weights_t, grid_weight_grad_t)
      CALL torch_tensor_grad(features%atomic_grid_weights_t, atomic_grid_weight_grad_t)
      CALL torch_tensor_data_ptr(grid_weight_grad_t, grid_weight_grad)
      CALL torch_tensor_data_ptr(atomic_grid_weight_grad_t, atomic_grid_weight_grad)
      IF (features%periodic_context) THEN
         ALLOCATE (periodic_atom_response(3, SIZE(particle_set)))
         CALL skala_gpw_periodic_partition_response(features, particle_set, cell, rho_r(1)%pw_grid, &
                                                    grid_weight_grad, atomic_grid_weight_grad, &
                                                    periodic_atom_response, periodic_strain_response)
         atom_force = atom_force + periodic_atom_response
         CALL torch_tensor_release(grid_weight_grad_t)
         CALL torch_tensor_release(atomic_grid_weight_grad_t)
         RETURN
      END IF

      natom = SIZE(particle_set)
      CPASSERT(SIZE(atom_force, 1) == 3)
      CPASSERT(SIZE(atom_force, 2) == natom)
      ALLOCATE (atom_coords_pbc(3, natom), included(natom), weights(natom), &
                dweights_datom(3, natom, natom), dweights_dstrain(3, 3, natom))
      DO iatom = 1, natom
         atom_coords_pbc(:, iatom) = pbc(particle_set(iatom)%r, cell, positive_range=.TRUE.)
      END DO

      bo = rho_r(1)%pw_grid%bounds_local
      local_row = 0
      DO k = bo(1, 3), bo(2, 3)
         DO j = bo(1, 2), bo(2, 2)
            DO i = bo(1, 1), bo(2, 1)
               local_row = local_row + 1
               grid_point = native_grid_coordinate(rho_r(1)%pw_grid, [i, j, k])
               CALL skala_gpw_smooth_partition_derivatives(grid_point, atom_coords_pbc, cell, &
                                                           weights, included, dweights_datom, &
                                                           dweights_dstrain)
               feature_begin = features%local_feature_offsets(local_row)
               feature_end = features%local_feature_offsets(local_row + 1) - 1
               CPASSERT(feature_end - feature_begin + 1 == COUNT(included))
               grid_base_weight = 0.0_dp
               DO feature_pos = feature_begin, feature_end
                  row = features%local_feature_rows(feature_pos)
                  grid_base_weight = grid_base_weight + features%grid_weights(row)
               END DO
               feature_pos = feature_begin
               DO iatom = 1, natom
                  IF (.NOT. included(iatom)) CYCLE
                  row = features%local_feature_rows(feature_pos)
                  weight_grad = grid_base_weight*grid_weight_grad(row) + &
                                rho_r(1)%pw_grid%dvol*atomic_grid_weight_grad(row)* &
                                smooth_partition_atomic_weight_scale_derivative(weights(iatom))
                  DO jatom = 1, natom
                     atom_force(:, jatom) = atom_force(:, jatom) + &
                                            weight_grad*dweights_datom(:, jatom, iatom)
                  END DO
                  feature_pos = feature_pos + 1
               END DO
               CPASSERT(feature_pos == feature_end + 1)
            END DO
         END DO
      END DO
      CPASSERT(local_row == features%nflat_local)

      DEALLOCATE (atom_coords_pbc, dweights_datom, dweights_dstrain, included, weights)
      CALL torch_tensor_release(grid_weight_grad_t)
      CALL torch_tensor_release(atomic_grid_weight_grad_t)

   END SUBROUTINE add_smooth_partition_force

! **************************************************************************************************
!> \brief Add the virial from SMOOTH native-grid atom partition weights.
!> \param virial_xc ...
!> \param features ...
!> \param particle_set ...
!> \param cell ...
!> \param rho_r ...
!> \param grid_weight_grad_t ...
!> \param atomic_grid_weight_grad_t ...
! **************************************************************************************************
   SUBROUTINE build_smooth_partition_virial(virial_xc, features, particle_set, cell, rho_r, &
                                            grid_weight_grad_t, atomic_grid_weight_grad_t)
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(INOUT)      :: virial_xc
      TYPE(skala_gpw_feature_type), INTENT(IN)           :: features
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: rho_r
      TYPE(torch_tensor_type), INTENT(INOUT)             :: grid_weight_grad_t, &
                                                            atomic_grid_weight_grad_t

      INTEGER                                            :: feature_begin, feature_end, feature_pos, &
                                                            i, iatom, idir, j, jdir, k, local_row, &
                                                            natom, row
      INTEGER, DIMENSION(2, 3)                           :: bo
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: included
      REAL(dp)                                           :: periodic_strain_response(3, 3)
      REAL(dp), ALLOCATABLE                              :: periodic_atom_response(:, :)
      REAL(KIND=dp)                                      :: grid_base_weight, tmp, weight_grad
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: weights
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: atom_coords_pbc
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: dweights_datom, dweights_dstrain
      REAL(KIND=dp), DIMENSION(3)                        :: grid_point
      REAL(KIND=dp), DIMENSION(:), POINTER               :: atomic_grid_weight_grad, grid_weight_grad

      NULLIFY (atomic_grid_weight_grad, grid_weight_grad)
      CALL torch_tensor_grad(features%grid_weights_t, grid_weight_grad_t)
      CALL torch_tensor_grad(features%atomic_grid_weights_t, atomic_grid_weight_grad_t)
      CALL torch_tensor_data_ptr(grid_weight_grad_t, grid_weight_grad)
      CALL torch_tensor_data_ptr(atomic_grid_weight_grad_t, atomic_grid_weight_grad)
      IF (features%periodic_context) THEN
         ALLOCATE (periodic_atom_response(3, SIZE(particle_set)))
         CALL skala_gpw_periodic_partition_response(features, particle_set, cell, rho_r(1)%pw_grid, &
                                                    grid_weight_grad, atomic_grid_weight_grad, &
                                                    periodic_atom_response, periodic_strain_response)
         virial_xc = virial_xc + periodic_strain_response
         CALL torch_tensor_release(grid_weight_grad_t)
         CALL torch_tensor_release(atomic_grid_weight_grad_t)
         RETURN
      END IF

      natom = SIZE(particle_set)
      ALLOCATE (atom_coords_pbc(3, natom), included(natom), weights(natom), &
                dweights_datom(3, natom, natom), dweights_dstrain(3, 3, natom))
      DO iatom = 1, natom
         atom_coords_pbc(:, iatom) = pbc(particle_set(iatom)%r, cell, positive_range=.TRUE.)
      END DO

      bo = rho_r(1)%pw_grid%bounds_local
      local_row = 0
      DO k = bo(1, 3), bo(2, 3)
         DO j = bo(1, 2), bo(2, 2)
            DO i = bo(1, 1), bo(2, 1)
               local_row = local_row + 1
               grid_point = native_grid_coordinate(rho_r(1)%pw_grid, [i, j, k])
               CALL skala_gpw_smooth_partition_derivatives(grid_point, atom_coords_pbc, cell, &
                                                           weights, included, dweights_datom, &
                                                           dweights_dstrain)
               feature_begin = features%local_feature_offsets(local_row)
               feature_end = features%local_feature_offsets(local_row + 1) - 1
               CPASSERT(feature_end - feature_begin + 1 == COUNT(included))
               grid_base_weight = 0.0_dp
               DO feature_pos = feature_begin, feature_end
                  row = features%local_feature_rows(feature_pos)
                  grid_base_weight = grid_base_weight + features%grid_weights(row)
               END DO
               feature_pos = feature_begin
               DO iatom = 1, natom
                  IF (.NOT. included(iatom)) CYCLE
                  row = features%local_feature_rows(feature_pos)
                  weight_grad = grid_base_weight*grid_weight_grad(row) + &
                                rho_r(1)%pw_grid%dvol*atomic_grid_weight_grad(row)* &
                                smooth_partition_atomic_weight_scale_derivative(weights(iatom))
                  DO idir = 1, 3
                     DO jdir = 1, idir
                        tmp = weight_grad*dweights_dstrain(idir, jdir, iatom)
                        virial_xc(jdir, idir) = virial_xc(jdir, idir) + tmp
                        IF (idir /= jdir) virial_xc(idir, jdir) = virial_xc(idir, jdir) + tmp
                     END DO
                  END DO
                  feature_pos = feature_pos + 1
               END DO
               CPASSERT(feature_pos == feature_end + 1)
            END DO
         END DO
      END DO
      CPASSERT(local_row == features%nflat_local)

      DEALLOCATE (atom_coords_pbc, dweights_datom, dweights_dstrain, included, weights)
      CALL torch_tensor_release(grid_weight_grad_t)
      CALL torch_tensor_release(atomic_grid_weight_grad_t)

   END SUBROUTINE build_smooth_partition_virial

! **************************************************************************************************
!> \brief Return the Cartesian coordinate of a regular GPW grid point.
!> \param pw_grid ...
!> \param index ...
!> \return ...
! **************************************************************************************************
   FUNCTION native_grid_coordinate(pw_grid, index) RESULT(coord)
      TYPE(pw_grid_type), POINTER                        :: pw_grid
      INTEGER, DIMENSION(3), INTENT(IN)                  :: index
      REAL(KIND=dp), DIMENSION(3)                        :: coord

      INTEGER, DIMENSION(3)                              :: relative_index

      relative_index = index - pw_grid%bounds(1, :)
      coord = REAL(relative_index(1), KIND=dp)*pw_grid%dh(:, 1) + &
              REAL(relative_index(2), KIND=dp)*pw_grid%dh(:, 2) + &
              REAL(relative_index(3), KIND=dp)*pw_grid%dh(:, 3)

   END FUNCTION native_grid_coordinate

! **************************************************************************************************
!> \brief Evaluate a rank-local atom chunk as multiple atom-contiguous Torch subchunks.
!> \param features ...
!> \param group ...
!> \param max_rows ...
!> \param compute_grads ...
!> \param exc ...
!> \param density_grad ...
!> \param grad_grad ...
!> \param kin_grad ...
!> \param collapse_spin_grads ...
!> \param density_direction Optional matching full-chunk primitive direction.
!> \param grad_direction ...
!> \param kin_direction ...
!> \param density_response New local-grid weighted covector response.
!> \param grad_response ...
!> \param kin_response ...
! **************************************************************************************************
   SUBROUTINE evaluate_atom_subchunks(features, group, max_rows, compute_grads, exc, &
                                      density_grad, grad_grad, kin_grad, collapse_spin_grads, &
                                      density_direction, grad_direction, kin_direction, &
                                      density_response, grad_response, kin_response)
      TYPE(skala_gpw_feature_type), INTENT(IN)           :: features

      CLASS(mp_comm_type), INTENT(IN)                    :: group
      INTEGER, INTENT(IN)                                :: max_rows
      LOGICAL, INTENT(IN)                                :: compute_grads, collapse_spin_grads
      REAL(KIND=dp), INTENT(OUT)                         :: exc
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: density_grad, kin_grad
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(OUT)                                     :: grad_grad

      REAL(dp), INTENT(IN), OPTIONAL                     :: density_direction(:, :), grad_direction(:, :, :), &
                                                            kin_direction(:, :)
      REAL(dp), ALLOCATABLE, INTENT(OUT), OPTIONAL        :: density_response(:, :), grad_response(:, :, :), &
                                                             kin_response(:, :)
      REAL(dp), ALLOCATABLE                              :: full_density_response(:, :), full_grad_response(:, :, :), &
                                                            full_kin_response(:, :)
      REAL(dp), POINTER                                  :: density_view(:, :), grad_view(:, :, :), kin_view(:, :)
      TYPE(torch_tensor_type)                            :: response_tensors(3)
      INTEGER                                            :: component, first, last, ns
      LOGICAL                                            :: need_response

      INTEGER                                            :: isubchunk, nroute_grad_per_point, &
                                                            nroute_recv_points, nroute_send_points, &
                                                            nsubchunks, phase_handle, subphase_handle
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: route_grad_return_recv_counts, &
                                                            route_grad_return_recv_displs, &
                                                            route_grad_return_send_counts, &
                                                            route_grad_return_send_displs, &
                                                            subchunk_atom_begin, &
                                                            subchunk_atom_count, &
                                                            subchunk_row_begin, &
                                                            subchunk_row_count
      REAL(KIND=dp)                                      :: subchunk_exc
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: recv_grad_buffer, send_grad_buffer
      TYPE(skala_gpw_feature_type)                       :: subchunk
      TYPE(torch_tensor_type)                            :: subchunk_exc_tensor

      CPASSERT(features%uses_atom_chunks)
      CPASSERT(max_rows > 0)
      CALL skala_gpw_atom_subchunk_layout(max_rows, subchunk_atom_begin, subchunk_atom_count, &
                                          subchunk_row_begin, subchunk_row_count)
      nsubchunks = SIZE(subchunk_atom_begin)
      need_response = PRESENT(density_response)
      CPASSERT(PRESENT(grad_response) .EQV. need_response)
      CPASSERT(PRESENT(kin_response) .EQV. need_response)
      CPASSERT(PRESENT(density_direction) .EQV. need_response)
      CPASSERT(PRESENT(grad_direction) .EQV. need_response)
      CPASSERT(PRESENT(kin_direction) .EQV. need_response)
      IF (need_response) THEN
         CPASSERT(compute_grads)
         ns = 2
         IF (features%uses_collapsed_rks_dynamic) ns = 1
         CPASSERT(ALL(SHAPE(density_direction) == [features%chunk_feature_count, ns]))
         CPASSERT(ALL(SHAPE(grad_direction) == [features%chunk_feature_count, 3, ns]))
         CPASSERT(ALL(SHAPE(kin_direction) == [features%chunk_feature_count, ns]))
         ALLOCATE (full_density_response(features%chunk_feature_count, ns), &
                   full_grad_response(features%chunk_feature_count, 3, ns), &
                   full_kin_response(features%chunk_feature_count, ns))
         full_density_response = 0.0_dp
         full_grad_response = 0.0_dp
         full_kin_response = 0.0_dp
      END IF

      exc = 0.0_dp
      IF (compute_grads) THEN
         CPASSERT(features%uses_atom_chunk_routing)
         nroute_recv_points = SUM(features%route_point_recv_counts)
         nroute_send_points = SIZE(features%route_send_local_rows)
         CPASSERT(SUM(features%route_point_send_counts) == nroute_send_points)
         CPASSERT(SIZE(features%route_chunk_offsets) == nroute_recv_points + 1)
         nroute_grad_per_point = ngrad_per_point
         IF (collapse_spin_grads) nroute_grad_per_point = ncollapsed_grad_per_point
         ALLOCATE (send_grad_buffer(MAX(1, nroute_grad_per_point*nroute_recv_points)), &
                   recv_grad_buffer(MAX(1, nroute_grad_per_point*nroute_send_points)), &
                   route_grad_return_send_counts(SIZE(features%route_point_recv_counts)), &
                   route_grad_return_send_displs(SIZE(features%route_point_recv_displs)), &
                   route_grad_return_recv_counts(SIZE(features%route_point_send_counts)), &
                   route_grad_return_recv_displs(SIZE(features%route_point_send_displs)))
         route_grad_return_send_counts(:) = &
            nroute_grad_per_point*features%route_point_recv_counts
         route_grad_return_send_displs(:) = &
            nroute_grad_per_point*features%route_point_recv_displs
         route_grad_return_recv_counts(:) = &
            nroute_grad_per_point*features%route_point_send_counts
         route_grad_return_recv_displs(:) = &
            nroute_grad_per_point*features%route_point_send_displs
      END IF

      CALL timeset("skala_gpw_atom_subchunks", phase_handle)
      DO isubchunk = 1, nsubchunks
         CALL timeset("skala_gpw_atom_subchunk_build", subphase_handle)
         CALL skala_gpw_feature_build_atom_subchunk_bounds(features, subchunk, &
                                                           subchunk_atom_begin(isubchunk), &
                                                           subchunk_atom_count(isubchunk), &
                                                           subchunk_row_begin(isubchunk), &
                                                           subchunk_row_count(isubchunk), &
                                                           compute_grads)
         CALL timestop(subphase_handle)
         CALL timeset("skala_gpw_atom_subchunk_forward", subphase_handle)
         CALL skala_torch_model_get_exc(cached_model, subchunk%inputs, &
                                        subchunk%grid_weights_t, subchunk_exc_tensor, &
                                        subchunk_exc)
         CALL timestop(subphase_handle)
         exc = exc + subchunk_exc
         IF (need_response) THEN
            first = subchunk_row_begin(isubchunk)
            last = first + subchunk_row_count(isubchunk) - 1
            CALL feature_hessian_product(subchunk, subchunk_exc_tensor, density_direction(first:last, :), &
                                         grad_direction(first:last, :, :), kin_direction(first:last, :), response_tensors)
            NULLIFY (density_view, grad_view, kin_view)
            CALL torch_tensor_data_ptr(response_tensors(1), density_view)
            CALL torch_tensor_data_ptr(response_tensors(2), grad_view)
            CALL torch_tensor_data_ptr(response_tensors(3), kin_view)
            full_density_response(first:last, :) = density_view
            full_grad_response(first:last, :, :) = grad_view
            full_kin_response(first:last, :) = kin_view
            NULLIFY (density_view, grad_view, kin_view)
            DO component = 1, 3
               CALL torch_tensor_release(response_tensors(component))
            END DO
         END IF
         IF (compute_grads) THEN
            CALL timeset("skala_gpw_atom_subchunk_backward", subphase_handle)
            CALL torch_tensor_backward_scalar(subchunk_exc_tensor)
            CALL timestop(subphase_handle)
         END IF
         CALL timeset("skala_gpw_atom_subchunk_release", subphase_handle)
         CALL torch_tensor_release(subchunk_exc_tensor)
         CALL skala_gpw_feature_release(subchunk)
         CALL timestop(subphase_handle)
      END DO
      IF (compute_grads .AND. features%chunk_feature_count > 0) THEN
         CALL timeset("skala_gpw_atom_subchunk_grad_pack", subphase_handle)
         CALL pack_atom_chunk_grads(features, send_grad_buffer, .TRUE., collapse_spin_grads)
         CALL timestop(subphase_handle)
      END IF
      CALL timestop(phase_handle)

      IF (compute_grads) THEN
         CALL timeset("skala_gpw_grad_route_comm", phase_handle)
         CALL group%alltoall(send_grad_buffer, route_grad_return_send_counts, &
                             route_grad_return_send_displs, recv_grad_buffer, &
                             route_grad_return_recv_counts, route_grad_return_recv_displs)
         CALL timestop(phase_handle)

         CALL timeset("skala_gpw_grad_route_scatter", phase_handle)
         CALL scatter_routed_atom_chunk_grads(features, recv_grad_buffer, collapse_spin_grads, &
                                              density_grad, grad_grad, kin_grad)
         CALL timestop(phase_handle)

         DEALLOCATE (recv_grad_buffer, route_grad_return_recv_counts, &
                     route_grad_return_recv_displs, route_grad_return_send_counts, &
                     route_grad_return_send_displs, send_grad_buffer)
      END IF
      IF (need_response) THEN
         IF (features%chunk_feature_count > 0) THEN
            CALL torch_tensor_from_array(response_tensors(1), full_density_response)
            CALL torch_tensor_from_array(response_tensors(2), full_grad_response)
            CALL torch_tensor_from_array(response_tensors(3), full_kin_response)
         END IF
         CALL fetch_and_gather_atom_chunk_grads(features, group, density_response, grad_response, kin_response, &
                                                response_tensors)
         IF (features%chunk_feature_count > 0) THEN
            DO component = 1, 3
               CALL torch_tensor_release(response_tensors(component))
            END DO
         END IF
      END IF
      DEALLOCATE (subchunk_atom_begin, subchunk_atom_count, subchunk_row_begin, subchunk_row_count)

   END SUBROUTINE evaluate_atom_subchunks

! **************************************************************************************************
!> \brief Joint model Hessian product with independent directions and retained ordinary gradients.
!> \param features ...
!> \param energy Weighted original scalar model energy.
!> \param density_direction ...
!> \param grad_direction ...
!> \param kin_direction ...
!> \param response Independent CPU response tensors, released by the caller.
! **************************************************************************************************
   SUBROUTINE feature_hessian_product(features, energy, density_direction, grad_direction, kin_direction, response)
      TYPE(skala_gpw_feature_type), INTENT(IN)           :: features
      TYPE(torch_tensor_type), INTENT(IN)                :: energy
      REAL(dp), INTENT(IN)                               :: density_direction(:, :), &
                                                            grad_direction(:, :, :), &
                                                            kin_direction(:, :)
      TYPE(torch_tensor_type), INTENT(INOUT)             :: response(3)

      INTEGER                                            :: component
      REAL(dp), ALLOCATABLE                              :: density(:, :), gradient(:, :, :), &
                                                            kinetic(:, :)
      TYPE(torch_tensor_type)                            :: directions(3), fields(3)

      fields = [features%density_t, features%grad_t, features%kin_t]
      ALLOCATE (density, SOURCE=density_direction)
      ALLOCATE (gradient, SOURCE=grad_direction)
      ALLOCATE (kinetic, SOURCE=kin_direction)
      CALL torch_tensor_from_array(directions(1), density)
      CALL torch_tensor_from_array(directions(2), gradient)
      CALL torch_tensor_from_array(directions(3), kinetic)
      CALL torch_tensor_hessian_product3(energy, fields, directions, response)
      DO component = 1, 3
         CALL torch_tensor_release(directions(component))
      END DO
   END SUBROUTINE feature_hessian_product

! **************************************************************************************************
!> \brief Select an automatic atom-subchunk row cap.
!> \param features ...
!> \param group ...
!> \return ...
! **************************************************************************************************
   FUNCTION auto_atom_chunk_max_rows(features, group) RESULT(max_rows)
      TYPE(skala_gpw_feature_type), INTENT(IN)           :: features

      CLASS(mp_comm_type), INTENT(IN)                    :: group
      INTEGER                                            :: max_rows

      INTEGER                                            :: local_rows_max, target_rows

      local_rows_max = features%chunk_feature_count
      CALL group%max(local_rows_max)
      IF (local_rows_max <= 0) THEN
         max_rows = 0
         RETURN
      END IF

      IF (group%num_pe > 1) THEN
         target_rows = CEILING(REAL(local_rows_max, KIND=dp)/2.0_dp)
         max_rows = atom_chunk_auto_row_quantum* &
                    ((target_rows + atom_chunk_auto_row_quantum - 1)/atom_chunk_auto_row_quantum)
      ELSE
         target_rows = NINT(REAL(local_rows_max, KIND=dp)/4.0_dp)
         max_rows = atom_chunk_auto_row_quantum* &
                    MAX(1, NINT(REAL(target_rows, KIND=dp)/ &
                                REAL(atom_chunk_auto_row_quantum, KIND=dp)))
      END IF
      max_rows = MAX(atom_chunk_auto_min_rows, MIN(atom_chunk_auto_max_rows, max_rows))

   END FUNCTION auto_atom_chunk_max_rows

! **************************************************************************************************
!> \brief Check the geometry associated with a selected density subgradient.
!> \param state ...
!> \param rho_r ...
!> \param particle_set ...
!> \param cell ...
!> \return ...
! **************************************************************************************************
   FUNCTION density_constraints_match_geometry(state, rho_r, particle_set, cell) RESULT(matches)
      TYPE(skala_density_constraints_type), INTENT(IN)   :: state
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: rho_r
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(cell_type), POINTER                           :: cell
      LOGICAL                                            :: matches

      INTEGER                                            :: i
      REAL(dp), ALLOCATABLE                              :: positions(:, :)

      matches = .FALSE.
      IF (skala_density_constraints_count(state, skala_fft_domain) == 0) RETURN
      ALLOCATE (positions(3, SIZE(particle_set)))
      DO i = 1, SIZE(particle_set)
         positions(:, i) = particle_set(i)%r
      END DO
      matches = skala_density_constraints_matches(state, SIZE(rho_r), rho_r(1)%pw_grid%npts, &
                                                  rho_r(1)%pw_grid%bounds, cell%hmat, positions)
   END FUNCTION density_constraints_match_geometry

! **************************************************************************************************
!> \brief Record the same local spin inputs packed by skala_gpw_feature_build.
!>        The XC rho set already includes the caller's core and density processing.
!>        Samples precede atom duplication, routing and nonlinear model features.
!> \param state ...
!> \param rho_set ...
!> \param rho_r ...
! **************************************************************************************************
   SUBROUTINE record_fft_descriptors(state, rho_set, rho_r)
      TYPE(skala_density_constraints_type), &
         INTENT(INOUT)                                   :: state
      TYPE(xc_rho_set_type), INTENT(IN)                  :: rho_set
      TYPE(pw_r3d_rs_type), POINTER                      :: rho_r(:)

      INTEGER                                            :: bo(2, 3), i, j, k, np, ns, row
      INTEGER, ALLOCATABLE                               :: points(:, :)
      REAL(dp), ALLOCATABLE                              :: density(:, :), kinetic(:, :)
      REAL(dp), POINTER                                  :: rho(:, :, :), rhoa(:, :, :), &
                                                            rhob(:, :, :), tau(:, :, :), &
                                                            tau_a(:, :, :), tau_b(:, :, :)

      IF (.NOT. skala_density_samples_active(state)) RETURN
      ns = SIZE(rho_r)
      bo = rho_r(1)%pw_grid%bounds_local
      np = PRODUCT(MAX(0, bo(2, :) - bo(1, :) + 1))
      ALLOCATE (points(3, np), density(np, ns), kinetic(np, ns))
      IF (ns == 1) THEN
         CALL xc_rho_set_get(rho_set, rho=rho, tau=tau)
      ELSE
         CPASSERT(ns == 2)
         CALL xc_rho_set_get(rho_set, rhoa=rhoa, rhob=rhob, tau_a=tau_a, tau_b=tau_b)
      END IF
      row = 0
      DO k = bo(1, 3), bo(2, 3)
         DO j = bo(1, 2), bo(2, 2)
            DO i = bo(1, 1), bo(2, 1)
               row = row + 1
               points(:, row) = [i, j, k]
               IF (ns == 1) THEN
                  density(row, 1) = 0.5_dp*rho(i, j, k)
                  kinetic(row, 1) = 0.5_dp*tau(i, j, k)
               ELSE
                  density(row, :) = [rhoa(i, j, k), rhob(i, j, k)]
                  kinetic(row, :) = [tau_a(i, j, k), tau_b(i, j, k)]
               END IF
            END DO
         END DO
      END DO
      CALL skala_density_samples_put(state, skala_fft_domain, points, density, kinetic)
   END SUBROUTINE record_fft_descriptors

! **************************************************************************************************
!> \brief Apply selected subgradients at registered physical density and kinetic-density targets.
!>        Every duplicated atomic feature row at an RKS point uses the same weight.
!>        All other covectors and the original energy remain unchanged.
!> \param state ...
!> \param features ...
!> \param rho_r ...
!> \param density_grad ...
!> \param kin_grad ...
! **************************************************************************************************
   SUBROUTINE apply_descriptor_constraints(state, features, rho_r, density_grad, kin_grad)
      TYPE(skala_density_constraints_type), &
         INTENT(INOUT)                                   :: state
      TYPE(skala_gpw_feature_type), INTENT(IN)           :: features
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: rho_r
      REAL(dp), INTENT(INOUT), TARGET                    :: density_grad(:, :), kin_grad(:, :)

      INTEGER                                            :: bo(2, 3), c, descriptor, e, first_spin, &
                                                            i, last_spin, n, nactive, nlocal, pos, &
                                                            row, spin
      INTEGER, ALLOCATABLE :: active(:, :), counts(:), descriptors(:), displs(:), &
         global_entries(:), local_entries(:), local_entry_counts(:), local_rows(:), owners(:), &
         points(:, :), spins(:)
      LOGICAL                                            :: branches(2), capture_jumps
      LOGICAL, ALLOCATABLE                               :: jump_owned(:)
      REAL(dp)                                           :: probes(2), tolerances(2)
      REAL(dp), ALLOCATABLE :: active_targets(:), constraint_error(:), jump_values(:), &
         negative(:, :), negative_weights(:), positive(:, :), targets(:)
      REAL(dp), POINTER                                  :: covector(:, :)

      CPASSERT(.NOT. features%uses_atom_chunks)
      CALL skala_density_constraints_get(state, points, spins, negative_weights, probes(1), tolerances(1), &
                                         descriptors, probes(2), tolerances(2), domain=skala_fft_domain, targets=targets)
      n = SIZE(spins)
      capture_jumps = skala_density_jumps_active(state)
      IF (capture_jumps) THEN
         ALLOCATE (jump_values(n), jump_owned(n))
         jump_values = 0.0_dp
         jump_owned = .FALSE.
      END IF
      bo = rho_r(1)%pw_grid%bounds_local
      ALLOCATE (local_rows(n), owners(n), local_entry_counts(n), constraint_error(n))
      local_rows = 0
      owners = 0
      local_entry_counts = 0
      constraint_error = 0.0_dp
      DO c = 1, n
         IF (ANY(points(:, c) < bo(1, :)) .OR. ANY(points(:, c) > bo(2, :))) CYCLE
         owners(c) = 1
         ! The feature index refers to a global atomic feature, not the local FFT point.
         row = points(1, c) - bo(1, 1) + 1 + (bo(2, 1) - bo(1, 1) + 1)* &
               (points(2, c) - bo(1, 2) + (bo(2, 2) - bo(1, 2) + 1)*(points(3, c) - bo(1, 3)))
         CPASSERT(row >= 1 .AND. row <= features%nflat_local)
         local_rows(c) = row
         first_spin = spins(c)
         last_spin = spins(c)
         IF (spins(c) == 0) THEN
            first_spin = 1
            last_spin = 2
         END IF
         DO pos = features%local_feature_offsets(row), features%local_feature_offsets(row + 1) - 1
            i = features%local_feature_rows(pos)
            IF (descriptors(c) == skala_density_descriptor) THEN
               constraint_error(c) = MAX(constraint_error(c), MAXVAL(ABS(features%density(i, first_spin:last_spin) - targets(c))))
            ELSE
               constraint_error(c) = MAX(constraint_error(c), MAXVAL(ABS(features%kin(i, first_spin:last_spin) - targets(c))))
            END IF
            local_entry_counts(c) = local_entry_counts(c) + 4*(last_spin - first_spin + 1)
         END DO
      END DO
      CALL rho_r(1)%pw_grid%para%group%sum(owners)
      CPASSERT(ALL(owners == 1))
      CALL rho_r(1)%pw_grid%para%group%max(constraint_error)
      ! Away from the selected target the original differentiable branch supplies the derivative.
      IF (ALL(constraint_error > tolerances(descriptors))) RETURN
      nlocal = SUM(local_entry_counts, MASK=constraint_error <= tolerances(descriptors))

      ALLOCATE (local_entries(nlocal), counts(rho_r(1)%pw_grid%para%group%num_pe), &
                displs(rho_r(1)%pw_grid%para%group%num_pe))
      e = 0
      DO c = 1, n
         IF (constraint_error(c) > tolerances(descriptors(c))) CYCLE
         row = local_rows(c)
         IF (row == 0) CYCLE
         first_spin = MAX(1, spins(c))
         last_spin = spins(c)
         IF (spins(c) == 0) last_spin = 2
         DO pos = features%local_feature_offsets(row), features%local_feature_offsets(row + 1) - 1
            DO spin = first_spin, last_spin
               local_entries(e + 1:e + 4) = [features%local_feature_rows(pos), spin, c, descriptors(c)]
               e = e + 4
            END DO
         END DO
      END DO
      CPASSERT(e == nlocal)
      CALL rho_r(1)%pw_grid%para%group%allgather(nlocal, counts)
      IF (SUM(counts) == 0) RETURN
      displs(1) = 0
      DO i = 2, SIZE(counts)
         displs(i) = displs(i - 1) + counts(i - 1)
      END DO
      ALLOCATE (global_entries(SUM(counts)))
      CALL rho_r(1)%pw_grid%para%group%allgatherv(local_entries, global_entries, counts, displs)
      DO descriptor = skala_density_descriptor, skala_kinetic_descriptor
         nactive = COUNT(global_entries(4::4) == descriptor)
         IF (nactive == 0) CYCLE
         ALLOCATE (active(2, nactive), active_targets(nactive))
         e = 0
         DO i = 1, SIZE(global_entries)/4
            IF (global_entries(4*i) /= descriptor) CYCLE
            e = e + 1
            active(:, e) = global_entries(4*i - 3:4*i - 2)
            active_targets(e) = targets(global_entries(4*i - 1))
         END DO
         branches = .FALSE.
         DO i = 1, SIZE(global_entries)/4
            IF (global_entries(4*i) /= descriptor) CYCLE
            c = global_entries(4*i - 1)
            branches(1) = branches(1) .OR. negative_weights(c) /= 0.0_dp
            branches(2) = branches(2) .OR. negative_weights(c) /= 1.0_dp
         END DO
         IF (capture_jumps) branches = .TRUE.
         IF (descriptor == skala_density_descriptor) THEN
            CALL skala_torch_model_descriptor_limits(cached_model, features%inputs, "density", features%density, &
                                                     active, probes(descriptor), negative, positive, branches=branches, &
                                                     targets=active_targets, adjacent=cached_model_fp64)
            covector => density_grad
         ELSE
            CALL skala_torch_model_descriptor_limits(cached_model, features%inputs, "kin", features%kin, &
                                                     active, probes(descriptor), negative, positive, branches=branches, &
                                                     targets=active_targets, adjacent=cached_model_fp64)
            covector => kin_grad
         END IF
         DO c = 1, n
            IF (descriptors(c) /= descriptor) CYCLE
            IF (constraint_error(c) > tolerances(descriptor)) CYCLE
            row = local_rows(c)
            IF (row == 0) CYCLE
            first_spin = MAX(1, spins(c))
            last_spin = spins(c)
            IF (spins(c) == 0) last_spin = 2
            covector(row, first_spin:last_spin) = 0.0_dp
            DO pos = features%local_feature_offsets(row), features%local_feature_offsets(row + 1) - 1
               i = features%local_feature_rows(pos)
               covector(row, first_spin:last_spin) = covector(row, first_spin:last_spin) + &
                                                     negative_weights(c)*negative(i, first_spin:last_spin) + &
                                                     (1.0_dp - negative_weights(c))*positive(i, first_spin:last_spin)
               IF (capture_jumps) THEN
                  jump_values(c) = jump_values(c) + &
                                   SUM(negative(i, first_spin:last_spin) - positive(i, first_spin:last_spin))
               END IF
            END DO
            IF (capture_jumps) jump_owned(c) = .TRUE.
         END DO
         NULLIFY (covector)
         DEALLOCATE (active, active_targets)
      END DO
      IF (capture_jumps) CALL skala_density_jumps_put(state, skala_fft_domain, jump_values, jump_owned)
   END SUBROUTINE apply_descriptor_constraints

! **************************************************************************************************
!> \brief Map full Torch feature gradients back to this rank's local grid order.
!> \param features ...
!> \param density_grad ...
!> \param grad_grad ...
!> \param kin_grad ...
!> \param field_response Optional borrowed weighted response tensors in the same feature layout.
! **************************************************************************************************
   SUBROUTINE fetch_local_feature_grads(features, density_grad, grad_grad, kin_grad, field_response)
      TYPE(skala_gpw_feature_type), INTENT(IN)           :: features
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: density_grad
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(OUT)                                     :: grad_grad
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: kin_grad
      TYPE(torch_tensor_type), INTENT(IN), OPTIONAL      :: field_response(3)

      INTEGER                                            :: feature_pos, i, j, k, local_row, row
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: density_grad_all, kin_grad_all
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: grad_grad_all
      TYPE(torch_tensor_type)                            :: density_grad_t, grad_grad_t, kin_grad_t

      NULLIFY (density_grad_all, grad_grad_all, kin_grad_all)
      CALL get_feature_grad_views(features, density_grad_t, grad_grad_t, kin_grad_t, &
                                  density_grad_all, grad_grad_all, kin_grad_all, field_response)
      CPASSERT(SIZE(density_grad_all, 1) == features%nflat)
      CPASSERT(SIZE(density_grad_all, 2) == 2)
      CPASSERT(SIZE(grad_grad_all, 1) == features%nflat)
      CPASSERT(SIZE(grad_grad_all, 2) == 3)
      CPASSERT(SIZE(grad_grad_all, 3) == 2)
      CPASSERT(SIZE(kin_grad_all, 1) == features%nflat)
      CPASSERT(SIZE(kin_grad_all, 2) == 2)

      ALLOCATE (density_grad(features%nflat_local, 2), &
                grad_grad(features%nflat_local, 3, 2), &
                kin_grad(features%nflat_local, 2))
      density_grad = 0.0_dp
      grad_grad = 0.0_dp
      kin_grad = 0.0_dp
      local_row = 0
      DO k = LBOUND(features%feature_index, 3), UBOUND(features%feature_index, 3)
         DO j = LBOUND(features%feature_index, 2), UBOUND(features%feature_index, 2)
            DO i = LBOUND(features%feature_index, 1), UBOUND(features%feature_index, 1)
               local_row = local_row + 1
               DO feature_pos = features%local_feature_offsets(local_row), &
                  features%local_feature_offsets(local_row + 1) - 1
                  row = features%local_feature_rows(feature_pos)
                  CPASSERT(row >= 1 .AND. row <= features%nflat)
                  density_grad(local_row, :) = density_grad(local_row, :) + &
                                               density_grad_all(row, :)
                  grad_grad(local_row, :, :) = grad_grad(local_row, :, :) + &
                                               grad_grad_all(row, :, :)
                  kin_grad(local_row, :) = kin_grad(local_row, :) + kin_grad_all(row, :)
               END DO
            END DO
         END DO
      END DO
      CPASSERT(local_row == features%nflat_local)

      IF (.NOT. PRESENT(field_response)) THEN
         CALL torch_tensor_release(density_grad_t)
         CALL torch_tensor_release(grad_grad_t)
         CALL torch_tensor_release(kin_grad_t)
      END IF

   END SUBROUTINE fetch_local_feature_grads

! **************************************************************************************************
!> \brief Pack atom-chunk Torch gradients into CP2K communication buffers.
!> \param features ...
!> \param TARGET ...
!> \param route_to_return_positions ...
!> \param collapse_spin_grads ...
!> \param field_response Optional borrowed weighted response tensors in the same feature layout.
! **************************************************************************************************
   SUBROUTINE pack_atom_chunk_grads(features, TARGET, route_to_return_positions, &
                                    collapse_spin_grads, field_response)
      TYPE(skala_gpw_feature_type), INTENT(IN)           :: features
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(INOUT)                                   :: target
      LOGICAL, INTENT(IN)                                :: route_to_return_positions
      LOGICAL, INTENT(IN), OPTIONAL                      :: collapse_spin_grads
      TYPE(torch_tensor_type), INTENT(IN), OPTIONAL      :: field_response(3)

      INTEGER                                            :: base, feature_pos, irow, &
                                                            ngrad_buffer_per_point, point_pos, &
                                                            target_points
      LOGICAL                                            :: my_collapse_spin_grads
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: chunk_density_grad, chunk_kin_grad
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: chunk_grad_grad
      TYPE(torch_tensor_type)                            :: density_grad_t, grad_grad_t, kin_grad_t

      my_collapse_spin_grads = .FALSE.
      IF (PRESENT(collapse_spin_grads)) my_collapse_spin_grads = collapse_spin_grads
      ngrad_buffer_per_point = ngrad_per_point
      IF (my_collapse_spin_grads) ngrad_buffer_per_point = ncollapsed_grad_per_point

      NULLIFY (chunk_density_grad, chunk_grad_grad, chunk_kin_grad)
      CALL get_feature_grad_views(features, density_grad_t, grad_grad_t, kin_grad_t, &
                                  chunk_density_grad, chunk_grad_grad, chunk_kin_grad, field_response)
      CPASSERT(MOD(SIZE(TARGET), ngrad_buffer_per_point) == 0)
      target_points = SIZE(TARGET)/ngrad_buffer_per_point
      CPASSERT(SIZE(chunk_density_grad, 1) == features%chunk_feature_count)
      CPASSERT(SIZE(chunk_grad_grad, 1) == features%chunk_feature_count)
      CPASSERT(SIZE(chunk_grad_grad, 2) == 3)
      CPASSERT(SIZE(chunk_kin_grad, 1) == features%chunk_feature_count)
      IF (features%uses_collapsed_rks_dynamic) THEN
         CPASSERT(my_collapse_spin_grads)
         CPASSERT(SIZE(chunk_density_grad, 2) == 1)
         CPASSERT(SIZE(chunk_grad_grad, 3) == 1)
         CPASSERT(SIZE(chunk_kin_grad, 2) == 1)
      ELSE
         CPASSERT(SIZE(chunk_density_grad, 2) == 2)
         CPASSERT(SIZE(chunk_grad_grad, 3) == 2)
         CPASSERT(SIZE(chunk_kin_grad, 2) == 2)
      END IF

      IF (route_to_return_positions) THEN
         CPASSERT(target_points == SIZE(features%route_chunk_offsets) - 1)
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP SHARED(chunk_density_grad, chunk_grad_grad, chunk_kin_grad, features, &
!$OMP        my_collapse_spin_grads, ngrad_buffer_per_point, TARGET, target_points) &
!$OMP PRIVATE(base, feature_pos, irow, point_pos)
         DO point_pos = 1, target_points
            base = ngrad_buffer_per_point*(point_pos - 1)
            TARGET(base + 1:base + ngrad_buffer_per_point) = 0.0_dp
            DO feature_pos = features%route_chunk_offsets(point_pos), &
               features%route_chunk_offsets(point_pos + 1) - 1
               irow = features%route_chunk_rows(feature_pos)
               CPASSERT(irow >= 1 .AND. irow <= features%chunk_feature_count)
               IF (my_collapse_spin_grads) THEN
                  IF (features%uses_collapsed_rks_dynamic) THEN
                     TARGET(base + 1) = TARGET(base + 1) + &
                                        0.5_dp*chunk_density_grad(irow, 1)
                     TARGET(base + 2) = TARGET(base + 2) + &
                                        0.5_dp*chunk_grad_grad(irow, 1, 1)
                     TARGET(base + 3) = TARGET(base + 3) + &
                                        0.5_dp*chunk_grad_grad(irow, 2, 1)
                     TARGET(base + 4) = TARGET(base + 4) + &
                                        0.5_dp*chunk_grad_grad(irow, 3, 1)
                     TARGET(base + 5) = TARGET(base + 5) + &
                                        0.5_dp*chunk_kin_grad(irow, 1)
                  ELSE
                     TARGET(base + 1) = TARGET(base + 1) + &
                                        0.5_dp*(chunk_density_grad(irow, 1) + &
                                                chunk_density_grad(irow, 2))
                     TARGET(base + 2) = TARGET(base + 2) + &
                                        0.5_dp*(chunk_grad_grad(irow, 1, 1) + &
                                                chunk_grad_grad(irow, 1, 2))
                     TARGET(base + 3) = TARGET(base + 3) + &
                                        0.5_dp*(chunk_grad_grad(irow, 2, 1) + &
                                                chunk_grad_grad(irow, 2, 2))
                     TARGET(base + 4) = TARGET(base + 4) + &
                                        0.5_dp*(chunk_grad_grad(irow, 3, 1) + &
                                                chunk_grad_grad(irow, 3, 2))
                     TARGET(base + 5) = TARGET(base + 5) + &
                                        0.5_dp*(chunk_kin_grad(irow, 1) + &
                                                chunk_kin_grad(irow, 2))
                  END IF
               ELSE
                  TARGET(base + 1:base + 2) = TARGET(base + 1:base + 2) + &
                                              chunk_density_grad(irow, :)
                  TARGET(base + 3) = TARGET(base + 3) + chunk_grad_grad(irow, 1, 1)
                  TARGET(base + 4) = TARGET(base + 4) + chunk_grad_grad(irow, 2, 1)
                  TARGET(base + 5) = TARGET(base + 5) + chunk_grad_grad(irow, 3, 1)
                  TARGET(base + 6) = TARGET(base + 6) + chunk_grad_grad(irow, 1, 2)
                  TARGET(base + 7) = TARGET(base + 7) + chunk_grad_grad(irow, 2, 2)
                  TARGET(base + 8) = TARGET(base + 8) + chunk_grad_grad(irow, 3, 2)
                  TARGET(base + 9:base + 10) = TARGET(base + 9:base + 10) + &
                                               chunk_kin_grad(irow, :)
               END IF
            END DO
         END DO
!$OMP END PARALLEL DO
      ELSE
         CPASSERT(target_points >= features%chunk_feature_count)
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP SHARED(chunk_density_grad, chunk_grad_grad, chunk_kin_grad, features, &
!$OMP        my_collapse_spin_grads, ngrad_buffer_per_point, TARGET) &
!$OMP PRIVATE(base, irow)
         DO irow = 1, features%chunk_feature_count
            base = ngrad_buffer_per_point*(irow - 1)
            IF (my_collapse_spin_grads) THEN
               IF (features%uses_collapsed_rks_dynamic) THEN
                  TARGET(base + 1) = 0.5_dp*chunk_density_grad(irow, 1)
                  TARGET(base + 2) = 0.5_dp*chunk_grad_grad(irow, 1, 1)
                  TARGET(base + 3) = 0.5_dp*chunk_grad_grad(irow, 2, 1)
                  TARGET(base + 4) = 0.5_dp*chunk_grad_grad(irow, 3, 1)
                  TARGET(base + 5) = 0.5_dp*chunk_kin_grad(irow, 1)
               ELSE
                  TARGET(base + 1) = 0.5_dp*(chunk_density_grad(irow, 1) + &
                                             chunk_density_grad(irow, 2))
                  TARGET(base + 2) = 0.5_dp*(chunk_grad_grad(irow, 1, 1) + &
                                             chunk_grad_grad(irow, 1, 2))
                  TARGET(base + 3) = 0.5_dp*(chunk_grad_grad(irow, 2, 1) + &
                                             chunk_grad_grad(irow, 2, 2))
                  TARGET(base + 4) = 0.5_dp*(chunk_grad_grad(irow, 3, 1) + &
                                             chunk_grad_grad(irow, 3, 2))
                  TARGET(base + 5) = 0.5_dp*(chunk_kin_grad(irow, 1) + &
                                             chunk_kin_grad(irow, 2))
               END IF
            ELSE
               TARGET(base + 1:base + 2) = chunk_density_grad(irow, :)
               TARGET(base + 3) = chunk_grad_grad(irow, 1, 1)
               TARGET(base + 4) = chunk_grad_grad(irow, 2, 1)
               TARGET(base + 5) = chunk_grad_grad(irow, 3, 1)
               TARGET(base + 6) = chunk_grad_grad(irow, 1, 2)
               TARGET(base + 7) = chunk_grad_grad(irow, 2, 2)
               TARGET(base + 8) = chunk_grad_grad(irow, 3, 2)
               TARGET(base + 9:base + 10) = chunk_kin_grad(irow, :)
            END IF
         END DO
!$OMP END PARALLEL DO
      END IF

      IF (.NOT. PRESENT(field_response)) THEN
         CALL torch_tensor_release(density_grad_t)
         CALL torch_tensor_release(grad_grad_t)
         CALL torch_tensor_release(kin_grad_t)
      END IF

   END SUBROUTINE pack_atom_chunk_grads

! **************************************************************************************************
!> \brief Scatter routed atom-chunk gradients into local grid-row order.
!> \param features ...
!> \param recv_grad_buffer ...
!> \param collapse_spin_grads ...
!> \param density_grad ...
!> \param grad_grad ...
!> \param kin_grad ...
! **************************************************************************************************
   SUBROUTINE scatter_routed_atom_chunk_grads(features, recv_grad_buffer, collapse_spin_grads, &
                                              density_grad, grad_grad, kin_grad)
      TYPE(skala_gpw_feature_type), INTENT(IN)           :: features
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: recv_grad_buffer
      LOGICAL, INTENT(IN)                                :: collapse_spin_grads
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: density_grad
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(OUT)                                     :: grad_grad
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: kin_grad

      INTEGER                                            :: base, local_row, nflat_local, &
                                                            nroute_grad_per_point, nroute_points, &
                                                            point_pos, row_route_pos

      nflat_local = features%nflat_local
      nroute_points = SIZE(features%route_send_local_rows)
      nroute_grad_per_point = ngrad_per_point
      IF (collapse_spin_grads) nroute_grad_per_point = ncollapsed_grad_per_point
      CPASSERT(SIZE(recv_grad_buffer) >= nroute_grad_per_point*nroute_points)
      ALLOCATE (density_grad(nflat_local, 2), grad_grad(nflat_local, 3, 2), &
                kin_grad(nflat_local, 2))

!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP SHARED(collapse_spin_grads, density_grad, features, grad_grad, kin_grad, nflat_local, &
!$OMP        nroute_grad_per_point, nroute_points, recv_grad_buffer) &
!$OMP PRIVATE(base, local_row, point_pos, row_route_pos)
      DO local_row = 1, nflat_local
         density_grad(local_row, :) = 0.0_dp
         grad_grad(local_row, :, :) = 0.0_dp
         kin_grad(local_row, :) = 0.0_dp
         DO row_route_pos = features%route_row_offsets(local_row), &
            features%route_row_offsets(local_row + 1) - 1
            point_pos = features%route_row_positions(row_route_pos)
            CPASSERT(point_pos >= 1 .AND. point_pos <= nroute_points)
            base = nroute_grad_per_point*(point_pos - 1)
            IF (collapse_spin_grads) THEN
               density_grad(local_row, :) = density_grad(local_row, :) + &
                                            recv_grad_buffer(base + 1)
               grad_grad(local_row, 1, :) = grad_grad(local_row, 1, :) + &
                                            recv_grad_buffer(base + 2)
               grad_grad(local_row, 2, :) = grad_grad(local_row, 2, :) + &
                                            recv_grad_buffer(base + 3)
               grad_grad(local_row, 3, :) = grad_grad(local_row, 3, :) + &
                                            recv_grad_buffer(base + 4)
               kin_grad(local_row, :) = kin_grad(local_row, :) + recv_grad_buffer(base + 5)
            ELSE
               density_grad(local_row, :) = density_grad(local_row, :) + &
                                            recv_grad_buffer(base + 1:base + 2)
               grad_grad(local_row, 1, 1) = grad_grad(local_row, 1, 1) + &
                                            recv_grad_buffer(base + 3)
               grad_grad(local_row, 2, 1) = grad_grad(local_row, 2, 1) + &
                                            recv_grad_buffer(base + 4)
               grad_grad(local_row, 3, 1) = grad_grad(local_row, 3, 1) + &
                                            recv_grad_buffer(base + 5)
               grad_grad(local_row, 1, 2) = grad_grad(local_row, 1, 2) + &
                                            recv_grad_buffer(base + 6)
               grad_grad(local_row, 2, 2) = grad_grad(local_row, 2, 2) + &
                                            recv_grad_buffer(base + 7)
               grad_grad(local_row, 3, 2) = grad_grad(local_row, 3, 2) + &
                                            recv_grad_buffer(base + 8)
               kin_grad(local_row, :) = kin_grad(local_row, :) + &
                                        recv_grad_buffer(base + 9:base + 10)
            END IF
         END DO
      END DO
!$OMP END PARALLEL DO

   END SUBROUTINE scatter_routed_atom_chunk_grads

! **************************************************************************************************
!> \brief Return CPU views of autograd outputs for the SKALA dynamic feature tensors.
!> \param features ...
!> \param density_grad_t ...
!> \param grad_grad_t ...
!> \param kin_grad_t ...
!> \param density_grad ...
!> \param grad_grad ...
!> \param kin_grad ...
!> \param field_response Optional borrowed CPU responses. Returned tensor owners remain empty in this case.
! **************************************************************************************************
   SUBROUTINE get_feature_grad_views(features, density_grad_t, grad_grad_t, kin_grad_t, &
                                     density_grad, grad_grad, kin_grad, field_response)
      TYPE(skala_gpw_feature_type), INTENT(IN)           :: features
      TYPE(torch_tensor_type), INTENT(INOUT)             :: density_grad_t, grad_grad_t, kin_grad_t
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: density_grad
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: grad_grad
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: kin_grad
      TYPE(torch_tensor_type), INTENT(IN), OPTIONAL      :: field_response(3)

      NULLIFY (density_grad, grad_grad, kin_grad)
      IF (PRESENT(field_response)) THEN
         CALL torch_tensor_data_ptr(field_response(1), density_grad)
         CALL torch_tensor_data_ptr(field_response(2), grad_grad)
         CALL torch_tensor_data_ptr(field_response(3), kin_grad)
      ELSE
         CALL torch_tensor_grad_batch3(features%density_t, features%grad_t, features%kin_t, &
                                       density_grad_t, grad_grad_t, kin_grad_t)
         CALL torch_tensor_data_ptr(density_grad_t, density_grad)
         CALL torch_tensor_data_ptr(grad_grad_t, grad_grad)
         CALL torch_tensor_data_ptr(kin_grad_t, kin_grad)
      END IF

   END SUBROUTINE get_feature_grad_views

! **************************************************************************************************
!> \brief Fetch atom-chunk gradients and route them back to their local grid owners.
!> \param features ...
!> \param group ...
!> \param density_grad ...
!> \param grad_grad ...
!> \param kin_grad ...
!> \param field_response Optional borrowed weighted response tensors in the same feature layout.
! **************************************************************************************************
   SUBROUTINE fetch_and_gather_atom_chunk_grads(features, group, density_grad, grad_grad, &
                                                kin_grad, field_response)
      TYPE(skala_gpw_feature_type), INTENT(IN)           :: features
      TYPE(torch_tensor_type), INTENT(IN), OPTIONAL      :: field_response(3)

      CLASS(mp_comm_type), INTENT(IN)                    :: group
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: density_grad, kin_grad
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(OUT)                                     :: grad_grad

      INTEGER                                            :: base, feature_pos, i, j, k, local_row, &
                                                            nflat_local, nroute_grad_per_point, &
                                                            nroute_recv_points, nroute_send_points, &
                                                            phase_handle, row
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: route_grad_return_recv_counts, &
                                                            route_grad_return_recv_displs, &
                                                            route_grad_return_send_counts, &
                                                            route_grad_return_send_displs
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: chunk_grad_buffer, global_grad_buffer, &
                                                            recv_grad_buffer, send_grad_buffer

      CPASSERT(features%uses_atom_chunks)

      nflat_local = features%nflat_local
      IF (features%uses_atom_chunk_routing) THEN
         nroute_recv_points = SUM(features%route_point_recv_counts)
         nroute_send_points = SIZE(features%route_send_local_rows)
         CPASSERT(SUM(features%route_point_send_counts) == nroute_send_points)
         CPASSERT(SIZE(features%route_chunk_offsets) == nroute_recv_points + 1)

         nroute_grad_per_point = ngrad_per_point
         IF (features%uses_collapsed_rks_dynamic) THEN
            nroute_grad_per_point = ncollapsed_grad_per_point
         END IF
         ALLOCATE (send_grad_buffer(MAX(1, nroute_grad_per_point*nroute_recv_points)), &
                   recv_grad_buffer(MAX(1, nroute_grad_per_point*nroute_send_points)), &
                   route_grad_return_send_counts(SIZE(features%route_point_recv_counts)), &
                   route_grad_return_send_displs(SIZE(features%route_point_recv_displs)), &
                   route_grad_return_recv_counts(SIZE(features%route_point_send_counts)), &
                   route_grad_return_recv_displs(SIZE(features%route_point_send_displs)))
         route_grad_return_send_counts(:) = &
            nroute_grad_per_point*features%route_point_recv_counts
         route_grad_return_send_displs(:) = &
            nroute_grad_per_point*features%route_point_recv_displs
         route_grad_return_recv_counts(:) = &
            nroute_grad_per_point*features%route_point_send_counts
         route_grad_return_recv_displs(:) = &
            nroute_grad_per_point*features%route_point_send_displs

         IF (features%chunk_feature_count > 0) THEN
            CALL timeset("skala_gpw_grad_torch_pack", phase_handle)
            CALL pack_atom_chunk_grads(features, send_grad_buffer, .TRUE., &
                                       features%uses_collapsed_rks_dynamic, field_response)
            CALL timestop(phase_handle)
         END IF

         CALL timeset("skala_gpw_grad_route_comm", phase_handle)
         CALL group%alltoall(send_grad_buffer, route_grad_return_send_counts, &
                             route_grad_return_send_displs, recv_grad_buffer, &
                             route_grad_return_recv_counts, route_grad_return_recv_displs)
         CALL timestop(phase_handle)

         CALL timeset("skala_gpw_grad_route_scatter", phase_handle)
         CALL scatter_routed_atom_chunk_grads(features, recv_grad_buffer, &
                                              features%uses_collapsed_rks_dynamic, &
                                              density_grad, grad_grad, kin_grad)
         CALL timestop(phase_handle)

         DEALLOCATE (recv_grad_buffer, route_grad_return_recv_counts, &
                     route_grad_return_recv_displs, route_grad_return_send_counts, &
                     route_grad_return_send_displs, send_grad_buffer)
      ELSE
         ALLOCATE (chunk_grad_buffer(MAX(1, ngrad_per_point*features%chunk_feature_count)), &
                   global_grad_buffer(ngrad_per_point*features%nflat))
         IF (features%chunk_feature_count > 0) THEN
            CALL timeset("skala_gpw_grad_torch_pack", phase_handle)
            CALL pack_atom_chunk_grads(features, chunk_grad_buffer, .FALSE., field_response=field_response)
            CALL timestop(phase_handle)
         END IF

         CALL timeset("skala_gpw_grad_allgatherv", phase_handle)
         CALL group%allgatherv(chunk_grad_buffer, global_grad_buffer, &
                               features%chunk_grad_counts, features%chunk_grad_displs)
         CALL timestop(phase_handle)

         CALL timeset("skala_gpw_grad_scatter", phase_handle)
         ALLOCATE (density_grad(nflat_local, 2), grad_grad(nflat_local, 3, 2), &
                   kin_grad(nflat_local, 2))
         density_grad = 0.0_dp
         grad_grad = 0.0_dp
         kin_grad = 0.0_dp
         local_row = 0
         DO k = LBOUND(features%feature_index, 3), UBOUND(features%feature_index, 3)
            DO j = LBOUND(features%feature_index, 2), UBOUND(features%feature_index, 2)
               DO i = LBOUND(features%feature_index, 1), UBOUND(features%feature_index, 1)
                  local_row = local_row + 1
                  DO feature_pos = features%local_feature_offsets(local_row), &
                     features%local_feature_offsets(local_row + 1) - 1
                     row = features%local_feature_rows(feature_pos)
                     CPASSERT(row >= 1 .AND. row <= features%nflat)
                     base = ngrad_per_point*(row - 1)
                     density_grad(local_row, :) = density_grad(local_row, :) + &
                                                  global_grad_buffer(base + 1:base + 2)
                     grad_grad(local_row, 1, 1) = grad_grad(local_row, 1, 1) + &
                                                  global_grad_buffer(base + 3)
                     grad_grad(local_row, 2, 1) = grad_grad(local_row, 2, 1) + &
                                                  global_grad_buffer(base + 4)
                     grad_grad(local_row, 3, 1) = grad_grad(local_row, 3, 1) + &
                                                  global_grad_buffer(base + 5)
                     grad_grad(local_row, 1, 2) = grad_grad(local_row, 1, 2) + &
                                                  global_grad_buffer(base + 6)
                     grad_grad(local_row, 2, 2) = grad_grad(local_row, 2, 2) + &
                                                  global_grad_buffer(base + 7)
                     grad_grad(local_row, 3, 2) = grad_grad(local_row, 3, 2) + &
                                                  global_grad_buffer(base + 8)
                     kin_grad(local_row, :) = kin_grad(local_row, :) + &
                                              global_grad_buffer(base + 9:base + 10)
                  END DO
               END DO
            END DO
         END DO
         CALL timestop(phase_handle)
         DEALLOCATE (chunk_grad_buffer, global_grad_buffer)

      END IF

   END SUBROUTINE fetch_and_gather_atom_chunk_grads

! **************************************************************************************************
!> \brief Build the native SKALA XC virial from feature gradients.
!> \param virial_xc ...
!> \param rho_set ...
!> \param rho_r ...
!> \param grad_grad ...
! **************************************************************************************************
   SUBROUTINE build_virial_from_feature_grads(virial_xc, rho_set, rho_r, grad_grad)
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(INOUT)      :: virial_xc
      TYPE(xc_rho_set_type), INTENT(IN)                  :: rho_set
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: rho_r
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: grad_grad

      INTEGER                                            :: i, idir, ipt, ispin, j, jdir, k, nspins
      INTEGER, DIMENSION(2, 3)                           :: bo
      REAL(KIND=dp)                                      :: grad_i, tmp
      TYPE(cp_3d_r_cp_type), DIMENSION(3)                :: drho, drhoa, drhob

      nspins = SIZE(rho_r)
      bo = rho_r(1)%pw_grid%bounds_local
      ipt = 0

      IF (nspins == 1) THEN
         CALL xc_rho_set_get(rho_set, drho=drho)
         DO k = bo(1, 3), bo(2, 3)
            DO j = bo(1, 2), bo(2, 2)
               DO i = bo(1, 1), bo(2, 1)
                  ipt = ipt + 1
                  DO idir = 1, 3
                     grad_i = 0.5_dp*(grad_grad(ipt, idir, 1) + grad_grad(ipt, idir, 2))
                     DO jdir = 1, 3
                        tmp = -grad_i*drho(jdir)%array(i, j, k)
                        virial_xc(jdir, idir) = virial_xc(jdir, idir) + tmp
                     END DO
                  END DO
               END DO
            END DO
         END DO
      ELSE
         CALL xc_rho_set_get(rho_set, drhoa=drhoa, drhob=drhob)
         DO k = bo(1, 3), bo(2, 3)
            DO j = bo(1, 2), bo(2, 2)
               DO i = bo(1, 1), bo(2, 1)
                  ipt = ipt + 1
                  DO idir = 1, 3
                     DO jdir = 1, 3
                        tmp = 0.0_dp
                        DO ispin = 1, 2
                           IF (ispin == 1) THEN
                              tmp = tmp - grad_grad(ipt, idir, ispin)*drhoa(jdir)%array(i, j, k)
                           ELSE
                              tmp = tmp - grad_grad(ipt, idir, ispin)*drhob(jdir)%array(i, j, k)
                           END IF
                        END DO
                        virial_xc(jdir, idir) = virial_xc(jdir, idir) + tmp
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END IF

   END SUBROUTINE build_virial_from_feature_grads

! **************************************************************************************************
!> \brief Print a native SKALA XC virial contribution for diagnostics.
!> \param label ...
!> \param delta ...
!> \param root_rank ...
! **************************************************************************************************
   SUBROUTINE print_virial_delta(label, delta, root_rank)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: delta
      LOGICAL, INTENT(IN)                                :: root_rank

      INTEGER                                            :: i, iw

      IF (.NOT. root_rank) RETURN
      iw = cp_logger_get_default_io_unit()
      IF (iw <= 0) RETURN
      WRITE (iw, "(T2,A,1X,A)") "SKALA_GPW| XC virial contribution", TRIM(label)
      DO i = 1, 3
         WRITE (iw, "(T2,A,1X,3ES20.10)") "SKALA_GPW|", delta(i, 1:3)
      END DO

   END SUBROUTINE print_virial_delta

! **************************************************************************************************
!> \brief Add explicit SKALA coordinate-feature contributions to the XC virial.
!> \param virial_xc ...
!> \param features ...
!> \param atom_coord_grad_t ...
!> \param grid_coord_grad_t ...
!> \param root_rank ...
!> \param print_components ...
! **************************************************************************************************
   SUBROUTINE build_static_coordinate_virial(virial_xc, features, atom_coord_grad_t, &
                                             grid_coord_grad_t, root_rank, print_components)
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(INOUT)      :: virial_xc
      TYPE(skala_gpw_feature_type), INTENT(IN)           :: features
      TYPE(torch_tensor_type), INTENT(INOUT)             :: atom_coord_grad_t, grid_coord_grad_t
      LOGICAL, INTENT(IN)                                :: root_rank
      LOGICAL, INTENT(IN), OPTIONAL                      :: print_components

      INTEGER                                            :: feature_pos, i, iatom, idir, iw, j, &
                                                            jdir, k, local_row, row
      LOGICAL                                            :: my_print_components
      REAL(KIND=dp)                                      :: tmp
      REAL(KIND=dp), DIMENSION(3)                        :: atom_grad_sum, grid_grad_sum
      REAL(KIND=dp), DIMENSION(3, 3)                     :: atom_virial, grid_virial
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: atom_coord_grad, grid_coord_grad

      my_print_components = .FALSE.
      IF (PRESENT(print_components)) my_print_components = print_components

      NULLIFY (atom_coord_grad, grid_coord_grad)
      CALL torch_tensor_grad(features%grid_coords_t, grid_coord_grad_t)
      CALL torch_tensor_data_ptr(grid_coord_grad_t, grid_coord_grad)
      CALL torch_tensor_data_ptr(atom_coord_grad_t, atom_coord_grad)

      grid_virial = 0.0_dp
      atom_virial = 0.0_dp
      grid_grad_sum = 0.0_dp
      atom_grad_sum = 0.0_dp
      local_row = 0
      DO k = LBOUND(features%feature_index, 3), UBOUND(features%feature_index, 3)
         DO j = LBOUND(features%feature_index, 2), UBOUND(features%feature_index, 2)
            DO i = LBOUND(features%feature_index, 1), UBOUND(features%feature_index, 1)
               local_row = local_row + 1
               DO feature_pos = features%local_feature_offsets(local_row), &
                  features%local_feature_offsets(local_row + 1) - 1
                  row = features%local_feature_rows(feature_pos)
                  grid_grad_sum(:) = grid_grad_sum(:) + grid_coord_grad(:, row)
                  DO idir = 1, 3
                     DO jdir = 1, 3
                        tmp = grid_coord_grad(idir, row)*features%grid_coords(jdir, row)
                        grid_virial(idir, jdir) = grid_virial(idir, jdir) + tmp
                        virial_xc(idir, jdir) = virial_xc(idir, jdir) + tmp
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
      CPASSERT(local_row == features%nflat_local)

      IF (root_rank) THEN
         DO iatom = 1, SIZE(features%coarse_0_atomic_coords, 2)
            atom_grad_sum(:) = atom_grad_sum(:) + atom_coord_grad(:, iatom)
            DO idir = 1, 3
               DO jdir = 1, 3
                  tmp = atom_coord_grad(idir, iatom)*features%coarse_0_atomic_coords(jdir, iatom)
                  atom_virial(idir, jdir) = atom_virial(idir, jdir) + tmp
                  virial_xc(idir, jdir) = virial_xc(idir, jdir) + tmp
               END DO
            END DO
         END DO
      END IF

      IF (my_print_components .AND. root_rank) THEN
         iw = cp_logger_get_default_io_unit()
         IF (iw > 0) THEN
            CALL print_virial_delta("static-grid", grid_virial, .TRUE.)
            CALL print_virial_delta("static-atom", atom_virial, .TRUE.)
            WRITE (iw, "(T2,A,1X,3ES20.10)") "SKALA_GPW| XC coordinate gradient grid sum", &
               grid_grad_sum
            WRITE (iw, "(T2,A,1X,3ES20.10)") "SKALA_GPW| XC coordinate gradient atom sum", &
               atom_grad_sum
            WRITE (iw, "(T2,A,1X,3ES20.10)") "SKALA_GPW| XC coordinate gradient total sum", &
               grid_grad_sum + atom_grad_sum
         END IF
      END IF

      CALL torch_tensor_release(grid_coord_grad_t)

   END SUBROUTINE build_static_coordinate_virial

! **************************************************************************************************
!> \brief Add residual SKALA weight-feature contributions to the XC virial.
!> \param virial_xc ...
!> \param features ...
!> \param exc ...
!> \param grid_weight_grad_t ...
!> \param atomic_grid_weight_grad_t ...
!> \param root_rank ...
!> \param print_components ...
! **************************************************************************************************
   SUBROUTINE build_weight_virial(virial_xc, features, exc, grid_weight_grad_t, &
                                  atomic_grid_weight_grad_t, root_rank, print_components)
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(INOUT)      :: virial_xc
      TYPE(skala_gpw_feature_type), INTENT(IN)           :: features
      REAL(KIND=dp), INTENT(IN)                          :: exc
      TYPE(torch_tensor_type), INTENT(INOUT)             :: grid_weight_grad_t, &
                                                            atomic_grid_weight_grad_t
      LOGICAL, INTENT(IN)                                :: root_rank
      LOGICAL, INTENT(IN), OPTIONAL                      :: print_components

      INTEGER                                            :: feature_pos, i, idir, iw, j, k, &
                                                            local_row, row
      LOGICAL                                            :: my_print_components
      REAL(KIND=dp)                                      :: atomic_tmp, exc_tmp, grid_tmp, tmp
      REAL(KIND=dp), DIMENSION(:), POINTER               :: atomic_grid_weight_grad, grid_weight_grad

      my_print_components = .FALSE.
      IF (PRESENT(print_components)) my_print_components = print_components

      NULLIFY (atomic_grid_weight_grad, grid_weight_grad)
      CALL torch_tensor_grad(features%grid_weights_t, grid_weight_grad_t)
      CALL torch_tensor_grad(features%atomic_grid_weights_t, atomic_grid_weight_grad_t)
      CALL torch_tensor_data_ptr(grid_weight_grad_t, grid_weight_grad)
      CALL torch_tensor_data_ptr(atomic_grid_weight_grad_t, atomic_grid_weight_grad)

      grid_tmp = 0.0_dp
      atomic_tmp = 0.0_dp
      local_row = 0
      DO k = LBOUND(features%feature_index, 3), UBOUND(features%feature_index, 3)
         DO j = LBOUND(features%feature_index, 2), UBOUND(features%feature_index, 2)
            DO i = LBOUND(features%feature_index, 1), UBOUND(features%feature_index, 1)
               local_row = local_row + 1
               DO feature_pos = features%local_feature_offsets(local_row), &
                  features%local_feature_offsets(local_row + 1) - 1
                  row = features%local_feature_rows(feature_pos)
                  grid_tmp = grid_tmp + grid_weight_grad(row)*features%grid_weights(row)
                  atomic_tmp = atomic_tmp + &
                               atomic_grid_weight_grad(row)*features%atomic_grid_weights(row)
               END DO
            END DO
         END DO
      END DO
      CPASSERT(local_row == features%nflat_local)
      exc_tmp = 0.0_dp
      IF (root_rank) exc_tmp = -exc
      tmp = grid_tmp + atomic_tmp + exc_tmp

      IF (my_print_components .AND. root_rank) THEN
         iw = cp_logger_get_default_io_unit()
         IF (iw > 0) THEN
            WRITE (iw, "(T2,A,1X,ES20.10)") "SKALA_GPW| XC virial weight grid", grid_tmp
            WRITE (iw, "(T2,A,1X,ES20.10)") "SKALA_GPW| XC virial weight atomic", atomic_tmp
            WRITE (iw, "(T2,A,1X,ES20.10)") "SKALA_GPW| XC virial weight final", exc_tmp
            WRITE (iw, "(T2,A,1X,ES20.10)") "SKALA_GPW| XC virial weight residual", tmp
         END IF
      END IF

      DO idir = 1, 3
         virial_xc(idir, idir) = virial_xc(idir, idir) + tmp
      END DO

      CALL torch_tensor_release(grid_weight_grad_t)
      CALL torch_tensor_release(atomic_grid_weight_grad_t)

   END SUBROUTINE build_weight_virial

! **************************************************************************************************
!> \brief Fill CP2K VXC real-space arrays from Torch feature gradients.
!> \param vxc_rho ...
!> \param vxc_tau ...
!> \param rho_r ...
!> \param pw_pool ...
!> \param density_grad ...
!> \param grad_grad ...
!> \param kin_grad ...
!> \param xc_deriv_method_id ...
!> \param global_grid_layout ...
! **************************************************************************************************
   SUBROUTINE build_vxc_from_feature_grads(vxc_rho, vxc_tau, rho_r, pw_pool, &
                                           density_grad, grad_grad, kin_grad, &
                                           xc_deriv_method_id, global_grid_layout)
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: vxc_rho, vxc_tau, rho_r
      TYPE(pw_pool_type), POINTER                        :: pw_pool
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: density_grad
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: grad_grad
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: kin_grad
      INTEGER, INTENT(IN)                                :: xc_deriv_method_id
      LOGICAL, INTENT(IN), OPTIONAL                      :: global_grid_layout

      INTEGER                                            :: i, ipt, ispin, j, k, nspins
      INTEGER, DIMENSION(2, 3)                           :: bo
      LOGICAL                                            :: my_global_grid_layout, valid_spin_shape
      REAL(KIND=dp)                                      :: dvol_inv
      TYPE(pw_c1d_gs_type)                               :: tmp_g, vxc_g
      TYPE(pw_r3d_rs_type), DIMENSION(3)                 :: grad_pw

      nspins = SIZE(rho_r)
      bo = rho_r(1)%pw_grid%bounds_local
      dvol_inv = 1.0_dp/rho_r(1)%pw_grid%dvol
      my_global_grid_layout = .FALSE.
      IF (PRESENT(global_grid_layout)) my_global_grid_layout = global_grid_layout
      IF (my_global_grid_layout) THEN
         CPASSERT(SIZE(density_grad, 1) == PRODUCT(rho_r(1)%pw_grid%npts))
      END IF
      CPASSERT(SIZE(kin_grad, 2) == SIZE(density_grad, 2))
      CPASSERT(SIZE(grad_grad, 3) == SIZE(density_grad, 2))
      CPASSERT(SIZE(grad_grad, 2) == 3)
      valid_spin_shape = SIZE(density_grad, 2) == nspins .OR. &
                         (nspins == 1 .AND. SIZE(density_grad, 2) == 2)
      CPASSERT(valid_spin_shape)

      ALLOCATE (vxc_rho(nspins), vxc_tau(nspins))
      DO ispin = 1, nspins
         CALL pw_pool%create_pw(vxc_rho(ispin))
         CALL pw_pool%create_pw(vxc_tau(ispin))
         CALL pw_zero(vxc_rho(ispin))
         CALL pw_zero(vxc_tau(ispin))
      END DO

      IF (xc_requires_tmp_g(xc_deriv_method_id) .OR. rho_r(1)%pw_grid%spherical) THEN
         CALL pw_pool%create_pw(vxc_g)
         IF (.NOT. rho_r(1)%pw_grid%spherical) CALL pw_pool%create_pw(tmp_g)
      END IF

      DO ispin = 1, nspins
         DO i = 1, 3
            CALL pw_pool%create_pw(grad_pw(i))
            CALL pw_zero(grad_pw(i))
         END DO

         ipt = 0
         DO k = bo(1, 3), bo(2, 3)
            DO j = bo(1, 2), bo(2, 2)
               DO i = bo(1, 1), bo(2, 1)
                  IF (my_global_grid_layout) THEN
                     ipt = 1 + i - rho_r(1)%pw_grid%bounds(1, 1) + &
                           rho_r(1)%pw_grid%npts(1)*( &
                           j - rho_r(1)%pw_grid%bounds(1, 2) + &
                           rho_r(1)%pw_grid%npts(2)*( &
                           k - rho_r(1)%pw_grid%bounds(1, 3)))
                  ELSE
                     ipt = ipt + 1
                  END IF
                  IF (nspins == 1) THEN
                     IF (SIZE(density_grad, 2) == 1) THEN
                        vxc_rho(1)%array(i, j, k) = 0.5_dp*dvol_inv*density_grad(ipt, 1)
                        vxc_tau(1)%array(i, j, k) = 0.5_dp*dvol_inv*kin_grad(ipt, 1)
                        grad_pw(1)%array(i, j, k) = 0.5_dp*dvol_inv*grad_grad(ipt, 1, 1)
                        grad_pw(2)%array(i, j, k) = 0.5_dp*dvol_inv*grad_grad(ipt, 2, 1)
                        grad_pw(3)%array(i, j, k) = 0.5_dp*dvol_inv*grad_grad(ipt, 3, 1)
                     ELSE
                        vxc_rho(1)%array(i, j, k) = 0.5_dp*dvol_inv* &
                                                    (density_grad(ipt, 1) + density_grad(ipt, 2))
                        vxc_tau(1)%array(i, j, k) = 0.5_dp*dvol_inv* &
                                                    (kin_grad(ipt, 1) + kin_grad(ipt, 2))
                        grad_pw(1)%array(i, j, k) = 0.5_dp*dvol_inv* &
                                                    (grad_grad(ipt, 1, 1) + grad_grad(ipt, 1, 2))
                        grad_pw(2)%array(i, j, k) = 0.5_dp*dvol_inv* &
                                                    (grad_grad(ipt, 2, 1) + grad_grad(ipt, 2, 2))
                        grad_pw(3)%array(i, j, k) = 0.5_dp*dvol_inv* &
                                                    (grad_grad(ipt, 3, 1) + grad_grad(ipt, 3, 2))
                     END IF
                  ELSE
                     vxc_rho(ispin)%array(i, j, k) = dvol_inv*density_grad(ipt, ispin)
                     vxc_tau(ispin)%array(i, j, k) = dvol_inv*kin_grad(ipt, ispin)
                     grad_pw(1)%array(i, j, k) = dvol_inv*grad_grad(ipt, 1, ispin)
                     grad_pw(2)%array(i, j, k) = dvol_inv*grad_grad(ipt, 2, ispin)
                     grad_pw(3)%array(i, j, k) = dvol_inv*grad_grad(ipt, 3, ispin)
                  END IF
               END DO
            END DO
         END DO

         DO i = 1, 3
            CALL pw_scale(grad_pw(i), -1.0_dp)
         END DO
         CALL xc_pw_divergence(xc_deriv_method_id, grad_pw, tmp_g, vxc_g, vxc_rho(ispin))

         DO i = 1, 3
            CALL pw_pool%give_back_pw(grad_pw(i))
         END DO
      END DO

      IF (ASSOCIATED(vxc_g%pw_grid)) CALL pw_pool%give_back_pw(vxc_g)
      IF (ASSOCIATED(tmp_g%pw_grid)) CALL pw_pool%give_back_pw(tmp_g)

   END SUBROUTINE build_vxc_from_feature_grads

! **************************************************************************************************
!> \brief Print optional diagnostics for the CP2K-native SKALA GPW feature block.
!> \param features ...
!> \param print_active ...
! **************************************************************************************************
   SUBROUTINE print_native_grid_diagnostics(features, print_active)
      TYPE(skala_gpw_feature_type), INTENT(IN)           :: features
      LOGICAL, INTENT(IN)                                :: print_active

      INTEGER                                            :: atom_rows_max, atom_rows_min, &
                                                            chunk_rows_max, chunk_rows_min, iw
      REAL(KIND=dp)                                      :: chunk_imbalance

      IF (.NOT. print_active) RETURN

      iw = cp_logger_get_default_io_unit()
      IF (iw <= 0) RETURN
      WRITE (UNIT=iw, FMT="(/,T2,A,1X,ES19.11)") &
         "SKALA_GPW| Native grid feature electrons", features%electron_count
      WRITE (UNIT=iw, FMT="(T2,A,1X,ES19.11)") &
         "SKALA_GPW| Native grid feature spin moment", features%spin_moment
      WRITE (UNIT=iw, FMT="(T2,A,1X,ES19.11)") &
         "SKALA_GPW| Native grid feature kinetic integral", features%kinetic_integral
      WRITE (UNIT=iw, FMT="(T2,A,1X,ES19.11)") &
         "SKALA_GPW| Native grid feature weight sum", features%grid_weight_sum
      IF (ALLOCATED(features%atomic_grid_sizes)) THEN
         atom_rows_min = INT(MINVAL(features%atomic_grid_sizes))
         atom_rows_max = INT(MAXVAL(features%atomic_grid_sizes))
         WRITE (UNIT=iw, FMT="(T2,A,1X,I0,1X,A,1X,I0,1X,A,1X,I0)") &
            "SKALA_GPW| Native grid atom row range", atom_rows_min, "to", &
            atom_rows_max, "sum", INT(SUM(features%atomic_grid_sizes))
      END IF
      IF (features%uses_atom_chunks) THEN
         WRITE (UNIT=iw, FMT="(T2,A,1X,I0,1X,A,1X,I0)") &
            "SKALA_GPW| Native grid atom chunk rows", features%chunk_feature_count, &
            "of", features%nflat
         IF (ALLOCATED(features%chunk_grad_counts)) THEN
            chunk_rows_min = MINVAL(features%chunk_grad_counts)/ngrad_per_point
            chunk_rows_max = MAXVAL(features%chunk_grad_counts)/ngrad_per_point
            chunk_imbalance = REAL(chunk_rows_max, KIND=dp)/REAL(MAX(1, chunk_rows_min), KIND=dp)
            WRITE (UNIT=iw, FMT="(T2,A,1X,I0,1X,A,1X,I0,1X,A,1X,ES12.5)") &
               "SKALA_GPW| Native grid atom chunk row range", chunk_rows_min, &
               "to", chunk_rows_max, "imbalance", chunk_imbalance
         END IF
      END IF

   END SUBROUTINE print_native_grid_diagnostics

! **************************************************************************************************
!> \brief Configure CUDA device selection for the native SKALA GPW Torch path.
!> \param use_cuda ...
!> \param requested_device ...
!> \param group ...
!> \return selected CUDA device, or -1 for CPU fallback/no visible CUDA device
! **************************************************************************************************
   FUNCTION configure_native_grid_cuda(use_cuda, requested_device, group) RESULT(selected_device)
      LOGICAL, INTENT(IN)                                :: use_cuda
      INTEGER, INTENT(IN)                                :: requested_device

      CLASS(mp_comm_type), INTENT(IN)                    :: group

      INTEGER                                            :: cuda_device_count, iw, pe, selected_device
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: selected_devices

      selected_device = -1

      IF (.NOT. use_cuda) RETURN

      IF (.NOT. torch_cuda_is_available()) THEN
         cuda_device_count = 0
      ELSE
         cuda_device_count = torch_cuda_device_count()
      END IF
      IF (cuda_device_count > 0) THEN
         IF (requested_device < 0) THEN
            selected_device = MOD(group%mepos, cuda_device_count)
         ELSE
            selected_device = requested_device
         END IF
      END IF
      IF (selected_device >= cuda_device_count) THEN
         CALL cp_abort(__LOCATION__, &
                       "GAUXC%NATIVE_GRID_CUDA_DEVICE selects a CUDA device outside the visible "// &
                       "Torch CUDA device range.")
      END IF
      IF (selected_device >= 0) CALL offload_set_chosen_device(selected_device)

      ALLOCATE (selected_devices(group%num_pe))
      CALL group%allgather(selected_device, selected_devices)

      IF (group%mepos /= 0) THEN
         DEALLOCATE (selected_devices)
         RETURN
      END IF
      IF (selected_device == logged_cuda_device .AND. &
          cuda_device_count == logged_cuda_device_count .AND. &
          group%num_pe == logged_cuda_nproc .AND. &
          requested_device == logged_cuda_request) THEN
         DEALLOCATE (selected_devices)
         RETURN
      END IF

      iw = cp_logger_get_default_io_unit()
      IF (iw <= 0) THEN
         DEALLOCATE (selected_devices)
         RETURN
      END IF
      IF (selected_device >= 0) THEN
         WRITE (UNIT=iw, FMT="(/,T2,A,1X,I0,1X,A,1X,I0,1X,A,1X,I0)") &
            "SKALA_GPW| Native grid Torch CUDA device", selected_device, &
            "of", cuda_device_count, "requested", requested_device
      ELSE
         WRITE (UNIT=iw, FMT="(/,T2,A)") &
            "SKALA_GPW| Native grid Torch CUDA requested, but no Torch CUDA device is visible"
      END IF
      WRITE (UNIT=iw, FMT="(T2,A)", ADVANCE="NO") &
         "SKALA_GPW| Native grid Torch CUDA rank devices"
      DO pe = 1, group%num_pe
         WRITE (UNIT=iw, FMT="(1X,I0,A,I0)", ADVANCE="NO") pe - 1, ":", selected_devices(pe)
      END DO
      WRITE (UNIT=iw, FMT=*)

      logged_cuda_device = selected_device
      logged_cuda_device_count = cuda_device_count
      logged_cuda_nproc = group%num_pe
      logged_cuda_request = requested_device
      DEALLOCATE (selected_devices)

   END FUNCTION configure_native_grid_cuda

! **************************************************************************************************
!> \brief Share the current model's spin-density clamp bounds, including ranks with no atom rows.
!> \param xc_section Current functional settings. A native evaluation must already have loaded the model.
!> \param group Communicator of the native density samples.
!> \return Physical density bounds. Empty for functionals using another evaluator.
! **************************************************************************************************
   FUNCTION native_skala_spin_floors(xc_section, group) RESULT(floors)
      TYPE(section_vals_type), INTENT(IN), POINTER       :: xc_section
      TYPE(mp_comm_type), INTENT(IN)                     :: group
      REAL(dp), ALLOCATABLE                              :: floors(:)

      CHARACTER(len=default_path_length)                 :: model_path
      INTEGER                                            :: count, owner
      LOGICAL                                            :: fp64
      REAL(dp), ALLOCATABLE                              :: local_floors(:)
      TYPE(section_vals_type), POINTER                   :: gauxc_section

      IF (.NOT. xc_section_uses_native_skala_evaluator(xc_section)) THEN
         ALLOCATE (floors(0))
         RETURN
      END IF
      CALL get_skala_model_path(xc_section, model_path)
      gauxc_section => get_gauxc_section(xc_section)
      CALL section_vals_val_get(gauxc_section, "NATIVE_GRID_MODEL_FP64", l_val=fp64)
      owner = HUGE(owner)
      IF (cached_model_loaded) THEN
         IF (TRIM(cached_model_path) == TRIM(model_path) .AND. (cached_model_fp64 .EQV. fp64)) THEN
            owner = group%mepos
         END IF
      END IF
      CALL group%min(owner)
      CPASSERT(owner < group%num_pe)
      count = 0
      IF (group%mepos == owner) THEN
         ALLOCATE (local_floors, SOURCE=skala_torch_model_density_floors(cached_model))
         count = SIZE(local_floors)
      END IF
      CALL group%max(count)
      ALLOCATE (floors(count), SOURCE=0.0_dp)
      IF (group%mepos == owner) floors(:) = local_floors
      IF (count > 0) CALL group%sum(floors)
   END FUNCTION native_skala_spin_floors

! **************************************************************************************************
!> \brief Load and cache the TorchScript SKALA model.
!> \param model_path ...
!> \param cuda_device ...
!> \param xc_section XC settings, including the native model precision.
! **************************************************************************************************
   SUBROUTINE ensure_model_loaded(model_path, cuda_device, xc_section)
      CHARACTER(len=*), INTENT(IN)                       :: model_path
      INTEGER, INTENT(IN)                                :: cuda_device
      TYPE(section_vals_type), POINTER                   :: xc_section

      LOGICAL                                            :: fp64
      TYPE(section_vals_type), POINTER                   :: gauxc_section

      gauxc_section => get_gauxc_section(xc_section)
      CALL section_vals_val_get(gauxc_section, "NATIVE_GRID_MODEL_FP64", l_val=fp64)

      IF (cached_model_loaded) THEN
         IF (TRIM(cached_model_path) == TRIM(model_path) .AND. &
             cached_model_cuda_device == cuda_device .AND. &
             (cached_model_fp64 .EQV. fp64)) RETURN
         CALL skala_torch_model_release(cached_model)
         cached_model_loaded = .FALSE.
      END IF

      CALL skala_torch_model_load(cached_model, TRIM(model_path), fp64=fp64)
      cached_model_path = model_path
      cached_model_cuda_device = cuda_device
      cached_model_fp64 = fp64
      cached_model_loaded = .TRUE.

   END SUBROUTINE ensure_model_loaded

! **************************************************************************************************
!> \brief Resolve the SKALA TorchScript model path from the GAUXC subsection.
!> \param xc_section ...
!> \param model_path ...
! **************************************************************************************************
   SUBROUTINE get_skala_model_path(xc_section, model_path)
      TYPE(section_vals_type), INTENT(IN), POINTER       :: xc_section
      CHARACTER(len=default_path_length), INTENT(OUT)    :: model_path

      CHARACTER(len=default_path_length)                 :: model_key
      INTEGER                                            :: env_status
      LOGICAL                                            :: native_grid_use_cuda
      TYPE(section_vals_type), POINTER                   :: gauxc_section

      gauxc_section => get_gauxc_section(xc_section)
      IF (.NOT. ASSOCIATED(gauxc_section)) THEN
         CPABORT("Native SKALA GPW requires an XC_FUNCTIONAL%GAUXC section")
      END IF

      CALL section_vals_val_get(gauxc_section, "MODEL", c_val=model_path)
      model_key = ADJUSTL(model_path)
      CALL uppercase(model_key)
      IF (TRIM(model_key) == "NONE" .OR. TRIM(model_key) == "") THEN
         CPABORT("Native SKALA GPW requires GAUXC%MODEL SKALA or a TorchScript model path")
      ELSE IF (TRIM(model_key) == "SKALA") THEN
         CALL section_vals_val_get(gauxc_section, "NATIVE_GRID_USE_CUDA", l_val=native_grid_use_cuda)
         IF (native_grid_use_cuda) THEN
            CALL GET_ENVIRONMENT_VARIABLE("GAUXC_SKALA_CUDA_MODEL", model_path, STATUS=env_status)
            IF (env_status == 0 .AND. LEN_TRIM(model_path) > 0) RETURN
         END IF
         CALL GET_ENVIRONMENT_VARIABLE("SKALA_MODEL", model_path, STATUS=env_status)
         IF (env_status /= 0 .OR. LEN_TRIM(model_path) == 0) THEN
            IF (native_grid_use_cuda) THEN
               CALL cp_abort(__LOCATION__, &
                             "MODEL SKALA CUDA path requires GAUXC_SKALA_CUDA_MODEL or SKALA_MODEL")
            ELSE
               CALL cp_abort(__LOCATION__, &
                             "MODEL SKALA requires the SKALA_MODEL environment variable")
            END IF
         END IF
      END IF

   END SUBROUTINE get_skala_model_path

! **************************************************************************************************
!> \brief Return the first GAUXC functional subsection, if present.
!> \param xc_section ...
!> \return ...
! **************************************************************************************************
   FUNCTION get_gauxc_section(xc_section) RESULT(gauxc_section)
      TYPE(section_vals_type), INTENT(IN), POINTER       :: xc_section
      TYPE(section_vals_type), POINTER                   :: gauxc_section

      INTEGER                                            :: ifun
      TYPE(section_vals_type), POINTER                   :: functionals, xc_fun

      NULLIFY (gauxc_section)
      IF (.NOT. ASSOCIATED(xc_section)) RETURN

      functionals => section_vals_get_subs_vals(xc_section, "XC_FUNCTIONAL")
      IF (.NOT. ASSOCIATED(functionals)) RETURN

      ifun = 0
      DO
         ifun = ifun + 1
         xc_fun => section_vals_get_subs_vals2(functionals, i_section=ifun)
         IF (.NOT. ASSOCIATED(xc_fun)) EXIT
         IF (xc_fun%section%name == "GAUXC") THEN
            gauxc_section => xc_fun
            EXIT
         END IF
      END DO

   END FUNCTION get_gauxc_section

END MODULE skala_gpw_functional
