!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Explicit inverse-orbital attachment to native Quickstep density and operator consumers.
!>        Scalar semilocal GPW and mixed XC on full meshes with original full-wave USPP data.
! **************************************************************************************************
MODULE qs_upf_inverse_scf
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE cell_types,                      ONLY: cell_type
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_release,&
                                              cp_cfm_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_add,&
                                              dbcsr_copy,&
                                              dbcsr_p_type,&
                                              dbcsr_release,&
                                              dbcsr_set,&
                                              dbcsr_type
   USE cp_dbcsr_operations,             ONLY: dbcsr_allocate_matrix_set
   USE cp_fm_struct,                    ONLY: cp_fm_struct_type
   USE hfx_types,                       ONLY: hfx_type
   USE input_constants,                 ONLY: sic_none,&
                                              uspp_method_direct,&
                                              uspp_method_inverse,&
                                              uspp_method_inverse_1c
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE kpoint_transitional,             ONLY: set_1d_pointer,&
                                              set_2d_pointer
   USE kpoint_types,                    ONLY: kpoint_native_operator_clear
   USE lebedev,                         ONLY: deallocate_lebedev_grids,&
                                              init_lebedev_grids
   USE mathlib,                         ONLY: invmat
   USE particle_types,                  ONLY: particle_type
   USE paw_types,                       ONLY: paw_potential_type
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_grid_types,                   ONLY: pw_grid_type
   USE pw_methods,                      ONLY: pw_integrate_function,&
                                              pw_scale,&
                                              pw_transfer
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_kpoint_density_state,         ONLY: kpoint_density_state_initialize
   USE qs_kpoint_operators,             ONLY: kpoint_operator_set_native,&
                                              kpoint_slot_ks,&
                                              kpoint_slot_s,&
                                              kpoint_slot_t
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_did_change,&
                                              qs_ks_env_type,&
                                              set_ks_env
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_set,&
                                              qs_rho_type
   USE qs_upf_inverse,                  ONLY: qs_upf_inverse_atomic_operators
   USE qs_upf_inverse_augmentation,     ONLY: qs_upf_inverse_augmentation_prepare
   USE qs_upf_inverse_background,       ONLY: qs_upf_inverse_background_collocate,&
                                              qs_upf_inverse_background_gaussian_geometry,&
                                              qs_upf_inverse_background_operator,&
                                              qs_upf_inverse_background_prepare
   USE qs_upf_inverse_core,             ONLY: qs_upf_inverse_core_collocate
   USE qs_upf_inverse_grid,             ONLY: qs_upf_inverse_grid_collocate,&
                                              qs_upf_inverse_grid_integrate,&
                                              qs_upf_inverse_grid_integrate_dense,&
                                              qs_upf_inverse_grid_metric,&
                                              qs_upf_inverse_grid_prepare,&
                                              qs_upf_inverse_grid_response,&
                                              qs_upf_inverse_matrix_gather,&
                                              qs_upf_inverse_matrix_scatter
   USE qs_upf_inverse_kpoints,          ONLY: qs_upf_inverse_kpoints_collocate,&
                                              qs_upf_inverse_kpoints_density,&
                                              qs_upf_inverse_kpoints_integrate,&
                                              qs_upf_inverse_kpoints_metric,&
                                              qs_upf_inverse_kpoints_prepare,&
                                              qs_upf_inverse_kpoints_release_operators
   USE qs_upf_inverse_paw,              ONLY: qs_upf_inverse_paw_core_prepare,&
                                              qs_upf_inverse_paw_evaluate,&
                                              qs_upf_inverse_paw_init,&
                                              qs_upf_inverse_paw_response
   USE qs_upf_inverse_types,            ONLY: qs_upf_inverse_augmentation_type,&
                                              qs_upf_inverse_grid_type,&
                                              qs_upf_inverse_scf_release,&
                                              qs_upf_inverse_scf_type
   USE qs_uspp_types,                   ONLY: qs_uspp_init,&
                                              qs_uspp_release,&
                                              qs_uspp_type
   USE upf_inverse,                     ONLY: upf_inverse_fields_init,&
                                              upf_inverse_kinetic_init
   USE upf_local_integrals,             ONLY: upf_local_fourier,&
                                              upf_local_init
   USE upf_radial_fields,               ONLY: upf_radial_fields_type
   USE upf_types,                       ONLY: upf_nlcc_kinetic
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: qs_upf_inverse_scf_prepare, qs_upf_inverse_scf_active
   PUBLIC :: qs_upf_inverse_scf_paw_composite
   PUBLIC :: qs_upf_inverse_scf_collocate, qs_upf_inverse_scf_integrate, qs_upf_inverse_scf_core
   PUBLIC :: qs_upf_inverse_scf_integrate_kinetic
   PUBLIC :: qs_upf_inverse_scf_initialize, qs_upf_inverse_scf_reset
   PUBLIC :: qs_upf_inverse_scf_install_metric
   PUBLIC :: qs_upf_inverse_scf_paw_evaluate, qs_upf_inverse_scf_add_paw
   PUBLIC :: qs_upf_inverse_scf_add_coulomb
   PUBLIC :: qs_upf_inverse_scf_integrate_background
   PUBLIC :: qs_upf_inverse_scf_traces
   PUBLIC :: qs_upf_inverse_scf_core_response
   PUBLIC :: qs_upf_inverse_scf_response
CONTAINS

! **************************************************************************************************
!> \brief Explicitly install one inverse frame, its charge metric, and core operator in native KS.
!>        Never inferred from full-wave availability; caller supplies the fixed image set.
!> \param ks_env Prepared native KS environment.
!> \param grid Native physical-cell PW grid.
!> \param images Fixed physical-cell images for the Gaussian and radial fields.
!> \param info Factorization status; nonzero leaves the reference KS environment unchanged.
!> \param atom_local Optional nonoverlapping atomic inverse with native Gaussian kinetic integrals.
!> \param direct Use original projectors and soft orbitals without reconstruction.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_scf_prepare(ks_env, grid, images, info, atom_local, direct)
      TYPE(qs_ks_env_type), INTENT(INOUT)                :: ks_env
      TYPE(pw_grid_type), POINTER                        :: grid
      INTEGER, INTENT(IN)                                :: images(:, :)
      INTEGER, INTENT(OUT)                               :: info
      LOGICAL, INTENT(IN), OPTIONAL                      :: atom_local, direct

      CHARACTER(LEN=256)                                 :: core_reason
      INTEGER                                            :: ikind, img
      LOGICAL                                            :: core_available, do_kpoints, use_atomic
      REAL(dp), ALLOCATABLE                              :: core_tau(:)
      TYPE(cp_cfm_type), ALLOCATABLE                     :: native_metric(:)
      TYPE(dbcsr_p_type), POINTER                        :: kinetic_k(:, :), metric(:), &
                                                            metric_k(:, :), overlap(:), &
                                                            overlap_k(:, :)
      TYPE(dbcsr_type)                                   :: imaginary
      TYPE(dft_control_type), POINTER                    :: control
      TYPE(qs_kind_type), POINTER                        :: kinds(:)
      TYPE(qs_upf_inverse_scf_type), POINTER             :: previous, state
      TYPE(qs_uspp_type), POINTER                        :: uspp

      CALL get_ks_env(ks_env, upf_inverse=previous, matrix_s_kp=overlap_k, uspp=uspp, &
                      dft_control=control, do_kpoints=do_kpoints, qs_kind_set=kinds)
      CPASSERT(.NOT. ASSOCIATED(previous))
      CPASSERT(ASSOCIATED(uspp) .AND. ASSOCIATED(grid))
      IF (.NOT. (control%qs_control%gpw .OR. control%qs_control%gapw_xc .OR. control%qs_control%gapw)) THEN
         CPABORT("Unsupported inverse SCF density representation")
      END IF
      IF (control%drho_by_collocation) THEN
         CPABORT("Inverse SCF needs spectral XC gradients")
      END IF
      IF (control%do_admm .OR. control%do_sccs .OR. control%qs_control%cdft .OR. &
        control%qs_control%rigpw .OR. control%qs_control%lrigpw .OR. control%dft_plus_u .OR. control%sic_method_id /= sic_none) THEN
         CPABORT("Unsupported inverse SCF density consumer")
      END IF
      ALLOCATE (state)
      ALLOCATE (state%paw(SIZE(kinds)))
      CALL init_lebedev_grids()
      DO ikind = 1, SIZE(kinds)
         IF (kinds(ikind)%ghost) CYCLE
         CPASSERT(ASSOCIATED(kinds(ikind)%sgp_potential))
         CPASSERT(ALLOCATED(kinds(ikind)%sgp_potential%upf_data))
         IF (kinds(ikind)%sgp_potential%upf_data%is_paw) THEN
            CALL qs_upf_inverse_paw_init(kinds(ikind)%sgp_potential%upf_data, &
                                         kinds(ikind)%sgp_potential%upf_projectors, kinds(ikind)%ngrid_ang, state%paw(ikind))
         END IF
         IF (.NOT. kinds(ikind)%sgp_potential%upf_data%is_ultrasoft .AND. &
             .NOT. kinds(ikind)%sgp_potential%upf_data%is_paw .AND. &
             TRIM(kinds(ikind)%sgp_potential%upf_data%pseudo_type) /= "NC") THEN
            CPABORT("Unsupported inverse UPF type")
         END IF
         IF (control%use_kinetic_energy_density) THEN
            IF (state%paw(ikind)%active) THEN
               IF (.NOT. state%paw(ikind)%xcgrid%core_kinetic%available) THEN
                  CPABORT(TRIM(state%paw(ikind)%xcgrid%core_kinetic%failure_reason))
               END IF
            ELSE
               CALL upf_nlcc_kinetic(kinds(ikind)%sgp_potential%upf_data, core_tau, core_available, core_reason)
               IF (.NOT. core_available) CPABORT(TRIM(core_reason))
            END IF
         END IF
      END DO
      CALL deallocate_lebedev_grids()
      IF (control%use_gauxc .AND. ANY(state%paw%active)) THEN
         IF (.NOT. control%use_kinetic_energy_density) THEN
            CPABORT("Inverse PAW learned XC requires the native kinetic-density evaluator")
         END IF
         ALLOCATE (state%paw_reconstruction(SIZE(kinds)))
         DO ikind = 1, SIZE(kinds)
            IF (.NOT. state%paw(ikind)%active) CYCLE
            CALL upf_inverse_kinetic_init(kinds(ikind)%sgp_potential%upf_data, state%paw_reconstruction(ikind))
         END DO
      END IF
      use_atomic = .FALSE.
      IF (PRESENT(atom_local)) use_atomic = atom_local
      IF (PRESENT(direct)) THEN
         IF (direct) use_atomic = .TRUE.
      END IF
      IF (use_atomic) THEN
         CALL qs_upf_inverse_background_prepare(ks_env, state%gaussian)
         CALL get_ks_env(ks_env, kinetic_kp=kinetic_k)
         CPASSERT(ASSOCIATED(kinetic_k))
         ALLOCATE (state%gaussian_kinetic(state%gaussian%nao, state%gaussian%nao, SIZE(state%gaussian%weights)))
         CALL qs_upf_inverse_background_operator(state%gaussian, kinetic_k(1, :), state%gaussian_kinetic)
      END IF
      IF (control%use_kinetic_energy_density) THEN
         ALLOCATE (state%kinetic_cache(SIZE(kinds)))
         DO ikind = 1, SIZE(kinds)
            IF (kinds(ikind)%ghost) CYCLE
            ! PAW uses smooth PS tau and a separate AE-minus-PS one-center functional.
            IF (state%paw(ikind)%active) CYCLE
            IF (.NOT. kinds(ikind)%sgp_potential%upf_data%is_ultrasoft) CYCLE
            CALL upf_inverse_kinetic_init(kinds(ikind)%sgp_potential%upf_data, state%kinetic_cache(ikind))
         END DO
      END IF
      IF (do_kpoints) THEN
         CALL qs_upf_inverse_kpoints_prepare(ks_env, grid, images, state%bloch, info, use_atomic, direct)
         IF (info /= 0) THEN
            CALL qs_upf_inverse_scf_release(state)
            RETURN
         END IF
         CALL qs_upf_inverse_kpoints_metric(state%bloch, native_metric, dense=state%metric_k)
         CALL qs_upf_inverse_kpoints_release_operators(native_metric)
         ! The image-indexed matrix set is only the direct metric's source identity.
         ! Physical matrix values are the complete native Bloch matrices, never these zeros.
         CALL dbcsr_allocate_matrix_set(metric_k, 1, SIZE(overlap_k, 2))
         DO img = 1, SIZE(overlap_k, 2)
            ALLOCATE (metric_k(1, img)%matrix)
            CALL dbcsr_copy(metric_k(1, img)%matrix, overlap_k(1, img)%matrix)
            CALL dbcsr_set(metric_k(1, img)%matrix, 0.0_dp)
         END DO
         CALL set_2d_pointer(state%metric, metric_k)
         state%enabled = .TRUE.
         CALL set_ks_env(ks_env, upf_inverse=state)
         CALL qs_upf_inverse_core_collocate(ks_env)
         CALL register_bloch(state, metric_k, kpoint_slot_s, state%metric_k)
         CALL qs_upf_inverse_scf_core(ks_env)
         CALL kpoint_density_state_initialize(state%bloch%kpoints, control%qs_control%uspp_method)
         CALL qs_ks_did_change(ks_env, s_mstruct_changed=.TRUE., rho_changed=.TRUE., potential_changed=.TRUE.)
         RETURN
      END IF
      CPASSERT(control%nimages == 1)
      CALL get_ks_env(ks_env, matrix_s=overlap)
      CALL qs_upf_inverse_grid_prepare(ks_env, grid, images, [0.0_dp, 0.0_dp, 0.0_dp], state%frame, info, use_atomic, direct)
      IF (info /= 0) THEN
         CALL qs_upf_inverse_scf_release(state)
         RETURN
      END IF
      CALL qs_upf_inverse_augmentation_prepare(ks_env, state%frame%ao, state%augmentation)
      ALLOCATE (metric(1))
      ALLOCATE (metric(1)%matrix)
      CALL qs_upf_inverse_grid_metric(state%frame, state%frame%ao%radial_q, overlap(1)%matrix, &
                                      metric(1)%matrix, imaginary)
      CALL check_real(state, metric(1)%matrix, imaginary)
      CALL dbcsr_release(imaginary)
      CALL set_1d_pointer(state%metric, metric)
      ! Retain all inverse-mediated AO fill in native density/KS allocation as well as in M.
      CALL dbcsr_copy(uspp%ao_pattern, metric(1)%matrix)
      CALL dbcsr_set(uspp%ao_pattern, 0.0_dp)
      state%enabled = .TRUE.
      CALL set_ks_env(ks_env, upf_inverse=state)
      CALL qs_upf_inverse_core_collocate(ks_env)
      CALL qs_upf_inverse_scf_core(ks_env)
      CALL qs_ks_did_change(ks_env, s_mstruct_changed=.TRUE., rho_changed=.TRUE., potential_changed=.TRUE.)
   END SUBROUTINE qs_upf_inverse_scf_prepare

! **************************************************************************************************
!> \brief Install the charge metric of a supplied density functional before native SCF starts.
!>        The same AO order and global k order as the defining inverse frames are required.
!>        Geometry changes must rebuild the functional, metric and allocation pattern together.
!> \param ks_env Freshly prepared inverse environment, before density or operator evaluations.
!> \param dense Hermitian per-k unit-potential adjoints of the complete electronic density.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_scf_install_metric(ks_env, dense)
      TYPE(qs_ks_env_type), INTENT(INOUT)                :: ks_env
      COMPLEX(dp), INTENT(IN)                            :: dense(:, :, :)

      INTEGER                                            :: k, n
      REAL(dp)                                           :: defect, tolerance
      TYPE(dbcsr_p_type), POINTER                        :: metric(:), metric_k(:, :)
      TYPE(dbcsr_type)                                   :: imaginary, real_part
      TYPE(qs_upf_inverse_scf_type), POINTER             :: state
      TYPE(qs_uspp_type), POINTER                        :: uspp

      CALL get_ks_env(ks_env, upf_inverse=state, uspp=uspp)
      CPASSERT(ASSOCIATED(state) .AND. ASSOCIATED(uspp))
      CALL check_geometry(ks_env, state)
      CPASSERT(state%density_evaluations == 0 .AND. state%operator_evaluations == 0)
      n = SIZE(dense, 1)
      CPASSERT(n > 0 .AND. SIZE(dense, 2) == n .AND. SIZE(dense, 3) > 0)
      DO k = 1, SIZE(dense, 3)
         tolerance = 128.0_dp*EPSILON(1.0_dp)*MAX(1.0_dp, MAXVAL(ABS(dense(:, :, k))))
         defect = MAXVAL(ABS(dense(:, :, k) - CONJG(TRANSPOSE(dense(:, :, k)))))
         CPASSERT(defect <= tolerance)
      END DO
      IF (ASSOCIATED(state%bloch%kpoints)) THEN
         CPASSERT(ALL(SHAPE(dense) == SHAPE(state%metric_k)))
         state%metric_k(:, :, :) = dense
         CALL get_ks_env(ks_env, matrix_metric_kp=metric_k)
         CALL register_bloch(state, metric_k, kpoint_slot_s, state%metric_k)
      ELSE
         CPASSERT(SIZE(dense, 3) == 1)
         CALL get_ks_env(ks_env, matrix_metric=metric)
         CALL qs_upf_inverse_matrix_scatter(state%frame%ao, dense(:, :, 1), metric(1)%matrix, real_part, imaginary)
         CALL check_real(state, real_part, imaginary)
         CALL dbcsr_copy(metric(1)%matrix, real_part)
         CALL dbcsr_add(uspp%ao_pattern, real_part, 1.0_dp, 1.0_dp)
         CALL dbcsr_set(uspp%ao_pattern, 0.0_dp)
         CALL dbcsr_release(real_part)
         CALL dbcsr_release(imaginary)
      END IF
      CALL qs_ks_did_change(ks_env, s_mstruct_changed=.TRUE., rho_changed=.TRUE., potential_changed=.TRUE.)
   END SUBROUTINE qs_upf_inverse_scf_install_metric

! **************************************************************************************************
!> \brief Add the complete common Coulomb AO operator once after fresh core and XC assembly.
!>        The caller must exclude the earlier PW/one-center Hartree operators from that assembly.
!>        Stored Coulomb operators must come from the current accepted native density evaluation.
!> \param ks_env Current native inverse state and newly assembled spin Hamiltonians.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_scf_add_coulomb(ks_env)
      TYPE(qs_ks_env_type), INTENT(INOUT)                :: ks_env

      INTEGER                                            :: spin
      TYPE(dbcsr_p_type), POINTER                        :: source(:, :)
      TYPE(dbcsr_type)                                   :: imaginary, real_part
      TYPE(qs_upf_inverse_scf_type), POINTER             :: state

      CALL get_ks_env(ks_env, upf_inverse=state, matrix_ks_kp=source)
      CPASSERT(ASSOCIATED(state) .AND. ASSOCIATED(source))
      CALL check_geometry(ks_env, state)
      CPASSERT(ASSOCIATED(state%coulomb))
      CPASSERT(state%coulomb%ready .AND. ALLOCATED(state%hartree_k))
      CPASSERT(state%coulomb_density_evaluation == state%density_evaluations)
      CPASSERT(state%coulomb_operator_evaluation < state%operator_evaluations)
      IF (ASSOCIATED(state%bloch%kpoints)) THEN
         CPASSERT(ALL(SHAPE(state%hartree_k) == SHAPE(state%core_k)))
         DO spin = 1, SIZE(state%ks_k, 4)
            state%ks_k(:, :, :, spin) = state%ks_k(:, :, :, spin) + state%hartree_k
         END DO
         CALL register_bloch_spin(state, source)
      ELSE
         CPASSERT(SIZE(source, 2) == 1 .AND. SIZE(state%hartree_k, 3) == 1)
         DO spin = 1, SIZE(source, 1)
            CALL qs_upf_inverse_matrix_scatter(state%frame%ao, state%hartree_k(:, :, 1), &
                                               source(spin, 1)%matrix, real_part, imaginary)
            CALL check_real(state, real_part, imaginary)
            CALL dbcsr_add(source(spin, 1)%matrix, real_part, 1.0_dp, 1.0_dp)
            CALL dbcsr_release(real_part)
            CALL dbcsr_release(imaginary)
         END DO
      END IF
      state%coulomb_operator_evaluation = state%operator_evaluations
   END SUBROUTINE qs_upf_inverse_scf_add_coulomb

! **************************************************************************************************
!> \brief Evaluate the original PAW one-center functional from the accepted native Gaussian density.
!>        All spin densities and global k weights enter the same atomic functional before differentiation.
!> \param ks_env Current inverse state and native density.
!> \param xc_section Activated native XC section, including its density and gradient cutoffs.
!> \param energies PAW Hartree, XC and frozen-core constant, each counted once per atom.
!> \param exchange Optional native HF contributions, including their physical interaction parameters.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_scf_paw_evaluate(ks_env, xc_section, energies, exchange)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(section_vals_type), POINTER                   :: xc_section
      REAL(dp), INTENT(OUT)                              :: energies(3)
      TYPE(hfx_type), INTENT(IN), OPTIONAL               :: exchange(:, :)

      COMPLEX(dp), ALLOCATABLE                           :: gamma(:, :), matrix(:, :, :)
      INTEGER                                            :: ikind, nao, nk, ns, spin
      REAL(dp)                                           :: density_cutoff, gradient_cutoff
      REAL(dp), ALLOCATABLE                              :: core_energy(:)
      TYPE(dbcsr_p_type), POINTER                        :: density(:, :)
      TYPE(dft_control_type), POINTER                    :: control
      TYPE(paw_potential_type), ALLOCATABLE              :: core_potential(:)
      TYPE(qs_kind_type), POINTER                        :: kinds(:)
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_upf_inverse_scf_type), POINTER             :: state
      TYPE(section_vals_type), POINTER                   :: functional

      energies = 0.0_dp
      CALL get_ks_env(ks_env, upf_inverse=state, rho=rho, dft_control=control, qs_kind_set=kinds)
      IF (.NOT. ASSOCIATED(state)) RETURN
      IF (.NOT. ANY(state%paw%active)) RETURN
      CALL check_geometry(ks_env, state)
      ns = control%nspins
      IF (PRESENT(exchange)) THEN
         ALLOCATE (core_energy(SIZE(kinds)), SOURCE=0.0_dp)
         ALLOCATE (core_potential(SIZE(kinds)))
         DO ikind = 1, SIZE(kinds)
            IF (.NOT. state%paw(ikind)%active) CYCLE
            ASSOCIATE (pot => kinds(ikind)%sgp_potential%upf_data, &
                       projectors => kinds(ikind)%sgp_potential%upf_projectors)
               ALLOCATE (core_potential(ikind)%matrix(projectors%nproj, projectors%nproj, 1))
               CALL qs_upf_inverse_paw_core_prepare(pot, projectors, state%paw(ikind), exchange, ns, &
                                                    core_energy(ikind), core_potential(ikind)%matrix(:, :, 1))
            END ASSOCIATE
         END DO
      END IF
      nk = 1
      IF (ASSOCIATED(state%bloch%kpoints)) THEN
         nk = SIZE(state%bloch%weights)
         nao = SIZE(state%bloch%frames(1)%ao%coefficients, 2)
      ELSE
         nao = SIZE(state%frame%ao%coefficients, 2)
         CALL qs_rho_get(rho, rho_ao_kp=density)
      END IF
      IF (ALLOCATED(state%paw_density)) DEALLOCATE (state%paw_density, state%paw_operator)
      ALLOCATE (state%paw_density(nao, nao, nk, ns), state%paw_operator(nao, nao, nk, ns))
      DO spin = 1, ns
         IF (ASSOCIATED(state%bloch%kpoints)) THEN
            CALL qs_upf_inverse_kpoints_density(state%bloch, spin, matrix)
            state%paw_density(:, :, :, spin) = matrix
         ELSE
            CALL qs_upf_inverse_matrix_gather(state%frame%ao, density(spin, 1)%matrix, gamma)
            state%paw_density(:, :, 1, spin) = gamma
         END IF
      END DO
      functional => section_vals_get_subs_vals(xc_section, "XC_FUNCTIONAL")
      CALL section_vals_val_get(xc_section, "DENSITY_CUTOFF", r_val=density_cutoff)
      CALL section_vals_val_get(xc_section, "GRADIENT_CUTOFF", r_val=gradient_cutoff)
      IF (ASSOCIATED(state%bloch%kpoints)) THEN
         CALL qs_upf_inverse_paw_evaluate(state%paw, state%bloch%frames, state%paw_density, state%bloch%weights, &
                                          functional, density_cutoff, gradient_cutoff, energies, &
                                          state%paw_operator, state%paw_potential, &
                                          semilocal_xc=.NOT. ALLOCATED(state%paw_reconstruction), &
                                          core_energy=core_energy, core_potential=core_potential)
      ELSE
         CALL qs_upf_inverse_paw_evaluate(state%paw, [state%frame], state%paw_density, [1.0_dp], &
                                          functional, density_cutoff, gradient_cutoff, energies, &
                                          state%paw_operator, state%paw_potential, &
                                          semilocal_xc=.NOT. ALLOCATED(state%paw_reconstruction), &
                                          core_energy=core_energy, core_potential=core_potential)
      END IF
      state%paw_density_evaluation = state%density_evaluations
   END SUBROUTINE qs_upf_inverse_scf_paw_evaluate

! **************************************************************************************************
!> \brief Add PAW operators once after the complete smooth and GAPW operators have been assembled.
!> \param ks_env Current inverse state and newly assembled native spin Hamiltonians.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_scf_add_paw(ks_env)
      TYPE(qs_ks_env_type), INTENT(INOUT)                :: ks_env

      INTEGER                                            :: spin
      TYPE(dbcsr_p_type), POINTER                        :: source(:, :)
      TYPE(dbcsr_type)                                   :: imaginary, real_part
      TYPE(qs_upf_inverse_scf_type), POINTER             :: state

      CALL get_ks_env(ks_env, upf_inverse=state, matrix_ks_kp=source)
      IF (.NOT. ASSOCIATED(state)) RETURN
      IF (.NOT. ALLOCATED(state%paw_operator)) RETURN
      CALL check_geometry(ks_env, state)
      CPASSERT(state%paw_density_evaluation == state%density_evaluations)
      IF (ASSOCIATED(state%bloch%kpoints)) THEN
         state%ks_k(:, :, :, :) = state%ks_k + state%paw_operator
         CALL register_bloch_spin(state, source)
      ELSE
         DO spin = 1, SIZE(source, 1)
            CALL qs_upf_inverse_matrix_scatter(state%frame%ao, state%paw_operator(:, :, 1, spin), &
                                               source(spin, 1)%matrix, real_part, imaginary)
            CALL check_real(state, real_part, imaginary)
            CALL dbcsr_add(source(spin, 1)%matrix, real_part, 1.0_dp, 1.0_dp)
            CALL dbcsr_release(real_part)
            CALL dbcsr_release(imaginary)
         END DO
      END IF
   END SUBROUTINE qs_upf_inverse_scf_add_paw

! **************************************************************************************************
!> \brief Test explicit attachment; full-wave data alone never activates this representation.
!> \param ks_env ...
!> \return ...
! **************************************************************************************************
   LOGICAL FUNCTION qs_upf_inverse_scf_active(ks_env) RESULT(active)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env

      TYPE(qs_upf_inverse_scf_type), POINTER             :: state

      CALL get_ks_env(ks_env, upf_inverse=state)
      active = ASSOCIATED(state)
      IF (active) active = state%enabled
   END FUNCTION qs_upf_inverse_scf_active

! **************************************************************************************************
!> \brief Test whether PAW supplies joint reconstructed valence and original AE core XC fields.
!> \param ks_env Native inverse environment.
!> \return True only for the explicitly prepared learned PAW composite.
! **************************************************************************************************
   LOGICAL FUNCTION qs_upf_inverse_scf_paw_composite(ks_env) RESULT(active)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env

      TYPE(qs_upf_inverse_scf_type), POINTER             :: state

      CALL get_ks_env(ks_env, upf_inverse=state)
      active = .FALSE.
      IF (.NOT. ASSOCIATED(state)) RETURN
      active = state%enabled .AND. ALLOCATED(state%paw_reconstruction)
   END FUNCTION qs_upf_inverse_scf_paw_composite

! **************************************************************************************************
!> \brief Build the original short-range local, Dion and transformed kinetic operator.
!>        The long-range ionic/core Coulomb term remains in native Poisson energy accounting.
!> \param ks_env ...
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_scf_core(ks_env)
      TYPE(qs_ks_env_type), INTENT(INOUT)                :: ks_env

      COMPLEX(dp), ALLOCATABLE                           :: atomic_t(:, :)
      INTEGER                                            :: axis
      TYPE(dbcsr_p_type), POINTER                        :: h(:), kinetic(:), overlap(:)
      TYPE(dbcsr_type)                                   :: imaginary, real_part
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_grid_type), POINTER                        :: grid
      TYPE(pw_pool_type), POINTER                        :: pool
      TYPE(pw_r3d_rs_type)                               :: local_r, unit_tau, zero_density, zeros(3)
      TYPE(qs_upf_inverse_scf_type), POINTER             :: state

      CALL get_ks_env(ks_env, upf_inverse=state)
      CPASSERT(ASSOCIATED(state))
      IF (ASSOCIATED(state%bloch%kpoints)) THEN
         CALL core_bloch(ks_env, state)
         RETURN
      END IF
      CALL get_ks_env(ks_env, upf_inverse=state, &
                      matrix_h=h, kinetic=kinetic, matrix_s=overlap, pw_env=pw_env)
      CPASSERT(ASSOCIATED(state))
      CALL check_geometry(ks_env, state)
      CALL pw_env_get(pw_env, auxbas_pw_pool=pool)
      grid => pool%pw_grid
      CALL local_r%create(grid)
      CALL unit_tau%create(grid)
      CALL zero_density%create(grid)
      unit_tau%array = 1.0_dp
      zero_density%array = 0.0_dp
      DO axis = 1, 3
         CALL zeros(axis)%create(grid)
         zeros(axis)%array = 0.0_dp
      END DO
      IF (state%frame%ao%atom_local) THEN
         CALL atomic_kinetic(state, state%frame, 1, atomic_t)
         CALL qs_upf_inverse_matrix_scatter(state%frame%ao, atomic_t, overlap(1)%matrix, real_part, imaginary)
         unit_tau%array = 0.0_dp
      ELSE
         CALL qs_upf_inverse_grid_integrate(state%frame, zero_density, zeros, unit_tau, overlap(1)%matrix, &
                                            real_part, imaginary)
      END IF
      CALL check_real(state, real_part, imaginary)
      CALL dbcsr_copy(kinetic(1)%matrix, real_part)
      CALL dbcsr_release(real_part)
      CALL dbcsr_release(imaginary)
      CALL local_potential(ks_env, grid, state%frame%ao%positions, local_r)
      CALL qs_upf_inverse_grid_integrate(state%frame, local_r, zeros, unit_tau, overlap(1)%matrix, &
                                         real_part, imaginary, projector_operator=state%frame%ao%dion, &
                                         augmentation=state%augmentation, fourier_charge=.TRUE.)
      CALL check_real(state, real_part, imaginary)
      CALL dbcsr_copy(h(1)%matrix, real_part)
      IF (state%frame%ao%atom_local) CALL dbcsr_add(h(1)%matrix, kinetic(1)%matrix, 1.0_dp, 1.0_dp)
      CALL dbcsr_release(real_part)
      CALL dbcsr_release(imaginary)
      CALL local_r%release()
      CALL unit_tau%release()
      CALL zero_density%release()
      DO axis = 1, 3
         CALL zeros(axis)%release()
      END DO
   END SUBROUTINE qs_upf_inverse_scf_core

! **************************************************************************************************
!> \brief Native spin density update with original Fourier Q and the same charge metric.
!> \param ks_env ...
!> \param rho_struct Native AO/grid density with occupations in its density matrices.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_scf_collocate(ks_env, rho_struct)
      TYPE(qs_ks_env_type), INTENT(IN), TARGET           :: ks_env
      TYPE(qs_rho_type), INTENT(INOUT)                   :: rho_struct

      COMPLEX(dp), ALLOCATABLE                           :: gamma(:, :), matrix(:, :, :)
      INTEGER                                            :: axis, k, spin
      REAL(dp), POINTER                                  :: total(:)
      TYPE(dbcsr_p_type), POINTER                        :: density(:, :)
      TYPE(pw_c1d_gs_type), POINTER                      :: rho_g(:), tau_g(:)
      TYPE(pw_r3d_rs_type)                               :: gradient(3), tau
      TYPE(pw_r3d_rs_type), POINTER                      :: rho_r(:), tau_r(:)
      TYPE(qs_ks_env_type), POINTER                      :: native_ks
      TYPE(qs_upf_inverse_scf_type), POINTER             :: state

      CALL get_ks_env(ks_env, upf_inverse=state)
      CPASSERT(ASSOCIATED(state))
      CALL check_geometry(ks_env, state)
      CALL qs_rho_get(rho_struct, rho_ao_kp=density, rho_r=rho_r, rho_g=rho_g, tot_rho_r=total, tau_r=tau_r, tau_g=tau_g)
      IF (ALLOCATED(state%kinetic_cache)) THEN
         CPASSERT(ASSOCIATED(tau_r) .AND. ASSOCIATED(tau_g))
         CPASSERT(SIZE(tau_r) == SIZE(density, 1) .AND. SIZE(tau_g) == SIZE(density, 1))
      END IF
      IF (.NOT. ASSOCIATED(state%bloch%kpoints)) THEN
         CPASSERT(SIZE(density, 2) == 1)
      END IF
      CALL tau%create(rho_r(1)%pw_grid)
      DO axis = 1, 3
         CALL gradient(axis)%create(rho_r(1)%pw_grid)
      END DO
      DO spin = 1, SIZE(density, 1)
         IF (state%gaussian_background) THEN
            CPASSERT(ASSOCIATED(state%coulomb))
            CPASSERT(state%coulomb%ready)
            IF (ASSOCIATED(state%bloch%kpoints)) THEN
               CALL qs_upf_inverse_kpoints_density(state%bloch, spin, matrix)
            ELSE
               CALL qs_upf_inverse_matrix_gather(state%frame%ao, density(spin, 1)%matrix, gamma)
               ALLOCATE (matrix(SIZE(gamma, 1), SIZE(gamma, 2), 1))
               matrix(:, :, 1) = gamma
               DEALLOCATE (gamma)
            END IF
            native_ks => ks_env
            IF (ALLOCATED(state%kinetic_cache)) THEN
               CALL qs_upf_inverse_background_collocate(state%coulomb%background, native_ks, matrix, &
                                                        rho_r(spin), rho_g(spin), &
                                                        kinetic=tau_r(spin), kinetic_gspace=tau_g(spin))
            ELSE
               CALL qs_upf_inverse_background_collocate(state%coulomb%background, native_ks, matrix, rho_r(spin), rho_g(spin))
            END IF
            ! The soft-grid integral excludes the center charge. Count the complete electronic source.
            total(spin) = 0.0_dp
            DO k = 1, SIZE(matrix, 3)
               total(spin) = total(spin) - state%coulomb%weights(k)* &
                             REAL(SUM(matrix(:, :, k)*CONJG(state%coulomb%metric(:, :, k))), dp)
            END DO
            DEALLOCATE (matrix)
            CYCLE
         ELSE IF (ASSOCIATED(state%bloch%kpoints)) THEN
            CALL qs_upf_inverse_kpoints_collocate(state%bloch, spin, rho_r(spin), gradient, tau, state%kinetic_cache)
         ELSE
            CALL qs_upf_inverse_grid_collocate(state%frame, density(spin, 1)%matrix, rho_r(spin), gradient, tau, &
                                               augmentation=state%augmentation, fourier_charge=.TRUE., kinetic=state%kinetic_cache)
         END IF
         IF (ALLOCATED(state%kinetic_cache)) THEN
            tau_r(spin)%array = tau%array
            CALL pw_transfer(tau_r(spin), tau_g(spin))
         END IF
         CALL pw_transfer(rho_r(spin), rho_g(spin))
         total(spin) = pw_integrate_function(rho_r(spin), isign=-1)
      END DO
      CALL qs_rho_set(rho_struct, rho_r_valid=.TRUE., rho_g_valid=.TRUE.)
      IF (ALLOCATED(state%kinetic_cache)) THEN
         CALL qs_rho_set(rho_struct, tau_r_valid=.TRUE., tau_g_valid=.TRUE.)
      END IF
      CALL tau%release()
      DO axis = 1, 3
         CALL gradient(axis)%release()
      END DO
      state%density_evaluations = state%density_evaluations + 1
   END SUBROUTINE qs_upf_inverse_scf_collocate

! **************************************************************************************************
!> \brief Integrate the native folded Hartree/XC potential through the same inverse density adjoint.
!> \param ks_env ...
!> \param potential Native potential, already including the spectral XC gradient adjoint.
!> \param hamiltonian Native spin Hamiltonian, incremented.
!> \param weighted Whether potential already contains the native dvol factor.
!> \param spin ...
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_scf_integrate(ks_env, potential, hamiltonian, weighted, spin)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: potential
      TYPE(dbcsr_type), INTENT(INOUT)                    :: hamiltonian
      LOGICAL, INTENT(IN)                                :: weighted
      INTEGER, INTENT(IN), OPTIONAL                      :: spin

      COMPLEX(dp), ALLOCATABLE                           :: dense(:, :, :)
      INTEGER                                            :: axis
      TYPE(cp_cfm_type), ALLOCATABLE                     :: native(:)
      TYPE(dbcsr_p_type), POINTER                        :: source(:, :)
      TYPE(dbcsr_type)                                   :: imaginary, real_part
      TYPE(pw_r3d_rs_type)                               :: unweighted, zero_tau, zeros(3)
      TYPE(qs_upf_inverse_scf_type), POINTER             :: state

      CALL get_ks_env(ks_env, upf_inverse=state)
      CPASSERT(ASSOCIATED(state))
      CALL check_geometry(ks_env, state)
      CALL unweighted%create(potential%pw_grid)
      unweighted%array = potential%array
      IF (weighted) CALL pw_scale(unweighted, 1.0_dp/potential%pw_grid%dvol)
      CALL zero_tau%create(potential%pw_grid)
      zero_tau%array = 0.0_dp
      DO axis = 1, 3
         CALL zeros(axis)%create(potential%pw_grid)
         zeros(axis)%array = 0.0_dp
      END DO
      IF (ASSOCIATED(state%bloch%kpoints)) THEN
         CPASSERT(PRESENT(spin))
         CALL qs_upf_inverse_kpoints_integrate(state%bloch, unweighted, zeros, zero_tau, native, dense)
         CALL qs_upf_inverse_kpoints_release_operators(native)
         state%ks_k(:, :, :, spin) = state%core_k + dense
         IF (ALLOCATED(state%xc%operators)) THEN
            state%ks_k(:, :, :, spin) = state%ks_k(:, :, :, spin) + state%xc%operators(:, :, :, spin)
         END IF
         CALL get_ks_env(ks_env, matrix_ks_kp=source)
         CALL register_bloch_spin(state, source)
      ELSE
         CALL qs_upf_inverse_grid_integrate(state%frame, unweighted, zeros, zero_tau, hamiltonian, &
                                            real_part, imaginary, augmentation=state%augmentation, fourier_charge=.TRUE.)
         CALL check_real(state, real_part, imaginary)
         CALL dbcsr_add(hamiltonian, real_part, 1.0_dp, 1.0_dp)
         CALL dbcsr_release(real_part)
         CALL dbcsr_release(imaginary)
         IF (ALLOCATED(state%xc%operators)) THEN
            CPASSERT(PRESENT(spin))
            CALL qs_upf_inverse_matrix_scatter(state%frame%ao, state%xc%operators(:, :, 1, spin), &
                                               hamiltonian, real_part, imaginary)
            CALL check_real(state, real_part, imaginary)
            CALL dbcsr_add(hamiltonian, real_part, 1.0_dp, 1.0_dp)
            CALL dbcsr_release(real_part)
            CALL dbcsr_release(imaginary)
         END IF
      END IF
      CALL unweighted%release()
      CALL zero_tau%release()
      DO axis = 1, 3
         CALL zeros(axis)%release()
      END DO
      state%operator_evaluations = state%operator_evaluations + 1
   END SUBROUTINE qs_upf_inverse_scf_integrate

! **************************************************************************************************
!> \brief Add the native kinetic XC adjoint through the same augmented tau map used by collocation.
!>        The bare USPP soft kinetic operator remains part of the separate core Hamiltonian.
!> \param ks_env Live inverse environment after scalar-potential integration.
!> \param potential Native potential conjugate to tau, without the grid volume element.
!> \param hamiltonian Native Gamma spin Hamiltonian, incremented.
!> \param spin Native spin index, including Bloch direct operators.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_scf_integrate_kinetic(ks_env, potential, hamiltonian, spin)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: potential
      TYPE(dbcsr_type), INTENT(INOUT)                    :: hamiltonian
      INTEGER, INTENT(IN)                                :: spin

      COMPLEX(dp), ALLOCATABLE                           :: dense(:, :, :)
      INTEGER                                            :: axis
      TYPE(cp_cfm_type), ALLOCATABLE                     :: native(:)
      TYPE(dbcsr_p_type), POINTER                        :: source(:, :)
      TYPE(dbcsr_type)                                   :: imaginary, real_part
      TYPE(pw_r3d_rs_type)                               :: zero, zeros(3)
      TYPE(qs_upf_inverse_scf_type), POINTER             :: state

      CALL get_ks_env(ks_env, upf_inverse=state)
      CPASSERT(ASSOCIATED(state))
      CPASSERT(ALLOCATED(state%kinetic_cache))
      CALL check_geometry(ks_env, state)
      CALL zero%create(potential%pw_grid)
      zero%array = 0.0_dp
      DO axis = 1, 3
         CALL zeros(axis)%create(potential%pw_grid)
         zeros(axis)%array = 0.0_dp
      END DO
      IF (ASSOCIATED(state%bloch%kpoints)) THEN
         CALL qs_upf_inverse_kpoints_integrate(state%bloch, zero, zeros, potential, native, dense, state%kinetic_cache)
         CALL qs_upf_inverse_kpoints_release_operators(native)
         state%ks_k(:, :, :, spin) = state%ks_k(:, :, :, spin) + dense
         CALL get_ks_env(ks_env, matrix_ks_kp=source)
         CALL register_bloch_spin(state, source)
      ELSE
         CALL qs_upf_inverse_grid_integrate(state%frame, zero, zeros, potential, hamiltonian, &
                                            real_part, imaginary, kinetic=state%kinetic_cache)
         CALL check_real(state, real_part, imaginary)
         CALL dbcsr_add(hamiltonian, real_part, 1.0_dp, 1.0_dp)
         CALL dbcsr_release(real_part)
         CALL dbcsr_release(imaginary)
      END IF
      CALL zero%release()
      DO axis = 1, 3
         CALL zeros(axis)%release()
      END DO
   END SUBROUTINE qs_upf_inverse_scf_integrate_kinetic

! **************************************************************************************************
!> \brief Insert the Gaussian soft-grid XC adjoint and the matching direct atomic XC operator.
!>        Dense per-k operators contain neither mesh weights nor a Hartree contribution.
!> \param ks_env Current inverse native environment.
!> \param dense Gaussian background transpose of the folded XC potential.
!> \param hamiltonian Native Gamma spin matrix, incremented; Bloch values use the direct provider.
!> \param spin Native spin index.
!> \param additive Add a background kinetic operator to the assembled Hamiltonian without repeating atomic XC.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_scf_integrate_background(ks_env, dense, hamiltonian, spin, additive)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      COMPLEX(dp), INTENT(IN)                            :: dense(:, :, :)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: hamiltonian
      INTEGER, INTENT(IN)                                :: spin
      LOGICAL, INTENT(IN), OPTIONAL                      :: additive

      COMPLEX(dp), ALLOCATABLE                           :: operator(:, :, :)
      LOGICAL                                            :: add_to_current
      TYPE(dbcsr_p_type), POINTER                        :: source(:, :)
      TYPE(dbcsr_type)                                   :: imaginary, real_part
      TYPE(qs_upf_inverse_scf_type), POINTER             :: state

      CALL get_ks_env(ks_env, upf_inverse=state)
      CPASSERT(ASSOCIATED(state))
      CPASSERT(state%gaussian_background .AND. ASSOCIATED(state%coulomb))
      CALL check_geometry(ks_env, state)
      add_to_current = .FALSE.
      IF (PRESENT(additive)) add_to_current = additive
      ALLOCATE (OPERATOR, SOURCE=dense)
      IF (.NOT. add_to_current .AND. ALLOCATED(state%xc%operators)) THEN
         OPERATOR(:, :, :) = OPERATOR + state%xc%operators(:, :, :, spin)
      END IF
      IF (ASSOCIATED(state%bloch%kpoints)) THEN
         CPASSERT(ALL(SHAPE(OPERATOR) == SHAPE(state%core_k)))
         IF (add_to_current) THEN
            state%ks_k(:, :, :, spin) = state%ks_k(:, :, :, spin) + OPERATOR
         ELSE
            state%ks_k(:, :, :, spin) = state%core_k + OPERATOR
         END IF
         CALL get_ks_env(ks_env, matrix_ks_kp=source)
         CALL register_bloch_spin(state, source)
      ELSE
         CPASSERT(SIZE(OPERATOR, 3) == 1)
         CALL qs_upf_inverse_matrix_scatter(state%frame%ao, OPERATOR(:, :, 1), hamiltonian, real_part, imaginary)
         CALL check_real(state, real_part, imaginary)
         CALL dbcsr_add(hamiltonian, real_part, 1.0_dp, 1.0_dp)
         CALL dbcsr_release(real_part)
         CALL dbcsr_release(imaginary)
      END IF
      IF (.NOT. add_to_current) state%operator_evaluations = state%operator_evaluations + 1
   END SUBROUTINE qs_upf_inverse_scf_integrate_background

! **************************************************************************************************
!> \brief Joint geometry response of the inverse core Hamiltonian and its charge metric.
!>        Differentiates the inverse kinetic correction, original Dion, short-range potentials and Q.
!>        The general frame includes soft kinetic energy. For an atomic frame, the native Gaussian
!>        kinetic response is contracted separately by qs_upf_inverse_scf_response.
!>        Gaussian ionic/core electrostatics remains in the native Poisson and core-energy terms.
!> \param ks_env Native kinds and their original radial potentials; geometry is not modified.
!> \param frame Prepared inverse fields at the geometry being differentiated.
!> \param augmentation Original charge kernels shared with the density map.
!> \param grid Native physical-cell grid of frame.
!> \param atom_response Cartesian response of every native atom.
!> \param cell_response Lattice-vector response, including its reciprocal and volume effects.
!> \param reduced_k_response Reduced Bloch-vector response, zero for cell stress.
!> \param core_response New replicated core-Hamiltonian response without occupation/k weights.
!> \param metric_response New response of the same radial-Q metric used by inverse SCF.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_scf_core_response(ks_env, frame, augmentation, grid, atom_response, cell_response, &
                                               reduced_k_response, core_response, metric_response)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(qs_upf_inverse_grid_type), INTENT(IN)         :: frame
      TYPE(qs_upf_inverse_augmentation_type), &
         INTENT(INOUT)                                   :: augmentation
      TYPE(pw_grid_type), POINTER                        :: grid
      REAL(dp), INTENT(IN)                               :: atom_response(:, :), &
                                                            cell_response(3, 3), &
                                                            reduced_k_response(3)
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: core_response(:, :), &
                                                            metric_response(:, :)

      COMPLEX(dp), ALLOCATABLE                           :: atomic_m(:, :), atomic_t(:, :), &
                                                            dm(:, :), dt(:, :)
      INTEGER                                            :: axis, info
      REAL(dp)                                           :: cell_inverse(3, 3), strain(3, 3)
      TYPE(pw_r3d_rs_type)                               :: local_r, local_response, unit_tau, zero, &
                                                            zeros(3)

      cell_inverse = frame%ao%cell
      CALL invmat(cell_inverse, info)
      CPASSERT(info == 0)
      strain = MATMUL(cell_response, cell_inverse)
      CALL local_r%create(grid)
      CALL local_response%create(grid)
      CALL unit_tau%create(grid)
      CALL zero%create(grid)
      unit_tau%array = 1.0_dp
      IF (frame%ao%atom_local) unit_tau%array = 0.0_dp
      zero%array = 0.0_dp
      DO axis = 1, 3
         CALL zeros(axis)%create(grid)
         zeros(axis)%array = 0.0_dp
      END DO
      CALL local_potential(ks_env, grid, frame%ao%positions, local_r, local_response, atom_response, strain)
      CALL qs_upf_inverse_grid_response(frame, augmentation, frame%ao%radial_q, local_r, zeros, unit_tau, &
                                        atom_response, cell_response, reduced_k_response, local_response, zeros, zero, &
                                        core_response, metric_response, frame%ao%dion)
      IF (frame%ao%atom_local) THEN
         CALL qs_upf_inverse_atomic_operators(frame%ao, frame%ao%radial_q, atomic_t, atomic_m, dt, dm, &
                                              atom_response, cell_response, reduced_k_response)
         core_response(:, :) = core_response + dt
      END IF
      CALL local_r%release()
      CALL local_response%release()
      CALL unit_tau%release()
      CALL zero%release()
      DO axis = 1, 3
         CALL zeros(axis)%release()
      END DO
   END SUBROUTINE qs_upf_inverse_scf_core_response

! **************************************************************************************************
!> \brief Contract joint inverse operator and metric responses with native stationary P and W.
!>        W must have been prepared by the native energy-weighted-density service after SCF.
!>        Each k-point mesh weight is applied once, after replicated AO contractions.
!>        Ionic/core, Poisson-kernel and explicit XC measure responses remain separate energy terms.
!> \param ks_env Converged native inverse SCF state and prepared energy-weighted densities.
!> \param vdensity Unweighted primitive spin potentials conjugate to density.
!> \param vgradient Unweighted primitive spin potentials conjugate to Cartesian gradients.
!> \param vtau Unweighted primitive spin potentials conjugate to smooth kinetic density.
!> \param atom_response Cartesian response of all native atoms.
!> \param cell_response Lattice-vector response, with fixed fractional grid rows.
!> \param reduced_k_response Response of all actual reduced k vectors, or one column at Gamma.
!> \param response Stationary electronic response, including the metric constraint.
!> \param metric_response Signed metric contribution, already included in response.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_scf_response(ks_env, vdensity, vgradient, vtau, atom_response, cell_response, &
                                          reduced_k_response, response, metric_response)
      TYPE(qs_ks_env_type), INTENT(INOUT), TARGET        :: ks_env
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: vdensity(:), vgradient(:, :), vtau(:)
      REAL(dp), INTENT(IN)                               :: atom_response(:, :), &
                                                            cell_response(3, 3), &
                                                            reduced_k_response(:, :)
      REAL(dp), INTENT(OUT)                              :: response, metric_response

      COMPLEX(dp), ALLOCATABLE :: core(:, :, :), density(:, :, :), energy_weighted(:, :, :), &
         matrix(:, :), metric(:, :, :), metric_part(:, :)
      INTEGER                                            :: active_potential, axis, component, ik, &
                                                            info, nk, nspin
      REAL(dp)                                           :: cell_inverse(3, 3), &
                                                            gaussian_virial(3, 3), paw_response, &
                                                            strain(3, 3), weight
      REAL(dp), ALLOCATABLE                              :: gaussian_gradients(:, :), &
                                                            relative_motion(:, :)
      TYPE(dbcsr_p_type), POINTER                        :: density_gamma(:, :), energy_gamma(:)
      TYPE(pw_grid_type), POINTER                        :: grid
      TYPE(pw_r3d_rs_type)                               :: zero, zeros(3)
      TYPE(qs_ks_env_type), POINTER                      :: native_ks
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_upf_inverse_scf_type), POINTER             :: state

      CALL get_ks_env(ks_env, upf_inverse=state, rho=rho)
      CPASSERT(ASSOCIATED(state))
      CALL check_geometry(ks_env, state)
      nspin = SIZE(vdensity)
      CPASSERT(ALL(SHAPE(vgradient) == [3, nspin]) .AND. SIZE(vtau) == nspin)
      grid => vdensity(1)%pw_grid
      CPASSERT(ASSOCIATED(grid))
      CALL zero%create(grid)
      zero%array = 0.0_dp
      DO axis = 1, 3
         CALL zeros(axis)%create(grid)
         zeros(axis)%array = 0.0_dp
      END DO
      nk = 1
      IF (ASSOCIATED(state%bloch%kpoints)) nk = SIZE(state%bloch%weights)
      CPASSERT(ALL(SHAPE(reduced_k_response) == [3, nk]))
      IF (ALLOCATED(state%gaussian_kinetic)) THEN
         CPASSERT(ALL(reduced_k_response == 0.0_dp))
         native_ks => ks_env
         ALLOCATE (gaussian_gradients(3, SIZE(atom_response, 2)), relative_motion(3, SIZE(atom_response, 2)))
         IF (ASSOCIATED(state%bloch%kpoints)) THEN
            cell_inverse = state%bloch%frames(1)%ao%cell
         ELSE
            cell_inverse = state%frame%ao%cell
         END IF
         CALL invmat(cell_inverse, info)
         CPASSERT(info == 0)
         strain = MATMUL(cell_response, cell_inverse)
         IF (ASSOCIATED(state%bloch%kpoints)) THEN
            relative_motion(:, :) = atom_response - MATMUL(strain, state%bloch%frames(1)%ao%positions)
         ELSE
            relative_motion(:, :) = atom_response - MATMUL(strain, state%frame%ao%positions)
         END IF
      END IF
      DO ik = 1, nk
         IF (ASSOCIATED(state%bloch%kpoints)) THEN
            CALL qs_upf_inverse_scf_core_response(ks_env, state%bloch%frames(ik), state%bloch%augmentation, grid, &
                                                  atom_response, cell_response, reduced_k_response(:, ik), matrix, metric_part)
         ELSE
            CALL qs_upf_inverse_scf_core_response(ks_env, state%frame, state%augmentation, grid, atom_response, &
                                                  cell_response, reduced_k_response(:, ik), matrix, metric_part)
         END IF
         IF (ik == 1) ALLOCATE (core(SIZE(matrix, 1), SIZE(matrix, 2), nk), metric(SIZE(matrix, 1), SIZE(matrix, 2), nk))
         core(:, :, ik) = matrix
         metric(:, :, ik) = metric_part
      END DO
      response = 0.0_dp
      metric_response = 0.0_dp
      IF (.NOT. ASSOCIATED(state%bloch%kpoints)) THEN
         CALL qs_rho_get(rho, rho_ao_kp=density_gamma)
         CALL get_ks_env(ks_env, matrix_w=energy_gamma)
         CPASSERT(ASSOCIATED(density_gamma) .AND. ASSOCIATED(energy_gamma))
         CPASSERT(SIZE(density_gamma, 1) == nspin .AND. SIZE(energy_gamma) == nspin)
         ALLOCATE (density(SIZE(core, 1), SIZE(core, 2), 1), energy_weighted(SIZE(core, 1), SIZE(core, 2), 1))
      END IF
      DO axis = 1, nspin
         active_potential = MERGE(1, 0, ANY(vdensity(axis)%array /= 0.0_dp) .OR. ANY(vtau(axis)%array /= 0.0_dp))
         DO component = 1, 3
            IF (ANY(vgradient(component, axis)%array /= 0.0_dp)) active_potential = 1
         END DO
         ! Empty or zero local slabs must follow the same collective response path as other ranks.
         IF (ASSOCIATED(state%bloch%kpoints)) THEN
            CALL state%bloch%frames(1)%ao%point_comm%max(active_potential)
         ELSE
            CALL state%frame%ao%point_comm%max(active_potential)
         END IF
         IF (ASSOCIATED(state%bloch%kpoints)) THEN
            CALL qs_upf_inverse_kpoints_density(state%bloch, axis, density)
            CALL qs_upf_inverse_kpoints_density(state%bloch, axis, energy_weighted, energy_weighted=.TRUE.)
         ELSE
            CALL qs_upf_inverse_matrix_gather(state%frame%ao, density_gamma(axis, 1)%matrix, matrix)
            density(:, :, 1) = matrix
            CALL qs_upf_inverse_matrix_gather(state%frame%ao, energy_gamma(axis)%matrix, matrix)
            energy_weighted(:, :, 1) = matrix
         END IF
         IF (ALLOCATED(state%gaussian_kinetic)) THEN
            CALL qs_upf_inverse_background_gaussian_geometry(state%gaussian, native_ks, density, .TRUE., &
                                                             gaussian_gradients, gaussian_virial)
            response = response + SUM(gaussian_gradients*relative_motion) - SUM(gaussian_virial*strain)
         END IF
         DO ik = 1, nk
            weight = 1.0_dp
            IF (ASSOCIATED(state%bloch%kpoints)) THEN
               weight = state%bloch%weights(ik)
            END IF
            IF (active_potential == 0) THEN
               matrix(:, :) = CMPLX(0.0_dp, 0.0_dp, dp)
            ELSE IF (ASSOCIATED(state%bloch%kpoints)) THEN
               CALL qs_upf_inverse_grid_response(state%bloch%frames(ik), state%bloch%augmentation, &
                                                 state%bloch%frames(ik)%ao%radial_q, &
                                                 vdensity(axis), vgradient(:, axis), vtau(axis), atom_response, cell_response, &
                                                 reduced_k_response(:, ik), zero, zeros, zero, matrix, metric_part, &
                                                 kinetic=state%kinetic_cache)
            ELSE
               CALL qs_upf_inverse_grid_response(state%frame, state%augmentation, state%frame%ao%radial_q, &
                                                 vdensity(axis), vgradient(:, axis), vtau(axis), atom_response, cell_response, &
                                                 reduced_k_response(:, ik), zero, zeros, zero, matrix, metric_part, &
                                                 kinetic=state%kinetic_cache)
            END IF
            response = response + weight*REAL(SUM(density(:, :, ik)*TRANSPOSE(core(:, :, ik) + matrix)), dp)
            metric_response = metric_response - weight*REAL(SUM(energy_weighted(:, :, ik)*TRANSPOSE(metric(:, :, ik))), dp)
         END DO
      END DO
      IF (ALLOCATED(state%paw_potential)) THEN
         CPASSERT(state%paw_density_evaluation == state%density_evaluations)
         IF (ASSOCIATED(state%bloch%kpoints)) THEN
            CALL qs_upf_inverse_paw_response(state%bloch%frames, state%paw_density, state%bloch%weights, &
                                             state%paw_potential, atom_response, cell_response, reduced_k_response, paw_response)
         ELSE
            CALL qs_upf_inverse_paw_response([state%frame], state%paw_density, [1.0_dp], &
                                             state%paw_potential, atom_response, cell_response, reduced_k_response, paw_response)
         END IF
         response = response + paw_response
      END IF
      response = response + metric_response
      CALL zero%release()
      DO axis = 1, 3
         CALL zeros(axis)%release()
      END DO
   END SUBROUTINE qs_upf_inverse_scf_response

! **************************************************************************************************
!> \brief Original short-range ionic potential and its optional fixed-fractional-row response.
!>        The same Fourier-to-real projection is used for the core operator and its derivative.
!> \param ks_env Native kind mapping and prepared original UPF radial potentials.
!> \param grid Native physical-cell grid.
!> \param positions Cartesian centers at the current inverse geometry.
!> \param local_r Allocated unweighted scalar potential, overwritten.
!> \param local_response Optional allocated response, including the inverse-volume derivative.
!> \param atom_response Required with local_response, Cartesian center displacements.
!> \param strain Required with local_response, cell_response times inverse cell.
! **************************************************************************************************
   SUBROUTINE local_potential(ks_env, grid, positions, local_r, local_response, atom_response, strain)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(pw_grid_type), POINTER                        :: grid
      REAL(dp), INTENT(IN)                               :: positions(:, :)
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: local_r
      TYPE(pw_r3d_rs_type), INTENT(INOUT), OPTIONAL      :: local_response
      REAL(dp), INTENT(IN), OPTIONAL                     :: atom_response(:, :), strain(3, 3)

      COMPLEX(dp), ALLOCATABLE                           :: center_derivative(:, :), ionic_g(:), &
                                                            reciprocal_derivative(:, :)
      INTEGER                                            :: atom, axis, ikind, ng
      REAL(dp)                                           :: volume_response
      REAL(dp), ALLOCATABLE                              :: g_response(:, :)
      TYPE(particle_type), POINTER                       :: atoms(:)
      TYPE(pw_c1d_gs_type)                               :: local_g, response_g
      TYPE(qs_kind_type), POINTER                        :: kinds(:)

      CALL get_ks_env(ks_env, qs_kind_set=kinds, particle_set=atoms)
      CPASSERT(ALL(SHAPE(positions) == [3, SIZE(atoms)]))
      CPASSERT(ASSOCIATED(local_r%pw_grid, grid))
      CALL local_g%create(grid)
      ng = SIZE(local_g%array)
      ALLOCATE (ionic_g(ng))
      local_g%array = CMPLX(0.0_dp, 0.0_dp, dp)
      IF (PRESENT(local_response)) THEN
         CPASSERT(PRESENT(atom_response) .AND. PRESENT(strain))
         CPASSERT(ALL(SHAPE(atom_response) == SHAPE(positions)))
         CPASSERT(ASSOCIATED(local_response%pw_grid, grid))
         CALL response_g%create(grid)
         response_g%array = CMPLX(0.0_dp, 0.0_dp, dp)
         ALLOCATE (center_derivative(3, ng), reciprocal_derivative(3, ng), g_response(3, ng))
         g_response(:, :) = -MATMUL(TRANSPOSE(strain), grid%g)
         volume_response = strain(1, 1) + strain(2, 2) + strain(3, 3)
      END IF
      DO atom = 1, SIZE(atoms)
         ikind = atoms(atom)%atomic_kind%kind_number
         IF (kinds(ikind)%ghost) CYCLE
         ASSOCIATE (potential => kinds(ikind)%sgp_potential)
            IF (potential%upf_local%radial%degree < 0) THEN
               CALL upf_local_init(potential%upf_data, potential%alpha_core_charge, 0, potential%upf_local)
            END IF
            IF (PRESENT(local_response)) THEN
               CALL upf_local_fourier(potential%upf_local, grid%g, positions(:, atom), grid%vol, ionic_g, &
                                      center_derivative, reciprocal_derivative)
               response_g%array = response_g%array - volume_response*ionic_g
               DO axis = 1, 3
                  response_g%array = response_g%array + atom_response(axis, atom)*center_derivative(axis, :) + &
                                     g_response(axis, :)*reciprocal_derivative(axis, :)
               END DO
            ELSE
               CALL upf_local_fourier(potential%upf_local, grid%g, positions(:, atom), grid%vol, ionic_g)
            END IF
         END ASSOCIATE
         local_g%array = local_g%array + ionic_g
      END DO
      CALL pw_transfer(local_g, local_r)
      IF (PRESENT(local_response)) THEN
         CALL pw_transfer(response_g, local_response)
         CALL response_g%release()
      END IF
      CALL local_g%release()
   END SUBROUTINE local_potential

! **************************************************************************************************
!> \brief Reject a stale inverse geometry instead of mixing it with rebuilt Gaussian operators.
!> \param ks_env ...
!> \param state ...
! **************************************************************************************************
   SUBROUTINE check_geometry(ks_env, state)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(qs_upf_inverse_scf_type), INTENT(IN)          :: state

      INTEGER                                            :: atom
      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_type), POINTER                       :: atoms(:)

      CALL get_ks_env(ks_env, cell=cell, particle_set=atoms)
      CPASSERT(state%enabled)
      IF (ASSOCIATED(state%bloch%kpoints)) THEN
         CPASSERT(ALL(cell%hmat == state%bloch%frames(1)%ao%cell))
         DO atom = 1, SIZE(atoms)
            CPASSERT(ALL(atoms(atom)%r == state%bloch%frames(1)%ao%positions(:, atom)))
         END DO
      ELSE
         CPASSERT(ALL(cell%hmat == state%frame%ao%cell))
         DO atom = 1, SIZE(atoms)
            CPASSERT(ALL(atoms(atom)%r == state%frame%ao%positions(:, atom)))
         END DO
      END IF
   END SUBROUTINE check_geometry

! **************************************************************************************************
!> \brief Verify that the Gamma operator can be stored in the real native KS representation.
!> \param state ...
!> \param real_part ...
!> \param imaginary ...
! **************************************************************************************************
   SUBROUTINE check_real(state, real_part, imaginary)
      TYPE(qs_upf_inverse_scf_type), INTENT(IN)          :: state
      TYPE(dbcsr_type), INTENT(IN)                       :: real_part, imaginary

      COMPLEX(dp), ALLOCATABLE                           :: dense(:, :)
      REAL(dp)                                           :: defect, tolerance

      CALL qs_upf_inverse_matrix_gather(state%frame%ao, real_part, dense, imaginary)
      defect = MAXVAL(ABS(AIMAG(dense)))
      tolerance = 128.0_dp*EPSILON(1.0_dp)*MAX(1.0_dp, MAXVAL(ABS(dense)))
      CPASSERT(defect <= tolerance)
   END SUBROUTINE check_real

! **************************************************************************************************
!> \brief Install spin-free complete Bloch operators through the direct native provider.
!> \param state ...
!> \param source Image-indexed source identity, not physical matrix values.
!> \param slot ...
!> \param dense ...
! **************************************************************************************************
   SUBROUTINE register_bloch(state, source, slot, dense)
      TYPE(qs_upf_inverse_scf_type), INTENT(IN)          :: state
      TYPE(dbcsr_p_type), POINTER                        :: source(:, :)
      INTEGER, INTENT(IN)                                :: slot
      COMPLEX(dp), INTENT(IN)                            :: dense(:, :, :)

      INTEGER                                            :: col, ik, local, row
      TYPE(cp_cfm_type), ALLOCATABLE                     :: native(:, :)
      TYPE(cp_fm_struct_type), POINTER                   :: layout

      CPASSERT(SIZE(source, 1) == 1)
      ALLOCATE (native(SIZE(state%bloch%kpoints%kp_env), 1))
      DO local = 1, SIZE(native, 1)
         ik = state%bloch%kp_range(1) + local - 1
         layout => state%bloch%kpoints%kp_env(local)%kpoint_env%pmat(1, 1)%matrix_struct
         CALL cp_cfm_create(native(local, 1), layout)
         native(local, 1)%local_data = CMPLX(0.0_dp, 0.0_dp, dp)
         DO col = 1, SIZE(layout%col_indices)
            DO row = 1, SIZE(layout%row_indices)
               native(local, 1)%local_data(row, col) = dense(layout%row_indices(row), layout%col_indices(col), ik)
            END DO
         END DO
      END DO
      CALL kpoint_operator_set_native(state%bloch%kpoints, source, slot, native)
      DO local = 1, SIZE(native, 1)
         CALL cp_cfm_release(native(local, 1))
      END DO
   END SUBROUTINE register_bloch

! **************************************************************************************************
!> \brief Refresh the complete spin-resolved native KS provider after folding Hartree/XC.
!> \param state ...
!> \param source ...
! **************************************************************************************************
   SUBROUTINE register_bloch_spin(state, source)
      TYPE(qs_upf_inverse_scf_type), INTENT(IN)          :: state
      TYPE(dbcsr_p_type), POINTER                        :: source(:, :)

      INTEGER                                            :: col, ik, local, row, spin
      TYPE(cp_cfm_type), ALLOCATABLE                     :: native(:, :)
      TYPE(cp_fm_struct_type), POINTER                   :: layout

      CPASSERT(SIZE(source, 1) == SIZE(state%ks_k, 4))
      ALLOCATE (native(SIZE(state%bloch%kpoints%kp_env), SIZE(source, 1)))
      DO spin = 1, SIZE(source, 1)
         DO local = 1, SIZE(native, 1)
            ik = state%bloch%kp_range(1) + local - 1
            layout => state%bloch%kpoints%kp_env(local)%kpoint_env%pmat(1, spin)%matrix_struct
            CALL cp_cfm_create(native(local, spin), layout)
            native(local, spin)%local_data = CMPLX(0.0_dp, 0.0_dp, dp)
            DO col = 1, SIZE(layout%col_indices)
               DO row = 1, SIZE(layout%row_indices)
                  native(local, spin)%local_data(row, col) = &
                     state%ks_k(layout%row_indices(row), layout%col_indices(col), ik, spin)
               END DO
            END DO
         END DO
      END DO
      CALL kpoint_operator_set_native(state%bloch%kpoints, source, kpoint_slot_ks, native)
      DO spin = 1, SIZE(native, 2)
         DO local = 1, SIZE(native, 1)
            CALL cp_cfm_release(native(local, spin))
         END DO
      END DO
   END SUBROUTINE register_bloch_spin

! **************************************************************************************************
!> \brief Native Gaussian kinetic energy plus the exact finite atomic reconstruction correction.
!>        The original density quadrature continues to define the installed charge metric.
!> \param state Native Gaussian kinetic operators captured before the inverse core replaced them.
!> \param frame Nonoverlapping atomic inverse at the current geometry and k point.
!> \param ik Global k-point index in the captured native Gaussian operators.
!> \param matrix New complete transformed kinetic operator.
! **************************************************************************************************
   SUBROUTINE atomic_kinetic(state, frame, ik, matrix)
      TYPE(qs_upf_inverse_scf_type), INTENT(IN)          :: state
      TYPE(qs_upf_inverse_grid_type), INTENT(IN)         :: frame
      INTEGER, INTENT(IN)                                :: ik
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: matrix(:, :)

      COMPLEX(dp), ALLOCATABLE                           :: metric(:, :)

      CPASSERT(frame%ao%atom_local .AND. ALLOCATED(state%gaussian_kinetic))
      CALL qs_upf_inverse_atomic_operators(frame%ao, frame%ao%radial_q, matrix, metric)
      matrix(:, :) = matrix + state%gaussian_kinetic(:, :, ik)
   END SUBROUTINE atomic_kinetic

! **************************************************************************************************
!> \brief Build fixed original core/kinetic Bloch operators on the native physical-cell grid.
!> \param ks_env ...
!> \param state ...
! **************************************************************************************************
   SUBROUTINE core_bloch(ks_env, state)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(qs_upf_inverse_scf_type), INTENT(INOUT)       :: state

      COMPLEX(dp), ALLOCATABLE                           :: matrix(:, :)
      INTEGER                                            :: axis, ik, n, nk, spin
      TYPE(dbcsr_p_type), POINTER                        :: kinetic(:, :), source(:, :)
      TYPE(dft_control_type), POINTER                    :: control
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_grid_type), POINTER                        :: grid
      TYPE(pw_pool_type), POINTER                        :: pool
      TYPE(pw_r3d_rs_type)                               :: local_r, unit_tau, zero, zeros(3)

      CALL check_geometry(ks_env, state)
      CALL get_ks_env(ks_env, dft_control=control, pw_env=pw_env, kinetic_kp=kinetic, matrix_ks_kp=source)
      ! KS storage is recreated after the core step when the AO structure changes.
      ! Do not retain an operator provider tied to that released source identity.
      IF (.NOT. ASSOCIATED(source)) THEN
         CALL kpoint_native_operator_clear(state%bloch%kpoints%native_operators(kpoint_slot_ks))
      END IF
      CALL pw_env_get(pw_env, auxbas_pw_pool=pool)
      grid => pool%pw_grid
      CALL local_r%create(grid)
      CALL zero%create(grid)
      CALL unit_tau%create(grid)
      zero%array = 0.0_dp
      unit_tau%array = 1.0_dp
      IF (state%bloch%frames(1)%ao%atom_local) unit_tau%array = 0.0_dp
      DO axis = 1, 3
         CALL zeros(axis)%create(grid)
         zeros(axis)%array = 0.0_dp
      END DO
      CALL local_potential(ks_env, grid, state%bloch%frames(1)%ao%positions, local_r)
      n = state%bloch%nao
      nk = SIZE(state%bloch%weights)
      IF (.NOT. ALLOCATED(state%core_k)) THEN
         ALLOCATE (state%core_k(n, n, nk), state%kinetic_k(n, n, nk), state%ks_k(n, n, nk, control%nspins))
      END IF
      DO ik = 1, nk
         IF (state%bloch%frames(ik)%ao%atom_local) THEN
            CALL atomic_kinetic(state, state%bloch%frames(ik), ik, matrix)
         ELSE
            CALL qs_upf_inverse_grid_integrate_dense(state%bloch%frames(ik), zero, zeros, unit_tau, matrix)
         END IF
         state%kinetic_k(:, :, ik) = matrix
         CALL qs_upf_inverse_grid_integrate_dense(state%bloch%frames(ik), local_r, zeros, unit_tau, matrix, &
                                                  projector_operator=state%bloch%frames(ik)%ao%dion, &
                                                  augmentation=state%bloch%augmentation, fourier_charge=.TRUE.)
         state%core_k(:, :, ik) = matrix
         IF (state%bloch%frames(ik)%ao%atom_local) THEN
            state%core_k(:, :, ik) = state%core_k(:, :, ik) + state%kinetic_k(:, :, ik)
         END IF
      END DO
      DO spin = 1, SIZE(state%ks_k, 4)
         state%ks_k(:, :, :, spin) = state%core_k
      END DO
      CALL register_bloch(state, kinetic, kpoint_slot_t, state%kinetic_k)
      IF (ASSOCIATED(source)) CALL register_bloch_spin(state, source)
      CALL local_r%release()
      CALL zero%release()
      CALL unit_tau%release()
      DO axis = 1, 3
         CALL zeros(axis)%release()
      END DO
   END SUBROUTINE core_bloch

! **************************************************************************************************
!> \brief Trace fixed core and kinetic operators against the same accepted inverse Bloch density.
!> \param ks_env ...
!> \param core ...
!> \param kinetic ...
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_scf_traces(ks_env, core, kinetic)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      REAL(dp), INTENT(OUT)                              :: core, kinetic

      COMPLEX(dp), ALLOCATABLE                           :: density(:, :, :)
      INTEGER                                            :: ik, spin
      TYPE(qs_upf_inverse_scf_type), POINTER             :: state

      CALL get_ks_env(ks_env, upf_inverse=state)
      CPASSERT(ASSOCIATED(state))
      CALL check_geometry(ks_env, state)
      CPASSERT(ASSOCIATED(state%bloch%kpoints))
      core = 0.0_dp
      kinetic = 0.0_dp
      DO spin = 1, SIZE(state%bloch%kpoints%native_density, 2)
         CALL qs_upf_inverse_kpoints_density(state%bloch, spin, density)
         DO ik = 1, SIZE(state%bloch%weights)
            core = core + state%bloch%weights(ik)*REAL(SUM(density(:, :, ik)*TRANSPOSE(state%core_k(:, :, ik))), dp)
            kinetic = kinetic + state%bloch%weights(ik)*REAL(SUM(density(:, :, ik)*TRANSPOSE(state%kinetic_k(:, :, ik))), dp)
         END DO
      END DO
   END SUBROUTINE qs_upf_inverse_scf_traces

! **************************************************************************************************
!> \brief Drop the metric, fields and allocation pattern together before a native geometry rebuild.
!> \param ks_env Native environment. The input kind selection remains available for rebuilding.
!> \param atom_local Optional previous formulation, retained across geometry updates by native SCF.
!> \param was_enabled Optional previous inverse activation status.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_scf_reset(ks_env, atom_local, was_enabled)
      TYPE(qs_ks_env_type), INTENT(INOUT)                :: ks_env
      LOGICAL, INTENT(OUT), OPTIONAL                     :: atom_local, was_enabled

      TYPE(qs_upf_inverse_scf_type), POINTER             :: state
      TYPE(qs_uspp_type), POINTER                        :: uspp

      CALL get_ks_env(ks_env, upf_inverse=state, uspp=uspp)
      IF (PRESENT(atom_local)) atom_local = .FALSE.
      IF (PRESENT(was_enabled)) was_enabled = .FALSE.
      IF (ASSOCIATED(state)) THEN
         IF (PRESENT(atom_local)) atom_local = ALLOCATED(state%gaussian_kinetic)
         IF (PRESENT(was_enabled)) was_enabled = state%enabled
      END IF
      CALL qs_upf_inverse_scf_release(state)
      CALL qs_uspp_release(uspp)
      CALL set_ks_env(ks_env, upf_inverse=state, uspp=uspp)
   END SUBROUTINE qs_upf_inverse_scf_reset

! **************************************************************************************************
!> \brief Activate the scalar inverse USPP representation before allocating SCF matrices.
!>        Image bounds cover the native basis radii and original radial support on the current cell.
!>        A geometry update rebuilds both the charge metric and all density/operator fields.
!>        USPP_METHOD selects soft orbitals, the general inverse, or the nonoverlapping atomic inverse.
!>        DIRECT and INVERSE_1C use native Gaussian kinetic integrals. INVERSE retains intercenter terms.
!> \param ks_env Ordinary Gaussian core and PW environment at the new geometry.
!> \param atom_local Optional explicit formulation. True rejects overlapping original supports.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_scf_initialize(ks_env, atom_local)
      TYPE(qs_ks_env_type), INTENT(INOUT)                :: ks_env
      LOGICAL, INTENT(IN), OPTIONAL                      :: atom_local

      INTEGER                                            :: atom, axis, ikind, image, info, ix, iy, &
                                                            iz, lower(3), upper(3)
      INTEGER, ALLOCATABLE                               :: images(:, :)
      LOGICAL                                            :: use_atomic, use_direct
      REAL(dp)                                           :: fractional(3), lower_real(3), radius, &
                                                            reach, upper_real(3)
      REAL(dp), ALLOCATABLE                              :: radial_radius(:)
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: control
      TYPE(gto_basis_set_type), POINTER                  :: basis
      TYPE(particle_type), POINTER                       :: atoms(:)
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: pool
      TYPE(qs_kind_type), POINTER                        :: kinds(:)
      TYPE(qs_uspp_type), POINTER                        :: uspp
      TYPE(upf_radial_fields_type)                       :: beta, delta

      CALL get_ks_env(ks_env, qs_kind_set=kinds, particle_set=atoms, cell=cell, pw_env=pw_env, dft_control=control)
      CALL get_ks_env(ks_env, uspp=uspp)
      CALL qs_uspp_init(uspp, kinds, atoms)
      CALL set_ks_env(ks_env, uspp=uspp)
      IF (.NOT. ASSOCIATED(uspp)) RETURN
      SELECT CASE (control%qs_control%uspp_method)
      CASE (uspp_method_direct, uspp_method_inverse, uspp_method_inverse_1c)
         use_direct = control%qs_control%uspp_method == uspp_method_direct
         use_atomic = control%qs_control%uspp_method /= uspp_method_inverse
      CASE DEFAULT
         CPABORT("Unknown USPP orbital representation")
      END SELECT
      IF (PRESENT(atom_local)) use_atomic = atom_local
      IF (use_direct) use_atomic = .TRUE.
      ALLOCATE (radial_radius(SIZE(kinds)))
      radial_radius = 0.0_dp
      DO ikind = 1, SIZE(kinds)
         IF (kinds(ikind)%ghost) CYCLE
         IF (.NOT. ASSOCIATED(kinds(ikind)%sgp_potential)) THEN
            CPABORT("Inverse USPP mixed potentials are not implemented")
         END IF
         IF (.NOT. ALLOCATED(kinds(ikind)%sgp_potential%upf_data)) THEN
            CPABORT("Inverse USPP mixed potentials are not implemented")
         END IF
         IF (kinds(ikind)%sgp_potential%upf_data%number_of_proj == 0) CYCLE
         CALL upf_inverse_fields_init(kinds(ikind)%sgp_potential%upf_data, beta, delta, use_direct)
         radial_radius(ikind) = MAX(beta%radius, delta%radius)
      END DO
      IF (control%qs_control%gapw .AND. ANY(cell%perd == 0)) THEN
         CPABORT("Inverse GAPW currently requires fully periodic electrostatics")
      END IF
      lower_real = 0.0_dp
      upper_real = 0.0_dp
      DO atom = 1, SIZE(atoms)
         ikind = atoms(atom)%atomic_kind%kind_number
         CALL get_qs_kind(kinds(ikind), basis_set=basis, basis_type='ORB')
         radius = basis%kind_radius
         IF (.NOT. kinds(ikind)%ghost) THEN
            radius = MAX(radius, radial_radius(ikind))
         END IF
         fractional = MATMUL(cell%h_inv, atoms(atom)%r)
         DO axis = 1, 3
            IF (cell%perd(axis) == 0) CYCLE
            reach = radius*SQRT(SUM(cell%h_inv(axis, :)**2))
            ! One additional cell keeps exponentially small Gaussian tails away from image boundaries.
            lower_real(axis) = MIN(lower_real(axis), -fractional(axis) - reach - 1.0_dp)
            upper_real(axis) = MAX(upper_real(axis), 1.0_dp - fractional(axis) + reach + 1.0_dp)
         END DO
      END DO
      lower = FLOOR(lower_real)
      upper = CEILING(upper_real)
      ALLOCATE (images(3, PRODUCT(upper - lower + 1)))
      image = 0
      DO iz = lower(3), upper(3)
         DO iy = lower(2), upper(2)
            DO ix = lower(1), upper(1)
               image = image + 1
               images(:, image) = [ix, iy, iz]
            END DO
         END DO
      END DO
      CALL pw_env_get(pw_env, auxbas_pw_pool=pool)
      CALL qs_upf_inverse_scf_prepare(ks_env, pool%pw_grid, images, info, use_atomic, use_direct)
      IF (info == -2) THEN
         CPABORT("Overlapping USPP reconstruction regions. Use USPP_METHOD INVERSE")
      END IF
      IF (info /= 0) CPABORT("Inverse USPP reconstruction is singular or ill-conditioned")
   END SUBROUTINE qs_upf_inverse_scf_initialize

END MODULE qs_upf_inverse_scf
