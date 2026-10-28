!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
!> \brief Fixed-geometry native electronic operator response for semilocal, Skala and RI hybrid functionals.
MODULE qs_upf_inverse_scf_response
   USE cell_types,                      ONLY: cell_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_create,&
                                              dbcsr_p_type,&
                                              dbcsr_release
   USE hfx_ri,                          ONLY: hfx_ri_update_ks
   USE hfx_ri_kp,                       ONLY: hfx_ri_update_ks_kp
   USE hfx_types,                       ONLY: hfx_type
   USE input_constants,                 ONLY: sic_none,&
                                              xc_vdw_fun_nonloc
   USE input_cp2k_hfx,                  ONLY: ri_pmat
   USE input_section_types,             ONLY: section_get_ival,&
                                              section_get_rval,&
                                              section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_get_subs_vals2,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_comm_type,&
                                              mp_para_env_type
   USE particle_types,                  ONLY: particle_type
   USE paw_types,                       ONLY: paw_potential_type
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_grid_types,                   ONLY: pw_grid_type
   USE pw_grids,                        ONLY: pw_grid_compare
   USE pw_methods,                      ONLY: pw_derive,&
                                              pw_transfer
   USE pw_poisson_types,                ONLY: pw_poisson_analytic,&
                                              pw_poisson_mt,&
                                              pw_poisson_multipole,&
                                              pw_poisson_periodic,&
                                              pw_poisson_type
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_upf_inverse_core,             ONLY: qs_upf_inverse_core_kinetic
   USE qs_upf_inverse_coulomb,          ONLY: qs_upf_inverse_coulomb_operator_response
   USE qs_upf_inverse_grid,             ONLY: qs_upf_inverse_matrix_gather,&
                                              qs_upf_inverse_matrix_scatter
   USE qs_upf_inverse_paw,              ONLY: qs_upf_inverse_paw_evaluate,&
                                              qs_upf_inverse_paw_xc_kernel
   USE qs_upf_inverse_scf_xc,           ONLY: qs_upf_inverse_scf_xc_model_response,&
                                              qs_upf_inverse_scf_xc_source_adjoint,&
                                              qs_upf_inverse_scf_xc_source_direction
   USE qs_upf_inverse_types,            ONLY: qs_upf_inverse_scf_type
   USE qs_upf_inverse_xc,               ONLY: qs_upf_inverse_xc_field_kernel,&
                                              qs_upf_inverse_xc_semilocal_kernel_available
   USE rs_pw_interface,                 ONLY: pw_real_transfer
   USE skala_density_constraints,       ONLY: skala_density_constraints_count,&
                                              skala_density_constraints_get,&
                                              skala_density_constraints_type,&
                                              skala_density_descriptor
   USE skala_gpw_functional,            ONLY: native_skala_gapw_atom_composite_requested,&
                                              native_skala_gapw_composite_reference,&
                                              native_skala_spin_floors,&
                                              skala_gpw_eval,&
                                              xc_section_uses_native_skala_evaluator
   USE xc_derivatives,                  ONLY: xc_functionals_get_needs,&
                                              xc_functionals_spin_floors
   USE xc_input_constants,              ONLY: xc_deriv_pw,&
                                              xc_rho_no_smooth
   USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: qs_upf_inverse_scf_electronic_response
   PUBLIC :: qs_upf_inverse_scf_density_tangent_available

CONTAINS

! **************************************************************************************************
!> \brief Recognize native positive spin-density floors with a linear density constraint.
!>        Their selected and default covectors differ only in the density normal. A kernel
!>        from the default branch can therefore be used between tangent projections, together
!>        with the current selected KS operator and the physical constraint curvature.
!>        Native model density clamps have the same common tangent. Kinetic and nonlinear
!>        selector boundaries retain numerical response.
!> \param qs_env Original evaluated inverse environment and its registered physical constraints.
!> \return Whether the selected surfaces admit this native density-tangent construction.
! **************************************************************************************************
   FUNCTION qs_upf_inverse_scf_density_tangent_available(qs_env) RESULT(available)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL                                            :: available

      INTEGER                                            :: i
      INTEGER, ALLOCATABLE                               :: descriptors(:), domains(:), &
                                                            points(:, :), spins(:)
      LOGICAL                                            :: active, kinetic_boundary, model
      REAL(dp)                                           :: kinetic_probe, kinetic_tolerance, probe, &
                                                            tolerance
      REAL(dp), ALLOCATABLE                              :: floors(:), targets(:), weights(:)
      TYPE(dft_control_type), POINTER                    :: control
      TYPE(mp_comm_type)                                 :: group
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(section_vals_type), POINTER                   :: component, functional, xc
      TYPE(skala_density_constraints_type), POINTER      :: constraints
      TYPE(xc_rho_cflags_type)                           :: needs

      available = .FALSE.
      CALL get_qs_env(qs_env, ks_env=ks_env, dft_control=control, para_env=para_env)
      IF (control%nspins < 1 .OR. control%nspins > 2) RETURN
      CALL get_ks_env(ks_env, skala_constraints=constraints)
      IF (.NOT. ASSOCIATED(constraints)) RETURN
      IF (skala_density_constraints_count(constraints) == 0) RETURN
      xc => section_vals_get_subs_vals(qs_env%input, "DFT%XC")
      model = xc_section_uses_native_skala_evaluator(xc)
      functional => section_vals_get_subs_vals(xc, "XC_FUNCTIONAL")
      IF (model) THEN
         ! The complete native model kernel is used only for its own density clamps.
         i = 0
         DO
            i = i + 1
            component => section_vals_get_subs_vals2(functional, i_section=i)
            IF (.NOT. ASSOCIATED(component)) EXIT
            CALL section_vals_val_get(component, "_SECTION_PARAMETERS_", l_val=active)
            IF (active .AND. component%section%name /= "GAUXC") RETURN
         END DO
         CALL group%set_handle(para_env%get_handle())
         ALLOCATE (floors, SOURCE=native_skala_spin_floors(xc, group))
      ELSE
         needs = xc_functionals_get_needs(functional, control%nspins == 2, kinetic_density_boundary=kinetic_boundary)
         IF (kinetic_boundary .AND. (needs%tau .OR. needs%tau_spin)) RETURN
         ALLOCATE (floors, &
                   SOURCE=xc_functionals_spin_floors(functional, control%nspins == 2, section_get_rval(xc, "DENSITY_CUTOFF")))
      END IF
      CALL skala_density_constraints_get(constraints, points, spins, weights, probe, tolerance, &
                                         descriptors=descriptors, domains=domains, targets=targets, &
                                         kinetic_probe=kinetic_probe, kinetic_tolerance=kinetic_tolerance)
      IF (ANY(descriptors /= skala_density_descriptor)) RETURN
      IF (control%nspins == 1) THEN
         IF (ANY(spins /= 0)) RETURN
      ELSE
         IF (ANY(spins < 1) .OR. ANY(spins > 2)) RETURN
      END IF
      DO i = 1, SIZE(targets)
         IF (.NOT. model .AND. targets(i) <= 0.0_dp) RETURN
         IF (ALL(floors /= targets(i))) RETURN
      END DO
      available = .TRUE.
   END FUNCTION qs_upf_inverse_scf_density_tangent_available

! **************************************************************************************************
!> \brief Apply the original electronic KS kernel to a joint spin/k occupation-matrix direction.
!>        Combine the original global XC response, exact atomic-minus-background response,
!>        the actual native source transposes, native Hartree and PAW one-center Hartree. Original
!>        local, nonlocal, kinetic and frozen-core one-body operators have zero density response.
!>        Their orbital-chart curvature is supplied separately by the caller using the actual KS
!>        operator. Fixed cores remain in the unperturbed XC arguments and have zero direction.
!>        All work fields are owned. Accepted densities, operators and Coulomb caches are unchanged.
!>        Selected nonsmooth surfaces need their separate common-tangent and normal responses.
!> \param qs_env Evaluated native inverse environment at fixed geometry.
!> \param matrices_direction Replicated unweighted Hermitian AO directions (AO,AO,k,spin).
!> \param operators_direction New complete electronic operator directions, without mesh weights.
!> \param available Whether every density-dependent contribution has a supplied kernel.
!> \param density_tangent Use only between projections onto recognized native density-floor tangents.
!>        This is not an unconstrained Hessian on a nonsmooth surface.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_scf_electronic_response(qs_env, matrices_direction, operators_direction, available, density_tangent)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      COMPLEX(dp), INTENT(IN)                            :: matrices_direction(:, :, :, :)
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT) :: operators_direction(:, :, :, :)
      LOGICAL, INTENT(OUT)                               :: available
      LOGICAL, INTENT(IN), OPTIONAL                      :: density_tangent

      COMPLEX(dp), ALLOCATABLE                           :: atomic(:, :, :, :), &
                                                            exchange(:, :, :, :), &
                                                            hartree(:, :, :), paw(:, :, :, :)
      INTEGER                                            :: n, nk, ns, selected, spin
      LOGICAL                                            :: atomic_xc, do_hfx, fine, model, &
                                                            rho_valid, tangent, tau_g_valid, &
                                                            tau_valid, use_tau, valid
      REAL(dp)                                           :: energies(3), energy, virial(3, 3)
      REAL(dp), ALLOCATABLE                              :: kweights(:)
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: control
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), POINTER                       :: particles(:)
      TYPE(paw_potential_type), ALLOCATABLE              :: paw_potential(:)
      TYPE(pw_c1d_gs_type)                               :: core_g, work_g
      TYPE(pw_c1d_gs_type), POINTER                      :: density_g(:), direction_g(:), &
                                                            native_core_g, native_g(:), &
                                                            native_tau_g(:)
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_grid_type), POINTER                        :: grid
      TYPE(pw_poisson_type), POINTER                     :: poisson
      TYPE(pw_pool_type), POINTER                        :: native_pool, pool
      TYPE(pw_r3d_rs_type)                               :: transferred_core
      TYPE(pw_r3d_rs_type), ALLOCATABLE                  :: combined_rho(:), combined_tau(:)
      TYPE(pw_r3d_rs_type), POINTER :: core_r, density(:), direction(:), first_rho(:), &
         first_tau(:), kinetic(:), kinetic_direction(:), native_r(:), native_tau(:), &
         native_weights, response_rho(:), response_tau(:), weights
      TYPE(pw_r3d_rs_type), TARGET                       :: transferred_weights
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_upf_inverse_scf_type), POINTER             :: state
      TYPE(section_vals_type), POINTER                   :: functional, hf, xc
      TYPE(skala_density_constraints_type), POINTER      :: constraints, response_constraints
      TYPE(xc_rho_cflags_type)                           :: needs

      available = .FALSE.
      tangent = .FALSE.
      IF (PRESENT(density_tangent)) tangent = density_tangent
      CALL get_qs_env(qs_env, ks_env=ks_env, rho=rho, dft_control=control, cell=cell, &
                      particle_set=particles, pw_env=pw_env, para_env=para_env)
      CALL get_ks_env(ks_env, upf_inverse=state, skala_constraints=constraints, &
                      rho_nlcc=core_r, rho_nlcc_g=native_core_g, xcint_weights=native_weights)
      IF (.NOT. ASSOCIATED(state)) RETURN
      IF (.NOT. state%enabled) RETURN
      xc => section_vals_get_subs_vals(qs_env%input, "DFT%XC")
      model = xc_section_uses_native_skala_evaluator(xc)
      IF (.NOT. model) THEN
         IF (.NOT. qs_upf_inverse_xc_semilocal_kernel_available(xc, control%nspins == 2)) RETURN
         IF (section_get_ival(xc, 'XC_GRID%XC_DERIV') /= xc_deriv_pw) RETURN
         IF (section_get_ival(xc, 'XC_GRID%XC_SMOOTH_RHO') /= xc_rho_no_smooth) RETURN
      END IF
      functional => section_vals_get_subs_vals(xc, 'XC_FUNCTIONAL')
      needs = xc_functionals_get_needs(functional, lsd=control%nspins == 2, calc_potential=.TRUE.)
      use_tau = model .OR. needs%tau .OR. needs%tau_spin
      IF (use_tau) THEN
         IF (.NOT. ALLOCATED(state%kinetic_cache)) RETURN
      END IF
      hf => section_vals_get_subs_vals(xc, "HF")
      CALL section_vals_get(hf, explicit=do_hfx)
      IF (control%do_admm .OR. control%sic_method_id /= sic_none) RETURN
      IF (control%dft_plus_u .OR. control%do_sccs .OR. control%do_fde .OR. control%do_pcc) RETURN
      IF (control%apply_embed_pot .OR. control%apply_dmfet_pot .OR. control%apply_external_density .OR. &
          control%apply_external_vxc .OR. control%qs_control%cdft .OR. control%correct_surf_dip .OR. qs_env%qmmm) RETURN
      IF (section_get_ival(xc, "VDW_POTENTIAL%POTENTIAL_TYPE") == xc_vdw_fun_nonloc) RETURN
      IF (native_skala_gapw_composite_reference(xc) .OR. native_skala_gapw_atom_composite_requested(xc)) RETURN
      selected = 0
      IF (ASSOCIATED(constraints)) selected = skala_density_constraints_count(constraints)
      CALL para_env%max(selected)
      IF (tangent) THEN
         IF (.NOT. qs_upf_inverse_scf_density_tangent_available(qs_env)) RETURN
      ELSE IF (selected /= 0) THEN
         RETURN
      END IF
      CALL qs_rho_get(rho, rho_r=native_r, rho_g=native_g, tau_r=native_tau, tau_g=native_tau_g, &
                      rho_g_valid=rho_valid, tau_r_valid=tau_valid, tau_g_valid=tau_g_valid)
      IF (.NOT. rho_valid) RETURN
      IF (use_tau) THEN
         IF (.NOT. tau_valid) RETURN
      END IF
      CPASSERT(ASSOCIATED(native_r) .AND. ASSOCIATED(native_g))
      IF (use_tau) THEN
         CPASSERT(ASSOCIATED(native_tau) .AND. ASSOCIATED(native_tau_g))
      END IF
      CALL pw_env_get(pw_env, xc_pw_pool=pool, auxbas_pw_pool=native_pool, poisson_env=poisson)
      IF (.NOT. ASSOCIATED(poisson%green_fft)) RETURN
      ! The response below uses the original Fourier influence function. Other
      ! Poisson solvers retain the existing numerical orbital curvature.
      SELECT CASE (poisson%parameters%solver)
      CASE (pw_poisson_periodic, pw_poisson_analytic, pw_poisson_mt, pw_poisson_multipole)
      CASE DEFAULT
         RETURN
      END SELECT
      IF (poisson%parameters%has_dielectric .OR. poisson%green_fft%p3m .OR. poisson%green_fft%sr_screening) RETURN
      ns = SIZE(native_r)
      n = SIZE(matrices_direction, 1)
      nk = SIZE(matrices_direction, 3)
      valid = ALL(SHAPE(matrices_direction) == [n, n, nk, ns])
      CPASSERT(valid)
      IF (do_hfx) THEN
         CALL native_exchange_response(qs_env, state, matrices_direction, exchange, valid)
         IF (.NOT. valid) RETURN
      END IF
      fine = .NOT. pw_grid_compare(native_pool%pw_grid, pool%pw_grid)
      grid => native_r(1)%pw_grid
      IF (fine) grid => pool%pw_grid
      atomic_xc = control%qs_control%gapw .OR. control%qs_control%gapw_xc .OR. &
                  ALLOCATED(state%paw_reconstruction)
      IF (atomic_xc) THEN
         CPASSERT(ASSOCIATED(state%xc_grid, grid))
         CPASSERT(ALLOCATED(state%xc_matrices))
      END IF
      ALLOCATE (density(ns), density_g(ns), kinetic(ns), direction(ns), direction_g(ns), kinetic_direction(ns))
      ALLOCATE (combined_rho(ns), combined_tau(ns))
      CALL work_g%create(grid)
      DO spin = 1, ns
         CALL density(spin)%create(grid)
         CALL density_g(spin)%create(grid)
         CALL kinetic(spin)%create(grid)
         kinetic(spin)%array = 0.0_dp
         CALL direction(spin)%create(grid)
         CALL direction_g(spin)%create(grid)
         CALL kinetic_direction(spin)%create(grid)
         CALL combined_rho(spin)%create(grid)
         CALL combined_tau(spin)%create(grid)
         CALL pw_real_transfer(native_g(spin), density_g(spin))
         IF (fine) THEN
            CALL pw_transfer(density_g(spin), density(spin))
            CALL pw_transfer(density(spin), density_g(spin))
            IF (use_tau) THEN
               IF (tau_g_valid) THEN
                  CALL pw_real_transfer(native_tau_g(spin), work_g)
               ELSE
                  BLOCK
                     TYPE(pw_c1d_gs_type) :: tau_aux
                     CALL tau_aux%create(native_tau(spin)%pw_grid)
                     CALL pw_transfer(native_tau(spin), tau_aux)
                     CALL pw_real_transfer(tau_aux, work_g)
                     CALL tau_aux%release()
                  END BLOCK
               END IF
               CALL pw_transfer(work_g, kinetic(spin))
            END IF
         ELSE
            density(spin)%array = native_r(spin)%array
            IF (use_tau) kinetic(spin)%array = native_tau(spin)%array
         END IF
      END DO
      IF (ASSOCIATED(core_r)) THEN
         CPASSERT(ASSOCIATED(native_core_g))
         CALL transferred_core%create(grid)
         CALL core_g%create(grid)
         CALL pw_real_transfer(native_core_g, core_g)
         IF (fine) THEN
            CALL pw_transfer(core_g, transferred_core)
            CALL pw_transfer(transferred_core, core_g)
         ELSE
            transferred_core%array = core_r%array
         END IF
         DO spin = 1, ns
            density(spin)%array = density(spin)%array + transferred_core%array
            density_g(spin)%array = density_g(spin)%array + core_g%array
         END DO
         CALL transferred_core%release()
         CALL core_g%release()
      END IF
      IF (use_tau) CALL qs_upf_inverse_core_kinetic(ks_env, kinetic, 1.0_dp)
      weights => native_weights
      IF (fine .AND. ASSOCIATED(native_weights)) THEN
         CALL transferred_weights%create(grid)
         CALL transfer_real(native_weights, transferred_weights)
         weights => transferred_weights
      END IF
      CALL qs_upf_inverse_scf_xc_source_direction(ks_env, matrices_direction, direction, kinetic_direction)
      DO spin = 1, ns
         CALL pw_transfer(direction(spin), direction_g(spin))
      END DO
      NULLIFY (first_rho, first_tau, response_rho, response_tau)
      IF (model) THEN
         response_constraints => constraints
         ! Keep the selected KS covector. Only the electronic kernel uses the default
         ! clamp branch between the caller's two physical tangent projections.
         IF (tangent) NULLIFY (response_constraints)
         CALL skala_gpw_eval(first_rho, first_tau, energy, density, density_g, kinetic, xc, weights, pool, &
                             particles, cell, .FALSE., virial, density_constraints=response_constraints, &
                             rho_direction=direction, rho_direction_g=direction_g, tau_direction=kinetic_direction, &
                             vxc_rho_response=response_rho, vxc_tau_response=response_tau)
      ELSE
         CALL semilocal_grid_response(xc, density, density_g, kinetic, direction, direction_g, kinetic_direction, &
                                      weights, pool, response_rho, response_tau)
      END IF
      IF (atomic_xc) THEN
         CALL qs_upf_inverse_scf_xc_model_response(ks_env, matrices_direction, direction, kinetic_direction, xc, &
                                                   atomic, combined_rho, combined_tau, density_tangent=tangent)
      ELSE
         DO spin = 1, ns
            combined_rho(spin)%array = 0.0_dp
            combined_tau(spin)%array = 0.0_dp
         END DO
      END IF
      DO spin = 1, ns
         combined_rho(spin)%array = combined_rho(spin)%array + response_rho(spin)%array
         combined_tau(spin)%array = combined_tau(spin)%array + response_tau(spin)%array
      END DO
      CALL qs_upf_inverse_scf_xc_source_adjoint(ks_env, combined_rho, combined_tau, operators_direction)
      IF (atomic_xc) operators_direction(:, :, :, :) = operators_direction + atomic
      ALLOCATE (hartree(n, n, nk))
      IF (state%gaussian_background) THEN
         CPASSERT(ASSOCIATED(state%coulomb))
         CALL qs_upf_inverse_coulomb_operator_response(qs_env, state, state%coulomb, &
                                                       SUM(matrices_direction, DIM=4), hartree)
         DO spin = 1, ns
            operators_direction(:, :, :, spin) = operators_direction(:, :, :, spin) + hartree
         END DO
      ELSE
         CALL native_hartree_response(ks_env, matrices_direction, atomic)
         operators_direction(:, :, :, :) = operators_direction + atomic
      END IF
      IF (ANY(state%paw%active)) THEN
         ALLOCATE (paw(n, n, nk, ns))
         functional => section_vals_get_subs_vals(xc, "XC_FUNCTIONAL")
         IF (ASSOCIATED(state%bloch%kpoints)) THEN
            ALLOCATE (kweights, SOURCE=state%bloch%weights)
            CALL qs_upf_inverse_paw_evaluate(state%paw, state%bloch%frames, matrices_direction, kweights, &
                                             functional, section_get_rval(xc, "DENSITY_CUTOFF"), &
                                             section_get_rval(xc, "GRADIENT_CUTOFF"), energies, paw, paw_potential, &
                                             semilocal_xc=.FALSE.)
         ELSE
            ALLOCATE (kweights(1), SOURCE=1.0_dp)
            CALL qs_upf_inverse_paw_evaluate(state%paw, [state%frame], matrices_direction, kweights, &
                                             functional, section_get_rval(xc, "DENSITY_CUTOFF"), &
                                             section_get_rval(xc, "GRADIENT_CUTOFF"), energies, paw, paw_potential, &
                                             semilocal_xc=.FALSE.)
         END IF
         operators_direction(:, :, :, :) = operators_direction + paw
         IF (.NOT. ALLOCATED(state%paw_reconstruction)) THEN
            CPASSERT(.NOT. model)
            CPASSERT(ALLOCATED(state%paw_density))
            CPASSERT(state%paw_density_evaluation == state%density_evaluations)
            IF (ASSOCIATED(state%bloch%kpoints)) THEN
               CALL qs_upf_inverse_paw_xc_kernel(state%paw, state%bloch%frames, state%paw_density, &
                                                 matrices_direction, kweights, xc, paw)
            ELSE
               CALL qs_upf_inverse_paw_xc_kernel(state%paw, [state%frame], state%paw_density, &
                                                 matrices_direction, kweights, xc, paw)
            END IF
            operators_direction(:, :, :, :) = operators_direction + paw
         END IF
      END IF
      DO spin = 1, ns
         IF (model) THEN
            CALL pool%give_back_pw(first_rho(spin))
            CALL pool%give_back_pw(first_tau(spin))
         END IF
         CALL pool%give_back_pw(response_rho(spin))
         CALL pool%give_back_pw(response_tau(spin))
         CALL density(spin)%release()
         CALL density_g(spin)%release()
         CALL kinetic(spin)%release()
         CALL direction(spin)%release()
         CALL direction_g(spin)%release()
         CALL kinetic_direction(spin)%release()
         CALL combined_rho(spin)%release()
         CALL combined_tau(spin)%release()
      END DO
      IF (do_hfx) operators_direction(:, :, :, :) = operators_direction + exchange
      IF (model) DEALLOCATE (first_rho, first_tau)
      DEALLOCATE (response_rho, response_tau, density, density_g, kinetic, &
                  direction, direction_g, kinetic_direction)
      IF (fine .AND. ASSOCIATED(native_weights)) CALL transferred_weights%release()
      CALL work_g%release()
      available = .TRUE.
   END SUBROUTINE qs_upf_inverse_scf_electronic_response

! **************************************************************************************************
!> \brief Fixed-geometry RI exchange response with owned density and operator workspaces.
!>        Reuse the original integrals and spin factors without changing the accepted SCF histories.
!>        Frozen-core one-body exchange has zero density response. Its orbital curvature remains
!>        in the caller's complete KS operator. Unprepared, MO, ACE and non-RI cases stay numerical.
!> \param qs_env Evaluated native environment and original HF components.
!> \param state Current USPP/PAW reconstruction and physical Bloch map.
!> \param direction Unweighted Hermitian AO density directions for every physical contribution.
!> \param operators New sum of the original HF components, without mesh weights.
!> \param available Whether every original HF component has a prepared RI-RHO response.
! **************************************************************************************************
   SUBROUTINE native_exchange_response(qs_env, state, direction, operators, available)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_upf_inverse_scf_type), INTENT(IN)          :: state
      COMPLEX(dp), INTENT(IN)                            :: direction(:, :, :, :)
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: operators(:, :, :, :)
      LOGICAL, INTENT(OUT)                               :: available

      COMPLEX(dp), ALLOCATABLE                           :: contribution(:, :, :, :), dense(:, :)
      INTEGER                                            :: component, image, ns, spin
      LOGICAL                                            :: bloch, imaginary
      REAL(dp)                                           :: unused_energy
      TYPE(dbcsr_p_type), POINTER                        :: imag_operator(:, :), imag_source(:, :), &
                                                            images(:, :), real_operator(:, :), &
                                                            real_source(:, :)
      TYPE(hfx_type), POINTER                            :: exchange(:, :)
      TYPE(qs_rho_type), POINTER                         :: rho

      available = .FALSE.
      CALL get_qs_env(qs_env, x_data=exchange, rho=rho)
      IF (.NOT. ASSOCIATED(exchange)) RETURN
      IF (SIZE(exchange, 1) == 0) RETURN
      bloch = ASSOCIATED(state%bloch%kpoints)
      DO component = 1, SIZE(exchange, 1)
         ASSOCIATE (hf => exchange(component, 1))
            IF (.NOT. hf%do_hfx_ri .OR. hf%use_ace) RETURN
            IF (.NOT. ASSOCIATED(hf%ri_data)) RETURN
            IF (hf%ri_data%flavor /= ri_pmat) RETURN
            IF (.NOT. ALLOCATED(hf%ri_data%rho_ao_t)) RETURN
            IF (.NOT. ALLOCATED(hf%ri_data%ks_t)) RETURN
            IF (.NOT. ALLOCATED(hf%ri_data%t_3c_int_ctr_2)) RETURN
         END ASSOCIATE
      END DO
      CALL qs_rho_get(rho, rho_ao_kp=images)
      ns = SIZE(direction, 4)
      ALLOCATE (operators, MOLD=direction)
      operators = CMPLX(0.0_dp, 0.0_dp, dp)
      ALLOCATE (real_operator(ns, SIZE(images, 2)))
      DO image = 1, SIZE(images, 2)
         DO spin = 1, ns
            ALLOCATE (real_operator(spin, image)%matrix)
            CALL dbcsr_create(real_operator(spin, image)%matrix, template=images(spin, image)%matrix)
         END DO
      END DO
      imaginary = ANY(AIMAG(direction) /= 0.0_dp)
      IF (.NOT. bloch) THEN
         CPASSERT(SIZE(direction, 3) == 1)
         ALLOCATE (real_source(ns, 1), imag_source(ns, 1), imag_operator(ns, 1))
         DO spin = 1, ns
            ALLOCATE (real_source(spin, 1)%matrix, imag_source(spin, 1)%matrix, imag_operator(spin, 1)%matrix)
            CALL qs_upf_inverse_matrix_scatter(state%frame%ao, direction(:, :, 1, spin), images(spin, 1)%matrix, &
                                               real_source(spin, 1)%matrix, imag_source(spin, 1)%matrix)
            CALL dbcsr_create(imag_operator(spin, 1)%matrix, template=imag_source(spin, 1)%matrix)
         END DO
      END IF
      DO component = 1, SIZE(exchange, 1)
         ASSOCIATE (hf => exchange(component, 1))
            IF (bloch) THEN
               CALL hfx_ri_update_ks_kp(qs_env, hf%ri_data, real_operator, unused_energy, images, .FALSE., ns, &
                                        hf%general_parameter%fraction, density_matrices=direction, &
                                        operators_response=contribution)
               operators(:, :, :, :) = operators + contribution
            ELSE
               CALL hfx_ri_update_ks(qs_env, hf%ri_data, real_operator, unused_energy, rho_ao=real_source, &
                                     geometry_did_change=.FALSE., nspins=ns, &
                                     hf_fraction=hf%general_parameter%fraction, response=.TRUE.)
               IF (imaginary) THEN
                  CALL hfx_ri_update_ks(qs_env, hf%ri_data, imag_operator, unused_energy, rho_ao=imag_source, &
                                        geometry_did_change=.FALSE., nspins=ns, &
                                        hf_fraction=hf%general_parameter%fraction, response=.TRUE.)
               END IF
            END IF
         END ASSOCIATE
      END DO
      IF (.NOT. bloch) THEN
         DO spin = 1, ns
            CALL qs_upf_inverse_matrix_gather(state%frame%ao, real_operator(spin, 1)%matrix, dense, &
                                              imag_operator(spin, 1)%matrix)
            operators(:, :, 1, spin) = dense
            CALL dbcsr_release(real_source(spin, 1)%matrix)
            CALL dbcsr_release(imag_source(spin, 1)%matrix)
            CALL dbcsr_release(imag_operator(spin, 1)%matrix)
            DEALLOCATE (real_source(spin, 1)%matrix, imag_source(spin, 1)%matrix, imag_operator(spin, 1)%matrix)
         END DO
         DEALLOCATE (real_source, imag_source, imag_operator)
      END IF
      DO image = 1, SIZE(real_operator, 2)
         DO spin = 1, ns
            CALL dbcsr_release(real_operator(spin, image)%matrix)
            DEALLOCATE (real_operator(spin, image)%matrix)
         END DO
      END DO
      DEALLOCATE (real_operator)
      available = .TRUE.
   END SUBROUTINE native_exchange_response

! **************************************************************************************************
!> \brief Native semilocal XC response using local points and the original Fourier gradient transpose.
!>        Independent density and kinetic-density directions use the original physical fields.
!>        Point covectors carry the original fixed XC integration weights before the spectral
!>        transpose. The returned unweighted scalar potentials enter the original AO source map.
!>        No global descriptor array or accepted native cache is modified.
!> \param xc Original complete XC section.
!> \param density Original spin densities, including the fixed core.
!> \param density_g Matching original reciprocal fields.
!> \param kinetic Original independent kinetic densities, zero when unused.
!> \param direction Spin-density directions, with fixed cores.
!> \param direction_g Matching reciprocal directions.
!> \param kinetic_direction Independent kinetic-density directions.
!> \param weights Original fixed XC integration weights, possibly disassociated.
!> \param pool Original XC-grid pool.
!> \param potential New pool-owned scalar density-potential responses.
!> \param kinetic_potential New pool-owned kinetic-potential responses.
! **************************************************************************************************
   SUBROUTINE semilocal_grid_response(xc, density, density_g, kinetic, direction, direction_g, kinetic_direction, &
                                      weights, pool, potential, kinetic_potential)
      TYPE(section_vals_type), POINTER                   :: xc
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: density(:)
      TYPE(pw_c1d_gs_type), INTENT(IN)                   :: density_g(:)
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: kinetic(:), direction(:)
      TYPE(pw_c1d_gs_type), INTENT(IN)                   :: direction_g(:)
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: kinetic_direction(:)
      TYPE(pw_r3d_rs_type), POINTER                      :: weights
      TYPE(pw_pool_type), POINTER                        :: pool
      TYPE(pw_r3d_rs_type), POINTER                      :: potential(:), kinetic_potential(:)

      INTEGER                                            :: axis, dims(3), np, ns, orders(3), spin
      REAL(dp), ALLOCATABLE                              :: dgrad(:, :, :), drho(:, :), dtau(:, :), &
                                                            grad(:, :, :), rho(:, :), tau(:, :), &
                                                            vg(:, :, :), vn(:, :), vt(:, :), w(:)
      TYPE(pw_c1d_gs_type)                               :: derivative_g
      TYPE(pw_grid_type), POINTER                        :: grid
      TYPE(pw_r3d_rs_type)                               :: derivative_r, gradient_r, work_r

      grid => density(1)%pw_grid
      np = SIZE(density(1)%array)
      ns = SIZE(density)
      dims = SHAPE(density(1)%array)
      ALLOCATE (rho(np, ns), grad(np, 3, ns), tau(np, ns), drho(np, ns), dgrad(np, 3, ns), dtau(np, ns), &
                vn(np, ns), vg(np, 3, ns), vt(np, ns), w(np))
      w = 1.0_dp
      IF (ASSOCIATED(weights)) THEN
         CPASSERT(pw_grid_compare(weights%pw_grid, grid))
         w(:) = RESHAPE(weights%array, [np])
      END IF
      CALL derivative_g%create(grid)
      CALL gradient_r%create(grid)
      CALL work_r%create(grid)
      CALL derivative_r%create(grid)
      DO spin = 1, ns
         CPASSERT(ASSOCIATED(density(spin)%pw_grid, grid))
         CPASSERT(ASSOCIATED(direction(spin)%pw_grid, grid))
         rho(:, spin) = RESHAPE(density(spin)%array, [np])
         tau(:, spin) = RESHAPE(kinetic(spin)%array, [np])
         drho(:, spin) = RESHAPE(direction(spin)%array, [np])
         dtau(:, spin) = RESHAPE(kinetic_direction(spin)%array, [np])
         DO axis = 1, 3
            orders = 0
            orders(axis) = 1
            derivative_g%array = density_g(spin)%array
            CALL pw_derive(derivative_g, orders)
            CALL pw_transfer(derivative_g, gradient_r)
            grad(:, axis, spin) = RESHAPE(gradient_r%array, [np])
            derivative_g%array = direction_g(spin)%array
            CALL pw_derive(derivative_g, orders)
            CALL pw_transfer(derivative_g, gradient_r)
            dgrad(:, axis, spin) = RESHAPE(gradient_r%array, [np])
         END DO
      END DO
      CALL qs_upf_inverse_xc_field_kernel(xc, rho, grad, tau, drho, dgrad, dtau, vn, vg, vt)
      ALLOCATE (potential(ns), kinetic_potential(ns))
      DO spin = 1, ns
         CALL pool%create_pw(potential(spin))
         CALL pool%create_pw(kinetic_potential(spin))
         potential(spin)%array = RESHAPE(w*vn(:, spin), dims)
         kinetic_potential(spin)%array = RESHAPE(w*vt(:, spin), dims)
         DO axis = 1, 3
            work_r%array = RESHAPE(w*vg(:, axis, spin), dims)
            CALL pw_transfer(work_r, derivative_g)
            orders = 0
            orders(axis) = 1
            CALL pw_derive(derivative_g, orders)
            CALL pw_transfer(derivative_g, derivative_r)
            potential(spin)%array = potential(spin)%array - derivative_r%array
         END DO
      END DO
      CALL derivative_g%release()
      CALL gradient_r%release()
      CALL work_r%release()
      CALL derivative_r%release()
   END SUBROUTINE semilocal_grid_response

! **************************************************************************************************
!> \brief Apply the native Fourier Hartree kernel and the same complete AO source transpose.
!> \param ks_env Native inverse environment.
!> \param matrices_direction Unweighted spin-resolved AO directions.
!> \param operators_direction Complete spin-independent Hartree response, replicated over spins.
! **************************************************************************************************
   SUBROUTINE native_hartree_response(ks_env, matrices_direction, operators_direction)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      COMPLEX(dp), INTENT(IN)                            :: matrices_direction(:, :, :, :)
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT) :: operators_direction(:, :, :, :)

      INTEGER                                            :: ns, spin
      TYPE(pw_c1d_gs_type)                               :: potential_g, source_g
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_poisson_type), POINTER                     :: poisson
      TYPE(pw_r3d_rs_type), ALLOCATABLE                  :: direction(:), kinetic_direction(:), &
                                                            potential(:), zero(:)
      TYPE(pw_r3d_rs_type), POINTER                      :: native_r(:)
      TYPE(qs_rho_type), POINTER                         :: rho

      CALL get_ks_env(ks_env, rho=rho, pw_env=pw_env)
      CALL qs_rho_get(rho, rho_r=native_r)
      CALL pw_env_get(pw_env, poisson_env=poisson)
      CPASSERT(ASSOCIATED(poisson%green_fft))
      ns = SIZE(native_r)
      ALLOCATE (direction(ns), kinetic_direction(ns), potential(ns), zero(ns))
      DO spin = 1, ns
         CALL direction(spin)%create(native_r(spin)%pw_grid)
         CALL kinetic_direction(spin)%create(native_r(spin)%pw_grid)
         CALL potential(spin)%create(native_r(spin)%pw_grid)
         CALL zero(spin)%create(native_r(spin)%pw_grid)
         zero(spin)%array = 0.0_dp
      END DO
      CALL qs_upf_inverse_scf_xc_source_direction(ks_env, matrices_direction, direction, kinetic_direction)
      DO spin = 2, ns
         direction(1)%array = direction(1)%array + direction(spin)%array
      END DO
      CALL source_g%create(native_r(1)%pw_grid)
      CALL potential_g%create(native_r(1)%pw_grid)
      CALL pw_transfer(direction(1), source_g)
      potential_g%array = poisson%green_fft%influence_fn%array*source_g%array
      CALL pw_transfer(potential_g, potential(1))
      DO spin = 2, ns
         potential(spin)%array = potential(1)%array
      END DO
      CALL qs_upf_inverse_scf_xc_source_adjoint(ks_env, potential, zero, operators_direction)
      DO spin = 1, ns
         CALL direction(spin)%release()
         CALL kinetic_direction(spin)%release()
         CALL potential(spin)%release()
         CALL zero(spin)%release()
      END DO
      CALL source_g%release()
      CALL potential_g%release()
   END SUBROUTINE native_hartree_response

! **************************************************************************************************
!> \brief Apply the original real/Fourier transfer to an owned real field.
!> \param source Native original field.
!> \param TARGET Created destination field.
! **************************************************************************************************
   SUBROUTINE transfer_real(source, TARGET)
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: source
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: target

      TYPE(pw_c1d_gs_type)                               :: source_g, target_g

      IF (pw_grid_compare(source%pw_grid, target%pw_grid)) THEN
         target%array = source%array
      ELSE
         CALL source_g%create(source%pw_grid)
         CALL target_g%create(target%pw_grid)
         CALL pw_transfer(source, source_g)
         CALL pw_real_transfer(source_g, target_g)
         CALL pw_transfer(target_g, TARGET)
         CALL source_g%release()
         CALL target_g%release()
      END IF
   END SUBROUTINE transfer_real

END MODULE qs_upf_inverse_scf_response
