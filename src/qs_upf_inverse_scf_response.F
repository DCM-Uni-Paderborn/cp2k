!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
!> \brief Fixed-geometry native electronic operator response for semilocal, Skala and RI hybrid functionals.
MODULE qs_upf_inverse_scf_response
   USE admm_types,                      ONLY: admm_type
   USE cell_types,                      ONLY: cell_type
   USE cp_control_types,                ONLY: dft_control_type
   USE input_constants,                 ONLY: sic_none,&
                                              xc_vdw_fun_nonloc
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
   USE pw_methods,                      ONLY: pw_transfer
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
   USE qs_upf_inverse_aux_operator,     ONLY: qs_upf_inverse_aux_electronic_response,&
                                              qs_upf_inverse_ri_response
   USE qs_upf_inverse_core,             ONLY: qs_upf_inverse_core_kinetic
   USE qs_upf_inverse_coulomb,          ONLY: qs_upf_inverse_coulomb_operator_response
   USE qs_upf_inverse_paw,              ONLY: qs_upf_inverse_paw_evaluate,&
                                              qs_upf_inverse_paw_xc_kernel
   USE qs_upf_inverse_scf_xc,           ONLY: qs_upf_inverse_scf_xc_model_response,&
                                              qs_upf_inverse_scf_xc_source_adjoint,&
                                              qs_upf_inverse_scf_xc_source_direction,&
                                              qs_upf_inverse_xc_transfer_real
   USE qs_upf_inverse_types,            ONLY: qs_upf_inverse_orbitals_type,&
                                              qs_upf_inverse_scf_type
   USE qs_upf_inverse_xc,               ONLY: qs_upf_inverse_xc_grid_kernel,&
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
   USE xc_input_constants,              ONLY: xc_deriv_pw
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
!> \param orbital_direction Optional certified coefficient and occupation responses of this density direction.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_scf_electronic_response(qs_env, matrices_direction, operators_direction, available, density_tangent, &
                                                     orbital_direction)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      COMPLEX(dp), INTENT(IN)                            :: matrices_direction(:, :, :, :)
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT) :: operators_direction(:, :, :, :)
      LOGICAL, INTENT(OUT)                               :: available
      LOGICAL, INTENT(IN), OPTIONAL                      :: density_tangent
      TYPE(qs_upf_inverse_orbitals_type), INTENT(IN), &
         OPTIONAL                                        :: orbital_direction

      COMPLEX(dp), ALLOCATABLE                           :: atomic(:, :, :, :), &
                                                            exchange(:, :, :, :), &
                                                            hartree(:, :, :), paw(:, :, :, :)
      INTEGER                                            :: n, nk, ns, selected, spin
      LOGICAL                                            :: atomic_xc, do_hfx, fine, model, &
                                                            rho_valid, tangent, tau_g_valid, &
                                                            tau_valid, use_tau, valid
      REAL(dp)                                           :: energies(3), energy, virial(3, 3)
      REAL(dp), ALLOCATABLE                              :: kweights(:)
      TYPE(admm_type), POINTER                           :: admm
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
      TYPE(qs_upf_inverse_orbitals_type), POINTER        :: orbital_factors
      TYPE(qs_upf_inverse_scf_type), POINTER             :: state
      TYPE(section_vals_type), POINTER                   :: functional, hf, xc
      TYPE(skala_density_constraints_type), POINTER      :: constraints, response_constraints
      TYPE(xc_rho_cflags_type)                           :: needs

      available = .FALSE.
      tangent = .FALSE.
      IF (PRESENT(density_tangent)) tangent = density_tangent
      CALL get_qs_env(qs_env, ks_env=ks_env, rho=rho, dft_control=control, cell=cell, &
                      particle_set=particles, pw_env=pw_env, para_env=para_env, admm_env=admm)
      CALL get_ks_env(ks_env, upf_inverse=state, skala_constraints=constraints, &
                      rho_nlcc=core_r, rho_nlcc_g=native_core_g, xcint_weights=native_weights)
      IF (.NOT. ASSOCIATED(state)) RETURN
      IF (.NOT. state%enabled) RETURN
      xc => section_vals_get_subs_vals(qs_env%input, "DFT%XC")
      IF (control%do_admm) THEN
         IF (.NOT. ASSOCIATED(admm)) RETURN
         xc => admm%xc_section_primary
      END IF
      model = xc_section_uses_native_skala_evaluator(xc)
      IF (.NOT. model) THEN
         IF (.NOT. qs_upf_inverse_xc_semilocal_kernel_available(xc, control%nspins == 2)) RETURN
         IF (section_get_ival(xc, 'XC_GRID%XC_DERIV') /= xc_deriv_pw) RETURN
      END IF
      functional => section_vals_get_subs_vals(xc, 'XC_FUNCTIONAL')
      needs = xc_functionals_get_needs(functional, lsd=control%nspins == 2, calc_potential=.TRUE.)
      use_tau = model .OR. needs%tau .OR. needs%tau_spin
      IF (use_tau) THEN
         IF (.NOT. ALLOCATED(state%kinetic_cache)) RETURN
      END IF
      hf => section_vals_get_subs_vals(xc, "HF")
      CALL section_vals_get(hf, explicit=do_hfx)
      IF (control%sic_method_id /= sic_none) RETURN
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
      IF (control%do_admm) THEN
         CALL qs_upf_inverse_aux_electronic_response(qs_env, matrices_direction, exchange, valid)
         IF (.NOT. valid) RETURN
         do_hfx = .TRUE.
      ELSE IF (do_hfx) THEN
         CALL qs_upf_inverse_ri_response(qs_env, state, rho, matrices_direction, exchange, valid)
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
         CALL qs_upf_inverse_xc_transfer_real(native_weights, transferred_weights)
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
         CALL qs_upf_inverse_xc_grid_kernel(xc, density, density_g, kinetic, direction, direction_g, kinetic_direction, &
                                            weights, pool, response_rho, response_tau)
      END IF
      IF (atomic_xc) THEN
         CALL qs_upf_inverse_scf_xc_model_response(ks_env, matrices_direction, direction, kinetic_direction, xc, &
                                                   atomic, combined_rho, combined_tau, density_tangent=tangent, &
                                                   orbital_direction=orbital_direction)
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
      IF (.NOT. state%gaussian_background .AND. .NOT. fine) THEN
         ! The XC and Hartree responses use the same native density map.
         ! Add their covectors before applying its linear AO transpose.
         work_g%array = direction_g(1)%array
         DO spin = 2, ns
            work_g%array = work_g%array + direction_g(spin)%array
         END DO
         work_g%array = poisson%green_fft%influence_fn%array*work_g%array
         CALL pw_transfer(work_g, direction(1))
         DO spin = 1, ns
            combined_rho(spin)%array = combined_rho(spin)%array + direction(1)%array
         END DO
      END IF
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
      ELSE IF (fine) THEN
         CALL native_hartree_response(ks_env, matrices_direction, atomic)
         operators_direction(:, :, :, :) = operators_direction + atomic
      END IF
      NULLIFY (orbital_factors)
      IF (state%orbital_factors%density_evaluation == state%density_evaluations) orbital_factors => state%orbital_factors
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
                                                 matrices_direction, kweights, xc, paw, orbital_factors, orbital_direction)
            ELSE
               CALL qs_upf_inverse_paw_xc_kernel(state%paw, [state%frame], state%paw_density, &
                                                 matrices_direction, kweights, xc, paw, orbital_factors, orbital_direction)
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

END MODULE qs_upf_inverse_scf_response
