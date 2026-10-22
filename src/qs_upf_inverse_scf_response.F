!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
!> \brief Fixed-geometry native electronic operator response for the original Skala functional.
MODULE qs_upf_inverse_scf_response
   USE cell_types,                      ONLY: cell_type
   USE cp_control_types,                ONLY: dft_control_type
   USE input_constants,                 ONLY: sic_none,&
                                              xc_vdw_fun_nonloc
   USE input_section_types,             ONLY: section_get_ival,&
                                              section_get_rval,&
                                              section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
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
   USE qs_upf_inverse_core,             ONLY: qs_upf_inverse_core_kinetic
   USE qs_upf_inverse_coulomb,          ONLY: qs_upf_inverse_coulomb_operator_response
   USE qs_upf_inverse_paw,              ONLY: qs_upf_inverse_paw_evaluate
   USE qs_upf_inverse_scf_xc,           ONLY: qs_upf_inverse_scf_xc_model_response,&
                                              qs_upf_inverse_scf_xc_source_adjoint,&
                                              qs_upf_inverse_scf_xc_source_direction
   USE qs_upf_inverse_types,            ONLY: qs_upf_inverse_scf_type
   USE skala_density_constraints,       ONLY: skala_density_constraints_count,&
                                              skala_density_constraints_type
   USE skala_gpw_functional,            ONLY: native_skala_gapw_atom_composite_requested,&
                                              native_skala_gapw_composite_reference,&
                                              skala_gpw_eval,&
                                              xc_section_uses_native_skala_evaluator
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: qs_upf_inverse_scf_electronic_response

CONTAINS

! **************************************************************************************************
!> \brief Apply the original electronic KS kernel to a joint spin/k occupation-matrix direction.
!>        Combine the complete global Skala model response, exact atomic-minus-background response,
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
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_scf_electronic_response(qs_env, matrices_direction, operators_direction, available)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      COMPLEX(dp), INTENT(IN)                            :: matrices_direction(:, :, :, :)
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT) :: operators_direction(:, :, :, :)
      LOGICAL, INTENT(OUT)                               :: available

      COMPLEX(dp), ALLOCATABLE                           :: atomic(:, :, :, :), hartree(:, :, :), &
                                                            paw(:, :, :, :)
      INTEGER                                            :: n, nk, ns, selected, spin
      LOGICAL                                            :: atomic_xc, do_hfx, fine, rho_valid, &
                                                            tau_g_valid, tau_valid, valid
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
      TYPE(skala_density_constraints_type), POINTER      :: constraints

      available = .FALSE.
      CALL get_qs_env(qs_env, ks_env=ks_env, rho=rho, dft_control=control, cell=cell, &
                      particle_set=particles, pw_env=pw_env, para_env=para_env)
      CALL get_ks_env(ks_env, upf_inverse=state, skala_constraints=constraints, &
                      rho_nlcc=core_r, rho_nlcc_g=native_core_g, xcint_weights=native_weights)
      IF (.NOT. ASSOCIATED(state)) RETURN
      IF (.NOT. state%enabled) RETURN
      xc => section_vals_get_subs_vals(qs_env%input, "DFT%XC")
      IF (.NOT. xc_section_uses_native_skala_evaluator(xc)) RETURN
      hf => section_vals_get_subs_vals(xc, "HF")
      CALL section_vals_get(hf, explicit=do_hfx)
      IF (do_hfx .OR. control%do_admm .OR. control%sic_method_id /= sic_none) RETURN
      IF (control%dft_plus_u .OR. control%do_sccs .OR. control%do_fde .OR. control%do_pcc) RETURN
      IF (control%apply_embed_pot .OR. control%apply_dmfet_pot .OR. control%apply_external_density .OR. &
          control%apply_external_vxc .OR. control%qs_control%cdft .OR. control%correct_surf_dip .OR. qs_env%qmmm) RETURN
      IF (section_get_ival(xc, "VDW_POTENTIAL%POTENTIAL_TYPE") == xc_vdw_fun_nonloc) RETURN
      IF (native_skala_gapw_composite_reference(xc) .OR. native_skala_gapw_atom_composite_requested(xc)) RETURN
      selected = 0
      IF (ASSOCIATED(constraints)) selected = skala_density_constraints_count(constraints)
      CALL para_env%max(selected)
      IF (selected /= 0) RETURN
      CALL qs_rho_get(rho, rho_r=native_r, rho_g=native_g, tau_r=native_tau, tau_g=native_tau_g, &
                      rho_g_valid=rho_valid, tau_r_valid=tau_valid, tau_g_valid=tau_g_valid)
      IF (.NOT. rho_valid .OR. .NOT. tau_valid) RETURN
      CPASSERT(ASSOCIATED(native_r) .AND. ASSOCIATED(native_g))
      CPASSERT(ASSOCIATED(native_tau) .AND. ASSOCIATED(native_tau_g))
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
         CALL direction(spin)%create(grid)
         CALL direction_g(spin)%create(grid)
         CALL kinetic_direction(spin)%create(grid)
         CALL combined_rho(spin)%create(grid)
         CALL combined_tau(spin)%create(grid)
         CALL pw_transfer(native_g(spin), density_g(spin))
         IF (fine) THEN
            CALL pw_transfer(density_g(spin), density(spin))
            CALL pw_transfer(density(spin), density_g(spin))
            IF (tau_g_valid) THEN
               CALL pw_transfer(native_tau_g(spin), work_g)
            ELSE
               BLOCK
                  TYPE(pw_c1d_gs_type) :: tau_aux
                  CALL tau_aux%create(native_tau(spin)%pw_grid)
                  CALL pw_transfer(native_tau(spin), tau_aux)
                  CALL pw_transfer(tau_aux, work_g)
                  CALL tau_aux%release()
               END BLOCK
            END IF
            CALL pw_transfer(work_g, kinetic(spin))
         ELSE
            density(spin)%array = native_r(spin)%array
            kinetic(spin)%array = native_tau(spin)%array
         END IF
      END DO
      IF (ASSOCIATED(core_r)) THEN
         CPASSERT(ASSOCIATED(native_core_g))
         CALL transferred_core%create(grid)
         CALL core_g%create(grid)
         CALL pw_transfer(native_core_g, core_g)
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
      CALL qs_upf_inverse_core_kinetic(ks_env, kinetic, 1.0_dp)
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
      CALL skala_gpw_eval(first_rho, first_tau, energy, density, density_g, kinetic, xc, weights, pool, &
                          particles, cell, .FALSE., virial, density_constraints=constraints, &
                          rho_direction=direction, rho_direction_g=direction_g, tau_direction=kinetic_direction, &
                          vxc_rho_response=response_rho, vxc_tau_response=response_tau)
      IF (atomic_xc) THEN
         CALL qs_upf_inverse_scf_xc_model_response(ks_env, matrices_direction, direction, kinetic_direction, xc, &
                                                   atomic, combined_rho, combined_tau)
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
      IF (atomic_xc) operators_direction = operators_direction + atomic
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
         operators_direction = operators_direction + atomic
      END IF
      IF (ANY(state%paw%active)) THEN
         ALLOCATE (paw(n, n, nk, ns))
         functional => section_vals_get_subs_vals(xc, "XC_FUNCTIONAL")
         IF (ASSOCIATED(state%bloch%kpoints)) THEN
            kweights = state%bloch%weights
            CALL qs_upf_inverse_paw_evaluate(state%paw, state%bloch%frames, matrices_direction, kweights, &
                                             functional, section_get_rval(xc, "DENSITY_CUTOFF"), &
                                             section_get_rval(xc, "GRADIENT_CUTOFF"), energies, paw, paw_potential, &
                                             semilocal_xc=.FALSE.)
         ELSE
            kweights = [1.0_dp]
            CALL qs_upf_inverse_paw_evaluate(state%paw, [state%frame], matrices_direction, kweights, &
                                             functional, section_get_rval(xc, "DENSITY_CUTOFF"), &
                                             section_get_rval(xc, "GRADIENT_CUTOFF"), energies, paw, paw_potential, &
                                             semilocal_xc=.FALSE.)
         END IF
         operators_direction = operators_direction + paw
      END IF
      DO spin = 1, ns
         CALL pool%give_back_pw(first_rho(spin))
         CALL pool%give_back_pw(first_tau(spin))
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
      DEALLOCATE (first_rho, first_tau, response_rho, response_tau, density, density_g, kinetic, &
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
         CALL pw_transfer(source_g, target_g)
         CALL pw_transfer(target_g, TARGET)
         CALL source_g%release()
         CALL target_g%release()
      END IF
   END SUBROUTINE transfer_real

END MODULE qs_upf_inverse_scf_response
