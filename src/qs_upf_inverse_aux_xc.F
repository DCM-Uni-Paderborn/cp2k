!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
!> \brief Auxiliary original-UPF XC energy and its complete physical AO derivative.
MODULE qs_upf_inverse_aux_xc
   USE accint_weights_forces,           ONLY: accint_weight_geometry
   USE admm_types,                      ONLY: admm_type
   USE cell_types,                      ONLY: cell_type
   USE cp_control_types,                ONLY: dft_control_type
   USE input_constants,                 ONLY: xc_none
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE particle_types,                  ONLY: particle_type
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: pw_r3d_rs_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type
   USE qs_rho_types,                    ONLY: qs_rho_get
   USE qs_upf_inverse_core,             ONLY: qs_upf_inverse_core_source_geometry
   USE qs_upf_inverse_paw,              ONLY: qs_upf_inverse_paw_evaluate,&
                                              qs_upf_inverse_paw_response
   USE qs_upf_inverse_scf_xc,           ONLY: qs_upf_inverse_scf_xc_response,&
                                              qs_upf_inverse_scf_xc_source_adjoint,&
                                              qs_upf_inverse_scf_xc_source_geometry
   USE qs_upf_inverse_types,            ONLY: qs_upf_inverse_scf_type
   USE qs_vxc,                          ONLY: qs_vxc_create
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: qs_upf_inverse_aux_xc_evaluate

CONTAINS

! **************************************************************************************************
!> \brief Evaluate the selected auxiliary XC functional on its own density and reconstruction.
!>        Retain the primary method's XC quadratures: the native grid reference and, where
!>        selected, the atom-centered original-UPF residual. The latter uses the auxiliary
!>        orbitals, original Q, prescribed cores and its complete grid and AO transposes.
!>        Native GAPW charge/OCE sources retain their separate charge and electrostatic role.
!>        PAW adds the selected AE-minus-PS XC term with original frozen cores, without Hartree
!>        or dataset reference energies. Operators retain original auxiliary AO coordinates.
!> \param qs_env Native auxiliary density and its independent original-UPF representation.
!> \param xc_section Actual selected auxiliary functional, including its exchange scaling.
!> \param grid_energy Native smooth-grid XC energy.
!> \param atom_energy Original reconstruction residual plus the selected PAW one-center XC.
!> \param operators Optional complete Hermitian AO derivative, without Brillouin-zone weights.
!> \param grid_virial Optional independently evaluated native grid XC virial.
!>        Atomic residual and frozen-core geometry contractions remain separate.
!> \param atom_response Optional total atomic motion, including affine cell motion.
!> \param cell_response Optional lattice derivative at fixed reduced k coordinates.
!> \param response Optional complete explicit XC energy derivative at fixed auxiliary AO density.
!>        Supply all three geometry arguments together. The result is globally reduced and
!>        includes moving sources, integration weights, XC coordinates and atomic corrections.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_aux_xc_evaluate(qs_env, xc_section, grid_energy, atom_energy, operators, grid_virial, &
                                             atom_response, cell_response, response)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: xc_section
      REAL(dp), INTENT(OUT)                              :: grid_energy, atom_energy
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT), OPTIONAL    :: operators(:, :, :, :)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: grid_virial(3, 3)
      REAL(dp), INTENT(IN), OPTIONAL                     :: atom_response(:, :), cell_response(3, 3)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: response

      INTEGER                                            :: functional_id, ns, spin
      LOGICAL                                            :: atomic, geometry, have_functional, valid
      REAL(dp)                                           :: density_cutoff, energies(3), &
                                                            gradient_cutoff, native_virial(3, 3)
      REAL(dp), ALLOCATABLE                              :: atom_gradient(:, :)
      TYPE(admm_type), POINTER                           :: admm
      TYPE(dft_control_type), POINTER                    :: control
      TYPE(particle_type), POINTER                       :: atoms(:)
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: pool
      TYPE(pw_r3d_rs_type)                               :: weight_potential
      TYPE(pw_r3d_rs_type), ALLOCATABLE                  :: kinetic_potential(:), potential(:)
      TYPE(pw_r3d_rs_type), POINTER                      :: density(:), v_rho(:), v_tau(:)
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_upf_inverse_scf_type), POINTER             :: state
      TYPE(section_vals_type), POINTER                   :: functional

      CALL get_qs_env(qs_env, ks_env=ks_env, admm_env=admm, dft_control=control)
      CPASSERT(ASSOCIATED(admm))
      CPASSERT(admm%upf_density%valid)
      state => admm%upf_inverse
      CPASSERT(ASSOCIATED(state))
      CPASSERT(state%basis_type == 'AUX_FIT')
      ns = SIZE(admm%upf_density%matrices, 4)
      geometry = PRESENT(response)
      CPASSERT(geometry .EQV. PRESENT(atom_response))
      CPASSERT(geometry .EQV. PRESENT(cell_response))
      functional => section_vals_get_subs_vals(xc_section, 'XC_FUNCTIONAL')
      CALL section_vals_val_get(xc_section, 'XC_FUNCTIONAL%_SECTION_PARAMETERS_', i_val=functional_id)
      have_functional = functional_id /= xc_none
      CALL get_ks_env(ks_env, pw_env=pw_env)
      CALL pw_env_get(pw_env, auxbas_pw_pool=pool)
      CALL qs_rho_get(admm%rho_aux_fit, rho_r=density)
      NULLIFY (v_rho, v_tau)
      IF (geometry) THEN
         CALL get_ks_env(ks_env, particle_set=atoms)
         valid = ALL(SHAPE(atom_response) == [3, SIZE(atoms)])
         CPASSERT(valid)
         ALLOCATE (atom_gradient(3, SIZE(atoms)), SOURCE=0.0_dp)
         CALL weight_potential%create(density(1)%pw_grid)
         CALL qs_vxc_create(ks_env, admm%rho_aux_fit, xc_section, v_rho, v_tau, grid_energy, &
                            just_energy=.FALSE., native_skala_atom_force=atom_gradient, upf_context=state, &
                            upf_density=admm%upf_density%matrices, upf_atomic_energy=atom_energy, &
                            upf_grid_virial=native_virial, upf_weight_potential=weight_potential)
         IF (PRESENT(grid_virial)) grid_virial = native_virial
      ELSE
         CALL qs_vxc_create(ks_env, admm%rho_aux_fit, xc_section, v_rho, v_tau, grid_energy, &
                            just_energy=.NOT. (PRESENT(operators) .OR. PRESENT(grid_virial)), upf_context=state, &
                            upf_density=admm%upf_density%matrices, upf_atomic_energy=atom_energy, upf_grid_virial=grid_virial)
      END IF

      IF (PRESENT(operators) .OR. geometry) THEN
         ALLOCATE (potential(ns), kinetic_potential(ns))
         DO spin = 1, ns
            CALL potential(spin)%create(density(spin)%pw_grid)
            CALL kinetic_potential(spin)%create(density(spin)%pw_grid)
            potential(spin)%array = 0.0_dp
            kinetic_potential(spin)%array = 0.0_dp
            IF (ASSOCIATED(v_rho)) THEN
               valid = ASSOCIATED(v_rho(spin)%pw_grid, density(spin)%pw_grid)
               CPASSERT(valid)
               potential(spin)%array = v_rho(spin)%array
            END IF
            IF (ASSOCIATED(v_tau)) THEN
               valid = ASSOCIATED(v_tau(spin)%pw_grid, density(spin)%pw_grid)
               CPASSERT(valid)
               kinetic_potential(spin)%array = v_tau(spin)%array
            END IF
         END DO
      END IF
      IF (PRESENT(operators)) THEN
         CALL qs_upf_inverse_scf_xc_source_adjoint(ks_env, potential, kinetic_potential, operators, context=state)
         ! GAPW_XC also uses the atom residual, although its charge remains a global grid source.
         atomic = control%qs_control%gapw .OR. control%qs_control%gapw_xc .OR. ALLOCATED(state%paw_reconstruction)
         IF (atomic .AND. have_functional) THEN
            CPASSERT(ALLOCATED(state%xc%operators))
            valid = ALL(SHAPE(state%xc%operators) == SHAPE(operators))
            CPASSERT(valid)
            operators(:, :, :, :) = operators + state%xc%operators
         END IF
      END IF

      IF (have_functional .AND. ANY(state%paw%active)) THEN
         CALL section_vals_val_get(xc_section, 'DENSITY_CUTOFF', r_val=density_cutoff)
         CALL section_vals_val_get(xc_section, 'GRADIENT_CUTOFF', r_val=gradient_cutoff)
         IF (ALLOCATED(state%paw_density)) DEALLOCATE (state%paw_density, state%paw_operator)
         ALLOCATE (state%paw_density, SOURCE=admm%upf_density%matrices)
         ALLOCATE (state%paw_operator, MOLD=admm%upf_density%matrices)
         IF (ASSOCIATED(state%bloch%kpoints)) THEN
            CALL qs_upf_inverse_paw_evaluate(state%paw, state%bloch%frames, state%paw_density, state%bloch%weights, &
                                             functional, density_cutoff, gradient_cutoff, energies, &
                                             state%paw_operator, state%paw_potential, &
                                             semilocal_xc=.NOT. ALLOCATED(state%paw_reconstruction), hartree=.FALSE.)
         ELSE
            CALL qs_upf_inverse_paw_evaluate(state%paw, [state%frame], state%paw_density, [1.0_dp], &
                                             functional, density_cutoff, gradient_cutoff, energies, &
                                             state%paw_operator, state%paw_potential, &
                                             semilocal_xc=.NOT. ALLOCATED(state%paw_reconstruction), hartree=.FALSE.)
         END IF
         atom_energy = atom_energy + energies(2)
         IF (PRESENT(operators)) operators(:, :, :, :) = operators + state%paw_operator
         state%paw_density_evaluation = state%density_evaluations
      END IF
      IF (geometry) THEN
         response = 0.0_dp
         IF (have_functional) THEN
            CALL xc_geometry(qs_env, state, potential, kinetic_potential, weight_potential, grid_energy, &
                             native_virial, atom_gradient, atom_response, cell_response, response)
         END IF
         CALL weight_potential%release()
      END IF
      IF (ALLOCATED(potential)) THEN
         DO spin = 1, ns
            CALL potential(spin)%release()
            CALL kinetic_potential(spin)%release()
         END DO
      END IF
      IF (ASSOCIATED(v_rho)) THEN
         DO spin = 1, SIZE(v_rho)
            CALL pool%give_back_pw(v_rho(spin))
         END DO
         DEALLOCATE (v_rho)
      END IF
      IF (ASSOCIATED(v_tau)) THEN
         DO spin = 1, SIZE(v_tau)
            CALL pool%give_back_pw(v_tau(spin))
         END DO
         DEALLOCATE (v_tau)
      END IF
   END SUBROUTINE qs_upf_inverse_aux_xc_evaluate

! **************************************************************************************************
!> \brief Differentiate the same auxiliary XC energy without changing accepted primary geometry data.
!>        AO occupations are fixed here. Projection and normalization responses are separate.
!> \param qs_env ...
!> \param state ...
!> \param potential ...
!> \param kinetic_potential ...
!> \param weight_potential Unweighted derivative with respect to the native GAPW grid weights.
!> \param grid_energy ...
!> \param grid_virial Local native XC tensor, excluding the ordinary volume measure.
!> \param atom_gradient Local explicit model energy gradient at fixed grid descriptors.
!> \param atom_response ...
!> \param cell_response ...
!> \param response Globally reduced positive energy directional derivative.
! **************************************************************************************************
   SUBROUTINE xc_geometry(qs_env, state, potential, kinetic_potential, weight_potential, grid_energy, grid_virial, &
                          atom_gradient, atom_response, cell_response, response)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_upf_inverse_scf_type), POINTER             :: state
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: potential(:), kinetic_potential(:), &
                                                            weight_potential
      REAL(dp), INTENT(IN)                               :: grid_energy, grid_virial(3, 3), &
                                                            atom_gradient(:, :), &
                                                            atom_response(:, :), &
                                                            cell_response(3, 3)
      REAL(dp), INTENT(OUT)                              :: response

      INTEGER                                            :: atom, nk
      REAL(dp)                                           :: partial, strain(3, 3), trace_strain
      REAL(dp), ALLOCATABLE                              :: reduced_k_response(:, :), &
                                                            relative_motion(:, :)
      TYPE(admm_type), POINTER                           :: admm
      TYPE(cell_type), POINTER                           :: cell
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), POINTER                       :: atoms(:)
      TYPE(qs_ks_env_type), POINTER                      :: ks_env

      CALL get_qs_env(qs_env, ks_env=ks_env, admm_env=admm, cell=cell, para_env=para_env, particle_set=atoms)
      strain = MATMUL(cell_response, cell%h_inv)
      trace_strain = strain(1, 1) + strain(2, 2) + strain(3, 3)
      ALLOCATE (relative_motion, SOURCE=atom_response)
      DO atom = 1, SIZE(atoms)
         relative_motion(:, atom) = atom_response(:, atom) - MATMUL(strain, atoms(atom)%r)
      END DO
      ! The native tensor already contains affine model coordinates, spectral-gradient
      ! response and the model's extra weight response relative to the ordinary measure.
      partial = SUM(grid_virial*strain) + SUM(atom_gradient*relative_motion)
      CALL para_env%sum(partial)
      response = trace_strain*grid_energy + partial
      CALL qs_upf_inverse_scf_xc_source_geometry(qs_env, admm%upf_density%matrices, potential, kinetic_potential, &
                                                 atom_response, cell_response, partial, context=state)
      response = response + partial
      CALL qs_upf_inverse_core_source_geometry(qs_env, potential, atom_response, cell_response, partial, &
                                               kinetic_potential=kinetic_potential, context=state)
      response = response + partial
      CALL accint_weight_geometry(qs_env, weight_potential, atom_response, cell_response, partial)
      response = response + partial
      nk = SIZE(admm%upf_density%matrices, 3)
      ALLOCATE (reduced_k_response(3, nk), SOURCE=0.0_dp)
      CALL qs_upf_inverse_scf_xc_response(ks_env, atom_response, cell_response, reduced_k_response, partial, context=state)
      response = response + partial
      IF (ANY(state%paw%active)) THEN
         IF (ASSOCIATED(state%bloch%kpoints)) THEN
            CALL qs_upf_inverse_paw_response(state%bloch%frames, state%paw_density, state%bloch%weights, &
                                             state%paw_potential, atom_response, cell_response, reduced_k_response, partial)
         ELSE
            CALL qs_upf_inverse_paw_response([state%frame], state%paw_density, [1.0_dp], state%paw_potential, &
                                             atom_response, cell_response, reduced_k_response, partial)
         END IF
         response = response + partial
      END IF
   END SUBROUTINE xc_geometry
END MODULE qs_upf_inverse_aux_xc
