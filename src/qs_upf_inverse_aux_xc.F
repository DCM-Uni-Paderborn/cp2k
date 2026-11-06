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
   USE input_constants,                 ONLY: xc_none,&
                                              xc_vdw_fun_nonloc
   USE input_section_types,             ONLY: section_get_ival,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE particle_types,                  ONLY: particle_type
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_grids,                        ONLY: pw_grid_compare
   USE pw_methods,                      ONLY: pw_integrate_function,&
                                              pw_transfer
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type
   USE qs_rho_types,                    ONLY: qs_rho_create,&
                                              qs_rho_get,&
                                              qs_rho_release,&
                                              qs_rho_set,&
                                              qs_rho_type
   USE qs_upf_inverse_core,             ONLY: qs_upf_inverse_core_kinetic,&
                                              qs_upf_inverse_core_source_geometry
   USE qs_upf_inverse_paw,              ONLY: qs_upf_inverse_paw_evaluate,&
                                              qs_upf_inverse_paw_response,&
                                              qs_upf_inverse_paw_xc_kernel
   USE qs_upf_inverse_scf_xc,           ONLY: qs_upf_inverse_scf_xc_model_response,&
                                              qs_upf_inverse_scf_xc_response,&
                                              qs_upf_inverse_scf_xc_source_adjoint,&
                                              qs_upf_inverse_scf_xc_source_direction,&
                                              qs_upf_inverse_scf_xc_source_geometry
   USE qs_upf_inverse_types,            ONLY: qs_upf_inverse_scf_type
   USE qs_upf_inverse_xc,               ONLY: qs_upf_inverse_xc_grid_kernel,&
                                              qs_upf_inverse_xc_semilocal_kernel_available
   USE qs_vxc,                          ONLY: qs_vxc_create
   USE rs_pw_interface,                 ONLY: pw_real_transfer
   USE xc_input_constants,              ONLY: xc_deriv_pw,&
                                              xc_rho_no_smooth
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: qs_upf_inverse_aux_xc_evaluate, qs_upf_inverse_aux_xc_rescale
   PUBLIC :: qs_upf_inverse_aux_xc_kernel, qs_upf_inverse_aux_xc_kernel_available

CONTAINS

! **************************************************************************************************
!> \brief Spin-resolved Merlot scaling of the selected original-UPF auxiliary exchange functional.
!>        ADMMP evaluates the unnormalized native density and multiplies exchange by eta**2.
!>        ADMMS evaluates the normalized native density and multiplies exchange by eta**(2/3).
!>        Keep the zero-valence exchange baseline unscaled and scale only each spin's
!>        valence-dependent increment. Original frozen-core fields are never scaled. The
!>        selected exchange must be spin separable, checked against the complete evaluation.
!>        This is the XC contribution only. HF and the normalization chain rule remain separate.
!> \param qs_env Native auxiliary density, spin charge ratios and original frozen-core context.
!> \param xc_section Selected signed auxiliary exchange correction.
!> \param grid_energy Scaled smooth-grid energy, summed over spins.
!> \param atom_energy Scaled original atomic exchange correction, summed over spins.
!> \param operators Derivative with respect to the actual native auxiliary density, without k weights.
!> \param eta_derivative Explicit derivative with respect to each spin's eta at fixed native density.
!> \param atom_response Optional total atomic motion, including affine cell motion.
!> \param cell_response Optional lattice derivative at fixed reduced k coordinates.
!> \param response Optional explicit geometry derivative at fixed native density and eta.
!> \param atom_gradient Optional complete Cartesian gradients using one XC evaluation per spin component.
!> \param strain_gradient Optional complete affine strain derivatives with the same prepared potentials.
!>        Do not combine these outputs with a single requested response direction.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_aux_xc_rescale(qs_env, xc_section, grid_energy, atom_energy, operators, eta_derivative, &
                                            atom_response, cell_response, response, atom_gradient, strain_gradient)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: xc_section
      REAL(dp), INTENT(OUT)                              :: grid_energy, atom_energy
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: operators(:, :, :, :)
      REAL(dp), INTENT(OUT)                              :: eta_derivative(:)
      REAL(dp), INTENT(IN), OPTIONAL                     :: atom_response(:, :), cell_response(3, 3)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: response, atom_gradient(:, :), &
                                                            strain_gradient(3, 3)

      COMPLEX(dp), ALLOCATABLE                           :: matrices(:, :, :, :), &
                                                            part_operator(:, :, :, :), &
                                                            restored_operator(:, :, :, :)
      INTEGER                                            :: ns, spin
      LOGICAL                                            :: geometry, valid
      REAL(dp) :: atom_part, atom_sum, atom_zero, eta, exponent, grid_part, grid_sum, grid_zero, &
         partial, response_zero, restored_atom, restored_grid, weight
      REAL(dp), ALLOCATABLE                              :: part_gradient(:, :), part_strain(:, :), &
                                                            zero_gradient(:, :), zero_strain(:, :)
      TYPE(admm_type), POINTER                           :: admm

      CALL get_qs_env(qs_env, admm_env=admm)
      CPASSERT(ASSOCIATED(admm))
      CPASSERT(admm%upf_density%valid)
      ns = SIZE(admm%upf_density%matrices, 4)
      CPASSERT(SIZE(eta_derivative) == ns)
      geometry = PRESENT(response)
      valid = .NOT. (geometry .AND. (PRESENT(atom_gradient) .OR. PRESENT(strain_gradient)))
      CPASSERT(valid)
      valid = geometry .EQV. PRESENT(atom_response)
      CPASSERT(valid)
      valid = geometry .EQV. PRESENT(cell_response)
      CPASSERT(valid)
      eta_derivative = 0.0_dp
      IF (.NOT. (admm%do_admmp .OR. admm%do_admms)) THEN
         CALL qs_upf_inverse_aux_xc_evaluate(qs_env, xc_section, grid_energy, atom_energy, operators, &
                                             atom_response=atom_response, cell_response=cell_response, response=response, &
                                             atom_gradient=atom_gradient, strain_gradient=strain_gradient)
         RETURN
      END IF
      IF (PRESENT(atom_gradient)) THEN
         ALLOCATE (zero_gradient, MOLD=atom_gradient)
         ALLOCATE (part_gradient, MOLD=atom_gradient)
      END IF
      IF (PRESENT(strain_gradient)) THEN
         ALLOCATE (zero_strain(3, 3), part_strain(3, 3))
      END IF
      CPASSERT(ns == 1 .OR. ns == 2)
      exponent = MERGE(2.0_dp, 2.0_dp/3.0_dp, admm%do_admmp)
      ALLOCATE (operators, MOLD=admm%upf_density%matrices)
      operators = CMPLX(0.0_dp, 0.0_dp, dp)
      grid_energy = 0.0_dp
      atom_energy = 0.0_dp
      IF (geometry) response = 0.0_dp
      grid_zero = 0.0_dp
      atom_zero = 0.0_dp
      response_zero = 0.0_dp
      ALLOCATE (matrices, MOLD=admm%upf_density%matrices)
      matrices = CMPLX(0.0_dp, 0.0_dp, dp)
      CALL xc_component(qs_env, xc_section, matrices, grid_zero, atom_zero, part_operator, &
                        atom_response, cell_response, response_zero, zero_gradient, zero_strain)
      grid_energy = grid_zero
      atom_energy = atom_zero
      IF (geometry) response = response_zero
      IF (PRESENT(atom_gradient)) atom_gradient(:, :) = zero_gradient
      IF (PRESENT(strain_gradient)) strain_gradient(:, :) = zero_strain
      grid_sum = grid_zero
      atom_sum = atom_zero
      DO spin = 1, ns
         matrices = CMPLX(0.0_dp, 0.0_dp, dp)
         matrices(:, :, :, spin) = admm%upf_density%matrices(:, :, :, spin)
         CALL xc_component(qs_env, xc_section, matrices, grid_part, atom_part, part_operator, &
                           atom_response, cell_response, partial, part_gradient, part_strain)
         grid_part = grid_part - grid_zero
         atom_part = atom_part - atom_zero
         grid_sum = grid_sum + grid_part
         atom_sum = atom_sum + atom_part
         eta = admm%gsi(spin)
         CPASSERT(eta > 0.0_dp)
         weight = eta**exponent
         grid_energy = grid_energy + weight*grid_part
         atom_energy = atom_energy + weight*atom_part
         operators(:, :, :, spin) = weight*part_operator(:, :, :, spin)
         eta_derivative(spin) = exponent*weight*(grid_part + atom_part)/eta
         IF (geometry) response = response + weight*(partial - response_zero)
         IF (PRESENT(atom_gradient)) atom_gradient(:, :) = atom_gradient + weight*(part_gradient - zero_gradient)
         IF (PRESENT(strain_gradient)) strain_gradient(:, :) = strain_gradient + weight*(part_strain - zero_strain)
      END DO
      ! Restore native-density XC scratch data and verify the spin decomposition.
      CALL qs_upf_inverse_aux_xc_evaluate(qs_env, xc_section, restored_grid, restored_atom, restored_operator)
      valid = ABS(grid_sum - restored_grid) < 1.0E-9_dp*MAX(1.0_dp, ABS(restored_grid)) .AND. &
              ABS(atom_sum - restored_atom) < 1.0E-9_dp*MAX(1.0_dp, ABS(restored_atom))
      IF (.NOT. valid) CPABORT('Auxiliary exchange is not spin separable')
   END SUBROUTINE qs_upf_inverse_aux_xc_rescale

! **************************************************************************************************
!> \brief Evaluate one independent XC component, including geometry only when requested.
!> \param qs_env ...
!> \param xc_section ...
!> \param matrices ...
!> \param grid_energy ...
!> \param atom_energy ...
!> \param operators ...
!> \param atom_response ...
!> \param cell_response ...
!> \param response Zero without a requested geometry direction.
!> \param atom_gradient ...
!> \param strain_gradient ...
! **************************************************************************************************
   SUBROUTINE xc_component(qs_env, xc_section, matrices, grid_energy, atom_energy, operators, &
                           atom_response, cell_response, response, atom_gradient, strain_gradient)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: xc_section
      COMPLEX(dp), INTENT(IN)                            :: matrices(:, :, :, :)
      REAL(dp), INTENT(OUT)                              :: grid_energy, atom_energy
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: operators(:, :, :, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: atom_response(:, :), cell_response(3, 3)
      REAL(dp), INTENT(OUT)                              :: response
      REAL(dp), INTENT(OUT), OPTIONAL                    :: atom_gradient(:, :), &
                                                            strain_gradient(3, 3)

      response = 0.0_dp
      IF (PRESENT(atom_response)) THEN
         CALL qs_upf_inverse_aux_xc_evaluate(qs_env, xc_section, grid_energy, atom_energy, operators, &
                                             atom_response=atom_response, cell_response=cell_response, response=response, &
                                             density_matrices=matrices)
      ELSE
         CALL qs_upf_inverse_aux_xc_evaluate(qs_env, xc_section, grid_energy, atom_energy, operators, &
                                             density_matrices=matrices, atom_gradient=atom_gradient, &
                                             strain_gradient=strain_gradient)
      END IF
   END SUBROUTINE xc_component

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
!> \param density_matrices Optional original-coordinate auxiliary density for a separate XC evaluation.
!>        Build owned grid and kinetic fields without changing the native auxiliary density,
!>        MOs or RI history. Atomic XC scratch data then refer to these supplied matrices.
!> \param atom_gradient Optional complete Cartesian derivatives at fixed auxiliary occupations.
!> \param strain_gradient Optional complete affine derivatives with atoms moving with the cell.
!>        The energy, native model potential and atomic XC fields are prepared once.
!>        Each geometry direction retains its original analytic source and measure derivatives.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_aux_xc_evaluate(qs_env, xc_section, grid_energy, atom_energy, operators, grid_virial, &
                                             atom_response, cell_response, response, density_matrices, &
                                             atom_gradient, strain_gradient)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: xc_section
      REAL(dp), INTENT(OUT)                              :: grid_energy, atom_energy
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT), OPTIONAL    :: operators(:, :, :, :)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: grid_virial(3, 3)
      REAL(dp), INTENT(IN), OPTIONAL                     :: atom_response(:, :), cell_response(3, 3)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: response
      COMPLEX(dp), INTENT(IN), OPTIONAL, TARGET          :: density_matrices(:, :, :, :)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: atom_gradient(:, :), &
                                                            strain_gradient(3, 3)

      COMPLEX(dp), POINTER                               :: matrices(:, :, :, :)
      INTEGER                                            :: atom, axis, functional_id, ns, other, &
                                                            spin
      LOGICAL                                            :: atomic, geometry, have_functional, valid
      REAL(dp)                                           :: density_cutoff, energies(3), &
                                                            gradient_cutoff, native_virial(3, 3), &
                                                            strain(3, 3)
      REAL(dp), ALLOCATABLE                              :: model_gradient(:, :), motion(:, :)
      TYPE(admm_type), POINTER                           :: admm
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: control
      TYPE(particle_type), POINTER                       :: atoms(:)
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: pool
      TYPE(pw_r3d_rs_type)                               :: weight_potential
      TYPE(pw_r3d_rs_type), ALLOCATABLE                  :: kinetic_potential(:), potential(:)
      TYPE(pw_r3d_rs_type), POINTER                      :: density(:), v_rho(:), v_tau(:)
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_rho_type), TARGET                          :: temporary_rho
      TYPE(qs_upf_inverse_scf_type), POINTER             :: state
      TYPE(section_vals_type), POINTER                   :: functional

      CALL get_qs_env(qs_env, ks_env=ks_env, admm_env=admm, dft_control=control)
      CPASSERT(ASSOCIATED(admm))
      CPASSERT(admm%upf_density%valid)
      state => admm%upf_inverse
      CPASSERT(ASSOCIATED(state))
      CPASSERT(state%basis_type == 'AUX_FIT')
      matrices => admm%upf_density%matrices
      rho => admm%rho_aux_fit
      ns = SIZE(admm%upf_density%matrices, 4)
      geometry = PRESENT(response)
      CPASSERT(geometry .EQV. PRESENT(atom_response))
      CPASSERT(geometry .EQV. PRESENT(cell_response))
      valid = .NOT. (geometry .AND. (PRESENT(atom_gradient) .OR. PRESENT(strain_gradient)))
      CPASSERT(valid)
      geometry = geometry .OR. PRESENT(atom_gradient) .OR. PRESENT(strain_gradient)
      functional => section_vals_get_subs_vals(xc_section, 'XC_FUNCTIONAL')
      CALL section_vals_val_get(xc_section, 'XC_FUNCTIONAL%_SECTION_PARAMETERS_', i_val=functional_id)
      have_functional = functional_id /= xc_none
      CALL get_ks_env(ks_env, pw_env=pw_env)
      CALL pw_env_get(pw_env, auxbas_pw_pool=pool)
      CALL qs_rho_get(rho, rho_r=density)
      IF (PRESENT(density_matrices)) THEN
         valid = ALL(SHAPE(density_matrices) == SHAPE(matrices))
         CPASSERT(valid)
         matrices => density_matrices
         CALL build_xc_rho(ks_env, state, matrices, density, temporary_rho)
         rho => temporary_rho
         CALL qs_rho_get(rho, rho_r=density)
         ! A supplied density cannot reuse orbital factors of the native density.
         state%orbital_factors%density_evaluation = -1
      END IF
      NULLIFY (v_rho, v_tau)
      IF (geometry) THEN
         CALL get_ks_env(ks_env, particle_set=atoms)
         IF (PRESENT(response)) THEN
            valid = ALL(SHAPE(atom_response) == [3, SIZE(atoms)])
            CPASSERT(valid)
         END IF
         IF (PRESENT(atom_gradient)) THEN
            valid = ALL(SHAPE(atom_gradient) == [3, SIZE(atoms)])
            CPASSERT(valid)
         END IF
         ALLOCATE (model_gradient(3, SIZE(atoms)), SOURCE=0.0_dp)
         CALL weight_potential%create(density(1)%pw_grid)
         CALL qs_vxc_create(ks_env, rho, xc_section, v_rho, v_tau, grid_energy, &
                            just_energy=.FALSE., native_skala_atom_force=model_gradient, upf_context=state, &
                            upf_density=matrices, upf_atomic_energy=atom_energy, &
                            upf_grid_virial=native_virial, upf_weight_potential=weight_potential)
         IF (PRESENT(grid_virial)) grid_virial = native_virial
      ELSE
         CALL qs_vxc_create(ks_env, rho, xc_section, v_rho, v_tau, grid_energy, &
                            just_energy=.NOT. (PRESENT(operators) .OR. PRESENT(grid_virial)), upf_context=state, &
                            upf_density=matrices, upf_atomic_energy=atom_energy, upf_grid_virial=grid_virial)
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
         ALLOCATE (state%paw_density, SOURCE=matrices)
         ALLOCATE (state%paw_operator, MOLD=matrices)
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
         IF (PRESENT(response)) THEN
            response = 0.0_dp
            IF (have_functional) THEN
               CALL xc_geometry(qs_env, state, matrices, potential, kinetic_potential, weight_potential, grid_energy, &
                                native_virial, model_gradient, atom_response, cell_response, response)
            END IF
         END IF
         IF (PRESENT(atom_gradient) .OR. PRESENT(strain_gradient)) THEN
            CALL get_qs_env(qs_env, cell=cell)
            ALLOCATE (motion(3, SIZE(atoms)))
         END IF
         IF (PRESENT(atom_gradient)) THEN
            atom_gradient = 0.0_dp
            IF (have_functional) THEN
               DO atom = 1, SIZE(atoms)
                  DO axis = 1, 3
                     motion = 0.0_dp
                     motion(axis, atom) = 1.0_dp
                     CALL xc_geometry(qs_env, state, matrices, potential, kinetic_potential, weight_potential, grid_energy, &
                                      native_virial, model_gradient, motion, 0.0_dp*cell%hmat, atom_gradient(axis, atom))
                  END DO
               END DO
            END IF
         END IF
         IF (PRESENT(strain_gradient)) THEN
            strain_gradient = 0.0_dp
            IF (have_functional) THEN
               DO other = 1, 3
                  DO axis = 1, 3
                     strain = 0.0_dp
                     strain(axis, other) = 1.0_dp
                     DO atom = 1, SIZE(atoms)
                        motion(:, atom) = MATMUL(strain, atoms(atom)%r)
                     END DO
                     CALL xc_geometry(qs_env, state, matrices, potential, kinetic_potential, weight_potential, grid_energy, &
                                      native_virial, model_gradient, motion, MATMUL(strain, cell%hmat), &
                                      strain_gradient(axis, other))
                  END DO
               END DO
            END IF
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
      IF (PRESENT(density_matrices)) THEN
         CALL qs_rho_release(temporary_rho)
         state%paw_density_evaluation = -1
      END IF
   END SUBROUTINE qs_upf_inverse_aux_xc_evaluate

! **************************************************************************************************
!> \brief Build the same native XC source map for independent auxiliary occupations.
!>        The grid is borrowed from the native owner. All field arrays belong to rho and
!>        are released with it. Original frozen-core fields remain in the shared context.
!> \param ks_env ...
!> \param state ...
!> \param matrices ...
!> \param reference Native auxiliary grid owners, which are not modified.
!> \param rho New owned XC fields and their actual valence grid integrals.
! **************************************************************************************************
   SUBROUTINE build_xc_rho(ks_env, state, matrices, reference, rho)
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_upf_inverse_scf_type), INTENT(INOUT), &
         TARGET                                          :: state
      COMPLEX(dp), INTENT(IN)                            :: matrices(:, :, :, :)
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: reference(:)
      TYPE(qs_rho_type), INTENT(OUT)                     :: rho

      INTEGER                                            :: ns, spin
      REAL(dp), POINTER                                  :: charge(:)
      TYPE(pw_c1d_gs_type), POINTER                      :: density_g(:), kinetic_g(:)
      TYPE(pw_r3d_rs_type), POINTER                      :: density_r(:), kinetic_r(:)

      ns = SIZE(matrices, 4)
      CPASSERT(SIZE(reference) == ns)
      CALL qs_rho_create(rho)
      ALLOCATE (density_r(ns), density_g(ns), kinetic_r(ns), kinetic_g(ns), charge(ns))
      DO spin = 1, ns
         CALL density_r(spin)%create(reference(spin)%pw_grid)
         CALL density_g(spin)%create(reference(spin)%pw_grid)
         CALL kinetic_r(spin)%create(reference(spin)%pw_grid)
         CALL kinetic_g(spin)%create(reference(spin)%pw_grid)
      END DO
      CALL qs_upf_inverse_scf_xc_source_direction(ks_env, matrices, density_r, kinetic_r, context=state)
      DO spin = 1, ns
         CALL pw_transfer(density_r(spin), density_g(spin))
         CALL pw_transfer(kinetic_r(spin), kinetic_g(spin))
         charge(spin) = pw_integrate_function(density_r(spin))
      END DO
      CALL qs_rho_set(rho, rho_r=density_r, rho_g=density_g, tau_r=kinetic_r, tau_g=kinetic_g, &
                      rho_r_valid=.TRUE., rho_g_valid=.TRUE., &
                      tau_r_valid=ALLOCATED(state%kinetic_cache), tau_g_valid=ALLOCATED(state%kinetic_cache), &
                      tot_rho_r=charge)
   END SUBROUTINE build_xc_rho

! **************************************************************************************************
!> \brief Differentiate the same auxiliary XC energy without changing accepted primary geometry data.
!>        AO occupations are fixed here. Projection and normalization responses are separate.
!> \param qs_env ...
!> \param state ...
!> \param matrices Original-coordinate auxiliary occupations defining this XC evaluation.
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
   SUBROUTINE xc_geometry(qs_env, state, matrices, potential, kinetic_potential, weight_potential, grid_energy, grid_virial, &
                          atom_gradient, atom_response, cell_response, response)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_upf_inverse_scf_type), POINTER             :: state
      COMPLEX(dp), INTENT(IN)                            :: matrices(:, :, :, :)
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
      TYPE(cell_type), POINTER                           :: cell
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), POINTER                       :: atoms(:)
      TYPE(qs_ks_env_type), POINTER                      :: ks_env

      CALL get_qs_env(qs_env, ks_env=ks_env, cell=cell, para_env=para_env, particle_set=atoms)
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
      CALL qs_upf_inverse_scf_xc_source_geometry(qs_env, matrices, potential, kinetic_potential, &
                                                 atom_response, cell_response, partial, context=state)
      response = response + partial
      CALL qs_upf_inverse_core_source_geometry(qs_env, potential, atom_response, cell_response, partial, &
                                               kinetic_potential=kinetic_potential, context=state)
      response = response + partial
      CALL accint_weight_geometry(qs_env, weight_potential, atom_response, cell_response, partial)
      response = response + partial
      nk = SIZE(matrices, 3)
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
! **************************************************************************************************
!> \brief Whether the auxiliary functional has a complete fixed-geometry semilocal kernel.
!> \param xc Selected auxiliary XC section.
!> \param ns Number of spin components.
!> \return Whether the selected semilocal functional and grid support the kernel.
! **************************************************************************************************
   FUNCTION qs_upf_inverse_aux_xc_kernel_available(xc, ns) RESULT(available)
      TYPE(section_vals_type), POINTER                   :: xc
      INTEGER, INTENT(IN)                                :: ns
      LOGICAL                                            :: available

      available = qs_upf_inverse_xc_semilocal_kernel_available(xc, ns == 2)
      available = available .AND. section_get_ival(xc, 'XC_GRID%XC_DERIV') == xc_deriv_pw .AND. &
                  section_get_ival(xc, 'XC_GRID%XC_SMOOTH_RHO') == xc_rho_no_smooth .AND. &
                  section_get_ival(xc, 'VDW_POTENTIAL%POTENTIAL_TYPE') /= xc_vdw_fun_nonloc
   END FUNCTION qs_upf_inverse_aux_xc_kernel_available

! **************************************************************************************************
!> \brief Differentiate the same auxiliary XC functional at its actual density and charge ratio.
!>        Retain the unscaled zero-valence core baseline for Merlot exchange. Density and
!>        explicit-ratio responses are total electronic variations. Primary XC is separate.
!> \param qs_env Native auxiliary density, reconstruction and charge ratios.
!> \param xc Same auxiliary XC section as the accepted energy.
!> \param direction Variation of the actual auxiliary AO density.
!> \param eta_direction Variation of each primary-to-auxiliary charge ratio.
!> \param operators Auxiliary XC derivative at fixed charge ratios.
!> \param operator_direction Total electronic variation of that XC derivative.
!> \param explicit_derivative Partial derivative of auxiliary XC with respect to each charge ratio.
!> \param explicit_direction Total electronic variation of those partial derivatives.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_aux_xc_kernel(qs_env, xc, direction, eta_direction, operators, operator_direction, &
                                           explicit_derivative, explicit_direction)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: xc
      COMPLEX(dp), INTENT(IN)                            :: direction(:, :, :, :)
      REAL(dp), INTENT(IN)                               :: eta_direction(:)
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: operators(:, :, :, :), &
                                                            operator_direction(:, :, :, :)
      REAL(dp), INTENT(OUT)                              :: explicit_derivative(:), &
                                                            explicit_direction(:)

      COMPLEX(dp), ALLOCATABLE                           :: matrices(:, :, :, :), &
                                                            part_direction(:, :, :, :), &
                                                            part_operator(:, :, :, :), &
                                                            part_response(:, :, :, :)
      INTEGER                                            :: ik, ns, saved_factors, spin
      LOGICAL                                            :: dimensions_valid
      REAL(dp)                                           :: atom_energy, energy, eta, exponent, &
                                                            grid_energy, pairing, weight, zero
      TYPE(admm_type), POINTER                           :: admm
      TYPE(qs_upf_inverse_scf_type)                      :: saved

      CALL get_qs_env(qs_env, admm_env=admm)
      CPASSERT(ASSOCIATED(admm))
      ns = SIZE(direction, 4)
      CPASSERT(ALL(SHAPE(direction) == SHAPE(admm%upf_density%matrices)))
      dimensions_valid = SIZE(eta_direction) == ns .AND. SIZE(explicit_derivative) == ns .AND. &
                         SIZE(explicit_direction) == ns
      CPASSERT(dimensions_valid)
      CPASSERT(qs_upf_inverse_aux_xc_kernel_available(xc, ns))
      saved_factors = admm%upf_inverse%orbital_factors%density_evaluation
      ! Keep accepted atomic XC fields exactly, including their compensated orbital gauge.
      ! Independent component evaluations own the temporary density-dependent scratch.
      ASSOCIATE (state => admm%upf_inverse)
         saved%xc = state%xc
         saved%xc_grid => state%xc_grid
         saved%paw_density_evaluation = state%paw_density_evaluation
         CALL MOVE_ALLOC(state%xc_matrices, saved%xc_matrices)
         CALL MOVE_ALLOC(state%xc_fields, saved%xc_fields)
         CALL MOVE_ALLOC(state%xc_kweights, saved%xc_kweights)
         CALL MOVE_ALLOC(state%paw_density, saved%paw_density)
         CALL MOVE_ALLOC(state%paw_operator, saved%paw_operator)
         CALL MOVE_ALLOC(state%paw_potential, saved%paw_potential)
      END ASSOCIATE
      explicit_derivative = 0.0_dp
      explicit_direction = 0.0_dp
      IF (.NOT. (admm%do_admmp .OR. admm%do_admms)) THEN
         CALL xc_component_kernel(qs_env, xc, admm%upf_density%matrices, direction, energy, operators, operator_direction)
      ELSE
         exponent = MERGE(2.0_dp, 2.0_dp/3.0_dp, admm%do_admmp)
         ALLOCATE (operators, MOLD=direction)
         ALLOCATE (operator_direction, MOLD=direction)
         operators = CMPLX(0.0_dp, 0.0_dp, dp)
         operator_direction = CMPLX(0.0_dp, 0.0_dp, dp)
         ALLOCATE (matrices, MOLD=direction)
         ALLOCATE (part_direction, MOLD=direction)
         matrices = CMPLX(0.0_dp, 0.0_dp, dp)
         CALL qs_upf_inverse_aux_xc_evaluate(qs_env, xc, grid_energy, atom_energy, density_matrices=matrices)
         zero = grid_energy + atom_energy
         DO spin = 1, ns
            matrices = CMPLX(0.0_dp, 0.0_dp, dp)
            part_direction = CMPLX(0.0_dp, 0.0_dp, dp)
            matrices(:, :, :, spin) = admm%upf_density%matrices(:, :, :, spin)
            part_direction(:, :, :, spin) = direction(:, :, :, spin)
            CALL xc_component_kernel(qs_env, xc, matrices, part_direction, energy, part_operator, part_response)
            eta = admm%gsi(spin)
            CPASSERT(eta > 0.0_dp)
            weight = eta**exponent
            operators(:, :, :, spin) = weight*part_operator(:, :, :, spin)
            operator_direction(:, :, :, spin) = weight*part_response(:, :, :, spin) + &
                                                exponent*weight*eta_direction(spin)/eta*part_operator(:, :, :, spin)
            pairing = 0.0_dp
            DO ik = 1, SIZE(direction, 3)
               pairing = pairing + admm%upf_aux_metric%weights(ik)* &
                         REAL(SUM(CONJG(part_operator(:, :, ik, spin))*direction(:, :, ik, spin)), dp)
            END DO
            explicit_derivative(spin) = exponent*weight/eta*(energy - zero)
            explicit_direction(spin) = exponent*(exponent - 1.0_dp)*weight/eta**2*eta_direction(spin)*(energy - zero) + &
                                       exponent*weight/eta*pairing
         END DO
      END IF
      ASSOCIATE (state => admm%upf_inverse)
         state%xc = saved%xc
         state%xc_grid => saved%xc_grid
         state%paw_density_evaluation = saved%paw_density_evaluation
         CALL MOVE_ALLOC(saved%xc_matrices, state%xc_matrices)
         CALL MOVE_ALLOC(saved%xc_fields, state%xc_fields)
         CALL MOVE_ALLOC(saved%xc_kweights, state%xc_kweights)
         CALL MOVE_ALLOC(saved%paw_density, state%paw_density)
         CALL MOVE_ALLOC(saved%paw_operator, state%paw_operator)
         CALL MOVE_ALLOC(saved%paw_potential, state%paw_potential)
         state%orbital_factors%density_evaluation = saved_factors
      END ASSOCIATE
   END SUBROUTINE qs_upf_inverse_aux_xc_kernel

! **************************************************************************************************
!> \brief Auxiliary XC kernel on owned global fields and the matching original atomic residual.
!>        Frozen cores enter only the base fields. There is no auxiliary Hartree functional.
!> \param qs_env Native geometry and independent auxiliary reconstruction.
!> \param xc Selected auxiliary XC section.
!> \param matrices Actual auxiliary AO density defining the base fields.
!> \param matrices_direction Variation of that auxiliary AO density.
!> \param energy Same grid and atomic XC energy as the native evaluation.
!> \param operators Complete auxiliary XC first derivative.
!> \param operators_direction Complete fixed-geometry XC operator response.
! **************************************************************************************************
   SUBROUTINE xc_component_kernel(qs_env, xc, matrices, matrices_direction, energy, operators, operators_direction)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: xc
      COMPLEX(dp), INTENT(IN)                            :: matrices(:, :, :, :), &
                                                            matrices_direction(:, :, :, :)
      REAL(dp), INTENT(OUT)                              :: energy
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT) :: operators(:, :, :, :), &
         operators_direction(:, :, :, :)

      COMPLEX(dp), ALLOCATABLE                           :: atomic(:, :, :, :), paw(:, :, :, :)
      INTEGER                                            :: functional_id, ns, spin
      LOGICAL                                            :: atomic_xc, fine, use_tau
      REAL(dp)                                           :: atom_energy, grid_energy
      REAL(dp), ALLOCATABLE                              :: kweights(:)
      TYPE(admm_type), POINTER                           :: admm
      TYPE(dft_control_type), POINTER                    :: control
      TYPE(pw_c1d_gs_type)                               :: core_g, weight_source_g, weight_target_g
      TYPE(pw_c1d_gs_type), ALLOCATABLE                  :: density_g(:), direction_g(:)
      TYPE(pw_c1d_gs_type), POINTER                      :: native_core_g
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: native_pool, pool
      TYPE(pw_r3d_rs_type)                               :: transferred_core
      TYPE(pw_r3d_rs_type), ALLOCATABLE                  :: combined(:), combined_tau(:), &
                                                            density(:), direction(:), kinetic(:), &
                                                            kinetic_direction(:)
      TYPE(pw_r3d_rs_type), POINTER                      :: core, native_weights, response(:), &
                                                            response_tau(:), weights
      TYPE(pw_r3d_rs_type), TARGET                       :: transferred_weights
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_upf_inverse_scf_type), POINTER             :: state

      CALL get_qs_env(qs_env, ks_env=ks_env, admm_env=admm, dft_control=control)
      state => admm%upf_inverse
      CALL get_ks_env(ks_env, pw_env=pw_env, rho_nlcc=core, rho_nlcc_g=native_core_g, xcint_weights=native_weights)
      CALL pw_env_get(pw_env, xc_pw_pool=pool, auxbas_pw_pool=native_pool)
      CALL section_vals_val_get(xc, 'XC_FUNCTIONAL%_SECTION_PARAMETERS_', i_val=functional_id)
      IF (functional_id == xc_none) THEN
         ALLOCATE (operators, MOLD=matrices)
         ALLOCATE (operators_direction, MOLD=matrices)
         operators = CMPLX(0.0_dp, 0.0_dp, dp)
         operators_direction = CMPLX(0.0_dp, 0.0_dp, dp)
         energy = 0.0_dp
         RETURN
      END IF
      CALL qs_upf_inverse_aux_xc_evaluate(qs_env, xc, grid_energy, atom_energy, operators, density_matrices=matrices)
      energy = grid_energy + atom_energy
      ns = SIZE(matrices, 4)
      use_tau = ALLOCATED(state%kinetic_cache)
      fine = .NOT. pw_grid_compare(native_pool%pw_grid, pool%pw_grid)
      ALLOCATE (density(ns), kinetic(ns), direction(ns), kinetic_direction(ns), combined(ns), combined_tau(ns), &
                density_g(ns), direction_g(ns))
      DO spin = 1, ns
         CALL density(spin)%create(pool%pw_grid)
         CALL density_g(spin)%create(pool%pw_grid)
         CALL kinetic(spin)%create(pool%pw_grid)
         CALL direction(spin)%create(pool%pw_grid)
         CALL direction_g(spin)%create(pool%pw_grid)
         CALL kinetic_direction(spin)%create(pool%pw_grid)
         CALL combined(spin)%create(pool%pw_grid)
         CALL combined_tau(spin)%create(pool%pw_grid)
      END DO
      CALL qs_upf_inverse_scf_xc_source_direction(ks_env, matrices, density, kinetic, context=state)
      CALL qs_upf_inverse_scf_xc_source_direction(ks_env, matrices_direction, direction, kinetic_direction, context=state)
      DO spin = 1, ns
         CALL pw_transfer(density(spin), density_g(spin))
         CALL pw_transfer(direction(spin), direction_g(spin))
      END DO
      IF (ASSOCIATED(core)) THEN
         CPASSERT(ASSOCIATED(native_core_g))
         CALL core_g%create(pool%pw_grid)
         CALL transferred_core%create(pool%pw_grid)
         CALL pw_real_transfer(native_core_g, core_g)
         IF (fine) THEN
            CALL pw_transfer(core_g, transferred_core)
            CALL pw_transfer(transferred_core, core_g)
         ELSE
            transferred_core%array = core%array
         END IF
         DO spin = 1, ns
            density(spin)%array = density(spin)%array + transferred_core%array
            density_g(spin)%array = density_g(spin)%array + core_g%array
         END DO
         CALL transferred_core%release()
         CALL core_g%release()
      END IF
      IF (use_tau) CALL qs_upf_inverse_core_kinetic(ks_env, kinetic, 1.0_dp, context=state)
      weights => native_weights
      IF (fine .AND. ASSOCIATED(native_weights)) THEN
         CALL weight_source_g%create(native_weights%pw_grid)
         CALL weight_target_g%create(pool%pw_grid)
         CALL transferred_weights%create(pool%pw_grid)
         CALL pw_transfer(native_weights, weight_source_g)
         CALL pw_real_transfer(weight_source_g, weight_target_g)
         CALL pw_transfer(weight_target_g, transferred_weights)
         weights => transferred_weights
         CALL weight_source_g%release()
         CALL weight_target_g%release()
      END IF
      NULLIFY (response, response_tau)
      CALL qs_upf_inverse_xc_grid_kernel(xc, density, density_g, kinetic, direction, direction_g, kinetic_direction, &
                                         weights, pool, response, response_tau)
      atomic_xc = control%qs_control%gapw .OR. control%qs_control%gapw_xc .OR. ALLOCATED(state%paw_reconstruction)
      IF (atomic_xc) THEN
         CALL qs_upf_inverse_scf_xc_model_response(ks_env, matrices_direction, direction, kinetic_direction, xc, &
                                                   atomic, combined, combined_tau, context=state)
      ELSE
         DO spin = 1, ns
            combined(spin)%array = 0.0_dp
            combined_tau(spin)%array = 0.0_dp
         END DO
      END IF
      DO spin = 1, ns
         combined(spin)%array = combined(spin)%array + response(spin)%array
         combined_tau(spin)%array = combined_tau(spin)%array + response_tau(spin)%array
      END DO
      CALL qs_upf_inverse_scf_xc_source_adjoint(ks_env, combined, combined_tau, operators_direction, context=state)
      IF (atomic_xc) operators_direction(:, :, :, :) = operators_direction + atomic
      IF (ANY(state%paw%active) .AND. .NOT. ALLOCATED(state%paw_reconstruction)) THEN
         ALLOCATE (paw, MOLD=matrices)
         IF (ASSOCIATED(state%bloch%kpoints)) THEN
            ALLOCATE (kweights, SOURCE=state%bloch%weights)
            CALL qs_upf_inverse_paw_xc_kernel(state%paw, state%bloch%frames, matrices, matrices_direction, kweights, xc, paw)
         ELSE
            CALL qs_upf_inverse_paw_xc_kernel(state%paw, [state%frame], matrices, matrices_direction, [1.0_dp], xc, paw)
         END IF
         operators_direction(:, :, :, :) = operators_direction + paw
      END IF
      DO spin = 1, ns
         CALL pool%give_back_pw(response(spin))
         CALL pool%give_back_pw(response_tau(spin))
         CALL density(spin)%release()
         CALL density_g(spin)%release()
         CALL kinetic(spin)%release()
         CALL direction(spin)%release()
         CALL direction_g(spin)%release()
         CALL kinetic_direction(spin)%release()
         CALL combined(spin)%release()
         CALL combined_tau(spin)%release()
      END DO
      DEALLOCATE (response, response_tau)
      IF (fine .AND. ASSOCIATED(native_weights)) CALL transferred_weights%release()
   END SUBROUTINE xc_component_kernel

END MODULE qs_upf_inverse_aux_xc
