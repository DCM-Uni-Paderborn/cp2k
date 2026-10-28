!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines to handle an external electrostatic field
!>        The external field can be generic and is provided by user input
! **************************************************************************************************
MODULE qs_external_potential
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type,&
                                              dbcsr_type
   USE cp_log_handling,                 ONLY: cp_to_string
   USE cp_realspace_grid_cube,          ONLY: cp_cube_to_pw
   USE fparser,                         ONLY: evalf,&
                                              evalfd,&
                                              finalizef,&
                                              initf,&
                                              parsef
   USE generic_function_input,          ONLY: parse_function_section
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp,&
                                              int_8
   USE kpoint_methods,                  ONLY: kpoint_star_operator_defect
   USE mathconstants,                   ONLY: fourpi
   USE mathlib,                         ONLY: invmat
   USE maxwell_solver_interface,        ONLY: maxwell_solver
   USE message_passing,                 ONLY: mp_para_env_type
   USE particle_types,                  ONLY: particle_type
   USE pw_grid_types,                   ONLY: PW_MODE_LOCAL
   USE pw_methods,                      ONLY: pw_integral_ab,&
                                              pw_zero
   USE pw_types,                        ONLY: pw_r3d_rs_type
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_upf_inverse_background,       ONLY: qs_upf_inverse_background_integrate
   USE qs_upf_inverse_coulomb,          ONLY: qs_upf_inverse_local_potential
   USE qs_upf_inverse_grid,             ONLY: qs_upf_inverse_grid_integrate_dense,&
                                              qs_upf_inverse_matrix_gather
   USE qs_upf_inverse_kpoints,          ONLY: qs_upf_inverse_kpoints_density
   USE qs_upf_inverse_scf,              ONLY: qs_upf_inverse_scf_integrate_background
   USE qs_upf_inverse_types,            ONLY: qs_upf_inverse_scf_type
   USE string_utilities,                ONLY: compress
   USE virial_types,                    ONLY: virial_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_external_potential'

! *** Public subroutines ***
   PUBLIC :: external_e_potential, &
             external_c_potential, &
             external_upf_potential, &
             interpolate_external_potential

CONTAINS

! **************************************************************************************************
!> \brief Electronic external-field corrections for the native USPP/PAW orbital representations.
!>        Ordinary grid density, operator and orbital response remain native. This adds signed
!>        GAPW centers and original PAW corrections through the same scalar-potential transpose.
!>        Geometry also restores the external grid measure and differentiates a prescribed
!>        Cartesian function on the moving grid. Cube/Maxwell samples remain fixed fractionally.
!> \param qs_env Current native field and density owners.
!> \param energy Correction energy, excluding effective point-ion/core counterterms.
!> \param hamiltonian Optional native spin matrix receiving the additive AO correction.
!> \param spin Spin index required with hamiltonian.
!> \param response Optional complete correction to the native folded-potential geometry derivative.
!> \param atom_response ...
!> \param cell_response ...
!> \param coefficients_response ...
! **************************************************************************************************
   SUBROUTINE external_upf_potential(qs_env, energy, hamiltonian, spin, response, atom_response, &
                                     cell_response, coefficients_response)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(dp), INTENT(OUT)                              :: energy
      TYPE(dbcsr_type), INTENT(INOUT), OPTIONAL          :: hamiltonian
      INTEGER, INTENT(IN), OPTIONAL                      :: spin
      REAL(dp), INTENT(OUT), OPTIONAL                    :: response
      REAL(dp), INTENT(IN), OPTIONAL                     :: atom_response(:, :), cell_response(3, 3)
      COMPLEX(dp), INTENT(IN), OPTIONAL                  :: coefficients_response(:, :, :)

      COMPLEX(dp), ALLOCATABLE                           :: operators(:, :, :)
      INTEGER                                            :: bounds(2, 3), changed, i, info, ispin, &
                                                            j, k, n
      REAL(dp)                                           :: inverse_cell(3, 3), strain(3, 3), &
                                                            trace_strain
      REAL(dp), ALLOCATABLE                              :: gradient(:, :), points(:, :)
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: control
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(pw_r3d_rs_type)                               :: field_response, potential
      TYPE(pw_r3d_rs_type), POINTER                      :: density(:), vee
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_upf_inverse_scf_type), POINTER             :: inverse
      TYPE(section_vals_type), POINTER                   :: input, section

      energy = 0.0_dp
      IF (PRESENT(response)) response = 0.0_dp
      CALL get_qs_env(qs_env, dft_control=control, ks_env=ks_env, vee=vee, cell=cell, rho=rho, para_env=para_env)
      IF (.NOT. control%apply_external_potential) RETURN
      CALL get_ks_env(ks_env, upf_inverse=inverse)
      IF (.NOT. ASSOCIATED(inverse)) RETURN
      IF (.NOT. inverse%enabled) RETURN
      CALL potential%create(vee%pw_grid)
      potential%array = -vee%array
      IF (PRESENT(response)) THEN
         CPASSERT(PRESENT(cell_response))
         CPASSERT(PRESENT(atom_response))
         CPASSERT(PRESENT(coefficients_response))
         inverse_cell = cell%hmat
         CALL invmat(inverse_cell, info)
         CPASSERT(info == 0)
         strain = MATMUL(cell_response, inverse_cell)
         trace_strain = strain(1, 1) + strain(2, 2) + strain(3, 3)
         CALL field_response%create(vee%pw_grid)
         field_response%array = 0.0_dp
         IF (ANY(strain /= 0.0_dp) .AND. .NOT. control%expot_control%read_from_cube .AND. &
             .NOT. control%expot_control%maxwell_solver) THEN
            bounds = vee%pw_grid%bounds_local
            n = PRODUCT(bounds(2, :) - bounds(1, :) + 1)
            ALLOCATE (points(3, n), gradient(3, n))
            n = 0
            DO k = bounds(1, 3), bounds(2, 3)
               DO j = bounds(1, 2), bounds(2, 2)
                  DO i = bounds(1, 1), bounds(2, 1)
                     n = n + 1
                     points(:, n) = MATMUL(vee%pw_grid%dh, REAL([i, j, k] - vee%pw_grid%bounds(1, :), dp))
                  END DO
               END DO
            END DO
            CALL get_qs_env(qs_env, input=input)
            section => section_vals_get_subs_vals(input, "DFT%EXTERNAL_POTENTIAL")
            CALL get_external_potential(points, section, dfunc=gradient, calc_derivatives=.TRUE.)
            n = 0
            DO k = bounds(1, 3), bounds(2, 3)
               DO j = bounds(1, 2), bounds(2, 2)
                  DO i = bounds(1, 1), bounds(2, 1)
                     n = n + 1
                     field_response%array(i, j, k) = -SUM(gradient(:, n)*MATMUL(strain, points(:, n)))
                  END DO
               END DO
            END DO
         END IF
         CALL qs_upf_inverse_local_potential(qs_env, potential, energy, response=response, atom_response=atom_response, &
                        cell_response=cell_response, coefficients_response=coefficients_response, potential_response=field_response)
         CALL qs_rho_get(rho, rho_r=density)
         DO ispin = 1, SIZE(density)
            response = response + trace_strain*pw_integral_ab(density(ispin), potential) + &
                       pw_integral_ab(density(ispin), field_response)
         END DO
         CALL field_response%release()
      ELSE
         changed = 1
         IF (ALLOCATED(inverse%external_operator) .AND. ALLOCATED(inverse%external_samples)) THEN
            IF (ALL(SHAPE(inverse%external_samples) == SHAPE(potential%array))) THEN
               IF (ALL(inverse%external_samples == potential%array)) changed = 0
            END IF
         END IF
         ! A change on any slab invalidates every rank before the collective atomic transpose.
         CALL para_env%max(changed)
         IF (changed /= 0) THEN
            CALL qs_upf_inverse_local_potential(qs_env, potential, energy, operators, &
                                                core_energy=inverse%external_core_energy)
            CALL MOVE_ALLOC(operators, inverse%external_operator)
            IF (ALLOCATED(inverse%external_samples)) DEALLOCATE (inverse%external_samples)
            ALLOCATE (inverse%external_samples, SOURCE=potential%array)
            CALL external_upf_check_symmetry(ks_env, potential, control%qs_control%eps_rho_gspace)
         END IF
         CALL external_upf_trace(ks_env, inverse%external_operator, energy)
         energy = energy + inverse%external_core_energy
         IF (PRESENT(hamiltonian)) THEN
            CPASSERT(PRESENT(spin))
            CALL qs_upf_inverse_scf_integrate_background(ks_env, inverse%external_operator, hamiltonian, spin, additive=.TRUE.)
         END IF
      END IF
      CALL potential%release()
   END SUBROUTINE external_upf_potential

! **************************************************************************************************
!> \brief Trace of a fixed scalar-field operator with the current native spin densities.
!>        Occupations are unweighted in each physical Bloch slice. Mesh weights enter once here.
!> \param ks_env Current native density and inverse frame owners.
!> \param operators Spin-free unweighted operator, indexed by physical k contribution.
!> \param energy Replicated electronic coupling, excluding the frozen-core constant.
! **************************************************************************************************
   SUBROUTINE external_upf_trace(ks_env, operators, energy)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      COMPLEX(dp), INTENT(IN)                            :: operators(:, :, :)
      REAL(dp), INTENT(OUT)                              :: energy

      COMPLEX(dp), ALLOCATABLE                           :: gamma(:, :), matrix(:, :, :)
      INTEGER                                            :: k, spin
      TYPE(dbcsr_p_type), POINTER                        :: density(:, :)
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_upf_inverse_scf_type), POINTER             :: inverse

      CALL get_ks_env(ks_env, rho=rho, upf_inverse=inverse)
      CALL qs_rho_get(rho, rho_ao_kp=density)
      energy = 0.0_dp
      DO spin = 1, SIZE(density, 1)
         IF (ASSOCIATED(inverse%bloch%kpoints)) THEN
            CALL qs_upf_inverse_kpoints_density(inverse%bloch, spin, matrix)
            DO k = 1, SIZE(matrix, 3)
               energy = energy + inverse%bloch%weights(k)*REAL(SUM(CONJG(operators(:, :, k))*matrix(:, :, k)), dp)
            END DO
         ELSE
            CALL qs_upf_inverse_matrix_gather(inverse%frame%ao, density(spin, 1)%matrix, gamma)
            energy = energy + REAL(SUM(CONJG(operators(:, :, 1))*gamma), dp)
         END IF
      END DO
   END SUBROUTINE external_upf_trace

! **************************************************************************************************
!> \brief Check the complete applied operator before averaging inequivalent Bloch contributions.
!>        Use the density-independent field alone, so an unconverged Hartree/XC potential cannot
!>        falsely reject a compatible applied field. The native AO adjoint includes antiunitary maps.
!> \param ks_env Prepared inverse geometry and native Gaussian density maps.
!> \param potential Scalar electron potential without the grid volume factor.
!> \param eps_grid Native integration accuracy used with a roundoff floor for the AO comparison.
! **************************************************************************************************
   SUBROUTINE external_upf_check_symmetry(ks_env, potential, eps_grid)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: potential
      REAL(dp), INTENT(IN)                               :: eps_grid

      COMPLEX(dp), ALLOCATABLE                           :: background(:, :, :), matrix(:, :), &
                                                            operators(:, :, :)
      INTEGER                                            :: axis, k
      LOGICAL                                            :: little_group
      REAL(dp)                                           :: defect, scale, tolerance
      TYPE(dft_control_type), POINTER                    :: control
      TYPE(pw_r3d_rs_type)                               :: zero, zeros(3)
      TYPE(qs_upf_inverse_scf_type), POINTER             :: inverse

      CALL get_ks_env(ks_env, upf_inverse=inverse, dft_control=control)
      IF (.NOT. ASSOCIATED(inverse%bloch%kpoints)) RETURN
      IF (inverse%bloch%kpoints%full_grid) RETURN
      little_group = .NOT. control%use_gauxc
      IF (SIZE(inverse%bloch%weights) == inverse%bloch%kpoints%nkp .AND. .NOT. little_group) RETURN
      ALLOCATE (operators, SOURCE=inverse%external_operator)
      IF (inverse%gaussian_background) THEN
         ALLOCATE (background, MOLD=inverse%coulomb%metric)
         CALL qs_upf_inverse_background_integrate(inverse%coulomb%background, potential=potential, &
                                                  ao_operators=background, ks_env=ks_env)
         operators = operators + background
      ELSE
         IF (inverse%native_gaussian) THEN
            ALLOCATE (background, MOLD=inverse%gaussian_kinetic)
            CALL qs_upf_inverse_background_integrate(inverse%gaussian, potential=potential, &
                                                     ao_operators=background, ks_env=ks_env)
            operators = operators + background
         END IF
         CALL zero%create(potential%pw_grid)
         zero%array = 0.0_dp
         DO axis = 1, 3
            CALL zeros(axis)%create(potential%pw_grid)
            zeros(axis)%array = 0.0_dp
         END DO
         DO k = 1, SIZE(operators, 3)
            CALL qs_upf_inverse_grid_integrate_dense(inverse%bloch%frames(k), potential, zeros, zero, matrix, &
                                                     augmentation=inverse%bloch%augmentation, fourier_charge=.TRUE., &
                                                     orbital_fields=.NOT. inverse%native_gaussian)
            operators(:, :, k) = operators(:, :, k) + matrix
         END DO
         CALL zero%release()
         DO axis = 1, 3
            CALL zeros(axis)%release()
         END DO
      END IF
      CALL kpoint_star_operator_defect(inverse%bloch%star, operators, defect, scale, little_group=little_group)
      tolerance = MAX(1024.0_dp*EPSILON(1.0_dp), 32.0_dp*eps_grid)*MAX(1.0_dp, scale)
      IF (defect > tolerance) THEN
         CALL cp_abort(__LOCATION__, &
                       "The external potential is incompatible with the k-point symmetry reduction. "// &
                       "Use FULL_GRID TRUE and SYMMETRY FALSE.")
      END IF
   END SUBROUTINE external_upf_check_symmetry

! **************************************************************************************************
!> \brief  Computes the external potential on the grid
!> \param qs_env ...
!> \date   12.2009
!> \author Teodoro Laino [tlaino]
! **************************************************************************************************
   SUBROUTINE external_e_potential(qs_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'external_e_potential'

      INTEGER                                            :: handle, i, j, k
      INTEGER(kind=int_8)                                :: npoints
      INTEGER, DIMENSION(2, 3)                           :: bo_global, bo_local
      REAL(kind=dp)                                      :: scaling_factor
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: efunc
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :)        :: grid_p
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_r3d_rs_type), POINTER                      :: v_ee
      TYPE(section_vals_type), POINTER                   :: ext_pot_section, input

      CALL timeset(routineN, handle)
      CALL get_qs_env(qs_env, dft_control=dft_control)
      IF (dft_control%apply_external_potential) THEN
         IF (dft_control%eval_external_potential) THEN
            CALL get_qs_env(qs_env, vee=v_ee)
            IF (dft_control%expot_control%maxwell_solver) THEN
               scaling_factor = dft_control%expot_control%scaling_factor
               CALL maxwell_solver(dft_control%maxwell_control, v_ee, &
                                   qs_env%sim_step, qs_env%sim_time, &
                                   scaling_factor)
               dft_control%eval_external_potential = .FALSE.
            ELSE IF (dft_control%expot_control%read_from_cube) THEN
               scaling_factor = dft_control%expot_control%scaling_factor
               CALL cp_cube_to_pw(v_ee, 'pot.cube', scaling_factor)
               dft_control%eval_external_potential = .FALSE.
            ELSE
               CALL get_qs_env(qs_env, input=input)
               ext_pot_section => section_vals_get_subs_vals(input, "DFT%EXTERNAL_POTENTIAL")

               CALL pw_zero(v_ee)

               bo_local = v_ee%pw_grid%bounds_local
               bo_global = v_ee%pw_grid%bounds

               npoints = INT(bo_local(2, 1) - bo_local(1, 1) + 1, kind=int_8)* &
                         INT(bo_local(2, 2) - bo_local(1, 2) + 1, kind=int_8)* &
                         INT(bo_local(2, 3) - bo_local(1, 3) + 1, kind=int_8)
               ALLOCATE (efunc(npoints))
               ALLOCATE (grid_p(3, npoints))
               npoints = 0
               DO k = bo_local(1, 3), bo_local(2, 3)
                  DO j = bo_local(1, 2), bo_local(2, 2)
                     DO i = bo_local(1, 1), bo_local(2, 1)
                        npoints = npoints + 1
                        grid_p(:, npoints) = MATMUL(v_ee%pw_grid%dh, REAL([i, j, k] - bo_global(1, :), dp))
                     END DO
                  END DO
               END DO

               CALL get_external_potential(grid_p, ext_pot_section, func=efunc)

               npoints = 0
               DO k = bo_local(1, 3), bo_local(2, 3)
                  DO j = bo_local(1, 2), bo_local(2, 2)
                     DO i = bo_local(1, 1), bo_local(2, 1)
                        npoints = npoints + 1
                        v_ee%array(i, j, k) = v_ee%array(i, j, k) + efunc(npoints)
                     END DO
                  END DO
               END DO

               DEALLOCATE (grid_p, efunc)

               dft_control%eval_external_potential = .FALSE.
            END IF
         END IF
      END IF
      CALL timestop(handle)
   END SUBROUTINE external_e_potential

! **************************************************************************************************
!> \brief  Computes the force and the energy due to the external potential on the cores
!> \param qs_env ...
!> \param calculate_forces ...
!> \date   12.2009
!> \author Teodoro Laino [tlaino]
! **************************************************************************************************
   SUBROUTINE external_c_potential(qs_env, calculate_forces)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, OPTIONAL                                  :: calculate_forces

      CHARACTER(len=*), PARAMETER :: routineN = 'external_c_potential'

      INTEGER                                            :: atom_a, handle, iatom, ikind, natom, &
                                                            nkind, nparticles
      INTEGER, DIMENSION(:), POINTER                     :: list
      LOGICAL                                            :: inverse_active, my_force, pot_on_grid
      REAL(KIND=dp)                                      :: ee_core_ener, zeff
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: efunc
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: dfunc, r
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_r3d_rs_type), POINTER                      :: v_ee
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_upf_inverse_scf_type), POINTER             :: inverse
      TYPE(section_vals_type), POINTER                   :: ext_pot_section, input
      TYPE(virial_type), POINTER                         :: virial

      CALL timeset(routineN, handle)
      NULLIFY (dft_control)

      CALL get_qs_env(qs_env=qs_env, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      energy=energy, &
                      particle_set=particle_set, &
                      input=input, &
                      cell=cell, &
                      dft_control=dft_control)
      CALL get_qs_env(qs_env, ks_env=ks_env, para_env=para_env, virial=virial)
      CALL get_ks_env(ks_env, upf_inverse=inverse)
      inverse_active = .FALSE.
      IF (ASSOCIATED(inverse)) inverse_active = inverse%enabled

      IF (dft_control%apply_external_potential .AND. .NOT. qs_env%mimic) THEN
         !ensure that external potential is loaded to grid
         IF (dft_control%eval_external_potential) THEN
            CALL external_e_potential(qs_env)
         END IF
         my_force = .FALSE.
         IF (PRESENT(calculate_forces)) my_force = calculate_forces
         ee_core_ener = 0.0_dp
         nkind = SIZE(atomic_kind_set)

         !check if external potential on grid has been loaded from a file instead of giving a function
         IF (dft_control%expot_control%read_from_cube .OR. &
             dft_control%expot_control%maxwell_solver) THEN
            CALL get_qs_env(qs_env, vee=v_ee)
            pot_on_grid = .TRUE.
         ELSE
            pot_on_grid = .FALSE.
            ext_pot_section => section_vals_get_subs_vals(input, "DFT%EXTERNAL_POTENTIAL")
         END IF

         nparticles = 0
         DO ikind = 1, SIZE(atomic_kind_set)
            CALL get_atomic_kind(atomic_kind_set(ikind), natom=natom)
            nparticles = nparticles + MAX(natom, 0)
         END DO

         ALLOCATE (efunc(nparticles))
         ALLOCATE (dfunc(3, nparticles), r(3, nparticles))

         nparticles = 0
         DO ikind = 1, SIZE(atomic_kind_set)
            CALL get_atomic_kind(atomic_kind_set(ikind), atom_list=list, natom=natom)

            DO iatom = 1, natom
               atom_a = list(iatom)
               nparticles = nparticles + 1
               !pbc returns r(i) in range [-cell%hmat(i,i)/2, cell%hmat(i,i)/2]
               !for periodic dimensions (assuming the cell is orthorombic).
               !This is not consistent with the potential on grid, where r(i) is
               !in range [0, cell%hmat(i,i)]
               !Use new pbc function with switch positive_range=.TRUE.
               r(:, nparticles) = pbc(particle_set(atom_a)%r(:), cell, positive_range=.TRUE.)
            END DO
         END DO

         !if potential is on grid, interpolate the value at r,
         !otherwise evaluate the given function
         IF (pot_on_grid) THEN
            DO iatom = 1, nparticles
               CALL interpolate_external_potential(r(:, iatom), v_ee, func=efunc(iatom), &
                                                   dfunc=dfunc(:, iatom), calc_derivatives=my_force)
            END DO
         ELSE
            CALL get_external_potential(r, ext_pot_section, func=efunc, dfunc=dfunc, calc_derivatives=my_force)
         END IF

         IF (my_force) THEN
            CALL get_qs_env(qs_env=qs_env, force=force)
         END IF

         nparticles = 0
         DO ikind = 1, SIZE(atomic_kind_set)
            CALL get_atomic_kind(atomic_kind_set(ikind), natom=natom)
            CALL get_qs_kind(qs_kind_set(ikind), zeff=zeff)
            IF (inverse_active) THEN
               IF (inverse%paw(ikind)%active) THEN
                  ! Neutral frozen core: its tabulated grid charge and point countercharge must agree.
                  zeff = zeff + fourpi*SUM(inverse%paw(ikind)%xcgrid%wr*inverse%paw(ikind)%xcgrid%core(:, 1))
               END IF
            END IF

            DO iatom = 1, natom
               nparticles = nparticles + 1

               ee_core_ener = ee_core_ener + zeff*efunc(nparticles)
               IF (my_force) THEN
                  force(ikind)%eev(1:3, iatom) = dfunc(1:3, nparticles)*zeff
                  IF (inverse_active .AND. .NOT. pot_on_grid .AND. virial%pv_availability .AND. &
                      .NOT. virial%pv_numer .AND. para_env%is_source()) THEN
                     DO atom_a = 1, 3
                        virial%pv_virial(:, atom_a) = virial%pv_virial(:, atom_a) - &
                                                      zeff*dfunc(:, nparticles)*r(atom_a, nparticles)
                     END DO
                  END IF
               END IF
            END DO
         END DO
         energy%ee_core = ee_core_ener

         DEALLOCATE (dfunc, r)
         DEALLOCATE (efunc)
      END IF
      CALL timestop(handle)
   END SUBROUTINE external_c_potential

! **************************************************************************************************
!> \brief  Low level function for computing the potential and the derivatives
!> \param r                position in realspace for each grid-point
!> \param ext_pot_section ...
!> \param func             external potential at r
!> \param dfunc            derivative of the external potential at r
!> \param calc_derivatives Whether to calculate dfunc
!> \date   12.2009
!> \par History
!>      12.2009            created [tlaino]
!>      11.2014            reading external cube file added [Juha Ritala & Matt Watkins]
!> \author Teodoro Laino [tlaino]
! **************************************************************************************************
   SUBROUTINE get_external_potential(r, ext_pot_section, func, dfunc, calc_derivatives)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: r
      TYPE(section_vals_type), POINTER                   :: ext_pot_section
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT), OPTIONAL :: func
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT), &
         OPTIONAL                                        :: dfunc
      LOGICAL, INTENT(IN), OPTIONAL                      :: calc_derivatives

      CHARACTER(len=*), PARAMETER :: routineN = 'get_external_potential'

      CHARACTER(LEN=default_path_length)                 :: coupling_function
      CHARACTER(LEN=default_string_length)               :: def_error, this_error
      CHARACTER(LEN=default_string_length), &
         ALLOCATABLE, DIMENSION(:)                       :: my_par
      INTEGER                                            :: handle, j
      INTEGER(kind=int_8)                                :: ipoint, npoints
      LOGICAL                                            :: check, my_force
      REAL(KIND=dp)                                      :: dedf, dx, err, lerr
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: my_val

      CALL timeset(routineN, handle)
      my_force = .FALSE.
      IF (PRESENT(calc_derivatives)) my_force = calc_derivatives
      check = PRESENT(dfunc) .EQV. PRESENT(calc_derivatives)
      CPASSERT(check)
      CALL parse_function_section(ext_pot_section, "FUNCTION", coupling_function, &
                                  my_par, my_val, i_rep_sec=1, &
                                  input_variables=["X", "Y", "Z"], dx=dx, error_limit=lerr)
      CALL initf(1)
      CALL parsef(1, TRIM(coupling_function), my_par)

      npoints = SIZE(r, 2, kind=int_8)

      DO ipoint = 1, npoints
         my_val(1) = r(1, ipoint)
         my_val(2) = r(2, ipoint)
         my_val(3) = r(3, ipoint)

         IF (PRESENT(func)) func(ipoint) = evalf(1, my_val)
         IF (my_force) THEN
            DO j = 1, 3
               dedf = evalfd(1, j, my_val, dx, err)
               IF (ABS(err) > lerr) THEN
                  WRITE (this_error, "(A,G12.6,A)") "(", err, ")"
                  WRITE (def_error, "(A,G12.6,A)") "(", lerr, ")"
                  CALL compress(this_error, .TRUE.)
                  CALL compress(def_error, .TRUE.)
                  CALL cp_warn(__LOCATION__, &
                               'ASSERTION (cond) failed at line '//cp_to_string(__LINE__)// &
                               ' Error '//TRIM(this_error)//' in computing numerical derivatives larger than'// &
                               TRIM(def_error)//' .')
               END IF
               dfunc(j, ipoint) = dedf
            END DO
         END IF
      END DO
      DEALLOCATE (my_par)
      DEALLOCATE (my_val)
      CALL finalizef()
      CALL timestop(handle)
   END SUBROUTINE get_external_potential

! **************************************************************************************************
!> \brief Trilinear external potential and the Cartesian derivative of the same interpolant.
!>        Periodic image wrapping and local/distributed grids use the full cell matrix.
!> \param r                 ...
!> \param grid             external potential pw grid, vee
!> \param func             value of vee at r
!> \param dfunc            derivatives of vee at r
!> \param calc_derivatives calc dfunc
! **************************************************************************************************
   SUBROUTINE interpolate_external_potential(r, grid, func, dfunc, calc_derivatives)
      REAL(KIND=dp), INTENT(IN)                          :: r(3)
      TYPE(pw_r3d_rs_type), POINTER                      :: grid
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: func, dfunc(3)
      LOGICAL, INTENT(IN), OPTIONAL                      :: calc_derivatives

      CHARACTER(len=*), PARAMETER :: routineN = 'interpolate_external_potential'

      INTEGER                                            :: axis, corner(3), grid_index(3), handle, &
                                                            i, info, j, k, lower(3), node(3), &
                                                            upper(3)
      LOGICAL                                            :: check, my_force
      REAL(dp)                                           :: buffer(0:1, 0:1, 0:1), derivative(3), &
                                                            grid_fraction(3), inverse(3, 3), &
                                                            point(3), value, weight(3)

      CALL timeset(routineN, handle)
      my_force = .FALSE.
      IF (PRESENT(calc_derivatives)) my_force = calc_derivatives
      check = PRESENT(dfunc) .EQV. PRESENT(calc_derivatives)
      CPASSERT(check)
      inverse = grid%pw_grid%dh
      CALL invmat(inverse, info)
      CPASSERT(info == 0)
      point = MATMUL(inverse, r)
      node = FLOOR(point)
      grid_fraction = point - REAL(node, dp)
      lower = grid%pw_grid%bounds(1, :)
      upper = grid%pw_grid%bounds(2, :)
      buffer = 0.0_dp
      ! Collect just the eight interpolation nodes, including wrapped local-grid boundaries.
      DO k = 0, 1
         DO j = 0, 1
            DO i = 0, 1
               grid_index = pbc_index(lower + node + [i, j, k], lower, upper)
               IF (ALL(grid_index >= grid%pw_grid%bounds_local(1, :)) .AND. &
                   ALL(grid_index <= grid%pw_grid%bounds_local(2, :))) THEN
                  buffer(i, j, k) = grid%array(grid_index(1), grid_index(2), grid_index(3))
               END IF
            END DO
         END DO
      END DO
      IF (grid%pw_grid%para%mode /= PW_MODE_LOCAL) CALL grid%pw_grid%para%group%sum(buffer)
      value = 0.0_dp
      derivative = 0.0_dp
      DO k = 0, 1
         DO j = 0, 1
            DO i = 0, 1
               corner = [i, j, k]
               weight = MERGE(grid_fraction, 1.0_dp - grid_fraction, corner == 1)
               value = value + buffer(i, j, k)*PRODUCT(weight)
               IF (.NOT. my_force) CYCLE
               DO axis = 1, 3
                  weight = MERGE(grid_fraction, 1.0_dp - grid_fraction, corner == 1)
                  weight(axis) = REAL(2*corner(axis) - 1, dp)
                  derivative(axis) = derivative(axis) + buffer(i, j, k)*PRODUCT(weight)
               END DO
            END DO
         END DO
      END DO
      IF (PRESENT(func)) func = value
      ! Differentiate the same trilinear energy and transform its grid derivative to Cartesian space.
      IF (my_force) dfunc = MATMUL(TRANSPOSE(inverse), derivative)
      CALL timestop(handle)
   END SUBROUTINE interpolate_external_potential

! **************************************************************************************************
!> \brief          get a correct value for possible out of bounds index using periodic
!>                  boundary conditions
!> \param i ...
!> \param lowbound ...
!> \param upbound ...
!> \return ...
! **************************************************************************************************
   ELEMENTAL FUNCTION pbc_index(i, lowbound, upbound)
      INTEGER, INTENT(IN)                                :: i, lowbound, upbound
      INTEGER                                            :: pbc_index

      pbc_index = lowbound + MODULO(i - lowbound, upbound - lowbound + 1)
   END FUNCTION pbc_index

END MODULE qs_external_potential
