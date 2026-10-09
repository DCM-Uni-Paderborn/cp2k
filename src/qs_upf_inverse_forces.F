!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Native forces and virial from the stationary inverse USPP density and metric responses.
!>        Ionic/core forces and explicit Poisson/XC virial terms retain their native consumers.
! **************************************************************************************************
MODULE qs_upf_inverse_forces
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE cell_types,                      ONLY: cell_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE particle_types,                  ONLY: particle_type
   USE pw_methods,                      ONLY: pw_integral_ab
   USE pw_types,                        ONLY: pw_r3d_rs_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_upf_inverse_background,       ONLY: qs_upf_inverse_background_geometry
   USE qs_upf_inverse_coulomb,          ONLY: qs_upf_inverse_coulomb_metric_response,&
                                              qs_upf_inverse_coulomb_response
   USE qs_upf_inverse_grid,             ONLY: qs_upf_inverse_matrix_gather
   USE qs_upf_inverse_kpoints,          ONLY: qs_upf_inverse_kpoints_density
   USE qs_upf_inverse_scf,              ONLY: qs_upf_inverse_scf_response
   USE qs_upf_inverse_scf_xc,           ONLY: qs_upf_inverse_scf_xc_response
   USE qs_upf_inverse_types,            ONLY: qs_upf_inverse_scf_type
   USE virial_types,                    ONLY: virial_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: qs_upf_inverse_force

CONTAINS

! **************************************************************************************************
!> \brief Add the stationary inverse electronic gradients to the native force and virial arrays.
!>        Replicated contractions contribute on one global rank, before the native final reduction.
!>        Remove the folded-potential operator measure before adding native functional volume terms.
!> \param qs_env Native converged environment with prepared P and energy-weighted W.
!> \param potential Folded native spin potentials, including their native dvol factor.
!> \param kinetic_potential Optional unweighted native spin potentials conjugate to augmented tau.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_force(qs_env, potential, kinetic_potential)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: potential(:)
      TYPE(pw_r3d_rs_type), INTENT(IN), OPTIONAL         :: kinetic_potential(:)

      COMPLEX(dp), ALLOCATABLE                           :: spin_density(:, :, :, :), &
                                                            total_density(:, :, :), &
                                                            total_weighted(:, :, :)
      INTEGER                                            :: atom, axis, index, kind, nk, nspin, &
                                                            other, spin
      INTEGER, ALLOCATABLE                               :: atom_of_kind(:), kind_of(:)
      LOGICAL                                            :: complete_coulomb
      REAL(dp)                                           :: field_measure, metric_response, &
                                                            response, source_virial(3, 3), &
                                                            strain(3, 3), work_virial(3, 3), &
                                                            xc_response
      REAL(dp), ALLOCATABLE                              :: atom_response(:, :), positions(:, :), &
                                                            reduced_k_response(:, :), &
                                                            source_forces(:, :), work_forces(:, :)
      TYPE(atomic_kind_type), POINTER                    :: atomic_kinds(:)
      TYPE(cell_type), POINTER                           :: cell
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), POINTER                       :: atoms(:)
      TYPE(pw_r3d_rs_type), ALLOCATABLE                  :: vdensity(:), vgradient(:, :), vtau(:)
      TYPE(pw_r3d_rs_type), POINTER                      :: density(:), kinetic_density(:)
      TYPE(qs_force_type), POINTER                       :: force(:)
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_upf_inverse_scf_type), POINTER             :: state
      TYPE(virial_type), POINTER                         :: virial

      CALL get_qs_env(qs_env, ks_env=ks_env, rho=rho, force=force, cell=cell, virial=virial, &
                      particle_set=atoms, atomic_kind_set=atomic_kinds, para_env=para_env)
      CALL get_ks_env(ks_env, upf_inverse=state)
      CPASSERT(ASSOCIATED(state) .AND. ASSOCIATED(force))
      complete_coulomb = ASSOCIATED(state%coulomb)
      CALL get_atomic_kind_set(atomic_kinds, atom_of_kind=atom_of_kind, kind_of=kind_of)
      nspin = SIZE(potential)
      IF (PRESENT(kinetic_potential)) THEN
         CPASSERT(SIZE(kinetic_potential) == nspin .AND. ALLOCATED(state%kinetic_cache))
      END IF
      nk = 1
      IF (ASSOCIATED(state%bloch%kpoints)) nk = SIZE(state%bloch%weights)
      ALLOCATE (atom_response(3, SIZE(atoms)), positions(3, SIZE(atoms)), reduced_k_response(3, nk), &
                vdensity(nspin), vgradient(3, nspin), vtau(nspin))
      reduced_k_response = 0.0_dp
      DO atom = 1, SIZE(atoms)
         positions(:, atom) = atoms(atom)%r
      END DO
      DO spin = 1, nspin
         CALL vdensity(spin)%create(potential(spin)%pw_grid)
         vdensity(spin)%array = potential(spin)%array/potential(spin)%pw_grid%dvol
         CALL vtau(spin)%create(potential(spin)%pw_grid)
         vtau(spin)%array = 0.0_dp
         IF (PRESENT(kinetic_potential)) THEN
            CPASSERT(ASSOCIATED(kinetic_potential(spin)%pw_grid, potential(spin)%pw_grid))
            vtau(spin)%array = kinetic_potential(spin)%array
         END IF
         DO axis = 1, 3
            CALL vgradient(axis, spin)%create(potential(spin)%pw_grid)
            vgradient(axis, spin)%array = 0.0_dp
         END DO
      END DO
      IF (complete_coulomb) THEN
         CALL stationary_matrices(ks_env, state, total_density, total_weighted, spin_density)
         ALLOCATE (source_forces(3, SIZE(atoms)), work_forces(3, SIZE(atoms)))
         source_forces = 0.0_dp
         source_virial = 0.0_dp
         DO spin = 1, nspin
            CALL qs_upf_inverse_background_geometry(state%coulomb%background, qs_env, &
                                                    spin_density(:, :, :, spin), vdensity(spin), &
                                                    work_forces, work_virial)
            source_forces(:, :) = source_forces + work_forces
            source_virial(:, :) = source_virial + work_virial
            IF (PRESENT(kinetic_potential)) THEN
               CALL qs_upf_inverse_background_geometry(state%coulomb%background, qs_env, &
                                                       spin_density(:, :, :, spin), vtau(spin), &
                                                       work_forces, work_virial, kinetic=.TRUE.)
               source_forces(:, :) = source_forces + work_forces
               source_virial(:, :) = source_virial + work_virial
            END IF
         END DO
      END IF
      DO atom = 1, SIZE(atoms)
         kind = kind_of(atom)
         index = atom_of_kind(atom)
         DO axis = 1, 3
            atom_response = 0.0_dp
            atom_response(axis, atom) = 1.0_dp
            IF (complete_coulomb) THEN
               CALL complete_response(qs_env, ks_env, state, total_density, total_weighted, vgradient, vgradient(1, :), &
                                      atom_response, 0.0_dp*cell%hmat, reduced_k_response, response, metric_response)
               response = response - SUM(source_forces*atom_response)
            ELSE
               CALL qs_upf_inverse_scf_response(ks_env, vdensity, vgradient, vtau, atom_response, &
                                                0.0_dp*cell%hmat, reduced_k_response, response, metric_response)
            END IF
            CALL qs_upf_inverse_scf_xc_response(ks_env, atom_response, 0.0_dp*cell%hmat, reduced_k_response, xc_response)
            response = response + xc_response
            IF (para_env%is_source()) THEN
               force(kind)%all_potential(axis, index) = force(kind)%all_potential(axis, index) + response - metric_response
               force(kind)%overlap(axis, index) = force(kind)%overlap(axis, index) + metric_response
            END IF
         END DO
      END DO
      IF (virial%pv_availability .AND. .NOT. virial%pv_numer) THEN
         IF (ANY(cell%perd == 0)) THEN
            CPABORT("Inverse analytic stress requires periodic electrostatics in all directions")
         END IF
         CALL qs_rho_get(rho, rho_r=density, tau_r=kinetic_density)
         CPASSERT(SIZE(density) == nspin)
         field_measure = 0.0_dp
         DO spin = 1, nspin
            field_measure = field_measure + pw_integral_ab(vdensity(spin), density(spin))
            IF (PRESENT(kinetic_potential)) THEN
               CPASSERT(ASSOCIATED(kinetic_density))
               field_measure = field_measure + pw_integral_ab(vtau(spin), kinetic_density(spin))
            END IF
         END DO
         DO axis = 1, 3
            DO other = 1, 3
               strain = 0.0_dp
               strain(axis, other) = 1.0_dp
               atom_response(:, :) = MATMUL(strain, positions)
               IF (complete_coulomb) THEN
                  CALL complete_response(qs_env, ks_env, state, total_density, total_weighted, vgradient, vgradient(1, :), &
                                         atom_response, MATMUL(strain, cell%hmat), reduced_k_response, response, metric_response)
                  ! The folded XC source has fixed fractional rows. Native XC owns its explicit measure.
                  response = response - SUM(source_virial*strain)
               ELSE
                  CALL qs_upf_inverse_scf_response(ks_env, vdensity, vgradient, vtau, atom_response, &
                                                   MATMUL(strain, cell%hmat), reduced_k_response, response, metric_response)
               END IF
               CALL qs_upf_inverse_scf_xc_response(ks_env, atom_response, MATMUL(strain, cell%hmat), &
                                                   reduced_k_response, xc_response)
               response = response + xc_response
               IF (axis == other .AND. .NOT. complete_coulomb) response = response - field_measure
               IF (para_env%is_source()) THEN
                  virial%pv_virial(axis, other) = virial%pv_virial(axis, other) - response
                  virial%pv_overlap(axis, other) = virial%pv_overlap(axis, other) - metric_response
               END IF
            END DO
         END DO
      END IF
      DO spin = 1, nspin
         CALL vdensity(spin)%release()
         CALL vtau(spin)%release()
         DO axis = 1, 3
            CALL vgradient(axis, spin)%release()
         END DO
      END DO
   END SUBROUTINE qs_upf_inverse_force

! **************************************************************************************************
!> \brief Spin-summed accepted and energy-weighted densities in the complete metric's global k order.
!> \param ks_env Native converged owners, including the prepared energy-weighted density service.
!> \param state Same defining inverse ordering as the installed complete charge metric.
!> \param density New Hermitian occupation matrices, without mesh weights.
!> \param weighted New Hermitian energy-weighted matrices, without mesh weights.
!> \param spin_density Same accepted occupation matrices retained separately for each XC spin source.
! **************************************************************************************************
   SUBROUTINE stationary_matrices(ks_env, state, density, weighted, spin_density)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(qs_upf_inverse_scf_type), INTENT(IN)          :: state
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: density(:, :, :), weighted(:, :, :), &
                                                            spin_density(:, :, :, :)

      COMPLEX(dp), ALLOCATABLE                           :: component(:, :, :), gamma(:, :)
      INTEGER                                            :: spin
      TYPE(dbcsr_p_type), POINTER                        :: density_gamma(:, :), weighted_gamma(:)
      TYPE(qs_rho_type), POINTER                         :: rho

      ALLOCATE (density, MOLD=state%coulomb%metric)
      ALLOCATE (weighted, MOLD=state%coulomb%metric)
      density = CMPLX(0.0_dp, 0.0_dp, dp)
      weighted = CMPLX(0.0_dp, 0.0_dp, dp)
      CALL get_ks_env(ks_env, rho=rho)
      CALL qs_rho_get(rho, rho_ao_kp=density_gamma)
      ALLOCATE (spin_density(SIZE(density, 1), SIZE(density, 2), SIZE(density, 3), SIZE(density_gamma, 1)))
      IF (.NOT. ASSOCIATED(state%bloch%kpoints)) THEN
         CALL get_ks_env(ks_env, matrix_w=weighted_gamma)
         CPASSERT(ASSOCIATED(weighted_gamma))
         CPASSERT(SIZE(weighted_gamma) == SIZE(density_gamma, 1))
      END IF
      DO spin = 1, SIZE(density_gamma, 1)
         IF (ASSOCIATED(state%bloch%kpoints)) THEN
            CALL qs_upf_inverse_kpoints_density(state%bloch, spin, component)
            spin_density(:, :, :, spin) = component
            density(:, :, :) = density + component
            CALL qs_upf_inverse_kpoints_density(state%bloch, spin, component, energy_weighted=.TRUE.)
            weighted(:, :, :) = weighted + component
         ELSE
            CALL qs_upf_inverse_matrix_gather(state%frame%ao, density_gamma(spin, 1)%matrix, gamma)
            spin_density(:, :, 1, spin) = gamma
            density(:, :, 1) = density(:, :, 1) + gamma
            CALL qs_upf_inverse_matrix_gather(state%frame%ao, weighted_gamma(spin)%matrix, gamma)
            weighted(:, :, 1) = weighted(:, :, 1) + gamma
         END IF
      END DO
   END SUBROUTINE stationary_matrices

! **************************************************************************************************
!> \brief Core, complete Coulomb and matching metric response of the stationary functional.
!>        The old radial-Q constraint and PW density potential are replaced. XC source and atom
!>        residual responses are added by the caller. Ionic self/overlap corrections remain native.
!> \param qs_env Native owners of the complete Coulomb source.
!> \param ks_env Converged inverse core and energy-weighted density owners.
!> \param state Complete metric, defining inverse and Coulomb source in the same native ordering.
!> \param density Spin-summed Hermitian occupations per global k, without mesh weights.
!> \param weighted Spin-summed Hermitian energy-weighted matrices in the same ordering.
!> \param zeros Zero primitive gradient potentials for the separate inverse core response.
!> \param zero Zero primitive density and tau potentials for that core response.
!> \param atom_response Physical atomic motion.
!> \param cell_response Physical lattice derivative.
!> \param reduced_k_response Fixed reduced k coordinates for force and stress.
!> \param response Core plus complete Coulomb and signed normalization response.
!> \param metric_response Signed complete Pulay contribution, already included in response.
! **************************************************************************************************
   SUBROUTINE complete_response(qs_env, ks_env, state, density, weighted, zeros, zero, atom_response, cell_response, &
                                reduced_k_response, response, metric_response)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_ks_env_type), INTENT(INOUT)                :: ks_env
      TYPE(qs_upf_inverse_scf_type), INTENT(IN), TARGET  :: state
      COMPLEX(dp), INTENT(IN)                            :: density(:, :, :), weighted(:, :, :)
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: zeros(:, :), zero(:)
      REAL(dp), INTENT(IN)                               :: atom_response(:, :), &
                                                            cell_response(3, 3), &
                                                            reduced_k_response(:, :)
      REAL(dp), INTENT(OUT)                              :: response, metric_response

      REAL(dp)                                           :: coulomb_response, radial_metric

      CPASSERT(ALL(reduced_k_response == 0.0_dp))
      CALL qs_upf_inverse_scf_response(ks_env, zero, zeros, zero, atom_response, cell_response, reduced_k_response, &
                                       response, radial_metric)
      CALL qs_upf_inverse_coulomb_metric_response(qs_env, state, state%coulomb, weighted, atom_response, cell_response, &
                                                  metric_response)
      metric_response = -metric_response
      response = response - radial_metric + metric_response
      CALL qs_upf_inverse_coulomb_response(qs_env, state, state%coulomb, density, atom_response, cell_response, coulomb_response)
      response = response + coulomb_response
   END SUBROUTINE complete_response

END MODULE qs_upf_inverse_forces
