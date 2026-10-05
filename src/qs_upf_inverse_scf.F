!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Explicit inverse-orbital attachment to native Quickstep density and operator consumers.
!>        Development entry for scalar semilocal GPW at Gamma, with original full-wave USPP data.
!>        K-point, GAPW, kinetic-augmentation and geometry-response consumers are extended jointly.
! **************************************************************************************************
MODULE qs_upf_inverse_scf
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
   USE input_constants,                 ONLY: sic_none
   USE kinds,                           ONLY: dp
   USE kpoint_transitional,             ONLY: set_1d_pointer,&
                                              set_2d_pointer
   USE kpoint_types,                    ONLY: kpoint_native_operator_clear
   USE mathlib,                         ONLY: invmat
   USE particle_types,                  ONLY: particle_type
   USE pw_grid_types,                   ONLY: pw_grid_type
   USE pw_methods,                      ONLY: pw_integrate_function,&
                                              pw_scale,&
                                              pw_transfer
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE qs_kind_types,                   ONLY: qs_kind_type
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
   USE qs_upf_inverse_augmentation,     ONLY: qs_upf_inverse_augmentation_prepare
   USE qs_upf_inverse_grid,             ONLY: qs_upf_inverse_grid_collocate,&
                                              qs_upf_inverse_grid_integrate,&
                                              qs_upf_inverse_grid_integrate_dense,&
                                              qs_upf_inverse_grid_metric,&
                                              qs_upf_inverse_grid_prepare,&
                                              qs_upf_inverse_grid_response,&
                                              qs_upf_inverse_matrix_gather
   USE qs_upf_inverse_kpoints,          ONLY: qs_upf_inverse_kpoints_collocate,&
                                              qs_upf_inverse_kpoints_density,&
                                              qs_upf_inverse_kpoints_integrate,&
                                              qs_upf_inverse_kpoints_metric,&
                                              qs_upf_inverse_kpoints_prepare,&
                                              qs_upf_inverse_kpoints_release_operators
   USE qs_upf_inverse_types,            ONLY: qs_upf_inverse_augmentation_type,&
                                              qs_upf_inverse_grid_type,&
                                              qs_upf_inverse_scf_release,&
                                              qs_upf_inverse_scf_type
   USE qs_uspp_types,                   ONLY: qs_uspp_type
   USE upf_local_integrals,             ONLY: upf_local_fourier,&
                                              upf_local_init
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: qs_upf_inverse_scf_prepare, qs_upf_inverse_scf_active
   PUBLIC :: qs_upf_inverse_scf_collocate, qs_upf_inverse_scf_integrate, qs_upf_inverse_scf_core
   PUBLIC :: qs_upf_inverse_scf_traces
   PUBLIC :: qs_upf_inverse_scf_core_response
CONTAINS

! **************************************************************************************************
!> \brief Explicitly install one inverse frame, its charge metric, and core operator in native KS.
!>        Never inferred from full-wave availability; caller supplies the fixed image set.
!> \param ks_env Prepared native KS environment.
!> \param grid Native physical-cell PW grid.
!> \param images Fixed physical-cell images for the Gaussian and radial fields.
!> \param info Factorization status; nonzero leaves the reference KS environment unchanged.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_scf_prepare(ks_env, grid, images, info)
      TYPE(qs_ks_env_type), INTENT(INOUT)                :: ks_env
      TYPE(pw_grid_type), POINTER                        :: grid
      INTEGER, INTENT(IN)                                :: images(:, :)
      INTEGER, INTENT(OUT)                               :: info

      INTEGER                                            :: ikind, img
      LOGICAL                                            :: do_kpoints
      TYPE(cp_cfm_type), ALLOCATABLE                     :: native_metric(:)
      TYPE(dbcsr_p_type), POINTER                        :: metric(:), metric_k(:, :), overlap(:), &
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
      IF (.NOT. control%qs_control%gpw .OR. control%qs_control%gapw .OR. control%qs_control%gapw_xc) THEN
         CPABORT("Inverse SCF GAPW reconstruction is not assembled")
      END IF
      IF (control%use_kinetic_energy_density) THEN
         CPABORT("Missing inverse kinetic augmentation for tau XC")
      END IF
      IF (control%drho_by_collocation) THEN
         CPABORT("Inverse SCF needs spectral XC gradients")
      END IF
      IF (control%do_admm .OR. control%do_sccs .OR. control%qs_control%cdft .OR. &
        control%qs_control%rigpw .OR. control%qs_control%lrigpw .OR. control%dft_plus_u .OR. control%sic_method_id /= sic_none) THEN
         CPABORT("Unsupported inverse SCF density consumer")
      END IF
      DO ikind = 1, SIZE(kinds)
         IF (kinds(ikind)%ghost) CYCLE
         CPASSERT(ASSOCIATED(kinds(ikind)%sgp_potential))
         CPASSERT(ALLOCATED(kinds(ikind)%sgp_potential%upf_data))
         IF (.NOT. kinds(ikind)%sgp_potential%upf_data%is_ultrasoft) THEN
            CPABORT("Missing mixed NC inverse nonlocal operator")
         END IF
      END DO
      ALLOCATE (state)
      IF (do_kpoints) THEN
         CALL qs_upf_inverse_kpoints_prepare(ks_env, grid, images, state%bloch, info)
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
         CALL register_bloch(state, metric_k, kpoint_slot_s, state%metric_k)
         CALL qs_upf_inverse_scf_core(ks_env)
         CALL kpoint_density_state_initialize(state%bloch%kpoints)
         CALL qs_ks_did_change(ks_env, s_mstruct_changed=.TRUE., rho_changed=.TRUE., potential_changed=.TRUE.)
         RETURN
      END IF
      CPASSERT(control%nimages == 1)
      CALL get_ks_env(ks_env, matrix_s=overlap)
      CALL qs_upf_inverse_grid_prepare(ks_env, grid, images, [0.0_dp, 0.0_dp, 0.0_dp], state%frame, info)
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
      CALL qs_upf_inverse_scf_core(ks_env)
      CALL qs_ks_did_change(ks_env, s_mstruct_changed=.TRUE., rho_changed=.TRUE., potential_changed=.TRUE.)
   END SUBROUTINE qs_upf_inverse_scf_prepare

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
!> \brief Build the original short-range local, Dion and soft kinetic operator in the inverse frame.
!>        The long-range ionic/core Coulomb term remains in native Poisson energy accounting.
!> \param ks_env ...
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_scf_core(ks_env)
      TYPE(qs_ks_env_type), INTENT(INOUT)                :: ks_env

      INTEGER                                            :: axis
      TYPE(dbcsr_p_type), POINTER                        :: h(:), kinetic(:), overlap(:)
      TYPE(dbcsr_type)                                   :: imaginary, real_part
      TYPE(pw_grid_type), POINTER                        :: grid
      TYPE(pw_r3d_rs_type)                               :: local_r, unit_tau, zero_density, zeros(3)
      TYPE(pw_r3d_rs_type), POINTER                      :: native_r(:)
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_upf_inverse_scf_type), POINTER             :: state

      CALL get_ks_env(ks_env, upf_inverse=state)
      CPASSERT(ASSOCIATED(state))
      IF (ASSOCIATED(state%bloch%kpoints)) THEN
         CALL core_bloch(ks_env, state)
         RETURN
      END IF
      CALL get_ks_env(ks_env, upf_inverse=state, &
                      matrix_h=h, kinetic=kinetic, matrix_s=overlap, rho=rho)
      CPASSERT(ASSOCIATED(state))
      CALL check_geometry(ks_env, state)
      CALL qs_rho_get(rho, rho_r=native_r)
      grid => native_r(1)%pw_grid
      CALL local_r%create(grid)
      CALL unit_tau%create(grid)
      CALL zero_density%create(grid)
      unit_tau%array = 1.0_dp
      zero_density%array = 0.0_dp
      DO axis = 1, 3
         CALL zeros(axis)%create(grid)
         zeros(axis)%array = 0.0_dp
      END DO
      CALL qs_upf_inverse_grid_integrate(state%frame, zero_density, zeros, unit_tau, overlap(1)%matrix, &
                                         real_part, imaginary)
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
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(qs_rho_type), INTENT(INOUT)                   :: rho_struct

      INTEGER                                            :: axis, spin
      REAL(dp), POINTER                                  :: total(:)
      TYPE(dbcsr_p_type), POINTER                        :: density(:, :)
      TYPE(pw_c1d_gs_type), POINTER                      :: rho_g(:)
      TYPE(pw_r3d_rs_type)                               :: gradient(3), tau
      TYPE(pw_r3d_rs_type), POINTER                      :: rho_r(:)
      TYPE(qs_upf_inverse_scf_type), POINTER             :: state

      CALL get_ks_env(ks_env, upf_inverse=state)
      CPASSERT(ASSOCIATED(state))
      CALL check_geometry(ks_env, state)
      CALL qs_rho_get(rho_struct, rho_ao_kp=density, rho_r=rho_r, rho_g=rho_g, tot_rho_r=total)
      IF (.NOT. ASSOCIATED(state%bloch%kpoints)) THEN
         CPASSERT(SIZE(density, 2) == 1)
      END IF
      CALL tau%create(rho_r(1)%pw_grid)
      DO axis = 1, 3
         CALL gradient(axis)%create(rho_r(1)%pw_grid)
      END DO
      DO spin = 1, SIZE(density, 1)
         IF (ASSOCIATED(state%bloch%kpoints)) THEN
            CALL qs_upf_inverse_kpoints_collocate(state%bloch, spin, rho_r(spin), gradient, tau)
         ELSE
            CALL qs_upf_inverse_grid_collocate(state%frame, density(spin, 1)%matrix, rho_r(spin), gradient, tau, &
                                               augmentation=state%augmentation, fourier_charge=.TRUE.)
         END IF
         CALL pw_transfer(rho_r(spin), rho_g(spin))
         total(spin) = pw_integrate_function(rho_r(spin), isign=-1)
      END DO
      CALL qs_rho_set(rho_struct, rho_r_valid=.TRUE., rho_g_valid=.TRUE.)
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
         CALL get_ks_env(ks_env, matrix_ks_kp=source)
         CALL register_bloch_spin(state, source)
      ELSE
         CALL qs_upf_inverse_grid_integrate(state%frame, unweighted, zeros, zero_tau, hamiltonian, &
                                            real_part, imaginary, augmentation=state%augmentation, fourier_charge=.TRUE.)
         CALL check_real(state, real_part, imaginary)
         CALL dbcsr_add(hamiltonian, real_part, 1.0_dp, 1.0_dp)
         CALL dbcsr_release(real_part)
         CALL dbcsr_release(imaginary)
      END IF
      CALL unweighted%release()
      CALL zero_tau%release()
      DO axis = 1, 3
         CALL zeros(axis)%release()
      END DO
      state%operator_evaluations = state%operator_evaluations + 1
   END SUBROUTINE qs_upf_inverse_scf_integrate

! **************************************************************************************************
!> \brief Joint geometry response of the inverse core Hamiltonian and its charge metric.
!>        Differentiates soft kinetic energy, original Dion, all short-range ionic potentials and Q.
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
      zero%array = 0.0_dp
      DO axis = 1, 3
         CALL zeros(axis)%create(grid)
         zeros(axis)%array = 0.0_dp
      END DO
      CALL local_potential(ks_env, grid, frame%ao%positions, local_r, local_response, atom_response, strain)
      CALL qs_upf_inverse_grid_response(frame, augmentation, frame%ao%radial_q, local_r, zeros, unit_tau, &
                                        atom_response, cell_response, reduced_k_response, local_response, zeros, zero, &
                                        core_response, metric_response, frame%ao%dion)
      CALL local_r%release()
      CALL local_response%release()
      CALL unit_tau%release()
      CALL zero%release()
      DO axis = 1, 3
         CALL zeros(axis)%release()
      END DO
   END SUBROUTINE qs_upf_inverse_scf_core_response

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
      TYPE(pw_grid_type), POINTER                        :: grid
      TYPE(pw_r3d_rs_type)                               :: local_r, unit_tau, zero, zeros(3)
      TYPE(pw_r3d_rs_type), POINTER                      :: native_r(:)
      TYPE(qs_rho_type), POINTER                         :: rho

      CALL check_geometry(ks_env, state)
      CALL get_ks_env(ks_env, rho=rho, kinetic_kp=kinetic, matrix_ks_kp=source)
      ! KS storage is recreated after the core step when the AO structure changes.
      ! Do not retain an operator provider tied to that released source identity.
      IF (.NOT. ASSOCIATED(source)) THEN
         CALL kpoint_native_operator_clear(state%bloch%kpoints%native_operators(kpoint_slot_ks))
      END IF
      CALL qs_rho_get(rho, rho_r=native_r)
      grid => native_r(1)%pw_grid
      CALL local_r%create(grid)
      CALL zero%create(grid)
      CALL unit_tau%create(grid)
      zero%array = 0.0_dp
      unit_tau%array = 1.0_dp
      DO axis = 1, 3
         CALL zeros(axis)%create(grid)
         zeros(axis)%array = 0.0_dp
      END DO
      CALL local_potential(ks_env, grid, state%bloch%frames(1)%ao%positions, local_r)
      n = state%bloch%nao
      nk = SIZE(state%bloch%weights)
      IF (.NOT. ALLOCATED(state%core_k)) THEN
         ALLOCATE (state%core_k(n, n, nk), state%kinetic_k(n, n, nk), state%ks_k(n, n, nk, SIZE(native_r)))
      END IF
      DO ik = 1, nk
         CALL qs_upf_inverse_grid_integrate_dense(state%bloch%frames(ik), zero, zeros, unit_tau, matrix)
         state%kinetic_k(:, :, ik) = matrix
         CALL qs_upf_inverse_grid_integrate_dense(state%bloch%frames(ik), local_r, zeros, unit_tau, matrix, &
                                                  projector_operator=state%bloch%frames(ik)%ao%dion, &
                                                  augmentation=state%bloch%augmentation, fourier_charge=.TRUE.)
         state%core_k(:, :, ik) = matrix
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

END MODULE qs_upf_inverse_scf
