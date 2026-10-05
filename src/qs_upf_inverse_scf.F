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
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_add,&
                                              dbcsr_copy,&
                                              dbcsr_p_type,&
                                              dbcsr_release,&
                                              dbcsr_set,&
                                              dbcsr_type
   USE input_constants,                 ONLY: sic_none
   USE kinds,                           ONLY: dp
   USE kpoint_transitional,             ONLY: set_1d_pointer
   USE particle_types,                  ONLY: particle_type
   USE pw_grid_types,                   ONLY: pw_grid_type
   USE pw_methods,                      ONLY: pw_integrate_function,&
                                              pw_scale,&
                                              pw_transfer
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE qs_kind_types,                   ONLY: qs_kind_type
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
                                              qs_upf_inverse_grid_metric,&
                                              qs_upf_inverse_grid_prepare,&
                                              qs_upf_inverse_matrix_gather
   USE qs_upf_inverse_types,            ONLY: qs_upf_inverse_scf_release,&
                                              qs_upf_inverse_scf_type
   USE qs_uspp_types,                   ONLY: qs_uspp_type
   USE upf_local_integrals,             ONLY: upf_local_fourier,&
                                              upf_local_init
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: qs_upf_inverse_scf_prepare, qs_upf_inverse_scf_active
   PUBLIC :: qs_upf_inverse_scf_collocate, qs_upf_inverse_scf_integrate, qs_upf_inverse_scf_core
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

      INTEGER                                            :: ikind
      LOGICAL                                            :: do_kpoints
      TYPE(dbcsr_p_type), POINTER                        :: metric(:), overlap(:)
      TYPE(dbcsr_type)                                   :: imaginary
      TYPE(dft_control_type), POINTER                    :: control
      TYPE(qs_kind_type), POINTER                        :: kinds(:)
      TYPE(qs_upf_inverse_scf_type), POINTER             :: previous, state
      TYPE(qs_uspp_type), POINTER                        :: uspp

      CALL get_ks_env(ks_env, upf_inverse=previous, matrix_s=overlap, uspp=uspp, &
                      dft_control=control, do_kpoints=do_kpoints, qs_kind_set=kinds)
      CPASSERT(.NOT. ASSOCIATED(previous))
      CPASSERT(ASSOCIATED(uspp) .AND. ASSOCIATED(grid))
      IF (do_kpoints .OR. control%nimages /= 1) THEN
         CPABORT("Inverse SCF k-point groups are not assembled")
      END IF
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
          control%qs_control%rigpw .OR. control%qs_control%lrigpw .OR. control%sic_method_id /= sic_none) THEN
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

      COMPLEX(dp), ALLOCATABLE                           :: ionic_g(:)
      INTEGER                                            :: atom, axis, ikind
      TYPE(dbcsr_p_type), POINTER                        :: h(:), kinetic(:), overlap(:)
      TYPE(dbcsr_type)                                   :: imaginary, real_part
      TYPE(particle_type), POINTER                       :: atoms(:)
      TYPE(pw_c1d_gs_type)                               :: local_g
      TYPE(pw_grid_type), POINTER                        :: grid
      TYPE(pw_r3d_rs_type)                               :: local_r, unit_tau, zero_density, zeros(3)
      TYPE(pw_r3d_rs_type), POINTER                      :: native_r(:)
      TYPE(qs_kind_type), POINTER                        :: kinds(:)
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_upf_inverse_scf_type), POINTER             :: state

      CALL get_ks_env(ks_env, upf_inverse=state, qs_kind_set=kinds, particle_set=atoms, &
                      matrix_h=h, kinetic=kinetic, matrix_s=overlap, rho=rho)
      CPASSERT(ASSOCIATED(state))
      CALL check_geometry(ks_env, state)
      CALL qs_rho_get(rho, rho_r=native_r)
      grid => native_r(1)%pw_grid
      CALL local_g%create(grid)
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
      ALLOCATE (ionic_g(SIZE(local_g%array)))
      local_g%array = CMPLX(0.0_dp, 0.0_dp, dp)
      DO atom = 1, SIZE(atoms)
         ikind = atoms(atom)%atomic_kind%kind_number
         IF (kinds(ikind)%ghost) CYCLE
         ASSOCIATE (potential => kinds(ikind)%sgp_potential)
            IF (potential%upf_local%radial%degree < 0) THEN
               CALL upf_local_init(potential%upf_data, potential%alpha_core_charge, 0, potential%upf_local)
            END IF
            CALL upf_local_fourier(potential%upf_local, grid%g, atoms(atom)%r, grid%vol, ionic_g)
         END ASSOCIATE
         local_g%array = local_g%array + ionic_g
      END DO
      CALL pw_transfer(local_g, local_r)
      CALL qs_upf_inverse_grid_integrate(state%frame, local_r, zeros, unit_tau, overlap(1)%matrix, &
                                         real_part, imaginary, projector_operator=state%frame%ao%dion, &
                                         augmentation=state%augmentation, fourier_charge=.TRUE.)
      CALL check_real(state, real_part, imaginary)
      CALL dbcsr_copy(h(1)%matrix, real_part)
      CALL dbcsr_release(real_part)
      CALL dbcsr_release(imaginary)
      CALL local_g%release()
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
      CPASSERT(SIZE(density, 2) == 1)
      CALL tau%create(rho_r(1)%pw_grid)
      DO axis = 1, 3
         CALL gradient(axis)%create(rho_r(1)%pw_grid)
      END DO
      DO spin = 1, SIZE(density, 1)
         CALL qs_upf_inverse_grid_collocate(state%frame, density(spin, 1)%matrix, rho_r(spin), gradient, tau, &
                                            augmentation=state%augmentation, fourier_charge=.TRUE.)
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
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_scf_integrate(ks_env, potential, hamiltonian, weighted)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: potential
      TYPE(dbcsr_type), INTENT(INOUT)                    :: hamiltonian
      LOGICAL, INTENT(IN)                                :: weighted

      INTEGER                                            :: axis
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
      CALL qs_upf_inverse_grid_integrate(state%frame, unweighted, zeros, zero_tau, hamiltonian, &
                                         real_part, imaginary, augmentation=state%augmentation, fourier_charge=.TRUE.)
      CALL check_real(state, real_part, imaginary)
      CALL dbcsr_add(hamiltonian, real_part, 1.0_dp, 1.0_dp)
      CALL dbcsr_release(real_part)
      CALL dbcsr_release(imaginary)
      CALL unweighted%release()
      CALL zero_tau%release()
      DO axis = 1, 3
         CALL zeros(axis)%release()
      END DO
      state%operator_evaluations = state%operator_evaluations + 1
   END SUBROUTINE qs_upf_inverse_scf_integrate

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
      CPASSERT(ALL(cell%hmat == state%frame%ao%cell))
      DO atom = 1, SIZE(atoms)
         CPASSERT(ALL(atoms(atom)%r == state%frame%ao%positions(:, atom)))
      END DO
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
END MODULE qs_upf_inverse_scf
