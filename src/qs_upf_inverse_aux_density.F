!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
!> \brief Native auxiliary valence density and its complete original-UPF charge transpose.
MODULE qs_upf_inverse_aux_density
   USE admm_types,                      ONLY: admm_type,&
                                              admm_upf_center_density_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: fourpi
   USE pw_methods,                      ONLY: pw_integrate_function,&
                                              pw_transfer
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_grid_atom,                    ONLY: grid_atom_type
   USE qs_harmonics_atom,               ONLY: harmonics_atom_type
   USE qs_kind_types,                   ONLY: get_qs_kind
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_set
   USE qs_upf_inverse_density,          ONLY: qs_upf_inverse_density_gapw,&
                                              qs_upf_inverse_density_gapw_integrate
   USE qs_upf_inverse_grid,             ONLY: qs_upf_inverse_matrix_gather
   USE qs_upf_inverse_kpoints,          ONLY: qs_upf_inverse_kpoints_density
   USE qs_upf_inverse_scf_xc,           ONLY: qs_upf_inverse_scf_xc_source_adjoint,&
                                              qs_upf_inverse_scf_xc_source_direction
   USE qs_upf_inverse_types,            ONLY: qs_upf_inverse_augmentation_type,&
                                              qs_upf_inverse_grid_type,&
                                              qs_upf_inverse_scf_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: qs_upf_inverse_aux_density_build, qs_upf_inverse_aux_density_integrate

CONTAINS

! **************************************************************************************************
!> \brief Build auxiliary grid and signed center sources from the same native occupation matrices.
!>        GPW/GAPW_XC retain the original inverse and Fourier-Q grid sources. Full GAPW retains
!>        the actual native soft Gaussian grid plus native Gaussian hard-minus-soft, inverse and
!>        original-Q center densities on their own radial/harmonic grids. Frozen cores are excluded.
!>        The stored matrices retain all physical Bloch contributions, not just Gaussian images.
!> \param qs_env Prepared native auxiliary density and independent original-UPF representation.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_aux_density_build(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      COMPLEX(dp), ALLOCATABLE                           :: gamma(:, :), matrix(:, :, :)
      INTEGER                                            :: center, ik, kind, n, nc, nh, nk, ns, spin
      LOGICAL                                            :: valid
      REAL(dp), ALLOCATABLE                              :: radial(:, :)
      REAL(dp), POINTER                                  :: grid_charge(:)
      TYPE(admm_type), POINTER                           :: admm
      TYPE(dbcsr_p_type), POINTER                        :: density_ao(:)
      TYPE(grid_atom_type), POINTER                      :: grid
      TYPE(harmonics_atom_type), POINTER                 :: harmonics
      TYPE(pw_c1d_gs_type), POINTER                      :: density_g(:), kinetic_g(:)
      TYPE(pw_r3d_rs_type), ALLOCATABLE                  :: temporary_kinetic(:)
      TYPE(pw_r3d_rs_type), POINTER                      :: density_r(:), kinetic_r(:)
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_upf_inverse_augmentation_type), POINTER    :: augmentation
      TYPE(qs_upf_inverse_grid_type), POINTER            :: frame
      TYPE(qs_upf_inverse_scf_type), POINTER             :: state

      CALL get_qs_env(qs_env, ks_env=ks_env, admm_env=admm)
      CPASSERT(ASSOCIATED(admm))
      CPASSERT(admm%upf_aux_metric%ready)
      CPASSERT(.NOT. admm%upf_aux_metric%ordered)
      state => admm%upf_inverse
      CPASSERT(ASSOCIATED(state))
      CALL qs_rho_get(admm%rho_aux_fit, rho_r=density_r, rho_g=density_g, tot_rho_r=grid_charge, &
                      tau_r=kinetic_r, tau_g=kinetic_g)
      CPASSERT(ASSOCIATED(density_r) .AND. ASSOCIATED(density_g))
      CPASSERT(ASSOCIATED(grid_charge))
      admm%upf_density%valid = .FALSE.
      ns = SIZE(density_r)
      n = SIZE(admm%upf_aux_metric%values, 1)
      nk = SIZE(admm%upf_aux_metric%weights)
      IF (.NOT. ALLOCATED(admm%upf_density%matrices)) THEN
         ALLOCATE (admm%upf_density%matrices(n, n, nk, ns), admm%upf_density%valence_charge(ns))
      END IF
      valid = ALL(SHAPE(admm%upf_density%matrices) == [n, n, nk, ns])
      CPASSERT(valid)
      CPASSERT(SIZE(admm%upf_density%valence_charge) == ns)
      IF (ASSOCIATED(state%bloch%kpoints)) THEN
         DO spin = 1, ns
            CALL qs_upf_inverse_kpoints_density(state%bloch, spin, matrix)
            admm%upf_density%matrices(:, :, :, spin) = matrix
         END DO
         augmentation => state%bloch%augmentation
      ELSE
         CPASSERT(nk == 1)
         CALL qs_rho_get(admm%rho_aux_fit, rho_ao=density_ao)
         DO spin = 1, ns
            CALL qs_upf_inverse_matrix_gather(state%frame%ao, density_ao(spin)%matrix, gamma)
            admm%upf_density%matrices(:, :, 1, spin) = gamma
         END DO
         augmentation => state%augmentation
      END IF

      IF (ALLOCATED(state%kinetic_cache) .AND. .NOT. ASSOCIATED(kinetic_r)) THEN
         CPASSERT(.NOT. ASSOCIATED(kinetic_g))
         ALLOCATE (kinetic_r(ns), kinetic_g(ns))
         DO spin = 1, ns
            CALL kinetic_r(spin)%create(density_r(spin)%pw_grid)
            CALL kinetic_g(spin)%create(density_g(spin)%pw_grid)
         END DO
         CALL qs_rho_set(admm%rho_aux_fit, tau_r=kinetic_r, tau_g=kinetic_g)
      END IF
      IF (ASSOCIATED(kinetic_r)) THEN
         CPASSERT(ASSOCIATED(kinetic_g))
         CALL qs_upf_inverse_scf_xc_source_direction(ks_env, admm%upf_density%matrices, density_r, kinetic_r, &
                                                     context=state)
         DO spin = 1, ns
            CALL pw_transfer(kinetic_r(spin), kinetic_g(spin))
         END DO
         CALL qs_rho_set(admm%rho_aux_fit, tau_r_valid=.TRUE., tau_g_valid=.TRUE.)
      ELSE
         ALLOCATE (temporary_kinetic(ns))
         DO spin = 1, ns
            CALL temporary_kinetic(spin)%create(density_r(spin)%pw_grid)
         END DO
         CALL qs_upf_inverse_scf_xc_source_direction(ks_env, admm%upf_density%matrices, density_r, temporary_kinetic, &
                                                     context=state)
         DO spin = 1, ns
            CALL temporary_kinetic(spin)%release()
         END DO
      END IF
      DO spin = 1, ns
         CALL pw_transfer(density_r(spin), density_g(spin))
         grid_charge(spin) = pw_integrate_function(density_r(spin))
         admm%upf_density%valence_charge(spin) = grid_charge(spin)
      END DO

      nc = 0
      IF (admm%upf_aux_metric%atomic) nc = SIZE(admm%upf_aux_metric%center_atoms)
      IF (.NOT. ALLOCATED(admm%upf_density%centers)) ALLOCATE (admm%upf_density%centers(nc))
      CPASSERT(SIZE(admm%upf_density%centers) == nc)
      DO center = 1, nc
         kind = admm%upf_aux_metric%center_kinds(center)
         CALL get_qs_kind(admm%upf_aux_metric%row_kinds(kind), grid_atom=grid, harmonics=harmonics)
         nh = harmonics%max_iso_not0
         IF (.NOT. ALLOCATED(admm%upf_density%centers(center)%values)) THEN
            ALLOCATE (admm%upf_density%centers(center)%values(grid%nr, nh, ns))
         END IF
         valid = ALL(SHAPE(admm%upf_density%centers(center)%values) == [grid%nr, nh, ns])
         CPASSERT(valid)
         ALLOCATE (radial(grid%nr, nh))
         admm%upf_density%centers(center)%values = 0.0_dp
         DO spin = 1, ns
            DO ik = 1, nk
               frame => auxiliary_frame(state, ik)
               CALL qs_upf_inverse_density_gapw(frame%ao, augmentation, admm%upf_aux_metric%inverse_centers(center), &
                                                admm%upf_aux_metric%row_kinds(kind), admm%upf_aux_metric%pairs(center), &
                                                admm%upf_density%matrices(:, :, ik, spin), radial, &
                                                atom_index=admm%upf_aux_metric%center_atoms(center))
               admm%upf_density%centers(center)%values(:, :, spin) = admm%upf_density%centers(center)%values(:, :, spin) + &
                                                                     admm%upf_aux_metric%weights(ik)*radial
            END DO
            admm%upf_density%valence_charge(spin) = admm%upf_density%valence_charge(spin) + &
                                                    SQRT(fourpi)*SUM(grid%wr*admm%upf_density%centers(center)%values(:, 1, spin))
         END DO
         DEALLOCATE (radial)
      END DO
      CALL qs_rho_set(admm%rho_aux_fit, rho_r_valid=.TRUE., rho_g_valid=.TRUE.)
      state%density_evaluations = state%density_evaluations + 1
      state%orbital_factors%density_evaluation = -1
      admm%upf_density%valid = .TRUE.
   END SUBROUTINE qs_upf_inverse_aux_density_build

! **************************************************************************************************
!> \brief Complete AO transpose of the auxiliary grid and signed native-center density sources.
!>        Grid and center covectors carry no mesh weights. Their quadrature measures enter once
!>        through the existing native transposes. Returned Hermitian AO operators carry no k weight.
!>        Atomic kinetic-density and PAW XC operators are separate functional contributions.
!> \param qs_env Same native auxiliary geometry and source owners as the density evaluation.
!> \param potential Grid density covectors on their actual XC grid.
!> \param kinetic_potential Matching grid kinetic-density covectors, zero without kinetic sources.
!> \param operators New globally reduced operators (AO,AO,physical k,spin).
!> \param centers Native radial/harmonic covectors, required for full GAPW.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_aux_density_integrate(qs_env, potential, kinetic_potential, operators, centers)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: potential(:), kinetic_potential(:)
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: operators(:, :, :, :)
      TYPE(admm_upf_center_density_type), INTENT(IN), &
         OPTIONAL                                        :: centers(:)

      COMPLEX(dp), ALLOCATABLE                           :: partial(:, :)
      INTEGER                                            :: center, ik, kind, n, ns, spin
      LOGICAL                                            :: valid
      TYPE(admm_type), POINTER                           :: admm
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_upf_inverse_augmentation_type), POINTER    :: augmentation
      TYPE(qs_upf_inverse_grid_type), POINTER            :: frame
      TYPE(qs_upf_inverse_scf_type), POINTER             :: state

      CALL get_qs_env(qs_env, ks_env=ks_env, admm_env=admm)
      CPASSERT(ASSOCIATED(admm))
      CPASSERT(admm%upf_density%valid)
      CPASSERT(admm%upf_aux_metric%ready)
      state => admm%upf_inverse
      CPASSERT(ASSOCIATED(state))
      ns = SIZE(potential)
      CALL qs_upf_inverse_scf_xc_source_adjoint(ks_env, potential, kinetic_potential, operators, context=state)
      IF (.NOT. admm%upf_aux_metric%atomic) RETURN
      CPASSERT(PRESENT(centers))
      CPASSERT(SIZE(centers) == SIZE(admm%upf_density%centers))
      n = SIZE(operators, 1)
      ALLOCATE (partial(n, n))
      augmentation => state%augmentation
      IF (ASSOCIATED(state%bloch%kpoints)) augmentation => state%bloch%augmentation
      DO center = 1, SIZE(centers)
         valid = ALL(SHAPE(centers(center)%values) == SHAPE(admm%upf_density%centers(center)%values))
         CPASSERT(valid)
         kind = admm%upf_aux_metric%center_kinds(center)
         DO spin = 1, ns
            DO ik = 1, SIZE(operators, 3)
               frame => auxiliary_frame(state, ik)
               CALL qs_upf_inverse_density_gapw_integrate(frame%ao, augmentation, &
                                                          admm%upf_aux_metric%inverse_centers(center), &
                                                          admm%upf_aux_metric%row_kinds(kind), &
                                                          admm%upf_aux_metric%pairs(center), &
                                                          centers(center)%values(:, :, spin), partial, &
                                                          atom_index=admm%upf_aux_metric%center_atoms(center))
               operators(:, :, ik, spin) = operators(:, :, ik, spin) + partial
            END DO
         END DO
      END DO
   END SUBROUTINE qs_upf_inverse_aux_density_integrate

! **************************************************************************************************
!> \brief Borrow the auxiliary frame for one physical Bloch contribution.
!> \param state Prepared auxiliary representation, kept alive by its native ADMM owner.
!> \param ik Physical Bloch index, unity at Gamma without an explicit mesh.
!> \return Original-coordinate inverse frame.
! **************************************************************************************************
   FUNCTION auxiliary_frame(state, ik) RESULT(frame)
      TYPE(qs_upf_inverse_scf_type), INTENT(IN), TARGET  :: state
      INTEGER, INTENT(IN)                                :: ik
      TYPE(qs_upf_inverse_grid_type), POINTER            :: frame

      frame => state%frame
      IF (ASSOCIATED(state%bloch%kpoints)) THEN
         CPASSERT(ik >= 1 .AND. ik <= SIZE(state%bloch%frames))
         frame => state%bloch%frames(ik)
      ELSE
         CPASSERT(ik == 1)
      END IF
   END FUNCTION auxiliary_frame
END MODULE qs_upf_inverse_aux_density
