!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
!> \brief Independent original-UPF basis representations and their native Gaussian references.
MODULE qs_upf_inverse_basis
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type
   USE cp_dbcsr_operations,             ONLY: dbcsr_deallocate_matrix_set
   USE kinds,                           ONLY: dp
   USE pw_grid_types,                   ONLY: pw_grid_type
   USE qs_environment_types,            ONLY: qs_environment_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_kinetic,                      ONLY: build_kinetic_matrix
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type
   USE qs_upf_inverse_augmentation,     ONLY: qs_upf_inverse_augmentation_prepare
   USE qs_upf_inverse_background,       ONLY: qs_upf_inverse_background_operator,&
                                              qs_upf_inverse_background_prepare
   USE qs_upf_inverse_grid,             ONLY: qs_upf_inverse_grid_prepare
   USE qs_upf_inverse_gth_local,        ONLY: qs_upf_inverse_gth_local_active
   USE qs_upf_inverse_kpoints,          ONLY: qs_upf_inverse_kpoints_prepare
   USE qs_upf_inverse_scf_xc,           ONLY: qs_upf_inverse_scf_xc_reference
   USE qs_upf_inverse_types,            ONLY: qs_upf_inverse_scf_release,&
                                              qs_upf_inverse_scf_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: qs_upf_inverse_basis_prepare

CONTAINS

! **************************************************************************************************
!> \brief Prepare an owned orbital or auxiliary basis without installing a native SCF state.
!>        Both bases use the original projectors, augmentation and selected XC reference data.
!>        The caller owns this context and must release it before its borrowed bases and grids.
!> \param ks_env Native geometry and basis owner. Density and operator caches remain unchanged.
!> \param grid Physical-cell source grid borrowed for the lifetime of this context.
!> \param images Fixed physical-cell images for Gaussian and radial fields.
!> \param state New owned representation. It remains disabled as a native SCF owner.
!> \param info Inverse factorization status. Failure releases the incomplete context.
!> \param atom_local Use independent atom-local inverse projectors.
!> \param direct Use original soft orbitals without reconstruction.
!> \param basis_type ORB by default, or AUX_FIT with its actual native auxiliary basis.
!> \param qs_env Matching native owner, required for auxiliary Gaussian references and task lists.
!> \param kinetic Override the primary functional's kinetic-density data requirement.
!> \param composite_paw Override its joint learned PAW valence reconstruction requirement.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_basis_prepare(ks_env, grid, images, state, info, atom_local, direct, &
                                           basis_type, qs_env, kinetic, composite_paw)
      TYPE(qs_ks_env_type), INTENT(IN), TARGET           :: ks_env
      TYPE(pw_grid_type), INTENT(IN), POINTER            :: grid
      INTEGER, INTENT(IN)                                :: images(:, :)
      TYPE(qs_upf_inverse_scf_type), POINTER             :: state
      INTEGER, INTENT(OUT)                               :: info
      LOGICAL, INTENT(IN), OPTIONAL                      :: atom_local, direct
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL             :: basis_type
      TYPE(qs_environment_type), OPTIONAL, POINTER       :: qs_env
      LOGICAL, INTENT(IN), OPTIONAL                      :: kinetic, composite_paw

      CHARACTER(LEN=16)                                  :: role
      INTEGER                                            :: ikind
      LOGICAL                                            :: auxiliary, do_kpoints, use_atomic, &
                                                            use_composite, use_kinetic
      TYPE(dbcsr_p_type), POINTER                        :: kinetic_k(:, :)
      TYPE(dft_control_type), POINTER                    :: control
      TYPE(qs_kind_type), POINTER                        :: kinds(:)
      TYPE(qs_ks_env_type), POINTER                      :: native_ks

      CPASSERT(.NOT. ASSOCIATED(state))
      CPASSERT(ASSOCIATED(grid))
      role = 'ORB'
      IF (PRESENT(basis_type)) role = basis_type
      CPASSERT(role == 'ORB' .OR. role == 'AUX_FIT')
      auxiliary = role == 'AUX_FIT'
      CALL get_ks_env(ks_env, dft_control=control, qs_kind_set=kinds, do_kpoints=do_kpoints)
      DO ikind = 1, SIZE(kinds)
         IF (kinds(ikind)%ghost) CYCLE
         IF (ASSOCIATED(kinds(ikind)%gth_potential)) THEN
            IF (kinds(ikind)%gth_potential%soc .OR. kinds(ikind)%gth_potential%lsdpot) THEN
               CPABORT("Mixed GTH needs scalar spin-independent data")
            END IF
            CYCLE
         END IF
         CPASSERT(ASSOCIATED(kinds(ikind)%sgp_potential))
         CPASSERT(ALLOCATED(kinds(ikind)%sgp_potential%upf_data))
         IF (.NOT. kinds(ikind)%sgp_potential%upf_data%is_ultrasoft .AND. &
             .NOT. kinds(ikind)%sgp_potential%upf_data%is_paw .AND. &
             TRIM(kinds(ikind)%sgp_potential%upf_data%pseudo_type) /= "NC") THEN
            CPABORT("Unsupported inverse UPF type")
         END IF
      END DO
      ALLOCATE (state)
      state%basis_type = role
      state%source_grid => grid
      state%native_gth_local = qs_upf_inverse_gth_local_active(kinds)
      use_kinetic = control%use_kinetic_energy_density
      IF (PRESENT(kinetic)) use_kinetic = kinetic
      use_composite = control%use_gauxc
      IF (PRESENT(composite_paw)) use_composite = composite_paw
      CALL qs_upf_inverse_scf_xc_reference(kinds, use_kinetic, use_composite, state)
      use_atomic = .FALSE.
      IF (PRESENT(atom_local)) use_atomic = atom_local
      IF (PRESENT(direct)) THEN
         IF (direct) use_atomic = .TRUE.
         state%native_gaussian = direct .AND. .NOT. control%qs_control%gapw
      END IF
      ! Auxiliary XC needs the same soft Gaussian source as full GAPW, without owning Hartree data.
      state%gaussian_background = auxiliary .AND. control%qs_control%gapw
      IF (use_atomic .OR. state%native_gth_local .OR. state%gaussian_background) THEN
         CALL qs_upf_inverse_background_prepare(ks_env, state%gaussian, soft_basis=control%qs_control%gapw, &
                                                qs_env=qs_env, basis_type=role)
      END IF
      IF (use_atomic) THEN
         IF (auxiliary) THEN
            native_ks => ks_env
            NULLIFY (kinetic_k)
            CALL build_kinetic_matrix(native_ks, matrixkp_t=kinetic_k, basis_type=role, &
                                      sab_nl=state%gaussian%neighbors)
         ELSE
            CALL get_ks_env(ks_env, kinetic_kp=kinetic_k)
         END IF
         CPASSERT(ASSOCIATED(kinetic_k))
         ALLOCATE (state%gaussian_kinetic(state%gaussian%nao, state%gaussian%nao, SIZE(state%gaussian%weights)))
         CALL qs_upf_inverse_background_operator(state%gaussian, kinetic_k(1, :), state%gaussian_kinetic)
         IF (auxiliary) CALL dbcsr_deallocate_matrix_set(kinetic_k)
      END IF
      IF (do_kpoints) THEN
         IF (state%native_gaussian) THEN
            CALL qs_upf_inverse_kpoints_prepare(ks_env, grid, images, state%bloch, info, use_atomic, direct, &
                                                store_fields=.FALSE., basis_type=role)
         ELSE
            CALL qs_upf_inverse_kpoints_prepare(ks_env, grid, images, state%bloch, info, use_atomic, direct, basis_type=role)
         END IF
      ELSE
         IF (state%native_gaussian) THEN
            CALL qs_upf_inverse_grid_prepare(ks_env, grid, images, [0.0_dp, 0.0_dp, 0.0_dp], state%frame, info, &
                                             use_atomic, direct, store_fields=.FALSE., basis_type=role)
         ELSE
            CALL qs_upf_inverse_grid_prepare(ks_env, grid, images, [0.0_dp, 0.0_dp, 0.0_dp], state%frame, info, &
                                             use_atomic, direct, basis_type=role)
         END IF
         IF (info == 0) CALL qs_upf_inverse_augmentation_prepare(ks_env, state%frame%ao, state%augmentation)
      END IF
      IF (info /= 0) CALL qs_upf_inverse_scf_release(state)
   END SUBROUTINE qs_upf_inverse_basis_prepare

END MODULE qs_upf_inverse_basis
