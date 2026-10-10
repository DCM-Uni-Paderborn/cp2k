!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Owned inverse USPP data, independent of the KS methods that consume it.
!>        Bases remain borrowed. Explicit SCF activation requires a prepared common frame.
! **************************************************************************************************
MODULE qs_upf_inverse_types
   USE basis_set_types,                 ONLY: gto_basis_set_p_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_deallocate_matrix,&
                                              dbcsr_p_type
   USE kinds,                           ONLY: dp
   USE kpoint_transitional,             ONLY: kpoint_transitional_release,&
                                              kpoint_transitional_type
   USE kpoint_types,                    ONLY: kp_transform_plan_type,&
                                              kpoint_clear_scf_smat,&
                                              kpoint_native_density_clear,&
                                              kpoint_native_operator_clear,&
                                              kpoint_type
   USE message_passing,                 ONLY: mp_comm_type
   USE paw_types,                       ONLY: paw_kind_type,&
                                              paw_potential_type
   USE pw_grid_types,                   ONLY: PW_MODE_LOCAL,&
                                              pw_grid_type
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE qs_upf_inverse_atom_grid,        ONLY: qs_upf_inverse_atom_grid_type
   USE qs_upf_inverse_gaussian_pairs,   ONLY: qs_upf_inverse_gaussian_pairs_type
   USE upf_inverse,                     ONLY: upf_inverse_atomic_type,&
                                              upf_inverse_kinetic_type,&
                                              upf_inverse_type
   USE upf_inverse_coulomb,             ONLY: upf_inverse_coulomb_center_type,&
                                              upf_inverse_coulomb_map_type
   USE upf_inverse_lattice,             ONLY: upf_inverse_lattice_type
   USE upf_radial_fields,               ONLY: upf_radial_fields_type
   USE uspp_augmentation,               ONLY: uspp_augmentation_type

   IMPLICIT NONE
   PRIVATE

   TYPE qs_upf_inverse_type
      TYPE(upf_inverse_type) :: inverse
      LOGICAL :: atom_local = .FALSE., direct = .FALSE.
      TYPE(upf_inverse_atomic_type), ALLOCATABLE :: atomic_cache(:)
      TYPE(mp_comm_type) :: point_comm
      ! Bases are borrowed from the live Quickstep environment; radial caches and fields are owned.
      TYPE(gto_basis_set_p_type), ALLOCATABLE :: bases(:)
      TYPE(upf_radial_fields_type), ALLOCATABLE :: beta_cache(:), delta_cache(:)
      INTEGER, ALLOCATABLE :: ao_first(:), ao_last(:), projector_first(:), projector_last(:), center_atoms(:), center_kinds(:)
      INTEGER, ALLOCATABLE :: images(:, :)
      INTEGER :: periodicity(3) = 0
      REAL(dp) :: cell(3, 3) = 0.0_dp, reduced_k(3) = 0.0_dp
      REAL(dp), ALLOCATABLE :: positions(:, :), points(:, :), weights(:)
      ! Supplied PP_Q and the reference metric's radial Q integral are separate data.
      COMPLEX(dp), ALLOCATABLE :: coefficients(:, :), q(:, :), radial_q(:, :), dion(:, :), beta(:, :), delta(:, :), &
                                  delta_gradient(:, :, :), gaussian(:, :), gaussian_gradient(:, :, :), &
                                  soft(:, :), soft_gradient(:, :, :)
   END TYPE qs_upf_inverse_type

   TYPE qs_upf_inverse_lattice_type
      LOGICAL :: ready = .FALSE.
      TYPE(upf_inverse_lattice_type) :: inverse
      ! B keeps its complete defining support, independently of the retained W domain.
      INTEGER, ALLOCATABLE :: overlap_images(:, :)
      COMPLEX(dp), ALLOCATABLE :: overlaps(:, :, :), coefficients(:, :, :)
      REAL(dp) :: tolerance = 0.0_dp, solve_tolerance = 0.0_dp, lower_bound = 0.0_dp
      REAL(dp) :: coefficient_error = HUGE(0.0_dp), bloch_error = HUGE(0.0_dp)
   END TYPE qs_upf_inverse_lattice_type

   TYPE qs_upf_inverse_grid_type
      TYPE(qs_upf_inverse_type) :: ao
      INTEGER :: bounds(2, 3) = 0, bounds_local(2, 3) = 0, npts(3) = 0, mode = PW_MODE_LOCAL
      REAL(dp) :: dh(3, 3) = 0.0_dp, dvol = 0.0_dp
   END TYPE qs_upf_inverse_grid_type

   TYPE qs_upf_inverse_augmentation_type
      TYPE(uspp_augmentation_type), ALLOCATABLE :: kinds(:)
   END TYPE qs_upf_inverse_augmentation_type

   TYPE qs_upf_inverse_kpoints_type
      ! The live initialized group environment is borrowed; all fields below are owned.
      TYPE(kpoint_type), POINTER :: kpoints => NULL()
      TYPE(qs_upf_inverse_grid_type), ALLOCATABLE :: frames(:)
      TYPE(qs_upf_inverse_augmentation_type) :: augmentation
      REAL(dp), ALLOCATABLE :: reduced_k(:, :), weights(:)
      INTEGER :: nao = 0, kp_range(2) = 0
   END TYPE qs_upf_inverse_kpoints_type

   TYPE qs_upf_inverse_xc_type
      LOGICAL :: kinetic = .FALSE.
      LOGICAL :: paw_reconstruction = .FALSE.
      LOGICAL :: weighted_adjoint = .FALSE.
      COMPLEX(dp), ALLOCATABLE :: soft(:, :, :), gradient(:, :, :, :), operators(:, :, :, :)
      REAL(dp), ALLOCATABLE :: density(:, :), density_gradient(:, :, :), energy_density(:), &
                               vdensity(:, :), vgradient(:, :, :), background_vdensity(:, :), &
                               background_vgradient(:, :, :), tau(:, :), vtau(:, :), background_vtau(:, :)
      REAL(dp), ALLOCATABLE :: point_derivative(:, :), weight_derivative(:), atom_derivative(:, :)
   END TYPE qs_upf_inverse_xc_type

   TYPE qs_upf_inverse_background_type
      INTEGER :: nao = 0
      LOGICAL :: bloch = .FALSE.
      TYPE(mp_comm_type) :: comm
      TYPE(kp_transform_plan_type) :: plan
      TYPE(dbcsr_p_type), POINTER :: density(:) => NULL(), potential(:) => NULL()
      REAL(dp), ALLOCATABLE :: reduced_k(:, :), weights(:)
   END TYPE qs_upf_inverse_background_type

   TYPE qs_upf_inverse_coulomb_type
      LOGICAL :: ready = .FALSE.
      TYPE(qs_upf_inverse_background_type) :: background
      TYPE(qs_upf_inverse_gaussian_pairs_type), ALLOCATABLE :: pairs(:)
      TYPE(upf_inverse_coulomb_center_type), ALLOCATABLE :: centers(:)
      INTEGER, ALLOCATABLE :: center_atoms(:), center_kinds(:), inverse_centers(:)
      TYPE(upf_inverse_coulomb_map_type) :: maps
      TYPE(pw_grid_type), POINTER :: grid => NULL()
      TYPE(pw_r3d_rs_type) :: source, potential
      TYPE(pw_c1d_gs_type) :: source_g, potential_g, ionic
      REAL(dp), ALLOCATABLE :: weights(:), multiplicity(:), poisson_weights(:), positions(:, :)
      REAL(dp) :: cell(3, 3) = 0.0_dp
      COMPLEX(dp), ALLOCATABLE :: metric(:, :, :)
   END TYPE qs_upf_inverse_coulomb_type

   TYPE qs_upf_inverse_scf_type
      LOGICAL :: enabled = .FALSE.
      LOGICAL :: gaussian_background = .FALSE.
      INTEGER :: density_evaluations = 0, operator_evaluations = 0
      TYPE(qs_upf_inverse_grid_type) :: frame
      TYPE(qs_upf_inverse_lattice_type) :: lattice
      TYPE(qs_upf_inverse_augmentation_type) :: augmentation
      TYPE(kpoint_transitional_type) :: metric
      TYPE(upf_inverse_kinetic_type), ALLOCATABLE :: kinetic_cache(:)
      ! PAW AE valence fields for one joint learned functional, separate from USPP kinetic kernels.
      TYPE(upf_inverse_kinetic_type), ALLOCATABLE :: paw_reconstruction(:)
      TYPE(paw_kind_type), ALLOCATABLE :: paw(:)
      TYPE(paw_potential_type), ALLOCATABLE :: paw_potential(:)
      COMPLEX(dp), ALLOCATABLE :: paw_operator(:, :, :, :), paw_density(:, :, :, :)
      INTEGER :: paw_density_evaluation = -1
      LOGICAL :: has_core_tau = .FALSE.
      TYPE(pw_r3d_rs_type) :: core_tau_r
      TYPE(pw_c1d_gs_type) :: core_tau_g
      TYPE(qs_upf_inverse_kpoints_type) :: bloch
      TYPE(qs_upf_inverse_background_type) :: gaussian
      COMPLEX(dp), ALLOCATABLE :: gaussian_kinetic(:, :, :)
      COMPLEX(dp), ALLOCATABLE :: metric_k(:, :, :), core_k(:, :, :), kinetic_k(:, :, :), ks_k(:, :, :, :)
      TYPE(qs_upf_inverse_atom_grid_type) :: atom_grid
      TYPE(qs_upf_inverse_xc_type) :: xc
      TYPE(qs_upf_inverse_coulomb_type), POINTER :: coulomb => NULL()
      COMPLEX(dp), ALLOCATABLE :: hartree_k(:, :, :)
      INTEGER :: coulomb_density_evaluation = -1, coulomb_operator_evaluation = -1
      TYPE(pw_grid_type), POINTER :: xc_grid => NULL()
      TYPE(upf_radial_fields_type), ALLOCATABLE :: core_cache(:), core_kinetic_cache(:)
      INTEGER, ALLOCATABLE :: core_images(:, :)
      REAL(dp), ALLOCATABLE :: xc_fields(:, :, :, :), xc_kweights(:)
      COMPLEX(dp), ALLOCATABLE :: xc_matrices(:, :, :, :)
   END TYPE qs_upf_inverse_scf_type

   PUBLIC :: qs_upf_inverse_type, qs_upf_inverse_grid_type, qs_upf_inverse_augmentation_type
   PUBLIC :: qs_upf_inverse_lattice_type
   PUBLIC :: qs_upf_inverse_kpoints_type, qs_upf_inverse_scf_type, qs_upf_inverse_scf_release
   PUBLIC :: qs_upf_inverse_xc_type
   PUBLIC :: qs_upf_inverse_background_type, qs_upf_inverse_background_release
   PUBLIC :: qs_upf_inverse_coulomb_type, qs_upf_inverse_coulomb_release

CONTAINS

! **************************************************************************************************
!> \brief Release the owned metric and fields before their borrowed KS basis owner is destroyed.
!> \param state Owned SCF state; pointer is nullified.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_scf_release(state)
      TYPE(qs_upf_inverse_scf_type), POINTER             :: state

      INTEGER                                            :: slot

      IF (.NOT. ASSOCIATED(state)) RETURN
      CALL qs_upf_inverse_background_release(state%gaussian)
      IF (state%has_core_tau) THEN
         CALL state%core_tau_r%release()
         CALL state%core_tau_g%release()
      END IF
      IF (ASSOCIATED(state%coulomb)) THEN
         CALL qs_upf_inverse_coulomb_release(state%coulomb)
         DEALLOCATE (state%coulomb)
      END IF
      IF (state%enabled .AND. ASSOCIATED(state%bloch%kpoints)) THEN
         DO slot = 1, SIZE(state%bloch%kpoints%native_operators)
            CALL kpoint_native_operator_clear(state%bloch%kpoints%native_operators(slot))
         END DO
         CALL kpoint_native_density_clear(state%bloch%kpoints)
         CALL kpoint_clear_scf_smat(state%bloch%kpoints)
      END IF
      CALL kpoint_transitional_release(state%metric)
      DEALLOCATE (state)
      NULLIFY (state)
   END SUBROUTINE qs_upf_inverse_scf_release

! **************************************************************************************************
!> \brief Release owned image matrices without changing the native density or operator state.
!> \param map ...
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_background_release(map)
      TYPE(qs_upf_inverse_background_type), &
         INTENT(INOUT)                                   :: map

      INTEGER                                            :: image

      IF (ASSOCIATED(map%density)) THEN
         DO image = 1, SIZE(map%density)
            CALL dbcsr_deallocate_matrix(map%density(image)%matrix)
            CALL dbcsr_deallocate_matrix(map%potential(image)%matrix)
         END DO
         DEALLOCATE (map%density, map%potential)
      END IF
      map = qs_upf_inverse_background_type()
   END SUBROUTINE qs_upf_inverse_background_release

! **************************************************************************************************
!> \brief Release owned native PW workspaces and image matrices before the environment is destroyed.
!> \param map ...
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_coulomb_release(map)
      TYPE(qs_upf_inverse_coulomb_type), INTENT(INOUT)   :: map

      IF (map%ready) THEN
         CALL qs_upf_inverse_background_release(map%background)
         CALL map%source%release()
         CALL map%potential%release()
         CALL map%source_g%release()
         CALL map%potential_g%release()
         CALL map%ionic%release()
      END IF
      map = qs_upf_inverse_coulomb_type()
   END SUBROUTINE qs_upf_inverse_coulomb_release

END MODULE qs_upf_inverse_types
