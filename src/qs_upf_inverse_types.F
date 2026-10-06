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
   USE kinds,                           ONLY: dp
   USE kpoint_transitional,             ONLY: kpoint_transitional_release,&
                                              kpoint_transitional_type
   USE kpoint_types,                    ONLY: kpoint_clear_scf_smat,&
                                              kpoint_native_density_clear,&
                                              kpoint_native_operator_clear,&
                                              kpoint_type
   USE message_passing,                 ONLY: mp_comm_type
   USE pw_grid_types,                   ONLY: PW_MODE_LOCAL,&
                                              pw_grid_type
   USE qs_upf_inverse_atom_grid,        ONLY: qs_upf_inverse_atom_grid_type
   USE upf_inverse,                     ONLY: upf_inverse_type
   USE upf_radial_fields,               ONLY: upf_radial_fields_type
   USE uspp_augmentation,               ONLY: uspp_augmentation_type

   IMPLICIT NONE
   PRIVATE

   TYPE qs_upf_inverse_type
      TYPE(upf_inverse_type) :: inverse
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
      COMPLEX(dp), ALLOCATABLE :: soft(:, :, :), gradient(:, :, :, :), operators(:, :, :, :)
      REAL(dp), ALLOCATABLE :: density(:, :), density_gradient(:, :, :), energy_density(:), &
                               vdensity(:, :), vgradient(:, :, :), background_vdensity(:, :), &
                               background_vgradient(:, :, :)
   END TYPE qs_upf_inverse_xc_type

   TYPE qs_upf_inverse_scf_type
      LOGICAL :: enabled = .FALSE.
      INTEGER :: density_evaluations = 0, operator_evaluations = 0
      TYPE(qs_upf_inverse_grid_type) :: frame
      TYPE(qs_upf_inverse_augmentation_type) :: augmentation
      TYPE(kpoint_transitional_type) :: metric
      TYPE(qs_upf_inverse_kpoints_type) :: bloch
      COMPLEX(dp), ALLOCATABLE :: metric_k(:, :, :), core_k(:, :, :), kinetic_k(:, :, :), ks_k(:, :, :, :)
      TYPE(qs_upf_inverse_atom_grid_type) :: atom_grid
      TYPE(qs_upf_inverse_xc_type) :: xc
      TYPE(pw_grid_type), POINTER :: xc_grid => NULL()
      TYPE(upf_radial_fields_type), ALLOCATABLE :: core_cache(:)
      INTEGER, ALLOCATABLE :: core_images(:, :)
      REAL(dp), ALLOCATABLE :: xc_fields(:, :, :, :), xc_kweights(:)
      COMPLEX(dp), ALLOCATABLE :: xc_matrices(:, :, :, :)
   END TYPE qs_upf_inverse_scf_type

   PUBLIC :: qs_upf_inverse_type, qs_upf_inverse_grid_type, qs_upf_inverse_augmentation_type
   PUBLIC :: qs_upf_inverse_kpoints_type, qs_upf_inverse_scf_type, qs_upf_inverse_scf_release
   PUBLIC :: qs_upf_inverse_xc_type

CONTAINS

! **************************************************************************************************
!> \brief Release the owned metric and fields before their borrowed KS basis owner is destroyed.
!> \param state Owned SCF state; pointer is nullified.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_scf_release(state)
      TYPE(qs_upf_inverse_scf_type), POINTER             :: state

      INTEGER                                            :: slot

      IF (.NOT. ASSOCIATED(state)) RETURN
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

END MODULE qs_upf_inverse_types
