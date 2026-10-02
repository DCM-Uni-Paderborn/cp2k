!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief PAW one-center cache types without dependencies on Kohn-Sham or XC evaluators.
! **************************************************************************************************
MODULE paw_types
   USE kinds,                           ONLY: dp
   USE paw_core_kinetic,                ONLY: paw_core_kinetic_type

   IMPLICIT NONE
   PRIVATE

   TYPE paw_onecenter_type
      INTEGER :: nproj = 0, npair = 0, lmax = -1
      INTEGER, ALLOCATABLE :: pair_index(:, :)
      REAL(dp), ALLOCATABLE :: gaunt(:, :, :), hartree(:, :, :)
      ! Lazy original-data frozen-core exchange; matrix is the unit-fraction spin operator.
      LOGICAL :: core_exchange_initialized = .FALSE., core_exchange_available = .FALSE.
      CHARACTER(LEN=256) :: core_exchange_reason = 'PAW core exchange has not been initialized'
      REAL(dp) :: core_exchange_energy = 0.0_dp
      REAL(dp), ALLOCATABLE :: core_exchange_matrix(:, :)
   END TYPE paw_onecenter_type

   TYPE paw_xc_type
      TYPE(paw_core_kinetic_type) :: core_kinetic
      INTEGER :: na = 0, nr = 0, nval = 0, npair = 0, lmax = -1
      REAL(dp), ALLOCATABLE :: r(:), wr(:), wa(:), directions(:, :), harmonics(:, :), angular_gradient(:, :, :)
      ! Final dimension: 1 is AE, 2 is compensated PS. Radial arrays contain true densities.
      REAL(dp), ALLOCATABLE :: radial(:, :, :, :), radial_gradient(:, :, :, :), core(:, :), core_gradient(:, :)
      ! Valence partial-wave kinetic density, 1/2 grad(phi_i).grad(phi_j), without charge compensation.
      REAL(dp), ALLOCATABLE :: kinetic_radial(:, :, :, :)
   END TYPE paw_xc_type

   TYPE paw_kind_type
      LOGICAL :: active = .FALSE.
      REAL(dp) :: core_energy = 0.0_dp
      TYPE(paw_onecenter_type) :: onecenter
      TYPE(paw_xc_type) :: xcgrid
   END TYPE paw_kind_type

   PUBLIC :: paw_onecenter_type, paw_xc_type, paw_kind_type

END MODULE paw_types
