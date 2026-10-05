!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Accepted and trial complete complex densities in the native Bloch variational space.
!>        No inverse transformation to a finite real-space image set is used for mixing.
! **************************************************************************************************
MODULE qs_kpoint_density_state
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_get_submatrix
   USE kinds,                           ONLY: dp
   USE kpoint_types,                    ONLY: kpoint_env_type,&
                                              kpoint_type
   USE qs_kpoint_operators,             ONLY: kpoint_slot_s
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: kpoint_density_state_active
   PUBLIC :: kpoint_density_state_initialize, kpoint_density_state_capture
   PUBLIC :: kpoint_density_state_accept, kpoint_density_state_delta, kpoint_density_state_restore
   PUBLIC :: kpoint_density_state_mulliken
CONTAINS

! **************************************************************************************************
!> \brief Test explicit complete Bloch density storage; null Gamma pointers are safe.
!> \param kpoints ...
!> \return ...
! **************************************************************************************************
   LOGICAL FUNCTION kpoint_density_state_active(kpoints) RESULT(active)
      TYPE(kpoint_type), POINTER                         :: kpoints

      active = .FALSE.
      IF (ASSOCIATED(kpoints)) active = ALLOCATED(kpoints%native_density)
   END FUNCTION kpoint_density_state_active

! **************************************************************************************************
!> \brief Explicitly attach group-local accepted/trial density storage to a native Bloch state.
!> \param kpoints ...
! **************************************************************************************************
   SUBROUTINE kpoint_density_state_initialize(kpoints)
      TYPE(kpoint_type), INTENT(INOUT)                   :: kpoints

      INTEGER                                            :: ik, nspin, spin

      CPASSERT(.NOT. ALLOCATED(kpoints%native_density))
      nspin = SIZE(kpoints%kp_env(1)%kpoint_env%pmat, 2)
      ALLOCATE (kpoints%native_density(SIZE(kpoints%kp_env), nspin), &
                kpoints%native_trial(SIZE(kpoints%kp_env), nspin))
      DO spin = 1, nspin
         DO ik = 1, SIZE(kpoints%kp_env)
            CALL cp_cfm_create(kpoints%native_density(ik, spin), &
                               kpoints%kp_env(ik)%kpoint_env%pmat(1, spin)%matrix_struct)
            CALL cp_cfm_create(kpoints%native_trial(ik, spin), &
                               kpoints%kp_env(ik)%kpoint_env%pmat(1, spin)%matrix_struct)
            kpoints%native_density(ik, spin)%local_data = CMPLX(0.0_dp, 0.0_dp, dp)
            kpoints%native_trial(ik, spin)%local_data = CMPLX(0.0_dp, 0.0_dp, dp)
         END DO
      END DO
      kpoints%native_trial_valid = .FALSE.
      ! A representation-specific generalized solve seeds P(k) in the physical metric.
      kpoints%native_density_seeded = .FALSE.
   END SUBROUTINE kpoint_density_state_initialize

! **************************************************************************************************
!> \brief Capture the newly diagonalized density while retaining the density that built H.
!> \param kpoints ...
! **************************************************************************************************
   SUBROUTINE kpoint_density_state_capture(kpoints)
      TYPE(kpoint_type), INTENT(INOUT)                   :: kpoints

      INTEGER                                            :: ik, spin
      TYPE(kpoint_env_type), POINTER                     :: kp

      CPASSERT(ALLOCATED(kpoints%native_density))
      DO spin = 1, SIZE(kpoints%native_trial, 2)
         DO ik = 1, SIZE(kpoints%native_trial, 1)
            kp => kpoints%kp_env(ik)%kpoint_env
            kpoints%native_trial(ik, spin)%local_data(:, :) = CMPLX(kp%pmat(1, spin)%local_data, 0.0_dp, dp)
            IF (.NOT. kpoints%use_real_wfn) THEN
               kpoints%native_trial(ik, spin)%local_data(:, :) = &
                  kpoints%native_trial(ik, spin)%local_data + CMPLX(0.0_dp, kp%pmat(2, spin)%local_data, dp)
            END IF
         END DO
      END DO
      kpoints%native_trial_valid = .TRUE.
   END SUBROUTINE kpoint_density_state_capture

! **************************************************************************************************
!> \brief Maximum complex AO density change, reduced over the world once without k weights.
!> \param kpoints ...
!> \param delta ...
! **************************************************************************************************
   SUBROUTINE kpoint_density_state_delta(kpoints, delta)
      TYPE(kpoint_type), INTENT(IN)                      :: kpoints
      REAL(dp), INTENT(OUT)                              :: delta

      INTEGER                                            :: ik, spin

      CPASSERT(kpoints%native_trial_valid)
      delta = 0.0_dp
      DO spin = 1, SIZE(kpoints%native_trial, 2)
         DO ik = 1, SIZE(kpoints%native_trial, 1)
            delta = MAX(delta, MAXVAL(ABS(kpoints%native_trial(ik, spin)%local_data - &
                                          kpoints%native_density(ik, spin)%local_data)))
         END DO
      END DO
      CALL kpoints%para_env%max(delta)
   END SUBROUTINE kpoint_density_state_delta

! **************************************************************************************************
!> \brief Accept a convex direct density step. The raw trial remains available for finalization.
!> \param kpoints ...
!> \param alpha ...
! **************************************************************************************************
   SUBROUTINE kpoint_density_state_accept(kpoints, alpha)
      TYPE(kpoint_type), INTENT(INOUT)                   :: kpoints
      REAL(dp), INTENT(IN)                               :: alpha

      INTEGER                                            :: ik, spin

      CPASSERT(kpoints%native_trial_valid)
      CPASSERT(alpha >= 0.0_dp .AND. alpha <= 1.0_dp)
      DO spin = 1, SIZE(kpoints%native_trial, 2)
         DO ik = 1, SIZE(kpoints%native_trial, 1)
            IF (alpha == 1.0_dp) THEN
               kpoints%native_density(ik, spin)%local_data(:, :) = kpoints%native_trial(ik, spin)%local_data
            ELSE
               kpoints%native_density(ik, spin)%local_data(:, :) = &
                  (1.0_dp - alpha)*kpoints%native_density(ik, spin)%local_data + &
                  alpha*kpoints%native_trial(ik, spin)%local_data
            END IF
         END DO
      END DO
      CALL kpoint_density_state_restore(kpoints)
   END SUBROUTINE kpoint_density_state_accept

! **************************************************************************************************
!> \brief Restore the accepted P(k) used by the next native commutator/DIIS evaluation.
!> \param kpoints ...
! **************************************************************************************************
   SUBROUTINE kpoint_density_state_restore(kpoints)
      TYPE(kpoint_type), INTENT(INOUT)                   :: kpoints

      INTEGER                                            :: ik, spin
      TYPE(kpoint_env_type), POINTER                     :: kp

      CPASSERT(ALLOCATED(kpoints%native_density))
      DO spin = 1, SIZE(kpoints%native_density, 2)
         DO ik = 1, SIZE(kpoints%native_density, 1)
            kp => kpoints%kp_env(ik)%kpoint_env
            kp%pmat(1, spin)%local_data(:, :) = REAL(kpoints%native_density(ik, spin)%local_data, dp)
            IF (.NOT. kpoints%use_real_wfn) THEN
               kp%pmat(2, spin)%local_data(:, :) = AIMAG(kpoints%native_density(ik, spin)%local_data)
            END IF
         END DO
      END DO
   END SUBROUTINE kpoint_density_state_restore

! **************************************************************************************************
!> \brief Mulliken populations from the accepted complex density and its registered charge metric.
!>        Each group source contributes its actual k-point weights once. The caller sums world ranks.
!> \param kpoints ...
!> \param populations Partial orbital populations, overwritten.
! **************************************************************************************************
   SUBROUTINE kpoint_density_state_mulliken(kpoints, populations)
      TYPE(kpoint_type), INTENT(IN)                      :: kpoints
      REAL(dp), INTENT(OUT)                              :: populations(:, :)

      COMPLEX(dp), ALLOCATABLE                           :: density(:, :), metric(:, :)
      INTEGER                                            :: ik, local, n, spin

      CPASSERT(ALLOCATED(kpoints%native_density))
      CPASSERT(ALLOCATED(kpoints%native_operators(kpoint_slot_s)%matrices))
      n = SIZE(populations, 1)
      CPASSERT(SIZE(populations, 2) == SIZE(kpoints%native_density, 2))
      ALLOCATE (density(n, n), metric(n, n))
      populations = 0.0_dp
      DO local = 1, SIZE(kpoints%native_density, 1)
         ik = kpoints%kp_range(1) + local - 1
         CALL cp_cfm_get_submatrix(kpoints%native_operators(kpoint_slot_s)%matrices(local, 1), metric)
         DO spin = 1, SIZE(populations, 2)
            CALL cp_cfm_get_submatrix(kpoints%native_density(local, spin), density)
            IF (kpoints%para_env_kp%is_source()) THEN
               populations(:, spin) = populations(:, spin) + &
                                      kpoints%wkp(ik)*SUM(REAL(density*CONJG(metric), dp), DIM=2)
            END IF
         END DO
      END DO
   END SUBROUTINE kpoint_density_state_mulliken
END MODULE qs_kpoint_density_state
