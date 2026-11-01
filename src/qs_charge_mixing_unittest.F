!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM qs_charge_mixing_unittest
   USE ieee_arithmetic, ONLY: ieee_is_finite
   USE kinds, ONLY: dp
   USE message_passing, ONLY: mp_para_env_type, &
                              mp_world_finalize, &
                              mp_world_init
   USE qs_charge_mixing, ONLY: broyden_vector_mixing, charge_mixing
   USE qs_density_mixing_types, ONLY: broyden_mixing_nr, &
                                      mixing_storage_release, &
                                      mixing_storage_type
   USE qs_density_mixing_types, ONLY: modified_broyden_mixing_nr
#if defined(__TBLITE)
   USE mctc_env, ONLY: error_type
   USE tblite_scf_mixer_type, ONLY: mixer_type
   USE tblite_scc_mixer, ONLY: new_cp2k_tblite_mixer
#endif

   IMPLICIT NONE

   REAL(KIND=dp), DIMENSION(1, 1) :: charges
   TYPE(mp_para_env_type), POINTER :: para_env
   TYPE(mixing_storage_type), POINTER :: mixing_store

   NULLIFY (mixing_store, para_env)
   ALLOCATE (para_env)
   CALL mp_world_init(para_env)

   ALLOCATE (mixing_store)
   mixing_store%ncall = 0
   mixing_store%nbuffer = 4
   mixing_store%n_simple_mix = 0
   mixing_store%nskip_mixing = 0
   mixing_store%alpha = 0.1_dp
   mixing_store%broy_w0 = 0.01_dp
   mixing_store%nat_local = 1
   mixing_store%max_shell = 1
   ALLOCATE (mixing_store%acharge(1, 1, 4), mixing_store%dacharge(1, 1, 4), &
             mixing_store%dfbroy(1, 1, 4), mixing_store%ubroy(1, 1, 4), &
             mixing_store%abroy(4, 4), mixing_store%wbroy(4), mixing_store%atlist(1))
   mixing_store%acharge = 0.0_dp
   mixing_store%dacharge = 0.0_dp
   mixing_store%dfbroy = 0.0_dp
   mixing_store%ubroy = 0.0_dp
   mixing_store%abroy = 0.0_dp
   mixing_store%wbroy = 0.0_dp
   mixing_store%atlist = 1

   charges = 1.0_dp
   CALL charge_mixing(broyden_mixing_nr, mixing_store, charges, para_env, 1)
   charges = 0.5_dp
   CALL charge_mixing(broyden_mixing_nr, mixing_store, charges, para_env, 2)

   ! Reproduce a degenerate Broyden history by supplying the previous residual again.
   charges(1, 1) = mixing_store%acharge(1, 1, 2) - mixing_store%dacharge(1, 1, 1)
   CALL charge_mixing(broyden_mixing_nr, mixing_store, charges, para_env, 3)

   IF (.NOT. ALL(ieee_is_finite(charges))) THEN
      ERROR STOP "Broyden mixing produced non-finite charges for identical residuals"
   END IF
   IF (.NOT. ALL(mixing_store%dfbroy(:, :, 1) == 0.0_dp)) THEN
      ERROR STOP "Degenerate Broyden direction was not discarded"
   END IF

   ! A history can be released before all vectors exist, as before the first SCF iteration.
   ALLOCATE (mixing_store%rhoin_buffer(2, 2), mixing_store%res_buffer(2, 2), &
             mixing_store%drho_buffer(2, 2), mixing_store%last_res(2), &
             mixing_store%rhoin(2), mixing_store%rhoin_old(2), &
             mixing_store%delta_res(2, 2), mixing_store%u_vec(2, 2), mixing_store%z_vec(2, 2))
   ALLOCATE (mixing_store%rhoin_buffer(1, 1)%cc(1), mixing_store%res_buffer(1, 1)%cc(1), &
             mixing_store%drho_buffer(1, 1)%cc(1), mixing_store%last_res(1)%cc(1), &
             mixing_store%rhoin(1)%cc(1), mixing_store%rhoin_old(1)%cc(1), &
             mixing_store%delta_res(1, 1)%cc(1), mixing_store%u_vec(1, 1)%cc(1), mixing_store%z_vec(1, 1)%cc(1))
   CALL mixing_storage_release(mixing_store)
   DEALLOCATE (mixing_store)
   CALL test_vector_mixing()
#if defined(__TBLITE)
   CALL test_tblite_history()
#endif
   DEALLOCATE (para_env)
   CALL mp_world_finalize()

CONTAINS

! **************************************************************************************************
!> \brief Distributed complex fixed points, wrapped history, empty ranks and SCF resets.
! **************************************************************************************************
   SUBROUTINE test_vector_mixing()
      COMPLEX(dp), ALLOCATABLE                           :: current(:), reference(:), trial(:)
      INTEGER                                            :: i, method, nlocal, step, variable
      INTEGER, ALLOCATABLE                               :: indices(:)
      REAL(dp)                                           :: error
      REAL(dp), ALLOCATABLE                              :: damping(:)
      TYPE(mixing_storage_type), POINTER                 :: store

      nlocal = 0
      DO variable = 1, 2
         IF (MOD(variable - 1, para_env%num_pe) == para_env%mepos) nlocal = nlocal + 1
      END DO
      ALLOCATE (current(nlocal), reference(nlocal), trial(nlocal), indices(nlocal), damping(nlocal))
      i = 0
      DO variable = 1, 2
         IF (MOD(variable - 1, para_env%num_pe) /= para_env%mepos) CYCLE
         i = i + 1
         indices(i) = variable
         damping(i) = 0.1_dp + 0.05_dp*REAL(variable, dp)
         reference(i) = CMPLX(REAL(variable, dp)/3.0_dp, 0.2_dp*REAL(3 - variable, dp), dp)
      END DO
      DO method = broyden_mixing_nr, modified_broyden_mixing_nr
         ALLOCATE (store)
         store%nbuffer = 4
         store%alpha = 0.2_dp
         store%nskip_mixing = 0
         store%n_simple_mix = 0
         store%broy_w0 = 1.e-10_dp
         IF (method == modified_broyden_mixing_nr) store%broy_w0 = 1.e-5_dp
         store%wc = 0.1_dp
         store%wmax = 10.0_dp
         current(:) = CMPLX(0.0_dp, 0.0_dp, dp)
         DO step = 1, 20
            DO i = 1, nlocal
               trial(i) = reference(i) + CMPLX(0.5_dp + 0.1_dp*indices(i), 0.1_dp, dp)*(current(i) - reference(i))
            END DO
            CALL broyden_vector_mixing(method, store, current, trial, para_env, step)
            IF (.NOT. ALL(ieee_is_finite(REAL(trial, dp))) .OR. .NOT. ALL(ieee_is_finite(AIMAG(trial)))) &
               ERROR STOP 'Nonfinite distributed complex Broyden step'
            current(:) = trial
         END DO
         error = 0.0_dp
         IF (nlocal > 0) error = MAXVAL(ABS(current - reference))
         CALL para_env%max(error)
         IF (error > 1.e-8_dp) ERROR STOP 'Distributed complex Broyden fixed point'
         IF (store%vector_count /= 5) ERROR STOP 'Complex Broyden history did not wrap'
         DO step = 1, 2
            current(:) = CMPLX(0.0_dp, 0.0_dp, dp)
            trial(:) = 0.5_dp*reference
            CALL broyden_vector_mixing(method, store, current, trial, para_env, step, damping)
            error = 0.0_dp
            IF (nlocal > 0) error = MAXVAL(ABS(trial - 0.5_dp*damping*reference))
            CALL para_env%max(error)
            IF (error > 1.e-12_dp .OR. store%vector_count /= step) ERROR STOP 'Complex Broyden SCF reset or zero direction'
         END DO
         CALL mixing_storage_release(store)
         DEALLOCATE (store)
      END DO
   END SUBROUTINE test_vector_mixing

#if defined(__TBLITE)
! **************************************************************************************************
!> \brief A wrapped Broyden history must match a fresh replay of the retained input/output pairs.
! **************************************************************************************************
   SUBROUTINE test_tblite_history()
      INTEGER, PARAMETER                                 :: nstep = 16, nvar = 4

      INTEGER                                            :: config, i, j, memory, step
      REAL(KIND=dp)                                      :: damping, max_weight, min_weight, omega0, &
                                                            weight_factor
      REAL(KIND=dp), DIMENSION(nvar)                     :: actual, expected
      REAL(KIND=dp), DIMENSION(nvar, nstep)              :: qin, qout

      CLASS(mixer_type), ALLOCATABLE                    :: mixer, replay
      TYPE(error_type), ALLOCATABLE                      :: error

      DO step = 1, nstep
         DO i = 1, nvar
            qin(i, step) = SIN(0.37_dp*REAL(i*step, dp))
            qout(i, step) = qin(i, step) + 0.05_dp*COS(0.21_dp*REAL(i + 2*step, dp))
         END DO
      END DO
      DO config = 1, 2
         damping = 0.4_dp
         omega0 = 0.01_dp
         min_weight = 1.0_dp
         max_weight = 100000.0_dp
         weight_factor = 0.01_dp
         IF (config == 2) THEN
            damping = 0.25_dp
            omega0 = 0.1_dp
            min_weight = 0.5_dp
            max_weight = 3.0_dp
            weight_factor = 0.2_dp
         END IF
         DO memory = 1, 4
            CALL new_cp2k_tblite_mixer(mixer, memory, nvar, damping, omega0, min_weight, max_weight, weight_factor)
            DO step = 1, nstep
               CALL mixer%set(qin(:, step))
               CALL mixer%diff(qout(:, step))
               CALL mixer%next(error)
               IF (ALLOCATED(error)) ERROR STOP "Wrapped tblite mixer failed"
               CALL mixer%get(actual)

               ! Replay at most memory + 1 pairs: this builds exactly the retained
               ! secant directions without ever overwriting a history slot.
               CALL new_cp2k_tblite_mixer(replay, memory, nvar, damping, omega0, min_weight, max_weight, weight_factor)
               DO j = MAX(1, step - memory), step
                  CALL replay%set(qin(:, j))
                  CALL replay%diff(qout(:, j))
                  CALL replay%next(error)
                  IF (ALLOCATED(error)) ERROR STOP "Replayed tblite mixer failed"
               END DO
               CALL replay%get(expected)
               IF (.NOT. ALL(ieee_is_finite(actual))) ERROR STOP "Non-finite wrapped tblite history"
               IF (MAXVAL(ABS(actual - expected)) > 1.0E-10_dp) THEN
                  ERROR STOP "Wrapped tblite history differs from an unwrapped replay"
               END IF
            END DO
         END DO
      END DO

   END SUBROUTINE test_tblite_history
#endif

END PROGRAM qs_charge_mixing_unittest
! vim: set ts=3 sw=3 tw=132 :
