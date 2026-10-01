!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM uspp_augmentation_unittest
   USE atom_upf,                        ONLY: atom_read_upf,&
                                              atom_release_upf,&
                                              atom_upfpot_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE libcp2k,                         ONLY: cp2k_finalize,&
                                              cp2k_init

   IMPLICIT NONE
   LOGICAL :: source_root_argument

   CHARACTER(LEN=*), PARAMETER :: fixture = 'uspp-augmentation-test.upf'
   CHARACTER(LEN=default_path_length) :: filename
   INTEGER :: unit, i, j, l, n, ns, pair, channels
   REAL(KIND=dp) :: error
   REAL(KIND=dp), ALLOCATABLE :: weights(:)
   TYPE(atom_upfpot_type) :: pot
   TYPE(cp_logger_type), POINTER :: logger

   CALL cp2k_init()
   logger => cp_get_default_logger()
   IF (logger%para_env%is_source()) THEN
      OPEN (NEWUNIT=unit, FILE=fixture, STATUS='REPLACE', ACTION='WRITE')
      WRITE (unit, '(A)') &
         '<UPF version="2.0.1">', &
         '<PP_HEADER element="Pt" pseudo_type="USPP" is_ultrasoft="true"', &
         ' z_valence="16" mesh_size="3" number_of_proj="2" l_max="1" />', &
         '<PP_MESH mesh="3">', &
         '<PP_R>', &
         '0.5 1.0 1.5', &
         '</PP_R>', &
         '<PP_RAB>', &
         '0.5 0.5 0.5', &
         '</PP_RAB>', &
         '</PP_MESH>', &
         '<PP_LOCAL type="real" size="3">', &
         '6 4 2', &
         '</PP_LOCAL>', &
         '<PP_NONLOCAL>', &
         '<PP_BETA.1 type="real" size="3" index="1" angular_momentum="0">', &
         '1 4 9', &
         '</PP_BETA.1>', &
         '<PP_BETA.2 type="real" size="3" index="2" angular_momentum="1">', &
         '2 -1 3', &
         '</PP_BETA.2>', &
         '<PP_DIJ>', &
         '2 0 0 8', &
         '</PP_DIJ>', &
         '<PP_AUGMENTATION q_with_l="T" nqf="0" nqlc="3">', &
         '<PP_Q>', &
         '2 0 0 3', &
         '</PP_Q>', &
         '<PP_QIJL.2.2.0 first_index="2" second_index="2" composite_index="3" angular_momentum="0">', &
         '-3 4 5', &
         '</PP_QIJL.2.2.0>', &
         '<PP_QIJL.1.1.0 type="real" size="3" first_index="1" second_index="1"', &
         ' composite_index="1" angular_momentum="0" is_null="false">', &
         '1 2 3', &
         '</PP_QIJL.1.1.0>', &
         '<PP_QIJL.1.2.1 first_index="1" second_index="2" composite_index="2" angular_momentum="1">', &
         '4 -5 6', &
         '</PP_QIJL.1.2.1>', &
         '<PP_QIJL.2.2.2 first_index="2" second_index="2" composite_index="3"', &
         ' angular_momentum="2" is_null="T" size="0"/>', &
         '</PP_AUGMENTATION>', &
         '</PP_NONLOCAL>', &
         '</UPF>'
      CLOSE (unit)
   END IF
   CALL logger%para_env%sync()
   CALL atom_read_upf(pot, fixture, read_augmentation=.TRUE.)
   IF (ANY(pot%econf(0:3) /= [1, 6, 9, 0])) ERROR STOP 'Wrong 62-electron core subtraction'
   IF (ANY(SHAPE(pot%qijl) /= [3, 3, 3]) .OR. LBOUND(pot%qijl, 3) /= 0) ERROR STOP 'Wrong augmentation shape'
   IF (MAXVAL(ABS(pot%qij - RESHAPE([2._dp, 0._dp, 0._dp, 3._dp], [2, 2]))) > 1.e-14_dp) THEN
      ERROR STOP 'PP_Q normalization changed'
   END IF
   IF (MAXVAL(ABS(pot%qijl(:, 1, 0) - [1._dp, 2._dp, 3._dp])) > 1.e-14_dp) ERROR STOP 'Wrong s-s augmentation'
   IF (MAXVAL(ABS(pot%qijl(:, 2, 1) - [4._dp, -5._dp, 6._dp])) > 1.e-14_dp) ERROR STOP 'Wrong s-p augmentation'
   IF (MAXVAL(ABS(pot%qijl(:, 3, 0) - [-3._dp, 4._dp, 5._dp])) > 1.e-14_dp) ERROR STOP 'Wrong p-p augmentation'
   IF (ANY(pot%qijl(:, :, 2) /= 0._dp)) ERROR STOP 'Null augmentation was not zero initialized'
   IF (MAXVAL(ABS(pot%beta(:, 2) - [4._dp, -1._dp, 2._dp])) > 1.e-14_dp) ERROR STOP 'USPP beta normalization changed'
   IF (ABS(pot%dion(2, 2) - 4._dp) > 1.e-14_dp) ERROR STOP 'USPP D is not in Hartree'
   CALL atom_release_upf(pot)
   IF (ALLOCATED(pot%qij) .OR. ALLOCATED(pot%qijl)) ERROR STOP 'Augmentation arrays were not released'
   CALL logger%para_env%sync()
   IF (logger%para_env%is_source()) THEN
      OPEN (NEWUNIT=unit, FILE=fixture, STATUS='OLD')
      CLOSE (unit, STATUS='DELETE')
   END IF

   ! Optional SSSP files: export invariants for an independent XML-reader comparison.
   DO i = 1, COMMAND_ARGUMENT_COUNT()
      CALL GET_COMMAND_ARGUMENT(i, filename)
      ! do_regtest.py passes the CP2K source root before any optional UPF audit files.
      IF (i == 1) THEN
         INQUIRE (FILE=TRIM(filename)//"/tests/UNIT_TESTS", EXIST=source_root_argument)
         IF (source_root_argument) CYCLE
      END IF
      CALL atom_read_upf(pot, TRIM(filename), read_augmentation=.TRUE.)
      IF (.NOT. ALL(ieee_is_finite(pot%qij)) .OR. .NOT. ALL(ieee_is_finite(pot%qijl))) THEN
         ERROR STOP 'Nonfinite USPP augmentation'
      END IF
      n = SIZE(pot%r)
      ALLOCATE (weights(n))
      weights = 0._dp
      ns = n - 1 + MODULO(n, 2)
      weights(1:ns:2) = 2._dp/3._dp
      weights(2:ns - 1:2) = 4._dp/3._dp
      weights(1) = 1._dp/3._dp
      weights(ns) = 1._dp/3._dp
      IF (ns < n) weights(n - 1:n) = weights(n - 1:n) + .5_dp
      weights(:) = weights*pot%rab
      error = 0._dp
      channels = 0
      DO j = 1, pot%number_of_proj
         DO l = 1, j
            pair = j*(j - 1)/2 + l
            channels = channels + MIN(pot%lbeta(j), pot%lbeta(l)) + 1
            IF (pot%lbeta(l) == pot%lbeta(j)) THEN
               error = MAX(error, ABS(SUM(weights*pot%qijl(:, pair, 0)) - pot%qij(l, j)))
            END IF
         END DO
      END DO
      IF (logger%para_env%is_source()) THEN
         WRITE (*, '(A,1X,A,1X,I0,4(1X,ES25.16))') 'USPP augmentation:', TRIM(filename), channels, &
            SUM(pot%qij), SUM(pot%qijl), SUM(ABS(pot%qijl)), error
      END IF
      DEALLOCATE (weights)
      CALL atom_release_upf(pot)
   END DO
   CALL cp2k_finalize()
END PROGRAM uspp_augmentation_unittest
