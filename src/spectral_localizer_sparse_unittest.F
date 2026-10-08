!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM spectral_localizer_sparse_unittest
   USE cp_dbcsr_api,                    ONLY: dbcsr_finalize_lib,&
                                              dbcsr_init_lib
   USE cp_log_handling,                 ONLY: cp_add_default_logger,&
                                              cp_logger_create,&
                                              cp_logger_release,&
                                              cp_logger_type,&
                                              cp_rm_default_logger
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE localizer_dbcsr_factor,          ONLY: localizer_factor_release,&
                                              localizer_factor_solve,&
                                              localizer_factor_type,&
                                              localizer_factorize
   USE localizer_sparse_inertia,        ONLY: sparse_inertia_available,&
                                              sparse_inertia_count,&
                                              sparse_inertia_finish,&
                                              sparse_inertia_start,&
                                              sparse_inertia_type,&
                                              sparse_metric_solve,&
                                              sparse_metric_start
   USE localizer_sparse_pfaffian,       ONLY: sparse_pfaffian_available,&
                                              sparse_pfaffian_sign,&
                                              sparse_z2_evaluate
   USE machine,                         ONLY: default_output_unit
   USE message_passing,                 ONLY: mp_para_env_type,&
                                              mp_world_finalize,&
                                              mp_world_init
   USE spectral_localizer,              ONLY: localizer_matrix_evaluate,&
                                              localizer_pfaffian_sign
   USE spectral_localizer_sparse,       ONLY: localizer_sparse_evaluate

   IMPLICIT NONE

   TYPE(mp_para_env_type), POINTER :: para_env
   TYPE(cp_logger_type), POINTER :: logger
   COMPLEX(KIND=dp), ALLOCATABLE :: matrix(:, :), metric(:, :), h(:, :)
   COMPLEX(KIND=dp) :: tx(2, 2), ty(2, 2), coordinate
   INTEGER :: length, n, i, j, x, y, io
   REAL(KIND=dp) :: factor
   LOGICAL :: pf_native = .FALSE.

   ALLOCATE (para_env)
   CALL mp_world_init(para_env)
   io = -1
   IF (para_env%is_source()) io = default_output_unit
   CALL cp_logger_create(logger, para_env=para_env, default_global_unit_nr=io, close_global_unit_on_dealloc=.FALSE.)
   CALL cp_add_default_logger(logger)
   CALL dbcsr_init_lib(para_env%get_handle(), io)
   CALL native_skew_tests()
   CALL pfaffian_tests(.TRUE.)
   IF (sparse_pfaffian_available()) CALL pfaffian_tests(.FALSE.)
   IF (sparse_inertia_available()) THEN
      CALL factor_tests()
      CALL metric_tests()
      CALL shifted_gap_tests()
      length = 5
      n = 2*length**2
      ALLOCATE (h(n, n), metric(n, n), matrix(2*n, 2*n))
      h = 0.0_dp
      metric = 0.0_dp
      matrix = 0.0_dp
      tx = RESHAPE([CMPLX(0.5_dp, 0.0_dp, dp), CMPLX(0.0_dp, -0.5_dp, dp), &
                    CMPLX(0.0_dp, -0.5_dp, dp), CMPLX(-0.5_dp, 0.0_dp, dp)], [2, 2])
      ty = RESHAPE([0.5_dp, 0.5_dp, -0.5_dp, -0.5_dp], [2, 2])
      DO y = 1, length
         DO x = 1, length
            i = 2*((y - 1)*length + x) - 1
            h(i, i) = -1.0_dp
            h(i + 1, i + 1) = 1.0_dp
            metric(i, i) = 1.0_dp
            metric(i + 1, i + 1) = 1.0_dp
            coordinate = 0.25_dp*CMPLX(x - 0.5_dp*(length + 1), -(y - 0.5_dp*(length + 1)), dp)
            matrix(i, n + i) = coordinate
            matrix(i + 1, n + i + 1) = coordinate
            IF (x < length) THEN
               j = i + 2
               h(i:i + 1, j:j + 1) = tx
               h(j:j + 1, i:i + 1) = CONJG(TRANSPOSE(tx))
            END IF
            IF (y < length) THEN
               j = i + 2*length
               h(i:i + 1, j:j + 1) = ty
               h(j:j + 1, i:i + 1) = CONJG(TRANSPOSE(ty))
            END IF
         END DO
      END DO
      matrix(:n, :n) = h
      matrix(n + 1:, n + 1:) = -h
      matrix(n + 1:, :n) = CONJG(TRANSPOSE(matrix(:n, n + 1:)))
      CALL compare(matrix, metric, -1, 'QWZ')
      CALL aii_tests(h, matrix, metric)
      DO i = 1, n
         factor = 0.7_dp + 0.1_dp*MOD(i, 5)
         matrix(i, :) = factor*matrix(i, :)
         matrix(:, i) = factor*matrix(:, i)
         matrix(n + i, :) = factor*matrix(n + i, :)
         matrix(:, n + i) = factor*matrix(:, n + i)
         metric(i, i) = factor**2
      END DO
      CALL compare(matrix, metric, -1, 'QWZ nonorthogonal')
      DEALLOCATE (h, metric, matrix)
      ALLOCATE (matrix(4, 4), metric(2, 2))
      metric = 0.0_dp
      metric(1, 1) = 1.0_dp
      metric(2, 2) = 1.0_dp
      matrix = 0.0_dp
      matrix(1, 1) = 1.e-6_dp
      matrix(2, 2) = 2.0_dp
      matrix(3, 3) = -matrix(1, 1)
      matrix(4, 4) = -2.0_dp
      CALL compare(matrix, metric, 0, 'near closing')
      matrix(1, 1) = 1.e-10_dp
      matrix(3, 3) = -matrix(1, 1)
      CALL compare(matrix, metric, 0, 'unresolved')
      matrix(1, 1) = 0.0_dp
      matrix(3, 3) = 0.0_dp
      CALL compare(matrix, metric, 0, 'singular')
      matrix = 0.0_dp
      matrix(1, 2) = CMPLX(0.0_dp, 1.0_dp, dp)
      matrix(2, 1) = CONJG(matrix(1, 2))
      matrix(3, 4) = matrix(1, 2)
      matrix(4, 3) = matrix(2, 1)
      CALL compare(matrix, metric, 0, 'complex zero-diagonal pivots')
      tx = RESHAPE([CMPLX(1.0_dp, 0.0_dp, dp), CMPLX(0.0_dp, 0.0_dp, dp), &
                    CMPLX(0.0_dp, 0.2_dp, dp), CMPLX(0.8_dp, 0.0_dp, dp)], [2, 2])
      matrix(:2, :2) = MATMUL(CONJG(TRANSPOSE(tx)), MATMUL(matrix(:2, :2), tx))
      matrix(3:, 3:) = matrix(:2, :2)
      metric(:, :) = MATMUL(CONJG(TRANSPOSE(tx)), tx)
      CALL compare(matrix, metric, 0, 'complex metric')
      metric = 0.0_dp
      metric(1, 1) = 1.0_dp
      metric(2, 2) = -1.0_dp
      CALL compare(matrix, metric, 0, 'indefinite metric')
      metric(2, 2) = 0.0_dp
      CALL compare(matrix, metric, 0, 'singular metric')
      metric(2, 2) = 1.e-12_dp
      CALL compare(matrix, metric, 0, 'ill-conditioned metric')
      DEALLOCATE (matrix, metric)
   END IF
   CALL dbcsr_finalize_lib()
   CALL cp_rm_default_logger()
   CALL cp_logger_release(logger)
   DEALLOCATE (para_env)
   CALL mp_world_finalize()

CONTAINS

! **************************************************************************************************
!> \brief Distributed skew pivots preserve permutation parity, scale and repeated solves.
! **************************************************************************************************
   SUBROUTINE native_skew_tests()
      INTEGER                                            :: case_number, count, dense_sign, &
                                                            dense_status, entry_index, i, j, k, n, &
                                                            negative, pass, sign, status, swap
      INTEGER, ALLOCATABLE                               :: cols(:), rows(:)
      REAL(dp)                                           :: scale
      REAL(dp), ALLOCATABLE                              :: a(:, :), exact(:, :), rhs(:, :), &
                                                            values(:), work(:)
      TYPE(localizer_factor_type)                        :: factor

      DO n = 2, 34, 8
         ALLOCATE (a(n, n), exact(n, 2), rhs(n, 2), work(n))
         DO i = 1, n
            exact(i, :) = [SIN(REAL(i, dp)), COS(REAL(i, dp))]
         END DO
         DO case_number = 1, 4
            a = 0.0_dp
            DO j = 2, n
               DO i = 1, j - 1
                  IF (MOD(7*i + 11*j + case_number, 3) /= 0) CYCLE
                  a(i, j) = SIN(REAL(17*i + 13*j + case_number, dp))
               END DO
            END DO
            DO i = 1, n - 1, 2
               a(i, i + 1) = a(i, i + 1) + 2.0_dp
            END DO
            a = a - TRANSPOSE(a)
            CALL localizer_pfaffian_sign(a, dense_sign, dense_status)
            IF (dense_status /= 0) ERROR STOP 'Invalid dense skew fixture'
            DO swap = 1, 2
               DO k = -1, 1
                  scale = 10.0_dp**(120*k)
                  DO pass = 1, 2
                     count = 0
                     entry_index = 0
                     DO j = 1, n
                        DO i = j + 1, n
                           IF (a(i, j) == 0.0_dp) CYCLE
                           entry_index = entry_index + 1
                           IF (MOD(entry_index - 1, para_env%num_pe) /= para_env%mepos) CYCLE
                           count = count + 1
                           IF (pass == 1) CYCLE
                           rows(count) = i
                           cols(count) = j
                           values(count) = scale*a(i, j)
                        END DO
                     END DO
                     IF (pass == 1) ALLOCATE (rows(count), cols(count), values(count))
                  END DO
                  CALL localizer_factorize(n, rows, cols, values, para_env, factor, negative, status, .TRUE., &
                                           skew=.TRUE., tolerance=1.e-10_dp, pf_sign=sign)
                  IF (status /= 0 .OR. sign /= dense_sign) ERROR STOP 'Native skew sign mismatch'
                  rhs = scale*MATMUL(a, exact)
                  CALL localizer_factor_solve(factor, rhs, para_env, status)
                  IF (status /= 0) ERROR STOP 'Native skew solve failed'
                  IF (MAXVAL(ABS(rhs - exact)) > 1.e-10_dp) ERROR STOP 'Native skew solution mismatch'
                  IF (MAXVAL(ABS(MATMUL(a, rhs - exact))) > 1.e-10_dp) ERROR STOP 'Native skew residual mismatch'
                  CALL localizer_factor_release(factor)
                  DEALLOCATE (rows, cols, values)
               END DO
               work = a(1, :)
               a(1, :) = a(2, :)
               a(2, :) = work
               work = a(:, 1)
               a(:, 1) = a(:, 2)
               a(:, 2) = work
               dense_sign = -dense_sign
            END DO
         END DO
         a(1, :) = 0.0_dp
         a(:, 1) = 0.0_dp
         count = n*(n - 1)/2
         ALLOCATE (rows(count), cols(count), values(count))
         count = 0
         DO j = 1, n
            DO i = j + 1, n
               count = count + 1
               rows(count) = i
               cols(count) = j
               values(count) = a(i, j)/REAL(para_env%num_pe, dp)
            END DO
         END DO
         CALL localizer_factorize(n, rows, cols, values, para_env, factor, negative, status, .TRUE., &
                                  skew=.TRUE., tolerance=1.e-10_dp, pf_sign=sign)
         IF (status /= 1 .OR. sign /= 0 .OR. factor%valid) ERROR STOP 'Singular native skew factor accepted'
         CALL localizer_factor_release(factor)
         DEALLOCATE (a, exact, rhs, work, rows, cols, values)
      END DO
      IF (para_env%is_source()) WRITE (*, '(A)') 'Native DBCSR skew signs, scales, permutations and solves passed'
   END SUBROUTINE native_skew_tests

! **************************************************************************************************
!> \brief Distributed duplicates and indefinite 2x2 factors must preserve inertia and solves.
! **************************************************************************************************
   SUBROUTINE factor_tests()
      INTEGER, PARAMETER                                 :: n = 34

      INTEGER                                            :: count, entry_index, i, j, k, negative, &
                                                            pass, status
      INTEGER, ALLOCATABLE                               :: cols(:), rows(:)
      REAL(dp)                                           :: a(n, n), exact(n, 2), rhs(n, 2)
      REAL(dp), ALLOCATABLE                              :: values(:)
      TYPE(localizer_factor_type)                        :: factor
      TYPE(sparse_inertia_type)                          :: state

      a = 0.0_dp
      DO i = 1, n - 1
         a(i + 1, i) = 1.0_dp + 0.01_dp*i
         a(i, i + 1) = a(i + 1, i)
      END DO
      DO j = 1, 2
         DO i = 1, n
            exact(i, j) = SIN(REAL(i + 7*j, dp))
         END DO
      END DO
      DO pass = 1, 2
         count = 0
         entry_index = 0
         DO i = 1, n - 1
            DO k = 1, 2
               entry_index = entry_index + 1
               IF (MOD(entry_index - 1, para_env%num_pe) /= para_env%mepos) CYCLE
               count = count + 1
               IF (pass == 1) CYCLE
               rows(count) = i + 1
               cols(count) = i
               values(count) = 0.5_dp*a(i + 1, i)
            END DO
         END DO
         IF (pass == 1) ALLOCATE (rows(count), cols(count), values(count))
      END DO
      CALL localizer_factorize(n, rows, cols, values, para_env, factor, negative, status, .TRUE.)
      IF (status /= 0 .OR. negative /= n/2) ERROR STOP 'Native indefinite factorization failed'
      rhs = MATMUL(a, exact)
      CALL localizer_factor_solve(factor, rhs, para_env, status)
      IF (status /= 0) ERROR STOP 'Native indefinite solve failed'
      IF (MAXVAL(ABS(rhs - exact)) > 1.e-12_dp) ERROR STOP 'Native indefinite solution mismatch'
      IF (MAXVAL(ABS(MATMUL(a, rhs - exact))) > 1.e-12_dp) ERROR STOP 'Native indefinite residual mismatch'
      CALL localizer_factor_release(factor)

      IF (para_env%is_source()) rows(1) = 0
      CALL sparse_inertia_start(state, n, rows, cols, para_env, status)
      IF (status /= -1) ERROR STOP 'Invalid rank-local coordinates accepted'
      IF (para_env%is_source()) rows(1) = 2
      CALL sparse_inertia_start(state, n, rows, cols, para_env, status)
      IF (status /= 0) ERROR STOP 'Distributed inertia initialization failed'
      IF (para_env%is_source()) state%initialized = .FALSE.
      CALL sparse_inertia_count(state, values, para_env, negative, status)
      IF (status /= -4) ERROR STOP 'Invalid rank-local inertia context accepted'
      IF (para_env%is_source()) state%initialized = .TRUE.
      CALL sparse_inertia_count(state, values, para_env, negative, status)
      IF (status /= 0 .OR. negative /= n/2) ERROR STOP 'Duplicate inertia mismatch'
      CALL sparse_inertia_finish(state)
      CALL sparse_inertia_start(state, n, rows, cols, para_env, status, use_native=.FALSE.)
      IF (sparse_inertia_available(use_native=.FALSE.)) THEN
         IF (status /= 0) ERROR STOP 'Optional MUMPS initialization failed'
         CALL sparse_inertia_count(state, values, para_env, negative, status)
         IF (status /= 0 .OR. negative /= n/2) ERROR STOP 'Native/MUMPS inertia mismatch'
         rhs = MATMUL(a, exact)
         CALL sparse_metric_solve(state, rhs, para_env, status)
         IF (status /= 0) ERROR STOP 'Optional MUMPS solve failed'
         IF (MAXVAL(ABS(rhs - exact)) > 1.e-12_dp) ERROR STOP 'Native/MUMPS solution mismatch'
         CALL sparse_inertia_finish(state)
      ELSE
         IF (status /= -5) ERROR STOP 'Unavailable MUMPS accepted'
      END IF
   END SUBROUTINE factor_tests

! **************************************************************************************************
!> \brief Shifted LDL brackets must track generalized eigenvalues, not pivot sizes.
! **************************************************************************************************
   SUBROUTINE shifted_gap_tests()
      COMPLEX(dp)                                        :: a(4, 4), change(2, 2), h(2, 2), s(2, 2)
      INTEGER                                            :: i
      REAL(dp)                                           :: crossing(5), exact_gap

      h = RESHAPE([1.0_dp, 20.0_dp, 20.0_dp, 399.0_dp], [2, 2])
      exact_gap = 1.0_dp/(SQRT(160004.0_dp) + 400.0_dp)
      s = 0.0_dp
      s(1, 1) = 2.0_dp
      s(2, 2) = 2.0_dp
      a = 0.0_dp
      a(:2, :2) = h
      a(3:, 3:) = -h
      CALL compare(a, s, 0, 'gap not pivot size', expected_gap=exact_gap)
      change = RESHAPE([0.5_dp, 0.0_dp, 0.2_dp, 3.0_dp], [2, 2])
      h = MATMUL(TRANSPOSE(change), MATMUL(h, change))
      s = MATMUL(TRANSPOSE(change), MATMUL(s, change))
      a(:2, :2) = h
      a(3:, 3:) = -h
      CALL compare(a, s, 0, 'gap not pivot size nonorthogonal', expected_gap=exact_gap)
      crossing = [-1.e-6_dp, -1.e-12_dp, 0.0_dp, 1.e-12_dp, 1.e-6_dp]
      s = 0.0_dp
      s(1, 1) = 2.0_dp
      s(2, 2) = 2.0_dp
      a = 0.0_dp
      a(2, 2) = 5.0_dp
      a(4, 4) = -5.0_dp
      DO i = 1, SIZE(crossing)
         a(1, 1) = crossing(i)
         a(3, 3) = -crossing(i)
         CALL compare(a, s, 0, 'two-sided gap crossing')
      END DO
   END SUBROUTINE shifted_gap_tests

! **************************************************************************************************
!> \brief Multiple right hand sides reuse validated metric factors, independent of rank count.
! **************************************************************************************************
   SUBROUTINE metric_tests()
      INTEGER                                            :: count, i, info, j, k, nfactor, pass
      INTEGER, ALLOCATABLE                               :: cols(:), rows(:)
      REAL(dp)                                           :: exact(6, 3), ratio, rhs(6, 3), s(6, 6)
      REAL(dp), ALLOCATABLE                              :: identity(:), values(:)
      TYPE(sparse_inertia_type)                          :: state

      s = 0.0_dp
      DO i = 1, 6
         s(i, i) = REAL(i + 1, dp)
         IF (i < 6) THEN
            s(i, i + 1) = 0.2_dp
            s(i + 1, i) = 0.2_dp
         END IF
         DO j = 1, 3
            exact(i, j) = SIN(REAL(i + 7*j, dp))
         END DO
      END DO
      DO pass = 1, 2
         k = 0
         count = 0
         DO j = 1, 6
            DO i = j, 6
               IF (s(i, j) == 0.0_dp) CYCLE
               k = k + 1
               IF (MOD(k - 1, para_env%num_pe) /= para_env%mepos) CYCLE
               count = count + 1
               IF (pass == 2) THEN
                  rows(count) = i
                  cols(count) = j
                  values(count) = s(i, j)
                  identity(count) = 0.0_dp
                  IF (i == j) identity(count) = 1.0_dp
               END IF
            END DO
         END DO
         IF (pass == 1) ALLOCATE (rows(count), cols(count), values(count), identity(count))
      END DO
      CALL sparse_metric_start(state, 6, rows, cols, values, identity, 1.e-10_dp, para_env, ratio, info)
      IF (info /= 0 .OR. ratio <= 1.e-10_dp) ERROR STOP 'Metric factorization failed'
      nfactor = state%factorizations
      DO i = 1, 2
         rhs = MATMUL(s, REAL(i, dp)*exact)
         CALL sparse_metric_solve(state, rhs, para_env, info)
         IF (info /= 0) ERROR STOP 'Repeated metric solution failed'
         IF (MAXVAL(ABS(rhs - REAL(i, dp)*exact)) > 1.e-12_dp) ERROR STOP 'Metric solution mismatch'
         IF (state%factorizations /= nfactor) ERROR STOP 'Unnecessary metric refactorization'
      END DO
      CALL sparse_inertia_finish(state)
      CALL sparse_metric_solve(state, rhs, para_env, info)
      IF (info /= -4) ERROR STOP 'Uninitialized metric factors accepted'
   END SUBROUTINE metric_tests

! **************************************************************************************************
!> \brief Compare distributed realified sparse inertia and gap with the complex dense oracle.
!> \param matrix Complex localizer
!> \param metric Physical metric (duplicated in the localizer pencil)
!> \param expected Expected resolved index
!> \param label Benchmark label
!> \param aii Test the time-reversal Z2 Pfaffian instead of the class A half-signature
!> \param expected_gap Optional analytic gap that must lie inside the returned bounds
! **************************************************************************************************
   SUBROUTINE compare(matrix, metric, expected, label, aii, expected_gap)
      COMPLEX(KIND=dp), INTENT(IN)                       :: matrix(:, :), metric(:, :)
      INTEGER, INTENT(IN)                                :: expected
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      LOGICAL, INTENT(IN), OPTIONAL                      :: aii
      REAL(dp), INTENT(IN), OPTIONAL                     :: expected_gap

      COMPLEX(KIND=dp)                                   :: av, bv
      INTEGER                                            :: c, col, dense_idx, dense_status, dim, &
                                                            idx, k, order, pass, physical, r, row, &
                                                            status
      INTEGER(int_8)                                     :: fill, memory_limit, peak_memory
      INTEGER(KIND=int_8)                                :: after, before, rate
      INTEGER, ALLOCATABLE                               :: cols(:), rows(:)
      LOGICAL                                            :: use_aii
      REAL(dp)                                           :: residual
      REAL(KIND=dp)                                      :: ar, br, dense_seconds, gap, lower, &
                                                            ratio, sparse_seconds, upper
      REAL(KIND=dp), ALLOCATABLE                         :: a(:), b(:), identity(:)
      TYPE(sparse_inertia_type)                          :: stats

      use_aii = .FALSE.
      IF (PRESENT(aii)) use_aii = aii
      dim = SIZE(matrix, 1)
      physical = dim/2
      order = 2*dim
      CALL SYSTEM_CLOCK(before, rate)
      IF (para_env%is_source()) THEN
         CALL localizer_matrix_evaluate(matrix, metric, use_aii, 1.e-10_dp, 1.e-10_dp, 1.e-8_dp, &
                                        dense_idx, gap, ratio, dense_status)
      END IF
      CALL para_env%bcast(dense_idx, para_env%source)
      CALL para_env%bcast(gap, para_env%source)
      CALL para_env%bcast(dense_status, para_env%source)
      CALL SYSTEM_CLOCK(after)
      dense_seconds = REAL(after - before, dp)/REAL(rate, dp)
      DO pass = 1, 2
         k = 0
         DO row = 1, order
            IF (MOD(row - 1, para_env%num_pe) /= para_env%mepos) CYCLE
            r = MOD(row - 1, dim) + 1
            DO col = 1, row
               c = MOD(col - 1, dim) + 1
               av = matrix(r, c)
               bv = 0.0_dp
               IF ((r - 1)/physical == (c - 1)/physical) THEN
                  bv = metric(MOD(r - 1, physical) + 1, MOD(c - 1, physical) + 1)
               END IF
               IF ((row - 1)/dim == (col - 1)/dim) THEN
                  ar = REAL(av, dp)
                  br = REAL(bv, dp)
               ELSE
                  ar = AIMAG(av)
                  br = AIMAG(bv)
               END IF
               IF (ar == 0.0_dp .AND. br == 0.0_dp .AND. row /= col) CYCLE
               k = k + 1
               IF (pass == 2) THEN
                  rows(k) = row
                  cols(k) = col
                  a(k) = ar
                  b(k) = br
                  identity(k) = 0.0_dp
                  IF (row == col) identity(k) = 1.0_dp
               END IF
            END DO
         END DO
         IF (pass == 1) ALLOCATE (rows(k), cols(k), a(k), b(k), identity(k))
      END DO
      CALL SYSTEM_CLOCK(before)
      CALL localizer_sparse_evaluate(order, rows, cols, a, b, identity, para_env, 1.e-10_dp, 1.e-8_dp, &
                                     idx, lower, upper, ratio, stats, status)
      IF (use_aii .AND. status == 0) THEN
         IF (idx /= 0) THEN
            status = -3
         ELSE
            IF (label == 'AII spin-mixed') THEN
               CALL sparse_z2_evaluate(order, rows, cols, a, b, para_env, 1.e-10_dp, ratio, lower, &
                                       idx, residual, fill, status, max_memory_bytes=1_int_8, &
                                       peak_memory_bytes=peak_memory)
               IF (status /= -6 .OR. idx /= 0 .OR. peak_memory /= 0) THEN
                  ERROR STOP 'Pfaffian memory failure not propagated collectively'
               END IF
               IF (para_env%num_pe > 1) THEN
                  memory_limit = INT(para_env%mepos + 1, int_8)
                  CALL sparse_z2_evaluate(order, rows, cols, a, b, para_env, 1.e-10_dp, ratio, lower, &
                                          idx, residual, fill, status, max_memory_bytes=memory_limit)
                  IF (status /= -1) ERROR STOP 'Inconsistent Pfaffian memory limits accepted'
               END IF
            END IF
            CALL sparse_z2_evaluate(order, rows, cols, a, b, para_env, 1.e-10_dp, ratio, lower, &
                                    idx, residual, fill, status)
            IF (para_env%is_source()) WRITE (*, '(A,ES14.5,I12)') TRIM(label)//' Pfaffian residual/fill ', residual, fill
         END IF
      END IF
      CALL SYSTEM_CLOCK(after)
      sparse_seconds = REAL(after - before, dp)/REAL(rate, dp)
      IF (status /= dense_status) THEN
         IF (para_env%is_source()) THEN
            WRITE (*, '(A,2I6,3ES16.7,I8)') TRIM(label)//' status(dense,sparse)/metric/gap/count ', &
               dense_status, status, ratio, lower, upper, stats%factorizations
            WRITE (*, '(A,2I12)') 'Sparse backend error code/detail: ', stats%error_code, stats%error_detail
         END IF
         ERROR STOP 'Sparse/dense resolution mismatch'
      END IF
      IF (status == 0) THEN
         IF (idx /= dense_idx .OR. idx /= expected) ERROR STOP 'Sparse/dense index mismatch'
         IF (gap < lower - 1.e-9_dp .OR. gap > upper + 1.e-9_dp) ERROR STOP 'Generalized gap bracket mismatch'
         IF (PRESENT(expected_gap)) THEN
            IF (expected_gap < lower .OR. expected_gap > upper) ERROR STOP 'Analytic gap outside shifted inertia bracket'
         END IF
      END IF
      IF (para_env%is_source()) THEN
         IF (stats%workspace_retries > 0) THEN
            WRITE (*, '(A,I6)') TRIM(label)//' workspace retries: ', stats%workspace_retries
         END IF
         WRITE (*, '(A,A,I8,2F12.5,3ES16.7,4I12,F12.3)') TRIM(label), &
            ' order/time(dense,sparse)/gap/nnz/factor/work/count/MB ', &
            order, dense_seconds, sparse_seconds, gap, lower, upper, stats%input_entries, &
            stats%factor_entries, stats%factor_workspace_entries, stats%factorizations, stats%peak_memory_mb
      END IF
   END SUBROUTINE compare

! **************************************************************************************************
!> \brief Compare nontrivial spin-mixed AII models and metric congruences with the dense oracle.
!> \param h Scalar QWZ Hamiltonian
!> \param local Scalar QWZ localizer
!> \param s Scalar metric
! **************************************************************************************************
   SUBROUTINE aii_tests(h, local, s)
      COMPLEX(dp), INTENT(IN)                            :: h(:, :), local(:, :), s(:, :)

      COMPLEX(dp), ALLOCATABLE                           :: a(:, :), change(:, :), spin_h(:, :), &
                                                            spin_s(:, :)
      INTEGER                                            :: i, n

      n = SIZE(h, 1)
      ALLOCATE (spin_h(2*n, 2*n), spin_s(2*n, 2*n), a(4*n, 4*n), change(2*n, 2*n))
      spin_h = 0.0_dp
      spin_s = 0.0_dp
      a = 0.0_dp
      spin_h(:n, :n) = h
      spin_h(n + 1:, n + 1:) = CONJG(h)
      spin_s(:n, :n) = s
      spin_s(n + 1:, n + 1:) = s
      DO i = 1, n, 2
         spin_h(i, n + i + 1) = CMPLX(0.1_dp, 0.05_dp, dp)
         spin_h(i + 1, n + i) = -spin_h(i, n + i + 1)
      END DO
      spin_h(n + 1:, :n) = CONJG(TRANSPOSE(spin_h(:n, n + 1:)))
      a(:2*n, :2*n) = spin_h
      a(2*n + 1:, 2*n + 1:) = -spin_h
      a(:n, 2*n + 1:3*n) = local(:n, n + 1:)
      a(n + 1:2*n, 3*n + 1:) = local(:n, n + 1:)
      a(2*n + 1:, :2*n) = CONJG(TRANSPOSE(a(:2*n, 2*n + 1:)))
      CALL compare(a, spin_s, 1, 'AII spin-mixed', .TRUE.)
      DO i = 1, 2*n, 2
         a(i, i) = a(i, i) - 2.0_dp
         a(i + 1, i + 1) = a(i + 1, i + 1) + 2.0_dp
      END DO
      a(2*n + 1:, 2*n + 1:) = -a(:2*n, :2*n)
      CALL compare(a, spin_s, 0, 'AII trivial QWZ pair', .TRUE.)
      a(:2*n, :2*n) = spin_h
      a(2*n + 1:, 2*n + 1:) = -spin_h
      change = 0.0_dp
      DO i = 1, n
         change(i, i) = CMPLX(0.8_dp + 0.02_dp*MOD(i, 5), 0.1_dp, dp)
         IF (i < n) change(i, i + 1) = CMPLX(0.02_dp, -0.01_dp, dp)
      END DO
      change(n + 1:, n + 1:) = CONJG(change(:n, :n))
      spin_s(:, :) = MATMUL(CONJG(TRANSPOSE(change)), MATMUL(spin_s, change))
      a(:2*n, :2*n) = MATMUL(CONJG(TRANSPOSE(change)), MATMUL(a(:2*n, :2*n), change))
      a(2*n + 1:, 2*n + 1:) = -a(:2*n, :2*n)
      a(:2*n, 2*n + 1:) = MATMUL(CONJG(TRANSPOSE(change)), MATMUL(a(:2*n, 2*n + 1:), change))
      a(2*n + 1:, :2*n) = CONJG(TRANSPOSE(a(:2*n, 2*n + 1:)))
      CALL compare(a, spin_s, 1, 'AII complex nonorthogonal', .TRUE.)
      spin_s(1, 1) = spin_s(1, 1) + 0.1_dp
      CALL compare(a, spin_s, 0, 'AII broken metric TR', .TRUE.)
      spin_s(1, 1) = spin_s(1, 1) - 0.1_dp
      a(1, 1) = a(1, 1) + 1.e-12_dp
      CALL compare(a, spin_s, 1, 'AII roundoff TR residual', .TRUE.)
      a(1, 1) = a(1, 1) + 0.1_dp
      CALL compare(a, spin_s, 0, 'AII broken TR', .TRUE.)
      a = 0.0_dp
      a(:2*n, :2*n) = spin_s
      a(2*n + 1:, 2*n + 1:) = -spin_s
      CALL compare(a, spin_s, 0, 'AII atomic reference', .TRUE.)
      CALL compare(1.e-10_dp*a, spin_s, 0, 'AII unresolved', .TRUE.)
      CALL compare(0.0_dp*a, spin_s, 0, 'AII singular', .TRUE.)
   END SUBROUTINE aii_tests

! **************************************************************************************************
!> \brief Exercise size parity, disconnected components, non-bipartite graphs and pivoting.
!> \param native Test native distributed factors instead of optional Tacho
! **************************************************************************************************
   SUBROUTINE pfaffian_tests(native)
      LOGICAL, INTENT(IN)                                :: native

      INTEGER                                            :: case_number, i, j, k, n, sign, status
      INTEGER(int_8)                                     :: fill, peak
      INTEGER, ALLOCATABLE                               :: cols(:), rows(:)
      REAL(dp)                                           :: residual
      REAL(dp), ALLOCATABLE                              :: a(:, :), blocks(:, :), lower(:, :), &
                                                            values(:)

      pf_native = native
      IF (.NOT. native .AND. .NOT. para_env%is_source()) RETURN
      DO n = 2, 24, 2
         ALLOCATE (a(n, n))
         DO case_number = 1, 4
            a = 0.0_dp
            DO j = 2, n
               DO i = 1, j - 1
                  IF (MOD(i*7 + j*11 + case_number, 3) /= 0) CYCLE
                  a(i, j) = SIN(REAL(i*17 + j*13 + case_number, dp))
               END DO
            END DO
            DO i = 1, n, 2
               a(i, i + 1) = a(i, i + 1) + 2.0_dp
            END DO
            a(:, :) = a - TRANSPOSE(a)
            CALL check_pfaffian(a)
            CALL check_pfaffian(1.e-120_dp*a)
            CALL check_pfaffian(1.e120_dp*a)
         END DO
         a = 0.0_dp
         DO i = 1, n/2
            a(i, n/2 + i) = REAL((-1)**i, dp)
            a(n/2 + i, i) = -a(i, n/2 + i)
         END DO
         CALL check_pfaffian(a)
         a(1, :) = 0.0_dp
         a(:, 1) = 0.0_dp
         CALL check_pfaffian(a)
         DEALLOCATE (a)
      END DO
      ! Exact rank-22 congruence with no zero rows. Probes b=A*x can hide its
      ! nullspace, whereas independent right-hand sides expose the failed solve.
      n = 24
      ALLOCATE (a(n, n), blocks(n, n), lower(n, n))
      lower = 0.0_dp
      blocks = 0.0_dp
      DO j = 1, n
         lower(j, j) = 1.0_dp
         DO i = j + 1, MIN(n, j + 4)
            lower(i, j) = REAL(MOD(17*i + 13*j + 3, 5) - 2, dp)/8.0_dp
         END DO
      END DO
      DO i = 1, n - 3, 2
         blocks(i, i + 1) = REAL((-1)**((i - 1)/2), dp)
         blocks(i + 1, i) = -blocks(i, i + 1)
      END DO
      a(:, :) = MATMUL(lower, MATMUL(blocks, TRANSPOSE(lower)))
      CALL check_pfaffian(a, expected_zero=.TRUE.)
      CALL check_pfaffian(1.e-120_dp*a, expected_zero=.TRUE.)
      CALL check_pfaffian(1.e120_dp*a, expected_zero=.TRUE.)
      DEALLOCATE (a, blocks, lower)
      ! Pf(A)=3 by exact recursive expansion. The fixed-front library cannot
      ! factor this connected graph even though its smallest singular value is 0.39.
      ALLOCATE (a(8, 8))
      a = 0.0_dp
      a(1, 2) = -1.0_dp
      a(1, 6) = 1.0_dp
      a(1, 8) = -1.0_dp
      a(2, 7) = -1.0_dp
      a(3, 4) = 1.0_dp
      a(3, 8) = -3.0_dp
      a(4, 5) = -1.0_dp
      a(5, 6) = 1.0_dp
      a(7, 8) = 1.0_dp
      a(:, :) = a - TRANSPOSE(a)
      CALL check_pfaffian(a)
      CALL check_pfaffian(1.e-120_dp*a)
      CALL check_pfaffian(1.e120_dp*a)
      CALL pfaffian_memory_test(a)
      DEALLOCATE (a)
      ! One-sided roundoff must not make the symbolic sparsity pattern asymmetric.
      ALLOCATE (a(32, 32))
      a = 0.0_dp
      DO i = 1, 31
         a(i, i + 1) = 1.0_dp
         a(i + 1, i) = -1.0_dp
      END DO
      a(2, 5) = 1.e-15_dp
      CALL check_pfaffian(a)
      CALL check_pfaffian(TRANSPOSE(a))
      CALL check_pfaffian(1.e-120_dp*a)
      CALL check_pfaffian(1.e120_dp*a)
      DEALLOCATE (a)
      ! Pf(A)=60 exactly. Integer abs in Tacho's pivot search misses nonzero
      ! Schur-complement candidates after the adapter normalizes the entries.
      ALLOCATE (a(8, 8))
      a = 0.0_dp
      a(1, 2:) = [5.0_dp, -4.0_dp, 4.0_dp, -3.0_dp, 1.0_dp, -4.0_dp, -4.0_dp]
      a(2, 3:) = [-3.0_dp, 5.0_dp, -5.0_dp, -4.0_dp, 4.0_dp, 2.0_dp]
      a(3, 4:) = [-2.0_dp, -3.0_dp, 2.0_dp, 2.0_dp, -4.0_dp]
      a(4, 5:) = [0.0_dp, 2.0_dp, 0.0_dp, 4.0_dp]
      a(5, 6:) = [-3.0_dp, -5.0_dp, -4.0_dp]
      a(6, 7:) = [-5.0_dp, 4.0_dp]
      a(7, 8) = -4.0_dp
      a(:, :) = a - TRANSPOSE(a)
      CALL localizer_pfaffian_sign(a, sign, status)
      IF (status /= 0 .OR. sign /= 1) ERROR STOP 'Dense pivot fixture failed'
      CALL check_pfaffian(a)
      CALL check_pfaffian(1.e-120_dp*a)
      CALL check_pfaffian(1.e120_dp*a)
      DEALLOCATE (a)
      ! A numerically singular matrix whose graph still has perfect matchings.
      ALLOCATE (a(4, 4))
      a = 0.0_dp
      a(1, 2) = 1.0_dp
      a(1, 3) = 1.0_dp
      a(2, 4) = 1.0_dp
      a(3, 4) = 1.0_dp
      a(:, :) = a - TRANSPOSE(a)
      CALL check_pfaffian(a)
      CALL evaluate_pfaffian(2, [1, 2], [2, 1], [1.0_dp, 1.0_dp], 1.e-10_dp, &
                             sign, residual, fill, status)
      IF (status /= -1 .OR. sign /= 0) ERROR STOP 'Non-skew Pfaffian input accepted'
      CALL evaluate_pfaffian(2, [1, 3], [2, 1], [1.0_dp, -1.0_dp], 1.e-10_dp, &
                             sign, residual, fill, status)
      IF (status /= -1 .OR. sign /= 0) ERROR STOP 'Invalid Pfaffian coordinates accepted'
      CALL evaluate_pfaffian(2, [1, 1, 2, 2], [2, 2, 1, 1], [2.0_dp, -1.0_dp, -2.0_dp, 1.0_dp], &
                             1.e-10_dp, sign, residual, fill, status)
      IF (status /= 0 .OR. sign /= 1) ERROR STOP 'Additive Pfaffian coordinates failed'
      CALL evaluate_pfaffian(4, [1, 2, 3, 4], [2, 1, 4, 3], [1.0_dp, -1.0_dp, 1.e-12_dp, -1.e-12_dp], &
                             1.e-10_dp, sign, residual, fill, status)
      IF (status /= 1 .OR. sign /= 0) ERROR STOP 'Unresolved Pfaffian pivot accepted'
      CALL evaluate_pfaffian(6, [1, 2, 3, 4, 5, 6, 3, 5, 3], [2, 1, 4, 3, 6, 5, 5, 3, 1], &
                             [1.0_dp, -1.0_dp, 2.0_dp, -2.0_dp, 3.0_dp, -3.0_dp, 0.5_dp, -0.5_dp, 1.e-12_dp], &
                             1.e-10_dp, sign, residual, fill, status)
      IF (status /= 0 .OR. sign /= 1) ERROR STOP 'One-sided roundoff between components failed'
      CALL evaluate_pfaffian(6, [1, 2, 3, 4, 5, 6, 3, 5, 3], [2, 1, 4, 3, 6, 5, 5, 3, 1], &
                             [1.0_dp, -1.0_dp, 2.0_dp, -2.0_dp, 3.0_dp, -3.0_dp, 0.5_dp, -0.5_dp, 1.e-6_dp], &
                             1.e-10_dp, sign, residual, fill, status)
      IF (status /= -1 .OR. sign /= 0) ERROR STOP 'Non-skew coupling between components accepted'
      IF (native) THEN
         ! The even nearest-neighbor chain has Pf(A)=1 by recursive expansion.
         ! Its dense matrix alone would need 8 MiB, exceeding this 2 MiB budget.
         n = 1024
         ALLOCATE (rows(2*(n - 1)), cols(2*(n - 1)), values(2*(n - 1)))
         DO i = 1, n - 1
            k = 2*i - 1
            rows(k:k + 1) = [i, i + 1]
            cols(k:k + 1) = [i + 1, i]
            values(k:k + 1) = [1.0_dp, -1.0_dp]
         END DO
         CALL evaluate_pfaffian(n, rows, cols, values, 1.e-10_dp, sign, residual, fill, status, &
                                max_memory_bytes=2_int_8*1024**2, peak_memory_bytes=peak)
         IF (status /= 0 .OR. sign /= 1 .OR. residual > 1.e-10_dp .OR. fill >= INT(n, int_8)**2) THEN
            ERROR STOP 'Distributed sparse chain Pfaffian failed'
         END IF
         IF (peak <= 0 .OR. peak > 2_int_8*1024**2) ERROR STOP 'Sparse chain numerical budget failed'
         IF (para_env%is_source()) WRITE (*, '(A,I0,ES14.5,2I12)') 'Native chain order/residual/fill/bytes ', &
            n, residual, fill, peak
         DEALLOCATE (rows, cols, values)
      END IF
      IF (para_env%is_source()) THEN
         WRITE (*, '(A,L2)') 'Sparse Pfaffian polynomial/permutation/scale/component tests passed, native=', native
      END IF
   END SUBROUTINE pfaffian_tests

! **************************************************************************************************
!> \brief Numerical buffer limits include both fronts and solve workspaces.
!> \param a Nonsingular connected fixture
! **************************************************************************************************
   SUBROUTINE pfaffian_memory_test(a)
      REAL(dp), INTENT(IN)                               :: a(:, :)

      INTEGER                                            :: i, j, k, n, sign, status
      INTEGER(int_8)                                     :: fill, limit, peak
      INTEGER, ALLOCATABLE                               :: cols(:), rows(:)
      REAL(dp)                                           :: residual
      REAL(dp), ALLOCATABLE                              :: values(:)

      n = SIZE(a, 1)
      ALLOCATE (rows(n*n), cols(n*n), values(n*n))
      k = 0
      DO j = 1, n
         DO i = 1, n
            k = k + 1
            rows(k) = i
            cols(k) = j
            values(k) = a(i, j)
         END DO
      END DO
      CALL evaluate_pfaffian(n, rows, cols, values, 1.e-10_dp, sign, residual, fill, status, &
                             peak_memory_bytes=peak)
      IF (status /= 0 .OR. peak <= 0_int_8) ERROR STOP 'Missing Pfaffian memory accounting'
      limit = peak
      CALL evaluate_pfaffian(n, rows, cols, values, 1.e-10_dp, sign, residual, fill, status, &
                             max_memory_bytes=limit, peak_memory_bytes=peak)
      IF (status /= 0 .OR. sign /= 1 .OR. peak /= limit) ERROR STOP 'Exact Pfaffian budget failed'
      CALL evaluate_pfaffian(n, rows, cols, values, 1.e-10_dp, sign, residual, fill, status, &
                             max_memory_bytes=limit - 1, peak_memory_bytes=peak)
      IF (status /= -6 .OR. sign /= 0 .OR. peak >= limit) ERROR STOP 'Pfaffian memory budget exceeded'
      CALL evaluate_pfaffian(n, rows, cols, values, 1.e-10_dp, sign, residual, fill, status, &
                             max_memory_bytes=0_int_8)
      IF (status /= -1 .OR. sign /= 0) ERROR STOP 'Invalid Pfaffian memory limit accepted'
   END SUBROUTINE pfaffian_memory_test

! **************************************************************************************************
!> \brief Compare sparse sign against the independent dense skew elimination, including odd swaps.
!> \param a Real skew matrix
!> \param expected_zero Override the dense sign for an analytically rank-deficient fixture
! **************************************************************************************************
   SUBROUTINE check_pfaffian(a, expected_zero)
      REAL(dp), INTENT(IN)                               :: a(:, :)
      LOGICAL, INTENT(IN), OPTIONAL                      :: expected_zero

      INTEGER                                            :: dense_sign, dense_status, i, j, k, n, &
                                                            pass, sign, status, swap
      INTEGER(int_8)                                     :: fill
      INTEGER, ALLOCATABLE                               :: cols(:), permutation(:), rows(:)
      REAL(dp)                                           :: residual
      REAL(dp), ALLOCATABLE                              :: values(:)

      n = SIZE(a, 1)
      ALLOCATE (permutation(n))
      permutation(:) = [(i, i=1, n)]
      CALL localizer_pfaffian_sign(a, dense_sign, dense_status)
      ! Scaling roundoff can leave a tiny nonzero last dense pivot. The known
      ! nullity of this fixture supplies the reference independently of it.
      IF (PRESENT(expected_zero)) THEN
         IF (expected_zero) THEN
            dense_sign = 0
            dense_status = -2
         END IF
      END IF
      DO swap = 1, 2
         DO pass = 1, 2
            k = 0
            DO j = 1, n
               DO i = 1, n
                  IF (a(permutation(i), permutation(j)) == 0.0_dp) CYCLE
                  k = k + 1
                  IF (pass == 2) THEN
                     rows(k) = i
                     cols(k) = j
                     values(k) = a(permutation(i), permutation(j))
                  END IF
               END DO
            END DO
            IF (pass == 1) ALLOCATE (rows(k), cols(k), values(k))
         END DO
         CALL evaluate_pfaffian(n, rows, cols, values, 1.e-10_dp, sign, residual, fill, status)
         IF ((dense_status == 0) .NEQV. (status == 0)) THEN
            WRITE (*, *) 'Pfaffian status: order, swap, dense, sparse, scale, residual = ', &
               n, swap, dense_status, status, MAXVAL(ABS(a)), residual
            ERROR STOP 'Sparse Pfaffian status mismatch'
         END IF
         IF (status == 0 .AND. sign /= dense_sign) ERROR STOP 'Sparse Pfaffian sign mismatch'
         IF (status /= 0 .AND. sign /= 0) ERROR STOP 'Invalid sparse Pfaffian sign returned'
         DEALLOCATE (rows, cols, values)
         permutation(:2) = [2, 1]
         dense_sign = -dense_sign
      END DO
   END SUBROUTINE check_pfaffian

! **************************************************************************************************
!> \brief Partition both triangles independently, including additive duplicates across MPI ranks.
!> \param order Matrix dimension
!> \param rows Full reference COO rows
!> \param cols Full reference COO columns
!> \param values Full reference COO values
!> \param tolerance Relative threshold
!> \param sign Pfaffian sign
!> \param residual Solve residual
!> \param fill Factor entries
!> \param status Solver status
!> \param max_memory_bytes Optional numerical buffer limit
!> \param peak_memory_bytes Optional accepted peak bytes
! **************************************************************************************************
   SUBROUTINE evaluate_pfaffian(order, rows, cols, values, tolerance, sign, residual, fill, status, &
                                max_memory_bytes, peak_memory_bytes)
      INTEGER, INTENT(IN)                                :: order, rows(:), cols(:)
      REAL(dp), INTENT(IN)                               :: values(:), tolerance
      INTEGER, INTENT(OUT)                               :: sign
      REAL(dp), INTENT(OUT)                              :: residual
      INTEGER(int_8), INTENT(OUT)                        :: fill
      INTEGER, INTENT(OUT)                               :: status
      INTEGER(int_8), INTENT(IN), OPTIONAL               :: max_memory_bytes
      INTEGER(int_8), INTENT(OUT), OPTIONAL              :: peak_memory_bytes

      INTEGER                                            :: i
      INTEGER, ALLOCATABLE                               :: owned(:)

      IF (pf_native) THEN
         owned = PACK([(i, i=1, SIZE(values))], &
                      [(MOD(i - 1, para_env%num_pe) == para_env%mepos, i=1, SIZE(values))])
         CALL sparse_pfaffian_sign(order, rows(owned), cols(owned), values(owned), tolerance, sign, residual, &
                                   fill, status, max_memory_bytes, peak_memory_bytes, para_env)
      ELSE
         CALL sparse_pfaffian_sign(order, rows, cols, values, tolerance, sign, residual, fill, status, &
                                   max_memory_bytes, peak_memory_bytes)
      END IF
   END SUBROUTINE evaluate_pfaffian

END PROGRAM spectral_localizer_sparse_unittest
