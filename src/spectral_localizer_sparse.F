!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Sparse class A index and metric-aware gap bounds from shifted inertia counts.
!> \author Thomas D. Kuehne, tkuehne@cp2k.org
! **************************************************************************************************
MODULE spectral_localizer_sparse
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE localizer_sparse_inertia,        ONLY: sparse_inertia_count,&
                                              sparse_inertia_finish,&
                                              sparse_inertia_start,&
                                              sparse_inertia_type
   USE message_passing,                 ONLY: mp_para_env_type

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: localizer_sparse_evaluate

CONTAINS

! **************************************************************************************************
!> \brief Analyze a realified Hermitian pencil without dense storage or a metric inverse square root.
!> \param order Realified order (twice the complex localizer dimension)
!> \param rows Distributed lower COO row indices, including identity entries
!> \param cols Distributed lower COO column indices
!> \param a Local covariant localizer values
!> \param b Local metric values on the union pattern
!> \param identity Local identity values on the same pattern
!> \param para_env Communicator
!> \param metric_tol Minimum numerical lower bound on the overlap condition ratio
!> \param gap_tol Absolute gap threshold
!> \param invariant Resolved half-signature of the complex localizer
!> \param gap_lower Lower bound on the nearest absolute generalized eigenvalue
!> \param gap_upper Upper bound on that eigenvalue
!> \param metric_ratio Conservative lower bound on lambda_min(B)/lambda_max(B)
!> \param statistics Factorization counts, workspace and input size (context is closed on return)
!> \param status 0 resolved, 1 unresolved, -1 input, -2 metric, -4 solver, -5 unavailable
! **************************************************************************************************
   SUBROUTINE localizer_sparse_evaluate(order, rows, cols, a, b, identity, para_env, metric_tol, gap_tol, &
                                        invariant, gap_lower, gap_upper, metric_ratio, statistics, status)
      INTEGER, INTENT(IN)                                :: order, rows(:), cols(:)
      REAL(KIND=dp), INTENT(IN)                          :: a(:), b(:), identity(:)
      TYPE(mp_para_env_type), POINTER                    :: para_env
      REAL(KIND=dp), INTENT(IN)                          :: metric_tol, gap_tol
      INTEGER, INTENT(OUT)                               :: invariant
      REAL(KIND=dp), INTENT(OUT)                         :: gap_lower, gap_upper, metric_ratio
      TYPE(sparse_inertia_type), INTENT(INOUT)           :: statistics
      INTEGER, INTENT(OUT)                               :: status

      INTEGER                                            :: i, info, invalid, iteration, &
                                                            metric_factorizations, metric_retries, &
                                                            negative, reference
      INTEGER(KIND=int_8)                                :: metric_entries, metric_workspace
      REAL(KIND=dp)                                      :: anorm, bnorm, metric_memory, metric_min, &
                                                            mid, resolution, threshold
      REAL(KIND=dp), ALLOCATABLE                         :: sums(:, :), values(:)

      invariant = 0
      gap_lower = 0.0_dp
      gap_upper = 0.0_dp
      metric_ratio = 0.0_dp
      metric_factorizations = 0
      metric_retries = 0
      metric_entries = 0
      metric_workspace = 0
      metric_memory = 0.0_dp
      invalid = 0
      IF (order < 4 .OR. MOD(order, 4) /= 0) invalid = 1
      IF (SIZE(rows) /= SIZE(cols) .OR. SIZE(rows) /= SIZE(a) .OR. &
          SIZE(rows) /= SIZE(b) .OR. SIZE(rows) /= SIZE(identity)) invalid = 1
      IF (.NOT. ALL(ieee_is_finite(a)) .OR. .NOT. ALL(ieee_is_finite(b))) invalid = 1
      IF (.NOT. ALL(ieee_is_finite(identity))) invalid = 1
      IF (.NOT. ALL(ieee_is_finite([metric_tol, gap_tol]))) invalid = 1
      IF (metric_tol <= 0.0_dp .OR. metric_tol >= 1.0_dp .OR. gap_tol < 0.0_dp) invalid = 1
      CALL para_env%max(invalid)
      status = -1
      IF (invalid /= 0) RETURN
      CALL sparse_inertia_start(statistics, order, rows, cols, para_env, status, positive_definite=.TRUE.)
      IF (status /= 0) THEN
         CALL sparse_inertia_finish(statistics)
         RETURN
      END IF
      ALLOCATE (sums(order, 2), values(SIZE(a)))
      sums = 0.0_dp
      DO i = 1, SIZE(a)
         sums(rows(i), :) = sums(rows(i), :) + [ABS(a(i)), ABS(b(i))]
         IF (rows(i) /= cols(i)) sums(cols(i), :) = sums(cols(i), :) + [ABS(a(i)), ABS(b(i))]
      END DO
      CALL para_env%sum(sums)
      anorm = MAXVAL(sums(:, 1))
      bnorm = MAXVAL(sums(:, 2))
      status = -2
      IF (bnorm <= TINY(1.0_dp) .OR. .NOT. ieee_is_finite(bnorm)) GOTO 900
      CALL sparse_inertia_count(statistics, b, para_env, negative, info)
      IF (info < 0) status = -4
      IF (info /= 0 .OR. negative /= 0) GOTO 900
      metric_min = 0.0_dp
      mid = 0.5_dp*bnorm
      ! A conservative bound needs only one successful shifted SPD factorization, not an eigenvalue.
      ! MUMPS uses non-pivoted LDL for SYM=1, so success alone does NOT imply positive definiteness.
      DO iteration = 1, 100
         mid = MAX(mid, metric_tol*bnorm*(1.0_dp + 16.0_dp*EPSILON(1.0_dp)))
         values = b - mid*identity
         CALL sparse_inertia_count(statistics, values, para_env, negative, info)
         IF (info < 0) THEN
            status = -4
            GOTO 900
         END IF
         IF (info == 0 .AND. negative == 0) THEN
            metric_min = mid
            EXIT
         END IF
         IF (mid <= metric_tol*bnorm*(1.0_dp + 16.0_dp*EPSILON(1.0_dp))) EXIT
         mid = 0.5_dp*mid
      END DO
      metric_ratio = metric_min/bnorm
      IF (metric_ratio <= metric_tol) GOTO 900
      metric_factorizations = statistics%factorizations
      metric_retries = statistics%workspace_retries
      metric_entries = statistics%factor_entries
      metric_workspace = statistics%factor_workspace_entries
      metric_memory = statistics%peak_memory_mb
      CALL sparse_inertia_finish(statistics)
      ! The localizer and its energy shifts are indefinite and must retain numerical pivoting.
      CALL sparse_inertia_start(statistics, order, rows, cols, para_env, status)
      IF (status /= 0) GOTO 900
      status = -1
      IF (.NOT. ieee_is_finite(anorm)) GOTO 900
      CALL sparse_inertia_count(statistics, a, para_env, reference, info)
      status = -4
      IF (info < 0) GOTO 900
      status = 1
      IF (info /= 0 .OR. MOD(reference, 2) /= 0) GOTO 900
      gap_upper = anorm/metric_min
      threshold = MAX(gap_tol, 100.0_dp*EPSILON(1.0_dp)*REAL(order, dp)*gap_upper)
      resolution = MAX(1.e-12_dp, 0.01_dp*gap_tol, 32.0_dp*EPSILON(1.0_dp)*gap_upper)
      DO iteration = 1, 100
         mid = 0.5_dp*(gap_lower + gap_upper)
         values = a - mid*b
         CALL sparse_inertia_count(statistics, values, para_env, negative, info)
         IF (info < 0) THEN
            status = -4
            GOTO 900
         END IF
         IF (info /= 0 .OR. negative /= reference) THEN
            gap_upper = mid
         ELSE
            values = a + mid*b
            CALL sparse_inertia_count(statistics, values, para_env, negative, info)
            IF (info < 0) THEN
               status = -4
               GOTO 900
            END IF
            IF (info /= 0 .OR. negative /= reference) THEN
               gap_upper = mid
            ELSE
               gap_lower = mid
            END IF
         END IF
         IF (gap_upper - gap_lower <= resolution) EXIT
      END DO
      IF (gap_upper - gap_lower > resolution .OR. gap_lower <= threshold) GOTO 900
      invariant = (order - 2*reference)/4
      status = 0
900   CONTINUE
      CALL sparse_inertia_finish(statistics)
      statistics%factorizations = statistics%factorizations + metric_factorizations
      statistics%workspace_retries = statistics%workspace_retries + metric_retries
      statistics%factor_entries = MAX(statistics%factor_entries, metric_entries)
      statistics%factor_workspace_entries = MAX(statistics%factor_workspace_entries, metric_workspace)
      statistics%peak_memory_mb = MAX(statistics%peak_memory_mb, metric_memory)
   END SUBROUTINE localizer_sparse_evaluate

END MODULE spectral_localizer_sparse
