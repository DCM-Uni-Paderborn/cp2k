!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Optional distributed sparse LDL inertia with reusable symbolic analysis.
!> \author Thomas D. Kuehne, tkuehne@cp2k.org
! **************************************************************************************************
MODULE localizer_sparse_inertia
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE message_passing,                 ONLY: mp_para_env_type

   IMPLICIT NONE
   PRIVATE
#if defined(__MUMPS)
   INCLUDE 'dmumps_struc.h'
#endif
   PUBLIC :: sparse_inertia_type, sparse_inertia_available, sparse_inertia_start, &
             sparse_inertia_count, sparse_inertia_finish
   PUBLIC :: sparse_metric_start, sparse_metric_solve

   TYPE sparse_inertia_type
#if defined(__MUMPS)
      TYPE(DMUMPS_STRUC) :: solver
#endif
      LOGICAL :: initialized = .FALSE.
      LOGICAL :: factor_valid = .FALSE.
      INTEGER :: order = 0, factorizations = 0, error_code = 0, error_detail = 0, workspace_retries = 0
      INTEGER(KIND=int_8) :: factor_entries = 0, input_entries = 0, factor_workspace_entries = 0
      REAL(KIND=dp) :: peak_memory_mb = 0.0_dp
   END TYPE sparse_inertia_type

CONTAINS

! **************************************************************************************************
!> \brief Whether the optional sparse inertia backend is linked.
!> \return Backend availability
! **************************************************************************************************
   LOGICAL FUNCTION sparse_inertia_available()
      sparse_inertia_available = .FALSE.
#if defined(__MUMPS)
      sparse_inertia_available = .TRUE.
#endif
   END FUNCTION sparse_inertia_available

! **************************************************************************************************
!> \brief Analyze one distributed, lower-triangular coordinate pattern. Duplicates are summed.
!> \param state Caller-owned factorization context
!> \param order Real matrix order
!> \param rows Local one-based row coordinates
!> \param cols Local one-based column coordinates
!> \param para_env Communicator for collective calls
!> \param status Zero on success; negative for invalid input or unavailable backend
!> \param positive_definite Use the non-pivoting SPD path, with negative-pivot checks still required
! **************************************************************************************************
   SUBROUTINE sparse_inertia_start(state, order, rows, cols, para_env, status, positive_definite)
      TYPE(sparse_inertia_type), INTENT(INOUT) :: state
      INTEGER, INTENT(IN) :: order, rows(:), cols(:)
      TYPE(mp_para_env_type), POINTER :: para_env
      INTEGER, INTENT(OUT) :: status
      LOGICAL, INTENT(IN), OPTIONAL :: positive_definite
#if defined(__MUMPS)
      INTEGER :: invalid

      invalid = 0
      IF (state%initialized .OR. order < 1 .OR. SIZE(rows) /= SIZE(cols)) invalid = 1
      IF (ANY(rows < 1) .OR. ANY(rows > order) .OR. ANY(cols < 1)) invalid = 1
      IF (SIZE(rows) == SIZE(cols)) THEN
         IF (ANY(cols > rows)) invalid = 1
      END IF
      CALL para_env%max(invalid)
      status = -1
      IF (invalid /= 0) RETURN
      state%solver%COMM = para_env%get_handle()
      state%solver%PAR = 1
      state%solver%SYM = 2
      IF (PRESENT(positive_definite)) THEN
         IF (positive_definite) state%solver%SYM = 1
      END IF
      state%solver%JOB = -1
      CALL dmumps(state%solver)
      state%error_code = state%solver%INFOG(1)
      state%error_detail = state%solver%INFOG(2)
      status = -4
      IF (state%error_code /= 0) RETURN
      state%initialized = .TRUE.
      state%order = order
      state%factorizations = 0
      state%workspace_retries = 0
      state%factor_entries = 0
      state%factor_workspace_entries = 0
      state%peak_memory_mb = 0.0_dp
      state%input_entries = INT(SIZE(rows), int_8)
      CALL para_env%sum(state%input_entries)
      state%solver%N = order
      state%solver%NZ_LOC = SIZE(rows)
      ALLOCATE (state%solver%IRN_LOC(SIZE(rows)), state%solver%JCN_LOC(SIZE(rows)), &
                state%solver%A_LOC(SIZE(rows)))
      state%solver%IRN_LOC = rows
      state%solver%JCN_LOC = cols
      state%solver%A_LOC = 0.0_dp
      state%solver%ICNTL(1:3) = -1
      state%solver%ICNTL(4) = 0
      state%solver%ICNTL(6) = 0
      state%solver%ICNTL(8) = 0
      ! The ScaLAPACK root-node LU does not report its inertia: force root LDL, retain other MPI work.
      state%solver%ICNTL(13) = 1
      state%solver%ICNTL(18) = 3
      state%solver%ICNTL(24) = 0
      state%solver%ICNTL(35) = 0
      state%solver%CNTL(1) = 0.01_dp
      state%solver%CNTL(4) = -1.0_dp
      state%solver%JOB = 1
      CALL dmumps(state%solver)
      state%error_code = state%solver%INFOG(1)
      state%error_detail = state%solver%INFOG(2)
      status = 0
      IF (state%error_code /= 0) status = -4
#else
      status = -5
#endif
   END SUBROUTINE sparse_inertia_start

! **************************************************************************************************
!> \brief Numerically factor the same pattern, with no static pivot perturbation or compression.
!> \param state Analyzed factorization context
!> \param values Local lower-triangle values
!> \param para_env Communicator for collective calls
!> \param negative Number of negative eigenvalues of the real matrix
!> \param status 0 success, 1 singular, -4 numerical failure, -5 backend unavailable
! **************************************************************************************************
   SUBROUTINE sparse_inertia_count(state, values, para_env, negative, status)
      TYPE(sparse_inertia_type), INTENT(INOUT) :: state
      REAL(KIND=dp), INTENT(IN) :: values(:)
      TYPE(mp_para_env_type), POINTER :: para_env
      INTEGER, INTENT(OUT) :: negative, status
#if defined(__MUMPS)
      INTEGER :: invalid, entries, attempt, i
      INTEGER, ALLOCATABLE :: row_support(:)

      negative = 0
      state%factor_valid = .FALSE.
      invalid = 0
      IF (.NOT. state%initialized) invalid = 1
      IF (state%initialized) THEN
         IF (SIZE(values) /= state%solver%NZ_LOC) invalid = 1
      END IF
      IF (.NOT. ALL(ieee_is_finite(values))) invalid = 1
      CALL para_env%max(invalid)
      status = -4
      IF (invalid /= 0) RETURN
      ! At an exact gap endpoint a shifted pencil may contain whole zero rows.
      ! Detect these before delayed pivots consume MUMPS's workspace. Additive
      ! duplicate cancellation may hide other zero rows, but cannot give a false one.
      ALLOCATE (row_support(state%order))
      row_support = 0
      DO i = 1, SIZE(values)
         IF (values(i) == 0.0_dp) CYCLE
         row_support(state%solver%IRN_LOC(i)) = 1
         row_support(state%solver%JCN_LOC(i)) = 1
      END DO
      CALL para_env%max(row_support)
      IF (ANY(row_support == 0)) THEN
         status = 1
         RETURN
      END IF
      state%solver%A_LOC = values
      state%solver%JOB = 2
      DO attempt = 1, 6
         CALL dmumps(state%solver)
         state%error_code = state%solver%INFOG(1)
         state%error_detail = state%solver%INFOG(2)
         state%factorizations = state%factorizations + 1
         ! Numerical fill can exceed the symbolic estimate; MUMPS permits retrying JOB=2.
         IF (.NOT. ANY(state%error_code == [-8, -9, -17, -20])) EXIT
         IF (attempt == 6 .OR. state%solver%ICNTL(14) >= 640) EXIT
         state%solver%ICNTL(14) = MIN(640, MAX(40, 2*state%solver%ICNTL(14)))
         state%workspace_retries = state%workspace_retries + 1
      END DO
      IF (state%solver%INFOG(1) == -10) THEN
         status = 1
         RETURN
      END IF
      IF (state%solver%INFOG(1) /= 0) RETURN
      negative = state%solver%INFOG(12)
      entries = state%solver%INFOG(9)
      IF (entries >= 0) THEN
         state%factor_workspace_entries = MAX(state%factor_workspace_entries, INT(entries, int_8))
      ELSE
         state%factor_workspace_entries = MAX(state%factor_workspace_entries, -INT(entries, int_8)*1000000_int_8)
      END IF
      entries = state%solver%INFOG(29)
      IF (entries >= 0) THEN
         state%factor_entries = MAX(state%factor_entries, INT(entries, int_8))
      ELSE
         state%factor_entries = MAX(state%factor_entries, -INT(entries, int_8)*1000000_int_8)
      END IF
      state%peak_memory_mb = MAX(state%peak_memory_mb, REAL(state%solver%INFOG(22), dp))
      state%factor_valid = .TRUE.
      status = 0
#else
      negative = 0
      status = -5
#endif
   END SUBROUTINE sparse_inertia_count

! **************************************************************************************************
!> \brief Validate a sparse SPD metric and retain its factors for repeated right hand sides.
!> \param state New factorization context
!> \param order Matrix dimension
!> \param rows Distributed lower COO rows
!> \param cols Distributed lower COO columns
!> \param values Additive metric entries
!> \param identity Identity on the same COO pattern
!> \param tolerance Lower condition-ratio threshold
!> \param para_env Communicator
!> \param ratio Conservative lower eigenvalue/row-norm ratio
!> \param status 0 success, -2 invalid metric, -4 solver failure, -5 unavailable
! **************************************************************************************************
   SUBROUTINE sparse_metric_start(state, order, rows, cols, values, identity, tolerance, para_env, ratio, status)
      TYPE(sparse_inertia_type), INTENT(INOUT)           :: state
      INTEGER, INTENT(IN)                                :: order, rows(:), cols(:)
      REAL(dp), INTENT(IN)                               :: values(:), identity(:), tolerance
      TYPE(mp_para_env_type), POINTER                    :: para_env
      REAL(dp), INTENT(OUT)                              :: ratio
      INTEGER, INTENT(OUT)                               :: status

      INTEGER                                            :: i, info, invalid, negative
      REAL(dp)                                           :: lower, upper
      REAL(dp), ALLOCATABLE                              :: shifted(:), sums(:)

      ratio = 0.0_dp
      invalid = 0
      IF (SIZE(values) /= SIZE(rows) .OR. SIZE(identity) /= SIZE(rows)) invalid = 1
      IF (.NOT. ieee_is_finite(tolerance)) invalid = 1
      IF (tolerance <= 0.0_dp .OR. tolerance >= 1.0_dp) invalid = 1
      IF (.NOT. ALL(ieee_is_finite(values)) .OR. .NOT. ALL(ieee_is_finite(identity))) invalid = 1
      CALL para_env%max(invalid)
      status = -2
      IF (invalid /= 0) RETURN
      CALL sparse_inertia_start(state, order, rows, cols, para_env, status, positive_definite=.TRUE.)
      IF (status /= 0) RETURN
      ALLOCATE (sums(order), shifted(SIZE(values)))
      sums = 0.0_dp
      DO i = 1, SIZE(values)
         sums(rows(i)) = sums(rows(i)) + ABS(values(i))
         IF (rows(i) /= cols(i)) sums(cols(i)) = sums(cols(i)) + ABS(values(i))
      END DO
      CALL para_env%sum(sums)
      upper = MAXVAL(sums)
      status = -2
      IF (.NOT. ieee_is_finite(upper) .OR. upper <= TINY(1.0_dp)) RETURN
      CALL sparse_inertia_count(state, values, para_env, negative, info)
      IF (info < 0) status = info
      IF (info /= 0 .OR. negative /= 0) RETURN
      lower = 0.5_dp*upper
      DO i = 1, 100
         lower = MAX(lower, tolerance*upper*(1.0_dp + 16.0_dp*EPSILON(lower)))
         shifted = values - lower*identity
         CALL sparse_inertia_count(state, shifted, para_env, negative, info)
         IF (info < 0) THEN
            status = info
            RETURN
         END IF
         IF (info == 0 .AND. negative == 0) THEN
            ratio = lower/upper
            EXIT
         END IF
         IF (lower <= tolerance*upper*(1.0_dp + 16.0_dp*EPSILON(lower))) RETURN
         lower = lower*0.5_dp
      END DO
      IF (ratio <= tolerance) RETURN
      CALL sparse_inertia_count(state, values, para_env, negative, status)
      IF (status == 0 .AND. negative /= 0) status = -2
   END SUBROUTINE sparse_metric_start

! **************************************************************************************************
!> \brief Collective repeated solves using existing MUMPS factors, without refactorization.
!> \param state Successfully factored metric
!> \param rhs Replicated right hand sides, replaced by replicated solutions
!> \param para_env Communicator (MUMPS host is communicator rank zero)
!> \param status 0 success, -4 failure, -5 unavailable
! **************************************************************************************************
   SUBROUTINE sparse_metric_solve(state, rhs, para_env, status)
      TYPE(sparse_inertia_type), INTENT(INOUT) :: state
      REAL(dp), CONTIGUOUS, INTENT(INOUT) :: rhs(:, :)
      TYPE(mp_para_env_type), POINTER :: para_env
      INTEGER, INTENT(OUT) :: status
#if defined(__MUMPS)
      INTEGER :: invalid, count
      invalid = 0
      IF (.NOT. state%factor_valid .OR. .NOT. state%initialized) invalid = 1
      IF (SIZE(rhs, 1) /= state%order .OR. SIZE(rhs, 2) < 1) invalid = 1
      IF (SIZE(rhs, KIND=int_8) > HUGE(count)) invalid = 1
      IF (.NOT. ALL(ieee_is_finite(rhs))) invalid = 1
      CALL para_env%max(invalid)
      status = -4
      IF (invalid /= 0) RETURN
      state%solver%NRHS = SIZE(rhs, 2)
      state%solver%LRHS = state%order
      state%solver%ICNTL(20) = 0
      state%solver%ICNTL(21) = 0
      IF (para_env%mepos == 0) THEN
         ALLOCATE (state%solver%RHS(SIZE(rhs)))
         state%solver%RHS = RESHAPE(rhs, [SIZE(rhs)])
      END IF
      state%solver%JOB = 3
      CALL dmumps(state%solver)
      state%error_code = state%solver%INFOG(1)
      state%error_detail = state%solver%INFOG(2)
      IF (para_env%mepos == 0) THEN
         rhs = RESHAPE(state%solver%RHS, SHAPE(rhs))
         DEALLOCATE (state%solver%RHS)
         NULLIFY (state%solver%RHS)
      END IF
      IF (state%error_code /= 0) RETURN
      CALL para_env%bcast(rhs, 0)
      IF (.NOT. ALL(ieee_is_finite(rhs))) RETURN
      status = 0
#else
      status = -5
#endif
   END SUBROUTINE sparse_metric_solve

! **************************************************************************************************
!> \brief Release factors and caller-owned coordinate buffers collectively.
!> \param state Factorization context, including a failed analysis
! **************************************************************************************************
   SUBROUTINE sparse_inertia_finish(state)
      TYPE(sparse_inertia_type), INTENT(INOUT) :: state
#if defined(__MUMPS)
      IF (.NOT. state%initialized) RETURN
      DEALLOCATE (state%solver%IRN_LOC, state%solver%JCN_LOC, state%solver%A_LOC)
      NULLIFY (state%solver%IRN_LOC, state%solver%JCN_LOC, state%solver%A_LOC)
      state%solver%JOB = -2
      CALL dmumps(state%solver)
#endif
      state%initialized = .FALSE.
      state%factor_valid = .FALSE.
   END SUBROUTINE sparse_inertia_finish

END MODULE localizer_sparse_inertia
