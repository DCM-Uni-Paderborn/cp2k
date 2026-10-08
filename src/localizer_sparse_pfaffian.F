!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Sparse Pfaffian evaluation using distributed DBCSR factors or optional Tacho.
!> \author Thomas D. Kuehne, tkuehne@cp2k.org
! **************************************************************************************************
MODULE localizer_sparse_pfaffian
#if defined(__HAS_IEEE_EXCEPTIONS)
   USE ieee_exceptions, ONLY: IEEE_ALL, ieee_get_halting_mode, ieee_set_halting_mode
#endif
   USE iso_c_binding, ONLY: c_double, &
                            c_int, &
                            c_long_long
   USE kinds, ONLY: dp, &
                    int_8
   USE ieee_arithmetic, ONLY: ieee_is_finite
   USE localizer_dbcsr_factor, ONLY: localizer_coo_norm, localizer_dbcsr_pfaffian
   USE message_passing, ONLY: mp_para_env_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: sparse_pfaffian_available, sparse_z2_evaluate, sparse_pfaffian_sign

#if defined(__TACHO)
   INTERFACE
      FUNCTION tacho_z2(order, nnz, rows, cols, a, b, tolerance, metric_ratio, gap, index, residual, fill, &
                        max_memory, peak_memory) RESULT(status) BIND(C, NAME='cp2k_tacho_z2_limited')
         IMPORT :: c_int, c_double, c_long_long
         INTEGER(c_int), VALUE :: order, nnz
         INTEGER(c_int), INTENT(IN) :: rows(*), cols(*)
         REAL(c_double), INTENT(IN) :: a(*), b(*)
         REAL(c_double), VALUE :: tolerance, metric_ratio, gap
         INTEGER(c_int), INTENT(OUT) :: index
         REAL(c_double), INTENT(OUT) :: residual
         INTEGER(c_long_long), INTENT(OUT) :: fill
         INTEGER(c_long_long), VALUE :: max_memory
         INTEGER(c_long_long), INTENT(OUT) :: peak_memory
         INTEGER(c_int) :: status
      END FUNCTION tacho_z2
      FUNCTION tacho_sign(n, nnz, rows, cols, values, tolerance, sign, residual, fill, max_memory, peak_memory) &
         RESULT(status) BIND(C, NAME='cp2k_tacho_pfaffian_limited')
         IMPORT :: c_int, c_double, c_long_long
         INTEGER(c_int), VALUE :: n, nnz
         INTEGER(c_int), INTENT(IN) :: rows(*), cols(*)
         REAL(c_double), INTENT(IN) :: values(*)
         REAL(c_double), VALUE :: tolerance
         INTEGER(c_int), INTENT(OUT) :: sign
         REAL(c_double), INTENT(OUT) :: residual
         INTEGER(c_long_long), INTENT(OUT) :: fill
         INTEGER(c_long_long), VALUE :: max_memory
         INTEGER(c_long_long), INTENT(OUT) :: peak_memory
         INTEGER(c_int) :: status
      END FUNCTION tacho_sign
   END INTERFACE
#endif

CONTAINS

! **************************************************************************************************
!> \brief Whether the optional Tacho skew-LDL backend is linked.
!> \return Backend availability
! **************************************************************************************************
   LOGICAL FUNCTION sparse_pfaffian_available()
      sparse_pfaffian_available = .FALSE.
#if defined(__TACHO)
      sparse_pfaffian_available = .TRUE.
#endif
   END FUNCTION sparse_pfaffian_available

! **************************************************************************************************
!> \brief Evaluate an unreplicated real skew COO matrix, primarily for independent solver tests.
!> \param order Matrix order
!> \param rows One-based COO rows, both triangles
!> \param cols One-based COO columns
!> \param values Additive COO entries
!> \param tolerance Relative factor residual and numerical pivot threshold
!> \param sign Pfaffian sign, zero on failure
!> \param residual Largest solve-probe relative residual
!> \param fill Sum of stored scalar panel entries over disconnected components
!> \param status 0 success, 1 unresolved, -1 input, -4 solver, -5 unavailable, -6 memory
!> \param max_memory_bytes Numerical factor/front/update/solve buffer limit; default 4 GiB
!> \param peak_memory_bytes Peak requested numerical buffer bytes (not whole-process memory)
!> \param para_env Use the native distributed DBCSR factorization when supplied
! **************************************************************************************************
   SUBROUTINE sparse_pfaffian_sign(order, rows, cols, values, tolerance, sign, residual, fill, status, &
                                   max_memory_bytes, peak_memory_bytes, para_env)
      INTEGER, INTENT(IN)                                :: order
      INTEGER, CONTIGUOUS, INTENT(IN)                    :: rows(:), cols(:)
      REAL(dp), CONTIGUOUS, INTENT(IN)                   :: values(:)
      REAL(dp), INTENT(IN)                               :: tolerance
      INTEGER, INTENT(OUT)                               :: sign
      REAL(dp), INTENT(OUT)                              :: residual
      INTEGER(int_8), INTENT(OUT)                        :: fill
      INTEGER, INTENT(OUT)                               :: status
      INTEGER(int_8), INTENT(IN), OPTIONAL                :: max_memory_bytes
      INTEGER(int_8), INTENT(OUT), OPTIONAL               :: peak_memory_bytes
      TYPE(mp_para_env_type), POINTER, OPTIONAL           :: para_env
      INTEGER(int_8)                                     :: limit, peak
#if defined(__TACHO) && defined(__HAS_IEEE_EXCEPTIONS)
      LOGICAL :: halt(SIZE(IEEE_ALL))
#endif

      sign = 0
      residual = 0.0_dp
      fill = 0
      status = -1
      limit = 4096_int_8*1024_int_8**2
      IF (PRESENT(max_memory_bytes)) limit = max_memory_bytes
      peak = 0_int_8
      IF (PRESENT(peak_memory_bytes)) peak_memory_bytes = peak
      IF (PRESENT(para_env)) THEN
         CALL localizer_dbcsr_pfaffian(order, rows, cols, values, para_env, tolerance, sign, residual, &
                                       fill, status, limit, peak)
         IF (PRESENT(peak_memory_bytes)) peak_memory_bytes = peak
         RETURN
      END IF
      IF (SIZE(rows) /= SIZE(cols) .OR. SIZE(rows) /= SIZE(values)) RETURN
      IF (limit <= 0_int_8) RETURN
#if defined(__TACHO)
#if defined(__HAS_IEEE_EXCEPTIONS)
      CALL ieee_get_halting_mode(IEEE_ALL, halt)
      CALL ieee_set_halting_mode(IEEE_ALL, .FALSE.)
#endif
      status = tacho_sign(order, SIZE(rows), rows, cols, values, tolerance, sign, residual, fill, limit, peak)
#if defined(__HAS_IEEE_EXCEPTIONS)
      CALL ieee_set_halting_mode(IEEE_ALL, halt)
#endif
#else
      MARK_USED(order)
      MARK_USED(tolerance)
      status = -5
#endif
      IF (PRESENT(peak_memory_bytes)) peak_memory_bytes = peak
   END SUBROUTINE sparse_pfaffian_sign

! **************************************************************************************************
!> \brief Resolve Z2 after the distributed inertia path has verified a positive metric and gap.
!> \param order Realified pencil order, eight times the scalar AO dimension
!> \param rows Distributed lower realified COO rows
!> \param cols Distributed lower realified COO columns
!> \param a Localizer values
!> \param b Metric values on the same union pattern
!> \param para_env Communicator
!> \param tolerance Relative structure and factor residual tolerance
!> \param metric_ratio Conservative numerical metric condition lower bound
!> \param gap Positive lower gap bracket from the independent inertia calculation
!> \param invariant Z2 index (valid only on success)
!> \param residual Largest factor solve-probe relative residual
!> \param fill Stored scalar panel entries (not total memory consumption)
!> \param status 0 success, 1 unresolved, -1 input, -3 structure, -4 solver, -5 unavailable, -6 memory
!> \param max_memory_bytes Numerical factor/front/update/solve buffer limit; default 4 GiB
!> \param peak_memory_bytes Peak requested numerical buffer bytes (not whole-process memory)
!> \param use_native Use distributed DBCSR skew factors (default true)
! **************************************************************************************************
   SUBROUTINE sparse_z2_evaluate(order, rows, cols, a, b, para_env, tolerance, metric_ratio, gap, &
                                 invariant, residual, fill, status, max_memory_bytes, peak_memory_bytes, use_native)
      INTEGER, INTENT(IN) :: order
      INTEGER, CONTIGUOUS, INTENT(IN) :: rows(:), cols(:)
      REAL(dp), CONTIGUOUS, INTENT(IN) :: a(:), b(:)
      REAL(dp), INTENT(IN) :: tolerance, metric_ratio, gap
      TYPE(mp_para_env_type), POINTER :: para_env
      INTEGER, INTENT(OUT) :: invariant, status
      REAL(dp), INTENT(OUT) :: residual
      INTEGER(int_8), INTENT(OUT) :: fill
      INTEGER(int_8), INTENT(IN), OPTIONAL :: max_memory_bytes
      INTEGER(int_8), INTENT(OUT), OPTIONAL :: peak_memory_bytes
      LOGICAL, INTENT(IN), OPTIONAL :: use_native
      INTEGER(int_8) :: limit, peak
      LOGICAL :: native
      INTEGER :: selection(2)
#if defined(__TACHO)
      INTEGER :: allocation_status, i, invalid, size_global
      INTEGER(int_8) :: limit_max, limit_min, total
      INTEGER, ALLOCATABLE :: counts(:), offsets(:), global_rows(:), global_cols(:)
      REAL(dp), ALLOCATABLE :: global_a(:), global_b(:)
#if defined(__HAS_IEEE_EXCEPTIONS)
      LOGICAL :: halt(SIZE(IEEE_ALL))
#endif
#endif

      invariant = 0
      residual = 0.0_dp
      fill = 0
      status = -5
      limit = 4096_int_8*1024_int_8**2
      IF (PRESENT(max_memory_bytes)) limit = max_memory_bytes
      peak = 0_int_8
      IF (PRESENT(peak_memory_bytes)) peak_memory_bytes = peak
      native = .TRUE.
      IF (PRESENT(use_native)) native = use_native
      selection = MERGE(1, 0, native)
      CALL para_env%min(selection(1))
      CALL para_env%max(selection(2))
      IF (selection(1) /= selection(2)) THEN
         status = -1
         RETURN
      END IF
      IF (native) THEN
         CALL native_z2_evaluate(order, rows, cols, a, b, para_env, tolerance, metric_ratio, gap, &
                                 invariant, residual, fill, status, limit, peak)
         IF (PRESENT(peak_memory_bytes)) peak_memory_bytes = peak
         RETURN
      END IF
#if defined(__TACHO)
      invalid = 0
      IF (SIZE(rows) /= SIZE(cols) .OR. SIZE(rows) /= SIZE(a) .OR. SIZE(rows) /= SIZE(b)) invalid = 1
      IF (order < 8 .OR. MOD(order, 8) /= 0 .OR. gap <= 0.0_dp .OR. metric_ratio <= 0.0_dp) invalid = 1
      limit_min = limit
      limit_max = limit
      CALL para_env%min(limit_min)
      CALL para_env%max(limit_max)
      IF (limit_min <= 0_int_8 .OR. limit_min /= limit_max) invalid = 1
      CALL para_env%max(invalid)
      status = -1
      IF (invalid /= 0) RETURN
      ALLOCATE (counts(para_env%num_pe), STAT=allocation_status)
      IF (allocation_status == 0) ALLOCATE (offsets(para_env%num_pe), STAT=allocation_status)
      invalid = MERGE(1, 0, allocation_status /= 0)
      CALL para_env%max(invalid)
      IF (invalid /= 0) THEN
         status = -6
         RETURN
      END IF
      counts = 0
      counts(para_env%mepos + 1) = SIZE(a)
      CALL para_env%sum(counts)
      total = SUM(INT(counts, int_8))
      IF (total > HUGE(size_global)) RETURN
      offsets(1) = 0
      DO i = 2, SIZE(offsets)
         offsets(i) = offsets(i - 1) + counts(i - 1)
      END DO
      size_global = 0
      IF (para_env%is_source()) size_global = INT(total)
      ALLOCATE (global_rows(size_global), STAT=allocation_status)
      IF (allocation_status == 0) ALLOCATE (global_cols(size_global), STAT=allocation_status)
      IF (allocation_status == 0) ALLOCATE (global_a(size_global), STAT=allocation_status)
      IF (allocation_status == 0) ALLOCATE (global_b(size_global), STAT=allocation_status)
      invalid = MERGE(1, 0, allocation_status /= 0)
      CALL para_env%max(invalid)
      IF (invalid /= 0) THEN
         status = -6
         RETURN
      END IF
      CALL para_env%gatherv(rows, global_rows, counts, offsets, para_env%source)
      CALL para_env%gatherv(cols, global_cols, counts, offsets, para_env%source)
      CALL para_env%gatherv(a, global_a, counts, offsets, para_env%source)
      CALL para_env%gatherv(b, global_b, counts, offsets, para_env%source)
      IF (para_env%is_source()) THEN
#if defined(__HAS_IEEE_EXCEPTIONS)
         CALL ieee_get_halting_mode(IEEE_ALL, halt)
         CALL ieee_set_halting_mode(IEEE_ALL, .FALSE.)
#endif
         status = tacho_z2(order, size_global, global_rows, global_cols, global_a, global_b, &
                           tolerance, metric_ratio, gap, invariant, residual, fill, limit, peak)
#if defined(__HAS_IEEE_EXCEPTIONS)
         CALL ieee_set_halting_mode(IEEE_ALL, halt)
#endif
      END IF
      CALL para_env%bcast(status, para_env%source)
      CALL para_env%bcast(invariant, para_env%source)
      CALL para_env%bcast(residual, para_env%source)
      CALL para_env%bcast(fill, para_env%source)
      CALL para_env%bcast(peak, para_env%source)
#else
      MARK_USED(order)
      MARK_USED(rows)
      MARK_USED(cols)
      MARK_USED(a)
      MARK_USED(b)
      MARK_USED(para_env)
      MARK_USED(tolerance)
      MARK_USED(metric_ratio)
      MARK_USED(gap)
#endif
      IF (PRESENT(peak_memory_bytes)) peak_memory_bytes = peak
   END SUBROUTINE sparse_z2_evaluate

! **************************************************************************************************
!> \brief Distributed AII transformation Q=(I-iC)/sqrt(2), followed by real skew LDL.
!> \param order Realified Hermitian pencil dimension
!> \param rows Distributed lower realified COO rows
!> \param cols Distributed lower realified COO columns
!> \param a Localizer entries
!> \param b Positive metric entries, verified independently by caller
!> \param para_env Communicator
!> \param tolerance Relative structure and solve tolerance
!> \param metric_ratio Lower metric eigenvalue divided by the original COO infinity-norm bound
!> \param gap Positive generalized gap lower bound
!> \param invariant Z2 relative to the fixed atomic reference orientation
!> \param residual Largest independent solve residual
!> \param fill Distributed factor scalar entries
!> \param status 0 success, 1 unresolved, -1 input, -3 AII structure, -6 memory limit
!> \param limit Per-rank numerical matrix/panel/work bytes, excluding COO and metadata
!> \param peak Peak accepted numerical buffer bytes
! **************************************************************************************************
   SUBROUTINE native_z2_evaluate(order, rows, cols, a, b, para_env, tolerance, metric_ratio, gap, &
                                 invariant, residual, fill, status, limit, peak)
      INTEGER, INTENT(IN) :: order, rows(:), cols(:)
      REAL(dp), INTENT(IN) :: a(:), b(:), tolerance, metric_ratio, gap
      TYPE(mp_para_env_type), POINTER :: para_env
      INTEGER, INTENT(OUT) :: invariant, status
      REAL(dp), INTENT(OUT) :: residual
      INTEGER(int_8), INTENT(OUT) :: fill, peak
      INTEGER(int_8), INTENT(IN) :: limit
      INTEGER :: dimensions(2), i, invalid, k, n, pf, quarter
      INTEGER(int_8) :: limits(2), current_peak
      INTEGER, ALLOCATABLE :: tr_rows(:), tr_cols(:), defect_rows(:), defect_cols(:)
      REAL(dp) :: a_scale, b_scale, error_a, error_b, metric_min, norm_a, norm_b, ranges(4, 2), symmetric_error
      REAL(dp), ALLOCATABLE :: real_a(:), imag_a(:), real_b(:), imag_b(:), bnorm(:), defect_values(:)
      COMPLEX(dp) :: av, bv

      invariant = 0
      residual = 0.0_dp
      fill = 0
      peak = 0
      status = -1
      invalid = 0
      IF (order < 8 .OR. MOD(order, 8) /= 0 .OR. limit <= 0) invalid = 1
      IF (SIZE(rows) /= SIZE(cols) .OR. SIZE(rows) /= SIZE(a) .OR. SIZE(rows) /= SIZE(b)) invalid = 1
      IF (SIZE(a) > HUGE(0)/16) invalid = 1
      IF (ANY(rows < 1) .OR. ANY(rows > order) .OR. ANY(cols < 1)) invalid = 1
      IF (SIZE(rows) == SIZE(cols)) THEN
         IF (ANY(cols > rows)) invalid = 1
      END IF
      IF (.NOT. ALL(ieee_is_finite(a)) .OR. .NOT. ALL(ieee_is_finite(b))) invalid = 1
      IF (.NOT. ALL(ieee_is_finite([tolerance, metric_ratio, gap]))) invalid = 1
      IF (tolerance <= 0.0_dp .OR. tolerance >= 1.0_dp .OR. &
          metric_ratio <= 0.0_dp .OR. metric_ratio > 1.0_dp .OR. gap <= 0.0_dp) invalid = 1
      CALL para_env%max(invalid)
      IF (invalid /= 0) RETURN
      dimensions = order
      limits = limit
      ranges(:, 1) = [tolerance, metric_ratio, gap, REAL(order, dp)]
      ranges(:, 2) = ranges(:, 1)
      CALL para_env%min(dimensions(1))
      CALL para_env%max(dimensions(2))
      CALL para_env%min(limits(1))
      CALL para_env%max(limits(2))
      CALL para_env%min(ranges(:, 1))
      CALL para_env%max(ranges(:, 2))
      IF (dimensions(1) /= dimensions(2) .OR. limits(1) /= limits(2) .OR. ANY(ranges(:, 1) /= ranges(:, 2))) RETURN
      status = -6
      IF (INT(order, int_8) > limit/8_int_8) RETURN
      a_scale = 0.0_dp
      b_scale = 0.0_dp
      IF (SIZE(a) > 0) THEN
         a_scale = MAXVAL(ABS(a))
         b_scale = MAXVAL(ABS(b))
      END IF
      CALL para_env%max(a_scale)
      CALL para_env%max(b_scale)
      status = 1
      IF (a_scale <= TINY(1.0_dp) .OR. b_scale <= TINY(1.0_dp)) RETURN
      n = order/2
      quarter = n/4
      ALLOCATE (bnorm(order))
      bnorm = 0.0_dp
      DO i = 1, SIZE(b)
         bnorm(rows(i)) = bnorm(rows(i)) + ABS(b(i)/b_scale)
         IF (rows(i) /= cols(i)) bnorm(cols(i)) = bnorm(cols(i)) + ABS(b(i)/b_scale)
      END DO
      CALL para_env%sum(bnorm)
      metric_min = metric_ratio*MAXVAL(bnorm)
      DEALLOCATE (bnorm)
      peak = 8_int_8*order
      ALLOCATE (tr_rows(8*SIZE(a)), tr_cols(8*SIZE(a)), real_a(8*SIZE(a)), imag_a(8*SIZE(a)), &
                real_b(8*SIZE(a)), imag_b(8*SIZE(a)))
      k = 0
      DO i = 1, SIZE(a)
         IF (rows(i) <= n) THEN
            av = CMPLX(a(i)/a_scale, 0.0_dp, dp)
            bv = CMPLX(b(i)/b_scale, 0.0_dp, dp)
            CALL transform_entry(rows(i), cols(i), av, bv)
            IF (rows(i) /= cols(i)) CALL transform_entry(cols(i), rows(i), av, bv)
         ELSE IF (cols(i) <= n) THEN
            av = CMPLX(0.0_dp, a(i)/a_scale, dp)
            bv = CMPLX(0.0_dp, b(i)/b_scale, dp)
            CALL transform_entry(rows(i) - n, cols(i), av, bv)
         END IF
      END DO
      CALL localizer_coo_norm(n, tr_rows(:k), tr_cols(:k), real_a(:k), para_env, norm_a, current_peak, status, limit)
      peak = MAX(peak, current_peak)
      IF (status /= 0) RETURN
      CALL localizer_coo_norm(n, tr_rows(:k), tr_cols(:k), imag_a(:k), para_env, error_a, current_peak, status, limit)
      peak = MAX(peak, current_peak)
      IF (status /= 0) RETURN
      CALL localizer_coo_norm(n, tr_rows(:k), tr_cols(:k), real_b(:k), para_env, norm_b, current_peak, status, limit)
      peak = MAX(peak, current_peak)
      IF (status /= 0) RETURN
      CALL localizer_coo_norm(n, tr_rows(:k), tr_cols(:k), imag_b(:k), para_env, error_b, current_peak, status, limit)
      peak = MAX(peak, current_peak)
      IF (status /= 0) RETURN
      ALLOCATE (defect_rows(2*k), defect_cols(2*k), defect_values(2*k))
      defect_rows(:k) = tr_rows(:k)
      defect_rows(k + 1:) = tr_cols(:k)
      defect_cols(:k) = tr_cols(:k)
      defect_cols(k + 1:) = tr_rows(:k)
      defect_values(:k) = 0.5_dp*real_a(:k)
      defect_values(k + 1:) = defect_values(:k)
      CALL localizer_coo_norm(n, defect_rows, defect_cols, defect_values, para_env, symmetric_error, &
                              current_peak, status, limit)
      peak = MAX(peak, current_peak)
      IF (status /= 0) RETURN
      error_a = error_a + symmetric_error
      status = -3
      IF (error_a > tolerance*norm_a .OR. error_b > tolerance*norm_b) RETURN
      status = 1
      IF (error_a > 0.0_dp) THEN
         ! Compare physical units logarithmically without overflow from disparate L/B scales.
         IF (LOG(error_a) + LOG(a_scale) >= LOG(0.25_dp) + LOG(gap) + LOG(metric_min) + LOG(b_scale)) RETURN
      END IF
      CALL localizer_dbcsr_pfaffian(n, tr_rows(:k), tr_cols(:k), real_a(:k), para_env, tolerance, pf, residual, &
                                    fill, status, limit, current_peak)
      peak = MAX(peak, current_peak)
      IF (status == 0) invariant = (1 - (-1)**quarter*pf)/2

   CONTAINS

! **************************************************************************************************
!> \brief Emit four coordinate contributions to i Q^H L Q and Q^H B Q.
!> \param row Complex input row
!> \param col Complex input column
!> \param a_value Normalized complex localizer entry
!> \param b_value Normalized complex metric entry
! **************************************************************************************************
      SUBROUTINE transform_entry(row, col, a_value, b_value)
      INTEGER, INTENT(IN)                                :: row, col
      COMPLEX(dp), INTENT(IN)                            :: a_value, b_value

      COMPLEX(dp)                                        :: factors(4), va, vb
      INTEGER                                            :: index, partner_col, partner_row
      REAL(dp)                                           :: sign_col, sign_row

         IF (ABS(a_value) == 0.0_dp .AND. ABS(b_value) == 0.0_dp) RETURN
         partner_row = (3 - (row - 1)/quarter)*quarter + MOD(row - 1, quarter) + 1
         partner_col = (3 - (col - 1)/quarter)*quarter + MOD(col - 1, quarter) + 1
         sign_row = -1.0_dp
         sign_col = -1.0_dp
         IF ((row - 1)/quarter == 0 .OR. (row - 1)/quarter == 3) sign_row = 1.0_dp
         IF ((col - 1)/quarter == 0 .OR. (col - 1)/quarter == 3) sign_col = 1.0_dp
         tr_rows(k + 1:k + 4) = [row, partner_row, partner_row, row]
         tr_cols(k + 1:k + 4) = [col, partner_col, col, partner_col]
         factors = [CMPLX(1.0_dp, 0.0_dp, dp), CMPLX(sign_row*sign_col, 0.0_dp, dp), &
                    CMPLX(0.0_dp, sign_row, dp), CMPLX(0.0_dp, -sign_col, dp)]
         DO index = 1, 4
            va = CMPLX(0.0_dp, 0.5_dp, dp)*a_value*factors(index)
            vb = 0.5_dp*b_value*factors(index)
            real_a(k + index) = REAL(va, dp)
            imag_a(k + index) = AIMAG(va)
            real_b(k + index) = REAL(vb, dp)
            imag_b(k + index) = AIMAG(vb)
         END DO
         k = k + 4
      END SUBROUTINE transform_entry
   END SUBROUTINE native_z2_evaluate

END MODULE localizer_sparse_pfaffian
