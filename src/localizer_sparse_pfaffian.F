!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Native sparse Pfaffian adapter; only coordinate data are gathered to the source rank.
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
   USE message_passing, ONLY: mp_para_env_type

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
! **************************************************************************************************
   SUBROUTINE sparse_pfaffian_sign(order, rows, cols, values, tolerance, sign, residual, fill, status, &
                                   max_memory_bytes, peak_memory_bytes)
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
! **************************************************************************************************
   SUBROUTINE sparse_z2_evaluate(order, rows, cols, a, b, para_env, tolerance, metric_ratio, gap, &
                                 invariant, residual, fill, status, max_memory_bytes, peak_memory_bytes)
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
      INTEGER(int_8) :: limit, peak
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
#endif
      IF (PRESENT(peak_memory_bytes)) peak_memory_bytes = peak
   END SUBROUTINE sparse_z2_evaluate

END MODULE localizer_sparse_pfaffian
