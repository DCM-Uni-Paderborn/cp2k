!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Distributed pivoted LDL factors without sparse-direct-solver dependencies.
!> \author Thomas D. Kuehne, tkuehne@cp2k.org
! **************************************************************************************************
MODULE localizer_dbcsr_factor
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_create, dbcsr_distribution_new, dbcsr_distribution_release, dbcsr_distribution_type, &
        dbcsr_finalize, dbcsr_get_block_p, dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, &
        dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_put_block, &
        dbcsr_release, dbcsr_type, dbcsr_type_no_symmetry
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE message_passing,                 ONLY: mp_para_env_type
   USE util,                            ONLY: sort

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: localizer_factor_type, localizer_factorize, localizer_factor_solve, localizer_factor_release
   PUBLIC :: localizer_skew_coo, localizer_dbcsr_pfaffian, localizer_coo_norm

   INTEGER, PARAMETER :: tile_size = 16

   TYPE factor_panel_type
      INTEGER :: pivot(2) = 0, width = 0
      REAL(dp) :: inverse(2, 2) = 0.0_dp
      INTEGER, ALLOCATABLE :: rows(:)
      REAL(dp), ALLOCATABLE :: values(:, :)
   END TYPE factor_panel_type

   TYPE localizer_factor_type
      TYPE(factor_panel_type), ALLOCATABLE :: panels(:)
      INTEGER :: order = 0, steps = 0
      ! Real-buffer estimate excludes COO routing buffers and DBCSR metadata.
      REAL(dp) :: scale = 0.0_dp, peak_memory_mb = 0.0_dp
      INTEGER(int_8) :: entries = 0, workspace_entries = 0
      INTEGER(int_8) :: peak_memory_bytes = 0, retained_entries = 0, memory_limit = HUGE(0_int_8)
      LOGICAL :: valid = .FALSE.
   END TYPE localizer_factor_type

CONTAINS

! **************************************************************************************************
!> \brief Complete-pivoting block LDL on a distributed full-storage DBCSR Schur complement.
!> \param order Matrix dimension
!> \param rows Distributed lower-triangle coordinates, with additive duplicates
!> \param cols Distributed lower-triangle coordinates
!> \param values Matrix entries
!> \param para_env Matrix communicator
!> \param factor Factors retained only when requested
!> \param negative Negative inertia count
!> \param status 0 success, 1 numerically singular, -1 invalid input, -4 numerical failure
!> \param retain Retain distributed elimination panels for repeated solves
!> \param skew Factor a real skew-symmetric matrix with oriented 2x2 pivots
!> \param tolerance Relative numerical pivot threshold
!> \param pf_sign Pfaffian sign for skew input, zero if unresolved
!> \param max_memory_bytes Per-rank numerical matrix/panel/work buffer limit, excluding COO and metadata
!> \note No fill is dropped and no pivot is perturbed. Only two matrix columns are replicated.
! **************************************************************************************************
   SUBROUTINE localizer_factorize(order, rows, cols, values, para_env, factor, negative, status, retain, &
                                  skew, tolerance, pf_sign, max_memory_bytes)
      INTEGER, INTENT(IN)                                :: order, rows(:), cols(:)
      REAL(dp), INTENT(IN)                               :: values(:)
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(localizer_factor_type), INTENT(INOUT)         :: factor
      INTEGER, INTENT(OUT)                               :: negative, status
      LOGICAL, INTENT(IN)                                :: retain
      LOGICAL, INTENT(IN), OPTIONAL                      :: skew
      REAL(dp), INTENT(IN), OPTIONAL                     :: tolerance
      INTEGER, INTENT(OUT), OPTIONAL                     :: pf_sign
      INTEGER(int_8), INTENT(IN), OPTIONAL               :: max_memory_bytes

      INTEGER                                            :: col, i, invalid, j, k, maximum(3), &
                                                            minimum(3), nrows, p(2), parity, &
                                                            remaining, row, sign_accum, step, width
      INTEGER(int_8)                                     :: candidate(2), key, limit_range(2), &
                                                            new_entries, panel_entries, reserved, &
                                                            stored
      LOGICAL                                            :: skew_mode
      LOGICAL, ALLOCATABLE                               :: active(:)
      REAL(dp)                                           :: alpha, bound, d(2, 2), determinant, &
                                                            eigenvalues(2), inverse(2, 2), &
                                                            maxima(2), scale, threshold_range(2)
      REAL(dp), ALLOCATABLE                              :: columns(:, :), multipliers(:, :)
      REAL(dp), POINTER                                  :: block(:, :)
      TYPE(dbcsr_iterator_type)                          :: iterator
      TYPE(dbcsr_type)                                   :: matrix

      CALL localizer_factor_release(factor)
      IF (PRESENT(max_memory_bytes)) factor%memory_limit = max_memory_bytes
      negative = 0
      IF (PRESENT(pf_sign)) pf_sign = 0
      status = -1
      skew_mode = .FALSE.
      IF (PRESENT(skew)) skew_mode = skew
      bound = 32.0_dp*EPSILON(1.0_dp)*REAL(order, dp)
      IF (PRESENT(tolerance)) bound = MAX(bound, tolerance)
      invalid = 0
      IF (order < 1 .OR. SIZE(rows) /= SIZE(cols) .OR. SIZE(rows) /= SIZE(values)) invalid = 1
      IF (skew_mode .AND. MOD(order, 2) /= 0) invalid = 1
      IF (PRESENT(tolerance)) THEN
         IF (.NOT. ieee_is_finite(tolerance)) invalid = 1
         IF (tolerance <= 0.0_dp .OR. tolerance >= 1.0_dp) invalid = 1
      END IF
      IF (ANY(rows < 1) .OR. ANY(rows > order) .OR. ANY(cols < 1)) invalid = 1
      IF (SIZE(rows) == SIZE(cols)) THEN
         IF (ANY(cols > rows)) invalid = 1
      END IF
      IF (.NOT. ALL(ieee_is_finite(values))) invalid = 1
      IF (factor%memory_limit <= 0) invalid = 1
      CALL para_env%max(invalid)
      IF (invalid /= 0) RETURN
      minimum = [order, MERGE(1, 0, retain), MERGE(1, 0, skew_mode)]
      maximum = minimum
      CALL para_env%min(minimum)
      CALL para_env%max(maximum)
      IF (ANY(minimum /= maximum)) RETURN
      threshold_range = bound
      CALL para_env%min(threshold_range(1))
      CALL para_env%max(threshold_range(2))
      IF (threshold_range(1) /= threshold_range(2)) RETURN
      limit_range = factor%memory_limit
      CALL para_env%min(limit_range(1))
      CALL para_env%max(limit_range(2))
      IF (limit_range(1) /= limit_range(2)) RETURN
      scale = 0.0_dp
      IF (SIZE(values) > 0) scale = MAXVAL(ABS(values))
      CALL para_env%max(scale)
      status = 1
      IF (scale <= TINY(1.0_dp)) RETURN
      factor%order = order
      factor%scale = scale
      reserved = 4_int_8*order + tile_size**2
      IF (retain) reserved = reserved + 4_int_8*order
      CALL import_matrix(order, rows, cols, values/scale, para_env, matrix, skew_mode, &
                         factor%memory_limit, reserved, stored, status)
      IF (status /= 0) THEN
         CALL dbcsr_release(matrix)
         RETURN
      END IF
      CALL reserve_buffers(stored + reserved, factor, para_env, status)
      IF (status /= 0) THEN
         CALL dbcsr_release(matrix)
         RETURN
      END IF
      ALLOCATE (active(order), columns(order, 2), multipliers(order, 2))
      IF (retain) ALLOCATE (factor%panels(order))
      active = .TRUE.
      remaining = order
      panel_entries = 0
      alpha = (1.0_dp + SQRT(17.0_dp))/8.0_dp
      step = 0
      sign_accum = 1
      DO WHILE (remaining > 0)
         maxima = 0.0_dp
         invalid = 0
         stored = 0
         CALL dbcsr_iterator_start(iterator, matrix)
         DO WHILE (dbcsr_iterator_blocks_left(iterator))
            CALL dbcsr_iterator_next_block(iterator, block=block, row_offset=row, col_offset=col)
            stored = stored + SIZE(block, KIND=int_8)
            IF (.NOT. ALL(ieee_is_finite(block))) invalid = 1
            DO j = 1, SIZE(block, 2)
               IF (.NOT. active(col + j - 1)) CYCLE
               DO i = 1, SIZE(block, 1)
                  IF (.NOT. active(row + i - 1)) CYCLE
                  IF (row + i < col + j) CYCLE
                  k = 2
                  IF (row + i == col + j) k = 1
                  maxima(k) = MAX(maxima(k), ABS(block(i, j)))
               END DO
            END DO
         END DO
         CALL dbcsr_iterator_stop(iterator)
         CALL reserve_buffers(stored + panel_entries + reserved, factor, para_env, status)
         CALL para_env%sum(stored)
         factor%workspace_entries = MAX(factor%workspace_entries, stored)
         CALL para_env%max(invalid)
         CALL para_env%max(maxima)
         status = -4
         IF (invalid /= 0) EXIT
         status = 1
         IF (MAXVAL(maxima) <= bound) EXIT
         IF (skew_mode .AND. maxima(1) > bound) THEN
            status = -1
            EXIT
         END IF
         ! Bunch-Parlett selection bounds multipliers and forces indefinite 2x2 pivots.
         width = 1
         IF (maxima(1) < alpha*maxima(2)) width = 2
         IF (skew_mode) width = 2
         candidate = HUGE(0_int_8)
         CALL dbcsr_iterator_start(iterator, matrix)
         DO WHILE (dbcsr_iterator_blocks_left(iterator))
            CALL dbcsr_iterator_next_block(iterator, block=block, row_offset=row, col_offset=col)
            DO j = 1, SIZE(block, 2)
               IF (.NOT. active(col + j - 1)) CYCLE
               DO i = 1, SIZE(block, 1)
                  IF (.NOT. active(row + i - 1)) CYCLE
                  IF (row + i == col + j) THEN
                     IF (ABS(block(i, j)) == maxima(1)) candidate(1) = MIN(candidate(1), INT(row + i - 1, int_8))
                  ELSE IF (row + i > col + j) THEN
                     key = INT(row + i - 1, int_8) + INT(order, int_8)*INT(col + j - 2, int_8)
                     IF (ABS(block(i, j)) == maxima(2)) candidate(2) = MIN(candidate(2), key)
                  END IF
               END DO
            END DO
         END DO
         CALL dbcsr_iterator_stop(iterator)
         CALL para_env%min(candidate)
         p = 0
         IF (width == 1) THEN
            p(1) = INT(candidate(1))
         ELSE
            p(1) = INT(MOD(candidate(2) - 1, INT(order, int_8))) + 1
            p(2) = INT((candidate(2) - 1)/INT(order, int_8)) + 1
         END IF
         CALL get_columns(matrix, p(:width), para_env, columns(:, :width))
         d = 0.0_dp
         inverse = 0.0_dp
         DO i = 1, width
            d(i, :width) = columns(p(i), :width)
         END DO
         IF (skew_mode) THEN
            IF (ABS(d(1, 2)) <= bound) EXIT
            IF (MAXVAL(ABS(d + TRANSPOSE(d))) > bound) THEN
               status = -4
               EXIT
            END IF
            inverse(1, 2) = -1.0_dp/d(1, 2)
            inverse(2, 1) = -inverse(1, 2)
            ! Removing pivots in this order permutes the active original indices.
            parity = COUNT(active(:p(1) - 1)) + COUNT(active(:p(2) - 1))
            IF (p(1) < p(2)) parity = parity - 1
            IF (MOD(parity, 2) /= 0) sign_accum = -sign_accum
            IF (d(1, 2) < 0.0_dp) sign_accum = -sign_accum
         ELSE IF (width == 1) THEN
            eigenvalues(1) = d(1, 1)
            IF (ABS(eigenvalues(1)) <= bound) EXIT
            inverse(1, 1) = 1.0_dp/d(1, 1)
         ELSE
            CALL dlae2(d(1, 1), d(1, 2), d(2, 2), eigenvalues(1), eigenvalues(2))
            IF (MINVAL(ABS(eigenvalues)) <= bound) EXIT
            scale = MAXVAL(ABS(d))
            determinant = (d(1, 1)/scale)*(d(2, 2)/scale) - (d(1, 2)/scale)**2
            inverse = RESHAPE([d(2, 2)/scale, -d(1, 2)/scale, -d(1, 2)/scale, d(1, 1)/scale], [2, 2])
            inverse = inverse/(determinant*scale)
         END IF
         status = -4
         IF (.NOT. ALL(ieee_is_finite(inverse))) EXIT
         IF (.NOT. skew_mode) negative = negative + COUNT(eigenvalues(:width) < 0.0_dp)
         active(p(:width)) = .FALSE.
         DO i = 1, order
            IF (.NOT. active(i)) columns(i, :width) = 0.0_dp
         END DO
         multipliers(:, :width) = MATMUL(columns(:, :width), inverse(:width, :width))
         step = step + 1
         nrows = 0
         DO i = 1, order
            IF (MOD((i - 1)/tile_size, para_env%num_pe) /= para_env%mepos) CYCLE
            IF (ANY(multipliers(i, :width) /= 0.0_dp)) nrows = nrows + 1
         END DO
         factor%entries = factor%entries + INT(nrows, int_8)*width
         IF (MOD((p(1) - 1)/tile_size, para_env%num_pe) == para_env%mepos) THEN
            factor%entries = factor%entries + width*(width + 1)/2
         END IF
         IF (retain) THEN
            CALL update_storage(matrix, order, columns(:, :width), multipliers(:, :width), para_env, new_entries)
            ! Reserve the complete fill update and panel before either allocation.
            CALL reserve_buffers(new_entries + panel_entries + INT(nrows, int_8)*width + reserved, &
                                 factor, para_env, status)
            IF (status /= 0) EXIT
            factor%panels(step)%width = width
            factor%panels(step)%pivot = p
            factor%panels(step)%inverse = inverse
            ALLOCATE (factor%panels(step)%rows(nrows), factor%panels(step)%values(nrows, width))
            nrows = 0
            DO i = 1, order
               IF (MOD((i - 1)/tile_size, para_env%num_pe) /= para_env%mepos) CYCLE
               IF (.NOT. ANY(multipliers(i, :width) /= 0.0_dp)) CYCLE
               nrows = nrows + 1
               factor%panels(step)%rows(nrows) = i
               factor%panels(step)%values(nrows, :) = multipliers(i, :width)
            END DO
            panel_entries = panel_entries + INT(nrows, int_8)*width
         ELSE
            CALL update_storage(matrix, order, columns(:, :width), multipliers(:, :width), para_env, new_entries)
            CALL reserve_buffers(new_entries + reserved, factor, para_env, status)
            IF (status /= 0) EXIT
         END IF
         CALL update_matrix(matrix, order, columns(:, :width), multipliers(:, :width), para_env, skew_mode)
         remaining = remaining - width
         status = 0
      END DO
      factor%steps = step
      IF (retain) factor%retained_entries = 4_int_8*order + panel_entries
      CALL para_env%sum(factor%entries)
      CALL para_env%max(factor%peak_memory_mb)
      CALL para_env%max(factor%peak_memory_bytes)
      CALL dbcsr_release(matrix)
      IF (remaining > 0) THEN
         factor%valid = .FALSE.
         RETURN
      END IF
      factor%valid = retain
      IF (skew_mode .AND. PRESENT(pf_sign)) pf_sign = sign_accum
      status = 0
   END SUBROUTINE localizer_factorize

! **************************************************************************************************
!> \brief Route additive COO entries directly to DBCSR row-tile owners, without gathering.
!> \param order Matrix order
!> \param rows Lower-triangle rows
!> \param cols Lower-triangle columns
!> \param values Normalized entries
!> \param para_env Communicator
!> \param matrix New distributed full-storage matrix
!> \param skew Mirror lower entries with a minus sign
!> \param limit Numerical matrix/work byte limit per rank
!> \param reserved Other numerical buffers reserved by caller
!> \param stored Number of local matrix scalar slots
!> \param status 0 success, -6 requested numerical storage exceeds limit
!> \param mirror_entries Mirror lower COO when true, otherwise import arbitrary full COO
! **************************************************************************************************
   SUBROUTINE import_matrix(order, rows, cols, values, para_env, matrix, skew, limit, reserved, stored, status, &
                            mirror_entries)
      INTEGER, INTENT(IN)                                :: order, rows(:), cols(:)
      REAL(dp), INTENT(IN)                               :: values(:)
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_type), INTENT(OUT)                      :: matrix
      LOGICAL, INTENT(IN)                                :: skew
      INTEGER(int_8), INTENT(IN)                         :: limit, reserved
      INTEGER(int_8), INTENT(OUT)                        :: stored
      INTEGER, INTENT(OUT)                               :: status
      LOGICAL, INTENT(IN), OPTIONAL                      :: mirror_entries

      INTEGER                                            :: col, dest, i, invalid, ir, jc, mirror, &
                                                            nb, nmirror, np, pass, row, slot
      INTEGER(int_8), ALLOCATABLE                        :: keys(:)
      INTEGER, ALLOCATABLE :: key_order(:), position(:), recv_count(:), recv_disp(:), &
         recv_index(:, :), send_count(:), send_disp(:), send_index(:, :)
      INTEGER, POINTER                                   :: col_dist(:), grid(:, :), row_dist(:), &
                                                            sizes(:)
      REAL(dp)                                           :: buffer(tile_size, tile_size)
      REAL(dp), ALLOCATABLE                              :: recv_value(:), send_value(:)
      TYPE(dbcsr_distribution_type)                      :: distribution

      np = para_env%num_pe
      nmirror = 2
      IF (PRESENT(mirror_entries)) THEN
         IF (.NOT. mirror_entries) nmirror = 1
      END IF
      nb = (order - 1)/tile_size + 1
      ALLOCATE (grid(0:np - 1, 0:0), row_dist(nb), col_dist(nb), sizes(nb))
      grid(:, 0) = [(i, i=0, np - 1)]
      row_dist = [(MOD(i - 1, np), i=1, nb)]
      col_dist = 0
      sizes = tile_size
      sizes(nb) = order - (nb - 1)*tile_size
      CALL dbcsr_distribution_new(distribution, group=para_env%get_handle(), pgrid=grid, &
                                  row_dist=row_dist, col_dist=col_dist)
      CALL dbcsr_create(matrix, 'LOCALIZER SCHUR', distribution, dbcsr_type_no_symmetry, sizes, sizes, &
                        mutable_work=.TRUE.)
      CALL dbcsr_distribution_release(distribution)
      ALLOCATE (send_count(0:np - 1), recv_count(0:np - 1), send_disp(0:np - 1), &
                recv_disp(0:np - 1), position(0:np - 1))
      send_count = 0
      send_disp = 0
      DO pass = 1, 2
         position = send_disp
         DO i = 1, SIZE(values)
            IF (values(i) == 0.0_dp) CYCLE
            DO mirror = 1, nmirror
               IF (mirror == 2 .AND. rows(i) == cols(i)) CYCLE
               row = rows(i)
               col = cols(i)
               IF (mirror == 2) THEN
                  row = cols(i)
                  col = rows(i)
               END IF
               dest = MOD((row - 1)/tile_size, np)
               IF (pass == 1) THEN
                  send_count(dest) = send_count(dest) + 1
               ELSE
                  position(dest) = position(dest) + 1
                  slot = position(dest)
                  send_index(:, slot) = [row, col]
                  send_value(slot) = values(i)
                  IF (skew .AND. mirror == 2) send_value(slot) = -values(i)
               END IF
            END DO
         END DO
         IF (pass == 1) THEN
            CALL para_env%alltoall(send_count, recv_count, 1)
            send_disp(0) = 0
            recv_disp(0) = 0
            DO i = 1, np - 1
               send_disp(i) = send_disp(i - 1) + send_count(i - 1)
               recv_disp(i) = recv_disp(i - 1) + recv_count(i - 1)
            END DO
            ALLOCATE (send_index(2, MAX(1, SUM(send_count))), recv_index(2, MAX(1, SUM(recv_count))), &
                      send_value(MAX(1, SUM(send_count))), recv_value(MAX(1, SUM(recv_count))))
         END IF
      END DO
      CALL para_env%alltoall(send_index, 2*send_count, 2*send_disp, recv_index, 2*recv_count, 2*recv_disp)
      CALL para_env%alltoall(send_value, send_count, send_disp, recv_value, recv_count, recv_disp)
      ALLOCATE (keys(SUM(recv_count)), key_order(SUM(recv_count)))
      DO i = 1, SIZE(keys)
         row = (recv_index(1, i) - 1)/tile_size + 1
         col = (recv_index(2, i) - 1)/tile_size + 1
         keys(i) = INT(row - 1, int_8)*nb + col
      END DO
      IF (SIZE(keys) > 0) CALL sort(keys, SIZE(keys), key_order)
      stored = 0
      DO i = 1, SIZE(keys)
         IF (i > 1) THEN
            IF (keys(i) == keys(i - 1)) CYCLE
         END IF
         row = INT((keys(i) - 1)/nb) + 1
         col = INT(MOD(keys(i) - 1, INT(nb, int_8))) + 1
         stored = stored + INT(sizes(row), int_8)*sizes(col)
      END DO
      invalid = MERGE(1, 0, stored + reserved > limit/8_int_8)
      CALL para_env%max(invalid)
      status = -6
      IF (invalid /= 0) THEN
         DEALLOCATE (grid, row_dist, col_dist, sizes)
         RETURN
      END IF
      DO i = 1, SUM(recv_count)
         row = (recv_index(1, i) - 1)/tile_size + 1
         col = (recv_index(2, i) - 1)/tile_size + 1
         ir = MOD(recv_index(1, i) - 1, tile_size) + 1
         jc = MOD(recv_index(2, i) - 1, tile_size) + 1
         buffer = 0.0_dp
         buffer(ir, jc) = recv_value(i)
         CALL dbcsr_put_block(matrix, row, col, buffer(:sizes(row), :sizes(col)), summation=.TRUE.)
      END DO
      CALL dbcsr_finalize(matrix)
      DEALLOCATE (grid, row_dist, col_dist, sizes)
      status = 0
   END SUBROUTINE import_matrix

! **************************************************************************************************
!> \brief Bound the spectral norm after distributed addition of arbitrary full COO entries.
!> \param order Scalar dimension, validated by caller
!> \param rows Distributed full COO rows, validated by caller
!> \param cols Distributed full COO columns
!> \param values Finite additive entries
!> \param para_env Communicator
!> \param norm Maximum of the assembled one- and infinity-norms
!> \param peak Peak accepted numerical matrix/work bytes
!> \param status 0 success, -1 nonfinite aggregate, -6 memory limit
!> \param limit Per-rank numerical buffer byte limit
! **************************************************************************************************
   SUBROUTINE localizer_coo_norm(order, rows, cols, values, para_env, norm, peak, status, limit)
      INTEGER, INTENT(IN)                                :: order, rows(:), cols(:)
      REAL(dp), INTENT(IN)                               :: values(:)
      TYPE(mp_para_env_type), POINTER                    :: para_env
      REAL(dp), INTENT(OUT)                              :: norm
      INTEGER(int_8), INTENT(OUT)                        :: peak
      INTEGER, INTENT(OUT)                               :: status
      INTEGER(int_8), INTENT(IN)                         :: limit

      INTEGER                                            :: col, invalid, row
      INTEGER(int_8)                                     :: reserved, stored
      REAL(dp), ALLOCATABLE                              :: col_norm(:), row_norm(:)
      REAL(dp), POINTER                                  :: block(:, :)
      TYPE(dbcsr_iterator_type)                          :: iterator
      TYPE(dbcsr_type)                                   :: matrix

      norm = 0.0_dp
      peak = 0
      reserved = 2_int_8*order + tile_size**2
      status = -6
      IF (reserved > limit/8_int_8) RETURN
      ALLOCATE (row_norm(order), col_norm(order))
      row_norm = 0.0_dp
      col_norm = 0.0_dp
      peak = 8_int_8*reserved
      CALL import_matrix(order, rows, cols, values, para_env, matrix, .FALSE., limit, reserved, stored, status, &
                         mirror_entries=.FALSE.)
      IF (status == 0) THEN
         peak = 8_int_8*(stored + reserved)
         invalid = 0
         CALL dbcsr_iterator_start(iterator, matrix)
         DO WHILE (dbcsr_iterator_blocks_left(iterator))
            CALL dbcsr_iterator_next_block(iterator, block=block, row_offset=row, col_offset=col)
            IF (.NOT. ALL(ieee_is_finite(block))) invalid = 1
            row_norm(row:row + SIZE(block, 1) - 1) = row_norm(row:row + SIZE(block, 1) - 1) + SUM(ABS(block), DIM=2)
            col_norm(col:col + SIZE(block, 2) - 1) = col_norm(col:col + SIZE(block, 2) - 1) + SUM(ABS(block), DIM=1)
         END DO
         CALL dbcsr_iterator_stop(iterator)
         CALL para_env%sum(row_norm)
         CALL para_env%sum(col_norm)
         CALL para_env%max(invalid)
         norm = MAX(MAXVAL(row_norm), MAXVAL(col_norm))
         IF (invalid /= 0 .OR. .NOT. ieee_is_finite(norm)) status = -1
      END IF
      CALL dbcsr_release(matrix)
      CALL para_env%max(peak)
   END SUBROUTINE localizer_coo_norm

! **************************************************************************************************
!> \brief Check collective numerical buffer budget before allocating matrix fill or factor panels.
!> \param entries Local numerical scalar slots
!> \param factor Factor context and accepted peak bytes
!> \param para_env Communicator
!> \param status 0 accepted, -6 rejected before allocation
! **************************************************************************************************
   SUBROUTINE reserve_buffers(entries, factor, para_env, status)
      INTEGER(int_8), INTENT(IN)                         :: entries
      TYPE(localizer_factor_type), INTENT(INOUT)         :: factor
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(OUT)                               :: status

      INTEGER                                            :: invalid

      invalid = MERGE(1, 0, entries > factor%memory_limit/8_int_8)
      CALL para_env%max(invalid)
      status = -6
      IF (invalid /= 0) RETURN
      factor%peak_memory_bytes = MAX(factor%peak_memory_bytes, 8_int_8*entries)
      factor%peak_memory_mb = REAL(factor%peak_memory_bytes, dp)/1048576.0_dp
      status = 0
   END SUBROUTINE reserve_buffers

! **************************************************************************************************
!> \brief Count local matrix slots after the impending update, without allocating new blocks.
!> \param matrix Schur complement
!> \param order Scalar dimension
!> \param columns Pivot columns
!> \param multipliers Update multipliers
!> \param para_env Row-tile communicator
!> \param entries Existing storage plus all required fill tiles
! **************************************************************************************************
   SUBROUTINE update_storage(matrix, order, columns, multipliers, para_env, entries)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      INTEGER, INTENT(IN)                                :: order
      REAL(dp), INTENT(IN)                               :: columns(:, :), multipliers(:, :)
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER(int_8), INTENT(OUT)                        :: entries

      INTEGER                                            :: col, finish, first, last, nb, row, start
      LOGICAL                                            :: found
      REAL(dp), POINTER                                  :: block(:, :)
      TYPE(dbcsr_iterator_type)                          :: iterator

      entries = 0
      CALL dbcsr_iterator_start(iterator, matrix)
      DO WHILE (dbcsr_iterator_blocks_left(iterator))
         CALL dbcsr_iterator_next_block(iterator, block=block)
         entries = entries + SIZE(block, KIND=int_8)
      END DO
      CALL dbcsr_iterator_stop(iterator)
      nb = (order - 1)/tile_size + 1
      DO row = 1, nb
         IF (MOD(row - 1, para_env%num_pe) /= para_env%mepos) CYCLE
         first = (row - 1)*tile_size + 1
         last = MIN(row*tile_size, order)
         IF (.NOT. ANY(multipliers(first:last, :) /= 0.0_dp)) CYCLE
         DO col = 1, nb
            start = (col - 1)*tile_size + 1
            finish = MIN(col*tile_size, order)
            IF (.NOT. ANY(columns(start:finish, :) /= 0.0_dp)) CYCLE
            CALL dbcsr_get_block_p(matrix, row, col, block, found=found)
            IF (.NOT. found) entries = entries + INT(last - first + 1, int_8)*(finish - start + 1)
         END DO
      END DO
   END SUBROUTINE update_storage

! **************************************************************************************************
!> \brief Validate additive full COO collectively and extract its projected skew lower triangle.
!> \param order Even scalar dimension
!> \param rows Distributed full COO rows
!> \param cols Distributed full COO columns
!> \param values Additive entries, neither triangle is assumed complete on any rank
!> \param para_env Communicator
!> \param tolerance Relative infinity-norm skew defect tolerance
!> \param lower_rows Owned lower-triangle coordinates
!> \param lower_cols Owned lower-triangle coordinates
!> \param lower_values Projected values normalized by a positive input scale
!> \param defect Infinity norm of discarded symmetric part, in the same normalization
!> \param peak Peak accepted numerical matrix/work bytes per rank
!> \param status 0 success, 1 zero matrix, -1 input or structure, -6 numerical memory limit
!> \param limit Numerical buffer byte limit, excluding COO routing and metadata
! **************************************************************************************************
   SUBROUTINE localizer_skew_coo(order, rows, cols, values, para_env, tolerance, lower_rows, lower_cols, &
                                 lower_values, defect, peak, status, limit)
      INTEGER, INTENT(IN)                                :: order, rows(:), cols(:)
      REAL(dp), INTENT(IN)                               :: values(:)
      TYPE(mp_para_env_type), POINTER                    :: para_env
      REAL(dp), INTENT(IN)                               :: tolerance
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: lower_rows(:), lower_cols(:)
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: lower_values(:)
      REAL(dp), INTENT(OUT)                              :: defect
      INTEGER(int_8), INTENT(OUT)                        :: peak
      INTEGER, INTENT(OUT)                               :: status
      INTEGER(int_8), INTENT(IN)                         :: limit

      INTEGER                                            :: col, dimensions(2), i, invalid, j, k, &
                                                            pass, row
      INTEGER(int_8)                                     :: limits(2), reserved, stored
      INTEGER, ALLOCATABLE                               :: canonical_cols(:), canonical_rows(:)
      REAL(dp)                                           :: norm, scale, tolerances(2)
      REAL(dp), ALLOCATABLE                              :: canonical(:), row_norm(:)
      REAL(dp), POINTER                                  :: block(:, :)
      TYPE(dbcsr_iterator_type)                          :: iterator
      TYPE(dbcsr_type)                                   :: matrix

      defect = 0.0_dp
      peak = 0
      status = -1
      invalid = 0
      IF (order < 2 .OR. MOD(order, 2) /= 0 .OR. limit <= 0) invalid = 1
      IF (SIZE(rows) /= SIZE(cols) .OR. SIZE(rows) /= SIZE(values)) invalid = 1
      IF (ANY(rows < 1) .OR. ANY(rows > order) .OR. ANY(cols < 1) .OR. ANY(cols > order)) invalid = 1
      IF (.NOT. ALL(ieee_is_finite(values)) .OR. .NOT. ieee_is_finite(tolerance)) invalid = 1
      IF (tolerance <= 0.0_dp .OR. tolerance >= 1.0_dp) invalid = 1
      CALL para_env%max(invalid)
      IF (invalid /= 0) RETURN
      dimensions = order
      limits = limit
      tolerances = tolerance
      CALL para_env%min(dimensions(1))
      CALL para_env%max(dimensions(2))
      CALL para_env%min(limits(1))
      CALL para_env%max(limits(2))
      CALL para_env%min(tolerances(1))
      CALL para_env%max(tolerances(2))
      IF (dimensions(1) /= dimensions(2) .OR. limits(1) /= limits(2) .OR. tolerances(1) /= tolerances(2)) RETURN
      scale = 0.0_dp
      IF (SIZE(values) > 0) scale = MAXVAL(ABS(values))
      CALL para_env%max(scale)
      status = 1
      IF (scale <= TINY(1.0_dp)) RETURN
      reserved = INT(order, int_8) + tile_size**2
      status = -6
      IF (reserved > limit/8_int_8) RETURN
      ALLOCATE (canonical_rows(SIZE(rows)), canonical_cols(SIZE(rows)), canonical(SIZE(values)), row_norm(order))
      peak = 8_int_8*reserved
      canonical_rows = MAX(rows, cols)
      canonical_cols = MIN(rows, cols)
      DO pass = 1, 2
         canonical = 0.5_dp*(values/scale)
         IF (pass == 1) THEN
            WHERE (rows == cols) canonical = values/scale
         ELSE
            WHERE (rows < cols) canonical = -canonical
            WHERE (rows == cols) canonical = 0.0_dp
         END IF
         CALL import_matrix(order, canonical_rows, canonical_cols, canonical, para_env, matrix, pass == 2, &
                            limit, reserved, stored, status)
         IF (status /= 0) THEN
            CALL dbcsr_release(matrix)
            EXIT
         END IF
         peak = MAX(peak, 8_int_8*(stored + reserved))
         row_norm = 0.0_dp
         k = 0
         invalid = 0
         CALL dbcsr_iterator_start(iterator, matrix)
         DO WHILE (dbcsr_iterator_blocks_left(iterator))
            CALL dbcsr_iterator_next_block(iterator, block=block, row_offset=row, col_offset=col)
            IF (.NOT. ALL(ieee_is_finite(block))) invalid = 1
            row_norm(row:row + SIZE(block, 1) - 1) = row_norm(row:row + SIZE(block, 1) - 1) + SUM(ABS(block), DIM=2)
            IF (pass == 1) CYCLE
            DO j = 1, SIZE(block, 2)
               DO i = 1, SIZE(block, 1)
                  IF (row + i > col + j .AND. block(i, j) /= 0.0_dp) k = k + 1
               END DO
            END DO
         END DO
         CALL dbcsr_iterator_stop(iterator)
         CALL para_env%sum(row_norm)
         CALL para_env%max(invalid)
         norm = MAXVAL(row_norm)
         IF (pass == 1) THEN
            defect = norm
         ELSE
            status = -1
            IF (invalid == 0 .AND. defect <= tolerance*norm) THEN
               status = 1
               IF (norm > TINY(1.0_dp)) THEN
                  ALLOCATE (lower_rows(k), lower_cols(k), lower_values(k))
                  k = 0
                  CALL dbcsr_iterator_start(iterator, matrix)
                  DO WHILE (dbcsr_iterator_blocks_left(iterator))
                     CALL dbcsr_iterator_next_block(iterator, block=block, row_offset=row, col_offset=col)
                     DO j = 1, SIZE(block, 2)
                        DO i = 1, SIZE(block, 1)
                           IF (row + i <= col + j .OR. block(i, j) == 0.0_dp) CYCLE
                           k = k + 1
                           lower_rows(k) = row + i - 1
                           lower_cols(k) = col + j - 1
                           lower_values(k) = block(i, j)
                        END DO
                     END DO
                  END DO
                  CALL dbcsr_iterator_stop(iterator)
                  status = 0
               END IF
            END IF
         END IF
         CALL dbcsr_release(matrix)
         IF (invalid /= 0) THEN
            status = -1
            EXIT
         END IF
      END DO
      CALL para_env%max(peak)
   END SUBROUTINE localizer_skew_coo

! **************************************************************************************************
!> \brief Distributed real Pfaffian sign with complete pivoting and independent solve probes.
!> \param order Even matrix dimension
!> \param rows Distributed additive full COO rows
!> \param cols Distributed additive full COO columns
!> \param values Real skew matrix entries
!> \param para_env Communicator
!> \param tolerance Structure, pivot and relative solve tolerance
!> \param sign Pfaffian sign, zero if unresolved or rejected
!> \param residual Largest independent solve relative residual
!> \param fill Stored factor scalar entries over all ranks
!> \param status 0 success, 1 unresolved, -1 input, -4 numerical failure, -6 memory limit
!> \param limit Per-rank numerical matrix/panel/work bytes, excluding COO and metadata
!> \param peak Peak accepted numerical buffer bytes
! **************************************************************************************************
   SUBROUTINE localizer_dbcsr_pfaffian(order, rows, cols, values, para_env, tolerance, sign, residual, &
                                       fill, status, limit, peak)
      INTEGER, INTENT(IN)                                :: order, rows(:), cols(:)
      REAL(dp), INTENT(IN)                               :: values(:)
      TYPE(mp_para_env_type), POINTER                    :: para_env
      REAL(dp), INTENT(IN)                               :: tolerance
      INTEGER, INTENT(OUT)                               :: sign
      REAL(dp), INTENT(OUT)                              :: residual
      INTEGER(int_8), INTENT(OUT)                        :: fill
      INTEGER, INTENT(OUT)                               :: status
      INTEGER(int_8), INTENT(IN)                         :: limit
      INTEGER(int_8), INTENT(OUT)                        :: peak

      INTEGER                                            :: i, j, negative, pf
      INTEGER, ALLOCATABLE                               :: lower_cols(:), lower_rows(:)
      REAL(dp)                                           :: defect, rhs_norm
      REAL(dp), ALLOCATABLE                              :: lower_values(:), probes(:, :), rhs(:, :)
      TYPE(localizer_factor_type)                        :: factor

      sign = 0
      residual = 0.0_dp
      fill = 0
      CALL localizer_skew_coo(order, rows, cols, values, para_env, tolerance, lower_rows, lower_cols, &
                              lower_values, defect, peak, status, limit)
      IF (status /= 0) RETURN
      CALL localizer_factorize(order, lower_rows, lower_cols, lower_values, para_env, factor, negative, status, &
                               retain=.TRUE., skew=.TRUE., tolerance=tolerance, pf_sign=pf, max_memory_bytes=limit)
      fill = factor%entries
      peak = MAX(peak, factor%peak_memory_bytes)
      IF (status == 0) THEN
         ! Include RHS, independent probes, solve update, and pivot workspace before allocation.
         CALL reserve_buffers(factor%retained_entries + 12_int_8*order + 8, factor, para_env, status)
         IF (status == 0) THEN
            ALLOCATE (rhs(order, 4), probes(order, 4))
            DO j = 1, 4
               DO i = 1, order
                  probes(i, j) = SIN(REAL(17*i + 13*j, dp)) + COS(REAL(11*i - 7*j, dp))
               END DO
            END DO
            rhs = probes
            CALL localizer_factor_solve(factor, rhs, para_env, status)
            IF (status == 0) THEN
               DO j = 1, 4
                  rhs_norm = SQRT(SUM(probes(:, j)**2))
                  probes(:, j) = 0.0_dp
                  DO i = 1, SIZE(lower_values)
                     probes(lower_rows(i), j) = probes(lower_rows(i), j) + lower_values(i)*rhs(lower_cols(i), j)
                     probes(lower_cols(i), j) = probes(lower_cols(i), j) - lower_values(i)*rhs(lower_rows(i), j)
                  END DO
                  CALL para_env%sum(probes(:, j))
                  DO i = 1, order
                     probes(i, j) = probes(i, j) - SIN(REAL(17*i + 13*j, dp)) - COS(REAL(11*i - 7*j, dp))
                  END DO
                  residual = MAX(residual, SQRT(SUM(probes(:, j)**2))/rhs_norm)
               END DO
               IF (ieee_is_finite(residual) .AND. residual <= tolerance) THEN
                  sign = pf
               ELSE
                  status = 1
               END IF
            END IF
         END IF
      END IF
      peak = MAX(peak, factor%peak_memory_bytes)
      CALL para_env%max(peak)
      CALL localizer_factor_release(factor)
   END SUBROUTINE localizer_dbcsr_pfaffian

! **************************************************************************************************
!> \brief Collect at most two pivot columns, not the global matrix.
!> \param matrix Distributed Schur complement
!> \param pivots Scalar column indices
!> \param para_env Communicator
!> \param columns Replicated pivot columns
! **************************************************************************************************
   SUBROUTINE get_columns(matrix, pivots, para_env, columns)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      INTEGER, INTENT(IN)                                :: pivots(:)
      TYPE(mp_para_env_type), POINTER                    :: para_env
      REAL(dp), INTENT(OUT)                              :: columns(:, :)

      INTEGER                                            :: col, k, row
      REAL(dp), POINTER                                  :: block(:, :)
      TYPE(dbcsr_iterator_type)                          :: iterator

      columns = 0.0_dp
      CALL dbcsr_iterator_start(iterator, matrix)
      DO WHILE (dbcsr_iterator_blocks_left(iterator))
         CALL dbcsr_iterator_next_block(iterator, block=block, row_offset=row, col_offset=col)
         DO k = 1, SIZE(pivots)
            IF (pivots(k) < col .OR. pivots(k) >= col + SIZE(block, 2)) CYCLE
            columns(row:row + SIZE(block, 1) - 1, k) = columns(row:row + SIZE(block, 1) - 1, k) + &
                                                       block(:, pivots(k) - col + 1)
         END DO
      END DO
      CALL dbcsr_iterator_stop(iterator)
      CALL para_env%sum(columns)
   END SUBROUTINE get_columns

! **************************************************************************************************
!> \brief Dense tile updates on distributed sparse storage, retaining all numerical fill.
!> \param matrix Schur complement
!> \param order Scalar order
!> \param columns Active pivot columns
!> \param multipliers Columns times inverse pivot block
!> \param para_env Communicator defining row-tile ownership
!> \param skew Use the skew Schur complement update
! **************************************************************************************************
   SUBROUTINE update_matrix(matrix, order, columns, multipliers, para_env, skew)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      INTEGER, INTENT(IN)                                :: order
      REAL(dp), INTENT(IN)                               :: columns(:, :), multipliers(:, :)
      TYPE(mp_para_env_type), POINTER                    :: para_env
      LOGICAL, INTENT(IN)                                :: skew

      INTEGER                                            :: col, finish, first, last, nb, nc, nr, &
                                                            row, start
      LOGICAL                                            :: found
      REAL(dp)                                           :: buffer(tile_size, tile_size)
      REAL(dp), POINTER                                  :: block(:, :)

      nb = (order - 1)/tile_size + 1
      DO row = 1, nb
         IF (MOD(row - 1, para_env%num_pe) /= para_env%mepos) CYCLE
         first = (row - 1)*tile_size + 1
         last = MIN(row*tile_size, order)
         IF (.NOT. ANY(multipliers(first:last, :) /= 0.0_dp)) CYCLE
         nr = last - first + 1
         DO col = 1, nb
            start = (col - 1)*tile_size + 1
            finish = MIN(col*tile_size, order)
            IF (.NOT. ANY(columns(start:finish, :) /= 0.0_dp)) CYCLE
            nc = finish - start + 1
            buffer(:nr, :nc) = -MATMUL(multipliers(first:last, :), TRANSPOSE(columns(start:finish, :)))
            IF (skew) buffer(:nr, :nc) = -buffer(:nr, :nc)
            CALL dbcsr_get_block_p(matrix, row, col, block, found=found)
            IF (found) THEN
               block = block + buffer(:nr, :nc)
            ELSE
               CALL dbcsr_put_block(matrix, row, col, buffer(:nr, :nc))
            END IF
         END DO
      END DO
      CALL dbcsr_finalize(matrix)
   END SUBROUTINE update_matrix

! **************************************************************************************************
!> \brief Repeated solves with distributed LDL panels and replicated right hand sides.
!> \param factor Valid retained factors
!> \param rhs Right hand sides, replaced by solutions
!> \param para_env Communicator
!> \param status 0 success, -4 invalid factors or nonfinite solution
! **************************************************************************************************
   SUBROUTINE localizer_factor_solve(factor, rhs, para_env, status)
      TYPE(localizer_factor_type), INTENT(IN)            :: factor
      REAL(dp), INTENT(INOUT)                            :: rhs(:, :)
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(OUT)                               :: status

      INTEGER                                            :: i, invalid, j, k, width
      REAL(dp), ALLOCATABLE                              :: pivot_rhs(:, :), update(:, :)

      status = -4
      invalid = 0
      IF (.NOT. factor%valid .OR. SIZE(rhs, 1) /= factor%order .OR. SIZE(rhs, 2) < 1) invalid = 1
      IF (.NOT. ALL(ieee_is_finite(rhs))) invalid = 1
      CALL para_env%max(invalid)
      IF (invalid /= 0) RETURN
      invalid = MERGE(1, 0, factor%retained_entries + &
                      2_int_8*SIZE(rhs, 2)*(INT(factor%order, int_8) + 1) > factor%memory_limit/8_int_8)
      CALL para_env%max(invalid)
      IF (invalid /= 0) THEN
         status = -6
         RETURN
      END IF
      ALLOCATE (update(SIZE(rhs, 1), SIZE(rhs, 2)), pivot_rhs(2, SIZE(rhs, 2)))
      rhs = rhs/factor%scale
      DO k = 1, factor%steps
         width = factor%panels(k)%width
         pivot_rhs(:width, :) = rhs(factor%panels(k)%pivot(:width), :)
         update = 0.0_dp
         DO i = 1, SIZE(factor%panels(k)%rows)
            DO j = 1, width
               update(factor%panels(k)%rows(i), :) = update(factor%panels(k)%rows(i), :) + &
                                                     factor%panels(k)%values(i, j)*pivot_rhs(j, :)
            END DO
         END DO
         CALL para_env%sum(update)
         rhs = rhs - update
      END DO
      DO k = 1, factor%steps
         width = factor%panels(k)%width
         pivot_rhs(:width, :) = MATMUL(factor%panels(k)%inverse(:width, :width), rhs(factor%panels(k)%pivot(:width), :))
         rhs(factor%panels(k)%pivot(:width), :) = pivot_rhs(:width, :)
      END DO
      DO k = factor%steps, 1, -1
         width = factor%panels(k)%width
         pivot_rhs = 0.0_dp
         DO i = 1, SIZE(factor%panels(k)%rows)
            pivot_rhs(:width, :) = pivot_rhs(:width, :) + &
                                   MATMUL(RESHAPE(factor%panels(k)%values(i, :), [width, 1]), &
                                          rhs(factor%panels(k)%rows(i):factor%panels(k)%rows(i), :))
         END DO
         CALL para_env%sum(pivot_rhs)
         rhs(factor%panels(k)%pivot(:width), :) = rhs(factor%panels(k)%pivot(:width), :) - pivot_rhs(:width, :)
      END DO
      invalid = 0
      IF (.NOT. ALL(ieee_is_finite(rhs))) invalid = 1
      CALL para_env%max(invalid)
      IF (invalid /= 0) RETURN
      status = 0
   END SUBROUTINE localizer_factor_solve

! **************************************************************************************************
!> \brief Release elimination panels, including a partially failed factorization.
!> \param factor Factor context
! **************************************************************************************************
   SUBROUTINE localizer_factor_release(factor)
      TYPE(localizer_factor_type), INTENT(INOUT)         :: factor

      IF (ALLOCATED(factor%panels)) DEALLOCATE (factor%panels)
      factor%order = 0
      factor%steps = 0
      factor%scale = 0.0_dp
      factor%peak_memory_mb = 0.0_dp
      factor%entries = 0
      factor%workspace_entries = 0
      factor%peak_memory_bytes = 0
      factor%retained_entries = 0
      factor%memory_limit = HUGE(0_int_8)
      factor%valid = .FALSE.
   END SUBROUTINE localizer_factor_release

END MODULE localizer_dbcsr_factor
