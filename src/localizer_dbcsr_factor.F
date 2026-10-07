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

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: localizer_factor_type, localizer_factorize, localizer_factor_solve, localizer_factor_release

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
!> \note No fill is dropped and no pivot is perturbed. Only two matrix columns are replicated.
! **************************************************************************************************
   SUBROUTINE localizer_factorize(order, rows, cols, values, para_env, factor, negative, status, retain)
      INTEGER, INTENT(IN)                                :: order, rows(:), cols(:)
      REAL(dp), INTENT(IN)                               :: values(:)
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(localizer_factor_type), INTENT(INOUT)         :: factor
      INTEGER, INTENT(OUT)                               :: negative, status
      LOGICAL, INTENT(IN)                                :: retain

      INTEGER                                            :: col, i, invalid, j, k, maximum(2), &
                                                            minimum(2), nrows, p(2), remaining, &
                                                            row, step, width
      INTEGER(int_8)                                     :: candidate(2), key, panel_entries, stored
      LOGICAL, ALLOCATABLE                               :: active(:)
      REAL(dp)                                           :: alpha, bound, d(2, 2), determinant, &
                                                            eigenvalues(2), inverse(2, 2), &
                                                            maxima(2), scale
      REAL(dp), ALLOCATABLE                              :: columns(:, :), multipliers(:, :)
      REAL(dp), POINTER                                  :: block(:, :)
      TYPE(dbcsr_iterator_type)                          :: iterator
      TYPE(dbcsr_type)                                   :: matrix

      CALL localizer_factor_release(factor)
      negative = 0
      status = -1
      invalid = 0
      IF (order < 1 .OR. SIZE(rows) /= SIZE(cols) .OR. SIZE(rows) /= SIZE(values)) invalid = 1
      IF (ANY(rows < 1) .OR. ANY(rows > order) .OR. ANY(cols < 1)) invalid = 1
      IF (SIZE(rows) == SIZE(cols)) THEN
         IF (ANY(cols > rows)) invalid = 1
      END IF
      IF (.NOT. ALL(ieee_is_finite(values))) invalid = 1
      CALL para_env%max(invalid)
      IF (invalid /= 0) RETURN
      minimum = [order, MERGE(1, 0, retain)]
      maximum = minimum
      CALL para_env%min(minimum)
      CALL para_env%max(maximum)
      IF (ANY(minimum /= maximum)) RETURN
      scale = 0.0_dp
      IF (SIZE(values) > 0) scale = MAXVAL(ABS(values))
      CALL para_env%max(scale)
      status = 1
      IF (scale <= TINY(1.0_dp)) RETURN
      factor%order = order
      factor%scale = scale
      CALL import_matrix(order, rows, cols, values/scale, para_env, matrix)
      ALLOCATE (active(order), columns(order, 2), multipliers(order, 2))
      IF (retain) ALLOCATE (factor%panels(order))
      active = .TRUE.
      remaining = order
      panel_entries = 0
      alpha = (1.0_dp + SQRT(17.0_dp))/8.0_dp
      bound = 32.0_dp*EPSILON(1.0_dp)*REAL(order, dp)
      step = 0
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
         factor%peak_memory_mb = MAX(factor%peak_memory_mb, &
                                     REAL(8_int_8*(stored + panel_entries + 4_int_8*order), dp)/1048576.0_dp)
         CALL para_env%sum(stored)
         factor%workspace_entries = MAX(factor%workspace_entries, stored)
         CALL para_env%max(invalid)
         CALL para_env%max(maxima)
         status = -4
         IF (invalid /= 0) EXIT
         status = 1
         IF (MAXVAL(maxima) <= bound) EXIT
         ! Bunch-Parlett selection bounds multipliers and forces indefinite 2x2 pivots.
         width = 1
         IF (maxima(1) < alpha*maxima(2)) width = 2
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
         IF (width == 1) THEN
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
         negative = negative + COUNT(eigenvalues(:width) < 0.0_dp)
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
         END IF
         CALL update_matrix(matrix, order, columns(:, :width), multipliers(:, :width), para_env)
         remaining = remaining - width
         status = 0
      END DO
      factor%steps = step
      CALL para_env%sum(factor%entries)
      CALL para_env%max(factor%peak_memory_mb)
      CALL dbcsr_release(matrix)
      IF (remaining > 0) THEN
         factor%valid = .FALSE.
         RETURN
      END IF
      factor%valid = retain
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
! **************************************************************************************************
   SUBROUTINE import_matrix(order, rows, cols, values, para_env, matrix)
      INTEGER, INTENT(IN)                                :: order, rows(:), cols(:)
      REAL(dp), INTENT(IN)                               :: values(:)
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_type), INTENT(OUT)                      :: matrix

      INTEGER                                            :: col, dest, i, ir, jc, mirror, nb, np, &
                                                            pass, row, slot
      INTEGER, ALLOCATABLE                               :: position(:), recv_count(:), &
                                                            recv_disp(:), recv_index(:, :), &
                                                            send_count(:), send_disp(:), &
                                                            send_index(:, :)
      INTEGER, POINTER                                   :: col_dist(:), grid(:, :), row_dist(:), &
                                                            sizes(:)
      REAL(dp)                                           :: buffer(tile_size, tile_size)
      REAL(dp), ALLOCATABLE                              :: recv_value(:), send_value(:)
      TYPE(dbcsr_distribution_type)                      :: distribution

      np = para_env%num_pe
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
            DO mirror = 1, 2
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
   END SUBROUTINE import_matrix

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
! **************************************************************************************************
   SUBROUTINE update_matrix(matrix, order, columns, multipliers, para_env)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      INTEGER, INTENT(IN)                                :: order
      REAL(dp), INTENT(IN)                               :: columns(:, :), multipliers(:, :)
      TYPE(mp_para_env_type), POINTER                    :: para_env

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
      factor%valid = .FALSE.
   END SUBROUTINE localizer_factor_release

END MODULE localizer_dbcsr_factor
