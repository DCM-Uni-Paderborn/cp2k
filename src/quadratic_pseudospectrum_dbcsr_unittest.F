!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Distributed non-Hermitian products and metric-aware quadratic translations.
!> \author Thomas D. Kuehne, tkuehne@cp2k.org
! **************************************************************************************************
PROGRAM quadratic_pseudospectrum_dbcsr_unittest
   USE bibliography,                    ONLY: add_all_references
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_create, dbcsr_distribution_new, dbcsr_distribution_release, dbcsr_distribution_type, &
        dbcsr_finalize, dbcsr_finalize_lib, dbcsr_get_block_p, dbcsr_init_lib, dbcsr_put_block, &
        dbcsr_type_no_symmetry
   USE cp_log_handling,                 ONLY: cp_add_default_logger,&
                                              cp_logger_create,&
                                              cp_logger_release,&
                                              cp_logger_type,&
                                              cp_rm_default_logger
   USE ieee_arithmetic,                 ONLY: ieee_quiet_nan,&
                                              ieee_value
   USE kinds,                           ONLY: dp
   USE localizer_sparse_inertia,        ONLY: sparse_inertia_available,&
                                              sparse_inertia_finish,&
                                              sparse_inertia_type,&
                                              sparse_metric_solve,&
                                              sparse_metric_start
   USE machine,                         ONLY: default_output_unit
   USE message_passing,                 ONLY: mp_para_env_type,&
                                              mp_world_finalize,&
                                              mp_world_init
   USE quadratic_pseudospectrum,        ONLY: quadratic_dense,&
                                              quadratic_result_type,&
                                              quadratic_solve
   USE reference_manager,               ONLY: remove_all_references
   USE spectral_localizer_dbcsr,        ONLY: localizer_dbcsr_apply,&
                                              localizer_dbcsr_check,&
                                              localizer_dbcsr_pencil,&
                                              localizer_dbcsr_release,&
                                              localizer_matrix_type

   IMPLICIT NONE
   INTEGER, PARAMETER :: n = 18, nb = 6, bs = 3
   TYPE(mp_para_env_type), POINTER :: para_env
   TYPE(cp_logger_type), POINTER :: logger
   TYPE(dbcsr_distribution_type) :: distribution
   TYPE(localizer_matrix_type) :: operators(2), metric, invalid_operator
   TYPE(sparse_inertia_type) :: factors
   TYPE(quadratic_result_type) :: reference, trial
   INTEGER, POINTER :: process_grid(:, :), row_dist(:), col_dist(:), block_size(:)
   INTEGER, ALLOCATABLE :: rows(:), cols(:)
   REAL(dp), ALLOCATABLE :: a(:), b(:), diagonal(:)
   REAL(dp), POINTER :: block(:, :)
   COMPLEX(dp) :: dense(n, n, 2), overlap(n, n), transform(n, n), cholesky(n, n), &
                  vectors(n, 3), image(n, 3), reverse(n, 3), query(2)
   REAL(dp) :: weights(2), ratio, error_forward, error_adjoint, error_metric, error_gap
   REAL(dp) :: error_hermiticity, expected_hermiticity
   INTEGER :: i, j, mode, io, order, status, nr, nc
   CHARACTER(LEN=32) :: argument
   LOGICAL :: use_sparse, column_grid, valid, found

   ALLOCATE (para_env)
   CALL mp_world_init(para_env)
   io = -1
   IF (para_env%is_source()) io = default_output_unit
   CALL cp_logger_create(logger, para_env=para_env, default_global_unit_nr=io, close_global_unit_on_dealloc=.FALSE.)
   CALL cp_add_default_logger(logger)
   CALL add_all_references()
   CALL dbcsr_init_lib(para_env%get_handle(), io)
   CALL get_command_argument(1, argument)
   column_grid = TRIM(argument) == '--column-grid'
   nr = para_env%num_pe
   nc = 1
   IF (column_grid) THEN
      nr = 1
      nc = para_env%num_pe
   END IF
   ALLOCATE (process_grid(0:nr - 1, 0:nc - 1), row_dist(nb), col_dist(nb), block_size(nb))
   process_grid = RESHAPE([(i, i=0, para_env%num_pe - 1)], [nr, nc])
   DO i = 1, nb
      row_dist(i) = MODULO(i - 1, nr)
      col_dist(i) = MODULO(i - 1, nc)
   END DO
   block_size = bs
   CALL dbcsr_distribution_new(distribution, group=para_env%get_handle(), pgrid=process_grid, &
                               row_dist=row_dist, col_dist=col_dist)
   dense = 0
   transform = 0
   DO i = 1, n
      transform(i, i) = 1 + 0.01_dp*i
      dense(i, i, 1) = 0.05_dp*i
      IF (i == n) CYCLE
      transform(i, i + 1) = CMPLX(0.13_dp, -0.07_dp, dp)
      dense(i, i + 1, 2) = CMPLX(0.9_dp, 0.03_dp*i, dp)
      dense(i, i + 1, 1) = CMPLX(-0.8_dp, 0.12_dp, dp)
      dense(i + 1, i, 1) = CONJG(dense(i, i + 1, 1))
   END DO
   overlap = MATMUL(CONJG(TRANSPOSE(transform)), transform)
   DO i = 1, 2
      dense(:, :, i) = MATMUL(CONJG(TRANSPOSE(transform)), MATMUL(dense(:, :, i), transform))
      CALL make_operator(dense(:, :, i), operators(i))
   END DO
   CALL make_operator(overlap, metric)
   CALL localizer_dbcsr_check(operators(1), para_env, 1.e-12_dp, valid, error_hermiticity)
   IF (.NOT. valid .OR. error_hermiticity > 1.e-12_dp) ERROR STOP 'Hermitian operator check failed'
   CALL localizer_dbcsr_check(operators(2), para_env, 1.e-12_dp, valid, error_hermiticity)
   expected_hermiticity = MAX(MAXVAL(ABS(REAL(dense(:, :, 2) - TRANSPOSE(dense(:, :, 2)), dp))), &
                              MAXVAL(ABS(AIMAG(dense(:, :, 2) + TRANSPOSE(dense(:, :, 2))))))/ &
                          MAX(MAXVAL(ABS(REAL(dense(:, :, 2), dp))), MAXVAL(ABS(AIMAG(dense(:, :, 2)))))
   IF (valid .OR. ABS(error_hermiticity - expected_hermiticity) > 1.e-13_dp) THEN
      ERROR STOP 'Non-Hermitian operator diagnostic failed'
   END IF
   CALL make_operator(overlap, invalid_operator)
   CALL dbcsr_get_block_p(invalid_operator%real_part, row=1, col=1, block=block, found=found)
   IF (found) block(1, 1) = ieee_value(0.0_dp, ieee_quiet_nan)
   CALL localizer_dbcsr_check(invalid_operator, para_env, 1.e-12_dp, valid, error_hermiticity)
   IF (valid .OR. error_hermiticity /= HUGE(1.0_dp)) ERROR STOP 'Nonfinite operator diagnostic failed'
   CALL localizer_dbcsr_release(invalid_operator)
   DO j = 1, 3
      DO i = 1, n
         vectors(i, j) = CMPLX(SIN(REAL(i*j, dp)), COS(REAL(i + j*j, dp)), dp)
      END DO
   END DO
   CALL localizer_dbcsr_apply(operators(2), vectors, image, para_env)
   CALL localizer_dbcsr_apply(operators(2), vectors, reverse, para_env, adjoint=.TRUE.)
   error_forward = MAXVAL(ABS(image - MATMUL(dense(:, :, 2), vectors)))
   error_adjoint = MAXVAL(ABS(reverse - MATMUL(CONJG(TRANSPOSE(dense(:, :, 2))), vectors)))
   IF (MAX(error_forward, error_adjoint) > 1.e-12_dp) ERROR STOP 'Distributed complex adjoint failed'
   query = [CMPLX(0.21_dp, 0.0_dp, dp), EXP(CMPLX(0.0_dp, 0.73_dp, dp))]
   weights = [1.0_dp, 1.1_dp]
   CALL quadratic_dense(dense, overlap, query, weights, 2, .FALSE., n, 100, 1.e-10_dp, &
                        1.e-10_dp, reference, status)
   IF (status /= 0) ERROR STOP 'Distributed quadratic oracle failed'
   cholesky = overlap
   CALL zpotrf('L', n, cholesky, n, status)
   IF (status /= 0) ERROR STOP 'Unit metric factor failed'
   DO mode = 0, 1
      use_sparse = mode == 1
      IF (use_sparse) THEN
         IF (.NOT. sparse_inertia_available()) CYCLE
         CALL localizer_dbcsr_pencil(metric, metric, para_env, order, rows, cols, a, b, diagonal)
         CALL sparse_metric_start(factors, order, rows, cols, b, diagonal, 1.e-10_dp, para_env, ratio, status)
         IF (status /= 0) ERROR STOP 'Distributed quadratic MUMPS metric failed'
      END IF
      CALL solve(vectors, image, status)
      IF (status /= 0) ERROR STOP 'Distributed metric solve failed'
      error_metric = MAXVAL(ABS(MATMUL(overlap, image) - vectors))
      IF (error_metric > 1.e-11_dp) ERROR STOP 'Distributed complex metric residual failed'
      CALL quadratic_solve(n, query, weights, 2, .TRUE., 10, 2000, 1.e-10_dp, &
                           apply, solve, trial, status, adjoint)
      IF (status /= 0) ERROR STOP 'Distributed quadratic iteration failed'
      error_gap = MAXVAL(ABS(trial%eigenvalues - reference%eigenvalues))
      IF (error_gap > 1.e-10_dp) ERROR STOP 'Distributed non-Hermitian spectrum mismatch'
      IF (MAXVAL(ABS(trial%complex_means - reference%complex_means)) > 1.e-8_dp) THEN
         ERROR STOP 'Distributed complex operator expectations failed'
      END IF
      IF (para_env%is_source()) PRINT '(A,I1,A,5ES15.6)', 'Quadratic DBCSR metric mode ', mode, &
         ' forward/adjoint/solve/spectrum/Ritz: ', error_forward, error_adjoint, error_metric, &
         error_gap, MAXVAL(trial%residuals)
      CALL quadratic_solve(n, query, weights, 2, .TRUE., 10, 2000, 1.e-10_dp, &
                           apply, solve, trial, status, apply)
      IF (status /= -4) ERROR STOP 'Incorrect non-Hermitian adjoint was accepted'
      IF (use_sparse) THEN
         CALL sparse_inertia_finish(factors)
         DEALLOCATE (rows, cols, a, b, diagonal)
      END IF
   END DO
   DO i = 1, 2
      CALL localizer_dbcsr_release(operators(i))
   END DO
   CALL localizer_dbcsr_release(metric)
   CALL dbcsr_distribution_release(distribution)
   DEALLOCATE (process_grid, row_dist, col_dist, block_size)
   CALL dbcsr_finalize_lib()
   CALL remove_all_references()
   CALL cp_rm_default_logger()
   CALL cp_logger_release(logger)
   DEALLOCATE (para_env)
   CALL mp_world_finalize()
CONTAINS
! **************************************************************************************************
!> \brief Fill only locally owned atom blocks; both complex components have full storage.
!> \param dense Dense oracle
!> \param matrix Distributed operator
! **************************************************************************************************
   SUBROUTINE make_operator(dense, matrix)
      COMPLEX(dp), INTENT(IN)                            :: dense(:, :)
      TYPE(localizer_matrix_type), INTENT(OUT)           :: matrix

      INTEGER                                            :: finish, first, i, j, last, start

      CALL dbcsr_create(matrix%real_part, 'QUADRATIC REAL', distribution, dbcsr_type_no_symmetry, block_size, block_size)
      CALL dbcsr_create(matrix%imaginary_part, 'QUADRATIC IMAGINARY', distribution, &
                        dbcsr_type_no_symmetry, block_size, block_size)
      DO i = 1, nb
         first = bs*(i - 1) + 1
         last = bs*i
         DO j = 1, nb
            IF (process_grid(row_dist(i), col_dist(j)) /= para_env%mepos) CYCLE
            start = bs*(j - 1) + 1
            finish = bs*j
            IF (MAXVAL(ABS(dense(first:last, start:finish))) == 0.0_dp) CYCLE
            CALL dbcsr_put_block(matrix%real_part, i, j, REAL(dense(first:last, start:finish), dp))
            CALL dbcsr_put_block(matrix%imaginary_part, i, j, AIMAG(dense(first:last, start:finish)))
         END DO
      END DO
      CALL dbcsr_finalize(matrix%real_part)
      CALL dbcsr_finalize(matrix%imaginary_part)
   END SUBROUTINE make_operator
! **************************************************************************************************
!> \brief Apply the metric or a distributed forward operator.
!> \param OPERATOR Zero for metric, positive observable index
!> \param v Input vectors
!> \param y Output vectors
!> \param status Return status
! **************************************************************************************************
   SUBROUTINE apply(OPERATOR, v, y, status)
      INTEGER, INTENT(IN)                                :: operator
      COMPLEX(dp), INTENT(IN)                            :: v(:, :)
      COMPLEX(dp), INTENT(OUT)                           :: y(:, :)
      INTEGER, INTENT(OUT)                               :: status

      IF (OPERATOR == 0) THEN
         CALL localizer_dbcsr_apply(metric, v, y, para_env)
      ELSE
         CALL localizer_dbcsr_apply(operators(OPERATOR), v, y, para_env)
      END IF
      status = 0
   END SUBROUTINE apply
! **************************************************************************************************
!> \brief Apply a distributed adjoint with its independently checked conjugation.
!> \param OPERATOR Observable index
!> \param v Input vectors
!> \param y Adjoint images
!> \param status Return status
! **************************************************************************************************
   SUBROUTINE adjoint(OPERATOR, v, y, status)
      INTEGER, INTENT(IN)                                :: operator
      COMPLEX(dp), INTENT(IN)                            :: v(:, :)
      COMPLEX(dp), INTENT(OUT)                           :: y(:, :)
      INTEGER, INTENT(OUT)                               :: status

      CALL localizer_dbcsr_apply(operators(OPERATOR), v, y, para_env, adjoint=.TRUE.)
      status = 0
   END SUBROUTINE adjoint
! **************************************************************************************************
!> \brief Solve the same complex metric with dense or sparse reusable factors.
!> \param v Right hand sides
!> \param y Solutions
!> \param status Return status
! **************************************************************************************************
   SUBROUTINE solve(v, y, status)
      COMPLEX(dp), INTENT(IN)                            :: v(:, :)
      COMPLEX(dp), INTENT(OUT)                           :: y(:, :)
      INTEGER, INTENT(OUT)                               :: status

      REAL(dp), ALLOCATABLE                              :: rhs(:, :)

      IF (use_sparse) THEN
         ALLOCATE (rhs(2*n, SIZE(v, 2)))
         rhs(:n, :) = REAL(v, dp)
         rhs(n + 1:, :) = AIMAG(v)
         CALL sparse_metric_solve(factors, rhs, para_env, status)
         y = CMPLX(rhs(:n, :), rhs(n + 1:, :), dp)
      ELSE
         y = v
         CALL zpotrs('L', n, SIZE(v, 2), cholesky, n, y, n, status)
      END IF
   END SUBROUTINE solve
END PROGRAM quadratic_pseudospectrum_dbcsr_unittest
