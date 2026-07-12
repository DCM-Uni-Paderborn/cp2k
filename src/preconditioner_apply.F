!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief computes preconditioners, and implements methods to apply them
!>      currently used in qs_ot
!> \par History
!>      - [UB] 2009-05-13 Adding stable approximate inverse (full and sparse)
!> \author Joost VandeVondele (09.2002)
! **************************************************************************************************
MODULE preconditioner_apply
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_gemm
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_get_info,&
                                              cp_cfm_get_submatrix,&
                                              cp_cfm_release,&
                                              cp_cfm_set_submatrix,&
                                              cp_cfm_to_fm,&
                                              cp_cfm_type,&
                                              cp_fm_to_cfm
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_copy, dbcsr_get_info, dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, &
        dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_multiply, &
        dbcsr_release, dbcsr_set, dbcsr_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr
   USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_restore
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE input_constants,                 ONLY: ot_precond_full_all,&
                                              ot_precond_full_kinetic,&
                                              ot_precond_full_single,&
                                              ot_precond_full_single_inverse,&
                                              ot_precond_s_inverse,&
                                              ot_precond_solver_direct,&
                                              ot_precond_solver_inv_chol,&
                                              ot_precond_solver_update
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: z_one,&
                                              z_zero
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE preconditioner_types,            ONLY: preconditioner_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'preconditioner_apply'

   PUBLIC :: apply_preconditioner_cfm_complex, apply_preconditioner_dbcsr, &
             apply_preconditioner_dbcsr_complex, apply_preconditioner_fm, &
             full_all_response_metric_dense

CONTAINS

! **************************************************************************************************
!> \brief Build the reduced column-resolved metric of a complex FULL_ALL inverse Hessian.
!> \param spectral_residual W^H R, with W the full generalized eigenvector matrix
!> \param full_evals full generalized eigenvalues associated with rows of W
!> \param occ_evals reference orbital energies associated with output columns
!> \param energy_gap lower bound applied to every spectral denominator
!> \param response_metric R^H P_j R for every output column j
!> \param column_scale optional additional positive scaling of P_j
! **************************************************************************************************
   SUBROUTINE full_all_response_metric_dense(spectral_residual, full_evals, occ_evals, energy_gap, &
                                             response_metric, column_scale)

      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: spectral_residual
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: full_evals, occ_evals
      REAL(KIND=dp), INTENT(IN)                          :: energy_gap
      COMPLEX(KIND=dp), DIMENSION(:, :, :), INTENT(OUT)  :: response_metric
      REAL(KIND=dp), DIMENSION(:), INTENT(IN), OPTIONAL  :: column_scale

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: weighted_residual
      INTEGER                                            :: imo, nfull, nmo
      REAL(KIND=dp)                                      :: scale
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: weights

      nfull = SIZE(spectral_residual, 1)
      nmo = SIZE(spectral_residual, 2)
      CPASSERT(SIZE(full_evals) == nfull)
      CPASSERT(SIZE(occ_evals) == nmo)
      CPASSERT(SIZE(response_metric, 1) == nmo)
      CPASSERT(SIZE(response_metric, 2) == nmo)
      CPASSERT(SIZE(response_metric, 3) == nmo)
      IF (PRESENT(column_scale)) THEN
         CPASSERT(SIZE(column_scale) == nmo)
      END IF

      ALLOCATE (weighted_residual(nfull, nmo), weights(nfull))
      DO imo = 1, nmo
         scale = 1.0_dp
         IF (PRESENT(column_scale)) scale = column_scale(imo)
         CPASSERT(scale >= 0.0_dp)
         weights(:) = scale/MAX(energy_gap, full_evals(:) - occ_evals(imo))
         weighted_residual(:, :) = &
            SPREAD(CMPLX(weights(:), 0.0_dp, KIND=dp), DIM=2, NCOPIES=nmo)*spectral_residual(:, :)
         response_metric(:, :, imo) = &
            MATMUL(CONJG(TRANSPOSE(spectral_residual)), weighted_residual)
         response_metric(:, :, imo) = 0.5_dp*(response_metric(:, :, imo) + &
                                              CONJG(TRANSPOSE(response_metric(:, :, imo))))
      END DO
      DEALLOCATE (weighted_residual, weights)

   END SUBROUTINE full_all_response_metric_dense

! **************************************************************************************************
!> \brief applies a previously created preconditioner to a full matrix
!> \param preconditioner_env ...
!> \param matrix_in ...
!> \param matrix_out ...
! **************************************************************************************************
   SUBROUTINE apply_preconditioner_fm(preconditioner_env, matrix_in, matrix_out)

      TYPE(preconditioner_type)                          :: preconditioner_env
      TYPE(cp_fm_type), INTENT(IN)                       :: matrix_in, matrix_out

      CHARACTER(len=*), PARAMETER :: routineN = 'apply_preconditioner_fm'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      SELECT CASE (preconditioner_env%in_use)
      CASE (0)
         CPABORT("No preconditioner in use")
      CASE (ot_precond_full_single)
         CALL apply_full_single(preconditioner_env, matrix_in, matrix_out)
      CASE (ot_precond_full_all)
         CALL apply_full_all(preconditioner_env, matrix_in, matrix_out)
      CASE (ot_precond_full_kinetic, ot_precond_full_single_inverse, ot_precond_s_inverse)
         SELECT CASE (preconditioner_env%solver)
         CASE (ot_precond_solver_inv_chol, ot_precond_solver_update)
            CALL apply_full_single(preconditioner_env, matrix_in, matrix_out)
         CASE (ot_precond_solver_direct)
            CALL apply_full_direct(preconditioner_env, matrix_in, matrix_out)
         CASE DEFAULT
            CPABORT("Solver not implemented")
         END SELECT
      CASE DEFAULT
         CPABORT("Unknown preconditioner")
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE apply_preconditioner_fm

! **************************************************************************************************
!> \brief ...
!> \param preconditioner_env ...
!> \param matrix_in ...
!> \param matrix_out ...
! **************************************************************************************************
   SUBROUTINE apply_preconditioner_dbcsr(preconditioner_env, matrix_in, matrix_out)

      TYPE(preconditioner_type)                          :: preconditioner_env
      TYPE(dbcsr_type)                                   :: matrix_in, matrix_out

      CHARACTER(len=*), PARAMETER :: routineN = 'apply_preconditioner_dbcsr'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      SELECT CASE (preconditioner_env%in_use)
      CASE (0)
         CPABORT("No preconditioner in use")
      CASE (ot_precond_full_single)
         CALL apply_single(preconditioner_env, matrix_in, matrix_out)
      CASE (ot_precond_full_all)
         CALL apply_all(preconditioner_env, matrix_in, matrix_out)
      CASE (ot_precond_full_kinetic, ot_precond_full_single_inverse, ot_precond_s_inverse)
         SELECT CASE (preconditioner_env%solver)
         CASE (ot_precond_solver_inv_chol, ot_precond_solver_update)
            CALL apply_single(preconditioner_env, matrix_in, matrix_out)
         CASE (ot_precond_solver_direct)
            CPABORT("Apply_full_direct not supported with ot")
            !CALL apply_full_direct(preconditioner_env, matrix_in, matrix_out)
         CASE DEFAULT
            CPABORT("Wrong solver")
         END SELECT
      CASE DEFAULT
         CPABORT("Wrong preconditioner")
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE apply_preconditioner_dbcsr

! **************************************************************************************************
!> \brief Apply the complex state-selective FULL_ALL spectral operator.
!>      The optional constrained branch is an exact reference for coupled
!>      REF/rotation models; the production path remains unconstrained.
!> \param preconditioner_env complex k-point preconditioner
!> \param matrix_in complex input channel
!> \param matrix_out complex output channel
!> \param matrix_constraint optional covariant constraint basis B=S*C
! **************************************************************************************************
   SUBROUTINE apply_preconditioner_cfm_complex(preconditioner_env, matrix_in, matrix_out, &
                                               matrix_constraint)

      TYPE(preconditioner_type)                          :: preconditioner_env
      TYPE(cp_cfm_type), INTENT(IN)                      :: matrix_in, matrix_out
      TYPE(cp_cfm_type), INTENT(IN), OPTIONAL            :: matrix_constraint

      CHARACTER(len=*), PARAMETER :: routineN = 'apply_preconditioner_cfm_complex'

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: constraint_dense, metric, &
                                                            metric_reference, rhs, rhs_reference, &
                                                            spectral_dense, weighted_constraint
      COMPLEX(KIND=dp), CONTIGUOUS, DIMENSION(:, :), &
         POINTER                                         :: local_data
      INTEGER                                            :: handle, i, info, j, k, n, ncol_local, &
                                                            nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp)                                      :: regularization, scale
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: metric_scale, weights
      TYPE(cp_cfm_type)                                  :: matrix_constraint_spectral, &
                                                            matrix_spectral

      CALL timeset(routineN, handle)

      CPASSERT(preconditioner_env%in_use == ot_precond_full_all)
      CPASSERT(ASSOCIATED(preconditioner_env%complex_fm))
      CPASSERT(ASSOCIATED(preconditioner_env%full_evals))
      CPASSERT(ASSOCIATED(preconditioner_env%occ_evals))
      CALL cp_cfm_get_info(matrix_in, nrow_global=n, ncol_global=k)
      CPASSERT(n == SIZE(preconditioner_env%full_evals))
      CPASSERT(k == SIZE(preconditioner_env%occ_evals))

      CALL cp_cfm_create(matrix_spectral, matrix_in%matrix_struct, &
                         name='complex FULL_ALL spectral input')
      CALL cp_cfm_gemm('C', 'N', n, k, n, z_one, preconditioner_env%complex_fm, &
                       matrix_in, z_zero, matrix_spectral)
      IF (PRESENT(matrix_constraint)) THEN
         CALL cp_cfm_get_info(matrix_constraint, nrow_global=nrow_local, ncol_global=ncol_local)
         CPASSERT(nrow_local == n)
         CPASSERT(ncol_local == k)
         CALL cp_cfm_create(matrix_constraint_spectral, matrix_in%matrix_struct, &
                            name='complex FULL_ALL spectral constraint')
         CALL cp_cfm_gemm('C', 'N', n, k, n, z_one, preconditioner_env%complex_fm, &
                          matrix_constraint, z_zero, matrix_constraint_spectral)

         ALLOCATE (constraint_dense(n, k), spectral_dense(n, k), weighted_constraint(n, k))
         ALLOCATE (metric(k, k), metric_reference(k, k), rhs(k, 1), rhs_reference(k, 1))
         ALLOCATE (metric_scale(k), weights(n))
         CALL cp_cfm_get_submatrix(matrix_constraint_spectral, constraint_dense)
         CALL cp_cfm_get_submatrix(matrix_spectral, spectral_dense)
         DO j = 1, k
            DO i = 1, n
               weights(i) = 1.0_dp/MAX(preconditioner_env%energy_gap, &
                                       preconditioner_env%full_evals(i) - &
                                       preconditioner_env%occ_evals(j))
            END DO
            weighted_constraint(:, :) = &
               SPREAD(CMPLX(weights(:), 0.0_dp, KIND=dp), DIM=2, NCOPIES=k)*constraint_dense(:, :)
            metric(:, :) = MATMUL(CONJG(TRANSPOSE(constraint_dense)), weighted_constraint)
            rhs(:, 1) = MATMUL(CONJG(TRANSPOSE(constraint_dense)), &
                               CMPLX(weights(:), 0.0_dp, KIND=dp)*spectral_dense(:, j))
            ! M=B^H P_j B is Hermitian positive definite.  Equilibrating its
            ! diagonal keeps the small dense solve stable when FULL_ALL spans
            ! eigenvalue gaps of very different size.
            metric(:, :) = 0.5_dp*(metric(:, :) + CONJG(TRANSPOSE(metric(:, :))))
            DO i = 1, k
               metric_scale(i) = 1.0_dp/SQRT(MAX(REAL(metric(i, i), KIND=dp), TINY(1.0_dp)))
            END DO
            DO i = 1, k
               metric(i, :) = metric_scale(i)*metric(i, :)*metric_scale(:)
            END DO
            rhs(:, 1) = metric_scale(:)*rhs(:, 1)
            metric_reference(:, :) = metric(:, :)
            rhs_reference(:, :) = rhs(:, :)
            CALL zposv('U', k, 1, metric, k, rhs, k, info)
            IF (info /= 0) THEN
               metric(:, :) = metric_reference(:, :)
               rhs(:, :) = rhs_reference(:, :)
               regularization = 100.0_dp*EPSILON(1.0_dp)*REAL(MAX(1, k), KIND=dp)
               DO i = 1, k
                  metric(i, i) = metric(i, i) + CMPLX(regularization, 0.0_dp, KIND=dp)
               END DO
               CALL zposv('U', k, 1, metric, k, rhs, k, info)
            END IF
            CPASSERT(info == 0)
            rhs(:, 1) = metric_scale(:)*rhs(:, 1)
            spectral_dense(:, j) = CMPLX(weights(:), 0.0_dp, KIND=dp)* &
                                   (spectral_dense(:, j) - MATMUL(constraint_dense, rhs(:, 1)))
         END DO
         CALL cp_cfm_set_submatrix(matrix_spectral, spectral_dense)
         DEALLOCATE (constraint_dense, metric, metric_reference, metric_scale, rhs, rhs_reference, &
                     spectral_dense, weighted_constraint, weights)
         CALL cp_cfm_release(matrix_constraint_spectral)
      ELSE
         CALL cp_cfm_get_info(matrix_spectral, nrow_local=nrow_local, ncol_local=ncol_local, &
                              row_indices=row_indices, col_indices=col_indices, local_data=local_data)
         DO j = 1, ncol_local
            DO i = 1, nrow_local
               scale = 1.0_dp/MAX(preconditioner_env%energy_gap, &
                                  preconditioner_env%full_evals(row_indices(i)) - &
                                  preconditioner_env%occ_evals(col_indices(j)))
               local_data(i, j) = scale*local_data(i, j)
            END DO
         END DO
      END IF
      CALL cp_cfm_gemm('N', 'N', n, k, n, z_one, preconditioner_env%complex_fm, &
                       matrix_spectral, z_zero, matrix_out)
      CALL cp_cfm_release(matrix_spectral)

      CALL timestop(handle)

   END SUBROUTINE apply_preconditioner_cfm_complex

! **************************************************************************************************
!> \brief Apply a complex state-selective FULL_ALL operator to paired real/imaginary DBCSR matrices.
!> \param preconditioner_env complex k-point preconditioner
!> \param matrix_in_re real input channel
!> \param matrix_in_im imaginary input channel
!> \param matrix_out_re real output channel
!> \param matrix_out_im imaginary output channel
!> \param matrix_constraint_re optional real part of B=S*C
!> \param matrix_constraint_im optional imaginary part of B=S*C
! **************************************************************************************************
   SUBROUTINE apply_preconditioner_dbcsr_complex(preconditioner_env, matrix_in_re, matrix_in_im, &
                                                 matrix_out_re, matrix_out_im, &
                                                 matrix_constraint_re, matrix_constraint_im)

      TYPE(preconditioner_type)                          :: preconditioner_env
      TYPE(dbcsr_type)                                   :: matrix_in_re, matrix_in_im, &
                                                            matrix_out_re, matrix_out_im
      TYPE(dbcsr_type), OPTIONAL                         :: matrix_constraint_re, &
                                                            matrix_constraint_im

      CHARACTER(len=*), PARAMETER :: routineN = 'apply_preconditioner_dbcsr_complex'

      INTEGER                                            :: handle, k, n
      TYPE(cp_cfm_type)                                  :: matrix_constraint, matrix_in, matrix_out
      TYPE(cp_fm_type) :: matrix_constraint_im_fm, matrix_constraint_re_fm, matrix_in_im_fm, &
         matrix_in_re_fm, matrix_out_im_fm, matrix_out_re_fm

      CALL timeset(routineN, handle)

      CPASSERT(preconditioner_env%in_use == ot_precond_full_all)
      CPASSERT(ASSOCIATED(preconditioner_env%complex_fm))
      CPASSERT(ASSOCIATED(preconditioner_env%full_evals))
      CPASSERT(ASSOCIATED(preconditioner_env%occ_evals))
      CALL dbcsr_get_info(matrix_in_re, nfullrows_total=n, nfullcols_total=k)
      CPASSERT(n == SIZE(preconditioner_env%full_evals))
      CPASSERT(k == SIZE(preconditioner_env%occ_evals))

      CALL cp_fm_create(matrix_in_re_fm, preconditioner_env%complex_fm%matrix_struct, &
                        nrow=n, ncol=k, name='complex FULL_ALL input real')
      CALL cp_fm_create(matrix_in_im_fm, preconditioner_env%complex_fm%matrix_struct, &
                        nrow=n, ncol=k, name='complex FULL_ALL input imaginary')
      CALL cp_fm_create(matrix_out_re_fm, preconditioner_env%complex_fm%matrix_struct, &
                        nrow=n, ncol=k, name='complex FULL_ALL output real')
      CALL cp_fm_create(matrix_out_im_fm, preconditioner_env%complex_fm%matrix_struct, &
                        nrow=n, ncol=k, name='complex FULL_ALL output imaginary')
      CALL copy_dbcsr_to_fm(matrix_in_re, matrix_in_re_fm)
      CALL copy_dbcsr_to_fm(matrix_in_im, matrix_in_im_fm)

      CALL cp_cfm_create(matrix_in, matrix_in_re_fm%matrix_struct, name='complex FULL_ALL input')
      CALL cp_cfm_create(matrix_out, matrix_in_re_fm%matrix_struct, name='complex FULL_ALL output')
      CALL cp_fm_to_cfm(matrix_in_re_fm, matrix_in_im_fm, matrix_in)

      IF (PRESENT(matrix_constraint_re) .OR. PRESENT(matrix_constraint_im)) THEN
         CPASSERT(PRESENT(matrix_constraint_re) .AND. PRESENT(matrix_constraint_im))
         CALL cp_fm_create(matrix_constraint_re_fm, preconditioner_env%complex_fm%matrix_struct, &
                           nrow=n, ncol=k, name='complex FULL_ALL constraint real')
         CALL cp_fm_create(matrix_constraint_im_fm, preconditioner_env%complex_fm%matrix_struct, &
                           nrow=n, ncol=k, name='complex FULL_ALL constraint imaginary')
         CALL copy_dbcsr_to_fm(matrix_constraint_re, matrix_constraint_re_fm)
         CALL copy_dbcsr_to_fm(matrix_constraint_im, matrix_constraint_im_fm)
         CALL cp_cfm_create(matrix_constraint, matrix_in_re_fm%matrix_struct, &
                            name='complex FULL_ALL constraint')
         CALL cp_fm_to_cfm(matrix_constraint_re_fm, matrix_constraint_im_fm, matrix_constraint)
         CALL apply_preconditioner_cfm_complex(preconditioner_env, matrix_in, matrix_out, &
                                               matrix_constraint)
         CALL cp_cfm_release(matrix_constraint)
         CALL cp_fm_release(matrix_constraint_im_fm)
         CALL cp_fm_release(matrix_constraint_re_fm)
      ELSE
         CALL apply_preconditioner_cfm_complex(preconditioner_env, matrix_in, matrix_out)
      END IF

      CALL cp_cfm_to_fm(matrix_out, matrix_out_re_fm, matrix_out_im_fm)
      CALL dbcsr_set(matrix_out_re, 0.0_dp)
      CALL dbcsr_set(matrix_out_im, 0.0_dp)
      CALL copy_fm_to_dbcsr(matrix_out_re_fm, matrix_out_re)
      CALL copy_fm_to_dbcsr(matrix_out_im_fm, matrix_out_im)

      CALL cp_cfm_release(matrix_out)
      CALL cp_cfm_release(matrix_in)
      CALL cp_fm_release(matrix_out_im_fm)
      CALL cp_fm_release(matrix_out_re_fm)
      CALL cp_fm_release(matrix_in_im_fm)
      CALL cp_fm_release(matrix_in_re_fm)

      CALL timestop(handle)

   END SUBROUTINE apply_preconditioner_dbcsr_complex

! **************************************************************************************************
!> \brief apply to full matrix, complete inversion has already been done
!> \param preconditioner_env ...
!> \param matrix_in ...
!> \param matrix_out ...
! **************************************************************************************************
   SUBROUTINE apply_full_single(preconditioner_env, matrix_in, matrix_out)

      TYPE(preconditioner_type)                          :: preconditioner_env
      TYPE(cp_fm_type), INTENT(IN)                       :: matrix_in, matrix_out

      CHARACTER(len=*), PARAMETER                        :: routineN = 'apply_full_single'

      INTEGER                                            :: handle, k, n

      CALL timeset(routineN, handle)

      CALL cp_fm_get_info(matrix_in, nrow_global=n, ncol_global=k)
      CALL parallel_gemm('N', 'N', n, k, n, 1.0_dp, preconditioner_env%fm, &
                         matrix_in, 0.0_dp, matrix_out)
      CALL timestop(handle)

   END SUBROUTINE apply_full_single

! **************************************************************************************************
!> \brief apply to dbcsr matrix, complete inversion has already been done
!> \param preconditioner_env ...
!> \param matrix_in ...
!> \param matrix_out ...
! **************************************************************************************************
   SUBROUTINE apply_single(preconditioner_env, matrix_in, matrix_out)

      TYPE(preconditioner_type)                          :: preconditioner_env
      TYPE(dbcsr_type)                                   :: matrix_in, matrix_out

      CHARACTER(len=*), PARAMETER                        :: routineN = 'apply_single'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (.NOT. ASSOCIATED(preconditioner_env%dbcsr_matrix)) &
         CPABORT("NOT ASSOCIATED preconditioner_env%dbcsr_matrix")
      CALL dbcsr_multiply('N', 'N', 1.0_dp, preconditioner_env%dbcsr_matrix, matrix_in, &
                          0.0_dp, matrix_out)

      CALL timestop(handle)

   END SUBROUTINE apply_single

! **************************************************************************************************
!> \brief preconditioner contains the factorization, application done by
!>        solving the linear system
!> \param preconditioner_env ...
!> \param matrix_in ...
!> \param matrix_out ...
! **************************************************************************************************
   SUBROUTINE apply_full_direct(preconditioner_env, matrix_in, matrix_out)

      TYPE(preconditioner_type)                          :: preconditioner_env
      TYPE(cp_fm_type), INTENT(IN)                       :: matrix_in, matrix_out

      CHARACTER(len=*), PARAMETER                        :: routineN = 'apply_full_direct'

      INTEGER                                            :: handle, k, n
      TYPE(cp_fm_type)                                   :: work

      CALL timeset(routineN, handle)

      CALL cp_fm_get_info(matrix_in, nrow_global=n, ncol_global=k)
      CALL cp_fm_create(work, matrix_in%matrix_struct, name="apply_full_single", &
                        use_sp=matrix_in%use_sp)
      CALL cp_fm_cholesky_restore(matrix_in, k, preconditioner_env%fm, work,&
           &                      "SOLVE", transa="T")
      CALL cp_fm_cholesky_restore(work, k, preconditioner_env%fm, matrix_out,&
           &                      "SOLVE", transa="N")
      CALL cp_fm_release(work)

      CALL timestop(handle)

   END SUBROUTINE apply_full_direct

! **************************************************************************************************
!> \brief full all to a full matrix
!> \param preconditioner_env ...
!> \param matrix_in ...
!> \param matrix_out ...
! **************************************************************************************************
   SUBROUTINE apply_full_all(preconditioner_env, matrix_in, matrix_out)

      TYPE(preconditioner_type)                          :: preconditioner_env
      TYPE(cp_fm_type), INTENT(IN)                       :: matrix_in, matrix_out

      CHARACTER(len=*), PARAMETER                        :: routineN = 'apply_full_all'

      INTEGER                                            :: handle, i, j, k, n, ncol_local, &
                                                            nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp)                                      :: dum
      REAL(KIND=dp), CONTIGUOUS, DIMENSION(:, :), &
         POINTER                                         :: local_data
      TYPE(cp_fm_type)                                   :: matrix_tmp

      CALL timeset(routineN, handle)

      CALL cp_fm_get_info(matrix_in, nrow_global=n, ncol_global=k)

      CALL cp_fm_create(matrix_tmp, matrix_in%matrix_struct, name="apply_full_all")
      CALL cp_fm_get_info(matrix_tmp, nrow_local=nrow_local, ncol_local=ncol_local, &
                          row_indices=row_indices, col_indices=col_indices, local_data=local_data)

      !
      CALL parallel_gemm('T', 'N', n, k, n, 1.0_dp, preconditioner_env%fm, &
                         matrix_in, 0.0_dp, matrix_tmp)

      ! do the right scaling
      DO j = 1, ncol_local
      DO i = 1, nrow_local
         dum = 1.0_dp/MAX(preconditioner_env%energy_gap, &
                          preconditioner_env%full_evals(row_indices(i)) - preconditioner_env%occ_evals(col_indices(j)))
         local_data(i, j) = local_data(i, j)*dum
      END DO
      END DO

      ! mult back
      CALL parallel_gemm('N', 'N', n, k, n, 1.0_dp, preconditioner_env%fm, &
                         matrix_tmp, 0.0_dp, matrix_out)

      CALL cp_fm_release(matrix_tmp)

      CALL timestop(handle)

   END SUBROUTINE apply_full_all

! **************************************************************************************************
!> \brief full all to a dbcsr matrix
!> \param preconditioner_env ...
!> \param matrix_in ...
!> \param matrix_out ...
! **************************************************************************************************
   SUBROUTINE apply_all(preconditioner_env, matrix_in, matrix_out)

      TYPE(preconditioner_type)                          :: preconditioner_env
      TYPE(dbcsr_type)                                   :: matrix_in, matrix_out

      CHARACTER(len=*), PARAMETER                        :: routineN = 'apply_all'

      INTEGER                                            :: col, col_offset, col_size, handle, i, j, &
                                                            row, row_offset, row_size
      REAL(KIND=dp)                                      :: dum
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: DATA
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_type)                                   :: matrix_tmp

      CALL timeset(routineN, handle)

      CALL dbcsr_copy(matrix_tmp, matrix_in, name="apply_full_all")
      CALL dbcsr_multiply('T', 'N', 1.0_dp, preconditioner_env%dbcsr_matrix, &
                          matrix_in, 0.0_dp, matrix_tmp)
      ! do the right scaling
      CALL dbcsr_iterator_start(iter, matrix_tmp)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, row, col, DATA, &
                                        row_size=row_size, col_size=col_size, &
                                        row_offset=row_offset, col_offset=col_offset)
         DO j = 1, col_size
         DO i = 1, row_size
            dum = 1.0_dp/MAX(preconditioner_env%energy_gap, &
                             preconditioner_env%full_evals(row_offset + i - 1) &
                             - preconditioner_env%occ_evals(col_offset + j - 1))
            DATA(i, j) = DATA(i, j)*dum
         END DO
         END DO
      END DO
      CALL dbcsr_iterator_stop(iter)

      ! mult back
      CALL dbcsr_multiply('N', 'N', 1.0_dp, preconditioner_env%dbcsr_matrix, &
                          matrix_tmp, 0.0_dp, matrix_out)
      CALL dbcsr_release(matrix_tmp)
      CALL timestop(handle)

   END SUBROUTINE apply_all

END MODULE preconditioner_apply
