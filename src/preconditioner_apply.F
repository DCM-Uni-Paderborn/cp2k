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
        dbcsr_add, dbcsr_copy, dbcsr_get_info, dbcsr_iterator_blocks_left, &
        dbcsr_iterator_next_block, dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, &
        dbcsr_multiply, dbcsr_release, dbcsr_set, dbcsr_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr
   USE cp_fm_basic_linalg,              ONLY: cp_fm_scale,&
                                              cp_fm_scale_and_add
   USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_restore
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE input_constants,                 ONLY: &
        ot_precond_fermi_low_rank, ot_precond_full_all, ot_precond_full_kinetic, &
        ot_precond_full_single, ot_precond_full_single_inverse, ot_precond_s_inverse, &
        ot_precond_solver_chebyshev, ot_precond_solver_direct, ot_precond_solver_inv_chol, &
        ot_precond_solver_update
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: pi,&
                                              z_one,&
                                              z_zero
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE preconditioner_types,            ONLY: preconditioner_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'preconditioner_apply'

   PUBLIC :: apply_preconditioner_cfm_complex, apply_preconditioner_dbcsr, &
             apply_preconditioner_dbcsr_complex, apply_preconditioner_fm, &
             full_all_response_metric_dbcsr_complex, full_all_response_metric_dense

CONTAINS

! **************************************************************************************************
!> \brief Chebyshev coefficients of 1/x on a positive spectral interval.
!> \param lambda_min ...
!> \param lambda_max ...
!> \param degree ...
!> \param coeff ...
! **************************************************************************************************
   SUBROUTINE chebyshev_inverse_coefficients(lambda_min, lambda_max, degree, coeff)

      REAL(KIND=dp), INTENT(IN)                          :: lambda_min, lambda_max
      INTEGER, INTENT(IN)                                :: degree
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(OUT)                                     :: coeff

      INTEGER                                            :: icoeff, iquad, nquad
      REAL(KIND=dp)                                      :: half_width, midpoint, theta, x

      CPASSERT(lambda_min > 0.0_dp .AND. lambda_max > lambda_min)
      CPASSERT(degree >= 0)
      midpoint = 0.5_dp*(lambda_max + lambda_min)
      half_width = 0.5_dp*(lambda_max - lambda_min)
      nquad = MAX(64, 4*(degree + 1))
      ALLOCATE (coeff(0:degree), source=0.0_dp)
      DO icoeff = 0, degree
         DO iquad = 1, nquad
            theta = pi*(REAL(iquad, KIND=dp) - 0.5_dp)/REAL(nquad, KIND=dp)
            x = midpoint + half_width*COS(theta)
            coeff(icoeff) = coeff(icoeff) + COS(REAL(icoeff, KIND=dp)*theta)/x
         END DO
         coeff(icoeff) = 2.0_dp*coeff(icoeff)/REAL(nquad, KIND=dp)
      END DO

   END SUBROUTINE chebyshev_inverse_coefficients

! **************************************************************************************************
!> \brief ...
!> \param preconditioner_env ...
!> \return ...
! **************************************************************************************************
   PURE LOGICAL FUNCTION is_spectral_preconditioner(preconditioner_env) RESULT(is_spectral)
      TYPE(preconditioner_type), INTENT(IN)              :: preconditioner_env

      is_spectral = preconditioner_env%in_use == ot_precond_full_all .OR. &
                    preconditioner_env%in_use == ot_precond_fermi_low_rank
   END FUNCTION is_spectral_preconditioner

! **************************************************************************************************
!> \brief Spectral inverse weight for exact FULL_ALL or its bounded Fermi-window approximation.
!> \param preconditioner_env ...
!> \param row ...
!> \param col ...
!> \return ...
! **************************************************************************************************
   PURE REAL(KIND=dp) FUNCTION spectral_inverse_weight(preconditioner_env, row, col) RESULT(weight)

      TYPE(preconditioner_type), INTENT(IN)              :: preconditioner_env
      INTEGER, INTENT(IN)                                :: row, col

      INTEGER                                            :: nocc
      REAL(KIND=dp)                                      :: delta

      delta = preconditioner_env%full_evals(row) - preconditioner_env%occ_evals(col)
      weight = 1.0_dp/MAX(preconditioner_env%energy_gap, delta)
      IF (preconditioner_env%in_use /= ot_precond_fermi_low_rank) RETURN

      nocc = SIZE(preconditioner_env%occ_evals)
      IF (row <= nocc .OR. row > nocc + preconditioner_env%spectral_rank .OR. &
          delta > preconditioner_env%spectral_window) THEN
         weight = preconditioner_env%spectral_base_scale
      END IF

   END FUNCTION spectral_inverse_weight

! **************************************************************************************************
!> \brief Build the reduced column-resolved metric of a complex FULL_ALL inverse Hessian.
!> \param spectral_residual W^H R, with W the full generalized eigenvector matrix
!> \param full_evals full generalized eigenvalues associated with rows of W
!> \param occ_evals reference orbital energies associated with output columns
!> \param energy_gap lower bound applied to every spectral denominator
!> \param response_metric R^H P_j R for every output column j
!> \param column_scale optional additional positive scaling of P_j
!> \param spectral_rank ...
!> \param spectral_window ...
!> \param spectral_base_scale ...
! **************************************************************************************************
   SUBROUTINE full_all_response_metric_dense(spectral_residual, full_evals, occ_evals, energy_gap, &
                                             response_metric, column_scale, spectral_rank, &
                                             spectral_window, spectral_base_scale)

      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: spectral_residual
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: full_evals, occ_evals
      REAL(KIND=dp), INTENT(IN)                          :: energy_gap
      COMPLEX(KIND=dp), DIMENSION(:, :, :), INTENT(OUT)  :: response_metric
      REAL(KIND=dp), DIMENSION(:), INTENT(IN), OPTIONAL  :: column_scale
      INTEGER, INTENT(IN), OPTIONAL                      :: spectral_rank
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: spectral_window, spectral_base_scale

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: weighted_residual
      INTEGER                                            :: ifull, imo, nfull, nmo
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
         IF (PRESENT(spectral_rank)) THEN
            CPASSERT(PRESENT(spectral_window) .AND. PRESENT(spectral_base_scale))
            DO ifull = 1, nfull
               IF (ifull <= nmo .OR. ifull > nmo + spectral_rank .OR. &
                   full_evals(ifull) - occ_evals(imo) > spectral_window) THEN
                  weights(ifull) = scale*spectral_base_scale
               END IF
            END DO
         END IF
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
!> \brief Build R^H P_j R directly from paired real/imaginary DBCSR residual matrices.
!>        FULL_ALL uses its column-resolved spectral weights.  FULL_SINGLE_INVERSE applies one
!>        common AO inverse to orbital gradients, but uses its associated state-selective
!>        spectral model for the coupled rotation/occupation response.
!> \param preconditioner_env complex orbital preconditioner
!> \param residual_re real part of the unweighted REF complement residual
!> \param residual_im imaginary part of the unweighted REF complement residual
!> \param response_metric reduced column-resolved FULL_ALL metric
!> \param column_scale optional additional positive scaling of P_j
! **************************************************************************************************
   SUBROUTINE full_all_response_metric_dbcsr_complex(preconditioner_env, residual_re, residual_im, &
                                                     response_metric, column_scale)

      TYPE(preconditioner_type), INTENT(IN)              :: preconditioner_env
      TYPE(dbcsr_type), INTENT(IN)                       :: residual_re, residual_im
      COMPLEX(KIND=dp), ALLOCATABLE, &
         DIMENSION(:, :, :), INTENT(INOUT)               :: response_metric
      REAL(KIND=dp), DIMENSION(:), INTENT(IN), OPTIONAL  :: column_scale

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: residual_dense, spectral_residual, &
                                                            weighted_residual
      INTEGER                                            :: imo, k, n
      LOGICAL                                            :: single_inverse, spectral_response
      REAL(KIND=dp)                                      :: scale
      TYPE(cp_cfm_type)                                  :: residual, residual_spectral
      TYPE(cp_cfm_type), POINTER                         :: response_basis
      TYPE(cp_fm_type)                                   :: residual_im_fm, residual_re_fm

      single_inverse = preconditioner_env%in_use == ot_precond_full_single_inverse
      spectral_response = is_spectral_preconditioner(preconditioner_env) .OR. &
                          (single_inverse .AND. ASSOCIATED(preconditioner_env%complex_response_fm))
      CPASSERT(single_inverse .OR. is_spectral_preconditioner(preconditioner_env))
      CPASSERT(ASSOCIATED(preconditioner_env%complex_fm))
      NULLIFY (response_basis)
      CALL dbcsr_get_info(residual_re, nfullrows_total=n, nfullcols_total=k)
      IF (spectral_response) THEN
         CPASSERT(ASSOCIATED(preconditioner_env%full_evals))
         CPASSERT(ASSOCIATED(preconditioner_env%occ_evals))
         CPASSERT(n == SIZE(preconditioner_env%full_evals))
         CPASSERT(k == SIZE(preconditioner_env%occ_evals))
         IF (single_inverse) THEN
            response_basis => preconditioner_env%complex_response_fm
         ELSE
            response_basis => preconditioner_env%complex_fm
         END IF
      ELSE
         CALL cp_cfm_get_info(preconditioner_env%complex_fm, nrow_global=imo)
         CPASSERT(n == imo)
      END IF
      IF (PRESENT(column_scale)) THEN
         CPASSERT(SIZE(column_scale) == k)
      END IF

      IF (ALLOCATED(response_metric)) THEN
         IF (ANY(SHAPE(response_metric) /= [k, k, k])) DEALLOCATE (response_metric)
      END IF
      IF (.NOT. ALLOCATED(response_metric)) ALLOCATE (response_metric(k, k, k))

      CALL cp_fm_create(residual_re_fm, preconditioner_env%complex_fm%matrix_struct, &
                        nrow=n, ncol=k, name='complex preconditioner response residual real')
      CALL cp_fm_create(residual_im_fm, preconditioner_env%complex_fm%matrix_struct, &
                        nrow=n, ncol=k, name='complex preconditioner response residual imaginary')
      CALL copy_dbcsr_to_fm(residual_re, residual_re_fm)
      CALL copy_dbcsr_to_fm(residual_im, residual_im_fm)
      CALL cp_cfm_create(residual, residual_re_fm%matrix_struct, &
                         name='complex preconditioner response residual')
      CALL cp_cfm_create(residual_spectral, residual_re_fm%matrix_struct, &
                         name='complex preconditioner weighted response residual')
      CALL cp_fm_to_cfm(residual_re_fm, residual_im_fm, residual)

      IF (.NOT. spectral_response) THEN
         ALLOCATE (residual_dense(n, k), weighted_residual(n, k))
         CALL cp_cfm_gemm('N', 'N', n, k, n, z_one, preconditioner_env%complex_fm, &
                          residual, z_zero, residual_spectral)
         CALL cp_cfm_get_submatrix(residual, residual_dense)
         CALL cp_cfm_get_submatrix(residual_spectral, weighted_residual)
         DO imo = 1, k
            scale = 1.0_dp
            IF (PRESENT(column_scale)) scale = column_scale(imo)
            CPASSERT(scale >= 0.0_dp)
            response_metric(:, :, imo) = scale* &
                                         MATMUL(CONJG(TRANSPOSE(residual_dense)), weighted_residual)
            response_metric(:, :, imo) = 0.5_dp*(response_metric(:, :, imo) + &
                                                 CONJG(TRANSPOSE(response_metric(:, :, imo))))
         END DO
         DEALLOCATE (weighted_residual, residual_dense)
      ELSE
         ALLOCATE (spectral_residual(n, k))
         CALL cp_cfm_gemm('C', 'N', n, k, n, z_one, response_basis, &
                          residual, z_zero, residual_spectral)
         CALL cp_cfm_get_submatrix(residual_spectral, spectral_residual)
         IF (.NOT. single_inverse .AND. &
             preconditioner_env%in_use == ot_precond_fermi_low_rank) THEN
            CALL full_all_response_metric_dense( &
               spectral_residual, preconditioner_env%full_evals, &
               preconditioner_env%occ_evals, preconditioner_env%energy_gap, &
               response_metric, column_scale, preconditioner_env%spectral_rank, &
               preconditioner_env%spectral_window, preconditioner_env%spectral_base_scale)
         ELSE
            CALL full_all_response_metric_dense( &
               spectral_residual, preconditioner_env%full_evals, &
               preconditioner_env%occ_evals, preconditioner_env%energy_gap, &
               response_metric, column_scale)
         END IF
         DEALLOCATE (spectral_residual)
      END IF

      CALL cp_cfm_release(residual_spectral)
      CALL cp_cfm_release(residual)
      CALL cp_fm_release(residual_im_fm)
      CALL cp_fm_release(residual_re_fm)

   END SUBROUTINE full_all_response_metric_dbcsr_complex

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
      CASE (ot_precond_full_all, ot_precond_fermi_low_rank)
         CALL apply_full_all(preconditioner_env, matrix_in, matrix_out)
      CASE (ot_precond_full_kinetic, ot_precond_full_single_inverse, ot_precond_s_inverse)
         SELECT CASE (preconditioner_env%solver)
         CASE (ot_precond_solver_inv_chol, ot_precond_solver_update)
            CALL apply_full_single(preconditioner_env, matrix_in, matrix_out)
         CASE (ot_precond_solver_direct)
            CALL apply_full_direct(preconditioner_env, matrix_in, matrix_out)
         CASE (ot_precond_solver_chebyshev)
            CALL apply_chebyshev_fm(preconditioner_env, matrix_in, matrix_out)
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
      CASE (ot_precond_full_all, ot_precond_fermi_low_rank)
         CALL apply_all(preconditioner_env, matrix_in, matrix_out)
      CASE (ot_precond_full_kinetic, ot_precond_full_single_inverse, ot_precond_s_inverse)
         SELECT CASE (preconditioner_env%solver)
         CASE (ot_precond_solver_inv_chol, ot_precond_solver_update)
            CALL apply_single(preconditioner_env, matrix_in, matrix_out)
         CASE (ot_precond_solver_direct)
            CPABORT("Apply_full_direct not supported with ot")
            !CALL apply_full_direct(preconditioner_env, matrix_in, matrix_out)
         CASE (ot_precond_solver_chebyshev)
            CALL apply_chebyshev_dbcsr(preconditioner_env, matrix_in, matrix_out)
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
                                                            nrow_local, rank_used
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp)                                      :: regularization, scale
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: metric_scale, weights
      TYPE(cp_cfm_type)                                  :: matrix_constraint_spectral, &
                                                            matrix_low_rank, matrix_spectral

      CALL timeset(routineN, handle)

      IF (preconditioner_env%in_use == ot_precond_full_single_inverse) THEN
         CPASSERT(ASSOCIATED(preconditioner_env%complex_fm))
         CPASSERT(.NOT. PRESENT(matrix_constraint))
         CALL cp_cfm_get_info(matrix_in, nrow_global=n, ncol_global=k)
         CALL cp_cfm_gemm('N', 'N', n, k, n, z_one, preconditioner_env%complex_fm, &
                          matrix_in, z_zero, matrix_out)
         CALL timestop(handle)
         RETURN
      END IF

      CPASSERT(is_spectral_preconditioner(preconditioner_env))
      CPASSERT(ASSOCIATED(preconditioner_env%complex_fm))
      CPASSERT(ASSOCIATED(preconditioner_env%full_evals))
      CPASSERT(ASSOCIATED(preconditioner_env%occ_evals))
      CALL cp_cfm_get_info(matrix_in, nrow_global=n, ncol_global=k)
      CPASSERT(n == SIZE(preconditioner_env%full_evals))
      CPASSERT(k == SIZE(preconditioner_env%occ_evals))

      IF (preconditioner_env%in_use == ot_precond_fermi_low_rank .AND. &
          .NOT. PRESENT(matrix_constraint)) THEN
         CPASSERT(ASSOCIATED(preconditioner_env%complex_base_fm))
         CALL cp_cfm_gemm('N', 'N', n, k, n, &
                          CMPLX(preconditioner_env%spectral_base_scale, 0.0_dp, KIND=dp), &
                          preconditioner_env%complex_base_fm, matrix_in, z_zero, matrix_out)
         rank_used = MIN(preconditioner_env%spectral_rank, n - k)
         IF (rank_used > 0) THEN
            CALL cp_cfm_create(matrix_low_rank, matrix_in%matrix_struct, nrow=rank_used, ncol=k, &
                               name='complex FERMI_LOW_RANK coefficients')
            CALL cp_cfm_gemm('C', 'N', rank_used, k, n, z_one, preconditioner_env%complex_fm, &
                             matrix_in, z_zero, matrix_low_rank, a_first_col=k + 1)
            CALL cp_cfm_get_info(matrix_low_rank, nrow_local=nrow_local, ncol_local=ncol_local, &
                                 row_indices=row_indices, col_indices=col_indices, local_data=local_data)
            DO j = 1, ncol_local
               DO i = 1, nrow_local
                  scale = spectral_inverse_weight(preconditioner_env, k + row_indices(i), &
                                                  col_indices(j)) - &
                          preconditioner_env%spectral_base_scale
                  local_data(i, j) = scale*local_data(i, j)
               END DO
            END DO
            CALL cp_cfm_gemm('N', 'N', n, k, rank_used, z_one, &
                             preconditioner_env%complex_fm, matrix_low_rank, z_one, matrix_out, &
                             a_first_col=k + 1)
            CALL cp_cfm_release(matrix_low_rank)
         END IF
         CALL timestop(handle)
         RETURN
      END IF

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
               weights(i) = spectral_inverse_weight(preconditioner_env, i, j)
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
               scale = spectral_inverse_weight(preconditioner_env, row_indices(i), col_indices(j))
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
!> \brief Apply a complex orbital preconditioner to paired real/imaginary DBCSR matrices.
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

      INTEGER                                            :: handle, k, n, npre
      TYPE(cp_cfm_type)                                  :: matrix_constraint, matrix_in, matrix_out
      TYPE(cp_fm_type) :: matrix_constraint_im_fm, matrix_constraint_re_fm, matrix_in_im_fm, &
         matrix_in_re_fm, matrix_out_im_fm, matrix_out_re_fm

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(preconditioner_env%complex_fm))
      CALL dbcsr_get_info(matrix_in_re, nfullrows_total=n, nfullcols_total=k)
      IF (preconditioner_env%in_use == ot_precond_full_single_inverse) THEN
         CALL cp_cfm_get_info(preconditioner_env%complex_fm, nrow_global=npre)
         CPASSERT(n == npre)
      ELSE
         CPASSERT(is_spectral_preconditioner(preconditioner_env))
         CPASSERT(ASSOCIATED(preconditioner_env%full_evals))
         CPASSERT(ASSOCIATED(preconditioner_env%occ_evals))
         CPASSERT(n == SIZE(preconditioner_env%full_evals))
         CPASSERT(k == SIZE(preconditioner_env%occ_evals))
      END IF

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

      IF (.NOT. ASSOCIATED(preconditioner_env%dbcsr_matrix)) THEN
         CPABORT("NOT ASSOCIATED preconditioner_env%dbcsr_matrix")
      END IF
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
      CALL cp_fm_create(work, matrix_in%matrix_struct, name="apply_full_single")
      CALL cp_fm_cholesky_restore(matrix_in, k, preconditioner_env%fm, work,&
           &                      "SOLVE", transa="T")
      CALL cp_fm_cholesky_restore(work, k, preconditioner_env%fm, matrix_out,&
           &                      "SOLVE", transa="N")
      CALL cp_fm_release(work)

      CALL timestop(handle)

   END SUBROUTINE apply_full_direct

! **************************************************************************************************
!> \brief Apply a Chebyshev approximation to the inverse of the stored dense SPD operator.
!> \param preconditioner_env ...
!> \param matrix_in ...
!> \param matrix_out ...
! **************************************************************************************************
   SUBROUTINE apply_chebyshev_fm(preconditioner_env, matrix_in, matrix_out)

      TYPE(preconditioner_type)                          :: preconditioner_env
      TYPE(cp_fm_type), INTENT(IN)                       :: matrix_in, matrix_out

      INTEGER                                            :: degree, icoeff, k, n
      REAL(KIND=dp)                                      :: half_width, midpoint
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: coeff
      TYPE(cp_fm_type)                                   :: t_next, t_now, t_previous

      CPASSERT(ASSOCIATED(preconditioner_env%fm))
      degree = preconditioner_env%polynomial_degree
      CALL chebyshev_inverse_coefficients(preconditioner_env%polynomial_min, &
                                          preconditioner_env%polynomial_max, degree, coeff)
      midpoint = 0.5_dp*(preconditioner_env%polynomial_max + preconditioner_env%polynomial_min)
      half_width = 0.5_dp*(preconditioner_env%polynomial_max - preconditioner_env%polynomial_min)
      CALL cp_fm_get_info(matrix_in, nrow_global=n, ncol_global=k)
      CALL cp_fm_create(t_previous, matrix_in%matrix_struct, name='Chebyshev previous vector')
      CALL cp_fm_create(t_now, matrix_in%matrix_struct, name='Chebyshev current vector')
      CALL cp_fm_create(t_next, matrix_in%matrix_struct, name='Chebyshev next vector')
      CALL cp_fm_to_fm(matrix_in, t_previous)
      CALL cp_fm_to_fm(matrix_in, matrix_out)
      CALL cp_fm_scale(0.5_dp*coeff(0), matrix_out)
      IF (degree >= 1) THEN
         CALL parallel_gemm('N', 'N', n, k, n, 1.0_dp/half_width, preconditioner_env%fm, &
                            matrix_in, 0.0_dp, t_now)
         CALL cp_fm_scale_and_add(1.0_dp, t_now, -midpoint/half_width, matrix_in)
         CALL cp_fm_scale_and_add(1.0_dp, matrix_out, coeff(1), t_now)
         DO icoeff = 2, degree
            CALL parallel_gemm('N', 'N', n, k, n, 2.0_dp/half_width, preconditioner_env%fm, &
                               t_now, 0.0_dp, t_next)
            CALL cp_fm_scale_and_add(1.0_dp, t_next, -2.0_dp*midpoint/half_width, t_now)
            CALL cp_fm_scale_and_add(1.0_dp, t_next, -1.0_dp, t_previous)
            CALL cp_fm_scale_and_add(1.0_dp, matrix_out, coeff(icoeff), t_next)
            CALL cp_fm_to_fm(t_now, t_previous)
            CALL cp_fm_to_fm(t_next, t_now)
         END DO
      END IF
      CALL cp_fm_release(t_next)
      CALL cp_fm_release(t_now)
      CALL cp_fm_release(t_previous)
      DEALLOCATE (coeff)

   END SUBROUTINE apply_chebyshev_fm

! **************************************************************************************************
!> \brief Apply a Chebyshev approximation to the inverse of the stored sparse SPD operator.
!> \param preconditioner_env ...
!> \param matrix_in ...
!> \param matrix_out ...
! **************************************************************************************************
   SUBROUTINE apply_chebyshev_dbcsr(preconditioner_env, matrix_in, matrix_out)

      TYPE(preconditioner_type)                          :: preconditioner_env
      TYPE(dbcsr_type)                                   :: matrix_in, matrix_out

      INTEGER                                            :: degree, icoeff
      REAL(KIND=dp)                                      :: half_width, midpoint
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: coeff
      TYPE(dbcsr_type)                                   :: t_next, t_now, t_previous

      CPASSERT(ASSOCIATED(preconditioner_env%dbcsr_matrix))
      degree = preconditioner_env%polynomial_degree
      CALL chebyshev_inverse_coefficients(preconditioner_env%polynomial_min, &
                                          preconditioner_env%polynomial_max, degree, coeff)
      midpoint = 0.5_dp*(preconditioner_env%polynomial_max + preconditioner_env%polynomial_min)
      half_width = 0.5_dp*(preconditioner_env%polynomial_max - preconditioner_env%polynomial_min)
      CALL dbcsr_copy(t_previous, matrix_in, name='Chebyshev previous vector')
      CALL dbcsr_set(matrix_out, 0.0_dp)
      CALL dbcsr_add(matrix_out, t_previous, 1.0_dp, 0.5_dp*coeff(0))
      IF (degree >= 1) THEN
         CALL dbcsr_copy(t_now, matrix_in, name='Chebyshev current vector')
         CALL dbcsr_multiply('N', 'N', 1.0_dp/half_width, preconditioner_env%dbcsr_matrix, &
                             matrix_in, -midpoint/half_width, t_now)
         CALL dbcsr_add(matrix_out, t_now, 1.0_dp, coeff(1))
         DO icoeff = 2, degree
            CALL dbcsr_copy(t_next, t_previous, name='Chebyshev next vector')
            CALL dbcsr_multiply('N', 'N', 2.0_dp/half_width, preconditioner_env%dbcsr_matrix, &
                                t_now, -1.0_dp, t_next)
            CALL dbcsr_add(t_next, t_now, 1.0_dp, -2.0_dp*midpoint/half_width)
            CALL dbcsr_add(matrix_out, t_next, 1.0_dp, coeff(icoeff))
            CALL dbcsr_release(t_previous)
            CALL dbcsr_copy(t_previous, t_now, name='Chebyshev previous vector')
            CALL dbcsr_release(t_now)
            CALL dbcsr_copy(t_now, t_next, name='Chebyshev current vector')
            CALL dbcsr_release(t_next)
         END DO
         CALL dbcsr_release(t_now)
      END IF
      CALL dbcsr_release(t_previous)
      DEALLOCATE (coeff)

   END SUBROUTINE apply_chebyshev_dbcsr

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
         dum = spectral_inverse_weight(preconditioner_env, row_indices(i), col_indices(j))
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
            dum = spectral_inverse_weight(preconditioner_env, row_offset + i - 1, col_offset + j - 1)
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
