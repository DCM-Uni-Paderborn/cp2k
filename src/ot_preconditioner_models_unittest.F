!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
PROGRAM ot_preconditioner_models_unittest
   USE fft_tools,                       ONLY: finalize_fft,&
                                              init_fft
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: twopi,&
                                              z_zero
   USE message_passing,                 ONLY: mp_comm_type,&
                                              mp_world_finalize,&
                                              mp_world_init
   USE ot_preconditioner_models,        ONLY: apply_chebyshev_inverse_dense,&
                                              apply_helmholtz_grid_inverse,&
                                              apply_two_level_inverse_dense
   USE preconditioner_apply,            ONLY: full_all_response_metric_dense

   IMPLICIT NONE

   COMPLEX(KIND=dp), DIMENSION(5, 5)                  :: matrix
   COMPLEX(KIND=dp), DIMENSION(5, 2)                  :: cheb_result, exact, rhs
   COMPLEX(KIND=dp), DIMENSION(5, 2)                  :: two_level_result
   COMPLEX(KIND=dp), DIMENSION(5, 2)                  :: prolongation
   COMPLEX(KIND=dp), DIMENSION(2, 2)                  :: coarse_matrix
   COMPLEX(KIND=dp), DIMENSION(2, 2, 2)               :: metric, metric_reference
   COMPLEX(KIND=dp), DIMENSION(6, 2)                  :: spectral_residual, weighted
   INTEGER                                            :: i, imo
   LOGICAL                                            :: coarse_applied, used_fft
   REAL(KIND=dp)                                      :: local_residual, two_level_residual
   REAL(KIND=dp)                                      :: weight
   REAL(KIND=dp), DIMENSION(6)                        :: full_evals
   REAL(KIND=dp), DIMENSION(2)                        :: occ_evals
   REAL(KIND=dp), DIMENSION(4, 4, 4)                  :: grid_result, grid_rhs
   REAL(KIND=dp), DIMENSION(3)                        :: cell_lengths
   TYPE(mp_comm_type)                                 :: mp_comm

   CALL mp_world_init(mp_comm)
   CALL init_fft("FFTSG", alltoall=.FALSE., fftsg_sizes=.TRUE., pool_limit=4, &
                 wisdom_file="", plan_style=1)

   matrix = z_zero
   rhs = z_zero
   DO i = 1, 5
      matrix(i, i) = CMPLX(REAL(i, KIND=dp), 0.0_dp, KIND=dp)
      rhs(i, 1) = CMPLX(SIN(0.3_dp*i), COS(0.2_dp*i), KIND=dp)
      rhs(i, 2) = CMPLX(COS(0.4_dp*i), -SIN(0.1_dp*i), KIND=dp)
      exact(i, :) = rhs(i, :)/REAL(i, KIND=dp)
   END DO
   CALL apply_chebyshev_inverse_dense(matrix, rhs, 1.0_dp, 5.0_dp, 32, cheb_result)
   IF (MAXVAL(ABS(cheb_result - exact)) > 2.0E-10_dp) THEN
      ERROR STOP "Chebyshev inverse model is inaccurate"
   END IF

   prolongation = z_zero
   prolongation(1, 1) = 1.0_dp
   prolongation(5, 2) = 1.0_dp
   coarse_matrix = z_zero
   coarse_matrix(1, 1) = matrix(1, 1)
   coarse_matrix(2, 2) = matrix(5, 5)
   CALL apply_chebyshev_inverse_dense(matrix, rhs, 1.0_dp, 5.0_dp, 0, cheb_result)
   CALL apply_two_level_inverse_dense(matrix, rhs, 1.0_dp, 5.0_dp, 0, prolongation, &
                                      coarse_matrix, two_level_result, coarse_applied)
   IF (.NOT. coarse_applied) ERROR STOP "Two-level coarse solve was not applied"
   local_residual = MAXVAL(ABS(rhs - MATMUL(matrix, cheb_result)))
   two_level_residual = MAXVAL(ABS(rhs - MATMUL(matrix, two_level_result)))
   IF (two_level_residual >= local_residual) THEN
      ERROR STOP "Two-level correction did not improve the polynomial model"
   END IF
   IF (MAXVAL(ABS(two_level_result([1, 5], :) - exact([1, 5], :))) > 2.0E-13_dp) THEN
      ERROR STOP "Two-level correction is not exact in its coarse subspace"
   END IF

   cell_lengths = [4.0_dp, 5.0_dp, 6.0_dp]
   DO i = 1, 4
      grid_rhs(i, :, :) = 3.25_dp + COS(twopi*REAL(i - 1, KIND=dp)/4.0_dp)
   END DO
   CALL apply_helmholtz_grid_inverse(grid_rhs, cell_lengths, 2.0_dp, grid_result, used_fft)
   IF (.NOT. used_fft) ERROR STOP "Coarse Helmholtz model did not use FFT"
   DO i = 1, 4
      grid_rhs(i, :, :) = 3.25_dp/2.0_dp + &
                          COS(twopi*REAL(i - 1, KIND=dp)/4.0_dp)/(2.0_dp + (twopi/cell_lengths(1))**2)
   END DO
   IF (MAXVAL(ABS(grid_result - grid_rhs)) > 2.0E-12_dp) THEN
      ERROR STOP "Coarse Helmholtz inverse differs from its Fourier reference"
   END IF

   full_evals = [-1.0_dp, -0.5_dp, 0.1_dp, 0.2_dp, 1.0_dp, 2.0_dp]
   occ_evals = [-0.5_dp, -0.4_dp]
   spectral_residual = z_zero
   DO i = 1, 6
      spectral_residual(i, 1) = CMPLX(0.1_dp*i, -0.03_dp*i, KIND=dp)
      spectral_residual(i, 2) = CMPLX(-0.04_dp*i, 0.07_dp*i, KIND=dp)
   END DO
   CALL full_all_response_metric_dense(spectral_residual, full_evals, occ_evals, 0.05_dp, &
                                       metric, spectral_rank=2, spectral_window=0.8_dp, &
                                       spectral_base_scale=1.25_dp)
   DO imo = 1, 2
      DO i = 1, 6
         IF (i > 2 .AND. i <= 4 .AND. full_evals(i) - occ_evals(imo) <= 0.8_dp) THEN
            weight = 1.0_dp/MAX(0.05_dp, full_evals(i) - occ_evals(imo))
         ELSE
            weight = 1.25_dp
         END IF
         weighted(i, :) = weight*spectral_residual(i, :)
      END DO
      metric_reference(:, :, imo) = MATMUL(CONJG(TRANSPOSE(spectral_residual)), weighted)
      metric_reference(:, :, imo) = 0.5_dp*(metric_reference(:, :, imo) + &
                                            CONJG(TRANSPOSE(metric_reference(:, :, imo))))
   END DO
   IF (MAXVAL(ABS(metric - metric_reference)) > 2.0E-13_dp) THEN
      ERROR STOP "Fermi low-rank response metric is inconsistent"
   END IF

   CALL finalize_fft(mp_comm, "")
   CALL mp_world_finalize()

END PROGRAM ot_preconditioner_models_unittest
