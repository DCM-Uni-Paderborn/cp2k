!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Reusable approximate inverse models for experimental OT preconditioners.
! **************************************************************************************************
MODULE ot_preconditioner_models
   USE fft_tools,                       ONLY: BWFFT,&
                                              FFT_RADIX_NEXT,&
                                              FWFFT,&
                                              fft3d,&
                                              fft_radix_operations,&
                                              fft_type
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: pi,&
                                              twopi
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ot_preconditioner_models'

   PUBLIC :: apply_chebyshev_inverse_dense, &
             apply_helmholtz_grid_inverse, &
             apply_two_level_inverse_dense

CONTAINS

! **************************************************************************************************
!> \brief Apply a Chebyshev approximation to A^-1 for a Hermitian positive-definite dense operator.
!> \param matrix ...
!> \param rhs ...
!> \param lambda_min ...
!> \param lambda_max ...
!> \param degree ...
!> \param RESULT ...
! **************************************************************************************************
   SUBROUTINE apply_chebyshev_inverse_dense(matrix, rhs, lambda_min, lambda_max, degree, RESULT)

      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: matrix, rhs
      REAL(KIND=dp), INTENT(IN)                          :: lambda_min, lambda_max
      INTEGER, INTENT(IN)                                :: degree
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(OUT)     :: result

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: t_next, t_now, t_previous
      INTEGER                                            :: icoeff, iquad, n, nquad
      REAL(KIND=dp)                                      :: half_width, midpoint, theta, x
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: coeff

      n = SIZE(matrix, 1)
      CPASSERT(SIZE(matrix, 2) == n)
      CPASSERT(SIZE(rhs, 1) == n)
      CPASSERT(ALL(SHAPE(RESULT) == SHAPE(rhs)))
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

      ALLOCATE (t_previous(SIZE(rhs, 1), SIZE(rhs, 2)), &
                t_now(SIZE(rhs, 1), SIZE(rhs, 2)), &
                t_next(SIZE(rhs, 1), SIZE(rhs, 2)))
      t_previous = rhs
      RESULT = 0.5_dp*coeff(0)*t_previous
      IF (degree >= 1) THEN
         t_now = (MATMUL(matrix, rhs) - midpoint*rhs)/half_width
         RESULT = RESULT + coeff(1)*t_now
         DO icoeff = 2, degree
            t_next = 2.0_dp*(MATMUL(matrix, t_now) - midpoint*t_now)/half_width - t_previous
            RESULT = RESULT + coeff(icoeff)*t_next
            t_previous = t_now
            t_now = t_next
         END DO
      END IF

      DEALLOCATE (coeff, t_next, t_now, t_previous)

   END SUBROUTINE apply_chebyshev_inverse_dense

! **************************************************************************************************
!> \brief Apply a local polynomial inverse followed by a Galerkin coarse-space correction.
!>        The columns of prolongation may represent long-wavelength grid modes, translated
!>        reference-cell orbitals, or any other compact collective subspace.  The correction
!>        is residual based, so the polynomial and coarse models are not double counted.
!> \param matrix ...
!> \param rhs ...
!> \param lambda_min ...
!> \param lambda_max ...
!> \param degree ...
!> \param prolongation ...
!> \param coarse_matrix ...
!> \param RESULT ...
!> \param coarse_applied ...
! **************************************************************************************************
   SUBROUTINE apply_two_level_inverse_dense(matrix, rhs, lambda_min, lambda_max, degree, &
                                            prolongation, coarse_matrix, RESULT, coarse_applied)

      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: matrix, rhs
      REAL(KIND=dp), INTENT(IN)                          :: lambda_min, lambda_max
      INTEGER, INTENT(IN)                                :: degree
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: prolongation, coarse_matrix
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(OUT)     :: result
      LOGICAL, INTENT(OUT), OPTIONAL                     :: coarse_applied

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: coarse_rhs, coarse_work, residual
      INTEGER                                            :: info, n, ncoarse, nrhs

      n = SIZE(matrix, 1)
      nrhs = SIZE(rhs, 2)
      ncoarse = SIZE(prolongation, 2)
      CPASSERT(SIZE(matrix, 2) == n)
      CPASSERT(SIZE(rhs, 1) == n)
      CPASSERT(ALL(SHAPE(RESULT) == SHAPE(rhs)))
      CPASSERT(SIZE(prolongation, 1) == n)
      CPASSERT(SIZE(coarse_matrix, 1) == ncoarse)
      CPASSERT(SIZE(coarse_matrix, 2) == ncoarse)

      CALL apply_chebyshev_inverse_dense(matrix, rhs, lambda_min, lambda_max, degree, RESULT)
      IF (ncoarse == 0) THEN
         IF (PRESENT(coarse_applied)) coarse_applied = .FALSE.
         RETURN
      END IF

      ALLOCATE (residual(n, nrhs), coarse_rhs(ncoarse, nrhs), coarse_work(ncoarse, ncoarse))
      residual = rhs - MATMUL(matrix, RESULT)
      coarse_rhs = MATMUL(CONJG(TRANSPOSE(prolongation)), residual)
      coarse_work = 0.5_dp*(coarse_matrix + CONJG(TRANSPOSE(coarse_matrix)))
      CALL zposv('U', ncoarse, nrhs, coarse_work, ncoarse, coarse_rhs, ncoarse, info)
      IF (info == 0) RESULT = RESULT + MATMUL(prolongation, coarse_rhs)
      IF (PRESENT(coarse_applied)) coarse_applied = info == 0
      DEALLOCATE (coarse_rhs, coarse_work, residual)

   END SUBROUTINE apply_two_level_inverse_dense

! **************************************************************************************************
!> \brief Apply (-nabla^2 + shift)^-1 on a periodic orthorhombic coarse grid.
!>        Unsupported serial FFT lengths use the positive diagonal fallback rhs/shift.
!> \param rhs ...
!> \param cell_lengths ...
!> \param shift ...
!> \param RESULT ...
!> \param used_fft ...
! **************************************************************************************************
   SUBROUTINE apply_helmholtz_grid_inverse(rhs, cell_lengths, shift, RESULT, used_fft)

      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: rhs
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: cell_lengths
      REAL(KIND=dp), INTENT(IN)                          :: shift
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(OUT)     :: result
      LOGICAL, INTENT(OUT), OPTIONAL                     :: used_fft

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)  :: spectral, work
      INTEGER                                            :: d, i, j, k, m(3), n(3), radix_length, &
                                                            stat
      LOGICAL                                            :: compatible
      REAL(KIND=dp)                                      :: denominator

      CPASSERT(ALL(SHAPE(RESULT) == SHAPE(rhs)))
      CPASSERT(ALL(cell_lengths > 0.0_dp))
      CPASSERT(shift > 0.0_dp)

      n = SHAPE(rhs)
      compatible = .TRUE.
      IF (fft_type /= 3) THEN
         DO d = 1, 3
            IF (n(d) < 3) compatible = .FALSE.
            IF (compatible) THEN
               CALL fft_radix_operations(n(d), radix_length, FFT_RADIX_NEXT)
               IF (radix_length /= n(d)) compatible = .FALSE.
            END IF
         END DO
      END IF
      IF (.NOT. compatible) THEN
         RESULT = rhs/shift
         IF (PRESENT(used_fft)) used_fft = .FALSE.
         RETURN
      END IF

      ALLOCATE (spectral(n(1), n(2), n(3)), work(n(1), n(2), n(3)))
      work = CMPLX(rhs, 0.0_dp, KIND=dp)
      CALL fft3d(FWFFT, n, work, spectral, status=stat)
      IF (stat == 0) THEN
         DO k = 1, n(3)
            m(3) = k - 1
            IF (m(3) > n(3)/2) m(3) = m(3) - n(3)
            DO j = 1, n(2)
               m(2) = j - 1
               IF (m(2) > n(2)/2) m(2) = m(2) - n(2)
               DO i = 1, n(1)
                  m(1) = i - 1
                  IF (m(1) > n(1)/2) m(1) = m(1) - n(1)
                  denominator = shift + SUM((twopi*REAL(m, KIND=dp)/cell_lengths)**2)
                  spectral(i, j, k) = spectral(i, j, k)/denominator
               END DO
            END DO
         END DO
         CALL fft3d(BWFFT, n, spectral, work, status=stat)
      END IF
      IF (stat == 0) THEN
         RESULT = REAL(work, KIND=dp)
      ELSE
         RESULT = rhs/shift
         compatible = .FALSE.
      END IF
      IF (PRESENT(used_fft)) used_fft = compatible
      DEALLOCATE (spectral, work)

   END SUBROUTINE apply_helmholtz_grid_inverse

END MODULE ot_preconditioner_models
