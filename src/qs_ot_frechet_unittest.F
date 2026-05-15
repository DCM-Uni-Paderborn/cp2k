!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief finite-difference checks for the complex OT rotation Frechet kernel
! **************************************************************************************************
PROGRAM qs_ot_frechet_unittest
   USE kinds,                           ONLY: dp
   USE qs_ot,                           ONLY: qs_ot_complex_exp_frechet_kernel

   IMPLICIT NONE

   INTEGER                                            :: nfail

   nfail = 0

   CALL check_offdiag("separated offdiag", 1.3_dp, -0.4_dp, &
                      CMPLX(0.7_dp, -0.2_dp, KIND=dp), 5.0E-12_dp, nfail)
   CALL check_offdiag("near-degenerate offdiag", 0.2_dp, 0.2000003_dp, &
                      CMPLX(-0.3_dp, 0.5_dp, KIND=dp), 5.0E-12_dp, nfail)
   CALL check_offdiag("degenerate offdiag", -0.6_dp, -0.6_dp, &
                      CMPLX(0.4_dp, 0.1_dp, KIND=dp), 5.0E-12_dp, nfail)
   CALL check_diag("diagonal", 0.8_dp, CMPLX(-0.2_dp, 0.35_dp, KIND=dp), &
                   5.0E-8_dp, nfail)

   IF (nfail > 0) ERROR STOP "qs_ot_frechet_unittest failed"

   PRINT *, "qs_ot_frechet_unittest: OK"

CONTAINS

! **************************************************************************************************
!> \brief Check L_exp(diag(-i*e1,-i*e2), E_12*z)
! **************************************************************************************************
   SUBROUTINE check_offdiag(label, e1, e2, z, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: e1, e2
      COMPLEX(KIND=dp), INTENT(IN)                       :: z
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      REAL(KIND=dp), PARAMETER                           :: eps = 1.0E-7_dp

      COMPLEX(KIND=dp)                                   :: fd, ref
      COMPLEX(KIND=dp), DIMENSION(2, 2)                  :: m0, meps, u0, ueps
      REAL(KIND=dp)                                      :: err

      m0(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      m0(1, 1) = CMPLX(0.0_dp, -e1, KIND=dp)
      m0(2, 2) = CMPLX(0.0_dp, -e2, KIND=dp)

      meps = m0
      meps(1, 2) = eps*z

      CALL dense_expm2(m0, u0)
      CALL dense_expm2(meps, ueps)

      fd = (ueps(1, 2) - u0(1, 2))/eps
      ref = qs_ot_complex_exp_frechet_kernel(e1, e2)*z
      err = ABS(fd - ref)

      IF (err > tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,ES12.5)') TRIM(label), "failed with error", err
      ELSE
         WRITE (*, '(A,1X,A,1X,ES12.5)') TRIM(label), "OK with error", err
      END IF

   END SUBROUTINE check_offdiag

! **************************************************************************************************
!> \brief Check L_exp(diag(-i*e,-i*e), E_11*z)
! **************************************************************************************************
   SUBROUTINE check_diag(label, e1, z, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: e1
      COMPLEX(KIND=dp), INTENT(IN)                       :: z
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      REAL(KIND=dp), PARAMETER                           :: eps = 1.0E-7_dp

      COMPLEX(KIND=dp)                                   :: fd, ref
      COMPLEX(KIND=dp), DIMENSION(2, 2)                  :: m0, meps, u0, ueps
      REAL(KIND=dp)                                      :: err

      m0(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      m0(1, 1) = CMPLX(0.0_dp, -e1, KIND=dp)
      m0(2, 2) = CMPLX(0.0_dp, -e1, KIND=dp)

      meps = m0
      meps(1, 1) = meps(1, 1) + eps*z

      CALL dense_expm2(m0, u0)
      CALL dense_expm2(meps, ueps)

      fd = (ueps(1, 1) - u0(1, 1))/eps
      ref = qs_ot_complex_exp_frechet_kernel(e1, e1)*z
      err = ABS(fd - ref)

      IF (err > tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,ES12.5)') TRIM(label), "failed with error", err
      ELSE
         WRITE (*, '(A,1X,A,1X,ES12.5)') TRIM(label), "OK with error", err
      END IF

   END SUBROUTINE check_diag

! **************************************************************************************************
!> \brief Small dense 2x2 matrix exponential via Taylor series.
! **************************************************************************************************
   SUBROUTINE dense_expm2(a, expa)
      COMPLEX(KIND=dp), DIMENSION(2, 2), INTENT(IN)      :: a
      COMPLEX(KIND=dp), DIMENSION(2, 2), INTENT(OUT)     :: expa

      COMPLEX(KIND=dp), DIMENSION(2, 2)                  :: term
      INTEGER                                            :: i, n

      expa(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      term(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      DO i = 1, 2
         expa(i, i) = CMPLX(1.0_dp, 0.0_dp, KIND=dp)
         term(i, i) = CMPLX(1.0_dp, 0.0_dp, KIND=dp)
      END DO

      DO n = 1, 80
         term = MATMUL(term, a)/REAL(n, KIND=dp)
         expa = expa + term
         IF (MAXVAL(ABS(term)) < 1.0E-20_dp) EXIT
      END DO

   END SUBROUTINE dense_expm2

END PROGRAM qs_ot_frechet_unittest
