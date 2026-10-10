!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculates Bessel functions
!> \note
!>      Functions adapted from netlib
!> \par History
!>      March-2006: Bessel Transform (JGH)
!> \author JGH (10-02-2001)
! **************************************************************************************************
MODULE bessel_lib

   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: fac,&
                                              pio2
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: bessk0, bessk1, bessk0_response, bessk1_response, bessel0

   REAL(dp), PARAMETER :: k0p(7) = [-0.57721566_dp, 0.42278420_dp, 0.23069756_dp, &
                                    0.3488590e-1_dp, 0.262698e-2_dp, 0.10750e-3_dp, 0.74e-5_dp]
   REAL(dp), PARAMETER :: k0q(7) = [1.25331414_dp, -0.7832358e-1_dp, 0.2189568e-1_dp, &
                                    -0.1062446e-1_dp, 0.587872e-2_dp, -0.251540e-2_dp, 0.53208e-3_dp]
   REAL(dp), PARAMETER :: k1p(7) = [1.0_dp, 0.15443144_dp, -0.67278579_dp, -0.18156897_dp, &
                                    -0.1919402e-1_dp, -0.110404e-2_dp, -0.4686e-4_dp]
   REAL(dp), PARAMETER :: k1q(7) = [1.25331414_dp, 0.23498619_dp, -0.3655620e-1_dp, &
                                    0.1504268e-1_dp, -0.780353e-2_dp, 0.325614e-2_dp, -0.68245e-3_dp]
   REAL(dp), PARAMETER :: i0p(7) = [1.0_dp, 3.5156229_dp, 3.0899424_dp, 1.2067492_dp, &
                                    0.2659732_dp, 0.360768e-1_dp, 0.45813e-2_dp]
   REAL(dp), PARAMETER :: i1p(7) = [0.5_dp, 0.87890594_dp, 0.51498869_dp, 0.15084934_dp, &
                                    0.2658733e-1_dp, 0.301532e-2_dp, 0.32411e-3_dp]

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param x must be positive
!> \return ...
! **************************************************************************************************
   ELEMENTAL FUNCTION bessk0(x)
      REAL(dp), INTENT(IN)                               :: x
      REAL(dp)                                           :: bessk0

      REAL(dp)                                           :: y

      IF (x < 2.0_dp) THEN
         y = x*x/4.0_dp
         bessk0 = -LOG(x/2.0_dp)*bessi0(x) + polynomial_value(k0p, y)
      ELSE
         y = 2.0_dp/x
         bessk0 = EXP(-x)/SQRT(x)*polynomial_value(k0q, y)
      END IF
   END FUNCTION bessk0

! **************************************************************************************************
!> \brief ...
!> \param x must be positive
!> \return ...
! **************************************************************************************************
   ELEMENTAL FUNCTION bessk1(x)
      REAL(dp), INTENT(IN)                               :: x
      REAL(dp)                                           :: bessk1

      REAL(dp)                                           :: y

      IF (x < 2.0_dp) THEN
         y = x*x/4.0_dp
         bessk1 = LOG(x/2.0_dp)*bessi1(x) + polynomial_value(k1p, y)/x
      ELSE
         y = 2.0_dp/x
         bessk1 = EXP(-x)/SQRT(x)*polynomial_value(k1q, y)
      END IF
   END FUNCTION bessk1

! **************************************************************************************************
!> \brief ...
!> \param x ...
!> \return ...
! **************************************************************************************************
   ELEMENTAL FUNCTION bessi0(x)

      REAL(KIND=dp), INTENT(IN)                          :: x
      REAL(KIND=dp)                                      :: bessi0

      REAL(KIND=dp), PARAMETER :: q1 = 0.39894228_dp, q2 = 0.1328592e-1_dp, q3 = 0.225319e-2_dp, &
         q4 = -0.157565e-2_dp, q5 = 0.916281e-2_dp, q6 = -0.2057706e-1_dp, q7 = 0.2635537e-1_dp, &
         q8 = -0.1647633e-1_dp, q9 = 0.392377e-2_dp

      REAL(KIND=dp)                                      :: ax, y

      IF (ABS(x) < 3.75_dp) THEN
         y = (x/3.75_dp)**2
         bessi0 = polynomial_value(i0p, y)
      ELSE
         ax = ABS(x)
         y = 3.75_dp/ax
         bessi0 = (EXP(ax)/SQRT(ax))*(q1 + y*(q2 + y* &
                                              (q3 + y*(q4 + y*(q5 + y*(q6 + y*(q7 + y* &
                                                                               (q8 + y*q9))))))))
      END IF

   END FUNCTION bessi0

! **************************************************************************************************
!> \brief ...
!> \param x ...
!> \return ...
! **************************************************************************************************
   ELEMENTAL FUNCTION bessi1(x)

      REAL(KIND=dp), INTENT(IN)                          :: x
      REAL(KIND=dp)                                      :: bessi1

      REAL(KIND=dp), PARAMETER :: q1 = 0.39894228_dp, q2 = -0.3988024e-1_dp, q3 = -0.362018e-2_dp, &
         q4 = 0.163801e-2_dp, q5 = -0.1031555e-1_dp, q6 = 0.2282967e-1_dp, q7 = -0.2895312e-1_dp, &
         q8 = 0.1787654e-1_dp, q9 = -0.420059e-2_dp

      REAL(KIND=dp)                                      :: ax, y

      IF (ABS(x) < 3.75_dp) THEN
         y = (x/3.75_dp)**2
         bessi1 = x*polynomial_value(i1p, y)
      ELSE
         ax = ABS(x)
         y = 3.75_dp/ax
         bessi1 = (EXP(ax)/SQRT(ax))*(q1 + y*(q2 + y* &
                                              (q3 + y*(q4 + y*(q5 + y*(q6 + y*(q7 + y* &
                                                                               (q8 + y*q9))))))))
         IF (x < 0.0_dp) bessi1 = -bessi1
      END IF

   END FUNCTION bessi1

! **************************************************************************************************
!> \brief Derivative of the same rational-polynomial approximation used by bessk0.
!> \param x Positive argument, away from the approximation's branch junction.
!> \return Derivative with respect to x.
! **************************************************************************************************
   ELEMENTAL FUNCTION bessk0_response(x) RESULT(response)
      REAL(dp), INTENT(IN)                               :: x
      REAL(dp)                                           :: response

      response = modified_bessel_response(x, k0p, k0q, i0p, 0)
   END FUNCTION bessk0_response

! **************************************************************************************************
!> \brief Derivative of the same rational-polynomial approximation used by bessk1.
!> \param x Positive argument, away from the approximation's branch junction.
!> \return Derivative with respect to x.
! **************************************************************************************************
   ELEMENTAL FUNCTION bessk1_response(x) RESULT(response)
      REAL(dp), INTENT(IN)                               :: x
      REAL(dp)                                           :: response

      response = modified_bessel_response(x, k1p, k1q, i1p, 1)
   END FUNCTION bessk1_response

! **************************************************************************************************
!> \brief Differentiate the retained small-argument and asymptotic approximations analytically.
!> \param x ...
!> \param p_coeff ...
!> \param q_coeff ...
!> \param i_coeff ...
!> \param order Zero or one.
!> \return ...
! **************************************************************************************************
   PURE FUNCTION modified_bessel_response(x, p_coeff, q_coeff, i_coeff, order) RESULT(response)
      REAL(dp), INTENT(IN)                               :: x, p_coeff(:), q_coeff(:), i_coeff(:)
      INTEGER, INTENT(IN)                                :: order
      REAL(dp)                                           :: response

      REAL(dp)                                           :: i_response, i_value, y

      IF (x < 2.0_dp) THEN
         y = (x/3.75_dp)**2
         i_value = polynomial_value(i_coeff, y)
         i_response = 2.0_dp*x/3.75_dp**2*polynomial_response(i_coeff, y)
         IF (order == 0) THEN
            response = -i_value/x - LOG(x/2.0_dp)*i_response + &
                       0.5_dp*x*polynomial_response(p_coeff, x*x/4.0_dp)
         ELSE
            i_response = i_value + x*i_response
            i_value = x*i_value
            response = i_value/x + LOG(x/2.0_dp)*i_response - &
                       polynomial_value(p_coeff, x*x/4.0_dp)/x**2 + &
                       0.5_dp*polynomial_response(p_coeff, x*x/4.0_dp)
         END IF
      ELSE
         y = 2.0_dp/x
         response = EXP(-x)/SQRT(x)*((-1.0_dp - 0.5_dp/x)*polynomial_value(q_coeff, y) - &
                                     y/x*polynomial_response(q_coeff, y))
      END IF
   END FUNCTION modified_bessel_response

! **************************************************************************************************
!> \brief Evaluate a polynomial using the original nested multiplication order.
!> \param coefficients Ascending powers.
!> \param x ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION polynomial_value(coefficients, x) RESULT(value)
      REAL(dp), INTENT(IN)                               :: coefficients(:), x
      REAL(dp)                                           :: value

      INTEGER                                            :: i

      value = coefficients(SIZE(coefficients))
      DO i = SIZE(coefficients) - 1, 1, -1
         value = coefficients(i) + x*value
      END DO
   END FUNCTION polynomial_value

! **************************************************************************************************
!> \brief Differentiate the same polynomial coefficients with nested multiplication.
!> \param coefficients Ascending powers.
!> \param x ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION polynomial_response(coefficients, x) RESULT(response)
      REAL(dp), INTENT(IN)                               :: coefficients(:), x
      REAL(dp)                                           :: response

      INTEGER                                            :: i

      response = REAL(SIZE(coefficients) - 1, dp)*coefficients(SIZE(coefficients))
      DO i = SIZE(coefficients) - 1, 2, -1
         response = REAL(i - 1, dp)*coefficients(i) + x*response
      END DO
   END FUNCTION polynomial_response

! **************************************************************************************************
!> \brief ...
!> \param x ...
!> \param l ...
!> \return ...
! **************************************************************************************************
   ELEMENTAL IMPURE FUNCTION bessel0(x, l)
      !
      ! Calculates spherical Bessel functions
      ! Abramowitz & Stegun using Formulas 10.1.2, 10.1.8, 10.1.9
      ! Adapted from P. Bloechl
      !
      REAL(KIND=dp), INTENT(IN)                          :: x
      INTEGER, INTENT(IN)                                :: l
      REAL(KIND=dp)                                      :: bessel0

      REAL(KIND=dp), PARAMETER                           :: tol = 1.e-12_dp

      INTEGER                                            :: i, ii, il, isvar, k
      REAL(KIND=dp)                                      :: arg, fact, xsq
      REAL(KIND=dp), DIMENSION(4)                        :: trig

      IF ((x*x - 0.5_dp) > l) THEN
         arg = x - pio2*REAL(l, KIND=dp)
         trig(1) = SIN(arg)/x
         trig(2) = COS(arg)/x
         trig(3) = -trig(1)
         trig(4) = -trig(2)
         bessel0 = trig(1)
         IF (l /= 0) THEN
            xsq = 0.5_dp/x
            fact = 1._dp
            DO k = 1, l
               ii = MOD(k, 4) + 1
               fact = fac(k + l)/fac(k)/fac(l - k)*xsq**k
               bessel0 = bessel0 + fact*trig(ii)
            END DO
         END IF
      ELSE
         ! Taylor expansion for small arguments
         isvar = 1
         DO il = 1, l
            isvar = isvar*(2*il + 1)
         END DO
         IF (l /= 0._dp) THEN
            fact = x**l/REAL(isvar, KIND=dp)
         ELSE
            fact = 1._dp/REAL(isvar, KIND=dp)
         END IF
         bessel0 = fact
         xsq = -0.5_dp*x*x
         isvar = 2*l + 1
         DO i = 1, 1000
            isvar = isvar + 2
            fact = fact*xsq/REAL(i*isvar, KIND=dp)
            bessel0 = bessel0 + fact
            IF (ABS(fact) < tol) EXIT
         END DO
         IF (ABS(fact) > tol) CPABORT("BESSEL0 NOT CONVERGED")
      END IF

   END FUNCTION bessel0

END MODULE bessel_lib

