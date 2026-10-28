!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Original Gaussian GTH PPL fields and Fourier coefficients.
!>        Match the polynomial terms passed to the native core_ppl integrator, including extended
!>        local Gaussian terms. The ionic Gaussian charge remains a separate native Poisson source.
! **************************************************************************************************
MODULE gth_local_fields
   USE external_potential_types,        ONLY: gth_potential_type
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: pi
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: gth_local_fields_type, gth_local_fields_init, gth_local_fields_evaluate, gth_local_fourier

   TYPE gth_local_fields_type
      INTEGER, ALLOCATABLE :: nterm(:)
      REAL(dp), ALLOCATABLE :: alpha(:), coefficients(:, :)
   END TYPE gth_local_fields_type

CONTAINS

! **************************************************************************************************
!> \brief Own the same exponents and already scaled coefficients as the native PPL integrator.
!> \param potential Original initialized GTH dataset, without fitting or additional truncation.
!> \param fields Owned Gaussian polynomial terms.
! **************************************************************************************************
   SUBROUTINE gth_local_fields_init(potential, fields)
      TYPE(gth_potential_type), INTENT(IN)               :: potential
      TYPE(gth_local_fields_type), INTENT(OUT)           :: fields

      INTEGER                                            :: i, maximum, n, nextra

      nextra = 0
      IF (potential%lpotextended) nextra = potential%nexp_lpot
      CPASSERT(nextra >= 0)
      maximum = MAX(1, potential%nexp_ppl)
      IF (nextra > 0) THEN
         CPASSERT(ASSOCIATED(potential%alpha_lpot))
         CPASSERT(ASSOCIATED(potential%nct_lpot))
         CPASSERT(ASSOCIATED(potential%cval_lpot))
         CPASSERT(SIZE(potential%alpha_lpot) >= nextra)
         CPASSERT(SIZE(potential%nct_lpot) >= nextra)
         CPASSERT(SIZE(potential%cval_lpot, 2) >= nextra)
         maximum = MAX(maximum, MAXVAL(potential%nct_lpot(:nextra)))
      END IF
      ALLOCATE (fields%nterm(1 + nextra), fields%alpha(1 + nextra), fields%coefficients(maximum, 1 + nextra))
      fields%coefficients = 0.0_dp
      fields%nterm(1) = potential%nexp_ppl
      fields%alpha(1) = potential%alpha_ppl
      n = fields%nterm(1)
      CPASSERT(n >= 0)
      IF (n > 0) THEN
         CPASSERT(ASSOCIATED(potential%cexp_ppl))
         CPASSERT(SIZE(potential%cexp_ppl) == n)
         fields%coefficients(:n, 1) = potential%cexp_ppl
      END IF
      DO i = 1, nextra
         n = potential%nct_lpot(i)
         CPASSERT(n >= 0)
         CPASSERT(SIZE(potential%cval_lpot, 1) >= n)
         fields%nterm(i + 1) = n
         fields%alpha(i + 1) = potential%alpha_lpot(i)
         fields%coefficients(:n, i + 1) = potential%cval_lpot(:n, i)
      END DO
      DO i = 1, SIZE(fields%nterm)
         IF (fields%nterm(i) > 0) THEN
            CPASSERT(fields%alpha(i) > 0.0_dp)
         END IF
      END DO
   END SUBROUTINE gth_local_fields_init

! **************************************************************************************************
!> \brief Evaluate sum_j exp(-alpha_j*r^2) sum_m C_mj r^(2m) and its point gradient.
!> \param fields Original owned local Gaussian terms.
!> \param displacement Point-minus-potential-center vectors, indexed (axis,point).
!> \param value New local potential, without an ionic Coulomb contribution.
!> \param gradient Optional point gradient, indexed (point,axis).
! **************************************************************************************************
   SUBROUTINE gth_local_fields_evaluate(fields, displacement, value, gradient)
      TYPE(gth_local_fields_type), INTENT(IN)            :: fields
      REAL(dp), INTENT(IN)                               :: displacement(:, :)
      REAL(dp), INTENT(OUT)                              :: value(:)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: gradient(:, :)

      INTEGER                                            :: exponent, m, point
      REAL(dp)                                           :: derivative, e, polynomial, r2

      CPASSERT(SIZE(displacement, 1) == 3)
      CPASSERT(SIZE(value) == SIZE(displacement, 2))
      IF (PRESENT(gradient)) THEN
         CPASSERT(SIZE(gradient, 1) == SIZE(value))
         CPASSERT(SIZE(gradient, 2) == 3)
         gradient = 0.0_dp
      END IF
      value = 0.0_dp
      DO exponent = 1, SIZE(fields%nterm)
         IF (fields%nterm(exponent) == 0) CYCLE
         DO point = 1, SIZE(value)
            r2 = SUM(displacement(:, point)**2)
            polynomial = fields%coefficients(fields%nterm(exponent), exponent)
            derivative = 0.0_dp
            DO m = fields%nterm(exponent) - 1, 1, -1
               derivative = derivative*r2 + polynomial
               polynomial = polynomial*r2 + fields%coefficients(m, exponent)
            END DO
            e = EXP(-fields%alpha(exponent)*r2)
            value(point) = value(point) + polynomial*e
            IF (PRESENT(gradient)) THEN
               gradient(point, :) = gradient(point, :) + &
                                    2.0_dp*e*(derivative - fields%alpha(exponent)*polynomial)*displacement(:, point)
            END IF
         END DO
      END DO
   END SUBROUTINE gth_local_fields_evaluate

! **************************************************************************************************
!> \brief Analytic Fourier transform with exp(-i*g.R)/volume, retaining G=0.
!>        Factorial-scaled Laguerre polynomials transform every original r^(2m) Gaussian term.
!> \param fields Original owned local Gaussian terms.
!> \param g Cartesian reciprocal vectors, indexed (axis,point).
!> \param center Physical potential center.
!> \param volume Positive physical cell volume.
!> \param potential New Fourier coefficients, including the original finite zero component.
!> \param center_derivative Optional Cartesian potential-center derivatives, axis first.
!> \param reciprocal_derivative Optional Cartesian reciprocal-vector derivatives, axis first.
! **************************************************************************************************
   SUBROUTINE gth_local_fourier(fields, g, center, volume, potential, center_derivative, reciprocal_derivative)
      TYPE(gth_local_fields_type), INTENT(IN)            :: fields
      REAL(dp), INTENT(IN)                               :: g(:, :), center(3), volume
      COMPLEX(dp), INTENT(OUT)                           :: potential(:)
      COMPLEX(dp), INTENT(OUT), OPTIONAL                 :: center_derivative(:, :), &
                                                            reciprocal_derivative(:, :)

      COMPLEX(dp)                                        :: phase
      INTEGER                                            :: exponent, m, maximum, point
      REAL(dp)                                           :: alpha, coefficient, derivative, dg(3), &
                                                            prefactor, total, value, x
      REAL(dp), ALLOCATABLE                              :: half(:), three_half(:)

      CPASSERT(SIZE(g, 1) == 3)
      CPASSERT(SIZE(potential) == SIZE(g, 2))
      CPASSERT(volume > 0.0_dp)
      IF (PRESENT(center_derivative)) THEN
         CPASSERT(SIZE(center_derivative, 1) == 3)
         CPASSERT(SIZE(center_derivative, 2) == SIZE(potential))
      END IF
      IF (PRESENT(reciprocal_derivative)) THEN
         CPASSERT(SIZE(reciprocal_derivative, 1) == 3)
         CPASSERT(SIZE(reciprocal_derivative, 2) == SIZE(potential))
      END IF
      maximum = MAX(0, MAXVAL(fields%nterm) - 1)
      ALLOCATE (half(0:maximum), three_half(0:maximum))
      DO point = 1, SIZE(potential)
         total = 0.0_dp
         dg = 0.0_dp
         DO exponent = 1, SIZE(fields%nterm)
            IF (fields%nterm(exponent) == 0) CYCLE
            alpha = fields%alpha(exponent)
            x = SUM(g(:, point)**2)/(4.0_dp*alpha)
            prefactor = (pi/alpha)**1.5_dp*EXP(-x)
            IF (prefactor == 0.0_dp) CYCLE
            CALL scaled_laguerre(0.5_dp, x, half)
            CALL scaled_laguerre(1.5_dp, x, three_half)
            coefficient = 1.0_dp
            DO m = 0, fields%nterm(exponent) - 1
               value = coefficient*fields%coefficients(m + 1, exponent)
               total = total + prefactor*value*half(m)
               derivative = -half(m)
               IF (m > 0) derivative = derivative - m*three_half(m - 1)
               dg = dg + prefactor*value*derivative*g(:, point)/(2.0_dp*alpha)
               coefficient = coefficient/alpha
            END DO
         END DO
         phase = EXP(CMPLX(0.0_dp, -DOT_PRODUCT(g(:, point), center), dp))
         potential(point) = phase*total/volume
         IF (PRESENT(center_derivative)) THEN
            center_derivative(:, point) = CMPLX(0.0_dp, -g(:, point), dp)*potential(point)
         END IF
         IF (PRESENT(reciprocal_derivative)) THEN
            reciprocal_derivative(:, point) = phase*dg/volume + CMPLX(0.0_dp, -center, dp)*potential(point)
         END IF
      END DO
   END SUBROUTINE gth_local_fourier

! **************************************************************************************************
!> \brief m! L_m^a(x), avoiding a factorial in the transform and its derivative.
!> \param a ...
!> \param x ...
!> \param value ...
! **************************************************************************************************
   SUBROUTINE scaled_laguerre(a, x, value)
      REAL(dp), INTENT(IN)                               :: a, x
      REAL(dp), INTENT(OUT)                              :: value(0:)

      INTEGER                                            :: m

      value(0) = 1.0_dp
      IF (UBOUND(value, 1) == 0) RETURN
      value(1) = 1.0_dp + a - x
      DO m = 2, UBOUND(value, 1)
         value(m) = (2*m - 1 + a - x)*value(m - 1) - (m - 1)*(m - 1 + a)*value(m - 2)
      END DO
   END SUBROUTINE scaled_laguerre
END MODULE gth_local_fields
