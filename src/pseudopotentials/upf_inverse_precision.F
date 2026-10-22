!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
!> \brief Compensated orbital expectations of signed point-density kernels.
MODULE upf_inverse_precision
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE iso_c_binding,                   ONLY: c_double
   USE kinds,                           ONLY: dp
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: upf_inverse_orbital_expectation

   TYPE compensated_real_type
      REAL(dp) :: hi = 0.0_dp, lo = 0.0_dp
   END TYPE compensated_real_type

   TYPE compensated_complex_type
      TYPE(compensated_real_type) :: re, im
   END TYPE compensated_complex_type

   INTERFACE
      FUNCTION fused_product_sum(a, b, c) RESULT(value) BIND(C, NAME="fma")
         IMPORT :: c_double
         REAL(c_double), VALUE, INTENT(IN) :: a, b, c
         REAL(c_double) :: value
      END FUNCTION fused_product_sum
   END INTERFACE

CONTAINS

! **************************************************************************************************
!> \brief Evaluate sum(k,i) w(k) f(i,k) (F(k) C(i,k))^dagger Q (F(k) C(i,k)).
!>        Keep the orbitals rather than forming a rounded occupation matrix first. Two FP64
!>        components retain products and cancellation between the original signed terms.
!>        Q is not required to be positive. No eigenvalue or density clipping is applied.
!> \param fields Point-field rows F, indexed (field,AO,k), including any projector transforms.
!> \param kernel Real symmetric point kernel Q in the same field-row order.
!> \param coefficients Original complex Gaussian MO coefficients, indexed (AO,MO,k).
!> \param occupations Original physical MO occupations, indexed (MO,k).
!> \param weights Original k weights, applied once.
!> \param fixed Additive fixed core contribution in the same physical units.
!> \return Original point expectation rounded to FP64 after the complete compensated sum.
! **************************************************************************************************
   FUNCTION upf_inverse_orbital_expectation(fields, kernel, coefficients, occupations, weights, fixed) RESULT(value)
      COMPLEX(dp), INTENT(IN)                            :: fields(:, :, :)
      REAL(dp), INTENT(IN)                               :: kernel(:, :)
      COMPLEX(dp), INTENT(IN)                            :: coefficients(:, :, :)
      REAL(dp), INTENT(IN)                               :: occupations(:, :), weights(:), fixed
      REAL(dp)                                           :: value

      INTEGER                                            :: a, b, i, ik, mu, nao, nf, nk, nmo
      LOGICAL                                            :: valid
      REAL(dp)                                           :: scale
      TYPE(compensated_complex_type)                     :: term
      TYPE(compensated_complex_type), ALLOCATABLE        :: applied(:), projected(:)
      TYPE(compensated_real_type)                        :: accumulator, occupation, point, work

      nf = SIZE(fields, 1)
      nao = SIZE(fields, 2)
      nk = SIZE(fields, 3)
      nmo = SIZE(coefficients, 2)
      CPASSERT(nf > 0 .AND. nao > 0 .AND. nk > 0)
      CPASSERT(ALL(SHAPE(kernel) == [nf, nf]))
      valid = SIZE(coefficients, 1) == nao .AND. SIZE(coefficients, 3) == nk
      CPASSERT(valid)
      valid = ALL(SHAPE(occupations) == [nmo, nk]) .AND. SIZE(weights) == nk
      CPASSERT(valid)
      CPASSERT(ALL(occupations >= 0.0_dp) .AND. ALL(weights >= 0.0_dp))
      valid = ALL(ieee_is_finite(kernel)) .AND. ALL(ieee_is_finite(occupations))
      CPASSERT(valid)
      CPASSERT(ALL(ieee_is_finite(weights)) .AND. ieee_is_finite(fixed))
      valid = ALL(ieee_is_finite(REAL(fields, dp))) .AND. ALL(ieee_is_finite(AIMAG(fields)))
      CPASSERT(valid)
      valid = ALL(ieee_is_finite(REAL(coefficients, dp))) .AND. ALL(ieee_is_finite(AIMAG(coefficients)))
      CPASSERT(valid)
      scale = MAX(MAXVAL(ABS(kernel)), TINY(0.0_dp))
      valid = MAXVAL(ABS(kernel - TRANSPOSE(kernel))) <= 128.0_dp*EPSILON(0.0_dp)*scale
      CPASSERT(valid)
      ALLOCATE (projected(nf), applied(nf))
      accumulator%hi = fixed
      DO ik = 1, nk
         IF (weights(ik) == 0.0_dp) CYCLE
         DO i = 1, nmo
            IF (occupations(i, ik) == 0.0_dp) CYCLE
            projected = compensated_complex_type()
            applied = compensated_complex_type()
            DO a = 1, nf
               DO mu = 1, nao
                  term = complex_product(fields(a, mu, ik), coefficients(mu, i, ik))
                  projected(a)%re = add(projected(a)%re, term%re)
                  projected(a)%im = add(projected(a)%im, term%im)
               END DO
            END DO
            DO a = 1, nf
               DO b = 1, nf
                  IF (kernel(a, b) == 0.0_dp) CYCLE
                  work%hi = kernel(a, b)
                  work%lo = 0.0_dp
                  applied(a)%re = add(applied(a)%re, multiply(work, projected(b)%re))
                  applied(a)%im = add(applied(a)%im, multiply(work, projected(b)%im))
               END DO
            END DO
            point = compensated_real_type()
            DO a = 1, nf
               point = add(point, multiply(projected(a)%re, applied(a)%re))
               point = add(point, multiply(projected(a)%im, applied(a)%im))
            END DO
            CALL PRODUCT(occupations(i, ik), weights(ik), occupation%hi, occupation%lo)
            accumulator = add(accumulator, multiply(occupation, point))
         END DO
      END DO
      value = accumulator%hi + accumulator%lo
      CPASSERT(ieee_is_finite(value))
   END FUNCTION upf_inverse_orbital_expectation

! **************************************************************************************************
!> \brief Error-free sum of two FP64 inputs, protected from reassociation across rounding steps.
!> \param a First input.
!> \param b Second input.
!> \param hi Rounded sum.
!> \param lo Its rounding remainder.
! **************************************************************************************************
   SUBROUTINE two_sum(a, b, hi, lo)
      REAL(dp), INTENT(IN)                               :: a, b
      REAL(dp), INTENT(OUT)                              :: hi, lo

      REAL(dp), VOLATILE                                 :: remainder_a, remainder_b, rounded, &
                                                            virtual_a, virtual_b

      rounded = a + b
      virtual_b = rounded - a
      virtual_a = rounded - virtual_b
      remainder_b = b - virtual_b
      remainder_a = a - virtual_a
      hi = rounded
      lo = remainder_a + remainder_b
   END SUBROUTINE two_sum

! **************************************************************************************************
!> \brief Product and its rounding remainder from the standard C fused multiply-add.
!> \param a First factor.
!> \param b Second factor.
!> \param hi Rounded product.
!> \param lo Its rounding remainder.
! **************************************************************************************************
   SUBROUTINE product(a, b, hi, lo)
      REAL(dp), INTENT(IN)                               :: a, b
      REAL(dp), INTENT(OUT)                              :: hi, lo

      REAL(dp), VOLATILE                                 :: rounded

      rounded = a*b
      hi = rounded
      lo = fused_product_sum(a, b, -rounded)
   END SUBROUTINE product

! **************************************************************************************************
!> \brief Add normalized two-component FP64 numbers.
!> \param a First input.
!> \param b Second input.
!> \return Compensated sum.
! **************************************************************************************************
   FUNCTION add(a, b) RESULT(value)
      TYPE(compensated_real_type), INTENT(IN)            :: a, b
      TYPE(compensated_real_type)                        :: value

      REAL(dp)                                           :: hi, lo

      CALL two_sum(a%hi, b%hi, hi, lo)
      lo = lo + a%lo + b%lo
      CALL two_sum(hi, lo, value%hi, value%lo)
   END FUNCTION add

! **************************************************************************************************
!> \brief Multiply normalized two-component FP64 numbers.
!> \param a First input.
!> \param b Second input.
!> \return Compensated product.
! **************************************************************************************************
   FUNCTION multiply(a, b) RESULT(value)
      TYPE(compensated_real_type), INTENT(IN)            :: a, b
      TYPE(compensated_real_type)                        :: value

      REAL(dp)                                           :: hi, lo

      CALL PRODUCT(a%hi, b%hi, hi, lo)
      lo = lo + a%hi*b%lo + a%lo*b%hi + a%lo*b%lo
      CALL two_sum(hi, lo, value%hi, value%lo)
   END FUNCTION multiply

! **************************************************************************************************
!> \brief Retain the product of two original complex FP64 field coefficients.
!> \param a First factor.
!> \param b Second factor.
!> \return Compensated real and imaginary components.
! **************************************************************************************************
   FUNCTION complex_product(a, b) RESULT(value)
      COMPLEX(dp), INTENT(IN)                            :: a, b
      TYPE(compensated_complex_type)                     :: value

      TYPE(compensated_real_type)                        :: first, second

      CALL PRODUCT(REAL(a, dp), REAL(b, dp), first%hi, first%lo)
      CALL PRODUCT(-AIMAG(a), AIMAG(b), second%hi, second%lo)
      value%re = add(first, second)
      CALL PRODUCT(REAL(a, dp), AIMAG(b), first%hi, first%lo)
      CALL PRODUCT(AIMAG(a), REAL(b, dp), second%hi, second%lo)
      value%im = add(first, second)
   END FUNCTION complex_product

END MODULE upf_inverse_precision
