!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
!> \brief Compensated orbital expectations of signed point-density kernels.
MODULE upf_inverse_precision
   USE compensated_arithmetic,          ONLY: add => twofold_add,&
                                              compensated_real_type => twofold_real_type,&
                                              multiply => twofold_multiply,&
                                              product => twofold_product
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: dp
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: upf_inverse_orbital_expectation, upf_inverse_orbital_projector_density, upf_inverse_projector_operator

   TYPE compensated_complex_type
      TYPE(compensated_real_type) :: re, im
   END TYPE compensated_complex_type

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
!> \brief Real projector occupation or its response, retaining the original orbital products.
!>        Two FP64 components preserve cancellation through subsequent atomic field contractions.
!>        Supplying any variation returns the matching signed occupation response instead.
!> \param fields Projector rows, indexed (projector,AO,k).
!> \param coefficients Original Gaussian orbital coefficients, indexed (AO,MO,k).
!> \param occupations Original physical occupations, indexed (MO,k).
!> \param weights Original normalized k weights, applied once.
!> \param density Real symmetric projector occupation or its signed response.
!> \param coefficients_direction Optional original coefficient variation at fixed geometry.
!> \param occupations_direction Optional signed occupation variation.
!> \param fields_direction Optional moving-projector variation in the same row and k order.
! **************************************************************************************************
   SUBROUTINE upf_inverse_orbital_projector_density(fields, coefficients, occupations, weights, density, &
                                                    coefficients_direction, occupations_direction, fields_direction)
      COMPLEX(dp), INTENT(IN)                            :: fields(:, :, :), coefficients(:, :, :)
      REAL(dp), INTENT(IN)                               :: occupations(:, :), weights(:)
      TYPE(compensated_real_type), INTENT(OUT)           :: density(:, :)
      COMPLEX(dp), INTENT(IN), OPTIONAL :: coefficients_direction(:, :, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: occupations_direction(:, :)
      COMPLEX(dp), INTENT(IN), OPTIONAL                  :: fields_direction(:, :, :)

      INTEGER                                            :: a, b, i, ik, mu, nao, nf, nk, nmo
      LOGICAL                                            :: assertion_valid, response, valid
      REAL(dp)                                           :: df
      TYPE(compensated_complex_type)                     :: term
      TYPE(compensated_complex_type), ALLOCATABLE        :: projected(:), tangent(:)
      TYPE(compensated_real_type)                        :: occupation, pair, pair_direction, &
                                                            response_weight

      nf = SIZE(fields, 1); nao = SIZE(fields, 2); nk = SIZE(fields, 3)
      nmo = SIZE(coefficients, 2)
      CPASSERT(ALL(SHAPE(density) == [nf, nf]))
      valid = ALL(SHAPE(coefficients) == [nao, nmo, nk]) .AND. &
              ALL(SHAPE(occupations) == [nmo, nk]) .AND. SIZE(weights) == nk
      CPASSERT(valid)
      assertion_valid = ALL(occupations >= 0.0_dp) .AND. ALL(weights >= 0.0_dp)
      CPASSERT(assertion_valid)
      assertion_valid = ALL(ieee_is_finite(occupations)) .AND. ALL(ieee_is_finite(weights))
      CPASSERT(assertion_valid)
      response = PRESENT(coefficients_direction) .OR. PRESENT(occupations_direction) .OR. PRESENT(fields_direction)
      IF (PRESENT(coefficients_direction)) THEN
         assertion_valid = ALL(SHAPE(coefficients_direction) == SHAPE(coefficients))
         CPASSERT(assertion_valid)
      END IF
      IF (PRESENT(occupations_direction)) THEN
         assertion_valid = ALL(SHAPE(occupations_direction) == SHAPE(occupations))
         CPASSERT(assertion_valid)
      END IF
      IF (PRESENT(fields_direction)) THEN
         assertion_valid = ALL(SHAPE(fields_direction) == SHAPE(fields))
         CPASSERT(assertion_valid)
      END IF
      ALLOCATE (projected(nf), tangent(nf))
      density = compensated_real_type()
      DO ik = 1, nk
         IF (weights(ik) == 0.0_dp) CYCLE
         DO i = 1, nmo
            df = 0.0_dp
            IF (PRESENT(occupations_direction)) df = occupations_direction(i, ik)
            IF (occupations(i, ik) == 0.0_dp .AND. df == 0.0_dp) CYCLE
            projected = compensated_complex_type()
            tangent = compensated_complex_type()
            DO a = 1, nf
               DO mu = 1, nao
                  term = complex_product(fields(a, mu, ik), coefficients(mu, i, ik))
                  projected(a)%re = add(projected(a)%re, term%re)
                  projected(a)%im = add(projected(a)%im, term%im)
                  IF (PRESENT(coefficients_direction)) THEN
                     term = complex_product(fields(a, mu, ik), coefficients_direction(mu, i, ik))
                     tangent(a)%re = add(tangent(a)%re, term%re)
                     tangent(a)%im = add(tangent(a)%im, term%im)
                  END IF
                  IF (PRESENT(fields_direction)) THEN
                     term = complex_product(fields_direction(a, mu, ik), coefficients(mu, i, ik))
                     tangent(a)%re = add(tangent(a)%re, term%re)
                     tangent(a)%im = add(tangent(a)%im, term%im)
                  END IF
               END DO
            END DO
            CALL PRODUCT(weights(ik), occupations(i, ik), occupation%hi, occupation%lo)
            CALL PRODUCT(weights(ik), df, response_weight%hi, response_weight%lo)
            DO b = 1, nf
               DO a = b, nf
                  pair = real_pair(projected(a), projected(b))
                  IF (response) THEN
                     pair_direction = add(real_pair(tangent(a), projected(b)), real_pair(projected(a), tangent(b)))
                     pair = add(multiply(occupation, pair_direction), multiply(response_weight, pair))
                  ELSE
                     pair = multiply(occupation, pair)
                  END IF
                  density(a, b) = add(density(a, b), pair)
               END DO
            END DO
         END DO
      END DO
      DO b = 1, nf
         DO a = b + 1, nf
            density(b, a) = density(a, b)
         END DO
      END DO
   END SUBROUTINE upf_inverse_orbital_projector_density

! **************************************************************************************************
!> \brief Original complex projector congruence, retaining an atomic kernel until final AO rounding.
!> \param fields Complex projector rows, indexed (projector,AO).
!> \param kernel Real symmetric atomic potential, including its retained product remainders.
!> \param ao_operator Hermitian Gaussian operator W^dagger kernel W.
! **************************************************************************************************
   SUBROUTINE upf_inverse_projector_operator(fields, kernel, ao_operator)
      COMPLEX(dp), INTENT(IN)                            :: fields(:, :)
      TYPE(compensated_real_type), INTENT(IN)            :: kernel(:, :)
      COMPLEX(dp), INTENT(OUT)                           :: ao_operator(:, :)

      INTEGER                                            :: a, b, i, ia, ib, j, n, nf
      INTEGER, ALLOCATABLE                               :: active(:)
      LOGICAL                                            :: assertion_valid
      TYPE(compensated_complex_type)                     :: value
      TYPE(compensated_complex_type), ALLOCATABLE        :: applied(:)
      TYPE(compensated_real_type)                        :: im, re

      nf = SIZE(fields, 1); n = SIZE(fields, 2)
      CPASSERT(ALL(SHAPE(kernel) == [nf, nf]))
      CPASSERT(ALL(SHAPE(ao_operator) == [n, n]))
      assertion_valid = MAXVAL(ABS(kernel%hi - TRANSPOSE(kernel%hi))) <= &
         128.0_dp*EPSILON(0.0_dp)*MAX(MAXVAL(ABS(kernel%hi)), TINY(0.0_dp))
      CPASSERT(assertion_valid)
      ALLOCATE (applied(nf))
      ALLOCATE (active(COUNT(ANY(fields /= CMPLX(0.0_dp, 0.0_dp, dp), DIM=1))))
      active(:) = PACK([(i, i=1, n)], ANY(fields /= CMPLX(0.0_dp, 0.0_dp, dp), DIM=1))
      ao_operator = CMPLX(0.0_dp, 0.0_dp, dp)
      DO ib = 1, SIZE(active)
         b = active(ib)
         applied = compensated_complex_type()
         DO j = 1, nf
            IF (fields(j, b) == CMPLX(0.0_dp, 0.0_dp, dp)) CYCLE
            re%hi = REAL(fields(j, b), dp); re%lo = 0.0_dp
            im%hi = AIMAG(fields(j, b)); im%lo = 0.0_dp
            DO i = 1, nf
               applied(i)%re = add(applied(i)%re, multiply(kernel(i, j), re))
               applied(i)%im = add(applied(i)%im, multiply(kernel(i, j), im))
            END DO
         END DO
         DO ia = 1, ib
            a = active(ia)
            value = compensated_complex_type()
            DO i = 1, nf
               re%hi = REAL(fields(i, a), dp); re%lo = 0.0_dp
               im%hi = AIMAG(fields(i, a)); im%lo = 0.0_dp
               value%re = add(value%re, add(multiply(re, applied(i)%re), multiply(im, applied(i)%im)))
               im%hi = -im%hi
               value%im = add(value%im, add(multiply(re, applied(i)%im), multiply(im, applied(i)%re)))
            END DO
            ao_operator(a, b) = CMPLX(value%re%hi + value%re%lo, value%im%hi + value%im%lo, dp)
            IF (a == b) ao_operator(a, b) = CMPLX(REAL(ao_operator(a, b), dp), 0.0_dp, dp)
            ao_operator(b, a) = CONJG(ao_operator(a, b))
         END DO
      END DO
   END SUBROUTINE upf_inverse_projector_operator

! **************************************************************************************************
!> \brief Real part of a times conjugate(b), retaining both components.
!> \param a First projector amplitude.
!> \param b Second projector amplitude.
!> \return Real Hermitian projector product.
! **************************************************************************************************
   FUNCTION real_pair(a, b) RESULT(value)
      TYPE(compensated_complex_type), INTENT(IN)         :: a, b
      TYPE(compensated_real_type)                        :: value

      value = add(multiply(a%re, b%re), multiply(a%im, b%im))
   END FUNCTION real_pair

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
