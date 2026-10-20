!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Common Gaussian representation of nonlinear core densities.
! **************************************************************************************************
MODULE qs_nlcc_methods
   USE external_potential_types,        ONLY: get_potential,&
                                              gth_potential_type,&
                                              sgp_potential_type
   USE kinds,                           ONLY: dp
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: get_nlcc_gaussians
   PUBLIC :: nlcc_gaussian_fields

CONTAINS

! **************************************************************************************************
!> \brief Return rho(r) = sum_i exp(-r**2/(2*alpha(i)**2)) *
!>        sum_j cval(j,i)*(r/alpha(i))**(2*j-2).
!>        The same representation is used for collocation and force integration.
!> \param qs_kind Atomic kind carrying a GTH or separable Gaussian potential.
!> \param alpha Gaussian widths; an empty array indicates no core correction.
!> \param nct Number of radial polynomial coefficients for each Gaussian.
!> \param cval Radial polynomial coefficients, including the density normalization.
! **************************************************************************************************
   SUBROUTINE get_nlcc_gaussians(qs_kind, alpha, nct, cval)
      TYPE(qs_kind_type), INTENT(IN)                     :: qs_kind
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(OUT)                                     :: alpha
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: nct
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: cval

      INTEGER                                            :: nexp
      INTEGER, DIMENSION(:), POINTER                     :: gth_nct
      LOGICAL                                            :: nlcc
      REAL(KIND=dp), DIMENSION(:), POINTER               :: a_nlcc, c_nlcc, gth_alpha
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: gth_cval
      TYPE(gth_potential_type), POINTER                  :: gth_potential
      TYPE(sgp_potential_type), POINTER                  :: sgp_potential

      CALL get_qs_kind(qs_kind, gth_potential=gth_potential, sgp_potential=sgp_potential)
      IF (ASSOCIATED(gth_potential)) THEN
         CALL get_potential(gth_potential, nlcc_present=nlcc, nexp_nlcc=nexp, &
                            alpha_nlcc=gth_alpha, nct_nlcc=gth_nct, cval_nlcc=gth_cval)
         IF (nlcc) THEN
            ALLOCATE (alpha(nexp), nct(nexp), cval(SIZE(gth_cval, 1), nexp))
            alpha(:) = gth_alpha(1:nexp)
            nct(:) = gth_nct(1:nexp)
            cval(:, :) = gth_cval(:, 1:nexp)
            RETURN
         END IF
      ELSE IF (ASSOCIATED(sgp_potential)) THEN
         CALL get_potential(sgp_potential, has_nlcc=nlcc, n_nlcc=nexp, &
                            a_nlcc=a_nlcc, c_nlcc=c_nlcc)
         IF (nlcc) THEN
            CPASSERT(ALL(a_nlcc(1:nexp) > 0.0_dp))
            ! SGP stores exponents in exp(-a*r**2), whereas GTH stores widths.
            ALLOCATE (alpha(nexp), nct(nexp), cval(1, nexp))
            alpha(:) = 1.0_dp/SQRT(2.0_dp*a_nlcc(1:nexp))
            nct = 1
            cval(1, :) = c_nlcc(1:nexp)
            RETURN
         END IF
      END IF
      ALLOCATE (alpha(0), nct(0), cval(0, 0))

   END SUBROUTINE get_nlcc_gaussians

! **************************************************************************************************
!> \brief Native nonlinear core Gaussians and their Cartesian derivatives at explicit points.
!> \param alpha Original Gaussian widths from get_nlcc_gaussians.
!> \param nct Number of original radial polynomial coefficients per Gaussian.
!> \param cval Original normalized coefficients of powers (r/alpha)**2.
!> \param displacement Point-minus-core-center vectors, indexed (axis,point).
!> \param density Overwritten positive core-density convention, without spin splitting.
!> \param gradient Overwritten Cartesian point gradient, indexed (point,axis).
!> \param hessian Optional Cartesian point Hessian, indexed (point,axis,axis).
! **************************************************************************************************
   SUBROUTINE nlcc_gaussian_fields(alpha, nct, cval, displacement, density, gradient, hessian)
      REAL(dp), INTENT(IN)                               :: alpha(:)
      INTEGER, INTENT(IN)                                :: nct(:)
      REAL(dp), INTENT(IN)                               :: cval(:, :), displacement(:, :)
      REAL(dp), INTENT(OUT)                              :: density(:), gradient(:, :)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: hessian(:, :, :)

      INTEGER                                            :: axis, exponent, m, point
      REAL(dp)                                           :: derivative, e, first, inverse_width2, &
                                                            polynomial, second, x

      CPASSERT(SIZE(nct) == SIZE(alpha))
      CPASSERT(SIZE(cval, 2) == SIZE(alpha))
      CPASSERT(SIZE(displacement, 1) == 3)
      CPASSERT(SIZE(displacement, 2) == SIZE(density))
      CPASSERT(ALL(SHAPE(gradient) == [SIZE(density), 3]))
      density = 0.0_dp
      gradient = 0.0_dp
      IF (PRESENT(hessian)) THEN
         CPASSERT(SIZE(hessian, 1) == SIZE(density))
         CPASSERT(SIZE(hessian, 2) == 3 .AND. SIZE(hessian, 3) == 3)
         hessian = 0.0_dp
      END IF
      DO exponent = 1, SIZE(alpha)
         CPASSERT(alpha(exponent) > 0.0_dp)
         CPASSERT(nct(exponent) >= 0)
         CPASSERT(nct(exponent) <= SIZE(cval, 1))
         IF (nct(exponent) == 0) CYCLE
         inverse_width2 = 1.0_dp/alpha(exponent)**2
         DO point = 1, SIZE(density)
            x = SUM(displacement(:, point)**2)*inverse_width2
            polynomial = cval(nct(exponent), exponent)
            derivative = 0.0_dp
            second = 0.0_dp
            DO m = nct(exponent) - 1, 1, -1
               second = second*x + 2.0_dp*derivative
               derivative = derivative*x + polynomial
               polynomial = polynomial*x + cval(m, exponent)
            END DO
            e = EXP(-0.5_dp*x)
            density(point) = density(point) + e*polynomial
            first = 2.0_dp*e*inverse_width2*(derivative - 0.5_dp*polynomial)
            gradient(point, :) = gradient(point, :) + first*displacement(:, point)
            IF (PRESENT(hessian)) THEN
               second = 4.0_dp*e*inverse_width2**2*(second - derivative + 0.25_dp*polynomial)
               DO axis = 1, 3
                  hessian(point, :, axis) = hessian(point, :, axis) + second*displacement(:, point)*displacement(axis, point)
                  hessian(point, axis, axis) = hessian(point, axis, axis) + first
               END DO
            END IF
         END DO
      END DO
   END SUBROUTINE nlcc_gaussian_fields

END MODULE qs_nlcc_methods
