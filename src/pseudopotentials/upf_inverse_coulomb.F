!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
!> \brief Radial multipole Coulomb map for signed inverse reconstruction densities.
!>        All maps use the same finite radial measure, including the diagonal shell term.
!>        Periodic background and interactions between distinct centers are assembled separately.
MODULE upf_inverse_coulomb
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: fourpi
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: upf_inverse_coulomb_radial, upf_inverse_coulomb_radial_response
   PUBLIC :: upf_inverse_coulomb_compensate, upf_inverse_coulomb_compensation_adjoint

CONTAINS

! **************************************************************************************************
!> \brief Apply the symmetric Coulomb kernel to real, orthonormal harmonic coefficients.
!> \param radii Positive radial nodes, in arbitrary order.
!> \param weights Radial measure including r**2 dr, without angular weights.
!> \param angular_momenta l of each real harmonic channel. Repeated l values are allowed.
!> \param density Signed radial density coefficients, indexed (radius,harmonic).
!> \param potential Unweighted Coulomb potential in the same harmonic order.
!> \param energy One half of the density/potential contraction with the same radial measure.
! **************************************************************************************************
   SUBROUTINE upf_inverse_coulomb_radial(radii, weights, angular_momenta, density, potential, energy)
      REAL(dp), INTENT(IN)                               :: radii(:), weights(:)
      INTEGER, INTENT(IN)                                :: angular_momenta(:)
      REAL(dp), INTENT(IN)                               :: density(:, :)
      REAL(dp), INTENT(OUT)                              :: potential(:, :), energy

      INTEGER                                            :: channel, i, j, l
      REAL(dp)                                           :: kernel

      CALL check_radial(radii, weights, angular_momenta, density)
      CPASSERT(ALL(SHAPE(potential) == SHAPE(density)))
      potential = 0.0_dp
      DO channel = 1, SIZE(angular_momenta)
         l = angular_momenta(channel)
         DO i = 1, SIZE(radii)
            DO j = 1, SIZE(radii)
               kernel = fourpi/REAL(2*l + 1, dp)/MAX(radii(i), radii(j))* &
                        (MIN(radii(i), radii(j))/MAX(radii(i), radii(j)))**l
               potential(i, channel) = potential(i, channel) + kernel*weights(j)*density(j, channel)
            END DO
         END DO
      END DO
      energy = 0.5_dp*SUM(SPREAD(weights, 2, SIZE(angular_momenta))*density*potential)
   END SUBROUTINE upf_inverse_coulomb_radial

! **************************************************************************************************
!> \brief Differentiate the same finite Coulomb energy, including the measure and kernel.
!>        Equal nodes must have equal radial motion so their common shell remains differentiable.
!> \param radii ...
!> \param weights ...
!> \param angular_momenta ...
!> \param density ...
!> \param radii_response Radial node motion. Zero for Cartesian grids moving with their centers.
!> \param weights_response Complete radial measure response.
!> \param density_response Signed density response, including the coupled inverse coefficients.
!> \param response Coulomb energy response.
! **************************************************************************************************
   SUBROUTINE upf_inverse_coulomb_radial_response(radii, weights, angular_momenta, density, &
                                                  radii_response, weights_response, density_response, response)
      REAL(dp), INTENT(IN)                               :: radii(:), weights(:)
      INTEGER, INTENT(IN)                                :: angular_momenta(:)
      REAL(dp), INTENT(IN)                               :: density(:, :), radii_response(:), &
                                                            weights_response(:), &
                                                            density_response(:, :)
      REAL(dp), INTENT(OUT)                              :: response

      INTEGER                                            :: channel, i, j, l, large, small
      REAL(dp)                                           :: derivative, energy, kernel
      REAL(dp), ALLOCATABLE                              :: potential(:, :)

      CPASSERT(SIZE(radii_response) == SIZE(radii))
      CPASSERT(SIZE(weights_response) == SIZE(weights))
      CPASSERT(ALL(SHAPE(density_response) == SHAPE(density)))
      ALLOCATE (potential(SIZE(density, 1), SIZE(density, 2)))
      CALL upf_inverse_coulomb_radial(radii, weights, angular_momenta, density, potential, energy)
      response = SUM((SPREAD(weights, 2, SIZE(angular_momenta))*density_response + &
                      SPREAD(weights_response, 2, SIZE(angular_momenta))*density)*potential)
      DO channel = 1, SIZE(angular_momenta)
         l = angular_momenta(channel)
         DO i = 1, SIZE(radii)
            DO j = 1, SIZE(radii)
               small = i
               large = j
               IF (radii(i) > radii(j)) THEN
                  small = j
                  large = i
               END IF
               IF (radii(i) == radii(j)) THEN
                  CPASSERT(radii_response(i) == radii_response(j))
               END IF
               kernel = fourpi/REAL(2*l + 1, dp)/radii(large)*(radii(small)/radii(large))**l
               derivative = kernel*(l*radii_response(small)/radii(small) - &
                                    (l + 1)*radii_response(large)/radii(large))
               response = response + 0.5_dp*weights(i)*weights(j)* &
                          density(i, channel)*density(j, channel)*derivative
            END DO
         END DO
      END DO
   END SUBROUTINE upf_inverse_coulomb_radial_response

! **************************************************************************************************
!> \brief Match each multipole with a supplied radial shape in the defining finite quadrature.
!>        The signed remainder has zero exterior moment in every represented harmonic channel.
!>        This identity applies to this radial representation, not to discarded angular channels.
!> \param radii ...
!> \param weights ...
!> \param angular_momenta ...
!> \param density ...
!> \param shapes Positive radial compensation shapes, indexed (radius,harmonic).
!> \param compensation Moment-matched compensation density in the same radial representation.
!> \param residual Signed density minus compensation.
!> \param moments Original exterior multipole moments.
! **************************************************************************************************
   SUBROUTINE upf_inverse_coulomb_compensate(radii, weights, angular_momenta, density, shapes, &
                                             compensation, residual, moments)
      REAL(dp), INTENT(IN)                               :: radii(:), weights(:)
      INTEGER, INTENT(IN)                                :: angular_momenta(:)
      REAL(dp), INTENT(IN)                               :: density(:, :), shapes(:, :)
      REAL(dp), INTENT(OUT)                              :: compensation(:, :), residual(:, :), &
                                                            moments(:)

      INTEGER                                            :: channel
      REAL(dp)                                           :: norm

      CALL check_radial(radii, weights, angular_momenta, density)
      CPASSERT(ALL(SHAPE(shapes) == SHAPE(density)))
      CPASSERT(ALL(SHAPE(compensation) == SHAPE(density)))
      CPASSERT(ALL(SHAPE(residual) == SHAPE(density)))
      CPASSERT(SIZE(moments) == SIZE(angular_momenta))
      CPASSERT(ALL(ieee_is_finite(shapes)) .AND. ALL(shapes >= 0.0_dp))
      DO channel = 1, SIZE(angular_momenta)
         norm = SUM(weights*radii**angular_momenta(channel)*shapes(:, channel))
         CPASSERT(norm > 0.0_dp)
         moments(channel) = SUM(weights*radii**angular_momenta(channel)*density(:, channel))
         compensation(:, channel) = moments(channel)/norm*shapes(:, channel)
      END DO
      residual = density - compensation
   END SUBROUTINE upf_inverse_coulomb_compensate

! **************************************************************************************************
!> \brief Apply the radial-measure adjoint of density minus its moment compensation.
!>        Compensation depends on the occupation density and contributes to the AO operator.
!> \param radii ...
!> \param weights ...
!> \param angular_momenta ...
!> \param shapes Same fixed radial compensation shapes as in the forward map.
!> \param potential Unweighted covector acting on the signed remainder.
!> \param adjoint Unweighted covector acting on the original density.
! **************************************************************************************************
   SUBROUTINE upf_inverse_coulomb_compensation_adjoint(radii, weights, angular_momenta, shapes, potential, adjoint)
      REAL(dp), INTENT(IN)                               :: radii(:), weights(:)
      INTEGER, INTENT(IN)                                :: angular_momenta(:)
      REAL(dp), INTENT(IN)                               :: shapes(:, :), potential(:, :)
      REAL(dp), INTENT(OUT)                              :: adjoint(:, :)

      INTEGER                                            :: channel
      REAL(dp)                                           :: norm

      CALL check_radial(radii, weights, angular_momenta, potential)
      CPASSERT(ALL(SHAPE(shapes) == SHAPE(potential)))
      CPASSERT(ALL(SHAPE(adjoint) == SHAPE(potential)))
      CPASSERT(ALL(ieee_is_finite(shapes)) .AND. ALL(shapes >= 0.0_dp))
      DO channel = 1, SIZE(angular_momenta)
         norm = SUM(weights*radii**angular_momenta(channel)*shapes(:, channel))
         CPASSERT(norm > 0.0_dp)
         adjoint(:, channel) = potential(:, channel) - radii**angular_momenta(channel)* &
                               SUM(weights*shapes(:, channel)*potential(:, channel))/norm
      END DO
   END SUBROUTINE upf_inverse_coulomb_compensation_adjoint

! **************************************************************************************************
!> \brief Validate the common real radial representation.
!> \param radii ...
!> \param weights ...
!> \param angular_momenta ...
!> \param density ...
! **************************************************************************************************
   SUBROUTINE check_radial(radii, weights, angular_momenta, density)
      REAL(dp), INTENT(IN)                               :: radii(:), weights(:)
      INTEGER, INTENT(IN)                                :: angular_momenta(:)
      REAL(dp), INTENT(IN)                               :: density(:, :)

      CPASSERT(SIZE(weights) == SIZE(radii))
      CPASSERT(ALL(SHAPE(density) == [SIZE(radii), SIZE(angular_momenta)]))
      CPASSERT(ALL(ieee_is_finite(radii)) .AND. ALL(radii > 0.0_dp))
      CPASSERT(ALL(ieee_is_finite(weights)) .AND. ALL(weights >= 0.0_dp))
      CPASSERT(ALL(ieee_is_finite(density)))
      CPASSERT(ALL(angular_momenta >= 0))
   END SUBROUTINE check_radial

END MODULE upf_inverse_coulomb
