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
   USE bessel_lib,                      ONLY: bessel0
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: fourpi
   USE orbital_pointers,                ONLY: indco,&
                                              init_orbital_pointers,&
                                              nco,&
                                              ncoset
   USE orbital_transformation_matrices, ONLY: init_spherical_harmonics,&
                                              orbtramat
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: upf_inverse_coulomb_radial, upf_inverse_coulomb_radial_response
   PUBLIC :: upf_inverse_coulomb_compensate, upf_inverse_coulomb_compensation_adjoint
   PUBLIC :: upf_inverse_coulomb_fourier_kernel

CONTAINS

! **************************************************************************************************
!> \brief Fourier map of the same radial shells, for coupling to a periodic background.
!>        rho_G = sum_ir,lm weights(ir)*kernel(G,ir,lm)*density(ir,lm).
!>        The reciprocal covector must use the caller's actual full/half-grid multiplicities.
!>        No reciprocal Coulomb cutoff or G=0 electrostatic convention is selected here.
!> \param radii Positive radial shell nodes.
!> \param angular_momenta l of each orthonormal real harmonic channel.
!> \param harmonic_indices Native harmonic indices, between l*l+1 and (l+1)**2.
!> \param gvecs Cartesian reciprocal vectors, including zero if required by the caller.
!> \param position Cartesian center of these radial shells.
!> \param volume Physical cell volume.
!> \param kernel Complex Fourier coefficients per unit radial measure and density.
!> \param response Optional joint geometry response of the same Fourier map.
!> \param radii_response Required with response.
!> \param gvecs_response Required with response, including reciprocal-vector strain.
!> \param position_response Required with response.
!> \param volume_response Required with response.
! **************************************************************************************************
   SUBROUTINE upf_inverse_coulomb_fourier_kernel(radii, angular_momenta, harmonic_indices, gvecs, position, volume, &
                                                 kernel, response, radii_response, gvecs_response, &
                                                 position_response, volume_response)
      REAL(dp), INTENT(IN)                               :: radii(:)
      INTEGER, INTENT(IN)                                :: angular_momenta(:), harmonic_indices(:)
      REAL(dp), INTENT(IN)                               :: gvecs(:, :), position(3), volume
      COMPLEX(dp), INTENT(OUT)                           :: kernel(:, :, :)
      COMPLEX(dp), INTENT(OUT), OPTIONAL                 :: response(:, :, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: radii_response(:), gvecs_response(:, :), &
                                                            position_response(3), volume_response

      COMPLEX(dp)                                        :: phase, phase_l, phase_response
      INTEGER                                            :: a, axis, channel, ig, ir, l, m, &
                                                            powers(3), reduced(3)
      REAL(dp)                                           :: angular, angular_gradient(3), &
                                                            direction(3), factor, gradient(3), jl, &
                                                            jnext, norm_g, radial

      CPASSERT(SIZE(harmonic_indices) == SIZE(angular_momenta))
      CPASSERT(SIZE(angular_momenta) > 0)
      CPASSERT(SIZE(gvecs, 1) == 3)
      IF (.NOT. ALL(SHAPE(kernel) == [SIZE(gvecs, 2), SIZE(radii), SIZE(angular_momenta)])) THEN
         CPABORT("Invalid radial Fourier kernel dimensions")
      END IF
      CPASSERT(ALL(ieee_is_finite(radii)) .AND. ALL(radii > 0.0_dp))
      CPASSERT(ALL(ieee_is_finite(gvecs)))
      CPASSERT(ALL(ieee_is_finite(position)))
      CPASSERT(ieee_is_finite(volume) .AND. volume > 0.0_dp)
      CPASSERT(ALL(angular_momenta >= 0 .AND. angular_momenta <= 20))
      CPASSERT(ALL(harmonic_indices > angular_momenta**2))
      CPASSERT(ALL(harmonic_indices <= (angular_momenta + 1)**2))
      IF (PRESENT(response)) THEN
         CPASSERT(ALL(SHAPE(response) == SHAPE(kernel)))
         CPASSERT(PRESENT(radii_response) .AND. PRESENT(gvecs_response))
         CPASSERT(PRESENT(position_response) .AND. PRESENT(volume_response))
         CPASSERT(SIZE(radii_response) == SIZE(radii))
         CPASSERT(ALL(SHAPE(gvecs_response) == SHAPE(gvecs)))
      END IF
      CALL init_orbital_pointers(MAXVAL(angular_momenta))
      CALL init_spherical_harmonics(MAXVAL(angular_momenta), -1)
      DO ig = 1, SIZE(gvecs, 2)
         norm_g = SQRT(SUM(gvecs(:, ig)**2))
         direction = 0.0_dp
         IF (norm_g > 0.0_dp) direction = gvecs(:, ig)/norm_g
         phase = fourpi/volume*EXP(CMPLX(0.0_dp, -DOT_PRODUCT(gvecs(:, ig), position), dp))
         IF (PRESENT(response)) THEN
            phase_response = CMPLX(-volume_response/volume, &
                                   -DOT_PRODUCT(gvecs_response(:, ig), position) - &
                                   DOT_PRODUCT(gvecs(:, ig), position_response), dp)
         END IF
         DO channel = 1, SIZE(angular_momenta)
            l = angular_momenta(channel)
            m = harmonic_indices(channel) - l*l
            phase_l = phase*CMPLX(0.0_dp, -1.0_dp, dp)**l
            angular = 0.0_dp
            angular_gradient = 0.0_dp
            DO a = 1, nco(l)
               powers = indco(:, ncoset(l - 1) + a)
               factor = orbtramat(l)%slm(m, a)
               angular = angular + factor*PRODUCT(direction**powers)
               IF (.NOT. PRESENT(response)) CYCLE
               DO axis = 1, 3
                  IF (powers(axis) == 0) CYCLE
                  reduced = powers
                  reduced(axis) = reduced(axis) - 1
                  angular_gradient(axis) = angular_gradient(axis) + &
                                           factor*REAL(powers(axis), dp)*PRODUCT(direction**reduced)
               END DO
            END DO
            DO ir = 1, SIZE(radii)
               jl = bessel0(norm_g*radii(ir), l)
               kernel(ig, ir, channel) = phase_l*angular*jl
               IF (.NOT. PRESENT(response)) CYCLE
               gradient = 0.0_dp
               radial = 0.0_dp
               IF (norm_g > 0.0_dp) THEN
                  jnext = bessel0(norm_g*radii(ir), l + 1)
                  gradient = jl/norm_g*angular_gradient - radii(ir)*jnext*angular*direction
                  radial = angular*(l/radii(ir)*jl - norm_g*jnext)
               ELSE IF (l == 1) THEN
                  ! The origin limit contains dipole derivatives even though rho_G itself is zero.
                  gradient = radii(ir)/3.0_dp*angular_gradient
               END IF
               response(ig, ir, channel) = kernel(ig, ir, channel)*phase_response + &
                                           phase_l*(DOT_PRODUCT(gradient, gvecs_response(:, ig)) + &
                                                    radial*radii_response(ir))
            END DO
         END DO
      END DO
   END SUBROUTINE upf_inverse_coulomb_fourier_kernel

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
