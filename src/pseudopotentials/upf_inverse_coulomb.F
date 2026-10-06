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
   USE mathconstants,                   ONLY: fourpi,&
                                              pi,&
                                              twopi
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
   PUBLIC :: upf_inverse_coulomb_pair_kernel

CONTAINS

! **************************************************************************************************
!> \brief Reciprocal Coulomb map between two translated harmonic radial-shell distributions.
!>        The kernel excludes both radial measures. Splitting the angular integral at the source
!>        shell avoids sampling a discontinuous radial derivative across an unsplit angular rule.
!>        The joint response differentiates the same finite rule, including its moving boundaries.
!> \param radii_a Positive target radii.
!> \param l_a Angular momenta of the target real harmonics.
!> \param indices_a Native harmonic indices of the target channels.
!> \param radii_b Positive source radii.
!> \param l_b Angular momenta of the source real harmonics.
!> \param indices_b Native harmonic indices of the source channels.
!> \param displacement Target center minus source center, including any lattice image.
!> \param kernel Coulomb map in (radius_a,harmonic_a,radius_b,harmonic_b) order.
!> \param response Optional joint geometry response.
!> \param dr_a Required with response.
!> \param dr_b Required with response.
!> \param ddisplacement Required with response.
! **************************************************************************************************
   SUBROUTINE upf_inverse_coulomb_pair_kernel(radii_a, l_a, indices_a, radii_b, l_b, indices_b, &
                                              displacement, kernel, response, dr_a, dr_b, ddisplacement)
      REAL(dp), INTENT(IN)                               :: radii_a(:)
      INTEGER, INTENT(IN)                                :: l_a(:), indices_a(:)
      REAL(dp), INTENT(IN)                               :: radii_b(:)
      INTEGER, INTENT(IN)                                :: l_b(:), indices_b(:)
      REAL(dp), INTENT(IN)                               :: displacement(3)
      REAL(dp), INTENT(OUT)                              :: kernel(:, :, :, :)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: response(:, :, :, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: dr_a(:), dr_b(:), ddisplacement(3)

      INTEGER                                            :: a, b, i, j
      REAL(dp)                                           :: da, db, delta(3), distance, value
      REAL(dp), ALLOCATABLE                              :: ab(:, :), ba(:, :), dab(:, :), dba(:, :)

      CPASSERT(ALL(radii_a > 0.0_dp) .AND. ALL(ieee_is_finite(radii_a)))
      CPASSERT(ALL(radii_b > 0.0_dp) .AND. ALL(ieee_is_finite(radii_b)))
      CPASSERT(ALL(ieee_is_finite(displacement)))
      CPASSERT(SIZE(l_a) > 0 .AND. SIZE(l_b) > 0)
      CPASSERT(SIZE(indices_a) == SIZE(l_a))
      CPASSERT(SIZE(indices_b) == SIZE(l_b))
      CPASSERT(ALL(l_a >= 0 .AND. l_a <= 20))
      CPASSERT(ALL(l_b >= 0 .AND. l_b <= 20))
      CPASSERT(ALL(indices_a > l_a**2 .AND. indices_a <= (l_a + 1)**2))
      CPASSERT(ALL(indices_b > l_b**2 .AND. indices_b <= (l_b + 1)**2))
      IF (.NOT. ALL(SHAPE(kernel) == [SIZE(radii_a), SIZE(l_a), SIZE(radii_b), SIZE(l_b)])) THEN
         CPABORT("Invalid translated shell kernel dimensions")
      END IF
      IF (PRESENT(response)) THEN
         CPASSERT(PRESENT(dr_a) .AND. PRESENT(dr_b))
         CPASSERT(PRESENT(ddisplacement))
         CPASSERT(SIZE(dr_a) == SIZE(radii_a))
         CPASSERT(SIZE(dr_b) == SIZE(radii_b))
         CPASSERT(ALL(SHAPE(response) == SHAPE(kernel)))
         response = 0.0_dp
      END IF
      kernel = 0.0_dp
      distance = SQRT(SUM(displacement**2))
      delta = 0.0_dp
      IF (PRESENT(response)) delta = ddisplacement
      IF (distance == 0.0_dp) THEN
         ! Same-center shells share their center motion. Distinct coincident centers need a limit.
         CPASSERT(ALL(delta == 0.0_dp))
         DO a = 1, SIZE(l_a)
            DO b = 1, SIZE(l_b)
               IF (indices_a(a) /= indices_b(b)) CYCLE
               DO i = 1, SIZE(radii_a)
                  DO j = 1, SIZE(radii_b)
                     value = fourpi/REAL(2*l_a(a) + 1, dp)/MAX(radii_a(i), radii_b(j))* &
                             (MIN(radii_a(i), radii_b(j))/MAX(radii_a(i), radii_b(j)))**l_a(a)
                     kernel(i, a, j, b) = value
                     IF (.NOT. PRESENT(response)) CYCLE
                     IF (radii_a(i) == radii_b(j)) THEN
                        CPASSERT(dr_a(i) == dr_b(j))
                     END IF
                     IF (radii_a(i) >= radii_b(j)) THEN
                        response(i, a, j, b) = value*(l_a(a)*dr_b(j)/radii_b(j) - &
                                                      (l_a(a) + 1)*dr_a(i)/radii_a(i))
                     ELSE
                        response(i, a, j, b) = value*(l_a(a)*dr_a(i)/radii_a(i) - &
                                                      (l_a(a) + 1)*dr_b(j)/radii_b(j))
                     END IF
                  END DO
               END DO
            END DO
         END DO
         RETURN
      END IF
      CALL init_orbital_pointers(MAX(MAXVAL(l_a), MAXVAL(l_b)))
      CALL init_spherical_harmonics(MAX(MAXVAL(l_a), MAXVAL(l_b)), -1)
      ALLOCATE (ab(SIZE(l_a), SIZE(l_b)), ba(SIZE(l_b), SIZE(l_a)), &
                dab(SIZE(l_a), SIZE(l_b)), dba(SIZE(l_b), SIZE(l_a)))
      DO i = 1, SIZE(radii_a)
         DO j = 1, SIZE(radii_b)
            da = 0.0_dp
            db = 0.0_dp
            IF (PRESENT(response)) THEN
               da = dr_a(i)
               db = dr_b(j)
            END IF
            CALL translated_shell(radii_a(i), radii_b(j), l_a, indices_a, l_b, indices_b, &
                                  displacement, da, db, delta, ab, dab)
            CALL translated_shell(radii_b(j), radii_a(i), l_b, indices_b, l_a, indices_a, &
                                  -displacement, db, da, -delta, ba, dba)
            kernel(i, :, j, :) = 0.5_dp*(ab + TRANSPOSE(ba))
            IF (PRESENT(response)) response(i, :, j, :) = 0.5_dp*(dab + TRANSPOSE(dba))
         END DO
      END DO
   END SUBROUTINE upf_inverse_coulomb_pair_kernel

! **************************************************************************************************
!> \brief Integrate a source shell's regular/irregular solid-harmonic potential on a target shell.
!>        Azimuthal integration is exact for the represented polynomial degree. Polar subintervals
!>        also limit the source-distance ratio to two in the irregular branch.
!> \param a ...
!> \param b ...
!> \param la ...
!> \param ia ...
!> \param lb ...
!> \param ib ...
!> \param d ...
!> \param da ...
!> \param db ...
!> \param dd ...
!> \param value ...
!> \param response ...
! **************************************************************************************************
   SUBROUTINE translated_shell(a, b, la, ia, lb, ib, d, da, db, dd, value, response)
      REAL(dp), INTENT(IN)                               :: a, b
      INTEGER, INTENT(IN)                                :: la(:), ia(:), lb(:), ib(:)
      REAL(dp), INTENT(IN)                               :: d(3), da, db, dd(3)
      REAL(dp), INTENT(OUT)                              :: value(:, :), response(:, :)

      INTEGER                                            :: axis, channel, iphi, iq, nphi, nq
      LOGICAL                                            :: regular
      REAL(dp) :: azimuth, bound, dbound, dex(3), dey(3), dez(3), dfactor, distance, dleft, &
         dlength, dmu, dright, dt, du(3), dweight, dx(3), ex(3), ey(3), ez(3), factor, left, &
         length, mu, norm_x, reference(3), right, scale, sine, t, tmax, u(3), weight, x(3)
      REAL(dp), ALLOCATABLE                              :: dya(:), dyb(:), ga(:, :), gb(:, :), &
                                                            nodes(:), weights(:), ya(:), yb(:)

      distance = SQRT(SUM(d**2))
      dt = DOT_PRODUCT(d, dd)/distance
      ez = d/distance
      dez = (dd - ez*dt)/distance
      reference = 0.0_dp
      axis = MINLOC(ABS(ez), DIM=1)
      reference(axis) = 1.0_dp
      ex = reference - ez(axis)*ez
      dex = -dez(axis)*ez - ez(axis)*dez
      length = SQRT(SUM(ex**2))
      dlength = DOT_PRODUCT(ex, dex)/length
      dex = (dex - ex*dlength/length)/length
      ex = ex/length
      ey = cross(ez, ex)
      dey = cross(dez, ex) + cross(ez, dex)
      nq = MAX(32, MAXVAL(la) + MAXVAL(lb) + 2)
      nphi = MAXVAL(la) + MAXVAL(lb) + 1
      ALLOCATE (nodes(nq), weights(nq), ya(SIZE(la)), yb(SIZE(lb)), &
                dya(SIZE(la)), dyb(SIZE(lb)), ga(3, SIZE(la)), gb(3, SIZE(lb)))
      CALL polar_rule(nodes, weights)
      value = 0.0_dp
      response = 0.0_dp
      tmax = distance + a
      left = -1.0_dp
      dleft = 0.0_dp
      regular = ABS(distance - a) < b
      DO WHILE (left < 1.0_dp)
         right = 1.0_dp
         dright = 0.0_dp
         t = SQRT(MAX(0.0_dp, (distance - a)**2 + 2.0_dp*distance*a*(left + 1.0_dp)))
         IF (regular .AND. b < tmax) THEN
            bound = b
            dbound = db
         ELSE IF (.NOT. regular .AND. 2.0_dp*t < tmax) THEN
            bound = 2.0_dp*t
            dbound = 2.0_dp*((distance - a)*(dt - da) + &
                             (dt*a + distance*da)*(left + 1.0_dp) + distance*a*dleft)/t
         ELSE
            bound = tmax
            dbound = dt + da
         END IF
         IF (bound < tmax) THEN
            right = ((bound - a)*(bound + a) - distance**2)/(2.0_dp*distance*a)
            dright = (bound*dbound - a*da - distance*dt)/(distance*a) - &
                     right*(dt/distance + da/a)
         END IF
         CPASSERT(right > left .AND. right <= 1.0_dp)
         DO iq = 1, nq
            mu = 0.5_dp*((1.0_dp - nodes(iq))*left + (1.0_dp + nodes(iq))*right)
            dmu = 0.5_dp*((1.0_dp - nodes(iq))*dleft + (1.0_dp + nodes(iq))*dright)
            weight = 0.5_dp*(right - left)*weights(iq)*twopi/REAL(nphi, dp)
            dweight = 0.5_dp*(dright - dleft)*weights(iq)*twopi/REAL(nphi, dp)
            sine = SQRT(MAX(0.0_dp, 1.0_dp - mu**2))
            CPASSERT(sine > 0.0_dp)
            DO iphi = 1, nphi
               azimuth = twopi*REAL(iphi - 1, dp)/REAL(nphi, dp)
               u = mu*ez + sine*(COS(azimuth)*ex + SIN(azimuth)*ey)
               du = dmu*ez + mu*dez - mu*dmu/sine*(COS(azimuth)*ex + SIN(azimuth)*ey) + &
                    sine*(COS(azimuth)*dex + SIN(azimuth)*dey)
               x = d + a*u
               dx = dd + da*u + a*du
               norm_x = SQRT(SUM(x**2))
               CALL solid_fields(u, la, ia, ya, ga)
               CALL solid_fields(x, lb, ib, yb, gb)
               dya(:) = MATMUL(du, ga)
               dyb(:) = MATMUL(dx, gb)
               DO channel = 1, SIZE(lb)
                  IF (norm_x < b) THEN
                     factor = fourpi/REAL(2*lb(channel) + 1, dp)/b**(lb(channel) + 1)
                     dfactor = -factor*(lb(channel) + 1)*db/b
                  ELSE
                     factor = fourpi/REAL(2*lb(channel) + 1, dp)* &
                              (b/norm_x)**lb(channel)/norm_x**(lb(channel) + 1)
                     dfactor = factor*(lb(channel)*db/b - &
                                       (2*lb(channel) + 1)*DOT_PRODUCT(x, dx)/norm_x**2)
                  END IF
                  scale = factor*yb(channel)
                  value(:, channel) = value(:, channel) + weight*ya*scale
                  response(:, channel) = response(:, channel) + dweight*ya*scale + &
                                         weight*(dya*scale + ya*(dfactor*yb(channel) + factor*dyb(channel)))
               END DO
            END DO
         END DO
         left = right
         dleft = dright
         regular = .FALSE.
      END DO
   END SUBROUTINE translated_shell

! **************************************************************************************************
!> \brief Native regular solid harmonics and Cartesian gradients, including their origin limits.
!> \param x ...
!> \param ls ...
!> \param indices ...
!> \param values ...
!> \param gradients ...
! **************************************************************************************************
   SUBROUTINE solid_fields(x, ls, indices, values, gradients)
      REAL(dp), INTENT(IN)                               :: x(3)
      INTEGER, INTENT(IN)                                :: ls(:), indices(:)
      REAL(dp), INTENT(OUT)                              :: values(:), gradients(:, :)

      INTEGER                                            :: axis, channel, i, l, m, powers(3), &
                                                            reduced(3)
      REAL(dp)                                           :: factor

      values = 0.0_dp
      gradients = 0.0_dp
      DO channel = 1, SIZE(ls)
         l = ls(channel)
         m = indices(channel) - l*l
         DO i = 1, nco(l)
            powers = indco(:, ncoset(l - 1) + i)
            factor = orbtramat(l)%slm(m, i)
            values(channel) = values(channel) + factor*PRODUCT(x**powers)
            DO axis = 1, 3
               IF (powers(axis) == 0) CYCLE
               reduced = powers
               reduced(axis) = reduced(axis) - 1
               gradients(axis, channel) = gradients(axis, channel) + &
                                          factor*powers(axis)*PRODUCT(x**reduced)
            END DO
         END DO
      END DO
   END SUBROUTINE solid_fields

! **************************************************************************************************
!> \brief Gauss-Legendre nodes and weights on (-1,1).
!> \param nodes ...
!> \param weights ...
! **************************************************************************************************
   SUBROUTINE polar_rule(nodes, weights)
      REAL(dp), INTENT(OUT)                              :: nodes(:), weights(:)

      INTEGER                                            :: i, iteration, j, n
      REAL(dp)                                           :: derivative, p, previous, saved, shift, x

      n = SIZE(nodes)
      CPASSERT(SIZE(weights) == n)
      DO i = 1, (n + 1)/2
         x = COS(pi*(REAL(i, dp) - 0.25_dp)/(REAL(n, dp) + 0.5_dp))
         DO iteration = 1, 20
            p = 1.0_dp
            previous = 0.0_dp
            DO j = 1, n
               saved = p
               p = ((2*j - 1)*x*p - (j - 1)*previous)/REAL(j, dp)
               previous = saved
            END DO
            derivative = n*(x*p - previous)/(x*x - 1.0_dp)
            shift = p/derivative
            x = x - shift
            IF (ABS(shift) < 2.0_dp*EPSILON(x)) EXIT
         END DO
         CPASSERT(iteration <= 20)
         nodes(i) = -x
         nodes(n + 1 - i) = x
         weights(i) = 2.0_dp/((1.0_dp - x*x)*derivative**2)
         weights(n + 1 - i) = weights(i)
      END DO
   END SUBROUTINE polar_rule

! **************************************************************************************************
!> \brief Cartesian vector product.
!> \param a ...
!> \param b ...
!> \return ...
! **************************************************************************************************
   FUNCTION cross(a, b) RESULT(c)
      REAL(dp), INTENT(IN)                               :: a(3), b(3)
      REAL(dp)                                           :: c(3)

      c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
   END FUNCTION cross

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
         CPASSERT(PRESENT(position_response))
         CPASSERT(PRESENT(volume_response))
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
