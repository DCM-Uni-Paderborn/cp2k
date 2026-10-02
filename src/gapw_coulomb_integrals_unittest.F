!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM gapw_coulomb_integrals_unittest
   USE gapw_coulomb_integrals,          ONLY: gapw_augmentation_integrals,&
                                              gapw_compensated_potential,&
                                              gapw_coulomb_init,&
                                              gapw_periodic_augmentation_integrals
   USE kinds,                           ONLY: dp
   USE lebedev,                         ONLY: deallocate_lebedev_grids,&
                                              get_number_of_lebedev_grid,&
                                              init_lebedev_grids,&
                                              lebedev_grid
   USE libcp2k,                         ONLY: cp2k_finalize,&
                                              cp2k_init
   USE mathconstants,                   ONLY: fourpi,&
                                              pi
   USE message_passing,                 ONLY: mp_comm_type,&
                                              mp_comm_world
   USE mt_util,                         ONLY: mt2d_kernel_strain
   USE pw_grid_types,                   ONLY: FULLSPACE,&
                                              HALFSPACE,&
                                              pw_grid_type
   USE pw_grids,                        ONLY: pw_grid_create,&
                                              pw_grid_release
   USE pw_methods,                      ONLY: do_accurate_sum,&
                                              do_standard_sum,&
                                              pw_integral_ab
   USE pw_poisson_methods,              ONLY: pw_poisson_set,&
                                              pw_poisson_solve
   USE pw_poisson_types,                ONLY: pw_poisson_mt,&
                                              pw_poisson_type
   USE pw_pool_types,                   ONLY: pw_pool_create
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r1d_gs_type
   USE spherical_harmonics,             ONLY: y_lm
   USE uspp_augmentation,               ONLY: uspp_augmentation_type

   IMPLICIT NONE
   REAL(dp) :: max_error

   CALL cp2k_init()
   CALL gapw_coulomb_init(6, 4)
   max_error = 0.0_dp
   CALL check_pw_integral_and_poisson()
   CALL check_mt2d_kernel()
   CALL check_multipoles()
   CALL check_augmentation()
   CALL check_periodic()
   PRINT '(A,ES16.7)', 'GAPW compensated Coulomb maximum scaled error: ', max_error
   CALL cp2k_finalize()

CONTAINS

! **************************************************************************************************
!> \brief Production PW multiplicity/volume oracle and auxiliary-density Poisson path equivalence.
! **************************************************************************************************
   SUBROUTINE check_pw_integral_and_poisson()
      REAL(dp), PARAMETER                                :: step = 5.e-4_dp

      COMPLEX(dp)                                        :: a_mode, b_mode
      INTEGER                                            :: axis, geometry, ig, sign_step, span, &
                                                            sumtype, test
      REAL(dp) :: cell(3, 3), derivative, energy, energy_dv, energy_v, energy_v_dv, expected, &
         expected_real, kernel, reference, shifted_g(3), strain_energy(4), &
         strain_steps(4) = [step, -step, 2.0_dp*step, -2.0_dp*step], stress(3, 3), &
         stress_dv(3, 3), stress_only(3, 3), stress_v(3, 3), stress_v_dv(3, 3)
      TYPE(mp_comm_type)                                 :: comm
      TYPE(pw_c1d_gs_type)                               :: a, b, gradient(3), gradient_v_dv(3), &
                                                            potential, potential_v_dv
      TYPE(pw_grid_type), POINTER                        :: grid
      TYPE(pw_poisson_type)                              :: poisson
      TYPE(pw_r1d_gs_type)                               :: a_real, b_real

      CALL comm%set_handle(mp_comm_world%get_handle())
      DO geometry = 1, 2
         cell = 0.0_dp
         cell(1, 1) = 1.0_dp
         cell(2, 2) = 1.0_dp
         cell(3, 3) = 7.0_dp
         IF (geometry == 2) THEN
            ! Same volume, with an appreciable truncation factor exp(-pi).
            DO axis = 1, 3
               cell(axis, axis) = 7.0_dp**(1.0_dp/3.0_dp)
            END DO
         END IF
         a_mode = CMPLX(0.4_dp, -0.2_dp, dp)
         b_mode = CMPLX(-0.7_dp, 0.9_dp, dp)
         NULLIFY (grid)
         DO span = FULLSPACE, HALFSPACE, -1
            ! Reciprocal inputs and outputs need no FFT. The production grid still
            ! supplies the distributed G0 ownership and half/full-space mappings.
            CALL pw_grid_create(grid, comm, cell, grid_span=span, npts=[5, 5, 5], odd=.TRUE., fft_usage=.FALSE.)
            CALL a%create(grid)
            CALL b%create(grid)
            CALL a_real%create(grid)
            CALL b_real%create(grid)
            CALL potential%create(grid)
            CALL potential_v_dv%create(grid)
            DO axis = 1, 3
               CALL gradient(axis)%create(grid)
               CALL gradient_v_dv(axis)%create(grid)
            END DO
            ALLOCATE (poisson%pw_pools(1))
            CALL pw_pool_create(poisson%pw_pools(1)%pool, grid)
            poisson%pw_level = 1
            poisson%parameters%solver = pw_poisson_mt
            poisson%parameters%periodic = [1, 1, 0]
            poisson%parameters%mt_alpha = 1.0_dp
            CALL pw_poisson_set(poisson, cell_hmat=cell, use_level=1)
            DO test = 1, 3
               ! Constant G0 only, G0 plus a Hermitian mode pair, and mode pair only.
               a%array = 0.0_dp
               b%array = 0.0_dp
               expected = 0.0_dp
               expected_real = 0.0_dp
               reference = 0.0_dp
               IF (test /= 3) THEN
                  IF (grid%have_g0) THEN
                     a%array(1) = 2.0_dp
                     b%array(1) = 3.0_dp
                  END IF
                  expected = 6.0_dp
                  expected_real = 6.0_dp
                  reference = 6.0_dp*mt2d_kernel_reference([0.0_dp, 0.0_dp, 0.0_dp], 3, cell(3, 3))
               END IF
               IF (test /= 1) THEN
                  DO ig = 1, grid%ngpts_cut_local
                     IF (ANY(grid%g_hat(2:3, ig) /= 0)) CYCLE
                     IF (grid%g_hat(1, ig) == 1) THEN
                        a%array(ig) = a_mode
                        b%array(ig) = b_mode
                     ELSE IF (grid%g_hat(1, ig) == -1) THEN
                        a%array(ig) = CONJG(a_mode)
                        b%array(ig) = CONJG(b_mode)
                     END IF
                  END DO
                  expected = expected + 2.0_dp*REAL(CONJG(a_mode)*b_mode, dp)
                  expected_real = expected_real + 2.0_dp*REAL(a_mode, dp)*REAL(b_mode, dp)
                  kernel = mt2d_kernel_reference([2.0_dp*pi/cell(1, 1), 0.0_dp, 0.0_dp], 3, cell(3, 3))
                  reference = reference + 2.0_dp*REAL(CONJG(a_mode)*b_mode, dp)*kernel
               END IF
               a_real%array = REAL(a%array, dp)
               b_real%array = REAL(b%array, dp)
               ! G0 must have weight one, the mode pair weight two, and both must
               ! receive V=7 for an integral. just_sum retains the unscaled oracle.
               DO sumtype = do_accurate_sum, do_standard_sum
                  CALL check_error(pw_integral_ab(a, b, sumtype=sumtype), 7.0_dp*expected, 3.e-14_dp)
                  CALL check_error(pw_integral_ab(a, b, sumtype=sumtype, just_sum=.TRUE.), expected, 3.e-14_dp)
                  CALL check_error(pw_integral_ab(a_real, b_real, sumtype=sumtype), 7.0_dp*expected_real, 3.e-14_dp)
                  CALL check_error(pw_integral_ab(a_real, b_real, sumtype=sumtype, just_sum=.TRUE.), &
                                   expected_real, 3.e-14_dp)
                  CALL check_error(pw_integral_ab(a, b_real, sumtype=sumtype), 7.0_dp*expected_real, 3.e-14_dp)
                  CALL check_error(pw_integral_ab(a_real, b, sumtype=sumtype), 7.0_dp*expected_real, 3.e-14_dp)
               END DO
               CALL pw_poisson_solve(poisson, a, energy, h_stress=stress, aux_density=b)
               CALL pw_poisson_solve(poisson, a, energy_dv, dvhartree=gradient, h_stress=stress_dv, aux_density=b)
               CALL pw_poisson_solve(poisson, a, energy_v, potential, h_stress=stress_v, aux_density=b)
               CALL pw_poisson_solve(poisson, a, energy_v_dv, potential_v_dv, dvhartree=gradient_v_dv, &
                                     h_stress=stress_v_dv, aux_density=b)
               CALL pw_poisson_solve(poisson, a, h_stress=stress_only, aux_density=b)
               CALL check_error(energy, 0.5_dp*7.0_dp*reference, 3.e-13_dp)
               CALL check_error(energy_dv, energy, 3.e-13_dp)
               CALL check_error(energy_v, energy, 3.e-13_dp)
               CALL check_error(energy_v_dv, energy, 3.e-13_dp)
               CALL check_error(MAXVAL(ABS(stress_dv - stress)), 0.0_dp, 3.e-11_dp)
               CALL check_error(MAXVAL(ABS(stress_v - stress)), 0.0_dp, 3.e-11_dp)
               CALL check_error(MAXVAL(ABS(stress_v_dv - stress)), 0.0_dp, 3.e-11_dp)
               CALL check_error(MAXVAL(ABS(stress_only - stress)), 0.0_dp, 3.e-11_dp)
               IF (test == 2) THEN
                  ! Fixed Fourier charges Q=V*rho give E(epsilon)=V/2/(1+epsilon)
                  ! times the strained kernel sum. This absolute production-virial
                  ! oracle includes both G0/volume and the nonzero-mode kernel term.
                  DO sign_step = 1, 4
                     shifted_g = [2.0_dp*pi/cell(1, 1)/(1.0_dp + strain_steps(sign_step)), 0.0_dp, 0.0_dp]
                     reference = 6.0_dp*mt2d_kernel_reference([0.0_dp, 0.0_dp, 0.0_dp], 3, cell(3, 3)) &
                                 + 2.0_dp*REAL(CONJG(a_mode)*b_mode, dp)*mt2d_kernel_reference(shifted_g, 3, cell(3, 3))
                     strain_energy(sign_step) = 0.5_dp*7.0_dp*reference/(1.0_dp + strain_steps(sign_step))
                  END DO
                  derivative = (8.0_dp*(strain_energy(1) - strain_energy(2)) &
                                - strain_energy(3) + strain_energy(4))/(12.0_dp*step)
                  CALL check_error(stress(1, 1) + derivative, 0.0_dp, 3.e-9_dp)
               END IF
               DO ig = 1, grid%ngpts_cut_local
                  kernel = mt2d_kernel_reference(grid%g(:, ig), 3, cell(3, 3))
                  CALL check_error(ABS(potential%array(ig) - kernel*a%array(ig)), 0.0_dp, 3.e-13_dp)
                  CALL check_error(ABS(potential_v_dv%array(ig) - kernel*a%array(ig)), 0.0_dp, 3.e-13_dp)
                  DO axis = 1, 3
                     CALL check_error(ABS(gradient(axis)%array(ig) &
                                          - CMPLX(0.0_dp, grid%g(axis, ig), dp)*kernel*a%array(ig)), 0.0_dp, 3.e-13_dp)
                     CALL check_error(ABS(gradient_v_dv(axis)%array(ig) &
                                          - CMPLX(0.0_dp, grid%g(axis, ig), dp)*kernel*a%array(ig)), 0.0_dp, 3.e-13_dp)
                  END DO
               END DO
            END DO
            CALL poisson%release()
            CALL a%release()
            CALL b%release()
            CALL a_real%release()
            CALL b_real%release()
            CALL potential%release()
            CALL potential_v_dv%release()
            DO axis = 1, 3
               CALL gradient(axis)%release()
               CALL gradient_v_dv(axis)%release()
            END DO
            CALL pw_grid_release(grid)
         END DO
      END DO
   END SUBROUTINE check_pw_integral_and_poisson

! **************************************************************************************************
!> \brief Independent affine-strain differences of the MT2D kernel and its volume contribution.
! **************************************************************************************************
   SUBROUTINE check_mt2d_kernel()
      REAL(dp), PARAMETER                                :: step = 3.e-4_dp

      COMPLEX(dp)                                        :: charge, charge_aux
      INTEGER                                            :: alpha, beta, ip, iz, jp, miller(3), &
                                                            mode, modes(3, 7), sign_step, test
      REAL(dp) :: bilinear, derivative, energy(4), g(3), h_inv(3, 3), kernel, lengths(3), &
         shifted_g(3), shifted_kernel(4), shifts(4) = [step, -step, 2.0_dp*step, -2.0_dp*step], &
         strain(3, 3), value, volume, volume_shifted

      ! Rotate the Cartesian vacuum axis. A shear within the periodic plane is
      ! permitted; only normal/periodic mixing is excluded by the MT2D kernel.
      DO iz = 1, 3
         ip = MOD(iz, 3) + 1
         jp = MOD(iz + 1, 3) + 1
         lengths = [9.0_dp, 13.0_dp, 17.0_dp]
         volume = PRODUCT(lengths)
         h_inv = 0.0_dp
         DO alpha = 1, 3
            h_inv(alpha, alpha) = 1.0_dp/lengths(alpha)
         END DO
         h_inv(ip, jp) = -1.4_dp/(lengths(ip)*lengths(jp))
         modes = 0
         ! G=0, both odd and even normal-only modes, and several modes with
         ! appreciable exponential truncation reject the periodic K**2 formula.
         modes(iz, 2:3) = [1, 2]
         modes(ip, 4) = 1
         modes(:, 5) = [1, 2, 1]
         modes(:, 6) = [2, -1, 2]
         modes(:, 7) = [-1, 2, -3]
         DO mode = 1, SIZE(modes, 2)
            miller = modes(:, mode)
            g = 2.0_dp*pi*MATMUL(TRANSPOSE(h_inv), REAL(miller, dp))
            CALL mt2d_kernel_strain(g, iz, lengths(iz), kernel, strain)
            CALL check_error(kernel, mt2d_kernel_reference(g, iz, lengths(iz)), 3.e-14_dp)
            CALL check_error(MAXVAL(ABS(strain(iz, :))), 0.0_dp, 0.0_dp)
            CALL check_error(MAXVAL(ABS(strain(:, iz))), 0.0_dp, 0.0_dp)
            IF (mode <= 3) CALL check_error(MAXVAL(ABS(strain)), 0.0_dp, 0.0_dp)
            DO alpha = 1, 3
               IF (alpha == iz) CYCLE
               DO beta = 1, 3
                  IF (beta == iz) CYCLE
                  DO sign_step = 1, 4
                     ! g'=(I+epsilon_ab)^(-T) g, evaluated independently of dK.
                     shifted_g = g
                     IF (alpha == beta) THEN
                        shifted_g(beta) = g(beta)/(1.0_dp + shifts(sign_step))
                     ELSE
                        shifted_g(beta) = g(beta) - shifts(sign_step)*g(alpha)
                     END IF
                     shifted_kernel(sign_step) = mt2d_kernel_reference(shifted_g, iz, lengths(iz))
                  END DO
                  derivative = (8.0_dp*(shifted_kernel(1) - shifted_kernel(2)) &
                                - shifted_kernel(3) + shifted_kernel(4))/(12.0_dp*step)
                  CALL check_error(strain(alpha, beta), derivative, 3.e-9_dp)
                  ! Fixed Fourier charges give rho_G=Q_G/volume. Check both a
                  ! self energy and an auxiliary-density bilinear energy, and
                  ! include G=0 separately with its unit reciprocal multiplicity.
                  charge = CMPLX(0.7_dp, -0.3_dp, dp)
                  DO test = 1, 2
                     charge_aux = charge
                     IF (test == 2) charge_aux = CMPLX(-0.4_dp, 0.6_dp, dp)
                     bilinear = REAL(CONJG(charge_aux)*charge, dp)
                     IF (mode /= 1) bilinear = 2.0_dp*bilinear
                     DO sign_step = 1, 4
                        volume_shifted = volume
                        IF (alpha == beta) volume_shifted = volume*(1.0_dp + shifts(sign_step))
                        energy(sign_step) = 0.5_dp*bilinear*shifted_kernel(sign_step)/volume_shifted
                     END DO
                     derivative = (8.0_dp*(energy(1) - energy(2)) - energy(3) + energy(4))/(12.0_dp*step)
                     value = -0.5_dp*bilinear*strain(alpha, beta)/volume
                     IF (alpha == beta) value = value + 0.5_dp*bilinear*kernel/volume
                     CALL check_error(value, -derivative, 3.e-11_dp)
                  END DO
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE check_mt2d_kernel

! **************************************************************************************************
!> \brief Original MT2D influence kernel, without using the analytic strain implementation.
!> \param g Strained Cartesian reciprocal vector.
!> \param iz Cartesian vacuum direction.
!> \param slab_size Fixed thickness of the slab cell.
!> \return MT2D influence kernel.
! **************************************************************************************************
   FUNCTION mt2d_kernel_reference(g, iz, slab_size) RESULT(kernel)
      REAL(dp), INTENT(IN)                               :: g(3)
      INTEGER, INTENT(IN)                                :: iz
      REAL(dp), INTENT(IN)                               :: slab_size
      REAL(dp)                                           :: kernel

      REAL(dp)                                           :: g2, g_parallel(3), radius

      g2 = SUM(g**2)
      IF (g2 == 0.0_dp) THEN
         kernel = pi*slab_size**2/2.0_dp
      ELSE
         g_parallel = g
         g_parallel(iz) = 0.0_dp
         radius = slab_size/2.0_dp
         kernel = fourpi/g2 - fourpi/g2*COS(g(iz)*radius)*EXP(-SQRT(SUM(g_parallel**2))*radius)
      END IF
   END FUNCTION mt2d_kernel_reference

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE check_multipoles()
      REAL(dp), PARAMETER                                :: alpha = 1.3_dp, eta = 0.55_dp, &
                                                            step = 1.e-3_dp

      INTEGER                                            :: axis, j, l, m, n, power, sample
      REAL(dp) :: angular(1), c, density, direction(3, 1), gm(3, 13), gmm(3, 13), gp(3, 13), &
         gpp(3, 13), grad(3, 13), hess(3, 3, 13), higher(13), point(3), points(3, 7), r, &
         reference, s, shifted(3), tail_difference, value(13), vm(13), vmm(13), vp(13), vpp(13)

      points(:, 1) = 0.0_dp
      points(:, 2) = [1.e-8_dp, -2.e-8_dp, 3.e-8_dp]
      points(:, 3) = [0.31_dp, -0.52_dp, 0.73_dp]
      points(:, 4) = [0.0_dp, 0.0_dp, 1.2_dp]
      points(:, 7) = [8.0_dp, -3.0_dp, 2.0_dp]
      DO l = 0, 6
         n = 2*l + 1
         DO power = 0, 3
            points(:, 5) = [SQRT(REAL(l + power + 2, dp)/eta)*(1.0_dp - 1.e-8_dp), 0.0_dp, 0.0_dp]
            points(:, 6) = [SQRT(REAL(l + power + 2, dp)/eta)*(1.0_dp + 1.e-8_dp), 0.0_dp, 0.0_dp]
            c = (eta/alpha)**(REAL(l, dp) + 1.5_dp)
            DO j = 1, power
               c = c*(REAL(l + j, dp) + 0.5_dp)/alpha
            END DO
            DO sample = 1, SIZE(points, 2)
               point = points(:, sample)
               s = SUM(point**2)
               r = SQRT(s)
               CALL gapw_compensated_potential(alpha, power, l, eta, point, value(:n), grad(:, :n), hess(:, :, :n))
               DO m = -l, l
                  angular = 0.0_dp
                  IF (r > 0.0_dp) THEN
                     direction(:, 1) = point/r
                     CALL y_lm(direction, angular, l, m)
                     angular = angular*r**l
                  ELSE IF (l == 0) THEN
                     angular = 1.0_dp/SQRT(fourpi)
                  END IF
                  density = (s**power*EXP(-alpha*s) - c*EXP(-eta*s))*angular(1)
                  reference = -fourpi*density
                  CALL check_error(hess(1, 1, m + l + 1) + hess(2, 2, m + l + 1) + &
                                   hess(3, 3, m + l + 1), reference, 2.e-11_dp)
               END DO
               DO axis = 1, 3
                  ! Fourth-order differences permit a larger step without magnifying
                  ! roundoff where the direct and complementary-tail formulas join.
                  shifted = point
                  shifted(axis) = shifted(axis) + step
                  CALL gapw_compensated_potential(alpha, power, l, eta, shifted, vp(:n), gp(:, :n))
                  shifted(axis) = shifted(axis) - 2.0_dp*step
                  CALL gapw_compensated_potential(alpha, power, l, eta, shifted, vm(:n), gm(:, :n))
                  shifted(axis) = point(axis) + 2.0_dp*step
                  CALL gapw_compensated_potential(alpha, power, l, eta, shifted, vpp(:n), gpp(:, :n))
                  shifted(axis) = point(axis) - 2.0_dp*step
                  CALL gapw_compensated_potential(alpha, power, l, eta, shifted, vmm(:n), gmm(:, :n))
                  DO m = 1, n
                     CALL check_error((8.0_dp*(vp(m) - vm(m)) - vpp(m) + vmm(m))/(12.0_dp*step), &
                                      grad(axis, m), 2.e-8_dp)
                     DO j = 1, 3
                        CALL check_error((8.0_dp*(gp(j, m) - gm(j, m)) - gpp(j, m) + gmm(j, m))/(12.0_dp*step), &
                                         hess(j, axis, m), 2.e-8_dp)
                     END DO
                  END DO
               END DO
               ! -d/d(alpha) raises the radial power, including the compensation moment.
               CALL gapw_compensated_potential(alpha + step, power, l, eta, point, vp(:n), gp(:, :n))
               CALL gapw_compensated_potential(alpha - step, power, l, eta, point, vm(:n), gm(:, :n))
               CALL gapw_compensated_potential(alpha + 2.0_dp*step, power, l, eta, point, vpp(:n), gp(:, :n))
               CALL gapw_compensated_potential(alpha - 2.0_dp*step, power, l, eta, point, vmm(:n), gm(:, :n))
               CALL gapw_compensated_potential(alpha, power + 1, l, eta, point, higher(:n), gp(:, :n))
               DO m = 1, n
                  CALL check_error(-(8.0_dp*(vp(m) - vm(m)) - vpp(m) + vmm(m))/(12.0_dp*step), &
                                   higher(m), 2.e-8_dp)
               END DO
            END DO
         END DO
      END DO
      ! Independent monopole solution, including the exponentially small far field.
      DO sample = 1, SIZE(points, 2)
         point = points(:, sample)
         r = SQRT(SUM(point**2))
         IF (r == 0.0_dp) THEN
            reference = SQRT(pi)/alpha*(1.0_dp - SQRT(eta/alpha))
         ELSE
            tail_difference = ERFC(SQRT(eta)*r) - ERFC(SQRT(alpha)*r)
            reference = pi/(2.0_dp*alpha**1.5_dp)*tail_difference/r
         END IF
         CALL gapw_compensated_potential(alpha, 0, 0, eta, point, value(:1), grad(:, :1))
         IF (r > 1.e-6_dp) THEN
            CALL check_error(value(1)/reference, 1.0_dp, 2.e-10_dp)
         ELSE
            ! Avoid cancellation in the independent erfc difference near r=0.
            reference = SQRT(pi)/alpha*(1.0_dp - SQRT(eta/alpha))
            CALL check_error(value(1), reference, 2.e-12_dp)
         END IF
      END DO
      CALL gapw_compensated_potential(alpha, 0, 6, alpha, point, value, grad, hess)
      IF (ANY(value /= 0.0_dp) .OR. ANY(grad /= 0.0_dp) .OR. ANY(hess /= 0.0_dp)) THEN
         ERROR STOP 'Identical compensated Gaussian did not cancel exactly'
      END IF
   END SUBROUTINE check_multipoles

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE check_augmentation()
      INTEGER, PARAMETER                                 :: nr = 801
      REAL(dp), PARAMETER                                :: a = 1.1_dp, b = 1.3_dp, dr = 0.01_dp, &
                                                            eta = 0.55_dp, step = 2.e-5_dp

      INTEGER                                            :: axis, ir, lg, sample
      REAL(dp) :: derivative(2, 2, 1, 3), dref(3), integral(2, 2, 1), local_error, minus(2, 2, 1), &
         plus(2, 2, 1), point(3), radial_weights(nr), reference, shifted(3)
      TYPE(uspp_augmentation_type)                       :: aug

      aug%nproj = 2
      aug%lmax = 2
      ALLOCATE (aug%r(nr), aug%weighted_q(nr, 3, 0:2), aug%radial_index(2), aug%lbeta(2), aug%gaunt(2, 2, 9))
      aug%radial_index(:) = [1, 2]
      aug%lbeta(:) = [0, 1]
      aug%gaunt = 0.0_dp
      aug%gaunt(1, 1, 1) = 1.0_dp/SQRT(fourpi)
      aug%gaunt(1, 2, 3) = 1.0_dp/SQRT(fourpi)
      aug%gaunt(2, 1, 3) = aug%gaunt(1, 2, 3)
      DO ir = 1, nr
         aug%r(ir) = REAL(ir - 1, dp)*dr
      END DO
      radial_weights = 2.0_dp*dr/3.0_dp
      radial_weights(2:nr - 1:2) = 4.0_dp*dr/3.0_dp
      radial_weights(1) = dr/3.0_dp
      radial_weights(nr) = dr/3.0_dp
      aug%weighted_q = 0.0_dp
      aug%weighted_q(:, 1, 0) = radial_weights*aug%r**2*EXP(-a*aug%r**2)
      aug%weighted_q(:, 2, 1) = aug%r*aug%weighted_q(:, 1, 0)
      CALL init_lebedev_grids()
      lg = get_number_of_lebedev_grid(n=194)
      DO sample = 1, 3
         point = REAL(sample - 1, dp)*[0.31_dp, -0.52_dp, 0.73_dp]
         CALL gapw_augmentation_integrals(aug, point, b, 0, 0, eta, lebedev_grid(lg)%r, &
                                          fourpi*lebedev_grid(lg)%w, integral, derivative)
         CALL gaussian_reference(point, a, b, eta, reference, dref)
         CALL check_error(integral(1, 1, 1), reference, 2.e-10_dp)
         ! s-p augmentation is sqrt(3)/(2*a) times the center derivative of s-s.
         CALL check_error(integral(1, 2, 1), SQRT(3.0_dp)/(2.0_dp*a)*dref(3), 2.e-10_dp)
         IF (ANY(integral(:, :, 1) /= TRANSPOSE(integral(:, :, 1)))) ERROR STOP 'Asymmetric Coulomb tensor'
         DO axis = 1, 3
            CALL check_error(derivative(1, 1, 1, axis), dref(axis), 2.e-10_dp)
            IF (sample /= 2) CYCLE
            shifted = point
            shifted(axis) = shifted(axis) + step
            CALL gapw_augmentation_integrals(aug, shifted, b, 0, 0, eta, lebedev_grid(lg)%r, &
                                             fourpi*lebedev_grid(lg)%w, plus)
            shifted(axis) = shifted(axis) - 2.0_dp*step
            CALL gapw_augmentation_integrals(aug, shifted, b, 0, 0, eta, lebedev_grid(lg)%r, &
                                             fourpi*lebedev_grid(lg)%w, minus)
            local_error = MAXVAL(ABS((plus(:, :, 1) - minus(:, :, 1))/(2.0_dp*step) - derivative(:, :, 1, axis)))
            CALL check_error(local_error, 0.0_dp, 2.e-9_dp)
         END DO
      END DO
      CALL deallocate_lebedev_grids()
   END SUBROUTINE check_augmentation

! **************************************************************************************************
!> \brief ...
!> \param point ...
!> \param a ...
!> \param b ...
!> \param eta ...
!> \param value ...
!> \param gradient ...
! **************************************************************************************************
   SUBROUTINE gaussian_reference(point, a, b, eta, value, gradient)
      REAL(dp), INTENT(IN)                               :: point(3), a, b, eta
      REAL(dp), INTENT(OUT)                              :: value, gradient(3)

      REAL(dp)                                           :: c, delta, r, radial, u, v

      r = SQRT(SUM(point**2))
      u = SQRT(a*b/(a + b))
      v = SQRT(a*eta/(a + eta))
      c = (pi*pi/(a*b))**1.5_dp/(fourpi*SQRT(fourpi))
      gradient = 0.0_dp
      IF (r == 0.0_dp) THEN
         value = c*2.0_dp/SQRT(pi)*(u - v)
      ELSE
         delta = ERF(u*r) - ERF(v*r)
         value = c*delta/r
         radial = c*(2.0_dp/SQRT(pi)*(u*EXP(-(u*r)**2) - v*EXP(-(v*r)**2))/r - delta/r**2)
         gradient = radial*point/r
      END IF
   END SUBROUTINE gaussian_reference

! **************************************************************************************************
!> \brief Triclinic lattice, wrapped representatives and strain against analytic Gaussian pairs.
! **************************************************************************************************
   SUBROUTINE check_periodic()
      INTEGER, PARAMETER                                 :: nr = 121
      REAL(dp), PARAMETER                                :: a = 2.0_dp, b = 1.3_dp, dr = 0.05_dp, &
                                                            eta = 0.55_dp, step = 2.e-5_dp

      INTEGER                                            :: axis, axis2, boundary, ir, ix, iy, iz, &
                                                            lg, periodicity(3), sample
      REAL(dp) :: boundary_point(3), changed(3, 3), dref(3), fractional(3), h_inv(3, 3), &
         hmat(3, 3), image_point(3), inverse(3, 3), minus(1, 1, 1, 13), plus(1, 1, 1, 13), &
         point(3), radial_weights(nr), reference(13), shifted(1, 1, 1, 13), tensor(1, 1, 1, 13), &
         value
      TYPE(uspp_augmentation_type)                       :: aug, compact

      aug%nproj = 1
      aug%lmax = 0
      ALLOCATE (aug%r(nr), aug%weighted_q(nr, 1, 0:0), aug%radial_index(1), aug%lbeta(1), aug%gaunt(1, 1, 1))
      aug%radial_index = 1
      aug%lbeta = 0
      aug%gaunt = 1.0_dp/SQRT(fourpi)
      DO ir = 1, nr
         aug%r(ir) = REAL(ir - 1, dp)*dr
      END DO
      radial_weights = 2.0_dp*dr/3.0_dp
      radial_weights(2:nr - 1:2) = 4.0_dp*dr/3.0_dp
      radial_weights([1, nr]) = dr/3.0_dp
      aug%weighted_q(:, 1, 0) = radial_weights*aug%r**2*EXP(-a*aug%r**2)
      hmat = RESHAPE([6.2_dp, 0.0_dp, 0.0_dp, 1.3_dp, 6.7_dp, 0.0_dp, -0.8_dp, 0.6_dp, 7.1_dp], [3, 3])
      CALL triangular_inverse(hmat, h_inv)
      point = [0.31_dp, -0.52_dp, 0.73_dp]
      fractional = MATMUL(h_inv, point)
      reference = 0.0_dp
      DO iz = -4, 4
         DO iy = -4, 4
            DO ix = -4, 4
               image_point = point + MATMUL(hmat, REAL([ix, iy, iz], dp))
               CALL gaussian_reference(image_point, a, b, eta, value, dref)
               reference(1) = reference(1) + value
               reference(2:4) = reference(2:4) + dref
               DO axis = 1, 3
                  DO axis2 = 1, 3
                     ir = 4 + axis + 3*(axis2 - 1)
                     reference(ir) = reference(ir) + dref(axis)*image_point(axis2)
                  END DO
               END DO
            END DO
         END DO
      END DO
      CALL init_lebedev_grids()
      lg = get_number_of_lebedev_grid(n=194)
      CALL gapw_periodic_augmentation_integrals(aug, point, b, 0, 0, eta, hmat, h_inv, &
                                                lebedev_grid(lg)%r, fourpi*lebedev_grid(lg)%w, 1.e-14_dp, tensor)
      DO ir = 1, 13
         CALL check_error(tensor(1, 1, 1, ir), reference(ir), 3.e-9_dp)
      END DO
      CALL gapw_periodic_augmentation_integrals(aug, point + MATMUL(hmat, REAL([2, -3, 1], dp)), b, 0, 0, &
                                                eta, hmat, h_inv, lebedev_grid(lg)%r, fourpi*lebedev_grid(lg)%w, 1.e-16_dp, shifted)
      CALL check_error(MAXVAL(ABS(shifted - tensor)), 0.0_dp, 2.e-12_dp)
      ! A diagonal strain and a shear, with fixed fractional centers and rigid radial functions.
      DO sample = 1, 2
         axis = 1
         axis2 = sample
         changed = hmat
         changed(axis, :) = changed(axis, :) + step*hmat(axis2, :)
         CALL triangular_inverse(changed, inverse)
         CALL gapw_periodic_augmentation_integrals(aug, MATMUL(changed, fractional), b, 0, 0, eta, changed, &
                                                   inverse, lebedev_grid(lg)%r, fourpi*lebedev_grid(lg)%w, 1.e-14_dp, plus)
         changed = hmat
         changed(axis, :) = changed(axis, :) - step*hmat(axis2, :)
         CALL triangular_inverse(changed, inverse)
         CALL gapw_periodic_augmentation_integrals(aug, MATMUL(changed, fractional), b, 0, 0, eta, changed, &
                                                   inverse, lebedev_grid(lg)%r, fourpi*lebedev_grid(lg)%w, 1.e-14_dp, minus)
         CALL check_error((plus(1, 1, 1, 1) - minus(1, 1, 1, 1))/(2.0_dp*step), &
                          tensor(1, 1, 1, 4 + axis + 3*(axis2 - 1)), 3.e-9_dp)
      END DO
      ! Isolated, wire and slab boundaries, including nonconsecutive periodic axes.
      ! Move beyond a nonperiodic face to detect accidental minimum-image wrapping.
      DO boundary = 0, 2
         SELECT CASE (boundary)
         CASE (0)
            periodicity = [0, 0, 0]
         CASE (1)
            periodicity = [0, 1, 0]
         CASE (2)
            periodicity = [1, 0, 1]
         END SELECT
         DO axis = 1, 3
            IF (periodicity(axis) == 0) EXIT
         END DO
         boundary_point = point + 0.55_dp*hmat(:, axis)
         reference = 0.0_dp
         DO iz = -4*periodicity(3), 4*periodicity(3)
            DO iy = -4*periodicity(2), 4*periodicity(2)
               DO ix = -4*periodicity(1), 4*periodicity(1)
                  image_point = boundary_point + MATMUL(hmat, REAL([ix, iy, iz], dp))
                  CALL gaussian_reference(image_point, a, b, eta, value, dref)
                  reference(1) = reference(1) + value
                  reference(2:4) = reference(2:4) + dref
                  DO axis = 1, 3
                     DO axis2 = 1, 3
                        ir = 4 + axis + 3*(axis2 - 1)
                        reference(ir) = reference(ir) + dref(axis)*image_point(axis2)
                     END DO
                  END DO
               END DO
            END DO
         END DO
         CALL gapw_periodic_augmentation_integrals(aug, boundary_point, b, 0, 0, eta, hmat, h_inv, &
                                                   lebedev_grid(lg)%r, fourpi*lebedev_grid(lg)%w, 1.e-14_dp, tensor, &
                                                   periodicity=periodicity)
         DO ir = 1, 13
            CALL check_error(tensor(1, 1, 1, ir), reference(ir), 3.e-9_dp)
         END DO
         boundary_point = boundary_point + MATMUL(hmat, REAL(periodicity*[2, -3, 1], dp))
         CALL gapw_periodic_augmentation_integrals(aug, boundary_point, b, 0, 0, eta, hmat, h_inv, &
                                                   lebedev_grid(lg)%r, fourpi*lebedev_grid(lg)%w, 1.e-14_dp, shifted, &
                                                   periodicity=periodicity)
         CALL check_error(MAXVAL(ABS(shifted - tensor)), 0.0_dp, 2.e-12_dp)
      END DO
      ! A zero radial suffix may be omitted from weighted_q while r retains the
      ! original mesh for spline construction. Check all thirteen tensor entries.
      aug%weighted_q(nr/2 + 1:, :, :) = 0.0_dp
      CALL gapw_periodic_augmentation_integrals(aug, point, b, 0, 0, eta, hmat, h_inv, &
                                                lebedev_grid(lg)%r, fourpi*lebedev_grid(lg)%w, 1.e-14_dp, tensor)
      compact = aug
      DEALLOCATE (compact%weighted_q)
      ALLOCATE (compact%weighted_q(nr/2, 1, 0:0))
      compact%weighted_q(:, :, :) = aug%weighted_q(:nr/2, :, :)
      CALL gapw_periodic_augmentation_integrals(compact, point, b, 0, 0, eta, hmat, h_inv, &
                                                lebedev_grid(lg)%r, fourpi*lebedev_grid(lg)%w, 1.e-14_dp, shifted)
      CALL check_error(MAXVAL(ABS(shifted - tensor)), 0.0_dp, 1.e-14_dp)
      DEALLOCATE (compact%weighted_q)
      ALLOCATE (compact%weighted_q(0, 1, 0:0))
      CALL gapw_periodic_augmentation_integrals(compact, point, b, 0, 0, eta, hmat, h_inv, &
                                                lebedev_grid(lg)%r, fourpi*lebedev_grid(lg)%w, 1.e-14_dp, shifted)
      CALL check_error(MAXVAL(ABS(shifted)), 0.0_dp, 0.0_dp)
      CALL deallocate_lebedev_grids()
   END SUBROUTINE check_periodic

! **************************************************************************************************
!> \brief Inverse of the upper-triangular test cell.
!> \param hmat Test cell.
!> \param inverse Inverse cell.
! **************************************************************************************************
   SUBROUTINE triangular_inverse(hmat, inverse)
      REAL(dp), INTENT(IN)                               :: hmat(3, 3)
      REAL(dp), INTENT(OUT)                              :: inverse(3, 3)

      INTEGER                                            :: i

      inverse = 0.0_dp
      DO i = 1, 3
         inverse(i, i) = 1.0_dp/hmat(i, i)
      END DO
      inverse(1, 2) = -hmat(1, 2)/(hmat(1, 1)*hmat(2, 2))
      inverse(2, 3) = -hmat(2, 3)/(hmat(2, 2)*hmat(3, 3))
      inverse(1, 3) = (hmat(1, 2)*hmat(2, 3)/hmat(2, 2) - hmat(1, 3))/(hmat(1, 1)*hmat(3, 3))
   END SUBROUTINE triangular_inverse

! **************************************************************************************************
!> \brief ...
!> \param value ...
!> \param reference ...
!> \param tolerance ...
! **************************************************************************************************
   SUBROUTINE check_error(value, reference, tolerance)
      REAL(dp), INTENT(IN)                               :: value, reference, tolerance

      REAL(dp)                                           :: error

      error = ABS(value - reference)/MAX(1.0_dp, ABS(reference))
      max_error = MAX(max_error, error)
      IF (.NOT. (error <= tolerance)) THEN
         PRINT *, 'Coulomb check: value, reference, error, tolerance:', value, reference, error, tolerance
         ERROR STOP 'Compensated Coulomb integral check failed'
      END IF
   END SUBROUTINE check_error

END PROGRAM gapw_coulomb_integrals_unittest
