!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Inverse reconstruction identities, conditioning and analytic field/operator responses.
! **************************************************************************************************
PROGRAM upf_inverse_unittest
   USE atom_upf,                        ONLY: atom_upfpot_type
   USE basis_set_types,                 ONLY: allocate_gto_basis_set,&
                                              deallocate_gto_basis_set,&
                                              gto_basis_set_type
   USE gaussian_basis_fields,           ONLY: gaussian_basis_fields_evaluate,&
                                              gaussian_basis_fields_periodic
   USE gaussian_coulomb_fields,         ONLY: gaussian_coulomb_fields_evaluate,&
                                              gaussian_coulomb_fields_init,&
                                              gaussian_coulomb_fields_prepare,&
                                              gaussian_coulomb_fields_type
   USE ieee_arithmetic,                 ONLY: ieee_quiet_nan,&
                                              ieee_value
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: pi
   USE orbital_pointers,                ONLY: init_orbital_pointers,&
                                              ncoset
   USE qs_grid_atom,                    ONLY: grid_atom_type
   USE qs_harmonics_atom,               ONLY: harmonics_atom_type
   USE qs_upf_inverse_density,          ONLY: qs_upf_inverse_density_radial,&
                                              qs_upf_inverse_density_radial_adjoint
   USE upf_core,                        ONLY: upf_core_fourier
   USE upf_inverse,                     ONLY: &
        upf_inverse_basis, upf_inverse_basis_response, upf_inverse_fields_init, upf_inverse_init, &
        upf_inverse_operator, upf_inverse_operator_response, upf_inverse_project, &
        upf_inverse_project_response, upf_inverse_type
   USE upf_inverse_coulomb,             ONLY: &
        upf_inverse_coulomb_center_type, upf_inverse_coulomb_compensate, &
        upf_inverse_coulomb_compensation_adjoint, upf_inverse_coulomb_fourier_kernel, &
        upf_inverse_coulomb_pair_kernel, upf_inverse_coulomb_periodic, &
        upf_inverse_coulomb_potential_type, upf_inverse_coulomb_radial, &
        upf_inverse_coulomb_radial_response
   USE upf_inverse_quadrature,          ONLY: upf_inverse_cross_density,&
                                              upf_inverse_density_residual,&
                                              upf_inverse_projector_density,&
                                              upf_inverse_semilocal_pair,&
                                              upf_inverse_semilocal_residual,&
                                              upf_inverse_semilocal_response,&
                                              upf_inverse_smooth_density
   USE upf_radial_fields,               ONLY: upf_radial_fields_evaluate,&
                                              upf_radial_fields_type

   IMPLICIT NONE
   REAL(dp), PARAMETER :: step = 1.e-5_dp
   REAL(dp) :: max_error = 0.0_dp
   CALL check_inverse()
   CALL check_projector_density()
   CALL check_density_residual()
   CALL check_radial_coulomb()
   CALL check_pair_coulomb()
   CALL check_periodic_coulomb()
   CALL check_center_angular()
   CALL check_fourier_coulomb()
   CALL check_core_fourier()
   CALL check_nc_fields()
   CALL check_gaussian_fields()
   WRITE (*, '(A)') 'upf_inverse_unittest PASS'
CONTAINS
! **************************************************************************************************
!> \brief Native hard/soft and inverse angular legs with a deliberately nonidentity Gram matrix.
! **************************************************************************************************
   SUBROUTINE check_center_angular()
      INTEGER                                            :: ir, side
      REAL(dp) :: adjoint(6), dh(2, 2), dpnt(6), ds(2, 2), exact, hard(2, 2), hard_adjoint(2, 2), &
         potential(2, 2), radial(2, 2), response, samples(2), shift, soft(2, 2), &
         soft_adjoint(2, 2), values(6), wr(2)
      TYPE(grid_atom_type)                               :: grid
      TYPE(harmonics_atom_type)                          :: harmonics

      grid%nr = 2
      grid%ng_sphere = 3
      ALLOCATE (grid%wa(3), harmonics%slm(3, 2))
      grid%wa = [1.1_dp, 0.7_dp, 1.4_dp]
      harmonics%slm(:, 1) = 1.0_dp/SQRT(4.0_dp*pi)
      harmonics%slm(:, 2) = [0.3_dp, -0.4_dp, 0.2_dp]
      values = [0.03_dp, -0.01_dp, 0.02_dp, 0.04_dp, 0.01_dp, -0.02_dp]
      dpnt = [0.011_dp, -0.003_dp, 0.007_dp, -0.004_dp, 0.009_dp, 0.005_dp]
      hard(:, 1) = [0.14_dp, 0.09_dp]
      hard(:, 2) = [-0.03_dp, 0.02_dp]
      soft = 0.3_dp*hard
      dh = 0.017_dp*hard
      ds = -0.023_dp*soft
      wr = [0.21_dp, 0.37_dp]
      potential(:, 1) = [0.08_dp, -0.03_dp]
      potential(:, 2) = [0.04_dp, 0.07_dp]
      CALL qs_upf_inverse_density_radial(grid, harmonics, values, radial, hard, soft)
      CALL qs_upf_inverse_density_radial_adjoint(grid, harmonics, potential, adjoint, hard_adjoint, soft_adjoint)
      exact = 0.0_dp
      response = 0.0_dp
      DO ir = 1, 2
         exact = exact + wr(ir)*SUM(grid%wa*adjoint((ir - 1)*3 + 1:ir*3)*values((ir - 1)*3 + 1:ir*3)) + &
                 wr(ir)*SUM(hard_adjoint(ir, :)*hard(ir, :) + soft_adjoint(ir, :)*soft(ir, :))
         response = response + wr(ir)*SUM(grid%wa*adjoint((ir - 1)*3 + 1:ir*3)*dpnt((ir - 1)*3 + 1:ir*3)) + &
                    wr(ir)*SUM(hard_adjoint(ir, :)*dh(ir, :) + soft_adjoint(ir, :)*ds(ir, :))
      END DO
      CALL check(ABS(SUM(SPREAD(wr, 2, 2)*radial*potential) - exact), 'complete finite angular adjoint', 2.e-14_dp)
      DO side = 1, 2
         shift = REAL(2*side - 3, dp)*step
         CALL qs_upf_inverse_density_radial(grid, harmonics, values + shift*dpnt, radial, &
                                            hard + shift*dh, soft + shift*ds)
         samples(side) = SUM(SPREAD(wr, 2, 2)*radial*potential)
      END DO
      CALL check(ABS((samples(2) - samples(1))/(2.0_dp*step) - response), 'all angular density legs', 2.e-11_dp)
      CALL qs_upf_inverse_density_radial(grid, harmonics, values, radial)
      CALL check(ABS(SUM(SPREAD(wr, 2, 2)*radial*potential) - &
                     SUM(wr*[SUM(grid%wa*adjoint(1:3)*values(1:3)), SUM(grid%wa*adjoint(4:6)*values(4:6))])), &
                 'inverse-only angular adjoint', 2.e-14_dp)
      DEALLOCATE (grid%wa, harmonics%slm)
   END SUBROUTINE check_center_angular

! **************************************************************************************************
!> \brief Periodic residual image sum against the shell theorem, both density legs and cell response.
! **************************************************************************************************
   SUBROUTINE check_periodic_coulomb()
      COMPLEX(dp)                                        :: adjoint(3), background(3), db(3)
      INTEGER                                            :: a, b, i, ix, iy, iz, j, side
      REAL(dp) :: cell(3, 3), dcell(3, 3), dg(3, 3), displacement(3), distance, dkw(3), dvolume, &
         energy, exact, g(3, 3), kernel, kw(3), norm, qa(2), qb(2), response, samples(2), shift, &
         trace, volume
      TYPE(upf_inverse_coulomb_center_type)              :: centers(2), changed(2), derivative(2)
      TYPE(upf_inverse_coulomb_potential_type), &
         ALLOCATABLE                                     :: potential(:)

      cell(:, 1) = [2.0_dp, 0.0_dp, 0.0_dp]
      cell(:, 2) = [0.25_dp, 1.9_dp, 0.0_dp]
      cell(:, 3) = [0.1_dp, 0.2_dp, 2.3_dp]
      dcell = 0.013_dp*cell
      volume = 2.0_dp*1.9_dp*2.3_dp
      dvolume = 0.039_dp*volume
      g(:, 1) = 0.0_dp
      g(:, 2) = pi*[1.0_dp, -0.25_dp/1.9_dp, (0.25_dp*0.2_dp/1.9_dp - 0.1_dp)/2.3_dp]
      g(:, 3) = -g(:, 2)
      dg = -0.013_dp*g
      background = [CMPLX(0.0_dp, 0.0_dp, dp), CMPLX(0.012_dp, -0.007_dp, dp), CMPLX(0.012_dp, 0.007_dp, dp)]
      db = 0.021_dp*background
      kw = 0.0_dp
      DO a = 1, 2
         ALLOCATE (centers(a)%l(1), centers(a)%indices(1), centers(a)%radii(2), centers(a)%weights(2), &
                   centers(a)%density(2, 1), centers(a)%shapes(2, 1))
         centers(a)%l(:) = [0]
         centers(a)%indices(:) = [1]
         centers(a)%radii(:) = [0.5_dp, 1.1_dp] + REAL(a - 1, dp)*0.1_dp
         centers(a)%weights(:) = [0.2_dp, 0.4_dp]
         centers(a)%density(:, :) = RESHAPE([0.7_dp, -0.15_dp]/REAL(a, dp), [2, 1])
         centers(a)%shapes(:, :) = RESHAPE([0.8_dp, 0.2_dp], [2, 1])
         centers(a)%position = REAL(a - 1, dp)*[0.7_dp, 0.4_dp, -0.3_dp]
         derivative(a) = centers(a)
         derivative(a)%radii(:) = [0.011_dp, -0.007_dp]
         derivative(a)%weights(:) = [0.003_dp, -0.001_dp]
         derivative(a)%density(:, 1) = REAL(a, dp)*[0.023_dp, -0.017_dp]
         derivative(a)%shapes(:, :) = 0.017_dp*centers(a)%shapes
         derivative(a)%position = REAL(a, dp)*[0.009_dp, -0.011_dp, 0.013_dp]
      END DO
      CALL upf_inverse_coulomb_periodic(centers, cell, g, volume, background, kw, energy, potential, adjoint)
      exact = 0.0_dp
      DO a = 1, 2
         norm = SUM(centers(a)%weights*centers(a)%shapes(:, 1))
         qa = centers(a)%weights*(centers(a)%density(:, 1) - centers(a)%shapes(:, 1)* &
                                  SUM(centers(a)%weights*centers(a)%density(:, 1))/norm)
         DO b = 1, 2
            norm = SUM(centers(b)%weights*centers(b)%shapes(:, 1))
            qb = centers(b)%weights*(centers(b)%density(:, 1) - centers(b)%shapes(:, 1)* &
                                     SUM(centers(b)%weights*centers(b)%density(:, 1))/norm)
            DO iz = -3, 3
               DO iy = -3, 3
                  DO ix = -3, 3
                     displacement = centers(a)%position - centers(b)%position - MATMUL(cell, REAL([ix, iy, iz], dp))
                     distance = SQRT(SUM(displacement**2))
                     DO i = 1, 2
                        DO j = 1, 2
                           IF (distance <= ABS(centers(a)%radii(i) - centers(b)%radii(j))) THEN
                              kernel = 4.0_dp*pi/MAX(centers(a)%radii(i), centers(b)%radii(j))
                           ELSE IF (distance >= centers(a)%radii(i) + centers(b)%radii(j)) THEN
                              kernel = 4.0_dp*pi/distance
                           ELSE
                              kernel = pi/(distance*centers(a)%radii(i)*centers(b)%radii(j))* &
                                       (2.0_dp*distance*(centers(a)%radii(i) + centers(b)%radii(j)) - distance**2 - &
                                        (centers(a)%radii(i) - centers(b)%radii(j))**2)
                           END IF
                           exact = exact + 0.5_dp*qa(i)*qb(j)*kernel
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
      CALL check(ABS(energy - exact), "Complete periodic residual shell sum", 3.e-12_dp)
      kw(2:3) = volume*4.0_dp*pi/SUM(g(:, 2)**2)
      dkw = 0.065_dp*kw
      CALL upf_inverse_coulomb_periodic(centers, cell, g, volume, background, kw, energy, potential, adjoint, &
                                        response=response, center_response=derivative, cell_response=dcell, &
                                        gvecs_response=dg, volume_response=dvolume, background_response=db, &
                                        poisson_weights_response=dkw)
      trace = REAL(SUM(CONJG(background)*adjoint), dp)
      DO a = 1, 2
         trace = trace + SUM(SPREAD(centers(a)%weights, 2, 1)*centers(a)%density*potential(a)%values)
      END DO
      CALL check(ABS(trace - 2.0_dp*energy), "Periodic Coulomb density and background adjoints", 3.e-12_dp)
      trace = REAL(SUM(CONJG(db)*adjoint), dp)
      DO a = 1, 2
         trace = trace + SUM(SPREAD(centers(a)%weights, 2, 1)*derivative(a)%density*potential(a)%values)
      END DO
      DO side = 1, 2
         shift = REAL(2*side - 3, dp)*step
         changed = centers
         DO a = 1, 2
            changed(a)%density(:, :) = centers(a)%density + shift*derivative(a)%density
         END DO
         CALL upf_inverse_coulomb_periodic(changed, cell, g, volume, background + shift*db, kw, &
                                           samples(side), potential, adjoint)
      END DO
      CALL check(ABS((samples(2) - samples(1))/(2.0_dp*step) - trace), "Periodic Coulomb density variation", 3.e-10_dp)
      DO side = 1, 2
         shift = REAL(2*side - 3, dp)*step
         changed = centers
         DO a = 1, 2
            changed(a)%position = centers(a)%position + shift*derivative(a)%position
            changed(a)%radii(:) = centers(a)%radii + shift*derivative(a)%radii
            changed(a)%weights(:) = centers(a)%weights + shift*derivative(a)%weights
            changed(a)%density(:, :) = centers(a)%density + shift*derivative(a)%density
            changed(a)%shapes(:, :) = centers(a)%shapes + shift*derivative(a)%shapes
         END DO
         CALL upf_inverse_coulomb_periodic(changed, cell + shift*dcell, g + shift*dg, volume + shift*dvolume, &
                                           background + shift*db, kw + shift*dkw, samples(side), potential, adjoint)
      END DO
      CALL check(ABS((samples(2) - samples(1))/(2.0_dp*step) - response), "Complete periodic Coulomb joint response", 3.e-9_dp)
   END SUBROUTINE check_periodic_coulomb
! **************************************************************************************************
!> \brief Translated shell theorem, overlap, reciprocity and joint geometry response.
! **************************************************************************************************
   SUBROUTINE check_pair_coulomb()
      INTEGER, PARAMETER                                 :: indices(4) = [1, 2, 3, 4], &
                                                            ls(4) = [0, 1, 1, 1]

      INTEGER                                            :: i, j, side
      REAL(dp) :: a(2), b(2), d(3), da(2), db(2), dd(3), distance, exact, kernel(2, 4, 2, 4), &
         response(2, 4, 2, 4), reverse(2, 4, 2, 4), samples(2, 4, 2, 4, 2), shift

      a = [0.7_dp, 1.3_dp]
      b = [0.5_dp, 1.1_dp]
      da = [0.013_dp, -0.007_dp]
      db = [-0.011_dp, 0.017_dp]
      dd = [0.019_dp, -0.023_dp, 0.009_dp]
      DO j = 1, 3
         d = REAL(j, dp)*[0.6_dp, -0.4_dp, 0.3_dp]
         distance = SQRT(SUM(d**2))
         CALL upf_inverse_coulomb_pair_kernel(a, ls, indices, b, ls, indices, d, &
                                              kernel, response, da, db, dd)
         CALL upf_inverse_coulomb_pair_kernel(b, ls, indices, a, ls, indices, -d, reverse)
         DO i = 1, 2
            DO side = 1, 2
               IF (distance >= a(i) + b(side)) THEN
                  exact = 4.0_dp*pi/distance
               ELSE IF (distance <= ABS(a(i) - b(side))) THEN
                  exact = 4.0_dp*pi/MAX(a(i), b(side))
               ELSE
                  exact = pi/(distance*a(i)*b(side))* &
                          (2.0_dp*distance*(a(i) + b(side)) - distance**2 - (a(i) - b(side))**2)
               END IF
               CALL check(ABS(kernel(i, 1, side, 1) - exact), "Translated monopole shells", 3.e-12_dp)
               CALL check(MAXVAL(ABS(kernel(i, :, side, :) - TRANSPOSE(reverse(side, :, i, :)))), &
                          "Translated Coulomb reciprocity", 1.e-12_dp)
            END DO
         END DO
         DO side = 1, 2
            shift = REAL(2*side - 3, dp)*step
            CALL upf_inverse_coulomb_pair_kernel(a + shift*da, ls, indices, b + shift*db, ls, indices, &
                                                 d + shift*dd, samples(:, :, :, :, side))
         END DO
         CALL check(MAXVAL(ABS((samples(:, :, :, :, 2) - samples(:, :, :, :, 1))/(2.0_dp*step) - response)), &
                    "Translated Coulomb joint response", 3.e-9_dp)
      END DO
      ! A source sphere passing through a target angular node is integrated in separate polar caps.
      d = [0.7_dp, 0.0_dp, 0.0_dp]
      CALL upf_inverse_coulomb_pair_kernel(a, ls, indices, b, ls, indices, d, &
                                           kernel, response, da, db, dd)
      DO side = 1, 2
         shift = REAL(2*side - 3, dp)*step
         CALL upf_inverse_coulomb_pair_kernel(a + shift*da, ls, indices, b + shift*db, ls, indices, &
                                              d + shift*dd, samples(:, :, :, :, side))
      END DO
      CALL check(MAXVAL(ABS((samples(:, :, :, :, 2) - samples(:, :, :, :, 1))/(2.0_dp*step) - response)), &
                 "Translated shell crossing response", 3.e-9_dp)
      d = [0.0_dp, 0.0_dp, a(2) + b(2)]
      CALL upf_inverse_coulomb_pair_kernel(a, ls, indices, b, ls, indices, d, &
                                           kernel, response, da, db, dd)
      DO side = 1, 2
         shift = REAL(2*side - 3, dp)*step
         CALL upf_inverse_coulomb_pair_kernel(a + shift*da, ls, indices, b + shift*db, ls, indices, &
                                              d + shift*dd, samples(:, :, :, :, side))
      END DO
      CALL check(MAXVAL(ABS((samples(:, :, :, :, 2) - samples(:, :, :, :, 1))/(2.0_dp*step) - response)), &
                 "Translated tangent shell response", 3.e-9_dp)
      d = 0.0_dp
      CALL upf_inverse_coulomb_pair_kernel(a, ls, indices, b, ls, indices, d, kernel)
      CALL check(ABS(kernel(1, 2, 1, 2) - 4.0_dp*pi/3.0_dp*b(1)/a(1)**2), "Concentric dipole shells")
      CALL check(ABS(kernel(1, 2, 1, 3)), "Concentric harmonic orthogonality")
   END SUBROUTINE check_pair_coulomb
! **************************************************************************************************
!> \brief Original core normalization, finite endpoint and joint phase/lattice/volume derivative.
! **************************************************************************************************
   SUBROUTINE check_core_fourier()
      COMPLEX(dp)                                        :: dg(3, 3), dr(3, 3), response(3), rho(3), &
                                                            samples(3, 2)
      INTEGER                                            :: side
      REAL(dp)                                           :: delta_g(3, 3), delta_r(3), g(3, 3), &
                                                            position(3), shift
      TYPE(atom_upfpot_type)                             :: pot

      ALLOCATE (pot%r(5), pot%rab(5), pot%rho_nlcc(5))
      pot%r(:) = [0.0_dp, 0.5_dp, 1.0_dp, 1.5_dp, 2.0_dp]
      pot%rab = 0.5_dp
      pot%rho_nlcc(:) = [1.0_dp, 0.8_dp, 0.5_dp, 0.2_dp, 0.1_dp]
      pot%core_correction = .TRUE.
      g(:, 1) = 0.0_dp
      g(:, 2) = [0.3_dp, -0.2_dp, 0.4_dp]
      g(:, 3) = -g(:, 2)
      position = [0.2_dp, -0.3_dp, 0.1_dp]
      delta_r = [0.02_dp, 0.01_dp, -0.03_dp]
      delta_g = 0.02_dp*g
      CALL upf_core_fourier(pot, g, position, 3.0_dp, rho, dr, dg)
      CALL check(ABS(3.0_dp*rho(1) - 4.0_dp*pi/6.0_dp* &
                     SUM([1.0_dp, 4.0_dp, 2.0_dp, 4.0_dp, 1.0_dp]*pot%r**2*pot%rho_nlcc)), "Original core zero mode", 1.e-12_dp)
      CALL check(ABS(rho(3) - CONJG(rho(2))), "Original core Fourier reality", 1.e-12_dp)
      response = SUM(dr*SPREAD(delta_r, 2, 3), DIM=1) + SUM(dg*delta_g, DIM=1) - 0.01_dp*rho
      DO side = 1, 2
         shift = REAL(2*side - 3, dp)*step
         CALL upf_core_fourier(pot, g + shift*delta_g, position + shift*delta_r, &
                               3.0_dp*(1.0_dp + 0.01_dp*shift), samples(:, side))
      END DO
      CALL check(MAXVAL(ABS((samples(:, 2) - samples(:, 1))/(2.0_dp*step) - response)), "Original core joint response", 2.e-9_dp)
      pot%core_correction = .FALSE.
      CALL upf_core_fourier(pot, g, position, 3.0_dp, rho, dr, dg)
      CALL check(MAXVAL(ABS(rho)) + MAXVAL(ABS(dr)) + MAXVAL(ABS(dg)), "Absent partial core", 0.0_dp)
   END SUBROUTINE check_core_fourier

! **************************************************************************************************
!> \brief Signed shell/PW coupling, reciprocal zero limit and joint geometry response.
! **************************************************************************************************
   SUBROUTINE check_fourier_coulomb()
      INTEGER, PARAMETER                                 :: ng = 4, nh = 3, nr = 2

      COMPLEX(dp)                                        :: derivative(ng, nr, nh), &
                                                            kernel(ng, nr, nh), other(ng, nr, nh), &
                                                            samples(ng, nr, nh, 2)
      INTEGER                                            :: l(nh), lm(nh), side
      REAL(dp)                                           :: dg(3, ng), dp_center(3), dr(nr), &
                                                            g(3, ng), position(3), radii(nr), &
                                                            shift, volume

      radii = [0.3_dp, 1.2_dp]
      l = [0, 1, 2]
      lm = [1, 3, 7]
      g(:, 1) = 0.0_dp
      g(:, 2) = [0.7_dp, 0.0_dp, 0.0_dp]
      g(:, 3) = [0.3_dp, -0.2_dp, 0.4_dp]
      g(:, 4) = -g(:, 3)
      position = [0.2_dp, -0.7_dp, 0.4_dp]
      volume = 7.0_dp
      dr = [0.07_dp, -0.03_dp]
      dg = 0.13_dp - 0.1_dp*g
      dp_center = [-0.04_dp, 0.08_dp, 0.02_dp]
      CALL upf_inverse_coulomb_fourier_kernel(radii, l, lm, g, position, volume, kernel, derivative, &
                                              dr, dg, dp_center, 0.3_dp)
      CALL check(MAXVAL(ABS(kernel(1, :, 1) - SQRT(4*pi)/volume)), 'shell monopole Fourier normalization')
      CALL check(MAXVAL(ABS(kernel(1, :, 2:))), 'zero higher shell multipoles at G=0', 0.0_dp)
      CALL check(MAXVAL(ABS(kernel(2, :, 1) - SQRT(4*pi)/volume* &
                            EXP(CMPLX(0.0_dp, -DOT_PRODUCT(g(:, 2), position), dp))* &
                            SIN(0.7_dp*radii)/(0.7_dp*radii))), 'analytic spherical shell form factor')
      CALL check(MAXVAL(ABS(kernel(3, :, :) - CONJG(kernel(4, :, :)))), 'real density reciprocal conjugacy')
      CALL check(MAXVAL(ABS(derivative(1, :, 2) - CMPLX(0.0_dp, &
                                               -SQRT(12*pi)/volume*radii*dg(3, 1)/3.0_dp, dp))), 'dipole Fourier derivative at G=0')
      DO side = 1, 2
         shift = step*(2*side - 3)
         CALL upf_inverse_coulomb_fourier_kernel(radii + shift*dr, l, lm, g + shift*dg, &
                                                 position + shift*dp_center, volume + shift*0.3_dp, samples(:, :, :, side))
      END DO
      CALL check(MAXVAL(ABS((samples(:, :, :, 2) - samples(:, :, :, 1))/(2*step) - derivative)), &
                 'shell Fourier radius center reciprocal volume response', 2.e-9_dp)
      CALL upf_inverse_coulomb_fourier_kernel(radii, l, lm, g, position, volume, other)
      CALL check(MAXVAL(ABS(kernel - other)), 'optional Fourier response preserves forward map', 0.0_dp)
   END SUBROUTINE check_fourier_coulomb

! **************************************************************************************************
!> \brief Signed Coulomb reciprocity, shell self term, moment compensation and energy adjoints.
! **************************************************************************************************
   SUBROUTINE check_radial_coulomb()
      INTEGER, PARAMETER                                 :: nh = 3, nr = 3

      INTEGER                                            :: channel, l(nh), side
      REAL(dp) :: adjoint(nr, nh), compensation(nr, nh), density(nr, nh), dr(nr), dw(nr), energy, &
         energy_other, moments(nh), one_density(1, 1), one_potential(1, 1), other(nr, nh), &
         potential(nr, nh), radii(nr), residual(nr, nh), response, samples(2), shapes(nr, nh), &
         shift, vp(nr, nh), weights(nr)

      ! An orthonormal monopole on one radial shell has E=2*pi*(w*rho)**2/r.
      one_density = 0.7_dp
      CALL upf_inverse_coulomb_radial([0.8_dp], [0.3_dp], [0], one_density, one_potential, energy)
      CALL check(ABS(energy - 2*pi*(0.3_dp*0.7_dp)**2/0.8_dp), 'Coulomb shell self energy')
      CALL upf_inverse_coulomb_radial([0.8_dp], [0.3_dp], [2], one_density, one_potential, energy)
      CALL check(ABS(energy - 2*pi*(0.3_dp*0.7_dp)**2/(5*0.8_dp)), 'Coulomb higher multipole shell')
      ! Unordered radii and signed, repeated-l channels are intentional.
      radii = [1.7_dp, 0.3_dp, 0.8_dp]
      weights = [0.4_dp, 0.2_dp, 0.7_dp]
      l = [0, 2, 2]
      density = RESHAPE([0.1_dp, -0.3_dp, 0.2_dp, -0.2_dp, 0.7_dp, -0.1_dp, &
                         0.3_dp, -0.2_dp, 0.4_dp], [nr, nh])
      other = 0.13_dp + 0.2_dp*density
      dr = [0.02_dp, -0.01_dp, 0.03_dp]
      dw = 0.03_dp*weights
      CALL upf_inverse_coulomb_radial(radii, weights, l, density, potential, energy)
      CALL upf_inverse_coulomb_radial(radii, weights, l, other, vp, energy_other)
      CALL check(ABS(SUM(SPREAD(weights, 2, nh)*(density*vp - other*potential))), 'signed Coulomb reciprocity')
      CALL check(MAX(0.0_dp, -energy), 'signed Coulomb nonnegative self energy', 0.0_dp)
      CALL upf_inverse_coulomb_radial_response(radii, weights, l, density, dr, dw, other, response)
      DO side = 1, 2
         shift = step*(2*side - 3)
         CALL upf_inverse_coulomb_radial(radii + shift*dr, weights + shift*dw, l, &
                                         density + shift*other, vp, samples(side))
      END DO
      CALL check(ABS((samples(2) - samples(1))/(2*step) - response), 'Coulomb density measure kernel response', 2.e-9_dp)
      DO channel = 1, nh
         shapes(:, channel) = radii**l(channel)*EXP(-0.6_dp*radii**2)
      END DO
      CALL upf_inverse_coulomb_compensate(radii, weights, l, density, shapes, compensation, residual, moments)
      CALL check(MAXVAL(ABS(density - compensation - residual)), 'Coulomb compensation density closure')
      DO channel = 1, nh
         CALL check(ABS(SUM(weights*radii**l(channel)*residual(:, channel))), 'zero compensated exterior moment')
         CALL check(ABS(moments(channel) - SUM(weights*radii**l(channel)*compensation(:, channel))), &
                    'original signed exterior moment retained')
      END DO
      CALL upf_inverse_coulomb_radial(radii, weights, l, residual, potential, energy)
      CALL upf_inverse_coulomb_compensation_adjoint(radii, weights, l, shapes, potential, adjoint)
      response = SUM(SPREAD(weights, 2, nh)*other*adjoint)
      DO side = 1, 2
         shift = step*(2*side - 3)
         CALL upf_inverse_coulomb_compensate(radii, weights, l, density + shift*other, shapes, &
                                             compensation, residual, moments)
         CALL upf_inverse_coulomb_radial(radii, weights, l, residual, vp, samples(side))
      END DO
      CALL check(ABS((samples(2) - samples(1))/(2*step) - response), 'occupation-dependent compensation adjoint', 2.e-9_dp)
      density = 0.0_dp
      CALL upf_inverse_coulomb_radial(radii, weights, l, density, potential, energy)
      CALL check(MAXVAL(ABS(potential)) + ABS(energy), 'exact zero NC Coulomb residual', 0.0_dp)
   END SUBROUTINE check_radial_coulomb

! **************************************************************************************************
!> \brief Gaussian/residual closure, adjoint and joint field/measure response on complex rows.
! **************************************************************************************************
   SUBROUTINE check_density_residual()
      INTEGER, PARAMETER                                 :: nb = 2, np = 3

      COMPLEX(dp) :: c(nb, nb), center_difference(np, nb), center_gradient(np, nb, 3), &
         dgr(np, nb, 3), dgs(np, nb, 3), dr(np, nb), ds(np, nb), gr(np, nb, 3), gs(np, nb, 3), &
         h(nb, nb), hr(nb, nb), hs(nb, nb), p(nb, nb), reference(np, nb), response(nb, nb), &
         soft(np, nb)
      INTEGER                                            :: axis, center, i, j, side
      REAL(dp) :: dw(np), gradient(np, 3), gradient_r(np, 3), gradient_s(np, 3), rho(np), &
         rho_r(np), rho_s(np), samples(2), shift, tau(np), tau_r(np), tau_s(np), value, vg(np, 3), &
         vn(np), vt(np), weights(np)

      c(:, 1) = [CMPLX(1.0_dp, 0.2_dp, dp), CMPLX(0.1_dp, -0.3_dp, dp)]
      c(:, 2) = [CMPLX(-0.4_dp, 0.3_dp, dp), CMPLX(0.2_dp, 0.5_dp, dp)]
      p = MATMUL(c, CONJG(TRANSPOSE(c)))
      DO j = 1, nb
         DO i = 1, np
            reference(i, j) = CMPLX(SIN(REAL(i + j, dp)), COS(REAL(i*j, dp)), dp)
            soft(i, j) = reference(i, j) + 0.1_dp*CONJG(reference(i, j))
            dr(i, j) = 0.07_dp*reference(i, j)
            ds(i, j) = 0.05_dp*soft(i, j)
            DO axis = 1, 3
               gr(i, j, axis) = CMPLX(COS(REAL(i + j + axis, dp)), SIN(REAL(i*j + axis, dp)), dp)
               gs(i, j, axis) = gr(i, j, axis) + 0.13_dp*CONJG(gr(i, j, axis))
               dgr(i, j, axis) = 0.03_dp*gr(i, j, axis)
               dgs(i, j, axis) = 0.09_dp*gs(i, j, axis)
            END DO
         END DO
      END DO
      weights = [0.3_dp, 0.5_dp, 0.7_dp]
      dw = 0.04_dp*weights
      vn = [0.2_dp, 0.7_dp, -0.3_dp]
      vt = [0.3_dp, 0.1_dp, 0.2_dp]
      vg = RESHAPE([0.1_dp, -0.2_dp, 0.3_dp, 0.4_dp, 0.5_dp, -0.6_dp, 0.7_dp, -0.8_dp, 0.9_dp], [np, 3])
      CALL upf_inverse_smooth_density(reference, gr, p, rho_r, gradient_r, tau_r)
      CALL upf_inverse_smooth_density(soft, gs, p, rho_s, gradient_s, tau_s)
      CALL upf_inverse_density_residual(reference, gr, soft, gs, p, rho, gradient, tau)
      CALL check(MAXVAL(ABS(rho_r + rho - rho_s)), 'Gaussian/residual density closure')
      CALL check(MAXVAL(ABS(gradient_r + gradient - gradient_s)), 'Gaussian/residual gradient closure')
      CALL check(MAXVAL(ABS(tau_r + tau - tau_s)), 'Gaussian/residual kinetic closure')
      CALL upf_inverse_semilocal_pair(reference, reference, gr, gr, weights, vn, vg, vt, hr)
      CALL upf_inverse_semilocal_pair(soft, soft, gs, gs, weights, vn, vg, vt, hs)
      CALL upf_inverse_semilocal_residual(reference, gr, soft, gs, weights, vn, vg, vt, h)
      CALL check(MAXVAL(ABS(hr + h - hs)), 'Gaussian/residual operator closure')
      value = SUM(weights*(vn*rho + SUM(vg*gradient, DIM=2) + vt*tau))
      CALL check(ABS(value - REAL(SUM(h*CONJG(p)), dp)), 'signed density operator adjoint')
      ! Overlapping center fields carry opposite off-diagonal terms that must not be discarded.
      rho_s = 0.0_dp
      gradient_s = 0.0_dp
      tau_s = 0.0_dp
      hs = CMPLX(0.0_dp, 0.0_dp, dp)
      DO center = 1, 2
         center_difference = (soft - reference)/2 + (2*center - 3)*0.2_dp*CONJG(reference)
         center_gradient = (gs - gr)/2 + (2*center - 3)*0.2_dp*CONJG(gr)
         CALL upf_inverse_cross_density(center_difference, soft + reference, center_gradient, gs + gr, &
                                        p, rho_r, gradient_r, tau_r)
         CALL upf_inverse_semilocal_pair(center_difference, soft + reference, center_gradient, gs + gr, &
                                         weights, vn, vg, vt, hr)
         hr = 0.5_dp*(hr + CONJG(TRANSPOSE(hr)))
         value = SUM(weights*(vn*rho_r + SUM(vg*gradient_r, DIM=2) + vt*tau_r))
         CALL check(ABS(value - REAL(SUM(hr*CONJG(p)), dp)), 'center residual operator adjoint')
         rho_s = rho_s + rho_r
         gradient_s = gradient_s + gradient_r
         tau_s = tau_s + tau_r
         hs = hs + hr
      END DO
      CALL check(MAXVAL(ABS(rho_s - rho)), 'overlapping center charge closure')
      CALL check(MAXVAL(ABS(gradient_s - gradient)), 'overlapping center gradient closure')
      CALL check(MAXVAL(ABS(tau_s - tau)), 'overlapping center kinetic closure')
      CALL check(MAXVAL(ABS(hs - h)), 'overlapping center operator closure')
      CALL upf_inverse_semilocal_response(soft - reference, soft + reference, gs - gr, gs + gr, &
                                          weights, vn, vg, vt, ds - dr, ds + dr, dgs - dgr, dgs + dgr, &
                                          dw, 0.05_dp*vn, 0.07_dp*vg, 0.09_dp*vt, response)
      DO side = 1, 2
         shift = step*(2*side - 3)
         CALL upf_inverse_density_residual(reference + shift*dr, gr + shift*dgr, soft + shift*ds, gs + shift*dgs, &
                                           p, rho, gradient, tau)
         samples(side) = SUM((weights + shift*dw)*((1.0_dp + 0.05_dp*shift)*vn*rho + &
                                                   (1.0_dp + 0.07_dp*shift)*SUM(vg*gradient, DIM=2) + &
                                                   (1.0_dp + 0.09_dp*shift)*vt*tau))
      END DO
      CALL check(ABS((samples(2) - samples(1))/(2*step) - REAL(SUM(response*CONJG(p)), dp)), &
                 'signed density geometry adjoint', 2.e-9_dp)
      CALL upf_inverse_density_residual(reference, gr, reference, gr, p, rho, gradient, tau)
      CALL check(MAXVAL(ABS(rho)) + MAXVAL(ABS(gradient)) + MAXVAL(ABS(tau)), 'exact zero NC residual', 0.0_dp)
      CALL upf_inverse_semilocal_residual(reference(:0, :), gr(:0, :, :), soft(:0, :), gs(:0, :, :), &
                                          weights(:0), vn(:0), vg(:0, :), vt(:0), h)
      CALL check(MAXVAL(ABS(h)), 'empty-row residual adjoint', 0.0_dp)
   END SUBROUTINE check_density_residual

! **************************************************************************************************
!> \brief Compare density congruence to independently projected occupied orbitals.
! **************************************************************************************************
   SUBROUTINE check_projector_density()
      COMPLEX(dp)                                        :: c(2, 2), density(2, 2), projected(3, 3), &
                                                            q(3, 3), reference(3, 3), w(3, 2), &
                                                            wc(3, 2)
      INTEGER                                            :: i, j

      DO j = 1, 2
         DO i = 1, 2
            c(i, j) = CMPLX(SIN(REAL(i + j, dp)), COS(REAL(i*j, dp)), dp)
         END DO
         DO i = 1, 3
            w(i, j) = 60.0_dp*CMPLX(COS(REAL(i + j, dp)), SIN(REAL(i*j, dp)), dp)
         END DO
      END DO
      density = MATMUL(c, CONJG(TRANSPOSE(c)))
      CALL upf_inverse_projector_density(w, density, projected)
      wc = MATMUL(w, c)
      reference = MATMUL(wc, CONJG(TRANSPOSE(wc)))
      CALL check(MAXVAL(ABS(projected - reference))/MAXVAL(ABS(reference)), 'projected occupied orbitals')
      CALL check(MAXVAL(ABS(projected - CONJG(TRANSPOSE(projected)))), 'Hermitian projector density', 0.0_dp)
      q = reference/MAXVAL(ABS(reference))
      CALL check(ABS(SUM(projected*CONJG(q)) - &
                     SUM(density*CONJG(MATMUL(CONJG(TRANSPOSE(w)), MATMUL(q, w)))))/MAXVAL(ABS(reference)), &
                 'projector density operator adjoint')
   END SUBROUTINE check_projector_density

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE check_inverse()
      INTEGER, PARAMETER                                 :: nb = 2, nf = 5, np = 3

      COMPLEX(KIND=dp) :: actual(nb, nb), b(np, nb), bd(np, np), chi(nf, nb), dchi(nf, nb), &
         dd(np, np), ddd(np, np), ddelta(nf, np), ddg(np, nb), delta(nf, np), dg(np, nb), &
         dgd(nb, np), dgg(nb, nb), donsite(np, np), dop(nf, nf), dprojector(np, nf), ds(nf, nb), &
         dw(np, nb), empty_b(0, nb), empty_delta(nf, 0), empty_q(0, 0), empty_w(0, nb), &
         gd(nb, np), gg(nb, nb), images(nf, nb, 2), onsite(np, np), op(nf, nf), projected(np, nb), &
         projector(np, nf), reference(nb, nb), response(nb, nb), sampled(nb, nb, 2), &
         shifted_chi(nf, nb), shifted_delta(nf, np), shifted_op(nf, nf), shifted_p(np, nf)
      COMPLEX(KIND=dp) :: shifted_q(np, np), soft(nf, nb), w(np, nb)
      INTEGER                                            :: i, info, j, phase, side
      REAL(KIND=dp)                                      :: shift
      TYPE(upf_inverse_type)                             :: inverse

      DO phase = 1, 3
         DO j = 1, nb
            DO i = 1, nf
               chi(i, j) = CMPLX(SIN(REAL(i + 2*j, dp)), (phase - 1)*COS(REAL(i*j, dp)), dp)
               dchi(i, j) = 0.07_dp*CMPLX(COS(REAL(i + j, dp)), (phase - 1)*SIN(REAL(2*i + j, dp)), dp)
            END DO
         END DO
         DO j = 1, np
            DO i = 1, nf
               delta(i, j) = 0.12_dp*CMPLX(COS(REAL(i*j, dp)), (phase - 1)*SIN(REAL(i + j, dp)), dp)
               ddelta(i, j) = 0.03_dp*CMPLX(SIN(REAL(i + 2*j, dp)), (phase - 1)*COS(REAL(i*j, dp)), dp)
               projector(j, i) = 0.15_dp*CMPLX(SIN(REAL(i + j, dp)), (phase - 1)*COS(REAL(2*i + j, dp)), dp)
               dprojector(j, i) = 0.02_dp*CMPLX(COS(REAL(i + j, dp)), (phase - 1)*SIN(REAL(i*j, dp)), dp)
            END DO
         END DO
         IF (phase == 3) THEN
            ! A mixed NC/USPP frame retains the NC bra, with no NC reconstruction.
            delta(:, np) = CMPLX(0.0_dp, 0.0_dp, dp)
            ddelta(:, np) = CMPLX(0.0_dp, 0.0_dp, dp)
         END IF
         ! A non-Hermitian operator exposes an accidental conjugation/symmetrization.
         DO j = 1, nf
            DO i = 1, nf
               op(i, j) = CMPLX(0.1_dp*i - 0.05_dp*j, (phase - 1)*0.03_dp*(i + j), dp)
               dop(i, j) = CMPLX(0.02_dp*(i + j), (phase - 1)*0.01_dp*(i - j), dp)
            END DO
         END DO
         onsite = op(1:np, 1:np)*0.3_dp
         donsite = dop(1:np, 1:np)*0.4_dp
         bd = MATMUL(projector, delta)
         b = MATMUL(projector, chi)
         CALL upf_inverse_init(bd, inverse, info)
         IF (info /= 0 .OR. .NOT. inverse%ready) ERROR STOP 'Inverse factorization failed'
         CALL upf_inverse_project(inverse, b, w)
         CALL upf_inverse_basis(chi, delta, w, soft)
         projected = MATMUL(projector, soft)
         CALL check(MAXVAL(ABS(projected - w)), 'inverse projector identity')
         CALL check(MAXVAL(ABS(soft + MATMUL(delta, projected) - chi)), 'full reconstruction round trip')
         CALL upf_inverse_project_response(inverse, w, MATMUL(dprojector, chi) + MATMUL(projector, dchi), &
                                           MATMUL(dprojector, delta) + MATMUL(projector, ddelta), dw)
         CALL upf_inverse_basis_response(dchi, delta, ddelta, w, dw, ds)
         gg = MATMUL(CONJG(TRANSPOSE(chi)), MATMUL(op, chi))
         gd = MATMUL(CONJG(TRANSPOSE(chi)), MATMUL(op, delta))
         dg = MATMUL(CONJG(TRANSPOSE(delta)), MATMUL(op, chi))
         dd = MATMUL(CONJG(TRANSPOSE(delta)), MATMUL(op, delta))
         CALL upf_inverse_operator(gg, gd, dg, dd, onsite, w, actual)
         reference = MATMUL(CONJG(TRANSPOSE(soft)), MATMUL(op, soft)) + &
                     MATMUL(CONJG(TRANSPOSE(projected)), MATMUL(onsite, projected))
         CALL check(MAXVAL(ABS(actual - reference)), 'inverse operator versus direct fields')
         dgg = MATMUL(CONJG(TRANSPOSE(dchi)), MATMUL(op, chi)) + &
               MATMUL(CONJG(TRANSPOSE(chi)), MATMUL(dop, chi) + MATMUL(op, dchi))
         dgd = MATMUL(CONJG(TRANSPOSE(dchi)), MATMUL(op, delta)) + &
               MATMUL(CONJG(TRANSPOSE(chi)), MATMUL(dop, delta) + MATMUL(op, ddelta))
         ddg = MATMUL(CONJG(TRANSPOSE(ddelta)), MATMUL(op, chi)) + &
               MATMUL(CONJG(TRANSPOSE(delta)), MATMUL(dop, chi) + MATMUL(op, dchi))
         ddd = MATMUL(CONJG(TRANSPOSE(ddelta)), MATMUL(op, delta)) + &
               MATMUL(CONJG(TRANSPOSE(delta)), MATMUL(dop, delta) + MATMUL(op, ddelta))
         CALL upf_inverse_operator_response(dgg, gd, dg, dd, onsite, w, dgd, ddg, ddd, donsite, dw, response)
         DO side = 1, 2
            shift = (2*side - 3)*step
            shifted_chi = chi + shift*dchi
            shifted_delta = delta + shift*ddelta
            shifted_p = projector + shift*dprojector
            shifted_op = op + shift*dop
            shifted_q = onsite + shift*donsite
            CALL upf_inverse_init(MATMUL(shifted_p, shifted_delta), inverse, info)
            IF (info /= 0) ERROR STOP 'Displaced inverse factorization failed'
            CALL upf_inverse_project(inverse, MATMUL(shifted_p, shifted_chi), w)
            CALL upf_inverse_basis(shifted_chi, shifted_delta, w, images(:, :, side))
            projected = MATMUL(shifted_p, images(:, :, side))
            sampled(:, :, side) = MATMUL(CONJG(TRANSPOSE(images(:, :, side))), &
                                         MATMUL(shifted_op, images(:, :, side))) + &
                                  MATMUL(CONJG(TRANSPOSE(projected)), MATMUL(shifted_q, projected))
         END DO
         CALL check(MAXVAL(ABS(ds - (images(:, :, 2) - images(:, :, 1))/(2*step))), 'inverse basis response', 1.e-9_dp)
         CALL check(MAXVAL(ABS(response - (sampled(:, :, 2) - sampled(:, :, 1))/(2*step))), &
                    'inverse operator response', 1.e-9_dp)
      END DO

      CALL upf_inverse_init(empty_q, inverse, info)
      IF (info /= 0 .OR. .NOT. inverse%ready) ERROR STOP 'Empty inverse failed'
      CALL upf_inverse_project(inverse, empty_b, empty_w)
      CALL upf_inverse_basis(chi, empty_delta, empty_w, soft)
      CALL check(MAXVAL(ABS(soft - chi)), 'identity without reconstruction projectors')
      bd = CMPLX(0.0_dp, 0.0_dp, dp)
      bd(2, 2) = -1.0_dp
      CALL upf_inverse_init(bd, inverse, info)
      IF (info /= 2 .OR. inverse%ready) ERROR STOP 'Singular inverse accepted'
      bd(2, 2) = -1.0_dp + 1.e-12_dp
      CALL upf_inverse_init(bd, inverse, info, minimum_rcond=1.e-10_dp)
      IF (info /= np + 1 .OR. inverse%ready) ERROR STOP 'Ill-conditioned inverse accepted'
      bd(1, 1) = CMPLX(ieee_value(0.0_dp, ieee_quiet_nan), 0.0_dp, dp)
      CALL upf_inverse_init(bd, inverse, info)
      IF (info /= np + 2 .OR. inverse%ready) ERROR STOP 'Nonfinite inverse accepted'
   END SUBROUTINE check_inverse

! **************************************************************************************************
!> \brief NC original projectors need neither reference waves nor augmentation data.
! **************************************************************************************************
   SUBROUTINE check_nc_fields()
      REAL(dp)                                           :: displacement(3), gradient(3, 1), value(1)
      TYPE(atom_upfpot_type)                             :: pot
      TYPE(upf_radial_fields_type)                       :: beta, delta

      pot%pseudo_type = "NC"
      ALLOCATE (pot%r(7), pot%beta(7, 1), pot%lbeta(1))
      pot%r(:) = [0.0_dp, 0.2_dp, 0.5_dp, 1.0_dp, 1.5_dp, 2.0_dp, 3.0_dp]
      pot%lbeta = 0
      pot%beta(:, 1) = EXP(-pot%r**2)
      pot%beta(6:, 1) = 0.0_dp
      CALL upf_inverse_fields_init(pot, beta, delta)
      displacement = [0.5_dp, 0.0_dp, 0.0_dp]
      CALL upf_radial_fields_evaluate(beta, displacement, value)
      CALL check(ABS(value(1) - EXP(-0.25_dp)/SQRT(4.0_dp*pi)), 'original NC projector')
      CALL upf_radial_fields_evaluate(delta, displacement, value, gradient)
      CALL check(MAXVAL(ABS(value)) + MAXVAL(ABS(gradient)), 'zero NC reconstruction')
   END SUBROUTINE check_nc_fields

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE check_gaussian_fields()
      INTEGER, PARAMETER                                 :: nb = 11, nf = 5

      COMPLEX(dp)                                        :: dg(nf, nb, 3), dv(nf, nb), g(nf, nb, 3), &
                                                            gss(nf, nb, 3, 2), other(nf, nb), &
                                                            v(nf, nb), vs(nf, nb, 2)
      INTEGER                                            :: axis, ico, image, images(3, 3), point, &
                                                            row, side, translated_images(3, 3)
      REAL(dp) :: bloch(3), cell(3, 3), center(3), cg(3, 8), cm(8), cp(8), cv(8), cx(3), dc(3, 3), &
         dcenter(3), dk(3), dpnt(3, nf), gradient(nf, nb, 3), gs(nf, nb, 3, 2), hh(nf, nb, 3, 3), &
         radius, samples(nf, nb, 2), shift, value(nf, nb), x(3, nf), xp(3, nf)
      TYPE(gaussian_coulomb_fields_type)                 :: coulomb
      TYPE(gto_basis_set_type), POINTER                  :: basis

      NULLIFY (basis)
      CALL init_orbital_pointers(3)
      CALL allocate_gto_basis_set(basis)
      basis%nset = 2
      basis%nsgf = nb
      ALLOCATE (basis%lmax(2), basis%npgf(2), basis%nsgf_set(2), basis%nshell(2), basis%l(2, 2), &
                basis%first_sgf(2, 2), basis%last_sgf(2, 2), basis%zet(2, 2), basis%sphi(40, nb))
      basis%lmax = [3, 1]
      basis%npgf = [2, 1]
      basis%nsgf_set = [8, 3]
      basis%nshell = [2, 1]
      basis%l = 0
      basis%l(:, 1) = [0, 3]
      basis%l(1, 2) = 1
      basis%first_sgf = 0
      basis%last_sgf = 0
      basis%first_sgf(:, 1) = [1, 2]
      basis%last_sgf(:, 1) = [1, 8]
      basis%first_sgf(1, 2) = 9
      basis%last_sgf(1, 2) = 11
      basis%zet = 0.0_dp
      basis%zet(:, 1) = [.4_dp, .9_dp]
      basis%zet(1, 2) = .7_dp
      basis%sphi = 0.0_dp
      basis%sphi(1, 1) = .6_dp
      basis%sphi(21, 1) = -.2_dp
      DO row = 2, 8
         DO ico = ncoset(2) + 1, ncoset(3)
            basis%sphi(ico, row) = .03_dp*SIN(REAL(ico*row, dp))
            basis%sphi(ico + 20, row) = -.02_dp*COS(REAL(ico + row, dp))
         END DO
      END DO
      DO row = 9, nb
         DO ico = 2, 4
            basis%sphi(ico, row) = .1_dp*SIN(REAL(ico*row, dp))
         END DO
      END DO
      x(:, 1) = 0.0_dp
      x(:, 2) = [.3_dp, 0.0_dp, 0.0_dp]
      x(:, 3) = [0.0_dp, -.7_dp, 0.0_dp]
      x(:, 4) = [.23_dp, .51_dp, -.4_dp]
      x(:, 5) = [-.6_dp, .2_dp, .8_dp]
      CALL gaussian_basis_fields_evaluate(basis, x, value, gradient, hh)
      CALL check(MAXVAL(ABS(value(:, 1) - (.6_dp*EXP(-.4_dp*SUM(x*x, DIM=1)) &
                                           - .2_dp*EXP(-.9_dp*SUM(x*x, DIM=1))))), 'native signed s contraction')
      DO axis = 1, 3
         DO side = 1, 2
            xp = x
            xp(axis, :) = xp(axis, :) + (2*side - 3)*step
            CALL gaussian_basis_fields_evaluate(basis, xp, samples(:, :, side), gs(:, :, :, side))
         END DO
         CALL check(MAXVAL(ABS(gradient(:, :, axis) - (samples(:, :, 2) - samples(:, :, 1))/(2*step))), &
                    'contracted Gaussian point gradient', 2.e-8_dp)
         CALL check(MAXVAL(ABS(hh(:, :, :, axis) - (gs(:, :, :, 2) - gs(:, :, :, 1))/(2*step))), &
                    'contracted Gaussian origin/axis Hessian', 2.e-7_dp)
      END DO
      cell = RESHAPE([2.3_dp, .1_dp, .2_dp, .2_dp, 2.7_dp, -.1_dp, .1_dp, .3_dp, 2.1_dp], [3, 3])
      center = [.1_dp, -.2_dp, .3_dp]
      images = RESHAPE([0, 0, 0, 1, -1, 0, -1, 0, 1], [3, 3])
      bloch = [.23_dp, -.17_dp, .11_dp]
      dc = .03_dp*cell
      dcenter = [.07_dp, -.02_dp, .01_dp]
      dpnt = .02_dp*x
      dk = [.03_dp, .02_dp, -.01_dp]
      CALL gaussian_basis_fields_periodic(basis, x, center, cell, images, bloch, v, g, dv, dg, dpnt, dcenter, dc, dk)
      DO side = 1, 2
         shift = (2*side - 3)*step
         CALL gaussian_basis_fields_periodic(basis, x + shift*dpnt, center + shift*dcenter, cell + shift*dc, &
                                             images, bloch + shift*dk, vs(:, :, side), gss(:, :, :, side))
      END DO
      CALL check(MAXVAL(ABS(dv - (vs(:, :, 2) - vs(:, :, 1))/(2*step))), 'native Gaussian Bloch joint response', 2.e-8_dp)
     CALL check(MAXVAL(ABS(dg - (gss(:, :, :, 2) - gss(:, :, :, 1))/(2*step))), 'native Gaussian Bloch gradient response', 2.e-7_dp)
      CALL gaussian_basis_fields_periodic(basis, x, center, cell, images, -bloch, other)
      CALL check(MAXVAL(ABS(other - CONJG(v))), 'native Gaussian Bloch time reversal')
      DO image = 1, 3
         translated_images(:, image) = images(:, image) - [1, 0, 0]
      END DO
      CALL gaussian_basis_fields_periodic(basis, x, center + cell(:, 1), cell, translated_images, bloch, other)
      CALL check(MAXVAL(ABS(other - EXP(CMPLX(0.0_dp, -2*pi*bloch(1), dp))*v)), 'native Gaussian lattice representative')
      CALL gaussian_coulomb_fields_init(3)
      CALL gaussian_coulomb_fields_prepare(basis, 1, coulomb)
      CALL deallocate_gto_basis_set(basis)
      DO point = 1, 3
         cx = [SQRT(12.0_dp/.4_dp)*REAL(point - 1, dp), 0.0_dp, 0.0_dp]
         CALL gaussian_coulomb_fields_evaluate(coulomb, cx, cv, cg)
         radius = SQRT(SUM(cx**2))
         IF (radius == 0.0_dp) THEN
            CALL check(ABS(cv(1) - 2*pi*(.6_dp/.4_dp - .2_dp/.9_dp)), 'Gaussian Coulomb origin')
         ELSE
            CALL check(ABS(cv(1) - pi**1.5_dp*(.6_dp*ERF(SQRT(.4_dp)*radius)/.4_dp**1.5_dp &
                                               - .2_dp*ERF(SQRT(.9_dp)*radius)/.9_dp**1.5_dp)/radius), &
                       'Signed Gaussian Coulomb tail')
         END IF
         DO axis = 1, 3
            cx(axis) = cx(axis) + step
            CALL gaussian_coulomb_fields_evaluate(coulomb, cx, cp)
            cx(axis) = cx(axis) - 2*step
            CALL gaussian_coulomb_fields_evaluate(coulomb, cx, cm)
            cx(axis) = cx(axis) + step
            CALL check(MAXVAL(ABS((cp - cm)/(2*step) - cg(axis, :))), 'Gaussian Coulomb gradient', 2.e-8_dp)
         END DO
      END DO
   END SUBROUTINE check_gaussian_fields

! **************************************************************************************************
!> \brief ...
!> \param value ...
!> \param label ...
!> \param tolerance ...
! **************************************************************************************************
   SUBROUTINE check(value, label, tolerance)
      REAL(KIND=dp), INTENT(IN)                          :: value
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: tolerance

      REAL(KIND=dp)                                      :: tol

      tol = 2.e-12_dp
      IF (PRESENT(tolerance)) tol = tolerance
      max_error = MAX(max_error, value)
      IF (value > tol) THEN
         WRITE (*, *) TRIM(label), value, tol
         ERROR STOP 'Inverse primitive mismatch'
      END IF
   END SUBROUTINE check
END PROGRAM upf_inverse_unittest
