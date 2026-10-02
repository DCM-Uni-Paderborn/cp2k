!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM uspp_hfx_integrals_unittest
   USE ai_coulomb,                      ONLY: coulomb2
   USE atom_upf,                        ONLY: atom_release_upf,&
                                              atom_upfpot_type
   USE basis_set_types,                 ONLY: allocate_gto_basis_set,&
                                              deallocate_gto_basis_set,&
                                              gto_basis_set_type
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: dp
   USE lebedev,                         ONLY: deallocate_lebedev_grids,&
                                              get_number_of_lebedev_grid,&
                                              init_lebedev_grids,&
                                              lebedev_grid
   USE libcp2k,                         ONLY: cp2k_finalize,&
                                              cp2k_init
   USE mathconstants,                   ONLY: fourpi,&
                                              pi
   USE orbital_pointers,                ONLY: coset,&
                                              indco,&
                                              ncoset
   USE orbital_transformation_matrices, ONLY: orbtramat
   USE upf_projector_integrals,         ONLY: upf_projector_init,&
                                              upf_projector_type
   USE uspp_augmentation,               ONLY: uspp_augmentation_init,&
                                              uspp_augmentation_type
   USE uspp_hfx_integrals,              ONLY: uspp_hfx_auxiliary_integrals,&
                                              uspp_hfx_integrals_init,&
                                              uspp_hfx_primitive_integrals

   IMPLICIT NONE
   REAL(dp), PARAMETER                                :: alpha = 0.75_dp, beta = 0.9_dp, step = 2.e-5_dp
   INTEGER, PARAMETER                                 :: np = 4, projector_axis(np) = [0, 2, 3, 1]
   INTEGER                                            :: axis, grid, l, rule, sample
   REAL(dp)                                           :: dg(np, np, 5, 3), dq(np, np, 5, 3), &
                                                         reference(np, np, 5), dref(np, np, 5, 3), &
                                                         integral(np, np, 5), minus(np, np, 5), &
                                                         plus(np, np, 5), displacement(3), shifted(3), &
                                                         max_error, closed_ss, kappa, t, f0, prefactor
   TYPE(uspp_augmentation_type)                       :: augmentation

   CALL cp2k_init()
   CALL uspp_hfx_integrals_init(6)
   CALL init_lebedev_grids()
   rule = get_number_of_lebedev_grid(n=110)
   max_error = 0.0_dp
   DO grid = 1, 2
      CALL make_augmentation(grid, augmentation)
      DO sample = 1, 2
         displacement = REAL(sample - 1, dp)*[0.31_dp, -0.52_dp, 0.73_dp]
         DO l = 0, 2
            CALL uspp_hfx_primitive_integrals(augmentation, alpha, 0, l, displacement, lebedev_grid(rule)%r, &
                                              fourpi*lebedev_grid(rule)%w, integral(:, :, :2*l + 1), &
                                              dg(:, :, :2*l + 1, :), dq(:, :, :2*l + 1, :))
            CALL primitive_reference(alpha, l, displacement, reference(:, :, :2*l + 1), dref(:, :, :2*l + 1, :))
            IF (.NOT. ALL(ieee_is_finite(integral(:, :, :2*l + 1))) .OR. &
                .NOT. ALL(ieee_is_finite(dg(:, :, :2*l + 1, :)))) ERROR STOP 'Nonfinite native RI-Q primitive'
            CALL check_error(MAXVAL(ABS(integral(:, :, :2*l + 1) - reference(:, :, :2*l + 1))), 3.e-8_dp)
            CALL check_error(MAXVAL(ABS(dg(:, :, :2*l + 1, :) - dref(:, :, :2*l + 1, :))), 3.e-8_dp)
            IF (ANY(dg(:, :, :2*l + 1, :) + dq(:, :, :2*l + 1, :) /= 0.0_dp)) &
               ERROR STOP 'Native RI-Q derivatives violate translational invariance'
            DO axis = 1, 2*l + 1
               IF (ANY(integral(:, :, axis) /= TRANSPOSE(integral(:, :, axis)))) &
                  ERROR STOP 'Native RI-Q integral lost projector-pair symmetry'
            END DO
            IF (l /= 0) CYCLE
            ! Closed Gaussian-Gaussian formula independent of both production paths.
            kappa = alpha*beta/(alpha + beta)
            t = kappa*SUM(displacement**2)
            f0 = 1.0_dp
            IF (t > 0.0_dp) f0 = SQRT(pi)*ERF(SQRT(t))/(2.0_dp*SQRT(t))
            prefactor = 2.0_dp*pi**2.5_dp/(alpha*beta*SQRT(alpha + beta))
            closed_ss = prefactor*f0/(fourpi*SQRT(fourpi))
            CALL check_error(ABS(integral(1, 1, 1) - closed_ss), 3.e-8_dp)
         END DO
         IF (grid /= 1 .OR. sample /= 2) CYCLE
         CALL uspp_hfx_primitive_integrals(augmentation, alpha, 0, 1, displacement, lebedev_grid(rule)%r, &
                                           fourpi*lebedev_grid(rule)%w, integral(:, :, :3), &
                                           augmentation_derivative=dq(:, :, :3, :))
         DO axis = 1, 3
            shifted = displacement
            shifted(axis) = shifted(axis) + step
            CALL uspp_hfx_primitive_integrals(augmentation, alpha, 0, 1, shifted, lebedev_grid(rule)%r, &
                                              fourpi*lebedev_grid(rule)%w, plus(:, :, :3))
            shifted(axis) = shifted(axis) - 2.0_dp*step
            CALL uspp_hfx_primitive_integrals(augmentation, alpha, 0, 1, shifted, lebedev_grid(rule)%r, &
                                              fourpi*lebedev_grid(rule)%w, minus(:, :, :3))
            CALL check_error(MAXVAL(ABS((plus(:, :, :3) - minus(:, :, :3))/(2.0_dp*step) - dq(:, :, :3, axis))), &
                             3.e-8_dp)
         END DO
      END DO
      IF (grid == 2) CALL check_auxiliary(augmentation, rule)
   END DO
   PRINT '(A,ES16.7)', 'Native USPP RI-Q Coulomb maximum absolute error: ', max_error
   CALL deallocate_lebedev_grids()
   CALL cp2k_finalize()

CONTAINS

! **************************************************************************************************
!> \brief Native s/p Q pairs equal products of Cartesian Gaussian polynomials.
!>        Linear and exponential meshes retain the original UPF quadrature and Jacobian.
! **************************************************************************************************
   SUBROUTINE make_augmentation(grid, augmentation)
      INTEGER, INTENT(IN)                                :: grid
      TYPE(uspp_augmentation_type), INTENT(OUT)          :: augmentation

      INTEGER                                            :: i, j, l, nr, pair
      REAL(dp)                                           :: dx
      TYPE(atom_upfpot_type)                             :: pot
      TYPE(upf_projector_type)                           :: projectors

      nr = 601
      IF (grid == 2) nr = 801
      pot%is_ultrasoft = .TRUE.
      pot%number_of_proj = 2
      ALLOCATE (pot%lbeta(2), pot%dion(2, 2), pot%qij(2, 2), pot%r(nr), pot%rab(nr), pot%beta(nr, 2), &
                pot%qijl(nr, 3, 0:2))
      pot%lbeta(:) = [0, 1]
      pot%dion = 0.0_dp
      pot%qij = 0.0_dp
      DO i = 1, nr
         IF (grid == 1) THEN
            pot%r(i) = 8.0_dp*REAL(i - 1, dp)/REAL(nr - 1, dp)
            pot%rab(i) = 8.0_dp/REAL(nr - 1, dp)
         ELSE
            dx = LOG(1.0_dp + 8.0_dp/0.001_dp)/REAL(nr - 1, dp)
            pot%r(i) = 0.001_dp*(EXP(REAL(i - 1, dp)*dx) - 1.0_dp)
            pot%rab(i) = dx*(pot%r(i) + 0.001_dp)
         END IF
      END DO
      pot%qijl = 0.0_dp
      DO j = 1, 2
         pot%beta(:, j) = pot%r**pot%lbeta(j)*EXP(-pot%r**2)
         DO i = 1, j
            pair = j*(j - 1)/2 + i
            DO l = ABS(pot%lbeta(i) - pot%lbeta(j)), pot%lbeta(i) + pot%lbeta(j), 2
               pot%qijl(:, pair, l) = pot%r**(pot%lbeta(i) + pot%lbeta(j) + 2)*EXP(-beta*pot%r**2)
            END DO
         END DO
      END DO
      CALL upf_projector_init(pot, 2, projectors)
      CALL uspp_augmentation_init(pot, projectors, augmentation)
      IF (augmentation%nproj /= np) ERROR STOP 'Wrong native RI-Q magnetic-projector count'
      CALL atom_release_upf(pot)
   END SUBROUTINE make_augmentation

! **************************************************************************************************
!> \brief Independent analytic Cartesian two-center Coulomb recurrence, with explicit Gaussian derivative.
!>        Q_s,s=g/(4*pi), Q_s,p_i=sqrt(3)*x_i*g/(4*pi), Q_p_i,p_j=3*x_i*x_j*g/(4*pi).
! **************************************************************************************************
   SUBROUTINE cartesian_reference(alpha, degree, displacement, values, derivatives)
      REAL(dp), INTENT(IN)                               :: alpha
      INTEGER, INTENT(IN)                                :: degree
      REAL(dp), INTENT(IN)                               :: displacement(3)
      REAL(dp), INTENT(OUT)                              :: values(:, :, :), derivatives(:, :, :, :)

      INTEGER                                            :: axis, c, i, p, powers(3), q, raised(3), t
      REAL(dp) :: f(0:degree + 3), factor, raw(ncoset(degree + 1), ncoset(2)), &
         scratch(ncoset(degree + 1), ncoset(2), degree + 4)

      raw = 0.0_dp
      CALL coulomb2(degree + 1, 1, [alpha], [0.0_dp], 0, 2, 1, [beta], [0.0_dp], 0, &
                    displacement, SUM(displacement**2), raw, scratch, f, screening=.FALSE.)
      DO q = 1, np
         DO p = 1, np
            powers = 0
            factor = 1.0_dp/fourpi
            IF (p > 1) THEN
               powers(projector_axis(p)) = powers(projector_axis(p)) + 1
               factor = factor*SQRT(3.0_dp)
            END IF
            IF (q > 1) THEN
               powers(projector_axis(q)) = powers(projector_axis(q)) + 1
               factor = factor*SQRT(3.0_dp)
            END IF
            c = coset(powers(1), powers(2), powers(3))
            DO t = 1, ncoset(degree)
               values(t, p, q) = factor*raw(t, c)
               DO axis = 1, 3
                  raised = indco(:, t)
                  raised(axis) = raised(axis) + 1
                  i = coset(raised(1), raised(2), raised(3))
                  derivatives(axis, t, p, q) = 2.0_dp*alpha*factor*raw(i, c)
                  raised = indco(:, t)
                  IF (raised(axis) == 0) CYCLE
                  raised(axis) = raised(axis) - 1
                  i = coset(raised(1), raised(2), raised(3))
                  derivatives(axis, t, p, q) = derivatives(axis, t, p, q) - &
                                               REAL(indco(axis, t), dp)*factor*raw(i, c)
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE cartesian_reference

! **************************************************************************************************
!> \brief Contract the independent Cartesian oracle to one normalized angular multipole.
! **************************************************************************************************
   SUBROUTINE primitive_reference(alpha, l, displacement, reference, derivative)
      REAL(dp), INTENT(IN)                               :: alpha
      INTEGER, INTENT(IN)                                :: l
      REAL(dp), INTENT(IN)                               :: displacement(3)
      REAL(dp), INTENT(OUT)                              :: reference(:, :, :), &
                                                            derivative(:, :, :, :)

      INTEGER                                            :: axis, m, t
      REAL(dp)                                           :: draw(3, ncoset(l), np, np), &
                                                            raw(ncoset(l), np, np)

      CALL cartesian_reference(alpha, l, displacement, raw, draw)
      reference = 0.0_dp
      derivative = 0.0_dp
      DO m = 1, 2*l + 1
         DO t = ncoset(l - 1) + 1, ncoset(l)
            reference(:, :, m) = reference(:, :, m) + orbtramat(l)%slm(m, t - ncoset(l - 1))*raw(t, :, :)
            DO axis = 1, 3
               derivative(:, :, m, axis) = derivative(:, :, m, axis) + &
                                           orbtramat(l)%slm(m, t - ncoset(l - 1))*draw(axis, t, :, :)
            END DO
         END DO
      END DO
   END SUBROUTINE primitive_reference

! **************************************************************************************************
!> \brief Full sphi contraction, mixed sets, Cartesian trace channel and both-center energy adjoints.
! **************************************************************************************************
   SUBROUTINE check_auxiliary(augmentation, rule)
      TYPE(uspp_augmentation_type), INTENT(IN)           :: augmentation
      INTEGER, INTENT(IN)                                :: rule

      INTEGER                                            :: axis, first, ip, iset, l, m, naux, p, q, &
                                                            row, t
      REAL(dp)                                           :: center_g(3), center_q(3), &
                                                            direction_g(3), direction_q(3), &
                                                            energy_derivative, energy_fd, factor
      REAL(dp), ALLOCATABLE :: adjoint(:, :, :), dg(:, :, :, :), dq(:, :, :, :), draw(:, :, :, :), &
         dref(:, :, :, :), integral(:, :, :), minus(:, :, :), plus(:, :, :), raw(:, :, :), &
         reference(:, :, :)
      TYPE(gto_basis_set_type), POINTER                  :: basis

      NULLIFY (basis)
      CALL allocate_gto_basis_set(basis)
      basis%nset = 2
      basis%nsgf = 11
      ALLOCATE (basis%lmax(2), basis%npgf(2), basis%nsgf_set(2), basis%zet(2, 2), basis%sphi(20, 11))
      basis%lmax(:) = [2, 0]
      basis%npgf(:) = [2, 1]
      basis%nsgf_set(:) = [10, 1]
      basis%zet(:, 1) = [0.75_dp, 1.6_dp]
      basis%zet(:, 2) = [2.3_dp, 0.0_dp]
      basis%sphi = 0.0_dp
      DO ip = 1, 2
         factor = 0.7_dp
         IF (ip == 2) factor = -0.23_dp
         DO l = 0, 2
            DO m = 1, 2*l + 1
               DO t = ncoset(l - 1) + 1, ncoset(l)
                  row = (ip - 1)*ncoset(2) + t
                  basis%sphi(row, l*l + m) = factor*REAL(l + 1, dp)*orbtramat(l)%slm(m, t - ncoset(l - 1))
               END DO
            END DO
         END DO
         ! A Cartesian x**2 function exercises the l=0 radial trace term as well as l=2.
         row = (ip - 1)*ncoset(2) + coset(2, 0, 0)
         basis%sphi(row, 10) = factor
      END DO
      basis%sphi(1, 11) = 0.83_dp
      center_g = [-0.21_dp, 0.12_dp, 0.35_dp]
      center_q = [0.10_dp, -0.40_dp, 1.08_dp]
      direction_g = [0.2_dp, -0.3_dp, 0.1_dp]
      direction_q = [-0.1_dp, 0.2_dp, 0.3_dp]
      DO iset = 1, 2
         naux = basis%nsgf_set(iset)
         ALLOCATE (integral(naux, np, np), plus(naux, np, np), minus(naux, np, np), reference(naux, np, np), &
                   adjoint(naux, np, np), dg(3, naux, np, np), dq(3, naux, np, np), dref(3, naux, np, np), &
                   raw(ncoset(basis%lmax(iset)), np, np), draw(3, ncoset(basis%lmax(iset)), np, np))
         CALL uspp_hfx_auxiliary_integrals(augmentation, basis, iset, center_g, center_q, lebedev_grid(rule)%r, &
                                           fourpi*lebedev_grid(rule)%w, integral, dg, dq)
         IF (.NOT. ALL(ieee_is_finite(integral)) .OR. .NOT. ALL(ieee_is_finite(dg))) &
            ERROR STOP 'Nonfinite contracted RI-Q integral'
         reference = 0.0_dp
         dref = 0.0_dp
         first = SUM(basis%nsgf_set(:iset - 1)) + 1
         DO ip = 1, basis%npgf(iset)
            CALL cartesian_reference(basis%zet(ip, iset), basis%lmax(iset), center_q - center_g, raw, draw)
            DO t = 1, ncoset(basis%lmax(iset))
               row = (ip - 1)*ncoset(basis%lmax(iset)) + t
               DO m = 1, naux
                  reference(m, :, :) = reference(m, :, :) + basis%sphi(row, first + m - 1)*raw(t, :, :)
                  DO axis = 1, 3
                     dref(axis, m, :, :) = dref(axis, m, :, :) + &
                                           basis%sphi(row, first + m - 1)*draw(axis, t, :, :)
                  END DO
               END DO
            END DO
         END DO
         CALL check_error(MAXVAL(ABS(integral - reference)), 5.e-8_dp)
         CALL check_error(MAXVAL(ABS(dg - dref)), 5.e-8_dp)
         IF (ANY(dg + dq /= 0.0_dp)) ERROR STOP 'Contracted RI-Q center gradients do not cancel'
         CALL uspp_hfx_auxiliary_integrals(augmentation, basis, iset, center_g + direction_g, center_q + direction_g, &
                                           lebedev_grid(rule)%r, fourpi*lebedev_grid(rule)%w, plus)
         CALL check_error(MAXVAL(ABS(plus - integral)), 2.e-13_dp)
         DO q = 1, np
            DO p = 1, np
               DO m = 1, naux
                  adjoint(m, p, q) = SIN(REAL(m + p + q, dp))
               END DO
            END DO
         END DO
         energy_derivative = 0.0_dp
         DO axis = 1, 3
            energy_derivative = energy_derivative + SUM(adjoint*(direction_g(axis)*dg(axis, :, :, :) + &
                                                                 direction_q(axis)*dq(axis, :, :, :)))
         END DO
         CALL uspp_hfx_auxiliary_integrals(augmentation, basis, iset, center_g + step*direction_g, &
                                           center_q + step*direction_q, lebedev_grid(rule)%r, &
                                           fourpi*lebedev_grid(rule)%w, plus)
         CALL uspp_hfx_auxiliary_integrals(augmentation, basis, iset, center_g - step*direction_g, &
                                           center_q - step*direction_q, lebedev_grid(rule)%r, &
                                           fourpi*lebedev_grid(rule)%w, minus)
         energy_fd = SUM(adjoint*(plus - minus))/(2.0_dp*step)
         CALL check_error(ABS(energy_fd - energy_derivative), 8.e-8_dp)
         DEALLOCATE (integral, plus, minus, reference, adjoint, dg, dq, dref, raw, draw)
      END DO
      CALL deallocate_gto_basis_set(basis)
   END SUBROUTINE check_auxiliary

! **************************************************************************************************
!> \brief Check an absolute error in the finite synthetic integral fixture.
! **************************************************************************************************
   SUBROUTINE check_error(error, tolerance)
      REAL(dp), INTENT(IN)                               :: error, tolerance

      max_error = MAX(max_error, error)
      IF (.NOT. ieee_is_finite(error) .OR. error > tolerance) THEN
         PRINT *, 'Native RI-Q error/tolerance:', error, tolerance
         ERROR STOP 'Native USPP Coulomb integral or derivative check failed'
      END IF
   END SUBROUTINE check_error

END PROGRAM uspp_hfx_integrals_unittest
