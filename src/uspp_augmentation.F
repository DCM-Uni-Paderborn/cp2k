!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Direct reciprocal-space reference for scalar UPF augmentation.
!>        Uses the complete tabulated radial mesh and CP2K real spherical harmonics.
!>        The overlap correction is the zero-frequency limit of this same transform.
! **************************************************************************************************
MODULE uspp_augmentation
   USE atom_upf,                        ONLY: atom_upfpot_type,&
                                              upf_radial_weights
   USE bessel_lib,                      ONLY: bessel0
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: fourpi
   USE orbital_pointers,                ONLY: indco,&
                                              init_orbital_pointers,&
                                              nco,&
                                              ncoset
   USE orbital_transformation_matrices, ONLY: init_spherical_harmonics,&
                                              orbtramat
   USE spherical_harmonics,             ONLY: y_lm
   USE splines,                         ONLY: spline3ders
   USE upf_projector_integrals,         ONLY: upf_projector_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'uspp_augmentation'

   TYPE uspp_augmentation_type
      INTEGER :: nproj = 0, lmax = -1
      INTEGER, ALLOCATABLE :: radial_index(:), angular_index(:), lbeta(:)
      REAL(KIND=dp), ALLOCATABLE :: r(:), weighted_q(:, :, :), gaunt(:, :, :), q(:, :)
      REAL(KIND=dp), ALLOCATABLE :: radial_x(:), radial_value(:, :, :), radial_slope(:, :, :)
   END TYPE uspp_augmentation_type

   PUBLIC :: uspp_augmentation_type, uspp_augmentation_init, uspp_augmentation_fourier, &
             uspp_augmentation_collocate, uspp_augmentation_integrate, &
             uspp_augmentation_prepare_realspace, uspp_augmentation_realspace, &
             uspp_augmentation_collocate_atom, uspp_augmentation_integrate_atom

CONTAINS

! **************************************************************************************************
!> \brief Prepare radial quadrature and Gaunt coefficients in the native beta convention.
!> \param pot Original scalar USPP with angular-resolved augmentation.
!> \param projectors Prepared native projectors, with Cartesian degree at least max(lbeta).
!> \param augmentation Immutable augmentation data, in the same magnetic-projector ordering.
! **************************************************************************************************
   SUBROUTINE uspp_augmentation_init(pot, projectors, augmentation)
      TYPE(atom_upfpot_type), INTENT(IN)                 :: pot
      TYPE(upf_projector_type), INTENT(IN)               :: projectors
      TYPE(uspp_augmentation_type), INTENT(OUT)          :: augmentation

      INTEGER                                            :: i, j, l, m, n, nc, p, pair, q, s
      REAL(KIND=dp), ALLOCATABLE                         :: weights(:)

      CPASSERT(.NOT. pot%has_so .AND. .NOT. pot%is_paw)
      CPASSERT(ALLOCATED(pot%qijl))
      CPASSERT(projectors%degree >= MAXVAL(pot%lbeta))
      CPASSERT(ALL(projectors%lbeta == pot%lbeta))
      n = SIZE(pot%r)
      CPASSERT(n >= 3 .AND. SIZE(pot%rab) == n)
      CPASSERT(ALL(pot%r >= 0.0_dp) .AND. ALL(pot%rab > 0.0_dp))
      CPASSERT(ALL(pot%r(2:) > pot%r(:n - 1)))
      nc = pot%number_of_proj*(pot%number_of_proj + 1)/2
      augmentation%nproj = projectors%nproj
      augmentation%lmax = 2*MAXVAL(pot%lbeta)
      CALL init_orbital_pointers(augmentation%lmax)
      CALL init_spherical_harmonics(augmentation%lmax, -1)
      CPASSERT(ALL(SHAPE(pot%qijl) == [n, nc, augmentation%lmax + 1]))
      CPASSERT(LBOUND(pot%qijl, 3) == 0)
      ALLOCATE (augmentation%r(n), weights(n))
      augmentation%r(:) = pot%r
      weights(:) = upf_radial_weights(pot)
      ALLOCATE (augmentation%weighted_q(n, nc, 0:augmentation%lmax))
      DO l = 0, augmentation%lmax
         DO pair = 1, nc
            augmentation%weighted_q(:, pair, l) = weights*pot%qijl(:, pair, l)
         END DO
      END DO
      n = augmentation%nproj
      ALLOCATE (augmentation%radial_index(n), augmentation%angular_index(n), &
                augmentation%lbeta(pot%number_of_proj))
      augmentation%radial_index(:) = projectors%radial_index
      augmentation%angular_index(:) = projectors%angular_index
      augmentation%lbeta(:) = pot%lbeta
      ALLOCATE (augmentation%gaunt(n, n, (augmentation%lmax + 1)**2), augmentation%q(n, n))
      augmentation%gaunt = 0.0_dp
      augmentation%q = 0.0_dp
      DO q = 1, n
         j = projectors%radial_index(q)
         l = pot%lbeta(j)
         m = projectors%angular_index(q) - l*l
         DO p = 1, n
            i = projectors%radial_index(p)
            ! Contract the existing analytic monomial/Y/Y integrals with Y_j's polynomial.
            DO s = 1, nco(l)
               augmentation%gaunt(p, q, :) = augmentation%gaunt(p, q, :) + &
                                             orbtramat(l)%slm(m, s)*projectors%gaunt(ncoset(l - 1) + s, &
                                                                           projectors%angular_index(p), :(augmentation%lmax + 1)**2)
            END DO
            pair = MAX(i, j)*(MAX(i, j) - 1)/2 + MIN(i, j)
            IF (projectors%angular_index(p) == projectors%angular_index(q)) THEN
               augmentation%q(p, q) = SUM(augmentation%weighted_q(:, pair, 0))
            END IF
         END DO
      END DO
   END SUBROUTINE uspp_augmentation_init

! **************************************************************************************************
!> \brief Prepare regular radial splines for atom-grid augmentation quadrature.
!>        UPF stores r**2 Q_ij^L(r). Interpolate Q_ij^L(r)/r**L in r**2 and multiply
!>        by solid harmonics, so values and Cartesian derivatives remain regular at the origin.
!>        The existing reciprocal-space quadrature and its overlap correction are unchanged.
!> \param pot Original scalar USPP used to initialize augmentation.
!> \param augmentation Prepared angular data; owns the optional real-space cache on return.
! **************************************************************************************************
   SUBROUTINE uspp_augmentation_prepare_realspace(pot, augmentation)
      TYPE(atom_upfpot_type), INTENT(IN)                 :: pot
      TYPE(uspp_augmentation_type), INTENT(INOUT)        :: augmentation

      INTEGER                                            :: first, l, n, nc, pair
      REAL(KIND=dp), ALLOCATABLE                         :: values(:), x(:)

      CPASSERT(ALLOCATED(pot%qijl) .AND. ALLOCATED(augmentation%r))
      CPASSERT(SIZE(pot%r) == SIZE(augmentation%r))
      CPASSERT(ALL(pot%r == augmentation%r))
      CPASSERT(ALL(pot%lbeta == augmentation%lbeta))
      IF (ALLOCATED(augmentation%radial_x)) RETURN
      first = 1
      IF (pot%r(1) == 0.0_dp) first = 2
      n = SIZE(pot%r) - first + 1
      CPASSERT(n >= 4)
      nc = SIZE(pot%qijl, 2)
      ALLOCATE (x(n), values(n), augmentation%radial_x(n + 1), &
                augmentation%radial_value(n + 1, nc, 0:augmentation%lmax), &
                augmentation%radial_slope(n + 1, nc, 0:augmentation%lmax))
      x(:) = pot%r(first:)**2
      augmentation%radial_x(1) = 0.0_dp
      augmentation%radial_x(2:) = x
      augmentation%radial_value = 0.0_dp
      augmentation%radial_slope = 0.0_dp
      DO l = 0, augmentation%lmax
         DO pair = 1, nc
            IF (ALL(pot%qijl(:, pair, l) == 0.0_dp)) CYCLE
            values(:) = pot%qijl(first:, pair, l)/pot%r(first:)**(l + 2)
            ! The first interval continues the endpoint cubic to x=0. A tabulated
            ! r**2 Q value at r=0 alone does not determine the regular radial limit.
            CALL spline3ders(x, values, augmentation%radial_x, &
                             ynew=augmentation%radial_value(:, pair, l), &
                             dynew=augmentation%radial_slope(:, pair, l))
         END DO
      END DO
   END SUBROUTINE uspp_augmentation_prepare_realspace

! **************************************************************************************************
!> \brief Evaluate the augmentation matrix and derivatives at a Cartesian displacement.
!>        Uses the same cubic for values, gradients and Hessians. No angular division
!>        by radius is needed, including at the center and on coordinate axes.
!> \param augmentation Data with an initialized real-space cache.
!> \param displacement Evaluation point minus augmentation center, in bohr.
!> \param qr Symmetric Q_ij(r), without radial Jacobian or quadrature weights.
!> \param gradient Optional Cartesian gradient of qr with respect to the evaluation point.
!> \param hessian Optional Cartesian Hessian of qr, also needed for GGA geometry derivatives.
! **************************************************************************************************
   SUBROUTINE uspp_augmentation_realspace(augmentation, displacement, qr, gradient, hessian)
      TYPE(uspp_augmentation_type), INTENT(IN)           :: augmentation
      REAL(KIND=dp), INTENT(IN)                          :: displacement(3)
      REAL(KIND=dp), INTENT(OUT)                         :: qr(:, :)
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: gradient(:, :, :), hessian(:, :, :, :)

      INTEGER                                            :: a, axis, axis2, hi, i, j, l, lm, lo, m, &
                                                            mid, p, pair, powers(3), q, reduced(3)
      REAL(KIND=dp) :: angular, angular_gradient(3), angular_hessian(3, 3), c0, c1, c2, c3, &
         factor, h, radial, radial_d, radial_dd, s, t, term, value, value_gradient(3), &
         value_hessian(3, 3)
      REAL(KIND=dp), ALLOCATABLE                         :: ddharmonics(:, :, :), dharmonics(:, :), &
                                                            harmonics(:)

      CPASSERT(ALLOCATED(augmentation%radial_x))
      CPASSERT(ALL(SHAPE(qr) == augmentation%nproj))
      qr = 0.0_dp
      IF (PRESENT(gradient)) THEN
         CPASSERT(ALL(SHAPE(gradient) == [augmentation%nproj, augmentation%nproj, 3]))
         gradient = 0.0_dp
      END IF
      IF (PRESENT(hessian)) THEN
         CPASSERT(ALL(SHAPE(hessian) == [augmentation%nproj, augmentation%nproj, 3, 3]))
         hessian = 0.0_dp
      END IF
      s = SUM(displacement**2)
      hi = SIZE(augmentation%radial_x)
      IF (s > augmentation%radial_x(hi)) RETURN
      lo = 1
      DO WHILE (hi - lo > 1)
         mid = (lo + hi)/2
         IF (s < augmentation%radial_x(mid)) THEN
            hi = mid
         ELSE
            lo = mid
         END IF
      END DO
      h = augmentation%radial_x(hi) - augmentation%radial_x(lo)
      t = (s - augmentation%radial_x(lo))/h
      ALLOCATE (harmonics((augmentation%lmax + 1)**2), &
                dharmonics(3, (augmentation%lmax + 1)**2), &
                ddharmonics(3, 3, (augmentation%lmax + 1)**2))
      harmonics = 0.0_dp
      dharmonics = 0.0_dp
      ddharmonics = 0.0_dp
      DO l = 0, augmentation%lmax
         DO m = -l, l
            lm = l*l + l + m + 1
            DO a = 1, nco(l)
               powers = indco(:, ncoset(l - 1) + a)
               factor = orbtramat(l)%slm(l + m + 1, a)
               harmonics(lm) = harmonics(lm) + factor*PRODUCT(displacement**powers)
               DO axis = 1, 3
                  IF (powers(axis) == 0) CYCLE
                  reduced = powers
                  reduced(axis) = reduced(axis) - 1
                  term = factor*REAL(powers(axis), dp)
                  dharmonics(axis, lm) = dharmonics(axis, lm) + term*PRODUCT(displacement**reduced)
                  DO axis2 = 1, 3
                     IF (reduced(axis2) == 0) CYCLE
                     reduced(axis2) = reduced(axis2) - 1
                     ddharmonics(axis, axis2, lm) = ddharmonics(axis, axis2, lm) + &
                                                    term*REAL(reduced(axis2) + 1, dp)*PRODUCT(displacement**reduced)
                     reduced(axis2) = reduced(axis2) + 1
                  END DO
               END DO
            END DO
         END DO
      END DO
      DO q = 1, augmentation%nproj
         j = augmentation%radial_index(q)
         DO p = 1, q
            i = augmentation%radial_index(p)
            pair = MAX(i, j)*(MAX(i, j) - 1)/2 + MIN(i, j)
            DO l = ABS(augmentation%lbeta(i) - augmentation%lbeta(j)), &
               augmentation%lbeta(i) + augmentation%lbeta(j), 2
               c0 = augmentation%radial_value(lo, pair, l)
               c1 = h*augmentation%radial_slope(lo, pair, l)
               c2 = 3.0_dp*(augmentation%radial_value(hi, pair, l) - c0) - &
                    2.0_dp*c1 - h*augmentation%radial_slope(hi, pair, l)
               c3 = 2.0_dp*(c0 - augmentation%radial_value(hi, pair, l)) + &
                    c1 + h*augmentation%radial_slope(hi, pair, l)
               radial = c0 + t*(c1 + t*(c2 + t*c3))
               radial_d = (c1 + t*(2.0_dp*c2 + 3.0_dp*t*c3))/h
               radial_dd = (2.0_dp*c2 + 6.0_dp*t*c3)/h**2
               angular = 0.0_dp
               angular_gradient = 0.0_dp
               angular_hessian = 0.0_dp
               DO lm = l*l + 1, (l + 1)**2
                  factor = augmentation%gaunt(p, q, lm)
                  angular = angular + factor*harmonics(lm)
                  angular_gradient = angular_gradient + factor*dharmonics(:, lm)
                  angular_hessian = angular_hessian + factor*ddharmonics(:, :, lm)
               END DO
               value = radial*angular
               qr(p, q) = qr(p, q) + value
               IF (PRESENT(gradient)) THEN
                  value_gradient = 2.0_dp*displacement*radial_d*angular + radial*angular_gradient
                  gradient(p, q, :) = gradient(p, q, :) + value_gradient
               END IF
               IF (PRESENT(hessian)) THEN
                  DO axis2 = 1, 3
                     DO axis = 1, 3
                        value_hessian(axis, axis2) = &
                           4.0_dp*displacement(axis)*displacement(axis2)*radial_dd*angular + &
                           2.0_dp*radial_d*(displacement(axis)*angular_gradient(axis2) + &
                                            displacement(axis2)*angular_gradient(axis)) + &
                           radial*angular_hessian(axis, axis2)
                     END DO
                     value_hessian(axis2, axis2) = value_hessian(axis2, axis2) + 2.0_dp*radial_d*angular
                  END DO
                  hessian(p, q, :, :) = hessian(p, q, :, :) + value_hessian
               END IF
            END DO
            qr(q, p) = qr(p, q)
            IF (PRESENT(gradient)) gradient(q, p, :) = gradient(p, q, :)
            IF (PRESENT(hessian)) hessian(q, p, :, :) = hessian(p, q, :, :)
         END DO
      END DO
   END SUBROUTINE uspp_augmentation_realspace

! **************************************************************************************************
!> \brief Add one center's augmentation density and gradient to arbitrary atom-grid points.
!> \param augmentation Data with an initialized real-space cache.
!> \param displacements Cartesian quadrature points relative to this center, in bohr.
!> \param density Real symmetric on-site density matrix, including off-diagonal entries.
!> \param rho Accumulated density values; no quadrature weights included.
!> \param gradient Accumulated Cartesian density gradients.
! **************************************************************************************************
   SUBROUTINE uspp_augmentation_collocate_atom(augmentation, displacements, density, rho, gradient)
      TYPE(uspp_augmentation_type), INTENT(IN)           :: augmentation
      REAL(KIND=dp), INTENT(IN)                          :: displacements(:, :), density(:, :)
      REAL(KIND=dp), INTENT(INOUT)                       :: rho(:), gradient(:, :)

      INTEGER                                            :: axis, n, point
      REAL(KIND=dp), ALLOCATABLE                         :: dqr(:, :, :), qr(:, :)

      n = augmentation%nproj
      CPASSERT(ALL(SHAPE(density) == n))
      CPASSERT(ALL(SHAPE(displacements) == [3, SIZE(rho)]))
      CPASSERT(ALL(SHAPE(gradient) == SHAPE(displacements)))
      ALLOCATE (qr(n, n), dqr(n, n, 3))
      DO point = 1, SIZE(rho)
         CALL uspp_augmentation_realspace(augmentation, displacements(:, point), qr, dqr)
         rho(point) = rho(point) + SUM(density*qr)
         DO axis = 1, 3
            gradient(axis, point) = gradient(axis, point) + SUM(density*dqr(:, :, axis))
         END DO
      END DO
   END SUBROUTINE uspp_augmentation_collocate_atom

! **************************************************************************************************
!> \brief Adjoint of the atom-grid density map, including the GGA gradient term.
!> \param augmentation Data with an initialized real-space cache.
!> \param displacements Cartesian quadrature points relative to this center, in bohr.
!> \param weights Quadrature weights, optionally signed for a hard/soft energy difference.
!> \param v_rho Derivative of the energy density with respect to density, without weights.
!> \param v_gradient Derivative with respect to the Cartesian density gradient, without weights.
!> \param matrix Resulting symmetric on-site operator.
!> \param density Optional on-site density, required for center_gradient.
!> \param center_gradient Optional energy gradient from translating this center at fixed grid points.
!>        Projector, moving-grid and partition-weight derivatives are separate contributions.
! **************************************************************************************************
   SUBROUTINE uspp_augmentation_integrate_atom(augmentation, displacements, weights, v_rho, v_gradient, &
                                               matrix, density, center_gradient)
      TYPE(uspp_augmentation_type), INTENT(IN)           :: augmentation
      REAL(KIND=dp), INTENT(IN)                          :: displacements(:, :), weights(:), &
                                                            v_rho(:), v_gradient(:, :)
      REAL(KIND=dp), INTENT(OUT)                         :: matrix(:, :)
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: density(:, :)
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: center_gradient(3)

      INTEGER                                            :: axis, axis2, n, point
      REAL(KIND=dp)                                      :: derivative
      REAL(KIND=dp), ALLOCATABLE                         :: ddqr(:, :, :, :), dqr(:, :, :), qr(:, :)

      n = augmentation%nproj
      CPASSERT(ALL(SHAPE(matrix) == n))
      CPASSERT(ALL(SHAPE(displacements) == [3, SIZE(weights)]))
      CPASSERT(SIZE(v_rho) == SIZE(weights))
      CPASSERT(ALL(SHAPE(v_gradient) == SHAPE(displacements)))
      ALLOCATE (qr(n, n), dqr(n, n, 3))
      matrix = 0.0_dp
      IF (PRESENT(center_gradient)) THEN
         CPASSERT(PRESENT(density))
         CPASSERT(ALL(SHAPE(density) == n))
         center_gradient = 0.0_dp
         ALLOCATE (ddqr(n, n, 3, 3))
      END IF
      DO point = 1, SIZE(weights)
         IF (PRESENT(center_gradient)) THEN
            CALL uspp_augmentation_realspace(augmentation, displacements(:, point), qr, dqr, ddqr)
         ELSE
            CALL uspp_augmentation_realspace(augmentation, displacements(:, point), qr, dqr)
         END IF
         matrix(:, :) = matrix + weights(point)*v_rho(point)*qr
         DO axis = 1, 3
            matrix(:, :) = matrix + weights(point)*v_gradient(axis, point)*dqr(:, :, axis)
         END DO
         IF (PRESENT(center_gradient)) THEN
            DO axis = 1, 3
               derivative = v_rho(point)*SUM(density*dqr(:, :, axis))
               DO axis2 = 1, 3
                  derivative = derivative + v_gradient(axis2, point)*SUM(density*ddqr(:, :, axis2, axis))
               END DO
               center_gradient(axis) = center_gradient(axis) - weights(point)*derivative
            END DO
         END IF
      END DO
   END SUBROUTINE uspp_augmentation_integrate_atom

! **************************************************************************************************
!> \brief Evaluate Q_ij(G) = integral Q_ij(r-R) exp(-i G.r) dr / volume.
!>        Direct radial quadrature is a reference without reciprocal interpolation.
!> \param augmentation Prepared data.
!> \param g Cartesian reciprocal vector in inverse bohr.
!> \param position Augmentation center in bohr.
!> \param volume Cell volume in bohr cubed.
!> \param qg Complex symmetric form-factor matrix. Q(-G) = conjugate(Q(G)).
!> \param strain Optional dQ/d epsilon_ab at fixed fractional center and reciprocal index.
!>        Includes the inverse-volume derivative; the center phase is strain invariant.
! **************************************************************************************************
   SUBROUTINE uspp_augmentation_fourier(augmentation, g, position, volume, qg, strain)
      TYPE(uspp_augmentation_type), INTENT(IN)           :: augmentation
      REAL(KIND=dp), INTENT(IN)                          :: g(3), position(3), volume
      COMPLEX(KIND=dp), INTENT(OUT)                      :: qg(:, :)
      COMPLEX(KIND=dp), INTENT(OUT), OPTIONAL            :: strain(:, :, :, :)

      COMPLEX(KIND=dp)                                   :: phase, phase_l
      INTEGER                                            :: a, axis, i, j, l, lm, m, p, pair, &
                                                            powers(3), q
      REAL(KIND=dp)                                      :: angular, angular_derivative(3), &
                                                            direction(3, 1), monomial, norm_g, &
                                                            value(1)
      REAL(KIND=dp), ALLOCATABLE                         :: harmonics(:), radial(:, :), &
                                                            radial_next(:, :), &
                                                            solid_gradient(:, :), work(:)

      CPASSERT(volume > 0.0_dp)
      CPASSERT(ALL(SHAPE(qg) == augmentation%nproj))
      IF (PRESENT(strain)) THEN
         CPASSERT(ALL(SHAPE(strain) == [augmentation%nproj, augmentation%nproj, 3, 3]))
         strain = CMPLX(0.0_dp, 0.0_dp, dp)
      END IF
      norm_g = SQRT(SUM(g*g))
      IF (norm_g == 0.0_dp) THEN
         qg(:, :) = CMPLX(augmentation%q/volume, 0.0_dp, KIND=dp)
         IF (PRESENT(strain)) THEN
            DO axis = 1, 3
               strain(:, :, axis, axis) = -qg
            END DO
         END IF
         RETURN
      END IF
      direction(:, 1) = g/norm_g
      ALLOCATE (work(SIZE(augmentation%r)), harmonics((augmentation%lmax + 1)**2))
      ALLOCATE (radial(SIZE(augmentation%weighted_q, 2), 0:augmentation%lmax))
      IF (PRESENT(strain)) THEN
         ALLOCATE (solid_gradient(3, (augmentation%lmax + 1)**2), &
                   radial_next(SIZE(radial, 1), 0:augmentation%lmax))
         solid_gradient = 0.0_dp
      END IF
      DO l = 0, augmentation%lmax
         IF (PRESENT(strain)) THEN
            work(:) = augmentation%r*bessel0(norm_g*augmentation%r, l + 1)
            radial_next(:, l) = MATMUL(work, augmentation%weighted_q(:, :, l))
         END IF
         work(:) = bessel0(norm_g*augmentation%r, l)
         radial(:, l) = MATMUL(work, augmentation%weighted_q(:, :, l))
         DO m = -l, l
            lm = l*l + l + m + 1
            CALL y_lm(direction, value, l, m)
            harmonics(lm) = value(1)
            IF (PRESENT(strain)) THEN
               ! Cartesian solid-harmonic derivatives avoid polar-axis singularities.
               DO a = 1, nco(l)
                  DO axis = 1, 3
                     powers = indco(:, ncoset(l - 1) + a)
                     IF (powers(axis) == 0) CYCLE
                     monomial = REAL(powers(axis), dp)
                     powers(axis) = powers(axis) - 1
                     solid_gradient(axis, lm) = solid_gradient(axis, lm) + &
                                                orbtramat(l)%slm(l + m + 1, a)*monomial*PRODUCT(direction(:, 1)**powers)
                  END DO
               END DO
            END IF
         END DO
      END DO
      phase = EXP(CMPLX(0.0_dp, -DOT_PRODUCT(g, position), KIND=dp))*fourpi/volume
      qg = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      DO q = 1, augmentation%nproj
         j = augmentation%radial_index(q)
         DO p = 1, q
            i = augmentation%radial_index(p)
            pair = MAX(i, j)*(MAX(i, j) - 1)/2 + MIN(i, j)
            DO l = ABS(augmentation%lbeta(i) - augmentation%lbeta(j)), &
               augmentation%lbeta(i) + augmentation%lbeta(j), 2
               angular = DOT_PRODUCT(augmentation%gaunt(p, q, l*l + 1:(l + 1)**2), &
                                     harmonics(l*l + 1:(l + 1)**2))
               phase_l = CMPLX(0.0_dp, -1.0_dp, KIND=dp)**l
               qg(p, q) = qg(p, q) + phase_l*angular*radial(pair, l)
               IF (PRESENT(strain)) THEN
                  angular_derivative = MATMUL(solid_gradient(:, l*l + 1:(l + 1)**2), &
                                              augmentation%gaunt(p, q, l*l + 1:(l + 1)**2))
                  ! d[j_l(Gr)Y_l(G)]/dG combines the l/G terms analytically.
                  angular_derivative = radial(pair, l)/norm_g*angular_derivative &
                                       - radial_next(pair, l)*angular*direction(:, 1)
                  DO axis = 1, 3
                     strain(p, q, :, axis) = strain(p, q, :, axis) - phase*phase_l*g*angular_derivative(axis)
                  END DO
               END IF
            END DO
            qg(p, q) = phase*qg(p, q)
            qg(q, p) = qg(p, q)
            IF (PRESENT(strain)) THEN
               DO axis = 1, 3
                  strain(p, q, axis, axis) = strain(p, q, axis, axis) - qg(p, q)
               END DO
               strain(q, p, :, :) = strain(p, q, :, :)
            END IF
         END DO
      END DO
   END SUBROUTINE uspp_augmentation_fourier

! **************************************************************************************************
!> \brief Add the augmentation density on a list of reciprocal vectors.
!> \param augmentation Prepared data.
!> \param gvecs Cartesian reciprocal vectors (3,ng).
!> \param position Center in bohr.
!> \param volume Cell volume.
!> \param density_matrix Real symmetric projector density matrix, with occupations included.
!> \param rho Fourier density coefficients (integral divided by volume), accumulated in place.
! **************************************************************************************************
   SUBROUTINE uspp_augmentation_collocate(augmentation, gvecs, position, volume, density_matrix, rho)
      TYPE(uspp_augmentation_type), INTENT(IN)           :: augmentation
      REAL(KIND=dp), INTENT(IN)                          :: gvecs(:, :), position(3), volume, &
                                                            density_matrix(:, :)
      COMPLEX(KIND=dp), INTENT(INOUT)                    :: rho(:)

      COMPLEX(KIND=dp), ALLOCATABLE                      :: qg(:, :)
      INTEGER                                            :: ig

      CPASSERT(SIZE(gvecs, 1) == 3 .AND. SIZE(gvecs, 2) == SIZE(rho))
      CPASSERT(ALL(SHAPE(density_matrix) == augmentation%nproj))
      ALLOCATE (qg(augmentation%nproj, augmentation%nproj))
      DO ig = 1, SIZE(rho)
         CALL uspp_augmentation_fourier(augmentation, gvecs(:, ig), position, volume, qg)
         rho(ig) = rho(ig) + SUM(density_matrix*qg)
      END DO
   END SUBROUTINE uspp_augmentation_collocate

! **************************************************************************************************
!> \brief Adjoint of augmentation collocation: integrate a real local potential against Q.
!> \param augmentation Prepared data.
!> \param gvecs Cartesian reciprocal vectors (3,ng).
!> \param position Center in bohr.
!> \param volume Cell volume.
!> \param vg Fourier potential coefficients (integral divided by volume).
!> \param multiplicity One for a complete reciprocal list; two for a retained +/-G pair.
!> \param matrix Real symmetric screened projector contribution, overwritten.
!> \param derivative Optional derivatives of matrix with respect to the augmentation center.
!> \param strain Optional density strain contribution at fixed Fourier potential and volume weight.
! **************************************************************************************************
   SUBROUTINE uspp_augmentation_integrate(augmentation, gvecs, position, volume, vg, multiplicity, matrix, derivative, strain)
      TYPE(uspp_augmentation_type), INTENT(IN)           :: augmentation
      REAL(KIND=dp), INTENT(IN)                          :: gvecs(:, :), position(3), volume
      COMPLEX(KIND=dp), INTENT(IN)                       :: vg(:)
      REAL(KIND=dp), INTENT(IN)                          :: multiplicity(:)
      REAL(KIND=dp), INTENT(OUT)                         :: matrix(:, :)
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: derivative(:, :, :), strain(:, :, :, :)

      COMPLEX(KIND=dp), ALLOCATABLE                      :: dq(:, :, :, :), qg(:, :)
      INTEGER                                            :: idir, ig, jdir

      CPASSERT(SIZE(gvecs, 1) == 3 .AND. SIZE(gvecs, 2) == SIZE(vg))
      CPASSERT(SIZE(multiplicity) == SIZE(vg))
      CPASSERT(ALL(SHAPE(matrix) == augmentation%nproj))
      ALLOCATE (qg(augmentation%nproj, augmentation%nproj))
      matrix = 0.0_dp
      IF (PRESENT(derivative)) THEN
         CPASSERT(ALL(SHAPE(derivative) == [augmentation%nproj, augmentation%nproj, 3]))
         derivative = 0.0_dp
      END IF
      IF (PRESENT(strain)) THEN
         CPASSERT(ALL(SHAPE(strain) == [augmentation%nproj, augmentation%nproj, 3, 3]))
         ALLOCATE (dq(augmentation%nproj, augmentation%nproj, 3, 3))
         strain = 0.0_dp
      END IF
      DO ig = 1, SIZE(vg)
         CALL uspp_augmentation_fourier(augmentation, gvecs(:, ig), position, volume, qg, dq)
         IF (PRESENT(strain)) THEN
            DO jdir = 1, 3
               DO idir = 1, 3
                  strain(:, :, idir, jdir) = strain(:, :, idir, jdir) + volume*multiplicity(ig)* &
                                             REAL(CONJG(dq(:, :, idir, jdir))*vg(ig), KIND=dp)
               END DO
            END DO
         END IF
         matrix(:, :) = matrix + volume*multiplicity(ig)*REAL(CONJG(qg)*vg(ig), KIND=dp)
         IF (PRESENT(derivative)) THEN
            DO idir = 1, 3
               derivative(:, :, idir) = derivative(:, :, idir) + volume*multiplicity(ig)* &
                                        REAL(CMPLX(0.0_dp, gvecs(idir, ig), dp)*CONJG(qg)*vg(ig), KIND=dp)
            END DO
         END IF
      END DO
   END SUBROUTINE uspp_augmentation_integrate

END MODULE uspp_augmentation
