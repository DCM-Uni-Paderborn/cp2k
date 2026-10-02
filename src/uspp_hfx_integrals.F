!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Native scalar USPP transition-charge integrals for isolated Coulomb RI-HFX.
!>        Integrates complete Q_ij against Gaussian Coulomb potentials, with both center derivatives.
!>        No GAPW compensation subtraction or fitting of the original augmentation is performed.
! **************************************************************************************************
MODULE uspp_hfx_integrals
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
   USE gapw_coulomb_integrals,          ONLY: gapw_coulomb_init,&
                                              gaussian_coulomb_radial
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: dfac,&
                                              fourpi
   USE orbital_pointers,                ONLY: indco,&
                                              nco,&
                                              ncoset
   USE orbital_transformation_matrices, ONLY: orbtramat
   USE spherical_harmonics,             ONLY: y_lm
   USE uspp_augmentation,               ONLY: uspp_augmentation_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'uspp_hfx_integrals'

   PUBLIC :: uspp_hfx_integrals_init, uspp_hfx_primitive_integrals, uspp_hfx_auxiliary_integrals

CONTAINS

! **************************************************************************************************
!> \brief Initialize shared angular and Boys tables before parallel integral evaluation.
!> \param maxl Largest Gaussian Cartesian degree or augmentation angular momentum.
! **************************************************************************************************
   SUBROUTINE uspp_hfx_integrals_init(maxl)
      INTEGER, INTENT(IN)                                :: maxl

      CPASSERT(maxl >= 0)
      CALL gapw_coulomb_init(maxl, maxl/2)
   END SUBROUTINE uspp_hfx_integrals_init

! **************************************************************************************************
!> \brief (Gaussian multipole | 1/r | Q_ij) on the original augmentation quadrature.
!> \param augmentation Native scalar Q radial/angular data, in magnetic-projector ordering.
!> \param alpha Positive Gaussian exponent, in bohr**(-2).
!> \param power Nonnegative additional radial power of r**2.
!> \param l Gaussian solid-harmonic angular momentum.
!> \param displacement Augmentation center minus Gaussian center, in bohr.
!> \param directions Unit vectors of the angular quadrature.
!> \param weights Angular weights, summing to 4*pi.
!> \param integrals Full projector-pair matrices, m=-l..l. Gaussian r**(2*power+l)Y_lm exp(-alpha*r**2)
!>        uses unit-normalized real Y_lm and includes no radial Gaussian normalization.
!> \param gaussian_derivative Optional d(integrals)/d(Gaussian center), Cartesian axis last.
!> \param augmentation_derivative Optional d(integrals)/d(Q center), Cartesian axis last.
! **************************************************************************************************
   SUBROUTINE uspp_hfx_primitive_integrals(augmentation, alpha, power, l, displacement, directions, weights, &
                                           integrals, gaussian_derivative, augmentation_derivative)
      TYPE(uspp_augmentation_type), INTENT(IN)           :: augmentation
      REAL(dp), INTENT(IN)                               :: alpha
      INTEGER, INTENT(IN)                                :: power, l
      REAL(dp), INTENT(IN)                               :: displacement(3), directions(:, :), &
                                                            weights(:)
      REAL(dp), INTENT(OUT)                              :: integrals(:, :, :)
      REAL(dp), INTENT(OUT), OPTIONAL :: gaussian_derivative(:, :, :, :), &
         augmentation_derivative(:, :, :, :)

      INTEGER                                            :: axis, i, ia, ir, j, la, lm, m, n, p, &
                                                            pair, q
      LOGICAL                                            :: derivatives, valid
      REAL(dp)                                           :: gradient(3, 2*l + 1), point(3), &
                                                            potential(2*l + 1), radial
      REAL(dp), ALLOCATABLE                              :: d_integrals(:, :, :, :), &
                                                            d_projected(:, :, :), harmonics(:, :), &
                                                            projected(:, :), values(:)

      valid = alpha > 0.0_dp .AND. ieee_is_finite(alpha) .AND. MIN(l, power) >= 0
      CPASSERT(valid)
      valid = ALL(ieee_is_finite(displacement))
      CPASSERT(valid)
      CPASSERT(ASSOCIATED(orbtramat))
      CPASSERT(UBOUND(orbtramat, 1) >= l)
      n = augmentation%nproj
      valid = n > 0 .AND. augmentation%lmax >= 0
      CPASSERT(valid)
      valid = ALL(SHAPE(integrals) == [n, n, 2*l + 1])
      CPASSERT(valid)
      valid = SIZE(directions, 1) == 3 .AND. SIZE(directions, 2) == SIZE(weights)
      CPASSERT(valid)
      valid = SIZE(weights) > 0 .AND. ALL(ieee_is_finite(weights)) .AND. ALL(ieee_is_finite(directions))
      CPASSERT(valid)
      valid = ABS(SUM(weights) - fourpi) < 1.e-10_dp
      CPASSERT(valid)
      valid = MAXVAL(ABS(SUM(directions**2, DIM=1) - 1.0_dp)) < 1.e-12_dp
      CPASSERT(valid)
      derivatives = PRESENT(gaussian_derivative) .OR. PRESENT(augmentation_derivative)
      IF (PRESENT(gaussian_derivative)) THEN
         valid = ALL(SHAPE(gaussian_derivative) == [n, n, 2*l + 1, 3])
         CPASSERT(valid)
      END IF
      IF (PRESENT(augmentation_derivative)) THEN
         valid = ALL(SHAPE(augmentation_derivative) == [n, n, 2*l + 1, 3])
         CPASSERT(valid)
      END IF
      ALLOCATE (harmonics(SIZE(weights), (augmentation%lmax + 1)**2), values(SIZE(weights)), &
                projected((augmentation%lmax + 1)**2, 2*l + 1))
      IF (derivatives) ALLOCATE (d_projected((augmentation%lmax + 1)**2, 2*l + 1, 3), &
                                 d_integrals(n, n, 2*l + 1, 3), SOURCE=0.0_dp)
      DO la = 0, augmentation%lmax
         DO m = -la, la
            lm = la*la + la + m + 1
            CALL y_lm(directions, values, la, m)
            harmonics(:, lm) = weights*values
         END DO
      END DO
      integrals = 0.0_dp
      DO ir = 1, SIZE(augmentation%weighted_q, 1)
         IF (ALL(augmentation%weighted_q(ir, :, :) == 0.0_dp)) CYCLE
         projected = 0.0_dp
         IF (derivatives) d_projected = 0.0_dp
         DO ia = 1, SIZE(weights)
            point = displacement + augmentation%r(ir)*directions(:, ia)
            CALL gaussian_potential(alpha, power, l, point, potential, gradient)
            DO m = 1, 2*l + 1
               projected(:, m) = projected(:, m) + harmonics(ia, :)*potential(m)
               IF (derivatives) THEN
                  DO axis = 1, 3
                     d_projected(:, m, axis) = d_projected(:, m, axis) + harmonics(ia, :)*gradient(axis, m)
                  END DO
               END IF
            END DO
         END DO
         DO q = 1, n
            j = augmentation%radial_index(q)
            DO p = 1, q
               i = augmentation%radial_index(p)
               pair = MAX(i, j)*(MAX(i, j) - 1)/2 + MIN(i, j)
               DO la = ABS(augmentation%lbeta(i) - augmentation%lbeta(j)), &
                  augmentation%lbeta(i) + augmentation%lbeta(j), 2
                  radial = augmentation%weighted_q(ir, pair, la)
                  IF (radial == 0.0_dp) CYCLE
                  DO lm = la*la + 1, (la + 1)**2
                     IF (augmentation%gaunt(p, q, lm) == 0.0_dp) CYCLE
                     integrals(p, q, :) = integrals(p, q, :) + &
                                          radial*augmentation%gaunt(p, q, lm)*projected(lm, :)
                     IF (derivatives) d_integrals(p, q, :, :) = d_integrals(p, q, :, :) + &
                                                                radial*augmentation%gaunt(p, q, lm)*d_projected(lm, :, :)
                  END DO
               END DO
            END DO
         END DO
      END DO
      DO q = 1, n
         DO p = 1, q - 1
            integrals(q, p, :) = integrals(p, q, :)
            IF (derivatives) d_integrals(q, p, :, :) = d_integrals(p, q, :, :)
         END DO
      END DO
      IF (PRESENT(gaussian_derivative)) gaussian_derivative = -d_integrals
      IF (PRESENT(augmentation_derivative)) augmentation_derivative = d_integrals
   END SUBROUTINE uspp_hfx_primitive_integrals

! **************************************************************************************************
!> \brief Native Q coupling to a complete contracted Gaussian auxiliary shell set.
!>        Uses the same Cartesian primitive ordering and sphi normalization as Gaussian RI/SAP.
!>        Homogeneous Cartesian monomials are expanded analytically in r**(2*power) solid harmonics;
!>        trace components are retained, so mixed angular-momentum sets need no special treatment.
!> \param augmentation Native scalar transition-charge data.
!> \param basis Ordinary RI Gaussian basis, with its existing normalized sphi coefficients.
!> \param iset Shell-set index in basis.
!> \param gaussian_center Gaussian RI atom center, in bohr.
!> \param augmentation_center Native Q atom center, in bohr.
!> \param directions Angular quadrature unit vectors.
!> \param weights Angular weights summing to 4*pi.
!> \param integrals (RI auxiliary | Coulomb | Q), shape (nsgf_set(iset), nproj, nproj).
!> \param gaussian_derivative Optional Cartesian center derivatives, axis first.
!> \param augmentation_derivative Optional Cartesian Q-center derivatives, axis first.
! **************************************************************************************************
   SUBROUTINE uspp_hfx_auxiliary_integrals(augmentation, basis, iset, gaussian_center, augmentation_center, &
                                           directions, weights, integrals, gaussian_derivative, augmentation_derivative)
      TYPE(uspp_augmentation_type), INTENT(IN)           :: augmentation
      TYPE(gto_basis_set_type), INTENT(IN)               :: basis
      INTEGER, INTENT(IN)                                :: iset
      REAL(dp), INTENT(IN)                               :: gaussian_center(3), &
                                                            augmentation_center(3), &
                                                            directions(:, :), weights(:)
      REAL(dp), INTENT(OUT)                              :: integrals(:, :, :)
      REAL(dp), INTENT(OUT), OPTIONAL :: gaussian_derivative(:, :, :, :), &
         augmentation_derivative(:, :, :, :)

      INTEGER                                            :: axis, degree, first, ip, l, last, m, n, &
                                                            naux, nc, nset, row, t
      INTEGER, POINTER                                   :: lmax(:), npgf(:), nsgf_set(:)
      LOGICAL                                            :: derivatives, valid
      REAL(dp)                                           :: displacement(3)
      REAL(dp), ALLOCATABLE :: coefficients(:, :), primitive(:, :, :), &
         primitive_derivative(:, :, :, :)
      REAL(dp), POINTER                                  :: sphi(:, :), zet(:, :)

      CALL get_gto_basis_set(basis, nset=nset, lmax=lmax, npgf=npgf, nsgf_set=nsgf_set, sphi=sphi, zet=zet)
      valid = iset >= 1 .AND. iset <= nset
      CPASSERT(valid)
      valid = ASSOCIATED(lmax) .AND. ASSOCIATED(npgf) .AND. ASSOCIATED(nsgf_set)
      CPASSERT(valid)
      valid = ASSOCIATED(sphi) .AND. ASSOCIATED(zet)
      CPASSERT(valid)
      valid = lmax(iset) >= 0 .AND. npgf(iset) > 0
      CPASSERT(valid)
      CPASSERT(ASSOCIATED(orbtramat))
      CPASSERT(UBOUND(orbtramat, 1) >= lmax(iset))
      n = augmentation%nproj
      naux = nsgf_set(iset)
      nc = ncoset(lmax(iset))
      first = SUM(nsgf_set(:iset - 1)) + 1
      last = first + naux - 1
      valid = naux > 0 .AND. ALL(SHAPE(integrals) == [naux, n, n])
      CPASSERT(valid)
      valid = SIZE(sphi, 1) >= nc*npgf(iset) .AND. SIZE(sphi, 2) >= last
      CPASSERT(valid)
      valid = SIZE(zet, 1) >= npgf(iset) .AND. SIZE(zet, 2) >= iset
      CPASSERT(valid)
      valid = ALL(ieee_is_finite(sphi(:nc*npgf(iset), first:last)))
      CPASSERT(valid)
      derivatives = PRESENT(gaussian_derivative) .OR. PRESENT(augmentation_derivative)
      IF (PRESENT(gaussian_derivative)) THEN
         valid = ALL(SHAPE(gaussian_derivative) == [3, naux, n, n])
         CPASSERT(valid)
      END IF
      IF (PRESENT(augmentation_derivative)) THEN
         valid = ALL(SHAPE(augmentation_derivative) == [3, naux, n, n])
         CPASSERT(valid)
      END IF
      displacement = augmentation_center - gaussian_center
      integrals = 0.0_dp
      IF (PRESENT(gaussian_derivative)) gaussian_derivative = 0.0_dp
      IF (PRESENT(augmentation_derivative)) augmentation_derivative = 0.0_dp
      DO ip = 1, npgf(iset)
         DO degree = 0, lmax(iset)
            DO l = MODULO(degree, 2), degree, 2
               ALLOCATE (coefficients(2*l + 1, naux), SOURCE=0.0_dp)
               DO t = ncoset(degree - 1) + 1, ncoset(degree)
                  row = (ip - 1)*nc + t
                  IF (ALL(sphi(row, first:last) == 0.0_dp)) CYCLE
                  DO m = 1, 2*l + 1
                     coefficients(m, :) = coefficients(m, :) + &
                                          monomial_harmonic(indco(:, t), l, m)*sphi(row, first:last)
                  END DO
               END DO
               IF (ANY(coefficients /= 0.0_dp)) THEN
                  ALLOCATE (primitive(n, n, 2*l + 1))
                  IF (derivatives) THEN
                     ALLOCATE (primitive_derivative(n, n, 2*l + 1, 3))
                     CALL uspp_hfx_primitive_integrals(augmentation, zet(ip, iset), (degree - l)/2, l, &
                                                       displacement, directions, weights, primitive, primitive_derivative)
                  ELSE
                     CALL uspp_hfx_primitive_integrals(augmentation, zet(ip, iset), (degree - l)/2, l, &
                                                       displacement, directions, weights, primitive)
                  END IF
                  DO t = 1, naux
                     DO m = 1, 2*l + 1
                        integrals(t, :, :) = integrals(t, :, :) + coefficients(m, t)*primitive(:, :, m)
                        IF (derivatives) THEN
                           DO axis = 1, 3
                              IF (PRESENT(gaussian_derivative)) gaussian_derivative(axis, t, :, :) = &
                                 gaussian_derivative(axis, t, :, :) + coefficients(m, t)*primitive_derivative(:, :, m, axis)
                              IF (PRESENT(augmentation_derivative)) augmentation_derivative(axis, t, :, :) = &
                                 augmentation_derivative(axis, t, :, :) - coefficients(m, t)*primitive_derivative(:, :, m, axis)
                           END DO
                        END IF
                     END DO
                  END DO
                  DEALLOCATE (primitive)
                  IF (derivatives) DEALLOCATE (primitive_derivative)
               END IF
               DEALLOCATE (coefficients)
            END DO
         END DO
      END DO
   END SUBROUTINE uspp_hfx_auxiliary_integrals

! **************************************************************************************************
!> \brief Exact projection of a Cartesian monomial onto unit-normalized real Y_lm on the sphere.
!>        The sphere moments retain all same-parity angular channels, including Cartesian traces.
!> \param powers Cartesian monomial powers.
!> \param l Angular momentum.
!> \param m CP2K harmonic index, 1..2*l+1.
!> \return Coefficient of r**(SUM(powers)-l) times the solid harmonic.
! **************************************************************************************************
   FUNCTION monomial_harmonic(powers, l, m) RESULT(coefficient)
      INTEGER, INTENT(IN)                                :: powers(3), l, m
      REAL(dp)                                           :: coefficient

      INTEGER                                            :: a, combined(3)

      coefficient = 0.0_dp
      DO a = 1, nco(l)
         combined = powers + indco(:, ncoset(l - 1) + a)
         IF (ANY(MODULO(combined, 2) /= 0)) CYCLE
         coefficient = coefficient + orbtramat(l)%slm(m, a)*fourpi* &
                       PRODUCT(dfac(combined - 1))/dfac(SUM(combined) + 1)
      END DO
   END FUNCTION monomial_harmonic

! **************************************************************************************************
!> \brief Full isolated Gaussian Coulomb potential and Cartesian evaluation-point gradient.
!>        In particular the asymptotic multipole is retained; the compensated tail is not used.
! **************************************************************************************************
   SUBROUTINE gaussian_potential(alpha, power, l, displacement, potential, gradient)
      REAL(dp), INTENT(IN)                               :: alpha
      INTEGER, INTENT(IN)                                :: power, l
      REAL(dp), INTENT(IN)                               :: displacement(3)
      REAL(dp), INTENT(OUT)                              :: potential(2*l + 1), gradient(3, 2*l + 1)

      INTEGER                                            :: a, axis, m, powers(3), reduced(3)
      REAL(dp)                                           :: angular, dangular(3), ddh, dh, factor, h

      CALL gaussian_coulomb_radial(alpha, power, l, SUM(displacement**2), .FALSE., h, dh, ddh)
      DO m = 1, 2*l + 1
         angular = 0.0_dp
         dangular = 0.0_dp
         DO a = 1, nco(l)
            powers = indco(:, ncoset(l - 1) + a)
            factor = orbtramat(l)%slm(m, a)
            IF (factor == 0.0_dp) CYCLE
            angular = angular + factor*PRODUCT(displacement**powers)
            DO axis = 1, 3
               IF (powers(axis) == 0) CYCLE
               reduced = powers
               reduced(axis) = reduced(axis) - 1
               dangular(axis) = dangular(axis) + factor*REAL(powers(axis), dp)*PRODUCT(displacement**reduced)
            END DO
         END DO
         potential(m) = h*angular
         gradient(:, m) = 2.0_dp*displacement*dh*angular + h*dangular
      END DO
   END SUBROUTINE gaussian_potential

END MODULE uspp_hfx_integrals
