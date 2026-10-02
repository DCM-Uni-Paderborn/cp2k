!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Original-data PAW frozen-core kinetic fields shared by smooth and one-center XC.
!>        The PS field is an auxiliary XC argument, not a physical pseudo core orbital.
! **************************************************************************************************
MODULE paw_core_kinetic
   USE atom_upf,                        ONLY: atom_upfpot_type,&
                                              upf_radial_weights
   USE bessel_lib,                      ONLY: bessel0
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: fourpi
   USE periodic_table,                  ONLY: get_ptable_info
   USE splines,                         ONLY: spline3ders
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'paw_core_kinetic'

   TYPE paw_core_kinetic_type
      LOGICAL :: available = .FALSE.
      CHARACTER(LEN=256) :: failure_reason = 'PAW core kinetic fields have not been initialized'
      REAL(dp) :: smoothing_radius = 0.0_dp, taper_radius = 0.0_dp, support_radius = 0.0_dp, anchor = 0.0_dp
      ! The AE samples retain the complete original mesh and amplitudes.
      REAL(dp), ALLOCATABLE :: radial_mesh(:), kinetic(:, :), weighted(:)
      ! Signed occupied R=u/r waves and their spline slopes use positive radii only.
      REAL(dp), ALLOCATABLE :: r(:), wave(:, :), slope(:, :)
      INTEGER, ALLOCATABLE :: angular(:), occupations(:)
   END TYPE paw_core_kinetic_type

   PUBLIC :: paw_core_kinetic_type, paw_core_kinetic_init, paw_core_kinetic_build_ae
   PUBLIC :: paw_core_kinetic_values, paw_core_kinetic_cartesian, paw_core_kinetic_fourier

CONTAINS

! **************************************************************************************************
!> \brief Cache a verified AE core and one explicit smooth, nonnegative auxiliary PS kinetic field.
!>        Inside the original augmentation radius, a quintic switch in (r/rc)^2 blends a
!>        Gaussian anchored at tau_AE(rc) into the original signed-orbital spline field.
!>        The final two mesh intervals carry a C2 taper; the full AE tail is retained.
!>        All consumers must use this same auxiliary definition and its derivatives.
!> \param pot Original scalar PAW UPF data.
!> \param cache Immutable radial cache; unavailable original cores retain a diagnostic.
! **************************************************************************************************
   SUBROUTINE paw_core_kinetic_init(pot, cache)
      TYPE(atom_upfpot_type), INTENT(IN)                 :: pot
      TYPE(paw_core_kinetic_type), INTENT(OUT)           :: cache

      INTEGER                                            :: first, i, j, nr, nshell
      INTEGER, ALLOCATABLE                               :: selected(:)
      LOGICAL                                            :: available
      REAL(dp)                                           :: derivative
      REAL(dp), ALLOCATABLE                              :: tau_ae(:)

      CALL paw_core_kinetic_build_ae(pot, tau_ae, available, cache%failure_reason, selected)
      IF (.NOT. available) RETURN
      nr = SIZE(pot%r)
      IF (SUM(selected) == 0) THEN
         cache%failure_reason = 'Core-free PAW cache contain nonzero PS core density'
         IF (ALLOCATED(pot%rho_nlcc)) THEN
            IF (ANY(pot%rho_nlcc /= 0.0_dp)) RETURN
         END IF
         ALLOCATE (cache%radial_mesh, SOURCE=pot%r)
         ALLOCATE (cache%kinetic(nr, 2), cache%weighted(nr), SOURCE=0.0_dp)
         cache%available = .TRUE.
         cache%failure_reason = ''
         RETURN
      END IF
      first = 1
      IF (pot%r(1) == 0.0_dp) first = 2
      cache%failure_reason = 'PS core tau needs an augmentation radius and two outer taper intervals'
      IF (nr - first + 1 < 4) RETURN
      IF (pot%paw%cutoff_index < first .OR. pot%paw%cutoff_index > nr - 2) RETURN
      cache%smoothing_radius = pot%r(pot%paw%cutoff_index)
      IF (cache%smoothing_radius <= 0.0_dp) RETURN
      cache%taper_radius = pot%r(nr - 2)
      cache%support_radius = pot%r(nr)
      nshell = COUNT(selected > 0)
      ALLOCATE (cache%radial_mesh, SOURCE=pot%r)
      ALLOCATE (cache%r, SOURCE=pot%r(first:))
      ALLOCATE (cache%wave(SIZE(cache%r), nshell), cache%slope(SIZE(cache%r), nshell), &
                cache%angular(nshell), cache%occupations(nshell), cache%kinetic(nr, 2))
      j = 0
      DO i = 1, SIZE(selected)
         IF (selected(i) == 0) CYCLE
         j = j + 1
         cache%wave(:, j) = pot%paw%gipaw_core_orbitals%wfc(first:, i)/cache%r
         CALL spline3ders(cache%r, cache%wave(:, j), cache%r, dynew=cache%slope(:, j))
         cache%angular(j) = pot%paw%gipaw_core_orbitals%angular_momenta(i)
         cache%occupations(j) = selected(i)
      END DO
      cache%kinetic(:, 1) = tau_ae
      CALL core_ae_value(cache, cache%smoothing_radius, cache%anchor, derivative)
      cache%available = .TRUE.
      CALL paw_core_kinetic_values(cache, cache%radial_mesh, cache%kinetic(:, 2))
      ALLOCATE (cache%weighted, SOURCE=fourpi*upf_radial_weights(pot)*pot%r**2*cache%kinetic(:, 2))
      cache%failure_reason = ''
   END SUBROUTINE paw_core_kinetic_init

! **************************************************************************************************
!> \brief Evaluate the shared PS kinetic field and radial derivative at arbitrary nonnegative radii.
!> \param cache Initialized original-data cache.
!> \param radii Evaluation radii, including the physical smooth origin.
!> \param tau Spin-summed kinetic density including the factor 1/2.
!> \param derivative Optional derivative with respect to radius.
! **************************************************************************************************
   SUBROUTINE paw_core_kinetic_values(cache, radii, tau, derivative)
      TYPE(paw_core_kinetic_type), INTENT(IN)            :: cache
      REAL(dp), INTENT(IN)                               :: radii(:)
      REAL(dp), INTENT(OUT)                              :: tau(:)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: derivative(:)

      INTEGER                                            :: i
      REAL(dp)                                           :: slope

      CPASSERT(cache%available)
      CPASSERT(SIZE(tau) == SIZE(radii))
      CPASSERT(ALL(ieee_is_finite(radii)) .AND. ALL(radii >= 0.0_dp))
      IF (PRESENT(derivative)) THEN
         CPASSERT(SIZE(derivative) == SIZE(radii))
      END IF
      DO i = 1, SIZE(radii)
         CALL core_ps_value(cache, radii(i), tau(i), slope)
         IF (PRESENT(derivative)) derivative(i) = slope
      END DO
   END SUBROUTINE paw_core_kinetic_values

! **************************************************************************************************
!> \brief Cartesian field and spatial gradient; center derivatives have the opposite sign.
!> \param cache Initialized original-data cache.
!> \param displacement Evaluation point minus the core center.
!> \param tau Spin-summed PS kinetic field.
!> \param gradient Optional spatial gradient, with the physical zero at the smooth origin.
! **************************************************************************************************
   SUBROUTINE paw_core_kinetic_cartesian(cache, displacement, tau, gradient)
      TYPE(paw_core_kinetic_type), INTENT(IN)            :: cache
      REAL(dp), INTENT(IN)                               :: displacement(3)
      REAL(dp), INTENT(OUT)                              :: tau
      REAL(dp), INTENT(OUT), OPTIONAL                    :: gradient(3)

      REAL(dp)                                           :: radius, slope

      CPASSERT(cache%available)
      CPASSERT(ALL(ieee_is_finite(displacement)))
      radius = SQRT(SUM(displacement**2))
      CALL core_ps_value(cache, radius, tau, slope)
      IF (PRESENT(gradient)) THEN
         gradient = 0.0_dp
         IF (radius > 0.0_dp) gradient = slope*displacement/radius
      END IF
   END SUBROUTINE paw_core_kinetic_cartesian

! **************************************************************************************************
!> \brief Fourier field and exact adjoints of its original-mesh discrete radial transform.
!>        Under affine strain g transforms contravariantly, position covariantly, and volume
!>        with the determinant; the center phase is invariant in the strain derivative.
!> \param cache Initialized original-data cache.
!> \param g Cartesian reciprocal vector, including zero or arbitrary complex-wave sampling.
!> \param position Cartesian core center.
!> \param volume Cell volume.
!> \param coefficient exp(-i g.R) times the radial transform divided by volume.
!> \param strain Optional derivative with respect to affine strain epsilon_ab.
!> \param center_derivative Optional derivative with respect to the core center at fixed g.
! **************************************************************************************************
   SUBROUTINE paw_core_kinetic_fourier(cache, g, position, volume, coefficient, strain, center_derivative)
      TYPE(paw_core_kinetic_type), INTENT(IN)            :: cache
      REAL(dp), INTENT(IN)                               :: g(3), position(3), volume
      COMPLEX(dp), INTENT(OUT)                           :: coefficient
      COMPLEX(dp), INTENT(OUT), OPTIONAL                 :: strain(3, 3), center_derivative(3)

      COMPLEX(dp)                                        :: phase
      INTEGER                                            :: a, b
      REAL(dp)                                           :: magnitude, next, radial

      CPASSERT(cache%available .AND. volume > 0.0_dp)
      CPASSERT(ALL(ieee_is_finite(g)))
      CPASSERT(ALL(ieee_is_finite(position)))
      CPASSERT(ieee_is_finite(volume))
      magnitude = SQRT(SUM(g**2))
      radial = SUM(cache%weighted*bessel0(magnitude*cache%radial_mesh, 0))
      phase = EXP(CMPLX(0.0_dp, -DOT_PRODUCT(g, position), KIND=dp))/volume
      coefficient = phase*radial
      IF (PRESENT(center_derivative)) center_derivative = CMPLX(0.0_dp, -g, KIND=dp)*coefficient
      IF (PRESENT(strain)) THEN
         strain = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
         IF (magnitude > 0.0_dp) THEN
            next = SUM(cache%weighted*cache%radial_mesh*bessel0(magnitude*cache%radial_mesh, 1))
            DO b = 1, 3
               DO a = 1, 3
                  strain(a, b) = phase*g(a)*g(b)*next/magnitude
               END DO
            END DO
         END IF
         DO a = 1, 3
            strain(a, a) = strain(a, a) - coefficient
         END DO
      END IF
   END SUBROUTINE paw_core_kinetic_fourier

! **************************************************************************************************
!> \brief Shared smooth auxiliary profile, evaluated without a singular angular quotient at zero.
!> \param cache Immutable original-data cache.
!> \param radius Nonnegative evaluation radius.
!> \param value PS kinetic field.
!> \param derivative Radial derivative.
! **************************************************************************************************
   SUBROUTINE core_ps_value(cache, radius, value, derivative)
      TYPE(paw_core_kinetic_type), INTENT(IN)            :: cache
      REAL(dp), INTENT(IN)                               :: radius
      REAL(dp), INTENT(OUT)                              :: value, derivative

      INTEGER                                            :: i
      REAL(dp)                                           :: angular, blend, complement, dblend, &
                                                            dquotient, dwindow, factor, gaussian, &
                                                            h, quotient, rsecond, rslope, rvalue, &
                                                            t, u, window

      value = 0.0_dp
      derivative = 0.0_dp
      IF (.NOT. ALLOCATED(cache%wave)) RETURN
      IF (radius >= cache%support_radius) RETURN
      IF (radius < cache%smoothing_radius) THEN
         u = (radius/cache%smoothing_radius)**2
         h = 10.0_dp - 15.0_dp*u + 6.0_dp*u**2
         blend = u**3*h
         dblend = 60.0_dp*radius/cache%smoothing_radius**2*u**2*(1.0_dp - u)**2
         ! S(u)/r^2 and its derivative are finite polynomials at r=0.
         quotient = u**2*h/cache%smoothing_radius**2
         dquotient = 2.0_dp*radius/cache%smoothing_radius**4*(2.0_dp*u*h + u**2*(-15.0_dp + 12.0_dp*u))
         complement = (1.0_dp - u)**3*(1.0_dp + 3.0_dp*u + 6.0_dp*u**2)
         gaussian = cache%anchor*EXP(1.0_dp - u)
         value = complement*gaussian
         derivative = -dblend*gaussian - complement*2.0_dp*radius/cache%smoothing_radius**2*gaussian
         DO i = 1, SIZE(cache%occupations)
            CALL core_wave_value(cache, i, radius, rvalue, rslope, rsecond)
            factor = 0.5_dp*REAL(cache%occupations(i), dp)/fourpi
            angular = REAL(cache%angular(i), dp)*REAL(cache%angular(i) + 1, dp)
            value = value + factor*(blend*rslope**2 + angular*quotient*rvalue**2)
            derivative = derivative + factor*(dblend*rslope**2 + 2.0_dp*blend*rslope*rsecond + &
                                              angular*(dquotient*rvalue**2 + 2.0_dp*quotient*rvalue*rslope))
         END DO
      ELSE
         CALL core_ae_value(cache, radius, value, derivative)
      END IF
      IF (radius > cache%taper_radius) THEN
         t = (radius - cache%taper_radius)/(cache%support_radius - cache%taper_radius)
         window = (1.0_dp - t)**3*(1.0_dp + 3.0_dp*t + 6.0_dp*t**2)
         dwindow = -30.0_dp*t**2*(1.0_dp - t)**2/(cache%support_radius - cache%taper_radius)
         derivative = window*derivative + dwindow*value
         value = window*value
      END IF
   END SUBROUTINE core_ps_value

! **************************************************************************************************
!> \brief Positive-radius AE field from the unchanged occupied signed-orbital splines.
!> \param cache Original-data cache.
!> \param radius Positive evaluation radius.
!> \param value AE kinetic field.
!> \param derivative Radial derivative.
! **************************************************************************************************
   SUBROUTINE core_ae_value(cache, radius, value, derivative)
      TYPE(paw_core_kinetic_type), INTENT(IN)            :: cache
      REAL(dp), INTENT(IN)                               :: radius
      REAL(dp), INTENT(OUT)                              :: value, derivative

      INTEGER                                            :: i
      REAL(dp)                                           :: angular, factor, radial, second, slope

      CPASSERT(radius > 0.0_dp)
      value = 0.0_dp
      derivative = 0.0_dp
      DO i = 1, SIZE(cache%occupations)
         CALL core_wave_value(cache, i, radius, radial, slope, second)
         factor = 0.5_dp*REAL(cache%occupations(i), dp)/fourpi
         angular = REAL(cache%angular(i), dp)*REAL(cache%angular(i) + 1, dp)
         value = value + factor*(slope**2 + angular*(radial/radius)**2)
         derivative = derivative + factor*(2.0_dp*slope*second + &
                                           angular*(2.0_dp*radial*slope/radius**2 - 2.0_dp*radial**2/radius**3))
      END DO
   END SUBROUTINE core_ae_value

! **************************************************************************************************
!> \brief Hermite evaluation of the original cubic signed R=u/r spline, including its adjoints.
!> \param cache Original-data cache.
!> \param shell Occupied shell index.
!> \param radius Evaluation radius, using the first positive interval for the smooth-origin blend.
!> \param value Signed radial wave.
!> \param derivative First radial derivative.
!> \param second Second radial derivative.
! **************************************************************************************************
   SUBROUTINE core_wave_value(cache, shell, radius, value, derivative, second)
      TYPE(paw_core_kinetic_type), INTENT(IN)            :: cache
      INTEGER, INTENT(IN)                                :: shell
      REAL(dp), INTENT(IN)                               :: radius
      REAL(dp), INTENT(OUT)                              :: value, derivative, second

      INTEGER                                            :: hi, lo, mid
      REAL(dp)                                           :: c0, c1, c2, c3, difference, h, right, t

      lo = 1
      hi = SIZE(cache%r)
      DO WHILE (hi - lo > 1)
         mid = (lo + hi)/2
         IF (radius < cache%r(mid)) THEN
            hi = mid
         ELSE
            lo = mid
         END IF
      END DO
      h = cache%r(hi) - cache%r(lo)
      t = (radius - cache%r(lo))/h
      c0 = cache%wave(lo, shell)
      c1 = h*cache%slope(lo, shell)
      right = h*cache%slope(hi, shell)
      difference = cache%wave(hi, shell) - c0
      c2 = 3.0_dp*difference - 2.0_dp*c1 - right
      c3 = -2.0_dp*difference + c1 + right
      value = c0 + t*(c1 + t*(c2 + t*c3))
      derivative = (c1 + t*(2.0_dp*c2 + 3.0_dp*t*c3))/h
      second = (2.0_dp*c2 + 6.0_dp*t*c3)/h**2
   END SUBROUTINE core_wave_value

! **************************************************************************************************
!> \brief Reconstruct AE frozen-core kinetic density from a uniquely verified closed GIPAW core.
!>        This helper does not provide PS core kinetic density or enable a tau-dependent XC method.
!> \param pot Original scalar PAW UPF data; orbital amplitudes are never normalized or modified.
!> \param tau_ae True spin-summed kinetic density, including 1/2, on the full original radial mesh.
!>        Allocated only on success. At r=0 the stored zero has zero cubature weight and is not
!>        a physical limit or a valid interpolation anchor; derivatives use only positive radii.
!> \param available Whether a unique complete closed core passed all independent data checks.
!> \param failure_reason Specific unavailable-data cause, empty on success.
!> \param occupations Optional selected spin-summed occupations in the original GIPAW index order.
! **************************************************************************************************
   SUBROUTINE paw_core_kinetic_build_ae(pot, tau_ae, available, failure_reason, occupations)
      TYPE(atom_upfpot_type), INTENT(IN)                 :: pot
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: tau_ae(:)
      LOGICAL, INTENT(OUT)                               :: available
      CHARACTER(LEN=*), INTENT(OUT)                      :: failure_reason
      INTEGER, ALLOCATABLE, INTENT(OUT), OPTIONAL         :: occupations(:)

      ! Original SSSP tables have maxima 1.45e-5 (norm), 3.45e-5 (charge),
      ! 1.34e-7 (q maximum) and 7.96e-8 (q radial L1). These are independent
      ! original-tabulation consistency checks, not decimal-rounding estimates.
      REAL(dp), PARAMETER :: norm_tolerance = 1.e-4_dp, charge_tolerance = 1.e-4_dp, &
                             maximum_tolerance = 1.e-6_dp, radial_tolerance = 1.e-6_dp

      INTEGER                                            :: first, i, j, l, matches, ncore, norb, nr, &
                                                            nuclear_charge, number_of_candidates, &
                                                            normalized_candidates
      INTEGER, ALLOCATABLE                               :: chosen(:), closed(:), selected(:), &
                                                            remaining(:)
      LOGICAL                                            :: found, maximum_passed, radial_passed
      REAL(dp)                                           :: charge, core_count, peak, radial_scale
      REAL(dp), ALLOCATABLE                              :: derivative(:), norm(:), q_core(:), &
                                                            q_orbital(:, :), radial(:), weights(:)

      available = .FALSE.
      failure_reason = 'AE core tau requires scalar PAW cache'
      IF (.NOT. pot%is_paw .OR. pot%has_so) RETURN
      IF (.NOT. ALLOCATED(pot%paw)) RETURN
      failure_reason = 'Missing AE core density or radial mesh'
      IF (.NOT. ALLOCATED(pot%r) .OR. .NOT. ALLOCATED(pot%rab)) RETURN
      IF (.NOT. ALLOCATED(pot%paw%rho_core)) RETURN
      nr = SIZE(pot%r)
      failure_reason = 'Invalid AE core density or radial mesh'
      IF (nr < 3) RETURN
      IF (SIZE(pot%rab) /= nr .OR. SIZE(pot%paw%rho_core) /= nr) RETURN
      IF (.NOT. ALL(ieee_is_finite(pot%r)) .OR. .NOT. ALL(ieee_is_finite(pot%rab))) RETURN
      IF (.NOT. ALL(ieee_is_finite(pot%paw%rho_core))) RETURN
      IF (ANY(pot%r < 0.0_dp) .OR. ANY(pot%rab <= 0.0_dp) .OR. ANY(pot%paw%rho_core < 0.0_dp)) RETURN
      IF (ANY(pot%r(2:) <= pot%r(:nr - 1))) RETURN
      failure_reason = 'Invalid PAW core electron count'
      CALL get_ptable_info(ADJUSTL(pot%symbol), number=nuclear_charge, found=found)
      IF (.NOT. found) RETURN
      IF (nuclear_charge < 1 .OR. .NOT. ieee_is_finite(pot%zion)) RETURN
      core_count = REAL(nuclear_charge, dp) - pot%zion
      IF (core_count < 0.0_dp .OR. core_count > REAL(nuclear_charge, dp)) RETURN
      ncore = NINT(core_count)
      IF (ABS(core_count - REAL(ncore, dp)) > 64.0_dp*EPSILON(1.0_dp)*REAL(nuclear_charge, dp)) RETURN

      IF (ncore == 0) THEN
         failure_reason = 'Core-free PAW cache contain nonzero AE core density'
         IF (ANY(pot%paw%rho_core /= 0.0_dp)) RETURN
         ALLOCATE (tau_ae(nr), SOURCE=0.0_dp)
         IF (PRESENT(occupations)) THEN
            norb = 0
            IF (ALLOCATED(pot%paw%gipaw_core_orbitals)) THEN
               IF (ALLOCATED(pot%paw%gipaw_core_orbitals%wfc)) norb = SIZE(pot%paw%gipaw_core_orbitals%wfc, 2)
            END IF
            ALLOCATE (occupations(norb), SOURCE=0)
         END IF
         available = .TRUE.
         failure_reason = ''
         RETURN
      END IF

      ALLOCATE (weights, SOURCE=upf_radial_weights(pot))
      ALLOCATE (q_core, SOURCE=fourpi*pot%r**2*pot%paw%rho_core)
      charge = SUM(weights*q_core)
      failure_reason = 'AE core charge does not match Z minus Z_valence'
      IF (.NOT. ieee_is_finite(charge)) RETURN
      IF (ABS(charge - REAL(ncore, dp)) > charge_tolerance) RETURN
      failure_reason = 'Missing GIPAW core orbital table'
      IF (.NOT. ALLOCATED(pot%paw%gipaw_core_orbitals)) RETURN
      ASSOCIATE (orbitals => pot%paw%gipaw_core_orbitals)
         IF (.NOT. ALLOCATED(orbitals%wfc) .OR. .NOT. ALLOCATED(orbitals%indices)) RETURN
         IF (.NOT. ALLOCATED(orbitals%principal_quantum_numbers) .OR. .NOT. ALLOCATED(orbitals%angular_momenta)) RETURN
         norb = SIZE(orbitals%wfc, 2)
         failure_reason = 'Invalid GIPAW core orbital table'
         IF (SIZE(orbitals%wfc, 1) /= nr .OR. SIZE(orbitals%indices) /= norb) RETURN
         IF (SIZE(orbitals%principal_quantum_numbers) /= norb .OR. SIZE(orbitals%angular_momenta) /= norb) RETURN
         IF (.NOT. ALL(ieee_is_finite(orbitals%wfc))) RETURN
         IF (ANY(orbitals%angular_momenta < 0)) RETURN
         IF (ANY(orbitals%principal_quantum_numbers <= orbitals%angular_momenta)) RETURN
         DO i = 1, norb
            IF (orbitals%indices(i) /= i) RETURN
         END DO
         ALLOCATE (closed(norb), selected(norb), chosen(norb), remaining(norb + 1), &
                   norm(norb), q_orbital(nr, norb))
         DO i = 1, norb
            l = orbitals%angular_momenta(i)
            ! A shell larger than this entire core cannot be selected. Capping
            ! its count also avoids integer overflow for invalid large l metadata.
            closed(i) = ncore + 1
            IF (l <= (ncore - 2)/4) closed(i) = 2*(2*l + 1)
            norm(i) = SUM(weights*orbitals%wfc(:, i)**2)
            q_orbital(:, i) = REAL(closed(i), dp)*orbitals%wfc(:, i)**2
         END DO
      END ASSOCIATE
      failure_reason = 'Nonfinite GIPAW core orbital norm'
      IF (.NOT. ALL(ieee_is_finite(norm)) .OR. .NOT. ALL(ieee_is_finite(q_orbital))) RETURN
      remaining(norb + 1) = 0
      DO i = norb, 1, -1
         remaining(i) = MIN(ncore, remaining(i + 1) + MIN(ncore, closed(i)))
      END DO
      failure_reason = 'Incomplete GIPAW core set: insufficient closed-shell electrons'
      IF (remaining(1) < ncore) RETURN
      selected(:) = 0
      chosen(:) = 0
      matches = 0
      number_of_candidates = 0
      normalized_candidates = 0
      maximum_passed = .FALSE.
      radial_passed = .FALSE.
      peak = MAXVAL(q_core)
      radial_scale = SUM(weights*ABS(q_core))
      CALL visit_core_subset(1, 0)
      IF (matches == 0) THEN
         IF (number_of_candidates == 0) THEN
            failure_reason = 'No closed-shell GIPAW subset matches the core electron count'
         ELSE IF (normalized_candidates == 0) THEN
            failure_reason = 'GIPAW core orbital normalization check failed'
         ELSE IF (.NOT. maximum_passed) THEN
            failure_reason = 'GIPAW core radial q maximum check failed'
         ELSE IF (.NOT. radial_passed) THEN
            failure_reason = 'GIPAW core radial q L1 check failed'
         ELSE
            failure_reason = 'No GIPAW core subset passes both radial q checks'
         END IF
         RETURN
      END IF
      failure_reason = 'Ambiguous GIPAW closed-core assignment'
      IF (matches > 1) RETURN

      first = 1
      IF (pot%r(1) == 0.0_dp) first = 2
      failure_reason = 'AE core tau requires at least four positive radial samples'
      IF (nr - first + 1 < 4) RETURN
      ALLOCATE (tau_ae(nr), SOURCE=0.0_dp)
      ALLOCATE (radial(nr - first + 1), derivative(nr - first + 1))
      DO i = 1, norb
         IF (chosen(i) == 0) CYCLE
         radial(:) = pot%paw%gipaw_core_orbitals%wfc(first:, i)/pot%r(first:)
         CALL spline3ders(pot%r(first:), radial, pot%r(first:), dynew=derivative)
         l = pot%paw%gipaw_core_orbitals%angular_momenta(i)
         tau_ae(first:) = tau_ae(first:) + 0.5_dp*REAL(chosen(i), dp)/fourpi* &
                          (derivative**2 + REAL(l, dp)*(REAL(l, dp) + 1.0_dp)*(radial/pot%r(first:))**2)
      END DO
      failure_reason = 'Nonfinite reconstructed AE core tau'
      IF (.NOT. ALL(ieee_is_finite(tau_ae))) THEN
         DEALLOCATE (tau_ae)
         RETURN
      END IF
      IF (PRESENT(occupations)) ALLOCATE (occupations, SOURCE=chosen)
      available = .TRUE.
      failure_reason = ''

   CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param index ...
!> \param electrons ...
! **************************************************************************************************
      RECURSIVE SUBROUTINE visit_core_subset(index, electrons)
      INTEGER, INTENT(IN)                                :: index, electrons

      REAL(dp)                                           :: delta(nr), maximum_error, radial_error

         IF (matches > 1 .OR. electrons > ncore) RETURN
         IF (electrons + remaining(index) < ncore) RETURN
         IF (index <= norb) THEN
            selected(index) = 0
            CALL visit_core_subset(index + 1, electrons)
            selected(index) = closed(index)
            CALL visit_core_subset(index + 1, electrons + closed(index))
            selected(index) = 0
            RETURN
         END IF
         number_of_candidates = number_of_candidates + 1
         IF (ANY(ABS(norm - 1.0_dp) > norm_tolerance .AND. selected > 0)) RETURN
         normalized_candidates = normalized_candidates + 1
         delta = -q_core
         DO j = 1, norb
            IF (selected(j) > 0) delta = delta + q_orbital(:, j)
         END DO
         maximum_error = MAXVAL(ABS(delta))/peak
         radial_error = SUM(weights*ABS(delta))/radial_scale
         maximum_passed = maximum_passed .OR. maximum_error <= maximum_tolerance
         radial_passed = radial_passed .OR. radial_error <= radial_tolerance
         IF (maximum_error > maximum_tolerance .OR. radial_error > radial_tolerance) RETURN
         matches = matches + 1
         chosen(:) = selected
      END SUBROUTINE visit_core_subset

   END SUBROUTINE paw_core_kinetic_build_ae

END MODULE paw_core_kinetic
