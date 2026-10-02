!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM paw_xc_unittest
   USE atom_paw,                        ONLY: atom_paw_add,&
                                              atom_paw_init,&
                                              atom_paw_type
   USE atom_types,                      ONLY: atom_basis_type,&
                                              atom_potential_type,&
                                              upf_pseudo
   USE atom_upf,                        ONLY: atom_read_upf,&
                                              atom_release_upf,&
                                              atom_upfpot_type
   USE atom_utils,                      ONLY: atom_density
   USE atom_xc,                         ONLY: atom_vxc_lda,&
                                              atom_vxc_lsd
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE input_cp2k_xc,                   ONLY: create_xc_fun_section,&
                                              create_xc_section
   USE input_section_types,             ONLY: section_release,&
                                              section_type,&
                                              section_vals_create,&
                                              section_vals_release,&
                                              section_vals_type,&
                                              section_vals_val_set
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE libcp2k,                         ONLY: cp2k_finalize,&
                                              cp2k_init
   USE mathconstants,                   ONLY: fourpi,&
                                              pi
   USE paw_core_kinetic,                ONLY: paw_core_kinetic_cartesian,&
                                              paw_core_kinetic_fourier,&
                                              paw_core_kinetic_init,&
                                              paw_core_kinetic_type,&
                                              paw_core_kinetic_values
   USE paw_onecenter,                   ONLY: paw_onecenter_init,&
                                              paw_onecenter_type
   USE paw_xc,                          ONLY: paw_xc_build_ae_core_tau,&
                                              paw_xc_collocate,&
                                              paw_xc_evaluate,&
                                              paw_xc_init,&
                                              paw_xc_kinetic_adjoint,&
                                              paw_xc_type
   USE upf_projector_integrals,         ONLY: upf_projector_init,&
                                              upf_projector_type
   USE xc_input_constants,              ONLY: xc_pbe_orig,&
                                              xc_pbe_sol

   IMPLICIT NONE
   INTEGER, PARAMETER :: n = 801
   REAL(dp), PARAMETER :: exponent = 0.9_dp, pseudo_scale = 0.8_dp**2 + 0.07_dp
   REAL(dp), PARAMETER :: core_exponents(4) = [1.3_dp, 0.55_dp, 0.8_dp, 0.25_dp]
   TYPE(atom_upfpot_type) :: pot
   TYPE(upf_projector_type) :: projectors
   TYPE(paw_onecenter_type) :: onecenter
   TYPE(paw_xc_type) :: grid
   TYPE(section_type), POINTER :: specification
   TYPE(section_vals_type), POINTER :: functional
   INTEGER :: i, j, l, pair

   CALL cp2k_init()
   pot%is_paw = .TRUE.
   pot%number_of_proj = 2
   ALLOCATE (pot%paw, pot%r(n), pot%rab(n), pot%lbeta(2), pot%beta(n, 2), &
             pot%dion(2, 2), pot%qijl(n, 3, 0:2), pot%rho_nlcc(n))
   ALLOCATE (pot%paw%ae_wfc(n, 2), pot%paw%ps_wfc(n, 2), pot%paw%rho_core(n))
   pot%paw%cutoff_index = n
   pot%lbeta(:) = [0, 1]
   pot%dion = 0.0_dp
   pot%qijl = 0.0_dp
   pot%rab = 8.0_dp/REAL(n - 1, dp)
   DO i = 1, n
      pot%r(i) = REAL(i - 1, dp)*pot%rab(i)
   END DO
   pot%paw%rho_core(:) = 0.2_dp*EXP(-exponent*pot%r**2)
   pot%rho_nlcc(:) = 0.1_dp*EXP(-exponent*pot%r**2)
   DO i = 1, 2
      pot%beta(:, i) = pot%r**pot%lbeta(i)*EXP(-pot%r**2)
      pot%paw%ae_wfc(:, i) = pot%r**(pot%lbeta(i) + 1)*EXP(-0.5_dp*exponent*pot%r**2)
      pot%paw%ps_wfc(:, i) = 0.8_dp*pot%paw%ae_wfc(:, i)
   END DO
   DO j = 1, 2
      DO i = 1, j
         pair = j*(j - 1)/2 + i
         DO l = ABS(pot%lbeta(i) - pot%lbeta(j)), pot%lbeta(i) + pot%lbeta(j), 2
            pot%qijl(:, pair, l) = 0.07_dp*pot%r**(pot%lbeta(i) + pot%lbeta(j) + 2)*EXP(-exponent*pot%r**2)
         END DO
      END DO
   END DO
   CALL upf_projector_init(pot, 1, projectors)
   CALL paw_onecenter_init(pot, projectors, onecenter)
   CALL paw_xc_init(pot, 50, grid)
   NULLIFY (specification, functional)
   CALL create_xc_fun_section(specification)
   CALL section_vals_create(functional, specification)
   CALL check_cartesian_fields()
   CALL check_kinetic_adjoint()
   CALL check_local_exchange()
   CALL check_gga()
   CALL check_meta_gga()
   CALL check_ae_core_tau()
   CALL check_shared_core_tau()
   CALL check_atomic_projection()
   CALL check_atomic_gradient_norm()
   CALL original_datasets()
   CALL section_vals_release(functional)
   CALL section_release(specification)
   CALL cp2k_finalize()

CONTAINS

! **************************************************************************************************
!> \brief Independent Gaussian core tau, complete-core assignment and unavailable-data guards.
! **************************************************************************************************
   SUBROUTINE check_ae_core_tau()
      CHARACTER(LEN=256)                                 :: reason
      INTEGER                                            :: first, i, origin, peak_index
      INTEGER, ALLOCATABLE                               :: occupations(:)
      LOGICAL                                            :: available
      REAL(dp)                                           :: amplitude, error, peak, positive_weight
      REAL(dp), ALLOCATABLE                              :: delta(:), original_u(:, :), q_core(:), &
                                                            reference(:), tau(:), weights(:)
      TYPE(atom_upfpot_type)                             :: baseline, changed, s_core

      error = 0.0_dp
      DO origin = 0, 1
         CALL make_core_tau_fixture(baseline, origin == 1)
         ALLOCATE (original_u, SOURCE=baseline%paw%gipaw_core_orbitals%wfc)
         CALL paw_xc_build_ae_core_tau(baseline, tau, available, reason, occupations)
         IF (.NOT. available) THEN
            PRINT *, 'Complete Gaussian PAW core unavailable:', TRIM(reason)
            ERROR STOP 'Complete Gaussian PAW core unavailable'
         END IF
         IF (LEN_TRIM(reason) /= 0) ERROR STOP 'Available PAW core retained an unavailable-data reason'
         IF (.NOT. ALLOCATED(tau)) ERROR STOP 'Available PAW core tau not allocated'
         IF (.NOT. ALLOCATED(occupations)) ERROR STOP 'Available PAW core occupations not allocated'
         IF (SIZE(tau) /= n) ERROR STOP 'PAW core tau truncated at the augmentation cutoff'
         IF (ANY(occupations /= [2, 2, 6, 0])) ERROR STOP 'GIPAW valence orbital was assigned to the frozen core'
         IF (ANY(baseline%paw%gipaw_core_orbitals%wfc /= original_u)) THEN
            ERROR STOP 'AE core tau changed original signed GIPAW orbital amplitudes'
         END IF
         IF (.NOT. ALL(ieee_is_finite(tau)) .OR. ANY(tau < 0.0_dp)) ERROR STOP 'Invalid Gaussian AE core tau'
         ALLOCATE (reference(n))
         reference(:) = 0.0_dp
         ! For normalized Gaussian radial orbitals the independent analytic
         ! s-shell factor is a**2*r**2; the p shell is 3-2*a*r**2+a**2*r**4.
         DO i = 1, 2
            amplitude = 2.0_dp*core_exponents(i)**1.5_dp/GAMMA(1.5_dp)
            reference(:) = reference + amplitude*core_exponents(i)**2*baseline%r**2* &
                           EXP(-core_exponents(i)*baseline%r**2)/fourpi
         END DO
         amplitude = 2.0_dp*core_exponents(3)**2.5_dp/GAMMA(2.5_dp)
         reference(:) = reference + 3.0_dp*amplitude*EXP(-core_exponents(3)*baseline%r**2)* &
                        (3.0_dp - 2.0_dp*core_exponents(3)*baseline%r**2 + &
                         core_exponents(3)**2*baseline%r**4)/fourpi
         first = 1
         IF (origin == 1) THEN
            first = 2
            ! A nonzero analytic p-shell origin is deliberately not replaced by
            ! the storage convention. This node is excluded from all derivatives.
            IF (tau(1) /= 0.0_dp .OR. reference(1) <= 0.0_dp) ERROR STOP 'PAW core tau origin convention changed'
         END IF
         error = MAX(error, MAXVAL(ABS(tau(first:) - reference(first:)))/MAXVAL(reference(first:)))
         IF (MAXVAL(tau(baseline%paw%cutoff_index + 1:)) <= 0.0_dp) ERROR STOP 'PAW core tau tail lost'
         DEALLOCATE (reference, original_u)
         CALL atom_release_upf(baseline)
      END DO
      WRITE (*, '(A,ES15.6)') 'PAW Gaussian AE core kinetic-density relative error:', error
      IF (error > 2.e-6_dp) ERROR STOP 'PAW AE core tau differs from independent Gaussian formula'

      CALL make_core_tau_fixture(baseline, .FALSE.)
      changed = baseline
      DEALLOCATE (changed%paw%gipaw_core_orbitals)
      CALL require_core_tau_unavailable(changed, 'Missing')

      changed = baseline
      ! Removing the p shell leaves only six supplied closed-shell electrons,
      ! while the independent AE core density and element still require ten.
      DEALLOCATE (changed%paw%gipaw_core_orbitals%wfc, changed%paw%gipaw_core_orbitals%indices, &
                  changed%paw%gipaw_core_orbitals%principal_quantum_numbers, &
                  changed%paw%gipaw_core_orbitals%angular_momenta)
      ALLOCATE (changed%paw%gipaw_core_orbitals%wfc(n, 3), SOURCE=baseline%paw%gipaw_core_orbitals%wfc(:, [1, 2, 4]))
      ALLOCATE (changed%paw%gipaw_core_orbitals%indices(3), SOURCE=[1, 2, 3])
      ALLOCATE (changed%paw%gipaw_core_orbitals%principal_quantum_numbers(3), SOURCE=[1, 2, 3])
      ALLOCATE (changed%paw%gipaw_core_orbitals%angular_momenta(3), SOURCE=[0, 0, 0])
      CALL require_core_tau_unavailable(changed, 'insufficient')

      changed = baseline
      changed%symbol = 'Be'
      changed%zion = 2.0_dp
      changed%paw%rho_core(:) = 2.0_dp*baseline%paw%gipaw_core_orbitals%wfc(:, 1)**2/(fourpi*baseline%r**2)
      changed%paw%gipaw_core_orbitals%wfc(:, 2) = -changed%paw%gipaw_core_orbitals%wfc(:, 1)
      CALL require_core_tau_unavailable(changed, 'Ambiguous')

      s_core = baseline
      s_core%symbol = 'Be'
      s_core%zion = 2.0_dp
      DEALLOCATE (s_core%paw%gipaw_core_orbitals%wfc, s_core%paw%gipaw_core_orbitals%indices, &
                  s_core%paw%gipaw_core_orbitals%principal_quantum_numbers, &
                  s_core%paw%gipaw_core_orbitals%angular_momenta)
      ALLOCATE (s_core%paw%gipaw_core_orbitals%wfc(n, 1), SOURCE=baseline%paw%gipaw_core_orbitals%wfc(:, 1:1))
      ALLOCATE (s_core%paw%gipaw_core_orbitals%indices(1), SOURCE=[1])
      ALLOCATE (s_core%paw%gipaw_core_orbitals%principal_quantum_numbers(1), SOURCE=[1])
      ALLOCATE (s_core%paw%gipaw_core_orbitals%angular_momenta(1), SOURCE=[0])
      s_core%paw%rho_core(:) = 2.0_dp*s_core%paw%gipaw_core_orbitals%wfc(:, 1)**2/(fourpi*s_core%r**2)
      changed = s_core
      changed%paw%gipaw_core_orbitals%wfc(:, :) = 1.01_dp*changed%paw%gipaw_core_orbitals%wfc
      CALL require_core_tau_unavailable(changed, 'normalization')
      changed = s_core
      changed%paw%rho_core(:) = 1.001_dp*changed%paw%rho_core
      CALL require_core_tau_unavailable(changed, 'charge')

      ALLOCATE (q_core, SOURCE=fourpi*s_core%r**2*s_core%paw%rho_core)
      peak = MAXVAL(q_core)
      peak_index = MAXLOC(q_core, DIM=1)
      changed = s_core
      ! A one-node perturbation violates the maximum check while its radial
      ! integral is small enough to satisfy both charge and L1 requirements.
      changed%paw%rho_core(peak_index) = changed%paw%rho_core(peak_index) + &
                                         3.e-6_dp*peak/(fourpi*changed%r(peak_index)**2)
      CALL require_core_tau_unavailable(changed, 'maximum')

      ! Independent Simpson weights on this known uniform, odd-sized fixture.
      ALLOCATE (weights(n), delta(n))
      weights(:) = (2.0_dp/3.0_dp)*s_core%rab
      weights(2:n - 1:2) = (4.0_dp/3.0_dp)*s_core%rab(2:n - 1:2)
      weights(1) = s_core%rab(1)/3.0_dp
      weights(n) = s_core%rab(n)/3.0_dp
      delta(:) = 0.0_dp
      WHERE (s_core%r >= 0.5_dp .AND. s_core%r <= 2.5_dp) delta(:) = -0.9e-6_dp*peak
      positive_weight = -SUM(weights*delta)/SUM(weights, MASK=s_core%r > 2.5_dp)
      WHERE (s_core%r > 2.5_dp) delta(:) = positive_weight
      changed = s_core
      changed%paw%rho_core(:) = changed%paw%rho_core + delta/(fourpi*changed%r**2)
      IF (ANY(changed%paw%rho_core < 0.0_dp)) ERROR STOP 'Invalid positive core-density L1 fixture'
      CALL require_core_tau_unavailable(changed, 'L1')

      changed = baseline
      changed%symbol = 'He'
      changed%zion = 2.0_dp
      changed%paw%rho_core(:) = 0.0_dp
      CALL paw_xc_build_ae_core_tau(changed, tau, available, reason, occupations)
      IF (.NOT. available) ERROR STOP 'Core-free He unavailable with GIPAW valence orbitals'
      IF (.NOT. ALLOCATED(tau) .OR. .NOT. ALLOCATED(occupations)) ERROR STOP 'Missing core-free He outputs'
      IF (ANY(tau /= 0.0_dp) .OR. ANY(occupations /= 0)) THEN
         ERROR STOP 'Core-free He acquired core tau from GIPAW valence orbitals'
      END IF
      DEALLOCATE (changed%paw%gipaw_core_orbitals)
      CALL paw_xc_build_ae_core_tau(changed, tau, available, reason)
      IF (.NOT. available) ERROR STOP 'Core-free PAW tau requires GIPAW orbitals'
      IF (.NOT. ALLOCATED(tau)) ERROR STOP 'Missing core-free PAW tau'
      IF (ANY(tau /= 0.0_dp)) ERROR STOP 'Core-free PAW tau is not zero on the full mesh'
      CALL atom_release_upf(baseline)
      CALL atom_release_upf(changed)
      CALL atom_release_upf(s_core)
   END SUBROUTINE check_ae_core_tau

! **************************************************************************************************
!> \brief Independent Gaussian smooth-core values and derivatives, including complex Fourier adjoints.
! **************************************************************************************************
   SUBROUTINE check_shared_core_tau()
      REAL(dp), PARAMETER                                :: step = 2.e-6_dp

      COMPLEX(dp)                                        :: center(3), cminus, coefficient, cplus, &
                                                            strain(3, 3)
      INTEGER                                            :: a, b, i, mode, origin
      REAL(dp) :: analytic, analytic_derivative, anchor, cartesian_error, derivative_error, error, &
         g(3), gm(3), gp(3), gradient(3), minus, plus, pm(3), point(3), position(3), pp(3), &
         radius(12), rc, reciprocal_error, shifted(3), slope(12), tau(12), unused, value, vm, vp
      TYPE(atom_upfpot_type)                             :: original
      TYPE(paw_core_kinetic_type)                        :: data

      error = 0.0_dp
      derivative_error = 0.0_dp
      cartesian_error = 0.0_dp
      reciprocal_error = 0.0_dp
      DO origin = 0, 1
         CALL make_core_tau_fixture(original, origin == 1)
         CALL paw_core_kinetic_init(original, DATA)
         IF (.NOT. data%available) ERROR STOP 'Verified Gaussian core has no shared kinetic fields'
         IF (SUM(data%occupations) /= 10) ERROR STOP 'Shared kinetic fields changed the closed-core selection'
         rc = data%smoothing_radius
         radius = [0.0_dp, 1.e-14_dp, 0.0003_dp, 0.237_dp, 0.914_dp, 1.783_dp, rc, rc + 0.031_dp, &
                   data%taper_radius, data%taper_radius + 0.007_dp, data%support_radius, data%support_radius + 0.1_dp]
         CALL paw_core_kinetic_values(DATA, radius, tau, slope)
         IF (.NOT. ALL(ieee_is_finite(tau)) .OR. ANY(tau < 0.0_dp)) ERROR STOP 'Invalid smooth core tau'
         CALL gaussian_core_kinetic(rc, anchor, unused)
         DO i = 1, SIZE(radius)
            CALL gaussian_core_kinetic(radius(i), analytic, analytic_derivative)
            CALL gaussian_core_blend(radius(i), rc, data%taper_radius, data%support_radius, anchor, &
                                     analytic, analytic_derivative)
            error = MAX(error, ABS(tau(i) - analytic))
            derivative_error = MAX(derivative_error, ABS(slope(i) - analytic_derivative))
         END DO
         IF (slope(1) /= 0.0_dp .OR. tau(11) /= 0.0_dp .OR. slope(11) /= 0.0_dp) THEN
            ERROR STOP 'Smooth core origin or compact support closure is incorrect'
         END IF
         point = [0.713_dp, -0.491_dp, 0.317_dp]
         CALL paw_core_kinetic_cartesian(DATA, point, value, gradient)
         DO a = 1, 3
            shifted = point
            shifted(a) = shifted(a) + step
            CALL paw_core_kinetic_cartesian(DATA, shifted, plus)
            shifted(a) = shifted(a) - 2.0_dp*step
            CALL paw_core_kinetic_cartesian(DATA, shifted, minus)
            cartesian_error = MAX(cartesian_error, ABS((plus - minus)/(2.0_dp*step) - gradient(a)))
         END DO
         position = [0.31_dp, -0.73_dp, 0.19_dp]
         DO mode = 0, 1
            g = REAL(mode, dp)*[0.71_dp, -0.37_dp, 1.13_dp]
            CALL paw_core_kinetic_fourier(DATA, g, position, 17.0_dp, coefficient, strain, center)
            IF (mode == 0) THEN
               IF (ABS(coefficient - SUM(data%weighted)/17.0_dp) > 1.e-14_dp) THEN
                  ERROR STOP 'Zero-G core kinetic normalization is incorrect'
               END IF
            ELSE
               IF (ABS(AIMAG(coefficient)) < 1.e-8_dp) ERROR STOP 'Core transform lost its complex center phase'
            END IF
            DO a = 1, 3
               pp = position
               pp(a) = pp(a) + step
               CALL paw_core_kinetic_fourier(DATA, g, pp, 17.0_dp, cplus)
               pp(a) = pp(a) - 2.0_dp*step
               CALL paw_core_kinetic_fourier(DATA, g, pp, 17.0_dp, cminus)
               reciprocal_error = MAX(reciprocal_error, ABS((cplus - cminus)/(2.0_dp*step) - center(a)))
               DO b = 1, 3
                  gp = g
                  gm = g
                  pp = position
                  pm = position
                  pp(a) = pp(a) + step*position(b)
                  pm(a) = pm(a) - step*position(b)
                  vp = 17.0_dp
                  vm = 17.0_dp
                  IF (a == b) THEN
                     gp(a) = g(a)/(1.0_dp + step)
                     gm(a) = g(a)/(1.0_dp - step)
                     vp = vp*(1.0_dp + step)
                     vm = vm*(1.0_dp - step)
                  ELSE
                     gp(b) = gp(b) - step*g(a)
                     gm(b) = gm(b) + step*g(a)
                  END IF
                  CALL paw_core_kinetic_fourier(DATA, gp, pp, vp, cplus)
                  CALL paw_core_kinetic_fourier(DATA, gm, pm, vm, cminus)
                  reciprocal_error = MAX(reciprocal_error, ABS((cplus - cminus)/(2.0_dp*step) - strain(a, b)))
               END DO
            END DO
         END DO
         CALL atom_release_upf(original)
      END DO
      WRITE (*, '(A,4ES15.6)') 'PAW shared core tau value/radial/Cartesian/Fourier errors:', &
         error, derivative_error, cartesian_error, reciprocal_error
      IF (error > 2.e-6_dp) ERROR STOP 'Smooth core tau differs from independent Gaussian values'
      IF (derivative_error > 2.e-5_dp) ERROR STOP 'Smooth core slope differs from independent Gaussian derivatives'
      IF (cartesian_error > 1.e-8_dp) ERROR STOP 'Core spatial gradient is not the field derivative'
      IF (reciprocal_error > 1.e-8_dp) ERROR STOP 'Complex core Fourier adjoints are not transform derivatives'
   END SUBROUTINE check_shared_core_tau

! **************************************************************************************************
!> \brief Closed-form independent s/p frozen-core kinetic density and radial derivative.
!> \param radius Evaluation radius.
!> \param value Spin-summed kinetic density.
!> \param derivative Radial derivative.
! **************************************************************************************************
   SUBROUTINE gaussian_core_kinetic(radius, value, derivative)
      REAL(dp), INTENT(IN)                               :: radius
      REAL(dp), INTENT(OUT)                              :: value, derivative

      INTEGER                                            :: i
      REAL(dp)                                           :: a, amplitude, factor, polynomial

      value = 0.0_dp
      derivative = 0.0_dp
      DO i = 1, 2
         a = core_exponents(i)
         amplitude = 2.0_dp*a**1.5_dp/GAMMA(1.5_dp)
         factor = amplitude*EXP(-a*radius**2)/fourpi
         value = value + factor*a**2*radius**2
         derivative = derivative + factor*a**2*(2.0_dp*radius - 2.0_dp*a*radius**3)
      END DO
      a = core_exponents(3)
      amplitude = 2.0_dp*a**2.5_dp/GAMMA(2.5_dp)
      factor = 3.0_dp*amplitude*EXP(-a*radius**2)/fourpi
      polynomial = 3.0_dp - 2.0_dp*a*radius**2 + a**2*radius**4
      value = value + factor*polynomial
      derivative = derivative + factor*(-4.0_dp*a*radius + 4.0_dp*a**2*radius**3 - 2.0_dp*a*radius*polynomial)
   END SUBROUTINE gaussian_core_kinetic

! **************************************************************************************************
!> \brief Independent analytic auxiliary profile used with the closed-form Gaussian core oracle.
!> \param radius Evaluation radius.
!> \param rc Augmentation radius.
!> \param taper Start of the compact outer taper.
!> \param support End of the radial support.
!> \param anchor AE core tau at rc.
!> \param value AE input, overwritten by the auxiliary PS value.
!> \param derivative AE derivative input, overwritten by the auxiliary PS derivative.
! **************************************************************************************************
   SUBROUTINE gaussian_core_blend(radius, rc, taper, support, anchor, value, derivative)
      REAL(dp), INTENT(IN)                               :: radius, rc, taper, support, anchor
      REAL(dp), INTENT(INOUT)                            :: value, derivative

      REAL(dp)                                           :: blend, dw, gaussian, slope, t, u, window

      IF (radius >= support) THEN
         value = 0.0_dp
         derivative = 0.0_dp
         RETURN
      END IF
      IF (radius < rc) THEN
         u = (radius/rc)**2
         blend = 10.0_dp*u**3 - 15.0_dp*u**4 + 6.0_dp*u**5
         slope = (30.0_dp*u**2 - 60.0_dp*u**3 + 30.0_dp*u**4)*2.0_dp*radius/rc**2
         gaussian = anchor*EXP(1.0_dp - u)
         derivative = blend*derivative + slope*(value - gaussian) - &
                      (1.0_dp - blend)*2.0_dp*radius/rc**2*gaussian
         value = blend*value + (1.0_dp - blend)*gaussian
      END IF
      IF (radius > taper) THEN
         t = (radius - taper)/(support - taper)
         window = 1.0_dp - 10.0_dp*t**3 + 15.0_dp*t**4 - 6.0_dp*t**5
         dw = (-30.0_dp*t**2 + 60.0_dp*t**3 - 30.0_dp*t**4)/(support - taper)
         derivative = window*derivative + dw*value
         value = window*value
      END IF
   END SUBROUTINE gaussian_core_blend

! **************************************************************************************************
!> \brief Build normalized Gaussian radial s/p cores plus an unassigned signed valence orbital.
!> \param potential Synthetic PAW data, without projector or PS kinetic-density assumptions.
!> \param zero_origin Whether the unused zero-radius node is present.
! **************************************************************************************************
   SUBROUTINE make_core_tau_fixture(potential, zero_origin)
      TYPE(atom_upfpot_type), INTENT(OUT)                :: potential
      LOGICAL, INTENT(IN)                                :: zero_origin

      INTEGER                                            :: i, ir, l
      REAL(dp)                                           :: amplitude, first_radius, occupation

      potential%symbol = 'Ar'
      potential%zion = 8.0_dp
      potential%is_paw = .TRUE.
      potential%mesh_size = n
      ALLOCATE (potential%paw, potential%r(n), potential%rab(n))
      ALLOCATE (potential%paw%rho_core(n), potential%paw%gipaw_core_orbitals)
      potential%paw%cutoff_index = n/2
      ASSOCIATE (orbitals => potential%paw%gipaw_core_orbitals)
         orbitals%data_format = 2
         ALLOCATE (orbitals%wfc(n, 4), orbitals%indices(4), orbitals%principal_quantum_numbers(4), &
                   orbitals%angular_momenta(4))
         orbitals%indices(:) = [1, 2, 3, 4]
         orbitals%principal_quantum_numbers(:) = [1, 2, 2, 3]
         orbitals%angular_momenta(:) = [0, 0, 1, 0]
         first_radius = MERGE(0.0_dp, 1.e-4_dp, zero_origin)
         potential%rab(:) = 8.0_dp/REAL(n - 1, dp)
         DO ir = 1, n
            potential%r(ir) = first_radius + REAL(ir - 1, dp)*potential%rab(ir)
         END DO
         potential%paw%rho_core(:) = 0.0_dp
         DO i = 1, 4
            l = orbitals%angular_momenta(i)
            amplitude = SQRT(2.0_dp*core_exponents(i)**(REAL(l, dp) + 1.5_dp)/GAMMA(REAL(l, dp) + 1.5_dp))
            orbitals%wfc(:, i) = amplitude*potential%r**(l + 1)*EXP(-0.5_dp*core_exponents(i)*potential%r**2)
            IF (i == 2) orbitals%wfc(:, i) = -orbitals%wfc(:, i)
            IF (i == 4) CYCLE
            occupation = REAL(2*(2*l + 1), dp)
            potential%paw%rho_core(:) = potential%paw%rho_core + occupation*amplitude**2* &
                                        potential%r**(2*l)*EXP(-core_exponents(i)*potential%r**2)/fourpi
         END DO
      END ASSOCIATE
   END SUBROUTINE make_core_tau_fixture

! **************************************************************************************************
!> \brief Require a specific unavailable-data diagnostic and no usable tau output.
!> \param potential Invalid or incomplete synthetic core data.
!> \param cause Expected diagnostic category.
! **************************************************************************************************
   SUBROUTINE require_core_tau_unavailable(potential, cause)
      TYPE(atom_upfpot_type), INTENT(IN)                 :: potential
      CHARACTER(LEN=*), INTENT(IN)                       :: cause

      CHARACTER(LEN=256)                                 :: reason
      LOGICAL                                            :: available
      REAL(dp), ALLOCATABLE                              :: tau(:)

      CALL paw_xc_build_ae_core_tau(potential, tau, available, reason)
      IF (available .OR. ALLOCATED(tau) .OR. INDEX(reason, cause) == 0) THEN
         PRINT *, 'PAW core tau availability/cause:', available, TRIM(reason), TRIM(cause)
         ERROR STOP 'Invalid PAW core data accepted or unavailable-data cause lost'
      END IF
   END SUBROUTINE require_core_tau_unavailable

! **************************************************************************************************
!> \brief TPSS energy derivatives, spin symmetry, and explicit frozen-core kinetic fields.
! **************************************************************************************************
   SUBROUTINE check_meta_gga()
      REAL(dp), PARAMETER                                :: step = 1.e-5_dp
#if defined(__LIBXC)
      INTEGER, PARAMETER                                 :: nvariant = 2
#else
      INTEGER, PARAMETER                                 :: nvariant = 1
#endif

      INTEGER                                            :: core_case, maxspin, ns, p, q, side, variant
      REAL(dp) :: check_energy, core_kinetic(n, 2), density(4, 4, 2), derivative, &
                  direction(4, 4, 2), energy, errors(2), minus, plus, potential(4, 4, 2), &
                  restricted(4, 4, 1), restricted_potential(4, 4, 1), scratch(4, 4, 2), spin_error
      TYPE(paw_xc_type)                                  :: xcgrid

      DO q = 1, 4
         DO p = 1, 4
            restricted(p, q, 1) = SIN(0.31_dp*p)*SIN(0.31_dp*q)
            IF (p == q) restricted(p, q, 1) = restricted(p, q, 1) + 0.3_dp
            direction(p, q, 1) = COS(0.23_dp*(p + q))
            direction(p, q, 2) = SIN(0.37_dp*(p + q))
         END DO
      END DO
      DO variant = 1, nvariant
         maxspin = variant
         CALL section_vals_release(functional)
         CALL section_vals_create(functional, specification)
         IF (variant == 1) THEN
            CALL section_vals_val_set(functional, 'TPSS%_SECTION_PARAMETERS_', l_val=.TRUE.)
         END IF
#if defined(__LIBXC)
         IF (variant == 2) THEN
            CALL section_vals_val_set(functional, 'MGGA_X_TPSS%_SECTION_PARAMETERS_', l_val=.TRUE.)
            CALL section_vals_val_set(functional, 'MGGA_C_TPSS%_SECTION_PARAMETERS_', l_val=.TRUE.)
         END IF
#endif
         DO core_case = 1, 2
            xcgrid = grid
            IF (core_case == 1) THEN
               xcgrid%core = 0.0_dp
               xcgrid%core_gradient = 0.0_dp
            END IF
            ! The synthetic spherical core is a Gaussian s orbital on each side.
            DO side = 1, 2
               core_kinetic(:, side) = 0.5_dp*exponent**2*grid%r**2*xcgrid%core(:, side)
            END DO
            CALL paw_xc_evaluate(xcgrid, onecenter, functional, restricted, 1.e-20_dp, 1.e-20_dp, &
                                 energy, restricted_potential, core_kinetic)
            IF (ALLOCATED(xcgrid%core_kinetic%kinetic)) DEALLOCATE (xcgrid%core_kinetic%kinetic)
            ALLOCATE (xcgrid%core_kinetic%kinetic, SOURCE=core_kinetic)
            xcgrid%core_kinetic%available = .TRUE.
            CALL paw_xc_evaluate(xcgrid, onecenter, functional, restricted, 1.e-20_dp, 1.e-20_dp, &
                                 check_energy, scratch(:, :, :1))
            IF (ABS(check_energy - energy) > 1.e-13_dp .OR. &
                MAXVAL(ABS(scratch(:, :, :1) - restricted_potential)) > 1.e-13_dp) THEN
               ERROR STOP 'Shared and explicit PAW core kinetic fields disagree'
            END IF
            ! An explicit caller field has precedence over the cached auxiliary profile.
            xcgrid%core_kinetic%kinetic(:, :) = 2.0_dp*core_kinetic + 0.1_dp
            CALL paw_xc_evaluate(xcgrid, onecenter, functional, restricted, 1.e-20_dp, 1.e-20_dp, &
                                 check_energy, scratch(:, :, :1), core_kinetic)
            IF (ABS(check_energy - energy) > 1.e-13_dp) ERROR STOP 'Explicit PAW core tau did not override the cache'
            xcgrid%core_kinetic%kinetic(:, :) = core_kinetic
            IF (core_case == 1) THEN
               CALL paw_xc_evaluate(xcgrid, onecenter, functional, restricted, 1.e-20_dp, 1.e-20_dp, &
                                    check_energy, scratch(:, :, :1))
               IF (ABS(check_energy - energy) > 1.e-13_dp) ERROR STOP 'Core-free PAW TPSS depends on optional zero core'
            END IF
            spin_error = 0.0_dp
            IF (maxspin == 2) THEN
               density(:, :, 1) = 0.5_dp*restricted(:, :, 1)
               density(:, :, 2) = density(:, :, 1)
               CALL paw_xc_evaluate(xcgrid, onecenter, functional, density, 1.e-20_dp, 1.e-20_dp, &
                                    check_energy, potential, core_kinetic)
               IF (ABS(check_energy - energy) > 2.e-10_dp) ERROR STOP 'Restricted and unpolarized spin PAW TPSS differ'
               CALL paw_xc_evaluate(xcgrid, onecenter, functional, density, 1.e-20_dp, 1.e-20_dp, &
                                    check_energy, scratch)
               IF (ABS(check_energy - energy) > 2.e-10_dp .OR. MAXVAL(ABS(scratch - potential)) > 1.e-13_dp) THEN
                  ERROR STOP 'Spin-split shared PAW core tau differs from the explicit field'
               END IF

               spin_error = MAXVAL(ABS(potential(:, :, 1) - restricted_potential(:, :, 1)))
            END IF
            errors = 0.0_dp
            DO ns = 1, maxspin
               density(:, :, 1) = restricted(:, :, 1)
               IF (ns == 2) THEN
                  density(:, :, 1) = 0.7_dp*restricted(:, :, 1)
                  density(:, :, 2) = 0.3_dp*restricted(:, :, 1)
               END IF
               CALL paw_xc_evaluate(xcgrid, onecenter, functional, density(:, :, :ns), 1.e-20_dp, 1.e-20_dp, &
                                    energy, potential(:, :, :ns), core_kinetic)
               CALL paw_xc_evaluate(xcgrid, onecenter, functional, density(:, :, :ns) + step*direction(:, :, :ns), &
                                    1.e-20_dp, 1.e-20_dp, plus, scratch(:, :, :ns), core_kinetic)
               CALL paw_xc_evaluate(xcgrid, onecenter, functional, density(:, :, :ns) - step*direction(:, :, :ns), &
                                    1.e-20_dp, 1.e-20_dp, minus, scratch(:, :, :ns), core_kinetic)
               derivative = SUM(potential(:, :, :ns)*direction(:, :, :ns))
               errors(ns) = ABS((plus - minus)/(2.0_dp*step) - derivative)
            END DO
            WRITE (*, '(A,2I2,2ES15.6)') 'PAW TPSS projector-density derivative errors:', core_case, maxspin, errors(:maxspin)
            WRITE (*, '(A,ES15.6)') 'PAW TPSS unpolarized spin derivative error:', spin_error
            IF (MAXVAL(errors) > 3.e-8_dp) ERROR STOP 'PAW TPSS XC matrix is not the energy derivative'
            IF (spin_error > 2.e-9_dp) ERROR STOP 'Restricted and unpolarized spin PAW TPSS derivatives differ'
            IF (maxspin == 2) THEN
               CALL paw_xc_evaluate(xcgrid, onecenter, functional, density(:, :, 2:1:-1), 1.e-20_dp, 1.e-20_dp, &
                                    check_energy, scratch, core_kinetic)
               IF (ABS(check_energy - energy) > 2.e-12_dp) ERROR STOP 'PAW TPSS spin exchange changed energy'
            IF (MAXVAL(ABS(scratch - potential(:, :, 2:1:-1))) > 2.e-11_dp) ERROR STOP 'PAW TPSS spin derivative permutation failed'
            END IF
         END DO
      END DO
      CALL section_vals_release(functional)
      CALL section_vals_create(functional, specification)
      CALL section_vals_val_set(functional, 'PBE%_SECTION_PARAMETERS_', l_val=.TRUE.)
   END SUBROUTINE check_meta_gga

! **************************************************************************************************
!> \brief Native TPSS and its first derivatives at stationary, positive-density points.
! **************************************************************************************************
   SUBROUTINE check_tpss_zero_gradient()
      REAL(dp), PARAMETER :: step = 1.e-5_dp, tau_scales(3) = [0.2_dp, 1.0_dp, 2.0_dp]

      INTEGER                                            :: axis, configuration
      REAL(dp)                                           :: energy, error, gradient(3), minus, plus, &
                                                            scratch(3), tau_unif, trial(3), x(3)

      error = 0.0_dp
      x(1) = 0.4_dp
      tau_unif = 0.3_dp*(3.0_dp*pi**2)**(2.0_dp/3.0_dp)*x(1)**(5.0_dp/3.0_dp)
      DO configuration = 1, SIZE(tau_scales)
         x(2) = 0.0_dp
         x(3) = tau_scales(configuration)*tau_unif
         CALL evaluate_tpss_point(x, energy, gradient)
         IF (ABS(gradient(2)) > 1.e-13_dp) ERROR STOP 'TPSS has a nonzero derivative at zero density gradient'
         DO axis = 1, 3, 2
            trial = x
            trial(axis) = x(axis) + step
            CALL evaluate_tpss_point(trial, plus, scratch)
            trial(axis) = x(axis) - step
            CALL evaluate_tpss_point(trial, minus, scratch)
            error = MAX(error, ABS((plus - minus)/(2.0_dp*step) - gradient(axis)))
         END DO
         trial = x
         trial(2) = 1.e-100_dp
         CALL evaluate_tpss_point(trial, plus, scratch)
         IF (ABS(plus - energy) > 1.e-13_dp .OR. ABS(scratch(2)) > 1.e-20_dp) THEN
            ERROR STOP 'TPSS does not approach its zero-gradient limit continuously'
         END IF
         x(2) = 1.e-4_dp
         CALL evaluate_tpss_point(x, energy, gradient)
         trial = x
         trial(2) = x(2) + step
         CALL evaluate_tpss_point(trial, plus, scratch)
         trial(2) = x(2) - step
         CALL evaluate_tpss_point(trial, minus, scratch)
         error = MAX(error, ABS((plus - minus)/(2.0_dp*step) - gradient(2)))
      END DO
      WRITE (*, '(A,ES15.6)') 'TPSS zero-gradient first-derivative error:', error
      IF (error > 3.e-8_dp) ERROR STOP 'TPSS zero-gradient potential is not the energy derivative'
   END SUBROUTINE check_tpss_zero_gradient

! **************************************************************************************************
!> \brief Evaluate native TPSS in independent density, gradient-norm, and kinetic-density variables.
!> \param x Density, gradient norm, and kinetic density, respectively.
!> \param energy Energy density.
!> \param gradient First derivatives with respect to x.
! **************************************************************************************************
   SUBROUTINE evaluate_tpss_point(x, energy, gradient)
      REAL(dp), INTENT(IN)                               :: x(3)
      REAL(dp), INTENT(OUT)                              :: energy, gradient(3)

      INTEGER, PARAMETER :: labels(3) = [deriv_rho, deriv_norm_drho, deriv_tau]

      INTEGER                                            :: bounds(2, 3), i
      REAL(dp), POINTER                                  :: values(:, :, :)
      TYPE(xc_derivative_set_type)                       :: derivatives
      TYPE(xc_derivative_type), POINTER                  :: derivative
      TYPE(xc_rho_cflags_type)                           :: needs
      TYPE(xc_rho_set_type)                              :: rho_set

      bounds = 1
      needs = xc_functionals_get_needs(functional, lsd=.FALSE., calc_potential=.TRUE.)
      CALL xc_rho_set_create(rho_set, bounds, rho_cutoff=1.e-20_dp, tau_cutoff=1.e-20_dp)
      CALL xc_rho_set_atom_update(rho_set, needs, 1, bounds)
      rho_set%rho = x(1)
      rho_set%norm_drho = x(2)
      rho_set%tau = x(3)
      rho_set%has%rho = .TRUE.
      rho_set%has%norm_drho = .TRUE.
      rho_set%has%tau = .TRUE.
      CALL xc_dset_create(derivatives, local_bounds=bounds)
      CALL xc_functionals_eval(functional, lsd=.FALSE., rho_set=rho_set, deriv_set=derivatives, deriv_order=1)
      derivative => xc_dset_get_derivative(derivatives, [INTEGER::])
      CALL xc_derivative_get(derivative, deriv_data=values)
      energy = values(1, 1, 1)
      DO i = 1, SIZE(labels)
         derivative => xc_dset_get_derivative(derivatives, [labels(i)])
         CALL xc_derivative_get(derivative, deriv_data=values)
         gradient(i) = values(1, 1, 1)
      END DO
      IF (.NOT. ieee_is_finite(energy) .OR. .NOT. ALL(ieee_is_finite(gradient))) THEN
         ERROR STOP 'TPSS energy or first derivative is nonfinite at a stationary density point'
      END IF
      CALL xc_dset_release(derivatives)
      CALL xc_rho_set_release(rho_set)
   END SUBROUTINE evaluate_tpss_point

! **************************************************************************************************
!> \brief Derivative of an independent quadratic kinetic-density functional, with one or two spins.
! **************************************************************************************************
   SUBROUTINE check_kinetic_adjoint()
      REAL(dp), PARAMETER                                :: step = 1.e-5_dp

      INTEGER                                            :: ir, ns, p, q, side, spin
      REAL(dp) :: density(4, 4, 2), derivative, direction(4, 4, 2), error, &
         gradient(3, grid%na, 2), minus, plus, potential(4, 4, 2), rho(grid%na, 2), &
         tau(grid%na, 2), weighted(grid%na, 2)

      error = 0.0_dp
      DO spin = 1, 2
         DO q = 1, 4
            DO p = 1, 4
               density(p, q, spin) = SIN(0.3_dp*p*spin)*SIN(0.3_dp*q*spin)
               IF (p == q) density(p, q, spin) = density(p, q, spin) + 0.2_dp
               direction(p, q, spin) = COS(0.23_dp*(p + q)*spin)
            END DO
         END DO
      END DO
      DO ns = 1, 2
         DO side = 1, 2
            DO ir = 51, 151, 50
               CALL paw_xc_collocate(grid, onecenter, density(:, :, :ns), ir, side, &
                                     rho(:, :ns), gradient(:, :, :ns), tau(:, :ns))
               DO spin = 1, ns
                  weighted(:, spin) = 2.0_dp*tau(:, spin)*grid%wa*grid%wr(ir)
               END DO
               CALL paw_xc_kinetic_adjoint(grid, onecenter, ir, side, weighted(:, :ns), potential(:, :, :ns))
               derivative = SUM(potential(:, :, :ns)*direction(:, :, :ns))
               CALL paw_xc_collocate(grid, onecenter, density(:, :, :ns) + step*direction(:, :, :ns), ir, side, &
                                     rho(:, :ns), gradient(:, :, :ns), tau(:, :ns))
               plus = SUM(MATMUL(grid%wa, tau(:, :ns)**2))*grid%wr(ir)
               CALL paw_xc_collocate(grid, onecenter, density(:, :, :ns) - step*direction(:, :, :ns), ir, side, &
                                     rho(:, :ns), gradient(:, :, :ns), tau(:, :ns))
               minus = SUM(MATMUL(grid%wa, tau(:, :ns)**2))*grid%wr(ir)
               error = MAX(error, ABS((plus - minus)/(2.0_dp*step) - derivative))
            END DO
         END DO
      END DO
      WRITE (*, '(A,ES15.6)') 'PAW kinetic-density adjoint error:', error
      IF (error > 1.e-10_dp) ERROR STOP 'PAW kinetic potential is not the energy derivative'
   END SUBROUTINE check_kinetic_adjoint

! **************************************************************************************************
!> \brief Atomic magnetic averaging against analytic p-shell exchange and AO directional derivatives.
! **************************************************************************************************
   SUBROUTINE check_atomic_projection()
      REAL(dp), PARAMETER                                :: step = 1.e-5_dp

      INTEGER                                            :: l, ns, spin
      REAL(dp) :: amplitude, derivative, direction(2, 2, 0:1, 2), energy, errors(3), exchange, &
         hartree, integral, jmat(2, 2, 0:1), matrix(2, 2, 0:1, 2), minus, occupation(2), &
         p(2, 2, 0:1, 2), plus, prefactor, reference, scratch(2, 2, 0:1, 2), xmat(2, 2, 0:1, 2)
      TYPE(atom_basis_type)                              :: basis
      TYPE(atom_paw_type)                                :: cache
      TYPE(atom_potential_type)                          :: potential
      TYPE(section_type), POINTER                        :: xc_spec
      TYPE(section_vals_type), POINTER                   :: xc

      potential%ppot_type = upf_pseudo
      potential%upf_pot = pot
      potential%upf_pot%rho_nlcc(:) = 0.0_dp
      potential%upf_pot%paw%rho_core(:) = 0.0_dp
      CALL atom_paw_init(potential, cache)
      basis%nbas(0:1) = 2
      ALLOCATE (basis%uspp)
      ALLOCATE (basis%uspp%lbeta, SOURCE=[0, 1])
      ALLOCATE (basis%uspp%beta(2, 2))
      basis%uspp%beta(:, 1) = [0.8_dp, -0.3_dp]
      basis%uspp%beta(:, 2) = [0.6_dp, 0.4_dp]
      NULLIFY (xc_spec, xc)
      CALL create_xc_section(xc_spec)
      CALL section_vals_create(xc, xc_spec)
      CALL section_vals_val_set(xc, 'DENSITY_CUTOFF', r_val=1.e-20_dp)
      CALL section_vals_val_set(xc, 'GRADIENT_CUTOFF', r_val=1.e-20_dp)
      CALL section_vals_val_set(xc, 'XC_FUNCTIONAL%XALPHA%_SECTION_PARAMETERS_', l_val=.TRUE.)
      ! Spherical p-shell density: N*beta**2*r**2*exp(-exponent*r**2)/(4*pi).
      integral = 0.5_dp*GAMMA(17.0_dp/6.0_dp)*(4.0_dp*exponent/3.0_dp)**(-17.0_dp/6.0_dp)
      errors = 0.0_dp
      DO ns = 1, 2
         occupation = [1.3_dp, 0.0_dp]
         IF (ns == 2) occupation = [0.9_dp, 0.4_dp]
         p = 0.0_dp
         p(1, 1, 1, :) = occupation
         jmat = 0.0_dp
         xmat = 0.0_dp
         hartree = 0.0_dp
         exchange = 0.0_dp
         IF (ns == 1) THEN
            CALL atom_paw_add(cache, basis, xc, p(:, :, :, 1), hartree, exchange, jmat, xmat(:, :, :, 1))
         ELSE
            CALL atom_paw_add(cache, basis, xc, p(:, :, :, 1), hartree, exchange, jmat, xmat(:, :, :, 1), &
                              pmatb=p(:, :, :, 2), xcmatb=xmat(:, :, :, 2))
         END IF
         prefactor = -0.75_dp*(3.0_dp/pi)**(1.0_dp/3.0_dp)*REAL(ns, dp)**(1.0_dp/3.0_dp)
         reference = 0.0_dp
         DO spin = 1, ns
            amplitude = occupation(spin)*basis%uspp%beta(1, 2)**2/fourpi
            energy = fourpi*prefactor*integral*amplitude**(4.0_dp/3.0_dp)*(1.0_dp - pseudo_scale**(4.0_dp/3.0_dp))
            reference = reference + energy
            derivative = (4.0_dp/3.0_dp)*energy/occupation(spin)
            errors(2) = MAX(errors(2), ABS(xmat(1, 1, 1, spin) - derivative))
         END DO
         errors(1) = MAX(errors(1), ABS(exchange - reference))
      END DO
      IF (MAXVAL(errors(1:2)) > 2.e-10_dp) ERROR STOP 'Atomic PAW magnetic averaging differs from analytic p-shell exchange'

      CALL section_vals_val_set(xc, 'XC_FUNCTIONAL%XALPHA%_SECTION_PARAMETERS_', l_val=.FALSE.)
      CALL section_vals_val_set(xc, 'XC_FUNCTIONAL%PBE%_SECTION_PARAMETERS_', l_val=.TRUE.)
      DO ns = 1, 2
         DO spin = 1, ns
            DO l = 0, 1
               p(:, :, l, spin) = RESHAPE([0.8_dp, 0.13_dp, 0.13_dp, 0.5_dp], [2, 2])/REAL(spin + l, dp)
               direction(:, :, l, spin) = RESHAPE([0.2_dp, -0.07_dp, -0.07_dp, 0.09_dp], [2, 2])*REAL(spin + l, dp)
            END DO
         END DO
         CALL atomic_correction(cache, basis, xc, p(:, :, :, :ns), energy, matrix(:, :, :, :ns))
         CALL atomic_correction(cache, basis, xc, p(:, :, :, :ns) + step*direction(:, :, :, :ns), &
                                plus, scratch(:, :, :, :ns))
         CALL atomic_correction(cache, basis, xc, p(:, :, :, :ns) - step*direction(:, :, :, :ns), &
                                minus, scratch(:, :, :, :ns))
         derivative = SUM(matrix(:, :, :, :ns)*direction(:, :, :, :ns))
         errors(3) = MAX(errors(3), ABS((plus - minus)/(2.0_dp*step) - derivative))
      END DO
      WRITE (*, '(A,3ES15.6)') 'Atomic PAW exchange/derivative/AO errors:', errors
      IF (errors(3) > 2.e-9_dp) ERROR STOP 'Atomic PAW operator is not the correction energy derivative'
      CALL section_vals_release(xc)
      CALL section_release(xc_spec)
      CALL atom_release_upf(potential%upf_pot)
   END SUBROUTINE check_atomic_projection

! **************************************************************************************************
!> \brief Atomic gradient norms and matrix pullbacks for same-sign and opposing spin radial slopes.
!>        This functional consistency check uses analytic Gaussian s/p fields, not a basis benchmark.
! **************************************************************************************************
   SUBROUTINE check_atomic_gradient_norm()
      REAL(dp), PARAMETER :: step = 2.e-6_dp
#if defined(__LIBXC)
      INTEGER, PARAMETER :: nvariant = 3
#else
      INTEGER, PARAMETER :: nvariant = 2
#endif
      INTEGER :: i, j, l, ns, variant
      REAL(dp) :: a, derivative, energy, plus, minus, errors(3), &
                  p(2, 2, 0:1, 2), trial(2, 2, 0:1, 2), direction(2, 2, 0:1, 2), &
                  matrix(2, 2, 0:1, 2), scratch(2, 2, 0:1, 2), slopes(n, 2)
      TYPE(atom_basis_type), POINTER :: basis
      TYPE(section_type), POINTER :: xc_spec
      TYPE(section_vals_type), POINTER :: xc

      ALLOCATE (basis)
      ALLOCATE (basis%grid)
      basis%nbas(0:1) = 2
      basis%grid%nr = n
      ALLOCATE (basis%grid%rad(n), basis%grid%rad2(n), basis%grid%wr(n))
      basis%grid%rad(:) = [(0.01_dp*(REAL(i, dp) - 0.5_dp), i=1, n)]
      basis%grid%rad2(:) = basis%grid%rad**2
      basis%grid%wr(:) = 0.01_dp*basis%grid%rad2
      ALLOCATE (basis%bf(n, 2, 0:1), basis%dbf(n, 2, 0:1))
      DO j = 1, 2
         a = 0.5_dp + 0.3_dp*REAL(j, dp)
         DO l = 0, 1
            basis%bf(:, j, l) = basis%grid%rad**l*EXP(-a*basis%grid%rad2)
            IF (l == 0) THEN
               basis%dbf(:, j, l) = -2.0_dp*a*basis%grid%rad*basis%bf(:, j, l)
            ELSE
               basis%dbf(:, j, l) = (1.0_dp - 2.0_dp*a*basis%grid%rad2)*EXP(-a*basis%grid%rad2)
            END IF
         END DO
      END DO
      p(:, :, 0, 1) = RESHAPE([0.8_dp, 0.03_dp, 0.03_dp, 0.3_dp], [2, 2])
      p(:, :, 1, 1) = 0.1_dp*p(:, :, 0, 1)
      p(:, :, 0, 2) = 0.1_dp*p(:, :, 0, 1)
      p(:, :, 1, 2) = 2.0_dp*p(:, :, 0, 1)
      CALL atom_density(slopes(:, 1), p(:, :, :, 1), basis, 1, typ="DER")
      CALL atom_density(slopes(:, 2), p(:, :, :, 2), basis, 1, typ="DER")
      IF (.NOT. ANY(slopes(:, 1)*slopes(:, 2) < 0.0_dp)) ERROR STOP 'Atomic norm fixture lacks opposing spin slopes'
      DO ns = 1, 2
         DO l = 0, 1
            direction(:, :, l, ns) = REAL(l + ns, dp)*RESHAPE([0.13_dp, -0.02_dp, -0.02_dp, 0.07_dp], [2, 2])
         END DO
      END DO
      NULLIFY (xc_spec, xc)
      CALL create_xc_section(xc_spec)
      errors = 0.0_dp
      DO variant = 1, nvariant
         CALL section_vals_create(xc, xc_spec)
         IF (variant == 1) CALL section_vals_val_set(xc, 'XC_FUNCTIONAL%PBE%_SECTION_PARAMETERS_', l_val=.TRUE.)
         IF (variant == 2) CALL section_vals_val_set(xc, 'XC_FUNCTIONAL%TPSS%_SECTION_PARAMETERS_', l_val=.TRUE.)
#if defined(__LIBXC)
         IF (variant == 3) THEN
            CALL section_vals_val_set(xc, 'XC_FUNCTIONAL%MGGA_X_TPSS%_SECTION_PARAMETERS_', l_val=.TRUE.)
            CALL section_vals_val_set(xc, 'XC_FUNCTIONAL%MGGA_C_TPSS%_SECTION_PARAMETERS_', l_val=.TRUE.)
         END IF
#endif
         DO ns = 1, 2
            IF (variant == 2 .AND. ns == 2) CYCLE
            trial = p
            IF (ns == 1) trial(:, :, :, 1) = SUM(p, DIM=4)
            CALL atomic_smooth_xc(basis, xc, trial(:, :, :, :ns), energy, matrix(:, :, :, :ns))
            IF (ABS(energy) < 1.e-8_dp) ERROR STOP 'Atomic norm test did not evaluate an XC functional'
            CALL atomic_smooth_xc(basis, xc, trial(:, :, :, :ns) + step*direction(:, :, :, :ns), &
                                  plus, scratch(:, :, :, :ns))
            CALL atomic_smooth_xc(basis, xc, trial(:, :, :, :ns) - step*direction(:, :, :, :ns), &
                                  minus, scratch(:, :, :, :ns))
            derivative = SUM(matrix(:, :, :, :ns)*direction(:, :, :, :ns))
            errors(variant) = MAX(errors(variant), ABS((plus - minus)/(2.0_dp*step) - derivative))
         END DO
         CALL section_vals_release(xc)
      END DO
      WRITE (*, '(A,3ES15.6)') 'Atomic gradient-norm PBE/TPSS/Libxc matrix errors:', errors
      IF (MAXVAL(errors) > 2.e-7_dp) ERROR STOP 'Atomic gradient-norm matrix is not the XC energy derivative'
      CALL section_release(xc_spec)
      DEALLOCATE (basis%bf, basis%dbf, basis%grid%rad, basis%grid%rad2, basis%grid%wr)
      DEALLOCATE (basis%grid)
      DEALLOCATE (basis)
   END SUBROUTINE check_atomic_gradient_norm

! **************************************************************************************************
!> \brief Smooth atomic functional and Gaussian matrix for either spin convention.
!> \param basis Analytic Gaussian radial fixture and its quadrature.
!> \param xc Atomic XC section.
!> \param density Physical spin density matrices.
!> \param energy Smooth XC energy.
!> \param matrix Its Gaussian density-matrix derivative.
! **************************************************************************************************
   SUBROUTINE atomic_smooth_xc(basis, xc, density, energy, matrix)
      TYPE(atom_basis_type), POINTER                     :: basis
      TYPE(section_vals_type), POINTER                   :: xc
      REAL(dp), INTENT(IN)                               :: density(:, :, 0:, :)
      REAL(dp), INTENT(OUT)                              :: energy, matrix(:, :, 0:, :)

      matrix = 0.0_dp
      IF (SIZE(density, 4) == 1) THEN
         CALL atom_vxc_lda(basis, density(:, :, :, 1), 1, xc, energy, matrix(:, :, :, 1))
      ELSE
         CALL atom_vxc_lsd(basis, density(:, :, :, 1), density(:, :, :, 2), 1, xc, energy, &
                           matrix(:, :, :, 1), matrix(:, :, :, 2))
      END IF
   END SUBROUTINE atomic_smooth_xc

! **************************************************************************************************
!> \brief Evaluate atomic correction and its complete AO derivative for either spin convention.
!> \param cache Prepared PAW data.
!> \param basis Radial basis with beta overlaps.
!> \param xc XC input section.
!> \param density Occupation-summed atomic density matrices.
!> \param energy Hartree plus XC correction energy.
!> \param matrix Derivative with respect to the atomic density matrices.
! **************************************************************************************************
   SUBROUTINE atomic_correction(cache, basis, xc, density, energy, matrix)
      TYPE(atom_paw_type), INTENT(IN)                    :: cache
      TYPE(atom_basis_type), INTENT(IN)                  :: basis
      TYPE(section_vals_type), POINTER                   :: xc
      REAL(dp), INTENT(IN)                               :: density(:, :, 0:, :)
      REAL(dp), INTENT(OUT)                              :: energy, matrix(:, :, 0:, :)

      INTEGER                                            :: spin
      REAL(dp) :: exchange, hartree, &
         coulomb(SIZE(density, 1), SIZE(density, 2), 0:UBOUND(density, 3))

      hartree = 0.0_dp
      exchange = 0.0_dp
      coulomb = 0.0_dp
      matrix = 0.0_dp
      IF (SIZE(density, 4) == 1) THEN
         CALL atom_paw_add(cache, basis, xc, density(:, :, :, 1), hartree, exchange, coulomb, matrix(:, :, :, 1))
      ELSE
         CALL atom_paw_add(cache, basis, xc, density(:, :, :, 1), hartree, exchange, coulomb, matrix(:, :, :, 1), &
                           pmatb=density(:, :, :, 2), xcmatb=matrix(:, :, :, 2))
      END IF
      energy = hartree + exchange
      DO spin = 1, SIZE(density, 4)
         matrix(:, :, :, spin) = matrix(:, :, :, spin) + coulomb
      END DO
   END SUBROUTINE atomic_correction

! **************************************************************************************************
!> \brief Optional original PAW files: polarized XC and its nonspherical projector-density derivative.
! **************************************************************************************************
   SUBROUTINE original_datasets()
      REAL(dp), PARAMETER                                :: step = 1.e-3_dp

      CHARACTER(LEN=default_path_length)                 :: filename
      INTEGER                                            :: argument, k, np, p, q, spin
      LOGICAL                                            :: source_root_argument
      REAL(dp)                                           :: derivative, energy, error, minus, plus
      REAL(dp), ALLOCATABLE                              :: density(:, :, :), direction(:, :, :), &
                                                            potential(:, :, :), scratch(:, :, :)
      TYPE(atom_upfpot_type)                             :: original
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(paw_onecenter_type)                           :: centered
      TYPE(paw_xc_type)                                  :: xcgrid
      TYPE(upf_projector_type)                           :: native

      logger => cp_get_default_logger()
      DO argument = 1, COMMAND_ARGUMENT_COUNT()
         CALL GET_COMMAND_ARGUMENT(argument, filename)
         IF (argument == 1) THEN
            INQUIRE (FILE=TRIM(filename)//'/tests/UNIT_TESTS', EXIST=source_root_argument)
            IF (source_root_argument) CYCLE
         END IF
         CALL atom_read_upf(original, TRIM(filename), read_paw=.TRUE.)
         CALL section_vals_val_set(functional, 'PBE%PARAMETRIZATION', &
                                   i_val=MERGE(xc_pbe_sol, xc_pbe_orig, INDEX(original%functional, 'PSX') > 0))
         CALL upf_projector_init(original, MAXVAL(original%lbeta), native)
         CALL paw_onecenter_init(original, native, centered)
         CALL paw_xc_init(original, 110, xcgrid)
         np = native%nproj
         ALLOCATE (density(np, np, 2), direction(np, np, 2), potential(np, np, 2), scratch(np, np, 2))
         density = 0.0_dp
         DO p = 1, np
            k = native%radial_index(p)
            density(p, p, 1) = MAX(0.2_dp, original%paw%occupations(k))/REAL(2*original%lbeta(k) + 1, dp)
         END DO
         DO q = 1, np
            DO p = 1, np
               DO spin = 1, 2
                  direction(p, q, spin) = 0.02_dp*COS(0.3_dp*(p + q)*spin)*SQRT(density(p, p, 1)*density(q, q, 1))
               END DO
            END DO
         END DO
         density(:, :, 2) = 0.35_dp*density(:, :, 1)
         density(:, :, 1) = 0.65_dp*density(:, :, 1)
         density(:, :, :) = density + direction
         CALL paw_xc_evaluate(xcgrid, centered, functional, density, 1.e-14_dp, 1.e-14_dp, energy, potential)
         IF (.NOT. ieee_is_finite(energy) .OR. .NOT. ALL(ieee_is_finite(potential))) THEN
            ERROR STOP 'Nonfinite original-dataset PAW XC correction'
         END IF
         CALL paw_xc_evaluate(xcgrid, centered, functional, density + step*direction, &
                              1.e-14_dp, 1.e-14_dp, plus, scratch)
         CALL paw_xc_evaluate(xcgrid, centered, functional, density - step*direction, &
                              1.e-14_dp, 1.e-14_dp, minus, scratch)
         derivative = SUM(potential*direction)
         error = ABS((plus - minus)/(2.0_dp*step) - derivative)/MAX(1.0_dp, ABS(derivative))
         IF (logger%para_env%is_source()) THEN
            WRITE (*, '(A,1X,A,3(1X,ES25.16))') 'PAW XC:', TRIM(filename), energy, derivative, error
         END IF
         IF (error > 5.e-7_dp) ERROR STOP 'Original-dataset PAW XC matrix is not the energy derivative'
         IF (xcgrid%core_kinetic%available) CALL original_core_tau_matrices(xcgrid, centered, density, direction)
         DEALLOCATE (density, direction, potential, scratch)
         CALL atom_release_upf(original)
      END DO
   END SUBROUTINE original_datasets

! **************************************************************************************************
!> \brief Original-dataset PAW TPSS matrix derivatives with the automatically supplied shared core.
!> \param xcgrid Original PAW radial cache.
!> \param centered Original one-center projection.
!> \param density Positive two-spin projector density.
!> \param direction Nonspherical two-spin variation.
! **************************************************************************************************
   SUBROUTINE original_core_tau_matrices(xcgrid, centered, density, direction)
      TYPE(paw_xc_type), INTENT(IN) :: xcgrid
      TYPE(paw_onecenter_type), INTENT(IN) :: centered
      REAL(dp), INTENT(IN) :: density(:, :, :), direction(:, :, :)
      REAL(dp), PARAMETER :: step = 1.e-4_dp
#if defined(__LIBXC)
      INTEGER, PARAMETER :: nvariant = 2
#else
      INTEGER, PARAMETER :: nvariant = 1
#endif
      INTEGER :: ns
      REAL(dp) :: energy, plus, minus, derivative, error
      REAL(dp) :: trial(SIZE(density, 1), SIZE(density, 2), 2), &
                  variation(SIZE(density, 1), SIZE(density, 2), 2), &
                  potential(SIZE(density, 1), SIZE(density, 2), 2), &
                  scratch(SIZE(density, 1), SIZE(density, 2), 2)
      TYPE(cp_logger_type), POINTER :: logger

      logger => cp_get_default_logger()
      DO ns = 1, nvariant
         CALL section_vals_release(functional)
         CALL section_vals_create(functional, specification)
         trial = density
         variation = direction
         IF (ns == 1) THEN
            CALL section_vals_val_set(functional, 'TPSS%_SECTION_PARAMETERS_', l_val=.TRUE.)
            trial(:, :, 1) = SUM(density, DIM=3)
            variation(:, :, 1) = SUM(direction, DIM=3)
         END IF
#if defined(__LIBXC)
         IF (ns == 2) THEN
            CALL section_vals_val_set(functional, 'MGGA_X_TPSS%_SECTION_PARAMETERS_', l_val=.TRUE.)
            CALL section_vals_val_set(functional, 'MGGA_C_TPSS%_SECTION_PARAMETERS_', l_val=.TRUE.)
         END IF
#endif
         CALL paw_xc_evaluate(xcgrid, centered, functional, trial(:, :, :ns), &
                              1.e-14_dp, 1.e-14_dp, energy, potential(:, :, :ns))
         CALL paw_xc_evaluate(xcgrid, centered, functional, trial(:, :, :ns) + step*variation(:, :, :ns), &
                              1.e-14_dp, 1.e-14_dp, plus, scratch(:, :, :ns))
         CALL paw_xc_evaluate(xcgrid, centered, functional, trial(:, :, :ns) - step*variation(:, :, :ns), &
                              1.e-14_dp, 1.e-14_dp, minus, scratch(:, :, :ns))
         derivative = SUM(potential(:, :, :ns)*variation(:, :, :ns))
         error = ABS((plus - minus)/(2.0_dp*step) - derivative)/MAX(1.0_dp, ABS(derivative))
         IF (logger%para_env%is_source()) WRITE (*, '(A,I2,3ES25.16)') &
            'Original shared-core PAW TPSS:', ns, energy, derivative, error
         IF (error > 5.e-7_dp) ERROR STOP 'Original shared-core PAW TPSS matrix is not the energy derivative'
      END DO
      CALL section_vals_release(functional)
      CALL section_vals_create(functional, specification)
      CALL section_vals_val_set(functional, 'PBE%_SECTION_PARAMETERS_', l_val=.TRUE.)
   END SUBROUTINE original_core_tau_matrices

! **************************************************************************************************
!> \brief Independent Cartesian s/p polynomial density and gradient at every angular direction.
! **************************************************************************************************
   SUBROUTINE check_cartesian_fields()
      INTEGER                                            :: a, axis, ir, p, q, shell, side, spin
      REAL(dp) :: basis(4), core, dbasis(4, 3), density(4, 4, 2), exact, exact_gradient(3), &
         factor, field_error, gradient(3, grid%na, 2), gradient_error, point(3), polynomial, &
         rho(grid%na, 2), scale, tau(grid%na, 2), tau_error, tau_reference, tau_scale, &
         wave_gradient(4)

      DO spin = 1, 2
         DO q = 1, 4
            DO p = 1, 4
               density(p, q, spin) = SIN(0.3_dp*p*spin)*SIN(0.3_dp*q*spin)
               IF (p == q) density(p, q, spin) = density(p, q, spin) + 0.2_dp
            END DO
         END DO
      END DO
      dbasis = 0.0_dp
      dbasis(2, 2) = SQRT(3.0_dp)
      dbasis(3, 3) = SQRT(3.0_dp)
      dbasis(4, 1) = SQRT(3.0_dp)
      field_error = 0.0_dp
      gradient_error = 0.0_dp
      tau_error = 0.0_dp
      DO shell = 1, 3
         ir = 1 + 50*shell
         DO side = 1, 2
            scale = MERGE(1.0_dp, pseudo_scale, side == 1)
            core = MERGE(0.2_dp, 0.1_dp, side == 1)/2.0_dp
            tau_scale = MERGE(1.0_dp, 0.8_dp**2, side == 1)
            CALL paw_xc_collocate(grid, onecenter, density, ir, side, rho, gradient, tau)
            DO a = 1, grid%na
               point = grid%r(ir)*grid%directions(:, a)
               basis = [1.0_dp, SQRT(3.0_dp)*point(2), SQRT(3.0_dp)*point(3), SQRT(3.0_dp)*point(1)]
               factor = EXP(-exponent*SUM(point**2))
               DO spin = 1, 2
                  polynomial = DOT_PRODUCT(basis, MATMUL(density(:, :, spin), basis))
                  exact = factor*(scale*polynomial/fourpi + core)
                  tau_reference = 0.0_dp
                  DO axis = 1, 3
                     exact_gradient(axis) = scale*factor*2.0_dp* &
                                            DOT_PRODUCT(dbasis(:, axis), MATMUL(density(:, :, spin), basis))/fourpi &
                                            - 2.0_dp*exponent*point(axis)*exact
                     wave_gradient = dbasis(:, axis) - exponent*point(axis)*basis
                     tau_reference = tau_reference + 0.5_dp*tau_scale*factor* &
                                     DOT_PRODUCT(wave_gradient, MATMUL(density(:, :, spin), wave_gradient))/fourpi
                  END DO
                  field_error = MAX(field_error, ABS(rho(a, spin) - exact))
                  gradient_error = MAX(gradient_error, MAXVAL(ABS(gradient(:, a, spin) - exact_gradient)))
                  tau_error = MAX(tau_error, ABS(tau(a, spin) - tau_reference))
               END DO
            END DO
         END DO
      END DO
      WRITE (*, '(A,2ES15.6)') 'PAW XC Cartesian density/gradient errors:', field_error, gradient_error
      IF (field_error > 2.e-13_dp) ERROR STOP 'PAW XC density differs from Cartesian Gaussian result'
      IF (gradient_error > 2.e-8_dp) ERROR STOP 'PAW XC gradient differs from Cartesian Gaussian result'
      WRITE (*, '(A,ES15.6)') 'PAW Cartesian partial-wave kinetic-density error:', tau_error
      IF (tau_error > 2.e-8_dp) ERROR STOP 'PAW kinetic density differs from Cartesian Gaussian result'
   END SUBROUTINE check_cartesian_fields

! **************************************************************************************************
!> \brief Closed-form Gaussian Dirac exchange with core and compensation, unpolarized and polarized.
! **************************************************************************************************
   SUBROUTINE check_local_exchange()
      INTEGER                                            :: nspin, spin
      REAL(dp)                                           :: a, b, density(4, 4, 2), derivative, &
                                                            energy, integral, occupation(2), &
                                                            potential(4, 4, 2), prefactor, &
                                                            reference

      CALL section_vals_val_set(functional, 'XALPHA%_SECTION_PARAMETERS_', l_val=.TRUE.)
      integral = (pi/(4.0_dp*exponent/3.0_dp))**1.5_dp
      DO nspin = 1, 2
         density = 0.0_dp
         occupation = [1.0_dp, 0.0_dp]
         IF (nspin == 2) occupation = [0.7_dp, 0.3_dp]
         DO spin = 1, nspin
            density(1, 1, spin) = occupation(spin)
         END DO
         CALL paw_xc_evaluate(grid, onecenter, functional, density(:, :, :nspin), 1.e-20_dp, 1.e-20_dp, &
                              energy, potential(:, :, :nspin))
         prefactor = -0.75_dp*(3.0_dp/pi)**(1.0_dp/3.0_dp)*REAL(nspin, dp)**(1.0_dp/3.0_dp)
         reference = 0.0_dp
         DO spin = 1, nspin
            a = occupation(spin)/fourpi + 0.2_dp/REAL(nspin, dp)
            b = occupation(spin)*pseudo_scale/fourpi + 0.1_dp/REAL(nspin, dp)
            reference = reference + prefactor*integral*(a**(4.0_dp/3.0_dp) - b**(4.0_dp/3.0_dp))
            derivative = prefactor*integral*(4.0_dp/3.0_dp)/fourpi*(a**(1.0_dp/3.0_dp) - &
                                                                    pseudo_scale*b**(1.0_dp/3.0_dp))
            IF (ABS(potential(1, 1, spin) - derivative) > 2.e-11_dp) ERROR STOP 'Wrong analytic PAW exchange derivative'
         END DO
         IF (ABS(energy - reference) > 2.e-11_dp) ERROR STOP 'Wrong analytic PAW exchange energy'
      END DO
      CALL section_vals_val_set(functional, 'XALPHA%_SECTION_PARAMETERS_', l_val=.FALSE.)
   END SUBROUTINE check_local_exchange

! **************************************************************************************************
!> \brief PBE/PBEsol variational checks, closed/open-shell equivalence and spin permutation.
! **************************************************************************************************
   SUBROUTINE check_gga()
      INTEGER                                            :: p, q, variant
      REAL(dp) :: check_energy, density(4, 4, 2), derivative, direction(4, 4, 2), energy, minus, &
         plus, potential(4, 4, 2), reference_energy, restricted(4, 4, 1), &
         restricted_potential(4, 4, 1), scratch(4, 4, 2), step

      CALL section_vals_val_set(functional, 'PBE%_SECTION_PARAMETERS_', l_val=.TRUE.)
      DO q = 1, 4
         DO p = 1, 4
            restricted(p, q, 1) = SIN(0.31_dp*p)*SIN(0.31_dp*q)
            IF (p == q) restricted(p, q, 1) = restricted(p, q, 1) + 0.3_dp
            direction(p, q, 1) = COS(0.23_dp*(p + q))
            direction(p, q, 2) = SIN(0.37_dp*(p + q))
         END DO
      END DO
      step = 1.e-5_dp
      reference_energy = 0.0_dp
      DO variant = 1, 2
         CALL section_vals_val_set(functional, 'PBE%PARAMETRIZATION', i_val=MERGE(xc_pbe_orig, xc_pbe_sol, variant == 1))
         CALL paw_xc_evaluate(grid, onecenter, functional, restricted, 1.e-20_dp, 1.e-20_dp, energy, restricted_potential)
         IF (variant == 2 .AND. ABS(energy - reference_energy) < 1.e-6_dp) ERROR STOP 'GGA parametrization did not affect PAW XC'
         reference_energy = energy
         density(:, :, 1) = 0.5_dp*restricted(:, :, 1)
         density(:, :, 2) = density(:, :, 1)
         CALL paw_xc_evaluate(grid, onecenter, functional, density, 1.e-20_dp, 1.e-20_dp, check_energy, potential)
         IF (ABS(check_energy - energy) > 2.e-11_dp) ERROR STOP 'Restricted and unpolarized spin PAW XC differ'
         IF (MAXVAL(ABS(potential(:, :, 1) - restricted_potential(:, :, 1))) > 2.e-10_dp) THEN
            ERROR STOP 'Restricted and unpolarized spin PAW XC derivatives differ'
         END IF
         density(:, :, 1) = 0.7_dp*restricted(:, :, 1)
         density(:, :, 2) = 0.3_dp*restricted(:, :, 1)
         CALL paw_xc_evaluate(grid, onecenter, functional, density, 1.e-20_dp, 1.e-20_dp, energy, potential)
         CALL paw_xc_evaluate(grid, onecenter, functional, density + step*direction, 1.e-20_dp, 1.e-20_dp, plus, scratch)
         CALL paw_xc_evaluate(grid, onecenter, functional, density - step*direction, 1.e-20_dp, 1.e-20_dp, minus, scratch)
         derivative = SUM(potential*direction)
   WRITE (*, '(A,I2,ES15.6)') 'PAW GGA projector-density derivative error:', variant, ABS((plus - minus)/(2.0_dp*step) - derivative)
         IF (ABS((plus - minus)/(2.0_dp*step) - derivative) > 3.e-8_dp) ERROR STOP 'PAW GGA XC matrix is not the energy derivative'
         CALL paw_xc_evaluate(grid, onecenter, functional, density(:, :, 2:1:-1), 1.e-20_dp, 1.e-20_dp, check_energy, scratch)
         IF (ABS(check_energy - energy) > 2.e-12_dp) ERROR STOP 'PAW XC spin exchange changed energy'
         IF (MAXVAL(ABS(scratch - potential(:, :, 2:1:-1))) > 2.e-11_dp) ERROR STOP 'PAW XC spin derivative permutation failed'
      END DO
   END SUBROUTINE check_gga

END PROGRAM paw_xc_unittest
