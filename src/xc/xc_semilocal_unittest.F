!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM xc_semilocal_unittest
   USE ieee_arithmetic, ONLY: ieee_is_finite
   USE input_cp2k_xc, ONLY: create_xc_fun_section
   USE input_section_types, ONLY: section_release, &
                                  section_type, &
                                  section_vals_create, &
                                  section_vals_release, &
                                  section_vals_type, &
                                  section_vals_val_set
   USE kinds, ONLY: dp
   USE libcp2k, ONLY: cp2k_finalize, &
                      cp2k_init
   USE mathconstants, ONLY: pi
   USE skala_density_constraints, ONLY: skala_density_constraints_type, skala_density_constraints_set, &
                                        skala_density_constraints_reset, skala_fft_domain, skala_last_sample_domain
   USE xc_atom, ONLY: xc_rho_set_atom_update
   USE xc_derivative_desc, ONLY: deriv_norm_drho, &
                                 deriv_norm_drhoa, &
                                 deriv_norm_drhob, &
                                 deriv_rho, &
                                 deriv_rhoa, &
                                 deriv_rhob, &
                                 deriv_tau
   USE xc_derivative_set_types, ONLY: xc_derivative_set_type, &
                                      xc_dset_create, &
                                      xc_dset_get_derivative, &
                                      xc_dset_release
   USE xc_derivative_types, ONLY: xc_derivative_get, &
                                  xc_derivative_type
   USE xc_derivatives, ONLY: xc_functionals_eval, &
                             xc_functionals_get_needs, xc_functionals_select_floor_limits, xc_functionals_spin_floors
   USE xc_input_constants, ONLY: xc_pbe_orig, &
                                 xc_pbe_rev, &
                                 xc_pbe_sol
   USE xc_pbe, ONLY: pbe_spin_density_floor
   USE xc_rho_cflags_types, ONLY: xc_rho_cflags_type
   USE xc_rho_set_types, ONLY: xc_rho_set_create, &
                               xc_rho_set_release, &
                               xc_rho_set_type
#if defined (__LIBXC)
   USE xc_derivative_desc, ONLY: deriv_gamma, deriv_gamma_aa, deriv_gamma_ab, deriv_gamma_bb, &
                                 deriv_tau_a, &
                                 deriv_tau_b
   USE xc_libxc, ONLY: libxc_density_floor, libxc_spin_polarized_eval, &
                       libxc_spin_unpolarized_eval
#endif

   IMPLICIT NONE
   INTEGER, PARAMETER :: labels(5) = [deriv_rhoa, deriv_rhob, deriv_norm_drho, deriv_norm_drhoa, deriv_norm_drhob]
   INTEGER, PARAMETER :: variants(3) = [xc_pbe_orig, xc_pbe_rev, xc_pbe_sol]
   TYPE(section_type), POINTER :: specification
   TYPE(section_vals_type), POINTER :: functional
   INTEGER :: axis, configuration, spin, variant
   REAL(dp) :: energy, gradient(5), hessian(5, 5), minus, plus, gm(5), gp(5), scratch(5, 5), &
               step, trial(5), x(5), energy_error, hessian_error

   CALL cp2k_init()
   NULLIFY (specification, functional)
   CALL create_xc_fun_section(specification)
   CALL section_vals_create(functional, specification)
   CALL section_vals_val_set(functional, 'PBE%_SECTION_PARAMETERS_', l_val=.TRUE.)
   energy_error = 0.0_dp
   hessian_error = 0.0_dp
   DO variant = 1, SIZE(variants)
      CALL section_vals_val_set(functional, 'PBE%PARAMETRIZATION', i_val=variants(variant))
      DO configuration = 1, 3
         DO spin = 1, 2
            IF (configuration == 1 .AND. spin == 2) CYCLE
            x = [0.4_dp, 0.2_dp, 0.3_dp, 0.2_dp, 0.25_dp]
            IF (configuration == 2) x(spin) = -0.1_dp
            IF (configuration == 3) x(spin) = 1.e-13_dp
            CALL evaluate(x, energy, gradient, hessian)
            DO axis = 1, 5
               step = 1.e-5_dp
               IF (configuration == 3 .AND. axis == spin) step = 1.e-14_dp
               trial = x
               trial(axis) = x(axis) + step
               CALL evaluate(trial, plus, gp, scratch)
               trial(axis) = x(axis) - step
               CALL evaluate(trial, minus, gm, scratch)
               energy_error = MAX(energy_error, ABS((plus - minus)/(2.0_dp*step) - gradient(axis)))
               hessian_error = MAX(hessian_error, MAXVAL(ABS((gp - gm)/(2.0_dp*step) - hessian(:, axis))))
            END DO
         END DO
      END DO
   END DO
   WRITE (*, '(A,2ES15.6)') 'PBE spin-floor first/second derivative errors:', energy_error, hessian_error
   IF (energy_error > 2.e-8_dp) ERROR STOP 'PBE potential is not the derivative of its clamped-density energy'
   IF (hessian_error > 2.e-7_dp) ERROR STOP 'PBE kernel is not the derivative of its clamped-density potential'
   CALL check_pbe_floor_limits()
   CALL check_xc_floor_selection()
   CALL section_vals_val_set(functional, 'PBE%_SECTION_PARAMETERS_', l_val=.FALSE.)
   CALL section_vals_val_set(functional, 'TPSS%_SECTION_PARAMETERS_', l_val=.TRUE.)
   CALL check_tpss_zero_gradient()
   CALL section_vals_val_set(functional, 'TPSS%_SECTION_PARAMETERS_', l_val=.FALSE.)
   CALL section_vals_release(functional)
   CALL section_vals_create(functional, specification)
   CALL section_vals_val_set(functional, 'XWPBE%_SECTION_PARAMETERS_', l_val=.TRUE.)
   CALL check_xwpbe_limits()
   CALL section_vals_val_set(functional, 'XWPBE%_SECTION_PARAMETERS_', l_val=.FALSE.)
#if defined (__LIBXC)
   CALL check_libxc_tau_projection()
   CALL check_libxc_c2_window()
   CALL check_libxc_c2_rounding()
   CALL check_libxc_individual_floor()
   CALL check_libxc_floor_limits()
   CALL check_libxc_spin_powers()
#endif
   CALL section_vals_release(functional)
   CALL section_release(specification)
   CALL cp2k_finalize()

CONTAINS

#if defined (__LIBXC)
! **************************************************************************************************
!> \brief Check meta-GGA potentials and kernels against independent differences of the projected energy.
! **************************************************************************************************
   SUBROUTINE check_libxc_tau_projection()
      CHARACTER(len=32), PARAMETER :: names(6) = [CHARACTER(len=32):: 'MGGA_X_TPSS', 'MGGA_C_TPSS',&
         'MGGA_X_R2SCAN', 'MGGA_C_R2SCAN', 'LDA_C_PW', 'GGA_C_PBE']
      REAL(dp), PARAMETER                                :: step = 1.e-5_dp

      INTEGER                                            :: axis, configuration, family, nconfig, &
                                                            nvar, spin_case
      REAL(dp)                                           :: energy, first_error, minus, plus, &
                                                            second_error
      REAL(dp), ALLOCATABLE                              :: gm(:), gp(:), gradient(:), &
                                                            hessian(:, :), scratch(:, :), &
                                                            trial(:), x(:)

      first_error = 0.0_dp
      second_error = 0.0_dp
      DO spin_case = 1, 2
         nvar = 3
         nconfig = 3
         IF (spin_case == 2) THEN
            nvar = 7
            nconfig = 6
         END IF
         ALLOCATE (x(nvar), trial(nvar), gradient(nvar), gp(nvar), gm(nvar), &
                   hessian(nvar, nvar), scratch(nvar, nvar))
         DO family = 1, SIZE(names)
            DO configuration = 1, nconfig
               IF (spin_case == 1) THEN
                  x(:) = [0.4_dp, 0.8_dp, 0.3_dp]
                  IF (configuration == 2) x(3) = 0.02_dp
                  IF (configuration == 3) x(3) = -0.02_dp
               ELSE
                  x(:) = [0.3_dp, 0.2_dp, 0.6_dp, 0.4_dp, 0.3_dp, 0.12_dp, 0.10_dp]
                  IF (configuration == 2) x(6:7) = [0.02_dp, 0.015_dp]
                  IF (configuration == 3) x(6) = -0.02_dp
                  IF (configuration == 4) x(6:7) = [-0.02_dp, -0.015_dp]
                  IF (configuration == 5) x(1) = -0.1_dp
                  IF (configuration == 6) x(2) = -0.1_dp
               END IF
               CALL evaluate_libxc_point(names(family), x, energy, gradient, hessian)
               IF (configuration >= 3) THEN
                  IF (INDEX(names(family), 'MGGA_X_') == 1 .AND. ABS(energy) < 1.e-4_dp) THEN
                     ERROR STOP 'Libxc discarded finite density at nonpositive tau'
                  END IF
                  trial(:) = x(:)
                  trial(nvar - spin_case + 1) = 1.e-10_dp
                  CALL evaluate_libxc_point(names(family), trial, plus, gp, scratch)
                  trial(nvar - spin_case + 1) = -1.e-10_dp
                  CALL evaluate_libxc_point(names(family), trial, minus, gm, scratch)
                  IF (ABS(plus - minus) > 1.e-9_dp) ERROR STOP 'Libxc energy is discontinuous through zero tau'
               END IF
               DO axis = 1, nvar
                  trial(:) = x(:)
                  trial(axis) = x(axis) + step
                  CALL evaluate_libxc_point(names(family), trial, plus, gp, scratch)
                  trial(axis) = x(axis) - step
                  CALL evaluate_libxc_point(names(family), trial, minus, gm, scratch)
                  first_error = MAX(first_error, &
                                    ABS((plus - minus)/(2.0_dp*step) - gradient(axis))/(1.0_dp + ABS(gradient(axis))))
                  second_error = MAX(second_error, &
                                     MAXVAL(ABS((gp - gm)/(2.0_dp*step) - hessian(:, axis))/ &
                                            (1.0_dp + ABS(hessian(:, axis)))))
               END DO
            END DO
         END DO
         DEALLOCATE (x, trial, gradient, gp, gm, hessian, scratch)
      END DO
      WRITE (*, '(A,2ES15.6)') 'Libxc projected-tau first/second derivative errors:', first_error, second_error
      IF (first_error > 3.e-8_dp) ERROR STOP 'Libxc projected-tau potential is not the energy derivative'
      IF (second_error > 2.e-6_dp) ERROR STOP 'Libxc projected-tau kernel is not the potential derivative'
   END SUBROUTINE check_libxc_tau_projection

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE check_libxc_c2_window()
      CHARACTER(len=32), PARAMETER :: names(4) = [CHARACTER(len=32):: 'MGGA_X_TPSS', 'MGGA_C_TPSS',&
         'MGGA_X_R2SCAN', 'MGGA_C_R2SCAN']
      REAL(dp), PARAMETER :: offsets(3) = [-0.6_dp, 0.2_dp, 0.6_dp], &
         relative_width = 64.0_dp*SQRT(EPSILON(0.0_dp))

      INTEGER                                            :: axis, configuration, family, nvar, &
                                                            spin_case
      REAL(dp)                                           :: em, emm, energy, ep, epp, first_error, &
                                                            first_fd, h, second_error, &
                                                            third_error, wa, wb, widtha, widthb
      REAL(dp), ALLOCATABLE :: gm(:), gmm(:), gp(:), gpp(:), gradient(:), hessian(:, :), hm(:, :), &
         hmm(:, :), hp(:, :), hpp(:, :), third(:, :, :), trial(:), x(:)

! Nonzero offsets keep the third-derivative differences well conditioned.

      first_error = 0.0_dp
      second_error = 0.0_dp
      third_error = 0.0_dp
      DO spin_case = 1, 2
         nvar = 3
         IF (spin_case == 2) nvar = 7
         ALLOCATE (x(nvar), trial(nvar), gradient(nvar), gm(nvar), gp(nvar), &
                   gmm(nvar), gpp(nvar), hessian(nvar, nvar), &
                   third(nvar, nvar, nvar), hm(nvar, nvar), hp(nvar, nvar), hmm(nvar, nvar), hpp(nvar, nvar))
         DO family = 1, SIZE(names)
            DO configuration = 1, SIZE(offsets)
               IF (spin_case == 1) THEN
                  x(:) = [0.4_dp, 0.64_dp, 0.3_dp]
                  wa = x(2)/(8.0_dp*x(1))
                  widtha = relative_width*wa
                  x(3) = wa + offsets(configuration)*widtha
               ELSE
                  x(:) = [0.3_dp, 0.2_dp, 0.055_dp, 0.16_dp, 0.09_dp, 0.12_dp, 0.10_dp]
                  wa = x(4)/(8.0_dp*x(1))
                  wb = x(5)/(8.0_dp*x(2))
                  widtha = relative_width*wa
                  widthb = relative_width*wb
                  x(6) = wa + offsets(configuration)*widtha
                  x(7) = wb - offsets(configuration)*widthb
               END IF
               ! The installed r2SCAN functions expose only first and second derivatives.
               IF (family <= 2) THEN
                  CALL evaluate_libxc_point(names(family), x, energy, gradient, hessian, third_all=third, contractions=.TRUE.)
               ELSE
                  CALL evaluate_libxc_point(names(family), x, energy, gradient, hessian, contractions=.TRUE.)
               END IF
               DO axis = 1, nvar
                  h = 0.01_dp*widtha
                  IF (spin_case == 1) THEN
                     IF (axis == 1) h = 0.01_dp*widtha/(wa/x(1))
                     IF (axis == 2) h = 0.01_dp*widtha/(wa/x(2))
                  ELSE
                     IF (axis == 1) h = 0.01_dp*widtha/(wa/x(1))
                     IF (axis == 2) h = 0.01_dp*widthb/(wb/x(2))
                     IF (axis == 3) h = 1.e-5_dp
                     IF (axis == 4) h = 0.01_dp*widtha/(wa/x(4))
                     IF (axis == 5) h = 0.01_dp*widthb/(wb/x(5))
                     IF (axis == 7) h = 0.01_dp*widthb
                  END IF
                  trial(:) = x(:); trial(axis) = x(axis) + h
                  CALL evaluate_libxc_point(names(family), trial, ep, gp, hp, contractions=.TRUE.)
                  trial(axis) = x(axis) - h
                  CALL evaluate_libxc_point(names(family), trial, em, gm, hm, contractions=.TRUE.)
                  trial(axis) = x(axis) + 2.0_dp*h
                  CALL evaluate_libxc_point(names(family), trial, epp, gpp, hpp, contractions=.TRUE.)
                  trial(axis) = x(axis) - 2.0_dp*h
                  CALL evaluate_libxc_point(names(family), trial, emm, gmm, hmm, contractions=.TRUE.)
                  first_fd = (8.0_dp*(ep - em) - (epp - emm))/(12.0_dp*h)
                  first_error = MAX(first_error, ABS(first_fd - gradient(axis))/(1.0_dp + ABS(gradient(axis))))
                  second_error = MAX(second_error, &
                      MAXVAL(ABS((8.0_dp*(gp - gm) - (gpp - gmm))/(12.0_dp*h) - hessian(:, axis))/(1.0_dp + ABS(hessian(:, axis)))))
                  IF (family <= 2) THEN
                     third_error = MAX(third_error, &
                                       MAXVAL(ABS((8.0_dp*(hp - hm) - (hpp - hmm))/(12.0_dp*h) - third(:, :, axis))/ &
                                              (1.0_dp + ABS(third(:, :, axis)))))
                  END IF
               END DO
            END DO
         END DO
         DEALLOCATE (x, trial, gradient, gm, gp, gmm, gpp, hessian, third, hm, hp, hmm, hpp)
      END DO
      WRITE (*, '(A,3ES15.6)') 'Libxc C2-window first/second/third derivative errors:', first_error, second_error, third_error
      IF (first_error > 2.e-6_dp) ERROR STOP 'C2-window energy Jacobian failed'
      IF (second_error > 2.e-5_dp) ERROR STOP 'C2-window potential Hessian failed'
      IF (third_error > 2.e-5_dp) ERROR STOP 'C2-window kernel third derivatives failed'
   END SUBROUTINE check_libxc_c2_window

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE check_libxc_c2_rounding()
      CHARACTER(len=32), PARAMETER :: names(4) = [CHARACTER(len=32):: 'MGGA_X_TPSS', 'MGGA_C_TPSS',&
         'MGGA_X_R2SCAN', 'MGGA_C_R2SCAN']
      REAL(dp), PARAMETER :: relative_width = 64.0_dp*SQRT(EPSILON(0.0_dp))

      INTEGER                                            :: channel, family, nvar, spin_case, &
                                                            tau_axis
      REAL(dp)                                           :: curve, energy, error, expected, slope, &
                                                            u, w(2), width(2)
      REAL(dp), ALLOCATABLE                              :: gp(:), gradient(:), hessian(:, :), &
                                                            hp(:, :), probe(:), x(:)

      error = 0.0_dp
      DO spin_case = 1, 2
         nvar = 3
         IF (spin_case == 2) nvar = 7
         ALLOCATE (x(nvar), probe(nvar), gradient(nvar), gp(nvar), hessian(nvar, nvar), hp(nvar, nvar))
         DO family = 1, SIZE(names)
            IF (spin_case == 1) THEN
               x(:) = [0.4_dp, 0.8_dp, 0.3_dp]
               w(1) = x(2)**2/(8.0_dp*x(1))
            ELSE
               x(:) = [0.3_dp, 0.2_dp, 0.6_dp, 0.4_dp, 0.3_dp, 0.12_dp, 0.10_dp]
               w(1) = x(4)**2/(8.0_dp*x(1))
               w(2) = x(5)**2/(8.0_dp*x(2))
            END IF
            probe(:) = x(:)
            DO channel = 1, spin_case
               tau_axis = nvar - spin_case + channel
               width(channel) = relative_width*w(channel)
               x(tau_axis) = w(channel) + 0.99999_dp*width(channel)
               ! The nearby unprojected point supplies raw Libxc derivatives.
               probe(tau_axis) = w(channel) + 1.00001_dp*width(channel)
            END DO
            CALL evaluate_libxc_point(names(family), x, energy, gradient, hessian)
            CALL evaluate_libxc_point(names(family), probe, energy, gp, hp)
            DO channel = 1, spin_case
               tau_axis = nvar - spin_case + channel
               u = (x(tau_axis) - w(channel))/width(channel)
               slope = (1.0_dp + u)**2*(2.0_dp - u)/4.0_dp
               curve = 0.75_dp*(1.0_dp - u*u)/width(channel)
               expected = hp(tau_axis, tau_axis)*slope*slope + gp(tau_axis)*curve
               error = MAX(error, ABS(hessian(tau_axis, tau_axis) - expected)/(1.0_dp + ABS(expected)))
            END DO
         END DO
         DEALLOCATE (x, probe, gradient, gp, hessian, hp)
      END DO
      WRITE (*, '(A,ES15.6)') 'Libxc rounded-window curvature error:', error
      IF (error > 2.e-6_dp) ERROR STOP 'Rounded energy value discarded C2 projection curvature'
   END SUBROUTINE check_libxc_c2_rounding

! **************************************************************************************************
!> \brief Check the original-density energy normalization below Libxc's individual spin floors.
! **************************************************************************************************
   SUBROUTINE check_libxc_individual_floor()
      CHARACTER(len=24), PARAMETER :: names(6) = [CHARACTER(len=24) :: 'LDA_C_PW', 'GGA_C_PBE', &
         'MGGA_X_TPSS', 'MGGA_C_TPSS', 'MGGA_X_R2SCAN', 'MGGA_C_R2SCAN']

      INTEGER                                            :: axis, configuration, f
      REAL(dp) :: energy, error_first, error_second, error_third, gm(7), gp(7), gradient(7), &
         hessian(7, 7), hm(7, 7), hp(7, 7), minus, plus, step, third(2, 2, 2), trial(7), x(7)

      error_first = 0.0_dp
      error_second = 0.0_dp
      error_third = 0.0_dp
      DO f = 1, SIZE(names)
         DO configuration = 1, 6
            x = [0.6_dp, 0.01_dp, 0.21_dp, 0.2_dp, 0.02_dp, 0.3_dp, 0.2_dp]
            SELECT CASE (configuration)
            CASE (2)
               x(2) = -0.01_dp
            CASE (3)
               x(1:2) = 0.03_dp
            CASE (4)
               x(1:2) = [0.01_dp, 0.6_dp]
            CASE (5)
               x(7) = -0.02_dp ! Compose the density floor with the kinetic-density projection.
            CASE (6)
               x(1:2) = [-0.01_dp, 0.6_dp]
            END SELECT
            IF (f <= 2) THEN
               CALL evaluate_libxc_point(TRIM(names(f)), x, energy, gradient, hessian, 0.05_dp, third)
            ELSE
               CALL evaluate_libxc_point(TRIM(names(f)), x, energy, gradient, hessian, 0.05_dp)
            END IF
            DO axis = 1, SIZE(x)
               step = 1.e-6_dp
               trial = x
               trial(axis) = x(axis) + step
               CALL evaluate_libxc_point(TRIM(names(f)), trial, plus, gp, hp, 0.05_dp)
               trial(axis) = x(axis) - step
               CALL evaluate_libxc_point(TRIM(names(f)), trial, minus, gm, hm, 0.05_dp)
               error_first = MAX(error_first, ABS((plus - minus)/(2.0_dp*step) - gradient(axis)))
               error_second = MAX(error_second, MAXVAL(ABS((gp - gm)/(2.0_dp*step) - hessian(:, axis))))
               IF (f <= 2 .AND. axis <= 2) THEN
                  error_third = MAX(error_third, MAXVAL(ABS((hp(1:2, 1:2) - hm(1:2, 1:2))/ &
                                                            (2.0_dp*step) - third(:, :, axis))))
               END IF
            END DO
         END DO
      END DO
      WRITE (*, '(A,3ES15.6)') 'Libxc individual spin-floor first/second/third derivative errors:', &
         error_first, error_second, error_third
      IF (error_first > 2.e-7_dp) ERROR STOP 'Libxc individual spin-floor potential is not an energy derivative'
      IF (error_second > 2.e-6_dp) ERROR STOP 'Libxc individual spin-floor kernel is not a potential derivative'
      IF (error_third > 2.e-5_dp) ERROR STOP 'Libxc individual spin-floor third derivatives are inconsistent'
   END SUBROUTINE check_libxc_individual_floor

! **************************************************************************************************
!> \brief Compare both floor limits with independent one-sided energy differences.
! **************************************************************************************************
   SUBROUTINE check_libxc_floor_limits()
      CHARACTER(len=24), PARAMETER :: &
         names(3) = [CHARACTER(len=24) :: 'MGGA_X_R2SCAN', 'MGGA_C_R2SCAN', 'GGA_C_PBE']
      REAL(dp), PARAMETER                                :: cutoff = 0.05_dp, step = 1.e-6_dp, &
                                                            weight = 0.37_dp

      INTEGER                                            :: family, side
      LOGICAL                                            :: spin_density_floor
      REAL(dp)                                           :: energy, error, first, gradient(7), &
                                                            hessian(7, 7), limits(7, 2), selected, &
                                                            trial(7), values(2), x(7)
      TYPE(section_vals_type), POINTER                   :: combination
      TYPE(xc_rho_cflags_type)                           :: needs

      NULLIFY (combination)
      CALL section_vals_create(combination, specification)
      CALL section_vals_val_set(combination, 'PBE%_SECTION_PARAMETERS_', l_val=.TRUE.)
      needs = xc_functionals_get_needs(combination, .TRUE., spin_density_floor=spin_density_floor)
      IF (spin_density_floor) ERROR STOP 'Native PBE incorrectly reports a Libxc spin floor'
      CALL check_spin_floor_targets(combination, .TRUE., [pbe_spin_density_floor])
      CALL section_vals_val_set(combination, 'GGA_C_PBE%_SECTION_PARAMETERS_', l_val=.TRUE.)
      needs = xc_functionals_get_needs(combination, .TRUE., spin_density_floor=spin_density_floor)
      IF (.NOT. spin_density_floor) ERROR STOP 'Active Libxc component lost its spin floor'
      CALL check_spin_floor_targets(combination, .TRUE., [pbe_spin_density_floor, libxc_density_floor(cutoff)])
      needs = xc_functionals_get_needs(combination, .FALSE., spin_density_floor=spin_density_floor)
      IF (spin_density_floor) ERROR STOP 'Unpolarized Libxc incorrectly reports an individual-spin floor'
      CALL check_spin_floor_targets(combination, .FALSE., [REAL(dp)::])
      CALL section_vals_val_set(combination, 'GGA_C_PBE%SCALE', r_val=0.0_dp)
      needs = xc_functionals_get_needs(combination, .TRUE., spin_density_floor=spin_density_floor)
      IF (spin_density_floor) ERROR STOP 'Zero-scaled Libxc component reports an active floor'
      CALL check_spin_floor_targets(combination, .TRUE., [pbe_spin_density_floor])
      CALL section_vals_val_set(combination, 'GGA_C_PBE%SCALE', r_val=1.0_dp)
      CALL section_vals_val_set(combination, 'GGA_C_PBE%_SECTION_PARAMETERS_', l_val=.FALSE.)
      needs = xc_functionals_get_needs(combination, .TRUE., spin_density_floor=spin_density_floor)
      IF (spin_density_floor) ERROR STOP 'Inactive Libxc component reports an active floor'
      CALL check_spin_floor_targets(combination, .TRUE., [pbe_spin_density_floor])
      CALL section_vals_val_set(combination, 'GGA_C_PBE%_SECTION_PARAMETERS_', l_val=.TRUE.)
      CALL section_vals_val_set(combination, 'PBE%_SECTION_PARAMETERS_', l_val=.FALSE.)
      CALL check_spin_floor_targets(combination, .TRUE., [libxc_density_floor(cutoff)])
      CALL section_vals_release(combination)
      error = 0.0_dp
      x = [0.4_dp, libxc_density_floor(cutoff), 0.7_dp, 0.6_dp, 0.2_dp, 0.3_dp, -0.02_dp]
      DO family = 1, SIZE(names)
         CALL evaluate_libxc_point(TRIM(names(family)), x, energy, gradient, hessian, cutoff)
         DO side = 1, 2
            CALL evaluate_libxc_point(TRIM(names(family)), x, selected, limits(:, side), hessian, cutoff, &
                                      floor_slopes=[-1.0_dp, REAL(side - 1, dp)])
            IF (ABS(selected - energy) > 1.e-14_dp) ERROR STOP 'A floor derivative selection changed the energy'
            trial = x
            trial(2) = x(2) + REAL(2*side - 3, dp)*step
            CALL evaluate_libxc_point(TRIM(names(family)), trial, values(1), gradient, hessian, cutoff)
            trial(2) = x(2) + REAL(2*side - 3, dp)*step*0.5_dp
            CALL evaluate_libxc_point(TRIM(names(family)), trial, values(2), gradient, hessian, cutoff)
            first = REAL(2*side - 3, dp)*(4.0_dp*values(2) - values(1) - 3.0_dp*energy)/step
            error = MAX(error, ABS(first - limits(2, side)))
            IF (ABS(first - limits(2, side)) > 1.e-7_dp) THEN
               WRITE (*, '(A,A,I3,3ES22.13)') 'Libxc floor limit ', TRIM(names(family)), side, &
                  first, limits(2, side), ABS(first - limits(2, side))
            END IF
         END DO
         CALL evaluate_libxc_point(TRIM(names(family)), x, selected, gradient, hessian, cutoff, &
                                   floor_slopes=[-1.0_dp, weight])
         IF (MAXVAL(ABS(gradient - ((1.0_dp - weight)*limits(:, 1) + weight*limits(:, 2)))) > 1.e-12_dp) THEN
            ERROR STOP 'Spin and kinetic density floors used different limiting weights'
         END IF
      END DO
      IF (error > 1.e-7_dp) ERROR STOP 'A selected spin-floor limit is not a one-sided energy derivative'
   END SUBROUTINE check_libxc_floor_limits

! **************************************************************************************************
!> \brief Check correlation potentials where a small spin power multiplies a large reduced tau.
! **************************************************************************************************
   SUBROUTINE check_libxc_spin_powers()
      CHARACTER(len=24), PARAMETER :: &
         names(2) = [CHARACTER(len=24) :: 'MGGA_C_TPSS', 'MGGA_C_R2SCAN']
      REAL(dp), PARAMETER :: &
         x(7) = [0.003_dp, -1.e-9_dp, 0.004_dp, 0.0040003_dp, 4.7e-7_dp, 0.000786_dp, 0.000112_dp]

      INTEGER                                            :: family, order
      REAL(dp)                                           :: energy, error, estimates(2), &
                                                            gradient(7), hessian(7, 7), minus, &
                                                            plus, scratch(7, 7), step, trial(7), &
                                                            values(7)

      error = 0.0_dp
      DO family = 1, SIZE(names)
         CALL evaluate_libxc_point(TRIM(names(family)), x, energy, gradient, hessian, 1.e-10_dp)
         DO order = 1, 2
            step = x(1)*1.e-3_dp/REAL(order, dp)
            trial = x
            trial(1) = x(1) + step
            CALL evaluate_libxc_point(TRIM(names(family)), trial, plus, values, scratch, 1.e-10_dp)
            trial(1) = x(1) - step
            CALL evaluate_libxc_point(TRIM(names(family)), trial, minus, values, scratch, 1.e-10_dp)
            estimates(order) = (plus - minus)/(2.0_dp*step)
         END DO
         error = MAX(error, ABS((4.0_dp*estimates(2) - estimates(1))/3.0_dp - gradient(1)))
      END DO
      WRITE (*, '(A,ES15.6)') 'Libxc small-spin-power potential error:', error
      IF (error > 2.e-7_dp) ERROR STOP 'Libxc small-spin-power potential is not an energy derivative'
   END SUBROUTINE check_libxc_spin_powers

! **************************************************************************************************
!> \brief Evaluate Libxc with independently varied density, gradient-norm and kinetic-density inputs.
!> \param name Functional name.
!> \param x Three restricted variables or seven spin variables.
!> \param energy Energy density.
!> \param gradient First derivatives with a preceding functional contribution removed.
!> \param hessian Second derivatives in the same independent variables.
!> \param density_cutoff Optional cutoff for the individual-spin-floor oracle.
!> \param third_rho Optional third derivatives in the two spin densities.
!> \param third_all Optional third derivatives in every independent variable.
!> \param contractions Use Libxc gradient contractions instead of CP2K gradient norms.
!> \param floor_slopes Optional first-derivative weights at an exact spin-density floor.
! **************************************************************************************************
   SUBROUTINE evaluate_libxc_point(name, x, energy, gradient, hessian, density_cutoff, third_rho, third_all, contractions, &
                                   floor_slopes)
      CHARACTER(len=*), INTENT(IN)                       :: name
      REAL(dp), INTENT(IN)                               :: x(:)
      REAL(dp), INTENT(OUT)                              :: energy, gradient(SIZE(x)), &
                                                            hessian(SIZE(x), SIZE(x))
      REAL(dp), INTENT(IN), OPTIONAL                     :: density_cutoff
      REAL(dp), INTENT(OUT), OPTIONAL :: third_rho(2, 2, 2), third_all(SIZE(x), SIZE(x), SIZE(x))
      LOGICAL, INTENT(IN), OPTIONAL                      :: contractions
      REAL(dp), INTENT(IN), OPTIONAL                     :: floor_slopes(2)

      INTEGER                                            :: bounds(2, 3), i, j, k, labels(7), &
                                                            nspins, order
      LOGICAL                                            :: use_contractions
      REAL(dp)                                           :: cutoff
      REAL(dp), POINTER                                  :: values(:, :, :)
      TYPE(section_vals_type), POINTER                   :: parameters
      TYPE(xc_derivative_set_type)                       :: derivatives
      TYPE(xc_derivative_type), POINTER                  :: derivative
      TYPE(xc_rho_cflags_type)                           :: needs
      TYPE(xc_rho_set_type)                              :: rho_set

      bounds = 1
      use_contractions = .FALSE.
      IF (PRESENT(contractions)) use_contractions = contractions
      labels = 0
      needs%norm_drho = .TRUE.
      IF (SIZE(x) == 3) THEN
         nspins = 1
         labels(1:3) = [deriv_rho, deriv_norm_drho, deriv_tau]
         needs%rho = .TRUE.
         needs%tau = .TRUE.
         needs%laplace_rho = .TRUE.
      ELSE
         nspins = 2
         labels = [deriv_rhoa, deriv_rhob, deriv_norm_drho, deriv_norm_drhoa, deriv_norm_drhob, deriv_tau_a, deriv_tau_b]
         needs%rho_spin = .TRUE.
         needs%norm_drho_spin = .TRUE.
         needs%tau_spin = .TRUE.
         needs%laplace_rho_spin = .TRUE.
      END IF
      IF (use_contractions) THEN
         IF (nspins == 1) THEN
            labels(2) = deriv_gamma
         ELSE
            labels(3:5) = [deriv_gamma_ab, deriv_gamma_aa, deriv_gamma_bb]
         END IF
      END IF
      cutoff = 1.e-20_dp
      IF (PRESENT(density_cutoff)) cutoff = density_cutoff
      CALL xc_rho_set_create(rho_set, bounds, rho_cutoff=cutoff, tau_cutoff=1.e-20_dp)
      CALL xc_rho_set_atom_update(rho_set, needs, nspins, bounds)
      IF (nspins == 1) THEN
         rho_set%rho = x(1)
         rho_set%norm_drho = x(2)
         rho_set%tau = x(3)
         rho_set%laplace_rho = 0.0_dp
         IF (use_contractions) rho_set%norm_drho = SQRT(x(2))
      ELSE
         rho_set%rhoa = x(1)
         rho_set%rhob = x(2)
         rho_set%norm_drho = x(3)
         rho_set%norm_drhoa = x(4)
         rho_set%norm_drhob = x(5)
         rho_set%tau_a = x(6)
         rho_set%tau_b = x(7)
         rho_set%laplace_rhoa = 0.0_dp
         rho_set%laplace_rhob = 0.0_dp
         IF (use_contractions) THEN
            rho_set%norm_drho = SQRT(2.0_dp*x(3) + x(4) + x(5))
            rho_set%norm_drhoa = SQRT(x(4))
            rho_set%norm_drhob = SQRT(x(5))
         END IF
      END IF
      rho_set%has = needs
      IF (PRESENT(floor_slopes)) THEN
         IF (nspins /= 2 .OR. PRESENT(third_rho) .OR. PRESENT(third_all)) ERROR STOP 'Invalid first-order floor control'
         ALLOCATE (rho_set%spin_floor_slopes(2, 1))
         rho_set%spin_floor_slopes(:, 1) = floor_slopes
      END IF
      CALL xc_dset_create(derivatives, local_bounds=bounds)
      DO i = 1, SIZE(x)
         derivative => xc_dset_get_derivative(derivatives, [labels(i)], allocate_deriv=.TRUE.)
         CALL xc_derivative_get(derivative, deriv_data=values)
         values = 0.125_dp
         DO j = i, SIZE(x)
            derivative => xc_dset_get_derivative(derivatives, [labels(i), labels(j)], allocate_deriv=.TRUE.)
            CALL xc_derivative_get(derivative, deriv_data=values)
            values = 0.125_dp
         END DO
      END DO
      order = 2
      IF (PRESENT(floor_slopes)) order = 1
      IF (PRESENT(third_rho) .OR. PRESENT(third_all)) THEN
         IF (PRESENT(third_rho) .AND. nspins /= 2) ERROR STOP 'Spin third-derivative oracle requires two spins'
         order = 3
         DO i = 1, SIZE(x)
            DO j = 1, SIZE(x)
               DO k = 1, SIZE(x)
                  derivative => xc_dset_get_derivative(derivatives, [labels(i), labels(j), labels(k)], allocate_deriv=.TRUE.)
                  CALL xc_derivative_get(derivative, deriv_data=values)
                  values = 0.125_dp
               END DO
            END DO
         END DO
      END IF
      NULLIFY (parameters)
      IF (nspins == 1) THEN
         CALL libxc_spin_unpolarized_eval(rho_set, derivatives, order, parameters, func_name_override=name)
      ELSE
         CALL libxc_spin_polarized_eval(rho_set, derivatives, order, parameters, func_name_override=name)
      END IF
      derivative => xc_dset_get_derivative(derivatives, [INTEGER::])
      CALL xc_derivative_get(derivative, deriv_data=values)
      energy = values(1, 1, 1)
      DO i = 1, SIZE(x)
         derivative => xc_dset_get_derivative(derivatives, [labels(i)])
         CALL xc_derivative_get(derivative, deriv_data=values)
         gradient(i) = values(1, 1, 1) - 0.125_dp
         DO j = i, SIZE(x)
            derivative => xc_dset_get_derivative(derivatives, [labels(i), labels(j)])
            CALL xc_derivative_get(derivative, deriv_data=values)
            hessian(i, j) = values(1, 1, 1) - 0.125_dp
            hessian(j, i) = hessian(i, j)
         END DO
      END DO
      IF (PRESENT(third_rho)) THEN
         DO i = 1, 2
            DO j = 1, 2
               DO k = 1, 2
                  derivative => xc_dset_get_derivative(derivatives, [labels(i), labels(j), labels(k)])
                  CALL xc_derivative_get(derivative, deriv_data=values)
                  third_rho(i, j, k) = values(1, 1, 1) - 0.125_dp
               END DO
            END DO
         END DO
      END IF
      IF (PRESENT(third_all)) THEN
         DO k = 1, SIZE(x)
            DO j = 1, SIZE(x)
               DO i = 1, SIZE(x)
                  derivative => xc_dset_get_derivative(derivatives, [labels(i), labels(j), labels(k)])
                  CALL xc_derivative_get(derivative, deriv_data=values)
                  third_all(i, j, k) = values(1, 1, 1) - 0.125_dp
               END DO
            END DO
         END DO
         IF (.NOT. ALL(ieee_is_finite(third_all))) ERROR STOP 'Nonfinite Libxc third derivatives'
      END IF
      IF (.NOT. ieee_is_finite(energy) .OR. .NOT. ALL(ieee_is_finite(gradient)) .OR. &
          .NOT. ALL(ieee_is_finite(hessian))) ERROR STOP 'Nonfinite Libxc projected-tau result'
      CALL xc_dset_release(derivatives)
      CALL xc_rho_set_release(rho_set)
   END SUBROUTINE evaluate_libxc_point
#endif

! **************************************************************************************************
!> \brief Check that actual physical floor values retain all active components exactly once.
!> \param combination Actual functional combination.
!> \param lsd Whether the physical density is spin polarized.
!> \param expected Distinct expected physical target values.
! **************************************************************************************************
   SUBROUTINE check_spin_floor_targets(combination, lsd, expected)
      TYPE(section_vals_type), POINTER                   :: combination
      LOGICAL, INTENT(IN)                                :: lsd
      REAL(dp), INTENT(IN)                               :: expected(:)

      INTEGER                                            :: i
      REAL(dp), ALLOCATABLE                              :: actual(:)

      ALLOCATE (actual, SOURCE=xc_functionals_spin_floors(combination, lsd, 0.05_dp))
      IF (SIZE(actual) /= SIZE(expected)) ERROR STOP 'Wrong number of component-dependent spin-density floors'
      DO i = 1, SIZE(expected)
         IF (COUNT(actual == expected(i)) /= 1) ERROR STOP 'Lost or duplicated an actual physical spin-density floor'
      END DO
   END SUBROUTINE check_spin_floor_targets

! **************************************************************************************************
!> \brief Check one-sided native PBE energy derivatives without moving its original density floor.
! **************************************************************************************************
   SUBROUTINE check_pbe_floor_limits()
      REAL(dp), PARAMETER                                :: level = pbe_spin_density_floor, &
                                                            step = level*1.e-3_dp, weight = 0.37_dp

      INTEGER                                            :: side, spin, variant
      LOGICAL                                            :: active
      REAL(dp)                                           :: energy, error, first, gradient(5), &
                                                            hessian(5, 5), limits(5, 2), selected, &
                                                            slopes(2), trial(5), values(2), x(5)
      TYPE(xc_rho_cflags_type)                           :: needs

      needs = xc_functionals_get_needs(functional, .TRUE., pbe_density_floor=active)
      IF (.NOT. active) ERROR STOP 'Native PBE lost its physical spin floor'
      CALL check_spin_floor_targets(functional, .TRUE., [pbe_spin_density_floor])
      needs = xc_functionals_get_needs(functional, .FALSE., pbe_density_floor=active)
      IF (active) ERROR STOP 'Unpolarized PBE reports an individual-spin floor'
      CALL check_spin_floor_targets(functional, .FALSE., [REAL(dp)::])
      CALL section_vals_val_set(functional, 'PBE%SCALE_X', r_val=0.0_dp)
      CALL section_vals_val_set(functional, 'PBE%SCALE_C', r_val=0.0_dp)
      needs = xc_functionals_get_needs(functional, .TRUE., pbe_density_floor=active)
      IF (active) ERROR STOP 'Zero-scaled PBE reports an active spin floor'
      CALL check_spin_floor_targets(functional, .TRUE., [REAL(dp)::])
      CALL section_vals_val_set(functional, 'PBE%SCALE_X', r_val=1.0_dp)
      CALL section_vals_val_set(functional, 'PBE%SCALE_C', r_val=1.0_dp)
      error = 0.0_dp
      DO variant = 1, SIZE(variants)
         CALL section_vals_val_set(functional, 'PBE%PARAMETRIZATION', i_val=variants(variant))
         DO spin = 1, 2
            x(1:2) = 2.0_dp*level
            x(spin) = level
            x(3:5) = [0.7_dp, 0.4_dp, 0.3_dp]*level**(4.0_dp/3.0_dp)
            CALL evaluate(x, energy, gradient, hessian, derivative_order=1)
            DO side = 1, 2
               slopes = -1.0_dp
               slopes(spin) = REAL(side - 1, dp)
               CALL evaluate(x, selected, limits(:, side), hessian, derivative_order=1, floor_slopes=slopes)
               IF (selected /= energy) ERROR STOP 'A native PBE floor selection changed the energy'
               trial = x
               trial(spin) = x(spin) + REAL(2*side - 3, dp)*step
               CALL evaluate(trial, values(1), gradient, hessian, derivative_order=1)
               trial(spin) = x(spin) + REAL(2*side - 3, dp)*step*0.5_dp
               CALL evaluate(trial, values(2), gradient, hessian, derivative_order=1)
               first = REAL(2*side - 3, dp)*(4.0_dp*values(2) - values(1) - 3.0_dp*energy)/step
               error = MAX(error, ABS(first - limits(spin, side)))
            END DO
            slopes(spin) = weight
            CALL evaluate(x, selected, gradient, hessian, derivative_order=1, floor_slopes=slopes)
            IF (MAXVAL(ABS(gradient - ((1.0_dp - weight)*limits(:, 1) + weight*limits(:, 2)))) > 1.e-12_dp) THEN
               ERROR STOP 'Native PBE did not preserve accumulated contributions and independent gradient derivatives'
            END IF
         END DO
      END DO
      WRITE (*, '(A,ES15.6)') 'Native PBE selected spin-floor energy derivative error:', error
      IF (error > 2.e-8_dp) ERROR STOP 'A native PBE floor limit is not a one-sided energy derivative'
   END SUBROUTINE check_pbe_floor_limits

! **************************************************************************************************
!> \brief Check component-specific floor routing, physical keys and clearing between atomic batches.
! **************************************************************************************************
   SUBROUTINE check_xc_floor_selection()
      INTEGER                                            :: bounds(2, 3), domain, families, family, &
                                                            keys(3, 4), point(3, 1)
      REAL(dp)                                           :: hmat(3, 3), level
      TYPE(section_vals_type), POINTER                   :: combination
      TYPE(skala_density_constraints_type)               :: constraints
      TYPE(xc_rho_cflags_type)                           :: needs
      TYPE(xc_rho_set_type)                              :: fields

      NULLIFY (combination)
      CALL section_vals_create(combination, specification)
      CALL section_vals_val_set(combination, 'PBE%_SECTION_PARAMETERS_', l_val=.TRUE.)
      families = 1
#if defined(__LIBXC)
      CALL section_vals_val_set(combination, 'GGA_C_PBE%_SECTION_PARAMETERS_', l_val=.TRUE.)
      families = 2
#endif
      bounds(:, 1) = [-1, 0]
      bounds(:, 2) = [2, 3]
      bounds(:, 3) = [4, 4]
      hmat = 0.0_dp
      hmat(1, 1) = 1.0_dp
      hmat(2, 2) = 1.0_dp
      hmat(3, 3) = 1.0_dp
      needs%rho_spin = .TRUE.
      CALL xc_rho_set_create(fields, bounds, rho_cutoff=0.05_dp)
      CALL xc_rho_set_atom_update(fields, needs, 2, bounds)
      fields%rhoa = 0.4_dp
      fields%rhob = 0.3_dp
      DO family = 1, families
         level = pbe_spin_density_floor
#if defined(__LIBXC)
         IF (family == 2) level = libxc_density_floor(fields%rho_cutoff)
#endif
         DO domain = skala_fft_domain, skala_last_sample_domain
            fields%rhob(0, 2, 4) = level
            point(:, 1) = [2, 4, 1]
            IF (domain == skala_fft_domain) point(:, 1) = [0, 2, 4]
            IF (domain > 2) point(:, 1) = [2, 4, 2]
            keys = 0
            keys(:, 2) = point(:, 1)
            CALL skala_density_constraints_set(constraints, point, [2], [0.37_dp], 1.e-14_dp, 1.e-20_dp, &
                                               2, [2, 2, 1], bounds, hmat, RESHAPE([REAL(dp)::0, 0, 0, 1, 0, 0], [3, 2]), &
                                               domains=[domain], atomic_sizes=[0, 6], targets=[level], &
                                               paw_sizes=RESHAPE([0, 0, 6, 3], [2, 2]))
            IF (domain == skala_fft_domain) THEN
               CALL xc_functionals_select_floor_limits(combination, constraints, fields, 2, domain)
            ELSE
               CALL xc_functionals_select_floor_limits(combination, constraints, fields, 2, domain, keys)
            END IF
            IF (family == 1) THEN
               IF (.NOT. ALLOCATED(fields%pbe_spin_floor_slopes)) ERROR STOP 'Missing native PBE floor derivative'
               IF (ALLOCATED(fields%spin_floor_slopes)) ERROR STOP 'PBE floor selected a Libxc derivative'
               IF (ABS(fields%pbe_spin_floor_slopes(2, 2) - 0.63_dp) > 1.e-15_dp) ERROR STOP 'Wrong PBE floor row'
               fields%pbe_spin_floor_slopes(2, 2) = -1.0_dp
               IF (ANY(fields%pbe_spin_floor_slopes /= -1.0_dp)) ERROR STOP 'PBE selected an unrelated sample'
            ELSE
               IF (.NOT. ALLOCATED(fields%spin_floor_slopes)) ERROR STOP 'Missing Libxc floor derivative'
               IF (ALLOCATED(fields%pbe_spin_floor_slopes)) ERROR STOP 'Libxc floor selected a PBE derivative'
               IF (ABS(fields%spin_floor_slopes(2, 2) - 0.63_dp) > 1.e-15_dp) ERROR STOP 'Wrong Libxc floor row'
               fields%spin_floor_slopes(2, 2) = -1.0_dp
               IF (ANY(fields%spin_floor_slopes /= -1.0_dp)) ERROR STOP 'Libxc selected an unrelated sample'
            END IF
            fields%rhob(0, 2, 4) = 2.0_dp*level
            CALL xc_functionals_select_floor_limits(combination, constraints, fields, 2, domain, keys)
            IF (ALLOCATED(fields%spin_floor_slopes) .OR. ALLOCATED(fields%pbe_spin_floor_slopes)) THEN
               ERROR STOP 'A stale floor derivative survived an off-surface batch'
            END IF
            CALL skala_density_constraints_reset(constraints)
         END DO
      END DO
      CALL xc_rho_set_release(fields)
      CALL section_vals_release(combination)
   END SUBROUTINE check_xc_floor_selection

! **************************************************************************************************
!> \brief Evaluate native PBE in its independent spin-density and gradient-norm variables.
!> \param x Two spin densities followed by total, alpha and beta gradient norms.
!> \param energy Energy density.
!> \param gradient First derivatives, with the accumulation of another functional preserved.
!> \param hessian Second derivatives in the same independent variables.
!> \param derivative_order ...
!> \param floor_slopes Optional first-derivative weights at the native PBE spin floor.
! **************************************************************************************************
   SUBROUTINE evaluate(x, energy, gradient, hessian, derivative_order, floor_slopes)
      REAL(dp), INTENT(IN)                               :: x(5)
      REAL(dp), INTENT(OUT)                              :: energy, gradient(5), hessian(5, 5)
      INTEGER, INTENT(IN), OPTIONAL                      :: derivative_order
      REAL(dp), INTENT(IN), OPTIONAL                     :: floor_slopes(2)

      INTEGER                                            :: bounds(2, 3), i, j, order
      REAL(dp), POINTER                                  :: values(:, :, :)
      TYPE(xc_derivative_set_type)                       :: derivatives
      TYPE(xc_derivative_type), POINTER                  :: derivative
      TYPE(xc_rho_cflags_type)                           :: needs
      TYPE(xc_rho_set_type)                              :: rho_set

      order = 2
      IF (PRESENT(derivative_order)) order = derivative_order
      bounds = 1
      needs = xc_functionals_get_needs(functional, lsd=.TRUE., calc_potential=.TRUE.)
      needs%norm_drho = .TRUE.
      CALL xc_rho_set_create(rho_set, bounds, rho_cutoff=1.e-20_dp)
      CALL xc_rho_set_atom_update(rho_set, needs, 2, bounds)
      rho_set%rhoa = x(1)
      rho_set%rhob = x(2)
      rho_set%norm_drho = x(3)
      rho_set%norm_drhoa = x(4)
      rho_set%norm_drhob = x(5)
      rho_set%has%rho_spin = .TRUE.
      rho_set%has%norm_drho = .TRUE.
      rho_set%has%norm_drho_spin = .TRUE.
      IF (PRESENT(floor_slopes)) THEN
         ALLOCATE (rho_set%pbe_spin_floor_slopes(2, 1))
         rho_set%pbe_spin_floor_slopes(:, 1) = floor_slopes
         rho_set%spin_floor_tolerance = 1.e-20_dp
      END IF
      CALL xc_dset_create(derivatives, local_bounds=bounds)
      ! Seed every first/second derivative to catch overwriting earlier functional contributions.
      DO i = 1, 5
         derivative => xc_dset_get_derivative(derivatives, [labels(i)], allocate_deriv=.TRUE.)
         CALL xc_derivative_get(derivative, deriv_data=values)
         values = 0.125_dp
         DO j = i, 5
            derivative => xc_dset_get_derivative(derivatives, [labels(i), labels(j)], allocate_deriv=.TRUE.)
            CALL xc_derivative_get(derivative, deriv_data=values)
            values = 0.125_dp
         END DO
      END DO
      CALL xc_functionals_eval(functional, lsd=.TRUE., rho_set=rho_set, deriv_set=derivatives, deriv_order=order)
      derivative => xc_dset_get_derivative(derivatives, [INTEGER::])
      CALL xc_derivative_get(derivative, deriv_data=values)
      energy = values(1, 1, 1)
      DO i = 1, 5
         derivative => xc_dset_get_derivative(derivatives, [labels(i)])
         CALL xc_derivative_get(derivative, deriv_data=values)
         gradient(i) = values(1, 1, 1) - 0.125_dp
         DO j = i, 5
            derivative => xc_dset_get_derivative(derivatives, [labels(i), labels(j)])
            CALL xc_derivative_get(derivative, deriv_data=values)
            hessian(i, j) = values(1, 1, 1) - 0.125_dp
            IF (order == 1) hessian(i, j) = 0.0_dp
            hessian(j, i) = hessian(i, j)
         END DO
      END DO
      CALL xc_dset_release(derivatives)
      CALL xc_rho_set_release(rho_set)
   END SUBROUTINE evaluate

! **************************************************************************************************
!> \brief Check native TPSS derivatives at stationary density points and below the von Weizsaecker bound.
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
            trial(:) = x(:)
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
      ! NLCC can give tau below |grad rho|^2/(8*rho). Check the projection Jacobian.
      x = [0.4_dp, 1.2_dp, 0.05_dp]
      CALL evaluate_tpss_point(x, energy, gradient)
      DO axis = 1, 3
         trial = x
         trial(axis) = x(axis) + step
         CALL evaluate_tpss_point(trial, plus, scratch)
         trial(axis) = x(axis) - step
         CALL evaluate_tpss_point(trial, minus, scratch)
         error = MAX(error, ABS((plus - minus)/(2.0_dp*step) - gradient(axis)))
      END DO
      ! Negative Fourier-interpolated tau must not discard a finite-density energy.
      DO configuration = 1, 2
         x = [0.4_dp, 1.2_dp, -0.05_dp]
         IF (configuration == 2) x(2) = 0.0_dp
         CALL evaluate_tpss_point(x, energy, gradient)
         IF (ABS(energy) < 1.e-4_dp .OR. gradient(3) /= 0.0_dp) THEN
            ERROR STOP 'TPSS negative-tau continuation is not finite and constant in tau'
         END IF
         DO axis = 1, 3
            IF (axis == 2 .AND. configuration == 2) CYCLE
            trial = x
            trial(axis) = x(axis) + step
            CALL evaluate_tpss_point(trial, plus, scratch)
            trial(axis) = x(axis) - step
            CALL evaluate_tpss_point(trial, minus, scratch)
            error = MAX(error, ABS((plus - minus)/(2.0_dp*step) - gradient(axis)))
         END DO
      END DO
      x = [0.4_dp, 1.2_dp, 0.0_dp]
      CALL evaluate_tpss_point(x, energy, gradient)
      trial = x
      trial(3) = 1.e-10_dp
      CALL evaluate_tpss_point(trial, plus, scratch)
      trial(3) = -1.e-10_dp
      CALL evaluate_tpss_point(trial, minus, scratch)
      IF (MAX(ABS(plus - energy), ABS(minus - energy)) > 1.e-9_dp) THEN
         ERROR STOP 'TPSS energy is discontinuous through zero kinetic density'
      END IF
      WRITE (*, '(A,ES15.6)') 'TPSS limiting-density first-derivative error:', error
      IF (error > 3.e-8_dp) ERROR STOP 'TPSS limiting-density potential is not the energy derivative'
   END SUBROUTINE check_tpss_zero_gradient

! **************************************************************************************************
!> \brief Check screened exchange limits and gradient scaling against energy differences.
! **************************************************************************************************
   SUBROUTINE check_xwpbe_limits()
      REAL(dp), PARAMETER :: reduced_gradients(3) = [0.0_dp, 0.2_dp, 12.0_dp], step = 1.e-5_dp

      INTEGER                                            :: axis, configuration, range_case
      REAL(dp)                                           :: energy, error, gr(3), gs(5), &
                                                            hessian(5, 5), minus, plus, &
                                                            scratchr(3), scratchs(5), trialr(3), &
                                                            trials(5), xr(3), xs(5)

      error = 0.0_dp
      CALL section_vals_val_set(functional, 'XWPBE%SCALE_X', r_val=-0.25_dp)
      CALL section_vals_val_set(functional, 'XWPBE%SCALE_X0', r_val=1.0_dp)
      DO range_case = 1, 2
         CALL section_vals_val_set(functional, 'XWPBE%OMEGA', r_val=0.11_dp*REAL(range_case - 1, dp))
         DO configuration = 1, SIZE(reduced_gradients)
            xr = [0.4_dp, 0.0_dp, 0.5_dp]
            xr(2) = reduced_gradients(configuration)*2.0_dp*(3.0_dp*pi**2)**(1.0_dp/3.0_dp)*xr(1)**(4.0_dp/3.0_dp)
            CALL evaluate_tpss_point(xr, energy, gr)
            DO axis = 1, 2
               IF (axis == 2 .AND. configuration == 1) CYCLE
               trialr = xr; trialr(axis) = xr(axis) + step
               CALL evaluate_tpss_point(trialr, plus, scratchr)
               trialr(axis) = xr(axis) - step
               CALL evaluate_tpss_point(trialr, minus, scratchr)
               error = MAX(error, ABS((plus - minus)/(2.0_dp*step) - gr(axis))/(1.0_dp + ABS(gr(axis))))
            END DO
            IF (configuration == 1) THEN
               IF (ABS(gr(2)) > 1.e-13_dp) ERROR STOP 'XWPBE zero-gradient potential is not stationary'
               trialr = xr; trialr(2) = 1.e-8_dp
               CALL evaluate_tpss_point(trialr, plus, scratchr)
               IF (ABS(plus - energy) > 1.e-12_dp) ERROR STOP 'XWPBE zero-gradient energy is discontinuous'
            END IF
            xs = [0.3_dp, 0.2_dp, 0.3_dp, 0.0_dp, 0.0_dp]
            xs(4) = reduced_gradients(configuration)*2.0_dp*(3.0_dp*pi**2)**(1.0_dp/3.0_dp)*(2.0_dp*xs(1))**(4.0_dp/3.0_dp)/2.0_dp
            xs(5) = reduced_gradients(configuration)*2.0_dp*(3.0_dp*pi**2)**(1.0_dp/3.0_dp)*(2.0_dp*xs(2))**(4.0_dp/3.0_dp)/2.0_dp
            CALL evaluate(xs, energy, gs, hessian, derivative_order=1)
            DO axis = 1, 5
               IF (axis >= 4 .AND. configuration == 1) CYCLE
               trials = xs; trials(axis) = xs(axis) + step
               CALL evaluate(trials, plus, scratchs, hessian, derivative_order=1)
               trials(axis) = xs(axis) - step
               CALL evaluate(trials, minus, scratchs, hessian, derivative_order=1)
               error = MAX(error, ABS((plus - minus)/(2.0_dp*step) - gs(axis))/(1.0_dp + ABS(gs(axis))))
            END DO
    IF (configuration == 1 .AND. MAXVAL(ABS(gs(4:5))) > 1.e-13_dp) ERROR STOP 'XWPBE spin zero-gradient potential is not stationary'
         END DO
      END DO
      WRITE (*, '(A,ES15.6)') 'XWPBE limiting/large-gradient first-derivative error:', error
      IF (error > 3.e-8_dp) ERROR STOP 'XWPBE potential is not the derivative of its scaled-gradient energy'
   END SUBROUTINE check_xwpbe_limits

! **************************************************************************************************
!> \brief Evaluate an unpolarized semilocal functional in density, gradient-norm and kinetic-density variables.
!> \param x ...
!> \param energy ...
!> \param gradient ...
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
      needs%tau = .TRUE.
      CALL xc_rho_set_create(rho_set, bounds, rho_cutoff=1.e-20_dp, tau_cutoff=1.e-20_dp)
      CALL xc_rho_set_atom_update(rho_set, needs, 1, bounds)
      rho_set%rho = x(1)
      rho_set%norm_drho = x(2)
      rho_set%tau = x(3)
      rho_set%has%rho = .TRUE.
      rho_set%has%norm_drho = .TRUE.
      rho_set%has%tau = .TRUE.
      CALL xc_dset_create(derivatives, local_bounds=bounds)
      DO i = 1, SIZE(labels)
         derivative => xc_dset_get_derivative(derivatives, [labels(i)], allocate_deriv=.TRUE.)
         CALL xc_derivative_get(derivative, deriv_data=values)
         values = 0.125_dp
      END DO
      CALL xc_functionals_eval(functional, lsd=.FALSE., rho_set=rho_set, deriv_set=derivatives, deriv_order=1)
      derivative => xc_dset_get_derivative(derivatives, [INTEGER::])
      CALL xc_derivative_get(derivative, deriv_data=values)
      energy = values(1, 1, 1)
      DO i = 1, SIZE(labels)
         derivative => xc_dset_get_derivative(derivatives, [labels(i)])
         CALL xc_derivative_get(derivative, deriv_data=values)
         gradient(i) = values(1, 1, 1) - 0.125_dp
      END DO
      IF (.NOT. ieee_is_finite(energy) .OR. .NOT. ALL(ieee_is_finite(gradient))) THEN
         ERROR STOP 'TPSS energy or first derivative is nonfinite at a stationary density point'
      END IF
      CALL xc_dset_release(derivatives)
      CALL xc_rho_set_release(rho_set)
   END SUBROUTINE evaluate_tpss_point

END PROGRAM xc_semilocal_unittest
