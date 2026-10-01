!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM xc_pbe_unittest
   USE input_cp2k_xc,                   ONLY: create_xc_fun_section
   USE input_section_types,             ONLY: section_release,&
                                              section_type,&
                                              section_vals_create,&
                                              section_vals_release,&
                                              section_vals_type,&
                                              section_vals_val_set
   USE kinds,                           ONLY: dp
   USE libcp2k,                         ONLY: cp2k_finalize,&
                                              cp2k_init
   USE xc_atom,                         ONLY: xc_rho_set_atom_update
   USE xc_derivative_desc,              ONLY: deriv_norm_drho,&
                                              deriv_norm_drhoa,&
                                              deriv_norm_drhob,&
                                              deriv_rhoa,&
                                              deriv_rhob
   USE xc_derivative_set_types,         ONLY: xc_derivative_set_type,&
                                              xc_dset_create,&
                                              xc_dset_get_derivative,&
                                              xc_dset_release
   USE xc_derivative_types,             ONLY: xc_derivative_get,&
                                              xc_derivative_type
   USE xc_derivatives,                  ONLY: xc_functionals_eval,&
                                              xc_functionals_get_needs
   USE xc_input_constants,              ONLY: xc_pbe_orig,&
                                              xc_pbe_rev,&
                                              xc_pbe_sol
   USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_type
   USE xc_rho_set_types,                ONLY: xc_rho_set_create,&
                                              xc_rho_set_release,&
                                              xc_rho_set_type

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
   CALL section_vals_release(functional)
   CALL section_release(specification)
   CALL cp2k_finalize()

CONTAINS

! **************************************************************************************************
!> \brief Evaluate native PBE in its independent spin-density and gradient-norm variables.
!> \param x Two spin densities followed by total, alpha and beta gradient norms.
!> \param energy Energy density.
!> \param gradient First derivatives, with the accumulation of another functional preserved.
!> \param hessian Second derivatives in the same independent variables.
! **************************************************************************************************
   SUBROUTINE evaluate(x, energy, gradient, hessian)
      REAL(dp), INTENT(IN)                               :: x(5)
      REAL(dp), INTENT(OUT)                              :: energy, gradient(5), hessian(5, 5)

      INTEGER                                            :: bounds(2, 3), i, j
      REAL(dp), POINTER                                  :: values(:, :, :)
      TYPE(xc_derivative_set_type)                       :: derivatives
      TYPE(xc_derivative_type), POINTER                  :: derivative
      TYPE(xc_rho_cflags_type)                           :: needs
      TYPE(xc_rho_set_type)                              :: rho_set

      bounds = 1
      needs = xc_functionals_get_needs(functional, lsd=.TRUE., calc_potential=.TRUE.)
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
      CALL xc_functionals_eval(functional, lsd=.TRUE., rho_set=rho_set, deriv_set=derivatives, deriv_order=2)
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
            hessian(j, i) = hessian(i, j)
         END DO
      END DO
      CALL xc_dset_release(derivatives)
      CALL xc_rho_set_release(rho_set)
   END SUBROUTINE evaluate

END PROGRAM xc_pbe_unittest
