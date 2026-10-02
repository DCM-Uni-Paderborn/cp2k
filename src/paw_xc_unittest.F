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
   USE paw_onecenter,                   ONLY: paw_onecenter_init,&
                                              paw_onecenter_type
   USE paw_xc,                          ONLY: paw_xc_collocate,&
                                              paw_xc_evaluate,&
                                              paw_xc_init,&
                                              paw_xc_kinetic_adjoint,&
                                              paw_xc_type
   USE upf_projector_integrals,         ONLY: upf_projector_init,&
                                              upf_projector_type
   USE xc_atom,                         ONLY: xc_rho_set_atom_update
   USE xc_derivative_desc,              ONLY: deriv_norm_drho,&
                                              deriv_rho,&
                                              deriv_tau
   USE xc_derivative_set_types,         ONLY: xc_derivative_set_type,&
                                              xc_dset_create,&
                                              xc_dset_get_derivative,&
                                              xc_dset_release
   USE xc_derivative_types,             ONLY: xc_derivative_get,&
                                              xc_derivative_type
   USE xc_derivatives,                  ONLY: xc_functionals_eval,&
                                              xc_functionals_get_needs
   USE xc_input_constants,              ONLY: xc_pbe_orig,&
                                              xc_pbe_sol
   USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_type
   USE xc_rho_set_types,                ONLY: xc_rho_set_create,&
                                              xc_rho_set_release,&
                                              xc_rho_set_type

   IMPLICIT NONE
   INTEGER, PARAMETER :: n = 801
   REAL(dp), PARAMETER :: exponent = 0.9_dp, pseudo_scale = 0.8_dp**2 + 0.07_dp
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
   CALL check_atomic_projection()
   CALL original_datasets()
   CALL section_vals_release(functional)
   CALL section_release(specification)
   CALL cp2k_finalize()

CONTAINS

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
            CALL check_tpss_zero_gradient()
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
         IF (ABS(plus - energy) > 1.e-13_dp .OR. ABS(scratch(2)) > 1.e-20_dp) &
            ERROR STOP 'TPSS does not approach its zero-gradient limit continuously'
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
      IF (.NOT. ieee_is_finite(energy) .OR. .NOT. ALL(ieee_is_finite(gradient))) &
         ERROR STOP 'TPSS energy or first derivative is nonfinite at a stationary density point'
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
         DEALLOCATE (density, direction, potential, scratch)
         CALL atom_release_upf(original)
      END DO
   END SUBROUTINE original_datasets

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
