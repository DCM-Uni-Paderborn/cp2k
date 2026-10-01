!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM paw_xc_unittest
   USE atom_upf,                        ONLY: atom_read_upf,&
                                              atom_release_upf,&
                                              atom_upfpot_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE input_cp2k_xc,                   ONLY: create_xc_fun_section
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
                                              paw_xc_type
   USE upf_projector_integrals,         ONLY: upf_projector_init,&
                                              upf_projector_type
   USE xc_input_constants,              ONLY: xc_pbe_orig,&
                                              xc_pbe_sol

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
   CALL check_local_exchange()
   CALL check_gga()
   CALL original_datasets()
   CALL section_vals_release(functional)
   CALL section_release(specification)
   CALL cp2k_finalize()

CONTAINS

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
         rho(grid%na, 2), scale

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
      DO shell = 1, 3
         ir = 1 + 50*shell
         DO side = 1, 2
            scale = MERGE(1.0_dp, pseudo_scale, side == 1)
            core = MERGE(0.2_dp, 0.1_dp, side == 1)/2.0_dp
            CALL paw_xc_collocate(grid, onecenter, density, ir, side, rho, gradient)
            DO a = 1, grid%na
               point = grid%r(ir)*grid%directions(:, a)
               basis = [1.0_dp, SQRT(3.0_dp)*point(2), SQRT(3.0_dp)*point(3), SQRT(3.0_dp)*point(1)]
               factor = EXP(-exponent*SUM(point**2))
               DO spin = 1, 2
                  polynomial = DOT_PRODUCT(basis, MATMUL(density(:, :, spin), basis))
                  exact = factor*(scale*polynomial/fourpi + core)
                  DO axis = 1, 3
                     exact_gradient(axis) = scale*factor*2.0_dp* &
                                            DOT_PRODUCT(dbasis(:, axis), MATMUL(density(:, :, spin), basis))/fourpi &
                                            - 2.0_dp*exponent*point(axis)*exact
                  END DO
                  field_error = MAX(field_error, ABS(rho(a, spin) - exact))
                  gradient_error = MAX(gradient_error, MAXVAL(ABS(gradient(:, a, spin) - exact_gradient)))
               END DO
            END DO
         END DO
      END DO
      WRITE (*, '(A,2ES15.6)') 'PAW XC Cartesian density/gradient errors:', field_error, gradient_error
      IF (field_error > 2.e-13_dp) ERROR STOP 'PAW XC density differs from Cartesian Gaussian result'
      IF (gradient_error > 2.e-8_dp) ERROR STOP 'PAW XC gradient differs from Cartesian Gaussian result'
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
