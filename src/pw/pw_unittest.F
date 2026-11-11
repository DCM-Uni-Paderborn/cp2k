!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM pw_unittest
   USE bessel_lib,                      ONLY: bessk1
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE fft_tools,                       ONLY: finalize_fft,&
                                              init_fft
   USE input_cp2k_xc,                   ONLY: create_xc_section
   USE input_section_types,             ONLY: section_release,&
                                              section_type,&
                                              section_vals_create,&
                                              section_vals_release,&
                                              section_vals_type,&
                                              section_vals_val_set
   USE kinds,                           ONLY: dp
   USE libcp2k,                         ONLY: cp2k_finalize,&
                                              cp2k_init
   USE mathconstants,                   ONLY: fourpi,&
                                              pi
   USE mathlib,                         ONLY: invmat
   USE message_passing,                 ONLY: mp_comm_self,&
                                              mp_comm_type,&
                                              mp_comm_world,&
                                              mp_para_env_type
   USE mt_util,                         ONLY: mt2d_kernel_strain
   USE pw_copy_all,                     ONLY: pw_copy_match
   USE pw_grid_types,                   ONLY: FULLSPACE,&
                                              HALFSPACE,&
                                              pw_grid_type
   USE pw_grids,                        ONLY: do_pw_grid_blocked_false,&
                                              pw_grid_create,&
                                              pw_grid_release
   USE pw_methods,                      ONLY: do_accurate_sum,&
                                              do_standard_sum,&
                                              pw_copy,&
                                              pw_integral_ab,&
                                              pw_real_fourier_aliases,&
                                              pw_transfer
   USE pw_poisson_methods,              ONLY: pw_green_strain_response,&
                                              pw_poisson_set,&
                                              pw_poisson_solve
   USE pw_poisson_types,                ONLY: greens_fn_type,&
                                              pw_green_create,&
                                              pw_green_release,&
                                              pw_poisson_analytic,&
                                              pw_poisson_mt,&
                                              pw_poisson_parameter_type,&
                                              pw_poisson_type
   USE pw_pool_types,                   ONLY: pw_pool_create,&
                                              pw_pool_release,&
                                              pw_pool_type
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r1d_gs_type,&
                                              pw_r3d_rs_type
   USE qs_external_potential,           ONLY: interpolate_external_potential
   USE qs_upf_inverse_scf_xc,           ONLY: qs_upf_inverse_xc_transfer_real
   USE qs_upf_inverse_xc,               ONLY: qs_upf_inverse_xc_grid_kernel,&
                                              qs_upf_inverse_xc_pw_adjoint,&
                                              qs_upf_inverse_xc_pw_fields
   USE realspace_grid_cube,             ONLY: cube_read_values,&
                                              cube_value_format,&
                                              cube_values_format
   USE rs_pw_interface,                 ONLY: pw_real_transfer
   USE skala_gpw_functional,            ONLY: build_vxc_from_feature_grads
   USE xc,                              ONLY: xc_vxc_pw_create
   USE xc_input_constants,              ONLY: xc_deriv_pw,&
                                              xc_rho_nn10,&
                                              xc_rho_nn50,&
                                              xc_rho_spline2_smooth,&
                                              xc_rho_spline3_smooth
   USE xc_util,                         ONLY: xc_pw_smooth

   IMPLICIT NONE

   INTEGER                                            :: mode
   REAL(dp)                                           :: max_error
   TYPE(cp_logger_type), POINTER                      :: logger
   TYPE(mp_para_env_type), POINTER                     :: para_env

   CALL cp2k_init()
   logger => cp_get_default_logger()
   para_env => logger%para_env
   CALL init_fft("FFTSG", alltoall=.FALSE., pool_limit=4, wisdom_file="", plan_style=1)
   CALL check_cube_values()
   DO mode = 0, 8
      CALL check_copy(mode)
   END DO
   CALL check_real_interpolation(.FALSE.)
   CALL check_real_interpolation(.TRUE.)
   max_error = 0.0_dp
   CALL check_upf_xc_smoothing()
   CALL check_pw_integral_and_poisson()
   CALL check_mt2d_kernel()
   CALL check_green_strain()
   PRINT '(A,ES16.7)', 'PW integral and Coulomb kernel maximum scaled error: ', max_error
   CALL finalize_fft(para_env, wisdom_file="")
   CALL cp2k_finalize()

CONTAINS

! **************************************************************************************************!
!> \brief Use the native XC filter in both directions of the mixed UPF grid map.
! **************************************************************************************************!
   SUBROUTINE check_upf_xc_smoothing()
      INTEGER                                            :: axis, i, j, k, mode, modes(4), nrows, &
                                                            point(3)
      REAL(dp)                                           :: cell(3, 3), phase, projected, sampled
      REAL(dp), ALLOCATABLE                              :: fields(:, :, :, :), &
                                                            reference(:, :, :, :), vgrad(:, :, :), &
                                                            vrho(:, :), vtau(:, :)
      TYPE(mp_comm_type)                                 :: comm
      TYPE(pw_grid_type), POINTER                        :: fine_grid, grid
      TYPE(pw_pool_type), POINTER                        :: pool
      TYPE(pw_r3d_rs_type)                               :: fine_density, fine_filtered, fine_point, &
                                                            kinetic(1), potential(1), smooth(1), &
                                                            tau_potential(1)
      TYPE(pw_r3d_rs_type), POINTER                      :: density_ptr(:), skala_density(:), &
                                                            skala_tau(:)
      TYPE(pw_r3d_rs_type), TARGET                       :: density(1)

      NULLIFY (fine_grid, grid, pool, skala_density, skala_tau)
      density_ptr => density
      CALL comm%set_handle(mp_comm_world%get_handle())
      cell = 0.0_dp
      DO i = 1, 3
         cell(i, i) = REAL(6 + i, dp)
      END DO
      CALL pw_grid_create(grid, comm, cell, grid_span=HALFSPACE, npts=[8, 6, 8], odd=.FALSE., fft_usage=.TRUE.)
      CALL pw_grid_create(fine_grid, comm, cell, grid_span=HALFSPACE, npts=[12, 10, 12], &
                          odd=.FALSE., fft_usage=.TRUE.)
      CALL fine_density%create(fine_grid)
      CALL fine_filtered%create(fine_grid)
      CALL fine_point%create(fine_grid)
      fine_point%array = 0.0_dp
      point = fine_grid%bounds(2, :)
      IF (ALL(point >= fine_grid%bounds_local(1, :)) .AND. ALL(point <= fine_grid%bounds_local(2, :))) THEN
         fine_point%array(point(1), point(2), point(3)) = 1.0_dp/fine_grid%dvol
      END IF
      CALL pw_pool_create(pool, grid)
      nrows = PRODUCT(grid%npts)
      ALLOCATE (vrho(nrows, 2), vgrad(nrows, 3, 2), vtau(nrows, 2))
      CALL density(1)%create(grid)
      CALL kinetic(1)%create(grid)
      CALL smooth(1)%create(grid)
      CALL potential(1)%create(grid)
      CALL tau_potential(1)%create(grid)
      DO k = grid%bounds_local(1, 3), grid%bounds_local(2, 3)
         DO j = grid%bounds_local(1, 2), grid%bounds_local(2, 2)
            DO i = grid%bounds_local(1, 1), grid%bounds_local(2, 1)
               phase = 2.0_dp*pi*(REAL(i, dp)/8.0_dp + REAL(j, dp)/6.0_dp + REAL(k, dp)/8.0_dp)
               density(1)%array(i, j, k) = 0.2_dp + 0.1_dp*SIN(phase)
               kinetic(1)%array(i, j, k) = 0.3_dp + 0.15_dp*COS(phase)
            END DO
         END DO
      END DO
      CALL qs_upf_inverse_xc_transfer_real(density(1), fine_density)
      modes = [xc_rho_nn10, xc_rho_nn50, xc_rho_spline2_smooth, xc_rho_spline3_smooth]
      DO mode = 1, SIZE(modes)
         CALL xc_pw_smooth(density(1), smooth(1), modes(mode))
         CALL qs_upf_inverse_xc_pw_fields(smooth, reference, kinetic)
         CALL qs_upf_inverse_xc_pw_fields(density, fields, kinetic, modes(mode), include_gradient=.FALSE.)
         CALL check_error(MAXVAL(ABS(fields(:, :, :, 1) - reference(:, :, :, 1))), 0.0_dp, 3.e-14_dp)
         CALL check_error(MAXVAL(ABS(fields(:, :, :, 2:4))), 0.0_dp, 0.0_dp)
         CALL check_error(MAXVAL(ABS(fields(:, :, :, 5) - reference(:, :, :, 5))), 0.0_dp, 0.0_dp)
         CALL qs_upf_inverse_xc_pw_fields(density, fields, kinetic, modes(mode))
         CALL check_error(MAXVAL(ABS(fields - reference)), 0.0_dp, 3.e-14_dp)
         CALL qs_upf_inverse_xc_pw_adjoint(fields, potential, tau_potential, modes(mode))
         projected = pw_integral_ab(density(1), potential(1)) + pw_integral_ab(kinetic(1), tau_potential(1))
         CALL check_error(projected, SUM(fields**2), 3.e-13_dp)
         ! The two equal model spin covectors contract into one restricted channel.
         vrho(:, :) = SPREAD(RESHAPE(fields(:, :, :, 1), [nrows]), 2, 2)
         vtau(:, :) = SPREAD(RESHAPE(fields(:, :, :, 5), [nrows]), 2, 2)
         DO axis = 1, 3
            vgrad(:, axis, :) = SPREAD(RESHAPE(fields(:, :, :, 1 + axis), [nrows]), 2, 2)
         END DO
         CALL build_vxc_from_feature_grads(skala_density, skala_tau, density_ptr, pool, &
                                           vrho, vgrad, vtau, xc_deriv_pw, &
                                           global_grid_layout=.TRUE., xc_smooth_id=modes(mode))
         CALL check_error(MAXVAL(ABS(skala_density(1)%array - potential(1)%array)), 0.0_dp, 3.e-13_dp)
         CALL check_error(MAXVAL(ABS(skala_tau(1)%array - tau_potential(1)%array)), 0.0_dp, 3.e-13_dp)
         CALL pool%give_back_pw(skala_density(1))
         CALL pool%give_back_pw(skala_tau(1))
         DEALLOCATE (skala_density, skala_tau)

         ! A point outside the coarse index range still has a coarse-grid adjoint.
         CALL xc_pw_smooth(fine_density, fine_filtered, modes(mode))
         sampled = pw_integral_ab(fine_filtered, fine_point)
         CALL xc_pw_smooth(fine_point, fine_filtered, modes(mode))
         CALL qs_upf_inverse_xc_transfer_real(fine_filtered, potential(1))
         CALL check_error(pw_integral_ab(density(1), potential(1)), sampled, 3.e-13_dp)
      END DO
      CALL check_smoothed_xc_kernel(density(1), kinetic(1), pool)
      CALL density(1)%release()
      CALL kinetic(1)%release()
      CALL smooth(1)%release()
      CALL potential(1)%release()
      CALL tau_potential(1)%release()
      CALL fine_density%release()
      CALL fine_filtered%release()
      CALL fine_point%release()
      CALL pw_grid_release(fine_grid)
      CALL pw_pool_release(pool)
      CALL pw_grid_release(grid)
   END SUBROUTINE check_upf_xc_smoothing

! **************************************************************************************************
!> \brief Compare the filtered semilocal kernel with differences of the native XC potential.
!> \param base Positive nonuniform source density.
!> \param perturbation Independent nonuniform direction.
!> \param pool Native PW pool of both fields.
! **************************************************************************************************
   SUBROUTINE check_smoothed_xc_kernel(base, perturbation, pool)
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: base, perturbation
      TYPE(pw_pool_type), POINTER                        :: pool

      INTEGER                                            :: functional, spin, trial
#if defined(__LIBXC)
      INTEGER, PARAMETER                                 :: functional_count = 2
#else
      INTEGER, PARAMETER                                 :: functional_count = 1
#endif
      REAL(dp), PARAMETER                                :: step = 1.0e-5_dp
      REAL(dp)                                           :: energy, error, shift, virial(3, 3)
      TYPE(pw_c1d_gs_type)                               :: direction_g(2)
      TYPE(pw_c1d_gs_type), POINTER                      :: rho_g(:)
      TYPE(pw_r3d_rs_type)                               :: direction(2), direction_tau(2), minus(2), minus_tau(2)
      TYPE(pw_r3d_rs_type), POINTER                      :: rho(:), tau(:), weights, v(:), vt(:), kernel(:), kernel_tau(:)
      TYPE(section_type), POINTER                        :: definition
      TYPE(section_vals_type), POINTER                   :: xc

      NULLIFY (definition, xc, weights, v, vt, kernel, kernel_tau)
      CALL create_xc_section(definition)
      CALL section_vals_create(xc, definition)
      CALL section_vals_val_set(xc, 'XC_GRID%XC_SMOOTH_RHO', i_val=xc_rho_nn10)
      ALLOCATE (rho(2), rho_g(2), tau(2))
      DO spin = 1, 2
         CALL pool%create_pw(rho(spin))
         CALL pool%create_pw(rho_g(spin))
         CALL pool%create_pw(tau(spin))
         CALL pool%create_pw(direction(spin))
         CALL pool%create_pw(direction_g(spin))
         CALL pool%create_pw(direction_tau(spin))
         CALL pool%create_pw(minus(spin))
         CALL pool%create_pw(minus_tau(spin))
         direction(spin)%array = perturbation%array/REAL(spin, dp)
         direction_tau(spin)%array = 0.15_dp*direction(spin)%array
         CALL pw_transfer(direction(spin), direction_g(spin))
      END DO
      DO functional = 1, functional_count
         CALL section_vals_val_set(xc, 'XC_FUNCTIONAL%PBE%_SECTION_PARAMETERS_', l_val=functional == 1)
         IF (functional == 2) THEN
            CALL section_vals_val_set(xc, 'XC_FUNCTIONAL%MGGA_X_R2SCAN%_SECTION_PARAMETERS_', l_val=.TRUE.)
            CALL section_vals_val_set(xc, 'XC_FUNCTIONAL%MGGA_C_R2SCAN%_SECTION_PARAMETERS_', l_val=.TRUE.)
         END IF
         DO trial = 0, 2
            shift = 0.0_dp
            IF (trial > 0) shift = REAL(2*trial - 3, dp)*step
            DO spin = 1, 2
               rho(spin)%array = base%array/REAL(spin, dp) + shift*direction(spin)%array
               tau(spin)%array = 0.12_dp + 0.25_dp*base%array/REAL(spin, dp) + shift*direction_tau(spin)%array
               CALL pw_transfer(rho(spin), rho_g(spin))
            END DO
            IF (trial == 0) THEN
               CALL qs_upf_inverse_xc_grid_kernel(xc, rho, rho_g, tau, direction, direction_g, direction_tau, &
                                                  weights, pool, kernel, kernel_tau)
               CYCLE
            END IF
            CALL xc_vxc_pw_create(v, vt, energy, rho, rho_g, tau, xc, weights, pool, .FALSE., virial, &
                                  project_density_potential=.FALSE.)
            DO spin = 1, 2
               IF (trial == 1) THEN
                  minus(spin)%array = v(spin)%array
                  minus_tau(spin)%array = 0.0_dp
                  IF (ASSOCIATED(vt)) minus_tau(spin)%array = vt(spin)%array
               ELSE
                  error = MAXVAL(ABS((v(spin)%array - minus(spin)%array)/(2.0_dp*step) - kernel(spin)%array))
                  CALL check_error(error, 0.0_dp, 3.e-8_dp)
                  IF (functional == 2 .AND. .NOT. ASSOCIATED(vt)) ERROR STOP 'Missing meta-GGA tau potential'
                  IF (ASSOCIATED(vt)) THEN
                     error = MAXVAL(ABS((vt(spin)%array - minus_tau(spin)%array)/(2.0_dp*step) - kernel_tau(spin)%array))
                     CALL check_error(error, 0.0_dp, 3.e-8_dp)
                  END IF
               END IF
               CALL pool%give_back_pw(v(spin))
               IF (ASSOCIATED(vt)) CALL pool%give_back_pw(vt(spin))
            END DO
            DEALLOCATE (v)
            IF (ASSOCIATED(vt)) DEALLOCATE (vt)
            NULLIFY (v, vt)
         END DO
         DO spin = 1, 2
            CALL pool%give_back_pw(kernel(spin))
            CALL pool%give_back_pw(kernel_tau(spin))
         END DO
         DEALLOCATE (kernel, kernel_tau)
      END DO
      DO spin = 1, 2
         CALL pool%give_back_pw(rho(spin))
         CALL pool%give_back_pw(rho_g(spin))
         CALL pool%give_back_pw(tau(spin))
         CALL pool%give_back_pw(direction(spin))
         CALL pool%give_back_pw(direction_g(spin))
         CALL pool%give_back_pw(direction_tau(spin))
         CALL pool%give_back_pw(minus(spin))
         CALL pool%give_back_pw(minus_tau(spin))
      END DO
      DEALLOCATE (rho, rho_g, tau)
      CALL section_vals_release(xc)
      CALL section_release(definition)
   END SUBROUTINE check_smoothed_xc_kernel

! **************************************************************************************************
!> \brief Preserve even-grid tensor-product Nyquist modes and the real-density/potential adjoint.
!> \param partial Keep one Nyquist extent unchanged while refining the other two directions.
! **************************************************************************************************
   SUBROUTINE check_real_interpolation(partial)
      LOGICAL, INTENT(IN)                                :: partial

      COMPLEX(dp)                                        :: phase
      INTEGER                                            :: alias, axes, fine_sizes(3), ig, nalias, &
                                                            nsplit
      LOGICAL                                            :: boundary(3)
      REAL(dp)                                           :: cell(3, 3), error, expected, left, &
                                                            point(3), right, vectors(3, 8)
      TYPE(mp_comm_type)                                 :: comm
      TYPE(pw_c1d_gs_type)                               :: coarse, fine, restricted, test
      TYPE(pw_grid_type), POINTER                        :: coarse_grid, fine_grid

      NULLIFY (coarse_grid, fine_grid)
      CALL comm%set_handle(mp_comm_world%get_handle())
      cell = 0.0_dp
      DO axes = 1, 3
         cell(axes, axes) = 7.0_dp
      END DO
      fine_sizes = [9, 9, 9]
      IF (partial) fine_sizes(1) = 6
      CALL pw_grid_create(fine_grid, comm, cell, grid_span=FULLSPACE, npts=fine_sizes, &
                          odd=.NOT. partial, fft_usage=.FALSE.)
      ! As in Quickstep, reference-aligned coarse grids use ray distribution.
      ! The fine grid may still use a blocked layout when its x extent is small.
      CALL pw_grid_create(coarse_grid, comm, cell, grid_span=FULLSPACE, npts=[6, 6, 6], &
                          odd=.FALSE., fft_usage=.FALSE., blocked=do_pw_grid_blocked_false, ref_grid=fine_grid)
      CALL coarse%create(coarse_grid)
      CALL fine%create(fine_grid)
      CALL test%create(fine_grid)
      CALL restricted%create(coarse_grid)
      ! Real atomic Fourier kernels need the same independent boundary signs.
      point = MATMUL(coarse_grid%dh, REAL([1, 2, -1], dp))
      DO ig = 1, SIZE(coarse%array)
         CALL pw_real_fourier_aliases(coarse_grid, ig, vectors, nalias)
         boundary = ABS(coarse_grid%g_hat(:, ig)) == 3
         IF (nalias /= 2**COUNT(boundary)) ERROR STOP "Missing real Fourier boundary aliases."
         phase = 0.0_dp
         DO alias = 1, nalias
            phase = phase + EXP(CMPLX(0.0_dp, DOT_PRODUCT(vectors(:, alias), point), dp))/REAL(nalias, dp)
         END DO
         IF (ABS(phase - EXP(CMPLX(0.0_dp, DOT_PRODUCT(coarse_grid%g(:, ig), point), dp))) > 3.e-14_dp) THEN
            ERROR STOP "Fourier boundary aliases changed the native grid phase."
         END IF
         IF (boundary(1) .AND. boundary(2)) THEN
            IF (ABS(SUM(vectors(1, :nalias)*vectors(2, :nalias))/REAL(nalias, dp)) > 3.e-14_dp) THEN
               ERROR STOP "A mixed odd kernel survived independent Nyquist reflection."
            END IF
         END IF
      END DO
      DO axes = 1, 3
         nsplit = axes
         IF (partial) nsplit = axes - 1
         ! The sampled product of one, two or three Nyquist cosines has a single
         ! coarse coefficient. Each refined direction splits into two signs.
         coarse%array = 0.0_dp
         DO ig = 1, SIZE(coarse%array)
            IF (ANY(coarse_grid%g_hat(:axes, ig) /= -3)) CYCLE
            IF (ANY(coarse_grid%g_hat(axes + 1:3, ig) /= 0)) CYCLE
            coarse%array(ig) = 1.0_dp
         END DO
         fine%array = (17.0_dp, -3.0_dp)
         CALL pw_real_transfer(coarse, fine)
         error = 0.0_dp
         DO ig = 1, SIZE(fine%array)
            expected = 0.0_dp
            IF (ALL(ABS(fine_grid%g_hat(:axes, ig)) == 3) .AND. &
                ALL(fine_grid%g_hat(axes + 1:3, ig) == 0)) THEN
               IF (.NOT. partial .OR. fine_grid%g_hat(1, ig) == -3) expected = 0.5_dp**nsplit
            END IF
            error = MAX(error, ABS(fine%array(ig) - expected))
         END DO
         CALL comm%max(error)
         IF (error > 3.e-14_dp) ERROR STOP "Real Fourier interpolation lost a Nyquist sign."
      END DO
      DO ig = 1, SIZE(coarse%array)
         coarse%array(ig) = CMPLX(COS(REAL(SUM(coarse_grid%g_hat(:, ig)), dp)), &
                                  SIN(REAL(DOT_PRODUCT([2, 3, 5], coarse_grid%g_hat(:, ig)), dp)), dp)
      END DO
      DO ig = 1, SIZE(test%array)
         test%array(ig) = CMPLX(SIN(REAL(SUM(fine_grid%g_hat(:, ig)), dp)), &
                                COS(REAL(DOT_PRODUCT([3, 5, 7], fine_grid%g_hat(:, ig)), dp)), dp)
      END DO
      fine%array = test%array
      CALL pw_real_transfer(coarse, fine, additive=.TRUE.)
      fine%array = fine%array - test%array
      CALL pw_real_transfer(test, restricted)
      left = REAL(SUM(CONJG(test%array)*fine%array), dp)
      right = REAL(SUM(CONJG(restricted%array)*coarse%array), dp)
      CALL comm%sum(left)
      CALL comm%sum(right)
      IF (ABS(left - right) > 1.e-13_dp*MAX(1.0_dp, ABS(left), ABS(right))) THEN
         ERROR STOP "Real Fourier density and potential transfers are not adjoint."
      END IF
      CALL coarse%release()
      CALL fine%release()
      CALL test%release()
      CALL restricted%release()
      CALL pw_grid_release(coarse_grid)
      CALL pw_grid_release(fine_grid)
   END SUBROUTINE check_real_interpolation

! **************************************************************************************************
!> \brief Preserve three-digit cube exponents and parse adjacent values across z-slices.
! **************************************************************************************************
   SUBROUTINE check_cube_values()
      CHARACTER(LEN=13)                                  :: value
      CHARACTER(LEN=79)                                  :: values
      CHARACTER(LEN=92)                                  :: slice_values
      INTEGER                                            :: axis, i, info, j, k, layout, sample, &
                                                            window
      REAL(dp) :: cell(3, 3), derivative(3), expected, FRACTION(3), gradient(3), inverse(3, 3), &
         minus, plus, point(3), position(3), sampled_value, shifted(3), step
      REAL(KIND=dp), DIMENSION(6)                        :: buffer, reference
      REAL(KIND=dp), DIMENSION(7)                        :: slice_buffer, slice_reference
      TYPE(mp_comm_type)                                 :: comm
      TYPE(pw_grid_type), POINTER                        :: grid
      TYPE(pw_r3d_rs_type), POINTER                      :: potential

      WRITE (value, cube_value_format) 3.33004E-101_dp
      IF (INDEX(value, "E-101") <= 0) THEN
         ERROR STOP "Cube value format must preserve three-digit negative exponents."
      END IF

      reference = [0.0_dp, -0.91458E-35_dp, 3.33004E-101_dp, &
                   -1.32487E-103_dp, 1.0_dp, -1.0_dp]
      WRITE (values(1:78), cube_values_format) reference
      values(79:79) = NEW_LINE("C")
      IF (INDEX(values, "E-101") <= 0) THEN
         ERROR STOP "Cube line format must preserve E-101."
      END IF
      IF (INDEX(values, "E-103") <= 0) THEN
         ERROR STOP "Cube line format must preserve E-103."
      END IF
      CALL cube_read_values(values, buffer)
      IF (MAXVAL(ABS(buffer - reference)) > 1.0E-12_dp) THEN
         ERROR STOP "Cube reader must parse adjacent explicit-exponent values."
      END IF

      slice_reference(1:6) = reference
      slice_reference(7) = 0.27183E-123_dp
      WRITE (slice_values(1:78), cube_values_format) slice_reference(1:6)
      slice_values(79:79) = NEW_LINE("C")
      WRITE (slice_values(80:92), cube_value_format) slice_reference(7)
      CALL cube_read_values(slice_values, slice_buffer)
      IF (MAXVAL(ABS(slice_buffer - slice_reference)) > 1.0E-12_dp) THEN
         ERROR STOP "Cube reader must parse multi-line z-slices."
      END IF
      ! Cube potentials use the same interpolated energy and Cartesian derivative.
      cell = RESHAPE([7.0_dp, 0.0_dp, 0.0_dp, 0.6_dp, 6.0_dp, 0.0_dp, 0.2_dp, 0.4_dp, 5.0_dp], [3, 3])
      DO layout = 1, 2
         IF (layout == 1) THEN
            CALL comm%set_handle(mp_comm_world%get_handle())
         ELSE
            CALL comm%set_handle(mp_comm_self%get_handle())
         END IF
         NULLIFY (grid)
         CALL pw_grid_create(grid, comm, cell, grid_span=FULLSPACE, npts=[6, 5, 4], odd=.FALSE., fft_usage=.FALSE.)
         ALLOCATE (potential)
         CALL potential%create(grid)
         DO k = grid%bounds_local(1, 3), grid%bounds_local(2, 3)
            DO j = grid%bounds_local(1, 2), grid%bounds_local(2, 2)
               DO i = grid%bounds_local(1, 1), grid%bounds_local(2, 1)
                  point = REAL([i, j, k] - grid%bounds(1, :), dp)
                  potential%array(i, j, k) = 0.7_dp + DOT_PRODUCT([0.2_dp, 0.3_dp, -0.1_dp], point**2) + &
                                             0.07_dp*point(1)*point(2) - 0.04_dp*point(2)*point(3) + 0.01_dp*PRODUCT(point)
               END DO
            END DO
         END DO
         inverse = grid%dh
         CALL invmat(inverse, info)
         IF (info /= 0) ERROR STOP "Invalid external-potential test cell."
         DO sample = 1, 2
            point = [1.23_dp, 1.34_dp, 1.45_dp]
            IF (sample == 2) point = [2.0_dp*grid%npts(1) + 0.61_dp, &
                                      -REAL(grid%npts(2), dp) + 1.32_dp, REAL(grid%npts(3), dp) - 0.45_dp]
            position = MATMUL(grid%dh, point)
            CALL interpolate_external_potential(position, potential, sampled_value, derivative, .TRUE.)
            IF (sample == 1) THEN
               fraction = point - REAL(FLOOR(point), dp)
               expected = 0.7_dp + DOT_PRODUCT([0.2_dp, 0.3_dp, -0.1_dp], point**2 + fraction*(1.0_dp - fraction)) + &
                          0.07_dp*point(1)*point(2) - 0.04_dp*point(2)*point(3) + 0.01_dp*PRODUCT(point)
               gradient = [0.2_dp*REAL(2*FLOOR(point(1)) + 1, dp) + 0.07_dp*point(2) + 0.01_dp*point(2)*point(3), &
                           0.3_dp*REAL(2*FLOOR(point(2)) + 1, dp) + 0.07_dp*point(1) - 0.04_dp*point(3) + &
                           0.01_dp*point(1)*point(3), &
                           -0.1_dp*REAL(2*FLOOR(point(3)) + 1, dp) - 0.04_dp*point(2) + 0.01_dp*point(1)*point(2)]
               gradient = MATMUL(TRANSPOSE(inverse), gradient)
               IF (ABS(sampled_value - expected) > 3.e-13_dp .OR. MAXVAL(ABS(derivative - gradient)) > 3.e-13_dp) THEN
                  ERROR STOP "External-potential interpolation lost the sheared-cell energy or gradient."
               END IF
            END IF
            DO window = 1, 2
               step = 1.e-4_dp/REAL(window, dp)
               DO axis = 1, 3
                  shifted = position
                  shifted(axis) = shifted(axis) + step
                  CALL interpolate_external_potential(shifted, potential, plus)
                  shifted(axis) = position(axis) - step
                  CALL interpolate_external_potential(shifted, potential, minus)
                  IF (ABS((plus - minus)/(2.0_dp*step) - derivative(axis)) > 2.e-8_dp) THEN
                     ERROR STOP "External-potential derivative differs from its interpolated energy."
                  END IF
               END DO
            END DO
         END DO
         CALL potential%release()
         DEALLOCATE (potential)
         CALL pw_grid_release(grid)
      END DO
   END SUBROUTINE check_cube_values

! **************************************************************************************************
!> \brief Compare identical, permuted, redistributed, partial and empty PW grids.
!> \param mode Layout variant; mode two changes only one rank's ordering.
!> \author Thomas D. Kuehne (tkuehne@cp2k.org)
! **************************************************************************************************
   SUBROUTINE check_copy(mode)
      INTEGER, INTENT(IN)                                :: mode

      COMPLEX(dp), PARAMETER :: untouched = (-999.0_dp, 17.0_dp)

      COMPLEX(dp)                                        :: expected
      INTEGER                                            :: dimensions(1), i, j, nsource, ntarget, &
                                                            owner, rank
      TYPE(pw_c1d_gs_type)                               :: destination, source
      TYPE(pw_grid_type), TARGET                         :: source_grid, target_grid

      rank = para_env%mepos
      nsource = 7
      IF (mode == 6 .OR. (mode >= 7 .AND. rank == 0)) nsource = 0
      owner = rank
      IF (mode == 3 .OR. mode == 8) owner = MOD(rank + 1, para_env%num_pe)
      ntarget = nsource
      IF (mode == 4) ntarget = nsource - 1
      IF (mode == 5) ntarget = nsource + 1
      IF (mode == 8) ntarget = MERGE(0, 7, owner == 0)
      ALLOCATE (source_grid%g_hat(3, nsource), target_grid%g_hat(3, ntarget))
      ALLOCATE (source%array(nsource), destination%array(ntarget))
      source%pw_grid => source_grid
      destination%pw_grid => target_grid
      dimensions = [para_env%num_pe]
      CALL source_grid%para%group%create(para_env, 1, dimensions)
      target_grid%para%group = source_grid%para%group
      source_grid%ngpts_cut_local = nsource
      target_grid%ngpts_cut_local = ntarget
      DO i = 1, nsource
         source_grid%g_hat(:, i) = [i, rank, 0]
         source%array(i) = CMPLX(i + 11*rank, -i - 3*rank, dp)
      END DO
      DO i = 1, ntarget
         j = i
         IF (mode == 1 .OR. (mode == 2 .AND. rank == 0)) j = ntarget + 1 - i
         target_grid%g_hat(:, i) = [j, owner, 0]
      END DO
      destination%array(:) = untouched
      CALL pw_copy_match(source, destination)
      DO i = 1, ntarget
         j = target_grid%g_hat(1, i)
         expected = CMPLX(j + 11*owner, -j - 3*owner, dp)
         IF (mode == 5 .AND. j > nsource) expected = untouched
         IF (ABS(destination%array(i) - expected) > 1.e-14_dp) ERROR STOP 'Incorrect PW copy'
      END DO
      IF (mode == 5) THEN
         CALL pw_copy(source, destination)
         IF (ANY(destination%array(:nsource) /= source%array) .OR. &
             destination%array(ntarget) /= 0.0_dp) ERROR STOP 'Incorrect PW overwrite'
      END IF
      DEALLOCATE (source%array, destination%array)
      CALL source_grid%para%group%free()
   END SUBROUTINE check_copy

! **************************************************************************************************
!> \brief Check PW multiplicity and volume normalization, and auxiliary-density Poisson routes.
! **************************************************************************************************
   SUBROUTINE check_pw_integral_and_poisson()
      REAL(dp), PARAMETER :: step = 5.e-4_dp, &
         strain_steps(4) = [step, -step, 2.0_dp*step, -2.0_dp*step]

      COMPLEX(dp)                                        :: a_mode, b_mode
      INTEGER                                            :: axis, ig, sign_step, span, sumtype, test
      REAL(dp) :: cell(3, 3), derivative, energy, energy_dv, energy_v, energy_v_dv, expected, &
         expected_real, kernel, reference, shifted_g(3), strain_energy(4), stress(3, 3), &
         stress_dv(3, 3), stress_only(3, 3), stress_v(3, 3), stress_v_dv(3, 3)
      TYPE(mp_comm_type)                                 :: comm
      TYPE(pw_c1d_gs_type)                               :: a, b, gradient(3), gradient_v_dv(3), &
                                                            potential, potential_v_dv
      TYPE(pw_grid_type), POINTER                        :: grid
      TYPE(pw_poisson_type)                              :: poisson
      TYPE(pw_r1d_gs_type)                               :: a_real, b_real

      CALL comm%set_handle(mp_comm_world%get_handle())
      ! Nonunit volume and an appreciable truncation factor exp(-pi).
      cell = 0.0_dp
      DO axis = 1, 3
         cell(axis, axis) = 7.0_dp**(1.0_dp/3.0_dp)
      END DO
      a_mode = CMPLX(0.4_dp, -0.2_dp, dp)
      b_mode = CMPLX(-0.7_dp, 0.9_dp, dp)
      NULLIFY (grid)
      DO span = FULLSPACE, HALFSPACE, -1
         ! Reciprocal inputs and outputs need no FFT. The production grid still
         ! supplies the distributed G0 ownership and half/full-space mappings.
         CALL pw_grid_create(grid, comm, cell, grid_span=span, npts=[5, 5, 5], odd=.TRUE., fft_usage=.FALSE.)
         CALL a%create(grid)
         CALL b%create(grid)
         CALL a_real%create(grid)
         CALL b_real%create(grid)
         CALL potential%create(grid)
         CALL potential_v_dv%create(grid)
         DO axis = 1, 3
            CALL gradient(axis)%create(grid)
            CALL gradient_v_dv(axis)%create(grid)
         END DO
         ALLOCATE (poisson%pw_pools(1))
         CALL pw_pool_create(poisson%pw_pools(1)%pool, grid)
         poisson%pw_level = 1
         poisson%parameters%solver = pw_poisson_mt
         poisson%parameters%periodic = [1, 1, 0]
         poisson%parameters%mt_alpha = 1.0_dp
         CALL pw_poisson_set(poisson, cell_hmat=cell, use_level=1)
         DO test = 1, 2
            ! Constant G0 only and G0 plus a Hermitian mode pair.
            a%array = 0.0_dp
            b%array = 0.0_dp
            expected = 0.0_dp
            expected_real = 0.0_dp
            reference = 0.0_dp
            IF (grid%have_g0) THEN
               a%array(1) = 2.0_dp
               b%array(1) = 3.0_dp
            END IF
            expected = 6.0_dp
            expected_real = 6.0_dp
            reference = 6.0_dp*mt2d_kernel_reference([0.0_dp, 0.0_dp, 0.0_dp], 3, cell(3, 3))
            IF (test /= 1) THEN
               DO ig = 1, grid%ngpts_cut_local
                  IF (ANY(grid%g_hat(2:3, ig) /= 0)) CYCLE
                  IF (grid%g_hat(1, ig) == 1) THEN
                     a%array(ig) = a_mode
                     b%array(ig) = b_mode
                  ELSE IF (grid%g_hat(1, ig) == -1) THEN
                     a%array(ig) = CONJG(a_mode)
                     b%array(ig) = CONJG(b_mode)
                  END IF
               END DO
               expected = expected + 2.0_dp*REAL(CONJG(a_mode)*b_mode, dp)
               expected_real = expected_real + 2.0_dp*REAL(a_mode, dp)*REAL(b_mode, dp)
               kernel = mt2d_kernel_reference([2.0_dp*pi/cell(1, 1), 0.0_dp, 0.0_dp], 3, cell(3, 3))
               reference = reference + 2.0_dp*REAL(CONJG(a_mode)*b_mode, dp)*kernel
            END IF
            a_real%array = REAL(a%array, dp)
            b_real%array = REAL(b%array, dp)
            ! G0 must have weight one, the mode pair weight two, and both must
            ! receive V=7 for an integral. just_sum retains the unscaled oracle.
            DO sumtype = do_accurate_sum, do_standard_sum
               CALL check_error(pw_integral_ab(a, b, sumtype=sumtype), 7.0_dp*expected, 3.e-14_dp)
               CALL check_error(pw_integral_ab(a, b, sumtype=sumtype, just_sum=.TRUE.), expected, 3.e-14_dp)
               CALL check_error(pw_integral_ab(a_real, b_real, sumtype=sumtype), 7.0_dp*expected_real, 3.e-14_dp)
               CALL check_error(pw_integral_ab(a_real, b_real, sumtype=sumtype, just_sum=.TRUE.), &
                                expected_real, 3.e-14_dp)
               CALL check_error(pw_integral_ab(a, b_real, sumtype=sumtype), 7.0_dp*expected_real, 3.e-14_dp)
               CALL check_error(pw_integral_ab(a_real, b, sumtype=sumtype), 7.0_dp*expected_real, 3.e-14_dp)
            END DO
            CALL pw_poisson_solve(poisson, a, energy, h_stress=stress, aux_density=b)
            CALL pw_poisson_solve(poisson, a, energy_dv, dvhartree=gradient, h_stress=stress_dv, aux_density=b)
            CALL pw_poisson_solve(poisson, a, energy_v, potential, h_stress=stress_v, aux_density=b)
            CALL pw_poisson_solve(poisson, a, energy_v_dv, potential_v_dv, dvhartree=gradient_v_dv, &
                                  h_stress=stress_v_dv, aux_density=b)
            CALL pw_poisson_solve(poisson, a, h_stress=stress_only, aux_density=b)
            CALL check_error(energy, 0.5_dp*7.0_dp*reference, 3.e-13_dp)
            CALL check_error(energy_dv, energy, 3.e-13_dp)
            CALL check_error(energy_v, energy, 3.e-13_dp)
            CALL check_error(energy_v_dv, energy, 3.e-13_dp)
            CALL check_error(MAXVAL(ABS(stress_dv - stress)), 0.0_dp, 3.e-11_dp)
            CALL check_error(MAXVAL(ABS(stress_v - stress)), 0.0_dp, 3.e-11_dp)
            CALL check_error(MAXVAL(ABS(stress_v_dv - stress)), 0.0_dp, 3.e-11_dp)
            CALL check_error(MAXVAL(ABS(stress_only - stress)), 0.0_dp, 3.e-11_dp)
            IF (test == 2) THEN
               ! Fixed Fourier charges Q=V*rho give E(epsilon)=V/2/(1+epsilon)
               ! times the strained kernel sum. This absolute production-virial
               ! oracle includes both G0/volume and the nonzero-mode kernel term.
               DO sign_step = 1, 4
                  shifted_g = [2.0_dp*pi/cell(1, 1)/(1.0_dp + strain_steps(sign_step)), 0.0_dp, 0.0_dp]
                  reference = 6.0_dp*mt2d_kernel_reference([0.0_dp, 0.0_dp, 0.0_dp], 3, cell(3, 3)) &
                              + 2.0_dp*REAL(CONJG(a_mode)*b_mode, dp)*mt2d_kernel_reference(shifted_g, 3, cell(3, 3))
                  strain_energy(sign_step) = 0.5_dp*7.0_dp*reference/(1.0_dp + strain_steps(sign_step))
               END DO
               derivative = (8.0_dp*(strain_energy(1) - strain_energy(2)) &
                             - strain_energy(3) + strain_energy(4))/(12.0_dp*step)
               CALL check_error(stress(1, 1) + derivative, 0.0_dp, 3.e-9_dp)
            END IF
            DO ig = 1, grid%ngpts_cut_local
               kernel = mt2d_kernel_reference(grid%g(:, ig), 3, cell(3, 3))
               CALL check_error(ABS(potential%array(ig) - kernel*a%array(ig)), 0.0_dp, 3.e-13_dp)
               CALL check_error(ABS(potential_v_dv%array(ig) - kernel*a%array(ig)), 0.0_dp, 3.e-13_dp)
               DO axis = 1, 3
                  CALL check_error(ABS(gradient(axis)%array(ig) &
                                       - CMPLX(0.0_dp, grid%g(axis, ig), dp)*kernel*a%array(ig)), 0.0_dp, 3.e-13_dp)
                  CALL check_error(ABS(gradient_v_dv(axis)%array(ig) &
                                       - CMPLX(0.0_dp, grid%g(axis, ig), dp)*kernel*a%array(ig)), 0.0_dp, 3.e-13_dp)
               END DO
            END DO
         END DO
         CALL poisson%release()
         CALL a%release()
         CALL b%release()
         CALL a_real%release()
         CALL b_real%release()
         CALL potential%release()
         CALL potential_v_dv%release()
         DO axis = 1, 3
            CALL gradient(axis)%release()
            CALL gradient_v_dv(axis)%release()
         END DO
         CALL pw_grid_release(grid)
      END DO
   END SUBROUTINE check_pw_integral_and_poisson

! **************************************************************************************************
!> \brief Compare native FFT kernels with newly constructed strained-cell kernels.
!>        Fixed reciprocal indices are matched across MPI layouts, including G0.
! **************************************************************************************************
   SUBROUTINE check_green_strain()
      REAL(dp), PARAMETER :: step = 3.e-4_dp, shifts(4) = [step, -step, 2.0_dp*step, -2.0_dp*step]

      INTEGER                                            :: alpha, axis, beta, ig, mask, sign_step, &
                                                            test
      REAL(dp) :: cell(3, 3), derivative, energy(2), norm, shifted_cell(3, 3), strain(3, 3), &
         strain_energy(4, 2), stress(3, 3, 2), volume_ratio
      REAL(dp), ALLOCATABLE                              :: response(:), shifted_values(:, :)
      TYPE(greens_fn_type)                               :: green, shifted_green
      TYPE(mp_comm_type)                                 :: comm
      TYPE(pw_c1d_gs_type)                               :: a, b, matched, weighted
      TYPE(pw_grid_type), POINTER                        :: dct_grid, grid, shifted_grid, super_grid
      TYPE(pw_poisson_parameter_type)                    :: parameters
      TYPE(pw_poisson_type)                              :: poisson
      TYPE(pw_pool_type), POINTER                        :: pool, shifted_pool

      CALL comm%set_handle(mp_comm_world%get_handle())
      NULLIFY (dct_grid, grid, shifted_grid, super_grid, pool, shifted_pool)
      ! x*K1(x) tends to one. The missing factor x in the I1 polynomial violates this limit.
      CALL check_error(1.e-6_dp*bessk1(1.e-6_dp), 1.0_dp, 2.e-11_dp)
      cell = 0.0_dp
      DO axis = 1, 3
         cell(axis, axis) = REAL(5 + 4*axis, dp)
      END DO
      DO mask = 0, 9
         parameters%solver = pw_poisson_analytic
         parameters%periodic = MERGE(1, 0, BTEST(mask, [0, 1, 2]))
         IF (mask == 8) THEN
            parameters%solver = pw_poisson_mt
            parameters%periodic = [1, 1, 0]
            parameters%mt_alpha = 1.0_dp
         ELSE IF (mask == 9) THEN
            parameters%solver = pw_poisson_mt
            parameters%periodic = 0
            parameters%mt_alpha = 1.0_dp
         END IF
         ! The isolated MT kernel uses the native FFT of its real-space long-range part.
         CALL pw_grid_create(grid, comm, cell, grid_span=HALFSPACE, npts=[5, 5, 5], odd=.TRUE., fft_usage=(mask == 9))
         CALL pw_pool_create(pool, grid)
         CALL pw_green_create(green, parameters, cell, pool, super_grid, dct_grid)
         CALL matched%create(grid)
         CALL a%create(grid)
         CALL b%create(grid)
         CALL weighted%create(grid)
         DO ig = 1, SIZE(a%array)
            norm = 1.0_dp + REAL(SUM(grid%g_hat(:, ig)**2), dp)
            a%array(ig) = CMPLX(1.0_dp, 0.3_dp*SIN(REAL(SUM(grid%g_hat(:, ig)), dp)), dp)/norm
            b%array(ig) = CMPLX(-0.2_dp, 0.5_dp*SIN(REAL(DOT_PRODUCT([1, 2, 3], grid%g_hat(:, ig)), dp)), dp)/norm
         END DO
         ALLOCATE (poisson%pw_pools(1))
         CALL pw_pool_create(poisson%pw_pools(1)%pool, grid)
         poisson%pw_level = 1
         poisson%parameters = parameters
         CALL pw_poisson_set(poisson, cell_hmat=cell, use_level=1)
         CALL pw_poisson_solve(poisson, a, energy(1), h_stress=stress(:, :, 1))
         CALL pw_poisson_solve(poisson, a, energy(2), h_stress=stress(:, :, 2), aux_density=b)
         weighted%array = green%influence_fn%array*a%array
         CALL check_error(energy(1), 0.5_dp*pw_integral_ab(a, weighted), 3.e-13_dp)
         CALL check_error(energy(2), 0.5_dp*pw_integral_ab(b, weighted), 3.e-13_dp)
         IF (mask == 0) THEN
            CALL pw_poisson_solve(poisson, a, energy(1), h_stress=stress(:, :, 1), greenfn=green%influence_fn)
            CALL check_error(energy(1), 0.5_dp*pw_integral_ab(a, weighted), 3.e-13_dp)
         END IF
         DO axis = 1, 3
            IF (parameters%periodic(axis) == 1) CYCLE
            CALL check_error(MAXVAL(ABS(stress(axis, :, :))), 0.0_dp, 0.0_dp)
            CALL check_error(MAXVAL(ABS(stress(:, axis, :))), 0.0_dp, 0.0_dp)
         END DO
         ALLOCATE (response(SIZE(matched%array)), shifted_values(SIZE(matched%array), 4))
         strain = 0.0_dp
         CALL pw_green_strain_response(green, strain, response)
         CALL check_error(MAXVAL(ABS(response)), 0.0_dp, 0.0_dp)
         DO alpha = 1, 3
            IF (parameters%periodic(alpha) == 0) CYCLE
            DO beta = 1, 3
               IF (parameters%periodic(beta) == 0) CYCLE
               strain = 0.0_dp
               strain(alpha, beta) = 1.0_dp
               CALL pw_green_strain_response(green, strain, response)
               DO sign_step = 1, 4
                  shifted_cell = cell + shifts(sign_step)*MATMUL(strain, cell)
                  CALL pw_grid_create(shifted_grid, comm, shifted_cell, grid_span=HALFSPACE, &
                                      npts=[5, 5, 5], odd=.TRUE., fft_usage=.FALSE.)
                  CALL pw_pool_create(shifted_pool, shifted_grid)
                  CALL pw_green_create(shifted_green, parameters, shifted_cell, shifted_pool, super_grid, dct_grid)
                  CALL pw_copy_match(shifted_green%influence_fn, matched)
                  shifted_values(:, sign_step) = REAL(matched%array, dp)
                  ! Hold Fourier charges V*rho fixed. The reciprocal integral
                  ! uses the original volume, so divide once by V_strained/V.
                  volume_ratio = 1.0_dp
                  IF (alpha == beta) volume_ratio = 1.0_dp + shifts(sign_step)
                  weighted%array = matched%array*a%array
                  strain_energy(sign_step, 1) = 0.5_dp*pw_integral_ab(a, weighted)/volume_ratio
                  strain_energy(sign_step, 2) = 0.5_dp*pw_integral_ab(b, weighted)/volume_ratio
                  CALL pw_green_release(shifted_green, shifted_pool)
                  CALL pw_pool_release(shifted_pool)
                  CALL pw_grid_release(shifted_grid)
               END DO
               DO test = 1, 2
                  derivative = (8.0_dp*(strain_energy(1, test) - strain_energy(2, test)) &
                                - strain_energy(3, test) + strain_energy(4, test))/(12.0_dp*step)
                  CALL check_error(stress(alpha, beta, test), -derivative, 3.e-9_dp)
               END DO
               DO ig = 1, SIZE(response)
                  derivative = (8.0_dp*(shifted_values(ig, 1) - shifted_values(ig, 2)) &
                                - shifted_values(ig, 3) + shifted_values(ig, 4))/(12.0_dp*step)
                  CALL check_error(response(ig), derivative, 3.e-9_dp)
               END DO
            END DO
         END DO
         DEALLOCATE (response, shifted_values)
         CALL matched%release()
         CALL a%release()
         CALL b%release()
         CALL weighted%release()
         CALL poisson%release()
         CALL pw_green_release(green, pool)
         CALL pw_pool_release(pool)
         CALL pw_grid_release(grid)
      END DO
   END SUBROUTINE check_green_strain

! **************************************************************************************************
!> \brief Compare the MT2D kernel and in-plane virial with independent strained-cell energies.
! **************************************************************************************************
   SUBROUTINE check_mt2d_kernel()
      REAL(dp), PARAMETER :: step = 3.e-4_dp, shifts(4) = [step, -step, 2.0_dp*step, -2.0_dp*step]

      INTEGER                                            :: alpha, beta, ip, iz, jp, miller(3), &
                                                            mode, modes(3, 5), sign_step
      REAL(dp)                                           :: derivative, g(3), h_inv(3, 3), kernel, &
                                                            lengths(3), shifted_g(3), &
                                                            shifted_kernel(4), strain(3, 3)

      ! Rotate the Cartesian vacuum axis. A shear within the periodic plane is
      ! permitted; only normal/periodic mixing is excluded by the MT2D kernel.
      DO iz = 1, 3
         ip = MOD(iz, 3) + 1
         jp = MOD(iz + 1, 3) + 1
         lengths = [9.0_dp, 13.0_dp, 17.0_dp]
         h_inv = 0.0_dp
         DO alpha = 1, 3
            h_inv(alpha, alpha) = 1.0_dp/lengths(alpha)
         END DO
         h_inv(ip, jp) = -1.4_dp/(lengths(ip)*lengths(jp))
         modes = 0
         ! G=0, both odd and even normal-only modes, and several modes with
         ! appreciable exponential truncation reject the periodic K**2 formula.
         modes(iz, 2:3) = [1, 2]
         modes(ip, 4) = 1
         modes(:, 5) = [1, 2, 1]
         DO mode = 1, SIZE(modes, 2)
            miller = modes(:, mode)
            g = 2.0_dp*pi*MATMUL(TRANSPOSE(h_inv), REAL(miller, dp))
            CALL mt2d_kernel_strain(g, iz, lengths(iz), kernel, strain)
            CALL check_error(kernel, mt2d_kernel_reference(g, iz, lengths(iz)), 3.e-14_dp)
            CALL check_error(MAXVAL(ABS(strain(iz, :))), 0.0_dp, 0.0_dp)
            CALL check_error(MAXVAL(ABS(strain(:, iz))), 0.0_dp, 0.0_dp)
            IF (mode <= 3) CALL check_error(MAXVAL(ABS(strain)), 0.0_dp, 0.0_dp)
            DO alpha = 1, 3
               IF (alpha == iz) CYCLE
               DO beta = 1, 3
                  IF (beta == iz) CYCLE
                  DO sign_step = 1, 4
                     ! g'=(I+epsilon_ab)^(-T) g, evaluated independently of dK.
                     shifted_g = g
                     IF (alpha == beta) THEN
                        shifted_g(beta) = g(beta)/(1.0_dp + shifts(sign_step))
                     ELSE
                        shifted_g(beta) = g(beta) - shifts(sign_step)*g(alpha)
                     END IF
                     shifted_kernel(sign_step) = mt2d_kernel_reference(shifted_g, iz, lengths(iz))
                  END DO
                  derivative = (8.0_dp*(shifted_kernel(1) - shifted_kernel(2)) &
                                - shifted_kernel(3) + shifted_kernel(4))/(12.0_dp*step)
                  CALL check_error(strain(alpha, beta), derivative, 3.e-9_dp)
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE check_mt2d_kernel

! **************************************************************************************************
!> \brief Evaluate the original MT2D influence kernel independently of its strain derivative.
!> \param g ...
!> \param iz ...
!> \param slab_size ...
!> \return ...
! **************************************************************************************************
   FUNCTION mt2d_kernel_reference(g, iz, slab_size) RESULT(kernel)
      REAL(dp), INTENT(IN)                               :: g(3)
      INTEGER, INTENT(IN)                                :: iz
      REAL(dp), INTENT(IN)                               :: slab_size
      REAL(dp)                                           :: kernel

      REAL(dp)                                           :: g2, g_parallel(3), radius

      g2 = SUM(g**2)
      IF (g2 == 0.0_dp) THEN
         kernel = pi*slab_size**2/2.0_dp
      ELSE
         g_parallel = g
         g_parallel(iz) = 0.0_dp
         radius = slab_size/2.0_dp
         kernel = fourpi/g2 - fourpi/g2*COS(g(iz)*radius)*EXP(-NORM2(g_parallel)*radius)
      END IF
   END FUNCTION mt2d_kernel_reference

! **************************************************************************************************
!> \brief Require the scaled error to be finite and within the supplied tolerance.
!> \param value ...
!> \param reference ...
!> \param tolerance ...
! **************************************************************************************************
   SUBROUTINE check_error(value, reference, tolerance)
      REAL(dp), INTENT(IN)                               :: value, reference, tolerance

      REAL(dp)                                           :: error

      error = ABS(value - reference)/MAX(1.0_dp, ABS(reference))
      max_error = MAX(max_error, error)
      IF (.NOT. (error <= tolerance)) THEN
         PRINT *, 'Coulomb check: value, reference, error, tolerance:', value, reference, error, tolerance
         ERROR STOP 'PW integral or MT2D kernel check failed'
      END IF
   END SUBROUTINE check_error

END PROGRAM pw_unittest
