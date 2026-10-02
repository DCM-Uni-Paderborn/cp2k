!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Scalar PAW one-center semilocal XC on radial and spherical quadrature grids.
!>        Evaluates AE minus compensated-PS densities, including the corresponding frozen cores.
! **************************************************************************************************
MODULE paw_xc
   USE atom_upf,                        ONLY: atom_upfpot_type,&
                                              upf_radial_weights
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE input_section_types,             ONLY: section_vals_type
   USE kinds,                           ONLY: dp
   USE lebedev,                         ONLY: get_number_of_lebedev_grid,&
                                              lebedev_grid
   USE mathconstants,                   ONLY: fourpi
   USE orbital_pointers,                ONLY: indco,&
                                              init_orbital_pointers,&
                                              nco,&
                                              ncoset
   USE orbital_transformation_matrices, ONLY: init_spherical_harmonics,&
                                              orbtramat
   USE paw_core_kinetic,                ONLY: paw_core_kinetic_init,&
                                              paw_xc_build_ae_core_tau => paw_core_kinetic_build_ae
   USE paw_onecenter,                   ONLY: paw_onecenter_type
   USE paw_types,                       ONLY: paw_xc_type
   USE splines,                         ONLY: spline3ders
   USE xc_atom,                         ONLY: fill_rho_set,&
                                              vxc_of_r_new,&
                                              xc_rho_set_atom_update
   USE xc_derivative_set_types,         ONLY: xc_derivative_set_type,&
                                              xc_dset_create,&
                                              xc_dset_release,&
                                              xc_dset_zero_all
   USE xc_derivatives,                  ONLY: xc_functionals_get_needs
   USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_type
   USE xc_rho_set_types,                ONLY: xc_rho_set_create,&
                                              xc_rho_set_release,&
                                              xc_rho_set_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'paw_xc'

   PUBLIC :: paw_xc_type, paw_xc_init, paw_xc_evaluate, paw_xc_collocate
   PUBLIC :: paw_xc_kinetic_adjoint
   PUBLIC :: paw_xc_build_ae_core_tau

CONTAINS

! **************************************************************************************************
!> \brief Prepare radial density products, gradients, frozen cores and angular quadrature.
!> \param pot Original scalar PAW UPF data.
!> \param nang Minimum number of angular points; the selected rule also resolves harmonic products.
!> \param grid Immutable per-kind one-center XC data.
! **************************************************************************************************
   SUBROUTINE paw_xc_init(pot, nang, grid)
      TYPE(atom_upfpot_type), INTENT(IN)                 :: pot
      INTEGER, INTENT(IN)                                :: nang
      TYPE(paw_xc_type), INTENT(OUT)                     :: grid

      INTEGER                                            :: a, axis, first, i, j, l, lm, m, nb, nr, &
                                                            ns, pair, powers(3), reduced(3), rule, &
                                                            s, side
      REAL(dp)                                           :: angular_factor, factor
      REAL(dp), ALLOCATABLE                              :: derivatives(:, :), values(:), waves(:, :)

      CPASSERT(pot%is_paw .AND. .NOT. pot%has_so .AND. ALLOCATED(pot%paw))
      CPASSERT(ALLOCATED(pot%qijl))
      nr = SIZE(pot%r)
      CALL paw_core_kinetic_init(pot, grid%core_kinetic)
      grid%nr = nr
      grid%nval = pot%paw%cutoff_index
      CPASSERT(grid%nval >= 4 .AND. grid%nval <= nr)
      nb = pot%number_of_proj
      grid%npair = nb*(nb + 1)/2
      grid%lmax = 2*MAXVAL(pot%lbeta)
      rule = MAX(get_number_of_lebedev_grid(n=nang), get_number_of_lebedev_grid(l=2*grid%lmax))
      grid%na = lebedev_grid(rule)%n
      ALLOCATE (grid%r(nr), grid%wr(nr), grid%wa(grid%na), grid%directions(3, grid%na), &
                grid%harmonics(grid%na, (grid%lmax + 1)**2), &
                grid%angular_gradient(3, grid%na, (grid%lmax + 1)**2))
      grid%r(:) = pot%r
      grid%wr(:) = upf_radial_weights(pot)*pot%r**2
      grid%wa(:) = fourpi*lebedev_grid(rule)%w
      grid%directions(:, :) = lebedev_grid(rule)%r
      CALL init_orbital_pointers(grid%lmax)
      CALL init_spherical_harmonics(grid%lmax, -1)
      grid%harmonics = 0.0_dp
      grid%angular_gradient = 0.0_dp
      DO l = 0, grid%lmax
         DO m = 1, 2*l + 1
            lm = l*l + m
            DO s = 1, nco(l)
               powers = indco(:, ncoset(l - 1) + s)
               factor = orbtramat(l)%slm(m, s)
               DO a = 1, grid%na
                  grid%harmonics(a, lm) = grid%harmonics(a, lm) + factor*PRODUCT(grid%directions(:, a)**powers)
                  DO axis = 1, 3
                     IF (powers(axis) == 0) CYCLE
                     reduced = powers
                     reduced(axis) = reduced(axis) - 1
                     grid%angular_gradient(axis, a, lm) = grid%angular_gradient(axis, a, lm) + &
                                                          factor*REAL(powers(axis), dp)*PRODUCT(grid%directions(:, a)**reduced)
                  END DO
               END DO
            END DO
            DO a = 1, grid%na
               ! Subtract the radial derivative of r**l*Y_lm to obtain the tangent gradient.
               grid%angular_gradient(:, a, lm) = grid%angular_gradient(:, a, lm) - &
                                                 REAL(l, dp)*grid%harmonics(a, lm)*grid%directions(:, a)
            END DO
         END DO
      END DO

      ALLOCATE (grid%radial(grid%nval, grid%npair, 0:grid%lmax, 2), &
                grid%radial_gradient(grid%nval, grid%npair, 0:grid%lmax, 2), &
                grid%core(nr, 2), grid%core_gradient(nr, 2))
      grid%radial = 0.0_dp
      grid%radial_gradient = 0.0_dp
      grid%core = 0.0_dp
      grid%core(:, 1) = pot%paw%rho_core
      IF (ALLOCATED(pot%rho_nlcc)) grid%core(:, 2) = pot%rho_nlcc
      DO side = 1, 2
         CALL spline3ders(pot%r, grid%core(:, side), pot%r, dynew=grid%core_gradient(:, side))
      END DO
      first = 1
      IF (pot%r(1) == 0.0_dp) first = 2
      ! Retain a few physical continuation samples for boundary derivatives, without
      ! differentiating an artificial jump to zero at the augmentation-sphere boundary.
      ns = MIN(nr, grid%nval + 3)
      CPASSERT(ns - first + 1 >= 4)
      ALLOCATE (values(ns))
      DO j = 1, nb
         DO i = 1, j
            pair = j*(j - 1)/2 + i
            DO l = ABS(pot%lbeta(i) - pot%lbeta(j)), pot%lbeta(i) + pot%lbeta(j), 2
               DO side = 1, 2
                  IF (side == 1) THEN
                     values(first:) = pot%paw%ae_wfc(first:ns, i)*pot%paw%ae_wfc(first:ns, j)/pot%r(first:ns)**2
                  ELSE
                     values(first:) = (pot%paw%ps_wfc(first:ns, i)*pot%paw%ps_wfc(first:ns, j) + &
                                       pot%qijl(first:ns, pair, l))/pot%r(first:ns)**2
                  END IF
                  CALL spline3ders(pot%r(first:ns), values(first:), pot%r(:grid%nval), &
                                   ynew=grid%radial(:, pair, l, side), &
                                   dynew=grid%radial_gradient(:, pair, l, side))
               END DO
            END DO
         END DO
      END DO
      ALLOCATE (grid%kinetic_radial(grid%nval, grid%npair, 0:grid%lmax, 2), &
                waves(grid%nval, nb), derivatives(grid%nval, nb))
      grid%kinetic_radial = 0.0_dp
      DO side = 1, 2
         DO i = 1, nb
            IF (side == 1) THEN
               values(first:) = pot%paw%ae_wfc(first:ns, i)/pot%r(first:ns)
            ELSE
               values(first:) = pot%paw%ps_wfc(first:ns, i)/pot%r(first:ns)
            END IF
            CALL spline3ders(pot%r(first:ns), values(first:), pot%r(:grid%nval), &
                             ynew=waves(:, i), dynew=derivatives(:, i))
         END DO
         DO j = 1, nb
            DO i = 1, j
               pair = j*(j - 1)/2 + i
               DO l = ABS(pot%lbeta(i) - pot%lbeta(j)), pot%lbeta(i) + pot%lbeta(j), 2
                  ! grad_Omega(Y_i).grad_Omega(Y_j) has the same Gaunt coefficients,
                  ! multiplied by [l_i(l_i+1)+l_j(l_j+1)-L(L+1)]/2.
                  angular_factor = 0.5_dp*REAL(pot%lbeta(i)*(pot%lbeta(i) + 1) + &
                                               pot%lbeta(j)*(pot%lbeta(j) + 1) - l*(l + 1), dp)
                  grid%kinetic_radial(first:, pair, l, side) = 0.5_dp*( &
                                                               derivatives(first:, i)*derivatives(first:, j) + &
                                                         angular_factor*waves(first:, i)*waves(first:, j)/pot%r(first:grid%nval)**2)
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE paw_xc_init

! **************************************************************************************************
!> \brief Contract a real symmetric spin-resolved projector density with the product Gaunt tensor.
!> \param onecenter Magnetic/radial mappings and angular integrals.
!> \param density Projector density matrix, one or two spin channels.
!> \param coefficients Radial-pair and real-harmonic coefficients for each spin.
! **************************************************************************************************
   SUBROUTINE density_coefficients(onecenter, density, coefficients)
      TYPE(paw_onecenter_type), INTENT(IN)               :: onecenter
      REAL(dp), INTENT(IN)                               :: density(:, :, :)
      REAL(dp), INTENT(OUT)                              :: coefficients(:, :, :)

      INTEGER                                            :: p, pair, q, spin

      CPASSERT(SIZE(density, 1) == onecenter%nproj .AND. SIZE(density, 2) == onecenter%nproj)
      CPASSERT(SIZE(density, 3) == 1 .OR. SIZE(density, 3) == 2)
      CPASSERT(ALL(SHAPE(coefficients) == [onecenter%npair, SIZE(onecenter%gaunt, 3), SIZE(density, 3)]))
      coefficients = 0.0_dp
      DO spin = 1, SIZE(density, 3)
         DO q = 1, onecenter%nproj
            DO p = 1, onecenter%nproj
               pair = onecenter%pair_index(p, q)
               coefficients(pair, :, spin) = coefficients(pair, :, spin) + onecenter%gaunt(p, q, :)*density(p, q, spin)
            END DO
         END DO
      END DO
   END SUBROUTINE density_coefficients

! **************************************************************************************************
!> \brief One radial shell of the spin density and Cartesian gradient, including frozen core.
!> \param grid Prepared radial and angular data.
!> \param coefficients Contracted projector-density coefficients.
!> \param ir Radial index with nonzero radius.
!> \param side AE (1) or compensated PS (2).
!> \param rho Angular point by spin density.
!> \param gradient Cartesian direction by angular point by spin gradient.
!> \param kinetic Optional valence partial-wave kinetic density; frozen core is not included.
! **************************************************************************************************
   SUBROUTINE collocate_coefficients(grid, coefficients, ir, side, rho, gradient, kinetic)
      TYPE(paw_xc_type), INTENT(IN)                      :: grid
      REAL(dp), INTENT(IN)                               :: coefficients(:, :, :)
      INTEGER, INTENT(IN)                                :: ir, side
      REAL(dp), INTENT(OUT)                              :: rho(:, :), gradient(:, :, :)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: kinetic(:, :)

      INTEGER                                            :: axis, l, lm, nlm, nspin, spin
      REAL(dp) :: radial_derivative(SIZE(coefficients, 2)), harmonics(SIZE(coefficients, 2)), &
         radial_field(grid%na)

      nspin = SIZE(coefficients, 3)
      nlm = SIZE(coefficients, 2)
      CPASSERT(side == 1 .OR. side == 2)
      CPASSERT(ir >= 1 .AND. ir <= grid%nr)
      CPASSERT(grid%r(ir) > 0.0_dp)
      CPASSERT(ALL(SHAPE(rho) == [grid%na, nspin]))
      CPASSERT(ALL(SHAPE(gradient) == [3, grid%na, nspin]))
      DO spin = 1, nspin
         harmonics = 0.0_dp
         radial_derivative = 0.0_dp
         IF (ir <= grid%nval) THEN
            DO lm = 1, nlm
               l = INT(SQRT(REAL(lm - 1, dp)))
               harmonics(lm) = SUM(coefficients(:, lm, spin)*grid%radial(ir, :, l, side))
               radial_derivative(lm) = SUM(coefficients(:, lm, spin)*grid%radial_gradient(ir, :, l, side))
            END DO
         END IF
         rho(:, spin) = MATMUL(grid%harmonics, harmonics) + grid%core(ir, side)/REAL(nspin, dp)
         radial_field = MATMUL(grid%harmonics, radial_derivative) + grid%core_gradient(ir, side)/REAL(nspin, dp)
         DO axis = 1, 3
            gradient(axis, :, spin) = grid%directions(axis, :)*radial_field + &
                                      MATMUL(grid%angular_gradient(axis, :, :), harmonics)/grid%r(ir)
         END DO
         IF (PRESENT(kinetic)) THEN
            CPASSERT(ALL(SHAPE(kinetic) == SHAPE(rho)))
            harmonics = 0.0_dp
            IF (ir <= grid%nval) THEN
               DO lm = 1, nlm
                  l = INT(SQRT(REAL(lm - 1, dp)))
                  harmonics(lm) = SUM(coefficients(:, lm, spin)*grid%kinetic_radial(ir, :, l, side))
               END DO
            END IF
            kinetic(:, spin) = MATMUL(grid%harmonics, harmonics)
         END IF
      END DO
   END SUBROUTINE collocate_coefficients

! **************************************************************************************************
!> \brief Public shell collocation for diagnostics and consumers needing reconstructed fields.
!> \param grid Prepared PAW XC grid.
!> \param onecenter Corresponding magnetic/radial mapping.
!> \param density Spin-resolved projector density.
!> \param ir Nonzero radial node.
!> \param side AE (1) or compensated PS (2).
!> \param rho Angular point by spin density including frozen core.
!> \param gradient Cartesian direction by angular point by spin gradient.
!> \param kinetic Optional valence partial-wave kinetic density, excluding frozen core and charge compensation.
! **************************************************************************************************
   SUBROUTINE paw_xc_collocate(grid, onecenter, density, ir, side, rho, gradient, kinetic)
      TYPE(paw_xc_type), INTENT(IN)                      :: grid
      TYPE(paw_onecenter_type), INTENT(IN)               :: onecenter
      REAL(dp), INTENT(IN)                               :: density(:, :, :)
      INTEGER, INTENT(IN)                                :: ir, side
      REAL(dp), INTENT(OUT)                              :: rho(:, :), gradient(:, :, :)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: kinetic(:, :)

      REAL(dp) :: coefficients(onecenter%npair, SIZE(onecenter%gaunt, 3), SIZE(density, 3))

      CPASSERT(grid%npair == onecenter%npair .AND. grid%lmax == onecenter%lmax)
      CALL density_coefficients(onecenter, density, coefficients)
      CALL collocate_coefficients(grid, coefficients, ir, side, rho, gradient, kinetic)
   END SUBROUTINE paw_xc_collocate

! **************************************************************************************************
!> \brief Adjoint of valence kinetic-density collocation on a single radial shell.
!> \param grid Prepared PAW radial and angular data.
!> \param onecenter Projector indices and Gaunt coefficients.
!> \param ir Nonzero radial node.
!> \param side AE (1) or PS (2), without charge compensation in either kinetic field.
!> \param weighted_potential Kinetic-density potential, including all integration weights.
!> \param potential Projector-density derivative from this shell, separately for each spin.
! **************************************************************************************************
   SUBROUTINE paw_xc_kinetic_adjoint(grid, onecenter, ir, side, weighted_potential, potential)
      TYPE(paw_xc_type), INTENT(IN)                      :: grid
      TYPE(paw_onecenter_type), INTENT(IN)               :: onecenter
      INTEGER, INTENT(IN)                                :: ir, side
      REAL(dp), INTENT(IN)                               :: weighted_potential(:, :)
      REAL(dp), INTENT(OUT)                              :: potential(:, :, :)

      INTEGER                                            :: l, lm, p, pair, q, spin
      REAL(dp) :: derivative(onecenter%npair, SIZE(onecenter%gaunt, 3)), &
         angular_value(SIZE(onecenter%gaunt, 3))

      CPASSERT(grid%npair == onecenter%npair .AND. grid%lmax == onecenter%lmax)
      CPASSERT(ir >= 1 .AND. ir <= grid%nr .AND. side >= 1 .AND. side <= 2)
      CPASSERT(SIZE(weighted_potential, 1) == grid%na)
      CPASSERT(ALL(SHAPE(potential) == [onecenter%nproj, onecenter%nproj, SIZE(weighted_potential, 2)]))
      potential = 0.0_dp
      IF (ir > grid%nval) RETURN
      DO spin = 1, SIZE(weighted_potential, 2)
         angular_value = MATMUL(TRANSPOSE(grid%harmonics), weighted_potential(:, spin))
         DO lm = 1, SIZE(angular_value)
            l = INT(SQRT(REAL(lm - 1, dp)))
            derivative(:, lm) = angular_value(lm)*grid%kinetic_radial(ir, :, l, side)
         END DO
         DO q = 1, onecenter%nproj
            DO p = 1, onecenter%nproj
               pair = onecenter%pair_index(p, q)
               potential(p, q, spin) = SUM(onecenter%gaunt(p, q, :)*derivative(pair, :))
            END DO
         END DO
      END DO
   END SUBROUTINE paw_xc_kinetic_adjoint

! **************************************************************************************************
!> \brief Semilocal PAW XC correction and its full projector-density derivative.
!> \param grid Prepared AE/PS radial and angular data.
!> \param onecenter Matching projector mapping and Gaunt tensor.
!> \param functional XC_FUNCTIONAL section, with activated functionals.
!> \param density Spin-resolved projector density, already summed over k points.
!> \param density_cutoff XC evaluator density cutoff.
!> \param gradient_cutoff XC evaluator gradient cutoff.
!> \param energy AE-minus-compensated-PS XC energy.
!> \param potential Derivative dE/d(density), separately for each spin.
!> \param core_kinetic Optional spin-summed frozen-core kinetic density on the UPF mesh, AE then PS.
!>        Required for tau-dependent XC when either frozen-core density is nonzero.
! **************************************************************************************************
   SUBROUTINE paw_xc_evaluate(grid, onecenter, functional, density, density_cutoff, gradient_cutoff, energy, potential, &
                              core_kinetic)
      TYPE(paw_xc_type), INTENT(IN)                      :: grid
      TYPE(paw_onecenter_type), INTENT(IN)               :: onecenter
      TYPE(section_vals_type), POINTER                   :: functional
      REAL(dp), INTENT(IN)                               :: density(:, :, :), density_cutoff, &
                                                            gradient_cutoff
      REAL(dp), INTENT(OUT)                              :: energy, potential(:, :, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: core_kinetic(:, :)

      INTEGER                                            :: axis, bounds(2, 3), first, ir, jr, l, &
                                                            lm, na, nlm, nr, nspin, p, pair, q, &
                                                            side, spin
      LOGICAL                                            :: lsd, use_tau
      REAL(dp)                                           :: batch_energy, sign
      REAL(dp), ALLOCATABLE :: angular_value(:), drho(:, :, :, :), kinetic_value(:), &
         radial_value(:), rho(:, :, :), tau(:, :, :), weights(:, :)
      REAL(dp), POINTER                                  :: vtau(:, :, :), vxc(:, :, :), &
                                                            vxg(:, :, :, :)
      REAL(dp) :: &
         coefficient_derivative(onecenter%npair, SIZE(onecenter%gaunt, 3), SIZE(density, 3)), &
         coefficients(onecenter%npair, SIZE(onecenter%gaunt, 3), SIZE(density, 3))
      TYPE(xc_derivative_set_type)                       :: derivatives
      TYPE(xc_rho_cflags_type)                           :: needs
      TYPE(xc_rho_set_type)                              :: rho_set

      CPASSERT(grid%npair == onecenter%npair .AND. grid%lmax == onecenter%lmax)
      CPASSERT(ALL(SHAPE(potential) == SHAPE(density)))
      CALL density_coefficients(onecenter, density, coefficients)
      nspin = SIZE(density, 3)
      lsd = nspin == 2
      na = grid%na
      nlm = SIZE(onecenter%gaunt, 3)
      needs = xc_functionals_get_needs(functional, lsd=lsd, calc_potential=.TRUE.)
      IF (needs%laplace_rho .OR. needs%laplace_rho_spin) THEN
         CPABORT("PAW one-center XC Laplacian dependence is not implemented")
      END IF
      use_tau = needs%tau .OR. needs%tau_spin
      IF (use_tau) THEN
         IF (PRESENT(core_kinetic)) THEN
            CPASSERT(ALL(SHAPE(core_kinetic) == [grid%nr, 2]))
            CPASSERT(ALL(ieee_is_finite(core_kinetic)) .AND. ALL(core_kinetic >= 0.0_dp))
         ELSE IF (.NOT. grid%core_kinetic%available .AND. ANY(grid%core /= 0.0_dp)) THEN
            CALL cp_abort(__LOCATION__, TRIM(grid%core_kinetic%failure_reason))
         END IF
      END IF
      bounds = 1
      bounds(2, 1) = na
      ! Batch radial shells to amortize the functional setup without changing
      ! their quadrature. Zero-weight padding keeps the final batch rectangular.
      nr = MIN(grid%nr, MAX(1, 2048/na))
      bounds(2, 2) = nr
      ALLOCATE (rho(na, nr, nspin), drho(4, na, nr, nspin), tau(na, nr, nspin), weights(na, nr), &
                vxc(na, nr, nspin), vxg(3, na, nr, nspin), vtau(na, nr, nspin), &
                angular_value(nlm), radial_value(nlm), kinetic_value(nlm))
      CALL xc_rho_set_create(rho_set, bounds, rho_cutoff=density_cutoff, drho_cutoff=gradient_cutoff, tau_cutoff=0.0_dp)
      CALL xc_rho_set_atom_update(rho_set, needs, nspin, bounds)
      CALL xc_dset_create(derivatives, local_bounds=bounds)
      energy = 0.0_dp
      coefficient_derivative = 0.0_dp
      DO side = 1, 2
         sign = REAL(3 - 2*side, dp)
         DO first = 1, grid%nr, nr
            rho = 0.0_dp
            drho = 0.0_dp
            tau = 0.0_dp
            weights = 0.0_dp
            DO jr = 1, nr
               ir = first + jr - 1
               IF (ir <= grid%nr) THEN
                  IF (grid%r(ir) > 0.0_dp) THEN
                     IF (use_tau) THEN
                        CALL collocate_coefficients(grid, coefficients, ir, side, rho(:, jr, :), &
                                                    drho(1:3, :, jr, :), tau(:, jr, :))
                        IF (PRESENT(core_kinetic)) THEN
                           tau(:, jr, :) = tau(:, jr, :) + core_kinetic(ir, side)/REAL(nspin, dp)
                        ELSE IF (grid%core_kinetic%available) THEN
                           tau(:, jr, :) = tau(:, jr, :) + grid%core_kinetic%kinetic(ir, side)/REAL(nspin, dp)
                        END IF
                     ELSE
                        CALL collocate_coefficients(grid, coefficients, ir, side, rho(:, jr, :), drho(1:3, :, jr, :))
                     END IF
                     drho(4, :, jr, :) = SQRT(SUM(drho(1:3, :, jr, :)**2, DIM=1))
                     weights(:, jr) = grid%wr(ir)*grid%wa
                  END IF
               END IF
               CALL fill_rho_set(rho_set, lsd, nspin, needs, rho, drho, tau, na, jr)
            END DO
            CALL xc_dset_zero_all(derivatives)
            vxc = 0.0_dp
            vxg = 0.0_dp
            vtau = 0.0_dp
            CALL vxc_of_r_new(functional, rho_set, derivatives, 1, needs, weights, lsd, na, nr, &
                              batch_energy, vxc, vxg, vtau)
            energy = energy + sign*batch_energy
            DO jr = 1, MIN(nr, grid%nval - first + 1)
               ir = first + jr - 1
               IF (grid%r(ir) == 0.0_dp) CYCLE
               DO spin = 1, nspin
                  angular_value(:) = MATMUL(TRANSPOSE(grid%harmonics), vxc(:, jr, spin))
                  radial_value = 0.0_dp
                  DO axis = 1, 3
                     angular_value(:) = angular_value + &
                                        MATMUL(TRANSPOSE(grid%angular_gradient(axis, :, :)), vxg(axis, :, jr, spin))/grid%r(ir)
                     radial_value(:) = radial_value + &
                                       MATMUL(TRANSPOSE(grid%harmonics), grid%directions(axis, :)*vxg(axis, :, jr, spin))
                  END DO
                  DO lm = 1, nlm
                     l = INT(SQRT(REAL(lm - 1, dp)))
                     coefficient_derivative(:, lm, spin) = coefficient_derivative(:, lm, spin) + sign*( &
                              angular_value(lm)*grid%radial(ir, :, l, side) + radial_value(lm)*grid%radial_gradient(ir, :, l, side))
                  END DO
                  IF (use_tau) THEN
                     ! vxc_of_r_new has already included radial and angular weights.
                     kinetic_value(:) = MATMUL(TRANSPOSE(grid%harmonics), vtau(:, jr, spin))
                     DO lm = 1, nlm
                        l = INT(SQRT(REAL(lm - 1, dp)))
                        coefficient_derivative(:, lm, spin) = coefficient_derivative(:, lm, spin) + &
                                                              sign*kinetic_value(lm)*grid%kinetic_radial(ir, :, l, side)
                     END DO
                  END IF
               END DO
            END DO
         END DO
      END DO
      potential = 0.0_dp
      DO spin = 1, nspin
         DO q = 1, onecenter%nproj
            DO p = 1, onecenter%nproj
               pair = onecenter%pair_index(p, q)
               potential(p, q, spin) = SUM(onecenter%gaunt(p, q, :)*coefficient_derivative(pair, :, spin))
            END DO
         END DO
      END DO
      CALL xc_dset_release(derivatives)
      CALL xc_rho_set_release(rho_set)
      DEALLOCATE (vxc, vxg, vtau)
   END SUBROUTINE paw_xc_evaluate

END MODULE paw_xc
