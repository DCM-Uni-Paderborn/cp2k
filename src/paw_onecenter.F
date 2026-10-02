!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief PAW valence Hartree and exchange corrections on the original partial-wave radial mesh.
!>        The pseudo one-center density includes the UPF compensation functions.
!>        Core densities belong to XC, while fixed kinetic and ionic terms are in PP_DIJ.
! **************************************************************************************************
MODULE paw_onecenter
   USE atom_upf,                        ONLY: atom_upfpot_type
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: fourpi
   USE paw_core_kinetic,                ONLY: paw_core_kinetic_build_ae
   USE paw_types,                       ONLY: paw_onecenter_type
   USE splines,                         ONLY: spline3ders
   USE upf_projector_integrals,         ONLY: upf_projector_angular_init,&
                                              upf_projector_product_gaunt,&
                                              upf_projector_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'paw_onecenter'

   PUBLIC :: paw_onecenter_type, paw_onecenter_init, paw_onecenter_hartree, &
             paw_onecenter_exchange, paw_onecenter_core_exchange_init, paw_radial_coulomb

CONTAINS

! **************************************************************************************************
!> \brief Cache the AE-minus-compensated-PS Coulomb matrices for each angular momentum.
!> \param pot Scalar PAW dataset including original r*phi partial waves and r**2*Q.
!> \param projectors Native radial/magnetic projector mapping and angular integrals.
!> \param onecenter Immutable data, independent of the orbital basis and METHOD.
! **************************************************************************************************
   SUBROUTINE paw_onecenter_init(pot, projectors, onecenter)
      TYPE(atom_upfpot_type), INTENT(IN)                 :: pot
      TYPE(upf_projector_type), INTENT(IN)               :: projectors
      TYPE(paw_onecenter_type), INTENT(OUT)              :: onecenter

      INTEGER                                            :: i, j, l, n, nb, nc, p, pair, q
      REAL(dp), ALLOCATABLE                              :: ae(:, :), ae_coulomb(:, :), ps(:, :), &
                                                            ps_coulomb(:, :)

      CPASSERT(pot%is_paw .AND. .NOT. pot%has_so)
      CPASSERT(ALLOCATED(pot%paw) .AND. ALLOCATED(pot%qijl))
      CPASSERT(ALL(pot%lbeta == projectors%lbeta))
      n = pot%paw%cutoff_index
      CPASSERT(n >= 4 .AND. n <= SIZE(pot%r))
      nb = pot%number_of_proj
      nc = nb*(nb + 1)/2
      onecenter%nproj = projectors%nproj
      onecenter%npair = nc
      onecenter%lmax = 2*MAXVAL(pot%lbeta)
      CALL upf_projector_product_gaunt(projectors, onecenter%gaunt)
      ALLOCATE (onecenter%pair_index(projectors%nproj, projectors%nproj))
      DO q = 1, projectors%nproj
         j = projectors%radial_index(q)
         DO p = 1, projectors%nproj
            i = projectors%radial_index(p)
            onecenter%pair_index(p, q) = MAX(i, j)*(MAX(i, j) - 1)/2 + MIN(i, j)
         END DO
      END DO
      ALLOCATE (onecenter%hartree(nc, nc, 0:onecenter%lmax), ae(n, nc), ps(n, nc), &
                ae_coulomb(nc, nc), ps_coulomb(nc, nc))
      DO l = 0, onecenter%lmax
         ae = 0.0_dp
         ps = 0.0_dp
         DO j = 1, nb
            DO i = 1, j
               IF (l < ABS(pot%lbeta(i) - pot%lbeta(j)) .OR. l > pot%lbeta(i) + pot%lbeta(j)) CYCLE
               IF (MODULO(l + pot%lbeta(i) + pot%lbeta(j), 2) /= 0) CYCLE
               pair = j*(j - 1)/2 + i
               ! Do not use unbound partial waves beyond the supplied augmentation sphere.
               ae(:, pair) = pot%paw%ae_wfc(:n, i)*pot%paw%ae_wfc(:n, j)
               ps(:, pair) = pot%paw%ps_wfc(:n, i)*pot%paw%ps_wfc(:n, j) + pot%qijl(:n, pair, l)
            END DO
         END DO
         CALL paw_radial_coulomb(pot%r(:n), ae, l, ae_coulomb)
         CALL paw_radial_coulomb(pot%r(:n), ps, l, ps_coulomb)
         onecenter%hartree(:, :, l) = ae_coulomb - ps_coulomb
      END DO
   END SUBROUTINE paw_onecenter_init

! **************************************************************************************************
!> \brief Original-data core-valence operator and core-core exchange reference for a closed PAW core.
!>        These are AE contributions: smooth RI exchange contains only valence transitions.
!>        Both use the complete supplied mesh, without truncating core tails at the augmentation sphere.
!>        The usual PAW core-locality and partial-wave completeness assumptions remain applicable.
!> \param pot Original scalar PAW data; neither core nor valence waves are normalized or fitted.
!> \param projectors Original valence magnetic ordering.
!> \param onecenter Prepared valence cache, with the lazy unit-fraction core exchange cache overwritten.
! **************************************************************************************************
   SUBROUTINE paw_onecenter_core_exchange_init(pot, projectors, onecenter)
      TYPE(atom_upfpot_type), INTENT(IN)                 :: pot
      TYPE(upf_projector_type), INTENT(IN)               :: projectors
      TYPE(paw_onecenter_type), INTENT(INOUT)            :: onecenter

      INTEGER                                            :: c, d, first_c, first_d, i, j, l, lc, ld, &
                                                            lm, nb, nr, nshell, p, q, shell, t
      INTEGER, ALLOCATABLE                               :: angular(:), selected(:), shells(:)
      LOGICAL                                            :: available
      REAL(dp)                                           :: angular_factor, core_fraction, &
                                                            multiplicity, scalar_coulomb(1, 1)
      REAL(dp), ALLOCATABLE                              :: coulomb(:, :), gaunt(:, :, :), &
                                                            products(:, :), scalar_product(:, :), &
                                                            tau_ae(:)
      TYPE(upf_projector_type)                           :: combined

      CPASSERT(onecenter%nproj == projectors%nproj)
      onecenter%core_exchange_initialized = .TRUE.
      onecenter%core_exchange_available = .FALSE.
      onecenter%core_exchange_energy = 0.0_dp
      IF (ALLOCATED(onecenter%core_exchange_matrix)) DEALLOCATE (onecenter%core_exchange_matrix)
      CALL paw_core_kinetic_build_ae(pot, tau_ae, available, onecenter%core_exchange_reason, selected)
      IF (.NOT. available) RETURN
      DEALLOCATE (tau_ae)
      ALLOCATE (onecenter%core_exchange_matrix(onecenter%nproj, onecenter%nproj), SOURCE=0.0_dp)
      nshell = COUNT(selected > 0)
      IF (nshell == 0) THEN
         onecenter%core_exchange_available = .TRUE.
         onecenter%core_exchange_reason = ''
         RETURN
      END IF
      nb = pot%number_of_proj
      nr = SIZE(pot%r)
      CPASSERT(ALL(pot%lbeta == projectors%lbeta))
      onecenter%core_exchange_reason = 'PAW core exchange requires complete original AE partial waves'
      IF (.NOT. ALLOCATED(pot%paw%ae_wfc)) RETURN
      IF (.NOT. ALL(SHAPE(pot%paw%ae_wfc) == [nr, nb])) RETURN
      IF (.NOT. ALL(ieee_is_finite(pot%paw%ae_wfc))) RETURN
      ALLOCATE (shells(nshell), angular(nb + nshell))
      angular(:nb) = pot%lbeta
      shell = 0
      DO i = 1, SIZE(selected)
         IF (selected(i) == 0) CYCLE
         shell = shell + 1
         shells(shell) = i
         angular(nb + shell) = pot%paw%gipaw_core_orbitals%angular_momenta(i)
      END DO
      CALL upf_projector_angular_init(angular, MAXVAL(angular), combined)
      CALL upf_projector_product_gaunt(combined, gaunt)
      CPASSERT(ALL(combined%radial_index(:onecenter%nproj) == projectors%radial_index))
      CPASSERT(ALL(combined%angular_index(:onecenter%nproj) == projectors%angular_index))
      ALLOCATE (products(nr, nb), coulomb(nb, nb), scalar_product(nr, 1))
      first_c = onecenter%nproj + 1
      DO shell = 1, nshell
         i = shells(shell)
         lc = angular(nb + shell)
         core_fraction = REAL(selected(i), dp)/REAL(2*(2*lc + 1), dp)
         ! Each closed-core magnetic state has this occupation in either spin.
         DO l = 0, MAXVAL(pot%lbeta) + lc
            products = 0.0_dp
            DO j = 1, nb
               IF (l < ABS(pot%lbeta(j) - lc) .OR. l > pot%lbeta(j) + lc) CYCLE
               IF (MODULO(l + pot%lbeta(j) + lc, 2) /= 0) CYCLE
               products(:, j) = pot%paw%ae_wfc(:, j)*pot%paw%gipaw_core_orbitals%wfc(:, i)
            END DO
            IF (.NOT. ANY(products /= 0.0_dp)) CYCLE
            CALL paw_radial_coulomb(pot%r, products, l, coulomb)
            DO q = 1, onecenter%nproj
               DO p = 1, q
                  angular_factor = 0.0_dp
                  DO c = first_c, first_c + 2*lc
                     DO lm = l*l + 1, (l + 1)**2
                        angular_factor = angular_factor + gaunt(p, c, lm)*gaunt(q, c, lm)
                     END DO
                  END DO
                  onecenter%core_exchange_matrix(p, q) = onecenter%core_exchange_matrix(p, q) - &
                                        core_fraction*angular_factor*coulomb(projectors%radial_index(p), projectors%radial_index(q))
               END DO
            END DO
         END DO
         first_d = onecenter%nproj + 1
         DO t = 1, shell
            j = shells(t)
            ld = angular(nb + t)
            ! Two physical spins cancel the exchange energy factor 1/2.
            ! The off-diagonal shell pair occurs in both ordered summations.
            multiplicity = MERGE(1.0_dp, 2.0_dp, t == shell)*core_fraction* &
                           REAL(selected(j), dp)/REAL(2*(2*ld + 1), dp)
            scalar_product(:, 1) = pot%paw%gipaw_core_orbitals%wfc(:, i)* &
                                   pot%paw%gipaw_core_orbitals%wfc(:, j)
            DO l = ABS(lc - ld), lc + ld, 2
               CALL paw_radial_coulomb(pot%r, scalar_product, l, scalar_coulomb)
               angular_factor = 0.0_dp
               DO d = first_d, first_d + 2*ld
                  DO c = first_c, first_c + 2*lc
                     DO lm = l*l + 1, (l + 1)**2
                        angular_factor = angular_factor + gaunt(c, d, lm)**2
                     END DO
                  END DO
               END DO
               onecenter%core_exchange_energy = onecenter%core_exchange_energy - &
                                                multiplicity*angular_factor*scalar_coulomb(1, 1)
            END DO
            first_d = first_d + 2*ld + 1
         END DO
         first_c = first_c + 2*lc + 1
      END DO
      DO q = 1, onecenter%nproj
         DO p = q + 1, onecenter%nproj
            onecenter%core_exchange_matrix(p, q) = onecenter%core_exchange_matrix(q, p)
         END DO
      END DO
      onecenter%core_exchange_reason = 'Nonfinite original-data PAW core exchange integral'
      IF (.NOT. ALL(ieee_is_finite(onecenter%core_exchange_matrix))) RETURN
      IF (.NOT. ieee_is_finite(onecenter%core_exchange_energy)) RETURN
      onecenter%core_exchange_available = .TRUE.
      onecenter%core_exchange_reason = ''
   END SUBROUTINE paw_onecenter_core_exchange_init

! **************************************************************************************************
!> \brief Evaluate the valence Hartree correction and its derivative in projector space.
!> \param onecenter Prepared PAW data.
!> \param density Real symmetric projector density, summed over spins and k points.
!> \param energy One-center Hartree AE-minus-PS energy in Hartree.
!> \param potential Derivative dE/d(density), with a full matrix contraction convention.
! **************************************************************************************************
   SUBROUTINE paw_onecenter_hartree(onecenter, density, energy, potential)
      TYPE(paw_onecenter_type), INTENT(IN)               :: onecenter
      REAL(dp), INTENT(IN)                               :: density(:, :)
      REAL(dp), INTENT(OUT)                              :: energy, potential(:, :)

      INTEGER                                            :: l, lm, p, pair, q
      REAL(dp)                                           :: derivative(onecenter%npair), &
                                                            coefficients(onecenter%npair)

      CPASSERT(ALL(SHAPE(density) == onecenter%nproj))
      CPASSERT(ALL(SHAPE(potential) == onecenter%nproj))
      energy = 0.0_dp
      potential = 0.0_dp
      DO lm = 1, SIZE(onecenter%gaunt, 3)
         l = INT(SQRT(REAL(lm - 1, dp)))
         coefficients = 0.0_dp
         DO q = 1, onecenter%nproj
            DO p = 1, onecenter%nproj
               pair = onecenter%pair_index(p, q)
               coefficients(pair) = coefficients(pair) + onecenter%gaunt(p, q, lm)*density(p, q)
            END DO
         END DO
         derivative = MATMUL(onecenter%hartree(:, :, l), coefficients)
         energy = energy + 0.5_dp*SUM(coefficients*derivative)
         DO q = 1, onecenter%nproj
            DO p = 1, onecenter%nproj
               pair = onecenter%pair_index(p, q)
               potential(p, q) = potential(p, q) + onecenter%gaunt(p, q, lm)*derivative(pair)
            END DO
         END DO
      END DO
   END SUBROUTINE paw_onecenter_hartree

! **************************************************************************************************
!> \brief Unscreened scalar PAW valence exact-exchange correction and projector derivative.
!>        Delta J_ij,kl is the AE-minus-compensated-PS pair Coulomb integral. Exchange
!>        crosses the density indices: E=-fraction/2 sum_sigma,ijkl D_ik D_jl Delta J_ij,kl.
!>        An initialized original-data core cache supplies linear core-valence and constant core-core terms.
!> \param onecenter Prepared radial Coulomb and magnetic Gaunt data.
!> \param density Real symmetric, physically spin-resolved projector densities.
!>        A spin-summed restricted density must be supplied as two identical half-density blocks.
!> \param fraction Exact-exchange fraction.
!> \param energy One-center exchange correction including an available initialized closed core, in Hartree.
!> \param potential Full-matrix derivative V_pq=-fraction sum_jk D_jk Delta J_pj,qk per spin.
! **************************************************************************************************
   SUBROUTINE paw_onecenter_exchange(onecenter, density, fraction, energy, potential)
      TYPE(paw_onecenter_type), INTENT(IN)               :: onecenter
      REAL(dp), INTENT(IN)                               :: density(:, :, :), fraction
      REAL(dp), INTENT(OUT)                              :: energy, potential(:, :, :)

      INTEGER                                            :: j, k, l, lm, p, pair_pj, pair_qk, q, spin
      LOGICAL                                            :: symmetric
      REAL(dp)                                           :: exchange, gaunt_pj, gaunt_qk

      CPASSERT(SIZE(density, 1) == onecenter%nproj .AND. SIZE(density, 2) == onecenter%nproj)
      CPASSERT(SIZE(density, 3) == 1 .OR. SIZE(density, 3) == 2)
      CPASSERT(ALL(SHAPE(potential) == SHAPE(density)))
      CPASSERT(ieee_is_finite(fraction) .AND. ALL(ieee_is_finite(density)))
      DO spin = 1, SIZE(density, 3)
         symmetric = MAXVAL(ABS(density(:, :, spin) - TRANSPOSE(density(:, :, spin)))) <= &
            64.0_dp*EPSILON(1.0_dp)*MAX(1.0_dp, MAXVAL(ABS(density(:, :, spin))))
         CPASSERT(symmetric)
      END DO
      potential(:, :, :) = 0.0_dp
      DO spin = 1, SIZE(density, 3)
         DO lm = 1, SIZE(onecenter%gaunt, 3)
            l = INT(SQRT(REAL(lm - 1, dp)))
            DO q = 1, onecenter%nproj
               DO p = 1, q
                  exchange = 0.0_dp
                  DO j = 1, onecenter%nproj
                     gaunt_pj = onecenter%gaunt(p, j, lm)
                     IF (gaunt_pj == 0.0_dp) CYCLE
                     pair_pj = onecenter%pair_index(p, j)
                     DO k = 1, onecenter%nproj
                        gaunt_qk = onecenter%gaunt(q, k, lm)
                        IF (gaunt_qk == 0.0_dp) CYCLE
                        pair_qk = onecenter%pair_index(q, k)
                        exchange = exchange + density(j, k, spin)*gaunt_pj*gaunt_qk* &
                                   onecenter%hartree(pair_pj, pair_qk, l)
                     END DO
                  END DO
                  potential(p, q, spin) = potential(p, q, spin) - fraction*exchange
               END DO
            END DO
         END DO
         DO q = 1, onecenter%nproj
            DO p = q + 1, onecenter%nproj
               potential(p, q, spin) = potential(q, p, spin)
            END DO
         END DO
      END DO
      energy = 0.5_dp*SUM(density*potential)
      IF (onecenter%core_exchange_initialized) THEN
         CPASSERT(onecenter%core_exchange_available)
         DO spin = 1, SIZE(density, 3)
            ! The core-valence contribution is linear, without the valence-valence factor 1/2.
            energy = energy + fraction*SUM(density(:, :, spin)*onecenter%core_exchange_matrix)
            potential(:, :, spin) = potential(:, :, spin) + fraction*onecenter%core_exchange_matrix
         END DO
         energy = energy + fraction*onecenter%core_exchange_energy
      END IF
   END SUBROUTINE paw_onecenter_exchange

! **************************************************************************************************
!> \brief Coulomb bilinear form for tabulated r**2*rho_LM radial channels.
!>        Integrate cubic interpolants over the two triangles r<s and s<r.
!>        This avoids applying a tensor-product radial rule across the kernel's diagonal cusp.
!>        Eight-point Gauss rules also integrate the within-interval cumulative moments.
!> \param r Increasing nonnegative radial coordinates. Functions vanish outside this interval.
!> \param products Radial density channels, including r**2 but no quadrature weight.
!> \param l Angular momentum.
!> \param coulomb Symmetric matrix including 4*pi/(2*l+1), without the energy factor 1/2.
! **************************************************************************************************
   SUBROUTINE paw_radial_coulomb(r, products, l, coulomb)
      REAL(dp), INTENT(IN)                               :: r(:), products(:, :)
      INTEGER, INTENT(IN)                                :: l
      REAL(dp), INTENT(OUT)                              :: coulomb(:, :)

      REAL(dp), PARAMETER :: nodes(8) = [0.019855071751231884_dp, 0.10166676129318663_dp, &
         0.23723379504183550_dp, 0.40828267875217510_dp, 0.59171732124782490_dp, &
         0.76276620495816450_dp, 0.89833323870681337_dp, 0.98014492824876812_dp], weights(8) = [ &
         0.05061426814518813_dp, 0.11119051722668724_dp, 0.15685332293894364_dp, &
         0.18134189168918099_dp, 0.18134189168918099_dp, 0.15685332293894364_dp, &
         0.11119051722668724_dp, 0.05061426814518813_dp]

      INTEGER                                            :: a, b, i, j, k, n, nc
      REAL(dp)                                           :: factor, h, radius, t, u
      REAL(dp), ALLOCATABLE                              :: c(:, :), cumulative(:), inner(:), &
                                                            panel(:), sample(:), slopes(:, :), &
                                                            value(:)

      n = SIZE(r)
      nc = SIZE(products, 2)
      CPASSERT(n >= 4 .AND. SIZE(products, 1) == n .AND. l >= 0)
      CPASSERT(ALL(SHAPE(coulomb) == nc))
      CPASSERT(r(1) >= 0.0_dp .AND. ALL(r(2:) > r(:n - 1)))
      ALLOCATE (slopes(n, nc), c(0:3, nc), cumulative(nc), panel(nc), inner(nc), value(nc), sample(nc))
      DO a = 1, nc
         CALL spline3ders(r, products(:, a), r, dynew=slopes(:, a))
      END DO
      coulomb = 0.0_dp
      cumulative = 0.0_dp
      DO i = 1, n - 1
         h = r(i + 1) - r(i)
         c(0, :) = products(i, :)
         c(1, :) = h*slopes(i, :)
         c(2, :) = 3.0_dp*(products(i + 1, :) - products(i, :)) - h*(2.0_dp*slopes(i, :) + slopes(i + 1, :))
         c(3, :) = 2.0_dp*(products(i, :) - products(i + 1, :)) + h*(slopes(i, :) + slopes(i + 1, :))
         panel = 0.0_dp
         DO j = 1, SIZE(nodes)
            t = nodes(j)
            radius = r(i) + h*t
            value(:) = c(0, :) + t*(c(1, :) + t*(c(2, :) + t*c(3, :)))
            ! The source moment up to this outer node includes a partial current interval.
            inner = 0.0_dp
            DO k = 1, SIZE(nodes)
               u = t*nodes(k)
               sample(:) = c(0, :) + u*(c(1, :) + u*(c(2, :) + u*c(3, :)))
               inner(:) = inner + weights(k)*(r(i) + h*u)**l*sample
            END DO
            inner(:) = cumulative + h*t*inner
            factor = h*weights(j)/radius**(l + 1)
            DO b = 1, nc
               DO a = 1, b
                  coulomb(a, b) = coulomb(a, b) + factor*(value(a)*inner(b) + inner(a)*value(b))
               END DO
            END DO
            panel(:) = panel + h*weights(j)*radius**l*value
         END DO
         cumulative(:) = cumulative + panel
      END DO
      coulomb = coulomb*fourpi/REAL(2*l + 1, dp)
      DO b = 1, nc
         DO a = b + 1, nc
            coulomb(a, b) = coulomb(b, a)
         END DO
      END DO
   END SUBROUTINE paw_radial_coulomb

END MODULE paw_onecenter
