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
   USE cp_array_sort,                   ONLY: cp_1d_r_sort
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: fourpi,&
                                              pi
   USE paw_core_kinetic,                ONLY: paw_core_kinetic_build_ae
   USE paw_types,                       ONLY: paw_core_exchange_type,&
                                              paw_onecenter_type
   USE splines,                         ONLY: spline3ders
   USE upf_projector_integrals,         ONLY: upf_projector_angular_init,&
                                              upf_projector_product_gaunt,&
                                              upf_projector_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'paw_onecenter'

   PUBLIC :: paw_onecenter_type, paw_onecenter_init, paw_onecenter_hartree, &
             paw_onecenter_exchange, paw_onecenter_core_exchange_init, paw_onecenter_core_exchange_prepare, &
             paw_onecenter_core_exchange, paw_radial_coulomb, paw_radial_truncated, paw_radial_short

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

      TYPE(paw_core_exchange_type)                       :: core

      CALL paw_core_exchange_build(pot, projectors, onecenter%nproj, core)
      onecenter%core_exchange_initialized = core%initialized
      onecenter%core_exchange_available = core%available
      onecenter%core_exchange_reason = core%reason
      onecenter%core_exchange_energy = core%energy
      CALL MOVE_ALLOC(core%matrix, onecenter%core_exchange_matrix)
   END SUBROUTINE paw_onecenter_core_exchange_init

! **************************************************************************************************
!> \brief Prepare or reuse original-data closed-core exchange for a specific truncated interaction.
!> \param pot Original scalar PAW data with complete core orbital tails.
!> \param projectors Original valence magnetic ordering.
!> \param onecenter Kind cache; the unscreened Coulomb data are retained separately.
!> \param cutoff_radius Fixed physical interaction radius, in bohr.
!> \param cache_index Index of the prepared interaction cache, including unavailable-data diagnostics.
!> \param omega Optional inverse erfc screening length instead of a cutoff_radius.
! **************************************************************************************************
   SUBROUTINE paw_onecenter_core_exchange_prepare(pot, projectors, onecenter, cutoff_radius, cache_index, omega)
      TYPE(atom_upfpot_type), INTENT(IN)                 :: pot
      TYPE(upf_projector_type), INTENT(IN)               :: projectors
      TYPE(paw_onecenter_type), INTENT(INOUT)            :: onecenter
      REAL(dp), INTENT(IN), OPTIONAL                     :: cutoff_radius
      INTEGER, INTENT(OUT)                               :: cache_index
      REAL(dp), INTENT(IN), OPTIONAL                     :: omega

      INTEGER                                            :: i, n
      TYPE(paw_core_exchange_type), ALLOCATABLE          :: extended(:)

      CPASSERT(PRESENT(cutoff_radius) .NEQV. PRESENT(omega))
      IF (PRESENT(cutoff_radius)) THEN
         CPASSERT(ieee_is_finite(cutoff_radius) .AND. cutoff_radius >= 0.0_dp)
      END IF
      IF (PRESENT(omega)) THEN
         CPASSERT(ieee_is_finite(omega) .AND. omega >= 0.0_dp)
      END IF
      n = 0
      IF (ALLOCATED(onecenter%interaction_core_exchange)) THEN
         n = SIZE(onecenter%interaction_core_exchange)
         DO i = 1, n
            IF (PRESENT(cutoff_radius)) THEN
               IF (onecenter%interaction_core_exchange(i)%cutoff_radius /= cutoff_radius) CYCLE
               IF (onecenter%interaction_core_exchange(i)%omega /= -1.0_dp) CYCLE
            ELSE
               IF (onecenter%interaction_core_exchange(i)%omega /= omega) CYCLE
               IF (onecenter%interaction_core_exchange(i)%cutoff_radius /= -1.0_dp) CYCLE
            END IF
            cache_index = i
            RETURN
         END DO
      END IF
      ALLOCATE (extended(n + 1))
      IF (n > 0) extended(:n) = onecenter%interaction_core_exchange
      CALL paw_core_exchange_build(pot, projectors, onecenter%nproj, extended(n + 1), cutoff_radius, omega)
      CALL MOVE_ALLOC(extended, onecenter%interaction_core_exchange)
      cache_index = n + 1
   END SUBROUTINE paw_onecenter_core_exchange_prepare

! **************************************************************************************************
!> \brief Build a unit-fraction AE core-valence matrix and core-core energy from original orbitals.
!> \param pot Original scalar PAW data, without fitted orbitals or shortened tails.
!> \param projectors Original valence magnetic ordering.
!> \param nproj Number of valence magnetic projectors.
!> \param core Prepared data or an explicit unavailable-data reason.
!> \param cutoff_radius Optional truncated interaction; absence selects Coulomb.
!> \param omega ...
! **************************************************************************************************
   SUBROUTINE paw_core_exchange_build(pot, projectors, nproj, core, cutoff_radius, omega)
      TYPE(atom_upfpot_type), INTENT(IN)                 :: pot
      TYPE(upf_projector_type), INTENT(IN)               :: projectors
      INTEGER, INTENT(IN)                                :: nproj
      TYPE(paw_core_exchange_type), INTENT(OUT)          :: core
      REAL(dp), INTENT(IN), OPTIONAL                     :: cutoff_radius, omega

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

      CPASSERT(nproj == projectors%nproj)
      core%initialized = .TRUE.
      core%available = .FALSE.
      core%energy = 0.0_dp
      CPASSERT(.NOT. (PRESENT(cutoff_radius) .AND. PRESENT(omega)))
      IF (PRESENT(cutoff_radius)) core%cutoff_radius = cutoff_radius
      IF (PRESENT(omega)) core%omega = omega
      CALL paw_core_kinetic_build_ae(pot, tau_ae, available, core%reason, selected)
      IF (.NOT. available) RETURN
      DEALLOCATE (tau_ae)
      ALLOCATE (core%matrix(nproj, nproj), SOURCE=0.0_dp)
      nshell = COUNT(selected > 0)
      IF (nshell == 0) THEN
         core%available = .TRUE.
         core%reason = ''
         RETURN
      END IF
      nb = pot%number_of_proj
      nr = SIZE(pot%r)
      CPASSERT(ALL(pot%lbeta == projectors%lbeta))
      core%reason = 'PAW core exchange requires complete original AE partial waves'
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
      CPASSERT(ALL(combined%radial_index(:nproj) == projectors%radial_index))
      CPASSERT(ALL(combined%angular_index(:nproj) == projectors%angular_index))
      ALLOCATE (products(nr, nb), coulomb(nb, nb), scalar_product(nr, 1))
      first_c = nproj + 1
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
            CALL paw_core_radial_integral(pot%r, products, l, coulomb, cutoff_radius, omega)
            DO q = 1, nproj
               DO p = 1, q
                  angular_factor = 0.0_dp
                  DO c = first_c, first_c + 2*lc
                     DO lm = l*l + 1, (l + 1)**2
                        angular_factor = angular_factor + gaunt(p, c, lm)*gaunt(q, c, lm)
                     END DO
                  END DO
                  core%matrix(p, q) = core%matrix(p, q) - &
                                      core_fraction*angular_factor*coulomb(projectors%radial_index(p), projectors%radial_index(q))
               END DO
            END DO
         END DO
         first_d = nproj + 1
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
               CALL paw_core_radial_integral(pot%r, scalar_product, l, scalar_coulomb, cutoff_radius, omega)
               angular_factor = 0.0_dp
               DO d = first_d, first_d + 2*ld
                  DO c = first_c, first_c + 2*lc
                     DO lm = l*l + 1, (l + 1)**2
                        angular_factor = angular_factor + gaunt(c, d, lm)**2
                     END DO
                  END DO
               END DO
               core%energy = core%energy - &
                             multiplicity*angular_factor*scalar_coulomb(1, 1)
            END DO
            first_d = first_d + 2*ld + 1
         END DO
         first_c = first_c + 2*lc + 1
      END DO
      DO q = 1, nproj
         DO p = q + 1, nproj
            core%matrix(p, q) = core%matrix(q, p)
         END DO
      END DO
      core%reason = 'Nonfinite original-data PAW core exchange integral'
      IF (.NOT. ALL(ieee_is_finite(core%matrix))) RETURN
      IF (.NOT. ieee_is_finite(core%energy)) RETURN
      core%available = .TRUE.
      core%reason = ''
   END SUBROUTINE paw_core_exchange_build

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
!> \brief Core-only PAW exchange for a globally reconstructed valence RI contraction.
!>        The valence exchange is already present in RI and must not be added again locally.
!> \param onecenter Prepared original-data core cache.
!> \param density Physical spin projector densities, including occupations and k-point weights.
!> \param fraction Exchange fraction for this interaction replica.
!> \param cache_index Positive truncated-cache index, or zero for the separate Coulomb cache.
!> \param energy Linear valence-core energy plus the core-core constant, counted once.
!> \param potential Full-matrix spin derivative, without a valence-valence factor 1/2.
! **************************************************************************************************
   SUBROUTINE paw_onecenter_core_exchange(onecenter, density, fraction, cache_index, energy, potential)
      TYPE(paw_onecenter_type), INTENT(IN)               :: onecenter
      REAL(dp), INTENT(IN)                               :: density(:, :, :), fraction
      INTEGER, INTENT(IN)                                :: cache_index
      REAL(dp), INTENT(OUT)                              :: energy, potential(:, :, :)

      INTEGER                                            :: spin

      CPASSERT(SIZE(density, 1) == onecenter%nproj .AND. SIZE(density, 2) == onecenter%nproj)
      CPASSERT(SIZE(density, 3) == 1 .OR. SIZE(density, 3) == 2)
      CPASSERT(ALL(SHAPE(potential) == SHAPE(density)))
      CPASSERT(ieee_is_finite(fraction) .AND. ALL(ieee_is_finite(density)))
      IF (cache_index == 0) THEN
         CPASSERT(onecenter%core_exchange_initialized .AND. onecenter%core_exchange_available)
         energy = fraction*onecenter%core_exchange_energy
         DO spin = 1, SIZE(density, 3)
            potential(:, :, spin) = fraction*onecenter%core_exchange_matrix
         END DO
      ELSE
         CPASSERT(ALLOCATED(onecenter%interaction_core_exchange))
         CPASSERT(cache_index > 0 .AND. cache_index <= SIZE(onecenter%interaction_core_exchange))
         ASSOCIATE (core => onecenter%interaction_core_exchange(cache_index))
            CPASSERT(core%initialized .AND. core%available)
            energy = fraction*core%energy
            DO spin = 1, SIZE(density, 3)
               potential(:, :, spin) = fraction*core%matrix
            END DO
         END ASSOCIATE
      END IF
      energy = energy + SUM(density*potential)
   END SUBROUTINE paw_onecenter_core_exchange

! **************************************************************************************************
!> \brief Select the original-data radial core interaction.
!> \param r Complete supplied radial mesh.
!> \param products Original orbital products, including r**2.
!> \param l Angular momentum.
!> \param integral Symmetric bilinear matrix.
!> \param cutoff_radius Optional truncated interaction; absence selects Coulomb.
!> \param omega ...
! **************************************************************************************************
   SUBROUTINE paw_core_radial_integral(r, products, l, integral, cutoff_radius, omega)
      REAL(dp), INTENT(IN)                               :: r(:), products(:, :)
      INTEGER, INTENT(IN)                                :: l
      REAL(dp), INTENT(OUT)                              :: integral(:, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: cutoff_radius, omega

      CPASSERT(.NOT. (PRESENT(cutoff_radius) .AND. PRESENT(omega)))
      IF (PRESENT(omega)) THEN
         CALL paw_radial_short(r, products, l, omega, integral)
      ELSE IF (PRESENT(cutoff_radius)) THEN
         CALL paw_radial_truncated(r, products, l, cutoff_radius, integral)
      ELSE
         CALL paw_radial_coulomb(r, products, l, integral)
      END IF
   END SUBROUTINE paw_core_radial_integral

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

! **************************************************************************************************
!> \brief Radial bilinear form for Theta(Rc-|r-s|)/|r-s| on the complete original mesh.
!>        The angular coefficient is integrated in the interpoint distance u. Its integrand
!>        is a polynomial of degree 2*l, so l+1 Gauss nodes suffice without a fixed angular cap.
!>        Radial triangles and the surfaces r+s=Rc and |r-s|=Rc are resolved separately.
!>        Cubic interpolation knots crossed by these surfaces also split the outer panels.
!> \param r Increasing nonnegative coordinates; channels vanish outside the supplied mesh.
!> \param products Radial channels including r**2, without quadrature weights.
!> \param l Angular momentum.
!> \param cutoff_radius Fixed physical radius in bohr, with zero selecting a zero interaction.
!> \param integral Symmetric matrix including 4*pi/(2*l+1), without the energy factor 1/2.
! **************************************************************************************************
   SUBROUTINE paw_radial_truncated(r, products, l, cutoff_radius, integral)
      REAL(dp), INTENT(IN)                               :: r(:), products(:, :)
      INTEGER, INTENT(IN)                                :: l
      REAL(dp), INTENT(IN)                               :: cutoff_radius
      REAL(dp), INTENT(OUT)                              :: integral(:, :)

      REAL(dp), PARAMETER :: nodes(8) = [0.019855071751231884_dp, 0.10166676129318663_dp, &
         0.23723379504183550_dp, 0.40828267875217510_dp, 0.59171732124782490_dp, &
         0.76276620495816450_dp, 0.89833323870681337_dp, 0.98014492824876812_dp], weights(8) = [ &
         0.05061426814518813_dp, 0.11119051722668724_dp, 0.15685332293894364_dp, &
         0.18134189168918099_dp, 0.18134189168918099_dp, 0.15685332293894364_dp, &
         0.11119051722668724_dp, 0.05061426814518813_dp]

      INTEGER                                            :: a, b, i, j, k, n, nc, panel, source_panel
      INTEGER, ALLOCATABLE                               :: order(:)
      REAL(dp)                                           :: factor, full_upper, h, kernel, left, &
                                                            partial_lower, radius, right, &
                                                            source_h, source_left, source_r, &
                                                            source_right, t
      REAL(dp), ALLOCATABLE :: angular_nodes(:), angular_weights(:), c(:, :, :), cumulative(:, :), &
         inner(:), knots(:), sample(:), slopes(:, :), value(:)

      n = SIZE(r)
      nc = SIZE(products, 2)
      CPASSERT(n >= 4 .AND. SIZE(products, 1) == n .AND. l >= 0)
      CPASSERT(ALL(SHAPE(integral) == nc))
      CPASSERT(r(1) >= 0.0_dp .AND. ALL(r(2:) > r(:n - 1)))
      CPASSERT(ieee_is_finite(cutoff_radius) .AND. cutoff_radius >= 0.0_dp)
      integral = 0.0_dp
      IF (cutoff_radius == 0.0_dp) RETURN
      IF (cutoff_radius >= 2.0_dp*r(n)) THEN
         CALL paw_radial_coulomb(r, products, l, integral)
         RETURN
      END IF
      ALLOCATE (slopes(n, nc), c(0:3, nc, n - 1), cumulative(nc, n), inner(nc), sample(nc), value(nc))
      ALLOCATE (knots(3*n + 2), order(3*n + 2), angular_nodes(l + 1), angular_weights(l + 1))
      CALL paw_angular_gauss(angular_nodes, angular_weights)
      DO a = 1, nc
         CALL spline3ders(r, products(:, a), r, dynew=slopes(:, a))
      END DO
      cumulative(:, 1) = 0.0_dp
      DO i = 1, n - 1
         h = r(i + 1) - r(i)
         c(0, :, i) = products(i, :)
         c(1, :, i) = h*slopes(i, :)
         c(2, :, i) = 3.0_dp*(products(i + 1, :) - products(i, :)) - h*(2.0_dp*slopes(i, :) + slopes(i + 1, :))
         c(3, :, i) = 2.0_dp*(products(i, :) - products(i + 1, :)) + h*(slopes(i, :) + slopes(i + 1, :))
         cumulative(:, i + 1) = cumulative(:, i)
         DO k = 1, SIZE(nodes)
            t = nodes(k)
            sample(:) = c(0, :, i) + t*(c(1, :, i) + t*(c(2, :, i) + t*c(3, :, i)))
            cumulative(:, i + 1) = cumulative(:, i + 1) + h*weights(k)*(r(i) + h*t)**l*sample
         END DO
      END DO
      knots(:n) = r
      knots(n + 1:2*n) = r + cutoff_radius
      knots(2*n + 1:3*n) = ABS(r - cutoff_radius)
      knots(3*n + 1:) = [0.5_dp*cutoff_radius, cutoff_radius]
      CALL cp_1d_r_sort(knots, SIZE(knots), order)
      DO panel = 1, SIZE(knots) - 1
         left = knots(panel)
         right = knots(panel + 1)
         IF (right <= left .OR. left < r(1) .OR. right > r(n)) CYCLE
         i = paw_radial_interval(r, left)
         h = r(i + 1) - r(i)
         DO j = 1, SIZE(nodes)
            radius = left + (right - left)*nodes(j)
            t = (radius - r(i))/h
            value(:) = c(0, :, i) + t*(c(1, :, i) + t*(c(2, :, i) + t*c(3, :, i)))
            ! The complete Coulomb triangle below r+s=Rc uses cumulative source moments.
            inner = 0.0_dp
            full_upper = MIN(radius, cutoff_radius - radius)
            IF (full_upper > r(1)) THEN
               source_panel = paw_radial_interval(r, full_upper)
               inner(:) = cumulative(:, source_panel)
               source_h = full_upper - r(source_panel)
               DO k = 1, SIZE(nodes)
                  source_r = r(source_panel) + source_h*nodes(k)
                  t = (source_r - r(source_panel))/(r(source_panel + 1) - r(source_panel))
                  sample(:) = c(0, :, source_panel) + t*(c(1, :, source_panel) + &
                                                         t*(c(2, :, source_panel) + t*c(3, :, source_panel)))
                  inner(:) = inner + source_h*weights(k)*source_r**l*sample
               END DO
               inner(:) = inner/radius**(l + 1)
            END IF
            ! For |r-s|<Rc<r+s, integrate only the resolved partial source panels.
            partial_lower = MAX(r(1), ABS(radius - cutoff_radius))
            IF (partial_lower < radius) THEN
               DO source_panel = paw_radial_interval(r, partial_lower), paw_radial_interval(r, radius)
                  source_left = MAX(r(source_panel), partial_lower)
                  source_right = MIN(r(source_panel + 1), radius)
                  source_h = source_right - source_left
                  IF (source_h <= 0.0_dp) CYCLE
                  DO k = 1, SIZE(nodes)
                     source_r = source_left + source_h*nodes(k)
                     t = (source_r - r(source_panel))/(r(source_panel + 1) - r(source_panel))
                     sample(:) = c(0, :, source_panel) + t*(c(1, :, source_panel) + &
                                                            t*(c(2, :, source_panel) + t*c(3, :, source_panel)))
                     kernel = paw_truncated_multipole(radius, source_r, cutoff_radius, l, angular_nodes, angular_weights)
                     inner(:) = inner + source_h*weights(k)*kernel*sample
                  END DO
               END DO
            END IF
            factor = (right - left)*weights(j)
            DO b = 1, nc
               DO a = 1, b
                  integral(a, b) = integral(a, b) + factor*(value(a)*inner(b) + inner(a)*value(b))
               END DO
            END DO
         END DO
      END DO
      integral = integral*fourpi/REAL(2*l + 1, dp)
      DO b = 1, nc
         DO a = b + 1, nc
            integral(a, b) = integral(b, a)
         END DO
      END DO
   END SUBROUTINE paw_radial_truncated

! **************************************************************************************************
!> \brief Original-mesh bilinear form for erfc(omega*|r-s|)/|r-s|.
!>        Integrate the two radial triangles, including each partial diagonal panel.
!>        Subpanels resolve the screening length; no original orbital tail is shortened.
!> \param r Increasing nonnegative radial coordinates.
!> \param products Radial channels including r**2, without quadrature weights.
!> \param l Angular momentum.
!> \param omega Nonnegative inverse screening length in bohr**-1; zero selects Coulomb.
!> \param integral Symmetric bilinear matrix including 4*pi/(2*l+1).
! **************************************************************************************************
   SUBROUTINE paw_radial_short(r, products, l, omega, integral)
      REAL(dp), INTENT(IN)                               :: r(:), products(:, :)
      INTEGER, INTENT(IN)                                :: l
      REAL(dp), INTENT(IN)                               :: omega
      REAL(dp), INTENT(OUT)                              :: integral(:, :)

      INTEGER                                            :: a, b, chunk, chunks, i, j, k, n, nc, &
                                                            outer, outers, source_panel
      REAL(dp)                                           :: factor, h, kernel, left, lower, &
                                                            nodes(8), radius, right, source_h, &
                                                            source_left, source_r, source_right, &
                                                            t, weights(8)
      REAL(dp), ALLOCATABLE                              :: angular_nodes(:), angular_weights(:), &
                                                            c(:, :, :), inner(:), sample(:), &
                                                            slopes(:, :), value(:)

      n = SIZE(r)
      nc = SIZE(products, 2)
      CPASSERT(n >= 4 .AND. SIZE(products, 1) == n .AND. l >= 0)
      CPASSERT(ALL(SHAPE(integral) == nc))
      CPASSERT(r(1) >= 0.0_dp .AND. ALL(r(2:) > r(:n - 1)))
      CPASSERT(ieee_is_finite(omega) .AND. omega >= 0.0_dp)
      IF (omega == 0.0_dp) THEN
         CALL paw_radial_coulomb(r, products, l, integral)
         RETURN
      END IF
      CPASSERT(omega*(r(n) - r(1)) < REAL(HUGE(0) - 1, dp))
      ALLOCATE (slopes(n, nc), c(0:3, nc, n - 1), inner(nc), sample(nc), value(nc))
      ALLOCATE (angular_nodes(MAX(16, l + 1)), angular_weights(MAX(16, l + 1)))
      CALL paw_angular_gauss(nodes, weights)
      CALL paw_angular_gauss(angular_nodes, angular_weights)
      DO a = 1, nc
         CALL spline3ders(r, products(:, a), r, dynew=slopes(:, a))
      END DO
      DO i = 1, n - 1
         h = r(i + 1) - r(i)
         c(0, :, i) = products(i, :)
         c(1, :, i) = h*slopes(i, :)
         c(2, :, i) = 3.0_dp*(products(i + 1, :) - products(i, :)) - h*(2.0_dp*slopes(i, :) + slopes(i + 1, :))
         c(3, :, i) = 2.0_dp*(products(i, :) - products(i + 1, :)) + h*(slopes(i, :) + slopes(i + 1, :))
      END DO
      integral = 0.0_dp
      DO i = 1, n - 1
         IF (.NOT. ANY(c(:, :, i) /= 0.0_dp)) CYCLE
         h = r(i + 1) - r(i)
         outers = MAX(1, CEILING(omega*h))
         DO outer = 1, outers
            left = r(i) + h*REAL(outer - 1, dp)/REAL(outers, dp)
            right = r(i) + h*REAL(outer, dp)/REAL(outers, dp)
            DO j = 1, SIZE(nodes)
               radius = left + (right - left)*nodes(j)
               t = (radius - r(i))/h
               value(:) = c(0, :, i) + t*(c(1, :, i) + t*(c(2, :, i) + t*c(3, :, i)))
               inner = 0.0_dp
               ! erfc(x) is below the smallest IEEE double at x>=28.
               lower = MAX(r(1), radius - 28.0_dp/omega)
               DO source_panel = paw_radial_interval(r, lower), i
                  IF (.NOT. ANY(c(:, :, source_panel) /= 0.0_dp)) CYCLE
                  source_left = MAX(lower, r(source_panel))
                  source_right = MIN(radius, r(source_panel + 1))
                  source_h = source_right - source_left
                  IF (source_h <= 0.0_dp) CYCLE
                  chunks = MAX(1, CEILING(omega*source_h))
                  DO chunk = 1, chunks
                     DO k = 1, SIZE(nodes)
                        source_r = source_left + source_h*(REAL(chunk - 1, dp) + nodes(k))/REAL(chunks, dp)
                        t = (source_r - r(source_panel))/(r(source_panel + 1) - r(source_panel))
                        sample(:) = c(0, :, source_panel) + t*(c(1, :, source_panel) + &
                                                            t*(c(2, :, source_panel) + t*c(3, :, source_panel)))
                        kernel = paw_short_multipole(radius, source_r, l, omega, angular_nodes, angular_weights)
                        inner(:) = inner + source_h/REAL(chunks, dp)*weights(k)*kernel*sample
                     END DO
                  END DO
               END DO
               factor = (right - left)*weights(j)
               DO b = 1, nc
                  DO a = 1, b
                     integral(a, b) = integral(a, b) + factor*(value(a)*inner(b) + inner(a)*value(b))
                  END DO
               END DO
            END DO
         END DO
      END DO
      integral = integral*fourpi/REAL(2*l + 1, dp)
      DO b = 1, nc
         DO a = b + 1, nc
            integral(a, b) = integral(b, a)
         END DO
      END DO
   END SUBROUTINE paw_radial_short

! **************************************************************************************************
!> \brief Short-range Legendre coefficient integrated in interpoint distance, for r>=s>0.
!>        This removes the coincident-angle Coulomb singularity before quadrature.
!> \param r Outer radius.
!> \param s Inner radius.
!> \param l Angular momentum.
!> \param omega Positive inverse screening length.
!> \param nodes Unit-interval Gauss nodes, with order at least l+1.
!> \param weights Corresponding weights.
!> \return Coefficient v_l without 4*pi/(2*l+1); retain its sign.
! **************************************************************************************************
   FUNCTION paw_short_multipole(r, s, l, omega, nodes, weights) RESULT(value)
      REAL(dp), INTENT(IN)                               :: r, s
      INTEGER, INTENT(IN)                                :: l
      REAL(dp), INTENT(IN)                               :: omega, nodes(:), weights(:)
      REAL(dp)                                           :: value

      INTEGER                                            :: chunk, chunks, i, k
      REAL(dp)                                           :: current, next, previous, span, u, z

      span = MIN(2.0_dp*s, MAX(0.0_dp, 28.0_dp/omega - (r - s)))
      value = 0.0_dp
      IF (span == 0.0_dp) RETURN
      chunks = MAX(1, CEILING(omega*span))
      DO chunk = 1, chunks
         DO i = 1, SIZE(nodes)
            u = span*(REAL(chunk - 1, dp) + nodes(i))/REAL(chunks, dp)
            z = 1.0_dp - u*(2.0_dp*(r - s) + u)/(2.0_dp*r*s)
            z = MIN(1.0_dp, MAX(-1.0_dp, z))
            previous = 1.0_dp
            current = z
            DO k = 2, l
               next = (REAL(2*k - 1, dp)*z*current - REAL(k - 1, dp)*previous)/REAL(k, dp)
               previous = current
               current = next
            END DO
            IF (l == 0) current = 1.0_dp
            value = value + weights(i)*current*ERFC(omega*(r - s + u))
         END DO
      END DO
      value = REAL(2*l + 1, dp)*span/REAL(chunks, dp)*value/(2.0_dp*r*s)
   END FUNCTION paw_short_multipole

! **************************************************************************************************
!> \brief Original mesh interval containing a supplied radius.
!> \param r Increasing radial coordinates.
!> \param x Radius inside the mesh.
!> \return Left endpoint index, including the final closed interval.
! **************************************************************************************************
   FUNCTION paw_radial_interval(r, x) RESULT(interval)
      REAL(dp), INTENT(IN)                               :: r(:), x
      INTEGER                                            :: interval

      INTEGER                                            :: high, low, middle

      low = 1
      high = SIZE(r)
      DO WHILE (high - low > 1)
         middle = (low + high)/2
         IF (r(middle) > x) THEN
            high = middle
         ELSE
            low = middle
         END IF
      END DO
      interval = low
   END FUNCTION paw_radial_interval

! **************************************************************************************************
!> \brief Partial truncated Legendre coefficient for the radial triangle r>=s.
!> \param r Outer radius.
!> \param s Inner radius.
!> \param cutoff_radius Interaction radius.
!> \param l Angular momentum.
!> \param nodes Unit-interval Gauss nodes of order l+1.
!> \param weights Corresponding weights.
!> \return Coefficient v_l, without 4*pi/(2*l+1); its sign must not be clipped.
! **************************************************************************************************
   FUNCTION paw_truncated_multipole(r, s, cutoff_radius, l, nodes, weights) RESULT(value)
      REAL(dp), INTENT(IN)                               :: r, s, cutoff_radius
      INTEGER, INTENT(IN)                                :: l
      REAL(dp), INTENT(IN)                               :: nodes(:), weights(:)
      REAL(dp)                                           :: value

      INTEGER                                            :: i, k
      REAL(dp)                                           :: current, next, previous, span, u, z

      span = MIN(2.0_dp*s, MAX(0.0_dp, cutoff_radius - r + s))
      value = 0.0_dp
      DO i = 1, SIZE(nodes)
         u = span*nodes(i)
         ! Factored difference of squares retains the angular coordinate for s<<r.
         z = 1.0_dp - u*(2.0_dp*(r - s) + u)/(2.0_dp*r*s)
         z = MIN(1.0_dp, MAX(-1.0_dp, z))
         previous = 1.0_dp
         current = z
         DO k = 2, l
            next = (REAL(2*k - 1, dp)*z*current - REAL(k - 1, dp)*previous)/REAL(k, dp)
            previous = current
            current = next
         END DO
         IF (l == 0) current = 1.0_dp
         value = value + weights(i)*current
      END DO
      value = REAL(2*l + 1, dp)*span*value/(2.0_dp*r*s)
   END FUNCTION paw_truncated_multipole

! **************************************************************************************************
!> \brief Generate a unit-interval Gauss rule with the order required by the angular polynomial.
!> \param nodes Gauss nodes; their array size determines the order.
!> \param weights Corresponding weights.
! **************************************************************************************************
   SUBROUTINE paw_angular_gauss(nodes, weights)
      REAL(dp), INTENT(OUT)                              :: nodes(:), weights(:)

      INTEGER                                            :: i, iteration, k, n
      REAL(dp)                                           :: current, derivative, next, previous, &
                                                            step, z

      n = SIZE(nodes)
      CPASSERT(n > 0 .AND. SIZE(weights) == n)
      DO i = 1, (n + 1)/2
         z = COS(pi*(REAL(i, dp) - 0.25_dp)/(REAL(n, dp) + 0.5_dp))
         DO iteration = 1, 100
            previous = 1.0_dp
            current = z
            DO k = 2, n
               next = (REAL(2*k - 1, dp)*z*current - REAL(k - 1, dp)*previous)/REAL(k, dp)
               previous = current
               current = next
            END DO
            derivative = REAL(n, dp)*(z*current - previous)/(z*z - 1.0_dp)
            step = current/derivative
            z = z - step
            IF (ABS(step) <= 4.0_dp*EPSILON(1.0_dp)) EXIT
         END DO
         CPASSERT(iteration <= 100)
         nodes(i) = 0.5_dp*(1.0_dp - z)
         nodes(n + 1 - i) = 0.5_dp*(1.0_dp + z)
         weights(i) = 1.0_dp/((1.0_dp - z*z)*derivative**2)
         weights(n + 1 - i) = weights(i)
      END DO
   END SUBROUTINE paw_angular_gauss

END MODULE paw_onecenter
