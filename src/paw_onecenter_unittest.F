!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM paw_onecenter_unittest
   USE atom_upf,                        ONLY: atom_read_upf,&
                                              atom_release_upf,&
                                              atom_upfpot_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE libcp2k,                         ONLY: cp2k_finalize,&
                                              cp2k_init
   USE mathconstants,                   ONLY: fourpi,&
                                              pi
   USE paw_onecenter,                   ONLY: paw_onecenter_exchange,&
                                              paw_onecenter_hartree,&
                                              paw_onecenter_init,&
                                              paw_onecenter_type,&
                                              paw_radial_coulomb
   USE upf_projector_integrals,         ONLY: upf_projector_init,&
                                              upf_projector_type

   IMPLICIT NONE

   CALL cp2k_init()
   CALL radial_integrals()
   CALL projector_functional()
   CALL original_datasets()
   CALL cp2k_finalize()

CONTAINS

! **************************************************************************************************
!> \brief Optional original PAW files: finite reference-occupation Hartree energy and matrix.
! **************************************************************************************************
   SUBROUTINE original_datasets()
      CHARACTER(LEN=default_path_length)                 :: filename
      INTEGER                                            :: argument, i, np, p
      LOGICAL                                            :: source_root_argument
      REAL(dp)                                           :: energy
      REAL(dp), ALLOCATABLE                              :: density(:, :), potential(:, :)
      TYPE(atom_upfpot_type)                             :: pot
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(paw_onecenter_type)                           :: onecenter
      TYPE(upf_projector_type)                           :: projectors

      logger => cp_get_default_logger()
      DO argument = 1, COMMAND_ARGUMENT_COUNT()
         CALL GET_COMMAND_ARGUMENT(argument, filename)
         IF (argument == 1) THEN
            INQUIRE (FILE=TRIM(filename)//'/tests/UNIT_TESTS', EXIST=source_root_argument)
            IF (source_root_argument) CYCLE
         END IF
         CALL atom_read_upf(pot, TRIM(filename), read_paw=.TRUE.)
         CALL upf_projector_init(pot, MAXVAL(pot%lbeta), projectors)
         CALL paw_onecenter_init(pot, projectors, onecenter)
         np = projectors%nproj
         ALLOCATE (density(np, np), potential(np, np))
         density = 0.0_dp
         DO p = 1, np
            i = projectors%radial_index(p)
            density(p, p) = pot%paw%occupations(i)/REAL(2*pot%lbeta(i) + 1, dp)
         END DO
         CALL paw_onecenter_hartree(onecenter, density, energy, potential)
         IF (.NOT. ieee_is_finite(energy) .OR. .NOT. ALL(ieee_is_finite(potential))) THEN
            ERROR STOP 'Nonfinite PAW reference-occupation Hartree correction'
         END IF
         IF (logger%para_env%is_source()) THEN
            WRITE (*, '(A,1X,A,2(1X,ES25.16))') 'PAW Hartree:', TRIM(filename), energy, MAXVAL(ABS(potential))
         END IF
         DEALLOCATE (density, potential)
         CALL atom_release_upf(pot)
      END DO
   END SUBROUTINE original_datasets

! **************************************************************************************************
!> \brief Analytic Coulomb integral for r**(l+2)*exp(-a*r**2) and exponent b.
!> \param a First exponent.
!> \param b Second exponent.
!> \param l Angular momentum.
!> \return Infinite-domain bilinear form including 4*pi/(2*l+1).
! **************************************************************************************************
   FUNCTION gaussian_coulomb(a, b, l) RESULT(value)
      REAL(dp), INTENT(IN)                               :: a, b
      INTEGER, INTENT(IN)                                :: l
      REAL(dp)                                           :: value

      value = pi*GAMMA(l + 1.5_dp)/(REAL(2*l + 1, dp)*a*b*(a + b)**(l + 0.5_dp))
   END FUNCTION gaussian_coulomb

! **************************************************************************************************
!> \brief Verify independent analytic values and radial convergence on linear and log meshes.
! **************************************************************************************************
   SUBROUTINE radial_integrals()
      INTEGER                                            :: i, j, l, mesh, n, refinement
      REAL(dp)                                           :: error, exact, exponents(2), &
                                                            integral(2, 2), previous
      REAL(dp), ALLOCATABLE                              :: products(:, :), r(:)

      exponents = [0.7_dp, 1.3_dp]
      DO mesh = 1, 2
         previous = HUGE(0.0_dp)
         DO refinement = 1, 3
            n = MERGE(64, 200, mesh == 1)*2**refinement + 1
            ALLOCATE (r(n), products(n, 2))
            DO i = 1, n
               IF (mesh == 1) THEN
                  r(i) = 10.0_dp*REAL(i - 1, dp)/REAL(n - 1, dp)
               ELSE
                  r(i) = 1.e-7_dp*EXP(LOG(1.e8_dp)*REAL(i - 1, dp)/REAL(n - 1, dp))
               END IF
            END DO
            error = 0.0_dp
            DO l = 0, 6
               DO i = 1, 2
                  products(:, i) = r**(l + 2)*EXP(-exponents(i)*r**2)
               END DO
               CALL paw_radial_coulomb(r, products, l, integral)
               DO j = 1, 2
                  DO i = 1, 2
                     exact = gaussian_coulomb(exponents(i), exponents(j), l)
                     error = MAX(error, ABS(integral(i, j)/exact - 1.0_dp))
                  END DO
               END DO
               IF (integral(1, 2) /= integral(2, 1)) ERROR STOP 'Asymmetric radial Coulomb matrix'
            END DO
            IF (refinement > 1 .AND. previous > 1.e-11_dp) THEN
               IF (error > 0.2_dp*previous) ERROR STOP 'PAW radial Coulomb integrals do not converge'
            END IF
            WRITE (*, '(A,2I6,ES15.6)') 'PAW radial Coulomb mesh/refinement/error:', mesh, refinement, error
            previous = error
            DEALLOCATE (r, products)
         END DO
         IF (error > 3.e-8_dp) ERROR STOP 'PAW radial Coulomb integral differs from analytic Gaussian result'
      END DO
   END SUBROUTINE radial_integrals

! **************************************************************************************************
!> \brief Analytic compensated s channel, nonspherical matrix derivatives and radial cutoff.
!> \brief Analytic compensated s channel, nonspherical matrix derivatives and radial cutoff.
!> \brief Analytic compensated s channel, nonspherical matrix derivatives and radial cutoff.
!> \brief Analytic compensated s channel, nonspherical matrix derivatives and radial cutoff.
!> \return ...
! **************************************************************************************************
   SUBROUTINE projector_functional()
      INTEGER, PARAMETER                                 :: cutoff = 449, n = 513, nb = 3

      INTEGER                                            :: i, j, l, np, p, pair, q
      REAL(dp)                                           :: ae_exponents(nb), directional, energy, &
                                                            exact, minus, plus, ps_exponents(nb), &
                                                            step
      REAL(dp), ALLOCATABLE                              :: density(:, :), direction(:, :), &
                                                            potential(:, :), scratch(:, :)
      TYPE(atom_upfpot_type)                             :: pot
      TYPE(paw_onecenter_type)                           :: changed_tail, onecenter
      TYPE(upf_projector_type)                           :: projectors

      pot%is_paw = .TRUE.
      pot%number_of_proj = nb
      ALLOCATE (pot%paw, pot%r(n), pot%rab(n), pot%lbeta(nb), pot%beta(n, nb), &
                pot%dion(nb, nb), pot%qijl(n, nb*(nb + 1)/2, 0:4))
      ALLOCATE (pot%paw%ae_wfc(n, nb), pot%paw%ps_wfc(n, nb))
      pot%paw%cutoff_index = cutoff
      pot%lbeta(:) = [0, 1, 2]
      pot%dion = 0.0_dp
      pot%qijl = 0.0_dp
      pot%rab = 8.0_dp/REAL(n - 1, dp)
      DO i = 1, n
         pot%r(i) = REAL(i - 1, dp)*pot%rab(i)
      END DO
      ae_exponents = [0.9_dp, 1.2_dp, 1.5_dp]
      ps_exponents = [0.7_dp, 1.0_dp, 1.3_dp]
      DO i = 1, nb
         pot%beta(:, i) = pot%r**pot%lbeta(i)*EXP(-pot%r**2)
         pot%paw%ae_wfc(:, i) = pot%r**(pot%lbeta(i) + 1)*EXP(-0.5_dp*ae_exponents(i)*pot%r**2)
         pot%paw%ps_wfc(:, i) = 0.8_dp*pot%r**(pot%lbeta(i) + 1)*EXP(-0.5_dp*ps_exponents(i)*pot%r**2)
      END DO
      DO j = 1, nb
         DO i = 1, j
            pair = j*(j - 1)/2 + i
            DO l = ABS(pot%lbeta(i) - pot%lbeta(j)), pot%lbeta(i) + pot%lbeta(j), 2
               pot%qijl(:, pair, l) = 0.07_dp*pot%r**(l + 2)*EXP(-2.0_dp*pot%r**2)
            END DO
         END DO
      END DO
      CALL upf_projector_init(pot, 2, projectors)
      CALL paw_onecenter_init(pot, projectors, onecenter)
      np = projectors%nproj
      ALLOCATE (density(np, np), direction(np, np), potential(np, np), scratch(np, np))
      density = 0.0_dp
      density(1, 1) = 1.0_dp
      CALL paw_onecenter_hartree(onecenter, density, energy, potential)
      exact = 0.5_dp/fourpi*(gaussian_coulomb(0.9_dp, 0.9_dp, 0) &
                             - 0.8_dp**4*gaussian_coulomb(0.7_dp, 0.7_dp, 0) &
                             - 2.0_dp*0.8_dp**2*0.07_dp*gaussian_coulomb(0.7_dp, 2.0_dp, 0) &
                             - 0.07_dp**2*gaussian_coulomb(2.0_dp, 2.0_dp, 0))
      IF (ABS(energy - exact) > 2.e-9_dp) ERROR STOP 'Wrong PAW compensation, Hartree prefactor or angular normalization'
      DO q = 1, np
         DO p = 1, np
            density(p, q) = SIN(0.31_dp*p)*SIN(0.31_dp*q) + 0.4_dp*COS(0.47_dp*p)*COS(0.47_dp*q)
            direction(p, q) = COS(0.23_dp*(p + q))
         END DO
      END DO
      CALL paw_onecenter_hartree(onecenter, density, energy, potential)
      IF (MAXVAL(ABS(potential - TRANSPOSE(potential))) > 1.e-13_dp) ERROR STOP 'Asymmetric PAW Hartree derivative'
      step = 1.e-5_dp
      CALL paw_onecenter_hartree(onecenter, density + step*direction, plus, scratch)
      CALL paw_onecenter_hartree(onecenter, density - step*direction, minus, scratch)
      directional = SUM(potential*direction)
      IF (ABS((plus - minus)/(2.0_dp*step) - directional) > 1.e-9_dp*MAX(1.0_dp, ABS(directional))) THEN
         ERROR STOP 'PAW projector-density derivative differs from energy finite difference'
      END IF
      CALL exchange_functional(onecenter, projectors, exact)
      pot%paw%ae_wfc(cutoff + 1:, :) = 1.e8_dp
      pot%paw%ps_wfc(cutoff + 1:, :) = -2.e8_dp
      pot%qijl(cutoff + 1:, :, :) = 1.e20_dp
      CALL paw_onecenter_init(pot, projectors, changed_tail)
      IF (ANY(changed_tail%hartree /= onecenter%hartree)) ERROR STOP 'PAW radial cutoff ignored'
   END SUBROUTINE projector_functional

! **************************************************************************************************
!> \brief Independent crossed-index exchange tensor, spin and full-matrix derivative checks.
!> \param onecenter Original s/p/d Coulomb cache.
!> \param projectors Independent radial-to-magnetic mapping for the four-index oracle.
!> \param s_hartree Analytic single-occupation compensated s-channel Hartree energy.
! **************************************************************************************************
   SUBROUTINE exchange_functional(onecenter, projectors, s_hartree)
      TYPE(paw_onecenter_type), INTENT(IN)               :: onecenter
      TYPE(upf_projector_type), INTENT(IN)               :: projectors
      REAL(dp), INTENT(IN)                               :: s_hartree

      REAL(dp), PARAMETER :: fractions(3) = [0.0_dp, 0.37_dp, 1.0_dp], step = 1.e-5_dp

      INTEGER                                            :: i, ij, j, k, kl, l, lm, np, p, q, spin, variant
      INTEGER, ALLOCATABLE                               :: permutation(:)
      REAL(dp)                                           :: directional, energy, fraction, minus, plus, &
                                                            reference_energy, scale, single_energy, &
                                                            tensor_error, uncrossed_energy
      REAL(dp), ALLOCATABLE                              :: coulomb(:, :, :, :), density(:, :, :), &
                                                            direction(:, :, :), permuted_density(:, :, :), &
                                                            potential(:, :, :), reference(:, :, :), scratch(:, :, :)
      TYPE(paw_onecenter_type)                           :: permuted

      np = onecenter%nproj
      ALLOCATE (coulomb(np, np, np, np), density(np, np, 2), direction(np, np, 2), &
                permuted_density(np, np, 2), potential(np, np, 2), reference(np, np, 2), &
                scratch(np, np, 2), permutation(np))
      coulomb(:, :, :, :) = 0.0_dp
      ! Build (ij|kl) before contracting densities, with independent radial-pair
      ! indices from the projector mapping rather than the production pair cache.
      DO l = 1, np
         DO k = 1, np
            kl = MAX(projectors%radial_index(k), projectors%radial_index(l))
            kl = kl*(kl - 1)/2 + MIN(projectors%radial_index(k), projectors%radial_index(l))
            DO j = 1, np
               DO i = 1, np
                  ij = MAX(projectors%radial_index(i), projectors%radial_index(j))
                  ij = ij*(ij - 1)/2 + MIN(projectors%radial_index(i), projectors%radial_index(j))
                  DO lm = 1, SIZE(onecenter%gaunt, 3)
                     p = INT(SQRT(REAL(lm - 1, dp)))
                     coulomb(i, j, k, l) = coulomb(i, j, k, l) + &
                                          onecenter%gaunt(i, j, lm)*onecenter%gaunt(k, l, lm)* &
                                          onecenter%hartree(ij, kl, p)
                  END DO
               END DO
            END DO
         END DO
      END DO
      tensor_error = 0.0_dp
      DO l = 1, np
         DO k = 1, np
            DO j = 1, np
               DO i = 1, np
                  tensor_error = MAX(tensor_error, ABS(coulomb(i, j, k, l) - coulomb(j, i, k, l)), &
                                     ABS(coulomb(i, j, k, l) - coulomb(i, j, l, k)), &
                                     ABS(coulomb(i, j, k, l) - coulomb(k, l, i, j)))
               END DO
            END DO
         END DO
      END DO
      IF (tensor_error > 1.e-13_dp*MAX(1.0_dp, MAXVAL(ABS(coulomb)))) &
         ERROR STOP 'PAW exchange Coulomb tensor permutations differ'

      density(:, :, :) = 0.0_dp
      density(1, 1, 1) = 1.0_dp
      CALL paw_onecenter_exchange(onecenter, density, fractions(2), energy, potential)
      IF (ABS(energy + fractions(2)*s_hartree) > 2.e-9_dp) &
         ERROR STOP 'Wrong PAW exchange fraction, sign, compensation or angular normalization'
      IF (ABS(potential(1, 1, 1) + 2.0_dp*fractions(2)*s_hartree) > 4.e-9_dp) &
         ERROR STOP 'Wrong PAW exchange single-spin matrix factor'

      DO spin = 1, 2
         DO q = 1, np
            DO p = 1, np
               density(p, q, spin) = SIN((0.17_dp + 0.13_dp*spin)*p)*SIN((0.17_dp + 0.13_dp*spin)*q) + &
                                     0.4_dp*COS((0.41_dp + 0.07_dp*spin)*p)*COS((0.41_dp + 0.07_dp*spin)*q)
            END DO
         END DO
      END DO
      DO variant = 1, SIZE(fractions)
         fraction = fractions(variant)
         reference(:, :, :) = 0.0_dp
         reference_energy = 0.0_dp
         uncrossed_energy = 0.0_dp
         DO spin = 1, 2
            DO l = 1, np
               DO k = 1, np
                  DO j = 1, np
                     DO i = 1, np
                        reference_energy = reference_energy - 0.5_dp*fraction* &
                                           density(i, k, spin)*density(j, l, spin)*coulomb(i, j, k, l)
                        uncrossed_energy = uncrossed_energy - 0.5_dp*fraction* &
                                          density(i, j, spin)*density(k, l, spin)*coulomb(i, j, k, l)
                        reference(i, k, spin) = reference(i, k, spin) - &
                                               fraction*density(j, l, spin)*coulomb(i, j, k, l)
                     END DO
                  END DO
               END DO
            END DO
         END DO
         CALL paw_onecenter_exchange(onecenter, density, fraction, energy, potential)
         scale = MAX(1.0_dp, ABS(reference_energy), MAXVAL(ABS(reference)))
         IF (ABS(energy - reference_energy) > 1.e-12_dp*scale .OR. &
             MAXVAL(ABS(potential - reference)) > 1.e-12_dp*scale) &
            ERROR STOP 'PAW exchange differs from independent crossed-index four-index contraction'
         IF (fraction > 0.0_dp) THEN
            IF (ABS(uncrossed_energy - reference_energy) < 1.e-6_dp*scale) &
               ERROR STOP 'PAW exchange fixture cannot distinguish Hartree and exchange density indices'
         END IF
         DO spin = 1, 2
            IF (ANY(potential(:, :, spin) /= TRANSPOSE(potential(:, :, spin)))) &
               ERROR STOP 'Asymmetric PAW exchange derivative'
            DO q = 1, np
               DO p = 1, np
                  direction(p, q, :) = 0.0_dp
                  direction(p, q, spin) = COS(0.23_dp*(p + q) + 0.19_dp*spin)
               END DO
            END DO
            CALL paw_onecenter_exchange(onecenter, density + step*direction, fraction, plus, scratch)
            CALL paw_onecenter_exchange(onecenter, density - step*direction, fraction, minus, scratch)
            directional = SUM(potential*direction)
            IF (ABS((plus - minus)/(2.0_dp*step) - directional) > 1.e-9_dp*scale) &
               ERROR STOP 'PAW exchange spin derivative differs from energy finite difference'
            direction(:, :, :) = 0.0_dp
            direction(2, 6, spin) = 0.7_dp
            direction(6, 2, spin) = 0.7_dp
            CALL paw_onecenter_exchange(onecenter, density + step*direction, fraction, plus, scratch)
            IF (ANY(scratch(:, :, 3 - spin) /= potential(:, :, 3 - spin))) &
               ERROR STOP 'PAW exchange coupled different spin density blocks'
            CALL paw_onecenter_exchange(onecenter, density - step*direction, fraction, minus, scratch)
            directional = 1.4_dp*potential(2, 6, spin)
            IF (ABS((plus - minus)/(2.0_dp*step) - directional) > 1.e-9_dp*scale) &
               ERROR STOP 'PAW exchange off-diagonal full-matrix factor differs from energy finite difference'
         END DO
      END DO

      ! Relabel radial/magnetic projectors, including mixed-l density entries.
      permuted = onecenter
      permutation(:) = [(np + 1 - p, p=1, np)]
      DO q = 1, np
         DO p = 1, np
            permuted%pair_index(p, q) = onecenter%pair_index(permutation(p), permutation(q))
            permuted%gaunt(p, q, :) = onecenter%gaunt(permutation(p), permutation(q), :)
            permuted_density(p, q, :) = density(permutation(p), permutation(q), :)
         END DO
      END DO
      CALL paw_onecenter_exchange(permuted, permuted_density, fraction, plus, scratch)
      IF (ABS(plus - energy) > 1.e-12_dp*scale) ERROR STOP 'PAW exchange depends on projector ordering'
      DO q = 1, np
         DO p = 1, np
            IF (MAXVAL(ABS(scratch(p, q, :) - potential(permutation(p), permutation(q), :))) > 1.e-12_dp*scale) &
               ERROR STOP 'PAW exchange matrix does not follow projector permutation'
         END DO
      END DO

      ! A single physical spin channel and two half-density restricted channels
      ! differ by exactly a factor two in both energy and per-spin potential.
      CALL paw_onecenter_exchange(onecenter, density(:, :, 1:1), fraction, single_energy, scratch(:, :, 1:1))
      density(:, :, 2) = 0.5_dp*density(:, :, 1)
      density(:, :, 1) = density(:, :, 2)
      CALL paw_onecenter_exchange(onecenter, density, fraction, energy, potential)
      IF (ABS(energy - 0.5_dp*single_energy) > 1.e-12_dp*scale .OR. &
          MAXVAL(ABS(potential(:, :, 1) - 0.5_dp*scratch(:, :, 1))) > 1.e-12_dp*scale .OR. &
          ANY(potential(:, :, 1) /= potential(:, :, 2))) &
         ERROR STOP 'PAW exchange restricted spin-summed conversion factor differs'
   END SUBROUTINE exchange_functional

END PROGRAM paw_onecenter_unittest
