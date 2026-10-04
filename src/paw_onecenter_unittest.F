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
   USE paw_onecenter,                   ONLY: &
        paw_onecenter_core_exchange, paw_onecenter_core_exchange_init, &
        paw_onecenter_core_exchange_prepare, paw_onecenter_exchange, paw_onecenter_hartree, &
        paw_onecenter_init, paw_onecenter_type, paw_radial_coulomb, paw_radial_short, &
        paw_radial_truncated
   USE upf_projector_integrals,         ONLY: upf_projector_init,&
                                              upf_projector_type

   IMPLICIT NONE

   CALL cp2k_init()
   CALL radial_integrals()
   CALL truncated_radial_integrals()
   CALL screened_radial_integrals()
   CALL projector_energy()
   CALL core_exchange_checks()
   CALL interaction_core_checks()
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
!> \brief Independent Fourier/Hermite reference for the truncated Gaussian multipole bilinear.
!> \param a First exponent.
!> \param b Second exponent.
!> \param l Angular momentum.
!> \param cutoff Interaction radius in bohr.
!> \return Infinite-domain truncated bilinear, with the same normalization as gaussian_coulomb.
! **************************************************************************************************
   FUNCTION gaussian_truncated(a, b, l, cutoff) RESULT(value)
      REAL(dp), INTENT(IN)                               :: a, b
      INTEGER, INTENT(IN)                                :: l
      REAL(dp), INTENT(IN)                               :: cutoff
      REAL(dp)                                           :: value

      INTEGER                                            :: k
      REAL(dp)                                           :: current, next, previous, x

      x = a*b*cutoff**2/(a + b)
      previous = 1.0_dp
      current = 2.0_dp*SQRT(x)
      DO k = 2, 2*l
         next = 2.0_dp*SQRT(x)*current - 2.0_dp*REAL(k - 1, dp)*previous
         previous = current
         current = next
      END DO
      IF (l == 0) current = 1.0_dp
      value = gaussian_coulomb(a, b, l)*(1.0_dp - REAL((-1)**l, dp)* &
                                         GAMMA(REAL(l + 1, dp))/GAMMA(REAL(2*l + 1, dp))*current*EXP(-x))
   END FUNCTION gaussian_truncated

! **************************************************************************************************
!> \brief Fourier Gaussian multipoles and a compact exact sphere on a logarithmic mesh.
!>        The l=1 reference exceeds its Coulomb value for a selected radius, so clipping is detectable.
! **************************************************************************************************
   SUBROUTINE truncated_radial_integrals()
      INTEGER, PARAMETER                                 :: angular(4) = [0, 1, 2, 8]
      REAL(dp), PARAMETER :: amplitudes(2) = [1.0_dp, 0.73_dp], cutoffs(2) = [0.7_dp, 2.3_dp], &
         exponents(2) = [0.7_dp, 1.1_dp]

      INTEGER                                            :: a, b, channel, i, l, mesh, n, variant
      REAL(dp)                                           :: error, exact, integral(2, 2), &
                                                            reference(2, 2), scale
      REAL(dp), ALLOCATABLE                              :: products(:, :), r(:)

      DO mesh = 1, 2
         n = MERGE(513, 17, mesh == 1)
         ALLOCATE (r(n), products(n, 2))
         DO i = 1, n
            IF (mesh == 1) THEN
               r(i) = 10.0_dp*REAL(i - 1, dp)/REAL(n - 1, dp)
            ELSE
               r(i) = 1.e-7_dp*EXP(LOG(1.e7_dp)*REAL(i - 1, dp)/REAL(n - 1, dp))
            END IF
         END DO
         error = 0.0_dp
         DO channel = 1, MERGE(SIZE(angular), 1, mesh == 1)
            l = angular(channel)
            DO a = 1, 2
               IF (mesh == 1) THEN
                  products(:, a) = r**(l + 2)*EXP(-exponents(a)*r**2)
               ELSE
                  products(:, a) = amplitudes(a)*3.0_dp*r**2/SQRT(fourpi)
               END IF
            END DO
            DO variant = 1, SIZE(cutoffs)
               CALL paw_radial_truncated(r, products, l, cutoffs(variant), integral)
               DO b = 1, 2
                  DO a = 1, 2
                     IF (mesh == 1) THEN
                        exact = gaussian_truncated(exponents(a), exponents(b), l, cutoffs(variant))
                        scale = MAX(1.0_dp, ABS(gaussian_coulomb(exponents(a), exponents(b), l)))
                     ELSE
                        exact = amplitudes(a)*amplitudes(b)*sphere_kernel(cutoffs(variant), -1.0_dp)
                        scale = 1.0_dp
                     END IF
                     error = MAX(error, ABS(integral(a, b) - exact)/scale)
                  END DO
               END DO
               IF (integral(1, 2) /= integral(2, 1)) ERROR STOP 'Asymmetric truncated PAW radial matrix'
               IF (l == 1 .AND. variant == 2) THEN
                  IF (integral(1, 2) <= gaussian_coulomb(exponents(1), exponents(2), l)) THEN
                     ERROR STOP 'Truncated multipole incorrectly bounded by Coulomb'
                  END IF
               END IF
            END DO
            CALL paw_radial_truncated(r, products, l, 0.0_dp, integral)
            IF (ANY(integral /= 0.0_dp)) ERROR STOP 'Zero-radius interaction is nonzero'
            CALL paw_radial_truncated(r, products, l, 2.0_dp*r(n), integral)
            CALL paw_radial_coulomb(r, products, l, reference)
            IF (ANY(integral /= reference)) ERROR STOP 'Saturated truncated radial kernel differs from Coulomb'
            DO b = 1, 2
               DO a = 1, 2
                  IF (mesh == 1) THEN
                     exact = gaussian_coulomb(exponents(a), exponents(b), l)
                  ELSE
                     exact = amplitudes(a)*amplitudes(b)*sphere_kernel(2.0_dp, -1.0_dp)
                  END IF
                  error = MAX(error, ABS(reference(a, b) - exact)/MAX(1.0_dp, ABS(exact)))
               END DO
            END DO
         END DO
         WRITE (*, '(A,I3,ES15.6)') 'PAW truncated radial mesh/error:', mesh, error
         IF (error > 3.e-8_dp) ERROR STOP 'Truncated radial PAW kernels differ from analytic references'
         DEALLOCATE (r, products)
      END DO
   END SUBROUTINE truncated_radial_integrals

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
! **************************************************************************************************
   SUBROUTINE projector_energy()
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
   END SUBROUTINE projector_energy

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

      INTEGER                                            :: i, ij, j, k, kl, l, lm, np, p, q, spin, &
                                                            variant
      INTEGER, ALLOCATABLE                               :: permutation(:)
      REAL(dp)                                           :: directional, energy, fraction, minus, &
                                                            plus, reference_energy, scale, &
                                                            single_energy, tensor_error, &
                                                            uncrossed_energy
      REAL(dp), ALLOCATABLE :: coulomb(:, :, :, :), density(:, :, :), direction(:, :, :), &
         permuted_density(:, :, :), potential(:, :, :), reference(:, :, :), scratch(:, :, :)
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
      IF (tensor_error > 1.e-13_dp*MAX(1.0_dp, MAXVAL(ABS(coulomb)))) THEN
         ERROR STOP 'PAW exchange Coulomb tensor permutations differ'
      END IF

      density(:, :, :) = 0.0_dp
      density(1, 1, 1) = 1.0_dp
      CALL paw_onecenter_exchange(onecenter, density, fractions(2), energy, potential)
      IF (ABS(energy + fractions(2)*s_hartree) > 2.e-9_dp) THEN
         ERROR STOP 'Wrong PAW exchange fraction, sign, compensation or angular normalization'
      END IF
      IF (ABS(potential(1, 1, 1) + 2.0_dp*fractions(2)*s_hartree) > 4.e-9_dp) THEN
         ERROR STOP 'Wrong PAW exchange single-spin matrix factor'
      END IF

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
             MAXVAL(ABS(potential - reference)) > 1.e-12_dp*scale) THEN
            ERROR STOP 'PAW exchange differs from independent crossed-index four-index contraction'
         END IF
         IF (fraction > 0.0_dp) THEN
            IF (ABS(uncrossed_energy - reference_energy) < 1.e-6_dp*scale) THEN
               ERROR STOP 'PAW exchange fixture cannot distinguish Hartree and exchange density indices'
            END IF
         END IF
         DO spin = 1, 2
            IF (ANY(potential(:, :, spin) /= TRANSPOSE(potential(:, :, spin)))) THEN
               ERROR STOP 'Asymmetric PAW exchange derivative'
            END IF
            DO q = 1, np
               DO p = 1, np
                  direction(p, q, :) = 0.0_dp
                  direction(p, q, spin) = COS(0.23_dp*(p + q) + 0.19_dp*spin)
               END DO
            END DO
            CALL paw_onecenter_exchange(onecenter, density + step*direction, fraction, plus, scratch)
            CALL paw_onecenter_exchange(onecenter, density - step*direction, fraction, minus, scratch)
            directional = SUM(potential*direction)
            IF (ABS((plus - minus)/(2.0_dp*step) - directional) > 1.e-9_dp*scale) THEN
               ERROR STOP 'PAW exchange spin derivative differs from energy finite difference'
            END IF
            direction(:, :, :) = 0.0_dp
            direction(2, 6, spin) = 0.7_dp
            direction(6, 2, spin) = 0.7_dp
            CALL paw_onecenter_exchange(onecenter, density + step*direction, fraction, plus, scratch)
            IF (ANY(scratch(:, :, 3 - spin) /= potential(:, :, 3 - spin))) THEN
               ERROR STOP 'PAW exchange coupled different spin density blocks'
            END IF
            CALL paw_onecenter_exchange(onecenter, density - step*direction, fraction, minus, scratch)
            directional = 1.4_dp*potential(2, 6, spin)
            IF (ABS((plus - minus)/(2.0_dp*step) - directional) > 1.e-9_dp*scale) THEN
               ERROR STOP 'PAW exchange off-diagonal full-matrix factor differs from energy finite difference'
            END IF
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
            IF (MAXVAL(ABS(scratch(p, q, :) - potential(permutation(p), permutation(q), :))) > 1.e-12_dp*scale) THEN
               ERROR STOP 'PAW exchange matrix does not follow projector permutation'
            END IF
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
          ANY(potential(:, :, 1) /= potential(:, :, 2))) THEN
         ERROR STOP 'PAW exchange restricted spin-summed conversion factor differs'
      END IF
   END SUBROUTINE exchange_functional

! **************************************************************************************************
!> \brief Independent Cartesian-Gaussian closed s/p core: full tails, amplitudes, spin and derivatives.
! **************************************************************************************************
   SUBROUTINE core_exchange_checks()
      INTEGER, PARAMETER                                 :: n = 1025, nb = 3, np = 5
      REAL(dp), PARAMETER :: amplitudes(nb) = [1.0_dp, 0.73_dp, 1.0_dp], &
         core_exponents(2) = [2.3_dp, 0.31_dp], fractions(3) = [0.0_dp, 0.37_dp, 1.0_dp], &
         step = 1.e-5_dp, valence_exponents(nb) = [0.43_dp, 0.83_dp, 0.67_dp]

      INTEGER                                            :: cache_index, i, j, spin, variant
      REAL(dp) :: a, b, density(np, np, 2), direction(np, np, 2), directional, energy, error, &
         exact, fraction, minus, plus, potential(np, np, 2), radial(n, 2), reference(np, np), &
         restricted(np, np, 2), scale, scratch(np, np, 2)
      TYPE(atom_upfpot_type)                             :: pot
      TYPE(paw_onecenter_type)                           :: onecenter, phase_changed, unavailable
      TYPE(upf_projector_type)                           :: projectors

      pot%is_paw = .TRUE.
      pot%symbol = 'Ne'
      pot%zion = 2.0_dp
      pot%number_of_proj = nb
      ALLOCATE (pot%paw, pot%r(n), pot%rab(n), pot%lbeta(nb), pot%beta(n, nb), &
                pot%dion(nb, nb), pot%qijl(n, nb*(nb + 1)/2, 0:2))
      ALLOCATE (pot%paw%ae_wfc(n, nb), pot%paw%ps_wfc(n, nb), pot%paw%rho_core(n))
      ALLOCATE (pot%paw%gipaw_core_orbitals)
      ASSOCIATE (core => pot%paw%gipaw_core_orbitals)
         ALLOCATE (core%wfc(n, 2), core%indices(2), core%principal_quantum_numbers(2), core%angular_momenta(2))
         core%indices(:) = [1, 2]
         core%principal_quantum_numbers(:) = [1, 2]
         core%angular_momenta(:) = [0, 1]
      END ASSOCIATE
      pot%paw%cutoff_index = 205
      pot%lbeta(:) = [0, 0, 1]
      pot%rab = 10.0_dp/REAL(n - 1, dp)
      pot%r(:) = [(REAL(i - 1, dp)*pot%rab(i), i=1, n)]
      pot%dion = 0.0_dp
      pot%qijl = 0.0_dp
      DO i = 1, 2
         a = core_exponents(i)
         radial(:, i) = SQRT(fourpi)*(2.0_dp*a/pi)**0.75_dp*EXP(-a*pot%r**2)
         IF (i == 2) radial(:, i) = radial(:, i)*2.0_dp*SQRT(a/3.0_dp)*pot%r
         pot%paw%gipaw_core_orbitals%wfc(:, i) = pot%r*radial(:, i)
      END DO
      pot%paw%rho_core(:) = (2.0_dp*radial(:, 1)**2 + 6.0_dp*radial(:, 2)**2)/fourpi
      DO i = 1, nb
         b = valence_exponents(i)
         pot%paw%ae_wfc(:, i) = amplitudes(i)*pot%r*SQRT(fourpi)*(2.0_dp*b/pi)**0.75_dp*EXP(-b*pot%r**2)
         IF (pot%lbeta(i) == 1) pot%paw%ae_wfc(:, i) = pot%paw%ae_wfc(:, i)*2.0_dp*SQRT(b/3.0_dp)*pot%r
         pot%beta(:, i) = pot%r**pot%lbeta(i)*EXP(-pot%r**2)
      END DO
      pot%paw%ps_wfc(:, :) = pot%paw%ae_wfc
      CALL upf_projector_init(pot, 1, projectors)
      CALL paw_onecenter_init(pot, projectors, onecenter)
      CALL paw_onecenter_core_exchange_init(pot, projectors, onecenter)
      IF (.NOT. onecenter%core_exchange_available) ERROR STOP 'Complete Gaussian s/p core exchange is unavailable'
      ! These formulas use normalized Cartesian Gaussians and their Fourier moments.
      ! They do not use production Gaunt coefficients or radial Coulomb quadrature.
      reference = 0.0_dp
      DO j = 1, 2
         DO i = 1, 2
            reference(i, j) = -amplitudes(i)*amplitudes(j)* &
                              (gaussian_cross(core_exponents(1), valence_exponents(i), valence_exponents(j)) + &
                               2.0_dp*core_exponents(2)/(2.0_dp*core_exponents(2) + &
                                                         valence_exponents(i) + valence_exponents(j))* &
                               gaussian_cross(core_exponents(2), valence_exponents(i), valence_exponents(j)))
         END DO
      END DO
      b = valence_exponents(3)
      a = core_exponents(1)
      scale = b/(3.0_dp*(a + b))*gaussian_cross(a, b, b)
      a = core_exponents(2)
      scale = scale + 11.0_dp*b*a/(3.0_dp*(a + b)**2)*gaussian_cross(a, b, b)
      DO i = 3, np
         reference(i, i) = -scale
      END DO
      a = core_exponents(1)
      b = core_exponents(2)
      exact = -(2.0_dp*SQRT(a/pi) + 2.0_dp*b/(a + b)*gaussian_cross(a, b, b) + &
                11.0_dp/4.0_dp*2.0_dp*SQRT(b/pi))
      error = MAX(MAXVAL(ABS(onecenter%core_exchange_matrix - reference)), ABS(onecenter%core_exchange_energy - exact))
      WRITE (*, '(A,ES15.6)') 'PAW original-core Cartesian-Gaussian exchange error:', error
      IF (error > 3.e-8_dp) ERROR STOP 'PAW closed-core exchange differs from independent Cartesian Gaussian integrals'
      DO spin = 1, 2
         DO j = 1, np
            DO i = 1, np
               density(i, j, spin) = SIN((0.19_dp + 0.11_dp*spin)*i)*SIN((0.19_dp + 0.11_dp*spin)*j)
               direction(i, j, spin) = COS((0.23_dp + 0.13_dp*spin)*(i + j))
            END DO
         END DO
      END DO
      DO variant = 1, SIZE(fractions)
         fraction = fractions(variant)
         CALL paw_onecenter_exchange(onecenter, density, fraction, energy, potential)
         IF (ABS(energy - fraction*(exact + SUM(density(:, :, 1)*reference) + &
                                    SUM(density(:, :, 2)*reference))) > 2.e-7_dp) THEN
            ERROR STOP 'PAW core exchange has wrong spin, fraction or linear energy factor'
         END IF
         DO spin = 1, 2
            IF (MAXVAL(ABS(potential(:, :, spin) - fraction*reference)) > 3.e-8_dp) THEN
               ERROR STOP 'PAW core-valence exchange has wrong full-matrix derivative'
            END IF
         END DO
         directional = SUM(potential*direction)
         CALL paw_onecenter_exchange(onecenter, density + step*direction, fraction, plus, scratch)
         CALL paw_onecenter_exchange(onecenter, density - step*direction, fraction, minus, scratch)
         IF (ABS((plus - minus)/(2.0_dp*step) - directional) > 1.e-9_dp*MAX(1.0_dp, ABS(directional))) THEN
            ERROR STOP 'PAW core-valence operator is not the energy derivative'
         END IF
         restricted(:, :, 1) = 0.5_dp*SUM(density, DIM=3)
         restricted(:, :, 2) = restricted(:, :, 1)
         CALL paw_onecenter_exchange(onecenter, restricted, fraction, plus, scratch)
         IF (ABS(plus - energy) > 1.e-13_dp) ERROR STOP 'PAW closed-core RKS and UKS contributions differ'
         CALL paw_onecenter_exchange(onecenter, 0.0_dp*density, fraction, plus, scratch)
         IF (ABS(plus - fraction*exact) > 3.e-8_dp) ERROR STOP 'PAW core-core exchange constant is omitted or repeated'
      END DO
      pot%paw%gipaw_core_orbitals%wfc(:, 2) = -pot%paw%gipaw_core_orbitals%wfc(:, 2)
      phase_changed = onecenter
      CALL paw_onecenter_core_exchange_init(pot, projectors, phase_changed)
      IF (.NOT. phase_changed%core_exchange_available) ERROR STOP 'Original core phase changes core availability'
      IF (MAXVAL(ABS(phase_changed%core_exchange_matrix - onecenter%core_exchange_matrix)) > 1.e-14_dp .OR. &
          ABS(phase_changed%core_exchange_energy - onecenter%core_exchange_energy) > 1.e-14_dp) THEN
         ERROR STOP 'PAW core exchange depends on an arbitrary original core orbital sign'
      END IF
      pot%paw%gipaw_core_orbitals%wfc(:, 2) = 0.9_dp*pot%paw%gipaw_core_orbitals%wfc(:, 2)
      unavailable = onecenter
      CALL paw_onecenter_core_exchange_init(pot, projectors, unavailable)
      IF (unavailable%core_exchange_available .OR. LEN_TRIM(unavailable%core_exchange_reason) == 0) THEN
         ERROR STOP 'Invalid original core data silently accepted for PAW exchange'
      END IF
      CALL paw_onecenter_core_exchange_prepare(pot, projectors, unavailable, 4.1_dp, cache_index)
      ASSOCIATE (core => unavailable%interaction_core_exchange(cache_index))
         IF (core%available .OR. LEN_TRIM(core%reason) == 0) THEN
            ERROR STOP 'Invalid original core data silently accepted for truncated PAW exchange'
         END IF
      END ASSOCIATE
      pot%zion = 10.0_dp
      pot%paw%rho_core = 0.0_dp
      CALL paw_onecenter_core_exchange_prepare(pot, projectors, unavailable, 4.3_dp, cache_index)
      ASSOCIATE (core => unavailable%interaction_core_exchange(cache_index))
         IF (.NOT. core%available) ERROR STOP 'Core-free PAW data unavailable for truncated exchange'
         IF (ANY(core%matrix /= 0.0_dp) .OR. core%energy /= 0.0_dp) THEN
            ERROR STOP 'Core-free truncated PAW exchange is nonzero'
         END IF
      END ASSOCIATE
      CALL atom_release_upf(pot)
   END SUBROUTINE core_exchange_checks

! **************************************************************************************************
!> \brief Compact uniform-sphere core: analytic kernels, radius caches, spin and density derivatives.
!>        Constant s orbitals have quadratic radial products, so nine nodes suffice without
!>        interpolation error. The complete domain extends beyond the augmentation sphere.
! **************************************************************************************************
   SUBROUTINE interaction_core_checks()
      INTEGER, PARAMETER                                 :: n = 9, np = 2
      REAL(dp), PARAMETER :: amplitudes(np) = [1.0_dp, 0.73_dp], &
         cutoffs(4) = [0.0_dp, 0.7_dp, 2.3_dp, 20.0_dp], fractions(3) = [0.0_dp, 0.37_dp, 1.0_dp], &
         omegas(2) = [0.17_dp, 0.61_dp], step = 1.e-5_dp

      INTEGER                                            :: cache_index, i, j, reused, spin, &
                                                            variant, weight
      REAL(dp) :: density(np, np, 2), direction(np, np, 2), energy, error, exact, fraction, minus, &
         plus, potential(np, np, 2), radius, reference(np, np), scratch(np, np, 2), screening
      TYPE(atom_upfpot_type)                             :: pot
      TYPE(paw_onecenter_type)                           :: onecenter, poisoned
      TYPE(upf_projector_type)                           :: projectors

      pot%is_paw = .TRUE.
      pot%symbol = 'Ne'
      pot%zion = 8.0_dp
      pot%number_of_proj = np
      ALLOCATE (pot%paw, pot%r(n), pot%rab(n), pot%lbeta(np), pot%beta(n, np), &
                pot%dion(np, np), pot%qijl(n, np*(np + 1)/2, 0:0))
      ALLOCATE (pot%paw%ae_wfc(n, np), pot%paw%ps_wfc(n, np), pot%paw%rho_core(n))
      ALLOCATE (pot%paw%gipaw_core_orbitals)
      ASSOCIATE (core => pot%paw%gipaw_core_orbitals)
         ALLOCATE (core%wfc(n, 1), core%indices(1), core%principal_quantum_numbers(1), core%angular_momenta(1))
         core%indices = 1
         core%principal_quantum_numbers = 1
         core%angular_momenta = 0
         core%wfc(:, 1) = [(SQRT(3.0_dp)*REAL(i - 1, dp)/REAL(n - 1, dp), i=1, n)]
      END ASSOCIATE
      pot%paw%cutoff_index = 5
      pot%rab = 1.0_dp/REAL(n - 1, dp)
      pot%r(:) = [(REAL(i - 1, dp)*pot%rab(i), i=1, n)]
      pot%lbeta = 0
      pot%dion = 0.0_dp
      pot%qijl = 0.0_dp
      pot%paw%rho_core = 6.0_dp/fourpi
      DO i = 1, np
         pot%paw%ae_wfc(:, i) = amplitudes(i)*pot%paw%gipaw_core_orbitals%wfc(:, 1)
         pot%beta(:, i) = EXP(-pot%r**2)
      END DO
      pot%paw%ps_wfc(:, :) = pot%paw%ae_wfc
      CALL upf_projector_init(pot, 0, projectors)
      CALL paw_onecenter_init(pot, projectors, onecenter)
      CALL paw_onecenter_core_exchange_init(pot, projectors, onecenter)
      IF (.NOT. onecenter%core_exchange_available) ERROR STOP 'Complete uniform-sphere core exchange is unavailable'
      DO spin = 1, 2
         DO j = 1, np
            DO i = 1, np
               density(i, j, spin) = SIN((0.19_dp + 0.11_dp*spin)*i)*SIN((0.19_dp + 0.11_dp*spin)*j)
               direction(i, j, spin) = COS((0.23_dp + 0.13_dp*spin)*(i + j))
            END DO
         END DO
      END DO
      error = 0.0_dp
      DO variant = 1, SIZE(cutoffs) + SIZE(omegas)
         IF (variant <= SIZE(cutoffs)) THEN
            screening = -1.0_dp
            radius = cutoffs(variant)
            CALL paw_onecenter_core_exchange_prepare(pot, projectors, onecenter, radius, cache_index)
            CALL paw_onecenter_core_exchange_prepare(pot, projectors, onecenter, radius, reused)
         ELSE
            radius = 1.0_dp
            screening = omegas(variant - SIZE(cutoffs))
            CALL paw_onecenter_core_exchange_prepare(pot, projectors, onecenter, cache_index=cache_index, omega=screening)
            CALL paw_onecenter_core_exchange_prepare(pot, projectors, onecenter, cache_index=reused, omega=screening)
         END IF
         IF (.NOT. onecenter%interaction_core_exchange(cache_index)%available) THEN
            ERROR STOP 'Complete uniform-sphere core is unavailable for the selected exchange kernel'
         END IF
         IF (reused /= cache_index .OR. SIZE(onecenter%interaction_core_exchange) /= variant) THEN
            ERROR STOP 'PAW core exchange kernel or parameter cache is inconsistent'
         END IF
         reference = 0.0_dp
         exact = -sphere_kernel(radius, screening)
         DO j = 1, np
            DO i = 1, np
               reference(i, j) = amplitudes(i)*amplitudes(j)*exact
            END DO
         END DO
         ASSOCIATE (core => onecenter%interaction_core_exchange(cache_index))
            error = MAX(error, MAXVAL(ABS(core%matrix - reference)), ABS(core%energy - exact))
         END ASSOCIATE
         ! A nonzero valence Coulomb correction must not enter the reconstructed-RI core-only action.
         poisoned = onecenter
         poisoned%hartree = 0.123_dp
         DO weight = 1, SIZE(fractions)
            fraction = fractions(weight)
            CALL paw_onecenter_core_exchange(poisoned, density, fraction, cache_index, energy, potential)
            IF (ABS(energy - fraction*(exact + SUM(density(:, :, 1)*reference) + &
                                       SUM(density(:, :, 2)*reference))) > 2.e-7_dp) THEN
               ERROR STOP 'PAW kernel core-only energy has wrong fraction, spin or duplicate valence terms'
            END IF
            DO spin = 1, 2
               IF (MAXVAL(ABS(potential(:, :, spin) - fraction*reference)) > 3.e-8_dp) THEN
                  ERROR STOP 'PAW kernel core exchange has wrong projector derivative'
               END IF
            END DO
            CALL paw_onecenter_core_exchange(poisoned, density + step*direction, fraction, cache_index, plus, scratch)
            CALL paw_onecenter_core_exchange(poisoned, density - step*direction, fraction, cache_index, minus, scratch)
            IF (ABS((plus - minus)/(2.0_dp*step) - SUM(potential*direction)) > &
                1.e-9_dp*MAX(1.0_dp, ABS(SUM(potential*direction)))) THEN
               ERROR STOP 'PAW kernel core operator is not the density derivative'
            END IF
         END DO
      END DO
     IF (SIZE(onecenter%interaction_core_exchange) /= SIZE(cutoffs) + SIZE(omegas)) ERROR STOP 'PAW kernel core radii share a cache'
      IF (MAXVAL(ABS(onecenter%interaction_core_exchange(4)%matrix - onecenter%core_exchange_matrix)) > 1.e-14_dp .OR. &
          ABS(onecenter%interaction_core_exchange(4)%energy - onecenter%core_exchange_energy) > 1.e-14_dp) THEN
         ERROR STOP 'PAW kernel core preparation replaced the separate Coulomb cache'
      END IF
      WRITE (*, '(A,ES15.6)') 'PAW kernel uniform-sphere exchange error:', error
      IF (error > 3.e-8_dp) ERROR STOP 'PAW kernel core exchange differs from sphere-overlap references'
      CALL atom_release_upf(pot)
   END SUBROUTINE interaction_core_checks

! **************************************************************************************************
!> \brief Crossed Coulomb integral of two normalized s-Gaussian transitions sharing exponent a.
!> \param a Common Gaussian exponent.
!> \param b First valence exponent.
!> \param c Second valence exponent.
!> \return Infinite-space Cartesian-Gaussian Coulomb integral in Hartree.
! **************************************************************************************************
   FUNCTION gaussian_cross(a, b, c) RESULT(value)
      REAL(dp), INTENT(IN)                               :: a, b, c
      REAL(dp)                                           :: value

      REAL(dp)                                           :: first, second

      first = (2.0_dp*SQRT(a*b)/(a + b))**1.5_dp
      second = (2.0_dp*SQRT(a*c)/(a + c))**1.5_dp
      value = first*second*2.0_dp/SQRT(pi)*SQRT((a + b)*(a + c)/(2.0_dp*a + b + c))
   END FUNCTION gaussian_cross

! **************************************************************************************************
!> \brief Independent overlap-volume integral for two normalized uniform unit spheres.
!>        Their overlap at distance d is pi*(4/3-d+d**3/12), for 0<=d<=2.
!> \param radius Truncation radius, capped at the sphere diameter.
!> \param omega Inverse screening length; a negative value selects the truncated reference.
!> \return Bilinear integral of the normalized densities, without a spin or energy factor.
! **************************************************************************************************
   FUNCTION sphere_kernel(radius, omega) RESULT(value)
      REAL(dp), INTENT(IN)                               :: radius, omega
      REAL(dp)                                           :: value

      INTEGER, PARAMETER                                 :: powers(3) = [1, 2, 4]

      INTEGER                                            :: i, m
      REAL(dp)                                           :: a, gaussian(0:5), moment(3), x

      IF (omega < 0.0_dp) THEN
         a = MIN(radius, 2.0_dp)
         value = 1.5_dp*a**2 - 0.75_dp*a**3 + 3.0_dp*a**5/80.0_dp
      ELSE IF (omega == 0.0_dp) THEN
         value = 1.2_dp
      ELSE
         x = 2.0_dp*omega
         gaussian(0) = 0.5_dp*SQRT(pi)*ERF(x)
         gaussian(1) = 0.5_dp*(1.0_dp - EXP(-x*x))
         DO m = 2, 5
            gaussian(m) = 0.5_dp*(REAL(m - 1, dp)*gaussian(m - 2) - x**(m - 1)*EXP(-x*x))
         END DO
         DO i = 1, SIZE(powers)
            m = powers(i) + 1
            moment(i) = (2.0_dp**m*ERFC(x) + 2.0_dp*gaussian(m)/(SQRT(pi)*omega**m))/REAL(m, dp)
         END DO
         value = 2.25_dp*(4.0_dp*moment(1)/3.0_dp - moment(2) + moment(3)/12.0_dp)
      END IF
   END FUNCTION sphere_kernel

! **************************************************************************************************
!> \brief Higher-angular-momentum short-range bilinear against an independent Fourier Gaussian integral.
! **************************************************************************************************
   SUBROUTINE screened_radial_integrals()
      INTEGER, PARAMETER                                 :: l = 8, n = 193
      REAL(dp), PARAMETER :: exponents(2) = [0.7_dp, 1.1_dp], omega = 0.11_dp

      INTEGER                                            :: a, b, i
      REAL(dp)                                           :: error, exact, integral(2, 2), &
                                                            products(n, 2), r(n), reference(2, 2), &
                                                            scale

      DO i = 1, n
         r(i) = 9.0_dp*REAL(i - 1, dp)/REAL(n - 1, dp)
      END DO
      DO a = 1, 2
         products(:, a) = r**(l + 2)*EXP(-exponents(a)*r**2)
      END DO
      CALL paw_radial_short(r, products, l, omega, integral)
      error = 0.0_dp
      DO b = 1, 2
         DO a = 1, 2
            scale = gaussian_coulomb(exponents(a), exponents(b), l)
            exact = scale*(1.0_dp - (omega**2/(exponents(a)*exponents(b)/(exponents(a) + exponents(b)) + &
                                               omega**2))**(REAL(l, dp) + 0.5_dp))
            error = MAX(error, ABS(integral(a, b) - exact)/MAX(1.0_dp, ABS(scale)))
         END DO
      END DO
      IF (integral(1, 2) /= integral(2, 1)) ERROR STOP 'Asymmetric screened PAW radial matrix'
      CALL paw_radial_short(r, products, l, 0.0_dp, integral)
      CALL paw_radial_coulomb(r, products, l, reference)
      IF (ANY(integral /= reference)) ERROR STOP 'Zero-screening PAW kernel differs from Coulomb'
      WRITE (*, '(A,ES15.6)') 'PAW screened Gaussian multipole error:', error
      IF (error > 3.e-6_dp) ERROR STOP 'Screened radial PAW multipole differs from Fourier Gaussian reference'
   END SUBROUTINE screened_radial_integrals

END PROGRAM paw_onecenter_unittest
