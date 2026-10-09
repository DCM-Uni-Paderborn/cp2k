!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
PROGRAM topology_band_unittest
   USE ieee_arithmetic,                 ONLY: ieee_quiet_nan,&
                                              ieee_value
   USE iso_fortran_env,                 ONLY: output_unit
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: twopi
   USE mathlib,                         ONLY: inv_3x3
   USE topology_band_representations,   ONLY: atomic_band_character_polynomial,&
                                              atomic_band_characters,&
                                              atomic_band_construct,&
                                              atomic_band_site_data,&
                                              atomic_band_site_match,&
                                              atomic_band_type,&
                                              atomic_character_compare,&
                                              atomic_character_polynomial_type
   USE topology_compatibility,          ONLY: little_group_data_type
   USE topology_corepresentations,      ONLY: corepresentation_multiplicities,&
                                              corepresentation_table_type,&
                                              projective_corepresentations
   USE topology_little_group,           ONLY: little_group_factor_system,&
                                              projective_characters
   USE topology_symmetry,               ONLY: character_multiplicities,&
                                              spin_rotation

   IMPLICIT NONE
   INTEGER :: io_unit
   INTEGER :: polynomial_bands = 0, polynomial_queries = 0
   REAL(dp) :: polynomial_residual = 0

   io_unit = output_unit
   CALL analytic_tests()
   CALL coordinate_tests()
   CALL polynomial_tests()
   CALL site_reference_tests()
   CALL screw_site_reference_tests()
   CALL magnetic_site_reference_tests()
   WRITE (io_unit, *) 'Fourier character bands/queries/largest residual:', &
      polynomial_bands, polynomial_queries, polynomial_residual
   WRITE (io_unit, *) 'Site-induced atomic band tests passed.'
CONTAINS

! **************************************************************************************************
!> \brief Compare quotient Fourier sums to raw induction on every torsion component and free samples.
!> \param band Prepared site-induced representations
!> \param polynomial Constructed Fourier form for further comparisons
! **************************************************************************************************
   SUBROUTINE check_polynomial(band, polynomial)
      TYPE(atomic_band_type), INTENT(IN)                 :: band
      TYPE(atomic_character_polynomial_type), &
         INTENT(OUT)                                     :: polynomial

      COMPLEX(dp), ALLOCATABLE                           :: chi(:, :), collected(:)
      INTEGER                                            :: g, i, j, m1, m2, m3, periods(3), row, &
                                                            sample, status
      INTEGER, ALLOCATABLE                               :: mapping(:)
      REAL(dp)                                           :: angle, error, k(3), y(3)

      CALL atomic_band_character_polynomial(band, polynomial, status)
      IF (status /= 0) ERROR STOP 'Fourier character construction failed'
      polynomial_bands = polynomial_bands + 1
      ALLOCATE (collected(SIZE(polynomial%dimensions)))
      DO g = 1, SIZE(polynomial%blocks)
         ASSOCIATE (block => polynomial%blocks(g))
            IF (ANY(block%period > 12)) ERROR STOP 'Unexpected crystallographic Smith period'
            periods = MAX(1, INT(block%period))
            DO m1 = 0, periods(1) - 1
            DO m2 = 0, periods(2) - 1
            DO m3 = 0, periods(3) - 1
            DO sample = 1, MERGE(3, 1, ANY(block%period == 0))
               y = REAL([m1, m2, m3], dp)/periods
               DO j = 1, 3
                  IF (block%period(j) == 0) y(j) = MODULO(SQRT(REAL(j + 1, dp))*(sample - 1), 1.0_dp)
               END DO
               k = MATMUL(TRANSPOSE(REAL(block%left, dp)), y)
               CALL atomic_band_characters(band, k, mapping, chi, 1.e-9_dp, status)
               IF (status /= 0) ERROR STOP 'Fourier raw character evaluation failed'
               row = FINDLOC(mapping, block%operation, DIM=1)
               IF (row == 0) ERROR STOP 'Smith dual sample not fixed by operation'
               collected = 0
               DO i = 1, SIZE(block%keys, 2)
                  angle = -twopi*DOT_PRODUCT(y, REAL(block%keys(:, i), dp))
                  collected = collected + CMPLX(COS(angle), SIN(angle), dp)*block%coefficients(i, :)
               END DO
               error = MAXVAL(ABS(collected - chi(row, :))/polynomial%dimensions)
               polynomial_residual = MAX(polynomial_residual, error)
               polynomial_queries = polynomial_queries + 1
               IF (error > 1.e-8_dp) ERROR STOP 'Quotient Fourier sum disagrees with raw Bloch induction'
            END DO
            END DO
            END DO
            END DO
         END ASSOCIATE
      END DO
   END SUBROUTINE check_polynomial

! **************************************************************************************************
!> \brief Inversion centers indistinguishable at Gamma differ on disconnected invariant components.
! **************************************************************************************************
   SUBROUTINE polynomial_tests()
      COMPLEX(dp), ALLOCATABLE                           :: gamma_a(:, :), gamma_b(:, :)
      INTEGER                                            :: i, j, n, r(3, 3, 2), status
      INTEGER, ALLOCATABLE                               :: mapping(:)
      LOGICAL                                            :: equal
      REAL(dp)                                           :: cell(3, 3), error, tau(3, 2)
      TYPE(atomic_band_type)                             :: band, empty
      TYPE(atomic_character_polynomial_type)             :: a, b, bad

      cell = 0
      r = 0
      tau = 0
      DO i = 1, 3
         cell(i, i) = 1
         r(i, i, 1) = 1
         r(i, i, 2) = -1
      END DO
      CALL atomic_band_construct(cell, r, tau, [.FALSE., .FALSE.], [0.0_dp, 0.0_dp, 0.0_dp], &
                                 .FALSE., band, 1.e-9_dp, error, status)
      IF (status /= 0) ERROR STOP 'Fourier inversion reference construction failed'
      CALL check_polynomial(band, a)
      CALL atomic_band_characters(band, [0.0_dp, 0.0_dp, 0.0_dp], mapping, gamma_a, 1.e-9_dp, status)
      IF (status /= 0) ERROR STOP 'Fourier inversion Gamma reference failed'
      CALL atomic_band_construct(cell, r, tau, [.FALSE., .FALSE.], [0.5_dp, 0.0_dp, 0.0_dp], &
                                 .FALSE., band, 1.e-9_dp, error, status)
      IF (status /= 0) ERROR STOP 'Fourier shifted inversion construction failed'
      CALL check_polynomial(band, b)
      CALL atomic_band_characters(band, [0.0_dp, 0.0_dp, 0.0_dp], mapping, gamma_b, 1.e-9_dp, status)
      IF (status /= 0) ERROR STOP 'Fourier shifted inversion Gamma failed'
      n = 0
      DO i = 1, SIZE(a%dimensions)
         CALL atomic_character_compare(a, i, a, i, equal, 1.e-9_dp, error, status)
         IF (status /= 0 .OR. .NOT. equal .OR. error /= 0) ERROR STOP 'Fourier self comparison failed'
         DO j = 1, SIZE(b%dimensions)
            IF (MAXVAL(ABS(gamma_a(:, i) - gamma_b(:, j))) > 1.e-9_dp) CYCLE
            n = n + 1
            CALL atomic_character_compare(a, i, b, j, equal, 1.e-9_dp, error, status)
            IF (status /= 0 .OR. equal) ERROR STOP 'Fourier comparison retained sampled inversion alias'
         END DO
      END DO
      IF (n /= 2) ERROR STOP 'Fourier inversion test lacks Gamma aliases'
      bad = a
      bad%blocks(1)%coefficients(1, 1) = CMPLX(ieee_value(0.0_dp, ieee_quiet_nan), 0.0_dp, dp)
      CALL atomic_character_compare(a, 1, bad, 1, equal, 1.e-9_dp, error, status)
      IF (status == 0 .OR. equal) ERROR STOP 'Fourier comparison accepted nonfinite coefficients'
      bad = a
      DEALLOCATE (bad%blocks(1)%keys)
      CALL atomic_character_compare(a, 1, bad, 1, equal, 1.e-9_dp, error, status)
      IF (status == 0 .OR. equal) ERROR STOP 'Fourier comparison accepted absent keys'
      CALL atomic_band_character_polynomial(empty, bad, status)
      IF (status == 0 .OR. ALLOCATED(bad%blocks)) ERROR STOP 'Fourier construction accepted empty band'
      WRITE (io_unit, *) 'Fourier invariant-torus alias and invalid-data tests passed.'
   END SUBROUTINE polynomial_tests

! **************************************************************************************************
!> \brief Distinguish equal-dimensional onsite C3 irreps with explicit spin and antiunitary conventions.
! **************************************************************************************************
   SUBROUTINE site_reference_tests()
      COMPLEX(dp)                                        :: chars(3, 3), gauged_spin(2, 2, 6), &
                                                            spin(2, 2, 6), theta(2, 2), z
      COMPLEX(dp), ALLOCATABLE                           :: bad(:, :), reference(:, :)
      INTEGER                                            :: g, i, inverse(3, 3), j, n, ni, order(6), &
                                                            r(3, 3, 6), s, shear(3, 3), status, tr
      INTEGER, ALLOCATABLE                               :: mapping(:), permuted(:)
      LOGICAL                                            :: anti(6)
      REAL(dp)                                           :: angle, bad_cart(3, 3, 6), cart(3, 3, 6), &
                                                            cell(3, 3), error, orient(3, 3), &
                                                            shift(3), tau(3, 6), transformed(3, 3)
      TYPE(atomic_band_type)                             :: band, empty

      cell = 0
      cell(:, 1) = [1.0_dp, 0.0_dp, 0.0_dp]
      cell(:, 2) = [-0.5_dp, SQRT(3.0_dp)/2, 0.0_dp]
      cell(3, 3) = 1
      r = 0
      DO i = 1, 3
         r(i, i, 1) = 1
      END DO
      r(:, :, 2) = RESHAPE([0, 1, 0, -1, -1, 0, 0, 0, 1], [3, 3])
      r(:, :, 3) = MATMUL(r(:, :, 2), r(:, :, 2))
      r(:, :, 4:) = r(:, :, :3)
      DO g = 1, 6
         cart(:, :, g) = MATMUL(MATMUL(cell, REAL(r(:, :, g), dp)), inv_3x3(cell))
      END DO
      angle = 0.73_dp
      orient = RESHAPE([COS(angle), 0.0_dp, -SIN(angle), 0.0_dp, 1.0_dp, 0.0_dp, &
                        SIN(angle), 0.0_dp, COS(angle)], [3, 3])
      shear = RESHAPE([1, 0, 0, 1, 1, 0, 0, 1, 1], [3, 3])
      inverse = NINT(inv_3x3(REAL(shear, dp)))
      transformed = MATMUL(MATMUL(orient, cell), REAL(shear, dp))
      shift = [0.137_dp, 0.219_dp, 0.317_dp]
      DO g = 1, 6
         r(:, :, g) = MATMUL(MATMUL(inverse, r(:, :, g)), shear)
         tau(:, g) = shift - MATMUL(REAL(r(:, :, g), dp), shift)
      END DO
      theta = 0
      theta(1, 2) = 1
      theta(2, 1) = -1
      DO s = 0, 1
         spin = 0
         DO g = 1, 3
            z = EXP(CMPLX(0.0_dp, -ACOS(-1.0_dp)*s*(g - 1)/3, dp))
            spin(1, 1, g) = z
            spin(2, 2, g) = CONJG(z)
            spin(:, :, g + 3) = MATMUL(spin(:, :, g), theta)
            DO j = 1, 3
               chars(g, j) = EXP(CMPLX(0.0_dp, ACOS(-1.0_dp)*(g - 1)*(2*(j - 1) + s)/3, dp))
            END DO
         END DO
         DO tr = 0, 1
            n = 3*(1 + tr)
            ni = 3 - tr
            anti = .FALSE.
            IF (tr == 1) anti(4:) = .TRUE.
            order = [1, 3, 2, 4, 6, 5]
            IF (tr == 1) order = [1, 4, 3, 6, 2, 5]
            ALLOCATE (reference(3, ni), bad(3, ni))
            IF (tr == 0) THEN
               reference(:, :) = chars
            ELSE IF (s == 0) THEN
               reference(:, 1) = chars(:, 1)
               reference(:, 2) = chars(:, 2) + chars(:, 3)
            ELSE
               reference(:, 1) = chars(:, 1) + chars(:, 3)
               reference(:, 2) = 2*chars(:, 2)
            END IF
            CALL atomic_band_construct(transformed, r(:, :, order(:n)), tau(:, order(:n)), anti(order(:n)), shift, &
                                       s == 1, band, 1.e-10_dp, error, status)
            IF (status /= 0) ERROR STOP 'Onsite C3 reference construction failed'
            CALL atomic_band_site_match(band, cart(:, :, :n), anti(:n), spin(:, :, :n), reference, orient, &
                                        mapping, 1.e-10_dp, error, status)
            IF (status /= 0) ERROR STOP 'Onsite C3 axes/spin/corepresentation matching failed'
            CALL atomic_band_site_match(band, cart(:, :, :n), anti(:n), spin(:, :, :n), reference(:, ni:1:-1), orient, &
                                        permuted, 1.e-10_dp, error, status)
            IF (status /= 0) ERROR STOP 'Permuted onsite reference table rejected'
            IF (ANY(permuted /= ni + 1 - mapping)) ERROR STOP 'Equal-dimensional onsite labels not distinguished'
            CALL atomic_band_site_match(band, cart(:, :, :n), anti(:n), spin(:, :, :n), reference, -orient, &
                                        permuted, 1.e-10_dp, error, status)
            IF (status /= 0) ERROR STOP 'Improper Cartesian reference axes rejected'
            IF (ANY(permuted /= mapping)) ERROR STOP 'Improper axes changed onsite labels'
            IF (s == 1) THEN
               gauged_spin = spin
               bad(:, :) = reference
               DO g = 2, n
                  z = EXP(CMPLX(0.0_dp, 0.219_dp*g, dp))
                  gauged_spin(:, :, g) = z*spin(:, :, g)
                  IF (g <= 3) bad(g, :) = z*reference(g, :)
               END DO
               CALL atomic_band_site_match(band, cart(:, :, :n), anti(:n), gauged_spin(:, :, :n), bad, orient, &
                                           permuted, 1.e-10_dp, error, status)
               IF (status /= 0) ERROR STOP 'Explicit onsite lift regauging rejected'
               IF (ANY(permuted /= mapping)) ERROR STOP 'Explicit lift gauge changed onsite labels'
               CALL atomic_band_site_match(band, cart(:, :, :n), anti(:n), gauged_spin(:, :, :n), reference, orient, &
                                           permuted, 1.e-10_dp, error, status)
               IF (status == 0 .OR. ALLOCATED(permuted)) ERROR STOP 'Unmatched onsite lift/character gauges accepted'
               gauged_spin = spin
               gauged_spin(:, :, 2) = TRANSPOSE(theta)
               CALL atomic_band_site_match(band, cart(:, :, :n), anti(:n), gauged_spin(:, :, :n), reference, orient, &
                                           permuted, 1.e-10_dp, error, status)
               IF (status == 0 .OR. ALLOCATED(permuted)) ERROR STOP 'Wrong onsite spin operation accepted'
               IF (tr == 1) THEN
                  gauged_spin = spin
                  gauged_spin(:, :, 4) = spin(:, :, 1)
                  CALL atomic_band_site_match(band, cart(:, :, :n), anti(:n), gauged_spin(:, :, :n), reference, orient, &
                                              permuted, 1.e-10_dp, error, status)
                  IF (status == 0 .OR. ALLOCATED(permuted)) ERROR STOP 'Missing physical Theta spin action accepted'
               END IF
            END IF
            bad_cart = cart
            bad_cart(:, :, 2) = cart(:, :, 1)
            CALL atomic_band_site_match(band, bad_cart(:, :, :n), anti(:n), spin(:, :, :n), reference, orient, &
                                        permuted, 1.e-10_dp, error, status)
            IF (status == 0 .OR. ALLOCATED(permuted)) ERROR STOP 'Duplicate Cartesian reference operation accepted'
            bad_cart(:, :, 1) = 0
            DO i = 1, 3
               bad_cart(i, i, 1) = 1
            END DO
            CALL atomic_band_site_match(band, cart(:, :, :n), anti(:n), spin(:, :, :n), reference, bad_cart(:, :, 1), &
                                        permuted, 1.e-10_dp, error, status)
            IF (status == 0 .OR. ALLOCATED(permuted)) ERROR STOP 'Wrong but orthogonal reference axes accepted'
            bad_cart(:, :, 1) = HUGE(1.0_dp)
            CALL atomic_band_site_match(band, cart(:, :, :n), anti(:n), spin(:, :, :n), reference, bad_cart(:, :, 1), &
                                        permuted, 1.e-10_dp, error, status)
            IF (status == 0 .OR. ALLOCATED(permuted)) ERROR STOP 'Overflowing reference axes accepted'
            gauged_spin = CMPLX(HUGE(1.0_dp), HUGE(1.0_dp), dp)
            CALL atomic_band_site_match(band, cart(:, :, :n), anti(:n), gauged_spin(:, :, :n), reference, orient, &
                                        permuted, 1.e-10_dp, error, status)
            IF (status == 0 .OR. ALLOCATED(permuted)) ERROR STOP 'Overflowing reference spin accepted'
            bad(:, :) = reference
            bad(:, 1) = reference(:, 2)
            CALL atomic_band_site_match(band, cart(:, :, :n), anti(:n), spin(:, :, :n), bad, orient, &
                                        mapping, 1.e-10_dp, error, status)
            IF (status == 0 .OR. ALLOCATED(mapping)) ERROR STOP 'Repeated onsite reference accepted'
            bad(:, :) = reference
            bad(2, 1) = ieee_value(0.0_dp, ieee_quiet_nan)
            CALL atomic_band_site_match(band, cart(:, :, :n), anti(:n), spin(:, :, :n), bad, orient, &
                                        mapping, 1.e-10_dp, error, status)
            IF (status == 0 .OR. ALLOCATED(mapping)) ERROR STOP 'Nonfinite onsite reference accepted'
            CALL atomic_band_site_match(band, cart(:, :, :n), anti(:n), spin(:, :, :n), reference, 0.0_dp*orient, &
                                        mapping, 1.e-10_dp, error, status)
            IF (status == 0) ERROR STOP 'Invalid onsite axes accepted'
            CALL atomic_band_site_match(empty, cart(:, :, :n), anti(:n), spin(:, :, :n), reference, orient, &
                                        mapping, 1.e-10_dp, error, status)
            IF (status == 0) ERROR STOP 'Unconstructed onsite representation accepted'
            DEALLOCATE (reference, bad)
         END DO
      END DO
      WRITE (io_unit, '(A)') 'Onsite reference axes, spin lifts and paired/doubled characters: passed.'
   END SUBROUTINE site_reference_tests

! **************************************************************************************************
!> \brief P6_3 has a two-site orbit with C3 onsite symmetry; retain its labels and screw Bloch phases.
! **************************************************************************************************
   SUBROUTINE screw_site_reference_tests()
      COMPLEX(dp)                                        :: chars(3, 3), lift(2, 2), phase, &
                                                            ref_spin(2, 2, 6), &
                                                            rotation_spin(2, 2), theta(2, 2), z
      COMPLEX(dp), ALLOCATABLE                           :: induced(:, :), onsite(:, :), &
                                                            reference(:, :)
      INTEGER                                            :: g, i, inverse(3, 3), j, n, &
                                                            native(3, 3, 12), nh, ni, order(12), &
                                                            p, query, r(3, 3, 12), s, shear(3, 3), &
                                                            status, tr
      INTEGER, ALLOCATABLE                               :: labels(:), mapping(:), site_mapping(:)
      LOGICAL                                            :: anti(12), site_anti(6)
      REAL(dp)                                           :: angle, axes(3, 3), cart(3, 3, 6), &
                                                            cell(3, 3), error, k(3), orient(3, 3), &
                                                            shift(3), tau(3, 12)
      REAL(dp), ALLOCATABLE                              :: positions(:, :)
      TYPE(atomic_band_type)                             :: band

      cell = 0
      cell(:, 1) = [1.0_dp, 0.0_dp, 0.0_dp]
      cell(:, 2) = [-0.5_dp, SQRT(3.0_dp)/2, 0.0_dp]
      cell(3, 3) = 1
      r = 0
      DO i = 1, 3
         r(i, i, 1) = 1
      END DO
      r(:, :, 2) = RESHAPE([1, 1, 0, -1, 0, 0, 0, 0, 1], [3, 3])
      DO g = 3, 6
         r(:, :, g) = MATMUL(r(:, :, g - 1), r(:, :, 2))
      END DO
      r(:, :, 7:) = r(:, :, :6)
      DO g = 1, 3
         cart(:, :, g) = MATMUL(MATMUL(cell, REAL(r(:, :, 2*g - 1), dp)), inv_3x3(cell))
      END DO
      cart(:, :, 4:) = cart(:, :, :3)
      angle = 0.73_dp
      orient = RESHAPE([COS(angle), 0.0_dp, -SIN(angle), 0.0_dp, 1.0_dp, 0.0_dp, &
                        SIN(angle), 0.0_dp, COS(angle)], [3, 3])
      CALL spin_rotation(orient, rotation_spin, 1.e-10_dp, status)
      IF (status /= 0) ERROR STOP 'Invalid screw-site reference axes'
      shear = RESHAPE([1, 0, 0, 1, 1, 0, 0, 1, 1], [3, 3])
      inverse = NINT(inv_3x3(REAL(shear, dp)))
      axes = MATMUL(MATMUL(orient, cell), REAL(shear, dp))
      shift = [0.137_dp, 0.219_dp, 0.317_dp]
      DO g = 1, 12
         native(:, :, g) = MATMUL(MATMUL(inverse, r(:, :, g)), shear)
         tau(:, g) = shift - MATMUL(REAL(native(:, :, g), dp), shift) + &
                     MATMUL(REAL(inverse, dp), [0.0_dp, 0.0_dp, 0.5_dp*MOD(g - 1, 6)])
      END DO
      theta = 0
      theta(1, 2) = 1
      theta(2, 1) = -1
      DO s = 0, 1
         ref_spin = 0
         DO g = 1, 3
            z = EXP(CMPLX(0.0_dp, -ACOS(-1.0_dp)*s*(g - 1)/3, dp))
            ref_spin(1, 1, g) = z
            ref_spin(2, 2, g) = CONJG(z)
            ref_spin(:, :, g + 3) = MATMUL(ref_spin(:, :, g), theta)
            DO j = 1, 3
               chars(g, j) = EXP(CMPLX(0.0_dp, ACOS(-1.0_dp)*(g - 1)*(2*(j - 1) + s)/3, dp))
            END DO
         END DO
         DO tr = 0, 1
            n = 6*(tr + 1)
            nh = 3*(tr + 1)
            ni = 3 - tr
            anti = .FALSE.
            site_anti = .FALSE.
            IF (tr == 1) THEN
               anti(7:) = .TRUE.
               site_anti(4:) = .TRUE.
            END IF
            order(:n) = [(n + 1 - g, g=1, n)]
            ALLOCATE (reference(3, ni))
            IF (tr == 0) THEN
               reference(:, :) = chars
            ELSE IF (s == 0) THEN
               reference(:, 1) = chars(:, 1)
               reference(:, 2) = chars(:, 2) + chars(:, 3)
            ELSE
               reference(:, 1) = chars(:, 1) + chars(:, 3)
               reference(:, 2) = 2*chars(:, 2)
            END IF
            CALL atomic_band_construct(axes, native(:, :, order(:n)), tau(:, order(:n)), anti(order(:n)), shift, &
                                       s == 1, band, 1.e-10_dp, error, status)
            IF (status /= 0) ERROR STOP 'P6_3 onsite reference construction failed'
            CALL atomic_band_site_data(band, positions, site_mapping, onsite)
            IF (SIZE(positions, 2) /= 2 .OR. SIZE(onsite, 1) /= 3) ERROR STOP 'P6_3 C3 stabilizer/orbit incorrect'
            CALL atomic_band_site_match(band, cart(:, :, :nh), site_anti(:nh), ref_spin(:, :, :nh), reference, orient, &
                                        labels, 1.e-10_dp, error, status)
            IF (status /= 0) ERROR STOP 'P6_3 C3 onsite labels not matched'
            DO query = 0, 2
               k = MATMUL(TRANSPOSE(REAL(shear, dp)), [0.0_dp, 0.0_dp, 0.25_dp*query])
               CALL atomic_band_characters(band, k, mapping, induced, 1.e-10_dp, status)
               IF (status /= 0) ERROR STOP 'P6_3 induced characters failed'
               DO i = 1, SIZE(mapping)
                  g = order(mapping(i))
                  IF (MOD(g, 2) == 0) THEN
                     IF (MAXVAL(ABS(induced(i, :))) > 1.e-10_dp) ERROR STOP 'Screw swapping sites has nonzero trace'
                     CYCLE
                  END IF
                  p = (g - 1)/2
                  phase = EXP(CMPLX(0.0_dp, -twopi*0.25_dp*query*p, dp))
                  IF (s == 1) THEN
                     CALL spin_rotation(MATMUL(MATMUL(orient, cart(:, :, p + 1)), TRANSPOSE(orient)), &
                                        lift, 1.e-10_dp, status)
                     IF (status /= 0) ERROR STOP 'Invalid induced spin convention'
                     lift = MATMUL(CONJG(TRANSPOSE(rotation_spin)), MATMUL(lift, rotation_spin))
                     phase = phase*SUM(CONJG(ref_spin(:, :, p + 1))*lift)/2
                  END IF
                  IF (MAXVAL(ABS(induced(i, :) - 2*phase*reference(p + 1, labels))) > 1.e-10_dp) THEN
                     ERROR STOP 'Named C3 onsite columns lost their induced screw phase'
                  END IF
               END DO
            END DO
            DEALLOCATE (reference)
         END DO
      END DO
      WRITE (io_unit, '(A)') 'P6_3 two-site C3 labels, scalar/spin/grey screw-phase induction: passed.'
   END SUBROUTINE screw_site_reference_tests

! **************************************************************************************************
!> \brief Distinguish an onsite antiunitary C2*T operation from a magnetic translation exchanging sites.
! **************************************************************************************************
   SUBROUTINE magnetic_site_reference_tests()
      COMPLEX(dp)                                        :: reference(1, 1), rotation(2, 2), &
                                                            spin(2, 2, 2), theta(2, 2)
      COMPLEX(dp), ALLOCATABLE                           :: characters(:, :)
      INTEGER                                            :: case, i, n, r(3, 3, 2), s, status
      INTEGER, ALLOCATABLE                               :: mapping(:), site_mapping(:)
      LOGICAL                                            :: anti(2), ref_anti(2)
      REAL(dp)                                           :: cart(3, 3, 2), cell(3, 3), error, &
                                                            tau(3, 2)
      REAL(dp), ALLOCATABLE                              :: positions(:, :)
      TYPE(atomic_band_type)                             :: band

      cell = 0
      spin = 0
      DO i = 1, 3
         cell(i, i) = 1
      END DO
      DO i = 1, 2
         spin(i, i, 1) = 1
      END DO
      theta = 0
      theta(1, 2) = 1
      theta(2, 1) = -1
      reference = 1
      anti = [.FALSE., .TRUE.]
      DO CASE = 1, 2
         r(:, :, 1) = NINT(cell)
         r(:, :, 2) = NINT(cell)
         tau = 0
         IF (CASE == 1) THEN
            r(1, 1, 2) = -1
            r(2, 2, 2) = -1
            n = 2
         ELSE
            tau(3, 2) = 0.5_dp
            n = 1
         END IF
         cart = REAL(r, dp)
         ref_anti = anti
         CALL spin_rotation(cart(:, :, 2), rotation, 1.e-10_dp, status)
         IF (status /= 0) ERROR STOP 'Invalid magnetic onsite reference rotation'
         spin(:, :, 2) = MATMUL(rotation, theta)
         DO s = 0, 1
            CALL atomic_band_construct(cell, r, tau, anti, [0.0_dp, 0.0_dp, 0.0_dp], s == 1, &
                                       band, 1.e-10_dp, error, status)
            IF (status /= 0) ERROR STOP 'Magnetic onsite reference construction failed'
            CALL atomic_band_site_data(band, positions, site_mapping, characters)
            IF (SIZE(positions, 2) /= 3 - n) ERROR STOP 'Magnetic onsite orbit incorrect'
            CALL atomic_band_site_match(band, cart(:, :, :n), ref_anti(:n), spin(:, :, :n), reference, cell, &
                                        mapping, 1.e-10_dp, error, status)
            IF (status /= 0) ERROR STOP 'Magnetic onsite reference matching failed'
            IF (SIZE(mapping) /= 1 .OR. mapping(1) /= 1) ERROR STOP 'Magnetic onsite label incorrect'
         END DO
      END DO
      WRITE (io_unit, '(A)') 'Magnetic C2*T stabilizer and site-exchanging half translation: passed.'
   END SUBROUTINE magnetic_site_reference_tests

! **************************************************************************************************
!> \brief Integer coset changes, passive origin shifts and nonorthogonal unimodular cell changes.
! **************************************************************************************************
   SUBROUTINE coordinate_tests()
      COMPLEX(dp), ALLOCATABLE                           :: chi(:, :), reference(:, :)
      INTEGER                                            :: found, g, i, info, inverse(3, 3), j, q, &
                                                            rotations(3, 3, 6), s, shear(3, 3), &
                                                            transformed_rotations(3, 3, 6), variant
      INTEGER, ALLOCATABLE                               :: mapping(:), transformed_mapping(:)
      LOGICAL                                            :: anti(6)
      LOGICAL, ALLOCATABLE                               :: used(:)
      REAL(dp) :: angle, cell(3, 3), error, k(3), origin(3), seed(3), shifts(3, 6), tau(3, 6), &
         transformed_cell(3, 3), transformed_k(3), transformed_seed(3), transformed_tau(3, 6)
      TYPE(atomic_band_type)                             :: band, changed

      cell = 0.0_dp
      rotations = 0
      DO i = 1, 3
         cell(i, i) = 1.0_dp
         rotations(i, i, 1) = 1
      END DO
      rotations(:, :, 2) = RESHAPE([0, 1, 0, 0, 0, 1, 1, 0, 0], [3, 3])
      rotations(:, :, 3) = MATMUL(rotations(:, :, 2), rotations(:, :, 2))
      rotations(:, :, 4:6) = rotations(:, :, 1:3)
      tau(:, :3) = 0.0_dp
      tau(:, 4:6) = 0.5_dp
      anti(:3) = .FALSE.
      anti(4:6) = .TRUE.
      seed = 0.125_dp
      origin = [0.123_dp, -0.217_dp, 0.319_dp]
      shear = RESHAPE([1, 0, 0, 1, 1, 0, 0, 1, 1], [3, 3])
      inverse = RESHAPE([1, 0, 0, -1, 1, 0, 1, -1, 1], [3, 3])
      shifts = 0.0_dp
      shifts(:, 2) = [1.0_dp, 0.0_dp, 0.0_dp]
      shifts(:, 3) = [0.0_dp, -1.0_dp, 2.0_dp]
      shifts(:, 5) = [1.0_dp, 2.0_dp, -1.0_dp]
      DO s = 0, 1
         CALL atomic_band_construct(cell, rotations, tau, anti, seed, s == 1, band, 1.e-9_dp, error, info)
         IF (info /= 0) ERROR STOP 'Coordinate reference site failed'
         DO variant = 1, 3
            transformed_cell = cell
            transformed_rotations = rotations
            transformed_tau = tau
            transformed_seed = seed
            SELECT CASE (variant)
            CASE (1)
               transformed_tau = tau + shifts
               transformed_seed = seed + [2.0_dp, -3.0_dp, 1.0_dp]
            CASE (2)
               transformed_seed = seed - origin
               DO g = 1, 6
                  transformed_tau(:, g) = tau(:, g) + MATMUL(REAL(rotations(:, :, g), dp), origin) - origin
               END DO
            CASE (3)
               transformed_cell = MATMUL(cell, REAL(shear, dp))
               transformed_seed = MATMUL(REAL(inverse, dp), seed)
               transformed_tau = MATMUL(REAL(inverse, dp), tau)
               DO g = 1, 6
                  transformed_rotations(:, :, g) = MATMUL(MATMUL(inverse, rotations(:, :, g)), shear)
               END DO
            END SELECT
            CALL atomic_band_construct(transformed_cell, transformed_rotations, transformed_tau, anti, transformed_seed, &
                                       s == 1, changed, 1.e-9_dp, error, info)
            IF (info /= 0) ERROR STOP 'Transformed site induction failed'
            DO q = 0, 2
               k = 0.25_dp*q
               transformed_k = k
               IF (variant == 3) transformed_k = MATMUL(REAL(TRANSPOSE(shear), dp), k)
               CALL atomic_band_characters(band, k, mapping, reference, 1.e-9_dp, info)
               IF (info /= 0) ERROR STOP 'Coordinate reference characters failed'
               CALL atomic_band_characters(changed, transformed_k, transformed_mapping, chi, 1.e-9_dp, info)
               IF (info /= 0) ERROR STOP 'Transformed characters failed'
               IF (ANY(SHAPE(chi) /= SHAPE(reference))) ERROR STOP 'Coordinate change altered dimensions'
               IF (ANY(mapping /= transformed_mapping)) ERROR STOP 'Coordinate change altered little group'
               IF (variant == 1) THEN
                  DO g = 1, SIZE(mapping)
                     angle = -twopi*DOT_PRODUCT(k, shifts(:, mapping(g)))
                     reference(g, :) = CMPLX(COS(angle), SIN(angle), dp)*reference(g, :)
                  END DO
               END IF
               ALLOCATE (used(SIZE(chi, 2)))
               used = .FALSE.
               DO i = 1, SIZE(reference, 2)
                  found = 0
                  DO j = 1, SIZE(chi, 2)
                     IF (used(j)) CYCLE
                     IF (MAXVAL(ABS(chi(:, j) - reference(:, i))) > 1.e-8_dp) CYCLE
                     found = j
                     EXIT
                  END DO
                  IF (found == 0) ERROR STOP 'Site induction depends on cell/origin/coset convention'
                  used(found) = .TRUE.
               END DO
               DEALLOCATE (used)
            END DO
         END DO
      END DO
      WRITE (io_unit, *) 'Atomic-band origin, cell and coset gauge checks passed.'
   END SUBROUTINE coordinate_tests

! **************************************************************************************************
!> \brief Analytic inversion, screw, quaternion, magnetic cosets, and malformed inputs.
! **************************************************************************************************
   SUBROUTINE analytic_tests()
      COMPLEX(dp), ALLOCATABLE                           :: chi(:, :), onsite(:, :), reference(:, :)
      INTEGER                                            :: a, i, info, j, n, perm(6), q, &
                                                            rotations(3, 3, 8), s
      INTEGER, ALLOCATABLE                               :: alternate(:), counts(:, :), mapping(:), &
                                                            site_map(:)
      LOGICAL                                            :: anti(8)
      REAL(dp)                                           :: cell(3, 3), error, k(3), parity, &
                                                            seed(3), tau(3, 8)
      REAL(dp), ALLOCATABLE                              :: positions(:, :)
      TYPE(atomic_band_type)                             :: band, other
      TYPE(little_group_data_type)                       :: point

      cell = 0.0_dp
      rotations = 0
      DO i = 1, 3
         cell(i, i) = 1.0_dp
         rotations(i, i, 1) = 1
      END DO
      rotations(:, :, 2) = -rotations(:, :, 1)
      tau = 0.0_dp
      anti = .FALSE.
      DO s = 0, 1
         rotations(:, :, 3:4) = rotations(:, :, 1:2)
         anti(3:4) = .TRUE.
         n = 2 + 2*s
         DO a = 0, 7
            DO i = 1, 3
               seed(i) = 0.5_dp*IBITS(a, i - 1, 1)
            END DO
            CALL atomic_band_construct(cell, rotations(:, :, :n), tau(:, :n), anti(:n), seed, s == 1, &
                                       band, 1.e-9_dp, error, info)
            IF (info /= 0) ERROR STOP 'Inversion site construction failed'
            CALL atomic_band_site_data(band, positions, site_map, onsite)
            IF (SIZE(positions, 2) /= 1 .OR. SIZE(onsite, 2) /= 2) ERROR STOP 'Inversion site rank failed'
            DO q = 0, 7
               DO i = 1, 3
                  k(i) = 0.5_dp*IBITS(q, i - 1, 1)
               END DO
               CALL atomic_band_characters(band, k, mapping, chi, 1.e-9_dp, info)
               IF (info /= 0) ERROR STOP 'Inversion band evaluation failed'
               parity = REAL((-1)**POPCNT(IAND(a, q)), dp)
               IF (MAXVAL(ABS(chi(1, :) - (s + 1))) > 1.e-8_dp) ERROR STOP 'Inversion atomic rank wrong'
               IF (MAXVAL(ABS(chi(2, :) - parity*onsite(2, :))) > 1.e-8_dp) ERROR STOP 'Inversion EBR phase wrong'
            END DO
         END DO
      END DO
      seed = [0.173_dp, 0.217_dp, 0.319_dp]
      CALL atomic_band_construct(cell, rotations(:, :, :4), tau(:, :4), anti(:4), seed, .TRUE., &
                                 band, 1.e-9_dp, error, info)
      IF (info /= 0) ERROR STOP 'Generic inversion orbit failed'
      CALL atomic_band_characters(band, [0.5_dp, 0.0_dp, 0.0_dp], mapping, chi, 1.e-9_dp, info)
      IF (info /= 0) ERROR STOP 'Generic inversion characters failed'
      IF (SIZE(chi, 2) /= 1) ERROR STOP 'Generic inversion onsite corep wrong'
      IF (MAXVAL(ABS(chi(:, 1) - [4, 0])) > 1.e-8_dp) ERROR STOP 'Generic inversion induced trace wrong'

      ! A screw has no fixed site: its induced trace vanishes, although its eigenvalues have k phases.
      rotations(:, :, 2) = rotations(:, :, 1)
      rotations(1, 1, 2) = -1
      rotations(2, 2, 2) = -1
      rotations(:, :, 3:4) = rotations(:, :, 1:2)
      tau(3, [2, 4]) = 0.5_dp
      DO s = 0, 1
         CALL atomic_band_construct(cell, rotations(:, :, :4), tau(:, :4), anti(:4), [0.0_dp, 0.0_dp, 0.0_dp], &
                                    s == 1, band, 1.e-9_dp, error, info)
         IF (info /= 0) ERROR STOP 'Screw band construction failed'
         CALL atomic_band_characters(band, [0.0_dp, 0.0_dp, 0.25_dp], mapping, chi, 1.e-9_dp, info)
         IF (info /= 0) ERROR STOP 'Screw band characters failed'
         IF (MAXVAL(ABS(chi(:, 1) - [2*(s + 1), 0])) > 1.e-8_dp) ERROR STOP 'Screw orbit rank/trace wrong'
      END DO

      tau = 0.0_dp
      rotations(:, :, 3) = rotations(:, :, 2)
      rotations(2, 2, 3) = 1
      rotations(3, 3, 3) = -1
      rotations(:, :, 4) = MATMUL(rotations(:, :, 2), rotations(:, :, 3))
      rotations(:, :, 5:8) = rotations(:, :, 1:4)
      anti(:4) = .FALSE.
      anti(5:8) = .TRUE.
      CALL atomic_band_construct(cell, rotations, tau, anti, [0.0_dp, 0.0_dp, 0.0_dp], .TRUE., &
                                 band, 1.e-9_dp, error, info)
      IF (info /= 0) ERROR STOP 'Quaternion site failed'
      CALL atomic_band_characters(band, [0.0_dp, 0.0_dp, 0.0_dp], mapping, chi, 1.e-9_dp, info)
      IF (info /= 0) ERROR STOP 'Quaternion character failed'
      IF (SIZE(chi, 2) /= 1 .OR. MAXVAL(ABS(chi(:, 1) - [2, 0, 0, 0])) > 1.e-8_dp) THEN
         ERROR STOP 'Quaternion onsite irrep incorrectly doubled'
      END IF

      ! C3 with an antiunitary half translation. Some orbit representatives are necessarily antiunitary.
      rotations(:, :, 2) = RESHAPE([0, 1, 0, 0, 0, 1, 1, 0, 0], [3, 3])
      rotations(:, :, 3) = MATMUL(rotations(:, :, 2), rotations(:, :, 2))
      rotations(:, :, 4:6) = rotations(:, :, 1:3)
      anti(:3) = .FALSE.
      anti(4:6) = .TRUE.
      tau = 0.0_dp
      tau(:, 4:6) = 0.5_dp
      DO s = 0, 1
         CALL atomic_band_construct(cell, rotations(:, :, :6), tau(:, :6), anti(:6), [0.0_dp, 0.0_dp, 0.0_dp], &
                                    s == 1, band, 1.e-9_dp, error, info)
         IF (info /= 0) ERROR STOP 'Antiunitary orbit construction failed'
         CALL atomic_band_characters(band, [0.0_dp, 0.0_dp, 0.0_dp], mapping, chi, 1.e-9_dp, info)
         IF (info /= 0) ERROR STOP 'Antiunitary orbit characters failed'
         IF (ANY(ABS(chi(1, :) - 2) > 1.e-8_dp)) ERROR STOP 'Magnetic orbit rank wrong'
         IF (MAXVAL(ABS(AIMAG(chi))) > 1.e-8_dp) ERROR STOP 'Antiunitary coset conjugation missing'
         IF (s == 0) THEN
            IF (COUNT(ABS(chi(2, :) + 1) < 1.e-8_dp) /= 2) ERROR STOP 'C3 conjugate orbital pair wrong'
            IF (COUNT(ABS(chi(2, :) - 2) < 1.e-8_dp) /= 1) ERROR STOP 'C3 real orbital wrong'
         END IF
         DO q = 0, 2
            k = 0.25_dp*q
            CALL point_characters(cell, rotations(:, :, :6), tau(:, :6), anti(:6), k, s == 1, band, point, counts)
         END DO
         ! Deliberately interleave antiunitary operations and move the identity from input index one.
         perm = [5, 2, 4, 1, 6, 3]
         CALL atomic_band_construct(cell, rotations(:, :, perm), tau(:, perm), anti(perm), [0.0_dp, 0.0_dp, 0.0_dp], &
                                    s == 1, other, 1.e-9_dp, error, info)
         IF (info /= 0) ERROR STOP 'Permuted magnetic orbit failed'
         CALL atomic_band_characters(other, [0.0_dp, 0.0_dp, 0.0_dp], alternate, reference, 1.e-9_dp, info)
         IF (info /= 0) ERROR STOP 'Permuted magnetic character failed'
         DO i = 1, SIZE(chi, 2)
            a = 0
            DO j = 1, SIZE(reference, 2)
               error = 0.0_dp
               DO q = 1, SIZE(mapping)
                  n = FINDLOC(perm(alternate), mapping(q), DIM=1)
                  IF (n == 0) ERROR STOP 'Permuted original operation mapping wrong'
                  error = MAX(error, ABS(chi(q, i) - reference(n, j)))
               END DO
               IF (error < 1.e-8_dp) a = a + 1
            END DO
            IF (a == 0) ERROR STOP 'Induction depends on coset representative choice'
         END DO
      END DO

      seed = 0.0_dp
      seed(1) = ieee_value(0.0_dp, ieee_quiet_nan)
      CALL atomic_band_construct(cell, rotations(:, :, :6), tau(:, :6), anti(:6), seed, .FALSE., &
                                 band, 1.e-9_dp, error, info)
      IF (info == 0) ERROR STOP 'Nonfinite site accepted'
      CALL atomic_band_site_data(band, positions, site_map, onsite)
      IF (ALLOCATED(positions) .OR. ALLOCATED(onsite)) ERROR STOP 'Failed construction retained data'
      CALL atomic_band_characters(band, [0.0_dp, 0.0_dp, 0.0_dp], mapping, chi, 1.e-9_dp, info)
      IF (info == 0 .OR. ALLOCATED(chi)) ERROR STOP 'Empty atomic band accepted'
      CALL atomic_band_characters(other, seed, mapping, chi, 1.e-9_dp, info)
      IF (info == 0 .OR. ALLOCATED(chi)) ERROR STOP 'Nonfinite reciprocal point accepted'
      CALL atomic_band_construct(cell, rotations(:, :, :2), tau(:, :2), anti(:2), [0.0_dp, 0.0_dp, 0.0_dp], &
                                 .FALSE., band, 1.e-9_dp, error, info)
      IF (info == 0) ERROR STOP 'Nonclosed site group accepted'
      WRITE (io_unit, *) 'Analytic atomic bands, antiunitary cosets and invalid-input checks passed.'
   END SUBROUTINE analytic_tests

! **************************************************************************************************
!> \brief Independently decompose every induced band into native little-group irreps/corepresentations.
!> \param cell Cartesian cell vectors
!> \param rotations Fractional direct-space operations
!> \param tau Translations of the supplied representatives
!> \param anti Antiunitary flags
!> \param k Fractional reciprocal vector
!> \param spinful Use spin-1/2 factors
!> \param band Prepared atomic reference
!> \param point Independent little-group data for subsequent compatibility checks
!> \param counts Unit-irrep multiplicities for every induced local-representation column
! **************************************************************************************************
   SUBROUTINE point_characters(cell, rotations, tau, anti, k, spinful, band, point, counts)
      REAL(dp), INTENT(IN)                               :: cell(3, 3)
      INTEGER, INTENT(IN)                                :: rotations(:, :, :)
      REAL(dp), INTENT(IN)                               :: tau(:, :)
      LOGICAL, INTENT(IN)                                :: anti(:)
      REAL(dp), INTENT(IN)                               :: k(3)
      LOGICAL, INTENT(IN)                                :: spinful
      TYPE(atomic_band_type), INTENT(IN)                 :: band
      TYPE(little_group_data_type), INTENT(OUT)          :: point
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: counts(:, :)

      COMPLEX(dp), ALLOCATABLE                           :: chi(:, :)
      INTEGER                                            :: g, i, info, n, nu
      INTEGER, ALLOCATABLE                               :: corep_counts(:), INDEX(:), &
                                                            induced_map(:), mapping(:), &
                                                            PRODUCT(:, :), unitary(:)
      REAL(dp)                                           :: error
      TYPE(corepresentation_table_type)                  :: coreps

      CALL little_group_factor_system(cell, rotations, tau, anti, k, spinful, &
                                      mapping, point%product, point%factors, 1.e-9_dp, info)
      IF (info /= 0) ERROR STOP 'Atomic-band test little group failed'
      n = SIZE(mapping)
      point%antiunitary = anti(mapping)
      point%cell = cell
      point%k = k
      point%spinful = spinful
      point%rotations = rotations(:, :, mapping)
      point%translations = tau(:, mapping)
      unitary = PACK([(g, g=1, n)],.NOT. point%antiunitary)
      nu = SIZE(unitary)
      ALLOCATE (INDEX(n), PRODUCT(nu, nu))
      index = 0
      DO g = 1, nu
         INDEX(unitary(g)) = g
      END DO
      DO g = 1, nu
         PRODUCT(:, g) = INDEX(point%product(unitary, unitary(g)))
      END DO
      CALL projective_characters(product, point%factors(unitary, unitary), point%characters, 1.e-8_dp, error, info)
      IF (info /= 0) ERROR STOP 'Atomic-band test little-group characters failed'
      CALL atomic_band_characters(band, k, induced_map, chi, 1.e-9_dp, info)
      IF (info /= 0) ERROR STOP 'Induced atomic characters failed'
      IF (SIZE(induced_map) /= nu) ERROR STOP 'Induced unitary group size wrong'
      ! Each call moves the identity to the first input position before selecting its little group.
      DO g = 1, nu
         INDEX(g) = FINDLOC(induced_map, mapping(unitary(g)), DIM=1)
         IF (INDEX(g) == 0) ERROR STOP 'Induced unitary member missing'
      END DO
      ALLOCATE (counts(SIZE(point%characters, 2), SIZE(chi, 2)))
      IF (ANY(point%antiunitary)) THEN
         CALL projective_corepresentations(point%product, point%factors, point%antiunitary, point%characters, &
                                           coreps, 1.e-8_dp, error, info)
         IF (info /= 0) ERROR STOP 'Induced band target corepresentations failed'
         ALLOCATE (corep_counts(SIZE(coreps%dimension)))
      END IF
      DO i = 1, SIZE(chi, 2)
         CALL character_multiplicities(chi(INDEX(:nu), i), point%characters, counts(:, i), 1.e-7_dp, info)
         IF (info /= 0) ERROR STOP 'Induced atomic band has invalid unitary multiplicities'
         IF (SUM(counts(:, i)*NINT(REAL(point%characters(1, :), dp))) /= NINT(REAL(chi(1, i), dp))) THEN
            ERROR STOP 'Induced atomic band rank lost'
         END IF
         IF (ANY(point%antiunitary)) THEN
            CALL corepresentation_multiplicities(coreps, counts(:, i), corep_counts, info)
            IF (info /= 0) ERROR STOP 'Induced atomic band violates antiunitary pairing'
         END IF
      END DO
      point%multiplicity = counts(:, 1)
   END SUBROUTINE point_characters

END PROGRAM topology_band_unittest
