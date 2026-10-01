!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
PROGRAM topology_reciprocal_unittest
   USE ieee_arithmetic,                 ONLY: ieee_quiet_nan,&
                                              ieee_value
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE mathlib,                         ONLY: inv_3x3
   USE topology_atomic_signatures,      ONLY: atomic_signature_construct,&
                                              atomic_signature_star_map,&
                                              atomic_signature_type
   USE topology_compatibility,          ONLY: little_group_compatibility,&
                                              little_group_edge_type
   USE topology_compatibility_lattice,  ONLY: compatibility_lattice_construct,&
                                              compatibility_lattice_type
   USE topology_integer_lattice,        ONLY: integer_lattice_certificate,&
                                              integer_lattice_factor,&
                                              integer_lattice_type,&
                                              integer_quotient_class
   USE topology_little_group,           ONLY: little_group_factor_system
   USE topology_reciprocal,             ONLY: reciprocal_refine,&
                                              reciprocal_sampling,&
                                              reciprocal_sampling_type,&
                                              reciprocal_strata,&
                                              reciprocal_stratum_type
   USE topology_star_transport,         ONLY: star_operation

   IMPLICIT NONE
   CALL analytic_tests()
   CALL refinement_tests()
   CALL star_analytic_tests()
   PRINT *, 'Reciprocal-stratum tests passed.'
CONTAINS

! **************************************************************************************************
!> \brief Complex screw characters, antiunitary conjugation and independent endpoint coset gauges.
! **************************************************************************************************
   SUBROUTINE star_analytic_tests()
      COMPLEX(dp)                                        :: saved
      INTEGER                                            :: a, b, c, d, g, i, ident(3, 3), info, &
                                                            operation, r(3, 3, 6), shift(3), spin
      INTEGER(int_8), ALLOCATABLE                        :: bad(:), indices(:)
      INTEGER, ALLOCATABLE                               :: permutation(:, :)
      LOGICAL                                            :: anti(6), compatible
      REAL(dp)                                           :: cell(3, 3), changed_tau(3, 6), error, &
                                                            k(3, 3), tau(3, 6)
      TYPE(atomic_signature_type)                        :: catalogue, shifted
      TYPE(compatibility_lattice_type)                   :: lattice

      ident = 0
      DO i = 1, 3
         ident(i, i) = 1
      END DO
      cell = REAL(ident, dp)
      cell(1, 2) = -0.5_dp
      cell(2, 2) = SQRT(3.0_dp)/2
      r(:, :, 1) = ident
      r(:, :, 2) = RESHAPE([0, 1, 0, -1, -1, 0, 0, 0, 1], [3, 3])
      r(:, :, 3) = MATMUL(r(:, :, 2), r(:, :, 2))
      r(:, :, 4:6) = r(:, :, 1:3)
      tau = 0.0_dp
      tau(3, 2:3) = [1.0_dp/3, 2.0_dp/3]
      tau(:, 4:6) = tau(:, 1:3)
      anti(:3) = .FALSE.
      anti(4:) = .TRUE.
      k = 0.0_dp
      k(3, :) = [0.25_dp, 0.75_dp, 1.25_dp]
      changed_tau = tau
      DO g = 2, 6
         changed_tau(:, g) = changed_tau(:, g) + REAL([MOD(g, 2), -MOD(g, 3), MOD(g + 1, 3)], dp)
      END DO
      DO spin = 0, 1
         CALL atomic_signature_construct(cell, r, tau, anti, k, spin == 1, catalogue, 1.e-8_dp, info)
         IF (info /= 0) ERROR STOP 'Screw star catalogue failed'
         CALL atomic_signature_construct(cell, r, changed_tau, anti, k, spin == 1, shifted, 1.e-8_dp, info)
         IF (info /= 0) ERROR STOP 'Coset-shifted screw star catalogue failed'
         CALL compatibility_lattice_construct(catalogue, [INTEGER ::], [INTEGER ::], &
                                              RESHAPE([INTEGER ::], [3, 0]), lattice, 1.e-8_dp, info, 64.0_dp)
         IF (info /= 0) ERROR STOP 'Star-only compatibility lattice failed'
         IF (SIZE(lattice%stars, 2) /= 2 .OR. SIZE(lattice%quotient%invariants) /= 3) THEN
            ERROR STOP 'Star-only independent character count wrong'
         END IF
         bad = INT(catalogue%matrix(:, 1), int_8)
         bad(catalogue%row_begin(2)) = bad(catalogue%row_begin(2)) + 1
         CALL integer_quotient_class(lattice%quotient, bad, compatible, indices, info)
         IF (info /= 0 .OR. compatible) ERROR STOP 'Star-inconsistent band signature accepted'
         CALL star_operation(catalogue%group, k(:, 1), k(:, 2), operation, shift, 1.e-8_dp, info)
         IF (info /= 0 .OR. .NOT. anti(operation)) ERROR STOP 'Missing antiunitary star operation'
         IF (ANY(shift /= [0, 0, -1])) ERROR STOP 'Star reciprocal shift lost'
         CALL atomic_signature_star_map(catalogue, 1, 2, operation, permutation, 1.e-8_dp, error, info)
         IF (info /= 0) ERROR STOP 'Complex screw/time-reversal transport failed'
         a = catalogue%row_begin(1)
         b = catalogue%row_begin(2) - 1
         c = catalogue%row_begin(2)
         d = catalogue%row_begin(3) - 1
         IF (ANY(MATMUL(permutation, catalogue%matrix(a:b, :)) /= catalogue%matrix(c:d, :))) THEN
            ERROR STOP 'Screw antiunitary transport disagrees with induced characters'
         END IF
         ! Change only the destination representative gauge, keeping the cached full group fixed.
         catalogue%points(2) = shifted%points(2)
         CALL atomic_signature_star_map(catalogue, 1, 2, operation, permutation, 1.e-8_dp, error, info)
         IF (info /= 0) ERROR STOP 'Independent endpoint coset gauge failed'
         IF (ANY(MATMUL(permutation, catalogue%matrix(a:b, :)) /= shifted%matrix(c:d, :))) THEN
            ERROR STOP 'Endpoint coset phase missing in star transport'
         END IF
         saved = catalogue%points(2)%little%factors(2, 2)
         catalogue%points(2)%little%factors(2, 2) = -saved
         CALL atomic_signature_star_map(catalogue, 1, 2, operation, permutation, 1.e-8_dp, error, info)
         IF (info == 0 .OR. ALLOCATED(permutation)) ERROR STOP 'Corrupt star factor accepted'
         catalogue%points(2)%little%factors(2, 2) = saved
         CALL star_operation(catalogue%group, k(:, 1), k(:, 3), operation, shift, 1.e-8_dp, info)
         IF (info /= 0 .OR. operation /= 1 .OR. ANY(shift /= [0, 0, -1])) ERROR STOP 'Reciprocal image not identified'
         CALL atomic_signature_star_map(catalogue, 1, 3, operation, permutation, 1.e-8_dp, error, info)
         IF (info /= 0) ERROR STOP 'Integer reciprocal image transport failed'
         CALL atomic_signature_star_map(catalogue, 1, 3, 0, permutation, 1.e-8_dp, error, info)
         IF (info == 0 .OR. ALLOCATED(permutation)) ERROR STOP 'Invalid star operation accepted'
         CALL star_operation(catalogue%group, k(:, 1), [0.123_dp, 0.222_dp, 0.111_dp], &
                             operation, shift, 1.e-8_dp, info)
         IF (info /= 1 .OR. operation /= 0) ERROR STOP 'Different reciprocal stars identified'
         CALL star_operation(catalogue%group, [HUGE(0.0_dp), 0.0_dp, 0.0_dp], k(:, 1), &
                             operation, shift, 1.e-8_dp, info)
         IF (info == 0 .OR. operation /= 0) ERROR STOP 'Unsafe star coordinate accepted'
         catalogue%points(3)%little%k(3) = HUGE(0.0_dp)
         CALL atomic_signature_star_map(catalogue, 1, 3, 1, permutation, 1.e-8_dp, error, info)
         IF (info == 0 .OR. ALLOCATED(permutation)) ERROR STOP 'Unsafe transport phase accepted'
         catalogue%points(3)%little%k = k(:, 3)
         catalogue%points(3)%little%characters(1, 1) = CMPLX(ieee_value(0.0_dp, ieee_quiet_nan), 0.0_dp, dp)
         CALL atomic_signature_star_map(catalogue, 1, 3, 1, permutation, 1.e-8_dp, error, info)
         IF (info == 0 .OR. ALLOCATED(permutation)) ERROR STOP 'Nonfinite star character accepted'
      END DO
      PRINT *, 'Star transport analytic, antiunitary and endpoint-gauge tests passed.'
   END SUBROUTINE star_analytic_tests

! **************************************************************************************************
!> \brief Preserve lifted paths, including negative and multiple reciprocal-boundary crossings.
! **************************************************************************************************
   SUBROUTINE refinement_tests()
      INTEGER                                            :: info, n
      TYPE(reciprocal_sampling_type)                     :: original, refined

      ALLOCATE (original%kpoints(3, 2), original%source(5), original%destination(5), original%shift(3, 5))
      original%kpoints(:, 1) = [0.13_dp, 0.29_dp, 0.77_dp]
      original%kpoints(:, 2) = [0.81_dp, 0.11_dp, 0.17_dp]
      original%source = [1, 2, 1, 2, 1]
      original%destination = [2, 1, 1, 2, 1]
      original%shift(:, 1) = [2, -1, 3]
      original%shift(:, 2) = [-2, 1, -3]
      original%shift(:, 3) = [1, 0, -1]
      original%shift(:, 4) = [-1, -2, 0]
      original%shift(:, 5) = 0
      DO n = 1, 7
         CALL check_refinement(original, n)
      END DO
      refined = original
      CALL reciprocal_refine(refined, 4, 3, 1.e-9_dp, info)
      IF (info /= -3) ERROR STOP 'Refinement point budget was not enforced'
      IF (SIZE(refined%kpoints, 2) /= 2 .OR. SIZE(refined%source) /= 5) ERROR STOP 'Budget failure mutated sampling'
      IF (ANY(refined%kpoints /= original%kpoints) .OR. ANY(refined%shift /= original%shift)) THEN
         ERROR STOP 'Budget failure changed path data'
      END IF
      CALL reciprocal_refine(refined, 0, 100, 1.e-9_dp, info)
      IF (info /= -1) ERROR STOP 'Zero subdivision count accepted'
      refined%source(1) = 0
      CALL reciprocal_refine(refined, 2, 100, 1.e-9_dp, info)
      IF (info /= -1) ERROR STOP 'Invalid source index accepted'
      refined = original
      refined%kpoints(1, 1) = ieee_value(0.0_dp, ieee_quiet_nan)
      CALL reciprocal_refine(refined, 2, 100, 1.e-9_dp, info)
      IF (info /= -1) ERROR STOP 'Nonfinite path accepted'
      refined = original
      refined%kpoints(1, 1) = HUGE(0.0_dp)
      CALL reciprocal_refine(refined, 2, 100, 1.e-9_dp, info)
      IF (info /= -1) ERROR STOP 'Unwrapped path accepted'
      PRINT *, 'Reciprocal path refinement tests passed.'
   END SUBROUTINE refinement_tests

! **************************************************************************************************
!> \brief Independently compare lifted endpoints, accumulated winding and uniform increments.
!> \param original Original sampling
!> \param intervals Subdivision count
! **************************************************************************************************
   SUBROUTINE check_refinement(original, intervals)
      TYPE(reciprocal_sampling_type), INTENT(IN)         :: original
      INTEGER, INTENT(IN)                                :: intervals

      INTEGER                                            :: e, info, j, k, p, q, winding(3)
      REAL(dp)                                           :: direction(3), expected(3), lifted(3), &
                                                            step(3)
      TYPE(reciprocal_sampling_type)                     :: refined

      refined = original
      CALL reciprocal_refine(refined, intervals, SIZE(original%kpoints, 2) + &
                             SIZE(original%source)*(intervals - 1), 1.e-9_dp, info)
      IF (info /= 0) ERROR STOP 'Path refinement failed'
      IF (SIZE(refined%source) /= intervals*SIZE(original%source)) ERROR STOP 'Refined edge count wrong'
      DO e = 1, SIZE(original%source)
         direction = original%kpoints(:, original%destination(e)) + REAL(original%shift(:, e), dp) - &
                     original%kpoints(:, original%source(e))
         winding = 0
         DO j = 1, intervals
            k = (e - 1)*intervals + j
            p = refined%source(k)
            q = refined%destination(k)
            IF (j == 1) THEN
               IF (p /= original%source(e)) ERROR STOP 'Refinement changed source'
            ELSE
               IF (p /= refined%destination(k - 1)) ERROR STOP 'Refined path is discontinuous'
            END IF
            winding = winding + refined%shift(:, k)
            step = refined%kpoints(:, q) + REAL(refined%shift(:, k), dp) - refined%kpoints(:, p)
            IF (MAXVAL(ABS(step - direction/intervals)) > 1.e-8_dp) ERROR STOP 'Refined step wrong'
            lifted = refined%kpoints(:, q) + REAL(winding, dp)
            expected = original%kpoints(:, original%source(e)) + REAL(j, dp)/intervals*direction
            IF (MAXVAL(ABS(lifted - expected)) > 1.e-8_dp) ERROR STOP 'Refined lifted point wrong'
         END DO
         IF (q /= original%destination(e)) ERROR STOP 'Refinement changed destination'
         IF (ANY(winding /= original%shift(:, e))) ERROR STOP 'Refinement changed total winding'
      END DO
   END SUBROUTINE check_refinement

! **************************************************************************************************
!> \brief Analytic symmorphic, nonsymmorphic and magnetic cases, including malformed input.
! **************************************************************************************************
   SUBROUTINE analytic_tests()
      INTEGER                                            :: g, hist(0:3), i, ident(3, 3), info, &
                                                            r(3, 3, 8)
      INTEGER(int_8)                                     :: comparisons
      LOGICAL                                            :: anti(8)
      REAL(dp)                                           :: cell(3, 3), tau(3, 8)
      TYPE(reciprocal_stratum_type), ALLOCATABLE         :: strata(:)

      ident = 0
      DO i = 1, 3
         ident(i, i) = 1
      END DO
      cell = REAL(ident, dp)
      r(:, :, 1) = ident
      tau = 0.0_dp
      anti = .FALSE.
      CALL reciprocal_strata(cell, r(:, :, :1), tau(:, :1), anti(:1), strata, 1.e-9_dp, info)
      IF (info /= 0) ERROR STOP 'P1 reciprocal strata failed'
      IF (SIZE(strata) /= 1 .OR. strata(1)%dimension /= 3) ERROR STOP 'P1 reciprocal stratum wrong'
      CALL check_geometry(cell, r(:, :, :1), tau(:, :1), anti(:1), strata)
      r(:, :, 2) = ident
      anti(2) = .TRUE.
      CALL reciprocal_strata(cell, r(:, :, :2), tau(:, :2), anti(:2), strata, 1.e-9_dp, info)
      IF (info /= 0) ERROR STOP 'Time-reversal reciprocal strata failed'
      IF (SIZE(strata) /= 9 .OR. COUNT(strata%dimension == 0) /= 8) ERROR STOP 'Eight TRIM were not retained'
      CALL check_geometry(cell, r(:, :, :2), tau(:, :2), anti(:2), strata)
      ! PT acts trivially on k but its full antiunitary stabilizer must not disappear.
      r(:, :, 2) = -ident
      CALL reciprocal_strata(cell, r(:, :, :2), tau(:, :2), anti(:2), strata, 1.e-9_dp, info)
      IF (info /= 0) ERROR STOP 'PT reciprocal strata failed'
      IF (SIZE(strata) /= 1) ERROR STOP 'PT generated false special k points'
      IF (SIZE(strata(1)%stabilizer) /= 2) ERROR STOP 'PT lost its antiunitary stabilizer'
      CALL check_geometry(cell, r(:, :, :2), tau(:, :2), anti(:2), strata)
      r(:, :, 2) = ident
      tau(:, 2) = [0.0_dp, 0.0_dp, 0.5_dp]
      CALL reciprocal_strata(cell, r(:, :, :2), tau(:, :2), anti(:2), strata, 1.e-9_dp, info)
      IF (info /= 0 .OR. COUNT(strata%dimension == 0) /= 8) ERROR STOP 'Magnetic half translation lost TRIM'
      CALL check_geometry(cell, r(:, :, :2), tau(:, :2), anti(:2), strata)
      anti = .FALSE.
      r(:, :, 2) = -ident
      r(3, 3, 2) = 1
      CALL reciprocal_strata(cell, r(:, :, :2), tau(:, :2), anti(:2), strata, 1.e-9_dp, info)
      IF (info /= 0) ERROR STOP 'Screw reciprocal strata failed'
      IF (SIZE(strata) /= 5 .OR. COUNT(strata%dimension == 1) /= 4) ERROR STOP 'Screw invariant lines lost'
      CALL check_geometry(cell, r(:, :, :2), tau(:, :2), anti(:2), strata)
      CALL sampling_test(cell, r(:, :, :2), tau(:, :2), anti(:2), comparisons)
      r(:, :, 2) = ident
      CALL reciprocal_strata(cell, r(:, :, :2), tau(:, :2), anti(:2), strata, 1.e-9_dp, info)
      IF (info /= 0) ERROR STOP 'Centered reciprocal strata failed'
      IF (SIZE(strata) /= 1 .OR. SIZE(strata(1)%stabilizer) /= 2) ERROR STOP 'Centering cosets lost'
      CALL check_geometry(cell, r(:, :, :2), tau(:, :2), anti(:2), strata)
      tau = 0.0_dp
      DO g = 1, 8
         r(:, :, g) = 0
         DO i = 1, 3
            r(i, i, g) = MERGE(-1, 1, BTEST(g - 1, i - 1))
         END DO
      END DO
      CALL reciprocal_strata(cell, r, tau, anti, strata, 1.e-9_dp, info)
      IF (info /= 0) ERROR STOP 'Orthorhombic reciprocal strata failed'
      DO i = 0, 3
         hist(i) = COUNT(strata%dimension == i)
      END DO
      IF (ANY(hist /= [8, 12, 6, 1])) ERROR STOP 'Orthorhombic reciprocal intersections incomplete'
      CALL check_geometry(cell, r, tau, anti, strata)
      CALL transformed_test(cell, r, tau, anti, strata)
      CALL sampling_test(cell, r, tau, anti, comparisons)
      tau(1, 2) = 0.317_dp
      CALL reciprocal_strata(cell, r, tau, anti, strata, 1.e-9_dp, info)
      IF (info == 0 .OR. ALLOCATED(strata)) ERROR STOP 'Invalid affine group accepted'
      CALL reciprocal_strata(cell, r(:, :, :0), tau(:, :0), anti(:0), strata, 1.e-9_dp, info)
      IF (info == 0 .OR. ALLOCATED(strata)) ERROR STOP 'Empty reciprocal group accepted'
      PRINT *, 'Reciprocal analytic, magnetic and invalid-data tests passed.'
   END SUBROUTINE analytic_tests

! **************************************************************************************************
!> \brief Check all original stabilizers, exact primitive cycles, and lifted incidence segments.
!> \param cell Direct cell
!> \param rotations Original rotations
!> \param translations Original translations
!> \param anti Original antiunitary flags
!> \param strata Generated reciprocal families
! **************************************************************************************************
   SUBROUTINE check_geometry(cell, rotations, translations, anti, strata)
      REAL(dp), INTENT(IN)                               :: cell(3, 3)
      INTEGER, INTENT(IN)                                :: rotations(:, :, :)
      REAL(dp), INTENT(IN)                               :: translations(:, :)
      LOGICAL, INTENT(IN)                                :: anti(:)
      TYPE(reciprocal_stratum_type), INTENT(IN)          :: strata(:)

      COMPLEX(dp), ALLOCATABLE                           :: factor(:, :)
      INTEGER                                            :: d, g, h, i, info, j, n, rank
      INTEGER(int_8), ALLOCATABLE                        :: diagonal(:), left(:, :), right(:, :)
      INTEGER, ALLOCATABLE                               :: mapping(:), products(:, :)
      REAL(dp)                                           :: delta(3), direction(3), image(3), sign
      TYPE(integer_lattice_type)                         :: integer_cycles

      n = SIZE(rotations, 3)
      DO i = 1, SIZE(strata)
         CALL little_group_factor_system(cell, rotations, translations, anti, strata(i)%k, .TRUE., &
                                         mapping, products, factor, 1.e-8_dp, info)
         IF (info /= 0) ERROR STOP 'Reciprocal seed has invalid full spinful little group'
         IF (SIZE(mapping) /= SIZE(strata(i)%stabilizer)) ERROR STOP 'Reciprocal stabilizer order mismatch'
         IF (SIZE(mapping)*strata(i)%star_size /= n) ERROR STOP 'Reciprocal orbit stabilizer mismatch'
         IF (ANY([(COUNT(strata(i)%stabilizer == mapping(g)), g=1, SIZE(mapping))] /= 1)) THEN
            ERROR STOP 'Reciprocal stabilizer membership mismatch'
         END IF
         d = strata(i)%dimension
         IF (ANY(SHAPE(strata(i)%cycles) /= [3, d])) ERROR STOP 'Reciprocal cycle dimension wrong'
         IF (d > 0) THEN
            CALL integer_lattice_factor(INT(strata(i)%cycles, int_8), integer_cycles, info)
            IF (info /= 0) ERROR STOP 'Reciprocal cycle factorization failed'
            CALL integer_lattice_certificate(integer_cycles, left, right, diagonal, rank)
            IF (rank /= d .OR. ANY(diagonal /= 1)) ERROR STOP 'Reciprocal cycles not primitive'
         END IF
         IF (ANY(MATMUL(INT(strata(i)%normal, int_8), INT(strata(i)%cycles, int_8)) /= 0)) THEN
            ERROR STOP 'Reciprocal cycle leaves closure'
         END IF
         DO j = 1, d
            direction = REAL(strata(i)%cycles(:, j), dp)
            DO h = 1, SIZE(mapping)
               g = mapping(h)
               delta = MATMUL(REAL(TRANSPOSE(rotations(:, :, g)), dp), direction) &
                       - MERGE(-1.0_dp, 1.0_dp, anti(g))*direction
               IF (MAXVAL(ABS(delta)) > 1.e-8_dp) ERROR STOP 'Cycle not fixed pointwise by stabilizer'
            END DO
         END DO
         DO j = 1, SIZE(strata(i)%specializations)
            h = strata(i)%specializations(j)
            g = strata(i)%image_operation(j)
            sign = MERGE(-1.0_dp, 1.0_dp, anti(g))
            image = sign*MATMUL(TRANSPOSE(inv_3x3(REAL(rotations(:, :, g), dp))), strata(h)%k) &
                    + REAL(strata(i)%image_shift(:, j), dp)
            direction = image - strata(i)%k
            IF (strata(h)%dimension >= strata(i)%dimension) ERROR STOP 'Incidence dimension does not decrease'
            delta = MATMUL(REAL(strata(i)%normal, dp), direction)
            IF (MAXVAL(ABS(delta)) > 1.e-8_dp) ERROR STOP 'Incidence left the lifted parent closure'
            DO h = 1, SIZE(mapping)
               g = mapping(h)
               delta = MATMUL(REAL(TRANSPOSE(rotations(:, :, g)), dp), direction) &
                       - MERGE(-1.0_dp, 1.0_dp, anti(g))*direction
               IF (MAXVAL(ABS(delta)) > 1.e-8_dp) ERROR STOP 'Incidence lost parent stabilizer'
            END DO
         END DO
      END DO
   END SUBROUTINE check_geometry

! **************************************************************************************************
!> \brief Independent rational-grid orbit/stabilizer coverage, using the direct transpose test.
!> \param rotations Original direct rotations
!> \param anti Original antiunitary flags
!> \param strata Generated reciprocal families
!> \param denominator Uniform reciprocal mesh size
! **************************************************************************************************
   SUBROUTINE grid_test(rotations, anti, strata, denominator)
      INTEGER, INTENT(IN)                                :: rotations(:, :, :)
      LOGICAL, INTENT(IN)                                :: anti(:)
      TYPE(reciprocal_stratum_type), INTENT(IN)          :: strata(:)
      INTEGER, INTENT(IN)                                :: denominator

      INTEGER                                            :: count_fixed, found, g, i, n, x, y, z
      LOGICAL                                            :: matched
      REAL(dp)                                           :: delta(3), image(3), k(3), sign

      n = SIZE(rotations, 3)
      DO z = 0, denominator - 1
         DO y = 0, denominator - 1
            DO x = 0, denominator - 1
               k = REAL([x, y, z], dp)/denominator
               count_fixed = 0
               DO g = 1, n
                  sign = MERGE(-1.0_dp, 1.0_dp, anti(g))
                  delta = MATMUL(REAL(TRANSPOSE(rotations(:, :, g)), dp), k) - sign*k
                  IF (MAXVAL(ABS(delta - ANINT(delta))) < 1.e-8_dp) count_fixed = count_fixed + 1
               END DO
               found = 0
               DO i = 1, SIZE(strata)
                  IF (SIZE(strata(i)%stabilizer) /= count_fixed) CYCLE
                  matched = .FALSE.
                  DO g = 1, n
                     ! Inverses run over the same group: W^T suffices for this orbit test.
                     sign = MERGE(-1.0_dp, 1.0_dp, anti(g))
                     image = sign*MATMUL(REAL(TRANSPOSE(rotations(:, :, g)), dp), k)
                     delta = MATMUL(REAL(strata(i)%normal, dp), image) - strata(i)%offset
                     IF (MAXVAL(ABS(delta - ANINT(delta))) < 1.e-8_dp) matched = .TRUE.
                  END DO
                  IF (matched) found = found + 1
               END DO
               IF (found /= 1) THEN
                  PRINT *, 'Grid point, stabilizer, matching strata:', k, count_fixed, found
                  ERROR STOP 'Reciprocal rational grid has missing or ambiguous stratum'
               END IF
            END DO
         END DO
      END DO
   END SUBROUTINE grid_test

! **************************************************************************************************
!> \brief Simultaneous unimodular shear, origin/coset shifts and operation permutation.
!> \param cell Original direct cell
!> \param rotations Original rotations
!> \param translations Original translations
!> \param anti Original antiunitary flags
!> \param reference Original reciprocal families
! **************************************************************************************************
   SUBROUTINE transformed_test(cell, rotations, translations, anti, reference)
      REAL(dp), INTENT(IN)                               :: cell(3, 3)
      INTEGER, INTENT(IN)                                :: rotations(:, :, :)
      REAL(dp), INTENT(IN)                               :: translations(:, :)
      LOGICAL, INTENT(IN)                                :: anti(:)
      TYPE(reciprocal_stratum_type), INTENT(IN)          :: reference(:)

      INTEGER                                            :: g, h, i, ident(3, 3), info, &
                                                            inverse(3, 3), j, n, shear(3, 3)
      LOGICAL                                            :: found, a(SIZE(anti))
      INTEGER                                            :: r(3, 3, SIZE(anti))
      REAL(dp)                                           :: delta(3), image(3), k(3), origin(3), &
                                                            tau(3, SIZE(anti))
      TYPE(reciprocal_stratum_type), ALLOCATABLE         :: changed(:)
      INTEGER :: counts(SIZE(reference)), correspondence(SIZE(reference))

      n = SIZE(anti)
      shear = RESHAPE([1, 0, 0, 1, 1, 0, 0, 1, 1], [3, 3])
      inverse = RESHAPE([1, 0, 0, -1, 1, 0, 1, -1, 1], [3, 3])
      ident = MATMUL(shear, inverse)
      origin = [0.123_dp, -0.217_dp, 0.319_dp]
      DO g = 1, n
         h = MODULO(g, n) + 1
         r(:, :, h) = MATMUL(MATMUL(inverse, rotations(:, :, g)), shear)
         a(h) = anti(g)
         tau(:, h) = MATMUL(REAL(inverse, dp), translations(:, g) &
                            + MATMUL(REAL(rotations(:, :, g) - ident, dp), origin))
         IF (ANY(rotations(:, :, g) /= ident) .OR. anti(g) &
             .OR. MAXVAL(ABS(translations(:, g))) > 1.e-9_dp) tau(:, h) = tau(:, h) + REAL([g, 1 - g, g], dp)
      END DO
      CALL reciprocal_strata(MATMUL(cell, REAL(shear, dp)), r, tau, a, changed, 1.e-8_dp, info)
      IF (info /= 0) ERROR STOP 'Transformed reciprocal strata failed'
      IF (SIZE(changed) /= SIZE(reference)) ERROR STOP 'Coordinate change lost reciprocal strata'
      CALL check_geometry(MATMUL(cell, REAL(shear, dp)), r, tau, a, changed)
      counts = 0
      correspondence = 0
      DO j = 1, SIZE(changed)
         k = MATMUL(REAL(TRANSPOSE(inverse), dp), changed(j)%k)
         DO i = 1, SIZE(reference)
            IF (reference(i)%dimension /= changed(j)%dimension) CYCLE
            IF (reference(i)%star_size /= changed(j)%star_size) CYCLE
            found = .FALSE.
            DO g = 1, n
               image = MERGE(-1.0_dp, 1.0_dp, anti(g))*MATMUL(REAL(TRANSPOSE(rotations(:, :, g)), dp), k)
               delta = MATMUL(REAL(reference(i)%normal, dp), image) - reference(i)%offset
               IF (MAXVAL(ABS(delta - ANINT(delta))) < 1.e-8_dp) found = .TRUE.
            END DO
            IF (.NOT. found) CYCLE
            IF (correspondence(j) /= 0) ERROR STOP 'Coordinate change is not a stratum injection'
            correspondence(j) = i
            counts(i) = counts(i) + 1
         END DO
      END DO
      IF (ANY(counts /= 1)) ERROR STOP 'Coordinate change is not a stratum bijection'
      DO j = 1, SIZE(changed)
         i = correspondence(j)
         IF (SIZE(changed(j)%specializations) /= SIZE(reference(i)%specializations)) THEN
            ERROR STOP 'Coordinate change lost incidence edges'
         END IF
         DO h = 1, SIZE(changed(j)%specializations)
            IF (.NOT. ANY(reference(i)%specializations == correspondence(changed(j)%specializations(h)))) THEN
               ERROR STOP 'Coordinate change altered incidence graph'
            END IF
         END DO
      END DO
   END SUBROUTINE transformed_test

! **************************************************************************************************
!> \brief Every induced atomic column obeys every generated incidence and primitive cycle relation.
!> \param cell Direct cell
!> \param rotations Original group rotations
!> \param translations Original group translations
!> \param anti Original antiunitary flags
!> \param comparisons Number of scalar/spinful atomic-column segment checks
! **************************************************************************************************
   SUBROUTINE sampling_test(cell, rotations, translations, anti, comparisons)
      REAL(dp), INTENT(IN)                               :: cell(3, 3)
      INTEGER, INTENT(IN)                                :: rotations(:, :, :)
      REAL(dp), INTENT(IN)                               :: translations(:, :)
      LOGICAL, INTENT(IN)                                :: anti(:)
      INTEGER(int_8), INTENT(OUT)                        :: comparisons

      INTEGER                                            :: a, b, c, d, e, g, h, info, p, q, spin
      INTEGER(int_8), ALLOCATABLE                        :: indices(:)
      INTEGER, ALLOCATABLE                               :: lc(:, :), left(:, :), rc(:, :), &
                                                            right(:, :)
      LOGICAL                                            :: compatible
      REAL(dp)                                           :: delta(3), direction(3), error
      TYPE(atomic_signature_type)                        :: catalogue
      TYPE(compatibility_lattice_type)                   :: changed, lattice
      TYPE(little_group_edge_type)                       :: edge
      TYPE(reciprocal_sampling_type)                     :: sampling

      comparisons = 0_int_8
      CALL reciprocal_sampling(cell, rotations, translations, anti, sampling, 1.e-9_dp, info)
      IF (info /= 0) ERROR STOP 'Reciprocal sampling graph failed'
      CALL check_refinement(sampling, 3)
      DO e = 1, SIZE(sampling%source)
         p = sampling%source(e)
         q = sampling%destination(e)
         direction = sampling%kpoints(:, q) + REAL(sampling%shift(:, e), dp) - sampling%kpoints(:, p)
         DO h = 1, SIZE(sampling%strata(p)%stabilizer)
            g = sampling%strata(p)%stabilizer(h)
            delta = MATMUL(REAL(TRANSPOSE(rotations(:, :, g)), dp), direction) &
                    - MERGE(-1.0_dp, 1.0_dp, anti(g))*direction
            IF (MAXVAL(ABS(delta)) > 1.e-8_dp) ERROR STOP 'Wrapped sampling lost lifted parent stabilizer'
         END DO
      END DO
      DO spin = 0, 1
         CALL atomic_signature_construct(cell, rotations, translations, anti, sampling%kpoints, spin == 1, &
                                         catalogue, 1.e-8_dp, info)
         IF (info /= 0) THEN
            PRINT *, 'Reciprocal atomic catalogue spin/status:', spin, info
            ERROR STOP 'Reciprocal atomic catalogue failed'
         END IF
         CALL compatibility_lattice_construct(catalogue, sampling%source, sampling%destination, sampling%shift, &
                                              lattice, 1.e-7_dp, info, 2048.0_dp)
         IF (info /= 0) THEN
            PRINT *, 'Sampling quotient spin/status:', spin, info
            ERROR STOP 'Reciprocal compatibility quotient failed'
         END IF
         CALL check_star_maps(catalogue, lattice, rotations, anti)
         PRINT *, 'Graph quotient spin/coordinates/rank/free/cyclic:', spin, SIZE(catalogue%matrix, 1), &
            SIZE(lattice%quotient%invariants), COUNT(lattice%quotient%invariants == 0), &
            COUNT(lattice%quotient%invariants > 1)
         PRINT *, 'Graph star relations:', SIZE(lattice%stars, 2)
         PRINT *, 'Graph cyclic orders:', PACK(lattice%quotient%invariants, lattice%quotient%invariants > 1)
         IF (ANY(MATMUL(lattice%constraints, INT(catalogue%matrix, int_8)) /= 0)) THEN
            ERROR STOP 'Atomic column violates assembled compatibility equations'
         END IF
         DO p = 1, SIZE(catalogue%matrix, 2)
            CALL integer_quotient_class(lattice%quotient, INT(catalogue%matrix(:, p), int_8), &
                                        compatible, indices, info)
            IF (info /= 0 .OR. .NOT. compatible) ERROR STOP 'Atomic column not compatible in joint lattice'
            IF (ANY(indices /= 0)) ERROR STOP 'Atomic column has nonzero joint quotient class'
         END DO
         IF (SIZE(anti) <= 4) THEN
            CALL compatibility_lattice_construct(catalogue, sampling%destination, sampling%source, -sampling%shift, &
                                                 changed, 1.e-7_dp, info, 2048.0_dp)
            IF (info /= 0) ERROR STOP 'Reversed compatibility graph failed'
            IF (SIZE(changed%quotient%invariants) /= SIZE(lattice%quotient%invariants)) THEN
               ERROR STOP 'Edge reversal changed compatibility rank'
            END IF
            IF (ANY(changed%quotient%invariants /= lattice%quotient%invariants)) THEN
               ERROR STOP 'Edge reversal changed quotient factors'
            END IF
            CALL compatibility_lattice_construct(catalogue, sampling%source, sampling%destination, sampling%shift, &
                                                 changed, 1.e-7_dp, info, 1.e-12_dp)
            IF (info /= -3 .OR. ALLOCATED(changed%constraints)) ERROR STOP 'Graph memory budget not enforced'
            q = FINDLOC(ANY(lattice%constraints /= 0, DIM=1), .TRUE., DIM=1)
            IF (q > 0) THEN
               catalogue%matrix(q, 1) = catalogue%matrix(q, 1) + 1
               CALL compatibility_lattice_construct(catalogue, sampling%source, sampling%destination, sampling%shift, &
                                                    changed, 1.e-7_dp, info, 2048.0_dp)
               IF (info /= -2 .OR. ALLOCATED(changed%constraints)) ERROR STOP 'Incompatible reference column accepted'
               catalogue%matrix(q, 1) = catalogue%matrix(q, 1) - 1
            END IF
         END IF
         DO e = 1, SIZE(sampling%source)
            p = sampling%source(e)
            q = sampling%destination(e)
            CALL little_group_compatibility(catalogue%points(p)%little, catalogue%points(q)%little, &
                                            sampling%shift(:, e), edge, 1.e-7_dp, error, info)
            IF (info /= 0) THEN
               PRINT *, 'Sampling spin/edge/status:', spin, e, info
               ERROR STOP 'Automatic reciprocal edge table failed'
            END IF
            a = catalogue%row_begin(p)
            b = catalogue%row_begin(p + 1) - 1
            c = catalogue%row_begin(q)
            d = catalogue%row_begin(q + 1) - 1
            ! Explicit allocation also avoids compiler-dependent automatic MATMUL resizing.
            ALLOCATE (left(SIZE(catalogue%points(p)%little%characters, 2), SIZE(catalogue%matrix, 2)), &
                      right(SIZE(catalogue%points(q)%little%characters, 2), SIZE(catalogue%matrix, 2)))
            IF (ALLOCATED(catalogue%points(p)%coreps%restriction)) THEN
               left(:, :) = MATMUL(catalogue%points(p)%coreps%restriction, catalogue%matrix(a:b, :))
            ELSE
               left(:, :) = catalogue%matrix(a:b, :)
            END IF
            IF (ALLOCATED(catalogue%points(q)%coreps%restriction)) THEN
               right(:, :) = MATMUL(catalogue%points(q)%coreps%restriction, catalogue%matrix(c:d, :))
            ELSE
               right(:, :) = catalogue%matrix(c:d, :)
            END IF
            ALLOCATE (lc(SIZE(edge%left_branch, 1), SIZE(left, 2)), &
                      rc(SIZE(edge%right_branch, 1), SIZE(right, 2)))
            lc(:, :) = MATMUL(edge%left_branch, left)
            rc(:, :) = MATMUL(edge%right_branch, right)
            IF (ANY(lc /= rc)) ERROR STOP 'Atomic column violates automatic unitary compatibility'
            DEALLOCATE (lc, rc)
            IF (ALLOCATED(edge%left_corep_branch)) THEN
               ALLOCATE (lc(SIZE(edge%left_corep_branch, 1), SIZE(left, 2)), &
                         rc(SIZE(edge%right_corep_branch, 1), SIZE(right, 2)))
               lc(:, :) = MATMUL(edge%left_corep_branch, catalogue%matrix(a:b, :))
               rc(:, :) = MATMUL(edge%right_corep_branch, catalogue%matrix(c:d, :))
               IF (ANY(lc /= rc)) ERROR STOP 'Atomic column violates automatic antiunitary compatibility'
               DEALLOCATE (lc, rc)
            END IF
            DEALLOCATE (left, right)
            comparisons = comparisons + SIZE(catalogue%matrix, 2)
         END DO
      END DO
   END SUBROUTINE sampling_test

! **************************************************************************************************
!> \brief Check atomic induction, inverse maps and all alternative operations on every retained star link.
!> \param catalogue Atomic reference and endpoint tables
!> \param lattice Joint graph with star identifications
!> \param rotations Original full-group rotations
!> \param anti Original operation time-reversal flags
! **************************************************************************************************
   SUBROUTINE check_star_maps(catalogue, lattice, rotations, anti)
      TYPE(atomic_signature_type), INTENT(IN)            :: catalogue
      TYPE(compatibility_lattice_type), INTENT(IN)       :: lattice
      INTEGER, INTENT(IN)                                :: rotations(:, :, :)
      LOGICAL, INTENT(IN)                                :: anti(:)

      INTEGER                                            :: a, b, c, d, e, g, h, info, p, q, shift(3)
      INTEGER, ALLOCATABLE                               :: other(:, :), permutation(:, :)
      REAL(dp)                                           :: delta(3), error

      DO e = 1, SIZE(lattice%stars, 2)
         p = lattice%stars(1, e)
         q = lattice%stars(2, e)
         g = lattice%stars(3, e)
         CALL atomic_signature_star_map(catalogue, p, q, g, permutation, 1.e-7_dp, error, info)
         IF (info /= 0) ERROR STOP 'Retained star map failed'
         a = catalogue%row_begin(p)
         b = catalogue%row_begin(p + 1) - 1
         c = catalogue%row_begin(q)
         d = catalogue%row_begin(q + 1) - 1
         IF (ANY(MATMUL(permutation, catalogue%matrix(a:b, :)) /= catalogue%matrix(c:d, :))) THEN
            ERROR STOP 'Star transport disagrees with independent atomic induction'
         END IF
         CALL star_operation(catalogue%group, catalogue%points(q)%little%k, catalogue%points(p)%little%k, &
                             g, shift, 1.e-7_dp, info)
         IF (info /= 0) ERROR STOP 'Inverse star operation missing'
         CALL atomic_signature_star_map(catalogue, q, p, g, other, 1.e-7_dp, error, info)
         IF (info /= 0) ERROR STOP 'Inverse star map failed'
         IF (ANY(other /= TRANSPOSE(permutation))) ERROR STOP 'Forward/reverse star maps do not cancel'
         DO h = 1, SIZE(anti)
            delta = MERGE(-1.0_dp, 1.0_dp, anti(h))* &
                    MATMUL(TRANSPOSE(inv_3x3(REAL(rotations(:, :, h), dp))), catalogue%points(p)%little%k) - &
                    catalogue%points(q)%little%k
            IF (MAXVAL(ABS(delta - ANINT(delta))) > 1.e-8_dp) CYCLE
            CALL atomic_signature_star_map(catalogue, p, q, h, other, 1.e-7_dp, error, info)
            IF (info /= 0) ERROR STOP 'Alternative star operation failed'
            IF (ANY(other /= permutation)) ERROR STOP 'Equivalent star operations give different corepresentation maps'
         END DO
      END DO
   END SUBROUTINE check_star_maps

END PROGRAM topology_reciprocal_unittest
