!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Transport projective little-group characters between symmetry-related reciprocal points.
!> \author Thomas D. Kuehne, tkuehne@cp2k.org
! **************************************************************************************************
MODULE topology_star_transport
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: twopi
   USE mathlib,                         ONLY: inv_3x3
   USE topology_compatibility,          ONLY: little_group_data_type
   USE topology_little_group,           ONLY: little_group_factor_system
   USE topology_symmetry,               ONLY: character_multiplicities

   IMPLICIT NONE
   PRIVATE

   TYPE, PUBLIC :: star_group_type
      PRIVATE
      REAL(dp) :: cell(3, 3) = 0.0_dp
      LOGICAL :: spinful = .FALSE.
      INTEGER, ALLOCATABLE :: rotations(:, :, :), PRODUCT(:, :), original(:)
      REAL(dp), ALLOCATABLE :: translations(:, :), reciprocal(:, :, :)
      LOGICAL, ALLOCATABLE :: antiunitary(:)
      COMPLEX(dp), ALLOCATABLE :: factors(:, :)
   END TYPE star_group_type
   PUBLIC :: star_group_construct, star_group_storage, star_operation, star_irrep_map
CONTAINS

! **************************************************************************************************
!> \brief Cache the validated full group and its spin factor at Gamma, retaining input-operation IDs.
!> \param cell Cartesian cell vectors as columns
!> \param rotations Full integer direct-space rotations modulo input-cell translations
!> \param translations Fractional coset representatives
!> \param antiunitary Physical time-reversal flags
!> \param spinful Include spin lifts and physical Theta squared
!> \param group Validated immutable group; empty on failure
!> \param tolerance Geometry/algebra tolerance
!> \param status Zero on success; negative invalid group/data
! **************************************************************************************************
   SUBROUTINE star_group_construct(cell, rotations, translations, antiunitary, spinful, group, tolerance, status)
      REAL(dp), INTENT(IN)                               :: cell(3, 3)
      INTEGER, INTENT(IN)                                :: rotations(:, :, :)
      REAL(dp), INTENT(IN)                               :: translations(:, :)
      LOGICAL, INTENT(IN)                                :: antiunitary(:), spinful
      TYPE(star_group_type), INTENT(OUT)                 :: group
      REAL(dp), INTENT(IN)                               :: tolerance
      INTEGER, INTENT(OUT)                               :: status

      INTEGER                                            :: g, n
      TYPE(star_group_type)                              :: work

      CALL little_group_factor_system(cell, rotations, translations, antiunitary, [0.0_dp, 0.0_dp, 0.0_dp], &
                                      spinful, work%original, work%product, work%factors, tolerance, status)
      IF (status /= 0) RETURN
      n = SIZE(work%original)
      status = -2
      IF (n /= SIZE(antiunitary)) RETURN
      work%cell = cell
      work%spinful = spinful
      work%rotations = rotations(:, :, work%original)
      work%translations = translations(:, work%original)
      work%antiunitary = antiunitary(work%original)
      ALLOCATE (work%reciprocal(3, 3, n))
      DO g = 1, n
         work%reciprocal(:, :, g) = MERGE(-1.0_dp, 1.0_dp, work%antiunitary(g))* &
                                    TRANSPOSE(inv_3x3(REAL(work%rotations(:, :, g), dp)))
      END DO
      group = work
      status = 0
   END SUBROUTINE star_group_construct

! **************************************************************************************************
!> \brief Estimated retained group storage, excluding descriptors.
!> \param group Cached full group
!> \return Storage in MiB
! **************************************************************************************************
   REAL(dp) FUNCTION star_group_storage(group) RESULT(memory)
      TYPE(star_group_type), INTENT(IN)                  :: group

      REAL(dp)                                           :: n

      memory = 0.0_dp
      IF (.NOT. ALLOCATED(group%original)) RETURN
      n = REAL(SIZE(group%original), dp)
      memory = (20*n*n + 140*n + 72)/1024.0_dp**2
   END FUNCTION star_group_storage

! **************************************************************************************************
!> \brief Find one operation taking source to target modulo a reciprocal lattice vector.
!> \param group Validated full group
!> \param source Fractional reciprocal source
!> \param TARGET Fractional reciprocal destination
!> \param operation Original input-operation index, or zero when not related
!> \param shift Integer vector satisfying sigma_s W_s^-T source = target + shift
!> \param tolerance Fractional tolerance
!> \param status 0 related; 1 distinct stars; -1 invalid data; -3 unsafe integer shift
! **************************************************************************************************
   SUBROUTINE star_operation(group, source, TARGET, operation, shift, tolerance, status)
      TYPE(star_group_type), INTENT(IN)                  :: group
      REAL(dp), INTENT(IN)                               :: source(3), TARGET(3)
      INTEGER, INTENT(OUT)                               :: operation, shift(3)
      REAL(dp), INTENT(IN)                               :: tolerance
      INTEGER, INTENT(OUT)                               :: status

      INTEGER                                            :: g
      REAL(dp)                                           :: delta(3)

      operation = 0
      shift = 0
      status = -1
      IF (.NOT. ALLOCATED(group%original)) RETURN
      IF (.NOT. ieee_is_finite(tolerance)) RETURN
      IF (tolerance <= 0.0_dp .OR. tolerance >= 0.01_dp) RETURN
      IF (.NOT. ALL(ieee_is_finite(source)) .OR. .NOT. ALL(ieee_is_finite(TARGET))) RETURN
      IF (MAXVAL(ABS(source)) > tolerance/(256*EPSILON(1.0_dp))/MAX(1.0_dp, MAXVAL(ABS(group%reciprocal)))) RETURN
      IF (MAXVAL(ABS(TARGET)) > tolerance/(256*EPSILON(1.0_dp))) RETURN
      DO g = 1, SIZE(group%original)
         delta = MATMUL(group%reciprocal(:, :, g), source) - TARGET
         IF (MAXVAL(ABS(delta - ANINT(delta))) > tolerance) CYCLE
         status = -3
         IF (MAXVAL(ABS(delta)) >= REAL(HUGE(0), dp)) RETURN
         operation = group%original(g)
         shift = NINT(delta)
         status = 0
         RETURN
      END DO
      status = 1
   END SUBROUTINE star_operation

! **************************************************************************************************
!> \brief Bijection of unitary irreps under a specified unitary/antiunitary star operation.
!> \param group Validated full group with Gamma spin factors
!> \param left Validated complete little group at the source
!> \param right Validated complete little group at the destination
!> \param operation Original input-operation index taking left k to right k modulo integers
!> \param permutation Destination-by-source unitary-irrep permutation, empty on failure
!> \param tolerance Fractional, factor and character tolerance
!> \param residual Largest conjugation, factor or character residual
!> \param status 0 success; -1 invalid data; -2 inconsistent transport
!> \note h*s and s*(s^-1*h*s) differ by a lattice translation L. The character multiplier is
!>       z(h,s)/z(s,s^-1*h*s)*exp(-2*pi*i*k_right.L), with source-character conjugation for anti_s.
!>       The full twisted factor law is checked, including antiunitary little-group operations.
! **************************************************************************************************
   SUBROUTINE star_irrep_map(group, left, right, operation, permutation, tolerance, residual, status)
      TYPE(star_group_type), INTENT(IN)                  :: group
      TYPE(little_group_data_type), INTENT(IN)           :: left, right
      INTEGER, INTENT(IN)                                :: operation
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: permutation(:, :)
      REAL(dp), INTENT(IN)                               :: tolerance
      REAL(dp), INTENT(OUT)                              :: residual
      INTEGER, INTENT(OUT)                               :: status

      COMPLEX(dp)                                        :: factor, multiplier
      COMPLEX(dp), ALLOCATABLE                           :: beta(:), chi(:)
      INTEGER                                            :: g, gh, h, i, j, n, ni, s, si, source, &
                                                            source_row
      INTEGER, ALLOCATABLE                               :: conjugate_index(:), il(:), ir(:), &
                                                            unit_left(:), unit_right(:), work(:, :)
      REAL(dp)                                           :: angle, delta(3)

      status = -1
      residual = HUGE(1.0_dp)
      IF (.NOT. ALLOCATED(group%original)) RETURN
      IF (.NOT. ieee_is_finite(tolerance)) RETURN
      IF (tolerance <= 0.0_dp .OR. tolerance >= 0.01_dp) RETURN
      s = FINDLOC(group%original, operation, DIM=1)
      IF (s == 0) RETURN
      CALL point_indices(group, left, il, tolerance, status)
      IF (status /= 0) RETURN
      CALL point_indices(group, right, ir, tolerance, status)
      IF (status /= 0) RETURN
      status = -2
      n = SIZE(il)
      IF (n /= SIZE(ir)) RETURN
      IF (MAXVAL(ABS(left%k)) > tolerance/(256*EPSILON(1.0_dp))/MAX(1.0_dp, MAXVAL(ABS(group%reciprocal)))) RETURN
      IF (MAXVAL(ABS(right%k)) > tolerance/(256*EPSILON(1.0_dp))) RETURN
      delta = MATMUL(group%reciprocal(:, :, s), left%k) - right%k
      residual = MAXVAL(ABS(delta - ANINT(delta)))
      IF (residual > tolerance) RETURN
      si = FINDLOC(group%product(:, s), 1, DIM=1)
      ALLOCATE (conjugate_index(n), beta(n), unit_left(n))
      unit_left = 0
      DO g = 1, n
         IF (.NOT. left%antiunitary(g)) unit_left(g) = COUNT(.NOT. left%antiunitary(:g))
      END DO
      unit_right = PACK([(g, g=1, n)],.NOT. right%antiunitary)
      DO g = 1, n
         h = ir(g)
         source = group%product(si, group%product(h, s))
         i = FINDLOC(il, source, DIM=1)
         IF (i == 0) RETURN
         conjugate_index(g) = i
         delta = right%translations(:, g) + MATMUL(REAL(right%rotations(:, :, g), dp), group%translations(:, s)) - &
                 group%translations(:, s) - MATMUL(REAL(group%rotations(:, :, s), dp), left%translations(:, i))
         residual = MAX(residual, MAXVAL(ABS(delta - ANINT(delta))))
         IF (residual > tolerance) RETURN
         IF (MAXVAL(ABS(delta)) > tolerance/(256*EPSILON(1.0_dp))/MAX(1.0_dp, MAXVAL(ABS(right%k)))) RETURN
         angle = -twopi*MODULO(DOT_PRODUCT(right%k, ANINT(delta)), 1.0_dp)
         beta(g) = group%factors(h, s)/group%factors(s, source)*CMPLX(COS(angle), SIN(angle), dp)
      END DO
      DO h = 1, n
         DO g = 1, n
            gh = right%product(g, h)
            IF (left%product(conjugate_index(g), conjugate_index(h)) /= conjugate_index(gh)) RETURN
            multiplier = beta(h)
            IF (right%antiunitary(g)) multiplier = CONJG(multiplier)
            factor = left%factors(conjugate_index(g), conjugate_index(h))
            IF (group%antiunitary(s)) factor = CONJG(factor)
            factor = beta(g)*multiplier/beta(gh)*factor
            residual = MAX(residual, ABS(factor - right%factors(g, h)))
         END DO
      END DO
      IF (residual > tolerance) RETURN
      ni = SIZE(left%characters, 2)
      IF (ni /= SIZE(right%characters, 2)) RETURN
      ALLOCATE (work(ni, ni), chi(SIZE(unit_right)))
      DO i = 1, ni
         DO j = 1, SIZE(unit_right)
            g = unit_right(j)
            source_row = unit_left(conjugate_index(g))
            IF (source_row == 0) RETURN
            multiplier = left%characters(source_row, i)
            IF (group%antiunitary(s)) multiplier = CONJG(multiplier)
            chi(j) = beta(g)*multiplier
         END DO
         CALL character_multiplicities(chi, right%characters, work(:, i), tolerance, status)
         IF (status /= 0) RETURN
         status = -2
         IF (SUM(work(:, i)) /= 1) RETURN
         j = FINDLOC(work(:, i), 1, DIM=1)
         residual = MAX(residual, MAXVAL(ABS(chi - right%characters(:, j))))
      END DO
      IF (ANY(SUM(work, DIM=2) /= 1) .OR. residual > tolerance) RETURN
      CALL MOVE_ALLOC(work, permutation)
      status = 0
   END SUBROUTINE star_irrep_map

! **************************************************************************************************
!> \brief Match endpoint cosets to the cached group, allowing integer representative changes.
!> \param group Validated full group
!> \param point Validated complete endpoint data
!> \param indices Matching Gamma-group indices
!> \param tolerance Geometry tolerance
!> \param status Zero on success
! **************************************************************************************************
   SUBROUTINE point_indices(group, point, indices, tolerance, status)
      TYPE(star_group_type), INTENT(IN)                  :: group
      TYPE(little_group_data_type), INTENT(IN)           :: point
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: indices(:)
      REAL(dp), INTENT(IN)                               :: tolerance
      INTEGER, INTENT(OUT)                               :: status

      INTEGER                                            :: g, h, n
      REAL(dp)                                           :: delta(3)

      status = -1
      IF (point%spinful .NEQV. group%spinful) RETURN
      IF (.NOT. ALL(ieee_is_finite(point%cell)) .OR. .NOT. ALL(ieee_is_finite(point%k))) RETURN
      IF (MAXVAL(ABS(point%cell - group%cell)) > tolerance*MAXVAL(ABS(group%cell))) RETURN
      IF (.NOT. ALLOCATED(point%product) .OR. .NOT. ALLOCATED(point%factors)) RETURN
      IF (.NOT. ALLOCATED(point%rotations) .OR. .NOT. ALLOCATED(point%translations)) RETURN
      IF (.NOT. ALLOCATED(point%antiunitary) .OR. .NOT. ALLOCATED(point%characters)) RETURN
      n = SIZE(point%antiunitary)
      IF (n < 1 .OR. ANY(SHAPE(point%product) /= [n, n]) .OR. ANY(SHAPE(point%factors) /= [n, n])) RETURN
      IF (ANY(SHAPE(point%rotations) /= [3, 3, n]) .OR. ANY(SHAPE(point%translations) /= [3, n])) RETURN
      IF (SIZE(point%characters, 1) /= COUNT(.NOT. point%antiunitary)) RETURN
      IF (SIZE(point%characters, 2) < 1) RETURN
      IF (ANY(point%product < 1) .OR. ANY(point%product > n)) RETURN
      IF (.NOT. ALL(ieee_is_finite(point%translations))) RETURN
      IF (MAXVAL(ABS(point%translations)) > tolerance/(256*EPSILON(1.0_dp)) &
          /MAX(1.0_dp, MAXVAL(ABS(REAL(group%rotations, dp))))) RETURN
      IF (.NOT. ALL(ieee_is_finite(REAL(point%factors, dp))) .OR. .NOT. ALL(ieee_is_finite(AIMAG(point%factors)))) RETURN
      IF (.NOT. ALL(ieee_is_finite(REAL(point%characters, dp))) .OR. &
          .NOT. ALL(ieee_is_finite(AIMAG(point%characters)))) RETURN
      ALLOCATE (indices(n))
      indices = 0
      DO g = 1, n
         DO h = 1, SIZE(group%original)
            IF (point%antiunitary(g) .NEQV. group%antiunitary(h)) CYCLE
            IF (ANY(point%rotations(:, :, g) /= group%rotations(:, :, h))) CYCLE
            delta = point%translations(:, g) - group%translations(:, h)
            IF (MAXVAL(ABS(delta - ANINT(delta))) > tolerance) CYCLE
            IF (indices(g) /= 0) RETURN
            indices(g) = h
         END DO
         IF (indices(g) == 0) RETURN
         IF (COUNT(indices == indices(g)) /= 1) RETURN
      END DO
      status = 0
   END SUBROUTINE point_indices
END MODULE topology_star_transport
