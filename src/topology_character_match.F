!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Character-level correspondence with explicit operation and projective-gauge provenance.
!> \author Thomas D. Kuehne, tkuehne@cp2k.org
! **************************************************************************************************
MODULE topology_character_match
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: dp
   USE topology_little_group,           ONLY: validate_factor_system

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: character_table_match
CONTAINS

! **************************************************************************************************
!> \brief Match complete character tables after checking the supplied group and factor-gauge maps.
!> \param product Target group product table, identity first
!> \param factor Target normalized twisted factor system
!> \param anti Target antiunitary flags
!> \param characters Target irreps or irreducible corepresentation restrictions on the unitary rows
!> \param reference_product Reference group product table, identity first
!> \param reference_factor Reference normalized twisted factor system
!> \param reference_anti Reference antiunitary flags
!> \param reference_characters Reference characters, unitary rows in ascending reference-group order
!> \param operations Reference operation for every target operation
!> \param phases Explicit lift convention D_target(g)=phase(g)*D_reference(operations(g))
!> \param mapping Reference column for every target column; unallocated on failure
!> \param tolerance Absolute algebra and dimension-scaled character tolerance
!> \param residual Largest algebra, orthogonality, rank-sum or character residual
!> \param status Zero on success; -1 malformed input; -2 algebra mismatch; -3 character mismatch
!> \note The spatial meaning of the operation map must be validated by the caller. Central
!>       projectors certify the unitary characters; in the antiunitary case the caller must
!>       additionally establish that the restrictions extend to irreducible corepresentations.
!>       A factor coboundary does not uniquely fix labels: one-dimensional character twists also
!>       preserve the factors. Phases are therefore supplied explicitly, never inferred here.
! **************************************************************************************************
   SUBROUTINE character_table_match(product, factor, anti, characters, reference_product, reference_factor, &
                                    reference_anti, reference_characters, operations, phases, mapping, &
                                    tolerance, residual, status)
      INTEGER, INTENT(IN)                                :: product(:, :)
      COMPLEX(dp), INTENT(IN)                            :: factor(:, :)
      LOGICAL, INTENT(IN)                                :: anti(:)
      COMPLEX(dp), INTENT(IN)                            :: characters(:, :)
      INTEGER, INTENT(IN)                                :: reference_product(:, :)
      COMPLEX(dp), INTENT(IN)                            :: reference_factor(:, :)
      LOGICAL, INTENT(IN)                                :: reference_anti(:)
      COMPLEX(dp), INTENT(IN)                            :: reference_characters(:, :)
      INTEGER, INTENT(IN)                                :: operations(:)
      COMPLEX(dp), INTENT(IN)                            :: phases(:)
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: mapping(:)
      REAL(dp), INTENT(IN)                               :: tolerance
      REAL(dp), INTENT(OUT)                              :: residual
      INTEGER, INTENT(OUT)                               :: status

      COMPLEX(dp)                                        :: expected, phase
      COMPLEX(dp), ALLOCATABLE                           :: aligned(:, :)
      INTEGER                                            :: g, h, i, j, n, ni, nu, q
      INTEGER, ALLOCATABLE                               :: indices(:), reference_rows(:), unitary(:)
      REAL(dp)                                           :: closest, error, scale

      status = -1
      residual = HUGE(1.0_dp)
      n = SIZE(anti)
      IF (n < 1 .OR. n > 1024 .OR. SIZE(reference_anti) /= n) RETURN
      IF (SIZE(operations) /= n .OR. SIZE(phases) /= n) RETURN
      IF (ANY(operations < 1) .OR. ANY(operations > n)) RETURN
      DO i = 1, n
         IF (COUNT(operations == i) /= 1) RETURN
      END DO
      IF (operations(1) /= 1 .OR. ANY(anti .NEQV. reference_anti(operations))) RETURN
      IF (.NOT. ALL(ieee_is_finite(REAL(phases, dp))) .OR. .NOT. ALL(ieee_is_finite(AIMAG(phases)))) RETURN
      IF (MAXVAL(ABS(REAL(phases, dp))) > 2 .OR. MAXVAL(ABS(AIMAG(phases))) > 2) RETURN
      IF (.NOT. ieee_is_finite(tolerance)) RETURN
      IF (tolerance <= 0 .OR. tolerance >= 0.01_dp) RETURN
      IF (MAXVAL(ABS(ABS(phases) - 1)) > tolerance .OR. ABS(phases(1) - 1) > tolerance) RETURN
      nu = COUNT(.NOT. anti)
      ni = SIZE(characters, 2)
      IF (ni < 1 .OR. ni > nu .OR. SIZE(characters, 1) /= nu) RETURN
      IF (ANY(SHAPE(reference_characters) /= [nu, ni])) RETURN
      IF (ANY(anti) .AND. n /= 2*nu) RETURN
      CALL validate_factor_system(product, factor, anti, tolerance, residual, status)
      IF (status /= 0) RETURN
      CALL validate_factor_system(reference_product, reference_factor, reference_anti, tolerance, error, status)
      IF (status /= 0) RETURN
      residual = MAX(residual, error)
      status = -2
      DO h = 1, n
         DO g = 1, n
            q = PRODUCT(g, h)
            IF (operations(q) /= reference_product(operations(g), operations(h))) RETURN
            phase = phases(h)
            IF (anti(g)) phase = CONJG(phase)
            expected = phases(g)*phase*CONJG(phases(q))*reference_factor(operations(g), operations(h))
            residual = MAX(residual, ABS(expected - factor(g, h)))
         END DO
      END DO
      IF (residual > tolerance) RETURN
      CALL table_checks(product, factor, anti, characters, tolerance, error, status)
      residual = MAX(residual, error)
      IF (status /= 0) RETURN
      CALL table_checks(reference_product, reference_factor, reference_anti, reference_characters, tolerance, error, status)
      residual = MAX(residual, error)
      IF (status /= 0) RETURN
      ALLOCATE (unitary(nu))
      unitary(:) = PACK([(g, g=1, n)],.NOT. anti)
      ALLOCATE (reference_rows(n), aligned(nu, ni), indices(ni))
      reference_rows = 0
      j = 0
      DO i = 1, n
         IF (reference_anti(i)) CYCLE
         j = j + 1
         reference_rows(i) = j
      END DO
      DO i = 1, nu
         g = unitary(i)
         aligned(i, :) = phases(g)*reference_characters(reference_rows(operations(g)), :)
      END DO
      indices = 0
      status = -3
      DO i = 1, ni
         scale = MAX(1.0_dp, ABS(characters(1, i)))
         closest = HUGE(1.0_dp)
         DO j = 1, ni
            error = MAXVAL(ABS(characters(:, i) - aligned(:, j)))/scale
            closest = MIN(closest, error)
            IF (error > tolerance) CYCLE
            IF (indices(i) /= 0) RETURN
            indices(i) = j
            residual = MAX(residual, error)
         END DO
         residual = MAX(residual, closest)
         IF (indices(i) == 0) RETURN
      END DO
      DO i = 1, ni
         IF (COUNT(indices == i) /= 1) RETURN
      END DO
      CALL MOVE_ALLOC(indices, mapping)
      status = 0
   END SUBROUTINE character_table_match

! **************************************************************************************************
!> \brief Check completeness, orthogonality and central projectors of unitary restrictions.
!> \param product Validated full group product table
!> \param factor Validated full group factor system
!> \param anti Antiunitary flags, permitting Wigner norms 2 and 4 in addition to norm 1
!> \param table Unitary characters or unitary restrictions of irreducible corepresentations
!> \param tolerance Absolute scaled tolerance
!> \param error Largest residual
!> \param status Zero on success, -3 inconsistent table
! **************************************************************************************************
   SUBROUTINE table_checks(product, factor, anti, table, tolerance, error, status)
      INTEGER, INTENT(IN)                                :: product(:, :)
      COMPLEX(dp), INTENT(IN)                            :: factor(:, :)
      LOGICAL, INTENT(IN)                                :: anti(:)
      COMPLEX(dp), INTENT(IN)                            :: table(:, :)
      REAL(dp), INTENT(IN)                               :: tolerance
      REAL(dp), INTENT(OUT)                              :: error
      INTEGER, INTENT(OUT)                               :: status

      COMPLEX(dp)                                        :: overlap
      COMPLEX(dp), ALLOCATABLE                           :: projector(:, :)
      INTEGER                                            :: d, g, h, i, j, n, norm, q
      INTEGER, ALLOCATABLE                               :: rows(:), unit_product(:, :), unitary(:)
      REAL(dp)                                           :: weight

      status = -3
      error = HUGE(1.0_dp)
      n = SIZE(table, 1)
      IF (.NOT. ALL(ieee_is_finite(REAL(table, dp))) .OR. .NOT. ALL(ieee_is_finite(AIMAG(table)))) RETURN
      IF (MAXVAL(ABS(REAL(table, dp))) > 2*n .OR. MAXVAL(ABS(AIMAG(table))) > 2*n) RETURN
      IF (MAXVAL(ABS(table)) > 2*n) RETURN
      ALLOCATE (rows(SIZE(anti)), unitary(n), unit_product(n, n), projector(n, n))
      unitary(:) = PACK([(g, g=1, SIZE(anti))],.NOT. anti)
      rows = 0
      rows(unitary) = [(g, g=1, n)]
      DO h = 1, n
         unit_product(:, h) = rows(PRODUCT(unitary, unitary(h)))
      END DO
      error = 0
      weight = 0
      DO i = 1, SIZE(table, 2)
         d = NINT(REAL(table(1, i), dp))
         IF (d < 1) RETURN
         error = MAX(error, ABS(table(1, i) - d)/d)
         overlap = SUM(CONJG(table(:, i))*table(:, i))/n
         norm = NINT(REAL(overlap, dp))
         IF (norm /= 1) THEN
            IF (.NOT. ANY(anti) .OR. (norm /= 2 .AND. norm /= 4)) RETURN
         END IF
         error = MAX(error, ABS(overlap - norm)/norm)
         weight = weight + REAL(d, dp)**2/norm
         DO j = 1, i - 1
            overlap = SUM(CONJG(table(:, j))*table(:, i))/n
            error = MAX(error, ABS(overlap))
         END DO
         ! Wigner normalization yields one projector also for paired/doubled restrictions.
         DO h = 1, n
            DO g = 1, n
               projector(unit_product(g, h), h) = &
                  REAL(d, dp)/(REAL(norm, dp)*n)*CONJG(table(g, i))*factor(unitary(g), unitary(h))
            END DO
         END DO
         error = MAX(error, MAXVAL(ABS(projector - CONJG(TRANSPOSE(projector)))))
         error = MAX(error, MAXVAL(ABS(MATMUL(projector, projector) - projector)))
         DO g = 1, n
            DO h = 1, n
               q = unit_product(g, h)
               error = MAX(error, MAXVAL(ABS(factor(unitary(g), unitary(h))*projector(h, :) - &
                                             projector(q, unit_product(g, :))*factor(unitary(g), unitary))))
            END DO
         END DO
         IF (error > tolerance) RETURN
      END DO
      error = MAX(error, ABS(weight/n - 1))
      IF (error > tolerance) RETURN
      status = 0
   END SUBROUTINE table_checks
END MODULE topology_character_match
