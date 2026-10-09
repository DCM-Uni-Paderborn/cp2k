!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Wigner classification and restrictions of projective antiunitary corepresentations.
!> \author Thomas D. Kuehne, tkuehne@cp2k.org
! **************************************************************************************************
MODULE topology_corepresentations
   USE kinds,                           ONLY: dp
   USE topology_little_group,           ONLY: validate_factor_system
   USE topology_symmetry,               ONLY: character_multiplicities

   IMPLICIT NONE
   PRIVATE

   TYPE, PUBLIC :: corepresentation_table_type
      INTEGER, ALLOCATABLE :: restriction(:, :), indicator(:), DIMENSION(:)
   END TYPE corepresentation_table_type

   PUBLIC :: projective_corepresentations, corepresentation_multiplicities
CONTAINS

! **************************************************************************************************
!> \brief Classify irreducible corepresentations through their unitary restrictions and Wigner test.
!> \param product Full group table with identity first
!> \param factor Normalized twisted factor system, including spin and nonsymmorphic phases
!> \param antiunitary Flags defining the unitary subgroup and its antiunitary coset
!> \param characters Validated complete unitary character table, rows in ascending full-group order
!> \param table Corepresentation dimensions, indicators (+1,-1,0), and unitary restrictions
!> \param tolerance Absolute factor/character tolerance
!> \param residual Largest factor, indicator or conjugate-character residual
!> \param status Zero on success, negative for invalid or inconsistent data
!> \note No matrix gauge or conventional irrep labels are assigned. Requires an antiunitary coset.
!>       Type +1 retains an irrep, -1 doubles it, and 0 pairs two inequivalent conjugate irreps.
!>       These signs are not eigenvalues of time reversal or insulating topological indices.
! **************************************************************************************************
   SUBROUTINE projective_corepresentations(product, factor, antiunitary, characters, table, tolerance, residual, status)
      INTEGER, INTENT(IN)                                :: product(:, :)
      COMPLEX(dp), INTENT(IN)                            :: factor(:, :)
      LOGICAL, INTENT(IN)                                :: antiunitary(:)
      COMPLEX(dp), INTENT(IN)                            :: characters(:, :)
      TYPE(corepresentation_table_type), INTENT(OUT)     :: table
      REAL(dp), INTENT(IN)                               :: tolerance
      REAL(dp), INTENT(OUT)                              :: residual
      INTEGER, INTENT(OUT)                               :: status

      COMPLEX(dp)                                        :: value
      COMPLEX(dp), ALLOCATABLE                           :: conjugated(:), regular(:)
      INTEGER                                            :: a, ainv, g, h, i, j, n, nc, ni, nu, q
      INTEGER, ALLOCATABLE                               :: branches(:, :), dimensions(:), &
                                                            indicators(:), multiplicities(:), &
                                                            partner(:), unitary(:), &
                                                            unitary_index(:)

      CALL validate_factor_system(product, factor, antiunitary, tolerance, residual, status)
      IF (status /= 0) RETURN
      status = -1
      n = SIZE(antiunitary)
      nu = COUNT(.NOT. antiunitary)
      ni = SIZE(characters, 2)
      IF (n /= 2*nu .OR. SIZE(characters, 1) /= nu .OR. ni < 1) RETURN
      ALLOCATE (unitary(nu), unitary_index(n), regular(nu), conjugated(nu), multiplicities(ni))
      unitary = PACK([(g, g=1, n)],.NOT. antiunitary)
      unitary_index = 0
      DO g = 1, nu
         unitary_index(unitary(g)) = g
      END DO
      regular = 0.0_dp
      regular(1) = nu
      CALL character_multiplicities(regular, characters, multiplicities, tolerance, status)
      IF (status /= 0) RETURN
      status = -2
      ALLOCATE (dimensions(ni), partner(ni), indicators(ni), branches(ni, ni))
      dimensions = NINT(REAL(characters(1, :), dp))
      IF (ANY(multiplicities /= dimensions)) RETURN
      a = FINDLOC(antiunitary, .TRUE., DIM=1)
      ainv = FINDLOC(PRODUCT(a, :), 1, DIM=1)
      status = -2
      DO i = 1, ni
         value = 0.0_dp
         DO g = 1, n
            IF (antiunitary(g)) value = value + factor(g, g)*characters(unitary_index(PRODUCT(g, g)), i)
         END DO
         value = value/REAL(nu, dp)
         q = NINT(REAL(value, dp))
         residual = MAX(residual, ABS(value - q))
         IF (ABS(q) > 1 .OR. residual > tolerance) RETURN
         indicators(i) = q
         ! h acting on a-transformed states is omega(h,a)/omega(a,a^-1 h a) times conjugate(D(a^-1 h a)).
         DO g = 1, nu
            h = PRODUCT(ainv, PRODUCT(unitary(g), a))
            conjugated(g) = factor(unitary(g), a)/factor(a, h)*CONJG(characters(unitary_index(h), i))
         END DO
         CALL character_multiplicities(conjugated, characters, multiplicities, tolerance, status)
         IF (status /= 0) RETURN
         status = -2
         IF (SUM(multiplicities) /= 1) RETURN
         j = FINDLOC(multiplicities, 1, DIM=1)
         partner(i) = j
         residual = MAX(residual, MAXVAL(ABS(conjugated - characters(:, j))))
         IF (residual > tolerance .OR. dimensions(i) /= dimensions(j)) RETURN
         IF ((q == 0) .EQV. (i == j)) RETURN
      END DO
      DO i = 1, ni
         j = partner(i)
         IF (partner(j) /= i .OR. indicators(i) /= indicators(j)) RETURN
      END DO
      branches = 0
      nc = 0
      DO i = 1, ni
         j = partner(i)
         IF (j < i) CYCLE
         nc = nc + 1
         SELECT CASE (indicators(i))
         CASE (1)
            branches(i, nc) = 1
         CASE (-1)
            branches(i, nc) = 2
         CASE (0)
            branches(i, nc) = 1
            branches(j, nc) = 1
         END SELECT
      END DO
      ALLOCATE (table%restriction(ni, nc), table%indicator(nc), table%dimension(nc))
      table%restriction = branches(:, :nc)
      table%dimension = MATMUL(dimensions, table%restriction)
      nc = 0
      DO i = 1, ni
         IF (partner(i) < i) CYCLE
         nc = nc + 1
         table%indicator(nc) = indicators(i)
      END DO
      status = 0
   END SUBROUTINE projective_corepresentations

! **************************************************************************************************
!> \brief Convert unitary irrep counts to nonnegative corepresentation counts without rounding.
!> \param table Validated irreducible corepresentation restrictions
!> \param unitary_counts Selected-band unitary multiplicities
!> \param counts Corepresentation multiplicities, invalid on failure
!> \param status Zero on success; -1 malformed table; -2 violated pairing or doubling
! **************************************************************************************************
   SUBROUTINE corepresentation_multiplicities(table, unitary_counts, counts, status)
      TYPE(corepresentation_table_type), INTENT(IN)      :: table
      INTEGER, INTENT(IN)                                :: unitary_counts(:)
      INTEGER, INTENT(OUT)                               :: counts(:), status

      INTEGER                                            :: i, j, n, ni, q

      status = -1
      counts = -1
      IF (.NOT. ALLOCATED(table%restriction) .OR. .NOT. ALLOCATED(table%indicator)) RETURN
      ni = SIZE(unitary_counts)
      n = SIZE(counts)
      IF (n < 1 .OR. ni < 1 .OR. SIZE(table%indicator) /= n) RETURN
      IF (ANY(SHAPE(table%restriction) /= [ni, n])) RETURN
      IF (ANY(unitary_counts < 0)) RETURN
      IF (ANY(table%restriction < 0) .OR. ANY(table%restriction > 2)) RETURN
      IF (ANY(COUNT(table%restriction > 0, DIM=2) /= 1)) RETURN
      DO j = 1, n
         SELECT CASE (table%indicator(j))
         CASE (1)
            IF (SUM(table%restriction(:, j)) /= 1) RETURN
         CASE (-1)
            IF (COUNT(table%restriction(:, j) == 2) /= 1 .OR. SUM(table%restriction(:, j)) /= 2) RETURN
         CASE (0)
            IF (COUNT(table%restriction(:, j) == 1) /= 2 .OR. SUM(table%restriction(:, j)) /= 2) RETURN
         CASE DEFAULT
            RETURN
         END SELECT
      END DO
      status = -2
      DO j = 1, n
         DO i = 1, ni
            q = table%restriction(i, j)
            IF (q == 0) CYCLE
            IF (MODULO(unitary_counts(i), q) /= 0) RETURN
            IF (counts(j) >= 0 .AND. counts(j) /= unitary_counts(i)/q) RETURN
            counts(j) = unitary_counts(i)/q
         END DO
      END DO
      status = 0
   END SUBROUTINE corepresentation_multiplicities
END MODULE topology_corepresentations
