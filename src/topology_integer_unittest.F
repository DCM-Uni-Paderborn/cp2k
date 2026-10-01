!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
PROGRAM topology_integer_unittest
   USE kinds,                           ONLY: int_8
   USE topology_integer_lattice,        ONLY: &
        integer_lattice_certificate, integer_lattice_factor, integer_lattice_product, &
        integer_lattice_quotient, integer_lattice_solve, integer_lattice_type, &
        integer_quotient_class, integer_quotient_type, nonnegative_integer_solve
   USE topology_tqc,                    ONLY: inversion_ebr_signature

   IMPLICIT NONE

   CALL certificate_tests()
   CALL inversion_tests(4)
   CALL inversion_tests(8)
   CALL nonnegative_tests()
   CALL rejection_tests()
   CALL quotient_tests()
   PRINT *, 'Integer atomic-signature decomposition tests passed.'
CONTAINS

! **************************************************************************************************
!> \brief Quotient certificates, free factors, nonsaturated constraints and inversion oracles.
! **************************************************************************************************
   SUBROUTINE quotient_tests()
      INTEGER                                            :: center, i, info, j, k, n, row, sample, &
                                                            sign, value
      INTEGER(int_8), ALLOCATABLE                        :: a(:, :), b(:), c(:, :), indices(:), &
                                                            verified(:, :), witness(:)
      LOGICAL                                            :: compatible, exists
      TYPE(integer_lattice_type)                         :: atomic
      TYPE(integer_quotient_type)                        :: quotient

      c = RESHAPE([2_int_8, -4_int_8, 0_int_8], [1, 3])
      a = RESHAPE([4_int_8, 2_int_8, 0_int_8], [3, 1])
      CALL integer_lattice_quotient(c, a, quotient, info)
      IF (info /= 0) ERROR STOP 'Free quotient failed'
      IF (COUNT(quotient%invariants == 0) /= 1 .OR. COUNT(quotient%invariants == 2) /= 1) THEN
         ERROR STOP 'Expected Z plus Z2 quotient'
      END IF
      CALL integer_quotient_class(quotient, [2_int_8, 1_int_8, 0_int_8], compatible, indices, info)
      IF (info /= 0 .OR. .NOT. compatible) ERROR STOP 'Saturated kernel lost integer solution'
      IF (SUM(ABS(indices)) /= 1) ERROR STOP 'Nontrivial Z2 class lost'
      CALL integer_quotient_class(quotient, [0_int_8, 0_int_8, -3_int_8], compatible, indices, info)
      IF (info /= 0 .OR. .NOT. compatible) ERROR STOP 'Free target rejected'
      IF (SUM(ABS(indices)) /= 3) ERROR STOP 'Free quotient coordinate reduced modulo zero'
      CALL integer_quotient_class(quotient, [1_int_8, 0_int_8, 0_int_8], compatible, indices, info)
      IF (info /= 0 .OR. compatible .OR. ALLOCATED(indices)) ERROR STOP 'Incompatible signature classified'
      a(1, 1) = 3
      CALL integer_lattice_quotient(c, a, quotient, info)
      IF (info /= -2 .OR. ALLOCATED(quotient%kernel)) ERROR STOP 'Incompatible atomic generator accepted'
      c = RESHAPE([1_int_8, 0_int_8, 0_int_8, 1_int_8], [2, 2])
      a = RESHAPE([0_int_8, 0_int_8], [2, 1])
      CALL integer_lattice_quotient(c, a, quotient, info)
      IF (info /= 0 .OR. SIZE(quotient%invariants) /= 0) ERROR STOP 'Zero kernel failed'
      CALL integer_quotient_class(quotient, [0_int_8, 0_int_8], compatible, indices, info)
      IF (info /= 0 .OR. .NOT. compatible .OR. SIZE(indices) /= 0) ERROR STOP 'Zero quotient class failed'
      DEALLOCATE (c, a)
      DO n = 4, 8, 4
         ALLOCATE (c(n - 1, 2*n), a(2*n, 2*n), b(2*n))
         c = 0
         a = 0
         DO i = 2, n
            c(i - 1, :2) = 1
            c(i - 1, 2*i - 1:2*i) = -1
         END DO
         DO center = 0, n - 1
            DO sign = 0, 1
               j = 2*center + sign + 1
               DO k = 0, n - 1
                  row = 2*k + MODULO(sign + POPCNT(IAND(center, k)), 2) + 1
                  a(row, j) = 1
               END DO
            END DO
         END DO
         CALL integer_lattice_quotient(c, a, quotient, info)
         IF (info /= 0) ERROR STOP 'Inversion compatibility quotient failed'
         IF (SIZE(quotient%invariants) /= n + 1 .OR. ANY(quotient%invariants == 0)) THEN
            ERROR STOP 'Inversion compatibility rank wrong'
         END IF
         IF (n == 4) THEN
            IF (COUNT(quotient%invariants == 2) /= 1 .OR. COUNT(quotient%invariants == 1) /= 4) THEN
               ERROR STOP '2D inversion quotient not Z2'
            END IF
         ELSE
            IF (COUNT(quotient%invariants == 2) /= 3 .OR. COUNT(quotient%invariants == 4) /= 1) THEN
               ERROR STOP '3D inversion quotient not Z2 cubed times Z4'
            END IF
         END IF
         CALL integer_lattice_product(c, quotient%kernel, verified, info)
         IF (info /= 0 .OR. ANY(verified /= 0)) ERROR STOP 'Quotient kernel certificate wrong'
         CALL integer_lattice_product(quotient%kernel, quotient%coordinates, verified, info)
         IF (info /= 0 .OR. ANY(verified /= a)) ERROR STOP 'Quotient generator certificate wrong'
         CALL integer_lattice_factor(a, atomic, info)
         IF (info /= 0) ERROR STOP 'Quotient oracle factorization failed'
         DO sample = 0, 3**n - 1
            value = sample
            DO k = 1, n
               b(2*k) = MODULO(value, 3)
               b(2*k - 1) = 2 - b(2*k)
               value = value/3
            END DO
            CALL integer_lattice_solve(atomic, b, exists, witness, info)
            IF (info /= 0) ERROR STOP 'Quotient signed membership oracle failed'
            CALL integer_quotient_class(quotient, b, compatible, indices, info)
            IF (info /= 0 .OR. .NOT. compatible) ERROR STOP 'Equal-rank inversion signature incompatible'
            IF (ALL(indices == 0) .NEQV. exists) ERROR STOP 'Quotient zero class disagrees with signed membership'
         END DO
         DEALLOCATE (c, a, b)
      END DO
      PRINT *, 'Compatibility quotient tests passed.'
   END SUBROUTINE quotient_tests

! **************************************************************************************************
!> \brief Check Smith invariants against gcds of all minors, including singular rectangular matrices.
! **************************************************************************************************
   SUBROUTINE certificate_tests()
      INTEGER                                            :: cols, i, ii, info, j, jj, k, m, n, rank, &
                                                            rows, sample
      INTEGER(int_8)                                     :: counts(0:4), divisors, expected, &
                                                            minor(4, 4), seed
      INTEGER(int_8), ALLOCATABLE                        :: a(:, :), b(:), diagonal(:), left(:, :), &
                                                            PRODUCT(:, :), right(:, :), &
                                                            witness(:), x(:)
      LOGICAL                                            :: exists
      TYPE(integer_lattice_type)                         :: factor

      seed = 19283
      counts = 0
      DO m = 1, 4
         DO n = 1, 4
            ALLOCATE (a(m, n), PRODUCT(m, n), x(n), b(m))
            DO sample = 1, 24
               DO j = 1, n
                  DO i = 1, m
                     seed = MODULO(48271_int_8*seed, 2147483647_int_8)
                     a(i, j) = MODULO(seed, 9_int_8) - 4
                  END DO
               END DO
               IF (MODULO(sample, 4) == 0) a(m, :) = 0
               IF (MODULO(sample, 6) == 0) a(:, n) = a(:, 1)
               IF (sample == 24) a = 0
               CALL integer_lattice_factor(a, factor, info)
               IF (info /= 0) ERROR STOP 'Small integer factorization failed'
               CALL integer_lattice_certificate(factor, left, right, diagonal, rank)
               PRODUCT(:, :) = MATMUL(MATMUL(left, a), right)
               DO k = 1, MIN(m, n)
                  PRODUCT(k, k) = PRODUCT(k, k) - diagonal(k)
               END DO
               IF (ANY(product /= 0)) ERROR STOP 'Smith certificate equation failed'
               IF (ABS(determinant(left)) /= 1 .OR. ABS(determinant(right)) /= 1) THEN
                  ERROR STOP 'Smith transformations not unimodular'
               END IF
               IF (COUNT(diagonal /= 0) /= rank .OR. ANY(diagonal < 0)) ERROR STOP 'Smith rank/sign incorrect'
               counts(rank) = counts(rank) + 1
               expected = 1
               DO k = 1, MIN(m, n)
                  divisors = 0
                  DO rows = 1, 2**m - 1
                     IF (POPCNT(rows) /= k) CYCLE
                     DO cols = 1, 2**n - 1
                        IF (POPCNT(cols) /= k) CYCLE
                        ii = 0
                        DO i = 1, m
                           IF (.NOT. BTEST(rows, i - 1)) CYCLE
                           ii = ii + 1
                           jj = 0
                           DO j = 1, n
                              IF (.NOT. BTEST(cols, j - 1)) CYCLE
                              jj = jj + 1
                              minor(ii, jj) = a(i, j)
                           END DO
                        END DO
                        divisors = gcd(divisors, ABS(determinant(minor(:k, :k))))
                     END DO
                  END DO
                  expected = expected*diagonal(k)
                  IF (divisors /= expected) ERROR STOP 'Smith invariant disagrees with determinantal divisor'
                  IF (k > 1 .AND. k <= rank) THEN
                     IF (MODULO(diagonal(k), diagonal(k - 1)) /= 0) ERROR STOP 'Smith divisibility failed'
                  END IF
               END DO
               x = [(INT(MODULO(i + sample, 5) - 2, int_8), i=1, n)]
               b(:) = MATMUL(a, x)
               CALL integer_lattice_solve(factor, b, exists, witness, info)
               IF (info /= 0 .OR. .NOT. exists) ERROR STOP 'Known signed witness not found'
               IF (ANY(MATMUL(a, witness) /= b)) ERROR STOP 'Signed witness incorrect'
            END DO
            DEALLOCATE (a, product, x, b)
         END DO
      END DO
      IF (ANY(counts == 0)) ERROR STOP 'Missing rank coverage'
      PRINT *, 'Smith minor/unimodularity certificates, ranks 0 through 4:', counts
   END SUBROUTINE certificate_tests

! **************************************************************************************************
!> \brief Independent inversion Fourier oracle, all odd-count patterns with two Kramers pairs.
!> \param n Number of TRIM, four or eight
! **************************************************************************************************
   SUBROUTINE inversion_tests(n)
      INTEGER, INTENT(IN)                                :: n

      INTEGER                                            :: center, column, info, k, parity, &
                                                            positive_count, row, sample, sign, &
                                                            signed_count, value
      INTEGER(int_8)                                     :: total_nodes, visited
      INTEGER(int_8), ALLOCATABLE                        :: witness(:)
      INTEGER(int_8)                                     :: b(2*n), a(2*n, 2*n)
      INTEGER                                            :: reference(2, n), odd(n)
      LOGICAL                                            :: exists, nonnegative_atomic, signed_atomic
      TYPE(integer_lattice_type)                         :: factor

      a = 0
      DO center = 0, n - 1
         DO sign = 0, 1
            column = 2*center + sign + 1
            DO k = 0, n - 1
               parity = MODULO(sign + POPCNT(IAND(center, k)), 2)
               row = 2*k + parity + 1
               a(row, column) = 1
            END DO
         END DO
      END DO
      CALL integer_lattice_factor(a, factor, info)
      IF (info /= 0) ERROR STOP 'Inversion signature factorization failed'
      signed_count = 0
      positive_count = 0
      total_nodes = 0
      DO sample = 0, 3**n - 1
         value = sample
         DO k = 1, n
            odd(k) = MODULO(value, 3)
            value = value/3
            b(2*k - 1:2*k) = [2 - odd(k), odd(k)]
         END DO
         CALL inversion_ebr_signature(odd, 2, signed_atomic, nonnegative_atomic, reference, info)
         IF (info /= 0) ERROR STOP 'Inversion analytic oracle failed'
         CALL integer_lattice_solve(factor, b, exists, witness, info)
         IF (info /= 0 .OR. (exists .NEQV. signed_atomic)) ERROR STOP 'Signed inversion oracle disagreement'
         IF (exists) signed_count = signed_count + 1
         CALL nonnegative_integer_solve(a, b, exists, witness, visited, info, 100000_int_8)
         IF (info /= 0 .OR. (exists .NEQV. nonnegative_atomic)) THEN
            PRINT *, 'Inversion size/pattern/status/nodes:', n, sample, info, visited
            ERROR STOP 'Nonnegative inversion oracle disagreement'
         END IF
         IF (exists) THEN
            positive_count = positive_count + 1
            IF (ANY(witness < 0) .OR. ANY(MATMUL(a, witness) /= b)) ERROR STOP 'Inversion witness incorrect'
         END IF
         total_nodes = total_nodes + visited
      END DO
      PRINT *, 'Inversion TRIM/patterns/signed/nonnegative/search nodes:', n, 3**n, signed_count, positive_count, total_nodes
   END SUBROUTINE inversion_tests

! **************************************************************************************************
!> \brief Compare complete nonnegative search with exhaustive coefficient enumeration for tiny matrices.
! **************************************************************************************************
   SUBROUTINE nonnegative_tests()
      INTEGER                                            :: b1, b2, cases, i, info, j, sample, &
                                                            value, x1, x2, x3
      INTEGER(int_8)                                     :: a(2, 3), b(2), visited, x(3)
      INTEGER(int_8), ALLOCATABLE                        :: witness(:)
      LOGICAL                                            :: exists, reference(0:4, 0:4)

      cases = 0
      DO sample = 0, 3**6 - 1
         value = sample
         DO j = 1, 3
            DO i = 1, 2
               a(i, j) = MODULO(value, 3)
               value = value/3
            END DO
         END DO
         reference = .FALSE.
         DO x1 = 0, 4
            DO x2 = 0, 4
               DO x3 = 0, 4
                  x = [x1, x2, x3]
                  b = MATMUL(a, x)
                  IF (ANY(b > 4)) CYCLE
                  reference(b(1), b(2)) = .TRUE.
               END DO
            END DO
         END DO
         DO b1 = 0, 4
            DO b2 = 0, 4
               b = [b1, b2]
               CALL nonnegative_integer_solve(a, b, exists, witness, visited, info, 10000_int_8)
               IF (info /= 0 .OR. (exists .NEQV. reference(b1, b2))) ERROR STOP 'Exhaustive search disagreement'
               IF (exists) THEN
                  IF (ANY(witness < 0) .OR. ANY(MATMUL(a, witness) /= b)) ERROR STOP 'Nonnegative witness incorrect'
               ELSE
                  IF (ALLOCATED(witness)) ERROR STOP 'Witness allocated without a solution'
               END IF
               cases = cases + 1
            END DO
         END DO
      END DO
      PRINT *, 'Independent bounded-coefficient cases:', cases
   END SUBROUTINE nonnegative_tests

! **************************************************************************************************
!> \brief Distinguish overflow/work limits and malformed data from proven nonmembership.
! **************************************************************************************************
   SUBROUTINE rejection_tests()
      INTEGER(int_8), PARAMETER                          :: limit = HUGE(0_int_8)/16

      INTEGER                                            :: info, rank
      INTEGER(int_8)                                     :: a(2, 2), b(2), visited
      INTEGER(int_8), ALLOCATABLE                        :: diagonal(:), left(:, :), right(:, :), &
                                                            witness(:)
      LOGICAL                                            :: exists
      TYPE(integer_lattice_type)                         :: factor

      a = 0
      a(1, 1) = 2
      a(2, 2) = 4
      b = [1, 0]
      CALL integer_lattice_factor(a, factor, info)
      IF (info /= 0) ERROR STOP 'Diagonal factor failed'
      CALL integer_lattice_solve(factor, b, exists, witness, info)
      IF (info /= 0 .OR. exists .OR. ALLOCATED(witness)) ERROR STOP 'Exact nonmembership failed'
      CALL integer_lattice_solve(factor, b(:1), exists, witness, info)
      IF (info /= -1 .OR. exists) ERROR STOP 'Wrong target dimension accepted'
      a(2, 2) = 0
      CALL integer_lattice_factor(a, factor, info)
      IF (info /= 0) ERROR STOP 'Singular factor failed'
      b = [2, 1]
      CALL integer_lattice_solve(factor, b, exists, witness, info)
      IF (info /= 0 .OR. exists) ERROR STOP 'Rank-deficient nonmembership failed'
      a(2, 2) = 4
      b = [2, 4]
      CALL nonnegative_integer_solve(a, b, exists, witness, visited, info, 0_int_8)
      IF (info /= 1 .OR. exists .OR. visited /= 0 .OR. ALLOCATED(witness)) ERROR STOP 'Zero budget misclassified'
      CALL nonnegative_integer_solve(a, b, exists, witness, visited, info, 1_int_8)
      IF (info /= 1 .OR. exists .OR. visited /= 1 .OR. ALLOCATED(witness)) ERROR STOP 'Work budget misclassified'
      CALL nonnegative_integer_solve(a, [-1_int_8, 0_int_8], exists, witness, visited, info, 10_int_8)
      IF (info /= -1 .OR. exists) ERROR STOP 'Negative target accepted'
      CALL nonnegative_integer_solve(-a, b, exists, witness, visited, info, 10_int_8)
      IF (info /= -1 .OR. exists) ERROR STOP 'Negative signature matrix accepted'
      CALL nonnegative_integer_solve(a, b, exists, witness, visited, info, -1_int_8)
      IF (info /= -1 .OR. exists) ERROR STOP 'Invalid search budget accepted'
      CALL integer_lattice_factor(a(:0, :), factor, info)
      IF (info /= -1) ERROR STOP 'Empty matrix accepted'
      CALL integer_lattice_solve(factor, b, exists, witness, info)
      IF (info /= -1 .OR. exists) ERROR STOP 'Invalid factor accepted'
      a = 1
      a(1, 1) = limit
      a(2, 2) = limit
      CALL integer_lattice_factor(a, factor, info)
      IF (info /= -3) ERROR STOP 'Intermediate coefficient growth not detected'
      CALL integer_lattice_certificate(factor, left, right, diagonal, rank)
      IF (rank /= -1 .OR. ALLOCATED(left) .OR. ALLOCATED(right) .OR. ALLOCATED(diagonal)) THEN
         ERROR STOP 'Failed factorization leaked a partial certificate'
      END IF
      a = 0
      a(1, 1) = limit
      a(2, 2) = limit
      CALL integer_lattice_factor(a, factor, info)
      IF (info /= 0) ERROR STOP 'Large diagonal unnecessarily rejected'
      b = limit
      CALL integer_lattice_solve(factor, b, exists, witness, info)
      IF (info /= 0 .OR. .NOT. exists) ERROR STOP 'Large valid target failed'
      a(1, 1) = HUGE(0_int_8)
      CALL integer_lattice_factor(a, factor, info)
      IF (info /= -3) ERROR STOP 'Unsafe input bound accepted'
      PRINT *, 'Integer invalid-input, growth and unresolved-budget tests passed.'
   END SUBROUTINE rejection_tests

! **************************************************************************************************
!> \brief Independent small-matrix determinant by Laplace expansion.
!> \param a Square integer matrix, order at most four
!> \return Determinant
! **************************************************************************************************
   RECURSIVE FUNCTION determinant(a) RESULT(value)
      INTEGER(int_8), INTENT(IN)                         :: a(:, :)
      INTEGER(int_8)                                     :: value

      INTEGER                                            :: j, n
      INTEGER(int_8) :: minor(SIZE(a, 1) - 1, SIZE(a, 1) - 1)

      n = SIZE(a, 1)
      IF (n == 1) THEN
         value = a(1, 1)
         RETURN
      END IF
      value = 0
      DO j = 1, n
         minor(:, :j - 1) = a(2:, :j - 1)
         minor(:, j:) = a(2:, j + 1:)
         value = value + (-1)**(j + 1)*a(1, j)*determinant(minor)
      END DO
   END FUNCTION determinant

! **************************************************************************************************
!> \brief Euclidean gcd for the independent determinantal-divisor oracle.
!> \param a First nonnegative value
!> \param b Second nonnegative value
!> \return Nonnegative gcd
! **************************************************************************************************
   FUNCTION gcd(a, b) RESULT(value)
      INTEGER(int_8), INTENT(IN)                         :: a, b
      INTEGER(int_8)                                     :: value

      INTEGER(int_8)                                     :: next, remainder

      value = a
      next = b
      DO WHILE (next /= 0)
         remainder = MODULO(value, next)
         value = next
         next = remainder
      END DO
   END FUNCTION gcd
END PROGRAM topology_integer_unittest
