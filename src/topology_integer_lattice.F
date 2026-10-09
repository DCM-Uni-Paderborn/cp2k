!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Checked integer-lattice membership and nonnegative atomic-signature decomposition.
!> \author Thomas D. Kuehne, tkuehne@cp2k.org
! **************************************************************************************************
MODULE topology_integer_lattice
   USE kinds,                           ONLY: int_8

   IMPLICIT NONE
   PRIVATE

   INTEGER(int_8), PARAMETER :: integer_limit = HUGE(0_int_8)/16
   TYPE, PUBLIC :: integer_lattice_type
      PRIVATE
      INTEGER :: rank = 0
      INTEGER(int_8), ALLOCATABLE :: matrix(:, :), left(:, :), right(:, :), diagonal(:)
   END TYPE integer_lattice_type

   TYPE, PUBLIC :: integer_quotient_type
      INTEGER(int_8), ALLOCATABLE :: kernel(:, :), coordinates(:, :), left(:, :), invariants(:)
      TYPE(integer_lattice_type), PRIVATE :: kernel_factor
   END TYPE integer_quotient_type

   PUBLIC :: integer_lattice_factor, integer_lattice_solve, integer_lattice_certificate
   PUBLIC :: nonnegative_integer_solve
   PUBLIC :: integer_lattice_quotient, integer_quotient_class, integer_lattice_product
CONTAINS

! **************************************************************************************************
!> \brief Checked integer product for lattice certificates, including empty dimensions.
!> \param a Left operand
!> \param b Right operand
!> \param c Product, unallocated on failure
!> \param status 0 success; -1 invalid shape; -3 arithmetic/work bound
! **************************************************************************************************
   SUBROUTINE integer_lattice_product(a, b, c, status)
      INTEGER(int_8), INTENT(IN)                         :: a(:, :), b(:, :)
      INTEGER(int_8), ALLOCATABLE, INTENT(OUT)           :: c(:, :)
      INTEGER, INTENT(OUT)                               :: status

      INTEGER(int_8), ALLOCATABLE                        :: work(:, :)

      status = -1
      IF (SIZE(a, 2) /= SIZE(b, 1)) RETURN
      status = -3
      IF (INT(SIZE(a, 1), int_8)*SIZE(b, 2) > 10000000_int_8) RETURN
      IF (ANY(a > integer_limit) .OR. ANY(a < -integer_limit)) RETURN
      IF (ANY(b > integer_limit) .OR. ANY(b < -integer_limit)) RETURN
      CALL checked_product(a, b, work, status)
      IF (status == 0) CALL MOVE_ALLOC(work, c)
   END SUBROUTINE integer_lattice_product

! **************************************************************************************************
!> \brief Exact quotient ker_Z(C)/im_Z(A), rejecting generators that violate C*A=0.
!> \param constraints Homogeneous compatibility equations C, at least one row/column
!> \param generators Atomic/reference columns A in the same integer coordinate basis
!> \param quotient Saturated kernel K, A=K*coordinates, Smith row map and quotient invariants
!> \param status 0 success; -1 invalid data; -2 inconsistent generators/certificate; -3 unresolved bound
!> \note Invariants 0 denote free Z factors, 1 trivial factors, and d>1 cyclic Z_d factors.
!>       This computes a supplied integer quotient, not physical completeness of either lattice.
! **************************************************************************************************
   SUBROUTINE integer_lattice_quotient(constraints, generators, quotient, status)
      INTEGER(int_8), INTENT(IN)                         :: constraints(:, :), generators(:, :)
      TYPE(integer_quotient_type), INTENT(OUT)           :: quotient
      INTEGER, INTENT(OUT)                               :: status

      INTEGER                                            :: j, n, q, rank
      INTEGER(int_8), ALLOCATABLE                        :: diagonal(:), left(:, :), right(:, :), &
                                                            values(:), verified(:, :)
      LOGICAL                                            :: exists
      TYPE(integer_lattice_type)                         :: factor
      TYPE(integer_quotient_type)                        :: work

      status = -1
      n = SIZE(constraints, 2)
      IF (n < 1 .OR. SIZE(constraints, 1) < 1) RETURN
      IF (SIZE(generators, 1) /= n .OR. SIZE(generators, 2) < 1) RETURN
      CALL integer_lattice_factor(constraints, factor, status)
      IF (status /= 0) RETURN
      CALL integer_lattice_certificate(factor, left, right, diagonal, rank)
      q = n - rank
      work%kernel = right(:, rank + 1:)
      CALL integer_lattice_product(constraints, work%kernel, verified, status)
      IF (status /= 0) RETURN
      status = -2
      IF (ANY(verified /= 0)) RETURN
      ALLOCATE (work%coordinates(q, SIZE(generators, 2)), work%invariants(q))
      IF (q == 0) THEN
         IF (ANY(generators /= 0)) RETURN
         ALLOCATE (work%left(0, 0))
      ELSE
         CALL integer_lattice_factor(work%kernel, work%kernel_factor, status)
         IF (status /= 0) RETURN
         DO j = 1, SIZE(generators, 2)
            CALL integer_lattice_solve(work%kernel_factor, generators(:, j), exists, values, status)
            IF (status /= 0) RETURN
            status = -2
            IF (.NOT. exists) RETURN
            work%coordinates(:, j) = values
         END DO
         CALL integer_lattice_factor(work%coordinates, factor, status)
         IF (status /= 0) RETURN
         CALL integer_lattice_certificate(factor, work%left, right, diagonal, rank)
         work%invariants = 0
         work%invariants(:SIZE(diagonal)) = diagonal
      END IF
      CALL integer_lattice_product(work%kernel, work%coordinates, verified, status)
      IF (status /= 0) RETURN
      status = -2
      IF (ANY(verified /= generators)) RETURN
      quotient = work
      status = 0
   END SUBROUTINE integer_lattice_quotient

! **************************************************************************************************
!> \brief Class of a target in a previously certified quotient, with compatibility separate from class.
!> \param quotient Valid output of integer_lattice_quotient
!> \param TARGET Signature in the original coordinate basis
!> \param compatible True if the target lies in the exact integer compatibility kernel
!> \param indices Residues modulo nonzero invariants, or signed free coordinates; empty on failure
!> \param status 0 decided; -1 invalid data; -3 unresolved arithmetic bound
!> \note All indices zero is equivalent to signed generator membership, not nonnegative membership.
! **************************************************************************************************
   SUBROUTINE integer_quotient_class(quotient, TARGET, compatible, indices, status)
      TYPE(integer_quotient_type), INTENT(IN)            :: quotient
      INTEGER(int_8), INTENT(IN)                         :: target(:)
      LOGICAL, INTENT(OUT)                               :: compatible
      INTEGER(int_8), ALLOCATABLE, INTENT(OUT)           :: indices(:)
      INTEGER, INTENT(OUT)                               :: status

      INTEGER                                            :: j, q
      INTEGER(int_8), ALLOCATABLE                        :: transformed(:, :), values(:)

      compatible = .FALSE.
      status = -1
      IF (.NOT. ALLOCATED(quotient%kernel) .OR. .NOT. ALLOCATED(quotient%invariants)) RETURN
      IF (SIZE(TARGET) /= SIZE(quotient%kernel, 1)) RETURN
      q = SIZE(quotient%kernel, 2)
      IF (q == 0) THEN
         status = 0
         compatible = ALL(TARGET == 0)
         IF (compatible) ALLOCATE (indices(0))
         RETURN
      END IF
      CALL integer_lattice_solve(quotient%kernel_factor, TARGET, compatible, values, status)
      IF (status /= 0 .OR. .NOT. compatible) RETURN
      CALL integer_lattice_product(quotient%left, RESHAPE(values, [q, 1]), transformed, status)
      IF (status /= 0) THEN
         compatible = .FALSE.
         RETURN
      END IF
      DO j = 1, q
         IF (quotient%invariants(j) > 0) transformed(j, 1) = MODULO(transformed(j, 1), quotient%invariants(j))
      END DO
      indices = transformed(:, 1)
   END SUBROUTINE integer_quotient_class

! **************************************************************************************************
!> \brief Smith form U*A*V=D, with unimodular U,V, positive diagonal and divisibility order.
!> \param matrix Integer matrix, rectangular and rank-deficient matrices allowed
!> \param factor Original matrix and verified exact certificate; empty on failure
!> \param status 0 success; -1 invalid dimensions; -2 certificate failure; -3 arithmetic/work limit
!> \note Euclidean row/column steps use no floating point. The complete certificate is multiplied
!>       out with checked arithmetic before acceptance. Coefficient growth can return unresolved.
! **************************************************************************************************
   SUBROUTINE integer_lattice_factor(matrix, factor, status)
      INTEGER(int_8), INTENT(IN)                         :: matrix(:, :)
      TYPE(integer_lattice_type), INTENT(OUT)            :: factor
      INTEGER, INTENT(OUT)                               :: status

      INTEGER                                            :: i, j, k, m, n, p, s
      INTEGER(int_8)                                     :: q, smallest, steps
      INTEGER(int_8), ALLOCATABLE                        :: a(:, :), test(:, :), verified(:, :)
      LOGICAL                                            :: restart
      TYPE(integer_lattice_type)                         :: work

      status = -1
      m = SIZE(matrix, 1)
      n = SIZE(matrix, 2)
      IF (MIN(m, n) < 1) RETURN
      status = -3
      IF (INT(m, int_8)**2 + INT(n, int_8)**2 + INT(m, int_8)*n > 10000000_int_8) RETURN
      IF (ANY(matrix > integer_limit) .OR. ANY(matrix < -integer_limit)) RETURN
      work%matrix = matrix
      a = matrix
      ALLOCATE (work%left(m, m), work%right(n, n), work%diagonal(MIN(m, n)))
      work%left = 0
      work%right = 0
      work%diagonal = 0
      DO i = 1, m
         work%left(i, i) = 1
      END DO
      DO i = 1, n
         work%right(i, i) = 1
      END DO
      steps = 0
      DO k = 1, MIN(m, n)
         smallest = integer_limit + 1
         p = 0
         s = 0
         DO j = k, n
            DO i = k, m
               IF (a(i, j) == 0 .OR. ABS(a(i, j)) >= smallest) CYCLE
               smallest = ABS(a(i, j))
               p = i
               s = j
            END DO
         END DO
         IF (p == 0) EXIT
         CALL exchange_rows(a, k, p)
         CALL exchange_rows(work%left, k, p)
         CALL exchange_columns(a, k, s)
         CALL exchange_columns(work%right, k, s)
         DO
            steps = steps + 1
            status = -3
            IF (steps > 1000000_int_8) RETURN
            restart = .FALSE.
            DO i = k + 1, m
               IF (a(i, k) == 0) CYCLE
               q = a(i, k)/a(k, k)
               CALL subtract_scaled(a(i, :), q, a(k, :), status)
               IF (status /= 0) RETURN
               CALL subtract_scaled(work%left(i, :), q, work%left(k, :), status)
               IF (status /= 0) RETURN
               IF (a(i, k) /= 0) THEN
                  CALL exchange_rows(a, k, i)
                  CALL exchange_rows(work%left, k, i)
                  restart = .TRUE.
                  EXIT
               END IF
            END DO
            IF (restart) CYCLE
            DO j = k + 1, n
               IF (a(k, j) == 0) CYCLE
               q = a(k, j)/a(k, k)
               CALL subtract_scaled(a(:, j), q, a(:, k), status)
               IF (status /= 0) RETURN
               CALL subtract_scaled(work%right(:, j), q, work%right(:, k), status)
               IF (status /= 0) RETURN
               IF (a(k, j) /= 0) THEN
                  CALL exchange_columns(a, k, j)
                  CALL exchange_columns(work%right, k, j)
                  restart = .TRUE.
                  EXIT
               END IF
            END DO
            IF (restart) CYCLE
            ! A pivot must divide the entire trailing block, not only its own row and column.
            p = 0
            DO j = k + 1, n
               DO i = k + 1, m
                  IF (MODULO(a(i, j), ABS(a(k, k))) == 0) CYCLE
                  p = i
                  EXIT
               END DO
               IF (p /= 0) EXIT
            END DO
            IF (p == 0) EXIT
            CALL subtract_scaled(a(k, :), -1_int_8, a(p, :), status)
            IF (status /= 0) RETURN
            CALL subtract_scaled(work%left(k, :), -1_int_8, work%left(p, :), status)
            IF (status /= 0) RETURN
         END DO
         IF (a(k, k) < 0) THEN
            a(k, :) = -a(k, :)
            work%left(k, :) = -work%left(k, :)
         END IF
         work%diagonal(k) = a(k, k)
         work%rank = k
      END DO
      CALL checked_product(work%left, matrix, test, status)
      IF (status /= 0) RETURN
      CALL checked_product(test, work%right, verified, status)
      IF (status /= 0) RETURN
      status = -2
      IF (ANY(verified /= a)) RETURN
      DO k = 1, MIN(m, n)
         a(k, k) = a(k, k) - work%diagonal(k)
      END DO
      IF (ANY(a /= 0)) RETURN
      DO k = 2, work%rank
         IF (MODULO(work%diagonal(k), work%diagonal(k - 1)) /= 0) RETURN
      END DO
      factor = work
      status = 0
   END SUBROUTINE integer_lattice_factor

! **************************************************************************************************
!> \brief Exact signed integer solution or exact nonmembership from a previously factored matrix.
!> \param factor Validated Smith form
!> \param TARGET Integer target vector
!> \param exists True only for a verified solution when status is zero
!> \param coefficients Signed witness, allocated only on success with exists true
!> \param status 0 decided; -1 invalid data; -2 witness failure; -3 arithmetic limit, not nonmembership
! **************************************************************************************************
   SUBROUTINE integer_lattice_solve(factor, TARGET, exists, coefficients, status)
      TYPE(integer_lattice_type), INTENT(IN)             :: factor
      INTEGER(int_8), INTENT(IN)                         :: target(:)
      LOGICAL, INTENT(OUT)                               :: exists
      INTEGER(int_8), ALLOCATABLE, INTENT(OUT)           :: coefficients(:)
      INTEGER, INTENT(OUT)                               :: status

      INTEGER                                            :: i, m, n
      INTEGER(int_8), ALLOCATABLE                        :: basis(:, :), solution(:, :), &
                                                            transformed(:, :), verified(:, :)

      exists = .FALSE.
      status = -1
      IF (.NOT. ALLOCATED(factor%matrix)) RETURN
      m = SIZE(factor%matrix, 1)
      n = SIZE(factor%matrix, 2)
      IF (SIZE(TARGET) /= m) RETURN
      status = -3
      IF (ANY(TARGET > integer_limit) .OR. ANY(TARGET < -integer_limit)) RETURN
      CALL checked_product(factor%left, RESHAPE(TARGET, [m, 1]), transformed, status)
      IF (status /= 0) RETURN
      IF (ANY(transformed(factor%rank + 1:, 1) /= 0)) RETURN
      ALLOCATE (basis(n, 1))
      basis = 0
      DO i = 1, factor%rank
         IF (MODULO(transformed(i, 1), factor%diagonal(i)) /= 0) RETURN
         basis(i, 1) = transformed(i, 1)/factor%diagonal(i)
      END DO
      CALL checked_product(factor%right, basis, solution, status)
      IF (status /= 0) RETURN
      CALL checked_product(factor%matrix, solution, verified, status)
      IF (status /= 0) RETURN
      status = -2
      IF (ANY(verified(:, 1) /= TARGET)) RETURN
      coefficients = solution(:, 1)
      exists = .TRUE.
      status = 0
   END SUBROUTINE integer_lattice_solve

! **************************************************************************************************
!> \brief Expose the integer certificate, without mutable access to the cached factorization.
!> \param factor Validated factorization
!> \param left Unimodular row transformation
!> \param right Unimodular column transformation
!> \param diagonal Smith diagonal of length min(rows,columns), including trailing zeros
!> \param rank Number of positive entries; -1 for an empty factorization
! **************************************************************************************************
   SUBROUTINE integer_lattice_certificate(factor, left, right, diagonal, rank)
      TYPE(integer_lattice_type), INTENT(IN)             :: factor
      INTEGER(int_8), ALLOCATABLE, INTENT(OUT)           :: left(:, :), right(:, :), diagonal(:)
      INTEGER, INTENT(OUT)                               :: rank

      rank = -1
      IF (.NOT. ALLOCATED(factor%matrix)) RETURN
      left = factor%left
      right = factor%right
      diagonal = factor%diagonal
      rank = factor%rank
   END SUBROUTINE integer_lattice_certificate

! **************************************************************************************************
!> \brief Exact finite search for nonnegative integer coefficients of a nonnegative signature matrix.
!> \param matrix Nonnegative atomic signature columns; zero columns are assigned coefficient zero
!> \param TARGET Nonnegative target counts
!> \param exists A verified witness exists, valid as a negative decision only when status is zero
!> \param coefficients Nonnegative witness, allocated only for a verified solution
!> \param visited Number of visited search nodes
!> \param status 0 decided; 1 unresolved search limit; -1 invalid; -2 witness failure; -3 arithmetic/work limit
!> \param max_nodes Required search budget; exhaustion must not be interpreted as an obstruction
!> \note Coordinate bounds and suffix row gcds prune the complete coefficient search. There is no
!>       floating-point linear program, rounding or assumption of a unique decomposition.
! **************************************************************************************************
   SUBROUTINE nonnegative_integer_solve(matrix, TARGET, exists, coefficients, visited, status, max_nodes)
      INTEGER(int_8), INTENT(IN) :: matrix(:, :), TARGET(:), max_nodes
      LOGICAL, INTENT(OUT) :: exists
      INTEGER(int_8), ALLOCATABLE, INTENT(OUT) :: coefficients(:)
      INTEGER(int_8), INTENT(OUT) :: visited
      INTEGER, INTENT(OUT) :: status
      INTEGER(int_8), ALLOCATABLE :: suffix(:, :), trial(:), verified(:, :)
      INTEGER :: m, n, i, j

      exists = .FALSE.
      visited = 0
      status = -1
      m = SIZE(matrix, 1)
      n = SIZE(matrix, 2)
      IF (MIN(m, n) < 1 .OR. SIZE(TARGET) /= m .OR. max_nodes < 0) RETURN
      IF (ANY(matrix < 0) .OR. ANY(TARGET < 0)) RETURN
      status = -3
      IF (ANY(matrix > integer_limit) .OR. ANY(TARGET > integer_limit)) RETURN
      IF (n > 4096 .OR. INT(m, int_8)*(n + 1_int_8) > 10000000_int_8) RETURN
      ALLOCATE (suffix(m, n + 1), trial(n))
      suffix(:, n + 1) = 0
      trial = 0
      DO j = n, 1, -1
         DO i = 1, m
            suffix(i, j) = integer_gcd(suffix(i, j + 1), matrix(i, j))
         END DO
      END DO
      status = 0
      CALL search(1, TARGET)
      IF (.NOT. exists .OR. status /= 0) RETURN
      CALL checked_product(matrix, RESHAPE(trial, [n, 1]), verified, status)
      IF (status /= 0) THEN
         exists = .FALSE.
         RETURN
      END IF
      status = -2
      exists = .FALSE.
      IF (ANY(verified(:, 1) /= TARGET)) RETURN
      coefficients = trial
      exists = .TRUE.
      status = 0
   CONTAINS
! **************************************************************************************************
!> \brief ...
!> \param column ...
!> \param remainder ...
! **************************************************************************************************
      RECURSIVE SUBROUTINE search(column, remainder)
      INTEGER, INTENT(IN)                                :: column
      INTEGER(int_8), INTENT(IN)                         :: remainder(:)

      INTEGER                                            :: row
      INTEGER(int_8)                                     :: bound, q, reduced(SIZE(remainder))

         IF (visited >= max_nodes) THEN
            status = 1
            RETURN
         END IF
         visited = visited + 1
         IF (ALL(remainder == 0)) THEN
            trial(column:) = 0
            exists = .TRUE.
            RETURN
         END IF
         IF (column > n) RETURN
         DO row = 1, m
            IF (suffix(row, column) == 0) THEN
               IF (remainder(row) /= 0) RETURN
            ELSE
               IF (MODULO(remainder(row), suffix(row, column)) /= 0) RETURN
            END IF
         END DO
         bound = integer_limit
         DO row = 1, m
            IF (matrix(row, column) == 0) CYCLE
            bound = MIN(bound, remainder(row)/matrix(row, column))
         END DO
         IF (ALL(matrix(:, column) == 0)) bound = 0
         DO q = bound, 0, -1
            trial(column) = q
            ! q*a never exceeds remainder componentwise, so this multiplication cannot overflow.
            reduced = remainder - q*matrix(:, column)
            CALL search(column + 1, reduced)
            IF (exists .OR. status /= 0) RETURN
         END DO
         trial(column) = 0
      END SUBROUTINE search
   END SUBROUTINE nonnegative_integer_solve

! **************************************************************************************************
!> \brief Integer x=x-q*y with conservative bounds checked before multiplication.
!> \param x Accumulator, potentially partially updated on failure
!> \param q Signed scalar
!> \param y Source vector
!> \param status Zero on success; -3 if a bound would be crossed
! **************************************************************************************************
   SUBROUTINE subtract_scaled(x, q, y, status)
      INTEGER(int_8), INTENT(INOUT) :: x(:)
      INTEGER(int_8), INTENT(IN) :: q, y(:)
      INTEGER, INTENT(OUT) :: status
      INTEGER(int_8) :: value
      INTEGER :: i
      status = -3
      DO i = 1, SIZE(x)
         IF (y(i) /= 0) THEN
            IF (ABS(q) > integer_limit/ABS(y(i))) RETURN
         END IF
         value = x(i) - q*y(i)
         IF (ABS(value) > integer_limit) RETURN
         x(i) = value
      END DO
      status = 0
   END SUBROUTINE subtract_scaled

! **************************************************************************************************
!> \brief Checked exact matrix product, including intermediate accumulation bounds.
!> \param a Left operand
!> \param b Right operand
!> \param c Product
!> \param status Zero on success; -3 on coefficient growth
! **************************************************************************************************
   SUBROUTINE checked_product(a, b, c, status)
      INTEGER(int_8), INTENT(IN)                         :: a(:, :), b(:, :)
      INTEGER(int_8), ALLOCATABLE, INTENT(OUT)           :: c(:, :)
      INTEGER, INTENT(OUT)                               :: status

      INTEGER                                            :: j, k

      ALLOCATE (c(SIZE(a, 1), SIZE(b, 2)))
      c = 0
      DO j = 1, SIZE(b, 2)
         DO k = 1, SIZE(a, 2)
            CALL subtract_scaled(c(:, j), -b(k, j), a(:, k), status)
            IF (status /= 0) RETURN
         END DO
      END DO
      status = 0
   END SUBROUTINE checked_product

! **************************************************************************************************
!> \brief Swap two integer matrix rows.
!> \param a Matrix
!> \param i First row
!> \param j Second row
! **************************************************************************************************
   SUBROUTINE exchange_rows(a, i, j)
      INTEGER(int_8), INTENT(INOUT)                      :: a(:, :)
      INTEGER, INTENT(IN)                                :: i, j

      INTEGER(int_8)                                     :: tmp(SIZE(a, 2))

      tmp = a(i, :)
      a(i, :) = a(j, :)
      a(j, :) = tmp
   END SUBROUTINE exchange_rows

! **************************************************************************************************
!> \brief Swap two integer matrix columns.
!> \param a Matrix
!> \param i First column
!> \param j Second column
! **************************************************************************************************
   SUBROUTINE exchange_columns(a, i, j)
      INTEGER(int_8), INTENT(INOUT)                      :: a(:, :)
      INTEGER, INTENT(IN)                                :: i, j

      INTEGER(int_8)                                     :: tmp(SIZE(a, 1))

      tmp = a(:, i)
      a(:, i) = a(:, j)
      a(:, j) = tmp
   END SUBROUTINE exchange_columns

! **************************************************************************************************
!> \brief Nonnegative greatest common divisor.
!> \param a First nonnegative integer
!> \param b Second nonnegative integer
!> \return gcd(a,b)
! **************************************************************************************************
   INTEGER(int_8) FUNCTION integer_gcd(a, b) RESULT(divisor)
      INTEGER(int_8), INTENT(IN)                         :: a, b

      INTEGER(int_8)                                     :: next, remainder

      divisor = a
      next = b
      DO WHILE (next /= 0)
         remainder = MODULO(divisor, next)
         divisor = next
         next = remainder
      END DO
   END FUNCTION integer_gcd
END MODULE topology_integer_lattice
