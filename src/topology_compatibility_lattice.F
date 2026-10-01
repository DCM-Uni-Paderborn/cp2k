!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Exact compatibility kernel and atomic quotient on a supplied reciprocal graph.
!> \author Thomas D. Kuehne, tkuehne@cp2k.org
! **************************************************************************************************
MODULE topology_compatibility_lattice
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE topology_atomic_signatures,      ONLY: atomic_signature_star_map,&
                                              atomic_signature_type
   USE topology_compatibility,          ONLY: little_group_compatibility,&
                                              little_group_edge_type
   USE topology_integer_lattice,        ONLY: integer_lattice_product,&
                                              integer_lattice_quotient,&
                                              integer_quotient_type
   USE topology_star_transport,         ONLY: star_operation

   IMPLICIT NONE
   PRIVATE

   TYPE, PUBLIC :: compatibility_lattice_type
      INTEGER(int_8), ALLOCATABLE :: constraints(:, :)
      INTEGER, ALLOCATABLE :: stars(:, :)
      TYPE(integer_quotient_type) :: quotient
   END TYPE compatibility_lattice_type
   PUBLIC :: compatibility_lattice_construct
CONTAINS

! **************************************************************************************************
!> \brief Assemble all edge restrictions in the catalogue's point-irrep/corepresentation basis.
!> \param catalogue Validated site-induced references and endpoint character tables
!> \param sources Directed source indices
!> \param destinations Directed destination indices
!> \param shifts Integer reciprocal target images, including nonzero self-edge winding
!> \param lattice Equations C*n=0, saturated integer kernel and quotient by atomic columns
!> \param tolerance Character and fractional geometry tolerance
!> \param status 0 success; -1 invalid input; -2 algebra/certificate failure; -3 unresolved resource limit
!> \param max_memory_mb Additional replicated storage estimate, excluding external-library workspace
!> \note Corepresentations are first restricted to unitary irreps; equality of unitary line counts
!>       is equivalent to equality of valid line corepresentation counts. A self-edge adds both
!>       endpoint restrictions into the same columns. Only duplicate/sign-equivalent equations
!>       are removed. Full-group star relations between supplied points are included; this does
!>       not establish complete reciprocal coverage or global band connectivity.
! **************************************************************************************************
   SUBROUTINE compatibility_lattice_construct(catalogue, sources, destinations, shifts, lattice, &
                                              tolerance, status, max_memory_mb)
      TYPE(atomic_signature_type), INTENT(IN)            :: catalogue
      INTEGER, INTENT(IN)                                :: sources(:), destinations(:), shifts(:, :)
      TYPE(compatibility_lattice_type), INTENT(OUT)      :: lattice
      REAL(dp), INTENT(IN)                               :: tolerance
      INTEGER, INTENT(OUT)                               :: status
      REAL(dp), INTENT(IN)                               :: max_memory_mb

      INTEGER                                            :: a, b, c, count, d, e, group_size, i, n, &
                                                            ne, np, ns, operation, p, q, shift(3)
      INTEGER(int_8)                                     :: max_rows
      INTEGER(int_8), ALLOCATABLE                        :: equations(:, :), left(:, :), &
                                                            right(:, :), row(:)
      INTEGER, ALLOCATABLE                               :: permutation(:, :), roots(:), stars(:, :)
      REAL(dp)                                           :: memory, residual
      TYPE(compatibility_lattice_type)                   :: work
      TYPE(little_group_edge_type)                       :: edge

      status = -1
      IF (.NOT. ieee_is_finite(max_memory_mb) .OR. .NOT. ieee_is_finite(tolerance)) RETURN
      IF (max_memory_mb <= 0.0_dp .OR. tolerance <= 0.0_dp .OR. tolerance >= 0.01_dp) RETURN
      IF (.NOT. ALLOCATED(catalogue%points) .OR. .NOT. ALLOCATED(catalogue%row_begin)) RETURN
      IF (.NOT. ALLOCATED(catalogue%matrix)) RETURN
      n = SIZE(catalogue%matrix, 1)
      np = SIZE(catalogue%points)
      ne = SIZE(sources)
      IF (n < 1 .OR. np < 1 .OR. SIZE(catalogue%row_begin) /= np + 1) RETURN
      IF (catalogue%row_begin(1) /= 1 .OR. catalogue%row_begin(np + 1) /= n + 1) RETURN
      DO p = 1, np
         IF (.NOT. ALLOCATED(catalogue%points(p)%dimensions)) RETURN
         IF (.NOT. ALLOCATED(catalogue%points(p)%little%product)) RETURN
         IF (catalogue%row_begin(p + 1) - catalogue%row_begin(p) /= SIZE(catalogue%points(p)%dimensions)) RETURN
      END DO
      IF (SIZE(destinations) /= ne .OR. ANY(SHAPE(shifts) /= [3, ne])) RETURN
      IF (ANY(sources < 1) .OR. ANY(sources > np)) RETURN
      IF (ANY(destinations < 1) .OR. ANY(destinations > np)) RETURN
      max_rows = 1_int_8 + n
      group_size = 0
      DO p = 1, np
         group_size = MAX(group_size, SIZE(catalogue%points(p)%little%product, 1))
      END DO
      DO e = 1, ne
         group_size = MAX(group_size, SIZE(catalogue%points(sources(e))%little%product, 1), &
                          SIZE(catalogue%points(destinations(e))%little%product, 1))
         max_rows = max_rows + MIN(SIZE(catalogue%points(sources(e))%little%product, 1), &
                                   SIZE(catalogue%points(destinations(e))%little%product, 1))
      END DO
      memory = (32.0_dp*REAL(max_rows, dp)*n + 512.0_dp*REAL(group_size, dp)**2)/1024.0_dp**2
      status = -3
      IF (max_rows > HUGE(0) .OR. memory > max_memory_mb) RETURN
      ALLOCATE (equations(INT(max_rows), n), row(n), roots(np), stars(6, np))
      equations = 0
      count = 0
      roots = 0
      ns = 0
      DO p = 1, np
         IF (roots(p) /= 0) CYCLE
         roots(p) = p
         DO q = p + 1, np
            IF (roots(q) /= 0) CYCLE
            CALL star_operation(catalogue%group, catalogue%points(p)%little%k, catalogue%points(q)%little%k, &
                                operation, shift, tolerance, status)
            IF (status == 1) CYCLE
            IF (status /= 0) RETURN
            CALL atomic_signature_star_map(catalogue, p, q, operation, permutation, tolerance, residual, status)
            IF (status /= 0) RETURN
            roots(q) = p
            ns = ns + 1
            stars(:, ns) = [p, q, operation, shift]
            a = catalogue%row_begin(p)
            b = catalogue%row_begin(p + 1) - 1
            c = catalogue%row_begin(q)
            DO i = 1, SIZE(permutation, 1)
               row = 0
               row(a:b) = INT(permutation(i, :), int_8)
               row(c + i - 1) = -1
               CALL append_equation(row, equations, count)
            END DO
         END DO
      END DO
      work%stars = stars(:, :ns)
      DEALLOCATE (roots, stars)
      DO e = 1, ne
         p = sources(e)
         q = destinations(e)
         IF (MAXVAL(ABS(catalogue%points(q)%little%k + REAL(shifts(:, e), dp) - &
                        catalogue%points(p)%little%k)) <= tolerance) CYCLE
         CALL little_group_compatibility(catalogue%points(p)%little, catalogue%points(q)%little, &
                                         shifts(:, e), edge, tolerance, residual, status)
         IF (status /= 0) RETURN
         IF (ALLOCATED(catalogue%points(p)%coreps%restriction)) THEN
            CALL integer_lattice_product(INT(edge%left_branch, int_8), &
                                         INT(catalogue%points(p)%coreps%restriction, int_8), left, status)
            IF (status /= 0) RETURN
         ELSE
            left = INT(edge%left_branch, int_8)
         END IF
         IF (ALLOCATED(catalogue%points(q)%coreps%restriction)) THEN
            CALL integer_lattice_product(INT(edge%right_branch, int_8), &
                                         INT(catalogue%points(q)%coreps%restriction, int_8), right, status)
            IF (status /= 0) RETURN
         ELSE
            right = INT(edge%right_branch, int_8)
         END IF
         a = catalogue%row_begin(p)
         b = catalogue%row_begin(p + 1) - 1
         c = catalogue%row_begin(q)
         d = catalogue%row_begin(q + 1) - 1
         status = -2
         IF (SIZE(left, 2) /= b - a + 1 .OR. SIZE(right, 2) /= d - c + 1) RETURN
         IF (SIZE(left, 1) /= SIZE(right, 1)) RETURN
         DO i = 1, SIZE(left, 1)
            row = 0
            row(a:b) = left(i, :)
            row(c:d) = row(c:d) - right(i, :)
            CALL append_equation(row, equations, count)
         END DO
      END DO
      work%constraints = equations(:MAX(1, count), :)
      DEALLOCATE (equations, row)
      ! Bound dense Smith certificates using the retained, deduplicated equations.
      memory = 96.0_dp*(REAL(MAX(1, count), dp)**2 + REAL(n, dp)**2 + &
                        3.0_dp*MAX(1, count)*n + REAL(SIZE(catalogue%matrix, 2), dp)**2)/1024.0_dp**2 + &
               512.0_dp*REAL(group_size, dp)**2/1024.0_dp**2
      status = -3
      IF (memory > max_memory_mb) RETURN
      CALL integer_lattice_quotient(work%constraints, INT(catalogue%matrix, int_8), work%quotient, status)
      IF (status /= 0) RETURN
      lattice = work
   END SUBROUTINE compatibility_lattice_construct

! **************************************************************************************************
!> \brief Retain a nonzero homogeneous equation once, normalized by gcd and leading sign.
!> \param row Candidate equation, normalized in place
!> \param equations Preallocated upper bound on rows
!> \param count Number of retained equations
! **************************************************************************************************
   SUBROUTINE append_equation(row, equations, count)
      INTEGER(int_8), INTENT(INOUT)                      :: row(:), equations(:, :)
      INTEGER, INTENT(INOUT)                             :: count

      INTEGER                                            :: first, j
      INTEGER(int_8)                                     :: divisor

      first = FINDLOC(row /= 0, .TRUE., DIM=1)
      IF (first == 0) RETURN
      divisor = 0
      DO j = 1, SIZE(row)
         divisor = gcd(divisor, ABS(row(j)))
      END DO
      IF (row(first) < 0) divisor = -divisor
      row = row/divisor
      DO j = 1, count
         IF (ALL(equations(j, :) == row)) RETURN
      END DO
      count = count + 1
      equations(count, :) = row
   END SUBROUTINE append_equation

! **************************************************************************************************
!> \brief Nonnegative gcd for canonical homogeneous equations.
!> \param a Nonnegative integer
!> \param b Nonnegative integer
!> \return Greatest common divisor
! **************************************************************************************************
   INTEGER(int_8) FUNCTION gcd(a, b) RESULT(value)
      INTEGER(int_8), INTENT(IN)                         :: a, b

      INTEGER(int_8)                                     :: other, remainder

      value = a
      other = b
      DO WHILE (other /= 0)
         remainder = MODULO(value, other)
         value = other
         other = remainder
      END DO
   END FUNCTION gcd
END MODULE topology_compatibility_lattice
