!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Generate connected site-symmetry families from periodic affine fixed-point equations.
!> \author Thomas D. Kuehne, tkuehne@cp2k.org
! **************************************************************************************************
MODULE topology_wyckoff
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE topology_little_group,           ONLY: little_group_factor_system

   IMPLICIT NONE
   PRIVATE

   !> B*x = offset modulo integers describes the closure of a site-symmetry stratum.
   TYPE, PUBLIC :: wyckoff_family_type
      INTEGER :: dimension = 3, multiplicity = 0
      INTEGER :: normal(3, 3) = 0
      REAL(dp) :: offset(3) = 0.0_dp, position(3) = 0.0_dp
      LOGICAL :: maximal_site_symmetry = .FALSE.
      INTEGER, ALLOCATABLE :: stabilizer(:), specializations(:)
   END TYPE wyckoff_family_type

   TYPE :: affine_torus_type
      INTEGER :: rank = 0, normal(3, 3) = 0
      REAL(dp) :: offset(3) = 0.0_dp
   END TYPE affine_torus_type

   PUBLIC :: wyckoff_families, wyckoff_match_affine
CONTAINS

! **************************************************************************************************
!> \brief Match a complete affine reference catalogue to generated families without choosing seeds.
!> \param families Complete output of wyckoff_families in the supplied group convention
!> \param rotations Complete fractional space-group rotations
!> \param translations Corresponding fractional translations
!> \param positions Reference affine origins; these need not be generic site positions
!> \param directions Integer columns spanning each reference family's free real parameters
!> \param multiplicities Reference orbit multiplicities in this cell, including centering
!> \param mapping Generated family index for each reference entry, unallocated on failure
!> \param operations Group operation mapping each reference closure to its generated family
!> \param tolerance Absolute fractional geometry tolerance
!> \param status Zero on success, -1 invalid input, -2 missing/ambiguous/nonbijective match
!> \note The caller must first put the reference in the same cell and origin convention.
!>       Rank, whole tangent space and affine offset are checked, not a sampled point.
!>       Names are external metadata: no standard letters or EBR labels are inferred here.
! **************************************************************************************************
   SUBROUTINE wyckoff_match_affine(families, rotations, translations, positions, directions, &
                                   multiplicities, mapping, operations, tolerance, status)
      TYPE(wyckoff_family_type), INTENT(IN)              :: families(:)
      INTEGER, INTENT(IN)                                :: rotations(:, :, :)
      REAL(dp), INTENT(IN)                               :: translations(:, :), positions(:, :)
      INTEGER, INTENT(IN)                                :: directions(:, :, :), multiplicities(:)
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: mapping(:), operations(:)
      REAL(dp), INTENT(IN)                               :: tolerance
      INTEGER, INTENT(OUT)                               :: status

      INTEGER                                            :: basis(3, 3), g, i, j, n, rank, rows(3), &
                                                            stat
      INTEGER(int_8)                                     :: tangent(3, 3)
      INTEGER, ALLOCATABLE                               :: indices(:), witnesses(:)
      REAL(dp)                                           :: delta(3), image(3)

      status = -1
      n = SIZE(families)
      IF (n < 1 .OR. SIZE(multiplicities) /= n) RETURN
      IF (ANY(SHAPE(positions) /= [3, n]) .OR. ANY(SHAPE(directions) /= [3, 3, n])) RETURN
      IF (SIZE(rotations, 1) /= 3 .OR. SIZE(rotations, 2) /= 3 .OR. SIZE(rotations, 3) < 1) RETURN
      IF (ANY(SHAPE(translations) /= [3, SIZE(rotations, 3)])) RETURN
      IF (.NOT. ieee_is_finite(tolerance)) RETURN
      IF (tolerance <= 0.0_dp .OR. tolerance >= 0.1_dp) RETURN
      IF (.NOT. ALL(ieee_is_finite(positions)) .OR. .NOT. ALL(ieee_is_finite(translations))) RETURN
      IF (MAXVAL(ABS(positions)) > 1.e6_dp .OR. MAXVAL(ABS(translations)) > 1.e6_dp) RETURN
      IF (MAXVAL(ABS(INT(rotations, int_8))) > 100000_int_8) RETURN
      IF (ANY(multiplicities < 1)) RETURN
      DO i = 1, n
         IF (families(i)%dimension < 0 .OR. families(i)%dimension > 3) RETURN
         IF (families(i)%multiplicity < 1 .OR. .NOT. ALL(ieee_is_finite(families(i)%offset))) RETURN
         IF (MAXVAL(ABS(families(i)%offset)) > 1.e6_dp) RETURN
         CALL saturated_normals(families(i)%normal, rows, rank, basis, stat)
         IF (stat /= 0 .OR. rank /= 3 - families(i)%dimension) RETURN
      END DO
      ALLOCATE (indices(n), witnesses(n))
      indices = 0
      witnesses = 0
      DO j = 1, n
         status = -1
         CALL saturated_normals(TRANSPOSE(directions(:, :, j)), rows, rank, basis, stat)
         IF (stat /= 0) RETURN
         status = -2
         DO i = 1, n
            IF (families(i)%dimension /= rank .OR. families(i)%multiplicity /= multiplicities(j)) CYCLE
            DO g = 1, SIZE(rotations, 3)
               tangent = MATMUL(INT(rotations(:, :, g), int_8), INT(directions(:, :, j), int_8))
               IF (ANY(MATMUL(INT(families(i)%normal, int_8), tangent) /= 0_int_8)) CYCLE
               image = MATMUL(REAL(rotations(:, :, g), dp), positions(:, j)) + translations(:, g)
               delta = MATMUL(REAL(families(i)%normal, dp), image) - families(i)%offset
               IF (MAXVAL(ABS(delta - ANINT(delta))) > tolerance) CYCLE
               IF (indices(j) /= 0) RETURN
               indices(j) = i
               witnesses(j) = g
               EXIT
            END DO
         END DO
         IF (indices(j) == 0) RETURN
      END DO
      DO i = 1, n
         IF (COUNT(indices == i) /= 1) RETURN
      END DO
      CALL MOVE_ALLOC(indices, mapping)
      CALL MOVE_ALLOC(witnesses, operations)
      status = 0
   END SUBROUTINE wyckoff_match_affine

! **************************************************************************************************
!> \brief Enumerate all connected fixed strata, modulo the supplied space-group action.
!> \param cell Cartesian lattice vectors as columns
!> \param rotations Complete integer fractional rotation group modulo cell translations
!> \param translations Fractional coset translations in the supplied cell and origin convention
!> \param antiunitary Antiunitary flags; these affect site stabilizers but not spatial fixed sets
!> \param families Generic seeds, dimensions, multiplicities, stabilizers and specialization indices
!> \param tolerance Absolute fractional geometry tolerance
!> \param status 0 success; -1/-2 invalid group; -3 arithmetic/work limit; -4 unstable generic seed
!> \note No coordinate sampling or Wyckoff database is used. All nonempty intersections of
!>       (W-I)*x=-tau modulo integers are retained before equivalent families are identified.
!>       Higher-symmetry subsets are excluded from each generic seed, not from its affine closure.
!>       Output has no standard letters and is not an elementary-band-representation catalogue.
!>       A work/arithmetic failure returns no partial list. Centering cosets are retained.
!>       Specializations list all higher-symmetry families in a closure, not just immediate covers.
!>       No specializations means a maximal site group, not necessarily an elementary induced band.
! **************************************************************************************************
   SUBROUTINE wyckoff_families(cell, rotations, translations, antiunitary, families, tolerance, status)
      REAL(dp), INTENT(IN)                               :: cell(3, 3)
      INTEGER, INTENT(IN)                                :: rotations(:, :, :)
      REAL(dp), INTENT(IN)                               :: translations(:, :)
      LOGICAL, INTENT(IN)                                :: antiunitary(:)
      TYPE(wyckoff_family_type), ALLOCATABLE, &
         INTENT(OUT)                                     :: families(:)
      REAL(dp), INTENT(IN)                               :: tolerance
      INTEGER, INTENT(OUT)                               :: status

      COMPLEX(dp), ALLOCATABLE                           :: factors(:, :)
      INTEGER                                            :: a(6, 3), basis(3, 3), g, i, j, nc, nf, &
                                                            ng, nr, old, r, rows(3), trial
      INTEGER, ALLOCATABLE                               :: mapping(:), products(:, :)
      LOGICAL                                            :: equivalent
      LOGICAL, ALLOCATABLE                               :: expected(:), fixed(:), specializes(:)
      REAL(dp) :: b(6), delta(3), direction(3, 3), image(3), mapped(3, 3), nullspace(3, 3), &
         point(3), right(3, 3), seed(3), tau(3), variation(3)
      TYPE(affine_torus_type), ALLOCATABLE               :: components(:)
      TYPE(wyckoff_family_type), ALLOCATABLE             :: work(:)

      CALL little_group_factor_system(cell, rotations, translations, antiunitary, [0.0_dp, 0.0_dp, 0.0_dp], &
                                      .FALSE., mapping, products, factors, tolerance, status)
      IF (status /= 0) RETURN
      ng = SIZE(mapping)
      ALLOCATE (components(32))
      nc = 1
      DO g = 1, ng
         a = 0
         a(4:6, :) = rotations(:, :, mapping(g))
         DO i = 1, 3
            a(i + 3, i) = a(i + 3, i) - 1
         END DO
         IF (ALL(a == 0)) CYCLE
         tau = MODULO(translations(:, mapping(g)), 1.0_dp)
         old = nc
         DO j = 1, old
            r = components(j)%rank
            IF (r == 3) CYCLE
            a(:3, :) = components(j)%normal
            b(:3) = components(j)%offset
            b(4:6) = -tau
            CALL saturated_normals(a, rows, nr, basis, status)
            IF (status /= 0) RETURN
            ! An intersection of the same dimension is either this connected torus or empty.
            IF (nr == r) CYCLE
            CALL intersect_components(a, b, rows, nr, basis, components, nc, tolerance, status)
            IF (status /= 0) RETURN
         END DO
      END DO

      ALLOCATE (work(nc), expected(ng), fixed(ng))
      nf = 0
      DO j = 1, nc
         r = components(j)%rank
         CALL affine_solution(components(j)%normal, [1, 2, 3], r, right, nullspace, status)
         IF (status /= 0) RETURN
         point = MATMUL(right, components(j)%offset)
         DO g = 1, ng
            a(:3, :) = rotations(:, :, mapping(g))
            DO i = 1, 3
               a(i, i) = a(i, i) - 1
            END DO
            delta = MATMUL(REAL(a(:3, :), dp), point) + translations(:, mapping(g))
            expected(g) = MAXVAL(ABS(delta - ANINT(delta))) <= tolerance &
                          .AND. MAXVAL(ABS(MATMUL(REAL(a(:3, :), dp), nullspace))) <= tolerance
         END DO
         DO trial = 1, 32
            variation = SQRT([2.0_dp, 3.0_dp, 5.0_dp]) + trial*SQRT([7.0_dp, 11.0_dp, 13.0_dp])
            seed = MODULO(point + MATMUL(nullspace, variation), 1.0_dp)
            DO g = 1, ng
               delta = MATMUL(REAL(rotations(:, :, mapping(g)), dp), seed) + translations(:, mapping(g)) - seed
               fixed(g) = MAXVAL(ABS(delta - ANINT(delta))) <= tolerance
            END DO
            IF (ALL(fixed .EQV. expected)) EXIT
         END DO
         status = -4
         IF (trial > 32 .OR. COUNT(fixed) == 0) RETURN
         IF (MODULO(ng, COUNT(fixed)) /= 0) RETURN
         equivalent = .FALSE.
         DO g = 1, ng
            image = MATMUL(REAL(rotations(:, :, mapping(g)), dp), seed) + translations(:, mapping(g))
            direction = MATMUL(REAL(rotations(:, :, mapping(g)), dp), nullspace)
            DO i = 1, nf
               IF (work(i)%dimension /= 3 - r) CYCLE
               IF (work(i)%multiplicity /= ng/COUNT(fixed)) CYCLE
               mapped = MATMUL(REAL(work(i)%normal, dp), direction)
               IF (MAXVAL(ABS(mapped)) > tolerance) CYCLE
               delta = MATMUL(REAL(work(i)%normal, dp), image) - work(i)%offset
               IF (MAXVAL(ABS(delta - ANINT(delta))) > tolerance) CYCLE
               equivalent = .TRUE.
               EXIT
            END DO
            IF (equivalent) EXIT
         END DO
         IF (equivalent) CYCLE
         nf = nf + 1
         work(nf)%normal = components(j)%normal
         work(nf)%offset = components(j)%offset
         work(nf)%position = seed
         work(nf)%dimension = 3 - r
         work(nf)%multiplicity = ng/COUNT(fixed)
         work(nf)%stabilizer = PACK(mapping, fixed)
      END DO
      ALLOCATE (specializes(nf))
      DO j = 1, nf
         specializes = .FALSE.
         DO i = 1, nf
            IF (work(i)%dimension >= work(j)%dimension) CYCLE
            CALL affine_solution(work(i)%normal, [1, 2, 3], 3 - work(i)%dimension, right, nullspace, status)
            IF (status /= 0) RETURN
            DO g = 1, ng
               image = MATMUL(REAL(rotations(:, :, mapping(g)), dp), work(i)%position) + translations(:, mapping(g))
               direction = MATMUL(REAL(rotations(:, :, mapping(g)), dp), nullspace)
               IF (MAXVAL(ABS(MATMUL(REAL(work(j)%normal, dp), direction))) > tolerance) CYCLE
               delta = MATMUL(REAL(work(j)%normal, dp), image) - work(j)%offset
               IF (MAXVAL(ABS(delta - ANINT(delta))) > tolerance) CYCLE
               specializes(i) = .TRUE.
               EXIT
            END DO
         END DO
         work(j)%specializations = PACK([(i, i=1, nf)], specializes)
         work(j)%maximal_site_symmetry = .NOT. ANY(specializes)
      END DO
      families = work(:nf)
      status = 0
   END SUBROUTINE wyckoff_families

! **************************************************************************************************
!> \brief Canonical saturated integer row lattice, rather than the possibly nonprimitive input lattice.
!> \param a Integer normals
!> \param rows Independent original rows
!> \param rank Row rank
!> \param basis Primitive normals of each connected solution torus, with zero unused rows
!> \param status Zero on success, -3 if exact intermediate integer arithmetic cannot be guaranteed
!> \note Saturation is essential: 2*x=0 modulo one has two components, each with normal x.
! **************************************************************************************************
   SUBROUTINE saturated_normals(a, rows, rank, basis, status)
      INTEGER, INTENT(IN)                                :: a(:, :)
      INTEGER, INTENT(OUT)                               :: rows(3), rank, basis(3, 3), status

      INTEGER                                            :: first, i
      INTEGER(int_8)                                     :: b(3, 3), d(3), divisor, u, v(3), w(3), z

      status = -3
      basis = 0
      rows = 0
      rank = 0
      IF (MAXVAL(ABS(INT(a, int_8))) > 100000_int_8) RETURN
      d = 0
      v = 0
      DO i = 1, SIZE(a, 1)
         w = INT(a(i, :), int_8)
         IF (ALL(w == 0)) CYCLE
         IF (rank == 0) THEN
            rows(1) = i
            rank = 1
            v = w
         ELSE IF (rank == 1) THEN
            d = [v(2)*w(3) - v(3)*w(2), v(3)*w(1) - v(1)*w(3), v(1)*w(2) - v(2)*w(1)]
            IF (ALL(d == 0)) CYCLE
            rows(2) = i
            rank = 2
         ELSE
            IF (SUM(d*w) == 0) CYCLE
            rows(3) = i
            rank = 3
            EXIT
         END IF
      END DO
      b = 0
      SELECT CASE (rank)
      CASE (1)
         CALL primitive_vector(v)
         b(1, :) = v
      CASE (2)
         CALL primitive_vector(d)
         IF (d(1) == 0 .AND. d(2) == 0) THEN
            b(1, 1) = 1
            b(2, 2) = 1
         ELSE
            CALL bezout(d(1), d(2), divisor, u, z)
            b(1, :) = [d(2)/divisor, -d(1)/divisor, 0_int_8]
            ! Together these span the full integer kernel of the primitive direction d.
            IF (MAXVAL(ABS([u, z])) > 100000_int_8 .OR. ABS(d(3)) > 100000_int_8) RETURN
            b(2, :) = [u*d(3), z*d(3), -divisor]
         END IF
         DO i = 1, 2
            first = FINDLOC(b(i, :) /= 0, .TRUE., DIM=1)
            IF (b(i, first) < 0) b(i, :) = -b(i, :)
         END DO
      CASE (3)
         DO i = 1, 3
            b(i, i) = 1
         END DO
      END SELECT
      IF (MAXVAL(ABS(b)) > 100000_int_8) RETURN
      basis = INT(b)
      status = 0
   END SUBROUTINE saturated_normals

! **************************************************************************************************
!> \brief Primitive vector with a deterministic sign.
!> \param v Nonzero integer vector
! **************************************************************************************************
   SUBROUTINE primitive_vector(v)
      INTEGER(int_8), INTENT(INOUT)                      :: v(3)

      INTEGER                                            :: i
      INTEGER(int_8)                                     :: divisor, next, u, z

      divisor = 0
      DO i = 1, 3
         CALL bezout(divisor, v(i), next, u, z)
         divisor = next
      END DO
      v = v/divisor
      i = FINDLOC(v /= 0, .TRUE., DIM=1)
      IF (v(i) < 0) v = -v
   END SUBROUTINE primitive_vector

! **************************************************************************************************
!> \brief Euclidean algorithm, with u*a+v*b=gcd(a,b)>=0.
!> \param a First integer
!> \param b Second integer
!> \param divisor Nonnegative greatest common divisor
!> \param u First coefficient
!> \param v Second coefficient
! **************************************************************************************************
   SUBROUTINE bezout(a, b, divisor, u, v)
      INTEGER(int_8), INTENT(IN)                         :: a, b
      INTEGER(int_8), INTENT(OUT)                        :: divisor, u, v

      INTEGER(int_8)                                     :: q, r, rr, s, ss, t, tmp, tt

      r = ABS(a)
      rr = ABS(b)
      s = 1
      ss = 0
      t = 0
      tt = 1
      DO WHILE (rr /= 0)
         q = r/rr
         tmp = r - q*rr
         r = rr
         rr = tmp
         tmp = s - q*ss
         s = ss
         ss = tmp
         tmp = t - q*tt
         t = tt
         tt = tmp
      END DO
      divisor = r
      u = s*SIGN(1_int_8, a)
      v = t*SIGN(1_int_8, b)
   END SUBROUTINE bezout

! **************************************************************************************************
!> \brief Right inverse and nullspace directions of at most three independent affine equations.
!> \param a Integer matrix
!> \param rows Independent rows, as returned by saturated_normals
!> \param rank Number of independent rows
!> \param right Particular solution is right(:,1:rank)*rhs(rows(1:rank))
!> \param nullspace Nullspace vectors, with zero unused columns
!> \param status Zero on success, -3 for numerically unstable elimination
! **************************************************************************************************
   SUBROUTINE affine_solution(a, rows, rank, right, nullspace, status)
      INTEGER, INTENT(IN)                                :: a(:, :), rows(3), rank
      REAL(dp), INTENT(OUT)                              :: right(3, 3), nullspace(3, 3)
      INTEGER, INTENT(OUT)                               :: status

      INTEGER                                            :: col, i, j, n, p, pivots(3)
      REAL(dp)                                           :: scale, tmp(6), work(3, 6)

      right = 0.0_dp
      nullspace = 0.0_dp
      work = 0.0_dp
      pivots = 0
      DO i = 1, rank
         work(i, :3) = REAL(a(rows(i), :), dp)
         work(i, 3 + i) = 1.0_dp
      END DO
      n = 0
      DO col = 1, 3
         IF (n == rank) EXIT
         p = MAXLOC(ABS(work(n + 1:rank, col)), DIM=1) + n
         IF (ABS(work(p, col)) < 64*EPSILON(1.0_dp)) CYCLE
         n = n + 1
         tmp = work(n, :)
         work(n, :) = work(p, :)
         work(p, :) = tmp
         scale = work(n, col)
         work(n, :) = work(n, :)/scale
         DO i = 1, rank
            IF (i == n) CYCLE
            scale = work(i, col)
            work(i, :) = work(i, :) - scale*work(n, :)
         END DO
         pivots(n) = col
      END DO
      status = -3
      IF (n /= rank) RETURN
      DO i = 1, rank
         right(pivots(i), :) = work(i, 4:6)
      END DO
      j = 0
      DO col = 1, 3
         IF (ANY(pivots == col)) CYCLE
         j = j + 1
         nullspace(col, j) = 1.0_dp
         DO i = 1, rank
            nullspace(pivots(i), j) = -work(i, col)
         END DO
      END DO
      status = 0
   END SUBROUTINE affine_solution

! **************************************************************************************************
!> \brief Split an affine congruence system into connected solution tori and append new components.
!> \param a Combined integer normals
!> \param b Combined offsets modulo one
!> \param rows Independent input rows
!> \param rank Input rank
!> \param basis Saturated normals from saturated_normals
!> \param components Unique components accumulated so far
!> \param nc Number of accumulated components
!> \param tolerance Geometry tolerance
!> \param status Zero on success, -3 for arithmetic or work limits
!> \note Enumerating integer right-hand sides over the unit cube covers every torus component;
!>       particular solutions need not lie inside the cube. Remaining rows reject inconsistent ones.
! **************************************************************************************************
   SUBROUTINE intersect_components(a, b, rows, rank, basis, components, nc, tolerance, status)
      INTEGER, INTENT(IN)                                :: a(:, :)
      REAL(dp), INTENT(IN)                               :: b(:)
      INTEGER, INTENT(IN)                                :: rows(3), rank, basis(3, 3)
      TYPE(affine_torus_type), ALLOCATABLE, &
         INTENT(INOUT)                                   :: components(:)
      INTEGER, INTENT(INOUT)                             :: nc
      REAL(dp), INTENT(IN)                               :: tolerance
      INTEGER, INTENT(OUT)                               :: status

      INTEGER                                            :: hi(3), i, j, lo(3), n1, n2, n3
      LOGICAL                                            :: found
      REAL(dp)                                           :: c(3), nullspace(3, 3), point(3), rhs(3), &
                                                            right(3, 3), shifts(3), delta(SIZE(b))
      TYPE(affine_torus_type), ALLOCATABLE               :: larger(:)

      CALL affine_solution(a, rows, rank, right, nullspace, status)
      IF (status /= 0) RETURN
      lo = 0
      hi = 0
      rhs = 0.0_dp
      DO i = 1, rank
         j = rows(i)
         rhs(i) = b(j)
         lo(i) = CEILING(SUM(MIN(a(j, :), 0)) - b(j) - tolerance)
         hi(i) = FLOOR(SUM(MAX(a(j, :), 0)) - b(j) + tolerance)
      END DO
      status = -3
      IF (PRODUCT(INT(MAX(0, hi - lo + 1), int_8)) > 1000000_int_8) RETURN
      DO n3 = lo(3), hi(3)
         DO n2 = lo(2), hi(2)
            DO n1 = lo(1), hi(1)
               point = MATMUL(right, rhs + REAL([n1, n2, n3], dp))
               delta = MATMUL(REAL(a, dp), point) - b
               IF (MAXVAL(ABS(delta - ANINT(delta))) > tolerance) CYCLE
               c = MODULO(MATMUL(REAL(basis, dp), point), 1.0_dp)
               found = .FALSE.
               DO i = 1, nc
                  IF (components(i)%rank /= rank) CYCLE
                  IF (ANY(components(i)%normal /= basis)) CYCLE
                  shifts = components(i)%offset - c
                  IF (MAXVAL(ABS(shifts - ANINT(shifts))) > tolerance) CYCLE
                  found = .TRUE.
                  EXIT
               END DO
               IF (found) CYCLE
               IF (nc == SIZE(components)) THEN
                  IF (nc >= 65536) RETURN
                  ALLOCATE (larger(2*nc))
                  larger(:nc) = components
                  CALL MOVE_ALLOC(larger, components)
               END IF
               nc = nc + 1
               components(nc)%rank = rank
               components(nc)%normal = basis
               components(nc)%offset = c
            END DO
         END DO
      END DO
      status = 0
   END SUBROUTINE intersect_components
END MODULE topology_wyckoff
