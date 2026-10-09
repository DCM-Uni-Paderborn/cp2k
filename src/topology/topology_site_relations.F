!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Site-specialization induction witnesses, independent of sampled reciprocal characters.
!> \author Thomas D. Kuehne, tkuehne@cp2k.org
! **************************************************************************************************
MODULE topology_site_relations
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE topology_band_representations,   ONLY: atomic_band_characters,&
                                              atomic_band_construct,&
                                              atomic_band_site_data,&
                                              atomic_band_type
   USE topology_integer_lattice,        ONLY: integer_lattice_factor,&
                                              integer_lattice_solve,&
                                              integer_lattice_type
   USE topology_wyckoff,                ONLY: wyckoff_family_type

   IMPLICIT NONE
   PRIVATE
   TYPE, PUBLIC :: site_relation_type
      INTEGER :: source = 0, destination = 0, operation = 0, shift(3) = 0
      REAL(dp) :: image(3) = 0.0_dp
      INTEGER, ALLOCATABLE :: coefficients(:, :)
   END TYPE site_relation_type
   TYPE, PUBLIC :: site_reduction_type
      INTEGER, ALLOCATABLE :: representative(:), expansion(:, :), witness_edge(:), witness_local(:)
   END TYPE site_reduction_type
   PUBLIC :: site_specialization, site_induction_reduce
CONTAINS

! **************************************************************************************************
!> \brief Compose only certified real-space equivalences and direct-sum induction witnesses.
!> \param relations Independently certified site-specialization induction matrices
!> \param first First global column for each family, including the final sentinel
!> \param dimensions Positive whole-band dimensions of the columns
!> \param maximal Maximal-site flags for the generated families
!> \param reduction Canonical equivalent columns and nonnegative expansion into unreduced generators
!> \param status Zero success, negative malformed/inconsistent/overflowing witness
!> \note No reciprocal signatures enter this construction. An unreduced generator is not by itself
!>       a completeness certificate for an elementary-band catalogue.
! **************************************************************************************************
   SUBROUTINE site_induction_reduce(relations, first, dimensions, maximal, reduction, status)
      TYPE(site_relation_type), INTENT(IN)               :: relations(:)
      INTEGER, INTENT(IN)                                :: first(:), dimensions(:)
      LOGICAL, INTENT(IN)                                :: maximal(:)
      TYPE(site_reduction_type), INTENT(OUT)             :: reduction
      INTEGER, INTENT(OUT)                               :: status

      INTEGER                                            :: a, b, c, d, e, i, j, k, n, ns, root, &
                                                            source, TARGET
      INTEGER(int_8)                                     :: rank_sum
      INTEGER(int_8), ALLOCATABLE                        :: column(:)
      INTEGER, ALLOCATABLE                               :: labels(:), order(:), preferred(:)
      LOGICAL, ALLOCATABLE                               :: at_maximum(:)
      TYPE(site_reduction_type)                          :: work

      status = -1
      n = SIZE(dimensions)
      ns = SIZE(maximal)
      IF (n < 1 .OR. SIZE(first) /= ns + 1) RETURN
      IF (first(1) /= 1 .OR. first(ns + 1) /= n + 1) RETURN
      IF (ANY(first(2:) <= first(:ns)) .OR. ANY(dimensions <= 0)) RETURN
      ALLOCATE (labels(n), preferred(n), order(n), at_maximum(n), column(n))
      labels = [(i, i=1, n)]
      DO i = 1, ns
         at_maximum(first(i):first(i + 1) - 1) = maximal(i)
      END DO
      DO e = 1, SIZE(relations)
         source = relations(e)%source
         TARGET = relations(e)%destination
         IF (MIN(source, TARGET) < 1 .OR. MAX(source, TARGET) > ns) RETURN
         IF (source == TARGET .OR. .NOT. ALLOCATED(relations(e)%coefficients)) RETURN
         a = first(source)
         b = first(source + 1) - 1
         c = first(TARGET)
         d = first(TARGET + 1) - 1
         IF (ANY(SHAPE(relations(e)%coefficients) /= [d - c + 1, b - a + 1])) RETURN
         IF (ANY(relations(e)%coefficients < 0)) RETURN
         DO j = 1, b - a + 1
            rank_sum = 0_int_8
            DO i = 1, d - c + 1
               rank_sum = rank_sum + INT(dimensions(c + i - 1), int_8)*INT(relations(e)%coefficients(i, j), int_8)
               IF (rank_sum > INT(dimensions(a + j - 1), int_8)) RETURN
            END DO
            IF (rank_sum /= INT(dimensions(a + j - 1), int_8)) RETURN
            IF (SUM(INT(relations(e)%coefficients(:, j), int_8)) /= 1) CYCLE
            i = FINDLOC(relations(e)%coefficients(:, j), 1, DIM=1)
            root = MIN(labels(a + j - 1), labels(c + i - 1))
            k = MAX(labels(a + j - 1), labels(c + i - 1))
            WHERE (labels == k) labels = root
         END DO
      END DO
      ! Prefer a maximal-site representative, without identifying columns by sampled characters.
      preferred = 0
      DO i = 1, n
         root = labels(i)
         IF (preferred(root) == 0) THEN
            preferred(root) = i
         ELSE IF (at_maximum(i) .AND. .NOT. at_maximum(preferred(root))) THEN
            preferred(root) = i
         END IF
      END DO
      ALLOCATE (work%representative(n), work%expansion(n, n), work%witness_edge(n), work%witness_local(n))
      work%representative = preferred(labels)
      work%expansion = 0
      work%witness_edge = 0
      work%witness_local = 0
      DO e = 1, SIZE(relations)
         a = first(relations(e)%source)
         DO j = 1, SIZE(relations(e)%coefficients, 2)
            IF (SUM(INT(relations(e)%coefficients(:, j), int_8)) <= 1) CYCLE
            root = work%representative(a + j - 1)
            IF (work%witness_edge(root) /= 0) CYCLE
            work%witness_edge(root) = e
            work%witness_local(root) = j
         END DO
      END DO
      order = [(i, i=1, n)]
      DO i = 2, n
         k = order(i)
         j = i - 1
         DO WHILE (j >= 1)
            IF (dimensions(order(j)) <= dimensions(k)) EXIT
            order(j + 1) = order(j)
            j = j - 1
         END DO
         order(j + 1) = k
      END DO
      ! Every summand of a nontrivial decomposition has strictly smaller positive band dimension.
      DO i = 1, n
         root = order(i)
         IF (work%representative(root) /= root) CYCLE
         e = work%witness_edge(root)
         IF (e == 0) THEN
            work%expansion(root, root) = 1
            CYCLE
         END IF
         j = work%witness_local(root)
         c = first(relations(e)%destination)
         column = 0_int_8
         DO k = 1, SIZE(relations(e)%coefficients, 1)
            a = relations(e)%coefficients(k, j)
            IF (a == 0) CYCLE
            b = work%representative(c + k - 1)
            status = -2
            IF (dimensions(b) >= dimensions(root)) RETURN
            column = column + INT(a, int_8)*INT(work%expansion(:, b), int_8)
            status = -3
            IF (ANY(column > INT(HUGE(0), int_8))) RETURN
         END DO
         work%expansion(:, root) = INT(column)
      END DO
      DO i = 1, n
         root = work%representative(i)
         work%expansion(:, i) = work%expansion(:, root)
         work%witness_edge(i) = work%witness_edge(root)
         work%witness_local(i) = work%witness_local(root)
      END DO
      status = -2
      DO j = 1, n
         rank_sum = 0_int_8
         DO i = 1, n
            rank_sum = rank_sum + INT(dimensions(i), int_8)*INT(work%expansion(i, j), int_8)
            IF (rank_sum > INT(dimensions(j), int_8)) RETURN
         END DO
         IF (rank_sum /= INT(dimensions(j), int_8)) RETURN
      END DO
      reduction = work
      status = 0
   END SUBROUTINE site_induction_reduce

! **************************************************************************************************
!> \brief Induce onsite (co)representations along a validated Wyckoff specialization.
!> \param cell Cartesian lattice vectors as columns
!> \param rotations Complete direct-space group rotations
!> \param translations Fractional coset representatives
!> \param antiunitary Full-group time-reversal flags
!> \param spinful Include physical spin factors
!> \param families Validated connected Wyckoff closures and generic seeds
!> \param bands Onsite representations constructed from exactly the same group and seeds
!> \param source Lower-symmetry family index
!> \param destination Higher-symmetry family in the source closure
!> \param relation Lifted endpoint and destination-by-source induction multiplicities
!> \param tolerance Geometry and character tolerance
!> \param status Zero success; negative inconsistent data or arithmetic/work failure
!> \note Finite site-group induction supplies a homotopy witness, not a fit of sampled k labels.
! **************************************************************************************************
   SUBROUTINE site_specialization(cell, rotations, translations, antiunitary, spinful, families, bands, &
                                  source, destination, relation, tolerance, status)
      REAL(dp), INTENT(IN)                               :: cell(3, 3)
      INTEGER, INTENT(IN)                                :: rotations(:, :, :)
      REAL(dp), INTENT(IN)                               :: translations(:, :)
      LOGICAL, INTENT(IN)                                :: antiunitary(:), spinful
      TYPE(wyckoff_family_type), INTENT(IN)              :: families(:)
      TYPE(atomic_band_type), INTENT(IN)                 :: bands(:)
      INTEGER, INTENT(IN)                                :: source, destination
      TYPE(site_relation_type), INTENT(OUT)              :: relation
      REAL(dp), INTENT(IN)                               :: tolerance
      INTEGER, INTENT(OUT)                               :: status

      COMPLEX(dp), ALLOCATABLE                           :: gamma_chars(:, :), local_chars(:, :), &
                                                            source_chars(:, :), target_chars(:, :)
      INTEGER                                            :: g, h, i, j, n, nt, orbit
      INTEGER(int_8), ALLOCATABLE                        :: lift(:)
      INTEGER, ALLOCATABLE :: local_counts(:, :), map_gamma(:), map_local(:), map_source(:), &
         map_target(:), members(:), source_permutation(:, :)
      LOGICAL                                            :: exists
      REAL(dp)                                           :: delta(3), error, x(3), y(3)
      REAL(dp), ALLOCATABLE                              :: positions(:, :), tau(:, :)
      TYPE(atomic_band_type)                             :: collision
      TYPE(integer_lattice_type)                         :: normal
      TYPE(site_relation_type)                           :: work

      status = -1
      n = SIZE(antiunitary)
      IF (.NOT. ieee_is_finite(tolerance)) RETURN
      IF (tolerance <= 0.0_dp .OR. tolerance >= 1.e-3_dp) RETURN
      IF (n < 1 .OR. ANY(SHAPE(rotations) /= [3, 3, n])) RETURN
      IF (ANY(SHAPE(translations) /= [3, n])) RETURN
      IF (.NOT. ALL(ieee_is_finite(cell)) .OR. .NOT. ALL(ieee_is_finite(translations))) RETURN
      IF (SIZE(bands) /= SIZE(families)) RETURN
      IF (MIN(source, destination) < 1 .OR. MAX(source, destination) > SIZE(families)) RETURN
      IF (.NOT. ALLOCATED(families(source)%specializations)) RETURN
      IF (.NOT. ALLOCATED(families(source)%stabilizer)) RETURN
      IF (.NOT. ALLOCATED(families(destination)%stabilizer)) RETURN
      IF (.NOT. ANY(families(source)%specializations == destination)) RETURN
      IF (ANY(families(source)%stabilizer < 1) .OR. ANY(families(source)%stabilizer > n)) RETURN
      IF (.NOT. ALL(ieee_is_finite(families(source)%position))) RETURN
      IF (.NOT. ALL(ieee_is_finite(families(destination)%position))) RETURN
      IF (MIN(families(source)%multiplicity, families(destination)%multiplicity) < 1) RETURN
      x = families(source)%position
      CALL integer_lattice_factor(INT(families(source)%normal, int_8), normal, status)
      IF (status /= 0) RETURN
      DO g = 1, n
         y = MATMUL(REAL(rotations(:, :, g), dp), families(destination)%position) + translations(:, g)
         delta = MATMUL(REAL(families(source)%normal, dp), x - y)
         IF (MAXVAL(ABS(delta - ANINT(delta))) > tolerance) CYCLE
         status = -3
         IF (MAXVAL(ABS(delta)) > REAL(HUGE(0), dp)) RETURN
         CALL integer_lattice_solve(normal, NINT(delta, int_8), exists, lift, status)
         IF (status /= 0) RETURN
         status = -2
         IF (.NOT. exists) RETURN
         status = -3
         IF (ANY(ABS(lift) > INT(HUGE(0), int_8))) RETURN
         y = y + REAL(lift, dp)
         ! A common lift fixes the entire straight displacement, not merely its periodic endpoints.
         DO i = 1, SIZE(families(source)%stabilizer)
            h = families(source)%stabilizer(i)
            delta = MATMUL(REAL(rotations(:, :, h), dp), x - y) - (x - y)
            IF (MAXVAL(ABS(delta)) > tolerance) EXIT
         END DO
         IF (i <= SIZE(families(source)%stabilizer)) CYCLE
         work%operation = g
         work%shift = INT(lift)
         EXIT
      END DO
      status = -2
      IF (work%operation == 0) RETURN
      work%image = y
      work%source = source
      work%destination = destination
      ALLOCATE (members(n), tau(3, n))
      nt = 0
      DO g = 1, n
         delta = MATMUL(REAL(rotations(:, :, g), dp), y) + translations(:, g) - y
         IF (MAXVAL(ABS(delta - ANINT(delta))) > tolerance) CYCLE
         nt = nt + 1
         members(nt) = g
         tau(:, nt) = y - MATMUL(REAL(rotations(:, :, g), dp), y)
      END DO
      IF (nt /= SIZE(families(destination)%stabilizer)) RETURN
      CALL atomic_band_construct(cell, rotations(:, :, members(:nt)), tau(:, :nt), antiunitary(members(:nt)), &
                                 x, spinful, collision, tolerance, error, status)
      IF (status /= 0) RETURN
      CALL atomic_band_site_data(collision, positions, map_local, local_chars)
      status = -2
      IF (SIZE(positions, 2)*families(destination)%multiplicity /= families(source)%multiplicity) RETURN
      CALL atomic_band_site_data(bands(source), positions, map_source, source_chars)
      IF (.NOT. ALLOCATED(map_source)) RETURN
      CALL site_counts(map_source, source_chars, members(map_local), local_chars, source_permutation, tolerance, status)
      IF (status /= 0) RETURN
      status = -2
      IF (SIZE(source_permutation, 1) /= SIZE(source_permutation, 2)) RETURN
      IF (ANY(SUM(source_permutation, DIM=1) /= 1) .OR. ANY(SUM(source_permutation, DIM=2) /= 1)) RETURN
      CALL atomic_band_characters(collision, [0.0_dp, 0.0_dp, 0.0_dp], map_gamma, gamma_chars, tolerance, status)
      IF (status /= 0) RETURN
      CALL atomic_band_site_data(bands(destination), positions, map_target, target_chars)
      status = -2
      IF (.NOT. ALLOCATED(map_target)) RETURN
      orbit = 0
      DO j = 1, SIZE(positions, 2)
         delta = y - positions(:, j)
         IF (MAXVAL(ABS(delta - ANINT(delta))) > tolerance) CYCLE
         orbit = j
         EXIT
      END DO
      status = -2
      IF (orbit == 0) RETURN
      CALL atomic_band_site_data(bands(destination), positions, map_target, target_chars, orbit)
      CALL site_counts(map_target, target_chars, members(map_gamma), gamma_chars, local_counts, tolerance, status)
      IF (status /= 0) RETURN
      ALLOCATE (work%coefficients(SIZE(local_counts, 1), SIZE(source_permutation, 1)))
      work%coefficients(:, :) = MATMUL(local_counts, TRANSPOSE(source_permutation))
      status = -2
      IF (ANY(MATMUL(NINT(REAL(target_chars(1, :), dp)), work%coefficients)*families(destination)%multiplicity &
              /= NINT(REAL(source_chars(1, :), dp))*families(source)%multiplicity)) RETURN
      relation = work
      status = 0
   END SUBROUTINE site_specialization

! **************************************************************************************************
!> \brief Decompose in validated full onsite corepresentation characters with their actual norms.
!> \param reference_map Original unitary operation IDs for reference rows
!> \param reference Validated onsite irrep/corepresentation characters
!> \param mapping Original unitary operation IDs for supplied rows
!> \param characters Characters to decompose, one representation per column
!> \param counts Nonnegative integer multiplicities, unallocated on failure
!> \param tolerance Absolute algebra tolerance
!> \param status Zero success, negative inconsistent mapping or characters
!> \note Wigner-paired/doubled corepresentations have character norms different from one.
! **************************************************************************************************
   SUBROUTINE site_counts(reference_map, reference, mapping, characters, counts, tolerance, status)
      INTEGER, INTENT(IN)                                :: reference_map(:)
      COMPLEX(dp), INTENT(IN)                            :: reference(:, :)
      INTEGER, INTENT(IN)                                :: mapping(:)
      COMPLEX(dp), INTENT(IN)                            :: characters(:, :)
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: counts(:, :)
      REAL(dp), INTENT(IN)                               :: tolerance
      INTEGER, INTENT(OUT)                               :: status

      COMPLEX(dp)                                        :: value
      COMPLEX(dp), ALLOCATABLE                           :: image(:)
      INTEGER                                            :: i, j, n
      INTEGER, ALLOCATABLE                               :: rows(:), work(:, :)
      REAL(dp)                                           :: norm

      status = -2
      n = SIZE(reference_map)
      IF (n < 1 .OR. SIZE(mapping) /= n) RETURN
      IF (SIZE(reference, 1) /= n .OR. SIZE(characters, 1) /= n) RETURN
      IF (MIN(SIZE(reference, 2), SIZE(characters, 2)) < 1) RETURN
      IF (.NOT. ALL(ieee_is_finite(REAL(reference, dp))) .OR. .NOT. ALL(ieee_is_finite(AIMAG(reference)))) RETURN
      IF (.NOT. ALL(ieee_is_finite(REAL(characters, dp))) .OR. .NOT. ALL(ieee_is_finite(AIMAG(characters)))) RETURN
      ALLOCATE (rows(n))
      DO i = 1, n
         rows(i) = FINDLOC(mapping, reference_map(i), DIM=1)
         IF (rows(i) == 0) RETURN
         IF (COUNT(mapping == reference_map(i)) /= 1) RETURN
         IF (COUNT(reference_map == reference_map(i)) /= 1) RETURN
      END DO
      ALLOCATE (work(SIZE(reference, 2), SIZE(characters, 2)), image(n))
      DO j = 1, SIZE(characters, 2)
         DO i = 1, SIZE(reference, 2)
            norm = SUM(ABS(reference(:, i))**2)
            IF (norm < 1.0_dp - tolerance) RETURN
            value = SUM(CONJG(reference(:, i))*characters(rows, j))/norm
            IF (ABS(AIMAG(value)) > tolerance) RETURN
            IF (ABS(REAL(value, dp)) >= REAL(HUGE(0), dp)) RETURN
            work(i, j) = NINT(REAL(value, dp))
            IF (work(i, j) < 0 .OR. ABS(value - work(i, j)) > tolerance) RETURN
         END DO
         image(:) = MATMUL(reference, CMPLX(work(:, j), 0, dp))
         IF (MAXVAL(ABS(image - characters(rows, j))) > tolerance) RETURN
      END DO
      CALL MOVE_ALLOC(work, counts)
      status = 0
   END SUBROUTINE site_counts
END MODULE topology_site_relations
