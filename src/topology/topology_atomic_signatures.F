!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Atomic signature matrices from all generated sites at an explicitly supplied set of k points.
!> \author Thomas D. Kuehne, tkuehne@cp2k.org
! **************************************************************************************************
MODULE topology_atomic_signatures
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE mathconstants,                   ONLY: twopi
   USE topology_band_representations,   ONLY: &
        atomic_band_character_polynomial, atomic_band_characters, atomic_band_construct, &
        atomic_band_site_data, atomic_band_site_snapshot, atomic_band_type, &
        atomic_character_compare, atomic_character_polynomial_type, atomic_site_data_type
   USE topology_compatibility,          ONLY: little_group_data_type
   USE topology_corepresentations,      ONLY: corepresentation_multiplicities,&
                                              corepresentation_table_type,&
                                              projective_corepresentations
   USE topology_little_group,           ONLY: little_group_factor_system,&
                                              projective_characters
   USE topology_site_relations,         ONLY: site_induction_reduce,&
                                              site_reduction_type,&
                                              site_relation_type,&
                                              site_specialization
   USE topology_star_transport,         ONLY: star_group_construct,&
                                              star_group_type,&
                                              star_irrep_map
   USE topology_symmetry,               ONLY: character_multiplicities
   USE topology_wyckoff,                ONLY: wyckoff_families,&
                                              wyckoff_family_type

   IMPLICIT NONE
   PRIVATE

   TYPE, PUBLIC :: atomic_signature_point_type
      TYPE(little_group_data_type) :: little
      TYPE(corepresentation_table_type) :: coreps
      INTEGER, ALLOCATABLE :: mapping(:), dimensions(:)
   END TYPE atomic_signature_point_type

   TYPE, PUBLIC :: atomic_signature_type
      TYPE(star_group_type) :: group
      TYPE(wyckoff_family_type), ALLOCATABLE :: families(:)
      TYPE(atomic_signature_point_type), ALLOCATABLE :: points(:)
      TYPE(atomic_site_data_type), ALLOCATABLE :: sites(:)
      COMPLEX(dp), ALLOCATABLE :: reference_spin(:, :, :)
      INTEGER, ALLOCATABLE :: matrix(:, :), row_begin(:), row_dimension(:)
      INTEGER, ALLOCATABLE :: column_site(:), column_local(:), band_dimension(:)
      ! Full invariant-torus unitary characters, not Bloch-bundle or site equivalence.
      INTEGER, ALLOCATABLE :: character_class(:)
      REAL(dp) :: character_residual = 0.0_dp
      TYPE(site_relation_type), ALLOCATABLE :: site_relations(:)
      TYPE(site_reduction_type) :: site_reduction
   END TYPE atomic_signature_type

   PUBLIC :: atomic_signature_construct, atomic_signature_decompose, atomic_signature_from_bands
   PUBLIC :: atomic_signature_star_map
CONTAINS

! **************************************************************************************************
!> \brief Star transport in the catalogue's full irrep/corepresentation coordinate basis.
!> \param catalogue Validated native reference data
!> \param source Source point index
!> \param destination Destination point index
!> \param operation Original full-group operation index
!> \param permutation Destination-by-source coordinate permutation; unallocated on failure
!> \param tolerance Algebra tolerance
!> \param residual Largest factor/character residual
!> \param status Zero on success; negative invalid/inconsistent data
! **************************************************************************************************
   SUBROUTINE atomic_signature_star_map(catalogue, source, destination, operation, permutation, tolerance, residual, status)
      TYPE(atomic_signature_type), INTENT(IN)            :: catalogue
      INTEGER, INTENT(IN)                                :: source, destination, operation
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: permutation(:, :)
      REAL(dp), INTENT(IN)                               :: tolerance
      REAL(dp), INTENT(OUT)                              :: residual
      INTEGER, INTENT(OUT)                               :: status

      INTEGER                                            :: j, n
      INTEGER, ALLOCATABLE                               :: restricted(:, :), unitary(:, :), &
                                                            work(:, :)

      status = -1
      residual = HUGE(1.0_dp)
      IF (.NOT. ALLOCATED(catalogue%points)) RETURN
      IF (MIN(source, destination) < 1 .OR. MAX(source, destination) > SIZE(catalogue%points)) RETURN
      CALL star_irrep_map(catalogue%group, catalogue%points(source)%little, catalogue%points(destination)%little, &
                          operation, unitary, tolerance, residual, status)
      IF (status /= 0) RETURN
      status = -2
      IF (ALLOCATED(catalogue%points(source)%coreps%restriction) .NEQV. &
          ALLOCATED(catalogue%points(destination)%coreps%restriction)) RETURN
      n = SIZE(catalogue%points(source)%dimensions)
      IF (SIZE(catalogue%points(destination)%dimensions) /= n) RETURN
      ALLOCATE (work(n, n))
      IF (ALLOCATED(catalogue%points(source)%coreps%restriction)) THEN
         ALLOCATE (restricted(SIZE(unitary, 1), n))
         restricted(:, :) = MATMUL(unitary, catalogue%points(source)%coreps%restriction)
         DO j = 1, n
            CALL corepresentation_multiplicities(catalogue%points(destination)%coreps, restricted(:, j), work(:, j), status)
            IF (status /= 0) RETURN
         END DO
      ELSE
         work(:, :) = unitary
      END IF
      status = -2
      IF (ANY(SUM(work, DIM=1) /= 1) .OR. ANY(SUM(work, DIM=2) /= 1)) RETURN
      IF (ANY(MATMUL(work, catalogue%points(source)%dimensions) /= catalogue%points(destination)%dimensions)) RETURN
      CALL MOVE_ALLOC(work, permutation)
      status = 0
   END SUBROUTINE atomic_signature_star_map

! **************************************************************************************************
!> \brief Build atomic integer signatures without eliminating sites or equivalent sampled columns.
!> \param cell Cartesian cell vectors as columns
!> \param rotations Complete group of integer direct-space rotations modulo lattice translations
!> \param translations Fractional coset representatives in the supplied cell/origin convention
!> \param antiunitary Physical time-reversal flags for the supplied operations
!> \param kpoints Explicit fractional reciprocal points, retained in input order
!> \param spinful Include double-group spin factors and physical Theta squared
!> \param catalogue Signature matrix, row tables and column provenance; empty on failure
!> \param tolerance Absolute geometry/character tolerance
!> \param status Zero on success; negative for invalid group, characters, dimensions or work limits
!> \param max_memory_mb Optional estimated matrix-storage budget in MiB, excluding geometric search
!> \note Rows are unitary irreps or full corepresentations when an antiunitary coset is present.
!>       Every column retains one site and onsite irrep/corepresentation across all k points.
!>       Includes nonmaximal sites and real-space induction witnesses; no band-bundle verdict is made.
!>       Supplied points need not cover all reciprocal strata. Centering sectors are not unfolded.
! **************************************************************************************************
   SUBROUTINE atomic_signature_construct(cell, rotations, translations, antiunitary, kpoints, spinful, &
                                         catalogue, tolerance, status, max_memory_mb)
      REAL(dp), INTENT(IN)                               :: cell(3, 3)
      INTEGER, INTENT(IN)                                :: rotations(:, :, :)
      REAL(dp), INTENT(IN)                               :: translations(:, :)
      LOGICAL, INTENT(IN)                                :: antiunitary(:)
      REAL(dp), INTENT(IN)                               :: kpoints(:, :)
      LOGICAL, INTENT(IN)                                :: spinful
      TYPE(atomic_signature_type), INTENT(OUT)           :: catalogue
      REAL(dp), INTENT(IN)                               :: tolerance
      INTEGER, INTENT(OUT)                               :: status
      REAL(dp), INTENT(IN), OPTIONAL                     :: max_memory_mb

      COMPLEX(dp), ALLOCATABLE                           :: characters(:, :)
      INTEGER                                            :: a, b, c, d, e, i, j, last, nc, np, nr, &
                                                            ns, p, row
      INTEGER, ALLOCATABLE                               :: first(:), mapping(:), site_dimensions(:)
      REAL(dp)                                           :: budget, group_size, memory, residual
      REAL(dp), ALLOCATABLE                              :: positions(:, :)
      TYPE(atomic_band_type), ALLOCATABLE                :: bands(:)
      TYPE(atomic_signature_type)                        :: work

      status = -1
      np = SIZE(kpoints, 2)
      IF (SIZE(kpoints, 1) /= 3 .OR. np == 0) RETURN
      IF (.NOT. ALL(ieee_is_finite(kpoints))) RETURN
      budget = HUGE(1.0_dp)
      IF (PRESENT(max_memory_mb)) budget = max_memory_mb
      IF (.NOT. ieee_is_finite(budget) .OR. budget <= 0.0_dp) RETURN
      group_size = REAL(SIZE(antiunitary), dp)
      memory = ((256.0_dp + 128.0_dp*np)*group_size**2 + 64.0_dp*group_size)/1024.0_dp**2
      status = -3
      IF (memory > budget) RETURN
      CALL star_group_construct(cell, rotations, translations, antiunitary, spinful, work%group, tolerance, status)
      IF (status /= 0) RETURN
      CALL wyckoff_families(cell, rotations, translations, antiunitary, work%families, tolerance, status)
      IF (status /= 0) RETURN
      ns = SIZE(work%families)
      memory = memory + (192.0_dp*group_size**2 + 256.0_dp*group_size + 512.0_dp)*ns/1024.0_dp**2
      status = -3
      IF (memory > budget) RETURN
      ALLOCATE (bands(ns), first(ns + 1), work%points(np), work%row_begin(np + 1), work%sites(ns))
      first(1) = 1
      DO i = 1, ns
         IF (i == 1) THEN
            CALL atomic_band_construct(cell, rotations, translations, antiunitary, work%families(i)%position, &
                                       spinful, bands(i), tolerance, residual, status, full_spin=work%reference_spin)
         ELSE
            CALL atomic_band_construct(cell, rotations, translations, antiunitary, work%families(i)%position, &
                                       spinful, bands(i), tolerance, residual, status)
         END IF
         IF (status /= 0) RETURN
         CALL atomic_band_site_data(bands(i), positions, mapping, characters)
         CALL atomic_band_site_snapshot(bands(i), work%sites(i))
         status = -2
         IF (SIZE(positions, 2) /= work%families(i)%multiplicity) RETURN
         first(i + 1) = first(i) + SIZE(characters, 2)
      END DO
      nc = first(ns + 1) - 1
      nr = 0
      DO i = 1, ns
         nr = nr + SIZE(work%families(i)%specializations)
      END DO
      memory = memory + (20.0_dp*REAL(nc, dp)**2 + 256.0_dp*nr + &
                         group_size*(512.0_dp*ns + 16.0_dp*nc) + 8.0_dp*nc)/1024.0_dp**2
      DO i = 1, ns
         DO j = 1, SIZE(work%families(i)%specializations)
            p = work%families(i)%specializations(j)
            memory = memory + 8.0_dp*(first(i + 1) - first(i))*(first(p + 1) - first(p))/1024.0_dp**2
         END DO
      END DO
      status = -3
      IF (memory > budget) RETURN
      ALLOCATE (work%column_site(nc), work%column_local(nc), work%band_dimension(nc))
      DO i = 1, ns
         CALL atomic_band_site_data(bands(i), positions, mapping, characters)
         site_dimensions = NINT(REAL(characters(1, :), dp))
         DO j = 1, SIZE(site_dimensions)
            p = first(i) + j - 1
            work%column_site(p) = i
            work%column_local(p) = j
            work%band_dimension(p) = work%families(i)%multiplicity*site_dimensions(j)
         END DO
      END DO
      ALLOCATE (work%site_relations(nr))
      e = 0
      DO i = 1, ns
         DO j = 1, SIZE(work%families(i)%specializations)
            p = work%families(i)%specializations(j)
            e = e + 1
            CALL site_specialization(cell, rotations, translations, antiunitary, spinful, work%families, bands, &
                                     i, p, work%site_relations(e), tolerance, status)
            IF (status /= 0) RETURN
         END DO
      END DO
      CALL site_induction_reduce(work%site_relations, first, work%band_dimension, work%families%maximal_site_symmetry, &
                                 work%site_reduction, status)
      IF (status /= 0) RETURN
      work%row_begin(1) = 1
      DO p = 1, np
         CALL signature_point(cell, rotations, translations, antiunitary, kpoints(:, p), spinful, &
                              work%points(p), tolerance, status)
         IF (status /= 0) RETURN
         work%row_begin(p + 1) = work%row_begin(p) + SIZE(work%points(p)%dimensions)
      END DO
      row = work%row_begin(np + 1) - 1
      status = -3
      IF (INT(row, int_8)*INT(nc, int_8) > 100000000_int_8) RETURN
      memory = memory + 16.0_dp*REAL(row, dp)*nc/1024.0_dp**2
      IF (memory > budget) RETURN
      ALLOCATE (work%matrix(row, nc), work%row_dimension(row))
      DO p = 1, np
         row = work%row_begin(p)
         last = work%row_begin(p + 1) - 1
         work%row_dimension(row:last) = work%points(p)%dimensions
         DO i = 1, ns
            CALL atomic_band_characters(bands(i), kpoints(:, p), mapping, characters, tolerance, status)
            IF (status /= 0) RETURN
            DO j = 1, SIZE(characters, 2)
               CALL atomic_signature_decompose(work%points(p), mapping, characters(:, j), &
                                               work%matrix(row:last, first(i) + j - 1), tolerance, status)
               IF (status /= 0) RETURN
            END DO
         END DO
         status = -2
         IF (ANY(MATMUL(work%row_dimension(row:last), work%matrix(row:last, :)) /= work%band_dimension)) RETURN
      END DO
      DO e = 1, nr
         i = work%site_relations(e)%source
         j = work%site_relations(e)%destination
         a = first(i)
         b = first(i + 1) - 1
         c = first(j)
         d = first(j + 1) - 1
         status = -2
         IF (ANY(MATMUL(work%matrix(:, c:d), work%site_relations(e)%coefficients) /= work%matrix(:, a:b))) RETURN
      END DO
      status = -2
      IF (ANY(MATMUL(work%matrix, work%site_reduction%expansion) /= work%matrix)) RETURN
      CALL character_classes(bands, work, tolerance, status)
      IF (status /= 0) RETURN
      catalogue = work
      status = 0
   END SUBROUTINE atomic_signature_construct

! **************************************************************************************************
!> \brief Classify unitary characters on the complete fixed tori without merging atomic columns.
!> \param bands Prepared onsite-induced representations in the same full-group convention
!> \param catalogue Sampled integer matrix and column provenance, updated with character classes
!> \param tolerance Dimension-scaled Fourier coefficient tolerance
!> \param status Zero on success; negative for invalid or nontransitive tolerance neighborhoods
!> \note Sampled signatures only prefilter candidate pairs. Neither site-reduction witnesses nor
!>       integer decompositions change, and character equality does not establish bundle equivalence.
! **************************************************************************************************
   SUBROUTINE character_classes(bands, catalogue, tolerance, status)
      TYPE(atomic_band_type), INTENT(IN)                 :: bands(:)
      TYPE(atomic_signature_type), INTENT(INOUT)         :: catalogue
      REAL(dp), INTENT(IN)                               :: tolerance
      INTEGER, INTENT(OUT)                               :: status

      INTEGER                                            :: a, b, i, j, n
      LOGICAL                                            :: equal
      LOGICAL, ALLOCATABLE                               :: equivalent(:, :)
      REAL(dp)                                           :: residual
      TYPE(atomic_character_polynomial_type), &
         ALLOCATABLE                                     :: polynomials(:)

      n = SIZE(catalogue%matrix, 2)
      ALLOCATE (polynomials(SIZE(bands)), equivalent(n, n), catalogue%character_class(n))
      DO i = 1, SIZE(bands)
         CALL atomic_band_character_polynomial(bands(i), polynomials(i), status)
         IF (status /= 0) RETURN
      END DO
      equivalent = .FALSE.
      catalogue%character_residual = 0
      DO j = 1, n
         equivalent(j, j) = .TRUE.
         b = catalogue%column_site(j)
         DO i = 1, j - 1
            IF (catalogue%band_dimension(i) /= catalogue%band_dimension(j)) CYCLE
            IF (ANY(catalogue%matrix(:, i) /= catalogue%matrix(:, j))) CYCLE
            a = catalogue%column_site(i)
            CALL atomic_character_compare(polynomials(a), catalogue%column_local(i), &
                                          polynomials(b), catalogue%column_local(j), equal, tolerance, residual, status)
            IF (status /= 0) RETURN
            equivalent(i, j) = equal
            equivalent(j, i) = equal
            IF (equal) catalogue%character_residual = MAX(catalogue%character_residual, residual)
         END DO
         catalogue%character_class(j) = FINDLOC(equivalent(:, j), .TRUE., DIM=1)
      END DO
      status = -2
      DO j = 1, n
         i = catalogue%character_class(j)
         IF (ANY(equivalent(:, j) .NEQV. equivalent(:, i))) RETURN
      END DO
      status = 0
   END SUBROUTINE character_classes

! **************************************************************************************************
!> \brief Build native little-group tables once per k point, sharing the compatibility convention.
!> \param cell Cartesian cell
!> \param rotations Complete group rotations
!> \param translations Fractional translations
!> \param antiunitary Antiunitary flags
!> \param k Reciprocal point
!> \param spinful Use double-group factors
!> \param point Prepared tables and original unitary operation indices
!> \param tolerance Absolute geometry/character tolerance
!> \param status Zero on success
! **************************************************************************************************
   SUBROUTINE signature_point(cell, rotations, translations, antiunitary, k, spinful, point, tolerance, status)
      REAL(dp), INTENT(IN)                               :: cell(3, 3)
      INTEGER, INTENT(IN)                                :: rotations(:, :, :)
      REAL(dp), INTENT(IN)                               :: translations(:, :)
      LOGICAL, INTENT(IN)                                :: antiunitary(:)
      REAL(dp), INTENT(IN)                               :: k(3)
      LOGICAL, INTENT(IN)                                :: spinful
      TYPE(atomic_signature_point_type), INTENT(OUT)     :: point
      REAL(dp), INTENT(IN)                               :: tolerance
      INTEGER, INTENT(OUT)                               :: status

      INTEGER                                            :: g, n, nu
      INTEGER, ALLOCATABLE                               :: index(:), mapping(:), products(:, :), &
                                                            unitary(:)
      REAL(dp)                                           :: residual

      CALL little_group_factor_system(cell, rotations, translations, antiunitary, k, spinful, mapping, &
                                      point%little%product, point%little%factors, tolerance, status)
      IF (status /= 0) RETURN
      point%little%cell = cell
      point%little%k = k
      point%little%spinful = spinful
      point%little%rotations = rotations(:, :, mapping)
      point%little%translations = translations(:, mapping)
      point%little%antiunitary = antiunitary(mapping)
      n = SIZE(mapping)
      unitary = PACK([(g, g=1, n)],.NOT. point%little%antiunitary)
      nu = SIZE(unitary)
      point%mapping = mapping(unitary)
      ALLOCATE (INDEX(n), products(nu, nu))
      index = 0
      DO g = 1, nu
         INDEX(unitary(g)) = g
      END DO
      DO g = 1, nu
         products(:, g) = INDEX(point%little%product(unitary, unitary(g)))
      END DO
      CALL projective_characters(products, point%little%factors(unitary, unitary), point%little%characters, &
                                 tolerance, residual, status)
      IF (status /= 0) RETURN
      ALLOCATE (point%little%multiplicity(SIZE(point%little%characters, 2)))
      point%little%multiplicity = 0
      IF (nu == n) THEN
         point%dimensions = NINT(REAL(point%little%characters(1, :), dp))
      ELSE
         CALL projective_corepresentations(point%little%product, point%little%factors, point%little%antiunitary, &
                                           point%little%characters, point%coreps, tolerance, residual, status)
         IF (status /= 0) RETURN
         point%dimensions = point%coreps%dimension
      END IF
   END SUBROUTINE signature_point

! **************************************************************************************************
!> \brief Decompose a unitary character vector using the catalogue's operation and irrep convention.
!> \param point Successfully prepared signature point
!> \param mapping Original unitary operation indices, in any order but with no duplicates/missing rows
!> \param chi Characters in exactly the same spin/Bloch coset gauge as point
!> \param counts Nonnegative irreducible/corepresentation counts; -1 on failure
!> \param tolerance Absolute character tolerance
!> \param status Zero on success; negative for invalid mapping, character or antiunitary pairing
!> \note Mapping refers to the original complete operation list used for construction. This routine
!>       checks decomposition, not physical subspace closure or symmetry of a Hamiltonian.
! **************************************************************************************************
   SUBROUTINE atomic_signature_decompose(point, mapping, chi, counts, tolerance, status)
      TYPE(atomic_signature_point_type), INTENT(IN)      :: point
      INTEGER, INTENT(IN)                                :: mapping(:)
      COMPLEX(dp), INTENT(IN)                            :: chi(:)
      INTEGER, INTENT(OUT)                               :: counts(:)
      REAL(dp), INTENT(IN)                               :: tolerance
      INTEGER, INTENT(OUT)                               :: status

      INTEGER                                            :: g, n
      INTEGER, ALLOCATABLE                               :: order(:), unitary_counts(:), work(:)

      status = -1
      counts = -1
      IF (.NOT. ALLOCATED(point%mapping) .OR. .NOT. ALLOCATED(point%dimensions)) RETURN
      IF (.NOT. ALLOCATED(point%little%characters)) RETURN
      n = SIZE(point%mapping)
      IF (SIZE(mapping) /= n .OR. SIZE(chi) /= n .OR. SIZE(counts) /= SIZE(point%dimensions)) RETURN
      ALLOCATE (order(n), unitary_counts(SIZE(point%little%characters, 2)), work(SIZE(counts)))
      DO g = 1, n
         IF (COUNT(mapping == point%mapping(g)) /= 1) RETURN
         order(g) = FINDLOC(mapping, point%mapping(g), DIM=1)
      END DO
      CALL character_multiplicities(chi(order), point%little%characters, unitary_counts, tolerance, status)
      IF (status /= 0) RETURN
      IF (ALLOCATED(point%coreps%restriction)) THEN
         CALL corepresentation_multiplicities(point%coreps, unitary_counts, work, status)
         IF (status /= 0) RETURN
      ELSE
         IF (SIZE(work) /= SIZE(unitary_counts)) THEN
            status = -1
            RETURN
         END IF
         work = unitary_counts
      END IF
      counts = work
      status = 0
   END SUBROUTINE atomic_signature_decompose

! **************************************************************************************************
!> \brief Map validated band characters to atomic rows, including exact coset/Bloch rephasing.
!> \param catalogue Atomic reference at the same ordered reciprocal points and cell
!> \param points Physically validated band representations; characters follow unitary operation order
!> \param TARGET Integer band signature in the catalogue's row convention; unallocated on failure
!> \param tolerance Geometry/factor/character tolerance; no gauge fitting is performed
!> \param status Zero on success, negative for inconsistent groups, factors, pairing or band ranks
!> \note All unitary and antiunitary operations must match bijectively modulo lattice translations.
!>       Cell/origin/spin conventions must agree. This routine does not check a Hamiltonian or gap.
! **************************************************************************************************
   SUBROUTINE atomic_signature_from_bands(catalogue, points, TARGET, tolerance, status)
      TYPE(atomic_signature_type), INTENT(IN)            :: catalogue
      TYPE(little_group_data_type), INTENT(IN)           :: points(:)
      INTEGER(int_8), ALLOCATABLE, INTENT(OUT)           :: target(:)
      REAL(dp), INTENT(IN)                               :: tolerance
      INTEGER, INTENT(OUT)                               :: status

      COMPLEX(dp)                                        :: expected, second_phase
      COMPLEX(dp), ALLOCATABLE                           :: character(:), old_character(:), phase(:)
      INTEGER                                            :: a, b, g, h, j, n, nu, p
      INTEGER(int_8)                                     :: first_rank, rank
      INTEGER(int_8), ALLOCATABLE                        :: work(:)
      INTEGER, ALLOCATABLE                               :: counts(:), inverse_unitary(:), order(:), &
                                                            unitary(:)
      REAL(dp)                                           :: angle, delta(3)

      status = -1
      IF (.NOT. ALLOCATED(catalogue%points) .OR. .NOT. ALLOCATED(catalogue%matrix)) RETURN
      IF (.NOT. ALLOCATED(catalogue%row_begin) .OR. .NOT. ALLOCATED(catalogue%row_dimension)) RETURN
      IF (.NOT. ieee_is_finite(tolerance)) RETURN
      IF (tolerance <= 0.0_dp .OR. tolerance >= 0.01_dp) RETURN
      IF (SIZE(points) /= SIZE(catalogue%points) .OR. SIZE(points) == 0) RETURN
      ALLOCATE (work(SIZE(catalogue%matrix, 1)))
      first_rank = -1
      DO p = 1, SIZE(points)
         ASSOCIATE (reference => catalogue%points(p)%little, band => points(p))
            status = -1
            IF (.NOT. ALL(ieee_is_finite(band%cell)) .OR. .NOT. ALL(ieee_is_finite(band%k))) RETURN
            IF (MAXVAL(ABS(band%cell - reference%cell)) > tolerance) RETURN
            IF (MAXVAL(ABS(band%k - reference%k)) > tolerance) RETURN
            IF (band%spinful .NEQV. reference%spinful) RETURN
            IF (.NOT. ALLOCATED(band%antiunitary) .OR. .NOT. ALLOCATED(band%rotations)) RETURN
            IF (.NOT. ALLOCATED(band%translations) .OR. .NOT. ALLOCATED(band%product)) RETURN
            IF (.NOT. ALLOCATED(band%factors) .OR. .NOT. ALLOCATED(band%characters)) RETURN
            IF (.NOT. ALLOCATED(band%multiplicity)) RETURN
            n = SIZE(reference%antiunitary)
            nu = COUNT(.NOT. reference%antiunitary)
            IF (SIZE(band%antiunitary) /= n .OR. COUNT(.NOT. band%antiunitary) /= nu) RETURN
            IF (ANY(SHAPE(band%rotations) /= [3, 3, n]) .OR. ANY(SHAPE(band%translations) /= [3, n])) RETURN
            IF (ANY(SHAPE(band%product) /= [n, n]) .OR. ANY(SHAPE(band%factors) /= [n, n])) RETURN
            IF (SIZE(band%characters, 1) /= nu .OR. SIZE(band%characters, 2) /= SIZE(band%multiplicity)) RETURN
            IF (ANY(band%multiplicity < 0) .OR. ANY(band%product < 1) .OR. ANY(band%product > n)) RETURN
            IF (.NOT. ALL(ieee_is_finite(band%translations))) RETURN
            IF (.NOT. ALL(ieee_is_finite(REAL(band%factors, dp)))) RETURN
            IF (.NOT. ALL(ieee_is_finite(AIMAG(band%factors)))) RETURN
            IF (.NOT. ALL(ieee_is_finite(REAL(band%characters, dp)))) RETURN
            IF (.NOT. ALL(ieee_is_finite(AIMAG(band%characters)))) RETURN
            ALLOCATE (order(n), phase(n), inverse_unitary(n), old_character(nu), CHARACTER(nu))
            order = 0
            DO g = 1, n
               DO h = 1, n
                  IF (reference%antiunitary(g) .NEQV. band%antiunitary(h)) CYCLE
                  IF (ANY(reference%rotations(:, :, g) /= band%rotations(:, :, h))) CYCLE
                  delta = reference%translations(:, g) - band%translations(:, h)
                  IF (MAXVAL(ABS(delta)) > tolerance/(32*EPSILON(1.0_dp))) RETURN
                  IF (MAXVAL(ABS(delta - ANINT(delta))) > tolerance) CYCLE
                  IF (order(g) /= 0) RETURN
                  order(g) = h
                  angle = -twopi*DOT_PRODUCT(reference%k, ANINT(delta))
                  phase(g) = CMPLX(COS(angle), SIN(angle), dp)
               END DO
               IF (order(g) == 0) RETURN
            END DO
            DO g = 1, n
               IF (COUNT(order == g) /= 1) RETURN
            END DO
            status = -2
            DO h = 1, n
               DO g = 1, n
                  j = reference%product(g, h)
                  IF (order(j) /= band%product(order(g), order(h))) RETURN
                  second_phase = phase(h)
                  IF (reference%antiunitary(g)) second_phase = CONJG(second_phase)
                  expected = phase(g)*second_phase/phase(j)*band%factors(order(g), order(h))
                  IF (ABS(expected - reference%factors(g, h)) > tolerance) RETURN
               END DO
            END DO
            unitary = PACK([(g, g=1, n)],.NOT. band%antiunitary)
            inverse_unitary = 0
            DO g = 1, nu
               inverse_unitary(unitary(g)) = g
            END DO
            old_character(:) = MATMUL(band%characters, CMPLX(band%multiplicity, 0, dp))
            unitary = PACK([(g, g=1, n)],.NOT. reference%antiunitary)
            DO g = 1, nu
               h = unitary(g)
               CHARACTER(g) = phase(h)*old_character(inverse_unitary(order(h)))
            END DO
         END ASSOCIATE
         a = catalogue%row_begin(p)
         b = catalogue%row_begin(p + 1) - 1
         ALLOCATE (counts(b - a + 1))
         CALL atomic_signature_decompose(catalogue%points(p), catalogue%points(p)%mapping, CHARACTER, counts, tolerance, status)
         IF (status /= 0) RETURN
         work(a:b) = INT(counts, int_8)
         rank = SUM(INT(catalogue%row_dimension(a:b), int_8)*work(a:b))
         IF (first_rank < 0) first_rank = rank
         status = -2
         IF (rank /= first_rank) RETURN
         DEALLOCATE (order, phase, inverse_unitary, old_character, CHARACTER, counts)
      END DO
      TARGET = work
      status = 0
   END SUBROUTINE atomic_signature_from_bands
END MODULE topology_atomic_signatures
