!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
PROGRAM topology_atomic_unittest
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE mathlib,                         ONLY: inv_3x3
   USE topology_atomic_signatures,      ONLY: atomic_signature_construct,&
                                              atomic_signature_decompose,&
                                              atomic_signature_from_bands,&
                                              atomic_signature_type
   USE topology_band_representations,   ONLY: atomic_band_site_snapshot,&
                                              atomic_band_type,&
                                              atomic_site_data_type,&
                                              atomic_site_storage
   USE topology_compatibility,          ONLY: little_group_compatibility,&
                                              little_group_data_type,&
                                              little_group_edge_type
   USE topology_integer_lattice,        ONLY: integer_lattice_certificate,&
                                              integer_lattice_factor,&
                                              integer_lattice_solve,&
                                              integer_lattice_type,&
                                              nonnegative_integer_solve
   USE topology_little_group,           ONLY: little_group_factor_system

   IMPLICIT NONE
   CALL analytic_tests()
   PRINT *, 'Atomic signature-matrix tests passed.'
CONTAINS

! **************************************************************************************************
!> \brief All inversion-center signatures, generic orbits, antiunitary pairing and invalid mappings.
! **************************************************************************************************
   SUBROUTINE analytic_tests()
      COMPLEX(dp)                                        :: chi(2)
      COMPLEX(dp), ALLOCATABLE                           :: image(:, :), reference(:, :)
      INTEGER                                            :: a, b, i, info, j, n, p, parity, perm(4), &
                                                            r(3, 3, 4), rr(3, 3, 4), s, site, spin
      INTEGER, ALLOCATABLE                               :: counts(:), reverse(:)
      LOGICAL                                            :: anti(4), used(17)
      REAL(dp)                                           :: angle, cell(3, 3), delta(3), &
                                                            points(3, 8), shifted(3, 4), tau(3, 4)
      TYPE(atomic_band_type)                             :: empty
      TYPE(atomic_signature_type)                        :: catalogue, changed
      TYPE(atomic_site_data_type)                        :: snapshot

      CALL atomic_band_site_snapshot(empty, snapshot)
      IF (ALLOCATED(snapshot%mapping) .OR. atomic_site_storage(snapshot) /= 0.0_dp) THEN
         ERROR STOP 'Empty onsite snapshot is not empty'
      END IF
      cell = 0.0_dp
      r = 0
      DO i = 1, 3
         cell(i, i) = 1.0_dp
         r(i, i, 1) = 1
      END DO
      r(:, :, 2) = -r(:, :, 1)
      r(:, :, 3:4) = r(:, :, 1:2)
      tau = 0.0_dp
      anti = [.FALSE., .FALSE., .TRUE., .TRUE.]
      DO p = 1, 8
         DO i = 1, 3
            points(i, p) = MERGE(0.5_dp, 0.0_dp, BTEST(p - 1, i - 1))
         END DO
      END DO
      DO spin = 0, 1
         n = 2 + 2*spin
         CALL atomic_signature_construct(cell, r(:, :, :n), tau(:, :n), anti(:n), points, spin == 1, &
                                         catalogue, 1.e-8_dp, info)
         IF (info /= 0) ERROR STOP 'Inversion signature construction failed'
         CALL check_onsite(catalogue, cell, r(:, :, :n), tau(:, :n), anti(:n))
         IF (SIZE(catalogue%matrix, 1) /= 16 .OR. SIZE(catalogue%matrix, 2) /= 17) THEN
            ERROR STOP 'Inversion signature matrix dimensions wrong'
         END IF
         IF (ANY(catalogue%character_class /= [(j, j=1, 17)])) ERROR STOP 'Distinct inversion character classes merged'
         CALL atomic_signature_construct(cell, r(:, :, :n), tau(:, :n), anti(:n), points(:, :1), spin == 1, &
                                         changed, 1.e-8_dp, info)
         IF (info /= 0) ERROR STOP 'Gamma-only atomic character classes failed'
         IF (ANY(changed%character_class /= [(j, j=1, 17)])) ERROR STOP 'Gamma aliases survive full-torus comparison'
         ALLOCATE (counts(2))
         DO j = 1, 17
            site = catalogue%column_site(j)
            s = 1 + spin
            IF (catalogue%families(site)%dimension == 3) THEN
               IF (catalogue%band_dimension(j) /= 2*s) ERROR STOP 'Generic inversion rank wrong'
               IF (ANY(catalogue%matrix(:, j) /= 1)) ERROR STOP 'Generic orbit must have both parities'
               CYCLE
            END IF
            IF (catalogue%band_dimension(j) /= s) ERROR STOP 'Centered inversion rank wrong'
            DO parity = -1, 1, 2
               chi = CMPLX([s, parity*s], 0, dp)
               CALL atomic_signature_decompose(catalogue%points(1), [1, 2], chi, counts, 1.e-8_dp, info)
               IF (info /= 0) ERROR STOP 'Analytic Gamma parity decomposition failed'
               IF (ALL(counts == catalogue%matrix(:2, j))) EXIT
            END DO
            IF (parity > 1) ERROR STOP 'Centered local parity not found'
            DO p = 1, 8
               delta = 2*catalogue%families(site)%position
               chi(1) = s
               chi(2) = parity*s*(-1.0_dp)**NINT(2*DOT_PRODUCT(points(:, p), delta))
               CALL atomic_signature_decompose(catalogue%points(p), [1, 2], chi, counts, 1.e-8_dp, info)
               IF (info /= 0) ERROR STOP 'Analytic TRIM parity decomposition failed'
               a = catalogue%row_begin(p)
               b = catalogue%row_begin(p + 1) - 1
               IF (ANY(counts /= catalogue%matrix(a:b, j))) ERROR STOP 'Inversion-center signature wrong'
            END DO
         END DO
         chi = CMPLX([1, 1], 0, dp)
         CALL atomic_signature_decompose(catalogue%points(1), [1, 1], chi, counts, 1.e-8_dp, info)
         IF (info == 0 .OR. ANY(counts /= -1)) ERROR STOP 'Duplicate character mapping accepted'
         CALL atomic_signature_decompose(catalogue%points(1), [1, 7], chi, counts, 1.e-8_dp, info)
         IF (info == 0 .OR. ANY(counts /= -1)) ERROR STOP 'Missing character mapping accepted'
         IF (spin == 1) THEN
            CALL atomic_signature_decompose(catalogue%points(1), [1, 2], chi, counts, 1.e-8_dp, info)
            IF (info == 0 .OR. ANY(counts /= -1)) ERROR STOP 'Odd Kramers multiplicity accepted'
         END IF
         DEALLOCATE (counts)
         CALL check_connections(catalogue)
         CALL check_integer_signatures(catalogue)
         CALL check_band_mapping(catalogue)

         ! Integer coset changes require row rephasing, not changes of column identity.
         perm = [2, 1, 4, 3]
         shifted = tau
         shifted(:, 2) = [1.0_dp, -2.0_dp, 3.0_dp]
         shifted(:, 4) = [-2.0_dp, 3.0_dp, 1.0_dp]
         rr = r(:, :, perm)
         CALL atomic_signature_construct(cell, rr(:, :, :n), shifted(:, perm(:n)), anti(perm(:n)), points, spin == 1, &
                                         changed, 1.e-8_dp, info)
         IF (info /= 0) ERROR STOP 'Permuted/rephased atomic signatures failed'
         CALL check_onsite(changed, cell, rr(:, :, :n), shifted(:, perm(:n)), anti(perm(:n)))
         used = .FALSE.
         DO j = 1, 17
            DO i = 1, 17
               IF (used(i)) CYCLE
               DO p = 1, 8
                  CALL physical_characters(catalogue, p, reference)
                  CALL physical_characters(changed, p, image)
                  reverse = [(FINDLOC(perm, a, DIM=1), a=1, 2)]
                  DO a = 1, 2
                     b = FINDLOC(changed%points(p)%mapping, reverse(a), DIM=1)
                     angle = -2*ACOS(-1.0_dp)*DOT_PRODUCT(points(:, p), shifted(:, a))
                     IF (ABS(image(b, i) - CMPLX(COS(angle), SIN(angle), dp)*reference(a, j)) > 1.e-7_dp) EXIT
                  END DO
                  IF (a <= 2) EXIT
               END DO
               IF (p <= 8) CYCLE
               used(i) = .TRUE.
               EXIT
            END DO
            IF (i > 17) ERROR STOP 'Atomic signature column lost under coset change'
         END DO
      END DO
      CALL atomic_signature_construct(cell, r(:, :, :2), tau(:, :2), anti(:2), points(:, :0), .FALSE., &
                                      catalogue, 1.e-8_dp, info)
      IF (info == 0 .OR. ALLOCATED(catalogue%matrix)) ERROR STOP 'Empty reciprocal set accepted'
      CALL atomic_signature_construct(cell, r(:, :, :2), tau(:, :2), anti(:2), points, .FALSE., &
                                      catalogue, 1.e-8_dp, info, 0.0001_dp)
      IF (info /= -3 .OR. ALLOCATED(catalogue%matrix)) ERROR STOP 'Atomic matrix memory budget ignored'
      PRINT *, 'Atomic signature analytic, gauge, and rejection tests passed.'
   END SUBROUTINE analytic_tests

! **************************************************************************************************
!> \brief Onsite snapshots retain full group membership, lattice shifts and physical spin factors.
!> \param catalogue Atomic reference catalogue
!> \param cell Input Cartesian cell
!> \param rotations Original direct-space rotations
!> \param translations Original fractional representatives
!> \param anti Original antiunitary flags
! **************************************************************************************************
   SUBROUTINE check_onsite(catalogue, cell, rotations, translations, anti)
      TYPE(atomic_signature_type), INTENT(IN)            :: catalogue
      REAL(dp), INTENT(IN)                               :: cell(3, 3)
      INTEGER, INTENT(IN)                                :: rotations(:, :, :)
      REAL(dp), INTENT(IN)                               :: translations(:, :)
      LOGICAL, INTENT(IN)                                :: anti(:)

      COMPLEX(dp)                                        :: composed(2, 2)
      INTEGER                                            :: g, h, i, j, k, n, p, s
      REAL(dp)                                           :: cart(3, 3), delta(3), seed(3)

      IF (SIZE(catalogue%sites) /= SIZE(catalogue%families)) ERROR STOP 'Missing onsite snapshots'
      DO s = 1, SIZE(catalogue%sites)
         ASSOCIATE (site => catalogue%sites(s))
            n = SIZE(site%mapping)
            IF (atomic_site_storage(site) <= 0.0_dp) ERROR STOP 'Onsite memory estimate is empty'
            IF (SIZE(site%characters, 1) /= COUNT(.NOT. site%antiunitary)) ERROR STOP 'Onsite character rows wrong'
            IF (SIZE(site%characters, 2) /= COUNT(catalogue%column_site == s)) ERROR STOP 'Onsite columns lost'
            seed = MODULO(catalogue%families(s)%position, 1.0_dp)
            DO g = 1, SIZE(anti)
               delta = MATMUL(REAL(rotations(:, :, g), dp), seed) + translations(:, g) - seed
               i = COUNT(site%mapping == g)
               IF (i /= MERGE(1, 0, MAXVAL(ABS(delta - ANINT(delta))) < 1.e-8_dp)) THEN
                  ERROR STOP 'Onsite snapshot loses stabilizer membership'
               END IF
            END DO
            DO i = 1, n
               g = site%mapping(i)
               IF (MAXVAL(ABS(site%spin(:, :, i) - catalogue%reference_spin(:, :, g))) > 1.e-12_dp) THEN
                  ERROR STOP 'Onsite and full-group lifts disagree'
               END IF
               IF (site%antiunitary(i) .NEQV. anti(g)) ERROR STOP 'Onsite antiunitary flag lost'
               delta = MATMUL(REAL(rotations(:, :, g), dp), seed) + translations(:, g) - seed
               IF (MAXVAL(ABS(delta - site%lattice_shifts(:, i))) > 1.e-8_dp) ERROR STOP 'Onsite lattice lift lost'
               cart = MATMUL(MATMUL(cell, REAL(rotations(:, :, g), dp)), inv_3x3(cell))
               IF (MAXVAL(ABS(cart - site%cartesian(:, :, i))) > 1.e-8_dp) ERROR STOP 'Onsite Cartesian axes lost'
               DO j = 1, n
                  h = site%mapping(j)
                  p = site%product(i, j)
                  k = site%mapping(p)
                  IF (ANY(MATMUL(rotations(:, :, g), rotations(:, :, h)) /= rotations(:, :, k))) THEN
                     ERROR STOP 'Onsite product operation lost'
                  END IF
                  composed = site%spin(:, :, j)
                  IF (site%antiunitary(i)) composed = CONJG(composed)
                  composed = MATMUL(site%spin(:, :, i), composed)
                  IF (.NOT. site%spinful) THEN
                     IF (ABS(site%factor(i, j) - 1.0_dp) > 1.e-8_dp) ERROR STOP 'Scalar onsite factor is not trivial'
                  ELSE IF (MAXVAL(ABS(composed - site%factor(i, j)*site%spin(:, :, p))) > 1.e-8_dp) THEN
                     ERROR STOP 'Onsite physical spin convention lost'
                  END IF
               END DO
            END DO
            DO j = 1, SIZE(catalogue%column_site)
               IF (catalogue%column_site(j) /= s) CYCLE
               i = catalogue%column_local(j)
               IF (ABS(site%characters(1, i)*catalogue%families(s)%multiplicity - catalogue%band_dimension(j)) > 1.e-8_dp) THEN
                  ERROR STOP 'Onsite snapshot does not describe its induced band column'
               END IF
            END DO
         END ASSOCIATE
      END DO
   END SUBROUTINE check_onsite

! **************************************************************************************************
!> \brief Recover unitary character vectors from signature rows, retaining corepresentation weights.
!> \param catalogue Generated catalogue
!> \param p Point index
!> \param characters Character columns in original unitary operation order
! **************************************************************************************************
   SUBROUTINE physical_characters(catalogue, p, characters)
      TYPE(atomic_signature_type), INTENT(IN)            :: catalogue
      INTEGER, INTENT(IN)                                :: p
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: characters(:, :)

      INTEGER, ALLOCATABLE                               :: counts(:, :)

      CALL unitary_counts(catalogue, p, counts)
      characters = MATMUL(catalogue%points(p)%little%characters, CMPLX(counts, 0, dp))
   END SUBROUTINE physical_characters

! **************************************************************************************************
!> \brief Expand signature corepresentation rows to unitary multiplicities for line subduction.
!> \param catalogue Catalogue
!> \param p Point index
!> \param counts Unit-only rows, identical columns to the catalogue
! **************************************************************************************************
   SUBROUTINE unitary_counts(catalogue, p, counts)
      TYPE(atomic_signature_type), INTENT(IN)            :: catalogue
      INTEGER, INTENT(IN)                                :: p
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: counts(:, :)

      INTEGER                                            :: a, b

      a = catalogue%row_begin(p)
      b = catalogue%row_begin(p + 1) - 1
      IF (ALLOCATED(catalogue%points(p)%coreps%restriction)) THEN
         counts = MATMUL(catalogue%points(p)%coreps%restriction, catalogue%matrix(a:b, :))
      ELSE
         counts = catalogue%matrix(a:b, :)
      END IF
   END SUBROUTINE unitary_counts

! **************************************************************************************************
!> \brief Every atomic column must satisfy directed segment restrictions, including a reciprocal loop.
!> \param catalogue Constructed catalogue
! **************************************************************************************************
   SUBROUTINE check_connections(catalogue)
      TYPE(atomic_signature_type), INTENT(IN)            :: catalogue

      COMPLEX(dp), ALLOCATABLE                           :: characters(:, :), chi(:)
      INTEGER                                            :: a, b, info, np, p, q, shift(3)
      INTEGER, ALLOCATABLE                               :: counts(:), expected(:), left(:, :), &
                                                            order(:), right(:, :), weights(:)
      REAL(dp)                                           :: residual
      TYPE(little_group_edge_type)                       :: edge

      np = SIZE(catalogue%points)
      CALL unitary_counts(catalogue, 1, left)
      DO p = 1, np + 1
         q = MODULO(p - 1, np) + 1
         shift = 0
         IF (p == np + 1) shift(3) = 1
         CALL little_group_compatibility(catalogue%points(1)%little, catalogue%points(q)%little, shift, edge, &
                                         1.e-7_dp, residual, info)
         IF (info /= 0) ERROR STOP 'Atomic catalogue line table failed'
         CALL unitary_counts(catalogue, q, right)
         IF (ANY(MATMUL(edge%left_branch, left) /= MATMUL(edge%right_branch, right))) THEN
            ERROR STOP 'Atomic signature column violates compatibility'
         END IF
      END DO
      weights = [(MODULO(p, 3), p=1, SIZE(catalogue%matrix, 2))]
      DO p = 1, np
         a = catalogue%row_begin(p)
         b = catalogue%row_begin(p + 1) - 1
         CALL physical_characters(catalogue, p, characters)
         ALLOCATE (chi(SIZE(characters, 1)), expected(b - a + 1), counts(b - a + 1))
         chi(:) = MATMUL(characters, CMPLX(weights, 0, dp))
         expected(:) = MATMUL(catalogue%matrix(a:b, :), weights)
         counts(:) = expected
         order = [(q, q=SIZE(chi), 1, -1)]
         CALL atomic_signature_decompose(catalogue%points(p), catalogue%points(p)%mapping(order), chi(order), &
                                         counts, 1.e-7_dp, info)
         IF (info /= 0 .OR. ANY(counts /= expected)) THEN
            PRINT *, 'Signature reconstruction point/status/expected/counts:', p, info, expected, counts
            ERROR STOP 'Atomic signature additivity/row mapping failed'
         END IF
         DEALLOCATE (chi, expected, counts)
      END DO
   END SUBROUTINE check_connections

! **************************************************************************************************
!> \brief Factor generated signature matrices and recover signed and nonnegative atomic witnesses.
!> \param catalogue Constructed catalogue, all retained reciprocal rows and onsite columns
! **************************************************************************************************
   SUBROUTINE check_integer_signatures(catalogue)
      TYPE(atomic_signature_type), INTENT(IN)            :: catalogue

      INTEGER                                            :: info, j, m, n, rank
      INTEGER(int_8)                                     :: visited
      INTEGER(int_8), ALLOCATABLE                        :: diagonal(:), left(:, :), matrix(:, :), &
                                                            right(:, :), TARGET(:), &
                                                            verified(:, :), weights(:), witness(:)
      LOGICAL                                            :: exists
      TYPE(integer_lattice_type)                         :: factor

      m = SIZE(catalogue%matrix, 1)
      n = SIZE(catalogue%matrix, 2)
      ALLOCATE (matrix(m, n), TARGET(m), weights(n), verified(m, n))
      matrix = INT(catalogue%matrix, int_8)
      CALL integer_lattice_factor(matrix, factor, info)
      IF (info /= 0) ERROR STOP 'Generated atomic Smith form failed'
      CALL integer_lattice_certificate(factor, left, right, diagonal, rank)
      verified(:, :) = MATMUL(MATMUL(left, matrix), right)
      DO j = 1, MIN(m, n)
         verified(j, j) = verified(j, j) - diagonal(j)
      END DO
      IF (ANY(verified /= 0) .OR. rank < 1) ERROR STOP 'Generated Smith certificate failed'
      weights = [(INT(MODULO(j, 5) - 2, int_8), j=1, n)]
      TARGET(:) = MATMUL(matrix, weights)
      CALL integer_lattice_solve(factor, TARGET, exists, witness, info)
      IF (info /= 0 .OR. .NOT. exists) ERROR STOP 'Generated signed signature not recovered'
      DO j = 1, n + 1
         IF (j <= n) THEN
            TARGET = matrix(:, j)
         ELSE
            TARGET = matrix(:, 1) + matrix(:, n) + matrix(:, (n + 1)/2)
         END IF
         CALL nonnegative_integer_solve(matrix, TARGET, exists, witness, visited, info, 1000000_int_8)
         IF (info /= 0 .OR. .NOT. exists) THEN
            PRINT *, 'Atomic integer dimensions/column/status/nodes:', m, n, j, info, visited
            ERROR STOP 'Generated nonnegative signature not recovered'
         END IF
      END DO
   END SUBROUTINE check_integer_signatures

! **************************************************************************************************
!> \brief Match band data with independently regenerated factors after coset and operation/irrep changes.
!> \param catalogue Constructed reference
! **************************************************************************************************
   SUBROUTINE check_band_mapping(catalogue)
      TYPE(atomic_signature_type), INTENT(IN)            :: catalogue

      INTEGER                                            :: g, h, i, info, n, nu, p
      INTEGER(int_8), ALLOCATABLE                        :: expected(:), TARGET(:)
      INTEGER, ALLOCATABLE                               :: counts(:, :), inverse(:), labels(:), &
                                                            mapping(:), permutation(:), &
                                                            unitary(:), weights(:)
      REAL(dp)                                           :: angle
      REAL(dp), ALLOCATABLE                              :: shifts(:, :), translations(:, :)
      TYPE(little_group_data_type)                       :: original
      TYPE(little_group_data_type), ALLOCATABLE          :: points(:)

      ALLOCATE (points(SIZE(catalogue%points)), weights(SIZE(catalogue%matrix, 2)), expected(SIZE(catalogue%matrix, 1)))
      weights = [(MODULO(i, 3), i=1, SIZE(weights))]
      expected(:) = MATMUL(INT(catalogue%matrix, int_8), INT(weights, int_8))
      DO p = 1, SIZE(points)
         points(p) = catalogue%points(p)%little
         CALL unitary_counts(catalogue, p, counts)
         points(p)%multiplicity(:) = MATMUL(counts, weights)
      END DO
      CALL atomic_signature_from_bands(catalogue, points, TARGET, 1.e-7_dp, info)
      IF (info /= 0) ERROR STOP 'Identity band-to-atomic mapping failed'
      IF (ANY(TARGET /= expected)) ERROR STOP 'Identity band signature incorrect'
      DO p = 1, SIZE(points)
         original = points(p)
         n = SIZE(original%antiunitary)
         nu = SIZE(original%characters, 1)
         ALLOCATE (translations(3, n), shifts(3, n), inverse(n))
         DO g = 1, n
            shifts(:, g) = REAL([MODULO(g, 3), -MODULO(g, 2), MODULO(g, 5)], dp)
         END DO
         shifts(:, 1) = 0.0_dp
         translations = original%translations + shifts
         permutation = [(i, i=n, 1, -1)]
         CALL little_group_factor_system(original%cell, original%rotations(:, :, permutation), &
                                         translations(:, permutation), original%antiunitary(permutation), original%k, &
                                         original%spinful, mapping, points(p)%product, points(p)%factors, 1.e-8_dp, info)
         IF (info /= 0) ERROR STOP 'Independent coset-shifted factor construction failed'
         permutation = permutation(mapping)
         points(p)%rotations = original%rotations(:, :, permutation)
         points(p)%translations = translations(:, permutation)
         points(p)%antiunitary = original%antiunitary(permutation)
         unitary = PACK([(g, g=1, n)],.NOT. original%antiunitary)
         inverse = 0
         DO i = 1, nu
            inverse(unitary(i)) = i
         END DO
         labels = [(i, i=SIZE(original%multiplicity), 1, -1)]
         points(p)%multiplicity = original%multiplicity(labels)
         unitary = PACK([(g, g=1, n)],.NOT. points(p)%antiunitary)
         DO i = 1, nu
            g = unitary(i)
            h = permutation(g)
            angle = -2*ACOS(-1.0_dp)*DOT_PRODUCT(original%k, shifts(:, h))
            points(p)%characters(i, :) = CMPLX(COS(angle), SIN(angle), dp)*original%characters(inverse(h), labels)
         END DO
         DEALLOCATE (translations, shifts, inverse)
      END DO
      CALL atomic_signature_from_bands(catalogue, points, TARGET, 1.e-7_dp, info)
      IF (info /= 0) THEN
         PRINT *, 'Rephased band mapping status:', info
         ERROR STOP 'Coset-shifted band mapping failed'
      END IF
      IF (ANY(TARGET /= expected)) ERROR STOP 'Rephased band signature incorrect'
      points(1)%factors(1, 1) = -points(1)%factors(1, 1)
      CALL atomic_signature_from_bands(catalogue, points, TARGET, 1.e-7_dp, info)
      IF (info /= -2 .OR. ALLOCATED(TARGET)) ERROR STOP 'Inconsistent band factors accepted'
      points(1)%factors(1, 1) = -points(1)%factors(1, 1)
      points(1)%translations(:, 1) = points(1)%translations(:, 1) + 0.031_dp
      CALL atomic_signature_from_bands(catalogue, points, TARGET, 1.e-7_dp, info)
      IF (info == 0 .OR. ALLOCATED(TARGET)) ERROR STOP 'Unmatched band operation accepted'
   END SUBROUTINE check_band_mapping

END PROGRAM topology_atomic_unittest
