!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
PROGRAM topology_site_unittest
   USE kinds,                           ONLY: dp
   USE topology_band_representations,   ONLY: atomic_band_characters,&
                                              atomic_band_construct,&
                                              atomic_band_site_data,&
                                              atomic_band_type
   USE topology_site_relations,         ONLY: site_induction_reduce,&
                                              site_reduction_type,&
                                              site_relation_type,&
                                              site_specialization
   USE topology_wyckoff,                ONLY: wyckoff_families,&
                                              wyckoff_family_type

   IMPLICIT NONE
   TYPE band_character_type
      INTEGER, ALLOCATABLE :: mapping(:)
      COMPLEX(dp), ALLOCATABLE :: characters(:, :)
   END TYPE band_character_type
   CALL analytic_tests()
   CALL reduction_tests()
   PRINT *, 'Site-specialization induction tests passed.'
CONTAINS

! **************************************************************************************************
!> \brief Transitive direct-sum witnesses, nonunique reductions and rejection of invalid integer data.
! **************************************************************************************************
   SUBROUTINE reduction_tests()
      INTEGER                                            :: status
      TYPE(site_reduction_type)                          :: reduced
      TYPE(site_relation_type)                           :: edges(3)

      edges(1)%source = 1
      edges(1)%destination = 2
      edges(1)%coefficients = RESHAPE([1], [1, 1])
      edges(2)%source = 1
      edges(2)%destination = 3
      edges(2)%coefficients = RESHAPE([1, 1], [2, 1])
      edges(3)%source = 3
      edges(3)%destination = 4
      edges(3)%coefficients = RESHAPE([1, 1], [1, 2])
      CALL site_induction_reduce(edges, [1, 2, 3, 5, 6], [2, 2, 1, 1, 1], &
                                 [.FALSE., .TRUE., .FALSE., .TRUE.], reduced, status)
      IF (status /= 0) ERROR STOP 'Synthetic transitive induction reduction failed'
      IF (ANY(reduced%representative /= [2, 2, 5, 5, 5])) ERROR STOP 'Maximal representatives wrong'
      IF (ANY(reduced%expansion(5, :) /= [2, 2, 1, 1, 1])) ERROR STOP 'Transitive direct sum wrong'
      IF (ANY(reduced%expansion(:4, :) /= 0)) ERROR STOP 'Composite generator survived expansion'
      IF (ANY(reduced%witness_edge(:2) /= 2)) ERROR STOP 'Transitive witness provenance lost'
      edges(2)%coefficients(1, 1) = -1
      CALL site_induction_reduce(edges, [1, 2, 3, 5, 6], [2, 2, 1, 1, 1], &
                                 [.FALSE., .TRUE., .FALSE., .TRUE.], reduced, status)
      IF (status == 0 .OR. ALLOCATED(reduced%expansion)) ERROR STOP 'Negative induction coefficient accepted'
      edges(2)%coefficients = HUGE(0)
      CALL site_induction_reduce(edges, [1, 2, 3, 5, 6], [2, 2, 1, 1, 1], &
                                 [.FALSE., .TRUE., .FALSE., .TRUE.], reduced, status)
      IF (status == 0 .OR. ALLOCATED(reduced%expansion)) ERROR STOP 'Invalid band dimension accepted'
      PRINT *, 'Transitive site-induction and malformed-witness tests passed.'
   END SUBROUTINE reduction_tests

! **************************************************************************************************
!> \brief Analytic inversion, mirror, shear/origin/coset conventions and magnetic half translations.
! **************************************************************************************************
   SUBROUTINE analytic_tests()
      INTEGER                                            :: composite, grey, i, ident(3, 3), &
                                                            inv(3, 3), n, r(3, 3, 4), relations, &
                                                            screw(3, 3, 8), shear(3, 3), spin, &
                                                            variant
      LOGICAL                                            :: anti(4)
      REAL(dp)                                           :: cell(3, 3), origin(3), shift(3, 8), &
                                                            tau(3, 4)

      ident = 0
      DO i = 1, 3
         ident(i, i) = 1
      END DO
      shear = RESHAPE([1, 0, 0, 1, 1, 0, 0, 1, 1], [3, 3])
      inv = RESHAPE([1, 0, 0, -1, 1, 0, 1, -1, 1], [3, 3])
      origin = [0.123_dp, -0.217_dp, 0.319_dp]
      DO variant = 1, 5
         r(:, :, 1) = ident
         r(:, :, 2) = -ident
         tau = 0.0_dp
         cell = REAL(ident, dp)
         IF (variant >= 2 .AND. variant <= 4) THEN
            r(:, :, 2) = ident
            r(3, 3, 2) = -1
         END IF
         IF (variant == 3) THEN
            cell = REAL(shear, dp)
            r(:, :, 2) = MATMUL(MATMUL(inv, r(:, :, 2)), shear)
         END IF
         IF (variant == 4) THEN
            tau(:, 2) = MATMUL(REAL(r(:, :, 2) - ident, dp), origin) + [1.0_dp, -2.0_dp, 3.0_dp]
         END IF
         IF (variant == 5) THEN
            r(:, :, 2) = ident
            tau(:, 2) = 0.5_dp
         END IF
         r(:, :, 3:4) = r(:, :, :2)
         tau(:, 3:4) = tau(:, :2)
         anti = [.FALSE., .FALSE., .TRUE., .TRUE.]
         DO grey = 1, 2
            n = 2*grey
            IF (variant == 5) THEN
               n = 2
               anti(2) = .TRUE.
            END IF
            DO spin = 0, 1
               CALL check_group(cell, r(:, :, :n), tau(:, :n), anti(:n), spin == 1, .TRUE., relations)
               IF (variant == 1 .AND. relations /= 8) ERROR STOP 'Inversion specialization count wrong'
               IF (variant >= 2 .AND. variant <= 4 .AND. relations /= 2) ERROR STOP 'Mirror specialization count wrong'
               IF (variant == 5 .AND. relations /= 0) ERROR STOP 'False magnetic specialization'
            END DO
            IF (variant == 5) EXIT
         END DO
      END DO
      ! P4_2/m has exceptional spinless S4 site representations, unlike its spinful counterpart.
      r(:, :, 1) = ident
      r(:, :, 2) = RESHAPE([0, 1, 0, -1, 0, 0, 0, 0, 1], [3, 3])
      r(:, :, 3) = MATMUL(r(:, :, 2), r(:, :, 2))
      r(:, :, 4) = MATMUL(r(:, :, 3), r(:, :, 2))
      shift = 0.0_dp
      DO i = 1, 4
         screw(:, :, 2*i - 1) = r(:, :, i)
         screw(:, :, 2*i) = -r(:, :, i)
         shift(3, 2*i - 1:2*i) = 0.5_dp*MODULO(i - 1, 2)
      END DO
      cell = REAL(ident, dp)
      DO spin = 0, 1
         CALL check_group(cell, RESHAPE([screw, screw], [3, 3, 16]), RESHAPE([shift, shift], [3, 16]), &
                          [(i > 8, i=1, 16)], &
                          spin == 1, .FALSE., relations, composite)
         IF (composite /= MERGE(2, 0, spin == 0)) ERROR STOP 'S4 maximal-site reduction wrong'
      END DO
      PRINT *, 'Site induction analytic, path, and convention checks passed.'
   END SUBROUTINE analytic_tests

! **************************************************************************************************
!> \brief Verify every geometric specialization against independent whole-band characters.
!> \param cell Invariant Cartesian cell
!> \param r Full direct-space rotations
!> \param tau Full fractional translations
!> \param anti Full time-reversal flags
!> \param spinful Physical spin representation
!> \param check_path Also independently construct intermediate real-space orbit positions
!> \param relations Number of tested specialization edges
!> \param composite Number of composite maximal-site representations
! **************************************************************************************************
   SUBROUTINE check_group(cell, r, tau, anti, spinful, check_path, relations, composite)
      REAL(dp), INTENT(IN)                               :: cell(3, 3)
      INTEGER, INTENT(IN)                                :: r(:, :, :)
      REAL(dp), INTENT(IN)                               :: tau(:, :)
      LOGICAL, INTENT(IN)                                :: anti(:), spinful, check_path
      INTEGER, INTENT(OUT)                               :: relations
      INTEGER, INTENT(OUT), OPTIONAL                     :: composite

      COMPLEX(dp), ALLOCATABLE                           :: image(:, :)
      INTEGER                                            :: a, b, e, j, nc, nr, p, status
      INTEGER, ALLOCATABLE                               :: dimensions(:), first(:)
      REAL(dp)                                           :: error, points(3, 8)
      TYPE(atomic_band_type), ALLOCATABLE                :: bands(:)
      TYPE(band_character_type), ALLOCATABLE             :: data(:, :)
      TYPE(site_reduction_type)                          :: reduction
      TYPE(site_relation_type)                           :: relation
      TYPE(site_relation_type), ALLOCATABLE              :: edges(:)
      TYPE(wyckoff_family_type), ALLOCATABLE             :: families(:)

      CALL wyckoff_families(cell, r, tau, anti, families, 1.e-8_dp, status)
      IF (status /= 0) ERROR STOP 'Site-induction families failed'
      points(:, 1) = 0.0_dp
      points(:, 2) = 0.5_dp
      points(:, 3) = [0.5_dp, 0.0_dp, 0.0_dp]
      points(:, 4) = [1.0_dp/3, 1.0_dp/3, 0.0_dp]
      points(:, 5) = [0.0_dp, 0.0_dp, 0.25_dp]
      points(:, 6) = [0.173_dp, 0.217_dp, 0.319_dp]
      points(:, 7) = [0.0_dp, 0.5_dp, 0.0_dp]
      points(:, 8) = [1.0_dp/3, 2.0_dp/3, 0.5_dp]
      ALLOCATE (bands(SIZE(families)), DATA(8, SIZE(families)), first(SIZE(families) + 1))
      first(1) = 1
      nr = 0
      DO a = 1, SIZE(families)
         CALL atomic_band_construct(cell, r, tau, anti, families(a)%position, spinful, bands(a), 1.e-8_dp, error, status)
         IF (status /= 0) ERROR STOP 'Site-induction band construction failed'
         DO p = 1, 8
            CALL atomic_band_characters(bands(a), points(:, p), DATA(p, a)%mapping, DATA(p, a)%characters, 1.e-8_dp, status)
            IF (status /= 0) ERROR STOP 'Site-induction band characters failed'
         END DO
         first(a + 1) = first(a) + SIZE(DATA(1, a)%characters, 2)
         nr = nr + SIZE(families(a)%specializations)
      END DO
      nc = first(SIZE(first)) - 1
      ALLOCATE (edges(nr), dimensions(nc))
      DO a = 1, SIZE(families)
         dimensions(first(a):first(a + 1) - 1) = NINT(REAL(DATA(1, a)%characters(1, :), dp))
      END DO
      relations = 0
      DO a = 1, SIZE(families)
         DO e = 1, SIZE(families(a)%specializations)
            b = families(a)%specializations(e)
            CALL site_specialization(cell, r, tau, anti, spinful, families, bands, a, b, relation, 1.e-8_dp, status)
            IF (status /= 0) THEN
               PRINT *, 'Site relation source/destination/status:', a, b, status
               PRINT *, 'Source and destination:', families(a)%position, families(b)%position
               ERROR STOP 'Site specialization failed'
            END IF
            DO p = 1, 8
               IF (ANY(DATA(p, a)%mapping /= DATA(p, b)%mapping)) ERROR STOP 'Global band operation mismatch'
               ALLOCATE (image(SIZE(DATA(p, a)%characters, 1), SIZE(DATA(p, a)%characters, 2)))
               image(:, :) = MATMUL(DATA(p, b)%characters, CMPLX(relation%coefficients, 0, dp))
               error = MAXVAL(ABS(image - DATA(p, a)%characters))
               IF (error > 1.e-7_dp) THEN
                  PRINT *, 'Site character source/destination/point/error:', a, b, p, error
                  ERROR STOP 'Induced character identity failed'
               END IF
               DEALLOCATE (image)
            END DO
            IF (check_path) CALL intermediate_test(cell, r, tau, anti, spinful, families(a)%position, &
                                                   relation%image, bands(a), points)
            relations = relations + 1
            edges(relations) = relation
         END DO
      END DO
      CALL site_induction_reduce(edges, first, dimensions, families%maximal_site_symmetry, reduction, status)
      IF (status /= 0) ERROR STOP 'Certified site-induction reduction failed'
      IF (PRESENT(composite)) composite = 0
      DO a = 1, SIZE(families)
         DO j = first(a), first(a + 1) - 1
            IF (PRESENT(composite)) THEN
               IF (families(a)%maximal_site_symmetry .AND. SUM(reduction%expansion(:, j)) > 1) composite = composite + 1
            END IF
            DO p = 1, 8
               ALLOCATE (image(SIZE(DATA(p, a)%characters, 1), 1))
               image = 0
               DO b = 1, SIZE(families)
                  image(:, 1) = image(:, 1) + MATMUL(DATA(p, b)%characters, &
                                                     CMPLX(reduction%expansion(first(b):first(b + 1) - 1, j), 0, dp))
               END DO
               IF (MAXVAL(ABS(image(:, 1) - DATA(p, a)%characters(:, j - first(a) + 1))) > 1.e-7_dp) THEN
                  ERROR STOP 'Transitive site-induction character identity failed'
               END IF
               DEALLOCATE (image)
            END DO
         END DO
      END DO
      ! Invalid and reversed edges must fail without returning a usable witness.
      CALL site_specialization(cell, r, tau, anti, spinful, families, bands, 0, 1, relation, 1.e-8_dp, status)
      IF (status == 0 .OR. ALLOCATED(relation%coefficients)) ERROR STOP 'Invalid site accepted'
      DO a = 1, SIZE(families)
         CALL site_specialization(cell, r, tau, anti, spinful, families, bands, a, a, relation, 1.e-8_dp, status)
         IF (status == 0 .OR. ALLOCATED(relation%coefficients)) ERROR STOP 'Self specialization accepted'
      END DO
   END SUBROUTINE check_group

! **************************************************************************************************
!> \brief Independently induce bands at intermediate points along the certified straight path.
!> \param cell Cartesian cell
!> \param r Full rotations
!> \param tau Full translations
!> \param anti Time-reversal flags
!> \param spinful Spin representation
!> \param source Generic source position
!> \param endpoint Lifted endpoint
!> \param band Source band representation
!> \param points Reciprocal sample points
! **************************************************************************************************
   SUBROUTINE intermediate_test(cell, r, tau, anti, spinful, source, endpoint, band, points)
      REAL(dp), INTENT(IN)                               :: cell(3, 3)
      INTEGER, INTENT(IN)                                :: r(:, :, :)
      REAL(dp), INTENT(IN)                               :: tau(:, :)
      LOGICAL, INTENT(IN)                                :: anti(:), spinful
      REAL(dp), INTENT(IN)                               :: source(3), endpoint(3)
      TYPE(atomic_band_type), INTENT(IN)                 :: band
      REAL(dp), INTENT(IN)                               :: points(:, :)

      COMPLEX(dp), ALLOCATABLE                           :: a(:, :), b(:, :)
      INTEGER                                            :: i, j, p, q, status
      INTEGER, ALLOCATABLE                               :: map_a(:), map_b(:), order(:)
      REAL(dp)                                           :: error, position(3)
      REAL(dp), ALLOCATABLE                              :: positions(:, :)
      TYPE(atomic_band_type)                             :: moved

      DO q = 1, 2
         position = source + REAL(q, dp)/3*(endpoint - source)
         CALL atomic_band_construct(cell, r, tau, anti, position, spinful, moved, 1.e-8_dp, error, status)
         IF (status /= 0) ERROR STOP 'Intermediate site construction failed'
         CALL atomic_band_site_data(band, positions, map_a, a)
         CALL atomic_band_site_data(moved, positions, map_b, b)
         IF (SIZE(map_a) /= SIZE(map_b)) ERROR STOP 'Intermediate path gained site symmetry'
         IF (ANY(map_a /= map_b)) ERROR STOP 'Intermediate site group changed'
         ALLOCATE (order(SIZE(a, 2)))
         DO i = 1, SIZE(a, 2)
            order(i) = 0
            DO j = 1, SIZE(b, 2)
               IF (MAXVAL(ABS(a(:, i) - b(:, j))) > 1.e-7_dp) CYCLE
               IF (order(i) /= 0) ERROR STOP 'Ambiguous intermediate irrep'
               order(i) = j
            END DO
            IF (order(i) == 0) ERROR STOP 'Intermediate irrep missing'
         END DO
         DO p = 1, SIZE(points, 2)
            CALL atomic_band_characters(band, points(:, p), map_a, a, 1.e-8_dp, status)
            IF (status /= 0) ERROR STOP 'Source path characters failed'
            CALL atomic_band_characters(moved, points(:, p), map_b, b, 1.e-8_dp, status)
            IF (status /= 0) ERROR STOP 'Intermediate path characters failed'
            IF (ANY(map_a /= map_b)) ERROR STOP 'Intermediate band operation mismatch'
            IF (MAXVAL(ABS(a - b(:, order))) > 1.e-7_dp) ERROR STOP 'Intermediate band characters changed'
         END DO
         DEALLOCATE (order)
      END DO
   END SUBROUTINE intermediate_test

END PROGRAM topology_site_unittest
