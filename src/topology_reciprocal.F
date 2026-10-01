!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Reciprocal fixed strata, incidence embeddings and primitive torus cycles.
!> \author Thomas D. Kuehne, tkuehne@cp2k.org
! **************************************************************************************************
MODULE topology_reciprocal
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE mathlib,                         ONLY: inv_3x3
   USE topology_integer_lattice,        ONLY: integer_lattice_certificate,&
                                              integer_lattice_factor,&
                                              integer_lattice_solve,&
                                              integer_lattice_type
   USE topology_little_group,           ONLY: little_group_factor_system
   USE topology_wyckoff,                ONLY: wyckoff_families,&
                                              wyckoff_family_type

   IMPLICIT NONE
   PRIVATE

   TYPE, PUBLIC :: reciprocal_stratum_type
      INTEGER :: dimension = 3, star_size = 0
      INTEGER :: normal(3, 3) = 0
      REAL(dp) :: offset(3) = 0.0_dp, k(3) = 0.0_dp
      INTEGER, ALLOCATABLE :: stabilizer(:), cycles(:, :), specializations(:)
      INTEGER, ALLOCATABLE :: image_operation(:), image_shift(:, :)
   END TYPE reciprocal_stratum_type

   TYPE, PUBLIC :: reciprocal_sampling_type
      TYPE(reciprocal_stratum_type), ALLOCATABLE :: strata(:)
      REAL(dp), ALLOCATABLE :: kpoints(:, :)
      INTEGER, ALLOCATABLE :: source(:), destination(:), shift(:, :)
   END TYPE reciprocal_sampling_type

   PUBLIC :: reciprocal_strata, reciprocal_sampling, reciprocal_refine
CONTAINS

! **************************************************************************************************
!> \brief Insert uniformly spaced interior samples, preserving every lifted directed segment.
!> \param sampling Valid output of reciprocal_sampling; unchanged on failure
!> \param intervals Number of subdivisions per segment, at least one
!> \param max_points Conservative point/work budget before point deduplication
!> \param tolerance Absolute fractional-coordinate matching tolerance
!> \param status Zero on success; -1 invalid input; -3 work budget exceeded
!> \note The shift of a subedge is the difference of the two endpoints' reciprocal-cell offsets.
!>       Discarding the source offset would change a path that crosses a reciprocal boundary.
! **************************************************************************************************
   SUBROUTINE reciprocal_refine(sampling, intervals, max_points, tolerance, status)
      TYPE(reciprocal_sampling_type), INTENT(INOUT)      :: sampling
      INTEGER, INTENT(IN)                                :: intervals, max_points
      REAL(dp), INTENT(IN)                               :: tolerance
      INTEGER, INTENT(OUT)                               :: status

      INTEGER                                            :: capacity, e, edge, j, ne, new_shift(3), &
                                                            np, old_shift(3), p, previous
      INTEGER, ALLOCATABLE                               :: destinations(:), shifts(:, :), sources(:)
      REAL(dp)                                           :: delta(3), direction(3), image(3), &
                                                            start(3), wrapped(3)
      REAL(dp), ALLOCATABLE                              :: points(:, :)

      status = -1
      IF (.NOT. ieee_is_finite(tolerance)) RETURN
      IF (intervals < 1 .OR. max_points < 1 .OR. tolerance <= 0.0_dp .OR. tolerance >= 0.01_dp) RETURN
      IF (.NOT. ALLOCATED(sampling%kpoints) .OR. .NOT. ALLOCATED(sampling%source)) RETURN
      IF (.NOT. ALLOCATED(sampling%destination) .OR. .NOT. ALLOCATED(sampling%shift)) RETURN
      np = SIZE(sampling%kpoints, 2)
      ne = SIZE(sampling%source)
      IF (np < 1 .OR. .NOT. ALL(ieee_is_finite(sampling%kpoints))) RETURN
      IF (ANY(sampling%kpoints < -tolerance) .OR. ANY(sampling%kpoints >= 1.0_dp + tolerance)) RETURN
      IF (SIZE(sampling%kpoints, 1) /= 3 .OR. SIZE(sampling%destination) /= ne) RETURN
      IF (ANY(SHAPE(sampling%shift) /= [3, ne])) RETURN
      IF (ANY(sampling%source < 1) .OR. ANY(sampling%source > np)) RETURN
      IF (ANY(sampling%destination < 1) .OR. ANY(sampling%destination > np)) RETURN
      status = -3
      IF (INT(np, int_8) + INT(ne, int_8)*(intervals - 1_int_8) > max_points) RETURN
      IF (INT(ne, int_8)*intervals > HUGE(0)) RETURN
      capacity = np + ne*(intervals - 1)
      ALLOCATE (points(3, capacity), sources(ne*intervals), destinations(ne*intervals), shifts(3, ne*intervals))
      points(:, :np) = sampling%kpoints
      edge = 0
      DO e = 1, ne
         previous = sampling%source(e)
         old_shift = 0
         start = sampling%kpoints(:, previous)
         direction = sampling%kpoints(:, sampling%destination(e)) + REAL(sampling%shift(:, e), dp) - start
         DO j = 1, intervals
            image = start + REAL(j, dp)/intervals*direction
            wrapped = MODULO(image, 1.0_dp)
            IF (j == intervals) THEN
               p = sampling%destination(e)
               new_shift = sampling%shift(:, e)
            ELSE
               DO p = 1, np
                  delta = points(:, p) - wrapped
                  IF (MAXVAL(ABS(delta - ANINT(delta))) <= tolerance) EXIT
               END DO
               IF (p > np) THEN
                  np = np + 1
                  points(:, np) = wrapped
               END IF
               IF (MAXVAL(ABS(image - points(:, p))) > REAL(HUGE(0)/2, dp)) RETURN
               new_shift = NINT(image - points(:, p))
            END IF
            IF (ANY(ABS(INT(new_shift, int_8) - INT(old_shift, int_8)) > HUGE(0))) RETURN
            edge = edge + 1
            sources(edge) = previous
            destinations(edge) = p
            shifts(:, edge) = new_shift - old_shift
            previous = p
            old_shift = new_shift
         END DO
      END DO
      sampling%kpoints = points(:, :np)
      CALL MOVE_ALLOC(sources, sampling%source)
      CALL MOVE_ALLOC(destinations, sampling%destination)
      CALL MOVE_ALLOC(shifts, sampling%shift)
      status = 0
   END SUBROUTINE reciprocal_refine

! **************************************************************************************************
!> \brief Generate explicit representative points and directed incidence/cycle segments.
!> \param cell Cartesian direct-lattice columns
!> \param rotations Complete integer direct rotations
!> \param translations Fractional coset translations
!> \param antiunitary Physical time-reversal flags
!> \param sampling Strata, wrapped points and directed edges with retained reciprocal shifts
!> \param tolerance Absolute geometry tolerance
!> \param status Zero on success, negative on failure; no partial output
!> \note The first points are the stratum seeds. Additional child images preserve incidence
!>       embeddings without requiring an irrep-label transport convention. Cycles are self-edges
!>       with nonzero integer shifts. These data can feed the explicit-segment compatibility code;
!>       they do not yet implement an open-cell complex or a global band-connectivity classifier.
! **************************************************************************************************
   SUBROUTINE reciprocal_sampling(cell, rotations, translations, antiunitary, sampling, tolerance, status)
      REAL(dp), INTENT(IN)                               :: cell(3, 3)
      INTEGER, INTENT(IN)                                :: rotations(:, :, :)
      REAL(dp), INTENT(IN)                               :: translations(:, :)
      LOGICAL, INTENT(IN)                                :: antiunitary(:)
      TYPE(reciprocal_sampling_type), INTENT(OUT)        :: sampling
      REAL(dp), INTENT(IN)                               :: tolerance
      INTEGER, INTENT(OUT)                               :: status

      INTEGER                                            :: edge, g, h, i, j, ne, nf, np, p, &
                                                            rotation(3, 3)
      REAL(dp)                                           :: delta(3), image(3), wrapped(3)
      REAL(dp), ALLOCATABLE                              :: points(:, :)
      TYPE(reciprocal_sampling_type)                     :: work

      CALL reciprocal_strata(cell, rotations, translations, antiunitary, work%strata, tolerance, status)
      IF (status /= 0) RETURN
      nf = SIZE(work%strata)
      ne = 0
      DO i = 1, nf
         ne = ne + SIZE(work%strata(i)%specializations) + work%strata(i)%dimension
      END DO
      ALLOCATE (points(3, nf + ne), work%source(ne), work%destination(ne), work%shift(3, ne))
      DO i = 1, nf
         points(:, i) = work%strata(i)%k
      END DO
      np = nf
      edge = 0
      DO i = 1, nf
         DO j = 1, SIZE(work%strata(i)%specializations)
            h = work%strata(i)%specializations(j)
            g = work%strata(i)%image_operation(j)
            CALL dual_rotation(rotations(:, :, g), rotation, status)
            IF (status /= 0) RETURN
            IF (antiunitary(g)) rotation = -rotation
            image = MATMUL(REAL(rotation, dp), work%strata(h)%k) + REAL(work%strata(i)%image_shift(:, j), dp)
            wrapped = MODULO(image, 1.0_dp)
            DO p = 1, np
               delta = points(:, p) - wrapped
               IF (MAXVAL(ABS(delta - ANINT(delta))) <= tolerance) EXIT
            END DO
            IF (p > np) THEN
               np = np + 1
               points(:, np) = wrapped
            END IF
            edge = edge + 1
            work%source(edge) = i
            work%destination(edge) = p
            status = -3
            IF (MAXVAL(ABS(image - points(:, p))) > REAL(HUGE(0), dp)) RETURN
            work%shift(:, edge) = NINT(image - points(:, p))
         END DO
         DO j = 1, work%strata(i)%dimension
            edge = edge + 1
            work%source(edge) = i
            work%destination(edge) = i
            work%shift(:, edge) = work%strata(i)%cycles(:, j)
         END DO
      END DO
      work%kpoints = points(:, :np)
      sampling = work
      status = 0
   END SUBROUTINE reciprocal_sampling

! **************************************************************************************************
!> \brief Generate connected fixed tori of the reciprocal action, modulo its point-group orbits.
!> \param cell Cartesian direct-lattice columns
!> \param rotations Complete integer direct-space rotation group modulo cell translations
!> \param translations Original fractional coset translations, retained for group validation
!> \param antiunitary Physical time-reversal flags
!> \param strata Generic seeds, full stabilizers, incidence embeddings and primitive cycles
!> \param tolerance Absolute geometry tolerance
!> \param status Zero on success; negative for malformed data, arithmetic limits or geometry failure
!> \note The effective geometric action is (-1)^a W^(-T). Identical actions are collapsed only
!>       for geometry; the returned stabilizers contain all original nonsymmorphic/centering cosets.
!>       N*k=offset modulo integers describes a connected closure; higher-symmetry subsets are
!>       excluded from its generic seed. Cycles span ker(N) over Z, not just over R.
!>       Each incidence gives an original operation and an integer shift taking a child seed into
!>       the parent's lifted affine closure. All higher-symmetry specializations are retained.
!>       This is not a Brillouin-zone polyhedron, a connected-cell decomposition of open strata,
!>       or a complete band-connectivity/EBR classification. Conventional-cell folding is retained.
!>       No partial result is returned on failure.
! **************************************************************************************************
   SUBROUTINE reciprocal_strata(cell, rotations, translations, antiunitary, strata, tolerance, status)
      REAL(dp), INTENT(IN)                               :: cell(3, 3)
      INTEGER, INTENT(IN)                                :: rotations(:, :, :)
      REAL(dp), INTENT(IN)                               :: translations(:, :)
      LOGICAL, INTENT(IN)                                :: antiunitary(:)
      TYPE(reciprocal_stratum_type), ALLOCATABLE, &
         INTENT(OUT)                                     :: strata(:)
      REAL(dp), INTENT(IN)                               :: tolerance
      INTEGER, INTENT(OUT)                               :: status

      COMPLEX(dp), ALLOCATABLE                           :: factors(:, :)
      INTEGER                                            :: c, g, h, i, j, n, nf, ng, rank
      INTEGER(int_8)                                     :: normal(3, 3), TARGET(3)
      INTEGER(int_8), ALLOCATABLE                        :: diagonal(:), left(:, :), right(:, :), &
                                                            shift(:)
      INTEGER, ALLOCATABLE                               :: dual(:, :, :), mapping(:), &
                                                            products(:, :), unique(:, :, :)
      LOGICAL                                            :: exists, found
      LOGICAL, ALLOCATABLE                               :: fixed(:)
      REAL(dp)                                           :: delta(3), image(3), reciprocal_cell(3, 3)
      TYPE(integer_lattice_type)                         :: integer_normal
      TYPE(reciprocal_stratum_type), ALLOCATABLE         :: work(:)
      TYPE(wyckoff_family_type), ALLOCATABLE             :: families(:)

      ! At Gamma every original coset participates, so this checks the entire supplied group.
      CALL little_group_factor_system(cell, rotations, translations, antiunitary, [0.0_dp, 0.0_dp, 0.0_dp], &
                                      .FALSE., mapping, products, factors, tolerance, status)
      IF (status /= 0) RETURN
      n = SIZE(mapping)
      ALLOCATE (dual(3, 3, n), unique(3, 3, n), fixed(n))
      ng = 0
      DO g = 1, n
         CALL dual_rotation(rotations(:, :, g), dual(:, :, g), status)
         IF (status /= 0) RETURN
         IF (antiunitary(g)) dual(:, :, g) = -dual(:, :, g)
         found = .FALSE.
         DO h = 1, ng
            IF (ALL(unique(:, :, h) == dual(:, :, g))) found = .TRUE.
         END DO
         IF (found) CYCLE
         ng = ng + 1
         unique(:, :, ng) = dual(:, :, g)
      END DO
      reciprocal_cell = TRANSPOSE(inv_3x3(cell/MAXVAL(ABS(cell))))
      CALL wyckoff_families(reciprocal_cell, unique(:, :, :ng), SPREAD([0.0_dp, 0.0_dp, 0.0_dp], 2, ng), &
                            [(.FALSE., g=1, ng)], families, tolerance, status)
      IF (status /= 0) RETURN
      nf = SIZE(families)
      ALLOCATE (work(nf))
      DO i = 1, nf
         work(i)%dimension = families(i)%dimension
         work(i)%star_size = families(i)%multiplicity
         work(i)%normal = families(i)%normal
         work(i)%offset = families(i)%offset
         work(i)%k = families(i)%position
         work(i)%specializations = families(i)%specializations
         DO g = 1, n
            delta = MATMUL(REAL(dual(:, :, g), dp), work(i)%k) - work(i)%k
            fixed(g) = MAXVAL(ABS(delta - ANINT(delta))) <= tolerance
         END DO
         status = -2
         IF (COUNT(fixed)*work(i)%star_size /= n) RETURN
         work(i)%stabilizer = PACK([(g, g=1, n)], fixed)
         normal = INT(work(i)%normal, int_8)
         CALL integer_lattice_factor(normal, integer_normal, status)
         IF (status /= 0) RETURN
         CALL integer_lattice_certificate(integer_normal, left, right, diagonal, rank)
         status = -2
         IF (rank /= 3 - work(i)%dimension .OR. ANY(diagonal(:rank) /= 1)) RETURN
         status = -3
         IF (ANY(ABS(right) > 100000_int_8)) RETURN
         work(i)%cycles = INT(right(:, rank + 1:))
      END DO

      DO i = 1, nf
         normal = INT(work(i)%normal, int_8)
         CALL integer_lattice_factor(normal, integer_normal, status)
         IF (status /= 0) RETURN
         c = SIZE(work(i)%specializations)
         ALLOCATE (work(i)%image_operation(c), work(i)%image_shift(3, c))
         DO j = 1, c
            h = work(i)%specializations(j)
            found = .FALSE.
            DO g = 1, n
               ! Both the child seed and every tangent direction must lie in the parent closure.
               IF (ANY(MATMUL(MATMUL(normal, INT(dual(:, :, g), int_8)), &
                              INT(work(h)%cycles, int_8)) /= 0)) CYCLE
               image = MATMUL(REAL(dual(:, :, g), dp), work(h)%k)
               delta = MATMUL(REAL(normal, dp), work(i)%k - image)
               IF (MAXVAL(ABS(delta - ANINT(delta))) > tolerance) CYCLE
               status = -3
               IF (MAXVAL(ABS(delta)) > REAL(HUGE(0), dp)) RETURN
               TARGET = NINT(delta, int_8)
               CALL integer_lattice_solve(integer_normal, TARGET, exists, shift, status)
               IF (status /= 0) RETURN
               status = -2
               IF (.NOT. exists) RETURN
               status = -3
               IF (ANY(ABS(shift) > INT(HUGE(0), int_8))) RETURN
               work(i)%image_operation(j) = g
               work(i)%image_shift(:, j) = INT(shift)
               found = .TRUE.
               EXIT
            END DO
            status = -2
            IF (.NOT. found) RETURN
         END DO
      END DO
      CALL MOVE_ALLOC(work, strata)
      status = 0
   END SUBROUTINE reciprocal_strata

! **************************************************************************************************
!> \brief Exact dual of a bounded unimodular three-dimensional rotation.
!> \param rotation Direct-space rotation
!> \param dual Inverse transpose
!> \param status Zero on success; -2 nonunimodular; -3 arithmetic bound exceeded
! **************************************************************************************************
   SUBROUTINE dual_rotation(rotation, dual, status)
      INTEGER, INTENT(IN)                                :: rotation(3, 3)
      INTEGER, INTENT(OUT)                               :: dual(3, 3), status

      INTEGER                                            :: i, i1, i2, j, j1, j2
      INTEGER(int_8)                                     :: cofactor(3, 3), determinant, r(3, 3)

      status = -3
      r = INT(rotation, int_8)
      IF (MAXVAL(ABS(r)) > 100000_int_8) RETURN
      DO j = 1, 3
         j1 = MODULO(j, 3) + 1
         j2 = MODULO(j + 1, 3) + 1
         DO i = 1, 3
            i1 = MODULO(i, 3) + 1
            i2 = MODULO(i + 1, 3) + 1
            cofactor(i, j) = r(i1, j1)*r(i2, j2) - r(i1, j2)*r(i2, j1)
         END DO
      END DO
      determinant = SUM(r(1, :)*cofactor(1, :))
      status = -2
      IF (ABS(determinant) /= 1_int_8) RETURN
      status = -3
      IF (MAXVAL(ABS(cofactor)) > 100000_int_8) RETURN
      dual = INT(cofactor/determinant)
      status = 0
   END SUBROUTINE dual_rotation
END MODULE topology_reciprocal
