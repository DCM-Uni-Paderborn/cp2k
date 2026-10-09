!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
PROGRAM topology_wyckoff_unittest
   USE ieee_arithmetic,                 ONLY: ieee_quiet_nan,&
                                              ieee_value
   USE kinds,                           ONLY: dp
   USE topology_wyckoff,                ONLY: wyckoff_families,&
                                              wyckoff_family_type,&
                                              wyckoff_match_affine

   IMPLICIT NONE
   CALL analytic_tests()
   CALL affine_match_tests()
   PRINT *, 'Wyckoff-family tests passed.'
CONTAINS

! **************************************************************************************************
!> \brief Whole-plane matching, nongeneric origins, duplicate catalogues and invalid numerical data.
! **************************************************************************************************
   SUBROUTINE affine_match_tests()
      INTEGER                                            :: dirs(3, 3, 3), i, ident(3, 3), mult(3), &
                                                            r(3, 3, 2), status, variant
      INTEGER, ALLOCATABLE                               :: mapping(:), operations(:)
      LOGICAL                                            :: expected
      REAL(dp)                                           :: cell(3, 3), positions(3, 3), tau(3, 2)
      TYPE(wyckoff_family_type), ALLOCATABLE             :: families(:)

      ident = 0
      DO i = 1, 3
         ident(i, i) = 1
      END DO
      cell = REAL(ident, dp)
      r(:, :, 1) = ident
      r(:, :, 2) = ident
      r(3, 3, 2) = -1
      tau = 0
      CALL wyckoff_families(cell, r, tau, [.FALSE., .FALSE.], families, 1.e-9_dp, status)
      IF (status /= 0) ERROR STOP 'Affine reference group failed'
      DO variant = 1, 7
         positions = 0
         positions(3, 1) = 0.5_dp
         dirs = 0
         dirs(:, :, 2) = ident
         dirs(:, :2, 1) = ident(:, :2)
         dirs(:, :2, 3) = ident(:, :2)
         mult = [1, 2, 1]
         SELECT CASE (variant)
         CASE (2)
            ! The point is shared, but the incorrect tangent must be rejected.
            dirs(:, 2, 3) = ident(:, 3)
         CASE (3)
            positions(:, 1) = positions(:, 3)
         CASE (4)
            positions(1, 2) = ieee_value(0.0_dp, ieee_quiet_nan)
         CASE (5)
            positions(1, 2) = HUGE(0.0_dp)
         CASE (6)
            mult(3) = 2
         CASE (7)
            dirs(:, :2, 1) = 2*dirs(:, :2, 1)
         END SELECT
         expected = variant == 1 .OR. variant == 7
         CALL wyckoff_match_affine(families, r, tau, positions, dirs, mult, &
                                   mapping, operations, 1.e-9_dp, status)
         IF ((status == 0) .NEQV. expected) ERROR STOP 'Affine family match status wrong'
         IF (ALLOCATED(mapping) .NEQV. expected) ERROR STOP 'Partial affine match escaped'
         IF (ALLOCATED(operations) .NEQV. expected) ERROR STOP 'Partial operation match escaped'
         IF (.NOT. expected) CYCLE
         IF (mapping(2) /= 1) ERROR STOP 'Nongeneric affine origin was used as a point sample'
         IF (mapping(1) == mapping(3)) ERROR STOP 'Distinct mirror planes identified'
      END DO
      PRINT *, 'Complete affine family matching and rejection tests passed.'
   END SUBROUTINE affine_match_tests

! **************************************************************************************************
!> \brief Analytic strata, origin/cell conventions, grey and magnetic groups, rejected malformed data.
! **************************************************************************************************
   SUBROUTINE analytic_tests()
      INTEGER :: changed(3, 3, 16), counts(27), found, g, hist(0:3), i, ident(3, 3), info, &
         inverse(3, 3), j, n, r(3, 3, 16), shear(3, 3), variant
      LOGICAL                                            :: anti(16)
      REAL(dp)                                           :: cell(3, 3), delta(3), image(3), &
                                                            lattice(3, 3), moved(3, 16), &
                                                            origin(3), tau(3, 16)
      TYPE(wyckoff_family_type), ALLOCATABLE             :: families(:), reference(:)

      ident = 0
      DO i = 1, 3
         ident(i, i) = 1
      END DO
      cell = REAL(ident, dp)
      r(:, :, 1) = ident
      tau = 0.0_dp
      anti = .FALSE.
      CALL wyckoff_families(cell, r(:, :, :1), tau(:, :1), anti(:1), families, 1.e-9_dp, info)
      IF (info /= 0) ERROR STOP 'Trivial Wyckoff construction failed'
      IF (SIZE(families) /= 1) ERROR STOP 'P1 should have only its general position'
      IF (families(1)%dimension /= 3 .OR. families(1)%multiplicity /= 1) ERROR STOP 'P1 wrong'
      IF (.NOT. families(1)%maximal_site_symmetry) ERROR STOP 'P1 site group must be maximal'
      r(:, :, 2) = -ident
      CALL wyckoff_families(cell, r(:, :, :2), tau(:, :2), anti(:2), families, 1.e-9_dp, info)
      IF (info /= 0) ERROR STOP 'Inversion Wyckoff construction failed'
      IF (SIZE(families) /= 9) ERROR STOP 'Eight distinct inversion centers were not preserved'
      IF (COUNT(families%dimension == 0) /= 8) ERROR STOP 'Inversion centers lost'
      IF (COUNT(families%maximal_site_symmetry) /= 8) ERROR STOP 'Inversion maximal site groups wrong'
      DO j = 2, SIZE(families)
         delta = 2*families(j)%position
         IF (MAXVAL(ABS(delta - ANINT(delta))) > 1.e-9_dp) ERROR STOP 'Invalid inversion center'
      END DO
      r(:, :, 2) = ident
      r(3, 3, 2) = -1
      CALL wyckoff_families(cell, r(:, :, :2), tau(:, :2), anti(:2), families, 1.e-9_dp, info)
      IF (info /= 0) ERROR STOP 'Mirror Wyckoff construction failed'
      IF (SIZE(families) /= 3 .OR. COUNT(families%dimension == 2) /= 2) ERROR STOP 'Two mirror planes lost'
      r(:, :, 2) = -ident
      r(3, 3, 2) = 1
      tau(3, 2) = 0.5_dp
      CALL wyckoff_families(cell, r(:, :, :2), tau(:, :2), anti(:2), families, 1.e-9_dp, info)
      IF (info /= 0) ERROR STOP 'Screw Wyckoff construction failed'
      IF (SIZE(families) /= 1 .OR. families(1)%multiplicity /= 2) ERROR STOP 'Screw has false fixed points'

      tau = 0.0_dp
      DO g = 1, 8
         r(:, :, g) = 0
         DO i = 1, 3
            r(i, i, g) = MERGE(-1, 1, BTEST(g - 1, i - 1))
         END DO
      END DO
      CALL wyckoff_families(cell, r(:, :, :8), tau(:, :8), anti(:8), reference, 1.e-9_dp, info)
      IF (info /= 0) ERROR STOP 'Pmmm reference failed'
      DO i = 0, 3
         hist(i) = COUNT(reference%dimension == i)
      END DO
      IF (ANY(hist /= [8, 12, 6, 1])) ERROR STOP 'Pmmm intersections incomplete'
      IF (COUNT(reference%maximal_site_symmetry) /= 8) ERROR STOP 'Pmmm maximal site groups wrong'
      DO j = 1, SIZE(reference)
         SELECT CASE (reference(j)%dimension)
         CASE (0)
            n = 0
         CASE (1)
            n = 2
         CASE (2)
            n = 8
         CASE (3)
            n = 26
         END SELECT
         IF (SIZE(reference(j)%specializations) /= n) ERROR STOP 'Pmmm specialization graph incomplete'
      END DO
      shear = RESHAPE([1, 0, 0, 1, 1, 0, 0, 1, 1], [3, 3])
      inverse = RESHAPE([1, 0, 0, -1, 1, 0, 1, -1, 1], [3, 3])
      origin = [0.123_dp, -0.217_dp, 0.319_dp]
      DO variant = 1, 4
         moved = 0.0_dp
         anti = .FALSE.
         lattice = cell
         n = 8
         DO g = 1, 8
            changed(:, :, g) = r(:, :, g)
            SELECT CASE (variant)
            CASE (1)
               moved(:, g) = MATMUL(REAL(r(:, :, g) - ident, dp), origin)
            CASE (2)
               changed(:, :, g) = MATMUL(MATMUL(inverse, r(:, :, g)), shear)
               lattice = MATMUL(cell, REAL(shear, dp))
            CASE (3)
               IF (g /= 1) moved(:, g) = REAL([g, 1 - g, 2*g], dp)
            CASE (4)
               changed(:, :, g + 8) = r(:, :, g)
               anti(g + 8) = .TRUE.
               n = 16
            END SELECT
         END DO
         CALL wyckoff_families(lattice, changed(:, :, :n), moved(:, :n), anti(:n), families, 1.e-9_dp, info)
         IF (info /= 0) ERROR STOP 'Transformed Wyckoff enumeration failed'
         IF (SIZE(families) /= SIZE(reference)) ERROR STOP 'Wyckoff count depends on conventions'
         counts = 0
         DO j = 1, SIZE(families)
            image = families(j)%position
            IF (variant == 1) image = image + origin
            IF (variant == 2) image = MATMUL(REAL(shear, dp), image)
            found = 0
            DO i = 1, SIZE(reference)
               IF (reference(i)%dimension /= families(j)%dimension) CYCLE
               IF (reference(i)%multiplicity /= families(j)%multiplicity) CYCLE
               IF (reference(i)%maximal_site_symmetry .NEQV. families(j)%maximal_site_symmetry) CYCLE
               IF (SIZE(reference(i)%specializations) /= SIZE(families(j)%specializations)) CYCLE
               delta = MATMUL(REAL(reference(i)%normal, dp), image) - reference(i)%offset
               IF (MAXVAL(ABS(delta - ANINT(delta))) > 1.e-8_dp) CYCLE
               found = found + 1
               counts(i) = counts(i) + 1
            END DO
            IF (found /= 1) ERROR STOP 'Transformed family is not unique'
         END DO
         IF (ANY(counts /= 1)) ERROR STOP 'Transformed family coverage is not bijective'
      END DO
      ! A half translation combined with time reversal has no additional spatial stabilizer.
      changed(:, :, 1) = ident
      changed(:, :, 2) = ident
      moved(:, 1) = 0.0_dp
      moved(:, 2) = 0.5_dp
      anti(:2) = [.FALSE., .TRUE.]
      CALL wyckoff_families(cell, changed(:, :, :2), moved(:, :2), anti(:2), families, 1.e-9_dp, info)
      IF (info /= 0) ERROR STOP 'Magnetic half-translation construction failed'
      IF (SIZE(families) /= 1 .OR. families(1)%multiplicity /= 2) ERROR STOP 'Magnetic fixed set wrong'
      moved(:, 2) = 0.3_dp
      CALL wyckoff_families(cell, changed(:, :, :2), moved(:, :2), anti(:2), families, 1.e-9_dp, info)
      IF (info == 0 .OR. ALLOCATED(families)) ERROR STOP 'Nonclosed group accepted'
      PRINT *, 'Wyckoff analytic, coordinate and magnetic tests passed.'
   END SUBROUTINE analytic_tests

END PROGRAM topology_wyckoff_unittest
