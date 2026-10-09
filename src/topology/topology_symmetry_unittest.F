!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
PROGRAM topology_symmetry_unittest
   USE ieee_arithmetic,                 ONLY: ieee_quiet_nan,&
                                              ieee_value
   USE iso_fortran_env,                 ONLY: output_unit
   USE kinds,                           ONLY: dp
   USE kpoint_mo_symmetry_methods,      ONLY: kpoint_cell_gauge_action
   USE kpoint_types,                    ONLY: kpoint_create,&
                                              kpoint_release,&
                                              kpoint_sym_create,&
                                              kpoint_type
   USE mathlib,                         ONLY: inv_3x3
   USE topology_compatibility,          ONLY: little_group_compatibility,&
                                              little_group_data_type,&
                                              little_group_edge_type
   USE topology_corepresentations,      ONLY: corepresentation_multiplicities,&
                                              corepresentation_table_type,&
                                              projective_corepresentations
   USE topology_little_group,           ONLY: little_group_factor_system,&
                                              projective_characters,&
                                              validate_factor_system
   USE topology_symmetry,               ONLY: character_multiplicities,&
                                              check_sewing_product,&
                                              inversion_representation,&
                                              spin_rotation,&
                                              subspace_sewing,&
                                              transform_sewing_link,&
                                              transform_spinor_frame
   USE topology_tqc,                    ONLY: inversion_ebr_signature,&
                                              inversion_indicators

   IMPLICIT NONE
   COMPLEX(KIND=dp) :: table(6, 3), group_chars(6), metric(2, 2), coeff(4, 4), phase(2), tmp(4)
   REAL(KIND=dp) :: error, energy_error, energies(4)
   INTEGER :: multiplicities(3), status, counts(2), odd(8), ebr(2, 8), a, k, sign, strong, weak(3), z4, io_unit
   LOGICAL :: signed_atomic, nonnegative_atomic

   io_unit = output_unit
   CALL sewing_tests()
   CALL spinor_tests()
   CALL cell_gauge_tests()
   CALL little_group_tests()
   CALL projective_character_tests()
   CALL compatibility_tests()
   CALL corepresentation_tests()
   ! Complete C3v/S3 character table per element, including a two-dimensional irrep.
   table(:, 1) = CMPLX([1, 1, 1, 1, 1, 1], 0, dp)
   table(:, 2) = CMPLX([1, 1, 1, -1, -1, -1], 0, dp)
   table(:, 3) = CMPLX([2, -1, -1, 0, 0, 0], 0, dp)
   group_chars = table(:, 1) + 2*table(:, 3)
   CALL character_multiplicities(group_chars, table, multiplicities, 1.e-10_dp, status)
   IF (status /= 0 .OR. ANY(multiplicities /= [1, 0, 2])) ERROR STOP 'C3v decomposition failed'
   group_chars(2) = group_chars(2) + 0.1_dp
   CALL character_multiplicities(group_chars, table, multiplicities, 1.e-10_dp, status)
   IF (status == 0) ERROR STOP 'Invalid character accepted'
   metric = CMPLX(0, 0, dp)
   metric(1, 1) = 2.0_dp
   metric(2, 2) = 0.5_dp
   coeff = CMPLX(0, 0, dp)
   coeff(1, 1) = SQRT(0.5_dp)
   coeff(2, 2) = SQRT(2.0_dp)
   coeff(3, 3) = SQRT(0.5_dp)
   coeff(4, 4) = SQRT(2.0_dp)
   phase = CMPLX([1, -1], 0, dp)
   energies = 0.0_dp
   CALL inversion_representation(metric, coeff, [1, 2], phase, energies, .TRUE., &
                                 1.e-9_dp, 1.e-9_dp, counts, error, energy_error, status)
   IF (status /= 0 .OR. ANY(counts /= [2, 2])) ERROR STOP 'Metric/spinor parity failed'
   tmp = coeff(:, 1)
   coeff(:, 1) = (tmp + CMPLX(0, 1, dp)*coeff(:, 2))/SQRT(2.0_dp)
   coeff(:, 2) = (CMPLX(0, 1, dp)*tmp + coeff(:, 2))/SQRT(2.0_dp)
   CALL inversion_representation(metric, coeff, [1, 2], phase, energies, .TRUE., &
                                 1.e-9_dp, 1.e-9_dp, counts, error, energy_error, status)
   IF (status /= 0 .OR. ANY(counts /= [2, 2])) ERROR STOP 'Degenerate-gauge parity failed'
   CALL inversion_representation(metric, coeff(:, :1), [1, 2], phase, energies(:1), .FALSE., &
                                 1.e-9_dp, 1.e-9_dp, counts, error, energy_error, status)
   IF (status == 0) ERROR STOP 'Nonclosed parity subspace accepted'
   CALL inversion_representation(metric, coeff(:, :2), [1, 2], phase, energies(:2), .TRUE., &
                                 1.e-9_dp, 1.e-9_dp, counts, error, energy_error, status)
   IF (status == 0) ERROR STOP 'Broken physical time reversal accepted'
   energies(2) = 1.0_dp
   CALL inversion_representation(metric, coeff, [1, 2], phase, energies, .TRUE., &
                                 1.e-9_dp, 1.e-9_dp, counts, error, energy_error, status)
   IF (status == 0) ERROR STOP 'Noncommuting eigenspace accepted'
   ! Scalar two-site inversion with a genuinely off-diagonal AO metric.
   metric = CMPLX(0.2_dp, 0, dp)
   metric(1, 1) = 1.0_dp
   metric(2, 2) = 1.0_dp
   coeff(:2, 1) = CMPLX([1.0_dp, 1.0_dp]/SQRT(2.4_dp), 0, dp)
   coeff(:2, 2) = CMPLX([1.0_dp, -1.0_dp]/SQRT(1.6_dp), 0, dp)
   phase = CMPLX(1, 0, dp)
   CALL inversion_representation(metric, coeff(:2, :2), [2, 1], phase, energies(:2), .FALSE., &
                                 1.e-9_dp, 1.e-9_dp, counts, error, energy_error, status)
   IF (status /= 0 .OR. ANY(counts /= [1, 1])) ERROR STOP 'Scalar two-site inversion failed'
   CALL inversion_representation(metric, coeff(:2, :1), [2, 1], -phase, energies(:1), .FALSE., &
                                 1.e-9_dp, 1.e-9_dp, counts, error, energy_error, status)
   IF (status /= 0 .OR. ANY(counts /= [0, 1])) ERROR STOP 'Zone-boundary inversion phase failed'
   metric(1, 2) = CMPLX(0.2_dp, 0.01_dp, dp)
   CALL inversion_representation(metric, coeff(:2, :2), [2, 1], phase, energies(:2), .FALSE., &
                                 1.e-9_dp, 1.e-9_dp, counts, error, energy_error, status)
   IF (status == 0) ERROR STOP 'Invalid metric accepted'
   metric(1, 2) = CMPLX(ieee_value(0.0_dp, ieee_quiet_nan), 0, dp)
   CALL inversion_representation(metric, coeff(:2, :2), [2, 1], phase, energies(:2), .FALSE., &
                                 1.e-9_dp, 1.e-9_dp, counts, error, energy_error, status)
   IF (status == 0) ERROR STOP 'Nonfinite metric accepted'
   DO a = 0, 7
      DO sign = -1, 1, 2
         DO k = 0, 7
            odd(k + 1) = (1 - sign*(-1)**POPCNT(IAND(a, k)))/2
         END DO
         CALL inversion_indicators(odd, 1, 3, strong, weak, z4, status)
         IF (status /= 0 .OR. strong /= 0 .OR. ANY(weak /= 0) .OR. z4 /= 0) ERROR STOP 'Atomic indicator failed'
         CALL inversion_ebr_signature(odd, 1, signed_atomic, nonnegative_atomic, ebr, status)
         IF (status /= 0 .OR. .NOT. nonnegative_atomic .OR. .NOT. signed_atomic) ERROR STOP 'Atomic EBR failed'
         IF (SUM(ebr) /= 1) ERROR STOP 'EBR rank failed'
      END DO
   END DO
   odd = 0
   odd(1) = 1
   CALL inversion_indicators(odd, 1, 3, strong, weak, z4, status)
   IF (status /= 0 .OR. strong /= 1 .OR. z4 /= 1) ERROR STOP 'Strong index failed'
   CALL inversion_ebr_signature(odd, 1, signed_atomic, nonnegative_atomic, ebr, status)
   IF (status /= 0 .OR. signed_atomic .OR. nonnegative_atomic) ERROR STOP 'Stable obstruction missed'
   odd(1) = 2
   CALL inversion_indicators(odd, 2, 3, strong, weak, z4, status)
   IF (status /= 0 .OR. strong /= 0 .OR. z4 /= 2) ERROR STOP 'Inversion Z4=2 missed'
   odd(1) = 4
   CALL inversion_ebr_signature(odd, 4, signed_atomic, nonnegative_atomic, ebr, status)
   IF (status /= 0 .OR. .NOT. signed_atomic .OR. nonnegative_atomic) ERROR STOP 'Signed-only EBR missed'
   odd = -1
   CALL inversion_ebr_signature(odd, 4, signed_atomic, nonnegative_atomic, ebr, status)
   IF (status == 0) ERROR STOP 'Invalid multiplicity accepted'
   CALL exhaustive_signatures(2)
   CALL exhaustive_signatures(3)
   WRITE (io_unit, *) 'Native characters, parity, inversion indicators and EBR tests passed.'
CONTAINS

! **************************************************************************************************
!> \brief Wigner types, magnetic translations, quaternion doublets and twisted-gauge invariance.
! **************************************************************************************************
   SUBROUTINE corepresentation_tests()
      COMPLEX(dp)                                        :: phase, phases(6)
      COMPLEX(dp), ALLOCATABLE                           :: characters(:, :), factors(:, :), &
                                                            rephased(:, :)
      INTEGER                                            :: boundary, g, h, i, info, inverse(6), nu, &
                                                            original(3, 3), perm(6), q, &
                                                            rotations(3, 3, 8), spin
      INTEGER, ALLOCATABLE                               :: counts(:), mapping(:), mult(:), &
                                                            products(:, :), reordered(:, :)
      LOGICAL                                            :: anti(8)
      REAL(dp)                                           :: cell(3, 3), error, k(3), tau(3, 8)
      TYPE(corepresentation_table_type)                  :: table

      cell = 0.0_dp
      rotations = 0
      tau = 0.0_dp
      anti = .FALSE.
      DO i = 1, 3
         cell(i, i) = 1.0_dp
         rotations(i, i, :) = 1
      END DO
      anti(2) = .TRUE.
      tau(3, 2) = 0.5_dp
      DO spin = 0, 1
         DO boundary = 0, 1
            k = [0.0_dp, 0.0_dp, 0.5_dp*boundary]
            CALL little_group_factor_system(cell, rotations(:, :, :2), tau(:, :2), anti(:2), k, spin == 1, &
                                            mapping, products, factors, 1.e-10_dp, info)
            IF (info /= 0) ERROR STOP 'Magnetic half-translation group failed'
            CALL projective_characters(products(:1, :1), factors(:1, :1), characters, 1.e-9_dp, error, info)
            IF (info /= 0) ERROR STOP 'Magnetic-translation unitary character failed'
            CALL projective_corepresentations(products, factors, anti(:2), characters, table, 1.e-9_dp, error, info)
            IF (info /= 0) ERROR STOP 'Magnetic-translation corep failed'
            q = (-1)**(spin + boundary)
            IF (ANY(table%indicator /= q)) ERROR STOP 'Spin/Bloch Wigner sign lost'
            IF (SIZE(table%indicator) /= 1 .OR. ANY(table%dimension /= (3 - q)/2)) ERROR STOP 'Corep doubling wrong'
            ALLOCATE (counts(1))
            CALL corepresentation_multiplicities(table, [(3 - q)/2], counts, info)
            IF (info /= 0 .OR. ANY(counts /= 1)) ERROR STOP 'Magnetic-translation counts failed'
            IF (q == -1) THEN
               CALL corepresentation_multiplicities(table, [1], counts, info)
               IF (info == 0) ERROR STOP 'Odd type-b multiplicity accepted'
            END IF
            DEALLOCATE (counts)
         END DO
      END DO
      ! The quaternion spin irrep already contains the Kramers doublet: no additional doubling.
      tau = 0.0_dp
      k = 0.0_dp
      DO g = 2, 4
         DO i = 1, 3
            rotations(i, i, g) = -1
         END DO
         rotations(g - 1, g - 1, g) = 1
      END DO
      rotations(:, :, 5:8) = rotations(:, :, :4)
      anti(:4) = .FALSE.
      anti(5:8) = .TRUE.
      CALL little_group_factor_system(cell, rotations, tau, anti, k, .TRUE., mapping, products, factors, 1.e-10_dp, info)
      IF (info /= 0) ERROR STOP 'Spin D2 grey group failed'
      CALL projective_characters(products(:4, :4), factors(:4, :4), characters, 1.e-9_dp, error, info)
      IF (info /= 0) ERROR STOP 'Quaternion unitary irrep failed'
      CALL projective_corepresentations(products, factors, anti, characters, table, 1.e-9_dp, error, info)
      IF (info /= 0) ERROR STOP 'Quaternion corep failed'
      IF (SIZE(table%indicator) /= 1 .OR. ANY(table%indicator /= 1)) ERROR STOP 'Quaternion Wigner type wrong'
      IF (ANY(table%dimension /= 2) .OR. ANY(table%restriction /= 1)) ERROR STOP 'Kramers doublet doubled twice'

      ! C3 plus spinless time reversal pairs its two complex one-dimensional irreps.
      cell(1, 2) = -0.5_dp
      cell(2, 2) = SQRT(3.0_dp)/2
      rotations(:, :, 2) = RESHAPE([0, 1, 0, -1, -1, 0, 0, 0, 1], [3, 3])
      rotations(:, :, 3) = MATMUL(rotations(:, :, 2), rotations(:, :, 2))
      rotations(:, :, 4:6) = rotations(:, :, :3)
      anti(:3) = .FALSE.
      anti(4:6) = .TRUE.
      CALL little_group_factor_system(cell, rotations(:, :, :6), tau(:, :6), anti(:6), k, .FALSE., &
                                      mapping, products, factors, 1.e-10_dp, info)
      IF (info /= 0) ERROR STOP 'C3 grey group failed'
      nu = 3
      CALL projective_characters(products(:nu, :nu), factors(:nu, :nu), characters, 1.e-9_dp, error, info)
      IF (info /= 0) ERROR STOP 'C3 unitary characters failed'
      CALL projective_corepresentations(products, factors, anti(:6), characters, table, 1.e-9_dp, error, info)
      IF (info /= 0) ERROR STOP 'C3 coreps failed'
      IF (COUNT(table%indicator == 1) /= 1 .OR. COUNT(table%indicator == 0) /= 1) ERROR STOP 'C3 conjugate pairing wrong'
      IF (SUM(table%dimension) /= 3) ERROR STOP 'C3 dimensions wrong'
      ALLOCATE (counts(2), mult(3))
      mult = 1
      CALL corepresentation_multiplicities(table, mult, counts, info)
      IF (info /= 0 .OR. ANY(counts /= 1)) ERROR STOP 'C3 multiplicities failed'
      q = FINDLOC(table%indicator, 0, DIM=1)
      i = FINDLOC(table%restriction(:, q), 1, DIM=1)
      mult(i) = 0
      CALL corepresentation_multiplicities(table, mult, counts, info)
      IF (info == 0) ERROR STOP 'Unpaired complex irrep accepted'
      original = MATMUL(table%restriction, TRANSPOSE(table%restriction))
      DO g = 1, 6
         phases(g) = EXP(CMPLX(0.0_dp, 0.17_dp*(g - 1)**2, dp))
      END DO
      rephased = factors
      DO h = 1, 6
         DO g = 1, 6
            phase = phases(h)
            IF (anti(g)) phase = CONJG(phase)
            rephased(g, h) = phases(g)*phase*CONJG(phases(products(g, h)))*factors(g, h)
         END DO
      END DO
      characters = SPREAD(phases(:3), 2, 3)*characters
      characters = characters(:, [3, 1, 2])
      CALL projective_corepresentations(products, rephased, anti(:6), characters, table, 1.e-9_dp, error, info)
      IF (info /= 0) ERROR STOP 'Twisted corep gauge failed'
      IF (ANY(MATMUL(table%restriction, TRANSPOSE(table%restriction)) /= original([3, 1, 2], [3, 1, 2]))) THEN
         ERROR STOP 'Corep pairing changed under phase/irrep gauge'
      END IF
      perm = [1, 4, 2, 5, 3, 6]
      inverse(perm) = [1, 2, 3, 4, 5, 6]
      reordered = products
      DO h = 1, 6
         DO g = 1, 6
            reordered(g, h) = inverse(products(perm(g), perm(h)))
         END DO
      END DO
      CALL projective_corepresentations(reordered, rephased(perm, perm), anti(perm), characters, table, &
                                        1.e-9_dp, error, info)
      IF (info /= 0) ERROR STOP 'Interleaved operation order failed'
      IF (ANY(MATMUL(table%restriction, TRANSPOSE(table%restriction)) /= original([3, 1, 2], [3, 1, 2]))) THEN
         ERROR STOP 'Corep pairing changed with operation order'
      END IF
      rephased(4, 4) = rephased(4, 4)*CMPLX(0.0_dp, 1.0_dp, dp)
      CALL projective_corepresentations(products, rephased, anti(:6), characters, table, 1.e-9_dp, error, info)
      IF (info == 0 .OR. ALLOCATED(table%restriction)) ERROR STOP 'Invalid antiunitary factor accepted'
      WRITE (io_unit, '(A)') 'Wigner corepresentations, pairing, magnetic translations and gauge tests passed.'
   END SUBROUTINE corepresentation_tests

! **************************************************************************************************
!> \brief Screw subduction, monodromy, representative gauges and non-line endpoint symmetries.
! **************************************************************************************************
   SUBROUTINE compatibility_tests()
      COMPLEX(dp)                                        :: phase, root
      INTEGER                                            :: i, info, rotations(3, 3, 2), spin
      REAL(dp)                                           :: angle, cell(3, 3), error, k(3), tau(3, 2)
      TYPE(little_group_data_type)                       :: left, right
      TYPE(little_group_edge_type)                       :: edge

      cell = 0.0_dp
      rotations = 0
      DO i = 1, 3
         cell(i, i) = 1.0_dp
         rotations(i, i, :) = 1
      END DO
      rotations(1, 1, 2) = -1
      rotations(2, 2, 2) = -1
      tau = 0.0_dp
      tau(3, 2) = 0.5_dp
      DO spin = 0, 1
         k = 0.0_dp
         CALL make_little_point(cell, rotations, tau, k, spin == 1, left)
         k(3) = 0.25_dp
         CALL make_little_point(cell, rotations, tau, k, spin == 1, right)
         root = 1.0_dp
         IF (spin == 1) root = CMPLX(0.0_dp, 1.0_dp, dp)
         angle = ACOS(-1.0_dp)/4
         phase = CMPLX(COS(angle), SIN(angle), dp)
         left%multiplicity(:) = MERGE(1, 0, REAL(left%characters(2, :)/root, dp) > 0.5_dp)
         right%multiplicity(:) = MERGE(1, 0, REAL(phase*right%characters(2, :)/root, dp) > 0.5_dp)
         IF (SUM(left%multiplicity) /= 1 .OR. SUM(right%multiplicity) /= 1) ERROR STOP 'Screw sector selection failed'
         CALL little_group_compatibility(left, right, [0, 0, 0], edge, 1.e-9_dp, error, info)
         IF (info /= 0) ERROR STOP 'Screw compatibility algebra failed'
         IF (SIZE(edge%left) /= 2 .OR. ANY(edge%left_count /= edge%right_count)) ERROR STOP 'Screw phases lost'
         right%multiplicity(:) = 1 - right%multiplicity
         CALL little_group_compatibility(left, right, [0, 0, 0], edge, 1.e-9_dp, error, info)
         IF (info /= 0) ERROR STOP 'Screw mismatch treated as invalid data'
         IF (ALL(edge%left_count == edge%right_count)) ERROR STOP 'Incompatible sectors missed'
         CALL little_group_compatibility(left, left, [0, 0, 1], edge, 1.e-9_dp, error, info)
         IF (info /= 0) ERROR STOP 'Screw monodromy algebra failed'
         IF (ALL(edge%left_count == edge%right_count)) ERROR STOP 'Reciprocal-image monodromy lost'
         left%multiplicity(:) = 1
         CALL little_group_compatibility(left, left, [0, 0, 1], edge, 1.e-9_dp, error, info)
         IF (info /= 0) ERROR STOP 'Complete screw pair rejected'
         IF (ANY(edge%left_count /= edge%right_count)) ERROR STOP 'Complete screw pair is incompatible'
      END DO
      k = 0.0_dp
      CALL make_little_point(cell, rotations, tau, k, .FALSE., left)
      left%multiplicity(:) = MERGE(1, 0, REAL(left%characters(2, :), dp) > 0.5_dp)
      ! Right-hand coset representatives may differ by an integer translation.
      tau(3, 2) = 1.5_dp
      k(3) = 0.25_dp
      CALL make_little_point(cell, rotations, tau, k, .FALSE., right)
      phase = EXP(CMPLX(0.0_dp, 3*ACOS(-1.0_dp)/4, dp))
      right%multiplicity(:) = MERGE(1, 0, REAL(phase*right%characters(2, :), dp) > 0.5_dp)
      right%characters(:, :) = right%characters(:, [2, 1])
      right%multiplicity(:) = right%multiplicity([2, 1])
      CALL little_group_compatibility(left, right, [0, 0, 0], edge, 1.e-9_dp, error, info)
      IF (info /= 0) ERROR STOP 'Representative/irrep-order gauge failed'
      IF (ANY(edge%left_count /= edge%right_count)) ERROR STOP 'Representative/irrep-order dependence'
      right%factors(2, 2) = -right%factors(2, 2)
      CALL little_group_compatibility(left, right, [0, 0, 0], edge, 1.e-9_dp, error, info)
      IF (info == 0) ERROR STOP 'Different line factor gauges accepted'
      tau(3, 2) = 0.5_dp
      k = [0.5_dp, 0.0_dp, 0.0_dp]
      CALL make_little_point(cell, rotations, tau, k, .FALSE., right)
      right%multiplicity(:) = MERGE(1, 0, REAL(right%characters(2, :), dp) < -0.5_dp)
      CALL little_group_compatibility(left, right, [0, 0, 0], edge, 1.e-9_dp, error, info)
      IF (info /= 0) ERROR STOP 'Endpoint-only symmetry test failed'
      IF (SIZE(edge%left) /= 1 .OR. ANY(edge%left_count /= edge%right_count)) THEN
         ERROR STOP 'Endpoint intersection substituted for segment stabilizer'
      END IF
      right%k(1) = ieee_value(0.0_dp, ieee_quiet_nan)
      CALL little_group_compatibility(left, right, [0, 0, 0], edge, 1.e-9_dp, error, info)
      IF (info == 0) ERROR STOP 'Nonfinite segment accepted'
      CALL degenerate_compatibility_test()
      CALL antiunitary_compatibility_tests()
      WRITE (io_unit, '(A)') 'Unitary band compatibility and nonsymmorphic monodromy: passed.'
   END SUBROUTINE compatibility_tests

! **************************************************************************************************
!> \brief A C4v doublet splits into one even and one odd branch on a mirror-invariant line.
! **************************************************************************************************
   SUBROUTINE degenerate_compatibility_test()
      INTEGER                                            :: g, i, info, mirror(3, 3), quarter(3, 3), &
                                                            rotations(3, 3, 8)
      REAL(dp)                                           :: cell(3, 3), error, tau(3, 8)
      TYPE(little_group_data_type)                       :: left, right
      TYPE(little_group_edge_type)                       :: edge

      cell = 0.0_dp
      rotations = 0
      tau = 0.0_dp
      DO i = 1, 3
         cell(i, i) = 1.0_dp
         rotations(i, i, 1) = 1
      END DO
      quarter = RESHAPE([0, 1, 0, -1, 0, 0, 0, 0, 1], [3, 3])
      mirror = rotations(:, :, 1)
      mirror(2, 2) = -1
      DO g = 2, 4
         rotations(:, :, g) = MATMUL(quarter, rotations(:, :, g - 1))
      END DO
      DO g = 1, 4
         rotations(:, :, g + 4) = MATMUL(rotations(:, :, g), mirror)
      END DO
      CALL make_little_point(cell, rotations, tau, [0.0_dp, 0.0_dp, 0.0_dp], .FALSE., left)
      CALL make_little_point(cell, rotations, tau, [0.25_dp, 0.0_dp, 0.0_dp], .FALSE., right)
      left%multiplicity = MERGE(1, 0, NINT(REAL(left%characters(1, :), dp)) == 2)
      IF (SUM(left%multiplicity) /= 1) ERROR STOP 'C4v doublet missing'
      right%multiplicity = 1
      CALL little_group_compatibility(left, right, [0, 0, 0], edge, 1.e-9_dp, error, info)
      IF (info /= 0) ERROR STOP 'C4v to mirror subduction failed'
      IF (SIZE(edge%left_count) /= 2) ERROR STOP 'Mirror irreps missing'
      IF (ANY(edge%left_count /= 1) .OR. ANY(edge%right_count /= 1)) ERROR STOP 'C4v doublet branching wrong'
      ! Equal rank alone is insufficient: two copies of one mirror branch are incompatible.
      right%multiplicity = [2, 0]
      CALL little_group_compatibility(left, right, [0, 0, 0], edge, 1.e-9_dp, error, info)
      IF (info /= 0) ERROR STOP 'Valid mirror mismatch rejected'
      IF (ALL(edge%left_count == edge%right_count)) ERROR STOP 'Mirror-sector obstruction missed'
   END SUBROUTINE degenerate_compatibility_test

! **************************************************************************************************
!> \brief PT line pairs, endpoint-only time reversal, glide Kramers planes and hourglass switching.
! **************************************************************************************************
   SUBROUTINE antiunitary_compatibility_tests()
      INTEGER                                            :: boundary, i, info, rotations(3, 3, 4), &
                                                            spin
      REAL(dp)                                           :: cell(3, 3), error, k(3), tau(3, 4)
      TYPE(little_group_data_type)                       :: left, right
      TYPE(little_group_edge_type)                       :: edge

      cell = 0.0_dp
      rotations = 0
      tau = 0.0_dp
      DO i = 1, 3
         cell(i, i) = 1.0_dp
         rotations(i, i, :) = 1
      END DO
      rotations(:, :, 2) = -rotations(:, :, 1)
      tau(:, 2) = [0.21_dp, 0.37_dp, 0.11_dp]
      DO spin = 0, 1
         CALL make_little_point(cell, rotations(:, :, :2), tau(:, :2), [0.173_dp, 0.217_dp, 0.319_dp], &
                                spin == 1, left, [.FALSE., .TRUE.])
         CALL make_little_point(cell, rotations(:, :, :2), tau(:, :2), [0.273_dp, 0.117_dp, 0.219_dp], &
                                spin == 1, right, [.FALSE., .TRUE.])
         left%multiplicity(:) = spin + 1
         right%multiplicity(:) = spin + 1
         CALL little_group_compatibility(left, right, [0, 0, 0], edge, 1.e-9_dp, error, info)
         IF (info /= 0) ERROR STOP 'PT line twisted rephasing failed'
         IF (ANY(edge%coreps%indicator /= (-1)**spin)) ERROR STOP 'PT line Wigner type wrong'
         IF (ANY(edge%left_corep_count /= 1) .OR. ANY(edge%right_corep_count /= 1)) ERROR STOP 'PT line pair lost'
         right%multiplicity(:) = 2*(spin + 1)
         CALL little_group_compatibility(left, right, [0, 0, 0], edge, 1.e-9_dp, error, info)
         IF (info /= 0) ERROR STOP 'Valid antiunitary mismatch rejected'
         IF (ALL(edge%left_corep_count == edge%right_corep_count)) ERROR STOP 'PT line mismatch lost'
         IF (spin == 1) THEN
            right%multiplicity(:) = 1
            CALL little_group_compatibility(left, right, [0, 0, 0], edge, 1.e-9_dp, error, info)
            IF (info == 0) ERROR STOP 'Unpaired PT endpoint accepted'
         END IF
      END DO
      ! Time reversal fixes both TRIM endpoints, but no interior point of the segment.
      rotations(:, :, 2) = rotations(:, :, 1)
      tau = 0.0_dp
      CALL make_little_point(cell, rotations(:, :, :2), tau(:, :2), [0.0_dp, 0.0_dp, 0.0_dp], &
                             .TRUE., left, [.FALSE., .TRUE.])
      CALL make_little_point(cell, rotations(:, :, :2), tau(:, :2), [0.5_dp, 0.0_dp, 0.0_dp], &
                             .TRUE., right, [.FALSE., .TRUE.])
      left%multiplicity(:) = 2
      right%multiplicity(:) = 2
      CALL little_group_compatibility(left, right, [0, 0, 0], edge, 1.e-9_dp, error, info)
      IF (info /= 0) ERROR STOP 'Endpoint time reversal subduction failed'
      IF (ALLOCATED(edge%coreps%indicator)) ERROR STOP 'Endpoint time reversal leaked onto line'
      IF (ANY(edge%left_corep_branch /= 2)) ERROR STOP 'Kramers pair did not split on unitary line'
      ! Antiunitary {m_z|x/2}T preserves k_z lines; its square is -1 on k_x=1/2.
      rotations(3, 3, 2) = -1
      tau(1, 2) = 0.5_dp
      DO spin = 0, 1
         DO boundary = 0, 1
            k = [0.5_dp*boundary, 0.0_dp, 0.0_dp]
            CALL make_little_point(cell, rotations(:, :, :2), tau(:, :2), k, spin == 1, left, [.FALSE., .TRUE.])
            k(3) = 0.25_dp
            CALL make_little_point(cell, rotations(:, :, :2), tau(:, :2), k, spin == 1, right, [.FALSE., .TRUE.])
            left%multiplicity(:) = 1 + boundary
            right%multiplicity(:) = 1 + boundary
            CALL little_group_compatibility(left, right, [0, 0, 0], edge, 1.e-9_dp, error, info)
            IF (info /= 0) ERROR STOP 'Antiunitary glide line failed'
            IF (ANY(edge%coreps%indicator /= (-1)**boundary)) ERROR STOP 'Glide Kramers plane lost'
         END DO
      END DO
      ! Spinful {m_y|z/2} plus T: opposite glide partners at Gamma, equal partners at Z.
      rotations(:, :, 2) = rotations(:, :, 1)
      rotations(2, 2, 2) = -1
      rotations(:, :, 4) = rotations(:, :, 2)
      tau = 0.0_dp
      tau(3, 2) = 0.5_dp
      tau(3, 4) = 0.5_dp
      CALL make_little_point(cell, rotations, tau, [0.0_dp, 0.0_dp, 0.0_dp], .TRUE., left, &
                             [.FALSE., .FALSE., .TRUE., .TRUE.])
      CALL make_little_point(cell, rotations, tau, [0.0_dp, 0.0_dp, 0.5_dp], .TRUE., right, &
                             [.FALSE., .FALSE., .TRUE., .TRUE.])
      left%multiplicity(:) = 1
      right%multiplicity(:) = [2, 0]
      CALL little_group_compatibility(left, right, [0, 0, 0], edge, 1.e-9_dp, error, info)
      IF (info /= 0) ERROR STOP 'Hourglass endpoint subduction failed'
      IF (ANY(edge%left_corep_branch /= 1)) ERROR STOP 'Gamma glide partner pair lost'
      IF (ANY(SUM(edge%right_corep_branch, DIM=1) /= 2)) ERROR STOP 'Z glide doublets lost'
      IF (ALL(edge%left_count == edge%right_count)) ERROR STOP 'Two-band hourglass obstruction missed'
      left%multiplicity(:) = 2
      right%multiplicity(:) = 2
      CALL little_group_compatibility(left, right, [0, 0, 0], edge, 1.e-9_dp, error, info)
      IF (info /= 0 .OR. ANY(edge%left_count /= edge%right_count)) ERROR STOP 'Four-band hourglass rejected'
      WRITE (io_unit, '(A)') 'Antiunitary segment subduction and hourglass connectivity: passed.'
   END SUBROUTINE antiunitary_compatibility_tests

! **************************************************************************************************
!> \brief Generate a complete unitary test point with one regular representation.
!> \param cell Lattice vectors as columns
!> \param rotations Full integer direct rotations
!> \param tau Full fractional translations
!> \param k Fractional reciprocal point
!> \param spinful Include the spin factor
!> \param point Validated point data
!> \param antiunitary Optional full-group flags, with unitary operations first
! **************************************************************************************************
   SUBROUTINE make_little_point(cell, rotations, tau, k, spinful, point, antiunitary)
      REAL(dp), INTENT(IN)                               :: cell(3, 3)
      INTEGER, INTENT(IN)                                :: rotations(:, :, :)
      REAL(dp), INTENT(IN)                               :: tau(:, :), k(3)
      LOGICAL, INTENT(IN)                                :: spinful
      TYPE(little_group_data_type), INTENT(OUT)          :: point
      LOGICAL, INTENT(IN), OPTIONAL                      :: antiunitary(:)

      INTEGER                                            :: info, n, nu
      INTEGER, ALLOCATABLE                               :: members(:)
      LOGICAL, ALLOCATABLE                               :: anti(:)
      REAL(dp)                                           :: residual

      ALLOCATE (anti(SIZE(tau, 2)))
      anti = .FALSE.
      IF (PRESENT(antiunitary)) anti(:) = antiunitary
      CALL little_group_factor_system(cell, rotations, tau, anti, k, spinful, &
                                      members, point%product, point%factors, 1.e-9_dp, info)
      IF (info /= 0) ERROR STOP 'Compatibility test point factors failed'
      nu = COUNT(.NOT. anti(members))
      CALL projective_characters(point%product(:nu, :nu), point%factors(:nu, :nu), point%characters, &
                                 1.e-9_dp, residual, info)
      IF (info /= 0) ERROR STOP 'Compatibility test point characters failed'
      n = SIZE(members)
      ALLOCATE (point%rotations(3, 3, n), point%translations(3, n), point%multiplicity(SIZE(point%characters, 2)))
      ALLOCATE (point%antiunitary(n))
      point%antiunitary(:) = anti(members)
      point%k = k
      point%cell = cell
      point%spinful = spinful
      point%rotations(:, :, :) = rotations(:, :, members)
      point%translations(:, :) = tau(:, members)
      point%multiplicity(:) = NINT(REAL(point%characters(1, :), dp))
   END SUBROUTINE make_little_point

! **************************************************************************************************
!> \brief Ordinary D2, its quaternion spin factor, and arbitrary coset-gauge character covariance.
! **************************************************************************************************
   SUBROUTINE projective_character_tests()
      COMPLEX(dp)                                        :: beta(4), regular(4)
      COMPLEX(dp), ALLOCATABLE                           :: characters(:, :), factor(:, :), &
                                                            reference(:, :)
      INTEGER                                            :: g, h, i, j, rotations(3, 3, 4), status
      INTEGER, ALLOCATABLE                               :: mapping(:), multiplicity(:), &
                                                            products(:, :)
      REAL(dp)                                           :: cell(3, 3), residual, translations(3, 4)

      cell = 0.0_dp
      rotations = 0
      translations = 0.0_dp
      DO i = 1, 3
         cell(i, i) = 1.0_dp
         rotations(i, i, 1) = 1
         DO j = 1, 3
            rotations(i, i, j + 1) = -1
         END DO
         rotations(i, i, i + 1) = 1
      END DO
      CALL little_group_factor_system(cell, rotations, translations, [.FALSE., .FALSE., .FALSE., .FALSE.], &
                                      [0.0_dp, 0.0_dp, 0.0_dp], .FALSE., mapping, products, factor, 1.e-10_dp, status)
      IF (status /= 0) ERROR STOP 'D2 scalar factors failed'
      CALL projective_characters(products, factor, characters, 1.e-9_dp, residual, status)
      IF (status /= 0) ERROR STOP 'D2 characters failed'
      IF (SIZE(characters, 2) /= 4) ERROR STOP 'D2 scalar irreps incomplete'
      IF (MAXVAL(ABS(characters(1, :) - 1.0_dp)) > 1.e-10_dp) ERROR STOP 'D2 scalar irrep dimension'
      IF (MAXVAL(ABS(characters**2 - 1.0_dp)) > 1.e-10_dp) ERROR STOP 'D2 scalar characters not +/-1'
      reference = characters
      beta(1) = 1.0_dp
      DO g = 2, 4
         beta(g) = EXP(CMPLX(0.0_dp, 0.317_dp*g, dp))
      END DO
      DO h = 1, 4
         DO g = 1, 4
            factor(g, h) = beta(g)*beta(h)/beta(products(g, h))
         END DO
      END DO
      CALL projective_characters(products, factor, characters, 1.e-9_dp, residual, status)
      IF (status /= 0) ERROR STOP 'Complex coboundary characters failed'
      ALLOCATE (multiplicity(4))
      DO i = 1, 4
         CALL character_multiplicities(beta*reference(:, i), characters, multiplicity, 1.e-9_dp, status)
         IF (status /= 0 .OR. SUM(multiplicity) /= 1) ERROR STOP 'Character coboundary covariance failed'
      END DO
      DEALLOCATE (multiplicity)
      CALL little_group_factor_system(cell, rotations, translations, [.FALSE., .FALSE., .FALSE., .FALSE.], &
                                      [0.0_dp, 0.0_dp, 0.0_dp], .TRUE., mapping, products, factor, 1.e-10_dp, status)
      IF (status /= 0) ERROR STOP 'D2 spin factors failed'
      CALL projective_characters(products, factor, characters, 1.e-9_dp, residual, status)
      IF (status /= 0) ERROR STOP 'Quaternion projective characters failed'
      IF (SIZE(characters, 2) /= 1) ERROR STOP 'Quaternion spin irrep count'
      regular = [2.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
      IF (MAXVAL(ABS(characters(:, 1) - regular)) > 1.e-10_dp) ERROR STOP 'Quaternion spin character'
      ALLOCATE (multiplicity(1))
      CALL character_multiplicities(2*regular, characters, multiplicity, 1.e-9_dp, status)
      IF (status /= 0 .OR. multiplicity(1) /= 2) ERROR STOP 'Quaternion regular decomposition'
      factor(2, 2) = 0.0_dp
      CALL projective_characters(products, factor, characters, 1.e-9_dp, residual, status)
      IF (status == 0 .OR. ALLOCATED(characters)) ERROR STOP 'Malformed factor character table accepted'
      WRITE (io_unit, *) 'Unitary projective characters and gauge covariance tests passed.'
   END SUBROUTINE projective_character_tests

! **************************************************************************************************
!> \brief Analytic screw, centering, antiunitary and representative-gauge factor tests.
! **************************************************************************************************
   SUBROUTINE little_group_tests()
      COMPLEX(dp)                                        :: ao(2, 2), beta(2), bloch, &
                                                            sewing(4, 4, 2), spin(2, 2)
      COMPLEX(dp), ALLOCATABLE                           :: characters(:, :), factors(:, :), &
                                                            original(:, :)
      INTEGER                                            :: g, h, i, ident(3, 3), j, &
                                                            rotations(3, 3, 4), status
      INTEGER, ALLOCATABLE                               :: mapping(:), products(:, :)
      LOGICAL                                            :: anti(4)
      REAL(dp)                                           :: cell(3, 3), error, k(3), origin(3), &
                                                            shifts(3, 4)

      cell = 0.0_dp
      ident = 0
      DO i = 1, 3
         cell(i, i) = 1.0_dp
         ident(i, i) = 1
      END DO
      rotations(:, :, 1) = ident
      rotations(:, :, 2) = ident
      rotations(1, 1, 2) = -1
      rotations(2, 2, 2) = -1
      shifts = 0.0_dp
      shifts(3, 2) = 0.5_dp
      anti = .FALSE.
      k = [0.0_dp, 0.0_dp, 0.25_dp]
      CALL little_group_factor_system(cell, rotations(:, :, :2), shifts(:, :2), anti(:2), k, .FALSE., &
                                      mapping, products, factors, 1.e-10_dp, status)
      IF (status /= 0) ERROR STOP 'Scalar screw little group failed'
      IF (ANY(mapping /= [1, 2]) .OR. products(2, 2) /= 1) ERROR STOP 'Screw product table failed'
      bloch = CMPLX(0.0_dp, -1.0_dp, dp)
      IF (ABS(factors(2, 2) - bloch) > 1.e-12_dp) ERROR STOP 'Screw translation phase failed'
      CALL projective_characters(products, factors, characters, 1.e-10_dp, error, status)
      IF (status /= 0) ERROR STOP 'Nonsymmorphic screw characters failed'
      IF (SIZE(characters, 2) /= 2) ERROR STOP 'Nonsymmorphic screw irrep count'
      IF (MAXVAL(ABS(characters(2, :)**2 - bloch)) > 1.e-10_dp) ERROR STOP 'Screw eigenvalues wrong'
      CALL little_group_factor_system(cell, rotations(:, :, :2), shifts(:, :2), anti(:2), k, .TRUE., &
                                      mapping, products, factors, 1.e-10_dp, status)
      IF (status /= 0) ERROR STOP 'Spinor screw little group failed'
      IF (ABS(factors(2, 2) + bloch) > 1.e-12_dp) ERROR STOP 'Screw spin/translation phase failed'
      ! Independent two-site orbital/spin action: the second site crosses the z cell boundary.
      ao = 0.0_dp
      ao(2, 1) = 1.0_dp
      ao(1, 2) = bloch
      CALL spin_rotation(REAL(rotations(:, :, 2), dp), spin, 1.e-12_dp, status)
      sewing = 0.0_dp
      DO i = 1, 4
         sewing(i, i, 1) = 1.0_dp
      END DO
      DO j = 1, 2
         DO i = 1, 2
            sewing(2*i - 1:2*i, 2*j - 1:2*j, 2) = spin(i, j)*ao
         END DO
      END DO
      CALL check_sewing_product(sewing, products, factors, anti(:2), 1.e-12_dp, error, status)
      IF (status /= 0) ERROR STOP 'Atom/spin sewing disagrees with little-group factors'
      original = factors
      CALL little_group_factor_system(cell, rotations(:, :, :2), shifts(:, :2), anti(:2), &
                                      k + [2.0_dp, -1.0_dp, 3.0_dp], .TRUE., &
                                      mapping, products, factors, 1.e-10_dp, status)
      IF (status /= 0) ERROR STOP 'Reciprocal-shifted little group failed'
      IF (MAXVAL(ABS(factors - original)) > 1.e-12_dp) ERROR STOP 'Factors not reciprocal periodic'
      shifts(3, 2) = shifts(3, 2) + 1.0_dp
      beta = [CMPLX(1.0_dp, 0.0_dp, dp), bloch]
      CALL little_group_factor_system(cell, rotations(:, :, :2), shifts(:, :2), anti(:2), k, .TRUE., &
                                      mapping, products, factors, 1.e-10_dp, status)
      IF (status /= 0) ERROR STOP 'Integer-shifted coset failed'
      DO h = 1, 2
         DO g = 1, 2
            IF (ABS(factors(g, h) - original(g, h)*beta(g)*beta(h)/beta(products(g, h))) > 1.e-12_dp) THEN
               ERROR STOP 'Coset representative coboundary failed'
            END IF
         END DO
      END DO
      shifts(3, 2) = 0.5_dp
      origin = [0.137_dp, 0.231_dp, 0.317_dp]
      shifts(:, 2) = shifts(:, 2) + MATMUL(REAL(rotations(:, :, 2) - ident, dp), origin)
      CALL little_group_factor_system(cell, rotations(:, :, 2:1:-1), shifts(:, 2:1:-1), anti(:2), k, .TRUE., &
                                      mapping, products, factors, 1.e-10_dp, status)
      IF (status /= 0) ERROR STOP 'Origin-shifted screw failed'
      IF (ANY(mapping /= [2, 1])) ERROR STOP 'Identity not moved to first position'
      IF (MAXVAL(ABS(factors - original)) > 1.e-12_dp) ERROR STOP 'Origin shift changed factor system'
      k(1) = 0.173_dp
      CALL little_group_factor_system(cell, rotations(:, :, :2), shifts(:, :2), anti(:2), k, .TRUE., &
                                      mapping, products, factors, 1.e-10_dp, status)
      IF (status /= 0) ERROR STOP 'Generic k selection failed'
      IF (SIZE(mapping) /= 1) ERROR STOP 'Non-little operation retained'

      ! Centering has the same rotation for distinct cosets and must not be collapsed.
      rotations(:, :, 2) = ident
      shifts = 0.0_dp
      shifts(3, 2) = 0.5_dp
      k = [0.0_dp, 0.0_dp, 0.25_dp]
      CALL little_group_factor_system(cell, rotations(:, :, :2), shifts(:, :2), anti(:2), k, .FALSE., &
                                      mapping, products, factors, 1.e-10_dp, status)
      IF (status /= 0) ERROR STOP 'Centering cosets failed'
      IF (SIZE(mapping) /= 2 .OR. ABS(factors(2, 2) - bloch) > 1.e-12_dp) THEN
         ERROR STOP 'Centering translation lost'
      END IF
      ! A magnetic half translation has (translation*Theta)^2=-exp(-2*pi*i*k.L).
      anti(2) = .TRUE.
      k(3) = 0.5_dp
      CALL little_group_factor_system(cell, rotations(:, :, :2), shifts(:, :2), anti(:2), k, .TRUE., &
                                      mapping, products, factors, 1.e-10_dp, status)
      IF (status /= 0) ERROR STOP 'Antiferromagnetic half translation failed'
      IF (ABS(factors(2, 2) - 1.0_dp) > 1.e-12_dp) ERROR STOP 'Magnetic half-translation square failed'
      shifts = 0.0_dp
      CALL little_group_factor_system(cell, rotations(:, :, :2), shifts(:, :2), anti(:2), k, .TRUE., &
                                      mapping, products, factors, 1.e-10_dp, status)
      IF (status /= 0) ERROR STOP 'Physical time reversal failed'
      IF (ABS(factors(2, 2) + 1.0_dp) > 1.e-12_dp) ERROR STOP 'Physical Theta square failed'
      CALL little_group_factor_system(cell, rotations(:, :, :2), shifts(:, :2), anti(:2), k, .FALSE., &
                                      mapping, products, factors, 1.e-10_dp, status)
      IF (status /= 0) ERROR STOP 'Scalar time reversal failed'
      IF (MAXVAL(ABS(factors - 1.0_dp)) > 1.e-12_dp) ERROR STOP 'Scalar Theta square failed'
      ! Antiunitary coboundaries need conjugation even when their phases are complex.
      beta(1) = 1.0_dp
      beta(2) = EXP(CMPLX(0.0_dp, 0.371_dp, dp))
      original = factors
      DO h = 1, 2
         DO g = 1, 2
            IF (anti(g)) THEN
               factors(g, h) = original(g, h)*beta(g)*CONJG(beta(h))/beta(products(g, h))
            ELSE
               factors(g, h) = original(g, h)*beta(g)*beta(h)/beta(products(g, h))
            END IF
         END DO
      END DO
      CALL validate_factor_system(products, factors, anti(:2), 1.e-12_dp, error, status)
      IF (status /= 0) ERROR STOP 'Twisted coboundary rejected'
      factors(2, 2) = EXP(CMPLX(0.0_dp, 0.25_dp, dp))
      CALL validate_factor_system(products, factors, anti(:2), 1.e-12_dp, error, status)
      IF (status == 0) ERROR STOP 'Untwisted antiunitary cocycle accepted'

      anti = .FALSE.
      shifts(3, 2) = 0.3_dp
      CALL little_group_factor_system(cell, rotations(:, :, :2), shifts(:, :2), anti(:2), k, .TRUE., &
                                      mapping, products, factors, 1.e-10_dp, status)
      IF (status == 0 .OR. ALLOCATED(mapping)) ERROR STOP 'Nonclosed translation set accepted'
      shifts = 0.0_dp
      CALL little_group_factor_system(cell, rotations(:, :, :2), shifts(:, :2), anti(:2), k, .TRUE., &
                                      mapping, products, factors, 1.e-10_dp, status)
      IF (status == 0) ERROR STOP 'Duplicate operation accepted'
      rotations(1, 1, 2) = 2
      CALL little_group_factor_system(cell, rotations(:, :, :2), shifts(:, :2), anti(:2), k, .FALSE., &
                                      mapping, products, factors, 1.e-10_dp, status)
      IF (status == 0) ERROR STOP 'Noncrystallographic operation accepted'
      cell(1, 1) = 0.0_dp
      CALL little_group_factor_system(cell, rotations(:, :, :1), shifts(:, :1), anti(:1), k, .FALSE., &
                                      mapping, products, factors, 1.e-10_dp, status)
      IF (status == 0) ERROR STOP 'Singular cell accepted'
      cell(1, 1) = ieee_value(0.0_dp, ieee_quiet_nan)
      CALL little_group_factor_system(cell, rotations(:, :, :1), shifts(:, :1), anti(:1), k, .FALSE., &
                                      mapping, products, factors, 1.e-10_dp, status)
      IF (status == 0) ERROR STOP 'Nonfinite cell accepted'
      WRITE (io_unit, *) 'Nonsymmorphic, spin and antiunitary little-group algebra tests passed.'
   END SUBROUTINE little_group_tests

! **************************************************************************************************
!> \brief Nonprimitive half-translation: two applications give the full-cell Bloch phase.
! **************************************************************************************************
   SUBROUTINE cell_gauge_tests()
      COMPLEX(dp)                                        :: action(2, 2), expected(2, 2), h(2, 2), &
                                                            translated(2, 2)
      INTEGER                                            :: i, status
      LOGICAL                                            :: antiunitary
      REAL(dp)                                           :: k(3), rotation(3, 3)
      TYPE(kpoint_type), POINTER                         :: points

      NULLIFY (points)
      CALL kpoint_create(points)
      points%nkp = 1
      ALLOCATE (points%xkp(3, 1), points%atype(2), points%ibrot(1), points%kind_rotmat(1, 1))
      k = [0.25_dp, 0.0_dp, 0.0_dp]
      points%xkp(:, 1) = k
      points%atype = 1
      points%ibrot = 1
      ALLOCATE (points%kind_rotmat(1, 1)%rmat(1, 1))
      points%kind_rotmat(1, 1)%rmat = 1.0_dp
      ALLOCATE (points%kp_sym(1))
      CALL kpoint_sym_create(points%kp_sym(1)%kpoint_sym)
      ASSOCIATE (sym => points%kp_sym(1)%kpoint_sym)
         sym%apply_symmetry = .TRUE.
         sym%nwred = 1
         ALLOCATE (sym%rot(3, 3, 1), sym%rotp(1), sym%xkp(3, 1), sym%f0(2, 1), &
                   sym%fcell_gauge(3, 2, 1), sym%kgphase(2, 1))
         sym%rot = 0.0_dp
         DO i = 1, 3
            sym%rot(i, i, 1) = 1.0_dp
         END DO
         sym%rotp = 1
         sym%xkp(:, 1) = k
         sym%f0(:, 1) = [2, 1]
         sym%fcell_gauge = 0
         sym%fcell_gauge(1, 2, 1) = 1
         sym%kgphase = 0.37_dp
         CALL kpoint_cell_gauge_action(points, 1, 1, k, action, rotation, antiunitary, status)
         IF (status /= 0 .OR. antiunitary) ERROR STOP 'Half-translation construction failed'
         expected = 0.0_dp
         expected(1, 2) = CMPLX(0.0_dp, -1.0_dp, dp)
         expected(2, 1) = 1.0_dp
         IF (MAXVAL(ABS(action - expected)) > 1.e-12_dp) ERROR STOP 'Cell-gauge image phase failed'
         translated = MATMUL(action, action)
         expected = 0.0_dp
         expected(1, 1) = CMPLX(0.0_dp, -1.0_dp, dp)
         expected(2, 2) = expected(1, 1)
         IF (MAXVAL(ABS(translated - expected)) > 1.e-12_dp) ERROR STOP 'Nonsymmorphic factor failed'
         h = 0.0_dp
         h(1, 2) = CMPLX(1.0_dp, -1.0_dp, dp)
         h(2, 1) = CONJG(h(1, 2))
         IF (MAXVAL(ABS(MATMUL(h, action) - MATMUL(action, h))) > 1.e-12_dp) THEN
            ERROR STOP 'Half-translation Hamiltonian covariance failed'
         END IF
         CALL kpoint_cell_gauge_action(points, 1, 1, k + [1.0_dp, 0.0_dp, 0.0_dp], &
                                       translated, rotation, antiunitary, status)
         IF (status /= 0 .OR. MAXVAL(ABS(action - translated)) > 1.e-12_dp) THEN
            ERROR STOP 'Reciprocal cell-gauge periodicity failed'
         END IF
         sym%f0(:, 1) = [1, 1]
         CALL kpoint_cell_gauge_action(points, 1, 1, k, action, rotation, antiunitary, status)
         IF (status == 0) ERROR STOP 'Nonbijective atom map accepted'
      END ASSOCIATE
      CALL kpoint_release(points)
   END SUBROUTINE cell_gauge_tests

! **************************************************************************************************
!> \brief All 48 cubic polar operations, double-group products, screw phases and Kramers conjugation.
! **************************************************************************************************
   SUBROUTINE spinor_tests()
      COMPLEX(dp) :: ao(2, 2), factor(2, 2), factors(48, 48), frame(4, 2), ident(2, 2), &
         image(4, 2), pauli(2, 2, 3), sewing(2, 2, 2), spin(2, 2, 48), twice(4, 2), work(2, 2)
      INTEGER                                            :: a, b, bits, count, g, gh, h, p1, p2, p3, &
                                                            PRODUCT(2, 2), products(48, 48), status
      REAL(dp)                                           :: det, error, pi, r(3, 3), &
                                                            rotations(3, 3, 48)

      ident = 0.0_dp
      ident(1, 1) = 1.0_dp
      ident(2, 2) = 1.0_dp
      pauli = 0.0_dp
      pauli(1, 2, 1) = 1.0_dp
      pauli(2, 1, 1) = 1.0_dp
      pauli(1, 2, 2) = CMPLX(0.0_dp, -1.0_dp, dp)
      pauli(2, 1, 2) = CMPLX(0.0_dp, 1.0_dp, dp)
      pauli(1, 1, 3) = 1.0_dp
      pauli(2, 2, 3) = -1.0_dp
      count = 0
      DO p1 = 1, 3
         DO p2 = 1, 3
            IF (p2 == p1) CYCLE
            p3 = 6 - p1 - p2
            DO bits = 0, 7
               count = count + 1
               r = 0.0_dp
               r(p1, 1) = REAL(1 - 2*IBITS(bits, 0, 1), dp)
               r(p2, 2) = REAL(1 - 2*IBITS(bits, 1, 1), dp)
               r(p3, 3) = REAL(1 - 2*IBITS(bits, 2, 1), dp)
               det = r(1, 1)*(r(2, 2)*r(3, 3) - r(2, 3)*r(3, 2)) &
                     - r(1, 2)*(r(2, 1)*r(3, 3) - r(2, 3)*r(3, 1)) &
                     + r(1, 3)*(r(2, 1)*r(3, 2) - r(2, 2)*r(3, 1))
               rotations(:, :, count) = r
               CALL spin_rotation(r, spin(:, :, count), 1.e-12_dp, status)
               IF (status /= 0) ERROR STOP 'Cubic spin rotation rejected'
               DO a = 1, 3
                  work = MATMUL(spin(:, :, count), MATMUL(pauli(:, :, a), CONJG(TRANSPOSE(spin(:, :, count)))))
                  DO b = 1, 3
                     work = work - det*r(b, a)*pauli(:, :, b)
                  END DO
                  IF (MAXVAL(ABS(work)) > 1.e-12_dp) ERROR STOP 'Axial Pauli covariance failed'
               END DO
            END DO
         END DO
      END DO
      DO g = 1, 48
         DO h = 1, 48
            r = MATMUL(rotations(:, :, g), rotations(:, :, h))
            DO gh = 1, 48
               IF (MAXVAL(ABS(r - rotations(:, :, gh))) < 1.e-12_dp) EXIT
            END DO
            IF (gh > 48) ERROR STOP 'Cubic product missing'
            products(g, h) = gh
            work = MATMUL(spin(:, :, g), MATMUL(spin(:, :, h), CONJG(TRANSPOSE(spin(:, :, gh)))))
            factors(g, h) = (work(1, 1) + work(2, 2))/2.0_dp
         END DO
      END DO
      CALL check_sewing_product(spin, products, factors, SPREAD(.FALSE., 1, 48), 1.e-12_dp, error, status)
      IF (status /= 0) ERROR STOP 'Double-group products failed'
      IF (MAXVAL(ABS(factors - 1.0_dp)) < 1.0_dp) ERROR STOP 'Spin double cover lost'
      CALL spin_lift_frame_tests(NINT(rotations))
      ! At a fixed k a nonsymmorphic operation can have a non-real Bloch factor.
      pi = ACOS(-1.0_dp)
      r = 0.0_dp
      r(1, 1) = COS(pi/3.0_dp)
      r(2, 2) = r(1, 1)
      r(1, 2) = -SIN(pi/3.0_dp)
      r(2, 1) = -r(1, 2)
      r(3, 3) = 1.0_dp
      CALL spin_rotation(r, sewing(:, :, 2), 1.e-12_dp, status)
      IF (status /= 0) ERROR STOP 'Hexagonal spin rotation rejected'
      work = ident
      DO a = 1, 6
         work = MATMUL(work, sewing(:, :, 2))
      END DO
      IF (MAXVAL(ABS(work + ident)) > 1.e-12_dp) ERROR STOP 'Hexagonal double cover failed'
      sewing(:, :, 1) = ident
      sewing(:, :, 2) = CMPLX(0.0_dp, -1.0_dp, dp)*pauli(:, :, 3)*EXP(CMPLX(0.0_dp, -pi/4, dp))
      product = RESHAPE([1, 2, 2, 1], [2, 2])
      factor = 1.0_dp
      factor(2, 2) = CMPLX(0.0_dp, 1.0_dp, dp)
      CALL check_sewing_product(sewing, product, factor, [.FALSE., .FALSE.], 1.e-12_dp, error, status)
      IF (status /= 0) ERROR STOP 'Nonsymmorphic Bloch factor failed'
      DO a = 0, 3
         CALL check_sewing_product(sewing, product, factor, [.FALSE., .FALSE.], 1.e-12_dp, error, status, [a, 4])
         IF (status /= 0 .OR. error > 1.e-12_dp) ERROR STOP 'Partitioned Bloch factor failed'
      END DO
      factor(2, 2) = -1.0_dp
      CALL check_sewing_product(sewing, product, factor, [.FALSE., .FALSE.], 1.e-12_dp, error, status, [0, 4])
      IF (status /= 0) ERROR STOP 'Untouched product partition failed'
      CALL check_sewing_product(sewing, product, factor, [.FALSE., .FALSE.], 1.e-12_dp, error, status, [1, 4])
      IF (status == 0) ERROR STOP 'Nonroot product error was missed'
      CALL check_sewing_product(sewing, product, factor, [.FALSE., .FALSE.], 1.e-12_dp, error, status)
      IF (status == 0) ERROR STOP 'Missing Bloch phase accepted'
      sewing(:, :, 2) = CMPLX(0.0_dp, 1.0_dp, dp)*pauli(:, :, 2)
      CALL check_sewing_product(sewing, product, factor, [.FALSE., .TRUE.], 1.e-12_dp, error, status)
      IF (status /= 0) ERROR STOP 'Kramers factor failed'
      DO a = 0, 3
         CALL check_sewing_product(sewing, product, factor, [.FALSE., .TRUE.], 1.e-12_dp, error, status, [a, 4])
         IF (status /= 0 .OR. error > 1.e-12_dp) ERROR STOP 'Partitioned antiunitary product failed'
      END DO
      CALL check_sewing_product(sewing, product, factor, [.FALSE., .TRUE.], 1.e-12_dp, error, status, [0, 0])
      IF (status == 0) ERROR STOP 'Empty worker partition accepted'
      CALL check_sewing_product(sewing, product, factor, [.FALSE., .TRUE.], 1.e-12_dp, error, status, [2, 2])
      IF (status == 0) ERROR STOP 'Invalid worker index accepted'
      frame = RESHAPE([CMPLX(1.0_dp, 0.2_dp, dp), CMPLX(-0.2_dp, 0.4_dp, dp), &
                       CMPLX(0.5_dp, 0.6_dp, dp), CMPLX(0.3_dp, -0.1_dp, dp), &
                       CMPLX(0.4_dp, -0.8_dp, dp), CMPLX(0.5_dp, 0.7_dp, dp), &
                       CMPLX(-0.3_dp, 0.2_dp, dp), CMPLX(0.8_dp, -0.6_dp, dp)], [4, 2])
      ao = ident
      CALL transform_spinor_frame(ao, ident, frame, .TRUE., image, 1.e-12_dp, status)
      IF (status /= 0) ERROR STOP 'Spinor time reversal failed'
      CALL transform_spinor_frame(ao, ident, image, .TRUE., twice, 1.e-12_dp, status)
      IF (status /= 0 .OR. MAXVAL(ABS(twice + frame)) > 1.e-12_dp) ERROR STOP 'Theta squared failed'
      r = 0.0_dp
      CALL spin_rotation(r, work, 1.e-12_dp, status)
      IF (status == 0) ERROR STOP 'Invalid Cartesian rotation accepted'
   END SUBROUTINE spinor_tests

! **************************************************************************************************
!> \brief Apply stored double-group lifts to complex frames in a primitive cubic cell.
!> \param cubic Complete cubic polar operations, including inversion
! **************************************************************************************************
   SUBROUTINE spin_lift_frame_tests(cubic)
      INTEGER, INTENT(IN)                                :: cubic(3, 3, 48)

      COMPLEX(dp)                                        :: frame(2, 2), image(2, 2), orbital(1, 1), &
                                                            spin(2, 2)
      COMPLEX(dp), ALLOCATABLE                           :: factors(:, :), lifts(:, :, :), &
                                                            sewing(:, :, :)
      INTEGER                                            :: g, h, ns, pass, status, time_reversal
      INTEGER, ALLOCATABLE                               :: mapping(:), PRODUCT(:, :), &
                                                            rotations(:, :, :)
      LOGICAL, ALLOCATABLE                               :: anti(:)
      REAL(dp)                                           :: cell(3, 3), error, inverse(3, 3), &
                                                            primitive(3, 3), scale

      primitive = RESHAPE([1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &
                           0.5_dp, 0.5_dp, 0.5_dp], [3, 3])
      inverse = inv_3x3(primitive)
      cell = primitive
      cell(1, :) = cell(1, :)*(1 + EPSILON(1.0_dp))
      cell(3, 2) = EPSILON(1.0_dp)/4
      orbital = 1
      frame(:, 1) = [CMPLX(1.0_dp, 0.0_dp, dp), CMPLX(0.0_dp, 1.0_dp, dp)]/SQRT(2.0_dp)
      frame(:, 2) = [CMPLX(0.0_dp, 1.0_dp, dp), CMPLX(1.0_dp, 0.0_dp, dp)]/SQRT(2.0_dp)
      DO time_reversal = 0, 1
         ns = 48*(time_reversal + 1)
         ALLOCATE (rotations(3, 3, ns), anti(ns), sewing(2, 2, ns))
         anti = .FALSE.
         DO g = 1, 48
            rotations(:, :, g) = NINT(MATMUL(MATMUL(inverse, REAL(cubic(:, :, g), dp)), primitive))
         END DO
         IF (time_reversal == 1) THEN
            rotations(:, :, 49:) = rotations(:, :, :48)
            anti(49:) = .TRUE.
         END IF
         DO pass = 1, 3
            scale = 4.2_dp*10.0_dp**(100*(pass - 2))
            CALL little_group_factor_system(scale*cell, rotations, SPREAD([0.0_dp, 0.0_dp, 0.0_dp], 2, ns), &
                                            anti, [0.0_dp, 0.0_dp, 0.0_dp], .TRUE., mapping, product, &
                                            factors, 1.e-10_dp, status, spin_lifts=lifts)
            IF (status /= 0) ERROR STOP 'Primitive cubic spin-lift algebra failed'
            DO g = 1, ns
               h = mapping(g)
               spin = CONJG(lifts(:, :, g))
               IF (anti(h)) THEN
                  spin(:, 1) = CONJG(lifts(:, 2, g))
                  spin(:, 2) = -CONJG(lifts(:, 1, g))
               END IF
               CALL transform_spinor_frame(orbital, spin, frame, anti(h), image, 1.e-10_dp, status)
               IF (status /= 0) ERROR STOP 'Stored spin-lift frame transformation failed'
               sewing(:, :, g) = MATMUL(CONJG(TRANSPOSE(frame)), image)
            END DO
            CALL check_sewing_product(sewing, product, factors, anti(mapping), 1.e-10_dp, error, status)
            IF (status /= 0) ERROR STOP 'Stored spin lifts violate the unitary/antiunitary frame algebra'
         END DO
         DEALLOCATE (rotations, anti, sewing)
      END DO
   END SUBROUTINE spin_lift_frame_tests

! **************************************************************************************************
!> \brief Compare the exact solver against independent enumeration of rank-two atomic sums.
!> \param DIMENSION Dimension of the inversion subgroup
! **************************************************************************************************
   SUBROUTINE exhaustive_signatures(DIMENSION)
      INTEGER, INTENT(IN)                                :: dimension

      INTEGER :: code, digit, i, ierr, indicator, indices(3), j, n, nu, x, y, &
         values(2**DIMENSION), signatures(2**DIMENSION, 2**(DIMENSION + 1)), &
         reconstructed(2**DIMENSION), coefficients(2, 2**DIMENSION)
      LOGICAL :: indicators_zero, nonnegative_ok, possible(0:3**(2**DIMENSION) - 1), signed_ok

      n = 2**DIMENSION
      DO x = 0, n - 1
         DO i = 0, n - 1
            signatures(i + 1, 2*x + 1) = (1 - (-1)**POPCNT(IAND(x, i)))/2
            signatures(i + 1, 2*x + 2) = 1 - signatures(i + 1, 2*x + 1)
         END DO
      END DO
      possible = .FALSE.
      DO x = 1, 2*n
         DO y = 1, 2*n
            code = SUM((signatures(:, x) + signatures(:, y))*[(3**i, i=0, n - 1)])
            possible(code) = .TRUE.
         END DO
      END DO
      DO code = 0, 3**n - 1
         digit = code
         DO i = 1, n
            values(i) = MOD(digit, 3)
            digit = digit/3
         END DO
         CALL inversion_ebr_signature(values, 2, signed_ok, nonnegative_ok, coefficients, ierr)
         IF (ierr /= 0 .OR. (nonnegative_ok .NEQV. possible(code))) ERROR STOP 'Enumerated EBR mismatch'
         CALL inversion_indicators(values, 2, DIMENSION, nu, indices, indicator, ierr)
         IF (ierr /= 0) ERROR STOP 'Enumerated indicator failure'
         indicators_zero = nu == 0
         IF (DIMENSION == 3) indicators_zero = indicator == 0 .AND. ALL(indices == 0)
         IF (signed_ok .NEQV. indicators_zero) ERROR STOP 'Indicator/atomic-lattice mismatch'
         IF (.NOT. nonnegative_ok) CYCLE
         reconstructed = 0
         DO j = 1, n
            reconstructed = reconstructed + coefficients(1, j)*signatures(:, 2*j - 1) + &
                            coefficients(2, j)*signatures(:, 2*j)
         END DO
         IF (ANY(reconstructed /= values) .OR. SUM(coefficients) /= 2) ERROR STOP 'EBR reconstruction failed'
      END DO
   END SUBROUTINE exhaustive_signatures

! **************************************************************************************************
!> \brief Complex gauges, AO metric, Bloch phase, Kramers sewing and antiunitary links.
! **************************************************************************************************
   SUBROUTINE sewing_tests()
      COMPLEX(dp) :: answer(2, 2), b(2, 2), basis(2, 2), exact(2, 2), image(2, 2), jtr(2, 2), &
         link(2, 2), right(2, 2), s(2, 2), single(1, 1), source(2, 2), TARGET(2, 2), u(2, 2), &
         v(2, 2), w(2, 2)
      INTEGER                                            :: i, info
      REAL(dp)                                           :: errors(4)

      s = 0.0_dp
      basis = 0.0_dp
      s(1, 1) = 2.0_dp
      s(2, 2) = 0.5_dp
      basis(1, 1) = 1/SQRT(2.0_dp)
      basis(2, 2) = SQRT(2.0_dp)
      u = RESHAPE([CMPLX(1, 0, dp), CMPLX(0, 1, dp), CMPLX(0, 1, dp), CMPLX(1, 0, dp)], [2, 2])/SQRT(2.0_dp)
      v = RESHAPE([0.8_dp, -0.6_dp, 0.6_dp, 0.8_dp], [2, 2])
      w = CMPLX(0.0_dp, -1.0_dp, dp)*v
      TARGET = MATMUL(basis, u)
      image = MATMUL(basis, w)
      CALL subspace_sewing(s, TARGET, image, b, 1.e-10_dp, errors, info)
      IF (info /= 0 .OR. MAXVAL(errors) > 1.e-12_dp) ERROR STOP 'Complex metric sewing failed'
      exact = MATMUL(CONJG(TRANSPOSE(u)), w)
      IF (MAXVAL(ABS(b - exact)) > 1.e-12_dp) ERROR STOP 'Degenerate/Bloch-phase sewing gauge failed'
      CALL subspace_sewing(s, TARGET(:, :1), image(:, :1), single, 1.e-10_dp, errors, info)
      IF (info /= -2) ERROR STOP 'Nonclosed sewing subspace accepted'
      CALL subspace_sewing(s, TARGET, 2*image, b, 1.e-10_dp, errors, info)
      IF (info /= -2) ERROR STOP 'Misnormalized transformed frame accepted'

      ! Two spin states with Theta=J K; arbitrary and DIFFERENT endpoint gauges.
      s = 0.0_dp
      s(1, 1) = 1.0_dp
      s(2, 2) = 1.0_dp
      jtr = RESHAPE([0.0_dp, -1.0_dp, 1.0_dp, 0.0_dp], [2, 2])
      source = u
      TARGET = v
      image = MATMUL(jtr, CONJG(source))
      CALL subspace_sewing(s(:1, :1), TARGET, image, b, 1.e-10_dp, errors, info)
      IF (info /= 0) ERROR STOP 'Spinor antiunitary sewing failed'
      CALL subspace_sewing(s(:1, :1), source, image, right, 1.e-10_dp, errors, info)
      IF (info /= 0) ERROR STOP 'Kramers self-sewing failed'
      IF (MAXVAL(ABS(right + TRANSPOSE(right))) > 1.e-12_dp) ERROR STOP 'Kramers sewing is not skew'
      exact = MATMUL(right, CONJG(right))
      DO i = 1, 2
         exact(i, i) = exact(i, i) + 1.0_dp
      END DO
      IF (MAXVAL(ABS(exact)) > 1.e-12_dp) ERROR STOP 'Spinful time reversal squared is not minus one'
      image = MATMUL(jtr, CONJG(v))
      CALL subspace_sewing(s(:1, :1), u, image, right, 1.e-10_dp, errors, info)
      IF (info /= 0) ERROR STOP 'Second endpoint antiunitary sewing failed'
      link = MATMUL(CONJG(TRANSPOSE(u)), v)
      CALL transform_sewing_link(link, b, right, .TRUE., answer, 1.e-10_dp, info)
      exact = MATMUL(CONJG(TRANSPOSE(v)), u)
      IF (info /= 0 .OR. MAXVAL(ABS(answer - exact)) > 1.e-12_dp) ERROR STOP 'Antiunitary link covariance failed'
      CALL subspace_sewing(s(:1, :1), v, u, b, 1.e-10_dp, errors, info)
      IF (info /= 0) ERROR STOP 'Unitary left sewing failed'
      CALL subspace_sewing(s(:1, :1), u, v, right, 1.e-10_dp, errors, info)
      IF (info /= 0) ERROR STOP 'Unitary right sewing failed'
      CALL transform_sewing_link(link, b, right, .FALSE., answer, 1.e-10_dp, info)
      IF (info /= 0 .OR. MAXVAL(ABS(answer - exact)) > 1.e-12_dp) ERROR STOP 'Unitary link covariance failed'
      CALL transform_sewing_link(link, 2*b, right, .FALSE., answer, 1.e-10_dp, info)
      IF (info /= -1) ERROR STOP 'Nonunitary link sewing accepted'
   END SUBROUTINE sewing_tests
END PROGRAM topology_symmetry_unittest
