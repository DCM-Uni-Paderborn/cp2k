!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
PROGRAM topology_character_unittest
   USE ieee_arithmetic,                 ONLY: ieee_quiet_nan,&
                                              ieee_value
   USE kinds,                           ONLY: dp
   USE mathlib,                         ONLY: inv_3x3
   USE topology_character_match,        ONLY: character_table_match
   USE topology_corepresentations,      ONLY: corepresentation_table_type,&
                                              projective_corepresentations
   USE topology_little_group,           ONLY: little_group_factor_system,&
                                              projective_characters

   IMPLICIT NONE
   CALL analytic_tests()
   CALL antiunitary_tests()
   CALL screw_lift_tests()
   PRINT *, 'Character correspondence tests passed.'
CONTAINS

! **************************************************************************************************
!> \brief Analytic ordinary/projective C3, label-changing character twists, and malformed input.
! **************************************************************************************************
   SUBROUTINE analytic_tests()
      COMPLEX(dp)                                        :: changed(3, 3), factor(3, 3), mix, &
                                                            other(3, 3), phase(3), ref(3, 3)
      COMPLEX(dp), ALLOCATABLE                           :: table(:, :)
      INTEGER                                            :: g, h, i, j, PRODUCT(3, 3), spin, status
      INTEGER, ALLOCATABLE                               :: mapping(:), original(:)
      REAL(dp)                                           :: angle, error

      DO spin = 0, 1
         DO h = 1, 3
            DO g = 1, 3
               PRODUCT(g, h) = MODULO(g + h - 2, 3) + 1
               factor(g, h) = 1
               IF (spin == 1 .AND. g + h - 2 >= 3) factor(g, h) = -1
            END DO
         END DO
         DO j = 1, 3
            DO g = 1, 3
               angle = ACOS(-1.0_dp)*(2*(j - 1) + spin)*(g - 1)/3
               ref(g, j) = CMPLX(COS(angle), SIN(angle), dp)
            END DO
         END DO
         CALL projective_characters(product, factor, table, 1.e-10_dp, error, status)
         IF (status /= 0) ERROR STOP 'Analytic C3 table failed'
         phase = 1
         CALL character_table_match(product, factor, [.FALSE., .FALSE., .FALSE.], table, &
                                    product, factor, [.FALSE., .FALSE., .FALSE.], ref, [1, 2, 3], &
                                    phase, mapping, 1.e-10_dp, error, status)
         IF (status /= 0) ERROR STOP 'Equal-dimensional analytic irreps not resolved'
         ALLOCATE (original(SIZE(mapping)))
         original(:) = mapping
         ! A one-dimensional twist leaves the factor system unchanged but permutes named irreps.
         DO g = 1, 3
            angle = 2*ACOS(-1.0_dp)*(g - 1)/3
            phase(g) = CMPLX(COS(angle), SIN(angle), dp)
         END DO
         CALL character_table_match(product, factor, [.FALSE., .FALSE., .FALSE.], table, &
                                    product, factor, [.FALSE., .FALSE., .FALSE.], ref, [1, 2, 3], &
                                    phase, mapping, 1.e-10_dp, error, status)
         IF (status /= 0) ERROR STOP 'Explicit character twist failed'
         IF (ANY(mapping == original)) ERROR STOP 'A factor system was mistaken for a unique label gauge'
         DO g = 1, 3
            angle = (g - 1)**2*SQRT(2.0_dp)
            phase(g) = CMPLX(COS(angle), SIN(angle), dp)
            changed(g, :) = phase(g)*ref(g, [3, 1, 2])
         END DO
         DO h = 1, 3
            DO g = 1, 3
               other(g, h) = phase(g)*phase(h)*CONJG(phase(PRODUCT(g, h)))*factor(g, h)
            END DO
         END DO
         CALL character_table_match(product, other, [.FALSE., .FALSE., .FALSE.], changed, &
                                    product, factor, [.FALSE., .FALSE., .FALSE.], ref, [1, 2, 3], &
                                    phase, mapping, 1.e-10_dp, error, status)
         IF (status /= 0) ERROR STOP 'Complex projective character rephasing failed'
         IF (ANY(mapping /= [3, 1, 2])) ERROR STOP 'Reference labels changed under reordering'
         DO i = 1, 6
            changed = ref
            phase = 1
            SELECT CASE (i)
            CASE (1)
               phase(2) = CMPLX(0, 1, dp)
            CASE (2)
               changed(:, 2) = changed(:, 1)
            CASE (3)
               changed(2, 1) = changed(2, 1) + 0.1_dp
            CASE (4)
               phase(2) = CMPLX(ieee_value(0.0_dp, ieee_quiet_nan), 0, dp)
            CASE (5)
               phase(2) = CMPLX(HUGE(1.0_dp), HUGE(1.0_dp), dp)
            CASE (6)
               changed(2, 1) = CMPLX(HUGE(1.0_dp), HUGE(1.0_dp), dp)
            END SELECT
            CALL character_table_match(product, factor, [.FALSE., .FALSE., .FALSE.], table, &
                                       product, factor, [.FALSE., .FALSE., .FALSE.], changed, [1, 2, 3], &
                                       phase, mapping, 1.e-10_dp, error, status)
            IF (status == 0 .OR. ALLOCATED(mapping)) ERROR STOP 'Invalid reference correspondence accepted'
         END DO
         ! This unitary column mixture preserves dimensions and orthogonality, but is not a character table.
         mix = CMPLX(COS(0.7_dp), SIN(0.7_dp), dp)
         changed = ref
         changed(:, 2) = 0.5_dp*((1 + mix)*ref(:, 2) + (1 - mix)*ref(:, 3))
         changed(:, 3) = 0.5_dp*((1 - mix)*ref(:, 2) + (1 + mix)*ref(:, 3))
         phase = 1
         CALL character_table_match(product, factor, [.FALSE., .FALSE., .FALSE.], changed, &
                                    product, factor, [.FALSE., .FALSE., .FALSE.], changed, [1, 2, 3], &
                                    phase, mapping, 1.e-10_dp, error, status)
         IF (status == 0 .OR. ALLOCATED(mapping)) ERROR STOP 'Orthogonality alone was accepted as a character certificate'
         DEALLOCATE (original)
      END DO
      PRINT *, 'Analytic C3, explicit lift conventions and malformed-reference tests passed.'
   END SUBROUTINE analytic_tests

! **************************************************************************************************
!> \brief Wigner paired/doubled restrictions with interleaved and rephased antiunitary operations.
! **************************************************************************************************
   SUBROUTINE antiunitary_tests()
      COMPLEX(dp)                                        :: composed(2, 2), phases(6), second
      COMPLEX(dp), ALLOCATABLE                           :: factor(:, :), lifts(:, :, :), &
                                                            other(:, :), ref(:, :), table(:, :), &
                                                            unit(:, :)
      INTEGER                                            :: g, h, i, inverse(6), order(6), &
                                                            r(3, 3, 6), spin, status
      INTEGER, ALLOCATABLE                               :: changed(:, :), mapping(:), &
                                                            operations(:), PRODUCT(:, :)
      LOGICAL                                            :: anti(6)
      REAL(dp)                                           :: angle, cell(3, 3), error, tau(3, 6)
      TYPE(corepresentation_table_type)                  :: coreps

      cell = 0
      cell(1, 1) = 1
      cell(1, 2) = -0.5_dp
      cell(2, 2) = SQRT(3.0_dp)/2
      cell(3, 3) = 1
      r = 0
      DO i = 1, 3
         r(i, i, 1) = 1
      END DO
      r(:, :, 2) = RESHAPE([0, 1, 0, -1, -1, 0, 0, 0, 1], [3, 3])
      r(:, :, 3) = MATMUL(r(:, :, 2), r(:, :, 2))
      r(:, :, 4:) = r(:, :, :3)
      anti = [.FALSE., .FALSE., .FALSE., .TRUE., .TRUE., .TRUE.]
      order = [1, 4, 3, 6, 2, 5]
      DO i = 1, 6
         inverse(order(i)) = i
         angle = (i - 1)**2*SQRT(3.0_dp)
         phases(i) = CMPLX(COS(angle), SIN(angle), dp)
      END DO
      tau = 0
      DO spin = 0, 1
         CALL little_group_factor_system(cell, r, tau, anti, [0.0_dp, 0.0_dp, 0.0_dp], &
                                         spin == 1, operations, product, factor, 1.e-10_dp, status, spin_lifts=lifts)
         IF (status /= 0) ERROR STOP 'Grey C3 algebra failed'
         IF (spin == 1) THEN
            DO h = 1, 6
               DO g = 1, 6
                  IF (anti(g)) THEN
                     composed = MATMUL(lifts(:, :, g), CONJG(lifts(:, :, h)))
                  ELSE
                     composed = MATMUL(lifts(:, :, g), lifts(:, :, h))
                  END IF
                  IF (MAXVAL(ABS(composed - factor(g, h)*lifts(:, :, PRODUCT(g, h)))) > 1.e-10_dp) THEN
                     ERROR STOP 'Returned spin lifts do not realize the twisted factor'
                  END IF
               END DO
            END DO
         END IF
         CALL projective_characters(PRODUCT(:3, :3), factor(:3, :3), unit, 1.e-10_dp, error, status)
         IF (status /= 0) ERROR STOP 'Grey C3 unitary table failed'
         CALL projective_corepresentations(product, factor, anti, unit, coreps, 1.e-10_dp, error, status)
         IF (status /= 0) ERROR STOP 'Grey C3 corepresentations failed'
         IF (.NOT. ANY(coreps%indicator == 0)) ERROR STOP 'Paired reference case missing'
         IF (spin == 1 .AND. .NOT. ANY(coreps%indicator == -1)) ERROR STOP 'Doubled reference case missing'
         ALLOCATE (ref(SIZE(unit, 1), SIZE(coreps%restriction, 2)))
         ref(:, :) = MATMUL(unit, CMPLX(coreps%restriction, 0, dp))
         ALLOCATE (other(6, 6), changed(6, 6), table(3, SIZE(ref, 2)))
         DO h = 1, 6
            DO g = 1, 6
               changed(g, h) = inverse(PRODUCT(order(g), order(h)))
               second = phases(h)
               IF (anti(order(g))) second = CONJG(second)
               other(g, h) = phases(g)*second*CONJG(phases(changed(g, h)))*factor(order(g), order(h))
            END DO
         END DO
         i = 0
         DO g = 1, 6
            IF (anti(order(g))) CYCLE
            i = i + 1
            table(i, :) = phases(g)*ref(order(g), :)
         END DO
         CALL character_table_match(changed, other, anti(order), table, product, factor, anti, ref, order, &
                                    phases, mapping, 1.e-10_dp, error, status)
         IF (status /= 0) ERROR STOP 'Twisted antiunitary character correspondence failed'
         IF (ANY(mapping /= [(i, i=1, SIZE(ref, 2))])) ERROR STOP 'Wigner column correspondence wrong'
         DEALLOCATE (other, changed, table, ref)
      END DO
      PRINT *, 'Paired/doubled corepresentations and interleaved antiunitary gauges passed.'
   END SUBROUTINE antiunitary_tests

! **************************************************************************************************
!> \brief Analytic body-centered three-screw boundary characters under rotated and rescaled cells.
! **************************************************************************************************
   SUBROUTINE screw_lift_tests()
      COMPLEX(dp)                                        :: frame(2, 2), imaginary, phases(4), &
                                                            ref(4, 4), ref_factor(4, 4), &
                                                            ref_spin(2, 2, 4), sigma(2, 2, 3)
      COMPLEX(dp), ALLOCATABLE                           :: factor(:, :), lifts(:, :, :), table(:, :)
      INTEGER                                            :: g, h, i, j, order(4), pass, q, &
                                                            ref_product(4, 4), rotations(3, 3, 4), &
                                                            signs(3), status
      INTEGER, ALLOCATABLE                               :: mapping(:), operations(:), PRODUCT(:, :)
      REAL(dp) :: angle, cell(3, 3), diagonal(3, 3), error, inverse(3, 3), k(3), orient(3, 3), &
         primitive(3, 3), scale, skew(3, 3), tau(3, 4)

      imaginary = CMPLX(0, 1, dp)
      primitive = RESHAPE([1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.5_dp, 0.5_dp, 0.5_dp], [3, 3])
      inverse = inv_3x3(primitive)
      angle = 0.73_dp
      skew = RESHAPE([0, 1, -1, -1, 0, 1, 1, -1, 0], [3, 3])/SQRT(3.0_dp)
      orient = (1 - COS(angle))/3 + SIN(angle)*skew
      DO i = 1, 3
         orient(i, i) = orient(i, i) + COS(angle)
      END DO
      cell = MATMUL(orient, primitive)
      sigma = 0
      sigma(1, 2, 1) = 1
      sigma(2, 1, 1) = 1
      sigma(1, 2, 2) = -imaginary
      sigma(2, 1, 2) = imaginary
      sigma(1, 1, 3) = 1
      sigma(2, 2, 3) = -1
      frame = -imaginary*SIN(angle/2)/SQRT(3.0_dp)*SUM(sigma, DIM=3)
      DO i = 1, 2
         frame(i, i) = frame(i, i) + COS(angle/2)
      END DO
      ref_spin = 0
      ref_spin(1, 1, 1) = 1
      ref_spin(2, 2, 1) = 1
      DO g = 1, 4
         signs = -1
         IF (g == 1) THEN
            signs = 1
         ELSE
            signs(5 - g) = 1
            ref_spin(:, :, g) = -imaginary*MATMUL(MATMUL(frame, sigma(:, :, 5 - g)), CONJG(TRANSPOSE(frame)))
         END IF
         diagonal = 0
         DO i = 1, 3
            diagonal(i, i) = signs(i)
         END DO
         rotations(:, :, g) = NINT(MATMUL(MATMUL(inverse, diagonal), primitive))
      END DO
      tau = 0
      tau(:, 2) = 0.5_dp*inverse(:, 2)
      tau(:, 3) = 0.5_dp*inverse(:, 1)
      tau(:, 4) = 0.5_dp*inverse(:, 3)
      k = -0.5_dp*SUM(primitive, DIM=1)
      DO j = 1, 4
         ref(1, j) = 1
         ref(2, j) = imaginary*(2*MOD(j - 1, 2) - 1)
         ref(3, j) = imaginary*(2*((j - 1)/2) - 1)
         ref(4, j) = -imaginary*ref(2, j)*ref(3, j)
      END DO
      DO h = 1, 4
         DO g = 1, 4
            q = IEOR(g - 1, h - 1) + 1
            ref_product(g, h) = q
            ref_factor(g, h) = ref(g, 1)*ref(h, 1)/ref(q, 1)
         END DO
      END DO
      DO pass = 1, 3
         scale = 10.0_dp**(100*(pass - 2))
         order = [1, 2, 3, 4]
         IF (pass /= 2) order = [3, 1, 4, 2]
         CALL little_group_factor_system(scale*cell, rotations(:, :, order), tau(:, order), &
                                         SPREAD(.FALSE., 1, 4), k, .TRUE., operations, product, factor, &
                                         1.e-10_dp, status, spin_lifts=lifts)
         IF (status /= 0) ERROR STOP 'Rotated three-screw boundary algebra failed'
         DO g = 1, 4
            i = order(operations(g))
            phases(g) = SIGN(1.0_dp, REAL(SUM(CONJG(ref_spin(:, :, i))*lifts(:, :, g)), dp))
            IF (MAXVAL(ABS(lifts(:, :, g) - phases(g)*ref_spin(:, :, i))) > 1.e-10_dp) THEN
               ERROR STOP 'Three-screw spin lift is not the analytic Cartesian lift'
            END IF
         END DO
         CALL projective_characters(product, factor, table, 1.e-10_dp, error, status)
         IF (status /= 0) ERROR STOP 'Three-screw projective characters failed'
         CALL character_table_match(product, factor, SPREAD(.FALSE., 1, 4), table, ref_product, ref_factor, &
                                    SPREAD(.FALSE., 1, 4), ref, order(operations), phases, mapping, 1.e-10_dp, error, status)
         IF (status /= 0) ERROR STOP 'Analytic boundary labels changed with the cell or spin-lift convention'
      END DO
      PRINT *, 'Rotated/rescaled three-screw boundary and actual spin-lift provenance passed.'
   END SUBROUTINE screw_lift_tests

END PROGRAM topology_character_unittest
