!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
PROGRAM kpsym_unittest
   USE cp_log_handling,                 ONLY: cp_logger_get_default_io_unit
   USE cryssym,                         ONLY: crys_sym_cartesian_halfturn,&
                                              crys_sym_gen,&
                                              csym_type,&
                                              kpoint_gen,&
                                              release_csym_type
   USE kinds,                           ONLY: dp
   USE kpoint_methods,                  ONLY: kpoint_star_create,&
                                              kpoint_star_density,&
                                              kpoint_star_operator,&
                                              kpoint_star_operator_defect,&
                                              kpoint_star_orbitals
   USE kpoint_mo_symmetry_methods,      ONLY: kpoint_cell_gauge_action
   USE kpoint_types,                    ONLY: kpoint_star_type,&
                                              kpoint_sym_create,&
                                              kpoint_sym_type,&
                                              kpoint_type
   USE libcp2k,                         ONLY: cp2k_finalize,&
                                              cp2k_init
   USE physcon,                         ONLY: bohr
   USE spglib_f08,                      ONLY: spg_get_major_version

   IMPLICIT NONE
   REAL(dp), PARAMETER :: ice_ix(3, 36) = RESHAPE([ &
                                                  -0.006830_dp, 0.340009_dp, 0.190161_dp, &
                                                  0.113981_dp, 0.164412_dp, 0.283330_dp, &
                                                  0.309270_dp, 0.370346_dp, 0.106284_dp, &
                                                  0.159991_dp, 0.493170_dp, 0.440161_dp, &
                                                  0.335588_dp, 0.613981_dp, 0.533330_dp, &
                                                  0.129654_dp, 0.809270_dp, 0.356284_dp, &
                                                  0.840009_dp, 0.506830_dp, -0.059839_dp, &
                                                  0.664412_dp, 0.386019_dp, 0.033330_dp, &
                                                  0.870346_dp, 0.190730_dp, 0.856284_dp, &
                                                  0.506830_dp, 0.840009_dp, 1.059839_dp, &
                                                  0.386019_dp, 0.664412_dp, 0.966670_dp, &
                                                  0.190730_dp, 0.870346_dp, 0.143716_dp, &
                                                  0.493170_dp, 0.159991_dp, 0.559839_dp, &
                                                  0.613981_dp, 0.335588_dp, 0.466670_dp, &
                                                  0.809270_dp, 0.129654_dp, 0.643716_dp, &
                                                  0.370346_dp, 0.309270_dp, -0.106284_dp, &
                                                  0.659991_dp, 1.006830_dp, 0.309839_dp, &
                                                  0.835588_dp, 0.886019_dp, 0.216670_dp, &
                                                  0.629654_dp, 0.690730_dp, 0.393716_dp, &
                                                  1.006830_dp, 0.659991_dp, 0.690161_dp, &
                                                  0.886019_dp, 0.835588_dp, 0.783330_dp, &
                                                  0.690730_dp, 0.629654_dp, 0.606284_dp, &
                                                  0.340009_dp, -0.006830_dp, 0.809839_dp, &
                                                  0.164412_dp, 0.113981_dp, 0.716670_dp, &
                                                  0.114088_dp, 0.310030_dp, 0.267366_dp, &
                                                  0.189970_dp, 0.614088_dp, 0.517366_dp, &
                                                  0.810030_dp, 0.385912_dp, 0.017365_dp, &
                                                  0.385912_dp, 0.810030_dp, 0.982634_dp, &
                                                  0.614088_dp, 0.189970_dp, 0.482634_dp, &
                                                  0.689970_dp, 0.885912_dp, 0.232634_dp, &
                                                  0.885912_dp, 0.689970_dp, 0.767366_dp, &
                                                  0.310030_dp, 0.114088_dp, 0.732634_dp, &
                                                  0.404256_dp, 0.404256_dp, 0.000000_dp, &
                                                  0.095744_dp, 0.904256_dp, 0.250000_dp, &
                                                  0.904256_dp, 0.095744_dp, 0.750000_dp, &
                                                  0.595744_dp, 0.595744_dp, 0.500000_dp], [3, 36])
   REAL(dp) :: cell(3, 3), coords(3, 36), halfturn_cell(3, 3), shifted(3, 36)
   INTEGER :: i, io_unit
   INTEGER :: types(36)

   CALL cp2k_init()
   CALL check_star_adjoint()
   io_unit = cp_logger_get_default_io_unit()
   cell = 0.0_dp
   cell(1, 1) = 6.784903_dp*bohr
   cell(2, 2) = 6.784903_dp*bohr
   cell(3, 3) = 6.807916_dp*bohr
   types(:24) = 1
   types(25:) = 2
   coords = ice_ix

   ! Rounded coordinates put seven operations within a few ulps of the boundary.
   CALL check_mesh(coords, cell, 1.0e-6_dp, 6, .FALSE., .TRUE.)
   CALL check_mesh(coords, cell, 1.0e-7_dp, 14, .FALSE., .TRUE.)
   ! FULL_GRID overrides both atomic symmetry and time-reversal reduction.
   CALL check_mesh(coords, cell, 1.0e-6_dp, 27, .FALSE., .FALSE., full_grid=.TRUE.)
   shifted = coords + SPREAD([0.137_dp, -0.231_dp, 0.319_dp], 2, 36)
   CALL check_mesh(shifted, cell, 1.0e-6_dp, 6, .FALSE., .TRUE.)
   CALL check_mesh(coords, 1.3_dp*cell, 1.0e-6_dp, 6, .FALSE., .TRUE.)
   DO i = 1, 24
      shifted(:, i) = coords(:, 25 - i)
   END DO
   DO i = 25, 36
      shifted(:, i) = coords(:, 61 - i)
   END DO
   CALL check_mesh(shifted, cell, 1.0e-6_dp, 6, .FALSE., .TRUE.)

   ! Restrict the complete orthorhombic group without adding little-group weights.
   halfturn_cell = 0.0_dp
   halfturn_cell(1, 1) = 7.0_dp
   halfturn_cell(2, 2) = 6.4_dp
   halfturn_cell(3, 3) = 6.7_dp
   CALL check_mesh(RESHAPE([0.0_dp, 0.0_dp, 0.0_dp], [3, 1]), halfturn_cell, &
                   1.e-8_dp, 8, .FALSE., .TRUE., cartesian_half_turns=.TRUE.)
   IF (spg_get_major_version() > 0) THEN
      CALL check_mesh(RESHAPE([0.0_dp, 0.0_dp, 0.0_dp], [3, 1]), halfturn_cell, &
                      1.e-8_dp, 8, .FALSE., .TRUE., cartesian_half_turns=.TRUE., &
                      use_spglib_reduction=.TRUE.)
      CALL check_mesh(RESHAPE([0.0_dp, 0.0_dp, 0.0_dp], [3, 1]), halfturn_cell, &
                      1.e-8_dp, 8, .FALSE., .TRUE., cartesian_half_turns=.TRUE., &
                      use_spglib_reduction=.TRUE., use_spglib_backend=.TRUE.)
   END IF

   ! Genuine tolerance-induced nonclosure must still trigger the orbit guard.
   coords = ice_ix
   coords(3, 27) = coords(3, 27) + 1.0e-6_dp
   coords(3, 33) = coords(3, 33) + 0.75e-6_dp
   coords(3, 36) = coords(3, 36) - 0.75e-6_dp
   CALL check_mesh(coords, cell, 1.0e-6_dp, 14, .TRUE., .FALSE.)
   IF (io_unit > 0) WRITE (io_unit, *) "K290 tolerance-boundary and orbit-guard tests passed"
   CALL cp2k_finalize()

CONTAINS

! **************************************************************************************************
!> \brief Check native symmetry phases, nonorthogonal AO rotations and the weighted trace adjoint.
!>        Include atom permutation, time reversal, unequal and zero representative weights,
!>        unused little-group metadata, and independence from legacy SCF phase calibration.
! **************************************************************************************************
   SUBROUTINE check_star_adjoint()
      COMPLEX(dp) :: action(6, 6), base(6, 6), coefficients(6, 2, 2), covariant(6, 6, 4), &
         density(6, 6, 2), expected(6, 2), inverse(6, 6), operators(6, 6, 4), reference(6, 6, 2)
      COMPLEX(dp), ALLOCATABLE                           :: expanded(:, :, :), orbitals(:, :, :), &
                                                            reduced(:, :, :), selected(:, :, :)
      INTEGER                                            :: i, ik, j, link, pass, status
      LOGICAL                                            :: antiunitary
      REAL(dp)                                           :: defect, lhs, metric_scale(6), rhs, &
                                                            rotation(3, 3), scale
      TYPE(kpoint_star_type)                             :: map
      TYPE(kpoint_sym_type), POINTER                     :: sym
      TYPE(kpoint_type), POINTER                         :: kpoint

      ALLOCATE (kpoint)
      kpoint%nkp = 2
      ALLOCATE (kpoint%xkp(3, 2), kpoint%wkp(2), kpoint%kp_sym(2), kpoint%atype(2), &
                kpoint%ibrot(2), kpoint%kind_rotmat(2, 1))
      kpoint%xkp(:, 1) = [0.13_dp, 0.21_dp, 0.07_dp]
      kpoint%xkp(:, 2) = [0.31_dp, 0.17_dp, 0.11_dp]
      kpoint%atype = 1
      kpoint%ibrot = [1, 2]
      DO i = 1, 2
         ALLOCATE (kpoint%kind_rotmat(i, 1)%rmat(3, 3))
         CALL kpoint_sym_create(kpoint%kp_sym(i)%kpoint_sym)
      END DO
      kpoint%kind_rotmat(1, 1)%rmat = 0.0_dp
      DO i = 1, 3
         kpoint%kind_rotmat(1, 1)%rmat(i, i) = 1.0_dp
      END DO
      ! A Cartesian rotation has a nonorthogonal representation in a scaled AO basis.
      kpoint%kind_rotmat(2, 1)%rmat = RESHAPE([0.0_dp, 0.5_dp, 0.0_dp, &
                                               -2.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [3, 3])
      sym => kpoint%kp_sym(1)%kpoint_sym
      sym%apply_symmetry = .TRUE.
      sym%nwght = 3
      sym%nwred = 4
      ALLOCATE (sym%rot(3, 3, 4), sym%rotp(4), sym%f0(2, 4), sym%fcell_gauge(3, 2, 4), &
                sym%phase_mode(4), sym%xkp(3, 4))
      DO i = 1, 4
         sym%rot(:, :, i) = kpoint%kind_rotmat(1, 1)%rmat
         sym%f0(:, i) = [1, 2]
         sym%xkp(:, i) = kpoint%xkp(:, 1)
      END DO
      sym%rot(:, :, 2) = RESHAPE([0.0_dp, 1.0_dp, 0.0_dp, &
                                  -1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [3, 3])
      sym%rot(:, :, 3) = sym%rot(:, :, 2)
      sym%rotp = [1, 2, -2, 1]
      sym%f0(:, 2) = [2, 1]
      sym%f0(:, 3) = [2, 1]
      sym%fcell_gauge = 0
      sym%fcell_gauge(:, 1, 2) = [1, 0, 0]
      sym%fcell_gauge(:, 2, 2) = [0, -1, 1]
      sym%fcell_gauge(:, :, 3) = sym%fcell_gauge(:, :, 2)
      sym%phase_mode = [1, 1, 2, 1]
      sym%xkp(:, 2) = [-0.21_dp, 0.13_dp, 0.07_dp]
      sym%xkp(:, 3) = -sym%xkp(:, 2)
      DO ik = 1, 2
         DO j = 1, 2
            DO i = 1, 6
               coefficients(i, j, ik) = CMPLX(SIN(0.17_dp*REAL(i + 2*j + 3*ik, dp)), &
                                              COS(0.23_dp*REAL(2*i - j + 4*ik, dp)), dp)
            END DO
         END DO
         density(:, :, ik) = MATMUL(coefficients(:, :, ik), CONJG(TRANSPOSE(coefficients(:, :, ik))))
      END DO
      DO link = 1, 4
         DO j = 1, 6
            DO i = 1, 6
               operators(i, j, link) = CMPLX(SIN(0.19_dp*REAL(i + 2*j + link, dp)), &
                                             COS(0.29_dp*REAL(2*i - j + 3*link, dp)), dp)
            END DO
         END DO
         operators(:, :, link) = 0.5_dp*(operators(:, :, link) + CONJG(TRANSPOSE(operators(:, :, link))))
      END DO
      DO pass = 1, 2
         kpoint%wkp = [0.37_dp, 0.63_dp]
         sym%phase_mode = 0
         IF (pass == 2) kpoint%wkp = [0.0_dp, 1.0_dp]
         IF (pass == 2) sym%phase_mode = [1, 2, 1, 2]
         CALL kpoint_star_create(kpoint, [3, 3], map)
         IF (SIZE(map%weights) /= 4 .OR. ANY(map%representative /= [1, 1, 1, 2])) THEN
            ERROR STOP 'Little-group metadata added spurious density contributions'
         END IF
         CALL kpoint_star_density(map, density, expanded)
         CALL kpoint_star_orbitals(map, coefficients, orbitals)
         DO link = 1, 4
            ik = map%representative(link)
            CALL kpoint_cell_gauge_action(kpoint, ik, map%operation(link), map%vectors(:, link), &
                                          action, rotation, antiunitary, status)
            IF (status /= 0) ERROR STOP 'The independent cell-gauge AO action is unavailable'
            IF (antiunitary) THEN
               expected = MATMUL(action, CONJG(coefficients(:, :, ik)))
            ELSE
               expected = MATMUL(action, coefficients(:, :, ik))
            END IF
            IF (.NOT. MAXVAL(ABS(orbitals(:, :, link) - expected)) <= 1.e-13_dp) THEN
               ERROR STOP 'Orbital symmetry changed an atom permutation, AO rotation or Bloch phase'
            END IF
            IF (.NOT. MAXVAL(ABS(expanded(:, :, link) - MATMUL(expected, CONJG(TRANSPOSE(expected))))) <= 1.e-12_dp) THEN
               ERROR STOP 'Native density symmetry differs from transformed orbital factors'
            END IF
         END DO
         ! Independently construct covariant operators using the scaled AO metric.
         metric_scale = [.25_dp, 1._dp, 1._dp, .25_dp, 1._dp, 1._dp]
         DO link = 1, SIZE(map%weights)
            ik = map%representative(link)
            CALL kpoint_cell_gauge_action(kpoint, ik, map%operation(link), map%vectors(:, link), &
                                          action, rotation, antiunitary, status)
            IF (status /= 0) ERROR STOP 'Covariance control has no independent AO action'
            DO j = 1, 6
               DO i = 1, 6
                  inverse(i, j) = CONJG(action(j, i))*metric_scale(j)/metric_scale(i)
               END DO
            END DO
            base = density(:, :, ik)
            IF (antiunitary) base = CONJG(base)
            covariant(:, :, link) = MATMUL(CONJG(TRANSPOSE(inverse)), MATMUL(base, inverse))
         END DO
         CALL kpoint_star_operator_defect(map, covariant, defect, scale)
         IF (defect > 1.e-12_dp*MAX(1._dp, scale)) ERROR STOP 'Covariant physical operators failed the star check'
         covariant(1, 1, 2) = covariant(1, 1, 2) + .1_dp
         CALL kpoint_star_operator_defect(map, covariant, defect, scale)
         IF (defect < 1.e-4_dp) ERROR STOP 'Broken operator covariance was hidden by star averaging'
         CALL kpoint_star_operator(map, operators, reduced)
         DO ik = 1, 2
            CALL kpoint_star_operator(map, operators, selected, representative=ik)
            IF (SIZE(selected, 3) /= 1) ERROR STOP 'A selected operator has more than one representative'
            IF (.NOT. MAXVAL(ABS(selected(:, :, 1) - reduced(:, :, ik))) <= 1.e-13_dp) THEN
               ERROR STOP 'A selected operator differs from the full irreducible adjoint'
            END IF
         END DO
         lhs = 0.0_dp
         rhs = 0.0_dp
         DO link = 1, 4
            lhs = lhs + map%weights(link)*REAL(SUM(CONJG(operators(:, :, link))*expanded(:, :, link)), dp)
         END DO
         DO ik = 1, 2
            rhs = rhs + kpoint%wkp(ik)*REAL(SUM(CONJG(reduced(:, :, ik))*density(:, :, ik)), dp)
         END DO
         IF (.NOT. ABS(lhs - rhs) <= 1.e-12_dp*MAX(1.0_dp, ABS(lhs), ABS(rhs))) THEN
            ERROR STOP 'Weighted k-point symmetry is not paired with its operator adjoint'
         END IF
         IF (pass == 1) reference = reduced
         IF (pass == 2) THEN
            IF (.NOT. MAXVAL(ABS(reduced - reference)) <= 1.e-13_dp) THEN
               ERROR STOP 'Zero representative weight changed an unweighted symmetry operator'
            END IF
         END IF
      END DO
      ! A fixed-point operation is usable independently of Brillouin-zone reduction.
      sym%apply_symmetry = .FALSE.
      sym%nwght = 1
      kpoint%xkp(:, 1) = 0.0_dp
      kpoint%wkp = [0.37_dp, 0.63_dp]
      sym%xkp = 0.0_dp
      CALL kpoint_cell_gauge_action(kpoint, 1, 2, kpoint%xkp(:, 1), action, rotation, antiunitary, status)
      IF (status /= 0 .OR. antiunitary) ERROR STOP 'Fixed-point AO action is unavailable'
      IF (MAXVAL(ABS(action(1:3, 4:6) - kpoint%kind_rotmat(2, 1)%rmat)) > 1.e-13_dp .OR. &
          MAXVAL(ABS(action(4:6, 1:3) - kpoint%kind_rotmat(2, 1)%rmat)) > 1.e-13_dp) THEN
         ERROR STOP 'Fixed-point AO action changed the native rotation or atom map'
      END IF
      CALL kpoint_star_create(kpoint, [3, 3], map)
      IF (SIZE(map%weights) /= 2) ERROR STOP 'Fixed-point operations added mesh contributions'
      IF (ANY(map%weights /= kpoint%wkp)) ERROR STOP 'Fixed-point operations changed mesh weights'
      CALL kpoint_star_density(map, density, expanded)
      IF (MAXVAL(ABS(expanded - density)) > 1.e-12_dp) THEN
         ERROR STOP 'Storing fixed-point operations changed the selected star density'
      END IF
      ! A complete reflection/conjugation group at Gamma adds no mesh points.
      kpoint%kind_rotmat(2, 1)%rmat(1, 2) = 2.0_dp
      sym%rot(:, :, 2) = RESHAPE([0.0_dp, 1.0_dp, 0.0_dp, &
                                  1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [3, 3])
      sym%rot(:, :, 3) = sym%rot(:, :, 2)
      sym%rotp(4) = -1
      covariant = 0.0_dp
      DO link = 1, 2
         DO i = 1, 6
            covariant(i, i, link) = metric_scale(i)
         END DO
      END DO
      CALL kpoint_star_operator_defect(map, covariant(:, :, :2), defect, scale, little_group=.TRUE.)
      IF (.NOT. defect <= 1.e-12_dp) ERROR STOP 'A compatible Gamma operator failed its complete-group check'
      covariant(1, 1, 1) = covariant(1, 1, 1) + .1_dp
      CALL kpoint_star_operator_defect(map, covariant(:, :, :2), defect, scale)
      IF (.NOT. defect <= 1.e-12_dp) ERROR STOP 'A fixed-point operation changed the distinct-star comparison'
      CALL kpoint_star_operator_defect(map, covariant(:, :, :2), defect, scale, little_group=.TRUE.)
      IF (.NOT. defect >= 1.e-4_dp) ERROR STOP 'A broken fixed-point operator escaped the complete-group check'
      covariant(1, 1, 1) = metric_scale(1)
      covariant(1, 2, 1) = CMPLX(0.0_dp, .1_dp, dp)
      covariant(2, 1, 1) = CONJG(covariant(1, 2, 1))
      CALL kpoint_star_operator_defect(map, covariant(:, :, :2), defect, scale, little_group=.TRUE.)
      IF (.NOT. defect >= 1.e-4_dp) ERROR STOP 'A broken conjugating fixed-point operator was accepted'
      kpoint%full_grid = .TRUE.
      CALL kpoint_star_operator_defect(map, covariant(:, :, :2), defect, scale, little_group=.TRUE.)
      IF (.NOT. defect <= 1.e-12_dp) ERROR STOP 'The full physical mesh acquired a fixed-point restriction'
      map = kpoint_star_type()
      DEALLOCATE (sym%rot, sym%rotp, sym%f0, sym%fcell_gauge, sym%phase_mode, sym%xkp)
      DO i = 1, 2
         DEALLOCATE (kpoint%kp_sym(i)%kpoint_sym, kpoint%kind_rotmat(i, 1)%rmat)
      END DO
      DEALLOCATE (kpoint%xkp, kpoint%wkp, kpoint%kp_sym, kpoint%atype, kpoint%ibrot, kpoint%kind_rotmat)
      DEALLOCATE (kpoint)
   END SUBROUTINE check_star_adjoint

! **************************************************************************************************
!> \brief Check mesh reduction and atomic operation consistency without an SCF calculation.
!> \param x fractional atomic coordinates
!> \param h cell matrix in bohr
!> \param delta K290 symmetry tolerance
!> \param expected expected irreducible k-point count
!> \param fallback expected inversion fallback
!> \param closed whether to check closure of the retained atomic operations
!> \param full_grid keep every mesh point, including time-reversed partners
!> \param cartesian_half_turns ...
!> \param use_spglib_reduction ...
!> \param use_spglib_backend ...
! **************************************************************************************************
   SUBROUTINE check_mesh(x, h, delta, expected, fallback, closed, full_grid, cartesian_half_turns, &
                         use_spglib_reduction, use_spglib_backend)
      REAL(dp), INTENT(IN)                               :: x(:, :), h(3, 3), delta
      INTEGER, INTENT(IN)                                :: expected
      LOGICAL, INTENT(IN)                                :: fallback, closed
      LOGICAL, INTENT(IN), OPTIONAL                      :: full_grid, cartesian_half_turns, &
                                                            use_spglib_reduction, &
                                                            use_spglib_backend

      INTEGER                                            :: a, i, j, k
      LOGICAL                                            :: found, spglib
      REAL(dp)                                           :: product(3, 3)
      TYPE(csym_type)                                    :: csym

      spglib = .FALSE.
      IF (PRESENT(use_spglib_reduction)) spglib = use_spglib_reduction
      CALL crys_sym_gen(csym, x, types(:SIZE(x, 2)), h, delta=delta, use_spglib=spglib)
      IF (PRESENT(full_grid)) THEN
         CALL kpoint_gen(csym, [3, 3, 3], symm=.TRUE., full_grid=full_grid, &
                         cartesian_half_turns=cartesian_half_turns, &
                         use_spglib_reduction=use_spglib_reduction, use_spglib_backend=use_spglib_backend)
      ELSE
         CALL kpoint_gen(csym, [3, 3, 3], symm=.TRUE., shift=[0.0_dp, 0.0_dp, 0.0_dp], &
                         cartesian_half_turns=cartesian_half_turns, &
                         use_spglib_reduction=use_spglib_reduction, use_spglib_backend=use_spglib_backend)
      END IF
      IF (csym%nkpoint /= expected) ERROR STOP "Incorrect irreducible mesh size"
      IF (csym%inversion_only .NEQV. fallback) ERROR STOP "Incorrect orbit-guard decision"
      IF (ABS(SUM(csym%wkpoint) - 27.0_dp) > 1.0e-12_dp) ERROR STOP "Incorrect weight sum"
      IF (ANY(csym%kplink <= 0)) ERROR STOP "Incomplete mesh mapping"
      IF (PRESENT(full_grid)) THEN
         IF (full_grid) THEN
            DO i = 1, 27
               IF (ANY(csym%kplink(:, i) /= i)) ERROR STOP "Full mesh has a reduced point mapping"
            END DO
            IF (ANY(ABS(csym%wkpoint - 1.0_dp) > 1.0e-12_dp)) ERROR STOP "Full mesh has reduced weights"
            IF (ANY(ABS(csym%xkpoint - csym%kpmesh) > 1.0e-12_dp)) ERROR STOP "Full mesh reordered points"
         END IF
      END IF
      IF (PRESENT(cartesian_half_turns)) THEN
         IF (cartesian_half_turns) THEN
            IF (csym%nrtot /= 4) ERROR STOP "Incorrect Cartesian half-turn group size"
            DO k = 1, csym%nrtot
               IF (.NOT. crys_sym_cartesian_halfturn(csym%rt(:, :, k))) ERROR STOP "Non-half-turn operation retained"
            END DO
         END IF
      END IF
      DO k = 1, csym%nrtot
         DO a = 1, SIZE(x, 2)
            IF (COUNT(csym%f0(:, k) == a) /= 1) ERROR STOP "Nonbijective atom mapping"
         END DO
         IF (ANY(types(csym%f0(:, k)) /= types(:SIZE(x, 2)))) ERROR STOP "Atom types not preserved"
      END DO
      IF (closed) THEN
         DO i = 1, csym%nrtot
            DO j = 1, csym%nrtot
               product = MATMUL(csym%rt(:, :, i), csym%rt(:, :, j))
               found = .FALSE.
               DO k = 1, csym%nrtot
                  IF (MAXVAL(ABS(product - csym%rt(:, :, k))) > 1.0e-12_dp) CYCLE
                  IF (ANY(csym%f0(csym%f0(:, j), i) /= csym%f0(:, k))) CYCLE
                  found = .TRUE.
               END DO
               IF (.NOT. found) ERROR STOP "Operation group not closed"
            END DO
         END DO
      END IF
      CALL release_csym_type(csym)
   END SUBROUTINE check_mesh
END PROGRAM kpsym_unittest
