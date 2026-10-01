!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Little-group products with nonsymmorphic, spin and antiunitary factor systems.
! **************************************************************************************************
MODULE topology_little_group
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE mathlib,                         ONLY: det_3x3,&
                                              inv_3x3
   USE topology_symmetry,               ONLY: spin_rotation

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: little_group_factor_system, validate_factor_system
   PUBLIC :: projective_characters
CONTAINS

! **************************************************************************************************
!> \brief Construct the little group modulo integer translations of the supplied cell.
!> \param cell Cartesian lattice vectors as columns, not necessarily primitive or orthogonal
!> \param rotations Integer direct-space operations x -> W*x + tau in fractional coordinates
!> \param translations Fractional tau; the representatives are retained, not reduced modulo one
!> \param antiunitary True for {W|tau}*Theta, with physical Theta^2=-1 for spinors
!> \param k Fractional reciprocal vector, with translation character exp(-2*pi*i*k.L)
!> \param spinful Include the physical spin-1/2 double-group factor
!> \param mapping Indices of retained input operations, identity first; unallocated on failure
!> \param product Product indices in mapping order; unallocated on failure
!> \param factor Sewing factors D_g*conjugate_if_anti_g(D_h)=factor(g,h)*D_gh
!> \param tolerance Absolute fractional-coordinate and algebra tolerance
!> \param status 0 success; -1 invalid geometry/data; -2 nonunique/nonclosed group or bad algebra
!> \param spin_lifts Optional actual spin lifts in mapping order, including Theta for antiunitary elements
!> \note This validates the selected little group, not the completeness of the supplied space group.
!>       Centering translations remain separate cosets. Characters for a conventional cell may
!>       therefore include folded translation sectors; this is not a primitive-cell EBR catalogue.
! **************************************************************************************************
   SUBROUTINE little_group_factor_system(cell, rotations, translations, antiunitary, k, spinful, &
                                         mapping, product, factor, tolerance, status, spin_lifts)
      REAL(dp), INTENT(IN)                               :: cell(3, 3)
      INTEGER, INTENT(IN)                                :: rotations(:, :, :)
      REAL(dp), INTENT(IN)                               :: translations(:, :)
      LOGICAL, INTENT(IN)                                :: antiunitary(:)
      REAL(dp), INTENT(IN)                               :: k(3)
      LOGICAL, INTENT(IN)                                :: spinful
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: mapping(:), PRODUCT(:, :)
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: factor(:, :)
      REAL(dp), INTENT(IN)                               :: tolerance
      INTEGER, INTENT(OUT)                               :: status
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT), OPTIONAL    :: spin_lifts(:, :, :)

      COMPLEX(dp)                                        :: composed(2, 2), phase, theta(2, 2)
      COMPLEX(dp), ALLOCATABLE                           :: factors(:, :), spin(:, :, :)
      INTEGER                                            :: g, gh, h, i, identity, j, ng, ns, s, t
      INTEGER(int_8)                                     :: rw(3, 3)
      INTEGER, ALLOCATABLE                               :: members(:), products(:, :), selected(:)
      REAL(dp)                                           :: angle, cart(3, 3), delta(3), error, &
                                                            ident(3, 3), inverse(3, 3), kp(3), &
                                                            lattice(3, 3), shift(3), spin_sign

      status = -1
      ns = SIZE(rotations, 3)
      IF (ns < 1 .OR. SIZE(rotations, 1) /= 3 .OR. SIZE(rotations, 2) /= 3) RETURN
      IF (ANY(SHAPE(translations) /= [3, ns]) .OR. SIZE(antiunitary) /= ns) RETURN
      IF (.NOT. ieee_is_finite(tolerance) .OR. tolerance <= 0.0_dp .OR. tolerance >= 0.01_dp) RETURN
      IF (.NOT. ALL(ieee_is_finite(cell)) .OR. .NOT. ALL(ieee_is_finite(k))) RETURN
      IF (.NOT. ALL(ieee_is_finite(translations))) RETURN
      IF (MAXVAL(ABS(cell)) < TINY(1.0_dp)) RETURN
      ! Keep the integer rotation products exact and fractional rounding numerically meaningful.
      IF (MAXVAL(ABS(INT(rotations, int_8))) > INT(SQRT(REAL(HUGE(0_int_8)/3, dp)), int_8)) RETURN
      IF (MAXVAL(ABS(translations)) > tolerance/(32*EPSILON(1.0_dp)) &
          /MAX(1.0_dp, MAXVAL(ABS(REAL(rotations, dp))))) RETURN
      IF (MAXVAL(ABS(k)) > tolerance/(32*EPSILON(1.0_dp))) RETURN
      lattice = cell/MAXVAL(ABS(cell))
      IF (ABS(det_3x3(lattice)) < EPSILON(1.0_dp)) RETURN
      inverse = inv_3x3(lattice)
      ident = 0.0_dp
      theta = 0.0_dp
      DO i = 1, 3
         ident(i, i) = 1.0_dp
      END DO
      theta(1, 2) = 1.0_dp
      theta(2, 1) = -1.0_dp
      kp = MODULO(k, 1.0_dp)
      ALLOCATE (selected(ns), spin(2, 2, ns))
      ng = 0
      identity = 0
      DO g = 1, ns
         cart = MATMUL(MATMUL(lattice, REAL(rotations(:, :, g), dp)), inverse)
         CALL spin_rotation(cart, spin(:, :, g), tolerance, s)
         IF (s /= 0) RETURN
         IF (antiunitary(g)) spin(:, :, g) = MATMUL(spin(:, :, g), theta)
         delta = MATMUL(REAL(TRANSPOSE(rotations(:, :, g)), dp), kp)
         IF (antiunitary(g)) THEN
            delta = delta + kp
         ELSE
            delta = delta - kp
         END IF
         IF (MAXVAL(ABS(delta - ANINT(delta))) > tolerance) CYCLE
         ng = ng + 1
         selected(ng) = g
         IF (.NOT. antiunitary(g) .AND. ALL(rotations(:, :, g) == NINT(ident)) &
             .AND. MAXVAL(ABS(translations(:, g))) <= tolerance) identity = ng
      END DO
      status = -2
      IF (identity == 0) RETURN
      ! A pure integer translation is not silently substituted for the identity representative.
      s = selected(1)
      selected(1) = selected(identity)
      selected(identity) = s
      members = selected(:ng)
      ALLOCATE (products(ng, ng), factors(ng, ng))
      DO h = 1, ng
         t = members(h)
         DO g = 1, ng
            s = members(g)
            rw = MATMUL(INT(rotations(:, :, s), int_8), INT(rotations(:, :, t), int_8))
            shift = translations(:, s) + MATMUL(REAL(rotations(:, :, s), dp), translations(:, t))
            gh = 0
            DO j = 1, ng
               i = members(j)
               IF (antiunitary(i) .NEQV. (antiunitary(s) .NEQV. antiunitary(t))) CYCLE
               IF (ANY(INT(rotations(:, :, i), int_8) /= rw)) CYCLE
               delta = shift - translations(:, i)
               IF (MAXVAL(ABS(delta - ANINT(delta))) > tolerance) CYCLE
               IF (gh /= 0) RETURN
               gh = j
            END DO
            IF (gh == 0) RETURN
            products(g, h) = gh
            delta = ANINT(shift - translations(:, members(gh)))
            angle = -2*ACOS(-1.0_dp)*MODULO(DOT_PRODUCT(kp, delta), 1.0_dp)
            phase = CMPLX(COS(angle), SIN(angle), dp)
            IF (spinful) THEN
               IF (antiunitary(s)) THEN
                  composed = MATMUL(spin(:, :, s), CONJG(spin(:, :, t)))
               ELSE
                  composed = MATMUL(spin(:, :, s), spin(:, :, t))
               END IF
               spin_sign = SIGN(1.0_dp, REAL(SUM(CONJG(spin(:, :, members(gh)))*composed), dp))
               IF (MAXVAL(ABS(composed - spin_sign*spin(:, :, members(gh)))) > tolerance) RETURN
               phase = spin_sign*phase
            END IF
            factors(g, h) = phase
         END DO
      END DO
      CALL validate_factor_system(products, factors, antiunitary(members), tolerance, error, status)
      IF (status /= 0) RETURN
      IF (PRESENT(spin_lifts)) THEN
         ALLOCATE (spin_lifts(2, 2, ng))
         spin_lifts(:, :, :) = spin(:, :, members)
      END IF
      CALL MOVE_ALLOC(members, mapping)
      CALL MOVE_ALLOC(products, product)
      CALL MOVE_ALLOC(factors, factor)
   END SUBROUTINE little_group_factor_system

! **************************************************************************************************
!> \brief Check a finite group and its normalized, possibly antiunitary, U(1) factor system.
!> \param product Product table, with identity at index one
!> \param factor Complex phases multiplying the product-operation representation
!> \param antiunitary Whether an operation conjugates scalars to its right
!> \param tolerance Absolute algebra tolerance
!> \param residual Largest normalization or twisted cocycle residual
!> \param status 0 success; -1 malformed data; -2 invalid group or factor algebra
! **************************************************************************************************
   SUBROUTINE validate_factor_system(product, factor, antiunitary, tolerance, residual, status)
      INTEGER, INTENT(IN)                                :: product(:, :)
      COMPLEX(dp), INTENT(IN)                            :: factor(:, :)
      LOGICAL, INTENT(IN)                                :: antiunitary(:)
      REAL(dp), INTENT(IN)                               :: tolerance
      REAL(dp), INTENT(OUT)                              :: residual
      INTEGER, INTENT(OUT)                               :: status

      COMPLEX(dp)                                        :: right
      INTEGER                                            :: g, gh, h, hj, j, n

      status = -1
      residual = HUGE(1.0_dp)
      n = SIZE(antiunitary)
      IF (n < 1 .OR. ANY(SHAPE(product) /= [n, n]) .OR. ANY(SHAPE(factor) /= [n, n])) RETURN
      IF (ANY(product < 1) .OR. ANY(product > n)) RETURN
      IF (.NOT. ieee_is_finite(tolerance) .OR. tolerance <= 0.0_dp .OR. tolerance >= 0.01_dp) RETURN
      IF (.NOT. ALL(ieee_is_finite(REAL(factor, dp))) .OR. .NOT. ALL(ieee_is_finite(AIMAG(factor)))) RETURN
      status = -2
      IF (antiunitary(1)) RETURN
      DO g = 1, n
         IF (PRODUCT(g, 1) /= g .OR. PRODUCT(1, g) /= g) RETURN
         IF (COUNT(PRODUCT(g, :) == 1) /= 1 .OR. COUNT(PRODUCT(:, g) == 1) /= 1) RETURN
         DO h = 1, n
            IF (antiunitary(PRODUCT(g, h)) .NEQV. (antiunitary(g) .NEQV. antiunitary(h))) RETURN
         END DO
      END DO
      residual = MAX(MAXVAL(ABS(ABS(factor) - 1.0_dp)), &
                     MAXVAL(ABS(factor(:, 1) - 1.0_dp)), MAXVAL(ABS(factor(1, :) - 1.0_dp)))
      IF (residual > tolerance) RETURN
      DO j = 1, n
         DO h = 1, n
            hj = PRODUCT(h, j)
            DO g = 1, n
               gh = PRODUCT(g, h)
               IF (PRODUCT(g, hj) /= PRODUCT(gh, j)) RETURN
               right = factor(h, j)
               IF (antiunitary(g)) right = CONJG(right)
               residual = MAX(residual, ABS(factor(g, h)*factor(gh, j) - right*factor(g, hj)))
            END DO
         END DO
      END DO
      IF (residual > tolerance) RETURN
      status = 0
   END SUBROUTINE validate_factor_system

! **************************************************************************************************
!> \brief Generate all unitary projective characters from a validated finite factor system.
!> \param product Group multiplication table, identity first
!> \param factor Normalized U(1) factor system; no antiunitary elements are permitted
!> \param table Characters per group element and irrep; unallocated on failure, order unspecified
!> \param tolerance Algebra, subspace and relative eigenvalue-clustering tolerance
!> \param residual Largest algebra, subspace or character-orthogonality residual
!> \param status 0 success; -1 invalid input; -2 invalid/incompletely resolved algebra; -3 LAPACK error
!> \note Right multiplication in the twisted group algebra commutes with its left regular action.
!>       Generic Hermitian right actions split it into invariant irreducible copies. Retry accidental
!>       coincidences with deterministic coefficients; never accept a partial character table.
!>       These are not named space-group irreps, magnetic corepresentations or EBR classifications.
! **************************************************************************************************
   SUBROUTINE projective_characters(product, factor, table, tolerance, residual, status)
      INTEGER, INTENT(IN)                                :: product(:, :)
      COMPLEX(dp), INTENT(IN)                            :: factor(:, :)
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: table(:, :)
      REAL(dp), INTENT(IN)                               :: tolerance
      REAL(dp), INTENT(OUT)                              :: residual
      INTEGER, INTENT(OUT)                               :: status

      COMPLEX(dp)                                        :: coefficient, overlap
      COMPLEX(dp), ALLOCATABLE                           :: block(:, :), CHARACTER(:), &
                                                            characters(:, :), image(:, :), &
                                                            regular(:), vectors(:, :), work(:)
      INTEGER                                            :: d, first, found, g, h, i, info, j, last, &
                                                            n, ni, trial
      INTEGER, ALLOCATABLE                               :: copies(:), dimensions(:)
      LOGICAL                                            :: valid
      REAL(dp)                                           :: algebra_error, scale
      REAL(dp), ALLOCATABLE                              :: eigenvalues(:), rwork(:)

      status = -1
      residual = HUGE(1.0_dp)
      n = SIZE(product, 1)
      IF (n < 1 .OR. n > 1024) RETURN
      CALL validate_factor_system(product, factor, SPREAD(.FALSE., 1, n), tolerance, algebra_error, status)
      IF (status /= 0) RETURN
      ALLOCATE (vectors(n, n), characters(n, n), CHARACTER(n), work(MAX(1, 2*n - 1)), &
                eigenvalues(n), rwork(MAX(1, 3*n - 2)), copies(n), dimensions(n), regular(n))
      DO trial = 1, 8
         vectors = 0.0_dp
         DO h = 1, n
            coefficient = CMPLX(SIN(REAL(h*trial, dp)*SQRT(2.0_dp)), &
                                COS(REAL(h*h + trial*h, dp)*SQRT(3.0_dp)), dp)/n
            DO g = 1, n
               vectors(PRODUCT(g, h), g) = vectors(PRODUCT(g, h), g) + coefficient*factor(g, h)
            END DO
         END DO
         vectors = 0.5_dp*(vectors + CONJG(TRANSPOSE(vectors)))
         CALL zheev('V', 'U', n, vectors, n, eigenvalues, work, SIZE(work), rwork, info)
         status = -3
         IF (info /= 0) RETURN
         scale = MAX(1.0_dp, MAXVAL(ABS(eigenvalues)))
         residual = algebra_error
         ni = 0
         first = 1
         valid = .TRUE.
         copies = 0
         dimensions = 0
         DO WHILE (first <= n)
            last = first
            DO WHILE (last < n)
               IF (eigenvalues(last + 1) - eigenvalues(first) > tolerance*scale) EXIT
               last = last + 1
            END DO
            d = last - first + 1
            ALLOCATE (image(n, d), block(d, d))
            DO g = 1, n
               DO h = 1, n
                  image(PRODUCT(g, h), :) = factor(g, h)*vectors(h, first:last)
               END DO
               block = MATMUL(CONJG(TRANSPOSE(vectors(:, first:last))), image)
               residual = MAX(residual, MAXVAL(ABS(image - MATMUL(vectors(:, first:last), block))))
               CHARACTER(g) = 0.0_dp
               DO i = 1, d
                  CHARACTER(g) = CHARACTER(g) + block(i, i)
               END DO
            END DO
            DEALLOCATE (image, block)
            residual = MAX(residual, ABS(SUM(ABS(CHARACTER)**2)/n - 1.0_dp))
            IF (residual > tolerance) THEN
               valid = .FALSE.
               EXIT
            END IF
            found = 0
            DO j = 1, ni
               overlap = SUM(CONJG(characters(:, j))*CHARACTER)/n
               IF (ABS(overlap - 1.0_dp) <= tolerance) THEN
                  IF (dimensions(j) /= d) valid = .FALSE.
                  found = j
                  residual = MAX(residual, MAXVAL(ABS(CHARACTER - characters(:, j)))/n)
               ELSE
                  residual = MAX(residual, ABS(overlap))
               END IF
            END DO
            IF (found == 0) THEN
               ni = ni + 1
               found = ni
               characters(:, ni) = CHARACTER
               dimensions(ni) = d
            END IF
            copies(found) = copies(found) + 1
            first = last + 1
         END DO
         status = -2
         IF (.NOT. valid .OR. residual > tolerance) CYCLE
         IF (SUM(dimensions(:ni)**2) /= n .OR. ANY(copies(:ni) /= dimensions(:ni))) CYCLE
         regular = MATMUL(characters(:, :ni), CMPLX(dimensions(:ni), 0, dp))
         regular(1) = regular(1) - n
         residual = MAX(residual, MAXVAL(ABS(regular))/n)
         IF (.NOT. ieee_is_finite(residual) .OR. residual > tolerance) CYCLE
         table = characters(:, :ni)
         status = 0
         RETURN
      END DO
   END SUBROUTINE projective_characters
END MODULE topology_little_group
