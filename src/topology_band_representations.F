!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Site-induced atomic band characters, with spin and antiunitary coset representatives.
!> \author Thomas D. Kuehne, tkuehne@cp2k.org
! **************************************************************************************************
MODULE topology_band_representations
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE mathconstants,                   ONLY: twopi
   USE mathlib,                         ONLY: inv_3x3
   USE topology_character_match,        ONLY: character_table_match
   USE topology_corepresentations,      ONLY: corepresentation_table_type,&
                                              projective_corepresentations
   USE topology_integer_lattice,        ONLY: integer_lattice_certificate,&
                                              integer_lattice_factor,&
                                              integer_lattice_product,&
                                              integer_lattice_type
   USE topology_little_group,           ONLY: little_group_factor_system,&
                                              projective_characters,&
                                              validate_factor_system
   USE topology_symmetry,               ONLY: spin_rotation

   IMPLICIT NONE
   PRIVATE

   !> Compact onsite provenance, without the orbit-sized induction cache.
   TYPE, PUBLIC :: atomic_site_data_type
      LOGICAL :: spinful = .FALSE.
      INTEGER, ALLOCATABLE :: mapping(:), PRODUCT(:, :)
      LOGICAL, ALLOCATABLE :: antiunitary(:)
      REAL(dp), ALLOCATABLE :: cartesian(:, :, :), lattice_shifts(:, :)
      COMPLEX(dp), ALLOCATABLE :: spin(:, :, :), factor(:, :), characters(:, :)
   END TYPE atomic_site_data_type

   !> Fourier characters on Fix(W^T) are indexed by Z^3/(W-I)Z^3.
   TYPE, PUBLIC :: atomic_character_block_type
      INTEGER :: operation = 0
      INTEGER(int_8) :: rotation(3, 3) = 0, left(3, 3) = 0, period(3) = 0
      INTEGER(int_8), ALLOCATABLE :: keys(:, :)
      COMPLEX(dp), ALLOCATABLE :: coefficients(:, :)
   END TYPE atomic_character_block_type

   TYPE, PUBLIC :: atomic_character_polynomial_type
      TYPE(atomic_character_block_type), ALLOCATABLE :: blocks(:)
      INTEGER, ALLOCATABLE :: dimensions(:)
   END TYPE atomic_character_polynomial_type

   !> Construct once per site; columns retain the same local representation at every k point.
   TYPE, PUBLIC :: atomic_band_type
      PRIVATE
      INTEGER, ALLOCATABLE :: mapping(:), rotations(:, :, :), fixed_unitary(:, :), dimensions(:), site_mapping(:)
      LOGICAL, ALLOCATABLE :: antiunitary(:), conjugate_site(:)
      REAL(dp), ALLOCATABLE :: positions(:, :), shifts(:, :, :)
      COMPLEX(dp), ALLOCATABLE :: site_characters(:, :), spin_factor(:, :)
      INTEGER, ALLOCATABLE :: site_product(:, :)
      LOGICAL :: spinful = .FALSE.
      LOGICAL, ALLOCATABLE :: site_antiunitary(:)
      REAL(dp), ALLOCATABLE :: site_cartesian(:, :, :)
      COMPLEX(dp), ALLOCATABLE :: site_factor(:, :), site_spin(:, :, :)
      INTEGER, ALLOCATABLE :: site_full_mapping(:)
      REAL(dp), ALLOCATABLE :: site_lattice_shifts(:, :)
   END TYPE atomic_band_type

   PUBLIC :: atomic_band_construct, atomic_band_characters, atomic_band_site_data
   PUBLIC :: atomic_band_site_match
   PUBLIC :: atomic_band_site_snapshot, atomic_site_storage
   PUBLIC :: atomic_band_character_polynomial, atomic_character_compare
CONTAINS

! **************************************************************************************************
!> \brief Collect finite Fourier coefficients on each complete unitary invariant reciprocal torus.
!> \param band Prepared induced representations
!> \param polynomial Sorted quotient-lattice keys and coefficients; empty on failure
!> \param status Zero success, negative for inconsistent data or integer-arithmetic bounds
!> \note For U(W-I)V=D, keys are U*L modulo each nonzero Smith divisor. Equal characters
!>       on Fix(W^T) have equal collected coefficients, including disconnected torsion components.
!>       Antiunitary orbit representatives still conjugate onsite characters. Antiunitary operators
!>       themselves have no basis-independent trace and are not assigned Fourier characters.
! **************************************************************************************************
   SUBROUTINE atomic_band_character_polynomial(band, polynomial, status)
      TYPE(atomic_band_type), INTENT(IN)                 :: band
      TYPE(atomic_character_polynomial_type), &
         INTENT(OUT)                                     :: polynomial
      INTEGER, INTENT(OUT)                               :: status

      COMPLEX(dp), ALLOCATABLE                           :: chi(:), coefficients(:, :)
      INTEGER                                            :: b, g, h, i, j, k, n, ni, no, nt, q, &
                                                            rank, row
      INTEGER(int_8)                                     :: key(3), matrix(3, 3)
      INTEGER(int_8), ALLOCATABLE                        :: diagonal(:), keys(:, :), left(:, :), &
                                                            right(:, :), shifts(:, :)
      INTEGER, ALLOCATABLE                               :: members(:), order(:)
      TYPE(atomic_character_polynomial_type)             :: work
      TYPE(integer_lattice_type)                         :: factor

      status = -1
      IF (.NOT. ALLOCATED(band%dimensions)) RETURN
      IF (MAXVAL(ABS(band%shifts)) > REAL(HUGE(0_int_8)/32, dp)) RETURN
      ni = SIZE(band%dimensions)
      no = SIZE(band%positions, 2)
      n = COUNT(.NOT. band%antiunitary)
      ALLOCATE (work%blocks(n), work%dimensions(ni), chi(ni))
      work%dimensions = no*band%dimensions
      b = 0
      DO q = 1, SIZE(band%mapping)
         g = FINDLOC(band%mapping, q, DIM=1)
         IF (g == 0) THEN
            status = -1
            RETURN
         END IF
         IF (band%antiunitary(g)) CYCLE
         b = b + 1
         work%blocks(b)%operation = q
         work%blocks(b)%rotation = INT(band%rotations(:, :, g), int_8)
         matrix = work%blocks(b)%rotation
         DO i = 1, 3
            matrix(i, i) = matrix(i, i) - 1_int_8
         END DO
         CALL integer_lattice_factor(matrix, factor, status)
         IF (status /= 0) RETURN
         CALL integer_lattice_certificate(factor, left, right, diagonal, rank)
         work%blocks(b)%left = left
         work%blocks(b)%period = diagonal
         members = PACK([(j, j=1, no)], band%fixed_unitary(g, :) /= 0)
         nt = SIZE(members)
         ALLOCATE (shifts(3, nt), order(nt), coefficients(nt, ni))
         shifts = NINT(band%shifts(:, g, members), int_8)
         CALL integer_lattice_product(left, shifts, keys, status)
         IF (status /= 0) RETURN
         DO i = 1, rank
            keys(i, :) = MODULO(keys(i, :), diagonal(i))
         END DO
         order = [(j, j=1, nt)]
         DO j = 2, nt
            k = order(j)
            i = j - 1
            DO WHILE (i >= 1)
               IF (key_compare(keys(:, order(i)), keys(:, k)) <= 0) EXIT
               order(i + 1) = order(i)
               i = i - 1
            END DO
            order(i + 1) = k
         END DO
         coefficients = 0
         row = 0
         key = 0
         DO i = 1, nt
            j = members(order(i))
            IF (row == 0 .OR. ANY(keys(:, order(i)) /= key)) THEN
               row = row + 1
               key = keys(:, order(i))
               shifts(:, row) = key
            END IF
            h = band%fixed_unitary(g, j)
            chi = band%site_characters(h, :)
            IF (band%conjugate_site(j)) chi = CONJG(chi)
            coefficients(row, :) = coefficients(row, :) + band%spin_factor(g, j)*chi
         END DO
         work%blocks(b)%keys = shifts(:, :row)
         work%blocks(b)%coefficients = coefficients(:row, :)
         DEALLOCATE (shifts, order, coefficients)
      END DO
      polynomial = work
      status = 0
   END SUBROUTINE atomic_band_character_polynomial

! **************************************************************************************************
!> \brief Compare complete invariant-torus character functions in the same group/lift convention.
!> \param a First constructed Fourier polynomial
!> \param ia First local column
!> \param b Second constructed Fourier polynomial
!> \param ib Second local column
!> \param equal Whether all collected coefficients agree within tolerance
!> \param tolerance Dimension-scaled absolute coefficient tolerance
!> \param residual Largest dimension-scaled coefficient difference
!> \param status Zero on success; negative for invalid data or different quotient conventions
!> \note This is unitary-character equality, not a site-intertwiner or Bloch-bundle equivalence.
!>       Both polynomials must come from the same cell, operation, coset and spin-lift convention.
! **************************************************************************************************
   SUBROUTINE atomic_character_compare(a, ia, b, ib, equal, tolerance, residual, status)
      TYPE(atomic_character_polynomial_type), INTENT(IN) :: a
      INTEGER, INTENT(IN)                                :: ia
      TYPE(atomic_character_polynomial_type), INTENT(IN) :: b
      INTEGER, INTENT(IN)                                :: ib
      LOGICAL, INTENT(OUT)                               :: equal
      REAL(dp), INTENT(IN)                               :: tolerance
      REAL(dp), INTENT(OUT)                              :: residual
      INTEGER, INTENT(OUT)                               :: status

      COMPLEX(dp)                                        :: difference
      INTEGER                                            :: g, i, j, na, nb, order
      REAL(dp)                                           :: scale

      status = -1
      equal = .FALSE.
      residual = HUGE(1.0_dp)
      IF (.NOT. ieee_is_finite(tolerance) .OR. tolerance <= 0 .OR. tolerance >= 0.01_dp) RETURN
      IF (.NOT. ALLOCATED(a%blocks) .OR. .NOT. ALLOCATED(b%blocks)) RETURN
      IF (.NOT. ALLOCATED(a%dimensions) .OR. .NOT. ALLOCATED(b%dimensions)) RETURN
      IF (ia < 1 .OR. ia > SIZE(a%dimensions) .OR. ib < 1 .OR. ib > SIZE(b%dimensions)) RETURN
      IF (.NOT. valid_polynomial(a) .OR. .NOT. valid_polynomial(b)) RETURN
      IF (SIZE(a%blocks) /= SIZE(b%blocks)) RETURN
      scale = REAL(MAX(a%dimensions(ia), b%dimensions(ib)), dp)
      IF (scale < 1) RETURN
      residual = ABS(REAL(a%dimensions(ia), dp) - REAL(b%dimensions(ib), dp))/scale
      DO g = 1, SIZE(a%blocks)
         IF (a%blocks(g)%operation /= b%blocks(g)%operation) RETURN
         IF (ANY(a%blocks(g)%rotation /= b%blocks(g)%rotation)) RETURN
         IF (ANY(a%blocks(g)%left /= b%blocks(g)%left)) RETURN
         IF (ANY(a%blocks(g)%period /= b%blocks(g)%period)) RETURN
         na = SIZE(a%blocks(g)%keys, 2)
         nb = SIZE(b%blocks(g)%keys, 2)
         i = 1
         j = 1
         DO WHILE (i <= na .OR. j <= nb)
            IF (i > na) THEN
               order = 1
            ELSE IF (j > nb) THEN
               order = -1
            ELSE
               order = key_compare(a%blocks(g)%keys(:, i), b%blocks(g)%keys(:, j))
            END IF
            difference = 0
            IF (order <= 0) THEN
               difference = a%blocks(g)%coefficients(i, ia)
               i = i + 1
            END IF
            IF (order >= 0) THEN
               difference = difference - b%blocks(g)%coefficients(j, ib)
               j = j + 1
            END IF
            residual = MAX(residual, ABS(difference)/scale)
         END DO
      END DO
      equal = residual <= tolerance
      status = 0
   END SUBROUTINE atomic_character_compare

! **************************************************************************************************
!> \brief Validate public polynomial storage before indexed comparison.
!> \param polynomial Fourier blocks to inspect
!> \return Whether shapes, finite coefficients and unique sorted keys are consistent
! **************************************************************************************************
   LOGICAL FUNCTION valid_polynomial(polynomial) RESULT(valid)
      TYPE(atomic_character_polynomial_type), INTENT(IN) :: polynomial

      INTEGER                                            :: g, j, last_operation, n

      valid = .FALSE.
      IF (.NOT. ALLOCATED(polynomial%dimensions) .OR. .NOT. ALLOCATED(polynomial%blocks)) RETURN
      IF (SIZE(polynomial%dimensions) == 0 .OR. SIZE(polynomial%blocks) == 0) RETURN
      IF (ANY(polynomial%dimensions < 1)) RETURN
      last_operation = 0
      DO g = 1, SIZE(polynomial%blocks)
         ASSOCIATE (block => polynomial%blocks(g))
            IF (block%operation <= last_operation .OR. ANY(block%period < 0)) RETURN
            last_operation = block%operation
            IF (.NOT. ALLOCATED(block%keys) .OR. .NOT. ALLOCATED(block%coefficients)) RETURN
            IF (SIZE(block%keys, 1) /= 3) RETURN
            n = SIZE(block%keys, 2)
            IF (SIZE(block%coefficients, 1) /= n) RETURN
            IF (SIZE(block%coefficients, 2) /= SIZE(polynomial%dimensions)) RETURN
            IF (.NOT. ALL(ieee_is_finite(REAL(block%coefficients, dp)))) RETURN
            IF (.NOT. ALL(ieee_is_finite(AIMAG(block%coefficients)))) RETURN
            DO j = 2, n
               IF (key_compare(block%keys(:, j - 1), block%keys(:, j)) >= 0) RETURN
            END DO
            DO j = 1, 3
               IF (block%period(j) == 0) CYCLE
               IF (ANY(block%keys(j, :) < 0) .OR. ANY(block%keys(j, :) >= block%period(j))) RETURN
            END DO
         END ASSOCIATE
      END DO
      valid = .TRUE.
   END FUNCTION valid_polynomial

! **************************************************************************************************
!> \brief Overflow-free lexicographic order of quotient-lattice keys.
!> \param a First key
!> \param b Second key
!> \return Negative, zero or positive according to key order
! **************************************************************************************************
   INTEGER FUNCTION key_compare(a, b) RESULT(order)
      INTEGER(int_8), INTENT(IN)                         :: a(3), b(3)

      INTEGER                                            :: i

      order = 0
      DO i = 1, 3
         IF (a(i) == b(i)) CYCLE
         order = MERGE(-1, 1, a(i) < b(i))
         RETURN
      END DO
   END FUNCTION key_compare

! **************************************************************************************************
!> \brief Induce all local irreps/corepresentations from the stabilizer of a supplied position.
!> \param cell Cartesian lattice vectors as columns
!> \param rotations Integer fractional direct-space rotations of a complete group modulo translations
!> \param translations Fractional translations, preserving the supplied coset representatives
!> \param antiunitary True for an operation followed by physical time reversal
!> \param position Fractional seed position, reduced modulo integer translations
!> \param spinful Include spin-1/2 double-group factors and Theta squared equal to minus one
!> \param band Prepared site-induced representations; empty on failure
!> \param tolerance Absolute geometry and representation tolerance
!> \param residual Largest orbit closure or local representation residual
!> \param status Zero on success, negative for invalid geometry or inconsistent algebra
!> \param full_spin Optional physical spin lifts in original full-group operation order; empty on failure
!> \note Site translations are adjusted to fix the seed exactly; onsite factors contain only spin.
!>       An antiunitary site stabilizer requires corepresentations, not arbitrary unitary irreps.
!>       These induced band representations need not be elementary. No Wyckoff completeness,
!>       primitive-cell EBR catalogue, canonical labels or topological classification is asserted.
! **************************************************************************************************
   SUBROUTINE atomic_band_construct(cell, rotations, translations, antiunitary, position, spinful, &
                                    band, tolerance, residual, status, full_spin)
      REAL(dp), INTENT(IN)                               :: cell(3, 3)
      INTEGER, INTENT(IN)                                :: rotations(:, :, :)
      REAL(dp), INTENT(IN)                               :: translations(:, :)
      LOGICAL, INTENT(IN)                                :: antiunitary(:)
      REAL(dp), INTENT(IN)                               :: position(3)
      LOGICAL, INTENT(IN)                                :: spinful
      TYPE(atomic_band_type), INTENT(OUT)                :: band
      REAL(dp), INTENT(IN)                               :: tolerance
      REAL(dp), INTENT(OUT)                              :: residual
      INTEGER, INTENT(OUT)                               :: status
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT), OPTIONAL    :: full_spin(:, :, :)

      COMPLEX(dp), ALLOCATABLE                           :: characters(:, :), factors(:, :), &
                                                            site_factor(:, :), spin_lifts(:, :, :)
      INTEGER                                            :: found, g, h, i, j, n, nh, no, nu, t
      INTEGER, ALLOCATABLE :: image(:), inverse(:), products(:, :), representatives(:), &
         site_index(:), site_members(:), site_product(:, :), unit_index(:), unit_product(:, :), &
         unitary(:)
      REAL(dp)                                           :: delta(3), error, lattice(3, 3), &
                                                            lattice_inverse(3, 3), seed(3), &
                                                            TARGET(3)
      REAL(dp), ALLOCATABLE                              :: positions(:, :)
      TYPE(atomic_band_type)                             :: work
      TYPE(corepresentation_table_type)                  :: coreps

      status = -1
      residual = HUGE(1.0_dp)
      IF (.NOT. ieee_is_finite(tolerance)) RETURN
      IF (tolerance <= 0.0_dp .OR. tolerance >= 0.01_dp) RETURN
      IF (.NOT. ALL(ieee_is_finite(position))) RETURN
      IF (MAXVAL(ABS(position)) > tolerance/(32*EPSILON(1.0_dp))) RETURN
      CALL little_group_factor_system(cell, rotations, translations, antiunitary, [0.0_dp, 0.0_dp, 0.0_dp], &
                                      spinful, work%mapping, products, factors, tolerance, status, spin_lifts=spin_lifts)
      IF (status /= 0) RETURN
      n = SIZE(work%mapping)
      work%rotations = rotations(:, :, work%mapping)
      work%antiunitary = antiunitary(work%mapping)
      seed = MODULO(position, 1.0_dp)
      ALLOCATE (positions(3, n), representatives(n), image(n), site_index(n), inverse(n))
      residual = 0.0_dp
      no = 0
      status = -2
      DO g = 1, n
         TARGET = MODULO(MATMUL(REAL(work%rotations(:, :, g), dp), seed) + translations(:, work%mapping(g)), 1.0_dp)
         found = 0
         DO j = 1, no
            delta = TARGET - positions(:, j)
            error = MAXVAL(ABS(delta - ANINT(delta)))
            IF (error > tolerance) CYCLE
            IF (found /= 0) RETURN
            found = j
            residual = MAX(residual, error)
         END DO
         IF (found == 0) THEN
            no = no + 1
            positions(:, no) = TARGET
            representatives(no) = g
            found = no
         END IF
         image(g) = found
         inverse(g) = FINDLOC(products(g, :), 1, DIM=1)
      END DO
      nh = COUNT(image == 1)
      IF (nh*no /= n) RETURN
      DO j = 1, no
         IF (COUNT(image == j) /= nh) RETURN
      END DO
      site_members = PACK([(g, g=1, n)], image == 1)
      site_index = 0
      DO h = 1, nh
         site_index(site_members(h)) = h
      END DO
      ALLOCATE (site_product(nh, nh))
      DO h = 1, nh
         site_product(:, h) = site_index(products(site_members, site_members(h)))
      END DO
      site_factor = factors(site_members, site_members)
      CALL validate_factor_system(site_product, site_factor, work%antiunitary(site_members), tolerance, error, status)
      IF (status /= 0) RETURN
      residual = MAX(residual, error)
      unitary = PACK([(h, h=1, nh)],.NOT. work%antiunitary(site_members))
      nu = SIZE(unitary)
      ALLOCATE (unit_index(nh), unit_product(nu, nu))
      unit_index = 0
      DO h = 1, nu
         unit_index(unitary(h)) = h
      END DO
      DO h = 1, nu
         unit_product(:, h) = unit_index(site_product(unitary, unitary(h)))
      END DO
      CALL projective_characters(unit_product, site_factor(unitary, unitary), characters, tolerance, error, status)
      IF (status /= 0) RETURN
      residual = MAX(residual, error)
      IF (nu /= nh) THEN
         CALL projective_corepresentations(site_product, site_factor, work%antiunitary(site_members), &
                                           characters, coreps, tolerance, error, status)
         IF (status /= 0) RETURN
         residual = MAX(residual, error)
         work%site_characters = MATMUL(characters, CMPLX(coreps%restriction, 0, dp))
      ELSE
         CALL MOVE_ALLOC(characters, work%site_characters)
      END IF
      work%dimensions = NINT(REAL(work%site_characters(1, :), dp))
      work%site_mapping = work%mapping(site_members(unitary))
      work%site_full_mapping = work%mapping(site_members)
      work%spinful = spinful
      ALLOCATE (work%site_product(nh, nh), work%site_factor(nh, nh), work%site_spin(2, 2, nh), &
                work%site_cartesian(3, 3, nh), work%site_antiunitary(nh), work%site_lattice_shifts(3, nh))
      work%site_product(:, :) = site_product
      work%site_factor(:, :) = site_factor
      work%site_spin(:, :, :) = spin_lifts(:, :, site_members)
      work%site_antiunitary(:) = work%antiunitary(site_members)
      lattice = cell/MAXVAL(ABS(cell))
      lattice_inverse = inv_3x3(lattice)
      DO h = 1, nh
         g = work%site_full_mapping(h)
         work%site_lattice_shifts(:, h) = ANINT(MATMUL(REAL(rotations(:, :, g), dp), seed) + translations(:, g) - seed)
         work%site_cartesian(:, :, h) = MATMUL(MATMUL(lattice, REAL(work%rotations(:, :, site_members(h)), dp)), &
                                               lattice_inverse)
      END DO
      work%positions = positions(:, :no)
      work%conjugate_site = work%antiunitary(representatives(:no))
      ALLOCATE (work%fixed_unitary(n, no), work%spin_factor(n, no), work%shifts(3, n, no))
      work%fixed_unitary = 0
      work%spin_factor = 0.0_dp
      work%shifts = 0.0_dp
      status = -2
      DO j = 1, no
         t = representatives(j)
         DO g = 1, n
            IF (work%antiunitary(g)) CYCLE
            TARGET = MATMUL(REAL(work%rotations(:, :, g), dp), positions(:, j)) + translations(:, work%mapping(g))
            found = 0
            DO i = 1, no
               delta = TARGET - positions(:, i)
               error = MAXVAL(ABS(delta - ANINT(delta)))
               IF (error > tolerance) CYCLE
               IF (found /= 0) RETURN
               found = i
               residual = MAX(residual, error)
            END DO
            IF (found == 0) RETURN
            IF (found /= j) CYCLE
            h = products(inverse(t), products(g, t))
            IF (site_index(h) == 0) RETURN
            i = unit_index(site_index(h))
            IF (i == 0) RETURN
            work%fixed_unitary(g, j) = i
            work%shifts(:, g, j) = ANINT(TARGET - positions(:, j))
            ! g*t_j = T_L*t_j*h; antiunitary t_j conjugates the local character of h.
            work%spin_factor(g, j) = factors(g, t)/factors(t, h)
         END DO
      END DO
      band = work
      IF (PRESENT(full_spin)) THEN
         ALLOCATE (full_spin(2, 2, n))
         full_spin(:, :, work%mapping) = spin_lifts
      END IF
      status = 0
   END SUBROUTINE atomic_band_construct

! **************************************************************************************************
!> \brief Evaluate induced characters on the unitary little group at a supplied k point.
!> \param band Successfully constructed site-induced representations
!> \param k Fractional reciprocal point in the cell convention used to construct band
!> \param mapping Original input operation indices of unitary little-group members, identity first
!> \param characters One column per local irrep/corepresentation; unallocated on failure
!> \param tolerance Absolute fractional reciprocal tolerance
!> \param status Zero on success, negative for invalid data
!> \note Characters sum only orbit sites fixed modulo an integer translation by each operation.
!>       The convention is exp(-2*pi*i*k.L), consistent with little_group_factor_system.
! **************************************************************************************************
   SUBROUTINE atomic_band_characters(band, k, mapping, characters, tolerance, status)
      TYPE(atomic_band_type), INTENT(IN)                 :: band
      REAL(dp), INTENT(IN)                               :: k(3)
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: mapping(:)
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: characters(:, :)
      REAL(dp), INTENT(IN)                               :: tolerance
      INTEGER, INTENT(OUT)                               :: status

      COMPLEX(dp)                                        :: phase
      COMPLEX(dp), ALLOCATABLE                           :: chi(:)
      INTEGER                                            :: g, h, i, j, n, nu
      INTEGER, ALLOCATABLE                               :: members(:)
      REAL(dp)                                           :: angle, delta(3), kp(3)

      status = -1
      IF (.NOT. ALLOCATED(band%mapping)) RETURN
      IF (.NOT. ieee_is_finite(tolerance)) RETURN
      IF (tolerance <= 0.0_dp .OR. tolerance >= 0.01_dp) RETURN
      IF (.NOT. ALL(ieee_is_finite(k))) RETURN
      IF (MAXVAL(ABS(k)) > tolerance/(32*EPSILON(1.0_dp))) RETURN
      kp = MODULO(k, 1.0_dp)
      n = SIZE(band%mapping)
      ALLOCATE (members(n))
      nu = 0
      DO g = 1, n
         IF (band%antiunitary(g)) CYCLE
         delta = MATMUL(REAL(TRANSPOSE(band%rotations(:, :, g)), dp), kp) - kp
         IF (MAXVAL(ABS(delta - ANINT(delta))) > tolerance) CYCLE
         nu = nu + 1
         members(nu) = g
      END DO
      ALLOCATE (mapping(nu), characters(nu, SIZE(band%dimensions)), chi(SIZE(band%dimensions)))
      mapping = band%mapping(members(:nu))
      characters = 0.0_dp
      DO i = 1, nu
         g = members(i)
         DO j = 1, SIZE(band%positions, 2)
            h = band%fixed_unitary(g, j)
            IF (h == 0) CYCLE
            angle = -twopi*MODULO(DOT_PRODUCT(kp, band%shifts(:, g, j)), 1.0_dp)
            phase = CMPLX(COS(angle), SIN(angle), dp)*band%spin_factor(g, j)
            chi = band%site_characters(h, :)
            IF (band%conjugate_site(j)) chi = CONJG(chi)
            characters(i, :) = characters(i, :) + phase*chi
         END DO
      END DO
      status = 0
   END SUBROUTINE atomic_band_characters

! **************************************************************************************************
!> \brief Get orbit and onsite unitary characters without exposing cached induction data.
!> \param band Constructed site-induced representations
!> \param positions Orbit positions in the supplied cell; unallocated if band is empty
!> \param mapping Original operation indices of the unitary site stabilizer
!> \param characters Local unitary irreps or restrictions of onsite corepresentations
!> \param orbit_index Optional orbit site, retaining the seed's representation-column labels
! **************************************************************************************************
   SUBROUTINE atomic_band_site_data(band, positions, mapping, characters, orbit_index)
      TYPE(atomic_band_type), INTENT(IN)                 :: band
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: positions(:, :)
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: mapping(:)
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: characters(:, :)
      INTEGER, INTENT(IN), OPTIONAL                      :: orbit_index

      INTEGER                                            :: g, h, i, j
      INTEGER, ALLOCATABLE                               :: members(:)

      IF (.NOT. ALLOCATED(band%dimensions)) RETURN
      j = 1
      IF (PRESENT(orbit_index)) j = orbit_index
      IF (j < 1 .OR. j > SIZE(band%positions, 2)) RETURN
      positions = band%positions
      IF (j == 1) THEN
         mapping = band%site_mapping
         characters = band%site_characters
      ELSE
         members = PACK([(g, g=1, SIZE(band%mapping))], band%fixed_unitary(:, j) /= 0)
         mapping = band%mapping(members)
         ALLOCATE (characters(SIZE(members), SIZE(band%dimensions)))
         DO i = 1, SIZE(members)
            g = members(i)
            h = band%fixed_unitary(g, j)
            characters(i, :) = band%site_characters(h, :)
            IF (band%conjugate_site(j)) characters(i, :) = CONJG(characters(i, :))
            characters(i, :) = band%spin_factor(g, j)*characters(i, :)
         END DO
      END IF
   END SUBROUTINE atomic_band_site_data

! **************************************************************************************************
!> \brief Copy the exact onsite convention used by the induced band columns.
!> \param band Constructed induced representations
!> \param site Compact reference data; unallocated components if band is empty
!> \note Mapping includes antiunitary operations; character rows contain only unitary operations.
!>       Lattice shifts obey W*q+tau=q+L for the stored seed reduced modulo the input lattice.
!>       The onsite action fixes q exactly and has only spin factors, not the Bloch factor of tau.
! **************************************************************************************************
   SUBROUTINE atomic_band_site_snapshot(band, site)
      TYPE(atomic_band_type), INTENT(IN)                 :: band
      TYPE(atomic_site_data_type), INTENT(OUT)           :: site

      IF (.NOT. ALLOCATED(band%site_product)) RETURN
      site%spinful = band%spinful
      site%mapping = band%site_full_mapping
      site%product = band%site_product
      site%antiunitary = band%site_antiunitary
      site%cartesian = band%site_cartesian
      site%lattice_shifts = band%site_lattice_shifts
      site%spin = band%site_spin
      site%factor = band%site_factor
      site%characters = band%site_characters
   END SUBROUTINE atomic_band_site_snapshot

! **************************************************************************************************
!> \brief Estimate retained onsite array storage in MiB, excluding allocation descriptors.
!> \param site Onsite snapshot, possibly empty
!> \return Array storage using at least eight bytes for integer and logical entries
! **************************************************************************************************
   FUNCTION atomic_site_storage(site) RESULT(memory)
      TYPE(atomic_site_data_type), INTENT(IN)            :: site
      REAL(dp)                                           :: memory

      memory = 0.0_dp
      IF (ALLOCATED(site%mapping)) memory = memory + 8.0_dp*SIZE(site%mapping)
      IF (ALLOCATED(site%product)) memory = memory + 8.0_dp*SIZE(site%product)
      IF (ALLOCATED(site%antiunitary)) memory = memory + 8.0_dp*SIZE(site%antiunitary)
      IF (ALLOCATED(site%cartesian)) memory = memory + 8.0_dp*SIZE(site%cartesian)
      IF (ALLOCATED(site%lattice_shifts)) memory = memory + 8.0_dp*SIZE(site%lattice_shifts)
      IF (ALLOCATED(site%spin)) memory = memory + 16.0_dp*SIZE(site%spin)
      IF (ALLOCATED(site%factor)) memory = memory + 16.0_dp*SIZE(site%factor)
      IF (ALLOCATED(site%characters)) memory = memory + 16.0_dp*SIZE(site%characters)
      memory = memory/1024.0_dp**2
   END FUNCTION atomic_site_storage

! **************************************************************************************************
!> \brief Identify the seed's onsite columns by explicit reference axes, spin lifts and characters.
!> \param band Constructed local representations, including full corepresentations when needed
!> \param reference_rotations Cartesian polar rotations of the reference site group, identity first
!> \param reference_antiunitary Antiunitary flags in reference operation order
!> \param reference_spin Spin lifts including i*sigma_y for antiunitary operations; ignored for scalar bands
!> \param reference_characters Complete unitary characters, or restrictions of complete corepresentations;
!>                            rows follow the unitary reference operations in ascending index order
!> \param orientation Orthogonal Cartesian map from the reference axes to the band's axes
!> \param mapping Native onsite column to reference column; unallocated on failure
!> \param tolerance Absolute geometry/spin tolerance and dimension-scaled character tolerance
!> \param residual Largest geometry, lift or character residual
!> \param status Zero on success; negative for invalid or inconsistent data
!> \param operation_mapping Optional reference operation for each native onsite operation
!> \param lift_phases Optional phase in D_native(g)=phase(g)*D_reference(operation_mapping(g))
!> \note The axes and spin lifts are supplied, not fitted from equal dimensions or induced k labels.
!>       This labels onsite representations, not elementary band representations or band connectivity.
! **************************************************************************************************
   SUBROUTINE atomic_band_site_match(band, reference_rotations, reference_antiunitary, reference_spin, &
                                     reference_characters, orientation, mapping, tolerance, residual, status, &
                                     operation_mapping, lift_phases)
      TYPE(atomic_band_type), INTENT(IN)                 :: band
      REAL(dp), INTENT(IN)                               :: reference_rotations(:, :, :)
      LOGICAL, INTENT(IN)                                :: reference_antiunitary(:)
      COMPLEX(dp), INTENT(IN)                            :: reference_spin(:, :, :), &
                                                            reference_characters(:, :)
      REAL(dp), INTENT(IN)                               :: orientation(3, 3)
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: mapping(:)
      REAL(dp), INTENT(IN)                               :: tolerance
      REAL(dp), INTENT(OUT)                              :: residual
      INTEGER, INTENT(OUT)                               :: status
      INTEGER, ALLOCATABLE, INTENT(OUT), OPTIONAL        :: operation_mapping(:)
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT), OPTIONAL    :: lift_phases(:)

      COMPLEX(dp)                                        :: composed(2, 2), orient_spin(2, 2), &
                                                            phase, rotated(2, 2)
      COMPLEX(dp), ALLOCATABLE                           :: phases(:), ref_factor(:, :)
      INTEGER                                            :: found, g, h, n, q
      INTEGER, ALLOCATABLE                               :: inverse(:), operations(:), &
                                                            ref_product(:, :)
      REAL(dp)                                           :: cart(3, 3), closest, error

      residual = HUGE(1.0_dp)
      status = -1
      IF (.NOT. ALLOCATED(band%site_product)) RETURN
      n = SIZE(band%site_antiunitary)
      IF (ANY(SHAPE(reference_rotations) /= [3, 3, n])) RETURN
      IF (SIZE(reference_antiunitary) /= n .OR. ANY(SHAPE(reference_spin) /= [2, 2, n])) RETURN
      IF (COUNT(.NOT. reference_antiunitary) /= SIZE(band%site_characters, 1)) RETURN
      IF (ANY(SHAPE(reference_characters) /= SHAPE(band%site_characters))) RETURN
      IF (.NOT. ieee_is_finite(tolerance) .OR. tolerance <= 0.0_dp .OR. tolerance >= 0.01_dp) RETURN
      IF (.NOT. ALL(ieee_is_finite(reference_rotations)) .OR. .NOT. ALL(ieee_is_finite(orientation))) RETURN
      IF (MAXVAL(ABS(reference_rotations)) > 2.0_dp .OR. MAXVAL(ABS(orientation)) > 2.0_dp) RETURN
      IF (.NOT. ALL(ieee_is_finite(REAL(reference_spin, dp))) .OR. &
          .NOT. ALL(ieee_is_finite(AIMAG(reference_spin)))) RETURN
      IF (MAXVAL(ABS(REAL(reference_spin, dp))) > 2.0_dp .OR. MAXVAL(ABS(AIMAG(reference_spin))) > 2.0_dp) RETURN
      CALL spin_rotation(orientation, orient_spin, tolerance, status)
      IF (status /= 0) RETURN
      ALLOCATE (operations(n), inverse(n), phases(n), ref_product(n, n), ref_factor(n, n))
      inverse = 0
      phases = 1.0_dp
      residual = 0.0_dp
      status = -2
      DO g = 1, n
         found = 0
         closest = HUGE(1.0_dp)
         DO h = 1, n
            IF (band%site_antiunitary(g) .NEQV. reference_antiunitary(h)) CYCLE
            cart = MATMUL(MATMUL(orientation, reference_rotations(:, :, h)), TRANSPOSE(orientation))
            error = MAXVAL(ABS(cart - band%site_cartesian(:, :, g)))
            closest = MIN(closest, error)
            IF (error > tolerance) CYCLE
            IF (found /= 0 .OR. inverse(h) /= 0) RETURN
            found = h
            residual = MAX(residual, error)
         END DO
         IF (found == 0) THEN
            residual = MAX(residual, closest)
            RETURN
         END IF
         operations(g) = found
         inverse(found) = g
         IF (band%spinful) THEN
            IF (band%site_antiunitary(g)) THEN
               rotated = MATMUL(MATMUL(orient_spin, reference_spin(:, :, found)), TRANSPOSE(orient_spin))
            ELSE
               rotated = MATMUL(MATMUL(orient_spin, reference_spin(:, :, found)), CONJG(TRANSPOSE(orient_spin)))
            END IF
            phase = SUM(CONJG(rotated)*band%site_spin(:, :, g))/2
            IF (ABS(phase) < 0.5_dp) RETURN
            phases(g) = phase/ABS(phase)
            error = MAXVAL(ABS(band%site_spin(:, :, g) - phases(g)*rotated))
            residual = MAX(residual, error)
            IF (error > tolerance) RETURN
         END IF
      END DO
      ref_factor = 1.0_dp
      DO h = 1, n
         DO g = 1, n
            q = operations(band%site_product(inverse(g), inverse(h)))
            ref_product(g, h) = q
            IF (.NOT. band%spinful) CYCLE
            IF (reference_antiunitary(g)) THEN
               composed = MATMUL(reference_spin(:, :, g), CONJG(reference_spin(:, :, h)))
            ELSE
               composed = MATMUL(reference_spin(:, :, g), reference_spin(:, :, h))
            END IF
            phase = SUM(CONJG(reference_spin(:, :, q))*composed)/2
            IF (ABS(phase) < 0.5_dp) RETURN
            ref_factor(g, h) = phase/ABS(phase)
            error = MAXVAL(ABS(composed - ref_factor(g, h)*reference_spin(:, :, q)))
            residual = MAX(residual, error)
            IF (error > tolerance) RETURN
         END DO
      END DO
      CALL character_table_match(band%site_product, band%site_factor, band%site_antiunitary, band%site_characters, &
                                 ref_product, ref_factor, reference_antiunitary, reference_characters, &
                                 operations, phases, mapping, tolerance, error, status)
      residual = MAX(residual, error)
      IF (status == 0) THEN
         IF (PRESENT(operation_mapping)) CALL MOVE_ALLOC(operations, operation_mapping)
         IF (PRESENT(lift_phases)) CALL MOVE_ALLOC(phases, lift_phases)
      END IF
   END SUBROUTINE atomic_band_site_match
END MODULE topology_band_representations
