!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Band compatibility along reciprocal-space segments, including antiunitary symmetries.
!> \author Thomas D. Kuehne, tkuehne@cp2k.org
! **************************************************************************************************
MODULE topology_compatibility
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: twopi
   USE topology_corepresentations,      ONLY: corepresentation_multiplicities,&
                                              corepresentation_table_type,&
                                              projective_corepresentations
   USE topology_little_group,           ONLY: projective_characters,&
                                              validate_factor_system
   USE topology_symmetry,               ONLY: character_multiplicities

   IMPLICIT NONE
   PRIVATE

   !> Validated little-group data; character rows follow the unitary operations in full-group order.
   TYPE, PUBLIC :: little_group_data_type
      REAL(dp) :: k(3) = 0.0_dp, cell(3, 3) = 0.0_dp
      LOGICAL :: spinful = .FALSE.
      LOGICAL, ALLOCATABLE :: antiunitary(:)
      INTEGER, ALLOCATABLE :: rotations(:, :, :), PRODUCT(:, :), multiplicity(:)
      REAL(dp), ALLOCATABLE :: translations(:, :)
      COMPLEX(dp), ALLOCATABLE :: factors(:, :), characters(:, :)
   END TYPE little_group_data_type

   !> Restriction matrices and band counts in a common line-character gauge; local irrep labels.
   TYPE, PUBLIC :: little_group_edge_type
      INTEGER, ALLOCATABLE :: left(:), right(:), left_branch(:, :), right_branch(:, :)
      INTEGER, ALLOCATABLE :: left_count(:), right_count(:)
      INTEGER, ALLOCATABLE :: anti_left(:), anti_right(:), left_corep_branch(:, :), right_corep_branch(:, :)
      INTEGER, ALLOCATABLE :: left_corep_count(:), right_corep_count(:)
      COMPLEX(dp), ALLOCATABLE :: characters(:, :)
      TYPE(corepresentation_table_type) :: coreps
   END TYPE little_group_edge_type

   PUBLIC :: little_group_compatibility, little_group_storage
CONTAINS

! **************************************************************************************************
!> \brief Restrict both endpoints to operations fixing the entire straight segment.
!> \param left Validated little-group data at the segment origin
!> \param right Validated little-group data at the stored destination point
!> \param reciprocal_shift Integer reciprocal image of the destination; it must not be discarded
!> \param edge Common line irreps, restrictions and band counts, allocated on success
!> \param tolerance Absolute fractional-coordinate, factor and character tolerance
!> \param residual Largest common-factor or character-generation residual
!> \param status Zero for valid results, including incompatible band counts; negative for bad data
!> \note A mismatch is a compatibility obstruction, not an algebra failure. Equal counts are
!>       necessary, not sufficient, for a gap along the segment or topological triviality.
!>       A common antiunitary coset adds corepresentation subduction to the unitary restrictions.
! **************************************************************************************************
   SUBROUTINE little_group_compatibility(left, right, reciprocal_shift, edge, tolerance, residual, status)
      TYPE(little_group_data_type), INTENT(IN)           :: left, right
      INTEGER, INTENT(IN)                                :: reciprocal_shift(3)
      TYPE(little_group_edge_type), INTENT(OUT)          :: edge
      REAL(dp), INTENT(IN)                               :: tolerance
      REAL(dp), INTENT(OUT)                              :: residual
      INTEGER, INTENT(OUT)                               :: status

      COMPLEX(dp)                                        :: phase_left, phase_right
      COMPLEX(dp), ALLOCATABLE                           :: chi(:), factors(:, :), other(:, :), &
                                                            phases(:, :)
      INTEGER                                            :: found, g, gh, h, i, j, n, ni, nl, nr, nu
      INTEGER, ALLOCATABLE :: il(:), inverse(:), ir(:), products(:, :), unit_index(:), &
         unit_left(:), unit_products(:, :), unit_right(:), unitary(:)
      LOGICAL, ALLOCATABLE                               :: anti(:), anti_left(:), anti_right(:)
      REAL(dp)                                           :: angle, delta(3), direction(3), error, &
                                                            target_k(3)

      status = -1
      residual = HUGE(1.0_dp)
      IF (.NOT. ieee_is_finite(tolerance)) RETURN
      IF (tolerance <= 0.0_dp .OR. tolerance >= 0.01_dp) RETURN
      IF (.NOT. valid_point(left) .OR. .NOT. valid_point(right)) RETURN
      IF (left%spinful .NEQV. right%spinful) RETURN
      IF (MAXVAL(ABS(left%cell - right%cell)) > tolerance*MAXVAL(ABS(left%cell))) RETURN
      nl = SIZE(left%product, 1)
      nr = SIZE(right%product, 1)
      ALLOCATE (anti_left(nl), anti_right(nr), unit_left(nl), unit_right(nr))
      anti_left = .FALSE.
      anti_right = .FALSE.
      IF (ALLOCATED(left%antiunitary)) anti_left(:) = left%antiunitary
      IF (ALLOCATED(right%antiunitary)) anti_right(:) = right%antiunitary
      unit_left(:) = 0
      unit_right(:) = 0
      DO g = 1, nl
         IF (.NOT. anti_left(g)) unit_left(g) = COUNT(.NOT. anti_left(:g))
      END DO
      DO g = 1, nr
         IF (.NOT. anti_right(g)) unit_right(g) = COUNT(.NOT. anti_right(:g))
      END DO
      CALL validate_factor_system(left%product, left%factors, anti_left, tolerance, error, status)
      IF (status /= 0) RETURN
      residual = error
      CALL validate_factor_system(right%product, right%factors, anti_right, tolerance, error, status)
      IF (status /= 0) RETURN
      residual = MAX(residual, error)
      status = -2
      target_k = right%k + REAL(reciprocal_shift, dp)
      direction = target_k - left%k
      ALLOCATE (il(nl), ir(nl), inverse(nl))
      inverse(:) = 0
      n = 0
      DO g = 1, nl
         delta = MERGE(-1.0_dp, 1.0_dp, anti_left(g))* &
                 MATMUL(REAL(TRANSPOSE(left%rotations(:, :, g)), dp), direction) - direction
         ! An integer change is not enough: the operation must fix every interior point.
         IF (MAXVAL(ABS(delta)) > tolerance) CYCLE
         found = 0
         DO h = 1, nr
            IF (anti_left(g) .NEQV. anti_right(h)) CYCLE
            IF (ANY(left%rotations(:, :, g) /= right%rotations(:, :, h))) CYCLE
            delta = left%translations(:, g) - right%translations(:, h)
            IF (MAXVAL(ABS(delta - ANINT(delta))) > tolerance) CYCLE
            IF (found /= 0) RETURN
            found = h
         END DO
         IF (found == 0) RETURN
         n = n + 1
         il(n) = g
         ir(n) = found
         inverse(g) = n
      END DO
      IF (n < 1) RETURN
      IF (il(1) /= 1 .OR. ir(1) /= 1) RETURN
      ALLOCATE (anti(n), unit_index(n))
      anti(:) = anti_left(il(:n))
      nu = COUNT(.NOT. anti)
      ALLOCATE (unitary(nu))
      unitary(:) = PACK([(g, g=1, n)],.NOT. anti)
      unit_index = 0
      DO g = 1, nu
         unit_index(unitary(g)) = g
      END DO
      ALLOCATE (products(n, n), factors(n, n), other(n, n), phases(n, 2), chi(nu))
      DO g = 1, n
         angle = twopi*DOT_PRODUCT(left%k, left%translations(:, il(g)))
         phases(g, 1) = CMPLX(COS(angle), SIN(angle), dp)
         angle = twopi*DOT_PRODUCT(target_k, right%translations(:, ir(g)))
         phases(g, 2) = CMPLX(COS(angle), SIN(angle), dp)
      END DO
      ! exp(+2*pi*i*k.tau_g) removes the continuously varying translational character.
      DO h = 1, n
         DO g = 1, n
            gh = inverse(left%product(il(g), il(h)))
            IF (gh == 0) RETURN
            IF (right%product(ir(g), ir(h)) /= ir(gh)) RETURN
            products(g, h) = gh
            phase_left = phases(h, 1)
            phase_right = phases(h, 2)
            IF (anti(g)) THEN
               phase_left = CONJG(phase_left)
               phase_right = CONJG(phase_right)
            END IF
            factors(g, h) = phases(g, 1)*phase_left*CONJG(phases(gh, 1))*left%factors(il(g), il(h))
            other(g, h) = phases(g, 2)*phase_right*CONJG(phases(gh, 2))*right%factors(ir(g), ir(h))
         END DO
      END DO
      residual = MAX(residual, MAXVAL(ABS(factors - other)))
      IF (residual > tolerance) RETURN
      ALLOCATE (unit_products(nu, nu))
      DO h = 1, nu
         unit_products(:, h) = unit_index(products(unitary, unitary(h)))
      END DO
      CALL projective_characters(unit_products, factors(unitary, unitary), edge%characters, tolerance, error, status)
      IF (status /= 0) RETURN
      residual = MAX(residual, error)
      ni = SIZE(edge%characters, 2)
      ALLOCATE (edge%left(nu), edge%right(nu), edge%left_count(ni), edge%right_count(ni))
      ALLOCATE (edge%left_branch(ni, SIZE(left%characters, 2)), edge%right_branch(ni, SIZE(right%characters, 2)))
      edge%left(:) = il(unitary)
      edge%right(:) = ir(unitary)
      DO i = 1, SIZE(left%characters, 2)
         chi(:) = phases(unitary, 1)*left%characters(unit_left(edge%left), i)
         CALL character_multiplicities(chi, edge%characters, edge%left_branch(:, i), tolerance, status)
         IF (status /= 0) RETURN
      END DO
      DO j = 1, SIZE(right%characters, 2)
         chi(:) = phases(unitary, 2)*right%characters(unit_right(edge%right), j)
         CALL character_multiplicities(chi, edge%characters, edge%right_branch(:, j), tolerance, status)
         IF (status /= 0) RETURN
      END DO
      edge%left_count(:) = MATMUL(edge%left_branch, left%multiplicity)
      edge%right_count(:) = MATMUL(edge%right_branch, right%multiplicity)
      IF (n > nu) THEN
         CALL projective_corepresentations(products, factors, anti, edge%characters, edge%coreps, &
                                           tolerance, error, status)
         IF (status /= 0) RETURN
         residual = MAX(residual, error)
         ALLOCATE (edge%anti_left(n - nu), edge%anti_right(n - nu))
         edge%anti_left(:) = PACK(il(:n), anti)
         edge%anti_right(:) = PACK(ir(:n), anti)
      END IF
      IF (ANY(anti_left)) THEN
         CALL corepresentation_branching(left, anti_left, edge%left_branch, edge%coreps, &
                                         edge%left_corep_branch, edge%left_corep_count, tolerance, error, status)
         IF (status /= 0) RETURN
         residual = MAX(residual, error)
      END IF
      IF (ANY(anti_right)) THEN
         CALL corepresentation_branching(right, anti_right, edge%right_branch, edge%coreps, &
                                         edge%right_corep_branch, edge%right_corep_count, tolerance, error, status)
         IF (status /= 0) RETURN
         residual = MAX(residual, error)
      END IF
      status = 0
   END SUBROUTINE little_group_compatibility

! **************************************************************************************************
!> \brief Subduce endpoint corepresentations into line corepresentations or unitary line irreps.
!> \param point Validated endpoint with an antiunitary coset
!> \param antiunitary Endpoint operation flags
!> \param unitary_branch Unitary endpoint-to-line branching matrix in the common line gauge
!> \param line Corepresentation table of the rephased line factor, unallocated for a unitary line
!> \param branch Line-by-endpoint corepresentation multiplicities
!> \param counts Selected-band line corepresentation counts
!> \param tolerance Character and factor tolerance
!> \param residual Endpoint classification residual
!> \param status Zero on success; negative for malformed data or noninteger restrictions
! **************************************************************************************************
   SUBROUTINE corepresentation_branching(point, antiunitary, unitary_branch, line, branch, counts, &
                                         tolerance, residual, status)
      TYPE(little_group_data_type), INTENT(IN)           :: point
      LOGICAL, INTENT(IN)                                :: antiunitary(:)
      INTEGER, INTENT(IN)                                :: unitary_branch(:, :)
      TYPE(corepresentation_table_type), INTENT(IN)      :: line
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: branch(:, :), counts(:)
      REAL(dp), INTENT(IN)                               :: tolerance
      REAL(dp), INTENT(OUT)                              :: residual
      INTEGER, INTENT(OUT)                               :: status

      INTEGER                                            :: i, n
      INTEGER, ALLOCATABLE                               :: endpoint_counts(:), unitary_counts(:)
      TYPE(corepresentation_table_type)                  :: endpoint

      CALL projective_corepresentations(point%product, point%factors, antiunitary, point%characters, endpoint, &
                                        tolerance, residual, status)
      IF (status /= 0) RETURN
      n = SIZE(endpoint%indicator)
      ALLOCATE (endpoint_counts(n), unitary_counts(SIZE(unitary_branch, 1)))
      CALL corepresentation_multiplicities(endpoint, point%multiplicity, endpoint_counts, status)
      IF (status /= 0) RETURN
      IF (.NOT. ALLOCATED(line%indicator)) THEN
         ALLOCATE (branch(SIZE(unitary_branch, 1), n), counts(SIZE(unitary_branch, 1)))
         branch(:, :) = MATMUL(unitary_branch, endpoint%restriction)
         counts(:) = MATMUL(branch, endpoint_counts)
         RETURN
      END IF
      ALLOCATE (branch(SIZE(line%indicator), n), counts(SIZE(line%indicator)))
      DO i = 1, n
         unitary_counts(:) = MATMUL(unitary_branch, endpoint%restriction(:, i))
         CALL corepresentation_multiplicities(line, unitary_counts, branch(:, i), status)
         IF (status /= 0) RETURN
      END DO
      counts(:) = MATMUL(branch, endpoint_counts)
   END SUBROUTINE corepresentation_branching

! **************************************************************************************************
!> \brief Structural checks; the producer must also validate the geometric group and band action.
!> \param point Group data with the identity first; absent antiunitary flags mean a unitary group
!> \return Whether the data have consistent finite dimensions and nonnegative band counts
! **************************************************************************************************
   LOGICAL FUNCTION valid_point(point)
      TYPE(little_group_data_type), INTENT(IN)           :: point

      INTEGER                                            :: n, ni, nu

      valid_point = .FALSE.
      IF (.NOT. ALLOCATED(point%rotations) .OR. .NOT. ALLOCATED(point%translations)) RETURN
      IF (.NOT. ALLOCATED(point%product) .OR. .NOT. ALLOCATED(point%factors)) RETURN
      IF (.NOT. ALLOCATED(point%characters) .OR. .NOT. ALLOCATED(point%multiplicity)) RETURN
      n = SIZE(point%product, 1)
      nu = n
      IF (ALLOCATED(point%antiunitary)) THEN
         IF (SIZE(point%antiunitary) /= n) RETURN
         nu = COUNT(.NOT. point%antiunitary)
      END IF
      ni = SIZE(point%multiplicity)
      IF (n < 1 .OR. ni < 1) RETURN
      IF (ANY(SHAPE(point%product) /= [n, n]) .OR. ANY(SHAPE(point%factors) /= [n, n])) RETURN
      IF (ANY(SHAPE(point%rotations) /= [3, 3, n]) .OR. ANY(SHAPE(point%translations) /= [3, n])) RETURN
      IF (ANY(SHAPE(point%characters) /= [nu, ni]) .OR. ANY(point%multiplicity < 0)) RETURN
      IF (.NOT. ALL(ieee_is_finite(point%k)) .OR. .NOT. ALL(ieee_is_finite(point%cell))) RETURN
      IF (MAXVAL(ABS(point%cell)) < TINY(1.0_dp)) RETURN
      IF (.NOT. ALL(ieee_is_finite(point%translations))) RETURN
      IF (.NOT. ALL(ieee_is_finite(REAL(point%characters, dp)))) RETURN
      IF (.NOT. ALL(ieee_is_finite(AIMAG(point%characters)))) RETURN
      valid_point = .TRUE.
   END FUNCTION valid_point

! **************************************************************************************************
!> \brief Storage estimate in MiB for retained point data, excluding descriptors and library workspaces.
!> \param point Allocated point data
!> \return Estimated retained memory
! **************************************************************************************************
   REAL(dp) FUNCTION little_group_storage(point) RESULT(memory)
      TYPE(little_group_data_type), INTENT(IN)           :: point

      memory = 0.0_dp
      IF (.NOT. ALLOCATED(point%product)) RETURN
      memory = (20.0_dp*REAL(SIZE(point%product), dp) + 60.0_dp*SIZE(point%rotations, 3) + &
                16.0_dp*REAL(SIZE(point%characters), dp) + 4.0_dp*SIZE(point%multiplicity))/1024.0_dp**2
      IF (ALLOCATED(point%antiunitary)) memory = memory + 4.0_dp*SIZE(point%antiunitary)/1024.0_dp**2
   END FUNCTION little_group_storage
END MODULE topology_compatibility
