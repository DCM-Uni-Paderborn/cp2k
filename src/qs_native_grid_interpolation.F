!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
!> \brief Shared native-grid interpolation, descriptor transposes and geometry responses.
MODULE qs_native_grid_interpolation
   USE cell_types,                      ONLY: cell_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE pw_grid_types,                   ONLY: pw_grid_type
   USE qs_native_grid_cache,            ONLY: native_grid_interp_npts,&
                                              native_grid_interp_offset_max,&
                                              native_grid_interp_offset_min,&
                                              native_grid_interpolation_stencil_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   ! A wrapped stencil can touch both end tiles and one adjacent interior tile when the final
   ! tile is shorter than the stencil. Three tiles per direction are therefore sufficient.
   INTEGER, PARAMETER, PRIVATE :: native_grid_adjoint_tile_edge = 64, &
                                  native_grid_adjoint_max_tiles_per_direction = 3, &
                                  native_grid_adjoint_max_bins_per_row = &
                                  native_grid_adjoint_max_tiles_per_direction**3

   PUBLIC :: gather_native_grid_field, create_native_grid_interpolation_stencil
   PUBLIC :: interpolate_native_grid_fields, interpolate_native_grid, interpolate_native_grid_gradient
   PUBLIC :: build_native_grid_adjoint_bins, native_grid_adjoint_tile_bounds, add_native_grid_fields_adjoint_tile
   PUBLIC :: interpolate_native_grid_batch, native_grid_batch_adjoint, native_grid_batch_response

   INTEGER, PARAMETER :: native_grid_batch_npts = native_grid_interp_npts + 1

   TYPE native_grid_batch_stencil_type
      INTEGER :: relative_index(native_grid_batch_npts, 3) = 0, npoints = 0
      REAL(dp) :: weight(native_grid_batch_npts, 3) = 0.0_dp
      LOGICAL :: valid(native_grid_batch_npts, 3) = .FALSE., active = .FALSE.
   END TYPE native_grid_batch_stencil_type
CONTAINS

! **************************************************************************************************
!> \brief Replicate a distributed real-space field for atom-grid interpolation.
!> \param local_values ...
!> \param pw_grid ...
!> \param group ...
!> \param global_values ...
! **************************************************************************************************
   SUBROUTINE gather_native_grid_field(local_values, pw_grid, group, global_values)
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN), &
         POINTER                                         :: local_values
      TYPE(pw_grid_type), INTENT(IN), POINTER            :: pw_grid
      TYPE(mp_para_env_type), INTENT(IN)                 :: group
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(OUT), &
         POINTER                                         :: global_values

      INTEGER, DIMENSION(2, 3)                           :: bo

      CPASSERT(ASSOCIATED(local_values))
      CPASSERT(ASSOCIATED(pw_grid))
      CPASSERT(.NOT. ASSOCIATED(global_values))
      bo = pw_grid%bounds_local
      ALLOCATE (global_values(pw_grid%bounds(1, 1):pw_grid%bounds(2, 1), &
                              pw_grid%bounds(1, 2):pw_grid%bounds(2, 2), &
                              pw_grid%bounds(1, 3):pw_grid%bounds(2, 3)))
      global_values = 0.0_dp
      global_values(bo(1, 1):bo(2, 1), bo(1, 2):bo(2, 2), bo(1, 3):bo(2, 3)) = &
         local_values
      CALL group%sum(global_values)

   END SUBROUTINE gather_native_grid_field

! **************************************************************************************************
!> \brief Build the tensor-product interpolation stencil for one Cartesian point.
!> \param stencil ...
!> \param pw_grid ...
!> \param cell ...
!> \param point ...
!> \param wrap_auxiliary_cell wrap all auxiliary-grid directions
!> \param indices_only skip interpolation weights when only grid indices are needed
!> \param derivative_weights Optional fractional-coordinate derivatives of the same stencil
! **************************************************************************************************
   SUBROUTINE create_native_grid_interpolation_stencil(stencil, pw_grid, cell, point, &
                                                       wrap_auxiliary_cell, indices_only, derivative_weights)
      TYPE(native_grid_interpolation_stencil_type), &
         INTENT(OUT)                                     :: stencil
      TYPE(pw_grid_type), INTENT(IN), POINTER            :: pw_grid
      TYPE(cell_type), INTENT(IN), POINTER               :: cell
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: point
      LOGICAL, INTENT(IN), OPTIONAL                      :: wrap_auxiliary_cell, indices_only
      REAL(dp), INTENT(OUT), OPTIONAL :: derivative_weights(native_grid_interp_npts, 3)

      INTEGER                                            :: idir, inode, relative_index
      INTEGER, DIMENSION(3)                              :: base
      LOGICAL                                            :: build_weights, wrap_grid
      REAL(KIND=dp), DIMENSION(3)                        :: fraction, relative

      CPASSERT(ASSOCIATED(pw_grid))
      CPASSERT(ASSOCIATED(cell))
      stencil%active = .FALSE.
      stencil%valid = .FALSE.
      stencil%weight = 0.0_dp
      IF (PRESENT(derivative_weights)) derivative_weights = 0.0_dp
      build_weights = .TRUE.
      IF (PRESENT(indices_only)) build_weights = .NOT. indices_only
      CPASSERT(build_weights .OR. .NOT. PRESENT(derivative_weights))
      wrap_grid = ANY(cell%perd == 0)
      IF (PRESENT(wrap_auxiliary_cell)) wrap_grid = wrap_auxiliary_cell
      relative = MATMUL(pw_grid%dh_inv, point)
      DO idir = 1, 3
         IF (cell%perd(idir) == 1 .OR. wrap_grid) THEN
            relative(idir) = MODULO(relative(idir), REAL(pw_grid%npts(idir), KIND=dp))
         ELSE IF (relative(idir) <= -REAL(native_grid_interp_offset_max, dp) .OR. &
                  relative(idir) >= REAL(pw_grid%npts(idir) - &
                                         native_grid_interp_offset_min, KIND=dp)) THEN
            RETURN
         END IF
         base(idir) = FLOOR(relative(idir))
         FRACTION(idir) = relative(idir) - REAL(base(idir), KIND=dp)
         IF (build_weights) THEN
            IF (PRESENT(derivative_weights)) THEN
               CALL native_grid_lagrange_weights(FRACTION(idir), stencil%weight(:, idir), derivative_weights(:, idir))
            ELSE
               CALL native_grid_lagrange_weights(FRACTION(idir), stencil%weight(:, idir))
            END IF
         END IF
         DO inode = 1, native_grid_interp_npts
            relative_index = base(idir) + native_grid_interp_offset_min + inode - 1
            IF (cell%perd(idir) == 1 .OR. wrap_grid) THEN
               relative_index = MODULO(relative_index, pw_grid%npts(idir))
               stencil%valid(inode, idir) = .TRUE.
            ELSE IF (relative_index >= 0 .AND. relative_index < pw_grid%npts(idir)) THEN
               stencil%valid(inode, idir) = .TRUE.
            END IF
            stencil%relative_index(inode, idir) = relative_index
         END DO
      END DO
      stencil%active = .TRUE.

   END SUBROUTINE create_native_grid_interpolation_stencil

! **************************************************************************************************
!> \brief Interpolate density, gradient, and kinetic-density fields in one stencil traversal.
!> \param density ...
!> \param grad_x ...
!> \param grad_y ...
!> \param grad_z ...
!> \param kin ...
!> \param stencil ...
!> \param density_value ...
!> \param grad_value ...
!> \param kin_value ...
! **************************************************************************************************
   SUBROUTINE interpolate_native_grid_fields(density, grad_x, grad_y, grad_z, kin, stencil, &
                                             density_value, grad_value, kin_value)
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: density, grad_x, grad_y, grad_z, kin
      TYPE(native_grid_interpolation_stencil_type), &
         INTENT(IN)                                      :: stencil
      REAL(KIND=dp), INTENT(OUT)                         :: density_value
      REAL(KIND=dp), DIMENSION(3), INTENT(OUT)           :: grad_value
      REAL(KIND=dp), INTENT(OUT)                         :: kin_value

      INTEGER                                            :: inode_x, inode_y, inode_z
      INTEGER, DIMENSION(3)                              :: index, lower_bound
      REAL(KIND=dp)                                      :: coefficient

      CPASSERT(ALL(SHAPE(grad_x) == SHAPE(density)))
      CPASSERT(ALL(SHAPE(grad_y) == SHAPE(density)))
      CPASSERT(ALL(SHAPE(grad_z) == SHAPE(density)))
      CPASSERT(ALL(SHAPE(kin) == SHAPE(density)))
      density_value = 0.0_dp
      grad_value = 0.0_dp
      kin_value = 0.0_dp
      IF (.NOT. stencil%active) RETURN
      lower_bound = LBOUND(density)
      DO inode_z = 1, native_grid_interp_npts
         IF (.NOT. stencil%valid(inode_z, 3)) CYCLE
         INDEX(3) = lower_bound(3) + stencil%relative_index(inode_z, 3)
         DO inode_y = 1, native_grid_interp_npts
            IF (.NOT. stencil%valid(inode_y, 2)) CYCLE
            INDEX(2) = lower_bound(2) + stencil%relative_index(inode_y, 2)
            DO inode_x = 1, native_grid_interp_npts
               IF (.NOT. stencil%valid(inode_x, 1)) CYCLE
               INDEX(1) = lower_bound(1) + stencil%relative_index(inode_x, 1)
               coefficient = stencil%weight(inode_x, 1)* &
                             stencil%weight(inode_y, 2)* &
                             stencil%weight(inode_z, 3)
               density_value = density_value + coefficient*density(INDEX(1), INDEX(2), INDEX(3))
               grad_value(1) = grad_value(1) + coefficient*grad_x(INDEX(1), INDEX(2), INDEX(3))
               grad_value(2) = grad_value(2) + coefficient*grad_y(INDEX(1), INDEX(2), INDEX(3))
               grad_value(3) = grad_value(3) + coefficient*grad_z(INDEX(1), INDEX(2), INDEX(3))
               kin_value = kin_value + coefficient*kin(INDEX(1), INDEX(2), INDEX(3))
            END DO
         END DO
      END DO

   END SUBROUTINE interpolate_native_grid_fields

! **************************************************************************************************
!> \brief Group atom-grid rows by the disjoint PW tiles touched by their interpolation stencils.
!> \param pw_grid ...
!> \param cell ...
!> \param points atom-grid coordinates
!> \param wrap_auxiliary_cell wrap all auxiliary-grid directions
!> \param tile_count number of tiles along each PW-grid direction
!> \param bin_offsets CSR offsets into bin_rows
!> \param bin_rows atom-grid rows touching each tile
! **************************************************************************************************
   SUBROUTINE build_native_grid_adjoint_bins(pw_grid, cell, points, wrap_auxiliary_cell, &
                                             tile_count, bin_offsets, bin_rows)
      TYPE(pw_grid_type), INTENT(IN), POINTER            :: pw_grid
      TYPE(cell_type), INTENT(IN), POINTER               :: cell
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: points
      LOGICAL, INTENT(IN)                                :: wrap_auxiliary_cell
      INTEGER, DIMENSION(3), INTENT(OUT)                 :: tile_count
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: bin_offsets, bin_rows

      INTEGER                                            :: candidate, ibin, idir, inode, irow, ix, &
                                                            iy, iz, nbin, ntouched, position
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: bin_counts, bin_cursor, row_bin_count
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: row_bins
      INTEGER, DIMENSION(3)                              :: touched_count
      INTEGER, DIMENSION(&
         native_grid_adjoint_max_tiles_per_direction, 3) :: touched_tiles
      TYPE(native_grid_interpolation_stencil_type)       :: stencil

      CPASSERT(ASSOCIATED(pw_grid))
      CPASSERT(ASSOCIATED(cell))
      CPASSERT(SIZE(points, 1) == 3)
      tile_count = (pw_grid%npts + native_grid_adjoint_tile_edge - 1)/ &
                   native_grid_adjoint_tile_edge
      nbin = PRODUCT(tile_count)
      ALLOCATE (row_bin_count(SIZE(points, 2)), &
                row_bins(native_grid_adjoint_max_bins_per_row, SIZE(points, 2)))
      row_bin_count = 0
      row_bins = 0
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
!$OMP PRIVATE(candidate, idir, inode, ix, iy, iz, ntouched, stencil, &
!$OMP         touched_count, touched_tiles) &
!$OMP SHARED(cell, points, pw_grid, row_bin_count, row_bins, tile_count, wrap_auxiliary_cell)
      DO irow = 1, SIZE(points, 2)
         CALL create_native_grid_interpolation_stencil( &
            stencil, pw_grid, cell, points(:, irow), wrap_auxiliary_cell, indices_only=.TRUE.)
         IF (.NOT. stencil%active) CYCLE
         touched_count = 0
         touched_tiles = 0
         DO idir = 1, 3
            DO inode = 1, native_grid_interp_npts
               IF (.NOT. stencil%valid(inode, idir)) CYCLE
               candidate = stencil%relative_index(inode, idir)/native_grid_adjoint_tile_edge + 1
               IF (ALL(touched_tiles(1:touched_count(idir), idir) /= candidate)) THEN
                  touched_count(idir) = touched_count(idir) + 1
                  CPASSERT(touched_count(idir) <= SIZE(touched_tiles, 1))
                  touched_tiles(touched_count(idir), idir) = candidate
               END IF
            END DO
         END DO
         ntouched = 0
         DO iz = 1, touched_count(3)
            DO iy = 1, touched_count(2)
               DO ix = 1, touched_count(1)
                  ntouched = ntouched + 1
                  CPASSERT(ntouched <= native_grid_adjoint_max_bins_per_row)
                  row_bins(ntouched, irow) = 1 + touched_tiles(ix, 1) - 1 + tile_count(1)*( &
                                             touched_tiles(iy, 2) - 1 + tile_count(2)*( &
                                             touched_tiles(iz, 3) - 1))
               END DO
            END DO
         END DO
         row_bin_count(irow) = ntouched
      END DO
!$OMP END PARALLEL DO

      ALLOCATE (bin_counts(nbin), bin_cursor(nbin), bin_offsets(nbin + 1))
      bin_counts = 0
      DO irow = 1, SIZE(points, 2)
         DO ibin = 1, row_bin_count(irow)
            bin_counts(row_bins(ibin, irow)) = bin_counts(row_bins(ibin, irow)) + 1
         END DO
      END DO
      bin_offsets(1) = 1
      DO ibin = 1, nbin
         bin_offsets(ibin + 1) = bin_offsets(ibin) + bin_counts(ibin)
      END DO
      ALLOCATE (bin_rows(bin_offsets(nbin + 1) - 1))
      bin_cursor(:) = bin_offsets(1:nbin)
      DO irow = 1, SIZE(points, 2)
         DO ibin = 1, row_bin_count(irow)
            candidate = row_bins(ibin, irow)
            position = bin_cursor(candidate)
            bin_rows(position) = irow
            bin_cursor(candidate) = position + 1
         END DO
      END DO
      DEALLOCATE (bin_counts, bin_cursor, row_bin_count, row_bins)

   END SUBROUTINE build_native_grid_adjoint_bins

! **************************************************************************************************
!> \brief Return the inclusive PW-grid bounds owned by one linear tile index.
!> \param pw_grid ...
!> \param tile_index linear tile index
!> \param tile_count number of tiles along each PW-grid direction
!> \param tile_lower zero-based lower grid index
!> \param tile_upper zero-based upper grid index
! **************************************************************************************************
   SUBROUTINE native_grid_adjoint_tile_bounds(pw_grid, tile_index, tile_count, &
                                              tile_lower, tile_upper)
      TYPE(pw_grid_type), INTENT(IN), POINTER            :: pw_grid
      INTEGER, INTENT(IN)                                :: tile_index
      INTEGER, DIMENSION(3), INTENT(IN)                  :: tile_count
      INTEGER, DIMENSION(3), INTENT(OUT)                 :: tile_lower, tile_upper

      INTEGER                                            :: linear_tile
      INTEGER, DIMENSION(3)                              :: expected_tile_count, tile_coord

      CPASSERT(ASSOCIATED(pw_grid))
      expected_tile_count = (pw_grid%npts + native_grid_adjoint_tile_edge - 1)/ &
                            native_grid_adjoint_tile_edge
      CPASSERT(ALL(tile_count == expected_tile_count))
      CPASSERT(tile_index >= 1 .AND. tile_index <= PRODUCT(tile_count))
      linear_tile = tile_index - 1
      tile_coord(1) = MOD(linear_tile, tile_count(1))
      linear_tile = linear_tile/tile_count(1)
      tile_coord(2) = MOD(linear_tile, tile_count(2))
      tile_coord(3) = linear_tile/tile_count(2)
      tile_lower = tile_coord*native_grid_adjoint_tile_edge
      tile_upper = MIN(tile_lower + native_grid_adjoint_tile_edge - 1, pw_grid%npts - 1)

   END SUBROUTINE native_grid_adjoint_tile_bounds

! **************************************************************************************************
!> \brief Apply the primitive-field interpolation transpose inside one disjoint PW tile.
!> \param density ...
!> \param grad ...
!> \param kin ...
!> \param pw_grid ...
!> \param stencil ...
!> \param density_value ...
!> \param grad_value ...
!> \param kin_value ...
!> \param nchannels number of spin channels accumulated by the adjoint
!> \param tile_lower zero-based lower grid index owned by this call
!> \param tile_upper zero-based upper grid index owned by this call
! **************************************************************************************************
   SUBROUTINE add_native_grid_fields_adjoint_tile( &
      density, grad, kin, pw_grid, stencil, density_value, grad_value, kin_value, nchannels, &
      tile_lower, tile_upper)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: density
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: grad
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: kin
      TYPE(pw_grid_type), INTENT(IN), POINTER            :: pw_grid
      TYPE(native_grid_interpolation_stencil_type), &
         INTENT(IN)                                      :: stencil
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: density_value
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: grad_value
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: kin_value
      INTEGER, INTENT(IN)                                :: nchannels
      INTEGER, DIMENSION(3), INTENT(IN)                  :: tile_lower, tile_upper

      INTEGER                                            :: idir, inode_x, inode_y, inode_z, ipt, &
                                                            ispin
      INTEGER, DIMENSION(3)                              :: relative_index
      REAL(KIND=dp)                                      :: coefficient

      CPASSERT(ASSOCIATED(pw_grid))
      CPASSERT(SIZE(density, 1) == PRODUCT(pw_grid%npts))
      CPASSERT(SIZE(density, 2) == nchannels)
      CPASSERT(SIZE(grad, 1) == PRODUCT(pw_grid%npts))
      CPASSERT(SIZE(grad, 2) == 3)
      CPASSERT(SIZE(grad, 3) == nchannels)
      CPASSERT(SIZE(kin, 1) == PRODUCT(pw_grid%npts))
      CPASSERT(SIZE(kin, 2) == nchannels)
      CPASSERT(SIZE(density_value) == 2)
      CPASSERT(ALL(SHAPE(grad_value) == [3, 2]))
      CPASSERT(SIZE(kin_value) == 2)
      CPASSERT(nchannels >= 1 .AND. nchannels <= 2)
      CPASSERT(ALL(tile_lower >= 0))
      CPASSERT(ALL(tile_upper >= tile_lower))
      CPASSERT(ALL(tile_upper < pw_grid%npts))
      IF (.NOT. stencil%active) RETURN
      DO inode_z = 1, native_grid_interp_npts
         IF (.NOT. stencil%valid(inode_z, 3)) CYCLE
         relative_index(3) = stencil%relative_index(inode_z, 3)
         IF (relative_index(3) < tile_lower(3) .OR. relative_index(3) > tile_upper(3)) CYCLE
         DO inode_y = 1, native_grid_interp_npts
            IF (.NOT. stencil%valid(inode_y, 2)) CYCLE
            relative_index(2) = stencil%relative_index(inode_y, 2)
            IF (relative_index(2) < tile_lower(2) .OR. relative_index(2) > tile_upper(2)) CYCLE
            DO inode_x = 1, native_grid_interp_npts
               IF (.NOT. stencil%valid(inode_x, 1)) CYCLE
               relative_index(1) = stencil%relative_index(inode_x, 1)
               IF (relative_index(1) < tile_lower(1) .OR. relative_index(1) > tile_upper(1)) CYCLE
               ipt = 1 + relative_index(1) + pw_grid%npts(1)*( &
                     relative_index(2) + pw_grid%npts(2)*relative_index(3))
               coefficient = stencil%weight(inode_x, 1)* &
                             stencil%weight(inode_y, 2)* &
                             stencil%weight(inode_z, 3)
               DO ispin = 1, nchannels
                  density(ipt, ispin) = density(ipt, ispin) + coefficient*density_value(ispin)
                  DO idir = 1, 3
                     grad(ipt, idir, ispin) = grad(ipt, idir, ispin) + &
                                              coefficient*grad_value(idir, ispin)
                  END DO
                  kin(ipt, ispin) = kin(ipt, ispin) + coefficient*kin_value(ispin)
               END DO
            END DO
         END DO
      END DO

   END SUBROUTINE add_native_grid_fields_adjoint_tile

! **************************************************************************************************
!> \brief Interpolate a replicated native-grid field at a Cartesian point.
!> \param values ...
!> \param pw_grid ...
!> \param cell ...
!> \param point ...
!> \param wrap_auxiliary_cell wrap all auxiliary-grid directions
!> \return ...
! **************************************************************************************************
   FUNCTION interpolate_native_grid(values, pw_grid, cell, point, wrap_auxiliary_cell) RESULT(value)
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: values
      TYPE(pw_grid_type), INTENT(IN), POINTER            :: pw_grid
      TYPE(cell_type), INTENT(IN), POINTER               :: cell
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: point
      LOGICAL, INTENT(IN), OPTIONAL                      :: wrap_auxiliary_cell
      REAL(KIND=dp)                                      :: value

      INTEGER                                            :: corner_x, corner_y, corner_z, idir
      INTEGER, DIMENSION(3)                              :: base, index, relative_index
      LOGICAL                                            :: wrap_grid
      REAL(KIND=dp)                                      :: coefficient
      REAL(KIND=dp), DIMENSION(3)                        :: fraction, relative
      REAL(KIND=dp), &
         DIMENSION(native_grid_interp_npts, 3)           :: weights

      CPASSERT(ASSOCIATED(pw_grid))
      CPASSERT(ASSOCIATED(cell))
      DO idir = 1, 3
         CPASSERT(SIZE(values, idir) == pw_grid%npts(idir))
      END DO

      wrap_grid = ANY(cell%perd == 0)
      IF (PRESENT(wrap_auxiliary_cell)) wrap_grid = wrap_auxiliary_cell
      relative = MATMUL(pw_grid%dh_inv, point)
      DO idir = 1, 3
         IF (cell%perd(idir) == 1 .OR. wrap_grid) THEN
            relative(idir) = MODULO(relative(idir), REAL(pw_grid%npts(idir), KIND=dp))
         ELSE IF (relative(idir) <= -REAL(native_grid_interp_offset_max, dp) .OR. &
                  relative(idir) >= REAL(pw_grid%npts(idir) - &
                                         native_grid_interp_offset_min, KIND=dp)) THEN
            value = 0.0_dp
            RETURN
         END IF
         base(idir) = FLOOR(relative(idir))
         FRACTION(idir) = relative(idir) - REAL(base(idir), KIND=dp)
         CALL native_grid_lagrange_weights(FRACTION(idir), weights(:, idir))
      END DO

      value = 0.0_dp
      DO corner_z = native_grid_interp_offset_min, native_grid_interp_offset_max
         DO corner_y = native_grid_interp_offset_min, native_grid_interp_offset_max
            DO corner_x = native_grid_interp_offset_min, native_grid_interp_offset_max
               relative_index = base + [corner_x, corner_y, corner_z]
               coefficient = weights(corner_x - native_grid_interp_offset_min + 1, 1)* &
                             weights(corner_y - native_grid_interp_offset_min + 1, 2)* &
                             weights(corner_z - native_grid_interp_offset_min + 1, 3)
               DO idir = 1, 3
                  IF (cell%perd(idir) == 1 .OR. wrap_grid) THEN
                     relative_index(idir) = MODULO(relative_index(idir), pw_grid%npts(idir))
                  ELSE IF (relative_index(idir) < 0 .OR. &
                           relative_index(idir) >= pw_grid%npts(idir)) THEN
                     coefficient = 0.0_dp
                  END IF
               END DO
               IF (coefficient == 0.0_dp) CYCLE
               index = LBOUND(values) + relative_index
               value = value + coefficient*values(INDEX(1), INDEX(2), INDEX(3))
            END DO
         END DO
      END DO

   END FUNCTION interpolate_native_grid

! **************************************************************************************************
!> \brief Return the Cartesian gradient of native-grid interpolation at one point.
!> \param values ...
!> \param pw_grid ...
!> \param cell ...
!> \param point ...
!> \param wrap_auxiliary_cell wrap all auxiliary-grid directions
!> \return ...
! **************************************************************************************************
   FUNCTION interpolate_native_grid_gradient(values, pw_grid, cell, point, &
                                             wrap_auxiliary_cell) RESULT(gradient)
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: values
      TYPE(pw_grid_type), INTENT(IN), POINTER            :: pw_grid
      TYPE(cell_type), INTENT(IN), POINTER               :: cell
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: point
      LOGICAL, INTENT(IN), OPTIONAL                      :: wrap_auxiliary_cell
      REAL(KIND=dp), DIMENSION(3)                        :: gradient

      INTEGER                                            :: corner_x, corner_y, corner_z, idir, jdir
      INTEGER, DIMENSION(3)                              :: base, corner_index, index, relative_index
      LOGICAL                                            :: wrap_grid
      REAL(KIND=dp)                                      :: coefficient
      REAL(KIND=dp), DIMENSION(3)                        :: fraction, gradient_relative, relative
      REAL(KIND=dp), &
         DIMENSION(native_grid_interp_npts, 3)           :: derivative_weights, weights

      CPASSERT(ASSOCIATED(pw_grid))
      CPASSERT(ASSOCIATED(cell))
      DO idir = 1, 3
         CPASSERT(SIZE(values, idir) == pw_grid%npts(idir))
      END DO

      wrap_grid = ANY(cell%perd == 0)
      IF (PRESENT(wrap_auxiliary_cell)) wrap_grid = wrap_auxiliary_cell
      relative = MATMUL(pw_grid%dh_inv, point)
      DO idir = 1, 3
         IF (cell%perd(idir) == 1 .OR. wrap_grid) THEN
            relative(idir) = MODULO(relative(idir), REAL(pw_grid%npts(idir), KIND=dp))
         ELSE IF (relative(idir) <= -REAL(native_grid_interp_offset_max, dp) .OR. &
                  relative(idir) >= REAL(pw_grid%npts(idir) - &
                                         native_grid_interp_offset_min, KIND=dp)) THEN
            gradient = 0.0_dp
            RETURN
         END IF
         base(idir) = FLOOR(relative(idir))
         FRACTION(idir) = relative(idir) - REAL(base(idir), KIND=dp)
         CALL native_grid_lagrange_weights( &
            FRACTION(idir), weights(:, idir), derivative_weights(:, idir))
      END DO

      gradient_relative = 0.0_dp
      DO corner_z = native_grid_interp_offset_min, native_grid_interp_offset_max
         DO corner_y = native_grid_interp_offset_min, native_grid_interp_offset_max
            DO corner_x = native_grid_interp_offset_min, native_grid_interp_offset_max
               relative_index = base + [corner_x, corner_y, corner_z]
               corner_index = [corner_x, corner_y, corner_z] - &
                              native_grid_interp_offset_min + 1
               DO idir = 1, 3
                  IF (cell%perd(idir) == 1 .OR. wrap_grid) THEN
                     relative_index(idir) = MODULO(relative_index(idir), pw_grid%npts(idir))
                  ELSE IF (relative_index(idir) < 0 .OR. &
                           relative_index(idir) >= pw_grid%npts(idir)) THEN
                     EXIT
                  END IF
               END DO
               IF (idir <= 3) CYCLE
               index = LBOUND(values) + relative_index
               DO idir = 1, 3
                  coefficient = 1.0_dp
                  DO jdir = 1, 3
                     IF (jdir == idir) THEN
                        coefficient = coefficient*derivative_weights( &
                                      corner_index(jdir), jdir)
                     ELSE
                        coefficient = coefficient*weights( &
                                      corner_index(jdir), jdir)
                     END IF
                  END DO
                  gradient_relative(idir) = gradient_relative(idir) + &
                                            coefficient*values(INDEX(1), INDEX(2), INDEX(3))
               END DO
            END DO
         END DO
      END DO
      gradient = MATMUL(TRANSPOSE(pw_grid%dh_inv), gradient_relative)

   END FUNCTION interpolate_native_grid_gradient

! **************************************************************************************************
!> \brief Build the value and optional derivative weights for the native-grid interpolation.
!> \param fraction fractional coordinate between two grid points
!> \param weights interpolation weights
!> \param derivative_weights optional derivatives with respect to fraction
! **************************************************************************************************
   SUBROUTINE native_grid_lagrange_weights(fraction, weights, derivative_weights)
      REAL(KIND=dp), INTENT(IN)                          :: fraction
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: weights
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT), OPTIONAL :: derivative_weights

      INTEGER                                            :: inode, jnode, knode, npoints
      REAL(KIND=dp)                                      :: denominator, derivative, numerator
      REAL(KIND=dp), DIMENSION(SIZE(weights))            :: nodes

      npoints = SIZE(weights)
      IF (PRESENT(derivative_weights)) THEN
         CPASSERT(SIZE(derivative_weights) == npoints)
      END IF
      DO inode = 1, npoints
         nodes(inode) = REAL(native_grid_interp_offset_min + inode - 1, dp)
      END DO
      DO inode = 1, npoints
         denominator = 1.0_dp
         numerator = 1.0_dp
         DO jnode = 1, npoints
            IF (jnode == inode) CYCLE
            denominator = denominator*(nodes(inode) - nodes(jnode))
            numerator = numerator*(fraction - nodes(jnode))
         END DO
         weights(inode) = numerator/denominator
         IF (PRESENT(derivative_weights)) THEN
            derivative = 0.0_dp
            DO knode = 1, npoints
               IF (knode == inode) CYCLE
               numerator = 1.0_dp
               DO jnode = 1, npoints
                  IF (jnode == inode .OR. jnode == knode) CYCLE
                  numerator = numerator*(fraction - nodes(jnode))
               END DO
               derivative = derivative + numerator/denominator
            END DO
            derivative_weights(inode) = derivative
         END IF
      END DO

   END SUBROUTINE native_grid_lagrange_weights

! **************************************************************************************************
!> \brief Select the native stencil or its C2 blend without changing existing native-grid consumers.
!> \param stencil ...
!> \param pw_grid ...
!> \param cell ...
!> \param point ...
!> \param wrap_auxiliary_cell ...
!> \param continuous ...
!> \param derivative_weights ...
! **************************************************************************************************
   SUBROUTINE create_native_grid_batch_stencil(stencil, pw_grid, cell, point, wrap_auxiliary_cell, &
                                               continuous, derivative_weights)
      TYPE(native_grid_batch_stencil_type), INTENT(OUT)  :: stencil
      TYPE(pw_grid_type), INTENT(IN), POINTER            :: pw_grid
      TYPE(cell_type), INTENT(IN), POINTER               :: cell
      REAL(dp), INTENT(IN)                               :: point(3)
      LOGICAL, INTENT(IN), OPTIONAL                      :: wrap_auxiliary_cell, continuous
      REAL(dp), INTENT(OUT), OPTIONAL :: derivative_weights(native_grid_batch_npts, 3)

      INTEGER                                            :: axis, base, inode, relative_index
      LOGICAL                                            :: smooth, wrap
      REAL(dp) :: blend, blend_derivative, fraction, left(native_grid_interp_npts), &
         left_derivative(native_grid_interp_npts), native_derivative(native_grid_interp_npts, 3), &
         relative(3), right(native_grid_interp_npts), right_derivative(native_grid_interp_npts), t
      TYPE(native_grid_interpolation_stencil_type)       :: native

      smooth = .FALSE.
      IF (PRESENT(continuous)) smooth = continuous
      IF (PRESENT(derivative_weights)) derivative_weights = 0.0_dp
      IF (.NOT. smooth) THEN
         IF (PRESENT(derivative_weights)) THEN
            CALL create_native_grid_interpolation_stencil(native, pw_grid, cell, point, &
                                                          wrap_auxiliary_cell, derivative_weights=native_derivative)
            derivative_weights(:native_grid_interp_npts, :) = native_derivative
         ELSE
            CALL create_native_grid_interpolation_stencil(native, pw_grid, cell, point, wrap_auxiliary_cell)
         END IF
         stencil%active = native%active
         stencil%npoints = native_grid_interp_npts
         stencil%relative_index(:native_grid_interp_npts, :) = native%relative_index
         stencil%valid(:native_grid_interp_npts, :) = native%valid
         stencil%weight(:native_grid_interp_npts, :) = native%weight
         RETURN
      END IF
      CPASSERT(ASSOCIATED(pw_grid))
      CPASSERT(ASSOCIATED(cell))
      wrap = ANY(cell%perd == 0)
      IF (PRESENT(wrap_auxiliary_cell)) wrap = wrap_auxiliary_cell
      relative = MATMUL(pw_grid%dh_inv, point)
      stencil%npoints = native_grid_batch_npts
      DO axis = 1, 3
         IF (cell%perd(axis) == 1 .OR. wrap) THEN
            relative(axis) = MODULO(relative(axis), REAL(pw_grid%npts(axis), dp))
         ELSE IF (relative(axis) <= -REAL(native_grid_interp_offset_max, dp) - 0.5_dp .OR. &
                  relative(axis) >= REAL(pw_grid%npts(axis) - native_grid_interp_offset_min, dp) + 0.5_dp) THEN
            RETURN
         END IF
         base = FLOOR(relative(axis) + 0.5_dp)
         fraction = relative(axis) - REAL(base, dp)
         ! Translating the first polynomial by one node reuses the same native Lagrange routine.
         IF (PRESENT(derivative_weights)) THEN
            CALL native_grid_lagrange_weights(fraction + 1.0_dp, left, left_derivative)
            CALL native_grid_lagrange_weights(fraction, right, right_derivative)
         ELSE
            CALL native_grid_lagrange_weights(fraction + 1.0_dp, left)
            CALL native_grid_lagrange_weights(fraction, right)
         END IF
         t = fraction + 0.5_dp
         blend = t**3*(10.0_dp + t*(-15.0_dp + 6.0_dp*t))
         stencil%weight(:native_grid_interp_npts, axis) = (1.0_dp - blend)*left
         stencil%weight(2:, axis) = stencil%weight(2:, axis) + blend*right
         IF (PRESENT(derivative_weights)) THEN
            blend_derivative = 30.0_dp*t**2*(1.0_dp - t)**2
            derivative_weights(:native_grid_interp_npts, axis) = &
               (1.0_dp - blend)*left_derivative - blend_derivative*left
            derivative_weights(2:, axis) = derivative_weights(2:, axis) + &
                                           blend*right_derivative + blend_derivative*right
         END IF
         DO inode = 1, native_grid_batch_npts
            relative_index = base + native_grid_interp_offset_min - 1 + inode - 1
            IF (cell%perd(axis) == 1 .OR. wrap) THEN
               relative_index = MODULO(relative_index, pw_grid%npts(axis))
               stencil%valid(inode, axis) = .TRUE.
            ELSE IF (relative_index >= 0 .AND. relative_index < pw_grid%npts(axis)) THEN
               stencil%valid(inode, axis) = .TRUE.
            END IF
            stencil%relative_index(inode, axis) = relative_index
         END DO
      END DO
      stencil%active = .TRUE.
   END SUBROUTINE create_native_grid_batch_stencil

! **************************************************************************************************
!> \brief Interpolate replicated native-grid descriptors to a local point batch.
!> \param fields Full grid fields, indexed (x,y,z,descriptor), with no quadrature weights.
!> \param pw_grid Native grid whose ordering is retained.
!> \param cell Periodicity used by the stencil.
!> \param points Local Cartesian points, indexed (axis,point), possibly empty.
!> \param values Point descriptors, indexed (point,descriptor).
!> \param wrap_auxiliary_cell Wrap all auxiliary-grid directions.
!> \param continuous Use the C2 blend for moving atom-grid points. The default retains the native stencil.
! **************************************************************************************************
   SUBROUTINE interpolate_native_grid_batch(fields, pw_grid, cell, points, values, wrap_auxiliary_cell, continuous)
      REAL(dp), INTENT(IN)                               :: fields(:, :, :, :)
      TYPE(pw_grid_type), INTENT(IN), POINTER            :: pw_grid
      TYPE(cell_type), INTENT(IN), POINTER               :: cell
      REAL(dp), INTENT(IN)                               :: points(:, :)
      REAL(dp), INTENT(OUT)                              :: values(:, :)
      LOGICAL, INTENT(IN), OPTIONAL                      :: wrap_auxiliary_cell, continuous

      INTEGER                                            :: index(3), ipoint, ix, iy, iz
      REAL(dp)                                           :: coefficient
      TYPE(native_grid_batch_stencil_type)               :: stencil

      CPASSERT(ASSOCIATED(pw_grid))
      CPASSERT(SIZE(fields, 1) == pw_grid%npts(1))
      CPASSERT(SIZE(fields, 2) == pw_grid%npts(2))
      CPASSERT(SIZE(fields, 3) == pw_grid%npts(3))
      CPASSERT(SIZE(points, 1) == 3)
      CPASSERT(ALL(SHAPE(values) == [SIZE(points, 2), SIZE(fields, 4)]))
      values = 0.0_dp
      DO ipoint = 1, SIZE(points, 2)
         CALL create_native_grid_batch_stencil(stencil, pw_grid, cell, points(:, ipoint), wrap_auxiliary_cell, continuous)
         IF (.NOT. stencil%active) CYCLE
         DO iz = 1, stencil%npoints
            IF (.NOT. stencil%valid(iz, 3)) CYCLE
            INDEX(3) = 1 + stencil%relative_index(iz, 3)
            DO iy = 1, stencil%npoints
               IF (.NOT. stencil%valid(iy, 2)) CYCLE
               INDEX(2) = 1 + stencil%relative_index(iy, 2)
               DO ix = 1, stencil%npoints
                  IF (.NOT. stencil%valid(ix, 1)) CYCLE
                  INDEX(1) = 1 + stencil%relative_index(ix, 1)
                  coefficient = stencil%weight(ix, 1)*stencil%weight(iy, 2)*stencil%weight(iz, 3)
                  values(ipoint, :) = values(ipoint, :) + coefficient*fields(INDEX(1), INDEX(2), INDEX(3), :)
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE interpolate_native_grid_batch

! **************************************************************************************************
!> \brief Accumulate the exact descriptor transpose of the same interpolation.
!> \param values Point covectors, indexed (point,descriptor), including any atom-grid measure.
!> \param pw_grid ...
!> \param cell ...
!> \param points ...
!> \param fields Full-grid Euclidean covectors, incremented. This routine has no MPI reduction.
!>               Reduce distributed point batches before restricting them to native grid owners.
!>               Divide by the PW volume element before treating them as unweighted potentials.
!> \param wrap_auxiliary_cell ...
!> \param continuous ...
! **************************************************************************************************
   SUBROUTINE native_grid_batch_adjoint(values, pw_grid, cell, points, fields, wrap_auxiliary_cell, continuous)
      REAL(dp), INTENT(IN)                               :: values(:, :)
      TYPE(pw_grid_type), INTENT(IN), POINTER            :: pw_grid
      TYPE(cell_type), INTENT(IN), POINTER               :: cell
      REAL(dp), INTENT(IN)                               :: points(:, :)
      REAL(dp), INTENT(INOUT)                            :: fields(:, :, :, :)
      LOGICAL, INTENT(IN), OPTIONAL                      :: wrap_auxiliary_cell, continuous

      INTEGER                                            :: index(3), ipoint, ix, iy, iz
      REAL(dp)                                           :: coefficient
      TYPE(native_grid_batch_stencil_type)               :: stencil

      CPASSERT(ASSOCIATED(pw_grid))
      CPASSERT(SIZE(fields, 1) == pw_grid%npts(1))
      CPASSERT(SIZE(fields, 2) == pw_grid%npts(2))
      CPASSERT(SIZE(fields, 3) == pw_grid%npts(3))
      CPASSERT(SIZE(points, 1) == 3)
      CPASSERT(ALL(SHAPE(values) == [SIZE(points, 2), SIZE(fields, 4)]))
      DO ipoint = 1, SIZE(points, 2)
         CALL create_native_grid_batch_stencil(stencil, pw_grid, cell, points(:, ipoint), wrap_auxiliary_cell, continuous)
         IF (.NOT. stencil%active) CYCLE
         DO iz = 1, stencil%npoints
            IF (.NOT. stencil%valid(iz, 3)) CYCLE
            INDEX(3) = 1 + stencil%relative_index(iz, 3)
            DO iy = 1, stencil%npoints
               IF (.NOT. stencil%valid(iy, 2)) CYCLE
               INDEX(2) = 1 + stencil%relative_index(iy, 2)
               DO ix = 1, stencil%npoints
                  IF (.NOT. stencil%valid(ix, 1)) CYCLE
                  INDEX(1) = 1 + stencil%relative_index(ix, 1)
                  coefficient = stencil%weight(ix, 1)*stencil%weight(iy, 2)*stencil%weight(iz, 3)
                  fields(INDEX(1), INDEX(2), INDEX(3), :) = fields(INDEX(1), INDEX(2), INDEX(3), :) + &
                                                            coefficient*values(ipoint, :)
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE native_grid_batch_adjoint

! **************************************************************************************************
!> \brief Joint source-field, point and PW-lattice response of descriptor interpolation.
!> \param fields ...
!> \param fields_response Source descriptor response at fixed native grid indices.
!> \param pw_grid ...
!> \param cell ...
!> \param points ...
!> \param point_response Motion of the Cartesian query points.
!> \param lattice_response Response of the PW lattice, dh(:,axis)*npts(axis), at fixed grid size.
!> \param response Point descriptor response, including both I*delta(fields) and delta(I)*fields.
!> \param wrap_auxiliary_cell ...
!> \param continuous ...
! **************************************************************************************************
   SUBROUTINE native_grid_batch_response(fields, fields_response, pw_grid, cell, points, &
                                         point_response, lattice_response, response, wrap_auxiliary_cell, continuous)
      REAL(dp), INTENT(IN)                               :: fields(:, :, :, :), &
                                                            fields_response(:, :, :, :)
      TYPE(pw_grid_type), INTENT(IN), POINTER            :: pw_grid
      TYPE(cell_type), INTENT(IN), POINTER               :: cell
      REAL(dp), INTENT(IN)                               :: points(:, :), point_response(:, :), &
                                                            lattice_response(3, 3)
      REAL(dp), INTENT(OUT)                              :: response(:, :)
      LOGICAL, INTENT(IN), OPTIONAL                      :: wrap_auxiliary_cell, continuous

      INTEGER                                            :: axis, INDEX(3), ipoint, ix, iy, iz, &
                                                            nodes(3), other
      REAL(dp) :: coefficient, derivative, derivative_weights(native_grid_batch_npts, 3), &
         relative_response(3), spacing_response(3, 3)
      TYPE(native_grid_batch_stencil_type)               :: stencil

      CPASSERT(ALL(SHAPE(fields_response) == SHAPE(fields)))
      CPASSERT(ALL(SHAPE(point_response) == SHAPE(points)))
      CALL interpolate_native_grid_batch(fields_response, pw_grid, cell, points, response, wrap_auxiliary_cell, continuous)
      DO axis = 1, 3
         spacing_response(:, axis) = lattice_response(:, axis)/REAL(pw_grid%npts(axis), dp)
      END DO
      DO ipoint = 1, SIZE(points, 2)
         CALL create_native_grid_batch_stencil(stencil, pw_grid, cell, points(:, ipoint), &
                                               wrap_auxiliary_cell, continuous, derivative_weights)
         IF (.NOT. stencil%active) CYCLE
         relative_response = MATMUL(pw_grid%dh_inv, point_response(:, ipoint) - &
                                    MATMUL(spacing_response, MATMUL(pw_grid%dh_inv, points(:, ipoint))))
         DO iz = 1, stencil%npoints
            IF (.NOT. stencil%valid(iz, 3)) CYCLE
            INDEX(3) = 1 + stencil%relative_index(iz, 3)
            DO iy = 1, stencil%npoints
               IF (.NOT. stencil%valid(iy, 2)) CYCLE
               INDEX(2) = 1 + stencil%relative_index(iy, 2)
               DO ix = 1, stencil%npoints
                  IF (.NOT. stencil%valid(ix, 1)) CYCLE
                  INDEX(1) = 1 + stencil%relative_index(ix, 1)
                  nodes = [ix, iy, iz]
                  coefficient = 0.0_dp
                  DO axis = 1, 3
                     derivative = derivative_weights(nodes(axis), axis)*relative_response(axis)
                     DO other = 1, 3
                        IF (other /= axis) derivative = derivative*stencil%weight(nodes(other), other)
                     END DO
                     coefficient = coefficient + derivative
                  END DO
                  response(ipoint, :) = response(ipoint, :) + coefficient*fields(INDEX(1), INDEX(2), INDEX(3), :)
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE native_grid_batch_response

END MODULE qs_native_grid_interpolation
