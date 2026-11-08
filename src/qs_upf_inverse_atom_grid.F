!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
!> \brief Atom-centered quadrature and its geometry response for inverse mixed XC.
MODULE qs_upf_inverse_atom_grid
   USE cell_types,                      ONLY: cell_type
   USE cryssym,                         ONLY: crys_sym_lattice_rotations
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE message_passing,                 ONLY: mp_comm_type
   USE qs_grid_atom,                    ONLY: grid_atom_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE skala_periodic_partition,        ONLY: build_periodic_atom_image_layout,&
                                              periodic_atom_image_partition_from_layout,&
                                              skala_periodic_partition_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   INTEGER(int_8), PARAMETER :: partition_cache_bytes = 64_int_8*1024_int_8**2

   TYPE qs_upf_inverse_atom_grid_ref_type
      TYPE(grid_atom_type), POINTER :: grid => NULL()
      REAL(dp), ALLOCATABLE :: directions(:, :), angular_weights(:)
   END TYPE qs_upf_inverse_atom_grid_ref_type

   TYPE qs_upf_inverse_atom_grid_type
      INTEGER, ALLOCATABLE :: atoms(:)
      ! Row within the complete atom quadrature, independent of MPI distribution.
      INTEGER, ALLOCATABLE :: atom_rows(:)
      REAL(dp), ALLOCATABLE :: points(:, :), weights(:), base_weights(:)
      ! The FFT background is integrated once per auxiliary cell also for isolated electrostatics.
      REAL(dp), ALLOCATABLE :: background_weights(:)
      REAL(dp) :: inverse_cell(3, 3) = 0.0_dp
      INTEGER :: background_axes(3) = 0, natom = 0
      TYPE(skala_periodic_partition_type) :: partition
      ! Atom derivatives followed by the three strain columns, for each local row.
      REAL(dp), ALLOCATABLE :: partition_derivatives(:, :, :)
      ! Native inverse frames and this quadrature share one geometry lifetime.
      COMPLEX(dp), ALLOCATABLE :: inverse_values(:, :, :), inverse_gradients(:, :, :, :)
   END TYPE qs_upf_inverse_atom_grid_type

   PUBLIC :: qs_upf_inverse_atom_grid_ref_type, qs_upf_inverse_atom_grid_type
   PUBLIC :: qs_upf_inverse_atom_grid_prepare, qs_upf_inverse_atom_grid_response
   PUBLIC :: qs_upf_reference_rotations

CONTAINS

! **************************************************************************************************
!> \brief Select the shared angular reference for periodic USPP and PAW quadratures.
!>        The reference cell is independent of mesh reduction and survives geometry updates.
!> \param cell Current physical cell, used to select the periodic domain.
!> \param reference Cartesian lattice vectors of the fixed angular reference.
!> \param kinds Original potentials of the physical kinds.
!> \param atom_xc Whether the USPP functional includes atom-centered XC quadrature.
!> \param rotations Reference point group, unallocated when no angular average is needed.
! **************************************************************************************************
   SUBROUTINE qs_upf_reference_rotations(cell, reference, kinds, atom_xc, rotations)
      TYPE(cell_type), INTENT(IN)                       :: cell
      REAL(dp), INTENT(IN)                               :: reference(3, 3)
      TYPE(qs_kind_type), INTENT(IN)                     :: kinds(:)
      LOGICAL, INTENT(IN)                                :: atom_xc
      REAL(dp), ALLOCATABLE, INTENT(OUT)                  :: rotations(:, :, :)

      INTEGER                                            :: kind

      IF (ANY(cell%perd /= 1)) RETURN
      DO kind = 1, SIZE(kinds)
         IF (kinds(kind)%ghost) CYCLE
         IF (.NOT. ASSOCIATED(kinds(kind)%sgp_potential)) CYCLE
         IF (.NOT. ALLOCATED(kinds(kind)%sgp_potential%upf_data)) CYCLE
         ASSOCIATE (pot => kinds(kind)%sgp_potential%upf_data)
            IF (.NOT. (pot%is_paw .OR. (atom_xc .AND. pot%is_ultrasoft))) CYCLE
         END ASSOCIATE
         CALL crys_sym_lattice_rotations(reference, rotations)
         RETURN
      END DO
   END SUBROUTINE qs_upf_reference_rotations

! **************************************************************************************************
!> \brief Build distributed atom rows from native radial/Lebedev quadratures.
!>        Cartesian radial offsets and radial support are independent of geometry. A C2 window
!>        tapers the final quarter of each finite support. This defines a finite residual quadrature,
!>        not an exact partition of all space. The periodic Becke weights use a smooth image window.
!> \param atom_coords Cartesian base-atom positions.
!> \param cell Periodic cell, also supporting nonperiodic directions.
!> \param atom_kinds Kind index of each atom.
!> \param grids Borrowed initialized native quadratures, indexed by kind.
!> \param support_radii Fixed positive support radius of each atom.
!> \param comm Communicator distributing rows in their global radial/angular order.
!> \param state Owned local quadrature and finite image layouts.
!> \param atom_blocks Keep complete atom blocks on one rank for nonlocal atomic functionals.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_atom_grid_prepare(atom_coords, cell, atom_kinds, grids, support_radii, comm, state, &
                                               atom_blocks)
      REAL(dp), INTENT(IN)                               :: atom_coords(:, :)
      TYPE(cell_type), POINTER                           :: cell
      INTEGER, INTENT(IN)                                :: atom_kinds(:)
      TYPE(qs_upf_inverse_atom_grid_ref_type), &
         INTENT(IN)                                      :: grids(:)
      REAL(dp), INTENT(IN)                               :: support_radii(:)
      TYPE(mp_comm_type), INTENT(IN)                     :: comm
      TYPE(qs_upf_inverse_atom_grid_type), INTENT(OUT)   :: state
      LOGICAL, INTENT(IN), OPTIONAL                      :: atom_blocks

      INTEGER                                            :: atom, atom_row, ia, ir, kind, nang, &
                                                            nlocal, row, target
      LOGICAL                                            :: complete_atoms, owned
      REAL(dp)                                           :: radius, taper, u, weight, window, &
                                                            window_gradient(3)
      TYPE(grid_atom_type), POINTER                      :: grid

      CPASSERT(ASSOCIATED(cell))
      CPASSERT(SIZE(atom_coords, 1) == 3)
      CPASSERT(SIZE(atom_coords, 2) == SIZE(atom_kinds))
      CPASSERT(SIZE(support_radii) == SIZE(atom_kinds))
      CPASSERT(ALL(support_radii > 0.0_dp))
      CPASSERT(ALL(atom_kinds >= 1 .AND. atom_kinds <= SIZE(grids)))
      complete_atoms = .FALSE.
      IF (PRESENT(atom_blocks)) complete_atoms = atom_blocks
      state%natom = SIZE(atom_kinds)
      state%inverse_cell = cell%h_inv
      IF (ANY(cell%perd /= 0)) THEN
         ! Keep every atom quadrature row inside the plateau of its own image window.
         CALL build_periodic_atom_image_layout(atom_coords, cell, 1, layout=state%partition, target_image=target, &
                                               cutoff_radius=2.0_dp*MAXVAL(support_radii))
      ELSE
         CALL build_periodic_atom_image_layout(atom_coords, cell, 1, layout=state%partition, target_image=target)
      END IF
      CPASSERT(target == 1)
      nlocal = 0
      row = 0
      DO atom = 1, SIZE(atom_kinds)
         kind = atom_kinds(atom)
         grid => grids(kind)%grid
         CPASSERT(ASSOCIATED(grid))
         CPASSERT(ASSOCIATED(grid%rad) .AND. ASSOCIATED(grid%weight))
         nang = grid%ng_sphere
         IF (ALLOCATED(grids(kind)%angular_weights)) nang = SIZE(grids(kind)%angular_weights)
         DO ir = 1, grid%nr
            IF (grid%rad(ir) >= support_radii(atom)) CYCLE
            DO ia = 1, nang
               row = row + 1
               owned = MODULO(row - 1, comm%num_pe) == comm%mepos
               IF (complete_atoms) owned = MODULO(atom - 1, comm%num_pe) == comm%mepos
               IF (owned) nlocal = nlocal + 1
            END DO
         END DO
      END DO
      ALLOCATE (state%atoms(nlocal), state%atom_rows(nlocal), state%points(3, nlocal), &
                state%weights(nlocal), state%base_weights(nlocal))
      row = 0
      nlocal = 0
      DO atom = 1, SIZE(atom_kinds)
         kind = atom_kinds(atom)
         grid => grids(kind)%grid
         nang = grid%ng_sphere
         IF (ALLOCATED(grids(kind)%angular_weights)) nang = SIZE(grids(kind)%angular_weights)
         atom_row = 0
         DO ir = 1, grid%nr
            radius = grid%rad(ir)
            IF (radius >= support_radii(atom)) CYCLE
            u = MAX(0.0_dp, (radius/support_radii(atom) - 0.75_dp)/0.25_dp)
            taper = (1.0_dp - u)**3*(1.0_dp + 3.0_dp*u + 6.0_dp*u*u)
            DO ia = 1, nang
               row = row + 1
               atom_row = atom_row + 1
               owned = MODULO(row - 1, comm%num_pe) == comm%mepos
               IF (complete_atoms) owned = MODULO(atom - 1, comm%num_pe) == comm%mepos
               IF (.NOT. owned) CYCLE
               nlocal = nlocal + 1
               state%atoms(nlocal) = atom
               state%atom_rows(nlocal) = atom_row
               IF (ALLOCATED(grids(kind)%angular_weights)) THEN
                  state%points(:, nlocal) = atom_coords(:, atom) + radius*grids(kind)%directions(:, ia)
                  state%base_weights(nlocal) = grid%wr(ir)*grids(kind)%angular_weights(ia)*taper
               ELSE
                  state%points(:, nlocal) = atom_coords(:, atom) + radius* &
                                            [grid%sin_pol(ia)*grid%cos_azi(ia), &
                                             grid%sin_pol(ia)*grid%sin_azi(ia), grid%cos_pol(ia)]
                  state%base_weights(nlocal) = grid%weight(ia, ir)*taper
               END IF
               CALL periodic_atom_image_partition_from_layout(state%points(:, nlocal), state%partition, atom, weight)
               state%weights(nlocal) = state%base_weights(nlocal)*weight
            END DO
         END DO
      END DO
      IF (ANY(cell%perd == 0)) THEN
         state%background_axes = 1 - cell%perd
         ALLOCATE (state%background_weights(nlocal))
         DO row = 1, nlocal
            CALL auxiliary_cell_window(MATMUL(state%inverse_cell, state%points(:, row)), &
                                       state%background_axes, window, window_gradient)
            state%background_weights(row) = state%weights(row)*window
         END DO
      END IF
   END SUBROUTINE qs_upf_inverse_atom_grid_prepare

! **************************************************************************************************
!> \brief Differentiate the same finite quadrature at fixed radial offsets and image membership.
!>        Atom grids move with their centers. There is no cell-volume Jacobian for Cartesian
!>        atom quadrature. Cell dependence enters through atom positions and lattice images.
!> \param state Prepared local quadrature, optionally retaining its partition derivatives.
!> \param atom_response Cartesian response of each base atom.
!> \param cell_response Response of the cell's lattice columns.
!> \param points_response Local point response, indexed (axis,row).
!> \param weights_response Local response of the complete quadrature measure.
!> \param background_weights_response Optional response of the once-per-cell FFT background measure.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_atom_grid_response(state, atom_response, cell_response, points_response, weights_response, &
                                                background_weights_response)
      TYPE(qs_upf_inverse_atom_grid_type), INTENT(INOUT) :: state
      REAL(dp), INTENT(IN)                               :: atom_response(:, :), cell_response(3, 3)
      REAL(dp), INTENT(OUT)                              :: points_response(:, :), &
                                                            weights_response(:)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: background_weights_response(:)

      INTEGER                                            :: atom, nrow, row, status
      INTEGER(int_8)                                     :: bytes_per_row
      LOGICAL                                            :: cached, store_derivatives
      REAL(dp)                                           :: dstrain(3, 3), fractional(3), &
                                                            fractional_response(3), strain(3, 3), &
                                                            value, weight, window, window_gradient(3)
      REAL(dp), ALLOCATABLE                              :: datom(:, :)

      CPASSERT(ALL(SHAPE(atom_response) == [3, state%natom]))
      CPASSERT(ALL(SHAPE(points_response) == SHAPE(state%points)))
      CPASSERT(SIZE(weights_response) == SIZE(state%weights))
      nrow = SIZE(state%atoms)
      cached = ALLOCATED(state%partition_derivatives)
      store_derivatives = .FALSE.
      bytes_per_row = (3_int_8*state%natom + 9)*STORAGE_SIZE(0.0_dp)/8
      IF (.NOT. cached .AND. nrow > 0 .AND. nrow <= partition_cache_bytes/bytes_per_row) THEN
         ALLOCATE (state%partition_derivatives(3, state%natom + 3, nrow), STAT=status)
         store_derivatives = status == 0
      END IF
      ALLOCATE (datom(3, state%natom))
      strain = MATMUL(cell_response, state%inverse_cell)
      DO row = 1, nrow
         atom = state%atoms(row)
         points_response(:, row) = atom_response(:, atom)
         IF (cached) THEN
            datom(:, :) = state%partition_derivatives(:, :state%natom, row)
            dstrain(:, :) = state%partition_derivatives(:, state%natom + 1:, row)
         ELSE
            CALL periodic_atom_image_partition_from_layout(state%points(:, row), state%partition, atom, weight, datom, dstrain)
            IF (store_derivatives) THEN
               state%partition_derivatives(:, :state%natom, row) = datom
               state%partition_derivatives(:, state%natom + 1:, row) = dstrain
            END IF
         END IF
         ! Translation invariance gives the point derivative as minus the sum of atom derivatives.
         value = SUM(datom*(atom_response - SPREAD(points_response(:, row), 2, state%natom))) + SUM(dstrain*strain)
         weights_response(row) = state%base_weights(row)*value
      END DO
      IF (PRESENT(background_weights_response)) THEN
         CPASSERT(SIZE(background_weights_response) == SIZE(state%weights))
         background_weights_response(:) = weights_response
         IF (ALLOCATED(state%background_weights)) THEN
            DO row = 1, SIZE(state%weights)
               fractional = MATMUL(state%inverse_cell, state%points(:, row))
               fractional_response = MATMUL(state%inverse_cell, points_response(:, row) - &
                                            MATMUL(cell_response, fractional))
               CALL auxiliary_cell_window(fractional, state%background_axes, window, window_gradient)
               background_weights_response(row) = weights_response(row)*window + &
                                                  state%weights(row)*DOT_PRODUCT(window_gradient, fractional_response)
            END DO
         END IF
      END IF
   END SUBROUTINE qs_upf_inverse_atom_grid_response

! **************************************************************************************************
!> \brief C2 partition of one FFT cell along the physically nonperiodic directions.
!>        Each one-dimensional factor is one on [1/4,3/4], vanishes outside [-1/4,5/4],
!>        and sums to one over integer translates. The product therefore counts a periodic
!>        auxiliary field once without imposing a discontinuous physical-cell boundary.
!> \param fractional Point in the defining cell's fractional coordinates, without wrapping.
!> \param axes Directions needing an auxiliary-cell partition.
!> \param window Product of the selected one-dimensional factors.
!> \param gradient Derivative with respect to the fractional point.
! **************************************************************************************************
   PURE SUBROUTINE auxiliary_cell_window(fractional, axes, window, gradient)
      REAL(dp), INTENT(IN)                               :: fractional(3)
      INTEGER, INTENT(IN)                                :: axes(3)
      REAL(dp), INTENT(OUT)                              :: window, gradient(3)

      INTEGER                                            :: axis, other
      REAL(dp)                                           :: derivatives(3), factors(3), sign, u

      factors = 1.0_dp
      derivatives = 0.0_dp
      DO axis = 1, 3
         IF (axes(axis) == 0) CYCLE
         IF (fractional(axis) < 0.25_dp) THEN
            u = 2.0_dp*(fractional(axis) + 0.25_dp)
            sign = 2.0_dp
         ELSE IF (fractional(axis) > 0.75_dp) THEN
            u = 2.0_dp*(1.25_dp - fractional(axis))
            sign = -2.0_dp
         ELSE
            CYCLE
         END IF
         IF (u <= 0.0_dp) THEN
            factors(axis) = 0.0_dp
         ELSE
            factors(axis) = u**3*(10.0_dp - 15.0_dp*u + 6.0_dp*u*u)
            derivatives(axis) = sign*30.0_dp*u*u*(1.0_dp - u)**2
         END IF
      END DO
      window = PRODUCT(factors)
      gradient = derivatives
      DO axis = 1, 3
         DO other = 1, 3
            IF (axis /= other) gradient(axis) = gradient(axis)*factors(other)
         END DO
      END DO
   END SUBROUTINE auxiliary_cell_window

END MODULE qs_upf_inverse_atom_grid
