!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
!> \brief Atom-centered quadrature and its geometry response for inverse mixed XC.
MODULE qs_upf_inverse_atom_grid
   USE cell_types,                      ONLY: cell_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_comm_type
   USE qs_atom_grid_partition,          ONLY: build_periodic_atom_image_layout,&
                                              periodic_atom_image_partition_from_layout
   USE qs_grid_atom,                    ONLY: grid_atom_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   TYPE qs_upf_inverse_atom_grid_ref_type
      TYPE(grid_atom_type), POINTER :: grid => NULL()
   END TYPE qs_upf_inverse_atom_grid_ref_type

   TYPE atom_image_layout_type
      INTEGER :: target_image = 0
      INTEGER, ALLOCATABLE :: atoms(:), shifts(:, :)
      REAL(dp), ALLOCATABLE :: coordinates(:, :)
   END TYPE atom_image_layout_type

   TYPE qs_upf_inverse_atom_grid_type
      INTEGER, ALLOCATABLE :: atoms(:)
      ! Row within the complete atom quadrature, independent of MPI distribution.
      INTEGER, ALLOCATABLE :: atom_rows(:)
      REAL(dp), ALLOCATABLE :: points(:, :), weights(:), base_weights(:)
      ! The FFT background is integrated once per auxiliary cell also for isolated electrostatics.
      REAL(dp), ALLOCATABLE :: background_weights(:)
      REAL(dp) :: inverse_cell(3, 3) = 0.0_dp
      INTEGER :: background_axes(3) = 0
      ! Native inverse frames and this quadrature share one geometry lifetime.
      COMPLEX(dp), ALLOCATABLE :: inverse_values(:, :, :), inverse_gradients(:, :, :, :)
      TYPE(atom_image_layout_type), ALLOCATABLE :: layouts(:)
   END TYPE qs_upf_inverse_atom_grid_type

   PUBLIC :: qs_upf_inverse_atom_grid_ref_type, qs_upf_inverse_atom_grid_type
   PUBLIC :: qs_upf_inverse_atom_grid_prepare, qs_upf_inverse_atom_grid_response

CONTAINS

! **************************************************************************************************
!> \brief Build distributed atom rows from native radial/Lebedev quadratures.
!>        Cartesian radial offsets and radial support are independent of geometry. A C2 window
!>        tapers the final quarter of each finite support. This defines a finite residual quadrature,
!>        not an exact partition of all space. Its finite image layout is shared by the atom's rows.
!> \param atom_coords Cartesian base-atom positions.
!> \param cell Periodic cell, also supporting nonperiodic directions.
!> \param atom_kinds Kind index of each atom.
!> \param grids Borrowed initialized native quadratures, indexed by kind.
!> \param support_radii Fixed positive support radius of each atom.
!> \param comm Communicator distributing rows in their global radial/angular order.
!> \param state Owned local quadrature and finite image layouts.
!> \param image_shell Optional common image shell for a specified finite approximation.
!> \param atom_blocks Keep complete atom blocks on one rank for nonlocal atomic functionals.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_atom_grid_prepare(atom_coords, cell, atom_kinds, grids, support_radii, comm, state, &
                                               image_shell, atom_blocks)
      REAL(dp), INTENT(IN)                               :: atom_coords(:, :)
      TYPE(cell_type), POINTER                           :: cell
      INTEGER, INTENT(IN)                                :: atom_kinds(:)
      TYPE(qs_upf_inverse_atom_grid_ref_type), &
         INTENT(IN)                                      :: grids(:)
      REAL(dp), INTENT(IN)                               :: support_radii(:)
      TYPE(mp_comm_type), INTENT(IN)                     :: comm
      TYPE(qs_upf_inverse_atom_grid_type), INTENT(OUT)   :: state
      INTEGER, INTENT(IN), OPTIONAL                      :: image_shell(3)
      LOGICAL, INTENT(IN), OPTIONAL                      :: atom_blocks

      INTEGER                                            :: atom, atom_row, ia, idir, ir, kind, &
                                                            nlocal, row, shell(3)
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
      ALLOCATE (state%layouts(SIZE(atom_kinds)))
      nlocal = 0
      row = 0
      DO atom = 1, SIZE(atom_kinds)
         kind = atom_kinds(atom)
         grid => grids(kind)%grid
         CPASSERT(ASSOCIATED(grid))
         CPASSERT(ASSOCIATED(grid%rad) .AND. ASSOCIATED(grid%weight))
         DO idir = 1, 3
            shell(idir) = CEILING(support_radii(atom)*SQRT(SUM(cell%h_inv(idir, :)**2))) + 1
         END DO
         IF (PRESENT(image_shell)) shell = image_shell
         CALL build_periodic_atom_image_layout(atom_coords, cell, atom, image_coords=state%layouts(atom)%coordinates, &
                                               target_image=state%layouts(atom)%target_image, &
                                               image_atom=state%layouts(atom)%atoms, image_shift=state%layouts(atom)%shifts, &
                                               image_shell=shell)
         DO ir = 1, grid%nr
            IF (grid%rad(ir) >= support_radii(atom)) CYCLE
            DO ia = 1, grid%ng_sphere
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
         grid => grids(atom_kinds(atom))%grid
         atom_row = 0
         DO ir = 1, grid%nr
            radius = grid%rad(ir)
            IF (radius >= support_radii(atom)) CYCLE
            u = MAX(0.0_dp, (radius/support_radii(atom) - 0.75_dp)/0.25_dp)
            taper = (1.0_dp - u)**3*(1.0_dp + 3.0_dp*u + 6.0_dp*u*u)
            DO ia = 1, grid%ng_sphere
               row = row + 1
               atom_row = atom_row + 1
               owned = MODULO(row - 1, comm%num_pe) == comm%mepos
               IF (complete_atoms) owned = MODULO(atom - 1, comm%num_pe) == comm%mepos
               IF (.NOT. owned) CYCLE
               nlocal = nlocal + 1
               state%atoms(nlocal) = atom
               state%atom_rows(nlocal) = atom_row
               state%points(:, nlocal) = atom_coords(:, atom) + radius* &
                                         [grid%sin_pol(ia)*grid%cos_azi(ia), &
                                          grid%sin_pol(ia)*grid%sin_azi(ia), grid%cos_pol(ia)]
               state%base_weights(nlocal) = grid%weight(ia, ir)*taper
               CALL periodic_atom_image_partition_from_layout(state%points(:, nlocal), &
                                                              state%layouts(atom)%coordinates, &
                                                              state%layouts(atom)%target_image, weight)
               state%weights(nlocal) = state%base_weights(nlocal)*weight
            END DO
         END DO
      END DO
      IF (ANY(cell%perd == 0)) THEN
         state%inverse_cell = cell%h_inv
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
!> \param state Prepared local quadrature.
!> \param atom_response Cartesian response of each base atom.
!> \param cell_response Response of the cell's lattice columns.
!> \param points_response Local point response, indexed (axis,row).
!> \param weights_response Local response of the complete quadrature measure.
!> \param background_weights_response Optional response of the once-per-cell FFT background measure.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_atom_grid_response(state, atom_response, cell_response, points_response, weights_response, &
                                                background_weights_response)
      TYPE(qs_upf_inverse_atom_grid_type), INTENT(IN)    :: state
      REAL(dp), INTENT(IN)                               :: atom_response(:, :), cell_response(3, 3)
      REAL(dp), INTENT(OUT)                              :: points_response(:, :), &
                                                            weights_response(:)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: background_weights_response(:)

      INTEGER                                            :: atom, image, row
      REAL(dp)                                           :: dpoint(3), fractional(3), &
                                                            fractional_response(3), &
                                                            image_response(3), value, weight, &
                                                            window, window_gradient(3)
      REAL(dp), ALLOCATABLE                              :: dimage(:, :)

      CPASSERT(ALL(SHAPE(atom_response) == [3, SIZE(state%layouts)]))
      CPASSERT(ALL(SHAPE(points_response) == SHAPE(state%points)))
      CPASSERT(SIZE(weights_response) == SIZE(state%weights))
      DO atom = 1, SIZE(state%layouts)
         ALLOCATE (dimage(3, SIZE(state%layouts(atom)%atoms)))
         DO row = 1, SIZE(state%atoms)
            IF (state%atoms(row) /= atom) CYCLE
            points_response(:, row) = atom_response(:, atom)
            CALL periodic_atom_image_partition_from_layout(state%points(:, row), state%layouts(atom)%coordinates, &
                                                           state%layouts(atom)%target_image, weight, dpoint, dimage)
            value = DOT_PRODUCT(dpoint, points_response(:, row))
            DO image = 1, SIZE(state%layouts(atom)%atoms)
               image_response = atom_response(:, state%layouts(atom)%atoms(image)) + &
                                MATMUL(cell_response, REAL(state%layouts(atom)%shifts(:, image), dp))
               value = value + DOT_PRODUCT(dimage(:, image), image_response)
            END DO
            weights_response(row) = state%base_weights(row)*value
         END DO
         DEALLOCATE (dimage)
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
