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
      REAL(dp), ALLOCATABLE :: points(:, :), weights(:), base_weights(:)
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
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_atom_grid_prepare(atom_coords, cell, atom_kinds, grids, support_radii, comm, state, image_shell)
      REAL(dp), INTENT(IN)                               :: atom_coords(:, :)
      TYPE(cell_type), POINTER                           :: cell
      INTEGER, INTENT(IN)                                :: atom_kinds(:)
      TYPE(qs_upf_inverse_atom_grid_ref_type), &
         INTENT(IN)                                      :: grids(:)
      REAL(dp), INTENT(IN)                               :: support_radii(:)
      TYPE(mp_comm_type), INTENT(IN)                     :: comm
      TYPE(qs_upf_inverse_atom_grid_type), INTENT(OUT)   :: state
      INTEGER, INTENT(IN), OPTIONAL                      :: image_shell(3)

      INTEGER                                            :: atom, ia, idir, ir, kind, nlocal, row, &
                                                            shell(3)
      REAL(dp)                                           :: radius, taper, u, weight
      TYPE(grid_atom_type), POINTER                      :: grid

      CPASSERT(ASSOCIATED(cell))
      CPASSERT(SIZE(atom_coords, 1) == 3)
      CPASSERT(SIZE(atom_coords, 2) == SIZE(atom_kinds))
      CPASSERT(SIZE(support_radii) == SIZE(atom_kinds))
      CPASSERT(ALL(support_radii > 0.0_dp))
      CPASSERT(ALL(atom_kinds >= 1 .AND. atom_kinds <= SIZE(grids)))
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
               IF (MODULO(row - 1, comm%num_pe) == comm%mepos) nlocal = nlocal + 1
            END DO
         END DO
      END DO
      ALLOCATE (state%atoms(nlocal), state%points(3, nlocal), state%weights(nlocal), state%base_weights(nlocal))
      row = 0
      nlocal = 0
      DO atom = 1, SIZE(atom_kinds)
         grid => grids(atom_kinds(atom))%grid
         DO ir = 1, grid%nr
            radius = grid%rad(ir)
            IF (radius >= support_radii(atom)) CYCLE
            u = MAX(0.0_dp, (radius/support_radii(atom) - 0.75_dp)/0.25_dp)
            taper = (1.0_dp - u)**3*(1.0_dp + 3.0_dp*u + 6.0_dp*u*u)
            DO ia = 1, grid%ng_sphere
               row = row + 1
               IF (MODULO(row - 1, comm%num_pe) /= comm%mepos) CYCLE
               nlocal = nlocal + 1
               state%atoms(nlocal) = atom
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
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_atom_grid_response(state, atom_response, cell_response, points_response, weights_response)
      TYPE(qs_upf_inverse_atom_grid_type), INTENT(IN)    :: state
      REAL(dp), INTENT(IN)                               :: atom_response(:, :), cell_response(3, 3)
      REAL(dp), INTENT(OUT)                              :: points_response(:, :), &
                                                            weights_response(:)

      INTEGER                                            :: atom, image, row
      REAL(dp)                                           :: dpoint(3), image_response(3), value, &
                                                            weight
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
   END SUBROUTINE qs_upf_inverse_atom_grid_response

END MODULE qs_upf_inverse_atom_grid
