!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM periodic_atom_partition_unittest
   USE cell_types,                      ONLY: cell_type
   USE kinds,                           ONLY: dp
   USE libcp2k,                         ONLY: cp2k_finalize,&
                                              cp2k_init
   USE mathlib,                         ONLY: invmat
   USE skala_gpw_features,              ONLY: periodic_atom_image_partition

   IMPLICIT NONE

   INTEGER :: atom, axis, component, info, replica, sample
   REAL(dp), PARAMETER :: step = 1.e-5_dp
   REAL(dp) :: atoms(3, 2), super_atoms(3, 6), point(3), hmat(3, 3), transform(3, 3), &
               shifted(3, 2), da(3, 2), ds(3, 3), super_da(3, 6), super_ds(3, 3), &
               total_ds(3, 3), super_total_ds(3, 3), weight, super_weight, plus, minus, &
               default_weight, explicit_weight, max_error, unmatched_error
   TYPE(cell_type), POINTER :: cell, supercell

   CALL cp2k_init()
   ALLOCATE (cell, supercell)
   cell%perd = 1
   supercell%perd = 1
   hmat = RESHAPE([3.0_dp, 0.0_dp, 0.0_dp, 0.2_dp, 3.1_dp, 0.0_dp, &
                   -0.1_dp, 0.3_dp, 2.9_dp], [3, 3])
   cell%hmat = hmat
   cell%h_inv = cell%hmat
   CALL invmat(cell%h_inv, info)
   IF (info /= 0) ERROR STOP "Singular partition test cell"
   supercell%hmat = hmat
   supercell%hmat(:, 1) = 3.0_dp*hmat(:, 1)
   supercell%h_inv = supercell%hmat
   CALL invmat(supercell%h_inv, info)
   IF (info /= 0) ERROR STOP "Singular partition test cell"
   atoms(:, 1) = [1.1_dp, 1.2_dp, 1.3_dp]
   atoms(:, 2) = [1.85_dp, 1.42_dp, 1.49_dp]
   DO replica = 0, 2
      super_atoms(:, 2*replica + 1:2*replica + 2) = atoms + &
                                                    SPREAD(REAL(replica, dp)*hmat(:, 1), 2, 2)
   END DO
   max_error = 0.0_dp
   unmatched_error = 0.0_dp
   DO sample = 1, 3
      point = atoms(:, 1) + [0.2_dp - 0.4_dp*sample, 0.15_dp*sample, 0.11_dp*sample]
      CALL periodic_atom_image_partition(point, atoms, cell, 1, default_weight)
      CALL periodic_atom_image_partition(point, atoms, cell, 1, explicit_weight, image_shell=[1, 1, 1])
      CALL check(default_weight - explicit_weight, 1.e-14_dp)
      CALL periodic_atom_image_partition(point, atoms, cell, 1, weight, da, ds, image_shell=[4, 1, 1])
      CALL periodic_atom_image_partition(point, super_atoms, supercell, 1, super_weight, super_da, super_ds)
      CALL check(weight - super_weight, 2.e-13_dp)
      unmatched_error = MAX(unmatched_error, ABS(default_weight - super_weight))
      DO atom = 1, 2
         DO axis = 1, 3
            CALL check(da(axis, atom) - SUM(super_da(axis, atom:6:2)), 2.e-12_dp)
            shifted = atoms
            shifted(axis, atom) = shifted(axis, atom) + step
            CALL periodic_atom_image_partition(point, shifted, cell, 1, plus, image_shell=[4, 1, 1])
            shifted(axis, atom) = shifted(axis, atom) - 2.0_dp*step
            CALL periodic_atom_image_partition(point, shifted, cell, 1, minus, image_shell=[4, 1, 1])
            CALL check(da(axis, atom) - (plus - minus)/(2.0_dp*step), 2.e-8_dp)
         END DO
      END DO
      ! The returned strain contains image translations only. Include affine motion
      ! of the base atoms and the target point for the physical strain derivative.
      DO axis = 1, 3
         DO component = 1, 3
            total_ds(axis, component) = ds(axis, component) + &
                                        SUM(da(axis, :)*(atoms(component, :) - point(component)))
            super_total_ds(axis, component) = super_ds(axis, component) + &
                                              SUM(super_da(axis, :)*(super_atoms(component, :) - point(component)))
            CALL check(total_ds(axis, component) - super_total_ds(axis, component), 2.e-12_dp)
            transform = 0.0_dp
            transform(1, 1) = 1.0_dp
            transform(2, 2) = 1.0_dp
            transform(3, 3) = 1.0_dp
            transform(axis, component) = transform(axis, component) + step
            cell%hmat = MATMUL(transform, hmat)
            cell%h_inv = cell%hmat
            CALL invmat(cell%h_inv, info)
            IF (info /= 0) ERROR STOP "Singular partition test cell"
            CALL periodic_atom_image_partition(MATMUL(transform, point), MATMUL(transform, atoms), &
                                               cell, 1, plus, image_shell=[4, 1, 1])
            transform(axis, component) = transform(axis, component) - 2.0_dp*step
            cell%hmat = MATMUL(transform, hmat)
            cell%h_inv = cell%hmat
            CALL invmat(cell%h_inv, info)
            IF (info /= 0) ERROR STOP "Singular partition test cell"
            CALL periodic_atom_image_partition(MATMUL(transform, point), MATMUL(transform, atoms), &
                                               cell, 1, minus, image_shell=[4, 1, 1])
            CALL check(total_ds(axis, component) - (plus - minus)/(2.0_dp*step), 2.e-8_dp)
         END DO
      END DO
      cell%hmat = hmat
      cell%h_inv = cell%hmat
      CALL invmat(cell%h_inv, info)
      IF (info /= 0) ERROR STOP "Singular partition test cell"
   END DO
   IF (unmatched_error < 1.e-5_dp) ERROR STOP "Partition fixture does not resolve image-shell sensitivity"
   PRINT '(A,ES16.7)', 'Periodic atom partition maximum absolute error: ', max_error
   DEALLOCATE (cell, supercell)
   CALL cp2k_finalize()

CONTAINS

! **************************************************************************************************
!> \brief Check partition equivalence or an independent finite-difference derivative.
!> \param error ...
!> \param tolerance ...
! **************************************************************************************************
   SUBROUTINE check(error, tolerance)
      REAL(dp), INTENT(IN)                               :: error, tolerance

      max_error = MAX(max_error, ABS(error))
      IF (ABS(error) > tolerance) THEN
         PRINT *, 'Periodic atom partition error and tolerance: ', error, tolerance
         ERROR STOP "Periodic atom partition check failed"
      END IF
   END SUBROUTINE check

END PROGRAM periodic_atom_partition_unittest
