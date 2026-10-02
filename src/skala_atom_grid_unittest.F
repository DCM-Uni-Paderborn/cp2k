!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM skala_atom_grid_unittest
   USE kinds,                           ONLY: dp
   USE libcp2k,                         ONLY: cp2k_finalize,&
                                              cp2k_init
   USE mathconstants,                   ONLY: pi
   USE mathlib,                         ONLY: det_3x3,&
                                              inv_3x3
   USE message_passing,                 ONLY: mp_comm_self,&
                                              mp_comm_type,&
                                              mp_comm_world
   USE skala_atom_grid,                 ONLY: skala_atom_grid_build,&
                                              skala_partition_init,&
                                              skala_partition_type,&
                                              skala_partition_weights,&
                                              skala_radial_grid

   IMPLICIT NONE
   TYPE(mp_comm_type) :: comm, single
   TYPE(skala_partition_type) :: partition, translated, displaced
   REAL(dp) :: lattice(3, 3), inverse(3, 3), centres(3, 2), shifted(3, 2), x(3), offset(3), &
               p, d, p2, d2, total, scalar(4), serial(4), coarse_error, error, &
               radius(12), radial_weight(12), t(12), weight(12), moments(24)
   REAL(dp), ALLOCATABLE :: points(:, :), ew(:), dw(:), bc(:, :)
   INTEGER, ALLOCATABLE :: bs(:)
   INTEGER :: i, j, k, a, point, level, nr, na, blocks, empty, descriptor_only

   CALL cp2k_init()
   comm = mp_comm_world
   CALL comm%init()
   single = mp_comm_self
   CALL single%init()
   CALL skala_radial_grid(12, 1.7_dp, radius, radial_weight)
   t = radius/(radius + 1.7_dp)
   weight = radial_weight*(1.0_dp - t)**2/(radius**2*1.7_dp)
   DO i = 0, 23
      moments(i + 1) = SUM(weight*t**i) - 1.0_dp/REAL(i + 1, dp)
   END DO
   IF (MAXVAL(ABS(moments)) > 2.0e-14_dp) ERROR STOP "Gauss-Legendre moments or radial Jacobian are incorrect"

   lattice = RESHAPE([8.0_dp, 0.0_dp, 0.0_dp, 1.3_dp, 7.0_dp, 0.0_dp, -0.7_dp, 0.8_dp, 9.0_dp], [3, 3])
   inverse = inv_3x3(lattice)
   centres(:, 1) = MATMUL(lattice, [0.13_dp, 0.31_dp, 0.21_dp])
   centres(:, 2) = MATMUL(lattice, [0.68_dp, 0.73_dp, 0.49_dp])
   CALL skala_partition_init(partition, lattice, centres, 1)
   shifted = centres
   shifted(:, 2) = centres(:, 2) + MATMUL(lattice, [2.0_dp, -1.0_dp, 1.0_dp])
   CALL skala_partition_init(translated, lattice, shifted, 1)
   shifted(:, 2) = centres(:, 2) + [0.2_dp, -0.3_dp, 0.1_dp]
   CALL skala_partition_init(displaced, lattice, shifted, 1)
   DO point = 1, 6
      x = MATMUL(lattice, [0.21_dp*point, 0.11_dp*point, 0.07_dp*point])
      IF (point == 5) x = centres(:, 1) + lattice(:, 1)
      IF (point == 6) x = centres(:, 1) + 2*lattice(:, 1) - [1.0e-7_dp, 0.0_dp, 0.0_dp]
      total = 0.0_dp
      DO a = 1, 2
         DO k = -3, 3
            DO j = -3, 3
               DO i = -3, 3
                  offset = x - centres(:, a) - MATMUL(lattice, REAL([i, j, k], dp))
                  CALL skala_partition_weights(partition, a, offset, p, d)
                  IF (p < 0.0_dp .OR. p > 1.0_dp .OR. d < 0.0_dp .OR. d > 1.0_dp) &
                     ERROR STOP "Invalid partition weight"
                  total = total + p
               END DO
            END DO
         END DO
      END DO
      IF (ABS(total - 1.0_dp) > 5.0e-13_dp) ERROR STOP "Periodic partition is not a partition of unity"
      offset = x - centres(:, 1)
      CALL skala_partition_weights(partition, 1, offset, p, d)
      CALL skala_partition_weights(translated, 1, offset, p2, d2)
      IF (ABS(p - p2) + ABS(d - d2) > 5.0e-13_dp) ERROR STOP "Rewrapping an atom changes the partition"
      CALL skala_partition_weights(displaced, 1, offset, p2, d2)
      IF (ABS(d - d2) > 5.0e-13_dp) ERROR STOP "Foreign atoms change the descriptor domain"
   END DO
   CALL skala_partition_weights(partition, 1, centres(:, 2) - centres(:, 1), p, d)
   IF (p /= 0.0_dp .OR. d <= 0.0_dp) ERROR STOP "A foreign nucleus must remain a descriptor-only target"

   coarse_error = HUGE(1.0_dp)
   DO level = 1, 2
      nr = 40
      na = 110
      IF (level == 2) THEN
         nr = 160
         na = 974
      END IF
      CALL skala_atom_grid_build(comm, lattice, centres, nr, na, 1, 1.0_dp, points, ew, dw, bs, bc)
      IF (SUM(bs) /= SIZE(ew) .OR. SIZE(dw) /= SIZE(ew)) ERROR STOP "Incomplete atom blocks"
      blocks = SIZE(bs)
      empty = 0
      IF (SIZE(ew) == 0) empty = 1
      descriptor_only = COUNT(ew == 0.0_dp .AND. dw > 0.0_dp)
      CALL comm%sum(blocks)
      CALL comm%sum(empty)
      CALL comm%sum(descriptor_only)
      IF (blocks /= 2 .OR. empty /= MAX(0, comm%num_pe - 2)) ERROR STOP "Atom-block ownership is not unique"
      scalar = integrals(points, ew)
      CALL comm%sum(scalar)
      ! The exact integrals are volume, 0, 0, and volume/2.
      serial = [ABS(det_3x3(lattice)), 0.0_dp, 0.0_dp, 0.5_dp*ABS(det_3x3(lattice))]
      error = MAXVAL(ABS(scalar - serial))/ABS(det_3x3(lattice))
      IF (comm%mepos == 0) WRITE (*, '(A,2I6,5ES16.7)') "Atom-grid integration (nr,na,integrals,error): ", nr, na, scalar, error
      IF (level == 1) coarse_error = error
      IF (level == 2 .AND. (error >= coarse_error .OR. error > 1.0e-6_dp)) &
         ERROR STOP "Periodic volume/Fourier quadrature did not converge"
      IF (level == 1) THEN
         CALL skala_atom_grid_build(single, lattice, centres, nr, na, 1, 1.0_dp, points, ew, dw, bs, bc)
         serial = integrals(points, ew)
         IF (MAXVAL(ABS(scalar - serial)) > 1.0e-10_dp) ERROR STOP "Distributed atom-grid integrals differ"
      END IF
   END DO
   shifted = centres
   CALL skala_radial_grid(12, 1.0_dp, radius, radial_weight)
   shifted(:, 2) = centres(:, 1) + [radius(6), 0.0_dp, 0.0_dp]
   CALL skala_atom_grid_build(comm, lattice, shifted, 12, 6, 1, 1.0_dp, points, ew, dw, bs, bc)
   descriptor_only = COUNT(ew == 0.0_dp .AND. dw > 0.0_dp)
   CALL comm%sum(descriptor_only)
   IF (descriptor_only == 0) ERROR STOP "Descriptor-only quadrature rows have been discarded"
   IF (comm%mepos == 0) WRITE (*, '(A)') "Periodic Skala atom-grid checks passed"
   CALL cp2k_finalize()
CONTAINS
   FUNCTION integrals(positions, weights) RESULT(values)
      REAL(dp), INTENT(IN)                               :: positions(:, :), weights(:)
      REAL(dp)                                           :: values(4)

      INTEGER                                            :: row
      REAL(dp)                                           :: f, s(3)

      values = 0.0_dp
      DO row = 1, SIZE(weights)
         s = MATMUL(inverse, positions(:, row))
         f = COS(2*pi*(s(1) - s(2)))
         values = values + weights(row)*[1.0_dp, COS(2*pi*s(1)), SIN(2*pi*s(2)), f*f]
      END DO
   END FUNCTION
END PROGRAM
