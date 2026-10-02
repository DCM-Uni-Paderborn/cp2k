!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Fixed periodic atom quadrature with independent energy and descriptor domains.
! **************************************************************************************************
MODULE skala_atom_grid
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: dp
   USE lebedev,                         ONLY: get_number_of_lebedev_grid,&
                                              lebedev_grid
   USE mathconstants,                   ONLY: fourpi
   USE mathlib,                         ONLY: det_3x3,&
                                              inv_3x3
   USE message_passing,                 ONLY: mp_comm_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'skala_atom_grid'

   TYPE seed_environment_type
      REAL(dp), ALLOCATABLE :: displacement(:, :), distance(:), membership(:)
   END TYPE

   TYPE skala_partition_type
      PRIVATE
      INTEGER :: shell = 1
      REAL(dp) :: lattice(3, 3), inverse(3, 3)
      REAL(dp), ALLOCATABLE :: centres(:, :)
      TYPE(seed_environment_type), ALLOCATABLE :: neighbours(:)
      TYPE(seed_environment_type) :: self_neighbours
   END TYPE

   PUBLIC :: skala_partition_type, skala_partition_init, skala_partition_weights, &
             skala_atom_grid_build, skala_radial_grid

CONTAINS

! **************************************************************************************************
!> \brief Cache translated Becke seed environments. Lattice vectors are columns, in bohr.
! **************************************************************************************************
   SUBROUTINE skala_partition_init(partition, lattice, centres, shell)
      TYPE(skala_partition_type), INTENT(OUT)            :: partition
      REAL(dp), INTENT(IN)                               :: lattice(3, 3), centres(:, :)
      INTEGER, INTENT(IN)                                :: shell

      INTEGER                                            :: a

      CPASSERT(shell >= 0)
      CPASSERT(SIZE(centres, 1) == 3 .AND. SIZE(centres, 2) > 0)
      CPASSERT(ALL(ieee_is_finite(lattice)) .AND. ALL(ieee_is_finite(centres)))
      CPASSERT(ABS(det_3x3(lattice)) > 1.0e-12_dp)
      partition%shell = shell
      partition%lattice = lattice
      partition%inverse = inv_3x3(lattice)
      partition%centres = centres
      ALLOCATE (partition%neighbours(SIZE(centres, 2)))
      DO a = 1, SIZE(centres, 2)
         CALL build_environment(partition, a, .FALSE., partition%neighbours(a))
      END DO
      CALL build_environment(partition, 1, .TRUE., partition%self_neighbours)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Each image has the same translated environment, not a finite-cluster edge environment.
! **************************************************************************************************
   SUBROUTINE build_environment(partition, atom, self_only, environment)
      TYPE(skala_partition_type), INTENT(IN)             :: partition
      INTEGER, INTENT(IN)                                :: atom
      LOGICAL, INTENT(IN)                                :: self_only
      TYPE(seed_environment_type), INTENT(OUT)           :: environment

      INTEGER                                            :: b, capacity, count, hi(3), i, j, k, lo(3)
      REAL(dp)                                           :: capacity_real, delta(3), distance, &
                                                            fractional(3), weight
      REAL(dp), ALLOCATABLE                              :: distances(:), offsets(:, :), weights(:)

      capacity_real = REAL(SIZE(partition%centres, 2), dp)*(2.0_dp*partition%shell + 3)**3
      CPASSERT(capacity_real < REAL(HUGE(capacity), dp))
      capacity = SIZE(partition%centres, 2)*(2*partition%shell + 3)**3
      ALLOCATE (offsets(3, capacity), weights(capacity), distances(capacity))
      count = 0
      DO b = 1, SIZE(partition%centres, 2)
         IF (self_only .AND. b /= atom) CYCLE
         delta = MATMUL(partition%inverse, partition%centres(:, b) - partition%centres(:, atom))
         lo = CEILING(-delta - (partition%shell + 1))
         hi = FLOOR(-delta + (partition%shell + 1))
         DO k = lo(3), hi(3)
            DO j = lo(2), hi(2)
               DO i = lo(1), hi(1)
                  IF (b == atom .AND. ALL([i, j, k] == 0)) CYCLE
                  fractional = delta + [i, j, k]
                  weight = compact_window(fractional, partition%shell)
                  IF (weight == 0.0_dp) CYCLE
                  distance = NORM2(MATMUL(partition%lattice, fractional))
                  IF (distance < 1.0e-10_dp) CPABORT("Coincident atoms in the periodic Skala quadrature")
                  count = count + 1
                  offsets(:, count) = MATMUL(partition%lattice, fractional)
                  distances(count) = distance
                  weights(count) = weight
               END DO
            END DO
         END DO
      END DO
      environment%displacement = offsets(:, :count)
      environment%distance = distances(:count)
      environment%membership = weights(:count)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Energy partition of atom at zero image and descriptor turn-on of its self-image partition.
!> \param partition Cached fixed geometry.
!> \param atom Target atom block.
!> \param offset Cartesian displacement from the target atom, in bohr.
!> \param energy All-atom, all-image partition weight.
!> \param descriptor Self-image domain window, independent of all other atoms.
! **************************************************************************************************
   SUBROUTINE skala_partition_weights(partition, atom, offset, energy, descriptor)
      TYPE(skala_partition_type), INTENT(IN)             :: partition
      INTEGER, INTENT(IN)                                :: atom
      REAL(dp), INTENT(IN)                               :: offset(3)
      REAL(dp), INTENT(OUT)                              :: energy, descriptor

      REAL(dp)                                           :: d, z

      CPASSERT(atom >= 1 .AND. atom <= SIZE(partition%centres, 2))
      CPASSERT(ALL(ieee_is_finite(offset)))
      energy = 0.0_dp
      descriptor = 0.0_dp
      IF (compact_window(MATMUL(partition%inverse, offset), partition%shell) == 0.0_dp) RETURN
      d = normalized_seed(partition, atom, offset, .TRUE.)
      IF (SIZE(partition%centres, 2) == 1) THEN
         energy = d
      ELSE
         energy = normalized_seed(partition, atom, offset, .FALSE.)
      END IF
      IF (d >= 1.0e-11_dp) THEN
         descriptor = 1.0_dp
      ELSE IF (d > 1.0e-12_dp) THEN
         z = (d - 1.0e-12_dp)/9.0e-12_dp
         descriptor = z**3*(10.0_dp - 15.0_dp*z + 6.0_dp*z*z)
      END IF
   END SUBROUTINE

! **************************************************************************************************
!> \brief Normalize all nonzero translated seeds using a scaled log sum.
! **************************************************************************************************
   FUNCTION normalized_seed(partition, atom, offset, self_only) RESULT(weight)
      TYPE(skala_partition_type), INTENT(IN)             :: partition
      INTEGER, INTENT(IN)                                :: atom
      REAL(dp), INTENT(IN)                               :: offset(3)
      LOGICAL, INTENT(IN)                                :: self_only
      REAL(dp)                                           :: weight

      INTEGER                                            :: b, hi(3), i, j, k, lo(3)
      REAL(dp)                                           :: fractional(3), largest, relative(3), &
                                                            TARGET, total, value, x(3)

      ! A finite sentinel avoids infinities in builds trapping floating-point exceptions.
      TARGET = -HUGE(1.0_dp)
      largest = TARGET
      total = 0.0_dp
      DO b = 1, SIZE(partition%centres, 2)
         IF (self_only .AND. b /= atom) CYCLE
         relative = offset + partition%centres(:, atom) - partition%centres(:, b)
         fractional = MATMUL(partition%inverse, relative)
         lo = CEILING(fractional - (partition%shell + 1))
         hi = FLOOR(fractional + (partition%shell + 1))
         DO k = lo(3), hi(3)
            DO j = lo(2), hi(2)
               DO i = lo(1), hi(1)
                  x = relative - MATMUL(partition%lattice, REAL([i, j, k], dp))
                  IF (self_only) THEN
                     value = log_seed(partition, partition%self_neighbours, x)
                  ELSE
                     value = log_seed(partition, partition%neighbours(b), x)
                  END IF
                  IF (b == atom .AND. ALL([i, j, k] == 0)) TARGET = value
                  IF (value == -HUGE(1.0_dp)) CYCLE
                  IF (value > largest) THEN
                     total = total*EXP(largest - value) + 1.0_dp
                     largest = value
                  ELSE
                     total = total + EXP(value - largest)
                  END IF
               END DO
            END DO
         END DO
      END DO
      IF (total <= 0.0_dp) CPABORT("Empty periodic Skala partition")
      weight = 0.0_dp
      IF (TARGET /= -HUGE(1.0_dp)) weight = EXP(TARGET - largest)/total
   END FUNCTION

! **************************************************************************************************
!> \brief Compact Becke seed with triple-cubic pair switches and smooth neighbour membership.
! **************************************************************************************************
   FUNCTION log_seed(partition, environment, x) RESULT(value)
      TYPE(skala_partition_type), INTENT(IN)             :: partition
      TYPE(seed_environment_type), INTENT(IN)            :: environment
      REAL(dp), INTENT(IN)                               :: x(3)
      REAL(dp)                                           :: value

      INTEGER                                            :: n, repetition
      REAL(dp)                                           :: c, factor, mu, radius

      value = -HUGE(1.0_dp)
      c = compact_window(MATMUL(partition%inverse, x), partition%shell)
      IF (c == 0.0_dp) RETURN
      value = LOG(c)
      radius = NORM2(x)
      DO n = 1, SIZE(environment%distance)
         mu = (radius - NORM2(x - environment%displacement(:, n)))/environment%distance(n)
         mu = MAX(-1.0_dp, MIN(1.0_dp, mu))
         DO repetition = 1, 3
            mu = 0.5_dp*mu*(3.0_dp - mu*mu)
         END DO
         c = environment%membership(n)
         factor = (1.0_dp - c) + c*MAX(0.0_dp, 0.5_dp*(1.0_dp - mu))
         IF (factor == 0.0_dp) THEN
            value = -HUGE(1.0_dp)
            RETURN
         END IF
         value = value + LOG(factor)
      END DO
   END FUNCTION

! **************************************************************************************************
!> \brief C2 compact window, evaluated without cancellation at the outer edge.
! **************************************************************************************************
   PURE FUNCTION compact_window(fractional, shell) RESULT(weight)
      REAL(dp), INTENT(IN)                               :: fractional(3)
      INTEGER, INTENT(IN)                                :: shell
      REAL(dp)                                           :: weight

      INTEGER                                            :: direction
      REAL(dp)                                           :: z

      weight = 1.0_dp
      DO direction = 1, 3
         z = MAX(0.0_dp, ABS(fractional(direction)) - shell)
         IF (z >= 1.0_dp) THEN
            weight = 0.0_dp
            RETURN
         END IF
         weight = weight*(1.0_dp - z)**3*(1.0_dp + 3.0_dp*z + 6.0_dp*z*z)
      END DO
   END FUNCTION

! **************************************************************************************************
!> \brief Gauss-Legendre radial rule on [0,infinity), including r^2 dr and a scale in bohr.
! **************************************************************************************************
   SUBROUTINE skala_radial_grid(n, scale, radius, weights)
      INTEGER, INTENT(IN)                                :: n
      REAL(dp), INTENT(IN)                               :: scale
      REAL(dp), INTENT(OUT)                              :: radius(n), weights(n)

      INTEGER                                            :: i, info
      REAL(dp)                                           :: eigenvectors(n, n), &
                                                            offdiag(MAX(1, n - 1)), t, &
                                                            work(MAX(1, 2*n - 2))

      CPASSERT(n >= 2 .AND. ieee_is_finite(scale) .AND. scale > 0.0_dp)
      ! Golub-Welsch with the Legendre Jacobi matrix. LAPACK supplies the nodes and weights.
      radius = 0.0_dp
      DO i = 1, n - 1
         offdiag(i) = REAL(i, dp)/SQRT(4.0_dp*REAL(i, dp)**2 - 1.0_dp)
      END DO
      CALL DSTEV('V', n, radius, offdiag, eigenvectors, n, work, info)
      CPASSERT(info == 0)
      DO i = 1, n
         t = 0.5_dp*(1.0_dp + radius(i))
         radius(i) = scale*t/(1.0_dp - t)
         weights(i) = eigenvectors(1, i)**2*radius(i)**2*scale/(1.0_dp - t)**2
      END DO
   END SUBROUTINE

! **************************************************************************************************
!> \brief Build Gauss-Legendre/Lebedev blocks, owned once across the communicator.
!> \note Geometry derivatives are not provided. Never drop a descriptor-only row.
! **************************************************************************************************
   SUBROUTINE skala_atom_grid_build(group, lattice, centres, nradial, nangular, shell, scale, &
                                    points, energy_weights, descriptor_weights, block_sizes, block_centres)
      TYPE(mp_comm_type), INTENT(IN)                     :: group
      REAL(dp), INTENT(IN)                               :: lattice(3, 3), centres(:, :)
      INTEGER, INTENT(IN)                                :: nradial, nangular, shell
      REAL(dp), INTENT(IN)                               :: scale
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: points(:, :), energy_weights(:), &
                                                            descriptor_weights(:)
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: block_sizes(:)
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: block_centres(:, :)

      INTEGER                                            :: a, angular_grid, b, ia, ir, max_rows, &
                                                            na, nb, row, start
      REAL(dp)                                           :: base_weight, d, offset(3), p
      REAL(dp), ALLOCATABLE                              :: descriptor(:), energy(:), &
                                                            positions(:, :), radial_weights(:), &
                                                            radii(:)
      TYPE(skala_partition_type)                         :: partition

      CPASSERT(nradial >= 2 .AND. nangular >= 6 .AND. nangular <= 974)
      CPASSERT(shell >= 0 .AND. scale > 0.0_dp .AND. ieee_is_finite(scale))
      CALL skala_partition_init(partition, lattice, centres, shell)
      angular_grid = get_number_of_lebedev_grid(n=nangular)
      na = lebedev_grid(angular_grid)%n
      IF (na /= nangular) CPABORT("Select an available Lebedev point count for the Skala atom grid")
      IF (ANY(lebedev_grid(angular_grid)%w < 0.0_dp)) THEN
         CPABORT("Skala descriptor weights require a positive Lebedev rule")
      END IF
      nb = 0
      DO a = group%mepos + 1, SIZE(centres, 2), group%num_pe
         nb = nb + 1
      END DO
      CPASSERT(REAL(nb, dp)*nradial*na < REAL(HUGE(max_rows), dp))
      max_rows = nb*nradial*na
      ALLOCATE (positions(3, max_rows), energy(max_rows), descriptor(max_rows), block_sizes(nb), block_centres(3, nb))
      ALLOCATE (radii(nradial), radial_weights(nradial))
      CALL skala_radial_grid(nradial, scale, radii, radial_weights)
      row = 0
      b = 0
      DO a = group%mepos + 1, SIZE(centres, 2), group%num_pe
         b = b + 1
         start = row
         block_centres(:, b) = centres(:, a)
         DO ir = 1, nradial
            DO ia = 1, na
               offset = radii(ir)*lebedev_grid(angular_grid)%r(:, ia)
               CALL skala_partition_weights(partition, a, offset, p, d)
               IF (p == 0.0_dp .AND. d == 0.0_dp) CYCLE
               base_weight = fourpi*radial_weights(ir)*lebedev_grid(angular_grid)%w(ia)
               row = row + 1
               positions(:, row) = centres(:, a) + offset
               energy(row) = base_weight*p
               descriptor(row) = base_weight*d
            END DO
         END DO
         block_sizes(b) = row - start
         IF (block_sizes(b) == 0) CPABORT("Empty Skala atom block. Refine the radial quadrature")
      END DO
      points = positions(:, :row)
      energy_weights = energy(:row)
      descriptor_weights = descriptor(:row)
   END SUBROUTINE
END MODULE skala_atom_grid
