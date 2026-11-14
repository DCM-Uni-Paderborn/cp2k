!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
!> \brief Radial multipole Coulomb map for signed inverse reconstruction densities.
!>        All maps use the same finite radial measure, including the diagonal shell term.
!>        Periodic background and interactions between distinct centers are assembled separately.
MODULE upf_inverse_coulomb
   USE bessel_lib,                      ONLY: bessel0
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kahan_sum,                       ONLY: accurate_sum
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: fourpi,&
                                              pi,&
                                              twopi
   USE mathlib,                         ONLY: invmat,&
                                              real_fourier_aliases
   USE message_passing,                 ONLY: mp_comm_type
   USE orbital_pointers,                ONLY: indco,&
                                              init_orbital_pointers,&
                                              nco,&
                                              ncoset
   USE orbital_transformation_matrices, ONLY: init_spherical_harmonics,&
                                              orbtramat
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   INTEGER, PARAMETER :: reflection_bits(3) = [1, 2, 4]

   TYPE upf_inverse_coulomb_center_type
      REAL(dp) :: position(3) = 0.0_dp
      INTEGER, ALLOCATABLE :: l(:), indices(:)
      REAL(dp), ALLOCATABLE :: radii(:), weights(:), density(:, :), shapes(:, :)
   END TYPE upf_inverse_coulomb_center_type

   TYPE upf_inverse_coulomb_potential_type
      REAL(dp), ALLOCATABLE :: values(:, :)
   END TYPE upf_inverse_coulomb_potential_type

   TYPE center_work_type
      REAL(dp), ALLOCATABLE :: compensation(:, :), residual(:, :), q(:, :), dq(:, :), potential(:, :)
   END TYPE center_work_type

   TYPE fourier_map_type
      COMPLEX(dp), ALLOCATABLE :: values(:, :, :)
   END TYPE fourier_map_type

   TYPE pair_map_type
      REAL(dp), ALLOCATABLE :: values(:, :, :, :)
   END TYPE pair_map_type

   TYPE harmonic_work_type
      INTEGER :: max_degree
      INTEGER, ALLOCATABLE :: first_term(:), powers(:, :)
      REAL(dp), ALLOCATABLE :: coefficients(:)
      REAL(dp), ALLOCATABLE :: values(:), response(:), gradient(:, :)
      REAL(dp), ALLOCATABLE :: sampled_values(:, :), sampled_response(:, :)
   END TYPE harmonic_work_type

   TYPE contained_shell_type
      LOGICAL :: ready(2) = .FALSE.
      REAL(dp) :: source_radius(2) = 0.0_dp
      REAL(dp), ALLOCATABLE :: values(:, :, :), response(:, :, :)
   END TYPE contained_shell_type

   TYPE upf_inverse_coulomb_map_type
      PRIVATE
      INTEGER :: periodicity(3) = [1, 1, 1], grid_shape(3) = 0
      REAL(dp) :: cell(3, 3) = 0.0_dp, volume = 0.0_dp
      REAL(dp), ALLOCATABLE :: gvecs(:, :)
      INTEGER, ALLOCATABLE :: alias_rows(:)
      TYPE(upf_inverse_coulomb_center_type), ALLOCATABLE :: centers(:)
      TYPE(fourier_map_type), ALLOCATABLE :: fourier(:)
      TYPE(pair_map_type), ALLOCATABLE :: pairs(:, :)
   END TYPE upf_inverse_coulomb_map_type

   PUBLIC :: upf_inverse_coulomb_center_type, upf_inverse_coulomb_potential_type
   PUBLIC :: upf_inverse_coulomb_periodic
   PUBLIC :: upf_inverse_coulomb_map_type, upf_inverse_coulomb_prepare
   PUBLIC :: upf_inverse_coulomb_density, upf_inverse_coulomb_density_adjoint
   PUBLIC :: upf_inverse_coulomb_radial, upf_inverse_coulomb_radial_response
   PUBLIC :: upf_inverse_coulomb_compensate, upf_inverse_coulomb_compensation_adjoint
   PUBLIC :: upf_inverse_coulomb_fourier_kernel
   PUBLIC :: upf_inverse_coulomb_pair_kernel
   PUBLIC :: upf_inverse_coulomb_pair_geometry

CONTAINS

! **************************************************************************************************
!> \brief Freeze the Fourier and complete periodic residual-pair maps at a fixed geometry.
!>        Opposite self images are included by reciprocity. Cartesian cells also admit independent pair reflections.
!>        Density and compensation moments remain variable.
!> \param centers Same center geometry and radial/harmonic layout as the subsequent functional.
!> \param cell ...
!> \param gvecs ...
!> \param volume ...
!> \param map Owned maps. Rebuild after changing centers, cell, reciprocal vectors or radial grids.
!> \param periodicity Optional lattice directions to repeat. Defaults to three periodic directions.
!> \param grid_shape Optional full-box real-grid dimensions. Zero or omission retains continuous coefficients.
!> \param comm Optional communicator sharing shell integrals. Center and cell data are replicated.
! **************************************************************************************************
   SUBROUTINE upf_inverse_coulomb_prepare(centers, cell, gvecs, volume, map, periodicity, grid_shape, comm)
      TYPE(upf_inverse_coulomb_center_type), INTENT(IN)  :: centers(:)
      REAL(dp), INTENT(IN)                               :: cell(3, 3), gvecs(:, :), volume
      TYPE(upf_inverse_coulomb_map_type), INTENT(OUT)    :: map
      INTEGER, INTENT(IN), OPTIONAL                      :: periodicity(3), grid_shape(3)
      TYPE(mp_comm_type), INTENT(IN), OPTIONAL           :: comm

      INTEGER                                            :: a, b, channel, h, i, image(3), info, j, &
                                                            k, limits(2, 3), mask, nrank, owner, &
                                                            r, rank, row
      INTEGER, ALLOCATABLE                               :: parity(:, :)
      INTEGER                                            :: miller(3)
      LOGICAL, ALLOCATABLE                               :: aliased(:)
      LOGICAL                                            :: cartesian, reflection_axes(3)
      REAL(dp)                                           :: displacement(3), image_center(3), inverse_cell(3, 3), &
                                                            multiplicity, radius, separation(3)
      REAL(dp), ALLOCATABLE                              :: kernel(:, :, :, :)

      CPASSERT(volume > 0.0_dp .AND. SIZE(gvecs, 1) == 3)
      IF (PRESENT(periodicity)) map%periodicity = periodicity
      IF (PRESENT(grid_shape)) map%grid_shape = grid_shape
      CPASSERT(ALL(map%grid_shape == 0) .OR. ALL(map%grid_shape > 0))
      CPASSERT(ALL(map%periodicity == 0 .OR. map%periodicity == 1))
      map%cell = cell
      map%volume = volume
      ALLOCATE (map%gvecs, SOURCE=gvecs)
      ALLOCATE (aliased(SIZE(gvecs, 2)), SOURCE=.FALSE.)
      IF (ALL(map%grid_shape > 0)) THEN
         DO row = 1, SIZE(gvecs, 2)
            miller = NINT(MATMUL(TRANSPOSE(cell), gvecs(:, row))/twopi)
            aliased(row) = ANY(MOD(map%grid_shape, 2) == 0 .AND. ABS(miller) == map%grid_shape/2)
         END DO
      END IF
      ALLOCATE (map%alias_rows, SOURCE=PACK([(row, row=1, SIZE(aliased))], aliased))
      ALLOCATE (map%centers, SOURCE=centers)
      ALLOCATE (map%fourier(SIZE(centers)), map%pairs(SIZE(centers), SIZE(centers)))
      inverse_cell = cell
      CALL invmat(inverse_cell, info)
      CPASSERT(info == 0)
      rank = 0
      nrank = 1
      IF (PRESENT(comm)) THEN
         rank = comm%mepos
         nrank = comm%num_pe
      END IF
      cartesian = cartesian_reflection_cell(cell)
      owner = -1
      DO a = 1, SIZE(centers)
         ALLOCATE (map%fourier(a)%values(SIZE(gvecs, 2), SIZE(centers(a)%radii), SIZE(centers(a)%l)))
         CALL upf_inverse_coulomb_fourier_kernel(centers(a)%radii, centers(a)%l, centers(a)%indices, &
                                                 gvecs, centers(a)%position, volume, map%fourier(a)%values, &
                                                 cell=cell, grid_shape=map%grid_shape)
         DO b = a, SIZE(centers)
            separation = centers(a)%position - centers(b)%position
            CALL pair_reflection_center(cell, separation, map%periodicity, cartesian, reflection_axes, image_center)
            IF (ANY(reflection_axes)) parity = image_pair_parity(centers(a), centers(b))
            radius = MAXVAL(centers(a)%radii) + MAXVAL(centers(b)%radii)
            CALL residual_image_limits(inverse_cell, separation, radius, map%periodicity, limits)
            DO k = limits(1, 3), limits(2, 3)
               DO j = limits(1, 2), limits(2, 2)
                  DO i = limits(1, 1), limits(2, 1)
                     image = [i, j, k]
                     IF (ANY(reflection_axes)) THEN
                        IF (ANY(reflection_axes .AND. image < image_center)) CYCLE
                     ELSE IF (a == b .AND. negative_image(i, j, k)) THEN
                        CYCLE
                     END IF
                     displacement = separation - MATMUL(cell, REAL(image, dp))
                     IF (SUM(displacement**2) >= radius**2) CYCLE
                     ! Every rank needs the block for reduction, but only intersecting shells need storage.
                     IF (.NOT. ALLOCATED(kernel)) THEN
                        ALLOCATE (kernel(SIZE(centers(a)%radii), SIZE(centers(a)%l), &
                                         SIZE(centers(b)%radii), SIZE(centers(b)%l)))
                        ALLOCATE (map%pairs(a, b)%values, MOLD=kernel)
                        map%pairs(a, b)%values = 0.0_dp
                     END IF
                     owner = MOD(owner + 1, nrank)
                     IF (owner /= rank) CYCLE
                     CALL upf_inverse_coulomb_pair_kernel(centers(a)%radii, centers(a)%l, centers(a)%indices, &
                                                          centers(b)%radii, centers(b)%l, centers(b)%indices, &
                                                          displacement, kernel)
                     mask = SUM(MERGE(reflection_bits, 0, reflection_axes .AND. image /= image_center))
                     IF (mask /= 0) THEN
                        ! The image orbit projects out every odd channel pair.
                        multiplicity = REAL(2**POPCNT(mask), dp)
                        DO h = 1, SIZE(kernel, 4)
                           DO channel = 1, SIZE(kernel, 2)
                              IF (IAND(parity(channel, h), mask) /= 0) CYCLE
                              IF (a /= b) THEN
                                 map%pairs(a, b)%values(:, channel, :, h) = &
                                    map%pairs(a, b)%values(:, channel, :, h) + multiplicity*kernel(:, channel, :, h)
                              ELSE
                                 ! Average opposite orientations to retain exact reciprocity.
                                 DO r = 1, SIZE(kernel, 3)
                                    map%pairs(a, b)%values(:, channel, r, h) = &
                                       map%pairs(a, b)%values(:, channel, r, h) + 0.5_dp*multiplicity* &
                                       (kernel(:, channel, r, h) + kernel(r, h, :, channel))
                                 END DO
                              END IF
                           END DO
                        END DO
                        CYCLE
                     END IF
                     map%pairs(a, b)%values(:, :, :, :) = map%pairs(a, b)%values + kernel
                     IF (a == b .AND. ANY(image /= 0)) THEN
                        ! Reversing the lattice image exchanges the two shell distributions.
                        DO h = 1, SIZE(kernel, 2)
                           DO r = 1, SIZE(kernel, 1)
                              map%pairs(a, b)%values(:, :, r, h) = map%pairs(a, b)%values(:, :, r, h) + kernel(r, h, :, :)
                           END DO
                        END DO
                     END IF
                  END DO
               END DO
            END DO
            IF (.NOT. ALLOCATED(kernel)) CYCLE
            IF (nrank > 1) CALL comm%sum(map%pairs(a, b)%values)
            IF (a /= b) THEN
               ALLOCATE (map%pairs(b, a)%values(SIZE(kernel, 3), SIZE(kernel, 4), SIZE(kernel, 1), SIZE(kernel, 2)))
               DO h = 1, SIZE(kernel, 2)
                  DO r = 1, SIZE(kernel, 1)
                     map%pairs(b, a)%values(:, :, r, h) = map%pairs(a, b)%values(r, h, :, :)
                  END DO
               END DO
            END IF
            DEALLOCATE (kernel)
         END DO
      END DO
   END SUBROUTINE upf_inverse_coulomb_prepare

! **************************************************************************************************
!> \brief Contract residual-pair geometry derivatives once at fixed radial shells and density.
!>        Keep every overlapping image and both orientations through pair reciprocity. Only
!>        center and cell covectors are retained, rather than a kernel for every motion.
!>        Density, measure and shape responses remain in the full Coulomb response.
!> \param centers Current replicated signed densities on the map's unchanged radial shells.
!> \param map Immutable geometry and periodic image convention.
!> \param center_gradient Derivative with respect to each Cartesian center position.
!> \param cell_gradient Derivative with respect to lattice vectors at fixed center positions.
!> \param comm Optional process group sharing the same replicated centers and map.
!>        Distributes Cartesian pair/image derivatives and reduces only the resulting covectors.
!> \param cell_derivative Whether to compute the cell gradient. Defaults to true.
! **************************************************************************************************
   SUBROUTINE upf_inverse_coulomb_pair_geometry(centers, map, center_gradient, cell_gradient, comm, cell_derivative)
      TYPE(upf_inverse_coulomb_center_type), INTENT(IN)  :: centers(:)
      TYPE(upf_inverse_coulomb_map_type), INTENT(IN)     :: map
      REAL(dp), INTENT(OUT)                              :: center_gradient(:, :), &
                                                            cell_gradient(3, 3)
      TYPE(mp_comm_type), INTENT(IN), OPTIONAL           :: comm
      LOGICAL, INTENT(IN), OPTIONAL                      :: cell_derivative

      INTEGER                                            :: a, axis, b, channel, column, h, i, &
                                                            image(3), info, j, k, limits(2, 3), &
                                                            mask, nr, nrank, owner, pair_parity, &
                                                            rank
      INTEGER, ALLOCATABLE                               :: parity(:, :)
      LOGICAL                                            :: cartesian, do_cell, reflection_axes(3)
      REAL(dp)                                           :: center_value, delta(3), displacement(3), &
                                                            image_center(3), inverse_cell(3, 3), multiplicity, &
                                                            radius, separation(3), value, values(3)
      REAL(dp), ALLOCATABLE                              :: dkernel(:, :, :, :), kernel(:, :, :, :), &
                                                            moments(:)
      TYPE(center_work_type), ALLOCATABLE                :: work(:)

      CALL check_map(map, centers, map%cell, map%gvecs, map%volume, map%periodicity)
      CPASSERT(ALL(SHAPE(center_gradient) == [3, SIZE(centers)]))
      center_gradient = 0.0_dp
      cell_gradient = 0.0_dp
      do_cell = .TRUE.
      IF (PRESENT(cell_derivative)) do_cell = cell_derivative
      IF (.NOT. do_cell .AND. SIZE(centers) <= 1) RETURN
      rank = 0
      nrank = 1
      IF (PRESENT(comm)) THEN
         rank = comm%mepos
         nrank = comm%num_pe
      END IF
      cartesian = cartesian_reflection_cell(map%cell)
      owner = -1
      inverse_cell = map%cell
      CALL invmat(inverse_cell, info)
      CPASSERT(info == 0)
      ALLOCATE (work(SIZE(centers)))
      DO a = 1, SIZE(centers)
         ALLOCATE (work(a)%compensation(SIZE(centers(a)%radii), SIZE(centers(a)%l)), &
                   work(a)%residual(SIZE(centers(a)%radii), SIZE(centers(a)%l)), &
                   work(a)%q(SIZE(centers(a)%radii), SIZE(centers(a)%l)), moments(SIZE(centers(a)%l)))
         CALL upf_inverse_coulomb_compensate(centers(a)%radii, centers(a)%weights, centers(a)%l, &
                                             centers(a)%density, centers(a)%shapes, work(a)%compensation, &
                                             work(a)%residual, moments)
         work(a)%q(:, :) = SPREAD(centers(a)%weights, 2, SIZE(centers(a)%l))*work(a)%residual
         DEALLOCATE (moments)
      END DO
      DO a = 1, SIZE(centers)
         DO b = a, SIZE(centers)
            IF (.NOT. ALLOCATED(map%pairs(a, b)%values)) CYCLE
            ! Moving a center and its own images together leaves their separation fixed.
            IF (a == b .AND. .NOT. do_cell) CYCLE
            separation = centers(a)%position - centers(b)%position
            CALL pair_reflection_center(map%cell, separation, map%periodicity, cartesian, reflection_axes, image_center)
            IF (ANY(reflection_axes)) parity = image_pair_parity(centers(a), centers(b))
            radius = MAXVAL(centers(a)%radii) + MAXVAL(centers(b)%radii)
            CALL residual_image_limits(inverse_cell, separation, radius, map%periodicity, limits)
            DO k = limits(1, 3), limits(2, 3)
               DO j = limits(1, 2), limits(2, 2)
                  DO i = limits(1, 1), limits(2, 1)
                     image = [i, j, k]
                     ! The two orientations of a self image give the same contracted cell derivative.
                     IF (a == b .AND. ALL(image == 0)) CYCLE
                     IF (ANY(reflection_axes)) THEN
                        IF (ANY(reflection_axes .AND. image < image_center)) CYCLE
                     ELSE IF (a == b .AND. negative_image(i, j, k)) THEN
                        CYCLE
                     END IF
                     displacement = separation - MATMUL(map%cell, REAL(image, dp))
                     IF (SUM(displacement**2) >= radius**2) CYCLE
                     mask = SUM(MERGE(reflection_bits, 0, reflection_axes .AND. image /= image_center))
                     DO axis = 1, 3
                        owner = MOD(owner + 1, nrank)
                        IF (owner /= rank) CYCLE
                        IF (.NOT. ALLOCATED(kernel)) THEN
                           ALLOCATE (kernel(SIZE(centers(a)%radii), SIZE(centers(a)%l), &
                                            SIZE(centers(b)%radii), SIZE(centers(b)%l)))
                           ALLOCATE (dkernel, MOLD=kernel)
                        END IF
                        delta = 0.0_dp
                        delta(axis) = 1.0_dp
                        CALL upf_inverse_coulomb_pair_kernel(centers(a)%radii, centers(a)%l, centers(a)%indices, &
                                                             centers(b)%radii, centers(b)%l, centers(b)%indices, &
                                                             displacement, kernel, dkernel, 0.0_dp*centers(a)%radii, &
                                                             0.0_dp*centers(b)%radii, delta)
                        IF (mask /= 0) THEN
                           ! d/dD_axis and the lattice index in column each add one reflection sign.
                           multiplicity = REAL(2**POPCNT(mask), dp)
                           IF (a == b) multiplicity = 0.5_dp*multiplicity
                           center_value = 0.0_dp
                           values = 0.0_dp
                           DO h = 1, SIZE(centers(b)%l)
                              DO channel = 1, SIZE(centers(a)%l)
                                 pair_parity = IEOR(parity(channel, h), reflection_bits(axis))
                                 value = DOT_PRODUCT(work(a)%q(:, channel), &
                                                     MATMUL(dkernel(:, channel, :, h), work(b)%q(:, h)))
                                 IF (IAND(pair_parity, mask) == 0) center_value = center_value + value
                                 IF (do_cell) THEN
                                    DO column = 1, 3
                                       IF (IAND(IEOR(pair_parity, reflection_bits(column)), mask) /= 0) CYCLE
                                       values(column) = values(column) + value
                                    END DO
                                 END IF
                              END DO
                           END DO
                           IF (a /= b) THEN
                              center_gradient(axis, a) = center_gradient(axis, a) + multiplicity*center_value
                              center_gradient(axis, b) = center_gradient(axis, b) - multiplicity*center_value
                           END IF
                           ! Reflected image indices are centered on the relative fractional coordinate.
                           IF (do_cell) cell_gradient(axis, :) = cell_gradient(axis, :) - multiplicity* &
                                                               (image_center*center_value + (REAL(image, dp) - image_center)*values)
                           CYCLE
                        END IF
                        value = 0.0_dp
                        DO channel = 1, SIZE(centers(a)%l)
                           DO nr = 1, SIZE(centers(a)%radii)
                              value = value + work(a)%q(nr, channel)* &
                                      SUM(dkernel(nr, channel, :, :)*work(b)%q)
                           END DO
                        END DO
                        IF (a /= b) THEN
                           center_gradient(axis, a) = center_gradient(axis, a) + value
                           center_gradient(axis, b) = center_gradient(axis, b) - value
                        END IF
                        IF (do_cell) cell_gradient(axis, :) = cell_gradient(axis, :) - value*REAL(image, dp)
                     END DO
                  END DO
               END DO
            END DO
            IF (ALLOCATED(kernel)) DEALLOCATE (kernel, dkernel)
         END DO
      END DO
      IF (PRESENT(comm)) THEN
         CALL comm%sum(center_gradient)
         IF (do_cell) CALL comm%sum(cell_gradient)
      END IF
   END SUBROUTINE upf_inverse_coulomb_pair_geometry

! **************************************************************************************************
!> \brief Fourier coefficients of the complete signed center density, before compensation.
!>        This is the finite-grid projection for density observables. It includes every supplied
!>        harmonic and the zero mode, without ions, frozen cores, or a Coulomb kernel.
!> \param centers Same radial/harmonic geometry as map, with the requested spin density.
!> \param cell Physical lattice vectors in columns.
!> \param gvecs Local reciprocal vectors of the observable's Fourier grid.
!> \param volume Cell volume defining the Fourier normalization.
!> \param density Overwritten local Fourier coefficients, without mesh multiplicity.
!> \param map Immutable geometry map shared with the Coulomb functional.
!> \param response Optional joint density derivative. The original map is not differentiated in place.
!> \param center_response Radial-node, measure, density and center motion, required with response.
!> \param gvecs_response Reciprocal-vector motion, required with response.
!> \param volume_response Cell-volume derivative, required with response.
!> \param cell_response Optional lattice derivative for the prepared real-grid boundary convention.
! **************************************************************************************************
   SUBROUTINE upf_inverse_coulomb_density(centers, cell, gvecs, volume, density, map, response, &
                                          center_response, gvecs_response, volume_response, cell_response)
      TYPE(upf_inverse_coulomb_center_type), INTENT(IN)  :: centers(:)
      REAL(dp), INTENT(IN)                               :: cell(3, 3), gvecs(:, :), volume
      COMPLEX(dp), INTENT(OUT)                           :: density(:)
      TYPE(upf_inverse_coulomb_map_type), INTENT(IN)     :: map
      COMPLEX(dp), INTENT(OUT), OPTIONAL                 :: response(:)
      TYPE(upf_inverse_coulomb_center_type), &
         INTENT(IN), OPTIONAL                            :: center_response(:)
      REAL(dp), INTENT(IN), OPTIONAL                     :: gvecs_response(:, :), volume_response, &
                                                            cell_response(3, 3)

      COMPLEX(dp), ALLOCATABLE                           :: dfourier(:, :, :), fourier(:, :, :), phase_response(:)
      INTEGER                                            :: a, channel, ir
      LOGICAL                                            :: derive_kernel, fixed_reciprocal, translation_only, valid_shape
      REAL(dp)                                           :: weight, weight_response

      CALL check_map(map, centers, cell, gvecs, volume, map%periodicity)
      CPASSERT(SIZE(density) == SIZE(gvecs, 2))
      density = CMPLX(0.0_dp, 0.0_dp, dp)
      fixed_reciprocal = .FALSE.
      IF (PRESENT(response)) THEN
         CPASSERT(PRESENT(center_response))
         CPASSERT(PRESENT(gvecs_response))
         CPASSERT(PRESENT(volume_response))
         CPASSERT(SIZE(response) == SIZE(density))
         CPASSERT(SIZE(center_response) == SIZE(centers))
         CPASSERT(ALL(SHAPE(gvecs_response) == SHAPE(gvecs)))
         response = CMPLX(0.0_dp, 0.0_dp, dp)
         fixed_reciprocal = ALL(gvecs_response == 0.0_dp) .AND. volume_response == 0.0_dp
         IF (PRESENT(cell_response)) fixed_reciprocal = fixed_reciprocal .AND. ALL(cell_response == 0.0_dp)
      END IF
      DO a = 1, SIZE(centers)
         CALL check_radial(centers(a)%radii, centers(a)%weights, centers(a)%l, centers(a)%density)
         IF (PRESENT(response)) THEN
            CALL check_radial(centers(a)%radii, centers(a)%weights, centers(a)%l, center_response(a)%density)
            valid_shape = SIZE(center_response(a)%weights) == SIZE(centers(a)%weights)
            CPASSERT(valid_shape)
            valid_shape = SIZE(center_response(a)%radii) == SIZE(centers(a)%radii)
            CPASSERT(valid_shape)
            derive_kernel = .TRUE.
            translation_only = .FALSE.
            IF (fixed_reciprocal .AND. ALL(center_response(a)%radii == 0.0_dp)) THEN
               derive_kernel = ANY(center_response(a)%position /= 0.0_dp)
               translation_only = derive_kernel
            END IF
            IF (translation_only) THEN
               CALL fourier_translation_response(centers(a), center_response(a)%position, map, phase_response, dfourier)
            ELSE IF (derive_kernel) THEN
               ALLOCATE (fourier, MOLD=map%fourier(a)%values)
               ALLOCATE (dfourier, MOLD=fourier)
               CALL upf_inverse_coulomb_fourier_kernel(centers(a)%radii, centers(a)%l, centers(a)%indices, &
                                                       gvecs, centers(a)%position, volume, fourier, dfourier, &
                                                       center_response(a)%radii, gvecs_response, &
                                                       center_response(a)%position, volume_response, &
                                                       cell=cell, grid_shape=map%grid_shape, cell_response=cell_response)
               DEALLOCATE (fourier)
            END IF
         END IF
         DO channel = 1, SIZE(centers(a)%l)
            DO ir = 1, SIZE(centers(a)%radii)
               density(:) = density + map%fourier(a)%values(:, ir, channel)* &
                            centers(a)%weights(ir)*centers(a)%density(ir, channel)
               IF (.NOT. PRESENT(response)) CYCLE
               weight = centers(a)%weights(ir)*centers(a)%density(ir, channel)
               weight_response = center_response(a)%weights(ir)*centers(a)%density(ir, channel) + &
                                 centers(a)%weights(ir)*center_response(a)%density(ir, channel)
               response(:) = response + map%fourier(a)%values(:, ir, channel)*weight_response
               IF (translation_only) THEN
                  response(:) = response + phase_response*map%fourier(a)%values(:, ir, channel)*weight
                  response(map%alias_rows) = response(map%alias_rows) + dfourier(:, ir, channel)*weight
               ELSE IF (derive_kernel) THEN
                  response(:) = response + dfourier(:, ir, channel)*weight
               END IF
            END DO
         END DO
         IF (ALLOCATED(phase_response)) DEALLOCATE (phase_response)
         IF (ALLOCATED(dfourier)) DEALLOCATE (dfourier)
      END DO
   END SUBROUTINE upf_inverse_coulomb_density

! **************************************************************************************************
!> \brief Derivative of the cached Fourier kernel for a pure center translation.
!> \param center Center and radial shell geometry of the cached coefficients.
!> \param displacement Center displacement at fixed lattice and radial nodes.
!> \param map Prepared Fourier geometry.
!> \param phase_response Phase derivative multiplier, zero on aliased rows.
!> \param alias_response Full kernel derivative on aliased rows.
! **************************************************************************************************
   SUBROUTINE fourier_translation_response(center, displacement, map, phase_response, alias_response)
      TYPE(upf_inverse_coulomb_center_type), INTENT(IN)  :: center
      REAL(dp), INTENT(IN)                               :: displacement(3)
      TYPE(upf_inverse_coulomb_map_type), INTENT(IN)     :: map
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: phase_response(:), alias_response(:, :, :)

      COMPLEX(dp), ALLOCATABLE                           :: fourier(:, :, :)
      INTEGER                                            :: nalias

      nalias = SIZE(map%alias_rows)
      ALLOCATE (phase_response, SOURCE=CMPLX(0.0_dp, -MATMUL(displacement, map%gvecs), dp))
      ALLOCATE (alias_response(nalias, SIZE(center%radii), SIZE(center%l)))
      IF (nalias == 0) RETURN
      ! Each Nyquist row averages phases of several distinct physical vectors.
      phase_response(map%alias_rows) = CMPLX(0.0_dp, 0.0_dp, dp)
      ALLOCATE (fourier, MOLD=alias_response)
      CALL upf_inverse_coulomb_fourier_kernel(center%radii, center%l, center%indices, &
                                              map%gvecs(:, map%alias_rows), center%position, map%volume, &
                                              fourier, alias_response, 0.0_dp*center%radii, &
                                              0.0_dp*map%gvecs(:, map%alias_rows), displacement, &
                                              0.0_dp, cell=map%cell, grid_shape=map%grid_shape)
   END SUBROUTINE fourier_translation_response

! **************************************************************************************************
!> \brief Transpose of the complete signed center-density Fourier projection for a local potential.
!>        The caller supplies the Fourier covector, including volume and multiplicity exactly once.
!>        No Coulomb kernel, compensation approximation, ionic or frozen-core term is inserted.
!> \param centers Same radial/harmonic geometry as map. The density values do not enter the transpose.
!> \param cell Physical lattice vectors in columns.
!> \param gvecs Local reciprocal vectors, in the same order as potential.
!> \param volume Cell volume defining the density Fourier coefficients.
!> \param potential Local Fourier covector for Re(sum(conjg(density)*potential)).
!> \param potentials Radial-measure adjoints, so dE=sum(weights*dn*potentials%values).
!> \param map Immutable geometry map shared with the density projection.
!> \param comm Optional communicator distributing Fourier rows. Radial potentials are replicated.
! **************************************************************************************************
   SUBROUTINE upf_inverse_coulomb_density_adjoint(centers, cell, gvecs, volume, potential, potentials, map, comm)
      TYPE(upf_inverse_coulomb_center_type), INTENT(IN)  :: centers(:)
      REAL(dp), INTENT(IN)                               :: cell(3, 3), gvecs(:, :), volume
      COMPLEX(dp), INTENT(IN)                            :: potential(:)
      TYPE(upf_inverse_coulomb_potential_type), &
         ALLOCATABLE, INTENT(OUT)                        :: potentials(:)
      TYPE(upf_inverse_coulomb_map_type), INTENT(IN)     :: map
      TYPE(mp_comm_type), INTENT(IN), OPTIONAL           :: comm

      INTEGER                                            :: a

      CALL check_map(map, centers, cell, gvecs, volume, map%periodicity)
      CPASSERT(SIZE(potential) == SIZE(gvecs, 2))
      ALLOCATE (potentials(SIZE(centers)))
      DO a = 1, SIZE(centers)
         ALLOCATE (potentials(a)%values(SIZE(centers(a)%radii), SIZE(centers(a)%l)))
         CALL fourier_density_adjoint(map%fourier(a)%values, potential, potentials(a)%values)
         IF (PRESENT(comm)) CALL comm%sum(potentials(a)%values)
      END DO
   END SUBROUTINE upf_inverse_coulomb_density_adjoint

! **************************************************************************************************
!> \brief Local Fourier contraction, shared by arbitrary potentials and the complete Coulomb adjoint.
!> \param fourier Density map without radial weights.
!> \param potential Fourier covector, with the caller's volume and multiplicity convention.
!> \param values Overwritten radial-measure adjoint.
! **************************************************************************************************
   SUBROUTINE fourier_density_adjoint(fourier, potential, values)
      COMPLEX(dp), INTENT(IN)                            :: fourier(:, :, :), potential(:)
      REAL(dp), INTENT(OUT)                              :: values(:, :)

      INTEGER                                            :: channel, ir
      LOGICAL                                            :: valid_shape

      CPASSERT(SIZE(fourier, 1) == SIZE(potential))
      valid_shape = SIZE(values, 1) == SIZE(fourier, 2) .AND. SIZE(values, 2) == SIZE(fourier, 3)
      CPASSERT(valid_shape)
      DO channel = 1, SIZE(values, 2)
         DO ir = 1, SIZE(values, 1)
            values(ir, channel) = REAL(SUM(CONJG(fourier(:, ir, channel))*potential), dp)
         END DO
      END DO
   END SUBROUTINE fourier_density_adjoint

! **************************************************************************************************
!> \brief Complete periodic Coulomb functional for a Fourier background and signed centered fields.
!>        Moment compensation uses the actual radial measure. Every overlapping residual image is
!>        retained, including self images. Separated residuals have exactly zero exterior potential.
!>        Gaussian hard/soft, reconstruction and original Q contributions belong in the supplied
!>        centered density. Ionic and Gaussian background charges belong in background.
!> \param centers Replicated harmonic densities, radial measures and compensation shapes.
!> \param cell Physical periodic cell, with lattice vectors in columns.
!> \param gvecs Local physical reciprocal vectors, in the same order as background.
!> \param volume Physical cell volume used in the Fourier coefficients.
!> \param background Fourier background before adding compensation charges.
!> \param poisson_weights Native Coulomb kernel times volume and Fourier multiplicity. The caller
!>        supplies the zero-mode convention. These weights are nonnegative.
!> \param energy Full Coulomb energy, excluding any separately defined ionic self energies.
!> \param potentials Optional radial-measure density adjoints, one for each supplied center.
!> \param background_adjoint Local Fourier covector, including the residual/background interaction.
!> \param comm Optional communicator distributing Fourier rows. Centered data are replicated.
!> \param response Optional joint energy derivative of this same finite representation.
!> \param center_response Radial, measure, density, shape and center derivatives, required with response.
!> \param cell_response Lattice derivative, required with response.
!> \param gvecs_response Physical reciprocal-vector derivative, required with response.
!> \param volume_response Cell-volume derivative, required with response.
!> \param background_response Background Fourier derivative, required with response.
!> \param poisson_weights_response Complete native Coulomb-weight derivative, required with response.
!> \param map Optional immutable maps for this exact geometry. Joint responses use the original
!>        differentiated kernels, so freezing a map never freezes the geometry derivative.
!> \param periodicity Optional lattice directions to repeat. Defaults to the map's periodicity
!>        when a map is supplied, otherwise to three periodic directions. The Fourier weights
!>        and their response must use the same electrostatic boundary conditions.
!> \param pair_geometry_response Optional previously contracted explicit residual-kernel response
!>        at this density and geometry. Used only with map and fixed radial shells. Radial motion
!>        retains the full original derivative. Density, measure and shape derivatives are unchanged.
! **************************************************************************************************
   SUBROUTINE upf_inverse_coulomb_periodic(centers, cell, gvecs, volume, background, poisson_weights, &
                                           energy, potentials, background_adjoint, comm, response, center_response, &
                                           cell_response, gvecs_response, volume_response, background_response, &
                                           poisson_weights_response, map, periodicity, pair_geometry_response)
      TYPE(upf_inverse_coulomb_center_type), INTENT(IN)  :: centers(:)
      REAL(dp), INTENT(IN)                               :: cell(3, 3), gvecs(:, :), volume
      COMPLEX(dp), INTENT(IN)                            :: background(:)
      REAL(dp), INTENT(IN)                               :: poisson_weights(:)
      REAL(dp), INTENT(OUT)                              :: energy
      TYPE(upf_inverse_coulomb_potential_type), &
         ALLOCATABLE, INTENT(OUT), OPTIONAL              :: potentials(:)
      COMPLEX(dp), INTENT(OUT)                           :: background_adjoint(:)
      TYPE(mp_comm_type), INTENT(IN), OPTIONAL           :: comm
      REAL(dp), INTENT(OUT), OPTIONAL                    :: response
      TYPE(upf_inverse_coulomb_center_type), &
         INTENT(IN), OPTIONAL                            :: center_response(:)
      REAL(dp), INTENT(IN), OPTIONAL                     :: cell_response(3, 3), &
                                                            gvecs_response(:, :), volume_response
      COMPLEX(dp), INTENT(IN), OPTIONAL                  :: background_response(:)
      REAL(dp), INTENT(IN), OPTIONAL                     :: poisson_weights_response(:)
      TYPE(upf_inverse_coulomb_map_type), INTENT(IN), TARGET, &
         OPTIONAL                                        :: map
      INTEGER, INTENT(IN), OPTIONAL                      :: periodicity(3)
      REAL(dp), INTENT(IN), OPTIONAL                     :: pair_geometry_response

      COMPLEX(dp), ALLOCATABLE                           :: dfourier(:, :, :), dg(:), dh(:), e(:), &
                                                            h(:), phase_response(:)
      COMPLEX(dp), CONTIGUOUS, POINTER                   :: kernel_values(:, :, :)
      INTEGER                                            :: a, b, channel, fourier_shape(3), i, &
                                                            image(3), info, j, limits(2, 3), na, &
                                                            nb, ng, nh, nr, REPEAT(3)
      LOGICAL                                            :: derive, derive_kernel, fixed_reciprocal, &
                                                            reuse_fourier, reuse_pair_geometry, translation_only
      REAL(dp)                                           :: ddisplacement(3), displacement(3), &
                                                            dmoment, dnorm, dpair, dpw, epair, &
                                                            inverse_cell(3, 3), moment, norm, &
                                                            radius
      REAL(dp), ALLOCATABLE :: adjoint(:, :), dc(:, :), de(:, :), df(:), dkernel(:, :, :, :), &
         f(:), kernel(:, :, :, :), moments(:), ve(:, :), vh(:, :)
      TYPE(center_work_type), ALLOCATABLE                :: work(:)
      TYPE(fourier_map_type), ALLOCATABLE, TARGET         :: scratch_fourier(:)
      TYPE(fourier_map_type), POINTER                     :: fourier(:)

      ng = SIZE(background)
      fourier_shape = 0
      IF (PRESENT(map)) fourier_shape = map%grid_shape
      repeat = [1, 1, 1]
      IF (PRESENT(map)) repeat = map%periodicity
      IF (PRESENT(periodicity)) repeat = periodicity
      CPASSERT(ALL(repeat == 0 .OR. repeat == 1))
      CPASSERT(ALL(SHAPE(gvecs) == [3, ng]))
      CPASSERT(SIZE(poisson_weights) == ng)
      CPASSERT(SIZE(background_adjoint) == ng)
      CPASSERT(ALL(ieee_is_finite(poisson_weights)))
      ! Truncated Coulomb kernels can have signed reciprocal coefficients.
      inverse_cell = cell
      CALL invmat(inverse_cell, info)
      CPASSERT(info == 0)
      derive = PRESENT(response)
      fixed_reciprocal = .FALSE.
      IF (derive) THEN
         CPASSERT(PRESENT(center_response) .AND. PRESENT(cell_response))
         CPASSERT(PRESENT(gvecs_response) .AND. PRESENT(volume_response))
         CPASSERT(PRESENT(background_response))
         CPASSERT(PRESENT(poisson_weights_response))
         CPASSERT(SIZE(center_response) == SIZE(centers))
         CPASSERT(ALL(SHAPE(gvecs_response) == SHAPE(gvecs)))
         CPASSERT(SIZE(background_response) == ng)
         CPASSERT(SIZE(poisson_weights_response) == ng)
         fixed_reciprocal = PRESENT(map) .AND. ALL(cell_response == 0.0_dp) .AND. &
                            ALL(gvecs_response == 0.0_dp) .AND. volume_response == 0.0_dp
      END IF
      IF (PRESENT(map)) CALL check_map(map, centers, cell, gvecs, volume, repeat)
      ! The reciprocal coefficients are immutable at fixed geometry. Borrow them
      ! for both the forward density and its adjoint without an SCF-sized copy.
      reuse_fourier = PRESENT(map) .AND. .NOT. derive
      IF (reuse_fourier) THEN
         fourier => map%fourier
      ELSE
         ALLOCATE (scratch_fourier(SIZE(centers)))
         fourier => scratch_fourier
      END IF
      reuse_pair_geometry = .FALSE.
      IF (derive .AND. PRESENT(map) .AND. PRESENT(pair_geometry_response)) THEN
         reuse_pair_geometry = .TRUE.
         DO a = 1, SIZE(centers)
            IF (ANY(center_response(a)%radii /= 0.0_dp)) reuse_pair_geometry = .FALSE.
         END DO
      END IF
      ALLOCATE (work(SIZE(centers)), h(ng), e(ng), dh(ng), dg(ng))
      h(:) = background
      e(:) = CMPLX(0.0_dp, 0.0_dp, dp)
      dh(:) = CMPLX(0.0_dp, 0.0_dp, dp)
      dg(:) = CMPLX(0.0_dp, 0.0_dp, dp)
      IF (derive) dh(:) = background_response
      DO a = 1, SIZE(centers)
         nr = SIZE(centers(a)%radii)
         nh = SIZE(centers(a)%l)
         ALLOCATE (work(a)%compensation(nr, nh), work(a)%residual(nr, nh), &
                   work(a)%q(nr, nh), work(a)%dq(nr, nh), work(a)%potential(nr, nh), &
                   moments(nh))
         derive_kernel = derive
         IF (fixed_reciprocal) THEN
            ! Density and compensation can respond even at a stationary center.
            derive_kernel = ANY(center_response(a)%position /= 0.0_dp) .OR. &
                            ANY(center_response(a)%radii /= 0.0_dp)
         END IF
         ! A translation changes only the phase of each unaliased coefficient.
         translation_only = .FALSE.
         IF (derive_kernel .AND. fixed_reciprocal) THEN
            translation_only = ALL(center_response(a)%radii == 0.0_dp)
         END IF
         IF ((derive_kernel .AND. .NOT. translation_only) .OR. (.NOT. derive .AND. .NOT. reuse_fourier)) THEN
            ALLOCATE (fourier(a)%values(ng, nr, nh))
         END IF
         CALL upf_inverse_coulomb_compensate(centers(a)%radii, centers(a)%weights, centers(a)%l, &
                                             centers(a)%density, centers(a)%shapes, work(a)%compensation, &
                                             work(a)%residual, moments)
         work(a)%q(:, :) = SPREAD(centers(a)%weights, 2, nh)*work(a)%residual
         work(a)%potential(:, :) = 0.0_dp
         work(a)%dq(:, :) = 0.0_dp
         IF (derive) THEN
            CALL check_radial(centers(a)%radii, centers(a)%weights, centers(a)%l, center_response(a)%density)
            CPASSERT(SIZE(center_response(a)%radii) == nr)
            CPASSERT(SIZE(center_response(a)%weights) == nr)
            CPASSERT(ALL(SHAPE(center_response(a)%shapes) == [nr, nh]))
            ALLOCATE (dc(nr, nh), de(nr, nh), f(nr), df(nr))
            DO channel = 1, nh
               f(:) = centers(a)%radii**centers(a)%l(channel)
               df(:) = f*centers(a)%l(channel)*center_response(a)%radii/centers(a)%radii
               norm = SUM(centers(a)%weights*f*centers(a)%shapes(:, channel))
               dnorm = SUM((center_response(a)%weights*f + centers(a)%weights*df)*centers(a)%shapes(:, channel) + &
                           centers(a)%weights*f*center_response(a)%shapes(:, channel))
               moment = moments(channel)
               dmoment = SUM((center_response(a)%weights*f + centers(a)%weights*df)*centers(a)%density(:, channel) + &
                             centers(a)%weights*f*center_response(a)%density(:, channel))
               dc(:, channel) = centers(a)%shapes(:, channel)*dmoment/norm + moment/norm* &
                                (center_response(a)%shapes(:, channel) - centers(a)%shapes(:, channel)*dnorm/norm)
            END DO
            de(:, :) = center_response(a)%density - dc
            work(a)%dq(:, :) = SPREAD(center_response(a)%weights, 2, nh)*work(a)%residual + &
                               SPREAD(centers(a)%weights, 2, nh)*de
            IF (translation_only) THEN
               CALL fourier_translation_response(centers(a), center_response(a)%position, map, phase_response, dfourier)
            ELSE IF (derive_kernel) THEN
               ALLOCATE (dfourier(ng, nr, nh))
               CALL upf_inverse_coulomb_fourier_kernel(centers(a)%radii, centers(a)%l, centers(a)%indices, &
                                                       gvecs, centers(a)%position, volume, fourier(a)%values, dfourier, &
                                                       center_response(a)%radii, gvecs_response, &
                                                       center_response(a)%position, volume_response, &
                                                       cell=cell, grid_shape=fourier_shape, cell_response=cell_response)
            END IF
         ELSE IF (.NOT. reuse_fourier) THEN
            CALL upf_inverse_coulomb_fourier_kernel(centers(a)%radii, centers(a)%l, centers(a)%indices, &
                                                    gvecs, centers(a)%position, volume, fourier(a)%values)
         END IF
         IF (ALLOCATED(fourier(a)%values)) THEN
            kernel_values => fourier(a)%values
         ELSE
            kernel_values => map%fourier(a)%values
         END IF
         DO channel = 1, nh
            DO i = 1, nr
               h(:) = h + kernel_values(:, i, channel)*centers(a)%weights(i)*work(a)%compensation(i, channel)
               e(:) = e + kernel_values(:, i, channel)*work(a)%q(i, channel)
               IF (.NOT. derive) CYCLE
               IF (translation_only) THEN
                  dh(:) = dh + phase_response*kernel_values(:, i, channel)*centers(a)%weights(i)*work(a)%compensation(i, channel)
                  dg(:) = dg + phase_response*kernel_values(:, i, channel)*work(a)%q(i, channel)
                  dh(map%alias_rows) = dh(map%alias_rows) + &
                                       dfourier(:, i, channel)*centers(a)%weights(i)*work(a)%compensation(i, channel)
                  dg(map%alias_rows) = dg(map%alias_rows) + dfourier(:, i, channel)*work(a)%q(i, channel)
               ELSE IF (derive_kernel) THEN
                  dh(:) = dh + dfourier(:, i, channel)*centers(a)%weights(i)*work(a)%compensation(i, channel)
                  dg(:) = dg + dfourier(:, i, channel)*work(a)%q(i, channel)
               END IF
               dh(:) = dh + kernel_values(:, i, channel)*(center_response(a)%weights(i)*work(a)%compensation(i, channel) + &
                                                         centers(a)%weights(i)*dc(i, channel))
               dg(:) = dg + kernel_values(:, i, channel)*work(a)%dq(i, channel)
            END DO
         END DO
         DEALLOCATE (moments)
         IF (derive) DEALLOCATE (dc, de, f, df)
         IF (ALLOCATED(phase_response)) DEALLOCATE (phase_response)
         IF (ALLOCATED(dfourier)) DEALLOCATE (dfourier)
         IF (.NOT. PRESENT(potentials) .AND. .NOT. reuse_fourier) THEN
            IF (ALLOCATED(fourier(a)%values)) DEALLOCATE (fourier(a)%values)
         END IF
      END DO
      ! Preserve small geometry differences when summing the reciprocal energy terms.
      energy = accurate_sum(poisson_weights*(0.5_dp*ABS(h)**2 + REAL(CONJG(h)*e, dp)))
      dpw = 0.0_dp
      IF (derive) THEN
         dpw = SUM(poisson_weights_response*(0.5_dp*ABS(h)**2 + REAL(CONJG(h)*e, dp)) + &
                   poisson_weights*REAL(CONJG(dh)*(h + e) + CONJG(h)*dg, dp))
      END IF
      IF (PRESENT(comm)) THEN
         CALL comm%sum(energy)
         IF (derive) CALL comm%sum(dpw)
      END IF
      background_adjoint = poisson_weights*(h + e)
      epair = 0.0_dp
      dpair = 0.0_dp
      DO a = 1, SIZE(centers)
         na = SIZE(centers(a)%radii)
         DO b = 1, SIZE(centers)
            IF (PRESENT(map)) THEN
               IF (.NOT. ALLOCATED(map%pairs(a, b)%values)) CYCLE
               ! The immutable map also supplies the primal and density response during motion.
               ! Only the explicit shell-kernel derivative needs to be integrated again.
               DO channel = 1, SIZE(centers(a)%l)
                  DO nr = 1, na
                     norm = SUM(map%pairs(a, b)%values(nr, channel, :, :)*work(b)%q)
                     work(a)%potential(nr, channel) = work(a)%potential(nr, channel) + norm
                     epair = epair + 0.5_dp*work(a)%q(nr, channel)*norm
                     IF (.NOT. derive) CYCLE
                     dpair = dpair + 0.5_dp*(work(a)%dq(nr, channel)*norm + work(a)%q(nr, channel)* &
                                             SUM(map%pairs(a, b)%values(nr, channel, :, :)*work(b)%dq))
                  END DO
               END DO
               IF (.NOT. derive .OR. reuse_pair_geometry) CYCLE
               ! The map already supplies both density covectors. Transposed pairs have
               ! the same joint geometry contraction, so integrate that derivative once.
               IF (b < a) CYCLE
            END IF
            nb = SIZE(centers(b)%radii)
            radius = MAXVAL(centers(a)%radii) + MAXVAL(centers(b)%radii)
            CALL residual_image_limits(inverse_cell, centers(a)%position - centers(b)%position, &
                                       radius, repeat, limits)
            ALLOCATE (kernel(na, SIZE(centers(a)%l), nb, SIZE(centers(b)%l)))
            IF (derive) ALLOCATE (dkernel, MOLD=kernel)
            DO j = limits(1, 3), limits(2, 3)
               DO i = limits(1, 2), limits(2, 2)
                  DO info = limits(1, 1), limits(2, 1)
                     image = [info, i, j]
                     displacement = centers(a)%position - centers(b)%position - MATMUL(cell, REAL(image, dp))
                     IF (SUM(displacement**2) >= radius**2) CYCLE
                     IF (derive) THEN
                        ddisplacement = center_response(a)%position - center_response(b)%position - &
                                        MATMUL(cell_response, REAL(image, dp))
                        IF (PRESENT(map)) THEN
                           IF (ALL(ddisplacement == 0.0_dp) .AND. ALL(center_response(a)%radii == 0.0_dp) .AND. &
                               ALL(center_response(b)%radii == 0.0_dp)) CYCLE
                        END IF
                        CALL upf_inverse_coulomb_pair_kernel(centers(a)%radii, centers(a)%l, centers(a)%indices, &
                                                             centers(b)%radii, centers(b)%l, centers(b)%indices, &
                                                             displacement, kernel, dkernel, center_response(a)%radii, &
                                                             center_response(b)%radii, ddisplacement)
                     ELSE
                        CALL upf_inverse_coulomb_pair_kernel(centers(a)%radii, centers(a)%l, centers(a)%indices, &
                                                             centers(b)%radii, centers(b)%l, centers(b)%indices, &
                                                             displacement, kernel)
                     END IF
                     DO channel = 1, SIZE(centers(a)%l)
                        DO nr = 1, na
                           IF (PRESENT(map)) THEN
                              dpair = dpair + MERGE(0.5_dp, 1.0_dp, a == b)*work(a)%q(nr, channel)* &
                                      SUM(dkernel(nr, channel, :, :)*work(b)%q)
                              CYCLE
                           END IF
                           norm = SUM(kernel(nr, channel, :, :)*work(b)%q)
                           work(a)%potential(nr, channel) = work(a)%potential(nr, channel) + norm
                           epair = epair + 0.5_dp*work(a)%q(nr, channel)*norm
                           IF (.NOT. derive) CYCLE
                           dpair = dpair + 0.5_dp*(work(a)%dq(nr, channel)*norm + work(a)%q(nr, channel)* &
                                                   SUM(kernel(nr, channel, :, :)*work(b)%dq + &
                                                       dkernel(nr, channel, :, :)*work(b)%q))
                        END DO
                     END DO
                  END DO
               END DO
            END DO
            DEALLOCATE (kernel)
            IF (derive) DEALLOCATE (dkernel)
         END DO
      END DO
      IF (reuse_pair_geometry) dpair = dpair + pair_geometry_response
      energy = energy + epair
      IF (derive) response = dpw + dpair
      IF (.NOT. PRESENT(potentials)) RETURN
      ALLOCATE (potentials(SIZE(centers)))
      DO a = 1, SIZE(centers)
         nr = SIZE(centers(a)%radii)
         nh = SIZE(centers(a)%l)
         ALLOCATE (potentials(a)%values(nr, nh), vh(nr, nh), ve(nr, nh), adjoint(nr, nh))
         IF (ALLOCATED(fourier(a)%values)) THEN
            kernel_values => fourier(a)%values
         ELSE
            kernel_values => map%fourier(a)%values
         END IF
         CALL fourier_density_adjoint(kernel_values, poisson_weights*h, vh)
         CALL fourier_density_adjoint(kernel_values, background_adjoint, ve)
         IF (PRESENT(comm)) THEN
            CALL comm%sum(vh)
            CALL comm%sum(ve)
         END IF
         CALL upf_inverse_coulomb_compensation_adjoint(centers(a)%radii, centers(a)%weights, centers(a)%l, &
                                                       centers(a)%shapes, vh + work(a)%potential, adjoint)
         potentials(a)%values(:, :) = adjoint + ve
         CALL upf_inverse_coulomb_compensation_adjoint(centers(a)%radii, centers(a)%weights, centers(a)%l, &
                                                       centers(a)%shapes, ve, adjoint)
         potentials(a)%values(:, :) = potentials(a)%values - adjoint
         DEALLOCATE (vh, ve, adjoint)
      END DO
   END SUBROUTINE upf_inverse_coulomb_periodic

! **************************************************************************************************
!> \brief Reject stale geometric maps rather than silently changing the Coulomb functional.
!> \param map ...
!> \param centers ...
!> \param cell ...
!> \param gvecs ...
!> \param volume ...
!> \param periodicity ...
! **************************************************************************************************
   SUBROUTINE check_map(map, centers, cell, gvecs, volume, periodicity)
      TYPE(upf_inverse_coulomb_map_type), INTENT(IN)     :: map
      TYPE(upf_inverse_coulomb_center_type), INTENT(IN)  :: centers(:)
      REAL(dp), INTENT(IN)                               :: cell(3, 3), gvecs(:, :), volume
      INTEGER, INTENT(IN)                                :: periodicity(3)

      INTEGER                                            :: a

      CPASSERT(ALLOCATED(map%centers) .AND. ALLOCATED(map%pairs))
      CPASSERT(SIZE(map%centers) == SIZE(centers))
      CPASSERT(ALL(map%cell == cell) .AND. map%volume == volume)
      CPASSERT(ALL(map%periodicity == periodicity))
      CPASSERT(ALL(SHAPE(map%gvecs) == SHAPE(gvecs)))
      CPASSERT(ALL(map%gvecs == gvecs))
      DO a = 1, SIZE(centers)
         CPASSERT(ALL(map%centers(a)%position == centers(a)%position))
         CPASSERT(SIZE(map%centers(a)%radii) == SIZE(centers(a)%radii))
         CPASSERT(SIZE(map%centers(a)%l) == SIZE(centers(a)%l))
         CPASSERT(ALL(map%centers(a)%radii == centers(a)%radii))
         CPASSERT(ALL(map%centers(a)%l == centers(a)%l))
         CPASSERT(ALL(map%centers(a)%indices == centers(a)%indices))
      END DO
   END SUBROUTINE check_map

! **************************************************************************************************
!> \brief Bound overlapping residual images, without repeating nonperiodic lattice directions.
!> \param inverse_cell ...
!> \param displacement ...
!> \param radius ...
!> \param periodicity ...
!> \param limits ...
! **************************************************************************************************
   SUBROUTINE residual_image_limits(inverse_cell, displacement, radius, periodicity, limits)
      REAL(dp), INTENT(IN)                               :: inverse_cell(3, 3), displacement(3), &
                                                            radius
      INTEGER, INTENT(IN)                                :: periodicity(3)
      INTEGER, INTENT(OUT)                               :: limits(2, 3)

      INTEGER                                            :: direction
      REAL(dp)                                           :: bound(3), fractional(3)

      fractional = MATMUL(inverse_cell, displacement)
      bound = radius*SQRT(SUM(inverse_cell**2, DIM=2))
      DO direction = 1, 3
         IF (periodicity(direction) == 1) THEN
            limits(1, direction) = CEILING(fractional(direction) - bound(direction))
            limits(2, direction) = FLOOR(fractional(direction) + bound(direction))
         ELSE
            limits(:, direction) = 0
         END IF
      END DO
   END SUBROUTINE residual_image_limits

! **************************************************************************************************
!> \brief Select one half of the opposite lattice images of a center.
!> \param i First lattice index.
!> \param j Second lattice index.
!> \param k Third lattice index.
!> \return True for the negative member of a nonzero image pair.
! **************************************************************************************************
   PURE LOGICAL FUNCTION negative_image(i, j, k) RESULT(negative)
      INTEGER, INTENT(IN)                                :: i, j, k

      negative = k < 0 .OR. (k == 0 .AND. (j < 0 .OR. (j == 0 .AND. i < 0)))
   END FUNCTION negative_image

! **************************************************************************************************
!> \brief Check whether independent Cartesian reflections preserve the lattice.
!> \param cell Lattice vectors in Cartesian columns.
!> \return True for a diagonal cell matrix.
! **************************************************************************************************
   PURE LOGICAL FUNCTION cartesian_reflection_cell(cell) RESULT(cartesian)
      REAL(dp), INTENT(IN)                               :: cell(3, 3)

      cartesian = ALL(cell([2, 3], 1) == 0.0_dp) .AND. ALL(cell([1, 3], 2) == 0.0_dp) .AND. &
                  ALL(cell([1, 2], 3) == 0.0_dp)
   END FUNCTION cartesian_reflection_cell

! **************************************************************************************************
!> \brief Find Cartesian axes where reflecting the displacement maps integer images to integer images.
!> \param cell Lattice vectors in Cartesian columns.
!> \param separation Difference of the two center positions, before subtracting an image.
!> \param periodicity Periodic lattice directions.
!> \param cartesian Whether the lattice is diagonal in the harmonic frame.
!> \param axes Axes with an exact integer or half-integer reflection center.
!> \param center Fractional position of the reflection center on each selected axis.
! **************************************************************************************************
   PURE SUBROUTINE pair_reflection_center(cell, separation, periodicity, cartesian, axes, center)
      REAL(dp), INTENT(IN)                               :: cell(3, 3), separation(3)
      INTEGER, INTENT(IN)                                :: periodicity(3)
      LOGICAL, INTENT(IN)                                :: cartesian
      LOGICAL, INTENT(OUT)                               :: axes(3)
      REAL(dp), INTENT(OUT)                              :: center(3)

      INTEGER                                            :: axis
      REAL(dp)                                           :: twice_center

      axes = .FALSE.
      center = 0.0_dp
      IF (.NOT. cartesian) RETURN
      DO axis = 1, 3
         IF (periodicity(axis) == 0) CYCLE
         twice_center = ANINT(2.0_dp*separation(axis)/cell(axis, axis))
         IF (2.0_dp*separation(axis) /= twice_center*cell(axis, axis)) CYCLE
         axes(axis) = .TRUE.
         center(axis) = 0.5_dp*twice_center
      END DO
   END SUBROUTINE pair_reflection_center

! **************************************************************************************************
!> \brief Encode the odd Cartesian reflection axes of each product of real solid harmonics.
!> \param first First native harmonic layout.
!> \param second Second native harmonic layout.
!> \return Bit masks for the two shell channels, with one bit per Cartesian axis.
! **************************************************************************************************
   FUNCTION image_pair_parity(first, second) RESULT(parity)
      TYPE(upf_inverse_coulomb_center_type), INTENT(IN)  :: first, second
      INTEGER :: parity(SIZE(first%l), SIZE(second%l))
      INTEGER :: first_bits(SIZE(first%l)), second_bits(SIZE(second%l))

      first_bits = harmonic_parity(first)
      second_bits = harmonic_parity(second)
      parity = IEOR(SPREAD(first_bits, 2, SIZE(second_bits)), SPREAD(second_bits, 1, SIZE(first_bits)))
   END FUNCTION image_pair_parity

! **************************************************************************************************
!> \brief Return the odd Cartesian reflection axes of the individual real solid harmonics.
!> \param center Native harmonic layout.
!> \return One bit mask per channel.
! **************************************************************************************************
   FUNCTION harmonic_parity(center) RESULT(bits)
      TYPE(upf_inverse_coulomb_center_type), INTENT(IN)  :: center
      INTEGER :: bits(SIZE(center%l))

      INTEGER                                            :: channel, l, monomial, row

      CALL init_orbital_pointers(MAXVAL(center%l))
      CALL init_spherical_harmonics(MAXVAL(center%l), -1)
      DO channel = 1, SIZE(center%l)
         l = center%l(channel)
         row = center%indices(channel) - l**2
         ! The nonzero monomials of each real solid harmonic have the same parity.
         monomial = MAXLOC(ABS(orbtramat(l)%slm(row, :)), DIM=1)
         bits(channel) = SUM(reflection_bits*MOD(indco(:, ncoset(l - 1) + monomial), 2))
      END DO
   END FUNCTION harmonic_parity

! **************************************************************************************************
!> \brief Reciprocal Coulomb map between two translated harmonic radial-shell distributions.
!>        The kernel excludes both radial measures. Splitting the angular integral at the source
!>        shell avoids sampling a discontinuous radial derivative across an unsplit angular rule.
!>        The joint response differentiates the same finite rule, including its moving boundaries.
!> \param radii_a Positive target radii.
!> \param l_a Angular momenta of the target real harmonics.
!> \param indices_a Native harmonic indices of the target channels.
!> \param radii_b Positive source radii.
!> \param l_b Angular momenta of the source real harmonics.
!> \param indices_b Native harmonic indices of the source channels.
!> \param displacement Target center minus source center, including any lattice image.
!> \param kernel Coulomb map in (radius_a,harmonic_a,radius_b,harmonic_b) order.
!> \param response Optional joint geometry response.
!> \param dr_a Required with response.
!> \param dr_b Required with response.
!> \param ddisplacement Required with response.
! **************************************************************************************************
   SUBROUTINE upf_inverse_coulomb_pair_kernel(radii_a, l_a, indices_a, radii_b, l_b, indices_b, &
                                              displacement, kernel, response, dr_a, dr_b, ddisplacement)
      REAL(dp), INTENT(IN)                               :: radii_a(:)
      INTEGER, INTENT(IN)                                :: l_a(:), indices_a(:)
      REAL(dp), INTENT(IN)                               :: radii_b(:)
      INTEGER, INTENT(IN)                                :: l_b(:), indices_b(:)
      REAL(dp), INTENT(IN)                               :: displacement(3)
      REAL(dp), INTENT(OUT)                              :: kernel(:, :, :, :)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: response(:, :, :, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: dr_a(:), dr_b(:), ddisplacement(3)

      INTEGER                                            :: a, b, i, j, nq, sample_count
      LOGICAL                                            :: separated
      REAL(dp)                                           :: da, db, delta(3), distance, factor, value
      REAL(dp), ALLOCATABLE                              :: ab(:, :), ba(:, :), dab(:, :), &
                                                            dba(:, :), dpoint(:, :), nodes(:), &
                                                            point(:, :), weights(:)
      TYPE(contained_shell_type), ALLOCATABLE            :: contained_a(:), contained_b(:)
      TYPE(harmonic_work_type)                           :: work_a, work_b

      CPASSERT(ALL(radii_a > 0.0_dp) .AND. ALL(ieee_is_finite(radii_a)))
      CPASSERT(ALL(radii_b > 0.0_dp) .AND. ALL(ieee_is_finite(radii_b)))
      CPASSERT(ALL(ieee_is_finite(displacement)))
      CPASSERT(SIZE(l_a) > 0 .AND. SIZE(l_b) > 0)
      CPASSERT(SIZE(indices_a) == SIZE(l_a))
      CPASSERT(SIZE(indices_b) == SIZE(l_b))
      CPASSERT(ALL(l_a >= 0 .AND. l_a <= 20))
      CPASSERT(ALL(l_b >= 0 .AND. l_b <= 20))
      CPASSERT(ALL(indices_a > l_a**2 .AND. indices_a <= (l_a + 1)**2))
      CPASSERT(ALL(indices_b > l_b**2 .AND. indices_b <= (l_b + 1)**2))
      IF (.NOT. ALL(SHAPE(kernel) == [SIZE(radii_a), SIZE(l_a), SIZE(radii_b), SIZE(l_b)])) THEN
         CPABORT("Invalid translated shell kernel dimensions")
      END IF
      IF (PRESENT(response)) THEN
         CPASSERT(PRESENT(dr_a) .AND. PRESENT(dr_b))
         CPASSERT(PRESENT(ddisplacement))
         CPASSERT(SIZE(dr_a) == SIZE(radii_a))
         CPASSERT(SIZE(dr_b) == SIZE(radii_b))
         CPASSERT(ALL(SHAPE(response) == SHAPE(kernel)))
         response = 0.0_dp
      END IF
      kernel = 0.0_dp
      distance = SQRT(SUM(displacement**2))
      delta = 0.0_dp
      IF (PRESENT(response)) delta = ddisplacement
      IF (distance == 0.0_dp) THEN
         ! Same-center shells share their center motion. Distinct coincident centers need a limit.
         CPASSERT(ALL(delta == 0.0_dp))
         DO a = 1, SIZE(l_a)
            DO b = 1, SIZE(l_b)
               IF (indices_a(a) /= indices_b(b)) CYCLE
               DO i = 1, SIZE(radii_a)
                  DO j = 1, SIZE(radii_b)
                     value = fourpi/REAL(2*l_a(a) + 1, dp)/MAX(radii_a(i), radii_b(j))* &
                             (MIN(radii_a(i), radii_b(j))/MAX(radii_a(i), radii_b(j)))**l_a(a)
                     kernel(i, a, j, b) = value
                     IF (.NOT. PRESENT(response)) CYCLE
                     IF (radii_a(i) == radii_b(j)) THEN
                        CPASSERT(dr_a(i) == dr_b(j))
                     END IF
                     IF (radii_a(i) >= radii_b(j)) THEN
                        response(i, a, j, b) = value*(l_a(a)*dr_b(j)/radii_b(j) - &
                                                      (l_a(a) + 1)*dr_a(i)/radii_a(i))
                     ELSE
                        response(i, a, j, b) = value*(l_a(a)*dr_a(i)/radii_a(i) - &
                                                      (l_a(a) + 1)*dr_b(j)/radii_b(j))
                     END IF
                  END DO
               END DO
            END DO
         END DO
         RETURN
      END IF
      CALL init_orbital_pointers(MAX(MAXVAL(l_a), MAXVAL(l_b)))
      CALL init_spherical_harmonics(MAX(MAXVAL(l_a), MAXVAL(l_b)), -1)
      ALLOCATE (ab(SIZE(l_a), SIZE(l_b)), ba(SIZE(l_b), SIZE(l_a)))
      IF (PRESENT(response)) ALLOCATE (dab(SIZE(l_a), SIZE(l_b)), dba(SIZE(l_b), SIZE(l_a)))
      nq = MAX(32, MAXVAL(l_a) + MAXVAL(l_b) + 2)
      sample_count = nq*(MAXVAL(l_a) + MAXVAL(l_b) + 1)
      CALL prepare_harmonic_work(l_a, indices_a, sample_count, PRESENT(response), work_a)
      CALL prepare_harmonic_work(l_b, indices_b, sample_count, PRESENT(response), work_b)
      ALLOCATE (nodes(nq), weights(nq))
      CALL polar_rule(nodes, weights)
      ALLOCATE (contained_a(SIZE(radii_a)), contained_b(SIZE(radii_b)))
      ! Separated harmonic shells have an exact point-multipole interaction.
      ! Higher orders retain the general angular rule, without discarding any channel.
      separated = MAXVAL(l_a) + MAXVAL(l_b) <= 8 .AND. distance > MINVAL(radii_a) + MINVAL(radii_b)
      IF (separated) THEN
         ALLOCATE (point(SIZE(l_a), SIZE(l_b)))
         IF (PRESENT(response)) THEN
            ALLOCATE (dpoint, MOLD=point)
            CALL separated_shell_map(l_a, indices_a, l_b, indices_b, displacement, point, dpoint, delta)
         ELSE
            CALL separated_shell_map(l_a, indices_a, l_b, indices_b, displacement, point)
         END IF
      END IF
      DO i = 1, SIZE(radii_a)
         DO j = 1, SIZE(radii_b)
            da = 0.0_dp
            db = 0.0_dp
            IF (PRESENT(response)) THEN
               da = dr_a(i)
               db = dr_b(j)
            END IF
            IF (separated .AND. distance > radii_a(i) + radii_b(j)) THEN
               DO b = 1, SIZE(l_b)
                  DO a = 1, SIZE(l_a)
                     factor = radii_a(i)**l_a(a)*radii_b(j)**l_b(b)
                     kernel(i, a, j, b) = factor*point(a, b)
                     IF (PRESENT(response)) response(i, a, j, b) = factor*(dpoint(a, b) + point(a, b)* &
                                                                           (l_a(a)*da/radii_a(i) + &
                                                                            l_b(b)*db/radii_b(j)))
                  END DO
               END DO
               CYCLE
            END IF
            IF (ABS(radii_a(i) - radii_b(j)) <= distance) CYCLE
            IF (PRESENT(response)) THEN
               CALL contained_shell(radii_a(i), radii_b(j), l_a, l_b, &
                                    displacement, da, db, delta, ab, nodes, weights, work_a, work_b, contained_a(i), dab)
               CALL contained_shell(radii_b(j), radii_a(i), l_b, l_a, &
                                    -displacement, db, da, -delta, ba, nodes, weights, work_b, work_a, contained_b(j), dba)
            ELSE
               CALL contained_shell(radii_a(i), radii_b(j), l_a, l_b, &
                                    displacement, da, db, delta, ab, nodes, weights, work_a, work_b, contained_a(i))
               CALL contained_shell(radii_b(j), radii_a(i), l_b, l_a, &
                                    -displacement, db, da, -delta, ba, nodes, weights, work_b, work_a, contained_b(j))
            END IF
            kernel(i, :, j, :) = 0.5_dp*(ab + TRANSPOSE(ba))
            IF (PRESENT(response)) response(i, :, j, :) = 0.5_dp*(dab + TRANSPOSE(dba))
         END DO
      END DO

      ! The remaining pairs have independent angular quadratures. Keep the shared
      ! contained-shell caches in the serial loop so their reference radii are unchanged.
!$OMP PARALLEL DEFAULT(NONE) IF(SIZE(radii_a) >= 32 .AND. SIZE(radii_b) >= 32) &
!$OMP SHARED(radii_a, radii_b, l_a, l_b, indices_a, indices_b, distance, displacement, delta, &
!$OMP        separated, kernel, response, dr_a, dr_b, sample_count, nodes, weights) &
!$OMP PRIVATE(i, j, da, db, ab, ba, dab, dba, work_a, work_b)
      CALL prepare_harmonic_work(l_a, indices_a, sample_count, PRESENT(response), work_a)
      CALL prepare_harmonic_work(l_b, indices_b, sample_count, PRESENT(response), work_b)
!$OMP DO SCHEDULE(DYNAMIC)
      DO i = 1, SIZE(radii_a)
         DO j = 1, SIZE(radii_b)
            IF (separated .AND. distance > radii_a(i) + radii_b(j)) CYCLE
            IF (ABS(radii_a(i) - radii_b(j)) > distance) CYCLE
            IF (PRESENT(response)) THEN
               CALL translated_shell(radii_a(i), radii_b(j), l_a, l_b, &
                                     displacement, dr_a(i), dr_b(j), delta, ab, nodes, weights, work_a, work_b, dab)
               CALL translated_shell(radii_b(j), radii_a(i), l_b, l_a, &
                                     -displacement, dr_b(j), dr_a(i), -delta, ba, nodes, weights, work_b, work_a, dba)
            ELSE
               CALL translated_shell(radii_a(i), radii_b(j), l_a, l_b, &
                                     displacement, 0.0_dp, 0.0_dp, delta, ab, nodes, weights, work_a, work_b)
               CALL translated_shell(radii_b(j), radii_a(i), l_b, l_a, &
                                     -displacement, 0.0_dp, 0.0_dp, -delta, ba, nodes, weights, work_b, work_a)
            END IF
            kernel(i, :, j, :) = 0.5_dp*(ab + TRANSPOSE(ba))
            IF (PRESENT(response)) response(i, :, j, :) = 0.5_dp*(dab + TRANSPOSE(dba))
         END DO
      END DO
!$OMP END DO
!$OMP END PARALLEL
   END SUBROUTINE upf_inverse_coulomb_pair_kernel

! **************************************************************************************************
!> \brief Reuse the same angular rule when one radial shell fully contains the other.
!>        At fixed target radius its nodes and boundaries do not depend on the source radius.
!>        The source potential scales as b**(-l-1) inside and b**l outside its shell.
!>        Cache the response at fixed source radius and add its explicit radial derivative.
!> \param a Target radius, fixed across calls using the same cache.
!> \param b Source radius.
!> \param la Target harmonic orders.
!> \param lb Source harmonic orders.
!> \param d Target center minus source center, fixed for this cache.
!> \param da Target radial motion, fixed for this cache.
!> \param db Source radial motion.
!> \param dd Center displacement response, fixed for this cache.
!> \param value Angular Coulomb map.
!> \param nodes Unchanged polar quadrature nodes.
!> \param weights Unchanged polar quadrature weights.
!> \param work_a Target harmonic scratch fields.
!> \param work_b Source harmonic scratch fields.
!> \param cache Owned maps for source shells enclosing or enclosed by the target.
!> \param response Optional joint response, requested consistently for this cache.
! **************************************************************************************************
   SUBROUTINE contained_shell(a, b, la, lb, d, da, db, dd, value, nodes, weights, work_a, work_b, cache, response)
      REAL(dp), INTENT(IN)                               :: a, b
      INTEGER, INTENT(IN)                                :: la(:), lb(:)
      REAL(dp), INTENT(IN)                               :: d(3), da, db, dd(3)
      REAL(dp), INTENT(OUT)                              :: value(:, :)
      REAL(dp), INTENT(IN)                               :: nodes(:), weights(:)
      TYPE(harmonic_work_type), INTENT(INOUT)            :: work_a, work_b
      TYPE(contained_shell_type), INTENT(INOUT)          :: cache
      REAL(dp), INTENT(OUT), OPTIONAL                    :: response(:, :)

      INTEGER                                            :: channel, power, side
      REAL(dp)                                           :: factor

      side = MERGE(1, 2, b > a)
      IF (.NOT. ALLOCATED(cache%values)) THEN
         ALLOCATE (cache%values(SIZE(la), SIZE(lb), 2))
         IF (PRESENT(response)) ALLOCATE (cache%response, MOLD=cache%values)
      END IF
      IF (.NOT. cache%ready(side)) THEN
         IF (PRESENT(response)) THEN
            CALL translated_shell(a, b, la, lb, d, da, 0.0_dp, dd, cache%values(:, :, side), &
                                  nodes, weights, work_a, work_b, cache%response(:, :, side))
         ELSE
            CALL translated_shell(a, b, la, lb, d, da, 0.0_dp, dd, cache%values(:, :, side), &
                                  nodes, weights, work_a, work_b)
         END IF
         cache%source_radius(side) = b
         cache%ready(side) = .TRUE.
      END IF
      DO channel = 1, SIZE(lb)
         IF (side == 1) THEN
            power = -lb(channel) - 1
            factor = (cache%source_radius(side)/b)**(-power)
         ELSE
            power = lb(channel)
            factor = (b/cache%source_radius(side))**power
         END IF
         value(:, channel) = factor*cache%values(:, channel, side)
         IF (PRESENT(response)) THEN
            response(:, channel) = factor*(cache%response(:, channel, side) + &
                                           power*db/b*cache%values(:, channel, side))
         END IF
      END DO
   END SUBROUTINE contained_shell

! **************************************************************************************************
!> \brief Point-multipole map for disjoint shells in the native real solid-harmonic convention.
!>        K_ab = (4*pi)^2*(-1)^lb*Ra(partial)*Rb(partial)(1/d)/[(2*la+1)!!*(2*lb+1)!!].
!>        Shell radii enter separately as a**la*b**lb. Cartesian Coulomb derivatives are exact.
!> \param la ...
!> \param ia ...
!> \param lb ...
!> \param ib ...
!> \param d Nonzero target-minus-source displacement.
!> \param values Kernel without the two shell-radius powers.
!> \param response Optional displacement response of the same kernel.
!> \param dd Required with response.
! **************************************************************************************************
   SUBROUTINE separated_shell_map(la, ia, lb, ib, d, values, response, dd)
      INTEGER, INTENT(IN)                                :: la(:), ia(:), lb(:), ib(:)
      REAL(dp), INTENT(IN)                               :: d(3)
      REAL(dp), INTENT(OUT)                              :: values(:, :)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: response(:, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: dd(3)

      INTEGER                                            :: a, axis, b, bx, by, bz, contractions, i, &
                                                            ix, iy, iz, j, n, order, powers(3), &
                                                            reduced(3), shifted(3)
      REAL(dp)                                           :: coefficient, direction(3), distance, &
                                                            factor, term
      REAL(dp), ALLOCATABLE                              :: factorial(:), odd_factorial(:), &
                                                            tensor(:, :, :)

      CPASSERT(PRESENT(response) .EQV. PRESENT(dd))
      n = MAXVAL(la) + MAXVAL(lb)
      IF (PRESENT(response)) n = n + 1
      distance = SQRT(SUM(d**2))
      CPASSERT(distance > 0.0_dp)
      direction = d/distance
      ALLOCATE (factorial(0:n + 1), odd_factorial(0:n + 1), tensor(0:n, 0:n, 0:n))
      factorial(0) = 1.0_dp
      odd_factorial(0) = 1.0_dp
      DO i = 1, n + 1
         factorial(i) = i*factorial(i - 1)
         odd_factorial(i) = (2*i - 1)*odd_factorial(i - 1)
      END DO
      tensor = 0.0_dp
      DO iz = 0, n
         DO iy = 0, n - iz
            DO ix = 0, n - iz - iy
               powers = [ix, iy, iz]
               order = SUM(powers)
               factor = PRODUCT(factorial(powers))/distance**(order + 1)
               DO bz = 0, iz/2
                  DO by = 0, iy/2
                     DO bx = 0, ix/2
                        contractions = bx + by + bz
                        reduced = powers - 2*[bx, by, bz]
                        term = REAL((-1)**(order - contractions), dp)*odd_factorial(order - contractions)* &
                               PRODUCT(direction**reduced)/(PRODUCT(factorial([bx, by, bz]))* &
                                                            PRODUCT(factorial(reduced))*2.0_dp**contractions)
                        tensor(ix, iy, iz) = tensor(ix, iy, iz) + factor*term
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
      values = 0.0_dp
      IF (PRESENT(response)) response = 0.0_dp
      DO b = 1, SIZE(lb)
         DO a = 1, SIZE(la)
            factor = fourpi**2*REAL((-1)**lb(b), dp)/(odd_factorial(la(a) + 1)*odd_factorial(lb(b) + 1))
            DO j = 1, nco(lb(b))
               DO i = 1, nco(la(a))
                  coefficient = factor*orbtramat(la(a))%slm(ia(a) - la(a)**2, i)* &
                                orbtramat(lb(b))%slm(ib(b) - lb(b)**2, j)
                  IF (coefficient == 0.0_dp) CYCLE
                  powers = indco(:, ncoset(la(a) - 1) + i) + indco(:, ncoset(lb(b) - 1) + j)
                  values(a, b) = values(a, b) + coefficient*tensor(powers(1), powers(2), powers(3))
                  IF (.NOT. PRESENT(response)) CYCLE
                  DO axis = 1, 3
                     shifted = powers
                     shifted(axis) = shifted(axis) + 1
                     response(a, b) = response(a, b) + &
                                      coefficient*dd(axis)*tensor(shifted(1), shifted(2), shifted(3))
                  END DO
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE separated_shell_map

! **************************************************************************************************
!> \brief Integrate a source shell's regular/irregular solid-harmonic potential on a target shell.
!>        Azimuthal integration is exact for the represented polynomial degree. Polar subintervals
!>        also limit the source-distance ratio to two in the irregular branch.
!> \param a ...
!> \param b ...
!> \param la ...
!> \param lb ...
!> \param d ...
!> \param da ...
!> \param db ...
!> \param dd ...
!> \param value ...
!> \param nodes Shared polar quadrature nodes, independent of the two shell radii.
!> \param weights Shared polar quadrature weights.
!> \param work_a Harmonic fields reused across target shells.
!> \param work_b Harmonic fields reused across source shells.
!> \param response ...
! **************************************************************************************************
   SUBROUTINE translated_shell(a, b, la, lb, d, da, db, dd, value, nodes, weights, work_a, work_b, response)
      REAL(dp), INTENT(IN)                               :: a, b
      INTEGER, INTENT(IN)                                :: la(:), lb(:)
      REAL(dp), INTENT(IN)                               :: d(3), da, db, dd(3)
      REAL(dp), CONTIGUOUS, INTENT(OUT)                  :: value(:, :)
      REAL(dp), INTENT(IN)                               :: nodes(:), weights(:)
      TYPE(harmonic_work_type), INTENT(INOUT)            :: work_a, work_b
      REAL(dp), CONTIGUOUS, INTENT(OUT), OPTIONAL        :: response(:, :)

      INTEGER                                            :: axis, channel, iphi, iq, max_order, &
                                                            min_order, nphi, nq, order, sample
      LOGICAL                                            :: active_order(0:20), derive, regular
      REAL(dp) :: azimuth, azimuth_directions(3, MAXVAL(la) + MAXVAL(lb) + 1), &
         azimuth_response(3, MAXVAL(la) + MAXVAL(lb) + 1), bound, dbound, dex(3), dey(3), dez(3), &
         dfactor, distance, dleft, dlength, dmu, dright, dt, du(3), dweight, dx(3), ex(3), ey(3), &
         ez(3), factor, left, length, mu, norm_x, radius_dot, reference(3), right, scale, sine, &
         source_factor(0:20), source_response(0:20), t, tmax, u(3), weight, x(3)

      EXTERNAL dgemm

      derive = PRESENT(response)
      distance = SQRT(SUM(d**2))
      ez = d/distance
      reference = 0.0_dp
      axis = MINLOC(ABS(ez), DIM=1)
      reference(axis) = 1.0_dp
      ex = reference - ez(axis)*ez
      length = SQRT(SUM(ex**2))
      IF (derive) THEN
         dt = DOT_PRODUCT(d, dd)/distance
         dez = (dd - ez*dt)/distance
         dex = -dez(axis)*ez - ez(axis)*dez
         dlength = DOT_PRODUCT(ex, dex)/length
         dex = (dex - ex*dlength/length)/length
      END IF
      ex = ex/length
      ey = cross(ez, ex)
      IF (derive) dey = cross(dez, ex) + cross(ez, dex)
      nq = SIZE(nodes)
      nphi = MAXVAL(la) + MAXVAL(lb) + 1
      DO iphi = 1, nphi
         azimuth = twopi*REAL(iphi - 1, dp)/REAL(nphi, dp)
         azimuth_directions(:, iphi) = COS(azimuth)*ex + SIN(azimuth)*ey
         IF (derive) azimuth_response(:, iphi) = COS(azimuth)*dex + SIN(azimuth)*dey
      END DO
      min_order = MINVAL(lb)
      max_order = MAXVAL(lb)
      active_order = .FALSE.
      DO channel = 1, SIZE(lb)
         active_order(lb(channel)) = .TRUE.
      END DO
      value = 0.0_dp
      IF (derive) response = 0.0_dp
      tmax = distance + a
      left = -1.0_dp
      dleft = 0.0_dp
      regular = ABS(distance - a) < b
      DO WHILE (left < 1.0_dp)
         right = 1.0_dp
         dright = 0.0_dp
         t = SQRT(MAX(0.0_dp, (distance - a)**2 + 2.0_dp*distance*a*(left + 1.0_dp)))
         IF (regular .AND. b < tmax) THEN
            bound = b
            IF (derive) dbound = db
         ELSE IF (.NOT. regular .AND. 2.0_dp*t < tmax) THEN
            bound = 2.0_dp*t
            IF (derive) dbound = 2.0_dp*((distance - a)*(dt - da) + &
                                         (dt*a + distance*da)*(left + 1.0_dp) + distance*a*dleft)/t
         ELSE
            bound = tmax
            IF (derive) dbound = dt + da
         END IF
         IF (bound < tmax) THEN
            right = ((bound - a)*(bound + a) - distance**2)/(2.0_dp*distance*a)
            IF (derive) dright = (bound*dbound - a*da - distance*dt)/(distance*a) - &
                                 right*(dt/distance + da/a)
         END IF
         ! At shell tangency, roundoff can leave a cap too narrow to place interior nodes.
         ! The two potentials agree on the shell, so its moving boundary terms cancel.
         IF (right - left <= 64.0_dp*EPSILON(1.0_dp)) THEN
            CPASSERT(regular)
            regular = .FALSE.
            CYCLE
         END IF
         IF (1.0_dp - right <= 64.0_dp*EPSILON(1.0_dp)) THEN
            right = 1.0_dp
            dright = 0.0_dp
         END IF
         CPASSERT(right > left .AND. right <= 1.0_dp)
         DO iq = 1, nq
            mu = 0.5_dp*((1.0_dp - nodes(iq))*left + (1.0_dp + nodes(iq))*right)
            IF (derive) dmu = 0.5_dp*((1.0_dp - nodes(iq))*dleft + (1.0_dp + nodes(iq))*dright)
            weight = 0.5_dp*(right - left)*weights(iq)*twopi/REAL(nphi, dp)
            IF (derive) dweight = 0.5_dp*(dright - dleft)*weights(iq)*twopi/REAL(nphi, dp)
            sine = SQRT(MAX(0.0_dp, 1.0_dp - mu**2))
            CPASSERT(sine > 0.0_dp)
            ! The source distance and its derivative are independent of the azimuth.
            norm_x = SQRT((distance - a)**2 + 2.0_dp*distance*a*(mu + 1.0_dp))
            IF (derive) radius_dot = (distance - a)*(dt - da) + &
                                     (dt*a + distance*da)*(mu + 1.0_dp) + distance*a*dmu
            DO order = min_order, max_order
               IF (.NOT. active_order(order)) CYCLE
               IF (norm_x < b) THEN
                  factor = fourpi/REAL(2*order + 1, dp)/b**(order + 1)
                  IF (derive) dfactor = -factor*(order + 1)*db/b
               ELSE
                  factor = fourpi/REAL(2*order + 1, dp)*(b/norm_x)**order/norm_x**(order + 1)
                  IF (derive) dfactor = factor*(order*db/b - (2*order + 1)*radius_dot/norm_x**2)
               END IF
               source_factor(order) = factor
               IF (derive) source_response(order) = dfactor
            END DO
            DO iphi = 1, nphi
               sample = (iq - 1)*nphi + iphi
               u = mu*ez + sine*azimuth_directions(:, iphi)
               x = d + a*u
               CALL solid_fields(u, work_a, derive)
               CALL solid_fields(x, work_b, derive)
               IF (derive) THEN
                  du = dmu*ez + mu*dez - mu*dmu/sine*azimuth_directions(:, iphi) + &
                       sine*azimuth_response(:, iphi)
                  dx = dd + da*u + a*du
                  DO channel = 1, SIZE(la)
                     work_a%response(channel) = DOT_PRODUCT(du, work_a%gradient(:, channel))
                  END DO
                  DO channel = 1, SIZE(lb)
                     work_b%response(channel) = DOT_PRODUCT(dx, work_b%gradient(:, channel))
                  END DO
               END IF
               work_a%sampled_values(:, sample) = work_a%values
               IF (derive) work_a%sampled_response(:, sample) = work_a%response
               DO channel = 1, SIZE(lb)
                  factor = source_factor(lb(channel))
                  scale = factor*work_b%values(channel)
                  work_b%sampled_values(channel, sample) = weight*scale
                  IF (derive) THEN
                     dfactor = source_response(lb(channel))
                     work_b%sampled_response(channel, sample) = dweight*scale + &
                        weight*(dfactor*work_b%values(channel) + factor*work_b%response(channel))
                  END IF
               END DO
            END DO
         END DO
         ! Each polar interval uses the original quadrature. Include the moving
         ! quadrature weight in the source derivative before applying the product rule.
         CALL dgemm('N', 'T', SIZE(la), SIZE(lb), nq*nphi, 1.0_dp, work_a%sampled_values, SIZE(la), &
                    work_b%sampled_values, SIZE(lb), 1.0_dp, value, SIZE(la))
         IF (derive) THEN
            CALL dgemm('N', 'T', SIZE(la), SIZE(lb), nq*nphi, 1.0_dp, work_a%sampled_response, SIZE(la), &
                       work_b%sampled_values, SIZE(lb), 1.0_dp, response, SIZE(la))
            CALL dgemm('N', 'T', SIZE(la), SIZE(lb), nq*nphi, 1.0_dp, work_a%sampled_values, SIZE(la), &
                       work_b%sampled_response, SIZE(lb), 1.0_dp, response, SIZE(la))
         END IF
         left = right
         dleft = dright
         regular = .FALSE.
      END DO
   END SUBROUTINE translated_shell

! **************************************************************************************************
!> \brief Retain the nonzero Cartesian terms in the native solid harmonics.
!>        The translated shells reuse these coefficients at every angular node.
!> \param ls ...
!> \param indices ...
!> \param sample_count Angular nodes in one polar interval.
!> \param derive Allocate scratch space for geometry derivatives.
!> \param work Prepared harmonic terms and scratch fields.
! **************************************************************************************************
   SUBROUTINE prepare_harmonic_work(ls, indices, sample_count, derive, work)
      INTEGER, INTENT(IN)                                :: ls(:), indices(:)
      INTEGER, INTENT(IN)                                :: sample_count
      LOGICAL, INTENT(IN)                                :: derive
      TYPE(harmonic_work_type), INTENT(OUT)              :: work

      INTEGER                                            :: channel, i, l, m, term
      REAL(dp)                                           :: coefficient

      work%max_degree = MAXVAL(ls)
      ALLOCATE (work%first_term(SIZE(ls) + 1), work%values(SIZE(ls)), work%sampled_values(SIZE(ls), sample_count))
      IF (derive) ALLOCATE (work%response(SIZE(ls)), work%gradient(3, SIZE(ls)), &
                            work%sampled_response(SIZE(ls), sample_count))
      work%first_term(1) = 1
      DO channel = 1, SIZE(ls)
         l = ls(channel)
         m = indices(channel) - l*l
         work%first_term(channel + 1) = work%first_term(channel) + COUNT(orbtramat(l)%slm(m, :) /= 0.0_dp)
      END DO
      term = work%first_term(SIZE(ls) + 1) - 1
      ALLOCATE (work%powers(3, term), work%coefficients(term))
      term = 0
      DO channel = 1, SIZE(ls)
         l = ls(channel)
         m = indices(channel) - l*l
         DO i = 1, nco(l)
            coefficient = orbtramat(l)%slm(m, i)
            IF (coefficient == 0.0_dp) CYCLE
            term = term + 1
            work%powers(:, term) = indco(:, ncoset(l - 1) + i)
            work%coefficients(term) = coefficient
         END DO
      END DO
   END SUBROUTINE prepare_harmonic_work

! **************************************************************************************************
!> \brief Native regular solid harmonics and Cartesian gradients, including their origin limits.
!> \param x ...
!> \param work Prepared Cartesian terms and output fields.
!> \param derive Evaluate Cartesian gradients.
! **************************************************************************************************
   SUBROUTINE solid_fields(x, work, derive)
      REAL(dp), INTENT(IN)                               :: x(3)
      TYPE(harmonic_work_type), INTENT(INOUT)            :: work
      LOGICAL, INTENT(IN)                                :: derive

      INTEGER                                            :: axis, channel, i, power, px, py, pz
      REAL(dp)                                           :: dx, dy, dz, factor, value, xp, xpowers(0:20, 3), yp, zp
      xpowers(0, :) = 1.0_dp
      DO axis = 1, 3
         DO power = 1, work%max_degree
            xpowers(power, axis) = x(axis)*xpowers(power - 1, axis)
         END DO
      END DO
      DO channel = 1, SIZE(work%values)
         value = 0.0_dp
         dx = 0.0_dp
         dy = 0.0_dp
         dz = 0.0_dp
         DO i = work%first_term(channel), work%first_term(channel + 1) - 1
            px = work%powers(1, i)
            py = work%powers(2, i)
            pz = work%powers(3, i)
            xp = xpowers(px, 1)
            yp = xpowers(py, 2)
            zp = xpowers(pz, 3)
            factor = work%coefficients(i)
            value = value + factor*xp*yp*zp
            IF (.NOT. derive) CYCLE
            IF (px > 0) dx = dx + factor*px*xpowers(px - 1, 1)*yp*zp
            IF (py > 0) dy = dy + factor*py*xp*xpowers(py - 1, 2)*zp
            IF (pz > 0) dz = dz + factor*pz*xp*yp*xpowers(pz - 1, 3)
         END DO
         work%values(channel) = value
         IF (derive) work%gradient(:, channel) = [dx, dy, dz]
      END DO
   END SUBROUTINE solid_fields

! **************************************************************************************************
!> \brief Gauss-Legendre nodes and weights on (-1,1).
!> \param nodes ...
!> \param weights ...
! **************************************************************************************************
   SUBROUTINE polar_rule(nodes, weights)
      REAL(dp), INTENT(OUT)                              :: nodes(:), weights(:)

      INTEGER                                            :: i, iteration, j, n
      REAL(dp)                                           :: derivative, p, previous, saved, shift, x

      n = SIZE(nodes)
      CPASSERT(SIZE(weights) == n)
      DO i = 1, (n + 1)/2
         x = COS(pi*(REAL(i, dp) - 0.25_dp)/(REAL(n, dp) + 0.5_dp))
         DO iteration = 1, 20
            p = 1.0_dp
            previous = 0.0_dp
            DO j = 1, n
               saved = p
               p = ((2*j - 1)*x*p - (j - 1)*previous)/REAL(j, dp)
               previous = saved
            END DO
            derivative = n*(x*p - previous)/(x*x - 1.0_dp)
            shift = p/derivative
            x = x - shift
            IF (ABS(shift) < 2.0_dp*EPSILON(x)) EXIT
         END DO
         CPASSERT(iteration <= 20)
         nodes(i) = -x
         nodes(n + 1 - i) = x
         weights(i) = 2.0_dp/((1.0_dp - x*x)*derivative**2)
         weights(n + 1 - i) = weights(i)
      END DO
   END SUBROUTINE polar_rule

! **************************************************************************************************
!> \brief Cartesian vector product.
!> \param a ...
!> \param b ...
!> \return ...
! **************************************************************************************************
   FUNCTION cross(a, b) RESULT(c)
      REAL(dp), INTENT(IN)                               :: a(3), b(3)
      REAL(dp)                                           :: c(3)

      c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
   END FUNCTION cross

! **************************************************************************************************
!> \brief Fourier map of the same radial shells, for coupling to a periodic background.
!>        rho_G = sum_ir,lm weights(ir)*kernel(G,ir,lm)*density(ir,lm).
!>        The reciprocal covector must use the caller's actual full/half-grid multiplicities.
!>        No reciprocal Coulomb cutoff or G=0 electrostatic convention is selected here.
!> \param radii Positive radial shell nodes.
!> \param angular_momenta l of each orthonormal real harmonic channel.
!> \param harmonic_indices Native harmonic indices, between l*l+1 and (l+1)**2.
!> \param gvecs Cartesian reciprocal vectors, including zero if required by the caller.
!> \param position Cartesian center of these radial shells.
!> \param volume Physical cell volume.
!> \param kernel Complex Fourier coefficients per unit radial measure and density.
!> \param response Optional joint geometry response of the same Fourier map.
!> \param radii_response Required with response.
!> \param gvecs_response Required with response, including reciprocal-vector strain.
!> \param position_response Required with response.
!> \param volume_response Required with response.
!> \param cell Lattice vectors in columns, required for a nonzero grid_shape.
!> \param grid_shape Optional full-box real-grid dimensions. Zero or omission gives continuous coefficients.
!> \param cell_response Optional lattice derivative of alias jumps. Omission holds the jumps fixed.
! **************************************************************************************************
   RECURSIVE SUBROUTINE upf_inverse_coulomb_fourier_kernel(radii, angular_momenta, harmonic_indices, gvecs, position, volume, &
                                                           kernel, response, radii_response, gvecs_response, &
                                                           position_response, volume_response, cell, grid_shape, cell_response)
      REAL(dp), INTENT(IN)                               :: radii(:)
      INTEGER, INTENT(IN)                                :: angular_momenta(:), harmonic_indices(:)
      REAL(dp), INTENT(IN)                               :: gvecs(:, :), position(3), volume
      COMPLEX(dp), INTENT(OUT)                           :: kernel(:, :, :)
      COMPLEX(dp), INTENT(OUT), OPTIONAL                 :: response(:, :, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: radii_response(:), gvecs_response(:, :), &
                                                            position_response(3), volume_response, &
                                                            cell(3, 3)
      INTEGER, INTENT(IN), OPTIONAL                      :: grid_shape(3)
      REAL(dp), INTENT(IN), OPTIONAL                     :: cell_response(3, 3)

      COMPLEX(dp)                                        :: phase, phase_l, phase_response
      COMPLEX(dp), ALLOCATABLE                           :: alias_kernel(:, :, :), &
                                                            alias_response(:, :, :)
      INTEGER                                            :: a, alias, axis, channel, count, ig, &
                                                            info, ir, l, last_order, m, miller(3), &
                                                            order, powers(3), reduced(3)
      LOGICAL :: bessel_ready(0:MAXVAL(angular_momenta) + 1), valid_lattice
      REAL(dp) :: alias_steps(3, 3), angular, angular_gradient(3), deformation(3, 3), &
         direction(3), factor, g_response(3, 1), gradient(3), inverse_cell(3, 3), jl, jnext, &
         lattice_error, lattice_index(3), norm_g, radial, vectors(3, 8), &
         bessel_values(SIZE(radii), 0:MAXVAL(angular_momenta) + 1)

      CPASSERT(SIZE(harmonic_indices) == SIZE(angular_momenta))
      CPASSERT(SIZE(angular_momenta) > 0)
      CPASSERT(SIZE(gvecs, 1) == 3)
      IF (.NOT. ALL(SHAPE(kernel) == [SIZE(gvecs, 2), SIZE(radii), SIZE(angular_momenta)])) THEN
         CPABORT("Invalid radial Fourier kernel dimensions")
      END IF
      CPASSERT(ALL(ieee_is_finite(radii)) .AND. ALL(radii > 0.0_dp))
      CPASSERT(ALL(ieee_is_finite(gvecs)))
      CPASSERT(ALL(ieee_is_finite(position)))
      CPASSERT(ieee_is_finite(volume) .AND. volume > 0.0_dp)
      CPASSERT(ALL(angular_momenta >= 0 .AND. angular_momenta <= 20))
      CPASSERT(ALL(harmonic_indices > angular_momenta**2))
      CPASSERT(ALL(harmonic_indices <= (angular_momenta + 1)**2))
      IF (PRESENT(response)) THEN
         CPASSERT(ALL(SHAPE(response) == SHAPE(kernel)))
         CPASSERT(PRESENT(radii_response) .AND. PRESENT(gvecs_response))
         CPASSERT(PRESENT(position_response))
         CPASSERT(PRESENT(volume_response))
         CPASSERT(SIZE(radii_response) == SIZE(radii))
         CPASSERT(ALL(SHAPE(gvecs_response) == SHAPE(gvecs)))
      END IF
      CALL init_orbital_pointers(MAXVAL(angular_momenta))
      CALL init_spherical_harmonics(MAXVAL(angular_momenta), -1)
      ! Reciprocal vectors have independent shell coefficients and derivatives.
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) IF(SIZE(gvecs, 2) > 64) &
!$OMP SHARED(radii, angular_momenta, harmonic_indices, gvecs, position, volume, kernel, response, &
!$OMP        radii_response, gvecs_response, position_response, volume_response, nco, ncoset, indco, orbtramat) &
!$OMP PRIVATE(ig, norm_g, bessel_ready, direction, phase, phase_response, channel, l, last_order, order, ir, m, &
!$OMP         phase_l, angular, angular_gradient, a, powers, factor, axis, reduced, jl, gradient, radial, jnext, bessel_values)
      DO ig = 1, SIZE(gvecs, 2)
         norm_g = SQRT(SUM(gvecs(:, ig)**2))
         bessel_ready = .FALSE.
         direction = 0.0_dp
         IF (norm_g > 0.0_dp) direction = gvecs(:, ig)/norm_g
         phase = fourpi/volume*EXP(CMPLX(0.0_dp, -DOT_PRODUCT(gvecs(:, ig), position), dp))
         IF (PRESENT(response)) THEN
            phase_response = CMPLX(-volume_response/volume, &
                                   -DOT_PRODUCT(gvecs_response(:, ig), position) - &
                                   DOT_PRODUCT(gvecs(:, ig), position_response), dp)
         END IF
         DO channel = 1, SIZE(angular_momenta)
            l = angular_momenta(channel)
            last_order = l
            IF (PRESENT(response) .AND. norm_g > 0.0_dp) last_order = l + 1
            ! Radial Bessel factors depend on l, not on the harmonic index m.
            DO order = l, last_order
               IF (bessel_ready(order)) CYCLE
               DO ir = 1, SIZE(radii)
                  bessel_values(ir, order) = bessel0(norm_g*radii(ir), order)
               END DO
               bessel_ready(order) = .TRUE.
            END DO
            m = harmonic_indices(channel) - l*l
            phase_l = phase*CMPLX(0.0_dp, -1.0_dp, dp)**l
            angular = 0.0_dp
            angular_gradient = 0.0_dp
            DO a = 1, nco(l)
               powers = indco(:, ncoset(l - 1) + a)
               factor = orbtramat(l)%slm(m, a)
               angular = angular + factor*PRODUCT(direction**powers)
               IF (.NOT. PRESENT(response)) CYCLE
               DO axis = 1, 3
                  IF (powers(axis) == 0) CYCLE
                  reduced = powers
                  reduced(axis) = reduced(axis) - 1
                  angular_gradient(axis) = angular_gradient(axis) + &
                                           factor*REAL(powers(axis), dp)*PRODUCT(direction**reduced)
               END DO
            END DO
            DO ir = 1, SIZE(radii)
               jl = bessel_values(ir, l)
               kernel(ig, ir, channel) = phase_l*angular*jl
               IF (.NOT. PRESENT(response)) CYCLE
               gradient = 0.0_dp
               radial = 0.0_dp
               IF (norm_g > 0.0_dp) THEN
                  jnext = bessel_values(ir, l + 1)
                  gradient = jl/norm_g*angular_gradient - radii(ir)*jnext*angular*direction
                  radial = angular*(l/radii(ir)*jl - norm_g*jnext)
               ELSE IF (l == 1) THEN
                  ! The origin limit contains dipole derivatives even though rho_G itself is zero.
                  gradient = radii(ir)/3.0_dp*angular_gradient
               END IF
               response(ig, ir, channel) = kernel(ig, ir, channel)*phase_response + &
                                           phase_l*(DOT_PRODUCT(gradient, gvecs_response(:, ig)) + &
                                                    radial*radii_response(ir))
            END DO
         END DO
      END DO
!$OMP END PARALLEL DO

      IF (.NOT. PRESENT(grid_shape)) RETURN
      IF (ALL(grid_shape == 0)) RETURN
      CPASSERT(ALL(grid_shape > 0))
      CPASSERT(PRESENT(cell))
      inverse_cell = cell
      CALL invmat(inverse_cell, info)
      CPASSERT(info == 0)
      alias_steps = twopi*TRANSPOSE(inverse_cell)
      DO axis = 1, 3
         alias_steps(:, axis) = alias_steps(:, axis)*REAL(grid_shape(axis), dp)
      END DO
      deformation = 0.0_dp
      IF (PRESENT(cell_response)) deformation = MATMUL(cell_response, inverse_cell)
      DO ig = 1, SIZE(gvecs, 2)
         lattice_index = MATMUL(TRANSPOSE(cell), gvecs(:, ig))/twopi
         miller = NINT(lattice_index)
         lattice_error = MAXVAL(ABS(lattice_index - REAL(miller, dp)))
         valid_lattice = lattice_error <= 128.0_dp*EPSILON(1.0_dp)*MAX(1.0_dp, MAXVAL(ABS(lattice_index)))
         CPASSERT(valid_lattice)
         CALL real_fourier_aliases(gvecs(:, ig), miller, grid_shape, alias_steps, vectors, count)
         IF (count == 1) CYCLE
         IF (.NOT. ALLOCATED(alias_kernel)) THEN
            ALLOCATE (alias_kernel(1, SIZE(radii), SIZE(angular_momenta)))
            IF (PRESENT(response)) ALLOCATE (alias_response, MOLD=alias_kernel)
         END IF
         kernel(ig, :, :) = kernel(ig, :, :)/REAL(count, dp)
         IF (PRESENT(response)) response(ig, :, :) = response(ig, :, :)/REAL(count, dp)
         DO alias = 2, count
            IF (PRESENT(response)) THEN
               ! The reciprocal alias jump follows the cell deformation.
               g_response(:, 1) = gvecs_response(:, ig) - &
                                  MATMUL(TRANSPOSE(deformation), vectors(:, alias) - gvecs(:, ig))
               CALL upf_inverse_coulomb_fourier_kernel(radii, angular_momenta, harmonic_indices, &
                                                       vectors(:, alias:alias), position, volume, alias_kernel, alias_response, &
                                                       radii_response, g_response, position_response, volume_response)
               response(ig, :, :) = response(ig, :, :) + alias_response(1, :, :)/REAL(count, dp)
            ELSE
               CALL upf_inverse_coulomb_fourier_kernel(radii, angular_momenta, harmonic_indices, &
                                                       vectors(:, alias:alias), position, volume, alias_kernel)
            END IF
            kernel(ig, :, :) = kernel(ig, :, :) + alias_kernel(1, :, :)/REAL(count, dp)
         END DO
      END DO
   END SUBROUTINE upf_inverse_coulomb_fourier_kernel

! **************************************************************************************************
!> \brief Apply the symmetric Coulomb kernel to real, orthonormal harmonic coefficients.
!> \param radii Positive radial nodes, in arbitrary order.
!> \param weights Radial measure including r**2 dr, without angular weights.
!> \param angular_momenta l of each real harmonic channel. Repeated l values are allowed.
!> \param density Signed radial density coefficients, indexed (radius,harmonic).
!> \param potential Unweighted Coulomb potential in the same harmonic order.
!> \param energy One half of the density/potential contraction with the same radial measure.
! **************************************************************************************************
   SUBROUTINE upf_inverse_coulomb_radial(radii, weights, angular_momenta, density, potential, energy)
      REAL(dp), INTENT(IN)                               :: radii(:), weights(:)
      INTEGER, INTENT(IN)                                :: angular_momenta(:)
      REAL(dp), INTENT(IN)                               :: density(:, :)
      REAL(dp), INTENT(OUT)                              :: potential(:, :), energy

      INTEGER                                            :: channel, i, j, l
      REAL(dp)                                           :: kernel

      CALL check_radial(radii, weights, angular_momenta, density)
      CPASSERT(ALL(SHAPE(potential) == SHAPE(density)))
      potential = 0.0_dp
      DO channel = 1, SIZE(angular_momenta)
         l = angular_momenta(channel)
         DO i = 1, SIZE(radii)
            DO j = 1, SIZE(radii)
               kernel = fourpi/REAL(2*l + 1, dp)/MAX(radii(i), radii(j))* &
                        (MIN(radii(i), radii(j))/MAX(radii(i), radii(j)))**l
               potential(i, channel) = potential(i, channel) + kernel*weights(j)*density(j, channel)
            END DO
         END DO
      END DO
      energy = 0.5_dp*SUM(SPREAD(weights, 2, SIZE(angular_momenta))*density*potential)
   END SUBROUTINE upf_inverse_coulomb_radial

! **************************************************************************************************
!> \brief Differentiate the same finite Coulomb energy, including the measure and kernel.
!>        Equal nodes must have equal radial motion so their common shell remains differentiable.
!> \param radii ...
!> \param weights ...
!> \param angular_momenta ...
!> \param density ...
!> \param radii_response Radial node motion. Zero for Cartesian grids moving with their centers.
!> \param weights_response Complete radial measure response.
!> \param density_response Signed density response, including the coupled inverse coefficients.
!> \param response Coulomb energy response.
! **************************************************************************************************
   SUBROUTINE upf_inverse_coulomb_radial_response(radii, weights, angular_momenta, density, &
                                                  radii_response, weights_response, density_response, response)
      REAL(dp), INTENT(IN)                               :: radii(:), weights(:)
      INTEGER, INTENT(IN)                                :: angular_momenta(:)
      REAL(dp), INTENT(IN)                               :: density(:, :), radii_response(:), &
                                                            weights_response(:), &
                                                            density_response(:, :)
      REAL(dp), INTENT(OUT)                              :: response

      INTEGER                                            :: channel, i, j, l, large, small
      REAL(dp)                                           :: derivative, energy, kernel
      REAL(dp), ALLOCATABLE                              :: potential(:, :)

      CPASSERT(SIZE(radii_response) == SIZE(radii))
      CPASSERT(SIZE(weights_response) == SIZE(weights))
      CPASSERT(ALL(SHAPE(density_response) == SHAPE(density)))
      ALLOCATE (potential(SIZE(density, 1), SIZE(density, 2)))
      CALL upf_inverse_coulomb_radial(radii, weights, angular_momenta, density, potential, energy)
      response = SUM((SPREAD(weights, 2, SIZE(angular_momenta))*density_response + &
                      SPREAD(weights_response, 2, SIZE(angular_momenta))*density)*potential)
      DO channel = 1, SIZE(angular_momenta)
         l = angular_momenta(channel)
         DO i = 1, SIZE(radii)
            DO j = 1, SIZE(radii)
               small = i
               large = j
               IF (radii(i) > radii(j)) THEN
                  small = j
                  large = i
               END IF
               IF (radii(i) == radii(j)) THEN
                  CPASSERT(radii_response(i) == radii_response(j))
               END IF
               kernel = fourpi/REAL(2*l + 1, dp)/radii(large)*(radii(small)/radii(large))**l
               derivative = kernel*(l*radii_response(small)/radii(small) - &
                                    (l + 1)*radii_response(large)/radii(large))
               response = response + 0.5_dp*weights(i)*weights(j)* &
                          density(i, channel)*density(j, channel)*derivative
            END DO
         END DO
      END DO
   END SUBROUTINE upf_inverse_coulomb_radial_response

! **************************************************************************************************
!> \brief Match each multipole with a supplied radial shape in the defining finite quadrature.
!>        The signed remainder has zero exterior moment in every represented harmonic channel.
!>        This identity applies to this radial representation, not to discarded angular channels.
!> \param radii ...
!> \param weights ...
!> \param angular_momenta ...
!> \param density ...
!> \param shapes Positive radial compensation shapes, indexed (radius,harmonic).
!> \param compensation Moment-matched compensation density in the same radial representation.
!> \param residual Signed density minus compensation.
!> \param moments Original exterior multipole moments.
! **************************************************************************************************
   SUBROUTINE upf_inverse_coulomb_compensate(radii, weights, angular_momenta, density, shapes, &
                                             compensation, residual, moments)
      REAL(dp), INTENT(IN)                               :: radii(:), weights(:)
      INTEGER, INTENT(IN)                                :: angular_momenta(:)
      REAL(dp), INTENT(IN)                               :: density(:, :), shapes(:, :)
      REAL(dp), INTENT(OUT)                              :: compensation(:, :), residual(:, :), &
                                                            moments(:)

      INTEGER                                            :: channel
      REAL(dp)                                           :: norm

      CALL check_radial(radii, weights, angular_momenta, density)
      CPASSERT(ALL(SHAPE(shapes) == SHAPE(density)))
      CPASSERT(ALL(SHAPE(compensation) == SHAPE(density)))
      CPASSERT(ALL(SHAPE(residual) == SHAPE(density)))
      CPASSERT(SIZE(moments) == SIZE(angular_momenta))
      CPASSERT(ALL(ieee_is_finite(shapes)) .AND. ALL(shapes >= 0.0_dp))
      DO channel = 1, SIZE(angular_momenta)
         norm = SUM(weights*radii**angular_momenta(channel)*shapes(:, channel))
         CPASSERT(norm > 0.0_dp)
         moments(channel) = SUM(weights*radii**angular_momenta(channel)*density(:, channel))
         compensation(:, channel) = moments(channel)/norm*shapes(:, channel)
      END DO
      residual = density - compensation
   END SUBROUTINE upf_inverse_coulomb_compensate

! **************************************************************************************************
!> \brief Apply the radial-measure adjoint of density minus its moment compensation.
!>        Compensation depends on the occupation density and contributes to the AO operator.
!> \param radii ...
!> \param weights ...
!> \param angular_momenta ...
!> \param shapes Same fixed radial compensation shapes as in the forward map.
!> \param potential Unweighted covector acting on the signed remainder.
!> \param adjoint Unweighted covector acting on the original density.
! **************************************************************************************************
   SUBROUTINE upf_inverse_coulomb_compensation_adjoint(radii, weights, angular_momenta, shapes, potential, adjoint)
      REAL(dp), INTENT(IN)                               :: radii(:), weights(:)
      INTEGER, INTENT(IN)                                :: angular_momenta(:)
      REAL(dp), INTENT(IN)                               :: shapes(:, :), potential(:, :)
      REAL(dp), INTENT(OUT)                              :: adjoint(:, :)

      INTEGER                                            :: channel
      REAL(dp)                                           :: norm

      CALL check_radial(radii, weights, angular_momenta, potential)
      CPASSERT(ALL(SHAPE(shapes) == SHAPE(potential)))
      CPASSERT(ALL(SHAPE(adjoint) == SHAPE(potential)))
      CPASSERT(ALL(ieee_is_finite(shapes)) .AND. ALL(shapes >= 0.0_dp))
      DO channel = 1, SIZE(angular_momenta)
         norm = SUM(weights*radii**angular_momenta(channel)*shapes(:, channel))
         CPASSERT(norm > 0.0_dp)
         adjoint(:, channel) = potential(:, channel) - radii**angular_momenta(channel)* &
                               SUM(weights*shapes(:, channel)*potential(:, channel))/norm
      END DO
   END SUBROUTINE upf_inverse_coulomb_compensation_adjoint

! **************************************************************************************************
!> \brief Validate the common real radial representation.
!> \param radii ...
!> \param weights ...
!> \param angular_momenta ...
!> \param density ...
! **************************************************************************************************
   SUBROUTINE check_radial(radii, weights, angular_momenta, density)
      REAL(dp), INTENT(IN)                               :: radii(:), weights(:)
      INTEGER, INTENT(IN)                                :: angular_momenta(:)
      REAL(dp), INTENT(IN)                               :: density(:, :)

      CPASSERT(SIZE(weights) == SIZE(radii))
      CPASSERT(SIZE(density, 1) == SIZE(radii))
      CPASSERT(SIZE(density, 2) == SIZE(angular_momenta))
      CPASSERT(ALL(ieee_is_finite(radii)) .AND. ALL(radii > 0.0_dp))
      CPASSERT(ALL(ieee_is_finite(weights)) .AND. ALL(weights >= 0.0_dp))
      CPASSERT(ALL(ieee_is_finite(density)))
      CPASSERT(ALL(angular_momenta >= 0))
   END SUBROUTINE check_radial

END MODULE upf_inverse_coulomb
