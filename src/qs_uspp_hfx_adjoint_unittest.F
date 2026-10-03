!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM qs_uspp_hfx_adjoint_unittest
   USE dbt_api,                         ONLY: &
        dbt_create, dbt_destroy, dbt_distribution_destroy, dbt_distribution_new, &
        dbt_distribution_type, dbt_finalize, dbt_get_block, dbt_get_stored_coordinates, &
        dbt_pgrid_create, dbt_pgrid_destroy, dbt_pgrid_type, dbt_put_block, dbt_type
   USE kinds,                           ONLY: dp
   USE libcp2k,                         ONLY: cp2k_finalize,&
                                              cp2k_init
   USE message_passing,                 ONLY: mp_comm_type,&
                                              mp_comm_world
   USE qs_tensors_types,                ONLY: create_3c_tensor
   USE qs_uspp_hfx,                     ONLY: qs_uspp_hfx_collect_adjoint,&
                                              qs_uspp_hfx_unstack_adjoint

   IMPLICIT NONE

   INTEGER, PARAMETER :: ri_sizes(4) = [1, 2, 1, 2], ao_sizes(3) = [1, 2, 1], &
                         images(5) = [3, 0, 1, 3, 2], original_image(3) = [2, 3, 1]
   INTEGER :: dims(3), group_case, data_case, ngroups, group_size, igroup, ind(3), sizes(3), &
              owner, image, first, slot, global_slot, i, j, k, a, b, c, group, mapped_image
   INTEGER, ALLOCATABLE :: d1(:), d2(:), d3(:), stack_dist3(:)
   LOGICAL :: found
   REAL(dp) :: error, expected, max_error
   REAL(dp), ALLOCATABLE :: block(:, :, :)
   TYPE(dbt_type) :: target(3), source(3), stack
   TYPE(dbt_distribution_type) :: stack_distribution
   TYPE(dbt_pgrid_type) :: world_grid, subgroup_grid
   TYPE(mp_comm_type) :: world, subgroup

   CALL cp2k_init()
   CALL world%set_handle(mp_comm_world%get_handle())
   max_error = 0.0_dp
   dims = [world%num_pe, 1, 1]
   CALL dbt_pgrid_create(world, dims, world_grid)
   DO group_case = 1, 3
      ngroups = 1
      IF (group_case == 2) ngroups = MIN(2, world%num_pe)
      IF (group_case == 3) ngroups = world%num_pe
      IF (group_case == 2 .AND. ngroups == 1) CYCLE
      IF (group_case == 3 .AND. ngroups <= 2) CYCLE
      IF (MODULO(world%num_pe, ngroups) /= 0) CYCLE
      group_size = world%num_pe/ngroups
      igroup = world%mepos/group_size
      CALL subgroup%from_split(world, igroup)
      dims = [1, group_size, 1]
      CALL dbt_pgrid_create(subgroup, dims, subgroup_grid)
      CALL create_3c_tensor(source(1), d1, d2, d3, subgroup_grid, ri_sizes, ao_sizes, ao_sizes, &
                            [1], [2, 3], name='Subgroup image adjoint')
      ALLOCATE (stack_dist3(2*SIZE(d3)))
      stack_dist3(:) = [d3, d3]
      CALL dbt_distribution_new(stack_distribution, subgroup_grid, d1, d2, stack_dist3)
      CALL dbt_create(stack, 'Batched image adjoint', stack_distribution, [1], [2, 3], &
                      ri_sizes, ao_sizes, [ao_sizes, ao_sizes])
      CALL dbt_distribution_destroy(stack_distribution)
      DEALLOCATE (d1, d2, d3, stack_dist3)
      CALL create_3c_tensor(TARGET(1), d1, d2, d3, world_grid, ri_sizes, ao_sizes, ao_sizes, &
                            [1], [2, 3], name='World image adjoint')
      DEALLOCATE (d1, d2, d3)
      DO image = 2, 3
         CALL dbt_create(source(1), source(image))
         CALL dbt_create(TARGET(1), TARGET(image))
      END DO
      DO data_case = 1, 3
         ! Seed only selected world blocks: gathering must preserve them and create missing blocks.
         DO image = 1, 3
            DO k = 1, SIZE(ao_sizes)
               DO j = 1, SIZE(ao_sizes)
                  DO i = 1, SIZE(ri_sizes)
                     ind = [i, j, k]
                     IF (MODULO(SUM(ind) + image, 3) /= 0) CYCLE
                     CALL dbt_get_stored_coordinates(TARGET(image), ind, owner)
                     IF (owner /= world%mepos) CYCLE
                     sizes = [ri_sizes(i), ao_sizes(j), ao_sizes(k)]
                     ALLOCATE (block(sizes(1), sizes(2), sizes(3)))
                     block = 0.025_dp*image
                     CALL dbt_put_block(TARGET(image), ind, sizes, block)
                     DEALLOCATE (block)
                  END DO
               END DO
            END DO
            CALL dbt_finalize(TARGET(image))
         END DO
         DO first = 1, SIZE(images), 2
            DO slot = 1, MIN(2, SIZE(images) - first + 1)
               global_slot = first + slot - 1
               DO k = 1, SIZE(ao_sizes)
                  DO j = 1, SIZE(ao_sizes)
                     DO i = 1, SIZE(ri_sizes)
                        IF (.NOT. active(data_case, igroup, [i, j, k])) CYCLE
                        ind = [i, j, (slot - 1)*SIZE(ao_sizes) + k]
                        CALL dbt_get_stored_coordinates(stack, ind, owner)
                        IF (owner /= subgroup%mepos) CYCLE
                        sizes = [ri_sizes(i), ao_sizes(j), ao_sizes(k)]
                        ALLOCATE (block(sizes(1), sizes(2), sizes(3)))
                        DO c = 1, sizes(3)
                           DO b = 1, sizes(2)
                              DO a = 1, sizes(1)
                                 block(a, b, c) = value(igroup, global_slot, [i, j, k], [a, b, c])
                              END DO
                           END DO
                        END DO
                        CALL dbt_put_block(stack, ind, sizes, block)
                        DEALLOCATE (block)
                     END DO
                  END DO
               END DO
            END DO
            CALL dbt_finalize(stack)
            ! Repeated accumulation must apply each independent prefactor exactly once.
            CALL qs_uspp_hfx_unstack_adjoint(source, stack, images, first, -0.375_dp)
            CALL qs_uspp_hfx_unstack_adjoint(source, stack, images, first, 0.125_dp)
            CALL clear_stack()
         END DO
         DO image = 1, 3
            CALL dbt_finalize(source(image))
            CALL qs_uspp_hfx_collect_adjoint(TARGET(original_image(image)), source(image), world)
         END DO
         DO image = 1, 3
            DO k = 1, SIZE(ao_sizes)
               DO j = 1, SIZE(ao_sizes)
                  DO i = 1, SIZE(ri_sizes)
                     ind = [i, j, k]
                     CALL dbt_get_stored_coordinates(TARGET(image), ind, owner)
                     IF (owner /= world%mepos) CYCLE
                     CALL dbt_get_block(TARGET(image), ind, block, found)
                     DO c = 1, ao_sizes(k)
                        DO b = 1, ao_sizes(j)
                           DO a = 1, ri_sizes(i)
                              expected = 0.0_dp
                              IF (MODULO(SUM(ind) + image, 3) == 0) expected = 0.025_dp*image
                              DO group = 0, ngroups - 1
                                 IF (.NOT. active(data_case, group, ind)) CYCLE
                                 DO global_slot = 1, SIZE(images)
                                    mapped_image = images(global_slot)
                                    IF (mapped_image == 0) CYCLE
                                    IF (original_image(mapped_image) /= image) CYCLE
                                    expected = expected - 0.25_dp*value(group, global_slot, ind, [a, b, c])
                                 END DO
                              END DO
                              error = ABS(expected)
                              IF (found) error = ABS(block(a, b, c) - expected)
                              max_error = MAX(max_error, error)
                              IF (.NOT. (error < 2.e-13_dp)) ERROR STOP 'Subgroup image adjoint mismatch'
                           END DO
                        END DO
                     END DO
                     IF (ALLOCATED(block)) DEALLOCATE (block)
                  END DO
               END DO
            END DO
         END DO
         CALL clear_images()
      END DO
      DO image = 1, 3
         CALL dbt_destroy(source(image))
         CALL dbt_destroy(TARGET(image))
      END DO
      CALL dbt_destroy(stack)
      CALL dbt_pgrid_destroy(subgroup_grid)
      CALL subgroup%free()
   END DO
   CALL dbt_pgrid_destroy(world_grid)
   CALL world%max(max_error)
   IF (world%mepos == 0) WRITE (*, '(A,ES16.7)') 'Image adjoint transport maximum error: ', max_error
   CALL cp2k_finalize()

CONTAINS

! **************************************************************************************************
!> \brief Independent sparse support, including empty groups, ranks and globally empty tensors.
!> \param configuration ...
!> \param group_index ...
!> \param index ...
!> \return ...
! **************************************************************************************************
   LOGICAL FUNCTION active(configuration, group_index, index)
      INTEGER, INTENT(IN)                                :: configuration, group_index, INDEX(3)

      active = configuration == 1 .AND. MODULO(INDEX(1) + 2*INDEX(2) + INDEX(3) + group_index, 4) /= 0
      IF (configuration == 2) active = group_index == 0 .AND. INDEX(1) == 1
   END FUNCTION active

! **************************************************************************************************
!> \brief Distinct nonsymmetric entries, independent of tensor layout and communicator ownership.
!> \param group_index ...
!> \param image_slot ...
!> \param index ...
!> \param ENTRY ...
!> \return ...
! **************************************************************************************************
   REAL(dp) FUNCTION value(group_index, image_slot, index, ENTRY)
      INTEGER, INTENT(IN)                                :: group_index, image_slot, INDEX(3), &
                                                            ENTRY(3)

      value = 0.013_dp*REAL(1000*group_index + 100*image_slot + 10*INDEX(1) + INDEX(2) &
                            - 2*INDEX(3) + 3*ENTRY(1) - ENTRY(2) + 2*ENTRY(3), dp)
   END FUNCTION value

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE clear_stack()
      USE dbt_api, ONLY: dbt_clear
      CALL dbt_clear(stack)
   END SUBROUTINE clear_stack

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE clear_images()
      USE dbt_api, ONLY: dbt_clear
      INTEGER                                            :: image

      DO image = 1, 3
         CALL dbt_clear(source(image))
         CALL dbt_clear(TARGET(image))
      END DO
   END SUBROUTINE clear_images

END PROGRAM qs_uspp_hfx_adjoint_unittest
