!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Distributed native USPP transition densities and their RI-HFX geometry adjoint.
!>        Coulomb-metric RI-RHO uses the ORB basis for isolated Coulomb or image-resolved exchange.
!>        Complete AO pairs are retained: the global RI contraction includes intercenter QQ.
! **************************************************************************************************
MODULE qs_uspp_hfx
   USE admm_types,                      ONLY: admm_type
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_p_type,&
                                              gto_basis_set_type
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_dbcsr_api,                    ONLY: dbcsr_create,&
                                              dbcsr_get_info,&
                                              dbcsr_release,&
                                              dbcsr_type
   USE dbt_api,                         ONLY: &
        dbt_clear, dbt_contract, dbt_copy, dbt_copy_matrix_to_tensor, dbt_copy_tensor_to_matrix, &
        dbt_create, dbt_destroy, dbt_finalize, dbt_get_block, dbt_get_info, &
        dbt_get_stored_coordinates, dbt_iterator_blocks_left, dbt_iterator_next_block, &
        dbt_iterator_start, dbt_iterator_stop, dbt_iterator_type, dbt_pgrid_type, dbt_put_block, &
        dbt_scale, dbt_type
   USE hfx_types,                       ONLY: hfx_ri_type
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE input_constants,                 ONLY: do_potential_coulomb,&
                                              do_potential_id,&
                                              do_potential_short,&
                                              do_potential_truncated
   USE input_cp2k_hfx,                  ONLY: ri_pmat
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_type
   USE lebedev,                         ONLY: get_number_of_lebedev_grid,&
                                              lebedev_grid
   USE libint_2c_3c,                    ONLY: libint_potential_type
   USE mathconstants,                   ONLY: fourpi
   USE message_passing,                 ONLY: mp_comm_type
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_integral_utils,               ONLY: basis_set_list_setup
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_tensors_types,                ONLY: create_3c_tensor
   USE qs_uspp_projectors,              ONLY: qs_uspp_projector_cells,&
                                              qs_uspp_projector_gradient,&
                                              qs_uspp_projector_matrix
   USE qs_uspp_types,                   ONLY: qs_uspp_prepare_paw_reconstruction,&
                                              qs_uspp_type
   USE sap_kind_types,                  ONLY: build_sap_ints_basis,&
                                              release_sap_int,&
                                              sap_int_type,&
                                              sap_sort
   USE uspp_augmentation,               ONLY: uspp_augmentation_type
   USE uspp_hfx_integrals,              ONLY: uspp_hfx_auxiliary_integrals,&
                                              uspp_hfx_integrals_init
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_uspp_hfx'
   PUBLIC :: qs_uspp_hfx_add_tensor, qs_uspp_hfx_gradient
   PUBLIC :: qs_uspp_hfx_transition_tensor, qs_uspp_hfx_transition_adjoint
   PUBLIC :: qs_uspp_hfx_native_tensor, qs_uspp_hfx_native_gradient
   PUBLIC :: qs_uspp_hfx_lattice_tensor, qs_uspp_hfx_lattice_gradient
   PUBLIC :: qs_uspp_hfx_add_kp_tensor, qs_uspp_hfx_kp_gradient
   PUBLIC :: qs_uspp_hfx_unstack_adjoint, qs_uspp_hfx_collect_adjoint

CONTAINS

! **************************************************************************************************
!> \brief Sum a batched image adjoint into its original image tensors on the same subgroup.
!> \param tensors Image adjoints with the repeated third-index distribution of stack.
!> \param stack Adjoint stacked along the third AO dimension.
!> \param images Image indices used by the forward stacking operation; zero denotes an absent image.
!> \param first First position in images represented by this batch.
!> \param pref Exchange and spin prefactor, included exactly once in the accumulated adjoint.
! **************************************************************************************************
   SUBROUTINE qs_uspp_hfx_unstack_adjoint(tensors, stack, images, first, pref)
      TYPE(dbt_type), INTENT(INOUT)                      :: tensors(:), stack
      INTEGER, INTENT(IN)                                :: images(:), first
      REAL(dp), INTENT(IN)                               :: pref

      INTEGER                                            :: image, ind(3), nblks(3), slot, TARGET(3)
      LOGICAL                                            :: found
      REAL(dp), ALLOCATABLE                              :: block(:, :, :)
      TYPE(dbt_iterator_type)                            :: iter

      CPASSERT(SIZE(tensors) > 0 .AND. first > 0)
      CALL dbt_get_info(tensors(1), nblks_total=nblks)
      CALL dbt_iterator_start(iter, stack)
      DO WHILE (dbt_iterator_blocks_left(iter))
         CALL dbt_iterator_next_block(iter, ind)
         slot = (ind(3) - 1)/nblks(3) + 1
         CPASSERT(first + slot - 1 <= SIZE(images))
         image = images(first + slot - 1)
         IF (image == 0) CYCLE
         CPASSERT(image > 0 .AND. image <= SIZE(tensors))
         TARGET = [ind(1), ind(2), MODULO(ind(3) - 1, nblks(3)) + 1]
         CALL dbt_get_block(stack, ind, block, found)
         CPASSERT(found)
         block(:, :, :) = pref*block(:, :, :)
         CALL dbt_put_block(tensors(image), TARGET, SHAPE(block), block, summation=.TRUE.)
         DEALLOCATE (block)
      END DO
      CALL dbt_iterator_stop(iter)
   END SUBROUTINE qs_uspp_hfx_unstack_adjoint

! **************************************************************************************************
!> \brief Sum subgroup-local image adjoints onto a world tensor without a global dense block map.
!> \param world_tensor World tensor, incremented with every independently computed contribution.
!> \param source Tensor on the caller's subgroup, with the same block sizes and index order.
!> \param comm World communicator shared by all subgroup ranks; call after their local loops finish.
! **************************************************************************************************
   SUBROUTINE qs_uspp_hfx_collect_adjoint(world_tensor, source, comm)
      TYPE(dbt_type), INTENT(INOUT)                      :: world_tensor, source

      CLASS(mp_comm_type), INTENT(IN)                    :: comm

      INTEGER                                            :: owner, iproc, i, ind(3), sizes(3), nblks(3), &
                                                            position, offset, length, nproc
      INTEGER, ALLOCATABLE :: send_blocks(:), recv_blocks(:), send_values(:), recv_values(:), &
                              send_ind_offset(:), recv_ind_offset(:), send_value_offset(:), recv_value_offset(:), &
                              next_ind(:), next_value(:), send_ind(:), recv_ind(:), b1(:), b2(:), b3(:)
      LOGICAL                                            :: found
      REAL(dp), ALLOCATABLE                              :: send_data(:), recv_data(:), block(:, :, :)
      TYPE(dbt_iterator_type)                            :: iter

      nproc = comm%num_pe
      ALLOCATE (send_blocks(0:nproc - 1), recv_blocks(0:nproc - 1), &
                send_values(0:nproc - 1), recv_values(0:nproc - 1), &
                send_ind_offset(0:nproc - 1), recv_ind_offset(0:nproc - 1), &
                send_value_offset(0:nproc - 1), recv_value_offset(0:nproc - 1), &
                next_ind(0:nproc - 1), next_value(0:nproc - 1))
      send_blocks = 0
      send_values = 0
      CALL dbt_iterator_start(iter, source)
      DO WHILE (dbt_iterator_blocks_left(iter))
         CALL dbt_iterator_next_block(iter, ind)
         CALL dbt_get_stored_coordinates(world_tensor, ind, owner)
         CALL dbt_get_block(source, ind, block, found)
         CPASSERT(found)
         CPASSERT(INT(send_values(owner), int_8) + SIZE(block, KIND=int_8) <= INT(HUGE(0), int_8))
         send_blocks(owner) = send_blocks(owner) + 1
         send_values(owner) = send_values(owner) + SIZE(block)
         DEALLOCATE (block)
      END DO
      CALL dbt_iterator_stop(iter)
      CALL comm%alltoall(send_blocks, recv_blocks, 1)
      CALL comm%alltoall(send_values, recv_values, 1)
      CPASSERT(SUM(INT(send_blocks, int_8)) <= INT(HUGE(0)/3, int_8))
      CPASSERT(SUM(INT(recv_blocks, int_8)) <= INT(HUGE(0)/3, int_8))
      CPASSERT(SUM(INT(send_values, int_8)) <= INT(HUGE(0), int_8))
      CPASSERT(SUM(INT(recv_values, int_8)) <= INT(HUGE(0), int_8))
      send_ind_offset(0) = 0
      recv_ind_offset(0) = 0
      send_value_offset(0) = 0
      recv_value_offset(0) = 0
      DO iproc = 1, nproc - 1
         send_ind_offset(iproc) = send_ind_offset(iproc - 1) + 3*send_blocks(iproc - 1)
         recv_ind_offset(iproc) = recv_ind_offset(iproc - 1) + 3*recv_blocks(iproc - 1)
         send_value_offset(iproc) = send_value_offset(iproc - 1) + send_values(iproc - 1)
         recv_value_offset(iproc) = recv_value_offset(iproc - 1) + recv_values(iproc - 1)
      END DO
      ALLOCATE (send_ind(3*SUM(send_blocks)), recv_ind(3*SUM(recv_blocks)), &
                send_data(SUM(send_values)), recv_data(SUM(recv_values)))
      next_ind(:) = send_ind_offset(:)
      next_value(:) = send_value_offset(:)
      CALL dbt_iterator_start(iter, source)
      DO WHILE (dbt_iterator_blocks_left(iter))
         CALL dbt_iterator_next_block(iter, ind)
         CALL dbt_get_stored_coordinates(world_tensor, ind, owner)
         CALL dbt_get_block(source, ind, block, found)
         CPASSERT(found)
         position = next_ind(owner)
         send_ind(position + 1:position + 3) = ind
         length = SIZE(block)
         offset = next_value(owner)
         send_data(offset + 1:offset + length) = RESHAPE(block, [length])
         next_ind(owner) = position + 3
         next_value(owner) = offset + length
         DEALLOCATE (block)
      END DO
      CALL dbt_iterator_stop(iter)
      CALL comm%alltoall(send_ind, 3*send_blocks, send_ind_offset, recv_ind, 3*recv_blocks, recv_ind_offset)
      CALL comm%alltoall(send_data, send_values, send_value_offset, recv_data, recv_values, recv_value_offset)
      CALL dbt_get_info(world_tensor, nblks_total=nblks)
      ALLOCATE (b1(nblks(1)), b2(nblks(2)), b3(nblks(3)))
      CALL dbt_get_info(world_tensor, blk_size_1=b1, blk_size_2=b2, blk_size_3=b3)
      DO iproc = 0, nproc - 1
         offset = recv_value_offset(iproc)
         DO i = 1, recv_blocks(iproc)
            position = recv_ind_offset(iproc) + 3*(i - 1)
            ind = recv_ind(position + 1:position + 3)
            sizes = [b1(ind(1)), b2(ind(2)), b3(ind(3))]
            length = PRODUCT(sizes)
            CALL dbt_put_block(world_tensor, ind, sizes, RESHAPE(recv_data(offset + 1:offset + length), sizes), &
                               summation=.TRUE.)
            offset = offset + length
         END DO
         CPASSERT(offset == recv_value_offset(iproc) + recv_values(iproc))
      END DO
      CALL dbt_finalize(world_tensor)
   END SUBROUTINE qs_uspp_hfx_collect_adjoint

! **************************************************************************************************
!> \brief Add full native augmentation after GG AO desymmetrization, before combined filtering.
!> \param qs_env Prepared Quickstep environment; an absent USPP state is a no-op.
!> \param ri_data Ordinary ORB/Coulomb RI-RHO configuration and process grid.
!> \param tensor Full (RI,AO,AO) tensor, possibly with split atomic blocks, incremented in place.
! **************************************************************************************************
   SUBROUTINE qs_uspp_hfx_add_tensor(qs_env, ri_data, tensor)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(hfx_ri_type), INTENT(IN)                      :: ri_data
      TYPE(dbt_type), INTENT(INOUT)                      :: tensor

      INTEGER, ALLOCATABLE                               :: nang(:)
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbt_type)                                     :: auxiliary
      TYPE(gto_basis_set_p_type), ALLOCATABLE            :: basis(:)
      TYPE(particle_type), POINTER                       :: particles(:)
      TYPE(qs_kind_type), POINTER                        :: kinds(:)
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_uspp_type), POINTER                        :: uspp

      CALL setup(qs_env, ri_data, ks_env, uspp, kinds, particles, basis, nang)
      IF (.NOT. ASSOCIATED(uspp)) RETURN
      CALL get_qs_env(qs_env, cell=cell)
      CALL qs_uspp_hfx_native_tensor(uspp, basis, particles, nang, ri_data%pgrid, auxiliary, &
                                     cell=cell, ri_metric=ri_data%ri_metric, periodic_sum=ANY(cell%perd /= 0), reconstructed_paw=.TRUE.)
      CALL qs_uspp_hfx_transition_tensor(uspp%overlaps, auxiliary, tensor)
      CALL dbt_destroy(auxiliary)
   END SUBROUTINE qs_uspp_hfx_add_tensor

! **************************************************************************************************
!> \brief Exact native geometry pullback of one complete first-AO batch of the HFX adjoint.
!>        Both B legs are differentiated separately; a batch need not be AO-transpose symmetric.
!> \param qs_env Prepared environment, including original unweighted B and SAP neighbor lists.
!> \param ri_data Same configuration and basis used for the forward tensor.
!> \param tensor_adjoint Opposite full transition tensor after density and RI-metric contraction.
!>        Its support must include native contributions, rather than only GG derivative blocks.
!> \param pref Energy-adjoint prefactor, -exchange_fraction*spin_factor, applied exactly once.
!> \param atom_gradient Local positive energy gradients (3,natom), incremented without MPI reduction.
!> \param cell_gradient Optional local physical affine-strain energy gradient, incremented in place.
! **************************************************************************************************
   SUBROUTINE qs_uspp_hfx_gradient(qs_env, ri_data, tensor_adjoint, pref, atom_gradient, cell_gradient)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(hfx_ri_type), INTENT(IN)                      :: ri_data
      TYPE(dbt_type), INTENT(INOUT)                      :: tensor_adjoint
      REAL(dp), INTENT(IN)                               :: pref
      REAL(dp), INTENT(INOUT)                            :: atom_gradient(:, :)
      REAL(dp), INTENT(INOUT), OPTIONAL                  :: cell_gradient(3, 3)

      INTEGER, ALLOCATABLE                               :: nang(:)
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_type)                                   :: overlap_gradient
      TYPE(dbt_type)                                     :: auxiliary, auxiliary_gradient
      TYPE(gto_basis_set_p_type), ALLOCATABLE            :: basis(:)
      TYPE(neighbor_list_set_p_type), POINTER            :: sap_ppnl(:)
      TYPE(particle_type), POINTER                       :: particles(:)
      TYPE(qs_kind_type), POINTER                        :: kinds(:)
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_uspp_type), POINTER                        :: uspp
      TYPE(sap_int_type), POINTER                        :: sap(:)

      CALL setup(qs_env, ri_data, ks_env, uspp, kinds, particles, basis, nang)
      IF (.NOT. ASSOCIATED(uspp)) RETURN
      IF (pref == 0.0_dp) RETURN
      CALL get_qs_env(qs_env, cell=cell)
      CALL qs_uspp_hfx_native_tensor(uspp, basis, particles, nang, ri_data%pgrid, auxiliary, &
                                     cell=cell, ri_metric=ri_data%ri_metric, periodic_sum=ANY(cell%perd /= 0), reconstructed_paw=.TRUE.)
      CALL qs_uspp_hfx_transition_adjoint(uspp%overlaps, auxiliary, tensor_adjoint, pref, &
                                          overlap_gradient, auxiliary_gradient)
      CALL projector_state(qs_env, ks_env, ri_data%orb_basis_type, uspp, sap_ppnl)
      CPASSERT(ASSOCIATED(sap_ppnl))
      ALLOCATE (sap(SIZE(kinds)**2))
      CALL build_sap_ints_basis(sap, sap_ppnl, kinds, ri_data%orb_basis_type, nder=1)
      CALL sap_sort(sap)
      CALL qs_uspp_projector_gradient(overlap_gradient, sap, uspp%center_atoms, atom_gradient, cell_gradient)
      CALL qs_uspp_hfx_native_gradient(uspp, basis, particles, nang, auxiliary_gradient, &
                                       atom_gradient, cell_gradient, cell=cell, ri_metric=ri_data%ri_metric, &
                                       periodic_sum=ANY(cell%perd /= 0), reconstructed_paw=.TRUE.)
      CALL release_sap_int(sap)
      CALL dbcsr_release(overlap_gradient)
      CALL dbt_destroy(auxiliary_gradient)
      CALL dbt_destroy(auxiliary)
   END SUBROUTINE qs_uspp_hfx_gradient

! **************************************************************************************************
!> \brief Add native image sums before the periodic RI metric and image reordering.
!> \param qs_env Prepared image-resolved projector environment.
!> \param ri_data Coulomb RI-RHO configuration with original RI and AO image maps.
!> \param tensors Full (stacked RI,AO,AO) tensors in the original AO-image order.
! **************************************************************************************************
   SUBROUTINE qs_uspp_hfx_add_kp_tensor(qs_env, ri_data, tensors)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(hfx_ri_type), INTENT(IN)                      :: ri_data
      TYPE(dbt_type), INTENT(INOUT)                      :: tensors(:)

      INTEGER, ALLOCATABLE                               :: ao_cells(:, :), nang(:), ri_cells(:, :)
      TYPE(cell_type), POINTER                           :: cell
      TYPE(gto_basis_set_p_type), ALLOCATABLE            :: basis(:)
      TYPE(particle_type), POINTER                       :: particles(:)
      TYPE(qs_kind_type), POINTER                        :: kinds(:)
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_uspp_type), POINTER                        :: uspp
      TYPE(sap_int_type), POINTER                        :: sap(:)

      CALL setup(qs_env, ri_data, ks_env, uspp, kinds, particles, basis, nang, with_images=.TRUE.)
      IF (.NOT. ASSOCIATED(uspp)) RETURN
      CALL prepare_lattice(qs_env, ks_env, ri_data, SIZE(kinds), 0, cell, sap, ri_cells, ao_cells)
      CALL qs_uspp_hfx_lattice_tensor(uspp%overlaps, sap, uspp, basis, particles, nang, ri_data%pgrid, &
                                      cell%hmat, ri_cells, ao_cells, tensors, cell=cell, reconstructed_paw=.TRUE.)
      CALL release_sap_int(sap)
   END SUBROUTINE qs_uspp_hfx_add_kp_tensor

! **************************************************************************************************
!> \brief Pull back complete world-communicator image adjoints to atoms and physical affine strain.
!> \param qs_env Same prepared environment as the forward map.
!> \param ri_data Same Coulomb RI-RHO configuration and image order.
!> \param tensor_adjoint Complete image adjoints after metric/bump/density contractions.
!> \param atom_gradient Local positive energy gradients, incremented without MPI reduction.
!> \param cell_gradient Optional physical affine-strain gradient, incremented without reduction.
!>        The image adjoints already contain their exchange/spin prefactor.
! **************************************************************************************************
   SUBROUTINE qs_uspp_hfx_kp_gradient(qs_env, ri_data, tensor_adjoint, atom_gradient, cell_gradient)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(hfx_ri_type), INTENT(IN)                      :: ri_data
      TYPE(dbt_type), INTENT(INOUT)                      :: tensor_adjoint(:)
      REAL(dp), INTENT(INOUT)                            :: atom_gradient(:, :)
      REAL(dp), INTENT(INOUT), OPTIONAL                  :: cell_gradient(3, 3)

      INTEGER, ALLOCATABLE                               :: ao_cells(:, :), nang(:), ri_cells(:, :)
      TYPE(cell_type), POINTER                           :: cell
      TYPE(gto_basis_set_p_type), ALLOCATABLE            :: basis(:)
      TYPE(particle_type), POINTER                       :: particles(:)
      TYPE(qs_kind_type), POINTER                        :: kinds(:)
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_uspp_type), POINTER                        :: uspp
      TYPE(sap_int_type), POINTER                        :: sap(:)

      CALL setup(qs_env, ri_data, ks_env, uspp, kinds, particles, basis, nang, with_images=.TRUE.)
      IF (.NOT. ASSOCIATED(uspp)) RETURN
      CALL prepare_lattice(qs_env, ks_env, ri_data, SIZE(kinds), 1, cell, sap, ri_cells, ao_cells)
      CALL qs_uspp_hfx_lattice_gradient(uspp%overlaps, sap, uspp, basis, particles, nang, ri_data%pgrid, &
                                        cell%hmat, ri_cells, ao_cells, tensor_adjoint, 1.0_dp, &
                                        atom_gradient, cell_gradient, cell=cell, reconstructed_paw=.TRUE.)
      CALL release_sap_int(sap)
   END SUBROUTINE qs_uspp_hfx_kp_gradient

! **************************************************************************************************
!> \brief Preserve the original integer image maps and SAP coordinate-folding convention.
!> \param qs_env ...
!> \param ks_env ...
!> \param ri_data ...
!> \param nkind ...
!> \param nder ...
!> \param cell ...
!> \param sap ...
!> \param ri_cells ...
!> \param ao_cells ...
! **************************************************************************************************
   SUBROUTINE prepare_lattice(qs_env, ks_env, ri_data, nkind, nder, cell, sap, ri_cells, ao_cells)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(hfx_ri_type), INTENT(IN)                      :: ri_data
      INTEGER, INTENT(IN)                                :: nkind, nder
      TYPE(cell_type), POINTER                           :: cell
      TYPE(sap_int_type), POINTER                        :: sap(:)
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: ri_cells(:, :), ao_cells(:, :)

      INTEGER, POINTER                                   :: index_to_cell(:, :)
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(neighbor_list_set_p_type), POINTER            :: sap_ppnl(:)
      TYPE(qs_kind_type), POINTER                        :: kinds(:)
      TYPE(qs_uspp_type), POINTER                        :: uspp

      CALL get_qs_env(qs_env, cell=cell, kpoints=kpoints, qs_kind_set=kinds)
      CPASSERT(ASSOCIATED(kpoints))
      CALL get_kpoint_info(kpoints, index_to_cell=index_to_cell)
      CPASSERT(ri_data%nimg <= SIZE(index_to_cell, 2))
      CPASSERT(ALL(ri_data%RI_cell_to_img > 0))
      CPASSERT(ALL(ri_data%RI_cell_to_img <= SIZE(index_to_cell, 2)))
      ALLOCATE (ri_cells(3, ri_data%ncell_RI), ao_cells(3, ri_data%nimg))
      ri_cells(:, :) = index_to_cell(:, ri_data%RI_cell_to_img)
      ao_cells(:, :) = index_to_cell(:, 1:ri_data%nimg)
      CALL projector_state(qs_env, ks_env, ri_data%orb_basis_type, uspp, sap_ppnl)
      CPASSERT(ASSOCIATED(sap_ppnl))
      ALLOCATE (sap(nkind*nkind))
      CALL build_sap_ints_basis(sap, sap_ppnl, kinds, ri_data%orb_basis_type, nder=nder)
      CALL sap_sort(sap)
   END SUBROUTINE prepare_lattice

! **************************************************************************************************
!> \brief Assemble all original projector-image pairs into unfiltered AO image tensors.
!>        T_Q(S,R) = sum_T B(T) C(S,T) B(T-R)^T, before the local RI inverse and bump.
!> \param template Matrix defining the SAP row distribution and communicator.
!> \param sap Original unweighted image-resolved projector overlaps.
!> \param uspp Original augmentation data and projector block ordering.
!> \param basis Original RI Gaussian bases.
!> \param particles Current RI atoms.
!> \param nang Per-kind angular controls.
!> \param pgrid Native tensor process grid on the template communicator.
!> \param cell_matrix Cartesian lattice vectors, including affine cell motion.
!> \param ri_cells Integer RI cells S in the first-dimension stacking order.
!> \param ao_cells Integer AO cells R in tensors order; no reciprocal phase or k weight here.
!> \param tensors Full (RI,AO,AO) image tensors, incremented without a Gaussian support restriction.
!> \param ri_domain Optional fixed local RI support (stacked RI atom,first AO atom).
!> \param cell Optional cell folding convention shared with the SAP labels.
! **************************************************************************************************
   SUBROUTINE qs_uspp_hfx_lattice_tensor(template, sap, uspp, basis, particles, nang, pgrid, &
                                         cell_matrix, ri_cells, ao_cells, tensors, ri_domain, cell, reconstructed_paw)
      TYPE(dbcsr_type), INTENT(IN)                       :: template
      TYPE(sap_int_type), INTENT(IN)                     :: sap(:)
      TYPE(qs_uspp_type), INTENT(IN)                     :: uspp
      TYPE(gto_basis_set_p_type), INTENT(IN)             :: basis(:)
      TYPE(particle_type), INTENT(IN)                    :: particles(:)
      INTEGER, INTENT(IN)                                :: nang(:)
      TYPE(dbt_pgrid_type), INTENT(IN)                   :: pgrid
      REAL(dp), INTENT(IN)                               :: cell_matrix(3, 3)
      INTEGER, INTENT(IN)                                :: ri_cells(:, :), ao_cells(:, :)
      TYPE(dbt_type), INTENT(INOUT)                      :: tensors(:)
      LOGICAL, INTENT(IN), OPTIONAL                      :: ri_domain(:, :)
      TYPE(cell_type), OPTIONAL, POINTER                 :: cell
      LOGICAL, INTENT(IN), OPTIONAL :: reconstructed_paw

      CALL lattice_action(template, sap, uspp, basis, particles, nang, pgrid, cell_matrix, &
                          ri_cells, ao_cells, tensors, ri_domain, cell, reconstructed_paw=reconstructed_paw)
   END SUBROUTINE qs_uspp_hfx_lattice_tensor

! **************************************************************************************************
!> \brief Exact geometry transpose of the complete projector-image sum on compatible process grids.
!>        Apply the same local RI mask before either B leg or the native C pullback. The caller
!>        supplies a complete image adjoint after its metric/bump/density contractions and performs
!>        MPI force reduction. Subgroup adjoints must first be lifted to this communicator.
!> \param template Same SAP distribution as the forward assembly.
!> \param sap Original overlaps with nder=1.
!> \param uspp Same original augmentation data.
!> \param basis Same original RI bases.
!> \param particles Current RI atoms.
!> \param nang Same angular controls.
!> \param pgrid Same native process grid and communicator.
!> \param cell_matrix Same Cartesian lattice vectors.
!> \param ri_cells Same RI-cell stack S.
!> \param ao_cells Same AO images R.
!> \param tensor_adjoint Complete image adjoints, unchanged.
!> \param pref Energy prefactor, applied exactly once.
!> \param atom_gradient Local positive energy gradients, incremented without MPI reduction.
!> \param cell_gradient Optional physical affine-strain gradient, incremented without reduction.
!> \param ri_domain Same optional fixed local RI support.
!> \param cell Same optional atom folding convention.
! **************************************************************************************************
   SUBROUTINE qs_uspp_hfx_lattice_gradient(template, sap, uspp, basis, particles, nang, pgrid, &
                                           cell_matrix, ri_cells, ao_cells, tensor_adjoint, pref, &
                                           atom_gradient, cell_gradient, ri_domain, cell, reconstructed_paw)
      TYPE(dbcsr_type), INTENT(IN)                       :: template
      TYPE(sap_int_type), INTENT(IN)                     :: sap(:)
      TYPE(qs_uspp_type), INTENT(IN)                     :: uspp
      TYPE(gto_basis_set_p_type), INTENT(IN)             :: basis(:)
      TYPE(particle_type), INTENT(IN)                    :: particles(:)
      INTEGER, INTENT(IN)                                :: nang(:)
      TYPE(dbt_pgrid_type), INTENT(IN)                   :: pgrid
      REAL(dp), INTENT(IN)                               :: cell_matrix(3, 3)
      INTEGER, INTENT(IN)                                :: ri_cells(:, :), ao_cells(:, :)
      TYPE(dbt_type), INTENT(INOUT)                      :: tensor_adjoint(:)
      REAL(dp), INTENT(IN)                               :: pref
      REAL(dp), INTENT(INOUT)                            :: atom_gradient(:, :)
      REAL(dp), INTENT(INOUT), OPTIONAL                  :: cell_gradient(3, 3)
      LOGICAL, INTENT(IN), OPTIONAL                      :: ri_domain(:, :)
      TYPE(cell_type), OPTIONAL, POINTER                 :: cell
      LOGICAL, INTENT(IN), OPTIONAL :: reconstructed_paw

      CALL lattice_action(template, sap, uspp, basis, particles, nang, pgrid, cell_matrix, &
                          ri_cells, ao_cells, tensor_adjoint, ri_domain, cell, pref, atom_gradient, cell_gradient, reconstructed_paw)
   END SUBROUTINE qs_uspp_hfx_lattice_gradient

! **************************************************************************************************
!> \brief Shared image-pair traversal for the forward map and its geometry transpose.
!> \param template ...
!> \param sap ...
!> \param uspp ...
!> \param basis ...
!> \param particles ...
!> \param nang ...
!> \param pgrid ...
!> \param cell_matrix ...
!> \param ri_cells ...
!> \param ao_cells ...
!> \param tensors ...
!> \param ri_domain ...
!> \param cell ...
!> \param pref ...
!> \param atom_gradient ...
!> \param cell_gradient ...
! **************************************************************************************************
   SUBROUTINE lattice_action(template, sap, uspp, basis, particles, nang, pgrid, cell_matrix, &
                             ri_cells, ao_cells, tensors, ri_domain, cell, pref, atom_gradient, cell_gradient, reconstructed_paw)
      TYPE(dbcsr_type), INTENT(IN)                       :: template
      TYPE(sap_int_type), INTENT(IN)                     :: sap(:)
      TYPE(qs_uspp_type), INTENT(IN)                     :: uspp
      TYPE(gto_basis_set_p_type), INTENT(IN)             :: basis(:)
      TYPE(particle_type), INTENT(IN)                    :: particles(:)
      INTEGER, INTENT(IN)                                :: nang(:)
      TYPE(dbt_pgrid_type), INTENT(IN)                   :: pgrid
      REAL(dp), INTENT(IN)                               :: cell_matrix(3, 3)
      INTEGER, INTENT(IN)                                :: ri_cells(:, :), ao_cells(:, :)
      TYPE(dbt_type), INTENT(INOUT)                      :: tensors(:)
      LOGICAL, INTENT(IN), OPTIONAL                      :: ri_domain(:, :)
      TYPE(cell_type), OPTIONAL, POINTER                 :: cell
      REAL(dp), INTENT(IN), OPTIONAL                     :: pref
      REAL(dp), INTENT(INOUT), OPTIONAL                  :: atom_gradient(:, :), cell_gradient(3, 3)
      LOGICAL, INTENT(IN), OPTIONAL :: reconstructed_paw

      INTEGER                                            :: image, nblks(3), t, u
      INTEGER, ALLOCATABLE                               :: d1(:), d2(:), d3(:), pp_cells(:, :), &
                                                            sizes_ri(:)
      INTEGER, POINTER                                   :: sizes_ao(:)
      REAL(dp)                                           :: qshift(3)
      REAL(dp), ALLOCATABLE                              :: shifts(:, :)
      TYPE(dbcsr_type)                                   :: gradient_left, gradient_right
      TYPE(dbcsr_type), ALLOCATABLE                      :: overlaps(:)
      TYPE(dbt_type)                                     :: auxiliary, auxiliary_gradient, &
                                                            auxiliary_total, full

      CPASSERT(SIZE(ri_cells, 1) == 3 .AND. SIZE(ri_cells, 2) > 0)
      CPASSERT(SIZE(ao_cells, 1) == 3 .AND. SIZE(ao_cells, 2) == SIZE(tensors))
      CPASSERT(ALL(ieee_is_finite(cell_matrix)))
      DO image = 1, SIZE(ao_cells, 2)
         CPASSERT(cell_column(ao_cells, ao_cells(:, image)) == image)
      END DO
      DO image = 1, SIZE(ri_cells, 2)
         CPASSERT(cell_column(ri_cells, ri_cells(:, image)) == image)
      END DO
      CPASSERT(PRESENT(pref) .EQV. PRESENT(atom_gradient))
      IF (PRESENT(cell_gradient)) THEN
         CPASSERT(PRESENT(atom_gradient))
      END IF
      CALL dbcsr_get_info(template, row_blk_size=sizes_ao)
      CPASSERT(SIZE(sizes_ao) == SIZE(particles))
      IF (PRESENT(ri_domain)) THEN
         CPASSERT(SIZE(ri_domain, 1) == SIZE(particles)*SIZE(ri_cells, 2))
         CPASSERT(SIZE(ri_domain, 2) == SIZE(particles))
      END IF
      CALL qs_uspp_projector_cells(template, sap, uspp%center_atoms, pp_cells)
      ALLOCATE (overlaps(SIZE(pp_cells, 2)))
      DO t = 1, SIZE(pp_cells, 2)
         CALL qs_uspp_projector_matrix(template, sap, uspp%center_atoms, uspp%projector_sizes, &
                                       overlaps(t), pp_cells(:, t))
      END DO
      ALLOCATE (shifts(3, SIZE(ri_cells, 2)))
      shifts(:, :) = MATMUL(cell_matrix, REAL(ri_cells, dp))
      DO t = 1, SIZE(pp_cells, 2)
         qshift = MATMUL(cell_matrix, REAL(pp_cells(:, t), dp))
         CALL qs_uspp_hfx_native_tensor(uspp, basis, particles, nang, pgrid, auxiliary, shifts, qshift, cell, reconstructed_paw=reconstructed_paw)
         CALL dbt_get_info(auxiliary, nblks_total=nblks)
         ALLOCATE (sizes_ri(nblks(1)))
         CALL dbt_get_info(auxiliary, blk_size_1=sizes_ri)
         CALL create_3c_tensor(full, d1, d2, d3, pgrid, sizes_ri, sizes_ao, sizes_ao, &
                               [1], [2, 3], name='Native lattice transition')
         DEALLOCATE (d1, d2, d3, sizes_ri)
         IF (PRESENT(atom_gradient)) CALL dbt_create(auxiliary, auxiliary_total)
         DO u = 1, SIZE(pp_cells, 2)
            image = cell_column(ao_cells, pp_cells(:, t) - pp_cells(:, u))
            IF (image == 0) CYCLE
            CALL dbt_clear(full)
            IF (PRESENT(atom_gradient)) THEN
               CALL dbt_copy(tensors(image), full)
               CALL apply_ri_domain(full, ri_domain)
               CALL qs_uspp_hfx_transition_adjoint(overlaps(t), auxiliary, full, pref, &
                                                   gradient_left, auxiliary_gradient, overlaps(u), gradient_right)
               CALL qs_uspp_projector_gradient(gradient_left, sap, uspp%center_atoms, atom_gradient, &
                                               cell_gradient, pp_cells(:, t))
               CALL qs_uspp_projector_gradient(gradient_right, sap, uspp%center_atoms, atom_gradient, &
                                               cell_gradient, pp_cells(:, u))
               CALL dbt_copy(auxiliary_gradient, auxiliary_total, summation=.TRUE., move_data=.TRUE.)
               CALL dbcsr_release(gradient_left)
               CALL dbcsr_release(gradient_right)
               CALL dbt_destroy(auxiliary_gradient)
            ELSE
               CALL qs_uspp_hfx_transition_tensor(overlaps(t), auxiliary, full, overlaps(u))
               CALL apply_ri_domain(full, ri_domain)
               CALL dbt_copy(full, tensors(image), summation=.TRUE., move_data=.TRUE.)
            END IF
         END DO
         IF (PRESENT(atom_gradient)) THEN
            CALL qs_uspp_hfx_native_gradient(uspp, basis, particles, nang, auxiliary_total, atom_gradient, &
                                             cell_gradient, shifts, qshift, cell, reconstructed_paw=reconstructed_paw)
            CALL dbt_destroy(auxiliary_total)
         END IF
         CALL dbt_destroy(full)
         CALL dbt_destroy(auxiliary)
      END DO
      DO t = 1, SIZE(overlaps)
         CALL dbcsr_release(overlaps(t))
      END DO
   END SUBROUTINE lattice_action

! **************************************************************************************************
!> \brief Locate an exact integer image; absent images remain outside the caller's finite domain.
!> \param cells Distinct integer image columns.
!> \param requested_cell Requested integer cell.
!> \return Matching column, or zero.
! **************************************************************************************************
   FUNCTION cell_column(cells, requested_cell) RESULT(index)
      INTEGER, INTENT(IN)                                :: cells(:, :), requested_cell(3)
      INTEGER                                            :: index

      INTEGER                                            :: i

      index = 0
      DO i = 1, SIZE(cells, 2)
         IF (ALL(cells(:, i) == requested_cell)) THEN
            CPASSERT(index == 0)
            index = i
         END IF
      END DO
   END FUNCTION cell_column

! **************************************************************************************************
!> \brief Apply the same fixed local RI mask to an unsplit forward tensor or its adjoint.
!>        No Gaussian derivative or AO-pair support mask is applied. The smooth RI bump is
!>        part of the subsequent local metric contraction, with its own geometry derivative.
!> \param tensor Unsplit (RI atom/cell,AO atom,AO atom) tensor.
!> \param domain Optional fixed support; absence retains every native block.
! **************************************************************************************************
   SUBROUTINE apply_ri_domain(tensor, domain)
      TYPE(dbt_type), INTENT(INOUT)                      :: tensor
      LOGICAL, INTENT(IN), OPTIONAL                      :: domain(:, :)

      INTEGER                                            :: ind(3)
      LOGICAL                                            :: found
      REAL(dp), ALLOCATABLE                              :: block(:, :, :)
      TYPE(dbt_iterator_type)                            :: iter

      IF (.NOT. PRESENT(domain)) RETURN
      CALL dbt_iterator_start(iter, tensor)
      DO WHILE (dbt_iterator_blocks_left(iter))
         CALL dbt_iterator_next_block(iter, ind)
         CPASSERT(ind(1) <= SIZE(domain, 1) .AND. ind(2) <= SIZE(domain, 2))
         IF (domain(ind(1), ind(2))) CYCLE
         CALL dbt_get_block(tensor, ind, block, found)
         CPASSERT(found)
         block = 0.0_dp
         CALL dbt_put_block(tensor, ind, SHAPE(block), block)
         DEALLOCATE (block)
      END DO
      CALL dbt_iterator_stop(iter)
   END SUBROUTINE apply_ri_domain

! **************************************************************************************************
!> \brief Distributed T += B C B^T with C(RI,P,P), B(AO,P); no packing or symmetry factor.
!>        C has only center-diagonal P blocks. T keeps all AO pairs and all RI centers.
!> \param overlaps Original rectangular distributed B matrix; never globally replicated.
!> \param auxiliary Complete native C, with unsplit RI and magnetic-projector blocks.
!> \param tensor Target full (RI,AO,AO), incremented after reblocking and redistribution.
!> \param right_overlaps Optional independent right B, e.g. B(AO^0,P^(T-R)) for an AO image R.
!>        The default is the same B on both legs. No AO symmetry or phase is imposed.
! **************************************************************************************************
   SUBROUTINE qs_uspp_hfx_transition_tensor(overlaps, auxiliary, tensor, right_overlaps)
      TYPE(dbcsr_type), INTENT(IN), TARGET               :: overlaps
      TYPE(dbt_type), INTENT(INOUT)                      :: auxiliary, tensor
      TYPE(dbcsr_type), INTENT(IN), OPTIONAL             :: right_overlaps

      TYPE(dbt_type)                                     :: b, full, right_b, u

      CALL transition_templates(overlaps, auxiliary, b, u, full)
      CALL dbt_copy_matrix_to_tensor(overlaps, b)
      CALL right_transition_tensor(overlaps, right_overlaps, right_b)
      CALL dbt_contract(1.0_dp, auxiliary, right_b, 0.0_dp, u, &
                        contract_1=[3], notcontract_1=[1, 2], &
                        contract_2=[2], notcontract_2=[1], map_1=[1, 2], map_2=[3], filter_eps=0.0_dp)
      CALL dbt_contract(1.0_dp, b, u, 0.0_dp, full, &
                        contract_1=[2], notcontract_1=[1], &
                        contract_2=[2], notcontract_2=[1, 3], map_1=[2], map_2=[1, 3], filter_eps=0.0_dp)
      CALL dbt_copy(full, tensor, summation=.TRUE., move_data=.TRUE.)
      CALL dbt_destroy(full)
      CALL dbt_destroy(u)
      CALL dbt_destroy(b)
      CALL dbt_destroy(right_b)
   END SUBROUTINE qs_uspp_hfx_transition_tensor

! **************************************************************************************************
!> \brief Distributed transpose of exactly the two forward contractions, including both B legs.
!>        Native C is a center-diagonal variable; off-center C adjoints are not retained.
!> \param overlaps Original B and output block ownership template.
!> \param auxiliary Original native C, unchanged.
!> \param tensor_adjoint Full or first-AO-batched (RI,AO,AO) opposite tensor, unchanged.
!> \param pref Prefactor applied once before the reverse contractions.
!> \param overlap_gradient New dE/dB matrix on B's original rectangular distribution.
!> \param auxiliary_gradient New dE/dC on C's original center-diagonal distribution.
!> \param right_overlaps Optional independent right B, as in the forward transition.
!> \param right_overlap_gradient New right-leg dE/dB; present iff right_overlaps is present.
!>        Otherwise both legs are accumulated into overlap_gradient, as in the Gamma path.
! **************************************************************************************************
   SUBROUTINE qs_uspp_hfx_transition_adjoint(overlaps, auxiliary, tensor_adjoint, pref, &
                                             overlap_gradient, auxiliary_gradient, right_overlaps, right_overlap_gradient)
      TYPE(dbcsr_type), INTENT(IN), TARGET               :: overlaps
      TYPE(dbt_type), INTENT(INOUT)                      :: auxiliary, tensor_adjoint
      REAL(dp), INTENT(IN)                               :: pref
      TYPE(dbcsr_type), INTENT(OUT)                      :: overlap_gradient
      TYPE(dbt_type), INTENT(OUT)                        :: auxiliary_gradient
      TYPE(dbcsr_type), INTENT(IN), OPTIONAL             :: right_overlaps
      TYPE(dbcsr_type), INTENT(OUT), OPTIONAL            :: right_overlap_gradient

      TYPE(dbt_type)                                     :: b, bar_b, bar_right_b, bar_u, full, &
                                                            right_b, u

      CPASSERT(PRESENT(right_overlaps) .EQV. PRESENT(right_overlap_gradient))
      CALL transition_templates(overlaps, auxiliary, b, u, full)
      CALL dbt_copy_matrix_to_tensor(overlaps, b)
      CALL right_transition_tensor(overlaps, right_overlaps, right_b)
      CALL dbt_contract(1.0_dp, auxiliary, right_b, 0.0_dp, u, &
                        contract_1=[3], notcontract_1=[1, 2], &
                        contract_2=[2], notcontract_2=[1], map_1=[1, 2], map_2=[3], filter_eps=0.0_dp)
      CALL dbt_copy(tensor_adjoint, full)
      CALL dbt_scale(full, pref)
      CALL dbt_create(u, bar_u)
      CALL dbt_create(b, bar_b)
      CALL dbt_create(right_b, bar_right_b)
      CALL dbt_contract(1.0_dp, full, u, 0.0_dp, bar_b, &
                        contract_1=[1, 3], notcontract_1=[2], &
                        contract_2=[1, 3], notcontract_2=[2], map_1=[1], map_2=[2], filter_eps=0.0_dp)
      CALL dbt_contract(1.0_dp, full, b, 0.0_dp, bar_u, &
                        contract_1=[2], notcontract_1=[1, 3], &
                        contract_2=[1], notcontract_2=[2], map_1=[1, 3], map_2=[2], filter_eps=0.0_dp)
      CALL dbt_contract(1.0_dp, bar_u, auxiliary, 0.0_dp, bar_right_b, &
                        contract_1=[1, 2], notcontract_1=[3], &
                        contract_2=[1, 2], notcontract_2=[3], map_1=[1], map_2=[2], filter_eps=0.0_dp)
      CALL dbt_create(auxiliary, auxiliary_gradient)
      CALL dbt_copy(auxiliary, auxiliary_gradient)
      CALL dbt_contract(1.0_dp, bar_u, right_b, 0.0_dp, auxiliary_gradient, &
                        contract_1=[3], notcontract_1=[1, 2], &
                        contract_2=[1], notcontract_2=[2], map_1=[1, 2], map_2=[3], filter_eps=0.0_dp, &
                        retain_sparsity=.TRUE.)
      ! bar_b retains the same matrix-derived block layout and ownership as the original B.
      IF (PRESENT(right_overlaps)) THEN
         CALL dbcsr_create(right_overlap_gradient, template=right_overlaps)
         CALL dbt_copy_tensor_to_matrix(bar_right_b, right_overlap_gradient)
      ELSE
         CALL dbt_copy(bar_right_b, bar_b, summation=.TRUE.)
      END IF
      CALL dbcsr_create(overlap_gradient, template=overlaps)
      CALL dbt_copy_tensor_to_matrix(bar_b, overlap_gradient)
      CALL dbt_destroy(bar_b)
      CALL dbt_destroy(bar_u)
      CALL dbt_destroy(full)
      CALL dbt_destroy(u)
      CALL dbt_destroy(b)
      CALL dbt_destroy(bar_right_b)
      CALL dbt_destroy(right_b)
   END SUBROUTINE qs_uspp_hfx_transition_adjoint

! **************************************************************************************************
!> \brief Copy the selected right projector leg, retaining its own sparse matrix distribution.
!> \param overlaps Left B, defining the common AO and magnetic-projector block ordering.
!> \param right_overlaps Optional independent right B; defaults to the left B.
!> \param right_b New matrix-derived tensor of the selected right leg.
! **************************************************************************************************
   SUBROUTINE right_transition_tensor(overlaps, right_overlaps, right_b)
      TYPE(dbcsr_type), INTENT(IN)                       :: overlaps
      TYPE(dbcsr_type), INTENT(IN), OPTIONAL             :: right_overlaps
      TYPE(dbt_type), INTENT(OUT)                        :: right_b

      INTEGER, POINTER                                   :: left_cols(:), left_rows(:), &
                                                            right_cols(:), right_rows(:)
      LOGICAL                                            :: compatible

      IF (PRESENT(right_overlaps)) THEN
         CALL dbcsr_get_info(overlaps, row_blk_size=left_rows, col_blk_size=left_cols)
         CALL dbcsr_get_info(right_overlaps, row_blk_size=right_rows, col_blk_size=right_cols)
         CPASSERT(SIZE(left_rows) == SIZE(right_rows) .AND. SIZE(left_cols) == SIZE(right_cols))
         compatible = ALL(left_rows == right_rows) .AND. ALL(left_cols == right_cols)
         CPASSERT(compatible)
         CALL dbt_create(right_overlaps, right_b)
         CALL dbt_copy_matrix_to_tensor(right_overlaps, right_b)
      ELSE
         CALL dbt_create(overlaps, right_b)
         CALL dbt_copy_matrix_to_tensor(overlaps, right_b)
      END IF
   END SUBROUTINE right_transition_tensor

! **************************************************************************************************
!> \brief Build all native (RI|Coulomb|Q) center-diagonal blocks on their unique DBT owners.
!>        Full Coulomb tails are evaluated for every RI atom; no GG neighbor list or screening.
!> \param uspp Current centers, original Q and magnetic-projector ordering.
!> \param basis RI basis per kind, with original normalized contraction coefficients.
!> \param particles Atomic positions and kind mapping, in RI-atom order.
!> \param nang Per-kind angular grid control, augmented by the polynomial angular minimum.
!> \param pgrid Three-dimensional tensor process grid.
!> \param auxiliary New C(RI,P,P), with atomwise RI and centerwise projector blocks.
!> \param ri_displacements Optional Cartesian RI-cell translations S, stacked in the RI dimension.
!> \param augmentation_displacement Optional Cartesian projector-cell translation T.
!>        The integral retains the full displacement (R_Q+T)-(R_RI+S).
!> \param cell Optional cell; fold atom coordinates to the same origin as the SAP image labels.
!> \param ri_metric RI operator, Coulomb by default.
!> \param periodic_sum Sum native Q images for the implicit Gamma engine.
! **************************************************************************************************
   SUBROUTINE qs_uspp_hfx_native_tensor(uspp, basis, particles, nang, pgrid, auxiliary, &
                                        ri_displacements, augmentation_displacement, cell, ri_metric, periodic_sum, reconstructed_paw)
      TYPE(qs_uspp_type), INTENT(IN), TARGET             :: uspp
      TYPE(gto_basis_set_p_type), INTENT(IN)             :: basis(:)
      TYPE(particle_type), INTENT(IN)                    :: particles(:)
      INTEGER, INTENT(IN)                                :: nang(:)
      TYPE(dbt_pgrid_type), INTENT(IN)                   :: pgrid
      TYPE(dbt_type), INTENT(OUT)                        :: auxiliary
      REAL(dp), INTENT(IN), OPTIONAL                     :: ri_displacements(:, :), &
                                                            augmentation_displacement(3)
      TYPE(cell_type), OPTIONAL, POINTER                 :: cell
      TYPE(libint_potential_type), INTENT(IN), OPTIONAL  :: ri_metric
      LOGICAL, INTENT(IN), OPTIONAL                      :: periodic_sum
      LOGICAL, INTENT(IN), OPTIONAL :: reconstructed_paw

      INTEGER                                            :: atom, center, first, ikind, iset, &
                                                            kind_q, n, nset, owner, ri_cell, row, &
                                                            rule
      INTEGER, ALLOCATABLE                               :: d1(:), d2(:), d3(:), sizes(:)
      INTEGER, POINTER                                   :: nsgf_set(:)
      REAL(dp)                                           :: q_position(3), qshift(3), ri_position(3)
      REAL(dp), ALLOCATABLE                              :: block(:, :, :), shifts(:, :)
      TYPE(gto_basis_set_type), POINTER                  :: ri_basis
      TYPE(uspp_augmentation_type), POINTER              :: augmentation

      CALL initialize_native(uspp, basis)
      CALL native_image_displacements(ri_displacements, augmentation_displacement, shifts, qshift)
      ALLOCATE (sizes(SIZE(particles)*SIZE(shifts, 2)))
      DO row = 1, SIZE(sizes)
         atom = MODULO(row - 1, SIZE(particles)) + 1
         CPASSERT(ASSOCIATED(particles(atom)%atomic_kind))
         ikind = particles(atom)%atomic_kind%kind_number
         CPASSERT(ikind >= 1 .AND. ikind <= SIZE(basis))
         CPASSERT(ASSOCIATED(basis(ikind)%gto_basis_set))
         CALL get_gto_basis_set(basis(ikind)%gto_basis_set, nsgf=sizes(row))
         CPASSERT(sizes(row) > 0)
      END DO
      CALL create_3c_tensor(auxiliary, d1, d2, d3, pgrid, sizes, uspp%projector_sizes, uspp%projector_sizes, &
                            [1], [2, 3], name="Native RI Q")
      DEALLOCATE (d1, d2, d3)
      DO center = 1, SIZE(uspp%center_atoms)
         kind_q = uspp%center_kinds(center)
         augmentation => native_augmentation(uspp, kind_q, reconstructed_paw)
         n = uspp%projector_sizes(center)
         q_position = uspp%position(:, center)
         IF (PRESENT(cell)) q_position = pbc(q_position, cell)
         q_position = q_position + qshift
         DO row = 1, SIZE(sizes)
            atom = MODULO(row - 1, SIZE(particles)) + 1
            ri_cell = (row - 1)/SIZE(particles) + 1
            ri_position = particles(atom)%r
            IF (PRESENT(cell)) ri_position = pbc(ri_position, cell)
            ri_position = ri_position + shifts(:, ri_cell)
            CALL dbt_get_stored_coordinates(auxiliary, [row, center, center], owner)
            IF (owner /= auxiliary%pgrid%mp_comm_2d%mepos) CYCLE
            ikind = particles(atom)%atomic_kind%kind_number
            ri_basis => basis(ikind)%gto_basis_set
            CALL get_gto_basis_set(ri_basis, nset=nset, nsgf_set=nsgf_set)
            rule = angular_rule(uspp, kind_q, ri_basis, nang(kind_q))
            ALLOCATE (block(sizes(row), n, n), SOURCE=0.0_dp)
            first = 1
            DO iset = 1, nset
               CALL native_auxiliary_block(augmentation, ri_basis, iset, ri_position, q_position, rule, &
                                           block(first:first + nsgf_set(iset) - 1, :, :), &
                                           cell=cell, ri_metric=ri_metric, periodic_sum=periodic_sum)
               first = first + nsgf_set(iset)
            END DO
            CALL dbt_put_block(auxiliary, [row, center, center], SHAPE(block), block)
            DEALLOCATE (block)
         END DO
      END DO
   END SUBROUTINE qs_uspp_hfx_native_tensor

! **************************************************************************************************
!> \brief Contract native C adjoints with RI- and Q-center derivatives on their unique owners.
!>        Native radial functions are fixed in Cartesian coordinates; strain moves their centers.
!> \param uspp Same original Q and centers as the forward tensor.
!> \param basis Same RI basis and normalized coefficients.
!> \param particles Current RI positions and kind mapping.
!> \param nang Same per-kind angular controls as the forward integral.
!> \param auxiliary_gradient Already scaled dE/dC, center-diagonal, in unsplit native layout.
!> \param atom_gradient Local energy gradients (3,natom), incremented without reduction.
!> \param cell_gradient Optional physical affine-strain energy gradient, incremented in place.
!> \param ri_displacements Same Cartesian RI-cell translations S and stacking as the forward C.
!> \param augmentation_displacement Same Cartesian projector-cell translation T.
!> \param cell Same optional cell and atom folding as the forward C.
!> \param ri_metric Same RI operator as the forward tensor.
!> \param periodic_sum Same implicit Gamma image sum as the forward tensor.
! **************************************************************************************************
   SUBROUTINE qs_uspp_hfx_native_gradient(uspp, basis, particles, nang, auxiliary_gradient, &
                                          atom_gradient, cell_gradient, ri_displacements, augmentation_displacement, cell, &
                                          ri_metric, periodic_sum, reconstructed_paw)
      TYPE(qs_uspp_type), INTENT(IN), TARGET             :: uspp
      TYPE(gto_basis_set_p_type), INTENT(IN)             :: basis(:)
      TYPE(particle_type), INTENT(IN)                    :: particles(:)
      INTEGER, INTENT(IN)                                :: nang(:)
      TYPE(dbt_type), INTENT(INOUT)                      :: auxiliary_gradient
      REAL(dp), INTENT(INOUT)                            :: atom_gradient(:, :)
      REAL(dp), INTENT(INOUT), OPTIONAL                  :: cell_gradient(3, 3)
      REAL(dp), INTENT(IN), OPTIONAL                     :: ri_displacements(:, :), &
                                                            augmentation_displacement(3)
      TYPE(cell_type), OPTIONAL, POINTER                 :: cell
      TYPE(libint_potential_type), INTENT(IN), OPTIONAL  :: ri_metric
      LOGICAL, INTENT(IN), OPTIONAL                      :: periodic_sum
      LOGICAL, INTENT(IN), OPTIONAL :: reconstructed_paw

      INTEGER                                            :: atom, axis, axis2, center, first, ikind, &
                                                            ind(3), iset, kind_q, n, nset, qatom, &
                                                            ri_cell, rule
      INTEGER, POINTER                                   :: nsgf_set(:)
      LOGICAL                                            :: found, valid
      REAL(dp)                                           :: gq(3), q_position(3), qshift(3), &
                                                            ri_position(3)
      REAL(dp), ALLOCATABLE                              :: adjoint(:, :, :), block(:, :, :), &
                                                            derivative(:, :, :, :), shifts(:, :), &
                                                            strain(:, :, :, :, :)
      TYPE(dbt_iterator_type)                            :: iter
      TYPE(gto_basis_set_type), POINTER                  :: ri_basis
      TYPE(uspp_augmentation_type), POINTER              :: augmentation

      valid = ALL(SHAPE(atom_gradient) == [3, SIZE(particles)])
      CPASSERT(valid)
      CALL initialize_native(uspp, basis)
      CALL native_image_displacements(ri_displacements, augmentation_displacement, shifts, qshift)
      CALL dbt_iterator_start(iter, auxiliary_gradient)
      DO WHILE (dbt_iterator_blocks_left(iter))
         CALL dbt_iterator_next_block(iter, ind)
         CPASSERT(ind(2) == ind(3))
         CPASSERT(ind(1) >= 1 .AND. ind(1) <= SIZE(particles)*SIZE(shifts, 2))
         atom = MODULO(ind(1) - 1, SIZE(particles)) + 1
         ri_cell = (ind(1) - 1)/SIZE(particles) + 1
         ri_position = particles(atom)%r
         IF (PRESENT(cell)) ri_position = pbc(ri_position, cell)
         ri_position = ri_position + shifts(:, ri_cell)
         center = ind(2)
         q_position = uspp%position(:, center)
         IF (PRESENT(cell)) q_position = pbc(q_position, cell)
         q_position = q_position + qshift
         qatom = uspp%center_atoms(center)
         kind_q = uspp%center_kinds(center)
         augmentation => native_augmentation(uspp, kind_q, reconstructed_paw)
         n = uspp%projector_sizes(center)
         ikind = particles(atom)%atomic_kind%kind_number
         ri_basis => basis(ikind)%gto_basis_set
         CALL get_gto_basis_set(ri_basis, nset=nset, nsgf_set=nsgf_set)
         rule = angular_rule(uspp, kind_q, ri_basis, nang(kind_q))
         CALL dbt_get_block(auxiliary_gradient, ind, adjoint, found)
         CPASSERT(found)
         gq = 0.0_dp
         first = 1
         DO iset = 1, nset
            ALLOCATE (block(nsgf_set(iset), n, n), derivative(3, nsgf_set(iset), n, n))
            IF (PRESENT(cell_gradient)) ALLOCATE (strain(3, 3, nsgf_set(iset), n, n))
            CALL native_auxiliary_block(augmentation, ri_basis, iset, ri_position, q_position, rule, block, &
                                        derivative, strain, cell, ri_metric, periodic_sum)
            DO axis = 1, 3
               gq(axis) = gq(axis) + SUM(adjoint(first:first + nsgf_set(iset) - 1, :, :)*derivative(axis, :, :, :))
               IF (PRESENT(cell_gradient)) THEN
                  DO axis2 = 1, 3
                     cell_gradient(axis, axis2) = cell_gradient(axis, axis2) + &
                                                  SUM(adjoint(first:first + nsgf_set(iset) - 1, :, :)*strain(axis, axis2, :, :, :))
                  END DO
               END IF
            END DO
            first = first + nsgf_set(iset)
            DEALLOCATE (block, derivative)
            IF (ALLOCATED(strain)) DEALLOCATE (strain)
         END DO
         atom_gradient(:, atom) = atom_gradient(:, atom) - gq
         atom_gradient(:, qatom) = atom_gradient(:, qatom) + gq
         DEALLOCATE (adjoint)
      END DO
      CALL dbt_iterator_stop(iter)
   END SUBROUTINE qs_uspp_hfx_native_gradient

! **************************************************************************************************
!> \brief Evaluate native RI coupling and the same center/strain adjoints over Q images.
!>        At Gamma B=sum_T B(T) is invariant under translation of the augmentation center;
!>        hence C=sum_T (RI^0|metric|Q^T) precedes both unweighted B legs. Image support uses
!>        the original Q support and RI Gaussian radius, independently of GG product lists.
!> \param augmentation Original Q radial/angular data.
!> \param basis Normalized RI basis.
!> \param iset RI shell set.
!> \param ri_position Folded RI center.
!> \param q_position Folded augmentation center.
!> \param rule Angular quadrature rule.
!> \param block Integral block.
!> \param derivative Optional Q-center derivative; the RI-center derivative is its negative.
!> \param strain Optional affine-strain derivative including each lattice translation.
!> \param cell Cell for periodic summation.
!> \param ri_metric RI metric; isolated Coulomb by default.
!> \param periodic_sum Sum the implicit Gamma lattice; false by default.
! **************************************************************************************************
   SUBROUTINE native_auxiliary_block(augmentation, basis, iset, ri_position, q_position, rule, &
                                     block, derivative, strain, cell, ri_metric, periodic_sum)
      TYPE(uspp_augmentation_type), INTENT(IN)           :: augmentation
      TYPE(gto_basis_set_type), INTENT(IN)               :: basis
      INTEGER, INTENT(IN)                                :: iset
      REAL(dp), INTENT(IN)                               :: ri_position(3), q_position(3)
      INTEGER, INTENT(IN)                                :: rule
      REAL(dp), INTENT(OUT)                              :: block(:, :, :)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: derivative(:, :, :, :), &
                                                            strain(:, :, :, :, :)
      TYPE(cell_type), OPTIONAL, POINTER                 :: cell
      TYPE(libint_potential_type), INTENT(IN), OPTIONAL  :: ri_metric
      LOGICAL, INTENT(IN), OPTIONAL                      :: periodic_sum

      INTEGER                                            :: axis, axis2, hi(3), ix, iy, iz, lo(3), &
                                                            metric
      LOGICAL                                            :: images
      REAL(dp)                                           :: delta(3), position(3), radius, width(3)
      REAL(dp), ALLOCATABLE                              :: gradient(:, :, :, :), values(:, :, :)

      metric = do_potential_coulomb
      IF (PRESENT(ri_metric)) metric = ri_metric%potential_type
      images = .FALSE.
      IF (PRESENT(periodic_sum)) images = periodic_sum
      lo = 0
      hi = 0
      IF (images) THEN
         CPASSERT(PRESENT(cell))
         CPASSERT(ASSOCIATED(cell))
         CPASSERT(metric == do_potential_id)
         CALL get_gto_basis_set(basis, kind_radius=radius)
         radius = radius + augmentation%realspace_radius
         delta = MATMUL(cell%h_inv, q_position - ri_position)
         width = radius*SQRT(SUM(cell%h_inv**2, DIM=2))
         WHERE (cell%perd /= 0)
            lo = CEILING(-delta - width)
            hi = FLOOR(-delta + width)
         END WHERE
      END IF
      block = 0.0_dp
      IF (PRESENT(derivative)) derivative = 0.0_dp
      IF (PRESENT(strain)) strain = 0.0_dp
      ALLOCATE (values(SIZE(block, 1), SIZE(block, 2), SIZE(block, 3)))
      IF (PRESENT(derivative) .OR. PRESENT(strain)) THEN
         ALLOCATE (gradient(3, SIZE(block, 1), SIZE(block, 2), SIZE(block, 3)))
      END IF
      DO iz = lo(3), hi(3)
         DO iy = lo(2), hi(2)
            DO ix = lo(1), hi(1)
               position = q_position
               IF (images) THEN
                  position = position + MATMUL(cell%hmat, REAL([ix, iy, iz], dp))
                  IF (NORM2(position - ri_position) > radius) CYCLE
               END IF
               CALL uspp_hfx_auxiliary_integrals(augmentation, basis, iset, ri_position, position, &
                                                 lebedev_grid(rule)%r, fourpi*lebedev_grid(rule)%w, values, &
                                                 augmentation_derivative=gradient, potential_type=metric)
               block = block + values
               IF (PRESENT(derivative)) derivative = derivative + gradient
               IF (PRESENT(strain)) THEN
                  DO axis2 = 1, 3
                     DO axis = 1, 3
                        strain(axis, axis2, :, :, :) = strain(axis, axis2, :, :, :) + &
                                                       gradient(axis, :, :, :)*(position(axis2) - ri_position(axis2))
                     END DO
                  END DO
               END IF
            END DO
         END DO
      END DO
   END SUBROUTINE native_auxiliary_block

! **************************************************************************************************
!> \brief Validate and retain Cartesian cell shifts, with an unchanged zero-image default.
!> \param ri_displacements Optional stack of RI-cell translations.
!> \param augmentation_displacement Optional translation of the augmentation center.
!> \param shifts Validated RI-cell translations, with a single zero cell by default.
!> \param qshift Validated augmentation translation, zero by default.
! **************************************************************************************************
   SUBROUTINE native_image_displacements(ri_displacements, augmentation_displacement, shifts, qshift)
      REAL(dp), INTENT(IN), OPTIONAL                     :: ri_displacements(:, :), &
                                                            augmentation_displacement(3)
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: shifts(:, :)
      REAL(dp), INTENT(OUT)                              :: qshift(3)

      IF (PRESENT(ri_displacements)) THEN
         CPASSERT(SIZE(ri_displacements, 1) == 3 .AND. SIZE(ri_displacements, 2) > 0)
         CPASSERT(ALL(ieee_is_finite(ri_displacements)))
         ALLOCATE (shifts, SOURCE=ri_displacements)
      ELSE
         ALLOCATE (shifts(3, 1), SOURCE=0.0_dp)
      END IF
      qshift = 0.0_dp
      IF (PRESENT(augmentation_displacement)) THEN
         CPASSERT(ALL(ieee_is_finite(augmentation_displacement)))
         qshift = augmentation_displacement
      END IF
   END SUBROUTINE native_image_displacements

! **************************************************************************************************
!> \brief Set up original RI basis and check the implemented production operator and orbital scope.
!> \param qs_env ...
!> \param ri_data ...
!> \param ks_env ...
!> \param uspp ...
!> \param kinds ...
!> \param particles ...
!> \param basis ...
!> \param nang ...
!> \param with_images Use explicit independent projector images for the k-point engine.
! **************************************************************************************************
   SUBROUTINE setup(qs_env, ri_data, ks_env, uspp, kinds, particles, basis, nang, with_images)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(hfx_ri_type), INTENT(IN)                      :: ri_data
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_uspp_type), POINTER                        :: uspp
      TYPE(qs_kind_type), POINTER                        :: kinds(:)
      TYPE(particle_type), POINTER                       :: particles(:)
      TYPE(gto_basis_set_p_type), ALLOCATABLE            :: basis(:)
      INTEGER, ALLOCATABLE                               :: nang(:)
      LOGICAL, INTENT(IN), OPTIONAL                      :: with_images

      INTEGER                                            :: ikind
      LOGICAL                                            :: do_kpoints, image_resolved, supported
      TYPE(cell_type), POINTER                           :: cell
      TYPE(neighbor_list_set_p_type), POINTER            :: sap_ppnl(:)

      CALL get_qs_env(qs_env, ks_env=ks_env)
      CALL projector_state(qs_env, ks_env, ri_data%orb_basis_type, uspp, sap_ppnl)
      IF (.NOT. ASSOCIATED(uspp)) RETURN
      CALL get_qs_env(qs_env, qs_kind_set=kinds, particle_set=particles, cell=cell)
      CALL get_ks_env(ks_env, do_kpoints=do_kpoints)
      image_resolved = .FALSE.
      IF (PRESENT(with_images)) image_resolved = with_images
      supported = ri_data%flavor == ri_pmat .AND. &
                  (ri_data%ri_metric%potential_type == do_potential_coulomb .OR. &
                   ri_data%ri_metric%potential_type == do_potential_id)
      IF (image_resolved) THEN
         supported = supported .AND. do_kpoints .AND. &
                     ri_data%ri_metric%potential_type == do_potential_coulomb .AND. &
                     (ri_data%hfx_pot%potential_type == do_potential_truncated .OR. &
                      ri_data%hfx_pot%potential_type == do_potential_short)
      ELSE
         supported = supported .AND. .NOT. do_kpoints
         IF (ANY(cell%perd /= 0)) THEN
            supported = supported .AND. ri_data%ri_metric%potential_type == do_potential_id .AND. &
                        (ri_data%hfx_pot%potential_type == do_potential_truncated .OR. &
                         ri_data%hfx_pot%potential_type == do_potential_short)
         ELSE
            supported = supported .AND. &
                        (ri_data%hfx_pot%potential_type == do_potential_coulomb .OR. &
                         ri_data%hfx_pot%potential_type == do_potential_short)
         END IF
      END IF
      IF (.NOT. supported) THEN
         CPABORT('Unsupported native Coulomb RI-HFX image setup')
      END IF
      CPASSERT(ASSOCIATED(uspp%overlaps))
      CPASSERT(ASSOCIATED(ri_data%pgrid))
      ALLOCATE (basis(SIZE(kinds)), nang(SIZE(kinds)))
      CALL qs_uspp_prepare_paw_reconstruction(uspp, kinds)
      CALL basis_set_list_setup(basis, ri_data%ri_basis_type, kinds)
      DO ikind = 1, SIZE(kinds)
         nang(ikind) = kinds(ikind)%ngrid_ang
      END DO
   END SUBROUTINE setup

! **************************************************************************************************
!> \brief Select the owned projector state of the actual RI orbital basis.
!> \param qs_env ...
!> \param ks_env ...
!> \param basis_type ...
!> \param uspp ...
!> \param neighbors ...
! **************************************************************************************************
   SUBROUTINE projector_state(qs_env, ks_env, basis_type, uspp, neighbors)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      CHARACTER(LEN=*), INTENT(IN)                       :: basis_type
      TYPE(qs_uspp_type), POINTER                        :: uspp
      TYPE(neighbor_list_set_p_type), POINTER            :: neighbors(:)

      TYPE(admm_type), POINTER                           :: admm_env

      CALL get_ks_env(ks_env, uspp=uspp, sap_ppnl=neighbors)
      IF (.NOT. ASSOCIATED(uspp)) RETURN
      IF (basis_type == 'ORB') RETURN
      CALL get_qs_env(qs_env, admm_env=admm_env)
      CPASSERT(ASSOCIATED(admm_env))
      CPASSERT(basis_type == admm_env%uspp_basis_type)
      CPASSERT(ASSOCIATED(admm_env%uspp_aux))
      uspp => admm_env%uspp_aux
      neighbors => admm_env%sap_aux_ppnl
   END SUBROUTINE projector_state

! **************************************************************************************************
!> \brief Prepare matrix-derived B and unpartitioned intermediates for exact DBT contractions.
!> \param overlaps ...
!> \param auxiliary ...
!> \param b ...
!> \param u ...
!> \param full ...
! **************************************************************************************************
   SUBROUTINE transition_templates(overlaps, auxiliary, b, u, full)
      TYPE(dbcsr_type), INTENT(IN)                       :: overlaps
      TYPE(dbt_type), INTENT(IN)                         :: auxiliary
      TYPE(dbt_type), INTENT(OUT)                        :: b, u, full

      INTEGER                                            :: nblks(3)
      INTEGER, ALLOCATABLE                               :: d1(:), d2(:), d3(:), sizes_p(:), &
                                                            sizes_p_right(:), sizes_ri(:)
      INTEGER, POINTER                                   :: sizes_ao(:), sizes_b_p(:)
      LOGICAL                                            :: compatible

      CALL dbcsr_get_info(overlaps, row_blk_size=sizes_ao, col_blk_size=sizes_b_p)
      CALL dbt_get_info(auxiliary, nblks_total=nblks)
      ALLOCATE (sizes_ri(nblks(1)), sizes_p(nblks(2)), sizes_p_right(nblks(3)))
      CALL dbt_get_info(auxiliary, blk_size_1=sizes_ri, blk_size_2=sizes_p, blk_size_3=sizes_p_right)
      compatible = nblks(2) == nblks(3) .AND. SIZE(sizes_p) == SIZE(sizes_b_p)
      CPASSERT(compatible)
      compatible = ALL(sizes_p == sizes_b_p)
      CPASSERT(compatible)
      compatible = ALL(sizes_p_right == sizes_b_p)
      CPASSERT(compatible)
      CALL dbt_create(overlaps, b)
      CALL create_3c_tensor(u, d1, d2, d3, auxiliary%pgrid, sizes_ri, sizes_p, sizes_ao, &
                            [1, 2], [3], name="Native RI P AO")
      DEALLOCATE (d1, d2, d3)
      CALL create_3c_tensor(full, d1, d2, d3, auxiliary%pgrid, sizes_ri, sizes_ao, sizes_ao, &
                            [1], [2, 3], name="Native RI AO AO")
      DEALLOCATE (d1, d2, d3)
   END SUBROUTINE transition_templates

! **************************************************************************************************
!> \brief Initialize shared Gaussian tables before serial native-block evaluation.
!> \param uspp ...
!> \param basis ...
! **************************************************************************************************
   SUBROUTINE initialize_native(uspp, basis)
      TYPE(qs_uspp_type), INTENT(IN)                     :: uspp
      TYPE(gto_basis_set_p_type), INTENT(IN)             :: basis(:)

      INTEGER                                            :: ikind, maxl
      INTEGER, POINTER                                   :: lmax(:)

      maxl = MAXVAL(uspp%augmentation%lmax)
      DO ikind = 1, SIZE(basis)
         IF (.NOT. ASSOCIATED(basis(ikind)%gto_basis_set)) CYCLE
         CALL get_gto_basis_set(basis(ikind)%gto_basis_set, lmax=lmax)
         maxl = MAX(maxl, MAXVAL(lmax))
      END DO
      CALL uspp_hfx_integrals_init(maxl)
   END SUBROUTINE initialize_native

! **************************************************************************************************
!> \brief Respect the per-kind angular control and the coincident-center polynomial minimum.
!> \param uspp ...
!> \param kind_q ...
!> \param basis ...
!> \param nang ...
!> \return ...
! **************************************************************************************************
   FUNCTION angular_rule(uspp, kind_q, basis, nang) RESULT(rule)
      TYPE(qs_uspp_type), INTENT(IN)                     :: uspp
      INTEGER, INTENT(IN)                                :: kind_q
      TYPE(gto_basis_set_type), INTENT(IN)               :: basis
      INTEGER, INTENT(IN)                                :: nang
      INTEGER                                            :: rule

      INTEGER, POINTER                                   :: lmax(:)

      CALL get_gto_basis_set(basis, lmax=lmax)
      rule = MAX(get_number_of_lebedev_grid(n=nang), &
                 get_number_of_lebedev_grid(l=uspp%augmentation(kind_q)%lmax + MAXVAL(lmax) + 1))
   END FUNCTION angular_rule

! **************************************************************************************************
!> \brief Select the same immutable pair density for a native RI forward map and its adjoint.
!> \param uspp Original compensated charge data and optional PAW valence reconstruction cache.
!> \param kind_q Native augmentation kind.
!> \param reconstructed_paw Optional explicit request; USPP kinds retain their original Q.
!> \return Pair density in the unchanged native radial/magnetic ordering.
! **************************************************************************************************
   FUNCTION native_augmentation(uspp, kind_q, reconstructed_paw) RESULT(augmentation)
      TYPE(qs_uspp_type), INTENT(IN), TARGET             :: uspp
      INTEGER, INTENT(IN)                                :: kind_q
      LOGICAL, INTENT(IN), OPTIONAL                      :: reconstructed_paw
      TYPE(uspp_augmentation_type), POINTER              :: augmentation

      augmentation => uspp%augmentation(kind_q)
      IF (.NOT. PRESENT(reconstructed_paw)) RETURN
      IF (.NOT. reconstructed_paw) RETURN
      CPASSERT(ALLOCATED(uspp%paw))
      IF (.NOT. uspp%paw(kind_q)%active) RETURN
      CPASSERT(ALLOCATED(uspp%paw_reconstruction))
      CPASSERT(uspp%paw_reconstruction(kind_q)%nproj == augmentation%nproj)
      CPASSERT(uspp%paw_reconstruction(kind_q)%lmax == augmentation%lmax)
      CPASSERT(ALL(uspp%paw_reconstruction(kind_q)%radial_index == augmentation%radial_index))
      CPASSERT(ALL(uspp%paw_reconstruction(kind_q)%angular_index == augmentation%angular_index))
      augmentation => uspp%paw_reconstruction(kind_q)
   END FUNCTION native_augmentation

END MODULE qs_uspp_hfx
