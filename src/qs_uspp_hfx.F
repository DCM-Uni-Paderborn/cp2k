!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Distributed native USPP transition densities and their RI-HFX geometry adjoint.
!>        The initial operator is isolated scalar Coulomb, Gamma, RI-RHO, and the ORB basis.
!>        Complete AO pairs are retained: the global RI contraction includes intercenter QQ.
! **************************************************************************************************
MODULE qs_uspp_hfx
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_p_type,&
                                              gto_basis_set_type
   USE cell_types,                      ONLY: cell_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_create,&
                                              dbcsr_get_info,&
                                              dbcsr_release,&
                                              dbcsr_type
   USE dbt_api,                         ONLY: &
        dbt_contract, dbt_copy, dbt_copy_matrix_to_tensor, dbt_copy_tensor_to_matrix, dbt_create, &
        dbt_destroy, dbt_get_block, dbt_get_info, dbt_get_stored_coordinates, &
        dbt_iterator_blocks_left, dbt_iterator_next_block, dbt_iterator_start, dbt_iterator_stop, &
        dbt_iterator_type, dbt_pgrid_type, dbt_put_block, dbt_scale, dbt_type
   USE hfx_types,                       ONLY: hfx_ri_type
   USE input_constants,                 ONLY: do_potential_coulomb
   USE input_cp2k_hfx,                  ONLY: ri_pmat
   USE kinds,                           ONLY: dp
   USE lebedev,                         ONLY: get_number_of_lebedev_grid,&
                                              lebedev_grid
   USE mathconstants,                   ONLY: fourpi
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_integral_utils,               ONLY: basis_set_list_setup
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_tensors_types,                ONLY: create_3c_tensor
   USE qs_uspp_projectors,              ONLY: qs_uspp_projector_gradient
   USE qs_uspp_types,                   ONLY: qs_uspp_type
   USE sap_kind_types,                  ONLY: build_sap_ints,&
                                              release_sap_int,&
                                              sap_int_type,&
                                              sap_sort
   USE uspp_hfx_integrals,              ONLY: uspp_hfx_auxiliary_integrals,&
                                              uspp_hfx_integrals_init
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_uspp_hfx'
   PUBLIC :: qs_uspp_hfx_add_tensor, qs_uspp_hfx_gradient
   PUBLIC :: qs_uspp_hfx_transition_tensor, qs_uspp_hfx_transition_adjoint
   PUBLIC :: qs_uspp_hfx_native_tensor, qs_uspp_hfx_native_gradient

CONTAINS

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
      TYPE(dbt_type)                                     :: auxiliary
      TYPE(gto_basis_set_p_type), ALLOCATABLE            :: basis(:)
      TYPE(particle_type), POINTER                       :: particles(:)
      TYPE(qs_kind_type), POINTER                        :: kinds(:)
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_uspp_type), POINTER                        :: uspp

      CALL setup(qs_env, ri_data, ks_env, uspp, kinds, particles, basis, nang)
      IF (.NOT. ASSOCIATED(uspp)) RETURN
      CALL qs_uspp_hfx_native_tensor(uspp, basis, particles, nang, ri_data%pgrid, auxiliary)
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
      CALL qs_uspp_hfx_native_tensor(uspp, basis, particles, nang, ri_data%pgrid, auxiliary)
      CALL qs_uspp_hfx_transition_adjoint(uspp%overlaps, auxiliary, tensor_adjoint, pref, &
                                          overlap_gradient, auxiliary_gradient)
      CALL get_ks_env(ks_env, sap_ppnl=sap_ppnl)
      CPASSERT(ASSOCIATED(sap_ppnl))
      ALLOCATE (sap(SIZE(kinds)**2))
      CALL build_sap_ints(sap, sap_ppnl, kinds, nder=1)
      CALL sap_sort(sap)
      CALL qs_uspp_projector_gradient(overlap_gradient, sap, uspp%center_atoms, atom_gradient, cell_gradient)
      CALL qs_uspp_hfx_native_gradient(uspp, basis, particles, nang, auxiliary_gradient, &
                                       atom_gradient, cell_gradient)
      CALL release_sap_int(sap)
      CALL dbcsr_release(overlap_gradient)
      CALL dbt_destroy(auxiliary_gradient)
      CALL dbt_destroy(auxiliary)
   END SUBROUTINE qs_uspp_hfx_gradient

! **************************************************************************************************
!> \brief Distributed T += B C B^T with C(RI,P,P), B(AO,P); no packing or symmetry factor.
!>        C has only center-diagonal P blocks. T keeps all AO pairs and all RI centers.
!> \param overlaps Original rectangular distributed B matrix; never globally replicated.
!> \param auxiliary Complete native C, with unsplit RI and magnetic-projector blocks.
!> \param tensor Target full (RI,AO,AO), incremented after reblocking and redistribution.
! **************************************************************************************************
   SUBROUTINE qs_uspp_hfx_transition_tensor(overlaps, auxiliary, tensor)
      TYPE(dbcsr_type), INTENT(IN), TARGET               :: overlaps
      TYPE(dbt_type), INTENT(INOUT)                      :: auxiliary, tensor

      TYPE(dbt_type)                                     :: b, full, u

      CALL transition_templates(overlaps, auxiliary, b, u, full)
      CALL dbt_copy_matrix_to_tensor(overlaps, b)
      CALL dbt_contract(1.0_dp, auxiliary, b, 0.0_dp, u, &
                        contract_1=[3], notcontract_1=[1, 2], &
                        contract_2=[2], notcontract_2=[1], map_1=[1, 2], map_2=[3], filter_eps=0.0_dp)
      CALL dbt_contract(1.0_dp, b, u, 0.0_dp, full, &
                        contract_1=[2], notcontract_1=[1], &
                        contract_2=[2], notcontract_2=[1, 3], map_1=[2], map_2=[1, 3], filter_eps=0.0_dp)
      CALL dbt_copy(full, tensor, summation=.TRUE., move_data=.TRUE.)
      CALL dbt_destroy(full)
      CALL dbt_destroy(u)
      CALL dbt_destroy(b)
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
! **************************************************************************************************
   SUBROUTINE qs_uspp_hfx_transition_adjoint(overlaps, auxiliary, tensor_adjoint, pref, &
                                             overlap_gradient, auxiliary_gradient)
      TYPE(dbcsr_type), INTENT(IN), TARGET               :: overlaps
      TYPE(dbt_type), INTENT(INOUT)                      :: auxiliary, tensor_adjoint
      REAL(dp), INTENT(IN)                               :: pref
      TYPE(dbcsr_type), INTENT(OUT)                      :: overlap_gradient
      TYPE(dbt_type), INTENT(OUT)                        :: auxiliary_gradient

      TYPE(dbt_type)                                     :: b, bar_b, bar_u, full, u

      CALL transition_templates(overlaps, auxiliary, b, u, full)
      CALL dbt_copy_matrix_to_tensor(overlaps, b)
      CALL dbt_contract(1.0_dp, auxiliary, b, 0.0_dp, u, &
                        contract_1=[3], notcontract_1=[1, 2], &
                        contract_2=[2], notcontract_2=[1], map_1=[1, 2], map_2=[3], filter_eps=0.0_dp)
      CALL dbt_copy(tensor_adjoint, full)
      CALL dbt_scale(full, pref)
      CALL dbt_create(u, bar_u)
      CALL dbt_create(b, bar_b)
      CALL dbt_contract(1.0_dp, full, u, 0.0_dp, bar_b, &
                        contract_1=[1, 3], notcontract_1=[2], &
                        contract_2=[1, 3], notcontract_2=[2], map_1=[1], map_2=[2], filter_eps=0.0_dp)
      CALL dbt_contract(1.0_dp, full, b, 0.0_dp, bar_u, &
                        contract_1=[2], notcontract_1=[1, 3], &
                        contract_2=[1], notcontract_2=[2], map_1=[1, 3], map_2=[2], filter_eps=0.0_dp)
      CALL dbt_contract(1.0_dp, bar_u, auxiliary, 1.0_dp, bar_b, &
                        contract_1=[1, 2], notcontract_1=[3], &
                        contract_2=[1, 2], notcontract_2=[3], map_1=[1], map_2=[2], filter_eps=0.0_dp)
      CALL dbt_create(auxiliary, auxiliary_gradient)
      CALL dbt_copy(auxiliary, auxiliary_gradient)
      CALL dbt_contract(1.0_dp, bar_u, b, 0.0_dp, auxiliary_gradient, &
                        contract_1=[3], notcontract_1=[1, 2], &
                        contract_2=[1], notcontract_2=[2], map_1=[1, 2], map_2=[3], filter_eps=0.0_dp, &
                        retain_sparsity=.TRUE.)
      ! bar_b retains the same matrix-derived block layout and ownership as the original B.
      CALL dbcsr_create(overlap_gradient, template=overlaps)
      CALL dbt_copy_tensor_to_matrix(bar_b, overlap_gradient)
      CALL dbt_destroy(bar_b)
      CALL dbt_destroy(bar_u)
      CALL dbt_destroy(full)
      CALL dbt_destroy(u)
      CALL dbt_destroy(b)
   END SUBROUTINE qs_uspp_hfx_transition_adjoint

! **************************************************************************************************
!> \brief Build all native (RI|Coulomb|Q) center-diagonal blocks on their unique DBT owners.
!>        Full Coulomb tails are evaluated for every RI atom; no GG neighbor list or screening.
!> \param uspp Current centers, original Q and magnetic-projector ordering.
!> \param basis RI basis per kind, with original normalized contraction coefficients.
!> \param particles Atomic positions and kind mapping, in RI-atom order.
!> \param nang Per-kind angular grid control, augmented by the polynomial angular minimum.
!> \param pgrid Three-dimensional tensor process grid.
!> \param auxiliary New C(RI,P,P), with atomwise RI and centerwise projector blocks.
! **************************************************************************************************
   SUBROUTINE qs_uspp_hfx_native_tensor(uspp, basis, particles, nang, pgrid, auxiliary)
      TYPE(qs_uspp_type), INTENT(IN)                     :: uspp
      TYPE(gto_basis_set_p_type), INTENT(IN)             :: basis(:)
      TYPE(particle_type), INTENT(IN)                    :: particles(:)
      INTEGER, INTENT(IN)                                :: nang(:)
      TYPE(dbt_pgrid_type), INTENT(IN)                   :: pgrid
      TYPE(dbt_type), INTENT(OUT)                        :: auxiliary

      INTEGER                                            :: atom, center, first, ikind, iset, &
                                                            kind_q, n, nset, owner, rule
      INTEGER, ALLOCATABLE                               :: d1(:), d2(:), d3(:), sizes(:)
      INTEGER, POINTER                                   :: nsgf_set(:)
      REAL(dp), ALLOCATABLE                              :: block(:, :, :)
      TYPE(gto_basis_set_type), POINTER                  :: ri_basis

      CALL initialize_native(uspp, basis)
      ALLOCATE (sizes(SIZE(particles)))
      DO atom = 1, SIZE(particles)
         CPASSERT(ASSOCIATED(particles(atom)%atomic_kind))
         ikind = particles(atom)%atomic_kind%kind_number
         CPASSERT(ikind >= 1 .AND. ikind <= SIZE(basis))
         CPASSERT(ASSOCIATED(basis(ikind)%gto_basis_set))
         CALL get_gto_basis_set(basis(ikind)%gto_basis_set, nsgf=sizes(atom))
         CPASSERT(sizes(atom) > 0)
      END DO
      CALL create_3c_tensor(auxiliary, d1, d2, d3, pgrid, sizes, uspp%projector_sizes, uspp%projector_sizes, &
                            [1], [2, 3], name="Native RI Q")
      DEALLOCATE (d1, d2, d3)
      DO center = 1, SIZE(uspp%center_atoms)
         kind_q = uspp%center_kinds(center)
         n = uspp%projector_sizes(center)
         DO atom = 1, SIZE(particles)
            CALL dbt_get_stored_coordinates(auxiliary, [atom, center, center], owner)
            IF (owner /= auxiliary%pgrid%mp_comm_2d%mepos) CYCLE
            ikind = particles(atom)%atomic_kind%kind_number
            ri_basis => basis(ikind)%gto_basis_set
            CALL get_gto_basis_set(ri_basis, nset=nset, nsgf_set=nsgf_set)
            rule = angular_rule(uspp, kind_q, ri_basis, nang(kind_q))
            ALLOCATE (block(sizes(atom), n, n), SOURCE=0.0_dp)
            first = 1
            DO iset = 1, nset
               CALL uspp_hfx_auxiliary_integrals(uspp%augmentation(kind_q), ri_basis, iset, particles(atom)%r, &
                                                 uspp%position(:, center), lebedev_grid(rule)%r, &
                                                 fourpi*lebedev_grid(rule)%w, block(first:first + nsgf_set(iset) - 1, :, :))
               first = first + nsgf_set(iset)
            END DO
            CALL dbt_put_block(auxiliary, [atom, center, center], SHAPE(block), block)
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
! **************************************************************************************************
   SUBROUTINE qs_uspp_hfx_native_gradient(uspp, basis, particles, nang, auxiliary_gradient, &
                                          atom_gradient, cell_gradient)
      TYPE(qs_uspp_type), INTENT(IN)                     :: uspp
      TYPE(gto_basis_set_p_type), INTENT(IN)             :: basis(:)
      TYPE(particle_type), INTENT(IN)                    :: particles(:)
      INTEGER, INTENT(IN)                                :: nang(:)
      TYPE(dbt_type), INTENT(INOUT)                      :: auxiliary_gradient
      REAL(dp), INTENT(INOUT)                            :: atom_gradient(:, :)
      REAL(dp), INTENT(INOUT), OPTIONAL                  :: cell_gradient(3, 3)

      INTEGER                                            :: atom, axis, center, first, ikind, &
                                                            ind(3), iset, kind_q, n, nset, qatom, &
                                                            rule
      INTEGER, POINTER                                   :: nsgf_set(:)
      LOGICAL                                            :: found, valid
      REAL(dp)                                           :: displacement(3), gq(3)
      REAL(dp), ALLOCATABLE                              :: adjoint(:, :, :), block(:, :, :), &
                                                            derivative(:, :, :, :)
      TYPE(dbt_iterator_type)                            :: iter
      TYPE(gto_basis_set_type), POINTER                  :: ri_basis

      valid = ALL(SHAPE(atom_gradient) == [3, SIZE(particles)])
      CPASSERT(valid)
      CALL initialize_native(uspp, basis)
      CALL dbt_iterator_start(iter, auxiliary_gradient)
      DO WHILE (dbt_iterator_blocks_left(iter))
         CALL dbt_iterator_next_block(iter, ind)
         CPASSERT(ind(2) == ind(3))
         atom = ind(1)
         center = ind(2)
         qatom = uspp%center_atoms(center)
         kind_q = uspp%center_kinds(center)
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
            CALL uspp_hfx_auxiliary_integrals(uspp%augmentation(kind_q), ri_basis, iset, particles(atom)%r, &
                                              uspp%position(:, center), lebedev_grid(rule)%r, &
                                              fourpi*lebedev_grid(rule)%w, block, augmentation_derivative=derivative)
            DO axis = 1, 3
               gq(axis) = gq(axis) + SUM(adjoint(first:first + nsgf_set(iset) - 1, :, :)*derivative(axis, :, :, :))
            END DO
            first = first + nsgf_set(iset)
            DEALLOCATE (block, derivative)
         END DO
         atom_gradient(:, atom) = atom_gradient(:, atom) - gq
         atom_gradient(:, qatom) = atom_gradient(:, qatom) + gq
         IF (PRESENT(cell_gradient)) THEN
            displacement = uspp%position(:, center) - particles(atom)%r
            DO axis = 1, 3
               cell_gradient(axis, :) = cell_gradient(axis, :) + gq(axis)*displacement
            END DO
         END IF
         DEALLOCATE (adjoint)
      END DO
      CALL dbt_iterator_stop(iter)
   END SUBROUTINE qs_uspp_hfx_native_gradient

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
! **************************************************************************************************
   SUBROUTINE setup(qs_env, ri_data, ks_env, uspp, kinds, particles, basis, nang)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(hfx_ri_type), INTENT(IN)                      :: ri_data
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_uspp_type), POINTER                        :: uspp
      TYPE(qs_kind_type), POINTER                        :: kinds(:)
      TYPE(particle_type), POINTER                       :: particles(:)
      TYPE(gto_basis_set_p_type), ALLOCATABLE            :: basis(:)
      INTEGER, ALLOCATABLE                               :: nang(:)

      INTEGER                                            :: ikind
      LOGICAL                                            :: do_kpoints, supported
      TYPE(cell_type), POINTER                           :: cell

      CALL get_qs_env(qs_env, ks_env=ks_env)
      CALL get_ks_env(ks_env, uspp=uspp)
      IF (.NOT. ASSOCIATED(uspp)) RETURN
      CALL get_qs_env(qs_env, qs_kind_set=kinds, particle_set=particles, cell=cell)
      CALL get_ks_env(ks_env, do_kpoints=do_kpoints)
      supported = ri_data%flavor == ri_pmat .AND. .NOT. do_kpoints .AND. ALL(cell%perd == 0) .AND. &
                  TRIM(ri_data%orb_basis_type) == 'ORB' .AND. ri_data%same_op .AND. &
                  ri_data%ri_metric%potential_type == do_potential_coulomb .AND. &
                  ri_data%hfx_pot%potential_type == do_potential_coulomb
      IF (.NOT. supported) THEN
         CPABORT('Native USPP HFX requires isolated Gamma Coulomb RI-RHO and ORB')
      END IF
      CPASSERT(ASSOCIATED(uspp%overlaps))
      CPASSERT(ASSOCIATED(ri_data%pgrid))
      ALLOCATE (basis(SIZE(kinds)), nang(SIZE(kinds)))
      CALL basis_set_list_setup(basis, ri_data%ri_basis_type, kinds)
      DO ikind = 1, SIZE(kinds)
         nang(ikind) = kinds(ikind)%ngrid_ang
      END DO
   END SUBROUTINE setup

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

END MODULE qs_uspp_hfx
