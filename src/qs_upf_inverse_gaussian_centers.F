!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
!> \brief Complete native Gaussian center expansions for two independent GAPW partitions.
MODULE qs_upf_inverse_gaussian_centers
   USE ao_util,                         ONLY: exp_radius
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
   USE cell_types,                      ONLY: cell_type
   USE cp_control_types,                ONLY: dft_control_type
   USE distribution_1d_types,           ONLY: distribution_1d_type
   USE distribution_2d_types,           ONLY: distribution_2d_type
   USE input_section_types,             ONLY: section_vals_val_get
   USE kinds,                           ONLY: dp
   USE molecule_types,                  ONLY: molecule_type
   USE particle_types,                  ONLY: particle_type
   USE paw_proj_set_types,              ONLY: allocate_paw_proj_set,&
                                              deallocate_paw_proj_set,&
                                              projectors
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type,&
                                              release_neighbor_list_sets
   USE qs_neighbor_lists,               ONLY: atom2d_build,&
                                              atom2d_cleanup,&
                                              build_neighbor_lists,&
                                              local_atoms_type,&
                                              pair_radius_setup
   USE qs_oce_methods,                  ONLY: build_oce_matrices
   USE qs_oce_types,                    ONLY: allocate_oce_set,&
                                              create_oce_set,&
                                              deallocate_oce_set,&
                                              oce_matrix_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   ! Kind copies borrow native bases and grids. Only the added projectors and OCE are owned.
   TYPE qs_upf_inverse_gaussian_centers_type
      TYPE(qs_kind_type), POINTER :: kinds(:) => NULL()
      TYPE(oce_matrix_type), POINTER :: oce => NULL()
      LOGICAL, ALLOCATABLE :: owns_projector(:)
   END TYPE qs_upf_inverse_gaussian_centers_type

   PUBLIC :: qs_upf_inverse_gaussian_centers_type
   PUBLIC :: qs_upf_inverse_gaussian_centers_prepare, qs_upf_inverse_gaussian_centers_release

CONTAINS

! **************************************************************************************************
!> \brief Expand both unchanged Gaussian bases at every center active in either GAPW partition.
!>        An originally unsplit basis has identical hard and soft expansions at an added center.
!>        No native basis, soft partition, projector or OCE owner is modified.
!> \param qs_env Native geometry, distribution and GAPW controls shared by both basis roles.
!> \param row_kinds Native row-basis kinds, including their original ORB and GAPW_1C bases.
!> \param column_kinds Independent column-basis kinds with the same atom-kind assignment.
!> \param row New owned row OCE and shallow kind copies, valid at this geometry.
!> \param column New independent column OCE and shallow kind copies.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_gaussian_centers_prepare(qs_env, row_kinds, column_kinds, row, column)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_kind_type), INTENT(IN)                     :: row_kinds(:), column_kinds(:)
      TYPE(qs_upf_inverse_gaussian_centers_type), &
         INTENT(OUT)                                     :: row, column

      INTEGER                                            :: kind
      LOGICAL, ALLOCATABLE                               :: centers(:)

      CPASSERT(SIZE(row_kinds) == SIZE(column_kinds))
      ALLOCATE (centers(SIZE(row_kinds)))
      DO kind = 1, SIZE(centers)
         centers(kind) = row_kinds(kind)%paw_atom .OR. column_kinds(kind)%paw_atom
      END DO
      CALL prepare_centers(qs_env, row_kinds, centers, row)
      CALL prepare_centers(qs_env, column_kinds, centers, column)
   END SUBROUTINE qs_upf_inverse_gaussian_centers_prepare

! **************************************************************************************************
!> \brief Release private OCE and added projectors, leaving every borrowed native field intact.
!> \param owner Prepared center expansion at an obsolete geometry.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_gaussian_centers_release(owner)
      TYPE(qs_upf_inverse_gaussian_centers_type), &
         INTENT(INOUT)                                   :: owner

      INTEGER                                            :: kind

      CALL deallocate_oce_set(owner%oce)
      IF (ASSOCIATED(owner%kinds)) THEN
         DO kind = 1, SIZE(owner%kinds)
            IF (owner%owns_projector(kind)) CALL deallocate_paw_proj_set(owner%kinds(kind)%paw_proj_set)
         END DO
         ! A native kind destructor would release the borrowed basis and grid pointers.
         DEALLOCATE (owner%kinds)
      END IF
      IF (ALLOCATED(owner%owns_projector)) DEALLOCATE (owner%owns_projector)
      owner = qs_upf_inverse_gaussian_centers_type()
   END SUBROUTINE qs_upf_inverse_gaussian_centers_release

! **************************************************************************************************
!> \brief Build native overlaps and their first geometry derivatives on the completed center set.
!> \param qs_env ...
!> \param kinds ...
!> \param centers Union of the two original Gaussian center sets.
!> \param owner ...
! **************************************************************************************************
   SUBROUTINE prepare_centers(qs_env, kinds, centers, owner)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_kind_type), INTENT(IN)                     :: kinds(:)
      LOGICAL, INTENT(IN)                                :: centers(:)
      TYPE(qs_upf_inverse_gaussian_centers_type), &
         INTENT(OUT)                                     :: owner

      INTEGER                                            :: kind, l, nkind, p
      LOGICAL, ALLOCATABLE                               :: present_basis(:)
      REAL(dp)                                           :: subcells
      REAL(dp), ALLOCATABLE                              :: basis_radius(:), center_radius(:), &
                                                            pair_radius(:, :)
      TYPE(atomic_kind_type), POINTER                    :: atomic_kinds(:)
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: control
      TYPE(distribution_1d_type), POINTER                :: distribution_1d
      TYPE(distribution_2d_type), POINTER                :: distribution_2d
      TYPE(gto_basis_set_type), POINTER                  :: basis, center_basis
      TYPE(local_atoms_type), ALLOCATABLE                :: atom2d(:)
      TYPE(molecule_type), POINTER                       :: molecules(:)
      TYPE(neighbor_list_set_p_type), POINTER            :: sap(:)
      TYPE(particle_type), POINTER                       :: particles(:)

      CALL get_qs_env(qs_env, dft_control=control, atomic_kind_set=atomic_kinds, cell=cell, &
                      distribution_2d=distribution_2d, local_particles=distribution_1d, &
                      particle_set=particles, molecule_set=molecules)
      nkind = SIZE(kinds)
      CPASSERT(SIZE(centers) == nkind)
      ALLOCATE (owner%kinds(nkind), owner%owns_projector(nkind))
      owner%kinds(:) = kinds
      owner%owns_projector = .FALSE.
      ALLOCATE (present_basis(nkind), basis_radius(nkind), center_radius(nkind), pair_radius(nkind, nkind))
      present_basis = .FALSE.
      basis_radius = 0.0_dp
      center_radius = 0.0_dp
      DO kind = 1, nkind
         CALL get_qs_kind(kinds(kind), basis_set=basis)
         IF (ASSOCIATED(basis)) THEN
            present_basis(kind) = .TRUE.
            CALL get_gto_basis_set(basis, kind_radius=basis_radius(kind))
         END IF
         IF (.NOT. centers(kind)) CYCLE
         IF (.NOT. kinds(kind)%paw_atom) THEN
            CALL get_qs_kind(kinds(kind), basis_set=center_basis, basis_type="GAPW_1C")
            CPASSERT(ASSOCIATED(basis) .AND. ASSOCIATED(center_basis))
            NULLIFY (owner%kinds(kind)%paw_proj_set)
            CALL allocate_paw_proj_set(owner%kinds(kind)%paw_proj_set)
            owner%owns_projector(kind) = .TRUE.
            CALL projectors(owner%kinds(kind)%paw_proj_set, center_basis, basis, owner%kinds(kind)%hard_radius, &
                            control%qs_control, owner%kinds(kind)%max_rad_local, qs_env%input, unsplit=.TRUE.)
            owner%kinds(kind)%paw_atom = .TRUE.
            ASSOCIATE (projector => owner%kinds(kind)%paw_proj_set)
               projector%rcprj = 0.0_dp
               projector%rzetprj = 0.0_dp
               DO l = 0, projector%maxl
                  DO p = 1, projector%nprj(l)
                     projector%rzetprj(p, l) = exp_radius(l, projector%zetprj(p, l), &
                                                          control%qs_control%eps_pgf_orb, 1.0_dp)
                     projector%rcprj = MAX(projector%rcprj, projector%rzetprj(p, l))
                  END DO
               END DO
            END ASSOCIATE
         END IF
         CPASSERT(ASSOCIATED(owner%kinds(kind)%paw_proj_set))
         center_radius(kind) = owner%kinds(kind)%paw_proj_set%rcprj
      END DO
      CALL pair_radius_setup(present_basis, centers, basis_radius, center_radius, pair_radius)
      CALL section_vals_val_get(qs_env%input, "DFT%SUBCELLS", r_val=subcells)
      ALLOCATE (atom2d(nkind))
      CALL atom2d_build(atom2d, distribution_1d, distribution_2d, atomic_kinds, molecules, .FALSE., particles)
      NULLIFY (sap)
      CALL build_neighbor_lists(sap, particles, atom2d, cell, pair_radius, &
                                subcells=subcells, operator_type="ABBA", nlname="MIXED-GAUSSIAN-CENTERS")
      CALL atom2d_cleanup(atom2d)
      CALL create_oce_set(owner%oce)
      CALL allocate_oce_set(owner%oce, nkind)
      CALL build_oce_matrices(owner%oce%intac, .TRUE., 1, owner%kinds, particles, sap, &
                              control%qs_control%gapw_control%eps_fit)
      CALL release_neighbor_list_sets(sap)
   END SUBROUTINE prepare_centers

END MODULE qs_upf_inverse_gaussian_centers
