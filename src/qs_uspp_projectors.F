!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Assemble the Gamma-point USPP projector matrix from Quickstep overlap lists.
! **************************************************************************************************
MODULE qs_uspp_projectors
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_create, dbcsr_distribution_get, dbcsr_distribution_new, dbcsr_distribution_release, &
        dbcsr_distribution_type, dbcsr_finalize, dbcsr_get_block_p, dbcsr_get_info, &
        dbcsr_get_stored_coordinates, dbcsr_reserve_blocks, dbcsr_set, dbcsr_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_comm_type
   USE sap_kind_types,                  ONLY: sap_int_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_uspp_projectors'
   PUBLIC :: qs_uspp_projector_matrix, qs_uspp_projector_gradient

CONTAINS

! **************************************************************************************************
!> \brief Sum all periodic images of <chi_mu|beta_i> into a rectangular DBCSR matrix.
!>        The row distribution must be that used to build the ABBA neighbor list.
!>        ABBA supplies every projector neighbor on each rank owning the AO's row or
!>        column. Selecting B's unique block owner therefore removes replicated work
!>        without discarding images. No filter or Gaussian-overlap sparsity is imposed.
!> \param template AO matrix whose atom blocks and row distribution generated sap_int.
!> \param sap_int Uncompressed overlaps from build_sap_ints (not H-weighted achint).
!> \param center_atoms Global atom indices for the USPP centers, in column-block order.
!> \param projector_sizes Magnetic-projector count per selected center.
!> \param overlaps New Gamma-point B matrix. Non-Gamma phases require a separate path.
! **************************************************************************************************
   SUBROUTINE qs_uspp_projector_matrix(template, sap_int, center_atoms, projector_sizes, overlaps)
      TYPE(dbcsr_type), INTENT(IN)                       :: template
      TYPE(sap_int_type), INTENT(IN)                     :: sap_int(:)
      INTEGER, INTENT(IN)                                :: center_atoms(:), projector_sizes(:)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: overlaps

      INTEGER                                            :: atom, center, i, j, k, npcols, owner, &
                                                            pass
      INTEGER, ALLOCATABLE                               :: atom_to_center(:)
      INTEGER, POINTER                                   :: col_dist(:), col_sizes(:), row_dist(:), &
                                                            row_sizes(:)
      LOGICAL                                            :: found
      REAL(KIND=dp), POINTER                             :: block(:, :)
      TYPE(dbcsr_distribution_type)                      :: dist, projector_dist
      TYPE(mp_comm_type)                                 :: group

      CALL dbcsr_get_info(template, distribution=dist, row_blk_size=row_sizes, group=group)
      CALL dbcsr_distribution_get(dist, row_dist=row_dist, npcols=npcols)
      CPASSERT(SIZE(center_atoms) == SIZE(projector_sizes))
      CPASSERT(SIZE(center_atoms) > 0 .AND. ALL(projector_sizes > 0))
      CPASSERT(ALL(center_atoms >= 1) .AND. ALL(center_atoms <= SIZE(row_sizes)))
      ALLOCATE (atom_to_center(SIZE(row_sizes)), col_sizes(SIZE(center_atoms)), col_dist(SIZE(center_atoms)))
      atom_to_center = 0
      DO i = 1, SIZE(center_atoms)
         CPASSERT(atom_to_center(center_atoms(i)) == 0)
         atom_to_center(center_atoms(i)) = i
         col_dist(i) = MODULO(i - 1, npcols)
      END DO
      col_sizes(:) = projector_sizes
      CALL dbcsr_distribution_new(projector_dist, template=dist, row_dist=row_dist, col_dist=col_dist)
      CALL dbcsr_create(overlaps, 'USPP Gaussian projector overlaps', projector_dist, 'N', row_sizes, col_sizes)
      CALL dbcsr_distribution_release(projector_dist)
      DEALLOCATE (col_sizes, col_dist)
      DO pass = 1, 2
      DO i = 1, SIZE(sap_int)
         IF (.NOT. ASSOCIATED(sap_int(i)%alist)) CYCLE
         DO j = 1, SIZE(sap_int(i)%alist)
            IF (.NOT. ASSOCIATED(sap_int(i)%alist(j)%clist)) CYCLE
            atom = sap_int(i)%alist(j)%aatom
            CPASSERT(atom >= 1 .AND. atom <= SIZE(row_sizes))
            DO k = 1, SIZE(sap_int(i)%alist(j)%clist)
               ! Unused entries have catom=0 and no integral storage.
               center = sap_int(i)%alist(j)%clist(k)%catom
               IF (center == 0) CYCLE
               CPASSERT(center >= 1 .AND. center <= SIZE(atom_to_center))
               center = atom_to_center(center)
               IF (center == 0) CYCLE
               CALL dbcsr_get_stored_coordinates(overlaps, atom, center, owner)
               IF (owner /= group%mepos) CYCLE
               CPASSERT(ASSOCIATED(sap_int(i)%alist(j)%clist(k)%acint))
               ASSOCIATE (integrals => sap_int(i)%alist(j)%clist(k)%acint)
                  CPASSERT(SIZE(integrals, 1) == row_sizes(atom))
                  CPASSERT(SIZE(integrals, 2) == projector_sizes(center))
                  CPASSERT(SIZE(integrals, 3) >= 1)
                  IF (pass == 1) THEN
                     CALL dbcsr_reserve_blocks(overlaps, rows=[atom], cols=[center])
                  ELSE
                     CALL dbcsr_get_block_p(overlaps, atom, center, block, found)
                     CPASSERT(found)
                     block(:, :) = block + integrals(:, :, 1)
                  END IF
               END ASSOCIATE
            END DO
         END DO
      END DO
      IF (pass == 1) THEN
         CALL dbcsr_finalize(overlaps)
         CALL dbcsr_set(overlaps, 0.0_dp)
      END IF
      END DO
   END SUBROUTINE qs_uspp_projector_matrix

! **************************************************************************************************
!> \brief Contract dE/dB with all image-resolved AO-center projector derivatives.
!>        Returns energy gradients, with equal and opposite AO/projector-center contributions.
!>        The unique owner of each rectangular block supplies its contribution; no MPI sum here.
!> \param gradient Distributed dE/dB, with the same row/column distribution as B.
!> \param sap_int Uncompressed nder=1 projector overlaps from the ABBA neighbor list.
!> \param center_atoms Global atom indices for the projector column blocks.
!> \param atom_gradient Local energy gradients (3,natom), incremented in place.
! **************************************************************************************************
   SUBROUTINE qs_uspp_projector_gradient(gradient, sap_int, center_atoms, atom_gradient)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: gradient
      TYPE(sap_int_type), INTENT(IN)                     :: sap_int(:)
      INTEGER, INTENT(IN)                                :: center_atoms(:)
      REAL(KIND=dp), INTENT(INOUT)                       :: atom_gradient(:, :)

      INTEGER                                            :: atom, catom, center, i, idir, j, k, owner
      INTEGER, ALLOCATABLE                               :: atom_to_center(:)
      LOGICAL                                            :: found
      REAL(KIND=dp)                                      :: value
      REAL(KIND=dp), POINTER                             :: block(:, :)
      TYPE(mp_comm_type)                                 :: group

      CPASSERT(SIZE(atom_gradient, 1) == 3)
      ALLOCATE (atom_to_center(SIZE(atom_gradient, 2)))
      atom_to_center = 0
      DO center = 1, SIZE(center_atoms)
         atom_to_center(center_atoms(center)) = center
      END DO
      CALL dbcsr_get_info(gradient, group=group)
      DO i = 1, SIZE(sap_int)
         IF (.NOT. ASSOCIATED(sap_int(i)%alist)) CYCLE
         DO j = 1, SIZE(sap_int(i)%alist)
            IF (.NOT. ASSOCIATED(sap_int(i)%alist(j)%clist)) CYCLE
            atom = sap_int(i)%alist(j)%aatom
            DO k = 1, SIZE(sap_int(i)%alist(j)%clist)
               catom = sap_int(i)%alist(j)%clist(k)%catom
               IF (catom == 0) CYCLE
               CPASSERT(catom >= 1 .AND. catom <= SIZE(atom_to_center))
               center = atom_to_center(catom)
               IF (center == 0) CYCLE
               CALL dbcsr_get_stored_coordinates(gradient, atom, center, owner)
               IF (owner /= group%mepos) CYCLE
               CALL dbcsr_get_block_p(gradient, atom, center, block, found)
               IF (.NOT. found) CYCLE
               CPASSERT(ASSOCIATED(sap_int(i)%alist(j)%clist(k)%acint))
               ASSOCIATE (integrals => sap_int(i)%alist(j)%clist(k)%acint)
                  CPASSERT(SIZE(integrals, 3) >= 4)
                  DO idir = 1, 3
                     value = SUM(block*integrals(:, :, idir + 1))
                     atom_gradient(idir, atom) = atom_gradient(idir, atom) + value
                     atom_gradient(idir, catom) = atom_gradient(idir, catom) - value
                  END DO
               END ASSOCIATE
            END DO
         END DO
      END DO
   END SUBROUTINE qs_uspp_projector_gradient

END MODULE qs_uspp_projectors
