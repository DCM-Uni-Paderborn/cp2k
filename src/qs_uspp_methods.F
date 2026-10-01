!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Connect native USPP state to the Quickstep geometry and neighbor-list lifecycle.
! **************************************************************************************************
MODULE qs_uspp_methods
   USE cp_dbcsr_api,                    ONLY: dbcsr_add,&
                                              dbcsr_set,&
                                              dbcsr_type
   USE kinds,                           ONLY: dp
   USE particle_types,                  ONLY: particle_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type,&
                                              set_ks_env
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_uspp_types,                   ONLY: qs_uspp_init,&
                                              qs_uspp_type,&
                                              qs_uspp_update
   USE sap_kind_types,                  ONLY: build_sap_ints,&
                                              release_sap_int,&
                                              sap_int_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_uspp_methods'
   PUBLIC :: qs_uspp_prepare, qs_uspp_init_ao_matrix

CONTAINS

! **************************************************************************************************
!> \brief Prepare USPP centers, unweighted projector overlaps, and the generalized metric.
!>        Called after rebuilding the Gaussian overlap. The NC path allocates nothing.
!> \param ks_env KS environment with current kinds, particles and ABBA projector neighbor list.
!> \param gaussian_overlap Current Gaussian overlap; never modified or aliased by the metric.
! **************************************************************************************************
   SUBROUTINE qs_uspp_prepare(ks_env, gaussian_overlap)
      TYPE(qs_ks_env_type), INTENT(INOUT)                :: ks_env
      TYPE(dbcsr_type), INTENT(IN)                       :: gaussian_overlap

      INTEGER                                            :: handle
      LOGICAL                                            :: do_kpoints
      TYPE(neighbor_list_set_p_type), POINTER            :: sap_ppnl(:)
      TYPE(particle_type), POINTER                       :: particle_set(:)
      TYPE(qs_kind_type), POINTER                        :: qs_kind_set(:)
      TYPE(qs_uspp_type), POINTER                        :: uspp
      TYPE(sap_int_type), POINTER                        :: sap_int(:)

      CALL get_ks_env(ks_env, uspp=uspp, qs_kind_set=qs_kind_set, particle_set=particle_set)
      CALL qs_uspp_init(uspp, qs_kind_set, particle_set)
      CALL set_ks_env(ks_env, uspp=uspp)
      IF (.NOT. ASSOCIATED(uspp)) RETURN
      CALL timeset('qs_uspp_prepare', handle)
      CALL get_ks_env(ks_env, do_kpoints=do_kpoints, sap_ppnl=sap_ppnl)
      CPASSERT(.NOT. do_kpoints)
      CPASSERT(ASSOCIATED(sap_ppnl))
      ALLOCATE (sap_int(SIZE(qs_kind_set)**2))
      CALL build_sap_ints(sap_int, sap_ppnl, qs_kind_set, nder=0)
      CALL qs_uspp_update(uspp, gaussian_overlap, sap_int, particle_set)
      CALL release_sap_int(sap_int)
      CALL timestop(handle)
   END SUBROUTINE qs_uspp_prepare

! **************************************************************************************************
!> \brief Initialize a matrix on the union of its ordinary and USPP projector-mediated blocks.
!>        Called after ordinary neighbor-list allocation, before assigning matrix values.
!> \param ks_env KS environment with current USPP geometry data.
!> \param matrix Newly allocated AO matrix, zeroed after expanding its block pattern.
! **************************************************************************************************
   SUBROUTINE qs_uspp_init_ao_matrix(ks_env, matrix)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix

      TYPE(qs_uspp_type), POINTER                        :: uspp

      CALL get_ks_env(ks_env, uspp=uspp)
      IF (.NOT. ASSOCIATED(uspp)) RETURN
      CPASSERT(ASSOCIATED(uspp%ao_pattern))
      CALL dbcsr_add(matrix, uspp%ao_pattern, 1.0_dp, 1.0_dp)
      CALL dbcsr_set(matrix, 0.0_dp)
   END SUBROUTINE qs_uspp_init_ao_matrix

END MODULE qs_uspp_methods
