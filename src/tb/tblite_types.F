!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief types for tblite
!> \author JVP
!> \history creation 09.2024
! **************************************************************************************************

MODULE tblite_types

#if defined(__TBLITE)
   USE mctc_io, ONLY: structure_type
#if defined(__TBLITE_GXTB)
   USE tblite_acp_cache, ONLY: acp_cache
   USE tblite_basis_cache, ONLY: basis_cache
   USE tblite_stream_plan, ONLY: cp2k_stream_plan_type
   USE tblite_bvk_plan, ONLY: cp2k_bvk_plan_type
   USE tblite_projector_cache, ONLY: cp2k_acp_projector_cache_type
   USE tblite_model_bridge, ONLY: cp2k_potential_mixer_type
#endif
   USE tblite_xtb_calculator, ONLY: xtb_calculator
   USE tblite_container_cache, ONLY: container_cache
   USE tblite_scc_mixer, ONLY: cp2k_tblite_broyden_mixer_type
   USE tblite_wavefunction_type, ONLY: wavefunction_type
   USE tblite_scf_potential, ONLY: potential_type
   USE tblite_param, ONLY: param_record
#endif

   USE cp_dbcsr_api, ONLY: dbcsr_p_type
   USE cp_dbcsr_operations, ONLY: dbcsr_deallocate_matrix_set
   USE input_constants, ONLY: tblite_gxtb_acceleration_auto, &
                              tblite_gxtb_acp_cache_enabled, &
                              tblite_gxtb_acp_mesh_dense, &
                              tblite_gxtb_exchange_complete_array, &
                              tblite_gxtb_forward_qualification_off, &
                              tblite_gxtb_gradient_dense, &
                              tblite_gxtb_image_batch_default, &
                              tblite_gxtb_star_dense, &
                              tblite_gxtb_transform_dense, &
                              tblite_mixer_damping_default, &
                              tblite_mixer_iterations_default, &
                              tblite_mixer_max_weight_default, &
                              tblite_mixer_min_weight_default, &
                              tblite_mixer_omega0_default, &
                              tblite_mixer_weight_factor_default
   USE kinds, ONLY: dp
   USE kpoint_lattice_fft, ONLY: kpoint_mesh_certificate_type

#include "../base/base_uses.f90"
   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'tblite_types'

   PUBLIC :: tblite_type, deallocate_tblite_type, allocate_tblite_type
#if defined(__TBLITE_GXTB)
   PUBLIC :: tb_gxtb_kmesh_map_type

   TYPE tb_gxtb_kmesh_map_type
      CHARACTER(LEN=32)                                :: scheme = ""
      INTEGER                                          :: nfull = 0, nred = 0
      LOGICAL                                          :: explicit_general = .FALSE.
      LOGICAL                                          :: full_grid = .FALSE., gamma_centered = .FALSE.
      LOGICAL                                          :: identity = .FALSE., inversion_only = .FALSE.
      INTEGER                                          :: symmetry_backend = 0, symmetry_method = 0
      INTEGER, DIMENSION(3)                            :: nmesh = 1
      REAL(KIND=dp), DIMENSION(3)                      :: kp_shift = 0.0_dp
      INTEGER, ALLOCATABLE, DIMENSION(:)               :: op_count, representative, source, &
                                                          star_member, star_offset, symop
      INTEGER, ALLOCATABLE, DIMENSION(:, :)            :: rotation_slot, symops
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)      :: reduced_xkp, xkp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)         :: reduced_weight, weight, fold_factor
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :) :: atom_phase
      ! Exact key for the consumer-owned whole-mesh Fourier certificate.
      ! This certifies only the reciprocal grid, never a geometry/charge kernel.
      TYPE(kpoint_mesh_certificate_type)               :: exchange_mesh_certificate
   END TYPE tb_gxtb_kmesh_map_type
#endif

   TYPE tblite_type

      LOGICAL                                            :: use_virial = .FALSE.
      INTEGER, ALLOCATABLE                               :: el_num(:)

      REAL(KIND=dp)                                      :: accuracy = 1.0_dp
      REAL(KIND=dp), DIMENSION(3, 3)                     :: sigma = -1.0_dp
      ! tblite contributions are currently replicated on every MPI rank.
      ! Force/virial assembly accounts for CP2K's subsequent MPI reductions.
      REAL(KIND=dp), ALLOCATABLE                         :: e_hal(:)
      REAL(KIND=dp), ALLOCATABLE                         :: e_rep(:)
      REAL(KIND=dp), ALLOCATABLE                         :: e_rep_sc(:)
      REAL(KIND=dp), ALLOCATABLE                         :: e_disp(:)
      REAL(KIND=dp), ALLOCATABLE                         :: e_scd(:)
      REAL(KIND=dp), ALLOCATABLE                         :: e_es(:)
      REAL(KIND=dp), ALLOCATABLE                         :: e_int(:)
      REAL(KIND=dp), ALLOCATABLE                         :: e_inc(:)
      REAL(KIND=dp), ALLOCATABLE                         :: e_exchange(:)
      REAL(KIND=dp), ALLOCATABLE                         :: e_spin(:)
      REAL(KIND=dp), ALLOCATABLE                         :: selfenergy(:)
      REAL(KIND=dp), ALLOCATABLE                         :: cn(:)

      REAL(KIND=dp), ALLOCATABLE                         :: grad(:, :)
      !derivatives w. r. t. coordination number
      REAL(KIND=dp), ALLOCATABLE                         :: dsedcn(:)
      REAL(KIND=dp), ALLOCATABLE                         :: dcndr(:, :, :)
      REAL(KIND=dp), ALLOCATABLE                         :: dcndL(:, :, :)

      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: dipbra => Null(), &
                                                            dipket => Null()
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: quadbra => Null(), &
                                                            quadket => Null()
      ! Spin-resolved density snapshot used by the UKS/LSD force path.
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rho_ao_kp_ref => Null()
      INTEGER                                            :: mixer_memory = tblite_mixer_iterations_default
      REAL(KIND=dp)                                      :: mixer_damping = tblite_mixer_damping_default
      REAL(KIND=dp)                                      :: mixer_max_weight = tblite_mixer_max_weight_default
      REAL(KIND=dp)                                      :: mixer_min_weight = tblite_mixer_min_weight_default
      REAL(KIND=dp)                                      :: mixer_omega0 = tblite_mixer_omega0_default
      REAL(KIND=dp)                                      :: mixer_weight_factor = tblite_mixer_weight_factor_default

#if defined(__TBLITE_GXTB)
      TYPE(basis_cache)                                  :: bcache = basis_cache()
      TYPE(acp_cache)                                    :: acache = acp_cache()
      TYPE(cp2k_acp_projector_cache_type)                 :: acp_projector_cache
      TYPE(cp2k_stream_plan_type) :: stream_plan
      TYPE(cp2k_bvk_plan_type)                            :: bvk_plan
      REAL(KIND=dp), ALLOCATABLE                         :: overlap_dense(:, :)
      TYPE(tb_gxtb_kmesh_map_type), ALLOCATABLE          :: gxtb_kmesh_map
      COMPLEX(KIND=dp), ALLOCATABLE                      :: acp_fock_kp(:, :, :)
      COMPLEX(KIND=dp), ALLOCATABLE                      :: qualification_acp_fock_full(:, :, :)
      COMPLEX(KIND=dp), ALLOCATABLE                      :: exchange_fock_kp(:, :, :, :)
      COMPLEX(KIND=dp), ALLOCATABLE                      :: mixed_fock_kp(:, :, :, :)
      COMPLEX(KIND=dp), ALLOCATABLE                      :: previous_density_kp(:, :, :, :)
      LOGICAL                                            :: gxtb_density_map_validated = .FALSE.
      LOGICAL                                            :: gxtb_transform_validated = .FALSE.
      LOGICAL                                            :: acp_fock_kp_valid = .FALSE.
      LOGICAL                                            :: exchange_fock_kp_valid = .FALSE.
      LOGICAL                                            :: gxtb_image_batch_explicit = .FALSE.
      LOGICAL                                            :: gxtb_acceleration_reported = .FALSE.
      LOGICAL                                            :: gxtb_auto_native_mixer = .FALSE.
      LOGICAL                                            :: mixed_fock_kp_valid = .FALSE.
      LOGICAL                                            :: previous_density_kp_valid = .FALSE.
      REAL(KIND=dp)                                      :: e_acp_kp = 0.0_dp
      INTEGER                                            :: gxtb_acceleration_mode = &
                                                            tblite_gxtb_acceleration_auto
      INTEGER                                            :: gxtb_exchange_stream_backend = &
                                                            tblite_gxtb_exchange_complete_array
      INTEGER                                            :: gxtb_forward_fock_qualification = &
                                                            tblite_gxtb_forward_qualification_off
      INTEGER                                            :: gxtb_gradient_contraction = tblite_gxtb_gradient_dense
      INTEGER                                            :: gxtb_fourier_transform_backend = tblite_gxtb_transform_dense
      INTEGER                                            :: gxtb_symmetry_star_contraction = tblite_gxtb_star_dense
      INTEGER                                            :: gxtb_image_batch_size = tblite_gxtb_image_batch_default
      INTEGER                                            :: gxtb_acp_projector_cache = tblite_gxtb_acp_cache_enabled
      INTEGER                                            :: gxtb_acp_mesh_contraction = tblite_gxtb_acp_mesh_dense
      INTEGER                                            :: potential_mixer_dimension = 0
      TYPE(cp2k_potential_mixer_type), ALLOCATABLE       :: potential_mixer
      TYPE(wavefunction_type)                            :: wfn_aux = wavefunction_type()
      TYPE(container_cache)                              :: ecache = container_cache()
      TYPE(container_cache)                              :: rcache = container_cache()
      TYPE(container_cache)                              :: scache = container_cache()
#endif
#if defined(__TBLITE)
      TYPE(structure_type)                               :: mol = structure_type()
      TYPE(xtb_calculator)                               :: calc = xtb_calculator()
      ! Geometry-dependent caches: Coulomb, dispersion and additional interactions.
      ! Refreshed by tb_init_ham, reused for SCC potentials, energies and gradients.
      TYPE(container_cache)                              :: cache = container_cache()
      TYPE(container_cache)                              :: dcache = container_cache()
      TYPE(container_cache)                              :: icache = container_cache()
      TYPE(wavefunction_type)                            :: wfn = wavefunction_type()
      TYPE(potential_type)                               :: pot = potential_type()
      TYPE(cp2k_tblite_broyden_mixer_type), ALLOCATABLE  :: mixer
      TYPE(param_record), ALLOCATABLE                    :: param
#endif

   END TYPE tblite_type

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param tb_tblite ...
! **************************************************************************************************
   SUBROUTINE allocate_tblite_type(tb_tblite)

      TYPE(tblite_type), POINTER                         :: tb_tblite

      CALL deallocate_tblite_type(tb_tblite)

      ALLOCATE (tb_tblite)

   END SUBROUTINE allocate_tblite_type

! **************************************************************************************************
!> \brief ...
!> \param tb_tblite ...
! **************************************************************************************************
   SUBROUTINE deallocate_tblite_type(tb_tblite)

      TYPE(tblite_type), POINTER                         :: tb_tblite

      IF (ASSOCIATED(tb_tblite)) THEN
         ! Release pointer-owned DBCSR objects; allocatable components are released
         ! automatically when the parent object is deallocated.
         IF (ASSOCIATED(tb_tblite%dipbra)) THEN
            CALL dbcsr_deallocate_matrix_set(tb_tblite%dipbra)
         END IF
         IF (ASSOCIATED(tb_tblite%dipket)) THEN
            CALL dbcsr_deallocate_matrix_set(tb_tblite%dipket)
         END IF
         IF (ASSOCIATED(tb_tblite%quadbra)) THEN
            CALL dbcsr_deallocate_matrix_set(tb_tblite%quadbra)
         END IF
         IF (ASSOCIATED(tb_tblite%quadket)) THEN
            CALL dbcsr_deallocate_matrix_set(tb_tblite%quadket)
         END IF
         IF (ASSOCIATED(tb_tblite%rho_ao_kp_ref)) THEN
            CALL dbcsr_deallocate_matrix_set(tb_tblite%rho_ao_kp_ref)
         END IF
         DEALLOCATE (tb_tblite)
         NULLIFY (tb_tblite)
      END IF

   END SUBROUTINE deallocate_tblite_type

END MODULE tblite_types
