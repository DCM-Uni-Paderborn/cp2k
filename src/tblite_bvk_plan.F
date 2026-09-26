!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Consumer-owned lifetime of geometry-static BvK exchange matrices.
!>        save_tblite still evaluates the model kernels and all responses. No charge-dependent
!>        onsite cache is retained here. Keys use values, not hashes or geometric tolerances.
! **************************************************************************************************
MODULE tblite_bvk_plan
#if defined(__TBLITE_GXTB)
   USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
   USE kinds, ONLY: dp
   USE kpoint_lattice_fft, ONLY: certify_kpoint_mesh, kpoint_mesh_certificate_type
   USE mctc_env, ONLY: error_type, fatal_error
   USE mctc_io, ONLY: structure_type
   USE tblite_exchange_cache, ONLY: exchange_cache
   USE tblite_exchange_fock, ONLY: exchange_bvk_kernel, exchange_fock
   USE tblite_exchange_contractions, ONLY: get_kfock_kmesh, get_kgrad_kmesh
   USE tblite_stream_plan, ONLY: cp2k_stream_plan_type
#if defined(__TBLITE_GXTB_TRANSFORM)
   USE tblite_exchange_fft, ONLY: cp2k_exchange_transform_type, prepare_exchange_transform, &
                                  release_exchange_transform
#endif
#endif
   IMPLICIT NONE
   PRIVATE
#if defined(__TBLITE_GXTB)
   PUBLIC :: cp2k_bvk_plan_type, prepare_bvk_kernel, evaluate_bvk_exchange, evaluate_bvk_gradient, &
             bvk_kernel_snapshot, bvk_plan_counts, invalidate_bvk_plan

   TYPE cp2k_bvk_plan_type
      PRIVATE
      LOGICAL :: valid = .FALSE.
      INTEGER :: builds = 0, hits = 0
      INTEGER, ALLOCATABLE :: integer_key(:)
      REAL(KIND=dp), ALLOCATABLE :: real_key(:)
      TYPE(exchange_bvk_kernel) :: kernel
#if defined(__TBLITE_GXTB_TRANSFORM)
      TYPE(cp2k_exchange_transform_type) :: transform
#endif
   END TYPE cp2k_bvk_plan_type

CONTAINS

! **************************************************************************************************
!> \brief Reuse only when every geometry, representation and static model input is unchanged.
!> \param plan CP2K-owned storage
!> \param exchange live model
!> \param mol current structure
!> \param nmesh BvK replication dimensions
!> \param reused true only for a previously validated exact key
!> \param error invalid descriptor or model kernel
! **************************************************************************************************
   SUBROUTINE prepare_bvk_kernel(plan, exchange, mol, nmesh, reused, error)
      TYPE(cp2k_bvk_plan_type), INTENT(INOUT)            :: plan
      TYPE(exchange_fock), INTENT(IN)                    :: exchange
      TYPE(structure_type), INTENT(IN)                   :: mol
      INTEGER, INTENT(IN)                                :: nmesh(3)
      LOGICAL, INTENT(OUT)                               :: reused
      TYPE(error_type), ALLOCATABLE, INTENT(OUT)         :: error

      INTEGER                                            :: axis, cell, cells, index0
      INTEGER, ALLOCATABLE                               :: integers(:)
      REAL(KIND=dp), ALLOCATABLE                         :: reals(:)

      reused = .FALSE.
      IF (ANY(nmesh < 1) .OR. ANY((.NOT. mol%periodic) .AND. nmesh /= 1)) THEN
         plan%valid = .FALSE.
         CALL fatal_error(error, 'Invalid periodic BvK dimensions')
         RETURN
      END IF
      cells = 1
      DO axis = 1, 3
         IF (cells > HUGE(cells)/nmesh(axis)) THEN
            plan%valid = .FALSE.
            CALL fatal_error(error, 'BvK cell count overflows the integer representation')
            RETURN
         END IF
         cells = cells*nmesh(axis)
      END DO
      integers = [mol%nat, mol%nid, nmesh, MERGE(1, 0, mol%periodic), SIZE(mol%id), mol%id, &
                  SIZE(mol%num), mol%num, SHAPE(mol%xyz), exchange%nao, exchange%nsh, exchange%maxsh, &
                  SIZE(exchange%nsh_id), exchange%nsh_id, SIZE(exchange%nao_sh), exchange%nao_sh, &
                  SIZE(exchange%ish_at), exchange%ish_at, SIZE(exchange%iao_sh), exchange%iao_sh, &
                  SHAPE(exchange%hubbard), SHAPE(exchange%offdiag_scale), SHAPE(exchange%onecxints), &
                  SHAPE(exchange%rad), SHAPE(exchange%kq), SHAPE(exchange%corr_scale), SHAPE(exchange%corr_rad)]
      reals = [RESHAPE(mol%xyz, [SIZE(mol%xyz)]), RESHAPE(mol%lattice, [9]), &
               exchange%frscale, exchange%onsite_frscale, exchange%omega, exchange%lrscale, &
               exchange%ondiag_scale, exchange%hubbard_exp, exchange%hubbard_exp_r0, exchange%gexp, exchange%corr_exp, &
               RESHAPE(exchange%hubbard, [SIZE(exchange%hubbard)]), &
               RESHAPE(exchange%offdiag_scale, [SIZE(exchange%offdiag_scale)]), &
               RESHAPE(exchange%onecxints, [SIZE(exchange%onecxints)]), &
               RESHAPE(exchange%rad, [SIZE(exchange%rad)]), RESHAPE(exchange%kq, [SIZE(exchange%kq)]), &
               RESHAPE(exchange%corr_scale, [SIZE(exchange%corr_scale)]), &
               RESHAPE(exchange%corr_rad, [SIZE(exchange%corr_rad)])]
      IF (.NOT. ALL(ieee_is_finite(reals))) THEN
         plan%valid = .FALSE.
         CALL fatal_error(error, 'Nonfinite geometry or model in BvK descriptor')
         RETURN
      END IF
      IF (plan%valid) THEN
         IF (SIZE(integers) == SIZE(plan%integer_key) .AND. SIZE(reals) == SIZE(plan%real_key)) THEN
            reused = ALL(integers == plan%integer_key) .AND. ALL(reals == plan%real_key)
         END IF
      END IF
      IF (reused) THEN
         plan%hits = plan%hits + 1
         RETURN
      END IF
      plan%valid = .FALSE.
      ! The provider owns Wigner-Seitz electrostatics and radial model formulas.
      ! Intent(out) on this unchanged public routine releases the old kernel first.
      CALL exchange%get_bvk_Kmatrix(mol, nmesh, plan%kernel)
      IF (.NOT. ALLOCATED(plan%kernel%reps) .OR. &
          .NOT. ALLOCATED(plan%kernel%g_mulliken_r) .OR. .NOT. ALLOCATED(plan%kernel%g_bocorr_r)) THEN
         CALL fatal_error(error, 'Incomplete BvK model kernel')
         RETURN
      END IF
      IF (ANY(SHAPE(plan%kernel%reps) /= [3, cells]) .OR. &
          ANY(SHAPE(plan%kernel%g_mulliken_r) /= [exchange%nsh, exchange%nsh, cells]) .OR. &
          ANY(SHAPE(plan%kernel%g_bocorr_r) /= [mol%nat, mol%nat, cells])) THEN
         CALL fatal_error(error, 'Inconsistent BvK model kernel dimensions')
         RETURN
      END IF
      IF (ANY(plan%kernel%nmesh /= nmesh)) THEN
         CALL fatal_error(error, 'Inconsistent BvK model mesh')
         RETURN
      END IF
      DO cell = 1, cells
         index0 = cell - 1
         IF (ANY(plan%kernel%reps(:, cell) /= &
                 [MODULO(index0, nmesh(1)), MODULO(index0/nmesh(1), nmesh(2)), &
                  index0/(nmesh(1)*nmesh(2))])) THEN
            CALL fatal_error(error, 'Noncanonical BvK model representatives')
            RETURN
         END IF
      END DO
      IF (.NOT. ALL(ieee_is_finite(plan%kernel%g_mulliken_r)) .OR. &
          .NOT. ALL(ieee_is_finite(plan%kernel%g_bocorr_r))) THEN
         CALL fatal_error(error, 'Nonfinite BvK model kernel')
         RETURN
      END IF
      plan%integer_key = integers
      plan%real_key = reals
      plan%builds = plan%builds + 1
      plan%valid = .TRUE.
   END SUBROUTINE prepare_bvk_kernel

! **************************************************************************************************
!> \brief Evaluate the unchanged model with a static kernel and the caller's current onsite state.
!> \param plan CP2K-owned geometry/model kernel
!> \param certificate CP2K-owned reciprocal-mesh proof
!> \param exchange current model
!> \param mol current structure
!> \param live current provider exchange cache; never retained in the plan
!> \param nmesh BvK dimensions
!> \param kfrac unweighted fractional k points
!> \param weights normalized full-grid weights
!> \param density spin-resolved unweighted Bloch density
!> \param overlap unweighted Bloch overlap
!> \param fock complete unweighted Fock response
!> \param vsh weighted shell potential
!> \param energy exchange energy per primitive cell
!> \param error invalid cache key or mesh
! **************************************************************************************************
   SUBROUTINE evaluate_bvk_exchange(plan, certificate, exchange, mol, live, nmesh, &
                                    kfrac, weights, density, overlap, fock, vsh, energy, error, use_consumer_transform)
      TYPE(cp2k_bvk_plan_type), INTENT(INOUT)            :: plan
      TYPE(kpoint_mesh_certificate_type), INTENT(INOUT)  :: certificate
      TYPE(exchange_fock), INTENT(IN)                    :: exchange
      TYPE(structure_type), INTENT(IN)                   :: mol
      TYPE(exchange_cache), INTENT(IN)                   :: live
      INTEGER, INTENT(IN)                                :: nmesh(3)
      REAL(KIND=dp), INTENT(IN)                          :: kfrac(:, :), weights(:)
      COMPLEX(KIND=dp), INTENT(IN)                       :: density(:, :, :, :), overlap(:, :, :)
      COMPLEX(KIND=dp), INTENT(OUT)                      :: fock(:, :, :, :)
      REAL(KIND=dp), INTENT(OUT)                         :: vsh(:), energy
      TYPE(error_type), ALLOCATABLE, INTENT(OUT)         :: error
      LOGICAL, INTENT(IN), OPTIONAL                      :: use_consumer_transform

      INTEGER                                            :: nao, nk, nspin
      LOGICAL                                            :: reused, use_transform, valid
      TYPE(cp2k_stream_plan_type)                        :: contraction_plan

      nk = SIZE(weights)
      nao = exchange%nao
      nspin = SIZE(density, 3)
      IF (nk < 1 .OR. ANY(SHAPE(kfrac) /= [3, nk]) .OR. &
          ANY(SHAPE(density) /= [nao, nao, nspin, nk]) .OR. &
          ANY(SHAPE(fock) /= [nao, nao, nspin, nk]) .OR. &
          ANY(SHAPE(overlap) /= [nao, nao, nk]) .OR. SIZE(vsh) /= exchange%nsh .OR. &
          (nspin /= 1 .AND. nspin /= 2)) THEN
         CALL fatal_error(error, 'Inconsistent BvK exchange array dimensions')
         RETURN
      END IF
      CALL prepare_bvk_kernel(plan, exchange, mol, nmesh, reused, error)
      IF (ALLOCATED(error)) RETURN
      CALL certify_kpoint_mesh(certificate, kfrac, nmesh, plan%kernel%reps, weights, valid, reused)
      IF (.NOT. valid) THEN
         CALL fatal_error(error, 'Incomplete or nonuniform BvK Fourier mesh')
         RETURN
      END IF
      ! Charge dependence and density/overlap are always evaluated from the current SCC state.
      use_transform = .TRUE.
      IF (PRESENT(use_consumer_transform)) use_transform = use_consumer_transform
#if defined(__TBLITE_GXTB_TRANSFORM)
      IF (use_transform) THEN
         CALL prepare_exchange_transform(plan%transform, kfrac, nmesh, plan%kernel%reps, weights)
         CALL get_kfock_kmesh(exchange, mol, live, contraction_plan, plan%kernel, kfrac, weights, &
                              density, overlap, fock, vsh, energy, transform=plan%transform)
         RETURN
      END IF
#endif
      CALL get_kfock_kmesh(exchange, mol, live, contraction_plan, plan%kernel, kfrac, weights, &
                           density, overlap, fock, vsh, energy)
   END SUBROUTINE evaluate_bvk_exchange

! **************************************************************************************************
!> \brief Reuse geometry kernels and CP2K Fourier maps for the exact contraction adjoint.
!>        Radial-kernel geometry derivatives, including the strain convention, remain in save_tblite.
! **************************************************************************************************
   SUBROUTINE evaluate_bvk_gradient(plan, certificate, exchange, mol, live, nmesh, &
                                    kfrac, weights, density, overlap, overlap_grad, gradient, sigma, error)
      TYPE(cp2k_bvk_plan_type), INTENT(INOUT)            :: plan
      TYPE(kpoint_mesh_certificate_type), INTENT(INOUT)  :: certificate
      TYPE(exchange_fock), INTENT(IN)                    :: exchange
      TYPE(structure_type), INTENT(IN)                   :: mol
      TYPE(exchange_cache), INTENT(IN)                   :: live
      INTEGER, INTENT(IN)                                :: nmesh(3)
      REAL(KIND=dp), INTENT(IN)                          :: kfrac(:, :), weights(:)
      COMPLEX(KIND=dp), INTENT(IN)                       :: density(:, :, :, :), overlap(:, :, :)
      COMPLEX(KIND=dp), INTENT(OUT)                      :: overlap_grad(:, :, :)
      REAL(KIND=dp), INTENT(OUT)                         :: gradient(:, :), sigma(:, :)
      TYPE(error_type), ALLOCATABLE, INTENT(OUT)         :: error

      INTEGER                                            :: nao, nk, nspin
      LOGICAL                                            :: reused, valid
      TYPE(cp2k_stream_plan_type)                        :: contraction_plan
      REAL(KIND=dp), ALLOCATABLE                         :: bocorr_grad(:, :, :), &
                                                            mulliken_grad(:, :, :)

      gradient = 0.0_dp
      sigma = 0.0_dp
      overlap_grad = CMPLX(0.0_dp, 0.0_dp, dp)
      nk = SIZE(weights)
      nao = exchange%nao
      nspin = SIZE(density, 3)
      IF (nk < 1 .OR. ANY(SHAPE(kfrac) /= [3, nk]) .OR. &
          ANY(SHAPE(density) /= [nao, nao, nspin, nk]) .OR. &
          ANY(SHAPE(overlap) /= [nao, nao, nk]) .OR. &
          ANY(SHAPE(overlap_grad) /= [nao, nao, nk]) .OR. &
          ANY(SHAPE(gradient) /= [3, mol%nat]) .OR. ANY(SHAPE(sigma) /= [3, 3]) .OR. &
          (nspin /= 1 .AND. nspin /= 2)) THEN
         CALL fatal_error(error, 'Inconsistent BvK gradient array dimensions')
         RETURN
      END IF
      CALL prepare_bvk_kernel(plan, exchange, mol, nmesh, reused, error)
      IF (ALLOCATED(error)) RETURN
      CALL certify_kpoint_mesh(certificate, kfrac, nmesh, plan%kernel%reps, weights, valid, reused)
      IF (.NOT. valid) THEN
         CALL fatal_error(error, 'Incomplete or nonuniform BvK Fourier mesh')
         RETURN
      END IF
      ALLOCATE (mulliken_grad(nao, nao, nk), bocorr_grad(nao, nao, nk))
#if defined(__TBLITE_GXTB_TRANSFORM)
      CALL prepare_exchange_transform(plan%transform, kfrac, nmesh, plan%kernel%reps, weights)
      CALL get_kgrad_kmesh(exchange, mol, live, contraction_plan, plan%kernel, kfrac, weights, &
                           density, overlap, overlap_grad, mulliken_grad, bocorr_grad, transform=plan%transform)
#else
      CALL get_kgrad_kmesh(exchange, mol, live, contraction_plan, plan%kernel, kfrac, weights, &
                           density, overlap, overlap_grad, mulliken_grad, bocorr_grad)
#endif
      CALL exchange%get_bvk_Kmatrix_derivs(mol, plan%kernel, mulliken_grad, bocorr_grad, gradient, sigma)
   END SUBROUTINE evaluate_bvk_gradient

! **************************************************************************************************
!> \brief Independent test/debug snapshot; production evaluation uses the private kernel in place.
!> \param plan CP2K-owned storage
!> \param kernel value copy, never an externally mutable alias
!> \param error no valid kernel available
! **************************************************************************************************
   SUBROUTINE bvk_kernel_snapshot(plan, kernel, error)
      TYPE(cp2k_bvk_plan_type), INTENT(IN)               :: plan
      TYPE(exchange_bvk_kernel), INTENT(OUT)             :: kernel
      TYPE(error_type), ALLOCATABLE, INTENT(OUT)         :: error

      IF (.NOT. plan%valid) THEN
         CALL fatal_error(error, 'No valid BvK kernel snapshot')
         RETURN
      END IF
      kernel = plan%kernel
   END SUBROUTINE bvk_kernel_snapshot

! **************************************************************************************************
!> \brief Invalidate reuse without discarding diagnostic counts.
!> \param plan CP2K-owned storage
! **************************************************************************************************
   SUBROUTINE invalidate_bvk_plan(plan)
      TYPE(cp2k_bvk_plan_type), INTENT(INOUT)            :: plan

      plan%valid = .FALSE.
#if defined(__TBLITE_GXTB_TRANSFORM)
      CALL release_exchange_transform(plan%transform)
#endif
   END SUBROUTINE invalidate_bvk_plan

! **************************************************************************************************
!> \brief Report actual rebuilds and exact-key hits.
!> \param plan CP2K-owned storage
!> \param builds successful provider evaluations
!> \param hits exact-key matches
!> \param valid reusable kernel
! **************************************************************************************************
   SUBROUTINE bvk_plan_counts(plan, builds, hits, valid)
      TYPE(cp2k_bvk_plan_type), INTENT(IN)               :: plan
      INTEGER, INTENT(OUT)                               :: builds, hits
      LOGICAL, INTENT(OUT)                               :: valid

      builds = plan%builds
      hits = plan%hits
      valid = plan%valid
   END SUBROUTINE bvk_plan_counts
#endif
END MODULE tblite_bvk_plan
