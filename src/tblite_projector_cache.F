!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Consumer-owned reuse of compact ACP projectors, not of integral or derivative kernels.
!>        An exact key includes effective (possibly charge-dependent) contraction coefficients.
!>        The provider continues to build all integrals on misses and all derivative integrals.
! **************************************************************************************************
MODULE tblite_projector_cache
#if defined(__TBLITE_GXTB)
   USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
   USE kinds, ONLY: dp
   USE mctc_io, ONLY: structure_type
   USE tblite_acp_cache, ONLY: acp_cache
   USE tblite_acp_type, ONLY: acp_type
   USE tblite_acp_images, ONLY: acp_image_type
   USE tblite_basis_cache, ONLY: basis_cache
   USE tblite_basis_type, ONLY: basis_type
#endif
   IMPLICIT NONE
   PRIVATE
#if defined(__TBLITE_GXTB)
   PUBLIC :: cp2k_acp_projector_cache_type, projector_cache_lookup, projector_cache_store, &
             projector_cache_bypass, projector_cache_counts

   TYPE cp2k_acp_projector_cache_type
      PRIVATE
      LOGICAL :: valid = .FALSE., pending = .FALSE.
      INTEGER :: hits = 0, misses = 0, bypasses = 0
      INTEGER, ALLOCATABLE :: integer_key(:)
      REAL(KIND=dp), ALLOCATABLE :: real_key(:)
      REAL(KIND=dp), ALLOCATABLE :: translation(:, :), overlap(:, :, :), level(:)
   END TYPE cp2k_acp_projector_cache_type

CONTAINS

! **************************************************************************************************
!> \brief Describe the actual basis representation consumed by the projector integral builder.
!> \param mol structure
!> \param basis valence or auxiliary basis
!> \param live current atom-resolved basis cache
!> \param integers exact shell layout and primitive descriptors
!> \param reals cutoffs, primitive exponents and effective contraction coefficients
! **************************************************************************************************
   SUBROUTINE basis_key(mol, basis, live, integers, reals)
      TYPE(structure_type), INTENT(IN)                   :: mol

      CLASS(basis_type), INTENT(IN) :: basis
      TYPE(basis_cache), INTENT(IN) :: live
      INTEGER, ALLOCATABLE, INTENT(OUT) :: integers(:)
      REAL(KIND=dp), ALLOCATABLE, INTENT(OUT) :: reals(:)
      INTEGER :: atom, species, shell, offset, total, count, base

      integers = [basis%nsh, basis%nao, basis%maxl, SIZE(basis%nsh_id), basis%nsh_id, &
                  SIZE(basis%nsh_at), basis%nsh_at, SIZE(basis%ish_at), basis%ish_at, &
                  SIZE(basis%iao_sh), basis%iao_sh, SIZE(basis%nao_sh), basis%nao_sh]
      base = SIZE(integers)
      total = 0
      count = 0
      DO atom = 1, mol%nat
         species = mol%id(atom)
         DO shell = 1, basis%nsh_id(species)
            total = total + basis%cgto(shell, species)%raw%nprim
            count = count + 1
         END DO
      END DO
      integers = [integers, SPREAD(0, 1, 2*count)]
      ALLOCATE (reals(1 + 2*total))
      reals(1) = basis%intcut
      offset = 1
      count = base
      DO atom = 1, mol%nat
         species = mol%id(atom)
         DO shell = 1, basis%nsh_id(species)
            ASSOCIATE (cgto => basis%cgto(shell, species)%raw)
               integers(count + 1:count + 2) = [cgto%ang, cgto%nprim]
               reals(offset + 1:offset + cgto%nprim) = cgto%alpha(1:cgto%nprim)
               CALL cgto%get_coeffs(live%cgto(shell, atom), &
                                    reals(offset + cgto%nprim + 1:offset + 2*cgto%nprim))
               offset = offset + 2*cgto%nprim
               count = count + 2
            END ASSOCIATE
         END DO
      END DO
   END SUBROUTINE basis_key

! **************************************************************************************************
!> \brief Look up compact projectors, refreshing live auxiliary-basis state even on a cache hit.
!> \param cache consumer-owned storage
!> \param mol current structure
!> \param basis valence basis
!> \param live current valence-basis cache
!> \param acp provider model
!> \param auxiliary live provider ACP cache, never replaced by an old derivative cache
!> \param translations ordered input translation list
!> \param images compact result on a hit
!> \param hit exact-key match with a reusable payload
! **************************************************************************************************
   SUBROUTINE projector_cache_lookup(cache, mol, basis, live, acp, auxiliary, translations, images, hit)
      TYPE(cp2k_acp_projector_cache_type), INTENT(INOUT) :: cache
      TYPE(structure_type), INTENT(IN)                   :: mol

      CLASS(basis_type), INTENT(IN) :: basis
      TYPE(basis_cache), INTENT(IN) :: live
      TYPE(acp_type), INTENT(IN) :: acp
      TYPE(acp_cache), INTENT(INOUT) :: auxiliary
      REAL(KIND=dp), INTENT(IN) :: translations(:, :)
      TYPE(acp_image_type), INTENT(OUT) :: images
      LOGICAL, INTENT(OUT) :: hit
      INTEGER, ALLOCATABLE :: integers(:), valence_i(:), auxiliary_i(:)
      REAL(KIND=dp), ALLOCATABLE :: reals(:), valence_r(:), auxiliary_r(:)
      INTEGER :: image, ao

      hit = .FALSE.
      CALL acp%update(mol, auxiliary)
      ! The minimal provider assumes an already allocated intermediate has the right shape.
      ! Make both cache hits and direct rebuilds safe across calculator/geometry-size changes.
      IF (ALLOCATED(auxiliary%scaled_pv_overlap)) THEN
         IF (ANY(SHAPE(auxiliary%scaled_pv_overlap) /= [acp%auxbas%nao, basis%nao])) &
            DEALLOCATE (auxiliary%scaled_pv_overlap)
      END IF
      CALL basis_key(mol, basis, live, valence_i, valence_r)
      CALL basis_key(mol, acp%auxbas, auxiliary%auxbas, auxiliary_i, auxiliary_r)
      integers = [mol%nat, SIZE(mol%id), mol%id, SIZE(mol%num), mol%num, &
                  MERGE(1, 0, mol%periodic), SHAPE(mol%xyz), SHAPE(translations), SHAPE(acp%levels), &
                  SIZE(valence_i), valence_i, SIZE(auxiliary_i), auxiliary_i]
      reals = [RESHAPE(mol%xyz, [SIZE(mol%xyz)]), RESHAPE(mol%lattice, [9]), &
               RESHAPE(translations, [SIZE(translations)]), RESHAPE(acp%levels, [SIZE(acp%levels)]), &
               valence_r, auxiliary_r]
      IF (.NOT. ALL(ieee_is_finite(reals))) THEN
         cache%valid = .FALSE.
         cache%pending = .FALSE.
         cache%misses = cache%misses + 1
         RETURN
      END IF
      IF (cache%valid) THEN
         IF (SIZE(cache%integer_key) == SIZE(integers) .AND. SIZE(cache%real_key) == SIZE(reals)) THEN
            hit = ALL(cache%integer_key == integers) .AND. ALL(cache%real_key == reals)
         END IF
      END IF
      IF (.NOT. hit) THEN
         ! Invalidate first: no stale payload can be accepted after a failed direct rebuild.
         cache%valid = .FALSE.
         cache%pending = .TRUE.
         cache%integer_key = integers
         cache%real_key = reals
         cache%misses = cache%misses + 1
         RETURN
      END IF
      cache%hits = cache%hits + 1
      images%projector_translation = cache%translation
      images%projector_overlap = cache%overlap
      images%level = cache%level
      ALLOCATE (images%scaled_projector_overlap, MOLD=images%projector_overlap)
      DO image = 1, SIZE(images%projector_overlap, 3)
         DO ao = 1, SIZE(images%projector_overlap, 2)
            images%scaled_projector_overlap(:, ao, image) = images%level*images%projector_overlap(:, ao, image)
         END DO
      END DO
      auxiliary%scaled_pv_overlap = SUM(images%scaled_projector_overlap, DIM=3)
   END SUBROUTINE projector_cache_lookup

! **************************************************************************************************
!> \brief Store the provider's compact payload after a miss; no H(R) or derivative arrays are retained.
!> \param cache consumer-owned storage with a pending exact key
!> \param images successful provider result
! **************************************************************************************************
   SUBROUTINE projector_cache_store(cache, images)
      TYPE(cp2k_acp_projector_cache_type), INTENT(INOUT) :: cache
      TYPE(acp_image_type), INTENT(IN)                   :: images

      cache%valid = .FALSE.
      IF (.NOT. cache%pending) RETURN
      cache%pending = .FALSE.
      IF (.NOT. ALLOCATED(images%projector_translation)) RETURN
      IF (.NOT. ALLOCATED(images%projector_overlap)) RETURN
      IF (.NOT. ALLOCATED(images%level)) RETURN
      IF (SIZE(images%projector_translation, 1) /= 3) RETURN
      IF (SIZE(images%projector_overlap, 1) /= SIZE(images%level)) RETURN
      IF (SIZE(images%projector_overlap, 3) /= SIZE(images%projector_translation, 2)) RETURN
      IF (.NOT. ALL(ieee_is_finite(images%projector_translation))) RETURN
      IF (.NOT. ALL(ieee_is_finite(images%projector_overlap))) RETURN
      IF (.NOT. ALL(ieee_is_finite(images%level))) RETURN
      cache%translation = images%projector_translation
      cache%overlap = images%projector_overlap
      cache%level = images%level
      cache%valid = .TRUE.
   END SUBROUTINE projector_cache_store

! **************************************************************************************************
!> \brief Disable reuse for the permanent direct oracle.
!> \param cache consumer-owned storage
! **************************************************************************************************
   SUBROUTINE projector_cache_bypass(cache)
      TYPE(cp2k_acp_projector_cache_type), INTENT(INOUT) :: cache

      cache%valid = .FALSE.
      cache%pending = .FALSE.
      cache%bypasses = cache%bypasses + 1
   END SUBROUTINE projector_cache_bypass

! **************************************************************************************************
!> \brief Diagnostic counters; callers cannot mutate the private cache keys.
!> \param cache consumer-owned storage
!> \param hits exact matches
!> \param misses new or changed keys
!> \param bypasses explicit direct requests
!> \param valid reusable compact payload
! **************************************************************************************************
   SUBROUTINE projector_cache_counts(cache, hits, misses, bypasses, valid)
      TYPE(cp2k_acp_projector_cache_type), INTENT(IN)    :: cache
      INTEGER, INTENT(OUT)                               :: hits, misses, bypasses
      LOGICAL, INTENT(OUT)                               :: valid

      hits = cache%hits
      misses = cache%misses
      bypasses = cache%bypasses
      valid = cache%valid
   END SUBROUTINE projector_cache_counts
#endif
END MODULE tblite_projector_cache
