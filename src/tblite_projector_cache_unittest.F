!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
PROGRAM tblite_projector_cache_unittest
#if defined(__TBLITE_GXTB)
   USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_quiet_nan, ieee_value
   USE kinds, ONLY: dp
   USE mctc_env, ONLY: error_type
   USE mctc_io, ONLY: new, structure_type
   USE tblite_acp_cache, ONLY: acp_cache
   USE tblite_acp_images, ONLY: acp_image_type
   USE tblite_basis_cache, ONLY: basis_cache
   USE tblite_model_bridge, ONLY: cp2k_build_acp_projectors, cp2k_update_basis
   USE tblite_projector_cache, ONLY: cp2k_acp_projector_cache_type, &
                                     projector_cache_lookup, projector_cache_store, &
                                     projector_cache_bypass, projector_cache_counts
   USE tblite_wavefunction_type, ONLY: wavefunction_type
   USE tblite_xtb_calculator, ONLY: xtb_calculator
   USE tblite_xtb_gxtb, ONLY: new_gxtb_calculator
   IMPLICIT NONE
   TYPE(structure_type) :: mol
   TYPE(xtb_calculator) :: calc
   TYPE(error_type), ALLOCATABLE :: error
   TYPE(basis_cache) :: basis
   TYPE(acp_cache) :: auxiliary, oracle_auxiliary
   TYPE(wavefunction_type) :: wfn_aux
   TYPE(cp2k_acp_projector_cache_type) :: cache
   TYPE(acp_image_type) :: images, oracle
   INTEGER :: hits, misses, bypasses, saved_int
   LOGICAL :: hit, valid
   REAL(KIND=dp) :: saved, lattice(3, 3), translations(3, 3), xyz(3, 3), column(3)

   lattice = 0.0_dp
   lattice(1, 1) = 15.0_dp
   lattice(2, 2) = 16.0_dp
   lattice(3, 3) = 17.0_dp
   xyz = RESHAPE([0.0_dp, 0.0_dp, 0.0_dp, 1.8_dp, 0.0_dp, 0.0_dp, &
                  -0.5_dp, 1.7_dp, 0.1_dp], [3, 3])
   CALL new(mol, [8, 1, 1], xyz, lattice=lattice, periodic=[.TRUE., .TRUE., .TRUE.])
   CALL new_gxtb_calculator(calc, mol, error, accuracy=0.1_dp)
   IF (ALLOCATED(error)) STOP 1
   CALL cp2k_update_basis(calc, mol, basis, error, grad=.TRUE., wfn_aux=wfn_aux)
   IF (ALLOCATED(error)) STOP 2
   translations(:, 1) = -lattice(:, 1)
   translations(:, 2) = 0.0_dp
   translations(:, 3) = lattice(:, 1)

   CALL prime()
   CALL projector_cache_counts(cache, hits, misses, bypasses, valid)
   CALL check(hits == 1 .AND. misses == 1 .AND. bypasses == 1 .AND. valid, 'initial counters')
   auxiliary%scaled_pv_overlap = -999.0_dp
   basis%cgto(1, 1)%dqeffdcn = 123.0_dp
   CALL lookup()
   CALL check(hit, 'repeat hit')
   CALL check(basis%cgto(1, 1)%dqeffdcn == 123.0_dp, 'live derivative cache is not restored from payload')
   CALL check(ALL(auxiliary%scaled_pv_overlap == SUM(oracle%scaled_projector_overlap, DIM=3)), &
              'scaled derivative intermediate restored')
   CALL compare()

   saved = mol%xyz(1, 2)
   mol%xyz(1, 2) = saved + 1.0E-13_dp
   CALL miss('tiny coordinate change')
   mol%xyz(1, 2) = saved
   CALL prime()
   saved = mol%lattice(1, 1)
   mol%lattice(1, 1) = saved + 1.0E-13_dp
   CALL miss('tiny lattice change')
   mol%lattice(1, 1) = saved
   CALL prime()
   mol%periodic(3) = .FALSE.
   CALL miss('partial periodicity')
   mol%periodic(3) = .TRUE.
   CALL prime()
   column = translations(:, 1)
   translations(:, 1) = translations(:, 3)
   translations(:, 3) = column
   CALL miss('translation order')
   column = translations(:, 1)
   translations(:, 1) = translations(:, 3)
   translations(:, 3) = column
   CALL prime()
   CALL projector_cache_lookup(cache, mol, calc%bas, basis, calc%acp, auxiliary, &
                               translations(:, 2:2), images, hit)
   CALL check(.NOT. hit, 'translation count')
   CALL prime()
   saved_int = mol%num(1)
   mol%num(1) = saved_int + 1
   CALL miss('species metadata')
   mol%num(1) = saved_int
   CALL prime()
   saved = calc%bas%intcut
   calc%bas%intcut = saved + 1.0E-12_dp
   CALL miss('valence cutoff')
   calc%bas%intcut = saved
   CALL prime()
   saved = calc%acp%auxbas%intcut
   calc%acp%auxbas%intcut = saved + 1.0E-12_dp
   CALL miss('auxiliary cutoff')
   calc%acp%auxbas%intcut = saved
   CALL prime()
   saved = calc%bas%cgto(1, 1)%raw%alpha(1)
   calc%bas%cgto(1, 1)%raw%alpha(1) = saved + 1.0E-12_dp
   CALL miss('valence primitive')
   calc%bas%cgto(1, 1)%raw%alpha(1) = saved
   CALL prime()
   saved = calc%acp%auxbas%cgto(1, 1)%raw%alpha(1)
   calc%acp%auxbas%cgto(1, 1)%raw%alpha(1) = saved + 1.0E-12_dp
   CALL miss('auxiliary primitive')
   calc%acp%auxbas%cgto(1, 1)%raw%alpha(1) = saved
   CALL prime()
   saved = basis%cgto(1, 1)%qeff
   basis%cgto(1, 1)%qeff = saved + 0.001_dp
   CALL miss('charge-dependent effective coefficients')
   basis%cgto(1, 1)%qeff = saved
   CALL prime()
   saved = basis%cgto(1, 1)%norm
   basis%cgto(1, 1)%norm = saved + 0.001_dp
   CALL miss('effective normalization')
   basis%cgto(1, 1)%norm = saved
   CALL prime()
   saved = calc%acp%levels(1, 1)
   calc%acp%levels(1, 1) = saved + 1.0E-12_dp
   CALL miss('projector level')
   calc%acp%levels(1, 1) = saved
   CALL prime()
   saved_int = calc%bas%iao_sh(1)
   calc%bas%iao_sh(1) = saved_int + 1
   CALL miss('AO layout')
   calc%bas%iao_sh(1) = saved_int
   CALL prime()
   saved = translations(1, 1)
   translations(1, 1) = ieee_value(0.0_dp, ieee_quiet_nan)
   CALL miss('nonfinite key')
   CALL projector_cache_store(cache, oracle)
   CALL projector_cache_counts(cache, hits, misses, bypasses, valid)
   CALL check(.NOT. valid, 'nonfinite key cannot be certified by storing an old payload')
   translations(1, 1) = saved
   CALL prime()
   CALL projector_cache_bypass(cache)
   CALL projector_cache_counts(cache, hits, misses, bypasses, valid)
   CALL check(.NOT. valid, 'direct bypass invalidates')
   CALL lookup()
   CALL check(.NOT. hit, 'post-bypass lookup must miss')
   WRITE (*, '(A)') 'CP2K ACP PROJECTOR CACHE: all exact-key and independent direct-builder checks passed'

CONTAINS
   SUBROUTINE check(condition, message)
      LOGICAL, INTENT(IN)                                :: condition
      CHARACTER(LEN=*), INTENT(IN)                       :: message

      IF (.NOT. condition) THEN
         WRITE (*, '(A,A)') 'FAILED: ', message
         STOP 3
      END IF
   END SUBROUTINE check

   SUBROUTINE lookup()
      CALL projector_cache_lookup(cache, mol, calc%bas, basis, calc%acp, auxiliary, translations, images, hit)
   END SUBROUTINE lookup

   SUBROUTINE miss(message)
      CHARACTER(LEN=*), INTENT(IN)                       :: message

      CALL lookup()
      CALL check(.NOT. hit, message)
      CALL projector_cache_counts(cache, hits, misses, bypasses, valid)
      CALL check(.NOT. valid, 'a miss invalidates old payload')
   END SUBROUTINE miss

   SUBROUTINE prime()
      CALL projector_cache_bypass(cache)
      CALL lookup()
      CALL check(.NOT. hit, 'fresh key must miss')
      ! The model bridge builds directly; cache ownership is exclusively in CP2K.
      ! A separate live auxiliary cache keeps this oracle independent of reuse.
      CALL cp2k_build_acp_projectors(calc, mol, basis, oracle_auxiliary, translations, oracle)
      CALL projector_cache_store(cache, oracle)
      CALL lookup()
      CALL check(hit, 'stored exact key must hit')
      CALL compare()
   END SUBROUTINE prime

   SUBROUTINE compare()
      CALL check(ALL(images%projector_translation == oracle%projector_translation), 'translations')
      CALL check(ALL(images%projector_overlap == oracle%projector_overlap), 'overlaps')
      CALL check(ALL(images%scaled_projector_overlap == oracle%scaled_projector_overlap), 'scaled overlaps')
      CALL check(ALL(images%level == oracle%level), 'AO projector levels')
      CALL check(.NOT. ALLOCATED(images%hamiltonian), 'compact-only cache')
   END SUBROUTINE compare
#else
   IMPLICIT NONE
   WRITE (*, '(A)') 'ACP cache test requires a g-xTB-enabled tblite provider'
#endif
END PROGRAM tblite_projector_cache_unittest
