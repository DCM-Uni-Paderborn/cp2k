!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
PROGRAM tblite_exchange_stream_unittest
#if defined(__TBLITE_GXTB_ACCELERATION)
   USE tblite_stream_plan, ONLY: cp2k_stream_plan_type
   USE, INTRINSIC :: iso_fortran_env, ONLY: int64
   USE mctc_env, ONLY: error_type, wp
   USE mctc_io, ONLY: new, structure_type
   USE tblite_container, ONLY: container_cache
   USE tblite_exchange_stream, ONLY: cp2k_exchange_stream_type, cp2k_exchange_partial_begin, &
                                    cp2k_exchange_partial_push, cp2k_exchange_partial_reverse_push, cp2k_exchange_partial_discard, &
                                     cp2k_exchange_partial_descriptor_type, cp2k_exchange_partial_forward_type, &
                                     cp2k_exchange_partial_reverse_type, cp2k_exchange_partial_export_forward, &
                                     cp2k_exchange_partial_export_reverse, &
                                     cp2k_exchange_transform_dense, cp2k_exchange_transform_separable_dft, &
                                     cp2k_exchange_transform_mixed_radix_fft, cp2k_exchange_stream_pair_tile_size, &
                                     cp2k_exchange_stream_image_batch_size
   USE tblite_wavefunction_type, ONLY: wavefunction_type, new_wavefunction
   USE tblite_xtb_calculator, ONLY: xtb_calculator
   USE tblite_xtb_gxtb, ONLY: new_gxtb_calculator
   IMPLICIT NONE
   TYPE(structure_type) :: mol
   TYPE(xtb_calculator) :: calc
   TYPE(container_cache) :: cache
   TYPE(cp2k_stream_plan_type) :: cache_plan
   TYPE(wavefunction_type) :: wfn
   TYPE(error_type), ALLOCATABLE :: error
   TYPE(cp2k_exchange_stream_type) :: stream
   TYPE(cp2k_exchange_partial_descriptor_type) :: descriptor
   TYPE(cp2k_exchange_partial_forward_type) :: forward_buffer
   TYPE(cp2k_exchange_partial_reverse_type) :: reverse_buffer
   INTEGER, ALLOCATABLE :: ownership(:)
   REAL(wp) :: lattice(3, 3), xyz(3, 3), kfrac(3, 3), weights(3)
   COMPLEX(wp), ALLOCATABLE :: density(:, :, :), overlap(:, :)
   INTEGER :: backend, direction, spin, a, b, nao, trial, ik, tile_size, nfields
   INTEGER, PARAMETER :: backends(3) = [cp2k_exchange_transform_dense, &
                                        cp2k_exchange_transform_separable_dft, cp2k_exchange_transform_mixed_radix_fft]
   INTEGER, PARAMETER :: invalid_indices(4) = [-1, 0, 4, HUGE(1)]
   INTEGER(int64) :: transaction(2)
   INTEGER(int64) :: bytes_per_pair

   CALL check(cp2k_exchange_stream_pair_tile_size(13, 216, 1, 3) == 168, 'NaCl forward tile')
   CALL check(cp2k_exchange_stream_pair_tile_size(13, 343, 2, 6) == 168, 'NaCl reverse UKS tile')
   CALL check(cp2k_exchange_stream_pair_tile_size(1, 343, 1, 3) == 1, 'scalar basis tile')
   CALL check(cp2k_exchange_stream_pair_tile_size(10000, 1000000, 2, 6) == 1, 'single pair lower bound')
   DO spin = 1, 2
      DO nfields = 3, 6, 3
         tile_size = cp2k_exchange_stream_pair_tile_size(10000, 512, spin, nfields)
         bytes_per_pair = 16_int64*512_int64*(INT(spin*nfields, int64) + 3_int64)
         CALL check(INT(tile_size, int64)*bytes_per_pair <= 64_int64*1024_int64**2, 'tile memory budget')
         CALL check(INT(tile_size + 1, int64)*bytes_per_pair > 64_int64*1024_int64**2, 'largest bounded tile')
         CALL check(tile_size < 10000**2, 'no complete AO-pair mesh')
      END DO
   END DO
   WRITE (*, '(A)') 'CP2K STREAM TILES: forward/reverse, scalar/NaCl/large basis, 64 MiB bounds passed'

   CALL check(cp2k_exchange_stream_image_batch_size(13, 216, 1) == 215, 'NaCl image batch')
   CALL check(cp2k_exchange_stream_image_batch_size(13, 343, 2) == 342, 'NaCl UKS image batch')
   CALL check(cp2k_exchange_stream_image_batch_size(1, 1, 1) == 1, 'Gamma image minimum')
   CALL check(cp2k_exchange_stream_image_batch_size(10000, 512, 2) == 1, 'large basis image minimum')
   DO spin = 1, 2
      tile_size = cp2k_exchange_stream_image_batch_size(256, 512, spin)
      bytes_per_pair = 16_int64*6_int64*256_int64**2*INT(spin, int64)
      CALL check(INT(tile_size, int64)*bytes_per_pair <= 64_int64*1024_int64**2, 'image memory budget')
      CALL check(INT(tile_size + 1, int64)*bytes_per_pair > 64_int64*1024_int64**2, 'largest bounded image batch')
      CALL check(tile_size < 512, 'no complete image mesh')
   END DO
   WRITE (*, '(A)') 'CP2K IMAGE BATCHES: forward/reverse, 64 MiB bounds passed'

   lattice = 0.0_wp
   lattice(1, 1) = 15.0_wp
   lattice(2, 2) = 16.0_wp
   lattice(3, 3) = 17.0_wp
   xyz = RESHAPE([0.0_wp, 0.0_wp, 0.0_wp, 1.8_wp, 0.0_wp, 0.0_wp, &
                  -0.5_wp, 1.7_wp, 0.1_wp], [3, 3])
   CALL new(mol, [8, 1, 1], xyz, lattice=lattice, periodic=[.TRUE., .TRUE., .TRUE.])
   CALL new_gxtb_calculator(calc, mol, error, accuracy=0.1_wp)
   CALL check(.NOT. ALLOCATED(error), 'calculator')
   nao = calc%bas%nao
   kfrac = 0.0_wp
   kfrac(1, :) = [0.0_wp, 1.0_wp/3.0_wp, 2.0_wp/3.0_wp]
   weights = 1.0_wp/3.0_wp
   transaction = [1_int64, 0_int64]
   DO spin = 1, 2
      CALL new_wavefunction(wfn, mol%nat, calc%bas%nsh, nao, spin, 0.0_wp)
      ALLOCATE (density(nao, nao, spin), overlap(nao, nao))
      DO a = 1, nao
         DO b = 1, nao
            overlap(a, b) = CMPLX(0.001_wp*(a + b), 0.0001_wp*(a - b), wp)
            density(a, b, :) = CMPLX(0.01_wp*(a + b), 0.001_wp*(a - b), wp)
         END DO
         overlap(a, a) = 1.0_wp
      END DO
      DO backend = 1, SIZE(backends)
         DO direction = 1, 2
            transaction(2) = transaction(2) + 1_int64
            CALL cp2k_exchange_partial_begin(stream, calc, mol, cache, cache_plan, wfn, [3, 1, 1], &
                                             kfrac, weights, 2, transaction, error, reverse_only=direction == 2, &
                                             transform_backend=backends(backend), transform_pair_batch_size=2)
            CALL check(.NOT. ALLOCATED(error), 'partial transaction')
            DO trial = 1, SIZE(invalid_indices)
               CALL push(invalid_indices(trial))
               CALL check(ALLOCATED(error), 'invalid index rejected before ownership access')
               CALL check(INDEX(error%message, 'index') > 0, 'specific index rejection')
               DEALLOCATE (error)
            END DO
            ! Rejection must not poison ownership or prevent subsequent valid input.
            DO ik = 3, 1, -1
               CALL push(ik)
               CALL check(.NOT. ALLOCATED(error), 'valid input after rejected index')
            END DO
            CALL push(1)
            CALL check(ALLOCATED(error), 'duplicate ownership rejected')
            DEALLOCATE (error)
            IF (direction == 1) THEN
               CALL cp2k_exchange_partial_export_forward(stream, descriptor, ownership, forward_buffer, error)
            ELSE
               CALL cp2k_exchange_partial_export_reverse(stream, descriptor, ownership, reverse_buffer, error)
            END IF
            CALL check(.NOT. ALLOCATED(error), 'export before contributor discard')
            CALL check(ALL(ownership == 1), 'invalid input leaves ownership intact')
            CALL cp2k_exchange_partial_discard(stream, error)
            CALL check(.NOT. ALLOCATED(error), 'contributor discard')
         END DO
      END DO
      DEALLOCATE (density, overlap)
   END DO
   WRITE (*, '(A)') 'CP2K STREAM BOUNDS: RKS/UKS, three transforms, forward/reverse, recovery passed'
CONTAINS
   SUBROUTINE push(index)
      INTEGER, INTENT(IN)                                :: index

      IF (direction == 1) THEN
         CALL cp2k_exchange_partial_push(stream, index, density, overlap, error)
      ELSE
         CALL cp2k_exchange_partial_reverse_push(stream, index, density, overlap, error)
      END IF
   END SUBROUTINE push

   SUBROUTINE check(condition, message)
      LOGICAL, INTENT(IN)                                :: condition
      CHARACTER(LEN=*), INTENT(IN)                       :: message

      IF (.NOT. condition) THEN
         WRITE (*, '(A)') 'FAIL: '//message
         IF (ALLOCATED(error)) WRITE (*, '(A)') error%message
         STOP 1
      END IF
   END SUBROUTINE check
#else
   WRITE (*, '(A)') 'CP2K STREAM BOUNDS: skipped (extended mathematical kernels unavailable)'
#endif
END PROGRAM tblite_exchange_stream_unittest
