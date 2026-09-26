!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                               !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                    !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                        !
!--------------------------------------------------------------------------------------------------!
PROGRAM tblite_exchange_derivative_unittest
#if defined(__TBLITE_GXTB_ACCELERATION)
   USE mctc_env, ONLY: error_type, wp
   USE mctc_io, ONLY: new, structure_type
   USE tblite_container, ONLY: container_cache
   USE tblite_cutoff, ONLY: get_lattice_points
   USE tblite_exchange_stream, ONLY: cp2k_exchange_kmesh, cp2k_exchange_kmesh_gradient, &
                                    cp2k_exchange_stream_type, cp2k_exchange_stream_begin, &
                                    cp2k_exchange_stream_batched, cp2k_exchange_stream_end, &
                                    cp2k_exchange_stream_reverse_push, cp2k_exchange_stream_reverse_transform_advance, &
                                    cp2k_exchange_stream_reverse_batch_apply, cp2k_exchange_stream_reverse_batch_pull, &
                                    cp2k_exchange_stream_reverse_batch_advance, cp2k_exchange_stream_reverse_result, &
                                    cp2k_exchange_transform_dense, cp2k_exchange_transform_separable_dft, &
                                    cp2k_exchange_transform_mixed_radix_fft
   USE tblite_stream_plan, ONLY: cp2k_stream_plan_type
   USE tblite_wavefunction_type, ONLY: new_wavefunction, wavefunction_type
   USE tblite_wignerseitz, ONLY: get_wignerseitz_pairs, get_wignerseitz_pair_weights, wignerseitz_threshold
   USE tblite_xtb_calculator, ONLY: xtb_calculator
   USE tblite_xtb_gxtb, ONLY: new_gxtb_calculator
   IMPLICIT NONE
   TYPE(structure_type) :: mol, reference
   TYPE(xtb_calculator) :: calc
   TYPE(container_cache) :: cache
   TYPE(cp2k_stream_plan_type) :: plan
   TYPE(wavefunction_type) :: wfn
   TYPE(error_type), ALLOCATABLE :: error
   INTEGER, PARAMETER :: mesh(3) = [3, 1, 1], nk = 3
   INTEGER, PARAMETER :: backends(3) = [cp2k_exchange_transform_dense, &
                                      cp2k_exchange_transform_separable_dft, cp2k_exchange_transform_mixed_radix_fft]
   REAL(wp), PARAMETER :: steps(2) = [2.0E-5_wp, 1.0E-5_wp]
   REAL(wp) :: lattice(3, 3), xyz(3, 2), kfrac(3, nk), weights(nk), h, plus, minus
   REAL(wp) :: gradient(3, 2), sigma(3, 3), numerical_gradient(3, 2), numerical_sigma(3, 3)
   REAL(wp) :: force_error, strain_error, batch_error, deformation(3, 3), unused_energy
   REAL(wp), ALLOCATABLE :: vsh(:)
   COMPLEX(wp), ALLOCATABLE :: density(:, :, :, :), overlap(:, :, :), fock(:, :, :, :), response(:, :, :)
   INTEGER :: position, spin, nao, i, j, k, channel, axis, atom, row, col, step, backend

   lattice = RESHAPE([4.0_wp, 0.2_wp, 0.1_wp, 0.0_wp, 8.0_wp, 0.3_wp, &
                     0.0_wp, 0.0_wp, 9.0_wp], [3, 3])
   xyz(:, 1) = 0.0_wp
   ! The perpendicular offset retains the first BvK face bisector in a skew cell.
   xyz(:, 2) = 1.5_wp*lattice(:, 1) + 0.12_wp*[-0.2_wp, 4.0_wp, 0.0_wp]
   kfrac = 0.0_wp
   kfrac(1, :) = [0.0_wp, 1.0_wp/3.0_wp, 2.0_wp/3.0_wp]
   weights = 1.0_wp/REAL(nk, wp)
   CALL new(reference, [1, 1], xyz, lattice=lattice, periodic=[.TRUE., .TRUE., .TRUE.])
   mol = reference
   CALL new_gxtb_calculator(calc, mol, error, accuracy=0.1_wp)
   CALL check(.NOT. ALLOCATED(error), 'calculator')
   nao = calc%bas%nao
   ALLOCATE (overlap(nao, nao, nk), response(nao, nao, nk), vsh(calc%bas%nsh))
   DO spin = 1, 2
      CALL new_wavefunction(wfn, mol%nat, calc%bas%nsh, nao, spin, 0.0_wp)
      ALLOCATE (density(nao, nao, spin, nk), fock(nao, nao, spin, nk))
      DO k = 1, nk
         DO i = 1, nao
            DO j = 1, nao
               overlap(i, j, k) = CMPLX(0.002_wp*(i + j), 0.0002_wp*(i - j)*k, wp)
               DO channel = 1, spin
                  density(i, j, channel, k) = CMPLX(0.02_wp*(i + j + channel), 0.002_wp*(i - j)*k, wp)
               END DO
            END DO
            overlap(i, i, k) = 1.0_wp
         END DO
      END DO
      DO position = -1, 1
         mol = reference
         mol%xyz(1, 2) = mol%xyz(1, 2) + REAL(position, wp)*0.005_wp
         CALL check_ws_face()
         unused_energy = energy(mol)
         CALL cp2k_exchange_kmesh_gradient(calc, mol, cache, plan, mesh, kfrac, weights, density, overlap, &
                                           response, gradient, sigma, error)
         CALL check(.NOT. ALLOCATED(error), 'complete gradient')
         DO backend = 1, SIZE(backends)
            CALL compare_batch(backends(backend), batch_error)
            WRITE (*, '(A,3I4,ES14.6)') 'WS batch residual: ', spin, position, backend, batch_error
            CALL check(batch_error < 1.0E-11_wp, 'batched/dense force, strain and overlap adjoints')
         END DO
         ! Hold AO density/overlap and model parameters fixed. This differentiates
         ! exactly the radial/image part returned above, not a full SCF force.
         DO step = 1, SIZE(steps)
            h = steps(step)
            DO atom = 1, mol%nat
               DO axis = 1, 3
                  xyz = mol%xyz
                  mol%xyz(axis, atom) = xyz(axis, atom) + h
                  plus = energy(mol)
                  mol%xyz(axis, atom) = xyz(axis, atom) - h
                  minus = energy(mol)
                  mol%xyz = xyz
                  numerical_gradient(axis, atom) = (plus - minus)/(2.0_wp*h)
               END DO
            END DO
            DO row = 1, 3
               DO col = 1, 3
                  deformation = 0.0_wp
                  DO axis = 1, 3
                     deformation(axis, axis) = 1.0_wp
                  END DO
                  xyz = mol%xyz
                  deformation(col, row) = deformation(col, row) + h
                  mol%xyz = MATMUL(deformation, xyz)
                  mol%lattice = MATMUL(deformation, lattice)
                  plus = energy(mol)
                  deformation(col, row) = deformation(col, row) - 2.0_wp*h
                  mol%xyz = MATMUL(deformation, xyz)
                  mol%lattice = MATMUL(deformation, lattice)
                  minus = energy(mol)
                  mol%xyz = xyz
                  mol%lattice = lattice
                  numerical_sigma(row, col) = (plus - minus)/(2.0_wp*h)
               END DO
            END DO
            force_error = MAXVAL(ABS(numerical_gradient - gradient))
            strain_error = MAXVAL(ABS(numerical_sigma - sigma))
            WRITE (*, '(A,3I4,2ES14.6)') 'WS finite-difference residuals: ', spin, position, step, &
               force_error, strain_error
            CALL check(force_error < 1.0E-7_wp, 'independent coordinate differences')
            CALL check(strain_error < 1.0E-7_wp, 'independent homogeneous-strain differences')
         END DO
      END DO
      DEALLOCATE (density, fock)
   END DO
   WRITE (*, '(A)') 'CP2K STREAM DERIVATIVES: RKS/UKS, smooth WS faces, batch adjoints and finite differences passed'
CONTAINS
   SUBROUTINE check_ws_face()
      REAL(wp) :: supercell(3, 3)
      REAL(wp), ALLOCATABLE :: translations(:, :), ws_weights(:), coordinate_response(:, :), strain_response(:, :, :)
      INTEGER, ALLOCATABLE :: indices(:)
      INTEGER :: nimages, direction

      supercell = mol%lattice
      DO direction = 1, 3
         supercell(:, direction) = supercell(:, direction)*REAL(mesh(direction), wp)
      END DO
      CALL get_lattice_points(mol%periodic, supercell, wignerseitz_threshold, translations)
      ALLOCATE (indices(SIZE(translations, 2)), ws_weights(SIZE(translations, 2)), &
                coordinate_response(3, SIZE(translations, 2)), strain_response(3, 3, SIZE(translations, 2)))
      CALL get_wignerseitz_pairs(nimages, translations, mol%xyz(:, 1) - mol%xyz(:, 2), indices)
      CALL check(nimages == 2, 'fixture really straddles a two-image WS face')
      CALL get_wignerseitz_pair_weights(translations, indices(:nimages), mol%xyz(:, 1) - mol%xyz(:, 2), &
                                       ws_weights, coordinate_response, strain_response)
      IF (position /= 0) THEN
         CALL check(MAXVAL(ABS(ws_weights(:nimages) - 0.5_wp)) > 1.0E-6_wp, 'unequal smooth image weights')
         CALL check(MAXVAL(ABS(coordinate_response(:, :nimages))) > 1.0E-6_wp, 'nonzero image-weight response')
      END IF
   END SUBROUTINE check_ws_face

   FUNCTION energy(structure) RESULT(value)
      TYPE(structure_type), INTENT(IN) :: structure
      REAL(wp) :: value

      CALL cp2k_exchange_kmesh(calc, structure, cache, plan, wfn, mesh, kfrac, weights, density, overlap, &
                               fock, vsh, value, error)
      CALL check(.NOT. ALLOCATED(error), 'displaced exchange energy')
   END FUNCTION energy

   SUBROUTINE compare_batch(transform, residual)
      INTEGER, INTENT(IN) :: transform
      REAL(wp), INTENT(OUT) :: residual
      TYPE(cp2k_exchange_stream_type) :: stream
      COMPLEX(wp) :: batched_response(nao, nao, nk)
      REAL(wp) :: batched_gradient(3, 2), batched_sigma(3, 3)
      INTEGER :: tile, tiles, ik, pairs
      LOGICAL :: done, tile_done

      pairs = MIN(2, nao*nao)
      tiles = (nao*nao + pairs - 1)/pairs
      IF (transform == cp2k_exchange_transform_dense) tiles = 1
      batched_response = 0.0_wp
      CALL cp2k_exchange_stream_begin(stream, calc, mol, cache, plan, wfn, mesh, kfrac, weights, error, &
                                       mode=cp2k_exchange_stream_batched, batch_size=1, reverse_only=.TRUE., &
                                       transform_backend=transform, transform_pair_batch_size=pairs)
      CALL check(.NOT. ALLOCATED(error), 'reverse batch begin')
      done = .FALSE.
      DO WHILE (.NOT. done)
         DO tile = 1, tiles
            DO ik = nk, 1, -1
               CALL cp2k_exchange_stream_reverse_push(stream, ik, density(:, :, :, ik), overlap(:, :, ik), error)
               CALL check(.NOT. ALLOCATED(error), 'reverse batch push')
            END DO
            CALL cp2k_exchange_stream_reverse_transform_advance(stream, tile_done, error)
            CALL check(.NOT. ALLOCATED(error), 'reverse transform')
            CALL check(tile_done .EQV. (tile == tiles), 'reverse tile count')
         END DO
         CALL cp2k_exchange_stream_reverse_batch_apply(stream, calc, mol, cache, plan, error)
         CALL check(.NOT. ALLOCATED(error), 'reverse batch apply')
         DO ik = 1, nk
            CALL cp2k_exchange_stream_reverse_batch_pull(stream, ik, density(:, :, :, ik), overlap(:, :, ik), &
                                                        batched_response(:, :, ik), error)
            CALL check(.NOT. ALLOCATED(error), 'reverse batch pull')
         END DO
         CALL cp2k_exchange_stream_reverse_batch_advance(stream, done, error)
         CALL check(.NOT. ALLOCATED(error), 'reverse batch advance')
      END DO
      CALL cp2k_exchange_stream_reverse_result(stream, calc, mol, cache, plan, batched_gradient, batched_sigma, error)
      CALL check(.NOT. ALLOCATED(error), 'reverse batch result')
      residual = MAX(MAXVAL(ABS(batched_gradient - gradient)), MAXVAL(ABS(batched_sigma - sigma)), &
                     MAXVAL(ABS(batched_response - response)))
      CALL cp2k_exchange_stream_end(stream, error)
      CALL check(.NOT. ALLOCATED(error), 'reverse batch end')
   END SUBROUTINE compare_batch

   SUBROUTINE check(condition, message)
      LOGICAL, INTENT(IN) :: condition
      CHARACTER(LEN=*), INTENT(IN) :: message

      IF (.NOT. condition) THEN
         WRITE (*, '(A)') 'FAIL: '//message
         IF (ALLOCATED(error)) WRITE (*, '(A)') error%message
         STOP 1
      END IF
   END SUBROUTINE check
#else
   WRITE (*, '(A)') 'CP2K STREAM DERIVATIVES: skipped (physical response interface unavailable)'
#endif
END PROGRAM tblite_exchange_derivative_unittest
