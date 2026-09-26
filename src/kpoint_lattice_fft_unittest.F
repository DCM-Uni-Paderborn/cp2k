!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
PROGRAM kpoint_lattice_fft_unittest
   USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_quiet_nan, ieee_value
   USE fft_lib, ONLY: fft_supports_arbitrary_lengths
   USE fft_tools, ONLY: finalize_fft, &
                        init_fft
   USE kinds, ONLY: dp
   USE kpoint_lattice_fft, ONLY: cell_to_k_grid_fft, &
                                 certify_kpoint_mesh, kpoint_mesh_certificate_type, &
                                 k_grid_to_cell_execute, &
                                 k_grid_to_cell_fft, &
                                 k_grid_to_cell_prepare, &
                                 k_grid_to_cell_release, &
                                 k_grid_to_cell_work_type, &
                                 regular_kpoint_grid
   USE lattice_preconditioner_operator, ONLY: &
      apply_lattice_inverse_dense, apply_lattice_state_inverse_dense, build_lattice_inverse, &
      build_lattice_local_correction, lattice_fft_auto_select, lattice_fft_cost_ratio, &
      lattice_fft_storage_ratio, lattice_grid, lattice_local_cost_ratio, &
      lattice_local_storage_ratio
   USE mathconstants, ONLY: gaussi, &
                            twopi, &
                            z_zero
   USE message_passing, ONLY: mp_comm_type, &
                              mp_world_finalize, &
                              mp_world_init

   IMPLICIT NONE

   INTEGER, PARAMETER                                  :: base_block = 2, base_cell = 8, base_n = 16, &
                                                          base_ncol = 3, ncell = 31, nkp = 64, &
                                                          nround = 16, nselected = 3
   INTEGER, PARAMETER                                  :: large_block = 1000, large_ncol = 1000, &
                                                          large_num_pe = 16
   COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)   :: inverse_k
   COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :) :: state_inverse_k
   COMPLEX(KIND=dp), DIMENSION(2, 3, nkp)              :: direct, transformed
   COMPLEX(KIND=dp), DIMENSION(2, 3, nround)           :: reconstructed, round_k
   COMPLEX(KIND=dp), DIMENSION(2, 3, nselected)         :: selected_transformed
   INTEGER                                             :: d, i, icell, info, j, jcell, ik, ix, iy, iz
   INTEGER, ALLOCATABLE, DIMENSION(:)                  :: correction_indices
   INTEGER, DIMENSION(3)                               :: base_dims, delta, large_dims
   INTEGER, DIMENSION(3, base_cell)                    :: base_index
   INTEGER, DIMENSION(3)                               :: nkp_grid
   INTEGER, DIMENSION(3, ncell)                        :: index_to_cell
   INTEGER, DIMENSION(3, nround)                       :: round_index
   INTEGER, DIMENSION(nselected)                       :: selected_kpoints
   LOGICAL                                             :: used_fft
   CHARACTER(LEN=16)                                   :: fft_backend
   REAL(KIND=dp), DIMENSION(3, 3)                      :: hmat
   REAL(KIND=dp)                                       :: base_residual_norm, corrected_projection_error, &
                                                          corrected_residual_norm, cost_ratio, &
                                                          expected_cost_ratio, imaginary_residual, &
                                                          large_block_real, large_ncell_real, &
                                                          large_ncol_real, projection_error
   REAL(KIND=dp), DIMENSION(3, base_cell)              :: base_xkp
   REAL(KIND=dp), DIMENSION(base_n, base_n)            :: base_overlap
   REAL(KIND=dp), DIMENSION(base_n, base_n)            :: defect_operator
   REAL(KIND=dp), DIMENSION(base_n, base_ncol)          :: base_input, base_input_rotated, &
                                                           base_output, base_output_rotated, &
                                                           base_state_output, defect_residual
   REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)         :: coarse_inverse, defect_base_output, &
                                                          defect_corrected_output, local_operator, &
                                                          selected_output, selected_rhs
   REAL(KIND=dp), DIMENSION(base_ncol, base_ncol)       :: column_rotation
   REAL(KIND=dp), DIMENSION(3, nkp)                    :: xkp, xkp_bad
   REAL(KIND=dp), DIMENSION(2, 3, ncell)               :: values_rs
   REAL(KIND=dp), DIMENSION(2, 3, nround)              :: round_rs
   TYPE(mp_comm_type)                                  :: mp_comm
   TYPE(k_grid_to_cell_work_type)                      :: inverse_work
   REAL(KIND=dp), DIMENSION(9)                         :: subset_weights
   COMPLEX(KIND=dp), DIMENSION(2, 3, ncell)            :: weighted_reference, weighted_result

   CALL mp_world_init(mp_comm)
   fft_backend = "FFTSG"
   IF (COMMAND_ARGUMENT_COUNT() > 0) CALL GET_COMMAND_ARGUMENT(1, fft_backend)
   CALL init_fft(TRIM(fft_backend), alltoall=.FALSE., pool_limit=4, &
                 wisdom_file="", plan_style=1)
   CALL check_mesh_certificate()

   nkp_grid = [4, 4, 4]
   hmat = RESHAPE([4.0_dp, 0.1_dp, 0.2_dp, &
                   0.0_dp, 5.0_dp, 0.3_dp, &
                   0.0_dp, 0.0_dp, 6.0_dp], [3, 3])

   ik = 0
   DO iz = 0, nkp_grid(3) - 1
      DO iy = 0, nkp_grid(2) - 1
         DO ix = 0, nkp_grid(1) - 1
            ik = ik + 1
            xkp(:, ik) = [(REAL(ix, dp) + 0.5_dp)/REAL(nkp_grid(1), dp), &
                          REAL(iy, dp)/REAL(nkp_grid(2), dp), &
                          (REAL(iz, dp) + 0.25_dp)/REAL(nkp_grid(3), dp)]
         END DO
      END DO
   END DO
   IF (.NOT. regular_kpoint_grid(xkp, nkp_grid)) ERROR STOP "Regular k-point grid not recognized"

   DO icell = 1, ncell
      index_to_cell(:, icell) = [MODULO(3*icell + 1, 11) - 5, &
                                 MODULO(5*icell + 2, 9) - 4, &
                                 MODULO(7*icell + 3, 7) - 3]
      DO d = 1, 3
         values_rs(1, d, icell) = SIN(0.13_dp*REAL(icell + 2*d, dp))
         values_rs(2, d, icell) = COS(0.17_dp*REAL(2*icell - d, dp))
      END DO
   END DO

   CALL cell_to_k_grid_fft(values_rs, index_to_cell, xkp, nkp_grid, transformed, used_fft)
   IF (.NOT. used_fft) ERROR STOP "Complete regular grid did not use lattice FFT"
   CALL direct_reference(values_rs, index_to_cell, xkp, direct)
   IF (MAXVAL(ABS(transformed - direct)) > 2.0E-12_dp) THEN
      ERROR STOP "Shifted lattice FFT differs from direct transform"
   END IF

   selected_kpoints = [2, 31, 64]
   CALL cell_to_k_grid_fft(values_rs, index_to_cell, xkp, nkp_grid, selected_transformed, &
                           used_fft, selected_kpoints=selected_kpoints)
   IF (.NOT. used_fft) ERROR STOP "Selected distributed output did not use lattice FFT"
   IF (MAXVAL(ABS(selected_transformed - direct(:, :, selected_kpoints))) > 2.0E-12_dp) THEN
      ERROR STOP "Selected lattice-FFT output differs from the full transform"
   END IF

   ! A reduced list is opt-in: no symmetry weights enter a forward transform.
   CALL cell_to_k_grid_fft(values_rs, index_to_cell, xkp(:, selected_kpoints), nkp_grid, &
                           selected_transformed, used_fft)
   IF (used_fft) ERROR STOP "Incomplete forward grid was accepted without opt-in"
   CALL cell_to_k_grid_fft(values_rs, index_to_cell, xkp(:, selected_kpoints), nkp_grid, &
                           selected_transformed, used_fft, allow_incomplete=.TRUE.)
   IF (.NOT. used_fft) ERROR STOP "Regular reduced forward grid did not use FFT"
   IF (MAXVAL(ABS(selected_transformed - direct(:, :, selected_kpoints))) > 2.0E-12_dp) THEN
      ERROR STOP "Reduced forward transform differs from the full-grid phase sum"
   END IF
   IF (regular_kpoint_grid(xkp(:, [2, 31, 2]), nkp_grid, allow_incomplete=.TRUE.)) THEN
      ERROR STOP "Duplicate reduced points were accepted"
   END IF

   CALL cell_to_k_grid_fft(values_rs, index_to_cell, xkp, nkp_grid, transformed, used_fft, 2, hmat)
   IF (.NOT. used_fft) ERROR STOP "Derivative did not use lattice FFT"
   CALL direct_reference(values_rs, index_to_cell, xkp, direct, 2, hmat)
   IF (MAXVAL(ABS(transformed - direct)) > 5.0E-12_dp) THEN
      ERROR STOP "Lattice-FFT derivative differs from direct transform"
   END IF
   CALL cell_to_k_grid_fft(values_rs, index_to_cell, xkp(:, selected_kpoints), nkp_grid, &
                           selected_transformed, used_fft, 2, hmat, allow_incomplete=.TRUE.)
   IF (.NOT. used_fft) ERROR STOP "Reduced forward derivative did not use FFT"
   IF (MAXVAL(ABS(selected_transformed - direct(:, :, selected_kpoints))) > 5.0E-12_dp) THEN
      ERROR STOP "Reduced forward derivative differs from direct transform"
   END IF

   CALL check_prime_grid()

   nkp_grid = [4, 4, 1]
   ik = 0
   DO iy = 0, nkp_grid(2) - 1
      DO ix = 0, nkp_grid(1) - 1
         ik = ik + 1
         xkp(:, ik) = [(REAL(ix, dp) + 0.5_dp)/REAL(nkp_grid(1), dp), &
                       REAL(iy, dp)/REAL(nkp_grid(2), dp), 0.0_dp]
      END DO
   END DO
   CALL cell_to_k_grid_fft(values_rs, index_to_cell, xkp(:, 1:ik), nkp_grid, &
                           transformed(:, :, 1:ik), used_fft)
   IF (.NOT. used_fft) ERROR STOP "Singleton grid direction did not use padded lattice FFT"
   CALL direct_reference(values_rs, index_to_cell, xkp(:, 1:ik), direct(:, :, 1:ik))
   IF (MAXVAL(ABS(transformed(:, :, 1:ik) - direct(:, :, 1:ik))) > 2.0E-12_dp) THEN
      ERROR STOP "Padded lattice FFT differs from direct transform"
   END IF

   icell = 0
   DO iy = 0, nkp_grid(2) - 1
      DO ix = 0, nkp_grid(1) - 1
         icell = icell + 1
         round_index(:, icell) = [ix, iy, 0]
         round_rs(:, :, icell) = values_rs(:, :, icell)
      END DO
   END DO
   CALL cell_to_k_grid_fft(round_rs, round_index, xkp(:, 1:ik), nkp_grid, round_k, used_fft)
   IF (.NOT. used_fft) ERROR STOP "Round-trip forward transform did not use lattice FFT"
   CALL k_grid_to_cell_fft(round_k, xkp(:, 1:ik), nkp_grid, round_index, reconstructed, used_fft)
   IF (.NOT. used_fft) ERROR STOP "Round-trip inverse transform did not use lattice FFT"
   IF (MAXVAL(ABS(reconstructed - CMPLX(round_rs, KIND=dp))) > 3.0E-12_dp) THEN
      ERROR STOP "Forward/inverse lattice-FFT round trip is inconsistent"
   END IF

   ! An incomplete shifted grid uses the supplied weights exactly, including
   ! negative translations. Reuse the preparation for separate matrix batches.
   subset_weights = [(REAL(ik, dp)/45.0_dp, ik=1, 9)]
   weighted_reference = z_zero
   DO icell = 1, ncell
      DO ik = 1, 9
         weighted_reference(:, :, icell) = weighted_reference(:, :, icell) + subset_weights(ik)* &
                                           EXP(-gaussi*twopi*SUM(xkp(:, ik)*REAL(index_to_cell(:, icell), dp)))*round_k(:, :, ik)
      END DO
   END DO
   CALL k_grid_to_cell_prepare(inverse_work, xkp(:, 1:9), nkp_grid, index_to_cell, &
                               weights=subset_weights, allow_incomplete=.TRUE.)
   DO i = 1, 2
      CALL k_grid_to_cell_execute(inverse_work, round_k(i:i, :, 1:9), weighted_result(i:i, :, :), used_fft)
      IF (.NOT. used_fft) ERROR STOP "Prepared incomplete grid did not use FFT"
   END DO
   IF (MAXVAL(ABS(weighted_result - weighted_reference)) > 3.0E-12_dp) THEN
      ERROR STOP "Prepared weighted inverse differs from explicit phase sum"
   END IF
   CALL k_grid_to_cell_fft(round_k(:, :, 1:9), xkp(:, 1:9), nkp_grid, index_to_cell, weighted_result, &
                           used_fft, weights=subset_weights, allow_incomplete=.TRUE.)
   IF (MAXVAL(ABS(weighted_result - weighted_reference)) > 3.0E-12_dp) THEN
      ERROR STOP "One-shot and prepared inverse transforms disagree"
   END IF
   CALL k_grid_to_cell_prepare(inverse_work, xkp(:, 1:9), nkp_grid, index_to_cell, &
                               weights=2.0_dp*subset_weights, allow_incomplete=.TRUE.)
   CALL k_grid_to_cell_execute(inverse_work, round_k(:, :, 1:9), weighted_result, used_fft)
   IF (MAXVAL(ABS(weighted_result - 2.0_dp*weighted_reference)) > 6.0E-12_dp) THEN
      ERROR STOP "Repreparing inverse transform retained stale weights"
   END IF
   CALL k_grid_to_cell_release(inverse_work)
   CALL k_grid_to_cell_release(inverse_work)

   nkp_grid = [4, 4, 4]

   xkp_bad = xkp
   xkp_bad(1, nkp) = xkp_bad(1, nkp) + 0.03125_dp
   CALL cell_to_k_grid_fft(values_rs, index_to_cell, xkp_bad, nkp_grid, transformed, used_fft)
   IF (used_fft) ERROR STOP "Irregular k-point list incorrectly used lattice FFT"
   CALL direct_reference(values_rs, index_to_cell, xkp_bad, direct)
   IF (MAXVAL(ABS(transformed - direct)) > 1.0E-14_dp) THEN
      ERROR STOP "Direct lattice-transform fallback is inconsistent"
   END IF

   ! Verify the block-circulant inverse against direct matrix application and right covariance.
   base_dims = [4, 2, 1]
   CALL lattice_grid(base_dims, base_index, base_xkp)
   base_overlap = 0.0_dp
   DO jcell = 1, base_cell
      DO icell = 1, base_cell
         delta = MODULO(base_index(:, icell) - base_index(:, jcell), base_dims)
         DO i = 1, base_block
            IF (ALL(delta == 0)) THEN
               base_overlap((icell - 1)*base_block + i, (jcell - 1)*base_block + i) = &
                  2.8_dp + 0.2_dp*REAL(i, KIND=dp)
            ELSE IF ((delta(1) == 1 .OR. delta(1) == base_dims(1) - 1) .AND. &
                     delta(2) == 0) THEN
               base_overlap((icell - 1)*base_block + i, (jcell - 1)*base_block + i) = -0.12_dp
            ELSE IF (delta(1) == 0 .AND. delta(2) == 1) THEN
               base_overlap((icell - 1)*base_block + i, (jcell - 1)*base_block + i) = -0.08_dp
            END IF
         END DO
         IF (ALL(delta == 0)) THEN
            base_overlap((icell - 1)*base_block + 1, (jcell - 1)*base_block + 2) = 0.15_dp
            base_overlap((icell - 1)*base_block + 2, (jcell - 1)*base_block + 1) = 0.15_dp
         END IF
      END DO
   END DO
   DO j = 1, base_ncol
      DO i = 1, base_n
         base_input(i, j) = SIN(0.11_dp*REAL(2*i + 3*j, KIND=dp))
      END DO
   END DO
   CALL build_lattice_inverse(base_overlap, base_dims, inverse_k, projection_error, used_fft, info)
   IF (info /= 0) ERROR STOP "Positive lattice overlap failed block Cholesky inversion"
   IF (.NOT. used_fft) ERROR STOP "Block-circulant overlap construction did not use lattice FFT"
   IF (projection_error > 2.0E-15_dp) ERROR STOP "Exact block-circulant overlap was changed by projection"
   IF (lattice_fft_auto_select(base_dims, base_block, base_ncol, 1, 1, &
                               projection_error, used_fft, info)) THEN
      ERROR STOP "Cost-aware lattice-FFT AUTO accepted a tiny eight-cell operator"
   END IF
   IF (.NOT. lattice_fft_auto_select(base_dims, 6, 8, 1, 2, &
                                     projection_error, used_fft, info)) THEN
      ERROR STOP "Cost-aware lattice-FFT AUTO rejected an economical eight-cell operator"
   END IF
   IF (lattice_fft_auto_select([4, 1, 1], 6, 8, 1, 2, &
                               projection_error, used_fft, info)) THEN
      ERROR STOP "Conservative lattice-FFT AUTO accepted fewer than eight cells"
   END IF
   IF (lattice_fft_auto_select(base_dims, 6, 8, 1, 2, 2.0E-3_dp, used_fft, info)) THEN
      ERROR STOP "Conservative lattice-FFT AUTO accepted a large projection error"
   END IF
   IF (lattice_fft_auto_select(base_dims, 6, 8, 17, 2, &
                               projection_error, used_fft, info)) THEN
      ERROR STOP "Storage-aware lattice-FFT AUTO accepted too many state operators"
   END IF
   IF (lattice_fft_cost_ratio(base_dims, base_block, base_ncol, 1) <= 0.8_dp) THEN
      ERROR STOP "Lattice-FFT cost estimate missed the tiny-batch overhead"
   END IF
   large_dims = [100, 100, 1]
   large_ncell_real = REAL(PRODUCT(large_dims), KIND=dp)
   large_block_real = REAL(large_block, KIND=dp)
   large_ncol_real = REAL(large_ncol, KIND=dp)
   expected_cost_ratio = ( &
                         large_ncell_real*large_block_real**2*large_ncol_real + &
                         4.0_dp*large_ncell_real*large_block_real*large_ncol_real* &
                         LOG(large_ncell_real)/LOG(2.0_dp) + &
                         8.0_dp*large_ncell_real*large_block_real + &
                         2.0_dp*large_ncell_real*large_block_real*large_ncol_real* &
                         LOG(REAL(large_num_pe + 1, KIND=dp))/LOG(2.0_dp))/ &
                         ((large_ncell_real*large_block_real)**2*large_ncol_real)
   cost_ratio = lattice_fft_cost_ratio(large_dims, large_block, large_ncol, large_num_pe)
   IF (ABS(cost_ratio - expected_cost_ratio) > 1.0E-14_dp*expected_cost_ratio) THEN
      ERROR STOP "Lattice-FFT cost estimate overflowed for a large problem"
   END IF
   IF (ABS(lattice_fft_storage_ratio(base_dims, 1) - 0.25_dp) > 1.0E-14_dp) THEN
      ERROR STOP "Lattice-FFT storage estimate is inconsistent"
   END IF
   IF (ABS(lattice_local_cost_ratio(base_dims, 1) - 0.28125_dp) > 1.0E-14_dp) THEN
      ERROR STOP "Lattice local-correction cost estimate is inconsistent"
   END IF
   IF (ABS(lattice_local_storage_ratio(base_dims, 3, 1) - 0.421875_dp) > 1.0E-14_dp) THEN
      ERROR STOP "Lattice local-correction storage estimate is inconsistent"
   END IF
   CALL apply_lattice_inverse_dense(inverse_k, base_dims, base_input, base_output, &
                                    used_fft, imaginary_residual)
   IF (.NOT. used_fft) ERROR STOP "Block-circulant inverse application did not use lattice FFT"
   IF (imaginary_residual > 3.0E-13_dp) ERROR STOP "Lattice inverse did not reconstruct a real matrix"
   IF (MAXVAL(ABS(MATMUL(base_overlap, base_output) - base_input)) > 2.0E-11_dp) THEN
      ERROR STOP "FFT block-circulant inverse differs from the dense operator inverse"
   END IF

   ALLOCATE (state_inverse_k(base_block, base_block, base_cell, base_ncol))
   DO j = 1, base_ncol
      state_inverse_k(:, :, :, j) = (0.5_dp + 0.25_dp*REAL(j, KIND=dp))*inverse_k
   END DO
   CALL apply_lattice_state_inverse_dense(state_inverse_k, base_dims, base_input, &
                                          base_state_output, used_fft, imaginary_residual)
   IF (.NOT. used_fft) ERROR STOP "State-dependent inverse application did not use lattice FFT"
   DO j = 1, base_ncol
      IF (MAXVAL(ABS(base_state_output(:, j) - &
                     (0.5_dp + 0.25_dp*REAL(j, KIND=dp))*base_output(:, j))) > 3.0E-12_dp) THEN
         ERROR STOP "State-dependent lattice inverse used the wrong column operator"
      END IF
   END DO
   DEALLOCATE (state_inverse_k)

   column_rotation = 0.0_dp
   column_rotation(1, 1) = COS(0.31_dp)
   column_rotation(1, 2) = -SIN(0.31_dp)
   column_rotation(2, 1) = SIN(0.31_dp)
   column_rotation(2, 2) = COS(0.31_dp)
   column_rotation(3, 3) = 1.0_dp
   base_input_rotated = MATMUL(base_input, column_rotation)
   CALL apply_lattice_inverse_dense(inverse_k, base_dims, base_input_rotated, &
                                    base_output_rotated, used_fft)
   IF (MAXVAL(ABS(base_output_rotated - MATMUL(base_output, column_rotation))) > 3.0E-12_dp) THEN
      ERROR STOP "FFT lattice base is not right-covariant under orbital rotations"
   END IF

   ! Add one localized positive defect and verify the balanced cell correction.
   defect_operator = base_overlap
   defect_operator(5, 5) = defect_operator(5, 5) + 0.90_dp
   defect_operator(6, 6) = defect_operator(6, 6) + 0.60_dp
   CALL build_lattice_inverse(defect_operator, base_dims, inverse_k, projection_error, used_fft, info)
   IF (info /= 0 .OR. .NOT. used_fft) ERROR STOP "Defective lattice operator failed construction"
   IF (projection_error <= 1.0E-3_dp) ERROR STOP "Localized defect did not break lattice translations"
   IF (lattice_fft_auto_select(base_dims, 6, 8, 1, 2, projection_error, used_fft, info)) THEN
      ERROR STOP "Lattice-FFT AUTO accepted a localized defect"
   END IF
   CALL build_lattice_local_correction( &
      defect_operator, base_dims, 1, correction_indices, local_operator, coarse_inverse, &
      corrected_projection_error, info)
   IF (info /= 0) ERROR STOP "Local correction coarse inverse failed"
   IF (ANY(correction_indices /= [5, 6])) ERROR STOP "Local correction selected the wrong defect cell"
   IF (corrected_projection_error >= projection_error) THEN
      ERROR STOP "Local correction did not reduce the non-circulant residual"
   END IF

   ALLOCATE (defect_base_output(base_n, base_ncol), defect_corrected_output(base_n, base_ncol))
   CALL apply_lattice_inverse_dense(inverse_k, base_dims, base_input, defect_base_output, used_fft)
   CALL apply_balanced_local_inverse(inverse_k, base_dims, correction_indices, local_operator, &
                                     coarse_inverse, base_input, defect_corrected_output)
   defect_residual = MATMUL(defect_operator, defect_base_output) - base_input
   base_residual_norm = NORM2(defect_residual)
   defect_residual = MATMUL(defect_operator, defect_corrected_output) - base_input
   corrected_residual_norm = NORM2(defect_residual)
   IF (corrected_residual_norm >= base_residual_norm) THEN
      ERROR STOP "Balanced local correction did not improve inverse application"
   END IF
   IF (SUM(base_input*defect_corrected_output) <= 0.0_dp) THEN
      ERROR STOP "Balanced local correction lost positive definiteness"
   END IF

   ALLOCATE (selected_rhs(base_n, SIZE(correction_indices)), &
             selected_output(base_n, SIZE(correction_indices)))
   selected_rhs(:, :) = local_operator
   CALL apply_balanced_local_inverse(inverse_k, base_dims, correction_indices, local_operator, &
                                     coarse_inverse, selected_rhs, selected_output)
   selected_rhs = 0.0_dp
   DO i = 1, SIZE(correction_indices)
      selected_rhs(correction_indices(i), i) = 1.0_dp
   END DO
   IF (MAXVAL(ABS(selected_output - selected_rhs)) > 3.0E-11_dp) THEN
      ERROR STOP "Balanced local correction is not exact in its selected subspace"
   END IF
   DEALLOCATE (coarse_inverse, correction_indices, defect_base_output, defect_corrected_output, &
               local_operator, selected_output, selected_rhs)
   DEALLOCATE (inverse_k)

   CALL finalize_fft(mp_comm, "")
   CALL mp_world_finalize()

CONTAINS

! **************************************************************************************************
!> \brief Verify exact-key reuse and invalidation independently of model kernels and FFT backend.
! **************************************************************************************************
   SUBROUTINE check_mesh_certificate()

      TYPE(kpoint_mesh_certificate_type)                 :: certificate
      INTEGER                                            :: ix, iy, ik
      INTEGER                                            :: cells(3, 6), changed_cells(3, 6)
      REAL(KIND=dp)                                      :: points(3, 6), changed_points(3, 6)
      REAL(KIND=dp)                                      :: weights(6), changed_weights(6)
      LOGICAL                                            :: valid, reused

      ik = 0
      DO iy = 0, 2
         DO ix = 0, 1
            ik = ik + 1
            cells(:, ik) = [ix, iy, 0]
            points(:, ik) = [REAL(ix, dp)/2.0_dp, REAL(iy, dp)/3.0_dp, 0.0_dp]
         END DO
      END DO
      weights = 1.0_dp/6.0_dp
      CALL certify_kpoint_mesh(certificate, points, [2, 3, 1], cells, weights, valid, reused)
      IF (.NOT. valid .OR. reused) ERROR STOP "Initial mesh certificate failed"
      CALL certify_kpoint_mesh(certificate, points, [2, 3, 1], cells, weights, valid, reused)
      IF (.NOT. valid .OR. .NOT. reused) ERROR STOP "Identical mesh certificate was not reused"

      changed_points = points(:, [6, 5, 4, 3, 2, 1]) + 0.125_dp
      CALL certify_kpoint_mesh(certificate, changed_points, [2, 3, 1], cells, weights, valid, reused)
      IF (.NOT. valid .OR. reused) ERROR STOP "Shifted/reordered mesh did not invalidate certificate"
      changed_points(1, 1) = changed_points(1, 1) + 1.0E-13_dp
      CALL certify_kpoint_mesh(certificate, changed_points, [2, 3, 1], cells, weights, valid, reused)
      IF (.NOT. valid .OR. reused) ERROR STOP "Near-equal mesh incorrectly reused certificate"
      CALL certify_kpoint_mesh(certificate, changed_points, [2, 3, 1], cells, weights, valid, reused)
      IF (.NOT. valid .OR. .NOT. reused) ERROR STOP "Recertified mesh did not become reusable"

      changed_cells = cells(:, [6, 5, 4, 3, 2, 1]) + 4
      CALL certify_kpoint_mesh(certificate, points, [2, 3, 1], changed_cells, weights, valid, reused)
      IF (.NOT. valid .OR. reused) ERROR STOP "Changed cell labels incorrectly reused certificate"
      changed_weights = weights
      changed_weights(1) = changed_weights(1) + 1.0E-14_dp
      changed_weights(2) = changed_weights(2) - 1.0E-14_dp
      CALL certify_kpoint_mesh(certificate, points, [2, 3, 1], changed_cells, changed_weights, valid, reused)
      IF (.NOT. valid .OR. reused) ERROR STOP "Near-equal weights incorrectly reused certificate"

      changed_points = points
      changed_points(:, 2) = changed_points(:, 1)
      CALL certify_kpoint_mesh(certificate, changed_points, [2, 3, 1], cells, weights, valid, reused)
      IF (valid .OR. reused) ERROR STOP "Duplicate k points received a certificate"
      changed_cells = cells
      changed_cells(:, 2) = changed_cells(:, 1)
      CALL certify_kpoint_mesh(certificate, points, [2, 3, 1], changed_cells, weights, valid, reused)
      IF (valid .OR. reused) ERROR STOP "Duplicate BvK representatives received a certificate"
      changed_weights = weights
      changed_weights(1) = changed_weights(1) + 0.01_dp
      changed_weights(2) = changed_weights(2) - 0.01_dp
      CALL certify_kpoint_mesh(certificate, points, [2, 3, 1], cells, changed_weights, valid, reused)
      IF (valid .OR. reused) ERROR STOP "Nonuniform weights received a certificate"
      CALL certify_kpoint_mesh(certificate, points, [2, 3, 1], cells, 2.0_dp*weights, valid, reused)
      IF (valid .OR. reused) ERROR STOP "Unnormalized weights received a certificate"

      CALL certify_kpoint_mesh(certificate, points, [2, 3, 1], cells, weights, valid, reused)
      IF (.NOT. valid .OR. reused) ERROR STOP "Invalid mesh did not invalidate the previous key"
      changed_points = points
      changed_points(1, 1) = ieee_value(0.0_dp, ieee_quiet_nan)
      CALL certify_kpoint_mesh(certificate, changed_points, [2, 3, 1], cells, weights, valid, reused)
      IF (valid .OR. reused) ERROR STOP "Non-finite point received a certificate"
      changed_weights = weights
      changed_weights(1) = ieee_value(0.0_dp, ieee_quiet_nan)
      CALL certify_kpoint_mesh(certificate, points, [2, 3, 1], cells, changed_weights, valid, reused)
      IF (valid .OR. reused) ERROR STOP "Non-finite weight received a certificate"

      CALL certify_kpoint_mesh(certificate, points, [2, 3, 1], cells, weights, valid, reused)
      IF (.NOT. valid .OR. reused) ERROR STOP "Non-finite inputs left a stale valid certificate"
      CALL certify_kpoint_mesh(certificate, points(:, 1:5), [2, 3, 1], cells, weights, valid, reused)
      IF (valid .OR. reused) ERROR STOP "Inconsistent point shape received a certificate"
      CALL certify_kpoint_mesh(certificate, points, [2, 3, 1], cells(:, 1:5), weights, valid, reused)
      IF (valid .OR. reused) ERROR STOP "Inconsistent cell shape received a certificate"
      CALL certify_kpoint_mesh(certificate, points, [2, 3, 1], cells, weights, valid, reused)
      IF (.NOT. valid .OR. reused) ERROR STOP "Invalid shape did not invalidate the certificate"
      CALL certify_kpoint_mesh(certificate, points, [0, 3, 1], cells, weights, valid, reused)
      IF (valid .OR. reused) ERROR STOP "Invalid grid dimensions received a certificate"
      CALL certify_kpoint_mesh(certificate, points, [2, 3, 1], cells, weights, valid, reused)
      IF (.NOT. valid .OR. reused) ERROR STOP "Invalid dimensions did not invalidate the certificate"

      ! A changed number of points must be safe without accessing stale-shaped cache arrays.
      CALL certify_kpoint_mesh(certificate, points(:, 1:1), [1, 1, 1], cells(:, 1:1), [1.0_dp], valid, reused)
      IF (.NOT. valid .OR. reused) ERROR STOP "Gamma mesh failed after larger-mesh certificate"
      CALL certify_kpoint_mesh(certificate, points, [2, 3, 1], cells, weights, valid, reused)
      IF (.NOT. valid .OR. reused) ERROR STOP "Larger mesh failed after Gamma certificate"

   END SUBROUTINE check_mesh_certificate

! **************************************************************************************************
!> \brief Prime-length full/reduced grids use FFTW, with an exact FFTSG fallback.
! **************************************************************************************************
   SUBROUTINE check_prime_grid()

      INTEGER, PARAMETER                               :: np = 343, ns = 3
      INTEGER                                           :: k, x, y, z, component
      INTEGER, DIMENSION(ns)                            :: selected
      REAL(KIND=dp), DIMENSION(3, np)                   :: points
      REAL(KIND=dp), DIMENSION(3, ns)                   :: bad_points
      REAL(KIND=dp), DIMENSION(ns)                      :: weights
      COMPLEX(KIND=dp), DIMENSION(2, 3, np)             :: exact, RESULT
      COMPLEX(KIND=dp), DIMENSION(2, 3, ns)             :: reduced
      COMPLEX(KIND=dp), DIMENSION(2, 3, ncell)          :: inverse, expected_inverse
      LOGICAL                                           :: fft_used, capable

      k = 0
      DO z = 0, 6
         DO y = 0, 6
            DO x = 0, 6
               k = k + 1
               points(:, k) = ([REAL(x, dp), REAL(y, dp), REAL(z, dp)] + &
                               [0.5_dp, 0.0_dp, 0.25_dp])/7.0_dp
            END DO
         END DO
      END DO
      selected = [2, 51, 319]
      capable = fft_supports_arbitrary_lengths()
      CALL cell_to_k_grid_fft(values_rs, index_to_cell, points, [7, 7, 7], RESULT, fft_used)
      IF (fft_used .NEQV. capable) ERROR STOP "Prime-grid forward path ignored backend capabilities"
      CALL direct_reference(values_rs, index_to_cell, points, exact)
      IF (MAXVAL(ABS(RESULT - exact)) > 3.0E-12_dp) ERROR STOP "Prime-grid forward transform differs"
      CALL cell_to_k_grid_fft(values_rs, index_to_cell, points(:, selected), [7, 7, 7], &
                              reduced, fft_used, allow_incomplete=.TRUE.)
      IF (fft_used .NEQV. capable) ERROR STOP "Prime reduced grid ignored backend capabilities"
      IF (MAXVAL(ABS(reduced - exact(:, :, selected))) > 3.0E-12_dp) THEN
         ERROR STOP "Prime reduced transform differs"
      END IF

      weights = [0.2_dp, 0.3_dp, 0.5_dp]
      expected_inverse = z_zero
      DO component = 1, ncell
         DO k = 1, ns
            expected_inverse(:, :, component) = expected_inverse(:, :, component) + weights(k)* &
                                                EXP(-gaussi*twopi*SUM(points(:, selected(k))* &
                                                                      REAL(index_to_cell(:, component), dp)))* &
                                                reduced(:, :, k)
         END DO
      END DO
      CALL k_grid_to_cell_fft(reduced, points(:, selected), [7, 7, 7], index_to_cell, inverse, &
                              fft_used, weights=weights, allow_incomplete=.TRUE.)
      IF (fft_used .NEQV. capable) ERROR STOP "Prime inverse path ignored backend capabilities"
      IF (MAXVAL(ABS(inverse - expected_inverse)) > 3.0E-12_dp) THEN
         ERROR STOP "Prime weighted inverse differs"
      END IF

      DO component = 1, 3
         CALL direct_reference(values_rs, index_to_cell, points, exact, component, hmat)
         CALL cell_to_k_grid_fft(values_rs, index_to_cell, points(:, selected), [7, 7, 7], &
                                 reduced, fft_used, component, hmat, allow_incomplete=.TRUE.)
         IF (fft_used .NEQV. capable) ERROR STOP "Prime derivative path ignored backend capabilities"
         IF (MAXVAL(ABS(reduced - exact(:, :, selected))) > 1.0E-11_dp) THEN
            ERROR STOP "Prime reduced derivative differs"
         END IF
      END DO

      bad_points = points(:, selected)
      bad_points(1, 2) = bad_points(1, 2) + 0.03_dp
      CALL cell_to_k_grid_fft(values_rs, index_to_cell, bad_points, [7, 7, 7], &
                              reduced, fft_used, allow_incomplete=.TRUE.)
      IF (fft_used) ERROR STOP "Irregular reduced list used FFT"
      CALL direct_reference(values_rs, index_to_cell, bad_points, exact(:, :, 1:ns))
      IF (MAXVAL(ABS(reduced - exact(:, :, 1:ns))) > 1.0E-14_dp) THEN
         ERROR STOP "Irregular reduced fallback differs"
      END IF

   END SUBROUTINE check_prime_grid

! **************************************************************************************************
!> \brief Apply the balanced two-level inverse used by the production local correction.
!> \param inverse_k block-circulant base inverse
!> \param lattice_dims replicated-cell dimensions
!> \param correction_indices selected coordinate basis
!> \param operator_columns columns A U
!> \param coarse_inverse inverse of U^T A U
!> \param input right-hand sides
!> \param output corrected inverse action
! **************************************************************************************************
   SUBROUTINE apply_balanced_local_inverse(inverse_k, lattice_dims, correction_indices, &
                                           operator_columns, coarse_inverse, input, output)

      COMPLEX(KIND=dp), DIMENSION(:, :, :), INTENT(IN)   :: inverse_k
      INTEGER, DIMENSION(3), INTENT(IN)                  :: lattice_dims
      INTEGER, DIMENSION(:), INTENT(IN)                  :: correction_indices
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: operator_columns, coarse_inverse, input
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: output

      INTEGER                                            :: i
      LOGICAL                                            :: used_fft
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: base_input, coarse_input, &
                                                            coarse_projected, coarse_residual, &
                                                            coarse_rhs

      ALLOCATE (base_input(SIZE(input, 1), SIZE(input, 2)), &
                coarse_input(SIZE(correction_indices), SIZE(input, 2)), &
                coarse_projected(SIZE(correction_indices), SIZE(input, 2)), &
                coarse_residual(SIZE(correction_indices), SIZE(input, 2)), &
                coarse_rhs(SIZE(correction_indices), SIZE(input, 2)))
      DO i = 1, SIZE(correction_indices)
         coarse_input(i, :) = input(correction_indices(i), :)
      END DO
      coarse_rhs(:, :) = MATMUL(coarse_inverse, coarse_input)
      base_input(:, :) = input - MATMUL(operator_columns, coarse_rhs)
      CALL apply_lattice_inverse_dense(inverse_k, lattice_dims, base_input, output, used_fft)
      IF (.NOT. used_fft) ERROR STOP "Balanced correction base did not use lattice FFT"
      coarse_projected(:, :) = MATMUL(TRANSPOSE(operator_columns), output)
      coarse_residual(:, :) = MATMUL(coarse_inverse, coarse_projected)
      output(correction_indices, :) = output(correction_indices, :) + coarse_rhs - coarse_residual
      DEALLOCATE (base_input, coarse_input, coarse_projected, coarse_residual, coarse_rhs)

   END SUBROUTINE apply_balanced_local_inverse

! **************************************************************************************************
!> \brief ...
!> \param values_rs ...
!> \param index_to_cell ...
!> \param xkp ...
!> \param values_k ...
!> \param deriv_direction ...
!> \param hmat ...
! **************************************************************************************************
   SUBROUTINE direct_reference(values_rs, index_to_cell, xkp, values_k, deriv_direction, hmat)

      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: values_rs
      INTEGER, DIMENSION(:, :), INTENT(IN)               :: index_to_cell
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: xkp
      COMPLEX(KIND=dp), DIMENSION(:, :, :), INTENT(OUT)  :: values_k
      INTEGER, INTENT(IN), OPTIONAL                      :: deriv_direction
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN), &
         OPTIONAL                                        :: hmat

      COMPLEX(KIND=dp)                                   :: factor
      INTEGER                                            :: icell, ik
      REAL(KIND=dp), DIMENSION(3)                        :: cell_index_real, cell_vector

      values_k = z_zero
      DO ik = 1, SIZE(xkp, 2)
         DO icell = 1, SIZE(values_rs, 3)
            factor = EXP(gaussi*twopi*SUM(xkp(:, ik)*REAL(index_to_cell(:, icell), KIND=dp)))
            IF (PRESENT(deriv_direction)) THEN
               cell_index_real = REAL(index_to_cell(:, icell), KIND=dp)
               cell_vector = MATMUL(hmat, cell_index_real)
               factor = factor*gaussi*cell_vector(deriv_direction)
            END IF
            values_k(:, :, ik) = values_k(:, :, ik) + factor*values_rs(:, :, icell)
         END DO
      END DO

   END SUBROUTINE direct_reference

END PROGRAM kpoint_lattice_fft_unittest
