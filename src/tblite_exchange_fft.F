!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
!> Consumer-owned complex BvK Fourier maps and scratch, including their exact adjoints.
MODULE tblite_exchange_fft
#if defined(__TBLITE_GXTB_TRANSFORM)
   USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
   USE fft_tools, ONLY: fft3d, fft_alloc, fft_dealloc, FWFFT, BWFFT
   USE kinds, ONLY: dp
   USE kpoint_lattice_fft, ONLY: regular_kpoint_grid, lattice_fft_shape
   USE mathconstants, ONLY: gaussi, twopi, z_zero
   USE tblite_transform_contract, ONLY: exchange_transform_type, exchange_k_to_r, &
                                            exchange_r_to_k, exchange_r_to_k_adjoint, exchange_k_to_r_adjoint
#include "./base/base_uses.f90"
#endif
   IMPLICIT NONE
   PRIVATE
#if defined(__TBLITE_GXTB_TRANSFORM)
   PUBLIC :: cp2k_exchange_transform_type, prepare_exchange_transform, release_exchange_transform, &
             exchange_transform_counts

   TYPE, EXTENDS(exchange_transform_type) :: cp2k_exchange_transform_type
      PRIVATE
      LOGICAL :: ready = .FALSE., use_fft = .FALSE., direct_requested = .FALSE.
      INTEGER :: nmesh(3) = 0, nfft(3) = 0, fft_calls = 0, direct_calls = 0
      INTEGER, ALLOCATABLE :: cells(:, :), grid_index(:, :), cell_index(:, :)
      REAL(KIND=dp), ALLOCATABLE :: points(:, :), weights(:)
      COMPLEX(KIND=dp), ALLOCATABLE :: cell_phase(:)
      COMPLEX(KIND=dp), CONTIGUOUS, POINTER :: fft_in(:, :, :) => NULL(), fft_out(:, :, :) => NULL()
   CONTAINS
      PROCEDURE :: apply => apply_exchange_transform
      FINAL :: finalize_exchange_transform
   END TYPE cp2k_exchange_transform_type

CONTAINS
!> Prepare an exact-key consumer map. Near-regular grids fall back to their actual phase sums.
   SUBROUTINE prepare_exchange_transform(self, points, nmesh, cells, weights, force_direct)
      TYPE(cp2k_exchange_transform_type), INTENT(INOUT)  :: self
      REAL(KIND=dp), INTENT(IN)                          :: points(:, :)
      INTEGER, INTENT(IN)                                :: nmesh(3), cells(:, :)
      REAL(KIND=dp), INTENT(IN)                          :: weights(:)
      LOGICAL, INTENT(IN), OPTIONAL                      :: force_direct

      INTEGER                                            :: axis, cell, k, nk
      LOGICAL                                            :: direct, same
      REAL(KIND=dp)                                      :: offset(3), scaled, tolerance

      direct = .FALSE.
      IF (PRESENT(force_direct)) direct = force_direct
      nk = SIZE(weights)
      CPASSERT(nk > 0 .AND. ALL(nmesh > 0))
      CPASSERT(ALL(SHAPE(points) == [3, nk]) .AND. ALL(SHAPE(cells) == [3, nk]))
      CPASSERT(ALL(ieee_is_finite(points)) .AND. ALL(ieee_is_finite(weights)))
      same = .FALSE.
      IF (self%ready) THEN
         IF (ALL(self%nmesh == nmesh) .AND. SIZE(self%weights) == nk) THEN
            same = ALL(self%points == points) .AND. ALL(self%cells == cells) .AND. &
                   ALL(self%weights == weights) .AND. (self%direct_requested .EQV. direct)
         END IF
      END IF
      IF (same) RETURN
      CALL release_exchange_transform(self)
      self%points = points
      self%cells = cells
      self%weights = weights
      self%nmesh = nmesh
      self%direct_requested = direct
      ALLOCATE (self%grid_index(3, nk))
      self%use_fft = regular_kpoint_grid(points, nmesh, self%grid_index, offset)
      ! The generic grid mapper tolerates numerical geometry noise. Model exchange uses
      ! a tighter machine-scale gate, retaining the original points if it fails.
      IF (self%use_fft) THEN
         DO k = 1, nk
            DO axis = 1, 3
               scaled = REAL(nmesh(axis), dp)*(points(axis, k) - offset(axis))
               tolerance = 64.0_dp*EPSILON(1.0_dp)*MAX(1.0_dp, ABS(scaled))
               IF (ABS(scaled - REAL(NINT(scaled), dp)) > tolerance) self%use_fft = .FALSE.
            END DO
         END DO
      END IF
      IF (self%use_fft) self%use_fft = lattice_fft_shape(nmesh, self%nfft, allow_arbitrary=.TRUE.)
      IF (direct) self%use_fft = .FALSE.
      IF (self%use_fft) THEN
         ALLOCATE (self%cell_index(3, nk), self%cell_phase(nk))
         DO cell = 1, nk
            DO axis = 1, 3
               self%cell_index(axis, cell) = MODULO(cells(axis, cell), nmesh(axis))* &
                                             (self%nfft(axis)/nmesh(axis)) + 1
            END DO
            self%cell_phase(cell) = EXP(-gaussi*twopi*SUM(offset*REAL(cells(:, cell), dp)))
         END DO
         CALL fft_alloc(self%fft_in, self%nfft)
         CALL fft_alloc(self%fft_out, self%nfft)
      END IF
      self%ready = .TRUE.
   END SUBROUTINE prepare_exchange_transform

!> Apply forward/inverse/adjoint maps outside OpenMP regions; scratch belongs to this instance.
   SUBROUTINE apply_exchange_transform(self, source, TARGET, operation)
      CLASS(cp2k_exchange_transform_type), INTENT(INOUT) :: self
      COMPLEX(KIND=dp), INTENT(IN) :: source(:, :)
      COMPLEX(KIND=dp), INTENT(OUT) :: target(:, :)
      INTEGER, INTENT(IN) :: operation
      INTEGER :: i, k, cell, nk, stat
      LOGICAL :: inverse, weighted
      COMPLEX(KIND=dp) :: factor

      CPASSERT(self%ready)
      nk = SIZE(self%weights)
      CPASSERT(ALL(SHAPE(source) == SHAPE(TARGET)) .AND. SIZE(source, 2) == nk)
      CPASSERT(operation >= 1 .AND. operation <= 4)
      inverse = operation == exchange_k_to_r .OR. operation == exchange_r_to_k_adjoint
      weighted = operation == exchange_k_to_r .OR. operation == exchange_k_to_r_adjoint
      stat = 0
      IF (self%use_fft) THEN
         DO i = 1, SIZE(source, 1)
            self%fft_in = z_zero
            IF (inverse) THEN
               DO k = 1, nk
                  factor = 1.0_dp
                  IF (weighted) factor = self%weights(k)
                  self%fft_in(self%grid_index(1, k), self%grid_index(2, k), self%grid_index(3, k)) = &
                     factor*source(i, k)
               END DO
               CALL fft3d(FWFFT, self%nfft, self%fft_in, self%fft_out, status=stat)
               IF (stat /= 0) EXIT
               DO cell = 1, nk
                  TARGET(i, cell) = REAL(PRODUCT(self%nfft), dp)*self%cell_phase(cell)* &
                                    self%fft_out(self%cell_index(1, cell), self%cell_index(2, cell), self%cell_index(3, cell))
               END DO
            ELSE
               DO cell = 1, nk
                  self%fft_in(self%cell_index(1, cell), self%cell_index(2, cell), self%cell_index(3, cell)) = &
                     CONJG(self%cell_phase(cell))*source(i, cell)
               END DO
               CALL fft3d(BWFFT, self%nfft, self%fft_in, self%fft_out, status=stat)
               IF (stat /= 0) EXIT
               DO k = 1, nk
                  factor = 1.0_dp
                  IF (weighted) factor = self%weights(k)
                  TARGET(i, k) = factor* &
                                 self%fft_out(self%grid_index(1, k), self%grid_index(2, k), self%grid_index(3, k))
               END DO
            END IF
         END DO
         IF (stat == 0) THEN
            self%fft_calls = self%fft_calls + 1
            RETURN
         END IF
         self%use_fft = .FALSE.
      END IF
      ! Whole-batch fallback replaces all partially computed outputs after a backend failure.
      TARGET = z_zero
      DO cell = 1, nk
         DO k = 1, nk
            factor = EXP(gaussi*twopi*DOT_PRODUCT(self%points(:, k), REAL(self%cells(:, cell), dp)))
            IF (inverse) factor = CONJG(factor)
            IF (weighted) factor = factor*self%weights(k)
            IF (inverse) THEN
               TARGET(:, cell) = TARGET(:, cell) + factor*source(:, k)
            ELSE
               TARGET(:, k) = TARGET(:, k) + factor*source(:, cell)
            END IF
         END DO
      END DO
      self%direct_calls = self%direct_calls + 1
   END SUBROUTINE apply_exchange_transform

!> Release scratch without resetting lifetime execution counters.
   SUBROUTINE release_exchange_transform(self)
      TYPE(cp2k_exchange_transform_type), INTENT(INOUT)  :: self

      IF (ASSOCIATED(self%fft_in)) CALL fft_dealloc(self%fft_in)
      IF (ASSOCIATED(self%fft_out)) CALL fft_dealloc(self%fft_out)
      IF (ALLOCATED(self%points)) DEALLOCATE (self%points)
      IF (ALLOCATED(self%cells)) DEALLOCATE (self%cells)
      IF (ALLOCATED(self%weights)) DEALLOCATE (self%weights)
      IF (ALLOCATED(self%grid_index)) DEALLOCATE (self%grid_index)
      IF (ALLOCATED(self%cell_index)) DEALLOCATE (self%cell_index)
      IF (ALLOCATED(self%cell_phase)) DEALLOCATE (self%cell_phase)
      self%ready = .FALSE.
      self%use_fft = .FALSE.
   END SUBROUTINE release_exchange_transform

   SUBROUTINE finalize_exchange_transform(self)
      TYPE(cp2k_exchange_transform_type), INTENT(INOUT)  :: self

      CALL release_exchange_transform(self)
   END SUBROUTINE finalize_exchange_transform

   SUBROUTINE exchange_transform_counts(self, fft_calls, direct_calls)
      TYPE(cp2k_exchange_transform_type), INTENT(IN)     :: self
      INTEGER, INTENT(OUT)                               :: fft_calls, direct_calls

      fft_calls = self%fft_calls
      direct_calls = self%direct_calls
   END SUBROUTINE exchange_transform_counts
#endif
END MODULE tblite_exchange_fft
