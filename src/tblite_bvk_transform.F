#if defined(__TBLITE_GXTB)
! CP2K-owned regular-grid transforms, relocated from the extended provider.
! Derived from save_tblite at aad26158b2e83f84b7b7d5b56f361ead4bc02a88.
! SPDX-Identifier: LGPL-3.0-or-later

!> Internal regular-grid Fourier transforms for BvK exchange meshes.
MODULE tblite_bvk_transform
   USE iso_fortran_env,                 ONLY: int64
   USE mctc_env,                        ONLY: wp

   IMPLICIT NONE
   PRIVATE

   !> Dense phase-table/BLAS reference backend.
   INTEGER, PARAMETER, PUBLIC :: bvk_transform_dense = 1
   !> Separable direct-DFT backend.  This is deliberately not called an FFT:
   !> it factorizes a three-dimensional transform into one-dimensional dense
   !> DFTs, reducing O(Nk**2) character products to O(Nk*sum(nmesh)).
   INTEGER, PARAMETER, PUBLIC :: bvk_transform_separable_dft = 2
   !> Internal mixed-radix Cooley--Tukey FFT backend.
   !>
   !> A length-n axis costs O(n*sum(radices)) operations.  This is O(n log n)
   !> when the prime radices are bounded; an unfactorable prime axis uses an
   !> explicit O(n**2) leaf DFT.  For B rows on an Nk-point grid the full cost
   !> is O(B*Nk*sum_i(sum(radices_i))).  No external FFT library is required.
   INTEGER, PARAMETER, PUBLIC :: bvk_transform_mixed_radix_fft = 3

   PUBLIC :: bvk_separable_k_to_r, bvk_separable_r_to_k
   PUBLIC :: bvk_mixed_radix_k_to_r, bvk_mixed_radix_r_to_k

CONTAINS

!> Weighted inverse transform from arbitrary input k order to representative
!> order.  With `normalized=.true.`, the uniform 1/Nk BZ normalization is
!> applied.  Setting it false gives the exact adjoint of the unnormalised
!> R-to-k transform.
   SUBROUTINE bvk_separable_k_to_r(mesh_k, mesh_r, nmesh, input_to_grid, reps, &
      & phase_roots, twist_phase, normalized)
      COMPLEX(wp), INTENT(in)                            :: mesh_k(:, :)
      COMPLEX(wp), INTENT(out)                           :: mesh_r(:, :)
      INTEGER, INTENT(in)                                :: nmesh(3), input_to_grid(:), reps(:, :)
      COMPLEX(wp), INTENT(in)                            :: phase_roots(:, :), twist_phase(:)
      LOGICAL, INTENT(in)                                :: normalized

      COMPLEX(wp), ALLOCATABLE                           :: work_a(:, :), work_b(:, :)
      INTEGER                                            :: icell, igrid, ik, nk
      REAL(wp)                                           :: scale

      nk = PRODUCT(nmesh)
      ALLOCATE (work_a(SIZE(mesh_k, 1), nk), work_b(SIZE(mesh_k, 1), nk), &
         & source=(0.0_wp, 0.0_wp))
      scale = 1.0_wp
      IF (normalized) scale = 1.0_wp/REAL(nk, wp)

      ! The external k order may be arbitrary.  Canonical grid order has the
      ! first reciprocal direction as its fastest varying index.
      DO ik = 1, nk
         work_a(:, input_to_grid(ik)) = scale*mesh_k(:, ik)
      END DO
      CALL transform_axis(work_a, work_b, nmesh, 1, phase_roots, .FALSE.)
      CALL transform_axis(work_b, work_a, nmesh, 2, phase_roots, .FALSE.)
      CALL transform_axis(work_a, work_b, nmesh, 3, phase_roots, .FALSE.)

      ! Representatives need not be returned in canonical grid order.  The
      ! common shifted-mesh character is applied only once after factorization.
      DO icell = 1, nk
         igrid = representative_grid_index(reps(:, icell), nmesh)
         mesh_r(:, icell) = CONJG(twist_phase(icell))*work_b(:, igrid)
      END DO
   END SUBROUTINE bvk_separable_k_to_r

!> Unnormalised forward transform from arbitrary representative order to the
!> original input k order.  This is the exact inverse/adjoint partner of the
!> normalized k-to-R transform; callers attach k weights to reverse-source
!> adjoints after this routine.
   SUBROUTINE bvk_separable_r_to_k(mesh_r, mesh_k, nmesh, input_to_grid, reps, &
      & phase_roots, twist_phase)
      COMPLEX(wp), INTENT(in)                            :: mesh_r(:, :)
      COMPLEX(wp), INTENT(out)                           :: mesh_k(:, :)
      INTEGER, INTENT(in)                                :: nmesh(3), input_to_grid(:), reps(:, :)
      COMPLEX(wp), INTENT(in)                            :: phase_roots(:, :), twist_phase(:)

      COMPLEX(wp), ALLOCATABLE                           :: work_a(:, :), work_b(:, :)
      INTEGER                                            :: icell, igrid, ik, nk

      nk = PRODUCT(nmesh)
      ALLOCATE (work_a(SIZE(mesh_r, 1), nk), work_b(SIZE(mesh_r, 1), nk), &
         & source=(0.0_wp, 0.0_wp))
      DO icell = 1, nk
         igrid = representative_grid_index(reps(:, icell), nmesh)
         work_a(:, igrid) = twist_phase(icell)*mesh_r(:, icell)
      END DO
      CALL transform_axis(work_a, work_b, nmesh, 1, phase_roots, .TRUE.)
      CALL transform_axis(work_b, work_a, nmesh, 2, phase_roots, .TRUE.)
      CALL transform_axis(work_a, work_b, nmesh, 3, phase_roots, .TRUE.)
      DO ik = 1, nk
         mesh_k(:, ik) = work_b(:, input_to_grid(ik))
      END DO
   END SUBROUTINE bvk_separable_r_to_k

!> Normalized inverse mixed-radix FFT from arbitrary input k order to
!> arbitrary representative order.  Setting `normalized` false evaluates the
!> unnormalised adjoint needed by the real reverse exchange map.
   SUBROUTINE bvk_mixed_radix_k_to_r(mesh_k, mesh_r, nmesh, input_to_grid, reps, &
      & phase_roots, twist_phase, normalized)
      COMPLEX(wp), INTENT(in)                            :: mesh_k(:, :)
      COMPLEX(wp), INTENT(out)                           :: mesh_r(:, :)
      INTEGER, INTENT(in)                                :: nmesh(3), input_to_grid(:), reps(:, :)
      COMPLEX(wp), INTENT(in)                            :: phase_roots(:, :), twist_phase(:)
      LOGICAL, INTENT(in)                                :: normalized

      COMPLEX(wp), ALLOCATABLE                           :: work_a(:, :), work_b(:, :)
      INTEGER                                            :: icell, igrid, ik, nk
      REAL(wp)                                           :: scale

      nk = PRODUCT(nmesh)
      ALLOCATE (work_a(SIZE(mesh_k, 1), nk), work_b(SIZE(mesh_k, 1), nk), &
         & source=(0.0_wp, 0.0_wp))
      scale = 1.0_wp
      IF (normalized) scale = 1.0_wp/REAL(nk, wp)

      DO ik = 1, nk
         work_a(:, input_to_grid(ik)) = scale*mesh_k(:, ik)
      END DO
      CALL transform_axis_mixed_radix(work_a, work_b, nmesh, 1, &
         & phase_roots, .FALSE.)
      CALL transform_axis_mixed_radix(work_b, work_a, nmesh, 2, &
         & phase_roots, .FALSE.)
      CALL transform_axis_mixed_radix(work_a, work_b, nmesh, 3, &
         & phase_roots, .FALSE.)

      DO icell = 1, nk
         igrid = representative_grid_index(reps(:, icell), nmesh)
         mesh_r(:, icell) = CONJG(twist_phase(icell))*work_b(:, igrid)
      END DO
   END SUBROUTINE bvk_mixed_radix_k_to_r

!> Unnormalised forward mixed-radix FFT from arbitrary representative order
!> to the original input k order.  Together with the normalized inverse this
!> is the exact inverse and complex/real-adjoint partner used by exchange.
   SUBROUTINE bvk_mixed_radix_r_to_k(mesh_r, mesh_k, nmesh, input_to_grid, reps, &
      & phase_roots, twist_phase)
      COMPLEX(wp), INTENT(in)                            :: mesh_r(:, :)
      COMPLEX(wp), INTENT(out)                           :: mesh_k(:, :)
      INTEGER, INTENT(in)                                :: nmesh(3), input_to_grid(:), reps(:, :)
      COMPLEX(wp), INTENT(in)                            :: phase_roots(:, :), twist_phase(:)

      COMPLEX(wp), ALLOCATABLE                           :: work_a(:, :), work_b(:, :)
      INTEGER                                            :: icell, igrid, ik, nk

      nk = PRODUCT(nmesh)
      ALLOCATE (work_a(SIZE(mesh_r, 1), nk), work_b(SIZE(mesh_r, 1), nk), &
         & source=(0.0_wp, 0.0_wp))
      DO icell = 1, nk
         igrid = representative_grid_index(reps(:, icell), nmesh)
         work_a(:, igrid) = twist_phase(icell)*mesh_r(:, icell)
      END DO
      CALL transform_axis_mixed_radix(work_a, work_b, nmesh, 1, &
         & phase_roots, .TRUE.)
      CALL transform_axis_mixed_radix(work_b, work_a, nmesh, 2, &
         & phase_roots, .TRUE.)
      CALL transform_axis_mixed_radix(work_a, work_b, nmesh, 3, &
         & phase_roots, .TRUE.)
      DO ik = 1, nk
         mesh_k(:, ik) = work_b(:, input_to_grid(ik))
      END DO
   END SUBROUTINE bvk_mixed_radix_r_to_k

!> Apply a mixed-radix FFT to every pencil of one canonical mesh axis.
!>
!> All independent data rows in a pencil are transformed together.  Only one
!> pencil is staged at a time, bounding the recursive scratch by the longest
!> mesh direction rather than by the complete three-dimensional mesh.
   SUBROUTINE transform_axis_mixed_radix(input, output, nmesh, idim, &
      & phase_roots, forward)
      COMPLEX(wp), INTENT(in)                            :: input(:, :)
      COMPLEX(wp), INTENT(out)                           :: output(:, :)
      INTEGER, INTENT(in)                                :: nmesh(3), idim
      COMPLEX(wp), INTENT(in)                            :: phase_roots(:, :)
      LOGICAL, INTENT(in)                                :: forward

      COMPLEX(wp), ALLOCATABLE                           :: pencil_in(:, :), pencil_out(:, :)
      INTEGER                                            :: coord, g(3), i1, i2, i3, index, n

      n = nmesh(idim)
      ALLOCATE (pencil_in(SIZE(input, 1), n), pencil_out(SIZE(input, 1), n))
      DO i3 = 0, MERGE(0, nmesh(3) - 1, idim == 3)
         DO i2 = 0, MERGE(0, nmesh(2) - 1, idim == 2)
            DO i1 = 0, MERGE(0, nmesh(1) - 1, idim == 1)
               g = [i1, i2, i3]
               DO coord = 0, n - 1
                  g(idim) = coord
                  index = 1 + g(1) + nmesh(1)*(g(2) + nmesh(2)*g(3))
                  pencil_in(:, coord + 1) = input(:, index)
               END DO
               CALL mixed_radix_fft_1d(pencil_in, pencil_out, n, 1, &
                  & phase_roots(:, idim), forward)
               DO coord = 0, n - 1
                  g(idim) = coord
                  index = 1 + g(1) + nmesh(1)*(g(2) + nmesh(2)*g(3))
                  output(:, index) = pencil_out(:, coord + 1)
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE transform_axis_mixed_radix

!> Recursive decimation-in-time mixed-radix Cooley--Tukey transform.
!>
!> `root_stride` embeds the current subtransform in the roots cached for the
!> full axis.  Composite lengths are always decomposed.  The direct routine is
!> reached only for a prime leaf, including small prime lengths.
   RECURSIVE SUBROUTINE mixed_radix_fft_1d(input, output, full_n, root_stride, &
      & roots, forward)
      COMPLEX(wp), INTENT(in)                            :: input(:, :)
      COMPLEX(wp), INTENT(out)                           :: output(:, :)
      INTEGER, INTENT(in)                                :: full_n, root_stride
      COMPLEX(wp), INTENT(in)                            :: roots(:)
      LOGICAL, INTENT(in)                                :: forward

      COMPLEX(wp)                                        :: phase
      COMPLEX(wp), ALLOCATABLE                           :: child_in(:, :), child_out(:, :), &
                                                            stage(:, :, :)
      INTEGER                                            :: j1, k, k1, m, n, q, radix, root
      INTEGER(int64)                                     :: root64

      n = SIZE(input, 2)
      IF (n == 1) THEN
         output = input
         RETURN
      END IF
      radix = smallest_factor(n)
      IF (radix == n) THEN
         CALL direct_prime_dft(input, output, full_n, root_stride, roots, forward)
         RETURN
      END IF

      m = n/radix
      ALLOCATE (child_in(SIZE(input, 1), m), child_out(SIZE(input, 1), m), &
         & stage(SIZE(input, 1), m, radix))
      DO j1 = 0, radix - 1
         child_in = input(:, j1 + 1:n:radix)
         CALL mixed_radix_fft_1d(child_in, child_out, full_n, &
            & root_stride*radix, roots, forward)
         stage(:, :, j1 + 1) = child_out
      END DO

      DO q = 0, m - 1
         DO k1 = 0, radix - 1
            k = q + m*k1
            output(:, k + 1) = (0.0_wp, 0.0_wp)
            DO j1 = 0, radix - 1
               root64 = MODULO(INT(root_stride, int64)*INT(j1, int64), &
                  & INT(full_n, int64))
               root64 = MODULO(root64*INT(k, int64), INT(full_n, int64))
               root = INT(root64)
               phase = roots(root + 1)
               IF (.NOT. forward) phase = CONJG(phase)
               output(:, k + 1) = output(:, k + 1) &
                  & + phase*stage(:, q + 1, j1 + 1)
            END DO
         END DO
      END DO
   END SUBROUTINE mixed_radix_fft_1d

!> Direct DFT for a prime Cooley--Tukey leaf.
   SUBROUTINE direct_prime_dft(input, output, full_n, root_stride, roots, forward)
      COMPLEX(wp), INTENT(in)                            :: input(:, :)
      COMPLEX(wp), INTENT(out)                           :: output(:, :)
      INTEGER, INTENT(in)                                :: full_n, root_stride
      COMPLEX(wp), INTENT(in)                            :: roots(:)
      LOGICAL, INTENT(in)                                :: forward

      COMPLEX(wp)                                        :: phase
      INTEGER                                            :: j, k, n, root
      INTEGER(int64)                                     :: root64

      n = SIZE(input, 2)
      output = (0.0_wp, 0.0_wp)
      DO k = 0, n - 1
         DO j = 0, n - 1
            root64 = MODULO(INT(root_stride, int64)*INT(j, int64), &
               & INT(full_n, int64))
            root64 = MODULO(root64*INT(k, int64), INT(full_n, int64))
            root = INT(root64)
            phase = roots(root + 1)
            IF (.NOT. forward) phase = CONJG(phase)
            output(:, k + 1) = output(:, k + 1) + phase*input(:, j + 1)
         END DO
      END DO
   END SUBROUTINE direct_prime_dft

!> Return the smallest prime factor without overflowing `factor*factor`.
   pure INTEGER FUNCTION smallest_factor(n) RESULT(factor)
      INTEGER, INTENT(in)                                :: n

      factor = 2
      DO WHILE (factor <= n/factor)
         IF (MODULO(n, factor) == 0) RETURN
         factor = factor + 1
      END DO
      factor = n
   END FUNCTION smallest_factor

!> Apply one one-dimensional DFT to all pencils of a canonical 3D mesh.
   SUBROUTINE transform_axis(input, output, nmesh, idim, phase_roots, forward)
      COMPLEX(wp), INTENT(in)                            :: input(:, :)
      COMPLEX(wp), INTENT(out)                           :: output(:, :)
      INTEGER, INTENT(in)                                :: nmesh(3), idim
      COMPLEX(wp), INTENT(in)                            :: phase_roots(:, :)
      LOGICAL, INTENT(in)                                :: forward

      COMPLEX(wp)                                        :: phase
      INTEGER                                            :: coord_in, coord_out, g(3), i1, i2, i3, &
                                                            idx_in, idx_out, n, root
      INTEGER(int64)                                     :: root64

      n = nmesh(idim)
      output = (0.0_wp, 0.0_wp)
      DO i3 = 0, nmesh(3) - 1
         DO i2 = 0, nmesh(2) - 1
            DO i1 = 0, nmesh(1) - 1
               g = [i1, i2, i3]
               coord_out = g(idim)
               idx_out = 1 + g(1) + nmesh(1)*(g(2) + nmesh(2)*g(3))
               DO coord_in = 0, n - 1
                  g(idim) = coord_in
                  idx_in = 1 + g(1) + nmesh(1)*(g(2) + nmesh(2)*g(3))
                  root64 = MODULO(INT(coord_in, int64)*INT(coord_out, int64), &
                     & INT(n, int64))
                  root = INT(root64)
                  phase = phase_roots(root + 1, idim)
                  IF (.NOT. forward) phase = CONJG(phase)
                  output(:, idx_out) = output(:, idx_out) + phase*input(:, idx_in)
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE transform_axis

   pure INTEGER FUNCTION representative_grid_index(rep, nmesh) RESULT(index)
      INTEGER, INTENT(in)                                :: rep(3), nmesh(3)

      INTEGER                                            :: grid(3)

      grid = MODULO(rep, nmesh)
      index = 1 + grid(1) + nmesh(1)*(grid(2) + nmesh(2)*grid(3))
   END FUNCTION representative_grid_index

END MODULE tblite_bvk_transform

#endif
