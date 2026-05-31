!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief finite-difference checks for the complex OT rotation Frechet kernel
! **************************************************************************************************
PROGRAM qs_ot_frechet_unittest
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_create, dbcsr_distribution_new, dbcsr_distribution_release, dbcsr_distribution_type, &
        dbcsr_finalize, dbcsr_finalize_lib, dbcsr_get_block_p, dbcsr_init_lib, dbcsr_put_block, &
        dbcsr_release, dbcsr_reserve_blocks, dbcsr_type, dbcsr_type_no_symmetry
   USE kinds,                           ONLY: dp
   USE machine,                         ONLY: default_output_unit
   USE mathlib,                         ONLY: diag_complex
   USE message_passing,                 ONLY: mp_comm_type,&
                                              mp_world_finalize,&
                                              mp_world_init
   USE qs_ot,                           ONLY: qs_ot_apply_complex_frechet_dbcsr,&
                                              qs_ot_complex_exp_frechet_kernel

   IMPLICIT NONE

   INTEGER                                            :: io_unit, mynode, nfail
   TYPE(mp_comm_type)                                 :: mp_comm

   nfail = 0
   CALL mp_world_init(mp_comm)
   mynode = mp_comm%mepos
   io_unit = -1
   IF (mynode == 0) io_unit = default_output_unit
   CALL dbcsr_init_lib(mp_comm%get_handle(), io_unit)

   CALL check_offdiag("separated offdiag", 1.3_dp, -0.4_dp, &
                      CMPLX(0.7_dp, -0.2_dp, KIND=dp), 5.0E-12_dp, nfail)
   CALL check_offdiag("near-degenerate offdiag", 0.2_dp, 0.2000003_dp, &
                      CMPLX(-0.3_dp, 0.5_dp, KIND=dp), 5.0E-12_dp, nfail)
   CALL check_offdiag("degenerate offdiag", -0.6_dp, -0.6_dp, &
                      CMPLX(0.4_dp, 0.1_dp, KIND=dp), 5.0E-12_dp, nfail)
   CALL check_diag("diagonal", 0.8_dp, CMPLX(-0.2_dp, 0.35_dp, KIND=dp), &
                   5.0E-8_dp, nfail)
   CALL check_sparse_dbcsr_pattern("sparse DBCSR-style block pattern", 5.0E-14_dp, nfail)
   CALL check_dbcsr_sparse_pattern("DBCSR sparse block pattern", 5.0E-14_dp, nfail, mp_comm)
   CALL check_weighted_kpoint_mermin_gradient("weighted k-point Mermin gradient", 1.0E-7_dp, nfail)
   CALL check_kpoint_mermin_energy_projection("fixed-N Mermin energy projection", 5.0E-11_dp, nfail)
   CALL check_fixed_n_preconditioner_gauge("fixed-N preconditioner gauge", 5.0E-15_dp, nfail)
   CALL check_complex_kpoint_mermin_coupled_gradient("complex k-point Mermin coupled gradient", &
                                                     2.0E-7_dp, nfail)
   CALL check_complex_fixed_n_schur_hessian("complex fixed-N Schur Hessian", &
                                            2.0E-5_dp, nfail)
   CALL check_complex_fixed_n_response_channel("complex fixed-N response channel", &
                                               2.0E-7_dp, nfail)
   CALL check_complex_projector_blended_response_channel("complex projector-blended response channel", &
                                                         5.0E-9_dp, nfail)
   CALL check_complex_projector_physical_response_channel("complex projector physical response channel", &
                                                          5.0E-9_dp, nfail)
   CALL check_mixed_fixed_n_tether_projector("mixed fixed-N/tether Schur projector", &
                                             5.0E-14_dp, nfail)
   CALL check_mixed_fixed_n_tether_cross_projector("mixed fixed-N/tether Schur cross projector", &
                                                   5.0E-14_dp, nfail)
   CALL check_mixed_fixed_n_asymmetric_response_schur("mixed fixed-N asymmetric response Schur", &
                                                      5.0E-14_dp, nfail)
   CALL check_mixed_fixed_n_tether_exact_response("mixed fixed-N/tether exact response", &
                                                  5.0E-12_dp, nfail)
   CALL check_mixed_fixed_n_tether_exact_schur("mixed fixed-N/tether exact Schur", &
                                               5.0E-11_dp, nfail)
   CALL check_complex_projector_rotation_gradient("complex projector rotation gradient", &
                                                  5.0E-8_dp, nfail)
   CALL check_complex_occ_projector_hessian("complex occupation-projector Hessian", &
                                            2.0E-6_dp, nfail)
   CALL check_complex_occ_projector_gauge_separation("complex occupation-projector gauge separation", &
                                                     1.0E-7_dp, nfail)
   CALL check_complex_nondiag_tether_gradient("complex non-diagonal tether gradient", &
                                              5.0E-9_dp, nfail)
   CALL check_complex_canonical_tether_gradient("complex canonical tether gradient", &
                                                5.0E-9_dp, nfail)

   IF (nfail == 0) PRINT *, "qs_ot_frechet_unittest: OK"

   CALL dbcsr_finalize_lib()
   CALL mp_world_finalize()

   IF (nfail > 0) ERROR STOP "qs_ot_frechet_unittest failed"

CONTAINS

! **************************************************************************************************
!> \brief Check L_exp(diag(-i*e1,-i*e2), E_12*z)
! **************************************************************************************************
   SUBROUTINE check_offdiag(label, e1, e2, z, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: e1, e2
      COMPLEX(KIND=dp), INTENT(IN)                       :: z
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      REAL(KIND=dp), PARAMETER                           :: eps = 1.0E-7_dp

      COMPLEX(KIND=dp)                                   :: fd, ref
      COMPLEX(KIND=dp), DIMENSION(2, 2)                  :: m0, meps, u0, ueps
      REAL(KIND=dp)                                      :: err

      m0(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      m0(1, 1) = CMPLX(0.0_dp, -e1, KIND=dp)
      m0(2, 2) = CMPLX(0.0_dp, -e2, KIND=dp)

      meps = m0
      meps(1, 2) = eps*z

      CALL dense_expm2(m0, u0)
      CALL dense_expm2(meps, ueps)

      fd = (ueps(1, 2) - u0(1, 2))/eps
      ref = qs_ot_complex_exp_frechet_kernel(e1, e2)*z
      err = ABS(fd - ref)

      IF (err > tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,ES12.5)') TRIM(label), "failed with error", err
      ELSE
         WRITE (*, '(A,1X,A,1X,ES12.5)') TRIM(label), "OK with error", err
      END IF

   END SUBROUTINE check_offdiag

! **************************************************************************************************
!> \brief Check L_exp(diag(-i*e,-i*e), E_11*z)
! **************************************************************************************************
   SUBROUTINE check_diag(label, e1, z, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: e1
      COMPLEX(KIND=dp), INTENT(IN)                       :: z
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      REAL(KIND=dp), PARAMETER                           :: eps = 1.0E-7_dp

      COMPLEX(KIND=dp)                                   :: fd, ref
      COMPLEX(KIND=dp), DIMENSION(2, 2)                  :: m0, meps, u0, ueps
      REAL(KIND=dp)                                      :: err

      m0(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      m0(1, 1) = CMPLX(0.0_dp, -e1, KIND=dp)
      m0(2, 2) = CMPLX(0.0_dp, -e1, KIND=dp)

      meps = m0
      meps(1, 1) = meps(1, 1) + eps*z

      CALL dense_expm2(m0, u0)
      CALL dense_expm2(meps, ueps)

      fd = (ueps(1, 1) - u0(1, 1))/eps
      ref = qs_ot_complex_exp_frechet_kernel(e1, e1)*z
      err = ABS(fd - ref)

      IF (err > tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,ES12.5)') TRIM(label), "failed with error", err
      ELSE
         WRITE (*, '(A,1X,A,1X,ES12.5)') TRIM(label), "OK with error", err
      END IF

   END SUBROUTINE check_diag

! **************************************************************************************************
!> \brief Check finite-generator Frechet application with mismatched Re/Im block patterns.
!>        This mirrors the DBCSR iterator path where real-only, imaginary-only,
!>        mixed, and absent sparse blocks must be handled consistently.
! **************************************************************************************************
   SUBROUTINE check_sparse_dbcsr_pattern(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 4, nb = 2

      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: dense_ref, sparse_out, z
      LOGICAL, DIMENSION(nb, nb)                         :: mask_im, mask_re
      REAL(KIND=dp)                                      :: err
      REAL(KIND=dp), DIMENSION(n)                        :: evals

      evals = [0.7_dp, -0.2_dp, 1.4_dp, -1.1_dp]
      z(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      z(1, 1) = CMPLX(0.3_dp, -0.4_dp, KIND=dp)
      z(1, 2) = CMPLX(-0.2_dp, 0.1_dp, KIND=dp)
      z(2, 1) = CMPLX(0.5_dp, 0.8_dp, KIND=dp)
      z(2, 2) = CMPLX(-0.7_dp, 0.2_dp, KIND=dp)
      z(1, 3) = CMPLX(0.0_dp, 0.6_dp, KIND=dp)
      z(1, 4) = CMPLX(0.0_dp, -0.3_dp, KIND=dp)
      z(2, 3) = CMPLX(0.0_dp, 0.4_dp, KIND=dp)
      z(2, 4) = CMPLX(0.0_dp, 0.2_dp, KIND=dp)
      z(3, 1) = CMPLX(-0.9_dp, 0.0_dp, KIND=dp)
      z(3, 2) = CMPLX(0.4_dp, 0.0_dp, KIND=dp)
      z(4, 1) = CMPLX(0.2_dp, 0.0_dp, KIND=dp)
      z(4, 2) = CMPLX(-0.1_dp, 0.0_dp, KIND=dp)

      mask_re(:, :) = .FALSE.
      mask_im(:, :) = .FALSE.
      mask_re(1, 1) = .TRUE.  ! real and imaginary block both present
      mask_im(1, 1) = .TRUE.
      mask_im(1, 2) = .TRUE.  ! imaginary-only block
      mask_re(2, 1) = .TRUE.  ! real-only block
      ! block (2,2) is absent in both patterns

      CALL apply_dense_frechet(evals, z, dense_ref)
      CALL apply_sparse_block_frechet(evals, z, mask_re, mask_im, sparse_out)
      err = MAXVAL(ABS(sparse_out - dense_ref))

      IF (err > tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,ES12.5)') TRIM(label), "failed with error", err
      ELSE
         WRITE (*, '(A,1X,A,1X,ES12.5)') TRIM(label), "OK with error", err
      END IF

   END SUBROUTINE check_sparse_dbcsr_pattern

! **************************************************************************************************
!> \brief Check the production DBCSR Frechet path with real/imaginary-only blocks.
! **************************************************************************************************
   SUBROUTINE check_dbcsr_sparse_pattern(label, tol, nfail, mp_comm)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail
      TYPE(mp_comm_type), INTENT(IN)                     :: mp_comm

      INTEGER, PARAMETER                                 :: n = 4

      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: dbcsr_out, dense_ref, z
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_size, col_dist, row_blk_size, &
                                                            row_dist
      INTEGER, DIMENSION(:, :), POINTER                  :: pgrid
      LOGICAL                                            :: found
      REAL(KIND=dp)                                      :: err
      REAL(KIND=dp), DIMENSION(2, 2)                     :: block
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block_p
      REAL(KIND=dp), DIMENSION(n)                        :: evals
      TYPE(dbcsr_distribution_type)                      :: dist
      TYPE(dbcsr_type)                                   :: inner_im, inner_re, outer_im, outer_re

      evals = [0.7_dp, -0.2_dp, 1.4_dp, -1.1_dp]
      CALL setup_sparse_pattern_data(z)
      CALL apply_dense_frechet(evals, z, dense_ref)

      ALLOCATE (pgrid(1, 1), row_dist(2), col_dist(2), row_blk_size(2), col_blk_size(2))
      pgrid(:, :) = 0
      row_dist(:) = 0
      col_dist(:) = 0
      row_blk_size(:) = 2
      col_blk_size(:) = 2
      CALL dbcsr_distribution_new(dist, group=mp_comm%get_handle(), pgrid=pgrid, &
                                  row_dist=row_dist, col_dist=col_dist)

      CALL dbcsr_create(inner_re, "inner_re", dist, dbcsr_type_no_symmetry, row_blk_size, col_blk_size)
      CALL dbcsr_create(inner_im, "inner_im", dist, dbcsr_type_no_symmetry, row_blk_size, col_blk_size)
      CALL dbcsr_reserve_blocks(inner_re, rows=[1, 2], cols=[1, 1])
      CALL dbcsr_reserve_blocks(inner_im, rows=[1, 1], cols=[1, 2])
      CALL dbcsr_finalize(inner_re)
      CALL dbcsr_finalize(inner_im)

      block = REAL(z(1:2, 1:2), KIND=dp)
      CALL dbcsr_put_block(inner_re, 1, 1, block)
      block = AIMAG(z(1:2, 1:2))
      CALL dbcsr_put_block(inner_im, 1, 1, block)
      block = REAL(z(3:4, 1:2), KIND=dp)
      CALL dbcsr_put_block(inner_re, 2, 1, block)
      block = AIMAG(z(1:2, 3:4))
      CALL dbcsr_put_block(inner_im, 1, 2, block)

      CALL qs_ot_apply_complex_frechet_dbcsr(evals, inner_re, inner_im, outer_re, outer_im)

      dbcsr_out(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      CALL dbcsr_get_block_p(outer_re, 1, 1, block_p, found)
      IF (found) dbcsr_out(1:2, 1:2) = dbcsr_out(1:2, 1:2) + CMPLX(block_p, 0.0_dp, KIND=dp)
      CALL dbcsr_get_block_p(outer_im, 1, 1, block_p, found)
      IF (found) dbcsr_out(1:2, 1:2) = dbcsr_out(1:2, 1:2) + CMPLX(0.0_dp, block_p, KIND=dp)
      CALL dbcsr_get_block_p(outer_re, 1, 2, block_p, found)
      IF (found) dbcsr_out(1:2, 3:4) = dbcsr_out(1:2, 3:4) + CMPLX(block_p, 0.0_dp, KIND=dp)
      CALL dbcsr_get_block_p(outer_im, 1, 2, block_p, found)
      IF (found) dbcsr_out(1:2, 3:4) = dbcsr_out(1:2, 3:4) + CMPLX(0.0_dp, block_p, KIND=dp)
      CALL dbcsr_get_block_p(outer_re, 2, 1, block_p, found)
      IF (found) dbcsr_out(3:4, 1:2) = dbcsr_out(3:4, 1:2) + CMPLX(block_p, 0.0_dp, KIND=dp)
      CALL dbcsr_get_block_p(outer_im, 2, 1, block_p, found)
      IF (found) dbcsr_out(3:4, 1:2) = dbcsr_out(3:4, 1:2) + CMPLX(0.0_dp, block_p, KIND=dp)
      CALL dbcsr_get_block_p(outer_re, 2, 2, block_p, found)
      IF (found) dbcsr_out(3:4, 3:4) = dbcsr_out(3:4, 3:4) + CMPLX(block_p, 0.0_dp, KIND=dp)
      CALL dbcsr_get_block_p(outer_im, 2, 2, block_p, found)
      IF (found) dbcsr_out(3:4, 3:4) = dbcsr_out(3:4, 3:4) + CMPLX(0.0_dp, block_p, KIND=dp)

      err = MAXVAL(ABS(dbcsr_out - dense_ref))
      CALL dbcsr_release(outer_re)
      CALL dbcsr_release(outer_im)

      CALL qs_ot_apply_complex_frechet_dbcsr(evals, inner_re, inner_im, outer_re, outer_im, &
                                             adjoint=.TRUE.)
      CALL apply_dense_frechet(evals, z, dense_ref, adjoint=.TRUE.)

      dbcsr_out(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      CALL dbcsr_get_block_p(outer_re, 1, 1, block_p, found)
      IF (found) dbcsr_out(1:2, 1:2) = dbcsr_out(1:2, 1:2) + CMPLX(block_p, 0.0_dp, KIND=dp)
      CALL dbcsr_get_block_p(outer_im, 1, 1, block_p, found)
      IF (found) dbcsr_out(1:2, 1:2) = dbcsr_out(1:2, 1:2) + CMPLX(0.0_dp, block_p, KIND=dp)
      CALL dbcsr_get_block_p(outer_re, 1, 2, block_p, found)
      IF (found) dbcsr_out(1:2, 3:4) = dbcsr_out(1:2, 3:4) + CMPLX(block_p, 0.0_dp, KIND=dp)
      CALL dbcsr_get_block_p(outer_im, 1, 2, block_p, found)
      IF (found) dbcsr_out(1:2, 3:4) = dbcsr_out(1:2, 3:4) + CMPLX(0.0_dp, block_p, KIND=dp)
      CALL dbcsr_get_block_p(outer_re, 2, 1, block_p, found)
      IF (found) dbcsr_out(3:4, 1:2) = dbcsr_out(3:4, 1:2) + CMPLX(block_p, 0.0_dp, KIND=dp)
      CALL dbcsr_get_block_p(outer_im, 2, 1, block_p, found)
      IF (found) dbcsr_out(3:4, 1:2) = dbcsr_out(3:4, 1:2) + CMPLX(0.0_dp, block_p, KIND=dp)
      CALL dbcsr_get_block_p(outer_re, 2, 2, block_p, found)
      IF (found) dbcsr_out(3:4, 3:4) = dbcsr_out(3:4, 3:4) + CMPLX(block_p, 0.0_dp, KIND=dp)
      CALL dbcsr_get_block_p(outer_im, 2, 2, block_p, found)
      IF (found) dbcsr_out(3:4, 3:4) = dbcsr_out(3:4, 3:4) + CMPLX(0.0_dp, block_p, KIND=dp)

      err = MAX(err, MAXVAL(ABS(dbcsr_out - dense_ref)))
      IF (err > tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,ES12.5)') TRIM(label), "failed with error", err
      ELSE
         WRITE (*, '(A,1X,A,1X,ES12.5)') TRIM(label), "OK with error", err
      END IF

      CALL dbcsr_release(inner_re)
      CALL dbcsr_release(inner_im)
      CALL dbcsr_release(outer_re)
      CALL dbcsr_release(outer_im)
      CALL dbcsr_distribution_release(dist)
      DEALLOCATE (pgrid, row_dist, col_dist, row_blk_size, col_blk_size)

   END SUBROUTINE check_dbcsr_sparse_pattern

! **************************************************************************************************
!> \brief Fill the dense companion to the sparse Frechet pattern tests.
! **************************************************************************************************
   SUBROUTINE setup_sparse_pattern_data(z)
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(OUT)     :: z

      z(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      z(1, 1) = CMPLX(0.3_dp, -0.4_dp, KIND=dp)
      z(1, 2) = CMPLX(-0.2_dp, 0.1_dp, KIND=dp)
      z(2, 1) = CMPLX(0.5_dp, 0.8_dp, KIND=dp)
      z(2, 2) = CMPLX(-0.7_dp, 0.2_dp, KIND=dp)
      z(1, 3) = CMPLX(0.0_dp, 0.6_dp, KIND=dp)
      z(1, 4) = CMPLX(0.0_dp, -0.3_dp, KIND=dp)
      z(2, 3) = CMPLX(0.0_dp, 0.4_dp, KIND=dp)
      z(2, 4) = CMPLX(0.0_dp, 0.2_dp, KIND=dp)
      z(3, 1) = CMPLX(-0.9_dp, 0.0_dp, KIND=dp)
      z(3, 2) = CMPLX(0.4_dp, 0.0_dp, KIND=dp)
      z(4, 1) = CMPLX(0.2_dp, 0.0_dp, KIND=dp)
      z(4, 2) = CMPLX(-0.1_dp, 0.0_dp, KIND=dp)

   END SUBROUTINE setup_sparse_pattern_data

! **************************************************************************************************
!> \brief Dense reference application of the complex exponential Frechet kernel.
! **************************************************************************************************
   SUBROUTINE apply_dense_frechet(evals, z, out, adjoint)
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: evals
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: z
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(OUT)     :: out
      LOGICAL, INTENT(IN), OPTIONAL                      :: adjoint

      COMPLEX(KIND=dp)                                   :: kernel
      INTEGER                                            :: i, j
      LOGICAL                                            :: use_adjoint

      use_adjoint = .FALSE.
      IF (PRESENT(adjoint)) use_adjoint = adjoint

      DO i = 1, SIZE(z, 1)
      DO j = 1, SIZE(z, 2)
         kernel = qs_ot_complex_exp_frechet_kernel(evals(i), evals(j))
         IF (use_adjoint) kernel = CONJG(kernel)
         out(i, j) = z(i, j)*kernel
      END DO
      END DO

   END SUBROUTINE apply_dense_frechet

! **************************************************************************************************
!> \brief Sparse block-pattern application equivalent to the DBCSR block iterator logic.
! **************************************************************************************************
   SUBROUTINE apply_sparse_block_frechet(evals, z, mask_re, mask_im, out)
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: evals
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: z
      LOGICAL, DIMENSION(:, :), INTENT(IN)               :: mask_re, mask_im
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(OUT)     :: out

      COMPLEX(KIND=dp)                                   :: cval
      INTEGER                                            :: bc, br, i, i1, i2, j, j1, j2

      out(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      DO br = 1, SIZE(mask_re, 1)
      DO bc = 1, SIZE(mask_re, 2)
         IF (.NOT. (mask_re(br, bc) .OR. mask_im(br, bc))) CYCLE
         i1 = (br - 1)*2 + 1
         i2 = br*2
         j1 = (bc - 1)*2 + 1
         j2 = bc*2
         DO i = i1, i2
         DO j = j1, j2
            cval = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
            IF (mask_re(br, bc)) cval = cval + CMPLX(REAL(z(i, j), KIND=dp), 0.0_dp, KIND=dp)
            IF (mask_im(br, bc)) cval = cval + CMPLX(0.0_dp, AIMAG(z(i, j)), KIND=dp)
            out(i, j) = cval*qs_ot_complex_exp_frechet_kernel(evals(i), evals(j))
         END DO
         END DO
      END DO
      END DO

   END SUBROUTINE apply_sparse_block_frechet

! **************************************************************************************************
!> \brief Check the weighted fixed-N Mermin gradient used for k-point OT energy variables.
! **************************************************************************************************
   SUBROUTINE check_weighted_kpoint_mermin_gradient(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: nkp = 2, nmo = 4
      REAL(KIND=dp), PARAMETER                           :: fd_step = 1.0E-6_dp, maxocc = 2.0_dp, &
                                                            nel = 5.0_dp, sigma = 0.025_dp

      INTEGER                                            :: ik, imo
      REAL(KIND=dp)                                      :: err, fd, ref
      REAL(KIND=dp), DIMENSION(nmo, nkp)                 :: eps, eps_m, eps_p, grad, h
      REAL(KIND=dp), DIMENSION(nkp)                      :: wk

      wk = [0.25_dp, 0.75_dp]
      eps(:, 1) = [-0.20_dp, -0.060_dp, -0.005_dp, 0.080_dp]
      eps(:, 2) = [-0.16_dp, -0.040_dp, 0.015_dp, 0.110_dp]
      h(:, 1) = [-0.19_dp, -0.045_dp, 0.020_dp, 0.100_dp]
      h(:, 2) = [-0.15_dp, -0.055_dp, 0.035_dp, 0.090_dp]

      CALL dense_kpoint_mermin_gradient(eps, h, wk, nel, sigma, maxocc, grad)

      err = 0.0_dp
      DO ik = 1, nkp
         DO imo = 1, nmo
            eps_p(:, :) = eps(:, :)
            eps_m(:, :) = eps(:, :)
            eps_p(imo, ik) = eps_p(imo, ik) + fd_step
            eps_m(imo, ik) = eps_m(imo, ik) - fd_step
            fd = (dense_kpoint_mermin_value(eps_p, h, wk, nel, sigma, maxocc) - &
                  dense_kpoint_mermin_value(eps_m, h, wk, nel, sigma, maxocc))/(2.0_dp*fd_step)
            ref = grad(imo, ik)
            err = MAX(err, ABS(fd - ref))
         END DO
      END DO

      IF (err > tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,ES12.5)') TRIM(label), "failed with error", err
      ELSE
         WRITE (*, '(A,1X,A,1X,ES12.5)') TRIM(label), "OK with error", err
      END IF

   END SUBROUTINE check_weighted_kpoint_mermin_gradient

! **************************************************************************************************
!> \brief Check the exact auxiliary-energy projection for the fixed-N Mermin subproblem.
! **************************************************************************************************
   SUBROUTINE check_kpoint_mermin_energy_projection(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: nkp = 2, nmo = 4
      REAL(KIND=dp), PARAMETER                           :: maxocc = 2.0_dp, nel = 4.35_dp, &
                                                            sigma = 0.030_dp

      REAL(KIND=dp)                                      :: err, gauge_shift
      REAL(KIND=dp), DIMENSION(nmo, nkp)                 :: eps, grad, h, projected
      REAL(KIND=dp), DIMENSION(nkp)                      :: wk

      wk = [0.35_dp, 0.65_dp]
      eps(:, 1) = [-0.34_dp, -0.07_dp, 0.025_dp, 0.18_dp]
      eps(:, 2) = [-0.29_dp, -0.02_dp, 0.075_dp, 0.22_dp]
      h(:, 1) = [-0.31_dp, -0.035_dp, 0.045_dp, 0.16_dp]
      h(:, 2) = [-0.27_dp, -0.055_dp, 0.055_dp, 0.20_dp]

      gauge_shift = SUM(eps(:, :) - h(:, :))/REAL(nkp*nmo, KIND=dp)
      projected(:, :) = h(:, :) + gauge_shift
      CALL dense_kpoint_mermin_gradient(projected, h, wk, nel, sigma, maxocc, grad)
      err = MAXVAL(ABS(grad))

      IF (err > tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,ES12.5)') TRIM(label), "failed with error", err
      ELSE
         WRITE (*, '(A,1X,A,1X,ES12.5)') TRIM(label), "OK with error", err
      END IF

   END SUBROUTINE check_kpoint_mermin_energy_projection

! **************************************************************************************************
!> \brief Check that the diagonal energy preconditioner removes the fixed-N gauge
!>        in the occupation-curvature metric, not by an arithmetic band average.
! **************************************************************************************************
   SUBROUTINE check_fixed_n_preconditioner_gauge(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 5
      REAL(KIND=dp), PARAMETER                           :: curvature_floor = 1.0_dp

      REAL(KIND=dp)                                      :: err, shift, wrong_err
      REAL(KIND=dp), DIMENSION(n)                        :: curv, grad, step, wrong_step

      curv = [0.15_dp, 0.8_dp, 2.4_dp, 6.0_dp, 11.0_dp]
      grad = [0.32_dp, -0.71_dp, 0.14_dp, 0.55_dp, -0.28_dp]

      step(:) = grad(:)/MAX(curv(:), curvature_floor)
      shift = SUM(curv(:)*step(:))/SUM(curv(:))
      step(:) = step(:) - shift
      err = ABS(SUM(curv(:)*step(:)))

      wrong_step(:) = grad(:)/MAX(curv(:), curvature_floor)
      wrong_step(:) = wrong_step(:) - SUM(wrong_step(:))/REAL(n, KIND=dp)
      wrong_err = ABS(SUM(curv(:)*wrong_step(:)))

      IF (err > tol .OR. wrong_err <= 1.0E-3_dp) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,2(1X,ES12.5))') TRIM(label), "failed with errors", err, wrong_err
      ELSE
         WRITE (*, '(A,1X,A,2(1X,ES12.5))') TRIM(label), "OK with errors", err, wrong_err
      END IF

   END SUBROUTINE check_fixed_n_preconditioner_gauge

! **************************************************************************************************
!> \brief Dense weighted k-point Mermin free-energy model.
! **************************************************************************************************
   FUNCTION dense_kpoint_mermin_value(eps, h, wk, nel, sigma, maxocc) RESULT(value)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: eps, h
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: wk
      REAL(KIND=dp), INTENT(IN)                          :: nel, sigma, maxocc
      REAL(KIND=dp)                                      :: value

      INTEGER                                            :: ik
      REAL(KIND=dp)                                      :: kTS, mu
      REAL(KIND=dp), DIMENSION(SIZE(eps, 1), SIZE(eps, 2)) :: occ

      CALL dense_fermi_kpoint(occ, mu, kTS, eps, nel, wk, sigma, maxocc)
      value = kTS
      DO ik = 1, SIZE(eps, 2)
         value = value + wk(ik)*SUM(occ(:, ik)*h(:, ik))
      END DO

   END FUNCTION dense_kpoint_mermin_value

! **************************************************************************************************
!> \brief Dense weighted fixed-N Mermin gradient.
! **************************************************************************************************
   SUBROUTINE dense_kpoint_mermin_gradient(eps, h, wk, nel, sigma, maxocc, grad)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: eps, h
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: wk
      REAL(KIND=dp), INTENT(IN)                          :: nel, sigma, maxocc
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: grad

      INTEGER                                            :: ik, imo
      REAL(KIND=dp)                                      :: gdotv, gsum, gvalue, mu, kTS, v
      REAL(KIND=dp), DIMENSION(SIZE(eps, 1), SIZE(eps, 2)) :: occ

      CALL dense_fermi_kpoint(occ, mu, kTS, eps, nel, wk, sigma, maxocc)

      gsum = 0.0_dp
      gdotv = 0.0_dp
      DO ik = 1, SIZE(eps, 2)
         DO imo = 1, SIZE(eps, 1)
            gvalue = maxocc*(occ(imo, ik)/maxocc)*(1.0_dp - occ(imo, ik)/maxocc)/sigma
            v = h(imo, ik) - eps(imo, ik) + mu
            gsum = gsum + wk(ik)*gvalue
            gdotv = gdotv + wk(ik)*gvalue*v
         END DO
      END DO

      grad(:, :) = 0.0_dp
      IF (ABS(gsum) > EPSILON(gsum)) THEN
         DO ik = 1, SIZE(eps, 2)
            DO imo = 1, SIZE(eps, 1)
               gvalue = maxocc*(occ(imo, ik)/maxocc)*(1.0_dp - occ(imo, ik)/maxocc)/sigma
               v = h(imo, ik) - eps(imo, ik) + mu
               grad(imo, ik) = wk(ik)*gvalue*(gdotv/gsum - v)
            END DO
         END DO
      END IF

   END SUBROUTINE dense_kpoint_mermin_gradient

! **************************************************************************************************
!> \brief Check the coupled complex rotation and auxiliary-energy Mermin gradient.
! **************************************************************************************************
   SUBROUTINE check_complex_kpoint_mermin_coupled_gradient(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 3, nkp = 2
      REAL(KIND=dp), PARAMETER                           :: fd_step = 1.0E-6_dp, maxocc = 2.0_dp, &
                                                            nel = 3.20_dp, sigma = 0.035_dp

      COMPLEX(KIND=dp), DIMENSION(n, n, nkp)             :: base
      INTEGER                                            :: i, ik, j
      REAL(KIND=dp)                                      :: err, fd, value_m, value_p
      REAL(KIND=dp), DIMENSION(n, nkp)                   :: eps, eps_m, eps_p, grad_eps
      REAL(KIND=dp), DIMENSION(nkp)                      :: wk
      REAL(KIND=dp), DIMENSION(n, n, nkp)                :: grad_im, grad_re, x_im, x_m_im, &
                                                            x_m_re, x_p_im, x_p_re, x_re

      wk = [0.40_dp, 0.60_dp]
      eps(:, 1) = [-0.23_dp, -0.035_dp, 0.090_dp]
      eps(:, 2) = [-0.19_dp, -0.010_dp, 0.125_dp]

      base(:, :, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      base(1, 1, 1) = CMPLX(-0.26_dp, 0.0_dp, KIND=dp)
      base(2, 2, 1) = CMPLX(-0.015_dp, 0.0_dp, KIND=dp)
      base(3, 3, 1) = CMPLX(0.14_dp, 0.0_dp, KIND=dp)
      base(1, 2, 1) = CMPLX(0.045_dp, -0.025_dp, KIND=dp)
      base(1, 3, 1) = CMPLX(-0.030_dp, 0.035_dp, KIND=dp)
      base(2, 3, 1) = CMPLX(0.020_dp, 0.030_dp, KIND=dp)
      base(2, 1, 1) = CONJG(base(1, 2, 1))
      base(3, 1, 1) = CONJG(base(1, 3, 1))
      base(3, 2, 1) = CONJG(base(2, 3, 1))

      base(1, 1, 2) = CMPLX(-0.22_dp, 0.0_dp, KIND=dp)
      base(2, 2, 2) = CMPLX(0.005_dp, 0.0_dp, KIND=dp)
      base(3, 3, 2) = CMPLX(0.18_dp, 0.0_dp, KIND=dp)
      base(1, 2, 2) = CMPLX(-0.035_dp, 0.015_dp, KIND=dp)
      base(1, 3, 2) = CMPLX(0.025_dp, 0.020_dp, KIND=dp)
      base(2, 3, 2) = CMPLX(0.030_dp, -0.025_dp, KIND=dp)
      base(2, 1, 2) = CONJG(base(1, 2, 2))
      base(3, 1, 2) = CONJG(base(1, 3, 2))
      base(3, 2, 2) = CONJG(base(2, 3, 2))

      x_re(:, :, :) = 0.0_dp
      x_re(1, 2, 1) = 0.055_dp
      x_re(1, 3, 1) = -0.025_dp
      x_re(2, 3, 1) = 0.035_dp
      x_re(1, 2, 2) = -0.040_dp
      x_re(1, 3, 2) = 0.030_dp
      x_re(2, 3, 2) = 0.020_dp
      DO ik = 1, nkp
         DO j = 2, n
            DO i = 1, j - 1
               x_re(j, i, ik) = -x_re(i, j, ik)
            END DO
         END DO
      END DO

      x_im(:, :, :) = 0.0_dp
      x_im(1, 1, 1) = 0.020_dp
      x_im(2, 2, 1) = -0.035_dp
      x_im(3, 3, 1) = 0.015_dp
      x_im(1, 2, 1) = -0.018_dp
      x_im(1, 3, 1) = 0.012_dp
      x_im(2, 3, 1) = 0.026_dp
      x_im(1, 1, 2) = -0.015_dp
      x_im(2, 2, 2) = 0.030_dp
      x_im(3, 3, 2) = -0.010_dp
      x_im(1, 2, 2) = 0.014_dp
      x_im(1, 3, 2) = -0.022_dp
      x_im(2, 3, 2) = 0.018_dp
      DO ik = 1, nkp
         DO j = 2, n
            DO i = 1, j - 1
               x_im(j, i, ik) = x_im(i, j, ik)
            END DO
         END DO
      END DO

      CALL dense_complex_kpoint_mermin_gradient(x_re, x_im, base, eps, wk, nel, sigma, maxocc, &
                                                grad_re, grad_im, grad_eps, value_p)

      err = 0.0_dp
      DO ik = 1, nkp
         DO i = 1, n
            eps_p(:, :) = eps(:, :)
            eps_m(:, :) = eps(:, :)
            eps_p(i, ik) = eps_p(i, ik) + fd_step
            eps_m(i, ik) = eps_m(i, ik) - fd_step
            fd = (dense_complex_kpoint_mermin_value(x_re, x_im, base, eps_p, wk, nel, sigma, maxocc) - &
                  dense_complex_kpoint_mermin_value(x_re, x_im, base, eps_m, wk, nel, sigma, maxocc))/ &
                 (2.0_dp*fd_step)
            err = MAX(err, ABS(fd - grad_eps(i, ik)))
         END DO
      END DO

      DO ik = 1, nkp
         DO j = 2, n
            DO i = 1, j - 1
               x_p_re(:, :, :) = x_re(:, :, :)
               x_m_re(:, :, :) = x_re(:, :, :)
               x_p_im(:, :, :) = x_im(:, :, :)
               x_m_im(:, :, :) = x_im(:, :, :)
               x_p_re(i, j, ik) = x_p_re(i, j, ik) + fd_step
               x_p_re(j, i, ik) = x_p_re(j, i, ik) - fd_step
               x_m_re(i, j, ik) = x_m_re(i, j, ik) - fd_step
               x_m_re(j, i, ik) = x_m_re(j, i, ik) + fd_step
               value_p = dense_complex_kpoint_mermin_value(x_p_re, x_p_im, base, eps, wk, nel, &
                                                           sigma, maxocc)
               value_m = dense_complex_kpoint_mermin_value(x_m_re, x_m_im, base, eps, wk, nel, &
                                                           sigma, maxocc)
               fd = (value_p - value_m)/(2.0_dp*fd_step)
               err = MAX(err, ABS(fd - grad_re(i, j, ik)))

               x_p_re(:, :, :) = x_re(:, :, :)
               x_m_re(:, :, :) = x_re(:, :, :)
               x_p_im(:, :, :) = x_im(:, :, :)
               x_m_im(:, :, :) = x_im(:, :, :)
               x_p_im(i, j, ik) = x_p_im(i, j, ik) + fd_step
               x_p_im(j, i, ik) = x_p_im(j, i, ik) + fd_step
               x_m_im(i, j, ik) = x_m_im(i, j, ik) - fd_step
               x_m_im(j, i, ik) = x_m_im(j, i, ik) - fd_step
               value_p = dense_complex_kpoint_mermin_value(x_p_re, x_p_im, base, eps, wk, nel, &
                                                           sigma, maxocc)
               value_m = dense_complex_kpoint_mermin_value(x_m_re, x_m_im, base, eps, wk, nel, &
                                                           sigma, maxocc)
               fd = (value_p - value_m)/(2.0_dp*fd_step)
               err = MAX(err, ABS(fd - grad_im(i, j, ik)))
            END DO
         END DO

         DO i = 1, n
            x_p_re(:, :, :) = x_re(:, :, :)
            x_m_re(:, :, :) = x_re(:, :, :)
            x_p_im(:, :, :) = x_im(:, :, :)
            x_m_im(:, :, :) = x_im(:, :, :)
            x_p_im(i, i, ik) = x_p_im(i, i, ik) + fd_step
            x_m_im(i, i, ik) = x_m_im(i, i, ik) - fd_step
            value_p = dense_complex_kpoint_mermin_value(x_p_re, x_p_im, base, eps, wk, nel, &
                                                        sigma, maxocc)
            value_m = dense_complex_kpoint_mermin_value(x_m_re, x_m_im, base, eps, wk, nel, &
                                                        sigma, maxocc)
            fd = (value_p - value_m)/(2.0_dp*fd_step)
            err = MAX(err, ABS(fd - 0.5_dp*grad_im(i, i, ik)))
         END DO
      END DO

      IF (err > tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,ES12.5)') TRIM(label), "failed with error", err
      ELSE
         WRITE (*, '(A,1X,A,1X,ES12.5)') TRIM(label), "OK with error", err
      END IF

   END SUBROUTINE check_complex_kpoint_mermin_coupled_gradient

! **************************************************************************************************
!> \brief Check the reduced fixed-N Schur Hessian against finite differences.
!>        The reduced problem uses eps=h(x), so the dense formula must include
!>        the fixed-occupation rotation Hessian and the fixed-N energy response
!>        projector, H_red = H_xx - R^T P_N R.
! **************************************************************************************************
   SUBROUTINE check_complex_fixed_n_schur_hessian(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 3, nkp = 1, nrot = n*(n - 1)
      REAL(KIND=dp), PARAMETER                           :: fd_step = 1.0E-5_dp, maxocc = 2.0_dp, &
                                                            nel = 3.15_dp, sigma = 0.045_dp

      COMPLEX(KIND=dp), DIMENSION(n, n, nkp)             :: base
      INTEGER                                            :: i, ivar, j, r, s
      INTEGER, DIMENSION(nrot)                           :: var_i, var_j, var_kind
      REAL(KIND=dp)                                      :: dense_term, err, kTS, mu, projected_term, &
                                                            rtd_r, rtd_s, value
      REAL(KIND=dp), DIMENSION(n, nkp)                   :: eps_base, eps_m, eps_p, grad_eps, h_base, &
                                                            h_m, h_p, occ
      REAL(KIND=dp), DIMENSION(nkp)                      :: wk
      REAL(KIND=dp), DIMENSION(n)                        :: curv, response_r, response_s
      REAL(KIND=dp), DIMENSION(nrot)                     :: grad_m_vec, grad_p_vec
      REAL(KIND=dp), DIMENSION(nrot, nrot)               :: hfixed, hred, schur
      REAL(KIND=dp), DIMENSION(n, n, nkp)                :: grad_im, grad_re, x_im, x_m_im, &
                                                            x_m_re, x_p_im, x_p_re, x_re

      wk = [1.0_dp]
      base(:, :, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      base(1, 1, 1) = CMPLX(-0.24_dp, 0.0_dp, KIND=dp)
      base(2, 2, 1) = CMPLX(0.01_dp, 0.0_dp, KIND=dp)
      base(3, 3, 1) = CMPLX(0.19_dp, 0.0_dp, KIND=dp)
      base(1, 2, 1) = CMPLX(0.045_dp, -0.025_dp, KIND=dp)
      base(1, 3, 1) = CMPLX(-0.025_dp, 0.030_dp, KIND=dp)
      base(2, 3, 1) = CMPLX(0.035_dp, 0.020_dp, KIND=dp)
      base(2, 1, 1) = CONJG(base(1, 2, 1))
      base(3, 1, 1) = CONJG(base(1, 3, 1))
      base(3, 2, 1) = CONJG(base(2, 3, 1))

      x_re(:, :, :) = 0.0_dp
      x_re(1, 2, 1) = 0.050_dp
      x_re(1, 3, 1) = -0.030_dp
      x_re(2, 3, 1) = 0.025_dp
      DO j = 2, n
         DO i = 1, j - 1
            x_re(j, i, 1) = -x_re(i, j, 1)
         END DO
      END DO

      x_im(:, :, :) = 0.0_dp
      x_im(1, 1, 1) = 0.015_dp
      x_im(2, 2, 1) = -0.025_dp
      x_im(3, 3, 1) = 0.010_dp
      x_im(1, 2, 1) = -0.020_dp
      x_im(1, 3, 1) = 0.015_dp
      x_im(2, 3, 1) = 0.018_dp
      DO j = 2, n
         DO i = 1, j - 1
            x_im(j, i, 1) = x_im(i, j, 1)
         END DO
      END DO

      ivar = 0
      DO j = 2, n
         DO i = 1, j - 1
            ivar = ivar + 1
            var_i(ivar) = i
            var_j(ivar) = j
            var_kind(ivar) = 1
            ivar = ivar + 1
            var_i(ivar) = i
            var_j(ivar) = j
            var_kind(ivar) = 2
         END DO
      END DO

      CALL dense_complex_kpoint_hdiag(x_re, x_im, base, h_base)
      eps_base(:, :) = h_base(:, :)
      CALL dense_fermi_kpoint(occ, mu, kTS, eps_base, nel, wk, sigma, maxocc)
      DO i = 1, n
         curv(i) = wk(1)*maxocc*(occ(i, 1)/maxocc)*(1.0_dp - occ(i, 1)/maxocc)/sigma
      END DO

      DO s = 1, nrot
         x_p_re(:, :, :) = x_re(:, :, :)
         x_m_re(:, :, :) = x_re(:, :, :)
         x_p_im(:, :, :) = x_im(:, :, :)
         x_m_im(:, :, :) = x_im(:, :, :)
         CALL perturb_complex_rotation_var(x_p_re, x_p_im, var_i(s), var_j(s), var_kind(s), fd_step)
         CALL perturb_complex_rotation_var(x_m_re, x_m_im, var_i(s), var_j(s), var_kind(s), -fd_step)

         CALL dense_complex_kpoint_mermin_gradient(x_p_re, x_p_im, base, eps_base, wk, nel, &
                                                   sigma, maxocc, grad_re, grad_im, grad_eps, value)
         CALL extract_rotation_gradient_vector(grad_re, grad_im, var_i, var_j, var_kind, grad_p_vec)
         CALL dense_complex_kpoint_mermin_gradient(x_m_re, x_m_im, base, eps_base, wk, nel, &
                                                   sigma, maxocc, grad_re, grad_im, grad_eps, value)
         CALL extract_rotation_gradient_vector(grad_re, grad_im, var_i, var_j, var_kind, grad_m_vec)
         hfixed(:, s) = (grad_p_vec(:) - grad_m_vec(:))/(2.0_dp*fd_step)

         CALL dense_complex_kpoint_hdiag(x_p_re, x_p_im, base, h_p)
         CALL dense_complex_kpoint_hdiag(x_m_re, x_m_im, base, h_m)
         eps_p(:, :) = h_p(:, :)
         eps_m(:, :) = h_m(:, :)
         response_s(:) = (h_p(:, 1) - h_m(:, 1))/(2.0_dp*fd_step)

         CALL dense_complex_kpoint_mermin_gradient(x_p_re, x_p_im, base, eps_p, wk, nel, &
                                                   sigma, maxocc, grad_re, grad_im, grad_eps, value)
         CALL extract_rotation_gradient_vector(grad_re, grad_im, var_i, var_j, var_kind, grad_p_vec)
         CALL dense_complex_kpoint_mermin_gradient(x_m_re, x_m_im, base, eps_m, wk, nel, &
                                                   sigma, maxocc, grad_re, grad_im, grad_eps, value)
         CALL extract_rotation_gradient_vector(grad_re, grad_im, var_i, var_j, var_kind, grad_m_vec)
         hred(:, s) = (grad_p_vec(:) - grad_m_vec(:))/(2.0_dp*fd_step)

         DO r = 1, nrot
            response_r(:) = 0.0_dp
            x_p_re(:, :, :) = x_re(:, :, :)
            x_m_re(:, :, :) = x_re(:, :, :)
            x_p_im(:, :, :) = x_im(:, :, :)
            x_m_im(:, :, :) = x_im(:, :, :)
            CALL perturb_complex_rotation_var(x_p_re, x_p_im, var_i(r), var_j(r), var_kind(r), fd_step)
            CALL perturb_complex_rotation_var(x_m_re, x_m_im, var_i(r), var_j(r), var_kind(r), -fd_step)
            CALL dense_complex_kpoint_hdiag(x_p_re, x_p_im, base, h_p)
            CALL dense_complex_kpoint_hdiag(x_m_re, x_m_im, base, h_m)
            response_r(:) = (h_p(:, 1) - h_m(:, 1))/(2.0_dp*fd_step)
            dense_term = DOT_PRODUCT(curv(:)*response_r(:), response_s(:))
            rtd_r = DOT_PRODUCT(curv(:), response_r(:))
            rtd_s = DOT_PRODUCT(curv(:), response_s(:))
            projected_term = dense_term - rtd_r*rtd_s/SUM(curv(:))
            schur(r, s) = hfixed(r, s) - projected_term
         END DO
      END DO

      err = MAXVAL(ABS(schur(:, :) - hred(:, :)))

      IF (err > tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,ES12.5)') TRIM(label), "failed with error", err
      ELSE
         WRITE (*, '(A,1X,A,1X,ES12.5)') TRIM(label), "OK with error", err
      END IF

   END SUBROUTINE check_complex_fixed_n_schur_hessian

! **************************************************************************************************
!> \brief Check the Schur response channel against finite fixed-N Mermin occupation response.
!>        The cross block is d grad_eps / d theta = -P_N R.  The production
!>        preconditioner uses the projected pseudo-inverse R = -P_N^+ d grad_eps/d theta,
!>        so null occupation-curvature modes cannot drive the rotation RHS.
! **************************************************************************************************
   SUBROUTINE check_complex_fixed_n_response_channel(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 3, nkp = 1, nrot = n*(n - 1)
      REAL(KIND=dp), PARAMETER                           :: fd_step = 1.0E-6_dp, maxocc = 2.0_dp, &
                                                            nel = 3.05_dp, sigma = 0.080_dp

      COMPLEX(KIND=dp), DIMENSION(n, n, nkp)             :: base
      INTEGER                                            :: i, ivar, j, s
      INTEGER, DIMENSION(nrot)                           :: var_i, var_j, var_kind
      REAL(KIND=dp)                                      :: curv_sum, err, kTS, mu, response_shift
      REAL(KIND=dp), DIMENSION(n, nkp)                   :: eps_base, grad_eps_m, grad_eps_p, &
                                                            h_base, h_m, h_p, occ
      REAL(KIND=dp), DIMENSION(nkp)                      :: wk
      REAL(KIND=dp), DIMENSION(n)                        :: cross_fd, curv, dh, projected_response, &
                                                            ref_cross, response_from_cross
      REAL(KIND=dp), DIMENSION(n, n, nkp)                :: x_im, x_m_im, x_m_re, x_p_im, x_p_re, x_re

      wk = [1.0_dp]
      base(:, :, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      base(1, 1, 1) = CMPLX(-0.20_dp, 0.0_dp, KIND=dp)
      base(2, 2, 1) = CMPLX(0.02_dp, 0.0_dp, KIND=dp)
      base(3, 3, 1) = CMPLX(0.16_dp, 0.0_dp, KIND=dp)
      base(1, 2, 1) = CMPLX(0.050_dp, -0.030_dp, KIND=dp)
      base(1, 3, 1) = CMPLX(-0.020_dp, 0.025_dp, KIND=dp)
      base(2, 3, 1) = CMPLX(0.030_dp, 0.035_dp, KIND=dp)
      base(2, 1, 1) = CONJG(base(1, 2, 1))
      base(3, 1, 1) = CONJG(base(1, 3, 1))
      base(3, 2, 1) = CONJG(base(2, 3, 1))

      x_re(:, :, :) = 0.0_dp
      x_re(1, 2, 1) = 0.045_dp
      x_re(1, 3, 1) = -0.020_dp
      x_re(2, 3, 1) = 0.030_dp
      DO j = 2, n
         DO i = 1, j - 1
            x_re(j, i, 1) = -x_re(i, j, 1)
         END DO
      END DO

      x_im(:, :, :) = 0.0_dp
      x_im(1, 1, 1) = 0.010_dp
      x_im(2, 2, 1) = -0.020_dp
      x_im(3, 3, 1) = 0.018_dp
      x_im(1, 2, 1) = -0.015_dp
      x_im(1, 3, 1) = 0.020_dp
      x_im(2, 3, 1) = 0.012_dp
      DO j = 2, n
         DO i = 1, j - 1
            x_im(j, i, 1) = x_im(i, j, 1)
         END DO
      END DO

      ivar = 0
      DO j = 2, n
         DO i = 1, j - 1
            ivar = ivar + 1
            var_i(ivar) = i
            var_j(ivar) = j
            var_kind(ivar) = 1
            ivar = ivar + 1
            var_i(ivar) = i
            var_j(ivar) = j
            var_kind(ivar) = 2
         END DO
      END DO

      CALL dense_complex_kpoint_hdiag(x_re, x_im, base, h_base)
      eps_base(:, 1) = h_base(:, 1) + [0.018_dp, -0.012_dp, 0.006_dp]
      CALL dense_fermi_kpoint(occ, mu, kTS, eps_base, nel, wk, sigma, maxocc)
      DO i = 1, n
         curv(i) = wk(1)*maxocc*(occ(i, 1)/maxocc)*(1.0_dp - occ(i, 1)/maxocc)/sigma
      END DO
      curv_sum = SUM(curv(:))

      err = 0.0_dp
      DO s = 1, nrot
         x_p_re(:, :, :) = x_re(:, :, :)
         x_m_re(:, :, :) = x_re(:, :, :)
         x_p_im(:, :, :) = x_im(:, :, :)
         x_m_im(:, :, :) = x_im(:, :, :)
         CALL perturb_complex_rotation_var(x_p_re, x_p_im, var_i(s), var_j(s), var_kind(s), fd_step)
         CALL perturb_complex_rotation_var(x_m_re, x_m_im, var_i(s), var_j(s), var_kind(s), -fd_step)

         CALL dense_complex_kpoint_hdiag(x_p_re, x_p_im, base, h_p)
         CALL dense_complex_kpoint_hdiag(x_m_re, x_m_im, base, h_m)
         CALL dense_kpoint_mermin_gradient(eps_base, h_p, wk, nel, sigma, maxocc, grad_eps_p)
         CALL dense_kpoint_mermin_gradient(eps_base, h_m, wk, nel, sigma, maxocc, grad_eps_m)

         dh(:) = (h_p(:, 1) - h_m(:, 1))/(2.0_dp*fd_step)
         response_shift = SUM(curv(:)*dh(:))/curv_sum
         projected_response(:) = dh(:) - response_shift
         ref_cross(:) = curv(:)*(response_shift - dh(:))
         cross_fd(:) = (grad_eps_p(:, 1) - grad_eps_m(:, 1))/(2.0_dp*fd_step)

         response_from_cross(:) = 0.0_dp
         DO i = 1, n
            IF (curv(i) > EPSILON(curv(i))) response_from_cross(i) = -cross_fd(i)/curv(i)
         END DO
         response_shift = SUM(curv(:)*response_from_cross(:))/curv_sum
         response_from_cross(:) = response_from_cross(:) - response_shift

         err = MAX(err, MAXVAL(ABS(cross_fd(:) - ref_cross(:))))
         err = MAX(err, MAXVAL(ABS(response_from_cross(:) - projected_response(:))))
      END DO

      IF (err > tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,ES12.5)') TRIM(label), "failed with error", err
      ELSE
         WRITE (*, '(A,1X,A,1X,ES12.5)') TRIM(label), "OK with error", err
      END IF

   END SUBROUTINE check_complex_fixed_n_response_channel

! **************************************************************************************************
!> \brief Check the response blend used by the finite complex split Schur model.
!>        The occupation-projector response is the derivative of the rotation
!>        gradient with respect to the band occupation weights.  For a purely
!>        diagonal projector this reduces to d[U^H H U]_ii/dtheta; with a
!>        spectral projector blend it also contains the finite complex projector
!>        rotation response.
! **************************************************************************************************
   SUBROUTINE check_complex_projector_blended_response_channel(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 3, nkp = 1, nrot = n*(n - 1)
      REAL(KIND=dp), PARAMETER                           :: active_tol = 1.0E-4_dp, &
                                                            fd_occ_step = 1.0E-6_dp, &
                                                            fd_step = 1.0E-6_dp, &
                                                            projector_weight = 0.35_dp

      COMPLEX(KIND=dp), DIMENSION(n, n, nkp)             :: base
      INTEGER                                            :: i, ivar, j, s
      INTEGER, DIMENSION(nrot)                           :: var_i, var_j, var_kind
      REAL(KIND=dp)                                      :: active_err, err
      REAL(KIND=dp), DIMENSION(n, nkp)                   :: h_m, h_p
      REAL(KIND=dp), DIMENSION(n)                        :: occ_scale, occ_scale_m, occ_scale_p
      REAL(KIND=dp), DIMENSION(nrot)                     :: grad_m_vec, grad_p_vec
      REAL(KIND=dp), DIMENSION(n, n)                     :: grad_im_2d, grad_re_2d
      REAL(KIND=dp), DIMENSION(n, nrot)                  :: diag_response, projector_response, &
                                                            response_fd, zero_response
      REAL(KIND=dp), DIMENSION(n, n, nkp)                :: x_im, x_m_im, x_m_re, x_p_im, x_p_re, x_re

      base(:, :, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      base(1, 1, 1) = CMPLX(-0.28_dp, 0.0_dp, KIND=dp)
      base(2, 2, 1) = CMPLX(0.03_dp, 0.0_dp, KIND=dp)
      base(3, 3, 1) = CMPLX(0.22_dp, 0.0_dp, KIND=dp)
      base(1, 2, 1) = CMPLX(0.060_dp, -0.020_dp, KIND=dp)
      base(1, 3, 1) = CMPLX(-0.025_dp, 0.035_dp, KIND=dp)
      base(2, 3, 1) = CMPLX(0.040_dp, 0.025_dp, KIND=dp)
      base(2, 1, 1) = CONJG(base(1, 2, 1))
      base(3, 1, 1) = CONJG(base(1, 3, 1))
      base(3, 2, 1) = CONJG(base(2, 3, 1))

      x_re(:, :, :) = 0.0_dp
      x_re(1, 2, 1) = 0.040_dp
      x_re(1, 3, 1) = -0.035_dp
      x_re(2, 3, 1) = 0.020_dp
      DO j = 2, n
         DO i = 1, j - 1
            x_re(j, i, 1) = -x_re(i, j, 1)
         END DO
      END DO

      x_im(:, :, :) = 0.0_dp
      x_im(1, 1, 1) = 0.018_dp
      x_im(2, 2, 1) = -0.022_dp
      x_im(3, 3, 1) = 0.014_dp
      x_im(1, 2, 1) = -0.018_dp
      x_im(1, 3, 1) = 0.016_dp
      x_im(2, 3, 1) = 0.026_dp
      DO j = 2, n
         DO i = 1, j - 1
            x_im(j, i, 1) = x_im(i, j, 1)
         END DO
      END DO

      ivar = 0
      DO j = 2, n
         DO i = 1, j - 1
            ivar = ivar + 1
            var_i(ivar) = i
            var_j(ivar) = j
            var_kind(ivar) = 1
            ivar = ivar + 1
            var_i(ivar) = i
            var_j(ivar) = j
            var_kind(ivar) = 2
         END DO
      END DO

      occ_scale(:) = [1.42_dp, 0.77_dp, 0.31_dp]

      CALL dense_blended_projector_response(x_re(:, :, 1), x_im(:, :, 1), base(:, :, 1), &
                                            0.0_dp, var_i, var_j, var_kind, zero_response)
      CALL dense_blended_projector_response(x_re(:, :, 1), x_im(:, :, 1), base(:, :, 1), &
                                            projector_weight, var_i, var_j, var_kind, projector_response)
      response_fd(:, :) = 0.0_dp
      err = 0.0_dp
      active_err = 0.0_dp
      DO s = 1, nrot
         x_p_re(:, :, :) = x_re(:, :, :)
         x_m_re(:, :, :) = x_re(:, :, :)
         x_p_im(:, :, :) = x_im(:, :, :)
         x_m_im(:, :, :) = x_im(:, :, :)
         CALL perturb_complex_rotation_var(x_p_re, x_p_im, var_i(s), var_j(s), var_kind(s), fd_step)
         CALL perturb_complex_rotation_var(x_m_re, x_m_im, var_i(s), var_j(s), var_kind(s), -fd_step)

         CALL dense_complex_kpoint_hdiag(x_p_re, x_p_im, base, h_p)
         CALL dense_complex_kpoint_hdiag(x_m_re, x_m_im, base, h_m)
         diag_response(:, s) = (h_p(:, 1) - h_m(:, 1))/(2.0_dp*fd_step)
         err = MAX(err, MAXVAL(ABS(zero_response(:, s) - diag_response(:, s))))
         active_err = MAX(active_err, MAXVAL(ABS(projector_response(:, s) - diag_response(:, s))))
      END DO

      DO i = 1, n
         occ_scale_p(:) = occ_scale(:)
         occ_scale_m(:) = occ_scale(:)
         occ_scale_p(i) = occ_scale_p(i) + fd_occ_step
         occ_scale_m(i) = occ_scale_m(i) - fd_occ_step
         CALL dense_blended_projector_rotation_gradient(x_re(:, :, 1), x_im(:, :, 1), base(:, :, 1), &
                                                        occ_scale_p, projector_weight, grad_re_2d, &
                                                        grad_im_2d, h_p(:, 1))
         CALL extract_rotation_gradient_vector_2d(grad_re_2d, grad_im_2d, var_i, var_j, &
                                                  var_kind, grad_p_vec)
         CALL dense_blended_projector_rotation_gradient(x_re(:, :, 1), x_im(:, :, 1), base(:, :, 1), &
                                                        occ_scale_m, projector_weight, grad_re_2d, &
                                                        grad_im_2d, h_m(:, 1))
         CALL extract_rotation_gradient_vector_2d(grad_re_2d, grad_im_2d, var_i, var_j, &
                                                  var_kind, grad_m_vec)
         response_fd(i, :) = 2.0_dp*(grad_p_vec(:) - grad_m_vec(:))/(2.0_dp*fd_occ_step)
      END DO

      err = MAX(err, MAXVAL(ABS(projector_response(:, :) - response_fd(:, :))))

      IF (err > tol .OR. active_err < active_tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,2(1X,ES12.5))') TRIM(label), "failed with error/active", err, active_err
      ELSE
         WRITE (*, '(A,1X,A,2(1X,ES12.5))') TRIM(label), "OK with error/active", err, active_err
      END IF

   END SUBROUTINE check_complex_projector_blended_response_channel

! **************************************************************************************************
!> \brief Check the finite-complex physical occupation response against eigenvalue differences.
!>        This is the gauge-stable companion to the band-projector basis response:
!>        the diagonal part remains the Rayleigh response, while the spectral part
!>        follows finite eigenvalue shifts rather than the instantaneous projector basis.
! **************************************************************************************************
   SUBROUTINE check_complex_projector_physical_response_channel(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 3, nkp = 1, nrot = n*(n - 1)
      REAL(KIND=dp), PARAMETER                           :: active_tol = 1.0E-5_dp, &
                                                            fd_step = 1.0E-6_dp, projector_weight = 0.70_dp

      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: chc_m, chc_p, occ_mat_m, occ_mat_p, u
      COMPLEX(KIND=dp), DIMENSION(n, n, nkp)             :: base
      INTEGER                                            :: i, ivar, j, s
      INTEGER, DIMENSION(nrot)                           :: var_i, var_j, var_kind
      REAL(KIND=dp)                                      :: active_err, err
      REAL(KIND=dp), DIMENSION(n)                        :: occ_curv, occ_scale
      REAL(KIND=dp), DIMENSION(n, nkp)                   :: diag_m, diag_p, eig_m, eig_p
      REAL(KIND=dp), DIMENSION(n, nrot)                  :: basis_response, physical_response, reference
      REAL(KIND=dp), DIMENSION(n, n, nkp)                :: x_im, x_m_im, x_m_re, x_p_im, x_p_re, x_re

      base(:, :, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      base(1, 1, 1) = CMPLX(-0.015_dp, 0.0_dp, KIND=dp)
      base(2, 2, 1) = CMPLX(0.011_dp, 0.0_dp, KIND=dp)
      base(3, 3, 1) = CMPLX(0.38_dp, 0.0_dp, KIND=dp)
      base(1, 2, 1) = CMPLX(0.020_dp, -0.035_dp, KIND=dp)
      base(1, 3, 1) = CMPLX(-0.045_dp, 0.015_dp, KIND=dp)
      base(2, 3, 1) = CMPLX(0.030_dp, 0.040_dp, KIND=dp)
      base(2, 1, 1) = CONJG(base(1, 2, 1))
      base(3, 1, 1) = CONJG(base(1, 3, 1))
      base(3, 2, 1) = CONJG(base(2, 3, 1))

      x_re(:, :, :) = 0.0_dp
      x_re(1, 2, 1) = 0.055_dp
      x_re(1, 3, 1) = -0.025_dp
      x_re(2, 3, 1) = 0.030_dp
      DO j = 2, n
         DO i = 1, j - 1
            x_re(j, i, 1) = -x_re(i, j, 1)
         END DO
      END DO

      x_im(:, :, :) = 0.0_dp
      x_im(1, 1, 1) = 0.012_dp
      x_im(2, 2, 1) = -0.018_dp
      x_im(3, 3, 1) = 0.021_dp
      x_im(1, 2, 1) = -0.028_dp
      x_im(1, 3, 1) = 0.019_dp
      x_im(2, 3, 1) = 0.024_dp
      DO j = 2, n
         DO i = 1, j - 1
            x_im(j, i, 1) = x_im(i, j, 1)
         END DO
      END DO

      ivar = 0
      DO j = 2, n
         DO i = 1, j - 1
            ivar = ivar + 1
            var_i(ivar) = i
            var_j(ivar) = j
            var_kind(ivar) = 1
            ivar = ivar + 1
            var_i(ivar) = i
            var_j(ivar) = j
            var_kind(ivar) = 2
         END DO
      END DO

      occ_curv(:) = [0.31_dp, 0.72_dp, 0.19_dp]
      occ_scale(:) = [1.42_dp, 0.77_dp, 0.31_dp]
      CALL dense_blended_projector_response(x_re(:, :, 1), x_im(:, :, 1), base(:, :, 1), &
                                            projector_weight, var_i, var_j, var_kind, basis_response)
      CALL dense_blended_projector_physical_response(x_re(:, :, 1), x_im(:, :, 1), base(:, :, 1), &
                                                     projector_weight, occ_scale, occ_curv, var_i, var_j, var_kind, &
                                                     physical_response)
      DO s = 1, nrot
         x_p_re(:, :, :) = x_re(:, :, :)
         x_m_re(:, :, :) = x_re(:, :, :)
         x_p_im(:, :, :) = x_im(:, :, :)
         x_m_im(:, :, :) = x_im(:, :, :)
         CALL perturb_complex_rotation_var(x_p_re, x_p_im, var_i(s), var_j(s), var_kind(s), fd_step)
         CALL perturb_complex_rotation_var(x_m_re, x_m_im, var_i(s), var_j(s), var_kind(s), -fd_step)
         CALL dense_complex_kpoint_hdiag(x_p_re, x_p_im, base, diag_p)
         CALL dense_complex_kpoint_hdiag(x_m_re, x_m_im, base, diag_m)
         CALL dense_complex_kpoint_eigvals(x_p_re, x_p_im, base, eig_p)
         CALL dense_complex_kpoint_eigvals(x_m_re, x_m_im, base, eig_m)
         CALL dense_rotation_u(x_p_re(:, :, 1), x_p_im(:, :, 1), u)
         chc_p(:, :) = MATMUL(CONJG(TRANSPOSE(u)), MATMUL(base(:, :, 1), u))
         CALL dense_rotation_u(x_m_re(:, :, 1), x_m_im(:, :, 1), u)
         chc_m(:, :) = MATMUL(CONJG(TRANSPOSE(u)), MATMUL(base(:, :, 1), u))
         CALL dense_spectral_occ_projector_from_chc(chc_p, occ_scale, occ_mat_p)
         CALL dense_spectral_occ_projector_from_chc(chc_m, occ_scale, occ_mat_m)
         reference(:, s) = (1.0_dp - projector_weight)*(diag_p(:, 1) - diag_m(:, 1))/(2.0_dp*fd_step) + &
                           projector_weight*(eig_p(:, 1) - eig_m(:, 1))/(2.0_dp*fd_step)
         DO i = 1, n
            reference(i, s) = reference(i, s) - projector_weight* &
                              0.5_dp*REAL(occ_mat_p(i, i) - occ_mat_m(i, i), KIND=dp)/ &
                              (2.0_dp*fd_step*occ_curv(i))
         END DO
      END DO

      err = MAXVAL(ABS(physical_response(:, :) - reference(:, :)))
      active_err = MAXVAL(ABS(physical_response(:, :) - basis_response(:, :)))

      IF (err > tol .OR. active_err < active_tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,2(1X,ES12.5))') TRIM(label), "failed with error/active", err, active_err
      ELSE
         WRITE (*, '(A,1X,A,2(1X,ES12.5))') TRIM(label), "OK with error/active", err, active_err
      END IF

   END SUBROUTINE check_complex_projector_physical_response_channel

! **************************************************************************************************
!> \brief Dense derivative of the blended projector rotation gradient with respect to occupations.
! **************************************************************************************************
   SUBROUTINE dense_blended_projector_response(x_re, x_im, base, projector_weight, &
                                               var_i, var_j, var_kind, response)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: x_re, x_im
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: base
      REAL(KIND=dp), INTENT(IN)                          :: projector_weight
      INTEGER, DIMENSION(:), INTENT(IN)                  :: var_i, var_j, var_kind
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: response

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: chc, evec, occ_mat, u
      INTEGER                                            :: band, i, j, n, r
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eps_diag, evals
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: grad_im, grad_re

      n = SIZE(x_re, 1)
      ALLOCATE (chc(n, n), evec(n, n), occ_mat(n, n), u(n, n), evals(n))
      ALLOCATE (grad_re(n, n), grad_im(n, n), eps_diag(n))

      CALL dense_rotation_u(x_re, x_im, u)
      chc(:, :) = MATMUL(CONJG(TRANSPOSE(u)), MATMUL(base, u))
      IF (projector_weight > 0.0_dp) CALL diag_complex(chc, evec, evals)

      DO band = 1, n
         occ_mat(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
         IF (projector_weight > 0.0_dp) THEN
            DO j = 1, n
               DO i = 1, n
                  occ_mat(i, j) = projector_weight*evec(i, band)*CONJG(evec(j, band))
               END DO
            END DO
         END IF
         occ_mat(band, band) = occ_mat(band, band) + &
                                CMPLX(1.0_dp - projector_weight, 0.0_dp, KIND=dp)
         CALL dense_projector_rotation_gradient(x_re, x_im, base, occ_mat, grad_re, grad_im, eps_diag)
         DO r = 1, SIZE(var_i)
            SELECT CASE (var_kind(r))
            CASE (1)
               response(band, r) = 2.0_dp*grad_re(var_i(r), var_j(r))
            CASE (2)
               response(band, r) = 2.0_dp*grad_im(var_i(r), var_j(r))
            END SELECT
         END DO
      END DO

      DEALLOCATE (chc, evec, occ_mat, u, evals, grad_re, grad_im, eps_diag)

   END SUBROUTINE dense_blended_projector_response

! **************************************************************************************************
!> \brief Dense finite-complex physical occupation response from Rayleigh/eigenvalue shifts.
! **************************************************************************************************
   SUBROUTINE dense_blended_projector_physical_response(x_re, x_im, base, projector_weight, occ_scale, &
                                                        occ_curv, var_i, var_j, var_kind, response)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: x_re, x_im
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: base
      REAL(KIND=dp), INTENT(IN)                          :: projector_weight
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: occ_curv, occ_scale
      INTEGER, DIMENSION(:), INTENT(IN)                  :: var_i, var_j, var_kind
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: response

      REAL(KIND=dp), PARAMETER                           :: fd_step = 1.0E-6_dp

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: chc_m, chc_p, evec, occ_mat_m, occ_mat_p, u
      INTEGER                                            :: band, n, s
      REAL(KIND=dp)                                      :: projector_response, raw_response
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: evals_m, evals_p
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: x_m_im, x_m_re, x_p_im, x_p_re

      n = SIZE(x_re, 1)
      ALLOCATE (chc_p(n, n), chc_m(n, n), evec(n, n), occ_mat_p(n, n), occ_mat_m(n, n), u(n, n), &
                evals_p(n), evals_m(n))
      ALLOCATE (x_p_re(n, n), x_p_im(n, n), x_m_re(n, n), x_m_im(n, n))

      DO s = 1, SIZE(var_i)
         x_p_re(:, :) = x_re(:, :)
         x_p_im(:, :) = x_im(:, :)
         x_m_re(:, :) = x_re(:, :)
         x_m_im(:, :) = x_im(:, :)
         CALL perturb_complex_rotation_var_2d(x_p_re, x_p_im, var_i(s), var_j(s), var_kind(s), fd_step)
         CALL perturb_complex_rotation_var_2d(x_m_re, x_m_im, var_i(s), var_j(s), var_kind(s), -fd_step)

         CALL dense_rotation_u(x_p_re, x_p_im, u)
         chc_p(:, :) = MATMUL(CONJG(TRANSPOSE(u)), MATMUL(base, u))
         CALL dense_rotation_u(x_m_re, x_m_im, u)
         chc_m(:, :) = MATMUL(CONJG(TRANSPOSE(u)), MATMUL(base, u))
         IF (projector_weight > 0.0_dp) THEN
            CALL diag_complex(chc_p, evec, evals_p)
            CALL diag_complex(chc_m, evec, evals_m)
            CALL dense_spectral_occ_projector_from_chc(chc_p, occ_scale, occ_mat_p)
            CALL dense_spectral_occ_projector_from_chc(chc_m, occ_scale, occ_mat_m)
         END IF

         DO band = 1, n
            raw_response = (REAL(chc_p(band, band), KIND=dp) - &
                            REAL(chc_m(band, band), KIND=dp))/(2.0_dp*fd_step)
            projector_response = raw_response
            IF (projector_weight > 0.0_dp) &
               projector_response = (evals_p(band) - evals_m(band))/(2.0_dp*fd_step)
            response(band, s) = (1.0_dp - projector_weight)*raw_response + &
                                projector_weight*projector_response
            IF (projector_weight > 0.0_dp .AND. ABS(occ_curv(band)) > EPSILON(occ_curv(band))) &
               response(band, s) = response(band, s) - projector_weight* &
                                   0.5_dp*REAL(occ_mat_p(band, band) - occ_mat_m(band, band), KIND=dp)/ &
                                   (2.0_dp*fd_step*occ_curv(band))
         END DO
      END DO

      DEALLOCATE (chc_p, chc_m, evec, occ_mat_p, occ_mat_m, u, evals_p, evals_m, &
                  x_p_re, x_p_im, x_m_re, x_m_im)

   END SUBROUTINE dense_blended_projector_physical_response

! **************************************************************************************************
!> \brief Dense gradient for the same blended projector model with arbitrary occupation weights.
! **************************************************************************************************
   SUBROUTINE dense_blended_projector_rotation_gradient(x_re, x_im, base, occ_scale, &
                                                        projector_weight, grad_re, grad_im, eps_diag)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: x_re, x_im
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: base
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: occ_scale
      REAL(KIND=dp), INTENT(IN)                          :: projector_weight
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: grad_re, grad_im
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: eps_diag

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: occ_mat
      INTEGER                                            :: i, n

      n = SIZE(x_re, 1)
      ALLOCATE (occ_mat(n, n))
      IF (projector_weight > 0.0_dp) THEN
         CALL dense_spectral_occ_projector(x_re, x_im, base, occ_scale, occ_mat)
         occ_mat(:, :) = projector_weight*occ_mat(:, :)
      ELSE
         occ_mat(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      END IF
      DO i = 1, n
         occ_mat(i, i) = occ_mat(i, i) + &
                         CMPLX((1.0_dp - projector_weight)*occ_scale(i), 0.0_dp, KIND=dp)
      END DO
      CALL dense_projector_rotation_gradient(x_re, x_im, base, occ_mat, grad_re, grad_im, eps_diag)
      DEALLOCATE (occ_mat)

   END SUBROUTINE dense_blended_projector_rotation_gradient

! **************************************************************************************************
!> \brief Check the split Schur projector for global fixed-N Mermin response plus channel-local tether.
! **************************************************************************************************
   SUBROUTINE check_mixed_fixed_n_tether_projector(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: ne = 5, nrot = 3, nchannels = 2

      INTEGER                                            :: first_energy, i, j, local_channel, nmo, r, s
      INTEGER, DIMENSION(nchannels)                      :: channel_nmo, e_offset
      REAL(KIND=dp)                                      :: err, occ_sum, projected_term, ref_term, &
                                                            shift_left, shift_right, weight_sum
      REAL(KIND=dp), DIMENSION(ne)                       :: occ_curv, tether_curv
      REAL(KIND=dp), DIMENSION(ne, ne)                   :: hproj
      REAL(KIND=dp), DIMENSION(ne, nrot)                 :: response
      REAL(KIND=dp), DIMENSION(nrot, nrot)               :: projected, reference

      channel_nmo = [3, 2]
      e_offset = [1, 4]
      occ_curv = [0.31_dp, 0.74_dp, 0.16_dp, 0.47_dp, 0.22_dp]
      tether_curv = [0.05_dp, 0.02_dp, 0.08_dp, 0.07_dp, 0.03_dp]
      response(:, 1) = [0.13_dp, -0.21_dp, 0.09_dp, 0.31_dp, -0.17_dp]
      response(:, 2) = [-0.05_dp, 0.18_dp, 0.24_dp, -0.11_dp, 0.07_dp]
      response(:, 3) = [0.22_dp, 0.04_dp, -0.19_dp, 0.15_dp, -0.08_dp]

      hproj(:, :) = 0.0_dp
      DO i = 1, ne
         hproj(i, i) = occ_curv(i) + tether_curv(i)
      END DO
      occ_sum = SUM(occ_curv)
      DO j = 1, ne
         DO i = 1, ne
            hproj(i, j) = hproj(i, j) - occ_curv(i)*occ_curv(j)/occ_sum
         END DO
      END DO
      DO local_channel = 1, nchannels
         nmo = channel_nmo(local_channel)
         first_energy = e_offset(local_channel)
         weight_sum = SUM(tether_curv(first_energy:first_energy + nmo - 1))
         DO j = first_energy, first_energy + nmo - 1
            DO i = first_energy, first_energy + nmo - 1
               hproj(i, j) = hproj(i, j) - tether_curv(i)*tether_curv(j)/weight_sum
            END DO
         END DO
      END DO

      DO s = 1, nrot
         DO r = 1, nrot
            ref_term = DOT_PRODUCT(response(:, r), MATMUL(hproj, response(:, s)))
            projected_term = DOT_PRODUCT((occ_curv(:) + tether_curv(:))*response(:, r), response(:, s))

            shift_left = DOT_PRODUCT(occ_curv, response(:, r))
            shift_right = DOT_PRODUCT(occ_curv, response(:, s))
            projected_term = projected_term - shift_left*shift_right/occ_sum

            DO local_channel = 1, nchannels
               nmo = channel_nmo(local_channel)
               first_energy = e_offset(local_channel)
               weight_sum = SUM(tether_curv(first_energy:first_energy + nmo - 1))
               shift_left = DOT_PRODUCT(tether_curv(first_energy:first_energy + nmo - 1), &
                                        response(first_energy:first_energy + nmo - 1, r))
               shift_right = DOT_PRODUCT(tether_curv(first_energy:first_energy + nmo - 1), &
                                         response(first_energy:first_energy + nmo - 1, s))
               projected_term = projected_term - shift_left*shift_right/weight_sum
            END DO

            reference(r, s) = ref_term
            projected(r, s) = projected_term
         END DO
      END DO

      err = MAXVAL(ABS(projected(:, :) - reference(:, :)))
      IF (err > tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,ES12.5)') TRIM(label), "failed with error", err
      ELSE
         WRITE (*, '(A,1X,A,1X,ES12.5)') TRIM(label), "OK with error", err
      END IF

   END SUBROUTINE check_mixed_fixed_n_tether_projector

! **************************************************************************************************
!> \brief Check split Schur projector with distinct occupation and tether response channels.
! **************************************************************************************************
   SUBROUTINE check_mixed_fixed_n_tether_cross_projector(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: ne = 5, nrot = 3, nchannels = 2

      INTEGER                                            :: first_energy, i, j, local_channel, nmo, r, s
      INTEGER, DIMENSION(nchannels)                      :: channel_nmo, e_offset
      REAL(KIND=dp)                                      :: err, occ_sum, projected_term, ref_term, &
                                                            shift_left, shift_right, weight_sum
      REAL(KIND=dp), DIMENSION(ne)                       :: occ_curv, tether_curv
      REAL(KIND=dp), DIMENSION(ne, ne)                   :: hocc, htether
      REAL(KIND=dp), DIMENSION(ne, nrot)                 :: response_occ, response_tether
      REAL(KIND=dp), DIMENSION(nrot, nrot)               :: projected, reference

      channel_nmo = [3, 2]
      e_offset = [1, 4]
      occ_curv = [0.31_dp, 0.74_dp, 0.16_dp, 0.47_dp, 0.22_dp]
      tether_curv = [0.05_dp, 0.02_dp, 0.08_dp, 0.07_dp, 0.03_dp]
      response_occ(:, 1) = [0.13_dp, -0.21_dp, 0.09_dp, 0.31_dp, -0.17_dp]
      response_occ(:, 2) = [-0.05_dp, 0.18_dp, 0.24_dp, -0.11_dp, 0.07_dp]
      response_occ(:, 3) = [0.22_dp, 0.04_dp, -0.19_dp, 0.15_dp, -0.08_dp]
      response_tether(:, 1) = [0.06_dp, -0.09_dp, 0.18_dp, 0.27_dp, -0.14_dp]
      response_tether(:, 2) = [-0.12_dp, 0.21_dp, 0.05_dp, -0.08_dp, 0.16_dp]
      response_tether(:, 3) = [0.19_dp, -0.02_dp, -0.11_dp, 0.23_dp, -0.04_dp]

      hocc(:, :) = 0.0_dp
      htether(:, :) = 0.0_dp
      DO i = 1, ne
         hocc(i, i) = occ_curv(i)
         htether(i, i) = tether_curv(i)
      END DO
      occ_sum = SUM(occ_curv)
      DO j = 1, ne
         DO i = 1, ne
            hocc(i, j) = hocc(i, j) - occ_curv(i)*occ_curv(j)/occ_sum
         END DO
      END DO
      DO local_channel = 1, nchannels
         nmo = channel_nmo(local_channel)
         first_energy = e_offset(local_channel)
         weight_sum = SUM(tether_curv(first_energy:first_energy + nmo - 1))
         DO j = first_energy, first_energy + nmo - 1
            DO i = first_energy, first_energy + nmo - 1
               htether(i, j) = htether(i, j) - tether_curv(i)*tether_curv(j)/weight_sum
            END DO
         END DO
      END DO

      DO s = 1, nrot
         DO r = 1, nrot
            ref_term = DOT_PRODUCT(response_occ(:, r), MATMUL(hocc, response_occ(:, s))) + &
                       DOT_PRODUCT(response_tether(:, r), MATMUL(htether, response_tether(:, s)))

            projected_term = DOT_PRODUCT(occ_curv(:)*response_occ(:, r), response_occ(:, s))
            shift_left = DOT_PRODUCT(occ_curv, response_occ(:, r))
            shift_right = DOT_PRODUCT(occ_curv, response_occ(:, s))
            projected_term = projected_term - shift_left*shift_right/occ_sum

            projected_term = projected_term + DOT_PRODUCT(tether_curv(:)*response_tether(:, r), &
                                                          response_tether(:, s))
            DO local_channel = 1, nchannels
               nmo = channel_nmo(local_channel)
               first_energy = e_offset(local_channel)
               weight_sum = SUM(tether_curv(first_energy:first_energy + nmo - 1))
               shift_left = DOT_PRODUCT(tether_curv(first_energy:first_energy + nmo - 1), &
                                        response_tether(first_energy:first_energy + nmo - 1, r))
               shift_right = DOT_PRODUCT(tether_curv(first_energy:first_energy + nmo - 1), &
                                         response_tether(first_energy:first_energy + nmo - 1, s))
               projected_term = projected_term - shift_left*shift_right/weight_sum
            END DO

            reference(r, s) = ref_term
            projected(r, s) = projected_term
         END DO
      END DO

      err = MAXVAL(ABS(projected(:, :) - reference(:, :)))
      IF (err > tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,ES12.5)') TRIM(label), "failed with error", err
      ELSE
         WRITE (*, '(A,1X,A,1X,ES12.5)') TRIM(label), "OK with error", err
      END IF

   END SUBROUTINE check_mixed_fixed_n_tether_cross_projector

! **************************************************************************************************
!> \brief Check the symmetrized Schur term for distinct rotation/energy responses.
! **************************************************************************************************
   SUBROUTINE check_mixed_fixed_n_asymmetric_response_schur(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: ne = 5, nrot = 3

      INTEGER                                            :: i, j, r, s
      REAL(KIND=dp)                                      :: active_err, err, occ_sum, shift, symmetry_err
      REAL(KIND=dp), DIMENSION(ne)                       :: occ_curv
      REAL(KIND=dp), DIMENSION(ne, ne)                   :: hocc
      REAL(KIND=dp), DIMENSION(ne, nrot)                 :: b_occ, b_raw, response_occ, response_occ_projected, &
                                                            response_raw, response_raw_projected
      REAL(KIND=dp), DIMENSION(nrot, nrot)               :: mixed, occ_only, reference

      occ_curv = [0.31_dp, 0.74_dp, 0.16_dp, 0.47_dp, 0.22_dp]
      response_occ(:, 1) = [0.13_dp, -0.21_dp, 0.09_dp, 0.31_dp, -0.17_dp]
      response_occ(:, 2) = [-0.05_dp, 0.18_dp, 0.24_dp, -0.11_dp, 0.07_dp]
      response_occ(:, 3) = [0.22_dp, 0.04_dp, -0.19_dp, 0.15_dp, -0.08_dp]
      response_raw(:, 1) = [0.10_dp, -0.16_dp, 0.14_dp, 0.22_dp, -0.11_dp]
      response_raw(:, 2) = [-0.01_dp, 0.15_dp, 0.20_dp, -0.03_dp, 0.02_dp]
      response_raw(:, 3) = [0.17_dp, 0.08_dp, -0.13_dp, 0.10_dp, -0.03_dp]

      hocc(:, :) = 0.0_dp
      DO i = 1, ne
         hocc(i, i) = occ_curv(i)
      END DO
      occ_sum = SUM(occ_curv)
      DO j = 1, ne
         DO i = 1, ne
            hocc(i, j) = hocc(i, j) - occ_curv(i)*occ_curv(j)/occ_sum
         END DO
      END DO

      response_occ_projected(:, :) = response_occ(:, :)
      response_raw_projected(:, :) = response_raw(:, :)
      DO r = 1, nrot
         shift = DOT_PRODUCT(occ_curv, response_occ_projected(:, r))/occ_sum
         response_occ_projected(:, r) = response_occ_projected(:, r) - shift
         shift = DOT_PRODUCT(occ_curv, response_raw_projected(:, r))/occ_sum
         response_raw_projected(:, r) = response_raw_projected(:, r) - shift
      END DO

      b_occ(:, :) = MATMUL(hocc, response_occ)
      b_raw(:, :) = MATMUL(hocc, response_raw)
      DO s = 1, nrot
         DO r = 1, nrot
            mixed(r, s) = 0.5_dp*(DOT_PRODUCT(b_occ(:, r), response_raw_projected(:, s)) + &
                                  DOT_PRODUCT(b_raw(:, r), response_occ_projected(:, s)))
            reference(r, s) = 0.5_dp*(DOT_PRODUCT(response_occ(:, r), MATMUL(hocc, response_raw(:, s))) + &
                                      DOT_PRODUCT(response_raw(:, r), MATMUL(hocc, response_occ(:, s))))
            occ_only(r, s) = DOT_PRODUCT(response_occ(:, r), MATMUL(hocc, response_occ(:, s)))
         END DO
      END DO

      symmetry_err = MAXVAL(ABS(mixed(:, :) - TRANSPOSE(mixed(:, :))))
      active_err = MAXVAL(ABS(mixed(:, :) - occ_only(:, :)))
      err = MAX(MAXVAL(ABS(mixed(:, :) - reference(:, :))), symmetry_err)

      IF (err > tol .OR. active_err < 1.0E-6_dp) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,2(1X,ES12.5))') TRIM(label), "failed with error/active", err, active_err
      ELSE
         WRITE (*, '(A,1X,A,2(1X,ES12.5))') TRIM(label), "OK with error/active", err, active_err
      END IF

   END SUBROUTINE check_mixed_fixed_n_asymmetric_response_schur

! **************************************************************************************************
!> \brief Check exact dense elimination for the mixed fixed-N/tether energy block.
! **************************************************************************************************
   SUBROUTINE check_mixed_fixed_n_tether_exact_response(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: ne = 6, nrot = 3, nchannels = 2

      INTEGER                                            :: first_energy, i, info, j, local_channel, nmo, r
      INTEGER, DIMENSION(ne)                             :: ipiv
      INTEGER, DIMENSION(nchannels)                      :: channel_nmo, e_offset
      REAL(KIND=dp)                                      :: active_err, channel_occ_sum, curv_max, err, &
                                                            gauge_weight, ref_value, shift, test_value, &
                                                            weight_sum
      REAL(KIND=dp), DIMENSION(ne)                       :: bq, occ_curv, tether_curv, x
      REAL(KIND=dp), DIMENSION(nrot)                     :: q
      REAL(KIND=dp), DIMENSION(ne, ne)                   :: h, h_ref
      REAL(KIND=dp), DIMENSION(ne, nrot)                 :: b, response_occ, response_tether, rhs
      REAL(KIND=dp), DIMENSION(nrot, nrot)               :: exact

      channel_nmo = [3, 3]
      e_offset = [1, 4]
      occ_curv = [0.31_dp, 0.74_dp, 0.16_dp, 0.0_dp, 0.0_dp, 0.0_dp]
      tether_curv = [0.05_dp, 0.02_dp, 0.08_dp, 0.07_dp, 0.03_dp, 0.0_dp]
      response_occ(:, 1) = [0.13_dp, -0.21_dp, 0.09_dp, 0.31_dp, -0.17_dp, 0.05_dp]
      response_occ(:, 2) = [-0.05_dp, 0.18_dp, 0.24_dp, -0.11_dp, 0.07_dp, -0.09_dp]
      response_occ(:, 3) = [0.22_dp, 0.04_dp, -0.19_dp, 0.15_dp, -0.08_dp, 0.12_dp]
      response_tether(:, 1) = [0.06_dp, -0.09_dp, 0.18_dp, 0.27_dp, -0.14_dp, 0.11_dp]
      response_tether(:, 2) = [-0.12_dp, 0.21_dp, 0.05_dp, -0.08_dp, 0.16_dp, -0.03_dp]
      response_tether(:, 3) = [0.19_dp, -0.02_dp, -0.11_dp, 0.23_dp, -0.04_dp, 0.08_dp]
      q = [0.21_dp, -0.13_dp, 0.17_dp]

      h(:, :) = 0.0_dp
      b(:, :) = 0.0_dp
      DO i = 1, ne
         h(i, i) = occ_curv(i) + tether_curv(i)
      END DO
      weight_sum = SUM(occ_curv)
      DO j = 1, ne
         DO i = 1, ne
            h(i, j) = h(i, j) - occ_curv(i)*occ_curv(j)/weight_sum
         END DO
      END DO
      DO r = 1, nrot
         shift = DOT_PRODUCT(occ_curv, response_occ(:, r))/weight_sum
         b(:, r) = b(:, r) + occ_curv(:)*(response_occ(:, r) - shift)
      END DO
      DO local_channel = 1, nchannels
         nmo = channel_nmo(local_channel)
         first_energy = e_offset(local_channel)
         weight_sum = SUM(tether_curv(first_energy:first_energy + nmo - 1))
         DO j = first_energy, first_energy + nmo - 1
            DO i = first_energy, first_energy + nmo - 1
               h(i, j) = h(i, j) - tether_curv(i)*tether_curv(j)/weight_sum
            END DO
         END DO
         DO r = 1, nrot
            shift = DOT_PRODUCT(tether_curv(first_energy:first_energy + nmo - 1), &
                                response_tether(first_energy:first_energy + nmo - 1, r))/weight_sum
            b(first_energy:first_energy + nmo - 1, r) = b(first_energy:first_energy + nmo - 1, r) + &
               tether_curv(first_energy:first_energy + nmo - 1)* &
               (response_tether(first_energy:first_energy + nmo - 1, r) - shift)
         END DO
      END DO

      curv_max = MAXVAL(occ_curv + tether_curv)
      gauge_weight = curv_max/REAL(ne, KIND=dp)
      h(:, :) = h(:, :) + gauge_weight
      DO local_channel = 1, nchannels
         nmo = channel_nmo(local_channel)
         first_energy = e_offset(local_channel)
         weight_sum = SUM(tether_curv(first_energy:first_energy + nmo - 1))
         channel_occ_sum = SUM(occ_curv(first_energy:first_energy + nmo - 1))
         IF (weight_sum <= EPSILON(weight_sum) .OR. channel_occ_sum > EPSILON(channel_occ_sum)) CYCLE
         DO j = first_energy, first_energy + nmo - 1
            DO i = first_energy, first_energy + nmo - 1
               h(i, j) = h(i, j) + gauge_weight
            END DO
         END DO
      END DO
      DO i = 1, ne
         IF (occ_curv(i) + tether_curv(i) > EPSILON(curv_max)) CYCLE
         h(i, :) = 0.0_dp
         h(:, i) = 0.0_dp
         h(i, i) = 1.0_dp
         b(i, :) = 0.0_dp
      END DO
      rhs(:, :) = b(:, :)
      h_ref(:, :) = h(:, :)
      CALL dgesv(ne, nrot, h, ne, ipiv, rhs, ne, info)
      IF (info /= 0) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,I6)') TRIM(label), "failed in dgesv", info
         RETURN
      END IF

      exact(:, :) = MATMUL(TRANSPOSE(b), rhs)
      active_err = ABS(exact(1, 1) - DOT_PRODUCT(response_occ(:, 1), occ_curv(:)*response_occ(:, 1)))
      bq(:) = MATMUL(b, q)
      x(:) = MATMUL(rhs, q)
      ref_value = -0.5_dp*DOT_PRODUCT(q, MATMUL(exact, q))
      test_value = 0.5_dp*DOT_PRODUCT(x, MATMUL(h_ref, x)) - DOT_PRODUCT(x, bq)
      err = MAX(MAXVAL(ABS(exact(:, :) - TRANSPOSE(exact(:, :)))), ABS(ref_value - test_value))

      IF (err > tol .OR. active_err < 1.0E-4_dp) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,2(1X,ES12.5))') TRIM(label), "failed with error/active", err, active_err
      ELSE
         WRITE (*, '(A,1X,A,2(1X,ES12.5))') TRIM(label), "OK with error/active", err, active_err
      END IF

   END SUBROUTINE check_mixed_fixed_n_tether_exact_response

! **************************************************************************************************
!> \brief Check exact energy-block elimination against the full augmented Schur system.
! **************************************************************************************************
   SUBROUTINE check_mixed_fixed_n_tether_exact_schur(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: ne = 6, nrot = 3, nchannels = 2, ntot = ne + nrot

      INTEGER                                            :: info
      INTEGER, DIMENSION(MAX(ne, ntot))                  :: ipiv
      INTEGER, DIMENSION(nchannels)                      :: channel_nmo, e_offset
      REAL(KIND=dp)                                      :: err
      REAL(KIND=dp), DIMENSION(ne)                       :: grad_e, z_aug, z_red
      REAL(KIND=dp), DIMENSION(nrot)                     :: grad_rot, q_aug, q_red
      REAL(KIND=dp), DIMENSION(ne, ne)                   :: h
      REAL(KIND=dp), DIMENSION(ne, nrot)                 :: b, response_occ, response_tether, rhs_energy
      REAL(KIND=dp), DIMENSION(nrot, nrot)               :: a, schur
      REAL(KIND=dp), DIMENSION(ntot, ntot)               :: aug
      REAL(KIND=dp), DIMENSION(ntot)                     :: rhs_aug

      channel_nmo = [3, 3]
      e_offset = [1, 4]
      grad_e = [0.11_dp, -0.07_dp, 0.04_dp, -0.09_dp, 0.13_dp, 0.03_dp]
      grad_rot = [0.17_dp, -0.19_dp, 0.08_dp]
      response_occ(:, 1) = [0.13_dp, -0.21_dp, 0.09_dp, 0.31_dp, -0.17_dp, 0.05_dp]
      response_occ(:, 2) = [-0.05_dp, 0.18_dp, 0.24_dp, -0.11_dp, 0.07_dp, -0.09_dp]
      response_occ(:, 3) = [0.22_dp, 0.04_dp, -0.19_dp, 0.15_dp, -0.08_dp, 0.12_dp]
      response_tether(:, 1) = [0.06_dp, -0.09_dp, 0.18_dp, 0.27_dp, -0.14_dp, 0.11_dp]
      response_tether(:, 2) = [-0.12_dp, 0.21_dp, 0.05_dp, -0.08_dp, 0.16_dp, -0.03_dp]
      response_tether(:, 3) = [0.19_dp, -0.02_dp, -0.11_dp, 0.23_dp, -0.04_dp, 0.08_dp]

      CALL build_test_mixed_energy_block(channel_nmo, e_offset, response_occ, response_tether, h, b)

      a(:, :) = RESHAPE([1.37_dp, 0.04_dp, -0.03_dp, &
                         0.04_dp, 1.11_dp, 0.05_dp, &
                         -0.03_dp, 0.05_dp, 0.93_dp], SHAPE(a))
      rhs_energy(:, 1:nrot) = b(:, 1:nrot)
      CALL dgesv(ne, nrot, h, ne, ipiv(1:ne), rhs_energy, ne, info)
      IF (info /= 0) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,I6)') TRIM(label), "failed energy dgesv", info
         RETURN
      END IF

      schur(:, :) = a(:, :) - MATMUL(TRANSPOSE(b), rhs_energy)
      z_red(:) = grad_e(:)
      CALL build_test_mixed_energy_block(channel_nmo, e_offset, response_occ, response_tether, h, b)
      CALL dgesv(ne, 1, h, ne, ipiv(1:ne), z_red, ne, info)
      IF (info /= 0) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,I6)') TRIM(label), "failed gradient dgesv", info
         RETURN
      END IF
      q_red(:) = grad_rot(:) + MATMUL(TRANSPOSE(rhs_energy), grad_e)
      CALL dgesv(nrot, 1, schur, nrot, ipiv(1:nrot), q_red, nrot, info)
      IF (info /= 0) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,I6)') TRIM(label), "failed Schur dgesv", info
         RETURN
      END IF
      z_red(:) = z_red(:) + MATMUL(rhs_energy, q_red)

      CALL build_test_mixed_energy_block(channel_nmo, e_offset, response_occ, response_tether, h, b)
      aug(:, :) = 0.0_dp
      aug(1:ne, 1:ne) = h(:, :)
      aug(1:ne, ne + 1:ntot) = -b(:, :)
      aug(ne + 1:ntot, 1:ne) = -TRANSPOSE(b)
      aug(ne + 1:ntot, ne + 1:ntot) = a(:, :)
      rhs_aug(1:ne) = grad_e(:)
      rhs_aug(ne + 1:ntot) = grad_rot(:)
      CALL dgesv(ntot, 1, aug, ntot, ipiv(1:ntot), rhs_aug, ntot, info)
      IF (info /= 0) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,I6)') TRIM(label), "failed augmented dgesv", info
         RETURN
      END IF
      z_aug(:) = rhs_aug(1:ne)
      q_aug(:) = rhs_aug(ne + 1:ntot)

      err = MAX(MAXVAL(ABS(z_red - z_aug)), MAXVAL(ABS(q_red - q_aug)))

      IF (err > tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,ES12.5)') TRIM(label), "failed with error", err
      ELSE
         WRITE (*, '(A,1X,A,1X,ES12.5)') TRIM(label), "OK with error", err
      END IF

   END SUBROUTINE check_mixed_fixed_n_tether_exact_schur

! **************************************************************************************************
!> \brief Build the dense mixed fixed-N/tether energy block used by exact-response tests.
! **************************************************************************************************
   SUBROUTINE build_test_mixed_energy_block(channel_nmo, e_offset, response_occ, response_tether, h, b)
      INTEGER, DIMENSION(:), INTENT(IN)                  :: channel_nmo, e_offset
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: response_occ, response_tether
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: h, b

      INTEGER                                            :: first_energy, i, j, local_channel, ne, nchannels, nmo, nrot, r
      REAL(KIND=dp)                                      :: channel_occ_sum, curv_max, gauge_weight, shift, weight_sum
      REAL(KIND=dp), DIMENSION(SIZE(h, 1))               :: occ_curv, tether_curv

      ne = SIZE(h, 1)
      nrot = SIZE(b, 2)
      nchannels = SIZE(channel_nmo)
      occ_curv = [0.31_dp, 0.74_dp, 0.16_dp, 0.0_dp, 0.0_dp, 0.0_dp]
      tether_curv = [0.05_dp, 0.02_dp, 0.08_dp, 0.07_dp, 0.03_dp, 0.0_dp]
      h(:, :) = 0.0_dp
      b(:, :) = 0.0_dp
      DO i = 1, ne
         h(i, i) = occ_curv(i) + tether_curv(i)
      END DO
      weight_sum = SUM(occ_curv)
      DO j = 1, ne
         DO i = 1, ne
            h(i, j) = h(i, j) - occ_curv(i)*occ_curv(j)/weight_sum
         END DO
      END DO
      DO r = 1, nrot
         shift = DOT_PRODUCT(occ_curv, response_occ(:, r))/weight_sum
         b(:, r) = b(:, r) + occ_curv(:)*(response_occ(:, r) - shift)
      END DO
      DO local_channel = 1, nchannels
         nmo = channel_nmo(local_channel)
         first_energy = e_offset(local_channel)
         weight_sum = SUM(tether_curv(first_energy:first_energy + nmo - 1))
         DO j = first_energy, first_energy + nmo - 1
            DO i = first_energy, first_energy + nmo - 1
               h(i, j) = h(i, j) - tether_curv(i)*tether_curv(j)/weight_sum
            END DO
         END DO
         DO r = 1, nrot
            shift = DOT_PRODUCT(tether_curv(first_energy:first_energy + nmo - 1), &
                                response_tether(first_energy:first_energy + nmo - 1, r))/weight_sum
            b(first_energy:first_energy + nmo - 1, r) = b(first_energy:first_energy + nmo - 1, r) + &
               tether_curv(first_energy:first_energy + nmo - 1)* &
               (response_tether(first_energy:first_energy + nmo - 1, r) - shift)
         END DO
      END DO

      curv_max = MAXVAL(occ_curv + tether_curv)
      gauge_weight = curv_max/REAL(ne, KIND=dp)
      h(:, :) = h(:, :) + gauge_weight
      DO local_channel = 1, nchannels
         nmo = channel_nmo(local_channel)
         first_energy = e_offset(local_channel)
         weight_sum = SUM(tether_curv(first_energy:first_energy + nmo - 1))
         channel_occ_sum = SUM(occ_curv(first_energy:first_energy + nmo - 1))
         IF (weight_sum <= EPSILON(weight_sum) .OR. channel_occ_sum > EPSILON(channel_occ_sum)) CYCLE
         DO j = first_energy, first_energy + nmo - 1
            DO i = first_energy, first_energy + nmo - 1
               h(i, j) = h(i, j) + gauge_weight
            END DO
         END DO
      END DO
      DO i = 1, ne
         IF (occ_curv(i) + tether_curv(i) > EPSILON(curv_max)) CYCLE
         h(i, :) = 0.0_dp
         h(:, i) = 0.0_dp
         h(i, i) = 1.0_dp
         b(i, :) = 0.0_dp
      END DO

   END SUBROUTINE build_test_mixed_energy_block

! **************************************************************************************************
!> \brief Diagonal of U^H H U for each dense complex k-point channel.
! **************************************************************************************************
   SUBROUTINE dense_complex_kpoint_hdiag(x_re, x_im, base, hdiag)
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: x_re, x_im
      COMPLEX(KIND=dp), DIMENSION(:, :, :), INTENT(IN)   :: base
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: hdiag

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: chc, u
      INTEGER                                            :: i, ik, n, nkp

      n = SIZE(x_re, 1)
      nkp = SIZE(x_re, 3)
      ALLOCATE (chc(n, n), u(n, n))
      DO ik = 1, nkp
         CALL dense_rotation_u(x_re(:, :, ik), x_im(:, :, ik), u)
         chc(:, :) = MATMUL(CONJG(TRANSPOSE(u)), MATMUL(base(:, :, ik), u))
         DO i = 1, n
            hdiag(i, ik) = REAL(chc(i, i), KIND=dp)
         END DO
      END DO
      DEALLOCATE (chc, u)

   END SUBROUTINE dense_complex_kpoint_hdiag

! **************************************************************************************************
!> \brief Eigenvalues of U^H H U for each dense complex k-point channel.
! **************************************************************************************************
   SUBROUTINE dense_complex_kpoint_eigvals(x_re, x_im, base, eigvals)
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: x_re, x_im
      COMPLEX(KIND=dp), DIMENSION(:, :, :), INTENT(IN)   :: base
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: eigvals

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: chc, evec, u
      INTEGER                                            :: ik, n, nkp

      n = SIZE(x_re, 1)
      nkp = SIZE(x_re, 3)
      ALLOCATE (chc(n, n), evec(n, n), u(n, n))
      DO ik = 1, nkp
         CALL dense_rotation_u(x_re(:, :, ik), x_im(:, :, ik), u)
         chc(:, :) = MATMUL(CONJG(TRANSPOSE(u)), MATMUL(base(:, :, ik), u))
         CALL diag_complex(chc, evec, eigvals(:, ik))
      END DO
      DEALLOCATE (chc, evec, u)

   END SUBROUTINE dense_complex_kpoint_eigvals

! **************************************************************************************************
!> \brief Perturb one off-diagonal real or imaginary anti-Hermitian rotation variable.
! **************************************************************************************************
   SUBROUTINE perturb_complex_rotation_var(x_re, x_im, i, j, kind, step)
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: x_re, x_im
      INTEGER, INTENT(IN)                                :: i, j, kind
      REAL(KIND=dp), INTENT(IN)                          :: step

      SELECT CASE (kind)
      CASE (1)
         x_re(i, j, 1) = x_re(i, j, 1) + step
         x_re(j, i, 1) = x_re(j, i, 1) - step
      CASE (2)
         x_im(i, j, 1) = x_im(i, j, 1) + step
         x_im(j, i, 1) = x_im(j, i, 1) + step
      END SELECT

   END SUBROUTINE perturb_complex_rotation_var

! **************************************************************************************************
!> \brief Perturb one packed complex rotation variable in a dense 2-D fixture.
! **************************************************************************************************
   SUBROUTINE perturb_complex_rotation_var_2d(x_re, x_im, i, j, kind, step)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: x_re, x_im
      INTEGER, INTENT(IN)                                :: i, j, kind
      REAL(KIND=dp), INTENT(IN)                          :: step

      SELECT CASE (kind)
      CASE (1)
         x_re(i, j) = x_re(i, j) + step
         x_re(j, i) = x_re(j, i) - step
      CASE (2)
         x_im(i, j) = x_im(i, j) + step
         x_im(j, i) = x_im(j, i) + step
      END SELECT

   END SUBROUTINE perturb_complex_rotation_var_2d

! **************************************************************************************************
!> \brief Pack off-diagonal real/imaginary rotation gradients into a vector.
! **************************************************************************************************
   SUBROUTINE extract_rotation_gradient_vector(grad_re, grad_im, var_i, var_j, var_kind, grad_vec)
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: grad_re, grad_im
      INTEGER, DIMENSION(:), INTENT(IN)                  :: var_i, var_j, var_kind
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: grad_vec

      INTEGER                                            :: r

      DO r = 1, SIZE(var_i)
         SELECT CASE (var_kind(r))
         CASE (1)
            grad_vec(r) = grad_re(var_i(r), var_j(r), 1)
         CASE (2)
            grad_vec(r) = grad_im(var_i(r), var_j(r), 1)
         END SELECT
      END DO

   END SUBROUTINE extract_rotation_gradient_vector

! **************************************************************************************************
!> \brief Pack off-diagonal real/imaginary rotation gradients from a dense 2-D fixture.
! **************************************************************************************************
   SUBROUTINE extract_rotation_gradient_vector_2d(grad_re, grad_im, var_i, var_j, var_kind, grad_vec)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: grad_re, grad_im
      INTEGER, DIMENSION(:), INTENT(IN)                  :: var_i, var_j, var_kind
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: grad_vec

      INTEGER                                            :: r

      DO r = 1, SIZE(var_i)
         SELECT CASE (var_kind(r))
         CASE (1)
            grad_vec(r) = grad_re(var_i(r), var_j(r))
         CASE (2)
            grad_vec(r) = grad_im(var_i(r), var_j(r))
         END SELECT
      END DO

   END SUBROUTINE extract_rotation_gradient_vector_2d

! **************************************************************************************************
!> \brief Dense complex k-point Mermin value with occupations from auxiliary energies.
! **************************************************************************************************
   FUNCTION dense_complex_kpoint_mermin_value(x_re, x_im, base, eps, wk, nel, sigma, maxocc) &
      RESULT(value)
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: x_re, x_im
      COMPLEX(KIND=dp), DIMENSION(:, :, :), INTENT(IN)   :: base
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: eps
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: wk
      REAL(KIND=dp), INTENT(IN)                          :: nel, sigma, maxocc
      REAL(KIND=dp)                                      :: value

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: chc, u
      INTEGER                                            :: i, ik, n, nkp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: h

      n = SIZE(eps, 1)
      nkp = SIZE(eps, 2)
      ALLOCATE (chc(n, n), u(n, n), h(n, nkp))
      DO ik = 1, nkp
         CALL dense_rotation_u(x_re(:, :, ik), x_im(:, :, ik), u)
         chc(:, :) = MATMUL(CONJG(TRANSPOSE(u)), MATMUL(base(:, :, ik), u))
         DO i = 1, n
            h(i, ik) = REAL(chc(i, i), KIND=dp)
         END DO
      END DO
      value = dense_kpoint_mermin_value(eps, h, wk, nel, sigma, maxocc)
      DEALLOCATE (chc, u, h)

   END FUNCTION dense_complex_kpoint_mermin_value

! **************************************************************************************************
!> \brief Dense coupled Mermin gradients for complex k-point rotations and auxiliary energies.
! **************************************************************************************************
   SUBROUTINE dense_complex_kpoint_mermin_gradient(x_re, x_im, base, eps, wk, nel, sigma, maxocc, &
                                                   grad_re, grad_im, grad_eps, value)
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: x_re, x_im
      COMPLEX(KIND=dp), DIMENSION(:, :, :), INTENT(IN)   :: base
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: eps
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: wk
      REAL(KIND=dp), INTENT(IN)                          :: nel, sigma, maxocc
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(OUT)     :: grad_re, grad_im
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: grad_eps
      REAL(KIND=dp), INTENT(OUT)                         :: value

      INTEGER                                            :: ik, n, nkp
      REAL(KIND=dp)                                      :: kTS, mu
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: occ_scale
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: h, occ

      n = SIZE(eps, 1)
      nkp = SIZE(eps, 2)
      ALLOCATE (h(n, nkp), occ(n, nkp), occ_scale(n))
      CALL dense_fermi_kpoint(occ, mu, kTS, eps, nel, wk, sigma, maxocc)
      DO ik = 1, nkp
         occ_scale(:) = 2.0_dp*wk(ik)*occ(:, ik)
         CALL dense_ref_rotation_gradient(x_re(:, :, ik), x_im(:, :, ik), base(:, :, ik), occ_scale, &
                                          grad_re(:, :, ik), grad_im(:, :, ik), h(:, ik))
      END DO
      CALL dense_kpoint_mermin_gradient(eps, h, wk, nel, sigma, maxocc, grad_eps)
      value = dense_kpoint_mermin_value(eps, h, wk, nel, sigma, maxocc)
      DEALLOCATE (h, occ, occ_scale)

   END SUBROUTINE dense_complex_kpoint_mermin_gradient

! **************************************************************************************************
!> \brief Check a non-diagonal occupation projector in the complex REF rotation gradient.
! **************************************************************************************************
   SUBROUTINE check_complex_projector_rotation_gradient(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 3
      REAL(KIND=dp), PARAMETER                           :: fd_step = 1.0E-6_dp

      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: base, occ_mat
      INTEGER                                            :: i, j
      REAL(KIND=dp)                                      :: err, fd, value_m, value_p
      REAL(KIND=dp), DIMENSION(n)                        :: eps_diag
      REAL(KIND=dp), DIMENSION(n, n)                     :: grad_im, grad_re, x_im, x_m_im, x_m_re, &
                                                            x_p_im, x_p_re, x_re

      base(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      base(1, 1) = CMPLX(-0.31_dp, 0.0_dp, KIND=dp)
      base(2, 2) = CMPLX(-0.04_dp, 0.0_dp, KIND=dp)
      base(3, 3) = CMPLX(0.16_dp, 0.0_dp, KIND=dp)
      base(1, 2) = CMPLX(0.05_dp, -0.02_dp, KIND=dp)
      base(1, 3) = CMPLX(-0.03_dp, 0.04_dp, KIND=dp)
      base(2, 3) = CMPLX(0.02_dp, 0.03_dp, KIND=dp)
      base(2, 1) = CONJG(base(1, 2))
      base(3, 1) = CONJG(base(1, 3))
      base(3, 2) = CONJG(base(2, 3))

      occ_mat(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      occ_mat(1, 1) = CMPLX(1.60_dp, 0.0_dp, KIND=dp)
      occ_mat(2, 2) = CMPLX(0.72_dp, 0.0_dp, KIND=dp)
      occ_mat(3, 3) = CMPLX(0.18_dp, 0.0_dp, KIND=dp)
      occ_mat(1, 2) = CMPLX(0.08_dp, -0.03_dp, KIND=dp)
      occ_mat(1, 3) = CMPLX(-0.04_dp, 0.02_dp, KIND=dp)
      occ_mat(2, 3) = CMPLX(0.03_dp, 0.05_dp, KIND=dp)
      occ_mat(2, 1) = CONJG(occ_mat(1, 2))
      occ_mat(3, 1) = CONJG(occ_mat(1, 3))
      occ_mat(3, 2) = CONJG(occ_mat(2, 3))

      x_re(:, :) = 0.0_dp
      x_re(1, 2) = 0.06_dp
      x_re(1, 3) = -0.03_dp
      x_re(2, 3) = 0.04_dp
      DO j = 2, n
         DO i = 1, j - 1
            x_re(j, i) = -x_re(i, j)
         END DO
      END DO

      x_im(:, :) = 0.0_dp
      x_im(1, 1) = 0.02_dp
      x_im(2, 2) = -0.04_dp
      x_im(3, 3) = 0.01_dp
      x_im(1, 2) = -0.02_dp
      x_im(1, 3) = 0.01_dp
      x_im(2, 3) = 0.03_dp
      DO j = 2, n
         DO i = 1, j - 1
            x_im(j, i) = x_im(i, j)
         END DO
      END DO

      CALL dense_projector_rotation_gradient(x_re, x_im, base, occ_mat, grad_re, grad_im, eps_diag)

      err = 0.0_dp
      DO j = 2, n
         DO i = 1, j - 1
            x_p_re(:, :) = x_re(:, :)
            x_m_re(:, :) = x_re(:, :)
            x_p_im(:, :) = x_im(:, :)
            x_m_im(:, :) = x_im(:, :)
            x_p_re(i, j) = x_p_re(i, j) + fd_step
            x_p_re(j, i) = x_p_re(j, i) - fd_step
            x_m_re(i, j) = x_m_re(i, j) - fd_step
            x_m_re(j, i) = x_m_re(j, i) + fd_step
            value_p = dense_projector_rotation_value(x_p_re, x_p_im, base, occ_mat)
            value_m = dense_projector_rotation_value(x_m_re, x_m_im, base, occ_mat)
            fd = (value_p - value_m)/(2.0_dp*fd_step)
            err = MAX(err, ABS(fd - grad_re(i, j)))

            x_p_re(:, :) = x_re(:, :)
            x_m_re(:, :) = x_re(:, :)
            x_p_im(:, :) = x_im(:, :)
            x_m_im(:, :) = x_im(:, :)
            x_p_im(i, j) = x_p_im(i, j) + fd_step
            x_p_im(j, i) = x_p_im(j, i) + fd_step
            x_m_im(i, j) = x_m_im(i, j) - fd_step
            x_m_im(j, i) = x_m_im(j, i) - fd_step
            value_p = dense_projector_rotation_value(x_p_re, x_p_im, base, occ_mat)
            value_m = dense_projector_rotation_value(x_m_re, x_m_im, base, occ_mat)
            fd = (value_p - value_m)/(2.0_dp*fd_step)
            err = MAX(err, ABS(fd - grad_im(i, j)))
         END DO
      END DO

      DO i = 1, n
         x_p_re(:, :) = x_re(:, :)
         x_m_re(:, :) = x_re(:, :)
         x_p_im(:, :) = x_im(:, :)
         x_m_im(:, :) = x_im(:, :)
         x_p_im(i, i) = x_p_im(i, i) + fd_step
         x_m_im(i, i) = x_m_im(i, i) - fd_step
         value_p = dense_projector_rotation_value(x_p_re, x_p_im, base, occ_mat)
         value_m = dense_projector_rotation_value(x_m_re, x_m_im, base, occ_mat)
         fd = (value_p - value_m)/(2.0_dp*fd_step)
         err = MAX(err, ABS(fd - 0.5_dp*grad_im(i, i)))
      END DO

      IF (err > tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,ES12.5)') TRIM(label), "failed with error", err
      ELSE
         WRITE (*, '(A,1X,A,1X,ES12.5)') TRIM(label), "OK with error", err
      END IF

   END SUBROUTINE check_complex_projector_rotation_gradient

! **************************************************************************************************
!> \brief Check that the finite rotation Hessian includes the moving occupation projector.
! **************************************************************************************************
   SUBROUTINE check_complex_occ_projector_hessian(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 3, nrot = n*(n - 1)
      REAL(KIND=dp), PARAMETER                           :: active_tol = 1.0E-7_dp, &
                                                            fd_step = 2.0E-5_dp

      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: base, chc, chc_m, chc_p, dchc, docc_diag, &
                                                            docc_mat, docc_offdiag, occ_mat, u
      INTEGER                                            :: i, j, r, s
      INTEGER, DIMENSION(nrot)                           :: var_i, var_j, var_kind
      REAL(KIND=dp)                                      :: active_diag_err, active_err, active_offdiag_err, err, fd, &
                                                            split_err, symmetry_err, value_m, value_p
      REAL(KIND=dp), DIMENSION(n)                        :: eps_diag, occ_scale
      REAL(KIND=dp), DIMENSION(nrot)                     :: grad_m_vec, grad_p_vec
      REAL(KIND=dp), DIMENSION(nrot, nrot)               :: hdiag, hfrechet, hfrozen, hoffdiag, hmove, hsplit
      REAL(KIND=dp), DIMENSION(n, n)                     :: grad_im, grad_re, x_im, x_m_im, x_m_re, &
                                                            x_p_im, x_p_re, x_re

      base(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      base(1, 1) = CMPLX(-0.37_dp, 0.0_dp, KIND=dp)
      base(2, 2) = CMPLX(-0.08_dp, 0.0_dp, KIND=dp)
      base(3, 3) = CMPLX(0.21_dp, 0.0_dp, KIND=dp)
      base(1, 2) = CMPLX(0.07_dp, -0.03_dp, KIND=dp)
      base(1, 3) = CMPLX(-0.04_dp, 0.05_dp, KIND=dp)
      base(2, 3) = CMPLX(0.03_dp, 0.04_dp, KIND=dp)
      base(2, 1) = CONJG(base(1, 2))
      base(3, 1) = CONJG(base(1, 3))
      base(3, 2) = CONJG(base(2, 3))

      occ_scale(:) = [1.54_dp, 0.83_dp, 0.25_dp]

      x_re(:, :) = 0.0_dp
      x_re(1, 2) = 0.05_dp
      x_re(1, 3) = -0.025_dp
      x_re(2, 3) = 0.035_dp
      DO j = 2, n
         DO i = 1, j - 1
            x_re(j, i) = -x_re(i, j)
         END DO
      END DO

      x_im(:, :) = 0.0_dp
      x_im(1, 1) = 0.015_dp
      x_im(2, 2) = -0.03_dp
      x_im(3, 3) = 0.02_dp
      x_im(1, 2) = -0.015_dp
      x_im(1, 3) = 0.02_dp
      x_im(2, 3) = 0.025_dp
      DO j = 2, n
         DO i = 1, j - 1
            x_im(j, i) = x_im(i, j)
         END DO
      END DO

      r = 0
      DO j = 2, n
         DO i = 1, j - 1
            r = r + 1
            var_i(r) = i
            var_j(r) = j
            var_kind(r) = 1
            r = r + 1
            var_i(r) = i
            var_j(r) = j
            var_kind(r) = 2
         END DO
      END DO

      CALL dense_spectral_projector_gradient(x_re, x_im, base, occ_scale, grad_re, grad_im, eps_diag)

      err = 0.0_dp
      hdiag(:, :) = 0.0_dp
      hoffdiag(:, :) = 0.0_dp
      DO s = 1, nrot
         x_p_re(:, :) = x_re(:, :)
         x_m_re(:, :) = x_re(:, :)
         x_p_im(:, :) = x_im(:, :)
         x_m_im(:, :) = x_im(:, :)
         CALL perturb_complex_rotation_var_2d(x_p_re, x_p_im, var_i(s), var_j(s), var_kind(s), fd_step)
         CALL perturb_complex_rotation_var_2d(x_m_re, x_m_im, var_i(s), var_j(s), var_kind(s), -fd_step)
         value_p = dense_spectral_projector_value(x_p_re, x_p_im, base, occ_scale)
         value_m = dense_spectral_projector_value(x_m_re, x_m_im, base, occ_scale)
         fd = (value_p - value_m)/(2.0_dp*fd_step)
         SELECT CASE (var_kind(s))
         CASE (1)
            err = MAX(err, ABS(fd - grad_re(var_i(s), var_j(s))))
         CASE (2)
            err = MAX(err, ABS(fd - grad_im(var_i(s), var_j(s))))
         END SELECT
      END DO

      CALL dense_rotation_u(x_re, x_im, u)
      chc(:, :) = MATMUL(CONJG(TRANSPOSE(u)), MATMUL(base, u))
      CALL dense_spectral_occ_projector(x_re, x_im, base, occ_scale, occ_mat)
      DO s = 1, nrot
         x_p_re(:, :) = x_re(:, :)
         x_m_re(:, :) = x_re(:, :)
         x_p_im(:, :) = x_im(:, :)
         x_m_im(:, :) = x_im(:, :)
         CALL perturb_complex_rotation_var_2d(x_p_re, x_p_im, var_i(s), var_j(s), var_kind(s), fd_step)
         CALL perturb_complex_rotation_var_2d(x_m_re, x_m_im, var_i(s), var_j(s), var_kind(s), -fd_step)

         CALL dense_spectral_projector_gradient(x_p_re, x_p_im, base, occ_scale, &
                                                grad_re, grad_im, eps_diag)
         CALL extract_rotation_gradient_vector_2d(grad_re, grad_im, var_i, var_j, var_kind, grad_p_vec)
         CALL dense_spectral_projector_gradient(x_m_re, x_m_im, base, occ_scale, &
                                                grad_re, grad_im, eps_diag)
         CALL extract_rotation_gradient_vector_2d(grad_re, grad_im, var_i, var_j, var_kind, grad_m_vec)
         hmove(:, s) = (grad_p_vec(:) - grad_m_vec(:))/(2.0_dp*fd_step)

         CALL dense_projector_rotation_gradient(x_p_re, x_p_im, base, occ_mat, grad_re, grad_im, eps_diag)
         CALL extract_rotation_gradient_vector_2d(grad_re, grad_im, var_i, var_j, var_kind, grad_p_vec)
         CALL dense_projector_rotation_gradient(x_m_re, x_m_im, base, occ_mat, grad_re, grad_im, eps_diag)
         CALL extract_rotation_gradient_vector_2d(grad_re, grad_im, var_i, var_j, var_kind, grad_m_vec)
         hfrozen(:, s) = (grad_p_vec(:) - grad_m_vec(:))/(2.0_dp*fd_step)

         CALL dense_rotation_u(x_p_re, x_p_im, u)
         chc_p(:, :) = MATMUL(CONJG(TRANSPOSE(u)), MATMUL(base, u))
         CALL dense_rotation_u(x_m_re, x_m_im, u)
         chc_m(:, :) = MATMUL(CONJG(TRANSPOSE(u)), MATMUL(base, u))
         dchc(:, :) = (chc_p(:, :) - chc_m(:, :))/(2.0_dp*fd_step)
         CALL dense_spectral_occ_projector_derivative_chc(chc, dchc, occ_scale, docc_mat)
         CALL dense_projector_rotation_gradient(x_re, x_im, base, docc_mat, grad_re, grad_im, eps_diag)
         CALL extract_rotation_gradient_vector_2d(grad_re, grad_im, var_i, var_j, var_kind, grad_p_vec)
         hfrechet(:, s) = hfrozen(:, s) + grad_p_vec(:)
         docc_diag(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
         docc_offdiag(:, :) = docc_mat(:, :)
         DO i = 1, n
            docc_diag(i, i) = docc_mat(i, i)
            docc_offdiag(i, i) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
         END DO
         CALL dense_projector_rotation_gradient(x_re, x_im, base, docc_diag, grad_re, grad_im, eps_diag)
         CALL extract_rotation_gradient_vector_2d(grad_re, grad_im, var_i, var_j, var_kind, grad_p_vec)
         hdiag(:, s) = grad_p_vec(:)
         CALL dense_projector_rotation_gradient(x_re, x_im, base, docc_offdiag, grad_re, grad_im, eps_diag)
         CALL extract_rotation_gradient_vector_2d(grad_re, grad_im, var_i, var_j, var_kind, grad_p_vec)
         hoffdiag(:, s) = grad_p_vec(:)
      END DO

      hsplit(:, :) = hfrozen(:, :) + hoffdiag(:, :)
      symmetry_err = MAXVAL(ABS(hmove(:, :) - TRANSPOSE(hmove(:, :))))
      active_err = MAXVAL(ABS(hmove(:, :) - hfrozen(:, :)))
      active_diag_err = MAXVAL(ABS(hdiag(:, :)))
      active_offdiag_err = MAXVAL(ABS(hoffdiag(:, :)))
      split_err = MAXVAL(ABS(hmove(:, :) - (hsplit(:, :) + hdiag(:, :))))
      err = MAX(err, symmetry_err)
      err = MAX(err, MAXVAL(ABS(hmove(:, :) - hfrechet(:, :))))
      err = MAX(err, split_err)

      IF (err > tol .OR. active_err < active_tol .OR. active_diag_err < active_tol .OR. &
          active_offdiag_err < active_tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), "failed with error/active/diag/offdiag", &
            err, active_err, active_diag_err, active_offdiag_err
      ELSE
         WRITE (*, '(A,1X,A,4(1X,ES12.5))') TRIM(label), "OK with error/active/diag/offdiag", &
            err, active_err, active_diag_err, active_offdiag_err
      END IF

   END SUBROUTINE check_complex_occ_projector_hessian

! **************************************************************************************************
!> \brief Separate projector-basis gauge noise from physical occupation response.
! **************************************************************************************************
   SUBROUTINE check_complex_occ_projector_gauge_separation(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 3
      REAL(KIND=dp), PARAMETER                           :: fd_step = 1.0E-8_dp, gap = 2.0E-5_dp, &
                                                            slope = -0.65_dp

      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: chc, dchc, docc_mat, fd_mat, occ_mat_m, &
                                                            occ_mat_p
      REAL(KIND=dp)                                      :: err, gauge_err
      REAL(KIND=dp), DIMENSION(n)                        :: occ_equal, occ_linear

      chc(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      chc(1, 1) = CMPLX(0.00_dp, 0.0_dp, KIND=dp)
      chc(2, 2) = CMPLX(gap, 0.0_dp, KIND=dp)
      chc(3, 3) = CMPLX(0.41_dp, 0.0_dp, KIND=dp)

      dchc(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      dchc(1, 2) = CMPLX(0.31_dp, -0.17_dp, KIND=dp)
      dchc(2, 1) = CONJG(dchc(1, 2))

      occ_equal(:) = [0.74_dp, 0.74_dp, 0.18_dp]
      CALL dense_spectral_occ_projector_derivative_chc(chc, dchc, occ_equal, docc_mat)
      CALL dense_spectral_occ_projector_from_chc(chc + fd_step*dchc, occ_equal, occ_mat_p)
      CALL dense_spectral_occ_projector_from_chc(chc - fd_step*dchc, occ_equal, occ_mat_m)
      fd_mat(:, :) = (occ_mat_p(:, :) - occ_mat_m(:, :))/(2.0_dp*fd_step)
      gauge_err = MAX(MAXVAL(ABS(docc_mat(:, :))), MAXVAL(ABS(fd_mat(:, :))))

      occ_linear(:) = occ_equal(:)
      occ_linear(2) = occ_linear(1) + slope*gap
      CALL dense_spectral_occ_projector_derivative_chc(chc, dchc, occ_linear, docc_mat)
      CALL dense_spectral_occ_projector_from_chc(chc + fd_step*dchc, occ_linear, occ_mat_p)
      CALL dense_spectral_occ_projector_from_chc(chc - fd_step*dchc, occ_linear, occ_mat_m)
      fd_mat(:, :) = (occ_mat_p(:, :) - occ_mat_m(:, :))/(2.0_dp*fd_step)
      err = MAXVAL(ABS(docc_mat(:, :) - fd_mat(:, :)))

      IF (err > tol .OR. gauge_err > tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,2(1X,ES12.5))') TRIM(label), "failed with error/gauge", err, gauge_err
      ELSE
         WRITE (*, '(A,1X,A,2(1X,ES12.5))') TRIM(label), "OK with error/gauge", err, gauge_err
      END IF

   END SUBROUTINE check_complex_occ_projector_gauge_separation

! **************************************************************************************************
!> \brief Dense value 1/2 sum_i f_i lambda_i[U^H H U].
! **************************************************************************************************
   FUNCTION dense_spectral_projector_value(x_re, x_im, base, occ_scale) RESULT(value)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: x_re, x_im
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: base
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: occ_scale
      REAL(KIND=dp)                                      :: value

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: chc, evec, u
      INTEGER                                            :: n
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: evals

      n = SIZE(x_re, 1)
      ALLOCATE (chc(n, n), evec(n, n), u(n, n), evals(n))
      CALL dense_rotation_u(x_re, x_im, u)
      chc(:, :) = MATMUL(CONJG(TRANSPOSE(u)), MATMUL(base, u))
      CALL diag_complex(chc, evec, evals)
      value = 0.5_dp*SUM(occ_scale(:)*evals(:))
      DEALLOCATE (chc, evec, u, evals)

   END FUNCTION dense_spectral_projector_value

! **************************************************************************************************
!> \brief Dense spectral occupation projector in the current rotated REF basis.
! **************************************************************************************************
   SUBROUTINE dense_spectral_occ_projector(x_re, x_im, base, occ_scale, occ_mat)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: x_re, x_im
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: base
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: occ_scale
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(OUT)     :: occ_mat

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: chc, u
      INTEGER                                            :: n

      n = SIZE(x_re, 1)
      ALLOCATE (chc(n, n), u(n, n))
      CALL dense_rotation_u(x_re, x_im, u)
      chc(:, :) = MATMUL(CONJG(TRANSPOSE(u)), MATMUL(base, u))
      CALL dense_spectral_occ_projector_from_chc(chc, occ_scale, occ_mat)
      DEALLOCATE (chc, u)

   END SUBROUTINE dense_spectral_occ_projector

! **************************************************************************************************
!> \brief Dense spectral occupation projector for a supplied Hermitian CHC matrix.
! **************************************************************************************************
   SUBROUTINE dense_spectral_occ_projector_from_chc(chc, occ_scale, occ_mat)
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: chc
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: occ_scale
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(OUT)     :: occ_mat

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: evec, weighted_evec
      INTEGER                                            :: i, n
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: evals

      n = SIZE(chc, 1)
      ALLOCATE (evec(n, n), weighted_evec(n, n), evals(n))
      CALL diag_complex(chc, evec, evals)
      weighted_evec(:, :) = evec(:, :)
      DO i = 1, n
         weighted_evec(:, i) = occ_scale(i)*weighted_evec(:, i)
      END DO
      occ_mat(:, :) = MATMUL(weighted_evec, CONJG(TRANSPOSE(evec)))
      DEALLOCATE (evec, weighted_evec, evals)

   END SUBROUTINE dense_spectral_occ_projector_from_chc

! **************************************************************************************************
!> \brief Fréchet derivative of the dense spectral occupation projector.
! **************************************************************************************************
   SUBROUTINE dense_spectral_occ_projector_derivative_chc(chc, dchc, occ_scale, docc_mat)
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: chc, dchc
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: occ_scale
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(OUT)     :: docc_mat

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: evec, kernel, work
      INTEGER                                            :: i, j, n
      REAL(KIND=dp)                                      :: coeff, denom, gap_tol, occ_tol, scale
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: evals

      n = SIZE(chc, 1)
      ALLOCATE (evec(n, n), kernel(n, n), work(n, n), evals(n))

      CALL diag_complex(chc, evec, evals)
      work(:, :) = MATMUL(CONJG(TRANSPOSE(evec)), MATMUL(dchc, evec))
      scale = MAX(1.0_dp, MAXVAL(ABS(evals(:))))
      gap_tol = SQRT(EPSILON(1.0_dp))*scale
      occ_tol = SQRT(EPSILON(1.0_dp))*MAX(1.0_dp, MAXVAL(ABS(occ_scale(:))))
      kernel(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      DO j = 1, n
         DO i = 1, n
            IF (i == j) CYCLE
            denom = evals(i) - evals(j)
            IF (ABS(denom) > gap_tol) THEN
               coeff = (occ_scale(i) - occ_scale(j))/denom
            ELSE IF (ABS(occ_scale(i) - occ_scale(j)) <= occ_tol) THEN
               coeff = 0.0_dp
            ELSE
               coeff = (occ_scale(i) - occ_scale(j))/SIGN(gap_tol, denom)
            END IF
            kernel(i, j) = coeff*work(i, j)
         END DO
      END DO

      docc_mat(:, :) = MATMUL(evec, MATMUL(kernel, CONJG(TRANSPOSE(evec))))
      docc_mat(:, :) = 0.5_dp*(docc_mat(:, :) + CONJG(TRANSPOSE(docc_mat(:, :))))
      DEALLOCATE (evec, kernel, work, evals)

   END SUBROUTINE dense_spectral_occ_projector_derivative_chc

! **************************************************************************************************
!> \brief Gradient of the dense spectral occupation-projector value.
! **************************************************************************************************
   SUBROUTINE dense_spectral_projector_gradient(x_re, x_im, base, occ_scale, grad_re, grad_im, eps_diag)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: x_re, x_im
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: base
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: occ_scale
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: grad_re, grad_im
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: eps_diag

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: occ_mat
      INTEGER                                            :: n

      n = SIZE(x_re, 1)
      ALLOCATE (occ_mat(n, n))
      CALL dense_spectral_occ_projector(x_re, x_im, base, occ_scale, occ_mat)
      CALL dense_projector_rotation_gradient(x_re, x_im, base, occ_mat, grad_re, grad_im, eps_diag)
      DEALLOCATE (occ_mat)

   END SUBROUTINE dense_spectral_projector_gradient

! **************************************************************************************************
!> \brief Dense value 1/2 Tr[P U^H H U] for a doubled Hermitian occupation projector.
! **************************************************************************************************
   FUNCTION dense_projector_rotation_value(x_re, x_im, base, occ_mat) RESULT(value)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: x_re, x_im
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: base, occ_mat
      REAL(KIND=dp)                                      :: value

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: chc, u
      INTEGER                                            :: n

      n = SIZE(x_re, 1)
      ALLOCATE (chc(n, n), u(n, n))
      CALL dense_rotation_u(x_re, x_im, u)
      chc(:, :) = MATMUL(CONJG(TRANSPOSE(u)), MATMUL(base, u))
      value = 0.5_dp*REAL(SUM(CONJG(occ_mat(:, :))*chc(:, :)), KIND=dp)
      DEALLOCATE (chc, u)

   END FUNCTION dense_projector_rotation_value

! **************************************************************************************************
!> \brief Dense gradient of 1/2 Tr[P U^H H U] for a doubled Hermitian occupation projector.
! **************************************************************************************************
   SUBROUTINE dense_projector_rotation_gradient(x_re, x_im, base, occ_mat, grad_re, grad_im, eps_diag)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: x_re, x_im
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: base, occ_mat
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: grad_re, grad_im
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: eps_diag

      COMPLEX(KIND=dp)                                   :: kernel
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: dedu, evec, frechet, generator, &
                                                            inner, outer, u, work
      INTEGER                                            :: i, j, n
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: evals

      n = SIZE(x_re, 1)
      ALLOCATE (generator(n, n), evec(n, n), u(n, n), work(n, n), dedu(n, n))
      ALLOCATE (inner(n, n), outer(n, n), frechet(n, n), evals(n))

      generator(:, :) = CMPLX(-x_im(:, :), x_re(:, :), KIND=dp)
      CALL diag_complex(generator, evec, evals)
      CALL dense_rotation_u_from_eigensystem(evec, evals, u)

      dedu(:, :) = MATMUL(base, MATMUL(u, occ_mat))

      work(:, :) = MATMUL(CONJG(TRANSPOSE(u)), MATMUL(base, u))
      DO i = 1, n
         eps_diag(i) = REAL(work(i, i), KIND=dp)
      END DO

      inner(:, :) = MATMUL(CONJG(TRANSPOSE(evec)), MATMUL(dedu, evec))
      DO j = 1, n
         DO i = 1, n
            kernel = CONJG(qs_ot_complex_exp_frechet_kernel(evals(i), evals(j)))
            outer(i, j) = inner(i, j)*kernel
         END DO
      END DO
      frechet(:, :) = MATMUL(evec, MATMUL(outer, CONJG(TRANSPOSE(evec))))

      grad_re(:, :) = REAL(frechet, KIND=dp) - TRANSPOSE(REAL(frechet, KIND=dp))
      grad_im(:, :) = AIMAG(frechet) + TRANSPOSE(AIMAG(frechet))

      DEALLOCATE (generator, evec, u, work, dedu, inner, outer, frechet, evals)

   END SUBROUTINE dense_projector_rotation_gradient

! **************************************************************************************************
!> \brief Check the complex k-point NONDIAG_ENERGY tether gradient.
! **************************************************************************************************
   SUBROUTINE check_complex_nondiag_tether_gradient(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 3
      REAL(KIND=dp), PARAMETER                           :: fd_step = 1.0E-6_dp, kweight = 0.75_dp, &
                                                            strength = 0.6_dp

      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: base
      INTEGER                                            :: i, j
      REAL(KIND=dp)                                      :: err, fd, value_m, value_p
      REAL(KIND=dp), DIMENSION(n)                        :: eps, grad_eps
      REAL(KIND=dp), DIMENSION(n, n)                     :: grad_im, grad_re, x_im, x_m_im, x_m_re, &
                                                            x_p_im, x_p_re, x_re

      base(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      base(1, 1) = CMPLX(-0.30_dp, 0.0_dp, KIND=dp)
      base(2, 2) = CMPLX(-0.02_dp, 0.0_dp, KIND=dp)
      base(3, 3) = CMPLX(0.12_dp, 0.0_dp, KIND=dp)
      base(1, 2) = CMPLX(0.04_dp, -0.03_dp, KIND=dp)
      base(2, 1) = CONJG(base(1, 2))
      base(1, 3) = CMPLX(-0.02_dp, 0.05_dp, KIND=dp)
      base(3, 1) = CONJG(base(1, 3))
      base(2, 3) = CMPLX(0.03_dp, 0.01_dp, KIND=dp)
      base(3, 2) = CONJG(base(2, 3))

      x_re(:, :) = 0.0_dp
      x_re(1, 2) = 0.08_dp
      x_re(2, 1) = -x_re(1, 2)
      x_re(1, 3) = -0.04_dp
      x_re(3, 1) = -x_re(1, 3)
      x_re(2, 3) = 0.02_dp
      x_re(3, 2) = -x_re(2, 3)

      x_im(:, :) = 0.0_dp
      x_im(1, 1) = 0.03_dp
      x_im(2, 2) = -0.05_dp
      x_im(3, 3) = 0.01_dp
      x_im(1, 2) = -0.02_dp
      x_im(2, 1) = x_im(1, 2)
      x_im(1, 3) = 0.01_dp
      x_im(3, 1) = x_im(1, 3)
      x_im(2, 3) = 0.04_dp
      x_im(3, 2) = x_im(2, 3)

      eps = [-0.24_dp, -0.01_dp, 0.09_dp]
      CALL dense_nondiag_tether_gradient(x_re, x_im, base, eps, strength, kweight, &
                                         grad_re, grad_im, grad_eps, value_p)

      err = 0.0_dp
      DO i = 1, n
         eps(i) = eps(i) + fd_step
         value_p = dense_nondiag_tether_value(x_re, x_im, base, eps, strength, kweight)
         eps(i) = eps(i) - 2.0_dp*fd_step
         value_m = dense_nondiag_tether_value(x_re, x_im, base, eps, strength, kweight)
         eps(i) = eps(i) + fd_step
         fd = (value_p - value_m)/(2.0_dp*fd_step)
         err = MAX(err, ABS(fd - grad_eps(i)))
      END DO

      DO j = 2, n
         DO i = 1, j - 1
            x_p_re(:, :) = x_re(:, :)
            x_m_re(:, :) = x_re(:, :)
            x_p_im(:, :) = x_im(:, :)
            x_m_im(:, :) = x_im(:, :)
            x_p_re(i, j) = x_p_re(i, j) + fd_step
            x_p_re(j, i) = x_p_re(j, i) - fd_step
            x_m_re(i, j) = x_m_re(i, j) - fd_step
            x_m_re(j, i) = x_m_re(j, i) + fd_step
            value_p = dense_nondiag_tether_value(x_p_re, x_p_im, base, eps, strength, kweight)
            value_m = dense_nondiag_tether_value(x_m_re, x_m_im, base, eps, strength, kweight)
            fd = (value_p - value_m)/(2.0_dp*fd_step)
            err = MAX(err, ABS(fd - grad_re(i, j)))

            x_p_re(:, :) = x_re(:, :)
            x_m_re(:, :) = x_re(:, :)
            x_p_im(:, :) = x_im(:, :)
            x_m_im(:, :) = x_im(:, :)
            x_p_im(i, j) = x_p_im(i, j) + fd_step
            x_p_im(j, i) = x_p_im(j, i) + fd_step
            x_m_im(i, j) = x_m_im(i, j) - fd_step
            x_m_im(j, i) = x_m_im(j, i) - fd_step
            value_p = dense_nondiag_tether_value(x_p_re, x_p_im, base, eps, strength, kweight)
            value_m = dense_nondiag_tether_value(x_m_re, x_m_im, base, eps, strength, kweight)
            fd = (value_p - value_m)/(2.0_dp*fd_step)
            err = MAX(err, ABS(fd - grad_im(i, j)))
         END DO
      END DO

      DO i = 1, n
         x_p_re(:, :) = x_re(:, :)
         x_m_re(:, :) = x_re(:, :)
         x_p_im(:, :) = x_im(:, :)
         x_m_im(:, :) = x_im(:, :)
         x_p_im(i, i) = x_p_im(i, i) + fd_step
         x_m_im(i, i) = x_m_im(i, i) - fd_step
         value_p = dense_nondiag_tether_value(x_p_re, x_p_im, base, eps, strength, kweight)
         value_m = dense_nondiag_tether_value(x_m_re, x_m_im, base, eps, strength, kweight)
         fd = (value_p - value_m)/(2.0_dp*fd_step)
         err = MAX(err, ABS(fd - 0.5_dp*grad_im(i, i)))
      END DO

      IF (err > tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,ES12.5)') TRIM(label), "failed with error", err
      ELSE
         WRITE (*, '(A,1X,A,1X,ES12.5)') TRIM(label), "OK with error", err
      END IF

   END SUBROUTINE check_complex_nondiag_tether_gradient

! **************************************************************************************************
!> \brief Check the gauge-projected complex canonical tether gradient.
! **************************************************************************************************
   SUBROUTINE check_complex_canonical_tether_gradient(label, tol, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: tol
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 3
      REAL(KIND=dp), PARAMETER                           :: fd_step = 1.0E-6_dp

      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: base
      INTEGER                                            :: i, j
      REAL(KIND=dp)                                      :: err, fd, value_m, value_p
      REAL(KIND=dp), DIMENSION(n)                        :: eps, grad_eps, weights
      REAL(KIND=dp), DIMENSION(n, n)                     :: grad_im, grad_re, x_im, x_m_im, x_m_re, &
                                                            x_p_im, x_p_re, x_re

      base(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      base(1, 1) = CMPLX(-0.31_dp, 0.0_dp, KIND=dp)
      base(2, 2) = CMPLX(-0.04_dp, 0.0_dp, KIND=dp)
      base(3, 3) = CMPLX(0.17_dp, 0.0_dp, KIND=dp)
      base(1, 2) = CMPLX(0.05_dp, -0.02_dp, KIND=dp)
      base(2, 1) = CONJG(base(1, 2))
      base(1, 3) = CMPLX(-0.03_dp, 0.04_dp, KIND=dp)
      base(3, 1) = CONJG(base(1, 3))
      base(2, 3) = CMPLX(0.02_dp, 0.03_dp, KIND=dp)
      base(3, 2) = CONJG(base(2, 3))

      x_re(:, :) = 0.0_dp
      x_re(1, 2) = 0.07_dp
      x_re(2, 1) = -x_re(1, 2)
      x_re(1, 3) = -0.03_dp
      x_re(3, 1) = -x_re(1, 3)
      x_re(2, 3) = 0.04_dp
      x_re(3, 2) = -x_re(2, 3)

      x_im(:, :) = 0.0_dp
      x_im(1, 1) = 0.02_dp
      x_im(2, 2) = -0.04_dp
      x_im(3, 3) = 0.03_dp
      x_im(1, 2) = -0.01_dp
      x_im(2, 1) = x_im(1, 2)
      x_im(1, 3) = 0.03_dp
      x_im(3, 1) = x_im(1, 3)
      x_im(2, 3) = 0.02_dp
      x_im(3, 2) = x_im(2, 3)

      eps = [-0.27_dp, -0.03_dp, 0.13_dp]
      weights = [0.42_dp, 0.73_dp, 0.28_dp]
      CALL dense_canonical_tether_gradient(x_re, x_im, base, eps, weights, &
                                           grad_re, grad_im, grad_eps, value_p)

      err = 0.0_dp
      DO i = 1, n
         eps(i) = eps(i) + fd_step
         value_p = dense_canonical_tether_value(x_re, x_im, base, eps, weights)
         eps(i) = eps(i) - 2.0_dp*fd_step
         value_m = dense_canonical_tether_value(x_re, x_im, base, eps, weights)
         eps(i) = eps(i) + fd_step
         fd = (value_p - value_m)/(2.0_dp*fd_step)
         err = MAX(err, ABS(fd - grad_eps(i)))
      END DO

      DO j = 2, n
         DO i = 1, j - 1
            x_p_re(:, :) = x_re(:, :)
            x_m_re(:, :) = x_re(:, :)
            x_p_im(:, :) = x_im(:, :)
            x_m_im(:, :) = x_im(:, :)
            x_p_re(i, j) = x_p_re(i, j) + fd_step
            x_p_re(j, i) = x_p_re(j, i) - fd_step
            x_m_re(i, j) = x_m_re(i, j) - fd_step
            x_m_re(j, i) = x_m_re(j, i) + fd_step
            value_p = dense_canonical_tether_value(x_p_re, x_p_im, base, eps, weights)
            value_m = dense_canonical_tether_value(x_m_re, x_m_im, base, eps, weights)
            fd = (value_p - value_m)/(2.0_dp*fd_step)
            err = MAX(err, ABS(fd - grad_re(i, j)))

            x_p_re(:, :) = x_re(:, :)
            x_m_re(:, :) = x_re(:, :)
            x_p_im(:, :) = x_im(:, :)
            x_m_im(:, :) = x_im(:, :)
            x_p_im(i, j) = x_p_im(i, j) + fd_step
            x_p_im(j, i) = x_p_im(j, i) + fd_step
            x_m_im(i, j) = x_m_im(i, j) - fd_step
            x_m_im(j, i) = x_m_im(j, i) - fd_step
            value_p = dense_canonical_tether_value(x_p_re, x_p_im, base, eps, weights)
            value_m = dense_canonical_tether_value(x_m_re, x_m_im, base, eps, weights)
            fd = (value_p - value_m)/(2.0_dp*fd_step)
            err = MAX(err, ABS(fd - grad_im(i, j)))
         END DO
      END DO

      DO i = 1, n
         x_p_re(:, :) = x_re(:, :)
         x_m_re(:, :) = x_re(:, :)
         x_p_im(:, :) = x_im(:, :)
         x_m_im(:, :) = x_im(:, :)
         x_p_im(i, i) = x_p_im(i, i) + fd_step
         x_m_im(i, i) = x_m_im(i, i) - fd_step
         value_p = dense_canonical_tether_value(x_p_re, x_p_im, base, eps, weights)
         value_m = dense_canonical_tether_value(x_m_re, x_m_im, base, eps, weights)
         fd = (value_p - value_m)/(2.0_dp*fd_step)
         err = MAX(err, ABS(fd - 0.5_dp*grad_im(i, i)))
      END DO

      IF (err > tol) THEN
         nfail = nfail + 1
         WRITE (*, '(A,1X,A,1X,ES12.5)') TRIM(label), "failed with error", err
      ELSE
         WRITE (*, '(A,1X,A,1X,ES12.5)') TRIM(label), "OK with error", err
      END IF

   END SUBROUTINE check_complex_canonical_tether_gradient

! **************************************************************************************************
!> \brief Dense gauge-projected complex canonical tether value.
! **************************************************************************************************
   FUNCTION dense_canonical_tether_value(x_re, x_im, base, eps, weights) RESULT(value)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: x_re, x_im
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: base
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: eps, weights
      REAL(KIND=dp)                                      :: value

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: chc, u
      INTEGER                                            :: i, n
      REAL(KIND=dp)                                      :: shift, weight_sum
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: delta

      n = SIZE(x_re, 1)
      ALLOCATE (chc(n, n), u(n, n), delta(n))
      CALL dense_rotation_u(x_re, x_im, u)
      chc(:, :) = MATMUL(CONJG(TRANSPOSE(u)), MATMUL(base, u))
      weight_sum = SUM(weights)
      shift = 0.0_dp
      DO i = 1, n
         delta(i) = REAL(chc(i, i), KIND=dp) - eps(i)
         shift = shift + weights(i)*delta(i)
      END DO
      IF (weight_sum > EPSILON(weight_sum)) shift = shift/weight_sum
      delta(:) = delta(:) - shift
      value = 0.5_dp*SUM(weights(:)*delta(:)**2)
      DEALLOCATE (chc, u, delta)

   END FUNCTION dense_canonical_tether_value

! **************************************************************************************************
!> \brief Dense gauge-projected complex canonical tether gradients.
! **************************************************************************************************
   SUBROUTINE dense_canonical_tether_gradient(x_re, x_im, base, eps, weights, &
                                              grad_re, grad_im, grad_eps, value)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: x_re, x_im
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: base
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: eps, weights
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: grad_re, grad_im
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: grad_eps
      REAL(KIND=dp), INTENT(OUT)                         :: value

      COMPLEX(KIND=dp)                                   :: kernel
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: chc, dedu, dmat, evec, frechet, &
                                                            generator, inner, outer, u
      INTEGER                                            :: i, j, n
      REAL(KIND=dp)                                      :: shift, weight_sum
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: delta, evals

      n = SIZE(x_re, 1)
      ALLOCATE (chc(n, n), dedu(n, n), dmat(n, n), evec(n, n), frechet(n, n), generator(n, n))
      ALLOCATE (inner(n, n), outer(n, n), u(n, n), delta(n), evals(n))

      generator(:, :) = CMPLX(-x_im(:, :), x_re(:, :), KIND=dp)
      CALL diag_complex(generator, evec, evals)
      CALL dense_rotation_u_from_eigensystem(evec, evals, u)

      chc(:, :) = MATMUL(CONJG(TRANSPOSE(u)), MATMUL(base, u))
      weight_sum = SUM(weights)
      shift = 0.0_dp
      DO i = 1, n
         delta(i) = REAL(chc(i, i), KIND=dp) - eps(i)
         shift = shift + weights(i)*delta(i)
      END DO
      IF (weight_sum > EPSILON(weight_sum)) shift = shift/weight_sum
      delta(:) = delta(:) - shift
      value = 0.5_dp*SUM(weights(:)*delta(:)**2)
      grad_eps(:) = -weights(:)*delta(:)

      dmat(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      DO i = 1, n
         dmat(i, i) = CMPLX(weights(i)*delta(i), 0.0_dp, KIND=dp)
      END DO
      dedu(:, :) = 2.0_dp*MATMUL(base, MATMUL(u, dmat))

      inner(:, :) = MATMUL(CONJG(TRANSPOSE(evec)), MATMUL(dedu, evec))
      DO j = 1, n
         DO i = 1, n
            kernel = CONJG(qs_ot_complex_exp_frechet_kernel(evals(i), evals(j)))
            outer(i, j) = inner(i, j)*kernel
         END DO
      END DO
      frechet(:, :) = MATMUL(evec, MATMUL(outer, CONJG(TRANSPOSE(evec))))

      grad_re(:, :) = REAL(frechet, KIND=dp) - TRANSPOSE(REAL(frechet, KIND=dp))
      grad_im(:, :) = AIMAG(frechet) + TRANSPOSE(AIMAG(frechet))

      DEALLOCATE (chc, dedu, dmat, evec, frechet, generator, inner, outer, u, delta, evals)

   END SUBROUTINE dense_canonical_tether_gradient

! **************************************************************************************************
!> \brief Dense complex NONDIAG_ENERGY tether value.
! **************************************************************************************************
   FUNCTION dense_nondiag_tether_value(x_re, x_im, base, eps, strength, kweight) RESULT(value)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: x_re, x_im
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: base
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: eps
      REAL(KIND=dp), INTENT(IN)                          :: strength, kweight
      REAL(KIND=dp)                                      :: value

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: dmat, u
      INTEGER                                            :: i, n

      n = SIZE(x_re, 1)
      ALLOCATE (dmat(n, n), u(n, n))
      CALL dense_rotation_u(x_re, x_im, u)
      dmat(:, :) = MATMUL(CONJG(TRANSPOSE(u)), MATMUL(base, u))
      DO i = 1, n
         dmat(i, i) = dmat(i, i) - eps(i)
      END DO
      value = 0.5_dp*strength*kweight*SUM(ABS(dmat(:, :))**2)
      DEALLOCATE (dmat, u)

   END FUNCTION dense_nondiag_tether_value

! **************************************************************************************************
!> \brief Dense complex NONDIAG_ENERGY tether gradients.
! **************************************************************************************************
   SUBROUTINE dense_nondiag_tether_gradient(x_re, x_im, base, eps, strength, kweight, &
                                            grad_re, grad_im, grad_eps, value)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: x_re, x_im
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: base
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: eps
      REAL(KIND=dp), INTENT(IN)                          :: strength, kweight
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: grad_re, grad_im
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: grad_eps
      REAL(KIND=dp), INTENT(OUT)                         :: value

      COMPLEX(KIND=dp)                                   :: kernel
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: dmat, dedu, evec, frechet, generator, &
                                                            inner, outer, u
      INTEGER                                            :: i, j, n
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: evals

      n = SIZE(x_re, 1)
      ALLOCATE (dmat(n, n), dedu(n, n), evec(n, n), frechet(n, n), generator(n, n))
      ALLOCATE (inner(n, n), outer(n, n), u(n, n), evals(n))

      generator(:, :) = CMPLX(-x_im(:, :), x_re(:, :), KIND=dp)
      CALL diag_complex(generator, evec, evals)
      CALL dense_rotation_u_from_eigensystem(evec, evals, u)

      dmat(:, :) = MATMUL(CONJG(TRANSPOSE(u)), MATMUL(base, u))
      DO i = 1, n
         dmat(i, i) = dmat(i, i) - eps(i)
      END DO
      value = 0.5_dp*strength*kweight*SUM(ABS(dmat(:, :))**2)
      DO i = 1, n
         grad_eps(i) = -strength*kweight*REAL(dmat(i, i), KIND=dp)
      END DO

      dedu(:, :) = 2.0_dp*strength*kweight*MATMUL(base, MATMUL(u, dmat))
      inner(:, :) = MATMUL(CONJG(TRANSPOSE(evec)), MATMUL(dedu, evec))
      DO j = 1, n
         DO i = 1, n
            kernel = CONJG(qs_ot_complex_exp_frechet_kernel(evals(i), evals(j)))
            outer(i, j) = inner(i, j)*kernel
         END DO
      END DO
      frechet(:, :) = MATMUL(evec, MATMUL(outer, CONJG(TRANSPOSE(evec))))

      grad_re(:, :) = REAL(frechet, KIND=dp) - TRANSPOSE(REAL(frechet, KIND=dp))
      grad_im(:, :) = AIMAG(frechet) + TRANSPOSE(AIMAG(frechet))

      DEALLOCATE (dmat, dedu, evec, frechet, generator, inner, outer, u, evals)

   END SUBROUTINE dense_nondiag_tether_gradient

! **************************************************************************************************
!> \brief Dense gradient of sum_i occ_scale_i [U^H H U]_ii.
! **************************************************************************************************
   SUBROUTINE dense_ref_rotation_gradient(x_re, x_im, base, occ_scale, grad_re, grad_im, eps_diag)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: x_re, x_im
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: base
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: occ_scale
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: grad_re, grad_im
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: eps_diag

      COMPLEX(KIND=dp)                                   :: kernel
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: dedu, evec, frechet, generator, &
                                                            inner, outer, u, work
      INTEGER                                            :: i, j, n
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: evals

      n = SIZE(x_re, 1)
      ALLOCATE (generator(n, n), evec(n, n), u(n, n), work(n, n), dedu(n, n))
      ALLOCATE (inner(n, n), outer(n, n), frechet(n, n), evals(n))

      generator(:, :) = CMPLX(-x_im(:, :), x_re(:, :), KIND=dp)
      CALL diag_complex(generator, evec, evals)
      CALL dense_rotation_u_from_eigensystem(evec, evals, u)

      work(:, :) = MATMUL(base, u)
      dedu(:, :) = work(:, :)
      DO j = 1, n
         dedu(:, j) = occ_scale(j)*dedu(:, j)
      END DO

      work(:, :) = MATMUL(CONJG(TRANSPOSE(u)), MATMUL(base, u))
      DO i = 1, n
         eps_diag(i) = REAL(work(i, i), KIND=dp)
      END DO

      inner(:, :) = MATMUL(CONJG(TRANSPOSE(evec)), MATMUL(dedu, evec))
      DO j = 1, n
         DO i = 1, n
            kernel = CONJG(qs_ot_complex_exp_frechet_kernel(evals(i), evals(j)))
            outer(i, j) = inner(i, j)*kernel
         END DO
      END DO
      frechet(:, :) = MATMUL(evec, MATMUL(outer, CONJG(TRANSPOSE(evec))))

      grad_re(:, :) = REAL(frechet, KIND=dp) - TRANSPOSE(REAL(frechet, KIND=dp))
      grad_im(:, :) = AIMAG(frechet) + TRANSPOSE(AIMAG(frechet))

      DEALLOCATE (generator, evec, u, work, dedu, inner, outer, frechet, evals)

   END SUBROUTINE dense_ref_rotation_gradient

! **************************************************************************************************
!> \brief Build U=exp(-i X) for the dense anti-Hermitian k-point generator convention.
! **************************************************************************************************
   SUBROUTINE dense_rotation_u(x_re, x_im, u)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: x_re, x_im
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(OUT)     :: u

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: evec, generator
      INTEGER                                            :: n
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: evals

      n = SIZE(x_re, 1)
      ALLOCATE (evec(n, n), generator(n, n), evals(n))
      generator(:, :) = CMPLX(-x_im(:, :), x_re(:, :), KIND=dp)
      CALL diag_complex(generator, evec, evals)
      CALL dense_rotation_u_from_eigensystem(evec, evals, u)
      DEALLOCATE (evec, generator, evals)

   END SUBROUTINE dense_rotation_u

! **************************************************************************************************
!> \brief Build U=V exp(-i evals) V^H.
! **************************************************************************************************
   SUBROUTINE dense_rotation_u_from_eigensystem(evec, evals, u)
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: evec
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: evals
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(OUT)     :: u

      INTEGER                                            :: j

      u(:, :) = evec(:, :)
      DO j = 1, SIZE(evals)
         u(:, j) = EXP(CMPLX(0.0_dp, -evals(j), KIND=dp))*u(:, j)
      END DO
      u(:, :) = MATMUL(u, CONJG(TRANSPOSE(evec)))

   END SUBROUTINE dense_rotation_u_from_eigensystem

! **************************************************************************************************
!> \brief Small local Fermi-Dirac k-point smearing helper for this standalone unit test.
! **************************************************************************************************
   SUBROUTINE dense_fermi_kpoint(occ, mu, kTS, eps, nel, wk, sigma, maxocc)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: occ
      REAL(KIND=dp), INTENT(OUT)                         :: mu, kTS
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: eps
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: wk
      REAL(KIND=dp), INTENT(IN)                          :: nel, sigma, maxocc

      REAL(KIND=dp), PARAMETER                           :: epsocc = 1.0E-12_dp

      INTEGER                                            :: iter
      REAL(KIND=dp)                                      :: de, mu_max, mu_min, n_now

      de = MAX(sigma*LOG((1.0_dp - epsocc)/epsocc), 0.5_dp)
      mu_min = MINVAL(eps) - de
      mu_max = MAXVAL(eps) + de
      DO iter = 1, 400
         mu = 0.5_dp*(mu_min + mu_max)
         CALL dense_fermi_kpoint_eval(occ, n_now, kTS, eps, mu, wk, sigma, maxocc)
         IF (ABS(n_now - nel) < nel*epsocc) EXIT
         IF (n_now <= nel) THEN
            mu_min = mu
         ELSE
            mu_max = mu
         END IF
      END DO
      mu = 0.5_dp*(mu_min + mu_max)
      CALL dense_fermi_kpoint_eval(occ, n_now, kTS, eps, mu, wk, sigma, maxocc)

   END SUBROUTINE dense_fermi_kpoint

! **************************************************************************************************
!> \brief Evaluate Fermi-Dirac occupations and -TS correction for a fixed chemical potential.
! **************************************************************************************************
   SUBROUTINE dense_fermi_kpoint_eval(occ, nel, kTS, eps, mu, wk, sigma, maxocc)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: occ
      REAL(KIND=dp), INTENT(OUT)                         :: nel, kTS
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: eps
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: wk
      REAL(KIND=dp), INTENT(IN)                          :: mu, sigma, maxocc

      INTEGER                                            :: ik, imo
      REAL(KIND=dp)                                      :: arg, term1, term2, tmp, tmp2, tmp3, &
                                                            tmp4, tmplog

      nel = 0.0_dp
      kTS = 0.0_dp
      DO ik = 1, SIZE(eps, 2)
         DO imo = 1, SIZE(eps, 1)
            IF (eps(imo, ik) > mu) THEN
               arg = -(eps(imo, ik) - mu)/sigma
               tmp = EXP(arg)
               tmp4 = tmp + 1.0_dp
               tmp2 = tmp/tmp4
               tmp3 = 1.0_dp/tmp4
               tmplog = -LOG(tmp4)
               term1 = tmp2*(arg + tmplog)
               term2 = tmp3*tmplog
            ELSE
               arg = (eps(imo, ik) - mu)/sigma
               tmp = EXP(arg)
               tmp4 = tmp + 1.0_dp
               tmp2 = 1.0_dp/tmp4
               tmp3 = tmp/tmp4
               tmplog = -LOG(tmp4)
               term1 = tmp2*tmplog
               term2 = tmp3*(arg + tmplog)
            END IF
            occ(imo, ik) = maxocc*tmp2
            nel = nel + wk(ik)*occ(imo, ik)
            kTS = kTS + sigma*maxocc*(term1 + term2)*wk(ik)
         END DO
      END DO

   END SUBROUTINE dense_fermi_kpoint_eval

! **************************************************************************************************
!> \brief Small dense 2x2 matrix exponential via Taylor series.
! **************************************************************************************************
   SUBROUTINE dense_expm2(a, expa)
      COMPLEX(KIND=dp), DIMENSION(2, 2), INTENT(IN)      :: a
      COMPLEX(KIND=dp), DIMENSION(2, 2), INTENT(OUT)     :: expa

      COMPLEX(KIND=dp), DIMENSION(2, 2)                  :: term
      INTEGER                                            :: i, n

      expa(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      term(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      DO i = 1, 2
         expa(i, i) = CMPLX(1.0_dp, 0.0_dp, KIND=dp)
         term(i, i) = CMPLX(1.0_dp, 0.0_dp, KIND=dp)
      END DO

      DO n = 1, 80
         term = MATMUL(term, a)/REAL(n, KIND=dp)
         expa = expa + term
         IF (MAXVAL(ABS(term)) < 1.0E-20_dp) EXIT
      END DO

   END SUBROUTINE dense_expm2

END PROGRAM qs_ot_frechet_unittest
