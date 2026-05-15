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
   SUBROUTINE apply_dense_frechet(evals, z, out)
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: evals
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: z
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(OUT)     :: out

      INTEGER                                            :: i, j

      DO i = 1, SIZE(z, 1)
      DO j = 1, SIZE(z, 2)
         out(i, j) = z(i, j)*qs_ot_complex_exp_frechet_kernel(evals(i), evals(j))
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
