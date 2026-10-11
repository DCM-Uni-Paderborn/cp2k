!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Coupled inverse reconstruction on explicit physical lattice images.
!>        K_R=delta_R0 I+C_R. Values outside the caller's solution domain are zero, without
!>        wrapping image indices. The same convolution and boundary define its response and
!>        adjoint. A complete finite-support residual includes the exterior equations.
! **************************************************************************************************
MODULE upf_inverse_lattice
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: twopi
   USE upf_inverse,                     ONLY: upf_inverse_init,&
                                              upf_inverse_project,&
                                              upf_inverse_type
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   TYPE upf_inverse_lattice_type
      TYPE(upf_inverse_type) :: diagonal
      INTEGER, ALLOCATABLE :: images(:, :), kernel_images(:, :), source_index(:, :)
      COMPLEX(dp), ALLOCATABLE :: coupling(:, :, :)
   END TYPE upf_inverse_lattice_type

   PUBLIC :: upf_inverse_lattice_type, upf_inverse_lattice_init, upf_inverse_lattice_project
   PUBLIC :: upf_inverse_lattice_response, upf_inverse_lattice_adjoint
   PUBLIC :: upf_inverse_lattice_residual, upf_inverse_lattice_bloch
   PUBLIC :: upf_inverse_lattice_response_residual
   PUBLIC :: upf_inverse_lattice_bound, upf_inverse_lattice_error_bounds

CONTAINS

! **************************************************************************************************
!> \brief Prepare the sparse convolution and its identical-cell LU preconditioner.
!>        Storage scales with the kernel and retained images, without a dense image-pair matrix.
!> \param kernel_images Distinct physical integer triples Q, including only supplied C_Q.
!> \param coupling Original non-Hermitian projector matrices C_Q.
!> \param images Distinct solution images R. Their order is arbitrary.
!> \param lattice New owned convolution with zero exterior coefficients.
!> \param info Original one-cell LU status. Failed preconditioners must not be used.
! **************************************************************************************************
   SUBROUTINE upf_inverse_lattice_init(kernel_images, coupling, images, lattice, info)
      INTEGER, INTENT(IN)                                :: kernel_images(:, :)
      COMPLEX(dp), INTENT(IN)                            :: coupling(:, :, :)
      INTEGER, INTENT(IN)                                :: images(:, :)
      TYPE(upf_inverse_lattice_type), INTENT(OUT)        :: lattice
      INTEGER, INTENT(OUT)                               :: info

      COMPLEX(dp), ALLOCATABLE                           :: diagonal(:, :)
      INTEGER                                            :: i, np, q, zero

      CALL check_images(kernel_images)
      CALL check_images(images)
      CPASSERT(SIZE(images, 2) > 0)
      np = SIZE(coupling, 1)
      CPASSERT(ALL(SHAPE(coupling) == [np, np, SIZE(kernel_images, 2)]))
      info = 0
      IF (.NOT. ALL(ieee_is_finite(REAL(coupling, dp))) .OR. .NOT. ALL(ieee_is_finite(AIMAG(coupling)))) THEN
         info = np + 2
         RETURN
      END IF
      ALLOCATE (lattice%images, SOURCE=images)
      ALLOCATE (lattice%kernel_images, SOURCE=kernel_images)
      ALLOCATE (lattice%coupling, SOURCE=coupling)
      ALLOCATE (lattice%source_index(SIZE(kernel_images, 2), SIZE(images, 2)))
      DO i = 1, SIZE(images, 2)
         DO q = 1, SIZE(kernel_images, 2)
            lattice%source_index(q, i) = find_image(images, images(:, i) - kernel_images(:, q))
         END DO
      END DO
      ALLOCATE (diagonal(np, np), SOURCE=CMPLX(0.0_dp, 0.0_dp, dp))
      zero = find_image(kernel_images, [0, 0, 0])
      IF (zero > 0) diagonal(:, :) = coupling(:, :, zero)
      CALL upf_inverse_init(diagonal, lattice%diagonal, info)
   END SUBROUTINE upf_inverse_lattice_init

! **************************************************************************************************
!> \brief Solve the coupled physical-image equations using restarted complex GMRES.
!>        Convergence uses the true unpreconditioned residual for each Gaussian column.
!>        A converged interior solve does not establish convergence of the image truncation.
!> \param lattice ...
!> \param overlaps B_R in the retained image order.
!> \param coefficients W_R, with the same projector, Gaussian and image order.
!> \param tolerance Relative Euclidean residual tolerance per right-hand side.
!> \param info Zero on success, one on exhaustion, two on Krylov breakdown.
!> \param max_iterations Optional maximum iterations per Gaussian column.
! **************************************************************************************************
   SUBROUTINE upf_inverse_lattice_project(lattice, overlaps, coefficients, tolerance, info, max_iterations)
      TYPE(upf_inverse_lattice_type), INTENT(IN)         :: lattice
      COMPLEX(dp), INTENT(IN)                            :: overlaps(:, :, :)
      COMPLEX(dp), INTENT(OUT)                           :: coefficients(:, :, :)
      REAL(dp), INTENT(IN)                               :: tolerance
      INTEGER, INTENT(OUT)                               :: info
      INTEGER, INTENT(IN), OPTIONAL                      :: max_iterations

      CALL solve(lattice, overlaps, coefficients, tolerance, info, .FALSE., max_iterations)
   END SUBROUTINE upf_inverse_lattice_project

! **************************************************************************************************
!> \brief Differentiate K*W=B on the unchanged physical domain, including all C_Q responses.
!> \param lattice ...
!> \param coefficients ...
!> \param overlaps_response Original dB_R, including its defining geometry and measure response.
!> \param coupling_response Original dC_Q in the kernel image order.
!> \param response ...
!> \param tolerance ...
!> \param info ...
! **************************************************************************************************
   SUBROUTINE upf_inverse_lattice_response(lattice, coefficients, overlaps_response, &
                                           coupling_response, response, tolerance, info)
      TYPE(upf_inverse_lattice_type), INTENT(IN)         :: lattice
      COMPLEX(dp), INTENT(IN)                            :: coefficients(:, :, :), &
                                                            overlaps_response(:, :, :), &
                                                            coupling_response(:, :, :)
      COMPLEX(dp), INTENT(OUT)                           :: response(:, :, :)
      REAL(dp), INTENT(IN)                               :: tolerance
      INTEGER, INTENT(OUT)                               :: info

      COMPLEX(dp), ALLOCATABLE                           :: rhs(:, :, :)
      INTEGER                                            :: i, q, source

      CPASSERT(ALL(SHAPE(coupling_response) == SHAPE(lattice%coupling)))
      CPASSERT(ALL(SHAPE(overlaps_response) == SHAPE(coefficients)))
      ALLOCATE (rhs, SOURCE=overlaps_response)
      DO i = 1, SIZE(lattice%images, 2)
         DO q = 1, SIZE(lattice%kernel_images, 2)
            source = lattice%source_index(q, i)
            IF (source == 0) CYCLE
            rhs(:, :, i) = rhs(:, :, i) - MATMUL(coupling_response(:, :, q), coefficients(:, :, source))
         END DO
      END DO
      CALL solve(lattice, rhs, response, tolerance, info, .FALSE.)
   END SUBROUTINE upf_inverse_lattice_response

! **************************************************************************************************
!> \brief Real energy adjoint of the same coupled solve, with Re sum(conjg(G)*dW) pairing.
!>        Shared C_Q blocks receive the sum from every physical image pair. No Hermitian,
!>        one-center, mesh-weight or occupation assumption is made.
!> \param lattice ...
!> \param coefficients ...
!> \param covector ...
!> \param overlaps_covector ...
!> \param coupling_covector ...
!> \param tolerance ...
!> \param info ...
! **************************************************************************************************
   SUBROUTINE upf_inverse_lattice_adjoint(lattice, coefficients, covector, overlaps_covector, &
                                          coupling_covector, tolerance, info)
      TYPE(upf_inverse_lattice_type), INTENT(IN)         :: lattice
      COMPLEX(dp), INTENT(IN)                            :: coefficients(:, :, :), covector(:, :, :)
      COMPLEX(dp), INTENT(OUT)                           :: overlaps_covector(:, :, :), &
                                                            coupling_covector(:, :, :)
      REAL(dp), INTENT(IN)                               :: tolerance
      INTEGER, INTENT(OUT)                               :: info

      INTEGER                                            :: i, q, source

      CPASSERT(ALL(SHAPE(coefficients) == SHAPE(covector)))
      CPASSERT(ALL(SHAPE(coupling_covector) == SHAPE(lattice%coupling)))
      CALL solve(lattice, covector, overlaps_covector, tolerance, info, .TRUE.)
      coupling_covector = CMPLX(0.0_dp, 0.0_dp, dp)
      IF (info /= 0) RETURN
      DO i = 1, SIZE(lattice%images, 2)
         DO q = 1, SIZE(lattice%kernel_images, 2)
            source = lattice%source_index(q, i)
            IF (source == 0) CYCLE
            coupling_covector(:, :, q) = coupling_covector(:, :, q) - &
                                         MATMUL(overlaps_covector(:, :, i), CONJG(TRANSPOSE(coefficients(:, :, source))))
         END DO
      END DO
   END SUBROUTINE upf_inverse_lattice_adjoint

! **************************************************************************************************
!> \brief Compute K*W-B on the full support of the supplied finite sequences, including exterior.
!>        Gaussian B_R can have infinite range. Its omitted tail must be bounded separately.
!>        The residual is not an inverse error bound without a bound on the infinite inverse.
!> \param lattice ...
!> \param coefficients Retained W_R, zero outside lattice%images.
!> \param overlap_images Independent complete support of the supplied finite B sequence.
!> \param overlaps ...
!> \param residual_images New union of W support, B support and all Q+R images.
!> \param residual New unscaled residual in this image order.
! **************************************************************************************************
   SUBROUTINE upf_inverse_lattice_residual(lattice, coefficients, overlap_images, overlaps, residual_images, residual)
      TYPE(upf_inverse_lattice_type), INTENT(IN)         :: lattice
      COMPLEX(dp), INTENT(IN)                           :: coefficients(:, :, :), overlaps(:, :, :)
      INTEGER, INTENT(IN)                                :: overlap_images(:, :)
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: residual_images(:, :)
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: residual(:, :, :)

      INTEGER                                            :: i, image(3), n, np, nr, q, TARGET
      INTEGER, ALLOCATABLE                               :: union_images(:, :)

      CALL check_images(overlap_images)
      np = lattice%diagonal%nproj
      nr = SIZE(lattice%images, 2)
      CPASSERT(ALL(SHAPE(coefficients) == [np, SIZE(overlaps, 2), nr]))
      CPASSERT(SIZE(overlaps, 1) == np)
      n = SIZE(overlap_images, 2)
      CPASSERT(SIZE(overlaps, 3) == n)
      ALLOCATE (union_images(3, nr*(1 + SIZE(lattice%kernel_images, 2)) + SIZE(overlap_images, 2)))
      n = nr
      union_images(:, :n) = lattice%images
      DO i = 1, SIZE(overlap_images, 2)
         CALL append_image(overlap_images(:, i))
      END DO
      DO i = 1, nr
         DO q = 1, SIZE(lattice%kernel_images, 2)
            image = lattice%images(:, i) + lattice%kernel_images(:, q)
            CALL append_image(image)
         END DO
      END DO
      ALLOCATE (residual_images, SOURCE=union_images(:, :n))
      ALLOCATE (residual(np, SIZE(overlaps, 2), n), SOURCE=CMPLX(0.0_dp, 0.0_dp, dp))
      residual(:, :, :nr) = coefficients
      DO i = 1, SIZE(overlap_images, 2)
         TARGET = find_image(residual_images, overlap_images(:, i))
         residual(:, :, TARGET) = residual(:, :, TARGET) - overlaps(:, :, i)
      END DO
      DO i = 1, nr
         DO q = 1, SIZE(lattice%kernel_images, 2)
            TARGET = find_image(residual_images, lattice%images(:, i) + lattice%kernel_images(:, q))
            residual(:, :, TARGET) = residual(:, :, TARGET) + MATMUL(lattice%coupling(:, :, q), coefficients(:, :, i))
         END DO
      END DO

   CONTAINS

      SUBROUTINE append_image(new_image)
      INTEGER, INTENT(IN)                                :: new_image(3)

         IF (find_image(union_images(:, :n), new_image) > 0) RETURN
         n = n + 1
         union_images(:, n) = new_image
      END SUBROUTINE append_image
   END SUBROUTINE upf_inverse_lattice_residual

! **************************************************************************************************
!> \brief Full interior and exterior residual K*dW+dC*W-dB of the coupled geometry response.
!>        A small forward residual alone does not control the derivative's omitted image tail.
!> \param lattice ...
!> \param coefficients ...
!> \param response ...
!> \param overlap_images Original support of the supplied finite dB sequence.
!> \param overlaps_response ...
!> \param coupling_response ...
!> \param residual_images ...
!> \param residual ...
! **************************************************************************************************
   SUBROUTINE upf_inverse_lattice_response_residual(lattice, coefficients, response, overlap_images, &
                                                    overlaps_response, coupling_response, residual_images, residual)
      TYPE(upf_inverse_lattice_type), INTENT(IN)         :: lattice
      COMPLEX(dp), INTENT(IN)                            :: coefficients(:, :, :), response(:, :, :), &
                                                            overlaps_response(:, :, :), coupling_response(:, :, :)
      INTEGER, INTENT(IN)                                :: overlap_images(:, :)
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: residual_images(:, :)
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: residual(:, :, :)

      INTEGER                                            :: i, q, TARGET

      CPASSERT(ALL(SHAPE(coefficients) == SHAPE(response)))
      CPASSERT(ALL(SHAPE(coupling_response) == SHAPE(lattice%coupling)))
      CALL upf_inverse_lattice_residual(lattice, response, overlap_images, overlaps_response, residual_images, residual)
      DO i = 1, SIZE(lattice%images, 2)
         DO q = 1, SIZE(lattice%kernel_images, 2)
            TARGET = find_image(residual_images, lattice%images(:, i) + lattice%kernel_images(:, q))
            residual(:, :, TARGET) = residual(:, :, TARGET) + MATMUL(coupling_response(:, :, q), coefficients(:, :, i))
         END DO
      END DO
   END SUBROUTINE upf_inverse_lattice_response_residual

! **************************************************************************************************
!> \brief Sufficient lower bound on the smallest singular value of K(k) for every Bloch k.
!>        For the computed M=(I+C_0)^(-1), bound I-K(k)M by the sum of Frobenius norms.
!>        The measured zero-cell inverse residual and a floating-point margin are retained.
!>        Failure of this sufficient condition does not establish singularity of K(k).
!> \param lattice Original finite convolution, including every supplied nonzero image.
!> \param lower_bound Positive uniform bound on success, otherwise zero.
!> \param info Zero on success, one for an inconclusive bound, two for an unavailable inverse.
! **************************************************************************************************
   SUBROUTINE upf_inverse_lattice_bound(lattice, lower_bound, info)
      TYPE(upf_inverse_lattice_type), INTENT(IN)         :: lattice
      REAL(dp), INTENT(OUT)                              :: lower_bound
      INTEGER, INTENT(OUT)                               :: info

      COMPLEX(dp), ALLOCATABLE                           :: diagonal(:, :), identity(:, :), &
                                                            inverse(:, :), PRODUCT(:, :)
      INTEGER                                            :: i, np, q, zero
      REAL(dp)                                           :: defect, inverse_norm, kernel_norm, &
                                                            roundoff

      lower_bound = 0.0_dp
      info = 2
      IF (.NOT. lattice%diagonal%ready) RETURN
      np = lattice%diagonal%nproj
      IF (np == 0) THEN
         lower_bound = 1.0_dp
         info = 0
         RETURN
      END IF
      ALLOCATE (identity(np, np), SOURCE=CMPLX(0.0_dp, 0.0_dp, dp))
      ALLOCATE (inverse(np, np), PRODUCT(np, np), diagonal(np, np))
      DO i = 1, np
         identity(i, i) = 1.0_dp
      END DO
      diagonal(:, :) = identity
      zero = find_image(lattice%kernel_images, [0, 0, 0])
      IF (zero > 0) diagonal(:, :) = diagonal + lattice%coupling(:, :, zero)
      CALL upf_inverse_project(lattice%diagonal, identity, inverse)
      inverse_norm = matrix_norm(inverse)
      PRODUCT(:, :) = identity - MATMUL(diagonal, inverse)
      defect = matrix_norm(product)
      kernel_norm = matrix_norm(diagonal)
      DO q = 1, SIZE(lattice%kernel_images, 2)
         IF (q == zero) CYCLE
         PRODUCT(:, :) = MATMUL(lattice%coupling(:, :, q), inverse)
         defect = defect + matrix_norm(product)
         kernel_norm = kernel_norm + matrix_norm(lattice%coupling(:, :, q))
      END DO
      roundoff = 32.0_dp*EPSILON(0.0_dp)*REAL(np, dp)* &
                 REAL(MAX(1, SIZE(lattice%kernel_images, 2)), dp)
      defect = defect + roundoff*(1.0_dp + kernel_norm*inverse_norm)
      info = 1
      IF (.NOT. ieee_is_finite(defect) .OR. .NOT. ieee_is_finite(inverse_norm)) RETURN
      IF (defect >= 1.0_dp .OR. inverse_norm <= 0.0_dp) RETURN
      lower_bound = (1.0_dp - defect)/inverse_norm
      info = 0
   END SUBROUTINE upf_inverse_lattice_bound

! **************************************************************************************************
!> \brief Bound omitted coefficients and their Bloch sums using the full exterior residual.
!>        The l2 sequence error uses the l2 residual. A uniform pointwise Bloch error instead
!>        uses its l1 norm. Responses also retain dC times the forward coefficient error.
!>        These bounds concern the supplied finite kernel. Omitted kernel data and quadrature
!>        errors require separate control. Margins account for ordinary arithmetic rounding,
!>        without providing interval-arithmetic certification.
!> \param lattice ...
!> \param coefficients Finite W, zero outside the retained domain.
!> \param overlap_images Complete support of the supplied B sequence, independent of W.
!> \param overlaps ...
!> \param coefficient_error Absolute l2 Frobenius bound for the full coefficient sequence.
!> \param bloch_error Absolute Frobenius bound for W(k) at every Bloch k.
!> \param lower_bound Uniform smallest-singular-value bound of the supplied K(k).
!> \param info Zero on success. Otherwise errors remain HUGE and no bound is established.
!> \param overlap_tail Optional proven omitted B bounds, ordered [l2,l1]. Absent means zero.
!> \param response Optional finite dW. All response arrays and error outputs must be supplied.
!> \param overlaps_response Optional dB in the overlap image order.
!> \param coupling_response Optional original dC in the kernel image order.
!> \param response_error Optional absolute l2 Frobenius bound for the full dW sequence.
!> \param response_bloch_error Optional uniform pointwise Bloch bound for dW(k).
!> \param response_overlap_tail Optional proven omitted dB bounds, ordered [l2,l1].
! **************************************************************************************************
   SUBROUTINE upf_inverse_lattice_error_bounds(lattice, coefficients, overlap_images, overlaps, &
                                               coefficient_error, bloch_error, lower_bound, info, overlap_tail, &
                                               response, overlaps_response, coupling_response, response_error, &
                                               response_bloch_error, response_overlap_tail)
      TYPE(upf_inverse_lattice_type), INTENT(IN)         :: lattice
      COMPLEX(dp), INTENT(IN)                            :: coefficients(:, :, :)
      INTEGER, INTENT(IN)                                :: overlap_images(:, :)
      COMPLEX(dp), INTENT(IN)                            :: overlaps(:, :, :)
      REAL(dp), INTENT(OUT)                              :: coefficient_error, bloch_error, &
                                                            lower_bound
      INTEGER, INTENT(OUT)                               :: info
      REAL(dp), INTENT(IN), OPTIONAL                     :: overlap_tail(2)
      COMPLEX(dp), INTENT(IN), OPTIONAL                  :: response(:, :, :), &
                                                            overlaps_response(:, :, :), &
                                                            coupling_response(:, :, :)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: response_error, response_bloch_error
      REAL(dp), INTENT(IN), OPTIONAL                     :: response_overlap_tail(2)

      COMPLEX(dp), ALLOCATABLE                           :: residual(:, :, :)
      INTEGER, ALLOCATABLE                               :: residual_images(:, :)
      LOGICAL                                            :: have_response
      REAL(dp)                                           :: b_norm, c_norm, db_norm, dc_norm, &
                                                            dw_norm, norm1, norm2, roundoff, w_norm

      have_response = PRESENT(response) .OR. PRESENT(overlaps_response) .OR. PRESENT(coupling_response) .OR. &
                      PRESENT(response_error) .OR. PRESENT(response_bloch_error) .OR. PRESENT(response_overlap_tail)
      IF (have_response) THEN
         IF (.NOT. (PRESENT(response) .AND. PRESENT(overlaps_response) .AND. PRESENT(coupling_response) .AND. &
                    PRESENT(response_error) .AND. PRESENT(response_bloch_error))) THEN
            CPABORT("Incomplete lattice response bound arguments")
         END IF
         response_error = HUGE(0.0_dp)
         response_bloch_error = HUGE(0.0_dp)
      END IF
      IF (PRESENT(overlap_tail)) THEN
         IF (.NOT. ALL(ieee_is_finite(overlap_tail)) .OR. ANY(overlap_tail < 0.0_dp)) THEN
            CPABORT("Invalid overlap tail bound")
         END IF
      END IF
      IF (PRESENT(response_overlap_tail)) THEN
         IF (.NOT. ALL(ieee_is_finite(response_overlap_tail)) .OR. ANY(response_overlap_tail < 0.0_dp)) THEN
            CPABORT("Invalid response overlap tail bound")
         END IF
      END IF
      coefficient_error = HUGE(0.0_dp)
      bloch_error = HUGE(0.0_dp)
      lower_bound = 0.0_dp
      info = 2
      IF (.NOT. ALL(ieee_is_finite(REAL(coefficients, dp))) .OR. &
          .NOT. ALL(ieee_is_finite(AIMAG(coefficients))) .OR. &
          .NOT. ALL(ieee_is_finite(REAL(overlaps, dp))) .OR. .NOT. ALL(ieee_is_finite(AIMAG(overlaps)))) RETURN
      IF (have_response) THEN
         IF (.NOT. ALL(ieee_is_finite(REAL(response, dp))) .OR. .NOT. ALL(ieee_is_finite(AIMAG(response))) .OR. &
             .NOT. ALL(ieee_is_finite(REAL(overlaps_response, dp))) .OR. &
             .NOT. ALL(ieee_is_finite(AIMAG(overlaps_response))) .OR. &
             .NOT. ALL(ieee_is_finite(REAL(coupling_response, dp))) .OR. &
             .NOT. ALL(ieee_is_finite(AIMAG(coupling_response)))) RETURN
      END IF
      CALL upf_inverse_lattice_bound(lattice, lower_bound, info)
      IF (info /= 0) RETURN
      CALL sequence_norms(lattice%coupling, norm2, c_norm)
      CALL sequence_norms(coefficients, norm2, w_norm)
      CALL sequence_norms(overlaps, norm2, b_norm)
      CALL upf_inverse_lattice_residual(lattice, coefficients, overlap_images, overlaps, residual_images, residual)
      CALL sequence_norms(residual, norm2, norm1)
      roundoff = 32.0_dp*EPSILON(0.0_dp)*REAL(MAX(1, lattice%diagonal%nproj), dp)* &
                 REAL(MAX(1, SIZE(lattice%kernel_images, 2) + SIZE(lattice%images, 2) + SIZE(overlap_images, 2)), dp)
      norm2 = norm2 + roundoff*((1.0_dp + c_norm)*w_norm + b_norm)
      norm1 = norm1 + roundoff*((1.0_dp + c_norm)*w_norm + b_norm)
      IF (PRESENT(overlap_tail)) THEN
         norm2 = norm2 + overlap_tail(1)
         norm1 = norm1 + overlap_tail(2)
      END IF
      coefficient_error = norm2/lower_bound
      bloch_error = norm1/lower_bound
      IF (have_response) THEN
         CALL sequence_norms(response, norm2, dw_norm)
         CALL sequence_norms(overlaps_response, norm2, db_norm)
         CALL sequence_norms(coupling_response, norm2, dc_norm)
         CALL upf_inverse_lattice_response_residual(lattice, coefficients, response, overlap_images, &
                                                    overlaps_response, coupling_response, residual_images, residual)
         CALL sequence_norms(residual, norm2, norm1)
         norm2 = norm2 + roundoff*((1.0_dp + c_norm)*dw_norm + dc_norm*w_norm + db_norm)
         norm1 = norm1 + roundoff*((1.0_dp + c_norm)*dw_norm + dc_norm*w_norm + db_norm)
         IF (PRESENT(response_overlap_tail)) THEN
            norm2 = norm2 + response_overlap_tail(1)
            norm1 = norm1 + response_overlap_tail(2)
         END IF
         response_error = (norm2 + dc_norm*coefficient_error)/lower_bound
         response_bloch_error = (norm1 + dc_norm*bloch_error)/lower_bound
      END IF
   END SUBROUTINE upf_inverse_lattice_error_bounds

! **************************************************************************************************
!> \brief Stable Frobenius norm and l1/l2 matrix-sequence norms, without squared large entries.
! **************************************************************************************************
   FUNCTION matrix_norm(matrix) RESULT(norm)
      COMPLEX(dp), INTENT(IN)                            :: matrix(:, :)
      REAL(dp)                                           :: norm

      INTEGER                                            :: column
      REAL(dp), EXTERNAL                                 :: dznrm2

      norm = 0.0_dp
      DO column = 1, SIZE(matrix, 2)
         norm = HYPOT(norm, dznrm2(SIZE(matrix, 1), matrix(:, column), 1))
      END DO
   END FUNCTION matrix_norm

   SUBROUTINE sequence_norms(SEQUENCE, norm2, norm1)
      COMPLEX(dp), INTENT(IN)                            :: sequence(:, :, :)
      REAL(dp), INTENT(OUT)                              :: norm2, norm1

      INTEGER                                            :: image
      REAL(dp)                                           :: norm

      norm2 = 0.0_dp
      norm1 = 0.0_dp
      DO image = 1, SIZE(SEQUENCE, 3)
         norm = matrix_norm(SEQUENCE(:, :, image))
         norm2 = HYPOT(norm2, norm)
         norm1 = norm1 + norm
      END DO
   END SUBROUTINE sequence_norms

! **************************************************************************************************
!> \brief Fourier sum of physical coefficients, without reciprocal mesh weights or image aliasing.
!> \param images ...
!> \param coefficients ...
!> \param reduced_k ...
!> \param bloch ...
! **************************************************************************************************
   SUBROUTINE upf_inverse_lattice_bloch(images, coefficients, reduced_k, bloch)
      INTEGER, INTENT(IN)                                :: images(:, :)
      COMPLEX(dp), INTENT(IN)                            :: coefficients(:, :, :)
      REAL(dp), INTENT(IN)                               :: reduced_k(3)
      COMPLEX(dp), INTENT(OUT)                           :: bloch(:, :)

      COMPLEX(dp)                                        :: phase
      INTEGER                                            :: i, nao, np, nr

      np = SIZE(coefficients, 1)
      nao = SIZE(coefficients, 2)
      nr = SIZE(coefficients, 3)
      CPASSERT(SIZE(images, 1) == 3)
      CPASSERT(SIZE(images, 2) == nr)
      CPASSERT(SIZE(bloch, 1) == np)
      CPASSERT(SIZE(bloch, 2) == nao)
      bloch = CMPLX(0.0_dp, 0.0_dp, dp)
      DO i = 1, SIZE(images, 2)
         phase = EXP(CMPLX(0.0_dp, twopi*DOT_PRODUCT(reduced_k, REAL(images(:, i), dp)), dp))
         bloch = bloch + phase*coefficients(:, :, i)
      END DO
   END SUBROUTINE upf_inverse_lattice_bloch

! **************************************************************************************************
!> \brief Columnwise left-preconditioned GMRES with double modified Gram-Schmidt.
!> \param lattice ...
!> \param rhs ...
!> \param solution ...
!> \param tolerance ...
!> \param info ...
!> \param adjoint ...
!> \param max_iterations ...
! **************************************************************************************************
   SUBROUTINE solve(lattice, rhs, solution, tolerance, info, adjoint, max_iterations)
      TYPE(upf_inverse_lattice_type), INTENT(IN)         :: lattice
      COMPLEX(dp), INTENT(IN)                            :: rhs(:, :, :)
      COMPLEX(dp), INTENT(OUT)                           :: solution(:, :, :)
      REAL(dp), INTENT(IN)                               :: tolerance
      INTEGER, INTENT(OUT)                               :: info
      LOGICAL, INTENT(IN)                                :: adjoint
      INTEGER, INTENT(IN), OPTIONAL                      :: max_iterations

      COMPLEX(dp)                                        :: product, rotated
      COMPLEX(dp), ALLOCATABLE                           :: b(:), g(:), h(:, :), r(:), s(:), &
                                                            v(:, :), x(:), y(:), z(:)
      INTEGER                                            :: column, i, iteration, j, limit, m, n, &
                                                            pass, used
      REAL(dp)                                           :: beta, bound, estimate_bound, norm_b, &
                                                            norm_z
      REAL(dp), ALLOCATABLE                              :: c(:)
      REAL(dp), EXTERNAL                                 :: dznrm2

      EXTERNAL zlartg

      CPASSERT(lattice%diagonal%ready)
      CPASSERT(tolerance > 0.0_dp .AND. tolerance < 1.0_dp)
      CPASSERT(SIZE(rhs, 1) == lattice%diagonal%nproj)
      CPASSERT(SIZE(rhs, 3) == SIZE(lattice%images, 2))
      CPASSERT(ALL(SHAPE(solution) == SHAPE(rhs)))
      info = 0
      solution = CMPLX(0.0_dp, 0.0_dp, dp)
      n = SIZE(rhs, 1)*SIZE(rhs, 3)
      IF (n == 0 .OR. SIZE(rhs, 2) == 0) RETURN
      limit = 500
      IF (PRESENT(max_iterations)) limit = max_iterations
      CPASSERT(limit > 0)
      m = MIN(30, n, limit)
      ALLOCATE (b(n), r(n), x(n), z(n), v(n, m + 1), h(m + 1, m), c(m), s(m), g(m + 1), y(m))
      DO column = 1, SIZE(rhs, 2)
         b(:) = RESHAPE(rhs(:, column, :), [n])
         norm_b = dznrm2(n, b, 1)
         IF (norm_b == 0.0_dp) CYCLE
         bound = tolerance*norm_b
         z(:) = b
         CALL precondition(lattice, z, adjoint)
         estimate_bound = tolerance*dznrm2(n, z, 1)
         x = CMPLX(0.0_dp, 0.0_dp, dp)
         iteration = 0
         DO
            CALL apply(lattice, x, r, adjoint)
            r(:) = b - r
            IF (dznrm2(n, r, 1) <= bound) EXIT
            IF (iteration >= limit) THEN
               info = 1
               EXIT
            END IF
            CALL precondition(lattice, r, adjoint)
            beta = dznrm2(n, r, 1)
            CPASSERT(beta > 0.0_dp)
            v(:, 1) = r/beta
            h = CMPLX(0.0_dp, 0.0_dp, dp)
            g = CMPLX(0.0_dp, 0.0_dp, dp)
            g(1) = CMPLX(beta, 0.0_dp, dp)
            used = 0
            DO j = 1, MIN(m, limit - iteration)
               iteration = iteration + 1
               CALL apply(lattice, v(:, j), z, adjoint)
               CALL precondition(lattice, z, adjoint)
               DO pass = 1, 2
                  DO i = 1, j
                     product = DOT_PRODUCT(v(:, i), z)
                     h(i, j) = h(i, j) + product
                     z(:) = z - product*v(:, i)
                  END DO
               END DO
               norm_z = dznrm2(n, z, 1)
               h(j + 1, j) = CMPLX(norm_z, 0.0_dp, dp)
               IF (norm_z > 0.0_dp) v(:, j + 1) = z/norm_z
               DO i = 1, j - 1
                  rotated = c(i)*h(i, j) + s(i)*h(i + 1, j)
                  h(i + 1, j) = -CONJG(s(i))*h(i, j) + c(i)*h(i + 1, j)
                  h(i, j) = rotated
               END DO
               CALL zlartg(h(j, j), h(j + 1, j), c(j), s(j), rotated)
               h(j, j) = rotated
               h(j + 1, j) = CMPLX(0.0_dp, 0.0_dp, dp)
               g(j + 1) = -CONJG(s(j))*g(j)
               g(j) = c(j)*g(j)
               IF (ABS(h(j, j)) == 0.0_dp) THEN
                  info = 2
                  EXIT
               END IF
               used = j
               IF (ABS(g(j + 1)) <= estimate_bound .OR. norm_z == 0.0_dp) EXIT
            END DO
            y(:used) = g(:used)
            DO i = used, 1, -1
               y(i) = (y(i) - SUM(h(i, i + 1:used)*y(i + 1:used)))/h(i, i)
            END DO
            x(:) = x + MATMUL(v(:, :used), y(:used))
            IF (info /= 0) EXIT
         END DO
         solution(:, column, :) = RESHAPE(x, [SIZE(rhs, 1), SIZE(rhs, 3)])
         IF (info /= 0) RETURN
      END DO
   END SUBROUTINE solve

! **************************************************************************************************
!> \brief Apply K or its exact conjugate transpose without wrapping physical indices.
!> \param lattice ...
!> \param vector ...
!> \param result ...
!> \param adjoint ...
! **************************************************************************************************
   SUBROUTINE apply(lattice, vector, RESULT, adjoint)
      TYPE(upf_inverse_lattice_type), INTENT(IN)         :: lattice
      COMPLEX(dp), INTENT(IN)                            :: vector(:)
      COMPLEX(dp), INTENT(OUT)                           :: result(:)
      LOGICAL, INTENT(IN)                                :: adjoint

      INTEGER                                            :: first, i, j, last, np, q, source_first, &
                                                            source_last

      EXTERNAL zgemv

      np = lattice%diagonal%nproj
      RESULT = vector
      IF (np == 0) RETURN
      DO i = 1, SIZE(lattice%images, 2)
         first = (i - 1)*np + 1
         last = i*np
         DO q = 1, SIZE(lattice%kernel_images, 2)
            j = lattice%source_index(q, i)
            IF (j == 0) CYCLE
            source_first = (j - 1)*np + 1
            source_last = j*np
            IF (adjoint) THEN
               CALL zgemv('C', np, np, CMPLX(1.0_dp, 0.0_dp, dp), lattice%coupling(:, :, q), np, &
                          vector(first:last), 1, CMPLX(1.0_dp, 0.0_dp, dp), RESULT(source_first:source_last), 1)
            ELSE
               CALL zgemv('N', np, np, CMPLX(1.0_dp, 0.0_dp, dp), lattice%coupling(:, :, q), np, &
                          vector(source_first:source_last), 1, CMPLX(1.0_dp, 0.0_dp, dp), RESULT(first:last), 1)
            END IF
         END DO
      END DO
   END SUBROUTINE apply

! **************************************************************************************************
!> \brief Reuse the one-cell original LU for every image, including the adjoint preconditioner.
!> \param lattice ...
!> \param vector ...
!> \param adjoint ...
! **************************************************************************************************
   SUBROUTINE precondition(lattice, vector, adjoint)
      TYPE(upf_inverse_lattice_type), INTENT(IN)         :: lattice
      COMPLEX(dp), INTENT(INOUT)                         :: vector(:)
      LOGICAL, INTENT(IN)                                :: adjoint

      CHARACTER                                          :: transpose
      INTEGER                                            :: info, np

      EXTERNAL zgetrs

      np = lattice%diagonal%nproj
      transpose = 'N'
      IF (adjoint) transpose = 'C'
      CALL zgetrs(transpose, np, SIZE(lattice%images, 2), lattice%diagonal%lu, np, &
                  lattice%diagonal%pivots, vector, np, info)
      CPASSERT(info == 0)
   END SUBROUTINE precondition

! **************************************************************************************************
!> \brief Check the physical image representation without imposing a rectangular mesh.
!> \param images ...
! **************************************************************************************************
   SUBROUTINE check_images(images)
      INTEGER, INTENT(IN)                                :: images(:, :)

      INTEGER                                            :: i

      CPASSERT(SIZE(images, 1) == 3)
      DO i = 1, SIZE(images, 2)
         CPASSERT(find_image(images(:, :i - 1), images(:, i)) == 0)
      END DO
   END SUBROUTINE check_images

! **************************************************************************************************
!> \brief Locate an exact physical image, returning zero for an absent triple.
!> \param images ...
!> \param image ...
!> \return ...
! **************************************************************************************************
   INTEGER FUNCTION find_image(images, image) RESULT(index)
      INTEGER, INTENT(IN)                                :: images(:, :), image(3)

      INTEGER                                            :: i

      index = 0
      DO i = 1, SIZE(images, 2)
         IF (.NOT. ALL(images(:, i) == image)) CYCLE
         index = i
         RETURN
      END DO
   END FUNCTION find_image
END MODULE upf_inverse_lattice
