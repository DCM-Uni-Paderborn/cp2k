!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Inverse of a specified finite-rank UPF reconstruction, including intercenter coupling.
!>        Delta contains original AE-minus-PS waves, P the original projector bras.
!>        For K=I+P*Delta and B=P*chi, W=K^-1*B gives U*chi=chi-Delta*W and P*U*chi=W.
!>        Complex arrays retain Bloch phases; no atomic block-diagonal approximation is made.
!>        This algebra neither supplies missing reference waves nor eliminates operator residuals.
! **************************************************************************************************
MODULE upf_inverse
   USE atom_upf,                        ONLY: atom_upfpot_type
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: dp
   USE upf_projector_integrals,         ONLY: upf_projector_init,&
                                              upf_projector_type
   USE upf_radial_fields,               ONLY: upf_radial_fields_init,&
                                              upf_radial_fields_type
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   TYPE upf_inverse_type
      INTEGER :: nproj = 0
      LOGICAL :: ready = .FALSE.
      REAL(KIND=dp) :: rcond = 0.0_dp
      COMPLEX(KIND=dp), ALLOCATABLE :: lu(:, :)
      INTEGER, ALLOCATABLE :: pivots(:)
   END TYPE upf_inverse_type

   PUBLIC :: upf_inverse_type, upf_inverse_init, upf_inverse_project, upf_inverse_project_response
   PUBLIC :: upf_inverse_project_adjoint
   PUBLIC :: upf_inverse_basis, upf_inverse_basis_response
   PUBLIC :: upf_inverse_operator, upf_inverse_operator_response
   PUBLIC :: upf_inverse_delta_init, upf_inverse_fields_init

CONTAINS

! **************************************************************************************************
!> \brief Prepare original AE-minus-PS radial overlaps in the native projector/harmonic order.
!>        The ordinary radial/Gaussian integrator supplies overlaps and AO-center derivatives.
!>        This is an overlap cache, not a replacement for the original nonlocal operator.
!> \param pot Original data with a complete full-wave pair for every beta projector.
!> \param degree Required Cartesian degree, including derivatives or moments.
!> \param delta Original radial quadrature and angular maps for (u_AE-u_PS)/r.
! **************************************************************************************************
   SUBROUTINE upf_inverse_delta_init(pot, degree, delta)
      TYPE(atom_upfpot_type), INTENT(IN)                 :: pot
      INTEGER, INTENT(IN)                                :: degree
      TYPE(upf_projector_type), INTENT(OUT)              :: delta

      TYPE(atom_upfpot_type)                             :: radial

      CALL inverse_radial_data(pot, radial)
      CALL upf_projector_init(radial, degree, delta)
   END SUBROUTINE upf_inverse_delta_init

! **************************************************************************************************
!> \brief Prepare original beta and reconstruction fields for common intercenter quadratures.
!>        Native radial channels and real harmonics agree with the Gaussian overlap cache.
!>        NC channels have zero reconstruction but retain their original projector bras.
!> \param pot Original UPF data, with complete full waves for reconstructing channels.
!> \param beta Original projector fields, without weights.
!> \param delta Original reconstruction differences, without weights.
! **************************************************************************************************
   SUBROUTINE upf_inverse_fields_init(pot, beta, delta)
      TYPE(atom_upfpot_type), INTENT(IN)                 :: pot
      TYPE(upf_radial_fields_type), INTENT(OUT)          :: beta, delta

      TYPE(atom_upfpot_type)                             :: radial

      CALL upf_radial_fields_init(pot%r, pot%beta, pot%lbeta, beta)
      IF (TRIM(pot%pseudo_type) == "NC" .AND. .NOT. pot%is_ultrasoft .AND. .NOT. pot%is_paw) THEN
         IF (pot%has_so) CPABORT("Spin-orbit UPF inverse fields are not supported")
         ALLOCATE (radial%beta(SIZE(pot%r), SIZE(pot%lbeta)))
         radial%beta = 0.0_dp
         CALL upf_radial_fields_init(pot%r, radial%beta, pot%lbeta, delta)
      ELSE
         CALL inverse_radial_data(pot, radial)
         CALL upf_radial_fields_init(radial%r, radial%beta, radial%lbeta, delta)
      END IF
   END SUBROUTINE upf_inverse_fields_init

! **************************************************************************************************
!> \brief Assemble original reconstruction differences for overlap and point-field caches.
!> \param pot Original UPF data.
!> \param radial Unweighted original differences in the native radial-channel order.
! **************************************************************************************************
   SUBROUTINE inverse_radial_data(pot, radial)
      TYPE(atom_upfpot_type), INTENT(IN)                 :: pot
      TYPE(atom_upfpot_type), INTENT(OUT)                :: radial

      INTEGER                                            :: i, n, nb

      IF (pot%has_so) THEN
         CPABORT("Spin-orbit UPF inverse fields are not supported")
      END IF
      IF (.NOT. ALLOCATED(pot%partial_waves)) THEN
         CPABORT("Inverse UPF requires original AE/PS full waves for every projector")
      END IF
      CPASSERT(ALLOCATED(pot%r) .AND. ALLOCATED(pot%rab))
      CPASSERT(ALLOCATED(pot%lbeta))
      IF (.NOT. (ALLOCATED(pot%partial_waves%ae_wfc) .AND. ALLOCATED(pot%partial_waves%ps_wfc))) THEN
         CPABORT("Incomplete original AE/PS full waves")
      END IF
      CPASSERT(ALLOCATED(pot%partial_waves%angular_momenta))
      n = SIZE(pot%r)
      nb = pot%number_of_proj
      CPASSERT(n >= 3 .AND. nb > 0 .AND. ALL(pot%r >= 0.0_dp))
      CPASSERT(SIZE(pot%rab) == n .AND. SIZE(pot%lbeta) == nb)
      CPASSERT(ALL(SHAPE(pot%partial_waves%ae_wfc) == [n, nb]))
      CPASSERT(ALL(SHAPE(pot%partial_waves%ps_wfc) == [n, nb]))
      CPASSERT(ALL(pot%partial_waves%angular_momenta == pot%lbeta))
      CPASSERT(ALL(ieee_is_finite(pot%partial_waves%ae_wfc)))
      CPASSERT(ALL(ieee_is_finite(pot%partial_waves%ps_wfc)))
      radial%number_of_proj = nb
      radial%has_so = pot%has_so
      ALLOCATE (radial%r, SOURCE=pot%r)
      ALLOCATE (radial%rab, SOURCE=pot%rab)
      ALLOCATE (radial%lbeta, SOURCE=pot%lbeta)
      ALLOCATE (radial%beta(n, nb), radial%dion(nb, nb))
      radial%dion(:, :) = 0.0_dp
      DO i = 1, nb
         radial%beta(:, i) = pot%partial_waves%ae_wfc(:, i) - pot%partial_waves%ps_wfc(:, i)
         IF (pot%r(1) == 0.0_dp) THEN
            radial%beta(2:, i) = radial%beta(2:, i)/pot%r(2:)
            ! Zero quadrature weight at r=0; this storage does not define a field continuation.
            radial%beta(1, i) = 0.0_dp
         ELSE
            radial%beta(:, i) = radial%beta(:, i)/pot%r
         END IF
      END DO
   END SUBROUTINE inverse_radial_data

! **************************************************************************************************
!> \brief Factor K with original projector normalization and all supplied intercenter blocks.
!> \param beta_delta P*Delta, including periodic images and Bloch phases.
!> \param inverse New factors; a failed construction cannot reuse an earlier factorization.
!> \param info Zero on success, LAPACK pivot on singularity, nproj+1 for poor conditioning,
!>             nproj+2 for nonfinite data. Failed factors must not be used.
!> \param minimum_rcond Optional condition limit; no regularization changes the reconstruction.
! **************************************************************************************************
   SUBROUTINE upf_inverse_init(beta_delta, inverse, info, minimum_rcond)
      COMPLEX(KIND=dp), INTENT(IN)                       :: beta_delta(:, :)
      TYPE(upf_inverse_type), INTENT(OUT)                :: inverse
      INTEGER, INTENT(OUT)                               :: info
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: minimum_rcond

      COMPLEX(KIND=dp), ALLOCATABLE                      :: work(:)
      INTEGER                                            :: i, n
      REAL(KIND=dp)                                      :: anorm, limit
      REAL(KIND=dp), ALLOCATABLE                         :: rwork(:)
      REAL(KIND=dp), EXTERNAL                            :: zlange

      EXTERNAL zgecon, zgetrf

      n = SIZE(beta_delta, 1)
      CPASSERT(SIZE(beta_delta, 2) == n)
      limit = 64.0_dp*EPSILON(0.0_dp)
      IF (PRESENT(minimum_rcond)) limit = minimum_rcond
      CPASSERT(ieee_is_finite(limit) .AND. limit >= 0.0_dp .AND. limit <= 1.0_dp)
      inverse%nproj = n
      info = 0
      IF (.NOT. ALL(ieee_is_finite(REAL(beta_delta, dp))) .OR. &
          .NOT. ALL(ieee_is_finite(AIMAG(beta_delta)))) THEN
         info = n + 2
         RETURN
      END IF
      ALLOCATE (inverse%lu(n, n), inverse%pivots(n))
      inverse%lu(:, :) = beta_delta
      DO i = 1, n
         inverse%lu(i, i) = inverse%lu(i, i) + 1.0_dp
      END DO
      IF (n == 0) THEN
         inverse%rcond = 1.0_dp
         inverse%ready = .TRUE.
         RETURN
      END IF
      ALLOCATE (work(2*n), rwork(2*n))
      anorm = zlange('1', n, n, inverse%lu, n, rwork)
      CALL zgetrf(n, n, inverse%lu, n, inverse%pivots, info)
      IF (info /= 0) RETURN
      CALL zgecon('1', n, inverse%lu, n, anorm, inverse%rcond, work, rwork, info)
      IF (info /= 0) RETURN
      IF (.NOT. ieee_is_finite(inverse%rcond)) THEN
         info = n + 2
         RETURN
      END IF
      IF (inverse%rcond < limit) THEN
         info = n + 1
         RETURN
      END IF
      inverse%ready = .TRUE.
   END SUBROUTINE upf_inverse_init

! **************************************************************************************************
!> \brief Solve for W=P*U*chi; do not replace original q by an orthogonalized metric.
!> \param inverse ...
!> \param beta_gaussian B=P*chi, with one Gaussian function per column.
!> \param coefficients W, using the already factored full K.
! **************************************************************************************************
   SUBROUTINE upf_inverse_project(inverse, beta_gaussian, coefficients)
      TYPE(upf_inverse_type), INTENT(IN)                 :: inverse
      COMPLEX(KIND=dp), INTENT(IN)                       :: beta_gaussian(:, :)
      COMPLEX(KIND=dp), INTENT(OUT)                      :: coefficients(:, :)

      INTEGER                                            :: info, n, nrhs

      EXTERNAL zgetrs

      CPASSERT(inverse%ready)
      n = inverse%nproj
      nrhs = SIZE(beta_gaussian, 2)
      CPASSERT(SIZE(beta_gaussian, 1) == n)
      CPASSERT(ALL(SHAPE(coefficients) == SHAPE(beta_gaussian)))
      coefficients = beta_gaussian
      IF (n == 0 .OR. nrhs == 0) RETURN
      CALL zgetrs('N', n, nrhs, inverse%lu, n, inverse%pivots, coefficients, n, info)
      CPASSERT(info == 0)
   END SUBROUTINE upf_inverse_project

! **************************************************************************************************
!> \brief Differentiate the same reconstruction: K*dW=dB-dK*W.
!>        Displacement/strain responses must include the original waves, projectors,
!>        their moving integration measure, and all relevant Bloch-image terms.
!> \param inverse ...
!> \param coefficients ...
!> \param beta_gaussian_response ...
!> \param beta_delta_response ...
!> \param response ...
! **************************************************************************************************
   SUBROUTINE upf_inverse_project_response(inverse, coefficients, beta_gaussian_response, &
                                           beta_delta_response, response)
      TYPE(upf_inverse_type), INTENT(IN)                 :: inverse
      COMPLEX(KIND=dp), INTENT(IN)                       :: coefficients(:, :), &
                                                            beta_gaussian_response(:, :), &
                                                            beta_delta_response(:, :)
      COMPLEX(KIND=dp), INTENT(OUT)                      :: response(:, :)

      COMPLEX(KIND=dp), ALLOCATABLE                      :: rhs(:, :)
      INTEGER                                            :: n

      n = inverse%nproj
      CPASSERT(SIZE(coefficients, 1) == n)
      CPASSERT(ALL(SHAPE(beta_delta_response) == [n, n]))
      CPASSERT(ALL(SHAPE(beta_gaussian_response) == SHAPE(coefficients)))
      ALLOCATE (rhs(n, SIZE(coefficients, 2)))
      rhs(:, :) = beta_gaussian_response - MATMUL(beta_delta_response, coefficients)
      CALL upf_inverse_project(inverse, rhs, response)
   END SUBROUTINE upf_inverse_project_response

! **************************************************************************************************
!> \brief Propagate a real energy covector through K*W=B using the existing full LU factors.
!>        The complex pairing is Re sum(conjg(G)*dW). Its reverse solve is K**H*GB=G,
!>        with GK=-GB*W**H. No symmetry of K, occupation factor or regularization is assumed.
!> \param inverse Original factored reconstruction matrix K.
!> \param coefficients Original solution W.
!> \param coefficients_covector Real energy derivative with respect to the complex W.
!> \param beta_gaussian_covector Covector with respect to B.
!> \param beta_delta_covector Covector with respect to K=I+P*Delta.
! **************************************************************************************************
   SUBROUTINE upf_inverse_project_adjoint(inverse, coefficients, coefficients_covector, &
                                          beta_gaussian_covector, beta_delta_covector)
      TYPE(upf_inverse_type), INTENT(IN)                 :: inverse
      COMPLEX(KIND=dp), INTENT(IN)                       :: coefficients(:, :), &
                                                            coefficients_covector(:, :)
      COMPLEX(KIND=dp), INTENT(OUT)                      :: beta_gaussian_covector(:, :), &
                                                            beta_delta_covector(:, :)

      INTEGER                                            :: info, n, nrhs

      EXTERNAL zgetrs

      CPASSERT(inverse%ready)
      n = inverse%nproj
      nrhs = SIZE(coefficients, 2)
      CPASSERT(SIZE(coefficients, 1) == n)
      CPASSERT(ALL(SHAPE(coefficients_covector) == SHAPE(coefficients)))
      CPASSERT(ALL(SHAPE(beta_gaussian_covector) == SHAPE(coefficients)))
      CPASSERT(ALL(SHAPE(beta_delta_covector) == [n, n]))
      beta_gaussian_covector = coefficients_covector
      beta_delta_covector = CMPLX(0.0_dp, 0.0_dp, dp)
      IF (n == 0 .OR. nrhs == 0) RETURN
      CALL zgetrs('C', n, nrhs, inverse%lu, n, inverse%pivots, beta_gaussian_covector, n, info)
      CPASSERT(info == 0)
      beta_delta_covector = -MATMUL(beta_gaussian_covector, CONJG(TRANSPOSE(coefficients)))
   END SUBROUTINE upf_inverse_project_adjoint

! **************************************************************************************************
!> \brief Evaluate inverse images on a field batch: U*chi=chi-Delta*W.
!>        Spatial derivatives use the corresponding derivative fields, with W fixed.
!> \param gaussian ...
!> \param delta ...
!> \param coefficients ...
!> \param soft ...
! **************************************************************************************************
   SUBROUTINE upf_inverse_basis(gaussian, delta, coefficients, soft)
      COMPLEX(KIND=dp), INTENT(IN)                       :: gaussian(:, :), delta(:, :), &
                                                            coefficients(:, :)
      COMPLEX(KIND=dp), INTENT(OUT)                      :: soft(:, :)

      CPASSERT(SIZE(delta, 1) == SIZE(gaussian, 1))
      CPASSERT(SIZE(delta, 2) == SIZE(coefficients, 1))
      CPASSERT(SIZE(coefficients, 2) == SIZE(gaussian, 2))
      CPASSERT(ALL(SHAPE(soft) == SHAPE(gaussian)))
      soft = gaussian - MATMUL(delta, coefficients)
   END SUBROUTINE upf_inverse_basis

! **************************************************************************************************
!> \brief Joint field response d(U*chi)=dchi-dDelta*W-Delta*dW.
!> \param gaussian_response ...
!> \param delta ...
!> \param delta_response ...
!> \param coefficients ...
!> \param coefficients_response ...
!> \param response ...
! **************************************************************************************************
   SUBROUTINE upf_inverse_basis_response(gaussian_response, delta, delta_response, &
                                         coefficients, coefficients_response, response)
      COMPLEX(KIND=dp), INTENT(IN) :: gaussian_response(:, :), delta(:, :), delta_response(:, :), &
         coefficients(:, :), coefficients_response(:, :)
      COMPLEX(KIND=dp), INTENT(OUT)                      :: response(:, :)

      CPASSERT(ALL(SHAPE(delta_response) == SHAPE(delta)))
      CPASSERT(ALL(SHAPE(coefficients_response) == SHAPE(coefficients)))
      CALL upf_inverse_basis(gaussian_response, delta_response, coefficients, response)
      response = response - MATMUL(delta, coefficients_response)
   END SUBROUTINE upf_inverse_basis_response

! **************************************************************************************************
!> \brief Pull back A+P^dagger*augmentation*P into reconstructed Gaussian variables.
!> \param gaussian_operator <chi|A|chi>.
!> \param gaussian_delta <chi|A|Delta>.
!> \param delta_gaussian <Delta|A|chi>; supplied separately for non-Hermitian operators.
!> \param delta_operator <Delta|A|Delta>.
!> \param augmentation Original projector-space contribution; q for the metric,
!>                     Dij+integral(v*Q) for nonlocal/local-augmentation terms.
!> \param coefficients W=K^-1*B. Original operators are retained, including residuals.
!> \param transformed <U*chi|A|U*chi>+W^dagger*augmentation*W.
! **************************************************************************************************
   SUBROUTINE upf_inverse_operator(gaussian_operator, gaussian_delta, delta_gaussian, delta_operator, &
                                   augmentation, coefficients, transformed)
      COMPLEX(KIND=dp), INTENT(IN) :: gaussian_operator(:, :), gaussian_delta(:, :), &
         delta_gaussian(:, :), delta_operator(:, :), augmentation(:, :), coefficients(:, :)
      COMPLEX(KIND=dp), INTENT(OUT)                      :: transformed(:, :)

      CALL check_operator_shapes(gaussian_operator, gaussian_delta, delta_gaussian, delta_operator, &
                                 augmentation, coefficients, SHAPE(transformed))
      transformed = gaussian_operator - MATMUL(gaussian_delta, coefficients) - &
                    MATMUL(CONJG(TRANSPOSE(coefficients)), delta_gaussian) + &
                    MATMUL(CONJG(TRANSPOSE(coefficients)), MATMUL(delta_operator + augmentation, coefficients))
   END SUBROUTINE upf_inverse_operator

! **************************************************************************************************
!> \brief Differentiate the full operator pullback, including both inverse-basis legs.
!>        Matrix responses already include their moving original fields and quadratures.
!>        The stationary energy also requires the transformed-metric Pulay term.
!> \param gaussian_response ...
!> \param gaussian_delta ...
!> \param delta_gaussian ...
!> \param delta_operator ...
!> \param augmentation ...
!> \param coefficients ...
!> \param gaussian_delta_response ...
!> \param delta_gaussian_response ...
!> \param delta_operator_response ...
!> \param augmentation_response ...
!> \param coefficients_response ...
!> \param response ...
! **************************************************************************************************
   SUBROUTINE upf_inverse_operator_response(gaussian_response, gaussian_delta, delta_gaussian, &
                                            delta_operator, augmentation, coefficients, &
                                            gaussian_delta_response, delta_gaussian_response, &
                                            delta_operator_response, augmentation_response, &
                                            coefficients_response, response)
      COMPLEX(KIND=dp), INTENT(IN) :: gaussian_response(:, :), gaussian_delta(:, :), &
         delta_gaussian(:, :), delta_operator(:, :), augmentation(:, :), coefficients(:, :), &
         gaussian_delta_response(:, :), delta_gaussian_response(:, :), &
         delta_operator_response(:, :), augmentation_response(:, :), coefficients_response(:, :)
      COMPLEX(KIND=dp), INTENT(OUT)                      :: response(:, :)

      CALL check_operator_shapes(gaussian_response, gaussian_delta, delta_gaussian, delta_operator, &
                                 augmentation, coefficients, SHAPE(response))
      CPASSERT(ALL(SHAPE(coefficients_response) == SHAPE(coefficients)))
      CALL upf_inverse_operator(gaussian_response, gaussian_delta_response, delta_gaussian_response, &
                                delta_operator_response, augmentation_response, coefficients, response)
      response = response - MATMUL(gaussian_delta, coefficients_response) - &
                 MATMUL(CONJG(TRANSPOSE(coefficients_response)), delta_gaussian) + &
                 MATMUL(CONJG(TRANSPOSE(coefficients_response)), MATMUL(delta_operator + augmentation, coefficients)) + &
                 MATMUL(CONJG(TRANSPOSE(coefficients)), MATMUL(delta_operator + augmentation, coefficients_response))
   END SUBROUTINE upf_inverse_operator_response

! **************************************************************************************************
!> \brief ...
!> \param gg ...
!> \param gd ...
!> \param dg ...
!> \param dd ...
!> \param augmentation ...
!> \param coefficients ...
!> \param output_shape ...
! **************************************************************************************************
   SUBROUTINE check_operator_shapes(gg, gd, dg, dd, augmentation, coefficients, output_shape)
      COMPLEX(KIND=dp), INTENT(IN)                       :: gg(:, :), gd(:, :), dg(:, :), dd(:, :), &
                                                            augmentation(:, :), coefficients(:, :)
      INTEGER, INTENT(IN)                                :: output_shape(2)

      INTEGER                                            :: m, n

      m = SIZE(coefficients, 1)
      n = SIZE(coefficients, 2)
      CPASSERT(ALL(SHAPE(gg) == [n, n]))
      CPASSERT(ALL(SHAPE(gd) == [n, m]))
      CPASSERT(ALL(SHAPE(dg) == [m, n]))
      CPASSERT(ALL(SHAPE(dd) == [m, m]))
      CPASSERT(ALL(SHAPE(augmentation) == [m, m]))
      CPASSERT(ALL(output_shape == [n, n]))
   END SUBROUTINE check_operator_shapes

END MODULE upf_inverse
