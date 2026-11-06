!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Contains methods used in the context of density fitting
!> \par History
!>      04.2008 created [Manuel Guidon]
!>      02.2013 moved from admm_methods
!> \author Manuel Guidon
! **************************************************************************************************
MODULE admm_utils
   USE admm_types,                      ONLY: admm_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_add,&
                                              dbcsr_copy,&
                                              dbcsr_create,&
                                              dbcsr_deallocate_matrix,&
                                              dbcsr_set,&
                                              dbcsr_type,&
                                              dbcsr_type_symmetric
   USE cp_dbcsr_operations,             ONLY: copy_fm_to_dbcsr
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE input_constants,                 ONLY: do_admm_purify_cauchy,&
                                              do_admm_purify_cauchy_subspace,&
                                              do_admm_purify_mo_diag,&
                                              do_admm_purify_mo_no_diag,&
                                              do_admm_purify_none
   USE kinds,                           ONLY: dp
   USE parallel_gemm_api,               ONLY: parallel_gemm
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: admm_inverse_sqrt_divdiff, admm_purified_operator, admm_charge_scaled_operator
   PUBLIC :: admm_charge_scaled_operator_response
   PUBLIC :: admm_correct_for_eigenvalues, &
             admm_uncorrect_for_eigenvalues

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'admm_utils'

!***

CONTAINS

! **************************************************************************************************
!> \brief Charge-ratio chain rule for F(R,eta), with R=D or R=eta*D in one spin sector.
!>        The actual native auxiliary density R is supplied. With the raw charge N_A=Tr(M_A D),
!>        lambda=[s'(eta) Tr(K_R D)+partial_eta F]/N_A and K_D=s K_R-eta lambda M_A.
!>        Traces include physical k weights once. Add lambda M_O in the primary operator and
!>        -eta lambda D / lambda P in the auxiliary / primary geometry metric covectors.
!> \param auxiliary_operator Unweighted Hermitian derivative K_R at fixed eta, in auxiliary physical-k order.
!> \param density Actual native auxiliary density R in the same coordinates.
!> \param auxiliary_metric Physical valence charge metric M_A, excluding frozen cores.
!> \param primary_density Primary physical density P.
!> \param primary_metric Primary valence charge metric M_O.
!> \param weights Physical normalized mesh weights.
!> \param eta Positive charge ratio N_O/N_A, shared over the full mesh.
!> \param normalized_density True for R=eta D. False for R=D with explicit exchange rescaling.
!> \param explicit_derivative Partial derivative of the complete energy with respect to eta at fixed R.
!> \param corrected Raw-density derivative K_D including the auxiliary charge normalization term.
!> \param lambda Primary charge-normalization multiplier.
! **************************************************************************************************
   SUBROUTINE admm_charge_scaled_operator(auxiliary_operator, density, auxiliary_metric, primary_density, primary_metric, &
                                          weights, eta, normalized_density, explicit_derivative, corrected, lambda)
      COMPLEX(dp), INTENT(IN) :: auxiliary_operator(:, :, :), density(:, :, :), &
         auxiliary_metric(:, :, :), primary_density(:, :, :), primary_metric(:, :, :)
      REAL(dp), INTENT(IN)                               :: weights(:), eta
      LOGICAL, INTENT(IN)                                :: normalized_density
      REAL(dp), INTENT(IN)                               :: explicit_derivative
      COMPLEX(dp), INTENT(OUT)                           :: corrected(:, :, :)
      REAL(dp), INTENT(OUT)                              :: lambda

      INTEGER                                            :: ik, na, nk, no
      LOGICAL                                            :: valid
      REAL(dp)                                           :: auxiliary_charge, numerator, &
                                                            primary_charge, scale, tolerance

      na = SIZE(density, 1)
      no = SIZE(primary_density, 1)
      nk = SIZE(weights)
      valid = ALL(SHAPE(auxiliary_operator) == [na, na, nk]) .AND. ALL(SHAPE(density) == [na, na, nk])
      CPASSERT(valid)
      valid = ALL(SHAPE(auxiliary_metric) == [na, na, nk]) .AND. ALL(SHAPE(corrected) == [na, na, nk])
      CPASSERT(valid)
      valid = ALL(SHAPE(primary_density) == [no, no, nk]) .AND. ALL(SHAPE(primary_metric) == [no, no, nk])
      CPASSERT(valid)
      valid = ieee_is_finite(eta) .AND. eta > 0.0_dp .AND. ieee_is_finite(explicit_derivative)
      CPASSERT(valid)
      CPASSERT(ALL(ieee_is_finite(weights)) .AND. ALL(weights >= 0.0_dp))
      scale = 1.0_dp
      IF (normalized_density) scale = eta
      auxiliary_charge = 0.0_dp
      primary_charge = 0.0_dp
      numerator = explicit_derivative
      DO ik = 1, nk
         auxiliary_charge = auxiliary_charge + weights(ik)* &
                            REAL(SUM(CONJG(auxiliary_metric(:, :, ik))*density(:, :, ik)), dp)/scale
         primary_charge = primary_charge + weights(ik)* &
                          REAL(SUM(CONJG(primary_metric(:, :, ik))*primary_density(:, :, ik)), dp)
         IF (normalized_density) numerator = numerator + weights(ik)* &
                                             REAL(SUM(CONJG(auxiliary_operator(:, :, ik))*density(:, :, ik)), dp)/scale
      END DO
      valid = ieee_is_finite(auxiliary_charge) .AND. ieee_is_finite(primary_charge) .AND. ieee_is_finite(numerator)
      CPASSERT(valid)
      lambda = 0.0_dp
      tolerance = 128.0_dp*EPSILON(1.0_dp)*MAX(1.0_dp, ABS(auxiliary_charge), ABS(primary_charge))
      IF (ABS(auxiliary_charge) > tolerance .OR. ABS(primary_charge) > tolerance) THEN
         CPASSERT(auxiliary_charge > 0.0_dp)
         lambda = numerator/auxiliary_charge
      END IF
      corrected(:, :, :) = scale*auxiliary_operator - eta*lambda*auxiliary_metric
   END SUBROUTINE admm_charge_scaled_operator

! **************************************************************************************************
!> \brief Divided difference of x**(-1/2) at positive metric eigenvalues, including equality.
!>        Rationalization avoids subtraction of nearly equal roots and a truncated Taylor series.
!> \param a First positive eigenvalue.
!> \param b Second positive eigenvalue.
!> \return Exact inverse-square-root derivative kernel, up to floating-point rounding.
! **************************************************************************************************
   PURE FUNCTION admm_inverse_sqrt_divdiff(a, b) RESULT(value)
      REAL(dp), INTENT(IN)                               :: a, b
      REAL(dp)                                           :: value

      REAL(dp)                                           :: root_a, root_b

      root_a = SQRT(a)
      root_b = SQRT(b)
      value = -1.0_dp/(root_a + root_b)/root_a/root_b
   END FUNCTION admm_inverse_sqrt_divdiff

! **************************************************************************************************
!> \brief Derivative of the MO_DIAG auxiliary density on the fixed-rank primary density manifold.
!>        Primary orbitals are orthonormal in the physical primary metric. Their positive
!>        occupations define the primary density. The auxiliary projection is orthonormalized
!>        with its actual Gram inverse square root. The returned Hermitian AO derivative also
!>        gives the occupation derivatives. Geometry includes the primary-metric response of
!>        the natural orbitals, in addition to the auxiliary and mixed metric responses.
!> \param projection Physical projection A = M_AA**(-1) M_AO.
!> \param primary_metric Physical primary normalization metric M_OO.
!> \param auxiliary_metric Physical auxiliary normalization metric M_AA.
!> \param coefficients Primary orbitals C with C**dagger M_OO C = I.
!> \param occupations Positive retained orbital occupations.
!> \param gram_vectors Eigenvectors of (A C)**dagger M_AA (A C).
!> \param gram_values Positive Gram eigenvalues in the same order.
!> \param inverse_root The same Gram inverse square root used to build the auxiliary density.
!> \param auxiliary_operator Hermitian derivative at fixed purified auxiliary density.
!> \param OPERATOR Hermitian derivative with respect to the primary AO density.
!> \param metric_inverse Optional inverse M_AA for geometry. Supply all four geometry arguments.
!> \param mixed_covector Optional derivative with respect to M_AO, with the factor two included.
!> \param auxiliary_covector Optional Hermitian derivative with respect to M_AA.
!> \param primary_covector Optional Hermitian derivative with respect to M_OO at fixed AO density.
! **************************************************************************************************
   SUBROUTINE admm_purified_operator(projection, primary_metric, auxiliary_metric, coefficients, occupations, &
                                     gram_vectors, gram_values, inverse_root, auxiliary_operator, OPERATOR, &
                                     metric_inverse, mixed_covector, auxiliary_covector, primary_covector)
      COMPLEX(dp), INTENT(IN)                            :: projection(:, :), primary_metric(:, :), &
                                                            auxiliary_metric(:, :), &
                                                            coefficients(:, :)
      REAL(dp), INTENT(IN)                               :: occupations(:)
      COMPLEX(dp), INTENT(IN)                            :: gram_vectors(:, :)
      REAL(dp), INTENT(IN)                               :: gram_values(:)
      COMPLEX(dp), INTENT(IN)                            :: inverse_root(:, :), &
                                                            auxiliary_operator(:, :)
      COMPLEX(dp), INTENT(OUT)                           :: operator(:, :)
      COMPLEX(dp), INTENT(IN), OPTIONAL                  :: metric_inverse(:, :)
      COMPLEX(dp), INTENT(OUT), OPTIONAL                 :: mixed_covector(:, :), &
                                                            auxiliary_covector(:, :), &
                                                            primary_covector(:, :)

      COMPLEX(dp), ALLOCATABLE :: delta(:, :), eigen_covector(:, :), gram_covector(:, :), &
         metric_orbitals(:, :), occupied_gradient(:, :), occupied_operator(:, :), &
         orbital_gradient(:, :), projected(:, :), purified(:, :), tangent(:, :), &
         weighted_gradient(:, :), z(:, :)
      INTEGER                                            :: i, j, m, na, no
      LOGICAL                                            :: geometry, valid

      no = SIZE(coefficients, 1)
      m = SIZE(coefficients, 2)
      na = SIZE(projection, 1)
      valid = ALL(SHAPE(projection) == [na, no]) .AND. &
              ALL(SHAPE(primary_metric) == [no, no]) .AND. &
              ALL(SHAPE(auxiliary_metric) == [na, na]) .AND. &
              ALL(SHAPE(auxiliary_operator) == [na, na]) .AND. &
              ALL(SHAPE(OPERATOR) == [no, no]) .AND. &
              ALL(SHAPE(gram_vectors) == [m, m]) .AND. &
              ALL(SHAPE(inverse_root) == [m, m]) .AND. &
              SIZE(occupations) == m .AND. SIZE(gram_values) == m
      CPASSERT(valid)
      geometry = PRESENT(metric_inverse) .OR. PRESENT(mixed_covector) .OR. &
                 PRESENT(auxiliary_covector) .OR. PRESENT(primary_covector)
      IF (geometry) THEN
         valid = PRESENT(metric_inverse) .AND. PRESENT(mixed_covector) .AND. &
                 PRESENT(auxiliary_covector) .AND. PRESENT(primary_covector)
         CPASSERT(valid)
         valid = ALL(SHAPE(metric_inverse) == [na, na]) .AND. &
                 ALL(SHAPE(mixed_covector) == [na, no]) .AND. &
                 ALL(SHAPE(auxiliary_covector) == [na, na]) .AND. &
                 ALL(SHAPE(primary_covector) == [no, no])
         CPASSERT(valid)
         mixed_covector = CMPLX(0.0_dp, 0.0_dp, dp)
         auxiliary_covector = CMPLX(0.0_dp, 0.0_dp, dp)
         primary_covector = CMPLX(0.0_dp, 0.0_dp, dp)
      END IF
      OPERATOR = CMPLX(0.0_dp, 0.0_dp, dp)
      IF (m == 0) RETURN
      CPASSERT(ALL(occupations > 0.0_dp) .AND. ALL(gram_values > 0.0_dp))
      ALLOCATE (projected, SOURCE=MATMUL(projection, coefficients))
      ALLOCATE (purified, SOURCE=MATMUL(projected, inverse_root))
      ALLOCATE (weighted_gradient, SOURCE=MATMUL(auxiliary_operator, purified))
      DO i = 1, m
         weighted_gradient(:, i) = weighted_gradient(:, i)*occupations(i)
      END DO
      ALLOCATE (gram_covector, SOURCE=MATMUL(CONJG(TRANSPOSE(projected)), weighted_gradient))
      gram_covector(:, :) = gram_covector + CONJG(TRANSPOSE(gram_covector))
      ALLOCATE (eigen_covector, SOURCE=MATMUL(CONJG(TRANSPOSE(gram_vectors)), MATMUL(gram_covector, gram_vectors)))
      DO j = 1, m
         DO i = 1, m
            eigen_covector(i, j) = eigen_covector(i, j)*admm_inverse_sqrt_divdiff(gram_values(i), gram_values(j))
         END DO
      END DO
      gram_covector(:, :) = MATMUL(gram_vectors, MATMUL(eigen_covector, CONJG(TRANSPOSE(gram_vectors))))
      weighted_gradient(:, :) = MATMUL(weighted_gradient, inverse_root) + &
                                MATMUL(auxiliary_metric, MATMUL(projected, gram_covector))
      ALLOCATE (orbital_gradient, SOURCE=MATMUL(CONJG(TRANSPOSE(projection)), weighted_gradient))
      ALLOCATE (occupied_gradient, SOURCE=MATMUL(CONJG(TRANSPOSE(coefficients)), orbital_gradient))
      ALLOCATE (occupied_operator, SOURCE=MATMUL(CONJG(TRANSPOSE(purified)), MATMUL(auxiliary_operator, purified)))
      ALLOCATE (metric_orbitals, SOURCE=MATMUL(primary_metric, coefficients))
      ALLOCATE (tangent, SOURCE=orbital_gradient - MATMUL(metric_orbitals, occupied_gradient))
      DO i = 1, m
         tangent(:, i) = tangent(:, i)/occupations(i)
      END DO
      OPERATOR(:, :) = MATMUL(metric_orbitals, MATMUL(occupied_operator, CONJG(TRANSPOSE(metric_orbitals)))) + &
                       MATMUL(tangent, CONJG(TRANSPOSE(metric_orbitals))) + &
                       MATMUL(metric_orbitals, CONJG(TRANSPOSE(tangent)))
      OPERATOR(:, :) = 0.5_dp*(OPERATOR + CONJG(TRANSPOSE(OPERATOR)))
      IF (geometry) THEN
         ALLOCATE (z, SOURCE=MATMUL(metric_inverse, MATMUL(weighted_gradient, CONJG(TRANSPOSE(coefficients)))))
         mixed_covector(:, :) = 2.0_dp*z
         auxiliary_covector(:, :) = MATMUL(projected, MATMUL(gram_covector, CONJG(TRANSPOSE(projected)))) - &
                                    MATMUL(z, CONJG(TRANSPOSE(projection))) - MATMUL(projection, CONJG(TRANSPOSE(z)))
         ALLOCATE (delta, SOURCE=occupied_operator)
         DO i = 1, m
            delta(:, i) = delta(:, i)*occupations(i)
         END DO
         delta(:, :) = delta - occupied_gradient
         delta(:, :) = 0.5_dp*(delta + CONJG(TRANSPOSE(delta)))
         primary_covector(:, :) = MATMUL(coefficients, MATMUL(delta, CONJG(TRANSPOSE(coefficients))))
      END IF
   END SUBROUTINE admm_purified_operator

! **************************************************************************************************
!> \brief ...
!> \param ispin ...
!> \param admm_env ...
!> \param ks_matrix ...
! **************************************************************************************************
   SUBROUTINE admm_correct_for_eigenvalues(ispin, admm_env, ks_matrix)
      INTEGER, INTENT(IN)                                :: ispin
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(dbcsr_type), POINTER                          :: ks_matrix

      INTEGER                                            :: nao_aux_fit, nao_orb
      TYPE(dbcsr_type), POINTER                          :: work

      nao_aux_fit = admm_env%nao_aux_fit
      nao_orb = admm_env%nao_orb

      IF (.NOT. admm_env%block_dm) THEN
         SELECT CASE (admm_env%purification_method)
         CASE (do_admm_purify_cauchy_subspace)
            !* remove what has been added and add the correction
            NULLIFY (work)
            ALLOCATE (work)
            CALL dbcsr_create(work, template=ks_matrix, name='work', matrix_type=dbcsr_type_symmetric)

            CALL dbcsr_copy(work, ks_matrix)
            CALL dbcsr_set(work, 0.0_dp)
            CALL copy_fm_to_dbcsr(admm_env%ks_to_be_merged(ispin), work, keep_sparsity=.TRUE.)

            CALL dbcsr_add(ks_matrix, work, 1.0_dp, -1.0_dp)

            ! ** calculate A^T*H_tilde*A
            CALL parallel_gemm('N', 'N', nao_aux_fit, nao_orb, nao_aux_fit, &
                               1.0_dp, admm_env%K(ispin), admm_env%A, 0.0_dp, &
                               admm_env%work_aux_orb)
            CALL parallel_gemm('T', 'N', nao_orb, nao_orb, nao_aux_fit, &
                               1.0_dp, admm_env%A, admm_env%work_aux_orb, 0.0_dp, &
                               admm_env%H_corr(ispin))

            CALL copy_fm_to_dbcsr(admm_env%H_corr(ispin), work, keep_sparsity=.TRUE.)

            CALL dbcsr_add(ks_matrix, work, 1.0_dp, 1.0_dp)
            CALL dbcsr_deallocate_matrix(work)

         CASE (do_admm_purify_mo_diag)
            !* remove what has been added and add the correction
            NULLIFY (work)
            ALLOCATE (work)
            CALL dbcsr_create(work, template=ks_matrix, name='work', matrix_type=dbcsr_type_symmetric)

            CALL dbcsr_copy(work, ks_matrix)
            CALL dbcsr_set(work, 0.0_dp)
            CALL copy_fm_to_dbcsr(admm_env%ks_to_be_merged(ispin), work, keep_sparsity=.TRUE.)

            ! ** calculate A^T*H_tilde*A
            CALL parallel_gemm('N', 'N', nao_aux_fit, nao_orb, nao_aux_fit, &
                               1.0_dp, admm_env%K(ispin), admm_env%A, 0.0_dp, &
                               admm_env%work_aux_orb)
            CALL parallel_gemm('T', 'N', nao_orb, nao_orb, nao_aux_fit, &
                               1.0_dp, admm_env%A, admm_env%work_aux_orb, 0.0_dp, &
                               admm_env%H_corr(ispin))

            CALL copy_fm_to_dbcsr(admm_env%H_corr(ispin), work, keep_sparsity=.TRUE.)

            CALL dbcsr_add(ks_matrix, work, 1.0_dp, 1.0_dp)
            CALL dbcsr_deallocate_matrix(work)

         CASE (do_admm_purify_mo_no_diag, do_admm_purify_none, do_admm_purify_cauchy)
            ! do nothing
         END SELECT
      END IF

   END SUBROUTINE admm_correct_for_eigenvalues

! **************************************************************************************************
!> \brief ...
!> \param ispin ...
!> \param admm_env ...
!> \param ks_matrix ...
! **************************************************************************************************
   SUBROUTINE admm_uncorrect_for_eigenvalues(ispin, admm_env, ks_matrix)
      INTEGER, INTENT(IN)                                :: ispin
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(dbcsr_type), POINTER                          :: ks_matrix

      INTEGER                                            :: nao_aux_fit, nao_orb
      TYPE(dbcsr_type), POINTER                          :: work

      nao_aux_fit = admm_env%nao_aux_fit
      nao_orb = admm_env%nao_orb

      IF (.NOT. admm_env%block_dm) THEN
         SELECT CASE (admm_env%purification_method)
         CASE (do_admm_purify_cauchy_subspace)
            !* remove what has been added and add the correction
            NULLIFY (work)
            ALLOCATE (work)
            CALL dbcsr_create(work, template=ks_matrix, name='work', matrix_type=dbcsr_type_symmetric)

            CALL dbcsr_copy(work, ks_matrix)
            CALL dbcsr_set(work, 0.0_dp)
            CALL copy_fm_to_dbcsr(admm_env%H_corr(ispin), work, keep_sparsity=.TRUE.)

            CALL dbcsr_add(ks_matrix, work, 1.0_dp, -1.0_dp)

            CALL copy_fm_to_dbcsr(admm_env%H_corr(ispin), work, keep_sparsity=.TRUE.)

            CALL dbcsr_set(work, 0.0_dp)
            CALL copy_fm_to_dbcsr(admm_env%ks_to_be_merged(ispin), work, keep_sparsity=.TRUE.)

            CALL dbcsr_add(ks_matrix, work, 1.0_dp, 1.0_dp)
            CALL dbcsr_deallocate_matrix(work)

         CASE (do_admm_purify_mo_diag)
            NULLIFY (work)
            ALLOCATE (work)
            CALL dbcsr_create(work, template=ks_matrix, name='work', matrix_type=dbcsr_type_symmetric)

            CALL dbcsr_copy(work, ks_matrix)
            CALL dbcsr_set(work, 0.0_dp)

            CALL copy_fm_to_dbcsr(admm_env%H_corr(ispin), work, keep_sparsity=.TRUE.)

            CALL dbcsr_add(ks_matrix, work, 1.0_dp, -1.0_dp)
            CALL dbcsr_deallocate_matrix(work)

         CASE (do_admm_purify_mo_no_diag, do_admm_purify_none, do_admm_purify_cauchy)
            ! do nothing
         END SELECT
      END IF
   END SUBROUTINE admm_uncorrect_for_eigenvalues

! **************************************************************************************************
!> \brief Fixed-geometry electronic response of the ADMM charge-ratio operator.
!>        Differentiate the actual native-density charge and numerator used above.
!>        The projection, both charge metrics and physical mesh weights stay fixed.
!>        The caller supplies the complete auxiliary operator response, including
!>        its explicit eta dependence, and the corresponding partial_eta response.
!> \param auxiliary_operator Derivative of the auxiliary energy at fixed eta.
!> \param operator_direction Total electronic variation of that derivative.
!> \param density Actual auxiliary density R, either D or eta D.
!> \param density_direction Variation of that same actual density.
!> \param auxiliary_metric Original auxiliary valence charge metric.
!> \param primary_density Original primary density for the charge-sector selection.
!> \param primary_metric Original primary valence charge metric.
!> \param weights Normalized physical k-point weights, applied once.
!> \param eta Original positive primary-to-raw-auxiliary charge ratio.
!> \param eta_direction Variation of that same charge ratio.
!> \param normalized_density True for R=eta D, false for R=D.
!> \param explicit_derivative Partial_eta of the auxiliary energy at fixed R.
!> \param explicit_direction Total variation of that partial derivative.
!> \param corrected_direction Raw-auxiliary operator response before projection.
!> \param lambda_direction Primary charge-multiplier response.
! **************************************************************************************************
   SUBROUTINE admm_charge_scaled_operator_response(auxiliary_operator, operator_direction, density, density_direction, &
                                                   auxiliary_metric, primary_density, primary_metric, weights, &
                                                   eta, eta_direction, normalized_density, explicit_derivative, &
                                                   explicit_direction, corrected_direction, lambda_direction)
      COMPLEX(dp), INTENT(IN) :: auxiliary_operator(:, :, :), operator_direction(:, :, :), &
         density(:, :, :), density_direction(:, :, :), auxiliary_metric(:, :, :), &
         primary_density(:, :, :), primary_metric(:, :, :)
      REAL(dp), INTENT(IN)                               :: weights(:), eta, eta_direction
      LOGICAL, INTENT(IN)                                :: normalized_density
      REAL(dp), INTENT(IN)                               :: explicit_derivative, explicit_direction
      COMPLEX(dp), INTENT(OUT)                           :: corrected_direction(:, :, :)
      REAL(dp), INTENT(OUT)                              :: lambda_direction

      COMPLEX(dp), ALLOCATABLE                           :: corrected(:, :, :)
      INTEGER                                            :: ik
      LOGICAL                                            :: valid
      REAL(dp) :: auxiliary_charge, charge_direction, lambda, numerator_direction, primary_charge, &
         scale, scale_direction, tolerance, trace_direction, trace_value

      valid = ALL(SHAPE(operator_direction) == SHAPE(auxiliary_operator)) .AND. &
              ALL(SHAPE(density_direction) == SHAPE(density)) .AND. &
              ALL(SHAPE(corrected_direction) == SHAPE(auxiliary_operator))
      CPASSERT(valid)
      CPASSERT(ieee_is_finite(eta_direction) .AND. ieee_is_finite(explicit_direction))
      ALLOCATE (corrected, MOLD=auxiliary_operator)
      CALL admm_charge_scaled_operator(auxiliary_operator, density, auxiliary_metric, primary_density, primary_metric, &
                                       weights, eta, normalized_density, explicit_derivative, corrected, lambda)
      scale = 1.0_dp
      scale_direction = 0.0_dp
      IF (normalized_density) THEN
         scale = eta
         scale_direction = eta_direction
      END IF
      auxiliary_charge = 0.0_dp
      charge_direction = 0.0_dp
      primary_charge = 0.0_dp
      numerator_direction = explicit_direction
      DO ik = 1, SIZE(weights)
         trace_value = REAL(SUM(CONJG(auxiliary_metric(:, :, ik))*density(:, :, ik)), dp)
         trace_direction = REAL(SUM(CONJG(auxiliary_metric(:, :, ik))*density_direction(:, :, ik)), dp)
         auxiliary_charge = auxiliary_charge + weights(ik)*trace_value/scale
         charge_direction = charge_direction + weights(ik)* &
                            (trace_direction/scale - trace_value*scale_direction/scale**2)
         primary_charge = primary_charge + weights(ik)* &
                          REAL(SUM(CONJG(primary_metric(:, :, ik))*primary_density(:, :, ik)), dp)
         IF (normalized_density) THEN
            trace_value = REAL(SUM(CONJG(auxiliary_operator(:, :, ik))*density(:, :, ik)), dp)
            trace_direction = REAL(SUM(CONJG(operator_direction(:, :, ik))*density(:, :, ik)) + &
                                   SUM(CONJG(auxiliary_operator(:, :, ik))*density_direction(:, :, ik)), dp)
            numerator_direction = numerator_direction + weights(ik)* &
                                  (trace_direction/scale - trace_value*scale_direction/scale**2)
         END IF
      END DO
      tolerance = 128.0_dp*EPSILON(1.0_dp)*MAX(1.0_dp, ABS(auxiliary_charge), ABS(primary_charge))
      lambda_direction = 0.0_dp
      IF (ABS(auxiliary_charge) > tolerance .OR. ABS(primary_charge) > tolerance) THEN
         CPASSERT(auxiliary_charge > 0.0_dp)
         lambda_direction = (numerator_direction - lambda*charge_direction)/auxiliary_charge
      END IF
      corrected_direction(:, :, :) = scale*operator_direction + scale_direction*auxiliary_operator - &
                                     (eta_direction*lambda + eta*lambda_direction)*auxiliary_metric
   END SUBROUTINE admm_charge_scaled_operator_response

END MODULE admm_utils
