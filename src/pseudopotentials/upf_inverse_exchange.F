!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Complex inverse transition integrals and their paired forward/reverse derivatives.
!>        Integrate s_left* s_right plus the original projector augmentation. Independent left
!>        and right fields allow different Bloch vectors and rectangular AO blocks. No Hermitian
!>        projection, symmetry packing or occupation factor is applied here.
! **************************************************************************************************
MODULE upf_inverse_exchange
   USE gaussian_coulomb_fields,         ONLY: gaussian_coulomb_fields_evaluate,&
                                              gaussian_coulomb_fields_type
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: fourpi
   USE message_passing,                 ONLY: mp_comm_type
   USE spherical_harmonics,             ONLY: y_lm
   USE uspp_augmentation,               ONLY: uspp_augmentation_type
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: upf_inverse_exchange_block, upf_inverse_exchange_response, upf_inverse_exchange_adjoint
   PUBLIC :: upf_inverse_exchange_contract
   PUBLIC :: upf_inverse_exchange_augmentation
   PUBLIC :: upf_inverse_exchange_kernel, upf_inverse_exchange_kernel_response
   PUBLIC :: upf_inverse_exchange_kernel_adjoint

CONTAINS

! **************************************************************************************************
!> \brief Effective RI interaction from the owner's fitting inverse and exchange interaction.
!>        The supplied inverse retains the owner's conditioning, local auxiliary selection and
!>        metric weighting. No independent inverse or regularization is introduced here.
!> \param fitting_inverse Hermitian fitting inverse A, prepared by the native RI owner.
!> \param exchange Hermitian exchange interaction V_x in the same auxiliary coordinates.
!> \param interaction New effective kernel A V_x A.
! **************************************************************************************************
   SUBROUTINE upf_inverse_exchange_kernel(fitting_inverse, exchange, interaction)
      COMPLEX(dp), INTENT(IN)                            :: fitting_inverse(:, :), exchange(:, :)
      COMPLEX(dp), INTENT(OUT)                           :: interaction(:, :)

      CALL check_hermitian(fitting_inverse)
      CALL check_hermitian(exchange)
      CPASSERT(ALL(SHAPE(exchange) == SHAPE(fitting_inverse)))
      CPASSERT(ALL(SHAPE(interaction) == SHAPE(fitting_inverse)))
      interaction(:, :) = MATMUL(MATMUL(fitting_inverse, exchange), fitting_inverse)
   END SUBROUTINE upf_inverse_exchange_kernel

! **************************************************************************************************
!> \brief Response of the effective RI kernel, including the owner's fitting inverse response.
!>        For a fixed auxiliary span dA=-A dV_f A. Local metric weighting, changing auxiliary
!>        support, or conditioning can supply additional terms in dA; retain the full owner map.
!> \param fitting_inverse ...
!> \param exchange ...
!> \param fitting_inverse_response Full native RI inverse response, including weighting.
!> \param exchange_response Full native exchange interaction response, including its kernel.
!> \param response New effective-kernel response.
! **************************************************************************************************
   SUBROUTINE upf_inverse_exchange_kernel_response(fitting_inverse, exchange, fitting_inverse_response, &
                                                   exchange_response, response)
      COMPLEX(dp), INTENT(IN)                            :: fitting_inverse(:, :), exchange(:, :), &
                                                            fitting_inverse_response(:, :), &
                                                            exchange_response(:, :)
      COMPLEX(dp), INTENT(OUT)                           :: response(:, :)

      COMPLEX(dp), ALLOCATABLE                           :: work(:, :)

      CALL check_hermitian(fitting_inverse)
      CALL check_hermitian(exchange)
      CALL check_hermitian(fitting_inverse_response)
      CALL check_hermitian(exchange_response)
      CPASSERT(ALL(SHAPE(exchange) == SHAPE(fitting_inverse)))
      CPASSERT(ALL(SHAPE(fitting_inverse_response) == SHAPE(fitting_inverse)))
      CPASSERT(ALL(SHAPE(exchange_response) == SHAPE(fitting_inverse)))
      CPASSERT(ALL(SHAPE(response) == SHAPE(fitting_inverse)))
      ALLOCATE (work(SIZE(response, 1), SIZE(response, 2)))
      work(:, :) = MATMUL(fitting_inverse, exchange)
      response(:, :) = MATMUL(work, fitting_inverse_response)
      work(:, :) = MATMUL(fitting_inverse_response, exchange)
      response = response + MATMUL(work, fitting_inverse)
      work(:, :) = MATMUL(fitting_inverse, exchange_response)
      response = response + MATMUL(work, fitting_inverse)
   END SUBROUTINE upf_inverse_exchange_kernel_response

! **************************************************************************************************
!> \brief Paired RI inverse and exchange-interaction adjoints for a Hermitian effective kernel.
!>        Adjoint pairing is Re sum(conjg(barX)*dX). The native RI owner transports the full
!>        inverse adjoint through its own weighted inverse, rather than assuming a global metric.
!> \param fitting_inverse ...
!> \param exchange ...
!> \param interaction_adjoint Complete Hermitian effective-kernel energy adjoint.
!> \param fitting_inverse_adjoint New fitting-inverse covector.
!> \param exchange_adjoint New exchange-interaction covector.
! **************************************************************************************************
   SUBROUTINE upf_inverse_exchange_kernel_adjoint(fitting_inverse, exchange, interaction_adjoint, &
                                                  fitting_inverse_adjoint, exchange_adjoint)
      COMPLEX(dp), INTENT(IN)                            :: fitting_inverse(:, :), exchange(:, :), &
                                                            interaction_adjoint(:, :)
      COMPLEX(dp), INTENT(OUT)                           :: fitting_inverse_adjoint(:, :), &
                                                            exchange_adjoint(:, :)

      COMPLEX(dp), ALLOCATABLE                           :: work(:, :)

      CALL check_hermitian(fitting_inverse)
      CALL check_hermitian(exchange)
      CALL check_hermitian(interaction_adjoint)
      CPASSERT(ALL(SHAPE(exchange) == SHAPE(fitting_inverse)))
      CPASSERT(ALL(SHAPE(interaction_adjoint) == SHAPE(fitting_inverse)))
      CPASSERT(ALL(SHAPE(fitting_inverse_adjoint) == SHAPE(fitting_inverse)))
      CPASSERT(ALL(SHAPE(exchange_adjoint) == SHAPE(fitting_inverse)))
      ALLOCATE (work(SIZE(exchange, 1), SIZE(exchange, 2)))
      work(:, :) = MATMUL(fitting_inverse, exchange)
      fitting_inverse_adjoint(:, :) = MATMUL(interaction_adjoint, work)
      work(:, :) = MATMUL(exchange, fitting_inverse)
      fitting_inverse_adjoint = fitting_inverse_adjoint + MATMUL(work, interaction_adjoint)
      work(:, :) = MATMUL(fitting_inverse, interaction_adjoint)
      exchange_adjoint(:, :) = MATMUL(work, fitting_inverse)
   END SUBROUTINE upf_inverse_exchange_kernel_adjoint

! **************************************************************************************************
!> \brief Contract one ordered, physical-spin RI exchange pair and all its energy adjoints.
!>        E = prefactor Re sum_AB interaction(A,B) Tr(P_left I_B P_right I_A^dagger).
!>        The Hermitian interaction is the effective fitted exchange kernel. For Coulomb fitting
!>        with a different exchange interaction it is V_fit^-1 V_exchange V_fit^-1, rather than
!>        just the inverse fitting metric. Transitions contain the auxiliary conjugation.
!>        A complex interaction therefore multiplies I_B, with I_A conjugated; transposing it
!>        changes the functional even though a real Gamma test would not detect that error.
!> \param transitions Complete rectangular (left AO, right AO, auxiliary) transitions.
!> \param interaction Complete Hermitian effective auxiliary interaction.
!> \param left_density Hermitian physical occupation matrix for one spin at the left k point.
!> \param right_density Independent right occupation matrix in the same spin.
!> \param prefactor Includes the exchange fraction, minus one-half, and both k weights once.
!> \param energy New contribution of this ordered spin/k pair.
!> \param left_operator New derivative conjugate to P_left, including the supplied prefactor.
!> \param right_operator New derivative conjugate to P_right, including the supplied prefactor.
!>        When both densities represent the same variable, add both derivatives. A native
!>        Hamiltonian also removes the target k weight from the accumulated energy derivative.
!> \param transition_adjoint Optional new complex transition adjoint for the real energy pairing.
!> \param interaction_adjoint Optional new effective-kernel adjoint, including metric geometry.
! **************************************************************************************************
   SUBROUTINE upf_inverse_exchange_contract(transitions, interaction, left_density, right_density, &
                                            prefactor, energy, left_operator, right_operator, &
                                            transition_adjoint, interaction_adjoint)
      COMPLEX(dp), INTENT(IN)                            :: transitions(:, :, :), interaction(:, :), &
                                                            left_density(:, :), right_density(:, :)
      REAL(dp), INTENT(IN)                               :: prefactor
      REAL(dp), INTENT(OUT)                              :: energy
      COMPLEX(dp), INTENT(OUT)                           :: left_operator(:, :), right_operator(:, :)
      COMPLEX(dp), INTENT(OUT), OPTIONAL                 :: transition_adjoint(:, :, :), &
                                                            interaction_adjoint(:, :)

      COMPLEX(dp), ALLOCATABLE                           :: fitted(:, :), left_work(:, :), &
                                                            right_work(:, :), work(:, :)
      INTEGER                                            :: a, b, naux, nl, nr

      nl = SIZE(transitions, 1)
      nr = SIZE(transitions, 2)
      naux = SIZE(transitions, 3)
      CPASSERT(ALL(SHAPE(left_density) == [nl, nl]))
      CPASSERT(ALL(SHAPE(right_density) == [nr, nr]))
      CPASSERT(ALL(SHAPE(interaction) == [naux, naux]))
      CPASSERT(ALL(SHAPE(left_operator) == [nl, nl]))
      CPASSERT(ALL(SHAPE(right_operator) == [nr, nr]))
      CALL check_hermitian(left_density)
      CALL check_hermitian(right_density)
      CALL check_hermitian(interaction)
      IF (PRESENT(transition_adjoint)) THEN
         CPASSERT(ALL(SHAPE(transition_adjoint) == SHAPE(transitions)))
      END IF
      IF (PRESENT(interaction_adjoint)) THEN
         CPASSERT(ALL(SHAPE(interaction_adjoint) == SHAPE(interaction)))
      END IF
      ALLOCATE (fitted(nl, nr), work(nl, nr), left_work(nl, nl), right_work(nr, nr))
      energy = 0.0_dp
      left_operator = CMPLX(0.0_dp, 0.0_dp, dp)
      right_operator = CMPLX(0.0_dp, 0.0_dp, dp)
      DO a = 1, naux
         fitted(:, :) = CMPLX(0.0_dp, 0.0_dp, dp)
         DO b = 1, naux
            fitted(:, :) = fitted + interaction(a, b)*transitions(:, :, b)
         END DO
         work(:, :) = MATMUL(fitted, right_density)
         left_work(:, :) = MATMUL(work, CONJG(TRANSPOSE(transitions(:, :, a))))
         left_operator = left_operator + prefactor*left_work
         energy = energy + prefactor*REAL(SUM(TRANSPOSE(left_density)*left_work), dp)
         work(:, :) = MATMUL(left_density, fitted)
         right_work(:, :) = MATMUL(CONJG(TRANSPOSE(transitions(:, :, a))), work)
         right_operator = right_operator + prefactor*right_work
         IF (PRESENT(transition_adjoint)) THEN
            transition_adjoint(:, :, a) = 2.0_dp*prefactor*MATMUL(work, right_density)
         END IF
         IF (PRESENT(interaction_adjoint)) THEN
            DO b = 1, naux
               work(:, :) = MATMUL(MATMUL(left_density, transitions(:, :, b)), right_density)
               interaction_adjoint(a, b) = prefactor*CONJG(SUM(work*CONJG(transitions(:, :, a))))
            END DO
         END IF
      END DO
      ! These variables are Hermitian. Remove only their contraction roundoff, without
      ! projecting the generally rectangular, non-Hermitian transitions themselves.
      left_operator = 0.5_dp*(left_operator + CONJG(TRANSPOSE(left_operator)))
      right_operator = 0.5_dp*(right_operator + CONJG(TRANSPOSE(right_operator)))
   END SUBROUTINE upf_inverse_exchange_contract

! **************************************************************************************************
!> \brief Original Q integrals against a native Gaussian Coulomb-potential cache.
!>        Original radial weights, angular kernels and complete nonzero tails are retained. The
!>        integral is in the native magnetic-projector and Gaussian-function coordinates, with
!>        no partial-wave product replacement or compensation subtraction. Image translations
!>        and Bloch phases belong to the caller.
!> \param fields Prepared native Gaussian auxiliary shell-set cache.
!> \param augmentation Original scalar projector augmentation with its radial quadrature.
!> \param displacement Augmentation center minus Gaussian center in bohr.
!> \param directions Original unit vectors of the angular quadrature.
!> \param weights Original angular weights, summing to 4*pi.
!> \param integrals New complete (projector, projector, auxiliary) block.
!> \param gaussian_derivative Optional new Gaussian-center derivative, Cartesian axis last.
!> \param augmentation_derivative Optional new augmentation-center derivative, Cartesian axis last.
! **************************************************************************************************
   SUBROUTINE upf_inverse_exchange_augmentation(fields, augmentation, displacement, directions, weights, &
                                                integrals, gaussian_derivative, augmentation_derivative)
      TYPE(gaussian_coulomb_fields_type), INTENT(IN)     :: fields
      TYPE(uspp_augmentation_type), INTENT(IN)           :: augmentation
      REAL(dp), INTENT(IN)                               :: displacement(3), directions(:, :), &
                                                            weights(:)
      REAL(dp), INTENT(OUT)                              :: integrals(:, :, :)
      REAL(dp), INTENT(OUT), OPTIONAL :: gaussian_derivative(:, :, :, :), &
         augmentation_derivative(:, :, :, :)

      INTEGER                                            :: axis, i, ia, ir, j, l, lm, m, n, naux, &
                                                            p, pair, q
      LOGICAL                                            :: derivatives
      REAL(dp)                                           :: point(3), radial
      REAL(dp), ALLOCATABLE :: d_integrals(:, :, :, :), d_projected(:, :, :), gradient(:, :), &
         harmonics(:, :), potential(:), projected(:, :), values(:)

      n = augmentation%nproj
      naux = fields%nfunction
      CPASSERT(ALL(SHAPE(integrals) == [n, n, naux]))
      integrals = 0.0_dp
      derivatives = PRESENT(gaussian_derivative) .OR. PRESENT(augmentation_derivative)
      IF (PRESENT(gaussian_derivative)) THEN
         CPASSERT(ALL(SHAPE(gaussian_derivative) == [n, n, naux, 3]))
         gaussian_derivative = 0.0_dp
      END IF
      IF (PRESENT(augmentation_derivative)) THEN
         CPASSERT(ALL(SHAPE(augmentation_derivative) == [n, n, naux, 3]))
         augmentation_derivative = 0.0_dp
      END IF
      IF (n == 0) RETURN
      CPASSERT(SIZE(directions, 1) == 3 .AND. SIZE(directions, 2) == SIZE(weights))
      CPASSERT(SIZE(weights) > 0)
      CPASSERT(ABS(SUM(weights) - fourpi) <= 1.0e-10_dp)
      CPASSERT(MAXVAL(ABS(SUM(directions**2, DIM=1) - 1.0_dp)) <= 1.0e-12_dp)
      ALLOCATE (harmonics(SIZE(weights), (augmentation%lmax + 1)**2), values(SIZE(weights)), &
                potential(naux), projected((augmentation%lmax + 1)**2, naux))
      IF (derivatives) THEN
         ALLOCATE (gradient(3, naux), d_projected((augmentation%lmax + 1)**2, naux, 3), &
                   d_integrals(n, n, naux, 3), SOURCE=0.0_dp)
      END IF
      DO l = 0, augmentation%lmax
         DO m = -l, l
            lm = l*l + l + m + 1
            CALL y_lm(directions, values, l, m)
            harmonics(:, lm) = weights*values
         END DO
      END DO
      DO ir = 1, SIZE(augmentation%weighted_q, 1)
         IF (ALL(augmentation%weighted_q(ir, :, :) == 0.0_dp)) CYCLE
         projected(:, :) = 0.0_dp
         IF (derivatives) d_projected(:, :, :) = 0.0_dp
         DO ia = 1, SIZE(weights)
            point = displacement + augmentation%r(ir)*directions(:, ia)
            IF (derivatives) THEN
               CALL gaussian_coulomb_fields_evaluate(fields, point, potential, gradient)
               DO axis = 1, 3
                  d_projected(:, :, axis) = d_projected(:, :, axis) + &
                                            SPREAD(harmonics(ia, :), 2, naux)*SPREAD(gradient(axis, :), 1, SIZE(harmonics, 2))
               END DO
            ELSE
               CALL gaussian_coulomb_fields_evaluate(fields, point, potential)
            END IF
            projected(:, :) = projected + SPREAD(harmonics(ia, :), 2, naux)*SPREAD(potential, 1, SIZE(harmonics, 2))
         END DO
         DO q = 1, n
            j = augmentation%radial_index(q)
            DO p = 1, q
               i = augmentation%radial_index(p)
               pair = MAX(i, j)*(MAX(i, j) - 1)/2 + MIN(i, j)
               DO l = ABS(augmentation%lbeta(i) - augmentation%lbeta(j)), augmentation%lbeta(i) + augmentation%lbeta(j), 2
                  radial = augmentation%weighted_q(ir, pair, l)
                  IF (radial == 0.0_dp) CYCLE
                  DO lm = l*l + 1, (l + 1)**2
                     integrals(p, q, :) = integrals(p, q, :) + radial*augmentation%gaunt(p, q, lm)*projected(lm, :)
                     IF (derivatives) THEN
                        d_integrals(p, q, :, :) = d_integrals(p, q, :, :) + &
                                                  radial*augmentation%gaunt(p, q, lm)*d_projected(lm, :, :)
                     END IF
                  END DO
               END DO
            END DO
         END DO
      END DO
      DO q = 1, n
         DO p = 1, q - 1
            integrals(q, p, :) = integrals(p, q, :)
            IF (derivatives) d_integrals(q, p, :, :) = d_integrals(p, q, :, :)
         END DO
      END DO
      IF (PRESENT(gaussian_derivative)) gaussian_derivative = -d_integrals
      IF (PRESENT(augmentation_derivative)) augmentation_derivative = d_integrals
   END SUBROUTINE upf_inverse_exchange_augmentation

! **************************************************************************************************
!> \brief I(A) = left^dagger diag(weights*potential(A)) right + W_left^dagger Q(A) W_right.
!>        Q(A) is the complete original augmentation contracted with the same auxiliary potential.
!>        In particular, Q is not replaced by an AE-minus-PS product. W is taken from the defining
!>        inverse frame, rather than recomputed on the auxiliary integration points.
!> \param left Local point rows of the left soft fields.
!> \param right Local point rows of the independent right soft fields.
!> \param weights Local original integration measure.
!> \param potential Local complex auxiliary potentials; includes the auxiliary conjugation/phase.
!> \param left_coefficients Original left W, replicated within the point communicator.
!> \param right_coefficients Original right W, replicated within the point communicator.
!> \param augmentation Complete integrated original Q, replicated within the point communicator.
!> \param integrals New complete rectangular (left AO, right AO, auxiliary) block.
!> \param comm Optional communicator partitioning point rows; empty ranks participate.
! **************************************************************************************************
   SUBROUTINE upf_inverse_exchange_block(left, right, weights, potential, left_coefficients, &
                                         right_coefficients, augmentation, integrals, comm)
      COMPLEX(dp), INTENT(IN)                            :: left(:, :), right(:, :)
      REAL(dp), INTENT(IN)                               :: weights(:)
      COMPLEX(dp), INTENT(IN)                            :: potential(:, :), &
                                                            left_coefficients(:, :), &
                                                            right_coefficients(:, :), &
                                                            augmentation(:, :, :)
      COMPLEX(dp), INTENT(OUT)                           :: integrals(:, :, :)
      TYPE(mp_comm_type), INTENT(IN), OPTIONAL           :: comm

      INTEGER                                            :: auxiliary

      CALL check_shapes(left, right, weights, potential, left_coefficients, right_coefficients, &
                        augmentation, SHAPE(integrals))
      DO auxiliary = 1, SIZE(potential, 2)
         CALL weighted_pair(left, right, weights*potential(:, auxiliary), integrals(:, :, auxiliary))
      END DO
      ! Only the point integral is distributed. Reducing the replicated Q contribution here
      ! would multiply it by the number of ranks, including ranks with no point rows.
      IF (PRESENT(comm)) CALL comm%sum(integrals)
      DO auxiliary = 1, SIZE(potential, 2)
         integrals(:, :, auxiliary) = integrals(:, :, auxiliary) + &
                                      MATMUL(CONJG(TRANSPOSE(left_coefficients)), &
                                             MATMUL(augmentation(:, :, auxiliary), right_coefficients))
      END DO
   END SUBROUTINE upf_inverse_exchange_block

! **************************************************************************************************
!> \brief Joint first response of the complete transition integral.
!>        Field responses include the coupled inverse response and all moving original fields.
!>        The auxiliary potential and integration measure also move. Q_response includes the
!>        original augmentation's center, auxiliary and periodic-image responses.
!> \param left ...
!> \param right ...
!> \param weights ...
!> \param potential ...
!> \param left_coefficients ...
!> \param right_coefficients ...
!> \param augmentation ...
!> \param left_response Joint response of the left soft fields.
!> \param right_response Joint response of the right soft fields.
!> \param weights_response Response of the local integration measure.
!> \param potential_response Response of the local complex auxiliary potentials.
!> \param left_coefficients_response dW_left from the original defining inverse frame.
!> \param right_coefficients_response dW_right from the original defining inverse frame.
!> \param augmentation_response Complete integrated original Q response, replicated as Q.
!> \param response New complete rectangular response block.
!> \param comm ...
! **************************************************************************************************
   SUBROUTINE upf_inverse_exchange_response(left, right, weights, potential, left_coefficients, &
                                            right_coefficients, augmentation, left_response, right_response, &
                                            weights_response, potential_response, left_coefficients_response, &
                                            right_coefficients_response, augmentation_response, response, comm)
      COMPLEX(dp), INTENT(IN)                            :: left(:, :), right(:, :)
      REAL(dp), INTENT(IN)                               :: weights(:)
      COMPLEX(dp), INTENT(IN) :: potential(:, :), left_coefficients(:, :), &
         right_coefficients(:, :), augmentation(:, :, :), left_response(:, :), right_response(:, :)
      REAL(dp), INTENT(IN)                               :: weights_response(:)
      COMPLEX(dp), INTENT(IN) :: potential_response(:, :), left_coefficients_response(:, :), &
         right_coefficients_response(:, :), augmentation_response(:, :, :)
      COMPLEX(dp), INTENT(OUT)                           :: response(:, :, :)
      TYPE(mp_comm_type), INTENT(IN), OPTIONAL           :: comm

      COMPLEX(dp), ALLOCATABLE                           :: work(:, :)
      INTEGER                                            :: auxiliary

      CALL check_shapes(left, right, weights, potential, left_coefficients, right_coefficients, &
                        augmentation, SHAPE(response))
      CPASSERT(ALL(SHAPE(left_response) == SHAPE(left)))
      CPASSERT(ALL(SHAPE(right_response) == SHAPE(right)))
      CPASSERT(SIZE(weights_response) == SIZE(weights))
      CPASSERT(ALL(SHAPE(potential_response) == SHAPE(potential)))
      CPASSERT(ALL(SHAPE(left_coefficients_response) == SHAPE(left_coefficients)))
      CPASSERT(ALL(SHAPE(right_coefficients_response) == SHAPE(right_coefficients)))
      CPASSERT(ALL(SHAPE(augmentation_response) == SHAPE(augmentation)))
      ALLOCATE (work(SIZE(left, 2), SIZE(right, 2)))
      DO auxiliary = 1, SIZE(potential, 2)
         CALL weighted_pair(left_response, right, weights*potential(:, auxiliary), response(:, :, auxiliary))
         CALL weighted_pair(left, right_response, weights*potential(:, auxiliary), work)
         response(:, :, auxiliary) = response(:, :, auxiliary) + work
         CALL weighted_pair(left, right, weights_response*potential(:, auxiliary) + &
                            weights*potential_response(:, auxiliary), work)
         response(:, :, auxiliary) = response(:, :, auxiliary) + work
      END DO
      IF (PRESENT(comm)) CALL comm%sum(response)
      DO auxiliary = 1, SIZE(potential, 2)
         response(:, :, auxiliary) = response(:, :, auxiliary) + &
                                     MATMUL(CONJG(TRANSPOSE(left_coefficients_response)), &
                                            MATMUL(augmentation(:, :, auxiliary), right_coefficients)) + &
                                     MATMUL(CONJG(TRANSPOSE(left_coefficients)), &
                                            MATMUL(augmentation_response(:, :, auxiliary), right_coefficients)) + &
                                     MATMUL(CONJG(TRANSPOSE(left_coefficients)), &
                                            MATMUL(augmentation(:, :, auxiliary), right_coefficients_response))
      END DO
   END SUBROUTINE upf_inverse_exchange_response

! **************************************************************************************************
!> \brief Reverse derivative for dE = Re sum(conjg(integral_adjoint)*dI), with no extra factor.
!>        Point adjoints remain local. The replicated W and Q adjoints are returned once and must
!>        not be summed over the point communicator. The caller reduces local point/measure
!>        contractions and adds replicated contractions once. Integral adjoints already include
!>        the exchange prefactor, occupations, k weights and any RI metric contractions.
!> \param left ...
!> \param right ...
!> \param weights ...
!> \param potential ...
!> \param left_coefficients ...
!> \param right_coefficients ...
!> \param augmentation ...
!> \param integral_adjoint Complete replicated adjoint of the rectangular transition integrals.
!> \param left_adjoint New local left field adjoint.
!> \param right_adjoint New local right field adjoint.
!> \param weight_adjoint New local real measure adjoint.
!> \param potential_adjoint New local complex auxiliary potential adjoint.
!> \param left_coefficients_adjoint New replicated original left W adjoint.
!> \param right_coefficients_adjoint New replicated original right W adjoint.
!> \param augmentation_adjoint New replicated complete original Q adjoint.
! **************************************************************************************************
   SUBROUTINE upf_inverse_exchange_adjoint(left, right, weights, potential, left_coefficients, &
                                           right_coefficients, augmentation, integral_adjoint, &
                                           left_adjoint, right_adjoint, weight_adjoint, potential_adjoint, &
                                           left_coefficients_adjoint, right_coefficients_adjoint, augmentation_adjoint)
      COMPLEX(dp), INTENT(IN)                            :: left(:, :), right(:, :)
      REAL(dp), INTENT(IN)                               :: weights(:)
      COMPLEX(dp), INTENT(IN) :: potential(:, :), left_coefficients(:, :), &
         right_coefficients(:, :), augmentation(:, :, :), integral_adjoint(:, :, :)
      COMPLEX(dp), INTENT(OUT)                           :: left_adjoint(:, :), right_adjoint(:, :)
      REAL(dp), INTENT(OUT)                              :: weight_adjoint(:)
      COMPLEX(dp), INTENT(OUT) :: potential_adjoint(:, :), left_coefficients_adjoint(:, :), &
         right_coefficients_adjoint(:, :), augmentation_adjoint(:, :, :)

      COMPLEX(dp), ALLOCATABLE                           :: kernel_adjoint(:), work(:, :)
      INTEGER                                            :: auxiliary

      CALL check_shapes(left, right, weights, potential, left_coefficients, right_coefficients, &
                        augmentation, SHAPE(integral_adjoint))
      CPASSERT(ALL(SHAPE(left_adjoint) == SHAPE(left)))
      CPASSERT(ALL(SHAPE(right_adjoint) == SHAPE(right)))
      CPASSERT(SIZE(weight_adjoint) == SIZE(weights))
      CPASSERT(ALL(SHAPE(potential_adjoint) == SHAPE(potential)))
      CPASSERT(ALL(SHAPE(left_coefficients_adjoint) == SHAPE(left_coefficients)))
      CPASSERT(ALL(SHAPE(right_coefficients_adjoint) == SHAPE(right_coefficients)))
      CPASSERT(ALL(SHAPE(augmentation_adjoint) == SHAPE(augmentation)))
      left_adjoint = CMPLX(0.0_dp, 0.0_dp, dp)
      right_adjoint = CMPLX(0.0_dp, 0.0_dp, dp)
      weight_adjoint = 0.0_dp
      left_coefficients_adjoint = CMPLX(0.0_dp, 0.0_dp, dp)
      right_coefficients_adjoint = CMPLX(0.0_dp, 0.0_dp, dp)
      ALLOCATE (work(SIZE(left, 1), SIZE(right, 2)), kernel_adjoint(SIZE(weights)))
      DO auxiliary = 1, SIZE(potential, 2)
         left_adjoint = left_adjoint + &
                        SPREAD(weights*potential(:, auxiliary), 2, SIZE(left, 2))* &
                        MATMUL(right, CONJG(TRANSPOSE(integral_adjoint(:, :, auxiliary))))
         work(:, :) = MATMUL(left, integral_adjoint(:, :, auxiliary))
         right_adjoint = right_adjoint + &
                         SPREAD(weights*CONJG(potential(:, auxiliary)), 2, SIZE(right, 2))*work
         kernel_adjoint(:) = SUM(work*CONJG(right), DIM=2)
         potential_adjoint(:, auxiliary) = weights*kernel_adjoint
         weight_adjoint = weight_adjoint + REAL(CONJG(kernel_adjoint)*potential(:, auxiliary), dp)
         left_coefficients_adjoint = left_coefficients_adjoint + &
                                     MATMUL(MATMUL(augmentation(:, :, auxiliary), right_coefficients), &
                                            CONJG(TRANSPOSE(integral_adjoint(:, :, auxiliary))))
         right_coefficients_adjoint = right_coefficients_adjoint + &
                                      MATMUL(CONJG(TRANSPOSE(augmentation(:, :, auxiliary))), &
                                             MATMUL(left_coefficients, integral_adjoint(:, :, auxiliary)))
         augmentation_adjoint(:, :, auxiliary) = &
            MATMUL(MATMUL(left_coefficients, integral_adjoint(:, :, auxiliary)), &
                   CONJG(TRANSPOSE(right_coefficients)))
      END DO
   END SUBROUTINE upf_inverse_exchange_adjoint

! **************************************************************************************************
!> \brief Integrate a complex multiplication kernel on local point rows.
!> \param left ...
!> \param right ...
!> \param kernel ...
!> \param integral ...
! **************************************************************************************************
   SUBROUTINE weighted_pair(left, right, kernel, integral)
      COMPLEX(dp), INTENT(IN)                            :: left(:, :), right(:, :), kernel(:)
      COMPLEX(dp), INTENT(OUT)                           :: integral(:, :)

      integral = MATMUL(CONJG(TRANSPOSE(left)), SPREAD(kernel, 2, SIZE(right, 2))*right)
   END SUBROUTINE weighted_pair

! **************************************************************************************************
!> \brief Check dimensions without imposing a density or symmetry convention on a transition.
!> \param left ...
!> \param right ...
!> \param weights ...
!> \param potential ...
!> \param left_coefficients ...
!> \param right_coefficients ...
!> \param augmentation ...
!> \param output_shape ...
! **************************************************************************************************
   SUBROUTINE check_shapes(left, right, weights, potential, left_coefficients, right_coefficients, &
                           augmentation, output_shape)
      COMPLEX(dp), INTENT(IN)                            :: left(:, :), right(:, :)
      REAL(dp), INTENT(IN)                               :: weights(:)
      COMPLEX(dp), INTENT(IN)                            :: potential(:, :), &
                                                            left_coefficients(:, :), &
                                                            right_coefficients(:, :), &
                                                            augmentation(:, :, :)
      INTEGER, INTENT(IN)                                :: output_shape(3)

      LOGICAL                                            :: valid

      CPASSERT(SIZE(left, 1) == SIZE(weights) .AND. SIZE(right, 1) == SIZE(weights))
      CPASSERT(SIZE(potential, 1) == SIZE(weights))
      CPASSERT(SIZE(left_coefficients, 2) == SIZE(left, 2))
      CPASSERT(SIZE(right_coefficients, 2) == SIZE(right, 2))
      valid = ALL(SHAPE(augmentation) == &
                  [SIZE(left_coefficients, 1), SIZE(right_coefficients, 1), SIZE(potential, 2)])
      CPASSERT(valid)
      CPASSERT(ALL(output_shape == [SIZE(left, 2), SIZE(right, 2), SIZE(potential, 2)]))
   END SUBROUTINE check_shapes

! **************************************************************************************************
!> \brief Check the Hermitian physical density or effective interaction, including empty blocks.
!> \param matrix ...
! **************************************************************************************************
   SUBROUTINE check_hermitian(matrix)
      COMPLEX(dp), INTENT(IN)                            :: matrix(:, :)

      REAL(dp)                                           :: error, scale

      CPASSERT(SIZE(matrix, 1) == SIZE(matrix, 2))
      IF (SIZE(matrix) == 0) RETURN
      scale = MAX(1.0_dp, MAXVAL(ABS(matrix)))
      error = MAXVAL(ABS(matrix - CONJG(TRANSPOSE(matrix))))
      CPASSERT(error <= 1.0e-12_dp*scale)
   END SUBROUTINE check_hermitian

END MODULE upf_inverse_exchange
