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
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_comm_type
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: upf_inverse_exchange_block, upf_inverse_exchange_response, upf_inverse_exchange_adjoint

CONTAINS

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

END MODULE upf_inverse_exchange
