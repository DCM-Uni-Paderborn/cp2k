!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
!> \brief Home-center inverse density and its finite Gaussian GAPW angular map.
MODULE qs_upf_inverse_density
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE gaussian_basis_fields,           ONLY: gaussian_basis_fields_periodic
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_comm_type
   USE paw_basis_types,                 ONLY: get_paw_basis_info
   USE qs_grid_atom,                    ONLY: grid_atom_type
   USE qs_harmonics_atom,               ONLY: harmonics_atom_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_upf_inverse,                  ONLY: qs_upf_inverse_center_fields,&
                                              qs_upf_inverse_evaluate
   USE qs_upf_inverse_gaussian,         ONLY: qs_upf_inverse_gaussian_integrate,&
                                              qs_upf_inverse_gaussian_radial
   USE qs_upf_inverse_gaussian_pairs,   ONLY: qs_upf_inverse_gaussian_pairs_density,&
                                              qs_upf_inverse_gaussian_pairs_integrate,&
                                              qs_upf_inverse_gaussian_pairs_type
   USE qs_upf_inverse_types,            ONLY: qs_upf_inverse_augmentation_type,&
                                              qs_upf_inverse_type
   USE upf_inverse,                     ONLY: upf_inverse_kinetic_evaluate,&
                                              upf_inverse_kinetic_type
   USE upf_inverse_quadrature,          ONLY: upf_inverse_cross_density,&
                                              upf_inverse_projector_density,&
                                              upf_inverse_semilocal_pair
   USE uspp_augmentation,               ONLY: uspp_augmentation_collocate_atom,&
                                              uspp_augmentation_integrate_atom
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: qs_upf_inverse_density_center, qs_upf_inverse_density_integrate
   PUBLIC :: qs_upf_inverse_density_radial, qs_upf_inverse_density_radial_adjoint
   PUBLIC :: qs_upf_inverse_density_gapw, qs_upf_inverse_density_gapw_integrate
   PUBLIC :: qs_upf_inverse_kinetic_center, qs_upf_inverse_kinetic_integrate
   PUBLIC :: qs_upf_inverse_kinetic_collocate, qs_upf_inverse_kinetic_operator

CONTAINS

! **************************************************************************************************
!> \brief Periodic valence kinetic augmentation and its complete fixed-matrix geometry response.
!>        Same-center kernels have no Bloch phase. Coupled inverse coefficients retain all k data.
!> \param frame Original defining inverse with the common finite image list.
!> \param kinetic Original kind-indexed AE/PS kernels. NC kinds have zero channels.
!> \param points Cartesian evaluation points, without partition weights.
!> \param matrix Hermitian AO occupations for one spin and k point, without mesh weights.
!> \param density Periodic kinetic augmentation, overwritten with factor one half.
!> \param response Optional fixed-matrix response, requiring all response inputs below.
!> \param point_response Cartesian motion of the evaluation points.
!> \param atom_response Cartesian motion of the original atoms.
!> \param cell_response Lattice response, including each translated image center.
!> \param coefficients_response Complete dW from the defining inverse.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_kinetic_collocate(frame, kinetic, points, matrix, density, response, &
                                               point_response, atom_response, cell_response, coefficients_response)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      TYPE(upf_inverse_kinetic_type), INTENT(IN)         :: kinetic(:)
      REAL(dp), INTENT(IN)                               :: points(:, :)
      COMPLEX(dp), INTENT(IN)                            :: matrix(:, :)
      REAL(dp), INTENT(OUT)                              :: density(:)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: response(:)
      REAL(dp), INTENT(IN), OPTIONAL                     :: point_response(:, :), &
                                                            atom_response(:, :), &
                                                            cell_response(3, 3)
      COMPLEX(dp), INTENT(IN), OPTIONAL                  :: coefficients_response(:, :)

      COMPLEX(dp), ALLOCATABLE                           :: projected(:, :), projected_response(:, :)
      INTEGER                                            :: atom, axis, center, first, ikind, image, &
                                                            last, n, nao, np, point
      LOGICAL                                            :: derive
      REAL(dp)                                           :: position(3), position_response(3)
      REAL(dp), ALLOCATABLE                              :: gradient(:, :, :), kernel(:, :)

      nao = SIZE(frame%coefficients, 2)
      np = SIZE(points, 2)
      CPASSERT(frame%atom_local .OR. frame%inverse%ready)
      CPASSERT(SIZE(points, 1) == 3)
      CPASSERT(ALL(SHAPE(matrix) == [nao, nao]) .AND. SIZE(density) == np)
      derive = PRESENT(response)
      CPASSERT(derive .EQV. PRESENT(point_response))
      CPASSERT(derive .EQV. PRESENT(atom_response))
      CPASSERT(derive .EQV. PRESENT(cell_response))
      CPASSERT(derive .EQV. PRESENT(coefficients_response))
      density = 0.0_dp
      IF (derive) THEN
         CPASSERT(SIZE(response) == np .AND. ALL(SHAPE(point_response) == SHAPE(points)))
         CPASSERT(ALL(SHAPE(atom_response) == SHAPE(frame%positions)))
         CPASSERT(ALL(SHAPE(coefficients_response) == SHAPE(frame%coefficients)))
         response = 0.0_dp
      END IF
      DO center = 1, SIZE(frame%center_atoms)
         ikind = frame%center_kinds(center)
         CPASSERT(ikind >= 1 .AND. ikind <= SIZE(kinetic))
         IF (kinetic(ikind)%delta%nproj == 0) CYCLE
         CALL kinetic_center_indices(frame, kinetic(ikind), center, points, atom, first, last)
         n = last - first + 1
         ALLOCATE (projected(n, n), kernel(n, n))
         projected(:, :) = MATMUL(frame%coefficients(first:last, :), &
                                  MATMUL(matrix, CONJG(TRANSPOSE(frame%coefficients(first:last, :)))))
         IF (derive) THEN
            ALLOCATE (projected_response(n, n), gradient(n, n, 3))
            projected_response(:, :) = MATMUL(coefficients_response(first:last, :), &
                                              MATMUL(matrix, CONJG(TRANSPOSE(frame%coefficients(first:last, :)))))
            projected_response(:, :) = projected_response + CONJG(TRANSPOSE(projected_response))
         END IF
         DO image = 1, SIZE(frame%images, 2)
            position = frame%positions(:, atom) + MATMUL(frame%cell, REAL(frame%images(:, image), dp))
            IF (derive) position_response = atom_response(:, atom) + &
                                            MATMUL(cell_response, REAL(frame%images(:, image), dp))
            DO point = 1, np
               IF (SUM((points(:, point) - position)**2) >= kinetic(ikind)%delta%radius**2) CYCLE
               IF (derive) THEN
                  CALL upf_inverse_kinetic_evaluate(kinetic(ikind), points(:, point) - position, kernel, gradient)
                  response(point) = response(point) + REAL(SUM(projected_response*kernel), dp)
                  DO axis = 1, 3
                     response(point) = response(point) + REAL(SUM(projected*gradient(:, :, axis)), dp)* &
                                       (point_response(axis, point) - position_response(axis))
                  END DO
               ELSE
                  CALL upf_inverse_kinetic_evaluate(kinetic(ikind), points(:, point) - position, kernel)
               END IF
               density(point) = density(point) + REAL(SUM(projected*kernel), dp)
            END DO
         END DO
         DEALLOCATE (projected, kernel)
         IF (derive) DEALLOCATE (projected_response, gradient)
      END DO
   END SUBROUTINE qs_upf_inverse_kinetic_collocate

! **************************************************************************************************
!> \brief Exact AO transpose of the same finite periodic kinetic-augmentation image sum.
!> \param frame ...
!> \param kinetic Original kind-indexed kernels.
!> \param points ...
!> \param weights Evaluation measure, including any partition, applied once.
!> \param potential Unweighted potential conjugate to the kinetic density.
!> \param ao_operator Hermitian AO operator without spin or mesh weights.
!> \param comm Optional communicator distributing point rows. Replicated points omit it.
!> \param ao_response Optional derivative of the same adjoint, requiring all response inputs.
!> \param point_response Cartesian motion of evaluation points.
!> \param atom_response Cartesian motion of original atoms.
!> \param cell_response Lattice response, including translated image centers.
!> \param coefficients_response Complete coupled inverse response dW.
!> \param weights_response Response of the integration measure.
!> \param potential_response Response of the unweighted kinetic potential.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_kinetic_operator(frame, kinetic, points, weights, potential, ao_operator, comm, &
                                              ao_response, point_response, atom_response, cell_response, &
                                              coefficients_response, weights_response, potential_response)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      TYPE(upf_inverse_kinetic_type), INTENT(IN)         :: kinetic(:)
      REAL(dp), INTENT(IN)                               :: points(:, :), weights(:), potential(:)
      COMPLEX(dp), INTENT(OUT)                           :: ao_operator(:, :)
      TYPE(mp_comm_type), INTENT(IN), OPTIONAL           :: comm
      COMPLEX(dp), INTENT(OUT), OPTIONAL                 :: ao_response(:, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: point_response(:, :), &
                                                            atom_response(:, :), &
                                                            cell_response(3, 3)
      COMPLEX(dp), INTENT(IN), OPTIONAL                  :: coefficients_response(:, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: weights_response(:), &
                                                            potential_response(:)

      INTEGER                                            :: atom, center, first, ikind, image, last, &
                                                            n, nao, np, point
      LOGICAL                                            :: derive
      REAL(dp)                                           :: position(3), position_response(3)
      REAL(dp), ALLOCATABLE                              :: kernel(:, :), kernel_gradient(:, :, :), &
                                                            OPERATOR(:, :), operator_response(:, :)

      nao = SIZE(frame%coefficients, 2)
      np = SIZE(points, 2)
      CPASSERT(frame%atom_local .OR. frame%inverse%ready)
      CPASSERT(SIZE(points, 1) == 3)
      CPASSERT(SIZE(weights) == np .AND. SIZE(potential) == np)
      CPASSERT(ALL(SHAPE(ao_operator) == [nao, nao]))
      derive = PRESENT(ao_response)
      CPASSERT(derive .EQV. PRESENT(point_response))
      CPASSERT(derive .EQV. PRESENT(atom_response))
      CPASSERT(derive .EQV. PRESENT(cell_response))
      CPASSERT(derive .EQV. PRESENT(coefficients_response))
      CPASSERT(derive .EQV. PRESENT(weights_response))
      CPASSERT(derive .EQV. PRESENT(potential_response))
      IF (derive) THEN
         CPASSERT(ALL(SHAPE(ao_response) == [nao, nao]))
         CPASSERT(ALL(SHAPE(point_response) == SHAPE(points)))
         CPASSERT(ALL(SHAPE(atom_response) == SHAPE(frame%positions)))
         CPASSERT(ALL(SHAPE(coefficients_response) == SHAPE(frame%coefficients)))
         CPASSERT(SIZE(weights_response) == np .AND. SIZE(potential_response) == np)
         ao_response = CMPLX(0.0_dp, 0.0_dp, dp)
      END IF
      ao_operator = CMPLX(0.0_dp, 0.0_dp, dp)
      DO center = 1, SIZE(frame%center_atoms)
         ikind = frame%center_kinds(center)
         CPASSERT(ikind >= 1 .AND. ikind <= SIZE(kinetic))
         IF (kinetic(ikind)%delta%nproj == 0) CYCLE
         CALL kinetic_center_indices(frame, kinetic(ikind), center, points, atom, first, last)
         n = last - first + 1
         ALLOCATE (kernel(n, n), OPERATOR(n, n))
         OPERATOR = 0.0_dp
         IF (derive) THEN
            ALLOCATE (kernel_gradient(n, n, 3), operator_response(n, n))
            operator_response = 0.0_dp
         END IF
         DO image = 1, SIZE(frame%images, 2)
            position = frame%positions(:, atom) + MATMUL(frame%cell, REAL(frame%images(:, image), dp))
            IF (derive) position_response = atom_response(:, atom) + &
                                            MATMUL(cell_response, REAL(frame%images(:, image), dp))
            DO point = 1, np
               IF (SUM((points(:, point) - position)**2) >= kinetic(ikind)%delta%radius**2) CYCLE
               IF (derive) THEN
                  CALL upf_inverse_kinetic_evaluate(kinetic(ikind), points(:, point) - position, kernel, kernel_gradient)
                  operator_response(:, :) = operator_response + &
                                    (weights_response(point)*potential(point) + weights(point)*potential_response(point))*kernel + &
                                            weights(point)*potential(point)*( &
                                            kernel_gradient(:, :, 1)*(point_response(1, point) - position_response(1)) + &
                                            kernel_gradient(:, :, 2)*(point_response(2, point) - position_response(2)) + &
                                            kernel_gradient(:, :, 3)*(point_response(3, point) - position_response(3)))
               ELSE
                  CALL upf_inverse_kinetic_evaluate(kinetic(ikind), points(:, point) - position, kernel)
               END IF
               OPERATOR(:, :) = OPERATOR + weights(point)*potential(point)*kernel
            END DO
         END DO
         ao_operator = ao_operator + MATMUL(CONJG(TRANSPOSE(frame%coefficients(first:last, :))), &
                                            MATMUL(OPERATOR, frame%coefficients(first:last, :)))
         IF (derive) THEN
            ao_response = ao_response + &
                          MATMUL(CONJG(TRANSPOSE(coefficients_response(first:last, :))), &
                                 MATMUL(OPERATOR, frame%coefficients(first:last, :))) + &
                          MATMUL(CONJG(TRANSPOSE(frame%coefficients(first:last, :))), &
                                 MATMUL(OPERATOR, coefficients_response(first:last, :))) + &
                          MATMUL(CONJG(TRANSPOSE(frame%coefficients(first:last, :))), &
                                 MATMUL(operator_response, frame%coefficients(first:last, :)))
            DEALLOCATE (kernel_gradient, operator_response)
         END IF
         DEALLOCATE (kernel, OPERATOR)
      END DO
      IF (PRESENT(comm)) CALL comm%sum(ao_operator)
      IF (PRESENT(comm) .AND. derive) CALL comm%sum(ao_response)
   END SUBROUTINE qs_upf_inverse_kinetic_operator

! **************************************************************************************************
!> \brief Home-center valence kinetic augmentation from the complete inverse projector density.
!>        Spin occupations and k weights are supplied by the caller. No core or charge Q is added.
!> \param frame Original defining inverse, including coupled coefficients W for this k point.
!> \param kinetic Independent reference-wave kernel of the selected original kind.
!> \param center Selected home center in the inverse frame.
!> \param points Cartesian points, with no partition or periodic image folding.
!> \param matrix Hermitian AO occupations for one spin and k point.
!> \param density Kinetic augmentation with the factor one half.
!> \param response Optional fixed-matrix response, requiring all three response inputs below.
!> \param point_response Cartesian motion of the evaluation points.
!> \param atom_response Cartesian motion of every original atom.
!> \param coefficients_response Complete dW, including any lattice or reduced-k response.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_kinetic_center(frame, kinetic, center, points, matrix, density, response, &
                                            point_response, atom_response, coefficients_response)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      TYPE(upf_inverse_kinetic_type), INTENT(IN)         :: kinetic
      INTEGER, INTENT(IN)                                :: center
      REAL(dp), INTENT(IN)                               :: points(:, :)
      COMPLEX(dp), INTENT(IN)                            :: matrix(:, :)
      REAL(dp), INTENT(OUT)                              :: density(:)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: response(:)
      REAL(dp), INTENT(IN), OPTIONAL                     :: point_response(:, :), atom_response(:, :)
      COMPLEX(dp), INTENT(IN), OPTIONAL                  :: coefficients_response(:, :)

      COMPLEX(dp), ALLOCATABLE                           :: projected(:, :), projected_response(:, :)
      INTEGER                                            :: atom, axis, first, last, n, nao, np, &
                                                            point
      LOGICAL                                            :: derive
      REAL(dp), ALLOCATABLE                              :: gradient(:, :, :), kernel(:, :)

      CALL kinetic_center_indices(frame, kinetic, center, points, atom, first, last)
      np = SIZE(points, 2)
      n = last - first + 1
      nao = SIZE(frame%coefficients, 2)
      CPASSERT(SIZE(density) == np)
      CPASSERT(ALL(SHAPE(matrix) == [nao, nao]))
      derive = PRESENT(response)
      CPASSERT(derive .EQV. PRESENT(point_response))
      CPASSERT(derive .EQV. PRESENT(atom_response))
      CPASSERT(derive .EQV. PRESENT(coefficients_response))
      ALLOCATE (projected(n, n), kernel(n, n))
    projected(:, :) = MATMUL(frame%coefficients(first:last, :), MATMUL(matrix, CONJG(TRANSPOSE(frame%coefficients(first:last, :)))))
      IF (derive) THEN
         CPASSERT(SIZE(response) == np .AND. ALL(SHAPE(point_response) == SHAPE(points)))
         CPASSERT(ALL(SHAPE(atom_response) == SHAPE(frame%positions)))
         CPASSERT(ALL(SHAPE(coefficients_response) == SHAPE(frame%coefficients)))
         ALLOCATE (projected_response(n, n), gradient(n, n, 3))
         projected_response(:, :) = MATMUL(coefficients_response(first:last, :), &
                                           MATMUL(matrix, CONJG(TRANSPOSE(frame%coefficients(first:last, :)))))
         projected_response(:, :) = projected_response + CONJG(TRANSPOSE(projected_response))
      END IF
      DO point = 1, np
         IF (derive) THEN
            CALL upf_inverse_kinetic_evaluate(kinetic, points(:, point) - frame%positions(:, atom), kernel, gradient)
            response(point) = REAL(SUM(projected_response*kernel), dp)
            DO axis = 1, 3
               response(point) = response(point) + REAL(SUM(projected*gradient(:, :, axis)), dp)* &
                                 (point_response(axis, point) - atom_response(axis, atom))
            END DO
         ELSE
            CALL upf_inverse_kinetic_evaluate(kinetic, points(:, point) - frame%positions(:, atom), kernel)
         END IF
         density(point) = REAL(SUM(projected*kernel), dp)
      END DO
   END SUBROUTINE qs_upf_inverse_kinetic_center

! **************************************************************************************************
!> \brief W^dagger K_tau[v_tau] W, the exact AO transpose of the home kinetic augmentation.
!> \param frame Same defining inverse as the density map.
!> \param kinetic Same original reference-wave kernel.
!> \param center Selected home center.
!> \param points Cartesian points, with no image folding or partition.
!> \param weights Point measure, applied once.
!> \param potential Unweighted real potential conjugate to the kinetic density.
!> \param ao_operator Hermitian AO operator, without spin or mesh weights.
!> \param comm Optional communicator distributing point rows. Replicated grids omit it.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_kinetic_integrate(frame, kinetic, center, points, weights, potential, ao_operator, comm)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      TYPE(upf_inverse_kinetic_type), INTENT(IN)         :: kinetic
      INTEGER, INTENT(IN)                                :: center
      REAL(dp), INTENT(IN)                               :: points(:, :), weights(:), potential(:)
      COMPLEX(dp), INTENT(OUT)                           :: ao_operator(:, :)
      TYPE(mp_comm_type), INTENT(IN), OPTIONAL           :: comm

      INTEGER                                            :: atom, first, last, n, nao, np, point
      REAL(dp), ALLOCATABLE                              :: kernel(:, :), OPERATOR(:, :)

      CALL kinetic_center_indices(frame, kinetic, center, points, atom, first, last)
      np = SIZE(points, 2)
      n = last - first + 1
      nao = SIZE(frame%coefficients, 2)
      CPASSERT(SIZE(weights) == np .AND. SIZE(potential) == np)
      CPASSERT(ALL(SHAPE(ao_operator) == [nao, nao]))
      ALLOCATE (kernel(n, n), OPERATOR(n, n))
      OPERATOR = 0.0_dp
      DO point = 1, np
         CALL upf_inverse_kinetic_evaluate(kinetic, points(:, point) - frame%positions(:, atom), kernel)
         OPERATOR(:, :) = OPERATOR + weights(point)*potential(point)*kernel
      END DO
      ao_operator = MATMUL(CONJG(TRANSPOSE(frame%coefficients(first:last, :))), &
                           MATMUL(OPERATOR, frame%coefficients(first:last, :)))
      IF (PRESENT(comm)) CALL comm%sum(ao_operator)
   END SUBROUTINE qs_upf_inverse_kinetic_integrate

! **************************************************************************************************
!> \brief Check the selected original projector block shared by kinetic density and its AO transpose.
!> \param frame ...
!> \param kinetic ...
!> \param center ...
!> \param points ...
!> \param atom ...
!> \param first ...
!> \param last ...
! **************************************************************************************************
   SUBROUTINE kinetic_center_indices(frame, kinetic, center, points, atom, first, last)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      TYPE(upf_inverse_kinetic_type), INTENT(IN)         :: kinetic
      INTEGER, INTENT(IN)                                :: center
      REAL(dp), INTENT(IN)                               :: points(:, :)
      INTEGER, INTENT(OUT)                               :: atom, first, last

      CPASSERT(center >= 1 .AND. center <= SIZE(frame%center_atoms))
      CPASSERT(SIZE(points, 1) == 3)
      atom = frame%center_atoms(center)
      CPASSERT(atom >= 1 .AND. atom <= SIZE(frame%positions, 2))
      first = frame%projector_first(center)
      last = frame%projector_last(center)
      CPASSERT(first >= 1 .AND. last <= SIZE(frame%coefficients, 1))
      CPASSERT(last - first + 1 == kinetic%delta%nproj)
   END SUBROUTINE kinetic_center_indices

! **************************************************************************************************
!> \brief Complete signed GAPW center density from the same accepted AO occupation matrix.
!>        Combines native Gaussian hard-minus-soft, inverse reconstruction and original UPF Q.
!>        Native atom grids have fixed Cartesian radial offsets. No partition or core is added.
!> \param frame Original defining inverse, including the full coupled W at this k point.
!> \param augmentation Original charge kernels.
!> \param center Selected home center in the inverse frame.
!> \param kind Borrowed native kind of this center, including its GAPW_1C basis and atom grid.
!> \param pairs Prepared native OCE pair map for this same center and geometry.
!> \param matrix Hermitian AO occupations for one spin and k point, without a mesh weight.
!> \param radial Complete signed harmonic density, without the radial measure.
!> \param response Optional fixed-matrix atom/cell/Bloch response of this complete center density.
!> \param atom_response Cartesian motion of all original atoms.
!> \param cell_response Physical lattice response.
!> \param reduced_k_response Reduced Bloch-vector response.
!> \param coefficients_response Complete dW from the original defining inverse quadrature.
!> \param atom_index Native Gaussian center, required when center is zero and no inverse correction is present.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_density_gapw(frame, augmentation, center, kind, pairs, matrix, radial, &
                                          response, atom_response, cell_response, reduced_k_response, &
                                          coefficients_response, atom_index)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      TYPE(qs_upf_inverse_augmentation_type), INTENT(IN) :: augmentation
      INTEGER, INTENT(IN)                                :: center
      TYPE(qs_kind_type), INTENT(IN)                     :: kind
      TYPE(qs_upf_inverse_gaussian_pairs_type), &
         INTENT(IN)                                      :: pairs
      COMPLEX(dp), INTENT(IN)                            :: matrix(:, :)
      REAL(dp), INTENT(OUT)                              :: radial(:, :)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: response(:, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: atom_response(:, :), &
                                                            cell_response(3, 3), &
                                                            reduced_k_response(3)
      COMPLEX(dp), INTENT(IN), OPTIONAL                  :: coefficients_response(:, :)
      INTEGER, INTENT(IN), OPTIONAL                      :: atom_index

      INTEGER                                            :: atom, n, nh, np
      LOGICAL                                            :: derive
      REAL(dp), ALLOCATABLE :: deh(:, :), des(:, :), dh(:, :), ds(:, :), dvalues(:), eh(:, :), &
         es(:, :), hard(:, :), point_response(:, :), points(:, :), soft(:, :), values(:), zero(:)
      TYPE(grid_atom_type), POINTER                      :: grid
      TYPE(harmonics_atom_type), POINTER                 :: harmonics

      CALL gapw_grid(frame, center, kind, grid, harmonics, n, points, atom_index)
      IF (center > 0) CALL check_center(frame, augmentation, center, points)
      nh = harmonics%max_iso_not0
      np = SIZE(points, 2)
      CPASSERT(ALL(SHAPE(radial) == [grid%nr, nh]))
      derive = PRESENT(response)
      CPASSERT(derive .EQV. PRESENT(atom_response))
      CPASSERT(derive .EQV. PRESENT(cell_response))
      CPASSERT(derive .EQV. PRESENT(reduced_k_response))
      CPASSERT(derive .EQV. PRESENT(coefficients_response))
      ALLOCATE (hard(n, n), soft(n, n), eh(grid%nr, nh), es(grid%nr, nh), values(np))
      IF (derive) THEN
         CPASSERT(ALL(SHAPE(response) == SHAPE(radial)))
         CPASSERT(ALL(SHAPE(atom_response) == SHAPE(frame%positions)))
         IF (PRESENT(atom_index)) THEN
            atom = atom_index
         ELSE
            atom = frame%center_atoms(center)
         END IF
         ALLOCATE (dh(n, n), ds(n, n), deh(grid%nr, nh), des(grid%nr, nh), &
                   dvalues(np), point_response(3, np), zero(grid%nr))
         point_response(:, :) = SPREAD(atom_response(:, atom), 2, np)
         zero = 0.0_dp
         CALL qs_upf_inverse_gaussian_pairs_density(pairs, frame%reduced_k, matrix, hard, soft, &
                                                    dh, ds, atom_response, cell_response, reduced_k_response)
         CALL qs_upf_inverse_gaussian_radial(kind, hard, eh, deh, dh, zero)
         CALL qs_upf_inverse_gaussian_radial(kind, soft, es, des, ds, zero)
         dvalues = 0.0_dp
         IF (center > 0) THEN
            CALL qs_upf_inverse_density_center(frame, augmentation, center, points, matrix, values, dvalues, &
                                               point_response, atom_response, cell_response, reduced_k_response, &
                                               coefficients_response)
         END IF
         CALL qs_upf_inverse_density_radial(grid, harmonics, dvalues, response, deh, des)
      ELSE
         CALL qs_upf_inverse_gaussian_pairs_density(pairs, frame%reduced_k, matrix, hard, soft)
         CALL qs_upf_inverse_gaussian_radial(kind, hard, eh)
         CALL qs_upf_inverse_gaussian_radial(kind, soft, es)
         values = 0.0_dp
         IF (center > 0) CALL qs_upf_inverse_density_center(frame, augmentation, center, points, matrix, values)
      END IF
      IF (center == 0) values = 0.0_dp
      CALL qs_upf_inverse_density_radial(grid, harmonics, values, radial, eh, es)
   END SUBROUTINE qs_upf_inverse_density_gapw

! **************************************************************************************************
!> \brief AO transpose of the complete finite GAPW center density, including every density leg.
!>        Setting the monopole covector to sqrt(4*pi) gives this center's charge metric correction.
!>        The complete metric must also include the Gaussian background's unit-potential adjoint.
!> \param frame ...
!> \param augmentation ...
!> \param center ...
!> \param kind ...
!> \param pairs ...
!> \param potential Harmonic radial-measure covector, with no radial weight or k weight included.
!> \param ao_operator Hermitian AO center contribution, without a mesh weight.
!> \param atom_index Native Gaussian center, required for center zero without an inverse correction.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_density_gapw_integrate(frame, augmentation, center, kind, pairs, potential, ao_operator, atom_index)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      TYPE(qs_upf_inverse_augmentation_type), INTENT(IN) :: augmentation
      INTEGER, INTENT(IN)                                :: center
      TYPE(qs_kind_type), INTENT(IN)                     :: kind
      TYPE(qs_upf_inverse_gaussian_pairs_type), &
         INTENT(IN)                                      :: pairs
      REAL(dp), INTENT(IN)                               :: potential(:, :)
      COMPLEX(dp), INTENT(OUT)                           :: ao_operator(:, :)
      INTEGER, INTENT(IN), OPTIONAL                      :: atom_index

      COMPLEX(dp), ALLOCATABLE                           :: gaussian(:, :)
      INTEGER                                            :: n, nh, np
      REAL(dp), ALLOCATABLE                              :: bh(:, :), bs(:, :), hard(:, :), &
                                                            points(:, :), soft(:, :), values(:), &
                                                            weights(:)
      TYPE(grid_atom_type), POINTER                      :: grid
      TYPE(harmonics_atom_type), POINTER                 :: harmonics

      CALL gapw_grid(frame, center, kind, grid, harmonics, n, points, atom_index)
      IF (center > 0) CALL check_center(frame, augmentation, center, points)
      nh = harmonics%max_iso_not0
      np = SIZE(points, 2)
      CPASSERT(ALL(SHAPE(potential) == [grid%nr, nh]))
      ALLOCATE (bh(n, n), bs(n, n), hard(grid%nr, nh), soft(grid%nr, nh), values(np), &
                weights(np), gaussian(SIZE(ao_operator, 1), SIZE(ao_operator, 2)))
      CALL qs_upf_inverse_density_radial_adjoint(grid, harmonics, potential, values, hard, soft)
      CALL qs_upf_inverse_gaussian_integrate(kind, hard, bh)
      CALL qs_upf_inverse_gaussian_integrate(kind, soft, bs)
      CALL qs_upf_inverse_gaussian_pairs_integrate(pairs, frame%reduced_k, bh, bs, gaussian)
      weights(:) = RESHAPE(grid%weight, [np])
      ! Atom rows are replicated. Only the Gaussian OCE pair transpose performs its native pair sum.
      ao_operator = CMPLX(0.0_dp, 0.0_dp, dp)
      IF (center > 0) THEN
         CALL qs_upf_inverse_density_integrate(frame, augmentation, center, points, weights, values, ao_operator)
      END IF
      ao_operator = ao_operator + gaussian
   END SUBROUTINE qs_upf_inverse_density_gapw_integrate

! **************************************************************************************************
!> \brief Borrow the same native layout for both sides of the complete GAPW density map.
!> \param frame ...
!> \param center ...
!> \param kind ...
!> \param grid Borrowed native radial/Lebedev grid.
!> \param harmonics Borrowed native harmonics and Gaunt coefficients.
!> \param n Number of compressed Gaussian primitives.
!> \param points Replicated Cartesian atom-grid points with angular rows fastest.
!> \param atom_index Explicit native Gaussian center. Center zero has no inverse or Q contribution.
! **************************************************************************************************
   SUBROUTINE gapw_grid(frame, center, kind, grid, harmonics, n, points, atom_index)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      INTEGER, INTENT(IN)                                :: center
      TYPE(qs_kind_type), INTENT(IN)                     :: kind
      TYPE(grid_atom_type), POINTER                      :: grid
      TYPE(harmonics_atom_type), POINTER                 :: harmonics
      INTEGER, INTENT(OUT)                               :: n
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: points(:, :)
      INTEGER, INTENT(IN), OPTIONAL                      :: atom_index

      INTEGER                                            :: atom, ia, ir, row
      TYPE(gto_basis_set_type), POINTER                  :: basis

      CPASSERT(center >= 0 .AND. center <= SIZE(frame%center_atoms))
      IF (PRESENT(atom_index)) THEN
         atom = atom_index
         IF (center > 0) THEN
            CPASSERT(atom == frame%center_atoms(center))
         END IF
      ELSE
         CPASSERT(center > 0)
         atom = frame%center_atoms(center)
      END IF
      CPASSERT(atom >= 1 .AND. atom <= SIZE(frame%positions, 2))
      CALL get_qs_kind(kind, basis_set=basis, basis_type="GAPW_1C", grid_atom=grid, harmonics=harmonics)
      CPASSERT(ASSOCIATED(basis) .AND. ASSOCIATED(grid) .AND. ASSOCIATED(harmonics))
      CPASSERT(ASSOCIATED(grid%rad) .AND. ASSOCIATED(grid%weight))
      CPASSERT(grid%nr > 0 .AND. grid%ng_sphere > 0)
      CALL get_paw_basis_info(basis, nsatbas=n)
      ALLOCATE (points(3, grid%nr*grid%ng_sphere))
      DO ir = 1, grid%nr
         DO ia = 1, grid%ng_sphere
            row = ia + (ir - 1)*grid%ng_sphere
            points(:, row) = frame%positions(:, atom) + grid%rad(ir)* &
                             [grid%sin_pol(ia)*grid%cos_azi(ia), grid%sin_pol(ia)*grid%sin_azi(ia), grid%cos_pol(ia)]
         END DO
      END DO
   END SUBROUTINE gapw_grid

! **************************************************************************************************
!> \brief Signed nR_a + Q_a from the full coupled inverse, without partition or image folding.
!>        nR_a = Re[d_a^* (s + chi) P], with d_a = -Delta_a W_a for the home image.
!>        Sum caller-supplied spin occupations and k weights outside this map exactly once.
!> \param frame Original inverse frame, including all centers and their coupled W.
!> \param augmentation Original charge kernels in the frame's magnetic-projector order.
!> \param center Selected home center, including a possible NC center with zero reconstruction/Q.
!> \param points Cartesian evaluation points, possibly an empty batch.
!> \param matrix Hermitian AO occupation density for one spin and k point.
!> \param density Signed home-center charge density, excluding Gaussian hard/soft and core.
!> \param response Optional joint density response at fixed matrix. All response inputs are required.
!> \param point_response Motion of the evaluation points.
!> \param atom_response Motion of every original atom.
!> \param cell_response Physical lattice response.
!> \param reduced_k_response Reduced Bloch-vector response.
!> \param coefficients_response Complete dW from the original defining inverse quadrature.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_density_center(frame, augmentation, center, points, matrix, density, response, &
                                            point_response, atom_response, cell_response, reduced_k_response, &
                                            coefficients_response)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      TYPE(qs_upf_inverse_augmentation_type), INTENT(IN) :: augmentation
      INTEGER, INTENT(IN)                                :: center
      REAL(dp), INTENT(IN)                               :: points(:, :)
      COMPLEX(dp), INTENT(IN)                            :: matrix(:, :)
      REAL(dp), INTENT(OUT)                              :: density(:)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: response(:)
      REAL(dp), INTENT(IN), OPTIONAL                     :: point_response(:, :), &
                                                            atom_response(:, :), &
                                                            cell_response(3, 3), &
                                                            reduced_k_response(3)
      COMPLEX(dp), INTENT(IN), OPTIONAL                  :: coefficients_response(:, :)

      COMPLEX(dp), ALLOCATABLE                           :: d(:, :), dd(:, :), de(:, :), e(:, :), &
                                                            gd(:, :, :), ge(:, :, :), &
                                                            projected(:, :), &
                                                            projected_response(:, :)
      INTEGER                                            :: atom, first, kind, last, n, np, nproj
      LOGICAL                                            :: derive
      REAL(dp), ALLOCATABLE                              :: dq(:), q(:), qgradient(:, :)

      CALL check_center(frame, augmentation, center, points)
      np = SIZE(points, 2)
      n = SIZE(frame%coefficients, 2)
      nproj = SIZE(frame%coefficients, 1)
      CPASSERT(SIZE(density) == np)
      CPASSERT(ALL(SHAPE(matrix) == [n, n]))
      derive = PRESENT(response)
      CPASSERT(derive .EQV. PRESENT(point_response))
      CPASSERT(derive .EQV. PRESENT(atom_response))
      CPASSERT(derive .EQV. PRESENT(cell_response))
      CPASSERT(derive .EQV. PRESENT(reduced_k_response))
      CPASSERT(derive .EQV. PRESENT(coefficients_response))
      ALLOCATE (d(np, n), e(np, n), gd(np, n, 3), ge(np, n, 3), &
                q(np), qgradient(3, np), projected(nproj, nproj))
      IF (derive) THEN
         CPASSERT(SIZE(response) == np)
         ALLOCATE (dd(np, n), de(np, n), dq(np))
         CALL center_fields(frame, center, points, d, e, gd, ge, dd, de, &
                            point_response=point_response, atom_response=atom_response, cell_response=cell_response, &
                            reduced_k_response=reduced_k_response, coefficients_response=coefficients_response)
         CALL upf_inverse_cross_density(dd, e, ao_density=matrix, density=response)
         CALL upf_inverse_cross_density(d, de, ao_density=matrix, density=dq)
         response = response + dq
      ELSE
         CALL center_fields(frame, center, points, d, e, gd, ge)
      END IF
      CALL upf_inverse_cross_density(d, e, ao_density=matrix, density=density)
      atom = frame%center_atoms(center)
      kind = frame%center_kinds(center)
      IF (augmentation%kinds(kind)%nproj == 0) RETURN
      first = frame%projector_first(center)
      last = frame%projector_last(center)
      CALL upf_inverse_projector_density(frame%coefficients, matrix, projected)
      q = 0.0_dp
      qgradient = 0.0_dp
      CALL uspp_augmentation_collocate_atom(augmentation%kinds(kind), &
                                            points - SPREAD(frame%positions(:, atom), 2, np), &
                                            REAL(projected(first:last, first:last), dp), q, qgradient)
      density = density + q
      IF (derive) THEN
         response = response + SUM(qgradient*(point_response - SPREAD(atom_response(:, atom), 2, np)), DIM=1)
         ALLOCATE (projected_response(nproj, nproj))
         projected_response(:, :) = MATMUL(coefficients_response, MATMUL(matrix, CONJG(TRANSPOSE(frame%coefficients))))
         projected_response(:, :) = projected_response + CONJG(TRANSPOSE(projected_response))
         q = 0.0_dp
         qgradient = 0.0_dp
         CALL uspp_augmentation_collocate_atom(augmentation%kinds(kind), &
                                               points - SPREAD(frame%positions(:, atom), 2, np), &
                                               REAL(projected_response(first:last, first:last), dp), q, qgradient)
         response = response + q
      END IF
   END SUBROUTINE qs_upf_inverse_density_center

! **************************************************************************************************
!> \brief Density adjoint with both inverse field legs and the original W^dagger Q W operator.
!>        A unit potential gives this center's contribution to the finite integrated charge metric.
!> \param frame ...
!> \param augmentation ...
!> \param center ...
!> \param points ...
!> \param weights Point measure, applied once.
!> \param potential Unweighted real density covector.
!> \param ao_operator Hermitian AO adjoint for one spin and k point, without a k weight.
!> \param comm Optional communicator distributing point rows. Replicated atom grids omit it.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_density_integrate(frame, augmentation, center, points, weights, potential, ao_operator, comm)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      TYPE(qs_upf_inverse_augmentation_type), INTENT(IN) :: augmentation
      INTEGER, INTENT(IN)                                :: center
      REAL(dp), INTENT(IN)                               :: points(:, :), weights(:), potential(:)
      COMPLEX(dp), INTENT(OUT)                           :: ao_operator(:, :)
      TYPE(mp_comm_type), INTENT(IN), OPTIONAL           :: comm

      COMPLEX(dp), ALLOCATABLE                           :: d(:, :), e(:, :), gd(:, :, :), &
                                                            ge(:, :, :), qoperator(:, :)
      INTEGER                                            :: atom, first, kind, last, n, np
      REAL(dp), ALLOCATABLE                              :: q(:, :), zero(:), zeros(:, :)

      CALL check_center(frame, augmentation, center, points)
      np = SIZE(points, 2)
      n = SIZE(frame%coefficients, 2)
      CPASSERT(SIZE(weights) == np .AND. SIZE(potential) == np)
      CPASSERT(ALL(SHAPE(ao_operator) == [n, n]))
      ALLOCATE (d(np, n), e(np, n), gd(np, n, 3), ge(np, n, 3), zero(np), zeros(np, 3))
      CALL center_fields(frame, center, points, d, e, gd, ge)
      zero = 0.0_dp
      zeros = 0.0_dp
      CALL upf_inverse_semilocal_pair(d, e, gd, ge, weights, potential, zeros, zero, ao_operator)
      ao_operator = 0.5_dp*(ao_operator + CONJG(TRANSPOSE(ao_operator)))
      atom = frame%center_atoms(center)
      kind = frame%center_kinds(center)
      IF (augmentation%kinds(kind)%nproj /= 0) THEN
         first = frame%projector_first(center)
         last = frame%projector_last(center)
         ALLOCATE (q(last - first + 1, last - first + 1), qoperator(last - first + 1, last - first + 1))
         CALL uspp_augmentation_integrate_atom(augmentation%kinds(kind), &
                                               points - SPREAD(frame%positions(:, atom), 2, np), &
                                               weights, potential, TRANSPOSE(zeros), q)
         qoperator(:, :) = CMPLX(q, 0.0_dp, dp)
         ao_operator = ao_operator + MATMUL(CONJG(TRANSPOSE(frame%coefficients(first:last, :))), &
                                            MATMUL(qoperator, frame%coefficients(first:last, :)))
      END IF
      IF (PRESENT(comm)) CALL comm%sum(ao_operator)
   END SUBROUTINE qs_upf_inverse_density_integrate

! **************************************************************************************************
!> \brief Native angular projection of nR + Q and optional Gaussian hard-minus-soft density.
!>        Use the same finite Lebedev map for every leg. Its Gram matrix is not assumed to be unity.
!>        The native harmonic order is l=0,1,... and m=-l,...,l, as in the Coulomb maps.
!> \param grid Borrowed native radial/angular quadrature, with complete replicated angular rows.
!> \param harmonics Borrowed native real harmonics on that same angular quadrature.
!> \param values Point density with angular rows fastest. No core or partition is added here.
!> \param radial New real harmonic density coefficients, without the radial measure.
!> \param gaussian_hard Optional native Gaussian hard harmonic coefficients.
!> \param gaussian_soft Matching native Gaussian soft harmonic coefficients.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_density_radial(grid, harmonics, values, radial, gaussian_hard, gaussian_soft)
      TYPE(grid_atom_type), INTENT(IN)                   :: grid
      TYPE(harmonics_atom_type), INTENT(IN)              :: harmonics
      REAL(dp), INTENT(IN)                               :: values(:)
      REAL(dp), INTENT(OUT)                              :: radial(:, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: gaussian_hard(:, :), gaussian_soft(:, :)

      INTEGER                                            :: first, ir, last, na, nh
      REAL(dp), ALLOCATABLE                              :: angular(:)

      na = grid%ng_sphere
      nh = SIZE(radial, 2)
      CALL check_radial(grid, harmonics, values, radial)
      CPASSERT(PRESENT(gaussian_hard) .EQV. PRESENT(gaussian_soft))
      IF (PRESENT(gaussian_hard)) THEN
         CPASSERT(ALL(SHAPE(gaussian_hard) == SHAPE(radial)))
         CPASSERT(ALL(SHAPE(gaussian_soft) == SHAPE(radial)))
      END IF
      ALLOCATE (angular(na))
      DO ir = 1, grid%nr
         first = (ir - 1)*na + 1
         last = ir*na
         angular(:) = values(first:last)
         IF (PRESENT(gaussian_hard)) THEN
            angular(:) = angular + MATMUL(harmonics%slm(:, 1:nh), gaussian_hard(ir, :) - gaussian_soft(ir, :))
         END IF
         radial(ir, :) = MATMUL(grid%wa*angular, harmonics%slm(:, 1:nh))
      END DO
   END SUBROUTINE qs_upf_inverse_density_radial

! **************************************************************************************************
!> \brief Exact transpose of the native finite angular projection for all three density legs.
!> \param grid ...
!> \param harmonics ...
!> \param potential Unweighted radial-measure covector, as returned by the Coulomb assembly.
!> \param values Unweighted point covector. Multiply by grid%weight once in point/AO integration.
!> \param gaussian_hard Optional native hard radial-measure covector, with the angular Gram map.
!> \param gaussian_soft Matching native soft radial-measure covector with opposite sign.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_density_radial_adjoint(grid, harmonics, potential, values, gaussian_hard, gaussian_soft)
      TYPE(grid_atom_type), INTENT(IN)                   :: grid
      TYPE(harmonics_atom_type), INTENT(IN)              :: harmonics
      REAL(dp), INTENT(IN)                               :: potential(:, :)
      REAL(dp), INTENT(OUT)                              :: values(:)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: gaussian_hard(:, :), gaussian_soft(:, :)

      INTEGER                                            :: first, ir, last, na, nh

      CALL check_radial(grid, harmonics, values, potential)
      CPASSERT(PRESENT(gaussian_hard) .EQV. PRESENT(gaussian_soft))
      IF (PRESENT(gaussian_hard)) THEN
         CPASSERT(ALL(SHAPE(gaussian_hard) == SHAPE(potential)))
         CPASSERT(ALL(SHAPE(gaussian_soft) == SHAPE(potential)))
      END IF
      na = grid%ng_sphere
      nh = SIZE(potential, 2)
      DO ir = 1, grid%nr
         first = (ir - 1)*na + 1
         last = ir*na
         values(first:last) = MATMUL(harmonics%slm(:, 1:nh), potential(ir, :))
         IF (PRESENT(gaussian_hard)) THEN
            gaussian_hard(ir, :) = MATMUL(grid%wa*values(first:last), harmonics%slm(:, 1:nh))
            gaussian_soft(ir, :) = -gaussian_hard(ir, :)
         END IF
      END DO
   END SUBROUTINE qs_upf_inverse_density_radial_adjoint

! **************************************************************************************************
!> \brief Full environmental field s+chi and a single home-center difference, including responses.
!> \param frame ...
!> \param center ...
!> \param points ...
!> \param d ...
!> \param e ...
!> \param gd ...
!> \param ge ...
!> \param dd ...
!> \param de ...
!> \param gdd ...
!> \param gde ...
!> \param point_response ...
!> \param atom_response ...
!> \param cell_response ...
!> \param reduced_k_response ...
!> \param coefficients_response ...
! **************************************************************************************************
   SUBROUTINE center_fields(frame, center, points, d, e, gd, ge, dd, de, gdd, gde, &
                            point_response, atom_response, cell_response, reduced_k_response, coefficients_response)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      INTEGER, INTENT(IN)                                :: center
      REAL(dp), INTENT(IN)                               :: points(:, :)
      COMPLEX(dp), INTENT(OUT)                           :: d(:, :), e(:, :), gd(:, :, :), &
                                                            ge(:, :, :)
      COMPLEX(dp), INTENT(OUT), OPTIONAL                 :: dd(:, :), de(:, :), gdd(:, :, :), &
                                                            gde(:, :, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: point_response(:, :), &
                                                            atom_response(:, :), &
                                                            cell_response(3, 3), &
                                                            reduced_k_response(3)
      COMPLEX(dp), INTENT(IN), OPTIONAL                  :: coefficients_response(:, :)

      COMPLEX(dp), ALLOCATABLE                           :: chi(:, :), dchi(:, :), dgc(:, :, :), &
                                                            gc(:, :, :)
      INTEGER                                            :: atom, first, last, n, np

      n = SIZE(frame%coefficients, 2)
      np = SIZE(points, 2)
      ALLOCATE (chi(np, n), gc(np, n, 3))
      CALL qs_upf_inverse_center_fields(frame, center, points, d, gd, dd, gdd, point_response, &
                                        atom_response, cell_response, reduced_k_response, coefficients_response, [0, 0, 0])
      CALL qs_upf_inverse_evaluate(frame, points, e, ge, de, gde, point_response, atom_response, &
                                   cell_response, reduced_k_response, coefficients_response)
      IF (PRESENT(de)) ALLOCATE (dchi(np, n))
      IF (PRESENT(gde)) ALLOCATE (dgc(np, n, 3))
      DO atom = 1, SIZE(frame%bases)
         first = frame%ao_first(atom)
         last = frame%ao_last(atom)
         IF (PRESENT(gde)) THEN
            CALL gaussian_basis_fields_periodic(frame%bases(atom)%gto_basis_set, points, frame%positions(:, atom), &
                                                frame%cell, frame%images, frame%reduced_k, &
                                                chi(:, first:last), gc(:, first:last, :), &
                                                dchi(:, first:last), dgc(:, first:last, :), &
                                                point_response, atom_response(:, atom), &
                                                cell_response, reduced_k_response)
         ELSE IF (PRESENT(de)) THEN
            CALL gaussian_basis_fields_periodic(frame%bases(atom)%gto_basis_set, points, frame%positions(:, atom), &
                                                frame%cell, frame%images, frame%reduced_k, chi(:, first:last), &
                                                gc(:, first:last, :), response=dchi(:, first:last), &
                                                point_response=point_response, center_response=atom_response(:, atom), &
                                                cell_response=cell_response, reduced_k_response=reduced_k_response)
         ELSE
            CALL gaussian_basis_fields_periodic(frame%bases(atom)%gto_basis_set, points, frame%positions(:, atom), &
                                                frame%cell, frame%images, frame%reduced_k, chi(:, first:last), gc(:, first:last, :))
         END IF
      END DO
      e = e + chi
      ge = ge + gc
      IF (PRESENT(de)) THEN
         de = de + dchi
         IF (PRESENT(gde)) gde = gde + dgc
      END IF
   END SUBROUTINE center_fields

! **************************************************************************************************
!> \brief Verify the original center and magnetic charge-kernel order.
!> \param frame ...
!> \param augmentation ...
!> \param center ...
!> \param points ...
! **************************************************************************************************
   SUBROUTINE check_center(frame, augmentation, center, points)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      TYPE(qs_upf_inverse_augmentation_type), INTENT(IN) :: augmentation
      INTEGER, INTENT(IN)                                :: center
      REAL(dp), INTENT(IN)                               :: points(:, :)

      INTEGER                                            :: kind, nproj

      CPASSERT(frame%atom_local .OR. frame%inverse%ready)
      CPASSERT(center >= 1 .AND. center <= SIZE(frame%center_atoms))
      CPASSERT(SIZE(points, 1) == 3)
      kind = frame%center_kinds(center)
      CPASSERT(kind >= 1 .AND. kind <= SIZE(augmentation%kinds))
      nproj = frame%projector_last(center) - frame%projector_first(center) + 1
      IF (augmentation%kinds(kind)%nproj /= 0) THEN
         CPASSERT(augmentation%kinds(kind)%nproj == nproj)
      END IF
   END SUBROUTINE check_center

! **************************************************************************************************
!> \brief Validate one complete replicated radial/angular field in native harmonic order.
!> \param grid ...
!> \param harmonics ...
!> \param values ...
!> \param radial ...
! **************************************************************************************************
   SUBROUTINE check_radial(grid, harmonics, values, radial)
      TYPE(grid_atom_type), INTENT(IN)                   :: grid
      TYPE(harmonics_atom_type), INTENT(IN)              :: harmonics
      REAL(dp), INTENT(IN)                               :: values(:), radial(:, :)

      CPASSERT(grid%nr >= 0 .AND. grid%ng_sphere > 0)
      CPASSERT(ASSOCIATED(grid%wa) .AND. ASSOCIATED(harmonics%slm))
      CPASSERT(SIZE(values) == grid%nr*grid%ng_sphere)
      CPASSERT(SIZE(radial, 1) == grid%nr)
      CPASSERT(SIZE(grid%wa) == grid%ng_sphere)
      CPASSERT(SIZE(harmonics%slm, 1) == grid%ng_sphere)
      CPASSERT(SIZE(radial, 2) <= SIZE(harmonics%slm, 2))
   END SUBROUTINE check_radial

END MODULE qs_upf_inverse_density
