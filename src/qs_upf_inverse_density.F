!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
!> \brief Home-center inverse density and its finite Gaussian GAPW angular map.
MODULE qs_upf_inverse_density
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE message_passing,                 ONLY: mp_comm_type
   USE orbital_pointers,                ONLY: indso
   USE paw_basis_types,                 ONLY: get_paw_basis_info
   USE qs_grid_atom,                    ONLY: grid_atom_type
   USE qs_harmonics_atom,               ONLY: harmonics_atom_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_upf_inverse,                  ONLY: qs_upf_inverse_block_size,&
                                              qs_upf_inverse_center_fields,&
                                              qs_upf_inverse_evaluate
   USE qs_upf_inverse_atom_grid,        ONLY: qs_upf_atom_angular_average
   USE qs_upf_inverse_augmentation,     ONLY: qs_upf_inverse_augmentation_has_support
   USE qs_upf_inverse_gaussian,         ONLY: qs_upf_inverse_gaussian_integrate,&
                                              qs_upf_inverse_gaussian_radial
   USE qs_upf_inverse_gaussian_pairs,   ONLY: qs_upf_inverse_gaussian_pairs_density,&
                                              qs_upf_inverse_gaussian_pairs_integrate,&
                                              qs_upf_inverse_gaussian_pairs_type
   USE qs_upf_inverse_types,            ONLY: qs_upf_inverse_augmentation_type,&
                                              qs_upf_inverse_density_cache_type,&
                                              qs_upf_inverse_type
   USE spherical_harmonics,             ONLY: y_lm
   USE upf_inverse,                     ONLY: upf_inverse_kinetic_evaluate,&
                                              upf_inverse_kinetic_type
   USE upf_inverse_quadrature,          ONLY: upf_inverse_cross_density,&
                                              upf_inverse_pair,&
                                              upf_inverse_pair_response,&
                                              upf_inverse_projector_density,&
                                              upf_inverse_semilocal_pair
   USE uspp_augmentation,               ONLY: uspp_augmentation_collocate_atom,&
                                              uspp_augmentation_integrate_atom
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   TYPE qs_upf_inverse_gapw_view_type
      TYPE(grid_atom_type) :: grid
      TYPE(harmonics_atom_type) :: harmonics
      REAL(dp), ALLOCATABLE :: angular_weights(:), weights(:, :), slm(:, :)
   END TYPE qs_upf_inverse_gapw_view_type

   PUBLIC :: qs_upf_inverse_gapw_grid, qs_upf_inverse_gapw_view_type
   PUBLIC :: qs_upf_inverse_density_center, qs_upf_inverse_density_integrate
   PUBLIC :: qs_upf_inverse_density_pair_integrate
   PUBLIC :: qs_upf_inverse_density_radial, qs_upf_inverse_density_radial_adjoint
   PUBLIC :: qs_upf_inverse_density_gapw, qs_upf_inverse_density_gapw_integrate
   PUBLIC :: qs_upf_inverse_density_cache_prepare
   PUBLIC :: qs_upf_inverse_kinetic_center, qs_upf_inverse_kinetic_integrate
   PUBLIC :: qs_upf_inverse_kinetic_collocate, qs_upf_inverse_kinetic_operator

CONTAINS

! **************************************************************************************************
!> \brief Store immutable home-center fields within the caller's remaining byte budget.
!>        An empty cache retains streamed evaluation. No density matrix or geometry response is stored.
!> \param frame Original defining inverse and borrowed live Gaussian basis owner.
!> \param center Home center in this frame. Zero and DIRECT need no reconstruction cache.
!> \param kind Native owner of the same replicated radial/angular quadrature.
!> \param remaining_bytes Shared budget for all centers and Bloch points in one Coulomb map.
!> \param cache Owned fields, released with that map before the frame or its bases change.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_density_cache_prepare(frame, center, kind, remaining_bytes, cache)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      INTEGER, INTENT(IN)                                :: center
      TYPE(qs_kind_type), INTENT(IN)                     :: kind
      INTEGER(int_8), INTENT(INOUT)                      :: remaining_bytes
      TYPE(qs_upf_inverse_density_cache_type), &
         INTENT(OUT)                                     :: cache

      INTEGER                                            :: n, nao, np
      INTEGER(int_8)                                     :: bytes, npoint
      REAL(dp), ALLOCATABLE                              :: points(:, :)
      TYPE(grid_atom_type), POINTER                      :: grid
      TYPE(harmonics_atom_type), POINTER                 :: harmonics
      TYPE(qs_upf_inverse_gapw_view_type), TARGET         :: view

      CPASSERT(remaining_bytes >= 0_int_8)
      IF (center == 0 .OR. frame%direct) RETURN
      IF (remaining_bytes == 0_int_8) RETURN
      CALL qs_upf_inverse_gapw_grid(frame, center, kind, grid, harmonics, n, points, view)
      CPASSERT(ASSOCIATED(grid))
      CPASSERT(grid%nr > 0 .AND. grid%ng_sphere > 0)
      npoint = INT(grid%nr, int_8)*grid%ng_sphere
      nao = SIZE(frame%coefficients, 2)
      bytes = 16_int_8*(8_int_8*npoint*nao + SIZE(frame%coefficients, KIND=int_8)) + &
              8_int_8*(3_int_8*npoint + SIZE(frame%positions, KIND=int_8) + 12_int_8) + &
              INT(STORAGE_SIZE(0)/8, int_8)*SIZE(frame%images, KIND=int_8)
      IF (bytes > remaining_bytes) RETURN
      np = SIZE(points, 2)
      ALLOCATE (cache%left(np, nao), cache%right(np, nao), &
                cache%left_gradient(np, nao, 3), cache%right_gradient(np, nao, 3))
      CALL center_fields(frame, center, points, cache%left, cache%right, cache%left_gradient, cache%right_gradient)
      cache%center = center
      cache%cell = frame%cell
      cache%reduced_k = frame%reduced_k
      ALLOCATE (cache%positions, SOURCE=frame%positions)
      ALLOCATE (cache%images, SOURCE=frame%images)
      ALLOCATE (cache%coefficients, SOURCE=frame%coefficients)
      CALL MOVE_ALLOC(points, cache%points)
      remaining_bytes = remaining_bytes - bytes
   END SUBROUTINE qs_upf_inverse_density_cache_prepare

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
      REAL(dp)                                           :: lower(3), position(3), position_response(3), radius, upper(3)
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
      lower = MINVAL(points, DIM=2)
      upper = MAXVAL(points, DIM=2)
      DO center = 1, SIZE(frame%center_atoms)
         ikind = frame%center_kinds(center)
         CPASSERT(ikind >= 1 .AND. ikind <= SIZE(kinetic))
         IF (kinetic(ikind)%delta%nproj == 0) CYCLE
         radius = kinetic(ikind)%delta%radius
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
            IF (.NOT. qs_upf_inverse_augmentation_has_support(position, radius, lower, upper)) CYCLE
            IF (derive) position_response = atom_response(:, atom) + &
                                            MATMUL(cell_response, REAL(frame%images(:, image), dp))
            DO point = 1, np
               IF (SUM((points(:, point) - position)**2) >= radius**2) CYCLE
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
      REAL(dp)                                           :: lower(3), position(3), position_response(3), radius, upper(3)
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
      lower = MINVAL(points, DIM=2)
      upper = MAXVAL(points, DIM=2)
      DO center = 1, SIZE(frame%center_atoms)
         ikind = frame%center_kinds(center)
         CPASSERT(ikind >= 1 .AND. ikind <= SIZE(kinetic))
         IF (kinetic(ikind)%delta%nproj == 0) CYCLE
         radius = kinetic(ikind)%delta%radius
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
            IF (.NOT. qs_upf_inverse_augmentation_has_support(position, radius, lower, upper)) CYCLE
            IF (derive) position_response = atom_response(:, atom) + &
                                            MATMUL(cell_response, REAL(frame%images(:, image), dp))
            DO point = 1, np
               IF (SUM((points(:, point) - position)**2) >= radius**2) CYCLE
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
!> \param radial Optional signed harmonic density, without the radial measure.
!> \param response Optional fixed-matrix atom/cell/Bloch response of this complete center density.
!> \param atom_response Cartesian motion of all original atoms.
!> \param cell_response Physical lattice response.
!> \param reduced_k_response Reduced Bloch-vector response.
!> \param coefficients_response Complete dW from the original defining inverse quadrature.
!> \param atom_index Native Gaussian center, required when center is zero and no inverse correction is present.
!> \param cache Optional immutable point fields. Geometry responses bypass this cache.
!> \param metric_matrix Optional energy-weighted AO density for the simultaneous normalization response.
!> \param metric_response Radial monopole of the response contracted with metric_matrix.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_density_gapw(frame, augmentation, center, kind, pairs, matrix, radial, &
                                          response, atom_response, cell_response, reduced_k_response, &
                                          coefficients_response, atom_index, cache, metric_matrix, metric_response)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      TYPE(qs_upf_inverse_augmentation_type), INTENT(IN) :: augmentation
      INTEGER, INTENT(IN)                                :: center
      TYPE(qs_kind_type), INTENT(IN)                     :: kind
      TYPE(qs_upf_inverse_gaussian_pairs_type), &
         INTENT(IN)                                      :: pairs
      COMPLEX(dp), INTENT(IN)                            :: matrix(:, :)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: radial(:, :), response(:, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: atom_response(:, :), &
                                                            cell_response(3, 3), &
                                                            reduced_k_response(3)
      COMPLEX(dp), INTENT(IN), OPTIONAL                  :: coefficients_response(:, :)
      INTEGER, INTENT(IN), OPTIONAL                      :: atom_index
      TYPE(qs_upf_inverse_density_cache_type), &
         INTENT(IN), OPTIONAL                            :: cache
      COMPLEX(dp), INTENT(IN), OPTIONAL                  :: metric_matrix(:, :)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: metric_response(:, :)

      INTEGER                                            :: atom, n, nh, np
      LOGICAL                                            :: derive
      REAL(dp), ALLOCATABLE :: deh(:, :), des(:, :), dh(:, :), ds(:, :), dvalues(:), eh(:, :), &
         es(:, :), hard(:, :), metric_values(:), point_response(:, :), points(:, :), soft(:, :), values(:)
      TYPE(grid_atom_type), POINTER                      :: grid
      TYPE(harmonics_atom_type), POINTER                 :: harmonics
      TYPE(qs_upf_inverse_gapw_view_type), TARGET         :: view

      CALL qs_upf_inverse_gapw_grid(frame, center, kind, grid, harmonics, n, points, view, atom_index)
      IF (center > 0) CALL check_center(frame, augmentation, center, points)
      nh = harmonics%max_iso_not0
      np = SIZE(points, 2)
      derive = PRESENT(response)
      CPASSERT(PRESENT(radial) .OR. derive)
      CPASSERT(derive .EQV. PRESENT(atom_response))
      CPASSERT(derive .EQV. PRESENT(cell_response))
      CPASSERT(derive .EQV. PRESENT(reduced_k_response))
      CPASSERT(derive .EQV. PRESENT(coefficients_response))
      CPASSERT(PRESENT(metric_matrix) .EQV. PRESENT(metric_response))
      IF (PRESENT(metric_response)) THEN
         CPASSERT(derive)
         CPASSERT(ALL(SHAPE(metric_matrix) == SHAPE(matrix)))
         CPASSERT(ALL(SHAPE(metric_response) == [grid%nr, 1]))
         ALLOCATE (metric_values(np), SOURCE=0.0_dp)
      END IF
      IF (PRESENT(radial)) THEN
         CPASSERT(ALL(SHAPE(radial) == [grid%nr, nh]))
         ALLOCATE (hard(n, n), soft(n, n), eh(grid%nr, nh), es(grid%nr, nh), values(np))
         values = 0.0_dp
      END IF
      IF (derive) THEN
         CPASSERT(ALL(SHAPE(response) == [grid%nr, nh]))
         CPASSERT(ALL(SHAPE(atom_response) == SHAPE(frame%positions)))
         IF (PRESENT(atom_index)) THEN
            atom = atom_index
         ELSE
            atom = frame%center_atoms(center)
         END IF
         ALLOCATE (dh(n, n), ds(n, n), deh(grid%nr, nh), des(grid%nr, nh), &
                   dvalues(np), point_response(3, np))
         point_response(:, :) = SPREAD(atom_response(:, atom), 2, np)
         CALL qs_upf_inverse_gaussian_pairs_density(pairs, frame%reduced_k, matrix, hard, soft, &
                                                    dh, ds, atom_response, cell_response, reduced_k_response)
         ! The physical radial nodes are fixed. The derivative acts only on the coefficients.
         CALL qs_upf_inverse_gaussian_radial(kind, dh, deh)
         CALL qs_upf_inverse_gaussian_radial(kind, ds, des)
         dvalues = 0.0_dp
         IF (center > 0) THEN
            CALL qs_upf_inverse_density_center(frame, augmentation, center, points, matrix, values, dvalues, &
                                               point_response, atom_response, cell_response, reduced_k_response, &
                                               coefficients_response, metric_matrix=metric_matrix, &
                                               metric_response=metric_values)
         END IF
         CALL qs_upf_inverse_density_radial(grid, harmonics, dvalues, response, deh, des)
         IF (PRESENT(metric_response)) THEN
            CALL qs_upf_inverse_gaussian_pairs_density(pairs, frame%reduced_k, metric_matrix, &
                                                       hard_response=dh, soft_response=ds, atom_response=atom_response, &
                                                       cell_response=cell_response, reduced_k_response=reduced_k_response)
            CALL qs_upf_inverse_gaussian_radial(kind, dh, deh)
            CALL qs_upf_inverse_gaussian_radial(kind, ds, des)
            CALL qs_upf_inverse_density_radial(grid, harmonics, metric_values, metric_response, deh, des)
         END IF
      ELSE
         CALL qs_upf_inverse_gaussian_pairs_density(pairs, frame%reduced_k, matrix, hard, soft)
         IF (center > 0) CALL qs_upf_inverse_density_center(frame, augmentation, center, points, matrix, values, cache=cache)
      END IF
      IF (PRESENT(radial)) THEN
         CALL qs_upf_inverse_gaussian_radial(kind, hard, eh)
         CALL qs_upf_inverse_gaussian_radial(kind, soft, es)
         CALL qs_upf_inverse_density_radial(grid, harmonics, values, radial, eh, es)
      END IF
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
!> \param cache Optional immutable point fields. Geometry responses bypass this cache.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_density_gapw_integrate(frame, augmentation, center, kind, pairs, potential, ao_operator, &
                                                    atom_index, cache)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      TYPE(qs_upf_inverse_augmentation_type), INTENT(IN) :: augmentation
      INTEGER, INTENT(IN)                                :: center
      TYPE(qs_kind_type), INTENT(IN)                     :: kind
      TYPE(qs_upf_inverse_gaussian_pairs_type), &
         INTENT(IN)                                      :: pairs
      REAL(dp), INTENT(IN)                               :: potential(:, :)
      COMPLEX(dp), INTENT(OUT)                           :: ao_operator(:, :)
      INTEGER, INTENT(IN), OPTIONAL                      :: atom_index
      TYPE(qs_upf_inverse_density_cache_type), &
         INTENT(IN), OPTIONAL                            :: cache

      COMPLEX(dp), ALLOCATABLE                           :: gaussian(:, :)
      INTEGER                                            :: n, nh, np
      REAL(dp), ALLOCATABLE                              :: bh(:, :), bs(:, :), hard(:, :), &
                                                            points(:, :), soft(:, :), values(:), &
                                                            weights(:)
      TYPE(grid_atom_type), POINTER                      :: grid
      TYPE(harmonics_atom_type), POINTER                 :: harmonics
      TYPE(qs_upf_inverse_gapw_view_type), TARGET         :: view

      CALL qs_upf_inverse_gapw_grid(frame, center, kind, grid, harmonics, n, points, view, atom_index)
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
         CALL qs_upf_inverse_density_integrate(frame, augmentation, center, points, weights, values, ao_operator, cache=cache)
      END IF
      ao_operator = ao_operator + gaussian
      ! The signed center density is real. Its Gaussian and UPF adjoints
      ! can cancel large terms in the two independently assembled triangles.
      ao_operator = 0.5_dp*(ao_operator + CONJG(TRANSPOSE(ao_operator)))
   END SUBROUTINE qs_upf_inverse_density_gapw_integrate

! **************************************************************************************************
!> \brief Borrow the same native layout for both sides of the complete GAPW density map.
!> \param frame ...
!> \param center ...
!> \param kind ...
!> \param grid Native or averaged atomic quadrature.
!> \param harmonics Harmonic values at the same angular nodes.
!> \param n Number of compressed Gaussian primitives.
!> \param points Replicated Cartesian atom-grid points with angular rows fastest.
!> \param view Storage for a rotated grid and its harmonic values.
!> \param atom_index Explicit native Gaussian center. Center zero has no inverse or Q contribution.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_gapw_grid(frame, center, kind, grid, harmonics, n, points, view, atom_index)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      INTEGER, INTENT(IN)                                :: center
      TYPE(qs_kind_type), INTENT(IN)                     :: kind
      TYPE(grid_atom_type), POINTER                      :: grid
      TYPE(harmonics_atom_type), POINTER                 :: harmonics
      INTEGER, INTENT(OUT)                               :: n
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: points(:, :)
      TYPE(qs_upf_inverse_gapw_view_type), INTENT(OUT), TARGET :: view
      INTEGER, INTENT(IN), OPTIONAL                      :: atom_index

      INTEGER                                            :: atom, first, ia, ir, iso, na, nh, row
      REAL(dp), ALLOCATABLE                              :: directions(:, :)
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
      IF (ALLOCATED(frame%center_rotations)) THEN
         CALL qs_upf_atom_angular_average(grid, frame%center_rotations, directions, view%angular_weights)
      END IF
      IF (ALLOCATED(directions)) THEN
         na = SIZE(directions, 2)
         nh = harmonics%max_iso_not0
         ALLOCATE (view%weights(na, grid%nr), view%slm(na, nh), points(3, na*grid%nr))
         DO ir = 1, grid%nr
            first = 1 + (ir - 1)*na
            points(:, first:first + na - 1) = SPREAD(frame%positions(:, atom), 2, na) + grid%rad(ir)*directions
            view%weights(:, ir) = grid%wr(ir)*view%angular_weights
         END DO
         DO iso = 1, nh
            CALL y_lm(directions, view%slm(:, iso), indso(1, iso), indso(2, iso))
         END DO
         view%grid%nr = grid%nr
         view%grid%ng_sphere = na
         view%grid%rad => grid%rad
         view%grid%wr => grid%wr
         view%grid%wa => view%angular_weights
         view%grid%weight => view%weights
         view%harmonics%max_iso_not0 = nh
         view%harmonics%slm => view%slm
         grid => view%grid
         harmonics => view%harmonics
         RETURN
      END IF
      ALLOCATE (points(3, grid%nr*grid%ng_sphere))
      DO ir = 1, grid%nr
         DO ia = 1, grid%ng_sphere
            row = ia + (ir - 1)*grid%ng_sphere
            points(:, row) = frame%positions(:, atom) + grid%rad(ir)* &
                             [grid%sin_pol(ia)*grid%cos_azi(ia), grid%sin_pol(ia)*grid%sin_azi(ia), grid%cos_pol(ia)]
         END DO
      END DO
   END SUBROUTINE qs_upf_inverse_gapw_grid

! **************************************************************************************************
!> \brief Signed nR_a + Q_a from the full coupled inverse, without partition or image folding.
!>        nR_a = Re[d_a^* (s + chi) P], with d_a = -Delta_a W_a for the home image.
!>        Sum caller-supplied spin occupations and k weights outside this map exactly once.
!> \param frame Original inverse frame, including all centers and their coupled W.
!> \param augmentation Original charge kernels in the frame's magnetic-projector order.
!> \param center Selected home center, including a possible NC center with zero reconstruction/Q.
!> \param points Cartesian evaluation points, possibly an empty batch.
!> \param matrix Hermitian AO occupation density for one spin and k point.
!> \param density Optional signed home-center charge density, excluding Gaussian hard/soft and core.
!> \param response Optional joint density response at fixed matrix. All response inputs are required.
!> \param point_response Motion of the evaluation points.
!> \param atom_response Motion of every original atom.
!> \param cell_response Physical lattice response.
!> \param reduced_k_response Reduced Bloch-vector response.
!> \param coefficients_response Complete dW from the original defining inverse quadrature.
!> \param cache Optional immutable point fields. Geometry responses bypass this cache.
!> \param metric_matrix Optional second AO density on a grid moving rigidly with this center.
!> \param metric_response Geometry response contracted with the second AO density.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_density_center(frame, augmentation, center, points, matrix, density, response, &
                                            point_response, atom_response, cell_response, reduced_k_response, &
                                            coefficients_response, cache, metric_matrix, metric_response)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      TYPE(qs_upf_inverse_augmentation_type), INTENT(IN) :: augmentation
      INTEGER, INTENT(IN)                                :: center
      REAL(dp), INTENT(IN)                               :: points(:, :)
      COMPLEX(dp), INTENT(IN)                            :: matrix(:, :)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: density(:), response(:)
      REAL(dp), INTENT(IN), OPTIONAL                     :: point_response(:, :), &
                                                            atom_response(:, :), &
                                                            cell_response(3, 3), &
                                                            reduced_k_response(3)
      COMPLEX(dp), INTENT(IN), OPTIONAL                  :: coefficients_response(:, :)
      TYPE(qs_upf_inverse_density_cache_type), &
         INTENT(IN), OPTIONAL                            :: cache
      COMPLEX(dp), INTENT(IN), OPTIONAL                  :: metric_matrix(:, :)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: metric_response(:)

      COMPLEX(dp), ALLOCATABLE                           :: d(:, :), dd(:, :), de(:, :), e(:, :), &
                                                            projected(:, :), &
                                                            projected_response(:, :)
      INTEGER                                            :: atom, first, kind, last, n, np, nproj
      LOGICAL                                            :: derive, moving_q
      REAL(dp), ALLOCATABLE                              :: displacements(:, :), dq(:), q(:), &
                                                            qgradient(:, :)

      CALL check_center(frame, augmentation, center, points)
      np = SIZE(points, 2)
      n = SIZE(frame%coefficients, 2)
      CPASSERT(ALL(SHAPE(matrix) == [n, n]))
      derive = PRESENT(response)
      CPASSERT(PRESENT(density) .OR. derive)
      CPASSERT(derive .EQV. PRESENT(point_response))
      CPASSERT(derive .EQV. PRESENT(atom_response))
      CPASSERT(derive .EQV. PRESENT(cell_response))
      CPASSERT(derive .EQV. PRESENT(reduced_k_response))
      CPASSERT(derive .EQV. PRESENT(coefficients_response))
      CPASSERT(PRESENT(metric_matrix) .EQV. PRESENT(metric_response))
      IF (PRESENT(metric_response)) THEN
         CPASSERT(derive)
         CPASSERT(ALL(SHAPE(metric_matrix) == SHAPE(matrix)))
         CPASSERT(SIZE(metric_response) == np)
         metric_response = 0.0_dp
      END IF
      IF (PRESENT(density)) THEN
         CPASSERT(SIZE(density) == np)
         density = 0.0_dp
      END IF
      IF (derive) THEN
         CPASSERT(SIZE(response) == np)
         CPASSERT(ALL(SHAPE(point_response) == SHAPE(points)))
         CPASSERT(ALL(SHAPE(atom_response) == SHAPE(frame%positions)))
         CPASSERT(ALL(SHAPE(coefficients_response) == SHAPE(frame%coefficients)))
         response = 0.0_dp
      END IF
      ! DIRECT has Delta=0. Its reconstruction density and every geometry response vanish.
      ! The original Q density and its projector and geometry responses remain below.
      IF (.NOT. frame%direct) THEN
         ALLOCATE (d(np, n), e(np, n))
         IF (derive) THEN
            ALLOCATE (dd(np, n), de(np, n), dq(np))
            CALL center_fields(frame, center, points, d, e, dd=dd, de=de, &
                               point_response=point_response, atom_response=atom_response, cell_response=cell_response, &
                               reduced_k_response=reduced_k_response, coefficients_response=coefficients_response)
            CALL contract_fields(matrix, response)
            IF (PRESENT(metric_response)) CALL contract_fields(metric_matrix, metric_response)
         ELSE
            CALL center_fields(frame, center, points, d, e, cache=cache)
         END IF
         IF (PRESENT(density)) CALL upf_inverse_cross_density(d, e, ao_density=matrix, density=density)
      END IF
      atom = frame%center_atoms(center)
      kind = frame%center_kinds(center)
      IF (augmentation%kinds(kind)%nproj == 0) RETURN
      first = frame%projector_first(center)
      last = frame%projector_last(center)
      nproj = last - first + 1
      ALLOCATE (q(np), displacements(3, np))
      displacements(:, :) = points - SPREAD(frame%positions(:, atom), 2, np)
      moving_q = .FALSE.
      IF (derive) moving_q = ANY(point_response /= SPREAD(atom_response(:, atom), 2, np))
      IF (PRESENT(metric_response)) THEN
         CPASSERT(.NOT. moving_q)
      END IF
      IF (PRESENT(density) .OR. moving_q) THEN
         ALLOCATE (projected(nproj, nproj))
         CALL upf_inverse_projector_density(frame%coefficients(first:last, :), matrix, projected)
         q = 0.0_dp
         IF (moving_q) ALLOCATE (qgradient(3, np), SOURCE=0.0_dp)
         CALL uspp_augmentation_collocate_atom(augmentation%kinds(kind), displacements, &
                                               REAL(projected, dp), q, qgradient)
         IF (PRESENT(density)) density = density + q
         IF (moving_q) response = response + SUM(qgradient*(point_response - SPREAD(atom_response(:, atom), 2, np)), DIM=1)
      END IF
      IF (derive) THEN
         ALLOCATE (projected_response(nproj, nproj))
         CALL contract_projectors(matrix, response)
         IF (PRESENT(metric_response)) CALL contract_projectors(metric_matrix, metric_response)
      END IF

   CONTAINS

!> \brief Contract the shared orbital fields and their derivatives with one AO matrix.
!> \param ao_density Hermitian density in the original Gaussian order.
!> \param field_response Pointwise reconstruction response.
      SUBROUTINE contract_fields(ao_density, field_response)
         COMPLEX(dp), INTENT(IN)                         :: ao_density(:, :)
         REAL(dp), INTENT(OUT)                           :: field_response(:)

         CALL upf_inverse_cross_density(dd, e, ao_density=ao_density, density=field_response)
         CALL upf_inverse_cross_density(d, de, ao_density=ao_density, density=dq)
         field_response = field_response + dq
      END SUBROUTINE contract_fields

!> \brief Add the derivative of the projected density at fixed radial grid coordinates.
!> \param ao_density Hermitian density in the original Gaussian order.
!> \param field_response Accumulated pointwise charge response.
      SUBROUTINE contract_projectors(ao_density, field_response)
         COMPLEX(dp), INTENT(IN)                         :: ao_density(:, :)
         REAL(dp), INTENT(INOUT)                         :: field_response(:)

         projected_response(:, :) = MATMUL(coefficients_response(first:last, :), &
                                           MATMUL(ao_density, CONJG(TRANSPOSE(frame%coefficients(first:last, :)))))
         projected_response(:, :) = projected_response + CONJG(TRANSPOSE(projected_response))
         q = 0.0_dp
         CALL uspp_augmentation_collocate_atom(augmentation%kinds(kind), displacements, &
                                               REAL(projected_response, dp), q)
         field_response = field_response + q
      END SUBROUTINE contract_projectors
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
!> \param cache Optional immutable point fields. Geometry responses bypass this cache.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_density_integrate(frame, augmentation, center, points, weights, potential, ao_operator, comm, cache)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      TYPE(qs_upf_inverse_augmentation_type), INTENT(IN) :: augmentation
      INTEGER, INTENT(IN)                                :: center
      REAL(dp), INTENT(IN)                               :: points(:, :), weights(:), potential(:)
      COMPLEX(dp), INTENT(OUT)                           :: ao_operator(:, :)
      TYPE(mp_comm_type), INTENT(IN), OPTIONAL           :: comm
      TYPE(qs_upf_inverse_density_cache_type), &
         INTENT(IN), OPTIONAL                            :: cache

      COMPLEX(dp), ALLOCATABLE                           :: d(:, :), e(:, :), gd(:, :, :), &
                                                            ge(:, :, :), qoperator(:, :)
      INTEGER                                            :: atom, first, kind, last, n, np
      REAL(dp), ALLOCATABLE                              :: q(:, :), zero(:), zeros(:, :)

      CALL check_center(frame, augmentation, center, points)
      np = SIZE(points, 2)
      n = SIZE(frame%coefficients, 2)
      CPASSERT(SIZE(weights) == np .AND. SIZE(potential) == np)
      CPASSERT(ALL(SHAPE(ao_operator) == [n, n]))
      ALLOCATE (zero(np), zeros(np, 3))
      zero = 0.0_dp
      zeros = 0.0_dp
      ao_operator = CMPLX(0.0_dp, 0.0_dp, dp)
      IF (.NOT. frame%direct) THEN
         ALLOCATE (d(np, n), e(np, n), gd(np, n, 3), ge(np, n, 3))
         CALL center_fields(frame, center, points, d, e, gd, ge, cache=cache)
         CALL upf_inverse_semilocal_pair(d, e, gd, ge, weights, potential, zeros, zero, ao_operator)
         ao_operator = 0.5_dp*(ao_operator + CONJG(TRANSPOSE(ao_operator)))
      END IF
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
!> \brief Independent-basis home-center charge residual and its complete geometry response.
!>        For each basis d_a is its home reconstruction difference and e=s+chi.
!>        The rectangular operator is (d_X^dagger e_Y+e_X^dagger d_Y)/2+W_X^dagger Q_a W_Y.
!>        Summing home centers gives the inverse-minus-Gaussian charge metric, including original Q.
!> \param frame Row basis with its original defining inverse and image domain.
!> \param right_frame Independent column basis for the same original atoms and UPF datasets.
!> \param augmentation Original home-center charge kernels shared by both bases.
!> \param center Same original center index and magnetic-projector order in both frames.
!> \param points Common Cartesian evaluation rows, without image folding or a new inverse fit.
!> \param weights Common point measure, applied once. Replicated atom grids omit comm.
!> \param potential Real point covector, without a measure or mesh weight.
!> \param ao_operator Rectangular row/column operator, with no Hermitian completion or mesh weight.
!> \param ao_response Optional joint response, requiring every remaining response argument.
!> \param point_response Motion of the common evaluation rows.
!> \param atom_response Motion of all original atoms.
!> \param cell_response Physical lattice response.
!> \param reduced_k_response Reduced Bloch-vector response.
!> \param coefficients_response Complete row dW from its original defining quadrature.
!> \param right_coefficients_response Complete independent column dW.
!> \param weights_response Response of the same point measure.
!> \param potential_response Response of the unweighted point covector.
!> \param comm Optional communicator distributing the common point rows.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_density_pair_integrate(frame, right_frame, augmentation, center, points, weights, &
                                                    potential, ao_operator, ao_response, point_response, atom_response, &
                                                    cell_response, reduced_k_response, coefficients_response, &
                                                    right_coefficients_response, weights_response, potential_response, comm)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame, right_frame
      TYPE(qs_upf_inverse_augmentation_type), INTENT(IN) :: augmentation
      INTEGER, INTENT(IN)                                :: center
      REAL(dp), INTENT(IN)                               :: points(:, :), weights(:), potential(:)
      COMPLEX(dp), INTENT(OUT)                           :: ao_operator(:, :)
      COMPLEX(dp), INTENT(OUT), OPTIONAL                 :: ao_response(:, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: point_response(:, :), &
                                                            atom_response(:, :), &
                                                            cell_response(3, 3), &
                                                            reduced_k_response(3)
      COMPLEX(dp), INTENT(IN), OPTIONAL :: coefficients_response(:, :), &
         right_coefficients_response(:, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: weights_response(:), &
                                                            potential_response(:)
      TYPE(mp_comm_type), INTENT(IN), OPTIONAL           :: comm

      COMPLEX(dp), ALLOCATABLE :: ddx(:, :), ddy(:, :), dex(:, :), dey(:, :), dx(:, :), dy(:, :), &
         ex(:, :), ey(:, :), partial(:, :)
      INTEGER                                            :: atom, block, first, hi, kind, last, lo, &
                                                            m, nl, np, nr
      LOGICAL                                            :: derive
      REAL(dp), ALLOCATABLE :: dq(:, :), measure(:), measure_response(:), one(:), q(:, :), &
         qpart(:, :), translated(:, :), vg(:, :), zero(:), zeros(:, :)

      CALL check_center(frame, augmentation, center, points)
      CALL check_center(right_frame, augmentation, center, points)
      CPASSERT(ALL(SHAPE(frame%positions) == SHAPE(right_frame%positions)))
      CPASSERT(ALL(frame%positions == right_frame%positions))
      CPASSERT(ALL(frame%cell == right_frame%cell))
      CPASSERT(ALL(frame%reduced_k == right_frame%reduced_k))
      CPASSERT(SIZE(frame%center_atoms) == SIZE(right_frame%center_atoms))
      CPASSERT(ALL(frame%center_atoms == right_frame%center_atoms))
      CPASSERT(ALL(frame%center_kinds == right_frame%center_kinds))
      CPASSERT(ALL(frame%projector_first == right_frame%projector_first))
      CPASSERT(ALL(frame%projector_last == right_frame%projector_last))
      nl = SIZE(frame%coefficients, 2)
      nr = SIZE(right_frame%coefficients, 2)
      np = SIZE(points, 2)
      CPASSERT(SIZE(weights) == np .AND. SIZE(potential) == np)
      CPASSERT(ALL(SHAPE(ao_operator) == [nl, nr]))
      derive = PRESENT(ao_response)
      CPASSERT(derive .EQV. PRESENT(point_response))
      CPASSERT(derive .EQV. PRESENT(atom_response))
      CPASSERT(derive .EQV. PRESENT(cell_response))
      CPASSERT(derive .EQV. PRESENT(reduced_k_response))
      CPASSERT(derive .EQV. PRESENT(coefficients_response))
      CPASSERT(derive .EQV. PRESENT(right_coefficients_response))
      CPASSERT(derive .EQV. PRESENT(weights_response))
      CPASSERT(derive .EQV. PRESENT(potential_response))
      ao_operator = CMPLX(0.0_dp, 0.0_dp, dp)
      IF (derive) THEN
         CPASSERT(ALL(SHAPE(ao_response) == [nl, nr]))
         CPASSERT(ALL(SHAPE(point_response) == SHAPE(points)))
         CPASSERT(ALL(SHAPE(atom_response) == SHAPE(frame%positions)))
         CPASSERT(ALL(SHAPE(coefficients_response) == SHAPE(frame%coefficients)))
         IF (.NOT. ALL(SHAPE(right_coefficients_response) == SHAPE(right_frame%coefficients))) THEN
            CPABORT("Invalid column inverse coefficient response")
         END IF
         CPASSERT(SIZE(weights_response) == np .AND. SIZE(potential_response) == np)
         ao_response = CMPLX(0.0_dp, 0.0_dp, dp)
      END IF
      block = MIN(qs_upf_inverse_block_size(frame), qs_upf_inverse_block_size(right_frame))
      ALLOCATE (partial(nl, nr), measure(block))
      IF (.NOT. frame%direct .OR. .NOT. right_frame%direct) THEN
         ALLOCATE (dx(block, nl), ex(block, nl), dy(block, nr), ey(block, nr))
         IF (derive) ALLOCATE (ddx(block, nl), dex(block, nl), ddy(block, nr), dey(block, nr), &
                               measure_response(block))
         DO lo = 1, np, block
            hi = MIN(np, lo + block - 1)
            m = hi - lo + 1
            measure(:m) = weights(lo:hi)*potential(lo:hi)
            IF (derive) THEN
               measure_response(:m) = weights_response(lo:hi)*potential(lo:hi) + &
                                      weights(lo:hi)*potential_response(lo:hi)
               CALL pair_fields(frame, center, points(:, lo:hi), dx(:m, :), ex(:m, :), &
                                ddx(:m, :), dex(:m, :), point_response(:, lo:hi), atom_response, cell_response, &
                                reduced_k_response, coefficients_response)
               CALL pair_fields(right_frame, center, points(:, lo:hi), dy(:m, :), ey(:m, :), &
                                ddy(:m, :), dey(:m, :), point_response(:, lo:hi), atom_response, cell_response, &
                                reduced_k_response, right_coefficients_response)
               CALL upf_inverse_pair_response(dx(:m, :), ey(:m, :), measure(:m), ddx(:m, :), dey(:m, :), &
                                              measure_response(:m), partial)
               ao_response(:, :) = ao_response + 0.5_dp*partial
               CALL upf_inverse_pair_response(ex(:m, :), dy(:m, :), measure(:m), dex(:m, :), ddy(:m, :), &
                                              measure_response(:m), partial)
               ao_response(:, :) = ao_response + 0.5_dp*partial
            ELSE
               CALL pair_fields(frame, center, points(:, lo:hi), dx(:m, :), ex(:m, :))
               CALL pair_fields(right_frame, center, points(:, lo:hi), dy(:m, :), ey(:m, :))
            END IF
            CALL upf_inverse_pair(dx(:m, :), ey(:m, :), measure(:m), partial)
            ao_operator(:, :) = ao_operator + 0.5_dp*partial
            CALL upf_inverse_pair(ex(:m, :), dy(:m, :), measure(:m), partial)
            ao_operator(:, :) = ao_operator + 0.5_dp*partial
         END DO
      END IF
      atom = frame%center_atoms(center)
      kind = frame%center_kinds(center)
      IF (augmentation%kinds(kind)%nproj /= 0) THEN
         first = frame%projector_first(center)
         last = frame%projector_last(center)
         ALLOCATE (q(last - first + 1, last - first + 1), translated(3, np), zeros(3, np))
         translated(:, :) = points - SPREAD(frame%positions(:, atom), 2, np)
         zeros = 0.0_dp
         CALL uspp_augmentation_integrate_atom(augmentation%kinds(kind), translated, weights, potential, zeros, q)
         ao_operator(:, :) = ao_operator + MATMUL(CONJG(TRANSPOSE(frame%coefficients(first:last, :))), &
                                                  MATMUL(CMPLX(q, 0.0_dp, dp), right_frame%coefficients(first:last, :)))
         IF (derive) THEN
            ALLOCATE (dq, MOLD=q)
            ALLOCATE (qpart, MOLD=q)
            ALLOCATE (one(np), zero(np), vg(3, np))
            one = 1.0_dp
            zero = 0.0_dp
            vg(:, :) = SPREAD(potential, 1, 3)*(point_response - SPREAD(atom_response(:, atom), 2, np))
            CALL uspp_augmentation_integrate_atom(augmentation%kinds(kind), translated, weights, zero, vg, dq)
            CALL uspp_augmentation_integrate_atom(augmentation%kinds(kind), translated, &
                                                  weights_response*potential + weights*potential_response, one, zeros, qpart)
            dq(:, :) = dq + qpart
            ao_response(:, :) = ao_response + &
                                MATMUL(CONJG(TRANSPOSE(coefficients_response(first:last, :))), &
                                       MATMUL(CMPLX(q, 0.0_dp, dp), right_frame%coefficients(first:last, :))) + &
                                MATMUL(CONJG(TRANSPOSE(frame%coefficients(first:last, :))), &
                                       MATMUL(CMPLX(q, 0.0_dp, dp), right_coefficients_response(first:last, :))) + &
                                MATMUL(CONJG(TRANSPOSE(frame%coefficients(first:last, :))), &
                                       MATMUL(CMPLX(dq, 0.0_dp, dp), right_frame%coefficients(first:last, :)))
         END IF
      END IF
      IF (PRESENT(comm)) THEN
         CALL comm%sum(ao_operator)
         IF (derive) CALL comm%sum(ao_response)
      END IF
   END SUBROUTINE qs_upf_inverse_density_pair_integrate

! **************************************************************************************************
!> \brief Home difference and s+chi for independent bases, including an unchanged DIRECT basis.
!> \param frame ...
!> \param center ...
!> \param points ...
!> \param d ...
!> \param e ...
!> \param dd ...
!> \param de ...
!> \param point_response ...
!> \param atom_response ...
!> \param cell_response ...
!> \param reduced_k_response ...
!> \param coefficients_response ...
! **************************************************************************************************
   SUBROUTINE pair_fields(frame, center, points, d, e, dd, de, point_response, atom_response, &
                          cell_response, reduced_k_response, coefficients_response)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      INTEGER, INTENT(IN)                                :: center
      REAL(dp), INTENT(IN)                               :: points(:, :)
      COMPLEX(dp), INTENT(OUT)                           :: d(:, :), e(:, :)
      COMPLEX(dp), INTENT(OUT), OPTIONAL                 :: dd(:, :), de(:, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: point_response(:, :), &
                                                            atom_response(:, :), &
                                                            cell_response(3, 3), &
                                                            reduced_k_response(3)
      COMPLEX(dp), INTENT(IN), OPTIONAL                  :: coefficients_response(:, :)

      IF (frame%direct) THEN
         d = CMPLX(0.0_dp, 0.0_dp, dp)
         IF (PRESENT(dd)) dd = CMPLX(0.0_dp, 0.0_dp, dp)
         CALL qs_upf_inverse_evaluate(frame, points, e, soft_response=de, point_response=point_response, &
                                      atom_response=atom_response, cell_response=cell_response, &
                                      reduced_k_response=reduced_k_response, coefficients_response=coefficients_response)
         e = 2.0_dp*e
         IF (PRESENT(de)) de = 2.0_dp*de
      ELSE
         CALL center_fields(frame, center, points, d, e, dd=dd, de=de, point_response=point_response, &
                            atom_response=atom_response, cell_response=cell_response, &
                            reduced_k_response=reduced_k_response, coefficients_response=coefficients_response)
      END IF
   END SUBROUTINE pair_fields

! **************************************************************************************************
!> \brief Native angular projection of nR + Q and optional Gaussian hard-minus-soft density.
!>        Use the same finite Lebedev map for every leg. Its Gram matrix is not assumed to be unity.
!>        The native harmonic order is l=0,1,... and m=-l,...,l, as in the Coulomb maps.
!> \param grid Borrowed native radial/angular quadrature, with complete replicated angular rows.
!> \param harmonics Borrowed native real harmonics on that same angular quadrature.
!> \param values Point density with angular rows fastest. No core or partition is added here.
!> \param radial Leading real harmonic density coefficients, without the radial measure.
!> \param gaussian_hard Optional native Gaussian hard harmonic coefficients.
!> \param gaussian_soft Matching native Gaussian soft harmonic coefficients.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_density_radial(grid, harmonics, values, radial, gaussian_hard, gaussian_soft)
      TYPE(grid_atom_type), INTENT(IN)                   :: grid
      TYPE(harmonics_atom_type), INTENT(IN)              :: harmonics
      REAL(dp), INTENT(IN)                               :: values(:)
      REAL(dp), INTENT(OUT)                              :: radial(:, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: gaussian_hard(:, :), gaussian_soft(:, :)

      INTEGER                                            :: first, ir, last, na, nh, nh_gaussian
      REAL(dp), ALLOCATABLE                              :: angular(:)

      na = grid%ng_sphere
      nh = SIZE(radial, 2)
      CALL check_radial(grid, harmonics, values, radial)
      CPASSERT(PRESENT(gaussian_hard) .EQV. PRESENT(gaussian_soft))
      IF (PRESENT(gaussian_hard)) THEN
         nh_gaussian = SIZE(gaussian_hard, 2)
         CPASSERT(nh_gaussian <= SIZE(harmonics%slm, 2))
         CPASSERT(SIZE(gaussian_hard, 1) == grid%nr)
         CPASSERT(ALL(SHAPE(gaussian_hard) == SHAPE(gaussian_soft)))
      END IF
      ALLOCATE (angular(na))
      DO ir = 1, grid%nr
         first = (ir - 1)*na + 1
         last = ir*na
         angular(:) = values(first:last)
         IF (PRESENT(gaussian_hard)) THEN
            ! All input harmonics contribute to the finite Gram map, even for monopole output.
            angular(:) = angular + MATMUL(harmonics%slm(:, 1:nh_gaussian), gaussian_hard(ir, :) - gaussian_soft(ir, :))
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
!> \param gd Optional spatial gradient of the home-center difference.
!> \param ge Matching spatial gradient of s+chi.
!> \param dd ...
!> \param de ...
!> \param point_response ...
!> \param atom_response ...
!> \param cell_response ...
!> \param reduced_k_response ...
!> \param coefficients_response ...
!> \param cache Optional immutable point fields. Geometry responses bypass this cache.
! **************************************************************************************************
   SUBROUTINE center_fields(frame, center, points, d, e, gd, ge, dd, de, &
                            point_response, atom_response, cell_response, reduced_k_response, coefficients_response, cache)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      INTEGER, INTENT(IN)                                :: center
      REAL(dp), INTENT(IN)                               :: points(:, :)
      COMPLEX(dp), INTENT(OUT)                           :: d(:, :), e(:, :)
      COMPLEX(dp), INTENT(OUT), OPTIONAL                 :: dd(:, :), de(:, :), gd(:, :, :), ge(:, :, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: point_response(:, :), &
                                                            atom_response(:, :), &
                                                            cell_response(3, 3), &
                                                            reduced_k_response(3)
      COMPLEX(dp), INTENT(IN), OPTIONAL                  :: coefficients_response(:, :)
      TYPE(qs_upf_inverse_density_cache_type), &
         INTENT(IN), OPTIONAL                            :: cache

      COMPLEX(dp), ALLOCATABLE                           :: chi(:, :), dchi(:, :), gc(:, :, :)
      INTEGER                                            :: n, np

      n = SIZE(frame%coefficients, 2)
      np = SIZE(points, 2)
      CPASSERT(PRESENT(gd) .EQV. PRESENT(ge))
      IF (PRESENT(cache)) THEN
         IF (ALLOCATED(cache%left) .AND. &
             .NOT. (PRESENT(dd) .OR. PRESENT(de) .OR. PRESENT(point_response) .OR. &
                    PRESENT(atom_response) .OR. PRESENT(cell_response) .OR. &
                    PRESENT(reduced_k_response) .OR. PRESENT(coefficients_response))) THEN
            CPASSERT(cache%center == center)
            CPASSERT(ALL(SHAPE(cache%points) == SHAPE(points)))
            CPASSERT(ALL(SHAPE(cache%positions) == SHAPE(frame%positions)))
            CPASSERT(ALL(SHAPE(cache%images) == SHAPE(frame%images)))
            CPASSERT(ALL(SHAPE(cache%coefficients) == SHAPE(frame%coefficients)))
            CPASSERT(ALL(cache%points == points))
            CPASSERT(ALL(cache%positions == frame%positions))
            CPASSERT(ALL(cache%cell == frame%cell))
            CPASSERT(ALL(cache%reduced_k == frame%reduced_k))
            CPASSERT(ALL(cache%images == frame%images))
            CPASSERT(ALL(cache%coefficients == frame%coefficients))
            d = cache%left
            e = cache%right
            IF (PRESENT(gd)) THEN
               gd = cache%left_gradient
               ge = cache%right_gradient
            END IF
            RETURN
         END IF
      END IF
      ALLOCATE (chi(np, n))
      IF (PRESENT(ge)) ALLOCATE (gc(np, n, 3))
      IF (PRESENT(de)) ALLOCATE (dchi(np, n))
      CALL qs_upf_inverse_center_fields(frame, center, points, d, gd, dd, point_response=point_response, &
                                        atom_response=atom_response, cell_response=cell_response, &
                                        reduced_k_response=reduced_k_response, coefficients_response=coefficients_response, &
                                        image_shift=[0, 0, 0])
      ! Reuse the original Gaussian fields and their joint response before inverse reconstruction.
      CALL qs_upf_inverse_evaluate(frame, points, e, ge, de, point_response=point_response, &
                                   atom_response=atom_response, cell_response=cell_response, &
                                   reduced_k_response=reduced_k_response, coefficients_response=coefficients_response, &
                                   gaussian=chi, gaussian_gradient=gc, gaussian_response=dchi)
      e = e + chi
      IF (PRESENT(ge)) ge = ge + gc
      IF (PRESENT(de)) de = de + dchi
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
