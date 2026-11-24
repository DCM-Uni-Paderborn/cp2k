!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Inverse USPP fields on native PW rows and real/imaginary Quickstep AO matrices.
!>        Point rows and projector contractions use the native PW communicator.
!>        AO columns and the returned dense operators remain replicated.
! **************************************************************************************************
MODULE qs_upf_inverse_grid
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_create, dbcsr_desymmetrize, dbcsr_finalize, dbcsr_get_info, &
        dbcsr_get_stored_coordinates, dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, &
        dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_put_block, &
        dbcsr_release, dbcsr_type, dbcsr_type_antisymmetric, dbcsr_type_no_symmetry, &
        dbcsr_type_symmetric
   USE gaussian_basis_fields,           ONLY: gaussian_basis_fields_periodic
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE mathlib,                         ONLY: invmat
   USE message_passing,                 ONLY: mp_comm_congruent,&
                                              mp_comm_ident,&
                                              mp_comm_self,&
                                              mp_comm_type
   USE pw_grid_types,                   ONLY: PW_MODE_DISTRIBUTED,&
                                              PW_MODE_LOCAL,&
                                              pw_grid_type
   USE pw_types,                        ONLY: pw_r3d_rs_type
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_upf_inverse,                  ONLY: max_point_field_bytes,&
                                              qs_upf_inverse_block_size,&
                                              qs_upf_inverse_evaluate,&
                                              qs_upf_inverse_from_ks,&
                                              qs_upf_inverse_release,&
                                              qs_upf_inverse_response,&
                                              qs_upf_inverse_type
   USE qs_upf_inverse_augmentation,     ONLY: qs_upf_inverse_augmentation_collocate,&
                                              qs_upf_inverse_augmentation_fourier_collocate,&
                                              qs_upf_inverse_augmentation_fourier_integrate,&
                                              qs_upf_inverse_augmentation_fourier_response,&
                                              qs_upf_inverse_augmentation_integrate,&
                                              qs_upf_inverse_augmentation_type
   USE qs_upf_inverse_density,          ONLY: qs_upf_inverse_kinetic_collocate,&
                                              qs_upf_inverse_kinetic_operator
   USE qs_upf_inverse_types,            ONLY: qs_upf_inverse_grid_type
   USE upf_inverse,                     ONLY: upf_inverse_kinetic_type
   USE upf_inverse_quadrature,          ONLY: &
        upf_inverse_pair, upf_inverse_pair_response, upf_inverse_projector_density, upf_inverse_projector_pair, &
        upf_inverse_semilocal_pair, upf_inverse_semilocal_response, upf_inverse_smooth_density
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   INTEGER, PARAMETER :: subdivisions = 2, interpolation_order = 8

   PUBLIC :: qs_upf_inverse_grid_type, qs_upf_inverse_grid_prepare, qs_upf_inverse_grid_release
   PUBLIC :: qs_upf_inverse_grid_collocate, qs_upf_inverse_grid_integrate
   PUBLIC :: qs_upf_inverse_grid_metric
   PUBLIC :: qs_upf_inverse_grid_collocate_dense, qs_upf_inverse_grid_integrate_dense, qs_upf_inverse_grid_metric_dense
   PUBLIC :: qs_upf_inverse_matrix_gather, qs_upf_inverse_matrix_scatter
   PUBLIC :: qs_upf_inverse_grid_response, qs_upf_inverse_grid_charge_type
   PUBLIC :: qs_upf_inverse_grid_gaussian_operator

   TYPE qs_upf_inverse_grid_charge_type
      COMPLEX(dp), ALLOCATABLE :: operator(:, :), derivative(:, :)
   END TYPE qs_upf_inverse_grid_charge_type

CONTAINS

! **************************************************************************************************
!> \brief Prepare original inverse fields with the native grid coordinates and volume element.
!>        A local PW grid is complete on each rank and uses SELF for point integrals.
!>        A distributed PW grid uses its own communicator, with each point counted once.
!> \param ks_env Prepared native KS environment; its ORB bases remain borrowed.
!> \param grid Native PW grid, including its auxiliary-cell geometry when applicable.
!> \param images Explicit images of the physical KS cell, held fixed for responses.
!> \param reduced_k Reduced Bloch vector of the physical KS cell.
!> \param frame New owned adapter; neither the grid nor the KS environment is modified.
!> \param info Inverse-factorization status.
!> \param atom_local Optional one-center reconstruction with native Gaussian projector integrals.
!> \param direct Use original projectors and soft orbitals without reconstruction.
!> \param store_fields Override the bounded default cache for dense reference or streamed evaluation.
!> \param basis_type Native basis role, ORB by default.
!> \param center_rotations Fixed reference group for the full GAPW atom integrals.
!> \param refine Use moving atomic quadrature for the compact orbital terms at Gamma.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_grid_prepare(ks_env, grid, images, reduced_k, frame, info, atom_local, direct, store_fields, &
                                          basis_type, center_rotations, refine)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(pw_grid_type), INTENT(IN)                     :: grid
      INTEGER, INTENT(IN)                                :: images(:, :)
      REAL(dp), INTENT(IN)                               :: reduced_k(3)
      TYPE(qs_upf_inverse_grid_type), INTENT(OUT)        :: frame
      INTEGER, INTENT(OUT)                               :: info
      LOGICAL, INTENT(IN), OPTIONAL                      :: atom_local, direct, store_fields
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL             :: basis_type
      REAL(dp), INTENT(IN), OPTIONAL                     :: center_rotations(:, :, :)

      LOGICAL, INTENT(IN), OPTIONAL                      :: refine

      INTEGER                                            :: ix, iy, iz, n, point
      REAL(dp), ALLOCATABLE                              :: points(:, :), weights(:)
      TYPE(mp_comm_type)                                 :: comm

      frame%bounds = grid%bounds
      frame%bounds_local = grid%bounds_local
      frame%npts = grid%npts
      frame%dh = grid%dh
      frame%dvol = grid%dvol
      frame%mode = grid%para%mode
      CPASSERT(frame%dvol > 0.0_dp .AND. ALL(frame%npts > 0))
      n = PRODUCT(MAX(0, frame%bounds_local(2, :) - frame%bounds_local(1, :) + 1))
      ! Reciprocal blocks need not have the size of the local real-space rectangle.
      CPASSERT(n == PRODUCT(MAX(0, grid%npts_local)))
      IF (frame%mode == PW_MODE_LOCAL) THEN
         CPASSERT(ALL(frame%bounds_local == frame%bounds))
         comm = mp_comm_self
      ELSE
         CPASSERT(frame%mode == PW_MODE_DISTRIBUTED)
         comm = grid%para%group%mp_comm_type
      END IF
      ALLOCATE (points(3, n), weights(n))
      weights = frame%dvol
      point = 0
      DO iz = frame%bounds_local(1, 3), frame%bounds_local(2, 3)
         DO iy = frame%bounds_local(1, 2), frame%bounds_local(2, 2)
            DO ix = frame%bounds_local(1, 1), frame%bounds_local(2, 1)
               point = point + 1
               points(:, point) = MATMUL(frame%dh, REAL([ix, iy, iz] - frame%bounds(1, :), dp))
            END DO
         END DO
      END DO
      CALL qs_upf_inverse_from_ks(ks_env, points, weights, images, reduced_k, comm, frame%ao, info, atom_local, direct, &
                                  store_fields, basis_type, center_rotations)
      IF (info /= 0) RETURN
      IF (PRESENT(refine)) THEN
         IF (refine .AND. .NOT. frame%ao%direct .AND. ALL(frame%ao%periodicity == 1)) THEN
            CALL reconstruction_quadrature(frame, points, weights)
            CALL MOVE_ALLOC(points, frame%orbital_points)
            CALL MOVE_ALLOC(weights, frame%orbital_weights)
            ALLOCATE (frame%row_center(SIZE(frame%orbital_weights)))
            frame%row_center(:n) = 0
            point = SIZE(frame%patch_center)
            frame%row_center(n + 1:n + point) = frame%patch_center
            frame%row_center(n + point + 1:) = -frame%patch_center
         END IF
      END IF
      IF (.NOT. PRESENT(store_fields)) CALL cache_gamma_fields(frame)
   END SUBROUTINE qs_upf_inverse_grid_prepare

!> \brief Build atom-centered midpoint samples and their local GPW interpolation stencils.
!> \param frame Prepared inverse frame that owns the interpolation stencils.
!> \param points Native grid followed by two copies of each atom-centered sample.
!> \param weights Signed integration weights for the orbital products.
   SUBROUTINE reconstruction_quadrature(frame, points, weights)
      TYPE(qs_upf_inverse_grid_type), INTENT(INOUT)      :: frame
      REAL(dp), ALLOCATABLE, INTENT(INOUT)               :: points(:, :), weights(:)

      INTEGER :: axis, base(3), box(3), center, count, extent(3), &
         grid_indices(interpolation_order, 3), i, info, j, kind, n, nc, owner(3), pass, q(3), &
         step(3), x, y, z
      REAL(dp)                                           :: coord(3), fraction, inverse_dh(3, 3), &
                                                            offset(3), origin(3), radius, width(3)
      REAL(dp), ALLOCATABLE                              :: refined(:, :)

      nc = SIZE(frame%ao%center_atoms)
      n = SIZE(weights)
      inverse_dh = frame%dh
      CALL invmat(inverse_dh, info)
      CPASSERT(info == 0)
      width = SQRT(SUM(inverse_dh**2, DIM=2))
      extent = subdivisions*frame%npts
      ALLOCATE (frame%stencil(interpolation_order, subdivisions, 3, nc), &
                frame%dstencil(interpolation_order, subdivisions, 3, nc))
      DO pass = 1, 2
         count = 0
         DO center = 1, nc
            kind = frame%ao%center_kinds(center)
            radius = frame%ao%delta_cache(kind)%radius
            IF (radius <= 0.0_dp) CYCLE
            origin = MATMUL(inverse_dh, frame%ao%positions(:, frame%ao%center_atoms(center)))
            box = CEILING(radius*width*subdivisions + 0.5_dp)
            IF (pass == 1) THEN
               DO axis = 1, 3
                  DO j = 1, subdivisions
                     fraction = MODULO(origin(axis) + (j - 0.5_dp)/subdivisions, 1.0_dp)
                     CALL cardinal_stencil(fraction, frame%stencil(:, j, axis, center), frame%dstencil(:, j, axis, center))
                  END DO
               END DO
            END IF
            DO z = -box(3), box(3)
               DO y = -box(2), box(2)
                  DO x = -box(1), box(1)
                     q = [x, y, z]
                     offset = (REAL(q, dp) + 0.5_dp)/subdivisions
                     IF (SUM(MATMUL(frame%dh, offset)**2) >= radius**2) CYCLE
                     IF (ANY(2*box + 1 > extent)) THEN
                        IF (previous_image(q, box, extent, frame%dh, radius)) CYCLE
                     END IF
                     coord = origin + offset
                     step = MODULO(q, subdivisions) + 1
                     ! Use the same fractional argument as the stencil, also at a grid plane.
                     base = FLOOR(origin + (REAL(step, dp) - 0.5_dp)/subdivisions) + (q - (step - 1))/subdivisions
                     DO axis = 1, 3
                        DO j = 1, interpolation_order
                           i = MODULO(base(axis) + j - interpolation_order/2, frame%npts(axis)) + frame%bounds(1, axis)
                           IF (i < frame%bounds_local(1, axis) .OR. i > frame%bounds_local(2, axis)) i = HUGE(0)
                           grid_indices(j, axis) = i
                        END DO
                     END DO
                     IF (ANY(ALL(grid_indices == HUGE(0), DIM=1))) CYCLE
                     count = count + 1
                     IF (pass == 1) CYCLE
                     frame%patch_index(:, :, count) = grid_indices
                     frame%patch_offset(:, count) = step
                     frame%patch_center(count) = center
                     refined(:, n + count) = MATMUL(frame%dh, MODULO(coord, REAL(frame%npts, dp)))
                     ! A unique coarse node owns each metric contribution across MPI quadrature halos.
                     owner = MODULO(base, frame%npts) + frame%bounds(1, :)
                     IF (ALL(owner >= frame%bounds_local(1, :)) .AND. ALL(owner <= frame%bounds_local(2, :))) THEN
                        frame%metric_weights(n + count) = frame%dvol/REAL(subdivisions**3, dp)
                     END IF
                  END DO
               END DO
            END DO
         END DO
         IF (pass == 1) THEN
            ALLOCATE (frame%patch_index(interpolation_order, 3, count), frame%patch_offset(3, count), frame%patch_center(count), &
                      refined(3, n + 2*count), frame%metric_weights(n + 2*count))
            refined(:, :n) = points
            frame%metric_weights = 0.0_dp
         END IF
      END DO
      refined(:, n + count + 1:) = refined(:, n + 1:n + count)
      CALL MOVE_ALLOC(refined, points)
      DEALLOCATE (weights)
      ALLOCATE (weights(n + 2*count))
      weights(:n) = frame%dvol
      weights(n + 1:n + count) = frame%dvol/REAL(subdivisions**3, dp)
      weights(n + count + 1:) = -weights(n + 1:n + count)
      frame%metric_weights(:n) = weights(:n)
      frame%metric_weights(n + count + 1:) = -frame%metric_weights(n + 1:n + count)
   END SUBROUTINE reconstruction_quadrature

!> \brief Keep each periodic fine node once when a reconstruction sphere intersects its own images.
!> \param q Fine-grid offset within the bounding box of the sphere.
!> \param box Maximum offsets in each direction.
!> \param extent Fine-grid periods.
!> \param dh Coarse-grid basis vectors.
!> \param radius Support of the radial reconstruction.
!> \return True if an earlier image contains the same quadrature point.
! **************************************************************************************************
   LOGICAL FUNCTION previous_image(q, box, extent, dh, radius) RESULT(duplicate)
      INTEGER, INTENT(IN)                                :: q(3), box(3), extent(3)
      REAL(dp), INTENT(IN)                               :: dh(3, 3), radius

      INTEGER                                            :: first(3), last(3), q2(3), x, y, z

      duplicate = .FALSE.
      first = -((box + q)/extent)
      last = (box - q)/extent
      DO z = first(3), MIN(0, last(3))
         DO y = first(2), last(2)
            DO x = first(1), last(1)
               IF (z == 0 .AND. (y > 0 .OR. (y == 0 .AND. x >= 0))) CYCLE
               q2 = q + extent*[x, y, z]
               IF (SUM(MATMUL(dh, (REAL(q2, dp) + 0.5_dp)/subdivisions)**2) < radius**2) THEN
                  duplicate = .TRUE.
                  RETURN
               END IF
            END DO
         END DO
      END DO
   END FUNCTION previous_image

! **************************************************************************************************
!> \brief Interpolate the resolved density spectrum with a moment-corrected cardinal B-spline.
!>        Values and two derivatives agree as a moving atom crosses a lattice plane.
!> \param t Fractional position between neighboring grid nodes.
!> \param weights Interpolation weights on the eight nearest nodes.
!> \param derivatives Derivative of the weights with respect to t.
! **************************************************************************************************
   SUBROUTINE cardinal_stencil(t, weights, derivatives)
      REAL(dp), INTENT(IN)                               :: t
      REAL(dp), INTENT(OUT) :: weights(interpolation_order), derivatives(interpolation_order)

      REAL(dp), PARAMETER :: second_moment = REAL(interpolation_order, dp)/24.0_dp

      INTEGER                                            :: cell, degree, j
      REAL(dp) :: b(-interpolation_order:3, 0:interpolation_order - 1), curvature, third

      ! Cardinal B-splines with their second moment removed. Low spatial frequencies
      ! are reproduced through cubic order while unresolved frequencies are attenuated.
      degree = interpolation_order - 1
      b = 0.0_dp
      b(0, 0) = 1.0_dp
      DO j = 1, degree
         DO cell = -j, 0
            b(cell, j) = ((t - cell)*b(cell, j - 1) + (cell + j + 1 - t)*b(cell + 1, j - 1))/j
         END DO
      END DO
      DO j = 1, interpolation_order
         cell = j - interpolation_order
         curvature = b(cell, degree - 2) - 2*b(cell + 1, degree - 2) + b(cell + 2, degree - 2)
         third = b(cell, degree - 3) - 3*b(cell + 1, degree - 3) + 3*b(cell + 2, degree - 3) - b(cell + 3, degree - 3)
         weights(j) = b(cell, degree) - second_moment*curvature
         derivatives(j) = b(cell, degree - 1) - b(cell + 1, degree - 1) - second_moment*third
      END DO
   END SUBROUTINE cardinal_stencil


! **************************************************************************************************
!> \brief Select the orbital quadrature, retaining the native samples for original UPF terms.
!> \param frame Prepared inverse frame.
!> \param points Orbital sampling coordinates.
!> \param weights Weights for contractions with interpolated potentials.
!> \param metric_weights Weights with unique ownership across overlapping MPI halos.
   SUBROUTINE orbital_samples(frame, points, weights, metric_weights)
      TYPE(qs_upf_inverse_grid_type), INTENT(IN), TARGET :: frame
      REAL(dp), INTENT(OUT), POINTER                     :: points(:, :), weights(:), &
                                                            metric_weights(:)

      points => frame%ao%points
      weights => frame%ao%weights
      metric_weights => frame%ao%weights
      IF (ALLOCATED(frame%orbital_points)) THEN
         points => frame%orbital_points
         weights => frame%orbital_weights
         metric_weights => frame%metric_weights
      END IF
   END SUBROUTINE orbital_samples

!> \brief Evaluate the same orbital partition for values, gradients and geometry derivatives.
!> \param frame Prepared inverse frame.
!> \param lo First sampling row.
!> \param hi Last sampling row.
!> \param soft Soft orbital samples in the chosen partition.
!> \param gradient Spatial gradients of the partitioned soft orbitals.
!> \param soft_response Joint derivative of the orbital samples.
!> \param gradient_response Joint derivative of their spatial gradients.
!> \param point_response Motion of the selected sampling points.
!> \param atom_response Atomic displacements.
!> \param cell_response Cell deformation.
!> \param reduced_k_response Change in the reduced Bloch vector.
!> \param coefficients_response Derivative of the original inverse coefficients.
   SUBROUTINE evaluate_rows(frame, lo, hi, soft, gradient, soft_response, gradient_response, &
                            point_response, atom_response, cell_response, reduced_k_response, coefficients_response)
      TYPE(qs_upf_inverse_grid_type), INTENT(IN)         :: frame
      INTEGER, INTENT(IN)                                :: lo, hi
      COMPLEX(dp), INTENT(OUT)                           :: soft(:, :)
      COMPLEX(dp), INTENT(OUT), OPTIONAL                 :: gradient(:, :, :), soft_response(:, :), &
                                                            gradient_response(:, :, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: point_response(:, :), &
                                                            atom_response(:, :), &
                                                            cell_response(3, 3), &
                                                            reduced_k_response(3)
      COMPLEX(dp), INTENT(IN), OPTIONAL                  :: coefficients_response(:, :)

      IF (ALLOCATED(frame%orbital_points)) THEN
         CALL qs_upf_inverse_evaluate(frame%ao, frame%orbital_points(:, lo:hi), soft, gradient, soft_response, &
                                      gradient_response, point_response, atom_response, cell_response, reduced_k_response, &
                                      coefficients_response, row_center=frame%row_center(lo:hi))
      ELSE
         CALL qs_upf_inverse_evaluate(frame%ao, frame%ao%points(:, lo:hi), soft, gradient, soft_response, &
                                      gradient_response, point_response, atom_response, cell_response, reduced_k_response, &
                                      coefficients_response)
      END IF
   END SUBROUTINE evaluate_rows

!> \brief Move each fine lattice with its atom relative to the strained native grid.
!> \param frame Prepared orbital quadrature.
!> \param strain Cartesian strain of the cell.
!> \param atom_response Atomic displacements.
!> \param motion Motion relative to the strained native grid in mesh coordinates.
!> \param dpoints Affine point displacements, updated to follow the atoms on the fine lattices.
   SUBROUTINE orbital_motion(frame, strain, atom_response, motion, dpoints)
      TYPE(qs_upf_inverse_grid_type), INTENT(IN)         :: frame
      REAL(dp), INTENT(IN)                               :: strain(3, 3), atom_response(:, :)
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: motion(:, :)
      REAL(dp), INTENT(INOUT)                            :: dpoints(:, :)

      INTEGER                                            :: atom, center, info, n, np, p
      REAL(dp)                                           :: inverse_dh(3, 3), relative(3)

      inverse_dh = frame%dh
      CALL invmat(inverse_dh, info)
      CPASSERT(info == 0)
      ALLOCATE (motion(3, SIZE(frame%ao%center_atoms)))
      DO center = 1, SIZE(frame%ao%center_atoms)
         atom = frame%ao%center_atoms(center)
         relative = atom_response(:, atom) - MATMUL(strain, frame%ao%positions(:, atom))
         motion(:, center) = MATMUL(inverse_dh, relative)
      END DO
      n = SIZE(frame%ao%weights)
      np = SIZE(frame%patch_center)
      DO p = 1, np
         center = frame%patch_center(p)
         dpoints(:, n + p) = dpoints(:, n + p) + MATMUL(frame%dh, motion(:, center))
         dpoints(:, n + np + p) = dpoints(:, n + p)
      END DO
   END SUBROUTINE orbital_motion

!> \brief Keep the real soft fields for repeated electronic operations on a Gamma grid.
!>        The full complex cache may exceed its budget even when these fields fit on their own.
!> \param frame Prepared grid and inverse coefficients at the current geometry.
! **************************************************************************************************
   SUBROUTINE cache_gamma_fields(frame)
      TYPE(qs_upf_inverse_grid_type), INTENT(INOUT)       :: frame

      COMPLEX(dp), ALLOCATABLE                           :: gradient(:, :, :), soft(:, :)
      INTEGER                                            :: block, first, last, m, nao, npoint, status

      IF ((frame%ao%fields_stored .AND. .NOT. ALLOCATED(frame%orbital_points)) .OR. ANY(frame%ao%reduced_k /= 0.0_dp)) RETURN
      nao = SIZE(frame%ao%coefficients, 2)
      npoint = SIZE(frame%ao%weights)
      IF (ALLOCATED(frame%orbital_weights)) npoint = SIZE(frame%orbital_weights)
      IF (npoint == 0 .OR. 32_int_8*nao*npoint > max_point_field_bytes) RETURN
      block = MIN(npoint, qs_upf_inverse_block_size(frame%ao))
      ALLOCATE (soft(block, nao), gradient(block, nao, 3), STAT=status)
      IF (status /= 0) RETURN
      ALLOCATE (frame%gamma_fields(npoint, nao, 4), STAT=status)
      IF (status /= 0) RETURN
      DO first = 1, npoint, block
         last = MIN(npoint, first + block - 1)
         m = last - first + 1
         CALL evaluate_rows(frame, first, last, soft(:m, :), gradient(:m, :, :))
         IF (ANY(AIMAG(soft(:m, :)) /= 0.0_dp) .OR. ANY(AIMAG(gradient(:m, :, :)) /= 0.0_dp)) THEN
            DEALLOCATE (frame%gamma_fields)
            RETURN
         END IF
         frame%gamma_fields(first:last, :, 1) = REAL(soft(:m, :), dp)
         frame%gamma_fields(first:last, :, 2:4) = REAL(gradient(:m, :, :), dp)
      END DO
   END SUBROUTINE cache_gamma_fields

! **************************************************************************************************
!> \brief Collocate one spin/k contribution and its projector density in the same inverse frame.
!>        The caller supplies occupations/k weights in P and optional original charge augmentation.
!>        Outputs are overwritten; this routine neither sums k points nor updates QS validity flags.
!> \param frame Prepared native-grid inverse frame.
!> \param density_real Native real AO density component in reconstructed Gaussian variables.
!> \param density Density on an allocated native PW grid; includes charge augmentation if supplied.
!> \param gradient Optional three physical Cartesian gradient components.
!> \param tau Smooth kinetic density, including the factor one half.
!> \param density_imag Optional native imaginary AO density component, with its actual sign.
!> \param projector_density Optional W P W^dagger for remaining original augmentation.
!> \param augmentation Optional owned original scalar USPP charge/gradient kernels.
!>        Tau remains smooth; supplying charge augmentation does not supply kinetic augmentation.
!> \param fourier_charge Optional native Fourier Q map, instead of point Q; default is false.
!>        Its zero mode equals the reference radial q; a supplied-q metric retains the explicit defect.
!> \param kinetic Optional independent original-wave valence kinetic kernels, indexed by kind.
!> \param orbital_fields Include soft-orbital fields. False returns only the requested augmentation.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_grid_collocate(frame, density_real, density, gradient, tau, density_imag, &
                                            projector_density, augmentation, fourier_charge, kinetic, orbital_fields)
      TYPE(qs_upf_inverse_grid_type), INTENT(IN)         :: frame
      TYPE(dbcsr_type), INTENT(IN)                       :: density_real
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: density, tau
      TYPE(pw_r3d_rs_type), INTENT(INOUT), OPTIONAL      :: gradient(3)
      TYPE(dbcsr_type), INTENT(IN), OPTIONAL             :: density_imag
      COMPLEX(dp), INTENT(OUT), OPTIONAL                 :: projector_density(:, :)
      TYPE(qs_upf_inverse_augmentation_type), &
         INTENT(INOUT), OPTIONAL                         :: augmentation
      LOGICAL, INTENT(IN), OPTIONAL                      :: fourier_charge
      TYPE(upf_inverse_kinetic_type), INTENT(IN), &
         OPTIONAL                                        :: kinetic(:)
      LOGICAL, INTENT(IN), OPTIONAL                      :: orbital_fields

      COMPLEX(dp), ALLOCATABLE                           :: p(:, :)

      CALL qs_upf_inverse_matrix_gather(frame%ao, density_real, p, density_imag)
      CALL qs_upf_inverse_grid_collocate_dense(frame, p, density, gradient, tau, &
                                               projector_density, augmentation, fourier_charge, kinetic, orbital_fields)
   END SUBROUTINE qs_upf_inverse_grid_collocate

! **************************************************************************************************
!> \brief The same native density map for a replicated complex AO density matrix.
!>        The caller supplies occupations and any k weight exactly once.
!> \param frame ...
!> \param p Hermitian AO density in the inverse variables.
!> \param density ...
!> \param gradient Optional Cartesian density gradient on the native grid.
!> \param tau Smooth kinetic density only.
!> \param projector_density ...
!> \param augmentation ...
!> \param fourier_charge ...
!> \param kinetic Optional independent valence kinetic augmentation; omitted for the bare soft kinetic operator.
!> \param orbital_fields ...
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_grid_collocate_dense(frame, p, density, gradient, tau, &
                                                  projector_density, augmentation, fourier_charge, kinetic, orbital_fields)
      TYPE(qs_upf_inverse_grid_type), INTENT(IN), TARGET :: frame
      COMPLEX(dp), INTENT(IN)                            :: p(:, :)
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: density, tau
      TYPE(pw_r3d_rs_type), INTENT(INOUT), OPTIONAL      :: gradient(3)
      COMPLEX(dp), INTENT(OUT), OPTIONAL                 :: projector_density(:, :)
      TYPE(qs_upf_inverse_augmentation_type), &
         INTENT(INOUT), OPTIONAL                         :: augmentation
      LOGICAL, INTENT(IN), OPTIONAL                      :: fourier_charge
      TYPE(upf_inverse_kinetic_type), INTENT(IN), &
         OPTIONAL                                        :: kinetic(:)
      LOGICAL, INTENT(IN), OPTIONAL                      :: orbital_fields

      COMPLEX(dp), ALLOCATABLE                           :: projected(:, :), sg(:, :, :), soft(:, :)
      INTEGER                                            :: axis, block, hi, lo, m, n, nao, nproj
      LOGICAL                                            :: use_fourier, use_orbitals
      REAL(dp), ALLOCATABLE                              :: dn(:, :), rho(:), t(:), t_augmentation(:)

      REAL(dp), POINTER :: points(:, :), weights(:), metric_weights(:)
      INTEGER :: nc

      CALL orbital_samples(frame, points, weights, metric_weights)
      nc = SIZE(frame%ao%weights)
      CPASSERT(frame%ao%atom_local .OR. frame%ao%inverse%ready)
      use_fourier = .FALSE.
      IF (PRESENT(fourier_charge)) use_fourier = fourier_charge
      use_orbitals = .TRUE.
      IF (PRESENT(orbital_fields)) use_orbitals = orbital_fields
      IF (use_fourier) THEN
         CPASSERT(PRESENT(augmentation))
      END IF
      CPASSERT(ALL(SHAPE(p) == SIZE(frame%ao%soft, 2)))
      n = SIZE(weights)
      nao = SIZE(frame%ao%coefficients, 2)
      block = qs_upf_inverse_block_size(frame%ao)
      ALLOCATE (rho(n), t(n))
      IF (PRESENT(gradient)) ALLOCATE (dn(n, 3), SOURCE=0.0_dp)
      rho = 0.0_dp
      t = 0.0_dp
      IF (use_orbitals) THEN
         ALLOCATE (soft(block, nao), sg(block, nao, 3))
         DO lo = 1, n, block
            hi = MIN(n, lo + block - 1)
            m = hi - lo + 1
            CALL point_fields(frame, lo, hi, soft(:m, :), sg(:m, :, :))
            ! The complete orbital density equals the Gaussian density plus its inverse residual.
            IF (PRESENT(gradient)) THEN
               CALL upf_inverse_smooth_density(soft(:m, :), sg(:m, :, :), p, rho(lo:hi), dn(lo:hi, :), t(lo:hi))
            ELSE
               CALL upf_inverse_smooth_density(soft(:m, :), sg(:m, :, :), p, rho(lo:hi), tau=t(lo:hi))
            END IF
         END DO
      END IF
      IF (PRESENT(projector_density) .OR. PRESENT(augmentation)) THEN
         nproj = SIZE(frame%ao%coefficients, 1)
         ALLOCATE (projected(nproj, nproj))
         CALL upf_inverse_projector_density(frame%ao%coefficients, p, projected, frame%ao%point_comm)
         IF (PRESENT(augmentation) .AND. .NOT. use_fourier) THEN
            IF (ALLOCATED(dn)) THEN
               CALL qs_upf_inverse_augmentation_collocate(frame%ao, augmentation, projected, rho(:nc), dn(:nc, :))
            ELSE
               CALL qs_upf_inverse_augmentation_collocate(frame%ao, augmentation, projected, rho(:nc))
            END IF
         END IF
         IF (PRESENT(projector_density)) THEN
            CPASSERT(ALL(SHAPE(projector_density) == [nproj, nproj]))
            projector_density = projected
         END IF
      END IF
      IF (PRESENT(kinetic)) THEN
         ALLOCATE (t_augmentation(nc))
         CALL qs_upf_inverse_kinetic_collocate(frame%ao, kinetic, frame%ao%points, p, t_augmentation)
         t(:nc) = t(:nc) + t_augmentation
      END IF
      CALL unpack_grid(frame, rho, density)
      CALL unpack_grid(frame, t, tau)
      IF (PRESENT(gradient)) THEN
         DO axis = 1, 3
            CALL unpack_grid(frame, dn(:, axis), gradient(axis))
         END DO
      END IF
      IF (use_fourier) THEN
         CALL qs_upf_inverse_augmentation_fourier_collocate(frame%ao, augmentation, projected, density, gradient)
      END IF
   END SUBROUTINE qs_upf_inverse_grid_collocate_dense

! **************************************************************************************************
!> \brief Integrate the density/gradient/tau adjoint into owned native real S / imaginary A matrices.
!>        Inverse-mediated AO fill is retained, independently of the template's stored blocks.
!> \param frame ...
!> \param vdensity Unweighted potential conjugate to smooth density.
!> \param vgradient Three unweighted potentials conjugate to physical grad(n), not sigma.
!> \param vtau Unweighted potential conjugate to smooth tau.
!> \param template Native AO block layout/distribution; its values and sparsity are not copied.
!> \param operator_real New owned real symmetric matrix; release it after use.
!> \param operator_imag New owned imaginary antisymmetric matrix; release it after use.
!> \param projector_operator Optional Hermitian projector operator.
!>        Missing charge or kinetic augmentation must not be inferred from the other kernel.
!> \param augmentation Optional original charge/gradient adjoint on the same native point rows.
!> \param fourier_charge Select the native Fourier adjoint of the same density map; default is false.
!> \param kinetic Optional independent valence kinetic augmentation; omitted for the bare soft kinetic operator.
!> \param orbital_fields ...
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_grid_integrate(frame, vdensity, vgradient, vtau, template, &
                            operator_real, operator_imag, projector_operator, augmentation, fourier_charge, kinetic, orbital_fields)
      TYPE(qs_upf_inverse_grid_type), INTENT(IN)         :: frame
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: vdensity, vgradient(3), vtau
      TYPE(dbcsr_type), INTENT(IN)                       :: template
      TYPE(dbcsr_type), INTENT(INOUT)                    :: operator_real, operator_imag
      COMPLEX(dp), INTENT(IN), OPTIONAL                  :: projector_operator(:, :)
      TYPE(qs_upf_inverse_augmentation_type), &
         INTENT(INOUT), OPTIONAL                         :: augmentation
      LOGICAL, INTENT(IN), OPTIONAL                      :: fourier_charge
      TYPE(upf_inverse_kinetic_type), INTENT(IN), &
         OPTIONAL                                        :: kinetic(:)
      LOGICAL, INTENT(IN), OPTIONAL                      :: orbital_fields

      COMPLEX(dp), ALLOCATABLE                           :: h(:, :)

      CALL qs_upf_inverse_grid_integrate_dense(frame, vdensity, vgradient, vtau, h, &
                                               projector_operator, augmentation, fourier_charge, kinetic, orbital_fields)
      CALL qs_upf_inverse_matrix_scatter(frame%ao, h, template, operator_real, operator_imag)
   END SUBROUTINE qs_upf_inverse_grid_integrate

! **************************************************************************************************
!> \brief Complex AO adjoint, before its native matrix distribution is selected.
!> \param frame ...
!> \param vdensity ...
!> \param vgradient ...
!> \param vtau ...
!> \param h New replicated Hermitian operator, without occupation or k weights.
!> \param projector_operator Optional Hermitian projector operator.
!> \param augmentation ...
!> \param fourier_charge ...
!> \param kinetic Optional independent valence kinetic augmentation; omitted for the bare soft kinetic operator.
!> \param orbital_fields ...
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_grid_integrate_dense(frame, vdensity, vgradient, vtau, h, &
                                                  projector_operator, augmentation, fourier_charge, kinetic, orbital_fields)
      TYPE(qs_upf_inverse_grid_type), INTENT(IN), TARGET :: frame
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: vdensity, vgradient(3), vtau
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: h(:, :)
      COMPLEX(dp), INTENT(IN), OPTIONAL                  :: projector_operator(:, :)
      TYPE(qs_upf_inverse_augmentation_type), &
         INTENT(INOUT), OPTIONAL                         :: augmentation
      LOGICAL, INTENT(IN), OPTIONAL                      :: fourier_charge
      TYPE(upf_inverse_kinetic_type), INTENT(IN), &
         OPTIONAL                                        :: kinetic(:)
      LOGICAL, INTENT(IN), OPTIONAL                      :: orbital_fields

      COMPLEX(dp), ALLOCATABLE                           :: projected(:, :), residual(:, :), sg(:, :, :), soft(:, :)
      INTEGER                                            :: axis, block, hi, lo, m, n, nao, nproj
      LOGICAL                                            :: use_fourier, use_gradients, use_orbitals
      REAL(dp), ALLOCATABLE                              :: dg(:, :), dn(:), dt(:)

      REAL(dp), POINTER :: points(:, :), weights(:), metric_weights(:)
      INTEGER :: nc

      CALL orbital_samples(frame, points, weights, metric_weights)
      nc = SIZE(frame%ao%weights)
      CPASSERT(frame%ao%atom_local .OR. frame%ao%inverse%ready)
      use_orbitals = .TRUE.
      IF (PRESENT(orbital_fields)) use_orbitals = orbital_fields
      use_fourier = .FALSE.
      IF (PRESENT(fourier_charge)) use_fourier = fourier_charge
      IF (use_fourier) THEN
         CPASSERT(PRESENT(augmentation))
      END IF
      n = SIZE(weights)
      nao = SIZE(frame%ao%soft, 2)
      block = qs_upf_inverse_block_size(frame%ao)
      ALLOCATE (dn(n), dg(n, 3), dt(n), h(nao, nao), residual(nao, nao))
      CALL pack_grid(frame, vdensity, dn)
      CALL pack_grid(frame, vtau, dt)
      DO axis = 1, 3
         CALL pack_grid(frame, vgradient(axis), dg(:, axis))
      END DO
      h = CMPLX(0.0_dp, 0.0_dp, dp)
      IF (use_orbitals) THEN
         use_gradients = ANY(dg /= 0.0_dp) .OR. ANY(dt /= 0.0_dp)
         ALLOCATE (soft(block, nao))
         IF (use_gradients) ALLOCATE (sg(block, nao, 3))
         DO lo = 1, n, block
            hi = MIN(n, lo + block - 1)
            m = hi - lo + 1
            IF (use_gradients) THEN
               CALL point_fields(frame, lo, hi, soft(:m, :), sg(:m, :, :))
               CALL upf_inverse_semilocal_pair(soft(:m, :), soft(:m, :), sg(:m, :, :), sg(:m, :, :), &
                                               weights(lo:hi), dn(lo:hi), dg(lo:hi, :), dt(lo:hi), residual)
            ELSE
               CALL point_fields(frame, lo, hi, soft(:m, :))
               CALL upf_inverse_pair(soft(:m, :), soft(:m, :), weights(lo:hi)*dn(lo:hi), residual)
            END IF
            h(:, :) = h + residual
         END DO
      END IF
      CALL frame%ao%point_comm%sum(h)
      IF (PRESENT(kinetic)) THEN
         CALL qs_upf_inverse_kinetic_operator(frame%ao, kinetic, frame%ao%points, frame%ao%weights, dt(:nc), &
                                              residual, frame%ao%point_comm)
         h(:, :) = h + residual
      END IF
      IF (PRESENT(augmentation)) THEN
         nproj = SIZE(frame%ao%coefficients, 1)
         ALLOCATE (projected(nproj, nproj))
         IF (use_fourier) THEN
            CALL qs_upf_inverse_augmentation_fourier_integrate(frame%ao, augmentation, vdensity, vgradient, projected)
         ELSE
            CALL qs_upf_inverse_augmentation_integrate(frame%ao, augmentation, dn(:nc), dg(:nc, :), projected)
         END IF
         CALL upf_inverse_projector_pair(frame%ao%coefficients, projected, frame%ao%coefficients, &
                                         residual, frame%ao%point_comm)
         h(:, :) = h + residual
      END IF
      IF (PRESENT(projector_operator)) THEN
         nproj = SIZE(frame%ao%coefficients, 1)
         CPASSERT(ALL(SHAPE(projector_operator) == [nproj, nproj]))
         CALL upf_inverse_projector_pair(frame%ao%coefficients, projector_operator, frame%ao%coefficients, &
                                         residual, frame%ao%point_comm)
         h(:, :) = h + residual
      END IF
      ! The variational adjoint acts on Hermitian density matrices. Independent
      ! contractions of its two triangles can lose symmetry when terms cancel.
      h(:, :) = 0.5_dp*(h + CONJG(TRANSPOSE(h)))
   END SUBROUTINE qs_upf_inverse_grid_integrate_dense

! **************************************************************************************************
!> \brief Full-Gaussian scalar-potential matrix on exactly the defining inverse grid and images.
!>        This is the Gaussian base alone, independent of W or the reconstruction. It permits a
!>        native analytic operator to replace its grid base while retaining the complete inverse
!>        and original Q contributions. Keep the same point ownership and integration measure.
!> \param frame Native grid coordinates, weights, full Gaussian bases and physical images.
!> \param potential Scalar potential at the defining grid points.
!> \param OPERATOR Optional new Gaussian matrix, without augmentation or reconstruction terms.
!> \param response Optional new joint geometry and Bloch response.
!> \param potential_response Required with response, the potential derivative at moving grid points.
!> \param atom_response Required with response, motion of the full Gaussian centers.
!> \param cell_response Required with response, motion of the cell vectors.
!> \param reduced_k_response Required with response, motion of the reduced Bloch vector.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_grid_gaussian_operator(frame, potential, OPERATOR, response, potential_response, &
                                                    atom_response, cell_response, reduced_k_response)
      TYPE(qs_upf_inverse_grid_type), INTENT(IN)         :: frame
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: potential
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT), OPTIONAL    :: operator(:, :), response(:, :)
      TYPE(pw_r3d_rs_type), INTENT(IN), OPTIONAL         :: potential_response
      REAL(dp), INTENT(IN), OPTIONAL                     :: atom_response(:, :), &
                                                            cell_response(3, 3), &
                                                            reduced_k_response(3)

      COMPLEX(dp), ALLOCATABLE                           :: gaussian(:, :), gaussian_response(:, :), &
                                                            partial(:, :), weighted(:, :), &
                                                            weighted_response(:, :)
      INTEGER                                            :: atom, block, first, hi, info, last, lo, &
                                                            m, n, nao
      REAL(dp)                                           :: cell_inverse(3, 3), strain(3, 3), &
                                                            volume_response
      REAL(dp), ALLOCATABLE                              :: points_response(:, :), values(:), &
                                                            values_response(:), weights_response(:)

      CPASSERT(PRESENT(OPERATOR) .OR. PRESENT(response))
      CPASSERT(SIZE(frame%ao%ao_last) > 0)
      nao = MAXVAL(frame%ao%ao_last)
      n = SIZE(frame%ao%weights)
      CPASSERT(SIZE(frame%ao%points, 2) == n)
      block = qs_upf_inverse_block_size(frame%ao)
      ALLOCATE (values(n), gaussian(block, nao), weighted(block, nao), partial(nao, nao))
      CALL pack_grid(frame, potential, values)
      IF (PRESENT(OPERATOR)) THEN
         ALLOCATE (OPERATOR(nao, nao), SOURCE=CMPLX(0.0_dp, 0.0_dp, dp))
      END IF
      IF (PRESENT(response)) THEN
         CPASSERT(PRESENT(potential_response))
         CPASSERT(PRESENT(atom_response))
         CPASSERT(PRESENT(cell_response))
         CPASSERT(PRESENT(reduced_k_response))
         CPASSERT(SIZE(atom_response, 1) == 3)
         CPASSERT(SIZE(atom_response, 2) == SIZE(frame%ao%bases))
         ALLOCATE (response(nao, nao), SOURCE=CMPLX(0.0_dp, 0.0_dp, dp))
         ALLOCATE (values_response(n), gaussian_response(block, nao), weighted_response(block, nao), &
                   points_response(3, block), weights_response(block))
         CALL pack_grid(frame, potential_response, values_response)
         cell_inverse = frame%ao%cell
         CALL invmat(cell_inverse, info)
         CPASSERT(info == 0)
         strain = MATMUL(cell_response, cell_inverse)
         volume_response = strain(1, 1) + strain(2, 2) + strain(3, 3)
      END IF
      DO lo = 1, n, block
         hi = MIN(n, lo + block - 1)
         m = hi - lo + 1
         IF (PRESENT(response)) THEN
            points_response(:, :m) = MATMUL(strain, frame%ao%points(:, lo:hi))
            weights_response(:m) = volume_response*frame%ao%weights(lo:hi)
         END IF
         IF (PRESENT(response) .OR. .NOT. frame%ao%fields_stored) THEN
            DO atom = 1, SIZE(frame%ao%bases)
               first = frame%ao%ao_first(atom)
               last = frame%ao%ao_last(atom)
               IF (PRESENT(response)) THEN
                  CALL gaussian_basis_fields_periodic(frame%ao%bases(atom)%gto_basis_set, frame%ao%points(:, lo:hi), &
                                                  frame%ao%positions(:, atom), frame%ao%cell, frame%ao%images, frame%ao%reduced_k, &
                                                      gaussian(:m, first:last), response=gaussian_response(:m, first:last), &
                                                    point_response=points_response(:, :m), center_response=atom_response(:, atom), &
                                                      cell_response=cell_response, reduced_k_response=reduced_k_response)
               ELSE
                  CALL gaussian_basis_fields_periodic(frame%ao%bases(atom)%gto_basis_set, frame%ao%points(:, lo:hi), &
                                                  frame%ao%positions(:, atom), frame%ao%cell, frame%ao%images, frame%ao%reduced_k, &
                                                      gaussian(:m, first:last))
               END IF
            END DO
         END IF
         IF (frame%ao%fields_stored) gaussian(:m, :) = frame%ao%gaussian(lo:hi, :)
         weighted(:m, :) = SPREAD(values(lo:hi), 2, nao)*gaussian(:m, :)
         IF (PRESENT(OPERATOR)) THEN
            CALL upf_inverse_pair(gaussian(:m, :), weighted(:m, :), frame%ao%weights(lo:hi), partial)
            OPERATOR(:, :) = OPERATOR + partial
         END IF
         IF (PRESENT(response)) THEN
            weighted_response(:m, :) = SPREAD(values(lo:hi), 2, nao)*gaussian_response(:m, :)
            CALL upf_inverse_pair_response(gaussian(:m, :), weighted(:m, :), frame%ao%weights(lo:hi), &
                                           gaussian_response(:m, :), weighted_response(:m, :), weights_response(:m), partial)
            response(:, :) = response + partial
            weighted_response(:m, :) = SPREAD(values_response(lo:hi), 2, nao)*gaussian(:m, :)
            CALL upf_inverse_pair(gaussian(:m, :), weighted_response(:m, :), frame%ao%weights(lo:hi), partial)
            response(:, :) = response + partial
         END IF
      END DO
      IF (PRESENT(OPERATOR)) CALL frame%ao%point_comm%sum(OPERATOR)
      IF (PRESENT(response)) CALL frame%ao%point_comm%sum(response)
   END SUBROUTINE qs_upf_inverse_grid_gaussian_operator

! **************************************************************************************************
!> \brief Assemble the original generalized metric in the inverse AO variables.
!>        The point rule is shared with the density/operator map; no metric residual is removed.
!> \param frame ...
!> \param projector_metric Explicit original projector metric (supplied q or reference radial q).
!>        The caller selects the comparison functional; this routine does not substitute either.
!> \param template ...
!> \param metric_real New owned real symmetric metric.
!> \param metric_imag New owned imaginary antisymmetric metric.
!> \param orbital_fields ...
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_grid_metric(frame, projector_metric, template, metric_real, metric_imag, orbital_fields)
      TYPE(qs_upf_inverse_grid_type), INTENT(IN)         :: frame
      COMPLEX(dp), INTENT(IN)                            :: projector_metric(:, :)
      TYPE(dbcsr_type), INTENT(IN)                       :: template
      TYPE(dbcsr_type), INTENT(INOUT)                    :: metric_real, metric_imag
      LOGICAL, INTENT(IN), OPTIONAL                      :: orbital_fields

      COMPLEX(dp), ALLOCATABLE                           :: metric(:, :)

      CALL qs_upf_inverse_grid_metric_dense(frame, projector_metric, metric, orbital_fields)
      CALL qs_upf_inverse_matrix_scatter(frame%ao, metric, template, metric_real, metric_imag)
   END SUBROUTINE qs_upf_inverse_grid_metric

! **************************************************************************************************
!> \brief Replicated complex inverse metric with the caller's explicit projector metric.
!> \param frame ...
!> \param projector_metric ...
!> \param metric New owned dense metric; radial Q and supplied PP_Q remain distinct.
!> \param orbital_fields ...
!> \param right_frame Optional independent right basis in the same physical projector and point frame.
!>        Each basis retains its own fixed Gaussian image set, including its geometry response.
!> \param response Optional joint geometry derivative of the rectangular metric.
!> \param atom_response Cartesian motion of every atom, required with response.
!> \param cell_response Lattice motion at fixed image membership, required with response.
!> \param reduced_k_response Reduced Bloch-vector motion, required with response.
!> \param point_response Motion of the shared defining quadrature, required with response.
!> \param weights_response Response of its measure, required with response.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_grid_metric_dense(frame, projector_metric, metric, orbital_fields, right_frame, response, &
                                               atom_response, cell_response, reduced_k_response, point_response, weights_response)
      TYPE(qs_upf_inverse_grid_type), INTENT(IN), TARGET :: frame
      COMPLEX(dp), INTENT(IN)                            :: projector_metric(:, :)
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: metric(:, :)
      LOGICAL, INTENT(IN), OPTIONAL                      :: orbital_fields
      TYPE(qs_upf_inverse_grid_type), INTENT(IN), &
         OPTIONAL, TARGET                                :: right_frame
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT), OPTIONAL    :: response(:, :)
      REAL(dp), INTENT(IN), OPTIONAL :: atom_response(:, :), cell_response(3, 3), &
         reduced_k_response(3), point_response(:, :), weights_response(:)

      COMPLEX(dp), ALLOCATABLE                           :: dr(:, :), ds(:, :), gradient(:, :, :), &
                                                            partial(:, :), right_gradient(:, :, :)
      COMPLEX(dp), ALLOCATABLE, TARGET                   :: dw(:, :), right_dw(:, :), &
                                                            right_soft(:, :), soft(:, :)
      COMPLEX(dp), POINTER :: right_coefficients_response(:, :), right_values(:, :)
      INTEGER                                            :: block, comparison, hi, lo, m, nao, &
                                                            nproj, nr
      LOGICAL                                            :: use_orbitals
      TYPE(qs_upf_inverse_grid_type), POINTER            :: right

      REAL(dp), POINTER :: points(:, :), weights(:), metric_weights(:)
      REAL(dp), ALLOCATABLE :: dpoints(:, :), dweights(:), motion(:, :)
      REAL(dp) :: inverse_cell(3, 3), strain(3, 3)
      INTEGER :: info

      CALL orbital_samples(frame, points, weights, metric_weights)
      use_orbitals = .TRUE.
      IF (PRESENT(orbital_fields)) use_orbitals = orbital_fields
      right => frame
      IF (PRESENT(right_frame)) THEN
         right => right_frame
         CPASSERT(ALLOCATED(frame%row_center) .EQV. ALLOCATED(right%row_center))
         IF (ALLOCATED(frame%row_center)) THEN
            CPASSERT(SIZE(frame%row_center) == SIZE(right%row_center))
            CPASSERT(ALL(frame%row_center == right%row_center))
         END IF
         comparison = frame%ao%point_comm%compare(right%ao%point_comm)
         CPASSERT(comparison == mp_comm_ident .OR. comparison == mp_comm_congruent)
         CPASSERT(ALL(SHAPE(frame%ao%points) == SHAPE(right%ao%points)))
         CPASSERT(ALL(SHAPE(frame%ao%positions) == SHAPE(right%ao%positions)))
         CPASSERT(ALL(frame%ao%points == right%ao%points))
         CPASSERT(ALL(frame%ao%weights == right%ao%weights))
         CPASSERT(ALL(frame%ao%positions == right%ao%positions))
         CPASSERT(ALL(frame%ao%cell == right%ao%cell))
         CPASSERT(ALL(frame%ao%reduced_k == right%ao%reduced_k))
         CPASSERT(SIZE(frame%ao%center_atoms) == SIZE(right%ao%center_atoms))
         CPASSERT(ALL(frame%ao%center_atoms == right%ao%center_atoms))
         CPASSERT(ALL(frame%ao%center_kinds == right%ao%center_kinds))
         CPASSERT(ALL(frame%ao%projector_first == right%ao%projector_first))
         CPASSERT(ALL(frame%ao%projector_last == right%ao%projector_last))
         CPASSERT(frame%ao%atom_local .EQV. right%ao%atom_local)
         CPASSERT(frame%ao%direct .EQV. right%ao%direct)
      END IF
      CPASSERT(frame%ao%atom_local .OR. frame%ao%inverse%ready)
      CPASSERT(right%ao%atom_local .OR. right%ao%inverse%ready)
      nao = SIZE(frame%ao%soft, 2)
      nr = SIZE(right%ao%soft, 2)
      nproj = SIZE(frame%ao%coefficients, 1)
      CPASSERT(ALL(SHAPE(projector_metric) == [nproj, nproj]))
      CPASSERT(SIZE(right%ao%coefficients, 1) == nproj)
      block = MIN(qs_upf_inverse_block_size(frame%ao), qs_upf_inverse_block_size(right%ao))
      ALLOCATE (metric(nao, nr), partial(nao, nr))
      metric = CMPLX(0.0_dp, 0.0_dp, dp)
      IF (PRESENT(response)) THEN
         CPASSERT(PRESENT(atom_response) .AND. PRESENT(cell_response))
         CPASSERT(PRESENT(reduced_k_response))
         CPASSERT(PRESENT(point_response) .AND. PRESENT(weights_response))
         ALLOCATE (response(nao, nr), dw(nproj, nao))
         response = CMPLX(0.0_dp, 0.0_dp, dp)
         CALL qs_upf_inverse_response(frame%ao, point_response, atom_response, cell_response, reduced_k_response, &
                                      weights_response, dw)
         right_coefficients_response => dw
         IF (PRESENT(right_frame)) THEN
            ALLOCATE (right_dw(nproj, nr))
            CALL qs_upf_inverse_response(right%ao, point_response, atom_response, cell_response, reduced_k_response, &
                                         weights_response, right_dw)
            right_coefficients_response => right_dw
         END IF
      END IF
      IF (use_orbitals) THEN
         ALLOCATE (soft(block, nao), gradient(block, nao, 3))
         right_values => soft
         IF (PRESENT(right_frame)) THEN
            ALLOCATE (right_soft(block, nr), right_gradient(block, nr, 3))
         END IF
         IF (PRESENT(response)) THEN
            ALLOCATE (ds(block, nao), dr(block, nr))
            IF (ALLOCATED(frame%orbital_points)) THEN
               inverse_cell = frame%ao%cell
               CALL invmat(inverse_cell, info)
               CPASSERT(info == 0)
               strain = MATMUL(cell_response, inverse_cell)
               ALLOCATE (dpoints, SOURCE=MATMUL(strain, points))
               ALLOCATE (dweights, SOURCE=(strain(1, 1) + strain(2, 2) + strain(3, 3))*metric_weights)
               CALL orbital_motion(frame, strain, atom_response, motion, dpoints)
            ELSE
               ALLOCATE (dpoints, SOURCE=point_response)
               ALLOCATE (dweights, SOURCE=weights_response)
            END IF
         END IF
         DO lo = 1, SIZE(weights), block
            hi = MIN(SIZE(weights), lo + block - 1)
            m = hi - lo + 1
            IF (PRESENT(response)) THEN
               CALL evaluate_rows(frame, lo, hi, soft(:m, :), gradient(:m, :, :), &
                                            soft_response=ds(:m, :), point_response=dpoints(:, lo:hi), &
                                            atom_response=atom_response, cell_response=cell_response, &
                                            reduced_k_response=reduced_k_response, coefficients_response=dw)
            ELSE
               CALL point_fields(frame, lo, hi, soft(:m, :), gradient(:m, :, :))
            END IF
            IF (PRESENT(right_frame)) THEN
               IF (PRESENT(response)) THEN
                  CALL evaluate_rows(right, lo, hi, right_soft(:m, :), &
                                               right_gradient(:m, :, :), soft_response=dr(:m, :), &
                                               point_response=dpoints(:, lo:hi), atom_response=atom_response, &
                                               cell_response=cell_response, reduced_k_response=reduced_k_response, &
                                               coefficients_response=right_coefficients_response)
               ELSE
                  CALL point_fields(right, lo, hi, right_soft(:m, :), right_gradient(:m, :, :))
               END IF
               right_values => right_soft
            ELSE IF (PRESENT(response)) THEN
               dr(:m, :) = ds(:m, :)
            END IF
            CALL upf_inverse_pair(soft(:m, :), right_values(:m, :), metric_weights(lo:hi), partial)
            metric(:, :) = metric + partial
            IF (PRESENT(response)) THEN
               CALL upf_inverse_pair_response(soft(:m, :), right_values(:m, :), metric_weights(lo:hi), &
                                              ds(:m, :), dr(:m, :), dweights(lo:hi), partial)
               response(:, :) = response + partial
            END IF
         END DO
      END IF
      CALL frame%ao%point_comm%sum(metric)
      CALL upf_inverse_projector_pair(frame%ao%coefficients, projector_metric, right%ao%coefficients, &
                                      partial, frame%ao%point_comm)
      metric(:, :) = metric + partial
      IF (PRESENT(response)) THEN
         CALL frame%ao%point_comm%sum(response)
         CALL upf_inverse_projector_pair(dw, projector_metric, right%ao%coefficients, partial, frame%ao%point_comm)
         response(:, :) = response + partial
         CALL upf_inverse_projector_pair(frame%ao%coefficients, projector_metric, right_coefficients_response, &
                                         partial, frame%ao%point_comm)
         response(:, :) = response + partial
      END IF
      IF (.NOT. PRESENT(right_frame)) THEN
         ! The same-basis metric is Hermitian. Independent contractions of its two
         ! triangles can differ by roundoff when the smooth and augmentation terms cancel.
         metric(:, :) = 0.5_dp*(metric + CONJG(TRANSPOSE(metric)))
         IF (PRESENT(response)) response(:, :) = 0.5_dp*(response + CONJG(TRANSPOSE(response)))
      END IF
   END SUBROUTINE qs_upf_inverse_grid_metric_dense

! **************************************************************************************************
!> \brief Joint native-grid geometry response of the inverse metric and original-Q semilocal operator.
!>        Fixed fractional grid rows, image membership and reciprocal indices define the derivative.
!>        Both moving inverse legs, their Cartesian gradients and the cell measure are differentiated.
!> \param frame ...
!> \param augmentation Original Fourier charge kernels; does not supply kinetic augmentation.
!> \param projector_metric Explicit original supplied-q or radial-Q metric, unchanged by geometry.
!> \param vdensity Unweighted conjugate potential.
!> \param vgradient Potentials conjugate to physical Cartesian grad(n).
!> \param vtau Potential conjugate to smooth tau, with the factor one half.
!> \param atom_response Cartesian response of every native atom.
!> \param cell_response Lattice-vector response.
!> \param reduced_k_response Response of the reduced Bloch vector, zero for physical cell stress.
!> \param vdensity_response Potential response at each moving grid row.
!> \param vgradient_response ...
!> \param vtau_response ...
!> \param operator_response New replicated complex operator response, without k/spin weights.
!> \param metric_response New metric response in the same inverse AO variables.
!> \param projector_operator Optional original geometry-independent operator, e.g. Dion.
!> \param kinetic Optional independent valence kinetic augmentation; omitted for the bare soft kinetic operator.
!> \param orbital_fields ...
!> \param known_coefficients_response Optional response for this defining frame and the same geometry variation.
!> \param saved_coefficients_response Optional owned coefficient response for reuse within this variation.
!> \param compute_metric Evaluate the metric response. If false, leave it unallocated.
!> \param charge Optional initially empty Fourier Q operator and derivative. Reuse only across
!>        k points at the same geometry, with identical potentials and geometry variation.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_grid_response(frame, augmentation, projector_metric, vdensity, vgradient, vtau, &
                                           atom_response, cell_response, reduced_k_response, vdensity_response, &
                                           vgradient_response, vtau_response, operator_response, metric_response, &
                                           projector_operator, kinetic, orbital_fields, &
                                           known_coefficients_response, saved_coefficients_response, compute_metric, charge)
      TYPE(qs_upf_inverse_grid_type), INTENT(IN), TARGET :: frame
      TYPE(qs_upf_inverse_augmentation_type), &
         INTENT(INOUT)                                   :: augmentation
      COMPLEX(dp), INTENT(IN)                            :: projector_metric(:, :)
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: vdensity, vgradient(3), vtau
      REAL(dp), INTENT(IN)                               :: atom_response(:, :), &
                                                            cell_response(3, 3), &
                                                            reduced_k_response(3)
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: vdensity_response, &
                                                            vgradient_response(3), vtau_response
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: operator_response(:, :), &
                                                            metric_response(:, :)
      COMPLEX(dp), INTENT(IN), OPTIONAL                  :: projector_operator(:, :)
      TYPE(upf_inverse_kinetic_type), INTENT(IN), &
         OPTIONAL                                        :: kinetic(:)
      LOGICAL, INTENT(IN), OPTIONAL                      :: orbital_fields
      COMPLEX(dp), INTENT(IN), OPTIONAL, TARGET :: known_coefficients_response(:, :)
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT), OPTIONAL :: saved_coefficients_response(:, :)
      LOGICAL, INTENT(IN), OPTIONAL                      :: compute_metric
      TYPE(qs_upf_inverse_grid_charge_type), &
         INTENT(INOUT), OPTIONAL, TARGET                 :: charge

      CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_upf_inverse_grid_response'

      COMPLEX(dp), ALLOCATABLE :: gradient(:, :, :), gradient_response(:, :, :), &
         kinetic_operator(:, :), kinetic_response(:, :), partial(:, :), projected(:, :), &
         soft(:, :), soft_response(:, :)
      COMPLEX(dp), ALLOCATABLE, TARGET                   :: owned_dw(:, :)
      COMPLEX(dp), POINTER                               :: dw(:, :)
      INTEGER                                            :: atom, axis, block, first, hi, info, last, lo, m, n, nao, &
                                                            nproj, timing_handle
      LOGICAL                                            :: cached_atomic_response, use_gradients, use_metric, use_orbitals
      REAL(dp)                                           :: cell_inverse(3, 3), strain(3, 3), &
                                                            volume_response
      REAL(dp), ALLOCATABLE                              :: dg(:, :), dn(:), dt(:), dvg(:, :), &
                                                            dvn(:), dvt(:), points_response(:, :), &
                                                            weights_response(:)
      TYPE(qs_upf_inverse_grid_charge_type), POINTER     :: charge_terms
      TYPE(qs_upf_inverse_grid_charge_type), TARGET      :: owned_charge

      REAL(dp), POINTER :: points(:, :), weights(:), metric_weights(:)
      INTEGER :: nc
      REAL(dp), ALLOCATABLE :: motion(:, :)

      CALL orbital_samples(frame, points, weights, metric_weights)
      nc = SIZE(frame%ao%weights)
      CALL timeset(routineN, timing_handle)
      use_orbitals = .TRUE.
      IF (PRESENT(orbital_fields)) use_orbitals = orbital_fields
      use_metric = .TRUE.
      IF (PRESENT(compute_metric)) use_metric = compute_metric
      CPASSERT(frame%ao%atom_local .OR. frame%ao%inverse%ready)
      n = SIZE(weights)
      nao = SIZE(frame%ao%soft, 2)
      nproj = SIZE(frame%ao%coefficients, 1)
      CPASSERT(ALL(SHAPE(projector_metric) == [nproj, nproj]))
      block = qs_upf_inverse_block_size(frame%ao)
      ALLOCATE (dn(n), dg(n, 3), dt(n), dvn(n), dvg(n, 3), dvt(n), points_response(3, n), weights_response(n), &
                operator_response(nao, nao))
      IF (use_metric) ALLOCATE (metric_response(nao, nao), SOURCE=CMPLX(0.0_dp, 0.0_dp, dp))
      cell_inverse = frame%ao%cell
      CALL invmat(cell_inverse, info)
      CPASSERT(info == 0)
      strain = MATMUL(cell_response, cell_inverse)
      volume_response = strain(1, 1) + strain(2, 2) + strain(3, 3)
      points_response(:, :) = MATMUL(strain, points)
      weights_response(:) = volume_response*weights
      IF (PRESENT(known_coefficients_response)) THEN
         CPASSERT(ALL(SHAPE(known_coefficients_response) == [nproj, nao]))
         dw => known_coefficients_response
      ELSE
         ALLOCATE (owned_dw(nproj, nao))
         dw => owned_dw
         CALL qs_upf_inverse_response(frame%ao, points_response(:, :nc), atom_response, cell_response, reduced_k_response, &
                                      weights_response(:nc), dw)
      END IF
      CALL pack_grid(frame, vdensity, dn)
      CALL pack_grid(frame, vtau, dt)
      CALL pack_grid(frame, vdensity_response, dvn)
      CALL pack_grid(frame, vtau_response, dvt)
      DO axis = 1, 3
         CALL pack_grid(frame, vgradient(axis), dg(:, axis))
         CALL pack_grid(frame, vgradient_response(axis), dvg(:, axis))
      END DO
      IF (ALLOCATED(frame%orbital_points)) THEN
         CALL orbital_motion(frame, strain, atom_response, motion, points_response)
         IF (ANY(motion /= 0.0_dp)) THEN
            CALL transfer_patch(frame, vdensity, gathered=dvn, motion=motion)
            CALL transfer_patch(frame, vtau, gathered=dvt, motion=motion)
            DO axis = 1, 3
               CALL transfer_patch(frame, vgradient(axis), gathered=dvg(:, axis), motion=motion)
            END DO
         END IF
      END IF
      operator_response = CMPLX(0.0_dp, 0.0_dp, dp)
      IF (use_orbitals) THEN
         use_gradients = ANY(dg /= 0.0_dp) .OR. ANY(dt /= 0.0_dp) .OR. &
                         ANY(dvg /= 0.0_dp) .OR. ANY(dvt /= 0.0_dp)
         ! For direct Gaussian orbitals at fixed grid and Bloch vector, moving an
         ! atom differentiates only its basis functions, with the opposite sign
         ! to their spatial gradient. Kinetic and gradient terms need Hessians.
         cached_atomic_response = frame%ao%direct .AND. frame%ao%fields_stored .AND. .NOT. use_gradients .AND. &
                                  ALL(cell_response == 0.0_dp) .AND. ALL(reduced_k_response == 0.0_dp)
         ALLOCATE (soft_response(block, nao), soft(block, nao), partial(nao, nao))
         IF (use_gradients) ALLOCATE (gradient_response(block, nao, 3), gradient(block, nao, 3))
         DO lo = 1, n, block
            hi = MIN(n, lo + block - 1)
            m = hi - lo + 1
            IF (.NOT. use_gradients) THEN
               IF (cached_atomic_response) THEN
                  soft(:m, :) = frame%ao%soft(lo:hi, :)
                  soft_response(:m, :) = CMPLX(0.0_dp, 0.0_dp, dp)
                  DO atom = 1, SIZE(frame%ao%ao_first)
                     first = frame%ao%ao_first(atom)
                     last = frame%ao%ao_last(atom)
                     DO axis = 1, 3
                        IF (atom_response(axis, atom) == 0.0_dp) CYCLE
                        soft_response(:m, first:last) = soft_response(:m, first:last) - &
                           atom_response(axis, atom)*frame%ao%soft_gradient(lo:hi, first:last, axis)
                     END DO
                  END DO
               ELSE
                  CALL evaluate_rows(frame, lo, hi, soft(:m, :), &
                                               soft_response=soft_response(:m, :), point_response=points_response(:, lo:hi), &
                                               atom_response=atom_response, cell_response=cell_response, &
                                               reduced_k_response=reduced_k_response, coefficients_response=dw)
               END IF
               ! Differentiate the product of the quadrature measure and the scalar potential.
               CALL upf_inverse_pair_response(soft(:m, :), soft(:m, :), weights(lo:hi)*dn(lo:hi), &
                                              soft_response(:m, :), soft_response(:m, :), &
                                              weights_response(lo:hi)*dn(lo:hi) + weights(lo:hi)*dvn(lo:hi), partial)
            ELSE
               IF (ANY(dg(lo:hi, :) /= 0.0_dp) .OR. ANY(dt(lo:hi) /= 0.0_dp)) THEN
                  CALL evaluate_rows(frame, lo, hi, soft(:m, :), gradient(:m, :, :), &
                                               soft_response(:m, :), gradient_response(:m, :, :), points_response(:, lo:hi), &
                                               atom_response, cell_response, reduced_k_response, dw)
               ELSE
                  ! Gradient responses enter only through the unperturbed gradient and tau potentials.
                  gradient_response(:m, :, :) = CMPLX(0.0_dp, 0.0_dp, dp)
                  CALL evaluate_rows(frame, lo, hi, soft(:m, :), gradient(:m, :, :), &
                                               soft_response=soft_response(:m, :), point_response=points_response(:, lo:hi), &
                                               atom_response=atom_response, cell_response=cell_response, &
                                               reduced_k_response=reduced_k_response, coefficients_response=dw)
               END IF
               CALL upf_inverse_semilocal_response(soft(:m, :), soft(:m, :), gradient(:m, :, :), gradient(:m, :, :), &
                                                   weights(lo:hi), dn(lo:hi), dg(lo:hi, :), dt(lo:hi), &
                                                   soft_response(:m, :), soft_response(:m, :), &
                                                   gradient_response(:m, :, :), gradient_response(:m, :, :), &
                                                   weights_response(lo:hi), dvn(lo:hi), dvg(lo:hi, :), dvt(lo:hi), partial)
            END IF
            operator_response(:, :) = operator_response + partial
            IF (use_metric) THEN
               CALL upf_inverse_pair_response(soft(:m, :), soft(:m, :), metric_weights(lo:hi), &
                                              soft_response(:m, :), soft_response(:m, :), &
                                              volume_response*metric_weights(lo:hi), partial)
               metric_response(:, :) = metric_response + partial
            END IF
         END DO
      END IF
      CALL frame%ao%point_comm%sum(operator_response)
      IF (use_metric) CALL frame%ao%point_comm%sum(metric_response)
      IF (PRESENT(kinetic)) THEN
         ALLOCATE (kinetic_operator(nao, nao), kinetic_response(nao, nao))
         CALL qs_upf_inverse_kinetic_operator(frame%ao, kinetic, frame%ao%points, frame%ao%weights, dt(:nc), &
                                              kinetic_operator, frame%ao%point_comm, &
                                              kinetic_response, points_response(:, :nc), atom_response, cell_response, &
                                              dw, weights_response(:nc), dvt(:nc))
         operator_response(:, :) = operator_response + kinetic_response
      END IF
      charge_terms => owned_charge
      IF (PRESENT(charge)) charge_terms => charge
      IF (.NOT. ALLOCATED(charge_terms%operator)) THEN
         ALLOCATE (charge_terms%operator(nproj, nproj), charge_terms%derivative(nproj, nproj))
         CALL qs_upf_inverse_augmentation_fourier_response(frame%ao, augmentation, vdensity, vgradient, atom_response, &
                                                           cell_response, vdensity_response, vgradient_response, &
                                                           charge_terms%derivative, charge_terms%operator)
      END IF
      CPASSERT(ALL(SHAPE(charge_terms%operator) == [nproj, nproj]))
      ALLOCATE (projected, SOURCE=charge_terms%operator)
      IF (PRESENT(projector_operator)) THEN
         CPASSERT(ALL(SHAPE(projector_operator) == [nproj, nproj]))
         projected(:, :) = projected + projector_operator
      END IF
      operator_response(:, :) = operator_response + MATMUL(CONJG(TRANSPOSE(frame%ao%coefficients)), &
                                                           MATMUL(charge_terms%derivative, frame%ao%coefficients))
      operator_response(:, :) = operator_response + MATMUL(CONJG(TRANSPOSE(dw)), &
                                                           MATMUL(projected, frame%ao%coefficients)) + &
                                MATMUL(CONJG(TRANSPOSE(frame%ao%coefficients)), MATMUL(projected, dw))
      IF (use_metric) THEN
         metric_response(:, :) = metric_response + MATMUL(CONJG(TRANSPOSE(dw)), &
                                                          MATMUL(projector_metric, frame%ao%coefficients)) + &
                                 MATMUL(CONJG(TRANSPOSE(frame%ao%coefficients)), MATMUL(projector_metric, dw))
         metric_response(:, :) = 0.5_dp*(metric_response + CONJG(TRANSPOSE(metric_response)))
      END IF
      operator_response(:, :) = 0.5_dp*(operator_response + CONJG(TRANSPOSE(operator_response)))
      IF (PRESENT(saved_coefficients_response)) THEN
         IF (ALLOCATED(owned_dw)) THEN
            CALL MOVE_ALLOC(owned_dw, saved_coefficients_response)
         ELSE
            ALLOCATE (saved_coefficients_response, SOURCE=dw)
         END IF
      END IF
      CALL timestop(timing_handle)
   END SUBROUTINE qs_upf_inverse_grid_response

! **************************************************************************************************
!> \brief Read or recompute one bounded block with the original W, images and point coordinates.
!> \param frame Original prepared inverse, with optional stored fields.
!> \param first First local row, inclusive.
!> \param last Last local row, inclusive.
!> \param soft Transformed orbital values.
!> \param gradient Optional physical spatial gradients.
! **************************************************************************************************
   SUBROUTINE point_fields(frame, first, last, soft, gradient)
      TYPE(qs_upf_inverse_grid_type), INTENT(IN)         :: frame
      INTEGER, INTENT(IN)                                :: first, last
      COMPLEX(dp), INTENT(OUT)                           :: soft(:, :)
      COMPLEX(dp), INTENT(OUT), OPTIONAL                 :: gradient(:, :, :)

      IF (frame%ao%fields_stored .AND. .NOT. ALLOCATED(frame%orbital_points)) THEN
         soft(:, :) = frame%ao%soft(first:last, :)
         IF (PRESENT(gradient)) gradient(:, :, :) = frame%ao%soft_gradient(first:last, :, :)
      ELSE IF (ALLOCATED(frame%gamma_fields)) THEN
         soft(:, :) = frame%gamma_fields(first:last, :, 1)
         IF (PRESENT(gradient)) gradient(:, :, :) = frame%gamma_fields(first:last, :, 2:4)
      ELSE
         CALL evaluate_rows(frame, first, last, soft, gradient)
      END IF
   END SUBROUTINE point_fields

! **************************************************************************************************
!> \brief Gather native AO components collectively, using native desymmetrization for S/A signs.
!>        All matrix-communicator ranks participate; the dense result is replicated.
!> \param frame Native atom/AO order used to validate both components.
!> \param matrix_real Native real AO component.
!> \param dense New owned complex dense matrix, real + i*imaginary.
!> \param matrix_imag Optional native imaginary AO component.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_matrix_gather(frame, matrix_real, dense, matrix_imag)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix_real
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: dense(:, :)
      TYPE(dbcsr_type), INTENT(IN), OPTIONAL             :: matrix_imag

      INTEGER                                            :: n
      REAL(dp), ALLOCATABLE                              :: component(:, :)

      n = frame%ao_last(SIZE(frame%ao_last))
      ALLOCATE (dense(n, n), component(n, n))
      CALL gather_component(frame, matrix_real, component)
      dense(:, :) = CMPLX(component, 0.0_dp, dp)
      IF (PRESENT(matrix_imag)) THEN
         CALL gather_component(frame, matrix_imag, component)
         dense(:, :) = dense + CMPLX(0.0_dp, component, dp)
      END IF
   END SUBROUTINE qs_upf_inverse_matrix_gather

! **************************************************************************************************
!> \brief Store every nonzero AO block on its native owner, including fill absent from the template.
!> \param frame ...
!> \param dense Hermitian complex operator in the native AO order.
!> \param template ...
!> \param matrix_real New owned symmetric real component.
!> \param matrix_imag New owned antisymmetric imaginary component.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_matrix_scatter(frame, dense, template, matrix_real, matrix_imag)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      COMPLEX(dp), INTENT(IN)                            :: dense(:, :)
      TYPE(dbcsr_type), INTENT(IN)                       :: template
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_real, matrix_imag

      INTEGER                                            :: c, i, j, nc, nr, owner, r
      REAL(dp)                                           :: tolerance
      TYPE(mp_comm_type)                                 :: comm

      CALL check_matrix(frame, template)
      nr = frame%ao_last(SIZE(frame%ao_last))
      CPASSERT(ALL(SHAPE(dense) == [nr, nr]))
      tolerance = 128.0_dp*EPSILON(1.0_dp)*MAX(1.0_dp, MAXVAL(ABS(dense)))
      CPASSERT(MAXVAL(ABS(dense - CONJG(TRANSPOSE(dense)))) <= tolerance)
      CALL dbcsr_get_info(template, group=comm)
      CALL dbcsr_create(matrix_real, name="Inverse UPF real", template=template, matrix_type=dbcsr_type_symmetric)
      CALL dbcsr_create(matrix_imag, name="Inverse UPF imaginary", template=template, matrix_type=dbcsr_type_antisymmetric)
      DO j = 1, SIZE(frame%ao_first)
         c = frame%ao_first(j)
         nc = frame%ao_last(j) - c + 1
         DO i = 1, j
            r = frame%ao_first(i)
            nr = frame%ao_last(i) - r + 1
            CALL dbcsr_get_stored_coordinates(matrix_real, i, j, owner)
            IF (owner /= comm%mepos) CYCLE
            IF (ANY(REAL(dense(r:r + nr - 1, c:c + nc - 1), dp) /= 0.0_dp)) THEN
               CALL dbcsr_put_block(matrix_real, i, j, REAL(dense(r:r + nr - 1, c:c + nc - 1), dp))
            END IF
            IF (ANY(AIMAG(dense(r:r + nr - 1, c:c + nc - 1)) /= 0.0_dp)) THEN
               CALL dbcsr_put_block(matrix_imag, i, j, AIMAG(dense(r:r + nr - 1, c:c + nc - 1)))
            END IF
         END DO
      END DO
      CALL dbcsr_finalize(matrix_real)
      CALL dbcsr_finalize(matrix_imag)
   END SUBROUTINE qs_upf_inverse_matrix_scatter

! **************************************************************************************************
!> \brief ...
!> \param frame ...
!> \param matrix ...
!> \param dense ...
! **************************************************************************************************
   SUBROUTINE gather_component(frame, matrix, dense)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      REAL(dp), INTENT(OUT)                              :: dense(:, :)

      INTEGER                                            :: c, nc, nr, r
      LOGICAL                                            :: transposed
      REAL(dp), POINTER                                  :: block(:, :)
      TYPE(dbcsr_iterator_type)                          :: iterator
      TYPE(dbcsr_type)                                   :: full
      TYPE(mp_comm_type)                                 :: comm

      CALL check_matrix(frame, matrix)
      CALL dbcsr_get_info(matrix, group=comm)
      CALL dbcsr_create(full, template=matrix, matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_desymmetrize(matrix, full)
      dense = 0.0_dp
      CALL dbcsr_iterator_start(iterator, full)
      DO WHILE (dbcsr_iterator_blocks_left(iterator))
         CALL dbcsr_iterator_next_block(iterator, block=block, row_offset=r, col_offset=c, transposed=transposed)
         nr = SIZE(block, 1)
         nc = SIZE(block, 2)
         IF (transposed) THEN
            dense(r:r + nc - 1, c:c + nr - 1) = TRANSPOSE(block)
         ELSE
            dense(r:r + nr - 1, c:c + nc - 1) = block
         END IF
      END DO
      CALL dbcsr_iterator_stop(iterator)
      CALL dbcsr_release(full)
      CALL comm%sum(dense)
   END SUBROUTINE gather_component

! **************************************************************************************************
!> \brief ...
!> \param frame ...
!> \param matrix ...
! **************************************************************************************************
   SUBROUTINE check_matrix(frame, matrix)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix

      INTEGER, POINTER                                   :: cols(:), rows(:)

      CALL dbcsr_get_info(matrix, row_blk_size=rows, col_blk_size=cols)
      CPASSERT(SIZE(rows) == SIZE(frame%ao_first) .AND. SIZE(cols) == SIZE(frame%ao_first))
      CPASSERT(ALL(rows == frame%ao_last - frame%ao_first + 1))
      CPASSERT(ALL(cols == rows))
   END SUBROUTINE check_matrix

! **************************************************************************************************
!> \brief ...
!> \param frame ...
!> \param field ...
! **************************************************************************************************
   SUBROUTINE check_grid(frame, field)
      TYPE(qs_upf_inverse_grid_type), INTENT(IN)         :: frame
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: field

      INTEGER                                            :: extents(3)

      CPASSERT(ASSOCIATED(field%array) .AND. ASSOCIATED(field%pw_grid))
      CPASSERT(ALL(field%pw_grid%bounds == frame%bounds))
      CPASSERT(ALL(field%pw_grid%bounds_local == frame%bounds_local))
      CPASSERT(ALL(field%pw_grid%npts == frame%npts))
      CPASSERT(field%pw_grid%para%mode == frame%mode)
      CPASSERT(ALL(field%pw_grid%dh == frame%dh) .AND. field%pw_grid%dvol == frame%dvol)
      extents = MAX(0, frame%bounds_local(2, :) - frame%bounds_local(1, :) + 1)
      CPASSERT(ALL(SHAPE(field%array) == extents))
      IF (SIZE(field%array) > 0) THEN
         CPASSERT(ALL(LBOUND(field%array) == frame%bounds_local(1, :)))
         CPASSERT(ALL(UBOUND(field%array) == frame%bounds_local(2, :)))
      END IF
      IF (frame%mode == PW_MODE_DISTRIBUTED) THEN
         CPASSERT(field%pw_grid%para%group%mp_comm_type == frame%ao%point_comm)
      END IF
   END SUBROUTINE check_grid

! **************************************************************************************************
!> \brief ...
!> \param frame ...
!> \param field ...
!> \param points ...
! **************************************************************************************************
   SUBROUTINE pack_grid(frame, field, points)
      TYPE(qs_upf_inverse_grid_type), INTENT(IN)         :: frame
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: field
      REAL(dp), INTENT(OUT)                              :: points(:)

      INTEGER                                            :: n

      CALL check_grid(frame, field)
      n = SIZE(field%array)
      points(:n) = RESHAPE(field%array, [n])
      points(n + 1:) = 0.0_dp
      IF (SIZE(points) > n) CALL transfer_patch(frame, field, gathered=points)
   END SUBROUTINE pack_grid

! **************************************************************************************************
!> \brief ...
!> \param frame ...
!> \param points ...
!> \param field ...
! **************************************************************************************************
   SUBROUTINE unpack_grid(frame, points, field)
      TYPE(qs_upf_inverse_grid_type), INTENT(IN)         :: frame
      REAL(dp), INTENT(IN)                               :: points(:)
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: field

      INTEGER                                            :: n

      CALL check_grid(frame, field)
      n = SIZE(field%array)
      field%array = RESHAPE(points(:n), SHAPE(field%array))
      IF (SIZE(points) > n) CALL transfer_patch(frame, field, scattered=points)
   END SUBROUTINE unpack_grid

!> \brief Interpolate a grid field onto atomic samples or deposit their weighted transpose.
!> \param frame Atomic quadrature and interpolation stencils.
!> \param field Native grid field read during interpolation or incremented during deposition.
!> \param gathered Interpolated potential or its response at the fine samples.
!> \param scattered Density samples deposited by the weighted transpose.
!> \param motion Optional motion of the interpolation stencils in mesh coordinates.
   SUBROUTINE transfer_patch(frame, field, gathered, scattered, motion)
      TYPE(qs_upf_inverse_grid_type), INTENT(IN)         :: frame
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: field
      REAL(dp), INTENT(INOUT), OPTIONAL                  :: gathered(:)
      REAL(dp), INTENT(IN), OPTIONAL                     :: scattered(:), motion(:, :)

      INTEGER                                            :: a, axis, b, c, center, i, j, k, n, np, &
                                                            offset, p
      REAL(dp)                                           :: dw(interpolation_order, 3), v, &
                                                            w(interpolation_order, 3), &
                                                            wx(interpolation_order), yz

      CPASSERT(PRESENT(gathered) .NEQV. PRESENT(scattered))
      IF (PRESENT(gathered)) THEN
         IF (ALL(field%array == 0.0_dp)) RETURN
      END IF
      n = SIZE(frame%ao%weights)
      np = SIZE(frame%patch_center)
      DO p = 1, np
         center = frame%patch_center(p)
         IF (PRESENT(motion)) THEN
            IF (ALL(motion(:, center) == 0.0_dp)) CYCLE
         END IF
         v = 0.0_dp
         IF (PRESENT(scattered)) THEN
            v = (scattered(n + p) - scattered(n + np + p))/REAL(subdivisions**3, dp)
            IF (v == 0.0_dp) CYCLE
         END IF
         DO axis = 1, 3
            offset = frame%patch_offset(axis, p)
            w(:, axis) = frame%stencil(:, offset, axis, center)
            IF (PRESENT(motion)) dw(:, axis) = frame%dstencil(:, offset, axis, center)*motion(axis, center)
         END DO
         DO k = 1, interpolation_order
            c = frame%patch_index(k, 3, p)
            IF (c == HUGE(0)) CYCLE
            DO j = 1, interpolation_order
               b = frame%patch_index(j, 2, p)
               IF (b == HUGE(0)) CYCLE
               yz = w(j, 2)*w(k, 3)
               IF (PRESENT(motion)) THEN
                  wx = dw(:, 1)*yz + w(:, 1)*(dw(j, 2)*w(k, 3) + w(j, 2)*dw(k, 3))
               ELSE
                  wx = w(:, 1)*yz
               END IF
               IF (PRESENT(scattered)) THEN
                  DO i = 1, interpolation_order
                     a = frame%patch_index(i, 1, p)
                     IF (a /= HUGE(0)) field%array(a, b, c) = field%array(a, b, c) + v*wx(i)
                  END DO
               ELSE
                  DO i = 1, interpolation_order
                     a = frame%patch_index(i, 1, p)
                     IF (a /= HUGE(0)) v = v + field%array(a, b, c)*wx(i)
                  END DO
               END IF
            END DO
         END DO
         IF (PRESENT(gathered)) THEN
            gathered(n + p) = gathered(n + p) + v
            gathered(n + np + p) = gathered(n + np + p) + v
         END IF
      END DO
   END SUBROUTINE transfer_patch

! **************************************************************************************************
!> \brief Release owned inverse fields before their borrowed KS basis owner is destroyed.
!> \param frame ...
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_grid_release(frame)
      TYPE(qs_upf_inverse_grid_type), INTENT(INOUT)      :: frame

      TYPE(qs_upf_inverse_grid_type)                     :: empty

      CALL qs_upf_inverse_release(frame%ao)
      frame = empty
   END SUBROUTINE qs_upf_inverse_grid_release

END MODULE qs_upf_inverse_grid
