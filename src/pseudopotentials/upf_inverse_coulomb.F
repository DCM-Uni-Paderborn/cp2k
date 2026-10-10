!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
!> \brief Radial multipole Coulomb map for signed inverse reconstruction densities.
!>        All maps use the same finite radial measure, including the diagonal shell term.
!>        Periodic background and interactions between distinct centers are assembled separately.
MODULE upf_inverse_coulomb
   USE bessel_lib,                      ONLY: bessel0
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: fourpi,&
                                              pi,&
                                              twopi
   USE mathlib,                         ONLY: invmat
   USE message_passing,                 ONLY: mp_comm_type
   USE orbital_pointers,                ONLY: indco,&
                                              init_orbital_pointers,&
                                              nco,&
                                              ncoset
   USE orbital_transformation_matrices, ONLY: init_spherical_harmonics,&
                                              orbtramat
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   TYPE upf_inverse_coulomb_center_type
      REAL(dp) :: position(3) = 0.0_dp
      INTEGER, ALLOCATABLE :: l(:), indices(:)
      REAL(dp), ALLOCATABLE :: radii(:), weights(:), density(:, :), shapes(:, :)
   END TYPE upf_inverse_coulomb_center_type

   TYPE upf_inverse_coulomb_potential_type
      REAL(dp), ALLOCATABLE :: values(:, :)
   END TYPE upf_inverse_coulomb_potential_type

   TYPE center_work_type
      REAL(dp), ALLOCATABLE :: compensation(:, :), residual(:, :), q(:, :), dq(:, :), potential(:, :)
      COMPLEX(dp), ALLOCATABLE :: fourier(:, :, :)
   END TYPE center_work_type

   TYPE fourier_map_type
      COMPLEX(dp), ALLOCATABLE :: values(:, :, :)
   END TYPE fourier_map_type

   TYPE pair_map_type
      REAL(dp), ALLOCATABLE :: values(:, :, :, :)
   END TYPE pair_map_type

   TYPE upf_inverse_coulomb_map_type
      PRIVATE
      INTEGER :: periodicity(3) = [1, 1, 1]
      REAL(dp) :: cell(3, 3) = 0.0_dp, volume = 0.0_dp
      REAL(dp), ALLOCATABLE :: gvecs(:, :)
      TYPE(upf_inverse_coulomb_center_type), ALLOCATABLE :: centers(:)
      TYPE(fourier_map_type), ALLOCATABLE :: fourier(:)
      TYPE(pair_map_type), ALLOCATABLE :: pairs(:, :)
   END TYPE upf_inverse_coulomb_map_type

   PUBLIC :: upf_inverse_coulomb_center_type, upf_inverse_coulomb_potential_type
   PUBLIC :: upf_inverse_coulomb_periodic
   PUBLIC :: upf_inverse_coulomb_map_type, upf_inverse_coulomb_prepare
   PUBLIC :: upf_inverse_coulomb_radial, upf_inverse_coulomb_radial_response
   PUBLIC :: upf_inverse_coulomb_compensate, upf_inverse_coulomb_compensation_adjoint
   PUBLIC :: upf_inverse_coulomb_fourier_kernel
   PUBLIC :: upf_inverse_coulomb_pair_kernel

CONTAINS

! **************************************************************************************************
!> \brief Freeze the Fourier and complete periodic residual-pair maps at a fixed geometry.
!>        Every overlapping image is retained. Density and compensation moments remain variable.
!> \param centers Same center geometry and radial/harmonic layout as the subsequent functional.
!> \param cell ...
!> \param gvecs ...
!> \param volume ...
!> \param map Owned maps. Rebuild after changing centers, cell, reciprocal vectors or radial grids.
!> \param periodicity Optional lattice directions to repeat. Defaults to three periodic directions.
! **************************************************************************************************
   SUBROUTINE upf_inverse_coulomb_prepare(centers, cell, gvecs, volume, map, periodicity)
      TYPE(upf_inverse_coulomb_center_type), INTENT(IN)  :: centers(:)
      REAL(dp), INTENT(IN)                               :: cell(3, 3), gvecs(:, :), volume
      TYPE(upf_inverse_coulomb_map_type), INTENT(OUT)    :: map
      INTEGER, INTENT(IN), OPTIONAL                      :: periodicity(3)

      INTEGER                                            :: a, b, h, i, image(3), info, j, k, &
                                                            limits(2, 3), r
      REAL(dp)                                           :: displacement(3), inverse_cell(3, 3), &
                                                            radius
      REAL(dp), ALLOCATABLE                              :: kernel(:, :, :, :)

      CPASSERT(volume > 0.0_dp .AND. SIZE(gvecs, 1) == 3)
      IF (PRESENT(periodicity)) map%periodicity = periodicity
      CPASSERT(ALL(map%periodicity == 0 .OR. map%periodicity == 1))
      map%cell = cell
      map%volume = volume
      ALLOCATE (map%gvecs, SOURCE=gvecs)
      ALLOCATE (map%centers, SOURCE=centers)
      ALLOCATE (map%fourier(SIZE(centers)), map%pairs(SIZE(centers), SIZE(centers)))
      inverse_cell = cell
      CALL invmat(inverse_cell, info)
      CPASSERT(info == 0)
      DO a = 1, SIZE(centers)
         ALLOCATE (map%fourier(a)%values(SIZE(gvecs, 2), SIZE(centers(a)%radii), SIZE(centers(a)%l)))
         CALL upf_inverse_coulomb_fourier_kernel(centers(a)%radii, centers(a)%l, centers(a)%indices, &
                                                 gvecs, centers(a)%position, volume, map%fourier(a)%values)
         DO b = a, SIZE(centers)
            ALLOCATE (kernel(SIZE(centers(a)%radii), SIZE(centers(a)%l), &
                             SIZE(centers(b)%radii), SIZE(centers(b)%l)))
            ALLOCATE (map%pairs(a, b)%values, MOLD=kernel)
            map%pairs(a, b)%values = 0.0_dp
            radius = MAXVAL(centers(a)%radii) + MAXVAL(centers(b)%radii)
            CALL residual_image_limits(inverse_cell, centers(a)%position - centers(b)%position, &
                                       radius, map%periodicity, limits)
            DO k = limits(1, 3), limits(2, 3)
               DO j = limits(1, 2), limits(2, 2)
                  DO i = limits(1, 1), limits(2, 1)
                     image = [i, j, k]
                     displacement = centers(a)%position - centers(b)%position - MATMUL(cell, REAL(image, dp))
                     IF (SUM(displacement**2) >= radius**2) CYCLE
                     CALL upf_inverse_coulomb_pair_kernel(centers(a)%radii, centers(a)%l, centers(a)%indices, &
                                                          centers(b)%radii, centers(b)%l, centers(b)%indices, &
                                                          displacement, kernel)
                     map%pairs(a, b)%values(:, :, :, :) = map%pairs(a, b)%values + kernel
                  END DO
               END DO
            END DO
            IF (a /= b) THEN
               ALLOCATE (map%pairs(b, a)%values(SIZE(kernel, 3), SIZE(kernel, 4), SIZE(kernel, 1), SIZE(kernel, 2)))
               DO h = 1, SIZE(kernel, 2)
                  DO r = 1, SIZE(kernel, 1)
                     map%pairs(b, a)%values(:, :, r, h) = map%pairs(a, b)%values(r, h, :, :)
                  END DO
               END DO
            END IF
            DEALLOCATE (kernel)
         END DO
      END DO
   END SUBROUTINE upf_inverse_coulomb_prepare

! **************************************************************************************************
!> \brief Complete periodic Coulomb functional for a Fourier background and signed centered fields.
!>        Moment compensation uses the actual radial measure. Every overlapping residual image is
!>        retained, including self images. Separated residuals have exactly zero exterior potential.
!>        Gaussian hard/soft, reconstruction and original Q contributions belong in the supplied
!>        centered density. Ionic and Gaussian background charges belong in background.
!> \param centers Replicated harmonic densities, radial measures and compensation shapes.
!> \param cell Physical periodic cell, with lattice vectors in columns.
!> \param gvecs Local physical reciprocal vectors, in the same order as background.
!> \param volume Physical cell volume used in the Fourier coefficients.
!> \param background Fourier background before adding compensation charges.
!> \param poisson_weights Native Coulomb kernel times volume and Fourier multiplicity. The caller
!>        supplies the zero-mode convention. These weights are nonnegative.
!> \param energy Full Coulomb energy, excluding any separately defined ionic self energies.
!> \param potentials Radial-measure density adjoints, one for each supplied center.
!> \param background_adjoint Local Fourier covector, including the residual/background interaction.
!> \param comm Optional communicator distributing Fourier rows. Centered data are replicated.
!> \param response Optional joint energy derivative of this same finite representation.
!> \param center_response Radial, measure, density, shape and center derivatives, required with response.
!> \param cell_response Lattice derivative, required with response.
!> \param gvecs_response Physical reciprocal-vector derivative, required with response.
!> \param volume_response Cell-volume derivative, required with response.
!> \param background_response Background Fourier derivative, required with response.
!> \param poisson_weights_response Complete native Coulomb-weight derivative, required with response.
!> \param map Optional immutable maps for this exact geometry. Joint responses use the original
!>        differentiated kernels, so freezing a map never freezes the geometry derivative.
!> \param periodicity Optional lattice directions to repeat. Defaults to the map's periodicity
!>        when a map is supplied, otherwise to three periodic directions. The Fourier weights
!>        and their response must use the same electrostatic boundary conditions.
! **************************************************************************************************
   SUBROUTINE upf_inverse_coulomb_periodic(centers, cell, gvecs, volume, background, poisson_weights, &
                                           energy, potentials, background_adjoint, comm, response, center_response, &
                                           cell_response, gvecs_response, volume_response, background_response, &
                                           poisson_weights_response, map, periodicity)
      TYPE(upf_inverse_coulomb_center_type), INTENT(IN)  :: centers(:)
      REAL(dp), INTENT(IN)                               :: cell(3, 3), gvecs(:, :), volume
      COMPLEX(dp), INTENT(IN)                            :: background(:)
      REAL(dp), INTENT(IN)                               :: poisson_weights(:)
      REAL(dp), INTENT(OUT)                              :: energy
      TYPE(upf_inverse_coulomb_potential_type), &
         ALLOCATABLE, INTENT(OUT)                        :: potentials(:)
      COMPLEX(dp), INTENT(OUT)                           :: background_adjoint(:)
      TYPE(mp_comm_type), INTENT(IN), OPTIONAL           :: comm
      REAL(dp), INTENT(OUT), OPTIONAL                    :: response
      TYPE(upf_inverse_coulomb_center_type), &
         INTENT(IN), OPTIONAL                            :: center_response(:)
      REAL(dp), INTENT(IN), OPTIONAL                     :: cell_response(3, 3), &
                                                            gvecs_response(:, :), volume_response
      COMPLEX(dp), INTENT(IN), OPTIONAL                  :: background_response(:)
      REAL(dp), INTENT(IN), OPTIONAL                     :: poisson_weights_response(:)
      TYPE(upf_inverse_coulomb_map_type), INTENT(IN), &
         OPTIONAL                                        :: map
      INTEGER, INTENT(IN), OPTIONAL                      :: periodicity(3)

      COMPLEX(dp), ALLOCATABLE                           :: dfourier(:, :, :), dg(:), dh(:), e(:), &
                                                            h(:)
      INTEGER                                            :: a, b, channel, i, image(3), info, j, &
                                                            limits(2, 3), na, nb, ng, nh, nr, &
                                                            REPEAT(3)
      LOGICAL                                            :: derive, fixed_pair
      REAL(dp)                                           :: ddisplacement(3), displacement(3), &
                                                            dmoment, dnorm, dpair, dpw, epair, &
                                                            inverse_cell(3, 3), moment, norm, &
                                                            radius
      REAL(dp), ALLOCATABLE :: adjoint(:, :), dc(:, :), de(:, :), df(:), dkernel(:, :, :, :), &
         f(:), kernel(:, :, :, :), moments(:), ve(:, :), vh(:, :)
      TYPE(center_work_type), ALLOCATABLE                :: work(:)

      ng = SIZE(background)
      repeat = [1, 1, 1]
      IF (PRESENT(map)) repeat = map%periodicity
      IF (PRESENT(periodicity)) repeat = periodicity
      CPASSERT(ALL(repeat == 0 .OR. repeat == 1))
      CPASSERT(ALL(SHAPE(gvecs) == [3, ng]))
      CPASSERT(SIZE(poisson_weights) == ng)
      CPASSERT(SIZE(background_adjoint) == ng)
      CPASSERT(ALL(ieee_is_finite(poisson_weights)))
      CPASSERT(ALL(poisson_weights >= 0.0_dp))
      inverse_cell = cell
      CALL invmat(inverse_cell, info)
      CPASSERT(info == 0)
      derive = PRESENT(response)
      IF (derive) THEN
         CPASSERT(PRESENT(center_response) .AND. PRESENT(cell_response))
         CPASSERT(PRESENT(gvecs_response) .AND. PRESENT(volume_response))
         CPASSERT(PRESENT(background_response))
         CPASSERT(PRESENT(poisson_weights_response))
         CPASSERT(SIZE(center_response) == SIZE(centers))
         CPASSERT(ALL(SHAPE(gvecs_response) == SHAPE(gvecs)))
         CPASSERT(SIZE(background_response) == ng)
         CPASSERT(SIZE(poisson_weights_response) == ng)
      END IF
      IF (PRESENT(map)) CALL check_map(map, centers, cell, gvecs, volume, repeat)
      ALLOCATE (work(SIZE(centers)), potentials(SIZE(centers)), h(ng), e(ng), dh(ng), dg(ng))
      h(:) = background
      e(:) = CMPLX(0.0_dp, 0.0_dp, dp)
      dh(:) = CMPLX(0.0_dp, 0.0_dp, dp)
      dg(:) = CMPLX(0.0_dp, 0.0_dp, dp)
      IF (derive) dh(:) = background_response
      DO a = 1, SIZE(centers)
         nr = SIZE(centers(a)%radii)
         nh = SIZE(centers(a)%l)
         ALLOCATE (work(a)%compensation(nr, nh), work(a)%residual(nr, nh), &
                   work(a)%q(nr, nh), work(a)%dq(nr, nh), work(a)%potential(nr, nh), &
                   work(a)%fourier(ng, nr, nh), moments(nh))
         CALL upf_inverse_coulomb_compensate(centers(a)%radii, centers(a)%weights, centers(a)%l, &
                                             centers(a)%density, centers(a)%shapes, work(a)%compensation, &
                                             work(a)%residual, moments)
         IF (PRESENT(map)) THEN
            work(a)%fourier(:, :, :) = map%fourier(a)%values
         ELSE
            CALL upf_inverse_coulomb_fourier_kernel(centers(a)%radii, centers(a)%l, centers(a)%indices, &
                                                    gvecs, centers(a)%position, volume, work(a)%fourier)
         END IF
         work(a)%q(:, :) = SPREAD(centers(a)%weights, 2, nh)*work(a)%residual
         work(a)%potential(:, :) = 0.0_dp
         work(a)%dq(:, :) = 0.0_dp
         IF (derive) THEN
            CALL check_radial(centers(a)%radii, centers(a)%weights, centers(a)%l, center_response(a)%density)
            CPASSERT(SIZE(center_response(a)%radii) == nr)
            CPASSERT(SIZE(center_response(a)%weights) == nr)
            CPASSERT(ALL(SHAPE(center_response(a)%shapes) == [nr, nh]))
            ALLOCATE (dc(nr, nh), de(nr, nh), f(nr), df(nr), dfourier(ng, nr, nh))
            DO channel = 1, nh
               f(:) = centers(a)%radii**centers(a)%l(channel)
               df(:) = f*centers(a)%l(channel)*center_response(a)%radii/centers(a)%radii
               norm = SUM(centers(a)%weights*f*centers(a)%shapes(:, channel))
               dnorm = SUM((center_response(a)%weights*f + centers(a)%weights*df)*centers(a)%shapes(:, channel) + &
                           centers(a)%weights*f*center_response(a)%shapes(:, channel))
               moment = moments(channel)
               dmoment = SUM((center_response(a)%weights*f + centers(a)%weights*df)*centers(a)%density(:, channel) + &
                             centers(a)%weights*f*center_response(a)%density(:, channel))
               dc(:, channel) = centers(a)%shapes(:, channel)*dmoment/norm + moment/norm* &
                                (center_response(a)%shapes(:, channel) - centers(a)%shapes(:, channel)*dnorm/norm)
            END DO
            de(:, :) = center_response(a)%density - dc
            work(a)%dq(:, :) = SPREAD(center_response(a)%weights, 2, nh)*work(a)%residual + &
                               SPREAD(centers(a)%weights, 2, nh)*de
            CALL upf_inverse_coulomb_fourier_kernel(centers(a)%radii, centers(a)%l, centers(a)%indices, &
                                                    gvecs, centers(a)%position, volume, work(a)%fourier, dfourier, &
                                                    center_response(a)%radii, gvecs_response, &
                                                    center_response(a)%position, volume_response)
         END IF
         DO channel = 1, nh
            DO i = 1, nr
               h(:) = h + work(a)%fourier(:, i, channel)*centers(a)%weights(i)*work(a)%compensation(i, channel)
               e(:) = e + work(a)%fourier(:, i, channel)*work(a)%q(i, channel)
               IF (.NOT. derive) CYCLE
               dh(:) = dh + dfourier(:, i, channel)*centers(a)%weights(i)*work(a)%compensation(i, channel) + &
                       work(a)%fourier(:, i, channel)*(center_response(a)%weights(i)*work(a)%compensation(i, channel) + &
                                                       centers(a)%weights(i)*dc(i, channel))
               dg(:) = dg + dfourier(:, i, channel)*work(a)%q(i, channel) + &
                       work(a)%fourier(:, i, channel)*work(a)%dq(i, channel)
            END DO
         END DO
         DEALLOCATE (moments)
         IF (derive) DEALLOCATE (dc, de, f, df, dfourier)
      END DO
      energy = SUM(poisson_weights*(0.5_dp*ABS(h)**2 + REAL(CONJG(h)*e, dp)))
      dpw = 0.0_dp
      IF (derive) THEN
         dpw = SUM(poisson_weights_response*(0.5_dp*ABS(h)**2 + REAL(CONJG(h)*e, dp)) + &
                   poisson_weights*REAL(CONJG(dh)*(h + e) + CONJG(h)*dg, dp))
      END IF
      IF (PRESENT(comm)) THEN
         CALL comm%sum(energy)
         IF (derive) CALL comm%sum(dpw)
      END IF
      background_adjoint = poisson_weights*(h + e)
      epair = 0.0_dp
      dpair = 0.0_dp
      DO a = 1, SIZE(centers)
         na = SIZE(centers(a)%radii)
         DO b = 1, SIZE(centers)
            IF (PRESENT(map)) THEN
               fixed_pair = .NOT. derive
               IF (derive) THEN
                  fixed_pair = ALL(cell_response == 0.0_dp) .AND. &
                               ALL(center_response(a)%position == center_response(b)%position) .AND. &
                               ALL(center_response(a)%radii == 0.0_dp) .AND. ALL(center_response(b)%radii == 0.0_dp)
               END IF
               IF (fixed_pair) THEN
                  DO channel = 1, SIZE(centers(a)%l)
                     DO nr = 1, na
                        norm = SUM(map%pairs(a, b)%values(nr, channel, :, :)*work(b)%q)
                        work(a)%potential(nr, channel) = work(a)%potential(nr, channel) + norm
                        epair = epair + 0.5_dp*work(a)%q(nr, channel)*norm
                        IF (.NOT. derive) CYCLE
                        dpair = dpair + 0.5_dp*(work(a)%dq(nr, channel)*norm + work(a)%q(nr, channel)* &
                                                SUM(map%pairs(a, b)%values(nr, channel, :, :)*work(b)%dq))
                     END DO
                  END DO
                  CYCLE
               END IF
            END IF
            nb = SIZE(centers(b)%radii)
            radius = MAXVAL(centers(a)%radii) + MAXVAL(centers(b)%radii)
            CALL residual_image_limits(inverse_cell, centers(a)%position - centers(b)%position, &
                                       radius, repeat, limits)
            ALLOCATE (kernel(na, SIZE(centers(a)%l), nb, SIZE(centers(b)%l)))
            IF (derive) ALLOCATE (dkernel, MOLD=kernel)
            DO j = limits(1, 3), limits(2, 3)
               DO i = limits(1, 2), limits(2, 2)
                  DO info = limits(1, 1), limits(2, 1)
                     image = [info, i, j]
                     displacement = centers(a)%position - centers(b)%position - MATMUL(cell, REAL(image, dp))
                     IF (SUM(displacement**2) >= radius**2) CYCLE
                     IF (derive) THEN
                        ddisplacement = center_response(a)%position - center_response(b)%position - &
                                        MATMUL(cell_response, REAL(image, dp))
                        CALL upf_inverse_coulomb_pair_kernel(centers(a)%radii, centers(a)%l, centers(a)%indices, &
                                                             centers(b)%radii, centers(b)%l, centers(b)%indices, &
                                                             displacement, kernel, dkernel, center_response(a)%radii, &
                                                             center_response(b)%radii, ddisplacement)
                     ELSE
                        CALL upf_inverse_coulomb_pair_kernel(centers(a)%radii, centers(a)%l, centers(a)%indices, &
                                                             centers(b)%radii, centers(b)%l, centers(b)%indices, &
                                                             displacement, kernel)
                     END IF
                     DO channel = 1, SIZE(centers(a)%l)
                        DO nr = 1, na
                           norm = SUM(kernel(nr, channel, :, :)*work(b)%q)
                           work(a)%potential(nr, channel) = work(a)%potential(nr, channel) + norm
                           epair = epair + 0.5_dp*work(a)%q(nr, channel)*norm
                           IF (.NOT. derive) CYCLE
                           dpair = dpair + 0.5_dp*(work(a)%dq(nr, channel)*norm + work(a)%q(nr, channel)* &
                                                   SUM(kernel(nr, channel, :, :)*work(b)%dq + &
                                                       dkernel(nr, channel, :, :)*work(b)%q))
                        END DO
                     END DO
                  END DO
               END DO
            END DO
            DEALLOCATE (kernel)
            IF (derive) DEALLOCATE (dkernel)
         END DO
      END DO
      energy = energy + epair
      IF (derive) response = dpw + dpair
      DO a = 1, SIZE(centers)
         nr = SIZE(centers(a)%radii)
         nh = SIZE(centers(a)%l)
         ALLOCATE (potentials(a)%values(nr, nh), vh(nr, nh), ve(nr, nh), adjoint(nr, nh))
         DO channel = 1, nh
            DO i = 1, nr
               vh(i, channel) = REAL(SUM(CONJG(work(a)%fourier(:, i, channel))*poisson_weights*h), dp)
               ve(i, channel) = REAL(SUM(CONJG(work(a)%fourier(:, i, channel))*background_adjoint), dp)
            END DO
         END DO
         IF (PRESENT(comm)) THEN
            CALL comm%sum(vh)
            CALL comm%sum(ve)
         END IF
         CALL upf_inverse_coulomb_compensation_adjoint(centers(a)%radii, centers(a)%weights, centers(a)%l, &
                                                       centers(a)%shapes, vh + work(a)%potential, adjoint)
         potentials(a)%values(:, :) = adjoint + ve
         CALL upf_inverse_coulomb_compensation_adjoint(centers(a)%radii, centers(a)%weights, centers(a)%l, &
                                                       centers(a)%shapes, ve, adjoint)
         potentials(a)%values(:, :) = potentials(a)%values - adjoint
         DEALLOCATE (vh, ve, adjoint)
      END DO
   END SUBROUTINE upf_inverse_coulomb_periodic

! **************************************************************************************************
!> \brief Reject stale geometric maps rather than silently changing the Coulomb functional.
!> \param map ...
!> \param centers ...
!> \param cell ...
!> \param gvecs ...
!> \param volume ...
!> \param periodicity ...
! **************************************************************************************************
   SUBROUTINE check_map(map, centers, cell, gvecs, volume, periodicity)
      TYPE(upf_inverse_coulomb_map_type), INTENT(IN)     :: map
      TYPE(upf_inverse_coulomb_center_type), INTENT(IN)  :: centers(:)
      REAL(dp), INTENT(IN)                               :: cell(3, 3), gvecs(:, :), volume
      INTEGER, INTENT(IN)                                :: periodicity(3)

      INTEGER                                            :: a

      CPASSERT(ALLOCATED(map%centers) .AND. ALLOCATED(map%pairs))
      CPASSERT(SIZE(map%centers) == SIZE(centers))
      CPASSERT(ALL(map%cell == cell) .AND. map%volume == volume)
      CPASSERT(ALL(map%periodicity == periodicity))
      CPASSERT(ALL(SHAPE(map%gvecs) == SHAPE(gvecs)))
      CPASSERT(ALL(map%gvecs == gvecs))
      DO a = 1, SIZE(centers)
         CPASSERT(ALL(map%centers(a)%position == centers(a)%position))
         CPASSERT(SIZE(map%centers(a)%radii) == SIZE(centers(a)%radii))
         CPASSERT(SIZE(map%centers(a)%l) == SIZE(centers(a)%l))
         CPASSERT(ALL(map%centers(a)%radii == centers(a)%radii))
         CPASSERT(ALL(map%centers(a)%l == centers(a)%l))
         CPASSERT(ALL(map%centers(a)%indices == centers(a)%indices))
      END DO
   END SUBROUTINE check_map

! **************************************************************************************************
!> \brief Bound overlapping residual images, without repeating nonperiodic lattice directions.
!> \param inverse_cell ...
!> \param displacement ...
!> \param radius ...
!> \param periodicity ...
!> \param limits ...
! **************************************************************************************************
   SUBROUTINE residual_image_limits(inverse_cell, displacement, radius, periodicity, limits)
      REAL(dp), INTENT(IN)                               :: inverse_cell(3, 3), displacement(3), &
                                                            radius
      INTEGER, INTENT(IN)                                :: periodicity(3)
      INTEGER, INTENT(OUT)                               :: limits(2, 3)

      INTEGER                                            :: direction
      REAL(dp)                                           :: bound(3), fractional(3)

      fractional = MATMUL(inverse_cell, displacement)
      bound = radius*SQRT(SUM(inverse_cell**2, DIM=2))
      DO direction = 1, 3
         IF (periodicity(direction) == 1) THEN
            limits(1, direction) = CEILING(fractional(direction) - bound(direction))
            limits(2, direction) = FLOOR(fractional(direction) + bound(direction))
         ELSE
            limits(:, direction) = 0
         END IF
      END DO
   END SUBROUTINE residual_image_limits

! **************************************************************************************************
!> \brief Reciprocal Coulomb map between two translated harmonic radial-shell distributions.
!>        The kernel excludes both radial measures. Splitting the angular integral at the source
!>        shell avoids sampling a discontinuous radial derivative across an unsplit angular rule.
!>        The joint response differentiates the same finite rule, including its moving boundaries.
!> \param radii_a Positive target radii.
!> \param l_a Angular momenta of the target real harmonics.
!> \param indices_a Native harmonic indices of the target channels.
!> \param radii_b Positive source radii.
!> \param l_b Angular momenta of the source real harmonics.
!> \param indices_b Native harmonic indices of the source channels.
!> \param displacement Target center minus source center, including any lattice image.
!> \param kernel Coulomb map in (radius_a,harmonic_a,radius_b,harmonic_b) order.
!> \param response Optional joint geometry response.
!> \param dr_a Required with response.
!> \param dr_b Required with response.
!> \param ddisplacement Required with response.
! **************************************************************************************************
   SUBROUTINE upf_inverse_coulomb_pair_kernel(radii_a, l_a, indices_a, radii_b, l_b, indices_b, &
                                              displacement, kernel, response, dr_a, dr_b, ddisplacement)
      REAL(dp), INTENT(IN)                               :: radii_a(:)
      INTEGER, INTENT(IN)                                :: l_a(:), indices_a(:)
      REAL(dp), INTENT(IN)                               :: radii_b(:)
      INTEGER, INTENT(IN)                                :: l_b(:), indices_b(:)
      REAL(dp), INTENT(IN)                               :: displacement(3)
      REAL(dp), INTENT(OUT)                              :: kernel(:, :, :, :)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: response(:, :, :, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: dr_a(:), dr_b(:), ddisplacement(3)

      INTEGER                                            :: a, b, i, j, nq
      LOGICAL                                            :: separated
      REAL(dp)                                           :: da, db, delta(3), distance, factor, value
      REAL(dp), ALLOCATABLE                              :: ab(:, :), ba(:, :), dab(:, :), &
                                                            dba(:, :), dpoint(:, :), nodes(:), &
                                                            point(:, :), weights(:)

      CPASSERT(ALL(radii_a > 0.0_dp) .AND. ALL(ieee_is_finite(radii_a)))
      CPASSERT(ALL(radii_b > 0.0_dp) .AND. ALL(ieee_is_finite(radii_b)))
      CPASSERT(ALL(ieee_is_finite(displacement)))
      CPASSERT(SIZE(l_a) > 0 .AND. SIZE(l_b) > 0)
      CPASSERT(SIZE(indices_a) == SIZE(l_a))
      CPASSERT(SIZE(indices_b) == SIZE(l_b))
      CPASSERT(ALL(l_a >= 0 .AND. l_a <= 20))
      CPASSERT(ALL(l_b >= 0 .AND. l_b <= 20))
      CPASSERT(ALL(indices_a > l_a**2 .AND. indices_a <= (l_a + 1)**2))
      CPASSERT(ALL(indices_b > l_b**2 .AND. indices_b <= (l_b + 1)**2))
      IF (.NOT. ALL(SHAPE(kernel) == [SIZE(radii_a), SIZE(l_a), SIZE(radii_b), SIZE(l_b)])) THEN
         CPABORT("Invalid translated shell kernel dimensions")
      END IF
      IF (PRESENT(response)) THEN
         CPASSERT(PRESENT(dr_a) .AND. PRESENT(dr_b))
         CPASSERT(PRESENT(ddisplacement))
         CPASSERT(SIZE(dr_a) == SIZE(radii_a))
         CPASSERT(SIZE(dr_b) == SIZE(radii_b))
         CPASSERT(ALL(SHAPE(response) == SHAPE(kernel)))
         response = 0.0_dp
      END IF
      kernel = 0.0_dp
      distance = SQRT(SUM(displacement**2))
      delta = 0.0_dp
      IF (PRESENT(response)) delta = ddisplacement
      IF (distance == 0.0_dp) THEN
         ! Same-center shells share their center motion. Distinct coincident centers need a limit.
         CPASSERT(ALL(delta == 0.0_dp))
         DO a = 1, SIZE(l_a)
            DO b = 1, SIZE(l_b)
               IF (indices_a(a) /= indices_b(b)) CYCLE
               DO i = 1, SIZE(radii_a)
                  DO j = 1, SIZE(radii_b)
                     value = fourpi/REAL(2*l_a(a) + 1, dp)/MAX(radii_a(i), radii_b(j))* &
                             (MIN(radii_a(i), radii_b(j))/MAX(radii_a(i), radii_b(j)))**l_a(a)
                     kernel(i, a, j, b) = value
                     IF (.NOT. PRESENT(response)) CYCLE
                     IF (radii_a(i) == radii_b(j)) THEN
                        CPASSERT(dr_a(i) == dr_b(j))
                     END IF
                     IF (radii_a(i) >= radii_b(j)) THEN
                        response(i, a, j, b) = value*(l_a(a)*dr_b(j)/radii_b(j) - &
                                                      (l_a(a) + 1)*dr_a(i)/radii_a(i))
                     ELSE
                        response(i, a, j, b) = value*(l_a(a)*dr_a(i)/radii_a(i) - &
                                                      (l_a(a) + 1)*dr_b(j)/radii_b(j))
                     END IF
                  END DO
               END DO
            END DO
         END DO
         RETURN
      END IF
      CALL init_orbital_pointers(MAX(MAXVAL(l_a), MAXVAL(l_b)))
      CALL init_spherical_harmonics(MAX(MAXVAL(l_a), MAXVAL(l_b)), -1)
      ALLOCATE (ab(SIZE(l_a), SIZE(l_b)), ba(SIZE(l_b), SIZE(l_a)))
      IF (PRESENT(response)) ALLOCATE (dab(SIZE(l_a), SIZE(l_b)), dba(SIZE(l_b), SIZE(l_a)))
      nq = MAX(32, MAXVAL(l_a) + MAXVAL(l_b) + 2)
      ALLOCATE (nodes(nq), weights(nq))
      CALL polar_rule(nodes, weights)
      ! Separated harmonic shells have an exact point-multipole interaction.
      ! Higher orders retain the general angular rule, without discarding any channel.
      separated = MAXVAL(l_a) + MAXVAL(l_b) <= 8 .AND. distance > MINVAL(radii_a) + MINVAL(radii_b)
      IF (separated) THEN
         ALLOCATE (point(SIZE(l_a), SIZE(l_b)))
         IF (PRESENT(response)) THEN
            ALLOCATE (dpoint, MOLD=point)
            CALL separated_shell_map(l_a, indices_a, l_b, indices_b, displacement, point, dpoint, delta)
         ELSE
            CALL separated_shell_map(l_a, indices_a, l_b, indices_b, displacement, point)
         END IF
      END IF
      DO i = 1, SIZE(radii_a)
         DO j = 1, SIZE(radii_b)
            da = 0.0_dp
            db = 0.0_dp
            IF (PRESENT(response)) THEN
               da = dr_a(i)
               db = dr_b(j)
            END IF
            IF (separated .AND. distance > radii_a(i) + radii_b(j)) THEN
               DO b = 1, SIZE(l_b)
                  DO a = 1, SIZE(l_a)
                     factor = radii_a(i)**l_a(a)*radii_b(j)**l_b(b)
                     kernel(i, a, j, b) = factor*point(a, b)
                     IF (PRESENT(response)) response(i, a, j, b) = factor*(dpoint(a, b) + point(a, b)* &
                                                                           (l_a(a)*da/radii_a(i) + &
                                                                            l_b(b)*db/radii_b(j)))
                  END DO
               END DO
               CYCLE
            END IF
            IF (PRESENT(response)) THEN
               CALL translated_shell(radii_a(i), radii_b(j), l_a, indices_a, l_b, indices_b, &
                                     displacement, da, db, delta, ab, nodes, weights, dab)
               CALL translated_shell(radii_b(j), radii_a(i), l_b, indices_b, l_a, indices_a, &
                                     -displacement, db, da, -delta, ba, nodes, weights, dba)
            ELSE
               CALL translated_shell(radii_a(i), radii_b(j), l_a, indices_a, l_b, indices_b, &
                                     displacement, da, db, delta, ab, nodes, weights)
               CALL translated_shell(radii_b(j), radii_a(i), l_b, indices_b, l_a, indices_a, &
                                     -displacement, db, da, -delta, ba, nodes, weights)
            END IF
            kernel(i, :, j, :) = 0.5_dp*(ab + TRANSPOSE(ba))
            IF (PRESENT(response)) response(i, :, j, :) = 0.5_dp*(dab + TRANSPOSE(dba))
         END DO
      END DO
   END SUBROUTINE upf_inverse_coulomb_pair_kernel

! **************************************************************************************************
!> \brief Point-multipole map for disjoint shells in the native real solid-harmonic convention.
!>        K_ab = (4*pi)^2*(-1)^lb*Ra(partial)*Rb(partial)(1/d)/[(2*la+1)!!*(2*lb+1)!!].
!>        Shell radii enter separately as a**la*b**lb. Cartesian Coulomb derivatives are exact.
!> \param la ...
!> \param ia ...
!> \param lb ...
!> \param ib ...
!> \param d Nonzero target-minus-source displacement.
!> \param values Kernel without the two shell-radius powers.
!> \param response Optional displacement response of the same kernel.
!> \param dd Required with response.
! **************************************************************************************************
   SUBROUTINE separated_shell_map(la, ia, lb, ib, d, values, response, dd)
      INTEGER, INTENT(IN)                                :: la(:), ia(:), lb(:), ib(:)
      REAL(dp), INTENT(IN)                               :: d(3)
      REAL(dp), INTENT(OUT)                              :: values(:, :)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: response(:, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: dd(3)

      INTEGER                                            :: a, axis, b, bx, by, bz, contractions, i, &
                                                            ix, iy, iz, j, n, order, powers(3), &
                                                            reduced(3), shifted(3)
      REAL(dp)                                           :: coefficient, direction(3), distance, &
                                                            factor, term
      REAL(dp), ALLOCATABLE                              :: factorial(:), odd_factorial(:), &
                                                            tensor(:, :, :)

      CPASSERT(PRESENT(response) .EQV. PRESENT(dd))
      n = MAXVAL(la) + MAXVAL(lb)
      IF (PRESENT(response)) n = n + 1
      distance = SQRT(SUM(d**2))
      CPASSERT(distance > 0.0_dp)
      direction = d/distance
      ALLOCATE (factorial(0:n + 1), odd_factorial(0:n + 1), tensor(0:n, 0:n, 0:n))
      factorial(0) = 1.0_dp
      odd_factorial(0) = 1.0_dp
      DO i = 1, n + 1
         factorial(i) = i*factorial(i - 1)
         odd_factorial(i) = (2*i - 1)*odd_factorial(i - 1)
      END DO
      tensor = 0.0_dp
      DO iz = 0, n
         DO iy = 0, n - iz
            DO ix = 0, n - iz - iy
               powers = [ix, iy, iz]
               order = SUM(powers)
               factor = PRODUCT(factorial(powers))/distance**(order + 1)
               DO bz = 0, iz/2
                  DO by = 0, iy/2
                     DO bx = 0, ix/2
                        contractions = bx + by + bz
                        reduced = powers - 2*[bx, by, bz]
                        term = REAL((-1)**(order - contractions), dp)*odd_factorial(order - contractions)* &
                               PRODUCT(direction**reduced)/(PRODUCT(factorial([bx, by, bz]))* &
                                                            PRODUCT(factorial(reduced))*2.0_dp**contractions)
                        tensor(ix, iy, iz) = tensor(ix, iy, iz) + factor*term
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
      values = 0.0_dp
      IF (PRESENT(response)) response = 0.0_dp
      DO b = 1, SIZE(lb)
         DO a = 1, SIZE(la)
            factor = fourpi**2*REAL((-1)**lb(b), dp)/(odd_factorial(la(a) + 1)*odd_factorial(lb(b) + 1))
            DO j = 1, nco(lb(b))
               DO i = 1, nco(la(a))
                  coefficient = factor*orbtramat(la(a))%slm(ia(a) - la(a)**2, i)* &
                                orbtramat(lb(b))%slm(ib(b) - lb(b)**2, j)
                  IF (coefficient == 0.0_dp) CYCLE
                  powers = indco(:, ncoset(la(a) - 1) + i) + indco(:, ncoset(lb(b) - 1) + j)
                  values(a, b) = values(a, b) + coefficient*tensor(powers(1), powers(2), powers(3))
                  IF (.NOT. PRESENT(response)) CYCLE
                  DO axis = 1, 3
                     shifted = powers
                     shifted(axis) = shifted(axis) + 1
                     response(a, b) = response(a, b) + &
                                      coefficient*dd(axis)*tensor(shifted(1), shifted(2), shifted(3))
                  END DO
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE separated_shell_map

! **************************************************************************************************
!> \brief Integrate a source shell's regular/irregular solid-harmonic potential on a target shell.
!>        Azimuthal integration is exact for the represented polynomial degree. Polar subintervals
!>        also limit the source-distance ratio to two in the irregular branch.
!> \param a ...
!> \param b ...
!> \param la ...
!> \param ia ...
!> \param lb ...
!> \param ib ...
!> \param d ...
!> \param da ...
!> \param db ...
!> \param dd ...
!> \param value ...
!> \param nodes Shared polar quadrature nodes, independent of the two shell radii.
!> \param weights Shared polar quadrature weights.
!> \param response ...
! **************************************************************************************************
   SUBROUTINE translated_shell(a, b, la, ia, lb, ib, d, da, db, dd, value, nodes, weights, response)
      REAL(dp), INTENT(IN)                               :: a, b
      INTEGER, INTENT(IN)                                :: la(:), ia(:), lb(:), ib(:)
      REAL(dp), INTENT(IN)                               :: d(3), da, db, dd(3)
      REAL(dp), INTENT(OUT)                              :: value(:, :)
      REAL(dp), INTENT(IN)                               :: nodes(:), weights(:)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: response(:, :)

      INTEGER                                            :: axis, channel, iphi, iq, nphi, nq
      LOGICAL                                            :: derive, regular
      REAL(dp) :: azimuth, bound, dbound, dex(3), dey(3), dez(3), dfactor, distance, dleft, &
         dlength, dmu, dright, dt, du(3), dweight, dx(3), ex(3), ey(3), ez(3), factor, left, &
         length, mu, norm_x, reference(3), right, scale, sine, t, tmax, u(3), weight, x(3)
      REAL(dp), ALLOCATABLE                              :: dya(:), dyb(:), ga(:, :), gb(:, :), &
                                                            ya(:), yb(:)

      derive = PRESENT(response)
      distance = SQRT(SUM(d**2))
      ez = d/distance
      reference = 0.0_dp
      axis = MINLOC(ABS(ez), DIM=1)
      reference(axis) = 1.0_dp
      ex = reference - ez(axis)*ez
      length = SQRT(SUM(ex**2))
      IF (derive) THEN
         dt = DOT_PRODUCT(d, dd)/distance
         dez = (dd - ez*dt)/distance
         dex = -dez(axis)*ez - ez(axis)*dez
         dlength = DOT_PRODUCT(ex, dex)/length
         dex = (dex - ex*dlength/length)/length
      END IF
      ex = ex/length
      ey = cross(ez, ex)
      IF (derive) dey = cross(dez, ex) + cross(ez, dex)
      nq = SIZE(nodes)
      nphi = MAXVAL(la) + MAXVAL(lb) + 1
      ALLOCATE (ya(SIZE(la)), yb(SIZE(lb)))
      IF (derive) ALLOCATE (dya(SIZE(la)), dyb(SIZE(lb)), ga(3, SIZE(la)), gb(3, SIZE(lb)))
      value = 0.0_dp
      IF (derive) response = 0.0_dp
      tmax = distance + a
      left = -1.0_dp
      dleft = 0.0_dp
      regular = ABS(distance - a) < b
      DO WHILE (left < 1.0_dp)
         right = 1.0_dp
         dright = 0.0_dp
         t = SQRT(MAX(0.0_dp, (distance - a)**2 + 2.0_dp*distance*a*(left + 1.0_dp)))
         IF (regular .AND. b < tmax) THEN
            bound = b
            IF (derive) dbound = db
         ELSE IF (.NOT. regular .AND. 2.0_dp*t < tmax) THEN
            bound = 2.0_dp*t
            IF (derive) dbound = 2.0_dp*((distance - a)*(dt - da) + &
                                         (dt*a + distance*da)*(left + 1.0_dp) + distance*a*dleft)/t
         ELSE
            bound = tmax
            IF (derive) dbound = dt + da
         END IF
         IF (bound < tmax) THEN
            right = ((bound - a)*(bound + a) - distance**2)/(2.0_dp*distance*a)
            IF (derive) dright = (bound*dbound - a*da - distance*dt)/(distance*a) - &
                                 right*(dt/distance + da/a)
         END IF
         ! At shell tangency, roundoff can leave a cap too narrow to place interior nodes.
         ! The two potentials agree on the shell, so its moving boundary terms cancel.
         IF (right - left <= 64.0_dp*EPSILON(1.0_dp)) THEN
            CPASSERT(regular)
            regular = .FALSE.
            CYCLE
         END IF
         IF (1.0_dp - right <= 64.0_dp*EPSILON(1.0_dp)) THEN
            right = 1.0_dp
            dright = 0.0_dp
         END IF
         CPASSERT(right > left .AND. right <= 1.0_dp)
         DO iq = 1, nq
            mu = 0.5_dp*((1.0_dp - nodes(iq))*left + (1.0_dp + nodes(iq))*right)
            IF (derive) dmu = 0.5_dp*((1.0_dp - nodes(iq))*dleft + (1.0_dp + nodes(iq))*dright)
            weight = 0.5_dp*(right - left)*weights(iq)*twopi/REAL(nphi, dp)
            IF (derive) dweight = 0.5_dp*(dright - dleft)*weights(iq)*twopi/REAL(nphi, dp)
            sine = SQRT(MAX(0.0_dp, 1.0_dp - mu**2))
            CPASSERT(sine > 0.0_dp)
            DO iphi = 1, nphi
               azimuth = twopi*REAL(iphi - 1, dp)/REAL(nphi, dp)
               u = mu*ez + sine*(COS(azimuth)*ex + SIN(azimuth)*ey)
               x = d + a*u
               norm_x = SQRT(SUM(x**2))
               IF (derive) THEN
                  du = dmu*ez + mu*dez - mu*dmu/sine*(COS(azimuth)*ex + SIN(azimuth)*ey) + &
                       sine*(COS(azimuth)*dex + SIN(azimuth)*dey)
                  dx = dd + da*u + a*du
                  CALL solid_fields(u, la, ia, ya, ga)
                  CALL solid_fields(x, lb, ib, yb, gb)
                  DO channel = 1, SIZE(la)
                     dya(channel) = DOT_PRODUCT(du, ga(:, channel))
                  END DO
                  DO channel = 1, SIZE(lb)
                     dyb(channel) = DOT_PRODUCT(dx, gb(:, channel))
                  END DO
               ELSE
                  CALL solid_fields(u, la, ia, ya)
                  CALL solid_fields(x, lb, ib, yb)
               END IF
               DO channel = 1, SIZE(lb)
                  IF (norm_x < b) THEN
                     factor = fourpi/REAL(2*lb(channel) + 1, dp)/b**(lb(channel) + 1)
                     IF (derive) dfactor = -factor*(lb(channel) + 1)*db/b
                  ELSE
                     factor = fourpi/REAL(2*lb(channel) + 1, dp)* &
                              (b/norm_x)**lb(channel)/norm_x**(lb(channel) + 1)
                     IF (derive) dfactor = factor*(lb(channel)*db/b - &
                                                   (2*lb(channel) + 1)*DOT_PRODUCT(x, dx)/norm_x**2)
                  END IF
                  scale = factor*yb(channel)
                  value(:, channel) = value(:, channel) + weight*ya*scale
                  IF (derive) response(:, channel) = response(:, channel) + dweight*ya*scale + &
                                                     weight*(dya*scale + ya*(dfactor*yb(channel) + factor*dyb(channel)))
               END DO
            END DO
         END DO
         left = right
         dleft = dright
         regular = .FALSE.
      END DO
   END SUBROUTINE translated_shell

! **************************************************************************************************
!> \brief Native regular solid harmonics and Cartesian gradients, including their origin limits.
!> \param x ...
!> \param ls ...
!> \param indices ...
!> \param values ...
!> \param gradients ...
! **************************************************************************************************
   SUBROUTINE solid_fields(x, ls, indices, values, gradients)
      REAL(dp), INTENT(IN)                               :: x(3)
      INTEGER, INTENT(IN)                                :: ls(:), indices(:)
      REAL(dp), INTENT(OUT)                              :: values(:)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: gradients(:, :)

      INTEGER                                            :: axis, channel, i, l, m, power, &
                                                            powers(3), reduced(3)
      REAL(dp)                                           :: factor, xpowers(0:20, 3)

      values = 0.0_dp
      IF (PRESENT(gradients)) gradients = 0.0_dp
      xpowers(0, :) = 1.0_dp
      DO axis = 1, 3
         DO power = 1, MAXVAL(ls)
            xpowers(power, axis) = x(axis)*xpowers(power - 1, axis)
         END DO
      END DO
      DO channel = 1, SIZE(ls)
         l = ls(channel)
         m = indices(channel) - l*l
         DO i = 1, nco(l)
            powers = indco(:, ncoset(l - 1) + i)
            factor = orbtramat(l)%slm(m, i)
            IF (factor == 0.0_dp) CYCLE
            values(channel) = values(channel) + factor* &
                              xpowers(powers(1), 1)*xpowers(powers(2), 2)*xpowers(powers(3), 3)
            IF (.NOT. PRESENT(gradients)) CYCLE
            DO axis = 1, 3
               IF (powers(axis) == 0) CYCLE
               reduced = powers
               reduced(axis) = reduced(axis) - 1
               gradients(axis, channel) = gradients(axis, channel) + &
                                          factor*powers(axis)*xpowers(reduced(1), 1)* &
                                          xpowers(reduced(2), 2)*xpowers(reduced(3), 3)
            END DO
         END DO
      END DO
   END SUBROUTINE solid_fields

! **************************************************************************************************
!> \brief Gauss-Legendre nodes and weights on (-1,1).
!> \param nodes ...
!> \param weights ...
! **************************************************************************************************
   SUBROUTINE polar_rule(nodes, weights)
      REAL(dp), INTENT(OUT)                              :: nodes(:), weights(:)

      INTEGER                                            :: i, iteration, j, n
      REAL(dp)                                           :: derivative, p, previous, saved, shift, x

      n = SIZE(nodes)
      CPASSERT(SIZE(weights) == n)
      DO i = 1, (n + 1)/2
         x = COS(pi*(REAL(i, dp) - 0.25_dp)/(REAL(n, dp) + 0.5_dp))
         DO iteration = 1, 20
            p = 1.0_dp
            previous = 0.0_dp
            DO j = 1, n
               saved = p
               p = ((2*j - 1)*x*p - (j - 1)*previous)/REAL(j, dp)
               previous = saved
            END DO
            derivative = n*(x*p - previous)/(x*x - 1.0_dp)
            shift = p/derivative
            x = x - shift
            IF (ABS(shift) < 2.0_dp*EPSILON(x)) EXIT
         END DO
         CPASSERT(iteration <= 20)
         nodes(i) = -x
         nodes(n + 1 - i) = x
         weights(i) = 2.0_dp/((1.0_dp - x*x)*derivative**2)
         weights(n + 1 - i) = weights(i)
      END DO
   END SUBROUTINE polar_rule

! **************************************************************************************************
!> \brief Cartesian vector product.
!> \param a ...
!> \param b ...
!> \return ...
! **************************************************************************************************
   FUNCTION cross(a, b) RESULT(c)
      REAL(dp), INTENT(IN)                               :: a(3), b(3)
      REAL(dp)                                           :: c(3)

      c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
   END FUNCTION cross

! **************************************************************************************************
!> \brief Fourier map of the same radial shells, for coupling to a periodic background.
!>        rho_G = sum_ir,lm weights(ir)*kernel(G,ir,lm)*density(ir,lm).
!>        The reciprocal covector must use the caller's actual full/half-grid multiplicities.
!>        No reciprocal Coulomb cutoff or G=0 electrostatic convention is selected here.
!> \param radii Positive radial shell nodes.
!> \param angular_momenta l of each orthonormal real harmonic channel.
!> \param harmonic_indices Native harmonic indices, between l*l+1 and (l+1)**2.
!> \param gvecs Cartesian reciprocal vectors, including zero if required by the caller.
!> \param position Cartesian center of these radial shells.
!> \param volume Physical cell volume.
!> \param kernel Complex Fourier coefficients per unit radial measure and density.
!> \param response Optional joint geometry response of the same Fourier map.
!> \param radii_response Required with response.
!> \param gvecs_response Required with response, including reciprocal-vector strain.
!> \param position_response Required with response.
!> \param volume_response Required with response.
! **************************************************************************************************
   SUBROUTINE upf_inverse_coulomb_fourier_kernel(radii, angular_momenta, harmonic_indices, gvecs, position, volume, &
                                                 kernel, response, radii_response, gvecs_response, &
                                                 position_response, volume_response)
      REAL(dp), INTENT(IN)                               :: radii(:)
      INTEGER, INTENT(IN)                                :: angular_momenta(:), harmonic_indices(:)
      REAL(dp), INTENT(IN)                               :: gvecs(:, :), position(3), volume
      COMPLEX(dp), INTENT(OUT)                           :: kernel(:, :, :)
      COMPLEX(dp), INTENT(OUT), OPTIONAL                 :: response(:, :, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: radii_response(:), gvecs_response(:, :), &
                                                            position_response(3), volume_response

      COMPLEX(dp)                                        :: phase, phase_l, phase_response
      INTEGER                                            :: a, axis, channel, ig, ir, l, m, &
                                                            powers(3), reduced(3)
      REAL(dp)                                           :: angular, angular_gradient(3), &
                                                            direction(3), factor, gradient(3), jl, &
                                                            jnext, norm_g, radial

      CPASSERT(SIZE(harmonic_indices) == SIZE(angular_momenta))
      CPASSERT(SIZE(angular_momenta) > 0)
      CPASSERT(SIZE(gvecs, 1) == 3)
      IF (.NOT. ALL(SHAPE(kernel) == [SIZE(gvecs, 2), SIZE(radii), SIZE(angular_momenta)])) THEN
         CPABORT("Invalid radial Fourier kernel dimensions")
      END IF
      CPASSERT(ALL(ieee_is_finite(radii)) .AND. ALL(radii > 0.0_dp))
      CPASSERT(ALL(ieee_is_finite(gvecs)))
      CPASSERT(ALL(ieee_is_finite(position)))
      CPASSERT(ieee_is_finite(volume) .AND. volume > 0.0_dp)
      CPASSERT(ALL(angular_momenta >= 0 .AND. angular_momenta <= 20))
      CPASSERT(ALL(harmonic_indices > angular_momenta**2))
      CPASSERT(ALL(harmonic_indices <= (angular_momenta + 1)**2))
      IF (PRESENT(response)) THEN
         CPASSERT(ALL(SHAPE(response) == SHAPE(kernel)))
         CPASSERT(PRESENT(radii_response) .AND. PRESENT(gvecs_response))
         CPASSERT(PRESENT(position_response))
         CPASSERT(PRESENT(volume_response))
         CPASSERT(SIZE(radii_response) == SIZE(radii))
         CPASSERT(ALL(SHAPE(gvecs_response) == SHAPE(gvecs)))
      END IF
      CALL init_orbital_pointers(MAXVAL(angular_momenta))
      CALL init_spherical_harmonics(MAXVAL(angular_momenta), -1)
      DO ig = 1, SIZE(gvecs, 2)
         norm_g = SQRT(SUM(gvecs(:, ig)**2))
         direction = 0.0_dp
         IF (norm_g > 0.0_dp) direction = gvecs(:, ig)/norm_g
         phase = fourpi/volume*EXP(CMPLX(0.0_dp, -DOT_PRODUCT(gvecs(:, ig), position), dp))
         IF (PRESENT(response)) THEN
            phase_response = CMPLX(-volume_response/volume, &
                                   -DOT_PRODUCT(gvecs_response(:, ig), position) - &
                                   DOT_PRODUCT(gvecs(:, ig), position_response), dp)
         END IF
         DO channel = 1, SIZE(angular_momenta)
            l = angular_momenta(channel)
            m = harmonic_indices(channel) - l*l
            phase_l = phase*CMPLX(0.0_dp, -1.0_dp, dp)**l
            angular = 0.0_dp
            angular_gradient = 0.0_dp
            DO a = 1, nco(l)
               powers = indco(:, ncoset(l - 1) + a)
               factor = orbtramat(l)%slm(m, a)
               angular = angular + factor*PRODUCT(direction**powers)
               IF (.NOT. PRESENT(response)) CYCLE
               DO axis = 1, 3
                  IF (powers(axis) == 0) CYCLE
                  reduced = powers
                  reduced(axis) = reduced(axis) - 1
                  angular_gradient(axis) = angular_gradient(axis) + &
                                           factor*REAL(powers(axis), dp)*PRODUCT(direction**reduced)
               END DO
            END DO
            DO ir = 1, SIZE(radii)
               jl = bessel0(norm_g*radii(ir), l)
               kernel(ig, ir, channel) = phase_l*angular*jl
               IF (.NOT. PRESENT(response)) CYCLE
               gradient = 0.0_dp
               radial = 0.0_dp
               IF (norm_g > 0.0_dp) THEN
                  jnext = bessel0(norm_g*radii(ir), l + 1)
                  gradient = jl/norm_g*angular_gradient - radii(ir)*jnext*angular*direction
                  radial = angular*(l/radii(ir)*jl - norm_g*jnext)
               ELSE IF (l == 1) THEN
                  ! The origin limit contains dipole derivatives even though rho_G itself is zero.
                  gradient = radii(ir)/3.0_dp*angular_gradient
               END IF
               response(ig, ir, channel) = kernel(ig, ir, channel)*phase_response + &
                                           phase_l*(DOT_PRODUCT(gradient, gvecs_response(:, ig)) + &
                                                    radial*radii_response(ir))
            END DO
         END DO
      END DO
   END SUBROUTINE upf_inverse_coulomb_fourier_kernel

! **************************************************************************************************
!> \brief Apply the symmetric Coulomb kernel to real, orthonormal harmonic coefficients.
!> \param radii Positive radial nodes, in arbitrary order.
!> \param weights Radial measure including r**2 dr, without angular weights.
!> \param angular_momenta l of each real harmonic channel. Repeated l values are allowed.
!> \param density Signed radial density coefficients, indexed (radius,harmonic).
!> \param potential Unweighted Coulomb potential in the same harmonic order.
!> \param energy One half of the density/potential contraction with the same radial measure.
! **************************************************************************************************
   SUBROUTINE upf_inverse_coulomb_radial(radii, weights, angular_momenta, density, potential, energy)
      REAL(dp), INTENT(IN)                               :: radii(:), weights(:)
      INTEGER, INTENT(IN)                                :: angular_momenta(:)
      REAL(dp), INTENT(IN)                               :: density(:, :)
      REAL(dp), INTENT(OUT)                              :: potential(:, :), energy

      INTEGER                                            :: channel, i, j, l
      REAL(dp)                                           :: kernel

      CALL check_radial(radii, weights, angular_momenta, density)
      CPASSERT(ALL(SHAPE(potential) == SHAPE(density)))
      potential = 0.0_dp
      DO channel = 1, SIZE(angular_momenta)
         l = angular_momenta(channel)
         DO i = 1, SIZE(radii)
            DO j = 1, SIZE(radii)
               kernel = fourpi/REAL(2*l + 1, dp)/MAX(radii(i), radii(j))* &
                        (MIN(radii(i), radii(j))/MAX(radii(i), radii(j)))**l
               potential(i, channel) = potential(i, channel) + kernel*weights(j)*density(j, channel)
            END DO
         END DO
      END DO
      energy = 0.5_dp*SUM(SPREAD(weights, 2, SIZE(angular_momenta))*density*potential)
   END SUBROUTINE upf_inverse_coulomb_radial

! **************************************************************************************************
!> \brief Differentiate the same finite Coulomb energy, including the measure and kernel.
!>        Equal nodes must have equal radial motion so their common shell remains differentiable.
!> \param radii ...
!> \param weights ...
!> \param angular_momenta ...
!> \param density ...
!> \param radii_response Radial node motion. Zero for Cartesian grids moving with their centers.
!> \param weights_response Complete radial measure response.
!> \param density_response Signed density response, including the coupled inverse coefficients.
!> \param response Coulomb energy response.
! **************************************************************************************************
   SUBROUTINE upf_inverse_coulomb_radial_response(radii, weights, angular_momenta, density, &
                                                  radii_response, weights_response, density_response, response)
      REAL(dp), INTENT(IN)                               :: radii(:), weights(:)
      INTEGER, INTENT(IN)                                :: angular_momenta(:)
      REAL(dp), INTENT(IN)                               :: density(:, :), radii_response(:), &
                                                            weights_response(:), &
                                                            density_response(:, :)
      REAL(dp), INTENT(OUT)                              :: response

      INTEGER                                            :: channel, i, j, l, large, small
      REAL(dp)                                           :: derivative, energy, kernel
      REAL(dp), ALLOCATABLE                              :: potential(:, :)

      CPASSERT(SIZE(radii_response) == SIZE(radii))
      CPASSERT(SIZE(weights_response) == SIZE(weights))
      CPASSERT(ALL(SHAPE(density_response) == SHAPE(density)))
      ALLOCATE (potential(SIZE(density, 1), SIZE(density, 2)))
      CALL upf_inverse_coulomb_radial(radii, weights, angular_momenta, density, potential, energy)
      response = SUM((SPREAD(weights, 2, SIZE(angular_momenta))*density_response + &
                      SPREAD(weights_response, 2, SIZE(angular_momenta))*density)*potential)
      DO channel = 1, SIZE(angular_momenta)
         l = angular_momenta(channel)
         DO i = 1, SIZE(radii)
            DO j = 1, SIZE(radii)
               small = i
               large = j
               IF (radii(i) > radii(j)) THEN
                  small = j
                  large = i
               END IF
               IF (radii(i) == radii(j)) THEN
                  CPASSERT(radii_response(i) == radii_response(j))
               END IF
               kernel = fourpi/REAL(2*l + 1, dp)/radii(large)*(radii(small)/radii(large))**l
               derivative = kernel*(l*radii_response(small)/radii(small) - &
                                    (l + 1)*radii_response(large)/radii(large))
               response = response + 0.5_dp*weights(i)*weights(j)* &
                          density(i, channel)*density(j, channel)*derivative
            END DO
         END DO
      END DO
   END SUBROUTINE upf_inverse_coulomb_radial_response

! **************************************************************************************************
!> \brief Match each multipole with a supplied radial shape in the defining finite quadrature.
!>        The signed remainder has zero exterior moment in every represented harmonic channel.
!>        This identity applies to this radial representation, not to discarded angular channels.
!> \param radii ...
!> \param weights ...
!> \param angular_momenta ...
!> \param density ...
!> \param shapes Positive radial compensation shapes, indexed (radius,harmonic).
!> \param compensation Moment-matched compensation density in the same radial representation.
!> \param residual Signed density minus compensation.
!> \param moments Original exterior multipole moments.
! **************************************************************************************************
   SUBROUTINE upf_inverse_coulomb_compensate(radii, weights, angular_momenta, density, shapes, &
                                             compensation, residual, moments)
      REAL(dp), INTENT(IN)                               :: radii(:), weights(:)
      INTEGER, INTENT(IN)                                :: angular_momenta(:)
      REAL(dp), INTENT(IN)                               :: density(:, :), shapes(:, :)
      REAL(dp), INTENT(OUT)                              :: compensation(:, :), residual(:, :), &
                                                            moments(:)

      INTEGER                                            :: channel
      REAL(dp)                                           :: norm

      CALL check_radial(radii, weights, angular_momenta, density)
      CPASSERT(ALL(SHAPE(shapes) == SHAPE(density)))
      CPASSERT(ALL(SHAPE(compensation) == SHAPE(density)))
      CPASSERT(ALL(SHAPE(residual) == SHAPE(density)))
      CPASSERT(SIZE(moments) == SIZE(angular_momenta))
      CPASSERT(ALL(ieee_is_finite(shapes)) .AND. ALL(shapes >= 0.0_dp))
      DO channel = 1, SIZE(angular_momenta)
         norm = SUM(weights*radii**angular_momenta(channel)*shapes(:, channel))
         CPASSERT(norm > 0.0_dp)
         moments(channel) = SUM(weights*radii**angular_momenta(channel)*density(:, channel))
         compensation(:, channel) = moments(channel)/norm*shapes(:, channel)
      END DO
      residual = density - compensation
   END SUBROUTINE upf_inverse_coulomb_compensate

! **************************************************************************************************
!> \brief Apply the radial-measure adjoint of density minus its moment compensation.
!>        Compensation depends on the occupation density and contributes to the AO operator.
!> \param radii ...
!> \param weights ...
!> \param angular_momenta ...
!> \param shapes Same fixed radial compensation shapes as in the forward map.
!> \param potential Unweighted covector acting on the signed remainder.
!> \param adjoint Unweighted covector acting on the original density.
! **************************************************************************************************
   SUBROUTINE upf_inverse_coulomb_compensation_adjoint(radii, weights, angular_momenta, shapes, potential, adjoint)
      REAL(dp), INTENT(IN)                               :: radii(:), weights(:)
      INTEGER, INTENT(IN)                                :: angular_momenta(:)
      REAL(dp), INTENT(IN)                               :: shapes(:, :), potential(:, :)
      REAL(dp), INTENT(OUT)                              :: adjoint(:, :)

      INTEGER                                            :: channel
      REAL(dp)                                           :: norm

      CALL check_radial(radii, weights, angular_momenta, potential)
      CPASSERT(ALL(SHAPE(shapes) == SHAPE(potential)))
      CPASSERT(ALL(SHAPE(adjoint) == SHAPE(potential)))
      CPASSERT(ALL(ieee_is_finite(shapes)) .AND. ALL(shapes >= 0.0_dp))
      DO channel = 1, SIZE(angular_momenta)
         norm = SUM(weights*radii**angular_momenta(channel)*shapes(:, channel))
         CPASSERT(norm > 0.0_dp)
         adjoint(:, channel) = potential(:, channel) - radii**angular_momenta(channel)* &
                               SUM(weights*shapes(:, channel)*potential(:, channel))/norm
      END DO
   END SUBROUTINE upf_inverse_coulomb_compensation_adjoint

! **************************************************************************************************
!> \brief Validate the common real radial representation.
!> \param radii ...
!> \param weights ...
!> \param angular_momenta ...
!> \param density ...
! **************************************************************************************************
   SUBROUTINE check_radial(radii, weights, angular_momenta, density)
      REAL(dp), INTENT(IN)                               :: radii(:), weights(:)
      INTEGER, INTENT(IN)                                :: angular_momenta(:)
      REAL(dp), INTENT(IN)                               :: density(:, :)

      CPASSERT(SIZE(weights) == SIZE(radii))
      CPASSERT(SIZE(density, 1) == SIZE(radii))
      CPASSERT(SIZE(density, 2) == SIZE(angular_momenta))
      CPASSERT(ALL(ieee_is_finite(radii)) .AND. ALL(radii > 0.0_dp))
      CPASSERT(ALL(ieee_is_finite(weights)) .AND. ALL(weights >= 0.0_dp))
      CPASSERT(ALL(ieee_is_finite(density)))
      CPASSERT(ALL(angular_momenta >= 0))
   END SUBROUTINE check_radial

END MODULE upf_inverse_coulomb
