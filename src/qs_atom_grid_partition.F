!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

!> \brief Shared finite-image atom partition and Cartesian geometry derivatives.
MODULE qs_atom_grid_partition
   USE cell_types,                      ONLY: cell_type
   USE kinds,                           ONLY: dp
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   REAL(dp), PARAMETER :: layout_tol = 1.0e-12_dp
   INTEGER, PARAMETER :: periodic_partition_image_shell = 1
   PUBLIC :: build_periodic_atom_image_layout, periodic_atom_image_partition
   PUBLIC :: periodic_atom_image_partition_from_layout, becke_shape, becke_shape_derivative
CONTAINS

! **************************************************************************************************
!> \brief Build screened logarithmic Becke weights for a set of periodic atom images.
!> \param image_coords ...
!> \param distances ...
!> \param log_weights ...
! **************************************************************************************************
   SUBROUTINE periodic_image_log_weights(image_coords, distances, log_weights)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: image_coords
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: distances
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: log_weights

      INTEGER                                            :: iimage, image_count, jimage, &
                                                            nearest_image
      REAL(KIND=dp)                                      :: factor, mu, rab, reference_log_weight, &
                                                            screen_log_weight
      REAL(KIND=dp), DIMENSION(3)                        :: pair

      image_count = SIZE(distances)
      CPASSERT(SIZE(image_coords, 1) == 3)
      CPASSERT(SIZE(image_coords, 2) == image_count)
      CPASSERT(SIZE(log_weights) == image_count)

      nearest_image = MINLOC(distances, DIM=1)
      log_weights = -HUGE(1.0_dp)
      reference_log_weight = 0.0_dp
      DO jimage = 1, image_count
         IF (jimage == nearest_image) CYCLE
         pair = image_coords(:, nearest_image) - image_coords(:, jimage)
         rab = SQRT(SUM(pair**2))
         IF (rab <= layout_tol) CYCLE
         mu = MAX(-1.0_dp, MIN(1.0_dp, &
                               (distances(nearest_image) - distances(jimage))/rab))
         factor = 0.5_dp*(1.0_dp - becke_shape(mu))
         CPASSERT(factor > TINY(1.0_dp))
         reference_log_weight = reference_log_weight + LOG(factor)
      END DO
      log_weights(nearest_image) = reference_log_weight

      ! Every omitted factor is at most one. A partial log product is therefore an upper bound
      ! on the final unnormalized weight. The normalization contains at least the nearest-image
      ! reference weight, so screening below this limit changes the complete partition by no more
      ! than one floating-point epsilon in total.
      screen_log_weight = reference_log_weight + &
                          LOG(EPSILON(1.0_dp)/REAL(image_count, KIND=dp))
      DO iimage = 1, image_count
         IF (iimage == nearest_image) CYCLE
         pair = image_coords(:, iimage) - image_coords(:, nearest_image)
         rab = SQRT(SUM(pair**2))
         IF (rab <= layout_tol) CYCLE
         mu = MAX(-1.0_dp, MIN(1.0_dp, &
                               (distances(iimage) - distances(nearest_image))/rab))
         factor = 0.5_dp*(1.0_dp - becke_shape(mu))
         IF (factor <= TINY(1.0_dp)) CYCLE
         log_weights(iimage) = LOG(factor)
         IF (log_weights(iimage) <= screen_log_weight) THEN
            log_weights(iimage) = -HUGE(1.0_dp)
            CYCLE
         END IF

         DO jimage = 1, image_count
            IF (jimage == iimage .OR. jimage == nearest_image) CYCLE
            pair = image_coords(:, iimage) - image_coords(:, jimage)
            rab = SQRT(SUM(pair**2))
            IF (rab <= layout_tol) CYCLE
            mu = MAX(-1.0_dp, MIN(1.0_dp, &
                                  (distances(iimage) - distances(jimage))/rab))
            factor = 0.5_dp*(1.0_dp - becke_shape(mu))
            IF (factor <= TINY(1.0_dp)) THEN
               log_weights(iimage) = -HUGE(1.0_dp)
               EXIT
            END IF
            log_weights(iimage) = log_weights(iimage) + LOG(factor)
            IF (log_weights(iimage) <= screen_log_weight) THEN
               log_weights(iimage) = -HUGE(1.0_dp)
               EXIT
            END IF
         END DO
      END DO

   END SUBROUTINE periodic_image_log_weights

! **************************************************************************************************
!> \brief Build screened normalized Becke weights, using direct products whenever they are safe.
!> \param image_coords ...
!> \param distances ...
!> \param normalized_weights ...
! **************************************************************************************************
   SUBROUTINE periodic_image_normalized_weights(image_coords, distances, normalized_weights)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: image_coords
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: distances
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: normalized_weights

      INTEGER                                            :: iimage, image_count, jimage, &
                                                            nearest_image
      LOGICAL                                            :: use_log_fallback
      REAL(KIND=dp)                                      :: factor, max_log_weight, mu, &
                                                            normalization, product, rab, &
                                                            reference_weight, screen_factor, &
                                                            screen_weight
      REAL(KIND=dp), DIMENSION(3)                        :: pair
      REAL(KIND=dp), DIMENSION(SIZE(distances))          :: log_weights, raw_weights

      image_count = SIZE(distances)
      CPASSERT(SIZE(image_coords, 1) == 3)
      CPASSERT(SIZE(image_coords, 2) == image_count)
      CPASSERT(SIZE(normalized_weights) == image_count)
      IF (image_count == 1) THEN
         normalized_weights(1) = 1.0_dp
         RETURN
      END IF

      nearest_image = MINLOC(distances, DIM=1)
      reference_weight = 1.0_dp
      use_log_fallback = .FALSE.
      DO jimage = 1, image_count
         IF (jimage == nearest_image) CYCLE
         pair = image_coords(:, nearest_image) - image_coords(:, jimage)
         rab = SQRT(SUM(pair**2))
         IF (rab <= layout_tol) CYCLE
         mu = MAX(-1.0_dp, MIN(1.0_dp, &
                               (distances(nearest_image) - distances(jimage))/rab))
         factor = 0.5_dp*(1.0_dp - becke_shape(mu))
         CPASSERT(factor > TINY(1.0_dp))
         IF (reference_weight <= TINY(1.0_dp)/factor) THEN
            use_log_fallback = .TRUE.
            EXIT
         END IF
         reference_weight = reference_weight*factor
      END DO

      screen_factor = EPSILON(1.0_dp)/REAL(image_count, KIND=dp)
      IF (.NOT. use_log_fallback) THEN
         use_log_fallback = reference_weight <= TINY(1.0_dp)/screen_factor
      END IF
      IF (use_log_fallback) THEN
         CALL periodic_image_log_weights(image_coords, distances, log_weights)
         max_log_weight = MAXVAL(log_weights)
         CPASSERT(max_log_weight > -0.5_dp*HUGE(1.0_dp))
         normalized_weights = 0.0_dp
         WHERE (log_weights > -0.5_dp*HUGE(1.0_dp))
            normalized_weights = EXP(log_weights - max_log_weight)
         END WHERE
         normalization = SUM(normalized_weights)
         CPASSERT(normalization > 0.0_dp)
         normalized_weights = normalized_weights/normalization
         RETURN
      END IF

      screen_weight = reference_weight*screen_factor
      raw_weights = 0.0_dp
      raw_weights(nearest_image) = reference_weight
      DO iimage = 1, image_count
         IF (iimage == nearest_image) CYCLE
         pair = image_coords(:, iimage) - image_coords(:, nearest_image)
         rab = SQRT(SUM(pair**2))
         IF (rab <= layout_tol) CYCLE
         mu = MAX(-1.0_dp, MIN(1.0_dp, &
                               (distances(iimage) - distances(nearest_image))/rab))
         factor = 0.5_dp*(1.0_dp - becke_shape(mu))
         IF (factor <= TINY(1.0_dp) .OR. factor <= screen_weight) CYCLE
         product = factor

         DO jimage = 1, image_count
            IF (jimage == iimage .OR. jimage == nearest_image) CYCLE
            pair = image_coords(:, iimage) - image_coords(:, jimage)
            rab = SQRT(SUM(pair**2))
            IF (rab <= layout_tol) CYCLE
            mu = MAX(-1.0_dp, MIN(1.0_dp, &
                                  (distances(iimage) - distances(jimage))/rab))
            factor = 0.5_dp*(1.0_dp - becke_shape(mu))
            IF (factor <= TINY(1.0_dp) .OR. product <= screen_weight/factor) THEN
               product = 0.0_dp
               EXIT
            END IF
            product = product*factor
         END DO
         raw_weights(iimage) = product
      END DO

      normalization = SUM(raw_weights)
      CPASSERT(normalization > 0.0_dp)
      normalized_weights = raw_weights/normalization

   END SUBROUTINE periodic_image_normalized_weights

! **************************************************************************************************
!> \brief Build the image coordinates shared by all points of one target-atom block.
!> \param atom_coords ...
!> \param cell ...
!> \param target_atom ...
!> \param image_periodicity optional image directions used by the auxiliary-grid partition
!> \param image_coords ...
!> \param target_image ...
!> \param image_atom Optional base-atom index of each image.
!> \param image_shift Optional lattice shift of each image.
!> \param image_shell Optional nonnegative shell extent along each lattice direction.
! **************************************************************************************************
   SUBROUTINE build_periodic_atom_image_layout(atom_coords, cell, target_atom, image_periodicity, &
                                               image_coords, target_image, image_atom, image_shift, image_shell)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: atom_coords
      TYPE(cell_type), POINTER                           :: cell
      INTEGER, INTENT(IN)                                :: target_atom
      INTEGER, DIMENSION(3), INTENT(IN), OPTIONAL        :: image_periodicity
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: image_coords
      INTEGER, INTENT(OUT)                               :: target_image
      INTEGER, ALLOCATABLE, INTENT(OUT), OPTIONAL        :: image_atom(:), image_shift(:, :)
      INTEGER, INTENT(IN), OPTIONAL                      :: image_shell(3)

      INTEGER                                            :: base_shift(3), iatom, idir, image_count, &
                                                            image_index, ix, ix_max, ix_min, iy, &
                                                            iy_max, iy_min, iz, iz_max, iz_min, &
                                                            natom
      INTEGER, DIMENSION(3)                              :: partition_periodicity, shell
      REAL(KIND=dp), DIMENSION(3)                        :: fractional

      natom = SIZE(atom_coords, 2)
      CPASSERT(ASSOCIATED(cell))
      CPASSERT(target_atom >= 1 .AND. target_atom <= natom)
      CPASSERT(SIZE(atom_coords, 1) == 3)
      partition_periodicity = cell%perd
      IF (PRESENT(image_periodicity)) partition_periodicity = image_periodicity
      CPASSERT(ALL(partition_periodicity == 0 .OR. partition_periodicity == 1))

      shell = periodic_partition_image_shell
      IF (PRESENT(image_shell)) shell = image_shell
      CPASSERT(ALL(shell >= 0))
      image_count = natom
      DO idir = 1, 3
         IF (partition_periodicity(idir) == 1) THEN
            image_count = image_count*(2*shell(idir) + 1)
         END IF
      END DO
      ALLOCATE (image_coords(3, image_count))
      IF (PRESENT(image_atom)) ALLOCATE (image_atom(image_count))
      IF (PRESENT(image_shift)) ALLOCATE (image_shift(3, image_count))

      image_index = 0
      target_image = 0
      DO iatom = 1, natom
         fractional = MATMUL(cell%h_inv, &
                             atom_coords(:, target_atom) - atom_coords(:, iatom))
         DO idir = 1, 3
            base_shift(idir) = partition_periodicity(idir)*NINT(fractional(idir))
         END DO
         ix_min = base_shift(1) - shell(1)*partition_periodicity(1)
         ix_max = base_shift(1) + shell(1)*partition_periodicity(1)
         iy_min = base_shift(2) - shell(2)*partition_periodicity(2)
         iy_max = base_shift(2) + shell(2)*partition_periodicity(2)
         iz_min = base_shift(3) - shell(3)*partition_periodicity(3)
         iz_max = base_shift(3) + shell(3)*partition_periodicity(3)
         DO iz = iz_min, iz_max
            DO iy = iy_min, iy_max
               DO ix = ix_min, ix_max
                  image_index = image_index + 1
                  IF (PRESENT(image_atom)) image_atom(image_index) = iatom
                  IF (PRESENT(image_shift)) image_shift(:, image_index) = [ix, iy, iz]
                  image_coords(:, image_index) = atom_coords(:, iatom) + &
                                                 MATMUL(cell%hmat, REAL([ix, iy, iz], KIND=dp))
                  IF (iatom == target_atom .AND. ALL([ix, iy, iz] == 0)) THEN
                     target_image = image_index
                  END IF
               END DO
            END DO
         END DO
      END DO
      CPASSERT(image_index == image_count)
      CPASSERT(target_image > 0)

   END SUBROUTINE build_periodic_atom_image_layout

! **************************************************************************************************
!> \brief Return an image-complete periodic atom weight for a prebuilt image layout.
!> \param grid_point ...
!> \param image_coords ...
!> \param target_image ...
!> \param weight ...
!> \param dweight_dpoint Optional Cartesian query derivative, supplied with dweight_dimage.
!> \param dweight_dimage Optional Cartesian derivative for each independent image coordinate.
! **************************************************************************************************
   SUBROUTINE periodic_atom_image_partition_from_layout(grid_point, image_coords, target_image, &
                                                        weight, dweight_dpoint, dweight_dimage)
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: grid_point
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: image_coords
      INTEGER, INTENT(IN)                                :: target_image
      REAL(KIND=dp), INTENT(OUT)                         :: weight
      REAL(dp), INTENT(OUT), OPTIONAL                    :: dweight_dpoint(3), dweight_dimage(:, :)

      INTEGER                                            :: active_count, active_image, &
                                                            active_index, candidate_image, iimage, &
                                                            image_count, jimage
      INTEGER, ALLOCATABLE                               :: active_images(:)
      REAL(dp) :: dist_diff, dmu_image_i(3), dmu_image_j(3), ds_dmu, ds_image_i(3), ds_image_j(3), &
         mu, mu_raw, one_minus_switch, pair(3), rab, switch, unit_pair(3)
      REAL(dp), ALLOCATABLE                              :: mean_log_derivative(:, :), &
                                                            target_log_derivative(:, :), &
                                                            unit_rvecs(:, :)
      REAL(KIND=dp), DIMENSION(SIZE(image_coords, 2))    :: distances, normalized_weights

      image_count = SIZE(image_coords, 2)
      CPASSERT(SIZE(image_coords, 1) == 3)
      CPASSERT(target_image >= 1 .AND. target_image <= image_count)
      CPASSERT(PRESENT(dweight_dpoint) .EQV. PRESENT(dweight_dimage))
      IF (PRESENT(dweight_dimage)) THEN
         CPASSERT(ALL(SHAPE(dweight_dimage) == [3, image_count]))
         dweight_dimage = 0.0_dp
         dweight_dpoint = 0.0_dp
      END IF
      IF (image_count == 1) THEN
         weight = 1.0_dp
         RETURN
      END IF

      weight = 0.0_dp
      DO iimage = 1, image_count
         distances(iimage) = SQRT(SUM((grid_point - image_coords(:, iimage))**2))
      END DO
      CALL periodic_image_normalized_weights(image_coords, distances, normalized_weights)
      weight = normalized_weights(target_image)
      IF (.NOT. PRESENT(dweight_dimage) .OR. weight <= 0.0_dp) RETURN
      ALLOCATE (active_images(image_count), unit_rvecs(3, image_count), &
                mean_log_derivative(3, image_count), target_log_derivative(3, image_count))
      DO iimage = 1, image_count
         IF (distances(iimage) > layout_tol) THEN
            unit_rvecs(:, iimage) = (grid_point - image_coords(:, iimage))/distances(iimage)
         ELSE
            unit_rvecs(:, iimage) = 0.0_dp
         END IF
      END DO
      active_count = 0
      DO iimage = 1, image_count
         IF (normalized_weights(iimage) > 0.0_dp) THEN
            active_count = active_count + 1
            active_images(active_count) = iimage
         END IF
      END DO
      CPASSERT(active_count > 0)
      ! Revisit the image pairs after normalization to contract the logarithmic
      ! weight derivatives directly, avoiding both an image_count-squared tensor
      ! and pairs for which both screened normalized weights vanish.
      mean_log_derivative = 0.0_dp
      target_log_derivative = 0.0_dp
      DO active_index = 1, active_count
         active_image = active_images(active_index)
         DO candidate_image = 1, image_count
            IF (candidate_image == active_image) CYCLE
            ! A pair of active images is visited by the lower-index image only.
            IF (normalized_weights(candidate_image) > 0.0_dp .AND. &
                candidate_image < active_image) CYCLE
            iimage = MIN(active_image, candidate_image)
            jimage = MAX(active_image, candidate_image)
            pair = image_coords(:, iimage) - image_coords(:, jimage)
            rab = SQRT(SUM(pair**2))
            IF (rab <= layout_tol) CYCLE
            unit_pair = pair/rab
            dist_diff = distances(iimage) - distances(jimage)
            mu_raw = dist_diff/rab
            mu = MAX(-1.0_dp, MIN(1.0_dp, mu_raw))
            switch = 0.5_dp*(1.0_dp - becke_shape(mu))
            one_minus_switch = 1.0_dp - switch
            IF (ABS(mu_raw) < 1.0_dp) THEN
               ds_dmu = -0.5_dp*becke_shape_derivative(mu)
            ELSE
               ds_dmu = 0.0_dp
            END IF
            IF (ABS(ds_dmu) <= 0.0_dp .OR. switch <= TINY(1.0_dp) .OR. &
                one_minus_switch <= TINY(1.0_dp)) CYCLE
            dmu_image_i = (-unit_rvecs(:, iimage)*rab - dist_diff*unit_pair)/rab**2
            dmu_image_j = (unit_rvecs(:, jimage)*rab + dist_diff*unit_pair)/rab**2
            ds_image_i = ds_dmu*dmu_image_i
            ds_image_j = ds_dmu*dmu_image_j

            mean_log_derivative(:, iimage) = mean_log_derivative(:, iimage) + &
                                             normalized_weights(iimage)*ds_image_i/switch - &
                                             normalized_weights(jimage)*ds_image_i/one_minus_switch
            mean_log_derivative(:, jimage) = mean_log_derivative(:, jimage) + &
                                             normalized_weights(iimage)*ds_image_j/switch - &
                                             normalized_weights(jimage)*ds_image_j/one_minus_switch
            IF (target_image == iimage) THEN
               target_log_derivative(:, iimage) = &
                  target_log_derivative(:, iimage) + ds_image_i/switch
               target_log_derivative(:, jimage) = &
                  target_log_derivative(:, jimage) + ds_image_j/switch
            ELSE IF (target_image == jimage) THEN
               target_log_derivative(:, iimage) = &
                  target_log_derivative(:, iimage) - ds_image_i/one_minus_switch
               target_log_derivative(:, jimage) = &
                  target_log_derivative(:, jimage) - ds_image_j/one_minus_switch
            END IF
         END DO
      END DO
      DO jimage = 1, image_count
         dweight_dimage(:, jimage) = weight*(target_log_derivative(:, jimage) - mean_log_derivative(:, jimage))
      END DO
      ! A common translation of query and images leaves the partition unchanged.
      dweight_dpoint = -SUM(dweight_dimage, DIM=2)

   END SUBROUTINE periodic_atom_image_partition_from_layout

! **************************************************************************************************
!> \brief Return the reference atom's weight and base-atom and image-strain derivatives.
!> \param grid_point Cartesian query, held fixed in the returned derivatives.
!> \param atom_coords Base-atom Cartesian coordinates.
!> \param cell Periodic cell.
!> \param target_atom Reference atom whose zero-shift image is selected.
!> \param weight Normalized image weight.
!> \param dweight_datom Optional derivatives for base atoms, with the query held fixed.
!> \param dweight_dstrain Optional derivative due only to lattice-image translations.
!> \param image_periodicity Optional periodic directions of the image layout.
! **************************************************************************************************
   SUBROUTINE periodic_atom_image_partition(grid_point, atom_coords, cell, target_atom, weight, &
                                            dweight_datom, dweight_dstrain, image_periodicity)
      REAL(dp), INTENT(IN)                               :: grid_point(3), atom_coords(:, :)
      TYPE(cell_type), POINTER                           :: cell
      INTEGER, INTENT(IN)                                :: target_atom
      REAL(dp), INTENT(OUT)                              :: weight
      REAL(dp), INTENT(OUT), OPTIONAL                    :: dweight_datom(:, :), &
                                                            dweight_dstrain(:, :)
      INTEGER, INTENT(IN), OPTIONAL                      :: image_periodicity(3)

      INTEGER                                            :: idir, jdir, jimage, target_image
      INTEGER, ALLOCATABLE                               :: image_atom(:), image_shift(:, :)
      REAL(dp)                                           :: dpoint(3), image_translation(3)
      REAL(dp), ALLOCATABLE                              :: dimage(:, :), image_coords(:, :)

      CPASSERT(PRESENT(dweight_datom) .EQV. PRESENT(dweight_dstrain))
      CALL build_periodic_atom_image_layout(atom_coords, cell, target_atom, image_periodicity, &
                                            image_coords, target_image, image_atom, image_shift)
      IF (.NOT. PRESENT(dweight_datom)) THEN
         CALL periodic_atom_image_partition_from_layout(grid_point, image_coords, target_image, weight)
         RETURN
      END IF
      CPASSERT(ALL(SHAPE(dweight_datom) == [3, SIZE(atom_coords, 2)]))
      CPASSERT(ALL(SHAPE(dweight_dstrain) == [3, 3]))
      ALLOCATE (dimage(3, SIZE(image_atom)))
      CALL periodic_atom_image_partition_from_layout(grid_point, image_coords, target_image, weight, dpoint, dimage)
      dweight_datom = 0.0_dp
      dweight_dstrain = 0.0_dp
      DO jimage = 1, SIZE(image_atom)
         dweight_datom(:, image_atom(jimage)) = dweight_datom(:, image_atom(jimage)) + dimage(:, jimage)
         image_translation = MATMUL(cell%hmat, REAL(image_shift(:, jimage), dp))
         DO idir = 1, 3
            DO jdir = 1, 3
               dweight_dstrain(idir, jdir) = dweight_dstrain(idir, jdir) + dimage(idir, jimage)*image_translation(jdir)
            END DO
         END DO
      END DO
   END SUBROUTINE periodic_atom_image_partition

! **************************************************************************************************
!> \brief Becke fuzzy-cell shape function.
!> \param mu Difference of distances divided by pair separation.
!> \return Threefold cubic switch.
! **************************************************************************************************
   PURE FUNCTION becke_shape(mu) RESULT(val)
      REAL(KIND=dp), INTENT(IN)                          :: mu
      REAL(KIND=dp)                                      :: val

      INTEGER                                            :: iter

      val = mu
      DO iter = 1, 3
         val = 0.5_dp*val*(3.0_dp - val*val)
      END DO

   END FUNCTION becke_shape

! **************************************************************************************************
!> \brief Derivative of the threefold Becke switch.
!> \param mu Difference of distances divided by pair separation.
!> \return Derivative with respect to mu.
! **************************************************************************************************
   PURE FUNCTION becke_shape_derivative(mu) RESULT(val)
      REAL(KIND=dp), INTENT(IN)                          :: mu
      REAL(KIND=dp)                                      :: val

      INTEGER                                            :: iter
      REAL(KIND=dp)                                      :: x

      x = mu
      val = 1.0_dp
      DO iter = 1, 3
         val = val*1.5_dp*(1.0_dp - x*x)
         x = 0.5_dp*x*(3.0_dp - x*x)
      END DO

   END FUNCTION becke_shape_derivative

END MODULE qs_atom_grid_partition
