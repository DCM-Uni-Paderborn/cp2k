!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Own the selected density subgradients of one KS environment.
! **************************************************************************************************
MODULE skala_density_constraints
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: dp
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'skala_density_constraints'

   PUBLIC :: skala_density_constraints_type, skala_density_constraints_reset, &
             skala_density_constraints_set, skala_density_constraints_count, &
             skala_density_constraints_matches, skala_density_constraints_get

   TYPE skala_density_constraints_type
      PRIVATE
      INTEGER :: nspins = 0, npts(3) = 0, bounds(2, 3) = 0
      INTEGER, ALLOCATABLE :: points(:, :), spins(:)
      REAL(dp) :: hmat(3, 3) = 0.0_dp, probe = 0.0_dp, tolerance = 0.0_dp
      REAL(dp), ALLOCATABLE :: positions(:, :), negative_weights(:)
   END TYPE skala_density_constraints_type

CONTAINS

! **************************************************************************************************
!> \brief Clear the selected constraints without changing their owner.
!> \param state ...
! **************************************************************************************************
   SUBROUTINE skala_density_constraints_reset(state)
      TYPE(skala_density_constraints_type), &
         INTENT(INOUT)                                   :: state

      IF (ALLOCATED(state%points)) DEALLOCATE (state%points)
      IF (ALLOCATED(state%spins)) DEALLOCATE (state%spins)
      IF (ALLOCATED(state%positions)) DEALLOCATE (state%positions)
      IF (ALLOCATED(state%negative_weights)) DEALLOCATE (state%negative_weights)
      state%nspins = 0
      state%npts = 0
      state%bounds = 0
      state%hmat = 0.0_dp
      state%probe = 0.0_dp
      state%tolerance = 0.0_dp
   END SUBROUTINE skala_density_constraints_reset

! **************************************************************************************************
!> \brief Store physical grid constraints and their selected convex weights.
!>        Spin zero denotes the paired alpha/beta density of an RKS calculation.
!>        The orbital solver supplies the weights. This routine does not select them.
!> \param state ...
!> \param points Global FFT grid indices, one column per physical constraint.
!> \param spins Zero for RKS, one or two for the physical UKS spin.
!> \param negative_weights Weight of the negative-density limiting derivative.
!> \param probe Positive branch displacement used by the original model.
!> \param tolerance Maximum absolute spin density for applying the selected derivative.
!> \param nspins ...
!> \param npts ...
!> \param bounds ...
!> \param hmat ...
!> \param positions ...
! **************************************************************************************************
   SUBROUTINE skala_density_constraints_set(state, points, spins, negative_weights, &
                                            probe, tolerance, nspins, npts, bounds, hmat, positions)
      TYPE(skala_density_constraints_type), &
         INTENT(INOUT)                                   :: state
      INTEGER, INTENT(IN)                                :: points(:, :), spins(:)
      REAL(dp), INTENT(IN)                               :: negative_weights(:), probe, tolerance
      INTEGER, INTENT(IN)                                :: nspins, npts(3), bounds(2, 3)
      REAL(dp), INTENT(IN)                               :: hmat(3, 3), positions(:, :)

      INTEGER                                            :: i, j

      CPASSERT(nspins == 1 .OR. nspins == 2)
      CPASSERT(SIZE(points, 1) == 3 .AND. SIZE(points, 2) == SIZE(spins))
      CPASSERT(SIZE(spins) == SIZE(negative_weights))
      CPASSERT(SIZE(positions, 1) == 3)
      CPASSERT(ALL(npts > 0) .AND. ALL(bounds(2, :) - bounds(1, :) + 1 == npts))
      CPASSERT(ieee_is_finite(probe) .AND. probe > 0.0_dp)
      CPASSERT(ieee_is_finite(tolerance) .AND. tolerance >= 0.0_dp)
      CPASSERT(ALL(ieee_is_finite(hmat)) .AND. ALL(ieee_is_finite(positions)))
      CPASSERT(ALL(ieee_is_finite(negative_weights)))
      CPASSERT(ALL(negative_weights >= 0.0_dp))
      CPASSERT(ALL(negative_weights <= 1.0_dp))
      IF (nspins == 1) THEN
         CPASSERT(ALL(spins == 0))
      ELSE
         CPASSERT(ALL(spins >= 1) .AND. ALL(spins <= 2))
      END IF
      DO i = 1, SIZE(spins)
         CPASSERT(ALL(points(:, i) >= bounds(1, :)))
         CPASSERT(ALL(points(:, i) <= bounds(2, :)))
         DO j = 1, i - 1
            IF (ALL(points(:, i) == points(:, j)) .AND. spins(i) == spins(j)) THEN
               CPABORT("Duplicate SKALA density constraint.")
            END IF
         END DO
      END DO

      CALL skala_density_constraints_reset(state)
      state%nspins = nspins
      state%npts = npts
      state%bounds = bounds
      state%hmat = hmat
      state%probe = probe
      state%tolerance = tolerance
      ALLOCATE (state%points(3, SIZE(spins)), state%spins(SIZE(spins)), &
                state%negative_weights(SIZE(spins)), state%positions(3, SIZE(positions, 2)))
      state%points(:, :) = points
      state%spins(:) = spins
      state%negative_weights(:) = negative_weights
      state%positions(:, :) = positions
   END SUBROUTINE skala_density_constraints_set

! **************************************************************************************************
!> \brief Number of selected physical constraints.
!> \param state ...
!> \return ...
! **************************************************************************************************
   FUNCTION skala_density_constraints_count(state) RESULT(n)
      TYPE(skala_density_constraints_type), INTENT(IN)   :: state
      INTEGER                                            :: n

      n = 0
      IF (ALLOCATED(state%spins)) n = SIZE(state%spins)
   END FUNCTION skala_density_constraints_count

! **************************************************************************************************
!> \brief Check that a selection belongs to this geometry and physical grid.
!> \param state ...
!> \param nspins ...
!> \param npts ...
!> \param bounds ...
!> \param hmat ...
!> \param positions ...
!> \return ...
! **************************************************************************************************
   FUNCTION skala_density_constraints_matches(state, nspins, npts, bounds, hmat, positions) RESULT(matches)
      TYPE(skala_density_constraints_type), INTENT(IN)   :: state
      INTEGER, INTENT(IN)                                :: nspins, npts(3), bounds(2, 3)
      REAL(dp), INTENT(IN)                               :: hmat(3, 3), positions(:, :)
      LOGICAL                                            :: matches

      matches = .FALSE.
      IF (skala_density_constraints_count(state) == 0) RETURN
      IF (state%nspins /= nspins .OR. ANY(state%npts /= npts) .OR. ANY(state%bounds /= bounds)) RETURN
      IF (ANY(state%hmat /= hmat)) RETURN
      IF (ANY(SHAPE(state%positions) /= SHAPE(positions))) RETURN
      matches = ALL(state%positions == positions)
   END FUNCTION skala_density_constraints_matches

! **************************************************************************************************
!> \brief Obtain owned copies of the selected physical constraints.
!> \param state ...
!> \param points ...
!> \param spins ...
!> \param negative_weights ...
!> \param probe ...
!> \param tolerance ...
! **************************************************************************************************
   SUBROUTINE skala_density_constraints_get(state, points, spins, negative_weights, probe, tolerance)
      TYPE(skala_density_constraints_type), INTENT(IN)   :: state
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: points(:, :), spins(:)
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: negative_weights(:)
      REAL(dp), INTENT(OUT)                              :: probe, tolerance

      INTEGER                                            :: n

      n = skala_density_constraints_count(state)
      CPASSERT(n > 0)
      ALLOCATE (points(3, n), spins(n), negative_weights(n))
      points(:, :) = state%points
      spins(:) = state%spins
      negative_weights(:) = state%negative_weights
      probe = state%probe
      tolerance = state%tolerance
   END SUBROUTINE skala_density_constraints_get

END MODULE skala_density_constraints
