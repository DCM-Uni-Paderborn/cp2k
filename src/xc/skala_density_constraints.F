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
             skala_density_constraints_matches, skala_density_constraints_get, &
             skala_density_constraints_balance, skala_density_crossings, skala_density_root

   ABSTRACT INTERFACE
      SUBROUTINE evaluate_physical_density(step, value, context)
         IMPORT :: dp
         REAL(dp), INTENT(IN) :: step
         REAL(dp), INTENT(OUT) :: value
         CLASS(*), INTENT(INOUT) :: context
      END SUBROUTINE evaluate_physical_density
   END INTERFACE

   TYPE skala_density_constraints_type
      PRIVATE
      INTEGER :: nspins = 0, npts(3) = 0, bounds(2, 3) = 0
      INTEGER, ALLOCATABLE :: points(:, :), spins(:)
      REAL(dp) :: hmat(3, 3) = 0.0_dp, probe = 0.0_dp, tolerance = 0.0_dp
      REAL(dp), ALLOCATABLE :: positions(:, :), negative_weights(:)
   END TYPE skala_density_constraints_type

CONTAINS

! **************************************************************************************************
!> \brief Detect physical density zeros bracketed by two orbital evaluations.
!>        The caller supplies matching physical point/spin layouts, before nonlinear features.
!>        A field that is zero at both ends supplies no sign-changing constraint.
!>        Endpoint zeros are candidates, not proof of an admissible stationary minimum.
!> \param first Physical spin densities at the first orbital state.
!> \param last Physical spin densities at the last orbital state.
!> \param crossings Candidate mask with the same physical point/spin layout.
! **************************************************************************************************
   SUBROUTINE skala_density_crossings(first, last, crossings)
      REAL(dp), INTENT(IN)                               :: first(:, :), last(:, :)
      LOGICAL, INTENT(OUT)                               :: crossings(:, :)

      CPASSERT(ALL(SHAPE(first) == SHAPE(last)))
      CPASSERT(ALL(SHAPE(first) == SHAPE(crossings)))
      CPASSERT(ALL(ieee_is_finite(first)))
      CPASSERT(ALL(ieee_is_finite(last)))
      crossings(:, :) = ((first < 0.0_dp .AND. last >= 0.0_dp) .OR. &
                         (first > 0.0_dp .AND. last <= 0.0_dp) .OR. &
                         (first == 0.0_dp .AND. last /= 0.0_dp))
   END SUBROUTINE skala_density_crossings

! **************************************************************************************************
!> \brief Locate a bracketed zero by reevaluating the physical density along an orbital curve.
!>        The callback retains the actual normalization, density map and frozen core.
!>        Only the absolute physical density tolerance establishes convergence.
!>        Unbracketed and unresolved roots return the best observed endpoint with failure.
!> \param evaluate Callback for one fixed physical point/spin on the caller's orbital curve.
!> \param context Owned orbital-curve data and native evaluation workspace.
!> \param first_step First endpoint of the supplied orbital-curve interval.
!> \param last_step Last endpoint, strictly larger than the first.
!> \param tolerance Maximum absolute physical spin density at the returned root.
!> \param step Best observed orbital-curve parameter.
!> \param value Its actual physical spin density.
!> \param converged Whether the physical density meets the requested tolerance.
!> \param max_iterations Maximum safeguarded interpolation steps, defaults to 100.
! **************************************************************************************************
   SUBROUTINE skala_density_root(evaluate, context, first_step, last_step, tolerance, step, value, &
                                 converged, max_iterations)
      PROCEDURE(evaluate_physical_density)                :: evaluate
      CLASS(*), INTENT(INOUT)                             :: context
      REAL(dp), INTENT(IN)                               :: first_step, last_step, tolerance
      REAL(dp), INTENT(OUT)                              :: step, value
      LOGICAL, INTENT(OUT)                               :: converged
      INTEGER, INTENT(IN), OPTIONAL                      :: max_iterations

      INTEGER                                            :: iteration, limit
      LOGICAL                                            :: bracketed
      REAL(dp)                                           :: first, first_value, fraction, last, last_value, &
                                                            scale, trial, trial_value

      CPASSERT(ieee_is_finite(first_step) .AND. ieee_is_finite(last_step))
      CPASSERT(last_step > first_step)
      CPASSERT(ieee_is_finite(tolerance) .AND. tolerance > 0.0_dp)
      limit = 100
      IF (PRESENT(max_iterations)) limit = max_iterations
      CPASSERT(limit > 0)
      first = first_step
      last = last_step
      CALL evaluate(first, first_value, context)
      CALL evaluate(last, last_value, context)
      CPASSERT(ieee_is_finite(first_value) .AND. ieee_is_finite(last_value))
      step = first
      value = first_value
      IF (ABS(last_value) < ABS(value)) THEN
         step = last
         value = last_value
      END IF
      converged = ABS(value) <= tolerance
      IF (converged) RETURN
      bracketed = (first_value < 0.0_dp .AND. last_value > 0.0_dp) .OR. &
                  (first_value > 0.0_dp .AND. last_value < 0.0_dp)
      IF (.NOT. bracketed) RETURN

      DO iteration = 1, limit
         scale = MAX(ABS(first_value), ABS(last_value))
         fraction = (first_value/scale)/(first_value/scale - last_value/scale)
         fraction = MAX(0.1_dp, MIN(0.9_dp, fraction))
         trial = (1.0_dp - fraction)*first + fraction*last
         IF (trial <= first .OR. trial >= last) EXIT
         CALL evaluate(trial, trial_value, context)
         CPASSERT(ieee_is_finite(trial_value))
         IF (ABS(trial_value) < ABS(value)) THEN
            step = trial
            value = trial_value
         END IF
         converged = ABS(value) <= tolerance
         IF (converged) RETURN
         IF ((trial_value < 0.0_dp) .EQV. (first_value < 0.0_dp)) THEN
            first = trial
            first_value = trial_value
         ELSE
            last = trial
            last_value = trial_value
         END IF
      END DO
   END SUBROUTINE skala_density_root

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

! **************************************************************************************************
!> \brief Select the minimum-norm convex orbital derivative at supplied density constraints.
!>        The columns are independent negative-minus-positive limiting gradient changes.
!>        All gradients must use the same orthonormal real orbital tangent coordinates.
!>        This solves a convex box problem, not the density roots or orbital stationarity.
!> \param state Physical constraints whose weights are updated only on convergence.
!> \param positive_gradient Complete orbital gradient with every constraint on its positive branch.
!> \param negative_changes One gradient change per physical constraint, in the state's order.
!> \param tolerance Absolute tolerance for projected gradients along normalized change columns.
!> \param balanced_gradient Complete orbital gradient using the state's resulting weights.
!> \param converged Whether all box optimality conditions meet tolerance.
!> \param max_iterations Maximum bounded minimization sweeps, defaults to 10000.
! **************************************************************************************************
   SUBROUTINE skala_density_constraints_balance(state, positive_gradient, negative_changes, tolerance, &
                                                balanced_gradient, converged, max_iterations)
      TYPE(skala_density_constraints_type), &
         INTENT(INOUT)                                   :: state
      REAL(dp), INTENT(IN)                               :: positive_gradient(:), &
                                                            negative_changes(:, :), tolerance
      REAL(dp), INTENT(OUT)                              :: balanced_gradient(:)
      LOGICAL, INTENT(OUT)                               :: converged
      INTEGER, INTENT(IN), OPTIONAL                      :: max_iterations

      INTEGER                                            :: c, iteration, limit, n
      REAL(dp)                                           :: derivative, error, model_scale, &
                                                            new_weight, scaled_tolerance
      REAL(dp), ALLOCATABLE                              :: changes(:, :), column_norms(:), &
                                                            gradient(:), weights(:)

      n = skala_density_constraints_count(state)
      CPASSERT(SIZE(negative_changes, 1) == SIZE(positive_gradient))
      CPASSERT(SIZE(negative_changes, 2) == n)
      CPASSERT(SIZE(balanced_gradient) == SIZE(positive_gradient))
      CPASSERT(ALL(ieee_is_finite(positive_gradient)))
      CPASSERT(ALL(ieee_is_finite(negative_changes)))
      CPASSERT(ieee_is_finite(tolerance) .AND. tolerance > 0.0_dp)
      limit = 10000
      IF (PRESENT(max_iterations)) limit = max_iterations
      CPASSERT(limit > 0)
      balanced_gradient(:) = positive_gradient
      converged = .TRUE.
      IF (n == 0 .OR. SIZE(positive_gradient) == 0) RETURN

      model_scale = MAX(MAXVAL(ABS(positive_gradient)), MAXVAL(ABS(negative_changes)))
      IF (model_scale == 0.0_dp) RETURN
      ALLOCATE (changes(SIZE(positive_gradient), n), gradient(SIZE(positive_gradient)), &
                column_norms(n), weights(n))
      changes(:, :) = negative_changes/model_scale
      weights(:) = state%negative_weights
      scaled_tolerance = 1.0_dp
      IF (model_scale > tolerance) scaled_tolerance = tolerance/model_scale
      DO c = 1, n
         column_norms(c) = DOT_PRODUCT(changes(:, c), changes(:, c))
      END DO

      converged = .FALSE.
      DO iteration = 1, limit
         gradient(:) = positive_gradient/model_scale + MATMUL(changes, weights)
         DO c = 1, n
            IF (column_norms(c) <= TINY(0.0_dp)) CYCLE
            derivative = DOT_PRODUCT(changes(:, c), gradient)
            new_weight = MAX(0.0_dp, MIN(1.0_dp, weights(c) - derivative/column_norms(c)))
            gradient(:) = gradient + (new_weight - weights(c))*changes(:, c)
            weights(c) = new_weight
         END DO
         ! Rebuild the residual before assessing the coupled optimality conditions.
         gradient(:) = positive_gradient/model_scale + MATMUL(changes, weights)
         CALL balance_free_weights(changes, column_norms, gradient, weights)
         gradient(:) = positive_gradient/model_scale + MATMUL(changes, weights)
         error = 0.0_dp
         DO c = 1, n
            IF (column_norms(c) <= TINY(0.0_dp)) CYCLE
            derivative = DOT_PRODUCT(changes(:, c), gradient)/SQRT(column_norms(c))
            IF (weights(c) == 0.0_dp) derivative = MIN(derivative, 0.0_dp)
            IF (weights(c) == 1.0_dp) derivative = MAX(derivative, 0.0_dp)
            error = MAX(error, ABS(derivative))
         END DO
         IF (error <= scaled_tolerance) THEN
            converged = .TRUE.
            EXIT
         END IF
      END DO
      IF (converged) state%negative_weights(:) = weights
      balanced_gradient(:) = positive_gradient + MATMUL(negative_changes, state%negative_weights)
   END SUBROUTINE skala_density_constraints_balance

! **************************************************************************************************
!> \brief Minimize the free-weight subproblem without forming a Gram matrix.
!>        A bounded SVD step resolves nearly dependent gradient changes and supplements
!>        the coordinate steps that release active bounds. Zero columns retain their weight.
!> \param changes ...
!> \param column_norms ...
!> \param gradient ...
!> \param weights ...
! **************************************************************************************************
   SUBROUTINE balance_free_weights(changes, column_norms, gradient, weights)
      REAL(dp), INTENT(IN)                               :: changes(:, :), column_norms(:), &
                                                            gradient(:)
      REAL(dp), INTENT(INOUT)                            :: weights(:)

      INTEGER                                            :: c, i, info, lwork, m, n, rank
      INTEGER, ALLOCATABLE                               :: free(:)
      REAL(dp)                                           :: fraction, query(1)
      REAL(dp), ALLOCATABLE                              :: matrix(:, :), residual(:), rhs(:), &
                                                            singular(:), work(:)

      n = COUNT(weights > 0.0_dp .AND. weights < 1.0_dp .AND. column_norms > TINY(0.0_dp))
      IF (n < 2) RETURN
      m = SIZE(gradient)
      ALLOCATE (free(n), matrix(m, n), rhs(MAX(m, n)), singular(MIN(m, n)), residual(m))
      i = 0
      DO c = 1, SIZE(weights)
         IF (weights(c) <= 0.0_dp .OR. weights(c) >= 1.0_dp .OR. column_norms(c) <= TINY(0.0_dp)) CYCLE
         i = i + 1
         free(i) = c
         matrix(:, i) = changes(:, c)
      END DO
      rhs = 0.0_dp
      rhs(1:m) = -gradient
      CALL dgelss(m, n, 1, matrix, m, rhs, MAX(m, n), singular, EPSILON(0.0_dp), rank, query, -1, info)
      IF (info /= 0) RETURN
      lwork = MAX(1, CEILING(query(1)))
      ALLOCATE (work(lwork))
      CALL dgelss(m, n, 1, matrix, m, rhs, MAX(m, n), singular, EPSILON(0.0_dp), rank, work, lwork, info)
      IF (info /= 0 .OR. .NOT. ALL(ieee_is_finite(rhs(1:n)))) RETURN
      fraction = 1.0_dp
      DO i = 1, n
         c = free(i)
         IF (rhs(i) > 1.0_dp - weights(c)) fraction = MIN(fraction, (1.0_dp - weights(c))/rhs(i))
         IF (rhs(i) < -weights(c)) fraction = MIN(fraction, -weights(c)/rhs(i))
      END DO
      residual(:) = gradient
      DO i = 1, n
         residual(:) = residual + fraction*rhs(i)*changes(:, free(i))
      END DO
      IF (DOT_PRODUCT(residual, residual) > DOT_PRODUCT(gradient, gradient)) RETURN
      DO i = 1, n
         c = free(i)
         weights(c) = MAX(0.0_dp, MIN(1.0_dp, weights(c) + fraction*rhs(i)))
      END DO
   END SUBROUTINE balance_free_weights

END MODULE skala_density_constraints
