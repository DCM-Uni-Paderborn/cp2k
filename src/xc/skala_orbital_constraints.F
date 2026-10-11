!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Solve orbital stationarity on independent physical descriptor constraints.
! **************************************************************************************************
MODULE skala_orbital_constraints
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: dp
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'skala_orbital_constraints'

   PUBLIC :: skala_constrained_stationarity, skala_cusp_admissibility
   PUBLIC :: skala_orbital_chart_type, skala_orbital_chart_create, skala_orbital_chart_release, &
             skala_orbital_chart_size, skala_orbital_chart_evaluate, skala_orbital_chart_derivative

   TYPE skala_orbital_chart_type
      PRIVATE
      COMPLEX(dp), ALLOCATABLE :: base(:, :), coefficients(:, :), inverse(:, :), rotation(:, :), &
                                  work(:, :), transformed(:, :)
      REAL(dp), ALLOCATABLE :: occupations(:), coordinate(:)
      INTEGER, ALLOCATABLE :: pairs(:, :), component(:), pivots(:)
      REAL(dp) :: metric_weight = 1.0_dp
      LOGICAL :: evaluated = .FALSE.
   END TYPE skala_orbital_chart_type

   ABSTRACT INTERFACE
      SUBROUTINE physical_constraints(coordinate, values, jacobian, context)
         IMPORT :: dp
         REAL(dp), INTENT(IN) :: coordinate(:)
         REAL(dp), INTENT(OUT) :: values(:), jacobian(:, :)
         CLASS(*), INTENT(INOUT) :: context
      END SUBROUTINE physical_constraints
      SUBROUTINE orbital_energy_gradient(coordinate, energy, gradient, context)
         IMPORT :: dp
         REAL(dp), INTENT(IN) :: coordinate(:)
         REAL(dp), INTENT(OUT) :: energy, gradient(:)
         CLASS(*), INTENT(INOUT) :: context
      END SUBROUTINE orbital_energy_gradient
   END INTERFACE

CONTAINS

! **************************************************************************************************
!> \brief Define a fixed-occupation orbital chart in the original charge metric.
!>        The supplied columns span the caller's retained orbital space. Cholesky normalization
!>        uses the supplied metric, which need not be the ordinary Gaussian overlap. Only pairs
!>        with different occupations define physical rotations. A complex chart retains their
!>        independent real and imaginary rotations. Create a new chart after a geometry change.
!> \param coefficients Reference columns, including the retained unoccupied space.
!> \param metric Original Hermitian charge-normalization operator.
!> \param occupations Fixed nonnegative occupations of the reference columns.
!> \param complex_orbitals Whether to include independent imaginary rotations.
!> \param chart Owned normalized reference, rotation and derivative storage.
!> \param valid Whether the reference Gram matrix has positive rank.
!> \param metric_weight Positive k-point weight for the physical occupied-orbital tangent norm.
! **************************************************************************************************
   SUBROUTINE skala_orbital_chart_create(coefficients, metric, occupations, complex_orbitals, chart, valid, metric_weight)
      COMPLEX(dp), INTENT(IN)                            :: coefficients(:, :), metric(:, :)
      REAL(dp), INTENT(IN)                               :: occupations(:)
      LOGICAL, INTENT(IN)                                :: complex_orbitals
      TYPE(skala_orbital_chart_type), INTENT(INOUT)      :: chart
      LOGICAL, INTENT(OUT)                               :: valid
      REAL(dp), INTENT(IN), OPTIONAL                     :: metric_weight

      COMPLEX(dp), ALLOCATABLE                           :: gram(:, :)
      INTEGER                                            :: a, b, count, i, info, n, rows

      CALL skala_orbital_chart_release(chart)
      rows = SIZE(coefficients, 1)
      n = SIZE(coefficients, 2)
      CPASSERT(rows > 0 .AND. n > 0 .AND. n <= rows)
      CPASSERT(SIZE(metric, 1) == rows .AND. SIZE(metric, 2) == rows)
      CPASSERT(SIZE(occupations) == n)
      CPASSERT(ALL(ieee_is_finite(REAL(coefficients, dp))))
      CPASSERT(ALL(ieee_is_finite(AIMAG(coefficients))))
      CPASSERT(ALL(ieee_is_finite(REAL(metric, dp))))
      CPASSERT(ALL(ieee_is_finite(AIMAG(metric))))
      CPASSERT(ALL(ieee_is_finite(occupations)) .AND. ALL(occupations >= 0.0_dp))
      IF (PRESENT(metric_weight)) chart%metric_weight = metric_weight
      CPASSERT(ieee_is_finite(chart%metric_weight))
      CPASSERT(chart%metric_weight > 0.0_dp)
      ALLOCATE (gram(n, n))
      gram(:, :) = MATMUL(CONJG(TRANSPOSE(coefficients)), MATMUL(metric, coefficients))
      CALL zpotrf('L', n, gram, n, info)
      valid = info == 0
      IF (.NOT. valid) RETURN
      count = 0
      DO b = 2, n
         DO a = 1, b - 1
            IF (occupations(a) == occupations(b)) CYCLE
            count = count + MERGE(2, 1, complex_orbitals)
         END DO
      END DO
      ALLOCATE (chart%base(rows, n), chart%coefficients(rows, n), chart%inverse(n, n), &
                chart%rotation(n, n), chart%work(n, n), chart%transformed(n, n), &
                chart%occupations(n), chart%coordinate(count), chart%pairs(2, count), &
                chart%component(count), chart%pivots(n))
      chart%base(:, :) = coefficients
      CALL ztrsm('R', 'L', 'C', 'N', rows, n, CMPLX(1.0_dp, 0.0_dp, dp), gram, n, chart%base, rows)
      chart%occupations(:) = occupations
      i = 0
      DO b = 2, n
         DO a = 1, b - 1
            IF (occupations(a) == occupations(b)) CYCLE
            i = i + 1
            chart%pairs(:, i) = [a, b]
            chart%component(i) = 1
            IF (complex_orbitals) THEN
               i = i + 1
               chart%pairs(:, i) = [a, b]
               chart%component(i) = 2
            END IF
         END DO
      END DO
   END SUBROUTINE skala_orbital_chart_create

! **************************************************************************************************
!> \brief Release the owned orbital chart, including its evaluated coordinates.
!> \param chart ...
! **************************************************************************************************
   SUBROUTINE skala_orbital_chart_release(chart)
      TYPE(skala_orbital_chart_type), INTENT(INOUT)      :: chart

      IF (ALLOCATED(chart%base)) THEN
         DEALLOCATE (chart%base, chart%coefficients, chart%inverse, chart%rotation, &
                     chart%work, chart%transformed, chart%occupations, chart%coordinate, &
                     chart%pairs, chart%component, chart%pivots)
      END IF
      chart%evaluated = .FALSE.
      chart%metric_weight = 1.0_dp
   END SUBROUTINE skala_orbital_chart_release

! **************************************************************************************************
!> \brief Number of independent real orbital coordinates at fixed occupations.
!> \param chart ...
!> \return ...
! **************************************************************************************************
   FUNCTION skala_orbital_chart_size(chart) RESULT(d)
      TYPE(skala_orbital_chart_type), INTENT(IN)         :: chart
      INTEGER                                            :: d

      d = 0
      IF (ALLOCATED(chart%coordinate)) d = SIZE(chart%coordinate)
   END FUNCTION skala_orbital_chart_size

! **************************************************************************************************
!> \brief Evaluate metric-normalized orbitals with a Cayley rotation.
!>        For the anti-Hermitian generator A, use L=(I-A/2)^(-1) and U=2L-I.
!>        The original metric is preserved for every finite coordinate vector. Repeated evaluation
!>        at identical coordinates reuses owned factors. No variable occupation or entropy update
!>        is performed by this fixed-occupation chart.
!> \param chart ...
!> \param coordinate ...
!> \param coefficients ...
! **************************************************************************************************
   SUBROUTINE skala_orbital_chart_evaluate(chart, coordinate, coefficients)
      TYPE(skala_orbital_chart_type), INTENT(INOUT)      :: chart
      REAL(dp), INTENT(IN)                               :: coordinate(:)
      COMPLEX(dp), INTENT(OUT)                           :: coefficients(:, :)

      COMPLEX(dp)                                        :: element
      INTEGER                                            :: a, b, i, info, n

      CPASSERT(ALLOCATED(chart%base))
      CPASSERT(SIZE(coordinate) == SIZE(chart%coordinate))
      CPASSERT(ALL(ieee_is_finite(coordinate)))
      CPASSERT(SIZE(coefficients, 1) == SIZE(chart%base, 1))
      CPASSERT(SIZE(coefficients, 2) == SIZE(chart%base, 2))
      IF (chart%evaluated) THEN
         IF (ALL(coordinate == chart%coordinate)) THEN
            coefficients(:, :) = chart%coefficients
            RETURN
         END IF
      END IF
      n = SIZE(chart%base, 2)
      chart%work = CMPLX(0.0_dp, 0.0_dp, dp)
      chart%inverse = CMPLX(0.0_dp, 0.0_dp, dp)
      DO i = 1, n
         chart%work(i, i) = 1.0_dp
         chart%inverse(i, i) = 1.0_dp
      END DO
      DO i = 1, SIZE(coordinate)
         a = chart%pairs(1, i)
         b = chart%pairs(2, i)
         element = CMPLX(coordinate(i), 0.0_dp, dp)
         IF (chart%component(i) == 2) element = CMPLX(0.0_dp, coordinate(i), dp)
         chart%work(a, b) = chart%work(a, b) - 0.5_dp*element
         chart%work(b, a) = chart%work(b, a) + 0.5_dp*CONJG(element)
      END DO
      CALL zgesv(n, n, chart%work, n, chart%pivots, chart%inverse, n, info)
      CPASSERT(info == 0)
      chart%rotation(:, :) = 2.0_dp*chart%inverse
      DO i = 1, n
         chart%rotation(i, i) = chart%rotation(i, i) - 1.0_dp
      END DO
      chart%coefficients(:, :) = MATMUL(chart%base, chart%rotation)
      chart%coordinate(:) = coordinate
      chart%evaluated = .TRUE.
      coefficients(:, :) = chart%coefficients
   END SUBROUTINE skala_orbital_chart_evaluate

! **************************************************************************************************
!> \brief Differentiate an original Hermitian AO expectation in the evaluated orbital chart.
!>        The coordinate derivative includes both factors of the Cayley inverse. Alternatively
!>        return the covector in orthonormal weighted occupied-orbital tangent directions at the
!>        current orbitals. This representation supplies physical cusp normals and convex weights
!>        without constructing a coordinate Gram matrix. The caller includes the k-point weight
!>        and any spin-density convention in the operator for the quantity being differentiated.
!> \param chart ...
!> \param operator ...
!> \param derivative ...
!> \param physical_tangent Whether to return the current orthonormal tangent covector.
! **************************************************************************************************
   SUBROUTINE skala_orbital_chart_derivative(chart, OPERATOR, derivative, physical_tangent)
      TYPE(skala_orbital_chart_type), INTENT(INOUT)      :: chart
      COMPLEX(dp), INTENT(IN)                            :: operator(:, :)
      REAL(dp), INTENT(OUT)                              :: derivative(:)
      LOGICAL, INTENT(IN), OPTIONAL                      :: physical_tangent

      INTEGER                                            :: a, b, i, n
      LOGICAL                                            :: physical
      REAL(dp)                                           :: value

      CPASSERT(chart%evaluated)
      CPASSERT(SIZE(OPERATOR, 1) == SIZE(chart%base, 1))
      CPASSERT(SIZE(OPERATOR, 2) == SIZE(chart%base, 1))
      CPASSERT(SIZE(derivative) == SIZE(chart%coordinate))
      physical = .FALSE.
      IF (PRESENT(physical_tangent)) physical = physical_tangent
      n = SIZE(chart%base, 2)
      IF (physical) THEN
         chart%transformed(:, :) = MATMUL(CONJG(TRANSPOSE(chart%coefficients)), &
                                          MATMUL(OPERATOR, chart%coefficients))
      ELSE
         chart%transformed(:, :) = MATMUL(CONJG(TRANSPOSE(chart%base)), MATMUL(OPERATOR, chart%base))
         chart%work(:, :) = MATMUL(CONJG(TRANSPOSE(chart%rotation)), chart%transformed)
         DO i = 1, n
            chart%work(i, :) = chart%occupations(i)*chart%work(i, :)
         END DO
         chart%transformed(:, :) = MATMUL(chart%inverse, MATMUL(chart%work, chart%inverse))
      END IF
      DO i = 1, SIZE(derivative)
         a = chart%pairs(1, i)
         b = chart%pairs(2, i)
         IF (physical) THEN
            IF (chart%component(i) == 1) THEN
               value = REAL(chart%transformed(a, b), dp)
            ELSE
               value = AIMAG(chart%transformed(a, b))
            END IF
            derivative(i) = 2.0_dp*(chart%occupations(b) - chart%occupations(a))*value/ &
                            SQRT(chart%metric_weight*(chart%occupations(a) + chart%occupations(b)))
         ELSE IF (chart%component(i) == 1) THEN
            derivative(i) = 2.0_dp*REAL(chart%transformed(b, a) - chart%transformed(a, b), dp)
         ELSE
            derivative(i) = -2.0_dp*AIMAG(chart%transformed(b, a) + chart%transformed(a, b))
         END IF
      END DO
   END SUBROUTINE skala_orbital_chart_derivative

! **************************************************************************************************
!> \brief Solve tangential orbital stationarity while reevaluating the physical constraints.
!>        Coordinates describe a differentiable chart of normalized orbitals. QR selects dependent
!>        coordinates, and physical Newton corrections restore every supplied constraint. The
!>        complete energy derivative supplies the implicit free-coordinate derivative. A bounded
!>        limited-memory BFGS iteration retains linear storage in the number of free coordinates.
!>        Convergence does not establish an admissible convex cusp or a physical minimum.
!> \param evaluate_constraints Original physical values and their coordinate derivatives.
!> \param evaluate_gradient Original energy and a consistently chosen limiting coordinate gradient.
!> \param context Owned native orbital and field evaluation data.
!> \param coordinate Initial chart coordinates, replaced only on convergence.
!> \param density_tolerances Absolute physical tolerance for each independent constraint.
!> \param gradient_tolerance Absolute tolerance for the norm of the implicit free derivative.
!> \param converged Whether physical constraints and tangential stationarity meet their tolerances.
!> \param max_iterations Maximum BFGS iterations, defaults to 200.
!> \param max_step Maximum free-coordinate step norm, defaults to 0.1.
! **************************************************************************************************
   SUBROUTINE skala_constrained_stationarity(evaluate_constraints, evaluate_gradient, context, coordinate, &
                                             density_tolerances, gradient_tolerance, converged, &
                                             max_iterations, max_step)
      PROCEDURE(physical_constraints) :: evaluate_constraints
      PROCEDURE(orbital_energy_gradient) :: evaluate_gradient
      CLASS(*), INTENT(INOUT) :: context
      REAL(dp), INTENT(INOUT) :: coordinate(:)
      REAL(dp), INTENT(IN) :: density_tolerances(:), gradient_tolerance
      LOGICAL, INTENT(OUT) :: converged
      INTEGER, INTENT(IN), OPTIONAL :: max_iterations
      REAL(dp), INTENT(IN), OPTIONAL :: max_step

      INTEGER :: d, history, i, iteration, j, limit, m, nf, shortening
      INTEGER, ALLOCATABLE :: dependent(:), free(:), pivots(:)
      LOGICAL :: accepted, valid
      REAL(dp) :: curvature, energy, fraction, radius, slope, trial_energy
      REAL(dp), ALLOCATABLE :: direction(:), g(:), jacobian(:, :), motion(:), &
                               old_gradient(:), reduced(:), scales(:), steps(:, :), trial(:), trial_reduced(:), &
                               values(:), x(:), changes(:, :)

      d = SIZE(coordinate)
      m = SIZE(density_tolerances)
      nf = d - m
      CPASSERT(d > 0 .AND. m <= d)
      CPASSERT(ALL(ieee_is_finite(coordinate)))
      CPASSERT(ALL(ieee_is_finite(density_tolerances)))
      CPASSERT(ALL(density_tolerances > 0.0_dp))
      CPASSERT(ieee_is_finite(gradient_tolerance))
      CPASSERT(gradient_tolerance > 0.0_dp)
      limit = 200
      IF (PRESENT(max_iterations)) limit = max_iterations
      radius = 0.1_dp
      IF (PRESENT(max_step)) radius = max_step
      CPASSERT(limit > 0)
      CPASSERT(ieee_is_finite(radius) .AND. radius > 0.0_dp)
      ALLOCATE (x(d), trial(d), values(m), jacobian(m, d), scales(m), dependent(m), free(nf), &
                g(d), reduced(nf), trial_reduced(nf), direction(nf), motion(nf), &
                old_gradient(nf), steps(nf, 7), changes(nf, 7), pivots(d))
      x(:) = coordinate
      CALL constraints(x)
      converged = .FALSE.
      IF (m > 0) THEN
         scales(:) = MAXVAL(ABS(jacobian), DIM=2)
         IF (ANY(scales <= TINY(0.0_dp))) RETURN
         CALL constraint_pivots(jacobian, scales, pivots, valid)
         IF (.NOT. valid) RETURN
         dependent(:) = pivots(:m)
      END IF
      j = 0
      DO i = 1, d
         IF (ANY(dependent == i)) CYCLE
         j = j + 1
         free(j) = i
      END DO
      CALL reduced_evaluation(x, energy, reduced, valid)
      IF (.NOT. valid) RETURN
      history = 0
      DO iteration = 1, limit
         IF (NORM2(reduced) <= gradient_tolerance) THEN
            coordinate(:) = x
            converged = .TRUE.
            RETURN
         END IF
         CALL bfgs_direction(reduced, steps, changes, history, direction)
         IF (DOT_PRODUCT(reduced, direction) >= 0.0_dp) THEN
            history = 0
            direction(:) = -reduced
         END IF
         IF (NORM2(direction) > radius) direction(:) = direction*(radius/NORM2(direction))
         slope = DOT_PRODUCT(reduced, direction)
         accepted = .FALSE.
         fraction = 1.0_dp
         DO shortening = 1, 24
            trial(:) = x
            trial(free) = x(free) + fraction*direction
            CALL reduced_evaluation(trial, trial_energy, trial_reduced, valid)
            IF (valid) THEN
               IF (NORM2(trial_reduced) < NORM2(reduced)) THEN
                  IF (trial_energy <= energy + 1.0e-4_dp*fraction*slope .OR. &
                      ABS(trial_energy - energy) <= 64.0_dp*EPSILON(0.0_dp)* &
                      MAX(ABS(energy), ABS(trial_energy), TINY(0.0_dp))) THEN
                     accepted = .TRUE.
                     EXIT
                  END IF
               END IF
            END IF
            fraction = 0.5_dp*fraction
         END DO
         IF (.NOT. accepted) EXIT
         motion(:) = trial(free) - x(free)
         old_gradient(:) = trial_reduced - reduced
         curvature = DOT_PRODUCT(motion, old_gradient)
         IF (curvature > 128.0_dp*EPSILON(0.0_dp)*NORM2(motion)*NORM2(old_gradient)) THEN
            IF (history == 7) THEN
               steps(:, 1:6) = steps(:, 2:7)
               changes(:, 1:6) = changes(:, 2:7)
               history = 6
            END IF
            history = history + 1
            steps(:, history) = motion
            changes(:, history) = old_gradient
         ELSE
            history = 0
         END IF
         x(:) = trial
         reduced(:) = trial_reduced
         energy = trial_energy
      END DO
      IF (NORM2(reduced) <= gradient_tolerance) THEN
         coordinate(:) = x
         converged = .TRUE.
      END IF
      ! Trial callbacks may update native caches; the caller's unchanged coordinates remain authoritative.
   CONTAINS
      SUBROUTINE constraints(c)
      REAL(dp), INTENT(IN)                               :: c(:)

         CALL evaluate_constraints(c, values, jacobian, context)
         CPASSERT(ALL(ieee_is_finite(values)))
         CPASSERT(ALL(ieee_is_finite(jacobian)))
      END SUBROUTINE constraints

      SUBROUTINE reduced_evaluation(c, e, rg, success)
      REAL(dp), INTENT(INOUT)                            :: c(:)
      REAL(dp), INTENT(OUT)                              :: e, rg(:)
      LOGICAL, INTENT(OUT)                               :: success

      INTEGER                                            :: correction, k
      LOGICAL                                            :: solved
      REAL(dp)                                           :: before, step_fraction
      REAL(dp), ALLOCATABLE                              :: candidate(:), matrix(:, :), &
                                                            multiplier(:, :), rhs(:, :)

         success = .FALSE.
         CALL constraints(c)
         IF (m > 0) THEN
            ALLOCATE (matrix(m, m), rhs(m, 1), candidate(d), multiplier(m, 1))
            DO correction = 1, 32
               IF (ALL(ABS(values) <= density_tolerances)) EXIT
               before = MAXVAL(ABS(values/scales))
               matrix(:, :) = jacobian(:, dependent)/SPREAD(scales, 2, m)
               rhs(:, 1) = -values/scales
               CALL solve_constraint_system(matrix, rhs, solved)
               IF (.NOT. solved) RETURN
               step_fraction = 1.0_dp
               DO k = 1, 24
                  candidate(:) = c
                  candidate(dependent) = c(dependent) + step_fraction*rhs(:, 1)
                  CALL constraints(candidate)
                  IF (ALL(ABS(values) <= density_tolerances)) EXIT
                  IF (MAXVAL(ABS(values/scales)) < before) EXIT
                  step_fraction = 0.5_dp*step_fraction
               END DO
               IF (k > 24) RETURN
               c(:) = candidate
            END DO
            IF (.NOT. ALL(ABS(values) <= density_tolerances)) RETURN
         END IF
         CALL evaluate_gradient(c, e, g, context)
         CPASSERT(ieee_is_finite(e))
         CPASSERT(ALL(ieee_is_finite(g)))
         rg(:) = g(free)
         IF (m > 0) THEN
            matrix(:, :) = TRANSPOSE(jacobian(:, dependent)/SPREAD(scales, 2, m))
            multiplier(:, 1) = g(dependent)
            CALL solve_constraint_system(matrix, multiplier, solved)
            IF (.NOT. solved) RETURN
            rg(:) = rg - MATMUL(TRANSPOSE(jacobian(:, free)/SPREAD(scales, 2, nf)), multiplier(:, 1))
         END IF
         success = .TRUE.
      END SUBROUTINE reduced_evaluation
   END SUBROUTINE skala_constrained_stationarity

! **************************************************************************************************
!> \brief Check the two limiting derivatives along each independent physical constraint normal.
!>        Gradients and constraint Jacobians use the same orthonormal real orbital tangent chart.
!>        Each negative change switches that sample alone from the all-positive branch. Dual
!>        normals retain the other constraints, so their density units do not affect the signs.
!>        This is a normal admissibility check, not proof of tangential stability or stationarity.
!> \param positive_gradient Complete all-positive limiting orbital derivative.
!> \param negative_changes Independent negative-minus-positive changes, one column per constraint.
!> \param jacobian Physical constraint derivatives, one row per constraint.
!> \param admissible Whether a convex cusp brackets zero along each positive physical normal.
!> \param release_side Minus one for a negative-side descent, plus one for a positive-side descent.
!> \param valid Whether the physical constraint normals are independent.
! **************************************************************************************************
   SUBROUTINE skala_cusp_admissibility(positive_gradient, negative_changes, jacobian, &
                                       admissible, release_side, valid)
      REAL(dp), INTENT(IN) :: positive_gradient(:), negative_changes(:, :), jacobian(:, :)
      LOGICAL, INTENT(OUT) :: admissible(:), valid
      INTEGER, INTENT(OUT) :: release_side(:)
      INTEGER :: d, i, info, lwork, m, rank
      REAL(dp) :: negative_slope, positive_slope, query(1)
      REAL(dp), ALLOCATABLE :: matrix(:, :), normals(:, :), scales(:), singular(:), work(:)

      d = SIZE(positive_gradient)
      m = SIZE(jacobian, 1)
      CPASSERT(m <= d .AND. SIZE(jacobian, 2) == d)
      CPASSERT(SIZE(negative_changes, 1) == d)
      CPASSERT(SIZE(negative_changes, 2) == m)
      CPASSERT(SIZE(admissible) == m .AND. SIZE(release_side) == m)
      CPASSERT(ALL(ieee_is_finite(positive_gradient)))
      CPASSERT(ALL(ieee_is_finite(negative_changes)))
      CPASSERT(ALL(ieee_is_finite(jacobian)))
      admissible = .FALSE.
      release_side = 0
      valid = .TRUE.
      IF (m == 0) RETURN
      valid = .FALSE.
      ALLOCATE (scales(m), matrix(m, d), normals(d, m), singular(m))
      scales(:) = MAXVAL(ABS(jacobian), DIM=2)
      IF (ANY(scales <= TINY(0.0_dp))) RETURN
      matrix(:, :) = jacobian/SPREAD(scales, 2, d)
      normals = 0.0_dp
      DO i = 1, m
         normals(i, i) = 1.0_dp
      END DO
      CALL dgelss(m, d, m, matrix, m, normals, d, singular, EPSILON(0.0_dp), rank, query, -1, info)
      IF (info /= 0) RETURN
      lwork = MAX(1, CEILING(query(1)))
      ALLOCATE (work(lwork))
      CALL dgelss(m, d, m, matrix, m, normals, d, singular, EPSILON(0.0_dp), rank, work, lwork, info)
      IF (info /= 0 .OR. rank /= m) RETURN
      DO i = 1, m
         normals(:, i) = normals(:, i)/NORM2(normals(:, i))
         positive_slope = DOT_PRODUCT(positive_gradient, normals(:, i))
         negative_slope = positive_slope + DOT_PRODUCT(negative_changes(:, i), normals(:, i))
         admissible(i) = negative_slope <= 0.0_dp .AND. positive_slope >= 0.0_dp
         IF (admissible(i)) CYCLE
         IF (-positive_slope >= negative_slope .AND. positive_slope < 0.0_dp) THEN
            release_side(i) = 1
         ELSE IF (negative_slope > 0.0_dp) THEN
            release_side(i) = -1
         END IF
      END DO
      valid = .TRUE.
   END SUBROUTINE skala_cusp_admissibility

! **************************************************************************************************
!> \brief Select dependent coordinates from the scaled physical constraint Jacobian.
!> \param jacobian ...
!> \param scales ...
!> \param pivots ...
!> \param success ...
! **************************************************************************************************
   SUBROUTINE constraint_pivots(jacobian, scales, pivots, success)
      REAL(dp), INTENT(IN)                               :: jacobian(:, :), scales(:)
      INTEGER, INTENT(OUT)                               :: pivots(:)
      LOGICAL, INTENT(OUT)                               :: success

      INTEGER                                            :: info, lwork, m, n
      REAL(dp)                                           :: query(1)
      REAL(dp), ALLOCATABLE                              :: matrix(:, :), tau(:), work(:)

      m = SIZE(jacobian, 1)
      n = SIZE(jacobian, 2)
      ALLOCATE (matrix(m, n), tau(m))
      matrix(:, :) = jacobian/SPREAD(scales, 2, n)
      pivots = 0
      CALL dgeqp3(m, n, matrix, m, pivots, tau, query, -1, info)
      IF (info /= 0) THEN
         success = .FALSE.
         RETURN
      END IF
      lwork = MAX(1, CEILING(query(1)))
      ALLOCATE (work(lwork))
      CALL dgeqp3(m, n, matrix, m, pivots, tau, work, lwork, info)
      success = info == 0
   END SUBROUTINE constraint_pivots

! **************************************************************************************************
!> \brief Solve a full-rank dependent-coordinate system without its Gram matrix.
!> \param matrix ...
!> \param rhs ...
!> \param success ...
! **************************************************************************************************
   SUBROUTINE solve_constraint_system(matrix, rhs, success)
      REAL(dp), INTENT(INOUT)                            :: matrix(:, :), rhs(:, :)
      LOGICAL, INTENT(OUT)                               :: success

      INTEGER                                            :: info, lwork, n, rank
      REAL(dp)                                           :: query(1)
      REAL(dp), ALLOCATABLE                              :: singular(:), work(:)

      n = SIZE(matrix, 1)
      ALLOCATE (singular(n))
      CALL dgelss(n, n, SIZE(rhs, 2), matrix, n, rhs, n, singular, EPSILON(0.0_dp), rank, query, -1, info)
      IF (info /= 0) THEN
         success = .FALSE.
         RETURN
      END IF
      lwork = MAX(1, CEILING(query(1)))
      ALLOCATE (work(lwork))
      CALL dgelss(n, n, SIZE(rhs, 2), matrix, n, rhs, n, singular, EPSILON(0.0_dp), rank, work, lwork, info)
      success = info == 0 .AND. rank == n .AND. ALL(ieee_is_finite(rhs))
   END SUBROUTINE solve_constraint_system

! **************************************************************************************************
!> \brief Apply the inverse limited-memory BFGS operator to the free orbital derivative.
!> \param gradient ...
!> \param steps ...
!> \param changes ...
!> \param history ...
!> \param direction ...
! **************************************************************************************************
   SUBROUTINE bfgs_direction(gradient, steps, changes, history, direction)
      REAL(dp), INTENT(IN)                               :: gradient(:), steps(:, :), changes(:, :)
      INTEGER, INTENT(IN)                                :: history
      REAL(dp), INTENT(OUT)                              :: direction(:)

      INTEGER                                            :: i
      REAL(dp)                                           :: alpha(7), beta, inverse_curvature(7), &
                                                            scale

      direction(:) = gradient
      DO i = history, 1, -1
         inverse_curvature(i) = 1.0_dp/DOT_PRODUCT(steps(:, i), changes(:, i))
         alpha(i) = inverse_curvature(i)*DOT_PRODUCT(steps(:, i), direction)
         direction(:) = direction - alpha(i)*changes(:, i)
      END DO
      scale = 1.0_dp
      IF (history > 0) scale = DOT_PRODUCT(steps(:, history), changes(:, history))/ &
                               DOT_PRODUCT(changes(:, history), changes(:, history))
      direction(:) = scale*direction
      DO i = 1, history
         beta = inverse_curvature(i)*DOT_PRODUCT(changes(:, i), direction)
         direction(:) = direction + steps(:, i)*(alpha(i) - beta)
      END DO
      direction(:) = -direction
   END SUBROUTINE bfgs_direction
END MODULE skala_orbital_constraints
