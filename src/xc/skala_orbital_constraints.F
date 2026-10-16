!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Solve orbital stationarity with physical descriptor constraints.
! **************************************************************************************************
MODULE skala_orbital_constraints
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: dp
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'skala_orbital_constraints'

   PUBLIC :: skala_constrained_stationarity, skala_cusp_admissibility
   PUBLIC :: skala_constraint_chart_type, skala_constraint_chart_create, skala_constraint_chart_restore
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

   TYPE skala_constraint_chart_type
      PRIVATE
      INTEGER, ALLOCATABLE :: rows(:), dependent(:), free(:)
      REAL(dp), ALLOCATABLE :: scales(:), tolerances(:)
   END TYPE skala_constraint_chart_type

   ABSTRACT INTERFACE
! **************************************************************************************************
!> \brief Evaluate physical constraints and their derivatives in orbital coordinates.
!> \param coordinate Complete orbital coordinates.
!> \param values Physical constraint residuals.
!> \param jacobian Derivative of each residual with respect to the coordinates.
!> \param context Borrowed problem state.
! **************************************************************************************************
      SUBROUTINE physical_constraints(coordinate, values, jacobian, context)
         IMPORT :: dp
         REAL(dp), INTENT(IN) :: coordinate(:)
         REAL(dp), INTENT(OUT) :: values(:), jacobian(:, :)
         CLASS(*), INTENT(INOUT) :: context
      END SUBROUTINE physical_constraints
! **************************************************************************************************
!> \brief Evaluate the complete energy and its orbital derivative.
!> \param coordinate Complete orbital coordinates.
!> \param energy Energy of the trial state.
!> \param gradient Energy derivative in the same coordinates.
!> \param context Borrowed problem state.
! **************************************************************************************************
      SUBROUTINE orbital_energy_gradient(coordinate, energy, gradient, context)
         IMPORT :: dp
         REAL(dp), INTENT(IN) :: coordinate(:)
         REAL(dp), INTENT(OUT) :: energy, gradient(:)
         CLASS(*), INTENT(INOUT) :: context
      END SUBROUTINE orbital_energy_gradient
! **************************************************************************************************
!> \brief Verify a stationary trial against an additional physical convergence criterion.
!> \param coordinate Complete orbital coordinates.
!> \param context Borrowed problem state.
!> \param converged Whether the physical criterion is satisfied.
! **************************************************************************************************
      SUBROUTINE orbital_stationarity(coordinate, context, converged)
         IMPORT :: dp
         REAL(dp), INTENT(IN) :: coordinate(:)
         CLASS(*), INTENT(INOUT) :: context
         LOGICAL, INTENT(OUT) :: converged
      END SUBROUTINE orbital_stationarity
! **************************************************************************************************
!> \brief Poll the caller's external stop condition between orbital trials.
!> \param context Borrowed problem state.
!> \param requested Whether the caller requests termination.
! **************************************************************************************************
      SUBROUTINE orbital_stop(context, requested)
         CLASS(*), INTENT(INOUT) :: context
         LOGICAL, INTENT(OUT) :: requested
      END SUBROUTINE orbital_stop
! **************************************************************************************************
!> \brief Check a trial segment before accepting a minimization step.
!> \param first Last accepted feasible coordinates.
!> \param last Proposed feasible coordinates.
!> \param curve Fixed independent coordinates and restoration of the original physical constraints.
!> \param context Borrowed problem state, including any discovered event.
!> \param interrupted Whether to return the accepted state for caller event handling.
! **************************************************************************************************
      SUBROUTINE orbital_segment(first, last, curve, context, interrupted)
         IMPORT :: dp, skala_constraint_chart_type
         REAL(dp), INTENT(IN) :: first(:), last(:)
         TYPE(skala_constraint_chart_type), INTENT(IN) :: curve
         CLASS(*), INTENT(INOUT) :: context
         LOGICAL, INTENT(OUT) :: interrupted
      END SUBROUTINE orbital_segment
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
!> \param valid Whether the requested orbital space has a positive charge metric.
!> \param metric_weight Positive k-point weight for the physical occupied-orbital tangent norm.
!> \param complete_space Append the missing zero-occupation directions of a positive charge metric.
!>        By default only the supplied orbital span is retained.
! **************************************************************************************************
   SUBROUTINE skala_orbital_chart_create(coefficients, metric, occupations, complex_orbitals, chart, valid, &
                                         metric_weight, complete_space)
      COMPLEX(dp), INTENT(IN)                            :: coefficients(:, :), metric(:, :)
      REAL(dp), INTENT(IN)                               :: occupations(:)
      LOGICAL, INTENT(IN)                                :: complex_orbitals
      TYPE(skala_orbital_chart_type), INTENT(INOUT)      :: chart
      LOGICAL, INTENT(OUT)                               :: valid
      REAL(dp), INTENT(IN), OPTIONAL                     :: metric_weight
      LOGICAL, INTENT(IN), OPTIONAL                      :: complete_space

      COMPLEX(dp), ALLOCATABLE                           :: gram(:, :)
      INTEGER                                            :: a, b, count, i, info, n, nreference, rows
      REAL(dp), ALLOCATABLE                              :: full_occupations(:)

      CALL skala_orbital_chart_release(chart)
      rows = SIZE(coefficients, 1)
      nreference = SIZE(coefficients, 2)
      n = nreference
      IF (PRESENT(complete_space)) THEN
         IF (complete_space) n = rows
      END IF
      valid = rows > 0 .AND. nreference > 0 .AND. nreference <= rows
      CPASSERT(valid)
      CPASSERT(SIZE(metric, 1) == rows .AND. SIZE(metric, 2) == rows)
      CPASSERT(SIZE(occupations) == nreference)
      CPASSERT(ALL(ieee_is_finite(REAL(coefficients, dp))))
      CPASSERT(ALL(ieee_is_finite(AIMAG(coefficients))))
      CPASSERT(ALL(ieee_is_finite(REAL(metric, dp))))
      CPASSERT(ALL(ieee_is_finite(AIMAG(metric))))
      CPASSERT(ALL(ieee_is_finite(occupations)) .AND. ALL(occupations >= 0.0_dp))
      IF (PRESENT(metric_weight)) chart%metric_weight = metric_weight
      CPASSERT(ieee_is_finite(chart%metric_weight))
      CPASSERT(chart%metric_weight > 0.0_dp)
      ALLOCATE (gram(nreference, nreference), full_occupations(n))
      full_occupations = 0.0_dp
      full_occupations(:nreference) = occupations
      gram(:, :) = MATMUL(CONJG(TRANSPOSE(coefficients)), MATMUL(metric, coefficients))
      CALL zpotrf('L', nreference, gram, nreference, info)
      valid = info == 0
      IF (.NOT. valid) RETURN
      count = 0
      DO b = 2, n
         DO a = 1, b - 1
            IF (full_occupations(a) == full_occupations(b)) CYCLE
            count = count + MERGE(2, 1, complex_orbitals)
         END DO
      END DO
      ALLOCATE (chart%base(rows, n), chart%coefficients(rows, n), chart%inverse(n, n), &
                chart%rotation(n, n), chart%work(n, n), chart%transformed(n, n), &
                chart%occupations(n), chart%coordinate(count), chart%pairs(2, count), &
                chart%component(count), chart%pivots(n))
      chart%base(:, :nreference) = coefficients
      CALL ztrsm('R', 'L', 'C', 'N', rows, nreference, CMPLX(1.0_dp, 0.0_dp, dp), &
                 gram, nreference, chart%base, rows)
      IF (n > nreference) THEN
         CALL complete_metric_space(chart%base, metric, nreference, valid)
         IF (.NOT. valid) THEN
            CALL skala_orbital_chart_release(chart)
            RETURN
         END IF
      END IF
      chart%occupations(:) = full_occupations
      i = 0
      DO b = 2, n
         DO a = 1, b - 1
            IF (full_occupations(a) == full_occupations(b)) CYCLE
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
!> \brief Complete normalized reference columns in the original positive charge metric.
!>        QR in the whitened space supplies only the missing complement. The existing reference
!>        columns are retained, including their phases, and no Gaussian-overlap substitution is made.
!> \param coefficients Normalized reference columns followed by storage for their complement.
!> \param metric Original Hermitian charge operator in the AO space.
!> \param nreference Number of supplied normalized columns.
!> \param valid Whether the full charge metric is positive definite.
! **************************************************************************************************
   SUBROUTINE complete_metric_space(coefficients, metric, nreference, valid)
      COMPLEX(dp), INTENT(INOUT)                         :: coefficients(:, :)
      COMPLEX(dp), INTENT(IN)                            :: metric(:, :)
      INTEGER, INTENT(IN)                                :: nreference
      LOGICAL, INTENT(OUT)                               :: valid

      COMPLEX(dp)                                        :: query(1)
      COMPLEX(dp), ALLOCATABLE                           :: q(:, :), tau(:), upper(:, :), work(:)
      INTEGER                                            :: info, n, workspace

      n = SIZE(coefficients, 1)
      CPASSERT(SIZE(coefficients, 2) == n)
      CPASSERT(nreference > 0 .AND. nreference < n)
      ALLOCATE (upper(n, n), q(n, n), tau(nreference))
      upper(:, :) = metric
      CALL zpotrf('U', n, upper, n, info)
      valid = info == 0
      IF (.NOT. valid) RETURN
      q = CMPLX(0.0_dp, 0.0_dp, dp)
      q(:, :nreference) = coefficients(:, :nreference)
      CALL ztrmm('L', 'U', 'N', 'N', n, nreference, CMPLX(1.0_dp, 0.0_dp, dp), upper, n, q, n)
      CALL zgeqrf(n, nreference, q, n, tau, query, -1, info)
      CPASSERT(info == 0)
      workspace = MAX(1, INT(REAL(query(1), dp)))
      ALLOCATE (work(workspace))
      CALL zgeqrf(n, nreference, q, n, tau, work, workspace, info)
      CPASSERT(info == 0)
      DEALLOCATE (work)
      CALL zungqr(n, n, nreference, q, n, tau, query, -1, info)
      CPASSERT(info == 0)
      workspace = MAX(1, INT(REAL(query(1), dp)))
      ALLOCATE (work(workspace))
      CALL zungqr(n, n, nreference, q, n, tau, work, workspace, info)
      CPASSERT(info == 0)
      coefficients(:, nreference + 1:n) = q(:, nreference + 1:n)
      CALL ztrsm('L', 'U', 'N', 'N', n, n - nreference, CMPLX(1.0_dp, 0.0_dp, dp), &
                 upper, n, coefficients(:, nreference + 1:n), n)
   END SUBROUTINE complete_metric_space

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
!> \param OPERATOR ...
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
!> \param coordinate Initial coordinates, replaced on convergence, a stop or a segment event.
!>        A stop or event returns the last accepted feasible trial. Other failures retain the initial state
!>        unless retain_progress is requested.
!> \param density_tolerances Absolute physical tolerance for every original constraint.
!> \param gradient_tolerance Absolute tolerance for the norm of the implicit free derivative.
!> \param converged Whether physical constraints and tangential stationarity meet their tolerances.
!> \param max_iterations Maximum BFGS iterations, defaults to 200.
!> \param max_step Maximum free-coordinate step norm, defaults to 0.1.
!> \param verify_stationarity Optional physical check, additional to the free-coordinate norm.
!> \param iterations Optional number of attempted BFGS iterations, bounded by max_iterations.
!> \param check_stop Optional external stop check between iterations and line-search trials.
!> \param stopped Whether an external stop was requested, independent of convergence.
!> \param check_segment Optional trial-segment check before step acceptance.
!> \param interrupted Whether that check requests event handling at the last accepted state.
!> \param retain_progress Return the last accepted feasible coordinates on failure, default false.
! **************************************************************************************************
   SUBROUTINE skala_constrained_stationarity(evaluate_constraints, evaluate_gradient, context, coordinate, &
                                             density_tolerances, gradient_tolerance, converged, &
                                             max_iterations, max_step, verify_stationarity, iterations, check_stop, stopped, &
                                             check_segment, interrupted, retain_progress)
      PROCEDURE(physical_constraints) :: evaluate_constraints
      PROCEDURE(orbital_energy_gradient) :: evaluate_gradient
      PROCEDURE(orbital_stationarity), OPTIONAL :: verify_stationarity
      PROCEDURE(orbital_stop), OPTIONAL :: check_stop
      PROCEDURE(orbital_segment), OPTIONAL :: check_segment
      CLASS(*), INTENT(INOUT) :: context
      REAL(dp), INTENT(INOUT) :: coordinate(:)
      REAL(dp), INTENT(IN) :: density_tolerances(:), gradient_tolerance
      LOGICAL, INTENT(OUT) :: converged
      INTEGER, INTENT(IN), OPTIONAL :: max_iterations
      REAL(dp), INTENT(IN), OPTIONAL :: max_step
      INTEGER, INTENT(OUT), OPTIONAL :: iterations
      LOGICAL, INTENT(OUT), OPTIONAL :: stopped
      LOGICAL, INTENT(OUT), OPTIONAL :: interrupted
      LOGICAL, INTENT(IN), OPTIONAL :: retain_progress

      INTEGER :: capacity, d, history, i, iteration, limit, m, nf, r, shortening
      TYPE(skala_constraint_chart_type) :: curve
      LOGICAL :: accepted, valid, event, newton_valid, newton_stopped
      REAL(dp) :: curvature, energy, fraction, radius, roundoff, slope, trial_energy, trial_slope
      REAL(dp), ALLOCATABLE :: direction(:), g(:), jacobian(:, :), motion(:), newton(:), &
                               old_gradient(:), reduced(:), steps(:, :), trial(:), trial_reduced(:), &
                               values(:), x(:), changes(:, :)

      d = SIZE(coordinate)
      IF (PRESENT(iterations)) iterations = 0
      IF (PRESENT(stopped)) stopped = .FALSE.
      IF (PRESENT(interrupted)) interrupted = .FALSE.
      m = SIZE(density_tolerances)
      CPASSERT(d > 0)
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
      ALLOCATE (x(d), trial(d), values(m), jacobian(m, d))
      x(:) = coordinate
      converged = .FALSE.
      CALL skala_constraint_chart_create(evaluate_constraints, context, x, density_tolerances, curve, valid)
      IF (.NOT. valid) RETURN
      r = SIZE(curve%dependent)
      nf = SIZE(curve%free)
      ! Use up to 256 KiB for curvature histories in small orbital spaces. Large
      ! spaces retain 32 entries and linear storage.
      capacity = MAX(7, MIN(nf, MAX(32, 16384/MAX(1, nf))))
      ALLOCATE (g(d), reduced(nf), trial_reduced(nf), direction(nf), motion(nf), &
                old_gradient(nf), steps(nf, capacity), changes(nf, capacity), newton(nf))
      CALL reduced_evaluation(x, energy, reduced, valid)
      IF (.NOT. valid) RETURN
      history = 0
      DO iteration = 1, limit
         IF (stop_requested()) RETURN
         IF (stationary(x, reduced)) THEN
            coordinate(:) = x
            converged = .TRUE.
            RETURN
         END IF
         IF (PRESENT(iterations)) iterations = iteration
         CALL bfgs_direction(reduced, steps, changes, history, direction)
         IF (DOT_PRODUCT(reduced, direction) >= 0.0_dp) THEN
            history = 0
            direction(:) = -reduced
         END IF
         roundoff = 64.0_dp*EPSILON(0.0_dp)*MAX(ABS(energy), TINY(0.0_dp))
         IF (history > 0 .AND. &
             -DOT_PRODUCT(reduced, direction) <= 32.0_dp*roundoff) THEN
            ! Resolve curvature below the energy resolution from reduced derivatives.
            ! Each perturbed trial restores the same physical constraints, retaining
            ! the curvature of the eliminated coordinates.
            CALL newton_direction(x, reduced, newton, newton_valid, newton_stopped)
            IF (newton_stopped) RETURN
            IF (newton_valid) direction(:) = newton
         END IF
         IF (NORM2(direction) > radius) direction(:) = direction*(radius/NORM2(direction))
         slope = DOT_PRODUCT(reduced, direction)
         accepted = .FALSE.
         fraction = 1.0_dp
         DO shortening = 1, DIGITS(1.0_dp)
            IF (stop_requested()) RETURN
            trial(:) = x
            trial(curve%free) = x(curve%free) + fraction*direction
            IF (ALL(trial(curve%free) == x(curve%free))) EXIT
            CALL reduced_evaluation(trial, trial_energy, trial_reduced, valid)
            IF (valid) THEN
               event = .FALSE.
               IF (PRESENT(check_segment)) CALL check_segment(x, trial, curve, context, event)
               IF (event) THEN
                  coordinate(:) = x
                  IF (PRESENT(interrupted)) interrupted = .TRUE.
                  RETURN
               END IF
               ! Descent may increase the gradient norm in a region of negative curvature.
               roundoff = 64.0_dp*EPSILON(0.0_dp)*MAX(ABS(energy), ABS(trial_energy), TINY(0.0_dp))
               trial_slope = DOT_PRODUCT(trial_reduced, direction)
               IF ((trial_energy < energy - roundoff .AND. &
                    trial_energy <= energy + 1.0e-4_dp*fraction*slope) .OR. &
                   (ABS(trial_energy - energy) <= roundoff .AND. &
                    (NORM2(trial_reduced) < NORM2(reduced) .OR. &
                     ABS(trial_slope) <= 0.9_dp*ABS(slope)))) THEN
                  ! Below energy resolution, directional curvature can establish a useful
                  ! step even when other gradient components grow in an anisotropic valley.
                  accepted = .TRUE.
                  EXIT
               END IF
            END IF
            fraction = 0.5_dp*fraction
         END DO
         IF (.NOT. accepted) THEN
            IF (history == 0) EXIT
            ! Retry with steepest descent when the stored curvature gives no acceptable step.
            history = 0
            CYCLE
         END IF
         motion(:) = trial(curve%free) - x(curve%free)
         old_gradient(:) = trial_reduced - reduced
         curvature = DOT_PRODUCT(motion, old_gradient)
         IF (curvature > 128.0_dp*EPSILON(0.0_dp)*NORM2(motion)*NORM2(old_gradient)) THEN
            IF (history == capacity) THEN
               steps(:, 1:capacity - 1) = steps(:, 2:capacity)
               changes(:, 1:capacity - 1) = changes(:, 2:capacity)
               history = capacity - 1
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
      IF (stop_requested()) RETURN
      IF (stationary(x, reduced)) THEN
         coordinate(:) = x
         converged = .TRUE.
      END IF
      IF (PRESENT(retain_progress)) THEN
         IF (retain_progress) coordinate(:) = x
      END IF
      ! Trial callbacks may update native caches; the returned coordinates remain authoritative.
   CONTAINS
! **************************************************************************************************
!> \brief Inexact Newton direction from central differences of the native derivatives.
!>        Truncated conjugate gradients use linear storage. Indefinite curvature returns
!>        the current descent direction, or failure if none was found. The outer line search
!>        retains its original energy and derivative checks.
!> \param center Accepted orbital coordinates.
!> \param derivative Derivative at the accepted state.
!> \param newton_step Candidate free-coordinate descent direction.
!> \param valid_direction Whether the candidate is a finite descent direction.
!> \param stopped_direction Whether an external stop interrupted a derivative evaluation.
! **************************************************************************************************
      SUBROUTINE newton_direction(center, derivative, newton_step, valid_direction, stopped_direction)
      REAL(dp), INTENT(IN)                               :: center(:), derivative(:)
      REAL(dp), INTENT(OUT)                              :: newton_step(:)
      LOGICAL, INTENT(OUT)                               :: valid_direction, stopped_direction

      INTEGER                                            :: iteration
      LOGICAL                                            :: valid_trial
      REAL(dp)                                           :: alpha, beta, curvature, energy_trial, &
                                                            next_residual2, norm_direction, &
                                                            residual2, step, TARGET
      REAL(dp), ALLOCATABLE                              :: gradient_minus(:), gradient_plus(:), &
                                                            minus(:), plus(:), PRODUCT(:), &
                                                            residual(:), search(:)

         valid_direction = .FALSE.
         stopped_direction = .FALSE.
         newton_step = 0.0_dp
         ALLOCATE (plus(d), minus(d), gradient_plus(nf), gradient_minus(nf), residual(nf), search(nf), PRODUCT(nf))
         residual(:) = -derivative
         search(:) = residual
         residual2 = DOT_PRODUCT(residual, residual)
         TARGET = MIN(0.1_dp, SQRT(NORM2(derivative)))*NORM2(derivative)
         DO iteration = 1, MIN(nf, 64)
            IF (stop_requested()) THEN
               stopped_direction = .TRUE.
               RETURN
            END IF
            norm_direction = NORM2(search)
            IF (norm_direction <= TINY(0.0_dp)) EXIT
            step = SQRT(EPSILON(0.0_dp))*MAX(1.0_dp, NORM2(center(curve%free)))/norm_direction
            plus(:) = center
            minus(:) = center
            plus(curve%free) = center(curve%free) + step*search
            minus(curve%free) = center(curve%free) - step*search
            CALL reduced_evaluation(plus, energy_trial, gradient_plus, valid_trial)
            IF (.NOT. valid_trial) RETURN
            IF (stop_requested()) THEN
               stopped_direction = .TRUE.
               RETURN
            END IF
            CALL reduced_evaluation(minus, energy_trial, gradient_minus, valid_trial)
            IF (.NOT. valid_trial) RETURN
            PRODUCT(:) = (gradient_plus - gradient_minus)/(2.0_dp*step)
            curvature = DOT_PRODUCT(search, product)
            IF (curvature <= 128.0_dp*EPSILON(0.0_dp)*norm_direction*NORM2(product)) EXIT
            alpha = residual2/curvature
            newton_step(:) = newton_step + alpha*search
            residual(:) = residual - alpha*product
            next_residual2 = DOT_PRODUCT(residual, residual)
            IF (SQRT(next_residual2) <= TARGET) EXIT
            beta = next_residual2/residual2
            search(:) = residual + beta*search
            residual2 = next_residual2
         END DO
         valid_direction = ALL(ieee_is_finite(newton_step)) .AND. DOT_PRODUCT(derivative, newton_step) < 0.0_dp
      END SUBROUTINE newton_direction
! **************************************************************************************************
!> \brief Preserve the last accepted feasible coordinates when the caller requests a stop.
!> \return Whether termination was requested.
! **************************************************************************************************
      FUNCTION stop_requested() RESULT(requested)
      LOGICAL                                            :: requested

         requested = .FALSE.
         IF (PRESENT(check_stop)) CALL check_stop(context, requested)
         IF (.NOT. requested) RETURN
         coordinate(:) = x
         IF (PRESENT(stopped)) stopped = .TRUE.
      END FUNCTION stop_requested
! **************************************************************************************************
!> \brief Require both the reduced derivative and any supplied physical convergence check.
!> \param c Complete trial coordinates.
!> \param rg Reduced energy derivative at the trial.
!> \return Whether the trial satisfies both criteria.
! **************************************************************************************************
      FUNCTION stationary(c, rg) RESULT(converged_trial)
      REAL(dp), INTENT(IN)                               :: c(:), rg(:)
      LOGICAL                                            :: converged_trial

         converged_trial = NORM2(rg) <= gradient_tolerance
         IF (.NOT. converged_trial) RETURN
         IF (PRESENT(verify_stationarity)) CALL verify_stationarity(c, context, converged_trial)
      END FUNCTION stationary
! **************************************************************************************************
!> \brief Restore the physical targets and eliminate the dependent-coordinate derivatives.
!> \param c Trial coordinates, including restored dependent coordinates.
!> \param e Complete energy at the restored trial.
!> \param rg Energy derivative in the independent coordinates.
!> \param success Whether the constraints and reduced derivative were evaluated successfully.
! **************************************************************************************************
      SUBROUTINE reduced_evaluation(c, e, rg, success)
      REAL(dp), INTENT(INOUT)                            :: c(:)
      REAL(dp), INTENT(OUT)                              :: e, rg(:)
      LOGICAL, INTENT(OUT)                               :: success

      LOGICAL                                            :: solved
      REAL(dp), ALLOCATABLE                              :: matrix(:, :), multiplier(:, :)

         CALL skala_constraint_chart_restore(curve, evaluate_constraints, context, c, values, jacobian, success)
         IF (.NOT. success) RETURN
         success = .FALSE.
         CALL evaluate_gradient(c, e, g, context)
         CPASSERT(ieee_is_finite(e))
         CPASSERT(ALL(ieee_is_finite(g)))
         rg(:) = g(curve%free)
         IF (r > 0) THEN
            ALLOCATE (matrix(r, r), multiplier(r, 1))
            matrix(:, :) = TRANSPOSE(jacobian(curve%rows(:r), curve%dependent)/SPREAD(curve%scales(curve%rows(:r)), 2, r))
            multiplier(:, 1) = g(curve%dependent)
            CALL solve_constraint_system(matrix, multiplier, solved)
            IF (.NOT. solved) RETURN
  rg(:) = rg - MATMUL(TRANSPOSE(jacobian(curve%rows(:r), curve%free)/SPREAD(curve%scales(curve%rows(:r)), 2, nf)), multiplier(:, 1))
         END IF
         success = .TRUE.
      END SUBROUTINE reduced_evaluation
   END SUBROUTINE skala_constrained_stationarity

! **************************************************************************************************
!> \brief Fix independent coordinates for restoration of simultaneous physical constraints.
!>        Rank selection uses all original physical residuals. Dependent samples retain their
!>        individual feasibility tests even when fewer rows determine the restoration.
!> \param evaluate_constraints Physical residuals and their orbital-coordinate Jacobian.
!> \param context Borrowed native problem state.
!> \param coordinate Reference orbital coordinates.
!> \param tolerances Absolute feasibility window of every original physical sample.
!> \param chart Fixed independent coordinates, dependent coordinates and physical row scales.
!> \param valid Whether the independent constraint rows and coordinates could be selected.
! **************************************************************************************************
   SUBROUTINE skala_constraint_chart_create(evaluate_constraints, context, coordinate, tolerances, chart, valid)
      PROCEDURE(physical_constraints) :: evaluate_constraints
      CLASS(*), INTENT(INOUT) :: context
      REAL(dp), INTENT(IN) :: coordinate(:), tolerances(:)
      TYPE(skala_constraint_chart_type), INTENT(OUT) :: chart
      LOGICAL, INTENT(OUT) :: valid

      INTEGER :: d, i, j, m, r
      INTEGER, ALLOCATABLE :: pivots(:)
      REAL(dp), ALLOCATABLE :: values(:), jacobian(:, :)

      d = SIZE(coordinate)
      m = SIZE(tolerances)
      CPASSERT(d > 0 .AND. ALL(ieee_is_finite(coordinate)))
      CPASSERT(ALL(ieee_is_finite(tolerances)) .AND. ALL(tolerances > 0.0_dp))
      ALLOCATE (chart%tolerances(m), chart%scales(m), chart%rows(m), pivots(d), values(m), jacobian(m, d))
      chart%tolerances(:) = tolerances
      CALL evaluate_constraints(coordinate, values, jacobian, context)
      CPASSERT(ALL(ieee_is_finite(values)))
      CPASSERT(ALL(ieee_is_finite(jacobian)))
      r = 0
      valid = .TRUE.
      IF (m > 0) THEN
         chart%scales(:) = MAX(MAXVAL(ABS(jacobian), DIM=2), tolerances)
         CALL independent_constraint_rows(jacobian, chart%scales, chart%rows, r, valid)
         IF (.NOT. valid) RETURN
         IF (r > 0) THEN
            CALL constraint_pivots(jacobian(chart%rows(:r), :), chart%scales(chart%rows(:r)), pivots, valid)
            IF (.NOT. valid) RETURN
         END IF
      END IF
      ALLOCATE (chart%dependent(r), chart%free(d - r))
      chart%dependent(:) = pivots(:r)
      j = 0
      DO i = 1, d
         IF (ANY(chart%dependent == i)) CYCLE
         j = j + 1
         chart%free(j) = i
      END DO
   END SUBROUTINE skala_constraint_chart_create

! **************************************************************************************************
!> \brief Restore physical constraints with the same dependent coordinates as orbital minimization.
!>        Free coordinates remain fixed. Interpolate them between feasible endpoints and restore
!>        the dependent coordinates to follow the minimizer's implicit constraint curve. Linear
!>        interpolation of all orbital coordinates need not satisfy the physical constraints.
!> \param chart Fixed constraint chart selected at minimization initialization.
!> \param evaluate_constraints Physical residuals and orbital-coordinate Jacobian.
!> \param context Borrowed native problem state.
!> \param coordinate Trial coordinates, replaced by the restored coordinates when successful.
!> \param values Original physical residuals at the returned trial.
!> \param jacobian Original physical Jacobian at the returned trial.
!> \param valid Whether every original sample meets its own absolute feasibility window.
! **************************************************************************************************
   SUBROUTINE skala_constraint_chart_restore(chart, evaluate_constraints, context, coordinate, values, jacobian, valid)
      TYPE(skala_constraint_chart_type), INTENT(IN)      :: chart

      PROCEDURE(physical_constraints) :: evaluate_constraints
      CLASS(*), INTENT(INOUT) :: context
      REAL(dp), INTENT(INOUT) :: coordinate(:)
      REAL(dp), INTENT(OUT) :: values(:), jacobian(:, :)
      LOGICAL, INTENT(OUT) :: valid

      INTEGER :: correction, d, k, m, r
      LOGICAL :: solved
      REAL(dp) :: before, fraction
      REAL(dp), ALLOCATABLE :: candidate(:), matrix(:, :), rhs(:, :)

      d = SIZE(coordinate)
      m = SIZE(chart%tolerances)
      r = SIZE(chart%dependent)
      CPASSERT(SIZE(chart%free) + r == d)
      CPASSERT(SIZE(values) == m .AND. ALL(SHAPE(jacobian) == [m, d]))
      valid = .FALSE.
      CALL evaluate_constraints(coordinate, values, jacobian, context)
      CPASSERT(ALL(ieee_is_finite(values)))
      CPASSERT(ALL(ieee_is_finite(jacobian)))
      IF (r > 0) THEN
         ALLOCATE (matrix(r, r), rhs(r, 1), candidate(d))
         DO correction = 1, 32
            IF (ALL(ABS(values) <= chart%tolerances)) EXIT
            before = MAXVAL(ABS(values/chart%scales))
            matrix(:, :) = jacobian(chart%rows(:r), chart%dependent)/SPREAD(chart%scales(chart%rows(:r)), 2, r)
            rhs(:, 1) = -values(chart%rows(:r))/chart%scales(chart%rows(:r))
            CALL solve_constraint_system(matrix, rhs, solved)
            IF (.NOT. solved) RETURN
            fraction = 1.0_dp
            DO k = 1, 24
               candidate(:) = coordinate
               candidate(chart%dependent) = coordinate(chart%dependent) + fraction*rhs(:, 1)
               CALL evaluate_constraints(candidate, values, jacobian, context)
               CPASSERT(ALL(ieee_is_finite(values)))
               CPASSERT(ALL(ieee_is_finite(jacobian)))
               IF (ALL(ABS(values) <= chart%tolerances)) EXIT
               IF (MAXVAL(ABS(values/chart%scales)) < before) EXIT
               fraction = 0.5_dp*fraction
            END DO
            IF (k > 24) RETURN
            coordinate(:) = candidate
         END DO
      END IF
      valid = ALL(ABS(values) <= chart%tolerances)
   END SUBROUTINE skala_constraint_chart_restore

! **************************************************************************************************
!> \brief Check attainable limiting derivatives along physical constraint normals.
!>        Gradients and constraint Jacobians use the same orthonormal real orbital tangent chart.
!>        Proportional normals share one variation. Retain every sample and switch their limits
!>        together, reversing the side of oppositely oriented samples. Dual normals retain the
!>        other independent groups. Their physical density units do not affect the signs.
!>        This is a normal admissibility check, not proof of tangential stability or stationarity.
!> \param positive_gradient Complete all-positive limiting orbital derivative.
!> \param negative_changes Independent negative-minus-positive changes, one column per constraint.
!> \param jacobian Physical constraint derivatives, one row per constraint.
!> \param admissible Whether a convex cusp brackets zero along each positive physical normal.
!> \param release_side Minus one for a negative-side descent, plus one for a positive-side descent.
!> \param valid Whether the nonzero normals form independent groups of proportional samples.
! **************************************************************************************************
   SUBROUTINE skala_cusp_admissibility(positive_gradient, negative_changes, jacobian, &
                                       admissible, release_side, valid)
      REAL(dp), INTENT(IN)                               :: positive_gradient(:), &
                                                            negative_changes(:, :), jacobian(:, :)
      LOGICAL, INTENT(OUT)                               :: admissible(:)
      INTEGER, INTENT(OUT)                               :: release_side(:)
      LOGICAL, INTENT(OUT)                               :: valid

      INTEGER                                            :: d, group, i, info, j, lwork, m, ngroups, &
                                                            rank, side
      INTEGER, ALLOCATABLE                               :: groups(:), orientation(:)
      REAL(dp)                                           :: negative_slope, positive_slope, &
                                                            query(1), roundoff, scale
      REAL(dp), ALLOCATABLE                              :: basis(:, :), matrix(:, :), &
                                                            normals(:, :), row(:), singular(:), &
                                                            work(:)

      d = SIZE(positive_gradient)
      m = SIZE(jacobian, 1)
      CPASSERT(SIZE(jacobian, 2) == d)
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
      IF (d == 0) RETURN
      ALLOCATE (basis(m, d), groups(m), orientation(m), row(d))
      ngroups = 0
      roundoff = 64.0_dp*EPSILON(0.0_dp)
      DO i = 1, m
         scale = MAXVAL(ABS(jacobian(i, :)))
         IF (scale <= TINY(0.0_dp)) RETURN
         row(:) = jacobian(i, :)/scale
         row(:) = row/NORM2(row)
         group = 0
         orientation(i) = 1
         DO j = 1, ngroups
            IF (NORM2(row - basis(j, :)) <= roundoff) THEN
               group = j
               EXIT
            ELSE IF (NORM2(row + basis(j, :)) <= roundoff) THEN
               group = j
               orientation(i) = -1
               EXIT
            END IF
         END DO
         IF (group == 0) THEN
            ngroups = ngroups + 1
            basis(ngroups, :) = row
            group = ngroups
         END IF
         groups(i) = group
      END DO
      IF (ngroups > d) RETURN
      ALLOCATE (matrix(ngroups, d), normals(d, ngroups), singular(ngroups))
      matrix(:, :) = basis(:ngroups, :)
      normals = 0.0_dp
      DO i = 1, ngroups
         normals(i, i) = 1.0_dp
      END DO
      CALL dgelss(ngroups, d, ngroups, matrix, ngroups, normals, d, singular, EPSILON(0.0_dp), rank, query, -1, info)
      IF (info /= 0) RETURN
      lwork = MAX(1, CEILING(query(1)))
      ALLOCATE (work(lwork))
      CALL dgelss(ngroups, d, ngroups, matrix, ngroups, normals, d, singular, EPSILON(0.0_dp), rank, work, lwork, info)
      IF (info /= 0 .OR. rank /= ngroups) RETURN
      DO group = 1, ngroups
         normals(:, group) = normals(:, group)/NORM2(normals(:, group))
         positive_slope = DOT_PRODUCT(positive_gradient, normals(:, group))
         negative_slope = positive_slope
         DO i = 1, m
            IF (groups(i) /= group) CYCLE
            IF (orientation(i) > 0) THEN
               negative_slope = negative_slope + DOT_PRODUCT(negative_changes(:, i), normals(:, group))
            ELSE
               positive_slope = positive_slope + DOT_PRODUCT(negative_changes(:, i), normals(:, group))
            END IF
         END DO
         WHERE (groups == group) admissible = negative_slope <= 0.0_dp .AND. positive_slope >= 0.0_dp
         IF (negative_slope <= 0.0_dp .AND. positive_slope >= 0.0_dp) CYCLE
         side = 0
         IF (-positive_slope >= negative_slope .AND. positive_slope < 0.0_dp) THEN
            side = 1
         ELSE IF (negative_slope > 0.0_dp) THEN
            side = -1
         END IF
         WHERE (groups == group) release_side = side*orientation
      END DO
      valid = .TRUE.
   END SUBROUTINE skala_cusp_admissibility

! **************************************************************************************************
!> \brief Select independent rows without discarding any physical residual checks.
!>        Column-pivoted QR of the transposed, row-scaled Jacobian determines its numerical rank.
!> \param jacobian All original physical constraint derivatives in the same orbital coordinates.
!> \param scales Positive physical row scales.
!> \param rows Permutation with the independent physical samples first.
!> \param rank Number of independent rows at the current regular point.
!> \param success Whether the rank-revealing factorization completed.
! **************************************************************************************************
   SUBROUTINE independent_constraint_rows(jacobian, scales, rows, rank, success)
      REAL(dp), INTENT(IN)                               :: jacobian(:, :), scales(:)
      INTEGER, INTENT(OUT)                               :: rows(:), rank
      LOGICAL, INTENT(OUT)                               :: success

      INTEGER                                            :: i, info, lwork, m, n
      REAL(dp)                                           :: query(1), threshold
      REAL(dp), ALLOCATABLE                              :: matrix(:, :), tau(:), work(:)

      m = SIZE(jacobian, 1)
      n = SIZE(jacobian, 2)
      CPASSERT(m > 0 .AND. n > 0 .AND. SIZE(rows) == m .AND. SIZE(scales) == m)
      CPASSERT(ALL(scales > 0.0_dp))
      ALLOCATE (matrix(n, m), tau(MIN(m, n)))
      matrix(:, :) = TRANSPOSE(jacobian/SPREAD(scales, 2, n))
      rows = 0
      rank = 0
      CALL dgeqp3(n, m, matrix, n, rows, tau, query, -1, info)
      success = info == 0
      IF (.NOT. success) RETURN
      lwork = MAX(1, CEILING(query(1)))
      ALLOCATE (work(lwork))
      CALL dgeqp3(n, m, matrix, n, rows, tau, work, lwork, info)
      success = info == 0
      IF (.NOT. success) RETURN
      threshold = 64.0_dp*EPSILON(0.0_dp)*ABS(matrix(1, 1))
      DO i = 1, MIN(m, n)
         IF (ABS(matrix(i, i)) <= threshold) EXIT
         rank = rank + 1
      END DO
   END SUBROUTINE independent_constraint_rows

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
      REAL(dp) :: beta, scale, inverse_curvature(SIZE(steps, 2)), alpha(SIZE(steps, 2))

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
