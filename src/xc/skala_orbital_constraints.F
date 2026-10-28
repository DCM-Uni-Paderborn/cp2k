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
             skala_orbital_chart_size, skala_orbital_chart_evaluate, skala_orbital_chart_derivative, &
             skala_orbital_chart_set_occupations, skala_orbital_chart_precondition, skala_orbital_chart_scales, &
             skala_orbital_chart_response_inverse, skala_orbital_chart_density_response, &
             skala_orbital_chart_operator_response

   PUBLIC :: skala_orbital_chart_set_reference
   PUBLIC :: skala_orbital_chart_restrict
   PUBLIC :: skala_orbital_chart_symmetry, skala_orbital_multiplets

   TYPE skala_orbital_chart_type
      PRIVATE
      COMPLEX(dp), ALLOCATABLE :: base(:, :), coefficients(:, :), inverse(:, :), rotation(:, :), &
                                  work(:, :), transformed(:, :)
      REAL(dp), ALLOCATABLE :: occupations(:), coordinate(:), coordinate_scale(:)
      INTEGER, ALLOCATABLE :: pairs(:, :), component(:), pivots(:)
      INTEGER, ALLOCATABLE :: occupation_groups(:)
      COMPLEX(dp), ALLOCATABLE :: symmetry(:, :, :)
      LOGICAL, ALLOCATABLE :: antiunitary(:)
      REAL(dp) :: metric_weight = 1.0_dp
      LOGICAL :: evaluated = .FALSE., variable_occupations = .FALSE.
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
! **************************************************************************************************
!> \brief Check a restored curvature probe for additional physical thresholds.
!> \param first Accepted feasible coordinates.
!> \param last Restored probe coordinates.
!> \param curve Current independent/dependent coordinate partition.
!> \param context Borrowed physical state and boundary admission settings.
!> \param smooth Whether no additional physical threshold is crossed.
!> \param interrupted Whether boundary or external-stop handling is required.
! **************************************************************************************************
      SUBROUTINE orbital_curvature_segment(first, last, curve, context, smooth, interrupted)
         IMPORT :: dp, skala_constraint_chart_type
         REAL(dp), INTENT(IN) :: first(:), last(:)
         TYPE(skala_constraint_chart_type), INTENT(IN) :: curve
         CLASS(*), INTENT(INOUT) :: context
         LOGICAL, INTENT(OUT) :: smooth, interrupted
      END SUBROUTINE orbital_curvature_segment
! **************************************************************************************************
!> \brief Apply a supplied full-coordinate energy Hessian at an unconstrained state.
!> \param coordinate Complete accepted orbital coordinates.
!> \param direction Full-coordinate direction.
!> \param context Borrowed problem state.
!> \param PRODUCT Full-coordinate Hessian product.
!> \param available Whether the complete physical response is supplied.
! **************************************************************************************************
      SUBROUTINE orbital_hessian_product(coordinate, direction, context, PRODUCT, available)
         IMPORT :: dp
         REAL(dp), INTENT(IN) :: coordinate(:), direction(:)
         CLASS(*), INTENT(INOUT) :: context
         REAL(dp), INTENT(OUT) :: PRODUCT(:)
         LOGICAL, INTENT(OUT) :: available
      END SUBROUTINE orbital_hessian_product
! **************************************************************************************************
!> \brief Apply the Lagrangian Hessian to a feasible full-coordinate direction.
!> \param coordinate Accepted coordinates on the original physical constraint curve.
!> \param direction Full tangent direction, including dependent-coordinate response.
!> \param multipliers Constraint coefficients in their original physical units.
!> \param context Borrowed problem state.
!> \param PRODUCT Hessian of E minus the multiplier-weighted constraints, applied to direction.
!> \param available Whether this complete constrained response is supplied.
! **************************************************************************************************
      SUBROUTINE orbital_lagrangian_hessian(coordinate, direction, multipliers, context, PRODUCT, available)
         IMPORT :: dp
         REAL(dp), INTENT(IN) :: coordinate(:), direction(:), multipliers(:)
         CLASS(*), INTENT(INOUT) :: context
         REAL(dp), INTENT(OUT) :: PRODUCT(:)
         LOGICAL, INTENT(OUT) :: available
      END SUBROUTINE orbital_lagrangian_hessian
! **************************************************************************************************
!> \brief Supply a positive inverse diagonal for the coupled orbital step.
!> \param coordinate Complete accepted coordinates.
!> \param context Borrowed problem state.
!> \param inverse_diagonal Current positive diagonal inverse.
!> \param available Whether this preconditioner is supplied.
! **************************************************************************************************
      SUBROUTINE orbital_preconditioner(coordinate, context, inverse_diagonal, available)
         IMPORT :: dp
         REAL(dp), INTENT(IN) :: coordinate(:)
         CLASS(*), INTENT(INOUT) :: context
         REAL(dp), INTENT(OUT) :: inverse_diagonal(:)
         LOGICAL, INTENT(OUT) :: available
      END SUBROUTINE orbital_preconditioner
! **************************************************************************************************
!> \brief Report stored accepted-state scalars without another physical evaluation.
!> \param iteration Attempted iteration in the current chart, zero at its reference.
!> \param energy Energy of the last accepted feasible coordinates.
!> \param gradient_norm Norm of the unscaled independent-coordinate covector.
!> \param moved Whether this iteration accepted a new step.
!> \param context Borrowed caller output context.
! **************************************************************************************************
      SUBROUTINE orbital_iteration_report(iteration, energy, gradient_norm, moved, context)
         IMPORT :: dp
         INTEGER, INTENT(IN) :: iteration
         REAL(dp), INTENT(IN) :: energy, gradient_norm
         LOGICAL, INTENT(IN) :: moved
         CLASS(*), INTENT(INOUT) :: context
      END SUBROUTINE orbital_iteration_report
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
!> \param variable_occupations ...
! **************************************************************************************************
   SUBROUTINE skala_orbital_chart_create(coefficients, metric, occupations, complex_orbitals, chart, valid, &
                                         metric_weight, complete_space, variable_occupations)
      COMPLEX(dp), INTENT(IN)                            :: coefficients(:, :), metric(:, :)
      REAL(dp), INTENT(IN)                               :: occupations(:)
      LOGICAL, INTENT(IN)                                :: complex_orbitals
      TYPE(skala_orbital_chart_type), INTENT(INOUT)      :: chart
      LOGICAL, INTENT(OUT)                               :: valid
      REAL(dp), INTENT(IN), OPTIONAL                     :: metric_weight
      LOGICAL, INTENT(IN), OPTIONAL                      :: complete_space, variable_occupations

      COMPLEX(dp), ALLOCATABLE                           :: gram(:, :)
      INTEGER                                            :: a, b, count, i, info, n, nreference, rows
      LOGICAL                                            :: varying
      REAL(dp), ALLOCATABLE                              :: full_occupations(:)

      CALL skala_orbital_chart_release(chart)
      rows = SIZE(coefficients, 1)
      nreference = SIZE(coefficients, 2)
      n = nreference
      varying = .FALSE.
      IF (PRESENT(variable_occupations)) varying = variable_occupations
      chart%variable_occupations = varying
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
            IF (varying) THEN
               IF (a > nreference .AND. b > nreference) CYCLE
            ELSE
               IF (full_occupations(a) == full_occupations(b)) CYCLE
            END IF
            count = count + MERGE(2, 1, complex_orbitals)
         END DO
      END DO
      ALLOCATE (chart%base(rows, n), chart%coefficients(rows, n), chart%inverse(n, n), &
                chart%rotation(n, n), chart%work(n, n), chart%transformed(n, n), &
                chart%occupations(n), chart%coordinate(count), chart%coordinate_scale(count), chart%pairs(2, count), &
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
      chart%coordinate_scale = 1.0_dp
      i = 0
      DO b = 2, n
         DO a = 1, b - 1
            IF (varying) THEN
               IF (a > nreference .AND. b > nreference) CYCLE
            ELSE
               IF (full_occupations(a) == full_occupations(b)) CYCLE
            END IF
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
!> \brief Change an initial orbital reference while retaining its occupation and rotation layout.
!>        The supplied reference spans the same complete space. The caller checks density and
!>        symmetry preservation. Apply only at zero coordinates before starting a new search.
!> \param chart Initial orbital chart, possibly already evaluated at zero coordinates.
!> \param coefficients New complete reference in the same AO order.
!> \param metric Original Hermitian charge metric.
!> \param valid Whether the new reference has a positive metric Gram matrix.
! **************************************************************************************************
   SUBROUTINE skala_orbital_chart_set_reference(chart, coefficients, metric, valid)
      TYPE(skala_orbital_chart_type), INTENT(INOUT)      :: chart
      COMPLEX(dp), INTENT(IN)                            :: coefficients(:, :), metric(:, :)
      LOGICAL, INTENT(OUT)                               :: valid

      COMPLEX(dp), ALLOCATABLE                           :: gram(:, :)
      INTEGER                                            :: info, n, rows

      CPASSERT(ALLOCATED(chart%base))
      CPASSERT(.NOT. ALLOCATED(chart%symmetry))
      IF (chart%evaluated) THEN
         CPASSERT(ALL(chart%coordinate == 0.0_dp))
      END IF
      CPASSERT(ALL(SHAPE(coefficients) == SHAPE(chart%base)))
      rows = SIZE(coefficients, 1)
      n = SIZE(coefficients, 2)
      CPASSERT(ALL(SHAPE(metric) == [rows, rows]))
      gram = MATMUL(CONJG(TRANSPOSE(coefficients)), MATMUL(metric, coefficients))
      CALL zpotrf('L', n, gram, n, info)
      valid = info == 0
      IF (.NOT. valid) RETURN
      chart%base = coefficients
      CALL ztrsm('R', 'L', 'C', 'N', rows, n, CMPLX(1.0_dp, 0.0_dp, dp), gram, n, chart%base, rows)
      chart%coordinate = 0.0_dp
      chart%evaluated = .FALSE.
   END SUBROUTINE skala_orbital_chart_set_reference

! **************************************************************************************************
!> \brief Restrict orbital rotations to supplied symmetry subspaces without changing occupations.
!>        The caller supplies a complete reference already adapted to these subspaces and checks
!>        its covariance in the original charge metric. Apply only before the first trial evaluation.
!>        Real generators additionally preserve conjugation in a compatible antiunitary frame.
!>        Density averaging, occupation rounding and a change of the reference are not performed.
!> \param chart Unevaluated metric-normalized orbital chart.
!> \param rotation_groups Positive subspace label for every reference column.
!> \param real_generators Keep only real rotations in the supplied reference frame.
! **************************************************************************************************
   SUBROUTINE skala_orbital_chart_restrict(chart, rotation_groups, real_generators)
      TYPE(skala_orbital_chart_type), INTENT(INOUT)      :: chart
      INTEGER, INTENT(IN)                                :: rotation_groups(:)
      LOGICAL, INTENT(IN)                                :: real_generators

      INTEGER                                            :: i, nkeep
      LOGICAL, ALLOCATABLE                               :: keep(:)

      CPASSERT(ALLOCATED(chart%base) .AND. ALLOCATED(chart%coordinate))
      CPASSERT(.NOT. ALLOCATED(chart%symmetry))
      CPASSERT(.NOT. chart%evaluated)
      CPASSERT(SIZE(rotation_groups) == SIZE(chart%base, 2))
      CPASSERT(ALL(rotation_groups > 0))
      ALLOCATE (keep(SIZE(chart%coordinate)))
      DO i = 1, SIZE(keep)
         keep(i) = rotation_groups(chart%pairs(1, i)) == rotation_groups(chart%pairs(2, i))
         IF (real_generators) keep(i) = keep(i) .AND. chart%component(i) == 1
      END DO
      nkeep = COUNT(keep)
      chart%pairs = RESHAPE(PACK(chart%pairs, SPREAD(keep, 1, 2)), [2, nkeep])
      chart%component = PACK(chart%component, keep)
      DEALLOCATE (chart%coordinate)
      ALLOCATE (chart%coordinate(nkeep), SOURCE=0.0_dp)
      chart%coordinate_scale = PACK(chart%coordinate_scale, keep)
   END SUBROUTINE skala_orbital_chart_restrict

! **************************************************************************************************
!> \brief Preserve a complete little group in the original metric-normalized reference.
!>        The real-linear group average acts on Cayley generators and its exact transpose acts
!>        on covectors. Coordinate storage stays quadratic rather than storing a dense matrix
!>        for every allowed generator. Projected coordinates may include null directions.
!>        The caller supplies the complete group, including its identity. Scalar projective phases
!>        cancel in the conjugation action. Initial occupations must commute with every action.
!>        Thermal occupations may vary independently between the supplied multiplets.
!> \param chart Initial orbital chart at zero coordinates.
!> \param actions Complete unitary parts of the unitary or antiunitary reference-space actions.
!> \param antiunitary Whether the corresponding action conjugates orbital coefficients.
!> \param valid Whether the original reference and occupation variables admit these actions.
!> \param occupation_groups Optional multiplet labels. Occupations are equal within each label.
! **************************************************************************************************
   SUBROUTINE skala_orbital_chart_symmetry(chart, actions, antiunitary, valid, occupation_groups)
      TYPE(skala_orbital_chart_type), INTENT(INOUT)      :: chart
      COMPLEX(dp), INTENT(IN)                            :: actions(:, :, :)
      LOGICAL, INTENT(IN)                                :: antiunitary(:)
      LOGICAL, INTENT(OUT)                               :: valid
      INTEGER, INTENT(IN), OPTIONAL                      :: occupation_groups(:)

      REAL(dp), PARAMETER                                :: tolerance = 3.e-11_dp

      COMPLEX(dp), ALLOCATABLE                           :: gram(:, :)
      INTEGER                                            :: i, j, n, op
      INTEGER, ALLOCATABLE                               :: groups(:)

      CPASSERT(ALLOCATED(chart%base))
      IF (chart%evaluated) THEN
         CPASSERT(ALL(chart%coordinate == 0.0_dp))
      END IF
      n = SIZE(chart%base, 2)
      CPASSERT(SIZE(actions, 1) == n .AND. SIZE(actions, 2) == n)
      CPASSERT(SIZE(actions, 3) == SIZE(antiunitary))
      valid = .FALSE.
      IF (SIZE(antiunitary) == 0) RETURN
      IF (.NOT. ALL(ieee_is_finite(REAL(actions, dp)))) RETURN
      IF (.NOT. ALL(ieee_is_finite(AIMAG(actions)))) RETURN
      groups = [(i, i=1, n)]
      IF (PRESENT(occupation_groups)) THEN
         CPASSERT(SIZE(occupation_groups) == n)
         DO i = 1, n
            groups(i) = MINVAL(PACK([(j, j=1, n)], occupation_groups == occupation_groups(i)))
         END DO
      END IF
      IF (chart%variable_occupations) THEN
         IF (MAXVAL(ABS(chart%occupations - chart%occupations(groups))) > tolerance) RETURN
      END IF
      ALLOCATE (gram(n, n))
      DO op = 1, SIZE(antiunitary)
         gram = MATMUL(CONJG(TRANSPOSE(actions(:, :, op))), actions(:, :, op))
         DO i = 1, n
            gram(i, i) = gram(i, i) - 1.0_dp
         END DO
         IF (MAXVAL(ABS(gram)) > tolerance) RETURN
         DO j = 1, n
            DO i = 1, n
               IF (i == j) CYCLE
               IF (ABS((chart%occupations(i) - chart%occupations(j))*actions(i, j, op)) > tolerance) RETURN
               IF (chart%variable_occupations) THEN
                  IF (groups(i) /= groups(j) .AND. ABS(actions(i, j, op)) > tolerance) RETURN
               END IF
            END DO
         END DO
      END DO
      chart%symmetry = actions
      chart%antiunitary = antiunitary
      IF (chart%variable_occupations) chart%occupation_groups = groups
      chart%evaluated = .FALSE.
      valid = .TRUE.
   END SUBROUTINE skala_orbital_chart_symmetry

! **************************************************************************************************
!> \brief Resolve irreducible thermal multiplets of a complete unitary and conjugating group.
!>        Preserve retained-band, occupation and auxiliary-energy sectors to symmetry tolerance.
!>        A nonscalar averaged Hermitian basis element splits equivalent copies. Repeating until
!>        every averaged basis element is scalar avoids tying independent copies accidentally.
!>        Permanently empty completion sectors need no irreducible decomposition or gauge change.
!> \param actions Reference-space actions, rotated into the resulting multiplet basis.
!> \param antiunitary Conjugating actions.
!> \param occupations Original complete occupation vector.
!> \param reference_energies Original complete auxiliary energies, in hartree.
!> \param nretained Number of thermally variable bands preceding the empty completion.
!> \param gauge Unitary reference rotation.
!> \param groups Minimal reference index of each irreducible multiplet.
!> \param energies New auxiliary energies, scalar inside each multiplet.
!> \param valid Whether the original reference admits the construction.
! **************************************************************************************************
   SUBROUTINE skala_orbital_multiplets(actions, antiunitary, occupations, reference_energies, nretained, &
                                       gauge, groups, energies, valid)
      COMPLEX(dp), INTENT(INOUT)                         :: actions(:, :, :)
      LOGICAL, INTENT(IN)                                :: antiunitary(:)
      REAL(dp), INTENT(IN)                               :: occupations(:), reference_energies(:)
      INTEGER, INTENT(IN)                                :: nretained
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: gauge(:, :)
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: groups(:)
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: energies(:)
      LOGICAL, INTENT(OUT)                               :: valid

      REAL(dp), PARAMETER                                :: tolerance = 3.e-11_dp

      COMPLEX(dp), ALLOCATABLE                           :: action(:, :), basis(:, :), label(:, :), &
                                                            source(:, :), term(:, :), work(:)
      INTEGER                                            :: first, i, j, m, mode, n, op, status
      INTEGER, ALLOCATABLE                               :: indices(:)
      LOGICAL                                            :: consistent, split
      REAL(dp)                                           :: scale
      REAL(dp), ALLOCATABLE                              :: real_work(:), values(:)

      valid = .FALSE.
      n = SIZE(occupations)
      consistent = SIZE(reference_energies) == n .AND. nretained >= 0 .AND. nretained <= n
      CPASSERT(consistent)
      consistent = ALL(SHAPE(actions) == [n, n, SIZE(antiunitary)])
      CPASSERT(consistent)
      IF (SIZE(antiunitary) == 0) RETURN
      IF (.NOT. ALL(ieee_is_finite(occupations)) .OR. .NOT. ALL(ieee_is_finite(reference_energies))) RETURN
      IF (.NOT. ALL(ieee_is_finite(REAL(actions, dp))) .OR. .NOT. ALL(ieee_is_finite(AIMAG(actions)))) RETURN
      ALLOCATE (groups(n), gauge(n, n))
      groups = 0
      gauge = 0.0_dp
      DO i = 1, n
         gauge(i, i) = 1.0_dp
         IF (groups(i) /= 0) CYCLE
         groups(i) = i
         DO j = i + 1, n
            IF (groups(j) /= 0) CYCLE
            IF ((i <= nretained) .NEQV. (j <= nretained)) CYCLE
            IF (ABS(occupations(i) - occupations(j)) > tolerance) CYCLE
            scale = MAX(1.0_dp, ABS(reference_energies(i)), ABS(reference_energies(j)))
            IF (ABS(reference_energies(i) - reference_energies(j)) > tolerance*scale) CYCLE
            groups(j) = i
         END DO
      END DO
      DO op = 1, SIZE(antiunitary)
         source = MATMUL(CONJG(TRANSPOSE(actions(:, :, op))), actions(:, :, op))
         DO i = 1, n
            source(i, i) = source(i, i) - 1.0_dp
         END DO
         IF (MAXVAL(ABS(source)) > tolerance) RETURN
         DO j = 1, n
            DO i = 1, n
               IF (groups(i) /= groups(j) .AND. ABS(actions(i, j, op)) > tolerance) RETURN
            END DO
         END DO
      END DO
      DO
         split = .FALSE.
         DO first = 1, n
            IF (groups(first) /= first) CYCLE
            IF (first > nretained) CYCLE
            indices = PACK([(i, i=1, n)], groups == first)
            m = SIZE(indices)
            IF (m <= 1) CYCLE
            ALLOCATE (basis(m, m), label(m, m), values(m), work(2*m), real_work(3*m))
            search: DO mode = 1, 3
               DO j = 1, m
                  DO i = 1, j
                     IF ((mode == 1) .NEQV. (i == j)) CYCLE
                     basis = 0.0_dp
                     SELECT CASE (mode)
                     CASE (1)
                        basis(i, i) = 1.0_dp
                     CASE (2)
                        basis(i, j) = 1.0_dp
                        basis(j, i) = 1.0_dp
                     CASE (3)
                        basis(i, j) = CMPLX(0.0_dp, 1.0_dp, dp)
                        basis(j, i) = CMPLX(0.0_dp, -1.0_dp, dp)
                     END SELECT
                     label = 0.0_dp
                     DO op = 1, SIZE(antiunitary)
                        action = actions(indices, indices, op)
                        term = basis
                        IF (antiunitary(op)) term = CONJG(term)
                        label = label + MATMUL(action, MATMUL(term, CONJG(TRANSPOSE(action))))/REAL(SIZE(antiunitary), dp)
                     END DO
                     label = 0.5_dp*(label + CONJG(TRANSPOSE(label)))
                     scale = SUM([(REAL(label(op, op), dp), op=1, m)])/REAL(m, dp)
                     DO op = 1, m
                        label(op, op) = label(op, op) - scale
                     END DO
                     IF (MAXVAL(ABS(label)) <= tolerance) CYCLE
                     CALL zheev('V', 'U', m, label, m, values, work, SIZE(work), real_work, status)
                     IF (status /= 0) RETURN
                     IF (MAXVAL(values(2:) - values(:m - 1)) <= tolerance) RETURN
                     gauge(:, indices) = MATMUL(gauge(:, indices), label)
                     DO op = 1, SIZE(antiunitary)
                        basis = label
                        IF (antiunitary(op)) basis = CONJG(basis)
                        source = actions(:, :, op)
                        source(:, indices) = MATMUL(source(:, indices), basis)
                        source(indices, :) = MATMUL(CONJG(TRANSPOSE(label)), source(indices, :))
                        actions(:, :, op) = source
                     END DO
                     groups(indices(1)) = indices(1)
                     DO op = 2, m
                        groups(indices(op)) = groups(indices(op - 1))
                        IF (values(op) - values(op - 1) > tolerance) groups(indices(op)) = indices(op)
                     END DO
                     split = .TRUE.
                     EXIT search
                  END DO
               END DO
            END DO search
            DEALLOCATE (basis, label, values, work, real_work)
            IF (split) EXIT
         END DO
         IF (.NOT. split) EXIT
      END DO
      energies = SUM(ABS(gauge)**2*SPREAD(reference_energies, 2, n), DIM=1)
      DO first = 1, n
         IF (groups(first) /= first) CYCLE
         indices = PACK([(i, i=1, n)], groups == first)
         energies(indices) = SUM(energies(indices))/REAL(SIZE(indices), dp)
      END DO
      source = MATMUL(gauge*SPREAD(energies, 1, n), CONJG(TRANSPOSE(gauge)))
      DO i = 1, n
         source(i, i) = source(i, i) - reference_energies(i)
      END DO
      IF (MAXVAL(ABS(source)) > tolerance*MAX(1.0_dp, MAXVAL(ABS(reference_energies)))) RETURN
      valid = .TRUE.
   END SUBROUTINE skala_orbital_multiplets

! **************************************************************************************************
!> \brief Apply the real-linear group average or its exact transpose in the reference frame.
!> \param actions Optional complete reference-space symmetry group.
!> \param antiunitary Whether the corresponding action conjugates coefficients.
!> \param matrix Reference-space matrix, replaced by its group average.
!> \param adjoint Use the transpose for the real Frobenius pairing.
! **************************************************************************************************
   SUBROUTINE orbital_symmetry_project(actions, antiunitary, matrix, adjoint)
      COMPLEX(dp), ALLOCATABLE, INTENT(IN)               :: actions(:, :, :)
      LOGICAL, ALLOCATABLE, INTENT(IN)                   :: antiunitary(:)
      COMPLEX(dp), INTENT(INOUT)                         :: matrix(:, :)
      LOGICAL, INTENT(IN)                                :: adjoint

      COMPLEX(dp), ALLOCATABLE                           :: projected(:, :), source(:, :), term(:, :)
      INTEGER                                            :: op

      IF (.NOT. ALLOCATED(actions)) RETURN
      ALLOCATE (projected, MOLD=matrix)
      projected = 0.0_dp
      DO op = 1, SIZE(antiunitary)
         IF (adjoint) THEN
            term = MATMUL(CONJG(TRANSPOSE(actions(:, :, op))), MATMUL(matrix, actions(:, :, op)))
            IF (antiunitary(op)) term = CONJG(term)
         ELSE
            source = matrix
            IF (antiunitary(op)) source = CONJG(source)
            term = MATMUL(actions(:, :, op), MATMUL(source, CONJG(TRANSPOSE(actions(:, :, op)))))
         END IF
         projected = projected + term/REAL(SIZE(antiunitary), dp)
      END DO
      matrix = projected
   END SUBROUTINE orbital_symmetry_project

! **************************************************************************************************
!> \brief Update physical occupations without changing the retained orbital coordinates.
!>        Variable-occupation charts retain equal-occupation rotations because later thermal
!>        updates can make them physical. Only rotations between permanently omitted columns
!>        are excluded. The original charge metric and Cayley factors stay unchanged.
!> \param chart Evaluated or unevaluated variable-occupation chart.
!> \param occupations Complete occupation vector, including zero occupations in its complement.
! **************************************************************************************************
   SUBROUTINE skala_orbital_chart_set_occupations(chart, occupations)
      TYPE(skala_orbital_chart_type), INTENT(INOUT)      :: chart
      REAL(dp), INTENT(IN)                               :: occupations(:)

      LOGICAL                                            :: invariant

      CPASSERT(chart%variable_occupations)
      CPASSERT(SIZE(occupations) == SIZE(chart%occupations))
      CPASSERT(ALL(ieee_is_finite(occupations)) .AND. ALL(occupations >= 0.0_dp))
      IF (ALLOCATED(chart%occupation_groups)) THEN
         invariant = MAXVAL(ABS(occupations - occupations(chart%occupation_groups))) <= 3.e-11_dp
         CPASSERT(invariant)
      END IF
      chart%occupations(:) = occupations
   END SUBROUTINE skala_orbital_chart_set_occupations

! **************************************************************************************************
!> \brief Freeze a positive diagonal approximation to the one-electron rotation Hessian.
!>        Energies are Rayleigh values in the normalized reference. The occupation metric and
!>        temperature regularize equal-occupation rotations. The map stays fixed during a search.
!> \param chart Variable-occupation chart at its reference coordinate.
!> \param energies Complete reference Rayleigh energies, in hartree.
!> \param temperature Positive kBT, in hartree.
! **************************************************************************************************
   SUBROUTINE skala_orbital_chart_precondition(chart, energies, temperature)
      TYPE(skala_orbital_chart_type), INTENT(INOUT)      :: chart
      REAL(dp), INTENT(IN)                               :: energies(:), temperature

      INTEGER                                            :: a, b, i
      REAL(dp)                                           :: curvature, floor, metric, peak

      CPASSERT(chart%variable_occupations)
      CPASSERT(SIZE(energies) == SIZE(chart%occupations))
      CPASSERT(ALL(ieee_is_finite(energies)))
      CPASSERT(ieee_is_finite(temperature) .AND. temperature > 0.0_dp)
      IF (chart%evaluated) THEN
         CPASSERT(ALL(chart%coordinate == 0.0_dp))
      END IF
      IF (SIZE(chart%coordinate_scale) == 0) RETURN
      peak = MAXVAL(chart%occupations)
      DO i = 1, SIZE(chart%coordinate_scale)
         a = chart%pairs(1, i)
         b = chart%pairs(2, i)
         ! Nearly empty reference pairs can acquire occupation during the same search.
         metric = chart%metric_weight*MAX(peak, chart%occupations(a) + chart%occupations(b))
         curvature = 2.0_dp*chart%metric_weight*ABS(chart%occupations(a) - chart%occupations(b))* &
                     MAX(ABS(energies(a) - energies(b)), temperature)
         chart%coordinate_scale(i) = MAX(curvature, temperature*metric)
      END DO
      floor = 64.0_dp*EPSILON(1.0_dp)*MAX(MAXVAL(chart%coordinate_scale), temperature*chart%metric_weight)
      WHERE (chart%coordinate_scale > 0.0_dp)
         chart%coordinate_scale = 1.0_dp/SQRT(MAX(chart%coordinate_scale, floor))
      ELSE WHERE
         chart%coordinate_scale = 1.0_dp
      END WHERE
      chart%evaluated = .FALSE.
   END SUBROUTINE skala_orbital_chart_precondition

! **************************************************************************************************
!> \brief Positive inverse of the current occupation-weighted one-electron rotation model.
!>        Exact pair coordinates coincide with physical rotations at the reference point.
!>        The original fixed coordinate scales remain part of the model curvature.
!> \param chart Evaluated variable-occupation chart.
!> \param energies Current orbital Rayleigh energies in hartree.
!> \param temperature Positive kBT in hartree.
!> \param inverse_diagonal Positive inverse curvature for each independent coordinate.
! **************************************************************************************************
   SUBROUTINE skala_orbital_chart_response_inverse(chart, energies, temperature, inverse_diagonal)
      TYPE(skala_orbital_chart_type), INTENT(IN)         :: chart
      REAL(dp), INTENT(IN)                               :: energies(:), temperature
      REAL(dp), INTENT(OUT)                              :: inverse_diagonal(:)

      INTEGER                                            :: a, b, i
      REAL(dp)                                           :: curvature, log_inverse, metric

      CPASSERT(chart%evaluated)
      CPASSERT(SIZE(energies) == SIZE(chart%occupations))
      CPASSERT(SIZE(inverse_diagonal) == SIZE(chart%coordinate_scale))
      CPASSERT(ALL(ieee_is_finite(energies)))
      CPASSERT(temperature > 0.0_dp)
      inverse_diagonal = 1.0_dp
      DO i = 1, SIZE(inverse_diagonal)
         a = chart%pairs(1, i)
         b = chart%pairs(2, i)
         metric = chart%metric_weight*(chart%occupations(a) + chart%occupations(b))
         curvature = 2.0_dp*chart%metric_weight*ABS(chart%occupations(a) - chart%occupations(b))* &
                     MAX(ABS(energies(a) - energies(b)), temperature)
         curvature = MAX(curvature, temperature*metric)
         IF (curvature <= 0.0_dp) CYCLE
         log_inverse = -LOG(curvature) - 2.0_dp*LOG(chart%coordinate_scale(i))
         inverse_diagonal(i) = EXP(MAX(LOG(TINY(1.0_dp)), &
                                       MIN(LOG(HUGE(1.0_dp)) - LOG(2.0_dp), log_inverse)))
      END DO
   END SUBROUTINE skala_orbital_chart_response_inverse

! **************************************************************************************************
!> \brief Read the fixed diagonal map from chart coordinates to Cayley generators.
!> \param chart Original chart and its unchanged reference map.
!> \param scales Positive scale of every independent real coordinate.
! **************************************************************************************************
   SUBROUTINE skala_orbital_chart_scales(chart, scales)
      TYPE(skala_orbital_chart_type), INTENT(IN)         :: chart
      REAL(dp), INTENT(OUT)                              :: scales(:)

      CPASSERT(SIZE(scales) == SIZE(chart%coordinate_scale))
      scales = chart%coordinate_scale
   END SUBROUTINE skala_orbital_chart_scales

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
                     chart%work, chart%transformed, chart%occupations, chart%coordinate, chart%coordinate_scale, &
                     chart%pairs, chart%component, chart%pivots)
      END IF
      IF (ALLOCATED(chart%symmetry)) DEALLOCATE (chart%symmetry, chart%antiunitary)
      IF (ALLOCATED(chart%occupation_groups)) DEALLOCATE (chart%occupation_groups)
      chart%evaluated = .FALSE.
      chart%metric_weight = 1.0_dp
      chart%variable_occupations = .FALSE.
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

      COMPLEX(dp), ALLOCATABLE                           :: generator(:, :)
      INTEGER                                            :: i, info, n

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
      CALL orbital_generator_direction(chart, coordinate, generator)
      chart%work = -0.5_dp*generator
      chart%inverse = CMPLX(0.0_dp, 0.0_dp, dp)
      DO i = 1, n
         chart%work(i, i) = chart%work(i, i) + 1.0_dp
         chart%inverse(i, i) = 1.0_dp
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

      INTEGER                                            :: a, b, i, j, n
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
         ! Form the occupation commutator before averaging to preserve the energy zero.
         DO j = 1, n
            DO i = 1, n
               IF (chart%occupations(i) + chart%occupations(j) == 0.0_dp) THEN
                  chart%transformed(i, j) = 0.0_dp
               ELSE
                  chart%transformed(i, j) = 2.0_dp*(chart%occupations(j) - chart%occupations(i))* &
                                            chart%transformed(i, j)/ &
                                            SQRT(chart%metric_weight*(chart%occupations(i) + chart%occupations(j)))
               END IF
            END DO
         END DO
      ELSE
         chart%transformed(:, :) = MATMUL(CONJG(TRANSPOSE(chart%base)), MATMUL(OPERATOR, chart%base))
         chart%work(:, :) = MATMUL(CONJG(TRANSPOSE(chart%rotation)), chart%transformed)
         DO i = 1, n
            chart%work(i, :) = chart%occupations(i)*chart%work(i, :)
         END DO
         chart%transformed(:, :) = MATMUL(chart%inverse, MATMUL(chart%work, chart%inverse))
         ! The Hermitian part has zero contraction with every anti-Hermitian generator.
         chart%transformed = 0.5_dp*(chart%transformed - CONJG(TRANSPOSE(chart%transformed)))
      END IF
      CALL orbital_symmetry_project(chart%symmetry, chart%antiunitary, chart%transformed, .TRUE.)
      DO i = 1, SIZE(derivative)
         a = chart%pairs(1, i)
         b = chart%pairs(2, i)
         IF (physical) THEN
            IF (chart%occupations(a) + chart%occupations(b) == 0.0_dp) THEN
               derivative(i) = 0.0_dp
               CYCLE
            END IF
            IF (chart%component(i) == 1) THEN
               value = REAL(chart%transformed(a, b), dp)
            ELSE
               value = AIMAG(chart%transformed(a, b))
            END IF
            derivative(i) = value
         ELSE IF (chart%component(i) == 1) THEN
            derivative(i) = 2.0_dp*REAL(chart%transformed(b, a) - chart%transformed(a, b), dp)
         ELSE
            derivative(i) = -2.0_dp*AIMAG(chart%transformed(b, a) + chart%transformed(a, b))
         END IF
         IF (.NOT. physical) derivative(i) = chart%coordinate_scale(i)*derivative(i)
      END DO
   END SUBROUTINE skala_orbital_chart_derivative

! **************************************************************************************************
!> \brief Differentiate the unweighted AO density in the evaluated original-metric Cayley chart.
!>        The charge metric and coordinate scales are fixed. The caller applies k weights once.
!>        Occupation directions are independent physical occupation changes, not band energies.
!> \param chart Evaluated orbital chart in the original normalization metric.
!> \param direction Real coordinate direction, including independent imaginary orbital rotations.
!> \param response Hermitian AO density response without a k-point weight.
!> \param occupations_direction Optional physical occupation direction for a variable-occupation chart.
!> \param orbital_direction Optional coefficient response for band-expectation derivatives.
! **************************************************************************************************
   SUBROUTINE skala_orbital_chart_density_response(chart, direction, response, occupations_direction, orbital_direction)
      TYPE(skala_orbital_chart_type), INTENT(IN)         :: chart
      REAL(dp), INTENT(IN)                               :: direction(:)
      COMPLEX(dp), INTENT(OUT)                           :: response(:, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: occupations_direction(:)
      COMPLEX(dp), INTENT(OUT), OPTIONAL                 :: orbital_direction(:, :)

      COMPLEX(dp), ALLOCATABLE                           :: coefficients_direction(:, :), &
                                                            generator(:, :), weighted(:, :)
      INTEGER                                            :: i, n, nao

      CPASSERT(chart%evaluated)
      n = SIZE(chart%base, 2)
      nao = SIZE(chart%base, 1)
      CPASSERT(ALL(SHAPE(response) == [nao, nao]))
      CALL orbital_generator_direction(chart, direction, generator)
      coefficients_direction = MATMUL(chart%base, MATMUL(chart%inverse, MATMUL(generator, chart%inverse)))
      IF (PRESENT(orbital_direction)) THEN
         CPASSERT(SIZE(orbital_direction, 1) == nao)
         CPASSERT(SIZE(orbital_direction, 2) == n)
         orbital_direction = coefficients_direction
      END IF
      weighted = coefficients_direction
      DO i = 1, n
         weighted(:, i) = chart%occupations(i)*weighted(:, i)
      END DO
      response = MATMUL(weighted, CONJG(TRANSPOSE(chart%coefficients)))
      response = response + CONJG(TRANSPOSE(response))
      IF (PRESENT(occupations_direction)) THEN
         CPASSERT(chart%variable_occupations)
         CPASSERT(SIZE(occupations_direction) == n)
         CPASSERT(ALL(ieee_is_finite(occupations_direction)))
         weighted = chart%coefficients
         DO i = 1, n
            weighted(:, i) = occupations_direction(i)*weighted(:, i)
         END DO
         response = response + MATMUL(weighted, CONJG(TRANSPOSE(chart%coefficients)))
      END IF
   END SUBROUTINE skala_orbital_chart_density_response

! **************************************************************************************************
!> \brief Differentiate the coordinate covector of an original Hermitian AO expectation.
!>        Retain both Cayley inverse responses and the curvature of the orbital map.
!>        The caller supplies identically weighted operators and their electronic response.
!>        Occupation and operator changes are retained together. The charge metric and
!>        coordinate scales stay fixed. Entropy and moving constraint normals remain separate.
!> \param chart Evaluated orbital chart.
!> \param direction Real coordinate direction.
!> \param OPERATOR Original Hermitian AO operator, including the caller's k-point weight.
!> \param operator_direction Response of the same weighted AO operator.
!> \param response Coordinate-covector response in the original fixed chart.
!> \param occupations_direction Optional physical occupation response.
! **************************************************************************************************
   SUBROUTINE skala_orbital_chart_operator_response(chart, direction, OPERATOR, operator_direction, response, &
                                                    occupations_direction)
      TYPE(skala_orbital_chart_type), INTENT(IN)         :: chart
      REAL(dp), INTENT(IN)                               :: direction(:)
      COMPLEX(dp), INTENT(IN)                            :: operator(:, :), operator_direction(:, :)
      REAL(dp), INTENT(OUT)                              :: response(:)
      REAL(dp), INTENT(IN), OPTIONAL                     :: occupations_direction(:)

      COMPLEX(dp), ALLOCATABLE :: generator(:, :), h0(:, :), h0_direction(:, :), &
         inverse_direction(:, :), transformed(:, :), transformed_direction(:, :), work(:, :)
      INTEGER                                            :: a, b, i, n, nao

      CPASSERT(chart%evaluated)
      n = SIZE(chart%base, 2)
      nao = SIZE(chart%base, 1)
      CPASSERT(ALL(SHAPE(OPERATOR) == [nao, nao]))
      CPASSERT(ALL(SHAPE(operator_direction) == [nao, nao]))
      CPASSERT(SIZE(response) == SIZE(chart%coordinate))
      CALL orbital_generator_direction(chart, direction, generator)
      inverse_direction = 0.5_dp*MATMUL(chart%inverse, MATMUL(generator, chart%inverse))
      h0 = MATMUL(CONJG(TRANSPOSE(chart%base)), MATMUL(OPERATOR, chart%base))
      h0_direction = MATMUL(CONJG(TRANSPOSE(chart%base)), MATMUL(operator_direction, chart%base))
      transformed = MATMUL(CONJG(TRANSPOSE(chart%rotation)), h0)
      transformed_direction = 2.0_dp*MATMUL(CONJG(TRANSPOSE(inverse_direction)), h0) + &
                              MATMUL(CONJG(TRANSPOSE(chart%rotation)), h0_direction)
      IF (PRESENT(occupations_direction)) THEN
         CPASSERT(chart%variable_occupations)
         CPASSERT(SIZE(occupations_direction) == n)
         CPASSERT(ALL(ieee_is_finite(occupations_direction)))
      END IF
      DO i = 1, n
         transformed_direction(i, :) = chart%occupations(i)*transformed_direction(i, :)
         IF (PRESENT(occupations_direction)) THEN
            transformed_direction(i, :) = transformed_direction(i, :) + occupations_direction(i)*transformed(i, :)
         END IF
         transformed(i, :) = chart%occupations(i)*transformed(i, :)
      END DO
      work = MATMUL(inverse_direction, MATMUL(transformed, chart%inverse)) + &
             MATMUL(chart%inverse, MATMUL(transformed_direction, chart%inverse)) + &
             MATMUL(chart%inverse, MATMUL(transformed, inverse_direction))
      CALL orbital_symmetry_project(chart%symmetry, chart%antiunitary, work, .TRUE.)
      DO i = 1, SIZE(response)
         a = chart%pairs(1, i)
         b = chart%pairs(2, i)
         IF (chart%component(i) == 1) THEN
            response(i) = 2.0_dp*REAL(work(b, a) - work(a, b), dp)
         ELSE
            response(i) = -2.0_dp*AIMAG(work(b, a) + work(a, b))
         END IF
         response(i) = chart%coordinate_scale(i)*response(i)
      END DO
   END SUBROUTINE skala_orbital_chart_operator_response

! **************************************************************************************************
!> \brief Form the anti-Hermitian generator direction with the fixed original coordinate scales.
!> \param chart ...
!> \param direction ...
!> \param generator ...
! **************************************************************************************************
   SUBROUTINE orbital_generator_direction(chart, direction, generator)
      TYPE(skala_orbital_chart_type), INTENT(IN)         :: chart
      REAL(dp), INTENT(IN)                               :: direction(:)
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: generator(:, :)

      COMPLEX(dp)                                        :: element
      INTEGER                                            :: a, b, i, n

      CPASSERT(SIZE(direction) == SIZE(chart%coordinate))
      CPASSERT(ALL(ieee_is_finite(direction)))
      n = SIZE(chart%base, 2)
      ALLOCATE (generator(n, n), SOURCE=CMPLX(0.0_dp, 0.0_dp, dp))
      DO i = 1, SIZE(direction)
         a = chart%pairs(1, i)
         b = chart%pairs(2, i)
         element = CMPLX(direction(i), 0.0_dp, dp)
         IF (chart%component(i) == 2) element = CMPLX(0.0_dp, direction(i), dp)
         element = chart%coordinate_scale(i)*element
         generator(a, b) = generator(a, b) + element
         generator(b, a) = generator(b, a) - CONJG(element)
      END DO
      CALL orbital_symmetry_project(chart%symmetry, chart%antiunitary, generator, .FALSE.)
   END SUBROUTINE orbital_generator_direction

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
!> \param check_stalled_segment Optional first-feasible-trial check after unsuccessful backtracking.
!>        Equal endpoints request an active-normal admissibility check before tangent minimization.
!>        or an unresolved energy change above the gradient tolerance.
!> \param coordinate_scales Optional constant diagonal map to original coordinates for stopping tests.
!>        Search directions use optimizer coordinates. The original unscaled gradient bound is retained.
!> \param curvature_steps Attempt coupled inexact Newton directions at every iteration, default false.
!> \param check_curvature_segment Optional branch check of every restored Hessian probe.
!> \param precondition Optional current positive diagonal inverse for coupled steps.
!>        The returned inverse stays fixed throughout each Newton solve.
!> \param apply_hessian Optional complete analytical response when no dependent constraints are present.
!>        An unavailable response retains the original numerical curvature and acceptance tests.
!> \param apply_lagrangian_hessian Optional analytical response on the original physical constraint curve.
!>        Constraint curvature and original-unit multipliers are included before tangent projection.
!> \param report_iteration Optional scalar-only report after reference evaluation and attempted steps.
! **************************************************************************************************
   SUBROUTINE skala_constrained_stationarity(evaluate_constraints, evaluate_gradient, context, coordinate, &
                                             density_tolerances, gradient_tolerance, converged, &
                                             max_iterations, max_step, verify_stationarity, iterations, check_stop, stopped, &
                                             check_segment, interrupted, retain_progress, check_stalled_segment, &
                                             coordinate_scales, curvature_steps, check_curvature_segment, precondition, &
                                             apply_hessian, apply_lagrangian_hessian, report_iteration)
      PROCEDURE(physical_constraints) :: evaluate_constraints
      PROCEDURE(orbital_energy_gradient) :: evaluate_gradient
      PROCEDURE(orbital_stationarity), OPTIONAL :: verify_stationarity
      PROCEDURE(orbital_stop), OPTIONAL :: check_stop
      PROCEDURE(orbital_segment), OPTIONAL :: check_segment
      PROCEDURE(orbital_segment), OPTIONAL :: check_stalled_segment
      PROCEDURE(orbital_curvature_segment), OPTIONAL :: check_curvature_segment
      PROCEDURE(orbital_preconditioner), OPTIONAL :: precondition
      PROCEDURE(orbital_hessian_product), OPTIONAL :: apply_hessian
      PROCEDURE(orbital_lagrangian_hessian), OPTIONAL :: apply_lagrangian_hessian
      PROCEDURE(orbital_iteration_report), OPTIONAL :: report_iteration
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
      REAL(dp), INTENT(IN), OPTIONAL :: coordinate_scales(:)
      LOGICAL, INTENT(IN), OPTIONAL :: curvature_steps

      INTEGER :: capacity, d, history, i, iteration, last_progress_check, limit, m, nf, r, shortening
      TYPE(skala_constraint_chart_type) :: curve
      LOGICAL :: accepted, valid, event, newton_valid, newton_stopped, have_segment, resolution_checked, &
                 coupled_curvature, have_curvature_check, curvature_interrupted, active_preconditioner, resolved_descent, &
                 curvature_stalled, newton_used
      REAL(dp) :: curvature, energy, fraction, progress_norm, radius, maximum_radius, roundoff, slope, trial_energy, trial_slope
      REAL(dp), ALLOCATABLE :: direction(:), g(:), jacobian(:, :), motion(:), newton(:), &
                               old_gradient(:), reduced(:), steps(:, :), trial(:), trial_reduced(:), &
                               values(:), x(:), changes(:, :), stalled_trial(:), &
                               inverse_diagonal(:), free_inverse(:)

      d = SIZE(coordinate)
      coupled_curvature = .FALSE.
      IF (PRESENT(curvature_steps)) coupled_curvature = curvature_steps
      have_curvature_check = PRESENT(check_curvature_segment)
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
      IF (PRESENT(coordinate_scales)) THEN
         CPASSERT(SIZE(coordinate_scales) == d)
         CPASSERT(ALL(ieee_is_finite(coordinate_scales)))
         CPASSERT(ALL(coordinate_scales > 0.0_dp))
      END IF
      limit = 200
      IF (PRESENT(max_iterations)) limit = max_iterations
      radius = 0.1_dp
      IF (PRESENT(max_step)) radius = max_step
      CPASSERT(limit > 0)
      CPASSERT(ieee_is_finite(radius) .AND. radius > 0.0_dp)
      maximum_radius = radius
      ALLOCATE (x(d), trial(d), values(m), jacobian(m, d), stalled_trial(d))
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
                old_gradient(nf), steps(nf, capacity), changes(nf, capacity), newton(nf), &
                inverse_diagonal(d), free_inverse(nf))
      CALL reduced_evaluation(x, energy, reduced, valid)
      IF (.NOT. valid) RETURN
      IF (PRESENT(report_iteration)) CALL report_iteration(0, energy, free_gradient_norm(reduced), .FALSE., context)
      history = 0
      progress_norm = NORM2(reduced)
      last_progress_check = 1
      DO iteration = 1, limit
         IF (stop_requested()) RETURN
         IF (stationary(x, reduced)) THEN
            coordinate(:) = x
            converged = .TRUE.
            RETURN
         END IF
         IF (PRESENT(iterations)) iterations = iteration
         IF (r > 0 .AND. PRESENT(check_stalled_segment)) THEN
            ! An active surface can stop being a normal minimum while its tangent
            ! still descends. Release it before spending steps on that tangent.
            event = .FALSE.
            CALL check_stalled_segment(x, x, curve, context, event)
            IF (event) THEN
               coordinate(:) = x
               IF (PRESENT(interrupted)) interrupted = .TRUE.
               RETURN
            END IF
            CALL reduced_evaluation(x, energy, reduced, valid)
            IF (.NOT. valid) RETURN
         END IF
         curvature_stalled = .FALSE.
         IF (iteration - last_progress_check >= MAX(8, capacity)) THEN
            curvature_stalled = NORM2(reduced) > 0.5_dp*progress_norm
            progress_norm = NORM2(reduced)
            last_progress_check = iteration
         END IF
         active_preconditioner = .FALSE.
         IF (coupled_curvature .AND. PRESENT(precondition)) THEN
            CALL precondition(x, context, inverse_diagonal, active_preconditioner)
            IF (active_preconditioner) THEN
               CPASSERT(ALL(ieee_is_finite(inverse_diagonal)))
               CPASSERT(ALL(inverse_diagonal > 0.0_dp))
               free_inverse = inverse_diagonal(curve%free)
               ! Nearly empty bands can yield an unresolvable inverse response. Bound
               ! only the search metric, leaving the original energy and stopping tests.
               free_inverse = MIN(free_inverse, MIN(1.0_dp, MINVAL(free_inverse))/SQRT(EPSILON(1.0_dp)))
            END IF
         END IF
         IF (active_preconditioner) THEN
            CALL bfgs_direction(reduced, steps, changes, history, direction, diagonal=free_inverse)
         ELSE
            CALL bfgs_direction(reduced, steps, changes, history, direction)
         END IF
         IF (DOT_PRODUCT(reduced, direction) >= 0.0_dp) THEN
            history = 0
            direction(:) = -reduced
         END IF
         roundoff = 64.0_dp*EPSILON(0.0_dp)*MAX(ABS(energy), TINY(0.0_dp))
         resolved_descent = .FALSE.
         newton_used = .FALSE.
         IF (coupled_curvature .OR. curvature_stalled .OR. (history > 0 .AND. &
                                                            -DOT_PRODUCT(reduced, direction) <= 32.0_dp*roundoff)) THEN
            IF (PRESENT(check_stalled_segment)) THEN
               ! Resolve a physical boundary on the feasible BFGS segment before
               ! estimating smooth curvature. The callback retains its admission tests.
               fraction = 1.0_dp
               IF (NORM2(direction) > radius) fraction = radius/NORM2(direction)
               trial(:) = x
               trial(curve%free) = x(curve%free) + fraction*direction
               IF (ANY(trial(curve%free) /= x(curve%free))) THEN
                  CALL reduced_evaluation(trial, trial_energy, trial_reduced, valid)
                  IF (valid) THEN
                     roundoff = 64.0_dp*EPSILON(0.0_dp)*MAX(ABS(energy), ABS(trial_energy), TINY(0.0_dp))
                     resolved_descent = trial_energy < energy - roundoff .AND. &
                                        trial_energy <= energy + 1.0e-4_dp*fraction*DOT_PRODUCT(reduced, direction)
                     IF (.NOT. resolved_descent) THEN
                        event = .FALSE.
                        CALL check_stalled_segment(x, trial, curve, context, event)
                        IF (event) THEN
                           coordinate(:) = x
                           IF (PRESENT(interrupted)) interrupted = .TRUE.
                           RETURN
                        END IF
                     END IF
                  END IF
                  ! A boundary scan can restore the native fields at a different point.
                  CALL reduced_evaluation(x, energy, reduced, valid)
                  IF (.NOT. valid) RETURN
               END IF
            END IF
            ! A resolved descent can increase the gradient in a region of negative curvature.
            ! Keep that original BFGS proposal for the ordinary acceptance checks below.
            ! Probe curvature only when that proposal has not already established descent.
            ! The extra probe is bounded to one per history-length interval. Every trial
            ! restores the same physical constraints and uses the original line search.
            IF (coupled_curvature .OR. .NOT. resolved_descent) THEN
               curvature_interrupted = .FALSE.
               CALL newton_direction(x, reduced, newton, newton_valid, newton_stopped)
               IF (curvature_interrupted) THEN
                  coordinate(:) = x
                  IF (PRESENT(interrupted)) interrupted = .TRUE.
                  RETURN
               END IF
               IF (newton_stopped) RETURN
               IF (newton_valid) THEN
                  direction(:) = newton
                  newton_used = .TRUE.
               END IF
            END IF
         END IF
         IF (NORM2(direction) > radius) direction(:) = direction*(radius/NORM2(direction))
         slope = DOT_PRODUCT(reduced, direction)
         accepted = .FALSE.
         fraction = 1.0_dp
         have_segment = .FALSE.
         resolution_checked = .FALSE.
         DO shortening = 1, DIGITS(1.0_dp)
            IF (stop_requested()) RETURN
            trial(:) = x
            trial(curve%free) = x(curve%free) + fraction*direction
            IF (ALL(trial(curve%free) == x(curve%free))) EXIT
            CALL reduced_evaluation(trial, trial_energy, trial_reduced, valid)
            IF (valid) THEN
               IF (.NOT. have_segment) THEN
                  stalled_trial(:) = trial
                  have_segment = .TRUE.
               END IF
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
               IF (PRESENT(check_stalled_segment) .AND. .NOT. resolution_checked .AND. &
                   ((curvature_stalled .AND. .NOT. resolved_descent) .OR. ABS(trial_energy - energy) <= roundoff .OR. &
                    (trial_slope > 0.9_dp*ABS(slope) .AND. trial_energy < energy - roundoff)) .AND. &
                   free_gradient_norm(trial_reduced) > gradient_tolerance) THEN
                  ! Inspect a reversal that fails the usual strong curvature condition.
                  ! A reversal satisfying it can be an ordinary smooth descent step.
                  resolution_checked = .TRUE.
                  event = .FALSE.
                  IF (trial_slope >= 0.0_dp .AND. trial_energy < energy - roundoff) THEN
                     CALL check_stalled_segment(x, trial, curve, context, event)
                  ELSE
                     CALL check_stalled_segment(x, stalled_trial, curve, context, event)
                  END IF
                  IF (event) THEN
                     coordinate(:) = x
                     IF (PRESENT(interrupted)) interrupted = .TRUE.
                     RETURN
                  END IF
                  ! A scan can leave native fields at x. Restore the actual trial and retain
                  ! its original energy/derivative acceptance tests even when no event is found.
                  CALL reduced_evaluation(trial, trial_energy, trial_reduced, valid)
                  IF (.NOT. valid) THEN
                     fraction = 0.5_dp*fraction
                     CYCLE
                  END IF
                  roundoff = 64.0_dp*EPSILON(0.0_dp)*MAX(ABS(energy), ABS(trial_energy), TINY(0.0_dp))
                  trial_slope = DOT_PRODUCT(trial_reduced, direction)
               END IF
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
            IF (PRESENT(report_iteration)) CALL report_iteration(iteration, energy, free_gradient_norm(reduced), .FALSE., context)
            IF (PRESENT(check_stalled_segment) .AND. have_segment) THEN
               ! A rejected wide trial can enclose a downhill cusp. Keep its original
               ! feasible endpoints rather than the last nearly coincident shortened trial.
               event = .FALSE.
               CALL check_stalled_segment(x, stalled_trial, curve, context, event)
               IF (event) THEN
                  coordinate(:) = x
                  IF (PRESENT(interrupted)) interrupted = .TRUE.
                  RETURN
               END IF
            END IF
            IF (history == 0) EXIT
            ! Retry with steepest descent when the stored curvature gives no acceptable step.
            history = 0
            CYCLE
         END IF
         motion(:) = trial(curve%free) - x(curve%free)
         IF (coupled_curvature .AND. newton_used) THEN
            ! A shortened step resolves a smaller neighborhood than the original CG sphere.
            ! Use that neighborhood in the next solve while retaining the user maximum.
            IF (fraction < 0.25_dp) THEN
               radius = MIN(radius, 2.0_dp*NORM2(motion))
            ELSE IF (fraction == 1.0_dp .AND. NORM2(motion) >= 0.9_dp*radius) THEN
               radius = MIN(maximum_radius, 2.0_dp*radius)
            END IF
         END IF
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
         IF (PRESENT(report_iteration)) CALL report_iteration(iteration, energy, free_gradient_norm(reduced), .TRUE., context)
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
!>        Truncated conjugate gradients reuse the positive BFGS inverse as a fixed preconditioner
!>        and retain linear storage. The constraint Jacobian and normal derivative are
!>        retained only for this accepted center. Indefinite curvature and an oversized CG iterate
!>        terminate at the current CG sphere, bounded by the user maximum. The outer line search
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
      REAL(dp)                                           :: alpha, beta, boundary, curvature, &
                                                            next_residual2, next_residual_dot, &
                                                            norm_direction, residual_dot, TARGET
      REAL(dp), ALLOCATABLE                              :: center_jacobian(:, :), center_normal(:), &
                                                            preconditioned(:), PRODUCT(:), &
                                                            residual(:), search(:)

         valid_direction = .FALSE.
         stopped_direction = .FALSE.
         newton_step = 0.0_dp
         ALLOCATE (residual(nf), search(nf), PRODUCT(nf), preconditioned(nf))
         residual(:) = -derivative
         IF (active_preconditioner) THEN
            CALL bfgs_direction(residual, steps, changes, history, preconditioned, diagonal=free_inverse)
         ELSE
            CALL bfgs_direction(residual, steps, changes, history, preconditioned)
         END IF
         preconditioned(:) = -preconditioned
         residual_dot = DOT_PRODUCT(residual, preconditioned)
         IF (.NOT. ieee_is_finite(residual_dot) .OR. residual_dot <= 0.0_dp) RETURN
         search(:) = preconditioned
         ! Match the inner solve to the one-percent resolution of numerical Hessian products.
         ! The original unscaled physical gradient still determines final convergence.
         TARGET = MAX(0.01_dp, MIN(0.1_dp, SQRT(NORM2(derivative))))*NORM2(derivative)
         ! Finite-precision CG can need more than nf steps after conjugacy is lost.
         DO iteration = 1, MIN(2*nf, 64)
            IF (stop_requested()) THEN
               stopped_direction = .TRUE.
               RETURN
            END IF
            norm_direction = NORM2(search)
            IF (norm_direction <= TINY(0.0_dp)) EXIT
            CALL reduced_hessian_product(center, search, PRODUCT, valid_trial, stopped_direction, &
                                         center_jacobian, center_normal)
            IF (.NOT. valid_trial .OR. stopped_direction) RETURN
            curvature = DOT_PRODUCT(search, product)
            boundary = boundary_length(newton_step, search, norm_direction)
            IF (curvature <= 128.0_dp*EPSILON(0.0_dp)*norm_direction*NORM2(product)) THEN
               newton_step = newton_step + boundary*(search/norm_direction)
               EXIT
            END IF
            alpha = residual_dot/curvature
            IF (alpha >= boundary/norm_direction) THEN
               newton_step = newton_step + boundary*(search/norm_direction)
               EXIT
            END IF
            newton_step(:) = newton_step + alpha*search
            residual(:) = residual - alpha*product
            next_residual2 = DOT_PRODUCT(residual, residual)
            IF (SQRT(next_residual2) <= TARGET) EXIT
            IF (active_preconditioner) THEN
               CALL bfgs_direction(residual, steps, changes, history, preconditioned, diagonal=free_inverse)
            ELSE
               CALL bfgs_direction(residual, steps, changes, history, preconditioned)
            END IF
            preconditioned(:) = -preconditioned
            next_residual_dot = DOT_PRODUCT(residual, preconditioned)
            IF (.NOT. ieee_is_finite(next_residual_dot) .OR. next_residual_dot <= 0.0_dp) EXIT
            beta = next_residual_dot/residual_dot
            search(:) = preconditioned + beta*search
            residual_dot = next_residual_dot
         END DO
         valid_direction = ALL(ieee_is_finite(newton_step)) .AND. DOT_PRODUCT(derivative, newton_step) < 0.0_dp
      END SUBROUTINE newton_direction
! **************************************************************************************************
!> \brief Distance along a normalized CG direction to the current CG sphere.
!>        Normalize by the radius and avoid cancellation at an outward-facing boundary.
!> \param step Current interior Newton step.
!> \param direction Current CG search direction.
!> \param norm_direction Positive norm of that direction.
!> \return Nonnegative distance to the sphere along the same direction.
! **************************************************************************************************
      FUNCTION boundary_length(step, direction, norm_direction) RESULT(length)
      REAL(dp), INTENT(IN)                               :: step(:), direction(:), norm_direction
      REAL(dp)                                           :: length

      REAL(dp)                                           :: projection, root, slack

         projection = DOT_PRODUCT(step/radius, direction/norm_direction)
         slack = MAX(0.0_dp, 1.0_dp - SUM((step/radius)**2))
         root = SQRT(projection**2 + slack)
         length = 0.0_dp
         IF (projection < 0.0_dp) THEN
            length = radius*(root - projection)
         ELSE IF (root + projection > 0.0_dp) THEN
            length = radius*(slack/(root + projection))
         END IF
      END FUNCTION boundary_length
! **************************************************************************************************
!> \brief Use the complete analytical unconstrained response when supplied, otherwise resolve native curvature.
!>        Dependent physical constraints retain their existing implicit-curve response.
!>        Refine central differences to one-percent agreement. If coordinate resolution prevents
!>        agreement, retain the estimate with the smallest relative change. This constructs a trial
!>        direction. The original energy, constraint and stationarity tests remain authoritative.
!> \param center Accepted orbital coordinates.
!> \param direction Free-coordinate direction.
!> \param PRODUCT Estimated reduced Hessian applied to direction.
!> \param valid_product Whether finite native differences supplied an estimate.
!> \param stopped_product Whether the caller requested a stop.
!> \param center_jacobian Accepted independent constraint rows, scaled once per Newton solve.
!> \param center_normal Complete normal energy derivative at the same center.
! **************************************************************************************************
      SUBROUTINE reduced_hessian_product(center, direction, PRODUCT, valid_product, stopped_product, &
                                         center_jacobian, center_normal)
      REAL(dp), INTENT(IN)                               :: center(:), direction(:)
      REAL(dp), INTENT(OUT)                              :: PRODUCT(:)
      LOGICAL, INTENT(OUT)                               :: valid_product, stopped_product
      REAL(dp), ALLOCATABLE, INTENT(INOUT)               :: center_jacobian(:, :), center_normal(:)

      INTEGER                                            :: refinement
      LOGICAL                                            :: have_previous, have_product, smooth, &
                                                            valid_trial
      REAL(dp)                                           :: best_error, difference, energy_trial, &
                                                            relative_error, scale, step
      REAL(dp), ALLOCATABLE :: candidate(:), full_direction(:), full_product(:), &
         gradient_minus(:), gradient_plus(:), matrix(:, :), minus(:), multipliers(:), plus(:), &
         previous(:), rhs(:, :)

         valid_product = .FALSE.
         stopped_product = .FALSE.
         PRODUCT = 0.0_dp
         scale = NORM2(direction)
         IF (scale <= TINY(0.0_dp)) RETURN
         IF (PRESENT(apply_hessian)) THEN
            IF (r == 0) THEN
               IF (stop_requested()) THEN
                  stopped_product = .TRUE.
                  RETURN
               END IF
               ALLOCATE (full_direction(d), full_product(d))
               full_direction = 0.0_dp
               full_direction(curve%free) = direction
               CALL apply_hessian(center, full_direction, context, full_product, valid_product)
               IF (valid_product) THEN
                  CPASSERT(ALL(ieee_is_finite(full_product)))
                  PRODUCT = full_product(curve%free)
                  RETURN
               END IF
            END IF
         END IF
         IF (PRESENT(apply_lagrangian_hessian) .AND. r > 0) THEN
            IF (stop_requested()) THEN
               stopped_product = .TRUE.
               RETURN
            END IF
            IF (.NOT. ALLOCATED(center_jacobian)) THEN
               CALL evaluate_constraints(center, values, jacobian, context)
               IF (ALL(ABS(values) <= curve%tolerances)) THEN
                  ALLOCATE (center_jacobian(r, d), center_normal(r))
                  center_jacobian = jacobian(curve%rows(:r), :)/SPREAD(curve%scales(curve%rows(:r)), 2, d)
                  ! Backtracking can leave a rejected full gradient cached. Refresh the
                  ! accepted center once, before any of its conjugate-gradient products.
                  CALL evaluate_gradient(center, energy_trial, g, context)
                  center_normal = g(curve%dependent)
               END IF
            END IF
            IF (ALLOCATED(center_jacobian)) THEN
               ALLOCATE (full_direction(d), full_product(d), matrix(r, r), rhs(r, 1), multipliers(m))
               matrix = center_jacobian(:, curve%dependent)
               rhs(:, 1) = -MATMUL(center_jacobian(:, curve%free), direction)
               CALL solve_constraint_system(matrix, rhs, valid_trial)
               IF (valid_trial) THEN
                  full_direction(curve%free) = direction
                  full_direction(curve%dependent) = rhs(:, 1)
                  matrix = TRANSPOSE(center_jacobian(:, curve%dependent))
                  rhs(:, 1) = center_normal
                  CALL solve_constraint_system(matrix, rhs, valid_trial)
                  IF (valid_trial) THEN
                     multipliers = 0.0_dp
                     multipliers(curve%rows(:r)) = rhs(:, 1)/curve%scales(curve%rows(:r))
                     CALL apply_lagrangian_hessian(center, full_direction, multipliers, context, &
                                                   full_product, valid_product)
                     IF (valid_product) THEN
                        CPASSERT(ALL(ieee_is_finite(full_product)))
                        matrix = TRANSPOSE(center_jacobian(:, curve%dependent))
                        rhs(:, 1) = full_product(curve%dependent)
                        CALL solve_constraint_system(matrix, rhs, valid_product)
                        IF (valid_product) THEN
                           PRODUCT = full_product(curve%free) - &
                                     MATMUL(TRANSPOSE(center_jacobian(:, curve%free)), rhs(:, 1))
                           RETURN
                        END IF
                     END IF
                  END IF
               END IF
            END IF
         END IF
         ! Central differences balance quadratic truncation against derivative roundoff.
         step = EPSILON(0.0_dp)**(1.0_dp/3.0_dp)*MAX(1.0_dp, NORM2(center(curve%free)))/scale
         best_error = HUGE(0.0_dp)
         have_previous = .FALSE.
         have_product = .FALSE.
         ALLOCATE (plus(d), minus(d), gradient_plus(nf), gradient_minus(nf), candidate(nf), previous(nf))
         DO refinement = 1, DIGITS(0.0_dp)
            plus(:) = center
            minus(:) = center
            plus(curve%free) = center(curve%free) + step*direction
            minus(curve%free) = center(curve%free) - step*direction
            IF (ALL(plus(curve%free) == center(curve%free)) .OR. &
                ALL(minus(curve%free) == center(curve%free))) EXIT
            IF (stop_requested()) THEN
               stopped_product = .TRUE.
               RETURN
            END IF
            CALL reduced_evaluation(plus, energy_trial, gradient_plus, valid_trial, center, smooth)
            IF (curvature_interrupted) RETURN
            IF (have_curvature_check) THEN
               IF (.NOT. smooth) THEN
                  have_previous = .FALSE.
                  have_product = .FALSE.
                  best_error = HUGE(0.0_dp)
                  step = 0.5_dp*step
                  CYCLE
               END IF
            END IF
            IF (.NOT. valid_trial) RETURN
            IF (stop_requested()) THEN
               stopped_product = .TRUE.
               RETURN
            END IF
            CALL reduced_evaluation(minus, energy_trial, gradient_minus, valid_trial, center, smooth)
            IF (curvature_interrupted) RETURN
            IF (have_curvature_check) THEN
               IF (.NOT. smooth) THEN
                  have_previous = .FALSE.
                  have_product = .FALSE.
                  best_error = HUGE(0.0_dp)
                  step = 0.5_dp*step
                  CYCLE
               END IF
            END IF
            IF (.NOT. valid_trial) RETURN
            candidate(:) = (gradient_plus - gradient_minus)/(2.0_dp*step)
            IF (.NOT. ALL(ieee_is_finite(candidate))) RETURN
            IF (.NOT. have_previous) THEN
               PRODUCT(:) = candidate
            ELSE
               have_product = .TRUE.
               difference = NORM2(candidate - previous)
               relative_error = difference/MAX(NORM2(candidate), NORM2(previous), TINY(0.0_dp))
               IF (relative_error < best_error) THEN
                  best_error = relative_error
                  PRODUCT(:) = candidate
               END IF
               IF (relative_error <= 0.01_dp) EXIT
            END IF
            previous(:) = candidate
            have_previous = .TRUE.
            step = 0.5_dp*step
         END DO
         valid_product = have_product
      END SUBROUTINE reduced_hessian_product
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

         converged_trial = free_gradient_norm(rg) <= gradient_tolerance
         IF (.NOT. converged_trial) RETURN
         IF (PRESENT(verify_stationarity)) CALL verify_stationarity(c, context, converged_trial)
      END FUNCTION stationary
! **************************************************************************************************
!> \brief Test stationarity in the original coordinates independently of preconditioning.
!> \param derivative Derivative in the optimizer's free coordinates.
!> \return Norm in the original unscaled free-coordinate units.
! **************************************************************************************************
      FUNCTION free_gradient_norm(derivative) RESULT(value)
      REAL(dp), INTENT(IN)                               :: derivative(:)
      REAL(dp)                                           :: value

         IF (PRESENT(coordinate_scales)) THEN
            value = NORM2(derivative/coordinate_scales(curve%free))
         ELSE
            value = NORM2(derivative)
         END IF
      END FUNCTION free_gradient_norm
! **************************************************************************************************
!> \brief Restore the physical targets and eliminate the dependent-coordinate derivatives.
!> \param c Trial coordinates, including restored dependent coordinates.
!> \param e Complete energy at the restored trial.
!> \param rg Energy derivative in the independent coordinates.
!> \param success Whether the constraints and reduced derivative were evaluated successfully.
!> \param probe_origin Optional center for a field-only branch check before requesting derivatives.
!> \param probe_smooth Whether that check permits this probe. No derivative is returned on rejection.
! **************************************************************************************************
      SUBROUTINE reduced_evaluation(c, e, rg, success, probe_origin, probe_smooth)
      REAL(dp), INTENT(INOUT)                            :: c(:)
      REAL(dp), INTENT(OUT)                              :: e, rg(:)
      LOGICAL, INTENT(OUT)                               :: success
      REAL(dp), INTENT(IN), OPTIONAL                     :: probe_origin(:)
      LOGICAL, INTENT(OUT), OPTIONAL                     :: probe_smooth

      LOGICAL                                            :: solved
      REAL(dp), ALLOCATABLE                              :: matrix(:, :), multiplier(:, :)

         IF (PRESENT(probe_smooth)) probe_smooth = .TRUE.
         CALL skala_constraint_chart_restore(curve, evaluate_constraints, context, c, values, jacobian, success)
         IF (.NOT. success) RETURN
         success = .FALSE.
         IF (PRESENT(probe_origin) .AND. have_curvature_check) THEN
            CPASSERT(PRESENT(probe_smooth))
            CALL check_curvature_segment(probe_origin, c, curve, context, probe_smooth, curvature_interrupted)
            IF (curvature_interrupted .OR. .NOT. probe_smooth) RETURN
         END IF
         ! A branch check can restore other native fields. Rebuild this trial's gradient afterwards.
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
            ! A small field residual can still displace a weak-normal implicit curve.
            ! Refine the target while representable Newton corrections reduce that residual.
            before = MAXVAL(ABS(values/chart%scales))
            IF (before == 0.0_dp) EXIT
            matrix(:, :) = jacobian(chart%rows(:r), chart%dependent)/SPREAD(chart%scales(chart%rows(:r)), 2, r)
            rhs(:, 1) = -values(chart%rows(:r))/chart%scales(chart%rows(:r))
            CALL solve_constraint_system(matrix, rhs, solved)
            IF (.NOT. solved) EXIT
            fraction = 1.0_dp
            DO k = 1, 24
               candidate(:) = coordinate
               candidate(chart%dependent) = coordinate(chart%dependent) + fraction*rhs(:, 1)
               IF (ALL(candidate == coordinate)) EXIT
               CALL evaluate_constraints(candidate, values, jacobian, context)
               CPASSERT(ALL(ieee_is_finite(values)))
               CPASSERT(ALL(ieee_is_finite(jacobian)))
               IF (MAXVAL(ABS(values/chart%scales)) < before) EXIT
               fraction = 0.5_dp*fraction
            END DO
            IF (k > 24 .OR. ALL(candidate == coordinate)) THEN
               ! Rejected trial callbacks must not leave fields at a different coordinate.
               CALL evaluate_constraints(coordinate, values, jacobian, context)
               EXIT
            END IF
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
!> \param group_ids Optional grouping of proportional normals, valid only on success.
! **************************************************************************************************
   SUBROUTINE skala_cusp_admissibility(positive_gradient, negative_changes, jacobian, &
                                       admissible, release_side, valid, group_ids)
      REAL(dp), INTENT(IN)                               :: positive_gradient(:), &
                                                            negative_changes(:, :), jacobian(:, :)
      LOGICAL, INTENT(OUT)                               :: admissible(:)
      INTEGER, INTENT(OUT)                               :: release_side(:)
      LOGICAL, INTENT(OUT)                               :: valid
      INTEGER, INTENT(OUT), OPTIONAL                     :: group_ids(:)

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
      IF (PRESENT(group_ids)) THEN
         CPASSERT(SIZE(group_ids) == m)
         group_ids = 0
      END IF
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
      IF (PRESENT(group_ids)) group_ids = groups
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
!> \param diagonal ...
! **************************************************************************************************
   SUBROUTINE bfgs_direction(gradient, steps, changes, history, direction, diagonal)
      REAL(dp), INTENT(IN)                               :: gradient(:), steps(:, :), changes(:, :)
      INTEGER, INTENT(IN)                                :: history
      REAL(dp), INTENT(OUT)                              :: direction(:)
      REAL(dp), INTENT(IN), OPTIONAL                     :: diagonal(:)

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
      IF (PRESENT(diagonal)) THEN
         CPASSERT(SIZE(diagonal) == SIZE(direction))
         direction(:) = diagonal*direction
      ELSE
         direction(:) = scale*direction
      END IF
      DO i = 1, history
         beta = inverse_curvature(i)*DOT_PRODUCT(changes(:, i), direction)
         direction(:) = direction + steps(:, i)*(alpha(i) - beta)
      END DO
      direction(:) = -direction
   END SUBROUTINE bfgs_direction
END MODULE skala_orbital_constraints
