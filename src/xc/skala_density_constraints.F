!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Own selected physical density and kinetic-density subgradients of one KS environment.
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
   PUBLIC :: skala_density_samples_begin, skala_density_samples_end, skala_density_samples_active, &
             skala_density_samples_put, skala_density_samples_get, skala_density_samples_only, &
             skala_density_samples_generation
   PUBLIC :: skala_density_constraints_update_weights, skala_density_constraints_clear_selection
   PUBLIC :: skala_density_jumps_begin, skala_density_jumps_active, &
             skala_density_jumps_put, skala_density_jumps_finish
   INTEGER, PARAMETER, PUBLIC :: skala_density_descriptor = 1, skala_kinetic_descriptor = 2
   INTEGER, PARAMETER, PUBLIC :: skala_fft_domain = 0, skala_atom_domain = 1, skala_background_domain = 2
   INTEGER, PARAMETER, PUBLIC :: skala_paw_ae_domain = 3, skala_paw_ps_domain = 4
   INTEGER, PARAMETER, PUBLIC :: skala_last_sample_domain = skala_paw_ps_domain

   ABSTRACT INTERFACE
! **************************************************************************************************
!> \brief ...
!> \param step ...
!> \param value ...
!> \param context ...
! **************************************************************************************************
      SUBROUTINE evaluate_physical_density(step, value, context)
         IMPORT :: dp
         REAL(dp), INTENT(IN) :: step
         REAL(dp), INTENT(OUT) :: value
         CLASS(*), INTENT(INOUT) :: context
      END SUBROUTINE evaluate_physical_density
   END INTERFACE

   TYPE skala_density_samples_type
      INTEGER, ALLOCATABLE :: points(:, :)
      REAL(dp), ALLOCATABLE :: density(:, :), kinetic(:, :)
   END TYPE skala_density_samples_type

   TYPE skala_density_constraints_type
      PRIVATE
      INTEGER :: nspins = 0, npts(3) = 0, bounds(2, 3) = 0
      INTEGER, ALLOCATABLE :: points(:, :), spins(:), descriptors(:), domains(:), atomic_sizes(:)
      INTEGER, ALLOCATABLE :: paw_sizes(:, :)
      REAL(dp) :: hmat(3, 3) = 0.0_dp, probe = 0.0_dp, tolerance = 0.0_dp
      REAL(dp) :: kinetic_probe = 0.0_dp, kinetic_tolerance = 0.0_dp
      REAL(dp), ALLOCATABLE :: positions(:, :), negative_weights(:), targets(:)
      LOGICAL :: sampling = .FALSE., descriptors_only = .FALSE.
      LOGICAL :: collecting_jumps = .FALSE.
      REAL(dp), ALLOCATABLE :: jump_values(:)
      INTEGER, ALLOCATABLE :: jump_owners(:)
      INTEGER :: sample_nspins = 0
      INTEGER :: sample_generation = 0
      TYPE(skala_density_samples_type) :: samples(skala_fft_domain:skala_last_sample_domain)
   END TYPE skala_density_constraints_type

CONTAINS

! **************************************************************************************************
!> \brief Collect primitive negative-minus-positive covectors during one complete KS evaluation.
!> \param state Registered physical constraints.
! **************************************************************************************************
   SUBROUTINE skala_density_jumps_begin(state)
      TYPE(skala_density_constraints_type), &
         INTENT(INOUT)                                   :: state

      INTEGER                                            :: n

      CALL clear_jumps(state)
      n = skala_density_constraints_count(state)
      ALLOCATE (state%jump_values(n), state%jump_owners(n))
      state%jump_values = 0.0_dp
      state%jump_owners = 0
      state%collecting_jumps = .TRUE.
   END SUBROUTINE skala_density_jumps_begin

! **************************************************************************************************
!> \brief Whether this native evaluation requests primitive covector differences.
!> \param state Registered physical constraints.
!> \return Whether collection is active.
! **************************************************************************************************
   FUNCTION skala_density_jumps_active(state) RESULT(active)
      TYPE(skala_density_constraints_type), INTENT(IN)   :: state
      LOGICAL                                            :: active

      active = state%collecting_jumps
   END FUNCTION skala_density_jumps_active

! **************************************************************************************************
!> \brief Record already weighted covectors only on their physical sample owner.
!> \param state Registered physical constraints.
!> \param domain Consumer domain.
!> \param values Covector differences in the domain's constraint order.
!> \param owned Locally evaluated physical samples in the same order.
! **************************************************************************************************
   SUBROUTINE skala_density_jumps_put(state, domain, values, owned)
      TYPE(skala_density_constraints_type), &
         INTENT(INOUT)                                   :: state
      INTEGER, INTENT(IN)                                :: domain
      REAL(dp), INTENT(IN)                               :: values(:)
      LOGICAL, INTENT(IN)                                :: owned(:)

      INTEGER                                            :: i, j

      IF (.NOT. state%collecting_jumps) RETURN
      CPASSERT(SIZE(values) == skala_density_constraints_count(state, domain))
      CPASSERT(SIZE(owned) == SIZE(values))
      j = 0
      DO i = 1, skala_density_constraints_count(state)
         IF (state%domains(i) /= domain) CYCLE
         j = j + 1
         IF (.NOT. owned(j)) CYCLE
         CPASSERT(ieee_is_finite(values(j)))
         state%jump_values(i) = state%jump_values(i) + values(j)
         state%jump_owners(i) = state%jump_owners(i) + 1
      END DO
   END SUBROUTINE skala_density_jumps_put

! **************************************************************************************************
!> \brief Consume this evaluation's local coefficients. The caller reduces values and owners.
!> \param state Registered physical constraints.
!> \param values Local covector differences, zero if collection was invalidated.
!> \param owners Number of local providers for each physical sample.
! **************************************************************************************************
   SUBROUTINE skala_density_jumps_finish(state, values, owners)
      TYPE(skala_density_constraints_type), &
         INTENT(INOUT)                                   :: state
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: values(:)
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: owners(:)

      INTEGER                                            :: n

      n = skala_density_constraints_count(state)
      ALLOCATE (values(n), owners(n))
      values = 0.0_dp
      owners = 0
      IF (state%collecting_jumps) THEN
         CPASSERT(SIZE(state%jump_values) == n)
         values = state%jump_values
         owners = state%jump_owners
      END IF
      CALL clear_jumps(state)
   END SUBROUTINE skala_density_jumps_finish

! **************************************************************************************************
!> \brief Invalidate transient coefficients when the selected physical samples change.
!> \param state Registered physical constraints.
! **************************************************************************************************
   SUBROUTINE clear_jumps(state)
      TYPE(skala_density_constraints_type), &
         INTENT(INOUT)                                   :: state

      IF (ALLOCATED(state%jump_values)) DEALLOCATE (state%jump_values, state%jump_owners)
      state%collecting_jumps = .FALSE.
   END SUBROUTINE clear_jumps

! **************************************************************************************************
!> \brief Start one native evaluation, discarding the previous physical descriptor samples.
!>        Consumers record their unchanged model inputs on their existing MPI owners.
!>        RKS stores one column containing the half density and half kinetic density.
!>        Every trial must begin a new sample. No orbital history or root is inferred here.
!> \param state ...
!> \param descriptors_only Record the native fields without evaluating XC, false by default.
! **************************************************************************************************
   SUBROUTINE skala_density_samples_begin(state, descriptors_only)
      TYPE(skala_density_constraints_type), &
         INTENT(INOUT)                                   :: state
      LOGICAL, INTENT(IN), OPTIONAL                      :: descriptors_only

      CALL clear_samples(state)
      state%descriptors_only = .FALSE.
      IF (PRESENT(descriptors_only)) state%descriptors_only = descriptors_only
      state%sampling = .TRUE.
   END SUBROUTINE skala_density_samples_begin

! **************************************************************************************************
!> \brief Freeze the descriptors consumed by the native evaluation.
!> \param state ...
! **************************************************************************************************
   SUBROUTINE skala_density_samples_end(state)
      TYPE(skala_density_constraints_type), &
         INTENT(INOUT)                                   :: state

      state%sampling = .FALSE.
      state%descriptors_only = .FALSE.
   END SUBROUTINE skala_density_samples_end

! **************************************************************************************************
!> \brief Identify the current descriptor collection independently of selected derivatives.
!>        Starting another collection or clearing samples changes this generation even when
!>        the orbital coordinates are unchanged. Empty distributed owners retain the same token.
!> \param state Native descriptor collection.
!> \return Collection generation, without interpreting local sample availability.
! **************************************************************************************************
   FUNCTION skala_density_samples_generation(state) RESULT(generation)
      TYPE(skala_density_constraints_type), INTENT(IN)   :: state
      INTEGER                                            :: generation

      generation = state%sample_generation
   END FUNCTION skala_density_samples_generation

! **************************************************************************************************
!> \brief Whether the native consumers should retain their current physical inputs.
!> \param state ...
!> \return ...
! **************************************************************************************************
   FUNCTION skala_density_samples_active(state) RESULT(active)
      TYPE(skala_density_constraints_type), INTENT(IN)   :: state
      LOGICAL                                            :: active

      active = state%sampling
   END FUNCTION skala_density_samples_active

! **************************************************************************************************
!> \brief Whether this traversal requests physical descriptors without an XC energy or operator.
!> \param state ...
!> \return ...
! **************************************************************************************************
   FUNCTION skala_density_samples_only(state) RESULT(ONLY)
      TYPE(skala_density_constraints_type), INTENT(IN)   :: state
      LOGICAL                                            :: only

      ONLY = state%sampling .AND. state%descriptors_only
   END FUNCTION skala_density_samples_only

! **************************************************************************************************
!> \brief Record unchanged local rho and tau inputs before nonlinear model features.
!>        Sampling is opt-in. Updating selected derivatives does not erase these inputs.
!> \param state ...
!> \param domain FFT, exact atomic, interpolated background, PAW AE or compensated PAW PS.
!> \param points FFT indices, (global atom, atomic row, 1), or PAW (global atom, radial node, angular node).
!> \param density Physical model spin density, one paired RKS or two UKS columns.
!> \param kinetic Optional physical kinetic density with the same spin convention.
!>        Omission records that the functional did not request this descriptor.
!> \param append Append a further chunk in this domain; otherwise replace its previous samples.
! **************************************************************************************************
   SUBROUTINE skala_density_samples_put(state, domain, points, density, kinetic, append)
      TYPE(skala_density_constraints_type), &
         INTENT(INOUT)                                   :: state
      INTEGER, INTENT(IN)                                :: domain, points(:, :)
      REAL(dp), INTENT(IN)                               :: density(:, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: kinetic(:, :)
      LOGICAL, INTENT(IN), OPTIONAL                      :: append

      INTEGER                                            :: first, nk, np, ns
      INTEGER, ALLOCATABLE                               :: combined_points(:, :)
      LOGICAL                                            :: valid
      REAL(dp), ALLOCATABLE                              :: combined_density(:, :), &
                                                            combined_kinetic(:, :)

      IF (.NOT. state%sampling) RETURN
      valid = domain >= skala_fft_domain .AND. domain <= skala_last_sample_domain
      CPASSERT(valid)
      np = SIZE(density, 1)
      ns = SIZE(density, 2)
      CPASSERT(ns == 1 .OR. ns == 2)
      CPASSERT(ALL(SHAPE(points) == [3, np]))
      CPASSERT(ALL(ieee_is_finite(density)))
      nk = 0
      IF (PRESENT(kinetic)) THEN
         CPASSERT(ALL(SHAPE(kinetic) == [np, ns]))
         CPASSERT(ALL(ieee_is_finite(kinetic)))
         nk = ns
      END IF
      IF (state%sample_nspins /= 0) THEN
         CPASSERT(state%sample_nspins == ns)
      END IF
      state%sample_nspins = ns
      IF (domain == skala_paw_ae_domain .OR. domain == skala_paw_ps_domain) THEN
         CPASSERT(ALL(points >= 1))
      ELSE IF (domain /= skala_fft_domain) THEN
         CPASSERT(ALL(points(1:2, :) >= 1) .AND. ALL(points(3, :) == 1))
      END IF
      ASSOCIATE (sample => state%samples(domain))
         first = 0
         IF (PRESENT(append)) THEN
            IF (append .AND. ALLOCATED(sample%points)) THEN
               valid = SIZE(sample%density, 2) == ns .AND. SIZE(sample%kinetic, 2) == nk
               CPASSERT(valid)
               first = SIZE(sample%points, 2)
            END IF
         END IF
         ALLOCATE (combined_points(3, first + np), combined_density(first + np, ns), combined_kinetic(first + np, nk))
         IF (first > 0) THEN
            combined_points(:, :first) = sample%points
            combined_density(:first, :) = sample%density
            combined_kinetic(:first, :) = sample%kinetic
         END IF
         combined_points(:, first + 1:) = points
         combined_density(first + 1:, :) = density
         IF (PRESENT(kinetic)) combined_kinetic(first + 1:, :) = kinetic
         CALL MOVE_ALLOC(combined_points, sample%points)
         CALL MOVE_ALLOC(combined_density, sample%density)
         CALL MOVE_ALLOC(combined_kinetic, sample%kinetic)
      END ASSOCIATE
   END SUBROUTINE skala_density_samples_put

! **************************************************************************************************
!> \brief Copy one domain's native descriptor samples without gathering MPI-owned rows.
!> \param state ...
!> \param domain ...
!> \param points ...
!> \param density ...
!> \param kinetic Native kinetic inputs, or zero columns when this descriptor was not requested.
!> \param available Whether this domain was evaluated since the last begin call.
! **************************************************************************************************
   SUBROUTINE skala_density_samples_get(state, domain, points, density, kinetic, available)
      TYPE(skala_density_constraints_type), INTENT(IN)   :: state
      INTEGER, INTENT(IN)                                :: domain
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: points(:, :)
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: density(:, :), kinetic(:, :)
      LOGICAL, INTENT(OUT)                               :: available

      INTEGER                                            :: nk, np, ns
      LOGICAL                                            :: valid

      valid = domain >= skala_fft_domain .AND. domain <= skala_last_sample_domain
      CPASSERT(valid)
      ASSOCIATE (sample => state%samples(domain))
         available = ALLOCATED(sample%points)
         np = 0
         ns = state%sample_nspins
         nk = 0
         IF (available) THEN
            np = SIZE(sample%points, 2)
            nk = SIZE(sample%kinetic, 2)
         END IF
         ALLOCATE (points(3, np), density(np, ns), kinetic(np, nk))
         IF (available) THEN
            points(:, :) = sample%points
            density(:, :) = sample%density
            kinetic(:, :) = sample%kinetic
         END IF
      END ASSOCIATE
   END SUBROUTINE skala_density_samples_get

! **************************************************************************************************
!> \brief Release current samples, including domains absent from the next evaluation.
!> \param state ...
! **************************************************************************************************
   SUBROUTINE clear_samples(state)
      TYPE(skala_density_constraints_type), &
         INTENT(INOUT)                                   :: state

      INTEGER                                            :: domain

      state%sample_generation = state%sample_generation + 1
      DO domain = skala_fft_domain, skala_last_sample_domain
         ASSOCIATE (sample => state%samples(domain))
            IF (ALLOCATED(sample%points)) DEALLOCATE (sample%points, sample%density, sample%kinetic)
         END ASSOCIATE
      END DO
      state%sampling = .FALSE.
      state%sample_nspins = 0
      IF (state%collecting_jumps) THEN
         state%jump_values = 0.0_dp
         state%jump_owners = 0
      END IF
   END SUBROUTINE clear_samples

! **************************************************************************************************
!> \brief Detect a physical density target bracketed by two orbital evaluations.
!>        The caller supplies matching physical point/spin layouts, before nonlinear features.
!>        A field equal to the target at both ends supplies no sign-changing constraint.
!>        Endpoint equality is a candidate, not proof of an admissible stationary minimum.
!> \param first Physical spin densities at the first orbital state.
!> \param last Physical spin densities at the last orbital state.
!> \param crossings Candidate mask with the same physical point/spin layout.
!> \param TARGET ...
! **************************************************************************************************
   SUBROUTINE skala_density_crossings(first, last, crossings, TARGET)
      REAL(dp), INTENT(IN)                               :: first(:, :), last(:, :)
      LOGICAL, INTENT(OUT)                               :: crossings(:, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: target

      REAL(dp)                                           :: level

      level = 0.0_dp
      IF (PRESENT(TARGET)) level = TARGET
      CPASSERT(ieee_is_finite(level))
      CPASSERT(ALL(SHAPE(first) == SHAPE(last)))
      CPASSERT(ALL(SHAPE(first) == SHAPE(crossings)))
      CPASSERT(ALL(ieee_is_finite(first)))
      CPASSERT(ALL(ieee_is_finite(last)))
      crossings(:, :) = ((first < level .AND. last >= level) .OR. &
                         (first > level .AND. last <= level) .OR. &
                         (first == level .AND. last /= level))
   END SUBROUTINE skala_density_crossings

! **************************************************************************************************
!> \brief Locate a bracketed target by reevaluating the physical density along an orbital curve.
!>        The callback retains the actual normalization, density map and frozen core.
!>        Only the absolute physical density tolerance establishes convergence.
!>        Unbracketed and unresolved roots return the best observed endpoint with failure.
!> \param evaluate Callback for one fixed physical point/spin on the caller's orbital curve.
!> \param context Owned orbital-curve data and native evaluation workspace.
!> \param first_step First endpoint of the supplied orbital-curve interval.
!> \param last_step Last endpoint, strictly larger than the first.
!> \param tolerance Maximum absolute difference from the physical target.
!> \param step Best observed orbital-curve parameter.
!> \param value Its actual physical spin density.
!> \param converged Whether the physical density meets the requested tolerance.
!> \param max_iterations Maximum safeguarded interpolation steps, defaults to 100.
!> \param TARGET ...
! **************************************************************************************************
   SUBROUTINE skala_density_root(evaluate, context, first_step, last_step, tolerance, step, value, &
                                 converged, max_iterations, TARGET)
      PROCEDURE(evaluate_physical_density)                :: evaluate
      CLASS(*), INTENT(INOUT)                             :: context
      REAL(dp), INTENT(IN)                               :: first_step, last_step, tolerance
      REAL(dp), INTENT(OUT)                              :: step, value
      LOGICAL, INTENT(OUT)                               :: converged
      INTEGER, INTENT(IN), OPTIONAL                      :: max_iterations
      REAL(dp), INTENT(IN), OPTIONAL                     :: target

      INTEGER                                            :: iteration, limit
      LOGICAL                                            :: bracketed
      REAL(dp)                                           :: first, first_value, fraction, last, last_value, level, &
                                                            scale, trial, trial_value

      level = 0.0_dp
      IF (PRESENT(TARGET)) level = TARGET
      CPASSERT(ieee_is_finite(level))
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
      IF (ABS(last_value - level) < ABS(value - level)) THEN
         step = last
         value = last_value
      END IF
      converged = ABS(value - level) <= tolerance
      IF (converged) RETURN
      bracketed = (first_value < level .AND. last_value > level) .OR. &
                  (first_value > level .AND. last_value < level)
      IF (.NOT. bracketed) RETURN

      DO iteration = 1, limit
         scale = MAX(ABS(first_value - level), ABS(last_value - level))
         fraction = ((first_value - level)/scale)/((first_value - level)/scale - (last_value - level)/scale)
         fraction = MAX(0.1_dp, MIN(0.9_dp, fraction))
         trial = (1.0_dp - fraction)*first + fraction*last
         IF (trial <= first .OR. trial >= last) EXIT
         CALL evaluate(trial, trial_value, context)
         CPASSERT(ieee_is_finite(trial_value))
         IF (ABS(trial_value - level) < ABS(value - level)) THEN
            step = trial
            value = trial_value
         END IF
         converged = ABS(value - level) <= tolerance
         IF (converged) RETURN
         IF ((trial_value < level) .EQV. (first_value < level)) THEN
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

      CALL skala_density_constraints_clear_selection(state)
      CALL clear_samples(state)
   END SUBROUTINE skala_density_constraints_reset

! **************************************************************************************************
!> \brief Clear selected derivatives while retaining this trial's physical samples.
!> \param state ...
! **************************************************************************************************
   SUBROUTINE skala_density_constraints_clear_selection(state)
      TYPE(skala_density_constraints_type), &
         INTENT(INOUT)                                   :: state

      CALL clear_jumps(state)
      IF (ALLOCATED(state%points)) DEALLOCATE (state%points)
      IF (ALLOCATED(state%spins)) DEALLOCATE (state%spins)
      IF (ALLOCATED(state%descriptors)) DEALLOCATE (state%descriptors)
      IF (ALLOCATED(state%domains)) DEALLOCATE (state%domains)
      IF (ALLOCATED(state%atomic_sizes)) DEALLOCATE (state%atomic_sizes)
      IF (ALLOCATED(state%paw_sizes)) DEALLOCATE (state%paw_sizes)
      IF (ALLOCATED(state%positions)) DEALLOCATE (state%positions)
      IF (ALLOCATED(state%negative_weights)) DEALLOCATE (state%negative_weights)
      IF (ALLOCATED(state%targets)) DEALLOCATE (state%targets)
      state%nspins = 0
      state%npts = 0
      state%bounds = 0
      state%hmat = 0.0_dp
      state%probe = 0.0_dp
      state%tolerance = 0.0_dp
      state%kinetic_probe = 0.0_dp
      state%kinetic_tolerance = 0.0_dp
   END SUBROUTINE skala_density_constraints_clear_selection

! **************************************************************************************************
!> \brief Store physical grid constraints and their selected convex weights.
!>        Spin zero denotes the paired alpha/beta density of an RKS calculation.
!>        The orbital solver supplies the weights. This routine does not select them.
!> \param state ...
!> \param points Global FFT indices or (global atom, atomic row, 1), one column per constraint.
!> \param spins Zero for RKS, one or two for the physical UKS spin.
!> \param negative_weights Weight of the negative-side limiting derivative.
!> \param probe Positive physical density displacement used by the original model.
!> \param tolerance Maximum absolute spin density for applying its selected derivative.
!> \param nspins ...
!> \param npts ...
!> \param bounds ...
!> \param hmat ...
!> \param positions ...
!> \param descriptors Physical density or kinetic density, defaults to density for every constraint.
!> \param kinetic_probe Positive physical kinetic-density displacement, required for kinetic constraints.
!> \param kinetic_tolerance Maximum absolute spin kinetic density, required for kinetic constraints.
!> \param domains FFT, exact atomic or interpolated atomic-background evaluation, defaults to FFT.
!> \param atomic_sizes Complete row count per global atom, required for atomic constraints.
!> \param targets Fixed physical descriptor values, zero by default.
!> \param paw_sizes Original radial and angular counts per global atom, required for PAW constraints.
! **************************************************************************************************
   SUBROUTINE skala_density_constraints_set(state, points, spins, negative_weights, &
                                            probe, tolerance, nspins, npts, bounds, hmat, positions, &
                                           descriptors, kinetic_probe, kinetic_tolerance, domains, atomic_sizes, targets, paw_sizes)
      TYPE(skala_density_constraints_type), &
         INTENT(INOUT)                                   :: state
      INTEGER, INTENT(IN)                                :: points(:, :), spins(:)
      REAL(dp), INTENT(IN)                               :: negative_weights(:), probe, tolerance
      INTEGER, INTENT(IN)                                :: nspins, npts(3), bounds(2, 3)
      REAL(dp), INTENT(IN)                               :: hmat(3, 3), positions(:, :)
      INTEGER, INTENT(IN), OPTIONAL                      :: descriptors(:)
      REAL(dp), INTENT(IN), OPTIONAL                     :: kinetic_probe, kinetic_tolerance
      INTEGER, INTENT(IN), OPTIONAL                      :: domains(:), atomic_sizes(:)
      REAL(dp), INTENT(IN), OPTIONAL                     :: targets(:)
      INTEGER, INTENT(IN), OPTIONAL                      :: paw_sizes(:, :)

      INTEGER                                            :: i, j
      INTEGER, ALLOCATABLE                               :: selected_descriptors(:), &
                                                            selected_domains(:)
      LOGICAL                                            :: valid

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
      IF (PRESENT(targets)) THEN
         CPASSERT(SIZE(targets) == SIZE(spins))
         CPASSERT(ALL(ieee_is_finite(targets)))
      END IF
      ALLOCATE (selected_descriptors(SIZE(spins)))
      ALLOCATE (selected_domains(SIZE(spins)))
      selected_descriptors = skala_density_descriptor
      selected_domains = skala_fft_domain
      IF (PRESENT(domains)) THEN
         CPASSERT(SIZE(domains) == SIZE(spins))
         valid = ALL(domains >= skala_fft_domain .AND. domains <= skala_last_sample_domain)
         CPASSERT(valid)
         selected_domains(:) = domains
      END IF
      IF (ANY(selected_domains == skala_atom_domain .OR. selected_domains == skala_background_domain)) THEN
         CPASSERT(PRESENT(atomic_sizes))
         CPASSERT(SIZE(atomic_sizes) == SIZE(positions, 2))
         CPASSERT(ALL(atomic_sizes >= 0))
      END IF
      IF (ANY(selected_domains == skala_paw_ae_domain .OR. selected_domains == skala_paw_ps_domain)) THEN
         CPASSERT(PRESENT(paw_sizes))
         CPASSERT(ALL(SHAPE(paw_sizes) == [2, SIZE(positions, 2)]))
         CPASSERT(ALL(paw_sizes >= 0))
      END IF
      IF (PRESENT(descriptors)) THEN
         CPASSERT(SIZE(descriptors) == SIZE(spins))
         valid = ALL(descriptors == skala_density_descriptor .OR. descriptors == skala_kinetic_descriptor)
         CPASSERT(valid)
         selected_descriptors(:) = descriptors
      END IF
      IF (ANY(selected_descriptors == skala_kinetic_descriptor)) THEN
         CPASSERT(PRESENT(kinetic_probe) .AND. PRESENT(kinetic_tolerance))
      END IF
      IF (PRESENT(kinetic_probe)) THEN
         CPASSERT(ieee_is_finite(kinetic_probe) .AND. kinetic_probe > 0.0_dp)
      END IF
      IF (PRESENT(kinetic_tolerance)) THEN
         valid = ieee_is_finite(kinetic_tolerance) .AND. kinetic_tolerance >= 0.0_dp
         CPASSERT(valid)
      END IF
      IF (nspins == 1) THEN
         CPASSERT(ALL(spins == 0))
      ELSE
         CPASSERT(ALL(spins >= 1) .AND. ALL(spins <= 2))
      END IF
      DO i = 1, SIZE(spins)
         IF (selected_domains(i) == skala_fft_domain) THEN
            CPASSERT(ALL(points(:, i) >= bounds(1, :)))
            CPASSERT(ALL(points(:, i) <= bounds(2, :)))
         ELSE IF (selected_domains(i) == skala_atom_domain .OR. selected_domains(i) == skala_background_domain) THEN
            valid = points(1, i) >= 1 .AND. points(1, i) <= SIZE(atomic_sizes)
            CPASSERT(valid)
            valid = points(2, i) >= 1 .AND. points(2, i) <= atomic_sizes(points(1, i))
            CPASSERT(valid)
            CPASSERT(points(3, i) == 1)
         ELSE
            valid = points(1, i) >= 1 .AND. points(1, i) <= SIZE(paw_sizes, 2)
            CPASSERT(valid)
            CPASSERT(ALL(points(2:3, i) >= 1))
            CPASSERT(ALL(points(2:3, i) <= paw_sizes(:, points(1, i))))
         END IF
         DO j = 1, i - 1
            IF (ALL(points(:, i) == points(:, j)) .AND. spins(i) == spins(j) .AND. &
                selected_descriptors(i) == selected_descriptors(j) .AND. &
                selected_domains(i) == selected_domains(j)) THEN
               CPABORT("Duplicate SKALA physical descriptor constraint.")
            END IF
         END DO
      END DO

      CALL skala_density_constraints_clear_selection(state)
      state%nspins = nspins
      state%npts = npts
      state%bounds = bounds
      state%hmat = hmat
      state%probe = probe
      state%tolerance = tolerance
      IF (PRESENT(kinetic_probe)) state%kinetic_probe = kinetic_probe
      IF (PRESENT(kinetic_tolerance)) state%kinetic_tolerance = kinetic_tolerance
      ALLOCATE (state%points(3, SIZE(spins)), state%spins(SIZE(spins)), &
                state%descriptors(SIZE(spins)), state%domains(SIZE(spins)), state%negative_weights(SIZE(spins)), &
                state%positions(3, SIZE(positions, 2)), state%targets(SIZE(spins)))
      state%points(:, :) = points
      state%spins(:) = spins
      state%descriptors(:) = selected_descriptors
      state%domains(:) = selected_domains
      IF (PRESENT(atomic_sizes)) THEN
         ALLOCATE (state%atomic_sizes(SIZE(atomic_sizes)))
         state%atomic_sizes(:) = atomic_sizes
      END IF
      IF (PRESENT(paw_sizes)) THEN
         ALLOCATE (state%paw_sizes(2, SIZE(paw_sizes, 2)))
         state%paw_sizes(:, :) = paw_sizes
      END IF
      state%negative_weights(:) = negative_weights
      state%targets = 0.0_dp
      IF (PRESENT(targets)) state%targets(:) = targets
      state%positions(:, :) = positions
   END SUBROUTINE skala_density_constraints_set

! **************************************************************************************************
!> \brief Number of selected physical constraints.
!> \param state ...
!> \param domain Optional evaluation domain whose constraints are counted.
!> \return ...
! **************************************************************************************************
   FUNCTION skala_density_constraints_count(state, domain) RESULT(n)
      TYPE(skala_density_constraints_type), INTENT(IN)   :: state
      INTEGER, INTENT(IN), OPTIONAL                      :: domain
      INTEGER                                            :: n

      n = 0
      IF (ALLOCATED(state%spins)) n = SIZE(state%spins)
      IF (PRESENT(domain) .AND. n > 0) n = COUNT(state%domains == domain)
   END FUNCTION skala_density_constraints_count

! **************************************************************************************************
!> \brief Check that a selection belongs to this geometry and physical grid.
!> \param state ...
!> \param nspins ...
!> \param npts ...
!> \param bounds ...
!> \param hmat ...
!> \param positions ...
!> \param atomic_sizes Optional current atomic layout, also required to validate atomic selections.
!> \param paw_sizes Optional original radial and angular layout per global atom.
!> \return ...
! **************************************************************************************************
   FUNCTION skala_density_constraints_matches(state, nspins, npts, bounds, hmat, positions, atomic_sizes, paw_sizes) RESULT(matches)
      TYPE(skala_density_constraints_type), INTENT(IN)   :: state
      INTEGER, INTENT(IN)                                :: nspins, npts(3), bounds(2, 3)
      REAL(dp), INTENT(IN)                               :: hmat(3, 3), positions(:, :)
      INTEGER, INTENT(IN), OPTIONAL                      :: atomic_sizes(:), paw_sizes(:, :)
      LOGICAL                                            :: matches

      matches = .FALSE.
      IF (skala_density_constraints_count(state) == 0) RETURN
      IF (state%nspins /= nspins .OR. ANY(state%npts /= npts) .OR. ANY(state%bounds /= bounds)) RETURN
      IF (ANY(state%hmat /= hmat)) RETURN
      IF (ANY(SHAPE(state%positions) /= SHAPE(positions))) RETURN
      IF (PRESENT(atomic_sizes)) THEN
         IF (.NOT. ALLOCATED(state%atomic_sizes)) RETURN
         IF (SIZE(state%atomic_sizes) /= SIZE(atomic_sizes)) RETURN
         IF (ANY(state%atomic_sizes /= atomic_sizes)) RETURN
      END IF
      IF (PRESENT(paw_sizes)) THEN
         IF (.NOT. ALLOCATED(state%paw_sizes)) RETURN
         IF (ANY(SHAPE(state%paw_sizes) /= SHAPE(paw_sizes))) RETURN
         IF (ANY(state%paw_sizes /= paw_sizes)) RETURN
      END IF
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
!> \param descriptors Owned physical descriptor identifiers, required for kinetic constraints.
!> \param kinetic_probe ...
!> \param kinetic_tolerance ...
!> \param domains Owned evaluation-domain identifiers, required when retrieving several domains.
!> \param domain Optional domain filter, preserving the original constraint order.
!> \param targets Fixed physical target of each selected descriptor.
! **************************************************************************************************
   SUBROUTINE skala_density_constraints_get(state, points, spins, negative_weights, probe, tolerance, &
                                            descriptors, kinetic_probe, kinetic_tolerance, domains, domain, targets)
      TYPE(skala_density_constraints_type), INTENT(IN)   :: state
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: points(:, :), spins(:)
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: negative_weights(:)
      REAL(dp), INTENT(OUT)                              :: probe, tolerance
      INTEGER, ALLOCATABLE, INTENT(OUT), OPTIONAL        :: descriptors(:)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: kinetic_probe, kinetic_tolerance
      INTEGER, ALLOCATABLE, INTENT(OUT), OPTIONAL        :: domains(:)
      INTEGER, INTENT(IN), OPTIONAL                      :: domain
      REAL(dp), ALLOCATABLE, INTENT(OUT), OPTIONAL       :: targets(:)

      INTEGER                                            :: i, n
      INTEGER, ALLOCATABLE                               :: selected(:)

      n = skala_density_constraints_count(state, domain)
      CPASSERT(n > 0)
      ALLOCATE (selected(n))
      IF (PRESENT(domain)) THEN
         selected(:) = PACK([(i, i=1, SIZE(state%spins))], state%domains == domain)
      ELSE
         selected(:) = [(i, i=1, n)]
         IF (ANY(state%domains /= skala_fft_domain)) THEN
            CPASSERT(PRESENT(domains))
         END IF
      END IF
      IF (ANY(state%descriptors(selected) == skala_kinetic_descriptor)) THEN
         CPASSERT(PRESENT(descriptors))
         CPASSERT(PRESENT(kinetic_probe) .AND. PRESENT(kinetic_tolerance))
      END IF
      ALLOCATE (points(3, n), spins(n), negative_weights(n))
      points(:, :) = state%points(:, selected)
      spins(:) = state%spins(selected)
      negative_weights(:) = state%negative_weights(selected)
      IF (PRESENT(targets)) THEN
         ALLOCATE (targets(n))
         targets(:) = state%targets(selected)
      ELSE IF (ANY(state%targets(selected) /= 0.0_dp)) THEN
         CPABORT("Nonzero physical descriptor constraints require explicit targets.")
      END IF
      probe = state%probe
      tolerance = state%tolerance
      IF (PRESENT(descriptors)) THEN
         ALLOCATE (descriptors(n))
         descriptors(:) = state%descriptors(selected)
      END IF
      IF (PRESENT(domains)) THEN
         ALLOCATE (domains(n))
         domains(:) = state%domains(selected)
      END IF
      IF (PRESENT(kinetic_probe)) kinetic_probe = state%kinetic_probe
      IF (PRESENT(kinetic_tolerance)) kinetic_tolerance = state%kinetic_tolerance
   END SUBROUTINE skala_density_constraints_get

! **************************************************************************************************
!> \brief Replace convex weights without changing the selected physical samples or geometry.
!>        A native caller must rebuild its KS operator after changing these weights.
!> \param state Registered physical constraints.
!> \param negative_weights Negative-side weights in the existing constraint order.
! **************************************************************************************************
   SUBROUTINE skala_density_constraints_update_weights(state, negative_weights)
      TYPE(skala_density_constraints_type), &
         INTENT(INOUT)                                   :: state
      REAL(dp), INTENT(IN)                               :: negative_weights(:)

      CPASSERT(SIZE(negative_weights) == skala_density_constraints_count(state))
      CPASSERT(ALL(ieee_is_finite(negative_weights)))
      CPASSERT(ALL(negative_weights >= 0.0_dp .AND. negative_weights <= 1.0_dp))
      IF (SIZE(negative_weights) == 0) RETURN
      CPASSERT(ALLOCATED(state%negative_weights))
      state%negative_weights(:) = negative_weights
   END SUBROUTINE skala_density_constraints_update_weights

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
