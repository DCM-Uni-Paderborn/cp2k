!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Small CP2K wrapper around the SKALA TorchScript functional protocol.
! **************************************************************************************************
MODULE skala_torch_api
#if defined (__HAS_IEEE_EXCEPTIONS)
   USE ieee_exceptions, ONLY: ieee_all, &
                              ieee_get_halting_mode, &
                              ieee_set_halting_mode
#endif
   USE kinds, ONLY: default_string_length, &
                    dp
   USE string_utilities, ONLY: uppercase
   USE torch_api, ONLY: &
      torch_dict_clone_detached, torch_dict_get, torch_dict_insert, torch_dict_release, torch_dict_type, &
      torch_model_disable_parameter_gradients, torch_model_feature_clamps, torch_model_forward_mol_tensor, &
      torch_model_freeze_preserving_method, torch_model_load_with_metadata, torch_model_release, &
      torch_model_normalize_gradient_squares, &
      torch_model_remap_device_constants, torch_model_promote_float32, torch_model_type, &
      torch_tensor_backward_scalar, torch_tensor_data_ptr, torch_tensor_from_array, torch_tensor_grad, &
      torch_tensor_item_double, torch_tensor_release, torch_tensor_to_device_leaf, torch_tensor_type, &
      torch_tensor_weighted_sum
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'skala_torch_api'

   PUBLIC :: skala_torch_model_type, skala_torch_model_load, skala_torch_model_release
   PUBLIC :: skala_torch_model_get_exc, skala_torch_model_get_exc_density
   PUBLIC :: skala_torch_model_density_limits, skala_torch_model_descriptor_limits
   PUBLIC :: skala_torch_model_needs_feature, skala_torch_model_protocol_version
   PUBLIC :: skala_torch_model_density_floors

   TYPE skala_torch_model_type
      PRIVATE
      INTEGER                                            :: protocol_version = -1
      CHARACTER(len=default_string_length), ALLOCATABLE, &
         DIMENSION(:)                                    :: features
      REAL(dp), ALLOCATABLE                              :: density_floors(:)
      TYPE(torch_model_type)                             :: torch_model
   END TYPE skala_torch_model_type

CONTAINS

! **************************************************************************************************
!> \brief Load a SKALA TorchScript model and its feature metadata.
!> \param model ...
!> \param filename ...
!> \param fp64 Promote float32 numerical state and casts; default true for differential consistency.
! **************************************************************************************************
   SUBROUTINE skala_torch_model_load(model, filename, fp64)
      TYPE(skala_torch_model_type), INTENT(INOUT)        :: model
      CHARACTER(len=*), INTENT(IN)                       :: filename
      LOGICAL, INTENT(IN), OPTIONAL                      :: fp64

      CHARACTER(:), ALLOCATABLE                          :: features_json, protocol_string
      INTEGER                                            :: ios
      LOGICAL                                            :: use_fp64

      CALL torch_model_load_with_metadata(model%torch_model, filename, &
                                          "protocol_version", protocol_string, &
                                          "features", features_json)
      CALL torch_model_remap_device_constants(model%torch_model)
      CALL torch_model_disable_parameter_gradients(model%torch_model)
      READ (protocol_string, *, IOSTAT=ios) model%protocol_version
      IF (ios /= 0) CPABORT("Could not parse SKALA TorchScript protocol_version metadata")
      IF (model%protocol_version /= 2) THEN
         CPABORT("Unsupported SKALA TorchScript protocol version")
      END IF

      CALL parse_feature_list(features_json, model%features)
      use_fp64 = .TRUE.
      IF (PRESENT(fp64)) use_fp64 = fp64
      IF (use_fp64) CALL torch_model_promote_float32(model%torch_model)
      ! Preserve the exported SKALA entry point while folding constant model state.
      CALL torch_model_freeze_preserving_method(model%torch_model, "get_exc_density")
      ! abs(|grad n|^2) is smooth in grad n. Its unreduced graph gives a false Hessian at zero.
      CALL torch_model_normalize_gradient_squares(model%torch_model, "get_exc_density")
      CALL torch_model_feature_clamps(model%torch_model, "get_exc_density", "density", model%density_floors)

   END SUBROUTINE skala_torch_model_load

! **************************************************************************************************
!> \brief Release a loaded SKALA TorchScript model.
!> \param model ...
! **************************************************************************************************
   SUBROUTINE skala_torch_model_release(model)
      TYPE(skala_torch_model_type), INTENT(INOUT)        :: model

      CALL torch_model_release(model%torch_model)
      IF (ALLOCATED(model%features)) DEALLOCATE (model%features)
      IF (ALLOCATED(model%density_floors)) DEALLOCATE (model%density_floors)
      model%protocol_version = -1

   END SUBROUTINE skala_torch_model_release

! **************************************************************************************************
!> \brief Return the loaded model's constant clamp bounds for each physical spin density.
!> \param model Loaded model.
!> \return Bounds in atomic density units, excluding clamps on transformed features or indices.
! **************************************************************************************************
   FUNCTION skala_torch_model_density_floors(model) RESULT(floors)
      TYPE(skala_torch_model_type), INTENT(IN)           :: model
      REAL(dp), ALLOCATABLE                              :: floors(:)

      CPASSERT(ALLOCATED(model%density_floors))
      ALLOCATE (floors(SIZE(model%density_floors)))
      floors(:) = model%density_floors
   END FUNCTION skala_torch_model_density_floors

! **************************************************************************************************
!> \brief Check whether a loaded SKALA model requests a feature.
!> \param model ...
!> \param feature ...
!> \return ...
! **************************************************************************************************
   FUNCTION skala_torch_model_needs_feature(model, feature) RESULT(needs_feature)
      TYPE(skala_torch_model_type), INTENT(IN)           :: model
      CHARACTER(len=*), INTENT(IN)                       :: feature
      LOGICAL                                            :: needs_feature

      CHARACTER(len=default_string_length)               :: feature_key, model_feature
      INTEGER                                            :: i

      feature_key = ADJUSTL(feature)
      CALL uppercase(feature_key)

      needs_feature = .FALSE.
      IF (.NOT. ALLOCATED(model%features)) RETURN

      DO i = 1, SIZE(model%features)
         model_feature = ADJUSTL(model%features(i))
         CALL uppercase(model_feature)
         IF (TRIM(model_feature) == TRIM(feature_key)) THEN
            needs_feature = .TRUE.
            RETURN
         END IF
      END DO

   END FUNCTION skala_torch_model_needs_feature

! **************************************************************************************************
!> \brief Return the loaded SKALA TorchScript protocol version.
!> \param model ...
!> \return ...
! **************************************************************************************************
   FUNCTION skala_torch_model_protocol_version(model) RESULT(protocol_version)
      TYPE(skala_torch_model_type), INTENT(IN)           :: model
      INTEGER                                            :: protocol_version

      protocol_version = model%protocol_version

   END FUNCTION skala_torch_model_protocol_version

! **************************************************************************************************
!> \brief Evaluate the SKALA exchange-correlation energy density.
!> \param model ...
!> \param inputs ...
!> \param exc_density ...
! **************************************************************************************************
   SUBROUTINE skala_torch_model_get_exc_density(model, inputs, exc_density)
      TYPE(skala_torch_model_type), INTENT(INOUT)        :: model
      TYPE(torch_dict_type), INTENT(IN)                  :: inputs
      TYPE(torch_tensor_type), INTENT(INOUT)             :: exc_density

#if defined (__HAS_IEEE_EXCEPTIONS)
      LOGICAL, DIMENSION(5)                              :: ieee_halt

      CALL ieee_get_halting_mode(IEEE_ALL, ieee_halt)
      CALL ieee_set_halting_mode(IEEE_ALL, .FALSE.)
#endif
      CALL torch_model_forward_mol_tensor(model%torch_model, "get_exc_density", inputs, exc_density)
#if defined (__HAS_IEEE_EXCEPTIONS)
      CALL ieee_set_halting_mode(IEEE_ALL, ieee_halt)
#endif

   END SUBROUTINE skala_torch_model_get_exc_density

! **************************************************************************************************
!> \brief Evaluate the weighted SKALA exchange-correlation energy.
!> \param model ...
!> \param inputs ...
!> \param grid_weights ...
!> \param exc_tensor ...
!> \param exc ...
! **************************************************************************************************
   SUBROUTINE skala_torch_model_get_exc(model, inputs, grid_weights, exc_tensor, exc)
      TYPE(skala_torch_model_type), INTENT(INOUT)        :: model
      TYPE(torch_dict_type), INTENT(IN)                  :: inputs
      TYPE(torch_tensor_type), INTENT(IN)                :: grid_weights
      TYPE(torch_tensor_type), INTENT(INOUT)             :: exc_tensor
      REAL(KIND=dp), INTENT(OUT)                         :: exc

      TYPE(torch_tensor_type)                            :: exc_density

#if defined (__HAS_IEEE_EXCEPTIONS)
      LOGICAL, DIMENSION(5)                              :: ieee_halt

      CALL ieee_get_halting_mode(IEEE_ALL, ieee_halt)
      CALL ieee_set_halting_mode(IEEE_ALL, .FALSE.)
#endif
      CALL torch_model_forward_mol_tensor(model%torch_model, "get_exc_density", inputs, exc_density)
      CALL torch_tensor_weighted_sum(exc_density, grid_weights, exc_tensor)
      CALL torch_tensor_release(exc_density)
      exc = torch_tensor_item_double(exc_tensor)
#if defined (__HAS_IEEE_EXCEPTIONS)
      CALL ieee_set_halting_mode(IEEE_ALL, ieee_halt)
#endif

   END SUBROUTINE skala_torch_model_get_exc

! **************************************************************************************************
!> \brief Evaluate both density covectors at selected zero-density samples.
!>        The original model and input tensors are unchanged. Both calculations have detached
!>        static inputs and a fresh density leaf, so they do not accumulate original input gradients.
!>        The caller must satisfy the selected density constraints before using their convex weights.
!> \param model Original loaded model.
!> \param inputs Complete original feature dictionary, including its integration weights.
!> \param density Original full two-spin density layout, with one row per model feature.
!> \param active One-based feature row and spin indices, with shape (2, number of selected samples).
!> \param probe Positive density magnitude for the one-sided limit. Refine it to verify the limit.
!> \param negative New owned weighted density covector on the negative side.
!> \param positive New owned weighted density covector on the positive side.
!> \param energies Complete discrete energies of the two limit evaluations.
!> \param targets Optional physical target per selected feature entry, zero by default.
!> \param adjacent Use adjacent double-precision branch values, false by default.
! **************************************************************************************************
   SUBROUTINE skala_torch_model_density_limits(model, inputs, density, active, probe, negative, positive, energies, &
                                               targets, adjacent)
      TYPE(skala_torch_model_type), INTENT(INOUT)        :: model
      TYPE(torch_dict_type), INTENT(IN)                  :: inputs
      REAL(dp), INTENT(IN)                               :: density(:, :)
      INTEGER, INTENT(IN)                                :: active(:, :)
      REAL(dp), INTENT(IN)                               :: probe
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: negative(:, :), positive(:, :)
      REAL(dp), INTENT(OUT)                              :: energies(2)
      REAL(dp), INTENT(IN), OPTIONAL                     :: targets(:)
      LOGICAL, INTENT(IN), OPTIONAL                      :: adjacent

      CALL skala_torch_model_descriptor_limits(model, inputs, "density", density, active, probe, &
                                               negative, positive, energies, targets=targets, adjacent=adjacent)
   END SUBROUTINE skala_torch_model_density_limits

! **************************************************************************************************
!> \brief Evaluate original-model one-sided covectors of a physical density or kinetic density.
!>        Other features are detached and unchanged. Original input gradients are not accumulated.
!>        Literal probes allow independent convergence checks. Adjacent double-precision values
!>        isolate nearby targets without crossing a second boundary. Zero uses the smallest normal
!>        signed values to preserve its side through device transfers. This mode requires model
!>        arithmetic that retains double-precision features.
!> \param model Original loaded model.
!> \param inputs Complete original feature dictionary, including its integration weights.
!> \param feature Dictionary key, density or kin.
!> \param values Original full two-spin physical descriptor layout.
!> \param active One-based feature row and spin indices, shape (2, number of selected samples).
!> \param probe Positive physical descriptor magnitude for literal probes.
!> \param negative New owned weighted descriptor covector on the negative side.
!> \param positive New owned weighted descriptor covector on the positive side.
!> \param energies Optional complete discrete energies. Requesting them evaluates both branches.
!> \param branches Optional negative/positive branch selection. Unrequested covectors are zero.
!> \param targets Optional physical target per selected feature entry, zero by default.
!> \param adjacent Use adjacent branch values instead of literal probes, false by default.
! **************************************************************************************************
   SUBROUTINE skala_torch_model_descriptor_limits(model, inputs, feature, values, active, probe, &
                                                  negative, positive, energies, branches, targets, adjacent)
      TYPE(skala_torch_model_type), INTENT(INOUT)        :: model
      TYPE(torch_dict_type), INTENT(IN)                  :: inputs
      CHARACTER(len=*), INTENT(IN)                       :: feature
      REAL(dp), INTENT(IN)                               :: values(:, :)
      INTEGER, INTENT(IN)                                :: active(:, :)
      REAL(dp), INTENT(IN)                               :: probe
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: negative(:, :), positive(:, :)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: energies(2)
      LOGICAL, INTENT(IN), OPTIONAL                      :: branches(2)
      REAL(dp), INTENT(IN), OPTIONAL                     :: targets(:)
      LOGICAL, INTENT(IN), OPTIONAL                      :: adjacent

      INTEGER                                            :: branch, sample
      LOGICAL                                            :: local_limit, needed(2)
      REAL(dp)                                           :: energy, TARGET
      REAL(dp), ALLOCATABLE                              :: branch_values(:, :)
      REAL(dp), POINTER                                  :: derivative(:, :)
      TYPE(torch_dict_type)                              :: branch_inputs
      TYPE(torch_tensor_type)                            :: derivative_t, descriptor_t, energy_t, &
                                                            weights_t

      CPASSERT(feature == "density" .OR. feature == "kin")
      CPASSERT(SIZE(values, 2) == 2)
      CPASSERT(SIZE(active, 1) == 2 .AND. SIZE(active, 2) > 0)
      CPASSERT(ALL(active(1, :) >= 1) .AND. ALL(active(1, :) <= SIZE(values, 1)))
      CPASSERT(ALL(active(2, :) >= 1) .AND. ALL(active(2, :) <= 2))
      CPASSERT(probe > 0.0_dp .AND. probe < HUGE(probe))
      IF (PRESENT(targets)) THEN
         CPASSERT(SIZE(targets) == SIZE(active, 2))
         CPASSERT(ALL(ABS(targets) < HUGE(probe)))
      END IF
      ALLOCATE (branch_values(SIZE(values, 1), 2))
      ALLOCATE (negative(SIZE(values, 1), 2), positive(SIZE(values, 1), 2))
      negative = 0.0_dp
      positive = 0.0_dp
      needed = .TRUE.
      local_limit = .FALSE.
      IF (PRESENT(adjacent)) local_limit = adjacent
      IF (PRESENT(branches)) needed = branches
      IF (PRESENT(energies)) needed = .TRUE.
      CPASSERT(ANY(needed))
      DO branch = 1, 2
         IF (.NOT. needed(branch)) CYCLE
         branch_values(:, :) = values
         DO sample = 1, SIZE(active, 2)
            TARGET = 0.0_dp
            IF (PRESENT(targets)) TARGET = targets(sample)
            IF (local_limit) THEN
               branch_values(active(1, sample), active(2, sample)) = descriptor_adjacent(TARGET, 2*branch - 3)
            ELSE
               branch_values(active(1, sample), active(2, sample)) = TARGET + REAL(2*branch - 3, dp)*probe
            END IF
         END DO
         CALL torch_dict_clone_detached(inputs, branch_inputs, feature)
         CALL torch_tensor_from_array(descriptor_t, branch_values, requires_grad=.TRUE.)
         CALL torch_tensor_to_device_leaf(descriptor_t, .TRUE.)
         CALL torch_dict_insert(branch_inputs, feature, descriptor_t)
         CALL torch_dict_get(branch_inputs, "grid_weights", weights_t)
         CALL skala_torch_model_get_exc(model, branch_inputs, weights_t, energy_t, energy)
         IF (PRESENT(energies)) energies(branch) = energy
         CALL torch_tensor_backward_scalar(energy_t)
         CALL torch_tensor_grad(descriptor_t, derivative_t)
         NULLIFY (derivative)
         CALL torch_tensor_data_ptr(derivative_t, derivative)
         CPASSERT(ALL(SHAPE(derivative) == SHAPE(values)))
         IF (branch == 1) THEN
            negative(:, :) = derivative
         ELSE
            positive(:, :) = derivative
         END IF
         CALL torch_tensor_release(derivative_t)
         CALL torch_tensor_release(energy_t)
         CALL torch_tensor_release(weights_t)
         CALL torch_dict_release(branch_inputs)
         CALL torch_tensor_release(descriptor_t)
      END DO
      DEALLOCATE (branch_values)
   END SUBROUTINE skala_torch_model_descriptor_limits

! **************************************************************************************************
!> \brief Select a representable side of a double-precision descriptor boundary.
!> \param center Physical boundary value.
!> \param side Negative or positive side, minus or plus one.
!> \return Adjacent value, or the smallest signed normal value at zero.
! **************************************************************************************************
   FUNCTION descriptor_adjacent(center, side) RESULT(value)
      REAL(dp), INTENT(IN)                               :: center
      INTEGER, INTENT(IN)                                :: side
      REAL(dp)                                           :: value

      CPASSERT(ABS(side) == 1)
      IF (center == 0.0_dp) THEN
         value = REAL(side, dp)*TINY(0.0_dp)
      ELSE
         value = NEAREST(center, REAL(side, dp))
      END IF
      CPASSERT(ABS(value) < HUGE(0.0_dp))
      IF (side < 0) THEN
         CPASSERT(value < center)
      ELSE
         CPASSERT(value > center)
      END IF
   END FUNCTION descriptor_adjacent

! **************************************************************************************************
!> \brief Parse a TorchScript extra_files JSON list of feature names.
!> \param features_json ...
!> \param features ...
! **************************************************************************************************
   SUBROUTINE parse_feature_list(features_json, features)
      CHARACTER(len=*), INTENT(IN)                       :: features_json
      CHARACTER(len=default_string_length), &
         ALLOCATABLE, DIMENSION(:), INTENT(OUT)          :: features

      INTEGER                                            :: end_pos, feature_count, i, pos, quote1, &
                                                            quote2, start_pos

      feature_count = 0
      pos = 1
      DO
         quote1 = INDEX(features_json(pos:), '"')
         IF (quote1 == 0) EXIT
         start_pos = pos + quote1
         quote2 = INDEX(features_json(start_pos:), '"')
         IF (quote2 == 0) EXIT
         feature_count = feature_count + 1
         pos = start_pos + quote2
      END DO

      IF (feature_count == 0) CPABORT("SKALA TorchScript model does not list any features")
      ALLOCATE (features(feature_count))
      features = ""

      pos = 1
      DO i = 1, feature_count
         quote1 = INDEX(features_json(pos:), '"')
         start_pos = pos + quote1
         quote2 = INDEX(features_json(start_pos:), '"')
         end_pos = start_pos + quote2 - 2
         features(i) = features_json(start_pos:end_pos)
         pos = start_pos + quote2
      END DO

   END SUBROUTINE parse_feature_list

END MODULE skala_torch_api
