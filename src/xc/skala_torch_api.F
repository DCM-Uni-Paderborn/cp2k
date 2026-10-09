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
   USE ieee_arithmetic, ONLY: ieee_is_finite
#if defined (__HAS_IEEE_EXCEPTIONS)
   USE ieee_exceptions, ONLY: ieee_all, &
                              ieee_get_halting_mode, &
                              ieee_set_halting_mode
#endif
   USE kinds, ONLY: default_string_length, &
                    dp
   USE string_utilities, ONLY: uppercase
   USE torch_api, ONLY: &
      torch_dict_type, torch_model_disable_parameter_gradients, torch_model_forward_mol_tensor, &
      torch_model_freeze_preserving_method, torch_model_load_with_metadata, &
      torch_model_release, torch_model_remap_device_constants, torch_model_promote_float32, &
      torch_model_prepare_higher_derivatives, &
      torch_model_type, &
      torch_dict_clone, torch_dict_get, torch_dict_insert, torch_dict_release, &
      torch_tensor_accumulate, torch_tensor_autograd, torch_tensor_cat, torch_tensor_from_array, &
      torch_tensor_hessian_vector, torch_tensor_linear_combination, torch_tensor_multiply, torch_tensor_narrow, &
      torch_tensor_to_device_leaf, torch_tensor_to_host, &
      torch_tensor_item_double, torch_tensor_release, torch_tensor_type, &
      torch_tensor_weighted_sum
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'skala_torch_api'

   PUBLIC :: skala_torch_model_type, skala_torch_model_load, skala_torch_model_release
   PUBLIC :: skala_torch_model_get_exc, skala_torch_model_get_exc_density
   PUBLIC :: skala_torch_model_spin_projected_eval
   PUBLIC :: skala_torch_model_needs_feature, skala_torch_model_protocol_version

   TYPE skala_torch_model_type
      PRIVATE
      INTEGER                                            :: protocol_version = -1
      CHARACTER(len=default_string_length), ALLOCATABLE, &
         DIMENSION(:)                                    :: features
      TYPE(torch_model_type)                             :: torch_model
   END TYPE skala_torch_model_type

CONTAINS

! **************************************************************************************************
!> \brief Load a SKALA TorchScript model and its feature metadata.
!> \param model ...
!> \param filename ...
!> \param fp64 Promote float32 numerical state and casts; default true for differential consistency.
!> \param higher_derivatives Use FP64 model arithmetic and stable higher input derivatives.
! **************************************************************************************************
   SUBROUTINE skala_torch_model_load(model, filename, fp64, higher_derivatives)
      TYPE(skala_torch_model_type), INTENT(INOUT)        :: model
      CHARACTER(len=*), INTENT(IN)                       :: filename
      LOGICAL, INTENT(IN), OPTIONAL                      :: fp64, higher_derivatives

      CHARACTER(:), ALLOCATABLE                          :: features_json, protocol_string
      INTEGER                                            :: ios
      LOGICAL                                            :: prepare_derivatives, use_fp64

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
      prepare_derivatives = .FALSE.
      IF (PRESENT(higher_derivatives)) prepare_derivatives = higher_derivatives
      IF (prepare_derivatives) THEN
         CALL torch_model_prepare_higher_derivatives(model%torch_model)
      ELSE IF (use_fp64) THEN
         CALL torch_model_promote_float32(model%torch_model)
      END IF
      ! Preserve the exported SKALA entry point while folding constant model state.
      CALL torch_model_freeze_preserving_method(model%torch_model, "get_exc_density")

   END SUBROUTINE skala_torch_model_load

! **************************************************************************************************
!> \brief Release a loaded SKALA TorchScript model.
!> \param model ...
! **************************************************************************************************
   SUBROUTINE skala_torch_model_release(model)
      TYPE(skala_torch_model_type), INTENT(INOUT)        :: model

      CALL torch_model_release(model%torch_model)
      IF (ALLOCATED(model%features)) DEALLOCATE (model%features)
      model%protocol_version = -1

   END SUBROUTINE skala_torch_model_release

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
!> \brief Streaming multicollinear extension of the unchanged collinear Skala model.
!> \param model Loaded Skala model.
!> \param inputs Molecular features with scalar/z/x/y primitive fields.
!> \param axes Unit angular directions in z/x/y order.
!> \param weights Angular weights summing to one.
!> \param variables Independent primitive and optional geometry tensors.
!> \param outputs Complete-energy covectors or their Hessian action.
!> \param exc Base-state energy.
!> \param status Zero on success. The caller must check failures collectively.
!> \param directions Optional joint input direction.
! **************************************************************************************************
   SUBROUTINE skala_torch_model_spin_projected_eval(model, inputs, axes, weights, variables, outputs, exc, status, directions)
      TYPE(skala_torch_model_type), INTENT(INOUT)        :: model
      TYPE(torch_dict_type), INTENT(IN)                  :: inputs
      REAL(KIND=dp), INTENT(IN)                          :: axes(:, :), weights(:)
      TYPE(torch_tensor_type), INTENT(IN)                :: variables(:)
      TYPE(torch_tensor_type), INTENT(INOUT)             :: outputs(:)
      REAL(KIND=dp), INTENT(OUT)                         :: exc
      INTEGER, INTENT(OUT)                               :: status
      TYPE(torch_tensor_type), INTENT(IN), OPTIONAL      :: directions(:)

      CHARACTER(len=7), PARAMETER :: keys(3) = [CHARACTER(len=7) :: "density", "grad", "kin"]

      INTEGER                                            :: f, i, j
      REAL(KIND=dp), ALLOCATABLE                         :: lambda_value(:)
      TYPE(torch_dict_type)                              :: projected
      TYPE(torch_tensor_type) :: channels(2), collinear, components(4), connected, effective, &
         energy, fields(3), grid_weights, lambda, scaled, scaling(1), spin, total, values, &
         derivatives(SIZE(variables))

      CPASSERT(SIZE(axes, 1) == 3)
      CPASSERT(SIZE(axes, 2) == SIZE(weights) .AND. SIZE(weights) > 0)
      CPASSERT(ALL(ieee_is_finite(axes)) .AND. ALL(ieee_is_finite(weights)))
      CPASSERT(ALL(weights >= 0.0_dp) .AND. ABS(SUM(weights) - 1.0_dp) < 1.0e-12_dp)
      CPASSERT(ALL(ABS(SUM(axes**2, DIM=1) - 1.0_dp) < 1.0e-12_dp))
      CPASSERT(SIZE(variables) > 0 .AND. SIZE(outputs) == SIZE(variables))
      IF (PRESENT(directions)) THEN
         CPASSERT(SIZE(directions) == SIZE(variables))
      END IF
      status = 0
      exc = 0.0_dp
      DO f = 1, SIZE(keys)
         CALL torch_dict_get(inputs, TRIM(keys(f)), fields(f), to_host=.FALSE.)
      END DO
      CALL torch_dict_get(inputs, "grid_weights", grid_weights, to_host=.FALSE.)
      ALLOCATE (lambda_value(1))
      lambda_value = 1.0_dp
      CALL torch_tensor_from_array(lambda, lambda_value)
      CALL torch_tensor_to_device_leaf(lambda, .TRUE.)

      ! Multicollinear extension averages E_col + dE_col/dlambda over spin directions.
      ! The common lambda scales all three spin fields (Phys. Rev. Research 5, 013036).
      DO i = 1, SIZE(weights)
         CALL torch_dict_clone(inputs, projected)
         DO f = 1, SIZE(keys)
            DO j = 1, 4
               CALL torch_tensor_narrow(fields(f), 0, j - 1, 1, components(j))
            END DO
            CALL torch_tensor_linear_combination(components(2:4), axes(:, i), spin)
            CALL torch_tensor_multiply(lambda, spin, scaled)
            CALL torch_tensor_linear_combination([components(1), scaled], [0.5_dp, 0.5_dp], channels(1))
            CALL torch_tensor_linear_combination([components(1), scaled], [0.5_dp, -0.5_dp], channels(2))
            CALL torch_tensor_cat(channels, collinear)
            CALL torch_dict_insert(projected, TRIM(keys(f)), collinear)
            DO j = 1, 4
               CALL torch_tensor_release(components(j))
            END DO
            DO j = 1, 2
               CALL torch_tensor_release(channels(j))
            END DO
            CALL torch_tensor_release(spin)
            CALL torch_tensor_release(scaled)
            CALL torch_tensor_release(collinear)
         END DO
         CALL torch_model_forward_mol_tensor(model%torch_model, "get_exc_density", projected, values, status)
         IF (status /= 0) EXIT
         CALL torch_tensor_weighted_sum(values, grid_weights, energy)
         CALL torch_tensor_linear_combination([energy, lambda], [1.0_dp, 0.0_dp], connected)
         CALL torch_tensor_autograd(connected, [lambda], scaling, .TRUE., status)
         IF (status /= 0) EXIT
         CALL torch_tensor_linear_combination([connected, scaling(1)], [1.0_dp, 1.0_dp], effective)
         CALL torch_tensor_accumulate(total, effective, weights(i))
         IF (PRESENT(directions)) THEN
            CALL torch_tensor_hessian_vector(effective, variables, directions, derivatives, status, to_host=.FALSE.)
         ELSE
            CALL torch_tensor_autograd(effective, variables, derivatives, .FALSE., status)
         END IF
         IF (status /= 0) EXIT
         DO j = 1, SIZE(variables)
            CALL torch_tensor_accumulate(outputs(j), derivatives(j), weights(i))
            CALL torch_tensor_release(derivatives(j))
         END DO
         CALL torch_tensor_release(values)
         CALL torch_tensor_release(energy)
         CALL torch_tensor_release(connected)
         CALL torch_tensor_release(scaling(1))
         CALL torch_tensor_release(effective)
         CALL torch_dict_release(projected)
      END DO
      IF (status == 0) THEN
         CALL torch_tensor_to_host(total, status)
         IF (status == 0) exc = torch_tensor_item_double(total)
         DO j = 1, SIZE(outputs)
            IF (status /= 0) EXIT
            CALL torch_tensor_to_host(outputs(j), status)
         END DO
      END IF
      IF (status /= 0) THEN
         exc = 0.0_dp
         DO j = 1, SIZE(outputs)
            CALL torch_tensor_release(outputs(j))
         END DO
      END IF
      DO j = 1, SIZE(derivatives)
         CALL torch_tensor_release(derivatives(j))
      END DO
      DO f = 1, SIZE(fields)
         CALL torch_tensor_release(fields(f))
      END DO
      CALL torch_tensor_release(values)
      CALL torch_tensor_release(energy)
      CALL torch_tensor_release(connected)
      CALL torch_tensor_release(scaling(1))
      CALL torch_tensor_release(effective)
      CALL torch_tensor_release(total)
      CALL torch_tensor_release(lambda)
      CALL torch_tensor_release(grid_weights)
      CALL torch_dict_release(projected)
   END SUBROUTINE skala_torch_model_spin_projected_eval

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
