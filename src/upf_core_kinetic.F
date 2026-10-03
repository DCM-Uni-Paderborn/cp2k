!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Explicitly normalized UPF model-core kinetic fields and their geometry adjoints.
!>        Retain the original radial samples, signs and quadrature.
! **************************************************************************************************
MODULE upf_core_kinetic
   USE atom_upf,                        ONLY: atom_upfpot_type,&
                                              upf_nlcc_kinetic,&
                                              upf_radial_weights
   USE bessel_lib,                      ONLY: bessel0
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: fourpi
   USE uspp_augmentation,               ONLY: upf_radial_scalar_init,&
                                              uspp_augmentation_realspace,&
                                              uspp_augmentation_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'upf_core_kinetic'

   TYPE upf_core_kinetic_type
      LOGICAL :: available = .FALSE.
      CHARACTER(LEN=256) :: failure_reason = 'UPF model-core kinetic field has not been initialized'
      REAL(dp) :: support_radius = 0.0_dp
      REAL(dp), ALLOCATABLE :: radial_mesh(:), weighted(:)
      TYPE(uspp_augmentation_type), ALLOCATABLE :: model_field
   END TYPE upf_core_kinetic_type

   PUBLIC :: upf_core_kinetic_type, upf_core_kinetic_init, upf_core_kinetic_values
   PUBLIC :: upf_core_kinetic_cartesian, upf_core_kinetic_fourier

CONTAINS

! **************************************************************************************************
!> \brief Prepare the explicitly supplied UPF model-core tau on its original mesh and quadrature.
!>        Charge augmentation and the generalized metric are independent of this fixed field.
!> \param pot Original UPF data, including the declared model-core kinetic convention.
!> \param cache Shared fixed-field value, Fourier and geometry-derivative cache.
! **************************************************************************************************
   SUBROUTINE upf_core_kinetic_init(pot, cache)
      TYPE(atom_upfpot_type), INTENT(IN)                 :: pot
      TYPE(upf_core_kinetic_type), INTENT(OUT)           :: cache

      LOGICAL                                            :: available
      REAL(dp), ALLOCATABLE                              :: kinetic(:)

      CPASSERT(.NOT. pot%is_paw)
      CALL upf_nlcc_kinetic(pot, kinetic, available, cache%failure_reason)
      IF (.NOT. available) RETURN
      ALLOCATE (cache%radial_mesh, SOURCE=pot%r)
      ALLOCATE (cache%weighted, SOURCE=fourpi*upf_radial_weights(pot)*pot%r**2*kinetic)
      ALLOCATE (cache%model_field)
      CALL upf_radial_scalar_init(pot, kinetic, cache%model_field, logarithmic_origin=.TRUE.)
      cache%support_radius = cache%model_field%realspace_radius
      cache%available = .TRUE.
   END SUBROUTINE upf_core_kinetic_init

! **************************************************************************************************
!> \brief Evaluate the shared PS kinetic field and radial derivative at arbitrary nonnegative radii.
!> \param cache Initialized original-data cache.
!> \param radii Evaluation radii, including the physical smooth origin.
!> \param tau Spin-summed kinetic density including the factor 1/2.
!> \param derivative Optional derivative with respect to radius.
! **************************************************************************************************
   SUBROUTINE upf_core_kinetic_values(cache, radii, tau, derivative)
      TYPE(upf_core_kinetic_type), INTENT(IN)            :: cache
      REAL(dp), INTENT(IN)                               :: radii(:)
      REAL(dp), INTENT(OUT)                              :: tau(:)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: derivative(:)

      INTEGER                                            :: i
      REAL(dp)                                           :: slope

      CPASSERT(cache%available)
      CPASSERT(SIZE(tau) == SIZE(radii))
      CPASSERT(ALL(ieee_is_finite(radii)) .AND. ALL(radii >= 0.0_dp))
      IF (PRESENT(derivative)) THEN
         CPASSERT(SIZE(derivative) == SIZE(radii))
      END IF
      DO i = 1, SIZE(radii)
         CALL core_value(cache, radii(i), tau(i), slope)
         IF (PRESENT(derivative)) derivative(i) = slope
      END DO
   END SUBROUTINE upf_core_kinetic_values

! **************************************************************************************************
!> \brief Cartesian field and spatial gradient; center derivatives have the opposite sign.
!> \param cache Initialized original-data cache.
!> \param displacement Evaluation point minus the core center.
!> \param tau Spin-summed PS kinetic field.
!> \param gradient Optional spatial gradient, with the physical zero at the smooth origin.
! **************************************************************************************************
   SUBROUTINE upf_core_kinetic_cartesian(cache, displacement, tau, gradient)
      TYPE(upf_core_kinetic_type), INTENT(IN)            :: cache
      REAL(dp), INTENT(IN)                               :: displacement(3)
      REAL(dp), INTENT(OUT)                              :: tau
      REAL(dp), INTENT(OUT), OPTIONAL                    :: gradient(3)

      REAL(dp)                                           :: radius, slope

      CPASSERT(cache%available)
      CPASSERT(ALL(ieee_is_finite(displacement)))
      radius = SQRT(SUM(displacement**2))
      CALL core_value(cache, radius, tau, slope)
      IF (PRESENT(gradient)) THEN
         gradient = 0.0_dp
         IF (radius > 0.0_dp) gradient = slope*displacement/radius
      END IF
   END SUBROUTINE upf_core_kinetic_cartesian

! **************************************************************************************************
!> \brief Fourier field and exact adjoints of its original-mesh discrete radial transform.
!>        Under affine strain g transforms contravariantly, position covariantly, and volume
!>        with the determinant; the center phase is invariant in the strain derivative.
!> \param cache Initialized original-data cache.
!> \param g Cartesian reciprocal vector, including zero or arbitrary complex-wave sampling.
!> \param position Cartesian core center.
!> \param volume Cell volume.
!> \param coefficient exp(-i g.R) times the radial transform divided by volume.
!> \param strain Optional derivative with respect to affine strain epsilon_ab.
!> \param center_derivative Optional derivative with respect to the core center at fixed g.
! **************************************************************************************************
   SUBROUTINE upf_core_kinetic_fourier(cache, g, position, volume, coefficient, strain, center_derivative)
      TYPE(upf_core_kinetic_type), INTENT(IN)            :: cache
      REAL(dp), INTENT(IN)                               :: g(3), position(3), volume
      COMPLEX(dp), INTENT(OUT)                           :: coefficient
      COMPLEX(dp), INTENT(OUT), OPTIONAL                 :: strain(3, 3), center_derivative(3)

      COMPLEX(dp)                                        :: phase
      INTEGER                                            :: a, b
      REAL(dp)                                           :: magnitude, next, radial

      CPASSERT(cache%available .AND. volume > 0.0_dp)
      CPASSERT(ALL(ieee_is_finite(g)))
      CPASSERT(ALL(ieee_is_finite(position)))
      CPASSERT(ieee_is_finite(volume))
      magnitude = SQRT(SUM(g**2))
      radial = SUM(cache%weighted*bessel0(magnitude*cache%radial_mesh, 0))
      phase = EXP(CMPLX(0.0_dp, -DOT_PRODUCT(g, position), KIND=dp))/volume
      coefficient = phase*radial
      IF (PRESENT(center_derivative)) center_derivative = CMPLX(0.0_dp, -g, KIND=dp)*coefficient
      IF (PRESENT(strain)) THEN
         strain = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
         IF (magnitude > 0.0_dp) THEN
            next = SUM(cache%weighted*cache%radial_mesh*bessel0(magnitude*cache%radial_mesh, 1))
            DO b = 1, 3
               DO a = 1, 3
                  strain(a, b) = phase*g(a)*g(b)*next/magnitude
               END DO
            END DO
         END IF
         DO a = 1, 3
            strain(a, a) = strain(a, a) - coefficient
         END DO
      END IF
   END SUBROUTINE upf_core_kinetic_fourier

! **************************************************************************************************
!> \brief Evaluate the common Cartesian interpolant on a radial ray.
!> \param cache Explicit original model-core data.
!> \param radius Nonnegative evaluation radius.
!> \param value Scalar kinetic density.
!> \param derivative Radial derivative of the same field.
! **************************************************************************************************
   SUBROUTINE core_value(cache, radius, value, derivative)
      TYPE(upf_core_kinetic_type), INTENT(IN)            :: cache
      REAL(dp), INTENT(IN)                               :: radius
      REAL(dp), INTENT(OUT)                              :: value, derivative

      REAL(dp)                                           :: gradient(1, 1, 3), model(1, 1)

      CPASSERT(ALLOCATED(cache%model_field))
      CALL uspp_augmentation_realspace(cache%model_field, [radius, 0.0_dp, 0.0_dp], model, gradient)
      value = model(1, 1)
      derivative = gradient(1, 1, 1)
   END SUBROUTINE core_value

END MODULE upf_core_kinetic
