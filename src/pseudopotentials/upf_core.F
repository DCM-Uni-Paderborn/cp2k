!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
!> \brief Fourier coefficients and geometry derivatives of the original UPF partial core density.
MODULE upf_core
   USE bessel_lib,                      ONLY: bessel0
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: fourpi
   USE upf_radial_fields,               ONLY: upf_radial_fields_init,&
                                              upf_radial_fields_type
   USE upf_types,                       ONLY: atom_upfpot_type,&
                                              upf_nlcc_kinetic,&
                                              upf_radial_weights
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   TYPE upf_core_fourier_cache_type
      PRIVATE
      LOGICAL :: ready = .FALSE.
      REAL(dp), ALLOCATABLE :: radial(:), slope(:)
   END TYPE upf_core_fourier_cache_type
   PUBLIC :: upf_core_fourier, upf_core_kinetic_fourier, upf_core_radial_init, upf_core_fourier_cache_type

CONTAINS

! **************************************************************************************************
!> \brief Prepare a spherical core field on its original finite radial domain.
!>        Keep a nonzero final sample and the existing smooth closure for an exact zero suffix.
!> \param r Original radial mesh.
!> \param values Spin-summed spherical core charge or kinetic density.
!> \param field Cartesian radial spline with the spherical-harmonic normalization removed.
! **************************************************************************************************
   SUBROUTINE upf_core_radial_init(r, values, field)
      REAL(dp), INTENT(IN)                               :: r(:), values(:)
      TYPE(upf_radial_fields_type), INTENT(OUT)          :: field

      CPASSERT(SIZE(r) == SIZE(values) .AND. SIZE(r) > 0)
      CALL upf_radial_fields_init(r, RESHAPE(SQRT(fourpi)*values, [SIZE(values), 1]), [0], field, &
                                  finite_table=values(SIZE(values)) /= 0.0_dp)
   END SUBROUTINE upf_core_radial_init

! **************************************************************************************************
!> \brief Integrate the supplied true spherical density with the original radial quadrature.
!>        The finite tabulated integration domain is retained, including its endpoint rule.
!>        This field contributes to XC only, not the valence charge or orbital metric.
!> \param pot Original scalar UPF dataset.
!> \param gvecs Native reciprocal vectors in Cartesian inverse bohr.
!> \param position Physical atom position in bohr.
!> \param volume Physical cell volume in cubic bohr.
!> \param density Fourier density coefficients, including 1/volume.
!> \param dposition Optional Cartesian center derivatives, axis first.
!> \param dg Optional reciprocal-vector derivatives at fixed center and volume, axis first.
!> \param cache Optional radial transforms for this field, grid and derivative selection only.
! **************************************************************************************************
   SUBROUTINE upf_core_fourier(pot, gvecs, position, volume, density, dposition, dg, cache)
      TYPE(atom_upfpot_type), INTENT(IN)                 :: pot
      REAL(dp), INTENT(IN)                               :: gvecs(:, :), position(3), volume
      COMPLEX(dp), INTENT(OUT)                           :: density(:)
      COMPLEX(dp), INTENT(OUT), OPTIONAL                 :: dposition(:, :), dg(:, :)
      TYPE(upf_core_fourier_cache_type), OPTIONAL, &
         INTENT(INOUT)                                   :: cache

      density = CMPLX(0.0_dp, 0.0_dp, dp)
      IF (PRESENT(dposition)) dposition = CMPLX(0.0_dp, 0.0_dp, dp)
      IF (PRESENT(dg)) dg = CMPLX(0.0_dp, 0.0_dp, dp)
      IF (.NOT. pot%core_correction) RETURN
      CPASSERT(ALLOCATED(pot%rho_nlcc))
      CALL core_field_fourier(pot, pot%rho_nlcc, gvecs, position, volume, density, dposition, dg, cache)
   END SUBROUTINE upf_core_fourier

! **************************************************************************************************
!> \brief Fourier image and derivatives of the independently supplied canonical core tau.
!> \param pot Original scalar UPF dataset with an explicit core tau convention.
!> \param gvecs Native reciprocal vectors in Cartesian inverse bohr.
!> \param position Physical atom position in bohr.
!> \param volume Physical cell volume in cubic bohr.
!> \param kinetic Fourier kinetic-density coefficients, including 1/volume.
!> \param dposition Optional Cartesian center derivatives, axis first.
!> \param dg Optional reciprocal-vector derivatives at fixed center and volume, axis first.
!> \param core_kinetic Optional independently prepared core tau on the original mesh, for PAW PS cores.
!> \param cache Optional radial transforms for this field, grid and derivative selection only.
! **************************************************************************************************
   SUBROUTINE upf_core_kinetic_fourier(pot, gvecs, position, volume, kinetic, dposition, dg, core_kinetic, cache)
      TYPE(atom_upfpot_type), INTENT(IN)                 :: pot
      REAL(dp), INTENT(IN)                               :: gvecs(:, :), position(3), volume
      COMPLEX(dp), INTENT(OUT)                           :: kinetic(:)
      COMPLEX(dp), INTENT(OUT), OPTIONAL                 :: dposition(:, :), dg(:, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: core_kinetic(:)
      TYPE(upf_core_fourier_cache_type), OPTIONAL, &
         INTENT(INOUT)                                   :: cache

      CHARACTER(LEN=256)                                 :: reason
      LOGICAL                                            :: available
      REAL(dp), ALLOCATABLE                              :: field(:)

      IF (PRESENT(core_kinetic)) THEN
         CPASSERT(SIZE(core_kinetic) == SIZE(pot%r))
         CALL core_field_fourier(pot, core_kinetic, gvecs, position, volume, kinetic, dposition, dg, cache)
      ELSE
         CALL upf_nlcc_kinetic(pot, field, available, reason)
         IF (.NOT. available) CPABORT(TRIM(reason))
         CALL core_field_fourier(pot, field, gvecs, position, volume, kinetic, dposition, dg, cache)
      END IF
   END SUBROUTINE upf_core_kinetic_fourier

! **************************************************************************************************
!> \brief Shared radial transform for independently specified spherical core fields.
!> \param pot Original radial mesh and quadrature.
!> \param field True spherical field without a radial prefactor.
!> \param gvecs ...
!> \param position ...
!> \param volume ...
!> \param density ...
!> \param dposition ...
!> \param dg ...
!> \param cache ...
! **************************************************************************************************
   SUBROUTINE core_field_fourier(pot, field, gvecs, position, volume, density, dposition, dg, cache)
      TYPE(atom_upfpot_type), INTENT(IN)                 :: pot
      REAL(dp), INTENT(IN)                               :: field(:), gvecs(:, :), position(3), &
                                                            volume
      COMPLEX(dp), INTENT(OUT)                           :: density(:)
      COMPLEX(dp), INTENT(OUT), OPTIONAL                 :: dposition(:, :), dg(:, :)
      TYPE(upf_core_fourier_cache_type), OPTIONAL, &
         INTENT(INOUT)                                   :: cache

      COMPLEX(dp)                                        :: phase
      INTEGER                                            :: ig, nr
      LOGICAL                                            :: reuse
      REAL(dp)                                           :: norm_g, radial, slope
      REAL(dp), ALLOCATABLE                              :: weighted(:)

      CPASSERT(volume > 0.0_dp)
      CPASSERT(SIZE(gvecs, 1) == 3 .AND. SIZE(gvecs, 2) == SIZE(density))
      IF (PRESENT(dposition)) THEN
         CPASSERT(ALL(SHAPE(dposition) == SHAPE(gvecs)))
      END IF
      IF (PRESENT(dg)) THEN
         CPASSERT(ALL(SHAPE(dg) == SHAPE(gvecs)))
      END IF
      density = CMPLX(0.0_dp, 0.0_dp, dp)
      IF (PRESENT(dposition)) dposition = CMPLX(0.0_dp, 0.0_dp, dp)
      IF (PRESENT(dg)) dg = CMPLX(0.0_dp, 0.0_dp, dp)
      reuse = .FALSE.
      IF (PRESENT(cache)) reuse = cache%ready
      IF (reuse) THEN
         IF (.NOT. ALLOCATED(cache%radial)) RETURN
         CPASSERT(SIZE(cache%radial) == SIZE(density))
         CPASSERT(PRESENT(dg) .EQV. ALLOCATED(cache%slope))
      ELSE
         nr = SIZE(pot%r)
         CPASSERT(SIZE(field) == nr)
         ALLOCATE (weighted(nr))
         weighted(:) = fourpi*upf_radial_weights(pot)*pot%r**2*field
         ! Form weights before omitting an exact zero suffix.
         DO WHILE (nr > 0)
            IF (weighted(nr) /= 0.0_dp) EXIT
            nr = nr - 1
         END DO
         IF (PRESENT(cache)) cache%ready = .TRUE.
         IF (nr == 0) RETURN
         IF (PRESENT(cache)) THEN
            ALLOCATE (cache%radial(SIZE(density)))
            IF (PRESENT(dg)) ALLOCATE (cache%slope(SIZE(density)), SOURCE=0.0_dp)
         END IF
      END IF
      DO ig = 1, SIZE(density)
         norm_g = SQRT(SUM(gvecs(:, ig)**2))
         phase = EXP(CMPLX(0.0_dp, -DOT_PRODUCT(gvecs(:, ig), position), dp))/volume
         IF (reuse) THEN
            radial = cache%radial(ig)
         ELSE
            radial = SUM(weighted(:nr)*bessel0(norm_g*pot%r(:nr), 0))
            IF (PRESENT(cache)) cache%radial(ig) = radial
         END IF
         density(ig) = phase*radial
         IF (PRESENT(dposition)) dposition(:, ig) = CMPLX(0.0_dp, -gvecs(:, ig), dp)*density(ig)
         IF (PRESENT(dg)) THEN
            dg(:, ig) = CMPLX(0.0_dp, -position, dp)*density(ig)
            IF (norm_g > 0.0_dp) THEN
               IF (reuse) THEN
                  slope = cache%slope(ig)
               ELSE
                  slope = -SUM(weighted(:nr)*pot%r(:nr)*bessel0(norm_g*pot%r(:nr), 1))
                  IF (PRESENT(cache)) cache%slope(ig) = slope
               END IF
               dg(:, ig) = dg(:, ig) + phase*slope*gvecs(:, ig)/norm_g
            END IF
         END IF
      END DO
   END SUBROUTINE core_field_fourier

END MODULE upf_core
