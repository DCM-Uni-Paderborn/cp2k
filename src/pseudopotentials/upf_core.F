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
   USE upf_types,                       ONLY: atom_upfpot_type,&
                                              upf_nlcc_kinetic,&
                                              upf_radial_weights
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: upf_core_fourier, upf_core_kinetic_fourier

CONTAINS

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
! **************************************************************************************************
   SUBROUTINE upf_core_fourier(pot, gvecs, position, volume, density, dposition, dg)
      TYPE(atom_upfpot_type), INTENT(IN)                 :: pot
      REAL(dp), INTENT(IN)                               :: gvecs(:, :), position(3), volume
      COMPLEX(dp), INTENT(OUT)                           :: density(:)
      COMPLEX(dp), INTENT(OUT), OPTIONAL                 :: dposition(:, :), dg(:, :)

      density = CMPLX(0.0_dp, 0.0_dp, dp)
      IF (PRESENT(dposition)) dposition = CMPLX(0.0_dp, 0.0_dp, dp)
      IF (PRESENT(dg)) dg = CMPLX(0.0_dp, 0.0_dp, dp)
      IF (.NOT. pot%core_correction) RETURN
      CPASSERT(ALLOCATED(pot%rho_nlcc))
      CALL core_field_fourier(pot, pot%rho_nlcc, gvecs, position, volume, density, dposition, dg)
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
! **************************************************************************************************
   SUBROUTINE upf_core_kinetic_fourier(pot, gvecs, position, volume, kinetic, dposition, dg, core_kinetic)
      TYPE(atom_upfpot_type), INTENT(IN)                 :: pot
      REAL(dp), INTENT(IN)                               :: gvecs(:, :), position(3), volume
      COMPLEX(dp), INTENT(OUT)                           :: kinetic(:)
      COMPLEX(dp), INTENT(OUT), OPTIONAL                 :: dposition(:, :), dg(:, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: core_kinetic(:)

      CHARACTER(LEN=256)                                 :: reason
      LOGICAL                                            :: available
      REAL(dp), ALLOCATABLE                              :: field(:)

      IF (PRESENT(core_kinetic)) THEN
         CPASSERT(SIZE(core_kinetic) == SIZE(pot%r))
         CALL core_field_fourier(pot, core_kinetic, gvecs, position, volume, kinetic, dposition, dg)
      ELSE
         CALL upf_nlcc_kinetic(pot, field, available, reason)
         IF (.NOT. available) CPABORT(TRIM(reason))
         CALL core_field_fourier(pot, field, gvecs, position, volume, kinetic, dposition, dg)
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
! **************************************************************************************************
   SUBROUTINE core_field_fourier(pot, field, gvecs, position, volume, density, dposition, dg)
      TYPE(atom_upfpot_type), INTENT(IN)                 :: pot
      REAL(dp), INTENT(IN)                               :: field(:), gvecs(:, :), position(3), &
                                                            volume
      COMPLEX(dp), INTENT(OUT)                           :: density(:)
      COMPLEX(dp), INTENT(OUT), OPTIONAL                 :: dposition(:, :), dg(:, :)

      COMPLEX(dp)                                        :: phase
      INTEGER                                            :: ig, nr
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
      nr = SIZE(pot%r)
      CPASSERT(SIZE(field) == nr)
      ALLOCATE (weighted(nr))
      weighted(:) = fourpi*upf_radial_weights(pot)*pot%r**2*field
      ! Form weights before omitting an exact zero suffix.
      DO WHILE (nr > 0)
         IF (weighted(nr) /= 0.0_dp) EXIT
         nr = nr - 1
      END DO
      IF (nr == 0) RETURN
      DO ig = 1, SIZE(density)
         norm_g = SQRT(SUM(gvecs(:, ig)**2))
         phase = EXP(CMPLX(0.0_dp, -DOT_PRODUCT(gvecs(:, ig), position), dp))/volume
         radial = SUM(weighted(:nr)*bessel0(norm_g*pot%r(:nr), 0))
         density(ig) = phase*radial
         IF (PRESENT(dposition)) dposition(:, ig) = CMPLX(0.0_dp, -gvecs(:, ig), dp)*density(ig)
         IF (PRESENT(dg)) THEN
            dg(:, ig) = CMPLX(0.0_dp, -position, dp)*density(ig)
            IF (norm_g > 0.0_dp) THEN
               slope = -SUM(weighted(:nr)*pot%r(:nr)*bessel0(norm_g*pot%r(:nr), 1))
               dg(:, ig) = dg(:, ig) + phase*slope*gvecs(:, ig)/norm_g
            END IF
         END IF
      END DO
   END SUBROUTINE core_field_fourier

END MODULE upf_core
