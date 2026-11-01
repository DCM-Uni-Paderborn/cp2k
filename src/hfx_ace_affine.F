!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Local complex ACE projectors and an affine frozen-operator energy.
!>        The interpolation factors include occupations. A mixed build density can be supplied
!>        independently of F F**H. Neither an orthogonal AO basis nor unit occupations are required.
!>        For a fixed compressed K,
!>        E_ref+Tr[K(D-D_ref)] has derivative K. Half of Tr[K D] does not.
!>        Compress valence-valence exchange separately from linear core exchange and XC.
MODULE hfx_ace_affine
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: dp
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: hfx_ace_affine_type, hfx_ace_affine_build, hfx_ace_affine_evaluate, hfx_ace_affine_release

   TYPE hfx_ace_affine_type
      COMPLEX(dp), ALLOCATABLE :: projectors(:, :)
      REAL(dp) :: offset = 0.0_dp
      LOGICAL :: ready = .FALSE.
   END TYPE hfx_ace_affine_type

CONTAINS

! **************************************************************************************************
!> \brief Build K_ACE=-Xi Xi**H with K_ACE F=K_ref F and an energy offset at the actual build density.
!>        With Y=K_ref F and G=-F**H Y=L L**H, Xi=Y L**(-H).
!>        Singular or inaccurate factorizations return an unusable cache. The caller must then
!>        evaluate full exchange. No eigenvalue shift or discarded occupied direction is used.
!> \param fock Complete unweighted Hermitian valence-valence exchange operator at the reference.
!> \param factors Occupied interpolation factors, including fractional occupations, in the same AO basis.
!> \param reference_energy Original valence-valence exchange energy at the build density, without k weight.
!> \param state New owned projectors and affine energy offset. A failed build releases all storage.
!> \param info Zero on success, a positive Cholesky failure, or -1 for excess factors / -2 for residual failure.
!> \param reference_density Optional actual mixed build density. Defaults to F F**H when absent.
!>        The supplied reference energy must correspond to this density, not necessarily to the factors.
! **************************************************************************************************
   SUBROUTINE hfx_ace_affine_build(fock, factors, reference_energy, state, info, reference_density)
      COMPLEX(dp), INTENT(IN)                            :: fock(:, :), factors(:, :)
      REAL(dp), INTENT(IN)                               :: reference_energy
      TYPE(hfx_ace_affine_type), INTENT(INOUT)           :: state
      INTEGER, INTENT(OUT)                               :: info
      COMPLEX(dp), INTENT(IN), OPTIONAL                  :: reference_density(:, :)

      COMPLEX(dp), PARAMETER :: one = CMPLX(1.0_dp, 0.0_dp, dp)

      COMPLEX(dp), ALLOCATABLE                           :: gram(:, :), image(:, :), overlaps(:, :), &
                                                            residual(:, :)
      INTEGER                                            :: n, rank
      LOGICAL                                            :: valid
      REAL(dp)                                           :: tolerance

      EXTERNAL zpotrf, ztrsm

      CALL hfx_ace_affine_release(state)
      n = SIZE(fock, 1)
      rank = SIZE(factors, 2)
      valid = n > 0 .AND. SIZE(fock, 2) == n .AND. SIZE(factors, 1) == n
      CPASSERT(valid)
      valid = ALL(ieee_is_finite(REAL(fock, dp))) .AND. ALL(ieee_is_finite(AIMAG(fock))) .AND. &
              ALL(ieee_is_finite(REAL(factors, dp))) .AND. ALL(ieee_is_finite(AIMAG(factors))) .AND. &
              ieee_is_finite(reference_energy)
      CPASSERT(valid)
      tolerance = 256.0_dp*EPSILON(1.0_dp)*REAL(n, dp)
      valid = MAXVAL(ABS(fock - CONJG(TRANSPOSE(fock)))) <= tolerance*MAX(1.0_dp, MAXVAL(ABS(fock)))
      CPASSERT(valid)
      IF (PRESENT(reference_density)) THEN
         valid = ALL(SHAPE(reference_density) == [n, n])
         CPASSERT(valid)
         valid = ALL(ieee_is_finite(REAL(reference_density, dp))) .AND. &
                 ALL(ieee_is_finite(AIMAG(reference_density)))
         CPASSERT(valid)
         valid = MAXVAL(ABS(reference_density - CONJG(TRANSPOSE(reference_density)))) <= &
            tolerance*MAX(1.0_dp, MAXVAL(ABS(reference_density)))
         CPASSERT(valid)
      END IF
      info = -1
      IF (rank > n) RETURN
      info = 0
      IF (rank == 0) THEN
         ALLOCATE (state%projectors(n, 0))
         state%offset = reference_energy
         state%ready = .TRUE.
         RETURN
      END IF
      ALLOCATE (image(n, rank), gram(rank, rank), overlaps(rank, rank), residual(n, rank))
      image(:, :) = MATMUL(fock, factors)
      gram(:, :) = -MATMUL(CONJG(TRANSPOSE(factors)), image)
      gram(:, :) = 0.5_dp*(gram + CONJG(TRANSPOSE(gram)))
      CALL zpotrf('L', rank, gram, rank, info)
      IF (info /= 0) RETURN
      ALLOCATE (state%projectors, SOURCE=image)
      CALL ztrsm('R', 'L', 'C', 'N', n, rank, one, gram, rank, state%projectors, n)
      overlaps(:, :) = MATMUL(CONJG(TRANSPOSE(state%projectors)), factors)
      residual(:, :) = -MATMUL(state%projectors, overlaps) - image
      IF (MAXVAL(ABS(residual)) > tolerance*MAX(1.0_dp, MAXVAL(ABS(image)))) THEN
         CALL hfx_ace_affine_release(state)
         info = -2
         RETURN
      END IF
      ! This offset preserves the independently supplied original reference energy.
      state%offset = reference_energy + SUM(ABS(overlaps)**2)
      IF (PRESENT(reference_density)) THEN
         image(:, :) = MATMUL(reference_density, state%projectors)
         state%offset = reference_energy + REAL(SUM(CONJG(state%projectors)*image), dp)
      END IF
      state%ready = .TRUE.
   END SUBROUTINE hfx_ace_affine_build

! **************************************************************************************************
!> \brief Evaluate the affine energy and its Hermitian AO derivative with fixed ACE projectors.
!>        Charge-ratio rescaling must multiply the complete affine energy, including its offset.
!>        Geometry derivatives require a fresh original exchange evaluation outside this routine.
!> \param state Previously built cache owned by this exchange context.
!> \param density Current Hermitian density in the same AO coordinates, without k weight.
!> \param energy Affine inner energy, equal to the original energy at the build density.
!> \param fock Its AO derivative. The physical mesh weight belongs to the caller.
! **************************************************************************************************
   SUBROUTINE hfx_ace_affine_evaluate(state, density, energy, fock)
      TYPE(hfx_ace_affine_type), INTENT(IN)              :: state
      COMPLEX(dp), INTENT(IN)                            :: density(:, :)
      REAL(dp), INTENT(OUT)                              :: energy
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: fock(:, :)

      INTEGER                                            :: n
      LOGICAL                                            :: valid

      CPASSERT(state%ready)
      CPASSERT(ALLOCATED(state%projectors))
      n = SIZE(state%projectors, 1)
      valid = ALL(SHAPE(density) == [n, n])
      CPASSERT(valid)
      ALLOCATE (fock(n, n))
      fock(:, :) = -MATMUL(state%projectors, CONJG(TRANSPOSE(state%projectors)))
      energy = state%offset + REAL(SUM(CONJG(fock)*density), dp)
   END SUBROUTINE hfx_ace_affine_evaluate

! **************************************************************************************************
!> \brief Release the projectors and reference energy of this cache only.
!> \param state Owned cache, also safe before initialization or after a failed build.
! **************************************************************************************************
   SUBROUTINE hfx_ace_affine_release(state)
      TYPE(hfx_ace_affine_type), INTENT(INOUT)           :: state

      IF (ALLOCATED(state%projectors)) DEALLOCATE (state%projectors)
      state%offset = 0.0_dp
      state%ready = .FALSE.
   END SUBROUTINE hfx_ace_affine_release

END MODULE hfx_ace_affine
