!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

!--------------------------------------------------------------------------------------------------!
!   Copyright (c) 2008, 2009, Joost VandeVondele and Manuel Guidon                                 !
!   All rights reserved.                                                                           !
!                                                                                                  !
!   Redistribution and use in source and binary forms, with or without                             !
!   modification, are permitted provided that the following conditions are met:                    !
!       * Redistributions of source code must retain the above copyright                           !
!         notice, this list of conditions and the following disclaimer.                            !
!       * Redistributions in binary form must reproduce the above copyright                        !
!         notice, this list of conditions and the following disclaimer in the                      !
!         documentation and/or other materials provided with the distribution.                     !
!                                                                                                  !
!   THIS SOFTWARE IS PROVIDED BY Joost VandeVondele and Manuel Guidon AS IS AND ANY                !
!   EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED                      !
!   WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE                         !
!   DISCLAIMED. IN NO EVENT SHALL Joost VandeVondele or Manuel Guidon BE LIABLE FOR ANY            !
!   DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES                     !
!   (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;                   !
!   LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND                    !
!   ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT                     !
!   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS                  !
!   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                                   !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief This module computes the basic integrals for the truncated coulomb operator
!>
!>          res(1) =G_0(R,T)= ((2*erf(sqrt(t))+erf(R-sqrt(t))-erf(R+sqrt(t)))/sqrt(t))
!>
!>        and up to 21 derivatives with respect to T
!>
!>          res(n+1)=(-1)**n d^n/dT^n G_0(R,T)
!>
!>        The function is only computed for values of R,T which fulfil
!>
!>          R**2 - 11.0_dp*R + 0.0_dp < T < R**2 + 11.0_dp*R + 50.0_dp where R>=0 T>=0
!>
!>        for T larger than the upper bound, 0 is returned
!>        (which is accurate at least up to 1.0E-16)
!>        while for T smaller than the lower bound, the caller is instructed
!>        to use the conventional gamma function instead
!>        (i.e. the limit of above expression for R to Infinity)
!>
!> \author Joost VandeVondele and Manuel Guidon
!> \par History
!>      Nov 2008, 2009 Joost VandeVondele and Manuel Guidon
!>      May 2019 A. Bussy: Added a get_maxl_init function to get current status of nderiv_init and
!>                moved the file to common (made it accessible from aobasis, same place as gamma.F).
!>      Oct 2025 M. Puligheddu: Added public qualifier to C0 to simplify reuse
! **************************************************************************************************
MODULE t_c_g0
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_comm_type
#include "../base/base_uses.f90"

   IMPLICIT NONE

   REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE, SAVE, PUBLIC :: C0

   PRIVATE

   PUBLIC :: t_c_g0_n, t_c_g0_response, init, free_C0, get_lmax_init

   INTEGER, PARAMETER :: degree = 13
   REAL(KIND=dp), PARAMETER :: target_error = 1.000000E-09_dp
   INTEGER, PARAMETER :: nderiv_max = 21
   INTEGER, SAVE      :: nderiv_init = -1
   INTEGER, SAVE      :: patches = -1

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param RES ...
!> \param use_gamma ...
!> \param R ...
!> \param T ...
!> \param NDERIV ...
! **************************************************************************************************
   SUBROUTINE t_c_g0_n(RES, use_gamma, R, T, NDERIV)
      REAL(dp), INTENT(OUT)                              :: RES(*)
      LOGICAL, INTENT(OUT)                               :: use_gamma
      REAL(dp), INTENT(IN)                               :: R, T
      INTEGER, INTENT(IN)                                :: NDERIV

      CALL evaluate_moments(RES, use_gamma, R, T, NDERIV)
   END SUBROUTINE t_c_g0_n

! **************************************************************************************************
!> \brief Analytic derivative of the original moment interpolants with respect to T.
!> \param RES ...
!> \param dres_dt ...
!> \param use_gamma ...
!> \param R ...
!> \param T ...
!> \param NDERIV ...
!> \note R and the coefficient table are held fixed. If use_gamma is true, the caller
!>       uses ordinary Boys moments and their exact recurrence instead.
! **************************************************************************************************
   SUBROUTINE t_c_g0_response(RES, dres_dt, use_gamma, R, T, NDERIV)
      REAL(dp), INTENT(OUT)                              :: RES(*), dres_dt(*)
      LOGICAL, INTENT(OUT)                               :: use_gamma
      REAL(dp), INTENT(IN)                               :: R, T
      INTEGER, INTENT(IN)                                :: NDERIV

      CALL evaluate_moments(RES, use_gamma, R, T, NDERIV, dres_dt)
   END SUBROUTINE t_c_g0_response

! **************************************************************************************************
!> \brief ...
!> \param x ...
!> \param center ...
!> \param factor ...
!> \param coordinate ...
!> \param derivative ...
! **************************************************************************************************
   PURE SUBROUTINE patch_coordinate(x, center, factor, coordinate, derivative)
      REAL(dp), INTENT(IN)                               :: x, center, factor
      REAL(dp), INTENT(OUT)                              :: coordinate, derivative

      coordinate = (2*x - center)*factor
      derivative = 2*factor
   END SUBROUTINE patch_coordinate

! **************************************************************************************************
!> \brief ...
!> \param RES ...
!> \param use_gamma ...
!> \param R ...
!> \param T ...
!> \param NDERIV ...
!> \param dres_dt ...
! **************************************************************************************************
   SUBROUTINE evaluate_moments(RES, use_gamma, R, T, NDERIV, dres_dt)
      REAL(KIND=dp), INTENT(OUT)                         :: RES(*)
      LOGICAL, INTENT(OUT)                               :: use_gamma
      REAL(KIND=dp), INTENT(IN)                          :: R, T
      INTEGER, INTENT(IN)                                :: NDERIV
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: dres_dt(*)

      INTEGER                                            :: patch
      REAL(KIND=dp)                                      :: lower, scale, TG1, TG2, upper, X1, X2

      IF (PRESENT(dres_dt)) dres_dt(1:NDERIV + 1) = 0._dp
      use_gamma = .FALSE.
      upper = R**2 + 11.0_dp*R + 50.0_dp
      lower = R**2 - 11.0_dp*R + 0.0_dp
      IF (T > upper) THEN
         RES(1:NDERIV + 1) = 0.0_dp
         RETURN
      END IF
      IF (R <= 11.0_dp) THEN
         X2 = R/11.0_dp
         upper = R**2 + 11.0_dp*R + 50.0_dp
         lower = 0.0_dp
         X1 = (T - lower)/(upper - lower)
         IF (X1 <= 0.500000000000000000E+00_dp) THEN
            IF (X2 <= 0.500000000000000000E+00_dp) THEN
               IF (X2 <= 0.250000000000000000E+00_dp) THEN
                  IF (X2 <= 0.125000000000000000E+00_dp) THEN
                     IF (X1 <= 0.250000000000000000E+00_dp) THEN
                        IF (X2 <= 0.625000000000000000E-01_dp) THEN
                           IF (X1 <= 0.125000000000000000E+00_dp) THEN
                              IF (X2 <= 0.312500000000000000E-01_dp) THEN
                                 IF (X1 <= 0.625000000000000000E-01_dp) THEN
                                    IF (X2 <= 0.156250000000000000E-01_dp) THEN
                                       IF (X1 <= 0.312500000000000000E-01_dp) THEN
                                          CALL patch_coordinate(X1, 0.03125_dp, 32.0_dp, TG1, scale)
                                          TG2 = (2*X2 - 0.156250000000000000E-01_dp)*0.640000000000000000E+02_dp
                                          patch = 1
                                       ELSE
                                          CALL patch_coordinate(X1, 0.09375_dp, 32.0_dp, TG1, scale)
                                          TG2 = (2*X2 - 0.156250000000000000E-01_dp)*0.640000000000000000E+02_dp
                                          patch = 2
                                       END IF
                                    ELSE
                                       IF (X1 <= 0.312500000000000000E-01_dp) THEN
                                          CALL patch_coordinate(X1, 0.03125_dp, 32.0_dp, TG1, scale)
                                          TG2 = (2*X2 - 0.468750000000000000E-01_dp)*0.640000000000000000E+02_dp
                                          patch = 3
                                       ELSE
                                          CALL patch_coordinate(X1, 0.09375_dp, 32.0_dp, TG1, scale)
                                          TG2 = (2*X2 - 0.468750000000000000E-01_dp)*0.640000000000000000E+02_dp
                                          patch = 4
                                       END IF
                                    END IF
                                 ELSE
                                    IF (X2 <= 0.156250000000000000E-01_dp) THEN
                                       CALL patch_coordinate(X1, 0.1875_dp, 16.0_dp, TG1, scale)
                                       TG2 = (2*X2 - 0.156250000000000000E-01_dp)*0.640000000000000000E+02_dp
                                       patch = 5
                                    ELSE
                                       CALL patch_coordinate(X1, 0.1875_dp, 16.0_dp, TG1, scale)
                                       TG2 = (2*X2 - 0.468750000000000000E-01_dp)*0.640000000000000000E+02_dp
                                       patch = 6
                                    END IF
                                 END IF
                              ELSE
                                 IF (X1 <= 0.625000000000000000E-01_dp) THEN
                                    IF (X2 <= 0.468750000000000000E-01_dp) THEN
                                       IF (X1 <= 0.312500000000000000E-01_dp) THEN
                                          CALL patch_coordinate(X1, 0.03125_dp, 32.0_dp, TG1, scale)
                                          TG2 = (2*X2 - 0.781250000000000000E-01_dp)*0.640000000000000000E+02_dp
                                          patch = 7
                                       ELSE
                                          CALL patch_coordinate(X1, 0.09375_dp, 32.0_dp, TG1, scale)
                                          TG2 = (2*X2 - 0.781250000000000000E-01_dp)*0.640000000000000000E+02_dp
                                          patch = 8
                                       END IF
                                    ELSE
                                       IF (X1 <= 0.312500000000000000E-01_dp) THEN
                                          CALL patch_coordinate(X1, 0.03125_dp, 32.0_dp, TG1, scale)
                                          TG2 = (2*X2 - 0.109375000000000000E+00_dp)*0.640000000000000000E+02_dp
                                          patch = 9
                                       ELSE
                                          CALL patch_coordinate(X1, 0.09375_dp, 32.0_dp, TG1, scale)
                                          TG2 = (2*X2 - 0.109375000000000000E+00_dp)*0.640000000000000000E+02_dp
                                          patch = 10
                                       END IF
                                    END IF
                                 ELSE
                                    IF (X2 <= 0.468750000000000000E-01_dp) THEN
                                       CALL patch_coordinate(X1, 0.1875_dp, 16.0_dp, TG1, scale)
                                       TG2 = (2*X2 - 0.781250000000000000E-01_dp)*0.640000000000000000E+02_dp
                                       patch = 11
                                    ELSE
                                       CALL patch_coordinate(X1, 0.1875_dp, 16.0_dp, TG1, scale)
                                       TG2 = (2*X2 - 0.109375000000000000E+00_dp)*0.640000000000000000E+02_dp
                                       patch = 12
                                    END IF
                                 END IF
                              END IF
                           ELSE
                              IF (X2 <= 0.312500000000000000E-01_dp) THEN
                                 IF (X1 <= 0.187500000000000000E+00_dp) THEN
                                    CALL patch_coordinate(X1, 0.3125_dp, 16.0_dp, TG1, scale)
                                    TG2 = (2*X2 - 0.312500000000000000E-01_dp)*0.320000000000000000E+02_dp
                                    patch = 13
                                 ELSE
                                    CALL patch_coordinate(X1, 0.4375_dp, 16.0_dp, TG1, scale)
                                    TG2 = (2*X2 - 0.312500000000000000E-01_dp)*0.320000000000000000E+02_dp
                                    patch = 14
                                 END IF
                              ELSE
                                 IF (X1 <= 0.187500000000000000E+00_dp) THEN
                                    CALL patch_coordinate(X1, 0.3125_dp, 16.0_dp, TG1, scale)
                                    TG2 = (2*X2 - 0.937500000000000000E-01_dp)*0.320000000000000000E+02_dp
                                    patch = 15
                                 ELSE
                                    CALL patch_coordinate(X1, 0.4375_dp, 16.0_dp, TG1, scale)
                                    TG2 = (2*X2 - 0.937500000000000000E-01_dp)*0.320000000000000000E+02_dp
                                    patch = 16
                                 END IF
                              END IF
                           END IF
                        ELSE
                           IF (X1 <= 0.125000000000000000E+00_dp) THEN
                              IF (X2 <= 0.937500000000000000E-01_dp) THEN
                                 IF (X1 <= 0.625000000000000000E-01_dp) THEN
                                    IF (X2 <= 0.781250000000000000E-01_dp) THEN
                                       IF (X1 <= 0.312500000000000000E-01_dp) THEN
                                          CALL patch_coordinate(X1, 0.03125_dp, 32.0_dp, TG1, scale)
                                          TG2 = (2*X2 - 0.140625000000000000E+00_dp)*0.640000000000000000E+02_dp
                                          patch = 17
                                       ELSE
                                          CALL patch_coordinate(X1, 0.09375_dp, 32.0_dp, TG1, scale)
                                          TG2 = (2*X2 - 0.140625000000000000E+00_dp)*0.640000000000000000E+02_dp
                                          patch = 18
                                       END IF
                                    ELSE
                                       IF (X1 <= 0.312500000000000000E-01_dp) THEN
                                          CALL patch_coordinate(X1, 0.03125_dp, 32.0_dp, TG1, scale)
                                          TG2 = (2*X2 - 0.171875000000000000E+00_dp)*0.640000000000000000E+02_dp
                                          patch = 19
                                       ELSE
                                          CALL patch_coordinate(X1, 0.09375_dp, 32.0_dp, TG1, scale)
                                          TG2 = (2*X2 - 0.171875000000000000E+00_dp)*0.640000000000000000E+02_dp
                                          patch = 20
                                       END IF
                                    END IF
                                 ELSE
                                    IF (X2 <= 0.781250000000000000E-01_dp) THEN
                                       CALL patch_coordinate(X1, 0.1875_dp, 16.0_dp, TG1, scale)
                                       TG2 = (2*X2 - 0.140625000000000000E+00_dp)*0.640000000000000000E+02_dp
                                       patch = 21
                                    ELSE
                                       CALL patch_coordinate(X1, 0.1875_dp, 16.0_dp, TG1, scale)
                                       TG2 = (2*X2 - 0.171875000000000000E+00_dp)*0.640000000000000000E+02_dp
                                       patch = 22
                                    END IF
                                 END IF
                              ELSE
                                 IF (X1 <= 0.625000000000000000E-01_dp) THEN
                                    IF (X2 <= 0.109375000000000000E+00_dp) THEN
                                       IF (X1 <= 0.312500000000000000E-01_dp) THEN
                                          CALL patch_coordinate(X1, 0.03125_dp, 32.0_dp, TG1, scale)
                                          TG2 = (2*X2 - 0.203125000000000000E+00_dp)*0.640000000000000000E+02_dp
                                          patch = 23
                                       ELSE
                                          CALL patch_coordinate(X1, 0.09375_dp, 32.0_dp, TG1, scale)
                                          TG2 = (2*X2 - 0.203125000000000000E+00_dp)*0.640000000000000000E+02_dp
                                          patch = 24
                                       END IF
                                    ELSE
                                       IF (X1 <= 0.312500000000000000E-01_dp) THEN
                                          CALL patch_coordinate(X1, 0.03125_dp, 32.0_dp, TG1, scale)
                                          TG2 = (2*X2 - 0.234375000000000000E+00_dp)*0.640000000000000000E+02_dp
                                          patch = 25
                                       ELSE
                                          CALL patch_coordinate(X1, 0.09375_dp, 32.0_dp, TG1, scale)
                                          TG2 = (2*X2 - 0.234375000000000000E+00_dp)*0.640000000000000000E+02_dp
                                          patch = 26
                                       END IF
                                    END IF
                                 ELSE
                                    IF (X2 <= 0.109375000000000000E+00_dp) THEN
                                       CALL patch_coordinate(X1, 0.1875_dp, 16.0_dp, TG1, scale)
                                       TG2 = (2*X2 - 0.203125000000000000E+00_dp)*0.640000000000000000E+02_dp
                                       patch = 27
                                    ELSE
                                       CALL patch_coordinate(X1, 0.1875_dp, 16.0_dp, TG1, scale)
                                       TG2 = (2*X2 - 0.234375000000000000E+00_dp)*0.640000000000000000E+02_dp
                                       patch = 28
                                    END IF
                                 END IF
                              END IF
                           ELSE
                              IF (X1 <= 0.187500000000000000E+00_dp) THEN
                                 IF (X2 <= 0.937500000000000000E-01_dp) THEN
                                    CALL patch_coordinate(X1, 0.3125_dp, 16.0_dp, TG1, scale)
                                    TG2 = (2*X2 - 0.156250000000000000E+00_dp)*0.320000000000000000E+02_dp
                                    patch = 29
                                 ELSE
                                    CALL patch_coordinate(X1, 0.3125_dp, 16.0_dp, TG1, scale)
                                    TG2 = (2*X2 - 0.218750000000000000E+00_dp)*0.320000000000000000E+02_dp
                                    patch = 30
                                 END IF
                              ELSE
                                 CALL patch_coordinate(X1, 0.4375_dp, 16.0_dp, TG1, scale)
                                 TG2 = (2*X2 - 0.187500000000000000E+00_dp)*0.160000000000000000E+02_dp
                                 patch = 31
                              END IF
                           END IF
                        END IF
                     ELSE
                        IF (X1 <= 0.375000000000000000E+00_dp) THEN
                           CALL patch_coordinate(X1, 0.625_dp, 8.0_dp, TG1, scale)
                           TG2 = (2*X2 - 0.125000000000000000E+00_dp)*0.800000000000000000E+01_dp
                           patch = 32
                        ELSE
                           CALL patch_coordinate(X1, 0.875_dp, 8.0_dp, TG1, scale)
                           TG2 = (2*X2 - 0.125000000000000000E+00_dp)*0.800000000000000000E+01_dp
                           patch = 33
                        END IF
                     END IF
                  ELSE
                     IF (X1 <= 0.250000000000000000E+00_dp) THEN
                        IF (X2 <= 0.187500000000000000E+00_dp) THEN
                           IF (X1 <= 0.125000000000000000E+00_dp) THEN
                              IF (X2 <= 0.156250000000000000E+00_dp) THEN
                                 IF (X1 <= 0.625000000000000000E-01_dp) THEN
                                    IF (X1 <= 0.312500000000000000E-01_dp) THEN
                                       IF (X2 <= 0.140625000000000000E+00_dp) THEN
                                          CALL patch_coordinate(X1, 0.03125_dp, 32.0_dp, TG1, scale)
                                          TG2 = (2*X2 - 0.265625000000000000E+00_dp)*0.640000000000000000E+02_dp
                                          patch = 34
                                       ELSE
                                          CALL patch_coordinate(X1, 0.03125_dp, 32.0_dp, TG1, scale)
                                          TG2 = (2*X2 - 0.296875000000000000E+00_dp)*0.640000000000000000E+02_dp
                                          patch = 35
                                       END IF
                                    ELSE
                                       CALL patch_coordinate(X1, 0.09375_dp, 32.0_dp, TG1, scale)
                                       TG2 = (2*X2 - 0.281250000000000000E+00_dp)*0.320000000000000000E+02_dp
                                       patch = 36
                                    END IF
                                 ELSE
                                    IF (X1 <= 0.937500000000000000E-01_dp) THEN
                                       CALL patch_coordinate(X1, 0.15625_dp, 32.0_dp, TG1, scale)
                                       TG2 = (2*X2 - 0.281250000000000000E+00_dp)*0.320000000000000000E+02_dp
                                       patch = 37
                                    ELSE
                                       CALL patch_coordinate(X1, 0.21875_dp, 32.0_dp, TG1, scale)
                                       TG2 = (2*X2 - 0.281250000000000000E+00_dp)*0.320000000000000000E+02_dp
                                       patch = 38
                                    END IF
                                 END IF
                              ELSE
                                 IF (X1 <= 0.625000000000000000E-01_dp) THEN
                                    IF (X1 <= 0.312500000000000000E-01_dp) THEN
                                       IF (X2 <= 0.171875000000000000E+00_dp) THEN
                                          CALL patch_coordinate(X1, 0.03125_dp, 32.0_dp, TG1, scale)
                                          TG2 = (2*X2 - 0.328125000000000000E+00_dp)*0.640000000000000000E+02_dp
                                          patch = 39
                                       ELSE
                                          CALL patch_coordinate(X1, 0.03125_dp, 32.0_dp, TG1, scale)
                                          TG2 = (2*X2 - 0.359375000000000000E+00_dp)*0.640000000000000000E+02_dp
                                          patch = 40
                                       END IF
                                    ELSE
                                       CALL patch_coordinate(X1, 0.09375_dp, 32.0_dp, TG1, scale)
                                       TG2 = (2*X2 - 0.343750000000000000E+00_dp)*0.320000000000000000E+02_dp
                                       patch = 41
                                    END IF
                                 ELSE
                                    CALL patch_coordinate(X1, 0.1875_dp, 16.0_dp, TG1, scale)
                                    TG2 = (2*X2 - 0.343750000000000000E+00_dp)*0.320000000000000000E+02_dp
                                    patch = 42
                                 END IF
                              END IF
                           ELSE
                              IF (X1 <= 0.187500000000000000E+00_dp) THEN
                                 CALL patch_coordinate(X1, 0.3125_dp, 16.0_dp, TG1, scale)
                                 TG2 = (2*X2 - 0.312500000000000000E+00_dp)*0.160000000000000000E+02_dp
                                 patch = 43
                              ELSE
                                 CALL patch_coordinate(X1, 0.4375_dp, 16.0_dp, TG1, scale)
                                 TG2 = (2*X2 - 0.312500000000000000E+00_dp)*0.160000000000000000E+02_dp
                                 patch = 44
                              END IF
                           END IF
                        ELSE
                           IF (X1 <= 0.125000000000000000E+00_dp) THEN
                              IF (X1 <= 0.625000000000000000E-01_dp) THEN
                                 IF (X2 <= 0.218750000000000000E+00_dp) THEN
                                    IF (X1 <= 0.312500000000000000E-01_dp) THEN
                                       IF (X2 <= 0.203125000000000000E+00_dp) THEN
                                          CALL patch_coordinate(X1, 0.03125_dp, 32.0_dp, TG1, scale)
                                          TG2 = (2*X2 - 0.390625000000000000E+00_dp)*0.640000000000000000E+02_dp
                                          patch = 45
                                       ELSE
                                          CALL patch_coordinate(X1, 0.03125_dp, 32.0_dp, TG1, scale)
                                          TG2 = (2*X2 - 0.421875000000000000E+00_dp)*0.640000000000000000E+02_dp
                                          patch = 46
                                       END IF
                                    ELSE
                                       CALL patch_coordinate(X1, 0.09375_dp, 32.0_dp, TG1, scale)
                                       TG2 = (2*X2 - 0.406250000000000000E+00_dp)*0.320000000000000000E+02_dp
                                       patch = 47
                                    END IF
                                 ELSE
                                    IF (X1 <= 0.312500000000000000E-01_dp) THEN
                                       IF (X2 <= 0.234375000000000000E+00_dp) THEN
                                          CALL patch_coordinate(X1, 0.03125_dp, 32.0_dp, TG1, scale)
                                          TG2 = (2*X2 - 0.453125000000000000E+00_dp)*0.640000000000000000E+02_dp
                                          patch = 48
                                       ELSE
                                          CALL patch_coordinate(X1, 0.03125_dp, 32.0_dp, TG1, scale)
                                          TG2 = (2*X2 - 0.484375000000000000E+00_dp)*0.640000000000000000E+02_dp
                                          patch = 49
                                       END IF
                                    ELSE
                                       CALL patch_coordinate(X1, 0.09375_dp, 32.0_dp, TG1, scale)
                                       TG2 = (2*X2 - 0.468750000000000000E+00_dp)*0.320000000000000000E+02_dp
                                       patch = 50
                                    END IF
                                 END IF
                              ELSE
                                 IF (X2 <= 0.218750000000000000E+00_dp) THEN
                                    CALL patch_coordinate(X1, 0.1875_dp, 16.0_dp, TG1, scale)
                                    TG2 = (2*X2 - 0.406250000000000000E+00_dp)*0.320000000000000000E+02_dp
                                    patch = 51
                                 ELSE
                                    CALL patch_coordinate(X1, 0.1875_dp, 16.0_dp, TG1, scale)
                                    TG2 = (2*X2 - 0.468750000000000000E+00_dp)*0.320000000000000000E+02_dp
                                    patch = 52
                                 END IF
                              END IF
                           ELSE
                              IF (X1 <= 0.187500000000000000E+00_dp) THEN
                                 CALL patch_coordinate(X1, 0.3125_dp, 16.0_dp, TG1, scale)
                                 TG2 = (2*X2 - 0.437500000000000000E+00_dp)*0.160000000000000000E+02_dp
                                 patch = 53
                              ELSE
                                 CALL patch_coordinate(X1, 0.4375_dp, 16.0_dp, TG1, scale)
                                 TG2 = (2*X2 - 0.437500000000000000E+00_dp)*0.160000000000000000E+02_dp
                                 patch = 54
                              END IF
                           END IF
                        END IF
                     ELSE
                        IF (X1 <= 0.375000000000000000E+00_dp) THEN
                           IF (X1 <= 0.312500000000000000E+00_dp) THEN
                              CALL patch_coordinate(X1, 0.5625_dp, 16.0_dp, TG1, scale)
                              TG2 = (2*X2 - 0.375000000000000000E+00_dp)*0.800000000000000000E+01_dp
                              patch = 55
                           ELSE
                              CALL patch_coordinate(X1, 0.6875_dp, 16.0_dp, TG1, scale)
                              TG2 = (2*X2 - 0.375000000000000000E+00_dp)*0.800000000000000000E+01_dp
                              patch = 56
                           END IF
                        ELSE
                           CALL patch_coordinate(X1, 0.875_dp, 8.0_dp, TG1, scale)
                           TG2 = (2*X2 - 0.375000000000000000E+00_dp)*0.800000000000000000E+01_dp
                           patch = 57
                        END IF
                     END IF
                  END IF
               ELSE
                  IF (X1 <= 0.250000000000000000E+00_dp) THEN
                     IF (X1 <= 0.125000000000000000E+00_dp) THEN
                        IF (X1 <= 0.625000000000000000E-01_dp) THEN
                           IF (X2 <= 0.375000000000000000E+00_dp) THEN
                              IF (X2 <= 0.312500000000000000E+00_dp) THEN
                                 IF (X1 <= 0.312500000000000000E-01_dp) THEN
                                    IF (X2 <= 0.281250000000000000E+00_dp) THEN
                                       CALL patch_coordinate(X1, 0.03125_dp, 32.0_dp, TG1, scale)
                                       TG2 = (2*X2 - 0.531250000000000000E+00_dp)*0.320000000000000000E+02_dp
                                       patch = 58
                                    ELSE
                                       CALL patch_coordinate(X1, 0.03125_dp, 32.0_dp, TG1, scale)
                                       TG2 = (2*X2 - 0.593750000000000000E+00_dp)*0.320000000000000000E+02_dp
                                       patch = 59
                                    END IF
                                 ELSE
                                    IF (X2 <= 0.281250000000000000E+00_dp) THEN
                                       CALL patch_coordinate(X1, 0.09375_dp, 32.0_dp, TG1, scale)
                                       TG2 = (2*X2 - 0.531250000000000000E+00_dp)*0.320000000000000000E+02_dp
                                       patch = 60
                                    ELSE
                                       CALL patch_coordinate(X1, 0.09375_dp, 32.0_dp, TG1, scale)
                                       TG2 = (2*X2 - 0.593750000000000000E+00_dp)*0.320000000000000000E+02_dp
                                       patch = 61
                                    END IF
                                 END IF
                              ELSE
                                 IF (X1 <= 0.312500000000000000E-01_dp) THEN
                                    IF (X2 <= 0.343750000000000000E+00_dp) THEN
                                       CALL patch_coordinate(X1, 0.03125_dp, 32.0_dp, TG1, scale)
                                       TG2 = (2*X2 - 0.656250000000000000E+00_dp)*0.320000000000000000E+02_dp
                                       patch = 62
                                    ELSE
                                       CALL patch_coordinate(X1, 0.03125_dp, 32.0_dp, TG1, scale)
                                       TG2 = (2*X2 - 0.718750000000000000E+00_dp)*0.320000000000000000E+02_dp
                                       patch = 63
                                    END IF
                                 ELSE
                                    CALL patch_coordinate(X1, 0.09375_dp, 32.0_dp, TG1, scale)
                                    TG2 = (2*X2 - 0.687500000000000000E+00_dp)*0.160000000000000000E+02_dp
                                    patch = 64
                                 END IF
                              END IF
                           ELSE
                              IF (X1 <= 0.312500000000000000E-01_dp) THEN
                                 IF (X2 <= 0.437500000000000000E+00_dp) THEN
                                    IF (X1 <= 0.156250000000000000E-01_dp) THEN
                                       CALL patch_coordinate(X1, 0.015625_dp, 64.0_dp, TG1, scale)
                                       TG2 = (2*X2 - 0.812500000000000000E+00_dp)*0.160000000000000000E+02_dp
                                       patch = 65
                                    ELSE
                                       CALL patch_coordinate(X1, 0.046875_dp, 64.0_dp, TG1, scale)
                                       TG2 = (2*X2 - 0.812500000000000000E+00_dp)*0.160000000000000000E+02_dp
                                       patch = 66
                                    END IF
                                 ELSE
                                    IF (X1 <= 0.156250000000000000E-01_dp) THEN
                                       CALL patch_coordinate(X1, 0.015625_dp, 64.0_dp, TG1, scale)
                                       TG2 = (2*X2 - 0.937500000000000000E+00_dp)*0.160000000000000000E+02_dp
                                       patch = 67
                                    ELSE
                                       CALL patch_coordinate(X1, 0.046875_dp, 64.0_dp, TG1, scale)
                                       TG2 = (2*X2 - 0.937500000000000000E+00_dp)*0.160000000000000000E+02_dp
                                       patch = 68
                                    END IF
                                 END IF
                              ELSE
                                 IF (X2 <= 0.437500000000000000E+00_dp) THEN
                                    CALL patch_coordinate(X1, 0.09375_dp, 32.0_dp, TG1, scale)
                                    TG2 = (2*X2 - 0.812500000000000000E+00_dp)*0.160000000000000000E+02_dp
                                    patch = 69
                                 ELSE
                                    CALL patch_coordinate(X1, 0.09375_dp, 32.0_dp, TG1, scale)
                                    TG2 = (2*X2 - 0.937500000000000000E+00_dp)*0.160000000000000000E+02_dp
                                    patch = 70
                                 END IF
                              END IF
                           END IF
                        ELSE
                           IF (X2 <= 0.375000000000000000E+00_dp) THEN
                              IF (X2 <= 0.312500000000000000E+00_dp) THEN
                                 IF (X1 <= 0.937500000000000000E-01_dp) THEN
                                    CALL patch_coordinate(X1, 0.15625_dp, 32.0_dp, TG1, scale)
                                    TG2 = (2*X2 - 0.562500000000000000E+00_dp)*0.160000000000000000E+02_dp
                                    patch = 71
                                 ELSE
                                    CALL patch_coordinate(X1, 0.21875_dp, 32.0_dp, TG1, scale)
                                    TG2 = (2*X2 - 0.562500000000000000E+00_dp)*0.160000000000000000E+02_dp
                                    patch = 72
                                 END IF
                              ELSE
                                 IF (X1 <= 0.937500000000000000E-01_dp) THEN
                                    CALL patch_coordinate(X1, 0.15625_dp, 32.0_dp, TG1, scale)
                                    TG2 = (2*X2 - 0.687500000000000000E+00_dp)*0.160000000000000000E+02_dp
                                    patch = 73
                                 ELSE
                                    CALL patch_coordinate(X1, 0.21875_dp, 32.0_dp, TG1, scale)
                                    TG2 = (2*X2 - 0.687500000000000000E+00_dp)*0.160000000000000000E+02_dp
                                    patch = 74
                                 END IF
                              END IF
                           ELSE
                              IF (X1 <= 0.937500000000000000E-01_dp) THEN
                                 IF (X2 <= 0.437500000000000000E+00_dp) THEN
                                    CALL patch_coordinate(X1, 0.15625_dp, 32.0_dp, TG1, scale)
                                    TG2 = (2*X2 - 0.812500000000000000E+00_dp)*0.160000000000000000E+02_dp
                                    patch = 75
                                 ELSE
                                    CALL patch_coordinate(X1, 0.15625_dp, 32.0_dp, TG1, scale)
                                    TG2 = (2*X2 - 0.937500000000000000E+00_dp)*0.160000000000000000E+02_dp
                                    patch = 76
                                 END IF
                              ELSE
                                 IF (X2 <= 0.437500000000000000E+00_dp) THEN
                                    CALL patch_coordinate(X1, 0.21875_dp, 32.0_dp, TG1, scale)
                                    TG2 = (2*X2 - 0.812500000000000000E+00_dp)*0.160000000000000000E+02_dp
                                    patch = 77
                                 ELSE
                                    CALL patch_coordinate(X1, 0.21875_dp, 32.0_dp, TG1, scale)
                                    TG2 = (2*X2 - 0.937500000000000000E+00_dp)*0.160000000000000000E+02_dp
                                    patch = 78
                                 END IF
                              END IF
                           END IF
                        END IF
                     ELSE
                        IF (X2 <= 0.375000000000000000E+00_dp) THEN
                           IF (X1 <= 0.187500000000000000E+00_dp) THEN
                              IF (X2 <= 0.312500000000000000E+00_dp) THEN
                                 CALL patch_coordinate(X1, 0.3125_dp, 16.0_dp, TG1, scale)
                                 TG2 = (2*X2 - 0.562500000000000000E+00_dp)*0.160000000000000000E+02_dp
                                 patch = 79
                              ELSE
                                 IF (X1 <= 0.156250000000000000E+00_dp) THEN
                                    CALL patch_coordinate(X1, 0.28125_dp, 32.0_dp, TG1, scale)
                                    TG2 = (2*X2 - 0.687500000000000000E+00_dp)*0.160000000000000000E+02_dp
                                    patch = 80
                                 ELSE
                                    CALL patch_coordinate(X1, 0.34375_dp, 32.0_dp, TG1, scale)
                                    TG2 = (2*X2 - 0.687500000000000000E+00_dp)*0.160000000000000000E+02_dp
                                    patch = 81
                                 END IF
                              END IF
                           ELSE
                              IF (X2 <= 0.312500000000000000E+00_dp) THEN
                                 CALL patch_coordinate(X1, 0.4375_dp, 16.0_dp, TG1, scale)
                                 TG2 = (2*X2 - 0.562500000000000000E+00_dp)*0.160000000000000000E+02_dp
                                 patch = 82
                              ELSE
                                 CALL patch_coordinate(X1, 0.4375_dp, 16.0_dp, TG1, scale)
                                 TG2 = (2*X2 - 0.687500000000000000E+00_dp)*0.160000000000000000E+02_dp
                                 patch = 83
                              END IF
                           END IF
                        ELSE
                           IF (X1 <= 0.187500000000000000E+00_dp) THEN
                              IF (X2 <= 0.437500000000000000E+00_dp) THEN
                                 CALL patch_coordinate(X1, 0.3125_dp, 16.0_dp, TG1, scale)
                                 TG2 = (2*X2 - 0.812500000000000000E+00_dp)*0.160000000000000000E+02_dp
                                 patch = 84
                              ELSE
                                 IF (X1 <= 0.156250000000000000E+00_dp) THEN
                                    CALL patch_coordinate(X1, 0.28125_dp, 32.0_dp, TG1, scale)
                                    TG2 = (2*X2 - 0.937500000000000000E+00_dp)*0.160000000000000000E+02_dp
                                    patch = 85
                                 ELSE
                                    CALL patch_coordinate(X1, 0.34375_dp, 32.0_dp, TG1, scale)
                                    TG2 = (2*X2 - 0.937500000000000000E+00_dp)*0.160000000000000000E+02_dp
                                    patch = 86
                                 END IF
                              END IF
                           ELSE
                              IF (X1 <= 0.218750000000000000E+00_dp) THEN
                                 CALL patch_coordinate(X1, 0.40625_dp, 32.0_dp, TG1, scale)
                                 TG2 = (2*X2 - 0.875000000000000000E+00_dp)*0.800000000000000000E+01_dp
                                 patch = 87
                              ELSE
                                 CALL patch_coordinate(X1, 0.46875_dp, 32.0_dp, TG1, scale)
                                 TG2 = (2*X2 - 0.875000000000000000E+00_dp)*0.800000000000000000E+01_dp
                                 patch = 88
                              END IF
                           END IF
                        END IF
                     END IF
                  ELSE
                     IF (X1 <= 0.375000000000000000E+00_dp) THEN
                        IF (X2 <= 0.375000000000000000E+00_dp) THEN
                           IF (X1 <= 0.312500000000000000E+00_dp) THEN
                              CALL patch_coordinate(X1, 0.5625_dp, 16.0_dp, TG1, scale)
                              TG2 = (2*X2 - 0.625000000000000000E+00_dp)*0.800000000000000000E+01_dp
                              patch = 89
                           ELSE
                              CALL patch_coordinate(X1, 0.6875_dp, 16.0_dp, TG1, scale)
                              TG2 = (2*X2 - 0.625000000000000000E+00_dp)*0.800000000000000000E+01_dp
                              patch = 90
                           END IF
                        ELSE
                           IF (X1 <= 0.312500000000000000E+00_dp) THEN
                              IF (X1 <= 0.281250000000000000E+00_dp) THEN
                                 CALL patch_coordinate(X1, 0.53125_dp, 32.0_dp, TG1, scale)
                                 TG2 = (2*X2 - 0.875000000000000000E+00_dp)*0.800000000000000000E+01_dp
                                 patch = 91
                              ELSE
                                 CALL patch_coordinate(X1, 0.59375_dp, 32.0_dp, TG1, scale)
                                 TG2 = (2*X2 - 0.875000000000000000E+00_dp)*0.800000000000000000E+01_dp
                                 patch = 92
                              END IF
                           ELSE
                              CALL patch_coordinate(X1, 0.6875_dp, 16.0_dp, TG1, scale)
                              TG2 = (2*X2 - 0.875000000000000000E+00_dp)*0.800000000000000000E+01_dp
                              patch = 93
                           END IF
                        END IF
                     ELSE
                        IF (X1 <= 0.437500000000000000E+00_dp) THEN
                           CALL patch_coordinate(X1, 0.8125_dp, 16.0_dp, TG1, scale)
                           TG2 = (2*X2 - 0.750000000000000000E+00_dp)*0.400000000000000000E+01_dp
                           patch = 94
                        ELSE
                           CALL patch_coordinate(X1, 0.9375_dp, 16.0_dp, TG1, scale)
                           TG2 = (2*X2 - 0.750000000000000000E+00_dp)*0.400000000000000000E+01_dp
                           patch = 95
                        END IF
                     END IF
                  END IF
               END IF
            ELSE
               IF (X1 <= 0.250000000000000000E+00_dp) THEN
                  IF (X1 <= 0.125000000000000000E+00_dp) THEN
                     IF (X1 <= 0.625000000000000000E-01_dp) THEN
                        IF (X1 <= 0.312500000000000000E-01_dp) THEN
                           IF (X1 <= 0.156250000000000000E-01_dp) THEN
                              IF (X1 <= 0.781250000000000000E-02_dp) THEN
                                 IF (X2 <= 0.750000000000000000E+00_dp) THEN
                                    CALL patch_coordinate(X1, 0.0078125_dp, 128.0_dp, TG1, scale)
                                    TG2 = (2*X2 - 0.125000000000000000E+01_dp)*0.400000000000000000E+01_dp
                                    patch = 96
                                 ELSE
                                    CALL patch_coordinate(X1, 0.0078125_dp, 128.0_dp, TG1, scale)
                                    TG2 = (2*X2 - 0.175000000000000000E+01_dp)*0.400000000000000000E+01_dp
                                    patch = 97
                                 END IF
                              ELSE
                                 IF (X2 <= 0.750000000000000000E+00_dp) THEN
                                    CALL patch_coordinate(X1, 0.0234375_dp, 128.0_dp, TG1, scale)
                                    TG2 = (2*X2 - 0.125000000000000000E+01_dp)*0.400000000000000000E+01_dp
                                    patch = 98
                                 ELSE
                                    CALL patch_coordinate(X1, 0.0234375_dp, 128.0_dp, TG1, scale)
                                    TG2 = (2*X2 - 0.175000000000000000E+01_dp)*0.400000000000000000E+01_dp
                                    patch = 99
                                 END IF
                              END IF
                           ELSE
                              IF (X2 <= 0.750000000000000000E+00_dp) THEN
                                 CALL patch_coordinate(X1, 0.046875_dp, 64.0_dp, TG1, scale)
                                 TG2 = (2*X2 - 0.125000000000000000E+01_dp)*0.400000000000000000E+01_dp
                                 patch = 100
                              ELSE
                                 CALL patch_coordinate(X1, 0.046875_dp, 64.0_dp, TG1, scale)
                                 TG2 = (2*X2 - 0.175000000000000000E+01_dp)*0.400000000000000000E+01_dp
                                 patch = 101
                              END IF
                           END IF
                        ELSE
                           IF (X2 <= 0.750000000000000000E+00_dp) THEN
                              IF (X2 <= 0.625000000000000000E+00_dp) THEN
                                 CALL patch_coordinate(X1, 0.09375_dp, 32.0_dp, TG1, scale)
                                 TG2 = (2*X2 - 0.112500000000000000E+01_dp)*0.800000000000000000E+01_dp
                                 patch = 102
                              ELSE
                                 CALL patch_coordinate(X1, 0.09375_dp, 32.0_dp, TG1, scale)
                                 TG2 = (2*X2 - 0.137500000000000000E+01_dp)*0.800000000000000000E+01_dp
                                 patch = 103
                              END IF
                           ELSE
                              IF (X1 <= 0.468750000000000000E-01_dp) THEN
                                 CALL patch_coordinate(X1, 0.078125_dp, 64.0_dp, TG1, scale)
                                 TG2 = (2*X2 - 0.175000000000000000E+01_dp)*0.400000000000000000E+01_dp
                                 patch = 104
                              ELSE
                                 CALL patch_coordinate(X1, 0.109375_dp, 64.0_dp, TG1, scale)
                                 TG2 = (2*X2 - 0.175000000000000000E+01_dp)*0.400000000000000000E+01_dp
                                 patch = 105
                              END IF
                           END IF
                        END IF
                     ELSE
                        IF (X2 <= 0.750000000000000000E+00_dp) THEN
                           IF (X2 <= 0.625000000000000000E+00_dp) THEN
                              IF (X1 <= 0.937500000000000000E-01_dp) THEN
                                 CALL patch_coordinate(X1, 0.15625_dp, 32.0_dp, TG1, scale)
                                 TG2 = (2*X2 - 0.112500000000000000E+01_dp)*0.800000000000000000E+01_dp
                                 patch = 106
                              ELSE
                                 CALL patch_coordinate(X1, 0.21875_dp, 32.0_dp, TG1, scale)
                                 TG2 = (2*X2 - 0.112500000000000000E+01_dp)*0.800000000000000000E+01_dp
                                 patch = 107
                              END IF
                           ELSE
                              IF (X1 <= 0.937500000000000000E-01_dp) THEN
                                 CALL patch_coordinate(X1, 0.15625_dp, 32.0_dp, TG1, scale)
                                 TG2 = (2*X2 - 0.137500000000000000E+01_dp)*0.800000000000000000E+01_dp
                                 patch = 108
                              ELSE
                                 CALL patch_coordinate(X1, 0.21875_dp, 32.0_dp, TG1, scale)
                                 TG2 = (2*X2 - 0.137500000000000000E+01_dp)*0.800000000000000000E+01_dp
                                 patch = 109
                              END IF
                           END IF
                        ELSE
                           IF (X1 <= 0.937500000000000000E-01_dp) THEN
                              CALL patch_coordinate(X1, 0.15625_dp, 32.0_dp, TG1, scale)
                              TG2 = (2*X2 - 0.175000000000000000E+01_dp)*0.400000000000000000E+01_dp
                              patch = 110
                           ELSE
                              CALL patch_coordinate(X1, 0.21875_dp, 32.0_dp, TG1, scale)
                              TG2 = (2*X2 - 0.175000000000000000E+01_dp)*0.400000000000000000E+01_dp
                              patch = 111
                           END IF
                        END IF
                     END IF
                  ELSE
                     IF (X2 <= 0.750000000000000000E+00_dp) THEN
                        IF (X2 <= 0.625000000000000000E+00_dp) THEN
                           IF (X1 <= 0.187500000000000000E+00_dp) THEN
                              IF (X1 <= 0.156250000000000000E+00_dp) THEN
                                 CALL patch_coordinate(X1, 0.28125_dp, 32.0_dp, TG1, scale)
                                 TG2 = (2*X2 - 0.112500000000000000E+01_dp)*0.800000000000000000E+01_dp
                                 patch = 112
                              ELSE
                                 CALL patch_coordinate(X1, 0.34375_dp, 32.0_dp, TG1, scale)
                                 TG2 = (2*X2 - 0.112500000000000000E+01_dp)*0.800000000000000000E+01_dp
                                 patch = 113
                              END IF
                           ELSE
                              IF (X1 <= 0.218750000000000000E+00_dp) THEN
                                 CALL patch_coordinate(X1, 0.40625_dp, 32.0_dp, TG1, scale)
                                 TG2 = (2*X2 - 0.112500000000000000E+01_dp)*0.800000000000000000E+01_dp
                                 patch = 114
                              ELSE
                                 CALL patch_coordinate(X1, 0.46875_dp, 32.0_dp, TG1, scale)
                                 TG2 = (2*X2 - 0.112500000000000000E+01_dp)*0.800000000000000000E+01_dp
                                 patch = 115
                              END IF
                           END IF
                        ELSE
                           IF (X1 <= 0.187500000000000000E+00_dp) THEN
                              IF (X1 <= 0.156250000000000000E+00_dp) THEN
                                 CALL patch_coordinate(X1, 0.28125_dp, 32.0_dp, TG1, scale)
                                 TG2 = (2*X2 - 0.137500000000000000E+01_dp)*0.800000000000000000E+01_dp
                                 patch = 116
                              ELSE
                                 CALL patch_coordinate(X1, 0.34375_dp, 32.0_dp, TG1, scale)
                                 TG2 = (2*X2 - 0.137500000000000000E+01_dp)*0.800000000000000000E+01_dp
                                 patch = 117
                              END IF
                           ELSE
                              IF (X1 <= 0.218750000000000000E+00_dp) THEN
                                 CALL patch_coordinate(X1, 0.40625_dp, 32.0_dp, TG1, scale)
                                 TG2 = (2*X2 - 0.137500000000000000E+01_dp)*0.800000000000000000E+01_dp
                                 patch = 118
                              ELSE
                                 CALL patch_coordinate(X1, 0.46875_dp, 32.0_dp, TG1, scale)
                                 TG2 = (2*X2 - 0.137500000000000000E+01_dp)*0.800000000000000000E+01_dp
                                 patch = 119
                              END IF
                           END IF
                        END IF
                     ELSE
                        IF (X1 <= 0.187500000000000000E+00_dp) THEN
                           IF (X2 <= 0.875000000000000000E+00_dp) THEN
                              CALL patch_coordinate(X1, 0.3125_dp, 16.0_dp, TG1, scale)
                              TG2 = (2*X2 - 0.162500000000000000E+01_dp)*0.800000000000000000E+01_dp
                              patch = 120
                           ELSE
                              CALL patch_coordinate(X1, 0.3125_dp, 16.0_dp, TG1, scale)
                              TG2 = (2*X2 - 0.187500000000000000E+01_dp)*0.800000000000000000E+01_dp
                              patch = 121
                           END IF
                        ELSE
                           IF (X2 <= 0.875000000000000000E+00_dp) THEN
                              IF (X1 <= 0.218750000000000000E+00_dp) THEN
                                 CALL patch_coordinate(X1, 0.40625_dp, 32.0_dp, TG1, scale)
                                 TG2 = (2*X2 - 0.162500000000000000E+01_dp)*0.800000000000000000E+01_dp
                                 patch = 122
                              ELSE
                                 CALL patch_coordinate(X1, 0.46875_dp, 32.0_dp, TG1, scale)
                                 TG2 = (2*X2 - 0.162500000000000000E+01_dp)*0.800000000000000000E+01_dp
                                 patch = 123
                              END IF
                           ELSE
                              CALL patch_coordinate(X1, 0.4375_dp, 16.0_dp, TG1, scale)
                              TG2 = (2*X2 - 0.187500000000000000E+01_dp)*0.800000000000000000E+01_dp
                              patch = 124
                           END IF
                        END IF
                     END IF
                  END IF
               ELSE
                  IF (X1 <= 0.375000000000000000E+00_dp) THEN
                     IF (X2 <= 0.750000000000000000E+00_dp) THEN
                        IF (X1 <= 0.312500000000000000E+00_dp) THEN
                           IF (X2 <= 0.625000000000000000E+00_dp) THEN
                              IF (X1 <= 0.281250000000000000E+00_dp) THEN
                                 CALL patch_coordinate(X1, 0.53125_dp, 32.0_dp, TG1, scale)
                                 TG2 = (2*X2 - 0.112500000000000000E+01_dp)*0.800000000000000000E+01_dp
                                 patch = 125
                              ELSE
                                 CALL patch_coordinate(X1, 0.59375_dp, 32.0_dp, TG1, scale)
                                 TG2 = (2*X2 - 0.112500000000000000E+01_dp)*0.800000000000000000E+01_dp
                                 patch = 126
                              END IF
                           ELSE
                              IF (X1 <= 0.281250000000000000E+00_dp) THEN
                                 CALL patch_coordinate(X1, 0.53125_dp, 32.0_dp, TG1, scale)
                                 TG2 = (2*X2 - 0.137500000000000000E+01_dp)*0.800000000000000000E+01_dp
                                 patch = 127
                              ELSE
                                 CALL patch_coordinate(X1, 0.59375_dp, 32.0_dp, TG1, scale)
                                 TG2 = (2*X2 - 0.137500000000000000E+01_dp)*0.800000000000000000E+01_dp
                                 patch = 128
                              END IF
                           END IF
                        ELSE
                           IF (X2 <= 0.625000000000000000E+00_dp) THEN
                              CALL patch_coordinate(X1, 0.6875_dp, 16.0_dp, TG1, scale)
                              TG2 = (2*X2 - 0.112500000000000000E+01_dp)*0.800000000000000000E+01_dp
                              patch = 129
                           ELSE
                              IF (X1 <= 0.343750000000000000E+00_dp) THEN
                                 CALL patch_coordinate(X1, 0.65625_dp, 32.0_dp, TG1, scale)
                                 TG2 = (2*X2 - 0.137500000000000000E+01_dp)*0.800000000000000000E+01_dp
                                 patch = 130
                              ELSE
                                 CALL patch_coordinate(X1, 0.71875_dp, 32.0_dp, TG1, scale)
                                 TG2 = (2*X2 - 0.137500000000000000E+01_dp)*0.800000000000000000E+01_dp
                                 patch = 131
                              END IF
                           END IF
                        END IF
                     ELSE
                        IF (X1 <= 0.312500000000000000E+00_dp) THEN
                           IF (X2 <= 0.875000000000000000E+00_dp) THEN
                              IF (X1 <= 0.281250000000000000E+00_dp) THEN
                                 CALL patch_coordinate(X1, 0.53125_dp, 32.0_dp, TG1, scale)
                                 TG2 = (2*X2 - 0.162500000000000000E+01_dp)*0.800000000000000000E+01_dp
                                 patch = 132
                              ELSE
                                 CALL patch_coordinate(X1, 0.59375_dp, 32.0_dp, TG1, scale)
                                 TG2 = (2*X2 - 0.162500000000000000E+01_dp)*0.800000000000000000E+01_dp
                                 patch = 133
                              END IF
                           ELSE
                              IF (X1 <= 0.281250000000000000E+00_dp) THEN
                                 CALL patch_coordinate(X1, 0.53125_dp, 32.0_dp, TG1, scale)
                                 TG2 = (2*X2 - 0.187500000000000000E+01_dp)*0.800000000000000000E+01_dp
                                 patch = 134
                              ELSE
                                 CALL patch_coordinate(X1, 0.59375_dp, 32.0_dp, TG1, scale)
                                 TG2 = (2*X2 - 0.187500000000000000E+01_dp)*0.800000000000000000E+01_dp
                                 patch = 135
                              END IF
                           END IF
                        ELSE
                           IF (X2 <= 0.875000000000000000E+00_dp) THEN
                              IF (X1 <= 0.343750000000000000E+00_dp) THEN
                                 CALL patch_coordinate(X1, 0.65625_dp, 32.0_dp, TG1, scale)
                                 TG2 = (2*X2 - 0.162500000000000000E+01_dp)*0.800000000000000000E+01_dp
                                 patch = 136
                              ELSE
                                 CALL patch_coordinate(X1, 0.71875_dp, 32.0_dp, TG1, scale)
                                 TG2 = (2*X2 - 0.162500000000000000E+01_dp)*0.800000000000000000E+01_dp
                                 patch = 137
                              END IF
                           ELSE
                              CALL patch_coordinate(X1, 0.6875_dp, 16.0_dp, TG1, scale)
                              TG2 = (2*X2 - 0.187500000000000000E+01_dp)*0.800000000000000000E+01_dp
                              patch = 138
                           END IF
                        END IF
                     END IF
                  ELSE
                     IF (X2 <= 0.750000000000000000E+00_dp) THEN
                        IF (X1 <= 0.437500000000000000E+00_dp) THEN
                           IF (X2 <= 0.625000000000000000E+00_dp) THEN
                              CALL patch_coordinate(X1, 0.8125_dp, 16.0_dp, TG1, scale)
                              TG2 = (2*X2 - 0.112500000000000000E+01_dp)*0.800000000000000000E+01_dp
                              patch = 139
                           ELSE
                              CALL patch_coordinate(X1, 0.8125_dp, 16.0_dp, TG1, scale)
                              TG2 = (2*X2 - 0.137500000000000000E+01_dp)*0.800000000000000000E+01_dp
                              patch = 140
                           END IF
                        ELSE
                           CALL patch_coordinate(X1, 0.9375_dp, 16.0_dp, TG1, scale)
                           TG2 = (2*X2 - 0.125000000000000000E+01_dp)*0.400000000000000000E+01_dp
                           patch = 141
                        END IF
                     ELSE
                        IF (X1 <= 0.437500000000000000E+00_dp) THEN
                           IF (X2 <= 0.875000000000000000E+00_dp) THEN
                              CALL patch_coordinate(X1, 0.8125_dp, 16.0_dp, TG1, scale)
                              TG2 = (2*X2 - 0.162500000000000000E+01_dp)*0.800000000000000000E+01_dp
                              patch = 142
                           ELSE
                              CALL patch_coordinate(X1, 0.8125_dp, 16.0_dp, TG1, scale)
                              TG2 = (2*X2 - 0.187500000000000000E+01_dp)*0.800000000000000000E+01_dp
                              patch = 143
                           END IF
                        ELSE
                           IF (X2 <= 0.875000000000000000E+00_dp) THEN
                              CALL patch_coordinate(X1, 0.9375_dp, 16.0_dp, TG1, scale)
                              TG2 = (2*X2 - 0.162500000000000000E+01_dp)*0.800000000000000000E+01_dp
                              patch = 144
                           ELSE
                              CALL patch_coordinate(X1, 0.9375_dp, 16.0_dp, TG1, scale)
                              TG2 = (2*X2 - 0.187500000000000000E+01_dp)*0.800000000000000000E+01_dp
                              patch = 145
                           END IF
                        END IF
                     END IF
                  END IF
               END IF
            END IF
         ELSE
            IF (X1 <= 0.750000000000000000E+00_dp) THEN
               IF (X2 <= 0.500000000000000000E+00_dp) THEN
                  IF (X1 <= 0.625000000000000000E+00_dp) THEN
                     IF (X2 <= 0.250000000000000000E+00_dp) THEN
                        CALL patch_coordinate(X1, 1.125_dp, 8.0_dp, TG1, scale)
                        TG2 = (2*X2 - 0.250000000000000000E+00_dp)*0.400000000000000000E+01_dp
                        patch = 146
                     ELSE
                        CALL patch_coordinate(X1, 1.125_dp, 8.0_dp, TG1, scale)
                        TG2 = (2*X2 - 0.750000000000000000E+00_dp)*0.400000000000000000E+01_dp
                        patch = 147
                     END IF
                  ELSE
                     CALL patch_coordinate(X1, 1.375_dp, 8.0_dp, TG1, scale)
                     TG2 = (2*X2 - 0.500000000000000000E+00_dp)*0.200000000000000000E+01_dp
                     patch = 148
                  END IF
               ELSE
                  IF (X1 <= 0.625000000000000000E+00_dp) THEN
                     IF (X2 <= 0.750000000000000000E+00_dp) THEN
                        IF (X1 <= 0.562500000000000000E+00_dp) THEN
                           CALL patch_coordinate(X1, 1.0625_dp, 16.0_dp, TG1, scale)
                           TG2 = (2*X2 - 0.125000000000000000E+01_dp)*0.400000000000000000E+01_dp
                           patch = 149
                        ELSE
                           CALL patch_coordinate(X1, 1.1875_dp, 16.0_dp, TG1, scale)
                           TG2 = (2*X2 - 0.125000000000000000E+01_dp)*0.400000000000000000E+01_dp
                           patch = 150
                        END IF
                     ELSE
                        IF (X1 <= 0.562500000000000000E+00_dp) THEN
                           CALL patch_coordinate(X1, 1.0625_dp, 16.0_dp, TG1, scale)
                           TG2 = (2*X2 - 0.175000000000000000E+01_dp)*0.400000000000000000E+01_dp
                           patch = 151
                        ELSE
                           CALL patch_coordinate(X1, 1.1875_dp, 16.0_dp, TG1, scale)
                           TG2 = (2*X2 - 0.175000000000000000E+01_dp)*0.400000000000000000E+01_dp
                           patch = 152
                        END IF
                     END IF
                  ELSE
                     IF (X1 <= 0.687500000000000000E+00_dp) THEN
                        CALL patch_coordinate(X1, 1.3125_dp, 16.0_dp, TG1, scale)
                        TG2 = (2*X2 - 0.150000000000000000E+01_dp)*0.200000000000000000E+01_dp
                        patch = 153
                     ELSE
                        CALL patch_coordinate(X1, 1.4375_dp, 16.0_dp, TG1, scale)
                        TG2 = (2*X2 - 0.150000000000000000E+01_dp)*0.200000000000000000E+01_dp
                        patch = 154
                     END IF
                  END IF
               END IF
            ELSE
               CALL patch_coordinate(X1, 1.75_dp, 4.0_dp, TG1, scale)
               TG2 = (2*X2 - 0.100000000000000000E+01_dp)*0.100000000000000000E+01_dp
               patch = 155
            END IF
         END IF
      ELSE
         IF (T < lower) THEN
            use_gamma = .TRUE.
            RETURN
         END IF
         X2 = 11.0_dp/R
         X1 = (T - lower)/(upper - lower)
         IF (X1 <= 0.500000000000000000E+00_dp) THEN
            IF (X1 <= 0.250000000000000000E+00_dp) THEN
               IF (X2 <= 0.500000000000000000E+00_dp) THEN
                  IF (X1 <= 0.125000000000000000E+00_dp) THEN
                     IF (X2 <= 0.250000000000000000E+00_dp) THEN
                        CALL patch_coordinate(X1, 0.125_dp, 8.0_dp, TG1, scale)
                        TG2 = (2*X2 - 0.250000000000000000E+00_dp)*0.400000000000000000E+01_dp
                        patch = 156
                     ELSE
                        CALL patch_coordinate(X1, 0.125_dp, 8.0_dp, TG1, scale)
                        TG2 = (2*X2 - 0.750000000000000000E+00_dp)*0.400000000000000000E+01_dp
                        patch = 157
                     END IF
                  ELSE
                     CALL patch_coordinate(X1, 0.375_dp, 8.0_dp, TG1, scale)
                     TG2 = (2*X2 - 0.500000000000000000E+00_dp)*0.200000000000000000E+01_dp
                     patch = 158
                  END IF
               ELSE
                  IF (X1 <= 0.125000000000000000E+00_dp) THEN
                     IF (X2 <= 0.750000000000000000E+00_dp) THEN
                        IF (X2 <= 0.625000000000000000E+00_dp) THEN
                           CALL patch_coordinate(X1, 0.125_dp, 8.0_dp, TG1, scale)
                           TG2 = (2*X2 - 0.112500000000000000E+01_dp)*0.800000000000000000E+01_dp
                           patch = 159
                        ELSE
                           IF (X1 <= 0.625000000000000000E-01_dp) THEN
                              CALL patch_coordinate(X1, 0.0625_dp, 16.0_dp, TG1, scale)
                              TG2 = (2*X2 - 0.137500000000000000E+01_dp)*0.800000000000000000E+01_dp
                              patch = 160
                           ELSE
                              CALL patch_coordinate(X1, 0.1875_dp, 16.0_dp, TG1, scale)
                              TG2 = (2*X2 - 0.137500000000000000E+01_dp)*0.800000000000000000E+01_dp
                              patch = 161
                           END IF
                        END IF
                     ELSE
                        IF (X1 <= 0.625000000000000000E-01_dp) THEN
                           IF (X2 <= 0.875000000000000000E+00_dp) THEN
                              IF (X1 <= 0.312500000000000000E-01_dp) THEN
                                 IF (X2 <= 0.812500000000000000E+00_dp) THEN
                                    CALL patch_coordinate(X1, 0.03125_dp, 32.0_dp, TG1, scale)
                                    TG2 = (2*X2 - 0.156250000000000000E+01_dp)*0.160000000000000000E+02_dp
                                    patch = 162
                                 ELSE
                                    CALL patch_coordinate(X1, 0.03125_dp, 32.0_dp, TG1, scale)
                                    TG2 = (2*X2 - 0.168750000000000000E+01_dp)*0.160000000000000000E+02_dp
                                    patch = 163
                                 END IF
                              ELSE
                                 CALL patch_coordinate(X1, 0.09375_dp, 32.0_dp, TG1, scale)
                                 TG2 = (2*X2 - 0.162500000000000000E+01_dp)*0.800000000000000000E+01_dp
                                 patch = 164
                              END IF
                           ELSE
                              IF (X1 <= 0.312500000000000000E-01_dp) THEN
                                 IF (X2 <= 0.937500000000000000E+00_dp) THEN
                                    IF (X1 <= 0.156250000000000000E-01_dp) THEN
                                       IF (X2 <= 0.906250000000000000E+00_dp) THEN
                                          CALL patch_coordinate(X1, 0.015625_dp, 64.0_dp, TG1, scale)
                                          TG2 = (2*X2 - 0.178125000000000000E+01_dp)*0.320000000000000000E+02_dp
                                          patch = 165
                                       ELSE
                                          CALL patch_coordinate(X1, 0.015625_dp, 64.0_dp, TG1, scale)
                                          TG2 = (2*X2 - 0.184375000000000000E+01_dp)*0.320000000000000000E+02_dp
                                          patch = 166
                                       END IF
                                    ELSE
                                       CALL patch_coordinate(X1, 0.046875_dp, 64.0_dp, TG1, scale)
                                       TG2 = (2*X2 - 0.181250000000000000E+01_dp)*0.160000000000000000E+02_dp
                                       patch = 167
                                    END IF
                                 ELSE
                                    IF (X1 <= 0.156250000000000000E-01_dp) THEN
                                       IF (X2 <= 0.968750000000000000E+00_dp) THEN
                                          IF (X1 <= 0.781250000000000000E-02_dp) THEN
                                             IF (X2 <= 0.953125000000000000E+00_dp) THEN
                                                CALL patch_coordinate(X1, 0.0078125_dp, 128.0_dp, TG1, scale)
                                                TG2 = (2*X2 - 0.189062500000000000E+01_dp)*0.640000000000000000E+02_dp
                                                patch = 168
                                             ELSE
                                                CALL patch_coordinate(X1, 0.0078125_dp, 128.0_dp, TG1, scale)
                                                TG2 = (2*X2 - 0.192187500000000000E+01_dp)*0.640000000000000000E+02_dp
                                                patch = 169
                                             END IF
                                          ELSE
                                             CALL patch_coordinate(X1, 0.0234375_dp, 128.0_dp, TG1, scale)
                                             TG2 = (2*X2 - 0.190625000000000000E+01_dp)*0.320000000000000000E+02_dp
                                             patch = 170
                                          END IF
                                       ELSE
                                          IF (X1 <= 0.781250000000000000E-02_dp) THEN
                                             IF (X2 <= 0.984375000000000000E+00_dp) THEN
                                                CALL patch_coordinate(X1, 0.0078125_dp, 128.0_dp, TG1, scale)
                                                TG2 = (2*X2 - 0.195312500000000000E+01_dp)*0.640000000000000000E+02_dp
                                                patch = 171
                                             ELSE
                                                CALL patch_coordinate(X1, 0.0078125_dp, 128.0_dp, TG1, scale)
                                                TG2 = (2*X2 - 0.198437500000000000E+01_dp)*0.640000000000000000E+02_dp
                                                patch = 172
                                             END IF
                                          ELSE
                                             IF (X2 <= 0.984375000000000000E+00_dp) THEN
                                                CALL patch_coordinate(X1, 0.0234375_dp, 128.0_dp, TG1, scale)
                                                TG2 = (2*X2 - 0.195312500000000000E+01_dp)*0.640000000000000000E+02_dp
                                                patch = 173
                                             ELSE
                                                CALL patch_coordinate(X1, 0.0234375_dp, 128.0_dp, TG1, scale)
                                                TG2 = (2*X2 - 0.198437500000000000E+01_dp)*0.640000000000000000E+02_dp
                                                patch = 174
                                             END IF
                                          END IF
                                       END IF
                                    ELSE
                                       IF (X2 <= 0.968750000000000000E+00_dp) THEN
                                          CALL patch_coordinate(X1, 0.046875_dp, 64.0_dp, TG1, scale)
                                          TG2 = (2*X2 - 0.190625000000000000E+01_dp)*0.320000000000000000E+02_dp
                                          patch = 175
                                       ELSE
                                          IF (X1 <= 0.234375000000000000E-01_dp) THEN
                                             CALL patch_coordinate(X1, 0.0390625_dp, 128.0_dp, TG1, scale)
                                             TG2 = (2*X2 - 0.196875000000000000E+01_dp)*0.320000000000000000E+02_dp
                                             patch = 176
                                          ELSE
                                             CALL patch_coordinate(X1, 0.0546875_dp, 128.0_dp, TG1, scale)
                                             TG2 = (2*X2 - 0.196875000000000000E+01_dp)*0.320000000000000000E+02_dp
                                             patch = 177
                                          END IF
                                       END IF
                                    END IF
                                 END IF
                              ELSE
                                 IF (X2 <= 0.937500000000000000E+00_dp) THEN
                                    CALL patch_coordinate(X1, 0.09375_dp, 32.0_dp, TG1, scale)
                                    TG2 = (2*X2 - 0.181250000000000000E+01_dp)*0.160000000000000000E+02_dp
                                    patch = 178
                                 ELSE
                                    IF (X1 <= 0.468750000000000000E-01_dp) THEN
                                       IF (X2 <= 0.968750000000000000E+00_dp) THEN
                                          CALL patch_coordinate(X1, 0.078125_dp, 64.0_dp, TG1, scale)
                                          TG2 = (2*X2 - 0.190625000000000000E+01_dp)*0.320000000000000000E+02_dp
                                          patch = 179
                                       ELSE
                                          CALL patch_coordinate(X1, 0.078125_dp, 64.0_dp, TG1, scale)
                                          TG2 = (2*X2 - 0.196875000000000000E+01_dp)*0.320000000000000000E+02_dp
                                          patch = 180
                                       END IF
                                    ELSE
                                       CALL patch_coordinate(X1, 0.109375_dp, 64.0_dp, TG1, scale)
                                       TG2 = (2*X2 - 0.193750000000000000E+01_dp)*0.160000000000000000E+02_dp
                                       patch = 181
                                    END IF
                                 END IF
                              END IF
                           END IF
                        ELSE
                           IF (X2 <= 0.875000000000000000E+00_dp) THEN
                              CALL patch_coordinate(X1, 0.1875_dp, 16.0_dp, TG1, scale)
                              TG2 = (2*X2 - 0.162500000000000000E+01_dp)*0.800000000000000000E+01_dp
                              patch = 182
                           ELSE
                              IF (X1 <= 0.937500000000000000E-01_dp) THEN
                                 IF (X2 <= 0.937500000000000000E+00_dp) THEN
                                    CALL patch_coordinate(X1, 0.15625_dp, 32.0_dp, TG1, scale)
                                    TG2 = (2*X2 - 0.181250000000000000E+01_dp)*0.160000000000000000E+02_dp
                                    patch = 183
                                 ELSE
                                    CALL patch_coordinate(X1, 0.15625_dp, 32.0_dp, TG1, scale)
                                    TG2 = (2*X2 - 0.193750000000000000E+01_dp)*0.160000000000000000E+02_dp
                                    patch = 184
                                 END IF
                              ELSE
                                 CALL patch_coordinate(X1, 0.21875_dp, 32.0_dp, TG1, scale)
                                 TG2 = (2*X2 - 0.187500000000000000E+01_dp)*0.800000000000000000E+01_dp
                                 patch = 185
                              END IF
                           END IF
                        END IF
                     END IF
                  ELSE
                     IF (X2 <= 0.750000000000000000E+00_dp) THEN
                        CALL patch_coordinate(X1, 0.375_dp, 8.0_dp, TG1, scale)
                        TG2 = (2*X2 - 0.125000000000000000E+01_dp)*0.400000000000000000E+01_dp
                        patch = 186
                     ELSE
                        IF (X1 <= 0.187500000000000000E+00_dp) THEN
                           IF (X2 <= 0.875000000000000000E+00_dp) THEN
                              CALL patch_coordinate(X1, 0.3125_dp, 16.0_dp, TG1, scale)
                              TG2 = (2*X2 - 0.162500000000000000E+01_dp)*0.800000000000000000E+01_dp
                              patch = 187
                           ELSE
                              CALL patch_coordinate(X1, 0.3125_dp, 16.0_dp, TG1, scale)
                              TG2 = (2*X2 - 0.187500000000000000E+01_dp)*0.800000000000000000E+01_dp
                              patch = 188
                           END IF
                        ELSE
                           CALL patch_coordinate(X1, 0.4375_dp, 16.0_dp, TG1, scale)
                           TG2 = (2*X2 - 0.175000000000000000E+01_dp)*0.400000000000000000E+01_dp
                           patch = 189
                        END IF
                     END IF
                  END IF
               END IF
            ELSE
               IF (X1 <= 0.375000000000000000E+00_dp) THEN
                  IF (X1 <= 0.312500000000000000E+00_dp) THEN
                     IF (X1 <= 0.281250000000000000E+00_dp) THEN
                        CALL patch_coordinate(X1, 0.53125_dp, 32.0_dp, TG1, scale)
                        TG2 = (2*X2 - 0.100000000000000000E+01_dp)*0.100000000000000000E+01_dp
                        patch = 190
                     ELSE
                        CALL patch_coordinate(X1, 0.59375_dp, 32.0_dp, TG1, scale)
                        TG2 = (2*X2 - 0.100000000000000000E+01_dp)*0.100000000000000000E+01_dp
                        patch = 191
                     END IF
                  ELSE
                     IF (X2 <= 0.500000000000000000E+00_dp) THEN
                        CALL patch_coordinate(X1, 0.6875_dp, 16.0_dp, TG1, scale)
                        TG2 = (2*X2 - 0.500000000000000000E+00_dp)*0.200000000000000000E+01_dp
                        patch = 192
                     ELSE
                        CALL patch_coordinate(X1, 0.6875_dp, 16.0_dp, TG1, scale)
                        TG2 = (2*X2 - 0.150000000000000000E+01_dp)*0.200000000000000000E+01_dp
                        patch = 193
                     END IF
                  END IF
               ELSE
                  IF (X1 <= 0.437500000000000000E+00_dp) THEN
                     IF (X2 <= 0.500000000000000000E+00_dp) THEN
                        CALL patch_coordinate(X1, 0.8125_dp, 16.0_dp, TG1, scale)
                        TG2 = (2*X2 - 0.500000000000000000E+00_dp)*0.200000000000000000E+01_dp
                        patch = 194
                     ELSE
                        IF (X1 <= 0.406250000000000000E+00_dp) THEN
                           CALL patch_coordinate(X1, 0.78125_dp, 32.0_dp, TG1, scale)
                           TG2 = (2*X2 - 0.150000000000000000E+01_dp)*0.200000000000000000E+01_dp
                           patch = 195
                        ELSE
                           CALL patch_coordinate(X1, 0.84375_dp, 32.0_dp, TG1, scale)
                           TG2 = (2*X2 - 0.150000000000000000E+01_dp)*0.200000000000000000E+01_dp
                           patch = 196
                        END IF
                     END IF
                  ELSE
                     IF (X2 <= 0.500000000000000000E+00_dp) THEN
                        CALL patch_coordinate(X1, 0.9375_dp, 16.0_dp, TG1, scale)
                        TG2 = (2*X2 - 0.500000000000000000E+00_dp)*0.200000000000000000E+01_dp
                        patch = 197
                     ELSE
                        CALL patch_coordinate(X1, 0.9375_dp, 16.0_dp, TG1, scale)
                        TG2 = (2*X2 - 0.150000000000000000E+01_dp)*0.200000000000000000E+01_dp
                        patch = 198
                     END IF
                  END IF
               END IF
            END IF
         ELSE
            IF (X1 <= 0.750000000000000000E+00_dp) THEN
               IF (X1 <= 0.625000000000000000E+00_dp) THEN
                  IF (X2 <= 0.500000000000000000E+00_dp) THEN
                     IF (X1 <= 0.562500000000000000E+00_dp) THEN
                        CALL patch_coordinate(X1, 1.0625_dp, 16.0_dp, TG1, scale)
                        TG2 = (2*X2 - 0.500000000000000000E+00_dp)*0.200000000000000000E+01_dp
                        patch = 199
                     ELSE
                        CALL patch_coordinate(X1, 1.1875_dp, 16.0_dp, TG1, scale)
                        TG2 = (2*X2 - 0.500000000000000000E+00_dp)*0.200000000000000000E+01_dp
                        patch = 200
                     END IF
                  ELSE
                     IF (X1 <= 0.562500000000000000E+00_dp) THEN
                        CALL patch_coordinate(X1, 1.0625_dp, 16.0_dp, TG1, scale)
                        TG2 = (2*X2 - 0.150000000000000000E+01_dp)*0.200000000000000000E+01_dp
                        patch = 201
                     ELSE
                        CALL patch_coordinate(X1, 1.1875_dp, 16.0_dp, TG1, scale)
                        TG2 = (2*X2 - 0.150000000000000000E+01_dp)*0.200000000000000000E+01_dp
                        patch = 202
                     END IF
                  END IF
               ELSE
                  IF (X2 <= 0.500000000000000000E+00_dp) THEN
                     IF (X1 <= 0.687500000000000000E+00_dp) THEN
                        CALL patch_coordinate(X1, 1.3125_dp, 16.0_dp, TG1, scale)
                        TG2 = (2*X2 - 0.500000000000000000E+00_dp)*0.200000000000000000E+01_dp
                        patch = 203
                     ELSE
                        CALL patch_coordinate(X1, 1.4375_dp, 16.0_dp, TG1, scale)
                        TG2 = (2*X2 - 0.500000000000000000E+00_dp)*0.200000000000000000E+01_dp
                        patch = 204
                     END IF
                  ELSE
                     CALL patch_coordinate(X1, 1.375_dp, 8.0_dp, TG1, scale)
                     TG2 = (2*X2 - 0.150000000000000000E+01_dp)*0.200000000000000000E+01_dp
                     patch = 205
                  END IF
               END IF
            ELSE
               IF (X1 <= 0.875000000000000000E+00_dp) THEN
                  CALL patch_coordinate(X1, 1.625_dp, 8.0_dp, TG1, scale)
                  TG2 = (2*X2 - 0.100000000000000000E+01_dp)*0.100000000000000000E+01_dp
                  patch = 206
               ELSE
                  CALL patch_coordinate(X1, 1.875_dp, 8.0_dp, TG1, scale)
                  TG2 = (2*X2 - 0.100000000000000000E+01_dp)*0.100000000000000000E+01_dp
                  patch = 207
               END IF
            END IF
         END IF
      END IF
      CALL PD2VAL(RES, NDERIV, TG1, TG2, C0(1, patch), dres_dt, scale/(upper - lower))
   END SUBROUTINE evaluate_moments

! **************************************************************************************************
!> \brief ...
!> \param Nder the number of derivatives that will actually be used
!> \param iunit contains the data file to initialize the table
!> \param mepos ...
!> \param group ...
! **************************************************************************************************
   SUBROUTINE init(Nder, iunit, mepos, group)
      INTEGER, INTENT(IN)                                :: Nder, iunit, mepos

      CLASS(mp_comm_type), INTENT(IN)                     :: group

      INTEGER                                            :: I
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: chunk

      patches = 207
      IF (Nder > nderiv_max) CPABORT("T_C_G0 init failed")
      nderiv_init = Nder
      IF (ALLOCATED(C0)) DEALLOCATE (C0)
      ! round up to a multiple of 32 to give some generous alignment for each C0
      ALLOCATE (C0(32*((31 + (Nder + 1)*(degree + 1)*(degree + 2)/2)/32), patches))
      ! init mpi'ed buffers to silence warnings under valgrind
      C0 = 1.0E99_dp
      IF (mepos == 0) THEN
         ALLOCATE (chunk((nderiv_max + 1)*(degree + 1)*(degree + 2)/2))
         DO I = 1, patches
            READ (iunit, *) chunk
            C0(1:(Nder + 1)*(degree + 1)*(degree + 2)/2, I) = chunk(1:(Nder + 1)*(degree + 1)*(degree + 2)/2)
         END DO
         DEALLOCATE (chunk)
      END IF
      CALL group%bcast(C0, 0)

   END SUBROUTINE init

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE free_C0()
      IF (ALLOCATED(C0)) DEALLOCATE (C0)
      nderiv_init = -1
   END SUBROUTINE free_C0

! **************************************************************************************************
!> \brief ...
!> \param RES ...
!> \param NDERIV ...
!> \param TG1 ...
!> \param TG2 ...
!> \param C0 ...
!> \param dres_dt ...
!> \param scale ...
! **************************************************************************************************
   SUBROUTINE PD2VAL(RES, NDERIV, TG1, TG2, C0, dres_dt, scale)
      REAL(KIND=dp), INTENT(OUT)                         :: res(*)
      INTEGER, INTENT(IN)                                :: NDERIV
      REAL(KIND=dp), INTENT(IN)                          :: TG1, TG2, C0(105, *)
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: dres_dt(*)
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: scale

      REAL(KIND=dp), PARAMETER :: SQRT2 = 1.4142135623730950488016887242096980785696718753_dp

      INTEGER                                            :: j, K
      REAL(KIND=dp)                                      :: dT1(0:13), T1(0:13), T2(0:13)

      T1(0) = 1.0_dp
      T2(0) = 1.0_dp
      T1(1) = SQRT2*TG1
      T2(1) = SQRT2*TG2
      T1(2) = 2*TG1*T1(1) - SQRT2
      T2(2) = 2*TG2*T2(1) - SQRT2
      T1(3) = 2*TG1*T1(2) - T1(1)
      T2(3) = 2*TG2*T2(2) - T2(1)
      T1(4) = 2*TG1*T1(3) - T1(2)
      T2(4) = 2*TG2*T2(3) - T2(2)
      T1(5) = 2*TG1*T1(4) - T1(3)
      T2(5) = 2*TG2*T2(4) - T2(3)
      T1(6) = 2*TG1*T1(5) - T1(4)
      T2(6) = 2*TG2*T2(5) - T2(4)
      T1(7) = 2*TG1*T1(6) - T1(5)
      T2(7) = 2*TG2*T2(6) - T2(5)
      T1(8) = 2*TG1*T1(7) - T1(6)
      T2(8) = 2*TG2*T2(7) - T2(6)
      T1(9) = 2*TG1*T1(8) - T1(7)
      T2(9) = 2*TG2*T2(8) - T2(7)
      T1(10) = 2*TG1*T1(9) - T1(8)
      T2(10) = 2*TG2*T2(9) - T2(8)
      T1(11) = 2*TG1*T1(10) - T1(9)
      T2(11) = 2*TG2*T2(10) - T2(9)
      T1(12) = 2*TG1*T1(11) - T1(10)
      T2(12) = 2*TG2*T2(11) - T2(10)
      T1(13) = 2*TG1*T1(12) - T1(11)
      T2(13) = 2*TG2*T2(12) - T2(11)
      DO K = 1, NDERIV + 1
         RES(K) = 0.0_dp
         RES(K) = RES(K) + DOT_PRODUCT(T1(0:13), C0(1:14, K))*T2(0)
         RES(K) = RES(K) + DOT_PRODUCT(T1(0:12), C0(15:27, K))*T2(1)
         RES(K) = RES(K) + DOT_PRODUCT(T1(0:11), C0(28:39, K))*T2(2)
         RES(K) = RES(K) + DOT_PRODUCT(T1(0:10), C0(40:50, K))*T2(3)
         RES(K) = RES(K) + DOT_PRODUCT(T1(0:9), C0(51:60, K))*T2(4)
         RES(K) = RES(K) + DOT_PRODUCT(T1(0:8), C0(61:69, K))*T2(5)
         RES(K) = RES(K) + DOT_PRODUCT(T1(0:7), C0(70:77, K))*T2(6)
         RES(K) = RES(K) + DOT_PRODUCT(T1(0:6), C0(78:84, K))*T2(7)
         RES(K) = RES(K) + DOT_PRODUCT(T1(0:5), C0(85:90, K))*T2(8)
         RES(K) = RES(K) + DOT_PRODUCT(T1(0:4), C0(91:95, K))*T2(9)
         RES(K) = RES(K) + DOT_PRODUCT(T1(0:3), C0(96:99, K))*T2(10)
         RES(K) = RES(K) + DOT_PRODUCT(T1(0:2), C0(100:102, K))*T2(11)
         RES(K) = RES(K) + DOT_PRODUCT(T1(0:1), C0(103:104, K))*T2(12)
         RES(K) = RES(K) + DOT_PRODUCT(T1(0:0), C0(105:105, K))*T2(13)
      END DO
      IF (PRESENT(dres_dt)) THEN
         CPASSERT(PRESENT(scale))
         dT1(0) = 0._dp
         dT1(1) = SQRT2
         dT1(2) = 2*T1(1) + 2*TG1*dT1(1)
         DO j = 3, 13
            dT1(j) = 2*T1(j - 1) + 2*TG1*dT1(j - 1) - dT1(j - 2)
         END DO
         DO K = 1, NDERIV + 1
            dres_dt(K) = 0.0_dp
            dres_dt(K) = dres_dt(K) + DOT_PRODUCT(dT1(0:13), C0(1:14, K))*T2(0)
            dres_dt(K) = dres_dt(K) + DOT_PRODUCT(dT1(0:12), C0(15:27, K))*T2(1)
            dres_dt(K) = dres_dt(K) + DOT_PRODUCT(dT1(0:11), C0(28:39, K))*T2(2)
            dres_dt(K) = dres_dt(K) + DOT_PRODUCT(dT1(0:10), C0(40:50, K))*T2(3)
            dres_dt(K) = dres_dt(K) + DOT_PRODUCT(dT1(0:9), C0(51:60, K))*T2(4)
            dres_dt(K) = dres_dt(K) + DOT_PRODUCT(dT1(0:8), C0(61:69, K))*T2(5)
            dres_dt(K) = dres_dt(K) + DOT_PRODUCT(dT1(0:7), C0(70:77, K))*T2(6)
            dres_dt(K) = dres_dt(K) + DOT_PRODUCT(dT1(0:6), C0(78:84, K))*T2(7)
            dres_dt(K) = dres_dt(K) + DOT_PRODUCT(dT1(0:5), C0(85:90, K))*T2(8)
            dres_dt(K) = dres_dt(K) + DOT_PRODUCT(dT1(0:4), C0(91:95, K))*T2(9)
            dres_dt(K) = dres_dt(K) + DOT_PRODUCT(dT1(0:3), C0(96:99, K))*T2(10)
            dres_dt(K) = dres_dt(K) + DOT_PRODUCT(dT1(0:2), C0(100:102, K))*T2(11)
            dres_dt(K) = dres_dt(K) + DOT_PRODUCT(dT1(0:1), C0(103:104, K))*T2(12)
            dres_dt(K) = dres_dt(K) + DOT_PRODUCT(dT1(0:0), C0(105:105, K))*T2(13)
            dres_dt(K) = dres_dt(K)*scale
         END DO
      END IF
   END SUBROUTINE PD2VAL

! **************************************************************************************************
!> \brief Returns the value of nderiv_init so that one can check if opening the potential file is
!>        worhtwhile
!> \return ...
!> \author A. Bussy, 05.2019
! **************************************************************************************************
   FUNCTION get_lmax_init() RESULT(lmax_init)

      INTEGER                                            :: lmax_init

      lmax_init = nderiv_init

   END FUNCTION get_lmax_init

END MODULE t_c_g0
