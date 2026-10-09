!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: LGPL-3.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
#if defined(__TBLITE_GXTB)
! CP2K-owned transform contract; derived from the integration callback API.
! SPDX-Identifier: LGPL-3.0-or-later

!> Internal CP2K Fourier transform interface; no provider callback API required.
MODULE tblite_transform_contract
   USE mctc_env,                        ONLY: wp

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: exchange_transform_type
   PUBLIC :: exchange_k_to_r, exchange_r_to_k, exchange_r_to_k_adjoint, exchange_k_to_r_adjoint

   INTEGER, PARAMETER :: exchange_k_to_r = 1, exchange_r_to_k = 2
   INTEGER, PARAMETER :: exchange_r_to_k_adjoint = 3, exchange_k_to_r_adjoint = 4

   !> The caller certifies the complete Fourier grid and owns its maps and scratch.
   !> No concurrent calls may share a mutable instance. Without this optional object,
   !> the existing explicit phase-sum implementation remains unchanged.
   TYPE, abstract :: exchange_transform_type
   CONTAINS
      PROCEDURE(apply_transform), deferred :: apply
   END TYPE exchange_transform_type

   abstract INTERFACE
      !> Apply one map to flattened complex AO-pair rows and grid columns.
      !> Callers pass their canonical matrix views, including for strided spin slices.
      !> k_to_r: X_R = sum_k w_k exp(-i k.R) X_k.
      !> r_to_k: X_k = sum_R exp(+i k.R) X_R.
      !> r_to_k_adjoint: X_R = sum_k exp(-i k.R) X_k (no weights).
      !> k_to_r_adjoint: X_k = w_k sum_R exp(+i k.R) X_R.
      !> Adjoint means the ordinary complex Euclidean inner product, not a BZ-weighted one.
! **************************************************************************************************
!> \brief ...
!> \param self ...
!> \param source ...
!> \param TARGET ...
!> \param operation ...
! **************************************************************************************************
      SUBROUTINE apply_transform(self, source, TARGET, operation)
         import :: exchange_transform_type, wp
         class(exchange_transform_type), INTENT(inout) :: self
         COMPLEX(wp), INTENT(in) :: source(:, :)
         COMPLEX(wp), INTENT(out) :: target(:, :)
         INTEGER, INTENT(in) :: operation
      END SUBROUTINE apply_transform
   END INTERFACE
END MODULE tblite_transform_contract

#endif
