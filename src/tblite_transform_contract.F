#if defined(__TBLITE_GXTB)
! CP2K-owned transform contract; derived from the integration callback API.
! SPDX-Identifier: LGPL-3.0-or-later

!> Internal CP2K Fourier transform interface; no provider callback API required.
module tblite_transform_contract
   use mctc_env, only : wp
   implicit none
   private
   public :: exchange_transform_type
   public :: exchange_k_to_r, exchange_r_to_k, exchange_r_to_k_adjoint, exchange_k_to_r_adjoint

   integer, parameter :: exchange_k_to_r = 1, exchange_r_to_k = 2
   integer, parameter :: exchange_r_to_k_adjoint = 3, exchange_k_to_r_adjoint = 4

   !> The caller certifies the complete Fourier grid and owns its maps and scratch.
   !> No concurrent calls may share a mutable instance. Without this optional object,
   !> the existing explicit phase-sum implementation remains unchanged.
   type, abstract :: exchange_transform_type
   contains
      procedure(apply_transform), deferred :: apply
   end type exchange_transform_type

   abstract interface
      !> Apply one map to flattened complex AO-pair rows and grid columns.
      !> Callers pass their canonical matrix views, including for strided spin slices.
      !> k_to_r: X_R = sum_k w_k exp(-i k.R) X_k.
      !> r_to_k: X_k = sum_R exp(+i k.R) X_R.
      !> r_to_k_adjoint: X_R = sum_k exp(-i k.R) X_k (no weights).
      !> k_to_r_adjoint: X_k = w_k sum_R exp(+i k.R) X_R.
      !> Adjoint means the ordinary complex Euclidean inner product, not a BZ-weighted one.
      subroutine apply_transform(self, source, target, operation)
         import :: exchange_transform_type, wp
         class(exchange_transform_type), intent(inout) :: self
         complex(wp), intent(in) :: source(:, :)
         complex(wp), intent(out) :: target(:, :)
         integer, intent(in) :: operation
      end subroutine apply_transform
   end interface
end module tblite_transform_contract

#endif
