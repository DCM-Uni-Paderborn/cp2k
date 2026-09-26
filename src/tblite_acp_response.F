! SPDX-Identifier: LGPL-3.0-or-later
! Derived from save_tblite ACP response orchestration at aad26158b2e83f84b7b7d5b56f361ead4bc02a88.
! CP2K owns density contractions, image labels and stream state. Gaussian
! integrals and their physical derivatives remain in the provider.
#if defined(__TBLITE_GXTB_ACCELERATION)
module tblite_acp_response
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io, only: structure_type
   use mctc_io_constants, only: pi
   use mctc_io_math, only: matinv_3x3
   use tblite_acp_type, only: acp_type
   use tblite_acp_images, only: acp_image_type, get_pv_overlap_deriv_images
   use tblite_acp_cache, only: acp_cache
   use tblite_basis_type, only: basis_type
   use tblite_basis_cache, only: basis_cache
   use tblite_blas, only: gemm
   implicit none
   private
   public :: acp_kmesh_gradient_stream_begin
   public :: acp_kmesh_gradient_stream_push
   public :: acp_kmesh_gradient_stream_end
   public :: acp_kmesh_gradient_stream_discard
   public :: acp_kmesh_gradient_stream_real_elements
   public :: acp_kmesh_gradient_stream_complex_elements
   public :: get_acp_projector_image_gradient
   type, public :: acp_kmesh_gradient_stream_type
      private
      logical :: active = .false.
      integer :: nao = 0
      integer :: naux = 0
      integer :: nk_expected = 0
      integer :: nk_pushed = 0
      integer :: nspin = 0
      integer :: ntr = 0
      real(wp) :: weight_sum = 0.0_wp
      real(wp), allocatable :: projector_cell(:, :)
      real(wp), allocatable :: response(:, :, :, :)
      complex(wp), allocatable :: csk(:, :)
      complex(wp), allocatable :: work(:, :)
   end type acp_kmesh_gradient_stream_type
contains
subroutine acp_kmesh_gradient_stream_begin(stream, mol, images, nspin, &
   & nk_expected, error)
   type(acp_kmesh_gradient_stream_type), intent(inout) :: stream
   type(structure_type), intent(in) :: mol
   type(acp_image_type), intent(in) :: images
   integer, intent(in) :: nspin, nk_expected
   type(error_type), allocatable, intent(out) :: error

   integer :: nao, naux, ntr

   if (stream%active) then
      call fatal_error(error, "ACP k-mesh gradient stream is already active")
      return
   end if
   call clear_acp_kmesh_gradient_stream(stream)
   if (.not.allocated(images%projector_translation) .or. &
      & .not.allocated(images%projector_overlap) .or. &
      & .not.allocated(images%scaled_projector_overlap)) then
      call fatal_error(error, &
         & "ACP k-mesh gradient stream requires projector images")
      return
   end if
   naux = size(images%projector_overlap, 1)
   nao = size(images%projector_overlap, 2)
   ntr = size(images%projector_overlap, 3)
   if (naux < 1 .or. nao < 1 .or. ntr < 1 .or. nspin < 1 .or. &
      & nk_expected < 1 .or. &
      & any(shape(images%scaled_projector_overlap) /= &
      & shape(images%projector_overlap)) .or. &
      & any(shape(images%projector_translation) /= [3, ntr])) then
      call fatal_error(error, &
         & "inconsistent ACP k-mesh gradient stream dimensions")
      return
   end if
   if (.not.all(ieee_is_finite(images%projector_translation)) .or. &
      & .not.all(ieee_is_finite(images%projector_overlap)) .or. &
      & .not.all(ieee_is_finite(images%scaled_projector_overlap))) then
      call fatal_error(error, &
         & "ACP k-mesh gradient stream projector data must be finite")
      return
   end if

   call get_acp_projector_cells(mol, images, stream%projector_cell, error)
   if (allocated(error)) return
   allocate(stream%response(naux, nao, ntr, nspin), source=0.0_wp)
   allocate(stream%csk(naux, nao), stream%work(naux, nao))
   stream%nao = nao
   stream%naux = naux
   stream%ntr = ntr
   stream%nspin = nspin
   stream%nk_expected = nk_expected
   stream%nk_pushed = 0
   stream%weight_sum = 0.0_wp
   stream%active = .true.
end subroutine acp_kmesh_gradient_stream_begin

subroutine acp_kmesh_gradient_stream_push(stream, images, ik, kfrac, weight, &
   & density, error)
   type(acp_kmesh_gradient_stream_type), intent(inout) :: stream
   type(acp_image_type), intent(in) :: images
   integer, intent(in) :: ik
   real(wp), intent(in) :: kfrac(:), weight
   complex(wp), intent(in) :: density(:, :, :)
   type(error_type), allocatable, intent(out) :: error

   integer :: itr, spin
   real(wp) :: angle
   complex(wp) :: phase

   if (.not.stream%active) then
      call fatal_error(error, "ACP k-mesh gradient stream is not active")
      return
   end if
   if (ik /= stream%nk_pushed+1 .or. ik > stream%nk_expected) then
      call fatal_error(error, &
         & "ACP k-mesh gradient stream k-point order changed")
      return
   end if
   if (.not.allocated(images%projector_translation) .or. &
      & .not.allocated(images%scaled_projector_overlap)) then
      call fatal_error(error, &
         & "ACP k-mesh gradient stream lost its projector images")
      return
   end if
   if (size(kfrac) /= 3 .or. &
      & any(shape(density) /= [stream%nao, stream%nao, stream%nspin]) .or. &
      & any(shape(images%scaled_projector_overlap) /= &
      & [stream%naux, stream%nao, stream%ntr]) .or. &
      & any(shape(images%projector_translation) /= [3, stream%ntr])) then
      call fatal_error(error, &
         & "inconsistent ACP k-mesh gradient stream block dimensions")
      return
   end if
   if (.not.all(ieee_is_finite(kfrac)) .or. &
      & .not.ieee_is_finite(weight) .or. weight < 0.0_wp .or. &
      & .not.all(ieee_is_finite(real(density, wp))) .or. &
      & .not.all(ieee_is_finite(aimag(density)))) then
      call fatal_error(error, "invalid ACP k-mesh gradient stream block")
      return
   end if

   stream%csk = (0.0_wp, 0.0_wp)
   do itr = 1, stream%ntr
      angle = -2.0_wp*pi*dot_product(kfrac, &
         & stream%projector_cell(:, itr))
      phase = exp(cmplx(0.0_wp, angle, wp))
      stream%csk = stream%csk + &
         & phase*images%scaled_projector_overlap(:, :, itr)
   end do
   do spin = 1, stream%nspin
      stream%work = matmul(stream%csk, density(:, :, spin))
      do itr = 1, stream%ntr
         angle = -2.0_wp*pi*dot_product(kfrac, &
            & stream%projector_cell(:, itr))
         phase = exp(cmplx(0.0_wp, angle, wp))
         stream%response(:, :, itr, spin) = &
            & stream%response(:, :, itr, spin) + &
            & 2.0_wp*weight*real(conjg(phase)*stream%work, wp)
      end do
   end do
   stream%nk_pushed = ik
   stream%weight_sum = stream%weight_sum+weight
end subroutine acp_kmesh_gradient_stream_push

subroutine acp_kmesh_gradient_stream_end(stream, mol, bas, bcache, acp, &
   & acache, images, dEdcnbas, dEdqbas, gradient, sigma, error)
   type(acp_kmesh_gradient_stream_type), intent(inout) :: stream
   type(structure_type), intent(in) :: mol
   class(basis_type), intent(in) :: bas
   type(basis_cache), intent(in) :: bcache
   type(acp_type), intent(in) :: acp
   type(acp_cache), intent(in) :: acache
   type(acp_image_type), intent(in) :: images
   real(wp), intent(inout) :: dEdcnbas(:), dEdqbas(:)
   real(wp), intent(inout) :: gradient(:, :), sigma(:, :)
   type(error_type), allocatable, intent(out) :: error

   if (.not.stream%active) then
      call fatal_error(error, "ACP k-mesh gradient stream is not active")
      return
   end if
   if (stream%nk_pushed /= stream%nk_expected) then
      call fatal_error(error, &
         & "ACP k-mesh gradient stream ended with missing k points")
      return
   end if
   if (abs(stream%weight_sum-1.0_wp) > 1.0e-12_wp) then
      call fatal_error(error, &
         & "ACP k-mesh gradient stream weights are not normalized")
      return
   end if
   if (.not.allocated(images%projector_translation)) then
      call fatal_error(error, &
         & "ACP k-mesh gradient stream lost its projector translations")
      return
   end if
   if (bas%nao /= stream%nao .or. acp%auxbas%nao /= stream%naux .or. &
      & any(shape(images%projector_translation) /= [3, stream%ntr]) .or. &
      & size(dEdcnbas) /= mol%nat .or. size(dEdqbas) /= mol%nat .or. &
      & any(shape(gradient) /= [3, mol%nat]) .or. &
      & any(shape(sigma) /= [3, 3])) then
      call fatal_error(error, &
         & "inconsistent ACP k-mesh gradient stream output dimensions")
      return
   end if

   call get_pv_overlap_deriv_images(mol, images%projector_translation, bas, &
      & bcache, acp, acache, stream%response, dEdcnbas, dEdqbas, gradient, &
      & sigma)
   call clear_acp_kmesh_gradient_stream(stream)
end subroutine acp_kmesh_gradient_stream_end

subroutine acp_kmesh_gradient_stream_discard(stream)
   type(acp_kmesh_gradient_stream_type), intent(inout) :: stream

   call clear_acp_kmesh_gradient_stream(stream)
end subroutine acp_kmesh_gradient_stream_discard

integer function acp_kmesh_gradient_stream_real_elements(stream) &
   & result(nelements)
   type(acp_kmesh_gradient_stream_type), intent(in) :: stream

   nelements = 0
   if (allocated(stream%projector_cell)) &
      & nelements = nelements+size(stream%projector_cell)
   if (allocated(stream%response)) nelements = nelements+size(stream%response)
end function acp_kmesh_gradient_stream_real_elements

integer function acp_kmesh_gradient_stream_complex_elements(stream) &
   & result(nelements)
   type(acp_kmesh_gradient_stream_type), intent(in) :: stream

   nelements = 0
   if (allocated(stream%csk)) nelements = nelements+size(stream%csk)
   if (allocated(stream%work)) nelements = nelements+size(stream%work)
end function acp_kmesh_gradient_stream_complex_elements

subroutine clear_acp_kmesh_gradient_stream(stream)
   type(acp_kmesh_gradient_stream_type), intent(inout) :: stream

   if (allocated(stream%projector_cell)) deallocate(stream%projector_cell)
   if (allocated(stream%response)) deallocate(stream%response)
   if (allocated(stream%csk)) deallocate(stream%csk)
   if (allocated(stream%work)) deallocate(stream%work)
   stream%active = .false.
   stream%nao = 0
   stream%naux = 0
   stream%nk_expected = 0
   stream%nk_pushed = 0
   stream%nspin = 0
   stream%ntr = 0
   stream%weight_sum = 0.0_wp
end subroutine clear_acp_kmesh_gradient_stream

subroutine get_acp_projector_image_gradient(mol, bas, bcache, acp, acache, &
   & images, density_translation, density, dEdcnbas, dEdqbas, gradient, &
   & sigma, error)
   type(structure_type), intent(in) :: mol
   class(basis_type), intent(in) :: bas
   type(basis_cache), intent(in) :: bcache
   type(acp_type), intent(in) :: acp
   type(acp_cache), intent(in) :: acache
   type(acp_image_type), intent(in) :: images
   real(wp), intent(in) :: density_translation(:, :)
   real(wp), intent(in) :: density(:, :, :, :)
   real(wp), intent(inout) :: dEdcnbas(:), dEdqbas(:)
   real(wp), intent(inout) :: gradient(:, :), sigma(:, :)
   type(error_type), allocatable, intent(out) :: error

   integer :: dtr, itr, jtr, ktr, nspin, ntr, spin
   integer, allocatable :: density_map(:)
   real(wp) :: diff(3), response_residual, scale, tol2
   real(wp), allocatable :: density_cell(:, :), projector_cell(:, :), &
      & response(:, :, :, :), response_oracle(:, :, :, :)
   character(len=256) :: message

   if (.not.allocated(images%projector_translation) .or. &
      & .not.allocated(images%projector_overlap) .or. &
      & .not.allocated(images%scaled_projector_overlap)) then
      call fatal_error(error, &
         & "ACP sparse image gradient requires compact projector images")
      return
   end if
   ntr = size(images%projector_translation, 2)
   nspin = size(density, 4)
   if (size(images%projector_overlap, 1) /= acp%auxbas%nao .or. &
      & size(images%projector_overlap, 2) /= bas%nao .or. &
      & size(images%projector_overlap, 3) /= ntr .or. &
      & any(shape(images%scaled_projector_overlap) /= &
      & shape(images%projector_overlap)) .or. &
      & any(shape(images%projector_translation) /= [3, ntr]) .or. &
      & size(density_translation, 1) /= 3 .or. &
      & size(density, 1) /= bas%nao .or. size(density, 2) /= bas%nao .or. &
      & size(density, 3) /= size(density_translation, 2) .or. &
      & size(density, 3) < 1 .or. nspin < 1 .or. &
      & size(dEdcnbas) /= mol%nat .or. size(dEdqbas) /= mol%nat .or. &
      & any(shape(gradient) /= [3, mol%nat]) .or. &
      & any(shape(sigma) /= [3, 3])) then
      call fatal_error(error, &
         & "inconsistent ACP sparse image gradient dimensions")
      return
   end if
   if (.not.all(ieee_is_finite(density_translation)) .or. &
      & .not.all(ieee_is_finite(density))) then
      call fatal_error(error, &
         & "ACP sparse image gradient inputs must be finite")
      return
   end if

   call get_acp_translation_cells(mol, images%projector_translation, &
      & projector_cell, error)
   if (allocated(error)) return
   call get_acp_translation_cells(mol, density_translation, density_cell, error)
   if (allocated(error)) return
   allocate(response(acp%auxbas%nao, bas%nao, ntr, nspin), source=0.0_wp)
   do itr = 1, ntr
      do jtr = 1, ntr
         diff = projector_cell(:, itr)-projector_cell(:, jtr)
         dtr = find_translation(density_cell, diff, 0.0_wp)
         if (dtr > 0) then
            do spin = 1, nspin
               call gemm(amat=images%scaled_projector_overlap(:, :, jtr), &
                  & bmat=density(:, :, dtr, spin), &
                  & cmat=response(:, :, itr, spin), beta=1.0_wp)
            end do
         end if

         diff = projector_cell(:, jtr)-projector_cell(:, itr)
         dtr = find_translation(density_cell, diff, 0.0_wp)
         if (dtr > 0) then
            do spin = 1, nspin
               call gemm(amat=images%scaled_projector_overlap(:, :, jtr), &
                  & bmat=density(:, :, dtr, spin), &
                  & cmat=response(:, :, itr, spin), transb='T', beta=1.0_wp)
            end do
         end if
      end do
   end do

   if (allocated(images%translation) .and. allocated(images%hamiltonian)) then
      scale = max(1.0_wp, maxval(abs(images%translation)), &
         & maxval(abs(density_translation)))
      tol2 = (2048.0_wp*epsilon(1.0_wp)*scale)**2
      allocate(density_map(size(images%translation, 2)), source=0)
      do ktr = 1, size(images%translation, 2)
         density_map(ktr) = find_translation(density_translation, &
            & images%translation(:, ktr), tol2)
         if (density_map(ktr) == 0) then
            call fatal_error(error, &
               & "ACP sparse qualification lacks a dense-oracle translation")
            return
         end if
      end do
      allocate(response_oracle, mold=response)
      response_oracle = 0.0_wp
      do itr = 1, ntr
         do jtr = 1, ntr
            diff = images%projector_translation(:, itr) &
               & - images%projector_translation(:, jtr)
            ktr = find_translation(images%translation, diff, tol2)
            dtr = density_map(ktr)
            do spin = 1, nspin
               call gemm(amat=images%scaled_projector_overlap(:, :, jtr), &
                  & bmat=density(:, :, dtr, spin), &
                  & cmat=response_oracle(:, :, itr, spin), beta=1.0_wp)
            end do
            diff = images%projector_translation(:, jtr) &
               & - images%projector_translation(:, itr)
            ktr = find_translation(images%translation, diff, tol2)
            dtr = density_map(ktr)
            do spin = 1, nspin
               call gemm(amat=images%scaled_projector_overlap(:, :, jtr), &
                  & bmat=density(:, :, dtr, spin), &
                  & cmat=response_oracle(:, :, itr, spin), transb='T', &
                  & beta=1.0_wp)
            end do
         end do
      end do
      if (.not.all(ieee_is_finite(response)) .or. &
         & .not.all(ieee_is_finite(response_oracle))) then
         call fatal_error(error, &
            & "ACP sparse qualification produced a non-finite response")
         return
      end if
      scale = max(1.0_wp, maxval(abs(response_oracle)))
      response_residual = maxval(abs(response-response_oracle))/scale
      write(*, '(a,es12.4)') &
         & " ACP_SPARSE_RESPONSE_QUALIFY residual=", response_residual
      if (response_residual > 1.0e-12_wp) then
         write(message, '(a,es24.16,a,es24.16)') &
            & "ACP sparse response exceeded its dense-oracle gate: relative=", &
            & response_residual, ", scale=", scale
         call fatal_error(error, trim(message))
         return
      end if
   end if

   call get_pv_overlap_deriv_images(mol, images%projector_translation, bas, &
      & bcache, acp, acache, response, dEdcnbas, dEdqbas, gradient, sigma)
end subroutine get_acp_projector_image_gradient

subroutine get_acp_projector_cells(mol, images, cells, error)
   type(structure_type), intent(in) :: mol
   type(acp_image_type), intent(in) :: images
   real(wp), allocatable, intent(out) :: cells(:, :)
   type(error_type), allocatable, intent(out) :: error

   call get_acp_translation_cells(mol, images%projector_translation, cells, &
      & error)
end subroutine get_acp_projector_cells

subroutine get_acp_translation_cells(mol, translation, cells, error)
   type(structure_type), intent(in) :: mol
   real(wp), intent(in) :: translation(:, :)
   real(wp), allocatable, intent(out) :: cells(:, :)
   type(error_type), allocatable, intent(out) :: error

   integer :: itr
   real(wp) :: scale
   real(wp) :: invlat(3, 3)

   allocate(cells(3, size(translation, 2)), source=0.0_wp)
   if (.not.any(mol%periodic)) return
   invlat = matinv_3x3(mol%lattice)
   do itr = 1, size(cells, 2)
      cells(:, itr) = matmul(invlat, translation(:, itr))
   end do
   scale = max(1.0_wp, maxval(abs(cells)))
   if (maxval(abs(cells-real(nint(cells), wp))) > &
      & 4096.0_wp*epsilon(1.0_wp)*scale) then
      call fatal_error(error, "ACP translation is not a lattice vector")
      return
   end if
   cells = real(nint(cells), wp)
end subroutine get_acp_translation_cells

pure function find_translation(trans, target, tol2) result(index)
   real(wp), intent(in) :: trans(:, :)
   real(wp), intent(in) :: target(:)
   real(wp), intent(in) :: tol2
   integer :: index
   integer :: itr

   index = 0
   do itr = 1, size(trans, 2)
      if (sum((trans(:, itr) - target)**2) <= tol2) then
         index = itr
         exit
      end if
   end do
end function find_translation
end module tblite_acp_response
#endif
