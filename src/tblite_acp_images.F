! SPDX-Identifier: LGPL-3.0-or-later
! Image/Bloch assembly moved from our save_tblite integration at 25456f4.
! The original Seidler ACP implementation, integral kernels, parameters,
! sparse Gamma contractions and BLAS/OpenMP policy remain in the library.
#if defined(__TBLITE_GXTB)
module tblite_acp_images
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io, only: structure_type
   use mctc_io_constants, only: pi
   use mctc_io_math, only: matinv_3x3
   use tblite_acp_type, only: acp_type
   use tblite_acp_cache, only: acp_cache
   use tblite_basis_cache, only: basis_cache
   use tblite_basis_type, only: basis_type
   use tblite_blas, only: gemm
   use tblite_integral_native_integrals, only: overlap_cgto, overlap_grad_cgto, msao
   implicit none
   private
   public :: get_acp_images, get_acp_kmesh, get_acp_kmesh_gradient
   public :: get_acp_image_gradient, get_pv_overlap_deriv_images

   !> Image representation of the separable ACP operator.
   !>
   !> `projector_overlap(:, :, T)` stores
   !> \f$C_T=\langle p_T|\mu_0\rangle\f$.  Consequently its Bloch
   !> representation is \f$C(k)=\sum_T C_T\exp(-i kT)\f$.  The Hamiltonian
   !> images use the convention expected by CP2K,
   !> \f$H(k)=\sum_R H_R\exp(+i kR)\f$, and are obtained from the exact
   !> finite-image convolution
   !> \f$H_R=\sum_{T-V=R}C_T^T L C_V\f$.
   type, public :: acp_image_type
      real(wp), allocatable :: projector_translation(:, :)
      real(wp), allocatable :: projector_overlap(:, :, :)
      real(wp), allocatable :: scaled_projector_overlap(:, :, :)
      real(wp), allocatable :: level(:)
      real(wp), allocatable :: translation(:, :)
      integer, allocatable :: inverse(:)
      integer :: origin = 0
      logical :: inversion_closed = .false.
      real(wp), allocatable :: hamiltonian(:, :, :)
   end type acp_image_type

contains

subroutine get_acp_images(mol, trans, bas, bcache, acp, acache, images, &
   & build_hamiltonian)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Lattice points used for the projector-valence overlap
   real(wp), intent(in) :: trans(:, :)
   !> Valence basis set information
   class(basis_type), intent(in) :: bas
   !> Valence basis set cache
   type(basis_cache), intent(in) :: bcache
   !> Atomic correction potential data
   type(acp_type), intent(in) :: acp
   !> Atomic correction potential cache
   type(acp_cache), intent(inout) :: acache
   !> Image-resolved separable ACP operator
   type(acp_image_type), intent(out) :: images
   !> Form the quadratic projector-translation difference set and H(R).
   !> Direct Bloch-space callers only need the compact projector images.
   logical, intent(in), optional :: build_hamiltonian

   integer :: iat, jat, isp, jsp, itr, jtr, ktr, ninput, ntr, nhtr
   integer :: ish, jproj, is, js, ii, jj, iao, jao, nao, ij
   logical :: with_hamiltonian
   real(wp) :: diff(3), level, r2, scale, tol2, vec(3)
   real(wp), allocatable :: htrans_work(:, :), projector_overlap(:, :, :), &
      & scaled_projector_overlap(:, :, :), stmp(:)

   call acp%update(mol, acache)

   ninput = size(trans, 2)
   allocate(projector_overlap(acp%auxbas%nao, bas%nao, ninput), &
      & scaled_projector_overlap(acp%auxbas%nao, bas%nao, ninput), &
      & images%level(acp%auxbas%nao), source=0.0_wp)

   if (.not. allocated(acache%scaled_pv_overlap)) then
      allocate(acache%scaled_pv_overlap(acp%auxbas%nao, bas%nao), &
         & source=0.0_wp)
   else
      acache%scaled_pv_overlap(:, :) = 0.0_wp
   end if

   ! Expand the shell levels once to the auxiliary-AO representation.
   do jat = 1, mol%nat
      jsp = mol%id(jat)
      js = acp%auxbas%ish_at(jat)
      do jproj = 1, acp%auxbas%nsh_id(jsp)
         jj = acp%auxbas%iao_sh(js+jproj)
         level = acp%levels(jproj, jsp)
         nao = msao(acp%auxbas%cgto(jproj, jsp)%raw%ang)
         images%level(jj+1:jj+nao) = level
      end do
   end do

   allocate(stmp(msao(bas%maxl)*msao(acp%auxbas%maxl)))

   ! C_T = <p_T|mu_0>.  This is the same integral and finite image set used
   ! by get_pv_overlap_3d, only without summing over T prematurely.
   do iat = 1, mol%nat
      isp = mol%id(iat)
      is = bas%ish_at(iat)
      do jat = 1, mol%nat
         jsp = mol%id(jat)
         js = acp%auxbas%ish_at(jat)
         do itr = 1, ninput
            vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat) - trans(:, itr)
            r2 = dot_product(vec, vec)
            do ish = 1, bas%nsh_id(isp)
               ii = bas%iao_sh(is+ish)
               do jproj = 1, acp%auxbas%nsh_id(jsp)
                  jj = acp%auxbas%iao_sh(js+jproj)
                  call overlap_cgto(acp%auxbas%cgto(jproj, jsp)%raw, &
                     & bas%cgto(ish, isp)%raw, acache%auxbas%cgto(jproj, jat), &
                     & bcache%cgto(ish, iat), r2, vec, acp%auxbas%intcut, stmp)

                  level = acp%levels(jproj, jsp)
                  nao = msao(acp%auxbas%cgto(jproj, jsp)%raw%ang)
                  do iao = 1, msao(bas%cgto(ish, isp)%raw%ang)
                     do jao = 1, nao
                        ij = jao + nao*(iao-1)
                        projector_overlap(jj+jao, ii+iao, itr) = &
                           & stmp(ij)
                        scaled_projector_overlap(jj+jao, ii+iao, itr) = &
                           & level*stmp(ij)
                     end do
                  end do
               end do
            end do
         end do
      end do
   end do

   ! CP2K may supply a wider image range than the ACP integral support.
   ! overlap_cgto marks images beyond that support with exact zeros; remove
   ! only those empty images so that the direct Bloch representation and the
   ! real-space convolution remain independent of harmless outer padding.
   ntr = 0
   do itr = 1, ninput
      if (any(projector_overlap(:, :, itr) /= 0.0_wp)) ntr = ntr + 1
   end do
   allocate(images%projector_translation(3, ntr), &
      & images%projector_overlap(acp%auxbas%nao, bas%nao, ntr), &
      & images%scaled_projector_overlap(acp%auxbas%nao, bas%nao, ntr))
   ktr = 0
   do itr = 1, ninput
      if (.not.any(projector_overlap(:, :, itr) /= 0.0_wp)) cycle
      ktr = ktr + 1
      images%projector_translation(:, ktr) = trans(:, itr)
      images%projector_overlap(:, :, ktr) = projector_overlap(:, :, itr)
      images%scaled_projector_overlap(:, :, ktr) = &
         & scaled_projector_overlap(:, :, itr)
   end do

   acache%scaled_pv_overlap = sum(images%scaled_projector_overlap, dim=3)

   with_hamiltonian = .true.
   if (present(build_hamiltonian)) with_hamiltonian = build_hamiltonian
   if (.not.with_hamiltonian) return

   ! First construct the unique translation-difference set.  The tolerance
   ! only absorbs floating-point noise from equivalent lattice differences.
   allocate(htrans_work(3, max(1, ntr*ntr)), source=0.0_wp)
   scale = max(1.0_wp, maxval(abs(trans)))
   tol2 = (1024.0_wp*epsilon(1.0_wp)*scale)**2
   nhtr = 0
   if (ntr == 0) then
      nhtr = 1
      htrans_work(:, nhtr) = 0.0_wp
   else
      do itr = 1, ntr
         do jtr = 1, ntr
            diff = images%projector_translation(:, itr) &
               & - images%projector_translation(:, jtr)
            ktr = find_translation(htrans_work(:, :nhtr), diff, tol2)
            if (ktr == 0) then
               nhtr = nhtr + 1
               htrans_work(:, nhtr) = diff
            end if
         end do
      end do
   end if

   allocate(images%translation(3, nhtr), &
      & images%hamiltonian(bas%nao, bas%nao, nhtr), source=0.0_wp)
   allocate(images%inverse(nhtr), source=0)
   images%translation = htrans_work(:, :nhtr)
   images%origin = find_translation(images%translation, &
      & [0.0_wp, 0.0_wp, 0.0_wp], tol2)
   images%inversion_closed = .true.
   do itr = 1, nhtr
      images%inverse(itr) = find_translation(images%translation, &
         & -images%translation(:, itr), tol2)
      images%inversion_closed = images%inversion_closed .and. &
         & images%inverse(itr) > 0
   end do

   ! H_R = sum_(T-V=R) C_T^T L C_V.  Each ordered translation pair is
   ! included, so H_R^T=H_-R and the Bloch sum is Hermitian at arbitrary k.
   do itr = 1, ntr
      do jtr = 1, ntr
         diff = images%projector_translation(:, itr) &
            & - images%projector_translation(:, jtr)
         ktr = find_translation(images%translation, diff, tol2)
         call gemm(amat=images%projector_overlap(:, :, itr), &
            & bmat=images%scaled_projector_overlap(:, :, jtr), &
            & cmat=images%hamiltonian(:, :, ktr), transa='T', beta=1.0_wp)
      end do
   end do

end subroutine get_acp_images


subroutine get_acp_kmesh(mol, images, kfrac, hamiltonian, error)
   type(structure_type), intent(in) :: mol
   type(acp_image_type), intent(in) :: images
   real(wp), intent(in) :: kfrac(:, :)
   complex(wp), intent(out) :: hamiltonian(:, :, :)
   type(error_type), allocatable, intent(out) :: error

   integer :: ik, itr, nao, nk, ntr
   real(wp) :: angle
   real(wp), allocatable :: projector_cell(:, :)
   complex(wp) :: phase
   complex(wp), allocatable :: ck(:, :), csk(:, :)

   if (.not.allocated(images%projector_translation) .or. &
      & .not.allocated(images%projector_overlap) .or. &
      & .not.allocated(images%scaled_projector_overlap)) then
      call fatal_error(error, "ACP k-mesh evaluation requires projector images")
      return
   end if
   nao = size(images%projector_overlap, 2)
   ntr = size(images%projector_overlap, 3)
   nk = size(kfrac, 2)
   if (size(kfrac, 1) /= 3 .or. nk < 1 .or. &
      & any(shape(images%scaled_projector_overlap) /= &
      & shape(images%projector_overlap)) .or. &
      & any(shape(images%projector_translation) /= [3, ntr]) .or. &
      & any(shape(hamiltonian) /= [nao, nao, nk])) then
      call fatal_error(error, "inconsistent ACP k-mesh dimensions")
      return
   end if
   if (.not.all(ieee_is_finite(kfrac))) then
      call fatal_error(error, "ACP k points must be finite")
      return
   end if

   call get_acp_projector_cells(mol, images, projector_cell, error)
   if (allocated(error)) return
   allocate(ck(size(images%projector_overlap, 1), nao), &
      & csk(size(images%projector_overlap, 1), nao))
   hamiltonian = (0.0_wp, 0.0_wp)
   do ik = 1, nk
      ck = (0.0_wp, 0.0_wp)
      csk = (0.0_wp, 0.0_wp)
      do itr = 1, ntr
         angle = -2.0_wp*pi*dot_product(kfrac(:, ik), projector_cell(:, itr))
         phase = exp(cmplx(0.0_wp, angle, wp))
         ck = ck + phase*images%projector_overlap(:, :, itr)
         csk = csk + phase*images%scaled_projector_overlap(:, :, itr)
      end do
      hamiltonian(:, :, ik) = matmul(conjg(transpose(ck)), csk)
   end do
end subroutine get_acp_kmesh


subroutine get_acp_kmesh_gradient(mol, bas, bcache, acp, acache, images, &
   & kfrac, weights, density, dEdcnbas, dEdqbas, gradient, sigma, error)
   type(structure_type), intent(in) :: mol
   class(basis_type), intent(in) :: bas
   type(basis_cache), intent(in) :: bcache
   type(acp_type), intent(in) :: acp
   type(acp_cache), intent(in) :: acache
   type(acp_image_type), intent(in) :: images
   real(wp), intent(in) :: kfrac(:, :), weights(:)
   complex(wp), intent(in) :: density(:, :, :, :)
   real(wp), intent(inout) :: dEdcnbas(:), dEdqbas(:)
   real(wp), intent(inout) :: gradient(:, :), sigma(:, :)
   type(error_type), allocatable, intent(out) :: error

   integer :: ik, itr, nao, nk, nspin, ntr, spin
   real(wp) :: angle
   real(wp), allocatable :: projector_cell(:, :), response(:, :, :, :)
   complex(wp) :: phase
   complex(wp), allocatable :: csk(:, :), work(:, :)

   if (.not.allocated(images%projector_translation) .or. &
      & .not.allocated(images%projector_overlap) .or. &
      & .not.allocated(images%scaled_projector_overlap)) then
      call fatal_error(error, "ACP k-mesh gradient requires projector images")
      return
   end if
   nao = size(images%projector_overlap, 2)
   ntr = size(images%projector_overlap, 3)
   nk = size(weights)
   nspin = size(density, 3)
   if (size(kfrac, 1) /= 3 .or. size(kfrac, 2) /= nk .or. nk < 1 .or. &
      & size(density, 1) /= nao .or. size(density, 2) /= nao .or. &
      & size(density, 4) /= nk .or. nspin < 1 .or. &
      & any(shape(images%scaled_projector_overlap) /= &
      & shape(images%projector_overlap)) .or. &
      & any(shape(images%projector_translation) /= [3, ntr])) then
      call fatal_error(error, "inconsistent ACP k-mesh gradient dimensions")
      return
   end if
   if (.not.all(ieee_is_finite(kfrac)) .or. &
      & .not.all(ieee_is_finite(weights)) .or. &
      & .not.all(ieee_is_finite(real(density, wp))) .or. &
      & .not.all(ieee_is_finite(aimag(density))) .or. &
      & any(weights < 0.0_wp) .or. &
      & abs(sum(weights)-1.0_wp) > 1.0e-12_wp) then
      call fatal_error(error, "invalid ACP k-mesh density or weights")
      return
   end if

   call get_acp_projector_cells(mol, images, projector_cell, error)
   if (allocated(error)) return
   allocate(response(size(images%projector_overlap, 1), nao, ntr, nspin), &
      & source=0.0_wp)
   allocate(csk(size(images%projector_overlap, 1), nao), &
      & work(size(images%projector_overlap, 1), nao))
   do ik = 1, nk
      csk = (0.0_wp, 0.0_wp)
      do itr = 1, ntr
         angle = -2.0_wp*pi*dot_product(kfrac(:, ik), projector_cell(:, itr))
         phase = exp(cmplx(0.0_wp, angle, wp))
         csk = csk + phase*images%scaled_projector_overlap(:, :, itr)
      end do
      do spin = 1, nspin
         work = matmul(csk, density(:, :, spin, ik))
         do itr = 1, ntr
            angle = -2.0_wp*pi*dot_product(kfrac(:, ik), projector_cell(:, itr))
            phase = exp(cmplx(0.0_wp, angle, wp))
            response(:, :, itr, spin) = response(:, :, itr, spin) + &
               & 2.0_wp*weights(ik)*real(conjg(phase)*work, wp)
         end do
      end do
   end do

   call get_pv_overlap_deriv_images(mol, images%projector_translation, bas, &
      & bcache, acp, acache, response, dEdcnbas, dEdqbas, gradient, sigma)
end subroutine get_acp_kmesh_gradient


subroutine get_acp_projector_cells(mol, images, cells, error)
   type(structure_type), intent(in) :: mol
   type(acp_image_type), intent(in) :: images
   real(wp), allocatable, intent(out) :: cells(:, :)
   type(error_type), allocatable, intent(out) :: error

   integer :: itr
   real(wp) :: scale
   real(wp) :: invlat(3, 3)

   allocate(cells(3, size(images%projector_translation, 2)), source=0.0_wp)
   if (.not.any(mol%periodic)) return
   invlat = matinv_3x3(mol%lattice)
   do itr = 1, size(cells, 2)
      cells(:, itr) = matmul(invlat, images%projector_translation(:, itr))
   end do
   scale = max(1.0_wp, maxval(abs(cells)))
   if (maxval(abs(cells-real(nint(cells), wp))) > &
      & 4096.0_wp*epsilon(1.0_wp)*scale) then
      call fatal_error(error, "ACP projector translation is not a lattice vector")
      return
   end if
   cells = real(nint(cells), wp)
end subroutine get_acp_projector_cells


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


subroutine get_acp_image_gradient(mol, bas, bcache, acp, acache, images, &
   & density_translation, density, dEdcnbas, dEdqbas, gradient, sigma, error)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Valence basis set information
   class(basis_type), intent(in) :: bas
   !> Valence basis set cache
   type(basis_cache), intent(in) :: bcache
   !> Atomic correction potential data
   type(acp_type), intent(in) :: acp
   !> Atomic correction potential cache
   type(acp_cache), intent(in) :: acache
   !> Projector and Hamiltonian image data at the current geometry
   type(acp_image_type), intent(in) :: images
   !> Cartesian translations of the density images
   real(wp), intent(in) :: density_translation(:, :)
   !> Density images `(AO, AO, image, spin)` dual to `images%hamiltonian`
   real(wp), intent(in) :: density(:, :, :, :)
   !> Derivative w.r.t. basis-set coordination numbers
   real(wp), intent(inout) :: dEdcnbas(:)
   !> Derivative w.r.t. basis-set atomic charges
   real(wp), intent(inout) :: dEdqbas(:)
   !> Cartesian ACP gradient
   real(wp), intent(inout) :: gradient(:, :)
   !> ACP strain derivative
   real(wp), intent(inout) :: sigma(:, :)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: itr, jtr, ktr, dtr, ntr, nspin, nearest
   integer, allocatable :: density_map(:)
   real(wp) :: diff(3), nearest2, scale, tol2
   real(wp), allocatable :: response(:, :, :, :)
   character(len=512) :: message

   if (.not. allocated(images%projector_translation) .or. &
      & .not. allocated(images%projector_overlap) .or. &
      & .not. allocated(images%scaled_projector_overlap) .or. &
      & .not. allocated(images%translation) .or. &
      & .not. allocated(images%hamiltonian)) then
      call fatal_error(error, "ACP image gradient requires initialized image data")
      return
   end if

   ntr = size(images%projector_translation, 2)
   nspin = size(density, 4)
   if (size(images%projector_overlap, 1) /= acp%auxbas%nao .or. &
      & size(images%projector_overlap, 2) /= bas%nao .or. &
      & size(images%projector_overlap, 3) /= ntr .or. &
      & any(shape(images%scaled_projector_overlap) /= &
      & shape(images%projector_overlap))) then
      call fatal_error(error, "inconsistent ACP projector image dimensions")
      return
   end if
   if (size(density, 1) /= bas%nao .or. size(density, 2) /= bas%nao .or. &
      & size(density, 3) /= size(density_translation, 2) .or. &
      & size(density, 3) < 1 .or. nspin < 1) then
      call fatal_error(error, "inconsistent ACP density image dimensions")
      return
   end if
   if (size(density_translation, 1) /= 3 .or. &
      & size(images%translation, 1) /= 3 .or. &
      & size(images%projector_translation, 1) /= 3) then
      call fatal_error(error, "ACP image translations must be three-dimensional")
      return
   end if
   if (size(images%hamiltonian, 1) /= bas%nao .or. &
      & size(images%hamiltonian, 2) /= bas%nao .or. &
      & size(images%hamiltonian, 3) /= size(images%translation, 2)) then
      call fatal_error(error, "inconsistent ACP Hamiltonian image dimensions")
      return
   end if
   if (size(dEdcnbas) /= mol%nat .or. size(dEdqbas) /= mol%nat .or. &
      & any(shape(gradient) /= [3, mol%nat]) .or. &
      & any(shape(sigma) /= [3, 3])) then
      call fatal_error(error, "inconsistent ACP image gradient output dimensions")
      return
   end if

   scale = max(1.0_wp, maxval(abs(images%translation)), &
      & maxval(abs(density_translation)))
   tol2 = (2048.0_wp*epsilon(1.0_wp)*scale)**2
   allocate(density_map(size(images%translation, 2)), source=0)
   do ktr = 1, size(images%translation, 2)
      density_map(ktr) = find_translation(density_translation, &
         & images%translation(:, ktr), tol2)
      if (density_map(ktr) == 0) then
         nearest = 0
         nearest2 = huge(1.0_wp)
         do dtr = 1, size(density_translation, 2)
            diff = density_translation(:, dtr) - images%translation(:, ktr)
            if (sum(diff**2) < nearest2) then
               nearest = dtr
               nearest2 = sum(diff**2)
            end if
         end do
         write(message, '(a,i0,a,3(es16.8,1x),a,i0,a,i0,a,es16.8)') &
            & "density image set does not cover ACP translation ", ktr, &
            & " at [", images%translation(:, ktr), "]; density images=", &
            & size(density_translation, 2), ", nearest=", nearest, &
            & ", distance=", sqrt(nearest2)
         call fatal_error(error, trim(message))
         return
      end if
   end do

   allocate(response(acp%auxbas%nao, bas%nao, ntr, nspin), source=0.0_wp)
   do itr = 1, ntr
      do jtr = 1, ntr
         ! dC_U from the left projector factor:
         !       L C_V P_(U-V)^T
         diff = images%projector_translation(:, itr) &
            & - images%projector_translation(:, jtr)
         ktr = find_translation(images%translation, diff, tol2)
         if (ktr == 0) then
            call fatal_error(error, "incomplete ACP Hamiltonian difference set")
            return
         end if
         dtr = density_map(ktr)
         do ktr = 1, nspin
            call gemm(amat=images%scaled_projector_overlap(:, :, jtr), &
               & bmat=density(:, :, dtr, ktr), &
               & cmat=response(:, :, itr, ktr), transb='T', beta=1.0_wp)
         end do

         ! dC_U from the right projector factor:
         !       L C_T P_(T-U)
         diff = images%projector_translation(:, jtr) &
            & - images%projector_translation(:, itr)
         ktr = find_translation(images%translation, diff, tol2)
         if (ktr == 0) then
            call fatal_error(error, "incomplete ACP Hamiltonian difference set")
            return
         end if
         dtr = density_map(ktr)
         do ktr = 1, nspin
            call gemm(amat=images%scaled_projector_overlap(:, :, jtr), &
               & bmat=density(:, :, dtr, ktr), &
               & cmat=response(:, :, itr, ktr), beta=1.0_wp)
         end do
      end do
   end do

   call get_pv_overlap_deriv_images(mol, images%projector_translation, bas, &
      & bcache, acp, acache, response, dEdcnbas, dEdqbas, gradient, sigma)
end subroutine get_acp_image_gradient


subroutine get_pv_overlap_deriv_images(mol, trans, bas, bcache, acp, acache, &
   & response, dEdcnbas, dEdqbas, gradient, sigma)
   type(structure_type), intent(in) :: mol
   real(wp), intent(in) :: trans(:, :)
   class(basis_type), intent(in) :: bas
   type(basis_cache), intent(in) :: bcache
   type(acp_type), intent(in) :: acp
   type(acp_cache), intent(in) :: acache
   !> Coefficient of every projector-overlap derivative `(aux, AO, T, spin)`
   real(wp), intent(in) :: response(:, :, :, :)
   real(wp), intent(inout) :: dEdcnbas(:)
   real(wp), intent(inout) :: dEdqbas(:)
   real(wp), intent(inout) :: gradient(:, :)
   real(wp), intent(inout) :: sigma(:, :)

   integer :: iat, jat, isp, jsp, spin, nspin, itr
   integer :: ish, jproj, is, js, ii, jj, iao, jao, nao, ij
   real(wp) :: r2, vec(3), dG(3), dcnbasi, dqbasi, tmp
   real(wp), allocatable :: stmp(:), dstmp(:, :)
   real(wp), allocatable :: dstmpdqeffi(:), dstmpdqeffj(:)
   logical :: compute_qeff_grad
   real(wp), allocatable :: dEdcnbas_local(:), dEdqbas_local(:), &
      & gradient_local(:, :), sigma_local(:, :)

   compute_qeff_grad = bas%charge_dependent .or. acp%auxbas%charge_dependent
   nspin = size(response, 4)
   allocate(stmp(msao(bas%maxl)*msao(acp%auxbas%maxl)), &
      & dstmp(3, msao(bas%maxl)*msao(acp%auxbas%maxl)))
   if (compute_qeff_grad) then
      allocate(dstmpdqeffi(msao(bas%maxl)*msao(acp%auxbas%maxl)), &
         & dstmpdqeffj(msao(bas%maxl)*msao(acp%auxbas%maxl)))
   end if

   do spin = 1, nspin
      !$omp parallel default(none) shared(dEdcnbas, dEdqbas, gradient, sigma, trans) &
      !$omp shared(mol, bas, bcache, acp, acache, response, compute_qeff_grad, spin) &
      !$omp private(iat, jat, itr, isp, jsp, is, js, ish, jproj, ii, jj, iao, jao, nao, ij) &
      !$omp private(r2, vec, stmp, dstmp, dstmpdqeffi, dstmpdqeffj, dcnbasi, dqbasi) &
      !$omp private(tmp, dG, dEdcnbas_local, dEdqbas_local, gradient_local, sigma_local)
      if (compute_qeff_grad) then
         allocate(dEdcnbas_local(size(dEdcnbas)), source=0.0_wp)
         allocate(dEdqbas_local(size(dEdqbas)), source=0.0_wp)
      end if
      allocate(gradient_local(size(gradient, 1), size(gradient, 2)), source=0.0_wp)
      allocate(sigma_local(size(sigma, 1), size(sigma, 2)), source=0.0_wp)

      !$omp do schedule(runtime)
      do iat = 1, mol%nat
         isp = mol%id(iat)
         is = bas%ish_at(iat)
         dcnbasi = 0.0_wp
         dqbasi = 0.0_wp
         do jat = 1, mol%nat
            jsp = mol%id(jat)
            js = acp%auxbas%ish_at(jat)
            do itr = 1, size(trans, 2)
               vec = mol%xyz(:, iat) - mol%xyz(:, jat) - trans(:, itr)
               r2 = dot_product(vec, vec)
               dG = 0.0_wp
               do ish = 1, bas%nsh_id(isp)
                  ii = bas%iao_sh(is+ish)
                  do jproj = 1, acp%auxbas%nsh_id(jsp)
                     jj = acp%auxbas%iao_sh(js+jproj)
                     call overlap_grad_cgto(acp%auxbas%cgto(jproj, jsp)%raw, &
                        & bas%cgto(ish, isp)%raw, acache%auxbas%cgto(jproj, jat), &
                        & bcache%cgto(ish, iat), r2, vec, acp%auxbas%intcut, &
                        & stmp, dstmp, dstmpdqeffj, dstmpdqeffi)

                     nao = msao(acp%auxbas%cgto(jproj, jsp)%raw%ang)
                     do iao = 1, msao(bas%cgto(ish, isp)%raw%ang)
                        do jao = 1, nao
                           ij = jao + nao*(iao-1)
                           dG = dG + response(jj+jao, ii+iao, itr, spin) &
                              & * dstmp(:, ij)
                           if (compute_qeff_grad) then
                              tmp = response(jj+jao, ii+iao, itr, spin) &
                                 & * dstmpdqeffi(ij)
                              dcnbasi = dcnbasi &
                                 & + bcache%cgto(ish, iat)%dqeffdcn*tmp
                              dqbasi = dqbasi &
                                 & + bcache%cgto(ish, iat)%dqeffdq*tmp
                           end if
                        end do
                     end do
                  end do
               end do

               gradient_local(:, iat) = gradient_local(:, iat) + dG
               gradient_local(:, jat) = gradient_local(:, jat) - dG
               sigma_local = sigma_local + 0.5_wp*(spread(vec, 1, 3) &
                  & * spread(dG, 2, 3) + spread(dG, 1, 3)*spread(vec, 2, 3))
            end do
         end do
         if (compute_qeff_grad) then
            dEdcnbas_local(iat) = dEdcnbas_local(iat) + dcnbasi
            dEdqbas_local(iat) = dEdqbas_local(iat) + dqbasi
         end if
      end do
      !$omp end do

      !$omp critical
      if (compute_qeff_grad) then
         dEdcnbas = dEdcnbas + dEdcnbas_local
         dEdqbas = dEdqbas + dEdqbas_local
      end if
      gradient = gradient + gradient_local
      sigma = sigma + sigma_local
      !$omp end critical
      if (compute_qeff_grad) then
         deallocate(dEdcnbas_local, dEdqbas_local)
      end if
      deallocate(gradient_local, sigma_local)
      !$omp end parallel
   end do
end subroutine get_pv_overlap_deriv_images

end module tblite_acp_images
#endif
