! This file is part of tblite.
! SPDX-Identifier: LGPL-3.0-or-later
!
! tblite is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! tblite is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with tblite.  If not, see <https://www.gnu.org/licenses/>.

! CP2K integration tests relocated from save_tblite; original assertions retained.
#if defined(__TBLITE_GXTB_ACCELERATION)
module tblite_model_bridge_tests
   use tblite_projector_cache, only: cp2k_acp_projector_cache_type, projector_cache_counts
   use tblite_stream_plan, only: cp2k_stream_plan_type
   !$ use omp_lib, only : omp_get_max_threads, omp_set_num_threads
   use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type, new
   use mctc_io_constants, only : pi
   use mctc_io_convert, only : aatoau
   use mctc_io_math, only : matinv_3x3
   use tblite_acp_cache, only : acp_cache
   use tblite_acp_type, only : get_acp, get_acp_gradient
   use tblite_adjlist, only : adjacency_list, new_adjacency_list
   use tblite_basis_cache, only : basis_cache
   use tblite_basis_type, only : basis_type, get_cutoff, get_pair_cutoff
   use tblite_context_type, only : context_type
   use tblite_exchange_stream, only : acp_image_type, cp2k_acp_image_gradient, &
      & cp2k_acp_projector_image_gradient, cp2k_acp_kmesh, &
      & cp2k_acp_kmesh_gradient, cp2k_acp_kmesh_gradient_stream_type, &
      & cp2k_acp_kmesh_gradient_stream_begin, &
      & cp2k_acp_kmesh_gradient_stream_push, &
      & cp2k_acp_kmesh_gradient_stream_end, &
      & cp2k_acp_kmesh_gradient_stream_discard, &
      & cp2k_acp_kmesh_gradient_stream_real_elements, &
      & cp2k_acp_kmesh_gradient_stream_complex_elements, &
      & cp2k_apply_basis_gradient, cp2k_build_gamma_h0, &
      & cp2k_build_acp_projectors, cp2k_build_image_acp, cp2k_build_image_h0, &
      & cp2k_get_acp_cutoff, &
      & cp2k_get_scc_mixer_dimension, cp2k_potential_mixer_type, &
      & cp2k_update_basis, hamiltonian_image_type, new_cp2k_potential_mixer
   use tblite_cutoff, only : get_lattice_points
   use tblite_wavefunction, only : wavefunction_type, new_wavefunction, eeqbc_guess
   use tblite_xtb_calculator, only : xtb_calculator, xtb_config
   use tblite_xtb_gxtb, only : new_gxtb_calculator
   use tblite_xtb_h0, only : get_selfenergy
   use tblite_xtb_singlepoint, only : xtb_singlepoint
   implicit none
   private

   public :: collect_model_bridge

   real(wp), parameter :: acc = 0.00001_wp
   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr1 = 1e5*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))
   real(wp), parameter :: kt = 3.166808578545117e-06_wp


contains

subroutine collect_model_bridge(testsuite)
   type(unittest_type), allocatable, intent(out) :: testsuite(:)
   testsuite = [ &
      new_unittest("image-h0-h2o-pbc", &
         & test_image_h0_h2o_pbc), &
      new_unittest("image-acp-gradient-h2o-pbc", &
         & test_image_acp_gradient_h2o_pbc), &
      new_unittest("kmesh-acp-h2o-pbc", &
         & test_kmesh_acp_h2o_pbc), &
      new_unittest("cp2k-acp-cutoff", &
         & test_cp2k_acp_cutoff), &
      new_unittest("cp2k-scc-mixer-dimension", &
         & test_cp2k_scc_mixer_dimension), &
      new_unittest("cp2k-potential-mixer", &
         & test_cp2k_potential_mixer), &
      new_unittest("cp2k-potential-mixer-uninitialized", &
         & test_cp2k_potential_mixer_uninitialized, should_fail=.true.), &
      new_unittest("cp2k-potential-mixer-invalid-dimension", &
         & test_cp2k_potential_mixer_invalid_dimension, should_fail=.true.), &
      new_unittest("cp2k-potential-mixer-size-mismatch", &
         & test_cp2k_potential_mixer_size_mismatch, should_fail=.true.), &
      new_unittest("cp2k-potential-mixer-nonfinite", &
         & test_cp2k_potential_mixer_nonfinite, should_fail=.true.) ]
end subroutine collect_model_bridge

subroutine test_image_h0_h2o_pbc(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(acp_cache) :: acache
   type(acp_image_type) :: acp_images, acp_images_ref
   type(basis_cache) :: bcache
   type(hamiltonian_image_type) :: images, images_serial
   integer, parameter :: num(3) = [1, 1, 8]
   real(wp), parameter :: xyz(3, 3) = aatoau * reshape([&
      & 0.200000_wp, 0.300000_wp, 0.400000_wp, &
      & 0.958602_wp, 0.300000_wp, 0.904284_wp, &
      & 0.200000_wp, 0.300000_wp, 1.400000_wp], [3, 3])
   real(wp), parameter :: lattice(3, 3) = aatoau * reshape([&
      & 7.800000_wp, 0.200000_wp, 0.100000_wp, &
      & 0.400000_wp, 8.300000_wp, 0.300000_wp, &
      & 0.100000_wp, 0.500000_wp, 9.100000_wp], [3, 3])
   logical, parameter :: periodic(3) = [.true., .true., .true.]
   real(wp), parameter :: qvec(3) = [0.17_wp, 0.13_wp, 0.11_wp]
   integer :: itr, nao, saved_threads
   complex(wp) :: phase
   complex(wp), allocatable :: ck(:, :), csk(:, :), hk(:, :), hk_acp(:, :), &
      & hk_acp_direct(:, :), sk(:, :)
   real(wp), allocatable :: acp_gamma(:, :), cn(:), dipole(:, :, :), &
      & hamiltonian(:, :), hamiltonian_total(:, :), lattr(:, :), &
      & lattr_acp(:, :), lattr_acp_ref(:, :), &
      & overlap(:, :), quadrupole(:, :, :), selfenergy(:)

   call new(mol, num, xyz, lattice=lattice, periodic=periodic)
   call new_gxtb_calculator(calc, mol, error, accuracy=acc)
   if (allocated(error)) return
   call cp2k_update_basis(calc, mol, bcache, error)
   if (allocated(error)) return

   allocate(cn(mol%nat), selfenergy(calc%bas%nsh))
   call calc%ncoord%get_cn(mol, cn)
   call get_selfenergy(calc%h0, mol%id, calc%bas%ish_at, calc%bas%nsh_id, &
      & cn=cn, selfenergy=selfenergy)
   call get_lattice_points(mol%periodic, mol%lattice, &
      & get_cutoff(calc%bas, acc), lattr)
   call get_lattice_points(mol%periodic, mol%lattice, &
      & cp2k_get_acp_cutoff(calc, acc), lattr_acp_ref)

   nao = calc%bas%nao
   allocate(overlap(nao, nao), dipole(3, nao, nao), &
      & quadrupole(6, nao, nao), hamiltonian(nao, nao))
   saved_threads = 1
   !$ saved_threads = omp_get_max_threads()
   !$ call omp_set_num_threads(1)
   call cp2k_build_image_h0(calc, mol, bcache, acc, selfenergy, lattr, images_serial)
   !$ call omp_set_num_threads(saved_threads)
   call cp2k_build_image_h0(calc, mol, bcache, acc, selfenergy, lattr, images, &
      & overlap, dipole, quadrupole, hamiltonian)

   call check(error, maxval(abs(images%hamiltonian-images_serial%hamiltonian)), 0.0_wp, thr=1.0e-12_wp)
   if (allocated(error)) return
   call check(error, maxval(abs(images%overlap-images_serial%overlap)), 0.0_wp, thr=1.0e-12_wp)
   if (allocated(error)) return
   call check(error, maxval(abs(images%dipole_bra-images_serial%dipole_bra)), 0.0_wp, thr=1.0e-12_wp)
   if (allocated(error)) return
   call check(error, maxval(abs(images%dipole_ket-images_serial%dipole_ket)), 0.0_wp, thr=1.0e-12_wp)
   if (allocated(error)) return
   call check(error, maxval(abs(images%quadrupole_bra-images_serial%quadrupole_bra)), 0.0_wp, thr=1.0e-12_wp)
   if (allocated(error)) return
   call check(error, maxval(abs(images%quadrupole_ket-images_serial%quadrupole_ket)), 0.0_wp, thr=1.0e-12_wp)
   if (allocated(error)) return

   call check(error, maxval(abs(sum(images%overlap, dim=3) - overlap)), &
      & 0.0_wp, thr=1.0e-11_wp)
   if (allocated(error)) then
      call test_failed(error, "Image overlap does not reproduce Gamma overlap")
      return
   end if
   call check(error, maxval(abs(sum(images%dipole_bra, dim=4) - dipole)), &
      & 0.0_wp, thr=1.0e-11_wp)
   if (allocated(error)) then
      call test_failed(error, "Image dipoles do not reproduce Gamma dipoles")
      return
   end if
   call check(error, maxval(abs(sum(images%quadrupole_bra, dim=4) - quadrupole)), &
      & 0.0_wp, thr=1.0e-11_wp)
   if (allocated(error)) then
      call test_failed(error, "Image quadrupoles do not reproduce Gamma quadrupoles")
      return
   end if
   call check(error, maxval(abs(sum(images%hamiltonian, dim=3) - hamiltonian)), &
      & 0.0_wp, thr=1.0e-11_wp)
   if (allocated(error)) then
      call test_failed(error, "Image H0 does not reproduce Gamma H0")
      return
   end if

   allocate(sk(nao, nao), hk(nao, nao), source=(0.0_wp, 0.0_wp))
   do itr = 1, size(images%translation, 2)
      phase = exp(cmplx(0.0_wp, dot_product(qvec, images%translation(:, itr)), wp))
      sk = sk + phase * images%overlap(:, :, itr)
      hk = hk + phase * images%hamiltonian(:, :, itr)
   end do
   call check(error, maxval(abs(sk - conjg(transpose(sk)))), &
      & 0.0_wp, thr=1.0e-11_wp)
   if (allocated(error)) then
      call test_failed(error, "Bloch-summed image overlap is not Hermitian")
      return
   end if
   call check(error, maxval(abs(hk - conjg(transpose(hk)))), &
      & 0.0_wp, thr=1.0e-11_wp)
   if (allocated(error)) then
      call test_failed(error, "Bloch-summed image H0 is not Hermitian")
      return
   end if

   ! ACP is a separable projector convolution, so its Hamiltonian image range
   ! is the difference set of the projector-overlap translations.
   allocate(hamiltonian_total(nao, nao), acp_gamma(nao, nao))
   call cp2k_build_gamma_h0(calc, mol, bcache, acache, acc, selfenergy, &
      & overlap, dipole, quadrupole, hamiltonian_total)
   acp_gamma = hamiltonian_total - hamiltonian
   allocate(lattr_acp(3, size(lattr_acp_ref, 2)+2))
   lattr_acp(:, :size(lattr_acp_ref, 2)) = lattr_acp_ref
   lattr_acp(:, size(lattr_acp_ref, 2)+1) = 20.0_wp*sum(lattice, dim=2)
   lattr_acp(:, size(lattr_acp_ref, 2)+2) = -lattr_acp(:, size(lattr_acp_ref, 2)+1)
   call cp2k_build_image_acp(calc, mol, bcache, acache, lattr_acp_ref, acp_images_ref)
   call cp2k_build_image_acp(calc, mol, bcache, acache, lattr_acp, acp_images)
   call check(error, size(acp_images%projector_translation, 2) <= &
      & size(lattr_acp, 2)-2)
   if (allocated(error)) then
      call test_failed(error, "Zero ACP projector translations were not pruned")
      return
   end if
   do itr = 1, size(acp_images%projector_overlap, 3)
      call check(error, any(acp_images%projector_overlap(:, :, itr) /= 0.0_wp))
      if (allocated(error)) then
         call test_failed(error, "Compact ACP image set contains an empty image")
         return
      end if
   end do
   call check(error, size(acp_images%projector_translation, 2), &
      & size(acp_images_ref%projector_translation, 2))
   if (.not. allocated(error)) then
      call check(error, size(acp_images%translation, 2), &
         & size(acp_images_ref%translation, 2))
   end if
   if (allocated(error)) then
      call test_failed(error, "Empty outer translations changed compact ACP support")
      return
   end if
   call check(error, max( &
      & maxval(abs(acp_images%projector_translation &
      & - acp_images_ref%projector_translation)), &
      & maxval(abs(acp_images%projector_overlap &
      & - acp_images_ref%projector_overlap)), &
      & maxval(abs(acp_images%translation-acp_images_ref%translation)), &
      & maxval(abs(acp_images%hamiltonian-acp_images_ref%hamiltonian))), &
      & 0.0_wp, thr=1.0e-14_wp)
   if (allocated(error)) then
      call test_failed(error, "Expanded zero input changed compact ACP images")
      return
   end if
   call check(error, &
      & maxval(abs(sum(acp_images%hamiltonian, dim=3) - acp_gamma)), &
      & 0.0_wp, thr=1.0e-11_wp)
   if (allocated(error)) then
      call test_failed(error, "Image ACP convolution does not reproduce Gamma ACP")
      return
   end if

   allocate(ck(size(acp_images%projector_overlap, 1), nao), &
      & csk(size(acp_images%projector_overlap, 1), nao), &
      & hk_acp(nao, nao), hk_acp_direct(nao, nao), &
      & source=(0.0_wp, 0.0_wp))
   do itr = 1, size(acp_images%projector_translation, 2)
      phase = exp(cmplx(0.0_wp, -dot_product(qvec, &
         & acp_images%projector_translation(:, itr)), wp))
      ck = ck + phase*acp_images%projector_overlap(:, :, itr)
      csk = csk + phase*acp_images%scaled_projector_overlap(:, :, itr)
   end do
   do itr = 1, size(acp_images%translation, 2)
      phase = exp(cmplx(0.0_wp, &
         & dot_product(qvec, acp_images%translation(:, itr)), wp))
      hk_acp = hk_acp + phase*acp_images%hamiltonian(:, :, itr)
   end do
   hk_acp_direct = matmul(conjg(transpose(ck)), csk)
   call check(error, maxval(abs(hk_acp - hk_acp_direct)), &
      & 0.0_wp, thr=1.0e-11_wp)
   if (allocated(error)) then
      call test_failed(error, "Image ACP does not reproduce direct Bloch product")
      return
   end if
   call check(error, maxval(abs(hk_acp - conjg(transpose(hk_acp)))), &
      & 0.0_wp, thr=1.0e-11_wp)
   if (allocated(error)) then
      call test_failed(error, "Bloch-summed image ACP is not Hermitian")
      return
   end if

end subroutine test_image_h0_h2o_pbc

subroutine test_image_acp_gradient_h2o_pbc(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol, mol_work
   type(xtb_calculator) :: calc
   type(acp_cache) :: acache, acache_work, acache_gamma
   type(acp_image_type) :: images, images_work
   type(basis_cache) :: bcache, bcache_work
   type(wavefunction_type) :: wfn, wfn_aux, wfn_aux_work
   type(adjacency_list) :: list
   integer, parameter :: num(3) = [1, 1, 8]
   real(wp), parameter :: xyz(3, 3) = aatoau * reshape([&
      & 0.200000_wp, 0.300000_wp, 0.400000_wp, &
      & 0.958602_wp, 0.300000_wp, 0.904284_wp, &
      & 0.200000_wp, 0.300000_wp, 1.400000_wp], [3, 3])
   real(wp), parameter :: lattice(3, 3) = aatoau * reshape([&
      & 7.800000_wp, 0.200000_wp, 0.100000_wp, &
      & 0.400000_wp, 8.300000_wp, 0.300000_wp, &
      & 0.100000_wp, 0.500000_wp, 9.100000_wp], [3, 3])
   logical, parameter :: periodic(3) = [.true., .true., .true.]
   real(wp), parameter :: qvec(3) = [0.17_wp, 0.13_wp, 0.11_wp]
   real(wp), parameter :: step = 1.0e-5_wp
   real(wp), parameter :: fd_thr = 2.0e-7_wp
   integer :: iao, jao, iat, ic, ish, isp, jc, itr, nao, nimg
   real(wp) :: arg, el, er, eps(3, 3), unity(3, 3)
   real(wp), allocatable :: amat(:, :), bmat(:, :), density(:, :, :, :), &
      & density_gamma(:, :, :, :), density_translation(:, :), &
      & density_translation_work(:, :), dEdcn(:), dEdcn_gamma(:), &
      & dEdcn_image(:), dEdq(:), dEdq_gamma(:), dEdq_image(:), &
      & gradient(:, :), gradient_full(:, :), gradient_gamma(:, :), &
      & gradient_image(:, :), hamiltonian_gamma(:, :), &
      & lattr(:, :), lattr_work(:, :), num_dEdcn(:), num_dEdq(:), &
      & numgrad(:, :), numgrad_full(:, :), numsigma(:, :), &
      & numsigma_full(:, :)
   real(wp) :: sigma(3, 3), sigma_full(3, 3), sigma_gamma(3, 3), &
      & sigma_image(3, 3)

   call new(mol, num, xyz, lattice=lattice, periodic=periodic)
   call new_gxtb_calculator(calc, mol, error, accuracy=acc)
   if (allocated(error)) return
   call cp2k_update_basis(calc, mol, bcache, error, grad=.true., &
      & wfn_aux=wfn_aux)
   if (allocated(error)) return
   call get_lattice_points(mol%periodic, mol%lattice, &
      & cp2k_get_acp_cutoff(calc, acc), lattr)
   call cp2k_build_image_acp(calc, mol, bcache, acache, lattr, images)

   nao = calc%bas%nao
   nimg = size(images%translation, 2)
   allocate(amat(nao, nao), bmat(nao, nao), &
      & density(nao, nao, nimg, 1), density_gamma(nao, nao, nimg, 1), &
      & density_translation(3, nimg))
   density_translation = images%translation

   ! One artificial complex-Hermitian k-point density, transformed to real
   ! images.  A is symmetric and B antisymmetric, hence
   ! P(-R)=P(R)^T exactly (up to the trigonometric roundoff).
   do iao = 1, nao
      do jao = 1, nao
         amat(iao, jao) = 0.01_wp*(sin(real(iao + 3*jao, wp)) &
            & + sin(real(jao + 3*iao, wp)))
         bmat(iao, jao) = 0.01_wp*(cos(real(iao + 2*jao, wp)) &
            & - cos(real(jao + 2*iao, wp)))
      end do
   end do
   do itr = 1, nimg
      arg = dot_product(qvec, density_translation(:, itr))
      density(:, :, itr, 1) = cos(arg)*amat + sin(arg)*bmat
      density_gamma(:, :, itr, 1) = amat
   end do

   allocate(dEdcn(mol%nat), dEdq(mol%nat), gradient(3, mol%nat), &
      & numgrad(3, mol%nat), numsigma(3, 3), source=0.0_wp)
   sigma = 0.0_wp
   call cp2k_acp_image_gradient(calc, mol, bcache, acache, images, &
      & density_translation, density, dEdcn, dEdq, gradient, sigma, error)
   if (allocated(error)) return

   ! Independently perturb the effective q-vSZP basis coordinate along its
   ! CN and charge derivatives.  Re-normalizing every affected CGTO makes
   ! these finite differences include both coefficient and normalization
   ! response for the nonconstant image density.
   allocate(num_dEdcn(mol%nat), num_dEdq(mol%nat), source=0.0_wp)
   do iat = 1, mol%nat
      isp = mol%id(iat)
      do ish = 1, calc%bas%nsh_id(isp)
         bcache%cgto(ish, iat)%qeff = bcache%cgto(ish, iat)%qeff &
            & + step*bcache%cgto(ish, iat)%dqeffdcn
         call calc%bas%cgto(ish, isp)%raw%get_normalization( &
            & bcache%cgto(ish, iat), .false.)
      end do
      call cp2k_build_image_acp(calc, mol, bcache, acache, lattr, images_work)
      er = contract_acp_image_energy(images_work, density_translation, density)
      do ish = 1, calc%bas%nsh_id(isp)
         bcache%cgto(ish, iat)%qeff = bcache%cgto(ish, iat)%qeff &
            & - 2.0_wp*step*bcache%cgto(ish, iat)%dqeffdcn
         call calc%bas%cgto(ish, isp)%raw%get_normalization( &
            & bcache%cgto(ish, iat), .false.)
      end do
      call cp2k_build_image_acp(calc, mol, bcache, acache, lattr, images_work)
      el = contract_acp_image_energy(images_work, density_translation, density)
      num_dEdcn(iat) = 0.5_wp*(er-el)/step
      do ish = 1, calc%bas%nsh_id(isp)
         bcache%cgto(ish, iat)%qeff = bcache%cgto(ish, iat)%qeff &
            & + step*bcache%cgto(ish, iat)%dqeffdcn
         call calc%bas%cgto(ish, isp)%raw%get_normalization( &
            & bcache%cgto(ish, iat), .false.)
      end do
   end do
   call check(error, maxval(abs(dEdcn-num_dEdcn)), 0.0_wp, thr=fd_thr)
   if (allocated(error)) then
      call test_failed(error, "Image ACP dE/dCN does not match finite differences")
      return
   end if

   do iat = 1, mol%nat
      isp = mol%id(iat)
      do ish = 1, calc%bas%nsh_id(isp)
         bcache%cgto(ish, iat)%qeff = bcache%cgto(ish, iat)%qeff &
            & + step*bcache%cgto(ish, iat)%dqeffdq
         call calc%bas%cgto(ish, isp)%raw%get_normalization( &
            & bcache%cgto(ish, iat), .false.)
      end do
      call cp2k_build_image_acp(calc, mol, bcache, acache, lattr, images_work)
      er = contract_acp_image_energy(images_work, density_translation, density)
      do ish = 1, calc%bas%nsh_id(isp)
         bcache%cgto(ish, iat)%qeff = bcache%cgto(ish, iat)%qeff &
            & - 2.0_wp*step*bcache%cgto(ish, iat)%dqeffdq
         call calc%bas%cgto(ish, isp)%raw%get_normalization( &
            & bcache%cgto(ish, iat), .false.)
      end do
      call cp2k_build_image_acp(calc, mol, bcache, acache, lattr, images_work)
      el = contract_acp_image_energy(images_work, density_translation, density)
      num_dEdq(iat) = 0.5_wp*(er-el)/step
      do ish = 1, calc%bas%nsh_id(isp)
         bcache%cgto(ish, iat)%qeff = bcache%cgto(ish, iat)%qeff &
            & + step*bcache%cgto(ish, iat)%dqeffdq
         call calc%bas%cgto(ish, isp)%raw%get_normalization( &
            & bcache%cgto(ish, iat), .false.)
      end do
   end do
   call check(error, maxval(abs(dEdq-num_dEdq)), 0.0_wp, thr=fd_thr)
   if (allocated(error)) then
      call test_failed(error, "Image ACP dE/dq does not match finite differences")
      return
   end if

   ! Restore derivative-bearing normalization caches for the remaining
   ! comparisons (the finite differences above only need their values).
   call calc%bas%update(mol, bcache, .true., wfn_aux)
   allocate(gradient_full(3, mol%nat), source=gradient)
   sigma_full = sigma
   call cp2k_apply_basis_gradient(calc, mol, wfn_aux, dEdcn, dEdq, &
      & gradient_full, sigma_full)

   ! Frozen-density Cartesian finite differences of the exact image energy.
   do iat = 1, mol%nat
      do ic = 1, 3
         mol_work = mol
         mol_work%xyz(ic, iat) = mol_work%xyz(ic, iat) + step
         call cp2k_build_image_acp(calc, mol_work, bcache, acache, lattr, &
            & images_work)
         er = contract_acp_image_energy(images_work, density_translation, density)

         mol_work = mol
         mol_work%xyz(ic, iat) = mol_work%xyz(ic, iat) - step
         call cp2k_build_image_acp(calc, mol_work, bcache, acache, lattr, &
            & images_work)
         el = contract_acp_image_energy(images_work, density_translation, density)
         numgrad(ic, iat) = 0.5_wp*(er-el)/step
      end do
   end do
   call check(error, maxval(abs(gradient-numgrad)), 0.0_wp, thr=fd_thr)
   if (allocated(error)) then
      call test_failed(error, "Image ACP Cartesian gradient does not match finite differences")
      return
   end if

   ! Repeat the same displacement finite difference with the q-vSZP basis
   ! and its external charge model rebuilt at every geometry.  This checks
   ! the full ACP -> basis CN/q -> EEQ Cartesian chain rule.
   allocate(numgrad_full(3, mol%nat), source=0.0_wp)
   do iat = 1, mol%nat
      do ic = 1, 3
         mol_work = mol
         mol_work%xyz(ic, iat) = mol_work%xyz(ic, iat) + step
         call cp2k_update_basis(calc, mol_work, bcache_work, error, &
            & wfn_aux=wfn_aux_work)
         if (allocated(error)) return
         call cp2k_build_image_acp(calc, mol_work, bcache_work, acache_work, &
            & lattr, images_work)
         er = contract_acp_image_energy(images_work, density_translation, density)

         mol_work = mol
         mol_work%xyz(ic, iat) = mol_work%xyz(ic, iat) - step
         call cp2k_update_basis(calc, mol_work, bcache_work, error, &
            & wfn_aux=wfn_aux_work)
         if (allocated(error)) return
         call cp2k_build_image_acp(calc, mol_work, bcache_work, acache_work, &
            & lattr, images_work)
         el = contract_acp_image_energy(images_work, density_translation, density)
         numgrad_full(ic, iat) = 0.5_wp*(er-el)/step
      end do
   end do
   call check(error, maxval(abs(gradient_full-numgrad_full)), 0.0_wp, &
      & thr=fd_thr)
   if (allocated(error)) then
      call test_failed(error, &
         & "Full image ACP q-vSZP Cartesian response does not match finite differences")
      return
   end if

   ! Frozen-density strain finite differences.  The image labels follow the
   ! strained Cartesian lattice while their matrix values stay fixed.
   unity = 0.0_wp
   do ic = 1, 3
      unity(ic, ic) = 1.0_wp
   end do
   allocate(lattr_work(3, size(lattr, 2)), &
      & density_translation_work(3, nimg))
   do ic = 1, 3
      do jc = 1, 3
         eps = unity
         eps(jc, ic) = eps(jc, ic) + step
         mol_work = mol
         mol_work%xyz = matmul(eps, mol%xyz)
         mol_work%lattice = matmul(eps, mol%lattice)
         lattr_work = matmul(eps, lattr)
         density_translation_work = matmul(eps, density_translation)
         call cp2k_build_image_acp(calc, mol_work, bcache, acache, lattr_work, &
            & images_work)
         er = contract_acp_image_energy(images_work, &
            & density_translation_work, density)

         eps = unity
         eps(jc, ic) = eps(jc, ic) - step
         mol_work = mol
         mol_work%xyz = matmul(eps, mol%xyz)
         mol_work%lattice = matmul(eps, mol%lattice)
         lattr_work = matmul(eps, lattr)
         density_translation_work = matmul(eps, density_translation)
         call cp2k_build_image_acp(calc, mol_work, bcache, acache, lattr_work, &
            & images_work)
         el = contract_acp_image_energy(images_work, &
            & density_translation_work, density)
         numsigma(jc, ic) = 0.5_wp*(er-el)/step
      end do
   end do
   numsigma = 0.5_wp*(numsigma + transpose(numsigma))
   call check(error, maxval(abs(sigma-numsigma)), 0.0_wp, thr=fd_thr)
   if (allocated(error)) then
      call test_failed(error, "Image ACP strain derivative does not match finite differences")
      return
   end if

   ! Full q-vSZP/EEQ strain response, with all Cartesian lattice labels and
   ! basis caches transformed/rebuilt at the displaced cell.
   allocate(numsigma_full(3, 3), source=0.0_wp)
   do ic = 1, 3
      do jc = 1, 3
         eps = unity
         eps(jc, ic) = eps(jc, ic) + step
         mol_work = mol
         mol_work%xyz = matmul(eps, mol%xyz)
         mol_work%lattice = matmul(eps, mol%lattice)
         lattr_work = matmul(eps, lattr)
         density_translation_work = matmul(eps, density_translation)
         call cp2k_update_basis(calc, mol_work, bcache_work, error, &
            & wfn_aux=wfn_aux_work)
         if (allocated(error)) return
         call cp2k_build_image_acp(calc, mol_work, bcache_work, acache_work, &
            & lattr_work, images_work)
         er = contract_acp_image_energy(images_work, &
            & density_translation_work, density)

         eps = unity
         eps(jc, ic) = eps(jc, ic) - step
         mol_work = mol
         mol_work%xyz = matmul(eps, mol%xyz)
         mol_work%lattice = matmul(eps, mol%lattice)
         lattr_work = matmul(eps, lattr)
         density_translation_work = matmul(eps, density_translation)
         call cp2k_update_basis(calc, mol_work, bcache_work, error, &
            & wfn_aux=wfn_aux_work)
         if (allocated(error)) return
         call cp2k_build_image_acp(calc, mol_work, bcache_work, acache_work, &
            & lattr_work, images_work)
         el = contract_acp_image_energy(images_work, &
            & density_translation_work, density)
         numsigma_full(jc, ic) = 0.5_wp*(er-el)/step
      end do
   end do
   numsigma_full = 0.5_wp*(numsigma_full + transpose(numsigma_full))
   call check(error, maxval(abs(sigma_full-numsigma_full)), 0.0_wp, &
      & thr=fd_thr)
   if (allocated(error)) then
      call test_failed(error, &
         & "Full image ACP q-vSZP strain response does not match finite differences")
      return
   end if

   ! A constant symmetric density on every image is precisely the Gamma
   ! contraction.  This also checks the q-vSZP dE/dCN and dE/dq channels
   ! against the established Gamma implementation.
   call cp2k_build_image_acp(calc, mol, bcache, acache, lattr, images)
   call new_adjacency_list(list, mol, lattr, cp2k_get_acp_cutoff(calc, acc), complete=.true.)
   ! Build the independent Gamma oracle through the library's public path.
   ! Its gradient needs the compact AO maps as well as projector overlaps;
   ! the image cache deliberately does not construct those Gamma-only data.
   call calc%acp%update(mol, acache_gamma)
   allocate(hamiltonian_gamma(nao, nao), source=0.0_wp)
   call get_acp(mol, lattr, list, calc%bas, bcache, calc%acp, &
      & acache_gamma, hamiltonian_gamma)
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, &
      & 0.0_wp, .false.)
   wfn%density(:, :, 1) = amat
   allocate(dEdcn_gamma(mol%nat), dEdq_gamma(mol%nat), &
      & gradient_gamma(3, mol%nat), dEdcn_image(mol%nat), &
      & dEdq_image(mol%nat), gradient_image(3, mol%nat), source=0.0_wp)
   sigma_gamma = 0.0_wp
   sigma_image = 0.0_wp
   call get_acp_gradient(mol, lattr, list, calc%bas, bcache, calc%acp, &
      & acache_gamma, wfn, dEdcn_gamma, dEdq_gamma, gradient_gamma, sigma_gamma)
   call cp2k_acp_image_gradient(calc, mol, bcache, acache, images, &
      & density_translation, density_gamma, dEdcn_image, dEdq_image, &
      & gradient_image, sigma_image, error)
   if (allocated(error)) return
   call check(error, max(maxval(abs(gradient_image-gradient_gamma)), &
      & maxval(abs(sigma_image-sigma_gamma)), &
      & maxval(abs(dEdcn_image-dEdcn_gamma)), &
      & maxval(abs(dEdq_image-dEdq_gamma))), 0.0_wp, thr=1.0e-11_wp)
   if (allocated(error)) then
      call test_failed(error, "Constant ACP density images do not reproduce Gamma response")
   end if
end subroutine test_image_acp_gradient_h2o_pbc

subroutine test_kmesh_acp_h2o_pbc(error)
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol, mol_changed
   type(xtb_calculator) :: calc
   type(acp_cache) :: acache
   type(cp2k_acp_projector_cache_type) :: projector_cache
   integer :: cache_hits, cache_misses, cache_bypasses
   logical :: cache_valid
   type(acp_image_type) :: images, projectors, projectors_changed
   type(cp2k_acp_kmesh_gradient_stream_type) :: gradient_stream
   type(basis_cache) :: bcache, bcache_changed
   type(wavefunction_type) :: wfn_aux, wfn_aux_changed
   integer, parameter :: num(3) = [1, 1, 8]
   real(wp), parameter :: xyz(3, 3) = aatoau * reshape([&
      & 0.200000_wp, 0.300000_wp, 0.400000_wp, &
      & 0.958602_wp, 0.300000_wp, 0.904284_wp, &
      & 0.200000_wp, 0.300000_wp, 1.400000_wp], [3, 3])
   real(wp), parameter :: lattice(3, 3) = aatoau * reshape([&
      & 7.800000_wp, 0.200000_wp, 0.100000_wp, &
      & 0.400000_wp, 8.300000_wp, 0.300000_wp, &
      & 0.100000_wp, 0.500000_wp, 9.100000_wp], [3, 3])
   logical, parameter :: periodic(3) = [.true., .true., .true.]
   real(wp), parameter :: kfrac(3, 2) = reshape([&
      & 0.17_wp, 0.13_wp, 0.11_wp, &
      &-0.17_wp,-0.13_wp,-0.11_wp], [3, 2])
   real(wp), parameter :: weights(2) = [0.5_wp, 0.5_wp]
   integer :: iao, ik, itr, jao, nao, nimg, origin
   real(wp) :: angle, energy_cached, energy_image, energy_k, saved_alpha, &
      & saved_level
   real(wp) :: sigma_cached(3, 3), sigma_image(3, 3), sigma_k(3, 3), &
      & sigma_sparse(3, 3), sigma_stream(3, 3)
   real(wp), allocatable :: amat(:, :), bmat(:, :), cells(:, :), &
      & dEdcn_cached(:), dEdcn_image(:), dEdcn_k(:), dEdq_cached(:), &
      & dEdcn_sparse(:), dEdcn_stream(:), dEdq_image(:), dEdq_k(:), &
      & dEdq_sparse(:), dEdq_stream(:), &
      & density_cp2k(:, :, :, :), density_image(:, :, :, :), &
      & density_sparse(:, :, :, :), &
      & density_zero(:, :, :, :), gradient_cached(:, :), &
      & gradient_image(:, :), gradient_k(:, :), gradient_sparse(:, :), &
      & gradient_streamed(:, :), lattr(:, :), lattr_changed(:, :), &
      & translation_sparse(:, :)
   complex(wp) :: phase
   complex(wp), allocatable :: density_k(:, :, :, :), h_compact(:, :, :), &
      & h_image(:, :, :), h_k(:, :, :)

   call new(mol, num, xyz, lattice=lattice, periodic=periodic)
   call new_gxtb_calculator(calc, mol, error, accuracy=acc)
   if (allocated(error)) return
   call cp2k_update_basis(calc, mol, bcache, error, grad=.true., &
      & wfn_aux=wfn_aux)
   if (allocated(error)) return
   call get_lattice_points(mol%periodic, mol%lattice, &
      & cp2k_get_acp_cutoff(calc, acc), lattr)
   call cp2k_build_image_acp(calc, mol, bcache, acache, lattr, images, projector_cache=projector_cache)
   call projector_cache_counts(projector_cache, cache_hits, cache_misses, cache_bypasses, cache_valid)
   call check(error, cache_misses, 1)
   if (.not.allocated(error)) call check(error, cache_hits, 0)
   if (allocated(error)) then
      call test_failed(error, "Initial ACP projector build did not miss the cache")
      return
   end if

   nao = calc%bas%nao
   nimg = size(images%translation, 2)
   allocate(amat(nao, nao), bmat(nao, nao), cells(3, nimg), &
      & density_cp2k(nao, nao, nimg, 1), &
      & density_image(nao, nao, nimg, 1), source=0.0_wp)
   allocate(density_k(nao, nao, 1, 2), h_compact(nao, nao, 2), &
      & h_image(nao, nao, 2), h_k(nao, nao, 2), &
      & source=cmplx(0.0_wp, 0.0_wp, wp))
   do iao = 1, nao
      do jao = 1, nao
         amat(iao, jao) = 0.01_wp*(sin(real(iao + 3*jao, wp)) &
            & + sin(real(jao + 3*iao, wp)))
         bmat(iao, jao) = 0.01_wp*(cos(real(iao + 2*jao, wp)) &
            & - cos(real(jao + 2*iao, wp)))
      end do
   end do
   density_k(:, :, 1, 1) = cmplx(amat, bmat, wp)
   density_k(:, :, 1, 2) = conjg(density_k(:, :, 1, 1))
   cells = matmul(matinv_3x3(mol%lattice), images%translation)
   do itr = 1, nimg
      angle = 2.0_wp*pi*dot_product(kfrac(:, 1), cells(:, itr))
      density_image(:, :, itr, 1) = cos(angle)*amat + sin(angle)*bmat
      ! CP2K's gathered AO/image matrices carry the opposite cell phase
      ! relative to the provider's direct image-energy convention.
      density_cp2k(:, :, itr, 1) = cos(angle)*amat - sin(angle)*bmat
   end do

   call cp2k_acp_kmesh(calc, mol, images, kfrac, h_k, error)
   if (allocated(error)) return
   do ik = 1, 2
      do itr = 1, nimg
         angle = 2.0_wp*pi*dot_product(kfrac(:, ik), cells(:, itr))
         phase = exp(cmplx(0.0_wp, angle, wp))
         h_image(:, :, ik) = h_image(:, :, ik) + &
            & phase*images%hamiltonian(:, :, itr)
      end do
   end do
   call check(error, maxval(abs(h_k-h_image)), 0.0_wp, thr=1.0e-11_wp)
   if (allocated(error)) then
      call test_failed(error, "Direct k-space ACP does not match image convolution")
      return
   end if
   call cp2k_build_acp_projectors(calc, mol, bcache, acache, lattr, projectors, projector_cache=projector_cache)
   call projector_cache_counts(projector_cache, cache_hits, cache_misses, cache_bypasses, cache_valid)
   call check(error, cache_misses, 1)
   if (.not.allocated(error)) call check(error, cache_hits, 1)
   if (allocated(error)) then
      call test_failed(error, "Identical ACP projector build did not hit the cache")
      return
   end if
   if (allocated(projectors%translation) .or. &
      & allocated(projectors%hamiltonian)) then
      call test_failed(error, &
         & "Compact ACP projector build materialized the Hamiltonian difference set")
      return
   end if
   call cp2k_acp_kmesh(calc, mol, projectors, kfrac, h_compact, error)
   if (allocated(error)) return
   if (any(h_compact /= h_k)) then
      call test_failed(error, "Cached ACP projector build changed the Fock matrix")
      return
   end if
   call check(error, maxval(abs(h_compact-h_k)), 0.0_wp, thr=1.0e-11_wp)
   if (allocated(error)) then
      call test_failed(error, &
         & "Compact ACP projector build changed the direct k-space Hamiltonian")
      return
   end if
   energy_k = 0.0_wp
   do ik = 1, 2
      energy_k = energy_k + weights(ik)*real(sum(conjg(density_k(:, :, 1, ik)) &
         & *h_k(:, :, ik)), wp)
   end do
   energy_cached = 0.0_wp
   do ik = 1, 2
      energy_cached = energy_cached + weights(ik)*real(sum( &
         & conjg(density_k(:, :, 1, ik))*h_compact(:, :, ik)), wp)
   end do
   if (energy_cached /= energy_k) then
      call test_failed(error, "Cached ACP projector build changed the energy")
      return
   end if
   energy_image = contract_acp_image_energy(images, images%translation, &
      & density_image)
   call check(error, energy_k, energy_image, thr=1.0e-11_wp)
   if (allocated(error)) then
      call test_failed(error, &
         & "Weighted k-space ACP energy does not match the complete image reference")
      return
   end if

   allocate(dEdcn_cached(mol%nat), dEdcn_image(mol%nat), dEdcn_k(mol%nat), &
      & dEdcn_sparse(mol%nat), dEdcn_stream(mol%nat), &
      & dEdq_cached(mol%nat), dEdq_image(mol%nat), dEdq_k(mol%nat), &
      & dEdq_sparse(mol%nat), dEdq_stream(mol%nat), &
      & gradient_cached(3, mol%nat), gradient_image(3, mol%nat), &
      & gradient_k(3, mol%nat), gradient_sparse(3, mol%nat), &
      & gradient_streamed(3, mol%nat), source=0.0_wp)
   sigma_cached = 0.0_wp
   sigma_image = 0.0_wp
   sigma_k = 0.0_wp
   sigma_sparse = 0.0_wp
   sigma_stream = 0.0_wp
   call cp2k_acp_image_gradient(calc, mol, bcache, acache, images, &
      & images%translation, density_image, dEdcn_image, dEdq_image, &
      & gradient_image, sigma_image, error)
   if (allocated(error)) return
   call cp2k_acp_projector_image_gradient(calc, mol, bcache, acache, images, &
      & images%translation, density_cp2k, dEdcn_sparse, dEdq_sparse, &
      & gradient_sparse, sigma_sparse, error)
   if (allocated(error)) return
   call cp2k_acp_kmesh_gradient(calc, mol, bcache, acache, images, &
      & kfrac, weights, density_k, dEdcn_k, dEdq_k, gradient_k, sigma_k, error)
   if (allocated(error)) return
   call cp2k_acp_kmesh_gradient(calc, mol, bcache, acache, projectors, &
      & kfrac, weights, density_k, dEdcn_cached, dEdq_cached, &
      & gradient_cached, sigma_cached, error)
   if (allocated(error)) return
   ! Independently evaluated OpenMP derivative reductions need not add their
   ! thread-local contributions in the same order. Keep a roundoff-level
   ! numerical gate; bitwise equality remains required for cached Fock data.
   call check(error, max(maxval(abs(dEdcn_cached-dEdcn_k)), &
      & maxval(abs(dEdq_cached-dEdq_k)), maxval(abs(gradient_cached-gradient_k)), &
      & maxval(abs(sigma_cached-sigma_k))), 0.0_wp, thr=16*epsilon(1.0_wp))
   if (allocated(error)) then
      write(*, '(a,4es24.15)') 'ACP cached response residuals: ', &
         & maxval(abs(dEdcn_cached-dEdcn_k)), maxval(abs(dEdq_cached-dEdq_k)), &
         & maxval(abs(gradient_cached-gradient_k)), maxval(abs(sigma_cached-sigma_k))
      call test_failed(error, &
         & "Cached ACP projector build changed the force or stress response")
      return
   end if
   call check(error, max(maxval(abs(dEdcn_k-dEdcn_image)), &
      & maxval(abs(dEdq_k-dEdq_image)), &
      & maxval(abs(gradient_k-gradient_image)), &
      & maxval(abs(sigma_k-sigma_image))), 0.0_wp, thr=1.0e-11_wp)
   if (allocated(error)) then
      call test_failed(error, &
         & "Direct complex k-space ACP reverse sweep does not match image response")
      return
   end if
   call check(error, max(maxval(abs(dEdcn_sparse-dEdcn_k)), &
      & maxval(abs(dEdq_sparse-dEdq_k)), &
      & maxval(abs(gradient_sparse-gradient_k)), &
      & maxval(abs(sigma_sparse-sigma_k))), 0.0_wp, thr=1.0e-11_wp)
   if (allocated(error)) then
      call test_failed(error, &
         & "Sparse CP2K-image ACP derivative changed the complex k-space response")
      return
   end if

   ! The bounded reverse stream must reproduce the whole-mesh reverse sweep
   ! while owning no AO-by-AO-by-k density or Hamiltonian payload.
   call cp2k_acp_kmesh_gradient_stream_begin(gradient_stream, calc, mol, &
      & projectors, 1, 2, error)
   if (allocated(error)) return
   call check(error, cp2k_acp_kmesh_gradient_stream_real_elements( &
      & gradient_stream), 3*size(projectors%projector_translation, 2) + &
      & calc%acp%auxbas%nao*nao*size(projectors%projector_translation, 2))
   if (.not.allocated(error)) then
      call check(error, cp2k_acp_kmesh_gradient_stream_complex_elements( &
         & gradient_stream), 2*calc%acp%auxbas%nao*nao)
   end if
   if (allocated(error)) then
      call test_failed(error, "ACP reverse stream storage is not mesh bounded")
      return
   end if
   do ik = 1, 2
      call cp2k_acp_kmesh_gradient_stream_push(gradient_stream, projectors, &
         & ik, kfrac(:, ik), weights(ik), density_k(:, :, :, ik), error)
      if (allocated(error)) return
   end do
   call cp2k_acp_kmesh_gradient_stream_end(gradient_stream, calc, mol, &
      & bcache, acache, projectors, dEdcn_stream, dEdq_stream, &
      & gradient_streamed, sigma_stream, error)
   if (allocated(error)) return
   call check(error, max(maxval(abs(dEdcn_stream-dEdcn_k)), &
      & maxval(abs(dEdq_stream-dEdq_k)), &
      & maxval(abs(gradient_streamed-gradient_k)), &
      & maxval(abs(sigma_stream-sigma_k))), 0.0_wp, thr=1.0e-11_wp)
   if (allocated(error)) then
      call test_failed(error, &
         & "Streamed ACP reverse sweep changed force or stress response")
      return
   end if

   ! The stream is fail-closed: no reordering, incomplete transaction or
   ! unnormalized irreducible-star weights may reach the derivative kernel.
   call cp2k_acp_kmesh_gradient_stream_begin(gradient_stream, calc, mol, &
      & projectors, 1, 2, error)
   if (allocated(error)) return
   call cp2k_acp_kmesh_gradient_stream_push(gradient_stream, projectors, &
      & 2, kfrac(:, 2), weights(2), density_k(:, :, :, 2), error)
   if (.not.allocated(error)) then
      call cp2k_acp_kmesh_gradient_stream_discard(gradient_stream)
      call test_failed(error, "Reordered ACP reverse-stream k point was accepted")
      return
   end if
   if (trim(error%message) /= &
      & "ACP k-mesh gradient stream k-point order changed") then
      deallocate(error)
      call cp2k_acp_kmesh_gradient_stream_discard(gradient_stream)
      call test_failed(error, "Unexpected ACP reverse-stream order diagnostic")
      return
   end if
   deallocate(error)
   call cp2k_acp_kmesh_gradient_stream_discard(gradient_stream)

   call cp2k_acp_kmesh_gradient_stream_begin(gradient_stream, calc, mol, &
      & projectors, 1, 2, error)
   if (allocated(error)) return
   call cp2k_acp_kmesh_gradient_stream_push(gradient_stream, projectors, &
      & 1, kfrac(:, 1), weights(1), density_k(:, :, :, 1), error)
   if (allocated(error)) return
   call cp2k_acp_kmesh_gradient_stream_end(gradient_stream, calc, mol, &
      & bcache, acache, projectors, dEdcn_stream, dEdq_stream, &
      & gradient_streamed, sigma_stream, error)
   if (.not.allocated(error)) then
      call cp2k_acp_kmesh_gradient_stream_discard(gradient_stream)
      call test_failed(error, "Incomplete ACP reverse stream was accepted")
      return
   end if
   if (trim(error%message) /= &
      & "ACP k-mesh gradient stream ended with missing k points") then
      deallocate(error)
      call cp2k_acp_kmesh_gradient_stream_discard(gradient_stream)
      call test_failed(error, "Unexpected incomplete ACP stream diagnostic")
      return
   end if
   deallocate(error)
   call cp2k_acp_kmesh_gradient_stream_discard(gradient_stream)

   call cp2k_acp_kmesh_gradient_stream_begin(gradient_stream, calc, mol, &
      & projectors, 1, 2, error)
   if (allocated(error)) return
   do ik = 1, 2
      call cp2k_acp_kmesh_gradient_stream_push(gradient_stream, projectors, &
         & ik, kfrac(:, ik), 0.5_wp*weights(ik), &
         & density_k(:, :, :, ik), error)
      if (allocated(error)) return
   end do
   call cp2k_acp_kmesh_gradient_stream_end(gradient_stream, calc, mol, &
      & bcache, acache, projectors, dEdcn_stream, dEdq_stream, &
      & gradient_streamed, sigma_stream, error)
   if (.not.allocated(error)) then
      call cp2k_acp_kmesh_gradient_stream_discard(gradient_stream)
      call test_failed(error, "Unnormalized ACP reverse stream was accepted")
      return
   end if
   if (trim(error%message) /= &
      & "ACP k-mesh gradient stream weights are not normalized") then
      deallocate(error)
      call cp2k_acp_kmesh_gradient_stream_discard(gradient_stream)
      call test_failed(error, "Unexpected ACP stream weight diagnostic")
      return
   end if
   deallocate(error)
   call cp2k_acp_kmesh_gradient_stream_discard(gradient_stream)

   ! A compact projector derivative may consume CP2K's sparse density-image
   ! support directly.  Missing translation differences are exact zero blocks,
   ! so this must equal the established zero-padded difference-set oracle.
   origin = images%origin
   if (origin < 1) then
      call test_failed(error, "ACP difference set has no origin image")
      return
   end if
   allocate(density_zero(nao, nao, nimg, 1), source=0.0_wp)
   allocate(density_sparse(nao, nao, 1, 1), source=0.0_wp)
   allocate(translation_sparse(3, 1))
   density_zero(:, :, origin, 1) = amat
   density_sparse(:, :, 1, 1) = amat
   translation_sparse(:, 1) = images%translation(:, origin)
   dEdcn_image = 0.0_wp
   dEdq_image = 0.0_wp
   gradient_image = 0.0_wp
   sigma_image = 0.0_wp
   dEdcn_sparse = 0.0_wp
   dEdq_sparse = 0.0_wp
   gradient_sparse = 0.0_wp
   sigma_sparse = 0.0_wp
   call cp2k_acp_image_gradient(calc, mol, bcache, acache, images, &
      & images%translation, density_zero, dEdcn_image, dEdq_image, &
      & gradient_image, sigma_image, error)
   if (allocated(error)) return
   call cp2k_acp_projector_image_gradient(calc, mol, bcache, acache, &
      & projectors, translation_sparse, density_sparse, dEdcn_sparse, &
      & dEdq_sparse, gradient_sparse, sigma_sparse, error)
   if (allocated(error)) return
   call check(error, max(maxval(abs(dEdcn_sparse-dEdcn_image)), &
      & maxval(abs(dEdq_sparse-dEdq_image)), &
      & maxval(abs(gradient_sparse-gradient_image)), &
      & maxval(abs(sigma_sparse-sigma_image))), 0.0_wp, thr=1.0e-11_wp)
   if (allocated(error)) then
      call test_failed(error, &
         & "Sparse ACP projector derivative changed the image response")
      return
   end if
   write(*, '(a,1x,i0,1x,a,1x,i0)') "ACP_STREAM_STATE_GATES", 3, &
      & "ACP_SPARSE_DENSITY_IMAGES", size(translation_sparse, 2)

   ! A Cartesian displacement must invalidate both the geometry and effective
   ! q-vSZP signature, even when the lattice translation set is unchanged.
   mol_changed = mol
   mol_changed%xyz(1, 1) = mol_changed%xyz(1, 1) + 1.0e-6_wp
   call cp2k_update_basis(calc, mol_changed, bcache_changed, error, grad=.true., &
      & wfn_aux=wfn_aux_changed)
   if (allocated(error)) return
   call cp2k_build_acp_projectors(calc, mol_changed, bcache_changed, acache, &
      & lattr, projectors_changed, projector_cache=projector_cache)
   call projector_cache_counts(projector_cache, cache_hits, cache_misses, cache_bypasses, cache_valid)
   call check(error, cache_misses, 2)
   if (allocated(error)) then
      call test_failed(error, "ACP projector cache survived a geometry change")
      return
   end if

   ! A cell deformation is an independent invalidation key, including the
   ! regenerated lattice-translation range used by the projector integrals.
   mol_changed = mol
   mol_changed%lattice(1, 1) = mol_changed%lattice(1, 1) + 1.0e-6_wp
   call cp2k_update_basis(calc, mol_changed, bcache_changed, error, grad=.true., &
      & wfn_aux=wfn_aux_changed)
   if (allocated(error)) return
   call get_lattice_points(mol_changed%periodic, mol_changed%lattice, &
      & cp2k_get_acp_cutoff(calc, acc), lattr_changed)
   call cp2k_build_acp_projectors(calc, mol_changed, bcache_changed, acache, &
      & lattr_changed, projectors_changed, projector_cache=projector_cache)
   call projector_cache_counts(projector_cache, cache_hits, cache_misses, cache_bypasses, cache_valid)
   call check(error, cache_misses, 3)
   if (allocated(error)) then
      call test_failed(error, "ACP projector cache survived a cell change")
      return
   end if

   ! Re-prime the original geometry, then perturb only the live effective
   ! q-vSZP contraction.  This isolates basis-cache invalidation from geometry.
   call cp2k_build_acp_projectors(calc, mol, bcache, acache, lattr, &
      & projectors_changed, projector_cache=projector_cache)
   call projector_cache_counts(projector_cache, cache_hits, cache_misses, cache_bypasses, cache_valid)
   call check(error, cache_misses, 4)
   if (allocated(error)) return
   bcache_changed = bcache
   bcache_changed%cgto(1, 1)%norm = &
      & bcache_changed%cgto(1, 1)%norm*(1.0_wp + 1.0e-6_wp)
   call cp2k_build_acp_projectors(calc, mol, bcache_changed, acache, lattr, &
      & projectors_changed, projector_cache=projector_cache)
   call projector_cache_counts(projector_cache, cache_hits, cache_misses, cache_bypasses, cache_valid)
   call check(error, cache_misses, 5)
   if (allocated(error)) then
      call test_failed(error, &
         & "ACP projector cache survived an effective-CGTO change")
      return
   end if

   ! Raw primitive metadata and ACP levels are independent parts of the key.
   call cp2k_build_acp_projectors(calc, mol, bcache, acache, lattr, &
      & projectors_changed, projector_cache=projector_cache)
   call projector_cache_counts(projector_cache, cache_hits, cache_misses, cache_bypasses, cache_valid)
   call check(error, cache_misses, 6)
   if (allocated(error)) return
   associate(cgto => calc%bas%cgto(1, mol%id(1))%raw)
      saved_alpha = cgto%alpha(1)
      cgto%alpha(1) = saved_alpha*(1.0_wp + 1.0e-6_wp)
   end associate
   call cp2k_build_acp_projectors(calc, mol, bcache, acache, lattr, &
      & projectors_changed, projector_cache=projector_cache)
   call projector_cache_counts(projector_cache, cache_hits, cache_misses, cache_bypasses, cache_valid)
   call check(error, cache_misses, 7)
   associate(cgto => calc%bas%cgto(1, mol%id(1))%raw)
      cgto%alpha(1) = saved_alpha
   end associate
   if (allocated(error)) then
      call test_failed(error, "ACP projector cache survived a primitive change")
      return
   end if

   call cp2k_build_acp_projectors(calc, mol, bcache, acache, lattr, &
      & projectors_changed, projector_cache=projector_cache)
   call projector_cache_counts(projector_cache, cache_hits, cache_misses, cache_bypasses, cache_valid)
   call check(error, cache_misses, 8)
   if (allocated(error)) return
   saved_level = calc%acp%levels(1, mol%id(1))
   calc%acp%levels(1, mol%id(1)) = saved_level + 1.0e-6_wp
   call cp2k_build_acp_projectors(calc, mol, bcache, acache, lattr, &
      & projectors_changed, projector_cache=projector_cache)
   call projector_cache_counts(projector_cache, cache_hits, cache_misses, cache_bypasses, cache_valid)
   call check(error, cache_misses, 9)
   calc%acp%levels(1, mol%id(1)) = saved_level
   if (allocated(error)) then
      call test_failed(error, "ACP projector cache survived an ACP-level change")
      return
   end if

   ! The explicit direct-build oracle must reproduce every ACP observable,
   ! invalidate the prior payload before rebuilding, and leave no reusable
   ! state behind.  The following normal build therefore has to miss once and
   ! only the subsequent identical build may hit.
   call cp2k_build_acp_projectors(calc, mol, bcache, acache, lattr, &
      & projectors_changed, use_projector_cache=.false., projector_cache=projector_cache)
   call projector_cache_counts(projector_cache, cache_hits, cache_misses, cache_bypasses, cache_valid)
   call check(error, cache_bypasses, 1)
   if (.not.allocated(error)) call check(error, .not.cache_valid)
   if (allocated(error)) then
      call test_failed(error, &
         & "Direct ACP oracle retained a stale projector cache")
      return
   end if
   call cp2k_acp_kmesh(calc, mol, projectors_changed, kfrac, h_compact, error)
   if (allocated(error)) return
   if (any(h_compact /= h_k)) then
      call test_failed(error, "Direct ACP oracle changed the Fock matrix")
      return
   end if
   energy_cached = 0.0_wp
   do ik = 1, 2
      energy_cached = energy_cached + weights(ik)*real(sum( &
         & conjg(density_k(:, :, 1, ik))*h_compact(:, :, ik)), wp)
   end do
   call check(error, energy_cached, energy_k, thr=1.0e-15_wp)
   if (allocated(error)) then
      call test_failed(error, "Direct ACP oracle changed the energy")
      return
   end if
   dEdcn_cached = 0.0_wp
   dEdq_cached = 0.0_wp
   gradient_cached = 0.0_wp
   sigma_cached = 0.0_wp
   call cp2k_acp_kmesh_gradient(calc, mol, bcache, acache, &
      & projectors_changed, kfrac, weights, density_k, dEdcn_cached, &
      & dEdq_cached, gradient_cached, sigma_cached, error)
   if (allocated(error)) return
   call check(error, max(maxval(abs(dEdcn_cached-dEdcn_k)), &
      & maxval(abs(dEdq_cached-dEdq_k)), maxval(abs(gradient_cached-gradient_k)), &
      & maxval(abs(sigma_cached-sigma_k))), 0.0_wp, thr=16*epsilon(1.0_wp))
   if (allocated(error)) then
      write(*, '(a,4es24.15)') 'ACP direct response residuals: ', &
         & maxval(abs(dEdcn_cached-dEdcn_k)), maxval(abs(dEdq_cached-dEdq_k)), &
         & maxval(abs(gradient_cached-gradient_k)), maxval(abs(sigma_cached-sigma_k))
      call test_failed(error, &
         & "Direct ACP oracle changed the force or stress response")
      return
   end if

   call cp2k_build_acp_projectors(calc, mol, bcache, acache, lattr, &
      & projectors_changed, projector_cache=projector_cache)
   call projector_cache_counts(projector_cache, cache_hits, cache_misses, cache_bypasses, cache_valid)
   call check(error, cache_misses, 10)
   if (.not.allocated(error)) call check(error, cache_valid)
   if (allocated(error)) then
      call test_failed(error, "ACP projector cache was stale after direct rebuild")
      return
   end if
   call cp2k_build_acp_projectors(calc, mol, bcache, acache, lattr, &
      & projectors_changed, projector_cache=projector_cache)
   call projector_cache_counts(projector_cache, cache_hits, cache_misses, cache_bypasses, cache_valid)
   call check(error, cache_hits, 2)
   if (allocated(error)) then
      call test_failed(error, "ACP projector cache did not recover after bypass")
      return
   end if
end subroutine test_kmesh_acp_h2o_pbc

subroutine test_cp2k_acp_cutoff(error)
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(basis_type) :: diffuse, compact
   real(wp) :: conservative, cutoff, expected, pair_cutoff, valence

   call get_structure(mol, "X23", "urea")
   call new_gxtb_calculator(calc, mol, error, accuracy=acc)
   if (allocated(error)) return

   valence = get_cutoff(calc%bas, acc)
   cutoff = cp2k_get_acp_cutoff(calc, acc)
   conservative = max(valence, get_cutoff(calc%acp%auxbas, acc))
   expected = max(valence, get_pair_cutoff(calc%bas, calc%acp%auxbas, acc))
   cutoff = cp2k_get_acp_cutoff(calc, acc)
   call check(error, cutoff, expected, thr=thr)
   if (allocated(error)) then
      call test_failed(error, &
         & "CP2K ACP cutoff does not cover valence and auxiliary bases")
      return
   end if
   ! For the current H/C/N/O g-xTB parameters the q-vSZP valence basis is
   ! already the more diffuse partner, so the mixed cutoff must not inflate
   ! the established valence range.
   call check(error, cutoff, valence, thr=thr)
   if (allocated(error)) then
      call test_failed(error, &
         & "H/C/N/O ACP cutoff unexpectedly exceeds the valence range")
      return
   end if
   call check(error, cutoff <= conservative)
   if (allocated(error)) then
      call test_failed(error, &
         & "mixed ACP cutoff exceeds the conservative same-basis bound")
      return
   end if

   ! Exercise the opposite exponent ordering independently of the current
   ! parametrization: a diffuse auxiliary basis must enlarge the required
   ! range, but never beyond the conservative same-basis maximum.
   compact%min_alpha = 1.0_wp
   diffuse%min_alpha = 0.25_wp
   valence = get_cutoff(compact, acc)
   pair_cutoff = get_pair_cutoff(compact, diffuse, acc)
   conservative = max(valence, get_cutoff(diffuse, acc))
   call check(error, pair_cutoff > valence .and. pair_cutoff <= conservative)
   if (allocated(error)) then
      call test_failed(error, "mixed Gaussian cutoff does not respect safe bounds")
   end if
end subroutine test_cp2k_acp_cutoff

subroutine test_cp2k_scc_mixer_dimension(error)
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   integer :: ndim

   call get_structure(mol, "MB16-43", "H2")
   call new_gxtb_calculator(calc, mol, error, accuracy=acc)
   if (allocated(error)) return

   ndim = cp2k_get_scc_mixer_dimension(mol, calc%bas, calc%variable_info())
   call check(error, ndim, calc%bas%nsh + 9*mol%nat)
   if (allocated(error)) then
      call test_failed(error, "CP2K mixer dimension includes CP2K-owned AO density")
   end if
end subroutine test_cp2k_scc_mixer_dimension

subroutine test_cp2k_potential_mixer(error)
   type(error_type), allocatable, intent(out) :: error

   type(cp2k_potential_mixer_type) :: mixer
   real(wp), parameter :: rms = sqrt(1.0_wp/3.0_wp)
   real(wp) :: mixed(3)

   call new_cp2k_potential_mixer(mixer, 3, error)
   if (allocated(error)) return

   call mixer%advance([0.0_wp, 0.0_wp, 0.0_wp], mixed, error)
   if (allocated(error)) return
   call check(error, maxval(abs(mixed)), 0.0_wp, thr=thr)
   if (allocated(error)) return
   call check(error, mixer%get_error(), 0.0_wp, thr=thr)
   if (allocated(error)) return

   call mixer%advance([1.0_wp, 0.0_wp, 0.0_wp], mixed, error)
   if (allocated(error)) return
   call check(error, maxval(abs(mixed-[0.2_wp, 0.0_wp, 0.0_wp])), &
      & 0.0_wp, thr=thr)
   if (allocated(error)) return
   call check(error, mixer%get_error(), rms, thr=thr)
   if (allocated(error)) return

   call mixer%advance([0.2_wp, 1.0_wp, 0.0_wp], mixed, error)
   if (allocated(error)) return
   call check(error, maxval(abs(mixed-[0.2_wp, 0.2_wp, 0.0_wp])), &
      & 0.0_wp, thr=thr)
   if (allocated(error)) return
   call check(error, mixer%get_error(), rms, thr=thr)
   if (allocated(error)) return

   ! The fourth Fock build switches to DIIS.  The two orthogonal unit
   ! residuals have equal Pulay weights, so the extrapolated vector is the
   ! arithmetic mean of the two corresponding raw Fock vectors.
   call mixer%advance([1.2_wp, 0.2_wp, 0.0_wp], mixed, error)
   if (allocated(error)) return
   call check(error, maxval(abs(mixed-[0.7_wp, 0.6_wp, 0.0_wp])), &
      & 0.0_wp, thr=thr)
   if (allocated(error)) return
   call check(error, mixer%get_error(), rms, thr=thr)
   if (allocated(error)) return

   call mixer%reset(3, error)
   if (allocated(error)) return
   call mixer%advance([2.0_wp, 3.0_wp, 4.0_wp], mixed, error)
   if (allocated(error)) return
   call check(error, maxval(abs(mixed-[2.0_wp, 3.0_wp, 4.0_wp])), &
      & 0.0_wp, thr=thr)
   if (allocated(error)) then
      call test_failed(error, "CP2K potential mixer reset did not clear its history")
   end if
end subroutine test_cp2k_potential_mixer

subroutine test_cp2k_potential_mixer_uninitialized(error)
   type(error_type), allocatable, intent(out) :: error

   type(cp2k_potential_mixer_type) :: mixer
   real(wp) :: mixed(1)

   call mixer%advance([0.0_wp], mixed, error)
end subroutine test_cp2k_potential_mixer_uninitialized

subroutine test_cp2k_potential_mixer_invalid_dimension(error)
   type(error_type), allocatable, intent(out) :: error

   type(cp2k_potential_mixer_type) :: mixer

   call new_cp2k_potential_mixer(mixer, 0, error)
end subroutine test_cp2k_potential_mixer_invalid_dimension

subroutine test_cp2k_potential_mixer_size_mismatch(error)
   type(error_type), allocatable, intent(out) :: error

   type(cp2k_potential_mixer_type) :: mixer
   real(wp) :: mixed(3)

   call new_cp2k_potential_mixer(mixer, 3, error)
   if (allocated(error)) return
   call mixer%advance([0.0_wp, 0.0_wp], mixed, error)
end subroutine test_cp2k_potential_mixer_size_mismatch

subroutine test_cp2k_potential_mixer_nonfinite(error)
   type(error_type), allocatable, intent(out) :: error

   type(cp2k_potential_mixer_type) :: mixer
   real(wp) :: mixed(1), raw(1)

   call new_cp2k_potential_mixer(mixer, 1, error)
   if (allocated(error)) return
   raw(1) = ieee_value(0.0_wp, ieee_quiet_nan)
   call mixer%advance(raw, mixed, error)
end subroutine test_cp2k_potential_mixer_nonfinite

function contract_acp_image_energy(images, density_translation, density) &
   & result(energy)
   type(acp_image_type), intent(in) :: images
   real(wp), intent(in) :: density_translation(:, :)
   real(wp), intent(in) :: density(:, :, :, :)
   real(wp) :: energy

   integer :: dtr, itr, jtr, spin
   real(wp) :: scale, tol2

   energy = 0.0_wp
   scale = max(1.0_wp, maxval(abs(images%translation)), &
      & maxval(abs(density_translation)))
   tol2 = (2048.0_wp*epsilon(1.0_wp)*scale)**2
   do itr = 1, size(images%translation, 2)
      dtr = 0
      do jtr = 1, size(density_translation, 2)
         if (sum((images%translation(:, itr) &
            & - density_translation(:, jtr))**2) <= tol2) then
            dtr = jtr
            exit
         end if
      end do
      if (dtr == 0) then
         energy = huge(energy)
         return
      end if
      do spin = 1, size(density, 4)
         energy = energy + sum(images%hamiltonian(:, :, itr) &
            & * density(:, :, dtr, spin))
      end do
   end do
end function contract_acp_image_energy

! These two existing fixtures are embedded to avoid adding mstore to CP2K.
subroutine get_structure(mol, collection, name)
   type(structure_type), intent(out) :: mol
   character(len=*), intent(in) :: collection, name
   if (collection == "X23" .and. name == "urea") then
      call urea(mol)
   else if (collection == "MB16-43" .and. name == "H2") then
      call H2(mol)
   else
      error stop "Unknown model-bridge test fixture"
   end if
end subroutine get_structure

! This file is part of mstore.
! SPDX-Identifier: Apache-2.0
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

subroutine urea(self)
   type(structure_type), intent(out) :: self
   integer, parameter :: nat = 16
   character(len=*), parameter :: sym(nat) = [character(len=4)::&
      & "H", "H", "H", "H", "H", "H", "H", "H", "C", "C", "N", "N", "N", "N", &
      & "O", "O"]
   real(wp), parameter :: xyz(3, nat) = reshape([&
      &  2.73311029549004_wp,  7.98560290858176_wp,  2.52731916558693_wp, &
      &  7.98560290858176_wp,  7.77168595812233_wp,  6.32245531053102_wp, &
      &  7.77168595812233_wp,  2.51919334503061_wp,  2.52731916558693_wp, &
      &  2.51919334503061_wp,  2.73311029549004_wp,  6.32245531053102_wp, &
      &  1.52274097774035_wp,  6.77523359083207_wp,  8.53343439213466_wp, &
      &  6.77523359083207_wp,  8.98205527587201_wp,  0.31615111141221_wp, &
      &  8.98205527587201_wp,  3.72956266278029_wp,  8.53343439213466_wp, &
      &  3.72956266278029_wp,  1.52274097774035_wp,  0.31615111141221_wp, &
      &  0.00000000000000_wp,  5.25249261309172_wp,  2.92453950999127_wp, &
      &  5.25249261309172_wp,  0.00000000000000_wp,  5.92523496612667_wp, &
      &  1.54107131713484_wp,  6.79337495765549_wp,  1.60305432044818_wp, &
      &  6.79337495765549_wp,  8.96372493647752_wp,  7.24653118309870_wp, &
      &  8.96372493647752_wp,  3.71142129595688_wp,  1.60305432044818_wp, &
      &  3.71142129595688_wp,  1.54107131713484_wp,  7.24653118309870_wp, &
      &  0.00000000000000_wp,  5.25249261309172_wp,  5.32978239466231_wp, &
      &  5.25249261309172_wp,  0.00000000000000_wp,  3.51999208145564_wp],&
      & shape(xyz))
   real(wp), parameter :: lattice(3, 3) = reshape([&
      & 10.50481250525348_wp,  0.00000000000000_wp,  0.00000000000000_wp, &
      &  0.00000000000000_wp, 10.50481250525348_wp,  0.00000000000000_wp, &
      &  0.00000000000000_wp,  0.00000000000000_wp,  8.84967734421642_wp],&
      & shape(lattice))
   call new(self, sym, xyz, lattice=lattice)
end subroutine urea

subroutine H2(self)
   type(structure_type), intent(out) :: self
   integer, parameter :: nat = 2
   character(len=*), parameter :: sym(nat) = [character(len=4)::&
      & "H", "H"]
   real(wp), parameter :: xyz(3, nat) = reshape([&
      &  0.00000000000000_wp,  0.00000000000000_wp, -0.70252931147690_wp, &
      &  0.00000000000000_wp,  0.00000000000000_wp,  0.70252931147690_wp],&
      & shape(xyz))
   call new(self, sym, xyz)
end subroutine H2

end module tblite_model_bridge_tests
#endif

program tblite_model_bridge_unittest
#if defined(__TBLITE_GXTB_ACCELERATION)
   use, intrinsic :: iso_fortran_env, only: output_unit
   use mctc_env_testing, only: run_testsuite
   use tblite_model_bridge_tests, only: collect_model_bridge
   implicit none
   integer :: failures
   failures = 0
   call run_testsuite(collect_model_bridge, output_unit, failures)
   if (failures /= 0) error stop "CP2K model-bridge tests failed"
   print '(a)', 'CP2K model-bridge integration tests passed'
#else
   implicit none
   print '(a)', 'CP2K model-bridge integration tests skipped: g-xTB response API unavailable'
#endif
end program tblite_model_bridge_unittest
