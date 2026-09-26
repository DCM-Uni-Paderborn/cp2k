! SPDX-Identifier: LGPL-3.0-or-later
! Derived from save_tblite mesh contractions at aad26158b2e83f84b7b7d5b56f361ead4bc02a88.
! CP2K owns mesh transforms, matrix contractions, buffering and their adjoints.
! Radial kernels and their geometry derivatives remain in the physical provider.
#if defined(__TBLITE_GXTB)
module tblite_exchange_contractions
   use mctc_env, only: wp
   use mctc_io, only: structure_type
   use mctc_io_constants, only: pi
   use tblite_exchange_cache, only: exchange_cache
   use tblite_exchange_fock, only: exchange_fock, exchange_bvk_kernel
   use tblite_stream_plan, only: cp2k_stream_plan_type
   use tblite_transform_contract, only: exchange_transform_type, exchange_k_to_r, &
      & exchange_r_to_k, exchange_r_to_k_adjoint, exchange_k_to_r_adjoint
   use tblite_bvk_transform, only: bvk_transform_dense, bvk_transform_separable_dft, &
      bvk_transform_mixed_radix_fft, bvk_separable_k_to_r, bvk_separable_r_to_k, &
      bvk_mixed_radix_k_to_r, bvk_mixed_radix_r_to_k
   implicit none
   private
   public :: get_kfock_kpoint, get_kgrad_kpoint
   public :: get_kfock_kmesh
   public :: get_kgrad_kmesh
   public :: get_kfock_stream_apply
   public :: get_kfock_batch_apply
   public :: get_kgrad_batch_apply
   public :: get_mulliken_derivs_direct, get_bocorr_derivs_direct
contains
pure logical function is_compact_transform_backend(backend) result(compact)
   integer, intent(in) :: backend

   compact = backend == bvk_transform_separable_dft &
      & .or. backend == bvk_transform_mixed_radix_fft
end function is_compact_transform_backend

pure subroutine matmul_complex_contiguous(amat, bmat, cmat)
   complex(wp), contiguous, intent(in) :: amat(:, :), bmat(:, :)
   complex(wp), contiguous, intent(out) :: cmat(:, :)

   cmat = matmul(amat, bmat)
end subroutine matmul_complex_contiguous

subroutine get_KFock_kmesh(self, mol, cache, plan, kernel, kfrac, weights, &
   & density, overlap, fock, vsh, energy, transform_backend, transform)
   class(exchange_fock), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(exchange_cache), intent(in) :: cache
   type(cp2k_stream_plan_type), intent(in), target :: plan
   type(exchange_bvk_kernel), intent(in) :: kernel
   real(wp), intent(in) :: kfrac(:, :), weights(:)
   complex(wp), intent(in) :: density(:, :, :, :), overlap(:, :, :)
   complex(wp), intent(out) :: fock(:, :, :, :)
   real(wp), intent(out) :: vsh(:), energy
   integer, intent(in), optional :: transform_backend
   class(exchange_transform_type), intent(inout), optional :: transform

   integer :: spin, ik, backend

   backend = bvk_transform_dense
   if (present(transform_backend)) backend = transform_backend

   if (present(transform)) then
      call build_KFock_kmesh(self, mol, density, overlap, kfrac, weights, &
         & kernel, cache%g_onsfx, cache%g_onsri, fock, .true., &
         & transform=transform)
   else if (is_compact_transform_backend(backend) .and. &
      & plan%plan_matches(self, mol, kernel, kfrac, weights, &
      & require_phases=.false.) .and. &
      & allocated(plan%bvk_phase_roots) .and. &
      & allocated(plan%bvk_twist_phase)) then
      call build_KFock_kmesh(self, mol, density, overlap, kfrac, weights, &
         & kernel, cache%g_onsfx, cache%g_onsri, fock, .true., &
         & transform_backend=backend, nmesh=plan%bvk_nmesh, &
         & input_to_grid=plan%bvk_input_to_grid, &
         & phase_roots=plan%bvk_phase_roots, &
         & twist_phase=plan%bvk_twist_phase, &
         & transform=transform)
   else if (plan%plan_matches(self, mol, kernel, kfrac, weights)) then
      call build_KFock_kmesh(self, mol, density, overlap, kfrac, weights, &
         & kernel, cache%g_onsfx, cache%g_onsri, fock, .true., &
         & plan%bvk_phase_forward, plan%bvk_phase_inverse, &
         & transform=transform)
   else
      call build_KFock_kmesh(self, mol, density, overlap, kfrac, weights, &
         & kernel, cache%g_onsfx, cache%g_onsri, fock, .true., &
         & transform=transform)
   end if

   energy = 0.0_wp
   do ik = 1, size(weights)
      do spin = 1, size(density, 3)
         energy = energy + 0.5_wp*weights(ik)*real(sum( &
            & conjg(density(:, :, spin, ik))*fock(:, :, spin, ik)), wp)
      end do
   end do

   call get_KFock_kmesh_vsh(self, mol, cache, weights, density, overlap, vsh)

end subroutine get_KFock_kmesh

subroutine get_KFock_kmesh_vsh(self, mol, cache, weights, density, overlap, vsh)
   class(exchange_fock), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(exchange_cache), intent(in) :: cache
   real(wp), intent(in) :: weights(:)
   complex(wp), intent(in) :: density(:, :, :, :), overlap(:, :, :)
   real(wp), intent(out) :: vsh(:)

   integer :: iao, iat, ik, iq, is, ish, izp, spin
   real(wp) :: spin_factor
   real(wp), allocatable :: grad_onsfx(:, :, :), grad_onsri(:, :)
   complex(wp), allocatable :: bq(:, :), bw(:, :), bt0(:, :), bu0(:, :), &
      & bv0(:, :), d0(:, :), dmat(:, :), p0(:, :), x0(:, :), xmat(:, :)
   complex(wp), allocatable :: bgp(:), bgsp(:), bgsps(:), diag2d(:), &
      & diagp(:), diagx(:)

   allocate(grad_onsfx(self%maxsh, self%maxsh, mol%nat), &
      & grad_onsri(self%maxsh, mol%nat), source=0.0_wp)
   allocate(bq(self%nao, self%nao), bw(self%nao, self%nao), &
      & bt0(self%nao, self%nao), bu0(self%nao, self%nao), &
      & bv0(self%nao, self%nao), d0(self%nao, self%nao), &
      & dmat(self%nao, self%nao), p0(self%nao, self%nao), &
      & x0(self%nao, self%nao), xmat(self%nao, self%nao), &
      & bgp(self%nao), bgsp(self%nao), bgsps(self%nao), &
      & diag2d(self%nao), diagp(self%nao), diagx(self%nao), &
      & source=(0.0_wp, 0.0_wp))

   spin_factor = 0.5_wp
   if (size(density, 3) > 1) spin_factor = 1.0_wp

   do spin = 1, size(density, 3)
      p0 = (0.0_wp, 0.0_wp)
      x0 = (0.0_wp, 0.0_wp)
      d0 = (0.0_wp, 0.0_wp)
      do ik = 1, size(weights)
         xmat = matmul(overlap(:, :, ik), density(:, :, spin, ik))
         dmat = 0.5_wp*matmul(xmat, overlap(:, :, ik))
         p0 = p0 + weights(ik)*density(:, :, spin, ik)
         x0 = x0 + weights(ik)*xmat
         d0 = d0 + weights(ik)*dmat
      end do

      bt0 = (0.0_wp, 0.0_wp)
      bu0 = (0.0_wp, 0.0_wp)
      bv0 = (0.0_wp, 0.0_wp)
      bgp = (0.0_wp, 0.0_wp)
      bgsp = (0.0_wp, 0.0_wp)
      bgsps = (0.0_wp, 0.0_wp)
      do ik = 1, size(weights)
         bq = -0.25_wp*spin_factor*weights(ik)*density(:, :, spin, ik)
         bu0 = bu0 + bq
         bw = matmul(bq, conjg(transpose(overlap(:, :, ik))))
         bt0 = bt0 + bw
         bv0 = bv0 + 0.5_wp*matmul( &
            & conjg(transpose(overlap(:, :, ik))), bw)
         do iao = 1, self%nao
            bgp(iao) = bgp(iao) + 0.25_wp*sum( &
               & conjg(overlap(:, iao, ik))*bw(:, iao))
            bgsp(iao) = bgsp(iao) + 0.5_wp*sum( &
               & conjg(overlap(:, iao, ik))*bq(:, iao))
            bgsps(iao) = bgsps(iao) + 0.25_wp*bq(iao, iao)
         end do
      end do

      call onsite_fx_parameter_gradient_complex(mol%nat, mol%id, &
         & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bt0, x0, &
         & 0.5_wp, grad_onsfx, adjoint_src=.true.)
      call onsite_fx_parameter_gradient_complex(mol%nat, mol%id, &
         & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bu0, d0, &
         & 0.5_wp, grad_onsfx)
      call onsite_fx_parameter_gradient_complex(mol%nat, mol%id, &
         & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bv0, p0, &
         & 0.5_wp, grad_onsfx)

      call onsite_ri_parameter_gradient_complex(mol%nat, mol%id, &
         & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bt0, x0, &
         & -2.0_wp, grad_onsri)
      call onsite_ri_parameter_gradient_complex(mol%nat, mol%id, &
         & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bu0, d0, &
         & -2.0_wp, grad_onsri)
      call onsite_ri_parameter_gradient_complex(mol%nat, mol%id, &
         & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bv0, p0, &
         & -2.0_wp, grad_onsri)

      do iao = 1, self%nao
         diagp(iao) = p0(iao, iao)
         diagx(iao) = x0(iao, iao)
         diag2d(iao) = 2.0_wp*d0(iao, iao)
      end do
      call onsite_fx_symv_parameter_gradient_complex(mol%nat, mol%id, &
         & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bgp, &
         & diagp, grad_onsfx)
      call onsite_fx_symv_parameter_gradient_complex(mol%nat, mol%id, &
         & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bgsp, &
         & diagx, grad_onsfx)
      call onsite_fx_symv_parameter_gradient_complex(mol%nat, mol%id, &
         & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bgsps, &
         & diag2d, grad_onsfx)
   end do

   vsh = 0.0_wp
   do iat = 1, mol%nat
      izp = mol%id(iat)
      is = self%ish_at(iat)
      do ish = 1, self%nsh_id(izp)
         iq = is + ish
         vsh(iq) = sum(grad_onsfx(:, :, iat) &
            & *cache%dgdq_onsfx(:, :, iq)) + sum(grad_onsri(:, iat) &
            & *cache%dgdq_onsri(:, iq))
      end do
   end do

end subroutine get_KFock_kmesh_vsh

subroutine build_KFock_kmesh(self, mol, density, overlap, kfrac, weights, &
   & kernel, g_onsfx, g_onsri, fock, use_bvk_kernel, &
   & phase_forward_cached, phase_inverse_cached, transform_backend, nmesh, &
   & input_to_grid, phase_roots, twist_phase, transform)
   class(exchange_fock), intent(in) :: self
   type(structure_type), intent(in) :: mol
   complex(wp), intent(in) :: density(:, :, :, :), overlap(:, :, :)
   real(wp), intent(in) :: kfrac(:, :), weights(:)
   type(exchange_bvk_kernel), intent(in) :: kernel
   real(wp), intent(in) :: g_onsfx(:, :, :), g_onsri(:, :)
   complex(wp), intent(out) :: fock(:, :, :, :)
   logical, intent(in) :: use_bvk_kernel
   complex(wp), contiguous, optional, target, intent(in) :: phase_forward_cached(:, :)
   complex(wp), contiguous, optional, target, intent(in) :: phase_inverse_cached(:, :)
   integer, intent(in), optional :: transform_backend, nmesh(3), &
      & input_to_grid(:)
   complex(wp), intent(in), optional :: phase_roots(:, :), twist_phase(:)
   class(exchange_transform_type), intent(inout), optional :: transform

   integer :: iao, backend, icell, ii, jj, ik, spin
   logical :: use_cached_phases, use_compact
   real(wp) :: angle, spin_factor
   complex(wp) :: phase, tmp
   complex(wp), allocatable, target :: amat(:, :, :), cmat(:, :, :)
   complex(wp), allocatable :: diagP(:), diagSP(:), diagSPS(:)
   complex(wp), allocatable :: gdiagP(:), gdiagSP(:), gdiagSPS(:)
   complex(wp), allocatable, target :: vmat(:, :, :)
   complex(wp), allocatable, target :: source_r(:, :, :), result_r(:, :, :)
   complex(wp), allocatable :: work(:, :)
   complex(wp), allocatable, target :: phase_forward_local(:, :), phase_inverse_local(:, :)
   complex(wp), contiguous, pointer :: phase_forward(:, :), phase_inverse(:, :)

   allocate(amat(self%nao, self%nao, size(weights)), &
      & cmat(self%nao, self%nao, size(weights)), &
      & vmat(self%nao, self%nao, size(weights)), &
      & source_r(self%nao, self%nao, size(kernel%reps, 2)), &
      & result_r(self%nao, self%nao, size(kernel%reps, 2)), &
      & work(self%nao, self%nao), diagP(self%nao), diagSP(self%nao), &
      & diagSPS(self%nao), gdiagP(self%nao), gdiagSP(self%nao), &
      & gdiagSPS(self%nao))
   backend = bvk_transform_dense
   if (present(transform_backend)) backend = transform_backend
   use_cached_phases = .false.
   use_compact = is_compact_transform_backend(backend) .and. present(nmesh) &
      & .and. present(input_to_grid) .and. present(phase_roots) &
      & .and. present(twist_phase)
   if (present(transform)) use_compact = .false.
   if (use_compact .or. present(transform)) then
      ! Compact-grid backends do not materialize an Nk-by-Nk phase table.
      ! Associate harmless scalar placeholders because the dense arguments of
      ! the common kernel entry point are intentionally retained as oracle.
      allocate(phase_forward_local(1, 1), phase_inverse_local(1, 1), &
         & source=(0.0_wp, 0.0_wp))
      phase_forward => phase_forward_local
      phase_inverse => phase_inverse_local
   else if (present(phase_forward_cached) .and. present(phase_inverse_cached)) then
      use_cached_phases = all(shape(phase_forward_cached) == &
         & [size(kernel%reps, 2), size(weights)]) .and. &
         & all(shape(phase_inverse_cached) == &
         & [size(weights), size(kernel%reps, 2)])
   end if
   if (.not.use_compact .and. use_cached_phases) then
      phase_forward => phase_forward_cached
      phase_inverse => phase_inverse_cached
   else if (.not.use_compact .and. .not.present(transform)) then
      allocate(phase_forward_local(size(kernel%reps, 2), size(weights)), &
         & phase_inverse_local(size(weights), size(kernel%reps, 2)))
      phase_forward => phase_forward_local
      phase_inverse => phase_inverse_local

      ! No ordering of the k points is assumed, and a common
      ! Monkhorst-Pack twist is retained explicitly in every phase.
      do ik = 1, size(weights)
         do icell = 1, size(kernel%reps, 2)
            angle = 2.0_wp*pi*dot_product(kfrac(:, ik), &
               & real(kernel%reps(:, icell), wp))
            phase = exp(cmplx(0.0_wp, angle, wp))
            phase_forward(icell, ik) = phase
            phase_inverse(ik, icell) = weights(ik)*conjg(phase)
         end do
      end do
   end if

   spin_factor = 0.5_wp
   if (size(density, 3) > 1) spin_factor = 1.0_wp
   fock = (0.0_wp, 0.0_wp)

   do spin = 1, size(density, 3)
      diagP = (0.0_wp, 0.0_wp)
      diagSP = (0.0_wp, 0.0_wp)
      diagSPS = (0.0_wp, 0.0_wp)
      do ik = 1, size(weights)
         amat(:, :, ik) = matmul(overlap(:, :, ik), density(:, :, spin, ik))
         cmat(:, :, ik) = 0.5_wp*matmul(amat(:, :, ik), overlap(:, :, ik))
         do iao = 1, self%nao
            diagP(iao) = diagP(iao) + weights(ik)*density(iao, iao, spin, ik)
            diagSP(iao) = diagSP(iao) + weights(ik)*amat(iao, iao, ik)
            diagSPS(iao) = diagSPS(iao) + 2.0_wp*weights(ik)*cmat(iao, iao, ik)
         end do
      end do
      call onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, g_onsfx, diagP, gdiagP)
      call onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, g_onsfx, diagSP, gdiagSP)
      call onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, g_onsfx, diagSPS, gdiagSPS)

      call apply_bvk_Kkernel_inplace(self, mol, kernel, phase_inverse, phase_forward, &
         & amat, source_r, result_r, &
         & g_onsfx, g_onsri, 1.0_wp, -4.0_wp, 0.5_wp, -2.0_wp, &
         & use_bvk_kernel, .true., backend, nmesh, input_to_grid, &
         & phase_roots, twist_phase, &
         & transform=transform)
      call apply_bvk_Kkernel_inplace(self, mol, kernel, phase_inverse, phase_forward, &
         & cmat, source_r, result_r, &
         & g_onsfx, g_onsri, 1.0_wp, 0.0_wp, 0.5_wp, -2.0_wp, &
         & use_bvk_kernel, .false., backend, nmesh, input_to_grid, &
         & phase_roots, twist_phase, &
         & transform=transform)

      ! A fixed-spin rank-four section is non-contiguous for UKS.  Copy it
      ! explicitly into the third reusable mesh buffer instead of relying on
      ! a compiler-created pack temporary at the contiguous kernel boundary.
      do ik = 1, size(weights)
         vmat(:, :, ik) = density(:, :, spin, ik)
      end do
      call apply_bvk_Kkernel_inplace(self, mol, kernel, phase_inverse, phase_forward, &
         & vmat, source_r, result_r, g_onsfx, g_onsri, 1.0_wp, 0.0_wp, &
         & 0.5_wp, -2.0_wp, use_bvk_kernel, .false., backend, nmesh, &
         & input_to_grid, phase_roots, twist_phase, &
         & transform=transform)

      do ik = 1, size(weights)
         work = amat(:, :, ik)
         do iao = 1, self%nao
            work(:, iao) = work(:, iao) &
               & + 0.25_wp*gdiagP(iao)*overlap(:, iao, ik)
         end do
         work = work + 0.5_wp*matmul(overlap(:, :, ik), vmat(:, :, ik))
         fock(:, :, spin, ik) = cmat(:, :, ik) &
            & + matmul(work, overlap(:, :, ik))
         do iao = 1, self%nao
            fock(:, iao, spin, ik) = fock(:, iao, spin, ik) &
               & + 0.5_wp*gdiagSP(iao)*overlap(:, iao, ik)
            fock(iao, iao, spin, ik) = fock(iao, iao, spin, ik) &
               & + 0.25_wp*gdiagSPS(iao)
         end do

         do ii = 1, self%nao
            fock(ii, ii, spin, ik) = cmplx(-0.5_wp*spin_factor &
               & *real(fock(ii, ii, spin, ik), wp), 0.0_wp, wp)
            do jj = 1, ii-1
               tmp = -0.25_wp*spin_factor*(fock(jj, ii, spin, ik) &
                  & + conjg(fock(ii, jj, spin, ik)))
               fock(jj, ii, spin, ik) = tmp
               fock(ii, jj, spin, ik) = conjg(tmp)
            end do
         end do
      end do
   end do

end subroutine build_KFock_kmesh

subroutine apply_bvk_Kkernel_realspace(self, mol, kernel, source_r, result_r, &
   & g_onsfx, g_onsri, mulliken_scale, bocorr_scale, onsite_fx_scale, &
   & onsite_ri_scale, use_bvk_kernel, adjoint_src)
   class(exchange_fock), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(exchange_bvk_kernel), intent(in) :: kernel
   complex(wp), intent(in) :: source_r(:, :, :)
   complex(wp), intent(out) :: result_r(:, :, :)
   real(wp), intent(in) :: g_onsfx(:, :, :), g_onsri(:, :)
   real(wp), intent(in) :: mulliken_scale, bocorr_scale
   real(wp), intent(in) :: onsite_fx_scale, onsite_ri_scale
   logical, intent(in) :: use_bvk_kernel, adjoint_src

   integer :: icell, origin

   result_r = (0.0_wp, 0.0_wp)
   origin = 0
   do icell = 1, size(kernel%reps, 2)
      if (all(kernel%reps(:, icell) == 0)) origin = icell
      if (use_bvk_kernel) then
         if (mulliken_scale /= 0.0_wp) &
            call shell_hadamard_add_complex(self%nsh, self%nao_sh, &
               & self%iao_sh, kernel%g_mulliken_r(:, :, icell), &
               & source_r(:, :, icell), mulliken_scale, &
               & result_r(:, :, icell))
         if (bocorr_scale /= 0.0_wp) &
            call atom_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
               & self%nao_sh, self%ish_at, self%iao_sh, &
               & kernel%g_bocorr_r(:, :, icell), source_r(:, :, icell), &
               & bocorr_scale, result_r(:, :, icell))
      end if
   end do

   if (origin > 0) then
      if (onsite_fx_scale /= 0.0_wp) &
         call onsite_fx_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
            & self%nao_sh, self%ish_at, self%iao_sh, g_onsfx, &
            & source_r(:, :, origin), onsite_fx_scale, &
            & result_r(:, :, origin), adjoint_src=adjoint_src)
      if (onsite_ri_scale /= 0.0_wp) &
         call onsite_ri_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
            & self%nao_sh, self%ish_at, self%iao_sh, g_onsri, &
            & source_r(:, :, origin), onsite_ri_scale, &
            & result_r(:, :, origin))
   end if
end subroutine apply_bvk_Kkernel_realspace

subroutine get_KFock_stream_apply(self, mol, cache, kernel, phase_forward, &
   & weights, amat_r, cmat_r, vmat_r, gdiagP, gdiagSP, gdiagSPS, vsh, energy)
   class(exchange_fock), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(exchange_cache), intent(in) :: cache
   type(exchange_bvk_kernel), intent(in) :: kernel
   complex(wp), intent(in) :: phase_forward(:, :)
   real(wp), intent(in) :: weights(:)
   complex(wp), intent(inout) :: amat_r(:, :, :, :), cmat_r(:, :, :, :), &
      & vmat_r(:, :, :, :)
   complex(wp), intent(out) :: gdiagP(:, :), gdiagSP(:, :), gdiagSPS(:, :)
   real(wp), intent(out) :: vsh(:), energy

   integer :: iao, iat, icell, iq, is, ish, izp, jcell, origin, spin
   real(wp) :: spin_factor
   complex(wp) :: coeff, energy_raw
   real(wp), allocatable :: grad_onsfx(:, :, :), grad_onsri(:, :)
   complex(wp), allocatable :: bgp(:), bgsp(:), bgsps(:), bmat(:, :), &
      & bt0(:, :), bu0(:, :), bv0(:, :), diag2d(:), diagp(:), diagx(:), &
      & result_r(:, :, :), result2_r(:, :, :)

   origin = 0
   do icell = 1, size(kernel%reps, 2)
      if (all(kernel%reps(:, icell) == 0)) origin = icell
   end do

   allocate(grad_onsfx(self%maxsh, self%maxsh, mol%nat), &
      & grad_onsri(self%maxsh, mol%nat), source=0.0_wp)
   allocate(bgp(self%nao), bgsp(self%nao), bgsps(self%nao), &
      & bmat(self%nao, self%nao), bt0(self%nao, self%nao), &
      & bu0(self%nao, self%nao), bv0(self%nao, self%nao), &
      & diag2d(self%nao), diagp(self%nao), diagx(self%nao), &
      & result_r(self%nao, self%nao, size(kernel%reps, 2)), &
      & result2_r(self%nao, self%nao, size(kernel%reps, 2)), &
      & source=(0.0_wp, 0.0_wp))

   energy = 0.0_wp
   vsh = 0.0_wp
   spin_factor = 0.5_wp
   if (size(amat_r, 4) > 1) spin_factor = 1.0_wp

   do spin = 1, size(amat_r, 4)
      do iao = 1, self%nao
         diagp(iao) = vmat_r(iao, iao, origin, spin)
         diagx(iao) = amat_r(iao, iao, origin, spin)
         diag2d(iao) = 2.0_wp*cmat_r(iao, iao, origin, spin)
      end do
      call onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, &
         & diagp, gdiagP(:, spin))
      call onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, &
         & diagx, gdiagSP(:, spin))
      call onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, &
         & diag2d, gdiagSPS(:, spin))

      bt0 = -0.25_wp*spin_factor*conjg( &
         & transpose(amat_r(:, :, origin, spin)))
      bu0 = -0.25_wp*spin_factor*vmat_r(:, :, origin, spin)
      bv0 = -0.25_wp*spin_factor*cmat_r(:, :, origin, spin)
      do iao = 1, self%nao
         bgp(iao) = -0.125_wp*spin_factor &
            & *cmat_r(iao, iao, origin, spin)
         bgsp(iao) = -0.125_wp*spin_factor &
            & *amat_r(iao, iao, origin, spin)
         bgsps(iao) = -0.0625_wp*spin_factor &
            & *vmat_r(iao, iao, origin, spin)
      end do

      call onsite_fx_parameter_gradient_complex(mol%nat, mol%id, &
         & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bt0, &
         & amat_r(:, :, origin, spin), 0.5_wp, grad_onsfx, &
         & adjoint_src=.true.)
      call onsite_fx_parameter_gradient_complex(mol%nat, mol%id, &
         & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bu0, &
         & cmat_r(:, :, origin, spin), 0.5_wp, grad_onsfx)
      call onsite_fx_parameter_gradient_complex(mol%nat, mol%id, &
         & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bv0, &
         & vmat_r(:, :, origin, spin), 0.5_wp, grad_onsfx)
      call onsite_ri_parameter_gradient_complex(mol%nat, mol%id, &
         & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bt0, &
         & amat_r(:, :, origin, spin), -2.0_wp, grad_onsri)
      call onsite_ri_parameter_gradient_complex(mol%nat, mol%id, &
         & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bu0, &
         & cmat_r(:, :, origin, spin), -2.0_wp, grad_onsri)
      call onsite_ri_parameter_gradient_complex(mol%nat, mol%id, &
         & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bv0, &
         & vmat_r(:, :, origin, spin), -2.0_wp, grad_onsri)
      call onsite_fx_symv_parameter_gradient_complex(mol%nat, mol%id, &
         & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bgp, &
         & diagp, grad_onsfx)
      call onsite_fx_symv_parameter_gradient_complex(mol%nat, mol%id, &
         & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bgsp, &
         & diagx, grad_onsfx)
      call onsite_fx_symv_parameter_gradient_complex(mol%nat, mol%id, &
         & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bgsps, &
         & diag2d, grad_onsfx)

      energy_raw = (0.0_wp, 0.0_wp)
      do iao = 1, self%nao
         energy_raw = energy_raw &
            & +0.5_wp*cmat_r(iao, iao, origin, spin)*gdiagP(iao, spin) &
            & +0.5_wp*conjg(amat_r(iao, iao, origin, spin)) &
            & *gdiagSP(iao, spin) &
            & +0.25_wp*vmat_r(iao, iao, origin, spin)*gdiagSPS(iao, spin)
      end do

      call apply_bvk_Kkernel_realspace(self, mol, kernel, &
         & amat_r(:, :, :, spin), result_r, cache%g_onsfx, cache%g_onsri, &
         & 1.0_wp, -4.0_wp, 0.5_wp, -2.0_wp, .true., .true.)
      do icell = 1, size(kernel%reps, 2)
         bmat = (0.0_wp, 0.0_wp)
         do jcell = 1, size(kernel%reps, 2)
            coeff = sum(weights*conjg(phase_forward(icell, :) &
               & *phase_forward(jcell, :)))
            bmat = bmat + coeff*conjg( &
               & transpose(amat_r(:, :, jcell, spin)))
         end do
         energy_raw = energy_raw + sum(conjg(bmat)*result_r(:, :, icell))
      end do
      amat_r(:, :, :, spin) = result_r

      call apply_bvk_Kkernel_realspace(self, mol, kernel, &
         & cmat_r(:, :, :, spin), result_r, cache%g_onsfx, cache%g_onsri, &
         & 1.0_wp, 0.0_wp, 0.5_wp, -2.0_wp, .true., .false.)
      energy_raw = energy_raw + sum(conjg(vmat_r(:, :, :, spin))*result_r)

      call apply_bvk_Kkernel_realspace(self, mol, kernel, &
         & vmat_r(:, :, :, spin), result2_r, cache%g_onsfx, &
         & cache%g_onsri, 1.0_wp, 0.0_wp, 0.5_wp, -2.0_wp, .true., .false.)
      energy_raw = energy_raw + sum(conjg(cmat_r(:, :, :, spin)) &
         & *result2_r)
      vmat_r(:, :, :, spin) = result2_r
      cmat_r(:, :, :, spin) = result_r
      energy = energy - 0.25_wp*spin_factor*real(energy_raw, wp)
   end do

   do iat = 1, mol%nat
      izp = mol%id(iat)
      is = self%ish_at(iat)
      do ish = 1, self%nsh_id(izp)
         iq = is + ish
         vsh(iq) = sum(grad_onsfx(:, :, iat) &
            & *cache%dgdq_onsfx(:, :, iq)) + sum(grad_onsri(:, iat) &
            & *cache%dgdq_onsri(:, iq))
      end do
   end do
end subroutine get_KFock_stream_apply

subroutine apply_bvk_Kkernel_image(self, mol, kernel, icell, source, result, &
   & g_onsfx, g_onsri, mulliken_scale, bocorr_scale, onsite_fx_scale, &
   & onsite_ri_scale, use_bvk_kernel, adjoint_src)
   class(exchange_fock), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(exchange_bvk_kernel), intent(in) :: kernel
   integer, intent(in) :: icell
   complex(wp), intent(in) :: source(:, :)
   complex(wp), intent(out) :: result(:, :)
   real(wp), intent(in) :: g_onsfx(:, :, :), g_onsri(:, :)
   real(wp), intent(in) :: mulliken_scale, bocorr_scale
   real(wp), intent(in) :: onsite_fx_scale, onsite_ri_scale
   logical, intent(in) :: use_bvk_kernel, adjoint_src

   result = (0.0_wp, 0.0_wp)
   if (use_bvk_kernel) then
      if (mulliken_scale /= 0.0_wp) &
         call shell_hadamard_add_complex(self%nsh, self%nao_sh, &
            & self%iao_sh, kernel%g_mulliken_r(:, :, icell), source, &
            & mulliken_scale, result)
      if (bocorr_scale /= 0.0_wp) &
         call atom_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
            & self%nao_sh, self%ish_at, self%iao_sh, &
            & kernel%g_bocorr_r(:, :, icell), source, bocorr_scale, result)
   end if
   if (all(kernel%reps(:, icell) == 0)) then
      if (onsite_fx_scale /= 0.0_wp) &
         call onsite_fx_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
            & self%nao_sh, self%ish_at, self%iao_sh, g_onsfx, source, &
            & onsite_fx_scale, result, adjoint_src=adjoint_src)
      if (onsite_ri_scale /= 0.0_wp) &
         call onsite_ri_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
            & self%nao_sh, self%ish_at, self%iao_sh, g_onsri, source, &
            & onsite_ri_scale, result)
   end if
end subroutine apply_bvk_Kkernel_image

subroutine get_KFock_batch_apply(self, mol, cache, kernel, image_index, &
   & inverse_coeff, amat_r, amat_partner_r, cmat_r, vmat_r, gdiagP, &
   & gdiagSP, gdiagSPS, vsh, energy)
   class(exchange_fock), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(exchange_cache), intent(in) :: cache
   type(exchange_bvk_kernel), intent(in) :: kernel
   integer, intent(in) :: image_index(:)
   complex(wp), intent(in) :: inverse_coeff(:)
   complex(wp), intent(inout) :: amat_r(:, :, :, :), &
      & amat_partner_r(:, :, :, :), cmat_r(:, :, :, :), &
      & vmat_r(:, :, :, :)
   complex(wp), intent(inout) :: gdiagP(:, :), gdiagSP(:, :), &
      & gdiagSPS(:, :)
   real(wp), intent(out) :: vsh(:), energy

   integer :: iao, iat, ibatch, icell, iq, is, ish, izp, origin_local, spin
   real(wp) :: spin_factor
   complex(wp) :: energy_raw
   real(wp), allocatable :: grad_onsfx(:, :, :), grad_onsri(:, :)
   complex(wp), allocatable :: bgp(:), bgsp(:), bgsps(:), bmat(:, :), &
      & bt0(:, :), bu0(:, :), bv0(:, :), diag2d(:), diagp(:), diagx(:), &
      & work(:, :)

   origin_local = 0
   do ibatch = 1, size(image_index)
      if (all(kernel%reps(:, image_index(ibatch)) == 0)) &
         & origin_local = ibatch
   end do

   allocate(grad_onsfx(self%maxsh, self%maxsh, mol%nat), &
      & grad_onsri(self%maxsh, mol%nat), source=0.0_wp)
   allocate(bgp(self%nao), bgsp(self%nao), bgsps(self%nao), &
      & bmat(self%nao, self%nao), bt0(self%nao, self%nao), &
      & bu0(self%nao, self%nao), bv0(self%nao, self%nao), &
      & diag2d(self%nao), diagp(self%nao), diagx(self%nao), &
      & work(self%nao, self%nao), source=(0.0_wp, 0.0_wp))

   energy = 0.0_wp
   vsh = 0.0_wp
   spin_factor = 0.5_wp
   if (size(amat_r, 4) > 1) spin_factor = 1.0_wp

   do spin = 1, size(amat_r, 4)
      energy_raw = (0.0_wp, 0.0_wp)
      if (origin_local > 0) then
         do iao = 1, self%nao
            diagp(iao) = vmat_r(iao, iao, origin_local, spin)
            diagx(iao) = amat_r(iao, iao, origin_local, spin)
            diag2d(iao) = 2.0_wp*cmat_r(iao, iao, origin_local, spin)
         end do
         call onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
            & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, &
            & diagp, gdiagP(:, spin))
         call onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
            & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, &
            & diagx, gdiagSP(:, spin))
         call onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
            & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, &
            & diag2d, gdiagSPS(:, spin))

         bt0 = -0.25_wp*spin_factor*conjg( &
            & transpose(amat_r(:, :, origin_local, spin)))
         bu0 = -0.25_wp*spin_factor*vmat_r(:, :, origin_local, spin)
         bv0 = -0.25_wp*spin_factor*cmat_r(:, :, origin_local, spin)
         do iao = 1, self%nao
            bgp(iao) = -0.125_wp*spin_factor &
               & *cmat_r(iao, iao, origin_local, spin)
            bgsp(iao) = -0.125_wp*spin_factor &
               & *amat_r(iao, iao, origin_local, spin)
            bgsps(iao) = -0.0625_wp*spin_factor &
               & *vmat_r(iao, iao, origin_local, spin)
         end do

         call onsite_fx_parameter_gradient_complex(mol%nat, mol%id, &
            & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bt0, &
            & amat_r(:, :, origin_local, spin), 0.5_wp, grad_onsfx, &
            & adjoint_src=.true.)
         call onsite_fx_parameter_gradient_complex(mol%nat, mol%id, &
            & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bu0, &
            & cmat_r(:, :, origin_local, spin), 0.5_wp, grad_onsfx)
         call onsite_fx_parameter_gradient_complex(mol%nat, mol%id, &
            & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bv0, &
            & vmat_r(:, :, origin_local, spin), 0.5_wp, grad_onsfx)
         call onsite_ri_parameter_gradient_complex(mol%nat, mol%id, &
            & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bt0, &
            & amat_r(:, :, origin_local, spin), -2.0_wp, grad_onsri)
         call onsite_ri_parameter_gradient_complex(mol%nat, mol%id, &
            & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bu0, &
            & cmat_r(:, :, origin_local, spin), -2.0_wp, grad_onsri)
         call onsite_ri_parameter_gradient_complex(mol%nat, mol%id, &
            & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bv0, &
            & vmat_r(:, :, origin_local, spin), -2.0_wp, grad_onsri)
         call onsite_fx_symv_parameter_gradient_complex(mol%nat, mol%id, &
            & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bgp, &
            & diagp, grad_onsfx)
         call onsite_fx_symv_parameter_gradient_complex(mol%nat, mol%id, &
            & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bgsp, &
            & diagx, grad_onsfx)
         call onsite_fx_symv_parameter_gradient_complex(mol%nat, mol%id, &
            & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bgsps, &
            & diag2d, grad_onsfx)

         do iao = 1, self%nao
            energy_raw = energy_raw &
               & +0.5_wp*cmat_r(iao, iao, origin_local, spin) &
               & *gdiagP(iao, spin) &
               & +0.5_wp*conjg(amat_r(iao, iao, origin_local, spin)) &
               & *gdiagSP(iao, spin) &
               & +0.25_wp*vmat_r(iao, iao, origin_local, spin) &
               & *gdiagSPS(iao, spin)
         end do
      end if

      do ibatch = 1, size(image_index)
         icell = image_index(ibatch)
         call apply_bvk_Kkernel_image(self, mol, kernel, icell, &
            & amat_r(:, :, ibatch, spin), work, cache%g_onsfx, &
            & cache%g_onsri, 1.0_wp, -4.0_wp, 0.5_wp, -2.0_wp, &
            & .true., .true.)
         bmat = inverse_coeff(ibatch)*conjg( &
            & transpose(amat_partner_r(:, :, ibatch, spin)))
         energy_raw = energy_raw + sum(conjg(bmat)*work)
         amat_r(:, :, ibatch, spin) = work

         call apply_bvk_Kkernel_image(self, mol, kernel, icell, &
            & cmat_r(:, :, ibatch, spin), work, cache%g_onsfx, &
            & cache%g_onsri, 1.0_wp, 0.0_wp, 0.5_wp, -2.0_wp, &
            & .true., .false.)
         energy_raw = energy_raw &
            & +sum(conjg(vmat_r(:, :, ibatch, spin))*work)
         amat_partner_r(:, :, ibatch, spin) = work

         call apply_bvk_Kkernel_image(self, mol, kernel, icell, &
            & vmat_r(:, :, ibatch, spin), work, cache%g_onsfx, &
            & cache%g_onsri, 1.0_wp, 0.0_wp, 0.5_wp, -2.0_wp, &
            & .true., .false.)
         energy_raw = energy_raw &
            & +sum(conjg(cmat_r(:, :, ibatch, spin))*work)
         cmat_r(:, :, ibatch, spin) = &
            & amat_partner_r(:, :, ibatch, spin)
         vmat_r(:, :, ibatch, spin) = work
      end do
      energy = energy - 0.25_wp*spin_factor*real(energy_raw, wp)
   end do

   do iat = 1, mol%nat
      izp = mol%id(iat)
      is = self%ish_at(iat)
      do ish = 1, self%nsh_id(izp)
         iq = is + ish
         vsh(iq) = sum(grad_onsfx(:, :, iat) &
            & *cache%dgdq_onsfx(:, :, iq)) + sum(grad_onsri(:, iat) &
            & *cache%dgdq_onsri(:, iq))
      end do
   end do
end subroutine get_KFock_batch_apply

subroutine reverse_bvk_Kkernel_image(self, mol, kernel, icell, source, &
   & result_adjoint, g_onsfx, g_onsri, mulliken_scale, bocorr_scale, &
   & onsite_fx_scale, onsite_ri_scale, use_bvk_kernel, adjoint_src, &
   & source_adjoint, mulliken_grad, bocorr_grad)
   class(exchange_fock), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(exchange_bvk_kernel), intent(in) :: kernel
   integer, intent(in) :: icell
   complex(wp), intent(in) :: source(:, :), result_adjoint(:, :)
   real(wp), intent(in) :: g_onsfx(:, :, :), g_onsri(:, :)
   real(wp), intent(in) :: mulliken_scale, bocorr_scale
   real(wp), intent(in) :: onsite_fx_scale, onsite_ri_scale
   logical, intent(in) :: use_bvk_kernel, adjoint_src
   complex(wp), intent(out) :: source_adjoint(:, :)
   real(wp), intent(inout) :: mulliken_grad(:, :), bocorr_grad(:, :)

   source_adjoint = (0.0_wp, 0.0_wp)
   if (use_bvk_kernel) then
      if (mulliken_scale /= 0.0_wp) then
         call shell_parameter_gradient_complex(self%nsh, self%nao_sh, &
            & self%iao_sh, result_adjoint, source, mulliken_scale, &
            & mulliken_grad)
         call shell_hadamard_add_complex(self%nsh, self%nao_sh, &
            & self%iao_sh, kernel%g_mulliken_r(:, :, icell), &
            & result_adjoint, mulliken_scale, source_adjoint)
      end if
      if (bocorr_scale /= 0.0_wp) then
         call atom_parameter_gradient_complex(mol%nat, mol%id, &
            & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, &
            & result_adjoint, source, bocorr_scale, bocorr_grad)
         call atom_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
            & self%nao_sh, self%ish_at, self%iao_sh, &
            & kernel%g_bocorr_r(:, :, icell), result_adjoint, &
            & bocorr_scale, source_adjoint)
      end if
   end if
   if (all(kernel%reps(:, icell) == 0)) then
      if (onsite_fx_scale /= 0.0_wp) &
         call onsite_fx_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
            & self%nao_sh, self%ish_at, self%iao_sh, g_onsfx, &
            & result_adjoint, onsite_fx_scale, source_adjoint, &
            & adjoint_src=adjoint_src)
      if (onsite_ri_scale /= 0.0_wp) &
         call onsite_ri_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
            & self%nao_sh, self%ish_at, self%iao_sh, g_onsri, &
            & result_adjoint, onsite_ri_scale, source_adjoint)
   end if
end subroutine reverse_bvk_Kkernel_image

subroutine get_KGrad_batch_apply(self, mol, cache, kernel, image_index, &
   & xmat_r, dmat_r, pmat_r, bmat_r, amat_r, vseed_r, gdiag_p, gdiag_sp, &
   & gdiag_sps, bdiag_sp, bdiag_sps, gradient, sigma)
   class(exchange_fock), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(exchange_cache), intent(in) :: cache
   type(exchange_bvk_kernel), intent(in) :: kernel
   integer, intent(in) :: image_index(:)
   complex(wp), intent(inout) :: xmat_r(:, :, :, :), &
      & dmat_r(:, :, :, :), pmat_r(:, :, :, :), bmat_r(:, :, :, :)
   complex(wp), intent(in) :: amat_r(:, :, :, :), vseed_r(:, :, :, :)
   complex(wp), intent(inout) :: gdiag_p(:, :), gdiag_sp(:, :), &
      & gdiag_sps(:, :), bdiag_sp(:, :), bdiag_sps(:, :)
   real(wp), intent(inout) :: gradient(:, :), sigma(:, :)

   integer :: ibatch, icell, spin
   real(wp), allocatable :: bocorr_grad(:, :, :), mulliken_grad(:, :, :)
   complex(wp), allocatable :: source_adjoint(:, :), work(:, :)

   allocate(mulliken_grad(self%nao, self%nao, size(image_index)), &
      & bocorr_grad(self%nao, self%nao, size(image_index)), source=0.0_wp)
   allocate(source_adjoint(self%nao, self%nao), &
      & work(self%nao, self%nao), source=(0.0_wp, 0.0_wp))

   do spin = 1, size(xmat_r, 4)
      do ibatch = 1, size(image_index)
         icell = image_index(ibatch)

         call apply_bvk_Kkernel_image(self, mol, kernel, icell, &
            & xmat_r(:, :, ibatch, spin), work, cache%g_onsfx, &
            & cache%g_onsri, 1.0_wp, -4.0_wp, 0.5_wp, -2.0_wp, &
            & .true., .true.)
         call reverse_bvk_Kkernel_image(self, mol, kernel, icell, &
            & xmat_r(:, :, ibatch, spin), bmat_r(:, :, ibatch, spin), &
            & cache%g_onsfx, cache%g_onsri, 1.0_wp, -4.0_wp, 0.5_wp, &
            & -2.0_wp, .true., .true., source_adjoint, &
            & mulliken_grad(:, :, ibatch), bocorr_grad(:, :, ibatch))
         xmat_r(:, :, ibatch, spin) = work
         bmat_r(:, :, ibatch, spin) = source_adjoint

         call reverse_bvk_Kkernel_image(self, mol, kernel, icell, &
            & dmat_r(:, :, ibatch, spin), amat_r(:, :, ibatch, spin), &
            & cache%g_onsfx, cache%g_onsri, 1.0_wp, 0.0_wp, 0.5_wp, &
            & -2.0_wp, .true., .false., source_adjoint, &
            & mulliken_grad(:, :, ibatch), bocorr_grad(:, :, ibatch))
         dmat_r(:, :, ibatch, spin) = source_adjoint

         call apply_bvk_Kkernel_image(self, mol, kernel, icell, &
            & pmat_r(:, :, ibatch, spin), work, cache%g_onsfx, &
            & cache%g_onsri, 1.0_wp, 0.0_wp, 0.5_wp, -2.0_wp, &
            & .true., .false.)
         call reverse_bvk_Kkernel_image(self, mol, kernel, icell, &
            & pmat_r(:, :, ibatch, spin), vseed_r(:, :, ibatch, spin), &
            & cache%g_onsfx, cache%g_onsri, 1.0_wp, 0.0_wp, 0.5_wp, &
            & -2.0_wp, .true., .false., source_adjoint, &
            & mulliken_grad(:, :, ibatch), bocorr_grad(:, :, ibatch))
         pmat_r(:, :, ibatch, spin) = work
      end do
   end do

   ! The three BZ-averaged diagonal maps belong to the unique origin image.
   if (any([(all(kernel%reps(:, image_index(ibatch)) == 0), &
      & ibatch=1, size(image_index))])) then
      do spin = 1, size(xmat_r, 4)
         call onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
            & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, &
            & gdiag_p(:, spin), work(:, 1))
         gdiag_p(:, spin) = work(:, 1)
         call onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
            & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, &
            & gdiag_sp(:, spin), work(:, 1))
         gdiag_sp(:, spin) = work(:, 1)
         call onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
            & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, &
            & gdiag_sps(:, spin), work(:, 1))
         gdiag_sps(:, spin) = work(:, 1)
         call onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
            & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, &
            & bdiag_sp(:, spin), work(:, 1))
         bdiag_sp(:, spin) = work(:, 1)
         call onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
            & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, &
            & bdiag_sps(:, spin), work(:, 1))
         bdiag_sps(:, spin) = work(:, 1)
      end do
   end if

   call get_kernel_image_batch_derivs(self, mol, kernel, image_index, &
      & mulliken_grad, bocorr_grad, gradient, sigma)
end subroutine get_KGrad_batch_apply

!> Reuse the public physical derivative with only the requested image labels.
!> No full response mesh or radial-kernel payload is copied. The provider
!> retains its Wigner--Seitz weights and their force/strain derivatives.
subroutine get_kernel_image_batch_derivs(self, mol, kernel, image_index, &
   & mulliken_grad, bocorr_grad, gradient, sigma)
   class(exchange_fock), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(exchange_bvk_kernel), intent(in) :: kernel
   integer, intent(in) :: image_index(:)
   real(wp), intent(in) :: mulliken_grad(:, :, :), bocorr_grad(:, :, :)
   real(wp), intent(inout) :: gradient(:, :), sigma(:, :)
   type(exchange_bvk_kernel) :: batch

   batch%nmesh = kernel%nmesh
   ! Materialize the vector-subscript section in an explicitly bounded array
   ! before passing its allocatable descriptor across the provider boundary.
   allocate(batch%reps(3, size(image_index)))
   batch%reps = kernel%reps(:, image_index)
   call self%get_bvk_Kmatrix_derivs(mol, batch, mulliken_grad, &
      & bocorr_grad, gradient, sigma)
end subroutine get_kernel_image_batch_derivs

subroutine apply_bvk_Kkernel_inplace(self, mol, kernel, phase_inverse, phase_forward, &
   & mesh_k, source_r, result_r, g_onsfx, g_onsri, mulliken_scale, &
   & bocorr_scale, onsite_fx_scale, onsite_ri_scale, use_bvk_kernel, &
   & adjoint_src, transform_backend, nmesh, input_to_grid, phase_roots, &
   & twist_phase, transform)
   class(exchange_fock), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(exchange_bvk_kernel), intent(in) :: kernel
   complex(wp), contiguous, intent(in) :: phase_inverse(:, :), phase_forward(:, :)
   complex(wp), contiguous, target, intent(inout) :: mesh_k(:, :, :)
   complex(wp), contiguous, target, intent(inout) :: source_r(:, :, :)
   complex(wp), contiguous, target, intent(inout) :: result_r(:, :, :)
   real(wp), intent(in) :: g_onsfx(:, :, :), g_onsri(:, :)
   real(wp), intent(in) :: mulliken_scale, bocorr_scale
   real(wp), intent(in) :: onsite_fx_scale, onsite_ri_scale
   logical, intent(in) :: use_bvk_kernel, adjoint_src
   integer, intent(in), optional :: transform_backend, nmesh(3), &
      & input_to_grid(:)
   complex(wp), intent(in), optional :: phase_roots(:, :), twist_phase(:)
   class(exchange_transform_type), intent(inout), optional :: transform

   integer :: backend, icell, npair, origin
   complex(wp), pointer, contiguous :: mesh_k_2d(:, :), source_r_2d(:, :), &
      & result_r_2d(:, :)

   npair = self%nao*self%nao
   mesh_k_2d(1:npair, 1:size(mesh_k, 3)) => mesh_k
   source_r_2d(1:npair, 1:size(source_r, 3)) => source_r
   result_r_2d(1:npair, 1:size(result_r, 3)) => result_r

   backend = bvk_transform_dense
   if (present(transform_backend)) backend = transform_backend

   ! Treat all AO pairs as independent rows of one transform.  The dense path
   ! uses one complex matrix multiplication; compact paths consume the cached
   ! grid ordering, roots, and twist without materializing phase tables.
   if (present(transform)) then
      call transform%apply(mesh_k_2d, source_r_2d, exchange_k_to_r)
   else if (is_compact_transform_backend(backend) .and. present(nmesh) &
      & .and. present(input_to_grid) .and. present(phase_roots) &
      & .and. present(twist_phase)) then
      select case (backend)
      case (bvk_transform_mixed_radix_fft)
         call bvk_mixed_radix_k_to_r(mesh_k_2d, source_r_2d, nmesh, &
            & input_to_grid, kernel%reps, phase_roots, twist_phase, .true.)
      case default
         call bvk_separable_k_to_r(mesh_k_2d, source_r_2d, nmesh, &
            & input_to_grid, kernel%reps, phase_roots, twist_phase, .true.)
      end select
   else
      call matmul_complex_contiguous(mesh_k_2d, phase_inverse, source_r_2d)
   end if
   result_r = (0.0_wp, 0.0_wp)

   origin = 0
   do icell = 1, size(kernel%reps, 2)
      if (all(kernel%reps(:, icell) == 0)) origin = icell
   end do

   if (use_bvk_kernel) then
      do icell = 1, size(kernel%reps, 2)
         if (mulliken_scale /= 0.0_wp) &
            call shell_hadamard_add_complex(self%nsh, self%nao_sh, &
               & self%iao_sh, kernel%g_mulliken_r(:, :, icell), &
               & source_r(:, :, icell), mulliken_scale, result_r(:, :, icell))
         if (bocorr_scale /= 0.0_wp) &
            call atom_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
               & self%nao_sh, self%ish_at, self%iao_sh, &
               & kernel%g_bocorr_r(:, :, icell), source_r(:, :, icell), &
               & bocorr_scale, result_r(:, :, icell))
      end do
   end if

   if (origin > 0) then
      if (onsite_fx_scale /= 0.0_wp) &
         call onsite_fx_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
            & self%nao_sh, self%ish_at, self%iao_sh, g_onsfx, &
            & source_r(:, :, origin), onsite_fx_scale, result_r(:, :, origin), &
            & adjoint_src=adjoint_src)
      if (onsite_ri_scale /= 0.0_wp) &
         call onsite_ri_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
            & self%nao_sh, self%ish_at, self%iao_sh, g_onsri, &
            & source_r(:, :, origin), onsite_ri_scale, result_r(:, :, origin))
   end if

   ! Reverse the same batched AO-pair map with the selected backend.
   if (present(transform)) then
      call transform%apply(result_r_2d, mesh_k_2d, exchange_r_to_k)
   else if (is_compact_transform_backend(backend) .and. present(nmesh) &
      & .and. present(input_to_grid) .and. present(phase_roots) &
      & .and. present(twist_phase)) then
      select case (backend)
      case (bvk_transform_mixed_radix_fft)
         call bvk_mixed_radix_r_to_k(result_r_2d, mesh_k_2d, nmesh, &
            & input_to_grid, kernel%reps, phase_roots, twist_phase)
      case default
         call bvk_separable_r_to_k(result_r_2d, mesh_k_2d, nmesh, &
            & input_to_grid, kernel%reps, phase_roots, twist_phase)
      end select
   else
      call matmul_complex_contiguous(result_r_2d, phase_forward, mesh_k_2d)
   end if

end subroutine apply_bvk_Kkernel_inplace

subroutine apply_bvk_Kkernel(self, mol, kernel, phase_inverse, phase_forward, &
   & source_k, g_onsfx, g_onsri, mulliken_scale, bocorr_scale, &
   & onsite_fx_scale, onsite_ri_scale, use_bvk_kernel, adjoint_src, result_k, &
   & transform_backend, nmesh, input_to_grid, phase_roots, twist_phase, transform)
   class(exchange_fock), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(exchange_bvk_kernel), intent(in) :: kernel
   complex(wp), intent(in) :: phase_inverse(:, :), phase_forward(:, :)
   complex(wp), contiguous, intent(in) :: source_k(:, :, :)
   real(wp), intent(in) :: g_onsfx(:, :, :), g_onsri(:, :)
   real(wp), intent(in) :: mulliken_scale, bocorr_scale
   real(wp), intent(in) :: onsite_fx_scale, onsite_ri_scale
   logical, intent(in) :: use_bvk_kernel, adjoint_src
   complex(wp), contiguous, target, intent(out) :: result_k(:, :, :)
   integer, intent(in), optional :: transform_backend, nmesh(3), &
      & input_to_grid(:)
   complex(wp), intent(in), optional :: phase_roots(:, :), twist_phase(:)
   class(exchange_transform_type), intent(inout), optional :: transform

   complex(wp), allocatable, target :: result_r(:, :, :), source_r(:, :, :)

   allocate(source_r(self%nao, self%nao, size(kernel%reps, 2)), &
      & result_r(self%nao, self%nao, size(kernel%reps, 2)))
   result_k = source_k
   call apply_bvk_Kkernel_inplace(self, mol, kernel, phase_inverse, &
      & phase_forward, result_k, source_r, result_r, g_onsfx, g_onsri, &
      & mulliken_scale, bocorr_scale, onsite_fx_scale, onsite_ri_scale, &
      & use_bvk_kernel, adjoint_src, transform_backend, nmesh, input_to_grid, &
      & phase_roots, twist_phase, &
         & transform=transform)

end subroutine apply_bvk_Kkernel

subroutine reverse_bvk_Kkernel(self, mol, kernel, weights, phase_inverse, &
   & phase_forward, source_k, result_adjoint_k, g_onsfx, g_onsri, &
   & mulliken_scale, bocorr_scale, onsite_fx_scale, onsite_ri_scale, &
   & use_bvk_kernel, adjoint_src, source_adjoint_k, mulliken_grad_r, &
   & bocorr_grad_r, transform_backend, nmesh, input_to_grid, phase_roots, &
   & twist_phase, transform)
   class(exchange_fock), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(exchange_bvk_kernel), intent(in) :: kernel
   real(wp), intent(in) :: weights(:)
   complex(wp), intent(in) :: phase_inverse(:, :), phase_forward(:, :)
   complex(wp), contiguous, target, intent(in) :: source_k(:, :, :)
   complex(wp), contiguous, target, intent(in) :: result_adjoint_k(:, :, :)
   real(wp), intent(in) :: g_onsfx(:, :, :), g_onsri(:, :)
   real(wp), intent(in) :: mulliken_scale, bocorr_scale
   real(wp), intent(in) :: onsite_fx_scale, onsite_ri_scale
   logical, intent(in) :: use_bvk_kernel, adjoint_src
   complex(wp), contiguous, target, intent(out) :: source_adjoint_k(:, :, :)
   real(wp), intent(inout) :: mulliken_grad_r(:, :, :), bocorr_grad_r(:, :, :)
   integer, intent(in), optional :: transform_backend, nmesh(3), &
      & input_to_grid(:)
   complex(wp), intent(in), optional :: phase_roots(:, :), twist_phase(:)
   class(exchange_transform_type), intent(inout), optional :: transform

   integer :: backend, icell, ik, npair, origin
   complex(wp), allocatable, target :: result_adjoint_r(:, :, :), &
      & source_adjoint_r(:, :, :), source_r(:, :, :)
   complex(wp), pointer :: result_adjoint_k_2d(:, :), &
      & result_adjoint_r_2d(:, :), source_adjoint_k_2d(:, :), &
      & source_adjoint_r_2d(:, :), source_k_2d(:, :), source_r_2d(:, :)

   allocate(source_r(self%nao, self%nao, size(kernel%reps, 2)), &
      & result_adjoint_r(self%nao, self%nao, size(kernel%reps, 2)), &
      & source_adjoint_r(self%nao, self%nao, size(kernel%reps, 2)), &
      & source=(0.0_wp, 0.0_wp))
   npair = self%nao*self%nao
   source_k_2d(1:npair, 1:size(source_k, 3)) => source_k
   source_r_2d(1:npair, 1:size(source_r, 3)) => source_r
   result_adjoint_k_2d(1:npair, 1:size(result_adjoint_k, 3)) &
      & => result_adjoint_k
   result_adjoint_r_2d(1:npair, 1:size(result_adjoint_r, 3)) &
      & => result_adjoint_r
   source_adjoint_r_2d(1:npair, 1:size(source_adjoint_r, 3)) &
      & => source_adjoint_r
   source_adjoint_k_2d(1:npair, 1:size(source_adjoint_k, 3)) &
      & => source_adjoint_k

   backend = bvk_transform_dense
   if (present(transform_backend)) backend = transform_backend
   if (present(transform)) then
      call transform%apply(source_k_2d, source_r_2d, exchange_k_to_r)
      call transform%apply(result_adjoint_k_2d, result_adjoint_r_2d, exchange_r_to_k_adjoint)
   else if (is_compact_transform_backend(backend) .and. present(nmesh) &
      & .and. present(input_to_grid) .and. present(phase_roots) &
      & .and. present(twist_phase)) then
      select case (backend)
      case (bvk_transform_mixed_radix_fft)
         call bvk_mixed_radix_k_to_r(source_k_2d, source_r_2d, nmesh, &
            & input_to_grid, kernel%reps, phase_roots, twist_phase, .true.)
         call bvk_mixed_radix_k_to_r(result_adjoint_k_2d, &
            & result_adjoint_r_2d, nmesh, input_to_grid, kernel%reps, &
            & phase_roots, twist_phase, .false.)
      case default
         call bvk_separable_k_to_r(source_k_2d, source_r_2d, nmesh, &
            & input_to_grid, kernel%reps, phase_roots, twist_phase, .true.)
         call bvk_separable_k_to_r(result_adjoint_k_2d, result_adjoint_r_2d, &
            & nmesh, input_to_grid, kernel%reps, phase_roots, twist_phase, &
            & .false.)
      end select
   else
      source_r_2d = matmul(source_k_2d, phase_inverse)
      result_adjoint_r_2d = matmul(result_adjoint_k_2d, &
         & conjg(transpose(phase_forward)))
   end if

   origin = 0
   do icell = 1, size(kernel%reps, 2)
      if (all(kernel%reps(:, icell) == 0)) origin = icell
   end do

   if (use_bvk_kernel) then
      do icell = 1, size(kernel%reps, 2)
         if (mulliken_scale /= 0.0_wp) then
            call shell_parameter_gradient_complex(self%nsh, self%nao_sh, &
               & self%iao_sh, result_adjoint_r(:, :, icell), &
               & source_r(:, :, icell), mulliken_scale, &
               & mulliken_grad_r(:, :, icell))
            call shell_hadamard_add_complex(self%nsh, self%nao_sh, &
               & self%iao_sh, kernel%g_mulliken_r(:, :, icell), &
               & result_adjoint_r(:, :, icell), mulliken_scale, &
               & source_adjoint_r(:, :, icell))
         end if
         if (bocorr_scale /= 0.0_wp) then
            call atom_parameter_gradient_complex(mol%nat, mol%id, &
               & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, &
               & result_adjoint_r(:, :, icell), source_r(:, :, icell), &
               & bocorr_scale, bocorr_grad_r(:, :, icell))
            call atom_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
               & self%nao_sh, self%ish_at, self%iao_sh, &
               & kernel%g_bocorr_r(:, :, icell), &
               & result_adjoint_r(:, :, icell), bocorr_scale, &
               & source_adjoint_r(:, :, icell))
         end if
      end do
   end if

   if (origin > 0) then
      if (onsite_fx_scale /= 0.0_wp) &
         call onsite_fx_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
            & self%nao_sh, self%ish_at, self%iao_sh, g_onsfx, &
            & result_adjoint_r(:, :, origin), onsite_fx_scale, &
            & source_adjoint_r(:, :, origin), adjoint_src=adjoint_src)
      if (onsite_ri_scale /= 0.0_wp) &
         call onsite_ri_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
            & self%nao_sh, self%ish_at, self%iao_sh, g_onsri, &
            & result_adjoint_r(:, :, origin), onsite_ri_scale, &
            & source_adjoint_r(:, :, origin))
   end if

   if (present(transform)) then
      call transform%apply(source_adjoint_r_2d, source_adjoint_k_2d, exchange_k_to_r_adjoint)
   else if (is_compact_transform_backend(backend) .and. present(nmesh) &
      & .and. present(input_to_grid) .and. present(phase_roots) &
      & .and. present(twist_phase)) then
      select case (backend)
      case (bvk_transform_mixed_radix_fft)
         call bvk_mixed_radix_r_to_k(source_adjoint_r_2d, &
            & source_adjoint_k_2d, nmesh, input_to_grid, kernel%reps, &
            & phase_roots, twist_phase)
      case default
         call bvk_separable_r_to_k(source_adjoint_r_2d, source_adjoint_k_2d, &
            & nmesh, input_to_grid, kernel%reps, phase_roots, twist_phase)
      end select
   else
      source_adjoint_k_2d = matmul(source_adjoint_r_2d, phase_forward)
   end if
   if (.not.present(transform)) then
      do ik = 1, size(weights)
         source_adjoint_k(:, :, ik) = weights(ik)*source_adjoint_k(:, :, ik)
      end do
   end if

end subroutine reverse_bvk_Kkernel

subroutine get_KGrad_kmesh(self, mol, cache, plan, kernel, kfrac, weights, density, &
   & overlap, overlap_grad, mulliken_grad_r, bocorr_grad_r, transform_backend, transform)
   class(exchange_fock), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(exchange_cache), target, intent(in) :: cache
   type(cp2k_stream_plan_type), intent(in), target :: plan
   type(exchange_bvk_kernel), intent(in) :: kernel
   real(wp), intent(in) :: kfrac(:, :), weights(:)
   complex(wp), intent(in) :: density(:, :, :, :), overlap(:, :, :)
   complex(wp), intent(out) :: overlap_grad(:, :, :)
   real(wp), intent(out) :: mulliken_grad_r(:, :, :), bocorr_grad_r(:, :, :)
   integer, intent(in), optional :: transform_backend
   class(exchange_transform_type), intent(inout), optional :: transform

   integer :: iao, backend, icell, ik, spin
   logical :: use_cached_phases, use_compact
   real(wp) :: angle, spin_factor
   complex(wp) :: phase
   complex(wp), allocatable :: amat(:, :, :), &
      & bdiag(:), bdiag_p(:), bdiag_sp(:), bdiag_sps(:), &
      & btmp(:, :, :), dmat(:, :, :), ddensity(:, :, :), &
      & gdiag_p(:), gdiag_sp(:), gdiag_sps(:), tmat(:, :, :), &
      & umat(:, :, :), vmat(:, :, :), work(:, :, :), xmat(:, :, :)
   complex(wp), allocatable, target :: phase_forward_local(:, :), phase_inverse_local(:, :)
   integer, allocatable :: plan_input_to_grid(:)
   complex(wp), allocatable :: plan_phase_roots(:, :), plan_twist_phase(:)
   complex(wp), contiguous, pointer :: phase_forward(:, :), phase_inverse(:, :)

   allocate(amat(self%nao, self%nao, size(weights)), &
      & btmp(self%nao, self%nao, size(weights)), &
      & dmat(self%nao, self%nao, size(weights)), &
      & ddensity(self%nao, self%nao, size(weights)), &
      & tmat(self%nao, self%nao, size(weights)), &
      & umat(self%nao, self%nao, size(weights)), &
      & vmat(self%nao, self%nao, size(weights)), &
      & work(self%nao, self%nao, size(weights)), &
      & xmat(self%nao, self%nao, size(weights)), &
      & bdiag(self%nao), bdiag_p(self%nao), bdiag_sp(self%nao), &
      & bdiag_sps(self%nao), gdiag_p(self%nao), gdiag_sp(self%nao), &
      & gdiag_sps(self%nao), source=(0.0_wp, 0.0_wp))
   backend = bvk_transform_dense
   if (present(transform_backend)) backend = transform_backend
   use_compact = is_compact_transform_backend(backend) .and. &
      & plan%plan_matches(self, mol, kernel, kfrac, weights, &
      & require_phases=.false.) .and. allocated(plan%bvk_phase_roots) &
      & .and. allocated(plan%bvk_twist_phase)
   if (present(transform)) use_compact = .false.
   if (present(transform)) then
      allocate(plan_input_to_grid(1), source=1)
      allocate(plan_phase_roots(1, 3), source=(1.0_wp, 0.0_wp))
      allocate(plan_twist_phase(1), source=(1.0_wp, 0.0_wp))
      allocate(phase_forward_local(1, 1), phase_inverse_local(1, 1), &
         & source=(0.0_wp, 0.0_wp))
      phase_forward => phase_forward_local
      phase_inverse => phase_inverse_local
      use_cached_phases = .false.
   else if (use_compact) then
      allocate(plan_input_to_grid, source=plan%bvk_input_to_grid)
      allocate(plan_phase_roots, source=plan%bvk_phase_roots)
      allocate(plan_twist_phase, source=plan%bvk_twist_phase)
      allocate(phase_forward_local(1, 1), phase_inverse_local(1, 1), &
         & source=(0.0_wp, 0.0_wp))
      phase_forward => phase_forward_local
      phase_inverse => phase_inverse_local
      use_cached_phases = .false.
   else
      backend = bvk_transform_dense
      allocate(plan_input_to_grid(1), source=1)
      allocate(plan_phase_roots(1, 3), source=(1.0_wp, 0.0_wp))
      allocate(plan_twist_phase(1), source=(1.0_wp, 0.0_wp))
      use_cached_phases = plan%plan_matches(self, mol, kernel, kfrac, weights)
   end if
   if (.not.use_compact .and. use_cached_phases) then
      phase_forward => plan%bvk_phase_forward
      phase_inverse => plan%bvk_phase_inverse
   else if (.not.use_compact .and. .not.present(transform)) then
      allocate(phase_forward_local(size(kernel%reps, 2), size(weights)), &
         & phase_inverse_local(size(weights), size(kernel%reps, 2)))
      phase_forward => phase_forward_local
      phase_inverse => phase_inverse_local
      do ik = 1, size(weights)
         do icell = 1, size(kernel%reps, 2)
            angle = 2.0_wp*pi*dot_product(kfrac(:, ik), &
               & real(kernel%reps(:, icell), wp))
            phase = exp(cmplx(0.0_wp, angle, wp))
            phase_forward(icell, ik) = phase
            phase_inverse(ik, icell) = weights(ik)*conjg(phase)
         end do
      end do
   end if

   spin_factor = 0.5_wp
   if (size(density, 3) > 1) spin_factor = 1.0_wp
   overlap_grad = (0.0_wp, 0.0_wp)
   mulliken_grad_r = 0.0_wp
   bocorr_grad_r = 0.0_wp

   do spin = 1, size(density, 3)
      bdiag_p = (0.0_wp, 0.0_wp)
      bdiag_sp = (0.0_wp, 0.0_wp)
      bdiag_sps = (0.0_wp, 0.0_wp)
      gdiag_p = (0.0_wp, 0.0_wp)
      gdiag_sp = (0.0_wp, 0.0_wp)
      gdiag_sps = (0.0_wp, 0.0_wp)

      do ik = 1, size(weights)
         xmat(:, :, ik) = matmul(overlap(:, :, ik), &
            & density(:, :, spin, ik))
         dmat(:, :, ik) = 0.5_wp*matmul(xmat(:, :, ik), overlap(:, :, ik))
         do iao = 1, self%nao
            gdiag_p(iao) = gdiag_p(iao) &
               & + weights(ik)*density(iao, iao, spin, ik)
            gdiag_sp(iao) = gdiag_sp(iao) &
               & + weights(ik)*xmat(iao, iao, ik)
            gdiag_sps(iao) = gdiag_sps(iao) &
               & + 2.0_wp*weights(ik)*dmat(iao, iao, ik)
         end do
      end do
      call onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, &
         & gdiag_p, bdiag_p)
      call onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, &
         & gdiag_sp, bdiag_sp)
      call onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, &
         & gdiag_sps, bdiag_sps)
      gdiag_p = bdiag_p
      gdiag_sp = bdiag_sp
      gdiag_sps = bdiag_sps
      bdiag_sp = (0.0_wp, 0.0_wp)
      bdiag_sps = (0.0_wp, 0.0_wp)

      call apply_bvk_Kkernel(self, mol, kernel, phase_inverse, phase_forward, &
         & xmat, cache%g_onsfx, cache%g_onsri, 1.0_wp, -4.0_wp, &
         & 0.5_wp, -2.0_wp, .true., .true., tmat, backend, plan%bvk_nmesh, &
         & plan_input_to_grid, plan_phase_roots, plan_twist_phase, &
         & transform=transform)
      call apply_bvk_Kkernel(self, mol, kernel, phase_inverse, phase_forward, &
         & dmat, cache%g_onsfx, cache%g_onsri, 1.0_wp, 0.0_wp, &
         & 0.5_wp, -2.0_wp, .true., .false., umat, backend, plan%bvk_nmesh, &
         & plan_input_to_grid, plan_phase_roots, plan_twist_phase, &
         & transform=transform)
      call apply_bvk_Kkernel(self, mol, kernel, phase_inverse, phase_forward, &
         & density(:, :, spin, :), cache%g_onsfx, cache%g_onsri, &
         & 1.0_wp, 0.0_wp, 0.5_wp, -2.0_wp, .true., .false., vmat, backend, &
         & plan%bvk_nmesh, plan_input_to_grid, plan_phase_roots, &
         & plan_twist_phase, &
         & transform=transform)

      do ik = 1, size(weights)
         work(:, :, ik) = tmat(:, :, ik)
         do iao = 1, self%nao
            work(:, iao, ik) = work(:, iao, ik) &
               & + 0.25_wp*gdiag_p(iao)*overlap(:, iao, ik)
         end do
         work(:, :, ik) = work(:, :, ik) &
            & + 0.5_wp*matmul(overlap(:, :, ik), vmat(:, :, ik))

         ! Seed E = -(spin_factor/4) sum_k w_k Re<P_k,A_k>.
         amat(:, :, ik) = -0.25_wp*spin_factor*weights(ik) &
            & *density(:, :, spin, ik)
         umat(:, :, ik) = amat(:, :, ik)
         btmp(:, :, ik) = matmul(amat(:, :, ik), &
            & conjg(transpose(overlap(:, :, ik))))
         overlap_grad(:, :, ik) = overlap_grad(:, :, ik) &
            & + matmul(conjg(transpose(work(:, :, ik))), amat(:, :, ik))
         do iao = 1, self%nao
            overlap_grad(:, iao, ik) = overlap_grad(:, iao, ik) &
               & + 0.5_wp*conjg(gdiag_sp(iao))*amat(:, iao, ik)
            bdiag_sp(iao) = bdiag_sp(iao) &
               & + 0.5_wp*dot_product(overlap(:, iao, ik), &
               & amat(:, iao, ik))
            bdiag_sps(iao) = bdiag_sps(iao) &
               & + 0.25_wp*amat(iao, iao, ik)
         end do

         do iao = 1, self%nao
            overlap_grad(:, iao, ik) = overlap_grad(:, iao, ik) &
               & + 0.25_wp*conjg(gdiag_p(iao))*btmp(:, iao, ik)
         end do
         overlap_grad(:, :, ik) = overlap_grad(:, :, ik) &
            & + 0.5_wp*matmul(btmp(:, :, ik), &
            & conjg(transpose(vmat(:, :, ik))))
         vmat(:, :, ik) = 0.5_wp*matmul( &
            & conjg(transpose(overlap(:, :, ik))), btmp(:, :, ik))
      end do

      call reverse_bvk_Kkernel(self, mol, kernel, weights, phase_inverse, &
         & phase_forward, dmat, umat, cache%g_onsfx, cache%g_onsri, &
         & 1.0_wp, 0.0_wp, 0.5_wp, -2.0_wp, .true., .false., tmat, &
         & mulliken_grad_r, bocorr_grad_r, backend, plan%bvk_nmesh, &
         & plan_input_to_grid, plan_phase_roots, plan_twist_phase, &
         & transform=transform)
      call reverse_bvk_Kkernel(self, mol, kernel, weights, phase_inverse, &
         & phase_forward, xmat, btmp, cache%g_onsfx, cache%g_onsri, &
         & 1.0_wp, -4.0_wp, 0.5_wp, -2.0_wp, .true., .true., umat, &
         & mulliken_grad_r, bocorr_grad_r, backend, plan%bvk_nmesh, &
         & plan_input_to_grid, plan_phase_roots, plan_twist_phase, &
         & transform=transform)
      call reverse_bvk_Kkernel(self, mol, kernel, weights, phase_inverse, &
         & phase_forward, density(:, :, spin, :), vmat, cache%g_onsfx, &
         & cache%g_onsri, 1.0_wp, 0.0_wp, 0.5_wp, -2.0_wp, .true., &
         & .false., ddensity, mulliken_grad_r, bocorr_grad_r, backend, &
         & plan%bvk_nmesh, plan_input_to_grid, plan_phase_roots, &
         & plan_twist_phase, &
         & transform=transform)

      ! Reverse the three BZ-averaged onsite diagonal maps.  diag(P) has no
      ! overlap response; its kernel response belongs to the shell potential.
      call onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, &
         & bdiag_sp, bdiag)
      do ik = 1, size(weights)
         do iao = 1, self%nao
            umat(iao, iao, ik) = umat(iao, iao, ik) &
               & + weights(ik)*bdiag(iao)
         end do
      end do
      call onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, &
         & bdiag_sps, bdiag)
      do ik = 1, size(weights)
         do iao = 1, self%nao
            umat(iao, :, ik) = umat(iao, :, ik) &
               & + weights(ik)*bdiag(iao)*overlap(iao, :, ik)
            overlap_grad(:, iao, ik) = overlap_grad(:, iao, ik) &
               & + weights(ik)*bdiag(iao)*conjg(xmat(iao, :, ik))
         end do

         ! D_k = X_k S_k/2 and X_k = S_k P_k.
         umat(:, :, ik) = umat(:, :, ik) &
            & + 0.5_wp*matmul(tmat(:, :, ik), &
            & conjg(transpose(overlap(:, :, ik))))
         overlap_grad(:, :, ik) = overlap_grad(:, :, ik) &
            & + 0.5_wp*matmul(conjg(transpose(xmat(:, :, ik))), &
            & tmat(:, :, ik))
         overlap_grad(:, :, ik) = overlap_grad(:, :, ik) &
            & + matmul(umat(:, :, ik), &
            & conjg(transpose(density(:, :, spin, ik))))
      end do
   end do

   do ik = 1, size(weights)
      overlap_grad(:, :, ik) = 0.5_wp*(overlap_grad(:, :, ik) &
         & + conjg(transpose(overlap_grad(:, :, ik))))/weights(ik)
   end do

end subroutine get_KGrad_kmesh

subroutine shell_hadamard_add_complex(nsh, nao_sh, iao_sh, g_sh, src, &
   & alpha, dst, adjoint_src)
   integer, intent(in) :: nsh, nao_sh(:), iao_sh(:)
   real(wp), intent(in) :: g_sh(:, :), alpha
   complex(wp), intent(in) :: src(:, :)
   complex(wp), intent(inout) :: dst(:, :)
   logical, intent(in), optional :: adjoint_src

   integer :: iao, ii, ish, jao, jj, jsh, ni, nj
   real(wp) :: scale
   logical :: adjoint

   adjoint = .false.
   if (present(adjoint_src)) adjoint = adjoint_src
   do ish = 1, nsh
      ii = iao_sh(ish)
      ni = nao_sh(ish)
      do jsh = 1, nsh
         jj = iao_sh(jsh)
         nj = nao_sh(jsh)
         scale = alpha*g_sh(jsh, ish)
         if (abs(scale) < epsilon(1.0_wp)) cycle
         if (.not. adjoint) then
            dst(jj+1:jj+nj, ii+1:ii+ni) = dst(jj+1:jj+nj, ii+1:ii+ni) &
               & + scale*src(jj+1:jj+nj, ii+1:ii+ni)
         else
            do iao = 1, ni
               do jao = 1, nj
                  dst(jj+jao, ii+iao) = dst(jj+jao, ii+iao) &
                     & + scale*conjg(src(ii+iao, jj+jao))
               end do
            end do
         end if
      end do
   end do
end subroutine shell_hadamard_add_complex

subroutine shell_parameter_gradient_complex(nsh, nao_sh, iao_sh, adjoint, &
   & source, alpha, gradient)
   integer, intent(in) :: nsh, nao_sh(:), iao_sh(:)
   complex(wp), intent(in) :: adjoint(:, :), source(:, :)
   real(wp), intent(in) :: alpha
   real(wp), intent(inout) :: gradient(:, :)

   integer :: iao, ii, ish, jao, jj, jsh, ni, nj

   do ish = 1, nsh
      ii = iao_sh(ish)
      ni = nao_sh(ish)
      do jsh = 1, nsh
         jj = iao_sh(jsh)
         nj = nao_sh(jsh)
         do iao = 1, ni
            do jao = 1, nj
               gradient(jj+jao, ii+iao) = gradient(jj+jao, ii+iao) &
                  & + alpha*real(conjg(adjoint(jj+jao, ii+iao)) &
                  & * source(jj+jao, ii+iao), wp)
            end do
         end do
      end do
   end do
end subroutine shell_parameter_gradient_complex

subroutine atom_hadamard_add_complex(nat, id, nsh_id, nao_sh, ish_at, iao_sh, &
   & g_at, src, alpha, dst, adjoint_src)
   integer, intent(in) :: nat, id(:), nsh_id(:), nao_sh(:), ish_at(:), iao_sh(:)
   real(wp), intent(in) :: g_at(:, :), alpha
   complex(wp), intent(in) :: src(:, :)
   complex(wp), intent(inout) :: dst(:, :)
   logical, intent(in), optional :: adjoint_src

   integer :: iao, iat, ii, is, ish, izp, jao, jat, jj, js, jsh, jzp, ni, nj
   real(wp) :: scale
   logical :: adjoint

   adjoint = .false.
   if (present(adjoint_src)) adjoint = adjoint_src
   do iat = 1, nat
      izp = id(iat)
      is = ish_at(iat)
      do jat = 1, nat
         jzp = id(jat)
         js = ish_at(jat)
         scale = alpha*g_at(jat, iat)
         if (abs(scale) < epsilon(1.0_wp)) cycle
         do ish = 1, nsh_id(izp)
            ii = iao_sh(is+ish)
            ni = nao_sh(is+ish)
            do jsh = 1, nsh_id(jzp)
               jj = iao_sh(js+jsh)
               nj = nao_sh(js+jsh)
               if (.not. adjoint) then
                  dst(jj+1:jj+nj, ii+1:ii+ni) = &
                     & dst(jj+1:jj+nj, ii+1:ii+ni) &
                     & + scale*src(jj+1:jj+nj, ii+1:ii+ni)
               else
                  do iao = 1, ni
                     do jao = 1, nj
                        dst(jj+jao, ii+iao) = dst(jj+jao, ii+iao) &
                           & + scale*conjg(src(ii+iao, jj+jao))
                     end do
                  end do
               end if
            end do
         end do
      end do
   end do
end subroutine atom_hadamard_add_complex

subroutine atom_parameter_gradient_complex(nat, id, nsh_id, nao_sh, ish_at, &
   & iao_sh, adjoint, source, alpha, gradient)
   integer, intent(in) :: nat, id(:), nsh_id(:), nao_sh(:), ish_at(:), iao_sh(:)
   complex(wp), intent(in) :: adjoint(:, :), source(:, :)
   real(wp), intent(in) :: alpha
   real(wp), intent(inout) :: gradient(:, :)

   integer :: iao, iat, ii, is, ish, izp, jao, jat, jj, js, jsh, jzp, ni, nj

   do iat = 1, nat
      izp = id(iat)
      is = ish_at(iat)
      do jat = 1, nat
         jzp = id(jat)
         js = ish_at(jat)
         do ish = 1, nsh_id(izp)
            ii = iao_sh(is+ish)
            ni = nao_sh(is+ish)
            do jsh = 1, nsh_id(jzp)
               jj = iao_sh(js+jsh)
               nj = nao_sh(js+jsh)
               do iao = 1, ni
                  do jao = 1, nj
                     gradient(jj+jao, ii+iao) = gradient(jj+jao, ii+iao) &
                        & + alpha*real(conjg(adjoint(jj+jao, ii+iao)) &
                        & * source(jj+jao, ii+iao), wp)
                  end do
               end do
            end do
         end do
      end do
   end do
end subroutine atom_parameter_gradient_complex

subroutine onsite_fx_hadamard_add_complex(nat, id, nsh_id, nao_sh, ish_at, &
   & iao_sh, g_onsfx, src, alpha, dst, adjoint_src)
   integer, intent(in) :: nat, id(:), nsh_id(:), nao_sh(:), ish_at(:), iao_sh(:)
   real(wp), intent(in) :: g_onsfx(:, :, :), alpha
   complex(wp), intent(in) :: src(:, :)
   complex(wp), intent(inout) :: dst(:, :)
   logical, intent(in), optional :: adjoint_src

   integer :: iao, iat, ii, is, ish, izp, jao, jj, jsh, ni, nj
   real(wp) :: scale
   logical :: adjoint

   adjoint = .false.
   if (present(adjoint_src)) adjoint = adjoint_src
   do iat = 1, nat
      izp = id(iat)
      is = ish_at(iat)
      do ish = 1, nsh_id(izp)
         ii = iao_sh(is+ish)
         ni = nao_sh(is+ish)
         do jsh = 1, nsh_id(izp)
            jj = iao_sh(is+jsh)
            nj = nao_sh(is+jsh)
            scale = alpha*g_onsfx(jsh, ish, iat)
            if (abs(scale) < epsilon(1.0_wp)) cycle
            if (.not. adjoint) then
               dst(jj+1:jj+nj, ii+1:ii+ni) = &
                  & dst(jj+1:jj+nj, ii+1:ii+ni) &
                  & + scale*src(jj+1:jj+nj, ii+1:ii+ni)
            else
               do iao = 1, ni
                  do jao = 1, nj
                     dst(jj+jao, ii+iao) = dst(jj+jao, ii+iao) &
                        & + scale*conjg(src(ii+iao, jj+jao))
                  end do
               end do
            end if
         end do
      end do
   end do
end subroutine onsite_fx_hadamard_add_complex

subroutine onsite_fx_symv_complex(nat, id, nsh_id, nao_sh, ish_at, iao_sh, &
   & g_onsfx, xvec, yvec)
   integer, intent(in) :: nat, id(:), nsh_id(:), nao_sh(:), ish_at(:), iao_sh(:)
   real(wp), intent(in) :: g_onsfx(:, :, :)
   complex(wp), intent(in) :: xvec(:)
   complex(wp), intent(out) :: yvec(:)

   integer :: iat, ii, is, ish, izp, jj, jsh, ni, nj
   complex(wp) :: shsum

   yvec = (0.0_wp, 0.0_wp)
   do iat = 1, nat
      izp = id(iat)
      is = ish_at(iat)
      do ish = 1, nsh_id(izp)
         ii = iao_sh(is+ish)
         ni = nao_sh(is+ish)
         shsum = sum(xvec(ii+1:ii+ni))
         if (abs(shsum) < epsilon(1.0_wp)) cycle
         do jsh = 1, nsh_id(izp)
            jj = iao_sh(is+jsh)
            nj = nao_sh(is+jsh)
            yvec(jj+1:jj+nj) = yvec(jj+1:jj+nj) &
               & + g_onsfx(jsh, ish, iat)*shsum
         end do
      end do
   end do
end subroutine onsite_fx_symv_complex

subroutine onsite_ri_hadamard_add_complex(nat, id, nsh_id, nao_sh, ish_at, &
   & iao_sh, g_onsri, src, alpha, dst)
   integer, intent(in) :: nat, id(:), nsh_id(:), nao_sh(:), ish_at(:), iao_sh(:)
   real(wp), intent(in) :: g_onsri(:, :), alpha
   complex(wp), intent(in) :: src(:, :)
   complex(wp), intent(inout) :: dst(:, :)

   integer :: iat, ii, is, ish, izp, ni
   real(wp) :: scale

   do iat = 1, nat
      izp = id(iat)
      is = ish_at(iat)
      do ish = 1, nsh_id(izp)
         ii = iao_sh(is+ish)
         ni = nao_sh(is+ish)
         scale = alpha*g_onsri(ish, iat)
         if (abs(scale) < epsilon(1.0_wp)) cycle
         dst(ii+1:ii+ni, ii+1:ii+ni) = &
            & dst(ii+1:ii+ni, ii+1:ii+ni) &
            & + scale*src(ii+1:ii+ni, ii+1:ii+ni)
      end do
   end do
end subroutine onsite_ri_hadamard_add_complex

subroutine onsite_fx_parameter_gradient_complex(nat, id, nsh_id, nao_sh, &
   & ish_at, iao_sh, adjoint, src, alpha, gradient, adjoint_src)
   integer, intent(in) :: nat, id(:), nsh_id(:), nao_sh(:), ish_at(:), iao_sh(:)
   complex(wp), intent(in) :: adjoint(:, :), src(:, :)
   real(wp), intent(in) :: alpha
   real(wp), intent(inout) :: gradient(:, :, :)
   logical, intent(in), optional :: adjoint_src

   integer :: iao, iat, ii, is, ish, izp, jao, jj, jsh, ni, nj
   real(wp) :: value
   logical :: take_adjoint

   take_adjoint = .false.
   if (present(adjoint_src)) take_adjoint = adjoint_src
   do iat = 1, nat
      izp = id(iat)
      is = ish_at(iat)
      do ish = 1, nsh_id(izp)
         ii = iao_sh(is+ish)
         ni = nao_sh(is+ish)
         do jsh = 1, nsh_id(izp)
            jj = iao_sh(is+jsh)
            nj = nao_sh(is+jsh)
            value = 0.0_wp
            if (.not.take_adjoint) then
               value = real(sum(conjg(adjoint(jj+1:jj+nj, ii+1:ii+ni)) &
                  & *src(jj+1:jj+nj, ii+1:ii+ni)), wp)
            else
               do iao = 1, ni
                  do jao = 1, nj
                     value = value + real(conjg(adjoint(jj+jao, ii+iao)) &
                        & *conjg(src(ii+iao, jj+jao)), wp)
                  end do
               end do
            end if
            gradient(jsh, ish, iat) = gradient(jsh, ish, iat) + alpha*value
         end do
      end do
   end do
end subroutine onsite_fx_parameter_gradient_complex

subroutine onsite_ri_parameter_gradient_complex(nat, id, nsh_id, nao_sh, &
   & ish_at, iao_sh, adjoint, src, alpha, gradient)
   integer, intent(in) :: nat, id(:), nsh_id(:), nao_sh(:), ish_at(:), iao_sh(:)
   complex(wp), intent(in) :: adjoint(:, :), src(:, :)
   real(wp), intent(in) :: alpha
   real(wp), intent(inout) :: gradient(:, :)

   integer :: iat, ii, is, ish, izp, ni

   do iat = 1, nat
      izp = id(iat)
      is = ish_at(iat)
      do ish = 1, nsh_id(izp)
         ii = iao_sh(is+ish)
         ni = nao_sh(is+ish)
         gradient(ish, iat) = gradient(ish, iat) + alpha*real(sum( &
            & conjg(adjoint(ii+1:ii+ni, ii+1:ii+ni)) &
            & *src(ii+1:ii+ni, ii+1:ii+ni)), wp)
      end do
   end do
end subroutine onsite_ri_parameter_gradient_complex

subroutine onsite_fx_symv_parameter_gradient_complex(nat, id, nsh_id, &
   & nao_sh, ish_at, iao_sh, adjoint, src, gradient)
   integer, intent(in) :: nat, id(:), nsh_id(:), nao_sh(:), ish_at(:), iao_sh(:)
   complex(wp), intent(in) :: adjoint(:), src(:)
   real(wp), intent(inout) :: gradient(:, :, :)

   integer :: iat, ii, is, ish, izp, jj, jsh, ni, nj
   complex(wp) :: adjoint_sum, source_sum

   do iat = 1, nat
      izp = id(iat)
      is = ish_at(iat)
      do ish = 1, nsh_id(izp)
         ii = iao_sh(is+ish)
         ni = nao_sh(is+ish)
         source_sum = sum(src(ii+1:ii+ni))
         do jsh = 1, nsh_id(izp)
            jj = iao_sh(is+jsh)
            nj = nao_sh(is+jsh)
            adjoint_sum = sum(adjoint(jj+1:jj+nj))
            gradient(jsh, ish, iat) = gradient(jsh, ish, iat) &
               & + real(conjg(adjoint_sum)*source_sum, wp)
         end do
      end do
   end do
end subroutine onsite_fx_symv_parameter_gradient_complex
subroutine get_KFock_kpoint(self, mol, cache, density, overlap, fock, vsh, energy)
   !> Instance of the Fock exchange container
   class(exchange_fock), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable exchange kernels and charge derivatives
   type(exchange_cache), intent(in) :: cache
   !> Complex Hermitian density matrix for all spin channels
   complex(wp), intent(in) :: density(:, :, :)
   !> Complex Hermitian overlap matrix
   complex(wp), intent(in) :: overlap(:, :)
   !> Complex Hermitian exchange Fock matrix
   complex(wp), intent(out) :: fock(:, :, :)
   !> Charge-dependent shell potential
   real(wp), intent(out) :: vsh(:)
   !> Optional exchange energy for this k-point block
   real(wp), intent(out), optional :: energy

   integer :: iat, izp, is, ish, iq, spin
   real(wp), allocatable :: dg_onsfx(:, :, :), dg_onsri(:, :)
   real(wp), allocatable :: zero_at(:, :), zero_sh(:, :)
   complex(wp), allocatable :: dfock(:, :, :)
   real(wp) :: tmp_energy

   call build_KFock_complex(self, mol, density, overlap, cache%g_mulliken, &
      & cache%g_bocorr, cache%g_onsfx, cache%g_onsri, fock)

   tmp_energy = 0.0_wp
   do spin = 1, size(density, 3)
      tmp_energy = tmp_energy + 0.5_wp*real(sum(conjg(density(:, :, spin)) &
         & * fock(:, :, spin)), wp)
   end do
   if (present(energy)) energy = tmp_energy

   allocate(zero_sh(self%nsh, self%nsh), zero_at(mol%nat, mol%nat), &
      & dg_onsfx(self%maxsh, self%maxsh, mol%nat), &
      & dg_onsri(self%maxsh, mol%nat), source=0.0_wp)
   allocate(dfock(self%nao, self%nao, size(density, 3)), &
      & source=(0.0_wp, 0.0_wp))
   vsh = 0.0_wp

   ! The onsite kernels are the only exchange kernels that depend explicitly
   ! on shell charges.  Apply their sparse derivative through the same linear
   ! Fock map, then contract with P.  This is intentionally an exact reference
   ! implementation; the contractions can be fused after validation in CP2K.
   do iat = 1, mol%nat
      izp = mol%id(iat)
      is = self%ish_at(iat)
      do ish = 1, self%nsh_id(izp)
         iq = is + ish
         dg_onsfx = 0.0_wp
         dg_onsri = 0.0_wp
         dg_onsfx(:, :, iat) = cache%dgdq_onsfx(:, :, iq)
         dg_onsri(:, iat) = cache%dgdq_onsri(:, iq)
         call build_KFock_complex(self, mol, density, overlap, zero_sh, &
            & zero_at, dg_onsfx, dg_onsri, dfock)
         do spin = 1, size(density, 3)
            vsh(iq) = vsh(iq) + 0.5_wp*real(sum( &
               & conjg(density(:, :, spin))*dfock(:, :, spin)), wp)
         end do
      end do
   end do

end subroutine get_KFock_kpoint

subroutine get_KGrad_kpoint(self, mol, cache, density, overlap, &
   & mulliken_grad, bocorr_grad, overlap_grad)
   class(exchange_fock), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(exchange_cache), intent(in) :: cache
   complex(wp), intent(in) :: density(:, :, :), overlap(:, :)
   real(wp), intent(out) :: mulliken_grad(:, :), bocorr_grad(:, :)
   complex(wp), intent(out) :: overlap_grad(:, :)

   integer :: iao, spin
   real(wp) :: spin_factor
   real(wp), allocatable :: mulliken_raw(:, :), bocorr_raw(:, :)
   complex(wp), allocatable :: amat(:, :), bamat(:, :), bdiag(:), &
      & bdiag_p(:), bdiag_sp(:), bdiag_sps(:), bmat(:, :), btmp(:, :), &
      & dmat(:, :), gdiag_p(:), gdiag_sp(:), gdiag_sps(:), tmat(:, :), &
      & t2mat(:, :), umat(:, :), vmat(:, :), xmat(:, :)

   allocate(amat(self%nao, self%nao), bamat(self%nao, self%nao), &
      & bmat(self%nao, self%nao), btmp(self%nao, self%nao), &
      & dmat(self%nao, self%nao), tmat(self%nao, self%nao), &
      & t2mat(self%nao, self%nao), umat(self%nao, self%nao), &
      & vmat(self%nao, self%nao), xmat(self%nao, self%nao), &
      & bdiag(self%nao), bdiag_p(self%nao), bdiag_sp(self%nao), &
      & bdiag_sps(self%nao), gdiag_p(self%nao), gdiag_sp(self%nao), &
      & gdiag_sps(self%nao), source=(0.0_wp, 0.0_wp))
   allocate(mulliken_raw(self%nao, self%nao), &
      & bocorr_raw(self%nao, self%nao), source=0.0_wp)

   spin_factor = 0.5_wp
   if (size(density, 3) > 1) spin_factor = 1.0_wp
   overlap_grad = (0.0_wp, 0.0_wp)

   do spin = 1, size(density, 3)
      ! Forward intermediates of build_KFock_complex, kept explicitly for the
      ! reverse sweep below.
      xmat = matmul(overlap, density(:, :, spin))
      dmat = 0.5_wp*matmul(xmat, overlap)
      call onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, &
         & [(density(iao, iao, spin), iao=1, self%nao)], gdiag_p)
      call onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, &
         & [(xmat(iao, iao), iao=1, self%nao)], gdiag_sp)
      call onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, &
         & [(2.0_wp*dmat(iao, iao), iao=1, self%nao)], gdiag_sps)

      tmat = (0.0_wp, 0.0_wp)
      call shell_hadamard_add_complex(self%nsh, self%nao_sh, self%iao_sh, &
         & cache%g_mulliken, xmat, 1.0_wp, tmat)
      call atom_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, cache%g_bocorr, xmat, &
         & -4.0_wp, tmat)
      call onsite_fx_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, xmat, &
         & 0.5_wp, tmat, adjoint_src=.true.)
      call onsite_ri_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsri, xmat, &
         & -2.0_wp, tmat)

      umat = (0.0_wp, 0.0_wp)
      call shell_hadamard_add_complex(self%nsh, self%nao_sh, self%iao_sh, &
         & cache%g_mulliken, dmat, 1.0_wp, umat)
      call onsite_fx_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, dmat, &
         & 0.5_wp, umat)
      call onsite_ri_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsri, dmat, &
         & -2.0_wp, umat)

      vmat = (0.0_wp, 0.0_wp)
      call shell_hadamard_add_complex(self%nsh, self%nao_sh, self%iao_sh, &
         & cache%g_mulliken, density(:, :, spin), 1.0_wp, vmat)
      call onsite_fx_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, &
         & density(:, :, spin), 0.5_wp, vmat)
      call onsite_ri_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsri, &
         & density(:, :, spin), -2.0_wp, vmat)

      t2mat = tmat + 0.5_wp*matmul(overlap, vmat)
      do iao = 1, self%nao
         t2mat(:, iao) = t2mat(:, iao) &
            & + 0.25_wp*gdiag_p(iao)*overlap(:, iao)
      end do
      amat = umat + matmul(t2mat, overlap)
      do iao = 1, self%nao
         amat(:, iao) = amat(:, iao) &
            & + 0.5_wp*gdiag_sp(iao)*overlap(:, iao)
         amat(iao, iao) = amat(iao, iao) + 0.25_wp*gdiag_sps(iao)
      end do

      ! Reverse sweep, seeded by E=-(spin_factor/4)*Re<P,A>.
      bamat = -0.25_wp*spin_factor*density(:, :, spin)
      bmat = bamat  ! adjoint of U
      btmp = matmul(bamat, conjg(transpose(overlap)))  ! adjoint of T2
      overlap_grad = overlap_grad + matmul(conjg(transpose(t2mat)), bamat)
      bdiag_sp = (0.0_wp, 0.0_wp)
      bdiag_sps = (0.0_wp, 0.0_wp)
      do iao = 1, self%nao
         overlap_grad(:, iao) = overlap_grad(:, iao) &
            & + 0.5_wp*conjg(gdiag_sp(iao))*bamat(:, iao)
         bdiag_sp(iao) = 0.5_wp*dot_product(overlap(:, iao), bamat(:, iao))
         bdiag_sps(iao) = 0.25_wp*bamat(iao, iao)
      end do

      ! T2 = T + 1/4 S diag(gdiagP) + 1/2 S V.
      tmat = btmp  ! adjoint of T
      bdiag_p = (0.0_wp, 0.0_wp)
      do iao = 1, self%nao
         overlap_grad(:, iao) = overlap_grad(:, iao) &
            & + 0.25_wp*conjg(gdiag_p(iao))*btmp(:, iao)
         bdiag_p(iao) = 0.25_wp*dot_product(overlap(:, iao), btmp(:, iao))
      end do
      overlap_grad = overlap_grad &
         & + 0.5_wp*matmul(btmp, conjg(transpose(vmat)))
      vmat = 0.5_wp*matmul(conjg(transpose(overlap)), btmp) ! adjoint of V

      ! U response: U = Gm(D) + 1/2 Gfx(D) - 2 Gri(D).
      dmat = (0.0_wp, 0.0_wp)
      call shell_hadamard_add_complex(self%nsh, self%nao_sh, self%iao_sh, &
         & cache%g_mulliken, bmat, 1.0_wp, dmat)
      call onsite_fx_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, bmat, &
         & 0.5_wp, dmat)
      call onsite_ri_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsri, bmat, &
         & -2.0_wp, dmat)
      call shell_parameter_gradient_complex(self%nsh, self%nao_sh, &
         & self%iao_sh, bmat, 0.5_wp*matmul(xmat, overlap), 1.0_wp, &
         & mulliken_raw)

      ! T response, including the adjoint of Gfx(X^H).
      bamat = (0.0_wp, 0.0_wp) ! adjoint of X
      call shell_hadamard_add_complex(self%nsh, self%nao_sh, self%iao_sh, &
         & cache%g_mulliken, tmat, 1.0_wp, bamat)
      call atom_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, cache%g_bocorr, tmat, &
         & -4.0_wp, bamat)
      call onsite_fx_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, tmat, &
         & 0.5_wp, bamat, adjoint_src=.true.)
      call onsite_ri_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsri, tmat, &
         & -2.0_wp, bamat)
      call shell_parameter_gradient_complex(self%nsh, self%nao_sh, &
         & self%iao_sh, tmat, xmat, 1.0_wp, mulliken_raw)
      call atom_parameter_gradient_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, tmat, xmat, -4.0_wp, &
         & bocorr_raw)

      ! V response contributes only to the Mulliken geometry kernel.
      call shell_parameter_gradient_complex(self%nsh, self%nao_sh, &
         & self%iao_sh, vmat, density(:, :, spin), 1.0_wp, mulliken_raw)

      ! Reverse the three onsite diagonal maps.
      call onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, &
         & bdiag_sp, bdiag)
      do iao = 1, self%nao
         bamat(iao, iao) = bamat(iao, iao) + bdiag(iao)
      end do
      call onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, &
         & bdiag_sps, bdiag)
      do iao = 1, self%nao
         bamat(iao, :) = bamat(iao, :) + bdiag(iao)*overlap(iao, :)
         overlap_grad(:, iao) = overlap_grad(:, iao) &
            & + bdiag(iao)*conjg(xmat(iao, :))
      end do

      ! D = 1/2 X S and X = S P.
      bamat = bamat + 0.5_wp*matmul(dmat, conjg(transpose(overlap)))
      overlap_grad = overlap_grad &
         & + 0.5_wp*matmul(conjg(transpose(xmat)), dmat)
      overlap_grad = overlap_grad &
         & + matmul(bamat, conjg(transpose(density(:, :, spin))))
   end do

   ! Only Hermitian overlap variations are physical at a full k point.
   overlap_grad = 0.5_wp*(overlap_grad+conjg(transpose(overlap_grad)))
   call symmetrize_shell_parameter_gradient(self%nsh, self%nao_sh, &
      & self%iao_sh, mulliken_raw, mulliken_grad)
   call symmetrize_atom_parameter_gradient(mol%nat, mol%id, self%nsh_id, &
      & self%nao_sh, self%ish_at, self%iao_sh, bocorr_raw, bocorr_grad)

end subroutine get_KGrad_kpoint

subroutine build_KFock_complex(self, mol, density, overlap, g_mulliken, &
   & g_bocorr, g_onsfx, g_onsri, fock)
   class(exchange_fock), intent(in) :: self
   type(structure_type), intent(in) :: mol
   complex(wp), intent(in) :: density(:, :, :), overlap(:, :)
   real(wp), intent(in) :: g_mulliken(:, :), g_bocorr(:, :)
   real(wp), intent(in) :: g_onsfx(:, :, :), g_onsri(:, :)
   complex(wp), intent(out) :: fock(:, :, :)

   integer :: iao, ii, jj, spin
   real(wp) :: spin_factor
   complex(wp) :: tmp
   complex(wp), allocatable :: diagP(:), diagSP(:), gdiagP(:), &
      & gdiagSP(:), gdiagSPS(:), tmpA(:, :), tmpB(:, :), tmpSPS(:)

   allocate(tmpA(self%nao, self%nao), tmpB(self%nao, self%nao), &
      & diagP(self%nao), diagSP(self%nao), gdiagP(self%nao), &
      & gdiagSP(self%nao), tmpSPS(self%nao), gdiagSPS(self%nao), &
      & source=(0.0_wp, 0.0_wp))

   spin_factor = 0.5_wp
   if (size(density, 3) > 1) spin_factor = 1.0_wp
   fock = (0.0_wp, 0.0_wp)

   do spin = 1, size(density, 3)
      ! A = S P
      tmpA = matmul(overlap, density(:, :, spin))

      do iao = 1, self%nao
         diagP(iao) = density(iao, iao, spin)
         diagSP(iao) = tmpA(iao, iao)
         tmpSPS(iao) = sum(tmpA(iao, :)*conjg(overlap(iao, :)))
      end do
      call onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, g_onsfx, diagP, gdiagP)

      call onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, g_onsfx, diagSP, gdiagSP)
      call onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, g_onsfx, tmpSPS, gdiagSPS)

      ! Initial 1/2 S P S intermediate.
      fock(:, :, spin) = 0.5_wp*matmul(tmpA, overlap)

      ! X = g*(S P), including the adjoint P S onsite contribution.
      tmpB = (0.0_wp, 0.0_wp)
      call shell_hadamard_add_complex(self%nsh, self%nao_sh, self%iao_sh, &
         & g_mulliken, tmpA, 1.0_wp, tmpB)
      call atom_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, g_bocorr, tmpA, &
         & -4.0_wp, tmpB)
      call onsite_fx_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, g_onsfx, tmpA, &
         & 0.5_wp, tmpB, adjoint_src=.true.)
      call onsite_ri_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, g_onsri, tmpA, &
         & -2.0_wp, tmpB)
      tmpA = tmpB

      ! g*(1/2 S P S)
      tmpB = (0.0_wp, 0.0_wp)
      call onsite_fx_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, g_onsfx, &
         & fock(:, :, spin), 0.5_wp, tmpB)
      call onsite_ri_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, g_onsri, &
         & fock(:, :, spin), -2.0_wp, tmpB)
      call shell_hadamard_add_complex(self%nsh, self%nao_sh, self%iao_sh, &
         & g_mulliken, fock(:, :, spin), 1.0_wp, tmpB)
      fock(:, :, spin) = tmpB

      ! g*P
      tmpB = (0.0_wp, 0.0_wp)
      call shell_hadamard_add_complex(self%nsh, self%nao_sh, self%iao_sh, &
         & g_mulliken, density(:, :, spin), 1.0_wp, tmpB)
      call onsite_fx_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, g_onsfx, &
         & density(:, :, spin), 0.5_wp, tmpB)
      call onsite_ri_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, g_onsri, &
         & density(:, :, spin), -2.0_wp, tmpB)

      do iao = 1, self%nao
         tmpA(:, iao) = tmpA(:, iao) &
            & + 0.25_wp*gdiagP(iao)*overlap(:, iao)
      end do
      tmpA = tmpA + 0.5_wp*matmul(overlap, tmpB)
      fock(:, :, spin) = fock(:, :, spin) + matmul(tmpA, overlap)

      do iao = 1, self%nao
         fock(:, iao, spin) = fock(:, iao, spin) &
            & + 0.5_wp*gdiagSP(iao)*overlap(:, iao)
         fock(iao, iao, spin) = fock(iao, iao, spin) &
            & + 0.25_wp*gdiagSPS(iao)
      end do

      ! Hermitian projection and exchange sign/spin scaling.
      do ii = 1, self%nao
         fock(ii, ii, spin) = cmplx(-0.5_wp*spin_factor &
            & * real(fock(ii, ii, spin), wp), 0.0_wp, wp)
         do jj = 1, ii-1
            tmp = -0.25_wp*spin_factor*(fock(jj, ii, spin) &
               & + conjg(fock(ii, jj, spin)))
            fock(jj, ii, spin) = tmp
            fock(ii, jj, spin) = conjg(tmp)
         end do
      end do
   end do

end subroutine build_KFock_complex

subroutine symmetrize_shell_parameter_gradient(nsh, nao_sh, iao_sh, raw, &
   & gradient)
   integer, intent(in) :: nsh, nao_sh(:), iao_sh(:)
   real(wp), intent(in) :: raw(:, :)
   real(wp), intent(out) :: gradient(:, :)

   integer :: iao, ii, ish, jao, jj, jsh, ni, nj
   real(wp) :: tmp

   gradient = raw
   do ish = 1, nsh
      ii = iao_sh(ish)
      ni = nao_sh(ish)
      do jsh = 1, ish-1
         jj = iao_sh(jsh)
         nj = nao_sh(jsh)
         do iao = 1, ni
            do jao = 1, nj
               tmp = raw(jj+jao, ii+iao) + raw(ii+iao, jj+jao)
               gradient(jj+jao, ii+iao) = tmp
               gradient(ii+iao, jj+jao) = tmp
            end do
         end do
      end do
   end do
end subroutine symmetrize_shell_parameter_gradient

subroutine symmetrize_atom_parameter_gradient(nat, id, nsh_id, nao_sh, ish_at, &
   & iao_sh, raw, gradient)
   integer, intent(in) :: nat, id(:), nsh_id(:), nao_sh(:), ish_at(:), iao_sh(:)
   real(wp), intent(in) :: raw(:, :)
   real(wp), intent(out) :: gradient(:, :)

   integer :: iao, iat, ii, is, ish, izp, jao, jat, jj, js, jsh, jzp, ni, nj
   real(wp) :: tmp

   gradient = raw
   do iat = 1, nat
      izp = id(iat)
      is = ish_at(iat)
      do jat = 1, iat-1
         jzp = id(jat)
         js = ish_at(jat)
         do ish = 1, nsh_id(izp)
            ii = iao_sh(is+ish)
            ni = nao_sh(is+ish)
            do jsh = 1, nsh_id(jzp)
               jj = iao_sh(js+jsh)
               nj = nao_sh(js+jsh)
               do iao = 1, ni
                  do jao = 1, nj
                     tmp = raw(jj+jao, ii+iao) + raw(ii+iao, jj+jao)
                     gradient(jj+jao, ii+iao) = tmp
                     gradient(ii+iao, jj+jao) = tmp
                  end do
               end do
            end do
         end do
      end do
   end do
end subroutine symmetrize_atom_parameter_gradient

!> Adapt direct symmetric-kernel adjoints to the original library convention.
!> Mulliken shell diagonals have a factor -2; other blocks a factor -1.
!> Bond-order atom diagonals have a factor 2; other blocks are unchanged.
!> Physical radial/WS derivative evaluation remains in the public library methods.
subroutine get_mulliken_derivs_direct(self, mol, cache, mulliken_grad, &
   & gradient, sigma)
   class(exchange_fock), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(exchange_cache), intent(inout) :: cache
   real(wp), contiguous, intent(in) :: mulliken_grad(:, :)
   real(wp), contiguous, intent(inout) :: gradient(:, :)
   real(wp), contiguous, intent(inout) :: sigma(:, :)

   integer :: ii, ish, ni
   real(wp), allocatable :: legacy_grad(:, :)

   allocate(legacy_grad(self%nao, self%nao), source=-mulliken_grad)
   do ish = 1, self%nsh
      ii = self%iao_sh(ish)
      ni = self%nao_sh(ish)
      legacy_grad(ii+1:ii+ni, ii+1:ii+ni) = &
         & -2.0_wp*mulliken_grad(ii+1:ii+ni, ii+1:ii+ni)
   end do
   call self%get_mulliken_derivs(mol, cache, legacy_grad, gradient, sigma)
end subroutine get_mulliken_derivs_direct


subroutine get_bocorr_derivs_direct(self, mol, cache, bocorr_grad, gradient, &
   & sigma)
   class(exchange_fock), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(exchange_cache), intent(inout) :: cache
   real(wp), contiguous, intent(in) :: bocorr_grad(:, :)
   real(wp), contiguous, intent(inout) :: gradient(:, :)
   real(wp), contiguous, intent(inout) :: sigma(:, :)

   integer :: iat, ii, is, izp, ni
   real(wp), allocatable :: legacy_grad(:, :)

   allocate(legacy_grad(self%nao, self%nao), source=bocorr_grad)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      is = self%ish_at(iat)
      ii = self%iao_sh(is+1)
      ni = sum(self%nao_sh(is+1:is+self%nsh_id(izp)))
      legacy_grad(ii+1:ii+ni, ii+1:ii+ni) = &
         & 2.0_wp*bocorr_grad(ii+1:ii+ni, ii+1:ii+ni)
   end do
   call self%get_bocorr_derivs(mol, cache, legacy_grad, gradient, sigma)
end subroutine get_bocorr_derivs_direct

end module tblite_exchange_contractions
#endif
