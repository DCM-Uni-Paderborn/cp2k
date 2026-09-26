! SPDX-Identifier: LGPL-3.0-or-later
! Derived from save_tblite consumer cache at aad26158b2e83f84b7b7d5b56f361ead4bc02a88.
! Persistent per-environment CP2K storage, distinct from the live physical SCC cache.
#if defined(__TBLITE_GXTB)
module tblite_stream_plan
   use mctc_env, only: wp
   use mctc_io, only: structure_type
   use tblite_exchange_fock, only: exchange_fock, exchange_bvk_kernel
   implicit none
   private
   public :: cp2k_stream_plan_type
   type :: cp2k_stream_plan_type
      !> Geometry- and mesh-static BvK kernel and Fourier plan
      logical :: bvk_plan_valid = .false.
      !> Number of successfully constructed BvK plans (diagnostic/test oracle)
      integer :: bvk_plan_builds = 0
      !> Exact exchange-model signature associated with the cached kernel
      logical :: bvk_model_valid = .false.
      integer :: bvk_model_nao = 0, bvk_model_nsh = 0, bvk_model_maxsh = 0
      integer, allocatable :: bvk_model_nsh_id(:), bvk_model_nao_sh(:)
      integer, allocatable :: bvk_model_ish_at(:), bvk_model_iao_sh(:)
      real(wp) :: bvk_model_frscale = 0.0_wp
      real(wp) :: bvk_model_onsite_frscale = 0.0_wp
      real(wp) :: bvk_model_omega = 0.0_wp
      real(wp) :: bvk_model_lrscale = 0.0_wp
      real(wp) :: bvk_model_ondiag_scale = 0.0_wp
      real(wp) :: bvk_model_hubbard_exp = 0.0_wp
      real(wp) :: bvk_model_hubbard_exp_r0 = 0.0_wp
      real(wp) :: bvk_model_gexp = 0.0_wp
      real(wp) :: bvk_model_corr_exp = 0.0_wp
      real(wp), allocatable :: bvk_model_hubbard(:, :, :, :)
      real(wp), allocatable :: bvk_model_onecxints(:, :, :)
      real(wp), allocatable :: bvk_model_offdiag_scale(:, :, :, :)
      real(wp), allocatable :: bvk_model_rad(:, :)
      real(wp), allocatable :: bvk_model_kq(:, :)
      real(wp), allocatable :: bvk_model_corr_scale(:, :)
      real(wp), allocatable :: bvk_model_corr_rad(:, :)
      logical :: bvk_periodic(3) = .false.
      integer :: bvk_nmesh(3) = 0
      integer, allocatable :: bvk_id(:)
      real(wp), allocatable :: bvk_xyz(:, :), bvk_lattice(:, :)
      real(wp), allocatable :: bvk_kfrac(:, :), bvk_weights(:)
      real(wp) :: bvk_twist(3) = 0.0_wp
      integer, allocatable :: bvk_input_to_grid(:), bvk_grid_to_input(:)
      !> Compact regular-grid Fourier plan.  `bvk_phase_roots(m+1,idim)`
      !> stores exp(2*pi*i*m/nmesh(idim)); `bvk_twist_phase` stores the
      !> twist character of each real-space representative.  Together these
      !> replace a dense Nk-by-Nk phase table in bounded image batching.
      complex(wp), allocatable :: bvk_phase_roots(:, :)
      complex(wp), allocatable :: bvk_twist_phase(:)
      complex(wp), allocatable :: bvk_phase_forward(:, :)
      complex(wp), allocatable :: bvk_phase_inverse(:, :)
      type(exchange_bvk_kernel) :: bvk_kernel
   contains
      procedure :: bvk_matches
      procedure :: set_bvk_signature
      procedure :: model_matches
      procedure :: set_model_signature
      procedure :: plan_matches
   end type cp2k_stream_plan_type
contains
!> Check whether the persistent BvK plan belongs to this exact geometry and mesh.
logical function bvk_matches(self, mol, nmesh, kfrac, weights, &
   & require_phases) result(matches)
   class(cp2k_stream_plan_type), intent(in) :: self
   type(structure_type), intent(in) :: mol
   integer, intent(in) :: nmesh(3)
   real(wp), intent(in) :: kfrac(:, :), weights(:)
   logical, intent(in), optional :: require_phases

   integer :: ik, nk
   logical :: need_phases

   matches = .false.
   need_phases = .true.
   if (present(require_phases)) need_phases = require_phases
   if (.not.self%bvk_plan_valid) return
   if (.not.self%bvk_model_valid) return
   nk = size(weights)
   if (any(self%bvk_nmesh /= nmesh)) return
   if (any(self%bvk_periodic .neqv. mol%periodic)) return
   if (.not.allocated(self%bvk_id) .or. .not.allocated(self%bvk_xyz) &
      & .or. .not.allocated(self%bvk_lattice) &
      & .or. .not.allocated(self%bvk_kfrac) &
      & .or. .not.allocated(self%bvk_weights) &
      & .or. .not.allocated(self%bvk_input_to_grid) &
      & .or. .not.allocated(self%bvk_grid_to_input) &
      & .or. .not.allocated(self%bvk_phase_roots) &
      & .or. .not.allocated(self%bvk_twist_phase) &
      & .or. .not.allocated(self%bvk_kernel%reps) &
      & .or. .not.allocated(self%bvk_kernel%g_mulliken_r) &
      & .or. .not.allocated(self%bvk_kernel%g_bocorr_r)) return
   if (need_phases) then
      if (.not.allocated(self%bvk_phase_forward) &
         & .or. .not.allocated(self%bvk_phase_inverse)) return
   end if
   if (size(self%bvk_id) /= size(mol%id)) return
   if (any(shape(self%bvk_xyz) /= shape(mol%xyz))) return
   if (any(shape(self%bvk_lattice) /= shape(mol%lattice))) return
   if (any(shape(self%bvk_kfrac) /= shape(kfrac))) return
   if (size(self%bvk_weights) /= size(weights)) return
   if (size(self%bvk_input_to_grid) /= nk) return
   if (size(self%bvk_grid_to_input) /= nk) return
   if (any(shape(self%bvk_phase_roots) /= [maxval(nmesh), 3])) return
   if (size(self%bvk_twist_phase) /= nk) return
   if (any(self%bvk_kernel%nmesh /= nmesh)) return
   if (any(shape(self%bvk_kernel%reps) /= [3, nk])) return
   if (size(self%bvk_kernel%g_mulliken_r, 1) /= &
      & size(self%bvk_kernel%g_mulliken_r, 2)) return
   if (size(self%bvk_kernel%g_mulliken_r, 3) /= nk) return
   if (any(shape(self%bvk_kernel%g_bocorr_r) /= [mol%nat, mol%nat, nk])) return
   if (need_phases) then
      if (any(shape(self%bvk_phase_forward) /= [nk, nk])) return
      if (any(shape(self%bvk_phase_inverse) /= [nk, nk])) return
   end if
   if (any(self%bvk_input_to_grid < 1) &
      & .or. any(self%bvk_input_to_grid > nk)) return
   if (any(self%bvk_grid_to_input < 1) &
      & .or. any(self%bvk_grid_to_input > nk)) return
   do ik = 1, nk
      if (self%bvk_grid_to_input(self%bvk_input_to_grid(ik)) /= ik) return
   end do
   if (any(self%bvk_id /= mol%id)) return
   if (any(self%bvk_xyz /= mol%xyz)) return
   if (any(self%bvk_lattice /= mol%lattice)) return
   if (any(self%bvk_kfrac /= kfrac)) return
   if (any(self%bvk_weights /= weights)) return
   matches = .true.
end function bvk_matches


!> Record the exact fingerprint after a BvK kernel and Fourier plan were validated.
subroutine set_bvk_signature(self, mol, nmesh, kfrac, weights)
   class(cp2k_stream_plan_type), intent(inout) :: self
   type(structure_type), intent(in) :: mol
   integer, intent(in) :: nmesh(3)
   real(wp), intent(in) :: kfrac(:, :), weights(:)

   self%bvk_plan_valid = .false.
   self%bvk_nmesh = nmesh
   self%bvk_periodic = mol%periodic
   if (allocated(self%bvk_id)) deallocate(self%bvk_id)
   if (allocated(self%bvk_xyz)) deallocate(self%bvk_xyz)
   if (allocated(self%bvk_lattice)) deallocate(self%bvk_lattice)
   if (allocated(self%bvk_kfrac)) deallocate(self%bvk_kfrac)
   if (allocated(self%bvk_weights)) deallocate(self%bvk_weights)
   allocate(self%bvk_id, source=mol%id)
   allocate(self%bvk_xyz, source=mol%xyz)
   allocate(self%bvk_lattice, source=mol%lattice)
   allocate(self%bvk_kfrac, source=kfrac)
   allocate(self%bvk_weights, source=weights)
   self%bvk_plan_builds = self%bvk_plan_builds + 1
   self%bvk_plan_valid = .true.
end subroutine set_bvk_signature


logical function model_matches(self, exchange) result(matches)
   class(exchange_fock), intent(in) :: exchange
   class(cp2k_stream_plan_type), intent(in) :: self

   matches = .false.
   if (.not.self%bvk_model_valid) return
   if (.not.allocated(self%bvk_model_nsh_id) &
      & .or. .not.allocated(self%bvk_model_nao_sh) &
      & .or. .not.allocated(self%bvk_model_ish_at) &
      & .or. .not.allocated(self%bvk_model_iao_sh) &
      & .or. .not.allocated(self%bvk_model_hubbard) &
      & .or. .not.allocated(self%bvk_model_onecxints) &
      & .or. .not.allocated(self%bvk_model_offdiag_scale) &
      & .or. .not.allocated(self%bvk_model_rad) &
      & .or. .not.allocated(self%bvk_model_kq) &
      & .or. .not.allocated(self%bvk_model_corr_scale) &
      & .or. .not.allocated(self%bvk_model_corr_rad)) return
   if (self%bvk_model_nao /= exchange%nao &
      & .or. self%bvk_model_nsh /= exchange%nsh &
      & .or. self%bvk_model_maxsh /= exchange%maxsh) return
   if (size(self%bvk_model_nsh_id) /= size(exchange%nsh_id) &
      & .or. size(self%bvk_model_nao_sh) /= size(exchange%nao_sh) &
      & .or. size(self%bvk_model_ish_at) /= size(exchange%ish_at) &
      & .or. size(self%bvk_model_iao_sh) /= size(exchange%iao_sh)) return
   if (any(shape(self%bvk_model_hubbard) /= shape(exchange%hubbard)) &
      & .or. any(shape(self%bvk_model_onecxints) /= shape(exchange%onecxints)) &
      & .or. any(shape(self%bvk_model_offdiag_scale) /= &
      & shape(exchange%offdiag_scale)) &
      & .or. any(shape(self%bvk_model_rad) /= shape(exchange%rad)) &
      & .or. any(shape(self%bvk_model_kq) /= shape(exchange%kq)) &
      & .or. any(shape(self%bvk_model_corr_scale) /= &
      & shape(exchange%corr_scale)) &
      & .or. any(shape(self%bvk_model_corr_rad) /= &
      & shape(exchange%corr_rad))) return
   if (any(self%bvk_model_nsh_id /= exchange%nsh_id) &
      & .or. any(self%bvk_model_nao_sh /= exchange%nao_sh) &
      & .or. any(self%bvk_model_ish_at /= exchange%ish_at) &
      & .or. any(self%bvk_model_iao_sh /= exchange%iao_sh)) return
   if (self%bvk_model_frscale /= exchange%frscale &
      & .or. self%bvk_model_onsite_frscale /= exchange%onsite_frscale &
      & .or. self%bvk_model_omega /= exchange%omega &
      & .or. self%bvk_model_lrscale /= exchange%lrscale &
      & .or. self%bvk_model_ondiag_scale /= exchange%ondiag_scale &
      & .or. self%bvk_model_hubbard_exp /= exchange%hubbard_exp &
      & .or. self%bvk_model_hubbard_exp_r0 /= exchange%hubbard_exp_r0 &
      & .or. self%bvk_model_gexp /= exchange%gexp &
      & .or. self%bvk_model_corr_exp /= exchange%corr_exp) return
   if (any(self%bvk_model_hubbard /= exchange%hubbard) &
      & .or. any(self%bvk_model_onecxints /= exchange%onecxints) &
      & .or. any(self%bvk_model_offdiag_scale /= exchange%offdiag_scale) &
      & .or. any(self%bvk_model_rad /= exchange%rad) &
      & .or. any(self%bvk_model_kq /= exchange%kq) &
      & .or. any(self%bvk_model_corr_scale /= exchange%corr_scale) &
      & .or. any(self%bvk_model_corr_rad /= exchange%corr_rad)) return
   matches = .true.
end function model_matches

subroutine set_model_signature(self, exchange)
   class(exchange_fock), intent(in) :: exchange
   class(cp2k_stream_plan_type), intent(inout) :: self

   self%bvk_model_valid = .false.
   self%bvk_model_nao = exchange%nao
   self%bvk_model_nsh = exchange%nsh
   self%bvk_model_maxsh = exchange%maxsh
   self%bvk_model_frscale = exchange%frscale
   self%bvk_model_onsite_frscale = exchange%onsite_frscale
   self%bvk_model_omega = exchange%omega
   self%bvk_model_lrscale = exchange%lrscale
   self%bvk_model_ondiag_scale = exchange%ondiag_scale
   self%bvk_model_hubbard_exp = exchange%hubbard_exp
   self%bvk_model_hubbard_exp_r0 = exchange%hubbard_exp_r0
   self%bvk_model_gexp = exchange%gexp
   self%bvk_model_corr_exp = exchange%corr_exp
   if (allocated(self%bvk_model_nsh_id)) deallocate(self%bvk_model_nsh_id)
   if (allocated(self%bvk_model_nao_sh)) deallocate(self%bvk_model_nao_sh)
   if (allocated(self%bvk_model_ish_at)) deallocate(self%bvk_model_ish_at)
   if (allocated(self%bvk_model_iao_sh)) deallocate(self%bvk_model_iao_sh)
   if (allocated(self%bvk_model_hubbard)) deallocate(self%bvk_model_hubbard)
   if (allocated(self%bvk_model_onecxints)) &
      & deallocate(self%bvk_model_onecxints)
   if (allocated(self%bvk_model_offdiag_scale)) &
      & deallocate(self%bvk_model_offdiag_scale)
   if (allocated(self%bvk_model_rad)) deallocate(self%bvk_model_rad)
   if (allocated(self%bvk_model_kq)) deallocate(self%bvk_model_kq)
   if (allocated(self%bvk_model_corr_scale)) &
      & deallocate(self%bvk_model_corr_scale)
   if (allocated(self%bvk_model_corr_rad)) deallocate(self%bvk_model_corr_rad)
   allocate(self%bvk_model_nsh_id, source=exchange%nsh_id)
   allocate(self%bvk_model_nao_sh, source=exchange%nao_sh)
   allocate(self%bvk_model_ish_at, source=exchange%ish_at)
   allocate(self%bvk_model_iao_sh, source=exchange%iao_sh)
   allocate(self%bvk_model_hubbard, source=exchange%hubbard)
   allocate(self%bvk_model_onecxints, source=exchange%onecxints)
   allocate(self%bvk_model_offdiag_scale, source=exchange%offdiag_scale)
   allocate(self%bvk_model_rad, source=exchange%rad)
   allocate(self%bvk_model_kq, source=exchange%kq)
   allocate(self%bvk_model_corr_scale, source=exchange%corr_scale)
   allocate(self%bvk_model_corr_rad, source=exchange%corr_rad)
   self%bvk_model_valid = .true.
end subroutine set_model_signature

logical function plan_matches(self, exchange, mol, kernel, kfrac, weights, &
   & require_phases) result(matches)
   class(exchange_fock), intent(in) :: exchange
   type(structure_type), intent(in) :: mol
   class(cp2k_stream_plan_type), intent(in) :: self
   type(exchange_bvk_kernel), intent(in) :: kernel
   real(wp), intent(in) :: kfrac(:, :), weights(:)
   logical, intent(in), optional :: require_phases

   logical :: need_phases

   matches = .false.
   need_phases = .true.
   if (present(require_phases)) need_phases = require_phases
   if (.not.self%bvk_matches(mol, kernel%nmesh, kfrac, weights, &
      & require_phases=need_phases)) return
   if (.not.self%model_matches(exchange)) return
   if (.not.allocated(kernel%reps)) return
   if (any(shape(kernel%reps) /= shape(self%bvk_kernel%reps))) return
   if (any(kernel%reps /= self%bvk_kernel%reps)) return
   matches = .true.
end function plan_matches

end module tblite_stream_plan
#endif
