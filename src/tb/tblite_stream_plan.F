!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: LGPL-3.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
! SPDX-Identifier: LGPL-3.0-or-later
! Derived from save_tblite consumer cache at aad26158b2e83f84b7b7d5b56f361ead4bc02a88.
! Persistent per-environment CP2K storage, distinct from the live physical SCC cache.
#if defined(__TBLITE_GXTB)
MODULE tblite_stream_plan
   USE mctc_env,                        ONLY: wp
   USE mctc_io,                         ONLY: structure_type
   USE tblite_exchange_fock,            ONLY: exchange_bvk_kernel,&
                                              exchange_fock

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: cp2k_stream_plan_type
   TYPE :: cp2k_stream_plan_type
      !> Geometry- and mesh-static BvK kernel and Fourier plan
      LOGICAL :: bvk_plan_valid = .FALSE.
      !> Number of successfully constructed BvK plans (diagnostic/test oracle)
      INTEGER :: bvk_plan_builds = 0
      !> Exact exchange-model signature associated with the cached kernel
      LOGICAL :: bvk_model_valid = .FALSE.
      INTEGER :: bvk_model_nao = 0, bvk_model_nsh = 0, bvk_model_maxsh = 0
      INTEGER, ALLOCATABLE :: bvk_model_nsh_id(:), bvk_model_nao_sh(:)
      INTEGER, ALLOCATABLE :: bvk_model_ish_at(:), bvk_model_iao_sh(:)
      REAL(wp) :: bvk_model_frscale = 0.0_wp
      REAL(wp) :: bvk_model_onsite_frscale = 0.0_wp
      REAL(wp) :: bvk_model_omega = 0.0_wp
      REAL(wp) :: bvk_model_lrscale = 0.0_wp
      REAL(wp) :: bvk_model_ondiag_scale = 0.0_wp
      REAL(wp) :: bvk_model_hubbard_exp = 0.0_wp
      REAL(wp) :: bvk_model_hubbard_exp_r0 = 0.0_wp
      REAL(wp) :: bvk_model_gexp = 0.0_wp
      REAL(wp) :: bvk_model_corr_exp = 0.0_wp
      REAL(wp), ALLOCATABLE :: bvk_model_hubbard(:, :, :, :)
      REAL(wp), ALLOCATABLE :: bvk_model_onecxints(:, :, :)
      REAL(wp), ALLOCATABLE :: bvk_model_offdiag_scale(:, :, :, :)
      REAL(wp), ALLOCATABLE :: bvk_model_rad(:, :)
      REAL(wp), ALLOCATABLE :: bvk_model_kq(:, :)
      REAL(wp), ALLOCATABLE :: bvk_model_corr_scale(:, :)
      REAL(wp), ALLOCATABLE :: bvk_model_corr_rad(:, :)
      LOGICAL :: bvk_periodic(3) = .FALSE.
      INTEGER :: bvk_nmesh(3) = 0
      INTEGER, ALLOCATABLE :: bvk_id(:)
      REAL(wp), ALLOCATABLE :: bvk_xyz(:, :), bvk_lattice(:, :)
      REAL(wp), ALLOCATABLE :: bvk_kfrac(:, :), bvk_weights(:)
      REAL(wp) :: bvk_twist(3) = 0.0_wp
      INTEGER, ALLOCATABLE :: bvk_input_to_grid(:), bvk_grid_to_input(:)
      !> Compact regular-grid Fourier plan.  `bvk_phase_roots(m+1,idim)`
      !> stores exp(2*pi*i*m/nmesh(idim)); `bvk_twist_phase` stores the
      !> twist character of each real-space representative.  Together these
      !> replace a dense Nk-by-Nk phase table in bounded image batching.
      COMPLEX(wp), ALLOCATABLE :: bvk_phase_roots(:, :)
      COMPLEX(wp), ALLOCATABLE :: bvk_twist_phase(:)
      COMPLEX(wp), ALLOCATABLE :: bvk_phase_forward(:, :)
      COMPLEX(wp), ALLOCATABLE :: bvk_phase_inverse(:, :)
      TYPE(exchange_bvk_kernel) :: bvk_kernel
   CONTAINS
      PROCEDURE :: bvk_matches
      PROCEDURE :: set_bvk_signature
      PROCEDURE :: model_matches
      PROCEDURE :: set_model_signature
      PROCEDURE :: plan_matches
   END TYPE cp2k_stream_plan_type
CONTAINS
! **************************************************************************************************
!> \brief ...
!> \param self ...
!> \param mol ...
!> \param nmesh ...
!> \param kfrac ...
!> \param weights ...
!> \param require_phases ...
!> \return ...
!> \note Check whether the persistent BvK plan belongs to this exact geometry and mesh.
! **************************************************************************************************
   LOGICAL FUNCTION bvk_matches(self, mol, nmesh, kfrac, weights, &
      & require_phases) RESULT(matches)
      class(cp2k_stream_plan_type), INTENT(in) :: self
      TYPE(structure_type), INTENT(in) :: mol
      INTEGER, INTENT(in) :: nmesh(3)
      REAL(wp), INTENT(in) :: kfrac(:, :), weights(:)
      LOGICAL, INTENT(in), OPTIONAL :: require_phases

      INTEGER :: ik, nk
      LOGICAL :: need_phases

      matches = .FALSE.
      need_phases = .TRUE.
      IF (PRESENT(require_phases)) need_phases = require_phases
      IF (.NOT. self%bvk_plan_valid) RETURN
      IF (.NOT. self%bvk_model_valid) RETURN
      nk = SIZE(weights)
      IF (ANY(self%bvk_nmesh /= nmesh)) RETURN
      IF (ANY(self%bvk_periodic .NEQV. mol%periodic)) RETURN
      IF (.NOT. ALLOCATED(self%bvk_id) .OR. .NOT. ALLOCATED(self%bvk_xyz) &
         & .OR. .NOT. ALLOCATED(self%bvk_lattice) &
         & .OR. .NOT. ALLOCATED(self%bvk_kfrac) &
         & .OR. .NOT. ALLOCATED(self%bvk_weights) &
         & .OR. .NOT. ALLOCATED(self%bvk_input_to_grid) &
         & .OR. .NOT. ALLOCATED(self%bvk_grid_to_input) &
         & .OR. .NOT. ALLOCATED(self%bvk_phase_roots) &
         & .OR. .NOT. ALLOCATED(self%bvk_twist_phase) &
         & .OR. .NOT. ALLOCATED(self%bvk_kernel%reps) &
         & .OR. .NOT. ALLOCATED(self%bvk_kernel%g_mulliken_r) &
         & .OR. .NOT. ALLOCATED(self%bvk_kernel%g_bocorr_r)) RETURN
      IF (need_phases) THEN
         IF (.NOT. ALLOCATED(self%bvk_phase_forward) &
            & .OR. .NOT. ALLOCATED(self%bvk_phase_inverse)) RETURN
      END IF
      IF (SIZE(self%bvk_id) /= SIZE(mol%id)) RETURN
      IF (ANY(SHAPE(self%bvk_xyz) /= SHAPE(mol%xyz))) RETURN
      IF (ANY(SHAPE(self%bvk_lattice) /= SHAPE(mol%lattice))) RETURN
      IF (ANY(SHAPE(self%bvk_kfrac) /= SHAPE(kfrac))) RETURN
      IF (SIZE(self%bvk_weights) /= SIZE(weights)) RETURN
      IF (SIZE(self%bvk_input_to_grid) /= nk) RETURN
      IF (SIZE(self%bvk_grid_to_input) /= nk) RETURN
      IF (ANY(SHAPE(self%bvk_phase_roots) /= [MAXVAL(nmesh), 3])) RETURN
      IF (SIZE(self%bvk_twist_phase) /= nk) RETURN
      IF (ANY(self%bvk_kernel%nmesh /= nmesh)) RETURN
      IF (ANY(SHAPE(self%bvk_kernel%reps) /= [3, nk])) RETURN
      IF (SIZE(self%bvk_kernel%g_mulliken_r, 1) /= &
         & SIZE(self%bvk_kernel%g_mulliken_r, 2)) RETURN
      IF (SIZE(self%bvk_kernel%g_mulliken_r, 3) /= nk) RETURN
      IF (ANY(SHAPE(self%bvk_kernel%g_bocorr_r) /= [mol%nat, mol%nat, nk])) RETURN
      IF (need_phases) THEN
         IF (ANY(SHAPE(self%bvk_phase_forward) /= [nk, nk])) RETURN
         IF (ANY(SHAPE(self%bvk_phase_inverse) /= [nk, nk])) RETURN
      END IF
      IF (ANY(self%bvk_input_to_grid < 1) &
         & .OR. ANY(self%bvk_input_to_grid > nk)) RETURN
      IF (ANY(self%bvk_grid_to_input < 1) &
         & .OR. ANY(self%bvk_grid_to_input > nk)) RETURN
      DO ik = 1, nk
         IF (self%bvk_grid_to_input(self%bvk_input_to_grid(ik)) /= ik) RETURN
      END DO
      IF (ANY(self%bvk_id /= mol%id)) RETURN
      IF (ANY(self%bvk_xyz /= mol%xyz)) RETURN
      IF (ANY(self%bvk_lattice /= mol%lattice)) RETURN
      IF (ANY(self%bvk_kfrac /= kfrac)) RETURN
      IF (ANY(self%bvk_weights /= weights)) RETURN
      matches = .TRUE.
   END FUNCTION bvk_matches

! **************************************************************************************************
!> \brief ...
!> \param self ...
!> \param mol ...
!> \param nmesh ...
!> \param kfrac ...
!> \param weights ...
!> \note Record the exact fingerprint after a BvK kernel and Fourier plan were validated.
! **************************************************************************************************
   SUBROUTINE set_bvk_signature(self, mol, nmesh, kfrac, weights)
      class(cp2k_stream_plan_type), INTENT(inout) :: self
      TYPE(structure_type), INTENT(in) :: mol
      INTEGER, INTENT(in) :: nmesh(3)
      REAL(wp), INTENT(in) :: kfrac(:, :), weights(:)

      self%bvk_plan_valid = .FALSE.
      self%bvk_nmesh = nmesh
      self%bvk_periodic = mol%periodic
      IF (ALLOCATED(self%bvk_id)) DEALLOCATE (self%bvk_id)
      IF (ALLOCATED(self%bvk_xyz)) DEALLOCATE (self%bvk_xyz)
      IF (ALLOCATED(self%bvk_lattice)) DEALLOCATE (self%bvk_lattice)
      IF (ALLOCATED(self%bvk_kfrac)) DEALLOCATE (self%bvk_kfrac)
      IF (ALLOCATED(self%bvk_weights)) DEALLOCATE (self%bvk_weights)
      ALLOCATE (self%bvk_id, source=mol%id)
      ALLOCATE (self%bvk_xyz, source=mol%xyz)
      ALLOCATE (self%bvk_lattice, source=mol%lattice)
      ALLOCATE (self%bvk_kfrac, source=kfrac)
      ALLOCATE (self%bvk_weights, source=weights)
      self%bvk_plan_builds = self%bvk_plan_builds + 1
      self%bvk_plan_valid = .TRUE.
   END SUBROUTINE set_bvk_signature

! **************************************************************************************************
!> \brief ...
!> \param self ...
!> \param exchange ...
!> \return ...
! **************************************************************************************************
   LOGICAL FUNCTION model_matches(self, exchange) RESULT(matches)
      class(exchange_fock), INTENT(in) :: exchange
      class(cp2k_stream_plan_type), INTENT(in) :: self

      matches = .FALSE.
      IF (.NOT. self%bvk_model_valid) RETURN
      IF (.NOT. ALLOCATED(self%bvk_model_nsh_id) &
         & .OR. .NOT. ALLOCATED(self%bvk_model_nao_sh) &
         & .OR. .NOT. ALLOCATED(self%bvk_model_ish_at) &
         & .OR. .NOT. ALLOCATED(self%bvk_model_iao_sh) &
         & .OR. .NOT. ALLOCATED(self%bvk_model_hubbard) &
         & .OR. .NOT. ALLOCATED(self%bvk_model_onecxints) &
         & .OR. .NOT. ALLOCATED(self%bvk_model_offdiag_scale) &
         & .OR. .NOT. ALLOCATED(self%bvk_model_rad) &
         & .OR. .NOT. ALLOCATED(self%bvk_model_kq) &
         & .OR. .NOT. ALLOCATED(self%bvk_model_corr_scale) &
         & .OR. .NOT. ALLOCATED(self%bvk_model_corr_rad)) RETURN
      IF (self%bvk_model_nao /= exchange%nao &
         & .OR. self%bvk_model_nsh /= exchange%nsh &
         & .OR. self%bvk_model_maxsh /= exchange%maxsh) RETURN
      IF (SIZE(self%bvk_model_nsh_id) /= SIZE(exchange%nsh_id) &
         & .OR. SIZE(self%bvk_model_nao_sh) /= SIZE(exchange%nao_sh) &
         & .OR. SIZE(self%bvk_model_ish_at) /= SIZE(exchange%ish_at) &
         & .OR. SIZE(self%bvk_model_iao_sh) /= SIZE(exchange%iao_sh)) RETURN
      IF (ANY(SHAPE(self%bvk_model_hubbard) /= SHAPE(exchange%hubbard)) &
         & .OR. ANY(SHAPE(self%bvk_model_onecxints) /= SHAPE(exchange%onecxints)) &
         & .OR. ANY(SHAPE(self%bvk_model_offdiag_scale) /= &
         & SHAPE(exchange%offdiag_scale)) &
         & .OR. ANY(SHAPE(self%bvk_model_rad) /= SHAPE(exchange%rad)) &
         & .OR. ANY(SHAPE(self%bvk_model_kq) /= SHAPE(exchange%kq)) &
         & .OR. ANY(SHAPE(self%bvk_model_corr_scale) /= &
         & SHAPE(exchange%corr_scale)) &
         & .OR. ANY(SHAPE(self%bvk_model_corr_rad) /= &
         & SHAPE(exchange%corr_rad))) RETURN
      IF (ANY(self%bvk_model_nsh_id /= exchange%nsh_id) &
         & .OR. ANY(self%bvk_model_nao_sh /= exchange%nao_sh) &
         & .OR. ANY(self%bvk_model_ish_at /= exchange%ish_at) &
         & .OR. ANY(self%bvk_model_iao_sh /= exchange%iao_sh)) RETURN
      IF (self%bvk_model_frscale /= exchange%frscale &
         & .OR. self%bvk_model_onsite_frscale /= exchange%onsite_frscale &
         & .OR. self%bvk_model_omega /= exchange%omega &
         & .OR. self%bvk_model_lrscale /= exchange%lrscale &
         & .OR. self%bvk_model_ondiag_scale /= exchange%ondiag_scale &
         & .OR. self%bvk_model_hubbard_exp /= exchange%hubbard_exp &
         & .OR. self%bvk_model_hubbard_exp_r0 /= exchange%hubbard_exp_r0 &
         & .OR. self%bvk_model_gexp /= exchange%gexp &
         & .OR. self%bvk_model_corr_exp /= exchange%corr_exp) RETURN
      IF (ANY(self%bvk_model_hubbard /= exchange%hubbard) &
         & .OR. ANY(self%bvk_model_onecxints /= exchange%onecxints) &
         & .OR. ANY(self%bvk_model_offdiag_scale /= exchange%offdiag_scale) &
         & .OR. ANY(self%bvk_model_rad /= exchange%rad) &
         & .OR. ANY(self%bvk_model_kq /= exchange%kq) &
         & .OR. ANY(self%bvk_model_corr_scale /= exchange%corr_scale) &
         & .OR. ANY(self%bvk_model_corr_rad /= exchange%corr_rad)) RETURN
      matches = .TRUE.
   END FUNCTION model_matches

! **************************************************************************************************
!> \brief ...
!> \param self ...
!> \param exchange ...
! **************************************************************************************************
   SUBROUTINE set_model_signature(self, exchange)
      class(exchange_fock), INTENT(in) :: exchange
      class(cp2k_stream_plan_type), INTENT(inout) :: self

      self%bvk_model_valid = .FALSE.
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
      IF (ALLOCATED(self%bvk_model_nsh_id)) DEALLOCATE (self%bvk_model_nsh_id)
      IF (ALLOCATED(self%bvk_model_nao_sh)) DEALLOCATE (self%bvk_model_nao_sh)
      IF (ALLOCATED(self%bvk_model_ish_at)) DEALLOCATE (self%bvk_model_ish_at)
      IF (ALLOCATED(self%bvk_model_iao_sh)) DEALLOCATE (self%bvk_model_iao_sh)
      IF (ALLOCATED(self%bvk_model_hubbard)) DEALLOCATE (self%bvk_model_hubbard)
      IF (ALLOCATED(self%bvk_model_onecxints)) THEN
         DEALLOCATE (self%bvk_model_onecxints)
      END IF
      IF (ALLOCATED(self%bvk_model_offdiag_scale)) THEN
         DEALLOCATE (self%bvk_model_offdiag_scale)
      END IF
      IF (ALLOCATED(self%bvk_model_rad)) DEALLOCATE (self%bvk_model_rad)
      IF (ALLOCATED(self%bvk_model_kq)) DEALLOCATE (self%bvk_model_kq)
      IF (ALLOCATED(self%bvk_model_corr_scale)) THEN
         DEALLOCATE (self%bvk_model_corr_scale)
      END IF
      IF (ALLOCATED(self%bvk_model_corr_rad)) DEALLOCATE (self%bvk_model_corr_rad)
      ALLOCATE (self%bvk_model_nsh_id, source=exchange%nsh_id)
      ALLOCATE (self%bvk_model_nao_sh, source=exchange%nao_sh)
      ALLOCATE (self%bvk_model_ish_at, source=exchange%ish_at)
      ALLOCATE (self%bvk_model_iao_sh, source=exchange%iao_sh)
      ALLOCATE (self%bvk_model_hubbard, source=exchange%hubbard)
      ALLOCATE (self%bvk_model_onecxints, source=exchange%onecxints)
      ALLOCATE (self%bvk_model_offdiag_scale, source=exchange%offdiag_scale)
      ALLOCATE (self%bvk_model_rad, source=exchange%rad)
      ALLOCATE (self%bvk_model_kq, source=exchange%kq)
      ALLOCATE (self%bvk_model_corr_scale, source=exchange%corr_scale)
      ALLOCATE (self%bvk_model_corr_rad, source=exchange%corr_rad)
      self%bvk_model_valid = .TRUE.
   END SUBROUTINE set_model_signature

! **************************************************************************************************
!> \brief ...
!> \param self ...
!> \param exchange ...
!> \param mol ...
!> \param kernel ...
!> \param kfrac ...
!> \param weights ...
!> \param require_phases ...
!> \return ...
! **************************************************************************************************
   LOGICAL FUNCTION plan_matches(self, exchange, mol, kernel, kfrac, weights, &
      & require_phases) RESULT(matches)
      class(exchange_fock), INTENT(in) :: exchange
      TYPE(structure_type), INTENT(in) :: mol
      class(cp2k_stream_plan_type), INTENT(in) :: self
      TYPE(exchange_bvk_kernel), INTENT(in) :: kernel
      REAL(wp), INTENT(in) :: kfrac(:, :), weights(:)
      LOGICAL, INTENT(in), OPTIONAL :: require_phases

      LOGICAL :: need_phases

      matches = .FALSE.
      need_phases = .TRUE.
      IF (PRESENT(require_phases)) need_phases = require_phases
      IF (.NOT. self%bvk_matches(mol, kernel%nmesh, kfrac, weights, &
         & require_phases=need_phases)) RETURN
      IF (.NOT. self%model_matches(exchange)) RETURN
      IF (.NOT. ALLOCATED(kernel%reps)) RETURN
      IF (ANY(SHAPE(kernel%reps) /= SHAPE(self%bvk_kernel%reps))) RETURN
      IF (ANY(kernel%reps /= self%bvk_kernel%reps)) RETURN
      matches = .TRUE.
   END FUNCTION plan_matches

END MODULE tblite_stream_plan
#endif
