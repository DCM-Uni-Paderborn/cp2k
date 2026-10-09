!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: LGPL-3.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
! SPDX-Identifier: LGPL-3.0-or-later
! Derived from save_tblite mesh contractions at aad26158b2e83f84b7b7d5b56f361ead4bc02a88.
! CP2K owns mesh transforms, matrix contractions, buffering and their adjoints.
! Radial kernels and their geometry derivatives remain in the physical provider.
#if defined(__TBLITE_GXTB)
MODULE tblite_exchange_contractions
   USE mctc_env,                        ONLY: wp
   USE mctc_io,                         ONLY: structure_type
   USE mctc_io_constants,               ONLY: pi
   USE tblite_bvk_transform,            ONLY: bvk_mixed_radix_k_to_r,&
                                              bvk_mixed_radix_r_to_k,&
                                              bvk_separable_k_to_r,&
                                              bvk_separable_r_to_k,&
                                              bvk_transform_dense,&
                                              bvk_transform_mixed_radix_fft,&
                                              bvk_transform_separable_dft
   USE tblite_exchange_cache,           ONLY: exchange_cache
   USE tblite_exchange_fock,            ONLY: exchange_bvk_kernel,&
                                              exchange_fock
   USE tblite_stream_plan,              ONLY: cp2k_stream_plan_type
   USE tblite_transform_contract,       ONLY: exchange_k_to_r,&
                                              exchange_k_to_r_adjoint,&
                                              exchange_r_to_k,&
                                              exchange_r_to_k_adjoint,&
                                              exchange_transform_type

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: get_kfock_kpoint, get_kgrad_kpoint
   PUBLIC :: get_kfock_kmesh
   PUBLIC :: get_kgrad_kmesh
   PUBLIC :: get_kfock_stream_apply
   PUBLIC :: get_kfock_batch_apply
   PUBLIC :: get_kgrad_batch_apply
   PUBLIC :: get_mulliken_derivs_direct, get_bocorr_derivs_direct
CONTAINS
! **************************************************************************************************
!> \brief ...
!> \param backend ...
!> \return ...
! **************************************************************************************************
   pure LOGICAL FUNCTION is_compact_transform_backend(backend) RESULT(compact)
      INTEGER, INTENT(in)                                :: backend

      compact = backend == bvk_transform_separable_dft &
         & .OR. backend == bvk_transform_mixed_radix_fft
   END FUNCTION is_compact_transform_backend

! **************************************************************************************************
!> \brief ...
!> \param amat ...
!> \param bmat ...
!> \param cmat ...
! **************************************************************************************************
   pure SUBROUTINE matmul_complex_contiguous(amat, bmat, cmat)
      COMPLEX(wp), contiguous, INTENT(in)                :: amat(:, :), bmat(:, :)
      COMPLEX(wp), contiguous, INTENT(out)               :: cmat(:, :)

      cmat = MATMUL(amat, bmat)
   END SUBROUTINE matmul_complex_contiguous

! **************************************************************************************************
!> \brief ...
!> \param self ...
!> \param mol ...
!> \param cache ...
!> \param plan ...
!> \param kernel ...
!> \param kfrac ...
!> \param weights ...
!> \param density ...
!> \param overlap ...
!> \param fock ...
!> \param vsh ...
!> \param energy ...
!> \param transform_backend ...
!> \param transform ...
! **************************************************************************************************
   SUBROUTINE get_KFock_kmesh(self, mol, cache, plan, kernel, kfrac, weights, &
      & density, overlap, fock, vsh, energy, transform_backend, transform)
      class(exchange_fock), INTENT(in) :: self
      TYPE(structure_type), INTENT(in) :: mol
      TYPE(exchange_cache), INTENT(in) :: cache
      TYPE(cp2k_stream_plan_type), INTENT(in), TARGET :: plan
      TYPE(exchange_bvk_kernel), INTENT(in) :: kernel
      REAL(wp), INTENT(in) :: kfrac(:, :), weights(:)
      COMPLEX(wp), INTENT(in) :: density(:, :, :, :), overlap(:, :, :)
      COMPLEX(wp), INTENT(out) :: fock(:, :, :, :)
      REAL(wp), INTENT(out) :: vsh(:), energy
      INTEGER, INTENT(in), OPTIONAL :: transform_backend
      class(exchange_transform_type), INTENT(inout), OPTIONAL :: transform

      INTEGER :: spin, ik, backend

      backend = bvk_transform_dense
      IF (PRESENT(transform_backend)) backend = transform_backend

      IF (PRESENT(transform)) THEN
         CALL build_KFock_kmesh(self, mol, density, overlap, kfrac, weights, &
            & kernel, cache%g_onsfx, cache%g_onsri, fock, .TRUE., &
            & transform=transform)
      ELSE IF (is_compact_transform_backend(backend) .AND. &
         & plan%plan_matches(self, mol, kernel, kfrac, weights, &
         & require_phases=.FALSE.) .AND. &
         & ALLOCATED(plan%bvk_phase_roots) .AND. &
         & ALLOCATED(plan%bvk_twist_phase)) THEN
         CALL build_KFock_kmesh(self, mol, density, overlap, kfrac, weights, &
            & kernel, cache%g_onsfx, cache%g_onsri, fock, .TRUE., &
            & transform_backend=backend, nmesh=plan%bvk_nmesh, &
            & input_to_grid=plan%bvk_input_to_grid, &
            & phase_roots=plan%bvk_phase_roots, &
            & twist_phase=plan%bvk_twist_phase, &
            & transform=transform)
      ELSE IF (plan%plan_matches(self, mol, kernel, kfrac, weights)) THEN
         CALL build_KFock_kmesh(self, mol, density, overlap, kfrac, weights, &
            & kernel, cache%g_onsfx, cache%g_onsri, fock, .TRUE., &
            & plan%bvk_phase_forward, plan%bvk_phase_inverse, &
            & transform=transform)
      ELSE
         CALL build_KFock_kmesh(self, mol, density, overlap, kfrac, weights, &
            & kernel, cache%g_onsfx, cache%g_onsri, fock, .TRUE., &
            & transform=transform)
      END IF

      energy = 0.0_wp
      DO ik = 1, SIZE(weights)
         DO spin = 1, SIZE(density, 3)
            energy = energy + 0.5_wp*weights(ik)*REAL(SUM( &
               & CONJG(density(:, :, spin, ik))*fock(:, :, spin, ik)), wp)
         END DO
      END DO

      CALL get_KFock_kmesh_vsh(self, mol, cache, weights, density, overlap, vsh)

   END SUBROUTINE get_KFock_kmesh

! **************************************************************************************************
!> \brief ...
!> \param self ...
!> \param mol ...
!> \param cache ...
!> \param weights ...
!> \param density ...
!> \param overlap ...
!> \param vsh ...
! **************************************************************************************************
   SUBROUTINE get_KFock_kmesh_vsh(self, mol, cache, weights, density, overlap, vsh)
      class(exchange_fock), INTENT(in) :: self
      TYPE(structure_type), INTENT(in) :: mol
      TYPE(exchange_cache), INTENT(in) :: cache
      REAL(wp), INTENT(in) :: weights(:)
      COMPLEX(wp), INTENT(in) :: density(:, :, :, :), overlap(:, :, :)
      REAL(wp), INTENT(out) :: vsh(:)

      INTEGER :: iao, iat, ik, iq, is, ish, izp, spin
      REAL(wp) :: spin_factor
      REAL(wp), ALLOCATABLE :: grad_onsfx(:, :, :), grad_onsri(:, :)
      COMPLEX(wp), ALLOCATABLE :: bq(:, :), bw(:, :), bt0(:, :), bu0(:, :), &
         & bv0(:, :), d0(:, :), dmat(:, :), p0(:, :), x0(:, :), xmat(:, :)
      COMPLEX(wp), ALLOCATABLE :: bgp(:), bgsp(:), bgsps(:), diag2d(:), &
         & diagp(:), diagx(:)

      ALLOCATE (grad_onsfx(self%maxsh, self%maxsh, mol%nat), &
         & grad_onsri(self%maxsh, mol%nat), source=0.0_wp)
      ALLOCATE (bq(self%nao, self%nao), bw(self%nao, self%nao), &
         & bt0(self%nao, self%nao), bu0(self%nao, self%nao), &
         & bv0(self%nao, self%nao), d0(self%nao, self%nao), &
         & dmat(self%nao, self%nao), p0(self%nao, self%nao), &
         & x0(self%nao, self%nao), xmat(self%nao, self%nao), &
         & bgp(self%nao), bgsp(self%nao), bgsps(self%nao), &
         & diag2d(self%nao), diagp(self%nao), diagx(self%nao), &
         & source=(0.0_wp, 0.0_wp))

      spin_factor = 0.5_wp
      IF (SIZE(density, 3) > 1) spin_factor = 1.0_wp

      DO spin = 1, SIZE(density, 3)
         p0 = (0.0_wp, 0.0_wp)
         x0 = (0.0_wp, 0.0_wp)
         d0 = (0.0_wp, 0.0_wp)
         DO ik = 1, SIZE(weights)
            xmat = MATMUL(overlap(:, :, ik), density(:, :, spin, ik))
            dmat = 0.5_wp*MATMUL(xmat, overlap(:, :, ik))
            p0 = p0 + weights(ik)*density(:, :, spin, ik)
            x0 = x0 + weights(ik)*xmat
            d0 = d0 + weights(ik)*dmat
         END DO

         bt0 = (0.0_wp, 0.0_wp)
         bu0 = (0.0_wp, 0.0_wp)
         bv0 = (0.0_wp, 0.0_wp)
         bgp = (0.0_wp, 0.0_wp)
         bgsp = (0.0_wp, 0.0_wp)
         bgsps = (0.0_wp, 0.0_wp)
         DO ik = 1, SIZE(weights)
            bq = -0.25_wp*spin_factor*weights(ik)*density(:, :, spin, ik)
            bu0 = bu0 + bq
            bw = MATMUL(bq, CONJG(TRANSPOSE(overlap(:, :, ik))))
            bt0 = bt0 + bw
            bv0 = bv0 + 0.5_wp*MATMUL( &
               & CONJG(TRANSPOSE(overlap(:, :, ik))), bw)
            DO iao = 1, self%nao
               bgp(iao) = bgp(iao) + 0.25_wp*SUM( &
                  & CONJG(overlap(:, iao, ik))*bw(:, iao))
               bgsp(iao) = bgsp(iao) + 0.5_wp*SUM( &
                  & CONJG(overlap(:, iao, ik))*bq(:, iao))
               bgsps(iao) = bgsps(iao) + 0.25_wp*bq(iao, iao)
            END DO
         END DO

         CALL onsite_fx_parameter_gradient_complex(mol%nat, mol%id, &
            & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bt0, x0, &
            & 0.5_wp, grad_onsfx, adjoint_src=.TRUE.)
         CALL onsite_fx_parameter_gradient_complex(mol%nat, mol%id, &
            & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bu0, d0, &
            & 0.5_wp, grad_onsfx)
         CALL onsite_fx_parameter_gradient_complex(mol%nat, mol%id, &
            & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bv0, p0, &
            & 0.5_wp, grad_onsfx)

         CALL onsite_ri_parameter_gradient_complex(mol%nat, mol%id, &
            & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bt0, x0, &
            & -2.0_wp, grad_onsri)
         CALL onsite_ri_parameter_gradient_complex(mol%nat, mol%id, &
            & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bu0, d0, &
            & -2.0_wp, grad_onsri)
         CALL onsite_ri_parameter_gradient_complex(mol%nat, mol%id, &
            & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bv0, p0, &
            & -2.0_wp, grad_onsri)

         DO iao = 1, self%nao
            diagp(iao) = p0(iao, iao)
            diagx(iao) = x0(iao, iao)
            diag2d(iao) = 2.0_wp*d0(iao, iao)
         END DO
         CALL onsite_fx_symv_parameter_gradient_complex(mol%nat, mol%id, &
            & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bgp, &
            & diagp, grad_onsfx)
         CALL onsite_fx_symv_parameter_gradient_complex(mol%nat, mol%id, &
            & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bgsp, &
            & diagx, grad_onsfx)
         CALL onsite_fx_symv_parameter_gradient_complex(mol%nat, mol%id, &
            & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bgsps, &
            & diag2d, grad_onsfx)
      END DO

      vsh = 0.0_wp
      DO iat = 1, mol%nat
         izp = mol%id(iat)
         is = self%ish_at(iat)
         DO ish = 1, self%nsh_id(izp)
            iq = is + ish
            vsh(iq) = SUM(grad_onsfx(:, :, iat) &
                         & *cache%dgdq_onsfx(:, :, iq)) + SUM(grad_onsri(:, iat) &
                                                             & *cache%dgdq_onsri(:, iq))
         END DO
      END DO

   END SUBROUTINE get_KFock_kmesh_vsh

! **************************************************************************************************
!> \brief ...
!> \param self ...
!> \param mol ...
!> \param density ...
!> \param overlap ...
!> \param kfrac ...
!> \param weights ...
!> \param kernel ...
!> \param g_onsfx ...
!> \param g_onsri ...
!> \param fock ...
!> \param use_bvk_kernel ...
!> \param phase_forward_cached ...
!> \param phase_inverse_cached ...
!> \param transform_backend ...
!> \param nmesh ...
!> \param input_to_grid ...
!> \param phase_roots ...
!> \param twist_phase ...
!> \param transform ...
! **************************************************************************************************
   SUBROUTINE build_KFock_kmesh(self, mol, density, overlap, kfrac, weights, &
      & kernel, g_onsfx, g_onsri, fock, use_bvk_kernel, &
      & phase_forward_cached, phase_inverse_cached, transform_backend, nmesh, &
      & input_to_grid, phase_roots, twist_phase, transform)
      class(exchange_fock), INTENT(in) :: self
      TYPE(structure_type), INTENT(in) :: mol
      COMPLEX(wp), INTENT(in) :: density(:, :, :, :), overlap(:, :, :)
      REAL(wp), INTENT(in) :: kfrac(:, :), weights(:)
      TYPE(exchange_bvk_kernel), INTENT(in) :: kernel
      REAL(wp), INTENT(in) :: g_onsfx(:, :, :), g_onsri(:, :)
      COMPLEX(wp), INTENT(out) :: fock(:, :, :, :)
      LOGICAL, INTENT(in) :: use_bvk_kernel
      COMPLEX(wp), contiguous, OPTIONAL, TARGET, INTENT(in) :: phase_forward_cached(:, :)
      COMPLEX(wp), contiguous, OPTIONAL, TARGET, INTENT(in) :: phase_inverse_cached(:, :)
      INTEGER, INTENT(in), OPTIONAL :: transform_backend, nmesh(3), &
         & input_to_grid(:)
      COMPLEX(wp), INTENT(in), OPTIONAL :: phase_roots(:, :), twist_phase(:)
      class(exchange_transform_type), INTENT(inout), OPTIONAL :: transform

      INTEGER :: iao, backend, icell, ii, jj, ik, spin
      LOGICAL :: use_cached_phases, use_compact
      REAL(wp) :: angle, spin_factor
      COMPLEX(wp) :: phase, tmp
      COMPLEX(wp), ALLOCATABLE, TARGET :: amat(:, :, :), cmat(:, :, :)
      COMPLEX(wp), ALLOCATABLE :: diagP(:), diagSP(:), diagSPS(:)
      COMPLEX(wp), ALLOCATABLE :: gdiagP(:), gdiagSP(:), gdiagSPS(:)
      COMPLEX(wp), ALLOCATABLE, TARGET :: vmat(:, :, :)
      COMPLEX(wp), ALLOCATABLE, TARGET :: source_r(:, :, :), result_r(:, :, :)
      COMPLEX(wp), ALLOCATABLE :: work(:, :)
      COMPLEX(wp), ALLOCATABLE, TARGET :: phase_forward_local(:, :), phase_inverse_local(:, :)
      COMPLEX(wp), contiguous, POINTER :: phase_forward(:, :), phase_inverse(:, :)

      ALLOCATE (amat(self%nao, self%nao, SIZE(weights)), &
         & cmat(self%nao, self%nao, SIZE(weights)), &
         & vmat(self%nao, self%nao, SIZE(weights)), &
         & source_r(self%nao, self%nao, SIZE(kernel%reps, 2)), &
         & result_r(self%nao, self%nao, SIZE(kernel%reps, 2)), &
         & work(self%nao, self%nao), diagP(self%nao), diagSP(self%nao), &
         & diagSPS(self%nao), gdiagP(self%nao), gdiagSP(self%nao), &
         & gdiagSPS(self%nao))
      backend = bvk_transform_dense
      IF (PRESENT(transform_backend)) backend = transform_backend
      use_cached_phases = .FALSE.
      use_compact = is_compact_transform_backend(backend) .AND. PRESENT(nmesh) &
         & .AND. PRESENT(input_to_grid) .AND. PRESENT(phase_roots) &
         & .AND. PRESENT(twist_phase)
      IF (PRESENT(transform)) use_compact = .FALSE.
      IF (use_compact .OR. PRESENT(transform)) THEN
         ! Compact-grid backends do not materialize an Nk-by-Nk phase table.
         ! Associate harmless scalar placeholders because the dense arguments of
         ! the common kernel entry point are intentionally retained as oracle.
         ALLOCATE (phase_forward_local(1, 1), phase_inverse_local(1, 1), &
            & source=(0.0_wp, 0.0_wp))
         phase_forward => phase_forward_local
         phase_inverse => phase_inverse_local
      ELSE IF (PRESENT(phase_forward_cached) .AND. PRESENT(phase_inverse_cached)) THEN
         use_cached_phases = ALL(SHAPE(phase_forward_cached) == &
            & [SIZE(kernel%reps, 2), SIZE(weights)]) .AND. &
            & ALL(SHAPE(phase_inverse_cached) == &
            & [SIZE(weights), SIZE(kernel%reps, 2)])
      END IF
      IF (.NOT. use_compact .AND. use_cached_phases) THEN
         phase_forward => phase_forward_cached
         phase_inverse => phase_inverse_cached
      ELSE IF (.NOT. use_compact .AND. .NOT. PRESENT(transform)) THEN
         ALLOCATE (phase_forward_local(SIZE(kernel%reps, 2), SIZE(weights)), &
            & phase_inverse_local(SIZE(weights), SIZE(kernel%reps, 2)))
         phase_forward => phase_forward_local
         phase_inverse => phase_inverse_local

         ! No ordering of the k points is assumed, and a common
         ! Monkhorst-Pack twist is retained explicitly in every phase.
         DO ik = 1, SIZE(weights)
            DO icell = 1, SIZE(kernel%reps, 2)
               angle = 2.0_wp*pi*DOT_PRODUCT(kfrac(:, ik), &
                  & REAL(kernel%reps(:, icell), wp))
               phase = EXP(CMPLX(0.0_wp, angle, wp))
               phase_forward(icell, ik) = phase
               phase_inverse(ik, icell) = weights(ik)*CONJG(phase)
            END DO
         END DO
      END IF

      spin_factor = 0.5_wp
      IF (SIZE(density, 3) > 1) spin_factor = 1.0_wp
      fock = (0.0_wp, 0.0_wp)

      DO spin = 1, SIZE(density, 3)
         diagP = (0.0_wp, 0.0_wp)
         diagSP = (0.0_wp, 0.0_wp)
         diagSPS = (0.0_wp, 0.0_wp)
         DO ik = 1, SIZE(weights)
            amat(:, :, ik) = MATMUL(overlap(:, :, ik), density(:, :, spin, ik))
            cmat(:, :, ik) = 0.5_wp*MATMUL(amat(:, :, ik), overlap(:, :, ik))
            DO iao = 1, self%nao
               diagP(iao) = diagP(iao) + weights(ik)*density(iao, iao, spin, ik)
               diagSP(iao) = diagSP(iao) + weights(ik)*amat(iao, iao, ik)
               diagSPS(iao) = diagSPS(iao) + 2.0_wp*weights(ik)*cmat(iao, iao, ik)
            END DO
         END DO
         CALL onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
            & self%nao_sh, self%ish_at, self%iao_sh, g_onsfx, diagP, gdiagP)
         CALL onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
            & self%nao_sh, self%ish_at, self%iao_sh, g_onsfx, diagSP, gdiagSP)
         CALL onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
            & self%nao_sh, self%ish_at, self%iao_sh, g_onsfx, diagSPS, gdiagSPS)

         CALL apply_bvk_Kkernel_inplace(self, mol, kernel, phase_inverse, phase_forward, &
            & amat, source_r, result_r, &
            & g_onsfx, g_onsri, 1.0_wp, -4.0_wp, 0.5_wp, -2.0_wp, &
            & use_bvk_kernel, .TRUE., backend, nmesh, input_to_grid, &
            & phase_roots, twist_phase, &
            & transform=transform)
         CALL apply_bvk_Kkernel_inplace(self, mol, kernel, phase_inverse, phase_forward, &
            & cmat, source_r, result_r, &
            & g_onsfx, g_onsri, 1.0_wp, 0.0_wp, 0.5_wp, -2.0_wp, &
            & use_bvk_kernel, .FALSE., backend, nmesh, input_to_grid, &
            & phase_roots, twist_phase, &
            & transform=transform)

         ! A fixed-spin rank-four section is non-contiguous for UKS.  Copy it
         ! explicitly into the third reusable mesh buffer instead of relying on
         ! a compiler-created pack temporary at the contiguous kernel boundary.
         DO ik = 1, SIZE(weights)
            vmat(:, :, ik) = density(:, :, spin, ik)
         END DO
         CALL apply_bvk_Kkernel_inplace(self, mol, kernel, phase_inverse, phase_forward, &
            & vmat, source_r, result_r, g_onsfx, g_onsri, 1.0_wp, 0.0_wp, &
            & 0.5_wp, -2.0_wp, use_bvk_kernel, .FALSE., backend, nmesh, &
            & input_to_grid, phase_roots, twist_phase, &
            & transform=transform)

         DO ik = 1, SIZE(weights)
            work = amat(:, :, ik)
            DO iao = 1, self%nao
               work(:, iao) = work(:, iao) &
                  & + 0.25_wp*gdiagP(iao)*overlap(:, iao, ik)
            END DO
            work = work + 0.5_wp*MATMUL(overlap(:, :, ik), vmat(:, :, ik))
            fock(:, :, spin, ik) = cmat(:, :, ik) &
               & + MATMUL(work, overlap(:, :, ik))
            DO iao = 1, self%nao
               fock(:, iao, spin, ik) = fock(:, iao, spin, ik) &
                  & + 0.5_wp*gdiagSP(iao)*overlap(:, iao, ik)
               fock(iao, iao, spin, ik) = fock(iao, iao, spin, ik) &
                  & + 0.25_wp*gdiagSPS(iao)
            END DO

            DO ii = 1, self%nao
               fock(ii, ii, spin, ik) = CMPLX(-0.5_wp*spin_factor &
                                             & *REAL(fock(ii, ii, spin, ik), wp), 0.0_wp, wp)
               DO jj = 1, ii - 1
                  tmp = -0.25_wp*spin_factor*(fock(jj, ii, spin, ik) &
                     & + CONJG(fock(ii, jj, spin, ik)))
                  fock(jj, ii, spin, ik) = tmp
                  fock(ii, jj, spin, ik) = CONJG(tmp)
               END DO
            END DO
         END DO
      END DO

   END SUBROUTINE build_KFock_kmesh

! **************************************************************************************************
!> \brief ...
!> \param self ...
!> \param mol ...
!> \param kernel ...
!> \param source_r ...
!> \param result_r ...
!> \param g_onsfx ...
!> \param g_onsri ...
!> \param mulliken_scale ...
!> \param bocorr_scale ...
!> \param onsite_fx_scale ...
!> \param onsite_ri_scale ...
!> \param use_bvk_kernel ...
!> \param adjoint_src ...
! **************************************************************************************************
   SUBROUTINE apply_bvk_Kkernel_realspace(self, mol, kernel, source_r, result_r, &
      & g_onsfx, g_onsri, mulliken_scale, bocorr_scale, onsite_fx_scale, &
      & onsite_ri_scale, use_bvk_kernel, adjoint_src)
      class(exchange_fock), INTENT(in) :: self
      TYPE(structure_type), INTENT(in) :: mol
      TYPE(exchange_bvk_kernel), INTENT(in) :: kernel
      COMPLEX(wp), INTENT(in) :: source_r(:, :, :)
      COMPLEX(wp), INTENT(out) :: result_r(:, :, :)
      REAL(wp), INTENT(in) :: g_onsfx(:, :, :), g_onsri(:, :)
      REAL(wp), INTENT(in) :: mulliken_scale, bocorr_scale
      REAL(wp), INTENT(in) :: onsite_fx_scale, onsite_ri_scale
      LOGICAL, INTENT(in) :: use_bvk_kernel, adjoint_src

      INTEGER :: icell, origin

      result_r = (0.0_wp, 0.0_wp)
      origin = 0
      DO icell = 1, SIZE(kernel%reps, 2)
         IF (ALL(kernel%reps(:, icell) == 0)) origin = icell
         IF (use_bvk_kernel) THEN
            IF (mulliken_scale /= 0.0_wp) THEN
               CALL shell_hadamard_add_complex(self%nsh, self%nao_sh, &
                   & self%iao_sh, kernel%g_mulliken_r(:, :, icell), &
                   & source_r(:, :, icell), mulliken_scale, &
                   & result_r(:, :, icell))
            END IF
            IF (bocorr_scale /= 0.0_wp) THEN
               CALL atom_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
                   & self%nao_sh, self%ish_at, self%iao_sh, &
                   & kernel%g_bocorr_r(:, :, icell), source_r(:, :, icell), &
                   & bocorr_scale, result_r(:, :, icell))
            END IF
         END IF
      END DO

      IF (origin > 0) THEN
         IF (onsite_fx_scale /= 0.0_wp) THEN
            CALL onsite_fx_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
                & self%nao_sh, self%ish_at, self%iao_sh, g_onsfx, &
                & source_r(:, :, origin), onsite_fx_scale, &
                & result_r(:, :, origin), adjoint_src=adjoint_src)
         END IF
         IF (onsite_ri_scale /= 0.0_wp) THEN
            CALL onsite_ri_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
                & self%nao_sh, self%ish_at, self%iao_sh, g_onsri, &
                & source_r(:, :, origin), onsite_ri_scale, &
                & result_r(:, :, origin))
         END IF
      END IF
   END SUBROUTINE apply_bvk_Kkernel_realspace

! **************************************************************************************************
!> \brief ...
!> \param self ...
!> \param mol ...
!> \param cache ...
!> \param kernel ...
!> \param phase_forward ...
!> \param weights ...
!> \param amat_r ...
!> \param cmat_r ...
!> \param vmat_r ...
!> \param gdiagP ...
!> \param gdiagSP ...
!> \param gdiagSPS ...
!> \param vsh ...
!> \param energy ...
! **************************************************************************************************
   SUBROUTINE get_KFock_stream_apply(self, mol, cache, kernel, phase_forward, &
      & weights, amat_r, cmat_r, vmat_r, gdiagP, gdiagSP, gdiagSPS, vsh, energy)
      class(exchange_fock), INTENT(in) :: self
      TYPE(structure_type), INTENT(in) :: mol
      TYPE(exchange_cache), INTENT(in) :: cache
      TYPE(exchange_bvk_kernel), INTENT(in) :: kernel
      COMPLEX(wp), INTENT(in) :: phase_forward(:, :)
      REAL(wp), INTENT(in) :: weights(:)
      COMPLEX(wp), INTENT(inout) :: amat_r(:, :, :, :), cmat_r(:, :, :, :), &
         & vmat_r(:, :, :, :)
      COMPLEX(wp), INTENT(out) :: gdiagP(:, :), gdiagSP(:, :), gdiagSPS(:, :)
      REAL(wp), INTENT(out) :: vsh(:), energy

      INTEGER :: iao, iat, icell, iq, is, ish, izp, jcell, origin, spin
      REAL(wp) :: spin_factor
      COMPLEX(wp) :: coeff, energy_raw
      REAL(wp), ALLOCATABLE :: grad_onsfx(:, :, :), grad_onsri(:, :)
      COMPLEX(wp), ALLOCATABLE :: bgp(:), bgsp(:), bgsps(:), bmat(:, :), &
         & bt0(:, :), bu0(:, :), bv0(:, :), diag2d(:), diagp(:), diagx(:), &
         & result_r(:, :, :), result2_r(:, :, :)

      origin = 0
      DO icell = 1, SIZE(kernel%reps, 2)
         IF (ALL(kernel%reps(:, icell) == 0)) origin = icell
      END DO

      ALLOCATE (grad_onsfx(self%maxsh, self%maxsh, mol%nat), &
         & grad_onsri(self%maxsh, mol%nat), source=0.0_wp)
      ALLOCATE (bgp(self%nao), bgsp(self%nao), bgsps(self%nao), &
         & bmat(self%nao, self%nao), bt0(self%nao, self%nao), &
         & bu0(self%nao, self%nao), bv0(self%nao, self%nao), &
         & diag2d(self%nao), diagp(self%nao), diagx(self%nao), &
         & result_r(self%nao, self%nao, SIZE(kernel%reps, 2)), &
         & result2_r(self%nao, self%nao, SIZE(kernel%reps, 2)), &
         & source=(0.0_wp, 0.0_wp))

      energy = 0.0_wp
      vsh = 0.0_wp
      spin_factor = 0.5_wp
      IF (SIZE(amat_r, 4) > 1) spin_factor = 1.0_wp

      DO spin = 1, SIZE(amat_r, 4)
         DO iao = 1, self%nao
            diagp(iao) = vmat_r(iao, iao, origin, spin)
            diagx(iao) = amat_r(iao, iao, origin, spin)
            diag2d(iao) = 2.0_wp*cmat_r(iao, iao, origin, spin)
         END DO
         CALL onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
            & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, &
            & diagp, gdiagP(:, spin))
         CALL onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
            & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, &
            & diagx, gdiagSP(:, spin))
         CALL onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
            & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, &
            & diag2d, gdiagSPS(:, spin))

         bt0 = -0.25_wp*spin_factor*CONJG( &
            & TRANSPOSE(amat_r(:, :, origin, spin)))
         bu0 = -0.25_wp*spin_factor*vmat_r(:, :, origin, spin)
         bv0 = -0.25_wp*spin_factor*cmat_r(:, :, origin, spin)
         DO iao = 1, self%nao
            bgp(iao) = -0.125_wp*spin_factor &
                      & *cmat_r(iao, iao, origin, spin)
            bgsp(iao) = -0.125_wp*spin_factor &
                       & *amat_r(iao, iao, origin, spin)
            bgsps(iao) = -0.0625_wp*spin_factor &
                        & *vmat_r(iao, iao, origin, spin)
         END DO

         CALL onsite_fx_parameter_gradient_complex(mol%nat, mol%id, &
            & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bt0, &
            & amat_r(:, :, origin, spin), 0.5_wp, grad_onsfx, &
            & adjoint_src=.TRUE.)
         CALL onsite_fx_parameter_gradient_complex(mol%nat, mol%id, &
            & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bu0, &
            & cmat_r(:, :, origin, spin), 0.5_wp, grad_onsfx)
         CALL onsite_fx_parameter_gradient_complex(mol%nat, mol%id, &
            & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bv0, &
            & vmat_r(:, :, origin, spin), 0.5_wp, grad_onsfx)
         CALL onsite_ri_parameter_gradient_complex(mol%nat, mol%id, &
            & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bt0, &
            & amat_r(:, :, origin, spin), -2.0_wp, grad_onsri)
         CALL onsite_ri_parameter_gradient_complex(mol%nat, mol%id, &
            & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bu0, &
            & cmat_r(:, :, origin, spin), -2.0_wp, grad_onsri)
         CALL onsite_ri_parameter_gradient_complex(mol%nat, mol%id, &
            & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bv0, &
            & vmat_r(:, :, origin, spin), -2.0_wp, grad_onsri)
         CALL onsite_fx_symv_parameter_gradient_complex(mol%nat, mol%id, &
            & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bgp, &
            & diagp, grad_onsfx)
         CALL onsite_fx_symv_parameter_gradient_complex(mol%nat, mol%id, &
            & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bgsp, &
            & diagx, grad_onsfx)
         CALL onsite_fx_symv_parameter_gradient_complex(mol%nat, mol%id, &
            & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bgsps, &
            & diag2d, grad_onsfx)

         energy_raw = (0.0_wp, 0.0_wp)
         DO iao = 1, self%nao
            energy_raw = energy_raw &
               & + 0.5_wp*cmat_r(iao, iao, origin, spin)*gdiagP(iao, spin) &
               & + 0.5_wp*CONJG(amat_r(iao, iao, origin, spin)) &
               & *gdiagSP(iao, spin) &
               & + 0.25_wp*vmat_r(iao, iao, origin, spin)*gdiagSPS(iao, spin)
         END DO

         CALL apply_bvk_Kkernel_realspace(self, mol, kernel, &
            & amat_r(:, :, :, spin), result_r, cache%g_onsfx, cache%g_onsri, &
            & 1.0_wp, -4.0_wp, 0.5_wp, -2.0_wp, .TRUE., .TRUE.)
         DO icell = 1, SIZE(kernel%reps, 2)
            bmat = (0.0_wp, 0.0_wp)
            DO jcell = 1, SIZE(kernel%reps, 2)
               coeff = SUM(weights*CONJG(phase_forward(icell, :) &
                                        & *phase_forward(jcell, :)))
               bmat = bmat + coeff*CONJG( &
                  & TRANSPOSE(amat_r(:, :, jcell, spin)))
            END DO
            energy_raw = energy_raw + SUM(CONJG(bmat)*result_r(:, :, icell))
         END DO
         amat_r(:, :, :, spin) = result_r

         CALL apply_bvk_Kkernel_realspace(self, mol, kernel, &
            & cmat_r(:, :, :, spin), result_r, cache%g_onsfx, cache%g_onsri, &
            & 1.0_wp, 0.0_wp, 0.5_wp, -2.0_wp, .TRUE., .FALSE.)
         energy_raw = energy_raw + SUM(CONJG(vmat_r(:, :, :, spin))*result_r)

         CALL apply_bvk_Kkernel_realspace(self, mol, kernel, &
            & vmat_r(:, :, :, spin), result2_r, cache%g_onsfx, &
            & cache%g_onsri, 1.0_wp, 0.0_wp, 0.5_wp, -2.0_wp, .TRUE., .FALSE.)
         energy_raw = energy_raw + SUM(CONJG(cmat_r(:, :, :, spin)) &
                                      & *result2_r)
         vmat_r(:, :, :, spin) = result2_r
         cmat_r(:, :, :, spin) = result_r
         energy = energy - 0.25_wp*spin_factor*REAL(energy_raw, wp)
      END DO

      DO iat = 1, mol%nat
         izp = mol%id(iat)
         is = self%ish_at(iat)
         DO ish = 1, self%nsh_id(izp)
            iq = is + ish
            vsh(iq) = SUM(grad_onsfx(:, :, iat) &
                         & *cache%dgdq_onsfx(:, :, iq)) + SUM(grad_onsri(:, iat) &
                                                             & *cache%dgdq_onsri(:, iq))
         END DO
      END DO
   END SUBROUTINE get_KFock_stream_apply

! **************************************************************************************************
!> \brief ...
!> \param self ...
!> \param mol ...
!> \param kernel ...
!> \param icell ...
!> \param source ...
!> \param RESULT ...
!> \param g_onsfx ...
!> \param g_onsri ...
!> \param mulliken_scale ...
!> \param bocorr_scale ...
!> \param onsite_fx_scale ...
!> \param onsite_ri_scale ...
!> \param use_bvk_kernel ...
!> \param adjoint_src ...
! **************************************************************************************************
   SUBROUTINE apply_bvk_Kkernel_image(self, mol, kernel, icell, source, RESULT, &
      & g_onsfx, g_onsri, mulliken_scale, bocorr_scale, onsite_fx_scale, &
      & onsite_ri_scale, use_bvk_kernel, adjoint_src)
      class(exchange_fock), INTENT(in) :: self
      TYPE(structure_type), INTENT(in) :: mol
      TYPE(exchange_bvk_kernel), INTENT(in) :: kernel
      INTEGER, INTENT(in) :: icell
      COMPLEX(wp), INTENT(in) :: source(:, :)
      COMPLEX(wp), INTENT(out) :: result(:, :)
      REAL(wp), INTENT(in) :: g_onsfx(:, :, :), g_onsri(:, :)
      REAL(wp), INTENT(in) :: mulliken_scale, bocorr_scale
      REAL(wp), INTENT(in) :: onsite_fx_scale, onsite_ri_scale
      LOGICAL, INTENT(in) :: use_bvk_kernel, adjoint_src

      RESULT = (0.0_wp, 0.0_wp)
      IF (use_bvk_kernel) THEN
         IF (mulliken_scale /= 0.0_wp) THEN
            CALL shell_hadamard_add_complex(self%nsh, self%nao_sh, &
                & self%iao_sh, kernel%g_mulliken_r(:, :, icell), source, &
                & mulliken_scale, RESULT)
         END IF
         IF (bocorr_scale /= 0.0_wp) THEN
            CALL atom_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
                & self%nao_sh, self%ish_at, self%iao_sh, &
                & kernel%g_bocorr_r(:, :, icell), source, bocorr_scale, RESULT)
         END IF
      END IF
      IF (ALL(kernel%reps(:, icell) == 0)) THEN
         IF (onsite_fx_scale /= 0.0_wp) THEN
            CALL onsite_fx_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
                & self%nao_sh, self%ish_at, self%iao_sh, g_onsfx, source, &
                & onsite_fx_scale, RESULT, adjoint_src=adjoint_src)
         END IF
         IF (onsite_ri_scale /= 0.0_wp) THEN
            CALL onsite_ri_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
                & self%nao_sh, self%ish_at, self%iao_sh, g_onsri, source, &
                & onsite_ri_scale, RESULT)
         END IF
      END IF
   END SUBROUTINE apply_bvk_Kkernel_image

! **************************************************************************************************
!> \brief ...
!> \param self ...
!> \param mol ...
!> \param cache ...
!> \param kernel ...
!> \param image_index ...
!> \param inverse_coeff ...
!> \param amat_r ...
!> \param amat_partner_r ...
!> \param cmat_r ...
!> \param vmat_r ...
!> \param gdiagP ...
!> \param gdiagSP ...
!> \param gdiagSPS ...
!> \param vsh ...
!> \param energy ...
! **************************************************************************************************
   SUBROUTINE get_KFock_batch_apply(self, mol, cache, kernel, image_index, &
      & inverse_coeff, amat_r, amat_partner_r, cmat_r, vmat_r, gdiagP, &
      & gdiagSP, gdiagSPS, vsh, energy)
      class(exchange_fock), INTENT(in) :: self
      TYPE(structure_type), INTENT(in) :: mol
      TYPE(exchange_cache), INTENT(in) :: cache
      TYPE(exchange_bvk_kernel), INTENT(in) :: kernel
      INTEGER, INTENT(in) :: image_index(:)
      COMPLEX(wp), INTENT(in) :: inverse_coeff(:)
      COMPLEX(wp), INTENT(inout) :: amat_r(:, :, :, :), &
         & amat_partner_r(:, :, :, :), cmat_r(:, :, :, :), &
         & vmat_r(:, :, :, :)
      COMPLEX(wp), INTENT(inout) :: gdiagP(:, :), gdiagSP(:, :), &
         & gdiagSPS(:, :)
      REAL(wp), INTENT(out) :: vsh(:), energy

      INTEGER :: iao, iat, ibatch, icell, iq, is, ish, izp, origin_local, spin
      REAL(wp) :: spin_factor
      COMPLEX(wp) :: energy_raw
      REAL(wp), ALLOCATABLE :: grad_onsfx(:, :, :), grad_onsri(:, :)
      COMPLEX(wp), ALLOCATABLE :: bgp(:), bgsp(:), bgsps(:), bmat(:, :), &
         & bt0(:, :), bu0(:, :), bv0(:, :), diag2d(:), diagp(:), diagx(:), &
         & work(:, :)

      origin_local = 0
      DO ibatch = 1, SIZE(image_index)
         IF (ALL(kernel%reps(:, image_index(ibatch)) == 0)) THEN
            origin_local = ibatch
         END IF
      END DO

      ALLOCATE (grad_onsfx(self%maxsh, self%maxsh, mol%nat), &
         & grad_onsri(self%maxsh, mol%nat), source=0.0_wp)
      ALLOCATE (bgp(self%nao), bgsp(self%nao), bgsps(self%nao), &
         & bmat(self%nao, self%nao), bt0(self%nao, self%nao), &
         & bu0(self%nao, self%nao), bv0(self%nao, self%nao), &
         & diag2d(self%nao), diagp(self%nao), diagx(self%nao), &
         & work(self%nao, self%nao), source=(0.0_wp, 0.0_wp))

      energy = 0.0_wp
      vsh = 0.0_wp
      spin_factor = 0.5_wp
      IF (SIZE(amat_r, 4) > 1) spin_factor = 1.0_wp

      DO spin = 1, SIZE(amat_r, 4)
         energy_raw = (0.0_wp, 0.0_wp)
         IF (origin_local > 0) THEN
            DO iao = 1, self%nao
               diagp(iao) = vmat_r(iao, iao, origin_local, spin)
               diagx(iao) = amat_r(iao, iao, origin_local, spin)
               diag2d(iao) = 2.0_wp*cmat_r(iao, iao, origin_local, spin)
            END DO
            CALL onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
               & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, &
               & diagp, gdiagP(:, spin))
            CALL onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
               & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, &
               & diagx, gdiagSP(:, spin))
            CALL onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
               & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, &
               & diag2d, gdiagSPS(:, spin))

            bt0 = -0.25_wp*spin_factor*CONJG( &
               & TRANSPOSE(amat_r(:, :, origin_local, spin)))
            bu0 = -0.25_wp*spin_factor*vmat_r(:, :, origin_local, spin)
            bv0 = -0.25_wp*spin_factor*cmat_r(:, :, origin_local, spin)
            DO iao = 1, self%nao
               bgp(iao) = -0.125_wp*spin_factor &
                         & *cmat_r(iao, iao, origin_local, spin)
               bgsp(iao) = -0.125_wp*spin_factor &
                          & *amat_r(iao, iao, origin_local, spin)
               bgsps(iao) = -0.0625_wp*spin_factor &
                           & *vmat_r(iao, iao, origin_local, spin)
            END DO

            CALL onsite_fx_parameter_gradient_complex(mol%nat, mol%id, &
               & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bt0, &
               & amat_r(:, :, origin_local, spin), 0.5_wp, grad_onsfx, &
               & adjoint_src=.TRUE.)
            CALL onsite_fx_parameter_gradient_complex(mol%nat, mol%id, &
               & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bu0, &
               & cmat_r(:, :, origin_local, spin), 0.5_wp, grad_onsfx)
            CALL onsite_fx_parameter_gradient_complex(mol%nat, mol%id, &
               & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bv0, &
               & vmat_r(:, :, origin_local, spin), 0.5_wp, grad_onsfx)
            CALL onsite_ri_parameter_gradient_complex(mol%nat, mol%id, &
               & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bt0, &
               & amat_r(:, :, origin_local, spin), -2.0_wp, grad_onsri)
            CALL onsite_ri_parameter_gradient_complex(mol%nat, mol%id, &
               & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bu0, &
               & cmat_r(:, :, origin_local, spin), -2.0_wp, grad_onsri)
            CALL onsite_ri_parameter_gradient_complex(mol%nat, mol%id, &
               & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bv0, &
               & vmat_r(:, :, origin_local, spin), -2.0_wp, grad_onsri)
            CALL onsite_fx_symv_parameter_gradient_complex(mol%nat, mol%id, &
               & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bgp, &
               & diagp, grad_onsfx)
            CALL onsite_fx_symv_parameter_gradient_complex(mol%nat, mol%id, &
               & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bgsp, &
               & diagx, grad_onsfx)
            CALL onsite_fx_symv_parameter_gradient_complex(mol%nat, mol%id, &
               & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bgsps, &
               & diag2d, grad_onsfx)

            DO iao = 1, self%nao
               energy_raw = energy_raw &
                  & + 0.5_wp*cmat_r(iao, iao, origin_local, spin) &
                  & *gdiagP(iao, spin) &
                  & + 0.5_wp*CONJG(amat_r(iao, iao, origin_local, spin)) &
                  & *gdiagSP(iao, spin) &
                  & + 0.25_wp*vmat_r(iao, iao, origin_local, spin) &
                  & *gdiagSPS(iao, spin)
            END DO
         END IF

         DO ibatch = 1, SIZE(image_index)
            icell = image_index(ibatch)
            CALL apply_bvk_Kkernel_image(self, mol, kernel, icell, &
               & amat_r(:, :, ibatch, spin), work, cache%g_onsfx, &
               & cache%g_onsri, 1.0_wp, -4.0_wp, 0.5_wp, -2.0_wp, &
               & .TRUE., .TRUE.)
            bmat = inverse_coeff(ibatch)*CONJG( &
               & TRANSPOSE(amat_partner_r(:, :, ibatch, spin)))
            energy_raw = energy_raw + SUM(CONJG(bmat)*work)
            amat_r(:, :, ibatch, spin) = work

            CALL apply_bvk_Kkernel_image(self, mol, kernel, icell, &
               & cmat_r(:, :, ibatch, spin), work, cache%g_onsfx, &
               & cache%g_onsri, 1.0_wp, 0.0_wp, 0.5_wp, -2.0_wp, &
               & .TRUE., .FALSE.)
            energy_raw = energy_raw &
               & + SUM(CONJG(vmat_r(:, :, ibatch, spin))*work)
            amat_partner_r(:, :, ibatch, spin) = work

            CALL apply_bvk_Kkernel_image(self, mol, kernel, icell, &
               & vmat_r(:, :, ibatch, spin), work, cache%g_onsfx, &
               & cache%g_onsri, 1.0_wp, 0.0_wp, 0.5_wp, -2.0_wp, &
               & .TRUE., .FALSE.)
            energy_raw = energy_raw &
               & + SUM(CONJG(cmat_r(:, :, ibatch, spin))*work)
            cmat_r(:, :, ibatch, spin) = &
               & amat_partner_r(:, :, ibatch, spin)
            vmat_r(:, :, ibatch, spin) = work
         END DO
         energy = energy - 0.25_wp*spin_factor*REAL(energy_raw, wp)
      END DO

      DO iat = 1, mol%nat
         izp = mol%id(iat)
         is = self%ish_at(iat)
         DO ish = 1, self%nsh_id(izp)
            iq = is + ish
            vsh(iq) = SUM(grad_onsfx(:, :, iat) &
                         & *cache%dgdq_onsfx(:, :, iq)) + SUM(grad_onsri(:, iat) &
                                                             & *cache%dgdq_onsri(:, iq))
         END DO
      END DO
   END SUBROUTINE get_KFock_batch_apply

! **************************************************************************************************
!> \brief ...
!> \param self ...
!> \param mol ...
!> \param kernel ...
!> \param icell ...
!> \param source ...
!> \param result_adjoint ...
!> \param g_onsfx ...
!> \param g_onsri ...
!> \param mulliken_scale ...
!> \param bocorr_scale ...
!> \param onsite_fx_scale ...
!> \param onsite_ri_scale ...
!> \param use_bvk_kernel ...
!> \param adjoint_src ...
!> \param source_adjoint ...
!> \param mulliken_grad ...
!> \param bocorr_grad ...
! **************************************************************************************************
   SUBROUTINE reverse_bvk_Kkernel_image(self, mol, kernel, icell, source, &
      & result_adjoint, g_onsfx, g_onsri, mulliken_scale, bocorr_scale, &
      & onsite_fx_scale, onsite_ri_scale, use_bvk_kernel, adjoint_src, &
      & source_adjoint, mulliken_grad, bocorr_grad)
      class(exchange_fock), INTENT(in) :: self
      TYPE(structure_type), INTENT(in) :: mol
      TYPE(exchange_bvk_kernel), INTENT(in) :: kernel
      INTEGER, INTENT(in) :: icell
      COMPLEX(wp), INTENT(in) :: source(:, :), result_adjoint(:, :)
      REAL(wp), INTENT(in) :: g_onsfx(:, :, :), g_onsri(:, :)
      REAL(wp), INTENT(in) :: mulliken_scale, bocorr_scale
      REAL(wp), INTENT(in) :: onsite_fx_scale, onsite_ri_scale
      LOGICAL, INTENT(in) :: use_bvk_kernel, adjoint_src
      COMPLEX(wp), INTENT(out) :: source_adjoint(:, :)
      REAL(wp), INTENT(inout) :: mulliken_grad(:, :), bocorr_grad(:, :)

      source_adjoint = (0.0_wp, 0.0_wp)
      IF (use_bvk_kernel) THEN
         IF (mulliken_scale /= 0.0_wp) THEN
            CALL shell_parameter_gradient_complex(self%nsh, self%nao_sh, &
               & self%iao_sh, result_adjoint, source, mulliken_scale, &
               & mulliken_grad)
            CALL shell_hadamard_add_complex(self%nsh, self%nao_sh, &
               & self%iao_sh, kernel%g_mulliken_r(:, :, icell), &
               & result_adjoint, mulliken_scale, source_adjoint)
         END IF
         IF (bocorr_scale /= 0.0_wp) THEN
            CALL atom_parameter_gradient_complex(mol%nat, mol%id, &
               & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, &
               & result_adjoint, source, bocorr_scale, bocorr_grad)
            CALL atom_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
               & self%nao_sh, self%ish_at, self%iao_sh, &
               & kernel%g_bocorr_r(:, :, icell), result_adjoint, &
               & bocorr_scale, source_adjoint)
         END IF
      END IF
      IF (ALL(kernel%reps(:, icell) == 0)) THEN
         IF (onsite_fx_scale /= 0.0_wp) THEN
            CALL onsite_fx_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
                & self%nao_sh, self%ish_at, self%iao_sh, g_onsfx, &
                & result_adjoint, onsite_fx_scale, source_adjoint, &
                & adjoint_src=adjoint_src)
         END IF
         IF (onsite_ri_scale /= 0.0_wp) THEN
            CALL onsite_ri_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
                & self%nao_sh, self%ish_at, self%iao_sh, g_onsri, &
                & result_adjoint, onsite_ri_scale, source_adjoint)
         END IF
      END IF
   END SUBROUTINE reverse_bvk_Kkernel_image

! **************************************************************************************************
!> \brief ...
!> \param self ...
!> \param mol ...
!> \param cache ...
!> \param kernel ...
!> \param image_index ...
!> \param xmat_r ...
!> \param dmat_r ...
!> \param pmat_r ...
!> \param bmat_r ...
!> \param amat_r ...
!> \param vseed_r ...
!> \param gdiag_p ...
!> \param gdiag_sp ...
!> \param gdiag_sps ...
!> \param bdiag_sp ...
!> \param bdiag_sps ...
!> \param gradient ...
!> \param sigma ...
! **************************************************************************************************
   SUBROUTINE get_KGrad_batch_apply(self, mol, cache, kernel, image_index, &
      & xmat_r, dmat_r, pmat_r, bmat_r, amat_r, vseed_r, gdiag_p, gdiag_sp, &
      & gdiag_sps, bdiag_sp, bdiag_sps, gradient, sigma)
      class(exchange_fock), INTENT(in) :: self
      TYPE(structure_type), INTENT(in) :: mol
      TYPE(exchange_cache), INTENT(in) :: cache
      TYPE(exchange_bvk_kernel), INTENT(in) :: kernel
      INTEGER, INTENT(in) :: image_index(:)
      COMPLEX(wp), INTENT(inout) :: xmat_r(:, :, :, :), &
         & dmat_r(:, :, :, :), pmat_r(:, :, :, :), bmat_r(:, :, :, :)
      COMPLEX(wp), INTENT(in) :: amat_r(:, :, :, :), vseed_r(:, :, :, :)
      COMPLEX(wp), INTENT(inout) :: gdiag_p(:, :), gdiag_sp(:, :), &
         & gdiag_sps(:, :), bdiag_sp(:, :), bdiag_sps(:, :)
      REAL(wp), INTENT(inout) :: gradient(:, :), sigma(:, :)

      INTEGER :: ibatch, icell, spin
      REAL(wp), ALLOCATABLE :: bocorr_grad(:, :, :), mulliken_grad(:, :, :)
      COMPLEX(wp), ALLOCATABLE :: source_adjoint(:, :), work(:, :)

      ALLOCATE (mulliken_grad(self%nao, self%nao, SIZE(image_index)), &
         & bocorr_grad(self%nao, self%nao, SIZE(image_index)), source=0.0_wp)
      ALLOCATE (source_adjoint(self%nao, self%nao), &
         & work(self%nao, self%nao), source=(0.0_wp, 0.0_wp))

      DO spin = 1, SIZE(xmat_r, 4)
         DO ibatch = 1, SIZE(image_index)
            icell = image_index(ibatch)

            CALL apply_bvk_Kkernel_image(self, mol, kernel, icell, &
               & xmat_r(:, :, ibatch, spin), work, cache%g_onsfx, &
               & cache%g_onsri, 1.0_wp, -4.0_wp, 0.5_wp, -2.0_wp, &
               & .TRUE., .TRUE.)
            CALL reverse_bvk_Kkernel_image(self, mol, kernel, icell, &
               & xmat_r(:, :, ibatch, spin), bmat_r(:, :, ibatch, spin), &
               & cache%g_onsfx, cache%g_onsri, 1.0_wp, -4.0_wp, 0.5_wp, &
               & -2.0_wp, .TRUE., .TRUE., source_adjoint, &
               & mulliken_grad(:, :, ibatch), bocorr_grad(:, :, ibatch))
            xmat_r(:, :, ibatch, spin) = work
            bmat_r(:, :, ibatch, spin) = source_adjoint

            CALL reverse_bvk_Kkernel_image(self, mol, kernel, icell, &
               & dmat_r(:, :, ibatch, spin), amat_r(:, :, ibatch, spin), &
               & cache%g_onsfx, cache%g_onsri, 1.0_wp, 0.0_wp, 0.5_wp, &
               & -2.0_wp, .TRUE., .FALSE., source_adjoint, &
               & mulliken_grad(:, :, ibatch), bocorr_grad(:, :, ibatch))
            dmat_r(:, :, ibatch, spin) = source_adjoint

            CALL apply_bvk_Kkernel_image(self, mol, kernel, icell, &
               & pmat_r(:, :, ibatch, spin), work, cache%g_onsfx, &
               & cache%g_onsri, 1.0_wp, 0.0_wp, 0.5_wp, -2.0_wp, &
               & .TRUE., .FALSE.)
            CALL reverse_bvk_Kkernel_image(self, mol, kernel, icell, &
               & pmat_r(:, :, ibatch, spin), vseed_r(:, :, ibatch, spin), &
               & cache%g_onsfx, cache%g_onsri, 1.0_wp, 0.0_wp, 0.5_wp, &
               & -2.0_wp, .TRUE., .FALSE., source_adjoint, &
               & mulliken_grad(:, :, ibatch), bocorr_grad(:, :, ibatch))
            pmat_r(:, :, ibatch, spin) = work
         END DO
      END DO

      ! The three BZ-averaged diagonal maps belong to the unique origin image.
      IF (ANY([(ALL(kernel%reps(:, image_index(ibatch)) == 0), &
         & ibatch=1, SIZE(image_index))])) THEN
         DO spin = 1, SIZE(xmat_r, 4)
            CALL onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
               & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, &
               & gdiag_p(:, spin), work(:, 1))
            gdiag_p(:, spin) = work(:, 1)
            CALL onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
               & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, &
               & gdiag_sp(:, spin), work(:, 1))
            gdiag_sp(:, spin) = work(:, 1)
            CALL onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
               & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, &
               & gdiag_sps(:, spin), work(:, 1))
            gdiag_sps(:, spin) = work(:, 1)
            CALL onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
               & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, &
               & bdiag_sp(:, spin), work(:, 1))
            bdiag_sp(:, spin) = work(:, 1)
            CALL onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
               & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, &
               & bdiag_sps(:, spin), work(:, 1))
            bdiag_sps(:, spin) = work(:, 1)
         END DO
      END IF

      CALL get_kernel_image_batch_derivs(self, mol, kernel, image_index, &
         & mulliken_grad, bocorr_grad, gradient, sigma)
   END SUBROUTINE get_KGrad_batch_apply

! **************************************************************************************************
!> \brief ...
!> \param self ...
!> \param mol ...
!> \param kernel ...
!> \param image_index ...
!> \param mulliken_grad ...
!> \param bocorr_grad ...
!> \param gradient ...
!> \param sigma ...
!> \note Reuse the public physical derivative with only the requested image labels.
!> \note No full response mesh or radial-kernel payload is copied. The provider
!> \note retains its Wigner--Seitz weights and their force/strain derivatives.
! **************************************************************************************************
   SUBROUTINE get_kernel_image_batch_derivs(self, mol, kernel, image_index, &
      & mulliken_grad, bocorr_grad, gradient, sigma)
      class(exchange_fock), INTENT(in) :: self
      TYPE(structure_type), INTENT(in) :: mol
      TYPE(exchange_bvk_kernel), INTENT(in) :: kernel
      INTEGER, INTENT(in) :: image_index(:)
      REAL(wp), INTENT(in) :: mulliken_grad(:, :, :), bocorr_grad(:, :, :)
      REAL(wp), INTENT(inout) :: gradient(:, :), sigma(:, :)
      TYPE(exchange_bvk_kernel) :: batch

      batch%nmesh = kernel%nmesh
      ! Materialize the vector-subscript section in an explicitly bounded array
      ! before passing its allocatable descriptor across the provider boundary.
      ALLOCATE (batch%reps(3, SIZE(image_index)))
      batch%reps = kernel%reps(:, image_index)
      CALL self%get_bvk_Kmatrix_derivs(mol, batch, mulliken_grad, &
         & bocorr_grad, gradient, sigma)
   END SUBROUTINE get_kernel_image_batch_derivs

! **************************************************************************************************
!> \brief ...
!> \param self ...
!> \param mol ...
!> \param kernel ...
!> \param phase_inverse ...
!> \param phase_forward ...
!> \param mesh_k ...
!> \param source_r ...
!> \param result_r ...
!> \param g_onsfx ...
!> \param g_onsri ...
!> \param mulliken_scale ...
!> \param bocorr_scale ...
!> \param onsite_fx_scale ...
!> \param onsite_ri_scale ...
!> \param use_bvk_kernel ...
!> \param adjoint_src ...
!> \param transform_backend ...
!> \param nmesh ...
!> \param input_to_grid ...
!> \param phase_roots ...
!> \param twist_phase ...
!> \param transform ...
! **************************************************************************************************
   SUBROUTINE apply_bvk_Kkernel_inplace(self, mol, kernel, phase_inverse, phase_forward, &
      & mesh_k, source_r, result_r, g_onsfx, g_onsri, mulliken_scale, &
      & bocorr_scale, onsite_fx_scale, onsite_ri_scale, use_bvk_kernel, &
      & adjoint_src, transform_backend, nmesh, input_to_grid, phase_roots, &
      & twist_phase, transform)
      class(exchange_fock), INTENT(in) :: self
      TYPE(structure_type), INTENT(in) :: mol
      TYPE(exchange_bvk_kernel), INTENT(in) :: kernel
      COMPLEX(wp), contiguous, INTENT(in) :: phase_inverse(:, :), phase_forward(:, :)
      COMPLEX(wp), contiguous, TARGET, INTENT(inout) :: mesh_k(:, :, :)
      COMPLEX(wp), contiguous, TARGET, INTENT(inout) :: source_r(:, :, :)
      COMPLEX(wp), contiguous, TARGET, INTENT(inout) :: result_r(:, :, :)
      REAL(wp), INTENT(in) :: g_onsfx(:, :, :), g_onsri(:, :)
      REAL(wp), INTENT(in) :: mulliken_scale, bocorr_scale
      REAL(wp), INTENT(in) :: onsite_fx_scale, onsite_ri_scale
      LOGICAL, INTENT(in) :: use_bvk_kernel, adjoint_src
      INTEGER, INTENT(in), OPTIONAL :: transform_backend, nmesh(3), &
         & input_to_grid(:)
      COMPLEX(wp), INTENT(in), OPTIONAL :: phase_roots(:, :), twist_phase(:)
      class(exchange_transform_type), INTENT(inout), OPTIONAL :: transform

      INTEGER :: backend, icell, npair, origin
      COMPLEX(wp), POINTER, contiguous :: mesh_k_2d(:, :), source_r_2d(:, :), &
         & result_r_2d(:, :)

      npair = self%nao*self%nao
      mesh_k_2d(1:npair, 1:SIZE(mesh_k, 3)) => mesh_k
      source_r_2d(1:npair, 1:SIZE(source_r, 3)) => source_r
      result_r_2d(1:npair, 1:SIZE(result_r, 3)) => result_r

      backend = bvk_transform_dense
      IF (PRESENT(transform_backend)) backend = transform_backend

      ! Treat all AO pairs as independent rows of one transform.  The dense path
      ! uses one complex matrix multiplication; compact paths consume the cached
      ! grid ordering, roots, and twist without materializing phase tables.
      IF (PRESENT(transform)) THEN
         CALL transform%apply(mesh_k_2d, source_r_2d, exchange_k_to_r)
      ELSE IF (is_compact_transform_backend(backend) .AND. PRESENT(nmesh) &
         & .AND. PRESENT(input_to_grid) .AND. PRESENT(phase_roots) &
         & .AND. PRESENT(twist_phase)) THEN
         SELECT CASE (backend)
         CASE (bvk_transform_mixed_radix_fft)
            CALL bvk_mixed_radix_k_to_r(mesh_k_2d, source_r_2d, nmesh, &
               & input_to_grid, kernel%reps, phase_roots, twist_phase, .TRUE.)
         CASE default
            CALL bvk_separable_k_to_r(mesh_k_2d, source_r_2d, nmesh, &
               & input_to_grid, kernel%reps, phase_roots, twist_phase, .TRUE.)
         END SELECT
      ELSE
         CALL matmul_complex_contiguous(mesh_k_2d, phase_inverse, source_r_2d)
      END IF
      result_r = (0.0_wp, 0.0_wp)

      origin = 0
      DO icell = 1, SIZE(kernel%reps, 2)
         IF (ALL(kernel%reps(:, icell) == 0)) origin = icell
      END DO

      IF (use_bvk_kernel) THEN
         DO icell = 1, SIZE(kernel%reps, 2)
            IF (mulliken_scale /= 0.0_wp) THEN
               CALL shell_hadamard_add_complex(self%nsh, self%nao_sh, &
                   & self%iao_sh, kernel%g_mulliken_r(:, :, icell), &
                   & source_r(:, :, icell), mulliken_scale, result_r(:, :, icell))
            END IF
            IF (bocorr_scale /= 0.0_wp) THEN
               CALL atom_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
                   & self%nao_sh, self%ish_at, self%iao_sh, &
                   & kernel%g_bocorr_r(:, :, icell), source_r(:, :, icell), &
                   & bocorr_scale, result_r(:, :, icell))
            END IF
         END DO
      END IF

      IF (origin > 0) THEN
         IF (onsite_fx_scale /= 0.0_wp) THEN
            CALL onsite_fx_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
                & self%nao_sh, self%ish_at, self%iao_sh, g_onsfx, &
                & source_r(:, :, origin), onsite_fx_scale, result_r(:, :, origin), &
                & adjoint_src=adjoint_src)
         END IF
         IF (onsite_ri_scale /= 0.0_wp) THEN
            CALL onsite_ri_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
                & self%nao_sh, self%ish_at, self%iao_sh, g_onsri, &
                & source_r(:, :, origin), onsite_ri_scale, result_r(:, :, origin))
         END IF
      END IF

      ! Reverse the same batched AO-pair map with the selected backend.
      IF (PRESENT(transform)) THEN
         CALL transform%apply(result_r_2d, mesh_k_2d, exchange_r_to_k)
      ELSE IF (is_compact_transform_backend(backend) .AND. PRESENT(nmesh) &
         & .AND. PRESENT(input_to_grid) .AND. PRESENT(phase_roots) &
         & .AND. PRESENT(twist_phase)) THEN
         SELECT CASE (backend)
         CASE (bvk_transform_mixed_radix_fft)
            CALL bvk_mixed_radix_r_to_k(result_r_2d, mesh_k_2d, nmesh, &
               & input_to_grid, kernel%reps, phase_roots, twist_phase)
         CASE default
            CALL bvk_separable_r_to_k(result_r_2d, mesh_k_2d, nmesh, &
               & input_to_grid, kernel%reps, phase_roots, twist_phase)
         END SELECT
      ELSE
         CALL matmul_complex_contiguous(result_r_2d, phase_forward, mesh_k_2d)
      END IF

   END SUBROUTINE apply_bvk_Kkernel_inplace

! **************************************************************************************************
!> \brief ...
!> \param self ...
!> \param mol ...
!> \param kernel ...
!> \param phase_inverse ...
!> \param phase_forward ...
!> \param source_k ...
!> \param g_onsfx ...
!> \param g_onsri ...
!> \param mulliken_scale ...
!> \param bocorr_scale ...
!> \param onsite_fx_scale ...
!> \param onsite_ri_scale ...
!> \param use_bvk_kernel ...
!> \param adjoint_src ...
!> \param result_k ...
!> \param transform_backend ...
!> \param nmesh ...
!> \param input_to_grid ...
!> \param phase_roots ...
!> \param twist_phase ...
!> \param transform ...
! **************************************************************************************************
   SUBROUTINE apply_bvk_Kkernel(self, mol, kernel, phase_inverse, phase_forward, &
      & source_k, g_onsfx, g_onsri, mulliken_scale, bocorr_scale, &
      & onsite_fx_scale, onsite_ri_scale, use_bvk_kernel, adjoint_src, result_k, &
      & transform_backend, nmesh, input_to_grid, phase_roots, twist_phase, transform)
      class(exchange_fock), INTENT(in) :: self
      TYPE(structure_type), INTENT(in) :: mol
      TYPE(exchange_bvk_kernel), INTENT(in) :: kernel
      COMPLEX(wp), INTENT(in) :: phase_inverse(:, :), phase_forward(:, :)
      COMPLEX(wp), contiguous, INTENT(in) :: source_k(:, :, :)
      REAL(wp), INTENT(in) :: g_onsfx(:, :, :), g_onsri(:, :)
      REAL(wp), INTENT(in) :: mulliken_scale, bocorr_scale
      REAL(wp), INTENT(in) :: onsite_fx_scale, onsite_ri_scale
      LOGICAL, INTENT(in) :: use_bvk_kernel, adjoint_src
      COMPLEX(wp), contiguous, TARGET, INTENT(out) :: result_k(:, :, :)
      INTEGER, INTENT(in), OPTIONAL :: transform_backend, nmesh(3), &
         & input_to_grid(:)
      COMPLEX(wp), INTENT(in), OPTIONAL :: phase_roots(:, :), twist_phase(:)
      class(exchange_transform_type), INTENT(inout), OPTIONAL :: transform

      COMPLEX(wp), ALLOCATABLE, TARGET :: result_r(:, :, :), source_r(:, :, :)

      ALLOCATE (source_r(self%nao, self%nao, SIZE(kernel%reps, 2)), &
         & result_r(self%nao, self%nao, SIZE(kernel%reps, 2)))
      result_k = source_k
      CALL apply_bvk_Kkernel_inplace(self, mol, kernel, phase_inverse, &
         & phase_forward, result_k, source_r, result_r, g_onsfx, g_onsri, &
         & mulliken_scale, bocorr_scale, onsite_fx_scale, onsite_ri_scale, &
         & use_bvk_kernel, adjoint_src, transform_backend, nmesh, input_to_grid, &
         & phase_roots, twist_phase, &
            & transform=transform)

   END SUBROUTINE apply_bvk_Kkernel

! **************************************************************************************************
!> \brief ...
!> \param self ...
!> \param mol ...
!> \param kernel ...
!> \param weights ...
!> \param phase_inverse ...
!> \param phase_forward ...
!> \param source_k ...
!> \param result_adjoint_k ...
!> \param g_onsfx ...
!> \param g_onsri ...
!> \param mulliken_scale ...
!> \param bocorr_scale ...
!> \param onsite_fx_scale ...
!> \param onsite_ri_scale ...
!> \param use_bvk_kernel ...
!> \param adjoint_src ...
!> \param source_adjoint_k ...
!> \param mulliken_grad_r ...
!> \param bocorr_grad_r ...
!> \param transform_backend ...
!> \param nmesh ...
!> \param input_to_grid ...
!> \param phase_roots ...
!> \param twist_phase ...
!> \param transform ...
! **************************************************************************************************
   SUBROUTINE reverse_bvk_Kkernel(self, mol, kernel, weights, phase_inverse, &
      & phase_forward, source_k, result_adjoint_k, g_onsfx, g_onsri, &
      & mulliken_scale, bocorr_scale, onsite_fx_scale, onsite_ri_scale, &
      & use_bvk_kernel, adjoint_src, source_adjoint_k, mulliken_grad_r, &
      & bocorr_grad_r, transform_backend, nmesh, input_to_grid, phase_roots, &
      & twist_phase, transform)
      class(exchange_fock), INTENT(in) :: self
      TYPE(structure_type), INTENT(in) :: mol
      TYPE(exchange_bvk_kernel), INTENT(in) :: kernel
      REAL(wp), INTENT(in) :: weights(:)
      COMPLEX(wp), INTENT(in) :: phase_inverse(:, :), phase_forward(:, :)
      COMPLEX(wp), contiguous, TARGET, INTENT(in) :: source_k(:, :, :)
      COMPLEX(wp), contiguous, TARGET, INTENT(in) :: result_adjoint_k(:, :, :)
      REAL(wp), INTENT(in) :: g_onsfx(:, :, :), g_onsri(:, :)
      REAL(wp), INTENT(in) :: mulliken_scale, bocorr_scale
      REAL(wp), INTENT(in) :: onsite_fx_scale, onsite_ri_scale
      LOGICAL, INTENT(in) :: use_bvk_kernel, adjoint_src
      COMPLEX(wp), contiguous, TARGET, INTENT(out) :: source_adjoint_k(:, :, :)
      REAL(wp), INTENT(inout) :: mulliken_grad_r(:, :, :), bocorr_grad_r(:, :, :)
      INTEGER, INTENT(in), OPTIONAL :: transform_backend, nmesh(3), &
         & input_to_grid(:)
      COMPLEX(wp), INTENT(in), OPTIONAL :: phase_roots(:, :), twist_phase(:)
      class(exchange_transform_type), INTENT(inout), OPTIONAL :: transform

      INTEGER :: backend, icell, ik, npair, origin
      COMPLEX(wp), ALLOCATABLE, TARGET :: result_adjoint_r(:, :, :), &
         & source_adjoint_r(:, :, :), source_r(:, :, :)
      COMPLEX(wp), POINTER :: result_adjoint_k_2d(:, :), &
         & result_adjoint_r_2d(:, :), source_adjoint_k_2d(:, :), &
         & source_adjoint_r_2d(:, :), source_k_2d(:, :), source_r_2d(:, :)

      ALLOCATE (source_r(self%nao, self%nao, SIZE(kernel%reps, 2)), &
         & result_adjoint_r(self%nao, self%nao, SIZE(kernel%reps, 2)), &
         & source_adjoint_r(self%nao, self%nao, SIZE(kernel%reps, 2)), &
         & source=(0.0_wp, 0.0_wp))
      npair = self%nao*self%nao
      source_k_2d(1:npair, 1:SIZE(source_k, 3)) => source_k
      source_r_2d(1:npair, 1:SIZE(source_r, 3)) => source_r
      result_adjoint_k_2d(1:npair, 1:SIZE(result_adjoint_k, 3)) &
         & => result_adjoint_k
      result_adjoint_r_2d(1:npair, 1:SIZE(result_adjoint_r, 3)) &
         & => result_adjoint_r
      source_adjoint_r_2d(1:npair, 1:SIZE(source_adjoint_r, 3)) &
         & => source_adjoint_r
      source_adjoint_k_2d(1:npair, 1:SIZE(source_adjoint_k, 3)) &
         & => source_adjoint_k

      backend = bvk_transform_dense
      IF (PRESENT(transform_backend)) backend = transform_backend
      IF (PRESENT(transform)) THEN
         CALL transform%apply(source_k_2d, source_r_2d, exchange_k_to_r)
         CALL transform%apply(result_adjoint_k_2d, result_adjoint_r_2d, exchange_r_to_k_adjoint)
      ELSE IF (is_compact_transform_backend(backend) .AND. PRESENT(nmesh) &
         & .AND. PRESENT(input_to_grid) .AND. PRESENT(phase_roots) &
         & .AND. PRESENT(twist_phase)) THEN
         SELECT CASE (backend)
         CASE (bvk_transform_mixed_radix_fft)
            CALL bvk_mixed_radix_k_to_r(source_k_2d, source_r_2d, nmesh, &
               & input_to_grid, kernel%reps, phase_roots, twist_phase, .TRUE.)
            CALL bvk_mixed_radix_k_to_r(result_adjoint_k_2d, &
               & result_adjoint_r_2d, nmesh, input_to_grid, kernel%reps, &
               & phase_roots, twist_phase, .FALSE.)
         CASE default
            CALL bvk_separable_k_to_r(source_k_2d, source_r_2d, nmesh, &
               & input_to_grid, kernel%reps, phase_roots, twist_phase, .TRUE.)
            CALL bvk_separable_k_to_r(result_adjoint_k_2d, result_adjoint_r_2d, &
               & nmesh, input_to_grid, kernel%reps, phase_roots, twist_phase, &
               & .FALSE.)
         END SELECT
      ELSE
         source_r_2d = MATMUL(source_k_2d, phase_inverse)
         result_adjoint_r_2d = MATMUL(result_adjoint_k_2d, &
            & CONJG(TRANSPOSE(phase_forward)))
      END IF

      origin = 0
      DO icell = 1, SIZE(kernel%reps, 2)
         IF (ALL(kernel%reps(:, icell) == 0)) origin = icell
      END DO

      IF (use_bvk_kernel) THEN
         DO icell = 1, SIZE(kernel%reps, 2)
            IF (mulliken_scale /= 0.0_wp) THEN
               CALL shell_parameter_gradient_complex(self%nsh, self%nao_sh, &
                  & self%iao_sh, result_adjoint_r(:, :, icell), &
                  & source_r(:, :, icell), mulliken_scale, &
                  & mulliken_grad_r(:, :, icell))
               CALL shell_hadamard_add_complex(self%nsh, self%nao_sh, &
                  & self%iao_sh, kernel%g_mulliken_r(:, :, icell), &
                  & result_adjoint_r(:, :, icell), mulliken_scale, &
                  & source_adjoint_r(:, :, icell))
            END IF
            IF (bocorr_scale /= 0.0_wp) THEN
               CALL atom_parameter_gradient_complex(mol%nat, mol%id, &
                  & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, &
                  & result_adjoint_r(:, :, icell), source_r(:, :, icell), &
                  & bocorr_scale, bocorr_grad_r(:, :, icell))
               CALL atom_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
                  & self%nao_sh, self%ish_at, self%iao_sh, &
                  & kernel%g_bocorr_r(:, :, icell), &
                  & result_adjoint_r(:, :, icell), bocorr_scale, &
                  & source_adjoint_r(:, :, icell))
            END IF
         END DO
      END IF

      IF (origin > 0) THEN
         IF (onsite_fx_scale /= 0.0_wp) THEN
            CALL onsite_fx_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
                & self%nao_sh, self%ish_at, self%iao_sh, g_onsfx, &
                & result_adjoint_r(:, :, origin), onsite_fx_scale, &
                & source_adjoint_r(:, :, origin), adjoint_src=adjoint_src)
         END IF
         IF (onsite_ri_scale /= 0.0_wp) THEN
            CALL onsite_ri_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
                & self%nao_sh, self%ish_at, self%iao_sh, g_onsri, &
                & result_adjoint_r(:, :, origin), onsite_ri_scale, &
                & source_adjoint_r(:, :, origin))
         END IF
      END IF

      IF (PRESENT(transform)) THEN
         CALL transform%apply(source_adjoint_r_2d, source_adjoint_k_2d, exchange_k_to_r_adjoint)
      ELSE IF (is_compact_transform_backend(backend) .AND. PRESENT(nmesh) &
         & .AND. PRESENT(input_to_grid) .AND. PRESENT(phase_roots) &
         & .AND. PRESENT(twist_phase)) THEN
         SELECT CASE (backend)
         CASE (bvk_transform_mixed_radix_fft)
            CALL bvk_mixed_radix_r_to_k(source_adjoint_r_2d, &
               & source_adjoint_k_2d, nmesh, input_to_grid, kernel%reps, &
               & phase_roots, twist_phase)
         CASE default
            CALL bvk_separable_r_to_k(source_adjoint_r_2d, source_adjoint_k_2d, &
               & nmesh, input_to_grid, kernel%reps, phase_roots, twist_phase)
         END SELECT
      ELSE
         source_adjoint_k_2d = MATMUL(source_adjoint_r_2d, phase_forward)
      END IF
      IF (.NOT. PRESENT(transform)) THEN
         DO ik = 1, SIZE(weights)
            source_adjoint_k(:, :, ik) = weights(ik)*source_adjoint_k(:, :, ik)
         END DO
      END IF

   END SUBROUTINE reverse_bvk_Kkernel

! **************************************************************************************************
!> \brief ...
!> \param self ...
!> \param mol ...
!> \param cache ...
!> \param plan ...
!> \param kernel ...
!> \param kfrac ...
!> \param weights ...
!> \param density ...
!> \param overlap ...
!> \param overlap_grad ...
!> \param mulliken_grad_r ...
!> \param bocorr_grad_r ...
!> \param transform_backend ...
!> \param transform ...
! **************************************************************************************************
   SUBROUTINE get_KGrad_kmesh(self, mol, cache, plan, kernel, kfrac, weights, density, &
      & overlap, overlap_grad, mulliken_grad_r, bocorr_grad_r, transform_backend, transform)
      class(exchange_fock), INTENT(in) :: self
      TYPE(structure_type), INTENT(in) :: mol
      TYPE(exchange_cache), TARGET, INTENT(in) :: cache
      TYPE(cp2k_stream_plan_type), INTENT(in), TARGET :: plan
      TYPE(exchange_bvk_kernel), INTENT(in) :: kernel
      REAL(wp), INTENT(in) :: kfrac(:, :), weights(:)
      COMPLEX(wp), INTENT(in) :: density(:, :, :, :), overlap(:, :, :)
      COMPLEX(wp), INTENT(out) :: overlap_grad(:, :, :)
      REAL(wp), INTENT(out) :: mulliken_grad_r(:, :, :), bocorr_grad_r(:, :, :)
      INTEGER, INTENT(in), OPTIONAL :: transform_backend
      class(exchange_transform_type), INTENT(inout), OPTIONAL :: transform

      INTEGER :: iao, backend, icell, ik, spin
      LOGICAL :: use_cached_phases, use_compact
      REAL(wp) :: angle, spin_factor
      COMPLEX(wp) :: phase
      COMPLEX(wp), ALLOCATABLE :: amat(:, :, :), &
         & bdiag(:), bdiag_p(:), bdiag_sp(:), bdiag_sps(:), &
         & btmp(:, :, :), dmat(:, :, :), ddensity(:, :, :), &
         & gdiag_p(:), gdiag_sp(:), gdiag_sps(:), tmat(:, :, :), &
         & umat(:, :, :), vmat(:, :, :), work(:, :, :), xmat(:, :, :)
      COMPLEX(wp), ALLOCATABLE, TARGET :: phase_forward_local(:, :), phase_inverse_local(:, :)
      INTEGER, ALLOCATABLE :: plan_input_to_grid(:)
      COMPLEX(wp), ALLOCATABLE :: plan_phase_roots(:, :), plan_twist_phase(:)
      COMPLEX(wp), contiguous, POINTER :: phase_forward(:, :), phase_inverse(:, :)

      ALLOCATE (amat(self%nao, self%nao, SIZE(weights)), &
         & btmp(self%nao, self%nao, SIZE(weights)), &
         & dmat(self%nao, self%nao, SIZE(weights)), &
         & ddensity(self%nao, self%nao, SIZE(weights)), &
         & tmat(self%nao, self%nao, SIZE(weights)), &
         & umat(self%nao, self%nao, SIZE(weights)), &
         & vmat(self%nao, self%nao, SIZE(weights)), &
         & work(self%nao, self%nao, SIZE(weights)), &
         & xmat(self%nao, self%nao, SIZE(weights)), &
         & bdiag(self%nao), bdiag_p(self%nao), bdiag_sp(self%nao), &
         & bdiag_sps(self%nao), gdiag_p(self%nao), gdiag_sp(self%nao), &
         & gdiag_sps(self%nao), source=(0.0_wp, 0.0_wp))
      backend = bvk_transform_dense
      IF (PRESENT(transform_backend)) backend = transform_backend
      use_compact = is_compact_transform_backend(backend) .AND. &
         & plan%plan_matches(self, mol, kernel, kfrac, weights, &
         & require_phases=.FALSE.) .AND. ALLOCATED(plan%bvk_phase_roots) &
         & .AND. ALLOCATED(plan%bvk_twist_phase)
      IF (PRESENT(transform)) use_compact = .FALSE.
      IF (PRESENT(transform)) THEN
         ALLOCATE (plan_input_to_grid(1), source=1)
         ALLOCATE (plan_phase_roots(1, 3), source=(1.0_wp, 0.0_wp))
         ALLOCATE (plan_twist_phase(1), source=(1.0_wp, 0.0_wp))
         ALLOCATE (phase_forward_local(1, 1), phase_inverse_local(1, 1), &
            & source=(0.0_wp, 0.0_wp))
         phase_forward => phase_forward_local
         phase_inverse => phase_inverse_local
         use_cached_phases = .FALSE.
      ELSE IF (use_compact) THEN
         ALLOCATE (plan_input_to_grid, source=plan%bvk_input_to_grid)
         ALLOCATE (plan_phase_roots, source=plan%bvk_phase_roots)
         ALLOCATE (plan_twist_phase, source=plan%bvk_twist_phase)
         ALLOCATE (phase_forward_local(1, 1), phase_inverse_local(1, 1), &
            & source=(0.0_wp, 0.0_wp))
         phase_forward => phase_forward_local
         phase_inverse => phase_inverse_local
         use_cached_phases = .FALSE.
      ELSE
         backend = bvk_transform_dense
         ALLOCATE (plan_input_to_grid(1), source=1)
         ALLOCATE (plan_phase_roots(1, 3), source=(1.0_wp, 0.0_wp))
         ALLOCATE (plan_twist_phase(1), source=(1.0_wp, 0.0_wp))
         use_cached_phases = plan%plan_matches(self, mol, kernel, kfrac, weights)
      END IF
      IF (.NOT. use_compact .AND. use_cached_phases) THEN
         phase_forward => plan%bvk_phase_forward
         phase_inverse => plan%bvk_phase_inverse
      ELSE IF (.NOT. use_compact .AND. .NOT. PRESENT(transform)) THEN
         ALLOCATE (phase_forward_local(SIZE(kernel%reps, 2), SIZE(weights)), &
            & phase_inverse_local(SIZE(weights), SIZE(kernel%reps, 2)))
         phase_forward => phase_forward_local
         phase_inverse => phase_inverse_local
         DO ik = 1, SIZE(weights)
            DO icell = 1, SIZE(kernel%reps, 2)
               angle = 2.0_wp*pi*DOT_PRODUCT(kfrac(:, ik), &
                  & REAL(kernel%reps(:, icell), wp))
               phase = EXP(CMPLX(0.0_wp, angle, wp))
               phase_forward(icell, ik) = phase
               phase_inverse(ik, icell) = weights(ik)*CONJG(phase)
            END DO
         END DO
      END IF

      spin_factor = 0.5_wp
      IF (SIZE(density, 3) > 1) spin_factor = 1.0_wp
      overlap_grad = (0.0_wp, 0.0_wp)
      mulliken_grad_r = 0.0_wp
      bocorr_grad_r = 0.0_wp

      DO spin = 1, SIZE(density, 3)
         bdiag_p = (0.0_wp, 0.0_wp)
         bdiag_sp = (0.0_wp, 0.0_wp)
         bdiag_sps = (0.0_wp, 0.0_wp)
         gdiag_p = (0.0_wp, 0.0_wp)
         gdiag_sp = (0.0_wp, 0.0_wp)
         gdiag_sps = (0.0_wp, 0.0_wp)

         DO ik = 1, SIZE(weights)
            xmat(:, :, ik) = MATMUL(overlap(:, :, ik), &
               & density(:, :, spin, ik))
            dmat(:, :, ik) = 0.5_wp*MATMUL(xmat(:, :, ik), overlap(:, :, ik))
            DO iao = 1, self%nao
               gdiag_p(iao) = gdiag_p(iao) &
                  & + weights(ik)*density(iao, iao, spin, ik)
               gdiag_sp(iao) = gdiag_sp(iao) &
                  & + weights(ik)*xmat(iao, iao, ik)
               gdiag_sps(iao) = gdiag_sps(iao) &
                  & + 2.0_wp*weights(ik)*dmat(iao, iao, ik)
            END DO
         END DO
         CALL onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
            & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, &
            & gdiag_p, bdiag_p)
         CALL onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
            & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, &
            & gdiag_sp, bdiag_sp)
         CALL onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
            & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, &
            & gdiag_sps, bdiag_sps)
         gdiag_p = bdiag_p
         gdiag_sp = bdiag_sp
         gdiag_sps = bdiag_sps
         bdiag_sp = (0.0_wp, 0.0_wp)
         bdiag_sps = (0.0_wp, 0.0_wp)

         CALL apply_bvk_Kkernel(self, mol, kernel, phase_inverse, phase_forward, &
            & xmat, cache%g_onsfx, cache%g_onsri, 1.0_wp, -4.0_wp, &
            & 0.5_wp, -2.0_wp, .TRUE., .TRUE., tmat, backend, plan%bvk_nmesh, &
            & plan_input_to_grid, plan_phase_roots, plan_twist_phase, &
            & transform=transform)
         CALL apply_bvk_Kkernel(self, mol, kernel, phase_inverse, phase_forward, &
            & dmat, cache%g_onsfx, cache%g_onsri, 1.0_wp, 0.0_wp, &
            & 0.5_wp, -2.0_wp, .TRUE., .FALSE., umat, backend, plan%bvk_nmesh, &
            & plan_input_to_grid, plan_phase_roots, plan_twist_phase, &
            & transform=transform)
         CALL apply_bvk_Kkernel(self, mol, kernel, phase_inverse, phase_forward, &
            & density(:, :, spin, :), cache%g_onsfx, cache%g_onsri, &
            & 1.0_wp, 0.0_wp, 0.5_wp, -2.0_wp, .TRUE., .FALSE., vmat, backend, &
            & plan%bvk_nmesh, plan_input_to_grid, plan_phase_roots, &
            & plan_twist_phase, &
            & transform=transform)

         DO ik = 1, SIZE(weights)
            work(:, :, ik) = tmat(:, :, ik)
            DO iao = 1, self%nao
               work(:, iao, ik) = work(:, iao, ik) &
                  & + 0.25_wp*gdiag_p(iao)*overlap(:, iao, ik)
            END DO
            work(:, :, ik) = work(:, :, ik) &
               & + 0.5_wp*MATMUL(overlap(:, :, ik), vmat(:, :, ik))

            ! Seed E = -(spin_factor/4) sum_k w_k Re<P_k,A_k>.
            amat(:, :, ik) = -0.25_wp*spin_factor*weights(ik) &
                            & *density(:, :, spin, ik)
            umat(:, :, ik) = amat(:, :, ik)
            btmp(:, :, ik) = MATMUL(amat(:, :, ik), &
               & CONJG(TRANSPOSE(overlap(:, :, ik))))
            overlap_grad(:, :, ik) = overlap_grad(:, :, ik) &
               & + MATMUL(CONJG(TRANSPOSE(work(:, :, ik))), amat(:, :, ik))
            DO iao = 1, self%nao
               overlap_grad(:, iao, ik) = overlap_grad(:, iao, ik) &
                  & + 0.5_wp*CONJG(gdiag_sp(iao))*amat(:, iao, ik)
               bdiag_sp(iao) = bdiag_sp(iao) &
                  & + 0.5_wp*DOT_PRODUCT(overlap(:, iao, ik), &
                  & amat(:, iao, ik))
               bdiag_sps(iao) = bdiag_sps(iao) &
                  & + 0.25_wp*amat(iao, iao, ik)
            END DO

            DO iao = 1, self%nao
               overlap_grad(:, iao, ik) = overlap_grad(:, iao, ik) &
                  & + 0.25_wp*CONJG(gdiag_p(iao))*btmp(:, iao, ik)
            END DO
            overlap_grad(:, :, ik) = overlap_grad(:, :, ik) &
               & + 0.5_wp*MATMUL(btmp(:, :, ik), &
               & CONJG(TRANSPOSE(vmat(:, :, ik))))
            vmat(:, :, ik) = 0.5_wp*MATMUL( &
               & CONJG(TRANSPOSE(overlap(:, :, ik))), btmp(:, :, ik))
         END DO

         CALL reverse_bvk_Kkernel(self, mol, kernel, weights, phase_inverse, &
            & phase_forward, dmat, umat, cache%g_onsfx, cache%g_onsri, &
            & 1.0_wp, 0.0_wp, 0.5_wp, -2.0_wp, .TRUE., .FALSE., tmat, &
            & mulliken_grad_r, bocorr_grad_r, backend, plan%bvk_nmesh, &
            & plan_input_to_grid, plan_phase_roots, plan_twist_phase, &
            & transform=transform)
         CALL reverse_bvk_Kkernel(self, mol, kernel, weights, phase_inverse, &
            & phase_forward, xmat, btmp, cache%g_onsfx, cache%g_onsri, &
            & 1.0_wp, -4.0_wp, 0.5_wp, -2.0_wp, .TRUE., .TRUE., umat, &
            & mulliken_grad_r, bocorr_grad_r, backend, plan%bvk_nmesh, &
            & plan_input_to_grid, plan_phase_roots, plan_twist_phase, &
            & transform=transform)
         CALL reverse_bvk_Kkernel(self, mol, kernel, weights, phase_inverse, &
            & phase_forward, density(:, :, spin, :), vmat, cache%g_onsfx, &
            & cache%g_onsri, 1.0_wp, 0.0_wp, 0.5_wp, -2.0_wp, .TRUE., &
            & .FALSE., ddensity, mulliken_grad_r, bocorr_grad_r, backend, &
            & plan%bvk_nmesh, plan_input_to_grid, plan_phase_roots, &
            & plan_twist_phase, &
            & transform=transform)

         ! Reverse the three BZ-averaged onsite diagonal maps.  diag(P) has no
         ! overlap response; its kernel response belongs to the shell potential.
         CALL onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
            & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, &
            & bdiag_sp, bdiag)
         DO ik = 1, SIZE(weights)
            DO iao = 1, self%nao
               umat(iao, iao, ik) = umat(iao, iao, ik) &
                  & + weights(ik)*bdiag(iao)
            END DO
         END DO
         CALL onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
            & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, &
            & bdiag_sps, bdiag)
         DO ik = 1, SIZE(weights)
            DO iao = 1, self%nao
               umat(iao, :, ik) = umat(iao, :, ik) &
                  & + weights(ik)*bdiag(iao)*overlap(iao, :, ik)
               overlap_grad(:, iao, ik) = overlap_grad(:, iao, ik) &
                  & + weights(ik)*bdiag(iao)*CONJG(xmat(iao, :, ik))
            END DO

            ! D_k = X_k S_k/2 and X_k = S_k P_k.
            umat(:, :, ik) = umat(:, :, ik) &
               & + 0.5_wp*MATMUL(tmat(:, :, ik), &
               & CONJG(TRANSPOSE(overlap(:, :, ik))))
            overlap_grad(:, :, ik) = overlap_grad(:, :, ik) &
               & + 0.5_wp*MATMUL(CONJG(TRANSPOSE(xmat(:, :, ik))), &
               & tmat(:, :, ik))
            overlap_grad(:, :, ik) = overlap_grad(:, :, ik) &
               & + MATMUL(umat(:, :, ik), &
               & CONJG(TRANSPOSE(density(:, :, spin, ik))))
         END DO
      END DO

      DO ik = 1, SIZE(weights)
         overlap_grad(:, :, ik) = 0.5_wp*(overlap_grad(:, :, ik) &
            & + CONJG(TRANSPOSE(overlap_grad(:, :, ik))))/weights(ik)
      END DO

   END SUBROUTINE get_KGrad_kmesh

! **************************************************************************************************
!> \brief ...
!> \param nsh ...
!> \param nao_sh ...
!> \param iao_sh ...
!> \param g_sh ...
!> \param src ...
!> \param alpha ...
!> \param dst ...
!> \param adjoint_src ...
! **************************************************************************************************
   SUBROUTINE shell_hadamard_add_complex(nsh, nao_sh, iao_sh, g_sh, src, &
      & alpha, dst, adjoint_src)
      INTEGER, INTENT(in)                                :: nsh, nao_sh(:), iao_sh(:)
      REAL(wp), INTENT(in)                               :: g_sh(:, :)
      COMPLEX(wp), INTENT(in)                            :: src(:, :)
      REAL(wp), INTENT(in)                               :: alpha
      COMPLEX(wp), INTENT(inout)                         :: dst(:, :)
      LOGICAL, INTENT(in), OPTIONAL                      :: adjoint_src

      INTEGER                                            :: iao, ii, ish, jao, jj, jsh, ni, nj
      LOGICAL                                            :: adjoint
      REAL(wp)                                           :: scale

      adjoint = .FALSE.
      IF (PRESENT(adjoint_src)) adjoint = adjoint_src
      DO ish = 1, nsh
         ii = iao_sh(ish)
         ni = nao_sh(ish)
         DO jsh = 1, nsh
            jj = iao_sh(jsh)
            nj = nao_sh(jsh)
            scale = alpha*g_sh(jsh, ish)
            IF (ABS(scale) < EPSILON(1.0_wp)) CYCLE
            IF (.NOT. adjoint) THEN
               dst(jj + 1:jj + nj, ii + 1:ii + ni) = dst(jj + 1:jj + nj, ii + 1:ii + ni) &
                  & + scale*src(jj + 1:jj + nj, ii + 1:ii + ni)
            ELSE
               DO iao = 1, ni
                  DO jao = 1, nj
                     dst(jj + jao, ii + iao) = dst(jj + jao, ii + iao) &
                        & + scale*CONJG(src(ii + iao, jj + jao))
                  END DO
               END DO
            END IF
         END DO
      END DO
   END SUBROUTINE shell_hadamard_add_complex

! **************************************************************************************************
!> \brief ...
!> \param nsh ...
!> \param nao_sh ...
!> \param iao_sh ...
!> \param adjoint ...
!> \param source ...
!> \param alpha ...
!> \param gradient ...
! **************************************************************************************************
   SUBROUTINE shell_parameter_gradient_complex(nsh, nao_sh, iao_sh, adjoint, &
      & source, alpha, gradient)
      INTEGER, INTENT(in)                                :: nsh, nao_sh(:), iao_sh(:)
      COMPLEX(wp), INTENT(in)                            :: adjoint(:, :), source(:, :)
      REAL(wp), INTENT(in)                               :: alpha
      REAL(wp), INTENT(inout)                            :: gradient(:, :)

      INTEGER                                            :: iao, ii, ish, jao, jj, jsh, ni, nj

      DO ish = 1, nsh
         ii = iao_sh(ish)
         ni = nao_sh(ish)
         DO jsh = 1, nsh
            jj = iao_sh(jsh)
            nj = nao_sh(jsh)
            DO iao = 1, ni
               DO jao = 1, nj
                  gradient(jj + jao, ii + iao) = gradient(jj + jao, ii + iao) &
                     & + alpha*REAL(CONJG(adjoint(jj + jao, ii + iao)) &
                     & *source(jj + jao, ii + iao), wp)
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE shell_parameter_gradient_complex

! **************************************************************************************************
!> \brief ...
!> \param nat ...
!> \param id ...
!> \param nsh_id ...
!> \param nao_sh ...
!> \param ish_at ...
!> \param iao_sh ...
!> \param g_at ...
!> \param src ...
!> \param alpha ...
!> \param dst ...
!> \param adjoint_src ...
! **************************************************************************************************
   SUBROUTINE atom_hadamard_add_complex(nat, id, nsh_id, nao_sh, ish_at, iao_sh, &
      & g_at, src, alpha, dst, adjoint_src)
      INTEGER, INTENT(in)                                :: nat, id(:), nsh_id(:), nao_sh(:), &
                                                            ish_at(:), iao_sh(:)
      REAL(wp), INTENT(in)                               :: g_at(:, :)
      COMPLEX(wp), INTENT(in)                            :: src(:, :)
      REAL(wp), INTENT(in)                               :: alpha
      COMPLEX(wp), INTENT(inout)                         :: dst(:, :)
      LOGICAL, INTENT(in), OPTIONAL                      :: adjoint_src

      INTEGER                                            :: iao, iat, ii, is, ish, izp, jao, jat, &
                                                            jj, js, jsh, jzp, ni, nj
      LOGICAL                                            :: adjoint
      REAL(wp)                                           :: scale

      adjoint = .FALSE.
      IF (PRESENT(adjoint_src)) adjoint = adjoint_src
      DO iat = 1, nat
         izp = id(iat)
         is = ish_at(iat)
         DO jat = 1, nat
            jzp = id(jat)
            js = ish_at(jat)
            scale = alpha*g_at(jat, iat)
            IF (ABS(scale) < EPSILON(1.0_wp)) CYCLE
            DO ish = 1, nsh_id(izp)
               ii = iao_sh(is + ish)
               ni = nao_sh(is + ish)
               DO jsh = 1, nsh_id(jzp)
                  jj = iao_sh(js + jsh)
                  nj = nao_sh(js + jsh)
                  IF (.NOT. adjoint) THEN
                     dst(jj + 1:jj + nj, ii + 1:ii + ni) = &
                        & dst(jj + 1:jj + nj, ii + 1:ii + ni) &
                        & + scale*src(jj + 1:jj + nj, ii + 1:ii + ni)
                  ELSE
                     DO iao = 1, ni
                        DO jao = 1, nj
                           dst(jj + jao, ii + iao) = dst(jj + jao, ii + iao) &
                              & + scale*CONJG(src(ii + iao, jj + jao))
                        END DO
                     END DO
                  END IF
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE atom_hadamard_add_complex

! **************************************************************************************************
!> \brief ...
!> \param nat ...
!> \param id ...
!> \param nsh_id ...
!> \param nao_sh ...
!> \param ish_at ...
!> \param iao_sh ...
!> \param adjoint ...
!> \param source ...
!> \param alpha ...
!> \param gradient ...
! **************************************************************************************************
   SUBROUTINE atom_parameter_gradient_complex(nat, id, nsh_id, nao_sh, ish_at, &
      & iao_sh, adjoint, source, alpha, gradient)
      INTEGER, INTENT(in)                                :: nat, id(:), nsh_id(:), nao_sh(:), &
                                                            ish_at(:), iao_sh(:)
      COMPLEX(wp), INTENT(in)                            :: adjoint(:, :), source(:, :)
      REAL(wp), INTENT(in)                               :: alpha
      REAL(wp), INTENT(inout)                            :: gradient(:, :)

      INTEGER                                            :: iao, iat, ii, is, ish, izp, jao, jat, &
                                                            jj, js, jsh, jzp, ni, nj

      DO iat = 1, nat
         izp = id(iat)
         is = ish_at(iat)
         DO jat = 1, nat
            jzp = id(jat)
            js = ish_at(jat)
            DO ish = 1, nsh_id(izp)
               ii = iao_sh(is + ish)
               ni = nao_sh(is + ish)
               DO jsh = 1, nsh_id(jzp)
                  jj = iao_sh(js + jsh)
                  nj = nao_sh(js + jsh)
                  DO iao = 1, ni
                     DO jao = 1, nj
                        gradient(jj + jao, ii + iao) = gradient(jj + jao, ii + iao) &
                           & + alpha*REAL(CONJG(adjoint(jj + jao, ii + iao)) &
                           & *source(jj + jao, ii + iao), wp)
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE atom_parameter_gradient_complex

! **************************************************************************************************
!> \brief ...
!> \param nat ...
!> \param id ...
!> \param nsh_id ...
!> \param nao_sh ...
!> \param ish_at ...
!> \param iao_sh ...
!> \param g_onsfx ...
!> \param src ...
!> \param alpha ...
!> \param dst ...
!> \param adjoint_src ...
! **************************************************************************************************
   SUBROUTINE onsite_fx_hadamard_add_complex(nat, id, nsh_id, nao_sh, ish_at, &
      & iao_sh, g_onsfx, src, alpha, dst, adjoint_src)
      INTEGER, INTENT(in)                                :: nat, id(:), nsh_id(:), nao_sh(:), &
                                                            ish_at(:), iao_sh(:)
      REAL(wp), INTENT(in)                               :: g_onsfx(:, :, :)
      COMPLEX(wp), INTENT(in)                            :: src(:, :)
      REAL(wp), INTENT(in)                               :: alpha
      COMPLEX(wp), INTENT(inout)                         :: dst(:, :)
      LOGICAL, INTENT(in), OPTIONAL                      :: adjoint_src

      INTEGER                                            :: iao, iat, ii, is, ish, izp, jao, jj, &
                                                            jsh, ni, nj
      LOGICAL                                            :: adjoint
      REAL(wp)                                           :: scale

      adjoint = .FALSE.
      IF (PRESENT(adjoint_src)) adjoint = adjoint_src
      DO iat = 1, nat
         izp = id(iat)
         is = ish_at(iat)
         DO ish = 1, nsh_id(izp)
            ii = iao_sh(is + ish)
            ni = nao_sh(is + ish)
            DO jsh = 1, nsh_id(izp)
               jj = iao_sh(is + jsh)
               nj = nao_sh(is + jsh)
               scale = alpha*g_onsfx(jsh, ish, iat)
               IF (ABS(scale) < EPSILON(1.0_wp)) CYCLE
               IF (.NOT. adjoint) THEN
                  dst(jj + 1:jj + nj, ii + 1:ii + ni) = &
                     & dst(jj + 1:jj + nj, ii + 1:ii + ni) &
                     & + scale*src(jj + 1:jj + nj, ii + 1:ii + ni)
               ELSE
                  DO iao = 1, ni
                     DO jao = 1, nj
                        dst(jj + jao, ii + iao) = dst(jj + jao, ii + iao) &
                           & + scale*CONJG(src(ii + iao, jj + jao))
                     END DO
                  END DO
               END IF
            END DO
         END DO
      END DO
   END SUBROUTINE onsite_fx_hadamard_add_complex

! **************************************************************************************************
!> \brief ...
!> \param nat ...
!> \param id ...
!> \param nsh_id ...
!> \param nao_sh ...
!> \param ish_at ...
!> \param iao_sh ...
!> \param g_onsfx ...
!> \param xvec ...
!> \param yvec ...
! **************************************************************************************************
   SUBROUTINE onsite_fx_symv_complex(nat, id, nsh_id, nao_sh, ish_at, iao_sh, &
      & g_onsfx, xvec, yvec)
      INTEGER, INTENT(in)                                :: nat, id(:), nsh_id(:), nao_sh(:), &
                                                            ish_at(:), iao_sh(:)
      REAL(wp), INTENT(in)                               :: g_onsfx(:, :, :)
      COMPLEX(wp), INTENT(in)                            :: xvec(:)
      COMPLEX(wp), INTENT(out)                           :: yvec(:)

      COMPLEX(wp)                                        :: shsum
      INTEGER                                            :: iat, ii, is, ish, izp, jj, jsh, ni, nj

      yvec = (0.0_wp, 0.0_wp)
      DO iat = 1, nat
         izp = id(iat)
         is = ish_at(iat)
         DO ish = 1, nsh_id(izp)
            ii = iao_sh(is + ish)
            ni = nao_sh(is + ish)
            shsum = SUM(xvec(ii + 1:ii + ni))
            IF (ABS(shsum) < EPSILON(1.0_wp)) CYCLE
            DO jsh = 1, nsh_id(izp)
               jj = iao_sh(is + jsh)
               nj = nao_sh(is + jsh)
               yvec(jj + 1:jj + nj) = yvec(jj + 1:jj + nj) &
                  & + g_onsfx(jsh, ish, iat)*shsum
            END DO
         END DO
      END DO
   END SUBROUTINE onsite_fx_symv_complex

! **************************************************************************************************
!> \brief ...
!> \param nat ...
!> \param id ...
!> \param nsh_id ...
!> \param nao_sh ...
!> \param ish_at ...
!> \param iao_sh ...
!> \param g_onsri ...
!> \param src ...
!> \param alpha ...
!> \param dst ...
! **************************************************************************************************
   SUBROUTINE onsite_ri_hadamard_add_complex(nat, id, nsh_id, nao_sh, ish_at, &
      & iao_sh, g_onsri, src, alpha, dst)
      INTEGER, INTENT(in)                                :: nat, id(:), nsh_id(:), nao_sh(:), &
                                                            ish_at(:), iao_sh(:)
      REAL(wp), INTENT(in)                               :: g_onsri(:, :)
      COMPLEX(wp), INTENT(in)                            :: src(:, :)
      REAL(wp), INTENT(in)                               :: alpha
      COMPLEX(wp), INTENT(inout)                         :: dst(:, :)

      INTEGER                                            :: iat, ii, is, ish, izp, ni
      REAL(wp)                                           :: scale

      DO iat = 1, nat
         izp = id(iat)
         is = ish_at(iat)
         DO ish = 1, nsh_id(izp)
            ii = iao_sh(is + ish)
            ni = nao_sh(is + ish)
            scale = alpha*g_onsri(ish, iat)
            IF (ABS(scale) < EPSILON(1.0_wp)) CYCLE
            dst(ii + 1:ii + ni, ii + 1:ii + ni) = &
               & dst(ii + 1:ii + ni, ii + 1:ii + ni) &
               & + scale*src(ii + 1:ii + ni, ii + 1:ii + ni)
         END DO
      END DO
   END SUBROUTINE onsite_ri_hadamard_add_complex

! **************************************************************************************************
!> \brief ...
!> \param nat ...
!> \param id ...
!> \param nsh_id ...
!> \param nao_sh ...
!> \param ish_at ...
!> \param iao_sh ...
!> \param adjoint ...
!> \param src ...
!> \param alpha ...
!> \param gradient ...
!> \param adjoint_src ...
! **************************************************************************************************
   SUBROUTINE onsite_fx_parameter_gradient_complex(nat, id, nsh_id, nao_sh, &
      & ish_at, iao_sh, adjoint, src, alpha, gradient, adjoint_src)
      INTEGER, INTENT(in)                                :: nat, id(:), nsh_id(:), nao_sh(:), &
                                                            ish_at(:), iao_sh(:)
      COMPLEX(wp), INTENT(in)                            :: adjoint(:, :), src(:, :)
      REAL(wp), INTENT(in)                               :: alpha
      REAL(wp), INTENT(inout)                            :: gradient(:, :, :)
      LOGICAL, INTENT(in), OPTIONAL                      :: adjoint_src

      INTEGER                                            :: iao, iat, ii, is, ish, izp, jao, jj, &
                                                            jsh, ni, nj
      LOGICAL                                            :: take_adjoint
      REAL(wp)                                           :: value

      take_adjoint = .FALSE.
      IF (PRESENT(adjoint_src)) take_adjoint = adjoint_src
      DO iat = 1, nat
         izp = id(iat)
         is = ish_at(iat)
         DO ish = 1, nsh_id(izp)
            ii = iao_sh(is + ish)
            ni = nao_sh(is + ish)
            DO jsh = 1, nsh_id(izp)
               jj = iao_sh(is + jsh)
               nj = nao_sh(is + jsh)
               value = 0.0_wp
               IF (.NOT. take_adjoint) THEN
                  value = REAL(SUM(CONJG(adjoint(jj + 1:jj + nj, ii + 1:ii + ni)) &
                                  & *src(jj + 1:jj + nj, ii + 1:ii + ni)), wp)
               ELSE
                  DO iao = 1, ni
                     DO jao = 1, nj
                        value = value + REAL(CONJG(adjoint(jj + jao, ii + iao)) &
                                            & *CONJG(src(ii + iao, jj + jao)), wp)
                     END DO
                  END DO
               END IF
               gradient(jsh, ish, iat) = gradient(jsh, ish, iat) + alpha*value
            END DO
         END DO
      END DO
   END SUBROUTINE onsite_fx_parameter_gradient_complex

! **************************************************************************************************
!> \brief ...
!> \param nat ...
!> \param id ...
!> \param nsh_id ...
!> \param nao_sh ...
!> \param ish_at ...
!> \param iao_sh ...
!> \param adjoint ...
!> \param src ...
!> \param alpha ...
!> \param gradient ...
! **************************************************************************************************
   SUBROUTINE onsite_ri_parameter_gradient_complex(nat, id, nsh_id, nao_sh, &
      & ish_at, iao_sh, adjoint, src, alpha, gradient)
      INTEGER, INTENT(in)                                :: nat, id(:), nsh_id(:), nao_sh(:), &
                                                            ish_at(:), iao_sh(:)
      COMPLEX(wp), INTENT(in)                            :: adjoint(:, :), src(:, :)
      REAL(wp), INTENT(in)                               :: alpha
      REAL(wp), INTENT(inout)                            :: gradient(:, :)

      INTEGER                                            :: iat, ii, is, ish, izp, ni

      DO iat = 1, nat
         izp = id(iat)
         is = ish_at(iat)
         DO ish = 1, nsh_id(izp)
            ii = iao_sh(is + ish)
            ni = nao_sh(is + ish)
            gradient(ish, iat) = gradient(ish, iat) + alpha*REAL(SUM( &
               & CONJG(adjoint(ii + 1:ii + ni, ii + 1:ii + ni)) &
               & *src(ii + 1:ii + ni, ii + 1:ii + ni)), wp)
         END DO
      END DO
   END SUBROUTINE onsite_ri_parameter_gradient_complex

! **************************************************************************************************
!> \brief ...
!> \param nat ...
!> \param id ...
!> \param nsh_id ...
!> \param nao_sh ...
!> \param ish_at ...
!> \param iao_sh ...
!> \param adjoint ...
!> \param src ...
!> \param gradient ...
! **************************************************************************************************
   SUBROUTINE onsite_fx_symv_parameter_gradient_complex(nat, id, nsh_id, &
      & nao_sh, ish_at, iao_sh, adjoint, src, gradient)
      INTEGER, INTENT(in)                                :: nat, id(:), nsh_id(:), nao_sh(:), &
                                                            ish_at(:), iao_sh(:)
      COMPLEX(wp), INTENT(in)                            :: adjoint(:), src(:)
      REAL(wp), INTENT(inout)                            :: gradient(:, :, :)

      COMPLEX(wp)                                        :: adjoint_sum, source_sum
      INTEGER                                            :: iat, ii, is, ish, izp, jj, jsh, ni, nj

      DO iat = 1, nat
         izp = id(iat)
         is = ish_at(iat)
         DO ish = 1, nsh_id(izp)
            ii = iao_sh(is + ish)
            ni = nao_sh(is + ish)
            source_sum = SUM(src(ii + 1:ii + ni))
            DO jsh = 1, nsh_id(izp)
               jj = iao_sh(is + jsh)
               nj = nao_sh(is + jsh)
               adjoint_sum = SUM(adjoint(jj + 1:jj + nj))
               gradient(jsh, ish, iat) = gradient(jsh, ish, iat) &
                  & + REAL(CONJG(adjoint_sum)*source_sum, wp)
            END DO
         END DO
      END DO
   END SUBROUTINE onsite_fx_symv_parameter_gradient_complex
! **************************************************************************************************
!> \brief ...
!> \param self ...
!> \param mol ...
!> \param cache ...
!> \param density ...
!> \param overlap ...
!> \param fock ...
!> \param vsh ...
!> \param energy ...
! **************************************************************************************************
   SUBROUTINE get_KFock_kpoint(self, mol, cache, density, overlap, fock, vsh, energy)
      !> Instance of the Fock exchange container
      class(exchange_fock), INTENT(in) :: self
      !> Molecular structure data
      TYPE(structure_type), INTENT(in) :: mol
      !> Reusable exchange kernels and charge derivatives
      TYPE(exchange_cache), INTENT(in) :: cache
      !> Complex Hermitian density matrix for all spin channels
      COMPLEX(wp), INTENT(in) :: density(:, :, :)
      !> Complex Hermitian overlap matrix
      COMPLEX(wp), INTENT(in) :: overlap(:, :)
      !> Complex Hermitian exchange Fock matrix
      COMPLEX(wp), INTENT(out) :: fock(:, :, :)
      !> Charge-dependent shell potential
      REAL(wp), INTENT(out) :: vsh(:)
      !> Optional exchange energy for this k-point block
      REAL(wp), INTENT(out), OPTIONAL :: energy

      INTEGER :: iat, izp, is, ish, iq, spin
      REAL(wp), ALLOCATABLE :: dg_onsfx(:, :, :), dg_onsri(:, :)
      REAL(wp), ALLOCATABLE :: zero_at(:, :), zero_sh(:, :)
      COMPLEX(wp), ALLOCATABLE :: dfock(:, :, :)
      REAL(wp) :: tmp_energy

      CALL build_KFock_complex(self, mol, density, overlap, cache%g_mulliken, &
         & cache%g_bocorr, cache%g_onsfx, cache%g_onsri, fock)

      tmp_energy = 0.0_wp
      DO spin = 1, SIZE(density, 3)
         tmp_energy = tmp_energy + 0.5_wp*REAL(SUM(CONJG(density(:, :, spin)) &
                                                  & *fock(:, :, spin)), wp)
      END DO
      IF (PRESENT(energy)) energy = tmp_energy

      ALLOCATE (zero_sh(self%nsh, self%nsh), zero_at(mol%nat, mol%nat), &
         & dg_onsfx(self%maxsh, self%maxsh, mol%nat), &
         & dg_onsri(self%maxsh, mol%nat), source=0.0_wp)
      ALLOCATE (dfock(self%nao, self%nao, SIZE(density, 3)), &
         & source=(0.0_wp, 0.0_wp))
      vsh = 0.0_wp

      ! The onsite kernels are the only exchange kernels that depend explicitly
      ! on shell charges.  Apply their sparse derivative through the same linear
      ! Fock map, then contract with P.  This is intentionally an exact reference
      ! implementation; the contractions can be fused after validation in CP2K.
      DO iat = 1, mol%nat
         izp = mol%id(iat)
         is = self%ish_at(iat)
         DO ish = 1, self%nsh_id(izp)
            iq = is + ish
            dg_onsfx = 0.0_wp
            dg_onsri = 0.0_wp
            dg_onsfx(:, :, iat) = cache%dgdq_onsfx(:, :, iq)
            dg_onsri(:, iat) = cache%dgdq_onsri(:, iq)
            CALL build_KFock_complex(self, mol, density, overlap, zero_sh, &
               & zero_at, dg_onsfx, dg_onsri, dfock)
            DO spin = 1, SIZE(density, 3)
               vsh(iq) = vsh(iq) + 0.5_wp*REAL(SUM( &
                  & CONJG(density(:, :, spin))*dfock(:, :, spin)), wp)
            END DO
         END DO
      END DO

   END SUBROUTINE get_KFock_kpoint

! **************************************************************************************************
!> \brief ...
!> \param self ...
!> \param mol ...
!> \param cache ...
!> \param density ...
!> \param overlap ...
!> \param mulliken_grad ...
!> \param bocorr_grad ...
!> \param overlap_grad ...
! **************************************************************************************************
   SUBROUTINE get_KGrad_kpoint(self, mol, cache, density, overlap, &
      & mulliken_grad, bocorr_grad, overlap_grad)
      class(exchange_fock), INTENT(in) :: self
      TYPE(structure_type), INTENT(in) :: mol
      TYPE(exchange_cache), INTENT(in) :: cache
      COMPLEX(wp), INTENT(in) :: density(:, :, :), overlap(:, :)
      REAL(wp), INTENT(out) :: mulliken_grad(:, :), bocorr_grad(:, :)
      COMPLEX(wp), INTENT(out) :: overlap_grad(:, :)

      INTEGER :: iao, spin
      REAL(wp) :: spin_factor
      REAL(wp), ALLOCATABLE :: mulliken_raw(:, :), bocorr_raw(:, :)
      COMPLEX(wp), ALLOCATABLE :: amat(:, :), bamat(:, :), bdiag(:), &
         & bdiag_p(:), bdiag_sp(:), bdiag_sps(:), bmat(:, :), btmp(:, :), &
         & dmat(:, :), gdiag_p(:), gdiag_sp(:), gdiag_sps(:), tmat(:, :), &
         & t2mat(:, :), umat(:, :), vmat(:, :), xmat(:, :)

      ALLOCATE (amat(self%nao, self%nao), bamat(self%nao, self%nao), &
         & bmat(self%nao, self%nao), btmp(self%nao, self%nao), &
         & dmat(self%nao, self%nao), tmat(self%nao, self%nao), &
         & t2mat(self%nao, self%nao), umat(self%nao, self%nao), &
         & vmat(self%nao, self%nao), xmat(self%nao, self%nao), &
         & bdiag(self%nao), bdiag_p(self%nao), bdiag_sp(self%nao), &
         & bdiag_sps(self%nao), gdiag_p(self%nao), gdiag_sp(self%nao), &
         & gdiag_sps(self%nao), source=(0.0_wp, 0.0_wp))
      ALLOCATE (mulliken_raw(self%nao, self%nao), &
         & bocorr_raw(self%nao, self%nao), source=0.0_wp)

      spin_factor = 0.5_wp
      IF (SIZE(density, 3) > 1) spin_factor = 1.0_wp
      overlap_grad = (0.0_wp, 0.0_wp)

      DO spin = 1, SIZE(density, 3)
         ! Forward intermediates of build_KFock_complex, kept explicitly for the
         ! reverse sweep below.
         xmat = MATMUL(overlap, density(:, :, spin))
         dmat = 0.5_wp*MATMUL(xmat, overlap)
         CALL onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
            & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, &
            & [(density(iao, iao, spin), iao=1, self%nao)], gdiag_p)
         CALL onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
            & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, &
            & [(xmat(iao, iao), iao=1, self%nao)], gdiag_sp)
         CALL onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
            & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, &
            & [(2.0_wp*dmat(iao, iao), iao=1, self%nao)], gdiag_sps)

         tmat = (0.0_wp, 0.0_wp)
         CALL shell_hadamard_add_complex(self%nsh, self%nao_sh, self%iao_sh, &
            & cache%g_mulliken, xmat, 1.0_wp, tmat)
         CALL atom_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
            & self%nao_sh, self%ish_at, self%iao_sh, cache%g_bocorr, xmat, &
            & -4.0_wp, tmat)
         CALL onsite_fx_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
            & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, xmat, &
            & 0.5_wp, tmat, adjoint_src=.TRUE.)
         CALL onsite_ri_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
            & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsri, xmat, &
            & -2.0_wp, tmat)

         umat = (0.0_wp, 0.0_wp)
         CALL shell_hadamard_add_complex(self%nsh, self%nao_sh, self%iao_sh, &
            & cache%g_mulliken, dmat, 1.0_wp, umat)
         CALL onsite_fx_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
            & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, dmat, &
            & 0.5_wp, umat)
         CALL onsite_ri_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
            & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsri, dmat, &
            & -2.0_wp, umat)

         vmat = (0.0_wp, 0.0_wp)
         CALL shell_hadamard_add_complex(self%nsh, self%nao_sh, self%iao_sh, &
            & cache%g_mulliken, density(:, :, spin), 1.0_wp, vmat)
         CALL onsite_fx_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
            & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, &
            & density(:, :, spin), 0.5_wp, vmat)
         CALL onsite_ri_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
            & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsri, &
            & density(:, :, spin), -2.0_wp, vmat)

         t2mat = tmat + 0.5_wp*MATMUL(overlap, vmat)
         DO iao = 1, self%nao
            t2mat(:, iao) = t2mat(:, iao) &
               & + 0.25_wp*gdiag_p(iao)*overlap(:, iao)
         END DO
         amat = umat + MATMUL(t2mat, overlap)
         DO iao = 1, self%nao
            amat(:, iao) = amat(:, iao) &
               & + 0.5_wp*gdiag_sp(iao)*overlap(:, iao)
            amat(iao, iao) = amat(iao, iao) + 0.25_wp*gdiag_sps(iao)
         END DO

         ! Reverse sweep, seeded by E=-(spin_factor/4)*Re<P,A>.
         bamat = -0.25_wp*spin_factor*density(:, :, spin)
         bmat = bamat  ! adjoint of U
         btmp = MATMUL(bamat, CONJG(TRANSPOSE(overlap)))  ! adjoint of T2
         overlap_grad = overlap_grad + MATMUL(CONJG(TRANSPOSE(t2mat)), bamat)
         bdiag_sp = (0.0_wp, 0.0_wp)
         bdiag_sps = (0.0_wp, 0.0_wp)
         DO iao = 1, self%nao
            overlap_grad(:, iao) = overlap_grad(:, iao) &
               & + 0.5_wp*CONJG(gdiag_sp(iao))*bamat(:, iao)
            bdiag_sp(iao) = 0.5_wp*DOT_PRODUCT(overlap(:, iao), bamat(:, iao))
            bdiag_sps(iao) = 0.25_wp*bamat(iao, iao)
         END DO

         ! T2 = T + 1/4 S diag(gdiagP) + 1/2 S V.
         tmat = btmp  ! adjoint of T
         bdiag_p = (0.0_wp, 0.0_wp)
         DO iao = 1, self%nao
            overlap_grad(:, iao) = overlap_grad(:, iao) &
               & + 0.25_wp*CONJG(gdiag_p(iao))*btmp(:, iao)
            bdiag_p(iao) = 0.25_wp*DOT_PRODUCT(overlap(:, iao), btmp(:, iao))
         END DO
         overlap_grad = overlap_grad &
            & + 0.5_wp*MATMUL(btmp, CONJG(TRANSPOSE(vmat)))
         vmat = 0.5_wp*MATMUL(CONJG(TRANSPOSE(overlap)), btmp) ! adjoint of V

         ! U response: U = Gm(D) + 1/2 Gfx(D) - 2 Gri(D).
         dmat = (0.0_wp, 0.0_wp)
         CALL shell_hadamard_add_complex(self%nsh, self%nao_sh, self%iao_sh, &
            & cache%g_mulliken, bmat, 1.0_wp, dmat)
         CALL onsite_fx_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
            & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, bmat, &
            & 0.5_wp, dmat)
         CALL onsite_ri_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
            & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsri, bmat, &
            & -2.0_wp, dmat)
         CALL shell_parameter_gradient_complex(self%nsh, self%nao_sh, &
            & self%iao_sh, bmat, 0.5_wp*MATMUL(xmat, overlap), 1.0_wp, &
            & mulliken_raw)

         ! T response, including the adjoint of Gfx(X^H).
         bamat = (0.0_wp, 0.0_wp) ! adjoint of X
         CALL shell_hadamard_add_complex(self%nsh, self%nao_sh, self%iao_sh, &
            & cache%g_mulliken, tmat, 1.0_wp, bamat)
         CALL atom_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
            & self%nao_sh, self%ish_at, self%iao_sh, cache%g_bocorr, tmat, &
            & -4.0_wp, bamat)
         CALL onsite_fx_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
            & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, tmat, &
            & 0.5_wp, bamat, adjoint_src=.TRUE.)
         CALL onsite_ri_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
            & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsri, tmat, &
            & -2.0_wp, bamat)
         CALL shell_parameter_gradient_complex(self%nsh, self%nao_sh, &
            & self%iao_sh, tmat, xmat, 1.0_wp, mulliken_raw)
         CALL atom_parameter_gradient_complex(mol%nat, mol%id, self%nsh_id, &
            & self%nao_sh, self%ish_at, self%iao_sh, tmat, xmat, -4.0_wp, &
            & bocorr_raw)

         ! V response contributes only to the Mulliken geometry kernel.
         CALL shell_parameter_gradient_complex(self%nsh, self%nao_sh, &
            & self%iao_sh, vmat, density(:, :, spin), 1.0_wp, mulliken_raw)

         ! Reverse the three onsite diagonal maps.
         CALL onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
            & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, &
            & bdiag_sp, bdiag)
         DO iao = 1, self%nao
            bamat(iao, iao) = bamat(iao, iao) + bdiag(iao)
         END DO
         CALL onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
            & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, &
            & bdiag_sps, bdiag)
         DO iao = 1, self%nao
            bamat(iao, :) = bamat(iao, :) + bdiag(iao)*overlap(iao, :)
            overlap_grad(:, iao) = overlap_grad(:, iao) &
               & + bdiag(iao)*CONJG(xmat(iao, :))
         END DO

         ! D = 1/2 X S and X = S P.
         bamat = bamat + 0.5_wp*MATMUL(dmat, CONJG(TRANSPOSE(overlap)))
         overlap_grad = overlap_grad &
            & + 0.5_wp*MATMUL(CONJG(TRANSPOSE(xmat)), dmat)
         overlap_grad = overlap_grad &
            & + MATMUL(bamat, CONJG(TRANSPOSE(density(:, :, spin))))
      END DO

      ! Only Hermitian overlap variations are physical at a full k point.
      overlap_grad = 0.5_wp*(overlap_grad + CONJG(TRANSPOSE(overlap_grad)))
      CALL symmetrize_shell_parameter_gradient(self%nsh, self%nao_sh, &
         & self%iao_sh, mulliken_raw, mulliken_grad)
      CALL symmetrize_atom_parameter_gradient(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, bocorr_raw, bocorr_grad)

   END SUBROUTINE get_KGrad_kpoint

! **************************************************************************************************
!> \brief ...
!> \param self ...
!> \param mol ...
!> \param density ...
!> \param overlap ...
!> \param g_mulliken ...
!> \param g_bocorr ...
!> \param g_onsfx ...
!> \param g_onsri ...
!> \param fock ...
! **************************************************************************************************
   SUBROUTINE build_KFock_complex(self, mol, density, overlap, g_mulliken, &
      & g_bocorr, g_onsfx, g_onsri, fock)
      class(exchange_fock), INTENT(in) :: self
      TYPE(structure_type), INTENT(in) :: mol
      COMPLEX(wp), INTENT(in) :: density(:, :, :), overlap(:, :)
      REAL(wp), INTENT(in) :: g_mulliken(:, :), g_bocorr(:, :)
      REAL(wp), INTENT(in) :: g_onsfx(:, :, :), g_onsri(:, :)
      COMPLEX(wp), INTENT(out) :: fock(:, :, :)

      INTEGER :: iao, ii, jj, spin
      REAL(wp) :: spin_factor
      COMPLEX(wp) :: tmp
      COMPLEX(wp), ALLOCATABLE :: diagP(:), diagSP(:), gdiagP(:), &
         & gdiagSP(:), gdiagSPS(:), tmpA(:, :), tmpB(:, :), tmpSPS(:)

      ALLOCATE (tmpA(self%nao, self%nao), tmpB(self%nao, self%nao), &
         & diagP(self%nao), diagSP(self%nao), gdiagP(self%nao), &
         & gdiagSP(self%nao), tmpSPS(self%nao), gdiagSPS(self%nao), &
         & source=(0.0_wp, 0.0_wp))

      spin_factor = 0.5_wp
      IF (SIZE(density, 3) > 1) spin_factor = 1.0_wp
      fock = (0.0_wp, 0.0_wp)

      DO spin = 1, SIZE(density, 3)
         ! A = S P
         tmpA = MATMUL(overlap, density(:, :, spin))

         DO iao = 1, self%nao
            diagP(iao) = density(iao, iao, spin)
            diagSP(iao) = tmpA(iao, iao)
            tmpSPS(iao) = SUM(tmpA(iao, :)*CONJG(overlap(iao, :)))
         END DO
         CALL onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
            & self%nao_sh, self%ish_at, self%iao_sh, g_onsfx, diagP, gdiagP)

         CALL onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
            & self%nao_sh, self%ish_at, self%iao_sh, g_onsfx, diagSP, gdiagSP)
         CALL onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
            & self%nao_sh, self%ish_at, self%iao_sh, g_onsfx, tmpSPS, gdiagSPS)

         ! Initial 1/2 S P S intermediate.
         fock(:, :, spin) = 0.5_wp*MATMUL(tmpA, overlap)

         ! X = g*(S P), including the adjoint P S onsite contribution.
         tmpB = (0.0_wp, 0.0_wp)
         CALL shell_hadamard_add_complex(self%nsh, self%nao_sh, self%iao_sh, &
            & g_mulliken, tmpA, 1.0_wp, tmpB)
         CALL atom_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
            & self%nao_sh, self%ish_at, self%iao_sh, g_bocorr, tmpA, &
            & -4.0_wp, tmpB)
         CALL onsite_fx_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
            & self%nao_sh, self%ish_at, self%iao_sh, g_onsfx, tmpA, &
            & 0.5_wp, tmpB, adjoint_src=.TRUE.)
         CALL onsite_ri_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
            & self%nao_sh, self%ish_at, self%iao_sh, g_onsri, tmpA, &
            & -2.0_wp, tmpB)
         tmpA = tmpB

         ! g*(1/2 S P S)
         tmpB = (0.0_wp, 0.0_wp)
         CALL onsite_fx_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
            & self%nao_sh, self%ish_at, self%iao_sh, g_onsfx, &
            & fock(:, :, spin), 0.5_wp, tmpB)
         CALL onsite_ri_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
            & self%nao_sh, self%ish_at, self%iao_sh, g_onsri, &
            & fock(:, :, spin), -2.0_wp, tmpB)
         CALL shell_hadamard_add_complex(self%nsh, self%nao_sh, self%iao_sh, &
            & g_mulliken, fock(:, :, spin), 1.0_wp, tmpB)
         fock(:, :, spin) = tmpB

         ! g*P
         tmpB = (0.0_wp, 0.0_wp)
         CALL shell_hadamard_add_complex(self%nsh, self%nao_sh, self%iao_sh, &
            & g_mulliken, density(:, :, spin), 1.0_wp, tmpB)
         CALL onsite_fx_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
            & self%nao_sh, self%ish_at, self%iao_sh, g_onsfx, &
            & density(:, :, spin), 0.5_wp, tmpB)
         CALL onsite_ri_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
            & self%nao_sh, self%ish_at, self%iao_sh, g_onsri, &
            & density(:, :, spin), -2.0_wp, tmpB)

         DO iao = 1, self%nao
            tmpA(:, iao) = tmpA(:, iao) &
               & + 0.25_wp*gdiagP(iao)*overlap(:, iao)
         END DO
         tmpA = tmpA + 0.5_wp*MATMUL(overlap, tmpB)
         fock(:, :, spin) = fock(:, :, spin) + MATMUL(tmpA, overlap)

         DO iao = 1, self%nao
            fock(:, iao, spin) = fock(:, iao, spin) &
               & + 0.5_wp*gdiagSP(iao)*overlap(:, iao)
            fock(iao, iao, spin) = fock(iao, iao, spin) &
               & + 0.25_wp*gdiagSPS(iao)
         END DO

         ! Hermitian projection and exchange sign/spin scaling.
         DO ii = 1, self%nao
            fock(ii, ii, spin) = CMPLX(-0.5_wp*spin_factor &
                                      & *REAL(fock(ii, ii, spin), wp), 0.0_wp, wp)
            DO jj = 1, ii - 1
               tmp = -0.25_wp*spin_factor*(fock(jj, ii, spin) &
                  & + CONJG(fock(ii, jj, spin)))
               fock(jj, ii, spin) = tmp
               fock(ii, jj, spin) = CONJG(tmp)
            END DO
         END DO
      END DO

   END SUBROUTINE build_KFock_complex

! **************************************************************************************************
!> \brief ...
!> \param nsh ...
!> \param nao_sh ...
!> \param iao_sh ...
!> \param raw ...
!> \param gradient ...
! **************************************************************************************************
   SUBROUTINE symmetrize_shell_parameter_gradient(nsh, nao_sh, iao_sh, raw, &
      & gradient)
      INTEGER, INTENT(in)                                :: nsh, nao_sh(:), iao_sh(:)
      REAL(wp), INTENT(in)                               :: raw(:, :)
      REAL(wp), INTENT(out)                              :: gradient(:, :)

      INTEGER                                            :: iao, ii, ish, jao, jj, jsh, ni, nj
      REAL(wp)                                           :: tmp

      gradient = raw
      DO ish = 1, nsh
         ii = iao_sh(ish)
         ni = nao_sh(ish)
         DO jsh = 1, ish - 1
            jj = iao_sh(jsh)
            nj = nao_sh(jsh)
            DO iao = 1, ni
               DO jao = 1, nj
                  tmp = raw(jj + jao, ii + iao) + raw(ii + iao, jj + jao)
                  gradient(jj + jao, ii + iao) = tmp
                  gradient(ii + iao, jj + jao) = tmp
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE symmetrize_shell_parameter_gradient

! **************************************************************************************************
!> \brief ...
!> \param nat ...
!> \param id ...
!> \param nsh_id ...
!> \param nao_sh ...
!> \param ish_at ...
!> \param iao_sh ...
!> \param raw ...
!> \param gradient ...
! **************************************************************************************************
   SUBROUTINE symmetrize_atom_parameter_gradient(nat, id, nsh_id, nao_sh, ish_at, &
      & iao_sh, raw, gradient)
      INTEGER, INTENT(in)                                :: nat, id(:), nsh_id(:), nao_sh(:), &
                                                            ish_at(:), iao_sh(:)
      REAL(wp), INTENT(in)                               :: raw(:, :)
      REAL(wp), INTENT(out)                              :: gradient(:, :)

      INTEGER                                            :: iao, iat, ii, is, ish, izp, jao, jat, &
                                                            jj, js, jsh, jzp, ni, nj
      REAL(wp)                                           :: tmp

      gradient = raw
      DO iat = 1, nat
         izp = id(iat)
         is = ish_at(iat)
         DO jat = 1, iat - 1
            jzp = id(jat)
            js = ish_at(jat)
            DO ish = 1, nsh_id(izp)
               ii = iao_sh(is + ish)
               ni = nao_sh(is + ish)
               DO jsh = 1, nsh_id(jzp)
                  jj = iao_sh(js + jsh)
                  nj = nao_sh(js + jsh)
                  DO iao = 1, ni
                     DO jao = 1, nj
                        tmp = raw(jj + jao, ii + iao) + raw(ii + iao, jj + jao)
                        gradient(jj + jao, ii + iao) = tmp
                        gradient(ii + iao, jj + jao) = tmp
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE symmetrize_atom_parameter_gradient

! **************************************************************************************************
!> \brief ...
!> \param self ...
!> \param mol ...
!> \param cache ...
!> \param mulliken_grad ...
!> \param gradient ...
!> \param sigma ...
!> \note Adapt direct symmetric-kernel adjoints to the original library convention.
!> \note Mulliken shell diagonals have a factor -2; other blocks a factor -1.
!> \note Bond-order atom diagonals have a factor 2; other blocks are unchanged.
!> \note Physical radial/WS derivative evaluation remains in the public library methods.
! **************************************************************************************************
   SUBROUTINE get_mulliken_derivs_direct(self, mol, cache, mulliken_grad, &
      & gradient, sigma)
      class(exchange_fock), INTENT(in) :: self
      TYPE(structure_type), INTENT(in) :: mol
      TYPE(exchange_cache), INTENT(inout) :: cache
      REAL(wp), contiguous, INTENT(in) :: mulliken_grad(:, :)
      REAL(wp), contiguous, INTENT(inout) :: gradient(:, :)
      REAL(wp), contiguous, INTENT(inout) :: sigma(:, :)

      INTEGER :: ii, ish, ni
      REAL(wp), ALLOCATABLE :: legacy_grad(:, :)

      ALLOCATE (legacy_grad(self%nao, self%nao), source=-mulliken_grad)
      DO ish = 1, self%nsh
         ii = self%iao_sh(ish)
         ni = self%nao_sh(ish)
         legacy_grad(ii + 1:ii + ni, ii + 1:ii + ni) = &
            & -2.0_wp*mulliken_grad(ii + 1:ii + ni, ii + 1:ii + ni)
      END DO
      CALL self%get_mulliken_derivs(mol, cache, legacy_grad, gradient, sigma)
   END SUBROUTINE get_mulliken_derivs_direct

! **************************************************************************************************
!> \brief ...
!> \param self ...
!> \param mol ...
!> \param cache ...
!> \param bocorr_grad ...
!> \param gradient ...
!> \param sigma ...
! **************************************************************************************************
   SUBROUTINE get_bocorr_derivs_direct(self, mol, cache, bocorr_grad, gradient, &
      & sigma)
      class(exchange_fock), INTENT(in) :: self
      TYPE(structure_type), INTENT(in) :: mol
      TYPE(exchange_cache), INTENT(inout) :: cache
      REAL(wp), contiguous, INTENT(in) :: bocorr_grad(:, :)
      REAL(wp), contiguous, INTENT(inout) :: gradient(:, :)
      REAL(wp), contiguous, INTENT(inout) :: sigma(:, :)

      INTEGER :: iat, ii, is, izp, ni
      REAL(wp), ALLOCATABLE :: legacy_grad(:, :)

      ALLOCATE (legacy_grad(self%nao, self%nao), source=bocorr_grad)
      DO iat = 1, mol%nat
         izp = mol%id(iat)
         is = self%ish_at(iat)
         ii = self%iao_sh(is + 1)
         ni = SUM(self%nao_sh(is + 1:is + self%nsh_id(izp)))
         legacy_grad(ii + 1:ii + ni, ii + 1:ii + ni) = &
            & 2.0_wp*bocorr_grad(ii + 1:ii + ni, ii + 1:ii + ni)
      END DO
      CALL self%get_bocorr_derivs(mol, cache, legacy_grad, gradient, sigma)
   END SUBROUTINE get_bocorr_derivs_direct

END MODULE tblite_exchange_contractions
#endif
