! Derived from save_tblite's CP2K compatibility helpers (25456f4adea68cbb98a4ac885983cdfe567288d8).
! SPDX-Identifier: LGPL-3.0-or-later

!> CP2K-owned model orchestration using the provider's public mathematical kernels.
!> Preserve native mixer scheduling, basis response and dense reference paths.
#if defined(__TBLITE_GXTB)
MODULE tblite_model_bridge
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE mctc_env,                        ONLY: error_type,&
                                              fatal_error,&
                                              wp
   USE mctc_io,                         ONLY: structure_type
   USE mctc_io_constants,               ONLY: pi
   USE multicharge,                     ONLY: get_charges
   USE tblite_acp_cache,                ONLY: acp_cache
   USE tblite_acp_type,                 ONLY: get_acp, get_acp_gradient
   USE tblite_acp_images,               ONLY: acp_image_type,&
                                              get_acp_image_gradient,&
                                              get_acp_images,&
                                              get_acp_kmesh,&
                                              get_acp_kmesh_gradient
   USE tblite_adjlist,                  ONLY: adjacency_list,&
                                              new_adjacency_list
   USE tblite_basis_cache,              ONLY: basis_cache
   USE tblite_basis_type,               ONLY: basis_type,&
                                              get_cutoff,&
                                              get_pair_cutoff
   USE tblite_blas,                     ONLY: gemv
   USE tblite_container_cache,          ONLY: container_cache
   USE tblite_cutoff,                   ONLY: get_lattice_points
   USE tblite_exchange_cache,           ONLY: exchange_cache
   USE tblite_exchange_contractions,    ONLY: get_mulliken_derivs_direct, get_bocorr_derivs_direct, &
                                              consumer_get_kfock_kpoint => get_kfock_kpoint, &
                                              consumer_get_kgrad_kpoint => get_kgrad_kpoint, &
                                              consumer_get_kfock_kmesh => get_kfock_kmesh, &
                                              consumer_get_kgrad_kmesh => get_kgrad_kmesh
   USE tblite_stream_plan,              ONLY: cp2k_stream_plan_type
   USE tblite_exchange_fock,            ONLY: exchange_bvk_kernel,&
                                              exchange_fock
   USE tblite_integral_native_integrals,ONLY: multipole_cgto,&
                                              multipole_grad_cgto
   USE tblite_integral_type,            ONLY: integral_type,&
                                              new_integral
   USE tblite_scf_info,                 ONLY: atom_resolved,&
                                              scf_info,&
                                              shell_resolved
   USE tblite_scf_mixer_cache,          ONLY: mixer_cache_container
   USE tblite_scf_mixer_diis,           ONLY: diis_input,&
                                              diis_mixer,&
                                              new_diis
   USE tblite_scf_mixer_input,          ONLY: mixer_mode
   USE tblite_scf_mixer_simple,         ONLY: new_simple,&
                                              simple_input,&
                                              simple_mixer
   USE tblite_scf_potential,            ONLY: add_pot_to_h1,&
                                              potential_type
   USE tblite_wavefunction_spin,        ONLY: magnet_to_updown,&
                                              updown_to_magnet
   USE tblite_wavefunction_type,        ONLY: new_wavefunction,&
                                              wavefunction_type
   USE tblite_xtb_calculator,           ONLY: xtb_calculator
   USE tblite_xtb_h0,                   ONLY: get_anisotropy,&
                                              get_anisotropy_gradient,&
                                              get_hamiltonian,&
                                              get_hamiltonian_gradient,&
                                              hamiltonian_image_type

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: cp2k_get_mixer_dimension, cp2k_get_scc_mixer_dimension
   PUBLIC :: cp2k_get_acp_cutoff
   PUBLIC :: cp2k_potential_mixer_type, new_cp2k_potential_mixer
   PUBLIC :: cp2k_set_scf_options
   PUBLIC :: cp2k_update_basis, cp2k_get_cgto
   PUBLIC :: cp2k_multipole_cgto, cp2k_multipole_grad_cgto
   PUBLIC :: cp2k_build_gamma_h0, cp2k_build_image_h0, cp2k_build_image_acp, &
      & cp2k_build_acp_projectors
   PUBLIC :: cp2k_prepare_exchange, cp2k_exchange_kpoint, cp2k_exchange_kmesh, &
      & cp2k_exchange_kmesh_gradient, cp2k_exchange_kpoint_gradient, &
      & cp2k_gamma_gradient
   PUBLIC :: cp2k_image_h0_gradient
   PUBLIC :: cp2k_acp_image_gradient, cp2k_acp_kmesh, &
      & cp2k_acp_kmesh_gradient, cp2k_apply_basis_gradient
   PUBLIC :: acp_image_type, hamiltonian_image_type

   !> CP2K-facing wrapper for g-xTB's native potential mixer.
   !>
   !> CP2K owns the packing of all real Fock degrees of freedom into one flat
   !> vector.  This wrapper deliberately owns only the mixer state and mirrors
   !> g-xTB's iterator schedule: simple mixing with damping 0.2 for the first
   !> three Fock builds, DIIS precollection on the third build, and DIIS from
   !> the fourth build onward.
   TYPE, PUBLIC :: cp2k_potential_mixer_type
      PRIVATE
      TYPE(simple_mixer) :: simple
      TYPE(diis_mixer) :: diis
      TYPE(mixer_cache_container) :: simple_cache
      TYPE(mixer_cache_container) :: diis_cache
      INTEGER :: ndim = 0
      INTEGER :: iteration = 0
      REAL(wp) :: residual = 0.0_wp
      LOGICAL :: initialized = .FALSE.
   CONTAINS
      PROCEDURE, PUBLIC :: reset => cp2k_potential_mixer_reset
      PROCEDURE, PUBLIC :: advance => cp2k_potential_mixer_advance
      PROCEDURE, PUBLIC :: get_error => cp2k_potential_mixer_get_error
   END TYPE cp2k_potential_mixer_type

CONTAINS

!> Return a real-space cutoff that covers both the valence basis and the
!> auxiliary ACP projectors used by g-xTB.
   pure FUNCTION cp2k_get_acp_cutoff(calc, accuracy) RESULT(cutoff)
      TYPE(xtb_calculator), INTENT(in)                   :: calc
      REAL(wp), INTENT(in)                               :: accuracy
      REAL(wp)                                           :: cutoff

      cutoff = get_cutoff(calc%bas, accuracy)
      IF (ALLOCATED(calc%acp)) THEN
         cutoff = MAX(cutoff, get_pair_cutoff(calc%bas, calc%acp%auxbas, accuracy))
      END IF
   END FUNCTION cp2k_get_acp_cutoff

!> Construct a CP2K-facing g-xTB potential mixer.
   SUBROUTINE new_cp2k_potential_mixer(self, ndim, error)
      TYPE(cp2k_potential_mixer_type), INTENT(out)       :: self
      INTEGER, INTENT(in)                                :: ndim
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      CALL self%reset(ndim, error)
   END SUBROUTINE new_cp2k_potential_mixer

!> Reset a CP2K-facing g-xTB potential mixer for a flat Fock vector.
   SUBROUTINE cp2k_potential_mixer_reset(self, ndim, error)
      class(cp2k_potential_mixer_type), INTENT(inout) :: self
      INTEGER, INTENT(in) :: ndim
      TYPE(error_type), ALLOCATABLE, INTENT(out) :: error

      IF (ALLOCATED(self%simple_cache%raw)) THEN
         CALL self%simple%cleanup(self%simple_cache)
         DEALLOCATE (self%simple_cache%raw)
      END IF
      IF (ALLOCATED(self%diis_cache%raw)) THEN
         CALL self%diis%cleanup(self%diis_cache)
         DEALLOCATE (self%diis_cache%raw)
      END IF

      self%ndim = 0
      self%iteration = 0
      self%residual = 0.0_wp
      self%initialized = .FALSE.

      IF (ndim <= 0) THEN
         CALL fatal_error(error, "CP2K potential mixer dimension must be positive")
         RETURN
      END IF

      CALL new_simple(self%simple, simple_input(mode=mixer_mode%potential, &
         & start=1, damp=0.2_wp))
      CALL new_diis(self%diis, diis_input(mode=mixer_mode%potential, start=4, &
         & precollect=1, memory=7, damp=1.0_wp, output_fraction=1.0_wp), error)
      IF (ALLOCATED(error)) RETURN

      self%ndim = ndim
      self%initialized = .TRUE.
   END SUBROUTINE cp2k_potential_mixer_reset

!> Advance g-xTB's native potential mixer from a raw to a mixed Fock vector.
   SUBROUTINE cp2k_potential_mixer_advance(self, raw_fock, mixed_fock, error)
      class(cp2k_potential_mixer_type), INTENT(inout) :: self
      REAL(wp), INTENT(in) :: raw_fock(:)
      REAL(wp), INTENT(out) :: mixed_fock(:)
      TYPE(error_type), ALLOCATABLE, INTENT(out) :: error

      INTEGER :: next_iteration

      mixed_fock = 0.0_wp
      IF (.NOT. self%initialized) THEN
         CALL fatal_error(error, "CP2K potential mixer must be reset before use")
         RETURN
      END IF
      IF (SIZE(raw_fock) /= self%ndim .OR. SIZE(mixed_fock) /= self%ndim) THEN
         CALL fatal_error(error, "CP2K potential mixer Fock vector dimension mismatch")
         RETURN
      END IF
      IF (.NOT. ALL(ieee_is_finite(raw_fock))) THEN
         CALL fatal_error(error, "CP2K potential mixer received a non-finite Fock vector")
         RETURN
      END IF

      next_iteration = self%iteration + 1
      SELECT CASE (next_iteration)
      CASE (1)
         CALL self%simple%update(self%simple_cache, self%ndim)
         self%simple_cache%raw%iset = 0
         CALL self%simple%set(self%simple_cache, raw_fock)
         self%simple_cache%raw%initialized = .TRUE.
         mixed_fock = raw_fock
         self%residual = 0.0_wp

      CASE (2)
         ! The DIIS cache is prepared one step before its first precollection.
         CALL self%diis%update(self%diis_cache, self%ndim)

         self%simple_cache%raw%idif = 0
         CALL self%simple%diff(self%simple_cache, raw_fock)
         self%residual = self%simple%get_residual_norm(self%simple_cache)
         CALL self%simple%collect(self%simple_cache)
         CALL self%simple%next(self%simple_cache, error)
         IF (ALLOCATED(error)) RETURN
         self%simple_cache%raw%iget = 0
         CALL self%simple%get(self%simple_cache, mixed_fock)

         CALL set_simple_input(self, mixed_fock)
         CALL set_diis_input(self, mixed_fock)

      CASE (3)
         self%simple_cache%raw%idif = 0
         CALL self%simple%diff(self%simple_cache, raw_fock)
         self%diis_cache%raw%idif = 0
         CALL self%diis%diff(self%diis_cache, raw_fock)
         self%residual = self%simple%get_residual_norm(self%simple_cache)

         CALL self%simple%collect(self%simple_cache)
         CALL self%diis%collect(self%diis_cache)
         CALL self%simple%next(self%simple_cache, error)
         IF (ALLOCATED(error)) RETURN
         self%simple_cache%raw%iget = 0
         CALL self%simple%get(self%simple_cache, mixed_fock)

         CALL set_simple_input(self, mixed_fock)
         CALL set_diis_input(self, mixed_fock)

      CASE default
         IF (next_iteration == 4) THEN
            ! This is the iterator's simple-to-DIIS transition point.
            CALL self%simple%cleanup(self%simple_cache)
            DEALLOCATE (self%simple_cache%raw)
            CALL self%diis%update(self%diis_cache, self%ndim)
         END IF

         self%diis_cache%raw%idif = 0
         CALL self%diis%diff(self%diis_cache, raw_fock)
         self%residual = self%diis%get_residual_norm(self%diis_cache)
         CALL self%diis%collect(self%diis_cache)
         CALL self%diis%next(self%diis_cache, error)
         IF (ALLOCATED(error)) RETURN
         self%diis_cache%raw%iget = 0
         CALL self%diis%get(self%diis_cache, mixed_fock)
         CALL set_diis_input(self, mixed_fock)
      END SELECT

      IF (.NOT. ieee_is_finite(self%residual)) THEN
         CALL fatal_error(error, "CP2K potential mixer produced a non-finite residual")
         RETURN
      END IF
      IF (.NOT. ALL(ieee_is_finite(mixed_fock))) THEN
         CALL fatal_error(error, "CP2K potential mixer produced a non-finite Fock vector")
         RETURN
      END IF
      self%iteration = next_iteration
   END SUBROUTINE cp2k_potential_mixer_advance

!> Return the RMS residual from the most recent raw Fock vector.
   pure FUNCTION cp2k_potential_mixer_get_error(self) RESULT(residual)
      class(cp2k_potential_mixer_type), INTENT(in) :: self
      REAL(wp) :: residual

      residual = self%residual
   END FUNCTION cp2k_potential_mixer_get_error

!> Store the current extrapolated vector in the simple mixer cache.
   SUBROUTINE set_simple_input(self, fock)
      class(cp2k_potential_mixer_type), INTENT(inout) :: self
      REAL(wp), INTENT(in) :: fock(:)

      self%simple_cache%raw%iset = 0
      CALL self%simple%set(self%simple_cache, fock)
      self%simple_cache%raw%initialized = .TRUE.
   END SUBROUTINE set_simple_input

!> Store the current extrapolated vector in the DIIS mixer cache.
   SUBROUTINE set_diis_input(self, fock)
      class(cp2k_potential_mixer_type), INTENT(inout) :: self
      REAL(wp), INTENT(in) :: fock(:)

      self%diis_cache%raw%iset = 0
      CALL self%diis%set(self%diis_cache, fock)
      self%diis_cache%raw%initialized = .TRUE.
   END SUBROUTINE set_diis_input

!> Return the CP2K-owned SCC mixer dimension without the spin multiplicity.
!>
!> The orbital density is deliberately not included here, even when a
!> contribution advertises an orbital-resolved density potential.  CP2K
!> diagonalizes and stores that density itself; only charges and multipoles
!> are exchanged with its SCC mixer through the compatibility interface.
   pure FUNCTION cp2k_get_scc_mixer_dimension(mol, bas, info) RESULT(ndim)
      TYPE(structure_type), INTENT(in)                   :: mol

      class(basis_type), INTENT(in) :: bas
      TYPE(scf_info), INTENT(in) :: info
      INTEGER :: ndim

      ndim = 0
      SELECT CASE (info%charge)
      CASE (atom_resolved)
         ndim = ndim + mol%nat
      CASE (shell_resolved)
         ndim = ndim + bas%nsh
      END SELECT
      IF (info%dipole == atom_resolved) ndim = ndim + 3*mol%nat
      IF (info%quadrupole == atom_resolved) ndim = ndim + 6*mol%nat
   END FUNCTION cp2k_get_scc_mixer_dimension

!> Backward-compatible name for the CP2K SCC mixer dimension.
   pure FUNCTION cp2k_get_mixer_dimension(mol, bas, info) RESULT(ndim)
      TYPE(structure_type), INTENT(in)                   :: mol

      class(basis_type), INTENT(in) :: bas
      TYPE(scf_info), INTENT(in) :: info
      INTEGER :: ndim

      ndim = cp2k_get_scc_mixer_dimension(mol, bas, info)
   END FUNCTION cp2k_get_mixer_dimension

!> Apply SCC settings using the iterator API used by save_tblite.
   SUBROUTINE cp2k_set_scf_options(calc, max_iter, damping)
      TYPE(xtb_calculator), INTENT(inout)                :: calc
      INTEGER, INTENT(in)                                :: max_iter
      REAL(wp), INTENT(in)                               :: damping

      INTEGER                                            :: imix

      calc%iterator%max_iter = max_iter
      IF (ALLOCATED(calc%iterator%mixers)) THEN
         DO imix = 1, SIZE(calc%iterator%mixers)
            calc%iterator%mixers(imix)%raw%damp = damping
         END DO
      END IF
   END SUBROUTINE cp2k_set_scf_options

!> Construct the immutable auxiliary charges and cache the effective basis.
   SUBROUTINE cp2k_update_basis(calc, mol, cache, error, grad, wfn_aux)
      TYPE(xtb_calculator), INTENT(inout)                :: calc
      TYPE(structure_type), INTENT(in)                   :: mol
      TYPE(basis_cache), INTENT(inout)                   :: cache
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error
      LOGICAL, INTENT(in), OPTIONAL                      :: grad
      TYPE(wavefunction_type), INTENT(inout), OPTIONAL   :: wfn_aux

      LOGICAL                                            :: do_grad
      TYPE(wavefunction_type)                            :: local_wfn_aux

      do_grad = .FALSE.
      IF (PRESENT(grad)) do_grad = grad

      IF (ALLOCATED(calc%charge_model)) THEN
         IF (PRESENT(wfn_aux)) THEN
            CALL new_wavefunction(wfn_aux, mol%nat, calc%bas%nsh, 0, 1, 0.0_wp, do_grad)
            IF (do_grad) THEN
               CALL get_charges(calc%charge_model, mol, error, wfn_aux%qat(:, 1), &
                  & dqdr=wfn_aux%dqatdr(:, :, :, 1), dqdL=wfn_aux%dqatdL(:, :, :, 1))
            ELSE
               CALL get_charges(calc%charge_model, mol, error, wfn_aux%qat(:, 1))
            END IF
            IF (ALLOCATED(error)) RETURN
            CALL calc%bas%update(mol, cache, do_grad, wfn_aux)
         ELSE
            CALL new_wavefunction(local_wfn_aux, mol%nat, calc%bas%nsh, 0, 1, 0.0_wp, do_grad)
            IF (do_grad) THEN
               CALL get_charges(calc%charge_model, mol, error, local_wfn_aux%qat(:, 1), &
                  & dqdr=local_wfn_aux%dqatdr(:, :, :, 1), &
                  & dqdL=local_wfn_aux%dqatdL(:, :, :, 1))
            ELSE
               CALL get_charges(calc%charge_model, mol, error, local_wfn_aux%qat(:, 1))
            END IF
            IF (ALLOCATED(error)) RETURN
            CALL calc%bas%update(mol, cache, do_grad, local_wfn_aux)
         END IF
      ELSE
         CALL calc%bas%update(mol, cache, do_grad)
      END IF
   END SUBROUTINE cp2k_update_basis

!> Build the complete Gamma-point overlap and one-electron g-xTB Hamiltonian.
!>
!> The dense matrices are deliberately returned in tblite's atom/AO ordering.
!> CP2K can scatter them into its block-sparse matrices while retaining the
!> element-wise basis metadata used for block dimensions and cutoffs.
   SUBROUTINE cp2k_build_gamma_h0(calc, mol, bcache, acache, accuracy, selfenergy, &
      & overlap, dipole, quadrupole, hamiltonian)
      TYPE(xtb_calculator), INTENT(inout)                :: calc
      TYPE(structure_type), INTENT(in)                   :: mol
      TYPE(basis_cache), INTENT(in)                      :: bcache
      TYPE(acp_cache), INTENT(inout)                     :: acache
      REAL(wp), INTENT(in)                               :: accuracy, selfenergy(:)
      REAL(wp), INTENT(out)                              :: overlap(:, :), dipole(:, :, :), &
                                                            quadrupole(:, :, :), hamiltonian(:, :)

      REAL(wp)                                           :: cutoff
      REAL(wp), ALLOCATABLE                              :: aniso_dip(:, :, :), lattr(:, :)
      TYPE(adjacency_list)                               :: list, list_acp

      ALLOCATE (aniso_dip(3, 3, mol%nat))
      cutoff = cp2k_get_acp_cutoff(calc, accuracy)
      CALL get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
      CALL get_anisotropy(calc%h0, mol, lattr, cutoff, aniso_dip)
      CALL new_adjacency_list(list, mol, lattr, cutoff)
      CALL get_hamiltonian(mol, lattr, list, calc%bas, bcache, calc%h0, &
         & selfenergy, aniso_dip, overlap, dipole, quadrupole, hamiltonian)

      IF (ALLOCATED(calc%acp)) THEN
         CALL calc%acp%update(mol, acache)
         ! H0 uses a half list, whereas the projector-centred ACP traversal
         ! requires every valence neighbour of each auxiliary atom.
         CALL new_adjacency_list(list_acp, mol, lattr, cutoff, complete=.TRUE.)
         CALL get_acp(mol, lattr, list_acp, calc%bas, bcache, calc%acp, acache, hamiltonian)
      END IF
   END SUBROUTINE cp2k_build_gamma_h0

!> Build image-resolved q-vSZP overlap, multipoles, and H0 matrices.
!>
!> This routine deliberately covers the one-electron H0 only.  The atomic
!> correction potential is a projector convolution and therefore requires a
!> separate image-resolved implementation before a complete k-point
!> Hamiltonian can be assembled.
   SUBROUTINE cp2k_build_image_h0(calc, mol, bcache, accuracy, selfenergy, &
      & translations, images, overlap_gamma, dipole_gamma, quadrupole_gamma, &
      & hamiltonian_gamma)
      TYPE(xtb_calculator), INTENT(inout)                :: calc
      TYPE(structure_type), INTENT(in)                   :: mol
      TYPE(basis_cache), INTENT(in)                      :: bcache
      REAL(wp), INTENT(in)                               :: accuracy, selfenergy(:), &
                                                            translations(:, :)
      TYPE(hamiltonian_image_type), INTENT(out)          :: images
      REAL(wp), INTENT(out), OPTIONAL                    :: overlap_gamma(:, :), &
                                                            dipole_gamma(:, :, :), &
                                                            quadrupole_gamma(:, :, :), &
                                                            hamiltonian_gamma(:, :)

      REAL(wp)                                           :: cutoff
      REAL(wp), ALLOCATABLE                              :: aniso_dip(:, :, :), dipole(:, :, :), &
                                                            hamiltonian(:, :), overlap(:, :), &
                                                            quadrupole(:, :, :)
      TYPE(adjacency_list)                               :: list

      ALLOCATE (aniso_dip(3, 3, mol%nat), overlap(calc%bas%nao, calc%bas%nao), &
         & dipole(3, calc%bas%nao, calc%bas%nao), &
         & quadrupole(6, calc%bas%nao, calc%bas%nao), &
         & hamiltonian(calc%bas%nao, calc%bas%nao))
      cutoff = get_cutoff(calc%bas, accuracy)
      CALL get_anisotropy(calc%h0, mol, translations, cutoff, aniso_dip)
      CALL new_adjacency_list(list, mol, translations, cutoff)
      CALL get_hamiltonian(mol, translations, list, calc%bas, bcache, calc%h0, &
         & selfenergy, aniso_dip, overlap, dipole, quadrupole, hamiltonian, images)

      IF (PRESENT(overlap_gamma)) overlap_gamma = overlap
      IF (PRESENT(dipole_gamma)) dipole_gamma = dipole
      IF (PRESENT(quadrupole_gamma)) quadrupole_gamma = quadrupole
      IF (PRESENT(hamiltonian_gamma)) hamiltonian_gamma = hamiltonian
   END SUBROUTINE cp2k_build_image_h0

!> Build the image-resolved separable ACP contribution.
!>
!> The returned Hamiltonian translations are generally a larger difference
!> set than the translations supplied for the projector integrals.  External
!> drivers must therefore scatter this set independently of the pair-local H0
!> image set.
   SUBROUTINE cp2k_build_image_acp(calc, mol, bcache, acache, translations, &
      & images, hamiltonian_gamma)
      TYPE(xtb_calculator), INTENT(inout)                :: calc
      TYPE(structure_type), INTENT(in)                   :: mol
      TYPE(basis_cache), INTENT(in)                      :: bcache
      TYPE(acp_cache), INTENT(inout)                     :: acache
      REAL(wp), INTENT(in)                               :: translations(:, :)
      TYPE(acp_image_type), INTENT(out)                  :: images
      REAL(wp), INTENT(out), OPTIONAL                    :: hamiltonian_gamma(:, :)

      CALL get_acp_images(mol, translations, calc%bas, bcache, calc%acp, acache, &
         & images)
      IF (PRESENT(hamiltonian_gamma)) THEN
         hamiltonian_gamma = SUM(images%hamiltonian, dim=3)
      END IF
   END SUBROUTINE cp2k_build_image_acp

!> Build only the compact projector-overlap images required by the direct
!> Bloch-space ACP forward and reverse sweeps.  Unlike
!> `cp2k_build_image_acp`, this deliberately does not materialize the
!> quadratic translation-difference set or H(R).
   SUBROUTINE cp2k_build_acp_projectors(calc, mol, bcache, acache, translations, &
      & images)
      TYPE(xtb_calculator), INTENT(inout)                :: calc
      TYPE(structure_type), INTENT(in)                   :: mol
      TYPE(basis_cache), INTENT(in)                      :: bcache
      TYPE(acp_cache), INTENT(inout)                     :: acache
      REAL(wp), INTENT(in)                               :: translations(:, :)
      TYPE(acp_image_type), INTENT(out)                  :: images

      CALL get_acp_images(mol, translations, calc%bas, bcache, calc%acp, acache, &
         & images, build_hamiltonian=.FALSE.)
   END SUBROUTINE cp2k_build_acp_projectors

!> Image-resolved ACP response for CP2K k-point density matrices.
!>
!> `density(:, :, R, spin)` must be in save_tblite's AO image orientation and
!> dual to the Hamiltonian returned by `cp2k_build_image_acp`, so that its
!> energy contribution is the element-wise contraction over equal images.
!> CP2K's irreducible-star weights must already be included in these real-space
!> densities.
   SUBROUTINE cp2k_acp_image_gradient(calc, mol, bcache, acache, images, &
      & density_translation, density, dEdcnbas, dEdqbas, gradient, sigma, error)
      TYPE(xtb_calculator), INTENT(in)                   :: calc
      TYPE(structure_type), INTENT(in)                   :: mol
      TYPE(basis_cache), INTENT(in)                      :: bcache
      TYPE(acp_cache), INTENT(in)                        :: acache
      TYPE(acp_image_type), INTENT(in)                   :: images
      REAL(wp), INTENT(in)                               :: density_translation(:, :), &
                                                            density(:, :, :, :)
      REAL(wp), INTENT(inout)                            :: dEdcnbas(:), dEdqbas(:), gradient(:, :), &
                                                            sigma(:, :)
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      IF (.NOT. ALLOCATED(calc%acp)) THEN
         CALL fatal_error(error, "g-xTB ACP image gradient requires an ACP model")
         RETURN
      END IF
      CALL get_acp_image_gradient(mol, calc%bas, bcache, calc%acp, acache, &
         & images, density_translation, density, dEdcnbas, dEdqbas, gradient, &
         & sigma, error)
   END SUBROUTINE cp2k_acp_image_gradient

!> Evaluate the separable ACP directly on CP2K Bloch points.
   SUBROUTINE cp2k_acp_kmesh(calc, mol, images, kfrac, hamiltonian, error)
      TYPE(xtb_calculator), INTENT(in)                   :: calc
      TYPE(structure_type), INTENT(in)                   :: mol
      TYPE(acp_image_type), INTENT(in)                   :: images
      REAL(wp), INTENT(in)                               :: kfrac(:, :)
      COMPLEX(wp), INTENT(out)                           :: hamiltonian(:, :, :)
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      IF (.NOT. ALLOCATED(calc%acp)) THEN
         CALL fatal_error(error, "g-xTB ACP container is not allocated")
         RETURN
      END IF
      CALL get_acp_kmesh(mol, images, kfrac, hamiltonian, error)
   END SUBROUTINE cp2k_acp_kmesh

!> Differentiate the complete weighted Bloch contraction of the ACP.
   SUBROUTINE cp2k_acp_kmesh_gradient(calc, mol, bcache, acache, images, &
      & kfrac, weights, density, dEdcnbas, dEdqbas, gradient, sigma, error)
      TYPE(xtb_calculator), INTENT(in)                   :: calc
      TYPE(structure_type), INTENT(in)                   :: mol
      TYPE(basis_cache), INTENT(in)                      :: bcache
      TYPE(acp_cache), INTENT(in)                        :: acache
      TYPE(acp_image_type), INTENT(in)                   :: images
      REAL(wp), INTENT(in)                               :: kfrac(:, :), weights(:)
      COMPLEX(wp), INTENT(in)                            :: density(:, :, :, :)
      REAL(wp), INTENT(inout)                            :: dEdcnbas(:), dEdqbas(:), gradient(:, :), &
                                                            sigma(:, :)
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      IF (.NOT. ALLOCATED(calc%acp)) THEN
         CALL fatal_error(error, "g-xTB ACP container is not allocated")
         RETURN
      END IF
      CALL get_acp_kmesh_gradient(mol, calc%bas, bcache, calc%acp, acache, &
         & images, kfrac, weights, density, dEdcnbas, dEdqbas, gradient, &
         & sigma, error)
   END SUBROUTINE cp2k_acp_kmesh_gradient

!> Apply a component's q-vSZP basis response and external-charge chain rule.
!>
!> Image-resolved one-electron components return their own additive
!> `dE/dCN_basis` and `dE/dq_basis`.  This helper maps exactly that component
!> onto Cartesian and strain derivatives.  Calling it separately for H0 and
!> ACP therefore remains additive and avoids applying either response twice.
   SUBROUTINE cp2k_apply_basis_gradient(calc, mol, wfn_aux, dEdcnbas, dEdqbas, &
      & gradient, sigma)
      TYPE(xtb_calculator), INTENT(in)                   :: calc
      TYPE(structure_type), INTENT(in)                   :: mol
      TYPE(wavefunction_type), INTENT(in)                :: wfn_aux
      REAL(wp), INTENT(in)                               :: dEdcnbas(:), dEdqbas(:)
      REAL(wp), INTENT(inout)                            :: gradient(:, :), sigma(:, :)

      REAL(wp), ALLOCATABLE                              :: dEdq(:)

      ALLOCATE (dEdq(mol%nat), source=0.0_wp)
      IF (calc%bas%charge_dependent) THEN
         CALL calc%bas%get_basis_gradient(mol, dEdcnbas, dEdqbas, dEdq, &
            & gradient, sigma)
      END IF
      IF (ALLOCATED(wfn_aux%dqatdr)) THEN
         CALL gemv(wfn_aux%dqatdr(:, :, :, 1), dEdq, gradient, beta=1.0_wp)
      END IF
      IF (ALLOCATED(wfn_aux%dqatdL)) THEN
         CALL gemv(wfn_aux%dqatdL(:, :, :, 1), dEdq, sigma, beta=1.0_wp)
      END IF
   END SUBROUTINE cp2k_apply_basis_gradient

!> Prepare geometry- and charge-dependent exchange kernels.
!>
!> This routine must be called once, serially, before evaluating the k-point
!> blocks of every SCC iteration.  Rebuilding here also makes changes of the
!> geometry or lattice visible to the Wigner-Seitz and bond-order kernels.
   SUBROUTINE cp2k_prepare_exchange(calc, mol, ecache, wfn, error)
      TYPE(xtb_calculator), INTENT(inout)                :: calc
      TYPE(structure_type), INTENT(in)                   :: mol
      TYPE(container_cache), INTENT(inout)               :: ecache
      TYPE(wavefunction_type), INTENT(in)                :: wfn
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      IF (.NOT. ALLOCATED(calc%exchange)) THEN
         CALL fatal_error(error, "g-xTB exchange container is not allocated")
         RETURN
      END IF

      CALL calc%exchange%update(mol, ecache)
      SELECT TYPE (cache => ecache%raw)
      TYPE is (exchange_cache)
         CALL calc%exchange%get_onsite_Kmatrix(mol, wfn, cache)
      class default
         CALL fatal_error(error, "unexpected g-xTB exchange cache type")
      END SELECT
   END SUBROUTINE cp2k_prepare_exchange

!> Evaluate the g-xTB exchange operator for one complex k-point block.
!>
!> The input P(k) must not contain a k-point or star weight.  The caller owns
!> all Brillouin-zone weights and symmetry expansion.  This routine returns
!> an unweighted Hermitian F(k), an atom/shell-resolved potential, and an
!> unweighted exchange energy.  For a reduced k mesh the shell potential must
!> be permuted with every star operation before the weighted star sum is made.
!> Call cp2k_prepare_exchange serially before entering the k-point loop; this
!> evaluator then only reads the shared cache and can be called concurrently.
   SUBROUTINE cp2k_exchange_kpoint(calc, mol, ecache, density, overlap, &
      & fock, vsh, energy, error)
      TYPE(xtb_calculator), INTENT(in)                   :: calc
      TYPE(structure_type), INTENT(in)                   :: mol
      TYPE(container_cache), INTENT(in)                  :: ecache
      COMPLEX(wp), INTENT(in)                            :: density(:, :, :), overlap(:, :)
      COMPLEX(wp), INTENT(out)                           :: fock(:, :, :)
      REAL(wp), INTENT(out)                              :: vsh(:), energy
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      IF (.NOT. ALLOCATED(calc%exchange)) THEN
         CALL fatal_error(error, "g-xTB exchange container is not allocated")
         RETURN
      END IF
      IF (.NOT. ALLOCATED(ecache%raw)) THEN
         CALL fatal_error(error, &
            & "k-point exchange cache is not prepared; call cp2k_prepare_exchange")
         RETURN
      END IF

      SELECT TYPE (cache => ecache%raw)
      TYPE is (exchange_cache)
         SELECT TYPE (exchange => calc%exchange)
         TYPE is (exchange_fock)
            CALL consumer_get_kfock_kpoint(exchange, mol, cache, density, overlap, fock, &
               & vsh, energy)
         class default
            CALL fatal_error(error, &
               & "k-point exchange requires the exchange_fock implementation")
         END SELECT
      class default
         CALL fatal_error(error, "unexpected g-xTB exchange cache type")
      END SELECT
   END SUBROUTINE cp2k_exchange_kpoint

!> Evaluate BvK-supercell-equivalent exchange on a complete regular k mesh.
!>
!> Density, overlap, and Fock contain unweighted Bloch blocks.  Energy and
!> shell potential are returned as complete Brillouin-zone contractions per
!> primitive cell and must not be weighted a second time by CP2K.
   SUBROUTINE cp2k_exchange_kmesh(calc, mol, ecache, wfn, nmesh, kfrac, weights, &
      & density, overlap, fock, vsh, energy, error)
      TYPE(xtb_calculator), INTENT(inout)                :: calc
      TYPE(structure_type), INTENT(in)                   :: mol
      TYPE(container_cache), INTENT(inout)               :: ecache
      TYPE(wavefunction_type), INTENT(in)                :: wfn
      INTEGER, INTENT(in)                                :: nmesh(3)
      REAL(wp), INTENT(in)                               :: kfrac(:, :), weights(:)
      COMPLEX(wp), INTENT(in)                            :: density(:, :, :, :), overlap(:, :, :)
      COMPLEX(wp), INTENT(out)                           :: fock(:, :, :, :)
      REAL(wp), INTENT(out)                              :: vsh(:), energy
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      COMPLEX(wp)                                        :: orthogonality
      INTEGER                                            :: icell, idim, jcell, nk
      REAL(wp)                                           :: target
      TYPE(exchange_bvk_kernel)                          :: kernel
      TYPE(cp2k_stream_plan_type)                        :: reference_plan

      nk = SIZE(weights)
      IF (ANY(nmesh < 1) .OR. PRODUCT(nmesh) /= nk) THEN
         CALL fatal_error(error, "g-xTB exchange requires a complete regular k mesh")
         RETURN
      END IF
      DO idim = 1, 3
         IF (.NOT. mol%periodic(idim) .AND. nmesh(idim) /= 1) THEN
            CALL fatal_error(error, &
               & "nonperiodic directions cannot be replicated by the g-xTB k mesh")
            RETURN
         END IF
      END DO
      IF (SIZE(kfrac, 1) /= 3 .OR. SIZE(kfrac, 2) /= nk &
         & .OR. SIZE(overlap, 3) /= nk .OR. SIZE(density, 4) /= nk &
         & .OR. SIZE(fock, 4) /= nk) THEN
         CALL fatal_error(error, "inconsistent g-xTB whole-mesh array dimensions")
         RETURN
      END IF
      IF (.NOT. ALL(ieee_is_finite(kfrac))) THEN
         CALL fatal_error(error, "g-xTB whole-mesh k points must be finite")
         RETURN
      END IF
      IF (.NOT. ALL(ieee_is_finite(weights)) .OR. &
         & ABS(SUM(weights) - 1.0_wp) > 1.0e-12_wp &
         & .OR. MAXVAL(ABS(weights - 1.0_wp/REAL(nk, wp))) > 1.0e-12_wp) THEN
         CALL fatal_error(error, "g-xTB whole-mesh exchange requires uniform weights")
         RETURN
      END IF

      CALL cp2k_prepare_exchange(calc, mol, ecache, wfn, error)
      IF (ALLOCATED(error)) RETURN
      SELECT TYPE (cache => ecache%raw)
      TYPE is (exchange_cache)
         SELECT TYPE (exchange => calc%exchange)
         TYPE is (exchange_fock)
            IF (SIZE(overlap, 1) /= exchange%nao &
               & .OR. SIZE(overlap, 2) /= exchange%nao &
               & .OR. SIZE(density, 1) /= exchange%nao &
               & .OR. SIZE(density, 2) /= exchange%nao &
               & .OR. SIZE(fock, 1) /= exchange%nao &
               & .OR. SIZE(fock, 2) /= exchange%nao &
               & .OR. (SIZE(density, 3) /= 1 .AND. SIZE(density, 3) /= 2) &
               & .OR. SIZE(density, 3) /= wfn%nspin &
               & .OR. SIZE(fock, 3) /= SIZE(density, 3) &
               & .OR. SIZE(vsh) /= exchange%nsh) THEN
               CALL fatal_error(error, &
                  & "inconsistent g-xTB whole-mesh AO, spin, or shell dimensions")
               RETURN
            END IF
            CALL exchange%get_bvk_Kmatrix(mol, nmesh, kernel)
            ! The actual fractional k points include any common MP twist.  Their
            ! phase matrix must nevertheless be unitary on the BvK representatives.
            DO icell = 1, nk
               DO jcell = 1, nk
                  orthogonality = SUM(weights*EXP(CMPLX(0.0_wp, 2.0_wp*pi* &
                     & MATMUL(TRANSPOSE(kfrac), REAL(kernel%reps(:, icell) &
                     & - kernel%reps(:, jcell), wp)), wp)))
                  TARGET = MERGE(1.0_wp, 0.0_wp, icell == jcell)
                  IF (ABS(orthogonality - TARGET) > 1.0e-10_wp) THEN
                     CALL fatal_error(error, &
                        & "g-xTB k points are not a complete Fourier mesh")
                     RETURN
                  END IF
               END DO
            END DO
            CALL consumer_get_kfock_kmesh(exchange, mol, cache, reference_plan, kernel, kfrac, weights, &
               & density, overlap, fock, vsh, energy)
         class default
            CALL fatal_error(error, &
               & "whole-mesh exchange requires the exchange_fock implementation")
         END SELECT
      class default
         CALL fatal_error(error, "unexpected g-xTB exchange cache type")
      END SELECT

   END SUBROUTINE cp2k_exchange_kmesh

!> Differentiate BvK-supercell-equivalent exchange on a complete k mesh.
!>
!> The overlap response is spin-summed, Hermitian, and unweighted:
!> ``dE = sum_k w_k Re sum_ij conjg(overlap_grad_ij(k))*dS_ij(k)``.
!> Atomic and positive homogeneous-strain derivatives already contain the
!> complete primitive-cell Brillouin-zone contraction and must not be weighted
!> again by CP2K.  The strain response has no inverse-volume factor.
   SUBROUTINE cp2k_exchange_kmesh_gradient(calc, mol, ecache, nmesh, kfrac, &
      & weights, density, overlap, overlap_grad, gradient, sigma, error)
      TYPE(xtb_calculator), INTENT(in)                   :: calc
      TYPE(structure_type), INTENT(in)                   :: mol
      TYPE(container_cache), INTENT(in)                  :: ecache
      INTEGER, INTENT(in)                                :: nmesh(3)
      REAL(wp), INTENT(in)                               :: kfrac(:, :), weights(:)
      COMPLEX(wp), INTENT(in)                            :: density(:, :, :, :), overlap(:, :, :)
      COMPLEX(wp), INTENT(out)                           :: overlap_grad(:, :, :)
      REAL(wp), INTENT(out)                              :: gradient(:, :), sigma(:, :)
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      COMPLEX(wp)                                        :: orthogonality
      INTEGER                                            :: icell, idim, jcell, nk
      REAL(wp)                                           :: target
      REAL(wp), ALLOCATABLE                              :: bocorr_grad_r(:, :, :), &
                                                            mulliken_grad_r(:, :, :)
      TYPE(exchange_bvk_kernel)                          :: kernel
      TYPE(cp2k_stream_plan_type)                        :: reference_plan

      gradient = 0.0_wp
      sigma = 0.0_wp
      overlap_grad = (0.0_wp, 0.0_wp)
      nk = SIZE(weights)

      IF (.NOT. ALLOCATED(calc%exchange)) THEN
         CALL fatal_error(error, "g-xTB exchange container is not allocated")
         RETURN
      END IF
      IF (.NOT. ALLOCATED(ecache%raw)) THEN
         CALL fatal_error(error, &
            & "k-point exchange cache is not prepared; call cp2k_prepare_exchange")
         RETURN
      END IF
      IF (ANY(nmesh < 1) .OR. PRODUCT(nmesh) /= nk) THEN
         CALL fatal_error(error, "g-xTB exchange requires a complete regular k mesh")
         RETURN
      END IF
      DO idim = 1, 3
         IF (.NOT. mol%periodic(idim) .AND. nmesh(idim) /= 1) THEN
            CALL fatal_error(error, &
               & "nonperiodic directions cannot be replicated by the g-xTB k mesh")
            RETURN
         END IF
      END DO
      IF (SIZE(kfrac, 1) /= 3 .OR. SIZE(kfrac, 2) /= nk &
         & .OR. SIZE(overlap, 3) /= nk .OR. SIZE(overlap_grad, 3) /= nk &
         & .OR. SIZE(density, 4) /= nk) THEN
         CALL fatal_error(error, &
            & "inconsistent g-xTB whole-mesh gradient array dimensions")
         RETURN
      END IF
      IF (.NOT. ALL(ieee_is_finite(kfrac))) THEN
         CALL fatal_error(error, "g-xTB whole-mesh k points must be finite")
         RETURN
      END IF
      IF (.NOT. ALL(ieee_is_finite(weights)) .OR. &
         & ABS(SUM(weights) - 1.0_wp) > 1.0e-12_wp &
         & .OR. MAXVAL(ABS(weights - 1.0_wp/REAL(nk, wp))) > 1.0e-12_wp) THEN
         CALL fatal_error(error, "g-xTB whole-mesh exchange requires uniform weights")
         RETURN
      END IF

      SELECT TYPE (cache => ecache%raw)
      TYPE is (exchange_cache)
         SELECT TYPE (exchange => calc%exchange)
         TYPE is (exchange_fock)
            IF (SIZE(overlap, 1) /= exchange%nao &
               & .OR. SIZE(overlap, 2) /= exchange%nao &
               & .OR. SIZE(overlap_grad, 1) /= exchange%nao &
               & .OR. SIZE(overlap_grad, 2) /= exchange%nao &
               & .OR. SIZE(density, 1) /= exchange%nao &
               & .OR. SIZE(density, 2) /= exchange%nao &
               & .OR. (SIZE(density, 3) /= 1 .AND. SIZE(density, 3) /= 2) &
               & .OR. SIZE(gradient, 1) /= 3 &
               & .OR. SIZE(gradient, 2) /= mol%nat &
               & .OR. SIZE(sigma, 1) /= 3 .OR. SIZE(sigma, 2) /= 3) THEN
               CALL fatal_error(error, &
                  & "inconsistent g-xTB whole-mesh gradient AO or atom dimensions")
               RETURN
            END IF

            CALL exchange%get_bvk_Kmatrix(mol, nmesh, kernel)
            DO icell = 1, nk
               DO jcell = 1, nk
                  orthogonality = SUM(weights*EXP(CMPLX(0.0_wp, 2.0_wp*pi* &
                     & MATMUL(TRANSPOSE(kfrac), REAL(kernel%reps(:, icell) &
                     & - kernel%reps(:, jcell), wp)), wp)))
                  TARGET = MERGE(1.0_wp, 0.0_wp, icell == jcell)
                  IF (ABS(orthogonality - TARGET) > 1.0e-10_wp) THEN
                     CALL fatal_error(error, &
                        & "g-xTB k points are not a complete Fourier mesh")
                     RETURN
                  END IF
               END DO
            END DO

            ALLOCATE (mulliken_grad_r(exchange%nao, exchange%nao, nk), &
               & bocorr_grad_r(exchange%nao, exchange%nao, nk))
            CALL consumer_get_kgrad_kmesh(exchange, mol, cache, reference_plan, kernel, kfrac, weights, &
               & density, overlap, overlap_grad, mulliken_grad_r, bocorr_grad_r)
            CALL exchange%get_bvk_Kmatrix_derivs(mol, kernel, mulliken_grad_r, &
               & bocorr_grad_r, gradient, sigma)
         class default
            CALL fatal_error(error, &
               & "whole-mesh exchange gradient requires exchange_fock")
         END SELECT
      class default
         CALL fatal_error(error, "unexpected g-xTB exchange cache type")
      END SELECT

   END SUBROUTINE cp2k_exchange_kmesh_gradient

!> Differentiate the g-xTB exchange functional for one complex k-point block.
!>
!> All returned quantities are unweighted.  For a reduced mesh the external
!> driver must expand `gradient` and `sigma` over the complete k star.  The
!> complex Hermitian `overlap_grad` transforms like a density matrix and is
!> intended to be folded to real-space overlap images by the same weighted
!> star transform as an ordinary k-point density.
   SUBROUTINE cp2k_exchange_kpoint_gradient(calc, mol, ecache, density, overlap, &
      & overlap_grad, gradient, sigma, error)
      TYPE(xtb_calculator), INTENT(in)                   :: calc
      TYPE(structure_type), INTENT(in)                   :: mol
      TYPE(container_cache), INTENT(inout)               :: ecache
      COMPLEX(wp), INTENT(in)                            :: density(:, :, :), overlap(:, :)
      COMPLEX(wp), INTENT(out)                           :: overlap_grad(:, :)
      REAL(wp), INTENT(out)                              :: gradient(:, :), sigma(:, :)
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      REAL(wp), ALLOCATABLE                              :: bocorr_grad(:, :), mulliken_grad(:, :)

      gradient = 0.0_wp
      sigma = 0.0_wp
      overlap_grad = (0.0_wp, 0.0_wp)
      IF (.NOT. ALLOCATED(calc%exchange)) THEN
         CALL fatal_error(error, "g-xTB exchange container is not allocated")
         RETURN
      END IF
      IF (.NOT. ALLOCATED(ecache%raw)) THEN
         CALL fatal_error(error, &
            & "k-point exchange cache is not prepared; call cp2k_prepare_exchange")
         RETURN
      END IF

      SELECT TYPE (cache => ecache%raw)
      TYPE is (exchange_cache)
         SELECT TYPE (exchange => calc%exchange)
         TYPE is (exchange_fock)
            ALLOCATE (mulliken_grad(exchange%nao, exchange%nao), &
               & bocorr_grad(exchange%nao, exchange%nao))
            CALL consumer_get_kgrad_kpoint(exchange, mol, cache, density, overlap, &
               & mulliken_grad, bocorr_grad, overlap_grad)
            CALL get_mulliken_derivs_direct(exchange, mol, cache, mulliken_grad, &
               & gradient, sigma)
            CALL get_bocorr_derivs_direct(exchange, mol, cache, bocorr_grad, &
               & gradient, sigma)
         class default
            CALL fatal_error(error, &
               & "k-point exchange gradient requires the exchange_fock implementation")
         END SELECT
      class default
         CALL fatal_error(error, "unexpected g-xTB exchange cache type")
      END SELECT
   END SUBROUTINE cp2k_exchange_kpoint_gradient

!> Complete Gamma-point H0/ACP/exchange and q-vSZP response for CP2K.
   SUBROUTINE cp2k_gamma_gradient(calc, mol, bcache, acache, ecache, wfn_aux, &
      & accuracy, selfenergy, dsedcn, overlap, pot, wfn, wdensity, gradient, sigma)
      TYPE(xtb_calculator), INTENT(inout)                :: calc
      TYPE(structure_type), INTENT(in)                   :: mol
      TYPE(basis_cache), INTENT(in)                      :: bcache
      TYPE(acp_cache), INTENT(inout)                     :: acache
      TYPE(container_cache), INTENT(inout)               :: ecache
      TYPE(wavefunction_type), INTENT(in)                :: wfn_aux
      REAL(wp), INTENT(in)                               :: accuracy, selfenergy(:), dsedcn(:), &
                                                            overlap(:, :)
      TYPE(potential_type), INTENT(in)                   :: pot
      TYPE(wavefunction_type), INTENT(inout)             :: wfn
      REAL(wp), INTENT(inout)                            :: wdensity(:, :, :)
      REAL(wp), INTENT(out)                              :: gradient(:, :), sigma(:, :)

      REAL(wp)                                           :: cutoff
      REAL(wp), ALLOCATABLE :: aniso_dip(:, :, :), cn(:), dcndL(:, :, :), dcndr(:, :, :), &
         dEdad(:, :, :), dEdcn(:), dEdcnbas(:), dEdq(:), dEdqbas(:), h1(:, :, :), lattr(:, :)
      TYPE(adjacency_list)                               :: list, list_acp
      TYPE(integral_type)                                :: ints
      TYPE(potential_type)                               :: pot_work

      gradient = 0.0_wp
      sigma = 0.0_wp
      ALLOCATE (aniso_dip(3, 3, mol%nat), dEdcn(mol%nat), dEdad(3, 3, mol%nat), &
         & dEdcnbas(mol%nat), dEdqbas(mol%nat), dEdq(mol%nat), source=0.0_wp)

      cutoff = cp2k_get_acp_cutoff(calc, accuracy)
      CALL get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
      CALL new_adjacency_list(list, mol, lattr, cutoff)
      CALL get_anisotropy(calc%h0, mol, lattr, cutoff, aniso_dip)

      ! CP2K constructs the potential contribution to its KS matrix directly.
      ! The native gradient kernel additionally expects the atom/shell shifts
      ! spread to pot%vao.  Prepare that representation on a private copy so the
      ! persistent CP2K potential is neither changed nor accumulated twice.
      pot_work = pot
      CALL new_integral(ints, SIZE(overlap, 1))
      ints%overlap = overlap
      ALLOCATE (h1(SIZE(overlap, 1), SIZE(overlap, 2), wfn%nspin))
      CALL add_pot_to_h1(calc%bas, ints, pot_work, h1)

      IF (ALLOCATED(calc%ncoord)) THEN
         ALLOCATE (cn(mol%nat), dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat))
         CALL calc%ncoord%get_cn(mol, cn, dcndr, dcndL)
      END IF

      IF (ALLOCATED(calc%exchange)) THEN
         CALL calc%exchange%get_gradient_w_overlap(mol, ecache, wfn, overlap, &
            & wdensity(:, :, 1), gradient, sigma)
      END IF

      IF (ALLOCATED(calc%acp)) THEN
         ! Match the complete ACP traversal used when its Gamma cache was
         ! built, retaining the separate half list for the H0 response below.
         CALL new_adjacency_list(list_acp, mol, lattr, cutoff, complete=.TRUE.)
         CALL get_acp_gradient(mol, lattr, list_acp, calc%bas, bcache, calc%acp, acache, &
            & wfn, dEdcnbas, dEdqbas, gradient, sigma)
      END IF

      CALL updown_to_magnet(wfn%density)
      CALL updown_to_magnet(wdensity)
      CALL get_hamiltonian_gradient(mol, lattr, list, calc%bas, bcache, calc%h0, &
         & selfenergy, dsedcn, aniso_dip, pot_work, wfn%density, wdensity, dEdcn, dEdad, &
         & dEdcnbas, dEdqbas, gradient, sigma)
      CALL magnet_to_updown(wfn%density)

      IF (ALLOCATED(dcndr)) CALL gemv(dcndr, dEdcn, gradient, beta=1.0_wp)
      IF (ALLOCATED(dcndL)) CALL gemv(dcndL, dEdcn, sigma, beta=1.0_wp)

      CALL cp2k_apply_basis_gradient(calc, mol, wfn_aux, dEdcnbas, dEdqbas, &
         & gradient, sigma)

      CALL get_anisotropy_gradient(calc%h0, mol, lattr, cutoff, dEdad, gradient, sigma)
   END SUBROUTINE cp2k_gamma_gradient

!> Image-resolved H0/Pulay and q-vSZP response for CP2K.
!>
!> The real-space matrices use layout `(nao,nao,ntranslation,nspin)` and
!> must use the exact translation ordering supplied here.  ACP and exchange
!> responses remain separate contributions.
   SUBROUTINE cp2k_image_h0_gradient(calc, mol, bcache, wfn_aux, accuracy, &
      & translations, selfenergy, dsedcn, overlap_gamma, pot, density_images, &
      & wdensity_images, gradient, sigma)
      TYPE(xtb_calculator), INTENT(inout)                :: calc
      TYPE(structure_type), INTENT(in)                   :: mol
      TYPE(basis_cache), INTENT(in)                      :: bcache
      TYPE(wavefunction_type), INTENT(in)                :: wfn_aux
      REAL(wp), INTENT(in)                               :: accuracy, translations(:, :), &
                                                            selfenergy(:), dsedcn(:), &
                                                            overlap_gamma(:, :)
      TYPE(potential_type), INTENT(in)                   :: pot
      REAL(wp), INTENT(in)                               :: density_images(:, :, :, :), &
                                                            wdensity_images(:, :, :, :)
      REAL(wp), INTENT(out)                              :: gradient(:, :), sigma(:, :)

      REAL(wp)                                           :: cutoff
      REAL(wp), ALLOCATABLE :: aniso(:, :, :), cn(:), dcndL(:, :, :), dcndr(:, :, :), &
         dEda(:, :, :), dEdcn(:), dEdcnbas(:), dEdq(:), dEdqbas(:), density_work(:, :, :, :), &
         h1(:, :, :), wdensity_work(:, :, :, :)
      TYPE(adjacency_list)                               :: list
      TYPE(integral_type)                                :: ints
      TYPE(potential_type)                               :: pot_work

      gradient = 0.0_wp
      sigma = 0.0_wp
      ALLOCATE (aniso(3, 3, mol%nat), dEdcn(mol%nat), dEda(3, 3, mol%nat), &
         & dEdcnbas(mol%nat), dEdqbas(mol%nat), dEdq(mol%nat), source=0.0_wp)

      cutoff = get_cutoff(calc%bas, accuracy)
      CALL new_adjacency_list(list, mol, translations, cutoff)
      CALL get_anisotropy(calc%h0, mol, translations, cutoff, aniso)

      pot_work = pot
      CALL new_integral(ints, SIZE(overlap_gamma, 1))
      ints%overlap = overlap_gamma
      ALLOCATE (h1(SIZE(overlap_gamma, 1), SIZE(overlap_gamma, 2), &
         & SIZE(density_images, 4)))
      CALL add_pot_to_h1(calc%bas, ints, pot_work, h1)

      IF (ALLOCATED(calc%ncoord)) THEN
         ALLOCATE (cn(mol%nat), dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat))
         CALL calc%ncoord%get_cn(mol, cn, dcndr, dcndL)
      END IF

      density_work = density_images
      wdensity_work = wdensity_images
      CALL updown_to_magnet(density_work)
      CALL updown_to_magnet(wdensity_work)
      CALL get_hamiltonian_gradient(mol, translations, list, calc%bas, bcache, &
         & calc%h0, selfenergy, dsedcn, aniso, pot_work, density_work, &
         & wdensity_work, dEdcn, dEda, dEdcnbas, dEdqbas, gradient, sigma)

      IF (ALLOCATED(dcndr)) CALL gemv(dcndr, dEdcn, gradient, beta=1.0_wp)
      IF (ALLOCATED(dcndL)) CALL gemv(dcndL, dEdcn, sigma, beta=1.0_wp)

      IF (calc%bas%charge_dependent) THEN
         CALL calc%bas%get_basis_gradient(mol, dEdcnbas, dEdqbas, dEdq, gradient, sigma)
      END IF
      IF (ALLOCATED(wfn_aux%dqatdr)) THEN
         CALL gemv(wfn_aux%dqatdr(:, :, :, 1), dEdq, gradient, beta=1.0_wp)
      END IF
      IF (ALLOCATED(wfn_aux%dqatdL)) THEN
         CALL gemv(wfn_aux%dqatdL(:, :, :, 1), dEdq, sigma, beta=1.0_wp)
      END IF

      CALL get_anisotropy_gradient(calc%h0, mol, translations, cutoff, dEda, &
         & gradient, sigma)
   END SUBROUTINE cp2k_image_h0_gradient

!> Extract one effective CGTO in the representation expected by CP2K.
   SUBROUTINE cp2k_get_cgto(bas, cache, ish, izp, iat, ang, nprim, alpha, coeff)
      class(basis_type), INTENT(in) :: bas
      TYPE(basis_cache), INTENT(in) :: cache
      INTEGER, INTENT(in) :: ish, izp, iat
      INTEGER, INTENT(out) :: ang, nprim
      REAL(wp), INTENT(out) :: alpha(:), coeff(:)

      associate (cgto => bas%cgto(ish, izp)%raw)
         ang = cgto%ang
         nprim = cgto%nprim
         alpha = 0.0_wp
         coeff = 0.0_wp
         alpha(:nprim) = cgto%alpha(:nprim)
         CALL cgto%get_coeffs(cache%cgto(ish, iat), coeff(:nprim))
      END associate
   END SUBROUTINE cp2k_get_cgto

!> Evaluate multipole integrals with save_tblite's explicit basis cache.
   pure SUBROUTINE cp2k_multipole_cgto(bas, cache, jsh, jzp, jat, ish, izp, iat, &
      & r2, vec, overlap, dpint, qpint)
      class(basis_type), INTENT(in) :: bas
      TYPE(basis_cache), INTENT(in) :: cache
      INTEGER, INTENT(in) :: jsh, jzp, jat, ish, izp, iat
      REAL(wp), INTENT(in) :: r2, vec(3)
      REAL(wp), INTENT(out) :: overlap(*), dpint(3, *), qpint(6, *)

      CALL multipole_cgto(bas%cgto(jsh, jzp)%raw, bas%cgto(ish, izp)%raw, &
         & cache%cgto(jsh, jat), cache%cgto(ish, iat), r2, vec, bas%intcut, &
         & overlap, dpint, qpint)
   END SUBROUTINE cp2k_multipole_cgto

!> Evaluate multipole derivatives with save_tblite's explicit basis cache.
   pure SUBROUTINE cp2k_multipole_grad_cgto(bas, cache, jsh, jzp, jat, ish, izp, iat, &
      & r2, vec, overlap, dpint, qpint, doverlap, ddpintj, dqpintj, ddpinti, dqpinti)
      class(basis_type), INTENT(in) :: bas
      TYPE(basis_cache), INTENT(in) :: cache
      INTEGER, INTENT(in) :: jsh, jzp, jat, ish, izp, iat
      REAL(wp), INTENT(in) :: r2, vec(3)
      REAL(wp), INTENT(out) :: overlap(*), dpint(3, *), qpint(6, *)
      REAL(wp), INTENT(out) :: doverlap(3, *)
      REAL(wp), INTENT(out) :: ddpintj(3, 3, *), dqpintj(3, 6, *)
      REAL(wp), INTENT(out) :: ddpinti(3, 3, *), dqpinti(3, 6, *)

      CALL multipole_grad_cgto(bas%cgto(jsh, jzp)%raw, bas%cgto(ish, izp)%raw, &
         & cache%cgto(jsh, jat), cache%cgto(ish, iat), r2, vec, bas%intcut, &
         & overlap, dpint, qpint, doverlap, ddpintj, dqpintj, ddpinti, dqpinti)
   END SUBROUTINE cp2k_multipole_grad_cgto

END MODULE tblite_model_bridge
#endif
