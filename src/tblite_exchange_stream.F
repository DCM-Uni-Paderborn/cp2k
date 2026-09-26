#if defined(__TBLITE_GXTB_ACCELERATION)
! Derived from save_tblite consumer orchestration at aad26158b2e83f84b7b7d5b56f361ead4bc02a88.
! SPDX-Identifier: LGPL-3.0-or-later

!> CP2K-owned batched and distributed exchange orchestration.
!> Public provider routines below are mathematical kernels; the common model
!> bridge is shared rather than duplicating mixer/basis orchestration.
MODULE tblite_exchange_stream
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE iso_fortran_env,                 ONLY: int64
   USE mctc_env,                        ONLY: error_type,&
                                              fatal_error,&
                                              wp
   USE mctc_io,                         ONLY: structure_type
   USE mctc_io_constants,               ONLY: pi
   USE tblite_acp_cache,                ONLY: acp_cache
   USE tblite_acp_images,               ONLY: acp_image_type
   USE tblite_acp_response,             ONLY: &
        acp_kmesh_gradient_stream_begin, &
        acp_kmesh_gradient_stream_complex_elements, acp_kmesh_gradient_stream_discard, &
        acp_kmesh_gradient_stream_end, acp_kmesh_gradient_stream_push, &
        acp_kmesh_gradient_stream_real_elements, &
        cp2k_acp_kmesh_gradient_stream_type => acp_kmesh_gradient_stream_type, &
        get_acp_projector_image_gradient
   USE tblite_projector_cache, ONLY: cp2k_acp_projector_cache_type, &
      projector_cache_lookup, projector_cache_store, projector_cache_bypass
   USE tblite_basis_cache,              ONLY: basis_cache
   USE tblite_bvk_transform,            ONLY: bvk_mixed_radix_k_to_r,&
                                              bvk_separable_k_to_r,&
                                              bvk_transform_dense,&
                                              bvk_transform_mixed_radix_fft,&
                                              bvk_transform_separable_dft
   USE tblite_container_cache,          ONLY: container_cache
   USE tblite_stream_plan, ONLY: cp2k_stream_plan_type
   USE tblite_exchange_cache,           ONLY: exchange_cache
   USE tblite_exchange_contractions, ONLY: &
      consumer_get_kfock_kmesh => get_kfock_kmesh, &
      consumer_get_kgrad_kmesh => get_kgrad_kmesh, &
      consumer_get_kfock_stream_apply => get_kfock_stream_apply, &
      consumer_get_kfock_batch_apply => get_kfock_batch_apply, &
      consumer_get_kgrad_batch_apply => get_kgrad_batch_apply
   USE tblite_exchange_fock,            ONLY: exchange_fock
   USE tblite_model_bridge,             ONLY: &
        direct_build_image_acp => cp2k_build_image_acp, &
        direct_build_acp_projectors => cp2k_build_acp_projectors, &
        cp2k_acp_image_gradient, cp2k_acp_kmesh, cp2k_acp_kmesh_gradient, &
        cp2k_apply_basis_gradient, cp2k_build_gamma_h0, cp2k_build_image_h0, cp2k_exchange_kpoint, &
        cp2k_exchange_kpoint_gradient, cp2k_gamma_gradient, cp2k_get_acp_cutoff, cp2k_get_cgto, &
        cp2k_get_mixer_dimension, cp2k_get_scc_mixer_dimension, cp2k_image_h0_gradient, &
        cp2k_multipole_cgto, cp2k_multipole_grad_cgto, cp2k_potential_mixer_type, &
        cp2k_prepare_exchange, cp2k_set_scf_options, cp2k_update_basis, new_cp2k_potential_mixer
   USE tblite_wavefunction_type,        ONLY: wavefunction_type
   USE tblite_xtb_calculator,           ONLY: xtb_calculator
   USE tblite_xtb_h0,                   ONLY: hamiltonian_image_type

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
   PUBLIC :: cp2k_exchange_stream_type, cp2k_exchange_stream_begin, &
      & cp2k_exchange_stream_push, cp2k_exchange_stream_apply, &
      & cp2k_exchange_stream_get, cp2k_exchange_stream_reverse_apply, &
      & cp2k_exchange_stream_reverse_get, cp2k_exchange_stream_end, &
      & cp2k_exchange_stream_reverse_begin, &
      & cp2k_exchange_stream_reverse_push, &
      & cp2k_exchange_stream_reverse_transform_advance, &
      & cp2k_exchange_stream_reverse_batch_apply, &
      & cp2k_exchange_stream_reverse_batch_pull, &
      & cp2k_exchange_stream_reverse_batch_advance, &
      & cp2k_exchange_stream_reverse_result, &
      & cp2k_exchange_stream_has_full_mesh_storage, &
      & cp2k_exchange_stream_batch_bounds, &
      & cp2k_exchange_stream_transform_advance, &
      & cp2k_exchange_stream_batch_apply, &
      & cp2k_exchange_stream_batch_pull, &
      & cp2k_exchange_stream_batch_advance, &
      & cp2k_exchange_stream_batch_result, &
      & cp2k_exchange_stream_peak_complex_elements, &
      & cp2k_exchange_stream_pair_tile_size, cp2k_exchange_stream_image_batch_size
   PUBLIC :: cp2k_exchange_partial_descriptor_type, &
      & cp2k_exchange_partial_forward_type, &
      & cp2k_exchange_partial_reverse_type, &
      & cp2k_exchange_partial_begin, cp2k_exchange_partial_push, &
      & cp2k_exchange_partial_export_forward, &
      & cp2k_exchange_partial_import_forward, &
      & cp2k_exchange_partial_forward_size, &
      & cp2k_exchange_partial_forward_pack, &
      & cp2k_exchange_partial_forward_pack_range, &
      & cp2k_exchange_partial_forward_unpack, &
      & cp2k_exchange_partial_forward_unpack_range, &
      & cp2k_exchange_partial_import_forward_range, &
      & cp2k_exchange_partial_descriptor_size, &
      & cp2k_exchange_partial_descriptor_pack, &
      & cp2k_exchange_partial_descriptor_unpack, &
      & cp2k_exchange_partial_reverse_push, &
      & cp2k_exchange_partial_export_reverse, &
      & cp2k_exchange_partial_import_reverse, &
      & cp2k_exchange_partial_reverse_size, &
      & cp2k_exchange_partial_reverse_pack, &
      & cp2k_exchange_partial_reverse_pack_range, &
      & cp2k_exchange_partial_reverse_unpack, &
      & cp2k_exchange_partial_reverse_unpack_range, &
      & cp2k_exchange_partial_import_reverse_range, &
      & cp2k_exchange_partial_descriptors_match, &
      & cp2k_exchange_partial_set_image_range, &
      & cp2k_exchange_partial_set_pull_ownership, &
      & cp2k_exchange_partial_contributor_advance, &
      & cp2k_exchange_partial_discard
   PUBLIC :: cp2k_image_h0_gradient
   PUBLIC :: cp2k_acp_image_gradient, cp2k_acp_projector_image_gradient, &
      & cp2k_acp_kmesh, cp2k_acp_kmesh_gradient, &
      & cp2k_acp_kmesh_gradient_stream_begin, &
      & cp2k_acp_kmesh_gradient_stream_push, &
      & cp2k_acp_kmesh_gradient_stream_end, &
      & cp2k_acp_kmesh_gradient_stream_discard, &
      & cp2k_acp_kmesh_gradient_stream_real_elements, &
      & cp2k_acp_kmesh_gradient_stream_complex_elements, &
      & cp2k_apply_basis_gradient
   PUBLIC :: acp_image_type, cp2k_acp_kmesh_gradient_stream_type, &
      & hamiltonian_image_type

   INTEGER, PARAMETER, PUBLIC :: cp2k_exchange_stream_reduced = 1
   INTEGER, PARAMETER, PUBLIC :: cp2k_exchange_stream_oracle = 2
   INTEGER, PARAMETER, PUBLIC :: cp2k_exchange_stream_batched = 3
   !> Dense Nk-by-Nk Fourier phase-table oracle.
   INTEGER, PARAMETER, PUBLIC :: cp2k_exchange_transform_dense = &
      & bvk_transform_dense
   !> Factorized direct DFT (not an FFT): O(Nk*sum(nmesh)) characters.
   INTEGER, PARAMETER, PUBLIC :: cp2k_exchange_transform_separable_dft = &
      & bvk_transform_separable_dft
   !> Internal mixed-radix Cooley--Tukey FFT with direct prime leaves.
   INTEGER, PARAMETER, PUBLIC :: cp2k_exchange_transform_mixed_radix_fft = &
      & bvk_transform_mixed_radix_fft
   INTEGER, PARAMETER, PUBLIC :: cp2k_exchange_partial_abi_version = 1
   INTEGER, PARAMETER, PUBLIC :: cp2k_exchange_partial_descriptor_words = 30

   !> CP2K-facing wrapper for g-xTB's native potential mixer.
   !>
   !> CP2K owns the packing of all real Fock degrees of freedom into one flat
   !> vector.  This wrapper deliberately owns only the mixer state and mirrors
   !> g-xTB's iterator schedule: simple mixing with damping 0.2 for the first
   !> three Fock builds, DIIS precollection on the third build, and DIIS from
   !> the fourth build onward.

   !> Stateful collector for Brillouin-zone-coupled exchange.
   !>
   !> The default reduced mode retains three complete weighted k-to-R image
   !> meshes per spin.  The explicitly selected batched mode instead bounds AO
   !> image storage and requires sequential k-block resubmission and Fock pulls
   !> for each batch; it also drops the dense Nk-by-Nk phase tables.  Oracle
   !> mode retains all Bloch blocks and dispatches the unchanged complete-mesh
   !> evaluator, including reverse mode.
   TYPE, PUBLIC :: cp2k_exchange_stream_type
      PRIVATE
      LOGICAL :: active = .FALSE.
      LOGICAL :: applied = .FALSE.
      LOGICAL :: reverse_applied = .FALSE.
      LOGICAL :: reverse_transaction = .FALSE.
      LOGICAL :: reverse_batch_applied = .FALSE.
      LOGICAL :: reverse_result_pulled = .FALSE.
      LOGICAL :: reverse_only = .FALSE.
      LOGICAL :: batch_applied = .FALSE.
      LOGICAL :: batch_result_pulled = .FALSE.
      LOGICAL :: partial_active = .FALSE.
      LOGICAL :: partial_forward_exported = .FALSE.
      LOGICAL :: partial_forward_imported = .FALSE.
      LOGICAL :: partial_reverse_exported = .FALSE.
      LOGICAL :: partial_reverse_imported = .FALSE.
      LOGICAL :: partial_image_distribution_set = .FALSE.
      LOGICAL :: partial_compact_image_storage = .FALSE.
      LOGICAL :: partial_pull_distribution_set = .FALSE.
      LOGICAL :: partial_result_owner = .FALSE.
      INTEGER :: nk = 0
      INTEGER :: nao = 0
      INTEGER :: nspin = 0
      INTEGER :: nsh = 0
      INTEGER :: mode = 0
      INTEGER :: transform_backend = bvk_transform_dense
      INTEGER :: transform_pair_batch_size = 0
      INTEGER :: transform_pair_first = 0
      INTEGER :: transform_pair_count = 0
      INTEGER :: transform_field_count = 0
      INTEGER :: batch_size = 0
      INTEGER :: batch_first = 0
      INTEGER :: batch_count = 0
      INTEGER :: batch_last = 0
      INTEGER :: image_first = 0
      INTEGER :: partial_image_first = 0
      INTEGER :: partial_image_last = 0
      INTEGER :: plan_builds = 0
      INTEGER :: partial_generation = 0
      INTEGER(int64) :: partial_transaction_id(2) = 0_int64
      INTEGER :: nmesh(3) = 0
      INTEGER(int64) :: peak_complex_elements = 0_int64
      INTEGER, ALLOCATABLE :: reps(:, :), inverse_image(:), input_to_grid(:)
      INTEGER(int64), ALLOCATABLE :: density_signature(:), &
         & overlap_signature(:)
      LOGICAL, ALLOCATABLE :: seen(:), forward_pulled(:), reverse_pulled(:), &
         & batch_pulled(:), block_signature_set(:), partial_forward_owned(:), &
         & partial_pull_mask(:)
      COMPLEX(wp), ALLOCATABLE :: inverse_coeff(:), phase_roots(:, :), &
         & twist_phase(:)
      REAL(wp), ALLOCATABLE :: kfrac(:, :), weights(:), vsh(:)
      REAL(wp), ALLOCATABLE :: reverse_gradient(:, :), reverse_sigma(:, :)
      REAL(wp), ALLOCATABLE :: g_onsfx(:, :, :), g_onsri(:, :), &
         & dgdq_onsfx(:, :, :), dgdq_onsri(:, :)
      COMPLEX(wp), ALLOCATABLE :: density(:, :, :, :), overlap(:, :, :), &
         & fock(:, :, :, :), overlap_adjoint(:, :, :)
      COMPLEX(wp), ALLOCATABLE :: amat_r(:, :, :, :), cmat_r(:, :, :, :), &
         & vmat_r(:, :, :, :), amat_partner_r(:, :, :, :), &
         & reverse_amat_r(:, :, :, :), reverse_vmat_r(:, :, :, :), &
         & gdiagP(:, :), gdiagSP(:, :), gdiagSPS(:, :)
      !> AO-pair tiles across the complete regular k mesh.  For nontrivial AO
      !> spaces the first extent is strictly smaller than nao**2, preventing a
      !> persistent complete Nk-by-nao**2 object in compact-transform streams.
      COMPLEX(wp), ALLOCATABLE :: transform_k(:, :, :, :)
      COMPLEX(wp), ALLOCATABLE :: reverse_bdiagSP(:, :), &
         & reverse_bdiagSPS(:, :)
      INTEGER(int64) :: kernel_signature(2) = 0_int64
      INTEGER(int64) :: phase_signature(2) = 0_int64
      INTEGER(int64) :: compact_phase_signature(2) = 0_int64
      INTEGER(int64) :: partial_state_signature(12) = 0_int64
      REAL(wp) :: twist(3) = 0.0_wp
      REAL(wp) :: energy = 0.0_wp
   END TYPE cp2k_exchange_stream_type

   !> Immutable identity of one pre-kernel partial-accumulation generation.
   !>
   !> Components are private so a host can copy and compare descriptors but
   !> cannot manufacture a compatible descriptor by changing individual
   !> fields.  Every import also rebuilds the descriptor from the live stream
   !> and validates the geometry/model/cache state before accepting payloads.
   TYPE, PUBLIC :: cp2k_exchange_partial_descriptor_type
      PRIVATE
      INTEGER :: abi_version = 0
      INTEGER :: scalar_bits = 0
      INTEGER :: phase = 0
      INTEGER :: generation = 0
      INTEGER(int64) :: transaction_id(2) = 0_int64
      INTEGER :: nk = 0
      INTEGER :: nao = 0
      INTEGER :: nspin = 0
      INTEGER :: nsh = 0
      INTEGER :: nmesh(3) = 0
      INTEGER :: batch_first = 0
      INTEGER :: batch_count = 0
      INTEGER :: batch_capacity = 0
      INTEGER :: image_first = 0
      INTEGER :: image_last = 0
      INTEGER(int64) :: state_signature(12) = 0_int64
   END TYPE cp2k_exchange_partial_descriptor_type

   !> Opaque contiguous complex payload for the additive forward k-to-R map.
   TYPE, PUBLIC :: cp2k_exchange_partial_forward_type
      PRIVATE
      INTEGER :: abi_version = 0
      INTEGER :: image_first = 0
      INTEGER :: image_count = 0
      INTEGER(int64) :: descriptor_signature(2) = 0_int64
      COMPLEX(wp), ALLOCATABLE :: values(:)
   END TYPE cp2k_exchange_partial_forward_type

   !> Opaque contiguous complex payload for the additive reverse source map.
   TYPE, PUBLIC :: cp2k_exchange_partial_reverse_type
      PRIVATE
      INTEGER :: abi_version = 0
      INTEGER :: image_first = 0
      INTEGER :: image_count = 0
      INTEGER(int64) :: descriptor_signature(2) = 0_int64
      COMPLEX(wp), ALLOCATABLE :: values(:)
   END TYPE cp2k_exchange_partial_reverse_type

CONTAINS

!> Choose a compact AO-pair tile independently of the real-space image batch.
!> Each tile resubmits the complete physical k mesh, so tying its size to a
!> small image batch repeats symmetry expansion and matrix products needlessly.
!> Limit the staged fields plus three transform work planes to 64 MiB, except
!> when even one pair exceeds that budget. Never retain a complete AO-pair mesh
!> for a non-scalar basis. The caller supplies positive dimensions and either
!> three forward fields or six reverse fields.
   PURE INTEGER FUNCTION cp2k_exchange_stream_pair_tile_size(nao, nk, nspin, nfields) RESULT(tile_size)
      INTEGER, INTENT(IN) :: nao, nk, nspin, nfields
      INTEGER(int64), PARAMETER :: target_bytes = 64_int64*1024_int64**2
      INTEGER(int64) :: bytes_per_pair, pair_limit

      pair_limit = MAX(1_int64, INT(nao, int64)**2 - 1_int64)
      bytes_per_pair = 16_int64*INT(nk, int64)*(INT(nspin, int64)*INT(nfields, int64) + 3_int64)
      tile_size = INT(MAX(1_int64, MIN(pair_limit, target_bytes/bytes_per_pair)))
   END FUNCTION cp2k_exchange_stream_pair_tile_size

!> Bound the six image fields needed by forward/reverse exchange to 64 MiB.
!> Small image batches repeatedly expand the same stars and transform the same
!> AO-pair tiles. Use the largest fitting batch, but never the complete image
!> mesh. One image is the irreducible minimum when the basis alone exceeds
!> the workspace target. Transform tiles have their own independent bound.
   PURE INTEGER FUNCTION cp2k_exchange_stream_image_batch_size(nao, nk, nspin) RESULT(batch_size)
      INTEGER, INTENT(IN) :: nao, nk, nspin
      INTEGER(int64), PARAMETER :: target_bytes = 64_int64*1024_int64**2
      INTEGER(int64) :: bytes_per_image, image_limit

      image_limit = MAX(1_int64, INT(nk, int64) - 1_int64)
      bytes_per_image = 16_int64*6_int64*INT(nao, int64)**2*INT(nspin, int64)
      batch_size = INT(MAX(1_int64, MIN(image_limit, target_bytes/bytes_per_image)))
   END FUNCTION cp2k_exchange_stream_image_batch_size

!> Whether a public exchange-transform selector is supported.
   pure LOGICAL FUNCTION valid_exchange_transform_backend(backend) RESULT(valid)
      INTEGER, INTENT(in)                                :: backend

      valid = backend == bvk_transform_dense &
         & .OR. backend == bvk_transform_separable_dft &
         & .OR. backend == bvk_transform_mixed_radix_fft
   END FUNCTION valid_exchange_transform_backend

!> Whether a selector consumes the compact regular-grid Fourier plan.
   pure LOGICAL FUNCTION compact_exchange_transform_backend(backend) &
      & RESULT(compact)
      INTEGER, INTENT(in)                                :: backend

      compact = backend == bvk_transform_separable_dft &
         & .OR. backend == bvk_transform_mixed_radix_fft
   END FUNCTION compact_exchange_transform_backend

!> Return a real-space cutoff that covers both the valence basis and the
!> auxiliary ACP projectors used by g-xTB.

!> Construct a CP2K-facing g-xTB potential mixer.

!> Reset a CP2K-facing g-xTB potential mixer for a flat Fock vector.

!> Advance g-xTB's native potential mixer from a raw to a mixed Fock vector.

!> Return the RMS residual from the most recent raw Fock vector.

!> Store the current extrapolated vector in the simple mixer cache.

!> Store the current extrapolated vector in the DIIS mixer cache.

!> Return the CP2K-owned SCC mixer dimension without the spin multiplicity.
!>
!> The orbital density is deliberately not included here, even when a
!> contribution advertises an orbital-resolved density potential.  CP2K
!> diagonalizes and stores that density itself; only charges and multipoles
!> are exchanged with its SCC mixer through the compatibility interface.

!> Backward-compatible name for the CP2K SCC mixer dimension.

!> Apply SCC settings using the iterator API used by save_tblite.

!> Construct the immutable auxiliary charges and cache the effective basis.

!> Build the complete Gamma-point overlap and one-electron g-xTB Hamiltonian.
!>
!> The dense matrices are deliberately returned in tblite's atom/AO ordering.
!> CP2K can scatter them into its block-sparse matrices while retaining the
!> element-wise basis metadata used for block dimensions and cutoffs.

!> Build image-resolved q-vSZP overlap, multipoles, and H0 matrices.
!>
!> This routine deliberately covers the one-electron H0 only.  The atomic
!> correction potential is a projector convolution and therefore requires a
!> separate image-resolved implementation before a complete k-point
!> Hamiltonian can be assembled.

!> Build the image-resolved separable ACP contribution.
!>
!> The returned Hamiltonian translations are generally a larger difference
!> set than the translations supplied for the projector integrals.  External
!> drivers must therefore scatter this set independently of the pair-local H0
!> image set.
   SUBROUTINE cp2k_build_image_acp(calc, mol, bcache, acache, translations, &
      & images, hamiltonian_gamma, use_projector_cache, projector_cache)
      TYPE(xtb_calculator), INTENT(inout)                :: calc
      TYPE(structure_type), INTENT(in)                   :: mol
      TYPE(basis_cache), INTENT(in)                      :: bcache
      TYPE(acp_cache), INTENT(inout)                     :: acache
      REAL(wp), INTENT(in)                               :: translations(:, :)
      TYPE(acp_image_type), INTENT(out)                  :: images
      REAL(wp), INTENT(out), OPTIONAL                    :: hamiltonian_gamma(:, :)
      LOGICAL, INTENT(in), OPTIONAL                      :: use_projector_cache
      TYPE(cp2k_acp_projector_cache_type), INTENT(INOUT), OPTIONAL :: projector_cache
      LOGICAL :: enabled, hit

      enabled = .TRUE.
      hit = .FALSE.
      IF (PRESENT(use_projector_cache)) enabled = use_projector_cache
      IF (PRESENT(projector_cache)) THEN
         IF (enabled) THEN
            CALL projector_cache_lookup(projector_cache, mol, calc%bas, bcache, &
                                        calc%acp, acache, translations, images, hit)
         ELSE
            CALL projector_cache_bypass(projector_cache)
         END IF
      END IF
      ! This explicit full-image oracle always forms H(R). Native compact
      ! calculations use the cached projector-only entry below instead.
      CALL direct_build_image_acp(calc, mol, bcache, acache, translations, images, hamiltonian_gamma)
      IF (PRESENT(projector_cache)) THEN
         IF (enabled .AND. .NOT. hit) CALL projector_cache_store(projector_cache, images)
      END IF
   END SUBROUTINE cp2k_build_image_acp

!> Build only the compact projector-overlap images required by the direct
!> Bloch-space ACP forward and reverse sweeps.  Unlike
!> `cp2k_build_image_acp`, this deliberately does not materialize the
!> quadratic translation-difference set or H(R).
   SUBROUTINE cp2k_build_acp_projectors(calc, mol, bcache, acache, translations, &
      & images, use_projector_cache, projector_cache)
      TYPE(xtb_calculator), INTENT(inout)                :: calc
      TYPE(structure_type), INTENT(in)                   :: mol
      TYPE(basis_cache), INTENT(in)                      :: bcache
      TYPE(acp_cache), INTENT(inout)                     :: acache
      REAL(wp), INTENT(in)                               :: translations(:, :)
      TYPE(acp_image_type), INTENT(out)                  :: images
      LOGICAL, INTENT(in), OPTIONAL                      :: use_projector_cache
      TYPE(cp2k_acp_projector_cache_type), INTENT(INOUT), OPTIONAL :: projector_cache
      LOGICAL :: enabled, hit

      enabled = .TRUE.
      IF (PRESENT(use_projector_cache)) enabled = use_projector_cache
      IF (PRESENT(projector_cache)) THEN
         IF (enabled) THEN
            CALL projector_cache_lookup(projector_cache, mol, calc%bas, bcache, &
                                        calc%acp, acache, translations, images, hit)
            IF (hit) RETURN
         ELSE
            CALL projector_cache_bypass(projector_cache)
         END IF
      END IF
      CALL direct_build_acp_projectors(calc, mol, bcache, acache, translations, images)
      IF (PRESENT(projector_cache)) THEN
         IF (enabled) CALL projector_cache_store(projector_cache, images)
      END IF
   END SUBROUTINE cp2k_build_acp_projectors

!> Image-resolved ACP response for CP2K k-point density matrices.
!>
!> `density(:, :, R, spin)` must be in save_tblite's AO image orientation and
!> dual to the Hamiltonian returned by `cp2k_build_image_acp`, so that its
!> energy contribution is the element-wise contraction over equal images.
!> CP2K's irreducible-star weights must already be included in these real-space
!> densities.

!> Contract compact ACP projector images with a sparse CP2K image density.
!>
!> Density translations absent from CP2K's overlap support contribute exact
!> zero blocks.  No quadratic ACP Hamiltonian difference set is materialized.
   SUBROUTINE cp2k_acp_projector_image_gradient(calc, mol, bcache, acache, &
      & images, density_translation, density, dEdcnbas, dEdqbas, gradient, &
      & sigma, error)
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
         CALL fatal_error(error, &
            & "g-xTB ACP projector-image gradient requires an ACP model")
         RETURN
      END IF
      CALL get_acp_projector_image_gradient(mol, calc%bas, bcache, calc%acp, &
         & acache, images, density_translation, density, dEdcnbas, dEdqbas, &
         & gradient, sigma, error)
   END SUBROUTINE cp2k_acp_projector_image_gradient

!> Evaluate the separable ACP directly on CP2K Bloch points.

!> Differentiate the complete weighted Bloch contraction of the ACP.

!> Begin a bounded exactly-once ACP k-mesh derivative stream.
   SUBROUTINE cp2k_acp_kmesh_gradient_stream_begin(stream, calc, mol, images, &
      & nspin, nk_expected, error)
      TYPE(cp2k_acp_kmesh_gradient_stream_type), &
         INTENT(inout)                                   :: stream
      TYPE(xtb_calculator), INTENT(in)                   :: calc
      TYPE(structure_type), INTENT(in)                   :: mol
      TYPE(acp_image_type), INTENT(in)                   :: images
      INTEGER, INTENT(in)                                :: nspin, nk_expected
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      IF (.NOT. ALLOCATED(calc%acp)) THEN
         CALL fatal_error(error, "g-xTB ACP container is not allocated")
         RETURN
      END IF
      CALL acp_kmesh_gradient_stream_begin(stream, mol, images, nspin, &
         & nk_expected, error)
   END SUBROUTINE cp2k_acp_kmesh_gradient_stream_begin

!> Add one ordered CP2K density block to an ACP derivative stream.
   SUBROUTINE cp2k_acp_kmesh_gradient_stream_push(stream, images, ik, kfrac, &
      & weight, density, error)
      TYPE(cp2k_acp_kmesh_gradient_stream_type), &
         INTENT(inout)                                   :: stream
      TYPE(acp_image_type), INTENT(in)                   :: images
      INTEGER, INTENT(in)                                :: ik
      REAL(wp), INTENT(in)                               :: kfrac(:), weight
      COMPLEX(wp), INTENT(in)                            :: density(:, :, :)
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      CALL acp_kmesh_gradient_stream_push(stream, images, ik, kfrac, weight, &
         & density, error)
   END SUBROUTINE cp2k_acp_kmesh_gradient_stream_push

!> Finish an ACP derivative stream and apply projector derivatives once.
   SUBROUTINE cp2k_acp_kmesh_gradient_stream_end(stream, calc, mol, bcache, &
      & acache, images, dEdcnbas, dEdqbas, gradient, sigma, error)
      TYPE(cp2k_acp_kmesh_gradient_stream_type), &
         INTENT(inout)                                   :: stream
      TYPE(xtb_calculator), INTENT(in)                   :: calc
      TYPE(structure_type), INTENT(in)                   :: mol
      TYPE(basis_cache), INTENT(in)                      :: bcache
      TYPE(acp_cache), INTENT(in)                        :: acache
      TYPE(acp_image_type), INTENT(in)                   :: images
      REAL(wp), INTENT(inout)                            :: dEdcnbas(:), dEdqbas(:), gradient(:, :), &
                                                            sigma(:, :)
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      IF (.NOT. ALLOCATED(calc%acp)) THEN
         CALL fatal_error(error, "g-xTB ACP container is not allocated")
         RETURN
      END IF
      CALL acp_kmesh_gradient_stream_end(stream, mol, calc%bas, bcache, &
         & calc%acp, acache, images, dEdcnbas, dEdqbas, gradient, sigma, error)
   END SUBROUTINE cp2k_acp_kmesh_gradient_stream_end

!> Discard an ACP derivative stream.
   SUBROUTINE cp2k_acp_kmesh_gradient_stream_discard(stream)
      TYPE(cp2k_acp_kmesh_gradient_stream_type), &
         INTENT(inout)                                   :: stream

      CALL acp_kmesh_gradient_stream_discard(stream)
   END SUBROUTINE cp2k_acp_kmesh_gradient_stream_discard

!> Count real scalar elements owned by an ACP derivative stream.
   INTEGER FUNCTION cp2k_acp_kmesh_gradient_stream_real_elements(stream) &
      & RESULT(nelements)
      TYPE(cp2k_acp_kmesh_gradient_stream_type), &
         INTENT(in)                                      :: stream

      nelements = acp_kmesh_gradient_stream_real_elements(stream)
   END FUNCTION cp2k_acp_kmesh_gradient_stream_real_elements

!> Count complex scalar elements owned by an ACP derivative stream.
   INTEGER FUNCTION cp2k_acp_kmesh_gradient_stream_complex_elements(stream) &
      & RESULT(nelements)
      TYPE(cp2k_acp_kmesh_gradient_stream_type), &
         INTENT(in)                                      :: stream

      nelements = acp_kmesh_gradient_stream_complex_elements(stream)
   END FUNCTION cp2k_acp_kmesh_gradient_stream_complex_elements

!> Apply a component's q-vSZP basis response and external-charge chain rule.
!>
!> Image-resolved one-electron components return their own additive
!> `dE/dCN_basis` and `dE/dq_basis`.  This helper maps exactly that component
!> onto Cartesian and strain derivatives.  Calling it separately for H0 and
!> ACP therefore remains additive and avoids applying either response twice.

!> Prepare geometry- and charge-dependent exchange kernels.
!>
!> This routine must be called once, serially, before evaluating the k-point
!> blocks of every SCC iteration.  Rebuilding here also makes changes of the
!> geometry or lattice visible to the Wigner-Seitz and bond-order kernels.

!> Evaluate the g-xTB exchange operator for one complex k-point block.
!>
!> The input P(k) must not contain a k-point or star weight.  The caller owns
!> all Brillouin-zone weights and symmetry expansion.  This routine returns
!> an unweighted Hermitian F(k), an atom/shell-resolved potential, and an
!> unweighted exchange energy.  For a reduced k mesh the shell potential must
!> be permuted with every star operation before the weighted star sum is made.
!> Call cp2k_prepare_exchange serially before entering the k-point loop; this
!> evaluator then only reads the shared cache and can be called concurrently.

!> Evaluate one regular-grid Fourier character from a compact immutable plan.
!>
!> The supplied representative may be canonical or a difference of two
!> representatives.  Reducing each integer product modulo its mesh extent
!> makes the result independent of representative wrapping.  No transcendental
!> operation and no dense Nk-by-Nk phase storage is required here.
   pure COMPLEX(wp) FUNCTION regular_grid_phase_forward(nmesh, input_to_grid, &
      & representative, phase_roots, twist_phase, ik) RESULT(phase)
      INTEGER, INTENT(in)                                :: nmesh(3), input_to_grid(:), &
                                                            representative(3)
      COMPLEX(wp), INTENT(in)                            :: phase_roots(:, :), twist_phase
      INTEGER, INTENT(in)                                :: ik

      INTEGER                                            :: grid(3), idim, igrid, root
      INTEGER(int64)                                     :: root64

      igrid = input_to_grid(ik) - 1
      grid(1) = MODULO(igrid, nmesh(1))
      grid(2) = MODULO(igrid/nmesh(1), nmesh(2))
      grid(3) = igrid/(nmesh(1)*nmesh(2))
      phase = twist_phase
      DO idim = 1, 3
         root64 = MODULO(INT(grid(idim), int64) &
            & *INT(MODULO(representative(idim), nmesh(idim)), int64), &
            & INT(nmesh(idim), int64))
         root = INT(root64)
         phase = phase*phase_roots(root + 1, idim)
      END DO
   END FUNCTION regular_grid_phase_forward

!> Verify the roundoff-exact finite-group grid required by factorized DFTs.
!> The dense oracle deliberately accepts the looser geometric plan tolerance;
!> a factorized transform must not silently project such perturbed points back
!> onto the nearest regular grid.
   pure LOGICAL FUNCTION machine_regular_bvk_plan(plan, nmesh, kfrac) &
      & RESULT(regular)
      TYPE(cp2k_stream_plan_type), INTENT(in)                   :: plan
      INTEGER, INTENT(in)                                :: nmesh(3)
      REAL(wp), INTENT(in)                               :: kfrac(:, :)

      INTEGER                                            :: grid(3), idim, igrid, ik
      REAL(wp)                                           :: canonical, periodic_delta

      regular = .FALSE.
      IF (.NOT. ALLOCATED(plan%bvk_input_to_grid)) RETURN
      IF (SIZE(plan%bvk_input_to_grid) /= SIZE(kfrac, 2)) RETURN
      DO ik = 1, SIZE(kfrac, 2)
         igrid = plan%bvk_input_to_grid(ik) - 1
         grid(1) = MODULO(igrid, nmesh(1))
         grid(2) = MODULO(igrid/nmesh(1), nmesh(2))
         grid(3) = igrid/(nmesh(1)*nmesh(2))
         DO idim = 1, 3
            canonical = plan%bvk_twist(idim) &
               & + REAL(grid(idim), wp)/REAL(nmesh(idim), wp)
            periodic_delta = MODULO(kfrac(idim, ik), 1.0_wp) &
               & - MODULO(canonical, 1.0_wp)
            periodic_delta = periodic_delta - NINT(periodic_delta)
            IF (ABS(periodic_delta) > 64.0_wp*EPSILON(1.0_wp)) RETURN
         END DO
      END DO
      regular = .TRUE.
   END FUNCTION machine_regular_bvk_plan

!> Build or reuse the geometry- and mesh-static BvK kernel/Fourier plan.
   SUBROUTINE prepare_bvk_plan(exchange, mol, plan, nmesh, kfrac, weights, &
      & error, retain_dense_phases)
      TYPE(exchange_fock), INTENT(in)                    :: exchange
      TYPE(structure_type), INTENT(in)                   :: mol
      TYPE(cp2k_stream_plan_type), INTENT(inout)                :: plan
      INTEGER, INTENT(in)                                :: nmesh(3)
      REAL(wp), INTENT(in)                               :: kfrac(:, :), weights(:)
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error
      LOGICAL, INTENT(in), OPTIONAL                      :: retain_dense_phases

      COMPLEX(wp)                                        :: orthogonality, phase
      INTEGER                                            :: grid(3), icell, idim, igrid, ik, jcell, &
                                                            nearest, nk, root
      LOGICAL                                            :: need_dense_phases, reuse_kernel
      LOGICAL, ALLOCATABLE                               :: seen(:)
      REAL(wp)                                           :: angle, scaled, TARGET
      REAL(wp), ALLOCATABLE                              :: kred(:, :)

      nk = SIZE(weights)
      reuse_kernel = .FALSE.
      need_dense_phases = .TRUE.
      IF (PRESENT(retain_dense_phases)) &
         & need_dense_phases = retain_dense_phases
      IF (plan%bvk_matches(mol, nmesh, kfrac, weights, &
         & require_phases=.FALSE.) &
         & .AND. plan%model_matches(exchange)) THEN
         IF (ANY(SHAPE(plan%bvk_kernel%g_mulliken_r) /= &
            & [exchange%nsh, exchange%nsh, nk])) THEN
            plan%bvk_plan_valid = .FALSE.
         ELSE IF (.NOT. need_dense_phases) THEN
            IF (ALLOCATED(plan%bvk_phase_forward)) &
               & DEALLOCATE (plan%bvk_phase_forward)
            IF (ALLOCATED(plan%bvk_phase_inverse)) &
               & DEALLOCATE (plan%bvk_phase_inverse)
            RETURN
         ELSE IF (ALLOCATED(plan%bvk_phase_forward) &
            & .AND. ALLOCATED(plan%bvk_phase_inverse)) THEN
            IF (ALL(SHAPE(plan%bvk_phase_forward) == [nk, nk]) &
               & .AND. ALL(SHAPE(plan%bvk_phase_inverse) == [nk, nk])) RETURN
            reuse_kernel = .TRUE.
         ELSE
            reuse_kernel = .TRUE.
         END IF
      END IF
      ALLOCATE (kred(3, nk))
      kred = MODULO(kfrac, 1.0_wp)
      IF (.NOT. reuse_kernel) THEN
         plan%bvk_plan_valid = .FALSE.

         IF (ALLOCATED(plan%bvk_input_to_grid)) &
            & DEALLOCATE (plan%bvk_input_to_grid)
         IF (ALLOCATED(plan%bvk_grid_to_input)) &
            & DEALLOCATE (plan%bvk_grid_to_input)
         ALLOCATE (plan%bvk_input_to_grid(nk), plan%bvk_grid_to_input(nk), source=0)
         ALLOCATE (seen(nk), source=.FALSE.)
         DO idim = 1, 3
            plan%bvk_twist(idim) = kred(idim, 1) &
               & - REAL(FLOOR(kred(idim, 1)*REAL(nmesh(idim), wp)), wp) &
               & /REAL(nmesh(idim), wp)
         END DO
         DO ik = 1, nk
            DO idim = 1, 3
               scaled = (kred(idim, ik) - plan%bvk_twist(idim)) &
                       & *REAL(nmesh(idim), wp)
               nearest = NINT(scaled)
               IF (ABS(scaled - REAL(nearest, wp)) > 1.0e-10_wp) THEN
                  CALL fatal_error(error, &
                     & "g-xTB k points do not share a regular-grid twist")
                  RETURN
               END IF
               grid(idim) = MODULO(nearest, nmesh(idim))
            END DO
            igrid = 1 + grid(1) + nmesh(1)*(grid(2) + nmesh(2)*grid(3))
            IF (seen(igrid)) THEN
               CALL fatal_error(error, &
                  & "g-xTB k-point mesh contains a duplicate grid point")
               RETURN
            END IF
            seen(igrid) = .TRUE.
            plan%bvk_input_to_grid(ik) = igrid
            plan%bvk_grid_to_input(igrid) = ik
         END DO
         IF (.NOT. ALL(seen)) THEN
            CALL fatal_error(error, "g-xTB k-point mesh is incomplete")
            RETURN
         END IF

         ! Construct the expensive image kernel only after the O(Nk) mesh proof.
         CALL exchange%get_bvk_Kmatrix(mol, nmesh, plan%bvk_kernel)
         IF (ANY(SHAPE(plan%bvk_kernel%reps) /= [3, nk]) &
            & .OR. ANY(SHAPE(plan%bvk_kernel%g_mulliken_r) /= &
            & [exchange%nsh, exchange%nsh, nk]) &
            & .OR. ANY(SHAPE(plan%bvk_kernel%g_bocorr_r) /= &
            & [mol%nat, mol%nat, nk])) THEN
            CALL fatal_error(error, "g-xTB BvK kernel has inconsistent dimensions")
            RETURN
         END IF

         ! A regular mesh needs only one roots-of-unity vector per direction and
         ! one twist character per image.  This O(Nk+3*max(nmesh)) plan is retained
         ! with the expensive BvK kernel and lets bounded batching evaluate every
         ! forward/inverse character without rebuilding dense Nk-by-Nk phases or
         ! calling EXP in its repeated push/pull loops.
         IF (ALLOCATED(plan%bvk_phase_roots)) &
            & DEALLOCATE (plan%bvk_phase_roots)
         IF (ALLOCATED(plan%bvk_twist_phase)) &
            & DEALLOCATE (plan%bvk_twist_phase)
         ALLOCATE (plan%bvk_phase_roots(MAXVAL(nmesh), 3), &
            & source=(1.0_wp, 0.0_wp))
         DO idim = 1, 3
            DO root = 0, nmesh(idim) - 1
               angle = 2.0_wp*pi*REAL(root, wp)/REAL(nmesh(idim), wp)
               plan%bvk_phase_roots(root + 1, idim) = &
                  & EXP(CMPLX(0.0_wp, angle, wp))
            END DO
         END DO
         ALLOCATE (plan%bvk_twist_phase(nk))
         DO icell = 1, nk
            angle = 2.0_wp*pi*DOT_PRODUCT(plan%bvk_twist, &
               & REAL(plan%bvk_kernel%reps(:, icell), wp))
            plan%bvk_twist_phase(icell) = EXP(CMPLX(0.0_wp, angle, wp))
         END DO
      END IF

      IF (ALLOCATED(plan%bvk_phase_forward)) &
         & DEALLOCATE (plan%bvk_phase_forward)
      IF (ALLOCATED(plan%bvk_phase_inverse)) &
         & DEALLOCATE (plan%bvk_phase_inverse)
      IF (need_dense_phases) THEN
         ALLOCATE (plan%bvk_phase_forward(SIZE(plan%bvk_kernel%reps, 2), nk), &
            & plan%bvk_phase_inverse(nk, SIZE(plan%bvk_kernel%reps, 2)))
         DO ik = 1, nk
            DO icell = 1, SIZE(plan%bvk_kernel%reps, 2)
               angle = 2.0_wp*pi*DOT_PRODUCT(kred(:, ik), &
                  & REAL(plan%bvk_kernel%reps(:, icell), wp))
               phase = EXP(CMPLX(0.0_wp, angle, wp))
               plan%bvk_phase_forward(icell, ik) = phase
               plan%bvk_phase_inverse(ik, icell) = weights(ik)*CONJG(phase)
            END DO
         END DO
      END IF

      ! Keep a dense phase-orthogonality oracle for small meshes.  The unique
      ! integer-grid map above is the O(Nk) production proof for larger meshes
      ! and is also the ordering descriptor consumed by the internal FFT backend.
      IF (nk <= 64) THEN
         DO icell = 1, nk
            DO jcell = 1, nk
               IF (need_dense_phases) THEN
                  orthogonality = SUM(plan%bvk_phase_forward(icell, :) &
                                     & *plan%bvk_phase_inverse(:, jcell))
               ELSE
                  orthogonality = (0.0_wp, 0.0_wp)
                  DO ik = 1, nk
                     orthogonality = orthogonality &
                        & + weights(ik)*regular_grid_phase_forward(nmesh, &
                        & plan%bvk_input_to_grid, &
                        & plan%bvk_kernel%reps(:, icell) &
                        & - plan%bvk_kernel%reps(:, jcell), &
                        & plan%bvk_phase_roots, &
                        & plan%bvk_twist_phase(icell) &
                        & *CONJG(plan%bvk_twist_phase(jcell)), ik)
                  END DO
               END IF
               TARGET = MERGE(1.0_wp, 0.0_wp, icell == jcell)
               IF (ABS(orthogonality - TARGET) > 1.0e-10_wp) THEN
                  CALL fatal_error(error, &
                     & "g-xTB k points are not a complete Fourier mesh")
                  RETURN
               END IF
            END DO
         END DO
      END IF
      IF (.NOT. reuse_kernel) THEN
         CALL plan%set_model_signature(exchange)
         CALL plan%set_bvk_signature(mol, nmesh, kfrac, weights)
      END IF
   END SUBROUTINE prepare_bvk_plan

!> Evaluate BvK-supercell-equivalent exchange on a complete regular k mesh.
!>
!> Density, overlap, and Fock contain unweighted Bloch blocks.  Energy and
!> shell potential are returned as complete Brillouin-zone contractions per
!> primitive cell and must not be weighted a second time by CP2K.
   SUBROUTINE cp2k_exchange_kmesh(calc, mol, ecache, plan, wfn, nmesh, kfrac, weights, &
      & density, overlap, fock, vsh, energy, error, transform_backend)
      TYPE(xtb_calculator), INTENT(inout)                :: calc
      TYPE(structure_type), INTENT(in)                   :: mol
      TYPE(container_cache), INTENT(inout)               :: ecache
      TYPE(cp2k_stream_plan_type), INTENT(inout) :: plan
      TYPE(wavefunction_type), INTENT(in)                :: wfn
      INTEGER, INTENT(in)                                :: nmesh(3)
      REAL(wp), INTENT(in)                               :: kfrac(:, :), weights(:)
      COMPLEX(wp), INTENT(in)                            :: density(:, :, :, :), overlap(:, :, :)
      COMPLEX(wp), INTENT(out)                           :: fock(:, :, :, :)
      REAL(wp), INTENT(out)                              :: vsh(:), energy
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error
      INTEGER, INTENT(in), OPTIONAL                      :: transform_backend

      INTEGER                                            :: backend, idim, nk

      backend = bvk_transform_dense
      IF (PRESENT(transform_backend)) backend = transform_backend
      IF (.NOT. valid_exchange_transform_backend(backend)) THEN
         CALL fatal_error(error, "unknown g-xTB k-to-R transform backend")
         RETURN
      END IF
      nk = SIZE(weights)
      IF (.NOT. complete_regular_mesh(nmesh, nk)) THEN
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
      IF (compact_exchange_transform_backend(backend) .AND. &
         & MAXVAL(ABS(weights - 1.0_wp/REAL(nk, wp))) > &
         & 64.0_wp*EPSILON(1.0_wp)) THEN
         IF (backend == bvk_transform_separable_dft) THEN
            CALL fatal_error(error, &
               & "separable-DFT g-xTB exchange requires machine-uniform weights")
         ELSE
            CALL fatal_error(error, &
               & "mixed-radix FFT g-xTB exchange requires machine-uniform weights")
         END IF
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
            CALL prepare_bvk_plan(exchange, mol, plan, nmesh, kfrac, weights, &
               & error, retain_dense_phases=backend == bvk_transform_dense)
            IF (ALLOCATED(error)) RETURN
            IF (compact_exchange_transform_backend(backend) .AND. &
               & .NOT. machine_regular_bvk_plan(plan, nmesh, kfrac)) THEN
               IF (backend == bvk_transform_separable_dft) THEN
                  CALL fatal_error(error, &
                     & "separable-DFT g-xTB exchange requires a machine-regular grid")
               ELSE
                  CALL fatal_error(error, &
                     & "mixed-radix FFT g-xTB exchange requires a machine-regular grid")
               END IF
               RETURN
            END IF
            CALL consumer_get_kfock_kmesh(exchange, mol, cache, plan, plan%bvk_kernel, &
               & kfrac, weights, density, overlap, fock, vsh, energy, backend)
         class default
            CALL fatal_error(error, &
               & "whole-mesh exchange requires the exchange_fock implementation")
         END SELECT
      class default
         CALL fatal_error(error, "unexpected g-xTB exchange cache type")
      END SELECT

   END SUBROUTINE cp2k_exchange_kmesh

!> Release all storage owned by a Brillouin-zone exchange stream.
   SUBROUTINE clear_exchange_stream(stream)
      TYPE(cp2k_exchange_stream_type), INTENT(inout)     :: stream

      IF (ALLOCATED(stream%seen)) DEALLOCATE (stream%seen)
      IF (ALLOCATED(stream%forward_pulled)) DEALLOCATE (stream%forward_pulled)
      IF (ALLOCATED(stream%reverse_pulled)) DEALLOCATE (stream%reverse_pulled)
      IF (ALLOCATED(stream%batch_pulled)) DEALLOCATE (stream%batch_pulled)
      IF (ALLOCATED(stream%partial_pull_mask)) &
         & DEALLOCATE (stream%partial_pull_mask)
      IF (ALLOCATED(stream%partial_forward_owned)) &
         & DEALLOCATE (stream%partial_forward_owned)
      IF (ALLOCATED(stream%block_signature_set)) &
         & DEALLOCATE (stream%block_signature_set)
      IF (ALLOCATED(stream%reps)) DEALLOCATE (stream%reps)
      IF (ALLOCATED(stream%inverse_image)) DEALLOCATE (stream%inverse_image)
      IF (ALLOCATED(stream%input_to_grid)) DEALLOCATE (stream%input_to_grid)
      IF (ALLOCATED(stream%density_signature)) &
         & DEALLOCATE (stream%density_signature)
      IF (ALLOCATED(stream%overlap_signature)) &
         & DEALLOCATE (stream%overlap_signature)
      IF (ALLOCATED(stream%inverse_coeff)) DEALLOCATE (stream%inverse_coeff)
      IF (ALLOCATED(stream%phase_roots)) DEALLOCATE (stream%phase_roots)
      IF (ALLOCATED(stream%twist_phase)) DEALLOCATE (stream%twist_phase)
      IF (ALLOCATED(stream%kfrac)) DEALLOCATE (stream%kfrac)
      IF (ALLOCATED(stream%weights)) DEALLOCATE (stream%weights)
      IF (ALLOCATED(stream%vsh)) DEALLOCATE (stream%vsh)
      IF (ALLOCATED(stream%reverse_gradient)) DEALLOCATE (stream%reverse_gradient)
      IF (ALLOCATED(stream%reverse_sigma)) DEALLOCATE (stream%reverse_sigma)
      IF (ALLOCATED(stream%g_onsfx)) DEALLOCATE (stream%g_onsfx)
      IF (ALLOCATED(stream%g_onsri)) DEALLOCATE (stream%g_onsri)
      IF (ALLOCATED(stream%dgdq_onsfx)) DEALLOCATE (stream%dgdq_onsfx)
      IF (ALLOCATED(stream%dgdq_onsri)) DEALLOCATE (stream%dgdq_onsri)
      IF (ALLOCATED(stream%density)) DEALLOCATE (stream%density)
      IF (ALLOCATED(stream%overlap)) DEALLOCATE (stream%overlap)
      IF (ALLOCATED(stream%fock)) DEALLOCATE (stream%fock)
      IF (ALLOCATED(stream%overlap_adjoint)) DEALLOCATE (stream%overlap_adjoint)
      IF (ALLOCATED(stream%transform_k)) DEALLOCATE (stream%transform_k)
      IF (ALLOCATED(stream%amat_r)) DEALLOCATE (stream%amat_r)
      IF (ALLOCATED(stream%cmat_r)) DEALLOCATE (stream%cmat_r)
      IF (ALLOCATED(stream%vmat_r)) DEALLOCATE (stream%vmat_r)
      IF (ALLOCATED(stream%amat_partner_r)) DEALLOCATE (stream%amat_partner_r)
      IF (ALLOCATED(stream%reverse_amat_r)) DEALLOCATE (stream%reverse_amat_r)
      IF (ALLOCATED(stream%reverse_vmat_r)) DEALLOCATE (stream%reverse_vmat_r)
      IF (ALLOCATED(stream%gdiagP)) DEALLOCATE (stream%gdiagP)
      IF (ALLOCATED(stream%gdiagSP)) DEALLOCATE (stream%gdiagSP)
      IF (ALLOCATED(stream%gdiagSPS)) DEALLOCATE (stream%gdiagSPS)
      IF (ALLOCATED(stream%reverse_bdiagSP)) &
         & DEALLOCATE (stream%reverse_bdiagSP)
      IF (ALLOCATED(stream%reverse_bdiagSPS)) &
         & DEALLOCATE (stream%reverse_bdiagSPS)
      stream%active = .FALSE.
      stream%applied = .FALSE.
      stream%reverse_applied = .FALSE.
      stream%reverse_transaction = .FALSE.
      stream%reverse_batch_applied = .FALSE.
      stream%reverse_result_pulled = .FALSE.
      stream%reverse_only = .FALSE.
      stream%batch_applied = .FALSE.
      stream%batch_result_pulled = .FALSE.
      stream%partial_active = .FALSE.
      stream%partial_forward_exported = .FALSE.
      stream%partial_forward_imported = .FALSE.
      stream%partial_reverse_exported = .FALSE.
      stream%partial_reverse_imported = .FALSE.
      stream%partial_image_distribution_set = .FALSE.
      stream%partial_compact_image_storage = .FALSE.
      stream%partial_pull_distribution_set = .FALSE.
      stream%partial_result_owner = .FALSE.
      stream%nk = 0
      stream%nao = 0
      stream%nspin = 0
      stream%nsh = 0
      stream%mode = 0
      stream%transform_backend = bvk_transform_dense
      stream%transform_pair_batch_size = 0
      stream%transform_pair_first = 0
      stream%transform_pair_count = 0
      stream%transform_field_count = 0
      stream%batch_size = 0
      stream%batch_first = 0
      stream%batch_count = 0
      stream%batch_last = 0
      stream%image_first = 0
      stream%partial_image_first = 0
      stream%partial_image_last = 0
      stream%plan_builds = 0
      stream%partial_generation = 0
      stream%partial_transaction_id = 0_int64
      stream%nmesh = 0
      stream%peak_complex_elements = 0_int64
      stream%kernel_signature = 0_int64
      stream%phase_signature = 0_int64
      stream%compact_phase_signature = 0_int64
      stream%partial_state_signature = 0_int64
      stream%twist = 0.0_wp
      stream%energy = 0.0_wp
   END SUBROUTINE clear_exchange_stream

!> Report whether a stream owns any complete AO-matrix Bloch or image mesh.
   LOGICAL FUNCTION cp2k_exchange_stream_has_full_mesh_storage(stream) &
      & RESULT(has_storage)
      TYPE(cp2k_exchange_stream_type), INTENT(in)        :: stream

      has_storage = ALLOCATED(stream%density) .OR. ALLOCATED(stream%overlap) &
         & .OR. ALLOCATED(stream%fock) .OR. ALLOCATED(stream%overlap_adjoint)
      IF (ALLOCATED(stream%amat_r)) has_storage = has_storage &
         & .OR. SIZE(stream%amat_r, 3) >= stream%nk
      IF (ALLOCATED(stream%amat_partner_r)) has_storage = has_storage &
         & .OR. SIZE(stream%amat_partner_r, 3) >= stream%nk
      IF (ALLOCATED(stream%cmat_r)) has_storage = has_storage &
         & .OR. SIZE(stream%cmat_r, 3) >= stream%nk
      IF (ALLOCATED(stream%vmat_r)) has_storage = has_storage &
         & .OR. SIZE(stream%vmat_r, 3) >= stream%nk
      IF (ALLOCATED(stream%transform_k)) has_storage = has_storage &
         & .OR. SIZE(stream%transform_k, 1) >= stream%nao*stream%nao
   END FUNCTION cp2k_exchange_stream_has_full_mesh_storage

!> Exact batched-mode peak number of explicitly allocated complex scalars.
!>
!> This includes both transaction-owned arrays and the largest explicit
!> push/apply/pull workspace.  It deliberately excludes immutable
!> geometry/model caches shared with other calls and compiler-created
!> expression temporaries, neither of which is owned by this transaction.
!> The query returns -1 outside an active batched transaction rather than
!> under-reporting another mode.
   INTEGER(int64) FUNCTION cp2k_exchange_stream_peak_complex_elements(stream) &
      & RESULT(nelements)
      TYPE(cp2k_exchange_stream_type), INTENT(in)        :: stream

      IF (stream%active .AND. stream%mode == cp2k_exchange_stream_batched) THEN
         nelements = stream%peak_complex_elements
      ELSE
         nelements = -1_int64
      END IF
   END FUNCTION cp2k_exchange_stream_peak_complex_elements

!> Count complex scalars currently owned by a stream object.
   INTEGER(int64) FUNCTION exchange_stream_complex_elements(stream) &
      & RESULT(nelements)
      TYPE(cp2k_exchange_stream_type), INTENT(in)        :: stream

      nelements = 0_int64
      IF (ALLOCATED(stream%density)) &
         & nelements = nelements + SIZE(stream%density, kind=int64)
      IF (ALLOCATED(stream%overlap)) &
         & nelements = nelements + SIZE(stream%overlap, kind=int64)
      IF (ALLOCATED(stream%fock)) &
         & nelements = nelements + SIZE(stream%fock, kind=int64)
      IF (ALLOCATED(stream%overlap_adjoint)) &
         & nelements = nelements + SIZE(stream%overlap_adjoint, kind=int64)
      IF (ALLOCATED(stream%transform_k)) &
         & nelements = nelements + SIZE(stream%transform_k, kind=int64)
      IF (ALLOCATED(stream%amat_r)) &
         & nelements = nelements + SIZE(stream%amat_r, kind=int64)
      IF (ALLOCATED(stream%amat_partner_r)) &
         & nelements = nelements + SIZE(stream%amat_partner_r, kind=int64)
      IF (ALLOCATED(stream%reverse_amat_r)) &
         & nelements = nelements + SIZE(stream%reverse_amat_r, kind=int64)
      IF (ALLOCATED(stream%reverse_vmat_r)) &
         & nelements = nelements + SIZE(stream%reverse_vmat_r, kind=int64)
      IF (ALLOCATED(stream%cmat_r)) &
         & nelements = nelements + SIZE(stream%cmat_r, kind=int64)
      IF (ALLOCATED(stream%vmat_r)) &
         & nelements = nelements + SIZE(stream%vmat_r, kind=int64)
      IF (ALLOCATED(stream%gdiagP)) &
         & nelements = nelements + SIZE(stream%gdiagP, kind=int64)
      IF (ALLOCATED(stream%gdiagSP)) &
         & nelements = nelements + SIZE(stream%gdiagSP, kind=int64)
      IF (ALLOCATED(stream%gdiagSPS)) &
         & nelements = nelements + SIZE(stream%gdiagSPS, kind=int64)
      IF (ALLOCATED(stream%reverse_bdiagSP)) &
         & nelements = nelements + SIZE(stream%reverse_bdiagSP, kind=int64)
      IF (ALLOCATED(stream%reverse_bdiagSPS)) &
         & nelements = nelements + SIZE(stream%reverse_bdiagSPS, kind=int64)
      IF (ALLOCATED(stream%inverse_coeff)) &
         & nelements = nelements + SIZE(stream%inverse_coeff, kind=int64)
      IF (ALLOCATED(stream%phase_roots)) &
         & nelements = nelements + SIZE(stream%phase_roots, kind=int64)
      IF (ALLOCATED(stream%twist_phase)) &
         & nelements = nelements + SIZE(stream%twist_phase, kind=int64)
   END FUNCTION exchange_stream_complex_elements

!> Mix one exact storage bit pattern into a deterministic transaction hash.
   pure INTEGER(int64) FUNCTION exchange_hash_mix(hash, bits) RESULT(mixed)
      INTEGER(int64), INTENT(in)                         :: hash, bits

      mixed = IEOR(ISHFTC(hash, 19), bits)
      mixed = IEOR(mixed, ISHFTC(bits, 41))
   END FUNCTION exchange_hash_mix

!> Fingerprint one complex AO matrix exactly at the storage-bit level.
   pure INTEGER(int64) FUNCTION exchange_hash_complex2(matrix) RESULT(hash)
      COMPLEX(wp), INTENT(in)                            :: matrix(:, :)

      INTEGER                                            :: ii, jj
      INTEGER(int64)                                     :: bits

      hash = IEOR(INT(SIZE(matrix, 1), int64), &
         & ISHFT(INT(SIZE(matrix, 2), int64), 32))
      DO jj = 1, SIZE(matrix, 2)
         DO ii = 1, SIZE(matrix, 1)
            bits = TRANSFER(REAL(matrix(ii, jj), wp), bits)
            hash = exchange_hash_mix(hash, bits)
            bits = TRANSFER(AIMAG(matrix(ii, jj)), bits)
            hash = exchange_hash_mix(hash, bits)
         END DO
      END DO
   END FUNCTION exchange_hash_complex2

!> Fingerprint one complex vector exactly at the storage-bit level.
   pure INTEGER(int64) FUNCTION exchange_hash_complex1(vector) RESULT(hash)
      COMPLEX(wp), INTENT(in)                            :: vector(:)

      INTEGER                                            :: ii
      INTEGER(int64)                                     :: bits

      hash = INT(SIZE(vector), int64)
      DO ii = 1, SIZE(vector)
         bits = TRANSFER(REAL(vector(ii), wp), bits)
         hash = exchange_hash_mix(hash, bits)
         bits = TRANSFER(AIMAG(vector(ii)), bits)
         hash = exchange_hash_mix(hash, bits)
      END DO
   END FUNCTION exchange_hash_complex1

!> Fingerprint one spin-resolved complex AO block.
   pure INTEGER(int64) FUNCTION exchange_hash_complex3(matrix) RESULT(hash)
      COMPLEX(wp), INTENT(in)                            :: matrix(:, :, :)

      INTEGER                                            :: spin

      hash = INT(SIZE(matrix, 3), int64)
      DO spin = 1, SIZE(matrix, 3)
         hash = exchange_hash_mix(hash, &
            & exchange_hash_complex2(matrix(:, :, spin)))
      END DO
   END FUNCTION exchange_hash_complex3

!> Fingerprint one real rank-three kernel array.
   pure INTEGER(int64) FUNCTION exchange_hash_real3(matrix) RESULT(hash)
      REAL(wp), INTENT(in)                               :: matrix(:, :, :)

      INTEGER                                            :: ii, jj, kk
      INTEGER(int64)                                     :: bits

      hash = IEOR(INT(SIZE(matrix, 1), int64), &
         & IEOR(ISHFT(INT(SIZE(matrix, 2), int64), 21), &
         & ISHFT(INT(SIZE(matrix, 3), int64), 42)))
      DO kk = 1, SIZE(matrix, 3)
         DO jj = 1, SIZE(matrix, 2)
            DO ii = 1, SIZE(matrix, 1)
               bits = TRANSFER(matrix(ii, jj, kk), bits)
               hash = exchange_hash_mix(hash, bits)
            END DO
         END DO
      END DO
   END FUNCTION exchange_hash_real3

!> Fingerprint the exact representative ordering of a BvK plan.
   pure INTEGER(int64) FUNCTION exchange_hash_integer2(matrix) RESULT(hash)
      INTEGER, INTENT(in)                                :: matrix(:, :)

      INTEGER                                            :: ii, jj

      hash = IEOR(INT(SIZE(matrix, 1), int64), &
         & ISHFT(INT(SIZE(matrix, 2), int64), 32))
      DO jj = 1, SIZE(matrix, 2)
         DO ii = 1, SIZE(matrix, 1)
            hash = exchange_hash_mix(hash, INT(matrix(ii, jj), int64))
         END DO
      END DO
   END FUNCTION exchange_hash_integer2

!> Fingerprint one real vector exactly at the storage-bit level.
   pure INTEGER(int64) FUNCTION exchange_hash_real1(vector) RESULT(hash)
      REAL(wp), INTENT(in)                               :: vector(:)

      INTEGER                                            :: ii
      INTEGER(int64)                                     :: bits

      hash = INT(SIZE(vector), int64)
      DO ii = 1, SIZE(vector)
         bits = TRANSFER(vector(ii), bits)
         hash = exchange_hash_mix(hash, bits)
      END DO
   END FUNCTION exchange_hash_real1

!> Fingerprint one real matrix exactly at the storage-bit level.
   pure INTEGER(int64) FUNCTION exchange_hash_real2(matrix) RESULT(hash)
      REAL(wp), INTENT(in)                               :: matrix(:, :)

      INTEGER                                            :: ii, jj
      INTEGER(int64)                                     :: bits

      hash = IEOR(INT(SIZE(matrix, 1), int64), &
         & ISHFT(INT(SIZE(matrix, 2), int64), 32))
      DO jj = 1, SIZE(matrix, 2)
         DO ii = 1, SIZE(matrix, 1)
            bits = TRANSFER(matrix(ii, jj), bits)
            hash = exchange_hash_mix(hash, bits)
         END DO
      END DO
   END FUNCTION exchange_hash_real2

!> Fingerprint one integer vector including its exact ordering.
   pure INTEGER(int64) FUNCTION exchange_hash_integer1(vector) RESULT(hash)
      INTEGER, INTENT(in)                                :: vector(:)

      INTEGER                                            :: ii

      hash = INT(SIZE(vector), int64)
      DO ii = 1, SIZE(vector)
         hash = exchange_hash_mix(hash, INT(vector(ii), int64))
      END DO
   END FUNCTION exchange_hash_integer1

!> Capture all immutable state that defines a partial k-to-R transaction.
   SUBROUTINE initialize_partial_state_signature(stream, mol)
      TYPE(cp2k_exchange_stream_type), INTENT(inout)     :: stream
      TYPE(structure_type), INTENT(in)                   :: mol

      INTEGER                                            :: periodic_integer(3)

      periodic_integer = MERGE(1, 0, mol%periodic)
      stream%partial_state_signature = 0_int64
      stream%partial_state_signature(1) = exchange_hash_real2(stream%kfrac)
      stream%partial_state_signature(2) = exchange_hash_real1(stream%weights)
      stream%partial_state_signature(3) = exchange_hash_integer2(stream%reps)
      stream%partial_state_signature(4) = &
         & exchange_hash_integer1(stream%input_to_grid)
      stream%partial_state_signature(5) = &
         & exchange_hash_complex2(stream%phase_roots)
      stream%partial_state_signature(6) = &
         & exchange_hash_complex1(stream%twist_phase)
      stream%partial_state_signature(7:8) = stream%kernel_signature
      stream%partial_state_signature(9) = exchange_hash_real3(stream%g_onsfx)
      stream%partial_state_signature(9) = exchange_hash_mix( &
         & stream%partial_state_signature(9), &
         & exchange_hash_real2(stream%g_onsri))
      stream%partial_state_signature(9) = exchange_hash_mix( &
         & stream%partial_state_signature(9), &
         & exchange_hash_real3(stream%dgdq_onsfx))
      stream%partial_state_signature(9) = exchange_hash_mix( &
         & stream%partial_state_signature(9), &
         & exchange_hash_real2(stream%dgdq_onsri))
      stream%partial_state_signature(10) = exchange_hash_real2(mol%xyz)
      stream%partial_state_signature(11) = exchange_hash_real2(mol%lattice)
      stream%partial_state_signature(12) = exchange_hash_integer1(mol%num)
      stream%partial_state_signature(12) = exchange_hash_mix( &
         & stream%partial_state_signature(12), &
         & exchange_hash_integer1(mol%id))
      stream%partial_state_signature(12) = exchange_hash_mix( &
         & stream%partial_state_signature(12), &
         & exchange_hash_integer1(periodic_integer))
      stream%partial_state_signature(12) = exchange_hash_mix( &
         & stream%partial_state_signature(12), &
         & exchange_hash_real1(stream%twist))
      stream%partial_state_signature(12) = exchange_hash_mix( &
         & stream%partial_state_signature(12), &
         & INT(stream%transform_backend, int64))
      stream%partial_state_signature(12) = exchange_hash_mix( &
         & stream%partial_state_signature(12), &
         & INT(stream%transform_pair_batch_size, int64))
   END SUBROUTINE initialize_partial_state_signature

!> Build the immutable public descriptor for the current stream generation.
   SUBROUTINE get_partial_descriptor(stream, phase, descriptor)
      TYPE(cp2k_exchange_stream_type), INTENT(in)        :: stream
      INTEGER, INTENT(in)                                :: phase
      TYPE(cp2k_exchange_partial_descriptor_type), &
         INTENT(out)                                     :: descriptor

      descriptor%abi_version = cp2k_exchange_partial_abi_version
      descriptor%scalar_bits = storage_size(0.0_wp)
      descriptor%phase = phase
      descriptor%generation = stream%partial_generation
      descriptor%transaction_id = stream%partial_transaction_id
      descriptor%nk = stream%nk
      descriptor%nao = stream%nao
      descriptor%nspin = stream%nspin
      descriptor%nsh = stream%nsh
      descriptor%nmesh = stream%nmesh
      descriptor%batch_first = stream%batch_first
      descriptor%batch_count = stream%batch_count
      descriptor%batch_capacity = stream%batch_size
      descriptor%image_first = stream%image_first
      descriptor%image_last = stream%batch_last
      descriptor%state_signature = stream%partial_state_signature
   END SUBROUTINE get_partial_descriptor

!> Return two independent storage-bit digests for buffer/descriptor binding.
   pure FUNCTION partial_descriptor_signature(descriptor) RESULT(signature)
      TYPE(cp2k_exchange_partial_descriptor_type), &
         INTENT(in)                                      :: descriptor
      INTEGER(int64)                                     :: signature(2)

      INTEGER                                            :: ii

      signature(1) = INT(descriptor%abi_version, int64)
      signature(1) = exchange_hash_mix(signature(1), &
         & INT(descriptor%scalar_bits, int64))
      signature(1) = exchange_hash_mix(signature(1), &
         & INT(descriptor%phase, int64))
      signature(1) = exchange_hash_mix(signature(1), &
         & INT(descriptor%generation, int64))
      DO ii = 1, SIZE(descriptor%transaction_id)
         signature(1) = exchange_hash_mix(signature(1), &
            & descriptor%transaction_id(ii))
      END DO
      signature(1) = exchange_hash_mix(signature(1), INT(descriptor%nk, int64))
      signature(1) = exchange_hash_mix(signature(1), INT(descriptor%nao, int64))
      signature(1) = exchange_hash_mix(signature(1), &
         & INT(descriptor%nspin, int64))
      signature(1) = exchange_hash_mix(signature(1), INT(descriptor%nsh, int64))
      DO ii = 1, 3
         signature(1) = exchange_hash_mix(signature(1), &
            & INT(descriptor%nmesh(ii), int64))
      END DO
      signature(1) = exchange_hash_mix(signature(1), &
         & INT(descriptor%batch_first, int64))
      signature(1) = exchange_hash_mix(signature(1), &
         & INT(descriptor%batch_count, int64))
      signature(1) = exchange_hash_mix(signature(1), &
         & INT(descriptor%batch_capacity, int64))
      signature(1) = exchange_hash_mix(signature(1), &
         & INT(descriptor%image_first, int64))
      signature(1) = exchange_hash_mix(signature(1), &
         & INT(descriptor%image_last, int64))
      signature(2) = NOT(signature(1))
      DO ii = 1, SIZE(descriptor%state_signature)
         signature(1) = exchange_hash_mix(signature(1), &
            & descriptor%state_signature(ii))
         signature(2) = exchange_hash_mix(signature(2), &
            & ISHFTC(descriptor%state_signature(ii), 17))
      END DO
   END FUNCTION partial_descriptor_signature

!> Compare every descriptor field rather than relying only on its digest.
   pure LOGICAL FUNCTION cp2k_exchange_partial_descriptors_match(left, right) &
      & RESULT(matches)
      TYPE(cp2k_exchange_partial_descriptor_type), &
         INTENT(in)                                      :: left, right

      matches = left%abi_version == right%abi_version &
         & .AND. left%scalar_bits == right%scalar_bits &
         & .AND. left%phase == right%phase &
         & .AND. left%generation == right%generation &
         & .AND. ALL(left%transaction_id == right%transaction_id) &
         & .AND. left%nk == right%nk &
         & .AND. left%nao == right%nao &
         & .AND. left%nspin == right%nspin &
         & .AND. left%nsh == right%nsh &
         & .AND. ALL(left%nmesh == right%nmesh) &
         & .AND. left%batch_first == right%batch_first &
         & .AND. left%batch_count == right%batch_count &
         & .AND. left%batch_capacity == right%batch_capacity &
         & .AND. left%image_first == right%image_first &
         & .AND. left%image_last == right%image_last &
         & .AND. ALL(left%state_signature == right%state_signature)
   END FUNCTION cp2k_exchange_partial_descriptors_match

!> Multiply nonnegative factors without exceeding a caller-supplied limit.
   pure SUBROUTINE checked_partial_product(factors, limit, value, ok)
      INTEGER(int64), INTENT(in)                         :: factors(:), limit
      INTEGER(int64), INTENT(out)                        :: value
      LOGICAL, INTENT(out)                               :: ok

      INTEGER                                            :: ii

      ok = .FALSE.
      value = 0_int64
      IF (limit < 0_int64 .OR. ANY(factors < 0_int64)) RETURN
      value = 1_int64
      DO ii = 1, SIZE(factors)
         IF (factors(ii) == 0_int64) THEN
            value = 0_int64
            ok = .TRUE.
            RETURN
         END IF
         IF (value > limit/factors(ii)) THEN
            value = 0_int64
            RETURN
         END IF
         value = value*factors(ii)
      END DO
      ok = .TRUE.
   END SUBROUTINE checked_partial_product

!> Check a regular-mesh extent without overflowing default-integer PRODUCT.
   pure LOGICAL FUNCTION complete_regular_mesh(nmesh, nk) RESULT(complete)
      INTEGER, INTENT(in)                                :: nmesh(3), nk

      INTEGER(int64)                                     :: mesh_size
      LOGICAL                                            :: product_ok

      complete = .FALSE.
      IF (nk < 1 .OR. ANY(nmesh < 1)) RETURN
      CALL checked_partial_product(INT(nmesh, int64), INT(nk, int64), &
         & mesh_size, product_ok)
      IF (.NOT. product_ok) RETURN
      complete = mesh_size == INT(nk, int64)
   END FUNCTION complete_regular_mesh

!> Return forward block/total extents that fit default-integer indexing.
   pure SUBROUTINE partial_forward_extents(nao, batch_count, nspin, nblock, total, &
      & ok)
      INTEGER, INTENT(in)                                :: nao, batch_count, nspin
      INTEGER(int64), INTENT(out)                        :: nblock, total
      LOGICAL, INTENT(out)                               :: ok

      INTEGER(int64), PARAMETER :: index_limit = INT(HUGE(0), int64)

      CALL checked_partial_product([INT(nao, int64), INT(nao, int64), &
         & INT(batch_count, int64), INT(nspin, int64)], index_limit, nblock, ok)
      total = 0_int64
      IF (.NOT. ok) RETURN
      IF (nblock > index_limit/4_int64) THEN
         ok = .FALSE.
         RETURN
      END IF
      total = 4_int64*nblock
   END SUBROUTINE partial_forward_extents

!> Return reverse block/onsite/total extents within default-integer indexing.
   pure SUBROUTINE partial_reverse_extents(nao, batch_count, nspin, nblock, &
      & nonsite, total, ok)
      INTEGER, INTENT(in)                                :: nao, batch_count, nspin
      INTEGER(int64), INTENT(out)                        :: nblock, nonsite, total
      LOGICAL, INTENT(out)                               :: ok

      INTEGER(int64), PARAMETER :: index_limit = INT(HUGE(0), int64)

      INTEGER(int64)                                     :: block_total, onsite_total

      CALL checked_partial_product([INT(nao, int64), INT(nao, int64), &
         & INT(batch_count, int64), INT(nspin, int64)], index_limit, nblock, ok)
      nonsite = 0_int64
      total = 0_int64
      IF (.NOT. ok) RETURN
      CALL checked_partial_product([INT(nao, int64), INT(nspin, int64)], &
         & index_limit, nonsite, ok)
      IF (.NOT. ok) RETURN
      IF (nblock > index_limit/6_int64 &
         & .OR. nonsite > index_limit/5_int64) THEN
         ok = .FALSE.
         RETURN
      END IF
      block_total = 6_int64*nblock
      onsite_total = 5_int64*nonsite
      IF (block_total > index_limit - onsite_total) THEN
         ok = .FALSE.
         RETURN
      END IF
      total = block_total + onsite_total
   END SUBROUTINE partial_reverse_extents

!> Validate one contiguous subset of the descriptor's current image batch.
   pure SUBROUTINE partial_payload_range(descriptor, image_first, image_count, &
      & image_offset, ok)
      TYPE(cp2k_exchange_partial_descriptor_type), &
         INTENT(in)                                      :: descriptor
      INTEGER, INTENT(in)                                :: image_first, image_count
      INTEGER, INTENT(out)                               :: image_offset
      LOGICAL, INTENT(out)                               :: ok

      image_offset = 0
      ok = .FALSE.
      IF (image_count < 0) RETURN
      IF (image_count == 0) THEN
         ok = image_first == 0
         RETURN
      END IF
      IF (image_first < descriptor%batch_first) RETURN
      image_offset = image_first - descriptor%batch_first
      IF (image_offset >= descriptor%batch_count) RETURN
      IF (image_count > descriptor%batch_count - image_offset) RETURN
      ok = .TRUE.
   END SUBROUTINE partial_payload_range

!> Validate all structural invariants of a canonical partial descriptor.
   pure LOGICAL FUNCTION partial_descriptor_is_valid(descriptor) RESULT(valid)
      TYPE(cp2k_exchange_partial_descriptor_type), &
         INTENT(in)                                      :: descriptor

      INTEGER(int64)                                     :: mesh_size
      LOGICAL                                            :: product_ok

      valid = .FALSE.
      IF (descriptor%abi_version /= cp2k_exchange_partial_abi_version) RETURN
      IF (descriptor%scalar_bits /= storage_size(0.0_wp)) RETURN
      IF (descriptor%phase /= 1 .AND. descriptor%phase /= 2) RETURN
      IF (descriptor%generation < 1) RETURN
      IF (ALL(descriptor%transaction_id == 0_int64)) RETURN
      IF (descriptor%nk < 1 .OR. descriptor%nao < 1 &
         & .OR. descriptor%nspin < 1 .OR. descriptor%nsh < 1) RETURN
      IF (ANY(descriptor%nmesh < 1)) RETURN
      CALL checked_partial_product(INT(descriptor%nmesh, int64), &
         & INT(descriptor%nk, int64), mesh_size, product_ok)
      IF (.NOT. product_ok) RETURN
      IF (mesh_size /= INT(descriptor%nk, int64)) RETURN
      IF (descriptor%image_first < 1 &
         & .OR. descriptor%image_last > descriptor%nk &
         & .OR. descriptor%image_first > descriptor%image_last) RETURN
      IF (descriptor%batch_capacity < 1 &
         & .OR. descriptor%batch_first < descriptor%image_first &
         & .OR. descriptor%batch_first > descriptor%image_last) RETURN
      IF (descriptor%batch_count /= MIN(descriptor%batch_capacity, &
         & descriptor%image_last - descriptor%batch_first + 1)) RETURN
      valid = .TRUE.
   END FUNCTION partial_descriptor_is_valid

!> Number of int64 words in the portable canonical descriptor encoding.
   INTEGER(int64) FUNCTION cp2k_exchange_partial_descriptor_size() &
      & RESULT(nwords)
      nwords = INT(cp2k_exchange_partial_descriptor_words, int64)
   END FUNCTION cp2k_exchange_partial_descriptor_size

!> Pack a descriptor into a portable, padding-free int64 representation.
!>
!> CP2K can allgather these words, unpack every rank's descriptor, and compare
!> every field with `cp2k_exchange_partial_descriptors_match`.  The internal
!> two-word payload digest is never sufficient for accepting a collective.
   SUBROUTINE cp2k_exchange_partial_descriptor_pack(descriptor, values, error)
      TYPE(cp2k_exchange_partial_descriptor_type), &
         INTENT(in)                                      :: descriptor
      INTEGER(int64), ALLOCATABLE, INTENT(out)           :: values(:)
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      IF (.NOT. partial_descriptor_is_valid(descriptor)) THEN
         CALL fatal_error(error, "invalid g-xTB exchange partial descriptor")
         RETURN
      END IF
      ALLOCATE (values(cp2k_exchange_partial_descriptor_words), source=0_int64)
      values(1) = INT(descriptor%abi_version, int64)
      values(2) = INT(descriptor%scalar_bits, int64)
      values(3) = INT(descriptor%phase, int64)
      values(4) = INT(descriptor%generation, int64)
      values(5:6) = descriptor%transaction_id
      values(7) = INT(descriptor%nk, int64)
      values(8) = INT(descriptor%nao, int64)
      values(9) = INT(descriptor%nspin, int64)
      values(10) = INT(descriptor%nsh, int64)
      values(11:13) = INT(descriptor%nmesh, int64)
      values(14) = INT(descriptor%batch_first, int64)
      values(15) = INT(descriptor%batch_count, int64)
      values(16) = INT(descriptor%batch_capacity, int64)
      values(17) = INT(descriptor%image_first, int64)
      values(18) = INT(descriptor%image_last, int64)
      values(19:30) = descriptor%state_signature
   END SUBROUTINE cp2k_exchange_partial_descriptor_pack

!> Unpack and structurally validate a portable canonical descriptor.
   SUBROUTINE cp2k_exchange_partial_descriptor_unpack(descriptor, values, error)
      TYPE(cp2k_exchange_partial_descriptor_type), &
         INTENT(out)                                     :: descriptor
      INTEGER(int64), INTENT(in)                         :: values(:)
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      INTEGER(int64), PARAMETER :: default_max = INT(HUGE(0), int64)

      INTEGER                                            :: ii

      IF (SIZE(values) /= cp2k_exchange_partial_descriptor_words) THEN
         CALL fatal_error(error, &
            & "inconsistent g-xTB exchange partial descriptor size")
         RETURN
      END IF
      DO ii = 1, 18
         IF (ii == 5 .OR. ii == 6) CYCLE
         IF (values(ii) > default_max .OR. values(ii) < -default_max) THEN
            CALL fatal_error(error, &
               & "out-of-range g-xTB exchange partial descriptor field")
            RETURN
         END IF
      END DO
      descriptor%abi_version = INT(values(1))
      descriptor%scalar_bits = INT(values(2))
      descriptor%phase = INT(values(3))
      descriptor%generation = INT(values(4))
      descriptor%transaction_id = values(5:6)
      descriptor%nk = INT(values(7))
      descriptor%nao = INT(values(8))
      descriptor%nspin = INT(values(9))
      descriptor%nsh = INT(values(10))
      descriptor%nmesh = INT(values(11:13))
      descriptor%batch_first = INT(values(14))
      descriptor%batch_count = INT(values(15))
      descriptor%batch_capacity = INT(values(16))
      descriptor%image_first = INT(values(17))
      descriptor%image_last = INT(values(18))
      descriptor%state_signature = values(19:30)
      IF (.NOT. partial_descriptor_is_valid(descriptor)) THEN
         CALL fatal_error(error, "invalid g-xTB exchange partial descriptor")
         RETURN
      END IF
   END SUBROUTINE cp2k_exchange_partial_descriptor_unpack

!> Reconstruct and Hermitize one Fock contribution from response images.
!>
!> The same projection/scaling kernel serves both a complete reduced image
!> mesh and one additive bounded batch.  `add_onsite` is true exactly for the
!> complete mesh or for the unique batch containing the origin image.
   SUBROUTINE reconstruct_stream_fock_block(overlap, phase_forward, add_onsite, &
      & amat_r, cmat_r, vmat_r, gdiagP, gdiagSP, gdiagSPS, fock)
      COMPLEX(wp), INTENT(in)                            :: overlap(:, :), phase_forward(:)
      LOGICAL, INTENT(in)                                :: add_onsite
      COMPLEX(wp), INTENT(in)                            :: amat_r(:, :, :, :), cmat_r(:, :, :, :), &
                                                            vmat_r(:, :, :, :), gdiagP(:, :), &
                                                            gdiagSP(:, :), gdiagSPS(:, :)
      COMPLEX(wp), INTENT(out)                           :: fock(:, :, :)

      COMPLEX(wp)                                        :: tmp
      COMPLEX(wp), ALLOCATABLE                           :: amat(:, :), cmat(:, :), vmat(:, :), &
                                                            work(:, :)
      INTEGER                                            :: iao, icell, ii, jj, spin
      REAL(wp)                                           :: spin_factor

      ALLOCATE (amat(SIZE(overlap, 1), SIZE(overlap, 2)), &
         & cmat(SIZE(overlap, 1), SIZE(overlap, 2)), &
         & vmat(SIZE(overlap, 1), SIZE(overlap, 2)), &
         & work(SIZE(overlap, 1), SIZE(overlap, 2)))
      spin_factor = 0.5_wp
      IF (SIZE(fock, 3) > 1) spin_factor = 1.0_wp
      DO spin = 1, SIZE(fock, 3)
         amat = (0.0_wp, 0.0_wp)
         cmat = (0.0_wp, 0.0_wp)
         vmat = (0.0_wp, 0.0_wp)
         DO icell = 1, SIZE(phase_forward)
            amat = amat + phase_forward(icell)*amat_r(:, :, icell, spin)
            cmat = cmat + phase_forward(icell)*cmat_r(:, :, icell, spin)
            vmat = vmat + phase_forward(icell)*vmat_r(:, :, icell, spin)
         END DO
         work = amat + 0.5_wp*MATMUL(overlap, vmat)
         IF (add_onsite) THEN
            DO iao = 1, SIZE(overlap, 1)
               work(:, iao) = work(:, iao) &
                  & + 0.25_wp*gdiagP(iao, spin)*overlap(:, iao)
            END DO
         END IF
         fock(:, :, spin) = cmat + MATMUL(work, overlap)
         IF (add_onsite) THEN
            DO iao = 1, SIZE(overlap, 1)
               fock(:, iao, spin) = fock(:, iao, spin) &
                  & + 0.5_wp*gdiagSP(iao, spin)*overlap(:, iao)
               fock(iao, iao, spin) = fock(iao, iao, spin) &
                  & + 0.25_wp*gdiagSPS(iao, spin)
            END DO
         END IF
         DO ii = 1, SIZE(overlap, 1)
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
   END SUBROUTINE reconstruct_stream_fock_block

!> Check the exact charge-dependent onsite state captured by `begin`.
   LOGICAL FUNCTION exchange_stream_onsite_matches(stream, cache) RESULT(matches)
      TYPE(cp2k_exchange_stream_type), INTENT(in)        :: stream
      TYPE(exchange_cache), INTENT(in)                   :: cache

      matches = .FALSE.
      IF (.NOT. ALLOCATED(stream%g_onsfx) &
         & .OR. .NOT. ALLOCATED(stream%g_onsri) &
         & .OR. .NOT. ALLOCATED(stream%dgdq_onsfx) &
         & .OR. .NOT. ALLOCATED(stream%dgdq_onsri) &
         & .OR. .NOT. ALLOCATED(cache%g_onsfx) &
         & .OR. .NOT. ALLOCATED(cache%g_onsri) &
         & .OR. .NOT. ALLOCATED(cache%dgdq_onsfx) &
         & .OR. .NOT. ALLOCATED(cache%dgdq_onsri)) RETURN
      IF (ANY(SHAPE(stream%g_onsfx) /= SHAPE(cache%g_onsfx)) &
         & .OR. ANY(SHAPE(stream%g_onsri) /= SHAPE(cache%g_onsri)) &
         & .OR. ANY(SHAPE(stream%dgdq_onsfx) /= SHAPE(cache%dgdq_onsfx)) &
         & .OR. ANY(SHAPE(stream%dgdq_onsri) /= SHAPE(cache%dgdq_onsri))) RETURN
      IF (ANY(stream%g_onsfx /= cache%g_onsfx) &
         & .OR. ANY(stream%g_onsri /= cache%g_onsri) &
         & .OR. ANY(stream%dgdq_onsfx /= cache%dgdq_onsfx) &
         & .OR. ANY(stream%dgdq_onsri /= cache%dgdq_onsri)) RETURN
      matches = .TRUE.
   END FUNCTION exchange_stream_onsite_matches

!> Validate the exact geometry, plan generation, model, and onsite snapshot.
   LOGICAL FUNCTION exchange_stream_state_matches(stream, calc, mol, ecache, plan) &
      & RESULT(matches)
      TYPE(cp2k_exchange_stream_type), INTENT(in)        :: stream
      TYPE(xtb_calculator), INTENT(in)                   :: calc
      TYPE(structure_type), INTENT(in)                   :: mol
      TYPE(container_cache), INTENT(in)                  :: ecache
      TYPE(cp2k_stream_plan_type), INTENT(in) :: plan

      matches = .FALSE.
      IF (.NOT. ALLOCATED(calc%exchange) .OR. .NOT. ALLOCATED(ecache%raw)) RETURN
      SELECT TYPE (cache => ecache%raw)
      TYPE is (exchange_cache)
         SELECT TYPE (exchange => calc%exchange)
         TYPE is (exchange_fock)
            IF (exchange%nao /= stream%nao .OR. exchange%nsh /= stream%nsh) RETURN
            IF (plan%bvk_plan_builds /= stream%plan_builds) RETURN
            IF (.NOT. plan%bvk_matches(mol, stream%nmesh, stream%kfrac, &
               & stream%weights, require_phases= &
               & stream%mode /= cp2k_exchange_stream_batched)) RETURN
            IF (.NOT. plan%model_matches(exchange)) RETURN
            IF (.NOT. exchange_stream_onsite_matches(stream, cache)) RETURN
            IF (.NOT. ALLOCATED(stream%input_to_grid) &
               & .OR. .NOT. ALLOCATED(plan%bvk_input_to_grid)) RETURN
            IF (SIZE(stream%input_to_grid) /= &
               & SIZE(plan%bvk_input_to_grid)) RETURN
            IF (ANY(stream%input_to_grid /= plan%bvk_input_to_grid)) RETURN
            IF (ANY(stream%twist /= plan%bvk_twist)) RETURN
            IF (ANY(SHAPE(stream%reps) /= SHAPE(plan%bvk_kernel%reps))) RETURN
            IF (ANY(stream%reps /= plan%bvk_kernel%reps)) RETURN
            IF (stream%kernel_signature(1) /= &
               & exchange_hash_real3(plan%bvk_kernel%g_mulliken_r)) RETURN
            IF (stream%kernel_signature(2) /= IEOR( &
               & exchange_hash_real3(plan%bvk_kernel%g_bocorr_r), &
               & exchange_hash_integer2(plan%bvk_kernel%reps))) RETURN
            IF (stream%mode /= cp2k_exchange_stream_batched) THEN
               IF (stream%phase_signature(1) /= &
                  & exchange_hash_complex2(plan%bvk_phase_forward)) RETURN
               IF (stream%phase_signature(2) /= &
                  & exchange_hash_complex2(plan%bvk_phase_inverse)) RETURN
            ELSE
               IF (.NOT. ALLOCATED(plan%bvk_phase_roots) &
                  & .OR. .NOT. ALLOCATED(plan%bvk_twist_phase)) RETURN
               IF (stream%compact_phase_signature(1) /= &
                  & exchange_hash_complex2(plan%bvk_phase_roots)) RETURN
               IF (stream%compact_phase_signature(2) /= &
                  & exchange_hash_complex1(plan%bvk_twist_phase)) RETURN
            END IF
            matches = .TRUE.
         END SELECT
      END SELECT
   END FUNCTION exchange_stream_state_matches

!> Begin collecting the unweighted Bloch blocks of one complete regular mesh.
!>
!> Only one transaction may be active in a stream object.  The expensive BvK
!> plan is prepared here, before blocks are accepted, so that `apply` can also
!> detect an intervening geometry, mesh, or model change.  The optional mode
!> selects the image-mesh reduced path, the full-mesh oracle, or the strictly
!> sequential bounded image-batch path; reduced mode is the default.  Batched
!> mode additionally requires a positive `batch_size`.  Its mandatory caller
!> sequence is `batch_bounds`, one `push` for every k point, `batch_apply`, one
!> `batch_pull` for every k point, and `batch_advance`; repeat until `done`,
!> then obtain exactly one `batch_result` and call `end`.  The caller owns and
!> zeroes the final Fock accumulators before the first batch.  The optional
!> `image_first,image_count` pair restricts a batched transaction to one
!> additive, contiguous image range.  Disjoint ranges may therefore be
!> evaluated concurrently by independent stream objects or MPI ranks and their
!> energy, shell potential, and Fock contributions summed by the caller.  A
!> complete result requires an exact, nonoverlapping cover of images 1..Nk.
!> Choosing `batch_size >= image_count` gives one K-block submission per
!> independent range.  A serial range larger than its bounded batch must still
!> resubmit K blocks: without retaining a full Bloch/image AO mesh, released
!> input information cannot contribute to a later image sum.
!> `reverse_only=.true.` starts a derivative-only transaction in reduced or
!> batched mode.  It deliberately skips the forward Fock/result protocol; the
!> caller immediately enters `reverse_push`/`reverse_batch_apply`/blockwise
!> `reverse_batch_pull`.  This is the preferred force/stress path when the SCF
!> Fock matrix has already been built elsewhere.  No concurrent use of one
!> stream object is supported.
   SUBROUTINE cp2k_exchange_stream_begin(stream, calc, mol, ecache, plan, wfn, nmesh, &
      & kfrac, weights, error, mode, batch_size, image_first, image_count, &
      & reverse_only, transform_backend, transform_pair_batch_size, &
      & local_image_first, local_image_count)
      TYPE(cp2k_exchange_stream_type), INTENT(inout)     :: stream
      TYPE(xtb_calculator), INTENT(inout)                :: calc
      TYPE(structure_type), INTENT(in)                   :: mol
      TYPE(container_cache), INTENT(inout)               :: ecache
      TYPE(cp2k_stream_plan_type), INTENT(inout) :: plan
      TYPE(wavefunction_type), INTENT(in)                :: wfn
      INTEGER, INTENT(in)                                :: nmesh(3)
      REAL(wp), INTENT(in)                               :: kfrac(:, :), weights(:)
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error
      INTEGER, INTENT(in), OPTIONAL                      :: mode, batch_size, image_first, &
                                                            image_count
      LOGICAL, INTENT(in), OPTIONAL                      :: reverse_only
      INTEGER, INTENT(in), OPTIONAL                      :: transform_backend, &
                                                            transform_pair_batch_size, &
                                                            local_image_first, local_image_count

      INTEGER :: grid(3), icell, idim, igrid, ik_oracle, inverse_cell, inverse_grid(3), nk, npair, &
         selected_batch_size, selected_image_count, selected_image_first, selected_image_last, &
         selected_mode, selected_pair_batch_size, selected_storage_capacity, &
         selected_transform_backend
      INTEGER, ALLOCATABLE                               :: rep_to_cell(:)
      LOGICAL                                            :: selected_reverse_only
      REAL(wp)                                           :: canonical, periodic_delta

      IF (stream%active) THEN
         CALL fatal_error(error, "g-xTB exchange stream is already active")
         RETURN
      END IF

      selected_mode = cp2k_exchange_stream_reduced
      IF (PRESENT(mode)) selected_mode = mode
      selected_transform_backend = bvk_transform_dense
      IF (PRESENT(transform_backend)) &
         & selected_transform_backend = transform_backend
      IF (.NOT. valid_exchange_transform_backend(selected_transform_backend)) THEN
         CALL fatal_error(error, "unknown streamed g-xTB transform backend")
         RETURN
      END IF
      selected_reverse_only = .FALSE.
      IF (PRESENT(reverse_only)) selected_reverse_only = reverse_only
      IF (selected_mode /= cp2k_exchange_stream_reduced &
         & .AND. selected_mode /= cp2k_exchange_stream_oracle &
         & .AND. selected_mode /= cp2k_exchange_stream_batched) THEN
         CALL fatal_error(error, "unknown g-xTB exchange stream mode")
         RETURN
      END IF
      IF (selected_reverse_only .AND. &
         & selected_mode == cp2k_exchange_stream_oracle) THEN
         CALL fatal_error(error, &
            & "reverse-only g-xTB exchange requires reduced or batched mode")
         RETURN
      END IF
      IF (selected_transform_backend /= bvk_transform_dense .AND. &
         & selected_mode /= cp2k_exchange_stream_batched) THEN
         CALL fatal_error(error, &
            & "compact streamed g-xTB transforms require batched mode")
         RETURN
      END IF
      selected_batch_size = 0
      IF (selected_mode == cp2k_exchange_stream_batched) THEN
         IF (.NOT. PRESENT(batch_size)) THEN
            CALL fatal_error(error, &
               & "batched g-xTB exchange stream requires a batch size")
            RETURN
         END IF
         IF (batch_size < 1) THEN
            CALL fatal_error(error, &
               & "g-xTB exchange stream batch size must be positive")
            RETURN
         END IF
         selected_batch_size = batch_size
      ELSE IF (PRESENT(batch_size)) THEN
         CALL fatal_error(error, &
            & "g-xTB exchange stream batch size requires batched mode")
         RETURN
      END IF
      selected_pair_batch_size = 0
      IF (selected_transform_backend /= bvk_transform_dense) THEN
         IF (PRESENT(transform_pair_batch_size)) THEN
            IF (transform_pair_batch_size < 1) THEN
               CALL fatal_error(error, &
                  & "g-xTB transform AO-pair batch size must be positive")
               RETURN
            END IF
            selected_pair_batch_size = transform_pair_batch_size
         ELSE
            selected_pair_batch_size = MAX(1, selected_batch_size)
         END IF
      END IF

      nk = SIZE(weights)
      selected_image_first = 1
      selected_image_count = nk
      selected_image_last = nk
      IF (selected_mode == cp2k_exchange_stream_batched) THEN
         IF (PRESENT(image_first) .NEQV. PRESENT(image_count)) THEN
            CALL fatal_error(error, &
               & "g-xTB exchange image range requires first and count")
            RETURN
         END IF
         IF (PRESENT(image_first)) THEN
            IF (image_first < 1 .OR. image_first > nk) THEN
               CALL fatal_error(error, "g-xTB exchange image range is invalid")
               RETURN
            END IF
            IF (image_count < 1 .OR. image_count > nk - image_first + 1) THEN
               CALL fatal_error(error, "g-xTB exchange image range is invalid")
               RETURN
            END IF
            selected_image_first = image_first
            selected_image_count = image_count
            selected_image_last = image_first + (image_count - 1)
         END IF
      ELSE IF (PRESENT(image_first) .OR. PRESENT(image_count)) THEN
         CALL fatal_error(error, &
            & "g-xTB exchange image range requires batched mode")
         RETURN
      END IF
      selected_storage_capacity = MIN(selected_batch_size, selected_image_count)
      IF (PRESENT(local_image_first) .NEQV. PRESENT(local_image_count)) THEN
         CALL fatal_error(error, &
            & "local g-xTB image storage requires first and count")
         RETURN
      END IF
      IF (PRESENT(local_image_first)) THEN
         IF (selected_mode /= cp2k_exchange_stream_batched) THEN
            CALL fatal_error(error, &
               & "local g-xTB image storage requires batched mode")
            RETURN
         END IF
         IF (local_image_count < 0) THEN
            CALL fatal_error(error, "negative local g-xTB image count")
            RETURN
         END IF
         IF (local_image_count == 0) THEN
            IF (local_image_first /= 0) THEN
               CALL fatal_error(error, &
                  & "empty local g-xTB image storage must start at zero")
               RETURN
            END IF
         ELSE
            IF (local_image_first < selected_image_first &
               & .OR. local_image_first > selected_image_last &
               & .OR. local_image_count > &
               & selected_image_last - local_image_first + 1) THEN
               CALL fatal_error(error, "local g-xTB image storage is invalid")
               RETURN
            END IF
         END IF
         selected_storage_capacity = MIN(selected_batch_size, local_image_count)
      END IF
      IF (.NOT. complete_regular_mesh(nmesh, nk)) THEN
         CALL fatal_error(error, &
            & "g-xTB exchange stream requires a complete regular k mesh")
         RETURN
      END IF
      DO idim = 1, 3
         IF (.NOT. mol%periodic(idim) .AND. nmesh(idim) /= 1) THEN
            CALL fatal_error(error, &
               & "nonperiodic directions cannot be replicated by the g-xTB k mesh")
            RETURN
         END IF
      END DO
      IF (SIZE(kfrac, 1) /= 3 .OR. SIZE(kfrac, 2) /= nk) THEN
         CALL fatal_error(error, &
            & "inconsistent g-xTB exchange stream mesh dimensions")
         RETURN
      END IF
      IF (.NOT. ALL(ieee_is_finite(kfrac))) THEN
         CALL fatal_error(error, "g-xTB exchange stream k points must be finite")
         RETURN
      END IF
      IF (.NOT. ALL(ieee_is_finite(weights)) .OR. &
         & ABS(SUM(weights) - 1.0_wp) > 1.0e-12_wp &
         & .OR. MAXVAL(ABS(weights - 1.0_wp/REAL(nk, wp))) > 1.0e-12_wp) THEN
         CALL fatal_error(error, &
            & "g-xTB exchange stream requires uniform k-point weights")
         RETURN
      END IF

      CALL cp2k_prepare_exchange(calc, mol, ecache, wfn, error)
      IF (ALLOCATED(error)) RETURN
      SELECT TYPE (cache => ecache%raw)
      TYPE is (exchange_cache)
         SELECT TYPE (exchange => calc%exchange)
         TYPE is (exchange_fock)
            IF (wfn%nspin /= 1 .AND. wfn%nspin /= 2) THEN
               CALL fatal_error(error, &
                  & "g-xTB exchange stream requires one or two spin channels")
               RETURN
            END IF
            CALL prepare_bvk_plan(exchange, mol, plan, nmesh, kfrac, weights, &
               & error, retain_dense_phases= &
               & selected_mode /= cp2k_exchange_stream_batched)
            IF (ALLOCATED(error)) RETURN
            IF (selected_mode == cp2k_exchange_stream_batched) THEN
               ! The inverse-image energy collapse is exact only for characters
               ! of a uniform finite Fourier group.  The general plan accepts a
               ! looser geometric tolerance, but bounded batching deliberately
               ! rejects anything beyond roundoff instead of silently dropping
               ! small off-pair Gram coefficients.
               IF (MAXVAL(ABS(weights - 1.0_wp/REAL(nk, wp))) > &
                  & 64.0_wp*EPSILON(1.0_wp)) THEN
                  CALL fatal_error(error, &
                     & "batched g-xTB exchange requires machine-uniform weights")
                  RETURN
               END IF
               DO ik_oracle = 1, nk
                  igrid = plan%bvk_input_to_grid(ik_oracle) - 1
                  grid(1) = MODULO(igrid, nmesh(1))
                  grid(2) = MODULO(igrid/nmesh(1), nmesh(2))
                  grid(3) = igrid/(nmesh(1)*nmesh(2))
                  DO idim = 1, 3
                     canonical = plan%bvk_twist(idim) &
                        & + REAL(grid(idim), wp)/REAL(nmesh(idim), wp)
                     periodic_delta = MODULO(kfrac(idim, ik_oracle), 1.0_wp) &
                        & - MODULO(canonical, 1.0_wp)
                     periodic_delta = periodic_delta - NINT(periodic_delta)
                     IF (ABS(periodic_delta) > &
                        & 64.0_wp*EPSILON(1.0_wp)) THEN
                        CALL fatal_error(error, &
                           & "batched g-xTB exchange requires a machine-regular grid")
                        RETURN
                     END IF
                  END DO
               END DO
            END IF
            IF (.NOT. ALLOCATED(cache%g_onsfx) &
               & .OR. .NOT. ALLOCATED(cache%g_onsri) &
               & .OR. .NOT. ALLOCATED(cache%dgdq_onsfx) &
               & .OR. .NOT. ALLOCATED(cache%dgdq_onsri)) THEN
               CALL fatal_error(error, &
                  & "g-xTB exchange stream onsite state is not prepared")
               RETURN
            END IF

            CALL clear_exchange_stream(stream)
            stream%nk = nk
            stream%nao = exchange%nao
            stream%nspin = wfn%nspin
            stream%nsh = exchange%nsh
            stream%mode = selected_mode
            stream%transform_backend = selected_transform_backend
            stream%reverse_only = selected_reverse_only
            stream%plan_builds = plan%bvk_plan_builds
            stream%nmesh = nmesh
            ALLOCATE (stream%seen(nk), source=.FALSE.)
            ALLOCATE (stream%forward_pulled(nk), source=.FALSE.)
            ALLOCATE (stream%reverse_pulled(nk), source=.FALSE.)
            ALLOCATE (stream%batch_pulled(nk), source=.FALSE.)
            ALLOCATE (stream%reps, source=plan%bvk_kernel%reps)
            ALLOCATE (stream%input_to_grid, source=plan%bvk_input_to_grid)
            ALLOCATE (stream%kfrac(3, nk), source=kfrac)
            ALLOCATE (stream%weights(nk), source=weights)
            stream%twist = plan%bvk_twist
            stream%kernel_signature(1) = &
               & exchange_hash_real3(plan%bvk_kernel%g_mulliken_r)
            stream%kernel_signature(2) = IEOR( &
               & exchange_hash_real3(plan%bvk_kernel%g_bocorr_r), &
               & exchange_hash_integer2(plan%bvk_kernel%reps))
            IF (selected_mode /= cp2k_exchange_stream_batched) THEN
               stream%phase_signature(1) = &
                  & exchange_hash_complex2(plan%bvk_phase_forward)
               stream%phase_signature(2) = &
                  & exchange_hash_complex2(plan%bvk_phase_inverse)
            ELSE
               ALLOCATE (stream%phase_roots, source=plan%bvk_phase_roots)
               ALLOCATE (stream%twist_phase, source=plan%bvk_twist_phase)
               stream%compact_phase_signature(1) = &
                  & exchange_hash_complex2(plan%bvk_phase_roots)
               stream%compact_phase_signature(2) = &
                  & exchange_hash_complex1(plan%bvk_twist_phase)
            END IF
            ALLOCATE (stream%g_onsfx, source=cache%g_onsfx)
            ALLOCATE (stream%g_onsri, source=cache%g_onsri)
            ALLOCATE (stream%dgdq_onsfx, source=cache%dgdq_onsfx)
            ALLOCATE (stream%dgdq_onsri, source=cache%dgdq_onsri)
            IF (selected_mode == cp2k_exchange_stream_oracle) THEN
               ALLOCATE (stream%density(exchange%nao, exchange%nao, wfn%nspin, nk), &
                  & stream%overlap(exchange%nao, exchange%nao, nk), &
                  & stream%fock(exchange%nao, exchange%nao, wfn%nspin, nk), &
                  & source=(0.0_wp, 0.0_wp))
               ALLOCATE (stream%overlap_adjoint(exchange%nao, exchange%nao, nk), &
                  & source=(0.0_wp, 0.0_wp))
            ELSE IF (selected_mode == cp2k_exchange_stream_reduced) THEN
               ALLOCATE (stream%amat_r(exchange%nao, exchange%nao, nk, wfn%nspin), &
                  & stream%cmat_r(exchange%nao, exchange%nao, nk, wfn%nspin), &
                  & stream%vmat_r(exchange%nao, exchange%nao, nk, wfn%nspin), &
                  & source=(0.0_wp, 0.0_wp))
               ALLOCATE (stream%gdiagP(exchange%nao, wfn%nspin), &
                  & stream%gdiagSP(exchange%nao, wfn%nspin), &
                  & stream%gdiagSPS(exchange%nao, wfn%nspin), &
                  & source=(0.0_wp, 0.0_wp))
            ELSE
               stream%batch_size = MIN(selected_batch_size, selected_image_count)
               stream%batch_first = selected_image_first
               stream%batch_count = MIN(stream%batch_size, selected_image_count)
               stream%batch_last = selected_image_last
               stream%image_first = selected_image_first
               IF (PRESENT(local_image_first)) THEN
                  stream%partial_image_first = local_image_first
                  IF (local_image_count > 0) THEN
                     stream%partial_image_last = &
                        & local_image_first + (local_image_count - 1)
                  ELSE
                     stream%partial_image_last = 0
                  END IF
                  stream%partial_image_distribution_set = .TRUE.
                  stream%partial_compact_image_storage = .TRUE.
               END IF
               ALLOCATE (stream%inverse_image(nk), source=0)
               ALLOCATE (stream%inverse_coeff(nk), &
                  & source=(0.0_wp, 0.0_wp))
               ALLOCATE (stream%density_signature(nk), source=0_int64)
               ALLOCATE (stream%overlap_signature(nk), source=0_int64)
               ALLOCATE (stream%block_signature_set(nk), source=.FALSE.)
               ALLOCATE (rep_to_cell(nk), source=0)
               DO icell = 1, nk
                  grid = MODULO(stream%reps(:, icell), nmesh)
                  IF (ANY(grid /= stream%reps(:, icell))) THEN
                     CALL clear_exchange_stream(stream)
                     CALL fatal_error(error, &
                        & "g-xTB image batch has a noncanonical representative")
                     RETURN
                  END IF
                  igrid = 1 + grid(1) &
                     & + nmesh(1)*(grid(2) + nmesh(2)*grid(3))
                  IF (rep_to_cell(igrid) /= 0) THEN
                     CALL clear_exchange_stream(stream)
                     CALL fatal_error(error, &
                        & "g-xTB image batch has duplicate representatives")
                     RETURN
                  END IF
                  rep_to_cell(igrid) = icell
               END DO
               IF (ANY(rep_to_cell == 0)) THEN
                  CALL clear_exchange_stream(stream)
                  CALL fatal_error(error, &
                     & "g-xTB image batch has missing representatives")
                  RETURN
               END IF
               DO icell = 1, nk
                  inverse_grid = MODULO(-stream%reps(:, icell), nmesh)
                  igrid = 1 + inverse_grid(1) &
                     & + nmesh(1)*(inverse_grid(2) + nmesh(2)*inverse_grid(3))
                  inverse_cell = rep_to_cell(igrid)
                  stream%inverse_image(icell) = inverse_cell
                  ! The grid character of R+(-R mod n) is unity for every k.
                  ! Only the supplied weight sum and cached twist character
                  ! remain in this exact inverse-image contraction coefficient.
                  stream%inverse_coeff(icell) = SUM(weights)*CONJG( &
                     & stream%twist_phase(icell) &
                     & *stream%twist_phase(inverse_cell))
               END DO
               ALLOCATE (stream%amat_r(exchange%nao, exchange%nao, &
                  & selected_storage_capacity, wfn%nspin), &
                  & stream%amat_partner_r(exchange%nao, exchange%nao, &
                  & selected_storage_capacity, wfn%nspin), &
                  & stream%cmat_r(exchange%nao, exchange%nao, &
                  & selected_storage_capacity, wfn%nspin), &
                  & stream%vmat_r(exchange%nao, exchange%nao, &
                  & selected_storage_capacity, wfn%nspin), &
                  & source=(0.0_wp, 0.0_wp))
               ALLOCATE (stream%gdiagP(exchange%nao, wfn%nspin), &
                  & stream%gdiagSP(exchange%nao, wfn%nspin), &
                  & stream%gdiagSPS(exchange%nao, wfn%nspin), &
                  & source=(0.0_wp, 0.0_wp))
               IF (selected_transform_backend /= bvk_transform_dense) THEN
                  npair = exchange%nao*exchange%nao
                  ! Keep the nontrivial compact path strictly below a complete
                  ! AO-pair mesh.  A one-function basis is already a scalar and
                  ! cannot be subdivided further.
                  IF (npair > 1) THEN
                     selected_pair_batch_size = MIN(selected_pair_batch_size, &
                        & npair - 1)
                  ELSE
                     selected_pair_batch_size = 1
                  END IF
                  stream%transform_pair_batch_size = selected_pair_batch_size
                  stream%transform_pair_first = 1
                  stream%transform_pair_count = MIN( &
                     & selected_pair_batch_size, npair)
                  stream%transform_field_count = MERGE(6, 3, &
                     & selected_reverse_only)
                  ALLOCATE (stream%transform_k(selected_pair_batch_size, nk, &
                     & wfn%nspin, stream%transform_field_count), &
                     & source=(0.0_wp, 0.0_wp))
               END IF
            END IF
            ALLOCATE (stream%vsh(exchange%nsh), source=0.0_wp)
            stream%energy = 0.0_wp
            stream%peak_complex_elements = &
               & exchange_stream_complex_elements(stream)
            IF (selected_mode == cp2k_exchange_stream_batched) THEN
               stream%peak_complex_elements = MAX( &
                  & stream%peak_complex_elements &
                  & + 2_int64*INT(exchange%nao, int64)**2, &
                  & stream%peak_complex_elements &
                  & + 5_int64*INT(exchange%nao, int64)**2 &
                  & + 6_int64*INT(exchange%nao, int64), &
                  & stream%peak_complex_elements &
                  & + INT(exchange%nao, int64)**2*INT(wfn%nspin, int64) &
                  & + 4_int64*INT(exchange%nao, int64)**2 &
                  & + INT(SIZE(stream%amat_r, 3), int64))
               IF (selected_transform_backend /= bvk_transform_dense) THEN
                  stream%peak_complex_elements = MAX( &
                     & stream%peak_complex_elements, &
                     & exchange_stream_complex_elements(stream) &
                     & + 3_int64*INT(stream%transform_pair_batch_size, int64) &
                     & *INT(nk, int64))
               END IF
            END IF
            stream%active = .TRUE.
            IF (selected_reverse_only) THEN
               stream%applied = .TRUE.
               stream%forward_pulled = .TRUE.
               stream%batch_result_pulled = .TRUE.
               CALL initialize_exchange_stream_reverse(stream, mol, error)
               IF (ALLOCATED(error)) THEN
                  CALL clear_exchange_stream(stream)
                  RETURN
               END IF
            END IF
         class default
            CALL fatal_error(error, &
               & "exchange stream requires the exchange_fock implementation")
         END SELECT
      class default
         CALL fatal_error(error, "unexpected g-xTB exchange cache type")
      END SELECT
   END SUBROUTINE cp2k_exchange_stream_begin

!> Start one MPI-free partial k-to-R accumulator transaction.
!>
!> The object is an ordinary bounded exchange stream and can therefore enter
!> the unchanged batch apply/pull protocol after a successful import.  Partial
!> contributors may be empty; global ownership is proved only at import.
!> `transaction_id` is a mandatory host-generated 128-bit identity.  CP2K
!> must create a fresh nonzero value for every logical exchange transaction
!> and broadcast the same value to all contributors; reuse would void stale-
!> payload protection and is therefore forbidden by the public host contract.
   SUBROUTINE cp2k_exchange_partial_begin(stream, calc, mol, ecache, plan, wfn, nmesh, &
      & kfrac, weights, batch_size, transaction_id, error, image_first, &
      & image_count, reverse_only, local_image_first, local_image_count, &
      & transform_backend, transform_pair_batch_size)
      TYPE(cp2k_exchange_stream_type), INTENT(inout)     :: stream
      TYPE(xtb_calculator), INTENT(inout)                :: calc
      TYPE(structure_type), INTENT(in)                   :: mol
      TYPE(container_cache), INTENT(inout)               :: ecache
      TYPE(cp2k_stream_plan_type), INTENT(inout) :: plan
      TYPE(wavefunction_type), INTENT(in)                :: wfn
      INTEGER, INTENT(in)                                :: nmesh(3)
      REAL(wp), INTENT(in)                               :: kfrac(:, :), weights(:)
      INTEGER, INTENT(in)                                :: batch_size
      INTEGER(int64), INTENT(in)                         :: transaction_id(2)
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error
      INTEGER, INTENT(in), OPTIONAL                      :: image_first, image_count
      LOGICAL, INTENT(in), OPTIONAL                      :: reverse_only
      INTEGER, INTENT(in), OPTIONAL                      :: local_image_first, local_image_count, &
                                                            transform_backend, &
                                                            transform_pair_batch_size

      IF (ALL(transaction_id == 0_int64)) THEN
         CALL fatal_error(error, &
            & "partial exchange requires a nonzero host transaction identity")
         RETURN
      END IF

      CALL cp2k_exchange_stream_begin(stream, calc, mol, ecache, plan, wfn, nmesh, &
         & kfrac, weights, error, mode=cp2k_exchange_stream_batched, &
         & batch_size=batch_size, image_first=image_first, &
         & image_count=image_count, reverse_only=reverse_only, &
         & local_image_first=local_image_first, &
         & local_image_count=local_image_count, &
         & transform_backend=transform_backend, &
         & transform_pair_batch_size=transform_pair_batch_size)
      IF (ALLOCATED(error)) RETURN
      stream%partial_active = .TRUE.
      stream%partial_generation = 1
      stream%partial_transaction_id = transaction_id
      ALLOCATE (stream%partial_forward_owned(stream%nk), source=.FALSE.)
      CALL initialize_partial_state_signature(stream, mol)
   END SUBROUTINE cp2k_exchange_partial_begin

!> Add one uniquely owned Bloch block to a forward partial accumulator.
   SUBROUTINE cp2k_exchange_partial_push(stream, ik, density, overlap, error)
      TYPE(cp2k_exchange_stream_type), INTENT(inout)     :: stream
      INTEGER, INTENT(in)                                :: ik
      COMPLEX(wp), INTENT(in)                            :: density(:, :, :), overlap(:, :)
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      IF (.NOT. stream%active .OR. .NOT. stream%partial_active) THEN
         CALL fatal_error(error, "g-xTB exchange partial stream is not active")
         RETURN
      END IF
      IF (stream%reverse_only) THEN
         CALL fatal_error(error, &
            & "forward partial push is invalid in a reverse-only transaction")
         RETURN
      END IF
      IF (stream%partial_compact_image_storage) THEN
         CALL fatal_error(error, &
            & "range-local g-xTB importer cannot accept forward contributors")
         RETURN
      END IF
      IF (stream%partial_forward_imported) THEN
         CALL fatal_error(error, "g-xTB exchange partial push followed import")
         RETURN
      END IF
      IF (stream%partial_forward_exported) THEN
         CALL fatal_error(error, "g-xTB exchange partial push followed export")
         RETURN
      END IF
      ! Fortran does not guarantee short-circuit evaluation below.
      IF (ik < 1 .OR. ik > stream%nk) THEN
         CALL fatal_error(error, "g-xTB exchange partial k-point index out of range")
         RETURN
      END IF
      IF (stream%transform_backend /= bvk_transform_dense &
         & .AND. stream%transform_pair_first > 1 &
         & .AND. .NOT. stream%partial_forward_owned(ik)) THEN
         CALL fatal_error(error, &
            & "compact partial forward k-point ownership changed between tiles")
         RETURN
      END IF
      CALL cp2k_exchange_stream_push(stream, ik, density, overlap, error)
      IF (.NOT. ALLOCATED(error)) stream%partial_forward_owned(ik) = .TRUE.
   END SUBROUTINE cp2k_exchange_partial_push

!> Export one additive forward payload and its local ownership vector.
   SUBROUTINE cp2k_exchange_partial_export_forward(stream, descriptor, &
      & owner_count, buffer, error)
      TYPE(cp2k_exchange_stream_type), INTENT(inout)     :: stream
      TYPE(cp2k_exchange_partial_descriptor_type), &
         INTENT(out)                                     :: descriptor
      INTEGER, ALLOCATABLE, INTENT(out)                  :: owner_count(:)
      TYPE(cp2k_exchange_partial_forward_type), &
         INTENT(out)                                     :: buffer
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      INTEGER                                            :: nblock_default, total_default
      INTEGER(int64)                                     :: first, last, nblock, total
      LOGICAL                                            :: extents_ok

      IF (.NOT. stream%active .OR. .NOT. stream%partial_active) THEN
         CALL fatal_error(error, "g-xTB exchange partial stream is not active")
         RETURN
      END IF
      IF (stream%reverse_only) THEN
         CALL fatal_error(error, &
            & "forward partial export is invalid in a reverse-only transaction")
         RETURN
      END IF
      IF (stream%partial_forward_imported) THEN
         CALL fatal_error(error, "g-xTB exchange forward payload was imported")
         RETURN
      END IF
      IF (stream%partial_forward_exported) THEN
         CALL fatal_error(error, "g-xTB exchange forward payload was exported")
         RETURN
      END IF
      IF (stream%applied .OR. stream%batch_applied) THEN
         CALL fatal_error(error, &
            & "g-xTB exchange forward payload export followed kernel apply")
         RETURN
      END IF
      IF (stream%partial_compact_image_storage) THEN
         CALL fatal_error(error, &
            & "range-local g-xTB importer cannot export a forward payload")
         RETURN
      END IF

      CALL get_partial_descriptor(stream, 1, descriptor)
      ALLOCATE (owner_count(stream%nk), source=0)
      WHERE (stream%seen)
         owner_count = 1
      END WHERE
      CALL partial_forward_extents(stream%nao, stream%batch_count, stream%nspin, &
         & nblock, total, extents_ok)
      IF (.NOT. extents_ok) THEN
         CALL fatal_error(error, &
            & "g-xTB exchange forward payload exceeds indexing limit")
         RETURN
      END IF
      nblock_default = INT(nblock)
      total_default = INT(total)
      ALLOCATE (buffer%values(total_default))
      first = 1_int64
      last = nblock
      buffer%values(first:last) = RESHAPE( &
         & stream%amat_r(:, :, 1:stream%batch_count, :), [nblock_default])
      first = last + 1_int64
      last = last + nblock
      buffer%values(first:last) = RESHAPE( &
         & stream%amat_partner_r(:, :, 1:stream%batch_count, :), [nblock_default])
      first = last + 1_int64
      last = last + nblock
      buffer%values(first:last) = RESHAPE( &
         & stream%cmat_r(:, :, 1:stream%batch_count, :), [nblock_default])
      first = last + 1_int64
      last = last + nblock
      buffer%values(first:last) = RESHAPE( &
         & stream%vmat_r(:, :, 1:stream%batch_count, :), [nblock_default])
      IF (.NOT. ALL(ieee_is_finite(REAL(buffer%values, wp))) &
         & .OR. .NOT. ALL(ieee_is_finite(AIMAG(buffer%values)))) THEN
         CALL fatal_error(error, "non-finite g-xTB exchange forward payload")
         RETURN
      END IF
      buffer%image_first = descriptor%batch_first
      buffer%image_count = descriptor%batch_count
      buffer%abi_version = cp2k_exchange_partial_abi_version
      buffer%descriptor_signature = partial_descriptor_signature(descriptor)
      stream%partial_forward_exported = .TRUE.
   END SUBROUTINE cp2k_exchange_partial_export_forward

!> Number of complex scalars in an exported forward payload.
!>
!> Construction fails closed before allocation if the count exceeds
!> `huge(default integer)`; a host may still chunk MPI reductions below that
!> provider-side indexing limit.
   INTEGER(int64) FUNCTION cp2k_exchange_partial_forward_size(buffer) &
      & RESULT(nelements)
      TYPE(cp2k_exchange_partial_forward_type), &
         INTENT(in)                                      :: buffer

      IF (buffer%abi_version == cp2k_exchange_partial_abi_version &
         & .AND. ALLOCATED(buffer%values)) THEN
         nelements = SIZE(buffer%values, kind=int64)
      ELSE
         nelements = -1_int64
      END IF
   END FUNCTION cp2k_exchange_partial_forward_size

!> Copy an opaque forward payload into contiguous host/MPI storage.
   SUBROUTINE cp2k_exchange_partial_forward_pack(buffer, values, error)
      TYPE(cp2k_exchange_partial_forward_type), &
         INTENT(in)                                      :: buffer
      COMPLEX(wp), ALLOCATABLE, INTENT(out)              :: values(:)
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      IF (buffer%abi_version /= cp2k_exchange_partial_abi_version &
         & .OR. .NOT. ALLOCATED(buffer%values)) THEN
         CALL fatal_error(error, "invalid g-xTB exchange forward payload")
         RETURN
      END IF
      ALLOCATE (values, source=buffer%values)
   END SUBROUTINE cp2k_exchange_partial_forward_pack

!> Copy only one contiguous image subset of a full forward payload.
   SUBROUTINE cp2k_exchange_partial_forward_pack_range(buffer, descriptor, &
      & image_first, image_count, values, error)
      TYPE(cp2k_exchange_partial_forward_type), &
         INTENT(in)                                      :: buffer
      TYPE(cp2k_exchange_partial_descriptor_type), &
         INTENT(in)                                      :: descriptor
      INTEGER, INTENT(in)                                :: image_first, image_count
      COMPLEX(wp), ALLOCATABLE, INTENT(out)              :: values(:)
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      INTEGER                                            :: field, image_offset, spin
      INTEGER(int64)                                     :: chunk, dst_first, dst_last, full_block, &
                                                            full_total, pair, range_block, &
                                                            range_total, src_first, src_last
      LOGICAL                                            :: extents_ok, range_ok

      IF (.NOT. partial_descriptor_is_valid(descriptor) &
         & .OR. descriptor%phase /= 1) THEN
         CALL fatal_error(error, "inconsistent g-xTB exchange forward descriptor")
         RETURN
      END IF
      CALL partial_payload_range(descriptor, image_first, image_count, &
         & image_offset, range_ok)
      CALL partial_forward_extents(descriptor%nao, descriptor%batch_count, &
         & descriptor%nspin, full_block, full_total, extents_ok)
      IF (.NOT. range_ok .OR. .NOT. extents_ok) THEN
         CALL fatal_error(error, "invalid g-xTB exchange forward payload range")
         RETURN
      END IF
      IF (buffer%abi_version /= cp2k_exchange_partial_abi_version &
         & .OR. .NOT. ALLOCATED(buffer%values) &
         & .OR. buffer%image_first /= descriptor%batch_first &
         & .OR. buffer%image_count /= descriptor%batch_count &
         & .OR. SIZE(buffer%values, kind=int64) /= full_total &
         & .OR. ANY(buffer%descriptor_signature /= &
         & partial_descriptor_signature(descriptor))) THEN
         CALL fatal_error(error, "inconsistent g-xTB exchange forward payload")
         RETURN
      END IF
      CALL partial_forward_extents(descriptor%nao, image_count, &
         & descriptor%nspin, range_block, range_total, extents_ok)
      IF (.NOT. extents_ok) THEN
         CALL fatal_error(error, &
            & "g-xTB exchange forward range exceeds indexing limit")
         RETURN
      END IF
      ALLOCATE (values(INT(range_total)))
      pair = INT(descriptor%nao, int64)**2
      chunk = pair*INT(image_count, int64)
      DO field = 0, 3
         DO spin = 0, descriptor%nspin - 1
            src_first = INT(field, int64)*full_block &
               & + INT(spin, int64)*pair*INT(descriptor%batch_count, int64) &
               & + INT(image_offset, int64)*pair + 1_int64
            src_last = src_first + chunk - 1_int64
            dst_first = INT(field, int64)*range_block &
               & + INT(spin, int64)*chunk + 1_int64
            dst_last = dst_first + chunk - 1_int64
            IF (chunk > 0_int64) values(dst_first:dst_last) = &
               & buffer%values(src_first:src_last)
         END DO
      END DO
      IF (.NOT. ALL(ieee_is_finite(REAL(values, wp))) &
         & .OR. .NOT. ALL(ieee_is_finite(AIMAG(values)))) THEN
         CALL fatal_error(error, "non-finite g-xTB exchange forward range")
         RETURN
      END IF
   END SUBROUTINE cp2k_exchange_partial_forward_pack_range

!> Recreate one opaque forward image-range payload after host reduction.
   SUBROUTINE cp2k_exchange_partial_forward_unpack_range(buffer, descriptor, &
      & image_first, image_count, values, error)
      TYPE(cp2k_exchange_partial_forward_type), &
         INTENT(out)                                     :: buffer
      TYPE(cp2k_exchange_partial_descriptor_type), &
         INTENT(in)                                      :: descriptor
      INTEGER, INTENT(in)                                :: image_first, image_count
      COMPLEX(wp), INTENT(in)                            :: values(:)
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      INTEGER                                            :: image_offset
      INTEGER(int64)                                     :: expected, nblock
      LOGICAL                                            :: extents_ok, range_ok

      IF (.NOT. partial_descriptor_is_valid(descriptor) &
         & .OR. descriptor%phase /= 1) THEN
         CALL fatal_error(error, "inconsistent g-xTB exchange forward descriptor")
         RETURN
      END IF
      CALL partial_payload_range(descriptor, image_first, image_count, &
         & image_offset, range_ok)
      CALL partial_forward_extents(descriptor%nao, image_count, descriptor%nspin, &
         & nblock, expected, extents_ok)
      IF (.NOT. range_ok .OR. .NOT. extents_ok &
         & .OR. SIZE(values, kind=int64) /= expected) THEN
         CALL fatal_error(error, "inconsistent g-xTB exchange forward range")
         RETURN
      END IF
      IF (.NOT. ALL(ieee_is_finite(REAL(values, wp))) &
         & .OR. .NOT. ALL(ieee_is_finite(AIMAG(values)))) THEN
         CALL fatal_error(error, "non-finite g-xTB exchange forward range")
         RETURN
      END IF
      ALLOCATE (buffer%values, source=values)
      buffer%image_first = image_first
      buffer%image_count = image_count
      buffer%abi_version = cp2k_exchange_partial_abi_version
      buffer%descriptor_signature = partial_descriptor_signature(descriptor)
   END SUBROUTINE cp2k_exchange_partial_forward_unpack_range

!> Recreate an opaque forward payload from a host-reduced contiguous buffer.
!>
!> `descriptor` must be the exact collectively validated companion descriptor
!> exported with `values`.  Rebinding saved raw values to a fresh descriptor is
!> forbidden by the host contract.
   SUBROUTINE cp2k_exchange_partial_forward_unpack(buffer, descriptor, values, &
      & error)
      TYPE(cp2k_exchange_partial_forward_type), &
         INTENT(out)                                     :: buffer
      TYPE(cp2k_exchange_partial_descriptor_type), &
         INTENT(in)                                      :: descriptor
      COMPLEX(wp), INTENT(in)                            :: values(:)
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      INTEGER(int64)                                     :: expected, nblock
      LOGICAL                                            :: extents_ok

      IF (.NOT. partial_descriptor_is_valid(descriptor) &
         & .OR. descriptor%phase /= 1) THEN
         CALL fatal_error(error, "inconsistent g-xTB exchange forward descriptor")
         RETURN
      END IF
      CALL partial_forward_extents(descriptor%nao, descriptor%batch_count, &
         & descriptor%nspin, nblock, expected, extents_ok)
      IF (.NOT. extents_ok) THEN
         CALL fatal_error(error, &
            & "g-xTB exchange forward payload exceeds indexing limit")
         RETURN
      END IF
      IF (SIZE(values, kind=int64) /= expected) THEN
         CALL fatal_error(error, "inconsistent g-xTB exchange forward payload")
         RETURN
      END IF
      IF (.NOT. ALL(ieee_is_finite(REAL(values, wp))) &
         & .OR. .NOT. ALL(ieee_is_finite(AIMAG(values)))) THEN
         CALL fatal_error(error, "non-finite g-xTB exchange forward payload")
         RETURN
      END IF
      ALLOCATE (buffer%values, source=values)
      buffer%image_first = descriptor%batch_first
      buffer%image_count = descriptor%batch_count
      buffer%abi_version = cp2k_exchange_partial_abi_version
      buffer%descriptor_signature = partial_descriptor_signature(descriptor)
   END SUBROUTINE cp2k_exchange_partial_forward_unpack

!> Import the globally summed forward sources before one nonlinear apply.
   SUBROUTINE cp2k_exchange_partial_import_forward(stream, descriptor, &
      & owner_count, buffer, calc, mol, ecache, plan, error)
      TYPE(cp2k_exchange_stream_type), INTENT(inout)     :: stream
      TYPE(cp2k_exchange_partial_descriptor_type), &
         INTENT(in)                                      :: descriptor
      INTEGER, INTENT(in)                                :: owner_count(:)
      TYPE(cp2k_exchange_partial_forward_type), &
         INTENT(in)                                      :: buffer
      TYPE(xtb_calculator), INTENT(in)                   :: calc
      TYPE(structure_type), INTENT(in)                   :: mol
      TYPE(container_cache), INTENT(in)                  :: ecache
      TYPE(cp2k_stream_plan_type), INTENT(in) :: plan
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      INTEGER(int64)                                     :: first, last, nblock, total
      LOGICAL                                            :: extents_ok
      TYPE(cp2k_exchange_partial_descriptor_type)        :: current

      IF (.NOT. stream%active .OR. .NOT. stream%partial_active) THEN
         CALL fatal_error(error, "g-xTB exchange partial stream is not active")
         RETURN
      END IF
      IF (stream%reverse_only) THEN
         CALL fatal_error(error, &
            & "forward partial import is invalid in a reverse-only transaction")
         RETURN
      END IF
      IF (stream%partial_forward_imported) THEN
         CALL fatal_error(error, "g-xTB exchange forward payload was imported")
         RETURN
      END IF
      IF (stream%applied .OR. stream%batch_applied) THEN
         CALL fatal_error(error, &
            & "g-xTB exchange forward payload import followed kernel apply")
         RETURN
      END IF
      IF (stream%partial_compact_image_storage) THEN
         CALL fatal_error(error, &
            & "range-local g-xTB importer requires a forward range payload")
         RETURN
      END IF
      IF (.NOT. exchange_stream_state_matches(stream, calc, mol, ecache, plan)) THEN
         CALL fatal_error(error, &
            & "g-xTB exchange partial state was invalidated before import")
         RETURN
      END IF
      CALL get_partial_descriptor(stream, 1, current)
      IF (.NOT. cp2k_exchange_partial_descriptors_match(current, descriptor)) THEN
         CALL fatal_error(error, "g-xTB exchange partial descriptor mismatch")
         RETURN
      END IF
      IF (SIZE(owner_count) /= stream%nk) THEN
         CALL fatal_error(error, "inconsistent g-xTB exchange owner-count size")
         RETURN
      END IF
      IF (ANY(owner_count < 0)) THEN
         CALL fatal_error(error, "negative g-xTB exchange owner count")
         RETURN
      END IF
      IF (ANY(owner_count > 1)) THEN
         CALL fatal_error(error, "duplicate g-xTB exchange k-point ownership")
         RETURN
      END IF
      IF (ANY(owner_count == 0)) THEN
         CALL fatal_error(error, "missing g-xTB exchange k-point ownership")
         RETURN
      END IF
      CALL partial_forward_extents(stream%nao, stream%batch_count, stream%nspin, &
         & nblock, total, extents_ok)
      IF (.NOT. extents_ok) THEN
         CALL fatal_error(error, &
            & "g-xTB exchange forward payload exceeds indexing limit")
         RETURN
      END IF
      IF (buffer%abi_version /= cp2k_exchange_partial_abi_version &
         & .OR. .NOT. ALLOCATED(buffer%values)) THEN
         CALL fatal_error(error, "inconsistent g-xTB exchange forward payload")
         RETURN
      END IF
      IF (SIZE(buffer%values, kind=int64) /= total &
         & .OR. buffer%image_first /= descriptor%batch_first &
         & .OR. buffer%image_count /= descriptor%batch_count &
         & .OR. ANY(buffer%descriptor_signature /= &
         & partial_descriptor_signature(descriptor))) THEN
         CALL fatal_error(error, "inconsistent g-xTB exchange forward payload")
         RETURN
      END IF
      IF (.NOT. ALL(ieee_is_finite(REAL(buffer%values, wp))) &
         & .OR. .NOT. ALL(ieee_is_finite(AIMAG(buffer%values)))) THEN
         CALL fatal_error(error, "non-finite g-xTB exchange forward payload")
         RETURN
      END IF

      first = 1_int64
      last = nblock
      stream%amat_r(:, :, 1:stream%batch_count, :) = &
         & RESHAPE(buffer%values(first:last), &
         & SHAPE(stream%amat_r(:, :, 1:stream%batch_count, :)))
      first = last + 1_int64
      last = last + nblock
      stream%amat_partner_r(:, :, 1:stream%batch_count, :) = &
         & RESHAPE(buffer%values(first:last), &
         & SHAPE(stream%amat_partner_r(:, :, 1:stream%batch_count, :)))
      first = last + 1_int64
      last = last + nblock
      stream%cmat_r(:, :, 1:stream%batch_count, :) = &
         & RESHAPE(buffer%values(first:last), &
         & SHAPE(stream%cmat_r(:, :, 1:stream%batch_count, :)))
      first = last + 1_int64
      last = last + nblock
      stream%vmat_r(:, :, 1:stream%batch_count, :) = &
         & RESHAPE(buffer%values(first:last), &
         & SHAPE(stream%vmat_r(:, :, 1:stream%batch_count, :)))
      stream%seen = .TRUE.
      IF (stream%transform_backend /= bvk_transform_dense) THEN
         stream%transform_pair_first = 0
         stream%transform_pair_count = 0
         stream%transform_k = (0.0_wp, 0.0_wp)
      END IF
      stream%partial_forward_imported = .TRUE.
   END SUBROUTINE cp2k_exchange_partial_import_forward

!> Import only this rank's compact intersection with the current image batch.
   SUBROUTINE cp2k_exchange_partial_import_forward_range(stream, descriptor, &
      & owner_count, buffer, calc, mol, ecache, plan, error)
      TYPE(cp2k_exchange_stream_type), INTENT(inout)     :: stream
      TYPE(cp2k_exchange_partial_descriptor_type), &
         INTENT(in)                                      :: descriptor
      INTEGER, INTENT(in)                                :: owner_count(:)
      TYPE(cp2k_exchange_partial_forward_type), &
         INTENT(in)                                      :: buffer
      TYPE(xtb_calculator), INTENT(in)                   :: calc
      TYPE(structure_type), INTENT(in)                   :: mol
      TYPE(container_cache), INTENT(in)                  :: ecache
      TYPE(cp2k_stream_plan_type), INTENT(in) :: plan
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      INTEGER                                            :: first_global, first_local, last_local, &
                                                            local_count
      INTEGER(int64)                                     :: first, last, nblock, total
      LOGICAL                                            :: extents_ok
      TYPE(cp2k_exchange_partial_descriptor_type)        :: current

      IF (.NOT. stream%active .OR. .NOT. stream%partial_active &
         & .OR. .NOT. stream%partial_compact_image_storage) THEN
         CALL fatal_error(error, &
            & "g-xTB exchange forward range importer is not active")
         RETURN
      END IF
      IF (stream%reverse_only) THEN
         CALL fatal_error(error, &
            & "forward range import is invalid in a reverse-only transaction")
         RETURN
      END IF
      IF (stream%partial_forward_imported) THEN
         CALL fatal_error(error, "g-xTB exchange forward payload was imported")
         RETURN
      END IF
      IF (stream%applied .OR. stream%batch_applied) THEN
         CALL fatal_error(error, &
            & "g-xTB exchange forward range import followed kernel apply")
         RETURN
      END IF
      IF (.NOT. exchange_stream_state_matches(stream, calc, mol, ecache, plan)) THEN
         CALL fatal_error(error, &
            & "g-xTB exchange partial state was invalidated before range import")
         RETURN
      END IF
      CALL get_partial_descriptor(stream, 1, current)
      IF (.NOT. cp2k_exchange_partial_descriptors_match(current, descriptor)) THEN
         CALL fatal_error(error, "g-xTB exchange partial descriptor mismatch")
         RETURN
      END IF
      IF (SIZE(owner_count) /= stream%nk .OR. ANY(owner_count /= 1)) THEN
         CALL fatal_error(error, &
            & "incomplete g-xTB exchange k-point ownership for range import")
         RETURN
      END IF
      CALL partial_local_batch_bounds(stream, first_local, last_local, &
         & local_count, first_global)
      CALL partial_forward_extents(stream%nao, local_count, stream%nspin, &
         & nblock, total, extents_ok)
      IF (.NOT. extents_ok .OR. local_count > SIZE(stream%amat_r, 3)) THEN
         CALL fatal_error(error, &
            & "g-xTB exchange forward range exceeds local storage")
         RETURN
      END IF
      IF (buffer%abi_version /= cp2k_exchange_partial_abi_version &
         & .OR. .NOT. ALLOCATED(buffer%values) &
         & .OR. buffer%image_first /= MERGE(first_global, 0, local_count > 0) &
         & .OR. buffer%image_count /= local_count &
         & .OR. SIZE(buffer%values, kind=int64) /= total &
         & .OR. ANY(buffer%descriptor_signature /= &
         & partial_descriptor_signature(descriptor))) THEN
         CALL fatal_error(error, "inconsistent g-xTB exchange forward range")
         RETURN
      END IF
      IF (.NOT. ALL(ieee_is_finite(REAL(buffer%values, wp))) &
         & .OR. .NOT. ALL(ieee_is_finite(AIMAG(buffer%values)))) THEN
         CALL fatal_error(error, "non-finite g-xTB exchange forward range")
         RETURN
      END IF

      IF (local_count > 0) THEN
         first = 1_int64
         last = nblock
         stream%amat_r(:, :, 1:local_count, :) = &
            & RESHAPE(buffer%values(first:last), &
            & SHAPE(stream%amat_r(:, :, 1:local_count, :)))
         first = last + 1_int64
         last = last + nblock
         stream%amat_partner_r(:, :, 1:local_count, :) = &
            & RESHAPE(buffer%values(first:last), &
            & SHAPE(stream%amat_partner_r(:, :, 1:local_count, :)))
         first = last + 1_int64
         last = last + nblock
         stream%cmat_r(:, :, 1:local_count, :) = &
            & RESHAPE(buffer%values(first:last), &
            & SHAPE(stream%cmat_r(:, :, 1:local_count, :)))
         first = last + 1_int64
         last = last + nblock
         stream%vmat_r(:, :, 1:local_count, :) = &
            & RESHAPE(buffer%values(first:last), &
            & SHAPE(stream%vmat_r(:, :, 1:local_count, :)))
      END IF
      stream%seen = .TRUE.
      IF (stream%transform_backend /= bvk_transform_dense) THEN
         stream%transform_pair_first = 0
         stream%transform_pair_count = 0
         stream%transform_k = (0.0_wp, 0.0_wp)
      END IF
      stream%partial_forward_imported = .TRUE.
   END SUBROUTINE cp2k_exchange_partial_import_forward_range

!> Add one uniquely owned Bloch block to a reverse partial accumulator.
   SUBROUTINE cp2k_exchange_partial_reverse_push(stream, ik, density, overlap, &
      & error)
      TYPE(cp2k_exchange_stream_type), INTENT(inout)     :: stream
      INTEGER, INTENT(in)                                :: ik
      COMPLEX(wp), INTENT(in)                            :: density(:, :, :), overlap(:, :)
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      IF (.NOT. stream%active .OR. .NOT. stream%partial_active &
         & .OR. .NOT. stream%reverse_transaction) THEN
         CALL fatal_error(error, &
            & "g-xTB exchange reverse partial stream is not active")
         RETURN
      END IF
      IF (stream%partial_reverse_imported) THEN
         CALL fatal_error(error, &
            & "g-xTB exchange reverse partial push followed import")
         RETURN
      END IF
      IF (stream%partial_compact_image_storage) THEN
         CALL fatal_error(error, &
            & "range-local g-xTB importer cannot accept reverse contributors")
         RETURN
      END IF
      IF (stream%partial_reverse_exported) THEN
         CALL fatal_error(error, &
            & "g-xTB exchange reverse partial push followed export")
         RETURN
      END IF
      ! Validate before indexing the ownership vector, also for dense transforms.
      IF (ik < 1 .OR. ik > stream%nk) THEN
         CALL fatal_error(error, "g-xTB exchange reverse partial k-point index out of range")
         RETURN
      END IF
      IF (stream%transform_backend /= bvk_transform_dense &
         & .AND. stream%transform_pair_first > 1 &
         & .AND. .NOT. stream%partial_forward_owned(ik)) THEN
         CALL fatal_error(error, &
            & "compact partial reverse k-point ownership changed between tiles")
         RETURN
      END IF
      CALL cp2k_exchange_stream_reverse_push(stream, ik, density, overlap, error)
      IF (.NOT. ALLOCATED(error)) stream%partial_forward_owned(ik) = .TRUE.
   END SUBROUTINE cp2k_exchange_partial_reverse_push

!> Export one additive reverse-source payload and local ownership vector.
   SUBROUTINE cp2k_exchange_partial_export_reverse(stream, descriptor, &
      & owner_count, buffer, error)
      TYPE(cp2k_exchange_stream_type), INTENT(inout)     :: stream
      TYPE(cp2k_exchange_partial_descriptor_type), &
         INTENT(out)                                     :: descriptor
      INTEGER, ALLOCATABLE, INTENT(out)                  :: owner_count(:)
      TYPE(cp2k_exchange_partial_reverse_type), &
         INTENT(out)                                     :: buffer
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      INTEGER                                            :: nblock_default, nonsite_default, &
                                                            total_default
      INTEGER(int64)                                     :: first, last, nblock, nonsite, total
      LOGICAL                                            :: extents_ok

      IF (.NOT. stream%active .OR. .NOT. stream%partial_active &
         & .OR. .NOT. stream%reverse_transaction) THEN
         CALL fatal_error(error, &
            & "g-xTB exchange reverse partial stream is not active")
         RETURN
      END IF
      IF (stream%partial_reverse_imported) THEN
         CALL fatal_error(error, "g-xTB exchange reverse payload was imported")
         RETURN
      END IF
      IF (stream%partial_reverse_exported) THEN
         CALL fatal_error(error, "g-xTB exchange reverse payload was exported")
         RETURN
      END IF
      IF (stream%reverse_batch_applied .OR. stream%reverse_applied) THEN
         CALL fatal_error(error, &
            & "g-xTB exchange reverse payload export followed kernel apply")
         RETURN
      END IF
      IF (stream%partial_compact_image_storage) THEN
         CALL fatal_error(error, &
            & "range-local g-xTB importer cannot export a reverse payload")
         RETURN
      END IF

      CALL get_partial_descriptor(stream, 2, descriptor)
      ALLOCATE (owner_count(stream%nk), source=0)
      WHERE (stream%seen)
         owner_count = 1
      END WHERE
      CALL partial_reverse_extents(stream%nao, stream%batch_count, stream%nspin, &
         & nblock, nonsite, total, extents_ok)
      IF (.NOT. extents_ok) THEN
         CALL fatal_error(error, &
            & "g-xTB exchange reverse payload exceeds indexing limit")
         RETURN
      END IF
      nblock_default = INT(nblock)
      nonsite_default = INT(nonsite)
      total_default = INT(total)
      ALLOCATE (buffer%values(total_default))
      first = 1_int64
      last = nblock
      buffer%values(first:last) = RESHAPE( &
         & stream%amat_r(:, :, 1:stream%batch_count, :), [nblock_default])
      first = last + 1_int64
      last = last + nblock
      buffer%values(first:last) = RESHAPE( &
         & stream%cmat_r(:, :, 1:stream%batch_count, :), [nblock_default])
      first = last + 1_int64
      last = last + nblock
      buffer%values(first:last) = RESHAPE( &
         & stream%vmat_r(:, :, 1:stream%batch_count, :), [nblock_default])
      first = last + 1_int64
      last = last + nblock
      buffer%values(first:last) = RESHAPE( &
         & stream%amat_partner_r(:, :, 1:stream%batch_count, :), [nblock_default])
      first = last + 1_int64
      last = last + nblock
      buffer%values(first:last) = RESHAPE( &
         & stream%reverse_amat_r(:, :, 1:stream%batch_count, :), [nblock_default])
      first = last + 1_int64
      last = last + nblock
      buffer%values(first:last) = RESHAPE( &
         & stream%reverse_vmat_r(:, :, 1:stream%batch_count, :), [nblock_default])
      first = last + 1_int64
      last = last + nonsite
      buffer%values(first:last) = RESHAPE(stream%gdiagP, [nonsite_default])
      first = last + 1_int64
      last = last + nonsite
      buffer%values(first:last) = RESHAPE(stream%gdiagSP, [nonsite_default])
      first = last + 1_int64
      last = last + nonsite
      buffer%values(first:last) = RESHAPE(stream%gdiagSPS, [nonsite_default])
      first = last + 1_int64
      last = last + nonsite
      buffer%values(first:last) = RESHAPE(stream%reverse_bdiagSP, [nonsite_default])
      first = last + 1_int64
      last = last + nonsite
      buffer%values(first:last) = RESHAPE(stream%reverse_bdiagSPS, &
         & [nonsite_default])
      IF (.NOT. ALL(ieee_is_finite(REAL(buffer%values, wp))) &
         & .OR. .NOT. ALL(ieee_is_finite(AIMAG(buffer%values)))) THEN
         CALL fatal_error(error, "non-finite g-xTB exchange reverse payload")
         RETURN
      END IF
      buffer%image_first = descriptor%batch_first
      buffer%image_count = descriptor%batch_count
      buffer%abi_version = cp2k_exchange_partial_abi_version
      buffer%descriptor_signature = partial_descriptor_signature(descriptor)
      stream%partial_reverse_exported = .TRUE.
   END SUBROUTINE cp2k_exchange_partial_export_reverse

!> Number of complex scalars in an exported reverse payload.
!>
!> Construction fails closed before allocation if the count exceeds
!> `huge(default integer)`; a host may still chunk MPI reductions below that
!> provider-side indexing limit.
   INTEGER(int64) FUNCTION cp2k_exchange_partial_reverse_size(buffer) &
      & RESULT(nelements)
      TYPE(cp2k_exchange_partial_reverse_type), &
         INTENT(in)                                      :: buffer

      IF (buffer%abi_version == cp2k_exchange_partial_abi_version &
         & .AND. ALLOCATED(buffer%values)) THEN
         nelements = SIZE(buffer%values, kind=int64)
      ELSE
         nelements = -1_int64
      END IF
   END FUNCTION cp2k_exchange_partial_reverse_size

!> Copy an opaque reverse payload into contiguous host/MPI storage.
   SUBROUTINE cp2k_exchange_partial_reverse_pack(buffer, values, error)
      TYPE(cp2k_exchange_partial_reverse_type), &
         INTENT(in)                                      :: buffer
      COMPLEX(wp), ALLOCATABLE, INTENT(out)              :: values(:)
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      IF (buffer%abi_version /= cp2k_exchange_partial_abi_version &
         & .OR. .NOT. ALLOCATED(buffer%values)) THEN
         CALL fatal_error(error, "invalid g-xTB exchange reverse payload")
         RETURN
      END IF
      ALLOCATE (values, source=buffer%values)
   END SUBROUTINE cp2k_exchange_partial_reverse_pack

!> Copy one image subset plus the small replicated onsite reverse source.
   SUBROUTINE cp2k_exchange_partial_reverse_pack_range(buffer, descriptor, &
      & image_first, image_count, values, error)
      TYPE(cp2k_exchange_partial_reverse_type), &
         INTENT(in)                                      :: buffer
      TYPE(cp2k_exchange_partial_descriptor_type), &
         INTENT(in)                                      :: descriptor
      INTEGER, INTENT(in)                                :: image_first, image_count
      COMPLEX(wp), ALLOCATABLE, INTENT(out)              :: values(:)
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      INTEGER                                            :: field, image_offset, spin
      INTEGER(int64) :: chunk, dst_first, dst_last, full_block, full_nonsite, full_total, &
         onsite_chunk, pair, range_block, range_nonsite, range_total, src_first, src_last
      LOGICAL                                            :: extents_ok, range_ok

      IF (.NOT. partial_descriptor_is_valid(descriptor) &
         & .OR. descriptor%phase /= 2) THEN
         CALL fatal_error(error, "inconsistent g-xTB exchange reverse descriptor")
         RETURN
      END IF
      CALL partial_payload_range(descriptor, image_first, image_count, &
         & image_offset, range_ok)
      CALL partial_reverse_extents(descriptor%nao, descriptor%batch_count, &
         & descriptor%nspin, full_block, full_nonsite, full_total, extents_ok)
      IF (.NOT. range_ok .OR. .NOT. extents_ok) THEN
         CALL fatal_error(error, "invalid g-xTB exchange reverse payload range")
         RETURN
      END IF
      IF (buffer%abi_version /= cp2k_exchange_partial_abi_version &
         & .OR. .NOT. ALLOCATED(buffer%values) &
         & .OR. buffer%image_first /= descriptor%batch_first &
         & .OR. buffer%image_count /= descriptor%batch_count &
         & .OR. SIZE(buffer%values, kind=int64) /= full_total &
         & .OR. ANY(buffer%descriptor_signature /= &
         & partial_descriptor_signature(descriptor))) THEN
         CALL fatal_error(error, "inconsistent g-xTB exchange reverse payload")
         RETURN
      END IF
      CALL partial_reverse_extents(descriptor%nao, image_count, descriptor%nspin, &
         & range_block, range_nonsite, range_total, extents_ok)
      IF (.NOT. extents_ok .OR. range_nonsite /= full_nonsite) THEN
         CALL fatal_error(error, &
            & "g-xTB exchange reverse range exceeds indexing limit")
         RETURN
      END IF
      ALLOCATE (values(INT(range_total)))
      pair = INT(descriptor%nao, int64)**2
      chunk = pair*INT(image_count, int64)
      DO field = 0, 5
         DO spin = 0, descriptor%nspin - 1
            src_first = INT(field, int64)*full_block &
               & + INT(spin, int64)*pair*INT(descriptor%batch_count, int64) &
               & + INT(image_offset, int64)*pair + 1_int64
            src_last = src_first + chunk - 1_int64
            dst_first = INT(field, int64)*range_block &
               & + INT(spin, int64)*chunk + 1_int64
            dst_last = dst_first + chunk - 1_int64
            IF (chunk > 0_int64) values(dst_first:dst_last) = &
               & buffer%values(src_first:src_last)
         END DO
      END DO
      onsite_chunk = 5_int64*full_nonsite
      IF (onsite_chunk > 0_int64) THEN
         src_first = 6_int64*full_block + 1_int64
         src_last = src_first + onsite_chunk - 1_int64
         dst_first = 6_int64*range_block + 1_int64
         dst_last = dst_first + onsite_chunk - 1_int64
         values(dst_first:dst_last) = buffer%values(src_first:src_last)
      END IF
      IF (.NOT. ALL(ieee_is_finite(REAL(values, wp))) &
         & .OR. .NOT. ALL(ieee_is_finite(AIMAG(values)))) THEN
         CALL fatal_error(error, "non-finite g-xTB exchange reverse range")
         RETURN
      END IF
   END SUBROUTINE cp2k_exchange_partial_reverse_pack_range

!> Recreate one opaque reverse image-range payload after host reduction.
   SUBROUTINE cp2k_exchange_partial_reverse_unpack_range(buffer, descriptor, &
      & image_first, image_count, values, error)
      TYPE(cp2k_exchange_partial_reverse_type), &
         INTENT(out)                                     :: buffer
      TYPE(cp2k_exchange_partial_descriptor_type), &
         INTENT(in)                                      :: descriptor
      INTEGER, INTENT(in)                                :: image_first, image_count
      COMPLEX(wp), INTENT(in)                            :: values(:)
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      INTEGER                                            :: image_offset
      INTEGER(int64)                                     :: expected, nblock, nonsite
      LOGICAL                                            :: extents_ok, range_ok

      IF (.NOT. partial_descriptor_is_valid(descriptor) &
         & .OR. descriptor%phase /= 2) THEN
         CALL fatal_error(error, "inconsistent g-xTB exchange reverse descriptor")
         RETURN
      END IF
      CALL partial_payload_range(descriptor, image_first, image_count, &
         & image_offset, range_ok)
      CALL partial_reverse_extents(descriptor%nao, image_count, descriptor%nspin, &
         & nblock, nonsite, expected, extents_ok)
      IF (.NOT. range_ok .OR. .NOT. extents_ok &
         & .OR. SIZE(values, kind=int64) /= expected) THEN
         CALL fatal_error(error, "inconsistent g-xTB exchange reverse range")
         RETURN
      END IF
      IF (.NOT. ALL(ieee_is_finite(REAL(values, wp))) &
         & .OR. .NOT. ALL(ieee_is_finite(AIMAG(values)))) THEN
         CALL fatal_error(error, "non-finite g-xTB exchange reverse range")
         RETURN
      END IF
      ALLOCATE (buffer%values, source=values)
      buffer%image_first = image_first
      buffer%image_count = image_count
      buffer%abi_version = cp2k_exchange_partial_abi_version
      buffer%descriptor_signature = partial_descriptor_signature(descriptor)
   END SUBROUTINE cp2k_exchange_partial_reverse_unpack_range

!> Recreate an opaque reverse payload from host-reduced contiguous storage.
!>
!> `descriptor` must be the exact collectively validated companion descriptor
!> exported with `values`.  Rebinding saved raw values to a fresh descriptor is
!> forbidden by the host contract.
   SUBROUTINE cp2k_exchange_partial_reverse_unpack(buffer, descriptor, values, &
      & error)
      TYPE(cp2k_exchange_partial_reverse_type), &
         INTENT(out)                                     :: buffer
      TYPE(cp2k_exchange_partial_descriptor_type), &
         INTENT(in)                                      :: descriptor
      COMPLEX(wp), INTENT(in)                            :: values(:)
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      INTEGER(int64)                                     :: expected, nblock, nonsite
      LOGICAL                                            :: extents_ok

      IF (.NOT. partial_descriptor_is_valid(descriptor) &
         & .OR. descriptor%phase /= 2) THEN
         CALL fatal_error(error, "inconsistent g-xTB exchange reverse descriptor")
         RETURN
      END IF
      CALL partial_reverse_extents(descriptor%nao, descriptor%batch_count, &
         & descriptor%nspin, nblock, nonsite, expected, extents_ok)
      IF (.NOT. extents_ok) THEN
         CALL fatal_error(error, &
            & "g-xTB exchange reverse payload exceeds indexing limit")
         RETURN
      END IF
      IF (SIZE(values, kind=int64) /= expected) THEN
         CALL fatal_error(error, "inconsistent g-xTB exchange reverse payload")
         RETURN
      END IF
      IF (.NOT. ALL(ieee_is_finite(REAL(values, wp))) &
         & .OR. .NOT. ALL(ieee_is_finite(AIMAG(values)))) THEN
         CALL fatal_error(error, "non-finite g-xTB exchange reverse payload")
         RETURN
      END IF
      ALLOCATE (buffer%values, source=values)
      buffer%image_first = descriptor%batch_first
      buffer%image_count = descriptor%batch_count
      buffer%abi_version = cp2k_exchange_partial_abi_version
      buffer%descriptor_signature = partial_descriptor_signature(descriptor)
   END SUBROUTINE cp2k_exchange_partial_reverse_unpack

!> Import globally summed reverse sources before the coupled adjoint kernel.
   SUBROUTINE cp2k_exchange_partial_import_reverse(stream, descriptor, &
      & owner_count, buffer, calc, mol, ecache, plan, error)
      TYPE(cp2k_exchange_stream_type), INTENT(inout)     :: stream
      TYPE(cp2k_exchange_partial_descriptor_type), &
         INTENT(in)                                      :: descriptor
      INTEGER, INTENT(in)                                :: owner_count(:)
      TYPE(cp2k_exchange_partial_reverse_type), &
         INTENT(in)                                      :: buffer
      TYPE(xtb_calculator), INTENT(in)                   :: calc
      TYPE(structure_type), INTENT(in)                   :: mol
      TYPE(container_cache), INTENT(in)                  :: ecache
      TYPE(cp2k_stream_plan_type), INTENT(in) :: plan
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      INTEGER(int64)                                     :: first, last, nblock, nonsite, total
      LOGICAL                                            :: extents_ok
      TYPE(cp2k_exchange_partial_descriptor_type)        :: current

      IF (.NOT. stream%active .OR. .NOT. stream%partial_active &
         & .OR. .NOT. stream%reverse_transaction) THEN
         CALL fatal_error(error, &
            & "g-xTB exchange reverse partial stream is not active")
         RETURN
      END IF
      IF (stream%partial_reverse_imported) THEN
         CALL fatal_error(error, "g-xTB exchange reverse payload was imported")
         RETURN
      END IF
      IF (stream%reverse_batch_applied .OR. stream%reverse_applied) THEN
         CALL fatal_error(error, &
            & "g-xTB exchange reverse payload import followed kernel apply")
         RETURN
      END IF
      IF (stream%partial_compact_image_storage) THEN
         CALL fatal_error(error, &
            & "range-local g-xTB importer requires a reverse range payload")
         RETURN
      END IF
      IF (.NOT. exchange_stream_state_matches(stream, calc, mol, ecache, plan)) THEN
         CALL fatal_error(error, &
            & "g-xTB exchange partial state was invalidated before import")
         RETURN
      END IF
      CALL get_partial_descriptor(stream, 2, current)
      IF (.NOT. cp2k_exchange_partial_descriptors_match(current, descriptor)) THEN
         CALL fatal_error(error, "g-xTB exchange partial descriptor mismatch")
         RETURN
      END IF
      IF (SIZE(owner_count) /= stream%nk) THEN
         CALL fatal_error(error, "inconsistent g-xTB exchange owner-count size")
         RETURN
      END IF
      IF (ANY(owner_count < 0)) THEN
         CALL fatal_error(error, "negative g-xTB exchange owner count")
         RETURN
      END IF
      IF (ANY(owner_count > 1)) THEN
         CALL fatal_error(error, "duplicate g-xTB exchange k-point ownership")
         RETURN
      END IF
      IF (ANY(owner_count == 0)) THEN
         CALL fatal_error(error, "missing g-xTB exchange k-point ownership")
         RETURN
      END IF
      CALL partial_reverse_extents(stream%nao, stream%batch_count, stream%nspin, &
         & nblock, nonsite, total, extents_ok)
      IF (.NOT. extents_ok) THEN
         CALL fatal_error(error, &
            & "g-xTB exchange reverse payload exceeds indexing limit")
         RETURN
      END IF
      IF (buffer%abi_version /= cp2k_exchange_partial_abi_version &
         & .OR. .NOT. ALLOCATED(buffer%values)) THEN
         CALL fatal_error(error, "inconsistent g-xTB exchange reverse payload")
         RETURN
      END IF
      IF (SIZE(buffer%values, kind=int64) /= total &
         & .OR. buffer%image_first /= descriptor%batch_first &
         & .OR. buffer%image_count /= descriptor%batch_count &
         & .OR. ANY(buffer%descriptor_signature /= &
         & partial_descriptor_signature(descriptor))) THEN
         CALL fatal_error(error, "inconsistent g-xTB exchange reverse payload")
         RETURN
      END IF
      IF (.NOT. ALL(ieee_is_finite(REAL(buffer%values, wp))) &
         & .OR. .NOT. ALL(ieee_is_finite(AIMAG(buffer%values)))) THEN
         CALL fatal_error(error, "non-finite g-xTB exchange reverse payload")
         RETURN
      END IF

      first = 1_int64
      last = nblock
      stream%amat_r(:, :, 1:stream%batch_count, :) = &
         & RESHAPE(buffer%values(first:last), &
         & SHAPE(stream%amat_r(:, :, 1:stream%batch_count, :)))
      first = last + 1_int64
      last = last + nblock
      stream%cmat_r(:, :, 1:stream%batch_count, :) = &
         & RESHAPE(buffer%values(first:last), &
         & SHAPE(stream%cmat_r(:, :, 1:stream%batch_count, :)))
      first = last + 1_int64
      last = last + nblock
      stream%vmat_r(:, :, 1:stream%batch_count, :) = &
         & RESHAPE(buffer%values(first:last), &
         & SHAPE(stream%vmat_r(:, :, 1:stream%batch_count, :)))
      first = last + 1_int64
      last = last + nblock
      stream%amat_partner_r(:, :, 1:stream%batch_count, :) = &
         & RESHAPE(buffer%values(first:last), &
         & SHAPE(stream%amat_partner_r(:, :, 1:stream%batch_count, :)))
      first = last + 1_int64
      last = last + nblock
      stream%reverse_amat_r(:, :, 1:stream%batch_count, :) = &
         & RESHAPE(buffer%values(first:last), &
         & SHAPE(stream%reverse_amat_r(:, :, 1:stream%batch_count, :)))
      first = last + 1_int64
      last = last + nblock
      stream%reverse_vmat_r(:, :, 1:stream%batch_count, :) = &
         & RESHAPE(buffer%values(first:last), &
         & SHAPE(stream%reverse_vmat_r(:, :, 1:stream%batch_count, :)))
      first = last + 1_int64
      last = last + nonsite
      stream%gdiagP = RESHAPE(buffer%values(first:last), SHAPE(stream%gdiagP))
      first = last + 1_int64
      last = last + nonsite
      stream%gdiagSP = RESHAPE(buffer%values(first:last), SHAPE(stream%gdiagSP))
      first = last + 1_int64
      last = last + nonsite
      stream%gdiagSPS = RESHAPE(buffer%values(first:last), SHAPE(stream%gdiagSPS))
      first = last + 1_int64
      last = last + nonsite
      stream%reverse_bdiagSP = RESHAPE(buffer%values(first:last), &
         & SHAPE(stream%reverse_bdiagSP))
      first = last + 1
      last = last + nonsite
      stream%reverse_bdiagSPS = RESHAPE(buffer%values(first:last), &
         & SHAPE(stream%reverse_bdiagSPS))
      stream%seen = .TRUE.
      IF (stream%transform_backend /= bvk_transform_dense) THEN
         stream%transform_pair_first = 0
         stream%transform_pair_count = 0
         stream%transform_k = (0.0_wp, 0.0_wp)
      END IF
      stream%partial_reverse_imported = .TRUE.
   END SUBROUTINE cp2k_exchange_partial_import_reverse

!> Import the compact local image source plus replicated onsite adjoint data.
   SUBROUTINE cp2k_exchange_partial_import_reverse_range(stream, descriptor, &
      & owner_count, buffer, calc, mol, ecache, plan, error)
      TYPE(cp2k_exchange_stream_type), INTENT(inout)     :: stream
      TYPE(cp2k_exchange_partial_descriptor_type), &
         INTENT(in)                                      :: descriptor
      INTEGER, INTENT(in)                                :: owner_count(:)
      TYPE(cp2k_exchange_partial_reverse_type), &
         INTENT(in)                                      :: buffer
      TYPE(xtb_calculator), INTENT(in)                   :: calc
      TYPE(structure_type), INTENT(in)                   :: mol
      TYPE(container_cache), INTENT(in)                  :: ecache
      TYPE(cp2k_stream_plan_type), INTENT(in) :: plan
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      INTEGER                                            :: first_global, first_local, last_local, &
                                                            local_count
      INTEGER(int64)                                     :: first, last, nblock, nonsite, total
      LOGICAL                                            :: extents_ok
      TYPE(cp2k_exchange_partial_descriptor_type)        :: current

      IF (.NOT. stream%active .OR. .NOT. stream%partial_active &
         & .OR. .NOT. stream%partial_compact_image_storage &
         & .OR. .NOT. stream%reverse_transaction) THEN
         CALL fatal_error(error, &
            & "g-xTB exchange reverse range importer is not active")
         RETURN
      END IF
      IF (stream%partial_reverse_imported) THEN
         CALL fatal_error(error, "g-xTB exchange reverse payload was imported")
         RETURN
      END IF
      IF (stream%reverse_batch_applied .OR. stream%reverse_applied) THEN
         CALL fatal_error(error, &
            & "g-xTB exchange reverse range import followed kernel apply")
         RETURN
      END IF
      IF (.NOT. exchange_stream_state_matches(stream, calc, mol, ecache, plan)) THEN
         CALL fatal_error(error, &
            & "g-xTB exchange partial state was invalidated before range import")
         RETURN
      END IF
      CALL get_partial_descriptor(stream, 2, current)
      IF (.NOT. cp2k_exchange_partial_descriptors_match(current, descriptor)) THEN
         CALL fatal_error(error, "g-xTB exchange partial descriptor mismatch")
         RETURN
      END IF
      IF (SIZE(owner_count) /= stream%nk .OR. ANY(owner_count /= 1)) THEN
         CALL fatal_error(error, &
            & "incomplete g-xTB exchange k-point ownership for range import")
         RETURN
      END IF
      CALL partial_local_batch_bounds(stream, first_local, last_local, &
         & local_count, first_global)
      CALL partial_reverse_extents(stream%nao, local_count, stream%nspin, &
         & nblock, nonsite, total, extents_ok)
      IF (.NOT. extents_ok .OR. local_count > SIZE(stream%amat_r, 3)) THEN
         CALL fatal_error(error, &
            & "g-xTB exchange reverse range exceeds local storage")
         RETURN
      END IF
      IF (buffer%abi_version /= cp2k_exchange_partial_abi_version &
         & .OR. .NOT. ALLOCATED(buffer%values) &
         & .OR. buffer%image_first /= MERGE(first_global, 0, local_count > 0) &
         & .OR. buffer%image_count /= local_count &
         & .OR. SIZE(buffer%values, kind=int64) /= total &
         & .OR. ANY(buffer%descriptor_signature /= &
         & partial_descriptor_signature(descriptor))) THEN
         CALL fatal_error(error, "inconsistent g-xTB exchange reverse range")
         RETURN
      END IF
      IF (.NOT. ALL(ieee_is_finite(REAL(buffer%values, wp))) &
         & .OR. .NOT. ALL(ieee_is_finite(AIMAG(buffer%values)))) THEN
         CALL fatal_error(error, "non-finite g-xTB exchange reverse range")
         RETURN
      END IF

      IF (local_count > 0) THEN
         first = 1_int64
         last = nblock
         stream%amat_r(:, :, 1:local_count, :) = &
            & RESHAPE(buffer%values(first:last), &
            & SHAPE(stream%amat_r(:, :, 1:local_count, :)))
         first = last + 1_int64
         last = last + nblock
         stream%cmat_r(:, :, 1:local_count, :) = &
            & RESHAPE(buffer%values(first:last), &
            & SHAPE(stream%cmat_r(:, :, 1:local_count, :)))
         first = last + 1_int64
         last = last + nblock
         stream%vmat_r(:, :, 1:local_count, :) = &
            & RESHAPE(buffer%values(first:last), &
            & SHAPE(stream%vmat_r(:, :, 1:local_count, :)))
         first = last + 1_int64
         last = last + nblock
         stream%amat_partner_r(:, :, 1:local_count, :) = &
            & RESHAPE(buffer%values(first:last), &
            & SHAPE(stream%amat_partner_r(:, :, 1:local_count, :)))
         first = last + 1_int64
         last = last + nblock
         stream%reverse_amat_r(:, :, 1:local_count, :) = &
            & RESHAPE(buffer%values(first:last), &
            & SHAPE(stream%reverse_amat_r(:, :, 1:local_count, :)))
         first = last + 1_int64
         last = last + nblock
         stream%reverse_vmat_r(:, :, 1:local_count, :) = &
            & RESHAPE(buffer%values(first:last), &
            & SHAPE(stream%reverse_vmat_r(:, :, 1:local_count, :)))
      END IF
      first = 6_int64*nblock + 1_int64
      last = first + nonsite - 1_int64
      stream%gdiagP = RESHAPE(buffer%values(first:last), SHAPE(stream%gdiagP))
      first = last + 1_int64
      last = last + nonsite
      stream%gdiagSP = RESHAPE(buffer%values(first:last), SHAPE(stream%gdiagSP))
      first = last + 1_int64
      last = last + nonsite
      stream%gdiagSPS = RESHAPE(buffer%values(first:last), SHAPE(stream%gdiagSPS))
      first = last + 1_int64
      last = last + nonsite
      stream%reverse_bdiagSP = RESHAPE(buffer%values(first:last), &
         & SHAPE(stream%reverse_bdiagSP))
      first = last + 1_int64
      last = last + nonsite
      stream%reverse_bdiagSPS = RESHAPE(buffer%values(first:last), &
         & SHAPE(stream%reverse_bdiagSPS))
      stream%seen = .TRUE.
      IF (stream%transform_backend /= bvk_transform_dense) THEN
         stream%transform_pair_first = 0
         stream%transform_pair_count = 0
         stream%transform_k = (0.0_wp, 0.0_wp)
      END IF
      stream%partial_reverse_imported = .TRUE.
   END SUBROUTINE cp2k_exchange_partial_import_reverse_range

!> Restrict one replicated importer to a disjoint contiguous image range.
!>
!> The host first reduces and imports the complete bounded pre-kernel source on
!> every participating importer.  It may then assign each importer one
!> contiguous subset of the global image mesh.  The nonlinear image kernel,
!> its inverse transform, and the direct force/stress response are evaluated
!> only for this subset.  A zero-length range is represented by `(0, 0)` and
!> permits more importers than BvK images.  A nonempty range must be contained
!> in the image window selected when the partial stream began.  The assignment
!> is fixed for the complete forward/reverse transaction; changing it between
!> phases requires a fresh reverse-only stream.  CP2K owns the collective
!> exact-cover proof and must sum the partial Fock, shell-potential, energy,
!> force, stress, and overlap-adjoint results over all importers.
!>
!> Image-range distribution and disjoint k-point pull ownership are mutually
!> exclusive: an image-local response contributes to every k point and must be
!> reduced before a k-point result is complete.
   SUBROUTINE cp2k_exchange_partial_set_image_range(stream, first_image, &
      & image_count, error)
      TYPE(cp2k_exchange_stream_type), INTENT(inout)     :: stream
      INTEGER, INTENT(in)                                :: first_image, image_count
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      IF (.NOT. stream%active .OR. .NOT. stream%partial_active) THEN
         CALL fatal_error(error, "g-xTB exchange partial stream is not active")
         RETURN
      END IF
      IF (stream%mode /= cp2k_exchange_stream_batched) THEN
         CALL fatal_error(error, &
            & "distributed image ownership requires batched exchange mode")
         RETURN
      END IF
      IF (stream%reverse_transaction) THEN
         IF (stream%reverse_batch_applied) THEN
            CALL fatal_error(error, &
               & "g-xTB exchange partial image range followed kernel apply")
            RETURN
         END IF
         IF (.NOT. stream%partial_reverse_imported) THEN
            CALL fatal_error(error, &
               & "reverse partial image range requires an imported payload")
            RETURN
         END IF
      ELSE
         IF (stream%batch_applied) THEN
            CALL fatal_error(error, &
               & "g-xTB exchange partial image range followed kernel apply")
            RETURN
         END IF
         IF (.NOT. stream%partial_forward_imported) THEN
            CALL fatal_error(error, &
               & "forward partial image range requires an imported payload")
            RETURN
         END IF
      END IF
      IF (stream%partial_generation /= 1) THEN
         CALL fatal_error(error, &
            & "g-xTB exchange partial image range followed an earlier batch")
         RETURN
      END IF
      IF (stream%partial_image_distribution_set) THEN
         CALL fatal_error(error, &
            & "g-xTB exchange partial image range was already set")
         RETURN
      END IF
      IF (stream%partial_pull_distribution_set) THEN
         CALL fatal_error(error, &
            & "distributed images cannot use disjoint k-point pulls")
         RETURN
      END IF
      IF (image_count < 0) THEN
         CALL fatal_error(error, "negative g-xTB exchange partial image count")
         RETURN
      END IF
      IF (image_count == 0) THEN
         IF (first_image /= 0) THEN
            CALL fatal_error(error, &
               & "empty g-xTB exchange partial image range must start at zero")
            RETURN
         END IF
         stream%partial_image_first = 0
         stream%partial_image_last = 0
      ELSE
         IF (first_image < stream%image_first &
            & .OR. first_image > stream%batch_last) THEN
            CALL fatal_error(error, &
               & "invalid g-xTB exchange partial image range")
            RETURN
         END IF
         IF (image_count > stream%batch_last - first_image + 1) THEN
            CALL fatal_error(error, &
               & "invalid g-xTB exchange partial image range")
            RETURN
         END IF
         stream%partial_image_first = first_image
         stream%partial_image_last = first_image + (image_count - 1)
      END IF
      stream%partial_image_distribution_set = .TRUE.
   END SUBROUTINE cp2k_exchange_partial_set_image_range

!> Return the importer-local contiguous slice of the current global batch.
   pure SUBROUTINE partial_local_batch_bounds(stream, first_local, last_local, &
      & local_count, first_global)
      TYPE(cp2k_exchange_stream_type), INTENT(in)        :: stream
      INTEGER, INTENT(out)                               :: first_local, last_local, local_count, &
                                                            first_global

      INTEGER                                            :: last_global

      IF (.NOT. stream%partial_image_distribution_set) THEN
         first_local = 1
         last_local = stream%batch_count
         local_count = stream%batch_count
         first_global = stream%batch_first
         RETURN
      END IF
      IF (stream%partial_image_first == 0) THEN
         first_local = 1
         last_local = 0
         local_count = 0
         first_global = 0
         RETURN
      END IF

      first_global = MAX(stream%batch_first, stream%partial_image_first)
      last_global = MIN(stream%batch_first + (stream%batch_count - 1), &
         & stream%partial_image_last)
      IF (last_global < first_global) THEN
         first_local = 1
         last_local = 0
         local_count = 0
         first_global = 0
      ELSE
         local_count = last_global - first_global + 1
         IF (stream%partial_compact_image_storage) THEN
            first_local = 1
            last_local = local_count
         ELSE
            first_local = first_global - stream%batch_first + 1
            last_local = last_global - stream%batch_first + 1
         END IF
      END IF
   END SUBROUTINE partial_local_batch_bounds

!> Select a disjoint distributed pull partition for replicated importers.
!>
!> `owner_rank(ik)` is the one-based contributor that reconstructs Bloch block
!> `ik`; values must cover the complete mesh and lie in
!> `1:contributor_count`.  Every importer must receive the same host-broadcast
!> owner map, while `contributor_rank` identifies the local importer.  Exactly
!> `result_owner_rank` may retrieve the replicated scalar shell-potential/
!> energy or direct force/stress result.  This routine is MPI-free: CP2K owns
!> the broadcast and must first establish exact descriptor agreement.
!>
!> The operation is optional.  Without it, the legacy designated-root path
!> remains strict and requires every Bloch block plus the scalar result to be
!> pulled from the single importer.
   SUBROUTINE cp2k_exchange_partial_set_pull_ownership(stream, owner_rank, &
      & contributor_rank, contributor_count, result_owner_rank, error)
      TYPE(cp2k_exchange_stream_type), INTENT(inout)     :: stream
      INTEGER, INTENT(in)                                :: owner_rank(:), contributor_rank, &
                                                            contributor_count, result_owner_rank
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      IF (.NOT. stream%active .OR. .NOT. stream%partial_active) THEN
         CALL fatal_error(error, "g-xTB exchange partial stream is not active")
         RETURN
      END IF
      IF (stream%mode /= cp2k_exchange_stream_batched) THEN
         CALL fatal_error(error, &
            & "distributed partial pulls require batched exchange mode")
         RETURN
      END IF
      IF (stream%partial_pull_distribution_set) THEN
         CALL fatal_error(error, &
            & "g-xTB exchange partial pull ownership was already set")
         RETURN
      END IF
      IF (stream%partial_image_distribution_set) THEN
         CALL fatal_error(error, &
            & "disjoint k-point pulls cannot use distributed images")
         RETURN
      END IF
      IF (stream%reverse_transaction) THEN
         IF (.NOT. stream%partial_reverse_imported) THEN
            CALL fatal_error(error, &
               & "reverse partial pull ownership requires an imported payload")
            RETURN
         END IF
         IF (ANY(stream%reverse_pulled)) THEN
            CALL fatal_error(error, &
               & "reverse partial pull ownership followed a pull")
            RETURN
         END IF
      ELSE
         IF (.NOT. stream%partial_forward_imported) THEN
            CALL fatal_error(error, &
               & "forward partial pull ownership requires an imported payload")
            RETURN
         END IF
         IF (ANY(stream%batch_pulled)) THEN
            CALL fatal_error(error, &
               & "forward partial pull ownership followed a pull")
            RETURN
         END IF
      END IF
      IF (SIZE(owner_rank) /= stream%nk) THEN
         CALL fatal_error(error, &
            & "inconsistent g-xTB exchange partial pull-owner size")
         RETURN
      END IF
      IF (contributor_count < 1 &
         & .OR. contributor_rank < 1 &
         & .OR. contributor_rank > contributor_count &
         & .OR. result_owner_rank < 1 &
         & .OR. result_owner_rank > contributor_count) THEN
         CALL fatal_error(error, &
            & "invalid g-xTB exchange partial contributor identity")
         RETURN
      END IF
      IF (ANY(owner_rank < 1) .OR. ANY(owner_rank > contributor_count)) THEN
         CALL fatal_error(error, &
            & "invalid g-xTB exchange partial pull-owner map")
         RETURN
      END IF

      ALLOCATE (stream%partial_pull_mask(stream%nk), &
         & source=owner_rank == contributor_rank)
      stream%partial_result_owner = contributor_rank == result_owner_rank
      stream%partial_pull_distribution_set = .TRUE.
   END SUBROUTINE cp2k_exchange_partial_set_pull_ownership

!> Reuse an exported contributor for the next sequential image batch.
!>
!> A contributor has no permission to run the nonlinear kernel.  This
!> operation therefore advances only its linear accumulation buffers after
!> the host has copied the exported payload.  The importing accumulator uses
!> the ordinary batch-advance operation after apply and pull; both paths then
!> arrive at the same next descriptor generation.
   SUBROUTINE cp2k_exchange_partial_contributor_advance(stream, done, error)
      TYPE(cp2k_exchange_stream_type), INTENT(inout)     :: stream
      LOGICAL, INTENT(out)                               :: done
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      LOGICAL                                            :: reverse_phase

      done = .FALSE.
      IF (.NOT. stream%active .OR. .NOT. stream%partial_active) THEN
         CALL fatal_error(error, "g-xTB exchange partial stream is not active")
         RETURN
      END IF
      IF (stream%mode /= cp2k_exchange_stream_batched) THEN
         CALL fatal_error(error, &
            & "g-xTB exchange partial stream is not in batched mode")
         RETURN
      END IF
      reverse_phase = stream%reverse_transaction
      IF (reverse_phase) THEN
         IF (.NOT. stream%partial_reverse_exported &
            & .OR. stream%partial_reverse_imported) THEN
            CALL fatal_error(error, &
               & "reverse contributor advance requires an unimported export")
            RETURN
         END IF
      ELSE
         IF (.NOT. stream%partial_forward_exported &
            & .OR. stream%partial_forward_imported) THEN
            CALL fatal_error(error, &
               & "forward contributor advance requires an unimported export")
            RETURN
         END IF
      END IF
      IF (stream%batch_count > stream%batch_last - stream%batch_first) THEN
         done = .TRUE.
         RETURN
      END IF

      stream%batch_first = stream%batch_first + stream%batch_count
      stream%batch_count = MIN(stream%batch_size, &
         & stream%batch_last - stream%batch_first + 1)
      stream%seen = .FALSE.
      IF (stream%transform_backend /= bvk_transform_dense) THEN
         stream%transform_pair_first = 1
         stream%transform_pair_count = MIN( &
            & stream%transform_pair_batch_size, stream%nao*stream%nao)
         stream%transform_k = (0.0_wp, 0.0_wp)
      END IF
      stream%partial_generation = stream%partial_generation + 1
      stream%amat_r = (0.0_wp, 0.0_wp)
      stream%cmat_r = (0.0_wp, 0.0_wp)
      stream%vmat_r = (0.0_wp, 0.0_wp)
      stream%amat_partner_r = (0.0_wp, 0.0_wp)
      IF (reverse_phase) THEN
         stream%partial_reverse_exported = .FALSE.
         stream%reverse_amat_r = (0.0_wp, 0.0_wp)
         stream%reverse_vmat_r = (0.0_wp, 0.0_wp)
         stream%gdiagP = (0.0_wp, 0.0_wp)
         stream%gdiagSP = (0.0_wp, 0.0_wp)
         stream%gdiagSPS = (0.0_wp, 0.0_wp)
         stream%reverse_bdiagSP = (0.0_wp, 0.0_wp)
         stream%reverse_bdiagSPS = (0.0_wp, 0.0_wp)
      ELSE
         stream%partial_forward_exported = .FALSE.
      END IF
   END SUBROUTINE cp2k_exchange_partial_contributor_advance

!> Abandon an exported contributor after the host has copied its payload.
!>
!> This is the only early-close operation.  It is intentionally restricted to
!> partial streams so ordinary stream transactions retain their strict result
!> and pull completion requirements.
   SUBROUTINE cp2k_exchange_partial_discard(stream, error)
      TYPE(cp2k_exchange_stream_type), INTENT(inout)     :: stream
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      IF (.NOT. stream%active .OR. .NOT. stream%partial_active) THEN
         CALL fatal_error(error, "g-xTB exchange partial stream is not active")
         RETURN
      END IF
      IF (.NOT. stream%partial_forward_exported .AND. &
         & .NOT. stream%partial_reverse_exported) THEN
         CALL fatal_error(error, &
            & "only an exported g-xTB exchange contributor can be discarded")
         RETURN
      END IF
      IF (stream%partial_forward_imported .OR. &
         & stream%partial_reverse_imported) THEN
         CALL fatal_error(error, &
            & "an imported g-xTB exchange accumulator cannot be discarded")
         RETURN
      END IF
      CALL clear_exchange_stream(stream)
   END SUBROUTINE cp2k_exchange_partial_discard

!> Initialize the bounded reverse state after mesh/model validation.
   SUBROUTINE initialize_exchange_stream_reverse(stream, mol, error)
      TYPE(cp2k_exchange_stream_type), INTENT(inout)     :: stream
      TYPE(structure_type), INTENT(in)                   :: mol
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      INTEGER                                            :: capacity

      IF (stream%reverse_transaction .OR. stream%reverse_applied) THEN
         CALL fatal_error(error, &
            & "g-xTB exchange stream reverse transaction already started")
         RETURN
      END IF
      IF (stream%mode /= cp2k_exchange_stream_reduced .AND. &
         & stream%mode /= cp2k_exchange_stream_batched) THEN
         CALL fatal_error(error, &
            & "streamed g-xTB exchange reverse requires reduced or batched mode")
         RETURN
      END IF
      IF (ALLOCATED(stream%partial_pull_mask)) &
         & DEALLOCATE (stream%partial_pull_mask)
      stream%partial_pull_distribution_set = .FALSE.
      stream%partial_result_owner = .FALSE.

      IF (stream%mode == cp2k_exchange_stream_reduced) THEN
         stream%batch_size = stream%nk
         stream%image_first = 1
         stream%batch_first = 1
         stream%batch_count = stream%nk
         stream%batch_last = stream%nk
      ELSE
         stream%batch_first = stream%image_first
         stream%batch_count = MIN(stream%batch_size, &
            & stream%batch_last - stream%image_first + 1)
      END IF
      capacity = SIZE(stream%amat_r, 3)
      IF (stream%transform_backend /= bvk_transform_dense) THEN
         stream%transform_pair_first = 1
         stream%transform_pair_count = MIN( &
            & stream%transform_pair_batch_size, stream%nao*stream%nao)
         stream%transform_field_count = 6
         IF (ALLOCATED(stream%transform_k)) DEALLOCATE (stream%transform_k)
         ALLOCATE (stream%transform_k(stream%transform_pair_batch_size, &
            & stream%nk, stream%nspin, stream%transform_field_count), &
            & source=(0.0_wp, 0.0_wp))
      END IF

      IF (.NOT. ALLOCATED(stream%amat_r)) &
         & ALLOCATE (stream%amat_r(stream%nao, stream%nao, capacity, &
         & stream%nspin))
      IF (.NOT. ALLOCATED(stream%cmat_r)) &
         & ALLOCATE (stream%cmat_r(stream%nao, stream%nao, capacity, &
         & stream%nspin))
      IF (.NOT. ALLOCATED(stream%vmat_r)) &
         & ALLOCATE (stream%vmat_r(stream%nao, stream%nao, capacity, &
         & stream%nspin))
      IF (.NOT. ALLOCATED(stream%amat_partner_r)) &
         & ALLOCATE (stream%amat_partner_r(stream%nao, stream%nao, capacity, &
         & stream%nspin))
      ALLOCATE (stream%reverse_amat_r(stream%nao, stream%nao, capacity, &
         & stream%nspin), stream%reverse_vmat_r(stream%nao, stream%nao, &
         & capacity, stream%nspin), source=(0.0_wp, 0.0_wp))
      IF (.NOT. ALLOCATED(stream%gdiagP)) &
         & ALLOCATE (stream%gdiagP(stream%nao, stream%nspin))
      IF (.NOT. ALLOCATED(stream%gdiagSP)) &
         & ALLOCATE (stream%gdiagSP(stream%nao, stream%nspin))
      IF (.NOT. ALLOCATED(stream%gdiagSPS)) &
         & ALLOCATE (stream%gdiagSPS(stream%nao, stream%nspin))
      ALLOCATE (stream%reverse_bdiagSP(stream%nao, stream%nspin), &
         & stream%reverse_bdiagSPS(stream%nao, stream%nspin), &
         & source=(0.0_wp, 0.0_wp))
      ALLOCATE (stream%reverse_gradient(3, mol%nat), &
         & stream%reverse_sigma(3, 3), source=0.0_wp)
      IF (.NOT. ALLOCATED(stream%density_signature)) &
         & ALLOCATE (stream%density_signature(stream%nk), source=0_int64)
      IF (.NOT. ALLOCATED(stream%overlap_signature)) &
         & ALLOCATE (stream%overlap_signature(stream%nk), source=0_int64)
      IF (.NOT. ALLOCATED(stream%block_signature_set)) &
         & ALLOCATE (stream%block_signature_set(stream%nk), source=.FALSE.)

      stream%amat_r = (0.0_wp, 0.0_wp)
      stream%cmat_r = (0.0_wp, 0.0_wp)
      stream%vmat_r = (0.0_wp, 0.0_wp)
      stream%amat_partner_r = (0.0_wp, 0.0_wp)
      stream%gdiagP = (0.0_wp, 0.0_wp)
      stream%gdiagSP = (0.0_wp, 0.0_wp)
      stream%gdiagSPS = (0.0_wp, 0.0_wp)
      stream%seen = .FALSE.
      stream%reverse_pulled = .FALSE.
      stream%forward_pulled = .TRUE.
      stream%batch_result_pulled = .TRUE.
      stream%reverse_batch_applied = .FALSE.
      stream%reverse_result_pulled = .FALSE.
      IF (stream%partial_active) THEN
         ! Forward partial importers can reconstruct Fock blocks for k points
         ! whose densities they never owned; those pulls establish only overlap
         ! signatures.  Preserve complete P/S signatures for locally owned
         ! forward blocks, but let reverse push/pull establish a complete pair
         ! for foreign blocks.  With the same stream objects, the collective set
         ! of retained owners therefore enforces forward-to-reverse identity for
         ! every k point without exporting additional signature metadata.
         WHERE (.NOT. stream%partial_forward_owned)
            stream%block_signature_set = .FALSE.
            stream%density_signature = 0_int64
            stream%overlap_signature = 0_int64
         END WHERE
         stream%partial_generation = 1
         stream%partial_reverse_exported = .FALSE.
         stream%partial_reverse_imported = .FALSE.
      END IF
      stream%reverse_transaction = .TRUE.
      stream%reverse_applied = .FALSE.
      stream%peak_complex_elements = MAX(stream%peak_complex_elements, &
         & exchange_stream_complex_elements(stream) &
         & + 10_int64*INT(stream%nao, int64)**2)
      IF (stream%transform_backend /= bvk_transform_dense) THEN
         stream%peak_complex_elements = MAX(stream%peak_complex_elements, &
            & exchange_stream_complex_elements(stream) &
            & + 3_int64*INT(stream%transform_pair_batch_size, int64) &
            & *INT(stream%nk, int64))
      END IF
   END SUBROUTINE initialize_exchange_stream_reverse

!> Enter a bounded reverse transaction after a reduced/batched forward apply.
!>
!> Forward Fock pulls are intentionally waived: the reverse sweep recomputes
!> only the image-local responses needed by the overlap adjoint and therefore
!> does not force a gradient-only caller to reconstruct unused Fock blocks.
   SUBROUTINE cp2k_exchange_stream_reverse_begin(stream, mol, error)
      TYPE(cp2k_exchange_stream_type), INTENT(inout)     :: stream
      TYPE(structure_type), INTENT(in)                   :: mol
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      IF (.NOT. stream%active) THEN
         CALL fatal_error(error, "g-xTB exchange stream is not active")
         RETURN
      END IF
      IF (.NOT. stream%applied) THEN
         CALL fatal_error(error, &
            & "g-xTB exchange stream reverse requires forward apply")
         RETURN
      END IF
      CALL initialize_exchange_stream_reverse(stream, mol, error)
   END SUBROUTINE cp2k_exchange_stream_reverse_begin

!> Push one unweighted density/overlap Bloch block into an active stream.
   SUBROUTINE cp2k_exchange_stream_push(stream, ik, density, overlap, error)
      TYPE(cp2k_exchange_stream_type), INTENT(inout)     :: stream
      INTEGER, INTENT(in)                                :: ik
      COMPLEX(wp), INTENT(in)                            :: density(:, :, :), overlap(:, :)
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      COMPLEX(wp)                                        :: phase
      COMPLEX(wp), ALLOCATABLE                           :: amat(:, :), cmat(:, :)
      INTEGER                                            :: iao, ibatch, icell, ilocal, &
                                                            inverse_cell, ipair, jao, spin
      INTEGER(int64)                                     :: density_signature, overlap_signature
      REAL(wp)                                           :: angle

      IF (.NOT. stream%active) THEN
         CALL fatal_error(error, "g-xTB exchange stream is not active")
         RETURN
      END IF
      IF (stream%partial_active .AND. stream%partial_forward_imported) THEN
         CALL fatal_error(error, "g-xTB exchange partial push followed import")
         RETURN
      END IF
      IF (stream%partial_active .AND. stream%partial_forward_exported) THEN
         CALL fatal_error(error, "g-xTB exchange partial push followed export")
         RETURN
      END IF
      IF (stream%applied) THEN
         CALL fatal_error(error, "g-xTB exchange stream was already applied")
         RETURN
      END IF
      IF (stream%mode == cp2k_exchange_stream_batched &
         & .AND. stream%batch_applied) THEN
         CALL fatal_error(error, &
            & "g-xTB exchange image batch was already applied")
         RETURN
      END IF
      IF (stream%transform_backend /= bvk_transform_dense .AND. &
         & stream%transform_pair_first < 1) THEN
         CALL fatal_error(error, &
            & "g-xTB compact forward transform batch is already complete")
         RETURN
      END IF
      IF (ik < 1 .OR. ik > stream%nk) THEN
         CALL fatal_error(error, "g-xTB exchange stream block index is out of range")
         RETURN
      END IF
      IF (stream%seen(ik)) THEN
         CALL fatal_error(error, "duplicate g-xTB exchange stream block")
         RETURN
      END IF
      IF (SIZE(density, 1) /= stream%nao &
         & .OR. SIZE(density, 2) /= stream%nao &
         & .OR. SIZE(density, 3) /= stream%nspin &
         & .OR. SIZE(overlap, 1) /= stream%nao &
         & .OR. SIZE(overlap, 2) /= stream%nao) THEN
         CALL fatal_error(error, &
            & "inconsistent g-xTB exchange stream block dimensions")
         RETURN
      END IF

      IF (stream%mode == cp2k_exchange_stream_batched) THEN
         density_signature = exchange_hash_complex3(density)
         overlap_signature = exchange_hash_complex2(overlap)
         IF (stream%block_signature_set(ik)) THEN
            IF (density_signature /= stream%density_signature(ik)) THEN
               CALL fatal_error(error, &
                  & "changed g-xTB exchange image-batch density block")
               RETURN
            END IF
            IF (overlap_signature /= stream%overlap_signature(ik)) THEN
               CALL fatal_error(error, &
                  & "changed g-xTB exchange image-batch overlap block")
               RETURN
            END IF
         ELSE
            stream%density_signature(ik) = density_signature
            stream%overlap_signature(ik) = overlap_signature
            stream%block_signature_set(ik) = .TRUE.
         END IF
      END IF

      IF (stream%mode == cp2k_exchange_stream_oracle) THEN
         stream%density(:, :, :, ik) = density
         stream%overlap(:, :, ik) = overlap
      ELSE IF (stream%mode == cp2k_exchange_stream_reduced) THEN
         ALLOCATE (amat(stream%nao, stream%nao), cmat(stream%nao, stream%nao))
         DO spin = 1, stream%nspin
            amat = MATMUL(overlap, density(:, :, spin))
            cmat = 0.5_wp*MATMUL(amat, overlap)
            DO icell = 1, stream%nk
               angle = 2.0_wp*pi*DOT_PRODUCT(stream%kfrac(:, ik), &
                  & REAL(stream%reps(:, icell), wp))
               phase = stream%weights(ik)*EXP(CMPLX(0.0_wp, -angle, wp))
               stream%amat_r(:, :, icell, spin) = &
                  & stream%amat_r(:, :, icell, spin) + phase*amat
               stream%cmat_r(:, :, icell, spin) = &
                  & stream%cmat_r(:, :, icell, spin) + phase*cmat
               stream%vmat_r(:, :, icell, spin) = &
                  & stream%vmat_r(:, :, icell, spin) + phase*density(:, :, spin)
            END DO
         END DO
      ELSE IF (stream%transform_backend == bvk_transform_dense) THEN
         ALLOCATE (amat(stream%nao, stream%nao), cmat(stream%nao, stream%nao))
         DO spin = 1, stream%nspin
            amat = MATMUL(overlap, density(:, :, spin))
            cmat = 0.5_wp*MATMUL(amat, overlap)
            DO ibatch = 1, stream%batch_count
               icell = stream%batch_first + (ibatch - 1)
               inverse_cell = stream%inverse_image(icell)
               phase = stream%weights(ik)*CONJG( &
                  & regular_grid_phase_forward(stream%nmesh, &
                  & stream%input_to_grid, stream%reps(:, icell), &
                  & stream%phase_roots, stream%twist_phase(icell), ik))
               stream%amat_r(:, :, ibatch, spin) = &
                  & stream%amat_r(:, :, ibatch, spin) + phase*amat
               stream%cmat_r(:, :, ibatch, spin) = &
                  & stream%cmat_r(:, :, ibatch, spin) + phase*cmat
               stream%vmat_r(:, :, ibatch, spin) = &
                  & stream%vmat_r(:, :, ibatch, spin) &
                  & + phase*density(:, :, spin)
               phase = stream%weights(ik)*CONJG( &
                  & regular_grid_phase_forward(stream%nmesh, &
                  & stream%input_to_grid, stream%reps(:, inverse_cell), &
                  & stream%phase_roots, stream%twist_phase(inverse_cell), ik))
               stream%amat_partner_r(:, :, ibatch, spin) = &
                  & stream%amat_partner_r(:, :, ibatch, spin) + phase*amat
            END DO
         END DO
      ELSE
         ALLOCATE (amat(stream%nao, stream%nao), cmat(stream%nao, stream%nao))
         DO spin = 1, stream%nspin
            amat = MATMUL(overlap, density(:, :, spin))
            cmat = 0.5_wp*MATMUL(amat, overlap)
            DO ilocal = 1, stream%transform_pair_count
               ipair = stream%transform_pair_first + (ilocal - 1)
               iao = MODULO(ipair - 1, stream%nao) + 1
               jao = (ipair - 1)/stream%nao + 1
               stream%transform_k(ilocal, ik, spin, 1) = amat(iao, jao)
               stream%transform_k(ilocal, ik, spin, 2) = cmat(iao, jao)
               stream%transform_k(ilocal, ik, spin, 3) = &
                  & density(iao, jao, spin)
            END DO
         END DO
      END IF
      stream%seen(ik) = .TRUE.
   END SUBROUTINE cp2k_exchange_stream_push

!> Apply exchange after all blocks have been collected.
!>
!> Reduced mode applies the kernels directly to the accumulated BvK images;
!> oracle mode dispatches the unchanged complete-mesh evaluator.
   SUBROUTINE cp2k_exchange_stream_apply(stream, calc, mol, ecache, plan, vsh, energy, &
      & error)
      TYPE(cp2k_exchange_stream_type), INTENT(inout)     :: stream
      TYPE(xtb_calculator), INTENT(in)                   :: calc
      TYPE(structure_type), INTENT(in)                   :: mol
      TYPE(container_cache), INTENT(inout)               :: ecache
      TYPE(cp2k_stream_plan_type), INTENT(inout) :: plan
      REAL(wp), INTENT(out)                              :: vsh(:), energy
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      IF (.NOT. stream%active) THEN
         CALL fatal_error(error, "g-xTB exchange stream is not active")
         RETURN
      END IF
      IF (stream%applied) THEN
         CALL fatal_error(error, "g-xTB exchange stream was already applied")
         RETURN
      END IF
      IF (stream%mode == cp2k_exchange_stream_batched) THEN
         CALL fatal_error(error, &
            & "batched g-xTB exchange requires batch apply")
         RETURN
      END IF
      IF (.NOT. ALL(stream%seen)) THEN
         CALL fatal_error(error, "g-xTB exchange stream has missing blocks")
         RETURN
      END IF
      IF (SIZE(vsh) /= stream%nsh) THEN
         CALL fatal_error(error, &
            & "inconsistent g-xTB exchange stream shell-potential dimensions")
         RETURN
      END IF
      IF (.NOT. ALLOCATED(calc%exchange) .OR. .NOT. ALLOCATED(ecache%raw)) THEN
         CALL fatal_error(error, &
            & "g-xTB exchange stream state was invalidated before apply")
         RETURN
      END IF
      IF (.NOT. exchange_stream_state_matches(stream, calc, mol, ecache, plan)) THEN
         CALL fatal_error(error, &
            & "g-xTB exchange stream state was invalidated before apply")
         RETURN
      END IF

      SELECT TYPE (cache => ecache%raw)
      TYPE is (exchange_cache)
         SELECT TYPE (exchange => calc%exchange)
         TYPE is (exchange_fock)
            IF (stream%mode == cp2k_exchange_stream_oracle) THEN
               CALL consumer_get_kfock_kmesh(exchange, mol, cache, plan, plan%bvk_kernel, &
                  & stream%kfrac, stream%weights, stream%density, stream%overlap, &
                  & stream%fock, stream%vsh, stream%energy)
            ELSE
               CALL consumer_get_kfock_stream_apply(exchange, mol, cache, &
                  & plan%bvk_kernel, plan%bvk_phase_forward, stream%weights, &
                  & stream%amat_r, stream%cmat_r, stream%vmat_r, &
                  & stream%gdiagP, stream%gdiagSP, stream%gdiagSPS, &
                  & stream%vsh, stream%energy)
            END IF
            vsh = stream%vsh
            energy = stream%energy
            stream%applied = .TRUE.
         class default
            CALL fatal_error(error, &
               & "exchange stream requires the exchange_fock implementation")
         END SELECT
      class default
         CALL fatal_error(error, "unexpected g-xTB exchange cache type")
      END SELECT
   END SUBROUTINE cp2k_exchange_stream_apply

!> Retrieve one unweighted Fock block after a successful stream application.
!>
!> Reduced mode requires the corresponding overlap block again and
!> reconstructs only this requested k point.  Oracle mode reads the retained
!> complete-mesh result and does not require the optional overlap argument.
   SUBROUTINE cp2k_exchange_stream_get(stream, ik, fock, error, overlap)
      TYPE(cp2k_exchange_stream_type), INTENT(inout)     :: stream
      INTEGER, INTENT(in)                                :: ik
      COMPLEX(wp), INTENT(out)                           :: fock(:, :, :)
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error
      COMPLEX(wp), INTENT(in), OPTIONAL                  :: overlap(:, :)

      COMPLEX(wp), ALLOCATABLE                           :: phase_forward(:)
      INTEGER                                            :: icell
      REAL(wp)                                           :: angle

      IF (.NOT. stream%active) THEN
         CALL fatal_error(error, "g-xTB exchange stream is not active")
         RETURN
      END IF
      IF (.NOT. stream%applied) THEN
         CALL fatal_error(error, "g-xTB exchange stream was not applied")
         RETURN
      END IF
      IF (stream%mode == cp2k_exchange_stream_batched) THEN
         CALL fatal_error(error, &
            & "batched g-xTB exchange requires batch pulls")
         RETURN
      END IF
      IF (ik < 1 .OR. ik > stream%nk) THEN
         CALL fatal_error(error, "g-xTB exchange stream block index is out of range")
         RETURN
      END IF
      IF (stream%forward_pulled(ik)) THEN
         CALL fatal_error(error, "duplicate g-xTB exchange stream forward pull")
         RETURN
      END IF
      IF (SIZE(fock, 1) /= stream%nao .OR. SIZE(fock, 2) /= stream%nao &
         & .OR. SIZE(fock, 3) /= stream%nspin) THEN
         CALL fatal_error(error, &
            & "inconsistent g-xTB exchange stream Fock dimensions")
         RETURN
      END IF

      IF (stream%mode == cp2k_exchange_stream_oracle) THEN
         fock = stream%fock(:, :, :, ik)
         stream%forward_pulled(ik) = .TRUE.
         RETURN
      END IF
      IF (.NOT. PRESENT(overlap)) THEN
         CALL fatal_error(error, &
            & "memory-reduced g-xTB exchange stream get requires overlap block")
         RETURN
      END IF
      IF (SIZE(overlap, 1) /= stream%nao .OR. SIZE(overlap, 2) /= stream%nao) THEN
         CALL fatal_error(error, &
            & "inconsistent g-xTB exchange stream get overlap dimensions")
         RETURN
      END IF

      ALLOCATE (phase_forward(stream%nk))
      DO icell = 1, stream%nk
         angle = 2.0_wp*pi*DOT_PRODUCT(stream%kfrac(:, ik), &
            & REAL(stream%reps(:, icell), wp))
         phase_forward(icell) = EXP(CMPLX(0.0_wp, angle, wp))
      END DO
      CALL reconstruct_stream_fock_block(overlap, phase_forward, .TRUE., &
         & stream%amat_r, stream%cmat_r, stream%vmat_r, stream%gdiagP, &
         & stream%gdiagSP, stream%gdiagSPS, fock)
      stream%forward_pulled(ik) = .TRUE.
   END SUBROUTINE cp2k_exchange_stream_get

!> Return the inclusive representative range of the current image batch.
!>
!> Batched exchange is a deliberately sequential caller-driven protocol:
!> every Bloch block is pushed once for this range, the range is applied, and
!> every forward or reverse contribution is pulled once before the matching
!> advance routine may reuse the bounded image storage.
   SUBROUTINE cp2k_exchange_stream_batch_bounds(stream, first_image, last_image, &
      & error)
      TYPE(cp2k_exchange_stream_type), INTENT(in)        :: stream
      INTEGER, INTENT(out)                               :: first_image, last_image
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      first_image = 0
      last_image = 0
      IF (.NOT. stream%active) THEN
         CALL fatal_error(error, "g-xTB exchange stream is not active")
         RETURN
      END IF
      IF (stream%mode /= cp2k_exchange_stream_batched) THEN
         CALL fatal_error(error, "g-xTB exchange stream is not in batched mode")
         RETURN
      END IF
      IF ((.NOT. stream%reverse_transaction .AND. stream%applied) .OR. &
         & (stream%reverse_transaction .AND. stream%reverse_applied)) THEN
         CALL fatal_error(error, "g-xTB exchange image batches are complete")
         RETURN
      END IF
      first_image = stream%batch_first
      last_image = stream%batch_first + (stream%batch_count - 1)
   END SUBROUTINE cp2k_exchange_stream_batch_bounds

!> Transform one AO-pair tile from k space to the complete representative
!> order.  Only the tile is materialized; the result is immediately scattered
!> into the current bounded image batch by the caller.
   SUBROUTINE compact_stream_k_to_r(stream, spin, field, RESULT)
      TYPE(cp2k_exchange_stream_type), INTENT(in)        :: stream
      INTEGER, INTENT(in)                                :: spin, field
      COMPLEX(wp), INTENT(out)                           :: result(:, :)

      SELECT CASE (stream%transform_backend)
      CASE (bvk_transform_separable_dft)
         CALL bvk_separable_k_to_r( &
            & stream%transform_k(1:stream%transform_pair_count, :, spin, field), &
            & RESULT, stream%nmesh, stream%input_to_grid, stream%reps, &
            & stream%phase_roots, stream%twist_phase, .TRUE.)
      CASE (bvk_transform_mixed_radix_fft)
         CALL bvk_mixed_radix_k_to_r( &
            & stream%transform_k(1:stream%transform_pair_count, :, spin, field), &
            & RESULT, stream%nmesh, stream%input_to_grid, stream%reps, &
            & stream%phase_roots, stream%twist_phase, .TRUE.)
      END SELECT
   END SUBROUTINE compact_stream_k_to_r

!> Finalize one compact forward AO-pair tile.
!>
!> The caller resubmits the complete k mesh until `done` is true.  Each pass
!> stages only `transform_pair_batch_size*Nk` scalars per field; transformed
!> values are scattered directly into the bounded image batch, so no complete
!> Nk-by-nao**2 object persists.
   SUBROUTINE cp2k_exchange_stream_transform_advance(stream, done, error)
      TYPE(cp2k_exchange_stream_type), INTENT(inout)     :: stream
      LOGICAL, INTENT(out)                               :: done
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      COMPLEX(wp), ALLOCATABLE                           :: transformed(:, :)
      INTEGER                                            :: iao, ibatch, icell, ilocal, &
                                                            inverse_cell, ipair, jao, spin

      done = .FALSE.
      IF (.NOT. stream%active .OR. stream%mode /= cp2k_exchange_stream_batched) THEN
         CALL fatal_error(error, &
            & "g-xTB transform advance requires an active batched stream")
         RETURN
      END IF
      IF (stream%applied .OR. stream%batch_applied) THEN
         CALL fatal_error(error, &
            & "g-xTB transform advance followed exchange application")
         RETURN
      END IF
      IF (stream%partial_active .AND. .NOT. stream%partial_forward_imported &
         & .AND. .NOT. stream%partial_compact_image_storage) THEN
         IF (ANY(stream%seen .NEQV. stream%partial_forward_owned)) THEN
            CALL fatal_error(error, &
               & "compact partial forward transform tile changed k ownership")
            RETURN
         END IF
      ELSE IF (.NOT. ALL(stream%seen)) THEN
         CALL fatal_error(error, "g-xTB transform tile has missing k blocks")
         RETURN
      END IF
      IF (stream%transform_backend == bvk_transform_dense) THEN
         done = .TRUE.
         RETURN
      END IF
      IF (stream%reverse_transaction .OR. stream%transform_field_count /= 3) THEN
         CALL fatal_error(error, "invalid g-xTB forward transform state")
         RETURN
      END IF
      IF (stream%transform_pair_first < 1) THEN
         CALL fatal_error(error, "g-xTB forward transform batch is already complete")
         RETURN
      END IF

      ALLOCATE (transformed(stream%transform_pair_count, stream%nk))
      DO spin = 1, stream%nspin
         CALL compact_stream_k_to_r(stream, spin, 1, transformed)
         DO ilocal = 1, stream%transform_pair_count
            ipair = stream%transform_pair_first + (ilocal - 1)
            iao = MODULO(ipair - 1, stream%nao) + 1
            jao = (ipair - 1)/stream%nao + 1
            DO ibatch = 1, stream%batch_count
               icell = stream%batch_first + (ibatch - 1)
               inverse_cell = stream%inverse_image(icell)
               stream%amat_r(iao, jao, ibatch, spin) = &
                  & transformed(ilocal, icell)
               stream%amat_partner_r(iao, jao, ibatch, spin) = &
                  & transformed(ilocal, inverse_cell)
            END DO
         END DO
         CALL compact_stream_k_to_r(stream, spin, 2, transformed)
         DO ilocal = 1, stream%transform_pair_count
            ipair = stream%transform_pair_first + (ilocal - 1)
            iao = MODULO(ipair - 1, stream%nao) + 1
            jao = (ipair - 1)/stream%nao + 1
            DO ibatch = 1, stream%batch_count
               icell = stream%batch_first + (ibatch - 1)
               stream%cmat_r(iao, jao, ibatch, spin) = &
                  & transformed(ilocal, icell)
            END DO
         END DO
         CALL compact_stream_k_to_r(stream, spin, 3, transformed)
         DO ilocal = 1, stream%transform_pair_count
            ipair = stream%transform_pair_first + (ilocal - 1)
            iao = MODULO(ipair - 1, stream%nao) + 1
            jao = (ipair - 1)/stream%nao + 1
            DO ibatch = 1, stream%batch_count
               icell = stream%batch_first + (ibatch - 1)
               stream%vmat_r(iao, jao, ibatch, spin) = &
                  & transformed(ilocal, icell)
            END DO
         END DO
      END DO

      IF (stream%transform_pair_count > &
         & stream%nao*stream%nao - stream%transform_pair_first) THEN
         stream%transform_pair_first = 0
         stream%transform_pair_count = 0
         done = .TRUE.
      ELSE
         stream%transform_pair_first = stream%transform_pair_first &
            & + stream%transform_pair_count
         stream%transform_pair_count = MIN(stream%transform_pair_batch_size, &
            & stream%nao*stream%nao - stream%transform_pair_first + 1)
         stream%seen = .FALSE.
         stream%transform_k = (0.0_wp, 0.0_wp)
      END IF
   END SUBROUTINE cp2k_exchange_stream_transform_advance

!> Apply the current bounded image batch after all k blocks were resubmitted.
   SUBROUTINE cp2k_exchange_stream_batch_apply(stream, calc, mol, ecache, plan, error)
      TYPE(cp2k_exchange_stream_type), INTENT(inout)     :: stream
      TYPE(xtb_calculator), INTENT(in)                   :: calc
      TYPE(structure_type), INTENT(in)                   :: mol
      TYPE(container_cache), INTENT(inout)               :: ecache
      TYPE(cp2k_stream_plan_type), INTENT(inout) :: plan
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      INTEGER                                            :: first_global, first_local, ibatch, &
                                                            last_image, last_local, local_count
      INTEGER, ALLOCATABLE                               :: image_index(:)
      REAL(wp)                                           :: energy_batch
      REAL(wp), ALLOCATABLE                              :: vsh_batch(:)

      IF (.NOT. stream%active) THEN
         CALL fatal_error(error, "g-xTB exchange stream is not active")
         RETURN
      END IF
      IF (stream%mode /= cp2k_exchange_stream_batched) THEN
         CALL fatal_error(error, "g-xTB exchange stream is not in batched mode")
         RETURN
      END IF
      IF (stream%applied) THEN
         CALL fatal_error(error, "g-xTB exchange image batches are complete")
         RETURN
      END IF
      IF (stream%batch_applied) THEN
         CALL fatal_error(error, "g-xTB exchange image batch was already applied")
         RETURN
      END IF
      IF (stream%partial_active .AND. &
         & .NOT. stream%partial_forward_imported) THEN
         CALL fatal_error(error, &
            & "g-xTB exchange partial forward payload was not imported")
         RETURN
      END IF
      IF (.NOT. ALL(stream%seen)) THEN
         CALL fatal_error(error, "g-xTB exchange image batch has missing blocks")
         RETURN
      END IF
      IF (stream%transform_backend /= bvk_transform_dense .AND. &
         & stream%transform_pair_first /= 0) THEN
         CALL fatal_error(error, &
            & "g-xTB compact forward transform batch is incomplete")
         RETURN
      END IF
      IF (.NOT. exchange_stream_state_matches(stream, calc, mol, ecache, plan)) THEN
         CALL fatal_error(error, &
            & "g-xTB exchange image batch state was invalidated before apply")
         RETURN
      END IF

      ALLOCATE (vsh_batch(stream%nsh), source=0.0_wp)
      CALL partial_local_batch_bounds(stream, first_local, last_local, &
         & local_count, first_global)
      energy_batch = 0.0_wp
      IF (local_count > 0) THEN
         ALLOCATE (image_index(local_count))
         image_index = [(first_global + (ibatch - first_local), &
            & ibatch=first_local, last_local)]
         last_image = first_global + (local_count - 1)
         SELECT TYPE (cache => ecache%raw)
         TYPE is (exchange_cache)
            SELECT TYPE (exchange => calc%exchange)
            TYPE is (exchange_fock)
               CALL consumer_get_kfock_batch_apply(exchange, mol, cache, &
                  & plan%bvk_kernel, image_index, &
                  & stream%inverse_coeff(first_global:last_image), &
                  & stream%amat_r(:, :, first_local:last_local, :), &
                  & stream%amat_partner_r(:, :, first_local:last_local, :), &
                  & stream%cmat_r(:, :, first_local:last_local, :), &
                  & stream%vmat_r(:, :, first_local:last_local, :), &
                  & stream%gdiagP, stream%gdiagSP, stream%gdiagSPS, &
                  & vsh_batch, energy_batch)
            class default
               CALL fatal_error(error, &
                  & "exchange image batch requires the exchange_fock implementation")
               RETURN
            END SELECT
         class default
            CALL fatal_error(error, "unexpected g-xTB exchange cache type")
            RETURN
         END SELECT
      END IF
      stream%vsh = stream%vsh + vsh_batch
      stream%energy = stream%energy + energy_batch
      stream%batch_applied = .TRUE.
   END SUBROUTINE cp2k_exchange_stream_batch_apply

!> Add one k-point Fock contribution from the current image batch.
!>
!> The caller owns the final Fock mesh and must initialize each accumulator
!> before the first batch.  This routine adds exactly one batch contribution;
!> every k point must be pulled exactly once per applied batch.
   SUBROUTINE cp2k_exchange_stream_batch_pull(stream, ik, overlap, fock, error)
      TYPE(cp2k_exchange_stream_type), INTENT(inout)     :: stream
      INTEGER, INTENT(in)                                :: ik
      COMPLEX(wp), INTENT(in)                            :: overlap(:, :)
      COMPLEX(wp), INTENT(inout)                         :: fock(:, :, :)
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      COMPLEX(wp), ALLOCATABLE                           :: fock_batch(:, :, :), phase_forward(:)
      INTEGER                                            :: first_global, first_local, ibatch, &
                                                            icell, last_local, local_count
      LOGICAL                                            :: add_onsite

      IF (.NOT. stream%active) THEN
         CALL fatal_error(error, "g-xTB exchange stream is not active")
         RETURN
      END IF
      IF (stream%mode /= cp2k_exchange_stream_batched) THEN
         CALL fatal_error(error, "g-xTB exchange stream is not in batched mode")
         RETURN
      END IF
      IF (.NOT. stream%batch_applied) THEN
         CALL fatal_error(error, "g-xTB exchange image batch was not applied")
         RETURN
      END IF
      IF (ik < 1 .OR. ik > stream%nk) THEN
         CALL fatal_error(error, "g-xTB exchange stream block index is out of range")
         RETURN
      END IF
      IF (stream%partial_pull_distribution_set) THEN
         IF (.NOT. stream%partial_pull_mask(ik)) THEN
            CALL fatal_error(error, "unowned g-xTB exchange image-batch pull")
            RETURN
         END IF
      END IF
      IF (stream%batch_pulled(ik)) THEN
         CALL fatal_error(error, "duplicate g-xTB exchange image-batch pull")
         RETURN
      END IF
      IF (SIZE(overlap, 1) /= stream%nao &
         & .OR. SIZE(overlap, 2) /= stream%nao &
         & .OR. SIZE(fock, 1) /= stream%nao &
         & .OR. SIZE(fock, 2) /= stream%nao &
         & .OR. SIZE(fock, 3) /= stream%nspin) THEN
         CALL fatal_error(error, &
            & "inconsistent g-xTB exchange image-batch pull dimensions")
         RETURN
      END IF
      IF (stream%partial_active .AND. &
         & .NOT. stream%block_signature_set(ik)) THEN
         stream%overlap_signature(ik) = exchange_hash_complex2(overlap)
         stream%block_signature_set(ik) = .TRUE.
      ELSE IF (exchange_hash_complex2(overlap) /= &
         & stream%overlap_signature(ik)) THEN
         CALL fatal_error(error, &
            & "changed g-xTB exchange image-batch pull overlap")
         RETURN
      END IF
      CALL partial_local_batch_bounds(stream, first_local, last_local, &
         & local_count, first_global)
      ALLOCATE (fock_batch(stream%nao, stream%nao, stream%nspin), &
         & source=(0.0_wp, 0.0_wp))
      add_onsite = .FALSE.
      IF (local_count > 0) THEN
         ALLOCATE (phase_forward(local_count))
         DO ibatch = first_local, last_local
            icell = first_global + (ibatch - first_local)
            phase_forward(ibatch - first_local + 1) = &
               & regular_grid_phase_forward(stream%nmesh, &
               & stream%input_to_grid, stream%reps(:, icell), &
               & stream%phase_roots, stream%twist_phase(icell), ik)
            IF (ALL(stream%reps(:, icell) == 0)) add_onsite = .TRUE.
         END DO
         CALL reconstruct_stream_fock_block(overlap, phase_forward, add_onsite, &
            & stream%amat_r(:, :, first_local:last_local, :), &
            & stream%cmat_r(:, :, first_local:last_local, :), &
            & stream%vmat_r(:, :, first_local:last_local, :), &
            & stream%gdiagP, stream%gdiagSP, stream%gdiagSPS, fock_batch)
      END IF
      fock = fock + fock_batch
      stream%batch_pulled(ik) = .TRUE.
   END SUBROUTINE cp2k_exchange_stream_batch_pull

!> Reuse the bounded storage for the next sequential image batch.
   SUBROUTINE cp2k_exchange_stream_batch_advance(stream, done, error)
      TYPE(cp2k_exchange_stream_type), INTENT(inout)     :: stream
      LOGICAL, INTENT(out)                               :: done
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      done = .FALSE.
      IF (.NOT. stream%active) THEN
         CALL fatal_error(error, "g-xTB exchange stream is not active")
         RETURN
      END IF
      IF (stream%mode /= cp2k_exchange_stream_batched) THEN
         CALL fatal_error(error, "g-xTB exchange stream is not in batched mode")
         RETURN
      END IF
      IF (stream%applied) THEN
         CALL fatal_error(error, "g-xTB exchange image batches are complete")
         RETURN
      END IF
      IF (.NOT. stream%batch_applied) THEN
         CALL fatal_error(error, "g-xTB exchange image batch was not applied")
         RETURN
      END IF
      IF (stream%partial_pull_distribution_set) THEN
         IF (ANY(stream%batch_pulled .NEQV. stream%partial_pull_mask)) THEN
            CALL fatal_error(error, &
               & "g-xTB exchange image batch violates pull ownership")
            RETURN
         END IF
      ELSE
         IF (.NOT. ALL(stream%batch_pulled)) THEN
            CALL fatal_error(error, "g-xTB exchange image batch has missing pulls")
            RETURN
         END IF
      END IF

      IF (stream%batch_count > stream%batch_last - stream%batch_first) THEN
         stream%applied = .TRUE.
         IF (stream%partial_pull_distribution_set &
            & .AND. .NOT. stream%partial_result_owner) &
            & stream%batch_result_pulled = .TRUE.
         done = .TRUE.
         RETURN
      END IF

      stream%batch_first = stream%batch_first + stream%batch_count
      stream%batch_count = MIN(stream%batch_size, &
         & stream%batch_last - stream%batch_first + 1)
      stream%seen = .FALSE.
      stream%batch_pulled = .FALSE.
      stream%batch_applied = .FALSE.
      IF (stream%transform_backend /= bvk_transform_dense) THEN
         stream%transform_pair_first = 1
         stream%transform_pair_count = MIN( &
            & stream%transform_pair_batch_size, stream%nao*stream%nao)
         stream%transform_k = (0.0_wp, 0.0_wp)
      END IF
      IF (stream%partial_active) THEN
         stream%partial_generation = stream%partial_generation + 1
         stream%partial_forward_exported = .FALSE.
         stream%partial_forward_imported = .FALSE.
      END IF
      stream%amat_r = (0.0_wp, 0.0_wp)
      stream%amat_partner_r = (0.0_wp, 0.0_wp)
      stream%cmat_r = (0.0_wp, 0.0_wp)
      stream%vmat_r = (0.0_wp, 0.0_wp)
   END SUBROUTINE cp2k_exchange_stream_batch_advance

!> Return the accumulated energy and shell potential exactly once after all
!> batches.  A successful result pull is mandatory before the stream can end.
   SUBROUTINE cp2k_exchange_stream_batch_result(stream, calc, mol, ecache, plan, vsh, &
      & energy, error)
      TYPE(cp2k_exchange_stream_type), INTENT(inout)     :: stream
      TYPE(xtb_calculator), INTENT(in)                   :: calc
      TYPE(structure_type), INTENT(in)                   :: mol
      TYPE(container_cache), INTENT(in)                  :: ecache
      TYPE(cp2k_stream_plan_type), INTENT(in) :: plan
      REAL(wp), INTENT(out)                              :: vsh(:), energy
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      IF (.NOT. stream%active) THEN
         CALL fatal_error(error, "g-xTB exchange stream is not active")
         RETURN
      END IF
      IF (stream%mode /= cp2k_exchange_stream_batched) THEN
         CALL fatal_error(error, "g-xTB exchange stream is not in batched mode")
         RETURN
      END IF
      IF (.NOT. stream%applied) THEN
         CALL fatal_error(error, "g-xTB exchange image batches are incomplete")
         RETURN
      END IF
      IF (stream%batch_result_pulled) THEN
         CALL fatal_error(error, &
            & "g-xTB exchange image-batch result was already pulled")
         RETURN
      END IF
      IF (SIZE(vsh) /= stream%nsh) THEN
         CALL fatal_error(error, &
            & "inconsistent g-xTB exchange stream shell-potential dimensions")
         RETURN
      END IF
      IF (.NOT. exchange_stream_state_matches(stream, calc, mol, ecache, plan)) THEN
         CALL fatal_error(error, &
            & "g-xTB exchange image batch state was invalidated before result")
         RETURN
      END IF
      vsh = stream%vsh
      energy = stream%energy
      stream%batch_result_pulled = .TRUE.
   END SUBROUTINE cp2k_exchange_stream_batch_result

!> Resubmit one unweighted density/overlap block to the current reverse batch.
   SUBROUTINE cp2k_exchange_stream_reverse_push(stream, ik, density, overlap, &
      & error)
      TYPE(cp2k_exchange_stream_type), INTENT(inout)     :: stream
      INTEGER, INTENT(in)                                :: ik
      COMPLEX(wp), INTENT(in)                            :: density(:, :, :), overlap(:, :)
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      COMPLEX(wp)                                        :: phase_adjoint, phase_source
      COMPLEX(wp), ALLOCATABLE                           :: amat(:, :), bmat(:, :), dmat(:, :), &
                                                            vseed(:, :), xmat(:, :)
      INTEGER                                            :: iao, ibatch, icell, ilocal, ipair, jao, &
                                                            spin
      INTEGER(int64)                                     :: density_signature, overlap_signature
      LOGICAL                                            :: has_origin
      REAL(wp)                                           :: angle, spin_factor

      IF (.NOT. stream%active .OR. .NOT. stream%reverse_transaction) THEN
         CALL fatal_error(error, &
            & "g-xTB exchange stream reverse transaction is not active")
         RETURN
      END IF
      IF (stream%partial_active .AND. stream%partial_reverse_imported) THEN
         CALL fatal_error(error, &
            & "g-xTB exchange reverse partial push followed import")
         RETURN
      END IF
      IF (stream%partial_active .AND. stream%partial_reverse_exported) THEN
         CALL fatal_error(error, &
            & "g-xTB exchange reverse partial push followed export")
         RETURN
      END IF
      IF (stream%reverse_applied) THEN
         CALL fatal_error(error, "g-xTB exchange stream reverse is complete")
         RETURN
      END IF
      IF (stream%reverse_batch_applied) THEN
         CALL fatal_error(error, &
            & "g-xTB exchange stream reverse batch was already applied")
         RETURN
      END IF
      IF (stream%transform_backend /= bvk_transform_dense .AND. &
         & stream%transform_pair_first < 1) THEN
         CALL fatal_error(error, &
            & "g-xTB compact reverse transform batch is already complete")
         RETURN
      END IF
      IF (ik < 1 .OR. ik > stream%nk) THEN
         CALL fatal_error(error, "g-xTB exchange stream block index is out of range")
         RETURN
      END IF
      IF (stream%seen(ik)) THEN
         CALL fatal_error(error, "duplicate g-xTB exchange reverse block")
         RETURN
      END IF
      IF (SIZE(density, 1) /= stream%nao &
         & .OR. SIZE(density, 2) /= stream%nao &
         & .OR. SIZE(density, 3) /= stream%nspin &
         & .OR. SIZE(overlap, 1) /= stream%nao &
         & .OR. SIZE(overlap, 2) /= stream%nao) THEN
         CALL fatal_error(error, &
            & "inconsistent g-xTB exchange reverse block dimensions")
         RETURN
      END IF

      density_signature = exchange_hash_complex3(density)
      overlap_signature = exchange_hash_complex2(overlap)
      IF (stream%block_signature_set(ik)) THEN
         IF (density_signature /= stream%density_signature(ik) &
            & .OR. overlap_signature /= stream%overlap_signature(ik)) THEN
            CALL fatal_error(error, "changed g-xTB exchange reverse block")
            RETURN
         END IF
      ELSE
         stream%density_signature(ik) = density_signature
         stream%overlap_signature(ik) = overlap_signature
         stream%block_signature_set(ik) = .TRUE.
      END IF

      ALLOCATE (amat(stream%nao, stream%nao), bmat(stream%nao, stream%nao), &
         & dmat(stream%nao, stream%nao), vseed(stream%nao, stream%nao), &
         & xmat(stream%nao, stream%nao))
      spin_factor = 0.5_wp
      IF (stream%nspin > 1) spin_factor = 1.0_wp
      has_origin = .FALSE.
      DO ibatch = 1, stream%batch_count
         icell = stream%batch_first + (ibatch - 1)
         IF (ALL(stream%reps(:, icell) == 0)) has_origin = .TRUE.
      END DO

      DO spin = 1, stream%nspin
         xmat = MATMUL(overlap, density(:, :, spin))
         dmat = 0.5_wp*MATMUL(xmat, overlap)
         IF (stream%transform_backend == bvk_transform_dense) THEN
            amat = -0.25_wp*spin_factor*stream%weights(ik) &
                  & *density(:, :, spin)
         ELSE
            amat = -0.25_wp*spin_factor*density(:, :, spin)
         END IF
         bmat = MATMUL(amat, CONJG(TRANSPOSE(overlap)))
         vseed = 0.5_wp*MATMUL(CONJG(TRANSPOSE(overlap)), bmat)

         IF (has_origin .AND. (stream%transform_backend == bvk_transform_dense &
            & .OR. stream%transform_pair_first == 1)) THEN
            DO iao = 1, stream%nao
               stream%gdiagP(iao, spin) = stream%gdiagP(iao, spin) &
                  & + stream%weights(ik)*density(iao, iao, spin)
               stream%gdiagSP(iao, spin) = stream%gdiagSP(iao, spin) &
                  & + stream%weights(ik)*xmat(iao, iao)
               stream%gdiagSPS(iao, spin) = stream%gdiagSPS(iao, spin) &
                  & + 2.0_wp*stream%weights(ik)*dmat(iao, iao)
               stream%reverse_bdiagSP(iao, spin) = &
                  & stream%reverse_bdiagSP(iao, spin) &
                  & + 0.5_wp*DOT_PRODUCT(overlap(:, iao), amat(:, iao)) &
                  & *MERGE(1.0_wp, stream%weights(ik), &
                  & stream%transform_backend == bvk_transform_dense)
               stream%reverse_bdiagSPS(iao, spin) = &
                  & stream%reverse_bdiagSPS(iao, spin) + 0.25_wp*amat(iao, iao) &
                  & *MERGE(1.0_wp, stream%weights(ik), &
                  & stream%transform_backend == bvk_transform_dense)
            END DO
         END IF

         IF (stream%transform_backend == bvk_transform_dense) THEN
            DO ibatch = 1, stream%batch_count
               icell = stream%batch_first + (ibatch - 1)
               angle = 2.0_wp*pi*DOT_PRODUCT(stream%kfrac(:, ik), &
                  & REAL(stream%reps(:, icell), wp))
               phase_adjoint = EXP(CMPLX(0.0_wp, -angle, wp))
               phase_source = stream%weights(ik)*phase_adjoint
               stream%amat_r(:, :, ibatch, spin) = &
                  & stream%amat_r(:, :, ibatch, spin) + phase_source*xmat
               stream%cmat_r(:, :, ibatch, spin) = &
                  & stream%cmat_r(:, :, ibatch, spin) + phase_source*dmat
               stream%vmat_r(:, :, ibatch, spin) = &
                  & stream%vmat_r(:, :, ibatch, spin) &
                  & + phase_source*density(:, :, spin)
               stream%amat_partner_r(:, :, ibatch, spin) = &
                  & stream%amat_partner_r(:, :, ibatch, spin) + phase_adjoint*bmat
               stream%reverse_amat_r(:, :, ibatch, spin) = &
                  & stream%reverse_amat_r(:, :, ibatch, spin) + phase_adjoint*amat
               stream%reverse_vmat_r(:, :, ibatch, spin) = &
                  & stream%reverse_vmat_r(:, :, ibatch, spin) + phase_adjoint*vseed
            END DO
         ELSE
            DO ilocal = 1, stream%transform_pair_count
               ipair = stream%transform_pair_first + (ilocal - 1)
               iao = MODULO(ipair - 1, stream%nao) + 1
               jao = (ipair - 1)/stream%nao + 1
               stream%transform_k(ilocal, ik, spin, 1) = xmat(iao, jao)
               stream%transform_k(ilocal, ik, spin, 2) = dmat(iao, jao)
               stream%transform_k(ilocal, ik, spin, 3) = &
                  & density(iao, jao, spin)
               stream%transform_k(ilocal, ik, spin, 4) = bmat(iao, jao)
               stream%transform_k(ilocal, ik, spin, 5) = amat(iao, jao)
               stream%transform_k(ilocal, ik, spin, 6) = vseed(iao, jao)
            END DO
         END IF
      END DO
      stream%seen(ik) = .TRUE.
   END SUBROUTINE cp2k_exchange_stream_reverse_push

!> Finalize one compact reverse AO-pair tile.
!>
!> All six reverse source maps use the same normalized regular-grid transform.
!> Fields that carried a k weight in the dense character sum are staged here
!> without that weight, so the transform normalization supplies it exactly
!> once.  The subsequent image-local reverse kernel and R-to-k pull are
!> therefore the real adjoint of the compact forward path.
   SUBROUTINE cp2k_exchange_stream_reverse_transform_advance(stream, done, error)
      TYPE(cp2k_exchange_stream_type), INTENT(inout)     :: stream
      LOGICAL, INTENT(out)                               :: done
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      COMPLEX(wp), ALLOCATABLE                           :: transformed(:, :)
      INTEGER                                            :: iao, ibatch, icell, ilocal, ipair, jao, &
                                                            spin

      done = .FALSE.
      IF (.NOT. stream%active .OR. .NOT. stream%reverse_transaction &
         & .OR. stream%mode /= cp2k_exchange_stream_batched) THEN
         CALL fatal_error(error, &
            & "g-xTB reverse transform advance requires a batched transaction")
         RETURN
      END IF
      IF (stream%reverse_applied .OR. stream%reverse_batch_applied) THEN
         CALL fatal_error(error, &
            & "g-xTB reverse transform advance followed kernel application")
         RETURN
      END IF
      IF (stream%partial_active .AND. .NOT. stream%partial_reverse_imported &
         & .AND. .NOT. stream%partial_compact_image_storage) THEN
         IF (ANY(stream%seen .NEQV. stream%partial_forward_owned)) THEN
            CALL fatal_error(error, &
               & "compact partial reverse transform tile changed k ownership")
            RETURN
         END IF
      ELSE IF (.NOT. ALL(stream%seen)) THEN
         CALL fatal_error(error, "g-xTB reverse transform tile has missing k blocks")
         RETURN
      END IF
      IF (stream%transform_backend == bvk_transform_dense) THEN
         done = .TRUE.
         RETURN
      END IF
      IF (stream%transform_field_count /= 6) THEN
         CALL fatal_error(error, "invalid g-xTB reverse transform state")
         RETURN
      END IF
      IF (stream%transform_pair_first < 1) THEN
         CALL fatal_error(error, "g-xTB reverse transform batch is already complete")
         RETURN
      END IF

      ALLOCATE (transformed(stream%transform_pair_count, stream%nk))
      DO spin = 1, stream%nspin
         CALL compact_stream_k_to_r(stream, spin, 1, transformed)
         DO ilocal = 1, stream%transform_pair_count
            ipair = stream%transform_pair_first + (ilocal - 1)
            iao = MODULO(ipair - 1, stream%nao) + 1
            jao = (ipair - 1)/stream%nao + 1
            DO ibatch = 1, stream%batch_count
               icell = stream%batch_first + (ibatch - 1)
               stream%amat_r(iao, jao, ibatch, spin) = &
                  & transformed(ilocal, icell)
            END DO
         END DO
         CALL compact_stream_k_to_r(stream, spin, 2, transformed)
         DO ilocal = 1, stream%transform_pair_count
            ipair = stream%transform_pair_first + (ilocal - 1)
            iao = MODULO(ipair - 1, stream%nao) + 1
            jao = (ipair - 1)/stream%nao + 1
            DO ibatch = 1, stream%batch_count
               icell = stream%batch_first + (ibatch - 1)
               stream%cmat_r(iao, jao, ibatch, spin) = &
                  & transformed(ilocal, icell)
            END DO
         END DO
         CALL compact_stream_k_to_r(stream, spin, 3, transformed)
         DO ilocal = 1, stream%transform_pair_count
            ipair = stream%transform_pair_first + (ilocal - 1)
            iao = MODULO(ipair - 1, stream%nao) + 1
            jao = (ipair - 1)/stream%nao + 1
            DO ibatch = 1, stream%batch_count
               icell = stream%batch_first + (ibatch - 1)
               stream%vmat_r(iao, jao, ibatch, spin) = &
                  & transformed(ilocal, icell)
            END DO
         END DO
         CALL compact_stream_k_to_r(stream, spin, 4, transformed)
         DO ilocal = 1, stream%transform_pair_count
            ipair = stream%transform_pair_first + (ilocal - 1)
            iao = MODULO(ipair - 1, stream%nao) + 1
            jao = (ipair - 1)/stream%nao + 1
            DO ibatch = 1, stream%batch_count
               icell = stream%batch_first + (ibatch - 1)
               stream%amat_partner_r(iao, jao, ibatch, spin) = &
                  & transformed(ilocal, icell)
            END DO
         END DO
         CALL compact_stream_k_to_r(stream, spin, 5, transformed)
         DO ilocal = 1, stream%transform_pair_count
            ipair = stream%transform_pair_first + (ilocal - 1)
            iao = MODULO(ipair - 1, stream%nao) + 1
            jao = (ipair - 1)/stream%nao + 1
            DO ibatch = 1, stream%batch_count
               icell = stream%batch_first + (ibatch - 1)
               stream%reverse_amat_r(iao, jao, ibatch, spin) = &
                  & transformed(ilocal, icell)
            END DO
         END DO
         CALL compact_stream_k_to_r(stream, spin, 6, transformed)
         DO ilocal = 1, stream%transform_pair_count
            ipair = stream%transform_pair_first + (ilocal - 1)
            iao = MODULO(ipair - 1, stream%nao) + 1
            jao = (ipair - 1)/stream%nao + 1
            DO ibatch = 1, stream%batch_count
               icell = stream%batch_first + (ibatch - 1)
               stream%reverse_vmat_r(iao, jao, ibatch, spin) = &
                  & transformed(ilocal, icell)
            END DO
         END DO
      END DO

      IF (stream%transform_pair_count > &
         & stream%nao*stream%nao - stream%transform_pair_first) THEN
         stream%transform_pair_first = 0
         stream%transform_pair_count = 0
         done = .TRUE.
      ELSE
         stream%transform_pair_first = stream%transform_pair_first &
            & + stream%transform_pair_count
         stream%transform_pair_count = MIN(stream%transform_pair_batch_size, &
            & stream%nao*stream%nao - stream%transform_pair_first + 1)
         stream%seen = .FALSE.
         stream%transform_k = (0.0_wp, 0.0_wp)
      END IF
   END SUBROUTINE cp2k_exchange_stream_reverse_transform_advance

!> Apply the current bounded reverse batch after all k blocks were resubmitted.
   SUBROUTINE cp2k_exchange_stream_reverse_batch_apply(stream, calc, mol, ecache, plan, &
      & error)
      TYPE(cp2k_exchange_stream_type), INTENT(inout)     :: stream
      TYPE(xtb_calculator), INTENT(in)                   :: calc
      TYPE(structure_type), INTENT(in)                   :: mol
      TYPE(container_cache), INTENT(inout)               :: ecache
      TYPE(cp2k_stream_plan_type), INTENT(inout) :: plan
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      INTEGER                                            :: first_global, first_local, ibatch, &
                                                            last_local, local_count
      INTEGER, ALLOCATABLE                               :: image_index(:)

      IF (.NOT. stream%active .OR. .NOT. stream%reverse_transaction) THEN
         CALL fatal_error(error, &
            & "g-xTB exchange stream reverse transaction is not active")
         RETURN
      END IF
      IF (stream%reverse_applied) THEN
         CALL fatal_error(error, "g-xTB exchange stream reverse is complete")
         RETURN
      END IF
      IF (stream%reverse_batch_applied) THEN
         CALL fatal_error(error, &
            & "g-xTB exchange stream reverse batch was already applied")
         RETURN
      END IF
      IF (stream%partial_active .AND. &
         & .NOT. stream%partial_reverse_imported) THEN
         CALL fatal_error(error, &
            & "g-xTB exchange partial reverse payload was not imported")
         RETURN
      END IF
      IF (.NOT. ALL(stream%seen)) THEN
         CALL fatal_error(error, "g-xTB exchange reverse batch has missing blocks")
         RETURN
      END IF
      IF (stream%transform_backend /= bvk_transform_dense .AND. &
         & stream%transform_pair_first /= 0) THEN
         CALL fatal_error(error, &
            & "g-xTB compact reverse transform batch is incomplete")
         RETURN
      END IF
      IF (.NOT. exchange_stream_state_matches(stream, calc, mol, ecache, plan)) THEN
         CALL fatal_error(error, &
            & "g-xTB exchange reverse state was invalidated before apply")
         RETURN
      END IF

      CALL partial_local_batch_bounds(stream, first_local, last_local, &
         & local_count, first_global)
      IF (local_count > 0) THEN
         ALLOCATE (image_index(local_count))
         image_index = [(first_global + (ibatch - first_local), &
            & ibatch=first_local, last_local)]
         SELECT TYPE (cache => ecache%raw)
         TYPE is (exchange_cache)
            SELECT TYPE (exchange => calc%exchange)
            TYPE is (exchange_fock)
               CALL consumer_get_kgrad_batch_apply(exchange, mol, cache, plan%bvk_kernel, &
                  & image_index, stream%amat_r(:, :, first_local:last_local, :), &
                  & stream%cmat_r(:, :, first_local:last_local, :), &
                  & stream%vmat_r(:, :, first_local:last_local, :), &
                  & stream%amat_partner_r(:, :, first_local:last_local, :), &
                  & stream%reverse_amat_r(:, :, first_local:last_local, :), &
                  & stream%reverse_vmat_r(:, :, first_local:last_local, :), &
                  & stream%gdiagP, stream%gdiagSP, stream%gdiagSPS, &
                  & stream%reverse_bdiagSP, stream%reverse_bdiagSPS, &
                  & stream%reverse_gradient, stream%reverse_sigma)
            class default
               CALL fatal_error(error, &
                  & "exchange reverse batch requires the exchange_fock implementation")
               RETURN
            END SELECT
         class default
            CALL fatal_error(error, "unexpected g-xTB exchange cache type")
            RETURN
         END SELECT
      END IF
      stream%reverse_pulled = .FALSE.
      stream%reverse_batch_applied = .TRUE.
   END SUBROUTINE cp2k_exchange_stream_reverse_batch_apply

!> Add one bounded-batch contribution to an unweighted overlap adjoint.
   SUBROUTINE cp2k_exchange_stream_reverse_batch_pull(stream, ik, density, &
      & overlap, overlap_adjoint, error)
      TYPE(cp2k_exchange_stream_type), INTENT(inout)     :: stream
      INTEGER, INTENT(in)                                :: ik
      COMPLEX(wp), INTENT(in)                            :: density(:, :, :), overlap(:, :)
      COMPLEX(wp), INTENT(inout)                         :: overlap_adjoint(:, :)
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      COMPLEX(wp)                                        :: phase
      COMPLEX(wp), ALLOCATABLE :: amat(:, :), bmat(:, :), dadj(:, :), raw(:, :), raw_spin(:, :), &
         tmat(:, :), vmat(:, :), work(:, :), xadj(:, :), xmat(:, :)
      INTEGER                                            :: first_global, first_local, iao, ibatch, &
                                                            icell, last_local, local_count, spin
      LOGICAL                                            :: has_origin
      REAL(wp)                                           :: angle, spin_factor

      IF (.NOT. stream%active .OR. .NOT. stream%reverse_transaction) THEN
         CALL fatal_error(error, &
            & "g-xTB exchange stream reverse transaction is not active")
         RETURN
      END IF
      IF (.NOT. stream%reverse_batch_applied) THEN
         CALL fatal_error(error, "g-xTB exchange reverse batch was not applied")
         RETURN
      END IF
      IF (ik < 1 .OR. ik > stream%nk) THEN
         CALL fatal_error(error, "g-xTB exchange stream block index is out of range")
         RETURN
      END IF
      IF (stream%partial_pull_distribution_set) THEN
         IF (.NOT. stream%partial_pull_mask(ik)) THEN
            CALL fatal_error(error, "unowned g-xTB exchange reverse-batch pull")
            RETURN
         END IF
      END IF
      IF (stream%reverse_pulled(ik)) THEN
         CALL fatal_error(error, "duplicate g-xTB exchange reverse-batch pull")
         RETURN
      END IF
      IF (SIZE(density, 1) /= stream%nao &
         & .OR. SIZE(density, 2) /= stream%nao &
         & .OR. SIZE(density, 3) /= stream%nspin &
         & .OR. SIZE(overlap, 1) /= stream%nao &
         & .OR. SIZE(overlap, 2) /= stream%nao &
         & .OR. SIZE(overlap_adjoint, 1) /= stream%nao &
         & .OR. SIZE(overlap_adjoint, 2) /= stream%nao) THEN
         CALL fatal_error(error, &
            & "inconsistent g-xTB exchange reverse pull dimensions")
         RETURN
      END IF
      IF (stream%partial_active .AND. &
         & .NOT. stream%block_signature_set(ik)) THEN
         stream%density_signature(ik) = exchange_hash_complex3(density)
         stream%overlap_signature(ik) = exchange_hash_complex2(overlap)
         stream%block_signature_set(ik) = .TRUE.
      ELSE IF (exchange_hash_complex3(density) /= &
         & stream%density_signature(ik) &
         & .OR. exchange_hash_complex2(overlap) /= &
         & stream%overlap_signature(ik)) THEN
         CALL fatal_error(error, "changed g-xTB exchange reverse pull block")
         RETURN
      END IF

      CALL partial_local_batch_bounds(stream, first_local, last_local, &
         & local_count, first_global)
      IF (local_count == 0) THEN
         stream%reverse_pulled(ik) = .TRUE.
         RETURN
      END IF

      ALLOCATE (amat(stream%nao, stream%nao), bmat(stream%nao, stream%nao), &
         & dadj(stream%nao, stream%nao), raw(stream%nao, stream%nao), &
         & raw_spin(stream%nao, stream%nao), tmat(stream%nao, stream%nao), &
         & vmat(stream%nao, stream%nao), &
         & work(stream%nao, stream%nao), xadj(stream%nao, stream%nao), &
         & xmat(stream%nao, stream%nao), source=(0.0_wp, 0.0_wp))
      spin_factor = 0.5_wp
      IF (stream%nspin > 1) spin_factor = 1.0_wp
      has_origin = .FALSE.
      DO ibatch = first_local, last_local
         icell = first_global + (ibatch - first_local)
         IF (ALL(stream%reps(:, icell) == 0)) has_origin = .TRUE.
      END DO

      raw = (0.0_wp, 0.0_wp)
      DO spin = 1, stream%nspin
         tmat = (0.0_wp, 0.0_wp)
         vmat = (0.0_wp, 0.0_wp)
         dadj = (0.0_wp, 0.0_wp)
         xadj = (0.0_wp, 0.0_wp)
         DO ibatch = first_local, last_local
            icell = first_global + (ibatch - first_local)
            angle = 2.0_wp*pi*DOT_PRODUCT(stream%kfrac(:, ik), &
               & REAL(stream%reps(:, icell), wp))
            phase = EXP(CMPLX(0.0_wp, angle, wp))
            tmat = tmat + phase*stream%amat_r(:, :, ibatch, spin)
            vmat = vmat + phase*stream%vmat_r(:, :, ibatch, spin)
            dadj = dadj + stream%weights(ik)*phase &
                  & *stream%cmat_r(:, :, ibatch, spin)
            xadj = xadj + stream%weights(ik)*phase &
                  & *stream%amat_partner_r(:, :, ibatch, spin)
         END DO

         xmat = MATMUL(overlap, density(:, :, spin))
         amat = -0.25_wp*spin_factor*stream%weights(ik) &
               & *density(:, :, spin)
         bmat = MATMUL(amat, CONJG(TRANSPOSE(overlap)))
         work = tmat + 0.5_wp*MATMUL(overlap, vmat)
         raw_spin = (0.0_wp, 0.0_wp)

         IF (has_origin) THEN
            DO iao = 1, stream%nao
               work(:, iao) = work(:, iao) &
                  & + 0.25_wp*stream%gdiagP(iao, spin)*overlap(:, iao)
            END DO
         END IF
         raw_spin = raw_spin + MATMUL(CONJG(TRANSPOSE(work)), amat)
         IF (has_origin) THEN
            DO iao = 1, stream%nao
               raw_spin(:, iao) = raw_spin(:, iao) &
                  & + 0.5_wp*CONJG(stream%gdiagSP(iao, spin))*amat(:, iao) &
                  & + 0.25_wp*CONJG(stream%gdiagP(iao, spin))*bmat(:, iao)
            END DO
         END IF
         raw_spin = raw_spin + 0.5_wp*MATMUL(bmat, CONJG(TRANSPOSE(vmat)))

         IF (has_origin) THEN
            DO iao = 1, stream%nao
               xadj(iao, iao) = xadj(iao, iao) &
                  & + stream%weights(ik)*stream%reverse_bdiagSP(iao, spin)
               xadj(iao, :) = xadj(iao, :) &
                  & + stream%weights(ik)*stream%reverse_bdiagSPS(iao, spin) &
                  & *overlap(iao, :)
               raw_spin(:, iao) = raw_spin(:, iao) &
                  & + stream%weights(ik)*stream%reverse_bdiagSPS(iao, spin) &
                  & *CONJG(xmat(iao, :))
            END DO
         END IF
         xadj = xadj + 0.5_wp*MATMUL(dadj, CONJG(TRANSPOSE(overlap)))
         raw_spin = raw_spin &
            & + 0.5_wp*MATMUL(CONJG(TRANSPOSE(xmat)), dadj) &
            & + MATMUL(xadj, CONJG(TRANSPOSE(density(:, :, spin))))
         raw = raw + raw_spin
      END DO

      overlap_adjoint = overlap_adjoint &
         & + 0.5_wp*(raw + CONJG(TRANSPOSE(raw)))/stream%weights(ik)
      stream%reverse_pulled(ik) = .TRUE.
   END SUBROUTINE cp2k_exchange_stream_reverse_batch_pull

!> Reuse bounded reverse storage for the next image batch.
   SUBROUTINE cp2k_exchange_stream_reverse_batch_advance(stream, done, error)
      TYPE(cp2k_exchange_stream_type), INTENT(inout)     :: stream
      LOGICAL, INTENT(out)                               :: done
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      done = .FALSE.
      IF (.NOT. stream%active .OR. .NOT. stream%reverse_transaction) THEN
         CALL fatal_error(error, &
            & "g-xTB exchange stream reverse transaction is not active")
         RETURN
      END IF
      IF (.NOT. stream%reverse_batch_applied) THEN
         CALL fatal_error(error, "g-xTB exchange reverse batch was not applied")
         RETURN
      END IF
      IF (stream%partial_pull_distribution_set) THEN
         IF (ANY(stream%reverse_pulled .NEQV. stream%partial_pull_mask)) THEN
            CALL fatal_error(error, &
               & "g-xTB exchange reverse batch violates pull ownership")
            RETURN
         END IF
      ELSE
         IF (.NOT. ALL(stream%reverse_pulled)) THEN
            CALL fatal_error(error, &
               & "g-xTB exchange reverse batch has missing pulls")
            RETURN
         END IF
      END IF
      IF (stream%batch_count > stream%batch_last - stream%batch_first) THEN
         stream%reverse_applied = .TRUE.
         IF (stream%partial_pull_distribution_set &
            & .AND. .NOT. stream%partial_result_owner) &
            & stream%reverse_result_pulled = .TRUE.
         done = .TRUE.
         RETURN
      END IF

      stream%batch_first = stream%batch_first + stream%batch_count
      stream%batch_count = MIN(stream%batch_size, &
         & stream%batch_last - stream%batch_first + 1)
      stream%seen = .FALSE.
      stream%reverse_pulled = .FALSE.
      stream%reverse_batch_applied = .FALSE.
      IF (stream%transform_backend /= bvk_transform_dense) THEN
         stream%transform_pair_first = 1
         stream%transform_pair_count = MIN( &
            & stream%transform_pair_batch_size, stream%nao*stream%nao)
         stream%transform_k = (0.0_wp, 0.0_wp)
      END IF
      IF (stream%partial_active) THEN
         stream%partial_generation = stream%partial_generation + 1
         stream%partial_reverse_exported = .FALSE.
         stream%partial_reverse_imported = .FALSE.
      END IF
      stream%amat_r = (0.0_wp, 0.0_wp)
      stream%cmat_r = (0.0_wp, 0.0_wp)
      stream%vmat_r = (0.0_wp, 0.0_wp)
      stream%amat_partner_r = (0.0_wp, 0.0_wp)
      stream%reverse_amat_r = (0.0_wp, 0.0_wp)
      stream%reverse_vmat_r = (0.0_wp, 0.0_wp)
      stream%gdiagP = (0.0_wp, 0.0_wp)
      stream%gdiagSP = (0.0_wp, 0.0_wp)
      stream%gdiagSPS = (0.0_wp, 0.0_wp)
      stream%reverse_bdiagSP = (0.0_wp, 0.0_wp)
      stream%reverse_bdiagSPS = (0.0_wp, 0.0_wp)
   END SUBROUTINE cp2k_exchange_stream_reverse_batch_advance

!> Return accumulated force and positive-strain responses exactly once.
   SUBROUTINE cp2k_exchange_stream_reverse_result(stream, calc, mol, ecache, plan, &
      & gradient, sigma, error)
      TYPE(cp2k_exchange_stream_type), INTENT(inout)     :: stream
      TYPE(xtb_calculator), INTENT(in)                   :: calc
      TYPE(structure_type), INTENT(in)                   :: mol
      TYPE(container_cache), INTENT(in)                  :: ecache
      TYPE(cp2k_stream_plan_type), INTENT(in) :: plan
      REAL(wp), INTENT(out)                              :: gradient(:, :), sigma(:, :)
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      IF (.NOT. stream%active .OR. .NOT. stream%reverse_transaction) THEN
         CALL fatal_error(error, &
            & "g-xTB exchange stream reverse transaction is not active")
         RETURN
      END IF
      IF (.NOT. stream%reverse_applied) THEN
         CALL fatal_error(error, "g-xTB exchange stream reverse is incomplete")
         RETURN
      END IF
      IF (stream%reverse_result_pulled) THEN
         CALL fatal_error(error, "g-xTB exchange reverse result was already pulled")
         RETURN
      END IF
      IF (SIZE(gradient, 1) /= 3 .OR. SIZE(gradient, 2) /= mol%nat &
         & .OR. SIZE(sigma, 1) /= 3 .OR. SIZE(sigma, 2) /= 3) THEN
         CALL fatal_error(error, &
            & "inconsistent g-xTB exchange reverse-result dimensions")
         RETURN
      END IF
      IF (.NOT. exchange_stream_state_matches(stream, calc, mol, ecache, plan)) THEN
         CALL fatal_error(error, &
            & "g-xTB exchange reverse state was invalidated before result")
         RETURN
      END IF
      gradient = stream%reverse_gradient
      sigma = stream%reverse_sigma
      stream%reverse_result_pulled = .TRUE.
   END SUBROUTINE cp2k_exchange_stream_reverse_result

!> Close a successfully applied exchange stream and release all storage.
   SUBROUTINE cp2k_exchange_stream_end(stream, error)
      TYPE(cp2k_exchange_stream_type), INTENT(inout)     :: stream
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      IF (.NOT. stream%active) THEN
         CALL fatal_error(error, "g-xTB exchange stream is not active")
         RETURN
      END IF
      IF (.NOT. stream%applied) THEN
         CALL clear_exchange_stream(stream)
         CALL fatal_error(error, &
            & "g-xTB exchange stream ended before successful apply")
         RETURN
      END IF
      IF (stream%mode == cp2k_exchange_stream_batched &
         & .AND. .NOT. stream%batch_result_pulled) THEN
         CALL clear_exchange_stream(stream)
         CALL fatal_error(error, &
            & "g-xTB exchange stream ended with missing batch result")
         RETURN
      END IF
      IF (stream%mode /= cp2k_exchange_stream_batched &
         & .AND. .NOT. ALL(stream%forward_pulled)) THEN
         CALL clear_exchange_stream(stream)
         CALL fatal_error(error, &
            & "g-xTB exchange stream ended with missing forward pulls")
         RETURN
      END IF
      IF (stream%reverse_transaction .AND. .NOT. stream%reverse_applied) THEN
         CALL clear_exchange_stream(stream)
         CALL fatal_error(error, &
            & "g-xTB exchange stream ended before reverse completion")
         RETURN
      END IF
      IF (stream%reverse_transaction .AND. &
         & .NOT. stream%reverse_result_pulled) THEN
         CALL clear_exchange_stream(stream)
         CALL fatal_error(error, &
            & "g-xTB exchange stream ended with missing reverse result")
         RETURN
      END IF
      IF (stream%reverse_applied) THEN
         IF (stream%partial_pull_distribution_set) THEN
            IF (ANY(stream%reverse_pulled .NEQV. &
               & stream%partial_pull_mask)) THEN
               CALL clear_exchange_stream(stream)
               CALL fatal_error(error, &
                  & "g-xTB exchange stream ended with invalid owned reverse pulls")
               RETURN
            END IF
         ELSE IF (.NOT. ALL(stream%reverse_pulled)) THEN
            CALL clear_exchange_stream(stream)
            CALL fatal_error(error, &
               & "g-xTB exchange stream ended with missing reverse pulls")
            RETURN
         END IF
      END IF
      CALL clear_exchange_stream(stream)
   END SUBROUTINE cp2k_exchange_stream_end

!> Differentiate BvK-supercell-equivalent exchange on a complete k mesh.
!>
!> The overlap response is spin-summed, Hermitian, and unweighted:
!> ``dE = sum_k w_k Re sum_ij conjg(overlap_grad_ij(k))*dS_ij(k)``.
!> Atomic and positive homogeneous-strain derivatives already contain the
!> complete primitive-cell Brillouin-zone contraction and must not be weighted
!> again by CP2K.  The strain response has no inverse-volume factor.
   SUBROUTINE cp2k_exchange_kmesh_gradient(calc, mol, ecache, plan, nmesh, kfrac, &
      & weights, density, overlap, overlap_grad, gradient, sigma, error, &
      & transform_backend)
      TYPE(xtb_calculator), INTENT(in)                   :: calc
      TYPE(structure_type), INTENT(in)                   :: mol
      TYPE(container_cache), INTENT(inout)               :: ecache
      TYPE(cp2k_stream_plan_type), INTENT(inout) :: plan
      INTEGER, INTENT(in)                                :: nmesh(3)
      REAL(wp), INTENT(in)                               :: kfrac(:, :), weights(:)
      COMPLEX(wp), INTENT(in)                            :: density(:, :, :, :), overlap(:, :, :)
      COMPLEX(wp), INTENT(out)                           :: overlap_grad(:, :, :)
      REAL(wp), INTENT(out)                              :: gradient(:, :), sigma(:, :)
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error
      INTEGER, INTENT(in), OPTIONAL                      :: transform_backend

      INTEGER                                            :: backend, idim, nk
      REAL(wp), ALLOCATABLE                              :: bocorr_grad_r(:, :, :), &
                                                            mulliken_grad_r(:, :, :)

      gradient = 0.0_wp
      sigma = 0.0_wp
      overlap_grad = (0.0_wp, 0.0_wp)
      backend = bvk_transform_dense
      IF (PRESENT(transform_backend)) backend = transform_backend
      IF (.NOT. valid_exchange_transform_backend(backend)) THEN
         CALL fatal_error(error, "unknown g-xTB k-to-R transform backend")
         RETURN
      END IF
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
      IF (.NOT. complete_regular_mesh(nmesh, nk)) THEN
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
      IF (compact_exchange_transform_backend(backend) .AND. &
         & MAXVAL(ABS(weights - 1.0_wp/REAL(nk, wp))) > &
         & 64.0_wp*EPSILON(1.0_wp)) THEN
         IF (backend == bvk_transform_separable_dft) THEN
            CALL fatal_error(error, &
               & "separable-DFT g-xTB exchange requires machine-uniform weights")
         ELSE
            CALL fatal_error(error, &
               & "mixed-radix FFT g-xTB exchange requires machine-uniform weights")
         END IF
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

            CALL prepare_bvk_plan(exchange, mol, plan, nmesh, kfrac, weights, &
               & error, retain_dense_phases=backend == bvk_transform_dense)
            IF (ALLOCATED(error)) RETURN
            IF (compact_exchange_transform_backend(backend) .AND. &
               & .NOT. machine_regular_bvk_plan(plan, nmesh, kfrac)) THEN
               IF (backend == bvk_transform_separable_dft) THEN
                  CALL fatal_error(error, &
                     & "separable-DFT g-xTB exchange requires a machine-regular grid")
               ELSE
                  CALL fatal_error(error, &
                     & "mixed-radix FFT g-xTB exchange requires a machine-regular grid")
               END IF
               RETURN
            END IF

            ALLOCATE (mulliken_grad_r(exchange%nao, exchange%nao, nk), &
               & bocorr_grad_r(exchange%nao, exchange%nao, nk))
            CALL consumer_get_kgrad_kmesh(exchange, mol, cache, plan, plan%bvk_kernel, &
               & kfrac, weights, density, overlap, overlap_grad, &
               & mulliken_grad_r, bocorr_grad_r, backend)
            CALL exchange%get_bvk_Kmatrix_derivs(mol, plan%bvk_kernel, &
               & mulliken_grad_r, &
               & bocorr_grad_r, gradient, sigma)
         class default
            CALL fatal_error(error, &
               & "whole-mesh exchange gradient requires exchange_fock")
         END SELECT
      class default
         CALL fatal_error(error, "unexpected g-xTB exchange cache type")
      END SELECT

   END SUBROUTINE cp2k_exchange_kmesh_gradient

!> Apply the exact complete-mesh reverse exchange evaluator once.
!>
!> This transition is legal only after the forward stream was applied.  The
!> complete primitive-cell atomic and positive-strain derivatives are returned
!> directly exactly once.  The unweighted overlap adjoints remain in the
!> stream and must subsequently be pulled blockwise with `reverse_get`.
   SUBROUTINE cp2k_exchange_stream_reverse_apply(stream, calc, mol, ecache, plan, &
      & gradient, sigma, error)
      TYPE(cp2k_exchange_stream_type), INTENT(inout)     :: stream
      TYPE(xtb_calculator), INTENT(in)                   :: calc
      TYPE(structure_type), INTENT(in)                   :: mol
      TYPE(container_cache), INTENT(inout)               :: ecache
      TYPE(cp2k_stream_plan_type), INTENT(inout) :: plan
      REAL(wp), INTENT(out)                              :: gradient(:, :), sigma(:, :)
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      IF (.NOT. stream%active) THEN
         CALL fatal_error(error, "g-xTB exchange stream is not active")
         RETURN
      END IF
      IF (.NOT. stream%applied) THEN
         CALL fatal_error(error, &
            & "g-xTB exchange stream reverse requires forward apply")
         RETURN
      END IF
      IF (stream%mode /= cp2k_exchange_stream_oracle) THEN
         CALL fatal_error(error, &
            & "g-xTB exchange stream reverse requires oracle mode")
         RETURN
      END IF
      IF (stream%reverse_applied) THEN
         CALL fatal_error(error, "g-xTB exchange stream reverse was already applied")
         RETURN
      END IF
      IF (SIZE(gradient, 1) /= 3 .OR. SIZE(gradient, 2) /= mol%nat &
         & .OR. SIZE(sigma, 1) /= 3 .OR. SIZE(sigma, 2) /= 3) THEN
         CALL fatal_error(error, &
            & "inconsistent g-xTB exchange stream reverse dimensions")
         RETURN
      END IF
      IF (.NOT. ALLOCATED(calc%exchange) .OR. .NOT. ALLOCATED(ecache%raw)) THEN
         CALL fatal_error(error, &
            & "g-xTB exchange stream state was invalidated before reverse apply")
         RETURN
      END IF
      IF (.NOT. exchange_stream_state_matches(stream, calc, mol, ecache, plan)) THEN
         CALL fatal_error(error, &
            & "g-xTB exchange stream state was invalidated before reverse apply")
         RETURN
      END IF

      CALL cp2k_exchange_kmesh_gradient(calc, mol, ecache, plan, stream%nmesh, &
         & stream%kfrac, stream%weights, stream%density, stream%overlap, &
         & stream%overlap_adjoint, gradient, sigma, error)
      IF (ALLOCATED(error)) RETURN
      stream%reverse_pulled = .FALSE.
      stream%reverse_applied = .TRUE.
   END SUBROUTINE cp2k_exchange_stream_reverse_apply

!> Pull one unweighted overlap-adjoint block from the reverse stream.
   SUBROUTINE cp2k_exchange_stream_reverse_get(stream, ik, overlap_adjoint, error)
      TYPE(cp2k_exchange_stream_type), INTENT(inout)     :: stream
      INTEGER, INTENT(in)                                :: ik
      COMPLEX(wp), INTENT(out)                           :: overlap_adjoint(:, :)
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      IF (.NOT. stream%active) THEN
         CALL fatal_error(error, "g-xTB exchange stream is not active")
         RETURN
      END IF
      IF (stream%mode /= cp2k_exchange_stream_oracle) THEN
         CALL fatal_error(error, &
            & "blockwise streamed reverse requires reverse-batch pull")
         RETURN
      END IF
      IF (.NOT. stream%reverse_applied) THEN
         CALL fatal_error(error, "g-xTB exchange stream reverse was not applied")
         RETURN
      END IF
      IF (ik < 1 .OR. ik > stream%nk) THEN
         CALL fatal_error(error, &
            & "g-xTB exchange stream reverse block index is out of range")
         RETURN
      END IF
      IF (stream%reverse_pulled(ik)) THEN
         CALL fatal_error(error, "duplicate g-xTB exchange stream reverse pull")
         RETURN
      END IF
      IF (SIZE(overlap_adjoint, 1) /= stream%nao &
         & .OR. SIZE(overlap_adjoint, 2) /= stream%nao) THEN
         CALL fatal_error(error, &
            & "inconsistent g-xTB exchange stream overlap-adjoint dimensions")
         RETURN
      END IF

      overlap_adjoint = stream%overlap_adjoint(:, :, ik)
      stream%reverse_pulled(ik) = .TRUE.
   END SUBROUTINE cp2k_exchange_stream_reverse_get

!> Differentiate the g-xTB exchange functional for one complex k-point block.
!>
!> All returned quantities are unweighted.  For a reduced mesh the external
!> driver must expand `gradient` and `sigma` over the complete k star.  The
!> complex Hermitian `overlap_grad` transforms like a density matrix and is
!> intended to be folded to real-space overlap images by the same weighted
!> star transform as an ordinary k-point density.

!> Complete Gamma-point H0/ACP/exchange and q-vSZP response for CP2K.

!> Image-resolved H0/Pulay and q-vSZP response for CP2K.
!>
!> The real-space matrices use layout `(nao,nao,ntranslation,nspin)` and
!> must use the exact translation ordering supplied here.  ACP and exchange
!> responses remain separate contributions.

!> Extract one effective CGTO in the representation expected by CP2K.

!> Evaluate multipole integrals with save_tblite's explicit basis cache.

!> Evaluate multipole derivatives with save_tblite's explicit basis cache.

END MODULE tblite_exchange_stream

#endif
