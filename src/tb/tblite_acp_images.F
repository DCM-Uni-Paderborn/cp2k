!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: LGPL-3.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
! SPDX-Identifier: LGPL-3.0-or-later
! Image/Bloch assembly moved from our save_tblite integration at 25456f4.
! The original Seidler ACP implementation, integral kernels, parameters,
! sparse Gamma contractions and BLAS/OpenMP policy remain in the library.
#if defined(__TBLITE_GXTB)
MODULE tblite_acp_images
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE mctc_env,                        ONLY: error_type,&
                                              fatal_error,&
                                              wp
   USE mctc_io,                         ONLY: structure_type
   USE mctc_io_constants,               ONLY: pi
   USE mctc_io_math,                    ONLY: matinv_3x3
   USE tblite_acp_cache,                ONLY: acp_cache
   USE tblite_acp_type,                 ONLY: acp_type
   USE tblite_basis_cache,              ONLY: basis_cache
   USE tblite_basis_type,               ONLY: basis_type
   USE tblite_blas,                     ONLY: gemm
   USE tblite_integral_native_integrals,ONLY: msao,&
                                              overlap_cgto,&
                                              overlap_grad_cgto

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: get_acp_images, get_acp_kmesh, get_acp_kmesh_gradient
   PUBLIC :: get_acp_image_gradient, get_pv_overlap_deriv_images

   !> Image representation of the separable ACP operator.
   !>
   !> `projector_overlap(:, :, T)` stores
   !> \f$C_T=\langle p_T|\mu_0\rangle\f$.  Consequently its Bloch
   !> representation is \f$C(k)=\sum_T C_T\exp(-i kT)\f$.  The Hamiltonian
   !> images use the convention expected by CP2K,
   !> \f$H(k)=\sum_R H_R\exp(+i kR)\f$, and are obtained from the exact
   !> finite-image convolution
   !> \f$H_R=\sum_{T-V=R}C_T^T L C_V\f$.
   TYPE, PUBLIC :: acp_image_type
      REAL(wp), ALLOCATABLE :: projector_translation(:, :)
      REAL(wp), ALLOCATABLE :: projector_overlap(:, :, :)
      REAL(wp), ALLOCATABLE :: scaled_projector_overlap(:, :, :)
      REAL(wp), ALLOCATABLE :: level(:)
      REAL(wp), ALLOCATABLE :: translation(:, :)
      INTEGER, ALLOCATABLE :: inverse(:)
      INTEGER :: origin = 0
      LOGICAL :: inversion_closed = .FALSE.
      REAL(wp), ALLOCATABLE :: hamiltonian(:, :, :)
   END TYPE acp_image_type

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param mol ...
!> \param trans ...
!> \param bas ...
!> \param bcache ...
!> \param acp ...
!> \param acache ...
!> \param images ...
!> \param build_hamiltonian ...
! **************************************************************************************************
   SUBROUTINE get_acp_images(mol, trans, bas, bcache, acp, acache, images, &
      & build_hamiltonian)
      !> Molecular structure data
      TYPE(structure_type), INTENT(in)                   :: mol
      REAL(wp), INTENT(in)                               :: trans(:, :)

!> Lattice points used for the projector-valence overlap
      !> Valence basis set information

      class(basis_type), INTENT(in) :: bas
      !> Valence basis set cache
      TYPE(basis_cache), INTENT(in) :: bcache
      !> Atomic correction potential data
      TYPE(acp_type), INTENT(in) :: acp
      !> Atomic correction potential cache
      TYPE(acp_cache), INTENT(inout) :: acache
      !> Image-resolved separable ACP operator
      TYPE(acp_image_type), INTENT(out) :: images
      !> Form the quadratic projector-translation difference set and H(R).
      !> Direct Bloch-space callers only need the compact projector images.
      LOGICAL, INTENT(in), OPTIONAL :: build_hamiltonian

      INTEGER :: iat, jat, isp, jsp, itr, jtr, ktr, ninput, ntr, nhtr
      INTEGER :: ish, jproj, is, js, ii, jj, iao, jao, nao, ij
      LOGICAL :: with_hamiltonian
      REAL(wp) :: diff(3), level, r2, scale, tol2, vec(3)
      REAL(wp), ALLOCATABLE :: htrans_work(:, :), projector_overlap(:, :, :), &
         & scaled_projector_overlap(:, :, :), stmp(:)

      CALL acp%update(mol, acache)

      ninput = SIZE(trans, 2)
      ALLOCATE (projector_overlap(acp%auxbas%nao, bas%nao, ninput), &
         & scaled_projector_overlap(acp%auxbas%nao, bas%nao, ninput), &
         & images%level(acp%auxbas%nao), source=0.0_wp)

      IF (.NOT. ALLOCATED(acache%scaled_pv_overlap)) THEN
         ALLOCATE (acache%scaled_pv_overlap(acp%auxbas%nao, bas%nao), &
            & source=0.0_wp)
      ELSE
         acache%scaled_pv_overlap(:, :) = 0.0_wp
      END IF

      ! Expand the shell levels once to the auxiliary-AO representation.
      DO jat = 1, mol%nat
         jsp = mol%id(jat)
         js = acp%auxbas%ish_at(jat)
         DO jproj = 1, acp%auxbas%nsh_id(jsp)
            jj = acp%auxbas%iao_sh(js + jproj)
            level = acp%levels(jproj, jsp)
            nao = msao(acp%auxbas%cgto(jproj, jsp)%raw%ang)
            images%level(jj + 1:jj + nao) = level
         END DO
      END DO

      ALLOCATE (stmp(msao(bas%maxl)*msao(acp%auxbas%maxl)))

      ! C_T = <p_T|mu_0>.  This is the same integral and finite image set used
      ! by get_pv_overlap_3d, only without summing over T prematurely.
      DO iat = 1, mol%nat
         isp = mol%id(iat)
         is = bas%ish_at(iat)
         DO jat = 1, mol%nat
            jsp = mol%id(jat)
            js = acp%auxbas%ish_at(jat)
            DO itr = 1, ninput
               vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat) - trans(:, itr)
               r2 = DOT_PRODUCT(vec, vec)
               DO ish = 1, bas%nsh_id(isp)
                  ii = bas%iao_sh(is + ish)
                  DO jproj = 1, acp%auxbas%nsh_id(jsp)
                     jj = acp%auxbas%iao_sh(js + jproj)
                     CALL overlap_cgto(acp%auxbas%cgto(jproj, jsp)%raw, &
                        & bas%cgto(ish, isp)%raw, acache%auxbas%cgto(jproj, jat), &
                        & bcache%cgto(ish, iat), r2, vec, acp%auxbas%intcut, stmp)

                     level = acp%levels(jproj, jsp)
                     nao = msao(acp%auxbas%cgto(jproj, jsp)%raw%ang)
                     DO iao = 1, msao(bas%cgto(ish, isp)%raw%ang)
                        DO jao = 1, nao
                           ij = jao + nao*(iao - 1)
                           projector_overlap(jj + jao, ii + iao, itr) = &
                              & stmp(ij)
                           scaled_projector_overlap(jj + jao, ii + iao, itr) = &
                              & level*stmp(ij)
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      ! CP2K may supply a wider image range than the ACP integral support.
      ! overlap_cgto marks images beyond that support with exact zeros; remove
      ! only those empty images so that the direct Bloch representation and the
      ! real-space convolution remain independent of harmless outer padding.
      ntr = 0
      DO itr = 1, ninput
         IF (ANY(projector_overlap(:, :, itr) /= 0.0_wp)) ntr = ntr + 1
      END DO
      ALLOCATE (images%projector_translation(3, ntr), &
         & images%projector_overlap(acp%auxbas%nao, bas%nao, ntr), &
         & images%scaled_projector_overlap(acp%auxbas%nao, bas%nao, ntr))
      ktr = 0
      DO itr = 1, ninput
         IF (.NOT. ANY(projector_overlap(:, :, itr) /= 0.0_wp)) CYCLE
         ktr = ktr + 1
         images%projector_translation(:, ktr) = trans(:, itr)
         images%projector_overlap(:, :, ktr) = projector_overlap(:, :, itr)
         images%scaled_projector_overlap(:, :, ktr) = &
            & scaled_projector_overlap(:, :, itr)
      END DO

      acache%scaled_pv_overlap = SUM(images%scaled_projector_overlap, dim=3)

      with_hamiltonian = .TRUE.
      IF (PRESENT(build_hamiltonian)) with_hamiltonian = build_hamiltonian
      IF (.NOT. with_hamiltonian) RETURN

      ! First construct the unique translation-difference set.  The tolerance
      ! only absorbs floating-point noise from equivalent lattice differences.
      ALLOCATE (htrans_work(3, MAX(1, ntr*ntr)), source=0.0_wp)
      scale = MAX(1.0_wp, MAXVAL(ABS(trans)))
      tol2 = (1024.0_wp*EPSILON(1.0_wp)*scale)**2
      nhtr = 0
      IF (ntr == 0) THEN
         nhtr = 1
         htrans_work(:, nhtr) = 0.0_wp
      ELSE
         DO itr = 1, ntr
            DO jtr = 1, ntr
               diff = images%projector_translation(:, itr) &
                  & - images%projector_translation(:, jtr)
               ktr = find_translation(htrans_work(:, :nhtr), diff, tol2)
               IF (ktr == 0) THEN
                  nhtr = nhtr + 1
                  htrans_work(:, nhtr) = diff
               END IF
            END DO
         END DO
      END IF

      ALLOCATE (images%translation(3, nhtr), &
         & images%hamiltonian(bas%nao, bas%nao, nhtr), source=0.0_wp)
      ALLOCATE (images%inverse(nhtr), source=0)
      images%translation = htrans_work(:, :nhtr)
      images%origin = find_translation(images%translation, &
         & [0.0_wp, 0.0_wp, 0.0_wp], tol2)
      images%inversion_closed = .TRUE.
      DO itr = 1, nhtr
         images%inverse(itr) = find_translation(images%translation, &
            & -images%translation(:, itr), tol2)
         images%inversion_closed = images%inversion_closed .AND. &
            & images%inverse(itr) > 0
      END DO

      ! H_R = sum_(T-V=R) C_T^T L C_V.  Each ordered translation pair is
      ! included, so H_R^T=H_-R and the Bloch sum is Hermitian at arbitrary k.
      DO itr = 1, ntr
         DO jtr = 1, ntr
            diff = images%projector_translation(:, itr) &
               & - images%projector_translation(:, jtr)
            ktr = find_translation(images%translation, diff, tol2)
            CALL gemm(amat=images%projector_overlap(:, :, itr), &
               & bmat=images%scaled_projector_overlap(:, :, jtr), &
               & cmat=images%hamiltonian(:, :, ktr), transa='T', beta=1.0_wp)
         END DO
      END DO

   END SUBROUTINE get_acp_images

! **************************************************************************************************
!> \brief ...
!> \param mol ...
!> \param images ...
!> \param kfrac ...
!> \param hamiltonian ...
!> \param error ...
! **************************************************************************************************
   SUBROUTINE get_acp_kmesh(mol, images, kfrac, hamiltonian, error)
      TYPE(structure_type), INTENT(in)                   :: mol
      TYPE(acp_image_type), INTENT(in)                   :: images
      REAL(wp), INTENT(in)                               :: kfrac(:, :)
      COMPLEX(wp), INTENT(out)                           :: hamiltonian(:, :, :)
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      COMPLEX(wp)                                        :: phase
      COMPLEX(wp), ALLOCATABLE                           :: ck(:, :), csk(:, :)
      INTEGER                                            :: ik, itr, nao, nk, ntr
      REAL(wp)                                           :: angle
      REAL(wp), ALLOCATABLE                              :: projector_cell(:, :)

      IF (.NOT. ALLOCATED(images%projector_translation) .OR. &
         & .NOT. ALLOCATED(images%projector_overlap) .OR. &
         & .NOT. ALLOCATED(images%scaled_projector_overlap)) THEN
         CALL fatal_error(error, "ACP k-mesh evaluation requires projector images")
         RETURN
      END IF
      nao = SIZE(images%projector_overlap, 2)
      ntr = SIZE(images%projector_overlap, 3)
      nk = SIZE(kfrac, 2)
      IF (SIZE(kfrac, 1) /= 3 .OR. nk < 1 .OR. &
         & ANY(SHAPE(images%scaled_projector_overlap) /= &
         & SHAPE(images%projector_overlap)) .OR. &
         & ANY(SHAPE(images%projector_translation) /= [3, ntr]) .OR. &
         & ANY(SHAPE(hamiltonian) /= [nao, nao, nk])) THEN
         CALL fatal_error(error, "inconsistent ACP k-mesh dimensions")
         RETURN
      END IF
      IF (.NOT. ALL(ieee_is_finite(kfrac))) THEN
         CALL fatal_error(error, "ACP k points must be finite")
         RETURN
      END IF

      CALL get_acp_projector_cells(mol, images, projector_cell, error)
      IF (ALLOCATED(error)) RETURN
      ALLOCATE (ck(SIZE(images%projector_overlap, 1), nao), &
         & csk(SIZE(images%projector_overlap, 1), nao))
      hamiltonian = (0.0_wp, 0.0_wp)
      DO ik = 1, nk
         ck = (0.0_wp, 0.0_wp)
         csk = (0.0_wp, 0.0_wp)
         DO itr = 1, ntr
            angle = -2.0_wp*pi*DOT_PRODUCT(kfrac(:, ik), projector_cell(:, itr))
            phase = EXP(CMPLX(0.0_wp, angle, wp))
            ck = ck + phase*images%projector_overlap(:, :, itr)
            csk = csk + phase*images%scaled_projector_overlap(:, :, itr)
         END DO
         hamiltonian(:, :, ik) = MATMUL(CONJG(TRANSPOSE(ck)), csk)
      END DO
   END SUBROUTINE get_acp_kmesh

! **************************************************************************************************
!> \brief ...
!> \param mol ...
!> \param bas ...
!> \param bcache ...
!> \param acp ...
!> \param acache ...
!> \param images ...
!> \param kfrac ...
!> \param weights ...
!> \param density ...
!> \param dEdcnbas ...
!> \param dEdqbas ...
!> \param gradient ...
!> \param sigma ...
!> \param error ...
! **************************************************************************************************
   SUBROUTINE get_acp_kmesh_gradient(mol, bas, bcache, acp, acache, images, &
      & kfrac, weights, density, dEdcnbas, dEdqbas, gradient, sigma, error)
      TYPE(structure_type), INTENT(in)                   :: mol

      class(basis_type), INTENT(in) :: bas
      TYPE(basis_cache), INTENT(in) :: bcache
      TYPE(acp_type), INTENT(in) :: acp
      TYPE(acp_cache), INTENT(in) :: acache
      TYPE(acp_image_type), INTENT(in) :: images
      REAL(wp), INTENT(in) :: kfrac(:, :), weights(:)
      COMPLEX(wp), INTENT(in) :: density(:, :, :, :)
      REAL(wp), INTENT(inout) :: dEdcnbas(:), dEdqbas(:)
      REAL(wp), INTENT(inout) :: gradient(:, :), sigma(:, :)
      TYPE(error_type), ALLOCATABLE, INTENT(out) :: error

      INTEGER :: ik, itr, nao, nk, nspin, ntr, spin
      REAL(wp) :: angle
      REAL(wp), ALLOCATABLE :: projector_cell(:, :), response(:, :, :, :)
      COMPLEX(wp) :: phase
      COMPLEX(wp), ALLOCATABLE :: csk(:, :), work(:, :)

      IF (.NOT. ALLOCATED(images%projector_translation) .OR. &
         & .NOT. ALLOCATED(images%projector_overlap) .OR. &
         & .NOT. ALLOCATED(images%scaled_projector_overlap)) THEN
         CALL fatal_error(error, "ACP k-mesh gradient requires projector images")
         RETURN
      END IF
      nao = SIZE(images%projector_overlap, 2)
      ntr = SIZE(images%projector_overlap, 3)
      nk = SIZE(weights)
      nspin = SIZE(density, 3)
      IF (SIZE(kfrac, 1) /= 3 .OR. SIZE(kfrac, 2) /= nk .OR. nk < 1 .OR. &
         & SIZE(density, 1) /= nao .OR. SIZE(density, 2) /= nao .OR. &
         & SIZE(density, 4) /= nk .OR. nspin < 1 .OR. &
         & ANY(SHAPE(images%scaled_projector_overlap) /= &
         & SHAPE(images%projector_overlap)) .OR. &
         & ANY(SHAPE(images%projector_translation) /= [3, ntr])) THEN
         CALL fatal_error(error, "inconsistent ACP k-mesh gradient dimensions")
         RETURN
      END IF
      IF (.NOT. ALL(ieee_is_finite(kfrac)) .OR. &
         & .NOT. ALL(ieee_is_finite(weights)) .OR. &
         & .NOT. ALL(ieee_is_finite(REAL(density, wp))) .OR. &
         & .NOT. ALL(ieee_is_finite(AIMAG(density))) .OR. &
         & ANY(weights < 0.0_wp) .OR. &
         & ABS(SUM(weights) - 1.0_wp) > 1.0e-12_wp) THEN
         CALL fatal_error(error, "invalid ACP k-mesh density or weights")
         RETURN
      END IF

      CALL get_acp_projector_cells(mol, images, projector_cell, error)
      IF (ALLOCATED(error)) RETURN
      ALLOCATE (response(SIZE(images%projector_overlap, 1), nao, ntr, nspin), &
         & source=0.0_wp)
      ALLOCATE (csk(SIZE(images%projector_overlap, 1), nao), &
         & work(SIZE(images%projector_overlap, 1), nao))
      DO ik = 1, nk
         csk = (0.0_wp, 0.0_wp)
         DO itr = 1, ntr
            angle = -2.0_wp*pi*DOT_PRODUCT(kfrac(:, ik), projector_cell(:, itr))
            phase = EXP(CMPLX(0.0_wp, angle, wp))
            csk = csk + phase*images%scaled_projector_overlap(:, :, itr)
         END DO
         DO spin = 1, nspin
            work = MATMUL(csk, density(:, :, spin, ik))
            DO itr = 1, ntr
               angle = -2.0_wp*pi*DOT_PRODUCT(kfrac(:, ik), projector_cell(:, itr))
               phase = EXP(CMPLX(0.0_wp, angle, wp))
               response(:, :, itr, spin) = response(:, :, itr, spin) + &
                  & 2.0_wp*weights(ik)*REAL(CONJG(phase)*work, wp)
            END DO
         END DO
      END DO

      CALL get_pv_overlap_deriv_images(mol, images%projector_translation, bas, &
         & bcache, acp, acache, response, dEdcnbas, dEdqbas, gradient, sigma)
   END SUBROUTINE get_acp_kmesh_gradient

! **************************************************************************************************
!> \brief ...
!> \param mol ...
!> \param images ...
!> \param cells ...
!> \param error ...
! **************************************************************************************************
   SUBROUTINE get_acp_projector_cells(mol, images, cells, error)
      TYPE(structure_type), INTENT(in)                   :: mol
      TYPE(acp_image_type), INTENT(in)                   :: images
      REAL(wp), ALLOCATABLE, INTENT(out)                 :: cells(:, :)
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      INTEGER                                            :: itr
      REAL(wp)                                           :: invlat(3, 3), scale

      ALLOCATE (cells(3, SIZE(images%projector_translation, 2)), source=0.0_wp)
      IF (.NOT. ANY(mol%periodic)) RETURN
      invlat = matinv_3x3(mol%lattice)
      DO itr = 1, SIZE(cells, 2)
         cells(:, itr) = MATMUL(invlat, images%projector_translation(:, itr))
      END DO
      scale = MAX(1.0_wp, MAXVAL(ABS(cells)))
      IF (MAXVAL(ABS(cells - REAL(NINT(cells), wp))) > &
         & 4096.0_wp*EPSILON(1.0_wp)*scale) THEN
         CALL fatal_error(error, "ACP projector translation is not a lattice vector")
         RETURN
      END IF
      cells = REAL(NINT(cells), wp)
   END SUBROUTINE get_acp_projector_cells

! **************************************************************************************************
!> \brief ...
!> \param trans ...
!> \param TARGET ...
!> \param tol2 ...
!> \return ...
! **************************************************************************************************
   pure FUNCTION find_translation(trans, TARGET, tol2) RESULT(index)
      REAL(wp), INTENT(in)                               :: trans(:, :), TARGET(:), tol2
      INTEGER                                            :: index

      INTEGER                                            :: itr

      index = 0
      DO itr = 1, SIZE(trans, 2)
         IF (SUM((trans(:, itr) - TARGET)**2) <= tol2) THEN
            index = itr
            EXIT
         END IF
      END DO
   END FUNCTION find_translation

! **************************************************************************************************
!> \brief ...
!> \param mol ...
!> \param bas ...
!> \param bcache ...
!> \param acp ...
!> \param acache ...
!> \param images ...
!> \param density_translation ...
!> \param density ...
!> \param dEdcnbas ...
!> \param dEdqbas ...
!> \param gradient ...
!> \param sigma ...
!> \param error ...
! **************************************************************************************************
   SUBROUTINE get_acp_image_gradient(mol, bas, bcache, acp, acache, images, &
      & density_translation, density, dEdcnbas, dEdqbas, gradient, sigma, error)
      !> Molecular structure data
      TYPE(structure_type), INTENT(in)                   :: mol

      !> Valence basis set information
      class(basis_type), INTENT(in) :: bas
      !> Valence basis set cache
      TYPE(basis_cache), INTENT(in) :: bcache
      !> Atomic correction potential data
      TYPE(acp_type), INTENT(in) :: acp
      !> Atomic correction potential cache
      TYPE(acp_cache), INTENT(in) :: acache
      !> Projector and Hamiltonian image data at the current geometry
      TYPE(acp_image_type), INTENT(in) :: images
      !> Cartesian translations of the density images
      REAL(wp), INTENT(in) :: density_translation(:, :)
      !> Density images `(AO, AO, image, spin)` dual to `images%hamiltonian`
      REAL(wp), INTENT(in) :: density(:, :, :, :)
      !> Derivative w.r.t. basis-set coordination numbers
      REAL(wp), INTENT(inout) :: dEdcnbas(:)
      !> Derivative w.r.t. basis-set atomic charges
      REAL(wp), INTENT(inout) :: dEdqbas(:)
      !> Cartesian ACP gradient
      REAL(wp), INTENT(inout) :: gradient(:, :)
      !> ACP strain derivative
      REAL(wp), INTENT(inout) :: sigma(:, :)
      !> Error handling
      TYPE(error_type), ALLOCATABLE, INTENT(out) :: error

      INTEGER :: itr, jtr, ktr, dtr, ntr, nspin, nearest
      INTEGER, ALLOCATABLE :: density_map(:)
      REAL(wp) :: diff(3), nearest2, scale, tol2
      REAL(wp), ALLOCATABLE :: response(:, :, :, :)
      CHARACTER(len=512) :: message

      IF (.NOT. ALLOCATED(images%projector_translation) .OR. &
         & .NOT. ALLOCATED(images%projector_overlap) .OR. &
         & .NOT. ALLOCATED(images%scaled_projector_overlap) .OR. &
         & .NOT. ALLOCATED(images%translation) .OR. &
         & .NOT. ALLOCATED(images%hamiltonian)) THEN
         CALL fatal_error(error, "ACP image gradient requires initialized image data")
         RETURN
      END IF

      ntr = SIZE(images%projector_translation, 2)
      nspin = SIZE(density, 4)
      IF (SIZE(images%projector_overlap, 1) /= acp%auxbas%nao .OR. &
         & SIZE(images%projector_overlap, 2) /= bas%nao .OR. &
         & SIZE(images%projector_overlap, 3) /= ntr .OR. &
         & ANY(SHAPE(images%scaled_projector_overlap) /= &
         & SHAPE(images%projector_overlap))) THEN
         CALL fatal_error(error, "inconsistent ACP projector image dimensions")
         RETURN
      END IF
      IF (SIZE(density, 1) /= bas%nao .OR. SIZE(density, 2) /= bas%nao .OR. &
         & SIZE(density, 3) /= SIZE(density_translation, 2) .OR. &
         & SIZE(density, 3) < 1 .OR. nspin < 1) THEN
         CALL fatal_error(error, "inconsistent ACP density image dimensions")
         RETURN
      END IF
      IF (SIZE(density_translation, 1) /= 3 .OR. &
         & SIZE(images%translation, 1) /= 3 .OR. &
         & SIZE(images%projector_translation, 1) /= 3) THEN
         CALL fatal_error(error, "ACP image translations must be three-dimensional")
         RETURN
      END IF
      IF (SIZE(images%hamiltonian, 1) /= bas%nao .OR. &
         & SIZE(images%hamiltonian, 2) /= bas%nao .OR. &
         & SIZE(images%hamiltonian, 3) /= SIZE(images%translation, 2)) THEN
         CALL fatal_error(error, "inconsistent ACP Hamiltonian image dimensions")
         RETURN
      END IF
      IF (SIZE(dEdcnbas) /= mol%nat .OR. SIZE(dEdqbas) /= mol%nat .OR. &
         & ANY(SHAPE(gradient) /= [3, mol%nat]) .OR. &
         & ANY(SHAPE(sigma) /= [3, 3])) THEN
         CALL fatal_error(error, "inconsistent ACP image gradient output dimensions")
         RETURN
      END IF

      scale = MAX(1.0_wp, MAXVAL(ABS(images%translation)), &
         & MAXVAL(ABS(density_translation)))
      tol2 = (2048.0_wp*EPSILON(1.0_wp)*scale)**2
      ALLOCATE (density_map(SIZE(images%translation, 2)), source=0)
      DO ktr = 1, SIZE(images%translation, 2)
         density_map(ktr) = find_translation(density_translation, &
            & images%translation(:, ktr), tol2)
         IF (density_map(ktr) == 0) THEN
            nearest = 0
            nearest2 = HUGE(1.0_wp)
            DO dtr = 1, SIZE(density_translation, 2)
               diff = density_translation(:, dtr) - images%translation(:, ktr)
               IF (SUM(diff**2) < nearest2) THEN
                  nearest = dtr
                  nearest2 = SUM(diff**2)
               END IF
            END DO
            WRITE (message, '(a,i0,a,3(es16.8,1x),a,i0,a,i0,a,es16.8)') &
               & "density image set does not cover ACP translation ", ktr, &
               & " at [", images%translation(:, ktr), "]; density images=", &
               & SIZE(density_translation, 2), ", nearest=", nearest, &
               & ", distance=", SQRT(nearest2)
            CALL fatal_error(error, TRIM(message))
            RETURN
         END IF
      END DO

      ALLOCATE (response(acp%auxbas%nao, bas%nao, ntr, nspin), source=0.0_wp)
      DO itr = 1, ntr
         DO jtr = 1, ntr
            ! dC_U from the left projector factor:
            !       L C_V P_(U-V)^T
            diff = images%projector_translation(:, itr) &
               & - images%projector_translation(:, jtr)
            ktr = find_translation(images%translation, diff, tol2)
            IF (ktr == 0) THEN
               CALL fatal_error(error, "incomplete ACP Hamiltonian difference set")
               RETURN
            END IF
            dtr = density_map(ktr)
            DO ktr = 1, nspin
               CALL gemm(amat=images%scaled_projector_overlap(:, :, jtr), &
                  & bmat=density(:, :, dtr, ktr), &
                  & cmat=response(:, :, itr, ktr), transb='T', beta=1.0_wp)
            END DO

            ! dC_U from the right projector factor:
            !       L C_T P_(T-U)
            diff = images%projector_translation(:, jtr) &
               & - images%projector_translation(:, itr)
            ktr = find_translation(images%translation, diff, tol2)
            IF (ktr == 0) THEN
               CALL fatal_error(error, "incomplete ACP Hamiltonian difference set")
               RETURN
            END IF
            dtr = density_map(ktr)
            DO ktr = 1, nspin
               CALL gemm(amat=images%scaled_projector_overlap(:, :, jtr), &
                  & bmat=density(:, :, dtr, ktr), &
                  & cmat=response(:, :, itr, ktr), beta=1.0_wp)
            END DO
         END DO
      END DO

      CALL get_pv_overlap_deriv_images(mol, images%projector_translation, bas, &
         & bcache, acp, acache, response, dEdcnbas, dEdqbas, gradient, sigma)
   END SUBROUTINE get_acp_image_gradient

! **************************************************************************************************
!> \brief ...
!> \param mol ...
!> \param trans ...
!> \param bas ...
!> \param bcache ...
!> \param acp ...
!> \param acache ...
!> \param response ...
!> \param dEdcnbas ...
!> \param dEdqbas ...
!> \param gradient ...
!> \param sigma ...
! **************************************************************************************************
   SUBROUTINE get_pv_overlap_deriv_images(mol, trans, bas, bcache, acp, acache, &
      & response, dEdcnbas, dEdqbas, gradient, sigma)
      TYPE(structure_type), INTENT(in)                   :: mol
      REAL(wp), INTENT(in)                               :: trans(:, :)

      class(basis_type), INTENT(in) :: bas
      TYPE(basis_cache), INTENT(in) :: bcache
      TYPE(acp_type), INTENT(in) :: acp
      TYPE(acp_cache), INTENT(in) :: acache
      !> Coefficient of every projector-overlap derivative `(aux, AO, T, spin)`
      REAL(wp), INTENT(in) :: response(:, :, :, :)
      REAL(wp), INTENT(inout) :: dEdcnbas(:)
      REAL(wp), INTENT(inout) :: dEdqbas(:)
      REAL(wp), INTENT(inout) :: gradient(:, :)
      REAL(wp), INTENT(inout) :: sigma(:, :)

      INTEGER :: iat, jat, isp, jsp, spin, nspin, itr
      INTEGER :: ish, jproj, is, js, ii, jj, iao, jao, nao, ij
      REAL(wp) :: r2, vec(3), dG(3), dcnbasi, dqbasi, tmp
      REAL(wp), ALLOCATABLE :: stmp(:), dstmp(:, :)
      REAL(wp), ALLOCATABLE :: dstmpdqeffi(:), dstmpdqeffj(:)
      LOGICAL :: compute_qeff_grad
      REAL(wp), ALLOCATABLE :: dEdcnbas_local(:), dEdqbas_local(:), &
         & gradient_local(:, :), sigma_local(:, :)

      compute_qeff_grad = bas%charge_dependent .OR. acp%auxbas%charge_dependent
      nspin = SIZE(response, 4)
      ALLOCATE (stmp(msao(bas%maxl)*msao(acp%auxbas%maxl)), &
         & dstmp(3, msao(bas%maxl)*msao(acp%auxbas%maxl)))
      IF (compute_qeff_grad) THEN
         ALLOCATE (dstmpdqeffi(msao(bas%maxl)*msao(acp%auxbas%maxl)), &
            & dstmpdqeffj(msao(bas%maxl)*msao(acp%auxbas%maxl)))
      END IF

      DO spin = 1, nspin
         !$OMP PARALLEL DEFAULT(NONE) SHARED(dEdcnbas, dEdqbas, gradient, sigma, trans) &
         !$OMP SHARED(mol, bas, bcache, acp, acache, response, compute_qeff_grad, spin) &
         !$OMP PRIVATE(iat, jat, itr, isp, jsp, is, js, ish, jproj, ii, jj, iao, jao, nao, ij) &
         !$OMP PRIVATE(r2, vec, stmp, dstmp, dstmpdqeffi, dstmpdqeffj, dcnbasi, dqbasi) &
         !$OMP PRIVATE(tmp, dG, dEdcnbas_local, dEdqbas_local, gradient_local, sigma_local)
         IF (compute_qeff_grad) THEN
            ALLOCATE (dEdcnbas_local(SIZE(dEdcnbas)), source=0.0_wp)
            ALLOCATE (dEdqbas_local(SIZE(dEdqbas)), source=0.0_wp)
         END IF
         ALLOCATE (gradient_local(SIZE(gradient, 1), SIZE(gradient, 2)), source=0.0_wp)
         ALLOCATE (sigma_local(SIZE(sigma, 1), SIZE(sigma, 2)), source=0.0_wp)

         !$OMP DO schedule(runtime)
         DO iat = 1, mol%nat
            isp = mol%id(iat)
            is = bas%ish_at(iat)
            dcnbasi = 0.0_wp
            dqbasi = 0.0_wp
            DO jat = 1, mol%nat
               jsp = mol%id(jat)
               js = acp%auxbas%ish_at(jat)
               DO itr = 1, SIZE(trans, 2)
                  vec = mol%xyz(:, iat) - mol%xyz(:, jat) - trans(:, itr)
                  r2 = DOT_PRODUCT(vec, vec)
                  dG = 0.0_wp
                  DO ish = 1, bas%nsh_id(isp)
                     ii = bas%iao_sh(is + ish)
                     DO jproj = 1, acp%auxbas%nsh_id(jsp)
                        jj = acp%auxbas%iao_sh(js + jproj)
                        CALL overlap_grad_cgto(acp%auxbas%cgto(jproj, jsp)%raw, &
                           & bas%cgto(ish, isp)%raw, acache%auxbas%cgto(jproj, jat), &
                           & bcache%cgto(ish, iat), r2, vec, acp%auxbas%intcut, &
                           & stmp, dstmp, dstmpdqeffj, dstmpdqeffi)

                        nao = msao(acp%auxbas%cgto(jproj, jsp)%raw%ang)
                        DO iao = 1, msao(bas%cgto(ish, isp)%raw%ang)
                           DO jao = 1, nao
                              ij = jao + nao*(iao - 1)
                              dG = dG + response(jj + jao, ii + iao, itr, spin) &
                                  & *dstmp(:, ij)
                              IF (compute_qeff_grad) THEN
                                 tmp = response(jj + jao, ii + iao, itr, spin) &
                                      & *dstmpdqeffi(ij)
                                 dcnbasi = dcnbasi &
                                    & + bcache%cgto(ish, iat)%dqeffdcn*tmp
                                 dqbasi = dqbasi &
                                    & + bcache%cgto(ish, iat)%dqeffdq*tmp
                              END IF
                           END DO
                        END DO
                     END DO
                  END DO

                  gradient_local(:, iat) = gradient_local(:, iat) + dG
                  gradient_local(:, jat) = gradient_local(:, jat) - dG
                  sigma_local = sigma_local + 0.5_wp*(SPREAD(vec, 1, 3) &
                                                     & *SPREAD(dG, 2, 3) + SPREAD(dG, 1, 3)*SPREAD(vec, 2, 3))
               END DO
            END DO
            IF (compute_qeff_grad) THEN
               dEdcnbas_local(iat) = dEdcnbas_local(iat) + dcnbasi
               dEdqbas_local(iat) = dEdqbas_local(iat) + dqbasi
            END IF
         END DO
         !$OMP END DO

         !$OMP CRITICAL
         IF (compute_qeff_grad) THEN
            dEdcnbas = dEdcnbas + dEdcnbas_local
            dEdqbas = dEdqbas + dEdqbas_local
         END IF
         gradient = gradient + gradient_local
         sigma = sigma + sigma_local
         !$OMP END CRITICAL
         IF (compute_qeff_grad) THEN
            DEALLOCATE (dEdcnbas_local, dEdqbas_local)
         END IF
         DEALLOCATE (gradient_local, sigma_local)
         !$OMP END PARALLEL
      END DO
   END SUBROUTINE get_pv_overlap_deriv_images

END MODULE tblite_acp_images
#endif
