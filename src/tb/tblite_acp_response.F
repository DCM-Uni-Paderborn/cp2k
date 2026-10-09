!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: LGPL-3.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
! SPDX-Identifier: LGPL-3.0-or-later
! Derived from save_tblite ACP response orchestration at aad26158b2e83f84b7b7d5b56f361ead4bc02a88.
! CP2K owns density contractions, image labels and stream state. Gaussian
! integrals and their physical derivatives remain in the provider.
#if defined(__TBLITE_GXTB_ACCELERATION)
MODULE tblite_acp_response
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE mctc_env,                        ONLY: error_type,&
                                              fatal_error,&
                                              wp
   USE mctc_io,                         ONLY: structure_type
   USE mctc_io_constants,               ONLY: pi
   USE mctc_io_math,                    ONLY: matinv_3x3
   USE tblite_acp_cache,                ONLY: acp_cache
   USE tblite_acp_images,               ONLY: acp_image_type,&
                                              get_pv_overlap_deriv_images
   USE tblite_acp_type,                 ONLY: acp_type
   USE tblite_basis_cache,              ONLY: basis_cache
   USE tblite_basis_type,               ONLY: basis_type
   USE tblite_blas,                     ONLY: gemm

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: acp_kmesh_gradient_stream_begin
   PUBLIC :: acp_kmesh_gradient_stream_push
   PUBLIC :: acp_kmesh_gradient_stream_end
   PUBLIC :: acp_kmesh_gradient_stream_discard
   PUBLIC :: acp_kmesh_gradient_stream_real_elements
   PUBLIC :: acp_kmesh_gradient_stream_complex_elements
   PUBLIC :: get_acp_projector_image_gradient
   TYPE, PUBLIC :: acp_kmesh_gradient_stream_type
      PRIVATE
      LOGICAL :: active = .FALSE.
      INTEGER :: nao = 0
      INTEGER :: naux = 0
      INTEGER :: nk_expected = 0
      INTEGER :: nk_pushed = 0
      INTEGER :: nspin = 0
      INTEGER :: ntr = 0
      REAL(wp) :: weight_sum = 0.0_wp
      REAL(wp), ALLOCATABLE :: projector_cell(:, :)
      REAL(wp), ALLOCATABLE :: response(:, :, :, :)
      COMPLEX(wp), ALLOCATABLE :: csk(:, :)
      COMPLEX(wp), ALLOCATABLE :: work(:, :)
   END TYPE acp_kmesh_gradient_stream_type
CONTAINS
! **************************************************************************************************
!> \brief ...
!> \param stream ...
!> \param mol ...
!> \param images ...
!> \param nspin ...
!> \param nk_expected ...
!> \param error ...
! **************************************************************************************************
   SUBROUTINE acp_kmesh_gradient_stream_begin(stream, mol, images, nspin, &
      & nk_expected, error)
      TYPE(acp_kmesh_gradient_stream_type), &
         INTENT(inout)                                   :: stream
      TYPE(structure_type), INTENT(in)                   :: mol
      TYPE(acp_image_type), INTENT(in)                   :: images
      INTEGER, INTENT(in)                                :: nspin, nk_expected
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      INTEGER                                            :: nao, naux, ntr

      IF (stream%active) THEN
         CALL fatal_error(error, "ACP k-mesh gradient stream is already active")
         RETURN
      END IF
      CALL clear_acp_kmesh_gradient_stream(stream)
      IF (.NOT. ALLOCATED(images%projector_translation) .OR. &
         & .NOT. ALLOCATED(images%projector_overlap) .OR. &
         & .NOT. ALLOCATED(images%scaled_projector_overlap)) THEN
         CALL fatal_error(error, &
            & "ACP k-mesh gradient stream requires projector images")
         RETURN
      END IF
      naux = SIZE(images%projector_overlap, 1)
      nao = SIZE(images%projector_overlap, 2)
      ntr = SIZE(images%projector_overlap, 3)
      IF (naux < 1 .OR. nao < 1 .OR. ntr < 1 .OR. nspin < 1 .OR. &
         & nk_expected < 1 .OR. &
         & ANY(SHAPE(images%scaled_projector_overlap) /= &
         & SHAPE(images%projector_overlap)) .OR. &
         & ANY(SHAPE(images%projector_translation) /= [3, ntr])) THEN
         CALL fatal_error(error, &
            & "inconsistent ACP k-mesh gradient stream dimensions")
         RETURN
      END IF
      IF (.NOT. ALL(ieee_is_finite(images%projector_translation)) .OR. &
         & .NOT. ALL(ieee_is_finite(images%projector_overlap)) .OR. &
         & .NOT. ALL(ieee_is_finite(images%scaled_projector_overlap))) THEN
         CALL fatal_error(error, &
            & "ACP k-mesh gradient stream projector data must be finite")
         RETURN
      END IF

      CALL get_acp_projector_cells(mol, images, stream%projector_cell, error)
      IF (ALLOCATED(error)) RETURN
      ALLOCATE (stream%response(naux, nao, ntr, nspin), source=0.0_wp)
      ALLOCATE (stream%csk(naux, nao), stream%work(naux, nao))
      stream%nao = nao
      stream%naux = naux
      stream%ntr = ntr
      stream%nspin = nspin
      stream%nk_expected = nk_expected
      stream%nk_pushed = 0
      stream%weight_sum = 0.0_wp
      stream%active = .TRUE.
   END SUBROUTINE acp_kmesh_gradient_stream_begin

! **************************************************************************************************
!> \brief ...
!> \param stream ...
!> \param images ...
!> \param ik ...
!> \param kfrac ...
!> \param weight ...
!> \param density ...
!> \param error ...
! **************************************************************************************************
   SUBROUTINE acp_kmesh_gradient_stream_push(stream, images, ik, kfrac, weight, &
      & density, error)
      TYPE(acp_kmesh_gradient_stream_type), &
         INTENT(inout)                                   :: stream
      TYPE(acp_image_type), INTENT(in)                   :: images
      INTEGER, INTENT(in)                                :: ik
      REAL(wp), INTENT(in)                               :: kfrac(:), weight
      COMPLEX(wp), INTENT(in)                            :: density(:, :, :)
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      COMPLEX(wp)                                        :: phase
      INTEGER                                            :: itr, spin
      REAL(wp)                                           :: angle

      IF (.NOT. stream%active) THEN
         CALL fatal_error(error, "ACP k-mesh gradient stream is not active")
         RETURN
      END IF
      IF (ik /= stream%nk_pushed + 1 .OR. ik > stream%nk_expected) THEN
         CALL fatal_error(error, &
            & "ACP k-mesh gradient stream k-point order changed")
         RETURN
      END IF
      IF (.NOT. ALLOCATED(images%projector_translation) .OR. &
         & .NOT. ALLOCATED(images%scaled_projector_overlap)) THEN
         CALL fatal_error(error, &
            & "ACP k-mesh gradient stream lost its projector images")
         RETURN
      END IF
      IF (SIZE(kfrac) /= 3 .OR. &
         & ANY(SHAPE(density) /= [stream%nao, stream%nao, stream%nspin]) .OR. &
         & ANY(SHAPE(images%scaled_projector_overlap) /= &
         & [stream%naux, stream%nao, stream%ntr]) .OR. &
         & ANY(SHAPE(images%projector_translation) /= [3, stream%ntr])) THEN
         CALL fatal_error(error, &
            & "inconsistent ACP k-mesh gradient stream block dimensions")
         RETURN
      END IF
      IF (.NOT. ALL(ieee_is_finite(kfrac)) .OR. &
         & .NOT. ieee_is_finite(weight) .OR. weight < 0.0_wp .OR. &
         & .NOT. ALL(ieee_is_finite(REAL(density, wp))) .OR. &
         & .NOT. ALL(ieee_is_finite(AIMAG(density)))) THEN
         CALL fatal_error(error, "invalid ACP k-mesh gradient stream block")
         RETURN
      END IF

      stream%csk = (0.0_wp, 0.0_wp)
      DO itr = 1, stream%ntr
         angle = -2.0_wp*pi*DOT_PRODUCT(kfrac, &
            & stream%projector_cell(:, itr))
         phase = EXP(CMPLX(0.0_wp, angle, wp))
         stream%csk = stream%csk + &
            & phase*images%scaled_projector_overlap(:, :, itr)
      END DO
      DO spin = 1, stream%nspin
         stream%work = MATMUL(stream%csk, density(:, :, spin))
         DO itr = 1, stream%ntr
            angle = -2.0_wp*pi*DOT_PRODUCT(kfrac, &
               & stream%projector_cell(:, itr))
            phase = EXP(CMPLX(0.0_wp, angle, wp))
            stream%response(:, :, itr, spin) = &
               & stream%response(:, :, itr, spin) + &
               & 2.0_wp*weight*REAL(CONJG(phase)*stream%work, wp)
         END DO
      END DO
      stream%nk_pushed = ik
      stream%weight_sum = stream%weight_sum + weight
   END SUBROUTINE acp_kmesh_gradient_stream_push

! **************************************************************************************************
!> \brief ...
!> \param stream ...
!> \param mol ...
!> \param bas ...
!> \param bcache ...
!> \param acp ...
!> \param acache ...
!> \param images ...
!> \param dEdcnbas ...
!> \param dEdqbas ...
!> \param gradient ...
!> \param sigma ...
!> \param error ...
! **************************************************************************************************
   SUBROUTINE acp_kmesh_gradient_stream_end(stream, mol, bas, bcache, acp, &
      & acache, images, dEdcnbas, dEdqbas, gradient, sigma, error)
      TYPE(acp_kmesh_gradient_stream_type), &
         INTENT(inout)                                   :: stream
      TYPE(structure_type), INTENT(in)                   :: mol

      class(basis_type), INTENT(in) :: bas
      TYPE(basis_cache), INTENT(in) :: bcache
      TYPE(acp_type), INTENT(in) :: acp
      TYPE(acp_cache), INTENT(in) :: acache
      TYPE(acp_image_type), INTENT(in) :: images
      REAL(wp), INTENT(inout) :: dEdcnbas(:), dEdqbas(:)
      REAL(wp), INTENT(inout) :: gradient(:, :), sigma(:, :)
      TYPE(error_type), ALLOCATABLE, INTENT(out) :: error

      IF (.NOT. stream%active) THEN
         CALL fatal_error(error, "ACP k-mesh gradient stream is not active")
         RETURN
      END IF
      IF (stream%nk_pushed /= stream%nk_expected) THEN
         CALL fatal_error(error, &
            & "ACP k-mesh gradient stream ended with missing k points")
         RETURN
      END IF
      IF (ABS(stream%weight_sum - 1.0_wp) > 1.0e-12_wp) THEN
         CALL fatal_error(error, &
            & "ACP k-mesh gradient stream weights are not normalized")
         RETURN
      END IF
      IF (.NOT. ALLOCATED(images%projector_translation)) THEN
         CALL fatal_error(error, &
            & "ACP k-mesh gradient stream lost its projector translations")
         RETURN
      END IF
      IF (bas%nao /= stream%nao .OR. acp%auxbas%nao /= stream%naux .OR. &
         & ANY(SHAPE(images%projector_translation) /= [3, stream%ntr]) .OR. &
         & SIZE(dEdcnbas) /= mol%nat .OR. SIZE(dEdqbas) /= mol%nat .OR. &
         & ANY(SHAPE(gradient) /= [3, mol%nat]) .OR. &
         & ANY(SHAPE(sigma) /= [3, 3])) THEN
         CALL fatal_error(error, &
            & "inconsistent ACP k-mesh gradient stream output dimensions")
         RETURN
      END IF

      CALL get_pv_overlap_deriv_images(mol, images%projector_translation, bas, &
         & bcache, acp, acache, stream%response, dEdcnbas, dEdqbas, gradient, &
         & sigma)
      CALL clear_acp_kmesh_gradient_stream(stream)
   END SUBROUTINE acp_kmesh_gradient_stream_end

! **************************************************************************************************
!> \brief ...
!> \param stream ...
! **************************************************************************************************
   SUBROUTINE acp_kmesh_gradient_stream_discard(stream)
      TYPE(acp_kmesh_gradient_stream_type), &
         INTENT(inout)                                   :: stream

      CALL clear_acp_kmesh_gradient_stream(stream)
   END SUBROUTINE acp_kmesh_gradient_stream_discard

! **************************************************************************************************
!> \brief ...
!> \param stream ...
!> \return ...
! **************************************************************************************************
   INTEGER FUNCTION acp_kmesh_gradient_stream_real_elements(stream) &
      & RESULT(nelements)
      TYPE(acp_kmesh_gradient_stream_type), INTENT(in)   :: stream

      nelements = 0
      IF (ALLOCATED(stream%projector_cell)) THEN
         nelements = nelements + SIZE(stream%projector_cell)
      END IF
      IF (ALLOCATED(stream%response)) nelements = nelements + SIZE(stream%response)
   END FUNCTION acp_kmesh_gradient_stream_real_elements

! **************************************************************************************************
!> \brief ...
!> \param stream ...
!> \return ...
! **************************************************************************************************
   INTEGER FUNCTION acp_kmesh_gradient_stream_complex_elements(stream) &
      & RESULT(nelements)
      TYPE(acp_kmesh_gradient_stream_type), INTENT(in)   :: stream

      nelements = 0
      IF (ALLOCATED(stream%csk)) nelements = nelements + SIZE(stream%csk)
      IF (ALLOCATED(stream%work)) nelements = nelements + SIZE(stream%work)
   END FUNCTION acp_kmesh_gradient_stream_complex_elements

! **************************************************************************************************
!> \brief ...
!> \param stream ...
! **************************************************************************************************
   SUBROUTINE clear_acp_kmesh_gradient_stream(stream)
      TYPE(acp_kmesh_gradient_stream_type), &
         INTENT(inout)                                   :: stream

      IF (ALLOCATED(stream%projector_cell)) DEALLOCATE (stream%projector_cell)
      IF (ALLOCATED(stream%response)) DEALLOCATE (stream%response)
      IF (ALLOCATED(stream%csk)) DEALLOCATE (stream%csk)
      IF (ALLOCATED(stream%work)) DEALLOCATE (stream%work)
      stream%active = .FALSE.
      stream%nao = 0
      stream%naux = 0
      stream%nk_expected = 0
      stream%nk_pushed = 0
      stream%nspin = 0
      stream%ntr = 0
      stream%weight_sum = 0.0_wp
   END SUBROUTINE clear_acp_kmesh_gradient_stream

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
   SUBROUTINE get_acp_projector_image_gradient(mol, bas, bcache, acp, acache, &
      & images, density_translation, density, dEdcnbas, dEdqbas, gradient, &
      & sigma, error)
      TYPE(structure_type), INTENT(in)                   :: mol

      class(basis_type), INTENT(in) :: bas
      TYPE(basis_cache), INTENT(in) :: bcache
      TYPE(acp_type), INTENT(in) :: acp
      TYPE(acp_cache), INTENT(in) :: acache
      TYPE(acp_image_type), INTENT(in) :: images
      REAL(wp), INTENT(in) :: density_translation(:, :)
      REAL(wp), INTENT(in) :: density(:, :, :, :)
      REAL(wp), INTENT(inout) :: dEdcnbas(:), dEdqbas(:)
      REAL(wp), INTENT(inout) :: gradient(:, :), sigma(:, :)
      TYPE(error_type), ALLOCATABLE, INTENT(out) :: error

      INTEGER :: dtr, itr, jtr, ktr, nspin, ntr, spin
      INTEGER, ALLOCATABLE :: density_map(:)
      REAL(wp) :: diff(3), response_residual, scale, tol2
      REAL(wp), ALLOCATABLE :: density_cell(:, :), projector_cell(:, :), &
         & response(:, :, :, :), response_oracle(:, :, :, :)
      CHARACTER(len=256) :: message

      IF (.NOT. ALLOCATED(images%projector_translation) .OR. &
         & .NOT. ALLOCATED(images%projector_overlap) .OR. &
         & .NOT. ALLOCATED(images%scaled_projector_overlap)) THEN
         CALL fatal_error(error, &
            & "ACP sparse image gradient requires compact projector images")
         RETURN
      END IF
      ntr = SIZE(images%projector_translation, 2)
      nspin = SIZE(density, 4)
      IF (SIZE(images%projector_overlap, 1) /= acp%auxbas%nao .OR. &
         & SIZE(images%projector_overlap, 2) /= bas%nao .OR. &
         & SIZE(images%projector_overlap, 3) /= ntr .OR. &
         & ANY(SHAPE(images%scaled_projector_overlap) /= &
         & SHAPE(images%projector_overlap)) .OR. &
         & ANY(SHAPE(images%projector_translation) /= [3, ntr]) .OR. &
         & SIZE(density_translation, 1) /= 3 .OR. &
         & SIZE(density, 1) /= bas%nao .OR. SIZE(density, 2) /= bas%nao .OR. &
         & SIZE(density, 3) /= SIZE(density_translation, 2) .OR. &
         & SIZE(density, 3) < 1 .OR. nspin < 1 .OR. &
         & SIZE(dEdcnbas) /= mol%nat .OR. SIZE(dEdqbas) /= mol%nat .OR. &
         & ANY(SHAPE(gradient) /= [3, mol%nat]) .OR. &
         & ANY(SHAPE(sigma) /= [3, 3])) THEN
         CALL fatal_error(error, &
            & "inconsistent ACP sparse image gradient dimensions")
         RETURN
      END IF
      IF (.NOT. ALL(ieee_is_finite(density_translation)) .OR. &
         & .NOT. ALL(ieee_is_finite(density))) THEN
         CALL fatal_error(error, &
            & "ACP sparse image gradient inputs must be finite")
         RETURN
      END IF

      CALL get_acp_translation_cells(mol, images%projector_translation, &
         & projector_cell, error)
      IF (ALLOCATED(error)) RETURN
      CALL get_acp_translation_cells(mol, density_translation, density_cell, error)
      IF (ALLOCATED(error)) RETURN
      ALLOCATE (response(acp%auxbas%nao, bas%nao, ntr, nspin), source=0.0_wp)
      DO itr = 1, ntr
         DO jtr = 1, ntr
            diff = projector_cell(:, itr) - projector_cell(:, jtr)
            dtr = find_translation(density_cell, diff, 0.0_wp)
            IF (dtr > 0) THEN
               DO spin = 1, nspin
                  CALL gemm(amat=images%scaled_projector_overlap(:, :, jtr), &
                     & bmat=density(:, :, dtr, spin), &
                     & cmat=response(:, :, itr, spin), beta=1.0_wp)
               END DO
            END IF

            diff = projector_cell(:, jtr) - projector_cell(:, itr)
            dtr = find_translation(density_cell, diff, 0.0_wp)
            IF (dtr > 0) THEN
               DO spin = 1, nspin
                  CALL gemm(amat=images%scaled_projector_overlap(:, :, jtr), &
                     & bmat=density(:, :, dtr, spin), &
                     & cmat=response(:, :, itr, spin), transb='T', beta=1.0_wp)
               END DO
            END IF
         END DO
      END DO

      IF (ALLOCATED(images%translation) .AND. ALLOCATED(images%hamiltonian)) THEN
         scale = MAX(1.0_wp, MAXVAL(ABS(images%translation)), &
            & MAXVAL(ABS(density_translation)))
         tol2 = (2048.0_wp*EPSILON(1.0_wp)*scale)**2
         ALLOCATE (density_map(SIZE(images%translation, 2)), source=0)
         DO ktr = 1, SIZE(images%translation, 2)
            density_map(ktr) = find_translation(density_translation, &
               & images%translation(:, ktr), tol2)
            IF (density_map(ktr) == 0) THEN
               CALL fatal_error(error, &
                  & "ACP sparse qualification lacks a dense-oracle translation")
               RETURN
            END IF
         END DO
         ALLOCATE (response_oracle, mold=response)
         response_oracle = 0.0_wp
         DO itr = 1, ntr
            DO jtr = 1, ntr
               diff = images%projector_translation(:, itr) &
                  & - images%projector_translation(:, jtr)
               ktr = find_translation(images%translation, diff, tol2)
               dtr = density_map(ktr)
               DO spin = 1, nspin
                  CALL gemm(amat=images%scaled_projector_overlap(:, :, jtr), &
                     & bmat=density(:, :, dtr, spin), &
                     & cmat=response_oracle(:, :, itr, spin), beta=1.0_wp)
               END DO
               diff = images%projector_translation(:, jtr) &
                  & - images%projector_translation(:, itr)
               ktr = find_translation(images%translation, diff, tol2)
               dtr = density_map(ktr)
               DO spin = 1, nspin
                  CALL gemm(amat=images%scaled_projector_overlap(:, :, jtr), &
                     & bmat=density(:, :, dtr, spin), &
                     & cmat=response_oracle(:, :, itr, spin), transb='T', &
                     & beta=1.0_wp)
               END DO
            END DO
         END DO
         IF (.NOT. ALL(ieee_is_finite(response)) .OR. &
            & .NOT. ALL(ieee_is_finite(response_oracle))) THEN
            CALL fatal_error(error, &
               & "ACP sparse qualification produced a non-finite response")
            RETURN
         END IF
         scale = MAX(1.0_wp, MAXVAL(ABS(response_oracle)))
         response_residual = MAXVAL(ABS(response - response_oracle))/scale
         WRITE (*, '(a,es12.4)') &
            & " ACP_SPARSE_RESPONSE_QUALIFY residual=", response_residual
         IF (response_residual > 1.0e-12_wp) THEN
            WRITE (message, '(a,es24.16,a,es24.16)') &
               & "ACP sparse response exceeded its dense-oracle gate: relative=", &
               & response_residual, ", scale=", scale
            CALL fatal_error(error, TRIM(message))
            RETURN
         END IF
      END IF

      CALL get_pv_overlap_deriv_images(mol, images%projector_translation, bas, &
         & bcache, acp, acache, response, dEdcnbas, dEdqbas, gradient, sigma)
   END SUBROUTINE get_acp_projector_image_gradient

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

      CALL get_acp_translation_cells(mol, images%projector_translation, cells, &
         & error)
   END SUBROUTINE get_acp_projector_cells

! **************************************************************************************************
!> \brief ...
!> \param mol ...
!> \param translation ...
!> \param cells ...
!> \param error ...
! **************************************************************************************************
   SUBROUTINE get_acp_translation_cells(mol, translation, cells, error)
      TYPE(structure_type), INTENT(in)                   :: mol
      REAL(wp), INTENT(in)                               :: translation(:, :)
      REAL(wp), ALLOCATABLE, INTENT(out)                 :: cells(:, :)
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      INTEGER                                            :: itr
      REAL(wp)                                           :: invlat(3, 3), scale

      ALLOCATE (cells(3, SIZE(translation, 2)), source=0.0_wp)
      IF (.NOT. ANY(mol%periodic)) RETURN
      invlat = matinv_3x3(mol%lattice)
      DO itr = 1, SIZE(cells, 2)
         cells(:, itr) = MATMUL(invlat, translation(:, itr))
      END DO
      scale = MAX(1.0_wp, MAXVAL(ABS(cells)))
      IF (MAXVAL(ABS(cells - REAL(NINT(cells), wp))) > &
         & 4096.0_wp*EPSILON(1.0_wp)*scale) THEN
         CALL fatal_error(error, "ACP translation is not a lattice vector")
         RETURN
      END IF
      cells = REAL(NINT(cells), wp)
   END SUBROUTINE get_acp_translation_cells

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
END MODULE tblite_acp_response
#endif
