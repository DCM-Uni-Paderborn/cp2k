!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Localized corrections to native Gaussian RI transition integrals.
!>        Each original difference field is paired with the complete soft-plus-Gaussian field.
!>        This retains products between different centers, without a spatial partition or taper.
!>        The RI owner supplies its Gaussian shell, AO block and image coordinates. Its metric,
!>        conditioning, occupations and exchange interaction are not changed here.
! **************************************************************************************************
MODULE qs_upf_inverse_hfx
   USE basis_set_types,                 ONLY: gto_basis_set_p_type,&
                                              gto_basis_set_type
   USE dbt_api,                         ONLY: &
        dbt_clear, dbt_create, dbt_destroy, dbt_finalize, dbt_get_block, dbt_get_info, &
        dbt_get_stored_coordinates, dbt_put_block, dbt_type
   USE gaussian_basis_fields,           ONLY: gaussian_basis_fields_periodic
   USE gaussian_coulomb_fields,         ONLY: gaussian_coulomb_fields_evaluate,&
                                              gaussian_coulomb_fields_init,&
                                              gaussian_coulomb_fields_periodic,&
                                              gaussian_coulomb_fields_prepare,&
                                              gaussian_coulomb_fields_screen,&
                                              gaussian_coulomb_fields_type
   USE input_constants,                 ONLY: do_potential_coulomb,&
                                              do_potential_long,&
                                              do_potential_short
   USE integral_library_types,          ONLY: coulomb_operator_type
   USE kinds,                           ONLY: dp
   USE lebedev,                         ONLY: get_number_of_lebedev_grid,&
                                              init_lebedev_grids,&
                                              lebedev_grid
   USE libint_2c_3c,                    ONLY: cutoff_screen_factor
   USE mathconstants,                   ONLY: fourpi,&
                                              twopi
   USE mathlib,                         ONLY: inv_3x3
   USE paw_onecenter,                   ONLY: paw_valence_augmentation_init
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_upf_inverse,                  ONLY: qs_upf_inverse_atomic_overlap
   USE qs_upf_inverse_types,            ONLY: qs_upf_inverse_type
   USE upf_inverse_exchange,            ONLY: upf_inverse_exchange_augmentation
   USE upf_radial_fields,               ONLY: upf_radial_fields_periodic
   USE upf_types,                       ONLY: upf_radial_weights
   USE uspp_augmentation,               ONLY: uspp_augmentation_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: qs_upf_inverse_hfx_center, qs_upf_inverse_hfx_tensor, qs_upf_inverse_hfx_contract
   PUBLIC :: qs_upf_inverse_hfx_center_adjoint, qs_upf_inverse_hfx_lattice_center
   PUBLIC :: qs_upf_inverse_hfx_lattice_tensor, qs_upf_inverse_hfx_lattice_gradient
   PUBLIC :: qs_upf_inverse_hfx_tensor_adjoint
   PUBLIC :: qs_upf_inverse_hfx_gamma_images

   ! Scratch fields are recomputed at every radial row and retained only within one center call.
   TYPE inverse_field_workspace_type
      COMPLEX(dp), ALLOCATABLE :: delta(:, :), ddelta(:, :), single_delta(:, :), single_ddelta(:, :), &
                                  phi(:, :), dphi(:, :)
   END TYPE inverse_field_workspace_type

   ! Relative-image projector overlaps do not depend on the auxiliary Gaussian shell.
   ! Each entry lives only within one physical image, AO block pair and geometry call.
   TYPE lattice_overlap_type
      COMPLEX(dp), ALLOCATABLE :: left(:, :), right(:, :), dleft(:, :), dright(:, :)
   END TYPE lattice_overlap_type

   ! Original Q integrals are independent of the AO block and its coefficients.
   ! Entries are local to one RI block, geometry and tensor/adjoint call.
   TYPE transition_augmentation_type
      REAL(dp), ALLOCATABLE :: values(:, :, :), gradient(:, :, :, :), cell_gradient(:, :, :, :, :)
   END TYPE transition_augmentation_type

CONTAINS

! **************************************************************************************************
!> \brief Auxiliary images covering the native erfc range and original localized transition support.
!>        Each original center contributes once in the reference cell. Its localized integral
!>        uses the summed auxiliary field. These bounds follow the same physical interaction
!>        range as native Gaussian integrals and retain original Q or PAW valence support.
!>        Forward and geometry contractions use the same integer images at a fixed geometry.
!> \param frame Original Gamma reconstruction and physical cell.
!> \param kinds Original radial datasets.
!> \param augmentation Supplied USPP Q, distinct from PAW valence transitions.
!> \param auxiliary_bases Native RI bases in atom order with their native screening radii.
!> \param metric Native fitting interaction and prescribed erfc screening range.
!> \param images Fixed auxiliary translations, including zero.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_hfx_gamma_images(frame, kinds, augmentation, auxiliary_bases, metric, images)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      TYPE(qs_kind_type), INTENT(IN)                     :: kinds(:)
      TYPE(uspp_augmentation_type), INTENT(IN)           :: augmentation(:)
      TYPE(gto_basis_set_p_type), INTENT(IN)             :: auxiliary_bases(:)
      TYPE(coulomb_operator_type), INTENT(IN)            :: metric
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: images(:, :)

      INTEGER                                            :: atom, axis, center, image, ix, iy, iz, &
                                                            kind, lower(3), upper(3)
      REAL(dp)                                           :: fractional(3), inverse_cell(3, 3), &
                                                            lower_real(3), margin, radius, reach, &
                                                            support, upper_real(3)
      TYPE(uspp_augmentation_type), ALLOCATABLE          :: paw_augmentation(:)

      CPASSERT(ALL(frame%reduced_k == 0.0_dp))
      CPASSERT(SIZE(auxiliary_bases) == SIZE(frame%positions, 2))
      CPASSERT(SIZE(augmentation) == SIZE(kinds))
      IF (ALL(frame%periodicity == 0)) THEN
         ALLOCATE (images(3, 1))
         images = 0
         RETURN
      END IF
      IF (metric%potential_type /= do_potential_short) THEN
         CPABORT("Periodic Gamma UPF fitting currently requires erfc")
      END IF
      CPASSERT(metric%omega > 0.0_dp .AND. metric%cutoff_radius > 0.0_dp)
      inverse_cell = inv_3x3(frame%cell)
      CALL prepare_paw_transitions(kinds, frame, paw_augmentation)
      lower_real = 0.0_dp
      upper_real = 0.0_dp
      DO center = 1, SIZE(frame%center_atoms)
         kind = frame%center_kinds(center)
         support = MAX(frame%delta_cache(kind)%radius, augmentation(kind)%realspace_radius, &
                       paw_augmentation(kind)%realspace_radius)
         DO atom = 1, SIZE(auxiliary_bases)
            radius = support + auxiliary_bases(atom)%gto_basis_set%kind_radius + &
                     cutoff_screen_factor*metric%cutoff_radius
            fractional = MATMUL(inverse_cell, &
                                frame%positions(:, frame%center_atoms(center)) - frame%positions(:, atom))
            DO axis = 1, 3
               IF (frame%periodicity(axis) == 0) CYCLE
               reach = radius*NORM2(inverse_cell(axis, :))
               margin = 32.0_dp*EPSILON(1.0_dp)*MAX(1.0_dp, ABS(fractional(axis)), reach)
               lower_real(axis) = MIN(lower_real(axis), fractional(axis) - reach - margin)
               upper_real(axis) = MAX(upper_real(axis), fractional(axis) + reach + margin)
            END DO
         END DO
      END DO
      lower = FLOOR(lower_real)
      upper = CEILING(upper_real)
      ALLOCATE (images(3, PRODUCT(upper - lower + 1)))
      image = 0
      DO iz = lower(3), upper(3)
         DO iy = lower(2), upper(2)
            DO ix = lower(1), upper(1)
               image = image + 1
               images(:, image) = [ix, iy, iz]
            END DO
         END DO
      END DO
   END SUBROUTINE qs_upf_inverse_hfx_gamma_images

! **************************************************************************************************
!> \brief Add original transition corrections to physical AO image tensors before local RI fitting.
!>        Auxiliary sites are stacked by cell and then by atom, in the native RI order. AO images
!>        retain the supplied order. Reconstruction images are summed before any owner filtering.
!>        The common frame supplies original radial fields and geometry, never its Bloch-summed W.
!>        Its defining lattice domain is supplied explicitly. This routine does not select or
!>        truncate images, change the fitting domain, or apply a Bloch phase or exchange weight.
!> \param frame Common DIRECT or INVERSE_1C frame with the original basis and lattice geometry.
!> \param kinds ...
!> \param augmentation ...
!> \param auxiliary_bases Native per-atom RI bases.
!> \param auxiliary_images Integer RI cells in native stacked order.
!> \param ao_images Integer right AO cells corresponding to tensors.
!> \param center_images Original reconstruction cells.
!> \param tensors Unfiltered (RI,AO0,AOR) tensors, including initially absent Gaussian blocks.
!> \param responses Optional matching integral response tensors.
!> \param atom_response Motion of the unshifted physical atoms and their RI centers.
!> \param cell_response Lattice-vector response at fixed integer cells and reduced k points.
!> \param metric Optional native RI interaction, Coulomb by default.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_hfx_lattice_tensor(frame, kinds, augmentation, auxiliary_bases, auxiliary_images, &
                                                ao_images, center_images, tensors, responses, atom_response, cell_response, metric)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      TYPE(qs_kind_type), INTENT(IN)                     :: kinds(:)
      TYPE(uspp_augmentation_type), INTENT(IN), TARGET   :: augmentation(:)
      TYPE(gto_basis_set_p_type), INTENT(IN)             :: auxiliary_bases(:)
      INTEGER, INTENT(IN)                                :: auxiliary_images(:, :), ao_images(:, :), &
                                                            center_images(:, :)
      TYPE(dbt_type), INTENT(INOUT)                      :: tensors(:)
      TYPE(dbt_type), INTENT(INOUT), OPTIONAL            :: responses(:)
      REAL(dp), INTENT(IN), OPTIONAL                     :: atom_response(:, :), cell_response(3, 3)
      TYPE(coulomb_operator_type), INTENT(IN), OPTIONAL  :: metric

      COMPLEX(dp), ALLOCATABLE                           :: factor_response(:), factors(:)
      INTEGER                                            :: atom, image, natom, nsite, site
      REAL(dp), ALLOCATABLE                              :: position_response(:, :), positions(:, :)
      TYPE(gto_basis_set_p_type), ALLOCATABLE            :: bases(:)

      CPASSERT(frame%atom_local)
      natom = SIZE(frame%positions, 2)
      CPASSERT(SIZE(auxiliary_bases) == natom)
      CPASSERT(SIZE(auxiliary_images, 1) == 3 .AND. SIZE(auxiliary_images, 2) > 0)
      CPASSERT(SIZE(ao_images, 1) == 3 .AND. SIZE(ao_images, 2) == SIZE(tensors))
      CPASSERT(SIZE(center_images, 1) == 3 .AND. SIZE(center_images, 2) > 0)
      CPASSERT(SIZE(tensors) > 0)
      IF (.NOT. (ALL(auxiliary_images*SPREAD(frame%periodicity, 2, SIZE(auxiliary_images, 2)) == auxiliary_images))) THEN
         CPABORT("Auxiliary images violate the periodicity")
      END IF
      CPASSERT(ALL(ao_images*SPREAD(frame%periodicity, 2, SIZE(ao_images, 2)) == ao_images))
      IF (PRESENT(responses)) THEN
         CPASSERT(SIZE(responses) == SIZE(tensors))
         CPASSERT(PRESENT(atom_response) .AND. PRESENT(cell_response))
         CPASSERT(ALL(SHAPE(atom_response) == SHAPE(frame%positions)))
      ELSE
         CPASSERT(.NOT. PRESENT(atom_response) .AND. .NOT. PRESENT(cell_response))
      END IF
      nsite = natom*SIZE(auxiliary_images, 2)
      ALLOCATE (bases(nsite), positions(3, nsite), factors(nsite))
      factors = CMPLX(1.0_dp, 0.0_dp, dp)
      IF (PRESENT(responses)) THEN
         ALLOCATE (position_response(3, nsite), factor_response(nsite))
         factor_response = CMPLX(0.0_dp, 0.0_dp, dp)
      END IF
      DO image = 1, SIZE(auxiliary_images, 2)
      DO atom = 1, natom
         site = (image - 1)*natom + atom
         bases(site)%gto_basis_set => auxiliary_bases(atom)%gto_basis_set
         positions(:, site) = frame%positions(:, atom) + MATMUL(frame%cell, REAL(auxiliary_images(:, image), dp))
         IF (PRESENT(responses)) THEN
            position_response(:, site) = atom_response(:, atom) + MATMUL(cell_response, REAL(auxiliary_images(:, image), dp))
         END IF
      END DO
      END DO
      DO image = 1, SIZE(tensors)
         IF (PRESENT(responses)) THEN
            CALL qs_upf_inverse_hfx_tensor(frame, frame, kinds, augmentation, bases, positions, factors, center_images, &
                                           tensors(image), response_real=responses(image), &
                                           atom_response=atom_response, cell_response=cell_response, &
                                           auxiliary_positions_response=position_response, &
                                           auxiliary_factors_response=factor_response, ao_image=ao_images(:, image), metric=metric)
         ELSE
            CALL qs_upf_inverse_hfx_tensor(frame, frame, kinds, augmentation, bases, positions, factors, center_images, &
                                           tensors(image), ao_image=ao_images(:, image), metric=metric)
         END IF
      END DO
   END SUBROUTINE qs_upf_inverse_hfx_lattice_tensor

! **************************************************************************************************
!> \brief Contract complete physical AO image covectors with original transition geometry responses.
!>        The native owner supplies the derivative before restricting Gaussian derivative blocks,
!>        including local RI inverse, bump, density, spin and exchange weights. The same cell/atom
!>        auxiliary ordering and complete center sums as in the forward tensor are retained.
!>        The analytic DIRECT or INVERSE_1C projection response is included. Fitting-metric and
!>        ordinary Gaussian derivatives remain native. Outputs are summed once and replicated
!>        on the supplied tensor communicator, which may be a native RI subgroup.
!> \param frame Original atom-local reconstruction and physical geometry.
!> \param kinds ...
!> \param augmentation ...
!> \param auxiliary_bases Native per-atom RI bases.
!> \param auxiliary_images RI cells in the forward cell/atom ordering.
!> \param ao_images Right AO cells in covector ordering.
!> \param center_images Original reconstruction cells used for the forward tensor.
!> \param covectors Complete native energy derivatives in (RI,AO0,AOR) order.
!> \param atom_gradient Replicated positive Cartesian energy gradients.
!> \param strain_gradient Optional replicated affine-strain gradient, with atoms moving with the cell.
!> \param metric Optional native RI interaction, Coulomb by default.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_hfx_lattice_gradient(frame, kinds, augmentation, auxiliary_bases, auxiliary_images, &
                                                  ao_images, center_images, covectors, atom_gradient, strain_gradient, metric)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      TYPE(qs_kind_type), INTENT(IN)                     :: kinds(:)
      TYPE(uspp_augmentation_type), INTENT(IN), TARGET   :: augmentation(:)
      TYPE(gto_basis_set_p_type), INTENT(IN)             :: auxiliary_bases(:)
      INTEGER, INTENT(IN)                                :: auxiliary_images(:, :), ao_images(:, :), &
                                                            center_images(:, :)
      TYPE(dbt_type), INTENT(INOUT)                      :: covectors(:)
      REAL(dp), INTENT(OUT)                              :: atom_gradient(:, :)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: strain_gradient(3, 3)
      TYPE(coulomb_operator_type), INTENT(IN), OPTIONAL  :: metric

      INTEGER                                            :: atom, axis, column
      INTEGER, ALLOCATABLE                               :: active(:)
      REAL(dp)                                           :: cell_response(3, 3)
      REAL(dp), ALLOCATABLE                              :: atom_response(:, :)

      CPASSERT(frame%atom_local)
      CPASSERT(ALL(SHAPE(atom_gradient) == SHAPE(frame%positions)))
      CPASSERT(SIZE(covectors) > 0 .AND. SIZE(ao_images, 2) == SIZE(covectors))
      ALLOCATE (active(SIZE(covectors)))
      CALL active_lattice_covectors(covectors, active)
      ALLOCATE (atom_response, MOLD=frame%positions)
      cell_response = 0.0_dp
      DO atom = 1, SIZE(frame%positions, 2)
         DO axis = 1, 3
            atom_response = 0.0_dp
            atom_response(axis, atom) = 1.0_dp
            CALL lattice_gradient_direction(frame, kinds, augmentation, auxiliary_bases, auxiliary_images, &
                                            ao_images, center_images, covectors, active, atom_response, cell_response, &
                                            atom_gradient(axis, atom), metric)
         END DO
      END DO
      IF (PRESENT(strain_gradient)) THEN
         DO column = 1, 3
            DO axis = 1, 3
               cell_response = 0.0_dp
               cell_response(axis, :) = frame%cell(column, :)
               atom_response = 0.0_dp
               atom_response(axis, :) = frame%positions(column, :)
               CALL lattice_gradient_direction(frame, kinds, augmentation, auxiliary_bases, auxiliary_images, &
                                               ao_images, center_images, covectors, active, atom_response, cell_response, &
                                               strain_gradient(axis, column), metric)
            END DO
         END DO
      END IF
   END SUBROUTINE qs_upf_inverse_hfx_lattice_gradient

! **************************************************************************************************
!> \brief Select physical images with a nonzero energy covector on any native tensor owner.
!>        Exact zeros are omitted without introducing another screening threshold.
!> \param covectors ...
!> \param active Replicated selection in the original physical AO image order.
! **************************************************************************************************
   SUBROUTINE active_lattice_covectors(covectors, active)
      TYPE(dbt_type), INTENT(INOUT)                      :: covectors(:)
      INTEGER, INTENT(OUT)                               :: active(:)

      INTEGER                                            :: i, image, j, k, nblks(3), owner
      LOGICAL                                            :: found
      REAL(dp), ALLOCATABLE, TARGET                      :: values(:, :, :)

      CPASSERT(SIZE(active) == SIZE(covectors) .AND. SIZE(active) > 0)
      active = 0
      DO image = 1, SIZE(covectors)
         CALL matching_tensor(covectors(1), covectors(image))
         CALL dbt_get_info(covectors(image), nblks_total=nblks)
         DO i = 1, nblks(1)
         DO j = 1, nblks(2)
         DO k = 1, nblks(3)
            IF (active(image) /= 0) CYCLE
            CALL dbt_get_stored_coordinates(covectors(image), [i, j, k], owner)
            IF (owner /= covectors(image)%pgrid%mp_comm_2d%mepos) CYCLE
            CALL dbt_get_block(covectors(image), [i, j, k], values, found)
            IF (found) THEN
               IF (ANY(values /= 0.0_dp)) active(image) = 1
            END IF
            IF (ALLOCATED(values)) DEALLOCATE (values)
         END DO
         END DO
         END DO
      END DO
      ! Every rank must take the same image branches before collective integral evaluation.
      CALL covectors(1)%pgrid%mp_comm_2d%max(active)
   END SUBROUTINE active_lattice_covectors

! **************************************************************************************************
!> \brief Pair one complete analytic lattice response with the native physical integral covectors.
!> \param frame ...
!> \param kinds ...
!> \param augmentation ...
!> \param auxiliary_bases ...
!> \param auxiliary_images ...
!> \param ao_images ...
!> \param center_images ...
!> \param covectors ...
!> \param active ...
!> \param atom_response ...
!> \param cell_response ...
!> \param derivative ...
!> \param metric ...
! **************************************************************************************************
   SUBROUTINE lattice_gradient_direction(frame, kinds, augmentation, auxiliary_bases, auxiliary_images, &
                                         ao_images, center_images, covectors, active, atom_response, cell_response, derivative, &
                                         metric)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      TYPE(qs_kind_type), INTENT(IN)                     :: kinds(:)
      TYPE(uspp_augmentation_type), INTENT(IN), TARGET   :: augmentation(:)
      TYPE(gto_basis_set_p_type), INTENT(IN)             :: auxiliary_bases(:)
      INTEGER, INTENT(IN)                                :: auxiliary_images(:, :), ao_images(:, :), &
                                                            center_images(:, :)
      TYPE(dbt_type), INTENT(INOUT)                      :: covectors(:)
      INTEGER, INTENT(IN)                                :: active(:)
      REAL(dp), INTENT(IN)                               :: atom_response(:, :), cell_response(3, 3)
      REAL(dp), INTENT(OUT)                              :: derivative
      TYPE(coulomb_operator_type), INTENT(IN), OPTIONAL  :: metric

      INTEGER                                            :: image
      REAL(dp)                                           :: contribution
      TYPE(dbt_type)                                     :: integrals(1), responses(1)

      CALL dbt_create(covectors(1), integrals(1))
      CALL dbt_create(covectors(1), responses(1))
      derivative = 0.0_dp
      CPASSERT(SIZE(active) == SIZE(covectors))
      DO image = 1, SIZE(covectors)
         IF (active(image) == 0) CYCLE
         CALL matching_tensor(covectors(1), covectors(image))
         CALL dbt_clear(integrals(1))
         CALL dbt_clear(responses(1))
         CALL qs_upf_inverse_hfx_lattice_tensor(frame, kinds, augmentation, auxiliary_bases, auxiliary_images, &
                                                ao_images(:, image:image), center_images, integrals, responses, &
                                                atom_response, cell_response, metric)
         CALL qs_upf_inverse_hfx_contract(covectors(image), responses(1), contribution)
         derivative = derivative + contribution
      END DO
      CALL dbt_destroy(integrals(1))
      CALL dbt_destroy(responses(1))
   END SUBROUTINE lattice_gradient_direction

! **************************************************************************************************
!> \brief Contract the RI owner's energy covector with a distributed integral response.
!>        For a complex tensor the real energy variation is Re sum(conjg(G)*dC). Both native
!>        real components are paired without an extra factor of two. The owner supplies G with
!>        its occupations, exchange fraction, mesh weights and fitting convention. Explicit
!>        responses of the fitting metric remain the owner's separate contribution.
!>        The tensors must have the same full block partition and ownership. Missing blocks
!>        contribute zero. The returned scalar is summed once on their shared communicator.
!> \param covector_real Real component of the owner's energy derivative with respect to C.
!> \param response_real Real component of the original integral response dC.
!> \param contraction Replicated real energy variation from these tensors.
!> \param covector_imaginary Optional imaginary component, paired with response_imaginary.
!> \param response_imaginary Optional imaginary integral response.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_hfx_contract(covector_real, response_real, contraction, &
                                          covector_imaginary, response_imaginary)
      TYPE(dbt_type), INTENT(INOUT)                      :: covector_real, response_real
      REAL(dp), INTENT(OUT)                              :: contraction
      TYPE(dbt_type), INTENT(INOUT), OPTIONAL            :: covector_imaginary, response_imaginary

      CPASSERT(PRESENT(covector_imaginary) .EQV. PRESENT(response_imaginary))
      CALL contract_real_tensors(covector_real, response_real, contraction)
      IF (PRESENT(covector_imaginary)) THEN
         BLOCK
            REAL(dp) :: imaginary_contraction

            CALL matching_tensor(covector_real, covector_imaginary)
            CALL contract_real_tensors(covector_imaginary, response_imaginary, imaginary_contraction)
            contraction = contraction + imaginary_contraction
         END BLOCK
      END IF
      CALL covector_real%pgrid%mp_comm_2d%sum(contraction)
   END SUBROUTINE qs_upf_inverse_hfx_contract

! **************************************************************************************************
!> \brief Local Frobenius pairing on native tensor owners, without a collective reduction.
!> \param covector ...
!> \param response ...
!> \param contraction ...
! **************************************************************************************************
   SUBROUTINE contract_real_tensors(covector, response, contraction)
      TYPE(dbt_type), INTENT(INOUT)                      :: covector, response
      REAL(dp), INTENT(OUT)                              :: contraction

      INTEGER                                            :: i, j, k, nblks(3), owner
      LOGICAL                                            :: found_covector, found_response
      REAL(dp), ALLOCATABLE, TARGET                      :: dc(:, :, :), g(:, :, :)

      CALL matching_tensor(covector, response)
      CALL dbt_get_info(covector, nblks_total=nblks)
      contraction = 0.0_dp
      DO i = 1, nblks(1)
      DO j = 1, nblks(2)
      DO k = 1, nblks(3)
         CALL dbt_get_stored_coordinates(covector, [i, j, k], owner)
         IF (owner /= covector%pgrid%mp_comm_2d%mepos) CYCLE
         CALL dbt_get_block(covector, [i, j, k], g, found_covector)
         IF (.NOT. found_covector) CYCLE
         CALL dbt_get_block(response, [i, j, k], dc, found_response)
         IF (found_response) contraction = contraction + SUM(g*dc)
      END DO
      END DO
      END DO
   END SUBROUTINE contract_real_tensors

! **************************************************************************************************
!> \brief Add inverse corrections to the RI owner's full distributed (RI,left AO,right AO) tensors.
!>        The supplied frames must be replicated on the tensor communicator. Their defining
!>        coefficients and responses are prepared collectively before this routine. Integral
!>        evaluation is local to each tensor owner, including blocks absent from the Gaussian
!>        sparsity pattern. The native owner applies its combined filter and fitting metric later.
!>        RI images, factors and reconstruction images are supplied explicitly by that owner.
!>        Supplying ao_image selects individual real-space AO cells 0 and R instead of Bloch
!>        orbitals. This requires a common DIRECT or INVERSE_1C frame. Its hard overlaps are
!>        evaluated for each reconstruction image, before the native local fitting operation.
!>        The complex tensor is represented by two native real tensors with identical ownership.
!>        No exchange fraction, occupation, mesh weight, symmetry factor or RI inverse is applied.
!> \param left_frame ...
!> \param right_frame ...
!> \param kinds Original radial data and native per-kind angular controls.
!> \param augmentation Original Q caches, indexed by native kind. PAW transitions use full AE/PS products.
!> \param auxiliary_bases Native RI bases, one per supplied auxiliary site/image.
!> \param auxiliary_positions Cartesian RI site/image positions, indexed (axis,site).
!> \param auxiliary_factors Owner's complex RI conjugation/phase/normalization factors, by site.
!> \param center_images Prescribed original reconstruction images, indexed (axis,image).
!> \param tensor_real Native full tensor to increment, retaining its block partition.
!> \param tensor_imaginary Matching imaginary part, required for nonzero complex corrections.
!> \param response_real Optional native response tensor to increment.
!> \param response_imaginary Matching imaginary response, if the complex tensor is requested.
!> \param atom_response ...
!> \param cell_response ...
!> \param left_k_response ...
!> \param right_k_response ...
!> \param left_coefficients_response Full response from the defining left frame.
!> \param right_coefficients_response Full response from the defining right frame.
!> \param auxiliary_positions_response Cartesian RI site/image motion.
!> \param auxiliary_factors_response Response of the owner's complex RI factors.
!> \param ao_image Optional integer right AO cell R for individual lattice integrals.
!> \param metric Optional native RI interaction, Coulomb by default.
!> \param auxiliary_images Optional fixed auxiliary Gaussian images, at fixed Gamma.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_hfx_tensor(left_frame, right_frame, kinds, augmentation, auxiliary_bases, &
                                        auxiliary_positions, auxiliary_factors, center_images, tensor_real, &
                                        tensor_imaginary, response_real, response_imaginary, atom_response, cell_response, &
                                        left_k_response, right_k_response, left_coefficients_response, &
                                        right_coefficients_response, auxiliary_positions_response, auxiliary_factors_response, &
                                        ao_image, metric, auxiliary_images)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: left_frame, right_frame
      TYPE(qs_kind_type), INTENT(IN)                     :: kinds(:)
      TYPE(uspp_augmentation_type), INTENT(IN), TARGET   :: augmentation(:)
      TYPE(gto_basis_set_p_type), INTENT(IN)             :: auxiliary_bases(:)
      REAL(dp), INTENT(IN)                               :: auxiliary_positions(:, :)
      COMPLEX(dp), INTENT(IN)                            :: auxiliary_factors(:)
      INTEGER, INTENT(IN)                                :: center_images(:, :)
      TYPE(dbt_type), INTENT(INOUT)                      :: tensor_real
      TYPE(dbt_type), INTENT(INOUT), OPTIONAL            :: tensor_imaginary, response_real, &
                                                            response_imaginary
      REAL(dp), INTENT(IN), OPTIONAL                     :: atom_response(:, :), &
                                                            cell_response(3, 3), &
                                                            left_k_response(3), right_k_response(3)
      COMPLEX(dp), INTENT(IN), OPTIONAL :: left_coefficients_response(:, :), &
         right_coefficients_response(:, :)
      REAL(dp), INTENT(IN), OPTIONAL :: auxiliary_positions_response(:, :)
      COMPLEX(dp), INTENT(IN), OPTIONAL                  :: auxiliary_factors_response(:)
      INTEGER, INTENT(IN), OPTIONAL                      :: ao_image(3)
      TYPE(coulomb_operator_type), INTENT(IN), OPTIONAL  :: metric
      INTEGER, INTENT(IN), OPTIONAL                      :: auxiliary_images(:, :)

      COMPLEX(dp), ALLOCATABLE                           :: correction(:, :, :), dc(:, :, :), &
                                                            dn(:, :, :), native(:, :, :)
      INTEGER :: block_i, block_j, block_k, center, first, image, inner, iset, kind, last, maxl, &
         nblks(3), ninner, nl, nouter, nr, outer, owner, rule, site, slot, start, STOP
      INTEGER, ALLOCATABLE                               :: auxiliary_start(:), sizes_i(:), &
                                                            sizes_j(:), sizes_k(:), start_i(:), &
                                                            start_j(:), start_k(:)
      LOGICAL                                            :: derivatives
      REAL(dp), ALLOCATABLE                              :: radial_weights(:)
      TYPE(gaussian_coulomb_fields_type)                 :: fields
      TYPE(lattice_overlap_type), ALLOCATABLE            :: overlaps(:, :)
      TYPE(transition_augmentation_type), ALLOCATABLE, &
         TARGET                                          :: qcache(:, :, :)
      TYPE(transition_augmentation_type), POINTER        :: qentry
      TYPE(uspp_augmentation_type), ALLOCATABLE, TARGET  :: paw_augmentation(:)
      TYPE(uspp_augmentation_type), POINTER              :: transition_augmentation

      CPASSERT(ALL(SHAPE(auxiliary_positions) == [3, SIZE(auxiliary_bases)]))
      CPASSERT(SIZE(auxiliary_factors) == SIZE(auxiliary_bases))
      CPASSERT(SIZE(center_images, 1) == 3 .AND. SIZE(center_images, 2) > 0)
      CPASSERT(SIZE(augmentation) == SIZE(kinds))
      IF (PRESENT(ao_image)) THEN
         CPASSERT(.NOT. PRESENT(auxiliary_images))
         CPASSERT(left_frame%atom_local .AND. right_frame%atom_local)
         CPASSERT(left_frame%direct .EQV. right_frame%direct)
         CPASSERT(ALL(left_frame%cell == right_frame%cell))
         CPASSERT(ALL(left_frame%positions == right_frame%positions))
         CPASSERT(ALL(left_frame%periodicity == right_frame%periodicity))
         CPASSERT(ALL(left_frame%ao_first == right_frame%ao_first))
         CPASSERT(ALL(left_frame%ao_last == right_frame%ao_last))
         CPASSERT(ALL(left_frame%center_atoms == right_frame%center_atoms))
         CPASSERT(ALL(left_frame%center_kinds == right_frame%center_kinds))
         CPASSERT(SIZE(left_frame%bases) == SIZE(right_frame%bases))
         DO site = 1, SIZE(left_frame%bases)
            IF (.NOT. (ASSOCIATED(left_frame%bases(site)%gto_basis_set, right_frame%bases(site)%gto_basis_set))) THEN
               CPABORT("Transition bases must refer to the same object")
            END IF
         END DO
         CPASSERT(.NOT. PRESENT(left_k_response) .AND. .NOT. PRESENT(right_k_response))
         IF (.NOT. (.NOT. PRESENT(left_coefficients_response) .AND. .NOT. PRESENT(right_coefficients_response))) THEN
            CPABORT("Coefficient responses require a varying transition")
         END IF
      END IF
      CALL prepare_paw_transitions(kinds, left_frame, paw_augmentation)
      CALL tensor_partition(left_frame, right_frame, auxiliary_bases, tensor_real, nblks, sizes_i, sizes_j, sizes_k, &
                            start_i, start_j, start_k, auxiliary_start, maxl)
      derivatives = PRESENT(response_real)
      IF (PRESENT(tensor_imaginary)) CALL matching_tensor(tensor_real, tensor_imaginary)
      IF (derivatives) THEN
         CALL matching_tensor(tensor_real, response_real)
         CPASSERT(PRESENT(response_imaginary) .EQV. PRESENT(tensor_imaginary))
         IF (PRESENT(response_imaginary)) CALL matching_tensor(tensor_real, response_imaginary)
         CPASSERT(PRESENT(atom_response) .AND. PRESENT(cell_response))
         IF (.NOT. PRESENT(ao_image)) THEN
            CPASSERT(PRESENT(left_k_response) .AND. PRESENT(right_k_response))
            IF (.NOT. (PRESENT(left_coefficients_response) .AND. PRESENT(right_coefficients_response))) THEN
               CPABORT("Both coefficient responses are required")
            END IF
         END IF
         IF (.NOT. (PRESENT(auxiliary_positions_response) .AND. PRESENT(auxiliary_factors_response))) THEN
            CPABORT("Both auxiliary responses are required")
         END IF
         CPASSERT(ALL(SHAPE(auxiliary_positions_response) == SHAPE(auxiliary_positions)))
         CPASSERT(SIZE(auxiliary_factors_response) == SIZE(auxiliary_factors))
      ELSE
         CPASSERT(.NOT. PRESENT(response_imaginary))
      END IF
      CALL init_lebedev_grids()
      CALL gaussian_coulomb_fields_init(maxl)
      nouter = nblks(2)*nblks(3)
      ninner = nblks(1)
      IF (.NOT. PRESENT(ao_image)) THEN
         nouter = nblks(1)
         ninner = nblks(2)*nblks(3)
      END IF
      DO outer = 1, nouter
         IF (PRESENT(ao_image)) THEN
            block_j = (outer - 1)/nblks(3) + 1
            block_k = MOD(outer - 1, nblks(3)) + 1
            ALLOCATE (overlaps(SIZE(left_frame%center_atoms), SIZE(center_images, 2)))
         ELSE
            block_i = outer
            ALLOCATE (qcache(SIZE(left_frame%center_atoms), SIZE(center_images, 2), sizes_i(block_i)))
         END IF
         DO inner = 1, ninner
            IF (PRESENT(ao_image)) THEN
               block_i = inner
            ELSE
               block_j = (inner - 1)/nblks(3) + 1
               block_k = MOD(inner - 1, nblks(3)) + 1
            END IF
            CALL dbt_get_stored_coordinates(tensor_real, [block_i, block_j, block_k], owner)
            IF (owner /= tensor_real%pgrid%mp_comm_2d%mepos) CYCLE
            nl = sizes_j(block_j)
            nr = sizes_k(block_k)
            ALLOCATE (native(sizes_i(block_i), nl, nr), SOURCE=CMPLX(0.0_dp, 0.0_dp, dp))
            IF (derivatives) ALLOCATE (dn(sizes_i(block_i), nl, nr), SOURCE=CMPLX(0.0_dp, 0.0_dp, dp))
            DO site = 1, SIZE(auxiliary_bases)
               IF (auxiliary_start(site + 1) <= start_i(block_i) .OR. &
                   auxiliary_start(site) >= start_i(block_i + 1)) CYCLE
               first = auxiliary_start(site)
               DO iset = 1, auxiliary_bases(site)%gto_basis_set%nset
                  last = first + auxiliary_bases(site)%gto_basis_set%nsgf_set(iset) - 1
                  start = MAX(first, start_i(block_i))
                  STOP = MIN(last, start_i(block_i + 1) - 1)
                  IF (start <= STOP) THEN
                     CALL prepare_auxiliary_fields(auxiliary_bases(site)%gto_basis_set, iset, fields, metric)
                     ALLOCATE (correction(nl, nr, fields%nfunction))
                     IF (derivatives) ALLOCATE (dc(nl, nr, fields%nfunction))
                     DO center = 1, SIZE(left_frame%center_atoms)
                        kind = left_frame%center_kinds(center)
                        CPASSERT(ASSOCIATED(kinds(kind)%sgp_potential))
                        CPASSERT(ALLOCATED(kinds(kind)%sgp_potential%upf_data))
                        transition_augmentation => augmentation(kind)
                        IF (paw_augmentation(kind)%nproj > 0) transition_augmentation => paw_augmentation(kind)
                        ASSOCIATE (pot => kinds(kind)%sgp_potential%upf_data)
                           ALLOCATE (radial_weights(SIZE(pot%r)))
                           radial_weights(:) = upf_radial_weights(pot)*pot%r**2
                           rule = MAX(get_number_of_lebedev_grid(n=kinds(kind)%ngrid_ang), &
                                      get_number_of_lebedev_grid(l=MAX(0, transition_augmentation%lmax) + maxl + 1))
                           DO image = 1, SIZE(center_images, 2)
                              NULLIFY (qentry)
                              IF (ALLOCATED(qcache)) THEN
                                 slot = start - start_i(block_i) + 1
                                 qentry => qcache(center, image, slot)
                              END IF
                              IF (derivatives) THEN
                                 IF (PRESENT(ao_image)) THEN
                                    CALL lattice_center_with_overlaps(left_frame, center, &
                                                                      [start_j(block_j), start_j(block_j + 1) - 1], &
                                                                      [start_k(block_k), start_k(block_k + 1) - 1], &
                                                                      ao_image, center_images(:, image), fields, &
                                                                      auxiliary_positions(:, site), transition_augmentation, &
                                                                      pot%r, radial_weights, lebedev_grid(rule)%r, &
                                                                      fourpi*lebedev_grid(rule)%w, correction, dc, &
                                                                      atom_response, cell_response, &
                                                                      auxiliary_positions_response(:, site), &
                                                                      overlaps=overlaps(center, image))
                                 ELSE
                                    CALL qs_upf_inverse_hfx_center(left_frame, right_frame, center, &
                                                                   [start_j(block_j), start_j(block_j + 1) - 1], &
                                                                   [start_k(block_k), start_k(block_k + 1) - 1], fields, &
                                                                   auxiliary_positions(:, site), transition_augmentation, &
                                                                   pot%r, radial_weights, lebedev_grid(rule)%r, &
                                                                   fourpi*lebedev_grid(rule)%w, correction, dc, &
                                                                  atom_response, cell_response, left_k_response, right_k_response, &
                                                                   left_coefficients_response, right_coefficients_response, &
                                                                   auxiliary_positions_response(:, site), center_images(:, image), &
                                                                   auxiliary_images, qentry)
                                 END IF
                              ELSE
                                 IF (PRESENT(ao_image)) THEN
                                    CALL lattice_center_with_overlaps(left_frame, center, &
                                                                      [start_j(block_j), start_j(block_j + 1) - 1], &
                                                                      [start_k(block_k), start_k(block_k + 1) - 1], &
                                                                      ao_image, center_images(:, image), fields, &
                                                                      auxiliary_positions(:, site), transition_augmentation, &
                                                                      pot%r, radial_weights, lebedev_grid(rule)%r, &
                                                                      fourpi*lebedev_grid(rule)%w, correction, &
                                                                      overlaps=overlaps(center, image))
                                 ELSE
                                    CALL qs_upf_inverse_hfx_center(left_frame, right_frame, center, &
                                                                   [start_j(block_j), start_j(block_j + 1) - 1], &
                                                                   [start_k(block_k), start_k(block_k + 1) - 1], fields, &
                                                                   auxiliary_positions(:, site), transition_augmentation, &
                                                                   pot%r, radial_weights, lebedev_grid(rule)%r, &
                                                                   fourpi*lebedev_grid(rule)%w, correction, &
                                                                   image_shift=center_images(:, image), &
                                                                   auxiliary_images=auxiliary_images, augmentation_integrals=qentry)
                                 END IF
                              END IF
                              native(start - start_i(block_i) + 1:STOP - start_i(block_i) + 1, :, :) = &
                                 native(start - start_i(block_i) + 1:STOP - start_i(block_i) + 1, :, :) + &
                                 auxiliary_factors(site)*RESHAPE(correction(:, :, start - first + 1:STOP - first + 1), &
                                                                 [STOP - start + 1, nl, nr], ORDER=[2, 3, 1])
                              IF (derivatives) THEN
                                 dn(start - start_i(block_i) + 1:STOP - start_i(block_i) + 1, :, :) = &
                                    dn(start - start_i(block_i) + 1:STOP - start_i(block_i) + 1, :, :) + &
                                    RESHAPE(auxiliary_factors(site)*dc(:, :, start - first + 1:STOP - first + 1) + &
                                            auxiliary_factors_response(site)*correction(:, :, start - first + 1:STOP - first + 1), &
                                            [STOP - start + 1, nl, nr], ORDER=[2, 3, 1])
                              END IF
                           END DO
                           DEALLOCATE (radial_weights)
                        END ASSOCIATE
                     END DO
                     DEALLOCATE (correction)
                     IF (derivatives) DEALLOCATE (dc)
                  END IF
                  first = last + 1
               END DO
            END DO
            CALL dbt_put_block(tensor_real, [block_i, block_j, block_k], SHAPE(native), REAL(native, dp), summation=.TRUE.)
            IF (PRESENT(tensor_imaginary)) THEN
               CALL dbt_put_block(tensor_imaginary, [block_i, block_j, block_k], SHAPE(native), AIMAG(native), summation=.TRUE.)
            ELSE
               CPASSERT(ALL(AIMAG(native) == 0.0_dp))
            END IF
            IF (derivatives) THEN
               CALL dbt_put_block(response_real, [block_i, block_j, block_k], SHAPE(dn), REAL(dn, dp), summation=.TRUE.)
               IF (PRESENT(response_imaginary)) THEN
                  CALL dbt_put_block(response_imaginary, [block_i, block_j, block_k], SHAPE(dn), AIMAG(dn), summation=.TRUE.)
               ELSE
                  CPASSERT(ALL(AIMAG(dn) == 0.0_dp))
               END IF
               DEALLOCATE (dn)
            END IF
            DEALLOCATE (native)
         END DO
         IF (ALLOCATED(overlaps)) DEALLOCATE (overlaps)
         IF (ALLOCATED(qcache)) DEALLOCATE (qcache)
      END DO
      CALL dbt_finalize(tensor_real)
      IF (PRESENT(tensor_imaginary)) CALL dbt_finalize(tensor_imaginary)
      IF (derivatives) CALL dbt_finalize(response_real)
      IF (PRESENT(response_imaginary)) CALL dbt_finalize(response_imaginary)
   END SUBROUTINE qs_upf_inverse_hfx_tensor

! **************************************************************************************************
!> \brief Return full coefficient covectors from the native RI owner's distributed energy derivative.
!>        Only stored covector blocks contribute. Both complex components retain their original
!>        ownership and AO/RI ordering. The conjugate auxiliary image factor returns the local
!>        integral covector. Local center/image covectors are summed once on the tensor communicator.
!>        Every rank participates, including ranks without blocks. Frames and output coefficients
!>        are replicated on that communicator. The defining inverse solve follows this reduction.
!>        The optional geometry contraction holds W fixed. The caller adds the defining-projection
!>        response from the returned coefficient covectors. It can reuse that inverse adjoint for
!>        different directions. Explicit auxiliary image-factor responses are included even when
!>        the factor itself is zero. The native fitting-metric derivative remains separate.
!>        Omit coefficient outputs when reusing their adjoint for additional geometry directions.
!> \param left_frame ...
!> \param right_frame ...
!> \param kinds ...
!> \param augmentation ...
!> \param auxiliary_bases ...
!> \param auxiliary_positions ...
!> \param auxiliary_factors ...
!> \param center_images ...
!> \param covector_real Native real energy derivative, ordered (RI,left AO,right AO).
!> \param left_coefficients_covector Optional full replicated derivative with respect to the left W.
!> \param right_coefficients_covector Optional full replicated derivative with respect to the right W.
!> \param covector_imaginary Optional matching imaginary energy derivative, otherwise zero.
!> \param geometry_response Optional replicated explicit geometry/Bloch derivative at fixed W.
!> \param atom_response Cartesian motion of every native atom and Gaussian center.
!> \param cell_response Cartesian lattice response.
!> \param left_k_response Left reduced Bloch-vector response.
!> \param right_k_response Right reduced Bloch-vector response.
!> \param auxiliary_positions_response Motion of every supplied RI site/image.
!> \param auxiliary_factors_response Response of the owner's complex RI image factors.
!> \param metric Optional native RI interaction, Coulomb by default.
!> \param auxiliary_images Optional fixed auxiliary Gaussian images, at fixed Gamma.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_hfx_tensor_adjoint(left_frame, right_frame, kinds, augmentation, auxiliary_bases, &
                                                auxiliary_positions, auxiliary_factors, center_images, covector_real, &
                                                left_coefficients_covector, right_coefficients_covector, covector_imaginary, &
                                                geometry_response, atom_response, cell_response, left_k_response, &
                                                right_k_response, auxiliary_positions_response, auxiliary_factors_response, &
                                                metric, auxiliary_images)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: left_frame, right_frame
      TYPE(qs_kind_type), INTENT(IN)                     :: kinds(:)
      TYPE(uspp_augmentation_type), INTENT(IN), TARGET   :: augmentation(:)
      TYPE(gto_basis_set_p_type), INTENT(IN)             :: auxiliary_bases(:)
      REAL(dp), INTENT(IN)                               :: auxiliary_positions(:, :)
      COMPLEX(dp), INTENT(IN)                            :: auxiliary_factors(:)
      INTEGER, INTENT(IN)                                :: center_images(:, :)
      TYPE(dbt_type), INTENT(INOUT)                      :: covector_real
      COMPLEX(dp), INTENT(OUT), OPTIONAL :: left_coefficients_covector(:, :), &
         right_coefficients_covector(:, :)
      TYPE(dbt_type), INTENT(INOUT), OPTIONAL            :: covector_imaginary
      REAL(dp), INTENT(OUT), OPTIONAL                    :: geometry_response
      REAL(dp), INTENT(IN), OPTIONAL :: atom_response(:, :), cell_response(3, 3), &
         left_k_response(3), right_k_response(3), auxiliary_positions_response(:, :)
      COMPLEX(dp), INTENT(IN), OPTIONAL                  :: auxiliary_factors_response(:)
      TYPE(coulomb_operator_type), INTENT(IN), OPTIONAL  :: metric
      INTEGER, INTENT(IN), OPTIONAL                      :: auxiliary_images(:, :)

      COMPLEX(dp), ALLOCATABLE :: correction(:, :, :), dc(:, :, :), factor_covector(:, :, :), &
         g(:, :, :), gl(:, :), gr(:, :), local(:, :, :), zero_left(:, :), zero_right(:, :)
      INTEGER                                            :: block_i, block_j, block_k, center, &
                                                            finish, first, image, iset, kind, &
                                                            last, maxl, nblks(3), nl, nr, owner, &
                                                            rule, site, slot, start
      INTEGER, ALLOCATABLE                               :: auxiliary_start(:), sizes_i(:), &
                                                            sizes_j(:), sizes_k(:), start_i(:), &
                                                            start_j(:), start_k(:)
      LOGICAL                                            :: coefficients, derivatives, &
                                                            found_imaginary, found_real
      REAL(dp), ALLOCATABLE                              :: radial_weights(:)
      REAL(dp), ALLOCATABLE, TARGET                      :: imaginary_block(:, :, :), &
                                                            real_block(:, :, :)
      TYPE(gaussian_coulomb_fields_type)                 :: fields
      TYPE(transition_augmentation_type), ALLOCATABLE, &
         TARGET                                          :: qcache(:, :, :)
      TYPE(uspp_augmentation_type), ALLOCATABLE, TARGET  :: paw_augmentation(:)
      TYPE(uspp_augmentation_type), POINTER              :: transition_augmentation

      CPASSERT(ALL(SHAPE(auxiliary_positions) == [3, SIZE(auxiliary_bases)]))
      CPASSERT(SIZE(auxiliary_factors) == SIZE(auxiliary_bases))
      CPASSERT(SIZE(center_images, 1) == 3 .AND. SIZE(center_images, 2) > 0)
      CPASSERT(SIZE(augmentation) == SIZE(kinds))
      CALL prepare_paw_transitions(kinds, left_frame, paw_augmentation)
      IF (PRESENT(left_coefficients_covector)) THEN
         CPASSERT(ALL(SHAPE(left_coefficients_covector) == SHAPE(left_frame%coefficients)))
         left_coefficients_covector = CMPLX(0.0_dp, 0.0_dp, dp)
      END IF
      IF (PRESENT(right_coefficients_covector)) THEN
         CPASSERT(ALL(SHAPE(right_coefficients_covector) == SHAPE(right_frame%coefficients)))
         right_coefficients_covector = CMPLX(0.0_dp, 0.0_dp, dp)
      END IF
      CALL tensor_partition(left_frame, right_frame, auxiliary_bases, covector_real, nblks, sizes_i, sizes_j, sizes_k, &
                            start_i, start_j, start_k, auxiliary_start, maxl)
      IF (PRESENT(covector_imaginary)) CALL matching_tensor(covector_real, covector_imaginary)
      coefficients = PRESENT(left_coefficients_covector) .OR. PRESENT(right_coefficients_covector)
      derivatives = PRESENT(geometry_response)
      CPASSERT(coefficients .OR. derivatives)
      IF (derivatives) THEN
         CPASSERT(PRESENT(atom_response) .AND. PRESENT(cell_response))
         CPASSERT(PRESENT(left_k_response) .AND. PRESENT(right_k_response))
         IF (.NOT. (PRESENT(auxiliary_positions_response) .AND. PRESENT(auxiliary_factors_response))) THEN
            CPABORT("Both auxiliary responses are required")
         END IF
         CPASSERT(ALL(SHAPE(atom_response) == SHAPE(left_frame%positions)))
         CPASSERT(ALL(SHAPE(auxiliary_positions_response) == SHAPE(auxiliary_positions)))
         CPASSERT(SIZE(auxiliary_factors_response) == SIZE(auxiliary_factors))
         ALLOCATE (zero_left, MOLD=left_frame%coefficients)
         ALLOCATE (zero_right, MOLD=right_frame%coefficients)
         zero_left = CMPLX(0.0_dp, 0.0_dp, dp)
         zero_right = CMPLX(0.0_dp, 0.0_dp, dp)
         geometry_response = 0.0_dp
      END IF
      CALL init_lebedev_grids()
      CALL gaussian_coulomb_fields_init(maxl)
      DO block_i = 1, nblks(1)
         ALLOCATE (qcache(SIZE(left_frame%center_atoms), SIZE(center_images, 2), sizes_i(block_i)))
         DO block_j = 1, nblks(2)
         DO block_k = 1, nblks(3)
            CALL dbt_get_stored_coordinates(covector_real, [block_i, block_j, block_k], owner)
            IF (owner /= covector_real%pgrid%mp_comm_2d%mepos) CYCLE
            CALL dbt_get_block(covector_real, [block_i, block_j, block_k], real_block, found_real)
            found_imaginary = .FALSE.
            IF (PRESENT(covector_imaginary)) THEN
               CALL dbt_get_block(covector_imaginary, [block_i, block_j, block_k], imaginary_block, found_imaginary)
            END IF
            IF (.NOT. (found_real .OR. found_imaginary)) CYCLE
            nl = sizes_j(block_j)
            nr = sizes_k(block_k)
            ALLOCATE (g(sizes_i(block_i), nl, nr), SOURCE=CMPLX(0.0_dp, 0.0_dp, dp))
            IF (found_real) g(:, :, :) = CMPLX(real_block, 0.0_dp, dp)
            IF (found_imaginary) g(:, :, :) = g + CMPLX(0.0_dp, imaginary_block, dp)
            IF (coefficients) THEN
               ALLOCATE (gl(SIZE(left_frame%coefficients, 1), nl), gr(SIZE(right_frame%coefficients, 1), nr))
            END IF
            DO site = 1, SIZE(auxiliary_bases)
               IF (auxiliary_start(site + 1) <= start_i(block_i) .OR. &
                   auxiliary_start(site) >= start_i(block_i + 1)) CYCLE
               first = auxiliary_start(site)
               DO iset = 1, auxiliary_bases(site)%gto_basis_set%nset
                  last = first + auxiliary_bases(site)%gto_basis_set%nsgf_set(iset) - 1
                  start = MAX(first, start_i(block_i))
                  finish = MIN(last, start_i(block_i + 1) - 1)
                  IF (start <= finish) THEN
                     CALL prepare_auxiliary_fields(auxiliary_bases(site)%gto_basis_set, iset, fields, metric)
                     ALLOCATE (local(nl, nr, fields%nfunction), SOURCE=CMPLX(0.0_dp, 0.0_dp, dp))
                     local(:, :, start - first + 1:finish - first + 1) = CONJG(auxiliary_factors(site))* &
                                                      RESHAPE(g(start - start_i(block_i) + 1:finish - start_i(block_i) + 1, :, :), &
                                                                                 [nl, nr, finish - start + 1], ORDER=[3, 1, 2])
                     IF (derivatives) THEN
                        ALLOCATE (factor_covector(nl, nr, fields%nfunction), SOURCE=CMPLX(0.0_dp, 0.0_dp, dp))
                        ALLOCATE (correction(nl, nr, fields%nfunction), dc(nl, nr, fields%nfunction))
                        factor_covector(:, :, start - first + 1:finish - first + 1) = CONJG(auxiliary_factors_response(site))* &
                                                      RESHAPE(g(start - start_i(block_i) + 1:finish - start_i(block_i) + 1, :, :), &
                                                                                      [nl, nr, finish - start + 1], ORDER=[3, 1, 2])
                     END IF
                     DO center = 1, SIZE(left_frame%center_atoms)
                        kind = left_frame%center_kinds(center)
                        CPASSERT(ASSOCIATED(kinds(kind)%sgp_potential))
                        CPASSERT(ALLOCATED(kinds(kind)%sgp_potential%upf_data))
                        transition_augmentation => augmentation(kind)
                        IF (paw_augmentation(kind)%nproj > 0) transition_augmentation => paw_augmentation(kind)
                        ASSOCIATE (pot => kinds(kind)%sgp_potential%upf_data)
                           ALLOCATE (radial_weights(SIZE(pot%r)))
                           radial_weights(:) = upf_radial_weights(pot)*pot%r**2
                           rule = MAX(get_number_of_lebedev_grid(n=kinds(kind)%ngrid_ang), &
                                      get_number_of_lebedev_grid(l=MAX(0, transition_augmentation%lmax) + maxl + 1))
                           DO image = 1, SIZE(center_images, 2)
                              slot = start - start_i(block_i) + 1
                              IF (coefficients) THEN
                                 CALL qs_upf_inverse_hfx_center_adjoint(left_frame, right_frame, center, &
                                                                        [start_j(block_j), start_j(block_j + 1) - 1], &
                                                                        [start_k(block_k), start_k(block_k + 1) - 1], fields, &
                                                                     auxiliary_positions(:, site), transition_augmentation, pot%r, &
                                                                        radial_weights, lebedev_grid(rule)%r, &
                                                                        fourpi*lebedev_grid(rule)%w, local, gl, gr, &
                                                                        center_images(:, image), auxiliary_images, &
                                                                        qcache(center, image, slot))
                              END IF
                              IF (derivatives) THEN
                                 CALL qs_upf_inverse_hfx_center(left_frame, right_frame, center, &
                                                                [start_j(block_j), start_j(block_j + 1) - 1], &
                                                                [start_k(block_k), start_k(block_k + 1) - 1], fields, &
                                                                auxiliary_positions(:, site), transition_augmentation, pot%r, &
                                                                radial_weights, lebedev_grid(rule)%r, fourpi*lebedev_grid(rule)%w, &
                                                                correction, dc, atom_response, cell_response, left_k_response, &
                                                                right_k_response, zero_left, zero_right, &
                                                                auxiliary_positions_response(:, site), center_images(:, image), &
                                                                auxiliary_images, qcache(center, image, slot))
                                 geometry_response = geometry_response + REAL(SUM(CONJG(local)*dc) + &
                                                                              SUM(CONJG(factor_covector)*correction), dp)
                              END IF
                              IF (PRESENT(left_coefficients_covector)) THEN
                                 left_coefficients_covector(:, start_j(block_j):start_j(block_j + 1) - 1) = &
                                    left_coefficients_covector(:, start_j(block_j):start_j(block_j + 1) - 1) + gl
                              END IF
                              IF (PRESENT(right_coefficients_covector)) THEN
                                 right_coefficients_covector(:, start_k(block_k):start_k(block_k + 1) - 1) = &
                                    right_coefficients_covector(:, start_k(block_k):start_k(block_k + 1) - 1) + gr
                              END IF
                           END DO
                           DEALLOCATE (radial_weights)
                        END ASSOCIATE
                     END DO
                     DEALLOCATE (local)
                     IF (derivatives) DEALLOCATE (factor_covector, correction, dc)
                  END IF
                  first = last + 1
               END DO
            END DO
            DEALLOCATE (g)
            IF (coefficients) DEALLOCATE (gl, gr)
         END DO
         END DO
         DEALLOCATE (qcache)
      END DO
      IF (PRESENT(left_coefficients_covector)) CALL covector_real%pgrid%mp_comm_2d%sum(left_coefficients_covector)
      IF (PRESENT(right_coefficients_covector)) CALL covector_real%pgrid%mp_comm_2d%sum(right_coefficients_covector)
      IF (derivatives) CALL covector_real%pgrid%mp_comm_2d%sum(geometry_response)
   END SUBROUTINE qs_upf_inverse_hfx_tensor_adjoint

! **************************************************************************************************
!> \brief Prepare the same auxiliary interaction as the native RI integral owner.
!> \param basis ...
!> \param iset ...
!> \param fields ...
!> \param metric Optional native fitting interaction.
! **************************************************************************************************
   SUBROUTINE prepare_auxiliary_fields(basis, iset, fields, metric)
      TYPE(gto_basis_set_type), INTENT(IN)               :: basis
      INTEGER, INTENT(IN)                                :: iset
      TYPE(gaussian_coulomb_fields_type), INTENT(OUT)    :: fields
      TYPE(coulomb_operator_type), INTENT(IN), OPTIONAL  :: metric

      CALL gaussian_coulomb_fields_prepare(basis, iset, fields)
      IF (.NOT. PRESENT(metric)) RETURN
      SELECT CASE (metric%potential_type)
      CASE (do_potential_coulomb)
      CASE (do_potential_long, do_potential_short)
         CALL gaussian_coulomb_fields_screen(fields, metric%omega, metric%potential_type == do_potential_short)
      CASE DEFAULT
         CPABORT("Inverse UPF fitting interaction is not implemented yet")
      END SELECT
   END SUBROUTINE prepare_auxiliary_fields

! **************************************************************************************************
!> \brief Original augmentation integrals and gradients at fixed geometry, images and RI shell.
!> \param fields ...
!> \param augmentation ...
!> \param displacement ...
!> \param directions ...
!> \param angular_weights ...
!> \param integrals ...
!> \param derivatives ...
!> \param cell ...
!> \param images ...
! **************************************************************************************************
   SUBROUTINE prepare_transition_augmentation(fields, augmentation, displacement, directions, angular_weights, &
                                              integrals, derivatives, cell, images)
      TYPE(gaussian_coulomb_fields_type), INTENT(IN)     :: fields
      TYPE(uspp_augmentation_type), INTENT(IN)           :: augmentation
      REAL(dp), INTENT(IN)                               :: displacement(3), directions(:, :), &
                                                            angular_weights(:)
      TYPE(transition_augmentation_type), INTENT(INOUT)  :: integrals
      LOGICAL, INTENT(IN)                                :: derivatives
      REAL(dp), INTENT(IN)                               :: cell(3, 3)
      INTEGER, INTENT(IN), OPTIONAL                      :: images(:, :)

      INTEGER                                            :: n, naux

      IF (ALLOCATED(integrals%values)) THEN
         IF (.NOT. derivatives .OR. ALLOCATED(integrals%gradient)) RETURN
      ELSE
         ALLOCATE (integrals%values(augmentation%nproj, augmentation%nproj, fields%nfunction))
      END IF
      n = augmentation%nproj
      naux = fields%nfunction
      IF (derivatives) THEN
         ALLOCATE (integrals%gradient(n, n, naux, 3))
         IF (PRESENT(images)) THEN
            ALLOCATE (integrals%cell_gradient(n, n, naux, 3, 3))
            CALL upf_inverse_exchange_augmentation(fields, augmentation, displacement, directions, angular_weights, &
                                                   integrals%values, augmentation_derivative=integrals%gradient, &
                                                   auxiliary_cell=cell, auxiliary_images=images, &
                                                   cell_derivative=integrals%cell_gradient)
         ELSE
            CALL upf_inverse_exchange_augmentation(fields, augmentation, displacement, directions, angular_weights, &
                                                   integrals%values, augmentation_derivative=integrals%gradient)
         END IF
      ELSE
         IF (PRESENT(images)) THEN
            CALL upf_inverse_exchange_augmentation(fields, augmentation, displacement, directions, angular_weights, &
                                                   integrals%values, auxiliary_cell=cell, auxiliary_images=images)
         ELSE
            CALL upf_inverse_exchange_augmentation(fields, augmentation, displacement, directions, angular_weights, &
                                                   integrals%values)
         END IF
      END IF
   END SUBROUTINE prepare_transition_augmentation

! **************************************************************************************************
!> \brief Prepare PAW valence transition kernels without changing original compensation data.
!> \param kinds Original datasets and prepared projectors.
!> \param frame Selected reconstruction centers.
!> \param augmentation Per-kind AE-minus-PS products. Non-PAW entries remain empty.
! **************************************************************************************************
   SUBROUTINE prepare_paw_transitions(kinds, frame, augmentation)
      TYPE(qs_kind_type), INTENT(IN)                     :: kinds(:)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      TYPE(uspp_augmentation_type), ALLOCATABLE, &
         INTENT(OUT)                                     :: augmentation(:)

      INTEGER                                            :: center, kind

      ALLOCATE (augmentation(SIZE(kinds)))
      DO center = 1, SIZE(frame%center_kinds)
         kind = frame%center_kinds(center)
         CPASSERT(ASSOCIATED(kinds(kind)%sgp_potential))
         CPASSERT(ALLOCATED(kinds(kind)%sgp_potential%upf_data))
         IF (.NOT. kinds(kind)%sgp_potential%upf_data%is_paw) CYCLE
         IF (augmentation(kind)%nproj > 0) CYCLE
         CALL paw_valence_augmentation_init(kinds(kind)%sgp_potential%upf_data, &
                                            kinds(kind)%sgp_potential%upf_projectors, augmentation(kind))
      END DO
   END SUBROUTINE prepare_paw_transitions

! **************************************************************************************************
!> \brief Common original RI-site and AO partitions for the forward tensor and its adjoint.
!> \param left_frame ...
!> \param right_frame ...
!> \param auxiliary_bases ...
!> \param tensor ...
!> \param nblks ...
!> \param sizes_i ...
!> \param sizes_j ...
!> \param sizes_k ...
!> \param start_i ...
!> \param start_j ...
!> \param start_k ...
!> \param auxiliary_start ...
!> \param maxl ...
! **************************************************************************************************
   SUBROUTINE tensor_partition(left_frame, right_frame, auxiliary_bases, tensor, nblks, sizes_i, sizes_j, sizes_k, &
                               start_i, start_j, start_k, auxiliary_start, maxl)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: left_frame, right_frame
      TYPE(gto_basis_set_p_type), INTENT(IN)             :: auxiliary_bases(:)
      TYPE(dbt_type), INTENT(IN)                         :: tensor
      INTEGER, INTENT(OUT)                               :: nblks(3)
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: sizes_i(:), sizes_j(:), sizes_k(:), &
                                                            start_i(:), start_j(:), start_k(:), &
                                                            auxiliary_start(:)
      INTEGER, INTENT(OUT)                               :: maxl

      INTEGER                                            :: i, j, k, site

      CALL dbt_get_info(tensor, nblks_total=nblks)
      ALLOCATE (sizes_i(nblks(1)), sizes_j(nblks(2)), sizes_k(nblks(3)), &
                start_i(nblks(1) + 1), start_j(nblks(2) + 1), start_k(nblks(3) + 1), &
                auxiliary_start(SIZE(auxiliary_bases) + 1))
      CALL dbt_get_info(tensor, blk_size_1=sizes_i, blk_size_2=sizes_j, blk_size_3=sizes_k)
      start_i(1) = 1
      start_j(1) = 1
      start_k(1) = 1
      DO i = 1, nblks(1)
         start_i(i + 1) = start_i(i) + sizes_i(i)
      END DO
      DO j = 1, nblks(2)
         start_j(j + 1) = start_j(j) + sizes_j(j)
      END DO
      DO k = 1, nblks(3)
         start_k(k + 1) = start_k(k) + sizes_k(k)
      END DO
      CPASSERT(start_j(nblks(2) + 1) - 1 == SIZE(left_frame%coefficients, 2))
      CPASSERT(start_k(nblks(3) + 1) - 1 == SIZE(right_frame%coefficients, 2))
      auxiliary_start(1) = 1
      maxl = 0
      DO site = 1, SIZE(auxiliary_bases)
         CPASSERT(ASSOCIATED(auxiliary_bases(site)%gto_basis_set))
         ASSOCIATE (basis => auxiliary_bases(site)%gto_basis_set)
            auxiliary_start(site + 1) = auxiliary_start(site) + basis%nsgf
            maxl = MAX(maxl, MAXVAL(basis%lmax))
         END ASSOCIATE
      END DO
      CPASSERT(auxiliary_start(SIZE(auxiliary_bases) + 1) == start_i(nblks(1) + 1))
   END SUBROUTINE tensor_partition

! **************************************************************************************************
!> \brief Check identical native block shapes and ownership for paired real tensors.
!> \param reference ...
!> \param TARGET ...
! **************************************************************************************************
   SUBROUTINE matching_tensor(reference, TARGET)
      TYPE(dbt_type), INTENT(IN)                         :: reference, TARGET

      INTEGER                                            :: i, j, k, nr(3), nt(3), owner_r, owner_t
      INTEGER, ALLOCATABLE                               :: ri(:), rj(:), rk(:), ti(:), tj(:), tk(:)

      CALL dbt_get_info(reference, nblks_total=nr)
      CALL dbt_get_info(TARGET, nblks_total=nt)
      CPASSERT(ALL(nr == nt))
      ALLOCATE (ri(nr(1)), rj(nr(2)), rk(nr(3)), ti(nt(1)), tj(nt(2)), tk(nt(3)))
      CALL dbt_get_info(reference, blk_size_1=ri, blk_size_2=rj, blk_size_3=rk)
      CALL dbt_get_info(TARGET, blk_size_1=ti, blk_size_2=tj, blk_size_3=tk)
      CPASSERT(ALL(ri == ti) .AND. ALL(rj == tj) .AND. ALL(rk == tk))
      CPASSERT(reference%pgrid%mp_comm_2d%mepos == target%pgrid%mp_comm_2d%mepos)
      CPASSERT(reference%pgrid%mp_comm_2d%num_pe == target%pgrid%mp_comm_2d%num_pe)
      DO i = 1, nr(1)
      DO j = 1, nr(2)
      DO k = 1, nr(3)
         CALL dbt_get_stored_coordinates(reference, [i, j, k], owner_r)
         CALL dbt_get_stored_coordinates(TARGET, [i, j, k], owner_t)
         CPASSERT(owner_r == owner_t)
      END DO
      END DO
      END DO
   END SUBROUTINE matching_tensor

! **************************************************************************************************
!> \brief One original center/image contribution to an RI shell and selected native AO block.
!>        Sum over original centers/images and add the ordinary Gaussian three-center integrals.
!>        The soft-minus-Gaussian product is
!>        1/2 [d_left^* (s_right+chi_right) + (s_left+chi_left)^* d_right].
!>        The original Q contribution is added independently. Radial rows are bounded by the
!>        exact support of Delta, including its full nonzero tail. No numerical cutoff is applied.
!> \param left_frame ...
!> \param right_frame ...
!> \param center Common native center index in both frames.
!> \param left_range Contiguous native left AO indices, inclusive.
!> \param right_range Contiguous native right AO indices, inclusive.
!> \param fields Owned Coulomb potentials of the native RI shell.
!> \param auxiliary_center Cartesian position of the selected RI Gaussian image.
!> \param augmentation Original Q in the native magnetic-projector order.
!> \param r Complete original radial mesh for Delta.
!> \param radial_weights Fixed original radial measure, including r**2.
!> \param directions Unit angular directions, indexed (axis,point).
!> \param angular_weights Angular weights, summing to 4*pi.
!> \param integral New rectangular (left AO,right AO,RI function) correction.
!> \param response Optional new joint geometry/Bloch response, with all response inputs.
!> \param atom_response Cartesian motion of every native atom, including Gaussian centers.
!> \param cell_response Cartesian lattice motion.
!> \param left_k_response ...
!> \param right_k_response ...
!> \param left_coefficients_response dW from the original defining left frame.
!> \param right_coefficients_response dW from the original defining right frame.
!> \param auxiliary_response Cartesian motion of the RI Gaussian image.
!> \param image_shift Original center's integer image, zero by default.
!> \param auxiliary_images Optional fixed auxiliary Gaussian images, at fixed Gamma.
!> \param augmentation_integrals Optional workspace valid for this center/image, geometry and RI shell.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_hfx_center(left_frame, right_frame, center, left_range, right_range, fields, &
                                        auxiliary_center, augmentation, r, radial_weights, directions, angular_weights, &
                                        integral, response, atom_response, cell_response, left_k_response, right_k_response, &
                                        left_coefficients_response, right_coefficients_response, auxiliary_response, image_shift, &
                                        auxiliary_images, augmentation_integrals)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: left_frame, right_frame
      INTEGER, INTENT(IN)                                :: center, left_range(2), right_range(2)
      TYPE(gaussian_coulomb_fields_type), INTENT(IN)     :: fields
      REAL(dp), INTENT(IN)                               :: auxiliary_center(3)
      TYPE(uspp_augmentation_type), INTENT(IN)           :: augmentation
      REAL(dp), INTENT(IN)                               :: r(:), radial_weights(:), &
                                                            directions(:, :), angular_weights(:)
      COMPLEX(dp), INTENT(OUT)                           :: integral(:, :, :)
      COMPLEX(dp), INTENT(OUT), OPTIONAL                 :: response(:, :, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: atom_response(:, :), &
                                                            cell_response(3, 3), &
                                                            left_k_response(3), right_k_response(3)
      COMPLEX(dp), INTENT(IN), OPTIONAL :: left_coefficients_response(:, :), &
         right_coefficients_response(:, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: auxiliary_response(3)
      INTEGER, INTENT(IN), OPTIONAL                      :: image_shift(3), auxiliary_images(:, :)
      TYPE(transition_augmentation_type), &
         INTENT(INOUT), OPTIONAL, TARGET                 :: augmentation_integrals

      COMPLEX(dp)                                        :: dphase, phase
      COMPLEX(dp), ALLOCATABLE :: block(:, :), cl(:, :), cr(:, :), dl(:, :), dq(:, :), dr(:, :), &
         fcl(:, :), fcr(:, :), fdl(:, :), fdr(:, :), q(:, :), work(:, :)
      INTEGER                                            :: atom, aux, axis, ia, ir, j, kind, nang, &
                                                            naux, nl, nr, p, projectors, shift(3), &
                                                            t
      LOGICAL                                            :: derivatives
      REAL(dp)                                           :: dposition(3), position(3), radius, scale
      REAL(dp), ALLOCATABLE :: dpoints(:, :), points(:, :), potential_cell_gradients(:, :, :, :), &
         potential_gradients(:, :, :), potential_responses(:, :), potentials(:, :), &
         weighted_potential(:), weights(:)
      TYPE(inverse_field_workspace_type)                 :: left_fields, right_fields
      TYPE(transition_augmentation_type), POINTER        :: qi
      TYPE(transition_augmentation_type), TARGET         :: temporary_augmentation

      IF (PRESENT(auxiliary_images)) THEN
         CPASSERT(SIZE(auxiliary_images, 1) == 3)
         CPASSERT(ALL(left_frame%reduced_k == 0.0_dp))
         CPASSERT(ALL(right_frame%reduced_k == 0.0_dp))
      END IF
      CPASSERT(center >= 1 .AND. center <= SIZE(left_frame%center_atoms))
      CPASSERT(center <= SIZE(right_frame%center_atoms))
      atom = left_frame%center_atoms(center)
      kind = left_frame%center_kinds(center)
      IF (.NOT. (right_frame%center_atoms(center) == atom .AND. right_frame%center_kinds(center) == kind)) THEN
         CPABORT("Inconsistent transition center")
      END IF
      CPASSERT(ALL(SHAPE(left_frame%positions) == SHAPE(right_frame%positions)))
      CPASSERT(ALL(left_frame%positions == right_frame%positions))
      CPASSERT(ALL(left_frame%cell == right_frame%cell))
      CPASSERT(ALL(left_frame%periodicity == right_frame%periodicity))
      CPASSERT(left_range(1) >= 1)
      CPASSERT(left_range(2) <= SIZE(left_frame%coefficients, 2))
      CPASSERT(right_range(1) >= 1)
      CPASSERT(right_range(2) <= SIZE(right_frame%coefficients, 2))
      nl = left_range(2) - left_range(1) + 1
      nr = right_range(2) - right_range(1) + 1
      nang = SIZE(angular_weights)
      naux = fields%nfunction
      CPASSERT(nl > 0 .AND. nr > 0 .AND. nang > 0 .AND. naux > 0)
      CPASSERT(ALL(SHAPE(integral) == [nl, nr, naux]))
      CPASSERT(SIZE(r) >= 3 .AND. SIZE(radial_weights) == SIZE(r))
      CPASSERT(ALL(r >= 0.0_dp) .AND. ALL(r(2:) > r(:SIZE(r) - 1)))
      CPASSERT(ALL(SHAPE(directions) == [3, nang]))
      CPASSERT(ABS(SUM(angular_weights) - fourpi) <= 1.0e-10_dp)
      CPASSERT(MAXVAL(ABS(SUM(directions**2, DIM=1) - 1.0_dp)) <= 1.0e-12_dp)
      projectors = left_frame%projector_last(center) - left_frame%projector_first(center) + 1
      IF (.NOT. (projectors == right_frame%projector_last(center) - right_frame%projector_first(center) + 1)) THEN
         CPABORT("Inconsistent transition projector count")
      END IF
      CPASSERT(augmentation%nproj == 0 .OR. augmentation%nproj == projectors)
      shift = 0
      IF (PRESENT(image_shift)) shift = image_shift
      CPASSERT(ALL(shift*left_frame%periodicity == shift))
      position = left_frame%positions(:, atom) + MATMUL(left_frame%cell, REAL(shift, dp))
      derivatives = PRESENT(response)
      IF (derivatives) THEN
         CPASSERT(ALL(SHAPE(response) == SHAPE(integral)))
         CPASSERT(PRESENT(atom_response) .AND. PRESENT(cell_response))
         CPASSERT(PRESENT(left_k_response) .AND. PRESENT(right_k_response))
         IF (PRESENT(auxiliary_images)) THEN
            CPASSERT(ALL(left_k_response == 0.0_dp))
            CPASSERT(ALL(right_k_response == 0.0_dp))
         END IF
         IF (.NOT. (PRESENT(left_coefficients_response) .AND. PRESENT(right_coefficients_response))) THEN
            CPABORT("Both coefficient responses are required")
         END IF
         CPASSERT(PRESENT(auxiliary_response))
         CPASSERT(ALL(SHAPE(atom_response) == SHAPE(left_frame%positions)))
         IF (.NOT. ALL(SHAPE(left_coefficients_response) == SHAPE(left_frame%coefficients))) THEN
            CPABORT("Inconsistent left coefficient response")
         END IF
         IF (.NOT. ALL(SHAPE(right_coefficients_response) == SHAPE(right_frame%coefficients))) THEN
            CPABORT("Inconsistent right coefficient response")
         END IF
         dposition = atom_response(:, atom) + MATMUL(cell_response, REAL(shift, dp))
         response = CMPLX(0.0_dp, 0.0_dp, dp)
      END IF
      integral = CMPLX(0.0_dp, 0.0_dp, dp)
      ALLOCATE (points(3, nang), potentials(nang, naux), weights(nang), weighted_potential(nang), work(nang, nr), block(nl, nr))
      IF (derivatives) THEN
         ALLOCATE (dpoints(3, nang), potential_gradients(3, naux, nang), potential_responses(nang, naux))
         IF (PRESENT(auxiliary_images)) ALLOCATE (potential_cell_gradients(3, 3, naux, nang))
         DO ia = 1, nang
            dpoints(:, ia) = dposition
         END DO
      END IF
      radius = MAX(left_frame%delta_cache(kind)%radius, right_frame%delta_cache(kind)%radius)
      CPASSERT(r(SIZE(r)) >= radius)
      DO ir = 1, SIZE(r)
         IF (r(ir) >= radius .OR. radial_weights(ir) == 0.0_dp) CYCLE
         DO ia = 1, nang
            points(:, ia) = position + r(ir)*directions(:, ia)
         END DO
         weights(:) = 0.5_dp*radial_weights(ir)*angular_weights
         CALL residual_fields(left_frame, center, left_range, points, cl, dl, fcl, fdl, &
                              dpoints, atom_response, cell_response, left_k_response, left_coefficients_response, &
                              derivatives, shift, left_fields)
         CALL residual_fields(right_frame, center, right_range, points, cr, dr, fcr, fdr, &
                              dpoints, atom_response, cell_response, right_k_response, right_coefficients_response, &
                              derivatives, shift, right_fields)
         DO ia = 1, nang
            IF (derivatives) THEN
               IF (PRESENT(auxiliary_images)) THEN
                  CALL gaussian_coulomb_fields_periodic(fields, points(:, ia) - auxiliary_center, left_frame%cell, &
                                                        auxiliary_images, potentials(ia, :), potential_gradients(:, :, ia), &
                                                        potential_cell_gradients(:, :, :, ia))
               ELSE
                  CALL gaussian_coulomb_fields_evaluate(fields, points(:, ia) - auxiliary_center, &
                                                        potentials(ia, :), potential_gradients(:, :, ia))
               END IF
               potential_responses(ia, :) = MATMUL(dposition - auxiliary_response, potential_gradients(:, :, ia))
               IF (PRESENT(auxiliary_images)) THEN
                  DO j = 1, 3
                     DO axis = 1, 3
                        potential_responses(ia, :) = potential_responses(ia, :) + &
                                                     cell_response(axis, j)*potential_cell_gradients(axis, j, :, ia)
                     END DO
                  END DO
               END IF
            ELSE
               IF (PRESENT(auxiliary_images)) THEN
                  CALL gaussian_coulomb_fields_periodic(fields, points(:, ia) - auxiliary_center, left_frame%cell, &
                                                        auxiliary_images, potentials(ia, :))
               ELSE
                  CALL gaussian_coulomb_fields_evaluate(fields, points(:, ia) - auxiliary_center, potentials(ia, :))
               END IF
            END IF
         END DO
         DO aux = 1, naux
            weighted_potential(:) = weights*potentials(:, aux)
            DO j = 1, nr
               work(:, j) = cr(:, j)*weighted_potential
            END DO
            block(:, :) = MATMUL(CONJG(TRANSPOSE(dl)), work)
            DO j = 1, nr
               work(:, j) = dr(:, j)*weighted_potential
            END DO
            integral(:, :, aux) = integral(:, :, aux) + block + MATMUL(CONJG(TRANSPOSE(cl)), work)
            IF (.NOT. derivatives) CYCLE
            block(:, :) = MATMUL(CONJG(TRANSPOSE(fcl)), work)
            DO j = 1, nr
               work(:, j) = fdr(:, j)*weighted_potential
            END DO
            block(:, :) = block + MATMUL(CONJG(TRANSPOSE(cl)), work)
            DO j = 1, nr
               work(:, j) = cr(:, j)*weighted_potential
            END DO
            block(:, :) = block + MATMUL(CONJG(TRANSPOSE(fdl)), work)
            DO j = 1, nr
               work(:, j) = fcr(:, j)*weighted_potential
            END DO
            block(:, :) = block + MATMUL(CONJG(TRANSPOSE(dl)), work)
            DO ia = 1, nang
               scale = weights(ia)*potential_responses(ia, aux)
               work(ia, :) = scale*cr(ia, :)
            END DO
            block(:, :) = block + MATMUL(CONJG(TRANSPOSE(dl)), work)
            DO ia = 1, nang
               scale = weights(ia)*potential_responses(ia, aux)
               work(ia, :) = scale*dr(ia, :)
            END DO
            response(:, :, aux) = response(:, :, aux) + block + MATMUL(CONJG(TRANSPOSE(cl)), work)
         END DO
      END DO
      IF (augmentation%nproj == 0) RETURN
      ALLOCATE (q(projectors, projectors), dq(projectors, projectors))
      qi => temporary_augmentation
      IF (PRESENT(augmentation_integrals)) qi => augmentation_integrals
      CALL prepare_transition_augmentation(fields, augmentation, position - auxiliary_center, directions, angular_weights, &
                                           qi, derivatives, left_frame%cell, auxiliary_images)
      phase = EXP(CMPLX(0.0_dp, twopi*DOT_PRODUCT(right_frame%reduced_k - left_frame%reduced_k, REAL(shift, dp)), dp))
      IF (derivatives) dphase = phase*CMPLX(0.0_dp, &
                                            twopi*DOT_PRODUCT(right_k_response - left_k_response, REAL(shift, dp)), dp)
      p = left_frame%projector_first(center)
      t = right_frame%projector_first(center)
      ASSOCIATE (wl => left_frame%coefficients(p:p + projectors - 1, left_range(1):left_range(2)), &
                 wr => right_frame%coefficients(t:t + projectors - 1, right_range(1):right_range(2)))
         DO aux = 1, naux
            q(:, :) = phase*qi%values(:, :, aux)
            integral(:, :, aux) = integral(:, :, aux) + MATMUL(CONJG(TRANSPOSE(wl)), MATMUL(q, wr))
            IF (.NOT. derivatives) CYCLE
            dq(:, :) = dphase*qi%values(:, :, aux)
            DO ia = 1, 3
               dq(:, :) = dq + phase*qi%gradient(:, :, aux, ia)*(dposition(ia) - auxiliary_response(ia))
            END DO
            IF (PRESENT(auxiliary_images)) THEN
               DO j = 1, 3
                  DO axis = 1, 3
                     dq(:, :) = dq + phase*qi%cell_gradient(:, :, aux, axis, j)*cell_response(axis, j)
                  END DO
               END DO
            END IF
            response(:, :, aux) = response(:, :, aux) + MATMUL(CONJG(TRANSPOSE(wl)), MATMUL(dq, wr)) + &
                           MATMUL(CONJG(TRANSPOSE(left_coefficients_response(p:p + projectors - 1, left_range(1):left_range(2)))), &
                                         MATMUL(q, wr)) + &
                                  MATMUL(CONJG(TRANSPOSE(wl)), &
                                        MATMUL(q, right_coefficients_response(t:t + projectors - 1, right_range(1):right_range(2))))
         END DO
      END ASSOCIATE
   END SUBROUTINE qs_upf_inverse_hfx_center

! **************************************************************************************************
!> \brief One unphased transition correction for Gaussian AO images 0 and R and center image T.
!>        DIRECT retains the original Q. INVERSE_1C also retains both Gaussian/Delta products
!>        and the Delta/Delta product. Its disjoint reconstruction domains are checked when
!>        preparing the frame. The general inverse requires its complete lattice transformation.
!>        The native RI owner must sum centers/images before fitting or filtering this correction.
!> \param frame Prepared original DIRECT or INVERSE_1C data and unshifted Gaussian centers.
!> \param center Original center index.
!> \param left_range Selected Gaussian AO columns in image 0.
!> \param right_range Selected Gaussian AO columns in image R.
!> \param ao_image Fixed integer right Gaussian image R.
!> \param center_image Fixed integer reconstruction image T.
!> \param fields Coulomb potentials of the native RI shell.
!> \param auxiliary_center Actual Cartesian RI image position.
!> \param augmentation Original transition augmentation in magnetic-projector order.
!> \param r Complete original radial mesh.
!> \param radial_weights Original radial measure including r**2.
!> \param directions Unit angular directions.
!> \param angular_weights Angular weights summing to 4*pi.
!> \param integral Rectangular AO/RI correction without Bloch phases or occupations.
!> \param response Optional complete atomic/cell/RI-image response at fixed image indices.
!> \param atom_response Motion of unshifted atom centers.
!> \param cell_response Motion of the lattice vectors.
!> \param auxiliary_response Motion of the actual RI image position.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_hfx_lattice_center(frame, center, left_range, right_range, ao_image, center_image, &
                                                fields, auxiliary_center, augmentation, r, radial_weights, &
                                                directions, angular_weights, integral, response, atom_response, &
                                                cell_response, auxiliary_response)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      INTEGER, INTENT(IN)                                :: center, left_range(2), right_range(2), &
                                                            ao_image(3), center_image(3)
      TYPE(gaussian_coulomb_fields_type), INTENT(IN)     :: fields
      REAL(dp), INTENT(IN)                               :: auxiliary_center(3)
      TYPE(uspp_augmentation_type), INTENT(IN)           :: augmentation
      REAL(dp), INTENT(IN)                               :: r(:), radial_weights(:), &
                                                            directions(:, :), angular_weights(:)
      COMPLEX(dp), INTENT(OUT)                           :: integral(:, :, :)
      COMPLEX(dp), INTENT(OUT), OPTIONAL                 :: response(:, :, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: atom_response(:, :), &
                                                            cell_response(3, 3), &
                                                            auxiliary_response(3)

      TYPE(lattice_overlap_type)                         :: overlaps

      CALL lattice_center_with_overlaps(frame, center, left_range, right_range, ao_image, center_image, &
                                        fields, auxiliary_center, augmentation, r, radial_weights, directions, &
                                        angular_weights, integral, response, atom_response, cell_response, &
                                        auxiliary_response, overlaps)
   END SUBROUTINE qs_upf_inverse_hfx_lattice_center

! **************************************************************************************************
!> \brief Evaluate one lattice center with reusable relative-image projector overlaps.
!>        The private workspace belongs to one AO block pair, physical image and geometry call.
!> \param overlaps ...
! **************************************************************************************************
   SUBROUTINE lattice_center_with_overlaps(frame, center, left_range, right_range, ao_image, center_image, &
                                           fields, auxiliary_center, augmentation, r, radial_weights, &
                                           directions, angular_weights, integral, response, atom_response, &
                                           cell_response, auxiliary_response, overlaps)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      INTEGER, INTENT(IN)                                :: center, left_range(2), right_range(2), &
                                                            ao_image(3), center_image(3)
      TYPE(gaussian_coulomb_fields_type), INTENT(IN)     :: fields
      REAL(dp), INTENT(IN)                               :: auxiliary_center(3)
      TYPE(uspp_augmentation_type), INTENT(IN)           :: augmentation
      REAL(dp), INTENT(IN)                               :: r(:), radial_weights(:), &
                                                            directions(:, :), angular_weights(:)
      COMPLEX(dp), INTENT(OUT)                           :: integral(:, :, :)
      COMPLEX(dp), INTENT(OUT), OPTIONAL                 :: response(:, :, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: atom_response(:, :), &
                                                            cell_response(3, 3), &
                                                            auxiliary_response(3)
      TYPE(lattice_overlap_type), INTENT(INOUT)          :: overlaps

      COMPLEX(dp), ALLOCATABLE :: block(:, :), cl(:, :), cr(:, :), delta(:, :), dl(:, :), &
         dq(:, :), dr(:, :), dwl(:, :), dwr(:, :), fcl(:, :), fcr(:, :), fdl(:, :), fdr(:, :), &
         q(:, :), wl(:, :), work(:, :), wr(:, :)
      INTEGER                                            :: atom, aux, ia, ir, j, kind, nang, nl, &
                                                            np, nr, zero(3, 1)
      LOGICAL                                            :: derivatives
      REAL(dp)                                           :: dposition(3), position(3), radius, scale
      REAL(dp), ALLOCATABLE :: dpoints(:, :), dqi(:, :, :, :), local_points(:, :), points(:, :), &
         potential_gradient(:, :, :), potentials(:, :), qi(:, :, :), weighted_potential(:), &
         weights(:)
      TYPE(inverse_field_workspace_type)                 :: left_fields, right_fields

      CPASSERT(frame%atom_local)
      CPASSERT(center >= 1 .AND. center <= SIZE(frame%center_atoms))
      CPASSERT(ALL(ao_image*frame%periodicity == ao_image))
      CPASSERT(ALL(center_image*frame%periodicity == center_image))
      atom = frame%center_atoms(center)
      kind = frame%center_kinds(center)
      np = frame%projector_last(center) - frame%projector_first(center) + 1
      nl = left_range(2) - left_range(1) + 1
      nr = right_range(2) - right_range(1) + 1
      nang = SIZE(angular_weights)
      CPASSERT(nl > 0 .AND. nr > 0 .AND. nang > 0 .AND. fields%nfunction > 0)
      CPASSERT(left_range(1) >= 1 .AND. left_range(2) <= SIZE(frame%coefficients, 2))
      CPASSERT(right_range(1) >= 1 .AND. right_range(2) <= SIZE(frame%coefficients, 2))
      CPASSERT(ALL(SHAPE(integral) == [nl, nr, fields%nfunction]))
      CPASSERT(augmentation%nproj == 0 .OR. augmentation%nproj == np)
      CPASSERT(SIZE(r) >= 3 .AND. SIZE(radial_weights) == SIZE(r))
      CPASSERT(ALL(r >= 0.0_dp) .AND. ALL(r(2:) > r(:SIZE(r) - 1)))
      CPASSERT(ALL(SHAPE(directions) == [3, nang]))
      CPASSERT(ABS(SUM(angular_weights) - fourpi) <= 1.0e-10_dp)
      CPASSERT(MAXVAL(ABS(SUM(directions**2, DIM=1) - 1.0_dp)) <= 1.0e-12_dp)
      derivatives = PRESENT(response)
      dposition = 0.0_dp
      position = frame%positions(:, atom) + MATMUL(frame%cell, REAL(center_image, dp))
      IF (derivatives) THEN
         CPASSERT(ALL(SHAPE(response) == SHAPE(integral)))
         IF (.NOT. (PRESENT(atom_response) .AND. PRESENT(cell_response) .AND. PRESENT(auxiliary_response))) THEN
            CPABORT("All lattice geometry responses are required")
         END IF
         CPASSERT(ALL(SHAPE(atom_response) == SHAPE(frame%positions)))
         dposition = atom_response(:, atom) + MATMUL(cell_response, REAL(center_image, dp))
         response = CMPLX(0.0_dp, 0.0_dp, dp)
      END IF
      IF (.NOT. ALLOCATED(overlaps%left)) THEN
         ALLOCATE (overlaps%left(np, nl), overlaps%right(np, nr))
         IF (derivatives) THEN
            ALLOCATE (overlaps%dleft(np, nl), overlaps%dright(np, nr))
            CALL lattice_coefficients(frame, center, left_range, -center_image, &
                                      overlaps%left, overlaps%dleft, atom_response, cell_response)
            CALL lattice_coefficients(frame, center, right_range, ao_image - center_image, &
                                      overlaps%right, overlaps%dright, atom_response, cell_response)
         ELSE
            CALL lattice_coefficients(frame, center, left_range, -center_image, overlaps%left)
            CALL lattice_coefficients(frame, center, right_range, ao_image - center_image, overlaps%right)
         END IF
      END IF
      CPASSERT(ALL(SHAPE(overlaps%left) == [np, nl]))
      CPASSERT(ALL(SHAPE(overlaps%right) == [np, nr]))
      ALLOCATE (wl(np, nl), wr(np, nr))
      wl(:, :) = overlaps%left
      wr(:, :) = overlaps%right
      IF (derivatives) THEN
         CPASSERT(ALLOCATED(overlaps%dleft) .AND. ALLOCATED(overlaps%dright))
         ALLOCATE (dwl(np, nl), dwr(np, nr))
         dwl(:, :) = overlaps%dleft
         dwr(:, :) = overlaps%dright
      END IF
      integral = CMPLX(0.0_dp, 0.0_dp, dp)
      radius = frame%delta_cache(kind)%radius
      CPASSERT(r(SIZE(r)) >= radius)
      IF (.NOT. frame%direct .AND. radius > 0.0_dp) THEN
         ALLOCATE (points(3, nang), local_points(3, nang), potentials(nang, fields%nfunction), &
                   weights(nang), weighted_potential(nang), delta(nang, np), work(nang, nr), &
                   block(nl, nr), cl(nang, nl), dl(nang, nl), cr(nang, nr), dr(nang, nr))
         IF (derivatives) THEN
            ALLOCATE (potential_gradient(3, fields%nfunction, nang), dpoints(3, nang), &
                      fcl(nang, nl), fdl(nang, nl), fcr(nang, nr), fdr(nang, nr))
            DO ia = 1, nang
               dpoints(:, ia) = dposition
            END DO
         END IF
         zero = 0
         DO ir = 1, SIZE(r)
            IF (r(ir) >= radius .OR. radial_weights(ir) == 0.0_dp) CYCLE
            local_points(:, :) = r(ir)*directions
            DO ia = 1, nang
               points(:, ia) = position + local_points(:, ia)
            END DO
            ! The integration coordinates move with the reconstruction center. Delta and
            ! its radial/angular measure therefore have no geometry response on this sphere.
            CALL upf_radial_fields_periodic(frame%delta_cache(kind), local_points, &
                                            [0.0_dp, 0.0_dp, 0.0_dp], frame%cell, zero, &
                                            [0.0_dp, 0.0_dp, 0.0_dp], delta)
            CALL lattice_residual_fields(frame, left_range, [0, 0, 0], points, delta, wl, cl, dl, &
                                         left_fields, fcl, fdl, dwl, atom_response, cell_response, dpoints, derivatives)
            CALL lattice_residual_fields(frame, right_range, ao_image, points, delta, wr, cr, dr, &
                                         right_fields, fcr, fdr, dwr, atom_response, cell_response, dpoints, derivatives)
            weights(:) = 0.5_dp*radial_weights(ir)*angular_weights
            DO ia = 1, nang
               IF (derivatives) THEN
                  CALL gaussian_coulomb_fields_evaluate(fields, points(:, ia) - auxiliary_center, &
                                                        potentials(ia, :), potential_gradient(:, :, ia))
               ELSE
                  CALL gaussian_coulomb_fields_evaluate(fields, points(:, ia) - auxiliary_center, potentials(ia, :))
               END IF
            END DO
            DO aux = 1, fields%nfunction
               weighted_potential(:) = weights*potentials(:, aux)
               DO j = 1, nr
                  work(:, j) = cr(:, j)*weighted_potential
               END DO
               block(:, :) = MATMUL(CONJG(TRANSPOSE(dl)), work)
               DO j = 1, nr
                  work(:, j) = dr(:, j)*weighted_potential
               END DO
               integral(:, :, aux) = integral(:, :, aux) + block + MATMUL(CONJG(TRANSPOSE(cl)), work)
               IF (.NOT. derivatives) CYCLE
               block(:, :) = MATMUL(CONJG(TRANSPOSE(fcl)), work)
               DO j = 1, nr
                  work(:, j) = fdr(:, j)*weighted_potential
               END DO
               block(:, :) = block + MATMUL(CONJG(TRANSPOSE(cl)), work)
               DO j = 1, nr
                  work(:, j) = cr(:, j)*weighted_potential
               END DO
               block(:, :) = block + MATMUL(CONJG(TRANSPOSE(fdl)), work)
               DO j = 1, nr
                  work(:, j) = fcr(:, j)*weighted_potential
               END DO
               block(:, :) = block + MATMUL(CONJG(TRANSPOSE(dl)), work)
               DO ia = 1, nang
                  scale = weights(ia)*DOT_PRODUCT(potential_gradient(:, aux, ia), dposition - auxiliary_response)
                  work(ia, :) = scale*cr(ia, :)
               END DO
               block(:, :) = block + MATMUL(CONJG(TRANSPOSE(dl)), work)
               DO ia = 1, nang
                  scale = weights(ia)*DOT_PRODUCT(potential_gradient(:, aux, ia), dposition - auxiliary_response)
                  work(ia, :) = scale*dr(ia, :)
               END DO
               response(:, :, aux) = response(:, :, aux) + block + MATMUL(CONJG(TRANSPOSE(cl)), work)
            END DO
         END DO
      END IF
      IF (augmentation%nproj == 0) RETURN
      ALLOCATE (qi(np, np, fields%nfunction), q(np, np))
      IF (derivatives) THEN
         ALLOCATE (dqi(np, np, fields%nfunction, 3), dq(np, np))
         CALL upf_inverse_exchange_augmentation(fields, augmentation, position - auxiliary_center, &
                                                directions, angular_weights, qi, augmentation_derivative=dqi)
      ELSE
         CALL upf_inverse_exchange_augmentation(fields, augmentation, position - auxiliary_center, &
                                                directions, angular_weights, qi)
      END IF
      DO aux = 1, fields%nfunction
         q(:, :) = CMPLX(qi(:, :, aux), 0.0_dp, dp)
         integral(:, :, aux) = integral(:, :, aux) + MATMUL(CONJG(TRANSPOSE(wl)), MATMUL(q, wr))
         IF (.NOT. derivatives) CYCLE
         dq = CMPLX(0.0_dp, 0.0_dp, dp)
         DO ia = 1, 3
            dq(:, :) = dq + dqi(:, :, aux, ia)*(dposition(ia) - auxiliary_response(ia))
         END DO
         response(:, :, aux) = response(:, :, aux) + MATMUL(CONJG(TRANSPOSE(wl)), MATMUL(dq, wr)) + &
                               MATMUL(CONJG(TRANSPOSE(dwl)), MATMUL(q, wr)) + &
                               MATMUL(CONJG(TRANSPOSE(wl)), MATMUL(q, dwr))
      END DO
   END SUBROUTINE lattice_center_with_overlaps

! **************************************************************************************************
!> \brief Hard-projector overlaps for one relative Gaussian image, without a reciprocal phase.
!> \param frame ...
!> \param center ...
!> \param ao_range ...
!> \param image Gaussian image relative to the projector's cell.
!> \param coefficients ...
!> \param response ...
!> \param atom_response ...
!> \param cell_response ...
! **************************************************************************************************
   SUBROUTINE lattice_coefficients(frame, center, ao_range, image, coefficients, response, atom_response, cell_response)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      INTEGER, INTENT(IN)                                :: center, ao_range(2), image(3)
      COMPLEX(dp), INTENT(OUT)                           :: coefficients(:, :)
      COMPLEX(dp), INTENT(OUT), OPTIONAL                 :: response(:, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: atom_response(:, :), cell_response(3, 3)

      COMPLEX(dp), ALLOCATABLE                           :: dfull(:, :), full(:, :)
      INTEGER                                            :: atom, b, first, hi, images(3, 1), kind, &
                                                            last, lo, nproj

      atom = frame%center_atoms(center)
      kind = frame%center_kinds(center)
      nproj = frame%atomic_cache(kind)%hard%nproj
      images(:, 1) = image
      coefficients = CMPLX(0.0_dp, 0.0_dp, dp)
      IF (PRESENT(response)) response = CMPLX(0.0_dp, 0.0_dp, dp)
      DO b = 1, SIZE(frame%bases)
         first = frame%ao_first(b)
         last = frame%ao_last(b)
         lo = MAX(first, ao_range(1))
         hi = MIN(last, ao_range(2))
         IF (hi < lo) CYCLE
         ALLOCATE (full(nproj, last - first + 1))
         IF (PRESENT(response)) THEN
            ALLOCATE (dfull(nproj, last - first + 1))
            CALL qs_upf_inverse_atomic_overlap(frame%bases(b)%gto_basis_set, frame%atomic_cache(kind)%hard, &
                                               frame%positions(:, atom), frame%positions(:, b), frame%cell, images, &
                                               [0.0_dp, 0.0_dp, 0.0_dp], full, dfull, &
                                               atom_response(:, atom), atom_response(:, b), cell_response)
            response(:, lo - ao_range(1) + 1:hi - ao_range(1) + 1) = dfull(:, lo - first + 1:hi - first + 1)
            DEALLOCATE (dfull)
         ELSE
            CALL qs_upf_inverse_atomic_overlap(frame%bases(b)%gto_basis_set, frame%atomic_cache(kind)%hard, &
                                               frame%positions(:, atom), frame%positions(:, b), frame%cell, images, &
                                               [0.0_dp, 0.0_dp, 0.0_dp], full)
         END IF
         coefficients(:, lo - ao_range(1) + 1:hi - ao_range(1) + 1) = full(:, lo - first + 1:hi - first + 1)
         DEALLOCATE (full)
      END DO
   END SUBROUTINE lattice_coefficients

! **************************************************************************************************
!> \brief Local INVERSE_1C soft-plus-Gaussian and difference fields for one Gaussian image.
!> \param frame ...
!> \param ao_range ...
!> \param image ...
!> \param points ...
!> \param delta Original difference field on the moving integration sphere.
!> \param coefficients Hard projector overlaps for this AO image.
!> \param complete ...
!> \param difference ...
!> \param workspace ...
!> \param complete_response ...
!> \param difference_response ...
!> \param coefficients_response ...
!> \param atom_response ...
!> \param cell_response ...
!> \param point_response ...
!> \param derivatives ...
! **************************************************************************************************
   SUBROUTINE lattice_residual_fields(frame, ao_range, image, points, delta, coefficients, complete, difference, &
                                      workspace, complete_response, difference_response, coefficients_response, &
                                      atom_response, cell_response, point_response, derivatives)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      INTEGER, INTENT(IN)                                :: ao_range(2), image(3)
      REAL(dp), INTENT(IN)                               :: points(:, :)
      COMPLEX(dp), INTENT(IN)                            :: delta(:, :), coefficients(:, :)
      COMPLEX(dp), INTENT(OUT)                           :: complete(:, :), difference(:, :)
      TYPE(inverse_field_workspace_type), INTENT(INOUT)  :: workspace
      COMPLEX(dp), ALLOCATABLE, INTENT(INOUT)            :: complete_response(:, :), &
                                                            difference_response(:, :), &
                                                            coefficients_response(:, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: atom_response(:, :), cell_response(3, 3)
      REAL(dp), ALLOCATABLE, INTENT(IN)                  :: point_response(:, :)
      LOGICAL, INTENT(IN)                                :: derivatives

      INTEGER                                            :: atom, first, hi, images(3, 1), last, lo, &
                                                            max_basis, np

      np = SIZE(points, 2)
      images(:, 1) = image
      IF (.NOT. ALLOCATED(workspace%phi)) THEN
         max_basis = MAXVAL(frame%ao_last - frame%ao_first + 1)
         ALLOCATE (workspace%phi(np, max_basis))
         IF (derivatives) ALLOCATE (workspace%dphi(np, max_basis))
      END IF
      complete = CMPLX(0.0_dp, 0.0_dp, dp)
      IF (derivatives) THEN
         complete_response = CMPLX(0.0_dp, 0.0_dp, dp)
      END IF
      DO atom = 1, SIZE(frame%bases)
         first = frame%ao_first(atom)
         last = frame%ao_last(atom)
         lo = MAX(first, ao_range(1))
         hi = MIN(last, ao_range(2))
         IF (hi < lo) CYCLE
         IF (derivatives) THEN
            CALL gaussian_basis_fields_periodic(frame%bases(atom)%gto_basis_set, points, frame%positions(:, atom), &
                                                frame%cell, images, [0.0_dp, 0.0_dp, 0.0_dp], &
                                                workspace%phi(:, :last - first + 1), &
                                                response=workspace%dphi(:, :last - first + 1), &
                                                point_response=point_response, center_response=atom_response(:, atom), &
                                                cell_response=cell_response)
            complete_response(:, lo - ao_range(1) + 1:hi - ao_range(1) + 1) = &
               2.0_dp*workspace%dphi(:, lo - first + 1:hi - first + 1)
         ELSE
            CALL gaussian_basis_fields_periodic(frame%bases(atom)%gto_basis_set, points, frame%positions(:, atom), &
                                                frame%cell, images, [0.0_dp, 0.0_dp, 0.0_dp], &
                                                workspace%phi(:, :last - first + 1))
         END IF
         complete(:, lo - ao_range(1) + 1:hi - ao_range(1) + 1) = 2.0_dp*workspace%phi(:, lo - first + 1:hi - first + 1)
      END DO
      difference(:, :) = -MATMUL(delta, coefficients)
      complete(:, :) = complete + difference
      IF (derivatives) THEN
         difference_response(:, :) = -MATMUL(delta, coefficients_response)
         complete_response(:, :) = complete_response + difference_response
      END IF
   END SUBROUTINE lattice_residual_fields

! **************************************************************************************************
!> \brief Return the coefficient covectors of one original center/image RI correction.
!>        The real energy variation is Re sum(conjg(G)*dC). Both coefficient factors and
!>        the complete intercenter soft fields are retained, together with the independent Q.
!>        Geometry, radial data and the auxiliary shell are fixed in this partial adjoint.
!>        Its coefficient covectors must subsequently pass through their defining inverse solve.
!>        This owner-local routine performs no MPI reduction and adds no occupation factors.
!> \param left_frame ...
!> \param right_frame ...
!> \param center ...
!> \param left_range ...
!> \param right_range ...
!> \param fields ...
!> \param auxiliary_center ...
!> \param augmentation ...
!> \param r ...
!> \param radial_weights ...
!> \param directions ...
!> \param angular_weights ...
!> \param covector Owner's complex derivative, ordered (left AO,right AO,RI function).
!> \param left_coefficients_covector All left projector rows and selected left AO columns.
!> \param right_coefficients_covector All right projector rows and selected right AO columns.
!> \param image_shift ...
!> \param auxiliary_images Optional fixed auxiliary Gaussian images, at fixed Gamma.
!> \param augmentation_integrals Optional workspace valid for this center/image, geometry and RI shell.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_hfx_center_adjoint(left_frame, right_frame, center, left_range, right_range, &
                                                fields, auxiliary_center, augmentation, r, radial_weights, &
                                                directions, angular_weights, covector, left_coefficients_covector, &
                                                right_coefficients_covector, image_shift, auxiliary_images, augmentation_integrals)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: left_frame, right_frame
      INTEGER, INTENT(IN)                                :: center, left_range(2), right_range(2)
      TYPE(gaussian_coulomb_fields_type), INTENT(IN)     :: fields
      REAL(dp), INTENT(IN)                               :: auxiliary_center(3)
      TYPE(uspp_augmentation_type), INTENT(IN)           :: augmentation
      REAL(dp), INTENT(IN)                               :: r(:), radial_weights(:), &
                                                            directions(:, :), angular_weights(:)
      COMPLEX(dp), INTENT(IN)                            :: covector(:, :, :)
      COMPLEX(dp), INTENT(OUT) :: left_coefficients_covector(:, :), &
         right_coefficients_covector(:, :)
      INTEGER, INTENT(IN), OPTIONAL                      :: image_shift(3), auxiliary_images(:, :)
      TYPE(transition_augmentation_type), &
         INTENT(INOUT), OPTIONAL, TARGET                 :: augmentation_integrals

      COMPLEX(dp)                                        :: phase
      COMPLEX(dp), ALLOCATABLE :: cl(:, :), cr(:, :), dl(:, :), dr(:, :), gcl(:, :), gcr(:, :), &
         gdl(:, :), gdr(:, :), q(:, :), unused_c(:, :), unused_d(:, :), work(:, :)
      INTEGER                                            :: atom, aux, ia, ir, j, kind, nang, nl, &
                                                            nr, p, projectors, shift(3), t
      REAL(dp)                                           :: position(3), radius
      REAL(dp), ALLOCATABLE                              :: points(:, :), potentials(:, :), &
                                                            unused_points(:, :), &
                                                            weighted_potential(:), weights(:)
      TYPE(inverse_field_workspace_type)                 :: left_fields, right_fields
      TYPE(transition_augmentation_type), POINTER        :: qi
      TYPE(transition_augmentation_type), TARGET         :: temporary_augmentation

      IF (PRESENT(auxiliary_images)) THEN
         CPASSERT(SIZE(auxiliary_images, 1) == 3)
         CPASSERT(ALL(left_frame%reduced_k == 0.0_dp))
         CPASSERT(ALL(right_frame%reduced_k == 0.0_dp))
      END IF
      CPASSERT(center >= 1 .AND. center <= SIZE(left_frame%center_atoms))
      CPASSERT(center <= SIZE(right_frame%center_atoms))
      atom = left_frame%center_atoms(center)
      kind = left_frame%center_kinds(center)
      IF (.NOT. (right_frame%center_atoms(center) == atom .AND. right_frame%center_kinds(center) == kind)) THEN
         CPABORT("Inconsistent transition center")
      END IF
      CPASSERT(ALL(SHAPE(left_frame%positions) == SHAPE(right_frame%positions)))
      CPASSERT(ALL(left_frame%positions == right_frame%positions))
      CPASSERT(ALL(left_frame%cell == right_frame%cell))
      CPASSERT(ALL(left_frame%periodicity == right_frame%periodicity))
      CPASSERT(left_range(1) >= 1 .AND. left_range(2) <= SIZE(left_frame%coefficients, 2))
      IF (.NOT. (right_range(1) >= 1 .AND. right_range(2) <= SIZE(right_frame%coefficients, 2))) THEN
         CPABORT("Invalid transition coefficient range")
      END IF
      nl = left_range(2) - left_range(1) + 1
      nr = right_range(2) - right_range(1) + 1
      nang = SIZE(angular_weights)
      CPASSERT(nl > 0 .AND. nr > 0 .AND. nang > 0 .AND. fields%nfunction > 0)
      CPASSERT(ALL(SHAPE(covector) == [nl, nr, fields%nfunction]))
      IF (.NOT. (ALL(SHAPE(left_coefficients_covector) == [SIZE(left_frame%coefficients, 1), nl]))) THEN
         CPABORT("Inconsistent left-coefficient covector shape")
      END IF
      IF (.NOT. (ALL(SHAPE(right_coefficients_covector) == [SIZE(right_frame%coefficients, 1), nr]))) THEN
         CPABORT("Inconsistent right-coefficient covector shape")
      END IF
      CPASSERT(SIZE(r) >= 3 .AND. SIZE(radial_weights) == SIZE(r))
      CPASSERT(ALL(r >= 0.0_dp) .AND. ALL(r(2:) > r(:SIZE(r) - 1)))
      CPASSERT(ALL(SHAPE(directions) == [3, nang]))
      CPASSERT(ABS(SUM(angular_weights) - fourpi) <= 1.0e-10_dp)
      CPASSERT(MAXVAL(ABS(SUM(directions**2, DIM=1) - 1.0_dp)) <= 1.0e-12_dp)
      projectors = left_frame%projector_last(center) - left_frame%projector_first(center) + 1
      IF (.NOT. (projectors == right_frame%projector_last(center) - right_frame%projector_first(center) + 1)) THEN
         CPABORT("Inconsistent transition projector count")
      END IF
      CPASSERT(augmentation%nproj == 0 .OR. augmentation%nproj == projectors)
      shift = 0
      IF (PRESENT(image_shift)) shift = image_shift
      CPASSERT(ALL(shift*left_frame%periodicity == shift))
      position = left_frame%positions(:, atom) + MATMUL(left_frame%cell, REAL(shift, dp))
      left_coefficients_covector = CMPLX(0.0_dp, 0.0_dp, dp)
      right_coefficients_covector = CMPLX(0.0_dp, 0.0_dp, dp)
      ALLOCATE (points(3, nang), potentials(nang, fields%nfunction), weights(nang), weighted_potential(nang), work(nang, nr), &
                gcl(nang, nl), gdl(nang, nl), gcr(nang, nr), gdr(nang, nr))
      radius = MAX(left_frame%delta_cache(kind)%radius, right_frame%delta_cache(kind)%radius)
      CPASSERT(r(SIZE(r)) >= radius)
      DO ir = 1, SIZE(r)
         IF (r(ir) >= radius .OR. radial_weights(ir) == 0.0_dp) CYCLE
         DO ia = 1, nang
            points(:, ia) = position + r(ir)*directions(:, ia)
         END DO
         weights(:) = 0.5_dp*radial_weights(ir)*angular_weights
         CALL residual_fields(left_frame, center, left_range, points, cl, dl, unused_c, unused_d, &
                              unused_points, derivatives=.FALSE., shift=shift, workspace=left_fields)
         CALL residual_fields(right_frame, center, right_range, points, cr, dr, unused_c, unused_d, &
                              unused_points, derivatives=.FALSE., shift=shift, workspace=right_fields)
         DO ia = 1, nang
            IF (PRESENT(auxiliary_images)) THEN
               CALL gaussian_coulomb_fields_periodic(fields, points(:, ia) - auxiliary_center, left_frame%cell, &
                                                     auxiliary_images, potentials(ia, :))
            ELSE
               CALL gaussian_coulomb_fields_evaluate(fields, points(:, ia) - auxiliary_center, potentials(ia, :))
            END IF
         END DO
         gcl = CMPLX(0.0_dp, 0.0_dp, dp)
         gdl = CMPLX(0.0_dp, 0.0_dp, dp)
         gcr = CMPLX(0.0_dp, 0.0_dp, dp)
         gdr = CMPLX(0.0_dp, 0.0_dp, dp)
         DO aux = 1, fields%nfunction
            weighted_potential(:) = weights*potentials(:, aux)
            DO j = 1, nr
               work(:, j) = cr(:, j)*weighted_potential
            END DO
            gdl(:, :) = gdl + MATMUL(work, CONJG(TRANSPOSE(covector(:, :, aux))))
            DO j = 1, nr
               work(:, j) = dr(:, j)*weighted_potential
            END DO
            gcl(:, :) = gcl + MATMUL(work, CONJG(TRANSPOSE(covector(:, :, aux))))
            work(:, :) = MATMUL(cl, covector(:, :, aux))
            DO j = 1, nr
               gdr(:, j) = gdr(:, j) + work(:, j)*weighted_potential
            END DO
            work(:, :) = MATMUL(dl, covector(:, :, aux))
            DO j = 1, nr
               gcr(:, j) = gcr(:, j) + work(:, j)*weighted_potential
            END DO
         END DO
         CALL residual_coefficients_adjoint(left_frame, center, left_fields, gcl, gdl, left_coefficients_covector)
         CALL residual_coefficients_adjoint(right_frame, center, right_fields, gcr, gdr, right_coefficients_covector)
      END DO
      IF (augmentation%nproj == 0) RETURN
      ALLOCATE (q(projectors, projectors))
      qi => temporary_augmentation
      IF (PRESENT(augmentation_integrals)) qi => augmentation_integrals
      CALL prepare_transition_augmentation(fields, augmentation, position - auxiliary_center, directions, angular_weights, &
                                           qi, .FALSE., left_frame%cell, auxiliary_images)
      phase = EXP(CMPLX(0.0_dp, twopi*DOT_PRODUCT(right_frame%reduced_k - left_frame%reduced_k, REAL(shift, dp)), dp))
      p = left_frame%projector_first(center)
      t = right_frame%projector_first(center)
      ASSOCIATE (wl => left_frame%coefficients(p:p + projectors - 1, left_range(1):left_range(2)), &
                 wr => right_frame%coefficients(t:t + projectors - 1, right_range(1):right_range(2)), &
                 gl => left_coefficients_covector(p:p + projectors - 1, :), &
                 gr => right_coefficients_covector(t:t + projectors - 1, :))
         DO aux = 1, fields%nfunction
            q(:, :) = phase*qi%values(:, :, aux)
            gl = gl + MATMUL(MATMUL(q, wr), CONJG(TRANSPOSE(covector(:, :, aux))))
            gr = gr + MATMUL(CONJG(TRANSPOSE(q)), MATMUL(wl, covector(:, :, aux)))
         END DO
      END ASSOCIATE
   END SUBROUTINE qs_upf_inverse_hfx_center_adjoint

! **************************************************************************************************
!> \brief Accumulate coefficient covectors through the complete and single-image difference fields.
!> \param frame ...
!> \param center ...
!> \param workspace Original complete and single-image radial fields from the forward evaluation.
!> \param complete_covector ...
!> \param difference_covector ...
!> \param coefficients_covector ...
! **************************************************************************************************
   SUBROUTINE residual_coefficients_adjoint(frame, center, workspace, complete_covector, difference_covector, &
                                            coefficients_covector)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      INTEGER, INTENT(IN)                                :: center
      TYPE(inverse_field_workspace_type), INTENT(IN)     :: workspace
      COMPLEX(dp), INTENT(IN)                            :: complete_covector(:, :), &
                                                            difference_covector(:, :)
      COMPLEX(dp), INTENT(INOUT)                         :: coefficients_covector(:, :)

      INTEGER                                            :: first, last, other

      ! Use the same radial fields that formed the forward values on this row.
      DO other = 1, SIZE(frame%center_atoms)
         first = frame%projector_first(other)
         last = frame%projector_last(other)
         coefficients_covector(first:last, :) = coefficients_covector(first:last, :) - &
                                                MATMUL(CONJG(TRANSPOSE(workspace%delta(:, first:last))), complete_covector)
         IF (other == center) THEN
            coefficients_covector(first:last, :) = coefficients_covector(first:last, :) - &
                                                   MATMUL(CONJG(TRANSPOSE(workspace%single_delta)), difference_covector)
         END IF
      END DO
   END SUBROUTINE residual_coefficients_adjoint

! **************************************************************************************************
!> \brief Complete soft-plus-Gaussian and single-center difference in a selected native AO block.
!> \param frame ...
!> \param center ...
!> \param ao_range ...
!> \param points ...
!> \param complete ...
!> \param difference ...
!> \param complete_response ...
!> \param difference_response ...
!> \param point_response ...
!> \param atom_response ...
!> \param cell_response ...
!> \param k_response ...
!> \param coefficients_response ...
!> \param derivatives ...
!> \param shift ...
!> \param workspace Scratch arrays reused across the radial rows of this center call.
! **************************************************************************************************
   SUBROUTINE residual_fields(frame, center, ao_range, points, complete, difference, complete_response, &
                              difference_response, point_response, atom_response, cell_response, k_response, &
                              coefficients_response, derivatives, shift, workspace)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      INTEGER, INTENT(IN)                                :: center, ao_range(2)
      REAL(dp), INTENT(IN)                               :: points(:, :)
      COMPLEX(dp), ALLOCATABLE, INTENT(INOUT)            :: complete(:, :), difference(:, :), &
                                                            complete_response(:, :), &
                                                            difference_response(:, :)
      REAL(dp), ALLOCATABLE, INTENT(IN)                  :: point_response(:, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: atom_response(:, :), &
                                                            cell_response(3, 3), k_response(3)
      COMPLEX(dp), INTENT(IN), OPTIONAL                  :: coefficients_response(:, :)
      LOGICAL, INTENT(IN)                                :: derivatives
      INTEGER, INTENT(IN)                                :: shift(3)
      TYPE(inverse_field_workspace_type), INTENT(INOUT)  :: workspace

      INTEGER                                            :: atom, first, kind, last, max_basis, n, &
                                                            np, nproj, other, selected_first, &
                                                            selected_last, single_image(3, 1)

      CPASSERT(frame%atom_local .OR. frame%inverse%ready)
      CPASSERT(SIZE(points, 1) == 3)
      CPASSERT(center >= 1 .AND. center <= SIZE(frame%center_atoms))
      CPASSERT(ao_range(1) >= 1 .AND. ao_range(2) <= SIZE(frame%coefficients, 2))
      CPASSERT(ALL(shift*frame%periodicity == shift))
      np = SIZE(points, 2)
      n = ao_range(2) - ao_range(1) + 1
      nproj = SIZE(frame%coefficients, 1)
      CPASSERT(n > 0)
      IF (.NOT. ALLOCATED(complete)) ALLOCATE (complete(np, n), difference(np, n))
      CPASSERT(ALL(SHAPE(complete) == [np, n]) .AND. ALL(SHAPE(difference) == [np, n]))
      first = frame%projector_first(center)
      last = frame%projector_last(center)
      IF (.NOT. ALLOCATED(workspace%delta)) THEN
         max_basis = MAXVAL(frame%ao_last - frame%ao_first + 1)
         ALLOCATE (workspace%delta(np, nproj), workspace%single_delta(np, last - first + 1), &
                   workspace%phi(np, max_basis))
      END IF
      CPASSERT(ALL(SHAPE(workspace%delta) == [np, nproj]))
      CPASSERT(ALL(SHAPE(workspace%single_delta) == [np, last - first + 1]))
      complete = CMPLX(0.0_dp, 0.0_dp, dp)
      IF (derivatives) THEN
         CPASSERT(ALLOCATED(point_response))
         CPASSERT(PRESENT(atom_response) .AND. PRESENT(cell_response))
         CPASSERT(PRESENT(k_response) .AND. PRESENT(coefficients_response))
         CPASSERT(ALL(SHAPE(point_response) == SHAPE(points)))
         CPASSERT(ALL(SHAPE(atom_response) == SHAPE(frame%positions)))
         CPASSERT(ALL(SHAPE(coefficients_response) == SHAPE(frame%coefficients)))
         IF (.NOT. ALLOCATED(complete_response)) ALLOCATE (complete_response(np, n), difference_response(np, n))
         IF (.NOT. ALLOCATED(workspace%ddelta)) THEN
            ALLOCATE (workspace%ddelta(np, nproj), workspace%single_ddelta(np, last - first + 1), &
                      workspace%dphi(np, SIZE(workspace%phi, 2)))
         END IF
         complete_response = CMPLX(0.0_dp, 0.0_dp, dp)
      END IF
      ! Each intersecting native Gaussian basis is evaluated once, without unused spatial gradients.
      DO atom = 1, SIZE(frame%bases)
         first = frame%ao_first(atom)
         last = frame%ao_last(atom)
         selected_first = MAX(first, ao_range(1))
         selected_last = MIN(last, ao_range(2))
         IF (selected_first > selected_last) CYCLE
         IF (derivatives) THEN
            CALL gaussian_basis_fields_periodic(frame%bases(atom)%gto_basis_set, points, frame%positions(:, atom), &
                                                frame%cell, frame%images, frame%reduced_k, workspace%phi(:, :last - first + 1), &
                                                response=workspace%dphi(:, :last - first + 1), &
                                                point_response=point_response, center_response=atom_response(:, atom), &
                                                cell_response=cell_response, reduced_k_response=k_response)
            complete_response(:, selected_first - ao_range(1) + 1:selected_last - ao_range(1) + 1) = &
               2.0_dp*workspace%dphi(:, selected_first - first + 1:selected_last - first + 1)
         ELSE
            CALL gaussian_basis_fields_periodic(frame%bases(atom)%gto_basis_set, points, frame%positions(:, atom), &
                                                frame%cell, frame%images, frame%reduced_k, workspace%phi(:, :last - first + 1))
         END IF
         complete(:, selected_first - ao_range(1) + 1:selected_last - ao_range(1) + 1) = &
            2.0_dp*workspace%phi(:, selected_first - first + 1:selected_last - first + 1)
      END DO
      ! All original centers and images contribute, including rows outside the selected AO centers.
      DO other = 1, SIZE(frame%center_atoms)
         atom = frame%center_atoms(other)
         kind = frame%center_kinds(other)
         first = frame%projector_first(other)
         last = frame%projector_last(other)
         IF (derivatives) THEN
            CALL upf_radial_fields_periodic(frame%delta_cache(kind), points, frame%positions(:, atom), &
                                            frame%cell, frame%images, frame%reduced_k, workspace%delta(:, first:last), &
                                            response=workspace%ddelta(:, first:last), points_response=point_response, &
                                            center_response=atom_response(:, atom), cell_response=cell_response, &
                                            reduced_k_response=k_response)
         ELSE
            CALL upf_radial_fields_periodic(frame%delta_cache(kind), points, frame%positions(:, atom), &
                                            frame%cell, frame%images, frame%reduced_k, workspace%delta(:, first:last))
         END IF
      END DO
      complete(:, :) = complete - MATMUL(workspace%delta, frame%coefficients(:, ao_range(1):ao_range(2)))
      IF (derivatives) THEN
         complete_response(:, :) = complete_response - &
                                   MATMUL(workspace%ddelta, frame%coefficients(:, ao_range(1):ao_range(2))) - &
                                   MATMUL(workspace%delta, coefficients_response(:, ao_range(1):ao_range(2)))
      END IF
      atom = frame%center_atoms(center)
      kind = frame%center_kinds(center)
      first = frame%projector_first(center)
      last = frame%projector_last(center)
      single_image(:, 1) = shift
      IF (derivatives) THEN
         CALL upf_radial_fields_periodic(frame%delta_cache(kind), points, frame%positions(:, atom), &
                                         frame%cell, single_image, frame%reduced_k, workspace%single_delta, &
                                         response=workspace%single_ddelta, points_response=point_response, &
                                         center_response=atom_response(:, atom), cell_response=cell_response, &
                                         reduced_k_response=k_response)
         difference_response(:, :) = -MATMUL(workspace%single_ddelta, &
                                             frame%coefficients(first:last, ao_range(1):ao_range(2))) - &
                                     MATMUL(workspace%single_delta, &
                                            coefficients_response(first:last, ao_range(1):ao_range(2)))
      ELSE
         CALL upf_radial_fields_periodic(frame%delta_cache(kind), points, frame%positions(:, atom), &
                                         frame%cell, single_image, frame%reduced_k, workspace%single_delta)
      END IF
      difference(:, :) = -MATMUL(workspace%single_delta, frame%coefficients(first:last, ao_range(1):ao_range(2)))
   END SUBROUTINE residual_fields

END MODULE qs_upf_inverse_hfx
