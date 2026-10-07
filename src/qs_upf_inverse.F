!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Inverse USPP point fields in the actual Quickstep atom/AO and magnetic-projector order.
!>        Point rows are distributed; full AO/projector columns are replicated in this adapter.
! **************************************************************************************************
MODULE qs_upf_inverse
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE cell_types,                      ONLY: cell_type
   USE gaussian_basis_fields,           ONLY: gaussian_basis_fields_periodic
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_comm_type
   USE particle_types,                  ONLY: particle_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type
   USE qs_upf_inverse_types,            ONLY: qs_upf_inverse_type
   USE qs_uspp_types,                   ONLY: qs_uspp_type
   USE upf_inverse,                     ONLY: upf_inverse_basis,&
                                              upf_inverse_basis_response,&
                                              upf_inverse_fields_init,&
                                              upf_inverse_init,&
                                              upf_inverse_project,&
                                              upf_inverse_project_response
   USE upf_inverse_exchange,            ONLY: upf_inverse_exchange_block,&
                                              upf_inverse_exchange_response
   USE upf_inverse_quadrature,          ONLY: upf_inverse_pair,&
                                              upf_inverse_pair_response
   USE upf_radial_fields,               ONLY: upf_radial_fields_periodic
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: qs_upf_inverse_type, qs_upf_inverse_prepare, qs_upf_inverse_from_ks
   PUBLIC :: qs_upf_inverse_response, qs_upf_inverse_release
   PUBLIC :: qs_upf_inverse_evaluate, qs_upf_inverse_center_fields

   PUBLIC :: qs_upf_inverse_transition, qs_upf_inverse_transition_response

CONTAINS

! **************************************************************************************************
!> \brief Read native geometry, ORB bases and original USPP operators from a prepared KS environment.
!>        Does not select an input METHOD or replace the reference SCF's soft-orbital state.
!> \param ks_env Live Quickstep KS environment with prepared USPP state.
!> \param points Local Cartesian quadrature points, indexed (axis,point).
!> \param weights Matching local measure; each point occurs once in comm.
!> \param images Explicit image set, fixed during a geometry response.
!> \param reduced_k Reduced Bloch vector for this frame.
!> \param comm Communicator distributing point rows, not a replicated AO process dimension.
!> \param frame New owned fields; native basis pointers remain borrowed from ks_env.
!> \param info Zero on success, otherwise the inverse-factorization status.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_from_ks(ks_env, points, weights, images, reduced_k, comm, frame, info)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      REAL(dp), INTENT(IN)                               :: points(:, :), weights(:)
      INTEGER, INTENT(IN)                                :: images(:, :)
      REAL(dp), INTENT(IN)                               :: reduced_k(3)
      TYPE(mp_comm_type), INTENT(IN)                     :: comm
      TYPE(qs_upf_inverse_type), INTENT(OUT)             :: frame
      INTEGER, INTENT(OUT)                               :: info

      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_type), POINTER                       :: particles(:)
      TYPE(qs_kind_type), POINTER                        :: kinds(:)
      TYPE(qs_uspp_type), POINTER                        :: uspp

      CALL get_ks_env(ks_env, qs_kind_set=kinds, particle_set=particles, cell=cell, uspp=uspp)
      CPASSERT(ASSOCIATED(kinds) .AND. ASSOCIATED(particles) .AND. ASSOCIATED(cell))
      CPASSERT(ASSOCIATED(uspp))
      CALL qs_upf_inverse_prepare(kinds, particles, cell, uspp, points, weights, images, reduced_k, comm, frame, info)
   END SUBROUTINE qs_upf_inverse_from_ks

! **************************************************************************************************
!> \brief Build the full inverse from native Quickstep basis columns and original radial fields.
!>        Ghost atoms retain their Gaussian columns. NC projectors have zero reconstruction
!>        and act on the same inverse soft fields, including coupling to USPP centers.
!>        Missing original full-wave pairs are rejected by the radial importer.
!> \param kinds Native Quickstep kinds; ORB bases and potential data are not modified.
!> \param particles Native atom order, also used for AO matrix blocks.
!> \param cell Native cell and periodicity flags.
!> \param uspp Prepared center selection and native magnetic augmentation operators.
!> \param points ...
!> \param weights ...
!> \param images ...
!> \param reduced_k ...
!> \param comm ...
!> \param frame ...
!> \param info ...
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_prepare(kinds, particles, cell, uspp, points, weights, images, reduced_k, comm, frame, info)
      TYPE(qs_kind_type), INTENT(IN)                     :: kinds(:)
      TYPE(particle_type), INTENT(IN)                    :: particles(:)
      TYPE(cell_type), INTENT(IN)                        :: cell
      TYPE(qs_uspp_type), INTENT(IN)                     :: uspp
      REAL(dp), INTENT(IN)                               :: points(:, :), weights(:)
      INTEGER, INTENT(IN)                                :: images(:, :)
      REAL(dp), INTENT(IN)                               :: reduced_k(3)
      TYPE(mp_comm_type), INTENT(IN)                     :: comm
      TYPE(qs_upf_inverse_type), INTENT(OUT)             :: frame
      INTEGER, INTENT(OUT)                               :: info

      COMPLEX(dp), ALLOCATABLE                           :: coupling(:, :), overlaps(:, :)
      INTEGER                                            :: atom, axis, center, first, ikind, last, &
                                                            nao, ncenter, npoint, nproj, p, q
      TYPE(gto_basis_set_type), POINTER                  :: basis
      LOGICAL                                            :: nc_kind(SIZE(kinds))

      CPASSERT(SIZE(points, 1) == 3 .AND. SIZE(points, 2) == SIZE(weights))
      CPASSERT(SIZE(images, 1) == 3 .AND. SIZE(images, 2) > 0)
      frame%periodicity = cell%perd
      DO axis = 1, 3
         IF (cell%perd(axis) == 0) THEN
            CPASSERT(ALL(images(axis, :) == 0))
            CPASSERT(reduced_k(axis) == 0.0_dp)
         END IF
      END DO
      CPASSERT(SIZE(uspp%center_atoms) > 0)
      npoint = SIZE(weights)
      frame%point_comm = comm
      frame%cell = cell%hmat
      frame%reduced_k = reduced_k
      ALLOCATE (frame%points, SOURCE=points)
      ALLOCATE (frame%weights, SOURCE=weights)
      ALLOCATE (frame%images, SOURCE=images)
      ALLOCATE (frame%positions(3, SIZE(particles)), frame%bases(SIZE(particles)), &
                frame%ao_first(SIZE(particles)), frame%ao_last(SIZE(particles)))
      nao = 0
      DO atom = 1, SIZE(particles)
         CPASSERT(ASSOCIATED(particles(atom)%atomic_kind))
         ikind = particles(atom)%atomic_kind%kind_number
         CPASSERT(ikind >= 1 .AND. ikind <= SIZE(kinds))
         CALL get_qs_kind(kinds(ikind), basis_set=basis, basis_type='ORB')
         CPASSERT(ASSOCIATED(basis))
         frame%bases(atom)%gto_basis_set => basis
         frame%ao_first(atom) = nao + 1
         nao = nao + basis%nsgf
         frame%ao_last(atom) = nao
         frame%positions(:, atom) = particles(atom)%r
      END DO
      ALLOCATE (frame%gaussian(npoint, nao), frame%gaussian_gradient(npoint, nao, 3), &
                frame%soft(npoint, nao), frame%soft_gradient(npoint, nao, 3))
      DO atom = 1, SIZE(particles)
         first = frame%ao_first(atom)
         last = frame%ao_last(atom)
         CALL gaussian_basis_fields_periodic(frame%bases(atom)%gto_basis_set, points, frame%positions(:, atom), &
                            frame%cell, images, reduced_k, frame%gaussian(:, first:last), frame%gaussian_gradient(:, first:last, :))
      END DO
      nc_kind = .FALSE.
      DO ikind = 1, SIZE(kinds)
         IF (kinds(ikind)%ghost) CYCLE
         IF (.NOT. ASSOCIATED(kinds(ikind)%sgp_potential)) CYCLE
         IF (.NOT. ALLOCATED(kinds(ikind)%sgp_potential%upf_data)) CYCLE
         ASSOCIATE (pot => kinds(ikind)%sgp_potential%upf_data)
            nc_kind(ikind) = TRIM(pot%pseudo_type) == "NC" .AND. &
                             .NOT. pot%is_ultrasoft .AND. .NOT. pot%is_paw .AND. pot%number_of_proj > 0
         END ASSOCIATE
      END DO
      ncenter = SIZE(uspp%center_atoms)
      DO atom = 1, SIZE(particles)
         IF (nc_kind(particles(atom)%atomic_kind%kind_number)) ncenter = ncenter + 1
      END DO
      ALLOCATE (frame%center_atoms(ncenter), frame%center_kinds(ncenter))
      frame%center_atoms(:SIZE(uspp%center_atoms)) = uspp%center_atoms
      frame%center_kinds(:SIZE(uspp%center_atoms)) = uspp%center_kinds
      center = SIZE(uspp%center_atoms)
      DO atom = 1, SIZE(particles)
         ikind = particles(atom)%atomic_kind%kind_number
         IF (.NOT. nc_kind(ikind)) CYCLE
         center = center + 1
         frame%center_atoms(center) = atom
         frame%center_kinds(center) = ikind
      END DO
      ALLOCATE (frame%projector_first(ncenter), frame%projector_last(ncenter), &
                frame%beta_cache(SIZE(kinds)), frame%delta_cache(SIZE(kinds)))
      nproj = 0
      DO center = 1, ncenter
         ikind = frame%center_kinds(center)
         atom = frame%center_atoms(center)
         CPASSERT(.NOT. kinds(ikind)%ghost)
         CPASSERT(particles(atom)%atomic_kind%kind_number == ikind)
         IF (frame%beta_cache(ikind)%nproj == 0) THEN
            CALL upf_inverse_fields_init(kinds(ikind)%sgp_potential%upf_data, &
                                         frame%beta_cache(ikind), frame%delta_cache(ikind))
         END IF
         p = kinds(ikind)%sgp_potential%upf_projectors%nproj
         CPASSERT(frame%beta_cache(ikind)%nproj == p)
         frame%projector_first(center) = nproj + 1
         nproj = nproj + frame%beta_cache(ikind)%nproj
         frame%projector_last(center) = nproj
      END DO
      ALLOCATE (frame%beta(npoint, nproj), frame%delta(npoint, nproj), frame%delta_gradient(npoint, nproj, 3), &
                frame%q(nproj, nproj), frame%radial_q(nproj, nproj), frame%dion(nproj, nproj), frame%coefficients(nproj, nao), &
                coupling(nproj, nproj), overlaps(nproj, nao))
      frame%q(:, :) = CMPLX(0.0_dp, 0.0_dp, dp)
      frame%radial_q(:, :) = CMPLX(0.0_dp, 0.0_dp, dp)
      frame%dion(:, :) = CMPLX(0.0_dp, 0.0_dp, dp)
      DO center = 1, ncenter
         ikind = frame%center_kinds(center)
         atom = frame%center_atoms(center)
         first = frame%projector_first(center)
         last = frame%projector_last(center)
         IF (.NOT. nc_kind(ikind)) frame%radial_q(first:last, first:last) = uspp%augmentation(ikind)%q
         ASSOCIATE (projectors => kinds(ikind)%sgp_potential%upf_projectors, &
                    pot => kinds(ikind)%sgp_potential%upf_data)
            IF (.NOT. nc_kind(ikind)) THEN
               DO q = 1, last - first + 1
                  DO p = 1, last - first + 1
                     IF (projectors%angular_index(p) /= projectors%angular_index(q)) CYCLE
                     frame%q(first + p - 1, first + q - 1) = pot%qij(projectors%radial_index(p), projectors%radial_index(q))
                  END DO
               END DO
            END IF
         END ASSOCIATE
         frame%dion(first:last, first:last) = kinds(ikind)%sgp_potential%upf_projectors%d
         CALL upf_radial_fields_periodic(frame%beta_cache(ikind), points, frame%positions(:, atom), &
                                         frame%cell, images, reduced_k, frame%beta(:, first:last))
         CALL upf_radial_fields_periodic(frame%delta_cache(ikind), points, frame%positions(:, atom), &
                                  frame%cell, images, reduced_k, frame%delta(:, first:last), frame%delta_gradient(:, first:last, :))
      END DO
      CALL upf_inverse_pair(frame%beta, frame%delta, weights, coupling, comm)
      CALL upf_inverse_pair(frame%beta, frame%gaussian, weights, overlaps, comm)
      CALL upf_inverse_init(coupling, frame%inverse, info)
      IF (info /= 0) RETURN
      CALL upf_inverse_project(frame%inverse, overlaps, frame%coefficients)
      CALL upf_inverse_basis(frame%gaussian, frame%delta, frame%coefficients, frame%soft)
      DO axis = 1, 3
         CALL upf_inverse_basis(frame%gaussian_gradient(:, :, axis), frame%delta_gradient(:, :, axis), &
                                frame%coefficients, frame%soft_gradient(:, :, axis))
      END DO
   END SUBROUTINE qs_upf_inverse_prepare

! **************************************************************************************************
!> \brief Evaluate another point batch with the existing inverse transformation.
!>        Point values and gradients use the original frame's W and fixed image set.
!>        A geometry response requires dW from that frame's defining quadrature, rather than
!>        a new inverse fitted to the evaluation points or their integration weights.
!> \param frame Prepared inverse frame, including borrowed live native bases.
!> \param points New local Cartesian points, indexed (axis,point), possibly empty.
!> \param soft Inverse-transformed basis values, overwritten.
!> \param gradient Their physical spatial gradients, overwritten.
!> \param soft_response Optional joint response, required with the other response arguments.
!> \param gradient_response Optional physical-gradient response.
!> \param point_response Optional Cartesian motion of the evaluation points.
!> \param atom_response Optional Cartesian motion of all native atoms.
!> \param cell_response Optional lattice response at fixed reduced Bloch convention.
!> \param reduced_k_response Optional reduced Bloch-vector response.
!> \param coefficients_response Optional dW from the original frame's quadrature response.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_evaluate(frame, points, soft, gradient, soft_response, gradient_response, &
                                      point_response, atom_response, cell_response, reduced_k_response, &
                                      coefficients_response)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      REAL(dp), INTENT(IN)                               :: points(:, :)
      COMPLEX(dp), INTENT(OUT)                           :: soft(:, :), gradient(:, :, :)
      COMPLEX(dp), INTENT(OUT), OPTIONAL                 :: soft_response(:, :), &
                                                            gradient_response(:, :, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: point_response(:, :), &
                                                            atom_response(:, :), &
                                                            cell_response(3, 3), &
                                                            reduced_k_response(3)
      COMPLEX(dp), INTENT(IN), OPTIONAL                  :: coefficients_response(:, :)

      COMPLEX(dp), ALLOCATABLE                           :: chi(:, :), dchi(:, :), ddelta(:, :), &
                                                            delta(:, :), dgc(:, :, :), &
                                                            dgd(:, :, :), gc(:, :, :), gd(:, :, :)
      INTEGER                                            :: atom, axis, center, first, ikind, last, &
                                                            nao, npoint, nproj
      LOGICAL                                            :: calculate_response

      CPASSERT(frame%inverse%ready)
      CPASSERT(SIZE(points, 1) == 3)
      npoint = SIZE(points, 2)
      nao = SIZE(frame%coefficients, 2)
      nproj = SIZE(frame%coefficients, 1)
      CPASSERT(ALL(SHAPE(soft) == [npoint, nao]))
      CPASSERT(ALL(SHAPE(gradient) == [npoint, nao, 3]))
      calculate_response = PRESENT(soft_response) .OR. PRESENT(gradient_response) .OR. &
                           PRESENT(point_response) .OR. PRESENT(atom_response) .OR. &
                           PRESENT(cell_response) .OR. PRESENT(reduced_k_response) .OR. &
                           PRESENT(coefficients_response)
      IF (calculate_response) THEN
         CPASSERT(PRESENT(soft_response) .AND. PRESENT(gradient_response))
         CPASSERT(PRESENT(point_response) .AND. PRESENT(atom_response))
         CPASSERT(PRESENT(cell_response) .AND. PRESENT(reduced_k_response))
         CPASSERT(PRESENT(coefficients_response))
         CPASSERT(ALL(SHAPE(soft_response) == [npoint, nao]))
         CPASSERT(ALL(SHAPE(gradient_response) == [npoint, nao, 3]))
         CPASSERT(ALL(SHAPE(point_response) == SHAPE(points)))
         CPASSERT(ALL(SHAPE(atom_response) == SHAPE(frame%positions)))
         CPASSERT(ALL(SHAPE(coefficients_response) == SHAPE(frame%coefficients)))
      END IF
      ALLOCATE (chi(npoint, nao), gc(npoint, nao, 3), delta(npoint, nproj), gd(npoint, nproj, 3))
      IF (calculate_response) THEN
         ALLOCATE (dchi(npoint, nao), dgc(npoint, nao, 3), ddelta(npoint, nproj), dgd(npoint, nproj, 3))
      END IF
      DO atom = 1, SIZE(frame%bases)
         first = frame%ao_first(atom)
         last = frame%ao_last(atom)
         IF (calculate_response) THEN
            CALL gaussian_basis_fields_periodic(frame%bases(atom)%gto_basis_set, points, frame%positions(:, atom), &
                                                frame%cell, frame%images, frame%reduced_k, chi(:, first:last), &
                                                gc(:, first:last, :), dchi(:, first:last), dgc(:, first:last, :), &
                                                point_response, atom_response(:, atom), cell_response, reduced_k_response)
         ELSE
            CALL gaussian_basis_fields_periodic(frame%bases(atom)%gto_basis_set, points, frame%positions(:, atom), &
                                                frame%cell, frame%images, frame%reduced_k, chi(:, first:last), gc(:, first:last, :))
         END IF
      END DO
      DO center = 1, SIZE(frame%center_atoms)
         atom = frame%center_atoms(center)
         ikind = frame%center_kinds(center)
         first = frame%projector_first(center)
         last = frame%projector_last(center)
         IF (calculate_response) THEN
            CALL upf_radial_fields_periodic(frame%delta_cache(ikind), points, frame%positions(:, atom), &
                                            frame%cell, frame%images, frame%reduced_k, delta(:, first:last), gd(:, first:last, :), &
                                            ddelta(:, first:last), dgd(:, first:last, :), point_response, atom_response(:, atom), &
                                            cell_response, reduced_k_response)
         ELSE
            CALL upf_radial_fields_periodic(frame%delta_cache(ikind), points, frame%positions(:, atom), &
                                            frame%cell, frame%images, frame%reduced_k, delta(:, first:last), gd(:, first:last, :))
         END IF
      END DO
      CALL upf_inverse_basis(chi, delta, frame%coefficients, soft)
      IF (calculate_response) THEN
         CALL upf_inverse_basis_response(dchi, delta, ddelta, frame%coefficients, coefficients_response, soft_response)
      END IF
      DO axis = 1, 3
         CALL upf_inverse_basis(gc(:, :, axis), gd(:, :, axis), frame%coefficients, gradient(:, :, axis))
         IF (calculate_response) THEN
            CALL upf_inverse_basis_response(dgc(:, :, axis), gd(:, :, axis), dgd(:, :, axis), frame%coefficients, &
                                            coefficients_response, gradient_response(:, :, axis))
         END IF
      END DO
   END SUBROUTINE qs_upf_inverse_evaluate

! **************************************************************************************************
!> \brief Complete complex transition integrals between two prepared inverse Bloch frames.
!>        The original inverse coefficients are reused on the auxiliary quadrature. Q integrals
!>        are supplied separately, retaining the original augmentation and its image phases.
!> \param left_frame Original left inverse frame.
!> \param right_frame Independent right frame, possibly at a different Bloch vector.
!> \param points Common local auxiliary integration points.
!> \param weights Original local integration measure.
!> \param potential Complex local auxiliary potentials, with the auxiliary conjugation/phase.
!> \param augmentation Complete replicated original Q integrals in the two projector coordinates.
!> \param integrals New complete rectangular transition block.
!> \param comm Optional communicator partitioning the auxiliary point rows.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_transition(left_frame, right_frame, points, weights, potential, &
                                        augmentation, integrals, comm)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: left_frame, right_frame
      REAL(dp), INTENT(IN)                               :: points(:, :), weights(:)
      COMPLEX(dp), INTENT(IN)                            :: potential(:, :), augmentation(:, :, :)
      COMPLEX(dp), INTENT(OUT)                           :: integrals(:, :, :)
      TYPE(mp_comm_type), INTENT(IN), OPTIONAL           :: comm

      COMPLEX(dp), ALLOCATABLE                           :: left(:, :), lg(:, :, :), rg(:, :, :), &
                                                            right(:, :)
      INTEGER                                            :: nl, npoint, nr

      npoint = SIZE(points, 2)
      nl = SIZE(left_frame%coefficients, 2)
      nr = SIZE(right_frame%coefficients, 2)
      ALLOCATE (left(npoint, nl), right(npoint, nr), lg(npoint, nl, 3), rg(npoint, nr, 3))
      CALL qs_upf_inverse_evaluate(left_frame, points, left, lg)
      CALL qs_upf_inverse_evaluate(right_frame, points, right, rg)
      CALL upf_inverse_exchange_block(left, right, weights, potential, left_frame%coefficients, &
                                      right_frame%coefficients, augmentation, integrals, comm)
   END SUBROUTINE qs_upf_inverse_transition

! **************************************************************************************************
!> \brief Joint geometry/Bloch response of complete inverse transition integrals.
!>        dW must come from each original frame's defining quadrature response, including the
!>        coupled inverse solve. Re-solving on the auxiliary quadrature changes the inverse map.
!> \param left_frame ...
!> \param right_frame ...
!> \param points ...
!> \param weights ...
!> \param potential ...
!> \param augmentation ...
!> \param point_response Motion of the common auxiliary quadrature points.
!> \param weights_response Response of the auxiliary integration measure.
!> \param potential_response Response of the complex auxiliary potential and its image phases.
!> \param augmentation_response Complete replicated original Q integral response.
!> \param atom_response Common nuclear and Gaussian-center motion in native atom order.
!> \param cell_response Common lattice response.
!> \param left_k_response Left reduced Bloch-vector response.
!> \param right_k_response Right reduced Bloch-vector response.
!> \param left_coefficients_response dW from the original left frame.
!> \param right_coefficients_response dW from the original right frame.
!> \param response New complete transition response block.
!> \param comm ...
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_transition_response(left_frame, right_frame, points, weights, potential, &
                                                 augmentation, point_response, weights_response, potential_response, &
                                                 augmentation_response, atom_response, cell_response, left_k_response, &
                                                 right_k_response, left_coefficients_response, &
                                                 right_coefficients_response, response, comm)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: left_frame, right_frame
      REAL(dp), INTENT(IN)                               :: points(:, :), weights(:)
      COMPLEX(dp), INTENT(IN)                            :: potential(:, :), augmentation(:, :, :)
      REAL(dp), INTENT(IN)                               :: point_response(:, :), weights_response(:)
      COMPLEX(dp), INTENT(IN)                            :: potential_response(:, :), &
                                                            augmentation_response(:, :, :)
      REAL(dp), INTENT(IN)                               :: atom_response(:, :), &
                                                            cell_response(3, 3), &
                                                            left_k_response(3), right_k_response(3)
      COMPLEX(dp), INTENT(IN) :: left_coefficients_response(:, :), &
         right_coefficients_response(:, :)
      COMPLEX(dp), INTENT(OUT)                           :: response(:, :, :)
      TYPE(mp_comm_type), INTENT(IN), OPTIONAL           :: comm

      COMPLEX(dp), ALLOCATABLE                           :: dl(:, :), dr(:, :), left(:, :), &
                                                            lg(:, :, :), rg(:, :, :), right(:, :)
      INTEGER                                            :: nl, npoint, nr

      npoint = SIZE(points, 2)
      nl = SIZE(left_frame%coefficients, 2)
      nr = SIZE(right_frame%coefficients, 2)
      ALLOCATE (left(npoint, nl), right(npoint, nr), lg(npoint, nl, 3), rg(npoint, nr, 3), dl(npoint, nl), dr(npoint, nr))
      CALL qs_upf_inverse_evaluate(left_frame, points, left, lg, soft_response=dl, point_response=point_response, &
                                   atom_response=atom_response, cell_response=cell_response, &
                                   reduced_k_response=left_k_response, coefficients_response=left_coefficients_response)
      CALL qs_upf_inverse_evaluate(right_frame, points, right, rg, soft_response=dr, point_response=point_response, &
                                   atom_response=atom_response, cell_response=cell_response, &
                                   reduced_k_response=right_k_response, coefficients_response=right_coefficients_response)
      CALL upf_inverse_exchange_response(left, right, weights, potential, left_frame%coefficients, &
                                         right_frame%coefficients, augmentation, dl, dr, weights_response, potential_response, &
                                         left_coefficients_response, right_coefficients_response, augmentation_response, &
                                         response, comm)
   END SUBROUTINE qs_upf_inverse_transition_response

! **************************************************************************************************
!> \brief Center contribution d_a=-Delta_a*W_a to the inverse-minus-Gaussian field.
!>        W_a is taken from the full coupled solve, never an isolated atomic reconstruction.
!>        Pair this difference with the complete soft-plus-Gaussian field for center residuals.
!>        No spatial partition is applied, including when original wave supports overlap.
!> \param frame Existing inverse and original radial cache.
!> \param center Index in the original frame's center/projector order.
!> \param points Local Cartesian points, possibly an empty batch.
!> \param difference Center difference in the full AO order.
!> \param gradient Physical Cartesian gradient of that center difference.
!> \param difference_response Optional joint field response, with all other response arguments.
!> \param gradient_response Optional physical-gradient response.
!> \param point_response Optional Cartesian evaluation-point motion.
!> \param atom_response Optional Cartesian motion of all native atoms.
!> \param cell_response Optional lattice motion.
!> \param reduced_k_response Optional reduced Bloch-vector motion.
!> \param coefficients_response Optional dW from the original defining quadrature.
!> \param image_shift Optional single lattice image. Use zero for an unpartitioned home-center density.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_center_fields(frame, center, points, difference, gradient, &
                                           difference_response, gradient_response, point_response, &
                                           atom_response, cell_response, reduced_k_response, coefficients_response, image_shift)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      INTEGER, INTENT(IN)                                :: center
      REAL(dp), INTENT(IN)                               :: points(:, :)
      COMPLEX(dp), INTENT(OUT)                           :: difference(:, :), gradient(:, :, :)
      COMPLEX(dp), INTENT(OUT), OPTIONAL                 :: difference_response(:, :), &
                                                            gradient_response(:, :, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: point_response(:, :), &
                                                            atom_response(:, :), &
                                                            cell_response(3, 3), &
                                                            reduced_k_response(3)
      COMPLEX(dp), INTENT(IN), OPTIONAL                  :: coefficients_response(:, :)
      INTEGER, INTENT(IN), OPTIONAL                      :: image_shift(3)

      COMPLEX(dp), ALLOCATABLE                           :: ddelta(:, :), delta(:, :), dgd(:, :, :), &
                                                            gd(:, :, :)
      INTEGER                                            :: atom, axis, first, ikind, last, nao, np, &
                                                            nproj
      INTEGER, ALLOCATABLE                               :: images(:, :)
      LOGICAL                                            :: calculate_response

      CPASSERT(frame%inverse%ready)
      CPASSERT(center >= 1 .AND. center <= SIZE(frame%center_atoms))
      CPASSERT(SIZE(points, 1) == 3)
      IF (PRESENT(image_shift)) THEN
         CPASSERT(ALL(image_shift*frame%periodicity == image_shift))
         ALLOCATE (images(3, 1))
         images(:, 1) = image_shift
      ELSE
         ALLOCATE (images(3, SIZE(frame%images, 2)))
         images(:, :) = frame%images
      END IF
      atom = frame%center_atoms(center)
      ikind = frame%center_kinds(center)
      first = frame%projector_first(center)
      last = frame%projector_last(center)
      np = SIZE(points, 2)
      nao = SIZE(frame%coefficients, 2)
      nproj = last - first + 1
      CPASSERT(ALL(SHAPE(difference) == [np, nao]))
      CPASSERT(ALL(SHAPE(gradient) == [np, nao, 3]))
      calculate_response = PRESENT(difference_response) .OR. PRESENT(gradient_response) .OR. &
                           PRESENT(point_response) .OR. PRESENT(atom_response) .OR. &
                           PRESENT(cell_response) .OR. PRESENT(reduced_k_response) .OR. &
                           PRESENT(coefficients_response)
      ALLOCATE (delta(np, nproj), gd(np, nproj, 3))
      IF (calculate_response) THEN
         CPASSERT(PRESENT(difference_response) .AND. PRESENT(gradient_response))
         CPASSERT(PRESENT(point_response) .AND. PRESENT(atom_response))
         CPASSERT(PRESENT(cell_response) .AND. PRESENT(reduced_k_response))
         CPASSERT(PRESENT(coefficients_response))
         CPASSERT(ALL(SHAPE(difference_response) == [np, nao]))
         CPASSERT(ALL(SHAPE(gradient_response) == [np, nao, 3]))
         CPASSERT(ALL(SHAPE(point_response) == SHAPE(points)))
         CPASSERT(ALL(SHAPE(atom_response) == SHAPE(frame%positions)))
         CPASSERT(ALL(SHAPE(coefficients_response) == SHAPE(frame%coefficients)))
         ALLOCATE (ddelta(np, nproj), dgd(np, nproj, 3))
         CALL upf_radial_fields_periodic(frame%delta_cache(ikind), points, frame%positions(:, atom), &
                                         frame%cell, images, frame%reduced_k, delta, gd, ddelta, dgd, &
                                         point_response, atom_response(:, atom), cell_response, reduced_k_response)
         difference_response = -MATMUL(ddelta, frame%coefficients(first:last, :)) - &
                               MATMUL(delta, coefficients_response(first:last, :))
      ELSE
         CALL upf_radial_fields_periodic(frame%delta_cache(ikind), points, frame%positions(:, atom), &
                                         frame%cell, images, frame%reduced_k, delta, gd)
      END IF
      difference = -MATMUL(delta, frame%coefficients(first:last, :))
      DO axis = 1, 3
         gradient(:, :, axis) = -MATMUL(gd(:, :, axis), frame%coefficients(first:last, :))
         IF (calculate_response) THEN
            gradient_response(:, :, axis) = -MATMUL(dgd(:, :, axis), frame%coefficients(first:last, :)) - &
                                            MATMUL(gd(:, :, axis), coefficients_response(first:last, :))
         END IF
      END DO
   END SUBROUTINE qs_upf_inverse_center_fields

! **************************************************************************************************
!> \brief Joint point/atom/cell/Bloch response in the same native AO and projector frame.
!> \param frame Prepared fields whose native bases remain alive and unchanged.
!> \param point_response Cartesian local point response.
!> \param atom_response Cartesian response of every native atom, including ghost/NC basis centers.
!> \param cell_response Lattice response at the same reduced Bloch convention.
!> \param reduced_k_response Reduced Bloch-vector response.
!> \param weights_response Response of the original local measure.
!> \param coefficients_response Response of W=K^-1 B.
!> \param soft_response Response of the inverse AO values.
!> \param gradient_response Response of their physical spatial gradients.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_response(frame, point_response, atom_response, cell_response, reduced_k_response, &
                                      weights_response, coefficients_response, soft_response, gradient_response)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      REAL(dp), INTENT(IN) :: point_response(:, :), atom_response(:, :), cell_response(3, 3), &
         reduced_k_response(3), weights_response(:)
      COMPLEX(dp), INTENT(OUT)                           :: coefficients_response(:, :), &
                                                            soft_response(:, :), &
                                                            gradient_response(:, :, :)

      COMPLEX(dp), ALLOCATABLE :: beta(:, :), chi(:, :), db(:, :), dbeta(:, :), dchi(:, :), &
         ddelta(:, :), delta(:, :), dgc(:, :, :), dgd(:, :, :), dk(:, :), gc(:, :, :), gd(:, :, :)
      INTEGER                                            :: atom, axis, center, first, ikind, last, &
                                                            nao, npoint, nproj
      LOGICAL                                            :: shared_response

      CPASSERT(frame%inverse%ready)
      CPASSERT(ALL(SHAPE(point_response) == SHAPE(frame%points)))
      CPASSERT(ALL(SHAPE(atom_response) == SHAPE(frame%positions)))
      CPASSERT(SIZE(weights_response) == SIZE(frame%weights))
      npoint = SIZE(frame%weights)
      nao = SIZE(frame%gaussian, 2)
      nproj = frame%inverse%nproj
      CPASSERT(ALL(SHAPE(coefficients_response) == [nproj, nao]))
      CPASSERT(ALL(SHAPE(soft_response) == [npoint, nao]))
      CPASSERT(ALL(SHAPE(gradient_response) == [npoint, nao, 3]))
      ALLOCATE (beta(npoint, nproj), delta(npoint, nproj), dbeta(npoint, nproj), ddelta(npoint, nproj), &
                gd(npoint, nproj, 3), dgd(npoint, nproj, 3), chi(npoint, nao), dchi(npoint, nao), &
                gc(npoint, nao, 3), dgc(npoint, nao, 3), dk(nproj, nproj), db(nproj, nao))
      shared_response = ANY(point_response /= 0.0_dp) .OR. ANY(cell_response /= 0.0_dp) .OR. &
                        ANY(reduced_k_response /= 0.0_dp)
      DO atom = 1, SIZE(frame%bases)
         first = frame%ao_first(atom)
         last = frame%ao_last(atom)
         IF (.NOT. shared_response .AND. ALL(atom_response(:, atom) == 0.0_dp)) THEN
            dchi(:, first:last) = CMPLX(0.0_dp, 0.0_dp, dp)
            dgc(:, first:last, :) = CMPLX(0.0_dp, 0.0_dp, dp)
            CYCLE
         END IF
         CALL gaussian_basis_fields_periodic(frame%bases(atom)%gto_basis_set, frame%points, frame%positions(:, atom), &
                                             frame%cell, frame%images, frame%reduced_k, chi(:, first:last), gc(:, first:last, :), &
                                             dchi(:, first:last), dgc(:, first:last, :), point_response, atom_response(:, atom), &
                                             cell_response, reduced_k_response)
      END DO
      DO center = 1, SIZE(frame%center_atoms)
         atom = frame%center_atoms(center)
         ikind = frame%center_kinds(center)
         first = frame%projector_first(center)
         last = frame%projector_last(center)
         IF (.NOT. shared_response .AND. ALL(atom_response(:, atom) == 0.0_dp)) THEN
            dbeta(:, first:last) = CMPLX(0.0_dp, 0.0_dp, dp)
            ddelta(:, first:last) = CMPLX(0.0_dp, 0.0_dp, dp)
            dgd(:, first:last, :) = CMPLX(0.0_dp, 0.0_dp, dp)
            CYCLE
         END IF
         CALL upf_radial_fields_periodic(frame%beta_cache(ikind), frame%points, frame%positions(:, atom), &
                                    frame%cell, frame%images, frame%reduced_k, beta(:, first:last), response=dbeta(:, first:last), &
                                         points_response=point_response, center_response=atom_response(:, atom), &
                                         cell_response=cell_response, reduced_k_response=reduced_k_response)
         CALL upf_radial_fields_periodic(frame%delta_cache(ikind), frame%points, frame%positions(:, atom), &
                                         frame%cell, frame%images, frame%reduced_k, delta(:, first:last), gd(:, first:last, :), &
                                         ddelta(:, first:last), dgd(:, first:last, :), point_response, atom_response(:, atom), &
                                         cell_response, reduced_k_response)
      END DO
      ! Retain the inverse-coefficient response, which couples all centers through the original fields.
      CALL upf_inverse_pair_response(frame%beta, frame%delta, frame%weights, dbeta, ddelta, &
                                     weights_response, dk, frame%point_comm)
      CALL upf_inverse_pair_response(frame%beta, frame%gaussian, frame%weights, dbeta, dchi, &
                                     weights_response, db, frame%point_comm)
      CALL upf_inverse_project_response(frame%inverse, frame%coefficients, db, dk, coefficients_response)
      CALL upf_inverse_basis_response(dchi, frame%delta, ddelta, frame%coefficients, &
                                      coefficients_response, soft_response)
      DO axis = 1, 3
         CALL upf_inverse_basis_response(dgc(:, :, axis), frame%delta_gradient(:, :, axis), dgd(:, :, axis), &
                                         frame%coefficients, coefficients_response, gradient_response(:, :, axis))
      END DO
   END SUBROUTINE qs_upf_inverse_response

! **************************************************************************************************
!> \brief Release owned fields and borrowed basis pointers before destroying their KS environment.
!> \param frame ...
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_release(frame)
      TYPE(qs_upf_inverse_type), INTENT(INOUT)           :: frame

      TYPE(qs_upf_inverse_type)                          :: empty

      frame = empty
   END SUBROUTINE qs_upf_inverse_release

END MODULE qs_upf_inverse
