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
   USE upf_inverse_quadrature,          ONLY: upf_inverse_pair,&
                                              upf_inverse_pair_response
   USE upf_radial_fields,               ONLY: upf_radial_fields_periodic
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: qs_upf_inverse_type, qs_upf_inverse_prepare, qs_upf_inverse_from_ks
   PUBLIC :: qs_upf_inverse_response, qs_upf_inverse_release

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
!>        Ghost/NC atoms retain their Gaussian columns; only the selected USPP centers reconstruct.
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
                                                            nao, npoint, nproj, p, q
      TYPE(gto_basis_set_type), POINTER                  :: basis

      CPASSERT(SIZE(points, 1) == 3 .AND. SIZE(points, 2) == SIZE(weights))
      CPASSERT(SIZE(images, 1) == 3 .AND. SIZE(images, 2) > 0)
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
      ALLOCATE (frame%center_atoms, SOURCE=uspp%center_atoms)
      ALLOCATE (frame%center_kinds, SOURCE=uspp%center_kinds)
      ALLOCATE (frame%projector_first(SIZE(uspp%center_atoms)), frame%projector_last(SIZE(uspp%center_atoms)), &
                frame%beta_cache(SIZE(kinds)), frame%delta_cache(SIZE(kinds)))
      nproj = 0
      DO center = 1, SIZE(uspp%center_atoms)
         ikind = uspp%center_kinds(center)
         atom = uspp%center_atoms(center)
         CPASSERT(.NOT. kinds(ikind)%ghost)
         CPASSERT(particles(atom)%atomic_kind%kind_number == ikind)
         IF (frame%beta_cache(ikind)%nproj == 0) THEN
            CALL upf_inverse_fields_init(kinds(ikind)%sgp_potential%upf_data, &
                                         frame%beta_cache(ikind), frame%delta_cache(ikind))
         END IF
         CPASSERT(frame%beta_cache(ikind)%nproj == uspp%projector_sizes(center))
         frame%projector_first(center) = nproj + 1
         nproj = nproj + uspp%projector_sizes(center)
         frame%projector_last(center) = nproj
      END DO
      ALLOCATE (frame%beta(npoint, nproj), frame%delta(npoint, nproj), frame%delta_gradient(npoint, nproj, 3), &
                frame%q(nproj, nproj), frame%radial_q(nproj, nproj), frame%dion(nproj, nproj), frame%coefficients(nproj, nao), &
                coupling(nproj, nproj), overlaps(nproj, nao))
      frame%q(:, :) = CMPLX(0.0_dp, 0.0_dp, dp)
      frame%radial_q(:, :) = CMPLX(0.0_dp, 0.0_dp, dp)
      frame%dion(:, :) = CMPLX(0.0_dp, 0.0_dp, dp)
      DO center = 1, SIZE(uspp%center_atoms)
         ikind = uspp%center_kinds(center)
         atom = uspp%center_atoms(center)
         first = frame%projector_first(center)
         last = frame%projector_last(center)
         frame%radial_q(first:last, first:last) = uspp%augmentation(ikind)%q
         ASSOCIATE (projectors => kinds(ikind)%sgp_potential%upf_projectors, &
                    pot => kinds(ikind)%sgp_potential%upf_data)
            DO q = 1, last - first + 1
               DO p = 1, last - first + 1
                  IF (projectors%angular_index(p) /= projectors%angular_index(q)) CYCLE
                  frame%q(first + p - 1, first + q - 1) = pot%qij(projectors%radial_index(p), projectors%radial_index(q))
               END DO
            END DO
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
      DO atom = 1, SIZE(frame%bases)
         first = frame%ao_first(atom)
         last = frame%ao_last(atom)
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
         CALL upf_radial_fields_periodic(frame%beta_cache(ikind), frame%points, frame%positions(:, atom), &
                                    frame%cell, frame%images, frame%reduced_k, beta(:, first:last), response=dbeta(:, first:last), &
                                         points_response=point_response, center_response=atom_response(:, atom), &
                                         cell_response=cell_response, reduced_k_response=reduced_k_response)
         CALL upf_radial_fields_periodic(frame%delta_cache(ikind), frame%points, frame%positions(:, atom), &
                                         frame%cell, frame%images, frame%reduced_k, delta(:, first:last), gd(:, first:last, :), &
                                         ddelta(:, first:last), dgd(:, first:last, :), point_response, atom_response(:, atom), &
                                         cell_response, reduced_k_response)
      END DO
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
