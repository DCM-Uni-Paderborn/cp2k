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
   USE mathconstants,                   ONLY: twopi
   USE message_passing,                 ONLY: mp_comm_type
   USE orbital_pointers,                ONLY: coset,&
                                              ncoset
   USE particle_types,                  ONLY: particle_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type
   USE qs_upf_inverse_types,            ONLY: qs_upf_inverse_type
   USE qs_uspp_types,                   ONLY: qs_uspp_type
   USE upf_inverse,                     ONLY: &
        upf_inverse_atomic_init, upf_inverse_atomic_residual, upf_inverse_basis, &
        upf_inverse_basis_response, upf_inverse_fields_init, upf_inverse_init, &
        upf_inverse_operator_response, upf_inverse_project, upf_inverse_project_adjoint, &
        upf_inverse_project_response
   USE upf_inverse_exchange,            ONLY: upf_inverse_exchange_block,&
                                              upf_inverse_exchange_contract,&
                                              upf_inverse_exchange_response
   USE upf_inverse_quadrature,          ONLY: upf_inverse_pair,&
                                              upf_inverse_pair_response
   USE upf_projector_integrals,         ONLY: upf_projector_overlap,&
                                              upf_projector_type
   USE upf_radial_fields,               ONLY: upf_radial_fields_periodic
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: qs_upf_inverse_type, qs_upf_inverse_prepare, qs_upf_inverse_from_ks
   PUBLIC :: qs_upf_inverse_response, qs_upf_inverse_release
   PUBLIC :: qs_upf_inverse_projection_adjoint, qs_upf_inverse_projection_contract
   PUBLIC :: qs_upf_inverse_evaluate, qs_upf_inverse_center_fields
   PUBLIC :: qs_upf_inverse_atomic_overlap
   PUBLIC :: qs_upf_inverse_atomic_operators
   PUBLIC :: qs_upf_inverse_transition, qs_upf_inverse_transition_response
   PUBLIC :: qs_upf_inverse_exchange_pair, qs_upf_inverse_exchange_pair_response

CONTAINS

! **************************************************************************************************
!> \brief Contract atom-local hard projectors with native Gaussian AOs and explicit Bloch images.
!>        Only the original radial quadrature is numerical. Gaussian Cartesian integrals and
!>        AO-center derivatives use the existing UPF integrator, without a uniform-grid solve.
!> \param basis Native contracted ORB basis and unchanged primitive normalization.
!> \param hard Atom-local transformed projectors, including enough derivative degree.
!> \param center Home projector center.
!> \param ao_center Home Gaussian center.
!> \param cell Cartesian lattice vectors.
!> \param images Fixed integer Gaussian images.
!> \param reduced_k Bloch vector, with phase exp(i*2*pi*k.image).
!> \param overlaps Hard projector bras against Gaussian AO columns, without occupations.
!> \param response Optional joint first response.
!> \param center_response Optional projector-center displacement.
!> \param ao_response Optional Gaussian-center displacement.
!> \param cell_response Optional lattice displacement.
!> \param reduced_k_response Optional reduced Bloch-vector displacement.
!> \param kinetic Optional projector overlaps with the Gaussian kinetic operator.
!> \param kinetic_response Optional joint first response of the kinetic overlaps.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_atomic_overlap(basis, hard, center, ao_center, cell, images, reduced_k, overlaps, &
                                            response, center_response, ao_response, cell_response, reduced_k_response, &
                                            kinetic, kinetic_response)
      TYPE(gto_basis_set_type), INTENT(IN)               :: basis
      TYPE(upf_projector_type), INTENT(IN)               :: hard
      REAL(dp), INTENT(IN)                               :: center(3), ao_center(3), cell(3, 3)
      INTEGER, INTENT(IN)                                :: images(:, :)
      REAL(dp), INTENT(IN)                               :: reduced_k(3)
      COMPLEX(dp), INTENT(OUT)                           :: overlaps(:, :)
      COMPLEX(dp), INTENT(OUT), OPTIONAL                 :: response(:, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: center_response(3), ao_response(3), &
                                                            cell_response(3, 3), &
                                                            reduced_k_response(3)
      COMPLEX(dp), INTENT(OUT), OPTIONAL                 :: kinetic(:, :), kinetic_response(:, :)

      COMPLEX(dp)                                        :: phase, phase_response
      INTEGER                                            :: axis, derivative_index, first, image, &
                                                            iset, last, nder, np, nshell, powers(3)
      REAL(dp)                                           :: direction(3), displacement(3), dk(3), &
                                                            shift(3)
      REAL(dp), ALLOCATABLE                              :: contracted(:, :), derivative(:, :), &
                                                            integrals(:, :, :)

      CPASSERT(SIZE(images, 1) == 3)
      CPASSERT(ALL(SHAPE(overlaps) == [hard%nproj, basis%nsgf]))
      nder = 0
      IF (PRESENT(response)) THEN
         CPASSERT(ALL(SHAPE(response) == SHAPE(overlaps)))
         response = CMPLX(0.0_dp, 0.0_dp, dp)
         nder = 1
      END IF
      IF (PRESENT(kinetic)) THEN
         CPASSERT(ALL(SHAPE(kinetic) == SHAPE(overlaps)))
         kinetic = CMPLX(0.0_dp, 0.0_dp, dp)
         nder = MAX(nder, 2)
      END IF
      IF (PRESENT(kinetic_response)) THEN
         CPASSERT(PRESENT(kinetic))
         CPASSERT(ALL(SHAPE(kinetic_response) == SHAPE(overlaps)))
         kinetic_response = CMPLX(0.0_dp, 0.0_dp, dp)
         nder = 3
      END IF
      overlaps = CMPLX(0.0_dp, 0.0_dp, dp)
      dk = 0.0_dp
      IF (PRESENT(reduced_k_response)) dk = reduced_k_response
      DO image = 1, SIZE(images, 2)
         shift = MATMUL(cell, REAL(images(:, image), dp))
         displacement = center - ao_center - shift
         direction = 0.0_dp
         IF (PRESENT(ao_response)) direction = direction + ao_response
         IF (PRESENT(center_response)) direction = direction - center_response
         IF (PRESENT(cell_response)) direction = direction + MATMUL(cell_response, REAL(images(:, image), dp))
         phase = EXP(CMPLX(0.0_dp, twopi*DOT_PRODUCT(reduced_k, REAL(images(:, image), dp)), dp))
         phase_response = CMPLX(0.0_dp, twopi*DOT_PRODUCT(dk, REAL(images(:, image), dp)), dp)*phase
         DO iset = 1, basis%nset
            nshell = basis%nshell(iset)
            first = MINVAL(basis%first_sgf(:nshell, iset))
            last = MAXVAL(basis%last_sgf(:nshell, iset))
            np = basis%npgf(iset)*ncoset(basis%lmax(iset))
            ALLOCATE (integrals(np, hard%nproj, ncoset(nder)), contracted(last - first + 1, hard%nproj))
            CALL upf_projector_overlap(hard, basis%lmax(iset), basis%zet(:basis%npgf(iset), iset), &
                                       displacement, nder, integrals)
            contracted(:, :) = MATMUL(TRANSPOSE(basis%sphi(:np, first:last)), integrals(:, :, 1))
            overlaps(:, first:last) = overlaps(:, first:last) + phase*TRANSPOSE(contracted)
            IF (PRESENT(response)) THEN
               response(:, first:last) = response(:, first:last) + phase_response*TRANSPOSE(contracted)
               ALLOCATE (derivative(last - first + 1, hard%nproj))
               derivative = 0.0_dp
               DO axis = 1, 3
                  derivative(:, :) = derivative + direction(axis)* &
                                     MATMUL(TRANSPOSE(basis%sphi(:np, first:last)), integrals(:, :, axis + 1))
               END DO
               response(:, first:last) = response(:, first:last) + phase*TRANSPOSE(derivative)
               DEALLOCATE (derivative)
            END IF
            IF (PRESENT(kinetic)) THEN
               contracted = 0.0_dp
               DO axis = 1, 3
                  powers = 0
                  powers(axis) = 2
                  derivative_index = coset(powers(1), powers(2), powers(3))
                  contracted(:, :) = contracted - 0.5_dp* &
                                     MATMUL(TRANSPOSE(basis%sphi(:np, first:last)), integrals(:, :, derivative_index))
               END DO
               kinetic(:, first:last) = kinetic(:, first:last) + phase*TRANSPOSE(contracted)
               IF (PRESENT(kinetic_response)) THEN
                  kinetic_response(:, first:last) = kinetic_response(:, first:last) + phase_response*TRANSPOSE(contracted)
                  ALLOCATE (derivative(last - first + 1, hard%nproj))
                  derivative = 0.0_dp
                  DO axis = 1, 3
                     DO derivative_index = 1, 3
                        powers = 0
                        powers(axis) = 2
                        powers(derivative_index) = powers(derivative_index) + 1
                        derivative(:, :) = derivative - 0.5_dp*direction(derivative_index)* &
                                           MATMUL(TRANSPOSE(basis%sphi(:np, first:last)), &
                                                  integrals(:, :, coset(powers(1), powers(2), powers(3))))
                     END DO
                  END DO
                  kinetic_response(:, first:last) = kinetic_response(:, first:last) + phase*TRANSPOSE(derivative)
                  DEALLOCATE (derivative)
               END IF
            END IF
            DEALLOCATE (integrals, contracted)
         END DO
      END DO
   END SUBROUTINE qs_upf_inverse_atomic_overlap

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
!> \param atom_local Optional nonoverlapping one-center inverse; the general reference is the default.
!> \param direct Use original projectors and soft orbitals without reconstruction.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_from_ks(ks_env, points, weights, images, reduced_k, comm, frame, info, atom_local, direct)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      REAL(dp), INTENT(IN)                               :: points(:, :), weights(:)
      INTEGER, INTENT(IN)                                :: images(:, :)
      REAL(dp), INTENT(IN)                               :: reduced_k(3)
      TYPE(mp_comm_type), INTENT(IN)                     :: comm
      TYPE(qs_upf_inverse_type), INTENT(OUT)             :: frame
      INTEGER, INTENT(OUT)                               :: info
      LOGICAL, INTENT(IN), OPTIONAL                      :: atom_local, direct

      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_type), POINTER                       :: particles(:)
      TYPE(qs_kind_type), POINTER                        :: kinds(:)
      TYPE(qs_uspp_type), POINTER                        :: uspp

      CALL get_ks_env(ks_env, qs_kind_set=kinds, particle_set=particles, cell=cell, uspp=uspp)
      CPASSERT(ASSOCIATED(kinds) .AND. ASSOCIATED(particles) .AND. ASSOCIATED(cell))
      CPASSERT(ASSOCIATED(uspp))
      CALL qs_upf_inverse_prepare(kinds, particles, cell, uspp, points, weights, images, reduced_k, &
                                  comm, frame, info, atom_local, direct)
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
!> \param atom_local Optional nonoverlapping one-center inverse. Status -2 reports overlapping support.
!> \param direct Use original projectors and soft orbitals without reconstruction.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_prepare(kinds, particles, cell, uspp, points, weights, images, reduced_k, &
                                     comm, frame, info, atom_local, direct)
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
      LOGICAL, INTENT(IN), OPTIONAL                      :: atom_local, direct

      COMPLEX(dp), ALLOCATABLE                           :: coupling(:, :), overlaps(:, :)
      INTEGER                                            :: atom, axis, center, degree, first, &
                                                            ikind, last, nao, ncenter, npoint, &
                                                            nproj, p, q
      LOGICAL                                            :: use_atomic
      TYPE(gto_basis_set_type), POINTER                  :: basis
      LOGICAL                                            :: nc_kind(SIZE(kinds))

      CPASSERT(SIZE(points, 1) == 3 .AND. SIZE(points, 2) == SIZE(weights))
      CPASSERT(SIZE(images, 1) == 3 .AND. SIZE(images, 2) > 0)
      use_atomic = .FALSE.
      IF (PRESENT(atom_local)) use_atomic = atom_local
      IF (PRESENT(direct)) frame%direct = direct
      IF (frame%direct) use_atomic = .TRUE.
      info = 0
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
      degree = 0
      DO atom = 1, SIZE(particles)
         CPASSERT(ASSOCIATED(particles(atom)%atomic_kind))
         ikind = particles(atom)%atomic_kind%kind_number
         CPASSERT(ikind >= 1 .AND. ikind <= SIZE(kinds))
         CALL get_qs_kind(kinds(ikind), basis_set=basis, basis_type='ORB')
         CPASSERT(ASSOCIATED(basis))
         degree = MAX(degree, MAXVAL(basis%lmax) + 3)
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
      IF (use_atomic) ALLOCATE (frame%atomic_cache(SIZE(kinds)))
      nproj = 0
      DO center = 1, ncenter
         ikind = frame%center_kinds(center)
         atom = frame%center_atoms(center)
         CPASSERT(.NOT. kinds(ikind)%ghost)
         CPASSERT(particles(atom)%atomic_kind%kind_number == ikind)
         IF (frame%beta_cache(ikind)%nproj == 0) THEN
            CALL upf_inverse_fields_init(kinds(ikind)%sgp_potential%upf_data, &
                                         frame%beta_cache(ikind), frame%delta_cache(ikind), frame%direct)
            IF (use_atomic) THEN
               CALL upf_inverse_atomic_init(kinds(ikind)%sgp_potential%upf_data, degree, &
                                            frame%atomic_cache(ikind), info, frame%direct)
               IF (info /= 0) RETURN
            END IF
         END IF
         p = kinds(ikind)%sgp_potential%upf_projectors%nproj
         CPASSERT(frame%beta_cache(ikind)%nproj == p)
         frame%projector_first(center) = nproj + 1
         nproj = nproj + frame%beta_cache(ikind)%nproj
         frame%projector_last(center) = nproj
      END DO
      ALLOCATE (frame%beta(npoint, nproj), frame%delta(npoint, nproj), frame%delta_gradient(npoint, nproj, 3), &
                frame%q(nproj, nproj), frame%radial_q(nproj, nproj), frame%dion(nproj, nproj), frame%coefficients(nproj, nao))
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
      IF (use_atomic) THEN
         IF (.NOT. frame%direct) CALL check_atomic_domains(frame, info)
         IF (info /= 0) RETURN
         ALLOCATE (overlaps(nproj, nao))
         CALL atomic_coefficients(frame, overlaps)
         frame%coefficients(:, :) = overlaps
         frame%atom_local = .TRUE.
      ELSE
         ALLOCATE (coupling(nproj, nproj), overlaps(nproj, nao))
         CALL upf_inverse_pair(frame%beta, frame%delta, weights, coupling, comm)
         CALL upf_inverse_pair(frame%beta, frame%gaussian, weights, overlaps, comm)
         CALL upf_inverse_init(coupling, frame%inverse, info)
         IF (info /= 0) RETURN
         CALL upf_inverse_project(frame%inverse, overlaps, frame%coefficients)
      END IF
      CALL upf_inverse_basis(frame%gaussian, frame%delta, frame%coefficients, frame%soft)
      DO axis = 1, 3
         CALL upf_inverse_basis(frame%gaussian_gradient(:, :, axis), frame%delta_gradient(:, :, axis), &
                                frame%coefficients, frame%soft_gradient(:, :, axis))
      END DO
   END SUBROUTINE qs_upf_inverse_prepare

! **************************************************************************************************
!> \brief Reject cross-center beta-Delta and Delta-Delta support in the local inverse.
!> \param frame Original native centers, fields and explicitly supplied periodic images.
!> \param info Zero for nonoverlap, -2 otherwise. No overlap term is silently discarded.
! **************************************************************************************************
   SUBROUTINE check_atomic_domains(frame, info)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      INTEGER, INTENT(OUT)                               :: info

      INTEGER                                            :: a, b, image, ka, kb
      REAL(dp)                                           :: displacement(3), ra, radius, rb

      info = 0
      DO a = 1, SIZE(frame%center_atoms)
         ka = frame%center_kinds(a)
         ra = frame%delta_cache(ka)%radius
         DO b = 1, SIZE(frame%center_atoms)
            kb = frame%center_kinds(b)
            rb = frame%delta_cache(kb)%radius
            IF (ra == 0.0_dp .AND. rb == 0.0_dp) CYCLE
            radius = ra + rb
            IF (rb > 0.0_dp) radius = MAX(radius, frame%beta_cache(ka)%radius + rb)
            IF (ra > 0.0_dp) radius = MAX(radius, ra + frame%beta_cache(kb)%radius)
            DO image = 1, SIZE(frame%images, 2)
               IF (a == b .AND. ALL(frame%images(:, image) == 0)) CYCLE
               displacement = frame%positions(:, frame%center_atoms(a)) - frame%positions(:, frame%center_atoms(b)) - &
                              MATMUL(frame%cell, REAL(frame%images(:, image), dp))
               IF (SUM(displacement**2) >= radius**2) CYCLE
               info = -2
               RETURN
            END DO
         END DO
      END DO
   END SUBROUTINE check_atomic_domains

! **************************************************************************************************
!> \brief Native atom-local hard projector coefficients and their common geometry response.
!>        The radial K is fixed. No grid measure, point motion or global inverse enters W.
!> \param frame Prepared original kinds, Gaussian bases, centers and physical-cell images.
!> \param overlaps New hard projector bras against all Gaussian columns.
!> \param response Optional joint response, required with all three geometry arguments.
!> \param atom_response Optional nuclear and Gaussian-center motion in the native atom order.
!> \param cell_response Optional lattice motion.
!> \param reduced_k_response Optional reduced Bloch-vector motion.
! **************************************************************************************************
   SUBROUTINE atomic_coefficients(frame, overlaps, response, atom_response, cell_response, reduced_k_response)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      COMPLEX(dp), INTENT(OUT)                           :: overlaps(:, :)
      COMPLEX(dp), INTENT(OUT), OPTIONAL                 :: response(:, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: atom_response(:, :), &
                                                            cell_response(3, 3), &
                                                            reduced_k_response(3)

      INTEGER                                            :: atom, center, first, ikind, last, p, q

      CPASSERT(ALLOCATED(frame%atomic_cache))
      IF (PRESENT(response)) THEN
         CPASSERT(PRESENT(atom_response))
         CPASSERT(PRESENT(cell_response) .AND. PRESENT(reduced_k_response))
      END IF
      DO center = 1, SIZE(frame%center_atoms)
         ikind = frame%center_kinds(center)
         atom = frame%center_atoms(center)
         p = frame%projector_first(center)
         q = frame%projector_last(center)
         DO first = 1, SIZE(frame%bases)
            last = frame%ao_last(first)
            IF (PRESENT(response)) THEN
               CALL qs_upf_inverse_atomic_overlap(frame%bases(first)%gto_basis_set, frame%atomic_cache(ikind)%hard, &
                                                  frame%positions(:, atom), frame%positions(:, first), frame%cell, frame%images, &
                                                  frame%reduced_k, overlaps(p:q, frame%ao_first(first):last), &
                                                  response(p:q, frame%ao_first(first):last), &
                                                 atom_response(:, atom), atom_response(:, first), cell_response, reduced_k_response)
            ELSE
               CALL qs_upf_inverse_atomic_overlap(frame%bases(first)%gto_basis_set, frame%atomic_cache(ikind)%hard, &
                                                  frame%positions(:, atom), frame%positions(:, first), frame%cell, frame%images, &
                                                  frame%reduced_k, overlaps(p:q, frame%ao_first(first):last))
            END IF
         END DO
      END DO
   END SUBROUTINE atomic_coefficients

! **************************************************************************************************
!> \brief Atomic corrections to native Gaussian kinetic and overlap matrices, retaining closure.
!>        The kinetic correction equals -W^dagger T W plus its finite-map residual.
!>        Adding the transformed original D gives the intended hard D-T representation.
!> \param frame Nonoverlapping atom-local inverse with native hard projector coefficients.
!> \param projector_metric Original metric selected by the functional, without substituting q_wave.
!> \param kinetic Owned correction to the ordinary Gaussian kinetic matrix.
!> \param metric Owned correction to the ordinary Gaussian overlap matrix.
!> \param kinetic_response Optional owned joint kinetic response, with metric_response and geometry.
!> \param metric_response Optional owned joint metric response.
!> \param atom_response Optional nuclear and Gaussian-center displacement.
!> \param cell_response Optional lattice displacement.
!> \param reduced_k_response Optional reduced Bloch-vector displacement.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_atomic_operators(frame, projector_metric, kinetic, metric, kinetic_response, metric_response, &
                                              atom_response, cell_response, reduced_k_response)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      COMPLEX(dp), INTENT(IN)                            :: projector_metric(:, :)
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: kinetic(:, :), metric(:, :)
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT), OPTIONAL    :: kinetic_response(:, :), &
                                                            metric_response(:, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: atom_response(:, :), &
                                                            cell_response(3, 3), &
                                                            reduced_k_response(3)

      COMPLEX(dp), ALLOCATABLE :: dl(:, :), dw(:, :), dz(:, :), laplacian(:, :), overlaps(:, :), &
         rs(:, :), rt(:, :), z(:, :), zero_ao(:, :), zero_projectors(:, :)
      INTEGER                                            :: atom, b, center, first, ikind, last, n, &
                                                            np, p, q
      LOGICAL                                            :: calculate_response

      CPASSERT(frame%atom_local .AND. ALLOCATED(frame%atomic_cache))
      n = SIZE(frame%coefficients, 2)
      np = SIZE(frame%coefficients, 1)
      CPASSERT(ALL(SHAPE(projector_metric) == [np, np]))
      ALLOCATE (kinetic(n, n), metric(n, n))
      kinetic = CMPLX(0.0_dp, 0.0_dp, dp)
      metric = CMPLX(0.0_dp, 0.0_dp, dp)
      calculate_response = PRESENT(kinetic_response) .OR. PRESENT(metric_response)
      IF (calculate_response) THEN
         CPASSERT(PRESENT(kinetic_response))
         CPASSERT(PRESENT(metric_response))
         CPASSERT(PRESENT(atom_response))
         CPASSERT(PRESENT(cell_response) .AND. PRESENT(reduced_k_response))
         ALLOCATE (kinetic_response(n, n), metric_response(n, n))
         kinetic_response = CMPLX(0.0_dp, 0.0_dp, dp)
         metric_response = CMPLX(0.0_dp, 0.0_dp, dp)
      END IF
      IF (frame%direct) RETURN
      IF (calculate_response) THEN
         ALLOCATE (dw(np, n), overlaps(np, n), zero_ao(n, n))
         zero_ao = CMPLX(0.0_dp, 0.0_dp, dp)
         CALL atomic_coefficients(frame, overlaps, dw, atom_response, cell_response, reduced_k_response)
      END IF
      DO center = 1, SIZE(frame%center_atoms)
         ikind = frame%center_kinds(center)
         atom = frame%center_atoms(center)
         p = frame%projector_first(center)
         q = frame%projector_last(center)
         ALLOCATE (z(q - p + 1, n), laplacian(q - p + 1, n))
         IF (calculate_response) THEN
            ALLOCATE (dz(q - p + 1, n), dl(q - p + 1, n), zero_projectors(q - p + 1, q - p + 1))
            zero_projectors = CMPLX(0.0_dp, 0.0_dp, dp)
         END IF
         DO b = 1, SIZE(frame%bases)
            first = frame%ao_first(b)
            last = frame%ao_last(b)
            IF (calculate_response) THEN
               CALL qs_upf_inverse_atomic_overlap(frame%bases(b)%gto_basis_set, frame%atomic_cache(ikind)%difference, &
                                                  frame%positions(:, atom), frame%positions(:, b), frame%cell, frame%images, &
                                                  frame%reduced_k, z(:, first:last), dz(:, first:last), &
                                                  atom_response(:, atom), atom_response(:, b), cell_response, reduced_k_response, &
                                                  laplacian(:, first:last), dl(:, first:last))
            ELSE
               CALL qs_upf_inverse_atomic_overlap(frame%bases(b)%gto_basis_set, frame%atomic_cache(ikind)%difference, &
                                                  frame%positions(:, atom), frame%positions(:, b), frame%cell, frame%images, &
                                                  frame%reduced_k, z(:, first:last), kinetic=laplacian(:, first:last))
            END IF
         END DO
         ASSOCIATE (atomic => frame%atomic_cache(ikind), w => frame%coefficients(p:q, :))
            CALL upf_inverse_atomic_residual(atomic, w, z, laplacian, rt, rs)
            kinetic(:, :) = kinetic + rt - MATMUL(CONJG(TRANSPOSE(w)), MATMUL(atomic%kinetic, w))
            metric(:, :) = metric + rs + MATMUL(CONJG(TRANSPOSE(w)), &
                                                MATMUL(projector_metric(p:q, p:q) - atomic%wave_overlap, w))
            IF (calculate_response) THEN
               CALL upf_inverse_operator_response(zero_ao, CONJG(TRANSPOSE(laplacian)), laplacian, &
                                                  CMPLX(atomic%difference_kinetic, 0.0_dp, dp), zero_projectors, w, &
                                                  CONJG(TRANSPOSE(dl)), dl, zero_projectors, zero_projectors, dw(p:q, :), rt)
               CALL upf_inverse_operator_response(zero_ao, CONJG(TRANSPOSE(z)), z, &
                                                  CMPLX(atomic%difference_overlap, 0.0_dp, dp), projector_metric(p:q, p:q), w, &
                                                  CONJG(TRANSPOSE(dz)), dz, zero_projectors, zero_projectors, dw(p:q, :), rs)
               kinetic_response(:, :) = kinetic_response + rt
               metric_response(:, :) = metric_response + rs
            END IF
         END ASSOCIATE
         DEALLOCATE (z, laplacian, rt, rs)
         IF (calculate_response) DEALLOCATE (dz, dl, zero_projectors)
      END DO
   END SUBROUTINE qs_upf_inverse_atomic_operators

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

      CPASSERT(frame%atom_local .OR. frame%inverse%ready)
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
         CPASSERT(PRESENT(soft_response))
         CPASSERT(PRESENT(point_response) .AND. PRESENT(atom_response))
         CPASSERT(PRESENT(cell_response) .AND. PRESENT(reduced_k_response))
         CPASSERT(PRESENT(coefficients_response))
         CPASSERT(ALL(SHAPE(soft_response) == [npoint, nao]))
         IF (PRESENT(gradient_response)) THEN
            CPASSERT(ALL(SHAPE(gradient_response) == [npoint, nao, 3]))
         END IF
         CPASSERT(ALL(SHAPE(point_response) == SHAPE(points)))
         CPASSERT(ALL(SHAPE(atom_response) == SHAPE(frame%positions)))
         CPASSERT(ALL(SHAPE(coefficients_response) == SHAPE(frame%coefficients)))
      END IF
      ALLOCATE (chi(npoint, nao), gc(npoint, nao, 3), delta(npoint, nproj), gd(npoint, nproj, 3))
      IF (calculate_response) THEN
         ALLOCATE (dchi(npoint, nao), ddelta(npoint, nproj))
         IF (PRESENT(gradient_response)) ALLOCATE (dgc(npoint, nao, 3), dgd(npoint, nproj, 3))
      END IF
      DO atom = 1, SIZE(frame%bases)
         first = frame%ao_first(atom)
         last = frame%ao_last(atom)
         IF (PRESENT(gradient_response)) THEN
            CALL gaussian_basis_fields_periodic(frame%bases(atom)%gto_basis_set, points, frame%positions(:, atom), &
                                                frame%cell, frame%images, frame%reduced_k, chi(:, first:last), &
                                                gc(:, first:last, :), dchi(:, first:last), dgc(:, first:last, :), &
                                                point_response, atom_response(:, atom), cell_response, reduced_k_response)
         ELSE IF (calculate_response) THEN
            CALL gaussian_basis_fields_periodic(frame%bases(atom)%gto_basis_set, points, frame%positions(:, atom), &
                                                frame%cell, frame%images, frame%reduced_k, chi(:, first:last), &
                                                gc(:, first:last, :), response=dchi(:, first:last), &
                                                point_response=point_response, center_response=atom_response(:, atom), &
                                                cell_response=cell_response, reduced_k_response=reduced_k_response)
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
         IF (PRESENT(gradient_response)) THEN
            CALL upf_radial_fields_periodic(frame%delta_cache(ikind), points, frame%positions(:, atom), &
                                            frame%cell, frame%images, frame%reduced_k, delta(:, first:last), gd(:, first:last, :), &
                                            ddelta(:, first:last), dgd(:, first:last, :), point_response, atom_response(:, atom), &
                                            cell_response, reduced_k_response)
         ELSE IF (calculate_response) THEN
            CALL upf_radial_fields_periodic(frame%delta_cache(ikind), points, frame%positions(:, atom), &
                                            frame%cell, frame%images, frame%reduced_k, delta(:, first:last), gd(:, first:last, :), &
                                            response=ddelta(:, first:last), points_response=point_response, &
                                            center_response=atom_response(:, atom), cell_response=cell_response, &
                                            reduced_k_response=reduced_k_response)
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
         IF (PRESENT(gradient_response)) THEN
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
!> \brief Ordered physical-spin exchange pair from complete native inverse-frame transitions.
!>        Original W is retained on the auxiliary quadrature. The RI owner supplies its complete
!>        effective interaction, including any local auxiliary selection and metric weighting.
!>        Point rows are distributed; the reduced transition, densities, energy and operators
!>        are replicated. The owner counts each ordered spin/k pair once.
!> \param left_frame ...
!> \param right_frame ...
!> \param points ...
!> \param weights ...
!> \param potential ...
!> \param augmentation ...
!> \param interaction Complete Hermitian effective RI kernel in the auxiliary order.
!> \param left_density Physical occupation matrix for one spin, without its k weight.
!> \param right_density Independent right physical occupation matrix of the same spin.
!> \param prefactor Minus one-half times exchange fraction and both Bloch weights.
!> \param energy New ordered-pair energy contribution.
!> \param left_operator New left occupation derivative, including the prefactor.
!> \param right_operator New right occupation derivative, including the prefactor.
!> \param transition_adjoint Optional new complete transition energy adjoint.
!> \param interaction_adjoint Optional new effective-kernel energy adjoint.
!> \param comm Auxiliary point communicator. Replicated energy derivatives are not reduced again.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_exchange_pair(left_frame, right_frame, points, weights, potential, &
                                           augmentation, interaction, left_density, right_density, prefactor, &
                                           energy, left_operator, right_operator, transition_adjoint, &
                                           interaction_adjoint, comm)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: left_frame, right_frame
      REAL(dp), INTENT(IN)                               :: points(:, :), weights(:)
      COMPLEX(dp), INTENT(IN)                            :: potential(:, :), augmentation(:, :, :), &
                                                            interaction(:, :), left_density(:, :), &
                                                            right_density(:, :)
      REAL(dp), INTENT(IN)                               :: prefactor
      REAL(dp), INTENT(OUT)                              :: energy
      COMPLEX(dp), INTENT(OUT)                           :: left_operator(:, :), right_operator(:, :)
      COMPLEX(dp), INTENT(OUT), OPTIONAL                 :: transition_adjoint(:, :, :), &
                                                            interaction_adjoint(:, :)
      TYPE(mp_comm_type), INTENT(IN), OPTIONAL           :: comm

      COMPLEX(dp), ALLOCATABLE                           :: transitions(:, :, :)

      ALLOCATE (transitions(SIZE(left_frame%coefficients, 2), SIZE(right_frame%coefficients, 2), SIZE(potential, 2)))
      CALL qs_upf_inverse_transition(left_frame, right_frame, points, weights, potential, &
                                     augmentation, transitions, comm)
      CALL upf_inverse_exchange_contract(transitions, interaction, left_density, right_density, prefactor, &
                                         energy, left_operator, right_operator, transition_adjoint, interaction_adjoint)
   END SUBROUTINE qs_upf_inverse_exchange_pair

! **************************************************************************************************
!> \brief Complete ordered-pair energy response from original-frame and native RI responses.
!>        Includes both occupation matrices, the effective fitting/exchange interaction, and
!>        the full transition response. In particular no auxiliary-grid solve replaces dW.
!> \param left_frame ...
!> \param right_frame ...
!> \param points ...
!> \param weights ...
!> \param potential ...
!> \param augmentation ...
!> \param interaction ...
!> \param left_density ...
!> \param right_density ...
!> \param prefactor Geometry-independent ordered spin/k prefactor.
!> \param point_response ...
!> \param weights_response ...
!> \param potential_response ...
!> \param augmentation_response ...
!> \param interaction_response Full RI owner's effective-kernel response.
!> \param left_density_response Full left occupation response in the native AO coordinates.
!> \param right_density_response Full right occupation response.
!> \param atom_response ...
!> \param cell_response ...
!> \param left_k_response ...
!> \param right_k_response ...
!> \param left_coefficients_response ...
!> \param right_coefficients_response ...
!> \param energy New unperturbed ordered-pair energy.
!> \param response New joint directional energy response, counted once by its owner.
!> \param comm ...
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_exchange_pair_response(left_frame, right_frame, points, weights, potential, &
                                                    augmentation, interaction, left_density, right_density, prefactor, &
                                                    point_response, weights_response, potential_response, &
                                                    augmentation_response, interaction_response, left_density_response, &
                                                    right_density_response, atom_response, cell_response, &
                                                    left_k_response, right_k_response, left_coefficients_response, &
                                                    right_coefficients_response, energy, response, comm)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: left_frame, right_frame
      REAL(dp), INTENT(IN)                               :: points(:, :), weights(:)
      COMPLEX(dp), INTENT(IN)                            :: potential(:, :), augmentation(:, :, :), &
                                                            interaction(:, :), left_density(:, :), &
                                                            right_density(:, :)
      REAL(dp), INTENT(IN)                               :: prefactor, point_response(:, :), &
                                                            weights_response(:)
      COMPLEX(dp), INTENT(IN) :: potential_response(:, :), augmentation_response(:, :, :), &
         interaction_response(:, :), left_density_response(:, :), right_density_response(:, :)
      REAL(dp), INTENT(IN)                               :: atom_response(:, :), &
                                                            cell_response(3, 3), &
                                                            left_k_response(3), right_k_response(3)
      COMPLEX(dp), INTENT(IN) :: left_coefficients_response(:, :), &
         right_coefficients_response(:, :)
      REAL(dp), INTENT(OUT)                              :: energy, response
      TYPE(mp_comm_type), INTENT(IN), OPTIONAL           :: comm

      COMPLEX(dp), ALLOCATABLE                           :: di(:, :, :), gi(:, :, :), gr(:, :), &
                                                            hl(:, :), hr(:, :)
      INTEGER                                            :: naux, nl, nr

      nl = SIZE(left_frame%coefficients, 2)
      nr = SIZE(right_frame%coefficients, 2)
      naux = SIZE(potential, 2)
      CPASSERT(ALL(SHAPE(left_density_response) == [nl, nl]))
      CPASSERT(ALL(SHAPE(right_density_response) == [nr, nr]))
      CPASSERT(ALL(SHAPE(interaction_response) == [naux, naux]))
      ALLOCATE (di(nl, nr, naux), gi(nl, nr, naux), gr(naux, naux), hl(nl, nl), hr(nr, nr))
      CALL qs_upf_inverse_exchange_pair(left_frame, right_frame, points, weights, potential, augmentation, &
                                        interaction, left_density, right_density, prefactor, energy, hl, hr, gi, gr, comm)
      CALL qs_upf_inverse_transition_response(left_frame, right_frame, points, weights, potential, augmentation, &
                                              point_response, weights_response, potential_response, augmentation_response, &
                                              atom_response, cell_response, left_k_response, right_k_response, &
                                              left_coefficients_response, right_coefficients_response, di, comm)
      response = REAL(SUM(CONJG(gi)*di) + SUM(CONJG(gr)*interaction_response) + &
                      SUM(CONJG(hl)*left_density_response) + SUM(CONJG(hr)*right_density_response), dp)
   END SUBROUTINE qs_upf_inverse_exchange_pair_response

! **************************************************************************************************
!> \brief Center contribution d_a=-Delta_a*W_a to the inverse-minus-Gaussian field.
!>        W_a uses the prepared local inverse or the full coupled reference solve.
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

      CPASSERT(frame%atom_local .OR. frame%inverse%ready)
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
         CPASSERT(PRESENT(difference_response))
         CPASSERT(PRESENT(point_response) .AND. PRESENT(atom_response))
         CPASSERT(PRESENT(cell_response) .AND. PRESENT(reduced_k_response))
         CPASSERT(PRESENT(coefficients_response))
         CPASSERT(ALL(SHAPE(difference_response) == [np, nao]))
         IF (PRESENT(gradient_response)) THEN
            CPASSERT(ALL(SHAPE(gradient_response) == [np, nao, 3]))
         END IF
         CPASSERT(ALL(SHAPE(point_response) == SHAPE(points)))
         CPASSERT(ALL(SHAPE(atom_response) == SHAPE(frame%positions)))
         CPASSERT(ALL(SHAPE(coefficients_response) == SHAPE(frame%coefficients)))
         ALLOCATE (ddelta(np, nproj))
         IF (PRESENT(gradient_response)) THEN
            ALLOCATE (dgd(np, nproj, 3))
            CALL upf_radial_fields_periodic(frame%delta_cache(ikind), points, frame%positions(:, atom), &
                                            frame%cell, images, frame%reduced_k, delta, gd, ddelta, dgd, &
                                            point_response, atom_response(:, atom), cell_response, reduced_k_response)
         ELSE
            CALL upf_radial_fields_periodic(frame%delta_cache(ikind), points, frame%positions(:, atom), &
                                            frame%cell, images, frame%reduced_k, delta, gd, response=ddelta, &
                                            points_response=point_response, center_response=atom_response(:, atom), &
                                            cell_response=cell_response, reduced_k_response=reduced_k_response)
         END IF
         difference_response = -MATMUL(ddelta, frame%coefficients(first:last, :)) - &
                               MATMUL(delta, coefficients_response(first:last, :))
      ELSE
         CALL upf_radial_fields_periodic(frame%delta_cache(ikind), points, frame%positions(:, atom), &
                                         frame%cell, images, frame%reduced_k, delta, gd)
      END IF
      difference = -MATMUL(delta, frame%coefficients(first:last, :))
      DO axis = 1, 3
         gradient(:, :, axis) = -MATMUL(gd(:, :, axis), frame%coefficients(first:last, :))
         IF (PRESENT(gradient_response)) THEN
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
!> \param soft_response Optional response of the inverse AO values.
!> \param gradient_response Optional response of their physical spatial gradients.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_response(frame, point_response, atom_response, cell_response, reduced_k_response, &
                                      weights_response, coefficients_response, soft_response, gradient_response)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      REAL(dp), INTENT(IN) :: point_response(:, :), atom_response(:, :), cell_response(3, 3), &
         reduced_k_response(3), weights_response(:)
      COMPLEX(dp), INTENT(OUT)                           :: coefficients_response(:, :)
      COMPLEX(dp), INTENT(OUT), OPTIONAL                 :: soft_response(:, :), &
                                                            gradient_response(:, :, :)

      COMPLEX(dp), ALLOCATABLE                           :: db(:, :), dbeta(:, :), dchi(:, :), &
                                                            ddelta(:, :), dgc(:, :, :), &
                                                            dgd(:, :, :), dk(:, :)
      INTEGER                                            :: axis, nao, npoint, nproj

      CPASSERT(frame%atom_local .OR. frame%inverse%ready)
      CPASSERT(ALL(SHAPE(point_response) == SHAPE(frame%points)))
      CPASSERT(ALL(SHAPE(atom_response) == SHAPE(frame%positions)))
      CPASSERT(SIZE(weights_response) == SIZE(frame%weights))
      npoint = SIZE(frame%weights)
      nao = SIZE(frame%gaussian, 2)
      nproj = SIZE(frame%coefficients, 1)
      CPASSERT(ALL(SHAPE(coefficients_response) == [nproj, nao]))
      IF (PRESENT(soft_response)) THEN
         CPASSERT(ALL(SHAPE(soft_response) == [npoint, nao]))
      END IF
      IF (PRESENT(gradient_response)) THEN
         CPASSERT(ALL(SHAPE(gradient_response) == [npoint, nao, 3]))
      END IF
      ALLOCATE (dbeta(npoint, nproj), ddelta(npoint, nproj), dchi(npoint, nao), dk(nproj, nproj), db(nproj, nao))
      IF (PRESENT(gradient_response)) THEN
         ALLOCATE (dgd(npoint, nproj, 3), dgc(npoint, nao, 3))
      END IF
      IF (PRESENT(gradient_response)) THEN
         CALL defining_responses(frame, point_response, atom_response, cell_response, reduced_k_response, &
                                 dbeta, ddelta, dchi, dgc, dgd)
      ELSE
         CALL defining_responses(frame, point_response, atom_response, cell_response, reduced_k_response, dbeta, ddelta, dchi)
      END IF
      IF (frame%atom_local) THEN
         CALL atomic_coefficients(frame, db, coefficients_response, atom_response, cell_response, reduced_k_response)
      ELSE
         CALL upf_inverse_pair_response(frame%beta, frame%delta, frame%weights, dbeta, ddelta, &
                                        weights_response, dk, frame%point_comm)
         CALL upf_inverse_pair_response(frame%beta, frame%gaussian, frame%weights, dbeta, dchi, &
                                        weights_response, db, frame%point_comm)
         CALL upf_inverse_project_response(frame%inverse, frame%coefficients, db, dk, coefficients_response)
      END IF
      IF (PRESENT(soft_response)) THEN
         CALL upf_inverse_basis_response(dchi, frame%delta, ddelta, frame%coefficients, &
                                         coefficients_response, soft_response)
      END IF
      IF (PRESENT(gradient_response)) THEN
         DO axis = 1, 3
            CALL upf_inverse_basis_response(dgc(:, :, axis), frame%delta_gradient(:, :, axis), dgd(:, :, axis), &
                                            frame%coefficients, coefficients_response, gradient_response(:, :, axis))
         END DO
      END IF
   END SUBROUTINE qs_upf_inverse_response

! **************************************************************************************************
!> \brief Responses of the original defining fields, shared by forward and adjoint projection.
!> \param frame ...
!> \param point_response ...
!> \param atom_response ...
!> \param cell_response ...
!> \param reduced_k_response ...
!> \param dbeta ...
!> \param ddelta ...
!> \param dchi ...
!> \param gaussian_gradient_response Optional physical Gaussian-gradient response.
!> \param delta_gradient_response Optional physical difference-gradient response.
! **************************************************************************************************
   SUBROUTINE defining_responses(frame, point_response, atom_response, cell_response, reduced_k_response, &
                                 dbeta, ddelta, dchi, gaussian_gradient_response, delta_gradient_response)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      REAL(dp), INTENT(IN)                               :: point_response(:, :), &
                                                            atom_response(:, :), &
                                                            cell_response(3, 3), &
                                                            reduced_k_response(3)
      COMPLEX(dp), INTENT(OUT)                           :: dbeta(:, :), ddelta(:, :), dchi(:, :)
      COMPLEX(dp), INTENT(OUT), OPTIONAL :: gaussian_gradient_response(:, :, :), &
         delta_gradient_response(:, :, :)

      COMPLEX(dp), ALLOCATABLE                           :: beta(:, :), chi(:, :), delta(:, :), &
                                                            gc(:, :, :), gd(:, :, :)
      INTEGER                                            :: atom, center, first, ikind, last, nao, &
                                                            npoint, nproj
      LOGICAL                                            :: shared_response

      npoint = SIZE(frame%weights)
      nao = SIZE(frame%coefficients, 2)
      nproj = SIZE(frame%coefficients, 1)
      CPASSERT(PRESENT(gaussian_gradient_response) .EQV. PRESENT(delta_gradient_response))
      CPASSERT(ALL(SHAPE(dbeta) == [npoint, nproj]) .AND. ALL(SHAPE(ddelta) == [npoint, nproj]))
      CPASSERT(ALL(SHAPE(dchi) == [npoint, nao]))
      ALLOCATE (beta(npoint, nproj), delta(npoint, nproj), chi(npoint, nao))
      IF (PRESENT(gaussian_gradient_response)) THEN
         ALLOCATE (gc(npoint, nao, 3), gd(npoint, nproj, 3))
      END IF
      shared_response = ANY(point_response /= 0.0_dp) .OR. ANY(cell_response /= 0.0_dp) .OR. &
                        ANY(reduced_k_response /= 0.0_dp)
      DO atom = 1, SIZE(frame%bases)
         first = frame%ao_first(atom)
         last = frame%ao_last(atom)
         IF (.NOT. shared_response .AND. ALL(atom_response(:, atom) == 0.0_dp)) THEN
            dchi(:, first:last) = CMPLX(0.0_dp, 0.0_dp, dp)
            IF (PRESENT(gaussian_gradient_response)) gaussian_gradient_response(:, first:last, :) = CMPLX(0.0_dp, 0.0_dp, dp)
            CYCLE
         END IF
         IF (PRESENT(gaussian_gradient_response)) THEN
            CALL gaussian_basis_fields_periodic(frame%bases(atom)%gto_basis_set, frame%points, frame%positions(:, atom), &
                                              frame%cell, frame%images, frame%reduced_k, chi(:, first:last), gc(:, first:last, :), &
                        dchi(:, first:last), gaussian_gradient_response(:, first:last, :), point_response, atom_response(:, atom), &
                                                cell_response, reduced_k_response)
         ELSE
            CALL gaussian_basis_fields_periodic(frame%bases(atom)%gto_basis_set, frame%points, frame%positions(:, atom), &
                                                frame%cell, frame%images, frame%reduced_k, chi(:, first:last), &
                                                response=dchi(:, first:last), point_response=point_response, &
                                                center_response=atom_response(:, atom), cell_response=cell_response, &
                                                reduced_k_response=reduced_k_response)
         END IF
      END DO
      DO center = 1, SIZE(frame%center_atoms)
         atom = frame%center_atoms(center)
         ikind = frame%center_kinds(center)
         first = frame%projector_first(center)
         last = frame%projector_last(center)
         IF (.NOT. shared_response .AND. ALL(atom_response(:, atom) == 0.0_dp)) THEN
            dbeta(:, first:last) = CMPLX(0.0_dp, 0.0_dp, dp)
            ddelta(:, first:last) = CMPLX(0.0_dp, 0.0_dp, dp)
            IF (PRESENT(gaussian_gradient_response)) delta_gradient_response(:, first:last, :) = CMPLX(0.0_dp, 0.0_dp, dp)
            CYCLE
         END IF
         CALL upf_radial_fields_periodic(frame%beta_cache(ikind), frame%points, frame%positions(:, atom), &
                                    frame%cell, frame%images, frame%reduced_k, beta(:, first:last), response=dbeta(:, first:last), &
                                         points_response=point_response, center_response=atom_response(:, atom), &
                                         cell_response=cell_response, reduced_k_response=reduced_k_response)
         IF (PRESENT(gaussian_gradient_response)) THEN
            CALL upf_radial_fields_periodic(frame%delta_cache(ikind), frame%points, frame%positions(:, atom), &
                                            frame%cell, frame%images, frame%reduced_k, delta(:, first:last), gd(:, first:last, :), &
                         ddelta(:, first:last), delta_gradient_response(:, first:last, :), point_response, atom_response(:, atom), &
                                            cell_response, reduced_k_response)
         ELSE
            CALL upf_radial_fields_periodic(frame%delta_cache(ikind), frame%points, frame%positions(:, atom), &
                                            frame%cell, frame%images, frame%reduced_k, delta(:, first:last), &
                                            response=ddelta(:, first:last), points_response=point_response, &
                                            center_response=atom_response(:, atom), cell_response=cell_response, &
                                            reduced_k_response=reduced_k_response)
         END IF
      END DO
   END SUBROUTINE defining_responses

! **************************************************************************************************
!> \brief Pull a replicated real-energy coefficient covector back through the defining projection.
!>        The RI owner sums all blocks and both orbital contributions before this call. The full
!>        inverse uses K**H*GB=GW and GK=-GB*W**H. Atom-local hard projectors are fixed radial
!>        operators, so their Gaussian-overlap covector is GW and their radial-K response is zero.
!> \param frame Original prepared frame in the native AO/projector order.
!> \param coefficients_covector Replicated GW with pairing Re sum(conjg(GW)*dW).
!> \param beta_gaussian_covector New GB for the original Gaussian overlaps.
!> \param beta_delta_covector New GK for the full coupling, zero for the atom-local formulation.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_projection_adjoint(frame, coefficients_covector, &
                                                beta_gaussian_covector, beta_delta_covector)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      COMPLEX(dp), INTENT(IN)                            :: coefficients_covector(:, :)
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: beta_gaussian_covector(:, :), &
                                                            beta_delta_covector(:, :)

      INTEGER                                            :: nproj

      CPASSERT(frame%atom_local .OR. frame%inverse%ready)
      CPASSERT(ALL(SHAPE(coefficients_covector) == SHAPE(frame%coefficients)))
      nproj = SIZE(frame%coefficients, 1)
      ALLOCATE (beta_gaussian_covector, MOLD=frame%coefficients)
      ALLOCATE (beta_delta_covector(nproj, nproj))
      IF (frame%atom_local) THEN
         beta_gaussian_covector(:, :) = coefficients_covector
         beta_delta_covector = CMPLX(0.0_dp, 0.0_dp, dp)
      ELSE
         CALL upf_inverse_project_adjoint(frame%inverse, frame%coefficients, coefficients_covector, &
                                          beta_gaussian_covector, beta_delta_covector)
      END IF
   END SUBROUTINE qs_upf_inverse_projection_adjoint

! **************************************************************************************************
!> \brief Contract defining-projection covectors with a physical geometry or Bloch response.
!>        GB and GK are prepared once, then reused for different directions without another LU
!>        solve. The general inverse retains every original field, image and measure term and
!>        reduces the local real scalar once. The atom-local analytic overlaps are replicated
!>        and need no additional reduction. Native occupation or mesh factors are already in G.
!> \param frame Original prepared frame, including its defining points and integration measure.
!> \param beta_gaussian_covector Replicated GB from qs_upf_inverse_projection_adjoint.
!> \param beta_delta_covector Replicated GK from the same frame.
!> \param point_response Motion of the original defining points.
!> \param atom_response Physical motion of all native atoms and Gaussian centers.
!> \param cell_response Cartesian lattice response.
!> \param reduced_k_response Response of the reduced Bloch vector.
!> \param weights_response Response of the original local measure.
!> \param response Replicated real-energy response through W, excluding explicit RI field terms.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_projection_contract(frame, beta_gaussian_covector, beta_delta_covector, &
                                                 point_response, atom_response, cell_response, reduced_k_response, &
                                                 weights_response, response)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      COMPLEX(dp), INTENT(IN)                            :: beta_gaussian_covector(:, :), &
                                                            beta_delta_covector(:, :)
      REAL(dp), INTENT(IN) :: point_response(:, :), atom_response(:, :), cell_response(3, 3), &
         reduced_k_response(3), weights_response(:)
      REAL(dp), INTENT(OUT)                              :: response

      COMPLEX(dp), ALLOCATABLE                           :: db(:, :), dbeta(:, :), dchi(:, :), &
                                                            ddelta(:, :), dk(:, :), overlaps(:, :)
      INTEGER                                            :: nao, npoint, nproj

      CPASSERT(frame%atom_local .OR. frame%inverse%ready)
      nproj = SIZE(frame%coefficients, 1)
      nao = SIZE(frame%coefficients, 2)
      npoint = SIZE(frame%weights)
      CPASSERT(ALL(SHAPE(beta_gaussian_covector) == [nproj, nao]))
      CPASSERT(ALL(SHAPE(beta_delta_covector) == [nproj, nproj]))
      CPASSERT(ALL(SHAPE(point_response) == SHAPE(frame%points)))
      CPASSERT(ALL(SHAPE(atom_response) == SHAPE(frame%positions)))
      CPASSERT(SIZE(weights_response) == npoint)
      ALLOCATE (db(nproj, nao))
      IF (frame%atom_local) THEN
         CPASSERT(ALL(beta_delta_covector == CMPLX(0.0_dp, 0.0_dp, dp)))
         ALLOCATE (overlaps(nproj, nao))
         CALL atomic_coefficients(frame, overlaps, db, atom_response, cell_response, reduced_k_response)
         response = REAL(SUM(CONJG(beta_gaussian_covector)*db), dp)
      ELSE
         ALLOCATE (dbeta(npoint, nproj), ddelta(npoint, nproj), dchi(npoint, nao), dk(nproj, nproj))
         CALL defining_responses(frame, point_response, atom_response, cell_response, reduced_k_response, dbeta, ddelta, dchi)
         CALL upf_inverse_pair_response(frame%beta, frame%gaussian, frame%weights, dbeta, dchi, weights_response, db)
         CALL upf_inverse_pair_response(frame%beta, frame%delta, frame%weights, dbeta, ddelta, weights_response, dk)
         response = REAL(SUM(CONJG(beta_gaussian_covector)*db) + SUM(CONJG(beta_delta_covector)*dk), dp)
         CALL frame%point_comm%sum(response)
      END IF
   END SUBROUTINE qs_upf_inverse_projection_contract

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
