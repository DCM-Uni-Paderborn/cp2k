!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief PAW one-center functional in the native inverse Gaussian representation.
!>        Fixed kinetic and ionic terms remain in the original PP_DIJ.
! **************************************************************************************************
MODULE qs_upf_inverse_paw
   USE hfx_types,                       ONLY: hfx_type
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE input_constants,                 ONLY: do_potential_coulomb,&
                                              do_potential_long,&
                                              do_potential_mix_cl,&
                                              do_potential_mix_cl_trunc,&
                                              do_potential_short,&
                                              do_potential_truncated
   USE input_section_types,             ONLY: section_vals_type
   USE kinds,                           ONLY: dp
   USE mathlib,                         ONLY: invmat
   USE message_passing,                 ONLY: mp_comm_type
   USE paw_onecenter,                   ONLY: paw_onecenter_core_exchange,&
                                              paw_onecenter_core_exchange_init,&
                                              paw_onecenter_core_exchange_prepare,&
                                              paw_onecenter_hartree,&
                                              paw_onecenter_init
   USE paw_types,                       ONLY: paw_kind_type,&
                                              paw_potential_type
   USE paw_xc,                          ONLY: paw_xc_evaluate,&
                                              paw_xc_init
   USE qs_upf_inverse,                  ONLY: qs_upf_inverse_response
   USE qs_upf_inverse_types,            ONLY: qs_upf_inverse_augmentation_type,&
                                              qs_upf_inverse_grid_type,&
                                              qs_upf_inverse_type
   USE upf_inverse,                     ONLY: upf_inverse_density_evaluate,&
                                              upf_inverse_kinetic_type
   USE upf_projector_integrals,         ONLY: upf_projector_type
   USE upf_types,                       ONLY: atom_upfpot_type
   USE uspp_augmentation,               ONLY: uspp_augmentation_realspace,&
                                              uspp_augmentation_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_upf_inverse_paw'

   PUBLIC :: qs_upf_inverse_paw_init, qs_upf_inverse_paw_evaluate, qs_upf_inverse_paw_response
   PUBLIC :: qs_upf_inverse_paw_core_prepare
   PUBLIC :: qs_upf_inverse_paw_collocate, qs_upf_inverse_paw_operator

CONTAINS

! **************************************************************************************************
!> \brief Closed-core exchange with the same physical interactions as the native valence exchange.
!>        Sum each HF contribution once, independently of its thread replicas. The valence-valence
!>        exchange is already evaluated through the reconstructed RI transitions.
!> \param pot Original scalar PAW data, including the complete signed core orbitals.
!> \param projectors Original magnetic-projector ordering.
!> \param kind Original-data caches, reused without changing their interactions.
!> \param exchange Native HF contributions and thread replicas, with parameters in atomic units.
!> \param nspins Number of physical valence spins, used for the native RI activity threshold.
!> \param energy Frozen core-core exchange constant, summed over the selected HF contributions.
!> \param potential Spin-independent core-valence projector operator, without an extra factor 1/2.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_paw_core_prepare(pot, projectors, kind, exchange, nspins, energy, potential)
      TYPE(atom_upfpot_type), INTENT(IN)                 :: pot
      TYPE(upf_projector_type), INTENT(IN)               :: projectors
      TYPE(paw_kind_type), INTENT(INOUT)                 :: kind
      TYPE(hfx_type), INTENT(IN)                         :: exchange(:, :)
      INTEGER, INTENT(IN)                                :: nspins
      REAL(dp), INTENT(OUT)                              :: energy, potential(:, :)

      INTEGER                                            :: cache, contribution, interaction(2), term
      REAL(dp)                                           :: coefficient(2), fraction, value
      REAL(dp), ALLOCATABLE                              :: matrix(:, :, :), zero_density(:, :, :)

      CPASSERT(pot%is_paw .AND. kind%active)
      CPASSERT(nspins == 1 .OR. nspins == 2)
      CPASSERT(SIZE(exchange, 2) > 0)
      CPASSERT(ALL(SHAPE(potential) == [projectors%nproj, projectors%nproj]))
      energy = 0.0_dp
      potential = 0.0_dp
      ALLOCATE (zero_density(projectors%nproj, projectors%nproj, 1), SOURCE=0.0_dp)
      ALLOCATE (matrix(projectors%nproj, projectors%nproj, 1))
      DO contribution = 1, SIZE(exchange, 1)
         ASSOCIATE (hf => exchange(contribution, 1))
            fraction = hf%general_parameter%fraction
            CPASSERT(ieee_is_finite(fraction))
            ! Match the spin prefactor used by hfx_ri_update_ks and its Fock contraction.
            IF (.NOT. MERGE(0.5_dp*fraction, fraction, nspins == 1) > EPSILON(0.0_dp)) CYCLE
            coefficient = 0.0_dp
            interaction = 0
            SELECT CASE (hf%potential_parameter%potential_type)
            CASE (do_potential_coulomb)
               coefficient(1) = fraction
            CASE (do_potential_short)
               coefficient(1) = fraction
               interaction(1) = 1
            CASE (do_potential_long)
               coefficient = [fraction, -fraction]
               interaction = [0, 1]
            CASE (do_potential_mix_cl)
               coefficient = fraction*[hf%potential_parameter%scale_coulomb + &
                                       hf%potential_parameter%scale_longrange, &
                                       -hf%potential_parameter%scale_longrange]
               interaction = [0, 1]
            CASE (do_potential_truncated)
               coefficient(1) = fraction
               interaction(1) = 2
            CASE (do_potential_mix_cl_trunc)
               coefficient = fraction*[hf%potential_parameter%scale_coulomb + &
                                       hf%potential_parameter%scale_longrange, &
                                       -hf%potential_parameter%scale_longrange]
               interaction = [2, 1]
            CASE DEFAULT
               CPABORT("Inverse PAW core exchange for this interaction is not implemented yet")
            END SELECT
            CPASSERT(ALL(ieee_is_finite(coefficient)))
            DO term = 1, 2
               IF (coefficient(term) == 0.0_dp) CYCLE
               SELECT CASE (interaction(term))
               CASE (0)
                  IF (.NOT. kind%onecenter%core_exchange_initialized) THEN
                     CALL paw_onecenter_core_exchange_init(pot, projectors, kind%onecenter)
                  END IF
                  IF (.NOT. kind%onecenter%core_exchange_available) THEN
                     CPABORT(TRIM(kind%onecenter%core_exchange_reason))
                  END IF
                  cache = 0
               CASE (1)
                  CALL paw_onecenter_core_exchange_prepare(pot, projectors, kind%onecenter, cache_index=cache, &
                                                           omega=hf%potential_parameter%omega)
               CASE (2)
                  CALL paw_onecenter_core_exchange_prepare(pot, projectors, kind%onecenter, cache_index=cache, &
                                                           cutoff_radius=hf%potential_parameter%cutoff_radius)
               END SELECT
               IF (cache > 0) THEN
                  IF (.NOT. kind%onecenter%interaction_core_exchange(cache)%available) THEN
                     CPABORT(TRIM(kind%onecenter%interaction_core_exchange(cache)%reason))
                  END IF
               END IF
               CALL paw_onecenter_core_exchange(kind%onecenter, zero_density, coefficient(term), cache, value, matrix)
               energy = energy + value
               potential = potential + matrix(:, :, 1)
            END DO
         END ASSOCIATE
      END DO
   END SUBROUTINE qs_upf_inverse_paw_core_prepare

! **************************************************************************************************
!> \brief PAW valence reconstruction relative to the compensated soft density, including its gradient.
!>        Adds AE-minus-PS partial-wave products minus the original compensation Q. No core is added.
!>        Same-center kernels are periodic without Bloch phases; the inverse coefficients carry k.
!> \param frame Original inverse frame and finite image list.
!> \param reconstruction Kind-indexed original AE/PS caches. Inactive kinds have zero channels.
!> \param augmentation Original compensation kernels in the same magnetic-projector order.
!> \param points Cartesian evaluation points without quadrature or k weights.
!> \param matrix Hermitian AO occupation matrix for one spin and k point.
!> \param density Reconstructed valence correction, overwritten.
!> \param gradient Cartesian gradient of the same correction, overwritten.
!> \param density_response Optional fixed-matrix geometry response, requiring all response inputs.
!> \param gradient_response Optional response of the Cartesian gradient.
!> \param point_response Cartesian motion of the points.
!> \param atom_response Cartesian motion of the original atoms.
!> \param cell_response Lattice response, including translated image centers.
!> \param coefficients_response Complete coupled inverse response dW.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_paw_collocate(frame, reconstruction, augmentation, points, matrix, density, gradient, &
                                           density_response, gradient_response, point_response, atom_response, &
                                           cell_response, coefficients_response)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      TYPE(upf_inverse_kinetic_type), INTENT(IN)         :: reconstruction(:)
      TYPE(qs_upf_inverse_augmentation_type), INTENT(IN) :: augmentation
      REAL(dp), INTENT(IN)                               :: points(:, :)
      COMPLEX(dp), INTENT(IN)                            :: matrix(:, :)
      REAL(dp), INTENT(OUT)                              :: density(:), gradient(:, :)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: density_response(:), &
                                                            gradient_response(:, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: point_response(:, :), &
                                                            atom_response(:, :), &
                                                            cell_response(3, 3)
      COMPLEX(dp), INTENT(IN), OPTIONAL                  :: coefficients_response(:, :)

      COMPLEX(dp), ALLOCATABLE                           :: projected(:, :), projected_response(:, :)
      INTEGER                                            :: atom, axis, center, first, ikind, image, &
                                                            last, n, nao, np, point, response_axis
      LOGICAL                                            :: derive
      REAL(dp)                                           :: motion(3), position(3), &
                                                            position_response(3), radius
      REAL(dp), ALLOCATABLE                              :: kernel(:, :), kernel_gradient(:, :, :), &
                                                            kernel_hessian(:, :, :, :)

      nao = SIZE(frame%coefficients, 2)
      np = SIZE(points, 2)
      CPASSERT(frame%atom_local .OR. frame%inverse%ready)
      CPASSERT(SIZE(points, 1) == 3 .AND. ALL(SHAPE(matrix) == [nao, nao]))
      CPASSERT(SIZE(density) == np .AND. ALL(SHAPE(gradient) == [np, 3]))
      CPASSERT(ALLOCATED(augmentation%kinds))
      CPASSERT(SIZE(reconstruction) == SIZE(augmentation%kinds))
      derive = PRESENT(density_response)
      CPASSERT(derive .EQV. PRESENT(gradient_response))
      CPASSERT(derive .EQV. PRESENT(point_response))
      CPASSERT(derive .EQV. PRESENT(atom_response))
      CPASSERT(derive .EQV. PRESENT(cell_response))
      CPASSERT(derive .EQV. PRESENT(coefficients_response))
      density = 0.0_dp
      gradient = 0.0_dp
      IF (derive) THEN
         CPASSERT(SIZE(density_response) == np)
         CPASSERT(ALL(SHAPE(gradient_response) == [np, 3]))
         CPASSERT(ALL(SHAPE(point_response) == SHAPE(points)))
         CPASSERT(ALL(SHAPE(atom_response) == SHAPE(frame%positions)))
         CPASSERT(ALL(SHAPE(coefficients_response) == SHAPE(frame%coefficients)))
         density_response = 0.0_dp
         gradient_response = 0.0_dp
      END IF
      DO center = 1, SIZE(frame%center_atoms)
         ikind = frame%center_kinds(center)
         CPASSERT(ikind > 0 .AND. ikind <= SIZE(reconstruction))
         IF (reconstruction(ikind)%delta%nproj == 0) CYCLE
         atom = frame%center_atoms(center)
         first = frame%projector_first(center)
         last = frame%projector_last(center)
         n = last - first + 1
         CPASSERT(n == reconstruction(ikind)%delta%nproj)
         CPASSERT(n == augmentation%kinds(ikind)%nproj)
         radius = MAX(reconstruction(ikind)%delta%radius, augmentation%kinds(ikind)%realspace_radius)
         ALLOCATE (projected(n, n), kernel(n, n), kernel_gradient(n, n, 3))
         projected(:, :) = MATMUL(frame%coefficients(first:last, :), &
                                  MATMUL(matrix, CONJG(TRANSPOSE(frame%coefficients(first:last, :)))))
         IF (derive) THEN
            ALLOCATE (projected_response(n, n), kernel_hessian(n, n, 3, 3))
            projected_response(:, :) = MATMUL(coefficients_response(first:last, :), &
                                              MATMUL(matrix, CONJG(TRANSPOSE(frame%coefficients(first:last, :)))))
            projected_response(:, :) = projected_response + CONJG(TRANSPOSE(projected_response))
         END IF
         DO image = 1, SIZE(frame%images, 2)
            position = frame%positions(:, atom) + MATMUL(frame%cell, REAL(frame%images(:, image), dp))
            IF (derive) position_response = atom_response(:, atom) + &
                                            MATMUL(cell_response, REAL(frame%images(:, image), dp))
            DO point = 1, np
               IF (SUM((points(:, point) - position)**2) >= radius**2) CYCLE
               IF (derive) THEN
                  CALL reconstruction_fields(reconstruction(ikind), augmentation%kinds(ikind), &
                                             points(:, point) - position, kernel, kernel_gradient, kernel_hessian)
                  motion = point_response(:, point) - position_response
                  density_response(point) = density_response(point) + REAL(SUM(projected_response*kernel), dp)
                  DO axis = 1, 3
                     density_response(point) = density_response(point) + REAL(SUM(projected*kernel_gradient(:, :, axis)), dp)* &
                                               motion(axis)
                     gradient_response(point, axis) = gradient_response(point, axis) + &
                                                      REAL(SUM(projected_response*kernel_gradient(:, :, axis)), dp)
                     DO response_axis = 1, 3
                        gradient_response(point, axis) = gradient_response(point, axis) + &
                                                         REAL(SUM(projected*kernel_hessian(:, :, axis, response_axis)), dp)* &
                                                         motion(response_axis)
                     END DO
                  END DO
               ELSE
                  CALL reconstruction_fields(reconstruction(ikind), augmentation%kinds(ikind), &
                                             points(:, point) - position, kernel, kernel_gradient)
               END IF
               density(point) = density(point) + REAL(SUM(projected*kernel), dp)
               DO axis = 1, 3
                  gradient(point, axis) = gradient(point, axis) + REAL(SUM(projected*kernel_gradient(:, :, axis)), dp)
               END DO
            END DO
         END DO
         DEALLOCATE (projected, kernel, kernel_gradient)
         IF (derive) DEALLOCATE (projected_response, kernel_hessian)
      END DO
   END SUBROUTINE qs_upf_inverse_paw_collocate

! **************************************************************************************************
!> \brief Exact AO transpose of the same PAW density and Cartesian-gradient reconstruction.
!> \param frame ...
!> \param reconstruction ...
!> \param augmentation ...
!> \param points ...
!> \param weights Integration measure, applied once. Weighted model derivatives use unit weights.
!> \param vdensity Potential conjugate to the reconstructed density.
!> \param vgradient Potential conjugate to its Cartesian gradient.
!> \param ao_operator Hermitian AO operator without spin or k weights.
!> \param comm Optional communicator distributing point rows. Replicated points omit it.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_paw_operator(frame, reconstruction, augmentation, points, weights, vdensity, vgradient, &
                                          ao_operator, comm)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      TYPE(upf_inverse_kinetic_type), INTENT(IN)         :: reconstruction(:)
      TYPE(qs_upf_inverse_augmentation_type), INTENT(IN) :: augmentation
      REAL(dp), INTENT(IN)                               :: points(:, :), weights(:), vdensity(:), &
                                                            vgradient(:, :)
      COMPLEX(dp), INTENT(OUT)                           :: ao_operator(:, :)
      TYPE(mp_comm_type), INTENT(IN), OPTIONAL           :: comm

      INTEGER                                            :: atom, axis, center, first, ikind, image, &
                                                            last, n, nao, np, point
      REAL(dp)                                           :: position(3), radius
      REAL(dp), ALLOCATABLE                              :: kernel(:, :), kernel_gradient(:, :, :), &
                                                            OPERATOR(:, :)

      nao = SIZE(frame%coefficients, 2)
      np = SIZE(points, 2)
      CPASSERT(frame%atom_local .OR. frame%inverse%ready)
      CPASSERT(SIZE(points, 1) == 3 .AND. SIZE(weights) == np .AND. SIZE(vdensity) == np)
      CPASSERT(ALL(SHAPE(vgradient) == [np, 3]) .AND. ALL(SHAPE(ao_operator) == [nao, nao]))
      CPASSERT(ALLOCATED(augmentation%kinds))
      CPASSERT(SIZE(reconstruction) == SIZE(augmentation%kinds))
      ao_operator = 0.0_dp
      DO center = 1, SIZE(frame%center_atoms)
         ikind = frame%center_kinds(center)
         CPASSERT(ikind > 0 .AND. ikind <= SIZE(reconstruction))
         IF (reconstruction(ikind)%delta%nproj == 0) CYCLE
         atom = frame%center_atoms(center)
         first = frame%projector_first(center)
         last = frame%projector_last(center)
         n = last - first + 1
         CPASSERT(n == reconstruction(ikind)%delta%nproj)
         CPASSERT(n == augmentation%kinds(ikind)%nproj)
         radius = MAX(reconstruction(ikind)%delta%radius, augmentation%kinds(ikind)%realspace_radius)
         ALLOCATE (kernel(n, n), kernel_gradient(n, n, 3), OPERATOR(n, n))
         OPERATOR(:, :) = 0.0_dp
         DO image = 1, SIZE(frame%images, 2)
            position = frame%positions(:, atom) + MATMUL(frame%cell, REAL(frame%images(:, image), dp))
            DO point = 1, np
               IF (SUM((points(:, point) - position)**2) >= radius**2) CYCLE
               CALL reconstruction_fields(reconstruction(ikind), augmentation%kinds(ikind), &
                                          points(:, point) - position, kernel, kernel_gradient)
               OPERATOR(:, :) = OPERATOR + weights(point)*vdensity(point)*kernel
               DO axis = 1, 3
                  OPERATOR(:, :) = OPERATOR + weights(point)*vgradient(point, axis)*kernel_gradient(:, :, axis)
               END DO
            END DO
         END DO
         ao_operator = ao_operator + MATMUL(CONJG(TRANSPOSE(frame%coefficients(first:last, :))), &
                                            MATMUL(OPERATOR, frame%coefficients(first:last, :)))
         DEALLOCATE (kernel, kernel_gradient, OPERATOR)
      END DO
      IF (PRESENT(comm)) CALL comm%sum(ao_operator)
   END SUBROUTINE qs_upf_inverse_paw_operator

! **************************************************************************************************
!> \brief Original AE-minus-PS products minus original compensation, with shared spatial derivatives.
!> \param reconstruction ...
!> \param augmentation ...
!> \param displacement ...
!> \param kernel ...
!> \param gradient ...
!> \param hessian ...
! **************************************************************************************************
   SUBROUTINE reconstruction_fields(reconstruction, augmentation, displacement, kernel, gradient, hessian)
      TYPE(upf_inverse_kinetic_type), INTENT(IN)         :: reconstruction
      TYPE(uspp_augmentation_type), INTENT(IN)           :: augmentation
      REAL(dp), INTENT(IN)                               :: displacement(3)
      REAL(dp), INTENT(OUT)                              :: kernel(:, :), gradient(:, :, :)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: hessian(:, :, :, :)

      INTEGER                                            :: n
      REAL(dp), ALLOCATABLE                              :: q(:, :), qgradient(:, :, :), &
                                                            qhessian(:, :, :, :)

      n = reconstruction%delta%nproj
      ALLOCATE (q(n, n), qgradient(n, n, 3))
      CALL upf_inverse_density_evaluate(reconstruction, displacement, kernel, gradient, hessian)
      IF (PRESENT(hessian)) THEN
         ALLOCATE (qhessian(n, n, 3, 3))
         CALL uspp_augmentation_realspace(augmentation, displacement, q, qgradient, qhessian)
         hessian = hessian - qhessian
      ELSE
         CALL uspp_augmentation_realspace(augmentation, displacement, q, qgradient)
      END IF
      kernel = kernel - q
      gradient = gradient - qgradient
   END SUBROUTINE reconstruction_fields

! **************************************************************************************************
!> \brief Prepare original PAW radial kernels without changing the supplied operators or waves.
!> \param pot Original scalar PAW dataset.
!> \param projectors Matching radial and magnetic projector mapping.
!> \param nang Requested atomic angular quadrature.
!> \param cache Owned per-kind PAW data.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_paw_init(pot, projectors, nang, cache)
      TYPE(atom_upfpot_type), INTENT(IN)                 :: pot
      TYPE(upf_projector_type), INTENT(IN)               :: projectors
      INTEGER, INTENT(IN)                                :: nang
      TYPE(paw_kind_type), INTENT(OUT)                   :: cache

      CPASSERT(pot%is_paw .AND. .NOT. pot%has_so)
      CPASSERT(ALLOCATED(pot%paw))
      CALL paw_onecenter_init(pot, projectors, cache%onecenter)
      CALL paw_xc_init(pot, nang, cache%xcgrid)
      cache%core_energy = pot%paw%core_energy
      cache%active = .TRUE.
   END SUBROUTINE qs_upf_inverse_paw_init

! **************************************************************************************************
!> \brief Evaluate PAW Hartree, XC and core energy and their Gaussian Hamiltonian derivative.
!> \param kinds Prepared per-kind PAW kernels. Inactive NC and USPP kinds contribute zero.
!> \param frames Native inverse frames for the complete global k mesh, identically ordered on all ranks.
!> \param density Hermitian Gaussian occupation matrices, AO by AO by global k point by spin.
!> \param kweights Normalized global k weights. They enter the projected density once.
!> \param functional Activated XC_FUNCTIONAL section.
!> \param density_cutoff XC density cutoff.
!> \param gradient_cutoff XC gradient cutoff.
!> \param energies One-center Hartree, XC and frozen-core constant, counted once per atom.
!> \param hamiltonian Unweighted Gaussian Hamiltonian correction for each global k point and spin.
!> \param potential Owned per-center projector potentials for the corresponding geometry response.
!> \param semilocal_xc Include atomic semilocal XC. A separate composite XC consumer may disable it.
!> \param core_energy Optional kind-indexed core-core exchange constants for the native HF interactions.
!> \param core_potential Optional kind-indexed core-valence operators, with one spin-independent matrix.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_paw_evaluate(kinds, frames, density, kweights, functional, density_cutoff, &
                                          gradient_cutoff, energies, hamiltonian, potential, semilocal_xc, &
                                          core_energy, core_potential)
      TYPE(paw_kind_type), INTENT(IN)                    :: kinds(:)
      TYPE(qs_upf_inverse_grid_type), INTENT(IN)         :: frames(:)
      COMPLEX(dp), INTENT(IN)                            :: density(:, :, :, :)
      REAL(dp), INTENT(IN)                               :: kweights(:)
      TYPE(section_vals_type), POINTER                   :: functional
      REAL(dp), INTENT(IN)                               :: density_cutoff, gradient_cutoff
      REAL(dp), INTENT(OUT)                              :: energies(3)
      COMPLEX(dp), INTENT(OUT)                           :: hamiltonian(:, :, :, :)
      TYPE(paw_potential_type), ALLOCATABLE, INTENT(OUT) :: potential(:)
      LOGICAL, INTENT(IN), OPTIONAL                      :: semilocal_xc
      REAL(dp), INTENT(IN), OPTIONAL                     :: core_energy(:)
      TYPE(paw_potential_type), INTENT(IN), OPTIONAL     :: core_potential(:)

      COMPLEX(dp), ALLOCATABLE                           :: projected(:, :)
      INTEGER                                            :: center, first, ik, ikind, last, n, nao, &
                                                            nk, nproj, ns, spin
      REAL(dp)                                           :: eh, ex, hermitian_error, hermitian_limit
      REAL(dp), ALLOCATABLE                              :: dm(:, :, :), vh(:, :), vxc(:, :, :)

      nk = SIZE(frames)
      CPASSERT(PRESENT(core_energy) .EQV. PRESENT(core_potential))
      IF (PRESENT(core_energy)) THEN
         CPASSERT(SIZE(core_energy) == SIZE(kinds))
         CPASSERT(SIZE(core_potential) == SIZE(kinds))
         CPASSERT(ALL(ieee_is_finite(core_energy)))
      END IF
      CPASSERT(nk > 0 .AND. SIZE(kweights) == nk)
      nao = SIZE(frames(1)%ao%coefficients, 2)
      nproj = SIZE(frames(1)%ao%coefficients, 1)
      ns = SIZE(density, 4)
      CPASSERT(ns == 1 .OR. ns == 2)
      CPASSERT(ALL(SHAPE(density) == [nao, nao, nk, ns]))
      CPASSERT(ALL(SHAPE(hamiltonian) == SHAPE(density)))
      CPASSERT(ALL(ieee_is_finite(kweights)) .AND. ALL(kweights >= 0.0_dp))
      CPASSERT(ABS(SUM(kweights) - 1.0_dp) < 1.e-12_dp)
      DO ik = 1, nk
         CPASSERT(SIZE(frames(ik)%ao%coefficients, 1) == nproj)
         CPASSERT(SIZE(frames(ik)%ao%coefficients, 2) == nao)
         CPASSERT(ALL(frames(ik)%ao%center_atoms == frames(1)%ao%center_atoms))
         CPASSERT(ALL(frames(ik)%ao%center_kinds == frames(1)%ao%center_kinds))
         IF (.NOT. (ALL(frames(ik)%ao%projector_first == frames(1)%ao%projector_first))) THEN
            CPABORT("Inconsistent inverse PAW frame or density")
         END IF
         CPASSERT(ALL(frames(ik)%ao%projector_last == frames(1)%ao%projector_last))
         DO spin = 1, ns
            CPASSERT(ALL(ieee_is_finite(REAL(density(:, :, ik, spin), dp))))
            CPASSERT(ALL(ieee_is_finite(AIMAG(density(:, :, ik, spin)))))
            hermitian_error = MAXVAL(ABS(density(:, :, ik, spin) - CONJG(TRANSPOSE(density(:, :, ik, spin)))))
            hermitian_limit = 1.e-12_dp*(1.0_dp + MAXVAL(ABS(density(:, :, ik, spin))))
            CPASSERT(hermitian_error <= hermitian_limit)
         END DO
      END DO
      ALLOCATE (potential(SIZE(frames(1)%ao%center_atoms)))
      energies = 0.0_dp
      hamiltonian = 0.0_dp
      DO center = 1, SIZE(potential)
         ikind = frames(1)%ao%center_kinds(center)
         CPASSERT(ikind > 0 .AND. ikind <= SIZE(kinds))
         IF (.NOT. kinds(ikind)%active) CYCLE
         first = frames(1)%ao%projector_first(center)
         last = frames(1)%ao%projector_last(center)
         n = last - first + 1
         CPASSERT(n == kinds(ikind)%onecenter%nproj)
         IF (PRESENT(core_potential)) THEN
            CPASSERT(ALLOCATED(core_potential(ikind)%matrix))
            CPASSERT(ALL(SHAPE(core_potential(ikind)%matrix) == [n, n, 1]))
            CPASSERT(ALL(ieee_is_finite(core_potential(ikind)%matrix)))
         END IF
         ALLOCATE (potential(center)%matrix(n, n, ns), SOURCE=0.0_dp)
         ! The global density and inverse coefficients are replicated. One rank owns
         ! each atomic functional, including the core constant, independently of k groups.
         IF (MODULO(center - 1, frames(1)%ao%point_comm%num_pe) == frames(1)%ao%point_comm%mepos) THEN
            ALLOCATE (dm(n, n, ns), projected(n, n), vh(n, n), vxc(n, n, ns))
            dm = 0.0_dp
            DO spin = 1, ns
               DO ik = 1, nk
                  projected(:, :) = MATMUL(frames(ik)%ao%coefficients(first:last, :), &
                                           MATMUL(density(:, :, ik, spin), &
                                                  CONJG(TRANSPOSE(frames(ik)%ao%coefficients(first:last, :)))))
                  dm(:, :, spin) = dm(:, :, spin) + kweights(ik)*REAL(projected, dp)
               END DO
               dm(:, :, spin) = 0.5_dp*(dm(:, :, spin) + TRANSPOSE(dm(:, :, spin)))
            END DO
            CALL paw_onecenter_hartree(kinds(ikind)%onecenter, SUM(dm, DIM=3), eh, vh)
            ex = 0.0_dp
            vxc = 0.0_dp
            IF (PRESENT(semilocal_xc)) THEN
               IF (semilocal_xc) THEN
                  CALL paw_xc_evaluate(kinds(ikind)%xcgrid, kinds(ikind)%onecenter, functional, dm, &
                                       density_cutoff, gradient_cutoff, ex, vxc)
               END IF
            ELSE
               CALL paw_xc_evaluate(kinds(ikind)%xcgrid, kinds(ikind)%onecenter, functional, dm, &
                                    density_cutoff, gradient_cutoff, ex, vxc)
            END IF
            IF (PRESENT(core_energy)) THEN
               ex = ex + core_energy(ikind) + SUM(SUM(dm, DIM=3)*core_potential(ikind)%matrix(:, :, 1))
               DO spin = 1, ns
                  vxc(:, :, spin) = vxc(:, :, spin) + core_potential(ikind)%matrix(:, :, 1)
               END DO
            END IF
            energies = energies + [eh, ex, kinds(ikind)%core_energy]
            DO spin = 1, ns
               potential(center)%matrix(:, :, spin) = vh + vxc(:, :, spin)
            END DO
            DEALLOCATE (dm, projected, vh, vxc)
         END IF
         CALL frames(1)%ao%point_comm%sum(potential(center)%matrix)
         DO spin = 1, ns
            DO ik = 1, nk
               hamiltonian(:, :, ik, spin) = hamiltonian(:, :, ik, spin) + &
                                             MATMUL(CONJG(TRANSPOSE(frames(ik)%ao%coefficients(first:last, :))), &
                                                    MATMUL(CMPLX(potential(center)%matrix(:, :, spin), KIND=dp), &
                                                           frames(ik)%ao%coefficients(first:last, :)))
            END DO
         END DO
      END DO
      CALL frames(1)%ao%point_comm%sum(energies)
   END SUBROUTINE qs_upf_inverse_paw_evaluate

! **************************************************************************************************
!> \brief Differentiate the same PAW functional through the native inverse transformation.
!> \param frames Native frames used for the energy and projector potential.
!> \param density Fixed Gaussian occupation matrices.
!> \param kweights Normalized global k weights.
!> \param potential PAW projector derivative at this density and geometry.
!> \param atom_response Cartesian displacements of the atomic positions.
!> \param cell_response Cell response with fixed fractional grid coordinates.
!> \param reduced_k_response Response of the reduced Bloch vectors, one per global k point.
!> \param derivative Energy response. Rigid atomic quadratures have no additional measure response.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_paw_response(frames, density, kweights, potential, atom_response, cell_response, &
                                          reduced_k_response, derivative)
      TYPE(qs_upf_inverse_grid_type), INTENT(IN)         :: frames(:)
      COMPLEX(dp), INTENT(IN)                            :: density(:, :, :, :)
      REAL(dp), INTENT(IN)                               :: kweights(:)
      TYPE(paw_potential_type), INTENT(IN)               :: potential(:)
      REAL(dp), INTENT(IN)                               :: atom_response(:, :), &
                                                            cell_response(3, 3), &
                                                            reduced_k_response(:, :)
      REAL(dp), INTENT(OUT)                              :: derivative

      COMPLEX(dp), ALLOCATABLE                           :: dw(:, :), projected(:, :)
      INTEGER                                            :: center, first, ik, info, last, nao, nk, &
                                                            nproj, ns, spin
      REAL(dp)                                           :: inverse_cell(3, 3), strain(3, 3), &
                                                            volume_response
      REAL(dp), ALLOCATABLE                              :: points_response(:, :), &
                                                            weights_response(:)

      nk = SIZE(frames)
      CPASSERT(nk > 0 .AND. SIZE(kweights) == nk)
      nao = SIZE(frames(1)%ao%coefficients, 2)
      nproj = SIZE(frames(1)%ao%coefficients, 1)
      ns = SIZE(density, 4)
      CPASSERT(ALL(SHAPE(density) == [nao, nao, nk, ns]))
      CPASSERT(SIZE(potential) == SIZE(frames(1)%ao%center_atoms))
      CPASSERT(ALL(SHAPE(reduced_k_response) == [3, nk]))
      CPASSERT(ALL(ieee_is_finite(kweights)) .AND. ALL(kweights >= 0.0_dp))
      CPASSERT(ABS(SUM(kweights) - 1.0_dp) < 1.e-12_dp)
      ALLOCATE (dw(nproj, nao))
      derivative = 0.0_dp
      DO ik = 1, nk
         CPASSERT(ALL(SHAPE(atom_response) == SHAPE(frames(ik)%ao%positions)))
         inverse_cell = frames(ik)%ao%cell
         CALL invmat(inverse_cell, info)
         CPASSERT(info == 0)
         strain = MATMUL(cell_response, inverse_cell)
         volume_response = strain(1, 1) + strain(2, 2) + strain(3, 3)
         ALLOCATE (points_response(3, SIZE(frames(ik)%ao%weights)), &
                   weights_response(SIZE(frames(ik)%ao%weights)))
         points_response(:, :) = MATMUL(strain, frames(ik)%ao%points)
         weights_response(:) = volume_response*frames(ik)%ao%weights
         CALL qs_upf_inverse_response(frames(ik)%ao, points_response, atom_response, cell_response, &
                                      reduced_k_response(:, ik), weights_response, dw)
         DEALLOCATE (points_response, weights_response)
         DO center = 1, SIZE(potential)
            IF (.NOT. ALLOCATED(potential(center)%matrix)) CYCLE
            first = frames(ik)%ao%projector_first(center)
            last = frames(ik)%ao%projector_last(center)
            CPASSERT(SIZE(potential(center)%matrix, 1) == last - first + 1)
            CPASSERT(SIZE(potential(center)%matrix, 2) == last - first + 1)
            CPASSERT(SIZE(potential(center)%matrix, 3) == ns)
            ALLOCATE (projected(last - first + 1, last - first + 1))
            DO spin = 1, ns
               projected(:, :) = MATMUL(dw(first:last, :), MATMUL(density(:, :, ik, spin), &
                                                                  CONJG(TRANSPOSE(frames(ik)%ao%coefficients(first:last, :)))))
               derivative = derivative + 2.0_dp*kweights(ik)* &
                            SUM(potential(center)%matrix(:, :, spin)*REAL(TRANSPOSE(projected), dp))
            END DO
            DEALLOCATE (projected)
         END DO
      END DO
   END SUBROUTINE qs_upf_inverse_paw_response

END MODULE qs_upf_inverse_paw
