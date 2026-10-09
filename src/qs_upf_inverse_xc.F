!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Mixed-quadrature XC residual in the reconstructed Gaussian representation.
!>        The exact point density uses the defining PW inverse and the original Q functions.
!>        The background is supplied by the caller with its interpolation adjoint retained.
! **************************************************************************************************
MODULE qs_upf_inverse_xc
   USE cell_types,                      ONLY: cell_type
   USE input_section_types,             ONLY: section_get_rval,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE message_passing,                 ONLY: mp_comm_type
   USE pw_grid_types,                   ONLY: PW_MODE_LOCAL,&
                                              pw_grid_type
   USE pw_methods,                      ONLY: pw_derive,&
                                              pw_transfer
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE qs_native_grid_interpolation,    ONLY: interpolate_native_grid_batch,&
                                              native_grid_batch_adjoint,&
                                              native_grid_batch_response
   USE qs_upf_inverse,                  ONLY: qs_upf_inverse_evaluate
   USE qs_upf_inverse_atom_grid,        ONLY: qs_upf_inverse_atom_grid_type
   USE qs_upf_inverse_augmentation,     ONLY: qs_upf_inverse_augmentation_collocate,&
                                              qs_upf_inverse_augmentation_integrate,&
                                              qs_upf_inverse_augmentation_response
   USE qs_upf_inverse_density,          ONLY: qs_upf_inverse_kinetic_collocate,&
                                              qs_upf_inverse_kinetic_operator
   USE qs_upf_inverse_paw,              ONLY: qs_upf_inverse_paw_collocate,&
                                              qs_upf_inverse_paw_operator
   USE qs_upf_inverse_types,            ONLY: qs_upf_inverse_augmentation_type,&
                                              qs_upf_inverse_type,&
                                              qs_upf_inverse_xc_type
   USE skala_gpw_functional,            ONLY: skala_gapw_atom_composite_energy,&
                                              xc_section_uses_native_skala_evaluator
   USE upf_inverse,                     ONLY: upf_inverse_kinetic_type
   USE upf_inverse_quadrature,          ONLY: upf_inverse_projector_density,&
                                              upf_inverse_semilocal_pair,&
                                              upf_inverse_semilocal_response,&
                                              upf_inverse_smooth_density
   USE xc_atom,                         ONLY: fill_rho_set,&
                                              vxc_of_r_new,&
                                              xc_rho_set_atom_update
   USE xc_derivative_set_types,         ONLY: xc_derivative_set_type,&
                                              xc_dset_create,&
                                              xc_dset_get_derivative,&
                                              xc_dset_release
   USE xc_derivative_types,             ONLY: xc_derivative_get,&
                                              xc_derivative_type
   USE xc_derivatives,                  ONLY: xc_functionals_get_needs
   USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_type
   USE xc_rho_set_types,                ONLY: xc_rho_set_create,&
                                              xc_rho_set_release,&
                                              xc_rho_set_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   TYPE qs_upf_inverse_frame_ref_type
      TYPE(qs_upf_inverse_type), POINTER :: frame => NULL()
   END TYPE qs_upf_inverse_frame_ref_type

   PUBLIC :: qs_upf_inverse_frame_ref_type, qs_upf_inverse_xc_type
   PUBLIC :: qs_upf_inverse_xc_evaluate, qs_upf_inverse_xc_response
   PUBLIC :: qs_upf_inverse_xc_grid_evaluate, qs_upf_inverse_xc_grid_response
   PUBLIC :: qs_upf_inverse_xc_pw_fields, qs_upf_inverse_xc_pw_adjoint, qs_upf_inverse_xc_pw_response

CONTAINS

! **************************************************************************************************
!> \brief Evaluate sum(w * (e_XC[n_exact] - e_XC[n_background])) and both density adjoints.
!>        Sum k-resolved densities before the nonlinear XC evaluation. Operators carry no k weight.
!>        The background covectors are negative and unweighted. The caller applies interpolation's
!>        transpose, the same partition weights and the native spectral-gradient transpose.
!> \param frames Borrowed defining-quadrature inverse frames, in global k order.
!> \param augmentation Original real-space charge augmentation.
!> \param matrices Unweighted occupation density matrices, indexed (AO,AO,k,spin).
!> \param kweights Normalized mesh weights, applied once to the point density.
!> \param points Local Cartesian evaluation points, possibly empty.
!> \param weights Local atom quadrature measure including partition weights.
!> \param background_density Interpolated total XC density, including frozen core if present.
!> \param background_gradient Matching physical Cartesian gradients, indexed (point,axis,spin).
!> \param core_density Original frozen-core contribution to the exact point density.
!> \param core_gradient Its physical Cartesian gradient, with the same spin convention.
!> \param xc_section Native XC section supplying the functional and cutoffs.
!> \param lsd Native spin-polarization convention.
!> \param comm Communicator distributing the evaluation points.
!> \param state New local fields and potentials plus globally reduced direct AO operators.
!> \param energy Globally reduced residual energy.
!> \param kinetic Optional original kind-indexed valence kinetic kernels, independent of charge Q.
!> \param background_tau Matching interpolated background kinetic density, including prescribed core.
!> \param core_tau Prescribed frozen-core kinetic density in the exact fields, with the same spin convention.
!> \param atom_grid Complete local atom blocks required by the nonlocal Skala evaluator.
!> \param paw_reconstruction Optional PAW-only reference caches for the composite AE valence fields.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_xc_evaluate(frames, augmentation, matrices, kweights, points, weights, &
                                         background_density, background_gradient, core_density, core_gradient, &
                                         xc_section, lsd, comm, state, energy, kinetic, background_tau, core_tau, atom_grid, &
                                         paw_reconstruction)
      TYPE(qs_upf_inverse_frame_ref_type), INTENT(IN)    :: frames(:)
      TYPE(qs_upf_inverse_augmentation_type), INTENT(IN) :: augmentation
      COMPLEX(dp), INTENT(IN)                            :: matrices(:, :, :, :)
      REAL(dp), INTENT(IN) :: kweights(:), points(:, :), weights(:), background_density(:, :), &
         background_gradient(:, :, :), core_density(:, :), core_gradient(:, :, :)
      TYPE(section_vals_type), POINTER                   :: xc_section
      LOGICAL, INTENT(IN)                                :: lsd
      TYPE(mp_comm_type), INTENT(IN)                     :: comm
      TYPE(qs_upf_inverse_xc_type), INTENT(OUT)          :: state
      REAL(dp), INTENT(OUT)                              :: energy
      TYPE(upf_inverse_kinetic_type), INTENT(IN), &
         OPTIONAL                                        :: kinetic(:)
      REAL(dp), INTENT(IN), OPTIONAL                     :: background_tau(:, :), core_tau(:, :)
      TYPE(qs_upf_inverse_atom_grid_type), INTENT(IN), &
         OPTIONAL                                        :: atom_grid
      TYPE(upf_inverse_kinetic_type), INTENT(IN), &
         OPTIONAL                                        :: paw_reconstruction(:)

      COMPLEX(dp), ALLOCATABLE                           :: kinetic_operator(:, :), projected(:, :), &
                                                            qoperator(:, :)
      INTEGER                                            :: ik, nao, nk, np, nproj, nspin, spin
      REAL(dp), ALLOCATABLE :: ebackground(:), gradient(:, :), integration_weights(:), &
         reconstruction_density(:), reconstruction_gradient(:, :), rho(:), tau(:), &
         tau_augmentation(:)
      TYPE(qs_upf_inverse_type), POINTER                 :: frame

      nk = SIZE(frames)
      CPASSERT(nk > 0 .AND. SIZE(kweights) == nk)
      CPASSERT(ALL(kweights >= 0.0_dp))
      CPASSERT(ABS(SUM(kweights) - 1.0_dp) <= 128.0_dp*EPSILON(1.0_dp)*nk)
      CPASSERT(ASSOCIATED(frames(1)%frame))
      nao = SIZE(frames(1)%frame%coefficients, 2)
      np = SIZE(points, 2)
      nspin = SIZE(matrices, 4)
      CPASSERT(ALL(SHAPE(matrices) == [nao, nao, nk, nspin]))
      CPASSERT(SIZE(points, 1) == 3 .AND. SIZE(weights) == np)
      CPASSERT(ALL(SHAPE(background_density) == [np, nspin]))
      CPASSERT(ALL(SHAPE(background_gradient) == [np, 3, nspin]))
      CPASSERT(ALL(SHAPE(core_density) == [np, nspin]))
      CPASSERT(ALL(SHAPE(core_gradient) == [np, 3, nspin]))
      state%kinetic = PRESENT(kinetic)
      state%paw_reconstruction = PRESENT(paw_reconstruction)
      CPASSERT(state%kinetic .EQV. PRESENT(background_tau))
      CPASSERT(state%kinetic .EQV. PRESENT(core_tau))
      IF (state%kinetic) THEN
         CPASSERT(ALL(SHAPE(background_tau) == [np, nspin]))
         CPASSERT(ALL(SHAPE(core_tau) == [np, nspin]))
      END IF
      ALLOCATE (state%soft(np, nao, nk), state%gradient(np, nao, 3, nk), &
                state%operators(nao, nao, nk, nspin), state%density(np, nspin), &
                state%density_gradient(np, 3, nspin), state%energy_density(np), &
                state%vdensity(np, nspin), state%vgradient(np, 3, nspin), &
                state%background_vdensity(np, nspin), state%background_vgradient(np, 3, nspin), &
                state%tau(np, nspin), state%vtau(np, nspin), state%background_vtau(np, nspin), &
                rho(np), gradient(np, 3), tau(np), tau_augmentation(np), ebackground(np))
      IF (state%paw_reconstruction) ALLOCATE (reconstruction_density(np), reconstruction_gradient(np, 3))
      state%density(:, :) = core_density
      state%density_gradient(:, :, :) = core_gradient
      state%tau = 0.0_dp
      IF (state%kinetic) state%tau(:, :) = core_tau
      DO ik = 1, nk
         CPASSERT(ASSOCIATED(frames(ik)%frame))
         frame => frames(ik)%frame
         CPASSERT(SIZE(frame%coefficients, 2) == nao)
         CALL qs_upf_inverse_evaluate(frame, points, state%soft(:, :, ik), state%gradient(:, :, :, ik))
         nproj = SIZE(frame%coefficients, 1)
         ALLOCATE (projected(nproj, nproj))
         DO spin = 1, nspin
            CALL upf_inverse_smooth_density(state%soft(:, :, ik), state%gradient(:, :, :, ik), &
                                            matrices(:, :, ik, spin), rho, gradient, tau)
            CALL upf_inverse_projector_density(frame%coefficients, matrices(:, :, ik, spin), projected)
            CALL qs_upf_inverse_augmentation_collocate(frame, augmentation, projected, rho, gradient, points)
            IF (state%paw_reconstruction) THEN
               CALL qs_upf_inverse_paw_collocate(frame, paw_reconstruction, augmentation, points, matrices(:, :, ik, spin), &
                                                 reconstruction_density, reconstruction_gradient)
               rho(:) = rho + reconstruction_density
               gradient(:, :) = gradient + reconstruction_gradient
            END IF
            state%density(:, spin) = state%density(:, spin) + kweights(ik)*rho
            state%density_gradient(:, :, spin) = state%density_gradient(:, :, spin) + kweights(ik)*gradient
            IF (state%kinetic) THEN
               CALL qs_upf_inverse_kinetic_collocate(frame, kinetic, points, matrices(:, :, ik, spin), tau_augmentation)
               IF (state%paw_reconstruction) THEN
                  CALL qs_upf_inverse_kinetic_collocate(frame, paw_reconstruction, points, matrices(:, :, ik, spin), &
                                                        reconstruction_density)
                  tau_augmentation(:) = tau_augmentation + reconstruction_density
               END IF
               state%tau(:, spin) = state%tau(:, spin) + kweights(ik)*(tau + tau_augmentation)
            END IF
         END DO
         DEALLOCATE (projected)
      END DO
      ALLOCATE (integration_weights(np))
      integration_weights(:) = weights
      IF (xc_section_uses_native_skala_evaluator(xc_section)) THEN
         CPASSERT(PRESENT(atom_grid) .AND. state%kinetic)
         CALL evaluate_skala_residual(xc_section, comm, state, background_density, background_gradient, &
                                      background_tau, points, weights, atom_grid, frames(1)%frame%positions, energy)
         ! Skala supplies derivatives of the complete weighted energy. Do not divide by
         ! quadrature weights, including at zero-weight rows used by its nonlocal descriptors.
         integration_weights = 1.0_dp
      ELSE
         CALL evaluate_xc_fields(xc_section, lsd, state%density, state%density_gradient, &
                                 state%energy_density, state%vdensity, state%vgradient, &
                                 kinetic_density=state%tau, vkinetic=state%vtau, kinetic_available=state%kinetic)
         CALL evaluate_xc_fields(xc_section, lsd, background_density, background_gradient, &
                                 ebackground, state%background_vdensity, state%background_vgradient, &
                                 kinetic_density=background_tau, vkinetic=state%background_vtau)
         state%energy_density(:) = state%energy_density - ebackground
         state%background_vdensity = -state%background_vdensity
         state%background_vgradient = -state%background_vgradient
         state%background_vtau = -state%background_vtau
         energy = SUM(weights*state%energy_density)
         CALL comm%sum(energy)
      END IF
      DO ik = 1, nk
         frame => frames(ik)%frame
         nproj = SIZE(frame%coefficients, 1)
         ALLOCATE (qoperator(nproj, nproj), kinetic_operator(nao, nao))
         DO spin = 1, nspin
            CALL upf_inverse_semilocal_pair(state%soft(:, :, ik), state%soft(:, :, ik), &
                                            state%gradient(:, :, :, ik), state%gradient(:, :, :, ik), integration_weights, &
                                            state%vdensity(:, spin), state%vgradient(:, :, spin), state%vtau(:, spin), &
                                            state%operators(:, :, ik, spin), comm)
            CALL qs_upf_inverse_augmentation_integrate(frame, augmentation, state%vdensity(:, spin), &
                                                       state%vgradient(:, :, spin), qoperator, points, integration_weights, comm)
            state%operators(:, :, ik, spin) = state%operators(:, :, ik, spin) + &
                                              MATMUL(CONJG(TRANSPOSE(frame%coefficients)), MATMUL(qoperator, frame%coefficients))
            IF (state%paw_reconstruction) THEN
               CALL qs_upf_inverse_paw_operator(frame, paw_reconstruction, augmentation, points, integration_weights, &
                                                state%vdensity(:, spin), state%vgradient(:, :, spin), kinetic_operator, comm)
               state%operators(:, :, ik, spin) = state%operators(:, :, ik, spin) + kinetic_operator
            END IF
            IF (state%kinetic) THEN
               CALL qs_upf_inverse_kinetic_operator(frame, kinetic, points, integration_weights, state%vtau(:, spin), &
                                                    kinetic_operator, comm)
               state%operators(:, :, ik, spin) = state%operators(:, :, ik, spin) + kinetic_operator
               IF (state%paw_reconstruction) THEN
                  CALL qs_upf_inverse_kinetic_operator(frame, paw_reconstruction, points, integration_weights, &
                                                       state%vtau(:, spin), kinetic_operator, comm)
                  state%operators(:, :, ik, spin) = state%operators(:, :, ik, spin) + kinetic_operator
               END IF
            END IF
         END DO
         DEALLOCATE (qoperator, kinetic_operator)
      END DO
   END SUBROUTINE qs_upf_inverse_xc_evaluate

! **************************************************************************************************
!> \brief Difference of complete Skala functionals on the same finite atom quadrature.
!>        Descriptor adjoints include integration weights. The radial base weights are fixed
!>        under geometry changes. Point, partition-weight and atom-coordinate responses remain.
!> \param xc_section ...
!> \param comm ...
!> \param state Exact fields and resulting weighted adjoints.
!> \param background_density ...
!> \param background_gradient ...
!> \param background_tau ...
!> \param points ...
!> \param weights ...
!> \param atom_grid Complete local atom blocks, including fixed radial base weights.
!> \param positions All original-coordinate atomic positions.
!> \param energy Globally reduced difference of the two complete functional evaluations.
! **************************************************************************************************
   SUBROUTINE evaluate_skala_residual(xc_section, comm, state, background_density, background_gradient, &
                                      background_tau, points, weights, atom_grid, positions, energy)
      TYPE(section_vals_type), POINTER                   :: xc_section
      TYPE(mp_comm_type), INTENT(IN)                     :: comm
      TYPE(qs_upf_inverse_xc_type), INTENT(INOUT)        :: state
      REAL(dp), INTENT(IN)                               :: background_density(:, :), &
                                                            background_gradient(:, :, :), &
                                                            background_tau(:, :), points(:, :), &
                                                            weights(:)
      TYPE(qs_upf_inverse_atom_grid_type), INTENT(IN)    :: atom_grid
      REAL(dp), INTENT(IN)                               :: positions(:, :)
      REAL(dp), INTENT(OUT)                              :: energy

      INTEGER                                            :: atom, local_atom, np, ns, row, spin
      INTEGER(int_8), ALLOCATABLE                        :: sizes(:)
      INTEGER, ALLOCATABLE                               :: local_atoms(:)
      REAL(dp)                                           :: background_energy, scale
      REAL(dp), ALLOCATABLE :: atom_coords(:, :), atom_derivative(:, :), &
         background_atom_derivative(:, :), background_base_derivative(:), &
         background_point_derivative(:, :), background_weight_derivative(:), base_derivative(:), &
         base_weights(:), coords(:, :), density(:, :), gradient(:, :, :), grid_weights(:), &
         point_derivative(:, :), tau(:, :), vdensity(:, :), vgradient(:, :, :), vtau(:, :), &
         weight_derivative(:)

      np = SIZE(weights)
      ns = SIZE(state%density, 2)
      CPASSERT(ns == 1 .OR. ns == 2)
      CPASSERT(SIZE(atom_grid%atoms) == np)
      CPASSERT(ALL(SHAPE(atom_grid%points) == SHAPE(points)))
      CPASSERT(SIZE(positions, 2) == SIZE(atom_grid%layouts))
      local_atom = 0
      DO atom = 1, SIZE(positions, 2)
         IF (ANY(atom_grid%atoms == atom)) local_atom = local_atom + 1
      END DO
      ALLOCATE (local_atoms(local_atom), sizes(local_atom), atom_coords(3, local_atom), &
                density(np, 2), gradient(np, 3, 2), tau(np, 2), coords(3, np), &
                grid_weights(np), base_weights(np))
      row = 0
      local_atom = 0
      DO atom = 1, SIZE(positions, 2)
         IF (.NOT. ANY(atom_grid%atoms == atom)) CYCLE
         local_atom = local_atom + 1
         local_atoms(local_atom) = atom
         sizes(local_atom) = INT(COUNT(atom_grid%atoms == atom), int_8)
         CPASSERT(ALL(atom_grid%atoms(row + 1:row + INT(sizes(local_atom))) == atom))
         row = row + INT(sizes(local_atom))
         atom_coords(:, local_atom) = positions(:, atom)
      END DO
      CPASSERT(row == np)
      coords(:, :) = points
      grid_weights(:) = weights
      base_weights(:) = atom_grid%base_weights
      scale = 1.0_dp
      IF (ns == 1) scale = 0.5_dp
      DO spin = 1, 2
         density(:, spin) = scale*state%density(:, MIN(spin, ns))
         gradient(:, :, spin) = scale*state%density_gradient(:, :, MIN(spin, ns))
         tau(:, spin) = scale*state%tau(:, MIN(spin, ns))
      END DO
      CALL skala_gapw_atom_composite_energy(xc_section, comm, density, gradient, tau, coords, grid_weights, &
                                            base_weights, sizes, atom_coords, energy, vdensity, vgradient, vtau, &
                                            point_derivative, weight_derivative, base_derivative, atom_derivative)
      DO spin = 1, ns
         IF (ns == 1) THEN
            state%vdensity(:, spin) = 0.5_dp*SUM(vdensity, DIM=2)
            state%vgradient(:, :, spin) = 0.5_dp*SUM(vgradient, DIM=3)
            state%vtau(:, spin) = 0.5_dp*SUM(vtau, DIM=2)
         ELSE
            state%vdensity(:, spin) = vdensity(:, spin)
            state%vgradient(:, :, spin) = vgradient(:, :, spin)
            state%vtau(:, spin) = vtau(:, spin)
         END IF
      END DO
      DO spin = 1, 2
         density(:, spin) = scale*background_density(:, MIN(spin, ns))
         gradient(:, :, spin) = scale*background_gradient(:, :, MIN(spin, ns))
         tau(:, spin) = scale*background_tau(:, MIN(spin, ns))
      END DO
      CALL skala_gapw_atom_composite_energy(xc_section, comm, density, gradient, tau, coords, grid_weights, &
                                            base_weights, sizes, atom_coords, background_energy, vdensity, vgradient, vtau, &
                                            background_point_derivative, background_weight_derivative, &
                                            background_base_derivative, background_atom_derivative)
      DO spin = 1, ns
         IF (ns == 1) THEN
            state%background_vdensity(:, spin) = -0.5_dp*SUM(vdensity, DIM=2)
            state%background_vgradient(:, :, spin) = -0.5_dp*SUM(vgradient, DIM=3)
            state%background_vtau(:, spin) = -0.5_dp*SUM(vtau, DIM=2)
         ELSE
            state%background_vdensity(:, spin) = -vdensity(:, spin)
            state%background_vgradient(:, :, spin) = -vgradient(:, :, spin)
            state%background_vtau(:, spin) = -vtau(:, spin)
         END IF
      END DO
      state%weighted_adjoint = .TRUE.
      state%energy_density = 0.0_dp
      energy = energy - background_energy
      ALLOCATE (state%point_derivative(3, np), state%weight_derivative(np), state%atom_derivative(3, SIZE(positions, 2)))
      state%point_derivative(:, :) = point_derivative - background_point_derivative
      state%weight_derivative(:) = weight_derivative - background_weight_derivative
      state%atom_derivative = 0.0_dp
      DO atom = 1, SIZE(local_atoms)
         state%atom_derivative(:, local_atoms(atom)) = atom_derivative(:, atom) - background_atom_derivative(:, atom)
      END DO
   END SUBROUTINE evaluate_skala_residual

! **************************************************************************************************
!> \brief Joint geometry response of the same XC residual at fixed occupation density matrices.
!>        Use e_exact-e_background for the quadrature response, rather than the potential-density
!>        contraction. The defining PW quadrature supplies dW, including its own measure response.
!> \param frames ...
!> \param augmentation ...
!> \param matrices ...
!> \param kweights Fixed mesh weights.
!> \param points ...
!> \param weights ...
!> \param state Fields from the unperturbed evaluation.
!> \param point_response Motion of the local atom-grid points.
!> \param weights_response Response of the atom quadrature and partition measure.
!> \param atom_response Motion of every Gaussian and projector center.
!> \param cell_response Lattice response.
!> \param reduced_k_response Reduced Bloch-vector response, indexed (axis,k).
!> \param coefficients_response dW from each frame's original defining quadrature.
!> \param background_response Complete interpolation and source-field density response.
!> \param background_gradient_response Matching Cartesian-gradient response.
!> \param core_response Original frozen-core density response, with the same spin convention.
!> \param core_gradient_response Matching frozen-core gradient response.
!> \param comm ...
!> \param response Globally reduced XC energy response.
!> \param kinetic Same independent kind-indexed kernels as in the unperturbed evaluation.
!> \param background_tau_response Complete background kinetic-density response.
!> \param core_tau_response Prescribed frozen-core kinetic-density response.
!> \param paw_reconstruction Same optional PAW-only reconstruction caches as in the evaluation.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_xc_response(frames, augmentation, matrices, kweights, points, weights, state, &
                                         point_response, weights_response, atom_response, cell_response, &
                                         reduced_k_response, coefficients_response, background_response, &
                                         background_gradient_response, core_response, core_gradient_response, comm, response, &
                                         kinetic, background_tau_response, core_tau_response, paw_reconstruction)
      TYPE(qs_upf_inverse_frame_ref_type), INTENT(IN)    :: frames(:)
      TYPE(qs_upf_inverse_augmentation_type), INTENT(IN) :: augmentation
      COMPLEX(dp), INTENT(IN)                            :: matrices(:, :, :, :)
      REAL(dp), INTENT(IN)                               :: kweights(:), points(:, :), weights(:)
      TYPE(qs_upf_inverse_xc_type), INTENT(IN)           :: state
      REAL(dp), INTENT(IN) :: point_response(:, :), weights_response(:), atom_response(:, :), &
         cell_response(3, 3), reduced_k_response(:, :)
      COMPLEX(dp), INTENT(IN)                            :: coefficients_response(:, :, :)
      REAL(dp), INTENT(IN) :: background_response(:, :), background_gradient_response(:, :, :), &
         core_response(:, :), core_gradient_response(:, :, :)
      TYPE(mp_comm_type), INTENT(IN)                     :: comm
      REAL(dp), INTENT(OUT)                              :: response
      TYPE(upf_inverse_kinetic_type), INTENT(IN), &
         OPTIONAL                                        :: kinetic(:)
      REAL(dp), INTENT(IN), OPTIONAL                     :: background_tau_response(:, :), &
                                                            core_tau_response(:, :)
      TYPE(upf_inverse_kinetic_type), INTENT(IN), &
         OPTIONAL                                        :: paw_reconstruction(:)

      COMPLEX(dp), ALLOCATABLE :: dgradient(:, :, :), dqoperator(:, :), dsoft(:, :), &
         gradient(:, :, :), matrix(:, :), qoperator(:, :), soft(:, :)
      INTEGER                                            :: ik, nao, np, nproj, spin
      REAL(dp)                                           :: kinetic_response
      REAL(dp), ALLOCATABLE :: dtau(:), integration_weights(:), reconstruction_density(:), &
         reconstruction_gradient(:, :), reconstruction_gradient_response(:, :), &
         reconstruction_response(:), tau(:), zero(:), zeros(:, :)
      TYPE(qs_upf_inverse_type), POINTER                 :: frame

      np = SIZE(points, 2)
      nao = SIZE(matrices, 1)
      CPASSERT(ALL(SHAPE(point_response) == SHAPE(points)))
      CPASSERT(SIZE(weights_response) == np .AND. SIZE(weights) == np)
      CPASSERT(ALL(SHAPE(background_response) == SHAPE(state%vdensity)))
      CPASSERT(ALL(SHAPE(core_response) == SHAPE(state%vdensity)))
      CPASSERT(ALL(SHAPE(background_gradient_response) == SHAPE(state%vgradient)))
      CPASSERT(ALL(SHAPE(core_gradient_response) == SHAPE(state%vgradient)))
      CPASSERT(ALL(SHAPE(reduced_k_response) == [3, SIZE(frames)]))
      CPASSERT(state%kinetic .EQV. PRESENT(kinetic))
      CPASSERT(state%kinetic .EQV. PRESENT(background_tau_response))
      CPASSERT(state%kinetic .EQV. PRESENT(core_tau_response))
      CPASSERT(state%paw_reconstruction .EQV. PRESENT(paw_reconstruction))
      IF (state%kinetic) THEN
         CPASSERT(ALL(SHAPE(background_tau_response) == SHAPE(state%vtau)))
         CPASSERT(ALL(SHAPE(core_tau_response) == SHAPE(state%vtau)))
      END IF
      ALLOCATE (soft(np, nao), gradient(np, nao, 3), dsoft(np, nao), dgradient(np, nao, 3), &
                matrix(nao, nao), tau(np), dtau(np), zero(np), zeros(np, 3))
      IF (state%paw_reconstruction) ALLOCATE (reconstruction_density(np), reconstruction_gradient(np, 3), &
                                              reconstruction_response(np), reconstruction_gradient_response(np, 3))
      zero = 0.0_dp
      zeros = 0.0_dp
      ALLOCATE (integration_weights(np))
      integration_weights(:) = weights
      IF (state%weighted_adjoint) THEN
         integration_weights = 1.0_dp
         response = SUM(state%weight_derivative*weights_response) + SUM(state%point_derivative*point_response) + &
                    SUM(state%atom_derivative*atom_response)
      ELSE
         response = SUM(weights_response*state%energy_density)
      END IF
      DO spin = 1, SIZE(matrices, 4)
         response = response + SUM(integration_weights*(state%vdensity(:, spin)*core_response(:, spin) + &
                                                      SUM(state%vgradient(:, :, spin)*core_gradient_response(:, :, spin), DIM=2) + &
                                                        state%background_vdensity(:, spin)*background_response(:, spin) + &
                                       SUM(state%background_vgradient(:, :, spin)*background_gradient_response(:, :, spin), DIM=2)))
         IF (state%kinetic) response = response + SUM(integration_weights*(state%vtau(:, spin)*core_tau_response(:, spin) + &
                                                                   state%background_vtau(:, spin)*background_tau_response(:, spin)))
      END DO
      CALL comm%sum(response)
      DO ik = 1, SIZE(frames)
         frame => frames(ik)%frame
         nproj = SIZE(frame%coefficients, 1)
         ALLOCATE (qoperator(nproj, nproj), dqoperator(nproj, nproj))
         CALL qs_upf_inverse_evaluate(frame, points, soft, gradient, dsoft, dgradient, point_response, atom_response, &
                                      cell_response, reduced_k_response(:, ik), coefficients_response(:, :, ik))
         DO spin = 1, SIZE(matrices, 4)
            CALL upf_inverse_semilocal_response(state%soft(:, :, ik), state%soft(:, :, ik), &
                                                state%gradient(:, :, :, ik), state%gradient(:, :, :, ik), integration_weights, &
                    state%vdensity(:, spin), state%vgradient(:, :, spin), state%vtau(:, spin), dsoft, dsoft, dgradient, dgradient, &
                                                zero, zero, zeros, zero, matrix, comm)
            CALL qs_upf_inverse_augmentation_integrate(frame, augmentation, state%vdensity(:, spin), &
                                                       state%vgradient(:, :, spin), qoperator, points, integration_weights, comm)
            CALL qs_upf_inverse_augmentation_response(frame, augmentation, state%vdensity(:, spin), &
                                     state%vgradient(:, :, spin), point_response, atom_response, cell_response, zero, zero, zeros, &
                                                      dqoperator, points, integration_weights, comm)
            matrix(:, :) = matrix + &
                           MATMUL(CONJG(TRANSPOSE(coefficients_response(:, :, ik))), MATMUL(qoperator, frame%coefficients)) + &
                           MATMUL(CONJG(TRANSPOSE(frame%coefficients)), MATMUL(dqoperator, frame%coefficients)) + &
                           MATMUL(CONJG(TRANSPOSE(frame%coefficients)), MATMUL(qoperator, coefficients_response(:, :, ik)))
            response = response + kweights(ik)*REAL(SUM(matrices(:, :, ik, spin)*TRANSPOSE(matrix)), dp)
            IF (state%paw_reconstruction) THEN
               CALL qs_upf_inverse_paw_collocate(frame, paw_reconstruction, augmentation, points, matrices(:, :, ik, spin), &
                                                 reconstruction_density, reconstruction_gradient, reconstruction_response, &
                                                 reconstruction_gradient_response, point_response, atom_response, cell_response, &
                                                 coefficients_response(:, :, ik))
               kinetic_response = SUM(integration_weights*(state%vdensity(:, spin)*reconstruction_response + &
                                                          SUM(state%vgradient(:, :, spin)*reconstruction_gradient_response, DIM=2)))
               CALL comm%sum(kinetic_response)
               response = response + kweights(ik)*kinetic_response
            END IF
            IF (state%kinetic) THEN
               CALL qs_upf_inverse_kinetic_collocate(frame, kinetic, points, matrices(:, :, ik, spin), tau, dtau, &
                                                     point_response, atom_response, cell_response, coefficients_response(:, :, ik))
               kinetic_response = SUM(integration_weights*state%vtau(:, spin)*dtau)
               CALL comm%sum(kinetic_response)
               response = response + kweights(ik)*kinetic_response
               IF (state%paw_reconstruction) THEN
                  CALL qs_upf_inverse_kinetic_collocate(frame, paw_reconstruction, points, matrices(:, :, ik, spin), tau, dtau, &
                                                      point_response, atom_response, cell_response, coefficients_response(:, :, ik))
                  kinetic_response = SUM(integration_weights*state%vtau(:, spin)*dtau)
                  CALL comm%sum(kinetic_response)
                  response = response + kweights(ik)*kinetic_response
               END IF
            END IF
         END DO
         DEALLOCATE (qoperator, dqoperator)
      END DO
   END SUBROUTINE qs_upf_inverse_xc_response

! **************************************************************************************************
!> \brief Native semilocal XC evaluation with unweighted density, gradient and kinetic covectors.
!> \param xc_section ...
!> \param lsd ...
!> \param density ...
!> \param gradient ...
!> \param energy_density ...
!> \param vdensity ...
!> \param vgradient ...
!> \param kinetic_density Optional independent kinetic density, with factor one half.
!> \param vkinetic Optional potential conjugate to kinetic density.
!> \param kinetic_available Optional explicit availability flag for allocated zero placeholders.
! **************************************************************************************************
   SUBROUTINE evaluate_xc_fields(xc_section, lsd, density, gradient, energy_density, vdensity, vgradient, &
                                 kinetic_density, vkinetic, kinetic_available)
      TYPE(section_vals_type), POINTER                   :: xc_section
      LOGICAL, INTENT(IN)                                :: lsd
      REAL(dp), INTENT(IN)                               :: density(:, :), gradient(:, :, :)
      REAL(dp), INTENT(OUT)                              :: energy_density(:), vdensity(:, :), &
                                                            vgradient(:, :, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: kinetic_density(:, :)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: vkinetic(:, :)
      LOGICAL, INTENT(IN), OPTIONAL                      :: kinetic_available

      INTEGER                                            :: bounds(2, 3), np, nspin, spin
      LOGICAL                                            :: have_kinetic
      REAL(dp)                                           :: energy
      REAL(dp), ALLOCATABLE                              :: drho(:, :, :, :), rho(:, :, :), &
                                                            tau(:, :, :), weights(:, :)
      REAL(dp), POINTER                                  :: derivative(:, :, :), vtau(:, :, :), &
                                                            vxc(:, :, :), vxg(:, :, :, :)
      TYPE(section_vals_type), POINTER                   :: functional
      TYPE(xc_derivative_set_type)                       :: derivatives
      TYPE(xc_derivative_type), POINTER                  :: energy_term
      TYPE(xc_rho_cflags_type)                           :: needs
      TYPE(xc_rho_set_type)                              :: fields

      np = SIZE(density, 1)
      nspin = SIZE(density, 2)
      CPASSERT((nspin == 1 .AND. .NOT. lsd) .OR. (nspin == 2 .AND. lsd))
      functional => section_vals_get_subs_vals(xc_section, "XC_FUNCTIONAL")
      needs = xc_functionals_get_needs(functional, lsd=lsd, calc_potential=.TRUE.)
      have_kinetic = PRESENT(kinetic_density)
      IF (PRESENT(kinetic_available)) have_kinetic = have_kinetic .AND. kinetic_available
      IF ((needs%tau .OR. needs%tau_spin) .AND. .NOT. have_kinetic) THEN
         CPABORT("Missing independent valence or core kinetic density for inverse XC")
      END IF
      IF (needs%laplace_rho .OR. needs%laplace_rho_spin) THEN
         CPABORT("Inverse mixed-quadrature XC needs independent Laplacian augmentation")
      END IF
      energy_density = 0.0_dp
      vdensity = 0.0_dp
      vgradient = 0.0_dp
      IF (PRESENT(vkinetic)) vkinetic = 0.0_dp
      IF (PRESENT(kinetic_density)) THEN
         CPASSERT(ALL(SHAPE(kinetic_density) == [np, nspin]))
      END IF
      IF (PRESENT(vkinetic)) THEN
         CPASSERT(ALL(SHAPE(vkinetic) == [np, nspin]))
      END IF
      IF (np == 0) RETURN
      bounds(1, :) = 1
      bounds(2, :) = [np, 1, 1]
      ALLOCATE (rho(np, 1, nspin), drho(4, np, 1, nspin), tau(np, 1, nspin), weights(np, 1))
      ALLOCATE (vxc(np, 1, nspin), vxg(3, np, 1, nspin), vtau(np, 1, nspin))
      NULLIFY (derivative)
      rho(:, 1, :) = density
      tau = 0.0_dp
      IF (have_kinetic) tau(:, 1, :) = kinetic_density
      weights = 1.0_dp
      vxc = 0.0_dp
      vxg = 0.0_dp
      vtau = 0.0_dp
      DO spin = 1, nspin
         drho(1:3, :, 1, spin) = TRANSPOSE(gradient(:, :, spin))
         drho(4, :, 1, spin) = SQRT(SUM(gradient(:, :, spin)**2, DIM=2))
      END DO
      CALL xc_rho_set_create(fields, bounds, rho_cutoff=section_get_rval(xc_section, "DENSITY_CUTOFF"), &
                             drho_cutoff=section_get_rval(xc_section, "GRADIENT_CUTOFF"), &
                             tau_cutoff=section_get_rval(xc_section, "TAU_CUTOFF"))
      CALL xc_rho_set_atom_update(fields, needs, nspin, bounds)
      CALL fill_rho_set(fields, lsd, nspin, needs, rho, drho, tau, np, 1)
      CALL xc_dset_create(derivatives, local_bounds=bounds)
      CALL vxc_of_r_new(functional, fields, derivatives, 1, needs, weights, lsd, np, 1, energy, vxc, vxg, vtau)
      energy_term => xc_dset_get_derivative(derivatives, [INTEGER::])
      IF (ASSOCIATED(energy_term)) THEN
         CALL xc_derivative_get(energy_term, deriv_data=derivative)
         energy_density = derivative(:, 1, 1)
      END IF
      vdensity = vxc(:, 1, :)
      IF (PRESENT(vkinetic)) vkinetic(:, :) = vtau(:, 1, :)
      DO spin = 1, nspin
         vgradient(:, :, spin) = TRANSPOSE(vxg(:, :, 1, spin))
      END DO
      CALL xc_dset_release(derivatives)
      CALL xc_rho_set_release(fields)
      DEALLOCATE (vxc, vxg, vtau)
   END SUBROUTINE evaluate_xc_fields

! **************************************************************************************************
!> \brief Evaluate the inverse XC residual against interpolated native PW descriptors.
!> \param frames ...
!> \param augmentation ...
!> \param matrices ...
!> \param kweights ...
!> \param points ...
!> \param weights ...
!> \param fields Replicated native XC density and physical Cartesian gradient including frozen core.
!>               Descriptors per spin are (density,gx,gy,gz), with tau last when kinetic is supplied.
!> \param grid Native source grid.
!> \param cell Native source periodicity.
!> \param core_density ...
!> \param core_gradient ...
!> \param xc_section ...
!> \param lsd ...
!> \param comm Point-batch communicator.
!> \param state ...
!> \param energy ...
!> \param grid_adjoint Globally reduced negative background descriptor transpose including atom weights.
!>                     Fold gradient covectors with the native spectral transpose and divide by the
!>                     PW volume element before integrating an unweighted scalar potential.
!> \param kinetic Optional independent kind-indexed valence kinetic kernels.
!> \param core_tau Prescribed exact frozen-core kinetic density when kinetic is supplied.
!> \param atom_grid Complete local atom blocks for nonlocal Skala evaluation.
!> \param paw_reconstruction Optional PAW-only reconstruction caches for the composite fields.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_xc_grid_evaluate(frames, augmentation, matrices, kweights, points, weights, &
                                              fields, grid, cell, core_density, core_gradient, xc_section, &
                                           lsd, comm, state, energy, grid_adjoint, kinetic, core_tau, atom_grid, paw_reconstruction)
      TYPE(qs_upf_inverse_frame_ref_type), INTENT(IN)    :: frames(:)
      TYPE(qs_upf_inverse_augmentation_type), INTENT(IN) :: augmentation
      COMPLEX(dp), INTENT(IN)                            :: matrices(:, :, :, :)
      REAL(dp), INTENT(IN)                               :: kweights(:), points(:, :), weights(:), &
                                                            fields(:, :, :, :)
      TYPE(pw_grid_type), INTENT(IN), POINTER            :: grid
      TYPE(cell_type), INTENT(IN), POINTER               :: cell
      REAL(dp), INTENT(IN)                               :: core_density(:, :), &
                                                            core_gradient(:, :, :)
      TYPE(section_vals_type), POINTER                   :: xc_section
      LOGICAL, INTENT(IN)                                :: lsd
      TYPE(mp_comm_type), INTENT(IN)                     :: comm
      TYPE(qs_upf_inverse_xc_type), INTENT(OUT)          :: state
      REAL(dp), INTENT(OUT)                              :: energy
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: grid_adjoint(:, :, :, :)
      TYPE(upf_inverse_kinetic_type), INTENT(IN), &
         OPTIONAL                                        :: kinetic(:)
      REAL(dp), INTENT(IN), OPTIONAL                     :: core_tau(:, :)
      TYPE(qs_upf_inverse_atom_grid_type), INTENT(IN), &
         OPTIONAL                                        :: atom_grid
      TYPE(upf_inverse_kinetic_type), INTENT(IN), &
         OPTIONAL                                        :: paw_reconstruction(:)

      INTEGER                                            :: axis, np, ns, spin, stride
      REAL(dp), ALLOCATABLE                              :: background(:, :), covectors(:, :), &
                                                            gradient(:, :, :), &
                                                            integration_weights(:), tau(:, :), &
                                                            values(:, :)

      np = SIZE(points, 2)
      ns = SIZE(matrices, 4)
      CPASSERT(PRESENT(kinetic) .EQV. PRESENT(core_tau))
      stride = 4
      IF (PRESENT(kinetic)) stride = 5
      CPASSERT(SIZE(fields, 4) == stride*ns)
      ALLOCATE (values(np, stride*ns), covectors(np, stride*ns), background(np, ns), gradient(np, 3, ns))
      IF (PRESENT(kinetic)) ALLOCATE (tau(np, ns))
      CALL interpolate_native_grid_batch(fields, grid, cell, points, values, wrap_auxiliary_cell=.TRUE., continuous=.TRUE.)
      DO spin = 1, ns
         background(:, spin) = values(:, stride*(spin - 1) + 1)
         DO axis = 1, 3
            gradient(:, axis, spin) = values(:, stride*(spin - 1) + 1 + axis)
         END DO
         IF (PRESENT(kinetic)) tau(:, spin) = values(:, stride*spin)
      END DO
      CALL qs_upf_inverse_xc_evaluate(frames, augmentation, matrices, kweights, points, weights, &
                                      background, gradient, core_density, core_gradient, xc_section, &
                                      lsd, comm, state, energy, kinetic, tau, core_tau, atom_grid, paw_reconstruction)
      ALLOCATE (integration_weights(np))
      integration_weights(:) = weights
      IF (state%weighted_adjoint) integration_weights = 1.0_dp
      DO spin = 1, ns
         covectors(:, stride*(spin - 1) + 1) = integration_weights*state%background_vdensity(:, spin)
         DO axis = 1, 3
            covectors(:, stride*(spin - 1) + 1 + axis) = integration_weights*state%background_vgradient(:, axis, spin)
         END DO
         IF (PRESENT(kinetic)) covectors(:, stride*spin) = integration_weights*state%background_vtau(:, spin)
      END DO
      ALLOCATE (grid_adjoint(SIZE(fields, 1), SIZE(fields, 2), SIZE(fields, 3), SIZE(fields, 4)))
      grid_adjoint = 0.0_dp
      CALL native_grid_batch_adjoint(covectors, grid, cell, points, grid_adjoint, wrap_auxiliary_cell=.TRUE., continuous=.TRUE.)
      CALL comm%sum(grid_adjoint)
   END SUBROUTINE qs_upf_inverse_xc_grid_evaluate

! **************************************************************************************************
!> \brief Joint inverse XC response including native source-field and interpolation geometry.
!> \param frames ...
!> \param augmentation ...
!> \param matrices ...
!> \param kweights ...
!> \param points ...
!> \param weights ...
!> \param state ...
!> \param fields ...
!> \param fields_response Native source descriptors at fixed grid indices, including spectral/core response.
!> \param grid ...
!> \param cell ...
!> \param point_response ...
!> \param weights_response ...
!> \param atom_response ...
!> \param cell_response Response of the common physical cell and auxiliary PW lattice.
!> \param reduced_k_response ...
!> \param coefficients_response ...
!> \param core_response ...
!> \param core_gradient_response ...
!> \param comm ...
!> \param response ...
!> \param kinetic Same optional valence kinetic kernels as in the unperturbed evaluation.
!> \param core_tau_response Prescribed frozen-core kinetic-density response.
!> \param paw_reconstruction Same optional PAW-only reconstruction caches as in the evaluation.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_xc_grid_response(frames, augmentation, matrices, kweights, points, weights, state, &
                                              fields, fields_response, grid, cell, point_response, weights_response, &
                                              atom_response, cell_response, reduced_k_response, coefficients_response, &
                                              core_response, core_gradient_response, comm, response, kinetic, core_tau_response, &
                                              paw_reconstruction)
      TYPE(qs_upf_inverse_frame_ref_type), INTENT(IN)    :: frames(:)
      TYPE(qs_upf_inverse_augmentation_type), INTENT(IN) :: augmentation
      COMPLEX(dp), INTENT(IN)                            :: matrices(:, :, :, :)
      REAL(dp), INTENT(IN)                               :: kweights(:), points(:, :), weights(:)
      TYPE(qs_upf_inverse_xc_type), INTENT(IN)           :: state
      REAL(dp), INTENT(IN)                               :: fields(:, :, :, :), &
                                                            fields_response(:, :, :, :)
      TYPE(pw_grid_type), INTENT(IN), POINTER            :: grid
      TYPE(cell_type), INTENT(IN), POINTER               :: cell
      REAL(dp), INTENT(IN) :: point_response(:, :), weights_response(:), atom_response(:, :), &
         cell_response(3, 3), reduced_k_response(:, :)
      COMPLEX(dp), INTENT(IN)                            :: coefficients_response(:, :, :)
      REAL(dp), INTENT(IN) :: core_response(:, :), core_gradient_response(:, :, :)
      TYPE(mp_comm_type), INTENT(IN)                     :: comm
      REAL(dp), INTENT(OUT)                              :: response
      TYPE(upf_inverse_kinetic_type), INTENT(IN), &
         OPTIONAL                                        :: kinetic(:)
      REAL(dp), INTENT(IN), OPTIONAL                     :: core_tau_response(:, :)
      TYPE(upf_inverse_kinetic_type), INTENT(IN), &
         OPTIONAL                                        :: paw_reconstruction(:)

      INTEGER                                            :: axis, np, ns, spin, stride
      REAL(dp)                                           :: lattice_error, lattice_tolerance
      REAL(dp), ALLOCATABLE                              :: background(:, :), gradient(:, :, :), &
                                                            tau(:, :), values(:, :)

      np = SIZE(points, 2)
      ns = SIZE(matrices, 4)
      CPASSERT(state%kinetic .EQV. PRESENT(kinetic))
      CPASSERT(state%kinetic .EQV. PRESENT(core_tau_response))
      stride = 4
      IF (state%kinetic) stride = 5
      CPASSERT(SIZE(fields, 4) == stride*ns)
      ! The auxiliary Fourier fields wrap in every direction, also for isolated electrostatics.
      ! Their interpolation response uses the same lattice as the physical integration points.
      lattice_error = MAXVAL(ABS(grid%dh*SPREAD(REAL(grid%npts, dp), 1, 3) - cell%hmat))
      lattice_tolerance = 128.0_dp*EPSILON(1.0_dp)*MAX(1.0_dp, MAXVAL(ABS(cell%hmat)))
      CPASSERT(lattice_error < lattice_tolerance)
      ALLOCATE (values(np, stride*ns), background(np, ns), gradient(np, 3, ns))
      IF (state%kinetic) ALLOCATE (tau(np, ns))
      CALL native_grid_batch_response(fields, fields_response, grid, cell, points, &
                                      point_response, cell_response, values, wrap_auxiliary_cell=.TRUE., continuous=.TRUE.)
      DO spin = 1, ns
         background(:, spin) = values(:, stride*(spin - 1) + 1)
         DO axis = 1, 3
            gradient(:, axis, spin) = values(:, stride*(spin - 1) + 1 + axis)
         END DO
         IF (state%kinetic) tau(:, spin) = values(:, stride*spin)
      END DO
      CALL qs_upf_inverse_xc_response(frames, augmentation, matrices, kweights, points, weights, state, &
                                      point_response, weights_response, atom_response, cell_response, &
                                      reduced_k_response, coefficients_response, background, gradient, &
                                      core_response, core_gradient_response, comm, response, kinetic, tau, core_tau_response, &
                                      paw_reconstruction)
   END SUBROUTINE qs_upf_inverse_xc_grid_response

! **************************************************************************************************
!> \brief Replicate the native density and its plane-wave Cartesian gradient for mixed XC.
!> \param density Native spin densities, including any prescribed frozen core.
!> \param fields Replicated descriptors, ordered as density,gx,gy,gz and optional tau per spin.
!>               The gradient is the native Fourier derivative of the supplied grid density.
!> \param kinetic Optional matching native kinetic densities, including prescribed core.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_xc_pw_fields(density, fields, kinetic)
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: density(:)
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: fields(:, :, :, :)
      TYPE(pw_r3d_rs_type), INTENT(IN), OPTIONAL         :: kinetic(:)

      INTEGER                                            :: axis, bounds(2, 3), orders(3), spin, &
                                                            stride
      TYPE(pw_c1d_gs_type)                               :: density_g, derivative_g
      TYPE(pw_grid_type), POINTER                        :: grid
      TYPE(pw_r3d_rs_type)                               :: work_r

      CPASSERT(SIZE(density) > 0)
      grid => density(1)%pw_grid
      CPASSERT(ASSOCIATED(grid))
      stride = 4
      IF (PRESENT(kinetic)) THEN
         CPASSERT(SIZE(kinetic) == SIZE(density))
         stride = 5
      END IF
      bounds = grid%bounds_local - SPREAD(grid%bounds(1, :) - 1, 1, 2)
      ALLOCATE (fields(grid%npts(1), grid%npts(2), grid%npts(3), stride*SIZE(density)))
      fields = 0.0_dp
      CALL density_g%create(grid)
      CALL derivative_g%create(grid)
      CALL work_r%create(grid)
      DO spin = 1, SIZE(density)
         CPASSERT(ASSOCIATED(density(spin)%pw_grid, grid))
         fields(bounds(1, 1):bounds(2, 1), bounds(1, 2):bounds(2, 2), &
                bounds(1, 3):bounds(2, 3), stride*(spin - 1) + 1) = density(spin)%array
         IF (PRESENT(kinetic)) THEN
            CPASSERT(ASSOCIATED(kinetic(spin)%pw_grid, grid))
            fields(bounds(1, 1):bounds(2, 1), bounds(1, 2):bounds(2, 2), &
                   bounds(1, 3):bounds(2, 3), stride*spin) = kinetic(spin)%array
         END IF
         CALL pw_transfer(density(spin), density_g)
         DO axis = 1, 3
            derivative_g%array(:) = density_g%array
            orders = 0
            orders(axis) = 1
            CALL pw_derive(derivative_g, orders)
            CALL pw_transfer(derivative_g, work_r)
            fields(bounds(1, 1):bounds(2, 1), bounds(1, 2):bounds(2, 2), &
                   bounds(1, 3):bounds(2, 3), stride*(spin - 1) + 1 + axis) = work_r%array
         END DO
      END DO
      IF (grid%para%mode /= PW_MODE_LOCAL) CALL grid%para%group%sum(fields)
      CALL density_g%release()
      CALL derivative_g%release()
      CALL work_r%release()
   END SUBROUTINE qs_upf_inverse_xc_pw_fields

! **************************************************************************************************
!> \brief Fold weighted density/gradient covectors into the native unweighted scalar potential.
!> \param fields Globally assembled Euclidean covectors from the atomic interpolation transpose.
!> \param potential Created native spin potentials, with the same grid as the source density.
!>                  Contains (lambda_n - div(lambda_g))/dvol on each owned native row.
!> \param kinetic_potential Optional created native kinetic potentials, receiving lambda_tau/dvol.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_xc_pw_adjoint(fields, potential, kinetic_potential)
      REAL(dp), INTENT(IN)                               :: fields(:, :, :, :)
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: potential(:)
      TYPE(pw_r3d_rs_type), INTENT(INOUT), OPTIONAL      :: kinetic_potential(:)

      INTEGER                                            :: axis, bounds(2, 3), orders(3), spin, &
                                                            stride
      TYPE(pw_c1d_gs_type)                               :: derivative_g
      TYPE(pw_grid_type), POINTER                        :: grid
      TYPE(pw_r3d_rs_type)                               :: derivative_r, work_r

      CPASSERT(SIZE(potential) > 0)
      grid => potential(1)%pw_grid
      CPASSERT(ASSOCIATED(grid))
      stride = 4
      IF (PRESENT(kinetic_potential)) THEN
         CPASSERT(SIZE(kinetic_potential) == SIZE(potential))
         stride = 5
      END IF
      CPASSERT(ALL(SHAPE(fields) == [grid%npts, stride*SIZE(potential)]))
      CPASSERT(grid%dvol > 0.0_dp)
      bounds = grid%bounds_local - SPREAD(grid%bounds(1, :) - 1, 1, 2)
      CALL derivative_g%create(grid)
      CALL work_r%create(grid)
      CALL derivative_r%create(grid)
      DO spin = 1, SIZE(potential)
         CPASSERT(ASSOCIATED(potential(spin)%pw_grid, grid))
         potential(spin)%array(:, :, :) = fields(bounds(1, 1):bounds(2, 1), bounds(1, 2):bounds(2, 2), &
                                                 bounds(1, 3):bounds(2, 3), stride*(spin - 1) + 1)
         IF (PRESENT(kinetic_potential)) THEN
            CPASSERT(ASSOCIATED(kinetic_potential(spin)%pw_grid, grid))
            kinetic_potential(spin)%array(:, :, :) = fields(bounds(1, 1):bounds(2, 1), bounds(1, 2):bounds(2, 2), &
                                                            bounds(1, 3):bounds(2, 3), stride*spin)/grid%dvol
         END IF
         DO axis = 1, 3
            work_r%array(:, :, :) = fields(bounds(1, 1):bounds(2, 1), bounds(1, 2):bounds(2, 2), &
                                           bounds(1, 3):bounds(2, 3), stride*(spin - 1) + 1 + axis)
            CALL pw_transfer(work_r, derivative_g)
            orders = 0
            orders(axis) = 1
            CALL pw_derive(derivative_g, orders)
            CALL pw_transfer(derivative_g, derivative_r)
            potential(spin)%array(:, :, :) = potential(spin)%array - derivative_r%array
         END DO
         potential(spin)%array(:, :, :) = potential(spin)%array/grid%dvol
      END DO
      CALL derivative_g%release()
      CALL work_r%release()
      CALL derivative_r%release()
   END SUBROUTINE qs_upf_inverse_xc_pw_adjoint

! **************************************************************************************************
!> \brief Response of native source descriptors at fixed grid topology and reciprocal indices.
!> \param density_response Spin-density response at fixed native real-grid indices, including core.
!> \param fields Current descriptors from qs_upf_inverse_xc_pw_fields.
!> \param lattice_response Response of the PW lattice, dh(:,axis)*npts(axis).
!> \param fields_response Replicated density and gradient response, including delta(G)=-strain^T*G.
!> \param kinetic_response Optional complete kinetic-density response at fixed native grid indices.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_xc_pw_response(density_response, fields, lattice_response, fields_response, kinetic_response)
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: density_response(:)
      REAL(dp), INTENT(IN)                               :: fields(:, :, :, :), &
                                                            lattice_response(3, 3)
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: fields_response(:, :, :, :)
      TYPE(pw_r3d_rs_type), INTENT(IN), OPTIONAL         :: kinetic_response(:)

      INTEGER                                            :: axis, other, spin, stride
      REAL(dp)                                           :: lattice_inverse(3, 3), strain(3, 3)
      TYPE(pw_grid_type), POINTER                        :: grid

      CALL qs_upf_inverse_xc_pw_fields(density_response, fields_response, kinetic_response)
      stride = 4
      IF (PRESENT(kinetic_response)) stride = 5
      CPASSERT(ALL(SHAPE(fields) == SHAPE(fields_response)))
      grid => density_response(1)%pw_grid
      lattice_inverse = grid%dh_inv/SPREAD(REAL(grid%npts, dp), 2, 3)
      strain = MATMUL(lattice_response, lattice_inverse)
      DO spin = 1, SIZE(density_response)
         DO axis = 1, 3
            DO other = 1, 3
               fields_response(:, :, :, stride*(spin - 1) + 1 + axis) = &
                  fields_response(:, :, :, stride*(spin - 1) + 1 + axis) - &
                  strain(other, axis)*fields(:, :, :, stride*(spin - 1) + 1 + other)
            END DO
         END DO
      END DO
   END SUBROUTINE qs_upf_inverse_xc_pw_response

END MODULE qs_upf_inverse_xc
