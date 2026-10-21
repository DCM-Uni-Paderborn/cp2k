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
!>        Keep complete scalar descriptors for nonlocal XC, but evaluate orbitals in bounded blocks.
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
   USE qs_upf_inverse,                  ONLY: qs_upf_inverse_block_size,&
                                              qs_upf_inverse_evaluate
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
   USE skala_density_constraints,       ONLY: skala_atom_domain,&
                                              skala_background_domain,&
                                              skala_density_constraints_count,&
                                              skala_density_constraints_matches,&
                                              skala_density_constraints_type,&
                                              skala_density_samples_active,&
                                              skala_density_samples_only,&
                                              skala_density_samples_put
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
   USE xc_derivatives,                  ONLY: xc_functionals_get_needs,&
                                              xc_functionals_select_floor_limits
   USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_type
   USE xc_rho_set_types,                ONLY: xc_rho_set_create,&
                                              xc_rho_set_release,&
                                              xc_rho_set_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   INTEGER(KIND=int_8), PARAMETER :: max_xc_field_bytes = 64_int_8*1024*1024

   TYPE qs_upf_inverse_frame_ref_type
      TYPE(qs_upf_inverse_type), POINTER :: frame => NULL()
   END TYPE qs_upf_inverse_frame_ref_type

   PUBLIC :: qs_upf_inverse_frame_ref_type, qs_upf_inverse_xc_type
   PUBLIC :: qs_upf_inverse_xc_evaluate, qs_upf_inverse_xc_response, qs_upf_inverse_xc_model_response
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
!> \param density_constraints Selected exact and background atomic descriptor derivatives for this geometry.
!> \param just_energy Evaluate the original energy without AO or background derivative assembly.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_xc_evaluate(frames, augmentation, matrices, kweights, points, weights, &
                                         background_density, background_gradient, core_density, core_gradient, &
                                         xc_section, lsd, comm, state, energy, kinetic, background_tau, core_tau, atom_grid, &
                                         paw_reconstruction, density_constraints, just_energy)
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
      TYPE(skala_density_constraints_type), OPTIONAL, &
         POINTER                                         :: density_constraints
      LOGICAL, INTENT(IN), OPTIONAL                      :: just_energy

      COMPLEX(dp), ALLOCATABLE                           :: cached_gradient(:, :, :, :), &
                                                            cached_soft(:, :, :)
      INTEGER                                            :: nao, nk, np, nspin
      LOGICAL                                            :: energy_only, geometry_fields, &
                                                            store_fields
      REAL(dp), ALLOCATABLE                              :: ebackground(:), integration_weights(:)

      energy_only = .FALSE.
      IF (PRESENT(just_energy)) energy_only = just_energy
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
      geometry_fields = .FALSE.
      IF (PRESENT(atom_grid)) THEN
         geometry_fields = ALLOCATED(atom_grid%inverse_values)
         CPASSERT(geometry_fields .EQV. ALLOCATED(atom_grid%inverse_gradients))
         IF (geometry_fields) THEN
            CPASSERT(ALL(SHAPE(atom_grid%inverse_values) == [np, nao, nk]))
            CPASSERT(ALL(SHAPE(atom_grid%inverse_gradients) == [np, nao, 3, nk]))
            CPASSERT(ALL(SHAPE(atom_grid%points) == SHAPE(points)))
            CPASSERT(ALL(atom_grid%points == points))
         END IF
      END IF
      ! Reuse only the orbital fields shared by the density and operator in this call.
      ! Larger point sets retain the original bounded-block evaluation.
      store_fields = .NOT. energy_only .AND. .NOT. geometry_fields .AND. np > 0 .AND. &
                     64_int_8*INT(np, int_8)*nao*nk <= max_xc_field_bytes
      IF (PRESENT(density_constraints)) THEN
         IF (ASSOCIATED(density_constraints)) THEN
            IF (skala_density_samples_only(density_constraints)) store_fields = .FALSE.
         END IF
      END IF
      IF (store_fields) ALLOCATE (cached_soft(np, nao, nk), cached_gradient(np, nao, 3, nk))
      ALLOCATE (state%operators(nao, nao, nk, nspin), state%density(np, nspin), &
                state%density_gradient(np, 3, nspin), state%energy_density(np), &
                state%vdensity(np, nspin), state%vgradient(np, 3, nspin), &
                state%background_vdensity(np, nspin), state%background_vgradient(np, 3, nspin), &
                state%tau(np, nspin), state%vtau(np, nspin), state%background_vtau(np, nspin), &
                ebackground(np))
      state%density(:, :) = core_density
      state%density_gradient(:, :, :) = core_gradient
      state%tau = 0.0_dp
      IF (state%kinetic) state%tau(:, :) = core_tau
      CALL collocate_xc_fields(frames, augmentation, matrices, kweights, points, state%density, &
                               state%density_gradient, state%tau, kinetic, paw_reconstruction, atom_grid, &
                               cached_soft, cached_gradient)
      IF (PRESENT(density_constraints)) THEN
         IF (ASSOCIATED(density_constraints)) THEN
            IF (skala_density_samples_only(density_constraints)) THEN
               energy = 0.0_dp
               state%operators = CMPLX(0.0_dp, 0.0_dp, dp)
               state%energy_density = 0.0_dp
               state%vdensity = 0.0_dp
               state%vgradient = 0.0_dp
               state%vtau = 0.0_dp
               state%background_vdensity = 0.0_dp
               state%background_vgradient = 0.0_dp
               state%background_vtau = 0.0_dp
               state%weighted_adjoint = xc_section_uses_native_skala_evaluator(xc_section)
               IF (state%weighted_adjoint) THEN
                  ALLOCATE (state%point_derivative(3, np), state%weight_derivative(np), &
                            state%atom_derivative(3, SIZE(frames(1)%frame%positions, 2)), SOURCE=0.0_dp)
               END IF
               RETURN
            END IF
         END IF
      END IF
      ALLOCATE (integration_weights(np))
      integration_weights(:) = weights
      IF (xc_section_uses_native_skala_evaluator(xc_section)) THEN
         CPASSERT(PRESENT(atom_grid) .AND. state%kinetic)
         CALL evaluate_skala_residual(xc_section, comm, state, background_density, background_gradient, &
                                      background_tau, points, weights, atom_grid, frames(1)%frame%positions, energy, &
                                      density_constraints, energy_only)
         ! Skala supplies derivatives of the complete weighted energy. Do not divide by
         ! quadrature weights, including at zero-weight rows used by its nonlocal descriptors.
         integration_weights = 1.0_dp
      ELSE
         CALL evaluate_xc_fields(xc_section, lsd, skala_atom_domain, state%density, state%density_gradient, &
                                 state%energy_density, state%vdensity, state%vgradient, &
                                 kinetic_density=state%tau, vkinetic=state%vtau, kinetic_available=state%kinetic, &
                                 density_constraints=density_constraints, atom_grid=atom_grid)
         CALL evaluate_xc_fields(xc_section, lsd, skala_background_domain, background_density, background_gradient, &
                                 ebackground, state%background_vdensity, state%background_vgradient, &
                                 kinetic_density=background_tau, vkinetic=state%background_vtau, &
                                 density_constraints=density_constraints, atom_grid=atom_grid)
         state%energy_density(:) = state%energy_density - ebackground
         state%background_vdensity = -state%background_vdensity
         state%background_vgradient = -state%background_vgradient
         state%background_vtau = -state%background_vtau
         energy = SUM(weights*state%energy_density)
         CALL comm%sum(energy)
      END IF
      IF (energy_only) THEN
         ! Native energy-only callers leave rho_changed set. The next KS or force
         ! evaluation must rebuild these operators from the complete functional.
         state%operators = CMPLX(0.0_dp, 0.0_dp, dp)
         RETURN
      END IF
      CALL integrate_xc_fields(frames, augmentation, points, integration_weights, state%vdensity, &
                               state%vgradient, state%vtau, comm, state%operators, kinetic, paw_reconstruction, &
                               atom_grid, cached_soft, cached_gradient)
   END SUBROUTINE qs_upf_inverse_xc_evaluate

! **************************************************************************************************
!> \brief Apply the complete atomic Skala XC kernel at fixed geometry and feature layout.
!>        Differentiate both exact and background functionals, retaining their signed adjoints.
!>        Frozen-core fields in state have zero occupation-matrix response. This is a branchwise
!>        electronic response, not a geometry derivative or a selected nonsmooth limiting covector.
!> \param frames ...
!> \param augmentation ...
!> \param matrices_direction Unweighted Hermitian occupation-matrix direction.
!> \param kweights Mesh weights applied once in the field map, not in the returned AO operators.
!> \param points ...
!> \param weights ...
!> \param state Exact primitive fields from the original unperturbed functional.
!> \param background_density ...
!> \param background_gradient ...
!> \param background_tau ...
!> \param background_density_direction ...
!> \param background_gradient_direction ...
!> \param background_tau_direction ...
!> \param xc_section ...
!> \param comm ...
!> \param response Primitive-field directions, weighted covector responses and direct AO kernel.
!>        The negative background covectors still require the caller's interpolation transpose.
!> \param kinetic ...
!> \param atom_grid ...
!> \param paw_reconstruction ...
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_xc_model_response(frames, augmentation, matrices_direction, kweights, points, weights, &
                                               state, background_density, background_gradient, background_tau, &
                                               background_density_direction, background_gradient_direction, &
                                               background_tau_direction, xc_section, comm, response, kinetic, &
                                               atom_grid, paw_reconstruction)
      TYPE(qs_upf_inverse_frame_ref_type), INTENT(IN)    :: frames(:)
      TYPE(qs_upf_inverse_augmentation_type), INTENT(IN) :: augmentation
      COMPLEX(dp), INTENT(IN)                            :: matrices_direction(:, :, :, :)
      REAL(dp), INTENT(IN)                               :: kweights(:), points(:, :), weights(:)
      TYPE(qs_upf_inverse_xc_type), INTENT(IN)           :: state
      REAL(dp), INTENT(IN) :: background_density(:, :), background_gradient(:, :, :), &
         background_tau(:, :), background_density_direction(:, :), &
         background_gradient_direction(:, :, :), background_tau_direction(:, :)
      TYPE(section_vals_type), POINTER                   :: xc_section
      TYPE(mp_comm_type), INTENT(IN)                     :: comm
      TYPE(qs_upf_inverse_xc_type), INTENT(OUT)          :: response
      TYPE(upf_inverse_kinetic_type), INTENT(IN)         :: kinetic(:)
      TYPE(qs_upf_inverse_atom_grid_type), INTENT(IN)    :: atom_grid
      TYPE(upf_inverse_kinetic_type), INTENT(IN), &
         OPTIONAL                                        :: paw_reconstruction(:)

      COMPLEX(dp), ALLOCATABLE                           :: cached_gradient(:, :, :, :), &
                                                            cached_soft(:, :, :)
      INTEGER                                            :: nao, nk, np, ns
      LOGICAL                                            :: geometry_fields
      REAL(dp)                                           :: energy
      REAL(dp), ALLOCATABLE                              :: integration_weights(:)
      TYPE(qs_upf_inverse_xc_type)                       :: fields

      CPASSERT(xc_section_uses_native_skala_evaluator(xc_section))
      CPASSERT(state%kinetic)
      CPASSERT(ALLOCATED(state%density))
      CPASSERT(ALLOCATED(state%density_gradient) .AND. ALLOCATED(state%tau))
      nk = SIZE(frames)
      CPASSERT(nk > 0 .AND. SIZE(kweights) == nk)
      CPASSERT(ALL(kweights >= 0.0_dp))
      CPASSERT(ABS(SUM(kweights) - 1.0_dp) <= 128.0_dp*EPSILON(1.0_dp)*nk)
      CPASSERT(ASSOCIATED(frames(1)%frame))
      nao = SIZE(frames(1)%frame%coefficients, 2)
      np = SIZE(points, 2)
      ns = SIZE(state%density, 2)
      CPASSERT(ns == 1 .OR. ns == 2)
      CPASSERT(ALL(SHAPE(matrices_direction) == [nao, nao, nk, ns]))
      CPASSERT(SIZE(points, 1) == 3 .AND. SIZE(weights) == np)
      CPASSERT(ALL(SHAPE(state%density) == [np, ns]))
      CPASSERT(ALL(SHAPE(state%density_gradient) == [np, 3, ns]))
      CPASSERT(ALL(SHAPE(state%tau) == [np, ns]))
      CPASSERT(ALL(SHAPE(background_density) == [np, ns]))
      CPASSERT(ALL(SHAPE(background_gradient) == [np, 3, ns]))
      CPASSERT(ALL(SHAPE(background_tau) == [np, ns]))
      CPASSERT(ALL(SHAPE(background_density_direction) == [np, ns]))
      CPASSERT(ALL(SHAPE(background_gradient_direction) == [np, 3, ns]))
      CPASSERT(ALL(SHAPE(background_tau_direction) == [np, ns]))
      CPASSERT(state%paw_reconstruction .EQV. PRESENT(paw_reconstruction))
      CPASSERT(ALL(SHAPE(atom_grid%points) == SHAPE(points)))
      CPASSERT(ALL(atom_grid%points == points))
      geometry_fields = ALLOCATED(atom_grid%inverse_values)
      CPASSERT(geometry_fields .EQV. ALLOCATED(atom_grid%inverse_gradients))
      IF (geometry_fields) THEN
         CPASSERT(ALL(SHAPE(atom_grid%inverse_values) == [np, nao, nk]))
         CPASSERT(ALL(SHAPE(atom_grid%inverse_gradients) == [np, nao, 3, nk]))
      ELSE IF (np > 0 .AND. 64_int_8*INT(np, int_8)*nao*nk <= max_xc_field_bytes) THEN
         ALLOCATE (cached_soft(np, nao, nk), cached_gradient(np, nao, 3, nk))
      END IF
      response%kinetic = .TRUE.
      response%paw_reconstruction = PRESENT(paw_reconstruction)
      response%weighted_adjoint = .TRUE.
      ALLOCATE (response%density(np, ns), response%density_gradient(np, 3, ns), response%tau(np, ns), &
                response%vdensity(np, ns), response%vgradient(np, 3, ns), response%vtau(np, ns), &
                response%background_vdensity(np, ns), response%background_vgradient(np, 3, ns), &
                response%background_vtau(np, ns), response%operators(nao, nao, nk, ns))
      response%density = 0.0_dp
      response%density_gradient = 0.0_dp
      response%tau = 0.0_dp
      CALL collocate_xc_fields(frames, augmentation, matrices_direction, kweights, points, response%density, &
                               response%density_gradient, response%tau, kinetic, paw_reconstruction, atom_grid, &
                               cached_soft, cached_gradient)
      fields%density = state%density
      fields%density_gradient = state%density_gradient
      fields%tau = state%tau
      ALLOCATE (fields%energy_density(np), fields%vdensity(np, ns), fields%vgradient(np, 3, ns), fields%vtau(np, ns), &
                fields%background_vdensity(np, ns), fields%background_vgradient(np, 3, ns), fields%background_vtau(np, ns))
      CALL evaluate_skala_residual(xc_section, comm, fields, background_density, background_gradient, &
                                   background_tau, points, weights, atom_grid, frames(1)%frame%positions, energy, &
                                   just_energy=.FALSE., field_response=response, &
                                   background_density_direction=background_density_direction, &
                                   background_gradient_direction=background_gradient_direction, &
                                   background_tau_direction=background_tau_direction)
      ALLOCATE (integration_weights(np), SOURCE=1.0_dp)
      CALL integrate_xc_fields(frames, augmentation, points, integration_weights, response%vdensity, &
                               response%vgradient, response%vtau, comm, response%operators, kinetic, paw_reconstruction, &
                               atom_grid, cached_soft, cached_gradient)
   END SUBROUTINE qs_upf_inverse_xc_model_response

! **************************************************************************************************
!> \brief Add the linear occupation-matrix contribution to the exact atomic primitive fields.
!>        The caller seeds frozen-core fields before this sum. Cache only fixed orbital fields.
!> \param frames ...
!> \param augmentation ...
!> \param matrices Unweighted Hermitian matrices, including signed response directions.
!> \param kweights Mesh weights applied once to the summed fields.
!> \param points ...
!> \param density ...
!> \param density_gradient ...
!> \param kinetic_density ...
!> \param kinetic ...
!> \param paw_reconstruction ...
!> \param atom_grid ...
!> \param cached_soft Optional storage shared with the following operator transpose.
!> \param cached_gradient ...
! **************************************************************************************************
   SUBROUTINE collocate_xc_fields(frames, augmentation, matrices, kweights, points, density, density_gradient, &
                                  kinetic_density, kinetic, paw_reconstruction, atom_grid, cached_soft, cached_gradient)
      TYPE(qs_upf_inverse_frame_ref_type), INTENT(IN)    :: frames(:)
      TYPE(qs_upf_inverse_augmentation_type), INTENT(IN) :: augmentation
      COMPLEX(dp), INTENT(IN)                            :: matrices(:, :, :, :)
      REAL(dp), INTENT(IN)                               :: kweights(:), points(:, :)
      REAL(dp), INTENT(INOUT)                            :: density(:, :), &
                                                            density_gradient(:, :, :), &
                                                            kinetic_density(:, :)
      TYPE(upf_inverse_kinetic_type), INTENT(IN), &
         OPTIONAL                                        :: kinetic(:), paw_reconstruction(:)
      TYPE(qs_upf_inverse_atom_grid_type), INTENT(IN), &
         OPTIONAL                                        :: atom_grid
      COMPLEX(dp), INTENT(OUT), OPTIONAL                 :: cached_soft(:, :, :), &
                                                            cached_gradient(:, :, :, :)

      COMPLEX(dp), ALLOCATABLE                           :: orbital_gradient(:, :, :), &
                                                            projected(:, :), soft(:, :)
      INTEGER                                            :: block, hi, ik, lo, m, nao, nk, np, &
                                                            nproj, nspin, spin
      LOGICAL                                            :: geometry_fields
      REAL(dp), ALLOCATABLE :: gradient(:, :), reconstruction_density(:), &
         reconstruction_gradient(:, :), rho(:), tau(:), tau_augmentation(:)
      TYPE(qs_upf_inverse_type), POINTER                 :: frame

      nk = SIZE(frames)
      nao = SIZE(matrices, 1)
      np = SIZE(points, 2)
      nspin = SIZE(matrices, 4)
      block = qs_upf_inverse_block_size(frames(1)%frame)
      geometry_fields = .FALSE.
      IF (PRESENT(atom_grid)) geometry_fields = ALLOCATED(atom_grid%inverse_values)
      CPASSERT(PRESENT(cached_soft) .EQV. PRESENT(cached_gradient))
      ALLOCATE (soft(block, nao), orbital_gradient(block, nao, 3), rho(np), gradient(np, 3), &
                tau(np), tau_augmentation(np))
      IF (PRESENT(paw_reconstruction)) ALLOCATE (reconstruction_density(np), reconstruction_gradient(np, 3))
      DO ik = 1, nk
         CPASSERT(ASSOCIATED(frames(ik)%frame))
         frame => frames(ik)%frame
         CPASSERT(SIZE(frame%coefficients, 2) == nao)
         DO lo = 1, np, block
            hi = MIN(np, lo + block - 1)
            m = hi - lo + 1
            IF (geometry_fields) THEN
               soft(:m, :) = atom_grid%inverse_values(lo:hi, :, ik)
               orbital_gradient(:m, :, :) = atom_grid%inverse_gradients(lo:hi, :, :, ik)
            ELSE
               CALL qs_upf_inverse_evaluate(frame, points(:, lo:hi), soft(:m, :), orbital_gradient(:m, :, :))
            END IF
            IF (PRESENT(cached_soft)) THEN
               cached_soft(lo:hi, :, ik) = soft(:m, :)
               cached_gradient(lo:hi, :, :, ik) = orbital_gradient(:m, :, :)
            END IF
            DO spin = 1, nspin
               CALL upf_inverse_smooth_density(soft(:m, :), orbital_gradient(:m, :, :), &
                                               matrices(:, :, ik, spin), rho(lo:hi), gradient(lo:hi, :), tau(lo:hi))
               density(lo:hi, spin) = density(lo:hi, spin) + kweights(ik)*rho(lo:hi)
               density_gradient(lo:hi, :, spin) = density_gradient(lo:hi, :, spin) + &
                                                  kweights(ik)*gradient(lo:hi, :)
               IF (PRESENT(kinetic)) kinetic_density(lo:hi, spin) = kinetic_density(lo:hi, spin) + kweights(ik)*tau(lo:hi)
            END DO
         END DO
         nproj = SIZE(frame%coefficients, 1)
         ALLOCATE (projected(nproj, nproj))
         DO spin = 1, nspin
            rho = 0.0_dp
            gradient = 0.0_dp
            CALL upf_inverse_projector_density(frame%coefficients, matrices(:, :, ik, spin), projected)
            CALL qs_upf_inverse_augmentation_collocate(frame, augmentation, projected, rho, gradient, points)
            IF (PRESENT(paw_reconstruction)) THEN
               CALL qs_upf_inverse_paw_collocate(frame, paw_reconstruction, augmentation, points, matrices(:, :, ik, spin), &
                                                 reconstruction_density, reconstruction_gradient)
               rho(:) = rho + reconstruction_density
               gradient(:, :) = gradient + reconstruction_gradient
            END IF
            density(:, spin) = density(:, spin) + kweights(ik)*rho
            density_gradient(:, :, spin) = density_gradient(:, :, spin) + kweights(ik)*gradient
            IF (PRESENT(kinetic)) THEN
               CALL qs_upf_inverse_kinetic_collocate(frame, kinetic, points, matrices(:, :, ik, spin), tau_augmentation)
               IF (PRESENT(paw_reconstruction)) THEN
                  CALL qs_upf_inverse_kinetic_collocate(frame, paw_reconstruction, points, matrices(:, :, ik, spin), &
                                                        reconstruction_density)
                  tau_augmentation(:) = tau_augmentation + reconstruction_density
               END IF
               kinetic_density(:, spin) = kinetic_density(:, spin) + kweights(ik)*tau_augmentation
            END IF
         END DO
         DEALLOCATE (projected)
      END DO
   END SUBROUTINE collocate_xc_fields

! **************************************************************************************************
!> \brief Hermitian transpose of the same exact atomic primitive-field map.
!>        Input covectors include quadrature weights when required. AO operators have no k weight.
!> \param frames ...
!> \param augmentation ...
!> \param points ...
!> \param integration_weights ...
!> \param vdensity ...
!> \param vgradient ...
!> \param vtau ...
!> \param comm ...
!> \param operators ...
!> \param kinetic ...
!> \param paw_reconstruction ...
!> \param atom_grid ...
!> \param cached_soft ...
!> \param cached_gradient ...
! **************************************************************************************************
   SUBROUTINE integrate_xc_fields(frames, augmentation, points, integration_weights, vdensity, vgradient, vtau, &
                                  comm, operators, kinetic, paw_reconstruction, atom_grid, cached_soft, cached_gradient)
      TYPE(qs_upf_inverse_frame_ref_type), INTENT(IN)    :: frames(:)
      TYPE(qs_upf_inverse_augmentation_type), INTENT(IN) :: augmentation
      REAL(dp), INTENT(IN)                               :: points(:, :), integration_weights(:), &
                                                            vdensity(:, :), vgradient(:, :, :), &
                                                            vtau(:, :)
      TYPE(mp_comm_type), INTENT(IN)                     :: comm
      COMPLEX(dp), INTENT(OUT)                           :: operators(:, :, :, :)
      TYPE(upf_inverse_kinetic_type), INTENT(IN), &
         OPTIONAL                                        :: kinetic(:), paw_reconstruction(:)
      TYPE(qs_upf_inverse_atom_grid_type), INTENT(IN), &
         OPTIONAL                                        :: atom_grid
      COMPLEX(dp), INTENT(IN), OPTIONAL                  :: cached_soft(:, :, :), &
                                                            cached_gradient(:, :, :, :)

      COMPLEX(dp), ALLOCATABLE                           :: kinetic_operator(:, :), &
                                                            orbital_gradient(:, :, :), &
                                                            partial(:, :), qoperator(:, :), &
                                                            soft(:, :)
      INTEGER                                            :: block, hi, ik, lo, m, nao, nk, np, &
                                                            nproj, nspin, spin
      LOGICAL                                            :: geometry_fields
      TYPE(qs_upf_inverse_type), POINTER                 :: frame

      nk = SIZE(frames)
      nao = SIZE(operators, 1)
      np = SIZE(points, 2)
      nspin = SIZE(operators, 4)
      block = qs_upf_inverse_block_size(frames(1)%frame)
      geometry_fields = .FALSE.
      IF (PRESENT(atom_grid)) geometry_fields = ALLOCATED(atom_grid%inverse_values)
      CPASSERT(PRESENT(cached_soft) .EQV. PRESENT(cached_gradient))
      ALLOCATE (soft(block, nao), orbital_gradient(block, nao, 3), partial(nao, nao))
      DO ik = 1, nk
         frame => frames(ik)%frame
         nproj = SIZE(frame%coefficients, 1)
         ALLOCATE (qoperator(nproj, nproj), kinetic_operator(nao, nao))
         operators(:, :, ik, :) = CMPLX(0.0_dp, 0.0_dp, dp)
         DO lo = 1, np, block
            hi = MIN(np, lo + block - 1)
            m = hi - lo + 1
            IF (geometry_fields) THEN
               soft(:m, :) = atom_grid%inverse_values(lo:hi, :, ik)
               orbital_gradient(:m, :, :) = atom_grid%inverse_gradients(lo:hi, :, :, ik)
            ELSE IF (PRESENT(cached_soft)) THEN
               soft(:m, :) = cached_soft(lo:hi, :, ik)
               orbital_gradient(:m, :, :) = cached_gradient(lo:hi, :, :, ik)
            ELSE
               CALL qs_upf_inverse_evaluate(frame, points(:, lo:hi), soft(:m, :), orbital_gradient(:m, :, :))
            END IF
            DO spin = 1, nspin
               CALL upf_inverse_semilocal_pair(soft(:m, :), soft(:m, :), &
                                               orbital_gradient(:m, :, :), orbital_gradient(:m, :, :), &
                                               integration_weights(lo:hi), vdensity(lo:hi, spin), &
                                               vgradient(lo:hi, :, spin), vtau(lo:hi, spin), partial)
               operators(:, :, ik, spin) = operators(:, :, ik, spin) + partial
            END DO
         END DO
         DO spin = 1, nspin
            ! Point owners can have different block counts, including none.
            CALL comm%sum(operators(:, :, ik, spin))
            CALL qs_upf_inverse_augmentation_integrate(frame, augmentation, vdensity(:, spin), &
                                                       vgradient(:, :, spin), qoperator, points, integration_weights, comm)
            operators(:, :, ik, spin) = operators(:, :, ik, spin) + &
                                        MATMUL(CONJG(TRANSPOSE(frame%coefficients)), MATMUL(qoperator, frame%coefficients))
            IF (PRESENT(paw_reconstruction)) THEN
               CALL qs_upf_inverse_paw_operator(frame, paw_reconstruction, augmentation, points, integration_weights, &
                                                vdensity(:, spin), vgradient(:, :, spin), kinetic_operator, comm)
               operators(:, :, ik, spin) = operators(:, :, ik, spin) + kinetic_operator
            END IF
            IF (PRESENT(kinetic)) THEN
               CALL qs_upf_inverse_kinetic_operator(frame, kinetic, points, integration_weights, vtau(:, spin), &
                                                    kinetic_operator, comm)
               operators(:, :, ik, spin) = operators(:, :, ik, spin) + kinetic_operator
               IF (PRESENT(paw_reconstruction)) THEN
                  CALL qs_upf_inverse_kinetic_operator(frame, paw_reconstruction, points, integration_weights, &
                                                       vtau(:, spin), kinetic_operator, comm)
                  operators(:, :, ik, spin) = operators(:, :, ik, spin) + kinetic_operator
               END IF
            END IF
         END DO
         DEALLOCATE (qoperator, kinetic_operator)
      END DO
   END SUBROUTINE integrate_xc_fields

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
!> \param density_constraints ...
!> \param just_energy Evaluate the original energy without AO or background derivative assembly.
!> \param field_response Optional exact primitive direction and resulting signed covector responses.
!> \param background_density_direction ...
!> \param background_gradient_direction ...
!> \param background_tau_direction ...
! **************************************************************************************************
   SUBROUTINE evaluate_skala_residual(xc_section, comm, state, background_density, background_gradient, &
                                      background_tau, points, weights, atom_grid, positions, energy, &
                                      density_constraints, just_energy, field_response, background_density_direction, &
                                      background_gradient_direction, background_tau_direction)
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
      TYPE(skala_density_constraints_type), OPTIONAL, &
         POINTER                                         :: density_constraints
      LOGICAL, INTENT(IN)                                :: just_energy
      TYPE(qs_upf_inverse_xc_type), INTENT(INOUT), &
         OPTIONAL                                        :: field_response
      REAL(dp), INTENT(IN), OPTIONAL :: background_density_direction(:, :), &
         background_gradient_direction(:, :, :), background_tau_direction(:, :)

      INTEGER                                            :: atom, local_atom, np, ns, row, spin
      INTEGER(int_8), ALLOCATABLE                        :: sizes(:)
      INTEGER, ALLOCATABLE                               :: local_atoms(:)
      REAL(dp)                                           :: background_energy, scale
      REAL(dp), ALLOCATABLE :: atom_coords(:, :), atom_derivative(:, :), &
         background_atom_derivative(:, :), background_base_derivative(:), &
         background_point_derivative(:, :), background_weight_derivative(:), base_derivative(:), &
         base_weights(:), coords(:, :), density(:, :), density_direction(:, :), &
         density_response(:, :), gradient(:, :, :), gradient_direction(:, :, :), &
         gradient_response(:, :, :), grid_weights(:), point_derivative(:, :), tau(:, :), &
         tau_direction(:, :), tau_response(:, :), vdensity(:, :), vgradient(:, :, :), vtau(:, :), &
         weight_derivative(:)

      CPASSERT(PRESENT(background_density_direction) .EQV. PRESENT(field_response))
      CPASSERT(PRESENT(background_gradient_direction) .EQV. PRESENT(field_response))
      CPASSERT(PRESENT(background_tau_direction) .EQV. PRESENT(field_response))
      IF (PRESENT(field_response)) THEN
         CPASSERT(.NOT. just_energy)
      END IF
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
      IF (PRESENT(field_response)) THEN
         ALLOCATE (density_direction(np, 2), gradient_direction(np, 3, 2), tau_direction(np, 2))
         DO spin = 1, 2
            density_direction(:, spin) = scale*field_response%density(:, MIN(spin, ns))
            gradient_direction(:, :, spin) = scale*field_response%density_gradient(:, :, MIN(spin, ns))
            tau_direction(:, spin) = scale*field_response%tau(:, MIN(spin, ns))
         END DO
      END IF
      IF (just_energy) THEN
         CALL skala_gapw_atom_composite_energy(xc_section, comm, density, gradient, tau, coords, grid_weights, &
                                               base_weights, sizes, atom_coords, energy)
      ELSE IF (PRESENT(field_response)) THEN
         CALL skala_gapw_atom_composite_energy(xc_section, comm, density, gradient, tau, coords, grid_weights, &
                                               base_weights, sizes, atom_coords, energy, vdensity, vgradient, vtau, &
                                               point_derivative, weight_derivative, base_derivative, atom_derivative, &
                                               density_constraints, skala_atom_domain, local_atoms, &
                                               density_direction, gradient_direction, tau_direction, &
                                               density_response, gradient_response, tau_response)
      ELSE
         CALL skala_gapw_atom_composite_energy(xc_section, comm, density, gradient, tau, coords, grid_weights, &
                                               base_weights, sizes, atom_coords, energy, vdensity, vgradient, vtau, &
                                               point_derivative, weight_derivative, base_derivative, atom_derivative, &
                                               density_constraints, skala_atom_domain, local_atoms)
      END IF
      IF (.NOT. just_energy) THEN
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
      END IF
      IF (PRESENT(field_response)) THEN
         DO spin = 1, ns
            IF (ns == 1) THEN
               field_response%vdensity(:, spin) = 0.5_dp*SUM(density_response, DIM=2)
               field_response%vgradient(:, :, spin) = 0.5_dp*SUM(gradient_response, DIM=3)
               field_response%vtau(:, spin) = 0.5_dp*SUM(tau_response, DIM=2)
            ELSE
               field_response%vdensity(:, spin) = density_response(:, spin)
               field_response%vgradient(:, :, spin) = gradient_response(:, :, spin)
               field_response%vtau(:, spin) = tau_response(:, spin)
            END IF
         END DO
         DO spin = 1, 2
            density_direction(:, spin) = scale*background_density_direction(:, MIN(spin, ns))
            gradient_direction(:, :, spin) = scale*background_gradient_direction(:, :, MIN(spin, ns))
            tau_direction(:, spin) = scale*background_tau_direction(:, MIN(spin, ns))
         END DO
      END IF
      DO spin = 1, 2
         density(:, spin) = scale*background_density(:, MIN(spin, ns))
         gradient(:, :, spin) = scale*background_gradient(:, :, MIN(spin, ns))
         tau(:, spin) = scale*background_tau(:, MIN(spin, ns))
      END DO
      IF (just_energy) THEN
         CALL skala_gapw_atom_composite_energy(xc_section, comm, density, gradient, tau, coords, grid_weights, &
                                               base_weights, sizes, atom_coords, background_energy)
      ELSE IF (PRESENT(field_response)) THEN
         CALL skala_gapw_atom_composite_energy(xc_section, comm, density, gradient, tau, coords, grid_weights, &
                                               base_weights, sizes, atom_coords, background_energy, vdensity, vgradient, vtau, &
                                               background_point_derivative, background_weight_derivative, &
                                               background_base_derivative, background_atom_derivative, &
                                               density_constraints, skala_background_domain, local_atoms, &
                                               density_direction, gradient_direction, tau_direction, &
                                               density_response, gradient_response, tau_response)
      ELSE
         CALL skala_gapw_atom_composite_energy(xc_section, comm, density, gradient, tau, coords, grid_weights, &
                                               base_weights, sizes, atom_coords, background_energy, vdensity, vgradient, vtau, &
                                               background_point_derivative, background_weight_derivative, &
                                               background_base_derivative, background_atom_derivative, &
                                               density_constraints, skala_background_domain, local_atoms)
      END IF
      IF (.NOT. just_energy) THEN
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
      END IF
      IF (PRESENT(field_response)) THEN
         DO spin = 1, ns
            IF (ns == 1) THEN
               field_response%background_vdensity(:, spin) = -0.5_dp*SUM(density_response, DIM=2)
               field_response%background_vgradient(:, :, spin) = -0.5_dp*SUM(gradient_response, DIM=3)
               field_response%background_vtau(:, spin) = -0.5_dp*SUM(tau_response, DIM=2)
            ELSE
               field_response%background_vdensity(:, spin) = -density_response(:, spin)
               field_response%background_vgradient(:, :, spin) = -gradient_response(:, :, spin)
               field_response%background_vtau(:, spin) = -tau_response(:, spin)
            END IF
         END DO
      END IF
      state%weighted_adjoint = .TRUE.
      state%energy_density = 0.0_dp
      energy = energy - background_energy
      ALLOCATE (state%point_derivative(3, np), state%weight_derivative(np), state%atom_derivative(3, SIZE(positions, 2)))
      IF (just_energy) THEN
         state%vdensity = 0.0_dp
         state%vgradient = 0.0_dp
         state%vtau = 0.0_dp
         state%background_vdensity = 0.0_dp
         state%background_vgradient = 0.0_dp
         state%background_vtau = 0.0_dp
         state%point_derivative = 0.0_dp
         state%weight_derivative = 0.0_dp
         state%atom_derivative = 0.0_dp
         RETURN
      END IF
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
      INTEGER                                            :: block, hi, ik, lo, m, nao, np, nproj, &
                                                            spin
      REAL(dp)                                           :: kinetic_response, local_response
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
      block = qs_upf_inverse_block_size(frames(1)%frame)
      ALLOCATE (soft(block, nao), gradient(block, nao, 3), dsoft(block, nao), dgradient(block, nao, 3), &
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
         local_response = 0.0_dp
         DO lo = 1, np, block
            hi = MIN(np, lo + block - 1)
            m = hi - lo + 1
            CALL qs_upf_inverse_evaluate(frame, points(:, lo:hi), soft(:m, :), gradient(:m, :, :), &
                                         dsoft(:m, :), dgradient(:m, :, :), point_response(:, lo:hi), atom_response, &
                                         cell_response, reduced_k_response(:, ik), coefficients_response(:, :, ik))
            DO spin = 1, SIZE(matrices, 4)
               CALL upf_inverse_semilocal_response(soft(:m, :), soft(:m, :), &
                                                   gradient(:m, :, :), gradient(:m, :, :), integration_weights(lo:hi), &
                                                   state%vdensity(lo:hi, spin), state%vgradient(lo:hi, :, spin), &
                                                   state%vtau(lo:hi, spin), dsoft(:m, :), dsoft(:m, :), &
                                                   dgradient(:m, :, :), dgradient(:m, :, :), &
                                                   zero(lo:hi), zero(lo:hi), zeros(lo:hi, :), zero(lo:hi), matrix)
               local_response = local_response + kweights(ik)*REAL(SUM(matrices(:, :, ik, spin)*TRANSPOSE(matrix)), dp)
            END DO
         END DO
         CALL comm%sum(local_response)
         response = response + local_response
         DO spin = 1, SIZE(matrices, 4)
            CALL qs_upf_inverse_augmentation_integrate(frame, augmentation, state%vdensity(:, spin), &
                                                       state%vgradient(:, :, spin), qoperator, points, integration_weights, comm)
            CALL qs_upf_inverse_augmentation_response(frame, augmentation, state%vdensity(:, spin), &
                                     state%vgradient(:, :, spin), point_response, atom_response, cell_response, zero, zero, zeros, &
                                                      dqoperator, points, integration_weights, comm)
            matrix(:, :) = &
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
!> \param domain Exact atomic or interpolated background evaluation domain.
!> \param density ...
!> \param gradient ...
!> \param energy_density ...
!> \param vdensity ...
!> \param vgradient ...
!> \param kinetic_density Optional independent kinetic density, with factor one half.
!> \param vkinetic Optional potential conjugate to kinetic density.
!> \param kinetic_available Optional explicit availability flag for allocated zero placeholders.
!> \param density_constraints Owned physical descriptor targets and first-derivative weights.
!> \param atom_grid Original global atom and row mapping for selected local samples.
! **************************************************************************************************
   SUBROUTINE evaluate_xc_fields(xc_section, lsd, domain, density, gradient, energy_density, vdensity, vgradient, &
                                 kinetic_density, vkinetic, kinetic_available, density_constraints, atom_grid)
      TYPE(section_vals_type), POINTER                   :: xc_section
      LOGICAL, INTENT(IN)                                :: lsd
      INTEGER, INTENT(IN)                                :: domain
      REAL(dp), INTENT(IN)                               :: density(:, :), gradient(:, :, :)
      REAL(dp), INTENT(OUT)                              :: energy_density(:), vdensity(:, :), &
                                                            vgradient(:, :, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: kinetic_density(:, :)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: vkinetic(:, :)
      LOGICAL, INTENT(IN), OPTIONAL                      :: kinetic_available
      TYPE(skala_density_constraints_type), OPTIONAL, &
         POINTER                                         :: density_constraints
      TYPE(qs_upf_inverse_atom_grid_type), INTENT(IN), &
         OPTIONAL                                        :: atom_grid

      INTEGER                                            :: bounds(2, 3), np, nspin, spin
      INTEGER, ALLOCATABLE                               :: keys(:, :)
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
      IF (PRESENT(density_constraints)) THEN
         IF (ASSOCIATED(density_constraints)) THEN
            IF (skala_density_constraints_count(density_constraints, domain) > 0) THEN
               CPASSERT(PRESENT(atom_grid))
               CPASSERT(SIZE(atom_grid%atoms) == np .AND. SIZE(atom_grid%atom_rows) == np)
               ALLOCATE (keys(3, np))
               keys(1, :) = atom_grid%atoms
               keys(2, :) = atom_grid%atom_rows
               keys(3, :) = 1
               CALL xc_functionals_select_floor_limits(functional, density_constraints, fields, nspin, domain, keys)
            END IF
         END IF
      END IF
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
!> \param density_constraints Owned selected derivatives, invalidated by geometry or atomic-layout changes.
!> \param just_energy Evaluate the original energy without AO or background derivative assembly.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_xc_grid_evaluate(frames, augmentation, matrices, kweights, points, weights, &
                                              fields, grid, cell, core_density, core_gradient, xc_section, &
                                              lsd, comm, state, energy, grid_adjoint, kinetic, core_tau, atom_grid, &
                                              paw_reconstruction, density_constraints, just_energy)
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
      TYPE(skala_density_constraints_type), OPTIONAL, &
         POINTER                                         :: density_constraints
      LOGICAL, INTENT(IN), OPTIONAL                      :: just_energy

      INTEGER                                            :: atom, axis, np, ns, spin, stride
      INTEGER, ALLOCATABLE                               :: atomic_sizes(:)
      REAL(dp), ALLOCATABLE                              :: background(:, :), covectors(:, :), &
                                                            gradient(:, :, :), &
                                                            integration_weights(:), tau(:, :), &
                                                            values(:, :)
      TYPE(skala_density_constraints_type), POINTER      :: selected_constraints

      np = SIZE(points, 2)
      ns = SIZE(matrices, 4)
      NULLIFY (selected_constraints)
      IF (PRESENT(density_constraints)) THEN
         IF (ASSOCIATED(density_constraints)) THEN
            IF (skala_density_constraints_count(density_constraints, skala_atom_domain) + &
                skala_density_constraints_count(density_constraints, skala_background_domain) > 0) THEN
               CPASSERT(PRESENT(atom_grid))
               ALLOCATE (atomic_sizes(SIZE(frames(1)%frame%positions, 2)))
               DO atom = 1, SIZE(atomic_sizes)
                  atomic_sizes(atom) = COUNT(atom_grid%atoms == atom)
               END DO
               CALL comm%sum(atomic_sizes)
               IF (skala_density_constraints_matches(density_constraints, ns, grid%npts, grid%bounds, &
                                                     cell%hmat, frames(1)%frame%positions, atomic_sizes)) THEN
                  selected_constraints => density_constraints
               END IF
            END IF
         END IF
      END IF
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
      IF (PRESENT(density_constraints)) THEN
         IF (ASSOCIATED(density_constraints)) THEN
            IF (skala_density_samples_only(density_constraints)) selected_constraints => density_constraints
         END IF
      END IF
      CALL qs_upf_inverse_xc_evaluate(frames, augmentation, matrices, kweights, points, weights, &
                                      background, gradient, core_density, core_gradient, xc_section, &
                                      lsd, comm, state, energy, kinetic, tau, core_tau, atom_grid, paw_reconstruction, &
                                      selected_constraints, just_energy)
      IF (PRESENT(density_constraints)) THEN
         IF (ASSOCIATED(density_constraints)) THEN
            IF (skala_density_samples_active(density_constraints)) THEN
               IF (PRESENT(atom_grid)) THEN
                  IF (state%kinetic) THEN
                     CALL record_atom_descriptors(density_constraints, atom_grid, state, background, tau)
                  ELSE
                     CALL record_atom_descriptors(density_constraints, atom_grid, state, background)
                  END IF
               END IF
               IF (skala_density_samples_only(density_constraints)) THEN
                  ! Pure descriptor sampling has zero XC covectors. No transpose is needed.
                  ALLOCATE (grid_adjoint(SIZE(fields, 1), SIZE(fields, 2), SIZE(fields, 3), SIZE(fields, 4)), SOURCE=0.0_dp)
                  RETURN
               END IF
            END IF
         END IF
      END IF
      IF (PRESENT(just_energy)) THEN
         IF (just_energy) THEN
            ALLOCATE (grid_adjoint(SIZE(fields, 1), SIZE(fields, 2), SIZE(fields, 3), SIZE(fields, 4)), SOURCE=0.0_dp)
            RETURN
         END IF
      END IF
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
!> \brief Retain the physical inputs of the exact and background atomic Skala evaluations.
!>        Recording uses the KS owner even when no subgradient has been selected yet.
!>        Core terms, original augmentation and paired-spin scaling are unchanged.
!> \param constraints ...
!> \param atom_grid ...
!> \param state ...
!> \param background_density ...
!> \param background_tau Optional kinetic inputs, present when the functional requests them.
! **************************************************************************************************
   SUBROUTINE record_atom_descriptors(constraints, atom_grid, state, background_density, background_tau)
      TYPE(skala_density_constraints_type), &
         INTENT(INOUT)                                   :: constraints
      TYPE(qs_upf_inverse_atom_grid_type), INTENT(IN)    :: atom_grid
      TYPE(qs_upf_inverse_xc_type), INTENT(IN)           :: state
      REAL(dp), INTENT(IN)                               :: background_density(:, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: background_tau(:, :)

      INTEGER, ALLOCATABLE                               :: indices(:, :)
      REAL(dp)                                           :: scale

      ALLOCATE (indices(3, SIZE(atom_grid%atoms)))
      indices(1, :) = atom_grid%atoms
      indices(2, :) = atom_grid%atom_rows
      indices(3, :) = 1
      scale = 1.0_dp
      IF (SIZE(state%density, 2) == 1) scale = 0.5_dp
      IF (state%kinetic) THEN
         CPASSERT(PRESENT(background_tau))
         CALL skala_density_samples_put(constraints, skala_atom_domain, indices, &
                                        scale*state%density, scale*state%tau)
         CALL skala_density_samples_put(constraints, skala_background_domain, indices, &
                                        scale*background_density, scale*background_tau)
      ELSE
         CALL skala_density_samples_put(constraints, skala_atom_domain, indices, scale*state%density)
         CALL skala_density_samples_put(constraints, skala_background_domain, indices, scale*background_density)
      END IF
   END SUBROUTINE record_atom_descriptors

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
