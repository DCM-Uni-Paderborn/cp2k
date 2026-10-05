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
   USE kinds,                           ONLY: dp
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
   USE qs_upf_inverse_augmentation,     ONLY: qs_upf_inverse_augmentation_collocate,&
                                              qs_upf_inverse_augmentation_integrate,&
                                              qs_upf_inverse_augmentation_response
   USE qs_upf_inverse_types,            ONLY: qs_upf_inverse_augmentation_type,&
                                              qs_upf_inverse_type
   USE upf_inverse_quadrature,          ONLY: upf_inverse_semilocal_pair,&
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

   TYPE qs_upf_inverse_xc_type
      COMPLEX(dp), ALLOCATABLE :: soft(:, :, :), gradient(:, :, :, :), operators(:, :, :, :)
      REAL(dp), ALLOCATABLE :: density(:, :), density_gradient(:, :, :), energy_density(:), &
                               vdensity(:, :), vgradient(:, :, :), background_vdensity(:, :), &
                               background_vgradient(:, :, :)
   END TYPE qs_upf_inverse_xc_type

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
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_xc_evaluate(frames, augmentation, matrices, kweights, points, weights, &
                                         background_density, background_gradient, core_density, core_gradient, &
                                         xc_section, lsd, comm, state, energy)
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

      COMPLEX(dp), ALLOCATABLE                           :: projected(:, :), qoperator(:, :)
      INTEGER                                            :: ik, nao, nk, np, nproj, nspin, spin
      REAL(dp), ALLOCATABLE                              :: ebackground(:), gradient(:, :), rho(:), &
                                                            tau(:), zero(:)
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
      ALLOCATE (state%soft(np, nao, nk), state%gradient(np, nao, 3, nk), &
                state%operators(nao, nao, nk, nspin), state%density(np, nspin), &
                state%density_gradient(np, 3, nspin), state%energy_density(np), &
                state%vdensity(np, nspin), state%vgradient(np, 3, nspin), &
                state%background_vdensity(np, nspin), state%background_vgradient(np, 3, nspin), &
                rho(np), gradient(np, 3), tau(np), zero(np), ebackground(np))
      zero = 0.0_dp
      state%density(:, :) = core_density
      state%density_gradient(:, :, :) = core_gradient
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
            projected(:, :) = MATMUL(frame%coefficients, MATMUL(matrices(:, :, ik, spin), &
                                                                CONJG(TRANSPOSE(frame%coefficients))))
            CALL qs_upf_inverse_augmentation_collocate(frame, augmentation, projected, rho, gradient, points)
            state%density(:, spin) = state%density(:, spin) + kweights(ik)*rho
            state%density_gradient(:, :, spin) = state%density_gradient(:, :, spin) + kweights(ik)*gradient
         END DO
         DEALLOCATE (projected)
      END DO
      CALL evaluate_xc_fields(xc_section, lsd, state%density, state%density_gradient, &
                              state%energy_density, state%vdensity, state%vgradient)
      CALL evaluate_xc_fields(xc_section, lsd, background_density, background_gradient, &
                              ebackground, state%background_vdensity, state%background_vgradient)
      state%energy_density(:) = state%energy_density - ebackground
      state%background_vdensity = -state%background_vdensity
      state%background_vgradient = -state%background_vgradient
      energy = SUM(weights*state%energy_density)
      CALL comm%sum(energy)
      DO ik = 1, nk
         frame => frames(ik)%frame
         nproj = SIZE(frame%coefficients, 1)
         ALLOCATE (qoperator(nproj, nproj))
         DO spin = 1, nspin
            CALL upf_inverse_semilocal_pair(state%soft(:, :, ik), state%soft(:, :, ik), &
                                            state%gradient(:, :, :, ik), state%gradient(:, :, :, ik), weights, &
                                            state%vdensity(:, spin), state%vgradient(:, :, spin), zero, &
                                            state%operators(:, :, ik, spin), comm)
            CALL qs_upf_inverse_augmentation_integrate(frame, augmentation, state%vdensity(:, spin), &
                                                       state%vgradient(:, :, spin), qoperator, points, weights, comm)
            state%operators(:, :, ik, spin) = state%operators(:, :, ik, spin) + &
                                              MATMUL(CONJG(TRANSPOSE(frame%coefficients)), MATMUL(qoperator, frame%coefficients))
         END DO
         DEALLOCATE (qoperator)
      END DO
   END SUBROUTINE qs_upf_inverse_xc_evaluate

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
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_xc_response(frames, augmentation, matrices, kweights, points, weights, state, &
                                         point_response, weights_response, atom_response, cell_response, &
                                         reduced_k_response, coefficients_response, background_response, &
                                         background_gradient_response, core_response, core_gradient_response, comm, response)
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

      COMPLEX(dp), ALLOCATABLE :: dgradient(:, :, :), dqoperator(:, :), dsoft(:, :), &
         gradient(:, :, :), matrix(:, :), qoperator(:, :), soft(:, :)
      INTEGER                                            :: ik, nao, np, nproj, spin
      REAL(dp), ALLOCATABLE                              :: zero(:), zeros(:, :)
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
      ALLOCATE (soft(np, nao), gradient(np, nao, 3), dsoft(np, nao), dgradient(np, nao, 3), &
                matrix(nao, nao), zero(np), zeros(np, 3))
      zero = 0.0_dp
      zeros = 0.0_dp
      response = SUM(weights_response*state%energy_density)
      DO spin = 1, SIZE(matrices, 4)
         response = response + SUM(weights*(state%vdensity(:, spin)*core_response(:, spin) + &
                                            SUM(state%vgradient(:, :, spin)*core_gradient_response(:, :, spin), DIM=2) + &
                                            state%background_vdensity(:, spin)*background_response(:, spin) + &
                                       SUM(state%background_vgradient(:, :, spin)*background_gradient_response(:, :, spin), DIM=2)))
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
                                                state%gradient(:, :, :, ik), state%gradient(:, :, :, ik), weights, &
                                   state%vdensity(:, spin), state%vgradient(:, :, spin), zero, dsoft, dsoft, dgradient, dgradient, &
                                                zero, zero, zeros, zero, matrix, comm)
            CALL qs_upf_inverse_augmentation_integrate(frame, augmentation, state%vdensity(:, spin), &
                                                       state%vgradient(:, :, spin), qoperator, points, weights, comm)
            CALL qs_upf_inverse_augmentation_response(frame, augmentation, state%vdensity(:, spin), &
                                     state%vgradient(:, :, spin), point_response, atom_response, cell_response, zero, zero, zeros, &
                                                      dqoperator, points, weights, comm)
            matrix(:, :) = matrix + &
                           MATMUL(CONJG(TRANSPOSE(coefficients_response(:, :, ik))), MATMUL(qoperator, frame%coefficients)) + &
                           MATMUL(CONJG(TRANSPOSE(frame%coefficients)), MATMUL(dqoperator, frame%coefficients)) + &
                           MATMUL(CONJG(TRANSPOSE(frame%coefficients)), MATMUL(qoperator, coefficients_response(:, :, ik)))
            response = response + kweights(ik)*REAL(SUM(matrices(:, :, ik, spin)*TRANSPOSE(matrix)), dp)
         END DO
         DEALLOCATE (qoperator, dqoperator)
      END DO
   END SUBROUTINE qs_upf_inverse_xc_response

! **************************************************************************************************
!> \brief Native LDA/GGA XC evaluation with unweighted Cartesian density covectors.
!> \param xc_section ...
!> \param lsd ...
!> \param density ...
!> \param gradient ...
!> \param energy_density ...
!> \param vdensity ...
!> \param vgradient ...
! **************************************************************************************************
   SUBROUTINE evaluate_xc_fields(xc_section, lsd, density, gradient, energy_density, vdensity, vgradient)
      TYPE(section_vals_type), POINTER                   :: xc_section
      LOGICAL, INTENT(IN)                                :: lsd
      REAL(dp), INTENT(IN)                               :: density(:, :), gradient(:, :, :)
      REAL(dp), INTENT(OUT)                              :: energy_density(:), vdensity(:, :), &
                                                            vgradient(:, :, :)

      INTEGER                                            :: bounds(2, 3), np, nspin, spin
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
      IF (needs%tau .OR. needs%tau_spin .OR. needs%laplace_rho .OR. needs%laplace_rho_spin) THEN
         CPABORT("Inverse mixed-quadrature XC needs independent kinetic or Laplacian augmentation")
      END IF
      energy_density = 0.0_dp
      vdensity = 0.0_dp
      vgradient = 0.0_dp
      IF (np == 0) RETURN
      bounds(1, :) = 1
      bounds(2, :) = [np, 1, 1]
      ALLOCATE (rho(np, 1, nspin), drho(4, np, 1, nspin), tau(np, 1, nspin), weights(np, 1))
      ALLOCATE (vxc(np, 1, nspin), vxg(3, np, 1, nspin))
      NULLIFY (vtau, derivative)
      rho(:, 1, :) = density
      tau = 0.0_dp
      weights = 1.0_dp
      vxc = 0.0_dp
      vxg = 0.0_dp
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
      DO spin = 1, nspin
         vgradient(:, :, spin) = TRANSPOSE(vxg(:, :, 1, spin))
      END DO
      CALL xc_dset_release(derivatives)
      CALL xc_rho_set_release(fields)
      DEALLOCATE (vxc, vxg)
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
!>               Four descriptors per spin, ordered (density,gx,gy,gz), with no quadrature weights.
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
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_xc_grid_evaluate(frames, augmentation, matrices, kweights, points, weights, &
                                              fields, grid, cell, core_density, core_gradient, xc_section, &
                                              lsd, comm, state, energy, grid_adjoint)
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

      INTEGER                                            :: axis, np, ns, spin
      REAL(dp), ALLOCATABLE                              :: background(:, :), covectors(:, :), &
                                                            gradient(:, :, :), values(:, :)

      np = SIZE(points, 2)
      ns = SIZE(matrices, 4)
      CPASSERT(SIZE(fields, 4) == 4*ns)
      ALLOCATE (values(np, 4*ns), covectors(np, 4*ns), background(np, ns), gradient(np, 3, ns))
      CALL interpolate_native_grid_batch(fields, grid, cell, points, values, wrap_auxiliary_cell=.TRUE., continuous=.TRUE.)
      DO spin = 1, ns
         background(:, spin) = values(:, 4*spin - 3)
         DO axis = 1, 3
            gradient(:, axis, spin) = values(:, 4*spin - 3 + axis)
         END DO
      END DO
      CALL qs_upf_inverse_xc_evaluate(frames, augmentation, matrices, kweights, points, weights, &
                                      background, gradient, core_density, core_gradient, xc_section, &
                                      lsd, comm, state, energy)
      DO spin = 1, ns
         covectors(:, 4*spin - 3) = weights*state%background_vdensity(:, spin)
         DO axis = 1, 3
            covectors(:, 4*spin - 3 + axis) = weights*state%background_vgradient(:, axis, spin)
         END DO
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
!> \param cell_response Response of the common physical and fully periodic PW lattice.
!> \param reduced_k_response ...
!> \param coefficients_response ...
!> \param core_response ...
!> \param core_gradient_response ...
!> \param comm ...
!> \param response ...
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_xc_grid_response(frames, augmentation, matrices, kweights, points, weights, state, &
                                              fields, fields_response, grid, cell, point_response, weights_response, &
                                              atom_response, cell_response, reduced_k_response, coefficients_response, &
                                              core_response, core_gradient_response, comm, response)
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

      INTEGER                                            :: axis, np, ns, spin
      REAL(dp)                                           :: lattice_error, lattice_tolerance
      REAL(dp), ALLOCATABLE                              :: background(:, :), gradient(:, :, :), &
                                                            values(:, :)

      np = SIZE(points, 2)
      ns = SIZE(matrices, 4)
      CPASSERT(SIZE(fields, 4) == 4*ns)
      CPASSERT(ALL(cell%perd == 1))
      lattice_error = MAXVAL(ABS(grid%dh*SPREAD(REAL(grid%npts, dp), 1, 3) - cell%hmat))
      lattice_tolerance = 128.0_dp*EPSILON(1.0_dp)*MAX(1.0_dp, MAXVAL(ABS(cell%hmat)))
      CPASSERT(lattice_error < lattice_tolerance)
      ALLOCATE (values(np, 4*ns), background(np, ns), gradient(np, 3, ns))
      CALL native_grid_batch_response(fields, fields_response, grid, cell, points, &
                                      point_response, cell_response, values, wrap_auxiliary_cell=.TRUE., continuous=.TRUE.)
      DO spin = 1, ns
         background(:, spin) = values(:, 4*spin - 3)
         DO axis = 1, 3
            gradient(:, axis, spin) = values(:, 4*spin - 3 + axis)
         END DO
      END DO
      CALL qs_upf_inverse_xc_response(frames, augmentation, matrices, kweights, points, weights, state, &
                                      point_response, weights_response, atom_response, cell_response, &
                                      reduced_k_response, coefficients_response, background, gradient, &
                                      core_response, core_gradient_response, comm, response)
   END SUBROUTINE qs_upf_inverse_xc_grid_response

! **************************************************************************************************
!> \brief Replicate the native density and its plane-wave Cartesian gradient for mixed XC.
!> \param density Native spin densities, including any prescribed frozen core.
!> \param fields Replicated descriptors (x,y,z,4*spin), ordered as density,gx,gy,gz.
!>               The gradient is the native Fourier derivative of the supplied grid density.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_xc_pw_fields(density, fields)
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: density(:)
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: fields(:, :, :, :)

      INTEGER                                            :: axis, bounds(2, 3), orders(3), spin
      TYPE(pw_c1d_gs_type)                               :: density_g, derivative_g
      TYPE(pw_grid_type), POINTER                        :: grid
      TYPE(pw_r3d_rs_type)                               :: work_r

      CPASSERT(SIZE(density) > 0)
      grid => density(1)%pw_grid
      CPASSERT(ASSOCIATED(grid))
      bounds = grid%bounds_local - SPREAD(grid%bounds(1, :) - 1, 1, 2)
      ALLOCATE (fields(grid%npts(1), grid%npts(2), grid%npts(3), 4*SIZE(density)))
      fields = 0.0_dp
      CALL density_g%create(grid)
      CALL derivative_g%create(grid)
      CALL work_r%create(grid)
      DO spin = 1, SIZE(density)
         CPASSERT(ASSOCIATED(density(spin)%pw_grid, grid))
         fields(bounds(1, 1):bounds(2, 1), bounds(1, 2):bounds(2, 2), &
                bounds(1, 3):bounds(2, 3), 4*spin - 3) = density(spin)%array
         CALL pw_transfer(density(spin), density_g)
         DO axis = 1, 3
            derivative_g%array(:) = density_g%array
            orders = 0
            orders(axis) = 1
            CALL pw_derive(derivative_g, orders)
            CALL pw_transfer(derivative_g, work_r)
            fields(bounds(1, 1):bounds(2, 1), bounds(1, 2):bounds(2, 2), &
                   bounds(1, 3):bounds(2, 3), 4*spin - 3 + axis) = work_r%array
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
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_xc_pw_adjoint(fields, potential)
      REAL(dp), INTENT(IN)                               :: fields(:, :, :, :)
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: potential(:)

      INTEGER                                            :: axis, bounds(2, 3), orders(3), spin
      TYPE(pw_c1d_gs_type)                               :: derivative_g
      TYPE(pw_grid_type), POINTER                        :: grid
      TYPE(pw_r3d_rs_type)                               :: derivative_r, work_r

      CPASSERT(SIZE(potential) > 0)
      grid => potential(1)%pw_grid
      CPASSERT(ASSOCIATED(grid))
      CPASSERT(ALL(SHAPE(fields) == [grid%npts, 4*SIZE(potential)]))
      CPASSERT(grid%dvol > 0.0_dp)
      bounds = grid%bounds_local - SPREAD(grid%bounds(1, :) - 1, 1, 2)
      CALL derivative_g%create(grid)
      CALL work_r%create(grid)
      CALL derivative_r%create(grid)
      DO spin = 1, SIZE(potential)
         CPASSERT(ASSOCIATED(potential(spin)%pw_grid, grid))
         potential(spin)%array(:, :, :) = fields(bounds(1, 1):bounds(2, 1), bounds(1, 2):bounds(2, 2), &
                                                 bounds(1, 3):bounds(2, 3), 4*spin - 3)
         DO axis = 1, 3
            work_r%array(:, :, :) = fields(bounds(1, 1):bounds(2, 1), bounds(1, 2):bounds(2, 2), &
                                           bounds(1, 3):bounds(2, 3), 4*spin - 3 + axis)
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
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_xc_pw_response(density_response, fields, lattice_response, fields_response)
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: density_response(:)
      REAL(dp), INTENT(IN)                               :: fields(:, :, :, :), &
                                                            lattice_response(3, 3)
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: fields_response(:, :, :, :)

      INTEGER                                            :: axis, other, spin
      REAL(dp)                                           :: lattice_inverse(3, 3), strain(3, 3)
      TYPE(pw_grid_type), POINTER                        :: grid

      CALL qs_upf_inverse_xc_pw_fields(density_response, fields_response)
      CPASSERT(ALL(SHAPE(fields) == SHAPE(fields_response)))
      grid => density_response(1)%pw_grid
      lattice_inverse = grid%dh_inv/SPREAD(REAL(grid%npts, dp), 2, 3)
      strain = MATMUL(lattice_response, lattice_inverse)
      DO spin = 1, SIZE(density_response)
         DO axis = 1, 3
            DO other = 1, 3
               fields_response(:, :, :, 4*spin - 3 + axis) = fields_response(:, :, :, 4*spin - 3 + axis) - &
                                                             strain(other, axis)*fields(:, :, :, 4*spin - 3 + other)
            END DO
         END DO
      END DO
   END SUBROUTINE qs_upf_inverse_xc_pw_response

END MODULE qs_upf_inverse_xc
