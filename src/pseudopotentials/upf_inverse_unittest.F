!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Inverse reconstruction identities, conditioning and analytic field/operator responses.
! **************************************************************************************************
PROGRAM upf_inverse_unittest
   USE atom_upf,                        ONLY: atom_upfpot_type
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE basis_set_container_types,       ONLY: add_basis_set_to_container
   USE basis_set_types,                 ONLY: allocate_gto_basis_set,&
                                              deallocate_gto_basis_set,&
                                              gto_basis_set_type
   USE cell_types,                      ONLY: cell_type
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE gaussian_basis_fields,           ONLY: gaussian_basis_fields_evaluate,&
                                              gaussian_basis_fields_periodic,&
                                              gaussian_basis_workspace_type
   USE gaussian_coulomb_fields,         ONLY: gaussian_coulomb_fields_evaluate,&
                                              gaussian_coulomb_fields_evaluate_batch,&
                                              gaussian_coulomb_fields_init,&
                                              gaussian_coulomb_fields_periodic,&
                                              gaussian_coulomb_fields_prepare,&
                                              gaussian_coulomb_fields_screen,&
                                              gaussian_coulomb_fields_truncate,&
                                              gaussian_coulomb_fields_type
   USE ieee_arithmetic,                 ONLY: ieee_quiet_nan,&
                                              ieee_value
   USE kinds,                           ONLY: dp
   USE lebedev,                         ONLY: get_number_of_lebedev_grid,&
                                              lebedev_grid
   USE mathconstants,                   ONLY: fourpi,&
                                              pi
   USE message_passing,                 ONLY: mp_comm_type,&
                                              mp_world_finalize,&
                                              mp_world_init
   USE orbital_pointers,                ONLY: init_orbital_pointers,&
                                              ncoset
   USE particle_types,                  ONLY: particle_type
   USE paw_core_kinetic,                ONLY: paw_core_kinetic_ae_cartesian,&
                                              paw_core_kinetic_type,&
                                              paw_core_radial_slopes,&
                                              paw_core_radial_value
   USE paw_xc,                          ONLY: paw_xc_ae_core_cartesian,&
                                              paw_xc_type
   USE qs_grid_atom,                    ONLY: grid_atom_type
   USE qs_harmonics_atom,               ONLY: harmonics_atom_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_upf_inverse,                  ONLY: qs_upf_inverse_atomic_operators,&
                                              qs_upf_inverse_atomic_overlap,&
                                              qs_upf_inverse_evaluate,&
                                              qs_upf_inverse_prepare,&
                                              qs_upf_inverse_response
   USE qs_upf_inverse_density,          ONLY: qs_upf_inverse_density_radial,&
                                              qs_upf_inverse_density_radial_adjoint,&
                                              qs_upf_inverse_kinetic_center,&
                                              qs_upf_inverse_kinetic_collocate,&
                                              qs_upf_inverse_kinetic_integrate,&
                                              qs_upf_inverse_kinetic_operator
   USE qs_upf_inverse_paw,              ONLY: qs_upf_inverse_paw_collocate,&
                                              qs_upf_inverse_paw_operator
   USE qs_upf_inverse_types,            ONLY: qs_upf_inverse_augmentation_type,&
                                              qs_upf_inverse_type
   USE qs_uspp_types,                   ONLY: qs_uspp_type
   USE skala_orbital_constraints,       ONLY: skala_constrained_stationarity,&
                                              skala_constraint_chart_create,&
                                              skala_constraint_chart_restore,&
                                              skala_constraint_chart_type,&
                                              skala_cusp_admissibility
   USE t_c_g0,                          ONLY: free_C0,&
                                              init_truncated => init
   USE upf_core,                        ONLY: upf_core_fourier,&
                                              upf_core_kinetic_fourier
   USE upf_inverse,                     ONLY: &
        upf_inverse_atomic_init, upf_inverse_atomic_residual, upf_inverse_atomic_type, &
        upf_inverse_basis, upf_inverse_basis_response, upf_inverse_density_evaluate, &
        upf_inverse_fields_init, upf_inverse_init, upf_inverse_kinetic_evaluate, &
        upf_inverse_kinetic_init, upf_inverse_kinetic_matrix, upf_inverse_kinetic_type, &
        upf_inverse_operator, upf_inverse_operator_response, upf_inverse_project, &
        upf_inverse_project_adjoint, upf_inverse_project_response, upf_inverse_type
   USE upf_inverse_coulomb,             ONLY: &
        upf_inverse_coulomb_center_type, upf_inverse_coulomb_compensate, &
        upf_inverse_coulomb_compensation_adjoint, upf_inverse_coulomb_density, &
        upf_inverse_coulomb_fourier_kernel, upf_inverse_coulomb_map_type, &
        upf_inverse_coulomb_pair_kernel, upf_inverse_coulomb_periodic, &
        upf_inverse_coulomb_potential_type, upf_inverse_coulomb_prepare, &
        upf_inverse_coulomb_radial, upf_inverse_coulomb_radial_response
   USE upf_inverse_lattice,             ONLY: &
        upf_inverse_lattice_adjoint, upf_inverse_lattice_bloch, upf_inverse_lattice_bound, &
        upf_inverse_lattice_error_bounds, upf_inverse_lattice_init, upf_inverse_lattice_project, &
        upf_inverse_lattice_residual, upf_inverse_lattice_response, &
        upf_inverse_lattice_response_residual, upf_inverse_lattice_type
   USE upf_inverse_quadrature,          ONLY: upf_inverse_cross_density,&
                                              upf_inverse_density_residual,&
                                              upf_inverse_pair,&
                                              upf_inverse_projector_density,&
                                              upf_inverse_semilocal_pair,&
                                              upf_inverse_semilocal_residual,&
                                              upf_inverse_semilocal_response,&
                                              upf_inverse_smooth_density
   USE upf_projector_integrals,         ONLY: upf_projector_init,&
                                              upf_projector_overlap,&
                                              upf_projector_type
   USE upf_radial_fields,               ONLY: upf_radial_fields_evaluate,&
                                              upf_radial_fields_periodic,&
                                              upf_radial_fields_type
   USE upf_types,                       ONLY: upf_radial_weights
   USE uspp_augmentation,               ONLY: uspp_augmentation_init,&
                                              uspp_augmentation_prepare_realspace

   IMPLICIT NONE
   TYPE roundoff_cusp_context
      REAL(dp) :: energy_zero = 1.e8_dp, last_coordinate = 0.0_dp
      INTEGER :: queries = 0
      LOGICAL :: interrupt = .TRUE.
   END TYPE roundoff_cusp_context
   REAL(dp), PARAMETER :: step = 1.e-5_dp
   REAL(dp), PARAMETER :: stationarity_tolerance = 1.e-8_dp
   REAL(dp) :: max_error = 0.0_dp
   INTEGER :: stop_polls = 0
   CALL check_inverse()
   CALL check_physical_stationarity()
   CALL check_constrained_curvature()
   CALL check_dependent_constraints()
   CALL check_cusp_normals()
   CALL check_lattice_inverse()
   CALL check_projector_density()
   CALL check_density_residual()
   CALL check_radial_coulomb()
   CALL check_pair_coulomb()
   CALL check_periodic_coulomb()
   CALL check_center_angular()
   CALL check_fourier_coulomb()
   CALL check_core_fourier()
   CALL check_paw_core_fields()
   CALL check_nc_fields()
   CALL check_kinetic_fields()
   CALL check_gaussian_fields()
   WRITE (*, '(A)') 'upf_inverse_unittest PASS'
CONTAINS
! **************************************************************************************************
!> \brief Resolve an anisotropic quadratic on a curved constraint below the energy resolution.
! **************************************************************************************************
   SUBROUTINE check_constrained_curvature()
      REAL(dp), PARAMETER :: tolerance = 1.e-11_dp, weights(3) = [1.0_dp, 0.01_dp, 0.0001_dp]

      INTEGER                                            :: i, iterations
      LOGICAL                                            :: converged, valid
      REAL(dp) :: anisotropic(9), coordinate(4), derivative(9), empty(0), energy, energy_zero, &
         implicit_derivative, jacobian(1, 4), midpoint(4), values(1)
      TYPE(skala_constraint_chart_type)                  :: curve

      coordinate(:3) = [0.01_dp, 0.02_dp, 0.03_dp]
      coordinate(4) = SUM(coordinate(:3)**2)
      energy_zero = 1.e8_dp
      CALL skala_constraint_chart_create(curved_constraints, energy_zero, coordinate, [1.e-16_dp], curve, valid)
      IF (.NOT. valid) ERROR STOP 'Could not define the physical constraint curve'
      ! The midpoint of two feasible orbital-coordinate endpoints is not generally feasible.
      midpoint(:3) = 1.5_dp*coordinate(:3)
      midpoint(4) = 2.5_dp*coordinate(4)
      IF (ABS(midpoint(4) - SUM(midpoint(:3)**2)) < 1.e-5_dp) THEN
         ERROR STOP 'Linear interpolation did not distinguish the curved physical surface'
      END IF
      CALL skala_constraint_chart_restore(curve, curved_constraints, energy_zero, midpoint, values, jacobian, valid)
      IF (.NOT. valid .OR. ABS(values(1)) > 1.e-16_dp) ERROR STOP 'Constraint curve lost physical feasibility'
      IF (ANY(midpoint(:3) /= 1.5_dp*coordinate(:3))) ERROR STOP 'Constraint restoration changed free coordinates'
      IF (ABS(midpoint(4) - SUM(midpoint(:3)**2)) > 1.e-16_dp) ERROR STOP 'Wrong restored physical curve'
      CALL skala_constrained_stationarity(curved_constraints, curved_energy_gradient, energy_zero, coordinate, &
                                          [1.e-16_dp], tolerance, converged, max_iterations=4, iterations=iterations)
      IF (.NOT. converged .OR. iterations > 4 .OR. NORM2(weights*coordinate(:3)) > tolerance) THEN
         ERROR STOP 'Implicit constrained curvature was not resolved below energy resolution'
      END IF
      IF (ABS(coordinate(4) - SUM(coordinate(:3)**2)) > 1.e-16_dp) THEN
         ERROR STOP 'Reduced-gradient curvature changed the physical constraint'
      END IF
      ! Feasibility in field units alone does not resolve a weak-normal orbital curve.
      coordinate = [0.2_dp, 0.0_dp, 0.0_dp, 0.04_dp + 1.e-10_dp]
      CALL skala_constraint_chart_create(scaled_curved_constraints, energy_zero, coordinate, [1.e-17_dp], curve, valid)
      IF (.NOT. valid) ERROR STOP 'Could not define the scaled physical constraint curve'
      CALL skala_constraint_chart_restore(curve, scaled_curved_constraints, energy_zero, coordinate, values, jacobian, valid)
      IF (.NOT. valid .OR. ABS(values(1)) > 1.e-17_dp) ERROR STOP 'Scaled constraint lost physical feasibility'
      ! E=1e4*(h**2-|u|**4)+|u|**2/2 has reduced derivative u on h=|u|**2.
      implicit_derivative = coordinate(1) + 4.e4_dp*coordinate(1)*(coordinate(4) - SUM(coordinate(:3)**2))
      IF (ABS(implicit_derivative - coordinate(1)) > 1.e-11_dp) THEN
         ERROR STOP 'Feasible field windows did not resolve implicit-curve derivatives'
      END IF
      ! Resolve a rotated quadratic before its energy changes reach roundoff.
      DO i = 1, SIZE(anisotropic)
         anisotropic(i) = 0.3_dp*COS(1.7_dp*REAL(i, dp))
      END DO
      energy_zero = 0.0_dp
      CALL skala_constrained_stationarity(empty_constraints, anisotropic_energy_gradient, energy_zero, &
                                          anisotropic, empty, tolerance, converged, max_iterations=200, iterations=iterations)
      CALL anisotropic_energy_gradient(anisotropic, energy, derivative, energy_zero)
      IF (.NOT. converged .OR. iterations > 200 .OR. NORM2(derivative) > tolerance .OR. NORM2(anisotropic) > 1.e-7_dp) THEN
         ERROR STOP 'Anisotropic curvature did not reach the original gradient tolerance'
      END IF
   END SUBROUTINE check_constrained_curvature

! **************************************************************************************************
!> \brief Positive quadratic with rotated axes, eigenvalues 1e-3 through 1e3 and exact minimum zero.
!> \param coordinate Trial coordinates.
!> \param energy Original quadratic energy, including its independent zero.
!> \param gradient Exact analytical gradient.
!> \param context Additive energy zero.
! **************************************************************************************************
   SUBROUTINE anisotropic_energy_gradient(coordinate, energy, gradient, context)
      REAL(dp), INTENT(IN)                               :: coordinate(:)
      REAL(dp), INTENT(OUT)                              :: energy, gradient(:)

      CLASS(*), INTENT(INOUT) :: context

      INTEGER :: i, j, n
      REAL(dp) :: axes(SIZE(coordinate), SIZE(coordinate)), eigenvalues(SIZE(coordinate)), &
                  matrix(SIZE(coordinate), SIZE(coordinate))

      n = SIZE(coordinate)
      DO j = 1, n
         eigenvalues(j) = 10.0_dp**(-3.0_dp + 6.0_dp*REAL(j - 1, dp)/REAL(n - 1, dp))
         DO i = 1, n
            axes(i, j) = SQRT(2.0_dp/REAL(n, dp))*COS(pi*(REAL(i, dp) - 0.5_dp)*REAL(j - 1, dp)/REAL(n, dp))
         END DO
      END DO
      axes(:, 1) = axes(:, 1)/SQRT(2.0_dp)
      matrix = MATMUL(axes*SPREAD(eigenvalues, 1, n), TRANSPOSE(axes))
      gradient = MATMUL(matrix, coordinate)
      SELECT TYPE (context)
      TYPE IS (REAL(dp))
         energy = context + 0.5_dp*DOT_PRODUCT(coordinate, gradient)
      CLASS DEFAULT
         ERROR STOP 'Invalid anisotropic quadratic context'
      END SELECT
   END SUBROUTINE anisotropic_energy_gradient

! **************************************************************************************************
!> \brief A curved surface h=sum(u**2), with u the first three coordinates.
!> \param coordinate Complete trial coordinates (u,h).
!> \param values Original constraint residual.
!> \param jacobian Derivative of the residual.
!> \param context Additive energy zero of the same analytical problem.
! **************************************************************************************************
   SUBROUTINE curved_constraints(coordinate, values, jacobian, context)
      REAL(dp), INTENT(IN)                               :: coordinate(:)
      REAL(dp), INTENT(OUT)                              :: values(:), jacobian(:, :)

      CLASS(*), INTENT(INOUT)                            :: context

      SELECT TYPE (context)
      TYPE IS (REAL(dp))
         values(1) = coordinate(4) - SUM(coordinate(:3)**2)
         jacobian(1, :3) = -2.0_dp*coordinate(:3)
         jacobian(1, 4) = 1.0_dp
      CLASS DEFAULT
         ERROR STOP 'Invalid curved-constraint context'
      END SELECT
   END SUBROUTINE curved_constraints

! **************************************************************************************************
!> \brief The same physical curve expressed in small field units.
!> \param coordinate Complete trial coordinates (u,h).
!> \param values Scaled physical residual.
!> \param jacobian Scaled physical Jacobian.
!> \param context Additive energy zero of the same analytical problem.
! **************************************************************************************************
   SUBROUTINE scaled_curved_constraints(coordinate, values, jacobian, context)
      REAL(dp), INTENT(IN)                               :: coordinate(:)
      REAL(dp), INTENT(OUT)                              :: values(:), jacobian(:, :)

      CLASS(*), INTENT(INOUT) :: context

      CALL curved_constraints(coordinate, values, jacobian, context)
      values = 1.e-8_dp*values
      jacobian = 1.e-8_dp*jacobian
   END SUBROUTINE scaled_curved_constraints

! **************************************************************************************************
!> \brief A reduced quadratic with a nonzero normal derivative and an independent energy zero.
!> \param coordinate Complete trial coordinates (u,h).
!> \param energy E0+sum(weights*u**2)/2+normal*(h-sum(u**2)).
!> \param gradient Complete energy derivative before eliminating h.
!> \param context Additive energy zero E0.
! **************************************************************************************************
   SUBROUTINE curved_energy_gradient(coordinate, energy, gradient, context)
      REAL(dp), INTENT(IN)                               :: coordinate(:)
      REAL(dp), INTENT(OUT)                              :: energy, gradient(:)

      CLASS(*), INTENT(INOUT)                            :: context
      REAL(dp), PARAMETER :: weights(3) = [1.0_dp, 0.01_dp, 0.0001_dp], normal = 0.25_dp

      SELECT TYPE (context)
      TYPE IS (REAL(dp))
         energy = context + 0.5_dp*SUM(weights*coordinate(:3)**2) &
                  + normal*(coordinate(4) - SUM(coordinate(:3)**2))
         gradient(:3) = weights*coordinate(:3) - 2.0_dp*normal*coordinate(:3)
         gradient(4) = normal
      CLASS DEFAULT
         ERROR STOP 'Invalid curved-energy context'
      END SELECT
   END SUBROUTINE curved_energy_gradient

! **************************************************************************************************
!> \brief Retain every physical residual when several samples share an independent condition.
! **************************************************************************************************
   SUBROUTINE check_dependent_constraints()
      REAL(dp), PARAMETER :: initial(2) = [0.2_dp, -0.1_dp], &
         tolerances(4) = [1.e-30_dp, 2.0_dp, 3.e-12_dp, 1.e-12_dp]

      INTEGER                                            :: case_index, problem
      LOGICAL                                            :: converged
      REAL(dp)                                           :: coordinate(2), jacobian(4, 2), values(4)

      DO case_index = 0, 2
         problem = case_index
         coordinate = initial
         CALL skala_constrained_stationarity(dependent_constraints, dependent_energy_gradient, problem, coordinate, &
                                             tolerances, stationarity_tolerance, converged, max_iterations=40)
         IF (problem == 0) THEN
            CALL dependent_constraints(coordinate, values, jacobian, problem)
            IF (.NOT. converged .OR. ANY(ABS(values) > tolerances)) THEN
               ERROR STOP 'Redundant physical samples lost their original feasibility checks'
            END IF
            IF (NORM2(coordinate - [1.0_dp, -1.0_dp]) > stationarity_tolerance) THEN
               ERROR STOP 'Rank reduction changed the constrained stationary point'
            END IF
         ELSE IF (converged .OR. ANY(coordinate /= initial)) THEN
            ERROR STOP 'An incompatible dependent or zero-normal constraint was discarded'
         END IF
      END DO
   END SUBROUTINE check_dependent_constraints

! **************************************************************************************************
!> \brief Supply three proportional samples with different units and a constant fourth sample.
!> \param coordinate Two independent orbital coordinates.
!> \param values All four original physical residuals.
!> \param jacobian Derivative of each residual.
!> \param context Zero for compatible samples, one or two for incompatible dependent or constant data.
! **************************************************************************************************
   SUBROUTINE dependent_constraints(coordinate, values, jacobian, context)
      REAL(dp), INTENT(IN)                               :: coordinate(:)
      REAL(dp), INTENT(OUT)                              :: values(:), jacobian(:, :)

      CLASS(*), INTENT(INOUT)                            :: context
      REAL(dp), PARAMETER                                :: factors(4) = [1.e-18_dp, 2.e12_dp, -3.0_dp, 0.0_dp]

      values(:) = factors*SUM(coordinate)
      jacobian(:, :) = SPREAD(factors, 2, SIZE(coordinate))
      SELECT TYPE (context)
      TYPE IS (INTEGER)
         IF (context == 1) values(2) = values(2) + 20.0_dp
         IF (context == 2) values(4) = 1.e-11_dp
      CLASS DEFAULT
         ERROR STOP 'Invalid dependent-constraint control'
      END SELECT
   END SUBROUTINE dependent_constraints

! **************************************************************************************************
!> \brief A quadratic with free minimum (2,0) and constrained minimum (1,-1).
!> \param coordinate Complete trial coordinates.
!> \param energy Complete energy.
!> \param gradient Energy derivative in the same coordinates.
!> \param context The same physical-residual control used by dependent_constraints.
! **************************************************************************************************
   SUBROUTINE dependent_energy_gradient(coordinate, energy, gradient, context)
      REAL(dp), INTENT(IN)                               :: coordinate(:)
      REAL(dp), INTENT(OUT)                              :: energy, gradient(:)

      CLASS(*), INTENT(INOUT)                            :: context

      SELECT TYPE (context)
      TYPE IS (INTEGER)
         IF (context < 0 .OR. context > 2) ERROR STOP 'Invalid dependent-constraint energy control'
      CLASS DEFAULT
         ERROR STOP 'Invalid dependent-constraint energy context'
      END SELECT
      gradient(:) = coordinate - [2.0_dp, 0.0_dp]
      energy = 0.5_dp*SUM(gradient**2)
   END SUBROUTINE dependent_energy_gradient

! **************************************************************************************************
!> \brief Retain coincident physical samples and reject unattainable independent branch choices.
! **************************************************************************************************
   SUBROUTINE check_cusp_normals()
      INTEGER                                            :: dependent_sides(3), groups(2), sides(2)
      LOGICAL                                            :: admissible(2), dependent_admissible(3), &
                                                            valid
      REAL(dp) :: changes(1, 2), coupled_changes(2, 2), coupled_jacobian(2, 2), &
         dependent_changes(2, 3), dependent_jacobian(3, 2), jacobian(2, 1)

      jacobian(:, 1) = [2.e-18_dp, 4.e12_dp]
      changes(1, :) = [-1.0_dp, -4.0_dp]
      CALL skala_cusp_admissibility([3.0_dp], changes, jacobian, admissible, sides, valid)
      IF (.NOT. valid .OR. .NOT. ALL(admissible) .OR. ANY(sides /= 0)) THEN
         ERROR STOP 'Coincident density normals were discarded or affected by their physical units'
      END IF
      jacobian(:, 1) = [1.0_dp, -1.0_dp]
      changes(1, :) = [-4.0_dp, 2.0_dp]
      CALL skala_cusp_admissibility([3.0_dp], changes, jacobian, admissible, sides, valid)
      IF (.NOT. valid .OR. .NOT. ALL(admissible)) ERROR STOP 'Opposite density normals did not switch together'
      ! Independent weights can cancel -1.5 + 2*lambda_1 + 2*lambda_2. Both attainable
      ! physical sides have slope +0.5, requiring negative-side descent of the first sample.
      changes = 2.0_dp
      CALL skala_cusp_admissibility([-1.5_dp], changes, jacobian, admissible, sides, valid)
      IF (.NOT. valid .OR. ANY(admissible) .OR. ANY(sides /= [-1, 1])) THEN
         ERROR STOP 'An unattainable convex cancellation was accepted as a physical minimum'
      END IF
      jacobian(:, 1) = [1.0_dp, 2.0_dp]
      changes(1, :) = [1.0_dp, 1.0_dp]
      CALL skala_cusp_admissibility([-1.0_dp], changes, jacobian, admissible, sides, valid)
      IF (.NOT. valid .OR. ANY(admissible) .OR. ANY(sides /= 1)) ERROR STOP 'A concave coincident cusp was accepted'
      ! Dual normals of x and x+y preserve the other sample. The first cusp brackets
      ! zero, while the second admits positive-side descent without releasing the first.
      coupled_jacobian(1, :) = [1.0_dp, 0.0_dp]
      coupled_jacobian(2, :) = [1.0_dp, 1.0_dp]
      coupled_changes(:, 1) = [-2.0_dp, 0.0_dp]
      coupled_changes(:, 2) = [0.0_dp, -1.0_dp]
      CALL skala_cusp_admissibility([-1.0_dp, -2.0_dp], coupled_changes, coupled_jacobian, &
                                    admissible, sides, valid, group_ids=groups)
      IF (.NOT. valid .OR. .NOT. admissible(1) .OR. admissible(2) .OR. &
          ANY(sides /= [0, 1]) .OR. ANY(groups /= [1, 2])) THEN
         ERROR STOP 'A coupled normal release changed an admissible independent sample'
      END IF
      dependent_jacobian(:, 1) = [1.0_dp, 0.0_dp, 1.0_dp]
      dependent_jacobian(:, 2) = [0.0_dp, 1.0_dp, 1.0_dp]
      dependent_changes = 0.0_dp
      CALL skala_cusp_admissibility([0.0_dp, 0.0_dp], dependent_changes, dependent_jacobian, &
                                    dependent_admissible, dependent_sides, valid)
      IF (valid .OR. ANY(dependent_admissible) .OR. ANY(dependent_sides /= 0)) THEN
         ERROR STOP 'Coupled dependent normals were declared independent'
      END IF
   END SUBROUTINE check_cusp_normals

! **************************************************************************************************
!> \brief A scaled physical coordinate must satisfy its own gradient bound and iteration budget.
! **************************************************************************************************
   SUBROUTINE check_physical_stationarity()
      INTEGER                                            :: iterations
      LOGICAL                                            :: converged, interrupted, stopped
      REAL(dp)                                           :: coordinate(1), cusp_coordinate(2), &
                                                            empty(0), pair(2), scale
      TYPE(roundoff_cusp_context)                        :: cusp

      scale = 20.0_dp
      coordinate = 1.e-6_dp
      CALL skala_constrained_stationarity(empty_constraints, scaled_energy_gradient, scale, coordinate, &
                                          empty, stationarity_tolerance, converged, max_iterations=8, iterations=iterations)
      IF (.NOT. converged .OR. iterations /= 0) ERROR STOP 'Initial coordinate derivative should satisfy its bound'
      IF (NORM2(coordinate)/scale <= stationarity_tolerance) ERROR STOP 'Physical gradient control is not distinct'
      CALL skala_constrained_stationarity(empty_constraints, scaled_energy_gradient, scale, coordinate, &
                                          empty, stationarity_tolerance, converged, max_iterations=1, &
                                          verify_stationarity=scaled_physical_gradient, iterations=iterations)
      IF (converged .OR. iterations /= 1 .OR. coordinate(1) /= 1.e-6_dp) THEN
         ERROR STOP 'Physical stationarity exceeded its budget or changed a failed trial'
      END IF
      CALL skala_constrained_stationarity(empty_constraints, scaled_energy_gradient, scale, coordinate, &
                                          empty, stationarity_tolerance, converged, max_iterations=8, &
                                          verify_stationarity=scaled_physical_gradient, iterations=iterations, &
                                          curvature_steps=.TRUE., precondition=scaled_inverse_curvature)
      IF (.NOT. converged .OR. iterations <= 0 .OR. iterations > 8 .OR. ABS(coordinate(1)) > 1.e-12_dp) THEN
         ERROR STOP 'Known physical minimum was not established'
      END IF
      coordinate = 1.0_dp
      CALL skala_constrained_stationarity(empty_constraints, scaled_energy_gradient, scale, coordinate, &
                                          empty, stationarity_tolerance, converged, max_iterations=1, &
                                          verify_stationarity=scaled_physical_gradient, iterations=iterations, &
                                          retain_progress=.TRUE.)
      IF (converged .OR. iterations /= 1) ERROR STOP 'Retained progress was confused with convergence'
      IF (ABS(coordinate(1) - (1.0_dp - 1.0_dp/scale**2)) > 1.e-14_dp) THEN
         ERROR STOP 'An incomplete minimization discarded its accepted state'
      END IF
      coordinate = 1.0_dp
      stop_polls = 0
      CALL skala_constrained_stationarity(empty_constraints, scaled_energy_gradient, scale, coordinate, &
                                          empty, stationarity_tolerance, converged, max_iterations=32, &
                                          verify_stationarity=scaled_physical_gradient, iterations=iterations, &
                                          check_stop=stop_after_first_step, stopped=stopped)
      IF (converged .OR. .NOT. stopped .OR. iterations /= 1) ERROR STOP 'Requested stop was lost'
      IF (ABS(coordinate(1) - (1.0_dp - 1.0_dp/scale**2)) > 1.e-14_dp) THEN
         ERROR STOP 'Requested stop discarded the accepted state'
      END IF
      CALL skala_constrained_stationarity(empty_constraints, scaled_energy_gradient, scale, coordinate, &
                                          empty, stationarity_tolerance, converged, max_iterations=32, &
                                          verify_stationarity=scaled_physical_gradient, stopped=stopped)
      IF (.NOT. converged .OR. stopped .OR. ABS(coordinate(1)) > 1.e-12_dp) THEN
         ERROR STOP 'Could not resume a stopped orbital minimization'
      END IF
      coordinate = 1.0_dp
      CALL skala_constrained_stationarity(empty_constraints, scaled_energy_gradient, scale, coordinate, &
                                          empty, stationarity_tolerance, converged, max_iterations=32, &
                                          verify_stationarity=scaled_physical_gradient, iterations=iterations, &
                                          check_segment=crossing_segment, interrupted=interrupted, stopped=stopped)
      IF (converged .OR. stopped .OR. .NOT. interrupted .OR. iterations /= 2) THEN
         ERROR STOP 'A crossing event was confused with convergence or an external stop'
      END IF
      IF (ABS(coordinate(1) - (1.0_dp - 1.0_dp/scale**2)) > 1.e-14_dp) THEN
         ERROR STOP 'A crossing event lost the last accepted feasible state'
      END IF
      scale = 1.0_dp
      coordinate = 0.0_dp
      CALL skala_constrained_stationarity(empty_constraints, absolute_energy_gradient, scale, coordinate, &
                                          empty, stationarity_tolerance, converged, max_iterations=4, &
                                          iterations=iterations, interrupted=interrupted, &
                                          check_stalled_segment=stalled_cusp_segment)
      IF (converged .OR. .NOT. interrupted .OR. iterations /= 1 .OR. coordinate(1) /= 0.0_dp) THEN
         ERROR STOP 'A nonsmooth minimum lost the original rejected trial or was declared differentiably stationary'
      END IF
      ! Detect the first directional reversal across a tilted cusp while its shallow
      ! tangential energy decrease is still resolved, before the smooth search stalls.
      scale = 1.e-3_dp
      cusp_coordinate = [0.7_dp, 0.5_dp]
      CALL skala_constrained_stationarity(empty_constraints, absolute_energy_gradient, scale, cusp_coordinate, &
                                          empty, stationarity_tolerance, converged, max_iterations=4, &
                                          iterations=iterations, interrupted=interrupted, &
                                          check_stalled_segment=crossing_segment)
      IF (converged .OR. .NOT. interrupted) THEN
         ERROR STOP 'A persistent cusp derivative exhausted the smooth search before boundary detection'
      END IF
      coordinate = 1.e-8_dp
      CALL skala_constrained_stationarity(empty_constraints, roundoff_cusp_gradient, cusp, coordinate, &
                                          empty, stationarity_tolerance, converged, max_iterations=1, &
                                          iterations=iterations, interrupted=interrupted, &
                                          check_stalled_segment=roundoff_cusp_segment)
      IF (converged .OR. .NOT. interrupted .OR. iterations /= 1 .OR. &
          coordinate(1) /= 1.e-8_dp .OR. cusp%queries /= 1) THEN
         ERROR STOP 'An unresolved energy change hid a physical cusp behind accepted shortened trials'
      END IF
      cusp%queries = 0
      cusp%interrupt = .FALSE.
      coordinate = 1.e-8_dp
      CALL skala_constrained_stationarity(empty_constraints, roundoff_cusp_gradient, cusp, coordinate, &
                                          empty, stationarity_tolerance, converged, max_iterations=1, &
                                          iterations=iterations, interrupted=interrupted, retain_progress=.TRUE., &
                                          check_stalled_segment=roundoff_cusp_segment)
      IF (converged .OR. interrupted .OR. iterations /= 1 .OR. cusp%queries /= 1 .OR. &
          coordinate(1) == 1.e-8_dp .OR. coordinate(1) /= cusp%last_coordinate) THEN
         ERROR STOP 'A rejected boundary query lost the accepted native trial state'
      END IF
      scale = 1.e-8_dp
      coordinate = 1.e-12_dp
      CALL skala_constrained_stationarity(empty_constraints, scaled_energy_gradient, scale, coordinate, &
                                          empty, stationarity_tolerance, converged, max_iterations=8, iterations=iterations)
      IF (.NOT. converged .OR. iterations > 8 .OR. NORM2(coordinate)/scale**2 > stationarity_tolerance) THEN
         ERROR STOP 'A steep quadratic exhausted its line search above coordinate resolution'
      END IF
      ! A numerical coordinate map must preserve the original derivative bound.
      scale = 1.e-6_dp
      coordinate = 0.5_dp*stationarity_tolerance*scale
      CALL skala_constrained_stationarity(empty_constraints, scaled_energy_gradient, scale, coordinate, &
                                          empty, stationarity_tolerance, converged, max_iterations=1, iterations=iterations, &
                                          coordinate_scales=[1.0_dp/scale])
      IF (.NOT. converged .OR. iterations /= 0) ERROR STOP 'Preconditioning changed the original stopping bound'
      scale = 1.e6_dp
      coordinate = 2.0_dp*stationarity_tolerance*scale
      CALL skala_constrained_stationarity(empty_constraints, scaled_energy_gradient, scale, coordinate, &
                                          empty, stationarity_tolerance, converged, max_iterations=1, iterations=iterations, &
                                          coordinate_scales=[1.0_dp/scale])
      IF (converged .OR. iterations /= 1) ERROR STOP 'Small scaled gradients hid a nonstationary original coordinate'
      scale = 1.0_dp
      coordinate = 1.0_dp
      CALL skala_constrained_stationarity(empty_constraints, scaled_energy_gradient, scale, coordinate, &
                                          empty, stationarity_tolerance, converged, max_iterations=1, &
                                          iterations=iterations, interrupted=interrupted, curvature_steps=.TRUE., &
                                          check_curvature_segment=interrupt_curvature_probe)
      IF (converged .OR. .NOT. interrupted .OR. iterations /= 1 .OR. coordinate(1) /= 1.0_dp) THEN
         ERROR STOP 'A coupled-curvature boundary lost the accepted point or bypassed its iteration budget'
      END IF
      scale = 1.0_dp
      pair = [0.3_dp, -0.2_dp]
      CALL skala_constrained_stationarity(empty_constraints, scaled_energy_gradient, scale, pair, &
                                          empty, stationarity_tolerance, converged, max_iterations=8, max_step=1.0_dp, &
                                          curvature_steps=.TRUE., precondition=extreme_inverse_curvature)
      IF (.NOT. converged .OR. NORM2(pair) > stationarity_tolerance) THEN
         ERROR STOP 'An unresolved inverse response prevented a known coupled minimum'
      END IF
   END SUBROUTINE check_physical_stationarity

! **************************************************************************************************
!> \brief Detect a known threshold along a scaled physical orbital trial.
!> \param first Last accepted coordinate.
!> \param last Proposed coordinate.
!> \param curve Fixed physical constraint chart.
!> \param context Positive coordinate scale.
!> \param interrupted Whether the trial crosses the threshold.
! **************************************************************************************************
   SUBROUTINE scaled_inverse_curvature(coordinate, context, inverse_diagonal, available)
      REAL(dp), INTENT(IN)                               :: coordinate(:)

      CLASS(*), INTENT(INOUT) :: context
      REAL(dp), INTENT(OUT) :: inverse_diagonal(:)
      LOGICAL, INTENT(OUT) :: available

      available = .TRUE.
      SELECT TYPE (context)
      TYPE IS (REAL(dp))
         inverse_diagonal = context**2
      CLASS DEFAULT
         ERROR STOP 'Invalid curvature-preconditioner context'
      END SELECT
   END SUBROUTINE scaled_inverse_curvature

   SUBROUTINE extreme_inverse_curvature(coordinate, context, inverse_diagonal, available)
      REAL(dp), INTENT(IN)                               :: coordinate(:)

      CLASS(*), INTENT(INOUT) :: context
      REAL(dp), INTENT(OUT) :: inverse_diagonal(:)
      LOGICAL, INTENT(OUT) :: available

      IF (SIZE(coordinate) /= 2) ERROR STOP 'Wrong coupled inverse-response control'
      inverse_diagonal = [1.0_dp, HUGE(1.0_dp)/2.0_dp]
      available = .TRUE.
   END SUBROUTINE extreme_inverse_curvature

   SUBROUTINE interrupt_curvature_probe(first, last, curve, context, smooth, interrupted)
      REAL(dp), INTENT(IN)                               :: first(:), last(:)
      TYPE(skala_constraint_chart_type), INTENT(IN)      :: curve

      CLASS(*), INTENT(INOUT) :: context
      LOGICAL, INTENT(OUT) :: smooth, interrupted

      smooth = .FALSE.
      interrupted = ANY(first /= last)
   END SUBROUTINE interrupt_curvature_probe

   SUBROUTINE crossing_segment(first, last, curve, context, interrupted)
      REAL(dp), INTENT(IN)                               :: first(:), last(:)
      TYPE(skala_constraint_chart_type), INTENT(IN)      :: curve

      CLASS(*), INTENT(INOUT)                            :: context
      LOGICAL, INTENT(OUT)                               :: interrupted

      SELECT TYPE (context)
      TYPE IS (REAL(dp))
         IF (context <= 0.0_dp) ERROR STOP 'Nonpositive crossing coordinate scale'
      CLASS DEFAULT
         ERROR STOP 'Wrong crossing coordinate scale'
      END SELECT
      IF (SIZE(first) == 2) THEN
         interrupted = (first(2) - 0.3_dp*first(1))*(last(2) - 0.3_dp*last(1)) <= 0.0_dp
      ELSE
         interrupted = first(1) > 0.9_dp .AND. last(1) <= 0.9_dp
      END IF
   END SUBROUTINE crossing_segment

! **************************************************************************************************
!> \brief Absolute-value energy whose ordinary derivative cannot vanish at its minimum.
!> \param coordinate One absolute-value coordinate or two with a tilted cusp and shallow tangential curvature.
!> \param energy Exact nonsmooth energy.
!> \param gradient Positive-side derivative at zero and the ordinary derivative elsewhere.
!> \param context Positive absolute-value scale in one dimension, tangential curvature in two.
! **************************************************************************************************
   SUBROUTINE absolute_energy_gradient(coordinate, energy, gradient, context)
      REAL(dp), INTENT(IN)                               :: coordinate(:)
      REAL(dp), INTENT(OUT)                              :: energy, gradient(:)

      CLASS(*), INTENT(INOUT) :: context

      SELECT TYPE (context)
      TYPE IS (REAL(dp))
         IF (SIZE(coordinate) == 2) THEN
            energy = 0.5_dp*context*coordinate(1)**2 + ABS(coordinate(2) - 0.3_dp*coordinate(1))
            gradient(2) = SIGN(1.0_dp, coordinate(2) - 0.3_dp*coordinate(1))
            gradient(1) = context*coordinate(1) - 0.3_dp*gradient(2)
         ELSE
            energy = context*ABS(coordinate(1))
            gradient(1) = SIGN(context, coordinate(1))
         END IF
      CLASS DEFAULT
         ERROR STOP 'Invalid absolute-value context'
      END SELECT
   END SUBROUTINE absolute_energy_gradient

! **************************************************************************************************
!> \brief A failed line search must retain its first feasible segment around the cusp.
!> \param first Accepted physical cusp coordinate.
!> \param last Original feasible trial, before backtracking collapses its interval.
!> \param curve Original physical constraint partition.
!> \param context Positive energy scale.
!> \param interrupted Return the original accepted coordinate for boundary treatment.
! **************************************************************************************************
   SUBROUTINE stalled_cusp_segment(first, last, curve, context, interrupted)
      REAL(dp), INTENT(IN)                               :: first(:), last(:)
      TYPE(skala_constraint_chart_type), INTENT(IN)      :: curve

      CLASS(*), INTENT(INOUT) :: context
      LOGICAL, INTENT(OUT) :: interrupted

      IF (first(1) /= 0.0_dp .OR. ABS(last(1) + 0.1_dp) > 1.e-14_dp) THEN
         ERROR STOP 'Stalled cusp search discarded its original feasible trial segment'
      END IF
      interrupted = .TRUE.
   END SUBROUTINE stalled_cusp_segment

! **************************************************************************************************
!> \brief A cusp near an independent energy zero, with decreasing gradient norm toward its minimum.
!> \param coordinate One real orbital coordinate.
!> \param energy E0+abs(x)+x**2/2.
!> \param gradient Ordinary one-sided derivative.
!> \param context Additive energy zero and last evaluated native coordinate.
! **************************************************************************************************
   SUBROUTINE roundoff_cusp_gradient(coordinate, energy, gradient, context)
      REAL(dp), INTENT(IN)                               :: coordinate(:)
      REAL(dp), INTENT(OUT)                              :: energy, gradient(:)

      CLASS(*), INTENT(INOUT)                            :: context

      SELECT TYPE (context)
      TYPE IS (roundoff_cusp_context)
         energy = context%energy_zero + ABS(coordinate(1)) + 0.5_dp*coordinate(1)**2
         gradient(1) = SIGN(1.0_dp, coordinate(1)) + coordinate(1)
         context%last_coordinate = coordinate(1)
      CLASS DEFAULT
         ERROR STOP 'Invalid unresolved-energy cusp context'
      END SELECT
   END SUBROUTINE roundoff_cusp_gradient

! **************************************************************************************************
!> \brief Preserve the original wide segment and model a query restoring the previous native state.
!> \param first Last accepted coordinate.
!> \param last First feasible trial before shortening its interval.
!> \param curve Original physical constraint partition.
!> \param context Query count, requested interruption and last evaluated native coordinate.
!> \param interrupted Whether to return the last accepted state for physical boundary treatment.
! **************************************************************************************************
   SUBROUTINE roundoff_cusp_segment(first, last, curve, context, interrupted)
      REAL(dp), INTENT(IN)                               :: first(:), last(:)
      TYPE(skala_constraint_chart_type), INTENT(IN)      :: curve

      CLASS(*), INTENT(INOUT)                            :: context
      LOGICAL, INTENT(OUT)                               :: interrupted

      SELECT TYPE (context)
      TYPE IS (roundoff_cusp_context)
         IF (first(1) /= 1.e-8_dp .OR. ABS(last(1) - first(1) + 0.1_dp) > 1.e-14_dp) THEN
            ERROR STOP 'An unresolved-energy cusp query lost its first feasible segment'
         END IF
         context%queries = context%queries + 1
         context%last_coordinate = first(1)
         interrupted = context%interrupt
      CLASS DEFAULT
         ERROR STOP 'Invalid unresolved-energy cusp context'
      END SELECT
   END SUBROUTINE roundoff_cusp_segment

! **************************************************************************************************
!> \brief Request termination after the first accepted line-search trial.
!> \param context Positive coordinate scale.
!> \param requested Whether one complete step has been accepted.
! **************************************************************************************************
   SUBROUTINE stop_after_first_step(context, requested)
      CLASS(*), INTENT(INOUT) :: context
      LOGICAL, INTENT(OUT) :: requested

      SELECT TYPE (context)
      TYPE IS (REAL(dp))
         IF (context <= 0.0_dp) ERROR STOP 'Nonpositive stopped coordinate scale'
      CLASS DEFAULT
         ERROR STOP 'Wrong stopped coordinate scale'
      END SELECT
      stop_polls = stop_polls + 1
      requested = stop_polls >= 3
   END SUBROUTINE stop_after_first_step

! **************************************************************************************************
!> \brief Supply an unconstrained coordinate problem.
!> \param coordinate Trial coordinates.
!> \param values Empty residuals.
!> \param jacobian Empty residual derivatives.
!> \param context Energy-problem data, unused by the empty constraints.
! **************************************************************************************************
   SUBROUTINE empty_constraints(coordinate, values, jacobian, context)
      REAL(dp), INTENT(IN)                               :: coordinate(:)
      REAL(dp), INTENT(OUT)                              :: values(:), jacobian(:, :)

      CLASS(*), INTENT(INOUT) :: context

      IF (SIZE(coordinate) < 1 .OR. SIZE(values) /= 0 .OR. ANY(SHAPE(jacobian) /= [0, SIZE(coordinate)])) THEN
         ERROR STOP 'Wrong unconstrained-coordinate control'
      END IF
      values = 0.0_dp
      jacobian = 0.0_dp
   END SUBROUTINE empty_constraints

! **************************************************************************************************
!> \brief Differentiate E=2+r^2/2 for the scaled coordinate x=scale*r.
!> \param coordinate Scaled trial coordinate.
!> \param energy Physical energy.
!> \param gradient Derivative with respect to x.
!> \param context Coordinate scale.
! **************************************************************************************************
   SUBROUTINE scaled_energy_gradient(coordinate, energy, gradient, context)
      REAL(dp), INTENT(IN)                               :: coordinate(:)
      REAL(dp), INTENT(OUT)                              :: energy, gradient(:)

      CLASS(*), INTENT(INOUT) :: context

      SELECT TYPE (context)
      TYPE IS (REAL(dp))
         gradient(:) = coordinate/context**2
         energy = 2.0_dp + 0.5_dp*DOT_PRODUCT(coordinate, gradient)
      CLASS DEFAULT
         ERROR STOP 'Wrong physical coordinate scale'
      END SELECT
   END SUBROUTINE scaled_energy_gradient

! **************************************************************************************************
!> \brief Check the derivative with respect to the original physical coordinate r.
!> \param coordinate Scaled trial coordinate.
!> \param context Coordinate scale.
!> \param converged Whether the physical energy derivative satisfies its bound.
! **************************************************************************************************
   SUBROUTINE scaled_physical_gradient(coordinate, context, converged)
      REAL(dp), INTENT(IN)                               :: coordinate(:)

      CLASS(*), INTENT(INOUT) :: context
      LOGICAL, INTENT(OUT) :: converged

      SELECT TYPE (context)
      TYPE IS (REAL(dp))
         converged = NORM2(coordinate)/context <= stationarity_tolerance
      CLASS DEFAULT
         ERROR STOP 'Wrong physical coordinate scale'
      END SELECT
   END SUBROUTINE scaled_physical_gradient

! **************************************************************************************************
!> \brief Physical-image inverse against a noncommuting matrix geometric sequence.
!>        Checks the exterior equation, Bloch limit, coupled response and real energy adjoint.
! **************************************************************************************************
   SUBROUTINE check_lattice_inverse()
      INTEGER, PARAMETER                                 :: nao = 2, np = 2, nr = 129, &
                                                            domains(3) = [3, 8, nr]

      COMPLEX(dp) :: b(np, nao, nr), bloch(np, nao), c(np, np, 3), c_padded(np, np, 16), &
         covector(np, nao, nr), db(np, nao, nr), db_bloch(np, nao), dc(np, np, 3), &
         dc_bloch(np, np), difference(np, nao, nr), exact(np, nao, nr), exact_bloch(np, nao), &
         exact_response(np, nao, nr + 1), gb(np, nao, nr), gc(np, np, 3), inverse(np, np), &
         k(np, np), response(np, nao, nr), samples(np, nao, nr, 2), solution(np, nao, nr)
      COMPLEX(dp), ALLOCATABLE                           :: residual(:, :, :), &
                                                            residual_samples(:, :, :, :), &
                                                            response_residual(:, :, :)
      INTEGER                                            :: i, images(3, nr), info, &
                                                            kernel_images(3, 3), &
                                                            padded_images(3, 16), padding, side, &
                                                            tail
      INTEGER, ALLOCATABLE                               :: residual_images(:, :), &
                                                            response_images(:, :)
      REAL(dp) :: bloch_bound, coefficient_bound, error, lower_bound, previous, previous_bound, &
         reduced_k(3), response_bloch_bound, response_bound, shift, tail_bound(2)
      TYPE(upf_inverse_lattice_type)                     :: lattice, shifted

      kernel_images(:, 1) = [0, 0, 0]
      kernel_images(:, 2) = [1, -2, 1]
      kernel_images(:, 3) = -kernel_images(:, 2)
      DO i = 1, nr
         images(:, i) = (i - 1)*kernel_images(:, 2)
      END DO
      c = CMPLX(0.0_dp, 0.0_dp, dp)
      c(:, :, 1) = RESHAPE([CMPLX(0.2_dp, 0.04_dp, dp), CMPLX(-0.1_dp, 0.03_dp, dp), &
                            CMPLX(0.07_dp, -0.08_dp, dp), CMPLX(-0.12_dp, 0.06_dp, dp)], [np, np])
      c(:, :, 2) = RESHAPE([CMPLX(0.12_dp, -0.03_dp, dp), CMPLX(0.08_dp, 0.01_dp, dp), &
                            CMPLX(-0.04_dp, 0.09_dp, dp), CMPLX(0.17_dp, 0.02_dp, dp)], [np, np])
      b = CMPLX(0.0_dp, 0.0_dp, dp)
      b(:, :, 1) = RESHAPE([CMPLX(0.8_dp, -0.3_dp, dp), CMPLX(-0.2_dp, 0.1_dp, dp), &
                            CMPLX(0.3_dp, 0.4_dp, dp), CMPLX(0.6_dp, -0.2_dp, dp)], [np, nao])
      k = c(:, :, 1)
      k(1, 1) = k(1, 1) + 1.0_dp
      k(2, 2) = k(2, 2) + 1.0_dp
      CALL invert_two(k, inverse)
      exact(:, :, 1) = MATMUL(inverse, b(:, :, 1))
      DO i = 2, nr
         exact(:, :, i) = -MATMUL(inverse, MATMUL(c(:, :, 2), exact(:, :, i - 1)))
      END DO
      previous = HUGE(0.0_dp)
      previous_bound = HUGE(0.0_dp)
      DO padding = 1, 3
         i = domains(padding)
         CALL upf_inverse_lattice_init(kernel_images, c, images(:, :i), lattice, info)
         CALL check(ABS(REAL(info, dp)), "Physical-image preconditioner")
         CALL upf_inverse_lattice_project(lattice, b(:, :, :i), solution(:, :, :i), 2.e-14_dp, info)
         CALL check(ABS(REAL(info, dp)), "Coupled physical-image solve")
         CALL check(MAXVAL(ABS(solution(:, :, :i) - exact(:, :, :i))), "Noncommuting inverse image sequence")
         CALL upf_inverse_lattice_residual(lattice, solution(:, :, :i), images(:, :1), b(:, :, :1), &
                                           residual_images, residual)
         tail = 0
         DO side = 1, SIZE(residual_images, 2)
            IF (ALL(residual_images(:, side) == i*kernel_images(:, 2))) tail = side
         END DO
         CALL check(MERGE(0.0_dp, 1.0_dp, tail > 0), "Exterior convolution equation retained", 0.0_dp)
         CALL check(MAXVAL(ABS(residual(:, :, tail) - MATMUL(c(:, :, 2), exact(:, :, i)))), &
                    "Independent exterior residual")
         error = MAXVAL(ABS(residual))
         CALL check(MERGE(0.0_dp, 1.0_dp, error < previous), "Image padding reduces the full residual", 0.0_dp)
         previous = error
         CALL upf_inverse_lattice_error_bounds(lattice, solution(:, :, :i), images(:, :1), b(:, :, :1), &
                                               coefficient_bound, bloch_bound, lower_bound, info)
         CALL check(ABS(REAL(info, dp)), "Uniform inverse bound established")
         difference = exact
         difference(:, :, :i) = difference(:, :, :i) - solution(:, :, :i)
         error = SQRT(SUM(ABS(difference)**2))
         CALL check(MAX(0.0_dp, error - coefficient_bound), "Full coefficient tail bounded", 0.0_dp)
         CALL check(MERGE(0.0_dp, 1.0_dp, bloch_bound < previous_bound), &
                    "Image padding reduces the uniform Bloch bound", 0.0_dp)
         previous_bound = bloch_bound
         reduced_k = [0.17_dp, 0.11_dp, -0.09_dp]
         CALL upf_inverse_lattice_bloch(kernel_images, c, reduced_k, k)
         k(1, 1) = k(1, 1) + 1.0_dp
         k(2, 2) = k(2, 2) + 1.0_dp
         CALL invert_two(k, inverse)
         exact_bloch = MATMUL(inverse, b(:, :, 1))
         CALL upf_inverse_lattice_bloch(images(:, :i), solution(:, :, :i), reduced_k, bloch)
         error = SQRT(SUM(ABS(bloch - exact_bloch)**2))
         CALL check(MAX(0.0_dp, error - bloch_bound), "Pointwise Bloch error bounded", 0.0_dp)
      END DO
      c_padded = CMPLX(0.0_dp, 0.0_dp, dp)
      c_padded(:, :, :3) = c
      padded_images(:, :3) = kernel_images
      DO i = 4, 16
         padded_images(:, i) = [i, 3*i, 0]
      END DO
      CALL upf_inverse_lattice_init(padded_images, c_padded, images, shifted, info)
      CALL check(ABS(REAL(info, dp)), "Zero kernel images preserve the preconditioner")
      CALL upf_inverse_lattice_project(shifted, b, samples(:, :, :, 1), 2.e-14_dp, info)
      CALL check(ABS(REAL(info, dp)), "Zero kernel images preserve the solve")
      CALL check(MAXVAL(ABS(samples(:, :, :, 1) - solution)), "Zero kernels preserve the physical operator")
      CALL upf_inverse_lattice_error_bounds(shifted, solution, images(:, :1), b(:, :, :1), &
                                            error, previous, lower_bound, info)
      CALL check(ABS(REAL(info, dp)), "Zero kernel images retain the uniform inverse bound")
      CALL check(ABS(error/coefficient_bound - 1.0_dp), "Zero kernels preserve the coefficient error bound", 1.e-9_dp)
      CALL check(ABS(previous/bloch_bound - 1.0_dp), "Zero kernels preserve the Bloch error bound", 1.e-9_dp)
      CALL upf_inverse_lattice_init(kernel_images(:, 3:1:-1), c(:, :, 3:1:-1), images(:, nr:1:-1), shifted, info)
      CALL check(ABS(REAL(info, dp)), "Reordered physical-image preconditioner")
      CALL upf_inverse_lattice_project(shifted, b(:, :, nr:1:-1), samples(:, :, :, 1), 2.e-14_dp, info)
      CALL check(ABS(REAL(info, dp)), "Reordered coupled image solve")
      CALL check(MAXVAL(ABS(samples(:, :, :, 1) - exact(:, :, nr:1:-1))), "Inverse preserves physical image order")
      DO side = 1, 2
         reduced_k = 0.0_dp
         IF (side == 2) reduced_k = [0.17_dp, 0.11_dp, -0.09_dp]
         CALL upf_inverse_lattice_bloch(kernel_images, c, reduced_k, k)
         k(1, 1) = k(1, 1) + 1.0_dp
         k(2, 2) = k(2, 2) + 1.0_dp
         CALL invert_two(k, inverse)
         exact_bloch = MATMUL(inverse, b(:, :, 1))
         CALL upf_inverse_lattice_bloch(images, solution, reduced_k, bloch)
         CALL check(MAXVAL(ABS(bloch - exact_bloch)), "Unaliased image sum agrees with the Bloch inverse")
         IF (side == 1) CALL check(MERGE(0.0_dp, 1.0_dp, MAXVAL(ABS(bloch - solution(:, :, 1))) > 0.01_dp), &
                                   "Gamma coefficients differ from the individual home image", 0.0_dp)
      END DO
      DO i = 1, nr
         db(:, :, i) = CMPLX(0.013_dp*COS(REAL(i, dp)), -0.009_dp*SIN(REAL(i, dp)), dp)
         covector(:, :, i) = CMPLX(0.07_dp*SIN(REAL(i, dp)), 0.11_dp*COS(REAL(i, dp)), dp)
         covector(1, 2, i) = CONJG(covector(1, 2, i))*0.3_dp
      END DO
      dc(:, :, 1) = c(:, :, 2)*0.2_dp
      dc(:, :, 2) = CONJG(TRANSPOSE(c(:, :, 1)))*0.11_dp
      dc(:, :, 3) = c(:, :, 1)*CMPLX(0.03_dp, -0.02_dp, dp)
      CALL upf_inverse_lattice_response(lattice, solution, db, dc, response, 2.e-14_dp, info)
      CALL check(ABS(REAL(info, dp)), "Coupled image response converged")
      DO side = 1, 2
         shift = REAL(2*side - 3, dp)*step
         CALL upf_inverse_lattice_init(kernel_images, c + shift*dc, images, shifted, info)
         CALL check(ABS(REAL(info, dp)), "Perturbed image preconditioner")
         CALL upf_inverse_lattice_project(shifted, b + shift*db, samples(:, :, :, side), 2.e-14_dp, info)
         CALL check(ABS(REAL(info, dp)), "Perturbed coupled image solve")
      END DO
      CALL check(MAXVAL(ABS((samples(:, :, :, 2) - samples(:, :, :, 1))/(2.0_dp*step) - response)), &
                 "Coupled image response against fresh inverse solves", 3.e-9_dp)
      CALL upf_inverse_lattice_response_residual(lattice, solution, response, images, db, dc, &
                                                 response_images, response_residual)
      CALL check(MAXVAL(ABS(response_residual(:, :, :nr))), "Interior coupled response equations")
      CALL upf_inverse_lattice_error_bounds(lattice, solution, images, b, coefficient_bound, bloch_bound, &
                                            lower_bound, info, response=response, overlaps_response=db, &
                                            coupling_response=dc, response_error=response_bound, &
                                            response_bloch_error=response_bloch_bound)
      CALL check(ABS(REAL(info, dp)), "Full coupled response bound established")
      k = c(:, :, 1)
      k(1, 1) = k(1, 1) + 1.0_dp
      k(2, 2) = k(2, 2) + 1.0_dp
      CALL invert_two(k, inverse)
      exact_response(:, :, 1) = -MATMUL(inverse, MATMUL(dc(:, :, 3), exact(:, :, 1)))
      DO i = 1, nr
         db_bloch = db(:, :, i) - MATMUL(dc(:, :, 1), exact(:, :, i))
         IF (i > 1) db_bloch = db_bloch - MATMUL(dc(:, :, 2), exact(:, :, i - 1))
         bloch = -MATMUL(inverse, MATMUL(c(:, :, 2), exact(:, :, i)))
         db_bloch = db_bloch - MATMUL(dc(:, :, 3), bloch)
         exact_response(:, :, i + 1) = MATMUL(inverse, db_bloch - MATMUL(c(:, :, 2), exact_response(:, :, i)))
      END DO
      difference = exact_response(:, :, 2:) - response
      error = SQRT(SUM(ABS(difference)**2) + SUM(ABS(exact_response(:, :, 1))**2))
      CALL check(MAX(0.0_dp, error - response_bound), "Full response coefficient tail bounded", 0.0_dp)
      CALL upf_inverse_lattice_bloch(kernel_images, c, reduced_k, k)
      k(1, 1) = k(1, 1) + 1.0_dp
      k(2, 2) = k(2, 2) + 1.0_dp
      CALL invert_two(k, inverse)
      exact_bloch = MATMUL(inverse, b(:, :, 1))
      CALL upf_inverse_lattice_bloch(images, db, reduced_k, db_bloch)
      CALL upf_inverse_lattice_bloch(kernel_images, dc, reduced_k, dc_bloch)
      exact_bloch = MATMUL(inverse, db_bloch - MATMUL(dc_bloch, exact_bloch))
      CALL upf_inverse_lattice_bloch(images, response, reduced_k, bloch)
      CALL check(MAX(0.0_dp, SQRT(SUM(ABS(bloch - exact_bloch)**2)) - response_bloch_bound), &
                 "Exterior Bloch response error bounded", 0.0_dp)
      tail_bound = [0.02_dp, 0.03_dp]
      error = coefficient_bound
      previous = bloch_bound
      CALL upf_inverse_lattice_error_bounds(lattice, solution, images, b, coefficient_bound, bloch_bound, &
                                            lower_bound, info, overlap_tail=tail_bound)
      CALL check(ABS(coefficient_bound - error - tail_bound(1)/lower_bound), "Omitted B l2 tail retained")
      CALL check(ABS(bloch_bound - previous - tail_bound(2)/lower_bound), "Omitted B l1 tail retained")
      ALLOCATE (residual_samples(np, nao, SIZE(response_images, 2), 2))
      DO side = 1, 2
         shift = REAL(2*side - 3, dp)*step
         CALL upf_inverse_lattice_init(kernel_images, c + shift*dc, images, shifted, info)
         CALL upf_inverse_lattice_residual(shifted, samples(:, :, :, side), images, b + shift*db, &
                                           residual_images, residual)
         CALL check(REAL(MAXVAL(ABS(residual_images - response_images)), dp), "Identical physical residual support", 0.0_dp)
         residual_samples(:, :, :, side) = residual
      END DO
      CALL check(MAXVAL(ABS((residual_samples(:, :, :, 2) - residual_samples(:, :, :, 1))/(2.0_dp*step) - response_residual)), &
                 "Full exterior response residual against fresh kernels", 3.e-9_dp)
      CALL upf_inverse_lattice_adjoint(lattice, solution, covector, gb, gc, 2.e-14_dp, info)
      CALL check(ABS(REAL(info, dp)), "Coupled image adjoint converged")
      CALL check(ABS(REAL(SUM(CONJG(covector)*response) - SUM(CONJG(gb)*db) - SUM(CONJG(gc)*dc), dp)), &
                 "Complete shared-kernel image adjoint")
      CALL upf_inverse_lattice_project(lattice, b, samples(:, :, :, 1), 2.e-14_dp, info, max_iterations=1)
      CALL check(ABS(REAL(info - 1, dp)), "Unconverged lattice solve reports failure", 0.0_dp)
      CALL upf_inverse_lattice_residual(lattice, samples(:, :, :, 1), images, b, residual_images, residual)
      CALL check(MERGE(0.0_dp, 1.0_dp, MAXVAL(ABS(residual(:, :, :nr))) > 1.e-4_dp), &
                 "Failed solve retains its true interior residual", 0.0_dp)
      c = CMPLX(0.0_dp, 0.0_dp, dp)
      c(1, 1, 2) = 0.8_dp
      c(2, 2, 2) = 0.8_dp
      CALL upf_inverse_lattice_init(kernel_images, c, images(:, :1), shifted, info)
      CALL upf_inverse_lattice_bound(shifted, lower_bound, info)
      CALL check(ABS(REAL(info - 1, dp)) + ABS(lower_bound), "Inconclusive bound reports no certificate", 0.0_dp)
      CALL upf_inverse_lattice_project(shifted, b(:, :, :1), samples(:, :, :1, 1), 2.e-14_dp, info)
      CALL check(ABS(REAL(info, dp)), "Inconclusive uniform bound permits the finite solve")
      CALL upf_inverse_lattice_error_bounds(shifted, samples(:, :, :1, 1), images(:, :1), b(:, :, :1), &
                                            coefficient_bound, bloch_bound, lower_bound, info)
      CALL check(MERGE(0.0_dp, 1.0_dp, info == 1 .AND. coefficient_bound == HUGE(0.0_dp) .AND. &
                       bloch_bound == HUGE(0.0_dp)), "Uncertified errors are not reported as zero", 0.0_dp)

   END SUBROUTINE check_lattice_inverse
! **************************************************************************************************
!> \brief Independent two-projector inverse for the lattice geometric-series control.
!> \param matrix ...
!> \param inverse_matrix ...
! **************************************************************************************************
   SUBROUTINE invert_two(matrix, inverse_matrix)
      COMPLEX(dp), INTENT(IN)                            :: matrix(2, 2)
      COMPLEX(dp), INTENT(OUT)                           :: inverse_matrix(2, 2)

      inverse_matrix = RESHAPE([matrix(2, 2), -matrix(2, 1), -matrix(1, 2), matrix(1, 1)], [2, 2])/ &
                       (matrix(1, 1)*matrix(2, 2) - matrix(1, 2)*matrix(2, 1))
   END SUBROUTINE invert_two
! **************************************************************************************************
!> \brief Cusped radial AE charge and signed-orbital kinetic fields with independent exact derivatives.
! **************************************************************************************************
   SUBROUTINE check_paw_core_fields()
      INTEGER                                            :: axis, order, point, side, support
      REAL(dp) :: curvature, density, derivative, exact(3), fd(5, 2), gradient(3), hessian(3, 3), &
         hessian_only(3, 3), motion(3), points(3, 3), radius, raw_slopes(6), samples(5, 2), &
         second, shift, slope, tau, tau_gradient(3), unit(3), value, x(3)
      TYPE(paw_core_kinetic_type)                        :: cache, empty
      TYPE(paw_xc_type)                                  :: grid

      ALLOCATE (grid%r(6), grid%core(6, 2), grid%core_gradient(6, 2))
      grid%r(:) = [0.01_dp, 0.5_dp, 1.0_dp, 1.5_dp, 2.0_dp, 2.5_dp]
      grid%core(:, 1) = MAX(2.0_dp - grid%r, 0.0_dp)**3
      grid%core_gradient(:, 1) = -3.0_dp*MAX(2.0_dp - grid%r, 0.0_dp)**2
      grid%core(:, 2) = 0.0_dp
      grid%core_gradient(:, 2) = 0.0_dp
      grid%has_ae_core = .TRUE.
      grid%core_support_index = 6
      ALLOCATE (cache%r, SOURCE=grid%r)
      ALLOCATE (cache%wave(6, 1), cache%slope(6, 1), cache%angular(1), cache%occupations(1), &
                cache%endpoint_curvature(1), cache%support_index(1))
      cache%endpoint_curvature = 0.0_dp
      cache%support_index = 6
      cache%wave(:, 1) = grid%core(:, 1)
      cache%slope(:, 1) = grid%core_gradient(:, 1)
      cache%angular = 0
      cache%occupations = 2
      cache%available = .TRUE.
      points(:, 1) = [0.0_dp, 0.0_dp, 0.7_dp]
      points(:, 2) = [0.31_dp, -0.42_dp, 0.53_dp]
      points(:, 3) = [0.003_dp, -0.004_dp, 0.0_dp]
      motion = [0.02_dp, -0.01_dp, 0.03_dp]
      DO point = 1, 3
         x = points(:, point)
         radius = SQRT(SUM(x**2))
         unit = x/radius
         slope = -3.0_dp*(2.0_dp - radius)**2
         CALL paw_xc_ae_core_cartesian(grid, x, density, gradient, hessian)
         CALL check(ABS(density - (2.0_dp - radius)**3), "Cusped AE core density")
         CALL check(MAXVAL(ABS(gradient - slope*unit)), "Cusped AE core gradient")
         DO axis = 1, 3
            exact = (6.0_dp*(2.0_dp - radius) - slope/radius)*unit*unit(axis)
            exact(axis) = exact(axis) + slope/radius
            CALL check(MAXVAL(ABS(hessian(:, axis) - exact)), "Cusped AE core Hessian", 3.e-11_dp)
         END DO
         CALL paw_xc_ae_core_cartesian(grid, x, density, hessian=hessian_only)
         CALL check(MAXVAL(ABS(hessian_only - hessian)), "AE core optional Hessian", 0.0_dp)
         CALL paw_core_kinetic_ae_cartesian(cache, x, tau, tau_gradient)
         CALL check(ABS(tau - 9.0_dp*(2.0_dp - radius)**4/fourpi), "Original signed AE core kinetic field")
         exact = -36.0_dp*(2.0_dp - radius)**3*unit/fourpi
         CALL check(MAXVAL(ABS(tau_gradient - exact)), "AE core kinetic gradient")
         DO order = 1, 2
            DO side = 1, 2
               shift = REAL(2*side - 3, dp)*step/REAL(order, dp)
               CALL paw_xc_ae_core_cartesian(grid, x + shift*motion, samples(1, side), samples(2:4, side))
               CALL paw_core_kinetic_ae_cartesian(cache, x + shift*motion, samples(5, side))
            END DO
            fd(:, order) = (samples(:, 2) - samples(:, 1))/(2.0_dp*step/REAL(order, dp))
         END DO
         exact = (4.0_dp*fd(2:4, 2) - fd(2:4, 1))/3.0_dp
         CALL check(ABS((4.0_dp*fd(1, 2) - fd(1, 1))/3.0_dp - DOT_PRODUCT(gradient, motion)), &
                    "Moving AE core density", 3.e-8_dp)
         CALL check(MAXVAL(ABS(exact - MATMUL(hessian, motion))), &
                    "Moving AE core gradient", 3.e-8_dp)
         CALL check(ABS((4.0_dp*fd(5, 2) - fd(5, 1))/3.0_dp - DOT_PRODUCT(tau_gradient, motion)), &
                    "Moving AE core kinetic field", 3.e-8_dp)
      END DO
      CALL paw_xc_ae_core_cartesian(grid, [3.0_dp, 0.0_dp, 0.0_dp], density, gradient, hessian)
      CALL paw_core_kinetic_ae_cartesian(cache, [3.0_dp, 0.0_dp, 0.0_dp], tau, tau_gradient)
      CALL check(ABS(density) + MAXVAL(ABS(gradient)) + MAXVAL(ABS(hessian)) + &
                 ABS(tau) + MAXVAL(ABS(tau_gradient)), "Verified original AE zero suffix", 0.0_dp)
      CALL paw_core_radial_slopes(grid%r, grid%core(:, 1), raw_slopes, curvature, support)
      CALL check(ABS(REAL(support - 5, dp)), "Core support starts at original zero suffix", 0.0_dp)
      CALL paw_core_radial_value(grid%r, grid%core(:, 1), raw_slopes, curvature, support, &
                                 grid%r(support - 1), value, derivative, second)
      CALL check(ABS(value - grid%core(support - 1, 1)) + ABS(derivative - raw_slopes(support - 1)) + &
                 ABS(second - curvature), "Core closure left value slope and curvature", 3.e-12_dp)
      DO point = 1, 6
         CALL paw_core_radial_value(grid%r, grid%core(:, 1), raw_slopes, curvature, support, &
                                    grid%r(point), value, derivative, second)
         CALL check(ABS(value - grid%core(point, 1)), "Unchanged original core samples", 3.e-12_dp)
      END DO
      CALL paw_core_radial_value(grid%r, grid%core(:, 1), raw_slopes, curvature, support, &
                                 grid%r(6), value, derivative, second)
      CALL check(ABS(value) + ABS(derivative) + ABS(second), "Core C2 zero endpoint", 0.0_dp)
      radius = 0.37_dp*grid%r(support - 1) + 0.63_dp*grid%r(support)
      CALL paw_core_radial_value(grid%r, grid%core(:, 1), raw_slopes, curvature, support, radius, value, derivative, second)
      DO side = 1, 2
         shift = REAL(2*side - 3, dp)*step
         CALL paw_core_radial_value(grid%r, grid%core(:, 1), raw_slopes, curvature, support, radius + shift, &
                                    samples(1, side), samples(2, side), samples(3, side))
      END DO
      CALL check(ABS((samples(1, 2) - samples(1, 1))/(2.0_dp*step) - derivative), "Core quintic first derivative", 3.e-8_dp)
      CALL check(ABS((samples(2, 2) - samples(2, 1))/(2.0_dp*step) - second), "Core quintic second derivative", 3.e-8_dp)
      grid%core = 0.0_dp
      grid%core_gradient = 0.0_dp
      grid%has_ae_core = .FALSE.
      empty%available = .TRUE.
      CALL paw_xc_ae_core_cartesian(grid, [0.0_dp, 0.0_dp, 0.0_dp], density, gradient, hessian)
      CALL paw_core_kinetic_ae_cartesian(empty, [0.0_dp, 0.0_dp, 0.0_dp], tau, tau_gradient)
      CALL check(ABS(density) + MAXVAL(ABS(gradient)) + MAXVAL(ABS(hessian)) + &
                 ABS(tau) + MAXVAL(ABS(tau_gradient)), "Core-free AE origin", 0.0_dp)
   END SUBROUTINE check_paw_core_fields
! **************************************************************************************************
!> \brief Atomic hard-projector orientation, full kinetic matrix and exact NC reduction.
! **************************************************************************************************
   SUBROUTINE check_atomic_projectors(comm)
      TYPE(mp_comm_type), INTENT(IN)                     :: comm

      INTEGER, PARAMETER                                 :: nb = 3, np = 5, nr = 33

      COMPLEX(dp) :: coefficient_overlap(np, 2), coefficient_response(np, 2), &
         delta_kinetic(np, 2), delta_overlap(np, 2), exact_operator(2, 2), identity(nb, nb), &
         kinetic_overlap(np, 2), kinetic_response(np, 2), kinetic_samples(np, 2, 2), &
         laplacian(np, 2), local_samples(2, 2, 2), metric_samples(2, 2, 2), mixed(2, 2), &
         overlap(np, 2), response(np, 2), rotation(nb, nb), samples(np, 2, 2), &
         shifted_overlap(np, 2)
      COMPLEX(dp), ALLOCATABLE :: kinetic_residual(:, :), local_kinetic(:, :), &
         local_kinetic_response(:, :), local_metric(:, :), local_metric_response(:, :), &
         metric_residual(:, :)
      INTEGER                                            :: axis, i, images(3, 3), info, ir, j, side
      REAL(dp) :: ao_center(3), cell(3, 3), center(3), dao(3), dc(3, 3), dcenter(3), delta, &
         directions(3, 6), displacement(3), dk(3), expected(np, np), full_kinetic(np, np), &
         hard_overlap(1, np, 4), kernel(np, np), norm_difference, original_overlap(1, np, 4), r, &
         reduced_k(3), shape
      REAL(dp), ALLOCATABLE                              :: radial_kinetic(:, :), weights(:)
      TYPE(atom_upfpot_type)                             :: pot
      TYPE(gto_basis_set_type), POINTER                  :: basis
      TYPE(qs_upf_inverse_type)                          :: frame, shifted_frame
      TYPE(upf_inverse_atomic_type)                      :: atomic
      TYPE(upf_inverse_kinetic_type)                     :: kinetic
      TYPE(upf_projector_type)                           :: original

      pot%pseudo_type = "USPP"
      pot%is_ultrasoft = .TRUE.
      pot%number_of_proj = nb
      ALLOCATE (pot%r(nr), pot%rab(nr), pot%lbeta(nb), pot%beta(nr, nb), pot%dion(nb, nb), pot%qij(nb, nb), pot%partial_waves)
      ALLOCATE (pot%partial_waves%ae_wfc(nr, nb), pot%partial_waves%ps_wfc(nr, nb), &
                pot%partial_waves%angular_momenta(nb))
      pot%r(:) = [(REAL(i - 1, dp)/16.0_dp, i=1, nr)]
      pot%rab = 1.0_dp/16.0_dp
      pot%lbeta(:) = [0, 0, 1]
      pot%partial_waves%angular_momenta(:) = pot%lbeta
      pot%dion = 0.0_dp
      pot%dion(1:2, 1:2) = RESHAPE([0.8_dp, -0.2_dp, -0.2_dp, 0.5_dp], [2, 2])
      pot%dion(3, 3) = 0.3_dp
      DO j = 1, nb
         DO ir = 1, nr
            r = pot%r(ir)
            shape = 0.0_dp
            IF (r < 1.5_dp) shape = (1.0_dp - (r/1.5_dp)**2)**3
            pot%beta(ir, j) = 0.3_dp*REAL(j + 1, dp)*r**pot%lbeta(j)*shape*(1.0_dp + 0.07_dp*j*r*r)
            pot%partial_waves%ps_wfc(ir, j) = r**(pot%lbeta(j) + 1)*(0.7_dp + 0.1_dp*j*r*r)
            delta = 0.2_dp*r**(pot%lbeta(j) + 1)*shape/REAL(j, dp)
            pot%partial_waves%ae_wfc(ir, j) = pot%partial_waves%ps_wfc(ir, j) + delta
         END DO
      END DO
      pot%qij = 0.0_dp
      DO j = 1, nb
         DO i = 1, nb
            IF (pot%lbeta(i) /= pot%lbeta(j)) CYCLE
            pot%qij(i, j) = SUM(upf_radial_weights(pot)*( &
                                pot%partial_waves%ae_wfc(:, i)*pot%partial_waves%ae_wfc(:, j) - &
                                pot%partial_waves%ps_wfc(:, i)*pot%partial_waves%ps_wfc(:, j))) + 0.01_dp*i*j
         END DO
      END DO
      CALL upf_inverse_atomic_init(pot, 3, atomic, info)
      IF (info /= 0 .OR. .NOT. atomic%inverse%ready) ERROR STOP 'Atomic inverse factorization failed'
      CALL upf_projector_init(pot, 1, original)
      identity = CMPLX(0.0_dp, 0.0_dp, dp)
      DO i = 1, nb
         identity(i, i) = CMPLX(1.0_dp, 0.0_dp, dp)
      END DO
      CALL upf_inverse_project(atomic%inverse, identity, rotation)
      displacement = [0.2_dp, -0.1_dp, 0.3_dp]
      CALL upf_projector_overlap(original, 0, [0.4_dp], displacement, 1, original_overlap)
      CALL upf_projector_overlap(atomic%hard, 0, [0.4_dp], displacement, 1, hard_overlap)
      DO axis = 1, 4
         CALL check(MAXVAL(ABS(hard_overlap(1, :2, axis) - &
                               MATMUL(REAL(rotation(:2, :2), dp), original_overlap(1, :2, axis)))), &
                    'Hard projector bras use K inverse including Gaussian derivatives')
      END DO
      CALL upf_inverse_kinetic_init(pot, kinetic)
      CALL upf_inverse_kinetic_matrix(pot, kinetic, radial_kinetic)
      directions = 0.0_dp
      DO axis = 1, 3
         directions(axis, 2*axis - 1) = 1.0_dp
         directions(axis, 2*axis) = -1.0_dp
      END DO
      ALLOCATE (weights, SOURCE=upf_radial_weights(pot)*pot%r**2)
      full_kinetic = 0.0_dp
      DO ir = 1, nr
         DO side = 1, 6
            CALL upf_inverse_kinetic_evaluate(kinetic, pot%r(ir)*directions(:, side), kernel)
            full_kinetic = full_kinetic + fourpi*weights(ir)*kernel/6.0_dp
         END DO
      END DO
      CALL check(MAXVAL(ABS(full_kinetic - atomic%kinetic)), 'Kinetic matrix angular integration and magnetic selection')
      IF (ABS(radial_kinetic(1, 2)) <= 1.e-6_dp) ERROR STOP 'Missing off-diagonal atomic kinetic correction'
      expected = original%d - full_kinetic
      CALL check(MAXVAL(ABS(atomic%hard%d - expected)), 'Reduced hard coupling is original D minus full kinetic matrix')
      expected = 0.0_dp
      DO j = 1, np
         DO i = 1, np
            IF (original%angular_index(i) /= original%angular_index(j)) CYCLE
            norm_difference = 0.01_dp*original%radial_index(i)*original%radial_index(j)
            expected(i, j) = norm_difference
         END DO
      END DO
      CALL check(MAXVAL(ABS(atomic%metric - expected)), 'Original norm residual is retained after hard reconstruction')
      NULLIFY (basis)
      CALL allocate_gto_basis_set(basis)
      basis%nset = 1
      basis%nsgf = 2
      ALLOCATE (basis%lmax(1), basis%npgf(1), basis%nsgf_set(1), basis%nshell(1), basis%l(1, 1), basis%first_sgf(1, 1), &
                basis%last_sgf(1, 1), basis%zet(2, 1), basis%sphi(2, 2))
      basis%lmax = 0
      basis%npgf = 2
      basis%nsgf_set = 2
      basis%nshell = 1
      basis%l = 0
      basis%first_sgf = 1
      basis%last_sgf = 2
      basis%zet(:, 1) = [0.4_dp, 0.9_dp]
      basis%sphi(:, :) = RESHAPE([0.6_dp, -0.2_dp, 0.3_dp, 0.8_dp], [2, 2])
      images = 0
      images(1, :) = [-1, 0, 1]
      cell = 0.0_dp
      DO axis = 1, 3
         cell(axis, axis) = 1.7_dp + 0.2_dp*axis
      END DO
      dc = RESHAPE([(0.002_dp*i, i=1, 9)], [3, 3])
      center = [0.2_dp, -0.1_dp, 0.3_dp]
      ao_center = [-0.1_dp, 0.3_dp, 0.2_dp]
      dcenter = [0.01_dp, -0.02_dp, 0.03_dp]
      dao = [-0.02_dp, 0.03_dp, 0.01_dp]
      reduced_k = [0.23_dp, -0.17_dp, 0.11_dp]
      dk = [0.004_dp, -0.002_dp, 0.003_dp]
      CALL qs_upf_inverse_atomic_overlap(basis, atomic%hard, center, ao_center, cell, images, reduced_k, &
                                         overlap, response, dcenter, dao, dc, dk, kinetic_overlap, kinetic_response)
      DO side = 1, 2
         delta = REAL(2*side - 3, dp)*step
         CALL qs_upf_inverse_atomic_overlap(basis, atomic%hard, center + delta*dcenter, ao_center + delta*dao, &
                                            cell + delta*dc, images, reduced_k + delta*dk, samples(:, :, side), &
                                            kinetic=kinetic_samples(:, :, side))
      END DO
      CALL check(MAXVAL(ABS((samples(:, :, 2) - samples(:, :, 1))/(2.0_dp*step) - response)), &
                 'Atomic hard projector joint atom cell and complex Bloch response', 3.e-9_dp)
      CALL check(MAXVAL(ABS((kinetic_samples(:, :, 2) - kinetic_samples(:, :, 1))/(2.0_dp*step) - kinetic_response)), &
                 'Atomic Gaussian kinetic projector joint geometry and Bloch response', 3.e-9_dp)
      laplacian = CMPLX(0.0_dp, 0.0_dp, dp)
      DO axis = 1, 3
         DO side = 1, 2
            displacement = ao_center
            displacement(axis) = displacement(axis) + REAL(2*side - 3, dp)*1.e-3_dp
            CALL qs_upf_inverse_atomic_overlap(basis, atomic%hard, center, displacement, cell, images, reduced_k, shifted_overlap)
            laplacian = laplacian + shifted_overlap - overlap
         END DO
      END DO
      CALL check(MAXVAL(ABS(-0.5_dp*laplacian/1.e-6_dp - kinetic_overlap)), &
                 'Gaussian kinetic overlap from independent AO-center second differences', 2.e-7_dp)
      CALL qs_upf_inverse_atomic_overlap(basis, atomic%difference, center, ao_center, cell, images, reduced_k, &
                                         delta_overlap, kinetic=delta_kinetic)
      CALL upf_inverse_atomic_residual(atomic, overlap, delta_overlap, delta_kinetic, kinetic_residual, metric_residual)
      mixed = MATMUL(CONJG(TRANSPOSE(delta_kinetic)), overlap)
      exact_operator = -mixed - CONJG(TRANSPOSE(mixed)) + &
                       MATMUL(CONJG(TRANSPOSE(overlap)), MATMUL(atomic%difference_kinetic + original%d, overlap))
      CALL check(MAXVAL(ABS(kinetic_residual + MATMUL(CONJG(TRANSPOSE(overlap)), MATMUL(atomic%hard%d, overlap)) - &
                            exact_operator)), 'Reduced hard D minus T plus closure residual equals the explicit kinetic pullback')
      mixed = MATMUL(CONJG(TRANSPOSE(delta_overlap)), overlap)
      exact_operator = -mixed - CONJG(TRANSPOSE(mixed)) + &
                       MATMUL(CONJG(TRANSPOSE(overlap)), &
                              MATMUL(atomic%difference_overlap + atomic%wave_overlap + atomic%metric, overlap))
      CALL check(MAXVAL(ABS(metric_residual + MATMUL(CONJG(TRANSPOSE(overlap)), MATMUL(atomic%metric, overlap)) - &
                            exact_operator)), 'Original q metric is preserved with finite local closure residual')
      IF (MAXVAL(ABS(kinetic_residual)) <= 1.e-6_dp) ERROR STOP 'Finite AO closure residual was removed'
      frame%atom_local = .TRUE.
      frame%cell = cell
      frame%reduced_k = reduced_k
      ALLOCATE (frame%atomic_cache(1), frame%bases(1), frame%positions(3, 1), &
                frame%center_atoms(1), frame%center_kinds(1), frame%projector_first(1), frame%projector_last(1), &
                frame%ao_first(1), frame%ao_last(1))
      ALLOCATE (frame%coefficients(np, 2), frame%q(np, np))
      ALLOCATE (frame%images, SOURCE=images)
      frame%atomic_cache(1) = atomic
      frame%bases(1)%gto_basis_set => basis
      frame%positions(:, 1) = center
      frame%center_atoms = 1
      frame%center_kinds = 1
      frame%projector_first = 1
      frame%projector_last = np
      frame%ao_first = 1
      frame%ao_last = 2
      CALL qs_upf_inverse_atomic_overlap(basis, atomic%hard, center, center, cell, images, reduced_k, frame%coefficients)
      frame%q(:, :) = CMPLX(atomic%wave_overlap + atomic%metric, 0.0_dp, dp)
      CALL qs_upf_inverse_atomic_operators(frame, frame%q, local_kinetic, local_metric)
      CALL qs_upf_inverse_atomic_overlap(basis, atomic%difference, center, center, cell, images, reduced_k, &
                                         delta_overlap, kinetic=delta_kinetic)
      mixed = MATMUL(CONJG(TRANSPOSE(delta_kinetic)), frame%coefficients)
      exact_operator = -mixed - CONJG(TRANSPOSE(mixed)) + &
                       MATMUL(CONJG(TRANSPOSE(frame%coefficients)), MATMUL(atomic%difference_kinetic, frame%coefficients))
      CALL check(MAXVAL(ABS(local_kinetic - exact_operator)), 'Native atomic kinetic correction preserves the finite-map remainder')
      mixed = MATMUL(CONJG(TRANSPOSE(delta_overlap)), frame%coefficients)
      exact_operator = -mixed - CONJG(TRANSPOSE(mixed)) + &
                       MATMUL(CONJG(TRANSPOSE(frame%coefficients)), &
                              MATMUL(atomic%difference_overlap + frame%q, frame%coefficients))
      CALL check(MAXVAL(ABS(local_metric - exact_operator)), 'Native atomic metric correction retains the supplied original q')
      CALL qs_upf_inverse_atomic_operators(frame, frame%q, local_kinetic, local_metric, &
                                           local_kinetic_response, local_metric_response, &
                                           RESHAPE(dcenter, [3, 1]), 0.0_dp*cell, 0.0_dp*reduced_k)
      CALL check(MAXVAL(ABS(local_kinetic_response)) + MAXVAL(ABS(local_metric_response)), &
                 'Rigid translation leaves atom-local kinetic and metric corrections unchanged', 0.0_dp)
      CALL qs_upf_inverse_atomic_overlap(basis, atomic%hard, center, center, cell, images, reduced_k, &
                                         coefficient_overlap, coefficient_response, dcenter, dcenter, dc, dk)
      CALL qs_upf_inverse_atomic_operators(frame, frame%q, local_kinetic, local_metric, &
                                           local_kinetic_response, local_metric_response, &
                                           RESHAPE(dcenter, [3, 1]), dc, dk, coefficient_response)
      local_samples(:, :, 1) = local_kinetic_response
      metric_samples(:, :, 1) = local_metric_response
      CALL qs_upf_inverse_atomic_operators(frame, frame%q, local_kinetic, local_metric, &
                                           local_kinetic_response, local_metric_response, &
                                           RESHAPE(dcenter, [3, 1]), dc, dk)
      CALL check(MAXVAL(ABS(local_kinetic_response - local_samples(:, :, 1))) + &
                 MAXVAL(ABS(local_metric_response - metric_samples(:, :, 1))), &
                 'Reused coefficients preserve nonzero lattice and Bloch operator responses', 0.0_dp)
      DO side = 1, 2
         delta = REAL(2*side - 3, dp)*step
         shifted_frame = frame
         shifted_frame%positions(:, 1) = center + delta*dcenter
         shifted_frame%cell = cell + delta*dc
         shifted_frame%reduced_k = reduced_k + delta*dk
         CALL qs_upf_inverse_atomic_overlap(basis, atomic%hard, shifted_frame%positions(:, 1), &
                                            shifted_frame%positions(:, 1), shifted_frame%cell, images, &
                                            shifted_frame%reduced_k, shifted_frame%coefficients)
         CALL qs_upf_inverse_atomic_operators(shifted_frame, shifted_frame%q, local_kinetic, local_metric)
         local_samples(:, :, side) = local_kinetic
         metric_samples(:, :, side) = local_metric
      END DO
      CALL check(MAXVAL(ABS((local_samples(:, :, 2) - local_samples(:, :, 1))/(2.0_dp*step) - &
                            local_kinetic_response)), 'Nonzero lattice and Bloch kinetic response after rigid motion', 3.e-9_dp)
      CALL check(MAXVAL(ABS((metric_samples(:, :, 2) - metric_samples(:, :, 1))/(2.0_dp*step) - &
                            local_metric_response)), 'Nonzero lattice and Bloch metric response after rigid motion', 3.e-9_dp)
      ! A second Gaussian center without a projector tests relative motion independently of lattice motion.
      DEALLOCATE (frame%bases, frame%positions, frame%ao_first, frame%ao_last, frame%coefficients)
      ALLOCATE (frame%bases(2), frame%positions(3, 2), frame%ao_first(2), frame%ao_last(2), frame%coefficients(np, 4))
      frame%bases(1)%gto_basis_set => basis
      frame%bases(2)%gto_basis_set => basis
      frame%positions(:, 1) = center
      frame%positions(:, 2) = ao_center
      frame%ao_first = [1, 3]
      frame%ao_last = [2, 4]
      CALL qs_upf_inverse_atomic_overlap(basis, atomic%hard, center, center, cell, images, reduced_k, &
                                         frame%coefficients(:, 1:2))
      CALL qs_upf_inverse_atomic_overlap(basis, atomic%hard, center, ao_center, cell, images, reduced_k, &
                                         frame%coefficients(:, 3:4))
      CALL qs_upf_inverse_atomic_operators(frame, frame%q, local_kinetic, local_metric, &
                                           local_kinetic_response, local_metric_response, &
                                           SPREAD(dcenter, 2, 2), 0.0_dp*cell, 0.0_dp*reduced_k)
      CALL check(MAXVAL(ABS(local_kinetic_response)) + MAXVAL(ABS(local_metric_response)), &
                 'Common translation with a separate Gaussian center', 0.0_dp)
      CALL qs_upf_inverse_atomic_operators(frame, frame%q, local_kinetic, local_metric, &
                                           local_kinetic_response, local_metric_response, &
                                           RESHAPE([dcenter, dao], [3, 2]), 0.0_dp*cell, 0.0_dp*reduced_k)
      DO side = 1, 2
         delta = REAL(2*side - 3, dp)*step
         shifted_frame = frame
         shifted_frame%positions(:, 1) = center + delta*dcenter
         shifted_frame%positions(:, 2) = ao_center + delta*dao
         CALL qs_upf_inverse_atomic_overlap(basis, atomic%hard, shifted_frame%positions(:, 1), &
                                            shifted_frame%positions(:, 1), cell, images, reduced_k, &
                                            shifted_frame%coefficients(:, 1:2))
         CALL qs_upf_inverse_atomic_overlap(basis, atomic%hard, shifted_frame%positions(:, 1), &
                                            shifted_frame%positions(:, 2), cell, images, reduced_k, &
                                            shifted_frame%coefficients(:, 3:4))
         CALL qs_upf_inverse_atomic_operators(shifted_frame, shifted_frame%q, local_kinetic, local_metric)
         local_samples(:, :, side) = local_kinetic(3:4, 3:4)
         metric_samples(:, :, side) = local_metric(3:4, 3:4)
      END DO
      CALL check(MAXVAL(ABS((local_samples(:, :, 2) - local_samples(:, :, 1))/(2.0_dp*step) - &
                            local_kinetic_response(3:4, 3:4))), 'Relative Gaussian-center kinetic response', 3.e-9_dp)
      CALL check(MAXVAL(ABS((metric_samples(:, :, 2) - metric_samples(:, :, 1))/(2.0_dp*step) - &
                            local_metric_response(3:4, 3:4))), 'Relative Gaussian-center metric response', 3.e-9_dp)
      CALL qs_upf_inverse_atomic_overlap(basis, atomic%hard, center, ao_center, cell, images, -reduced_k, samples(:, :, 1))
      CALL check(MAXVAL(ABS(samples(:, :, 1) - CONJG(overlap))), 'Atomic hard projector time reversal')
      CALL check_nc_atomic_coupling(pot, basis, frame%q, comm)
      CALL deallocate_gto_basis_set(basis)
      pot%pseudo_type = "NC"
      pot%is_ultrasoft = .FALSE.
      CALL upf_inverse_atomic_init(pot, 1, atomic, info)
      IF (info /= 0 .OR. .NOT. atomic%inverse%ready) ERROR STOP 'NC atomic inverse factorization failed'
      ! Compact storage may remove only an exactly zero radial suffix.
      ir = SIZE(atomic%hard%r)
      IF (ir > SIZE(original%r)) ERROR STOP 'NC projector gained radial quadrature points'
      CALL check(MAXVAL(ABS(atomic%hard%r - original%r(:ir))), 'Original NC radial quadrature is retained', 0.0_dp)
      IF (ir < SIZE(original%r)) THEN
         CALL check(MAXVAL(ABS(original%weighted_beta(ir + 1:, :))), 'Omitted NC quadrature rows are exactly zero', 0.0_dp)
      END IF
      CALL check(MAXVAL(ABS(atomic%hard%weighted_beta - original%weighted_beta(:ir, :))) + &
                 MAXVAL(ABS(atomic%hard%d - original%d)) + MAXVAL(ABS(atomic%kinetic)) + MAXVAL(ABS(atomic%metric)), &
                 'Atomic inverse reduces exactly to original NC projectors', 0.0_dp)
   END SUBROUTINE check_atomic_projectors

! **************************************************************************************************
!> \brief Mixed NC projectors retain their overlap with a neighboring periodic reconstruction.
! **************************************************************************************************
   SUBROUTINE check_nc_atomic_coupling(pot, basis, metric, comm)
      TYPE(atom_upfpot_type), INTENT(IN)                 :: pot
      TYPE(gto_basis_set_type), POINTER                  :: basis
      COMPLEX(dp), INTENT(IN)                            :: metric(:, :)
      TYPE(mp_comm_type), INTENT(IN)                     :: comm

      COMPLEX(dp), ALLOCATABLE :: beta(:, :), chi(:, :), correction(:, :), density(:, :), &
         dw(:, :), gradient(:, :, :), op(:, :), original(:, :), samples(:, :, :), soft(:, :)
      INTEGER                                            :: i, images(3, 3), info, ir, j, np, &
                                                            npoint, row, rule, side
      REAL(dp)                                           :: atom_response(3, 2), dc(3, 3), &
                                                            derivative, dk(3), empty_points(3, 0), &
                                                            empty_weights(0), energies(2), k(3), &
                                                            shift
      REAL(dp), ALLOCATABLE                              :: points(:, :), radial_weights(:), &
                                                            weights(:)
      TYPE(atomic_kind_type), TARGET                     :: atomic_kinds(2)
      TYPE(cell_type)                                    :: cell, shifted_cell
      TYPE(particle_type)                                :: atoms(2)
      TYPE(qs_kind_type)                                 :: kinds(2)
      TYPE(qs_upf_inverse_type)                          :: frame, shifted
      TYPE(qs_uspp_type)                                 :: uspp

      np = SIZE(metric, 1)
      DO i = 1, 2
         atomic_kinds(i)%kind_number = i
         atoms(i)%atomic_kind => atomic_kinds(i)
         CALL add_basis_set_to_container(kinds(i)%basis_sets, basis, 'ORB')
         ALLOCATE (kinds(i)%sgp_potential)
         ALLOCATE (kinds(i)%sgp_potential%upf_data, SOURCE=pot)
         CALL upf_projector_init(pot, 3, kinds(i)%sgp_potential%upf_projectors)
      END DO
      ASSOCIATE (nc => kinds(2)%sgp_potential%upf_data)
         nc%pseudo_type = 'NC'
         nc%is_ultrasoft = .FALSE.
         nc%qij = 0.0_dp
         DEALLOCATE (nc%partial_waves)
      END ASSOCIATE
      ALLOCATE (uspp%center_atoms(1), uspp%center_kinds(1), uspp%augmentation(2))
      uspp%center_atoms = 1
      uspp%center_kinds = 1
      ALLOCATE (uspp%augmentation(1)%q, SOURCE=REAL(metric, dp))
      cell%hmat = 0.0_dp
      DO i = 1, 3
         cell%hmat(i, i) = 3.6_dp + 0.2_dp*i
      END DO
      cell%perd = 1
      atoms(1)%r = [0.2_dp, -0.1_dp, 0.3_dp]
      atoms(2)%r = atoms(1)%r + [3.4_dp, 0.3_dp, 0.1_dp]
      images = 0
      images(1, :) = [-1, 0, 1]
      k = [0.23_dp, -0.17_dp, 0.11_dp]
      dk = [0.004_dp, -0.002_dp, 0.003_dp]
      dc = RESHAPE([(0.002_dp*i, i=1, 9)], [3, 3])
      atom_response(:, 1) = [0.01_dp, -0.02_dp, 0.03_dp]
      atom_response(:, 2) = [-0.02_dp, 0.03_dp, 0.01_dp]
      CALL qs_upf_inverse_prepare(kinds, atoms, cell, uspp, empty_points, empty_weights, images, k, &
                                  comm, frame, info, atom_local=.TRUE., store_fields=.FALSE.)
      IF (info /= 0) ERROR STOP 'Mixed NC atomic frame rejected'
      ALLOCATE (original(np, 4), dw(2*np, 4), samples(2*np, 4, 2), density(4, 4), op(4, 4))
      DO i = 1, 2
         CALL qs_upf_inverse_atomic_overlap(basis, frame%atomic_cache(2)%hard, atoms(2)%r, atoms(i)%r, &
                                            cell%hmat, images, k, original(:, 2*i - 1:2*i))
      END DO
      IF (MAXVAL(ABS(original - frame%coefficients(np + 1:, :))) < 1.e-6_dp) ERROR STOP 'NC coupling not exercised'
      rule = get_number_of_lebedev_grid(n=kinds(1)%ngrid_ang)
      npoint = SIZE(pot%r)*SIZE(lebedev_grid(rule)%w)
      ALLOCATE (points(3, npoint), weights(npoint), beta(npoint, np), chi(npoint, 4), soft(npoint, 4), &
                gradient(npoint, 4, 3), correction(np, 4))
      radial_weights = upf_radial_weights(pot)*pot%r**2
      row = 0
      DO ir = 1, SIZE(pot%r)
         DO j = 1, SIZE(lebedev_grid(rule)%w)
            row = row + 1
            points(:, row) = atoms(1)%r + pot%r(ir)*lebedev_grid(rule)%r(:, j)
            weights(row) = radial_weights(ir)*fourpi*lebedev_grid(rule)%w(j)
         END DO
      END DO
      CALL upf_radial_fields_periodic(frame%beta_cache(2), points, atoms(2)%r, cell%hmat, images, k, beta)
      DO i = 1, 2
         CALL gaussian_basis_fields_periodic(basis, points, atoms(i)%r, cell%hmat, images, k, chi(:, 2*i - 1:2*i))
      END DO
      CALL qs_upf_inverse_evaluate(frame, points, soft, gradient)
      CALL upf_inverse_pair(beta, chi - soft, weights, correction)
      CALL check(MAXVAL(ABS(original - frame%coefficients(np + 1:, :) - correction)), &
                 'NC bra agrees with the reconstructed-orbital difference integral')
      CALL qs_upf_inverse_response(frame, empty_points, atom_response, dc, dk, empty_weights, dw)
      DO j = 1, 4
         DO i = 1, 4
            density(i, j) = CMPLX(0.1_dp*COS(REAL(i - j, dp)), 0.05_dp*SIN(REAL(i - j, dp)), dp)
         END DO
         density(j, j) = density(j, j) + 0.8_dp
      END DO
      derivative = 2.0_dp*REAL(SUM(CONJG(dw(np + 1:, :))* &
                                   MATMUL(frame%atomic_cache(2)%hard%d, &
                                          MATMUL(frame%coefficients(np + 1:, :), density))), dp)
      DO side = 1, 2
         shift = REAL(2*side - 3, dp)*step
         DO i = 1, 2
            atoms(i)%r = frame%positions(:, i) + shift*atom_response(:, i)
         END DO
         shifted_cell = cell
         shifted_cell%hmat = cell%hmat + shift*dc
         CALL qs_upf_inverse_prepare(kinds, atoms, shifted_cell, uspp, empty_points, empty_weights, images, k + shift*dk, &
                                     comm, shifted, info, atom_local=.TRUE., store_fields=.FALSE.)
         IF (info /= 0) ERROR STOP 'Displaced mixed NC atomic frame rejected'
         samples(:, :, side) = shifted%coefficients
         op = MATMUL(CONJG(TRANSPOSE(shifted%coefficients(np + 1:, :))), &
                     MATMUL(frame%atomic_cache(2)%hard%d, shifted%coefficients(np + 1:, :)))
         energies(side) = REAL(SUM(density*TRANSPOSE(op)), dp)
      END DO
      CALL check(MAXVAL(ABS((samples(:, :, 2) - samples(:, :, 1))/(2.0_dp*step) - dw)), &
                 'Mixed NC atomic/cell/Bloch coefficient response', 3.e-9_dp)
      CALL check(ABS((energies(2) - energies(1))/(2.0_dp*step) - derivative), &
                 'Both transformed NC nonlocal factors contribute to the geometry derivative', 3.e-9_dp)
      DO i = 1, 2
         atoms(i)%r = frame%positions(:, i)
      END DO
      CALL qs_upf_inverse_prepare(kinds, atoms, cell, uspp, empty_points, empty_weights, images, -k, &
                                  comm, shifted, info, atom_local=.TRUE., store_fields=.FALSE.)
      IF (info /= 0) ERROR STOP 'Time-reversed mixed NC frame rejected'
      CALL check(MAXVAL(ABS(shifted%coefficients - CONJG(frame%coefficients))), 'Mixed NC time reversal')
      kinds(2)%sgp_potential%upf_data = pot
      DEALLOCATE (uspp%center_atoms, uspp%center_kinds)
      ALLOCATE (uspp%center_atoms, SOURCE=[1, 2])
      ALLOCATE (uspp%center_kinds, SOURCE=[1, 2])
      ALLOCATE (uspp%augmentation(2)%q, SOURCE=REAL(metric, dp))
      CALL qs_upf_inverse_prepare(kinds, atoms, cell, uspp, empty_points, empty_weights, images, k, &
                                  comm, shifted, info, atom_local=.TRUE., store_fields=.FALSE.)
      IF (info /= -2) ERROR STOP 'Active reconstruction overlap was accepted'
      DO i = 1, 2
         DEALLOCATE (kinds(i)%sgp_potential)
      END DO
   END SUBROUTINE check_nc_atomic_coupling

! **************************************************************************************************
!> \brief Independent kinetic augmentation, noncompact PS support and its Cartesian response.
! **************************************************************************************************
   SUBROUTINE check_kinetic_fields()
      INTEGER, PARAMETER                                 :: n = 9, nr = 9

      COMPLEX(dp)                                        :: ao(2, 2), ao_direction(2, 2), &
                                                            coefficients(n, 2), dc(n, 2), &
                                                            matrix(2, 2)
      INTEGER                                            :: axis, i, image, j, l, pair, point, side
      REAL(dp) :: atoms_response(3, 1), cell(3, 3), cell_response(3, 3), charge_gradient(3, 3), &
         charge_gradient_response(3, 3), charge_gradient_samples(3, 3, 2), d(3, n), delta, &
         density(3), energy_samples(2), expected, gradient(n, n, 3), gradient_samples(n, n, 3, 2), &
         hard_value(n), hessian(n, n, 3, 3), kernel(n, n), p(3, n), points(3, 3), &
         points_response(3, 3), position(3), potential(3), reference(3), response(3), &
         samples(n, n, 2), shifted(3), soft_value(n), value(n), vg(3, 3), weights(3)
      TYPE(atom_upfpot_type)                             :: pot
      TYPE(qs_upf_inverse_augmentation_type)             :: augmentation
      TYPE(qs_upf_inverse_type)                          :: frame
      TYPE(upf_inverse_kinetic_type)                     :: kinetic, kinetics(2)
      TYPE(upf_projector_type)                           :: projectors

      pot%number_of_proj = 3
      ALLOCATE (pot%r(nr), pot%rab(nr), pot%lbeta(3), pot%partial_waves)
      ALLOCATE (pot%partial_waves%ae_wfc(nr, 3), pot%partial_waves%ps_wfc(nr, 3), &
                pot%partial_waves%angular_momenta(3))
      pot%lbeta(:) = [0, 1, 2]
      pot%partial_waves%angular_momenta(:) = pot%lbeta
      pot%r(:) = [(0.25_dp*REAL(i - 1, dp), i=1, nr)]
      pot%rab = 0.25_dp
      DO j = 1, 3
         DO i = 1, nr
            pot%partial_waves%ps_wfc(i, j) = pot%r(i)**j*(1.0_dp + 0.1_dp*pot%r(i)**2 - 0.01_dp*pot%r(i)**4)
            delta = 0.0_dp
            IF (pot%r(i) < 1.5_dp) delta = 0.2_dp*pot%r(i)**j*(1.0_dp - (pot%r(i)/1.5_dp)**2)**3
            pot%partial_waves%ae_wfc(i, j) = pot%partial_waves%ps_wfc(i, j) + delta
         END DO
      END DO
      CALL upf_inverse_kinetic_init(pot, kinetic)
      CALL upf_radial_fields_evaluate(kinetic%pseudo, [2.0_dp, 0.0_dp, 0.0_dp], value)
      CALL check(ABS(value(1) - 1.24_dp/SQRT(4.0_dp*pi)), 'Noncompact original PS endpoint retained')
      DO point = 1, 4
         SELECT CASE (point)
         CASE (1)
            position = 0.0_dp
         CASE (2)
            position = [0.4_dp, 0.0_dp, 0.0_dp]
         CASE (3)
            position = [0.3_dp, -0.2_dp, 0.7_dp]
         CASE (4)
            position = [1.45_dp, 0.0_dp, 0.0_dp]
         END SELECT
         CALL upf_inverse_kinetic_evaluate(kinetic, position, kernel, gradient)
         CALL upf_radial_fields_evaluate(kinetic%pseudo, position, value, p)
         CALL upf_radial_fields_evaluate(kinetic%delta, position, value, d)
         DO j = 1, n
            DO i = 1, n
               delta = 0.5_dp*(DOT_PRODUCT(p(:, i) + d(:, i), p(:, j) + d(:, j)) - DOT_PRODUCT(p(:, i), p(:, j)))
               CALL check(ABS(kernel(i, j) - delta), 'Original AE minus PS kinetic products')
            END DO
         END DO
         DO axis = 1, 3
            DO side = 1, 2
               position(axis) = position(axis) + REAL(2*side - 3, dp)*step
               CALL upf_inverse_kinetic_evaluate(kinetic, position, samples(:, :, side))
               position(axis) = position(axis) - REAL(2*side - 3, dp)*step
            END DO
            CALL check(MAXVAL(ABS((samples(:, :, 2) - samples(:, :, 1))/(2.0_dp*step) - gradient(:, :, axis))), &
                       'Independent kinetic augmentation Cartesian derivative', 3.e-8_dp)
         END DO
         CALL upf_inverse_density_evaluate(kinetic, position, kernel, gradient, hessian)
         CALL upf_radial_fields_evaluate(kinetic%pseudo, position, soft_value)
         CALL upf_radial_fields_evaluate(kinetic%delta, position, hard_value)
         hard_value = hard_value + soft_value
         DO j = 1, n
            DO i = 1, n
               CALL check(ABS(kernel(i, j) - hard_value(i)*hard_value(j) + soft_value(i)*soft_value(j)), &
                          'Hard minus soft density product kernel')
            END DO
         END DO
         DO axis = 1, 3
            DO side = 1, 2
               shifted = position
               shifted(axis) = shifted(axis) + REAL(2*side - 3, dp)*step
               CALL upf_inverse_density_evaluate(kinetic, shifted, samples(:, :, side), gradient_samples(:, :, :, side))
            END DO
            CALL check(MAXVAL(ABS((samples(:, :, 2) - samples(:, :, 1))/(2.0_dp*step) - gradient(:, :, axis))), &
                       'Hard minus soft density kernel Cartesian derivative', 3.e-8_dp)
            CALL check(MAXVAL(ABS((gradient_samples(:, :, :, 2) - gradient_samples(:, :, :, 1))/(2.0_dp*step) - &
                                  hessian(:, :, :, axis))), 'Hard minus soft density kernel Cartesian Hessian', 3.e-8_dp)
         END DO
         CALL upf_inverse_density_evaluate(kinetic, position, samples(:, :, 1), hessian=hessian)
         CALL check(MAXVAL(ABS(samples(:, :, 1) - kernel)), 'Density kernel with Hessian only')
      END DO
      CALL upf_inverse_density_evaluate(kinetic, [4.0_dp, 0.0_dp, 0.0_dp], kernel, gradient, hessian)
      CALL check(MAXVAL(ABS(kernel)) + MAXVAL(ABS(gradient)) + MAXVAL(ABS(hessian)), &
                 'Compact density difference outside PS table', 0.0_dp)
      CALL upf_inverse_kinetic_evaluate(kinetic, [4.0_dp, 0.0_dp, 0.0_dp], kernel, gradient)
      CALL check(MAXVAL(ABS(kernel)) + MAXVAL(ABS(gradient)), 'Compact kinetic difference outside PS table', 0.0_dp)
      ALLOCATE (frame%coefficients(n, 2), frame%positions(3, 1), &
                frame%center_atoms(1), frame%projector_first(1), frame%projector_last(1))
      frame%center_atoms = 1
      frame%projector_first = 1
      frame%projector_last = n
      frame%positions = 0.0_dp
      DO j = 1, 2
         DO i = 1, n
            coefficients(i, j) = CMPLX(0.02_dp*REAL(i + j, dp), -0.01_dp*REAL(i - j, dp), dp)
            dc(i, j) = CMPLX(0.003_dp*REAL(i + 2*j, dp), 0.002_dp*REAL(2*i - j, dp), dp)
         END DO
      END DO
      frame%coefficients(:, :) = coefficients
      ao = RESHAPE([CMPLX(1.0_dp, 0.0_dp, dp), CMPLX(0.02_dp, -0.03_dp, dp), &
                    CMPLX(0.02_dp, 0.03_dp, dp), CMPLX(0.7_dp, 0.0_dp, dp)], [2, 2])
      ao_direction = RESHAPE([CMPLX(0.01_dp, 0.0_dp, dp), CMPLX(-0.02_dp, 0.03_dp, dp), &
                              CMPLX(-0.02_dp, -0.03_dp, dp), CMPLX(-0.04_dp, 0.0_dp, dp)], [2, 2])
      points(:, 1) = [0.2_dp, -0.1_dp, 0.3_dp]
      points(:, 2) = [0.5_dp, 0.2_dp, -0.4_dp]
      points(:, 3) = [1.3_dp, 0.0_dp, 0.0_dp]
      weights = [0.3_dp, 0.4_dp, 0.2_dp]
      potential = [0.7_dp, -0.3_dp, 1.1_dp]
      CALL qs_upf_inverse_kinetic_center(frame, kinetic, 1, points, ao, density)
      CALL qs_upf_inverse_kinetic_integrate(frame, kinetic, 1, points, weights, potential, matrix)
      CALL check(ABS(SUM(weights*potential*density) - REAL(SUM(ao*CONJG(matrix)), dp)), 'Inverse kinetic AO adjoint')
      expected = REAL(SUM(ao_direction*CONJG(matrix)), dp)
      DO side = 1, 2
         delta = REAL(2*side - 3, dp)*step
         CALL qs_upf_inverse_kinetic_center(frame, kinetic, 1, points, ao + delta*ao_direction, shifted)
         energy_samples(side) = SUM(weights*potential*shifted)
      END DO
      CALL check(ABS((energy_samples(2) - energy_samples(1))/(2.0_dp*step) - expected), 'Inverse kinetic AO variation', 3.e-9_dp)
      points_response = 0.02_dp - 0.01_dp*points
      atoms_response(:, 1) = [0.01_dp, -0.02_dp, 0.03_dp]
      CALL qs_upf_inverse_kinetic_center(frame, kinetic, 1, points, ao, density, response, &
                                         points_response, atoms_response, dc)
      DO side = 1, 2
         delta = REAL(2*side - 3, dp)*step
         frame%coefficients(:, :) = coefficients + delta*dc
         frame%positions(:, :) = delta*atoms_response
         CALL qs_upf_inverse_kinetic_center(frame, kinetic, 1, points + delta*points_response, ao, shifted)
         energy_samples(side) = SUM(weights*potential*shifted)
      END DO
      CALL check(ABS((energy_samples(2) - energy_samples(1))/(2.0_dp*step) - SUM(weights*potential*response)), &
                 'Inverse kinetic center point and coefficient response', 3.e-8_dp)
      CALL qs_upf_inverse_kinetic_integrate(frame, kinetic, 1, points(:, :0), weights(:0), potential(:0), matrix)
      CALL check(MAXVAL(ABS(matrix)), 'Empty kinetic point batch', 0.0_dp)
      ALLOCATE (frame%center_kinds(1), frame%images(3, 3))
      frame%center_kinds = 1
      frame%images = 0
      frame%images(1, :) = [-1, 0, 1]
      frame%inverse%ready = .TRUE.
      cell = 0.0_dp
      DO axis = 1, 3
         cell(axis, axis) = 0.9_dp
      END DO
      cell_response = RESHAPE([(0.001_dp*REAL(i, dp), i=1, 9)], [3, 3])
      frame%cell = cell
      frame%coefficients(:, :) = coefficients
      frame%positions = 0.0_dp
      kinetics(1) = kinetic
      CALL qs_upf_inverse_kinetic_collocate(frame, kinetics, points, ao, density, response, &
                                            points_response, atoms_response, cell_response, dc)
      CALL qs_upf_inverse_kinetic_operator(frame, kinetics, points, weights, potential, matrix)
      CALL check(ABS(SUM(weights*potential*density) - REAL(SUM(ao*CONJG(matrix)), dp)), 'Periodic kinetic AO adjoint')
      reference = 0.0_dp
      DO image = 1, 3
         position = MATMUL(cell, REAL(frame%images(:, image), dp))
         CALL qs_upf_inverse_kinetic_center(frame, kinetic, 1, points - SPREAD(position, 2, 3), ao, shifted)
         reference = reference + shifted
      END DO
      CALL check(MAXVAL(ABS(density - reference)), 'Periodic kinetic home image sum')
      DO side = 1, 2
         delta = REAL(2*side - 3, dp)*step
         frame%coefficients(:, :) = coefficients + delta*dc
         frame%positions(:, :) = delta*atoms_response
         frame%cell = cell + delta*cell_response
         CALL qs_upf_inverse_kinetic_collocate(frame, kinetics, points + delta*points_response, ao, shifted)
         energy_samples(side) = SUM(weights*potential*shifted)
      END DO
      CALL check(ABS((energy_samples(2) - energy_samples(1))/(2.0_dp*step) - SUM(weights*potential*response)), &
                 'Periodic kinetic image cell and coefficient response', 3.e-8_dp)
      frame%center_kinds = 2
      CALL qs_upf_inverse_kinetic_collocate(frame, kinetics, points, ao, density)
      CALL qs_upf_inverse_kinetic_operator(frame, kinetics, points, weights, potential, matrix)
      CALL check(MAXVAL(ABS(density)) + MAXVAL(ABS(matrix)), 'Zero NC kinetic augmentation', 0.0_dp)
      ! Independent compensation functions deliberately differ from the reference-wave products.
      ALLOCATE (pot%beta(nr, 3), pot%dion(3, 3), pot%qijl(nr, 6, 0:4), augmentation%kinds(2))
      pot%beta = 0.0_dp
      pot%dion = 0.0_dp
      pot%qijl = 0.0_dp
      DO i = 1, nr
         IF (pot%r(i) >= 1.75_dp) CYCLE
         DO l = 0, 4
            DO pair = 1, 6
               pot%qijl(i, pair, l) = 0.01_dp*REAL(pair, dp)*pot%r(i)**(l + 2)*(1 - (pot%r(i)/1.75_dp)**2)**3
            END DO
         END DO
      END DO
      CALL upf_projector_init(pot, 2, projectors)
      CALL uspp_augmentation_init(pot, projectors, augmentation%kinds(1))
      CALL uspp_augmentation_prepare_realspace(pot, augmentation%kinds(1))
      frame%center_kinds = 1
      frame%coefficients(:, :) = coefficients
      frame%positions = 0.0_dp
      frame%cell = cell
      vg = RESHAPE([(0.01_dp*REAL(i - 5, dp), i=1, 9)], [3, 3])
      CALL qs_upf_inverse_paw_collocate(frame, kinetics, augmentation, points, ao, density, charge_gradient, &
                                        response, charge_gradient_response, points_response, atoms_response, cell_response, dc)
      CALL qs_upf_inverse_paw_operator(frame, kinetics, augmentation, points, weights, potential, vg, matrix)
      CALL check(ABS(SUM(weights*(potential*density + SUM(vg*charge_gradient, DIM=2))) - REAL(SUM(ao*CONJG(matrix)), dp)), &
                 'PAW reconstructed charge and gradient AO adjoint')
      DO side = 1, 2
         delta = REAL(2*side - 3, dp)*step
         frame%coefficients(:, :) = coefficients + delta*dc
         frame%positions(:, :) = delta*atoms_response
         frame%cell = cell + delta*cell_response
         CALL qs_upf_inverse_paw_collocate(frame, kinetics, augmentation, points + delta*points_response, ao, &
                                           shifted, charge_gradient_samples(:, :, side))
         energy_samples(side) = SUM(weights*(potential*shifted + SUM(vg*charge_gradient_samples(:, :, side), DIM=2)))
      END DO
      CALL check(MAXVAL(ABS((charge_gradient_samples(:, :, 2) - charge_gradient_samples(:, :, 1))/(2*step) - &
                            charge_gradient_response)), 'PAW reconstructed density gradient geometry response', 3.e-8_dp)
      CALL check(ABS((energy_samples(2) - energy_samples(1))/(2*step) - &
                     SUM(weights*(potential*response + SUM(vg*charge_gradient_response, DIM=2)))), &
                 'PAW reconstructed density periodic image and inverse response', 3.e-8_dp)
      CALL qs_upf_inverse_paw_operator(frame, kinetics, augmentation, points(:, :0), weights(:0), potential(:0), vg(:0, :), matrix)
      CALL check(MAXVAL(ABS(matrix)), 'Empty PAW reconstruction adjoint', 0.0_dp)
      frame%center_kinds = 2
      CALL qs_upf_inverse_paw_collocate(frame, kinetics, augmentation, points, ao, density, charge_gradient)
      CALL qs_upf_inverse_paw_operator(frame, kinetics, augmentation, points, weights, potential, vg, matrix)
      CALL check(MAXVAL(ABS(density)) + MAXVAL(ABS(charge_gradient)) + MAXVAL(ABS(matrix)), &
                 'Inactive PAW kind has no reconstruction', 0.0_dp)
   END SUBROUTINE check_kinetic_fields

! **************************************************************************************************
!> \brief Native hard/soft and inverse angular legs with a deliberately nonidentity Gram matrix.
! **************************************************************************************************
   SUBROUTINE check_center_angular()
      INTEGER                                            :: ir, side
      REAL(dp) :: adjoint(6), dh(2, 2), dpnt(6), ds(2, 2), exact, hard(2, 2), hard_adjoint(2, 2), &
         potential(2, 2), radial(2, 2), response, samples(2), shift, soft(2, 2), &
         soft_adjoint(2, 2), values(6), wr(2)
      TYPE(grid_atom_type)                               :: grid
      TYPE(harmonics_atom_type)                          :: harmonics

      grid%nr = 2
      grid%ng_sphere = 3
      ALLOCATE (grid%wa(3), harmonics%slm(3, 2))
      grid%wa = [1.1_dp, 0.7_dp, 1.4_dp]
      harmonics%slm(:, 1) = 1.0_dp/SQRT(4.0_dp*pi)
      harmonics%slm(:, 2) = [0.3_dp, -0.4_dp, 0.2_dp]
      values = [0.03_dp, -0.01_dp, 0.02_dp, 0.04_dp, 0.01_dp, -0.02_dp]
      dpnt = [0.011_dp, -0.003_dp, 0.007_dp, -0.004_dp, 0.009_dp, 0.005_dp]
      hard(:, 1) = [0.14_dp, 0.09_dp]
      hard(:, 2) = [-0.03_dp, 0.02_dp]
      soft = 0.3_dp*hard
      dh = 0.017_dp*hard
      ds = -0.023_dp*soft
      wr = [0.21_dp, 0.37_dp]
      potential(:, 1) = [0.08_dp, -0.03_dp]
      potential(:, 2) = [0.04_dp, 0.07_dp]
      CALL qs_upf_inverse_density_radial(grid, harmonics, values, radial, hard, soft)
      CALL qs_upf_inverse_density_radial_adjoint(grid, harmonics, potential, adjoint, hard_adjoint, soft_adjoint)
      exact = 0.0_dp
      response = 0.0_dp
      DO ir = 1, 2
         exact = exact + wr(ir)*SUM(grid%wa*adjoint((ir - 1)*3 + 1:ir*3)*values((ir - 1)*3 + 1:ir*3)) + &
                 wr(ir)*SUM(hard_adjoint(ir, :)*hard(ir, :) + soft_adjoint(ir, :)*soft(ir, :))
         response = response + wr(ir)*SUM(grid%wa*adjoint((ir - 1)*3 + 1:ir*3)*dpnt((ir - 1)*3 + 1:ir*3)) + &
                    wr(ir)*SUM(hard_adjoint(ir, :)*dh(ir, :) + soft_adjoint(ir, :)*ds(ir, :))
      END DO
      CALL check(ABS(SUM(SPREAD(wr, 2, 2)*radial*potential) - exact), 'complete finite angular adjoint', 2.e-14_dp)
      DO side = 1, 2
         shift = REAL(2*side - 3, dp)*step
         CALL qs_upf_inverse_density_radial(grid, harmonics, values + shift*dpnt, radial, &
                                            hard + shift*dh, soft + shift*ds)
         samples(side) = SUM(SPREAD(wr, 2, 2)*radial*potential)
      END DO
      CALL check(ABS((samples(2) - samples(1))/(2.0_dp*step) - response), 'all angular density legs', 2.e-11_dp)
      CALL qs_upf_inverse_density_radial(grid, harmonics, values, radial)
      CALL check(ABS(SUM(SPREAD(wr, 2, 2)*radial*potential) - &
                     SUM(wr*[SUM(grid%wa*adjoint(1:3)*values(1:3)), SUM(grid%wa*adjoint(4:6)*values(4:6))])), &
                 'inverse-only angular adjoint', 2.e-14_dp)
      DEALLOCATE (grid%wa, harmonics%slm)
   END SUBROUTINE check_center_angular

! **************************************************************************************************
!> \brief Residual image sums against the shell theorem, density adjoints and cell response.
! **************************************************************************************************
   SUBROUTINE check_periodic_coulomb()
      COMPLEX(dp)                                        :: adjoint(3), background(3), &
                                                            cached_adjoint(3), db(3), &
                                                            density_g(3), exact_density(3)
      INTEGER                                            :: a, b, boundary, i, ix, iy, iz, j, &
                                                            periodicity(3), side
      REAL(dp) :: cached_energy, cached_response, cell(3, 3), dcell(3, 3), dg(3, 3), &
         displacement(3), distance, dkw(3), dvolume, energy, exact, g(3, 3), kernel, kw(3), norm, &
         qa(2), qb(2), response, samples(2), shift, trace, volume
      TYPE(upf_inverse_coulomb_center_type)              :: centers(2), changed(2), derivative(2)
      TYPE(upf_inverse_coulomb_map_type), ALLOCATABLE    :: map
      TYPE(upf_inverse_coulomb_potential_type), &
         ALLOCATABLE                                     :: cached_potential(:), potential(:)

      cell(:, 1) = [2.0_dp, 0.0_dp, 0.0_dp]
      cell(:, 2) = [0.25_dp, 1.9_dp, 0.0_dp]
      cell(:, 3) = [0.1_dp, 0.2_dp, 2.3_dp]
      dcell = 0.013_dp*cell
      volume = 2.0_dp*1.9_dp*2.3_dp
      dvolume = 0.039_dp*volume
      g(:, 1) = 0.0_dp
      g(:, 2) = pi*[1.0_dp, -0.25_dp/1.9_dp, (0.25_dp*0.2_dp/1.9_dp - 0.1_dp)/2.3_dp]
      g(:, 3) = -g(:, 2)
      dg = -0.013_dp*g
      background = [CMPLX(0.0_dp, 0.0_dp, dp), CMPLX(0.012_dp, -0.007_dp, dp), CMPLX(0.012_dp, 0.007_dp, dp)]
      db = 0.021_dp*background
      kw = 0.0_dp
      DO a = 1, 2
         ALLOCATE (centers(a)%l(1), centers(a)%indices(1), centers(a)%radii(2), centers(a)%weights(2), &
                   centers(a)%density(2, 1), centers(a)%shapes(2, 1))
         centers(a)%l(:) = [0]
         centers(a)%indices(:) = [1]
         centers(a)%radii(:) = [0.5_dp, 1.1_dp] + REAL(a - 1, dp)*0.1_dp
         centers(a)%weights(:) = [0.2_dp, 0.4_dp]
         centers(a)%density(:, :) = RESHAPE([0.7_dp, -0.15_dp]/REAL(a, dp), [2, 1])
         centers(a)%shapes(:, :) = RESHAPE([0.8_dp, 0.2_dp], [2, 1])
         centers(a)%position = REAL(a - 1, dp)*[0.7_dp, 0.4_dp, -0.3_dp]
         derivative(a) = centers(a)
         derivative(a)%radii(:) = [0.011_dp, -0.007_dp]
         derivative(a)%weights(:) = [0.003_dp, -0.001_dp]
         derivative(a)%density(:, 1) = REAL(a, dp)*[0.023_dp, -0.017_dp]
         derivative(a)%shapes(:, :) = 0.017_dp*centers(a)%shapes
         derivative(a)%position = REAL(a, dp)*[0.009_dp, -0.011_dp, 0.013_dp]
      END DO
      ALLOCATE (map)
      DO boundary = 1, 4
         periodicity = 0
         IF (boundary > 1) periodicity(:boundary - 1) = 1
         kw = 0.0_dp
         CALL upf_inverse_coulomb_periodic(centers, cell, g, volume, background, kw, energy, potential, adjoint, &
                                           periodicity=periodicity)
         exact = 0.0_dp
         DO a = 1, 2
            norm = SUM(centers(a)%weights*centers(a)%shapes(:, 1))
            qa = centers(a)%weights*(centers(a)%density(:, 1) - centers(a)%shapes(:, 1)* &
                                     SUM(centers(a)%weights*centers(a)%density(:, 1))/norm)
            DO b = 1, 2
               norm = SUM(centers(b)%weights*centers(b)%shapes(:, 1))
               qb = centers(b)%weights*(centers(b)%density(:, 1) - centers(b)%shapes(:, 1)* &
                                        SUM(centers(b)%weights*centers(b)%density(:, 1))/norm)
               DO iz = -3, 3
                  DO iy = -3, 3
                     DO ix = -3, 3
                        IF (ANY([ix, iy, iz]*(1 - periodicity) /= 0)) CYCLE
                        displacement = centers(a)%position - centers(b)%position - MATMUL(cell, REAL([ix, iy, iz], dp))
                        distance = SQRT(SUM(displacement**2))
                        DO i = 1, 2
                           DO j = 1, 2
                              IF (distance <= ABS(centers(a)%radii(i) - centers(b)%radii(j))) THEN
                                 kernel = 4.0_dp*pi/MAX(centers(a)%radii(i), centers(b)%radii(j))
                              ELSE IF (distance >= centers(a)%radii(i) + centers(b)%radii(j)) THEN
                                 kernel = 4.0_dp*pi/distance
                              ELSE
                                 kernel = pi/(distance*centers(a)%radii(i)*centers(b)%radii(j))* &
                                          (2.0_dp*distance*(centers(a)%radii(i) + centers(b)%radii(j)) - distance**2 - &
                                           (centers(a)%radii(i) - centers(b)%radii(j))**2)
                              END IF
                              exact = exact + 0.5_dp*qa(i)*qb(j)*kernel
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
         CALL check(ABS(energy - exact), "Complete periodic residual shell sum", 3.e-12_dp)
         kw(2:3) = volume*4.0_dp*pi/SUM(g(:, 2)**2)
         ! A wire kernel can have negative Fourier coefficients. Its adjoints remain unchanged.
         IF (boundary == 2) kw = -kw
         dkw = 0.065_dp*kw
         CALL upf_inverse_coulomb_periodic(centers, cell, g, volume, background, kw, energy, potential, adjoint, &
                                           response=response, center_response=derivative, cell_response=dcell, &
                                           gvecs_response=dg, volume_response=dvolume, background_response=db, &
                                           poisson_weights_response=dkw, periodicity=periodicity)
         CALL upf_inverse_coulomb_prepare(centers, cell, g, volume, map, periodicity)
         CALL upf_inverse_coulomb_density(centers, cell, g, volume, density_g, map)
         exact_density = CMPLX(0.0_dp, 0.0_dp, dp)
         DO a = 1, SIZE(centers)
            exact_density(1) = exact_density(1) + SQRT(fourpi)/volume* &
                               SUM(centers(a)%weights*centers(a)%density(:, 1))
            DO i = 2, SIZE(g, 2)
               distance = SQRT(SUM(g(:, i)**2))
               exact_density(i) = exact_density(i) + SQRT(fourpi)/volume* &
                                  EXP(CMPLX(0.0_dp, -DOT_PRODUCT(g(:, i), centers(a)%position), dp))* &
                                  SUM(centers(a)%weights*centers(a)%density(:, 1)* &
                                      SIN(distance*centers(a)%radii)/(distance*centers(a)%radii))
            END DO
         END DO
         CALL check(MAXVAL(ABS(density_g - exact_density)), "Complete signed density shell form factors", 3.e-15_dp)
         CALL check(ABS(density_g(2) - CONJG(density_g(3))), "Complete signed density Fourier conjugacy", 3.e-15_dp)
         CALL upf_inverse_coulomb_periodic(centers, cell, g, volume, background, kw, cached_energy, &
                                           cached_potential, cached_adjoint, map=map)
         CALL check(ABS(cached_energy - energy), "Cached periodic Coulomb energy", 3.e-12_dp)
         CALL check(MAXVAL(ABS(cached_adjoint - adjoint)), "Cached periodic Coulomb background adjoint", 3.e-12_dp)
         DO a = 1, SIZE(centers)
            CALL check(MAXVAL(ABS(cached_potential(a)%values - potential(a)%values)), &
                       "Cached periodic Coulomb center adjoint", 3.e-12_dp)
         END DO
         CALL upf_inverse_coulomb_periodic(centers, cell, g, volume, background, kw, cached_energy, &
                                           cached_potential, cached_adjoint, response=cached_response, &
                                           center_response=derivative, cell_response=dcell, gvecs_response=dg, &
                                           volume_response=dvolume, background_response=db, poisson_weights_response=dkw, map=map)
         CALL check(ABS(cached_response - response), "Cached periodic Coulomb joint response", 3.e-12_dp)
         trace = REAL(SUM(CONJG(background)*adjoint), dp)
         DO a = 1, 2
            trace = trace + SUM(SPREAD(centers(a)%weights, 2, 1)*centers(a)%density*potential(a)%values)
         END DO
         CALL check(ABS(trace - 2.0_dp*energy), "Periodic Coulomb density and background adjoints", 3.e-12_dp)
         trace = REAL(SUM(CONJG(db)*adjoint), dp)
         DO a = 1, 2
            trace = trace + SUM(SPREAD(centers(a)%weights, 2, 1)*derivative(a)%density*potential(a)%values)
         END DO
         changed = derivative
         DO a = 1, SIZE(changed)
            changed(a)%position = 0.0_dp
            changed(a)%radii = 0.0_dp
            changed(a)%weights = 0.0_dp
            changed(a)%shapes = 0.0_dp
         END DO
         CALL upf_inverse_coulomb_periodic(centers, cell, g, volume, background, kw, cached_energy, &
                                           cached_potential, cached_adjoint, response=cached_response, &
                                           center_response=changed, cell_response=0.0_dp*dcell, gvecs_response=0.0_dp*dg, &
                                           volume_response=0.0_dp, background_response=db, poisson_weights_response=0.0_dp*dkw, &
                                           map=map)
         CALL check(ABS(cached_response - trace), "Cached fixed-shell density response", 3.e-12_dp)
         DO side = 1, 2
            changed = derivative
            DO a = 1, SIZE(changed)
               changed(a)%radii = 0.0_dp
            END DO
            IF (side == 1) changed(2)%position = changed(1)%position
            CALL upf_inverse_coulomb_periodic(centers, cell, g, volume, background, kw, exact, &
                                              potential, adjoint, response=trace, center_response=changed, &
                                              cell_response=REAL(side - 1, dp)*dcell, gvecs_response=dg, &
                                              volume_response=dvolume, background_response=db, &
                                              poisson_weights_response=dkw, periodicity=periodicity)
            CALL upf_inverse_coulomb_periodic(centers, cell, g, volume, background, kw, cached_energy, &
                                              cached_potential, cached_adjoint, response=cached_response, &
                                              center_response=changed, cell_response=REAL(side - 1, dp)*dcell, &
                                              gvecs_response=dg, volume_response=dvolume, background_response=db, &
                                              poisson_weights_response=dkw, map=map)
            CALL check(ABS(cached_energy - exact), "Cached moving-center Coulomb energy", 3.e-12_dp)
            CALL check(ABS(cached_response - trace), "Cached moving-center Coulomb response", 3.e-12_dp)
            DO a = 1, SIZE(centers)
               CALL check(MAXVAL(ABS(cached_potential(a)%values - potential(a)%values)), &
                          "Cached moving-center Coulomb adjoint", 3.e-12_dp)
            END DO
         END DO
         trace = REAL(SUM(CONJG(db)*adjoint), dp)
         DO a = 1, 2
            trace = trace + SUM(SPREAD(centers(a)%weights, 2, 1)*derivative(a)%density*potential(a)%values)
         END DO
         DO side = 1, 2
            shift = REAL(2*side - 3, dp)*step
            changed = centers
            DO a = 1, 2
               changed(a)%density(:, :) = centers(a)%density + shift*derivative(a)%density
            END DO
            CALL upf_inverse_coulomb_periodic(changed, cell, g, volume, background + shift*db, kw, &
                                              samples(side), potential, adjoint, periodicity=periodicity)
            CALL upf_inverse_coulomb_periodic(changed, cell, g, volume, background + shift*db, kw, &
                                              cached_energy, cached_potential, cached_adjoint, map=map)
            CALL check(ABS(cached_energy - samples(side)), "Cached periodic Coulomb density update", 3.e-12_dp)
         END DO
         CALL check(ABS((samples(2) - samples(1))/(2.0_dp*step) - trace), "Periodic Coulomb density variation", 3.e-10_dp)
         DO side = 1, 2
            shift = REAL(2*side - 3, dp)*step
            changed = centers
            DO a = 1, 2
               changed(a)%position = centers(a)%position + shift*derivative(a)%position
               changed(a)%radii(:) = centers(a)%radii + shift*derivative(a)%radii
               changed(a)%weights(:) = centers(a)%weights + shift*derivative(a)%weights
               changed(a)%density(:, :) = centers(a)%density + shift*derivative(a)%density
               changed(a)%shapes(:, :) = centers(a)%shapes + shift*derivative(a)%shapes
            END DO
            CALL upf_inverse_coulomb_periodic(changed, cell + shift*dcell, g + shift*dg, volume + shift*dvolume, &
                                              background + shift*db, kw + shift*dkw, samples(side), potential, adjoint, &
                                              periodicity=periodicity)
         END DO
         CALL check(ABS((samples(2) - samples(1))/(2.0_dp*step) - response), "Complete periodic Coulomb joint response", 3.e-9_dp)
      END DO
   END SUBROUTINE check_periodic_coulomb
! **************************************************************************************************
!> \brief Translated shell theorem, overlap, reciprocity and joint geometry response.
! **************************************************************************************************
   SUBROUTINE check_pair_coulomb()
      INTEGER, PARAMETER                                 :: indices(4) = [1, 2, 3, 4], &
                                                            ls(4) = [0, 1, 1, 1]

      INTEGER                                            :: i, j, side
      REAL(dp) :: a(2), b(2), d(3), da(2), db(2), dd(3), distance, exact, kernel(2, 4, 2, 4), &
         response(2, 4, 2, 4), reverse(2, 4, 2, 4), samples(2, 4, 2, 4, 2), shift

      a = [0.7_dp, 1.3_dp]
      b = [0.5_dp, 1.1_dp]
      da = [0.013_dp, -0.007_dp]
      db = [-0.011_dp, 0.017_dp]
      dd = [0.019_dp, -0.023_dp, 0.009_dp]
      DO j = 1, 3
         d = REAL(j, dp)*[0.6_dp, -0.4_dp, 0.3_dp]
         distance = SQRT(SUM(d**2))
         CALL upf_inverse_coulomb_pair_kernel(a, ls, indices, b, ls, indices, d, &
                                              kernel, response, da, db, dd)
         CALL upf_inverse_coulomb_pair_kernel(b, ls, indices, a, ls, indices, -d, reverse)
         DO i = 1, 2
            DO side = 1, 2
               IF (distance >= a(i) + b(side)) THEN
                  exact = 4.0_dp*pi/distance
               ELSE IF (distance <= ABS(a(i) - b(side))) THEN
                  exact = 4.0_dp*pi/MAX(a(i), b(side))
               ELSE
                  exact = pi/(distance*a(i)*b(side))* &
                          (2.0_dp*distance*(a(i) + b(side)) - distance**2 - (a(i) - b(side))**2)
               END IF
               CALL check(ABS(kernel(i, 1, side, 1) - exact), "Translated monopole shells", 3.e-12_dp)
               CALL check(MAXVAL(ABS(kernel(i, :, side, :) - TRANSPOSE(reverse(side, :, i, :)))), &
                          "Translated Coulomb reciprocity", 1.e-12_dp)
            END DO
         END DO
         DO side = 1, 2
            shift = REAL(2*side - 3, dp)*step
            CALL upf_inverse_coulomb_pair_kernel(a + shift*da, ls, indices, b + shift*db, ls, indices, &
                                                 d + shift*dd, samples(:, :, :, :, side))
         END DO
         CALL check(MAXVAL(ABS((samples(:, :, :, :, 2) - samples(:, :, :, :, 1))/(2.0_dp*step) - response)), &
                    "Translated Coulomb joint response", 3.e-9_dp)
      END DO
      ! A source sphere passing through a target angular node is integrated in separate polar caps.
      d = [0.7_dp, 0.0_dp, 0.0_dp]
      CALL upf_inverse_coulomb_pair_kernel(a, ls, indices, b, ls, indices, d, &
                                           kernel, response, da, db, dd)
      DO side = 1, 2
         shift = REAL(2*side - 3, dp)*step
         CALL upf_inverse_coulomb_pair_kernel(a + shift*da, ls, indices, b + shift*db, ls, indices, &
                                              d + shift*dd, samples(:, :, :, :, side))
      END DO
      CALL check(MAXVAL(ABS((samples(:, :, :, :, 2) - samples(:, :, :, :, 1))/(2.0_dp*step) - response)), &
                 "Translated shell crossing response", 3.e-9_dp)
      d = [0.0_dp, 0.0_dp, a(2) + b(2)]
      CALL upf_inverse_coulomb_pair_kernel(a, ls, indices, b, ls, indices, d, &
                                           kernel, response, da, db, dd)
      DO side = 1, 2
         shift = REAL(2*side - 3, dp)*step
         CALL upf_inverse_coulomb_pair_kernel(a + shift*da, ls, indices, b + shift*db, ls, indices, &
                                              d + shift*dd, samples(:, :, :, :, side))
      END DO
      CALL check(MAXVAL(ABS((samples(:, :, :, :, 2) - samples(:, :, :, :, 1))/(2.0_dp*step) - response)), &
                 "Translated tangent shell response", 3.e-9_dp)
      d = 0.0_dp
      CALL upf_inverse_coulomb_pair_kernel(a, ls, indices, b, ls, indices, d, kernel)
      CALL check(ABS(kernel(1, 2, 1, 2) - 4.0_dp*pi/3.0_dp*b(1)/a(1)**2), "Concentric dipole shells")
      CALL check(ABS(kernel(1, 2, 1, 3)), "Concentric harmonic orthogonality")
   END SUBROUTINE check_pair_coulomb
! **************************************************************************************************
!> \brief Original core normalization, finite endpoint and joint phase/lattice/volume derivative.
! **************************************************************************************************
   SUBROUTINE check_core_fourier()
      COMPLEX(dp)                                        :: dg(3, 3), dr(3, 3), response(3), rho(3), &
                                                            samples(3, 2)
      INTEGER                                            :: side
      REAL(dp)                                           :: delta_g(3, 3), delta_r(3), g(3, 3), &
                                                            position(3), shift
      TYPE(atom_upfpot_type)                             :: pot

      ALLOCATE (pot%r(5), pot%rab(5), pot%rho_nlcc(5), pot%tau_nlcc(5))
      pot%r(:) = [0.0_dp, 0.5_dp, 1.0_dp, 1.5_dp, 2.0_dp]
      pot%rab = 0.5_dp
      pot%rho_nlcc(:) = [1.0_dp, 0.8_dp, 0.5_dp, 0.2_dp, 0.1_dp]
      pot%core_correction = .TRUE.
      pot%tau_nlcc(:) = [0.2_dp, 0.5_dp, 0.8_dp, 0.3_dp, 0.05_dp]
      pot%tau_nlcc_convention = 'HARTREE_PER_VOLUME'
      g(:, 1) = 0.0_dp
      g(:, 2) = [0.3_dp, -0.2_dp, 0.4_dp]
      g(:, 3) = -g(:, 2)
      position = [0.2_dp, -0.3_dp, 0.1_dp]
      delta_r = [0.02_dp, 0.01_dp, -0.03_dp]
      delta_g = 0.02_dp*g
      CALL upf_core_fourier(pot, g, position, 3.0_dp, rho, dr, dg)
      CALL check(ABS(3.0_dp*rho(1) - 4.0_dp*pi/6.0_dp* &
                     SUM([1.0_dp, 4.0_dp, 2.0_dp, 4.0_dp, 1.0_dp]*pot%r**2*pot%rho_nlcc)), "Original core zero mode", 1.e-12_dp)
      CALL check(ABS(rho(3) - CONJG(rho(2))), "Original core Fourier reality", 1.e-12_dp)
      response = SUM(dr*SPREAD(delta_r, 2, 3), DIM=1) + SUM(dg*delta_g, DIM=1) - 0.01_dp*rho
      DO side = 1, 2
         shift = REAL(2*side - 3, dp)*step
         CALL upf_core_fourier(pot, g + shift*delta_g, position + shift*delta_r, &
                               3.0_dp*(1.0_dp + 0.01_dp*shift), samples(:, side))
      END DO
      CALL check(MAXVAL(ABS((samples(:, 2) - samples(:, 1))/(2.0_dp*step) - response)), "Original core joint response", 2.e-9_dp)
      CALL upf_core_kinetic_fourier(pot, g, position, 3.0_dp, rho, dr, dg)
      CALL check(ABS(3.0_dp*rho(1) - 4.0_dp*pi/6.0_dp* &
                     SUM([1.0_dp, 4.0_dp, 2.0_dp, 4.0_dp, 1.0_dp]*pot%r**2*pot%tau_nlcc)), "Core tau zero mode", 1.e-12_dp)
      response = SUM(dr*SPREAD(delta_r, 2, 3), DIM=1) + SUM(dg*delta_g, DIM=1) - 0.01_dp*rho
      DO side = 1, 2
         shift = REAL(2*side - 3, dp)*step
         CALL upf_core_kinetic_fourier(pot, g + shift*delta_g, position + shift*delta_r, &
                                       3.0_dp*(1.0_dp + 0.01_dp*shift), samples(:, side))
      END DO
      CALL check(MAXVAL(ABS((samples(:, 2) - samples(:, 1))/(2.0_dp*step) - response)), "Core tau joint response", 2.e-9_dp)
      pot%core_correction = .FALSE.
      CALL upf_core_fourier(pot, g, position, 3.0_dp, rho, dr, dg)
      CALL check(MAXVAL(ABS(rho)) + MAXVAL(ABS(dr)) + MAXVAL(ABS(dg)), "Absent partial core", 0.0_dp)
   END SUBROUTINE check_core_fourier

! **************************************************************************************************
!> \brief Signed shell/PW coupling, reciprocal zero limit and joint geometry response.
! **************************************************************************************************
   SUBROUTINE check_fourier_coulomb()
      INTEGER, PARAMETER                                 :: ng = 4, nh = 3, nr = 2

      COMPLEX(dp)                                        :: derivative(ng, nr, nh), &
                                                            kernel(ng, nr, nh), other(ng, nr, nh), &
                                                            samples(ng, nr, nh, 2)
      INTEGER                                            :: l(nh), lm(nh), side
      REAL(dp)                                           :: dg(3, ng), dp_center(3), dr(nr), &
                                                            g(3, ng), position(3), radii(nr), &
                                                            shift, volume

      radii = [0.3_dp, 1.2_dp]
      l = [0, 1, 2]
      lm = [1, 3, 7]
      g(:, 1) = 0.0_dp
      g(:, 2) = [0.7_dp, 0.0_dp, 0.0_dp]
      g(:, 3) = [0.3_dp, -0.2_dp, 0.4_dp]
      g(:, 4) = -g(:, 3)
      position = [0.2_dp, -0.7_dp, 0.4_dp]
      volume = 7.0_dp
      dr = [0.07_dp, -0.03_dp]
      dg = 0.13_dp - 0.1_dp*g
      dp_center = [-0.04_dp, 0.08_dp, 0.02_dp]
      CALL upf_inverse_coulomb_fourier_kernel(radii, l, lm, g, position, volume, kernel, derivative, &
                                              dr, dg, dp_center, 0.3_dp)
      CALL check(MAXVAL(ABS(kernel(1, :, 1) - SQRT(4*pi)/volume)), 'shell monopole Fourier normalization')
      CALL check(MAXVAL(ABS(kernel(1, :, 2:))), 'zero higher shell multipoles at G=0', 0.0_dp)
      CALL check(MAXVAL(ABS(kernel(2, :, 1) - SQRT(4*pi)/volume* &
                            EXP(CMPLX(0.0_dp, -DOT_PRODUCT(g(:, 2), position), dp))* &
                            SIN(0.7_dp*radii)/(0.7_dp*radii))), 'analytic spherical shell form factor')
      CALL check(MAXVAL(ABS(kernel(3, :, :) - CONJG(kernel(4, :, :)))), 'real density reciprocal conjugacy')
      CALL check(MAXVAL(ABS(derivative(1, :, 2) - CMPLX(0.0_dp, &
                                               -SQRT(12*pi)/volume*radii*dg(3, 1)/3.0_dp, dp))), 'dipole Fourier derivative at G=0')
      DO side = 1, 2
         shift = step*(2*side - 3)
         CALL upf_inverse_coulomb_fourier_kernel(radii + shift*dr, l, lm, g + shift*dg, &
                                                 position + shift*dp_center, volume + shift*0.3_dp, samples(:, :, :, side))
      END DO
      CALL check(MAXVAL(ABS((samples(:, :, :, 2) - samples(:, :, :, 1))/(2*step) - derivative)), &
                 'shell Fourier radius center reciprocal volume response', 2.e-9_dp)
      CALL upf_inverse_coulomb_fourier_kernel(radii, l, lm, g, position, volume, other)
      CALL check(MAXVAL(ABS(kernel - other)), 'optional Fourier response preserves forward map', 0.0_dp)
   END SUBROUTINE check_fourier_coulomb

! **************************************************************************************************
!> \brief Signed Coulomb reciprocity, shell self term, moment compensation and energy adjoints.
! **************************************************************************************************
   SUBROUTINE check_radial_coulomb()
      INTEGER, PARAMETER                                 :: nh = 3, nr = 3

      INTEGER                                            :: channel, l(nh), side
      REAL(dp) :: adjoint(nr, nh), compensation(nr, nh), density(nr, nh), dr(nr), dw(nr), energy, &
         energy_other, moments(nh), one_density(1, 1), one_potential(1, 1), other(nr, nh), &
         potential(nr, nh), radii(nr), residual(nr, nh), response, samples(2), shapes(nr, nh), &
         shift, vp(nr, nh), weights(nr)

      ! An orthonormal monopole on one radial shell has E=2*pi*(w*rho)**2/r.
      one_density = 0.7_dp
      CALL upf_inverse_coulomb_radial([0.8_dp], [0.3_dp], [0], one_density, one_potential, energy)
      CALL check(ABS(energy - 2*pi*(0.3_dp*0.7_dp)**2/0.8_dp), 'Coulomb shell self energy')
      CALL upf_inverse_coulomb_radial([0.8_dp], [0.3_dp], [2], one_density, one_potential, energy)
      CALL check(ABS(energy - 2*pi*(0.3_dp*0.7_dp)**2/(5*0.8_dp)), 'Coulomb higher multipole shell')
      ! Unordered radii and signed, repeated-l channels are intentional.
      radii = [1.7_dp, 0.3_dp, 0.8_dp]
      weights = [0.4_dp, 0.2_dp, 0.7_dp]
      l = [0, 2, 2]
      density = RESHAPE([0.1_dp, -0.3_dp, 0.2_dp, -0.2_dp, 0.7_dp, -0.1_dp, &
                         0.3_dp, -0.2_dp, 0.4_dp], [nr, nh])
      other = 0.13_dp + 0.2_dp*density
      dr = [0.02_dp, -0.01_dp, 0.03_dp]
      dw = 0.03_dp*weights
      CALL upf_inverse_coulomb_radial(radii, weights, l, density, potential, energy)
      CALL upf_inverse_coulomb_radial(radii, weights, l, other, vp, energy_other)
      CALL check(ABS(SUM(SPREAD(weights, 2, nh)*(density*vp - other*potential))), 'signed Coulomb reciprocity')
      CALL check(MAX(0.0_dp, -energy), 'signed Coulomb nonnegative self energy', 0.0_dp)
      CALL upf_inverse_coulomb_radial_response(radii, weights, l, density, dr, dw, other, response)
      DO side = 1, 2
         shift = step*(2*side - 3)
         CALL upf_inverse_coulomb_radial(radii + shift*dr, weights + shift*dw, l, &
                                         density + shift*other, vp, samples(side))
      END DO
      CALL check(ABS((samples(2) - samples(1))/(2*step) - response), 'Coulomb density measure kernel response', 2.e-9_dp)
      DO channel = 1, nh
         shapes(:, channel) = radii**l(channel)*EXP(-0.6_dp*radii**2)
      END DO
      CALL upf_inverse_coulomb_compensate(radii, weights, l, density, shapes, compensation, residual, moments)
      CALL check(MAXVAL(ABS(density - compensation - residual)), 'Coulomb compensation density closure')
      DO channel = 1, nh
         CALL check(ABS(SUM(weights*radii**l(channel)*residual(:, channel))), 'zero compensated exterior moment')
         CALL check(ABS(moments(channel) - SUM(weights*radii**l(channel)*compensation(:, channel))), &
                    'original signed exterior moment retained')
      END DO
      CALL upf_inverse_coulomb_radial(radii, weights, l, residual, potential, energy)
      CALL upf_inverse_coulomb_compensation_adjoint(radii, weights, l, shapes, potential, adjoint)
      response = SUM(SPREAD(weights, 2, nh)*other*adjoint)
      DO side = 1, 2
         shift = step*(2*side - 3)
         CALL upf_inverse_coulomb_compensate(radii, weights, l, density + shift*other, shapes, &
                                             compensation, residual, moments)
         CALL upf_inverse_coulomb_radial(radii, weights, l, residual, vp, samples(side))
      END DO
      CALL check(ABS((samples(2) - samples(1))/(2*step) - response), 'occupation-dependent compensation adjoint', 2.e-9_dp)
      density = 0.0_dp
      CALL upf_inverse_coulomb_radial(radii, weights, l, density, potential, energy)
      CALL check(MAXVAL(ABS(potential)) + ABS(energy), 'exact zero NC Coulomb residual', 0.0_dp)
   END SUBROUTINE check_radial_coulomb

! **************************************************************************************************
!> \brief Gaussian/residual closure, adjoint and joint field/measure response on complex rows.
! **************************************************************************************************
   SUBROUTINE check_density_residual()
      INTEGER, PARAMETER                                 :: nb = 2, np = 3

      COMPLEX(dp) :: c(nb, nb), center_difference(np, nb), center_gradient(np, nb, 3), &
         dgr(np, nb, 3), dgs(np, nb, 3), dr(np, nb), ds(np, nb), gr(np, nb, 3), gs(np, nb, 3), &
         h(nb, nb), hr(nb, nb), hs(nb, nb), p(nb, nb), reference(np, nb), response(nb, nb), &
         soft(np, nb)
      INTEGER                                            :: axis, center, i, j, side
      REAL(dp) :: dw(np), gradient(np, 3), gradient_r(np, 3), gradient_s(np, 3), rho(np), &
         rho_r(np), rho_s(np), samples(2), shift, tau(np), tau_r(np), tau_s(np), value, vg(np, 3), &
         vn(np), vt(np), weights(np)

      c(:, 1) = [CMPLX(1.0_dp, 0.2_dp, dp), CMPLX(0.1_dp, -0.3_dp, dp)]
      c(:, 2) = [CMPLX(-0.4_dp, 0.3_dp, dp), CMPLX(0.2_dp, 0.5_dp, dp)]
      p = MATMUL(c, CONJG(TRANSPOSE(c)))
      DO j = 1, nb
         DO i = 1, np
            reference(i, j) = CMPLX(SIN(REAL(i + j, dp)), COS(REAL(i*j, dp)), dp)
            soft(i, j) = reference(i, j) + 0.1_dp*CONJG(reference(i, j))
            dr(i, j) = 0.07_dp*reference(i, j)
            ds(i, j) = 0.05_dp*soft(i, j)
            DO axis = 1, 3
               gr(i, j, axis) = CMPLX(COS(REAL(i + j + axis, dp)), SIN(REAL(i*j + axis, dp)), dp)
               gs(i, j, axis) = gr(i, j, axis) + 0.13_dp*CONJG(gr(i, j, axis))
               dgr(i, j, axis) = 0.03_dp*gr(i, j, axis)
               dgs(i, j, axis) = 0.09_dp*gs(i, j, axis)
            END DO
         END DO
      END DO
      weights = [0.3_dp, 0.5_dp, 0.7_dp]
      dw = 0.04_dp*weights
      vn = [0.2_dp, 0.7_dp, -0.3_dp]
      vt = [0.3_dp, 0.1_dp, 0.2_dp]
      vg = RESHAPE([0.1_dp, -0.2_dp, 0.3_dp, 0.4_dp, 0.5_dp, -0.6_dp, 0.7_dp, -0.8_dp, 0.9_dp], [np, 3])
      CALL upf_inverse_smooth_density(reference, gr, p, rho_r, gradient_r, tau_r)
      CALL upf_inverse_smooth_density(soft, gs, p, rho_s, gradient_s, tau_s)
      CALL upf_inverse_density_residual(reference, gr, soft, gs, p, rho, gradient, tau)
      CALL check(MAXVAL(ABS(rho_r + rho - rho_s)), 'Gaussian/residual density closure')
      CALL check(MAXVAL(ABS(gradient_r + gradient - gradient_s)), 'Gaussian/residual gradient closure')
      CALL check(MAXVAL(ABS(tau_r + tau - tau_s)), 'Gaussian/residual kinetic closure')
      CALL upf_inverse_semilocal_pair(reference, reference, gr, gr, weights, vn, vg, vt, hr)
      CALL upf_inverse_semilocal_pair(soft, soft, gs, gs, weights, vn, vg, vt, hs)
      CALL upf_inverse_semilocal_residual(reference, gr, soft, gs, weights, vn, vg, vt, h)
      CALL check(MAXVAL(ABS(hr + h - hs)), 'Gaussian/residual operator closure')
      value = SUM(weights*(vn*rho + SUM(vg*gradient, DIM=2) + vt*tau))
      CALL check(ABS(value - REAL(SUM(h*CONJG(p)), dp)), 'signed density operator adjoint')
      ! Overlapping center fields carry opposite off-diagonal terms that must not be discarded.
      rho_s = 0.0_dp
      gradient_s = 0.0_dp
      tau_s = 0.0_dp
      hs = CMPLX(0.0_dp, 0.0_dp, dp)
      DO center = 1, 2
         center_difference = (soft - reference)/2 + (2*center - 3)*0.2_dp*CONJG(reference)
         center_gradient = (gs - gr)/2 + (2*center - 3)*0.2_dp*CONJG(gr)
         CALL upf_inverse_cross_density(center_difference, soft + reference, center_gradient, gs + gr, &
                                        p, rho_r, gradient_r, tau_r)
         CALL upf_inverse_semilocal_pair(center_difference, soft + reference, center_gradient, gs + gr, &
                                         weights, vn, vg, vt, hr)
         hr = 0.5_dp*(hr + CONJG(TRANSPOSE(hr)))
         value = SUM(weights*(vn*rho_r + SUM(vg*gradient_r, DIM=2) + vt*tau_r))
         CALL check(ABS(value - REAL(SUM(hr*CONJG(p)), dp)), 'center residual operator adjoint')
         rho_s = rho_s + rho_r
         gradient_s = gradient_s + gradient_r
         tau_s = tau_s + tau_r
         hs = hs + hr
      END DO
      CALL check(MAXVAL(ABS(rho_s - rho)), 'overlapping center charge closure')
      CALL check(MAXVAL(ABS(gradient_s - gradient)), 'overlapping center gradient closure')
      CALL check(MAXVAL(ABS(tau_s - tau)), 'overlapping center kinetic closure')
      CALL check(MAXVAL(ABS(hs - h)), 'overlapping center operator closure')
      CALL upf_inverse_semilocal_response(soft - reference, soft + reference, gs - gr, gs + gr, &
                                          weights, vn, vg, vt, ds - dr, ds + dr, dgs - dgr, dgs + dgr, &
                                          dw, 0.05_dp*vn, 0.07_dp*vg, 0.09_dp*vt, response)
      DO side = 1, 2
         shift = step*(2*side - 3)
         CALL upf_inverse_density_residual(reference + shift*dr, gr + shift*dgr, soft + shift*ds, gs + shift*dgs, &
                                           p, rho, gradient, tau)
         samples(side) = SUM((weights + shift*dw)*((1.0_dp + 0.05_dp*shift)*vn*rho + &
                                                   (1.0_dp + 0.07_dp*shift)*SUM(vg*gradient, DIM=2) + &
                                                   (1.0_dp + 0.09_dp*shift)*vt*tau))
      END DO
      CALL check(ABS((samples(2) - samples(1))/(2*step) - REAL(SUM(response*CONJG(p)), dp)), &
                 'signed density geometry adjoint', 2.e-9_dp)
      CALL upf_inverse_density_residual(reference, gr, reference, gr, p, rho, gradient, tau)
      CALL check(MAXVAL(ABS(rho)) + MAXVAL(ABS(gradient)) + MAXVAL(ABS(tau)), 'exact zero NC residual', 0.0_dp)
      CALL upf_inverse_semilocal_residual(reference(:0, :), gr(:0, :, :), soft(:0, :), gs(:0, :, :), &
                                          weights(:0), vn(:0), vg(:0, :), vt(:0), h)
      CALL check(MAXVAL(ABS(h)), 'empty-row residual adjoint', 0.0_dp)
   END SUBROUTINE check_density_residual

! **************************************************************************************************
!> \brief Compare density congruence to independently projected occupied orbitals.
! **************************************************************************************************
   SUBROUTINE check_projector_density()
      COMPLEX(dp)                                        :: c(2, 2), density(2, 2), projected(3, 3), &
                                                            q(3, 3), reference(3, 3), w(3, 2), &
                                                            wc(3, 2)
      INTEGER                                            :: i, j

      DO j = 1, 2
         DO i = 1, 2
            c(i, j) = CMPLX(SIN(REAL(i + j, dp)), COS(REAL(i*j, dp)), dp)
         END DO
         DO i = 1, 3
            w(i, j) = 60.0_dp*CMPLX(COS(REAL(i + j, dp)), SIN(REAL(i*j, dp)), dp)
         END DO
      END DO
      density = MATMUL(c, CONJG(TRANSPOSE(c)))
      CALL upf_inverse_projector_density(w, density, projected)
      wc = MATMUL(w, c)
      reference = MATMUL(wc, CONJG(TRANSPOSE(wc)))
      CALL check(MAXVAL(ABS(projected - reference))/MAXVAL(ABS(reference)), 'projected occupied orbitals')
      CALL check(MAXVAL(ABS(projected - CONJG(TRANSPOSE(projected)))), 'Hermitian projector density', 0.0_dp)
      q = reference/MAXVAL(ABS(reference))
      CALL check(ABS(SUM(projected*CONJG(q)) - &
                     SUM(density*CONJG(MATMUL(CONJG(TRANSPOSE(w)), MATMUL(q, w)))))/MAXVAL(ABS(reference)), &
                 'projector density operator adjoint')
   END SUBROUTINE check_projector_density

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE check_inverse()
      INTEGER, PARAMETER                                 :: nb = 2, nf = 5, np = 3

      COMPLEX(KIND=dp) :: actual(nb, nb), b(np, nb), bd(np, np), chi(nf, nb), db(np, nb), &
         dchi(nf, nb), dd(np, np), ddd(np, np), ddelta(nf, np), ddg(np, nb), delta(nf, np), &
         dg(np, nb), dgd(nb, np), dgg(nb, nb), dk(np, np), donsite(np, np), dop(nf, nf), &
         dprojector(np, nf), ds(nf, nb), dw(np, nb), empty_b(0, nb), empty_delta(nf, 0), &
         empty_q(0, 0), empty_w(0, nb), gb(np, nb), gd(nb, np), gg(nb, nb), gk(np, np), &
         gw(np, nb), images(nf, nb, 2), onsite(np, np), op(nf, nf), projected(np, nb), &
         projector(np, nf), reference(nb, nb), response(nb, nb), sampled(nb, nb, 2), &
         shifted_chi(nf, nb)
      COMPLEX(KIND=dp) :: shifted_delta(nf, np), shifted_op(nf, nf), shifted_p(np, nf), &
         shifted_q(np, np), soft(nf, nb), w(np, nb)
      INTEGER                                            :: i, info, j, phase, side
      REAL(KIND=dp)                                      :: objective(2), pairing, shift
      TYPE(upf_inverse_type)                             :: inverse

      DO phase = 1, 3
         DO j = 1, nb
            DO i = 1, nf
               chi(i, j) = CMPLX(SIN(REAL(i + 2*j, dp)), (phase - 1)*COS(REAL(i*j, dp)), dp)
               dchi(i, j) = 0.07_dp*CMPLX(COS(REAL(i + j, dp)), (phase - 1)*SIN(REAL(2*i + j, dp)), dp)
            END DO
         END DO
         DO j = 1, np
            DO i = 1, nf
               delta(i, j) = 0.12_dp*CMPLX(COS(REAL(i*j, dp)), (phase - 1)*SIN(REAL(i + j, dp)), dp)
               ddelta(i, j) = 0.03_dp*CMPLX(SIN(REAL(i + 2*j, dp)), (phase - 1)*COS(REAL(i*j, dp)), dp)
               projector(j, i) = 0.15_dp*CMPLX(SIN(REAL(i + j, dp)), (phase - 1)*COS(REAL(2*i + j, dp)), dp)
               dprojector(j, i) = 0.02_dp*CMPLX(COS(REAL(i + j, dp)), (phase - 1)*SIN(REAL(i*j, dp)), dp)
            END DO
         END DO
         IF (phase == 3) THEN
            ! A mixed NC/USPP frame retains the NC bra, with no NC reconstruction.
            delta(:, np) = CMPLX(0.0_dp, 0.0_dp, dp)
            ddelta(:, np) = CMPLX(0.0_dp, 0.0_dp, dp)
         END IF
         ! A non-Hermitian operator exposes an accidental conjugation/symmetrization.
         DO j = 1, nf
            DO i = 1, nf
               op(i, j) = CMPLX(0.1_dp*i - 0.05_dp*j, (phase - 1)*0.03_dp*(i + j), dp)
               dop(i, j) = CMPLX(0.02_dp*(i + j), (phase - 1)*0.01_dp*(i - j), dp)
            END DO
         END DO
         onsite = op(1:np, 1:np)*0.3_dp
         donsite = dop(1:np, 1:np)*0.4_dp
         bd = MATMUL(projector, delta)
         b = MATMUL(projector, chi)
         CALL upf_inverse_init(bd, inverse, info)
         IF (info /= 0 .OR. .NOT. inverse%ready) ERROR STOP 'Inverse factorization failed'
         CALL upf_inverse_project(inverse, b, w)
         CALL upf_inverse_basis(chi, delta, w, soft)
         projected = MATMUL(projector, soft)
         CALL check(MAXVAL(ABS(projected - w)), 'inverse projector identity')
         CALL check(MAXVAL(ABS(soft + MATMUL(delta, projected) - chi)), 'full reconstruction round trip')
         db = MATMUL(dprojector, chi) + MATMUL(projector, dchi)
         dk = MATMUL(dprojector, delta) + MATMUL(projector, ddelta)
         CALL upf_inverse_project_response(inverse, w, db, dk, dw)
         gw = chi(1:np, :) + CMPLX(0.13_dp, -0.09_dp, dp)*b
         CALL upf_inverse_project_adjoint(inverse, w, gw, gb, gk)
         pairing = REAL(SUM(CONJG(gb)*db) + SUM(CONJG(gk)*dk), dp)
         CALL check(ABS(pairing - REAL(SUM(CONJG(gw)*dw), dp)), 'inverse solve forward/adjoint pairing')
         CALL check(MAXVAL(ABS(gb + MATMUL(CONJG(TRANSPOSE(bd)), gb) - gw)), 'inverse conjugate solve')
         CALL upf_inverse_basis_response(dchi, delta, ddelta, w, dw, ds)
         gg = MATMUL(CONJG(TRANSPOSE(chi)), MATMUL(op, chi))
         gd = MATMUL(CONJG(TRANSPOSE(chi)), MATMUL(op, delta))
         dg = MATMUL(CONJG(TRANSPOSE(delta)), MATMUL(op, chi))
         dd = MATMUL(CONJG(TRANSPOSE(delta)), MATMUL(op, delta))
         CALL upf_inverse_operator(gg, gd, dg, dd, onsite, w, actual)
         reference = MATMUL(CONJG(TRANSPOSE(soft)), MATMUL(op, soft)) + &
                     MATMUL(CONJG(TRANSPOSE(projected)), MATMUL(onsite, projected))
         CALL check(MAXVAL(ABS(actual - reference)), 'inverse operator versus direct fields')
         dgg = MATMUL(CONJG(TRANSPOSE(dchi)), MATMUL(op, chi)) + &
               MATMUL(CONJG(TRANSPOSE(chi)), MATMUL(dop, chi) + MATMUL(op, dchi))
         dgd = MATMUL(CONJG(TRANSPOSE(dchi)), MATMUL(op, delta)) + &
               MATMUL(CONJG(TRANSPOSE(chi)), MATMUL(dop, delta) + MATMUL(op, ddelta))
         ddg = MATMUL(CONJG(TRANSPOSE(ddelta)), MATMUL(op, chi)) + &
               MATMUL(CONJG(TRANSPOSE(delta)), MATMUL(dop, chi) + MATMUL(op, dchi))
         ddd = MATMUL(CONJG(TRANSPOSE(ddelta)), MATMUL(op, delta)) + &
               MATMUL(CONJG(TRANSPOSE(delta)), MATMUL(dop, delta) + MATMUL(op, ddelta))
         CALL upf_inverse_operator_response(dgg, gd, dg, dd, onsite, w, dgd, ddg, ddd, donsite, dw, response)
         DO side = 1, 2
            shift = (2*side - 3)*step
            shifted_chi = chi + shift*dchi
            shifted_delta = delta + shift*ddelta
            shifted_p = projector + shift*dprojector
            shifted_op = op + shift*dop
            shifted_q = onsite + shift*donsite
            CALL upf_inverse_init(MATMUL(shifted_p, shifted_delta), inverse, info)
            IF (info /= 0) ERROR STOP 'Displaced inverse factorization failed'
            CALL upf_inverse_project(inverse, MATMUL(shifted_p, shifted_chi), w)
            objective(side) = REAL(SUM(CONJG(gw)*w), dp)
            CALL upf_inverse_basis(shifted_chi, shifted_delta, w, images(:, :, side))
            projected = MATMUL(shifted_p, images(:, :, side))
            sampled(:, :, side) = MATMUL(CONJG(TRANSPOSE(images(:, :, side))), &
                                         MATMUL(shifted_op, images(:, :, side))) + &
                                  MATMUL(CONJG(TRANSPOSE(projected)), MATMUL(shifted_q, projected))
         END DO
         CALL check(MAXVAL(ABS(ds - (images(:, :, 2) - images(:, :, 1))/(2*step))), 'inverse basis response', 1.e-9_dp)
         CALL check(MAXVAL(ABS(response - (sampled(:, :, 2) - sampled(:, :, 1))/(2*step))), &
                    'inverse operator response', 1.e-9_dp)
         CALL check(ABS(pairing - (objective(2) - objective(1))/(2*step)), 'inverse solve energy adjoint', 1.e-9_dp)
      END DO

      CALL upf_inverse_init(empty_q, inverse, info)
      IF (info /= 0 .OR. .NOT. inverse%ready) ERROR STOP 'Empty inverse failed'
      CALL upf_inverse_project(inverse, empty_b, empty_w)
      CALL upf_inverse_project_adjoint(inverse, empty_b, empty_b, empty_w, empty_q)
      CALL upf_inverse_basis(chi, empty_delta, empty_w, soft)
      CALL check(MAXVAL(ABS(soft - chi)), 'identity without reconstruction projectors')
      bd = CMPLX(0.0_dp, 0.0_dp, dp)
      bd(2, 2) = -1.0_dp
      CALL upf_inverse_init(bd, inverse, info)
      IF (info /= 2 .OR. inverse%ready) ERROR STOP 'Singular inverse accepted'
      bd(2, 2) = -1.0_dp + 1.e-12_dp
      CALL upf_inverse_init(bd, inverse, info, minimum_rcond=1.e-10_dp)
      IF (info /= np + 1 .OR. inverse%ready) ERROR STOP 'Ill-conditioned inverse accepted'
      bd(1, 1) = CMPLX(ieee_value(0.0_dp, ieee_quiet_nan), 0.0_dp, dp)
      CALL upf_inverse_init(bd, inverse, info)
      IF (info /= np + 2 .OR. inverse%ready) ERROR STOP 'Nonfinite inverse accepted'
   END SUBROUTINE check_inverse

! **************************************************************************************************
!> \brief NC original projectors need neither reference waves nor augmentation data.
! **************************************************************************************************
   SUBROUTINE check_nc_fields()
      REAL(dp)                                           :: displacement(3), gradient(3, 1), value(1)
      TYPE(atom_upfpot_type)                             :: pot
      TYPE(upf_radial_fields_type)                       :: beta, delta

      pot%pseudo_type = "NC"
      ALLOCATE (pot%r(7), pot%beta(7, 1), pot%lbeta(1))
      pot%r(:) = [0.0_dp, 0.2_dp, 0.5_dp, 1.0_dp, 1.5_dp, 2.0_dp, 3.0_dp]
      pot%lbeta = 0
      pot%beta(:, 1) = EXP(-pot%r**2)
      pot%beta(6:, 1) = 0.0_dp
      CALL upf_inverse_fields_init(pot, beta, delta)
      displacement = [0.5_dp, 0.0_dp, 0.0_dp]
      CALL upf_radial_fields_evaluate(beta, displacement, value)
      CALL check(ABS(value(1) - EXP(-0.25_dp)/SQRT(4.0_dp*pi)), 'original NC projector')
      CALL upf_radial_fields_evaluate(delta, displacement, value, gradient)
      CALL check(MAXVAL(ABS(value)) + MAXVAL(ABS(gradient)), 'zero NC reconstruction')
   END SUBROUTINE check_nc_fields

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE check_gaussian_fields()
      INTEGER, PARAMETER                                 :: nb = 11, nf = 5

      COMPLEX(dp) :: box_cell_gradient(nf, nb, 3, 3), dg(nf, nb, 3), dv(nf, nb), g(nf, nb, 3), &
         gss(nf, nb, 3, 2), other(nf, nb), phase, reference_cell_gradient(nf, nb, 3, 3), &
         small(2, nb), v(nf, nb), vs(nf, nb, 2)
      INTEGER :: axis, box_images(3, 27), column, ico, image, images(3, 3), ix, iy, iz, kernel, &
         point, row, side, translated_images(3, 3), unit
      REAL(dp) :: batch_gradient(3, 8, nf), batch_values(nf, 8), bloch(3), box_cell(3, 3), &
         cell(3, 3), center(3), cg(3, 8), cm(8), cp(8), cv(8), cx(3), dc(3, 3), dcenter(3), dk(3), &
         dpnt(3, nf), expected, gradient(nf, nb, 3), gs(nf, nb, 3, 2), hh(nf, nb, 3, 3), &
         lattice_gradient(3, 3, 8), lg(3, 8), lv(8), omega, radius, samples(nf, nb, 2), sg(3, 8), &
         shift, shifted_cell(3, 3), sv(8), value(nf, nb), width(2), x(3, nf), xp(3, nf)
      TYPE(gaussian_basis_workspace_type)                :: workspace
      TYPE(gaussian_coulomb_fields_type)                 :: coulomb, identity, long_range, &
                                                            periodic_fields, short_range, truncated
      TYPE(gto_basis_set_type), POINTER                  :: basis
      TYPE(mp_comm_type)                                 :: comm

      NULLIFY (basis)
      CALL init_orbital_pointers(3)
      CALL allocate_gto_basis_set(basis)
      basis%nset = 2
      basis%nsgf = nb
      ALLOCATE (basis%lmax(2), basis%npgf(2), basis%nsgf_set(2), basis%nshell(2), basis%l(2, 2), &
                basis%first_sgf(2, 2), basis%last_sgf(2, 2), basis%zet(2, 2), basis%sphi(40, nb))
      basis%lmax = [3, 1]
      basis%npgf = [2, 1]
      basis%nsgf_set = [8, 3]
      basis%nshell = [2, 1]
      basis%l = 0
      basis%l(:, 1) = [0, 3]
      basis%l(1, 2) = 1
      basis%first_sgf = 0
      basis%last_sgf = 0
      basis%first_sgf(:, 1) = [1, 2]
      basis%last_sgf(:, 1) = [1, 8]
      basis%first_sgf(1, 2) = 9
      basis%last_sgf(1, 2) = 11
      basis%zet = 0.0_dp
      basis%zet(:, 1) = [.4_dp, .9_dp]
      basis%zet(1, 2) = .7_dp
      basis%sphi = 0.0_dp
      basis%sphi(1, 1) = .6_dp
      basis%sphi(21, 1) = -.2_dp
      DO row = 2, 8
         DO ico = ncoset(2) + 1, ncoset(3)
            basis%sphi(ico, row) = .03_dp*SIN(REAL(ico*row, dp))
            basis%sphi(ico + 20, row) = -.02_dp*COS(REAL(ico + row, dp))
         END DO
      END DO
      DO row = 9, nb
         DO ico = 2, 4
            basis%sphi(ico, row) = .1_dp*SIN(REAL(ico*row, dp))
         END DO
      END DO
      x(:, 1) = 0.0_dp
      x(:, 2) = [.3_dp, 0.0_dp, 0.0_dp]
      x(:, 3) = [0.0_dp, -.7_dp, 0.0_dp]
      x(:, 4) = [.23_dp, .51_dp, -.4_dp]
      x(:, 5) = [-.6_dp, .2_dp, .8_dp]
      CALL gaussian_basis_fields_evaluate(basis, x, value, gradient, hh)
      CALL check(MAXVAL(ABS(value(:, 1) - (.6_dp*EXP(-.4_dp*SUM(x*x, DIM=1)) &
                                           - .2_dp*EXP(-.9_dp*SUM(x*x, DIM=1))))), 'native signed s contraction')
      DO axis = 1, 3
         DO side = 1, 2
            xp = x
            xp(axis, :) = xp(axis, :) + (2*side - 3)*step
            CALL gaussian_basis_fields_evaluate(basis, xp, samples(:, :, side), gs(:, :, :, side))
         END DO
         CALL check(MAXVAL(ABS(gradient(:, :, axis) - (samples(:, :, 2) - samples(:, :, 1))/(2*step))), &
                    'contracted Gaussian point gradient', 2.e-8_dp)
         CALL check(MAXVAL(ABS(hh(:, :, :, axis) - (gs(:, :, :, 2) - gs(:, :, :, 1))/(2*step))), &
                    'contracted Gaussian origin/axis Hessian', 2.e-7_dp)
      END DO
      cell = RESHAPE([2.3_dp, .1_dp, .2_dp, .2_dp, 2.7_dp, -.1_dp, .1_dp, .3_dp, 2.1_dp], [3, 3])
      center = [.1_dp, -.2_dp, .3_dp]
      images = RESHAPE([0, 0, 0, 1, -1, 0, -1, 0, 1], [3, 3])
      bloch = [.23_dp, -.17_dp, .11_dp]
      dc = .03_dp*cell
      dcenter = [.07_dp, -.02_dp, .01_dp]
      dpnt = .02_dp*x
      dk = [.03_dp, .02_dp, -.01_dp]
      CALL gaussian_basis_fields_periodic(basis, x, center, cell, images, bloch, v, g, dv, dg, dpnt, dcenter, dc, dk)
      CALL gaussian_basis_fields_periodic(basis, x, center, cell, images, bloch, vs(:, :, 1), gss(:, :, :, 1), &
                                          vs(:, :, 2), gss(:, :, :, 2), dpnt, dcenter, dc, dk, workspace=workspace)
      CALL check(MAX(MAXVAL(ABS(vs(:, :, 1) - v)), MAXVAL(ABS(vs(:, :, 2) - dv))), &
                 'reused Gaussian value and joint response')
      CALL check(MAX(MAXVAL(ABS(gss(:, :, :, 1) - g)), MAXVAL(ABS(gss(:, :, :, 2) - dg))), &
                 'reused Gaussian gradient and Hessian response')
      CALL gaussian_basis_fields_periodic(basis, x(:, :2), center, cell, images, bloch, small, workspace=workspace)
      CALL check(MAXVAL(ABS(small - v(:2, :))), 'resized Gaussian value-only workspace')
      DO side = 1, 2
         shift = (2*side - 3)*step
         CALL gaussian_basis_fields_periodic(basis, x + shift*dpnt, center + shift*dcenter, cell + shift*dc, &
                                             images, bloch + shift*dk, vs(:, :, side), gss(:, :, :, side))
      END DO
      CALL check(MAXVAL(ABS(dv - (vs(:, :, 2) - vs(:, :, 1))/(2*step))), 'native Gaussian Bloch joint response', 2.e-8_dp)
     CALL check(MAXVAL(ABS(dg - (gss(:, :, :, 2) - gss(:, :, :, 1))/(2*step))), 'native Gaussian Bloch gradient response', 2.e-7_dp)
      CALL gaussian_basis_fields_periodic(basis, x, center, cell, images, -bloch, other, workspace=workspace)
      CALL check(MAXVAL(ABS(other - CONJG(v))), 'native Gaussian Bloch time reversal')
      DO image = 1, 3
         translated_images(:, image) = images(:, image) - [1, 0, 0]
      END DO
      CALL gaussian_basis_fields_periodic(basis, x, center + cell(:, 1), cell, translated_images, bloch, other)
      CALL check(MAXVAL(ABS(other - EXP(CMPLX(0.0_dp, -2*pi*bloch(1), dp))*v)), 'native Gaussian lattice representative')
      CALL gaussian_coulomb_fields_init(3)
      CALL gaussian_coulomb_fields_prepare(basis, 1, coulomb)
      CALL gaussian_coulomb_fields_prepare(basis, 1, identity, identity=.TRUE.)
      expected = 0.0_dp
      DO point = 1, nf
         CALL gaussian_coulomb_fields_evaluate(identity, x(:, point), cv, cg)
         expected = MAX(expected, MAXVAL(ABS(cv - value(point, :8))), &
                        MAXVAL(ABS(cg - TRANSPOSE(gradient(point, :8, :)))))
      END DO
      CALL check(expected, 'overlap RI field and independent Cartesian Gaussian gradient')
      CALL gaussian_coulomb_fields_evaluate_batch(identity, x, batch_values, batch_gradient)
      CALL check(MAXVAL(ABS(batch_values - value(:, :8))), 'batched overlap RI fields')
      DO point = 1, nf
         CALL check(MAXVAL(ABS(batch_gradient(:, :, point) - TRANSPOSE(gradient(point, :8, :)))), &
                    'batched overlap RI field gradients')
      END DO
      CALL gaussian_coulomb_fields_evaluate_batch(identity, x, batch_values)
      CALL check(MAXVAL(ABS(batch_values - value(:, :8))), 'batched value-only overlap RI fields')
      box_cell = 0.0_dp
      DO axis = 1, 3
         box_cell(axis, axis) = 2.0_dp + 0.3_dp*REAL(axis, dp)
      END DO
      image = 0
      DO iz = -1, 1
      DO iy = -1, 1
      DO ix = -1, 1
         image = image + 1
         box_images(:, image) = [ix, iy, iz]
      END DO
      END DO
      END DO
      CALL gaussian_basis_fields_periodic(basis, x, center, box_cell, box_images, bloch, v, g, &
                                          cell_gradient=box_cell_gradient)
      other = CMPLX(0.0_dp, 0.0_dp, dp)
      dg = CMPLX(0.0_dp, 0.0_dp, dp)
      reference_cell_gradient = CMPLX(0.0_dp, 0.0_dp, dp)
      dv = CMPLX(0.0_dp, 0.0_dp, dp)
      DO image = 1, SIZE(box_images, 2)
         xp = x - SPREAD(center + MATMUL(box_cell, REAL(box_images(:, image), dp)), 2, nf)
         CALL gaussian_basis_fields_evaluate(basis, xp, value, gradient)
         phase = EXP(CMPLX(0.0_dp, 2*pi*DOT_PRODUCT(bloch, REAL(box_images(:, image), dp)), dp))
         other = other + phase*value
         dg = dg + phase*gradient
         dv = dv + CMPLX(0.0_dp, 2*pi*DOT_PRODUCT(dk, REAL(box_images(:, image), dp)), dp)*phase*value
         DO axis = 1, 3
            dv = dv + phase*gradient(:, :, axis)*SPREAD(dpnt(axis, :) - dcenter(axis) - &
                                                        DOT_PRODUCT(dc(axis, :), REAL(box_images(:, image), dp)), 2, nb)
         END DO
         DO column = 1, 3
            reference_cell_gradient(:, :, :, column) = reference_cell_gradient(:, :, :, column) - &
                                                       phase*gradient*REAL(box_images(column, image), dp)
         END DO
      END DO
      CALL check(MAXVAL(ABS(v - other)) + MAXVAL(ABS(g - dg)), 'factored Gaussian value and spatial gradient')
      CALL check(MAXVAL(ABS(box_cell_gradient - reference_cell_gradient)), 'factored Gaussian full cell gradient')
      CALL gaussian_basis_fields_periodic(basis, x, center, box_cell, box_images, bloch, v, response=other, &
                                          point_response=dpnt, center_response=dcenter, cell_response=dc, reduced_k_response=dk)
      CALL check(MAXVAL(ABS(other - dv)), 'factored Gaussian joint point center cell and Bloch response')
      dv = CMPLX(0.0_dp, 0.0_dp, dp)
      DO axis = 1, 3
         dv = dv - dcenter(axis)*dg(:, :, axis)
      END DO
      CALL gaussian_basis_fields_periodic(basis, x, center, box_cell, box_images, bloch, v, &
                                          response=other, center_response=dcenter)
      CALL check(MAXVAL(ABS(other - dv)), 'factored Gaussian center-only response')
      DO column = 1, 3
         DO axis = 1, 3
            DO side = 1, 2
               shifted_cell = box_cell
               shifted_cell(axis, column) = shifted_cell(axis, column) + REAL(2*side - 3, dp)*step
               CALL gaussian_basis_fields_periodic(basis, x, center, shifted_cell, box_images, bloch, vs(:, :, side))
            END DO
            CALL check(MAXVAL(ABS((vs(:, :, 2) - vs(:, :, 1))/(2*step) - box_cell_gradient(:, :, axis, column))), &
                       'factored Gaussian diagonal and shear cell differences', 2.e-8_dp)
         END DO
      END DO
      CALL deallocate_gto_basis_set(basis)
      omega = .31_dp
      width = SQRT([.4_dp, .9_dp]*omega**2/([.4_dp, .9_dp] + omega**2))
      long_range = coulomb
      short_range = coulomb
      CALL gaussian_coulomb_fields_screen(long_range, omega, .FALSE.)
      CALL gaussian_coulomb_fields_screen(short_range, omega, .TRUE.)
      DO point = 1, 3
         cx = [SQRT(12.0_dp/.4_dp)*REAL(point - 1, dp), 0.0_dp, 0.0_dp]
         CALL gaussian_coulomb_fields_evaluate(coulomb, cx, cv, cg)
         CALL gaussian_coulomb_fields_evaluate(long_range, cx, lv, lg)
         CALL gaussian_coulomb_fields_evaluate(short_range, cx, sv, sg)
         CALL check(MAXVAL(ABS(lv + sv - cv)), 'Gaussian screened potential partition')
         CALL check(MAXVAL(ABS(lg + sg - cg)), 'Gaussian screened gradient partition')
         radius = SQRT(SUM(cx**2))
         IF (radius == 0.0_dp) THEN
            CALL check(ABS(cv(1) - 2*pi*(.6_dp/.4_dp - .2_dp/.9_dp)), 'Gaussian Coulomb origin')
            expected = 2*pi*(.6_dp*width(1)/.4_dp**1.5_dp - .2_dp*width(2)/.9_dp**1.5_dp)
         ELSE
            CALL check(ABS(cv(1) - pi**1.5_dp*(.6_dp*ERF(SQRT(.4_dp)*radius)/.4_dp**1.5_dp &
                                               - .2_dp*ERF(SQRT(.9_dp)*radius)/.9_dp**1.5_dp)/radius), &
                       'Signed Gaussian Coulomb tail')
            expected = pi**1.5_dp*(.6_dp*ERF(width(1)*radius)/.4_dp**1.5_dp &
                                   - .2_dp*ERF(width(2)*radius)/.9_dp**1.5_dp)/radius
         END IF
         CALL check(ABS(lv(1) - expected), 'Signed Gaussian erf potential')
         DO axis = 1, 3
            cx(axis) = cx(axis) + step
            CALL gaussian_coulomb_fields_evaluate(coulomb, cx, cp)
            cx(axis) = cx(axis) - 2*step
            CALL gaussian_coulomb_fields_evaluate(coulomb, cx, cm)
            cx(axis) = cx(axis) + step
            CALL check(MAXVAL(ABS((cp - cm)/(2*step) - cg(axis, :))), 'Gaussian Coulomb gradient', 2.e-8_dp)
            cx(axis) = cx(axis) + step
            CALL gaussian_coulomb_fields_evaluate(long_range, cx, cp)
            cx(axis) = cx(axis) - 2*step
            CALL gaussian_coulomb_fields_evaluate(long_range, cx, cm)
            cx(axis) = cx(axis) + step
            CALL check(MAXVAL(ABS((cp - cm)/(2*step) - lg(axis, :))), 'Gaussian erf gradient', 2.e-8_dp)
            cx(axis) = cx(axis) + step
            CALL gaussian_coulomb_fields_evaluate(short_range, cx, cp)
            cx(axis) = cx(axis) - 2*step
            CALL gaussian_coulomb_fields_evaluate(short_range, cx, cm)
            cx(axis) = cx(axis) + step
            CALL check(MAXVAL(ABS((cp - cm)/(2*step) - sg(axis, :))), 'Gaussian erfc gradient', 2.e-8_dp)
         END DO
      END DO
      CALL mp_world_init(comm)
      CALL check_atomic_projectors(comm)
      unit = -1
      IF (comm%mepos == 0) CALL open_file(file_name='t_c_g.dat', unit_number=unit)
      CALL init_truncated(3, unit, comm%mepos, comm)
      IF (comm%mepos == 0) CALL close_file(unit)
      truncated = coulomb
      CALL gaussian_coulomb_fields_truncate(truncated, 2.3_dp)
      DO kernel = 1, 3
         periodic_fields = coulomb
         IF (kernel == 2) periodic_fields = short_range
         IF (kernel == 3) periodic_fields = truncated
         CALL gaussian_coulomb_fields_evaluate_batch(periodic_fields, x, batch_values, batch_gradient)
         DO point = 1, nf
            CALL gaussian_coulomb_fields_evaluate(periodic_fields, x(:, point), cv, cg)
            CALL check(MAXVAL(ABS(batch_values(point, :) - cv)) + MAXVAL(ABS(batch_gradient(:, :, point) - cg)), &
                       'batched interaction fields and gradients')
         END DO
      END DO
      DO point = 1, 4
         cx = REAL(point - 1, dp)*[1.15_dp, .13_dp, .07_dp]
         radius = NORM2(cx)
         CALL gaussian_coulomb_fields_evaluate(truncated, cx, sv, sg)
         IF (radius == 0.0_dp) THEN
            expected = 2*pi*(.6_dp*(1 - EXP(-.4_dp*2.3_dp**2))/.4_dp - &
                             .2_dp*(1 - EXP(-.9_dp*2.3_dp**2))/.9_dp)
         ELSE
            expected = pi**1.5_dp/(2*radius)* &
                       (.6_dp*(2*ERF(SQRT(.4_dp)*radius) + ERF(SQRT(.4_dp)*(2.3_dp - radius)) - &
                               ERF(SQRT(.4_dp)*(2.3_dp + radius)))/.4_dp**1.5_dp - &
                        .2_dp*(2*ERF(SQRT(.9_dp)*radius) + ERF(SQRT(.9_dp)*(2.3_dp - radius)) - &
                               ERF(SQRT(.9_dp)*(2.3_dp + radius)))/.9_dp**1.5_dp)
         END IF
         CALL check(ABS(sv(1) - expected), 'Signed truncated Gaussian potential', 3.e-8_dp)
         DO axis = 1, 3
            cx(axis) = cx(axis) + step
            CALL gaussian_coulomb_fields_evaluate(truncated, cx, cp)
            cx(axis) = cx(axis) - 2*step
            CALL gaussian_coulomb_fields_evaluate(truncated, cx, cm)
            cx(axis) = cx(axis) + step
            CALL check(MAXVAL(ABS((cp - cm)/(2*step) - sg(axis, :))), 'Truncated Gaussian gradient', 2.e-8_dp)
         END DO
      END DO
      long_range = coulomb
      CALL gaussian_coulomb_fields_truncate(long_range, 100.0_dp)
      CALL gaussian_coulomb_fields_evaluate(long_range, cx, lv, lg)
      CALL gaussian_coulomb_fields_evaluate(coulomb, cx, cv, cg)
      CALL check(MAXVAL(ABS(lv - cv)) + MAXVAL(ABS(lg - cg)), 'Truncated Gaussian full-Coulomb limit', 3.e-12_dp)
      DO kernel = 1, 2
         periodic_fields = short_range
         IF (kernel == 2) periodic_fields = truncated
         cx = [.23_dp, -.41_dp, .37_dp]
         CALL gaussian_coulomb_fields_periodic(periodic_fields, cx, cell, images, sv, sg, lattice_gradient)
         cv = 0.0_dp
         DO image = 1, SIZE(images, 2)
            CALL gaussian_coulomb_fields_evaluate(periodic_fields, cx - MATMUL(cell, REAL(images(:, image), dp)), cp)
            cv = cv + cp
         END DO
         CALL check(MAXVAL(ABS(sv - cv)), 'Gaussian auxiliary image ordering')
         DO column = 1, 3
            DO axis = 1, 3
               shifted_cell = cell
               shifted_cell(axis, column) = shifted_cell(axis, column) + step
               CALL gaussian_coulomb_fields_periodic(periodic_fields, cx, shifted_cell, images, cp)
               shifted_cell(axis, column) = shifted_cell(axis, column) - 2*step
               CALL gaussian_coulomb_fields_periodic(periodic_fields, cx, shifted_cell, images, cm)
               CALL check(MAXVAL(ABS((cp - cm)/(2*step) - lattice_gradient(axis, column, :))), &
                          'Gaussian auxiliary cell gradient', 2.e-8_dp)
            END DO
         END DO
         CALL gaussian_coulomb_fields_periodic(periodic_fields, cx + cell(:, 1), cell, images + &
                                               SPREAD([1, 0, 0], 2, SIZE(images, 2)), cp)
         CALL check(MAXVAL(ABS(cp - sv)), 'Gaussian auxiliary lattice representative')
         CALL gaussian_coulomb_fields_periodic(periodic_fields, cx, cell, images(:, :0), cp, lg, lattice_gradient)
         CALL check(MAXVAL(ABS(cp)) + MAXVAL(ABS(lg)) + MAXVAL(ABS(lattice_gradient)), 'Empty auxiliary image sum')
      END DO
      long_range = coulomb
      short_range = coulomb
      CALL gaussian_coulomb_fields_screen(long_range, 0.0_dp, .FALSE.)
      CALL gaussian_coulomb_fields_screen(short_range, 0.0_dp, .TRUE.)
      CALL gaussian_coulomb_fields_evaluate(long_range, cx, lv, lg)
      CALL gaussian_coulomb_fields_evaluate(short_range, cx, sv, sg)
      CALL gaussian_coulomb_fields_evaluate(coulomb, cx, cv, cg)
      CALL check(MAXVAL(ABS(lv)) + MAXVAL(ABS(lg)), 'Zero-screening erf field')
      CALL check(MAXVAL(ABS(sv - cv)) + MAXVAL(ABS(sg - cg)), 'Zero-screening erfc field')
      CALL free_C0()
      CALL mp_world_finalize()
   END SUBROUTINE check_gaussian_fields

! **************************************************************************************************
!> \brief ...
!> \param value ...
!> \param label ...
!> \param tolerance ...
! **************************************************************************************************
   SUBROUTINE check(value, label, tolerance)
      REAL(KIND=dp), INTENT(IN)                          :: value
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: tolerance

      REAL(KIND=dp)                                      :: tol

      tol = 2.e-12_dp
      IF (PRESENT(tolerance)) tol = tolerance
      max_error = MAX(max_error, value)
      IF (value > tol) THEN
         WRITE (*, *) TRIM(label), value, tol
         ERROR STOP 'Inverse primitive mismatch'
      END IF
   END SUBROUTINE check
END PROGRAM upf_inverse_unittest
