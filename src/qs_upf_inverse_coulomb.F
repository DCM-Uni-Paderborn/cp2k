!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
!> \brief Native periodic GAPW Coulomb energy, AO adjoint and electronic charge metric.
MODULE qs_upf_inverse_coulomb
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type,&
                                              dbcsr_type
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE mathconstants,                   ONLY: fourpi,&
                                              pi
   USE mathlib,                         ONLY: invmat
   USE message_passing,                 ONLY: mp_para_env_type
   USE mt_util,                         ONLY: MT0D,&
                                              MT2D
   USE orbital_pointers,                ONLY: indso
   USE particle_types,                  ONLY: particle_type
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_grid_types,                   ONLY: FULLSPACE,&
                                              HALFSPACE,&
                                              PW_MODE_LOCAL,&
                                              pw_grid_type
   USE pw_methods,                      ONLY: pw_integral_ab,&
                                              pw_integrate_function,&
                                              pw_real_fourier_aliases,&
                                              pw_scale,&
                                              pw_transfer
   USE pw_poisson_methods,              ONLY: pw_green_strain_response,&
                                              pw_poisson_rebuild
   USE pw_poisson_types,                ONLY: ANALYTIC0D,&
                                              ANALYTIC1D,&
                                              ANALYTIC2D,&
                                              PERIODIC3D,&
                                              pw_poisson_type
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE qs_collocate_density,            ONLY: calculate_rho_core
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_grid_atom,                    ONLY: grid_atom_type
   USE qs_harmonics_atom,               ONLY: harmonics_atom_type
   USE qs_integrate_potential_single,   ONLY: integrate_v_core_rspace
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_oce_methods,                  ONLY: build_oce_matrices
   USE qs_oce_types,                    ONLY: allocate_oce_set,&
                                              create_oce_set,&
                                              deallocate_oce_set,&
                                              oce_matrix_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_upf_inverse,                  ONLY: qs_upf_inverse_response
   USE qs_upf_inverse_background,       ONLY: qs_upf_inverse_background_collocate,&
                                              qs_upf_inverse_background_geometry,&
                                              qs_upf_inverse_background_integrate,&
                                              qs_upf_inverse_background_prepare,&
                                              qs_upf_inverse_background_release,&
                                              qs_upf_inverse_background_scratch
   USE qs_upf_inverse_density,          ONLY: qs_upf_inverse_density_cache_prepare,&
                                              qs_upf_inverse_density_gapw,&
                                              qs_upf_inverse_density_gapw_integrate
   USE qs_upf_inverse_gaussian_pairs,   ONLY: qs_upf_inverse_gaussian_pairs_has_derivatives,&
                                              qs_upf_inverse_gaussian_pairs_prepare
   USE qs_upf_inverse_grid,             ONLY: qs_upf_inverse_grid_charge_type,&
                                              qs_upf_inverse_grid_integrate_dense,&
                                              qs_upf_inverse_grid_response,&
                                              qs_upf_inverse_matrix_gather
   USE qs_upf_inverse_kpoints,          ONLY: qs_upf_inverse_kpoints_density
   USE qs_upf_inverse_paw,              ONLY: qs_upf_inverse_paw_density_adjoint,&
                                              qs_upf_inverse_paw_density_fourier
   USE qs_upf_inverse_scf,              ONLY: qs_upf_inverse_scf_install_metric,&
                                              qs_upf_inverse_scf_integrate_background
   USE qs_upf_inverse_types,            ONLY: qs_upf_inverse_augmentation_type,&
                                              qs_upf_inverse_background_type,&
                                              qs_upf_inverse_coulomb_release,&
                                              qs_upf_inverse_coulomb_type,&
                                              qs_upf_inverse_scf_type,&
                                              qs_upf_inverse_type
   USE upf_inverse_coulomb,             ONLY: upf_inverse_coulomb_center_type,&
                                              upf_inverse_coulomb_density,&
                                              upf_inverse_coulomb_density_adjoint,&
                                              upf_inverse_coulomb_pair_geometry,&
                                              upf_inverse_coulomb_periodic,&
                                              upf_inverse_coulomb_potential_type,&
                                              upf_inverse_coulomb_prepare
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: qs_upf_inverse_coulomb_type, qs_upf_inverse_coulomb_prepare
   PUBLIC :: qs_upf_inverse_coulomb_release, qs_upf_inverse_coulomb_evaluate, qs_upf_inverse_coulomb_metric
   PUBLIC :: qs_upf_inverse_coulomb_install_metric
   PUBLIC :: qs_upf_inverse_coulomb_attach, qs_upf_inverse_coulomb_update
   PUBLIC :: qs_upf_inverse_coulomb_xc_integrate
   PUBLIC :: qs_upf_inverse_coulomb_ionic_geometry
   PUBLIC :: qs_upf_inverse_coulomb_pair_geometry
   PUBLIC :: qs_upf_inverse_coulomb_response, qs_upf_inverse_coulomb_metric_response
   PUBLIC :: qs_upf_inverse_coulomb_metric_geometry, qs_upf_inverse_coulomb_center_metric_response
   PUBLIC :: qs_upf_inverse_coulomb_operator_response
   PUBLIC :: qs_upf_inverse_coulomb_density
   PUBLIC :: qs_upf_inverse_coulomb_charges
   PUBLIC :: qs_upf_inverse_paw_density, qs_upf_inverse_electronic_density
   PUBLIC :: qs_upf_inverse_local_potential, qs_upf_inverse_local_field

CONTAINS

! **************************************************************************************************
!> \brief Couple a scalar grid potential to density corrections missing from the native grid source.
!>        Full GAPW adds every signed center. PAW additionally adds AE-minus-PS-minus-Q and the
!>        frozen AE core once. The ordinary grid contribution and ionic counterterms remain with
!>        the caller. The same Fourier map supplies energy, its AO transpose and joint geometry
!>        response. All ranks participate; Bloch weights enter occupations once, never operators.
!> \param qs_env Current native density, inverse frames and original atomic data.
!> \param potential Unweighted scalar field on the native auxiliary grid.
!> \param energy Correction energy, replicated over the grid communicator.
!> \param operators Optional spin-independent AO operators, indexed by global k without weights.
!> \param response Optional joint atom/cell energy derivative at fixed AO occupations.
!> \param atom_response Physical Cartesian atomic motion, required with response.
!> \param cell_response Lattice derivative with vectors in columns, required with response.
!> \param coefficients_response Defining inverse coefficient responses, required with response.
!> \param potential_response Optional derivative of the sampled potential on the moving grid.
!> \param core_energy Optional full frozen-core energy, supplied with operators and counted once.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_local_potential(qs_env, potential, energy, operators, response, atom_response, &
                                             cell_response, coefficients_response, potential_response, core_energy)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: potential
      REAL(dp), INTENT(OUT)                              :: energy
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT), OPTIONAL    :: operators(:, :, :)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: response
      REAL(dp), INTENT(IN), OPTIONAL                     :: atom_response(:, :), cell_response(3, 3)
      COMPLEX(dp), INTENT(IN), OPTIONAL                  :: coefficients_response(:, :, :)
      TYPE(pw_r3d_rs_type), INTENT(IN), OPTIONAL         :: potential_response
      REAL(dp), INTENT(OUT), OPTIONAL                    :: core_energy

      COMPLEX(dp), ALLOCATABLE :: ao_operator(:, :), component(:, :, :), correction(:), &
         covector(:), derivative(:), gamma(:, :), matrix(:, :, :), paw_density(:, :), &
         paw_response(:, :), projected(:, :, :), projected_response(:, :, :)
      INTEGER                                            :: a, center, first, info, k, kind, last, &
                                                            nao, nk, spin
      LOGICAL                                            :: valid_shape
      REAL(dp)                                           :: core_part, inverse_cell(3, 3), &
                                                            strain(3, 3), trace_strain
      REAL(dp), ALLOCATABLE                              :: g_response(:, :), multiplicity(:), &
                                                            paw_operator(:, :), radial(:, :), &
                                                            weights(:)
      TYPE(dbcsr_p_type), POINTER                        :: ao_density(:, :)
      TYPE(pw_c1d_gs_type)                               :: field, field_response
      TYPE(qs_kind_type), POINTER                        :: kinds(:)
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_upf_inverse_augmentation_type), POINTER    :: augmentation
      TYPE(qs_upf_inverse_scf_type), POINTER             :: state
      TYPE(qs_upf_inverse_type), POINTER                 :: frame
      TYPE(upf_inverse_coulomb_center_type), ALLOCATABLE :: center_response(:), centers(:)
      TYPE(upf_inverse_coulomb_potential_type), &
         ALLOCATABLE                                     :: radial_potential(:)

      energy = 0.0_dp
      IF (PRESENT(response)) response = 0.0_dp
      IF (PRESENT(core_energy)) THEN
         CPASSERT(PRESENT(operators))
         core_energy = 0.0_dp
      END IF
      CALL get_qs_env(qs_env, ks_env=ks_env, rho=rho, qs_kind_set=kinds)
      CALL get_ks_env(ks_env, upf_inverse=state)
      IF (.NOT. ASSOCIATED(state)) RETURN
      IF (.NOT. state%enabled) RETURN
      CALL select_frame(state, 1, frame, augmentation)
      nao = SIZE(frame%coefficients, 2)
      nk = 1
      IF (ASSOCIATED(state%bloch%kpoints)) nk = SIZE(state%bloch%weights)
      ALLOCATE (matrix(nao, nao, nk), SOURCE=CMPLX(0.0_dp, 0.0_dp, dp))
      ALLOCATE (weights(nk), SOURCE=1.0_dp)
      IF (ASSOCIATED(state%bloch%kpoints)) weights(:) = state%bloch%weights
      CALL qs_rho_get(rho, rho_ao_kp=ao_density)
      DO spin = 1, SIZE(ao_density, 1)
         IF (ASSOCIATED(state%bloch%kpoints)) THEN
            CALL qs_upf_inverse_kpoints_density(state%bloch, spin, component)
            matrix(:, :, :) = matrix + component
         ELSE
            CALL qs_upf_inverse_matrix_gather(state%frame%ao, ao_density(spin, 1)%matrix, gamma)
            matrix(:, :, 1) = matrix(:, :, 1) + gamma
         END IF
      END DO
      CALL field%create(potential%pw_grid)
      CALL pw_transfer(potential, field)
      ALLOCATE (multiplicity(SIZE(field%array)), SOURCE=1.0_dp)
      IF (field%pw_grid%grid_span == HALFSPACE) THEN
         multiplicity = 2.0_dp
         IF (field%pw_grid%have_g0) multiplicity(1) = 1.0_dp
      END IF
      ALLOCATE (covector, SOURCE=field%pw_grid%vol*multiplicity*field%array)
      ALLOCATE (correction(SIZE(covector)), SOURCE=CMPLX(0.0_dp, 0.0_dp, dp))
      IF (PRESENT(operators)) THEN
         ALLOCATE (operators(nao, nao, nk), SOURCE=CMPLX(0.0_dp, 0.0_dp, dp))
         ALLOCATE (ao_operator(nao, nao))
      END IF
      IF (PRESENT(response)) THEN
         CPASSERT(PRESENT(atom_response))
         CPASSERT(PRESENT(cell_response))
         CPASSERT(PRESENT(coefficients_response))
         valid_shape = ALL(SHAPE(coefficients_response) == [SIZE(frame%coefficients, 1), nao, nk])
         CPASSERT(valid_shape)
         inverse_cell = frame%cell
         CALL invmat(inverse_cell, info)
         CPASSERT(info == 0)
         strain = MATMUL(cell_response, inverse_cell)
         trace_strain = strain(1, 1) + strain(2, 2) + strain(3, 3)
         ALLOCATE (g_response, SOURCE=-MATMUL(TRANSPOSE(strain), field%pw_grid%g))
         ALLOCATE (derivative, MOLD=correction)
         derivative = CMPLX(0.0_dp, 0.0_dp, dp)
      END IF
      IF (ASSOCIATED(state%coulomb)) THEN
         ASSOCIATE (map => state%coulomb)
            CPASSERT(ASSOCIATED(field%pw_grid, map%grid))
            ALLOCATE (centers, SOURCE=map%centers)
            DO a = 1, SIZE(centers)
               centers(a)%density = 0.0_dp
               ALLOCATE (radial, MOLD=centers(a)%density)
               kind = map%center_kinds(a)
               DO k = 1, nk
                  CALL select_frame(state, k, frame, augmentation)
                  CALL qs_upf_inverse_density_gapw(frame, augmentation, map%inverse_centers(a), kinds(kind), &
                                                   map%pairs(a), matrix(:, :, k), radial, &
                                                   atom_index=map%center_atoms(a), cache=map%fields(a, k))
                  centers(a)%density(:, :) = centers(a)%density + weights(k)*radial
               END DO
               DEALLOCATE (radial)
            END DO
            IF (PRESENT(response)) THEN
               CALL prepare_pair_responses(qs_env, state, map)
               CALL center_density_response(qs_env, state, map, matrix, atom_response, cell_response, &
                                            center_response, coefficients_response)
               CALL upf_inverse_coulomb_density(centers, map%cell, map%grid%g, map%grid%vol, correction, map%maps, &
                                                derivative, center_response, g_response, trace_strain*map%grid%vol, &
                                                cell_response=cell_response)
            ELSE
               CALL upf_inverse_coulomb_density(centers, map%cell, map%grid%g, map%grid%vol, correction, map%maps)
            END IF
            IF (PRESENT(operators)) THEN
               IF (map%grid%para%mode == PW_MODE_LOCAL) THEN
                  CALL upf_inverse_coulomb_density_adjoint(centers, map%cell, map%grid%g, map%grid%vol, &
                                                           covector, radial_potential, map%maps)
               ELSE
                  CALL upf_inverse_coulomb_density_adjoint(centers, map%cell, map%grid%g, map%grid%vol, &
                                                           covector, radial_potential, map%maps, map%grid%para%group%mp_comm_type)
               END IF
               DO k = 1, nk
                  CALL select_frame(state, k, frame, augmentation)
                  DO a = 1, SIZE(centers)
                     kind = map%center_kinds(a)
                     CALL qs_upf_inverse_density_gapw_integrate(frame, augmentation, map%inverse_centers(a), kinds(kind), &
                                                                map%pairs(a), radial_potential(a)%values, ao_operator, &
                                                                map%center_atoms(a), map%fields(a, k))
                     operators(:, :, k) = operators(:, :, k) + ao_operator
                  END DO
               END DO
            END IF
         END ASSOCIATE
      END IF
      CALL select_frame(state, 1, frame, augmentation)
      ALLOCATE (paw_density(SIZE(covector), 1))
      IF (PRESENT(response)) ALLOCATE (paw_response, MOLD=paw_density)
      DO center = 1, SIZE(frame%center_kinds)
         kind = frame%center_kinds(center)
         IF (.NOT. state%paw(kind)%active) CYCLE
         first = frame%projector_first(center)
         last = frame%projector_last(center)
         ALLOCATE (projected(last - first + 1, last - first + 1, 1), SOURCE=CMPLX(0.0_dp, 0.0_dp, dp))
         IF (PRESENT(response)) ALLOCATE (projected_response, MOLD=projected)
         IF (PRESENT(response)) projected_response = CMPLX(0.0_dp, 0.0_dp, dp)
         DO k = 1, nk
            CALL select_frame(state, k, frame, augmentation)
            projected(:, :, 1) = projected(:, :, 1) + weights(k)*MATMUL(frame%coefficients(first:last, :), &
                                                       MATMUL(matrix(:, :, k), CONJG(TRANSPOSE(frame%coefficients(first:last, :)))))
            IF (PRESENT(response)) THEN
               projected_response(:, :, 1) = projected_response(:, :, 1) + weights(k)* &
                                             (MATMUL(coefficients_response(first:last, :, k), &
                                                   MATMUL(matrix(:, :, k), CONJG(TRANSPOSE(frame%coefficients(first:last, :))))) + &
                                              MATMUL(frame%coefficients(first:last, :), &
                                                MATMUL(matrix(:, :, k), CONJG(TRANSPOSE(coefficients_response(first:last, :, k))))))
            END IF
         END DO
         IF (PRESENT(response)) THEN
            CALL qs_upf_inverse_paw_density_fourier(state%paw(kind), projected, field%pw_grid%g, &
                                     frame%positions(:, frame%center_atoms(center)), field%pw_grid%vol, paw_density, paw_response, &
                     projected_response, g_response, atom_response(:, frame%center_atoms(center)), trace_strain*field%pw_grid%vol, &
                                                    pw_grid=field%pw_grid, cell_response=cell_response)
            derivative(:) = derivative + paw_response(:, 1)
            DEALLOCATE (projected_response)
         ELSE
            CALL qs_upf_inverse_paw_density_fourier(state%paw(kind), projected, field%pw_grid%g, &
                                                   frame%positions(:, frame%center_atoms(center)), field%pw_grid%vol, paw_density, &
                                                    pw_grid=field%pw_grid)
         END IF
         correction(:) = correction + paw_density(:, 1)
         IF (PRESENT(operators)) THEN
            ALLOCATE (paw_operator(last - first + 1, last - first + 1))
            IF (field%pw_grid%para%mode == PW_MODE_LOCAL) THEN
               CALL qs_upf_inverse_paw_density_adjoint(state%paw(kind), field%pw_grid%g, &
                             frame%positions(:, frame%center_atoms(center)), field%pw_grid%vol, covector, paw_operator, core_part, &
                                                       pw_grid=field%pw_grid)
            ELSE
               CALL qs_upf_inverse_paw_density_adjoint(state%paw(kind), field%pw_grid%g, &
                             frame%positions(:, frame%center_atoms(center)), field%pw_grid%vol, covector, paw_operator, core_part, &
                                                       field%pw_grid%para%group%mp_comm_type, pw_grid=field%pw_grid)
            END IF
            IF (PRESENT(core_energy)) core_energy = core_energy + core_part
            DO k = 1, nk
               CALL select_frame(state, k, frame, augmentation)
               operators(:, :, k) = operators(:, :, k) + MATMUL(CONJG(TRANSPOSE(frame%coefficients(first:last, :))), &
                                                                MATMUL(paw_operator, frame%coefficients(first:last, :)))
            END DO
            DEALLOCATE (paw_operator)
         END IF
         DEALLOCATE (projected)
      END DO
      energy = REAL(SUM(CONJG(correction)*covector), dp)
      IF (PRESENT(response)) THEN
         response = REAL(SUM(CONJG(derivative)*covector), dp) + trace_strain*energy
         IF (PRESENT(potential_response)) THEN
            CALL field_response%create(field%pw_grid)
            CALL pw_transfer(potential_response, field_response)
            response = response + field%pw_grid%vol*REAL(SUM(CONJG(correction)*multiplicity*field_response%array), dp)
            CALL field_response%release()
         END IF
      END IF
      IF (field%pw_grid%para%mode /= PW_MODE_LOCAL) THEN
         CALL field%pw_grid%para%group%mp_comm_type%sum(energy)
         IF (PRESENT(response)) CALL field%pw_grid%para%group%mp_comm_type%sum(response)
      END IF
      CALL field%release()
   END SUBROUTINE qs_upf_inverse_local_potential

! **************************************************************************************************
!> \brief Couple a scalar field to the complete original-UPF electronic density.
!>        Add the native global source and its transpose to the signed-center and PAW correction.
!>        The frozen AE core is included once. Ionic counterterms remain with the caller.
!>        Geometry is differentiated at fixed AO occupations and fixed reduced Bloch vectors.
!> \param qs_env Current density, defining inverse coefficients and native Gaussian task lists.
!> \param potential Scalar field on the native auxiliary grid, without volume weights.
!> \param energy Complete electronic field energy, including the frozen PAW core.
!> \param operators Optional unweighted spin-independent AO operators in global k order.
!> \param response Optional joint atomic and cell derivative at fixed AO occupations.
!> \param atom_response Physical atomic motion, required with response.
!> \param cell_response Physical cell motion, required with response.
!> \param coefficients_response Defining inverse coefficient motion, required with response.
!> \param potential_response Optional derivative of the field sampled on the moving grid.
!> \param core_energy Optional frozen-core field energy, supplied with operators.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_local_field(qs_env, potential, energy, operators, response, atom_response, &
                                         cell_response, coefficients_response, potential_response, core_energy)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: potential
      REAL(dp), INTENT(OUT)                              :: energy
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT), OPTIONAL    :: operators(:, :, :)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: response
      REAL(dp), INTENT(IN), OPTIONAL                     :: atom_response(:, :), cell_response(3, 3)
      COMPLEX(dp), INTENT(IN), OPTIONAL                  :: coefficients_response(:, :, :)
      TYPE(pw_r3d_rs_type), INTENT(IN), OPTIONAL         :: potential_response
      REAL(dp), INTENT(OUT), OPTIONAL                    :: core_energy

      COMPLEX(dp), ALLOCATABLE                           :: density(:, :, :), &
                                                            global_operator(:, :, :), &
                                                            matrix(:, :), metric_response(:, :), &
                                                            part(:, :, :)
      INTEGER                                            :: axis, info, k, nao, nk, spin
      REAL(dp)                                           :: background_energy, inverse_cell(3, 3), &
                                                            strain(3, 3), trace_strain, &
                                                            virial(3, 3)
      REAL(dp), ALLOCATABLE                              :: forces(:, :), relative_motion(:, :), &
                                                            weights(:)
      TYPE(dbcsr_p_type), POINTER                        :: ao_density(:, :)
      TYPE(pw_r3d_rs_type)                               :: field_response, zero, zeros(3)
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_upf_inverse_augmentation_type), POINTER    :: augmentation
      TYPE(qs_upf_inverse_background_type), POINTER      :: background
      TYPE(qs_upf_inverse_grid_charge_type)              :: charge
      TYPE(qs_upf_inverse_scf_type), POINTER             :: state
      TYPE(qs_upf_inverse_type), POINTER                 :: frame

      CALL qs_upf_inverse_local_potential(qs_env, potential, energy, operators, response, atom_response, &
                                          cell_response, coefficients_response, potential_response, core_energy)
      CALL get_qs_env(qs_env, ks_env=ks_env, rho=rho)
      CALL get_ks_env(ks_env, upf_inverse=state)
      IF (.NOT. ASSOCIATED(state)) RETURN
      IF (.NOT. state%enabled) RETURN
      CALL select_frame(state, 1, frame, augmentation)
      nao = SIZE(frame%coefficients, 2)
      nk = 1
      IF (ASSOCIATED(state%bloch%kpoints)) nk = SIZE(state%bloch%weights)
      ALLOCATE (density(nao, nao, nk), global_operator(nao, nao, nk), weights(nk))
      density = CMPLX(0.0_dp, 0.0_dp, dp)
      global_operator = CMPLX(0.0_dp, 0.0_dp, dp)
      weights = 1.0_dp
      IF (ASSOCIATED(state%bloch%kpoints)) weights(:) = state%bloch%weights
      CALL qs_rho_get(rho, rho_ao_kp=ao_density)
      DO spin = 1, SIZE(ao_density, 1)
         IF (ASSOCIATED(state%bloch%kpoints)) THEN
            CALL qs_upf_inverse_kpoints_density(state%bloch, spin, part)
            density(:, :, :) = density + part
         ELSE
            CALL qs_upf_inverse_matrix_gather(frame, ao_density(spin, 1)%matrix, matrix)
            density(:, :, 1) = density(:, :, 1) + matrix
         END IF
      END DO
      NULLIFY (background)
      IF (state%gaussian_background) THEN
         background => state%coulomb%background
      ELSE IF (state%native_gaussian) THEN
         background => state%gaussian
      END IF
      background_energy = 0.0_dp
      IF (ASSOCIATED(background)) THEN
         CALL qs_upf_inverse_background_integrate(background, potential=potential, &
                                                  ao_operators=global_operator, ks_env=ks_env)
         DO k = 1, nk
            background_energy = background_energy + weights(k)* &
                                REAL(SUM(CONJG(global_operator(:, :, k))*density(:, :, k)), dp)
         END DO
      END IF
      CALL zero%create(potential%pw_grid)
      zero%array = 0.0_dp
      DO axis = 1, 3
         CALL zeros(axis)%create(potential%pw_grid)
         zeros(axis)%array = 0.0_dp
      END DO
      IF (.NOT. state%gaussian_background) THEN
         DO k = 1, nk
            IF (ASSOCIATED(state%bloch%kpoints)) THEN
               CALL qs_upf_inverse_grid_integrate_dense(state%bloch%frames(k), potential, zeros, zero, matrix, &
                                                        augmentation=state%bloch%augmentation, fourier_charge=.TRUE., &
                                                        orbital_fields=.NOT. state%native_gaussian)
            ELSE
               CALL qs_upf_inverse_grid_integrate_dense(state%frame, potential, zeros, zero, matrix, &
                                                        augmentation=state%augmentation, fourier_charge=.TRUE., &
                                                        orbital_fields=.NOT. state%native_gaussian)
            END IF
            global_operator(:, :, k) = global_operator(:, :, k) + matrix
         END DO
      END IF
      DO k = 1, nk
         energy = energy + weights(k)*REAL(SUM(CONJG(global_operator(:, :, k))*density(:, :, k)), dp)
      END DO
      IF (PRESENT(operators)) operators(:, :, :) = operators + global_operator
      IF (PRESENT(response)) THEN
         CALL field_response%create(potential%pw_grid)
         field_response%array = 0.0_dp
         IF (PRESENT(potential_response)) field_response%array = potential_response%array
         inverse_cell = frame%cell
         CALL invmat(inverse_cell, info)
         CPASSERT(info == 0)
         strain = MATMUL(cell_response, inverse_cell)
         trace_strain = strain(1, 1) + strain(2, 2) + strain(3, 3)
         IF (ASSOCIATED(background)) THEN
            ALLOCATE (forces(3, SIZE(atom_response, 2)), relative_motion(3, SIZE(atom_response, 2)))
            relative_motion(:, :) = atom_response - MATMUL(strain, frame%positions)
            CALL qs_upf_inverse_background_geometry(background, matrix=density, potential=potential, &
                                                    forces=forces, virial=virial, ks_env=ks_env)
            response = response - SUM(forces*relative_motion) - SUM(virial*strain) + trace_strain*background_energy
            CALL qs_upf_inverse_background_integrate(background, potential=field_response, &
                                                     ao_operators=global_operator, ks_env=ks_env)
            DO k = 1, nk
               response = response + weights(k)*REAL(SUM(CONJG(global_operator(:, :, k))*density(:, :, k)), dp)
            END DO
         END IF
         IF (.NOT. state%gaussian_background) THEN
            DO k = 1, nk
               IF (ASSOCIATED(state%bloch%kpoints)) THEN
                  CALL qs_upf_inverse_grid_response(state%bloch%frames(k), state%bloch%augmentation, &
                                                    state%bloch%frames(k)%ao%radial_q, potential, zeros, zero, &
                                                    atom_response, cell_response, [0.0_dp, 0.0_dp, 0.0_dp], &
                                                    field_response, zeros, zero, matrix, metric_response, &
                                                    orbital_fields=.NOT. state%native_gaussian, &
                                                    known_coefficients_response=coefficients_response(:, :, k), &
                                                    compute_metric=.FALSE., charge=charge)
               ELSE
                  CALL qs_upf_inverse_grid_response(state%frame, state%augmentation, frame%radial_q, potential, zeros, zero, &
                                                    atom_response, cell_response, [0.0_dp, 0.0_dp, 0.0_dp], &
                                                    field_response, zeros, zero, matrix, metric_response, &
                                                    orbital_fields=.NOT. state%native_gaussian, &
                                                    known_coefficients_response=coefficients_response(:, :, k), &
                                                    compute_metric=.FALSE.)
               END IF
               response = response + weights(k)*REAL(SUM(CONJG(matrix)*density(:, :, k)), dp)
            END DO
         END IF
         CALL field_response%release()
      END IF
      CALL zero%release()
      DO axis = 1, 3
         CALL zeros(axis)%release()
      END DO
   END SUBROUTINE qs_upf_inverse_local_field

! **************************************************************************************************
!> \brief Complete spin valence densities projected onto the existing auxiliary grid.
!>        Combine the accepted soft Gaussian density with every signed center contribution,
!>        including original Q. The SCF's soft fields and Coulomb state are not replaced.
!>        All ranks must participate. Atomic densities are replicated, Fourier rows retain
!>        their original distribution, and Bloch weights enter once. The result has the
!>        auxiliary grid's spatial resolution and excludes ions and frozen-core densities.
!> \param qs_env Native owner of the accepted occupation matrices and full-GAPW source map.
!> \param density Owned real-space spin fields. Unallocated when no inverse GAPW map is active.
!> \param monopole_exponent Optional inverse-square width for the native spherical hard-density
!>        approximation. Each complete center charge is preserved, with higher multipoles omitted.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_coulomb_density(qs_env, density, monopole_exponent)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(pw_r3d_rs_type), ALLOCATABLE, INTENT(OUT)     :: density(:)
      REAL(dp), INTENT(IN), OPTIONAL                     :: monopole_exponent

      INTEGER                                            :: spin
      TYPE(pw_c1d_gs_type), ALLOCATABLE                  :: density_g(:)

      CALL coulomb_density_fourier(qs_env, density_g, monopole_exponent)
      IF (.NOT. ALLOCATED(density_g)) RETURN
      ALLOCATE (density(SIZE(density_g)))
      DO spin = 1, SIZE(density)
         CALL density(spin)%create(density_g(spin)%pw_grid)
         CALL pw_transfer(density_g(spin), density(spin))
         CALL density_g(spin)%release()
      END DO
   END SUBROUTINE qs_upf_inverse_coulomb_density

! **************************************************************************************************
!> \brief Same complete valence source before its observable real-space FFT projection.
!> \param qs_env ...
!> \param density Owned Fourier spin fields. Unallocated when no inverse GAPW map is active.
!> \param monopole_exponent ...
! **************************************************************************************************
   SUBROUTINE coulomb_density_fourier(qs_env, density, monopole_exponent)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(pw_c1d_gs_type), ALLOCATABLE, INTENT(OUT)     :: density(:)
      REAL(dp), INTENT(IN), OPTIONAL                     :: monopole_exponent

      COMPLEX(dp), ALLOCATABLE                           :: center_g(:), gamma(:, :), matrix(:, :, :)
      INTEGER                                            :: a, k, kind, spin
      LOGICAL                                            :: rho_g_valid
      REAL(dp)                                           :: charge
      REAL(dp), ALLOCATABLE                              :: radial(:, :)
      TYPE(dbcsr_p_type), POINTER                        :: matrices(:, :)
      TYPE(pw_c1d_gs_type)                               :: physical_g
      TYPE(pw_c1d_gs_type), POINTER                      :: soft_g(:)
      TYPE(qs_kind_type), POINTER                        :: kinds(:)
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_upf_inverse_augmentation_type), POINTER    :: augmentation
      TYPE(qs_upf_inverse_scf_type), POINTER             :: state
      TYPE(qs_upf_inverse_type), POINTER                 :: frame
      TYPE(upf_inverse_coulomb_center_type), ALLOCATABLE :: centers(:)

      CALL get_qs_env(qs_env, ks_env=ks_env, rho=rho, qs_kind_set=kinds)
      IF (.NOT. ASSOCIATED(ks_env)) RETURN
      CALL get_ks_env(ks_env, upf_inverse=state)
      IF (.NOT. ASSOCIATED(state)) RETURN
      IF (.NOT. state%enabled .OR. .NOT. state%gaussian_background) RETURN
      CPASSERT(ASSOCIATED(state%coulomb))
      CPASSERT(state%coulomb%ready)
      CALL qs_rho_get(rho, rho_ao_kp=matrices, rho_g=soft_g, rho_g_valid=rho_g_valid)
      CPASSERT(rho_g_valid)
      CPASSERT(SIZE(soft_g) == SIZE(matrices, 1))
      ALLOCATE (density(SIZE(soft_g)))
      ALLOCATE (centers, SOURCE=state%coulomb%centers)
      CALL physical_g%create(state%coulomb%grid)
      ALLOCATE (center_g(SIZE(physical_g%array)))
      DO spin = 1, SIZE(soft_g)
         IF (ASSOCIATED(state%bloch%kpoints)) THEN
            CALL qs_upf_inverse_kpoints_density(state%bloch, spin, matrix)
         ELSE
            CALL qs_upf_inverse_matrix_gather(state%frame%ao, matrices(spin, 1)%matrix, gamma)
            ALLOCATE (matrix(SIZE(gamma, 1), SIZE(gamma, 2), 1))
            matrix(:, :, 1) = gamma
            DEALLOCATE (gamma)
         END IF
         DO a = 1, SIZE(centers)
            ALLOCATE (radial, MOLD=centers(a)%density)
            centers(a)%density = 0.0_dp
            kind = state%coulomb%center_kinds(a)
            DO k = 1, SIZE(matrix, 3)
               CALL select_frame(state, k, frame, augmentation)
               CALL qs_upf_inverse_density_gapw(frame, augmentation, state%coulomb%inverse_centers(a), kinds(kind), &
                                                state%coulomb%pairs(a), matrix(:, :, k), radial, &
                                                atom_index=state%coulomb%center_atoms(a), cache=state%coulomb%fields(a, k))
               centers(a)%density(:, :) = centers(a)%density + state%coulomb%weights(k)*radial
            END DO
            DEALLOCATE (radial)
         END DO
         IF (PRESENT(monopole_exponent)) THEN
            CPASSERT(monopole_exponent > 0.0_dp)
            center_g = CMPLX(0.0_dp, 0.0_dp, dp)
            DO a = 1, SIZE(centers)
               charge = SQRT(fourpi)*SUM(centers(a)%weights*centers(a)%density(:, 1))
               center_g(:) = center_g + charge/state%coulomb%grid%vol* &
                             EXP(CMPLX(-SUM(state%coulomb%grid%g**2, DIM=1)/(4.0_dp*monopole_exponent), &
                                       -MATMUL(TRANSPOSE(state%coulomb%grid%g), centers(a)%position), dp))
            END DO
         ELSE
            CALL upf_inverse_coulomb_density(centers, state%coulomb%cell, state%coulomb%grid%g, &
                                             state%coulomb%grid%vol, center_g, state%coulomb%maps)
         END IF
         physical_g%array = soft_g(spin)%array + center_g
         CALL density(spin)%create(state%coulomb%grid)
         density(spin)%array = physical_g%array
         DEALLOCATE (matrix)
      END DO
      CALL physical_g%release()
   END SUBROUTINE coulomb_density_fourier

! **************************************************************************************************
!> \brief Reconstructed PAW AE density for the existing X-ray Fourier convention.
!>        Combine the accepted compensated valence source with the original AE-minus-PS-minus-Q
!>        products and frozen AE core. Neither the accepted SCF density nor its Coulomb map changes.
!> \param qs_env Native owner of current occupations, inverse frames and original PAW data.
!> \param density Owned by the caller; returned coefficients are integrals, without 1/volume.
!> \param q_max Maximum magnitude of the requested reciprocal vectors.
!> \param rho_hard Signed electronic charge added to the compensated valence source.
!> \param rho_soft Signed charge of the complete compensated valence source.
!> \param active False when no inverse UPF state is active, leaving the legacy caller in charge.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_paw_density(qs_env, density, q_max, rho_hard, rho_soft, active)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(pw_c1d_gs_type), INTENT(INOUT)                :: density
      REAL(dp), INTENT(IN)                               :: q_max
      REAL(dp), INTENT(OUT)                              :: rho_hard, rho_soft
      LOGICAL, INTENT(OUT)                               :: active

      INTEGER                                            :: center, spin
      TYPE(pw_c1d_gs_type), ALLOCATABLE                  :: complete(:)
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_upf_inverse_augmentation_type), POINTER    :: augmentation
      TYPE(qs_upf_inverse_scf_type), POINTER             :: state
      TYPE(qs_upf_inverse_type), POINTER                 :: frame

      active = .FALSE.
      rho_hard = 0.0_dp
      rho_soft = 0.0_dp
      CALL get_qs_env(qs_env, ks_env=ks_env)
      IF (.NOT. ASSOCIATED(ks_env)) RETURN
      CALL get_ks_env(ks_env, upf_inverse=state)
      IF (.NOT. ASSOCIATED(state)) RETURN
      IF (.NOT. state%enabled) RETURN
      active = .TRUE.
      CALL select_frame(state, 1, frame, augmentation)
      DO center = 1, SIZE(frame%center_kinds)
         IF (.NOT. state%paw(frame%center_kinds(center))%active) THEN
            CPABORT("All-electron X-ray density for inverse UPF requires PAW data for every atom")
         END IF
      END DO
      CALL electronic_density_fourier(qs_env, density%pw_grid, complete, rho_soft, .FALSE., q_max)
      density%array = CMPLX(0.0_dp, 0.0_dp, dp)
      DO spin = 1, SIZE(complete)
         density%array = density%array + complete(spin)%array
         CALL complete(spin)%release()
      END DO
      rho_hard = pw_integrate_function(density, isign=-1) - rho_soft
      CALL pw_scale(density, density%pw_grid%vol)
   END SUBROUTINE qs_upf_inverse_paw_density

! **************************************************************************************************
!> \brief Complete electronic spin densities conjugate to the native scalar-potential service.
!>        Retain the original USPP augmentation and full-GAPW centers. PAW also contributes its
!>        AE-minus-PS-minus-Q density and the frozen AE core, split over physical spin channels.
!>        These are electron-number densities, not total electrostatic sources. Ionic charges
!>        and the PAW core countercharge remain separate. NLCC model cores are XC inputs only.
!>        Native grid projection is retained without normalization to the charge metric.
!> \param qs_env Native owner of the accepted occupation matrices and original atomic data.
!> \param density Owned Fourier spin fields in electron-number per volume units.
!> \param active False with unallocated density when no original-UPF state is enabled.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_electronic_density(qs_env, density, active)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(pw_c1d_gs_type), ALLOCATABLE, INTENT(OUT)     :: density(:)
      LOGICAL, INTENT(OUT)                               :: active

      REAL(dp)                                           :: soft_charge
      TYPE(pw_c1d_gs_type), POINTER                      :: soft(:)
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_upf_inverse_scf_type), POINTER             :: state

      active = .FALSE.
      CALL get_qs_env(qs_env, ks_env=ks_env, rho=rho)
      IF (.NOT. ASSOCIATED(ks_env)) RETURN
      CALL get_ks_env(ks_env, upf_inverse=state)
      IF (.NOT. ASSOCIATED(state)) RETURN
      IF (.NOT. state%enabled) RETURN
      CALL qs_rho_get(rho, rho_g=soft)
      CPASSERT(ASSOCIATED(soft))
      CALL electronic_density_fourier(qs_env, soft(1)%pw_grid, density, soft_charge, .TRUE.)
      active = .TRUE.
   END SUBROUTINE qs_upf_inverse_electronic_density

! **************************************************************************************************
!> \brief Shared Fourier assembly with the requested native-grid or X-ray convention.
!> \param qs_env ...
!> \param grid ...
!> \param density ...
!> \param soft_charge Signed compensated-valence charge before the PAW AE reconstruction.
!> \param grid_convention Retain native Nyquist projection for a scalar grid-potential consumer.
!> \param q_max Optional PAW correction cutoff used by the existing X-ray observable.
! **************************************************************************************************
   SUBROUTINE electronic_density_fourier(qs_env, grid, density, soft_charge, grid_convention, q_max)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(pw_grid_type), POINTER                        :: grid
      TYPE(pw_c1d_gs_type), ALLOCATABLE, INTENT(OUT)     :: density(:)
      REAL(dp), INTENT(OUT)                              :: soft_charge
      LOGICAL, INTENT(IN)                                :: grid_convention
      REAL(dp), INTENT(IN), OPTIONAL                     :: q_max

      COMPLEX(dp), ALLOCATABLE                           :: correction(:, :), gamma(:, :), &
                                                            matrices(:, :, :, :), matrix(:, :, :), &
                                                            projected(:, :, :)
      INTEGER                                            :: center, first, k, kind, last, nao, nk, &
                                                            ns, spin
      LOGICAL                                            :: valid
      REAL(dp), ALLOCATABLE                              :: weights(:)
      TYPE(dbcsr_p_type), POINTER                        :: ao_density(:, :)
      TYPE(pw_c1d_gs_type), ALLOCATABLE                  :: complete(:)
      TYPE(pw_c1d_gs_type), POINTER                      :: soft(:)
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_upf_inverse_augmentation_type), POINTER    :: augmentation
      TYPE(qs_upf_inverse_scf_type), POINTER             :: state
      TYPE(qs_upf_inverse_type), POINTER                 :: frame

      soft_charge = 0.0_dp
      CALL get_qs_env(qs_env, ks_env=ks_env, rho=rho)
      CALL get_ks_env(ks_env, upf_inverse=state)
      CPASSERT(ASSOCIATED(state))
      CPASSERT(state%enabled)
      CALL select_frame(state, 1, frame, augmentation)
      CALL qs_rho_get(rho, rho_ao_kp=ao_density, rho_g=soft, rho_g_valid=valid)
      CPASSERT(valid .AND. ASSOCIATED(soft))
      ns = SIZE(soft)
      nao = SIZE(frame%coefficients, 2)
      nk = 1
      IF (ASSOCIATED(state%bloch%kpoints)) nk = SIZE(state%bloch%weights)
      CALL coulomb_density_fourier(qs_env, complete)
      ALLOCATE (density(ns))
      DO spin = 1, ns
         CALL density(spin)%create(grid)
         IF (ALLOCATED(complete)) THEN
            CALL pw_transfer(complete(spin), density(spin))
            CALL complete(spin)%release()
         ELSE
            CALL pw_transfer(soft(spin), density(spin))
         END IF
         soft_charge = soft_charge + pw_integrate_function(density(spin), isign=-1)
      END DO
      IF (.NOT. ANY(state%paw%active)) RETURN
      ALLOCATE (matrices(nao, nao, nk, ns), weights(nk))
      weights = 1.0_dp
      IF (ASSOCIATED(state%bloch%kpoints)) weights(:) = state%bloch%weights
      DO spin = 1, ns
         IF (ASSOCIATED(state%bloch%kpoints)) THEN
            CALL qs_upf_inverse_kpoints_density(state%bloch, spin, matrix)
            matrices(:, :, :, spin) = matrix
            DEALLOCATE (matrix)
         ELSE
            CALL qs_upf_inverse_matrix_gather(state%frame%ao, ao_density(spin, 1)%matrix, gamma)
            matrices(:, :, 1, spin) = gamma
            DEALLOCATE (gamma)
         END IF
      END DO
      ALLOCATE (correction(SIZE(density(1)%array), ns))
      DO center = 1, SIZE(frame%center_kinds)
         kind = frame%center_kinds(center)
         IF (.NOT. state%paw(kind)%active) CYCLE
         first = frame%projector_first(center)
         last = frame%projector_last(center)
         ALLOCATE (projected(last - first + 1, last - first + 1, ns), SOURCE=CMPLX(0.0_dp, 0.0_dp, dp))
         DO k = 1, nk
            CALL select_frame(state, k, frame, augmentation)
            DO spin = 1, ns
               projected(:, :, spin) = projected(:, :, spin) + weights(k)* &
                                       MATMUL(frame%coefficients(first:last, :), &
                                              MATMUL(matrices(:, :, k, spin), &
                                                     CONJG(TRANSPOSE(frame%coefficients(first:last, :)))))
            END DO
         END DO
         IF (grid_convention) THEN
            CALL qs_upf_inverse_paw_density_fourier(state%paw(kind), projected, grid%g, &
                                                    frame%positions(:, frame%center_atoms(center)), grid%vol, &
                                                    correction, pw_grid=grid)
         ELSE
            CALL qs_upf_inverse_paw_density_fourier(state%paw(kind), projected, grid%g, &
                                                    frame%positions(:, frame%center_atoms(center)), grid%vol, correction)
         END IF
         DO spin = 1, ns
            IF (PRESENT(q_max)) THEN
               WHERE (grid%gsq <= q_max*q_max) density(spin)%array = density(spin)%array + correction(:, spin)
            ELSE
               density(spin)%array = density(spin)%array + correction(:, spin)
            END IF
         END DO
         DEALLOCATE (projected)
      END DO
   END SUBROUTINE electronic_density_fourier

! **************************************************************************************************
!> \brief Complete electronic metric counts and the original ionic Gaussian grid charge.
!>        All ranks participate in the ionic zero-mode reduction. No density field is rebuilt.
!> \param qs_env Native owner of the accepted inverse GAPW source.
!> \param electronic Positive spin electron counts. Unallocated when this source is inactive.
!> \param ionic Positive ionic Gaussian charge, excluding any NLCC model-core density.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_coulomb_charges(qs_env, electronic, ionic)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: electronic(:)
      REAL(dp), INTENT(OUT)                              :: ionic

      REAL(dp), POINTER                                  :: total(:)
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_upf_inverse_scf_type), POINTER             :: state

      ionic = 0.0_dp
      CALL get_qs_env(qs_env, ks_env=ks_env, rho=rho)
      IF (.NOT. ASSOCIATED(ks_env)) RETURN
      CALL get_ks_env(ks_env, upf_inverse=state)
      IF (.NOT. ASSOCIATED(state)) RETURN
      IF (.NOT. state%enabled .OR. .NOT. state%gaussian_background) RETURN
      CPASSERT(ASSOCIATED(state%coulomb))
      CPASSERT(state%coulomb%ready)
      CALL qs_rho_get(rho, tot_rho_r=total)
      ALLOCATE (electronic, SOURCE=-total)
      IF (state%coulomb%grid%have_g0) ionic = -state%coulomb%grid%vol*REAL(state%coulomb%ionic%array(1), dp)
      IF (state%coulomb%grid%para%mode /= PW_MODE_LOCAL) THEN
         CALL state%coulomb%grid%para%group%mp_comm_type%sum(ionic)
      END IF
   END SUBROUTINE qs_upf_inverse_coulomb_charges

! **************************************************************************************************
!> \brief Geometry adjoint of every native Gaussian ionic charge on the actual auxiliary grid.
!>        Differentiates calculate_rho_core with its collocation extent and kind selection.
!>        Native force, virial and atomic-energy owners are unchanged. Grid measure and Fourier
!>        kernel geometry remain the caller's responsibility, as for the electronic background.
!> \param qs_env Native owner of the same ionic charges and auxiliary grid.
!> \param potential Unweighted real-space covector conjugate to the negative ionic density.
!> \param forces Replicated physical atomic forces from the ionic grid source.
!> \param virial Replicated source virial, including the current grid measure.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_coulomb_ionic_geometry(qs_env, potential, forces, virial)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: potential
      REAL(dp), INTENT(OUT)                              :: forces(:, :), virial(3, 3)

      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(pw_r3d_rs_type)                               :: weighted

      CALL get_qs_env(qs_env, para_env=para_env)
      CALL weighted%create(potential%pw_grid)
      weighted%array = potential%pw_grid%dvol*potential%array
      CALL integrate_v_core_rspace(weighted, qs_env, force_response=forces, &
                                   virial_response=virial, only_nopaw=.FALSE.)
      CALL para_env%sum(forces)
      CALL para_env%sum(virial)
      forces(:, :) = -forces
      CALL weighted%release()
   END SUBROUTINE qs_upf_inverse_coulomb_ionic_geometry

! **************************************************************************************************
!> \brief Analytic Fourier coefficients of the periodic Gaussian ionic charges and their adjoint.
!>        Each term represents a physical wave vector folded into the real FFT coefficient.
!> \param grid Auxiliary reciprocal grid.
!> \param index Local reciprocal coefficient index.
!> \param position Cartesian atomic center.
!> \param alpha Gaussian exponent.
!> \param charge Integrated Gaussian charge divided by the cell volume.
!> \param vectors Physical wave vectors.
!> \param terms Fourier coefficients including the averaging weights.
!> \param count Number of physical wave vectors for this coefficient.
! **************************************************************************************************
   SUBROUTINE ionic_fourier_terms(grid, index, position, alpha, charge, vectors, terms, count)
      TYPE(pw_grid_type), INTENT(IN)                     :: grid
      INTEGER, INTENT(IN)                                :: index
      REAL(dp), INTENT(IN)                               :: position(3), alpha, charge
      REAL(dp), INTENT(OUT)                              :: vectors(3, 8)
      COMPLEX(dp), INTENT(OUT)                           :: terms(8)
      INTEGER, INTENT(OUT)                               :: count

      INTEGER                                            :: j
      REAL(dp)                                           :: g(3)

      CALL pw_real_fourier_aliases(grid, index, vectors, count)
      DO j = 1, count
         g = vectors(:, j)
         terms(j) = charge/REAL(count, dp)*EXP(CMPLX(-DOT_PRODUCT(g, g)/(4.0_dp*alpha), -DOT_PRODUCT(g, position), dp))
      END DO
   END SUBROUTINE ionic_fourier_terms

! **************************************************************************************************
!> \brief Construct the periodic ionic source without sampling narrow Gaussians in real space.
!> \param atoms Native atoms at the prepared geometry.
!> \param kinds Native Gaussian charge definitions.
!> \param density Fourier source including the negative ionic charge.
! **************************************************************************************************
   SUBROUTINE ionic_fourier_source(atoms, kinds, density)
      TYPE(particle_type), INTENT(IN)                    :: atoms(:)
      TYPE(qs_kind_type), INTENT(IN)                     :: kinds(:)
      TYPE(pw_c1d_gs_type), INTENT(INOUT)                :: density

      COMPLEX(dp)                                        :: terms(8)
      INTEGER                                            :: atom, count, i, kind
      REAL(dp)                                           :: alpha, charge, coefficient, vectors(3, 8)
      TYPE(pw_grid_type), POINTER                        :: grid

      grid => density%pw_grid
      density%array = 0.0_dp
      DO atom = 1, SIZE(atoms)
         kind = atoms(atom)%atomic_kind%kind_number
         CALL get_qs_kind(kinds(kind), alpha_core_charge=alpha, ccore_charge=coefficient)
         IF (alpha == 0.0_dp .OR. coefficient == 0.0_dp) CYCLE
         charge = -coefficient*(pi/alpha)**1.5_dp/grid%vol
         DO i = 1, SIZE(density%array)
            CALL ionic_fourier_terms(grid, i, atoms(atom)%r, alpha, charge, vectors, terms, count)
            density%array(i) = density%array(i) + SUM(terms(:count))
         END DO
      END DO
   END SUBROUTINE ionic_fourier_source

! **************************************************************************************************
!> \brief Differentiate the same ionic Fourier coefficients at fixed fractional centers for strain.
!> \param qs_env Owner of the atomic positions and charge definitions.
!> \param map Prepared reciprocal potential and FFT weights.
!> \param forces Negative derivative with respect to Cartesian displacements.
!> \param virial Negative derivative with respect to affine strain, excluding the integral measure.
! **************************************************************************************************
   SUBROUTINE ionic_fourier_geometry(qs_env, map, forces, virial)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_upf_inverse_coulomb_type), INTENT(IN)      :: map
      REAL(dp), INTENT(OUT)                              :: forces(:, :), virial(3, 3)

      COMPLEX(dp)                                        :: covector, term, terms(8)
      INTEGER                                            :: atom, axis, count, i, j, kind
      REAL(dp)                                           :: alpha, charge, coefficient, g(3), &
                                                            vectors(3, 8)
      TYPE(particle_type), POINTER                       :: atoms(:)
      TYPE(pw_grid_type), POINTER                        :: grid
      TYPE(qs_kind_type), POINTER                        :: kinds(:)

      CALL get_qs_env(qs_env, qs_kind_set=kinds, particle_set=atoms)
      grid => map%grid
      forces = 0.0_dp
      virial = 0.0_dp
      DO atom = 1, SIZE(atoms)
         kind = atoms(atom)%atomic_kind%kind_number
         CALL get_qs_kind(kinds(kind), alpha_core_charge=alpha, ccore_charge=coefficient)
         IF (alpha == 0.0_dp .OR. coefficient == 0.0_dp) CYCLE
         charge = -coefficient*(pi/alpha)**1.5_dp/grid%vol
         DO i = 1, SIZE(map%potential_g%array)
            covector = grid%vol*map%multiplicity(i)*CONJG(map%potential_g%array(i))
            CALL ionic_fourier_terms(grid, i, atoms(atom)%r, alpha, charge, vectors, terms, count)
            DO j = 1, count
               g = vectors(:, j)
               term = covector*terms(j)
               forces(:, atom) = forces(:, atom) - AIMAG(term)*g
               DO axis = 1, 3
                  virial(:, axis) = virial(:, axis) - REAL(term, dp)*g*g(axis)/(2.0_dp*alpha)
                  virial(axis, axis) = virial(axis, axis) + REAL(term, dp)
               END DO
            END DO
         END DO
      END DO
      IF (grid%para%mode /= PW_MODE_LOCAL) THEN
         CALL grid%para%group%mp_comm_type%sum(forces)
         CALL grid%para%group%mp_comm_type%sum(virial)
      END IF
   END SUBROUTINE ionic_fourier_geometry

! **************************************************************************************************
!> \brief Own the complete Coulomb functional and install its metric before the first SCF density.
!>        The caller must connect the same source and AO operator to native energy accounting.
!>        Native GAPW energy initialization attaches this source before allocating new KS matrices.
!> \param qs_env Fresh native inverse environment at the current geometry.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_coulomb_attach(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_upf_inverse_scf_type), POINTER             :: state

      CALL get_qs_env(qs_env, ks_env=ks_env)
      CALL get_ks_env(ks_env, upf_inverse=state)
      CPASSERT(ASSOCIATED(state))
      CPASSERT(state%density_evaluations == 0 .AND. state%operator_evaluations == 0)
      IF (.NOT. ASSOCIATED(state%coulomb)) ALLOCATE (state%coulomb)
      CALL qs_upf_inverse_coulomb_prepare(qs_env, state, state%coulomb)
      CALL qs_upf_inverse_coulomb_install_metric(qs_env, state%coulomb)
      state%gaussian_background = .TRUE.
      IF (ALLOCATED(state%hartree_k)) DEALLOCATE (state%hartree_k)
      ALLOCATE (state%hartree_k, MOLD=state%coulomb%metric)
      state%hartree_k(:, :, :) = CMPLX(0.0_dp, 0.0_dp, dp)
      state%coulomb_density_evaluation = -1
      state%coulomb_operator_evaluation = -1
   END SUBROUTINE qs_upf_inverse_coulomb_attach

! **************************************************************************************************
!> \brief Return the folded XC background potential through the same native Gaussian source map.
!>        The atomic XC adjoint is added separately, and the common Hartree operator only afterward.
!> \param qs_env Current native owners of the Gaussian source and inverse state.
!> \param potential Folded XC grid potential, including the native volume element if weighted.
!> \param hamiltonian Native spin matrix to receive the matching XC contribution.
!> \param weighted Whether the native grid measure is already included.
!> \param spin Native spin index.
!> \param kinetic Fold the soft kinetic-density potential into the existing Hamiltonian.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_coulomb_xc_integrate(qs_env, potential, hamiltonian, weighted, spin, kinetic)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: potential
      TYPE(dbcsr_type), INTENT(INOUT)                    :: hamiltonian
      LOGICAL, INTENT(IN)                                :: weighted
      INTEGER, INTENT(IN)                                :: spin
      LOGICAL, INTENT(IN), OPTIONAL                      :: kinetic

      COMPLEX(dp), ALLOCATABLE                           :: dense(:, :, :)
      LOGICAL                                            :: do_kinetic
      TYPE(pw_r3d_rs_type)                               :: unweighted
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_upf_inverse_scf_type), POINTER             :: state

      CALL get_qs_env(qs_env, ks_env=ks_env)
      CALL get_ks_env(ks_env, upf_inverse=state)
      CPASSERT(ASSOCIATED(state))
      CPASSERT(state%gaussian_background .AND. ASSOCIATED(state%coulomb))
      CPASSERT(state%coulomb%ready)
      do_kinetic = .FALSE.
      IF (PRESENT(kinetic)) do_kinetic = kinetic
      CALL unweighted%create(potential%pw_grid)
      unweighted%array = potential%array
      IF (weighted) CALL pw_scale(unweighted, 1.0_dp/potential%pw_grid%dvol)
      ALLOCATE (dense, MOLD=state%coulomb%metric)
      CALL qs_upf_inverse_background_integrate(state%coulomb%background, qs_env, unweighted, dense, do_kinetic)
      CALL qs_upf_inverse_scf_integrate_background(ks_env, dense, hamiltonian, spin, additive=do_kinetic)
      CALL unweighted%release()
   END SUBROUTINE qs_upf_inverse_coulomb_xc_integrate

! **************************************************************************************************
!> \brief Evaluate native Coulomb energy and store its common AO operator from accepted spin densities.
!>        Spin densities are summed before applying mesh weights once in the common source map.
!>        The per-k AO operators have no mesh weight and are shared by all spins.
!> \param qs_env Current native inverse environment and accepted density state.
!> \param energy Complete electronic/ionic Gaussian Coulomb energy. Native ionic corrections stay separate.
!> \param charge Optional complete electronic charge, excluding ions and frozen partial cores.
!> \param electronic_energy Optional electron-only Hartree energy from the same quadratic functional.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_coulomb_update(qs_env, energy, charge, electronic_energy)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(dp), INTENT(OUT)                              :: energy
      REAL(dp), INTENT(OUT), OPTIONAL                    :: charge, electronic_energy

      COMPLEX(dp), ALLOCATABLE                           :: component(:, :, :), gamma(:, :), &
                                                            matrix(:, :, :), operators(:, :, :)
      INTEGER                                            :: k, spin
      REAL(dp)                                           :: ionic_energy
      TYPE(dbcsr_p_type), POINTER                        :: density(:, :)
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_upf_inverse_scf_type), POINTER             :: state

      CALL get_qs_env(qs_env, ks_env=ks_env, rho=rho)
      CALL get_ks_env(ks_env, upf_inverse=state)
      CPASSERT(ASSOCIATED(state))
      CPASSERT(ASSOCIATED(state%coulomb))
      CPASSERT(state%coulomb%ready)
      CALL qs_rho_get(rho, rho_ao_kp=density)
      ALLOCATE (matrix, MOLD=state%coulomb%metric)
      matrix(:, :, :) = CMPLX(0.0_dp, 0.0_dp, dp)
      DO spin = 1, SIZE(density, 1)
         IF (ASSOCIATED(state%bloch%kpoints)) THEN
            CALL qs_upf_inverse_kpoints_density(state%bloch, spin, component)
            matrix(:, :, :) = matrix + component
         ELSE
            CPASSERT(SIZE(density, 2) == 1)
            CALL qs_upf_inverse_matrix_gather(state%frame%ao, density(spin, 1)%matrix, gamma)
            matrix(:, :, 1) = matrix(:, :, 1) + gamma
         END IF
      END DO
      IF (.NOT. ALLOCATED(state%hartree_k)) ALLOCATE (state%hartree_k, MOLD=matrix)
      ALLOCATE (operators, MOLD=matrix)
      CALL qs_upf_inverse_coulomb_evaluate(qs_env, state, state%coulomb, matrix, energy, operators, charge)
      IF (PRESENT(electronic_energy)) THEN
         ionic_energy = 0.5_dp*SUM(state%coulomb%poisson_weights*ABS(state%coulomb%ionic%array)**2)
         IF (state%coulomb%grid%para%mode /= PW_MODE_LOCAL) THEN
            CALL state%coulomb%grid%para%group%mp_comm_type%sum(ionic_energy)
         END IF
         electronic_energy = ionic_energy - energy
         DO k = 1, SIZE(matrix, 3)
            electronic_energy = electronic_energy + &
                                state%coulomb%weights(k)*REAL(SUM(matrix(:, :, k)*CONJG(operators(:, :, k))), dp)
         END DO
      END IF
      state%hartree_k(:, :, :) = operators
      state%coulomb_density_evaluation = state%density_evaluations
   END SUBROUTINE qs_upf_inverse_coulomb_update

! **************************************************************************************************
!> \brief Prepare all native sources, their finite metric and the density-independent Coulomb maps.
!>        Rebuild after native geometry, OCE, task-list, quadrature or defining-inverse changes.
!>        The ion source includes every kind and is separate from frozen partial-core XC density.
!>        Compensation shapes use the native ionic Gaussian exponent and the actual finite radial moments.
!> \param qs_env Live native environment with GAPW soft task lists and atom grids.
!> \param state Original defining inverse frames at the same geometry and unreduced mesh.
!> \param map Owned assembly. An earlier prepared assembly is released before rebuilding.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_coulomb_prepare(qs_env, state, map)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_upf_inverse_scf_type), INTENT(IN), TARGET  :: state
      TYPE(qs_upf_inverse_coulomb_type), INTENT(INOUT)   :: map

      COMPLEX(dp), ALLOCATABLE                           :: ao_operator(:, :)
      INTEGER                                            :: a, atom, b, channel, grid_shape(3), k, &
                                                            kind, nao, ncenter, nh, nr
      INTEGER(int_8)                                     :: remaining_bytes
      LOGICAL                                            :: paw_atom
      REAL(dp)                                           :: alpha, tolerance, total
      REAL(dp), ALLOCATABLE                              :: unit(:, :)
      TYPE(grid_atom_type), POINTER                      :: grid_atom
      TYPE(harmonics_atom_type), POINTER                 :: harmonics
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_set_p_type), POINTER            :: sab(:)
      TYPE(oce_matrix_type), POINTER                     :: oce
      TYPE(particle_type), POINTER                       :: atoms(:)
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_poisson_type), POINTER                     :: poisson
      TYPE(pw_pool_type), POINTER                        :: pool
      TYPE(qs_kind_type), POINTER                        :: kinds(:)
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_upf_inverse_augmentation_type), POINTER    :: augmentation
      TYPE(qs_upf_inverse_type), POINTER                 :: frame

      CPASSERT(state%enabled)
      IF (map%ready) CALL qs_upf_inverse_coulomb_release(map)
      CALL get_qs_env(qs_env, ks_env=ks_env, qs_kind_set=kinds, oce=oce, sab_orb=sab, &
                      para_env=para_env, pw_env=pw_env, particle_set=atoms)
      CPASSERT(ASSOCIATED(oce))
      CALL pw_env_get(pw_env, auxbas_pw_pool=pool, poisson_env=poisson)
      map%grid => pool%pw_grid
      CPASSERT(map%grid%vol > 0.0_dp)
      CPASSERT(map%grid%grid_span == FULLSPACE .OR. map%grid%grid_span == HALFSPACE)
      CALL map%source%create(map%grid)
      CALL map%potential%create(map%grid)
      CALL map%source_g%create(map%grid)
      CALL map%potential_g%create(map%grid)
      CALL map%ionic%create(map%grid)
      CALL pw_poisson_rebuild(poisson, map%source_g)
      CPASSERT(ASSOCIATED(poisson%green_fft))
      IF (poisson%parameters%has_dielectric .OR. &
          poisson%green_fft%p3m .OR. poisson%green_fft%sr_screening) THEN
         CPABORT("Inverse GAPW requires an unscreened Coulomb kernel without a dielectric")
      END IF
      SELECT CASE (poisson%green_fft%method)
      CASE (PERIODIC3D, ANALYTIC0D, ANALYTIC1D, ANALYTIC2D, MT0D, MT2D)
      CASE DEFAULT
         CPABORT("Inverse GAPW Coulomb for this Poisson solver is not implemented")
      END SELECT
      IF (poisson%green_fft%method == PERIODIC3D) THEN
         CALL ionic_fourier_source(atoms, kinds, map%ionic)
      ELSE
         CALL calculate_rho_core(map%ionic, total, qs_env, only_nopaw=.FALSE.)
      END IF
      ALLOCATE (map%multiplicity(SIZE(map%ionic%array)), map%poisson_weights(SIZE(map%ionic%array)))
      map%multiplicity = 1.0_dp
      IF (map%grid%grid_span == HALFSPACE) THEN
         map%multiplicity = 2.0_dp
         IF (map%grid%have_g0) map%multiplicity(1) = 1.0_dp
      END IF
      ! The inversion-symmetric MT kernel is real up to roundoff in its numerical FFT.
      tolerance = 64.0_dp*EPSILON(1.0_dp)*MAX(1.0_dp, MAXVAL(ABS(REAL(poisson%green_fft%influence_fn%array, dp))))
      CPASSERT(MAXVAL(ABS(AIMAG(poisson%green_fft%influence_fn%array))) <= tolerance)
      map%poisson_weights(:) = map%grid%vol*map%multiplicity*REAL(poisson%green_fft%influence_fn%array, dp)
      IF (ASSOCIATED(state%bloch%kpoints)) THEN
         ALLOCATE (map%weights, SOURCE=state%bloch%weights)
      ELSE
         ALLOCATE (map%weights(1))
         map%weights = 1.0_dp
      END IF
      CALL select_frame(state, 1, frame, augmentation)
      CPASSERT(ALL(frame%periodicity == poisson%parameters%periodic))
      nao = SIZE(frame%gaussian, 2)
      map%cell = frame%cell
      ! Wire and slab kernels use Cartesian isolated axes. Shear within a periodic plane is allowed.
      tolerance = 64.0_dp*EPSILON(1.0_dp)*MAXVAL(ABS(map%cell))
      DO a = 1, 3
         DO b = 1, 3
            IF (frame%periodicity(a) == frame%periodicity(b)) CYCLE
            IF (ABS(map%cell(a, b)) > tolerance) THEN
               CPABORT("Inverse GAPW requires Cartesian isolated axes")
            END IF
         END DO
      END DO
      ALLOCATE (map%positions, SOURCE=frame%positions)
      ncenter = SIZE(frame%center_atoms)
      DO atom = 1, SIZE(atoms)
         IF (ANY(frame%center_atoms == atom)) CYCLE
         kind = atoms(atom)%atomic_kind%kind_number
         CALL get_qs_kind(kinds(kind), paw_atom=paw_atom)
         IF (paw_atom) ncenter = ncenter + 1
      END DO
      ALLOCATE (map%centers(ncenter), map%pairs(ncenter), map%center_atoms(ncenter), &
                map%center_kinds(ncenter), map%inverse_centers(ncenter), map%fields(ncenter, SIZE(map%weights)))
      ! Bound replicated point fields across all centers and Bloch points in this map.
      remaining_bytes = 64_int_8*1024_int_8*1024_int_8
      ncenter = SIZE(frame%center_atoms)
      map%center_atoms(:ncenter) = frame%center_atoms
      map%center_kinds(:ncenter) = frame%center_kinds
      map%inverse_centers(:ncenter) = [(a, a=1, ncenter)]
      ! Every native softened Gaussian center needs its hard-minus-soft source, including ghosts.
      DO atom = 1, SIZE(atoms)
         IF (ANY(frame%center_atoms == atom)) CYCLE
         kind = atoms(atom)%atomic_kind%kind_number
         CALL get_qs_kind(kinds(kind), paw_atom=paw_atom)
         IF (.NOT. paw_atom) CYCLE
         ncenter = ncenter + 1
         map%center_atoms(ncenter) = atom
         map%center_kinds(ncenter) = kind
         map%inverse_centers(ncenter) = 0
      END DO
      CALL qs_upf_inverse_background_prepare(ks_env, map%background)
      ALLOCATE (map%metric(nao, nao, SIZE(map%weights)), ao_operator(nao, nao))
      map%potential%array = 1.0_dp
      CALL qs_upf_inverse_background_integrate(map%background, qs_env, map%potential, map%metric)
      DO a = 1, SIZE(map%centers)
         kind = map%center_kinds(a)
         CALL get_qs_kind(kinds(kind), grid_atom=grid_atom, harmonics=harmonics, alpha_core_charge=alpha)
         CPASSERT(alpha > 0.0_dp)
         CPASSERT(ASSOCIATED(grid_atom) .AND. ASSOCIATED(harmonics))
         nr = grid_atom%nr
         nh = harmonics%max_iso_not0
         ASSOCIATE (center => map%centers(a))
            center%position = frame%positions(:, map%center_atoms(a))
            ALLOCATE (center%radii, SOURCE=grid_atom%rad)
            ALLOCATE (center%weights, SOURCE=grid_atom%wr)
            ALLOCATE (center%l(nh), center%indices(nh), center%density(nr, nh), center%shapes(nr, nh))
            center%density = 0.0_dp
            DO channel = 1, nh
               center%l(channel) = indso(1, channel)
               center%indices(channel) = channel
               center%shapes(:, channel) = center%radii**center%l(channel)*EXP(-alpha*center%radii**2)
            END DO
         END ASSOCIATE
         CALL qs_upf_inverse_gaussian_pairs_prepare(oce, sab, kinds, map%center_atoms(a), kind, &
                                                 frame%ao_first, frame%ao_last, para_env, map%pairs(a), frame%positions, frame%cell)
         ALLOCATE (unit(nr, nh))
         unit = 0.0_dp
         unit(:, 1) = SQRT(fourpi)
         DO k = 1, SIZE(map%weights)
            CALL select_frame(state, k, frame, augmentation)
            CALL qs_upf_inverse_density_cache_prepare(frame, map%inverse_centers(a), kinds(kind), remaining_bytes, &
                                                      map%fields(a, k))
            CALL qs_upf_inverse_density_gapw_integrate(frame, augmentation, map%inverse_centers(a), kinds(kind), &
                                                       map%pairs(a), unit, ao_operator, map%center_atoms(a), map%fields(a, k))
            map%metric(:, :, k) = map%metric(:, :, k) + ao_operator
         END DO
         DEALLOCATE (unit)
      END DO
      grid_shape = 0
      IF (map%grid%grid_span == FULLSPACE .AND. .NOT. map%grid%spherical) grid_shape = map%grid%npts
      CALL upf_inverse_coulomb_prepare(map%centers, map%cell, map%grid%g, map%grid%vol, map%maps, &
                                       frame%periodicity, grid_shape, comm=para_env%mp_comm_type)
      map%ready = .TRUE.
   END SUBROUTINE qs_upf_inverse_coulomb_prepare

! **************************************************************************************************
!> \brief Evaluate the complete native periodic functional and its common Hartree AO operator.
!>        Includes Gaussian background, ions, every hard-minus-soft center, nR and original Q.
!>        Ionic self and overlap terms remain separate native energy terms. NLCC is excluded.
!> \param qs_env ...
!> \param state ...
!> \param map ...
!> \param matrix Hermitian AO occupations summed over spins, indexed by global k without mesh weights.
!> \param energy Coulomb energy including ionic interaction and the native ionic Gaussian self term.
!> \param ao_operators Same Hartree AO operator for every spin. Mesh weights belong to the outer trace.
!> \param charge Optional integrated electronic charge. Ionic and frozen-core charges are excluded.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_coulomb_evaluate(qs_env, state, map, matrix, energy, ao_operators, charge)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_upf_inverse_scf_type), INTENT(IN), TARGET  :: state
      TYPE(qs_upf_inverse_coulomb_type), INTENT(INOUT)   :: map
      COMPLEX(dp), INTENT(IN)                            :: matrix(:, :, :)
      REAL(dp), INTENT(OUT)                              :: energy
      COMPLEX(dp), INTENT(OUT)                           :: ao_operators(:, :, :)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: charge

      COMPLEX(dp), ALLOCATABLE                           :: ao_operator(:, :)
      INTEGER                                            :: a, k, kind
      TYPE(qs_kind_type), POINTER                        :: kinds(:)
      TYPE(qs_upf_inverse_augmentation_type), POINTER    :: augmentation
      TYPE(qs_upf_inverse_type), POINTER                 :: frame
      TYPE(upf_inverse_coulomb_potential_type), &
         ALLOCATABLE                                     :: potentials(:)

      CPASSERT(map%ready)
      CPASSERT(ALL(SHAPE(ao_operators) == SHAPE(map%metric)))
      CALL evaluate_coulomb_sources(qs_env, state, map, matrix, energy, potentials, charge)
      CALL get_qs_env(qs_env, qs_kind_set=kinds)
      ALLOCATE (ao_operator(SIZE(matrix, 1), SIZE(matrix, 2)))
      CALL qs_upf_inverse_background_integrate(map%background, qs_env, map%potential, ao_operators)
      DO k = 1, SIZE(map%weights)
         CALL select_frame(state, k, frame, augmentation)
         DO a = 1, SIZE(map%centers)
            kind = map%center_kinds(a)
            CALL qs_upf_inverse_density_gapw_integrate(frame, augmentation, map%inverse_centers(a), kinds(kind), &
                                                       map%pairs(a), potentials(a)%values, ao_operator, map%center_atoms(a), &
                                                       map%fields(a, k))
            ao_operators(:, :, k) = ao_operators(:, :, k) + ao_operator
         END DO
      END DO
   END SUBROUTINE qs_upf_inverse_coulomb_evaluate

! **************************************************************************************************
!> \brief Refresh the Coulomb sources and potentials for one spin-summed density.
!>        Geometry derivatives use these fields directly. The AO transpose is separate.
!> \param qs_env ...
!> \param state ...
!> \param map ...
!> \param matrix Hermitian AO occupation matrices, without k-point weights.
!> \param energy Coulomb energy with the same ionic convention as the AO evaluation.
!> \param potentials Radial Coulomb potentials for the center densities.
!> \param charge Optional integrated electronic charge.
! **************************************************************************************************
   SUBROUTINE evaluate_coulomb_sources(qs_env, state, map, matrix, energy, potentials, charge)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_upf_inverse_scf_type), INTENT(IN), TARGET  :: state
      TYPE(qs_upf_inverse_coulomb_type), INTENT(INOUT)   :: map
      COMPLEX(dp), INTENT(IN)                            :: matrix(:, :, :)
      REAL(dp), INTENT(OUT)                              :: energy
      TYPE(upf_inverse_coulomb_potential_type), &
         ALLOCATABLE, INTENT(OUT)                         :: potentials(:)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: charge

      COMPLEX(dp), ALLOCATABLE                           :: adjoint(:), background(:)
      INTEGER                                            :: a, k, kind
      REAL(dp)                                           :: electronic_charge, total
      REAL(dp), ALLOCATABLE                              :: radial(:, :)
      TYPE(qs_kind_type), POINTER                        :: kinds(:)
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_upf_inverse_augmentation_type), POINTER    :: augmentation
      TYPE(qs_upf_inverse_type), POINTER                 :: frame

      CPASSERT(map%ready)
      CPASSERT(ALL(SHAPE(matrix) == SHAPE(map%metric)))
      CALL select_frame(state, 1, frame, augmentation)
      CPASSERT(ALL(frame%cell == map%cell))
      CPASSERT(ALL(frame%positions == map%positions))
      CALL get_qs_env(qs_env, ks_env=ks_env, qs_kind_set=kinds)
      CALL qs_upf_inverse_background_collocate(map%background, ks_env, matrix, map%source, map%source_g, total)
      electronic_charge = -total
      DO a = 1, SIZE(map%centers)
         ALLOCATE (radial, MOLD=map%centers(a)%density)
         map%centers(a)%density = 0.0_dp
         DO k = 1, SIZE(map%weights)
            CALL select_frame(state, k, frame, augmentation)
            kind = map%center_kinds(a)
            CALL qs_upf_inverse_density_gapw(frame, augmentation, map%inverse_centers(a), kinds(kind), &
                                             map%pairs(a), matrix(:, :, k), radial, atom_index=map%center_atoms(a), &
                                             cache=map%fields(a, k))
            map%centers(a)%density(:, :) = map%centers(a)%density + map%weights(k)*radial
         END DO
         electronic_charge = electronic_charge + SQRT(fourpi)*SUM(map%centers(a)%weights*map%centers(a)%density(:, 1))
         DEALLOCATE (radial)
      END DO
      IF (PRESENT(charge)) charge = electronic_charge
      ALLOCATE (background, SOURCE=map%source_g%array + map%ionic%array)
      ALLOCATE (adjoint(SIZE(background)))
      IF (map%grid%para%mode == PW_MODE_LOCAL) THEN
         CALL upf_inverse_coulomb_periodic(map%centers, map%cell, map%grid%g, map%grid%vol, &
                                           background, map%poisson_weights, energy, potentials, adjoint, map=map%maps)
      ELSE
         CALL upf_inverse_coulomb_periodic(map%centers, map%cell, map%grid%g, map%grid%vol, &
                                           background, map%poisson_weights, energy, potentials, adjoint, &
                                           comm=map%grid%para%group%mp_comm_type, map=map%maps)
      END IF
      ! Convert the Fourier covector to the native unweighted potential before its grid transpose.
      map%potential_g%array = adjoint/(map%grid%vol*map%multiplicity)
      CALL pw_transfer(map%potential_g, map%potential)
   END SUBROUTINE evaluate_coulomb_sources

! **************************************************************************************************
!> \brief Apply the complete fixed-geometry GAPW Coulomb kernel to an electronic AO direction.
!>        Use the original smooth Gaussian source, signed center densities, compensation maps,
!>        complete center interactions and native Poisson weights. Fixed ionic fields have zero
!>        response. Independent matrix and PW workspaces preserve every accepted source/cache.
!>        The caller sums spin directions. Mesh weights enter source densities once, and the
!>        returned common Hartree operators are unweighted in the original global k-point order.
!> \param qs_env Native environment owning the original Gaussian task lists and atomic grids.
!> \param state Original inverse frames at the same fixed geometry.
!> \param map Immutable prepared Coulomb map, including its native boundary-condition kernel.
!> \param matrix_direction Replicated unweighted Hermitian spin-summed AO direction (AO,AO,k).
!> \param operators_direction Overwritten complete Hartree operator direction (AO,AO,k).
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_coulomb_operator_response(qs_env, state, map, matrix_direction, operators_direction)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_upf_inverse_scf_type), INTENT(IN), TARGET  :: state
      TYPE(qs_upf_inverse_coulomb_type), INTENT(IN)      :: map
      COMPLEX(dp), INTENT(IN)                            :: matrix_direction(:, :, :)
      COMPLEX(dp), INTENT(OUT)                           :: operators_direction(:, :, :)

      COMPLEX(dp), ALLOCATABLE                           :: adjoint(:), ao_operator(:, :)
      INTEGER                                            :: a, k, kind
      LOGICAL                                            :: valid
      REAL(dp)                                           :: quadratic_energy
      REAL(dp), ALLOCATABLE                              :: radial(:, :)
      TYPE(pw_c1d_gs_type)                               :: potential_g, source_g
      TYPE(pw_r3d_rs_type)                               :: potential, source
      TYPE(qs_kind_type), POINTER                        :: kinds(:)
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_upf_inverse_augmentation_type), POINTER    :: augmentation
      TYPE(qs_upf_inverse_background_type)               :: background
      TYPE(qs_upf_inverse_type), POINTER                 :: frame
      TYPE(upf_inverse_coulomb_center_type), ALLOCATABLE :: centers(:)
      TYPE(upf_inverse_coulomb_potential_type), &
         ALLOCATABLE                                     :: potentials(:)

      CPASSERT(map%ready)
      valid = ALL(SHAPE(matrix_direction) == SHAPE(map%metric))
      CPASSERT(valid)
      valid = ALL(SHAPE(operators_direction) == SHAPE(map%metric))
      CPASSERT(valid)
      CALL select_frame(state, 1, frame, augmentation)
      CPASSERT(ALL(frame%cell == map%cell))
      CPASSERT(ALL(frame%positions == map%positions))
      CALL get_qs_env(qs_env, ks_env=ks_env, qs_kind_set=kinds)
      CALL qs_upf_inverse_background_scratch(map%background, background)
      CALL source%create(map%grid)
      CALL source_g%create(map%grid)
      CALL potential%create(map%grid)
      CALL potential_g%create(map%grid)
      CALL qs_upf_inverse_background_collocate(background, ks_env, matrix_direction, source, source_g)
      ALLOCATE (centers, SOURCE=map%centers)
      DO a = 1, SIZE(centers)
         ALLOCATE (radial, MOLD=centers(a)%density)
         centers(a)%density = 0.0_dp
         kind = map%center_kinds(a)
         DO k = 1, SIZE(map%weights)
            CALL select_frame(state, k, frame, augmentation)
            CALL qs_upf_inverse_density_gapw(frame, augmentation, map%inverse_centers(a), kinds(kind), &
                                             map%pairs(a), matrix_direction(:, :, k), radial, &
                                             atom_index=map%center_atoms(a), cache=map%fields(a, k))
            centers(a)%density(:, :) = centers(a)%density + map%weights(k)*radial
         END DO
         DEALLOCATE (radial)
      END DO
      ALLOCATE (adjoint(SIZE(source_g%array)))
      ALLOCATE (ao_operator(SIZE(matrix_direction, 1), SIZE(matrix_direction, 2)))
      IF (map%grid%para%mode == PW_MODE_LOCAL) THEN
         CALL upf_inverse_coulomb_periodic(centers, map%cell, map%grid%g, map%grid%vol, &
                                           source_g%array, map%poisson_weights, quadratic_energy, potentials, &
                                           adjoint, map=map%maps)
      ELSE
         CALL upf_inverse_coulomb_periodic(centers, map%cell, map%grid%g, map%grid%vol, &
                                           source_g%array, map%poisson_weights, quadratic_energy, potentials, &
                                           adjoint, comm=map%grid%para%group%mp_comm_type, map=map%maps)
      END IF
      potential_g%array = adjoint/(map%grid%vol*map%multiplicity)
      CALL pw_transfer(potential_g, potential)
      CALL qs_upf_inverse_background_integrate(background, qs_env, potential, operators_direction)
      DO k = 1, SIZE(map%weights)
         CALL select_frame(state, k, frame, augmentation)
         DO a = 1, SIZE(centers)
            kind = map%center_kinds(a)
            CALL qs_upf_inverse_density_gapw_integrate(frame, augmentation, map%inverse_centers(a), kinds(kind), &
                                                       map%pairs(a), potentials(a)%values, ao_operator, &
                                                       atom_index=map%center_atoms(a), cache=map%fields(a, k))
            operators_direction(:, :, k) = operators_direction(:, :, k) + ao_operator
         END DO
      END DO
      CALL qs_upf_inverse_background_release(background)
      CALL source%release()
      CALL source_g%release()
      CALL potential%release()
      CALL potential_g%release()
   END SUBROUTINE qs_upf_inverse_coulomb_operator_response

! **************************************************************************************************
!> \brief Geometry covectors of the smooth Gaussian Coulomb sources at the prepared density.
!> \param qs_env Native Gaussian and ionic owners.
!> \param map Prepared source and potential.
!> \param matrix Fixed spin-summed AO occupations.
!> \param forces Negative Cartesian energy gradients.
!> \param virial Native virial, including both electronic and ionic sources.
! **************************************************************************************************
   SUBROUTINE coulomb_source_geometry(qs_env, map, matrix, forces, virial)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_upf_inverse_coulomb_type), INTENT(INOUT)   :: map
      COMPLEX(dp), INTENT(IN)                            :: matrix(:, :, :)
      REAL(dp), INTENT(OUT)                              :: forces(:, :), virial(3, 3)

      REAL(dp)                                           :: ionic_virial(3, 3)
      REAL(dp), ALLOCATABLE                              :: ionic_forces(:, :)
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_poisson_type), POINTER                     :: poisson

      ALLOCATE (ionic_forces, MOLD=forces)
      CALL qs_upf_inverse_background_geometry(map%background, qs_env, matrix, map%potential, forces, virial)
      CALL get_qs_env(qs_env, pw_env=pw_env)
      CALL pw_env_get(pw_env, poisson_env=poisson)
      IF (poisson%green_fft%method == PERIODIC3D) THEN
         CALL ionic_fourier_geometry(qs_env, map, ionic_forces, ionic_virial)
      ELSE
         CALL qs_upf_inverse_coulomb_ionic_geometry(qs_env, map%potential, ionic_forces, ionic_virial)
      END IF
      forces(:, :) = forces + ionic_forces
      virial(:, :) = virial + ionic_virial
   END SUBROUTINE coulomb_source_geometry

! **************************************************************************************************
!> \brief Joint atom/cell derivative of the complete native Coulomb functional at fixed AO occupations.
!>        Includes the actual Gaussian electron/ion grid sources, complete center densities, coupled
!>        inverse coefficients, compensation, periodic pair kernels and native Coulomb measure.
!>        Reduced Bloch coordinates and native Cartesian atom-grid offsets remain fixed.
!> \param qs_env Native owners at the same geometry as map.
!> \param state Original defining inverse frames, unchanged by this directional derivative.
!> \param map Complete source map. Its primal source/potential scratch is updated, without changing P.
!> \param matrix Spin-summed Hermitian occupations per global k, without mesh weights.
!> \param atom_response Physical Cartesian motion of all original atoms.
!> \param cell_response Physical lattice derivative, with lattice vectors in columns.
!> \param response Coulomb energy derivative, excluding separate native ionic self/overlap corrections.
!> \param charge_response Optional derivative of the same complete electronic charge metric contraction.
!> \param coefficients_response Optional defining-frame responses at fixed reduced k for this geometry variation.
!> \param pair_geometry_response Optional residual-kernel contraction captured at this accepted density.
!> \param source_forces Optional smooth-source forces prepared with the pair geometry at this density.
!> \param source_virial Corresponding smooth-source virial. Both arguments must be supplied together.
!> \param metric_matrix Optional energy-weighted densities in the same AO and k ordering.
!> \param metric_response One-center charge-metric derivative at metric_matrix, before the Pulay sign.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_coulomb_response(qs_env, state, map, matrix, atom_response, cell_response, &
                                              response, charge_response, coefficients_response, pair_geometry_response, &
                                              source_forces, source_virial, metric_matrix, metric_response)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_upf_inverse_scf_type), INTENT(IN), TARGET  :: state
      TYPE(qs_upf_inverse_coulomb_type), INTENT(INOUT)   :: map
      COMPLEX(dp), INTENT(IN)                            :: matrix(:, :, :)
      REAL(dp), INTENT(IN)                               :: atom_response(:, :), cell_response(3, 3)
      REAL(dp), INTENT(OUT)                              :: response
      REAL(dp), INTENT(OUT), OPTIONAL                    :: charge_response
      COMPLEX(dp), INTENT(IN), OPTIONAL, TARGET          :: coefficients_response(:, :, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: pair_geometry_response, source_forces(:, :), &
                                                            source_virial(3, 3)
      COMPLEX(dp), INTENT(IN), OPTIONAL, TARGET          :: metric_matrix(:, :, :)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: metric_response

      COMPLEX(dp), ALLOCATABLE                           :: adjoint(:), background(:), &
                                                            background_response(:)
      INTEGER                                            :: a, info
      REAL(dp)                                           :: energy, inverse_cell(3, 3), &
                                                            kernel_response, &
                                                            strain(3, 3), trace_strain, &
                                                            virial(3, 3)
      REAL(dp), ALLOCATABLE                              :: forces(:, :), gvecs_response(:, :), &
                                                            poisson_weights_response(:)
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_poisson_type), POINTER                     :: poisson
      TYPE(pw_r3d_rs_type)                               :: unit
      TYPE(upf_inverse_coulomb_center_type), ALLOCATABLE :: centers_response(:)
      TYPE(upf_inverse_coulomb_potential_type), &
         ALLOCATABLE                                     :: potentials(:)

      CPASSERT(map%ready)
      CPASSERT(ALL(SHAPE(atom_response) == SHAPE(map%positions)))
      CALL prepare_pair_responses(qs_env, state, map)
      CPASSERT(PRESENT(source_forces) .EQV. PRESENT(source_virial))
      inverse_cell = map%cell
      CALL invmat(inverse_cell, info)
      CPASSERT(info == 0)
      strain = MATMUL(cell_response, inverse_cell)
      trace_strain = strain(1, 1) + strain(2, 2) + strain(3, 3)
      ! Fourier coefficients use a grid average. Its source derivative has no volume-measure term.
      IF (PRESENT(source_forces)) THEN
         response = -SUM(source_forces*(atom_response - MATMUL(strain, map%positions))) - SUM(source_virial*strain)
      ELSE
         CALL evaluate_coulomb_sources(qs_env, state, map, matrix, energy, potentials)
         ALLOCATE (forces, MOLD=atom_response)
         CALL coulomb_source_geometry(qs_env, map, matrix, forces, virial)
         response = -SUM(forces*(atom_response - MATMUL(strain, map%positions))) - SUM(virial*strain)
      END IF
      IF (PRESENT(charge_response)) THEN
         IF (.NOT. ALLOCATED(forces)) ALLOCATE (forces, MOLD=atom_response)
         CALL unit%create(map%grid)
         unit%array = 1.0_dp
         CALL qs_upf_inverse_background_geometry(map%background, qs_env, matrix, unit, forces, virial)
         charge_response = -SUM(forces*(atom_response - MATMUL(strain, map%positions))) - &
                           SUM(virial*strain) + trace_strain*pw_integral_ab(map%source, unit)
         CALL unit%release()
      END IF
      CALL center_density_response(qs_env, state, map, matrix, atom_response, cell_response, centers_response, &
                                   coefficients_response, metric_matrix, metric_response)
      IF (PRESENT(charge_response)) THEN
         DO a = 1, SIZE(centers_response)
            charge_response = charge_response + SQRT(fourpi)* &
                              SUM(map%centers(a)%weights*centers_response(a)%density(:, 1))
         END DO
      END IF
      ALLOCATE (gvecs_response(3, SIZE(map%grid%g, 2)), poisson_weights_response(SIZE(map%poisson_weights)))
      gvecs_response(:, :) = -MATMUL(TRANSPOSE(strain), map%grid%g)
      CALL get_qs_env(qs_env, pw_env=pw_env)
      CALL pw_env_get(pw_env, poisson_env=poisson)
      CALL pw_green_strain_response(poisson%green_fft, strain, poisson_weights_response)
      poisson_weights_response(:) = trace_strain*map%poisson_weights + &
                                    map%grid%vol*map%multiplicity*poisson_weights_response
      ALLOCATE (background, SOURCE=map%source_g%array + map%ionic%array)
      ALLOCATE (adjoint(SIZE(background)), background_response(SIZE(background)))
      background_response = CMPLX(0.0_dp, 0.0_dp, dp)
      IF (map%grid%para%mode == PW_MODE_LOCAL) THEN
         CALL upf_inverse_coulomb_periodic(map%centers, map%cell, map%grid%g, map%grid%vol, &
                                           background, map%poisson_weights, energy, background_adjoint=adjoint, &
                                           response=kernel_response, center_response=centers_response, &
                                           cell_response=cell_response, &
                                           gvecs_response=gvecs_response, volume_response=trace_strain*map%grid%vol, &
                                           background_response=background_response, &
                                           poisson_weights_response=poisson_weights_response, map=map%maps, &
                                           pair_geometry_response=pair_geometry_response)
      ELSE
         CALL upf_inverse_coulomb_periodic(map%centers, map%cell, map%grid%g, map%grid%vol, &
                                           background, map%poisson_weights, energy, background_adjoint=adjoint, &
                                           comm=map%grid%para%group%mp_comm_type, response=kernel_response, &
                                           center_response=centers_response, cell_response=cell_response, &
                                           gvecs_response=gvecs_response, volume_response=trace_strain*map%grid%vol, &
                                           background_response=background_response, &
                                           poisson_weights_response=poisson_weights_response, map=map%maps, &
                                           pair_geometry_response=pair_geometry_response)
      END IF
      response = response + kernel_response
   END SUBROUTINE qs_upf_inverse_coulomb_response

! **************************************************************************************************
!> \brief Capture residual-kernel covectors for all motions of the same accepted density.
!>        Native center responses keep their radial shells fixed. Other density and source
!>        responses continue to use the complete original directional derivative.
!>        Replicated center fields share pair/image derivatives across the native process group.
!> \param qs_env Native source owners.
!> \param state Original defining inverse frames.
!> \param map Complete Coulomb source, refreshed from matrix before contraction.
!> \param matrix Accepted spin-summed Hermitian matrices in global k order.
!> \param atom_gradient Covector for physical Cartesian atom positions.
!> \param cell_gradient Covector for lattice vectors at fixed atom positions.
!> \param source_forces Optional negative Cartesian gradients of the smooth sources at these occupations.
!> \param source_virial Corresponding smooth-source virial.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_coulomb_pair_geometry(qs_env, state, map, matrix, atom_gradient, cell_gradient, &
                                                   source_forces, source_virial)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_upf_inverse_scf_type), INTENT(IN), TARGET  :: state
      TYPE(qs_upf_inverse_coulomb_type), INTENT(INOUT)   :: map
      COMPLEX(dp), INTENT(IN)                            :: matrix(:, :, :)
      REAL(dp), INTENT(OUT)                              :: atom_gradient(:, :), cell_gradient(3, 3)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: source_forces(:, :), source_virial(3, 3)

      INTEGER                                            :: a
      REAL(dp)                                           :: energy
      REAL(dp), ALLOCATABLE                              :: center_gradient(:, :)
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(upf_inverse_coulomb_potential_type), &
         ALLOCATABLE                                     :: potentials(:)

      CPASSERT(map%ready)
      CPASSERT(ALL(SHAPE(atom_gradient) == SHAPE(map%positions)))
      CPASSERT(PRESENT(source_forces) .EQV. PRESENT(source_virial))
      CALL evaluate_coulomb_sources(qs_env, state, map, matrix, energy, potentials)
      IF (PRESENT(source_forces)) CALL coulomb_source_geometry(qs_env, map, matrix, source_forces, source_virial)
      CALL get_qs_env(qs_env, para_env=para_env)
      ALLOCATE (center_gradient(3, SIZE(map%centers)))
      CALL upf_inverse_coulomb_pair_geometry(map%centers, map%maps, center_gradient, cell_gradient, &
                                             comm=para_env%mp_comm_type)
      atom_gradient = 0.0_dp
      DO a = 1, SIZE(map%centers)
         atom_gradient(:, map%center_atoms(a)) = atom_gradient(:, map%center_atoms(a)) + center_gradient(:, a)
      END DO
   END SUBROUTINE qs_upf_inverse_coulomb_pair_geometry

! **************************************************************************************************
!> \brief Geometry derivative of the complete electronic charge metric at fixed Hermitian matrix.
!>        The matrix may be an energy-weighted density for the stationary Pulay contraction.
!>        No ionic/core charges or Coulomb energy enter this linear functional. Native density,
!>        installed metric and accepted Coulomb potential owners are unchanged.
!> \param qs_env Current native Gaussian source owners.
!> \param state Defining inverse frames in the installed metric's global AO and k ordering.
!> \param map Complete source layout at the same geometry.
!> \param matrix Hermitian per-k matrices, without mesh weights, summed over spins if appropriate.
!> \param atom_response Physical atomic motion.
!> \param cell_response Physical lattice derivative at fixed reduced Bloch coordinates.
!> \param response Derivative of sum_k w_k Tr(matrix_k M_k), without the Pulay minus sign.
!> \param coefficients_response Optional defining-frame responses at fixed reduced k for this geometry variation.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_coulomb_metric_response(qs_env, state, map, matrix, atom_response, cell_response, response, &
                                                     coefficients_response)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_upf_inverse_scf_type), INTENT(IN), TARGET  :: state
      TYPE(qs_upf_inverse_coulomb_type), INTENT(INOUT)   :: map
      COMPLEX(dp), INTENT(IN)                            :: matrix(:, :, :)
      REAL(dp), INTENT(IN)                               :: atom_response(:, :), cell_response(3, 3)
      REAL(dp), INTENT(OUT)                              :: response
      COMPLEX(dp), INTENT(IN), OPTIONAL, TARGET          :: coefficients_response(:, :, :)

      REAL(dp)                                           :: atom_gradient(3, SIZE(atom_response, 2)), &
                                                            cell_gradient(3, 3)

      CPASSERT(ALL(SHAPE(atom_response) == SHAPE(map%positions)))
      CALL qs_upf_inverse_coulomb_metric_geometry(qs_env, map, matrix, atom_gradient, cell_gradient)
      CALL qs_upf_inverse_coulomb_center_metric_response(qs_env, state, map, matrix, atom_response, cell_response, &
                                                         response, coefficients_response)
      response = response + SUM(atom_gradient*atom_response) + SUM(cell_gradient*cell_response)
   END SUBROUTINE qs_upf_inverse_coulomb_metric_response

! **************************************************************************************************
!> \brief Smooth charge-metric gradients shared by all variations at a fixed AO matrix.
!> \param qs_env Current native Gaussian source owners.
!> \param map Complete source layout at the same geometry.
!> \param matrix Hermitian per-k matrices, without mesh weights.
!> \param atom_gradient Derivative with respect to Cartesian atom positions, without the Pulay sign.
!> \param cell_gradient Derivative with respect to lattice vectors at fixed atom positions.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_coulomb_metric_geometry(qs_env, map, matrix, atom_gradient, cell_gradient)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_upf_inverse_coulomb_type), INTENT(INOUT)   :: map
      COMPLEX(dp), INTENT(IN)                            :: matrix(:, :, :)
      REAL(dp), INTENT(OUT)                              :: atom_gradient(:, :), cell_gradient(3, 3)

      INTEGER                                            :: axis, info
      REAL(dp)                                           :: charge, inverse_cell(3, 3), virial(3, 3)
      TYPE(pw_c1d_gs_type)                               :: background_g
      TYPE(pw_r3d_rs_type)                               :: background, unit
      TYPE(qs_ks_env_type), POINTER                      :: ks_env

      CPASSERT(map%ready)
      CPASSERT(ALL(SHAPE(matrix) == SHAPE(map%metric)))
      CPASSERT(ALL(SHAPE(atom_gradient) == SHAPE(map%positions)))
      CALL get_qs_env(qs_env, ks_env=ks_env)
      inverse_cell = map%cell
      CALL invmat(inverse_cell, info)
      CPASSERT(info == 0)
      CALL background%create(map%grid)
      CALL background_g%create(map%grid)
      CALL unit%create(map%grid)
      unit%array = 1.0_dp
      CALL qs_upf_inverse_background_collocate(map%background, ks_env, matrix, background, background_g)
      charge = pw_integral_ab(background, unit)
      CALL qs_upf_inverse_background_geometry(map%background, qs_env, matrix, unit, atom_gradient, virial)
      ! The native virial uses affine motion. Convert to lattice derivatives at fixed positions.
      virial = MATMUL(atom_gradient, TRANSPOSE(map%positions)) - virial
      DO axis = 1, 3
         virial(axis, axis) = virial(axis, axis) + charge
      END DO
      cell_gradient = MATMUL(virial, TRANSPOSE(inverse_cell))
      atom_gradient = -atom_gradient
      CALL unit%release()
      CALL background%release()
      CALL background_g%release()
   END SUBROUTINE qs_upf_inverse_coulomb_metric_geometry

! **************************************************************************************************
!> \brief Response of the one-center part of the charge metric, without the Pulay sign.
!> \param qs_env Native owners of the original radial fields.
!> \param state Defining inverse frames.
!> \param map Complete source layout at the same geometry.
!> \param matrix Hermitian per-k matrices, without mesh weights.
!> \param atom_response Physical atomic motion.
!> \param cell_response Physical lattice derivative at fixed reduced Bloch coordinates.
!> \param response Derivative of the integrated one-center electronic density.
!> \param coefficients_response Optional defining-frame responses for this variation.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_coulomb_center_metric_response(qs_env, state, map, matrix, atom_response, cell_response, &
                                                            response, coefficients_response)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_upf_inverse_scf_type), INTENT(IN), TARGET  :: state
      TYPE(qs_upf_inverse_coulomb_type), INTENT(INOUT)   :: map
      COMPLEX(dp), INTENT(IN)                            :: matrix(:, :, :)
      REAL(dp), INTENT(IN)                               :: atom_response(:, :), cell_response(3, 3)
      REAL(dp), INTENT(OUT)                              :: response
      COMPLEX(dp), INTENT(IN), OPTIONAL, TARGET          :: coefficients_response(:, :, :)

      INTEGER                                            :: a
      TYPE(upf_inverse_coulomb_center_type), ALLOCATABLE :: centers_response(:)

      CPASSERT(map%ready)
      CPASSERT(ALL(SHAPE(matrix) == SHAPE(map%metric)))
      CPASSERT(ALL(SHAPE(atom_response) == SHAPE(map%positions)))
      CALL prepare_pair_responses(qs_env, state, map)
      CALL center_density_response(qs_env, state, map, matrix, atom_response, cell_response, centers_response, &
                                   coefficients_response)
      response = 0.0_dp
      DO a = 1, SIZE(centers_response)
         response = response + SQRT(fourpi)*SUM(map%centers(a)%weights*centers_response(a)%density(:, 1))
      END DO
   END SUBROUTINE qs_upf_inverse_coulomb_center_metric_response

! **************************************************************************************************
!> \brief Complete center density derivative at fixed native radial nodes and measures.
!> \param qs_env Native kinds with their original radial and Gaussian one-center fields.
!> \param state Original defining inverse frames.
!> \param map Current complete source layout and derivative pair coefficients.
!> \param matrix Hermitian spin-summed matrices per global k, without mesh weights.
!> \param atom_response Physical atomic motion.
!> \param cell_response Physical lattice derivative.
!> \param centers_response New center layout with density and position responses.
!> \param known_coefficients_response Optional defining-frame responses for the same variation and global k order.
!> \param metric_matrix Optional second set of per-k matrices for the integrated one-center response.
!> \param metric_response Optional integrated one-center response evaluated with metric_matrix.
! **************************************************************************************************
   SUBROUTINE center_density_response(qs_env, state, map, matrix, atom_response, cell_response, centers_response, &
                                      known_coefficients_response, metric_matrix, metric_response)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_upf_inverse_scf_type), INTENT(IN), TARGET  :: state
      TYPE(qs_upf_inverse_coulomb_type), INTENT(IN)      :: map
      COMPLEX(dp), INTENT(IN)                            :: matrix(:, :, :)
      REAL(dp), INTENT(IN)                               :: atom_response(:, :), cell_response(3, 3)
      TYPE(upf_inverse_coulomb_center_type), &
         ALLOCATABLE, INTENT(OUT)                        :: centers_response(:)
      COMPLEX(dp), INTENT(IN), OPTIONAL, TARGET :: known_coefficients_response(:, :, :)
      COMPLEX(dp), INTENT(IN), OPTIONAL, TARGET          :: metric_matrix(:, :, :)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: metric_response

      COMPLEX(dp), ALLOCATABLE, TARGET                   :: owned_coefficients(:, :)
      COMPLEX(dp), POINTER                               :: coefficients_response(:, :), metric_k(:, :)
      INTEGER                                            :: a, info, k, kind, nao, nk, nproj
      REAL(dp)                                           :: inverse_cell(3, 3), strain(3, 3), &
                                                            trace_strain
      REAL(dp), ALLOCATABLE                              :: metric_radial(:, :), point_response(:, :), &
                                                            radial_response(:, :), &
                                                            weights_response(:)
      TYPE(qs_kind_type), POINTER                        :: kinds(:)
      TYPE(qs_upf_inverse_augmentation_type), POINTER    :: augmentation
      TYPE(qs_upf_inverse_type), POINTER                 :: frame

      CPASSERT(PRESENT(metric_matrix) .EQV. PRESENT(metric_response))
      NULLIFY (metric_k)
      IF (PRESENT(metric_response)) THEN
         CPASSERT(ALL(SHAPE(metric_matrix) == SHAPE(matrix)))
         metric_response = 0.0_dp
      END IF
      CALL get_qs_env(qs_env, qs_kind_set=kinds)
      inverse_cell = map%cell
      CALL invmat(inverse_cell, info)
      CPASSERT(info == 0)
      strain = MATMUL(cell_response, inverse_cell)
      trace_strain = strain(1, 1) + strain(2, 2) + strain(3, 3)
      ALLOCATE (centers_response, SOURCE=map%centers)
      CALL select_frame(state, 1, frame, augmentation)
      DO a = 1, SIZE(centers_response)
         centers_response(a)%position = atom_response(:, map%center_atoms(a))
         centers_response(a)%radii = 0.0_dp
         centers_response(a)%weights = 0.0_dp
         centers_response(a)%shapes = 0.0_dp
         centers_response(a)%density = 0.0_dp
      END DO
      nao = SIZE(frame%coefficients, 2)
      nproj = SIZE(frame%coefficients, 1)
      nk = SIZE(map%weights)
      IF (PRESENT(known_coefficients_response)) THEN
         CPASSERT(ALL(SHAPE(known_coefficients_response) == [nproj, nao, nk]))
      ELSE
         ALLOCATE (owned_coefficients, MOLD=frame%coefficients)
      END IF
      DO k = 1, SIZE(map%weights)
         CALL select_frame(state, k, frame, augmentation)
         IF (PRESENT(metric_matrix)) metric_k => metric_matrix(:, :, k)
         IF (PRESENT(known_coefficients_response)) THEN
            coefficients_response => known_coefficients_response(:, :, k)
         ELSE
            coefficients_response => owned_coefficients
            ALLOCATE (point_response(3, SIZE(frame%weights)), weights_response(SIZE(frame%weights)))
            point_response(:, :) = MATMUL(strain, frame%points)
            weights_response(:) = trace_strain*frame%weights
            CALL qs_upf_inverse_response(frame, point_response, atom_response, cell_response, &
                                         [0.0_dp, 0.0_dp, 0.0_dp], weights_response, coefficients_response)
            DEALLOCATE (point_response, weights_response)
         END IF
         DO a = 1, SIZE(centers_response)
            kind = map%center_kinds(a)
            ALLOCATE (radial_response, MOLD=map%centers(a)%density)
            IF (PRESENT(metric_response)) ALLOCATE (metric_radial(SIZE(radial_response, 1), 1))
            CALL qs_upf_inverse_density_gapw(frame, augmentation, map%inverse_centers(a), kinds(kind), &
                                             map%pairs(a), matrix(:, :, k), response=radial_response, &
                                             atom_response=atom_response, cell_response=cell_response, &
                                             reduced_k_response=[0.0_dp, 0.0_dp, 0.0_dp], &
                                             coefficients_response=coefficients_response, atom_index=map%center_atoms(a), &
                                             metric_matrix=metric_k, metric_response=metric_radial)
            centers_response(a)%density(:, :) = centers_response(a)%density + map%weights(k)*radial_response
            IF (PRESENT(metric_response)) THEN
               metric_response = metric_response + map%weights(k)*SUM(map%centers(a)%weights*metric_radial(:, 1))
               DEALLOCATE (metric_radial)
            END IF
            DEALLOCATE (radial_response)
         END DO
      END DO
      IF (PRESENT(metric_response)) metric_response = SQRT(fourpi)*metric_response
   END SUBROUTINE center_density_response

! **************************************************************************************************
!> \brief Prepare missing native OCE derivatives in private pair caches at the same geometry.
!>        Energy-only native OCE owners need not contain derivatives. Their primal tables and
!>        ownership are unchanged while the same native builder supplies the response coefficients.
!> \param qs_env Native kinds, atoms, overlap lists and OCE neighbor layout.
!> \param state Original inverse AO and center ordering.
!> \param map Complete source map at the same current geometry.
! **************************************************************************************************
   SUBROUTINE prepare_pair_responses(qs_env, state, map)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_upf_inverse_scf_type), INTENT(IN), TARGET  :: state
      TYPE(qs_upf_inverse_coulomb_type), INTENT(INOUT)   :: map

      INTEGER                                            :: a, derivative_order, kind
      LOGICAL                                            :: needed
      TYPE(dft_control_type), POINTER                    :: control
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_set_p_type), POINTER            :: sab(:), sap(:)
      TYPE(oce_matrix_type), POINTER                     :: oce
      TYPE(particle_type), POINTER                       :: atoms(:)
      TYPE(qs_kind_type), POINTER                        :: kinds(:)
      TYPE(qs_upf_inverse_augmentation_type), POINTER    :: augmentation
      TYPE(qs_upf_inverse_type), POINTER                 :: frame

      needed = .FALSE.
      DO a = 1, SIZE(map%pairs)
         needed = needed .OR. .NOT. qs_upf_inverse_gaussian_pairs_has_derivatives(map%pairs(a))
      END DO
      IF (.NOT. needed) RETURN
      CALL get_qs_env(qs_env, dft_control=control, para_env=para_env, qs_kind_set=kinds, &
                      particle_set=atoms, sab_orb=sab, sap_oce=sap)
      CPASSERT(ASSOCIATED(sap))
      NULLIFY (oce)
      CALL create_oce_set(oce)
      CALL allocate_oce_set(oce, SIZE(kinds))
      derivative_order = 1
      CALL build_oce_matrices(oce%intac, .TRUE., derivative_order, kinds, atoms, sap, &
                              control%qs_control%gapw_control%eps_fit)
      CALL select_frame(state, 1, frame, augmentation)
      DO a = 1, SIZE(map%pairs)
         IF (qs_upf_inverse_gaussian_pairs_has_derivatives(map%pairs(a))) CYCLE
         kind = map%center_kinds(a)
         CALL qs_upf_inverse_gaussian_pairs_prepare(oce, sab, kinds, map%center_atoms(a), kind, &
                                                 frame%ao_first, frame%ao_last, para_env, map%pairs(a), frame%positions, frame%cell)
      END DO
      CALL deallocate_oce_set(oce)
   END SUBROUTINE prepare_pair_responses

! **************************************************************************************************
!> \brief Return the metric of exactly the electronic sources used in the Coulomb functional.
!> \param map ...
!> \param metric Overwritten per-k matrices. A unit source potential defines their normalization.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_coulomb_metric(map, metric)
      TYPE(qs_upf_inverse_coulomb_type), INTENT(IN)      :: map
      COMPLEX(dp), INTENT(OUT)                           :: metric(:, :, :)

      CPASSERT(map%ready)
      CPASSERT(ALL(SHAPE(metric) == SHAPE(map%metric)))
      metric(:, :, :) = map%metric
   END SUBROUTINE qs_upf_inverse_coulomb_metric

! **************************************************************************************************
!> \brief Install this complete electronic charge metric before the first native SCF density.
!>        The caller must use this same complete density and its operators throughout SCF.
!>        Native energy initialization must attach the source before the first SCF density.
!> \param qs_env Fresh inverse environment at the geometry used to prepare the Coulomb map.
!> \param map Complete native density functional, including all center monopoles.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_coulomb_install_metric(qs_env, map)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_upf_inverse_coulomb_type), INTENT(IN)      :: map

      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_upf_inverse_augmentation_type), POINTER    :: augmentation
      TYPE(qs_upf_inverse_scf_type), POINTER             :: state
      TYPE(qs_upf_inverse_type), POINTER                 :: frame

      CPASSERT(map%ready)
      CALL get_qs_env(qs_env, ks_env=ks_env)
      CALL get_ks_env(ks_env, upf_inverse=state)
      CPASSERT(ASSOCIATED(state))
      CALL select_frame(state, 1, frame, augmentation)
      CPASSERT(ALL(frame%cell == map%cell))
      CPASSERT(ALL(frame%positions == map%positions))
      CALL qs_upf_inverse_scf_install_metric(ks_env, map%metric)
   END SUBROUTINE qs_upf_inverse_coulomb_install_metric

! **************************************************************************************************
!> \brief Select a borrowed original inverse frame without modifying accepted density state.
!> \param state ...
!> \param k ...
!> \param frame ...
!> \param augmentation ...
! **************************************************************************************************
   SUBROUTINE select_frame(state, k, frame, augmentation)
      TYPE(qs_upf_inverse_scf_type), INTENT(IN), TARGET  :: state
      INTEGER, INTENT(IN)                                :: k
      TYPE(qs_upf_inverse_type), POINTER                 :: frame
      TYPE(qs_upf_inverse_augmentation_type), POINTER    :: augmentation

      IF (ASSOCIATED(state%bloch%kpoints)) THEN
         CPASSERT(k >= 1 .AND. k <= SIZE(state%bloch%frames))
         frame => state%bloch%frames(k)%ao
         augmentation => state%bloch%augmentation
      ELSE
         CPASSERT(k == 1)
         frame => state%frame%ao
         augmentation => state%augmentation
      END IF
   END SUBROUTINE select_frame

END MODULE qs_upf_inverse_coulomb
