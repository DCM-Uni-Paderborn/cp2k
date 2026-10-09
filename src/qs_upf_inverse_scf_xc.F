!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
!> \brief Native SCF consumers of the inverse atom-centered XC residual.
MODULE qs_upf_inverse_scf_xc
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE cell_types,                      ONLY: cell_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type
   USE input_section_types,             ONLY: section_get_ival,&
                                              section_vals_type
   USE kinds,                           ONLY: dp
   USE lebedev,                         ONLY: deallocate_lebedev_grids,&
                                              get_number_of_lebedev_grid,&
                                              init_lebedev_grids,&
                                              lebedev_grid
   USE mathconstants,                   ONLY: fourpi
   USE message_passing,                 ONLY: mp_comm_type,&
                                              mp_para_env_type
   USE particle_types,                  ONLY: particle_type
   USE paw_core_kinetic,                ONLY: paw_core_kinetic_ae_cartesian,&
                                              paw_core_kinetic_cartesian
   USE paw_xc,                          ONLY: paw_xc_ae_core_cartesian
   USE pw_types,                        ONLY: pw_r3d_rs_type
   USE qs_grid_atom,                    ONLY: allocate_grid_atom,&
                                              create_grid_atom,&
                                              grid_atom_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_upf_inverse,                  ONLY: qs_upf_inverse_response
   USE qs_upf_inverse_atom_grid,        ONLY: qs_upf_inverse_atom_grid_prepare,&
                                              qs_upf_inverse_atom_grid_ref_type,&
                                              qs_upf_inverse_atom_grid_response
   USE qs_upf_inverse_grid,             ONLY: qs_upf_inverse_matrix_gather
   USE qs_upf_inverse_kpoints,          ONLY: qs_upf_inverse_kpoints_density
   USE qs_upf_inverse_types,            ONLY: qs_upf_inverse_augmentation_type,&
                                              qs_upf_inverse_scf_type,&
                                              qs_upf_inverse_type
   USE qs_upf_inverse_xc,               ONLY: qs_upf_inverse_frame_ref_type,&
                                              qs_upf_inverse_xc_grid_evaluate,&
                                              qs_upf_inverse_xc_grid_response,&
                                              qs_upf_inverse_xc_pw_adjoint,&
                                              qs_upf_inverse_xc_pw_fields
   USE skala_gpw_functional,            ONLY: xc_section_uses_native_skala_evaluator
   USE upf_radial_fields,               ONLY: upf_radial_fields_init,&
                                              upf_radial_fields_periodic
   USE upf_types,                       ONLY: upf_nlcc_kinetic
   USE xc_input_constants,              ONLY: xc_deriv_pw,&
                                              xc_rho_no_smooth
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: qs_upf_inverse_scf_xc_evaluate, qs_upf_inverse_scf_xc_response

CONTAINS

! **************************************************************************************************
!> \brief Evaluate the same atom residual for native energy, direct AO and background operators.
!> \param ks_env Current native inverse environment, supplying the accepted occupation matrices.
!> \param density Actual native XC spin densities, including the native frozen-core background.
!> \param xc_section Functional and density cutoffs shared with the native grid evaluation.
!> \param energy Globally reduced atom residual, to be stored as the native one-center XC energy.
!> \param potential Optional native unweighted XC potentials to receive the background transpose.
!> \param kinetic Optional native kinetic density, including the independently supplied core.
!> \param kinetic_potential Optional unweighted kinetic XC potentials receiving the same transpose.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_scf_xc_evaluate(ks_env, density, xc_section, energy, potential, kinetic, kinetic_potential)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: density(:)
      TYPE(section_vals_type), POINTER                   :: xc_section
      REAL(dp), INTENT(OUT)                              :: energy
      TYPE(pw_r3d_rs_type), INTENT(INOUT), OPTIONAL      :: potential(:)
      TYPE(pw_r3d_rs_type), INTENT(IN), OPTIONAL         :: kinetic(:)
      TYPE(pw_r3d_rs_type), INTENT(INOUT), OPTIONAL      :: kinetic_potential(:)

      COMPLEX(dp), ALLOCATABLE                           :: gamma(:, :), matrix(:, :, :)
      INTEGER                                            :: nao, nk, np, ns, spin
      REAL(dp), ALLOCATABLE                              :: adjoint(:, :, :, :), core(:, :), &
                                                            core_gradient(:, :, :), core_tau(:, :)
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_p_type), POINTER                        :: density_ao(:, :)
      TYPE(mp_comm_type)                                 :: comm
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(pw_r3d_rs_type), ALLOCATABLE                  :: background_kinetic(:), &
                                                            background_potential(:)
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_upf_inverse_augmentation_type), POINTER    :: augmentation
      TYPE(qs_upf_inverse_frame_ref_type), ALLOCATABLE   :: frames(:)
      TYPE(qs_upf_inverse_scf_type), POINTER             :: state

      CALL get_ks_env(ks_env, upf_inverse=state, rho=rho, cell=cell, para_env=para_env)
      CALL comm%set_handle(para_env%get_handle())
      CPASSERT(ASSOCIATED(state) .AND. state%enabled)
      IF (ALLOCATED(state%paw_reconstruction)) THEN
         IF (.NOT. xc_section_uses_native_skala_evaluator(xc_section)) THEN
            CPABORT("Inverse PAW learned XC requires the native Skala evaluator")
         END IF
      END IF
      CPASSERT(PRESENT(kinetic) .EQV. ALLOCATED(state%kinetic_cache))
      IF (PRESENT(potential)) THEN
         CPASSERT(PRESENT(kinetic) .EQV. PRESENT(kinetic_potential))
      END IF
      IF (section_get_ival(xc_section, "XC_GRID%XC_DERIV") /= xc_deriv_pw .OR. &
          section_get_ival(xc_section, "XC_GRID%XC_SMOOTH_RHO") /= xc_rho_no_smooth) THEN
         CPABORT("Inverse mixed XC currently requires XC_DERIV PW and XC_SMOOTH_RHO NONE")
      END IF
      CALL frame_refs(state, frames, augmentation)
      nk = SIZE(frames)
      nao = SIZE(frames(1)%frame%coefficients, 2)
      ns = SIZE(density)
      IF (PRESENT(kinetic)) THEN
         CPASSERT(SIZE(kinetic) == ns)
      END IF
      IF (PRESENT(kinetic_potential)) THEN
         CPASSERT(SIZE(kinetic_potential) == ns)
      END IF
      IF (.NOT. ALLOCATED(state%atom_grid%points)) CALL prepare_atom_grid(ks_env, state, frames(1)%frame, &
                                                                          xc_section_uses_native_skala_evaluator(xc_section))
      np = SIZE(state%atom_grid%weights)
      IF (ALLOCATED(state%xc_matrices)) DEALLOCATE (state%xc_matrices)
      ALLOCATE (state%xc_matrices(nao, nao, nk, ns))
      IF (.NOT. ASSOCIATED(state%bloch%kpoints)) CALL qs_rho_get(rho, rho_ao_kp=density_ao)
      DO spin = 1, ns
         IF (ASSOCIATED(state%bloch%kpoints)) THEN
            CALL qs_upf_inverse_kpoints_density(state%bloch, spin, matrix)
            state%xc_matrices(:, :, :, spin) = matrix
         ELSE
            CALL qs_upf_inverse_matrix_gather(state%frame%ao, density_ao(spin, 1)%matrix, gamma)
            state%xc_matrices(:, :, 1, spin) = gamma
         END IF
      END DO
      state%xc_grid => density(1)%pw_grid
      CALL qs_upf_inverse_xc_pw_fields(density, state%xc_fields, kinetic)
      ALLOCATE (core(np, ns), core_gradient(np, 3, ns))
      IF (PRESENT(kinetic)) ALLOCATE (core_tau(np, ns))
      CALL core_fields(ks_env, state, core, core_gradient, kinetic=core_tau)
      CALL qs_upf_inverse_xc_grid_evaluate(frames, augmentation, state%xc_matrices, state%xc_kweights, &
                                           state%atom_grid%points, state%atom_grid%weights, state%xc_fields, &
                                           state%xc_grid, cell, core, core_gradient, xc_section, ns == 2, &
                                           comm, state%xc, energy, adjoint, state%kinetic_cache, core_tau, state%atom_grid, &
                                           state%paw_reconstruction)
      IF (PRESENT(potential)) THEN
         CPASSERT(SIZE(potential) == ns)
         ALLOCATE (background_potential(ns))
         IF (PRESENT(kinetic_potential)) ALLOCATE (background_kinetic(ns))
         DO spin = 1, ns
            CALL background_potential(spin)%create(state%xc_grid)
            IF (PRESENT(kinetic_potential)) CALL background_kinetic(spin)%create(state%xc_grid)
         END DO
         CALL qs_upf_inverse_xc_pw_adjoint(adjoint, background_potential, background_kinetic)
         DO spin = 1, ns
            CPASSERT(ASSOCIATED(potential(spin)%pw_grid, state%xc_grid))
            potential(spin)%array = potential(spin)%array + background_potential(spin)%array
            CALL background_potential(spin)%release()
            IF (PRESENT(kinetic_potential)) THEN
               CPASSERT(ASSOCIATED(kinetic_potential(spin)%pw_grid, state%xc_grid))
               kinetic_potential(spin)%array = kinetic_potential(spin)%array + background_kinetic(spin)%array
               CALL background_kinetic(spin)%release()
            END IF
         END DO
      END IF
   END SUBROUTINE qs_upf_inverse_scf_xc_evaluate

! **************************************************************************************************
!> \brief Atom XC geometry response after the native folded source-density response is accounted for.
!>        The folded background potential already differentiates the native valence/core source.
!>        Here retain its spectral-gradient lattice response and the moving interpolation, original
!>        core, inverse transformation and atom partition. No PW measure Jacobian is added here.
!> \param ks_env ...
!> \param atom_response ...
!> \param cell_response ...
!> \param reduced_k_response ...
!> \param response Globally reduced additional energy gradient.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_scf_xc_response(ks_env, atom_response, cell_response, reduced_k_response, response)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      REAL(dp), INTENT(IN)                               :: atom_response(:, :), &
                                                            cell_response(3, 3), &
                                                            reduced_k_response(:, :)
      REAL(dp), INTENT(OUT)                              :: response

      COMPLEX(dp), ALLOCATABLE                           :: coefficients(:, :, :)
      INTEGER                                            :: axis, ik, nk, np, ns, other, stride
      REAL(dp)                                           :: strain(3, 3)
      REAL(dp), ALLOCATABLE :: core(:, :), core_gradient(:, :, :), core_tau(:, :), &
         defining_points(:, :), defining_weights(:), fields(:, :, :, :), points(:, :), weights(:)
      TYPE(cell_type), POINTER                           :: cell
      TYPE(mp_comm_type)                                 :: comm
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(qs_upf_inverse_augmentation_type), POINTER    :: augmentation
      TYPE(qs_upf_inverse_frame_ref_type), ALLOCATABLE   :: frames(:)
      TYPE(qs_upf_inverse_scf_type), POINTER             :: state
      TYPE(qs_upf_inverse_type), POINTER                 :: frame

      response = 0.0_dp
      CALL get_ks_env(ks_env, upf_inverse=state, cell=cell, para_env=para_env)
      CALL comm%set_handle(para_env%get_handle())
      CPASSERT(ASSOCIATED(state))
      IF (.NOT. ASSOCIATED(state%xc_grid)) RETURN
      CALL frame_refs(state, frames, augmentation)
      frame => frames(1)%frame
      nk = SIZE(frames)
      np = SIZE(state%atom_grid%weights)
      ns = SIZE(state%xc_matrices, 4)
      stride = 4
      IF (state%xc%kinetic) stride = 5
      ALLOCATE (points(3, np), weights(np), core(np, ns), core_gradient(np, 3, ns), &
                coefficients(SIZE(frame%coefficients, 1), SIZE(frame%coefficients, 2), nk))
      IF (state%xc%kinetic) ALLOCATE (core_tau(np, ns))
      CALL qs_upf_inverse_atom_grid_response(state%atom_grid, atom_response, cell_response, points, weights)
      strain = MATMUL(cell_response, cell%h_inv)
      DO ik = 1, nk
         frame => frames(ik)%frame
         ALLOCATE (defining_points(3, SIZE(frame%weights)), defining_weights(SIZE(frame%weights)))
         defining_points(:, :) = MATMUL(strain, frame%points)
         defining_weights(:) = (strain(1, 1) + strain(2, 2) + strain(3, 3))*frame%weights
         CALL qs_upf_inverse_response(frame, defining_points, atom_response, cell_response, reduced_k_response(:, ik), &
                                      defining_weights, coefficients(:, :, ik))
         DEALLOCATE (defining_points, defining_weights)
      END DO
      ALLOCATE (fields, MOLD=state%xc_fields)
      fields = 0.0_dp
      DO ik = 1, ns
         DO axis = 1, 3
            DO other = 1, 3
               fields(:, :, :, stride*(ik - 1) + 1 + axis) = fields(:, :, :, stride*(ik - 1) + 1 + axis) - &
                                                           strain(other, axis)*state%xc_fields(:, :, :, stride*(ik - 1) + 1 + other)
            END DO
         END DO
      END DO
      CALL core_fields(ks_env, state, core, core_gradient, points, atom_response, cell_response, core_tau)
      CALL qs_upf_inverse_xc_grid_response(frames, augmentation, state%xc_matrices, state%xc_kweights, &
                                           state%atom_grid%points, state%atom_grid%weights, state%xc, state%xc_fields, &
                                           fields, state%xc_grid, cell, points, weights, atom_response, cell_response, &
                                           reduced_k_response, coefficients, core, core_gradient, comm, response, &
                                           state%kinetic_cache, core_tau, state%paw_reconstruction)
   END SUBROUTINE qs_upf_inverse_scf_xc_response

! **************************************************************************************************
!> \brief Borrow each defining inverse frame in the same global k order used by the native density.
!> \param state ...
!> \param frames ...
!> \param augmentation ...
! **************************************************************************************************
   SUBROUTINE frame_refs(state, frames, augmentation)
      TYPE(qs_upf_inverse_scf_type), POINTER             :: state
      TYPE(qs_upf_inverse_frame_ref_type), ALLOCATABLE, &
         INTENT(OUT)                                     :: frames(:)
      TYPE(qs_upf_inverse_augmentation_type), POINTER    :: augmentation

      INTEGER                                            :: ik, nk

      nk = 1
      IF (ASSOCIATED(state%bloch%kpoints)) nk = SIZE(state%bloch%weights)
      ALLOCATE (frames(nk))
      IF (.NOT. ALLOCATED(state%xc_kweights)) ALLOCATE (state%xc_kweights(nk))
      IF (ASSOCIATED(state%bloch%kpoints)) THEN
         DO ik = 1, nk
            frames(ik)%frame => state%bloch%frames(ik)%ao
         END DO
         state%xc_kweights(:) = state%bloch%weights
         augmentation => state%bloch%augmentation
      ELSE
         frames(1)%frame => state%frame%ao
         state%xc_kweights = 1.0_dp
         augmentation => state%augmentation
      END IF
   END SUBROUTINE frame_refs

! **************************************************************************************************
!> \brief Reuse native radial/Lebedev grids with fixed support from basis and original radial fields.
!> \param ks_env ...
!> \param state ...
!> \param frame ...
!> \param atom_blocks Keep complete local feature blocks for a nonlocal functional.
! **************************************************************************************************
   SUBROUTINE prepare_atom_grid(ks_env, state, frame, atom_blocks)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(qs_upf_inverse_scf_type), INTENT(INOUT)       :: state
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      LOGICAL, INTENT(IN)                                :: atom_blocks

      CHARACTER(LEN=256)                                 :: core_reason
      INTEGER                                            :: atom, ikind, ll, na, nr, quadrature
      INTEGER, ALLOCATABLE                               :: kind_of(:)
      LOGICAL                                            :: core_available
      REAL(dp)                                           :: ae_radius, core_radius
      REAL(dp), ALLOCATABLE                              :: core_tau(:), radial(:, :), support(:)
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: control
      TYPE(grid_atom_type), POINTER                      :: grid
      TYPE(gto_basis_set_type), POINTER                  :: basis
      TYPE(mp_comm_type)                                 :: comm
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), POINTER                       :: atoms(:)
      TYPE(qs_kind_type), POINTER                        :: kinds(:)
      TYPE(qs_upf_inverse_atom_grid_ref_type), &
         ALLOCATABLE                                     :: grids(:)

      CALL get_ks_env(ks_env, qs_kind_set=kinds, particle_set=atoms, cell=cell, para_env=para_env, &
                      dft_control=control)
      CALL comm%set_handle(para_env%get_handle())
      quadrature = control%qs_control%gapw_control%quadrature
      ALLOCATE (grids(SIZE(kinds)), state%core_cache(SIZE(kinds)), state%core_kinetic_cache(SIZE(kinds)), &
                kind_of(SIZE(atoms)), support(SIZE(atoms)))
      CALL init_lebedev_grids()
      DO ikind = 1, SIZE(kinds)
         CALL get_qs_kind(kinds(ikind), grid_atom=grid, ngrid_ang=na, ngrid_rad=nr)
         IF (.NOT. ASSOCIATED(grid)) THEN
            CALL allocate_grid_atom(kinds(ikind)%grid_atom)
            grid => kinds(ikind)%grid_atom
         END IF
         IF (.NOT. ASSOCIATED(grid%weight)) THEN
            ll = get_number_of_lebedev_grid(n=na)
            na = lebedev_grid(ll)%n
            grid%nr = nr
            grid%ng_sphere = na
            CALL create_grid_atom(grid, nr, na, 0, ll, quadrature)
         END IF
         grids(ikind)%grid => grid
         IF (kinds(ikind)%ghost) CYCLE
         ASSOCIATE (pot => kinds(ikind)%sgp_potential%upf_data)
            IF (ALLOCATED(pot%rho_nlcc)) THEN
               IF (ANY(pot%rho_nlcc /= 0.0_dp)) THEN
                  ALLOCATE (radial(SIZE(pot%r), 1))
                  radial(:, 1) = SQRT(fourpi)*pot%rho_nlcc
                  CALL upf_radial_fields_init(pot%r, radial, [0], state%core_cache(ikind))
                  DEALLOCATE (radial)
               END IF
            END IF
            IF (control%use_kinetic_energy_density) THEN
               IF (state%paw(ikind)%active) THEN
                  IF (ALLOCATED(core_tau)) DEALLOCATE (core_tau)
                  ALLOCATE (core_tau, SOURCE=state%paw(ikind)%xcgrid%core_kinetic%kinetic(:, 2))
               ELSE
                  CALL upf_nlcc_kinetic(pot, core_tau, core_available, core_reason)
                  IF (.NOT. core_available) CPABORT(TRIM(core_reason))
               END IF
               IF (ANY(core_tau /= 0.0_dp)) THEN
                  ALLOCATE (radial(SIZE(pot%r), 1))
                  radial(:, 1) = SQRT(fourpi)*core_tau
                  CALL upf_radial_fields_init(pot%r, radial, [0], state%core_kinetic_cache(ikind))
                  DEALLOCATE (radial)
               END IF
            END IF
         END ASSOCIATE
      END DO
      CALL deallocate_lebedev_grids()
      core_radius = MAX(MAXVAL(state%core_cache%radius), MAXVAL(state%core_kinetic_cache%radius))
      DO atom = 1, SIZE(atoms)
         ikind = atoms(atom)%atomic_kind%kind_number
         kind_of(atom) = ikind
         CALL get_qs_kind(kinds(ikind), basis_set=basis, basis_type='ORB')
         support(atom) = MAX(basis%kind_radius, frame%beta_cache(ikind)%radius, &
                             frame%delta_cache(ikind)%radius, state%core_cache(ikind)%radius, &
                             state%core_kinetic_cache(ikind)%radius)
         IF (ALLOCATED(state%paw_reconstruction) .AND. state%paw(ikind)%active) THEN
            ae_radius = state%paw(ikind)%xcgrid%core_kinetic%ae_support_radius
            IF (state%paw(ikind)%xcgrid%has_ae_core) THEN
               ae_radius = MAX(ae_radius, state%paw(ikind)%xcgrid%r(state%paw(ikind)%xcgrid%core_support_index))
            END IF
            support(atom) = MAX(support(atom), ae_radius)
            core_radius = MAX(core_radius, ae_radius)
         END IF
         IF (ASSOCIATED(state%bloch%kpoints)) THEN
            support(atom) = MAX(support(atom), state%bloch%augmentation%kinds(ikind)%realspace_radius)
         ELSE
            support(atom) = MAX(support(atom), state%augmentation%kinds(ikind)%realspace_radius)
         END IF
         support(atom) = support(atom)/0.75_dp
      END DO
      CALL qs_upf_inverse_atom_grid_prepare(frame%positions, cell, kind_of, grids, support, comm, state%atom_grid, &
                                            atom_blocks=atom_blocks)
      ! The core image bounds also cover atom rows outside the defining physical-cell grid.
      CALL prepare_core_images(frame%positions, support, cell, &
                               core_radius, state%core_images)
   END SUBROUTINE prepare_atom_grid

! **************************************************************************************************
!> \brief Enumerate core images reaching any supported atom row, at fixed topology during responses.
!> \param positions ...
!> \param support ...
!> \param cell ...
!> \param radius ...
!> \param images ...
! **************************************************************************************************
   SUBROUTINE prepare_core_images(positions, support, cell, radius, images)
      REAL(dp), INTENT(IN)                               :: positions(:, :), support(:)
      TYPE(cell_type), POINTER                           :: cell
      REAL(dp), INTENT(IN)                               :: radius
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: images(:, :)

      INTEGER                                            :: axis, i, j, k, lower(3), row, upper(3)
      REAL(dp)                                           :: reach
      REAL(dp), ALLOCATABLE                              :: fractional(:, :)

      ALLOCATE (fractional(3, SIZE(positions, 2)))
      fractional(:, :) = MATMUL(cell%h_inv, positions)
      DO axis = 1, 3
         reach = (MAXVAL(support) + radius)*SQRT(SUM(cell%h_inv(axis, :)**2))
         upper(axis) = CEILING(MAXVAL(fractional(axis, :)) - MINVAL(fractional(axis, :)) + reach)
         lower(axis) = -upper(axis)
      END DO
      ALLOCATE (images(3, PRODUCT(upper - lower + 1)))
      row = 0
      DO k = lower(3), upper(3)
         DO j = lower(2), upper(2)
            DO i = lower(1), upper(1)
               row = row + 1
               images(:, row) = [i, j, k]
            END DO
         END DO
      END DO
   END SUBROUTINE prepare_core_images

! **************************************************************************************************
!> \brief Original frozen-core point fields or their complete moving-point geometry response.
!> \param ks_env ...
!> \param state ...
!> \param density ...
!> \param gradient ...
!> \param point_response ...
!> \param atom_response ...
!> \param cell_response ...
!> \param kinetic Optional independently supplied core kinetic density or its geometry response.
! **************************************************************************************************
   SUBROUTINE core_fields(ks_env, state, density, gradient, point_response, atom_response, cell_response, kinetic)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(qs_upf_inverse_scf_type), INTENT(IN)          :: state
      REAL(dp), INTENT(OUT)                              :: density(:, :), gradient(:, :, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: point_response(:, :), &
                                                            atom_response(:, :), &
                                                            cell_response(3, 3)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: kinetic(:, :)

      COMPLEX(dp), ALLOCATABLE                           :: dg(:, :, :), dv(:, :), spatial(:, :, :), &
                                                            value(:, :)
      INTEGER                                            :: atom, ikind, image, np, point, spin
      LOGICAL                                            :: reconstruct_core
      REAL(dp)                                           :: displacement(3), motion(3), rho, &
                                                            rho_gradient(3), rho_hessian(3, 3), &
                                                            tau, tau_gradient(3)
      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_type), POINTER                       :: atoms(:)

      CALL get_ks_env(ks_env, particle_set=atoms, cell=cell)
      np = SIZE(density, 1)
      density = 0.0_dp
      gradient = 0.0_dp
      IF (PRESENT(kinetic)) THEN
         CPASSERT(ALL(SHAPE(kinetic) == SHAPE(density)))
         kinetic = 0.0_dp
      END IF
      ALLOCATE (value(np, 1), spatial(np, 1, 3), dv(np, 1), dg(np, 1, 3))
      DO atom = 1, SIZE(atoms)
         ikind = atoms(atom)%atomic_kind%kind_number
         reconstruct_core = ALLOCATED(state%paw_reconstruction) .AND. state%paw(ikind)%active
         IF (state%core_cache(ikind)%nproj == 0 .AND. .NOT. reconstruct_core) CYCLE
         IF (reconstruct_core) THEN
            value = CMPLX(0.0_dp, 0.0_dp, dp)
            spatial = CMPLX(0.0_dp, 0.0_dp, dp)
            DO image = 1, SIZE(state%core_images, 2)
               DO point = 1, np
                  displacement = state%atom_grid%points(:, point) - atoms(atom)%r - &
                                 MATMUL(cell%hmat, REAL(state%core_images(:, image), dp))
                  CALL paw_xc_ae_core_cartesian(state%paw(ikind)%xcgrid, displacement, rho, rho_gradient, rho_hessian)
                  IF (PRESENT(point_response)) THEN
                     CPASSERT(PRESENT(atom_response) .AND. PRESENT(cell_response))
                     motion = point_response(:, point) - atom_response(:, atom) - &
                              MATMUL(cell_response, REAL(state%core_images(:, image), dp))
                     rho = DOT_PRODUCT(rho_gradient, motion)
                     rho_gradient = MATMUL(rho_hessian, motion)
                  END IF
                  value(point, 1) = value(point, 1) + rho
                  spatial(point, 1, :) = spatial(point, 1, :) + rho_gradient
               END DO
            END DO
         ELSE IF (PRESENT(point_response)) THEN
            CPASSERT(PRESENT(atom_response) .AND. PRESENT(cell_response))
            CALL upf_radial_fields_periodic(state%core_cache(ikind), state%atom_grid%points, atoms(atom)%r, &
                                            cell%hmat, state%core_images, [0.0_dp, 0.0_dp, 0.0_dp], value, spatial, dv, dg, &
                                            point_response, atom_response(:, atom), cell_response)
            value(:, :) = dv
            spatial(:, :, :) = dg
         ELSE
            CALL upf_radial_fields_periodic(state%core_cache(ikind), state%atom_grid%points, atoms(atom)%r, &
                                            cell%hmat, state%core_images, [0.0_dp, 0.0_dp, 0.0_dp], value, spatial)
         END IF
         DO spin = 1, SIZE(density, 2)
            density(:, spin) = density(:, spin) + REAL(value(:, 1), dp)/REAL(SIZE(density, 2), dp)
            gradient(:, :, spin) = gradient(:, :, spin) + REAL(spatial(:, 1, :), dp)/REAL(SIZE(density, 2), dp)
         END DO
      END DO
      IF (.NOT. PRESENT(kinetic)) RETURN
      DO atom = 1, SIZE(atoms)
         ikind = atoms(atom)%atomic_kind%kind_number
         reconstruct_core = ALLOCATED(state%paw_reconstruction) .AND. state%paw(ikind)%active
         IF (state%core_kinetic_cache(ikind)%nproj == 0 .AND. .NOT. reconstruct_core) CYCLE
         IF (state%paw(ikind)%active) THEN
            value = CMPLX(0.0_dp, 0.0_dp, dp)
            DO image = 1, SIZE(state%core_images, 2)
               DO point = 1, np
                  displacement = state%atom_grid%points(:, point) - atoms(atom)%r - &
                                 MATMUL(cell%hmat, REAL(state%core_images(:, image), dp))
                  IF (reconstruct_core) THEN
                     CALL paw_core_kinetic_ae_cartesian(state%paw(ikind)%xcgrid%core_kinetic, displacement, tau, tau_gradient)
                  ELSE
                     CALL paw_core_kinetic_cartesian(state%paw(ikind)%xcgrid%core_kinetic, displacement, tau, tau_gradient)
                  END IF
                  IF (PRESENT(point_response)) THEN
                     CPASSERT(PRESENT(atom_response) .AND. PRESENT(cell_response))
                     motion = point_response(:, point) - atom_response(:, atom) - &
                              MATMUL(cell_response, REAL(state%core_images(:, image), dp))
                     tau = DOT_PRODUCT(tau_gradient, motion)
                  END IF
                  value(point, 1) = value(point, 1) + tau
               END DO
            END DO
         ELSE IF (PRESENT(point_response)) THEN
            CPASSERT(PRESENT(atom_response) .AND. PRESENT(cell_response))
            CALL upf_radial_fields_periodic(state%core_kinetic_cache(ikind), state%atom_grid%points, atoms(atom)%r, &
                                            cell%hmat, state%core_images, [0.0_dp, 0.0_dp, 0.0_dp], value, spatial, dv, dg, &
                                            point_response, atom_response(:, atom), cell_response)
            value(:, :) = dv
         ELSE
            CALL upf_radial_fields_periodic(state%core_kinetic_cache(ikind), state%atom_grid%points, atoms(atom)%r, &
                                            cell%hmat, state%core_images, [0.0_dp, 0.0_dp, 0.0_dp], value, spatial)
         END IF
         DO spin = 1, SIZE(kinetic, 2)
            kinetic(:, spin) = kinetic(:, spin) + REAL(value(:, 1), dp)/REAL(SIZE(kinetic, 2), dp)
         END DO
      END DO
   END SUBROUTINE core_fields

END MODULE qs_upf_inverse_scf_xc
