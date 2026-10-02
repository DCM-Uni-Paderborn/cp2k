!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM skala_sirius_unittest
   USE cell_methods, ONLY: cell_create
   USE cell_types, ONLY: cell_release, &
                         cell_type
   USE ieee_arithmetic, ONLY: ieee_is_finite
   USE input_cp2k_xc, ONLY: create_xc_section
   USE input_section_types, ONLY: section_release, &
                                  section_type, &
                                  section_vals_create, &
                                  section_vals_release, &
                                  section_vals_type, &
                                  section_vals_val_set
   USE kinds, ONLY: default_path_length, &
                    default_string_length, &
                    dp
   USE libcp2k, ONLY: cp2k_finalize, &
                      cp2k_init
   USE message_passing, ONLY: mp_comm_type, &
                              mp_comm_world
   USE particle_types, ONLY: particle_type
   USE sirius, ONLY: &
      sirius_add_atom, sirius_add_atom_type, sirius_add_atom_type_radial_function, &
      sirius_add_atom_type_aw_descriptor, sirius_add_atom_type_lo_descriptor, &
      sirius_context_handler, sirius_create_context, sirius_create_ground_state, &
      sirius_create_kset, sirius_find_ground_state, sirius_free_handler, &
      sirius_generate_effective_potential, sirius_get_energy, sirius_ground_state_handler, &
      sirius_import_parameters, sirius_initialize_context, sirius_kpoint_set_handler, &
      sirius_set_atom_type_dion, sirius_set_atom_type_paw, sirius_set_atom_type_radial_grid, &
      sirius_set_atom_type_radial_grid_inf, &
      sirius_set_lattice_vectors, sirius_set_parameters
#if defined(__SIRIUS_XC_DERIVATIVE)
   USE sirius, ONLY: sirius_check_external_xc_derivative
#endif
   USE sirius_interface, ONLY: cp_sirius_init
   USE skala_sirius, ONLY: skala_sirius_register_points, &
                           skala_sirius_release, &
                           skala_sirius_state_type

   IMPLICIT NONE

   INTEGER, PARAMETER :: nrad = 400, nmt = 801, nfree = 3000
   TYPE(mp_comm_type) :: comm
   TYPE(cell_type), POINTER :: cell => NULL()
   TYPE(particle_type), POINTER :: atoms(:) => NULL()
   TYPE(section_type), POINTER :: definition => NULL()
   TYPE(section_vals_type), POINTER :: xc => NULL()
   TYPE(skala_sirius_state_type), POINTER :: state => NULL()
   TYPE(sirius_context_handler) :: ctx
   TYPE(sirius_kpoint_set_handler) :: ks
   TYPE(sirius_ground_state_handler) :: gs
   CHARACTER(len=default_path_length) :: model, paw_setup
   CHARACTER(len=32) :: side_arg, step_arg
   INTEGER :: i, l, mode, mag, nr, nb, x, y, z, ierr, side, first_mode, last_mode
   REAL(dp) :: lattice(3, 3), r(nrad), v(nrad), rho(nrad), orbital(nrad), dion(1, 1), &
               energy, exc, repeat_exc, kpoint(3, 1), reference(3, 2), fd_step, &
               rmt(nmt), rfree(nfree), rhofree(nfree)
#if defined(__SIRIUS_XC_DERIVATIVE)
   REAL(dp) :: derivatives(3, 2), errors(2, 2)
#endif
   REAL(dp), ALLOCATABLE :: points(:, :), weights(:), descriptor_weights(:), centres(:, :)
   INTEGER, ALLOCATABLE :: sizes(:)
   LOGICAL :: converged

   CALL get_command_argument(1, model)
   CALL get_command_argument(2, paw_setup)
   CALL get_command_argument(3, side_arg)
   CALL get_command_argument(4, step_arg)
   fd_step = 0.0_dp
   IF (LEN_TRIM(step_arg) > 0) THEN
      READ (step_arg, *, IOSTAT=ierr) fd_step
      IF (ierr /= 0) ERROR STOP "Invalid orbital derivative step"
      IF (.NOT. ieee_is_finite(fd_step) .OR. fd_step <= 0.0_dp) ERROR STOP "Orbital derivative step must be positive"
#if !defined(__SIRIUS_XC_DERIVATIVE)
      ERROR STOP "SIRIUS lacks the external XC orbital derivative diagnostic"
#endif
   END IF
   side = 6
   IF (LEN_TRIM(side_arg) > 0) THEN
      READ (side_arg, *, IOSTAT=ierr) side
      IF (ierr /= 0) ERROR STOP "Invalid quadrature dimension"
   END IF
   IF (side < 4 .OR. side > 64) ERROR STOP "Quadrature dimension must be between 4 and 64"
   first_mode = 0
   last_mode = 1
   IF (LEN_TRIM(paw_setup) > 0) first_mode = 1
   IF (TRIM(paw_setup) == "LAPW") THEN
      first_mode = 2
      last_mode = 2
   END IF
   reference = 0.0_dp
   IF (LEN_TRIM(model) == 0 .OR. LEN_TRIM(model) > default_string_length) &
      ERROR STOP "Pass the test model path, or SKALA with GAUXC_SKALA_MODEL set"
   CALL cp2k_init()
   CALL cp_sirius_init()
   comm = mp_comm_world
   CALL comm%init()
   IF (comm%mepos == 0) WRITE (*, '(A,I0,2A)') "Points per block: ", side**3, "; PAW setup: ", TRIM(paw_setup)
   CALL create_xc_section(definition)
   CALL section_vals_create(xc, definition)
   CALL section_vals_val_set(xc, "XC_FUNCTIONAL%GAUXC%MODEL", c_val=TRIM(model))
   CALL section_vals_val_set(xc, "XC_FUNCTIONAL%GAUXC%NATIVE_GRID", l_val=.TRUE.)
   lattice = 0.0_dp
   lattice(1, 1) = 8.0_dp
   lattice(2, 2) = 9.0_dp
   lattice(3, 3) = 10.0_dp
   CALL cell_create(cell, hmat=lattice, periodic=[1, 1, 1])
   ALLOCATE (atoms(1))
   atoms(1)%r = [1.6_dp, 2.7_dp, 4.0_dp]
   nr = 0
   nb = 0
   IF (comm%mepos == 0) THEN
      nr = side**3
      nb = 1
   END IF
   ALLOCATE (points(3, nr), weights(nr), descriptor_weights(nr), centres(3, nb), sizes(nb))
   IF (nb > 0) THEN
      i = 0
      DO z = 0, side - 1
         DO y = 0, side - 1
            DO x = 0, side - 1
               i = i + 1
               points(:, i) = [8*(0.037_dp + REAL(x, dp)/side), 9*(0.029_dp + REAL(y, dp)/side), &
                               10*(0.041_dp + REAL(z, dp)/side)]
            END DO
         END DO
      END DO
      sizes = nr
      centres(:, 1) = atoms(1)%r
   END IF
   weights = cell%deth/REAL(side**3, dp)
   descriptor_weights = weights
   IF (nr > 0) weights(1) = 0.0_dp
   kpoint(:, 1) = [0.07_dp, 0.03_dp, 0.02_dp]

   DO mode = first_mode, last_mode
      DO mag = 0, 1
         CALL sirius_create_context(comm%get_handle(), ctx)
         IF (mode == 2) THEN
            CALL sirius_import_parameters(ctx, &
                                          '{"parameters":{"electronic_structure_method":"full_potential_lapwlo",'// &
                                          '"use_symmetry":false,"pw_cutoff":5.0,"gk_cutoff":1.5,"num_bands":4,'// &
                                          '"xc_functionals":[],"smearing_width":0.01,"auto_rmt":0,'// &
                                          '"core_relativity":"none","valence_relativity":"none",'// &
                                          '"lmax_apw":3,"lmax_rho":6,"lmax_pot":6},'// &
                                          '"control":{"mpi_grid_dims":[1,1],"processing_unit":"CPU","verbosity":1},'// &
                                          '"mixer":{"type":"linear","beta":0.2,"use_hartree":false},'// &
                                          '"iterative_solver":{"type":"exact","min_tolerance":1e-11}}')
         ELSE
            CALL sirius_import_parameters(ctx, &
                                          '{"parameters":{"electronic_structure_method":"pseudopotential",'// &
                                          '"use_symmetry":false,"pw_cutoff":8.0,"gk_cutoff":2.6,"num_bands":4,'// &
                                          '"xc_functionals":[],"smearing_width":0.01},'// &
                                          '"control":{"mpi_grid_dims":[1,1],"processing_unit":"CPU","verbosity":1},'// &
                                          '"mixer":{"type":"linear","beta":0.2,"use_hartree":false},'// &
                                          '"iterative_solver":{"type":"exact","min_tolerance":1e-11,"num_steps":80}}')
         END IF
         CALL sirius_set_parameters(ctx, num_mag_dims=mag)
         CALL sirius_set_lattice_vectors(ctx, lattice(:, 1), lattice(:, 2), lattice(:, 3))
         IF (mode == 2) THEN
            CALL sirius_add_atom_type(ctx, "He", zn=2)
            DO i = 1, nmt
               rmt(i) = 1.0e-6_dp + (1.8_dp - 1.0e-6_dp)*(REAL(i - 1, dp)/REAL(nmt - 1, dp))**2
            END DO
            DO i = 1, nfree
               rfree(i) = 1.0e-7_dp + (20.0_dp - 1.0e-7_dp)*(REAL(i - 1, dp)/REAL(nfree - 1, dp))**3
            END DO
            rhofree = 16.0_dp/ACOS(-1.0_dp)*EXP(-4*rfree)
            CALL sirius_set_atom_type_radial_grid(ctx, "He", nmt, rmt)
            CALL sirius_set_atom_type_radial_grid_inf(ctx, "He", nfree, rfree)
            CALL sirius_add_atom_type_radial_function(ctx, "He", "ae_rho", rhofree, nfree)
            DO l = 0, 3
               CALL sirius_add_atom_type_aw_descriptor(ctx, "He", l + 1, l, -0.3_dp, 0, .FALSE.)
               CALL sirius_add_atom_type_aw_descriptor(ctx, "He", l + 1, l, -0.3_dp, 1, .FALSE.)
            END DO
            CALL sirius_add_atom_type_lo_descriptor(ctx, "He", 1, 1, 0, -0.7_dp, 0, .FALSE.)
            CALL sirius_add_atom_type_lo_descriptor(ctx, "He", 1, 1, 0, -0.7_dp, 1, .FALSE.)
         ELSE IF (mode == 1 .AND. LEN_TRIM(paw_setup) > 0) THEN
            CALL sirius_add_atom_type(ctx, "He", fname=TRIM(paw_setup))
         ELSE
            CALL sirius_add_atom_type(ctx, "He", zn=2)
            DO i = 1, nrad
               r(i) = 20.0_dp*(REAL(i - 1, dp)/REAL(nrad - 1, dp))**2
               IF (mode == 1) r(i) = 0.001_dp + 2.999_dp*(REAL(i - 1, dp)/REAL(nrad - 1, dp))**2
            END DO
            v = -2.0_dp/SQRT(1.0_dp + r**2)
            rho = 4.0_dp*EXP(-r**2)*r
            orbital = EXP(-r)
            CALL sirius_set_atom_type_radial_grid(ctx, "He", nrad, r)
            CALL sirius_add_atom_type_radial_function(ctx, "He", "vloc", v, nrad)
            CALL sirius_add_atom_type_radial_function(ctx, "He", "ps_rho_total", rho, nrad)
            CALL sirius_add_atom_type_radial_function(ctx, "He", "ps_atomic_wf", orbital, nrad, n=1, l=0)
            dion = 0.0_dp
            IF (mode == 1) THEN
               ! Deliberately synthetic setup: tests the coupled implementation, not PAW accuracy.
               CALL sirius_add_atom_type_radial_function(ctx, "He", "beta", 0.2_dp*r*EXP(-2*r**2), nrad, l=0)
               CALL sirius_add_atom_type_radial_function(ctx, "He", "ae_paw_wf", r*EXP(-r**2), nrad)
               CALL sirius_add_atom_type_radial_function(ctx, "He", "ps_paw_wf", r*EXP(-0.8_dp*r**2), nrad)
               CALL sirius_add_atom_type_radial_function(ctx, "He", "q_aug", 0.05_dp*r**2*EXP(-2*r**2), nrad, &
                                                         l=0, idxrf1=1, idxrf2=1)
               CALL sirius_add_atom_type_radial_function(ctx, "He", "ae_paw_core", 0.02_dp*EXP(-3*r**2), nrad)
               CALL sirius_add_atom_type_radial_function(ctx, "He", "ae_paw_core_tau", 0.03_dp*EXP(-3*r**2), nrad)
               CALL sirius_add_atom_type_radial_function(ctx, "He", "ps_rho_core", 0.01_dp*EXP(-2*r**2), nrad)
               CALL sirius_set_atom_type_dion(ctx, "He", 1, dion)
               CALL sirius_set_atom_type_paw(ctx, "He", 0.0_dp, [2.0_dp], 1)
            ELSE
               CALL sirius_set_atom_type_dion(ctx, "He", 0, dion(:0, :0))
            END IF
         END IF
         CALL sirius_add_atom(ctx, "He", [0.2_dp, 0.3_dp, 0.4_dp])
         CALL skala_sirius_register_points(ctx, xc, atoms, cell, points, weights, descriptor_weights, sizes, centres, state)
         CALL sirius_initialize_context(ctx)
         CALL sirius_create_kset(ctx, 1, kpoint, [1.0_dp], .TRUE., ks)
         CALL sirius_create_ground_state(ks, gs)
         converged = .FALSE.
         CALL sirius_find_ground_state(gs, density_tol=2.0e-6_dp, energy_tol=2.0e-7_dp, iter_solver_tol=1.0e-9_dp, &
                                       max_niter=150, save_state=.FALSE., converged=converged, error_code=ierr)
         IF (ierr /= 0) ERROR STOP "SIRIUS Skala SCF failed"
         CALL sirius_get_energy(gs, "total", energy)
         CALL sirius_get_energy(gs, "exc", exc)
         IF (comm%mepos == 0) WRITE (*, '(A,2I3,L3,2ES22.12)') "Point SCF (method 0=PW 1=PAW 2=LAPW,mag,converged,E,Exc): ", &
            mode, mag, converged, energy, exc
         IF (.NOT. converged .OR. .NOT. ieee_is_finite(energy) .OR. .NOT. ieee_is_finite(exc)) &
            ERROR STOP "SIRIUS Skala point SCF did not converge with finite energies"
         IF (ABS(exc) < 1.0e-6_dp) ERROR STOP "The model must produce a nontrivial XC energy"
         reference(mode + 1, mag + 1) = energy
         CALL sirius_generate_effective_potential(gs)
         CALL sirius_get_energy(gs, "exc", repeat_exc)
         IF (ABS(exc - repeat_exc) > 1.0e-7_dp) ERROR STOP "Unchanged density changed the point XC energy"
#if defined(__SIRIUS_XC_DERIVATIVE)
         IF (fd_step > 0.0_dp) THEN
            CALL sirius_check_external_xc_derivative(gs, fd_step, derivatives, error_code=ierr)
            IF (ierr /= 0) ERROR STOP "SIRIUS Skala orbital derivative check failed"
            DO i = 1, 2
               errors(:, i) = ABS(derivatives(2:3, i) - derivatives(1, i))
               IF (comm%mepos == 0) WRITE (*, '(A,I0,6ES19.10)') "XC orbital derivative (direction,h,H,FDh,FDh/2,errors): ", &
                  i, fd_step, derivatives(:, i), errors(:, i)
            END DO
            IF (.NOT. ALL(ieee_is_finite(derivatives))) ERROR STOP "Nonfinite XC orbital derivatives"
            ! Float32 model energy roundoff is amplified by division by the finite-difference step.
            IF (MAXVAL(errors) > 2.0e-5_dp) ERROR STOP "XC orbital derivative mismatch"
            IF (MAXVAL(ABS(derivatives(1, :))) < 1.0e-5_dp) ERROR STOP "Vacuous orbital derivative check"
            CALL sirius_get_energy(gs, "exc", repeat_exc)
            IF (exc /= repeat_exc) ERROR STOP "Orbital diagnostic changed the ground-state XC energy"
            CALL sirius_generate_effective_potential(gs)
            CALL sirius_get_energy(gs, "exc", repeat_exc)
            IF (ABS(exc - repeat_exc) > 1.0e-7_dp) ERROR STOP "Orbital diagnostic did not restore the callback"
         END IF
#endif
         CALL sirius_free_handler(gs)
         CALL sirius_free_handler(ks)
         CALL sirius_free_handler(ctx)
         CALL skala_sirius_release(state)
         IF (ASSOCIATED(state)) ERROR STOP "Callback state was not released"
      END DO
   END DO
   IF (ANY(ABS(reference(:, 1) - reference(:, 2)) > 2.0e-5_dp)) &
      ERROR STOP "Scalar and unpolarized collinear SCF energies differ"
   DEALLOCATE (points, weights, descriptor_weights, centres, sizes, atoms)
   CALL cell_release(cell)
   CALL section_vals_release(xc)
   CALL section_release(definition)
   CALL cp2k_finalize()
END PROGRAM skala_sirius_unittest
