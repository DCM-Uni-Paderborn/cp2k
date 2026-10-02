!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Connect SIRIUS primitive XC fields to the CP2K Skala evaluator.
! **************************************************************************************************
MODULE skala_sirius
   USE ISO_C_BINDING, ONLY: C_ASSOCIATED, C_DOUBLE, C_F_POINTER, C_FUNLOC, C_INT, C_LOC, C_PTR
   USE cell_types, ONLY: cell_type
   USE input_section_types, ONLY: section_vals_type, section_vals_val_get
   USE kinds, ONLY: dp
   USE message_passing, ONLY: mp_comm_type
   USE particle_types, ONLY: particle_type
   USE skala_atom_grid, ONLY: skala_atom_grid_build, skala_atom_grid_pullback
   USE skala_gpw_features, ONLY: skala_grid_layout_type
   USE skala_gpw_functional, ONLY: skala_grid_eval, skala_point_eval
#if defined(__SIRIUS)
   USE sirius, ONLY: sirius_context_handler
#if defined(__SIRIUS_XC_CALLBACK) && defined(__LIBTORCH)
   USE sirius, ONLY: sirius_set_xc_callback
#endif
#if defined(__SIRIUS_XC_POINT_CALLBACK) && defined(__LIBTORCH)
   USE sirius, ONLY: sirius_set_xc_point_callback
#endif
#if defined(__SIRIUS_XC_POINT_GEOMETRY) && defined(__LIBTORCH)
   USE sirius, ONLY: sirius_set_xc_point_geometry_callbacks
#endif
#endif
#include "./base/base_uses.f90"
   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'skala_sirius'
   PUBLIC :: skala_sirius_state_type, skala_sirius_release, skala_sirius_update
#if defined(__SIRIUS)
   PUBLIC :: skala_sirius_register, skala_sirius_register_points, skala_sirius_register_atom_grid
#endif

   TYPE, BIND(C) :: callback_token_type
      INTEGER(C_INT) :: id = 0
   END TYPE

   TYPE skala_sirius_state_type
      PRIVATE
      TYPE(callback_token_type) :: token
      TYPE(cell_type), POINTER :: cell => NULL()
      TYPE(particle_type), POINTER :: particles(:) => NULL()
      TYPE(section_vals_type), POINTER :: xc_section => NULL()
      INTEGER :: nradial = 0, nangular = 0, shell = 0
      REAL(dp) :: scale = 0.0_dp
      REAL(dp) :: grid_lattice(3, 3) = 0.0_dp
      REAL(dp), ALLOCATABLE :: grid_centres(:, :)
      TYPE(skala_sirius_state_type), POINTER :: next => NULL()
   END TYPE

   ! Tokens cross the C ABI; the Fortran pointers remain in this per-context registry.
   TYPE(skala_sirius_state_type), POINTER, SAVE :: states => NULL()
   INTEGER(C_INT), SAVE :: next_id = 0

CONTAINS

#if defined(__SIRIUS)
! **************************************************************************************************
!> \brief Construct a fixed periodic atom grid and register the joint field callback.
! **************************************************************************************************
   SUBROUTINE skala_sirius_register_atom_grid(context, xc_section, particles, cell, group, grid_section, state)
      TYPE(sirius_context_handler), INTENT(IN)           :: context
      TYPE(section_vals_type), POINTER                   :: xc_section
      TYPE(particle_type), POINTER                       :: particles(:)
      TYPE(cell_type), POINTER                           :: cell
      TYPE(mp_comm_type), INTENT(IN)                     :: group
      TYPE(section_vals_type), POINTER                   :: grid_section
      TYPE(skala_sirius_state_type), POINTER             :: state

      INTEGER                                            :: a, error, na, nr, shell
      INTEGER, ALLOCATABLE                               :: sizes(:)
      REAL(dp)                                           :: scale
      REAL(dp), ALLOCATABLE                              :: block_centres(:, :), centres(:, :), &
                                                            descriptor(:), energy(:), points(:, :)

      IF (ANY(cell%perd /= 1)) THEN
         CPABORT("SIRIUS Skala atom quadrature requires a three-dimensional periodic cell")
      END IF
      CALL section_vals_val_get(grid_section, "N_RADIAL", i_val=nr)
      CALL section_vals_val_get(grid_section, "N_ANGULAR", i_val=na)
      CALL section_vals_val_get(grid_section, "IMAGE_SHELL", i_val=shell)
      CALL section_vals_val_get(grid_section, "RADIAL_SCALE", r_val=scale)
      ALLOCATE (centres(3, SIZE(particles)))
      DO a = 1, SIZE(particles)
         centres(:, a) = particles(a)%r
      END DO
      CALL skala_atom_grid_build(group, cell%hmat, centres, nr, na, shell, scale, &
                                 points, energy, descriptor, sizes, block_centres)
      CALL skala_sirius_register_points(context, xc_section, particles, cell, points, energy, descriptor, &
                                        sizes, block_centres, state)
      state%nradial = nr
      state%nangular = na
      state%shell = shell
      state%scale = scale
      state%grid_lattice = cell%hmat
      state%grid_centres = centres
#if defined(__SIRIUS_XC_POINT_GEOMETRY) && defined(__LIBTORCH)
      CALL sirius_set_xc_point_geometry_callbacks(context, C_FUNLOC(evaluate_skala_point_geometry), &
                                                  C_FUNLOC(pullback_skala_quadrature), C_LOC(state%token), error)
      IF (error /= 0) CPABORT("Could not register the SIRIUS Skala geometry callbacks")
#else
      MARK_USED(error)
#endif
   END SUBROUTINE

! **************************************************************************************************
!> \brief Register Skala before initializing the SIRIUS context.
!> \param context SIRIUS context.
!> \param xc_section XC input, owned by the PWDFT environment.
!> \param particles Particle set, owned by the PWDFT subsystem.
!> \param cell Cell, owned by the PWDFT subsystem.
!> \param state Callback state, released after destroying the SIRIUS context.
! **************************************************************************************************
   SUBROUTINE skala_sirius_register(context, xc_section, particles, cell, state)
      TYPE(sirius_context_handler), INTENT(IN) :: context
      TYPE(section_vals_type), POINTER :: xc_section
      TYPE(particle_type), POINTER :: particles(:)
      TYPE(cell_type), POINTER :: cell
      TYPE(skala_sirius_state_type), POINTER :: state
#if defined(__SIRIUS_XC_CALLBACK) && defined(__LIBTORCH)
      INTEGER :: error

      CALL create_state(xc_section, particles, cell, state)
      CALL sirius_set_xc_callback(context, C_FUNLOC(evaluate_skala_fields), C_LOC(state%token), error)
      IF (error /= 0) CPABORT("Could not register the SIRIUS Skala callback")
#else
      MARK_USED(context)
      MARK_USED(xc_section)
      MARK_USED(particles)
      MARK_USED(cell)
      MARK_USED(state)
      CPABORT("Skala with SIRIUS requires LibTorch and SIRIUS with the primitive-field XC callback API")
#endif
   END SUBROUTINE

! **************************************************************************************************
!> \brief Register a fixed atom-block quadrature for joint reconstructed Skala fields.
!> \param context SIRIUS context, before initialization.
!> \param xc_section XC input, owned by the PWDFT environment.
!> \param particles Particle set, owned by the PWDFT subsystem.
!> \param cell Cell, owned by the PWDFT subsystem.
!> \param points Cartesian points in bohr, owned once across the context communicator.
!> \param energy_weights Energy-partition quadrature weights.
!> \param descriptor_weights Nonlocal descriptor-domain quadrature weights.
!> \param atom_sizes Number of points in each complete local atom block.
!> \param atom_coords Cartesian centres of these blocks in bohr.
!> \param state Callback state, released after destroying the SIRIUS context.
! **************************************************************************************************
   SUBROUTINE skala_sirius_register_points(context, xc_section, particles, cell, points, energy_weights, &
                                           descriptor_weights, atom_sizes, atom_coords, state)
      TYPE(sirius_context_handler), INTENT(IN) :: context
      TYPE(section_vals_type), POINTER :: xc_section
      TYPE(particle_type), POINTER :: particles(:)
      TYPE(cell_type), POINTER :: cell
      REAL(dp), INTENT(IN) :: points(:, :), energy_weights(:), descriptor_weights(:), atom_coords(:, :)
      INTEGER, INTENT(IN) :: atom_sizes(:)
      TYPE(skala_sirius_state_type), POINTER :: state
#if defined(__SIRIUS_XC_POINT_CALLBACK) && defined(__LIBTORCH)
      INTEGER :: nr, nb, error

      nr = SIZE(energy_weights)
      nb = SIZE(atom_sizes)
      CPASSERT(ALL(SHAPE(points) == [3, nr]))
      CPASSERT(ALL(SHAPE(atom_coords) == [3, nb]))
      CPASSERT(SIZE(descriptor_weights) == nr)
      CALL create_state(xc_section, particles, cell, state)
      CALL sirius_set_xc_point_callback(context, C_FUNLOC(evaluate_skala_points), C_LOC(state%token), nr, nb, &
                                        points, energy_weights, descriptor_weights, atom_sizes, atom_coords, error)
      IF (error /= 0) CPABORT("Could not register the SIRIUS Skala point callback")
#else
      MARK_USED(context)
      MARK_USED(xc_section)
      MARK_USED(particles)
      MARK_USED(cell)
      MARK_USED(points)
      MARK_USED(energy_weights)
      MARK_USED(descriptor_weights)
      MARK_USED(atom_sizes)
      MARK_USED(atom_coords)
      MARK_USED(state)
      CPABORT("Reconstructed Skala fields require LibTorch and the SIRIUS point XC callback API")
#endif
   END SUBROUTINE
#endif

! **************************************************************************************************
!> \brief Allocate the caller-owned callback state and its C token.
! **************************************************************************************************
   SUBROUTINE create_state(xc_section, particles, cell, state)
      TYPE(section_vals_type), POINTER                   :: xc_section
      TYPE(particle_type), POINTER                       :: particles(:)
      TYPE(cell_type), POINTER                           :: cell
      TYPE(skala_sirius_state_type), POINTER             :: state

      CPASSERT(.NOT. ASSOCIATED(state))
      CPASSERT(next_id < HUGE(next_id))
      ALLOCATE (state)
      next_id = next_id + 1
      state%token%id = next_id
      state%xc_section => xc_section
      CALL skala_sirius_update(state, particles, cell)
      state%next => states
      states => state
   END SUBROUTINE

! **************************************************************************************************
!> \brief Refresh borrowed geometry pointers when the PWDFT subsystem changes.
!> \param state Callback state.
!> \param particles Current particle set.
!> \param cell Current cell.
! **************************************************************************************************
   SUBROUTINE skala_sirius_update(state, particles, cell)
      TYPE(skala_sirius_state_type), POINTER             :: state
      TYPE(particle_type), POINTER                       :: particles(:)
      TYPE(cell_type), POINTER                           :: cell

      IF (.NOT. ASSOCIATED(state)) RETURN
      CPASSERT(ASSOCIATED(particles))
      CPASSERT(ASSOCIATED(cell))
      state%particles => particles
      state%cell => cell
   END SUBROUTINE

! **************************************************************************************************
!> \brief Release callback state after the associated SIRIUS context has been destroyed.
!> \param state Callback state.
! **************************************************************************************************
   SUBROUTINE skala_sirius_release(state)
      TYPE(skala_sirius_state_type), POINTER             :: state

      TYPE(skala_sirius_state_type), POINTER             :: current, previous

      IF (.NOT. ASSOCIATED(state)) RETURN
      current => states
      NULLIFY (previous)
      DO WHILE (ASSOCIATED(current))
         IF (ASSOCIATED(current, state)) EXIT
         previous => current
         current => current%next
      END DO
      CPASSERT(ASSOCIATED(current))
      IF (ASSOCIATED(previous)) THEN
         previous%next => current%next
      ELSE
         states => current%next
      END IF
      DEALLOCATE (state)
      NULLIFY (state)
   END SUBROUTINE

#if defined(__SIRIUS_XC_CALLBACK) && defined(__LIBTORCH)
! **************************************************************************************************
!> \brief SIRIUS C callback. One FFT communicator partitions a complete regular grid.
! **************************************************************************************************
   SUBROUTINE evaluate_skala_fields(user_data, nr, ns, dims, z_offset, fcomm_context, fcomm_fft, &
                                    lattice, weight, rho, gradient, tau, energy, vrho, vgradient, vtau, status) BIND(C)
      TYPE(C_PTR), VALUE                                 :: user_data
      INTEGER(C_INT), VALUE                              :: nr, ns
      INTEGER(C_INT), INTENT(IN)                         :: dims(3)
      INTEGER(C_INT), VALUE                              :: z_offset, fcomm_context, fcomm_fft
      REAL(C_DOUBLE), INTENT(IN)                         :: lattice(3, 3)
      REAL(C_DOUBLE), VALUE                              :: weight
      REAL(C_DOUBLE), INTENT(IN)                         :: rho(nr, ns), gradient(nr, 3, ns), &
                                                            tau(nr, ns)
      REAL(C_DOUBLE), INTENT(OUT)                        :: energy, vrho(nr, ns), &
                                                            vgradient(nr, 3, ns), vtau(nr, ns)
      INTEGER(C_INT), INTENT(OUT)                        :: status

      INTEGER                                            :: i, invalid, nz
      TYPE(callback_token_type), POINTER                 :: token
      TYPE(mp_comm_type), TARGET                         :: context_comm, fft_comm
      TYPE(skala_grid_layout_type)                       :: grid
      TYPE(skala_sirius_state_type), POINTER             :: state

      CALL fft_comm%set_handle(INT(fcomm_fft))
      CALL context_comm%set_handle(INT(fcomm_context))
      energy = 0.0_dp
      vrho = 0.0_dp
      vgradient = 0.0_dp
      vtau = 0.0_dp
      status = 1
      NULLIFY (state)
      IF (C_ASSOCIATED(user_data)) THEN
         CALL C_F_POINTER(user_data, token)
         state => states
         DO WHILE (ASSOCIATED(state))
            IF (state%token%id == token%id) EXIT
            state => state%next
         END DO
      END IF
      invalid = 0
      IF (.NOT. ASSOCIATED(state)) invalid = 1
      IF (nr < 0 .OR. ALL(ns /= [1, 2]) .OR. ANY(dims <= 0) .OR. weight <= 0.0_dp) invalid = 1
      CALL context_comm%max(invalid)
      IF (invalid /= 0) RETURN
      nz = nr/(dims(1)*dims(2))
      IF (nr /= dims(1)*dims(2)*nz .OR. z_offset < 0 .OR. z_offset + nz > dims(3)) invalid = 1
      IF (MAXVAL(ABS(lattice - state%cell%hmat)) > 1.0e-10_dp) invalid = 1
      CALL context_comm%max(invalid)
      IF (invalid /= 0) RETURN

      grid%group => fft_comm
      grid%npts = dims
      grid%bounds(1, :) = 0
      grid%bounds(2, :) = dims - 1
      grid%bounds_local = grid%bounds
      grid%bounds_local(:, 3) = [z_offset, z_offset + nz - 1]
      grid%ngpts_local = nr
      grid%dvol = weight
      DO i = 1, 3
         grid%dh(:, i) = lattice(:, i)/REAL(dims(i), dp)
      END DO
      CALL skala_grid_eval(grid, rho, gradient, tau, state%xc_section, state%particles, state%cell, &
                           energy, vrho, vgradient, vtau)
      status = 0
   END SUBROUTINE
#endif

#if defined(__SIRIUS_XC_POINT_CALLBACK) && defined(__LIBTORCH)
! **************************************************************************************************
!> \brief Return raw discrete covectors for joint fields, without dividing by quadrature weights.
! **************************************************************************************************
   SUBROUTINE evaluate_skala_points(user_data, nr, ns, nb, fcomm, coordinates, energy_weights, &
                                    descriptor_weights, atom_sizes, atom_coords, rho, gradient, tau, &
                                    energy, drho, dgradient, dtau, status) BIND(C)
      TYPE(C_PTR), VALUE                                 :: user_data
      INTEGER(C_INT), VALUE                              :: nr, ns, nb, fcomm
      REAL(C_DOUBLE), INTENT(IN)                         :: coordinates(3, nr), energy_weights(nr), &
                                                            descriptor_weights(nr)
      INTEGER(C_INT), INTENT(IN)                         :: atom_sizes(nb)
      REAL(C_DOUBLE), INTENT(IN)                         :: atom_coords(3, nb), rho(nr, ns), &
                                                            gradient(nr, 3, ns), tau(nr, ns)
      REAL(C_DOUBLE), INTENT(OUT)                        :: energy, drho(nr, ns), &
                                                            dgradient(nr, 3, ns), dtau(nr, ns)
      INTEGER(C_INT), INTENT(OUT)                        :: status

      INTEGER                                            :: invalid
      TYPE(callback_token_type), POINTER                 :: token
      TYPE(mp_comm_type)                                 :: comm
      TYPE(skala_sirius_state_type), POINTER             :: state

      CALL comm%set_handle(INT(fcomm))
      energy = 0.0_dp
      drho = 0.0_dp
      dgradient = 0.0_dp
      dtau = 0.0_dp
      status = 1
      NULLIFY (state)
      IF (C_ASSOCIATED(user_data)) THEN
         CALL C_F_POINTER(user_data, token)
         state => states
         DO WHILE (ASSOCIATED(state))
            IF (state%token%id == token%id) EXIT
            state => state%next
         END DO
      END IF
      invalid = 0
      IF (.NOT. ASSOCIATED(state)) invalid = 1
      IF (nr < 0 .OR. nb < 0 .OR. ALL(ns /= [1, 2])) invalid = 1
      CALL comm%max(invalid)
      IF (invalid /= 0) RETURN
      CALL skala_point_eval(comm, coordinates, energy_weights, descriptor_weights, atom_sizes, atom_coords, &
                            rho, gradient, tau, state%xc_section, energy, drho, dgradient, dtau)
      status = 0
   END SUBROUTINE
#endif

#if defined(__SIRIUS_XC_POINT_GEOMETRY) && defined(__LIBTORCH)
! **************************************************************************************************
!> \brief Return the registry entry without exporting Fortran pointers across the C ABI.
! **************************************************************************************************
   FUNCTION geometry_state(user_data) RESULT(state)
      TYPE(C_PTR), VALUE                                 :: user_data
      TYPE(skala_sirius_state_type), POINTER             :: state

      TYPE(callback_token_type), POINTER                 :: token

      NULLIFY (state)
      IF (.NOT. C_ASSOCIATED(user_data)) RETURN
      CALL C_F_POINTER(user_data, token)
      state => states
      DO WHILE (ASSOCIATED(state))
         IF (state%token%id == token%id) EXIT
         state => state%next
      END DO
   END FUNCTION

! **************************************************************************************************
!> \brief Differentiate the model's independent fields, coordinates and quadrature weights.
! **************************************************************************************************
   SUBROUTINE evaluate_skala_point_geometry(user_data, nr, ns, nb, fcomm, coordinates, energy_weights, &
                                            descriptor_weights, atom_sizes, atom_coords, rho, gradient, tau, &
                                            energy, drho, dgradient, dtau, dpoints, denergy, ddescriptor, &
                                            dcentres, status) BIND(C)
      TYPE(C_PTR), VALUE                                 :: user_data
      INTEGER(C_INT), VALUE                              :: nr, ns, nb, fcomm
      REAL(C_DOUBLE), INTENT(IN)                         :: coordinates(3, nr), energy_weights(nr), &
                                                            descriptor_weights(nr)
      INTEGER(C_INT), INTENT(IN)                         :: atom_sizes(nb)
      REAL(C_DOUBLE), INTENT(IN)                         :: atom_coords(3, nb), rho(nr, ns), &
                                                            gradient(nr, 3, ns), tau(nr, ns)
      REAL(C_DOUBLE), INTENT(OUT) :: energy, drho(nr, ns), dgradient(nr, 3, ns), dtau(nr, ns), &
         dpoints(3, nr), denergy(nr), ddescriptor(nr), dcentres(3, nb)
      INTEGER(C_INT), INTENT(OUT)                        :: status

      INTEGER                                            :: invalid
      TYPE(mp_comm_type)                                 :: comm
      TYPE(skala_sirius_state_type), POINTER             :: state

      CALL comm%set_handle(INT(fcomm))
      status = 1
      state => geometry_state(user_data)
      invalid = 0
      IF (.NOT. ASSOCIATED(state)) invalid = 1
      IF (nr < 0 .OR. nb < 0 .OR. ALL(ns /= [1, 2])) invalid = 1
      CALL comm%max(invalid)
      IF (invalid /= 0) RETURN
      CALL skala_point_eval(comm, coordinates, energy_weights, descriptor_weights, atom_sizes, atom_coords, &
                            rho, gradient, tau, state%xc_section, energy, drho, dgradient, dtau, &
                            dpoints, dcentres, denergy, ddescriptor)
      status = 0
   END SUBROUTINE

! **************************************************************************************************
!> \brief Pull back combined model and field-sampling covectors through the periodic atom grid.
! **************************************************************************************************
   SUBROUTINE pullback_skala_quadrature(user_data, nr, nb, natom, fcomm, atom_sizes, dpoints, &
                                        denergy, ddescriptor, dcentres, datoms, dlattice, status) BIND(C)
      TYPE(C_PTR), VALUE                                 :: user_data
      INTEGER(C_INT), VALUE                              :: nr, nb, natom, fcomm
      INTEGER(C_INT), INTENT(IN)                         :: atom_sizes(nb)
      REAL(C_DOUBLE), INTENT(IN)                         :: dpoints(3, nr), denergy(nr), &
                                                            ddescriptor(nr), dcentres(3, nb)
      REAL(C_DOUBLE), INTENT(OUT)                        :: datoms(3, natom), dlattice(3, 3)
      INTEGER(C_INT), INTENT(OUT)                        :: status

      INTEGER                                            :: a, invalid
      REAL(dp)                                           :: centres(3, natom)
      TYPE(mp_comm_type)                                 :: comm
      TYPE(skala_sirius_state_type), POINTER             :: state

      CALL comm%set_handle(INT(fcomm))
      status = 1
      datoms = 0.0_dp
      dlattice = 0.0_dp
      state => geometry_state(user_data)
      invalid = 0
      IF (.NOT. ASSOCIATED(state)) invalid = 1
      IF (nr < 0 .OR. nb < 0 .OR. natom < 1) invalid = 1
      CALL comm%max(invalid)
      IF (invalid /= 0) RETURN
      IF (SIZE(state%particles) /= natom .OR. state%nradial < 2) invalid = 1
      CALL comm%max(invalid)
      IF (invalid /= 0) RETURN
      DO a = 1, natom
         centres(:, a) = state%particles(a)%r
      END DO
      IF (MAXVAL(ABS(centres - state%grid_centres)) > 1.0e-12_dp .OR. &
          MAXVAL(ABS(state%cell%hmat - state%grid_lattice)) > 1.0e-12_dp) invalid = 1
      CALL comm%max(invalid)
      IF (invalid /= 0) RETURN
      CALL skala_atom_grid_pullback(comm, state%cell%hmat, centres, state%nradial, state%nangular, state%shell, &
                                    state%scale, atom_sizes, dpoints, denergy, ddescriptor, dcentres, datoms, dlattice)
      status = 0
   END SUBROUTINE
#endif
END MODULE skala_sirius
