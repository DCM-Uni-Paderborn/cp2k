!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
!> \brief Evaluate fixed-occupation native orbital trials and physical Skala roots.
MODULE qs_skala_scf
   USE cell_types,                      ONLY: cell_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE input_section_types,             ONLY: section_get_ival
   USE kinds,                           ONLY: dp
   USE kpoint_types,                    ONLY: kpoint_env_type,&
                                              kpoint_type
   USE message_passing,                 ONLY: mp_comm_type,&
                                              mp_para_env_type
   USE particle_types,                  ONLY: particle_type
   USE pw_grid_types,                   ONLY: pw_grid_type
   USE pw_types,                        ONLY: pw_r3d_rs_type
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_ks_methods,                   ONLY: qs_ks_update_qs_env
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_did_change,&
                                              qs_ks_env_type
   USE qs_mo_types,                     ONLY: mo_set_type
   USE qs_native_grid_interpolation,    ONLY: native_grid_batch_adjoint
   USE qs_rho_methods,                  ONLY: qs_rho_update_rho
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_set,&
                                              qs_rho_type
   USE qs_skala_orbitals,               ONLY: qs_skala_orbitals_create,&
                                              qs_skala_orbitals_evaluate,&
                                              qs_skala_orbitals_install_gamma,&
                                              qs_skala_orbitals_install_kpoint,&
                                              qs_skala_orbitals_release,&
                                              qs_skala_orbitals_type
   USE qs_upf_inverse_augmentation,     ONLY: qs_upf_inverse_augmentation_integrate
   USE qs_upf_inverse_background,       ONLY: qs_upf_inverse_background_integrate
   USE qs_upf_inverse_density,          ONLY: qs_upf_inverse_kinetic_operator
   USE qs_upf_inverse_grid,             ONLY: qs_upf_inverse_grid_integrate_dense,&
                                              qs_upf_inverse_matrix_gather
   USE qs_upf_inverse_paw,              ONLY: qs_upf_inverse_paw_operator
   USE qs_upf_inverse_types,            ONLY: qs_upf_inverse_augmentation_type,&
                                              qs_upf_inverse_scf_type,&
                                              qs_upf_inverse_type
   USE skala_density_constraints,       ONLY: &
        skala_background_domain, skala_density_constraints_type, skala_density_descriptor, &
        skala_density_crossings, skala_density_root, skala_density_samples_begin, skala_density_samples_end, &
        skala_density_samples_get, skala_fft_domain, skala_kinetic_descriptor
   USE skala_orbital_constraints,       ONLY: skala_orbital_chart_derivative,&
                                              skala_orbital_chart_size
   USE upf_inverse_quadrature,          ONLY: upf_inverse_semilocal_pair
   USE xc_input_constants,              ONLY: xc_rho_no_smooth
   USE xc_util,                         ONLY: xc_pw_smooth
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: qs_skala_scf_type, qs_skala_scf_create, qs_skala_scf_release, qs_skala_scf_size
   PUBLIC :: qs_skala_scf_evaluate, qs_skala_scf_value, qs_skala_scf_root
   PUBLIC :: qs_skala_scf_point_derivative
   PUBLIC :: qs_skala_scf_crossings

   TYPE qs_skala_scf_type
      PRIVATE
      TYPE(qs_environment_type), POINTER :: qs_env => NULL()
      TYPE(qs_ks_env_type), POINTER :: ks_env => NULL()
      TYPE(qs_rho_type), POINTER :: rho => NULL()
      TYPE(qs_energy_type), POINTER :: energy => NULL()
      TYPE(qs_upf_inverse_scf_type), POINTER :: inverse => NULL()
      TYPE(skala_density_constraints_type), POINTER :: constraints => NULL()
      TYPE(kpoint_type), POINTER :: kpoints => NULL()
      TYPE(mp_para_env_type), POINTER :: para_env => NULL()
      TYPE(qs_skala_orbitals_type), ALLOCATABLE :: orbitals(:, :)
      INTEGER, ALLOCATABLE :: dimensions(:, :), offsets(:, :)
      REAL(dp), ALLOCATABLE :: coordinate(:), positions(:, :)
      REAL(dp) :: hmat(3, 3) = 0.0_dp
      LOGICAL :: bloch = .FALSE., evaluated = .FALSE.
   END TYPE qs_skala_scf_type

   TYPE physical_root_context
      TYPE(qs_skala_scf_type), POINTER :: state => NULL()
      REAL(dp), ALLOCATABLE :: first(:), direction(:)
      INTEGER :: domain, point(3), spin, descriptor
   END TYPE physical_root_context

   TYPE physical_domain_samples
      INTEGER, ALLOCATABLE :: points(:, :)
      REAL(dp), ALLOCATABLE :: density(:, :), kinetic(:, :)
      LOGICAL :: available = .FALSE.
   END TYPE physical_domain_samples

CONTAINS

! **************************************************************************************************
!> \brief Borrow the native QS state and complete each occupied block in its original charge metric.
!>        Collective over QS and each owning k group. Zero-weight and empty-spin blocks have no
!>        coordinates. Geometry, occupations and k weights stay fixed throughout this state.
!> \param state Owned orbital charts, borrowing a QS environment that must outlive this state.
!> \param qs_env Native USPP or PAW environment with initialized MOs and normalization operators.
!> \param valid Whether every active block has a positive complete charge metric.
! **************************************************************************************************
   SUBROUTINE qs_skala_scf_create(state, qs_env, valid)
      TYPE(qs_skala_scf_type), INTENT(INOUT)            :: state
      TYPE(qs_environment_type), POINTER               :: qs_env
      LOGICAL, INTENT(OUT)                              :: valid

      COMPLEX(dp), ALLOCATABLE                         :: metric(:, :)
      INTEGER                                          :: i, ik, local, nk, nl, ns, offset, spin, failed
      LOGICAL                                          :: block_valid, complex_orbitals, owner
      REAL(dp)                                         :: weight
      TYPE(cell_type), POINTER                         :: cell
      TYPE(dbcsr_p_type), POINTER                      :: metrics(:)
      TYPE(kpoint_env_type), POINTER                   :: kp
      TYPE(mo_set_type), POINTER                       :: mos(:), native
      TYPE(particle_type), POINTER                     :: particles(:)

      CALL qs_skala_scf_release(state)
      CPASSERT(ASSOCIATED(qs_env))
      state%qs_env => qs_env
      CALL get_qs_env(qs_env, ks_env=state%ks_env, rho=state%rho, energy=state%energy, &
                      kpoints=state%kpoints, mos=mos, para_env=state%para_env, cell=cell, particle_set=particles)
      CALL get_ks_env(state%ks_env, upf_inverse=state%inverse, skala_constraints=state%constraints)
      valid = ASSOCIATED(state%inverse) .AND. ASSOCIATED(state%constraints)
      IF (valid) valid = state%inverse%enabled
      IF (.NOT. valid) THEN
         CALL qs_skala_scf_release(state)
         RETURN
      END IF
      state%bloch = ASSOCIATED(state%inverse%bloch%kpoints)
      nk = 1
      nl = 1
      ns = SIZE(mos)
      owner = state%para_env%mepos == 0
      complex_orbitals = .FALSE.
      IF (state%bloch) THEN
         nk = SIZE(state%kpoints%wkp)
         nl = SIZE(state%kpoints%kp_env)
         ns = SIZE(state%kpoints%kp_env(1)%kpoint_env%mos)
         owner = state%kpoints%para_env_kp%mepos == 0
         complex_orbitals = .NOT. state%kpoints%use_real_wfn
      ELSE
         CALL get_ks_env(state%ks_env, matrix_metric=metrics)
      END IF
      ALLOCATE (state%orbitals(nl, ns), state%dimensions(nk, ns), state%offsets(nk, ns), &
                state%positions(3, SIZE(particles)))
      state%hmat(:, :) = cell%hmat
      DO i = 1, SIZE(particles)
         state%positions(:, i) = particles(i)%r
      END DO
      state%dimensions = 0
      failed = 0
      DO spin = 1, ns
         DO local = 1, nl
            ik = 1
            weight = 1.0_dp
            native => mos(spin)
            IF (state%bloch) THEN
               kp => state%kpoints%kp_env(local)%kpoint_env
               ik = kp%nkpoint
               native => kp%mos(spin)
               weight = state%kpoints%wkp(ik)
            END IF
            CPASSERT(ieee_is_finite(weight) .AND. weight >= 0.0_dp)
            IF (weight == 0.0_dp .OR. native%nmo == 0) CYCLE
            IF (ALL(native%occupation_numbers(:native%nmo) == 0.0_dp)) CYCLE
            IF (state%bloch) THEN
               IF (ALLOCATED(metric)) DEALLOCATE (metric)
               ALLOCATE (metric(native%nao, native%nao))
               metric(:, :) = state%inverse%metric_k(:, :, ik)
            ELSE
               CALL qs_upf_inverse_matrix_gather(state%inverse%frame%ao, metrics(1)%matrix, metric)
            END IF
            CALL qs_skala_orbitals_create(state%orbitals(local, spin), native, metric, &
                                           complex_orbitals, weight, block_valid)
            IF (.NOT. block_valid) failed = 1
            IF (block_valid .AND. owner) THEN
               state%dimensions(ik, spin) = skala_orbital_chart_size(state%orbitals(local, spin)%chart)
            END IF
         END DO
      END DO
      CALL state%para_env%max(failed)
      valid = failed == 0
      IF (.NOT. valid) THEN
         CALL qs_skala_scf_release(state)
         RETURN
      END IF
      CALL state%para_env%sum(state%dimensions)
      offset = 0
      DO spin = 1, ns
         DO ik = 1, nk
            state%offsets(ik, spin) = offset
            offset = offset + state%dimensions(ik, spin)
         END DO
      END DO
      ALLOCATE (state%coordinate(offset))
      state%coordinate = 0.0_dp
   END SUBROUTINE qs_skala_scf_create

! **************************************************************************************************
!> \brief Release owned charts and buffers without releasing or resetting the borrowed native QS state.
!> \param state ...
! **************************************************************************************************
   SUBROUTINE qs_skala_scf_release(state)
      TYPE(qs_skala_scf_type), INTENT(INOUT)            :: state

      INTEGER                                          :: local, spin

      IF (ALLOCATED(state%orbitals)) THEN
         DO spin = 1, SIZE(state%orbitals, 2)
            DO local = 1, SIZE(state%orbitals, 1)
               CALL qs_skala_orbitals_release(state%orbitals(local, spin))
            END DO
         END DO
         DEALLOCATE (state%orbitals, state%dimensions, state%offsets, state%positions)
      END IF
      IF (ALLOCATED(state%coordinate)) DEALLOCATE (state%coordinate)
      NULLIFY (state%qs_env, state%ks_env, state%rho, state%energy, state%inverse, state%constraints, &
               state%kpoints, state%para_env)
      state%hmat = 0.0_dp
      state%bloch = .FALSE.
      state%evaluated = .FALSE.
   END SUBROUTINE qs_skala_scf_release

! **************************************************************************************************
!> \brief Number of global independent fixed-occupation orbital coordinates.
!> \param state ...
!> \return ...
! **************************************************************************************************
   FUNCTION qs_skala_scf_size(state) RESULT(n)
      TYPE(qs_skala_scf_type), INTENT(IN)               :: state
      INTEGER                                          :: n

      CPASSERT(ALLOCATED(state%coordinate))
      n = SIZE(state%coordinate)
   END FUNCTION qs_skala_scf_size

! **************************************************************************************************
!> \brief Install all owning orbital blocks, rebuild rho/tau and evaluate the complete native energy.
!>        Samples are the actual local model inputs of this trial, including atomic/background maps.
!>        Gradients are globally ordered and include each k weight once. Reevaluate the accepted
!>        coordinate after rejecting a trial; native MOs, P and fields then describe the same state.
!> \param state ...
!> \param coordinate Replicated global coordinates, in spin then global k-point order.
!> \param energy Complete native energy at fixed occupations.
!> \param gradient Optional derivative in the supplied Cayley coordinates.
!> \param physical_gradient Optional derivative in the orthonormal physical orbital tangent.
! **************************************************************************************************
   SUBROUTINE qs_skala_scf_evaluate(state, coordinate, energy, gradient, physical_gradient)
      TYPE(qs_skala_scf_type), INTENT(INOUT)            :: state
      REAL(dp), INTENT(IN)                              :: coordinate(:)
      REAL(dp), INTENT(OUT)                             :: energy
      REAL(dp), INTENT(OUT), OPTIONAL                   :: gradient(:), physical_gradient(:)

      COMPLEX(dp), ALLOCATABLE                         :: hamiltonian(:, :)
      INTEGER                                          :: first, ik, last, local, spin
      LOGICAL                                          :: owner, need_gradient, valid
      REAL(dp), ALLOCATABLE                            :: block_gradient(:)
      TYPE(cell_type), POINTER                         :: cell
      TYPE(dbcsr_p_type), POINTER                      :: densities(:, :), matrices(:, :)
      TYPE(particle_type), POINTER                     :: particles(:)

      CPASSERT(ASSOCIATED(state%qs_env))
      CPASSERT(SIZE(coordinate) == qs_skala_scf_size(state))
      CPASSERT(ALL(ieee_is_finite(coordinate)))
      CALL get_qs_env(state%qs_env, cell=cell, particle_set=particles)
      CPASSERT(ALL(cell%hmat == state%hmat))
      CPASSERT(SIZE(particles) == SIZE(state%positions, 2))
      DO local = 1, SIZE(particles)
         CPASSERT(ALL(particles(local)%r == state%positions(:, local)))
      END DO
      CALL qs_rho_get(state%rho, rho_ao_kp=densities)
      DO spin = 1, SIZE(state%orbitals, 2)
         DO local = 1, SIZE(state%orbitals, 1)
            IF (.NOT. ASSOCIATED(state%orbitals(local, spin)%native)) CYCLE
            valid = ALL(state%orbitals(local, spin)%native%occupation_numbers(: &
                        state%orbitals(local, spin)%native%nmo) == &
                        state%orbitals(local, spin)%occupations(:state%orbitals(local, spin)%native%nmo))
            CPASSERT(valid)
            ik = 1
            IF (state%bloch) THEN
               ik = state%kpoints%kp_env(local)%kpoint_env%nkpoint
               CPASSERT(state%kpoints%wkp(ik) == state%orbitals(local, spin)%weight)
            END IF
            first = state%offsets(ik, spin) + 1
            last = first + state%dimensions(ik, spin) - 1
            CALL qs_skala_orbitals_evaluate(state%orbitals(local, spin), coordinate(first:last))
            IF (state%bloch) THEN
               CALL qs_skala_orbitals_install_kpoint(state%orbitals(local, spin), state%kpoints, local, spin)
            ELSE
               CALL qs_skala_orbitals_install_gamma(state%orbitals(local, spin), &
                                                      state%inverse%frame%ao, densities(spin, 1)%matrix)
            END IF
         END DO
      END DO
      CALL qs_rho_set(state%rho, rho_r_valid=.FALSE., rho_g_valid=.FALSE., tau_r_valid=.FALSE., tau_g_valid=.FALSE.)
      CALL qs_rho_update_rho(state%rho, state%qs_env)
      CALL skala_density_samples_begin(state%constraints)
      CALL qs_ks_did_change(state%ks_env, rho_changed=.TRUE.)
      need_gradient = PRESENT(gradient) .OR. PRESENT(physical_gradient)
      CALL qs_ks_update_qs_env(state%qs_env, calculate_forces=.FALSE., &
                               just_energy=.NOT. need_gradient, print_active=.FALSE.)
      CALL skala_density_samples_end(state%constraints)
      energy = state%energy%total
      IF (PRESENT(gradient)) THEN
         CPASSERT(SIZE(gradient) == SIZE(coordinate))
         gradient = 0.0_dp
      END IF
      IF (PRESENT(physical_gradient)) THEN
         CPASSERT(SIZE(physical_gradient) == SIZE(coordinate))
         physical_gradient = 0.0_dp
      END IF
      IF (need_gradient) THEN
         owner = state%para_env%mepos == 0
         IF (state%bloch) THEN
            owner = state%kpoints%para_env_kp%mepos == 0
         ELSE
            CALL get_ks_env(state%ks_env, matrix_ks_kp=matrices)
         END IF
         DO spin = 1, SIZE(state%orbitals, 2)
            DO local = 1, SIZE(state%orbitals, 1)
               IF (.NOT. ASSOCIATED(state%orbitals(local, spin)%native)) CYCLE
               ik = 1
               IF (state%bloch) ik = state%kpoints%kp_env(local)%kpoint_env%nkpoint
               first = state%offsets(ik, spin) + 1
               last = first + state%dimensions(ik, spin) - 1
               IF (state%bloch) THEN
                  IF (ALLOCATED(hamiltonian)) DEALLOCATE (hamiltonian)
                  ALLOCATE (hamiltonian(SIZE(state%inverse%ks_k, 1), SIZE(state%inverse%ks_k, 2)))
                  hamiltonian(:, :) = state%inverse%ks_k(:, :, ik, spin)
               ELSE
                  CALL qs_upf_inverse_matrix_gather(state%inverse%frame%ao, matrices(spin, 1)%matrix, hamiltonian)
               END IF
               IF (.NOT. owner) CYCLE
               ALLOCATE (block_gradient(last - first + 1))
               hamiltonian = state%orbitals(local, spin)%weight*hamiltonian
               IF (PRESENT(gradient)) THEN
                  CALL skala_orbital_chart_derivative(state%orbitals(local, spin)%chart, hamiltonian, block_gradient)
                  gradient(first:last) = block_gradient
               END IF
               IF (PRESENT(physical_gradient)) THEN
                  CALL skala_orbital_chart_derivative(state%orbitals(local, spin)%chart, hamiltonian, &
                                                       block_gradient, physical_tangent=.TRUE.)
                  physical_gradient(first:last) = block_gradient
               END IF
               DEALLOCATE (block_gradient)
            END DO
         END DO
         IF (PRESENT(gradient)) CALL state%para_env%sum(gradient)
         IF (PRESENT(physical_gradient)) CALL state%para_env%sum(physical_gradient)
      END IF
      state%coordinate(:) = coordinate
      state%evaluated = .TRUE.
   END SUBROUTINE qs_skala_scf_evaluate

! **************************************************************************************************
!> \brief Read one actual physical model input, requiring exactly one distributed owner.
!> \param state ...
!> \param domain FFT, exact atomic or interpolated atomic background.
!> \param point Global FFT index or (global atom, row within that atom, 1).
!> \param spin Zero for the paired RKS input, one or two for UKS.
!> \param descriptor Density or kinetic density.
!> \param value Physical model input, including native core and density processing.
! **************************************************************************************************
   SUBROUTINE qs_skala_scf_value(state, domain, point, spin, descriptor, value)
      TYPE(qs_skala_scf_type), INTENT(IN)               :: state
      INTEGER, INTENT(IN)                               :: domain, point(3), spin, descriptor
      REAL(dp), INTENT(OUT)                             :: value

      INTEGER                                          :: column, i, owners
      INTEGER, ALLOCATABLE                             :: points(:, :)
      LOGICAL                                          :: available, valid
      REAL(dp), ALLOCATABLE                            :: density(:, :), kinetic(:, :)

      CPASSERT(state%evaluated)
      valid = domain >= skala_fft_domain .AND. domain <= skala_background_domain
      CPASSERT(valid)
      valid = descriptor == skala_density_descriptor .OR. descriptor == skala_kinetic_descriptor
      CPASSERT(valid)
      CALL skala_density_samples_get(state%constraints, domain, points, density, kinetic, available)
      owners = 0
      value = 0.0_dp
      IF (available) THEN
         column = spin
         IF (SIZE(density, 2) == 1) THEN
            CPASSERT(spin == 0)
            column = 1
         END IF
         CPASSERT(column >= 1 .AND. column <= SIZE(density, 2))
         DO i = 1, SIZE(points, 2)
            IF (.NOT. ALL(points(:, i) == point)) CYCLE
            owners = owners + 1
            IF (descriptor == skala_density_descriptor) THEN
               value = density(i, column)
            ELSE
               value = kinetic(i, column)
            END IF
         END DO
      END IF
      CALL state%para_env%sum(owners)
      CPASSERT(owners == 1)
      CALL state%para_env%sum(value)
   END SUBROUTINE qs_skala_scf_value

! **************************************************************************************************
!> \brief Differentiate one actual physical rho/tau input without changing the native KS operator.
!>        Core terms are constant at fixed geometry. The paired-spin factor and each k weight are
!>        applied once. Atomic background inputs use the transpose of their actual interpolation,
!>        without the negative sign belonging to the background energy contribution.
!> \param state Evaluated native orbital state, with current physical model samples.
!> \param domain FFT, exact atomic or interpolated atomic background.
!> \param point Global FFT index or (global atom, row within that atom, 1).
!> \param spin Zero for the paired RKS input, one or two for UKS.
!> \param descriptor Density or kinetic density.
!> \param derivative Global derivative in the current Cayley coordinates.
!> \param physical_derivative Optional covector in the orthonormal physical orbital tangent.
! **************************************************************************************************
   SUBROUTINE qs_skala_scf_point_derivative(state, domain, point, spin, descriptor, derivative, physical_derivative)
      TYPE(qs_skala_scf_type), INTENT(INOUT)             :: state
      INTEGER, INTENT(IN)                                :: domain, point(3), spin, descriptor
      REAL(dp), INTENT(OUT)                              :: derivative(:)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: physical_derivative(:)

      COMPLEX(dp), ALLOCATABLE                           :: matrix(:, :), operators(:, :, :)
      INTEGER                                            :: first, ik, last, local, ns, selected_spin
      LOGICAL                                            :: owner
      REAL(dp)                                           :: scale, value

      CALL qs_skala_scf_value(state, domain, point, spin, descriptor, value)
      CPASSERT(SIZE(derivative) == qs_skala_scf_size(state))
      derivative = 0.0_dp
      IF (PRESENT(physical_derivative)) THEN
         CPASSERT(SIZE(physical_derivative) == SIZE(derivative))
         physical_derivative = 0.0_dp
      END IF
      ns = SIZE(state%orbitals, 2)
      selected_spin = spin
      scale = 1.0_dp
      IF (ns == 1) THEN
         selected_spin = 1
         scale = 0.5_dp
      END IF
      IF (domain == skala_fft_domain .OR. domain == skala_background_domain) THEN
         CALL grid_point_operators(state, domain, point, descriptor, operators)
      ELSE
         CALL atom_point_operators(state, point, descriptor, operators)
      END IF
      owner = state%para_env%mepos == 0
      IF (state%bloch) owner = state%kpoints%para_env_kp%mepos == 0
      ALLOCATE (matrix(SIZE(operators, 1), SIZE(operators, 2)))
      DO local = 1, SIZE(state%orbitals, 1)
         IF (.NOT. ASSOCIATED(state%orbitals(local, selected_spin)%native)) CYCLE
         ik = 1
         IF (state%bloch) ik = state%kpoints%kp_env(local)%kpoint_env%nkpoint
         first = state%offsets(ik, selected_spin) + 1
         last = first + state%dimensions(ik, selected_spin) - 1
         IF (.NOT. owner) CYCLE
         matrix(:, :) = scale*state%orbitals(local, selected_spin)%weight*operators(:, :, ik)
         CALL skala_orbital_chart_derivative(state%orbitals(local, selected_spin)%chart, &
                                             matrix, derivative(first:last))
         IF (PRESENT(physical_derivative)) THEN
            CALL skala_orbital_chart_derivative(state%orbitals(local, selected_spin)%chart, &
                                                matrix, physical_derivative(first:last), physical_tangent=.TRUE.)
         END IF
      END DO
      CALL state%para_env%sum(derivative)
      IF (PRESENT(physical_derivative)) CALL state%para_env%sum(physical_derivative)
   END SUBROUTINE qs_skala_scf_point_derivative

   SUBROUTINE grid_point_operators(state, domain, point, descriptor, operators)
      TYPE(qs_skala_scf_type), INTENT(INOUT)             :: state
      INTEGER, INTENT(IN)                                :: domain, point(3), descriptor
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: operators(:, :, :)

      COMPLEX(dp), ALLOCATABLE                           :: matrix(:, :)
      INTEGER                                            :: axis, bounds(2, 3), i, ik, nao, nk, row, &
                                                            smooth
      REAL(dp), ALLOCATABLE                              :: fields(:, :, :, :), values(:, :)
      TYPE(cell_type), POINTER                           :: cell
      TYPE(pw_grid_type), POINTER                        :: grid
      TYPE(pw_r3d_rs_type)                               :: density, filtered, gradient(3), kinetic
      TYPE(pw_r3d_rs_type), POINTER                      :: rho_r(:)
      TYPE(qs_upf_inverse_augmentation_type), POINTER    :: augmentation

      CALL qs_rho_get(state%rho, rho_r=rho_r)
      grid => rho_r(1)%pw_grid
      IF (domain == skala_background_domain) grid => state%inverse%xc_grid
      CPASSERT(ASSOCIATED(grid))
      CALL density%create(grid)
      CALL kinetic%create(grid)
      density%array = 0.0_dp
      kinetic%array = 0.0_dp
      DO axis = 1, 3
         CALL gradient(axis)%create(grid)
         gradient(axis)%array = 0.0_dp
      END DO
      IF (domain == skala_fft_domain) THEN
         IF (ALL(point >= grid%bounds_local(1, :)) .AND. ALL(point <= grid%bounds_local(2, :))) THEN
            IF (descriptor == skala_density_descriptor) THEN
               density%array(point(1), point(2), point(3)) = 1.0_dp/grid%dvol
            ELSE
               kinetic%array(point(1), point(2), point(3)) = 1.0_dp/grid%dvol
            END IF
         END IF
         smooth = section_get_ival(state%qs_env%input, "DFT%XC%XC_GRID%XC_SMOOTH_RHO")
         IF (descriptor == skala_density_descriptor .AND. smooth /= xc_rho_no_smooth) THEN
            CALL filtered%create(grid)
            ! The native smoothing stencils are symmetric on the auxiliary PW grid.
            CALL xc_pw_smooth(density, filtered, smooth)
            density%array(:, :, :) = filtered%array
            CALL filtered%release()
         END IF
      ELSE
         ALLOCATE (fields(grid%npts(1), grid%npts(2), grid%npts(3), 1), &
                   values(SIZE(state%inverse%atom_grid%atoms), 1))
         fields = 0.0_dp
         values = 0.0_dp
         row = 0
         DO i = 1, SIZE(values, 1)
            IF (state%inverse%atom_grid%atoms(i) /= point(1)) CYCLE
            row = row + 1
            IF (row == point(2)) values(i, 1) = 1.0_dp
         END DO
         CALL get_qs_env(state%qs_env, cell=cell)
         CALL native_grid_batch_adjoint(values, grid, cell, state%inverse%atom_grid%points, fields, &
                                        wrap_auxiliary_cell=.TRUE., continuous=.TRUE.)
         CALL state%para_env%sum(fields)
         bounds = grid%bounds_local - SPREAD(grid%bounds(1, :) - 1, 1, 2)
         IF (descriptor == skala_density_descriptor) THEN
            density%array(:, :, :) = fields(bounds(1, 1):bounds(2, 1), bounds(1, 2):bounds(2, 2), &
                                            bounds(1, 3):bounds(2, 3), 1)/grid%dvol
         ELSE
            kinetic%array(:, :, :) = fields(bounds(1, 1):bounds(2, 1), bounds(1, 2):bounds(2, 2), &
                                            bounds(1, 3):bounds(2, 3), 1)/grid%dvol
         END IF
      END IF
      nk = SIZE(state%dimensions, 1)
      IF (state%bloch) THEN
         nao = SIZE(state%inverse%bloch%frames(1)%ao%coefficients, 2)
      ELSE
         nao = SIZE(state%inverse%frame%ao%coefficients, 2)
      END IF
      ALLOCATE (operators(nao, nao, nk))
      IF (state%inverse%gaussian_background) THEN
         IF (descriptor == skala_density_descriptor) THEN
            CALL qs_upf_inverse_background_integrate(state%inverse%coulomb%background, state%qs_env, density, operators)
         ELSE
            CALL qs_upf_inverse_background_integrate(state%inverse%coulomb%background, state%qs_env, kinetic, operators, &
                                                     kinetic=.TRUE.)
         END IF
      ELSE
         augmentation => state%inverse%augmentation
         IF (state%bloch) augmentation => state%inverse%bloch%augmentation
         DO ik = 1, nk
            IF (state%bloch) THEN
               CALL qs_upf_inverse_grid_integrate_dense(state%inverse%bloch%frames(ik), density, gradient, kinetic, &
                                                        matrix, augmentation=augmentation, fourier_charge=.TRUE., &
                                                        kinetic=state%inverse%kinetic_cache)
            ELSE
               CALL qs_upf_inverse_grid_integrate_dense(state%inverse%frame, density, gradient, kinetic, matrix, &
                                                        augmentation=augmentation, fourier_charge=.TRUE., &
                                                        kinetic=state%inverse%kinetic_cache)
            END IF
            operators(:, :, ik) = matrix
         END DO
      END IF
      CALL density%release()
      CALL kinetic%release()
      DO axis = 1, 3
         CALL gradient(axis)%release()
      END DO
   END SUBROUTINE grid_point_operators

   SUBROUTINE atom_point_operators(state, point, descriptor, operators)
      TYPE(qs_skala_scf_type), INTENT(INOUT)             :: state
      INTEGER, INTENT(IN)                                :: point(3), descriptor
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: operators(:, :, :)

      COMPLEX(dp), ALLOCATABLE                           :: correction(:, :), q(:, :)
      INTEGER                                            :: i, ik, nao, nk, np, row
      REAL(dp), ALLOCATABLE                              :: density(:), gradient(:, :), kinetic(:), &
                                                            weights(:)
      TYPE(mp_comm_type)                                 :: comm
      TYPE(qs_upf_inverse_augmentation_type), POINTER    :: augmentation
      TYPE(qs_upf_inverse_type), POINTER                 :: frame

      CPASSERT(point(3) == 1)
      CPASSERT(ALLOCATED(state%inverse%xc%soft))
      np = SIZE(state%inverse%atom_grid%atoms)
      nk = SIZE(state%dimensions, 1)
      nao = SIZE(state%inverse%xc%soft, 2)
      ALLOCATE (density(np), kinetic(np), gradient(np, 3), weights(np), operators(nao, nao, nk), correction(nao, nao))
      density = 0.0_dp
      kinetic = 0.0_dp
      gradient = 0.0_dp
      weights = 1.0_dp
      row = 0
      DO i = 1, np
         IF (state%inverse%atom_grid%atoms(i) /= point(1)) CYCLE
         row = row + 1
         IF (row /= point(2)) CYCLE
         IF (descriptor == skala_density_descriptor) THEN
            density(i) = 1.0_dp
         ELSE
            kinetic(i) = 1.0_dp
         END IF
      END DO
      CALL comm%set_handle(state%para_env%get_handle())
      augmentation => state%inverse%augmentation
      IF (state%bloch) augmentation => state%inverse%bloch%augmentation
      DO ik = 1, nk
         frame => state%inverse%frame%ao
         IF (state%bloch) frame => state%inverse%bloch%frames(ik)%ao
         CALL upf_inverse_semilocal_pair(state%inverse%xc%soft(:, :, ik), state%inverse%xc%soft(:, :, ik), &
                                         state%inverse%xc%gradient(:, :, :, ik), state%inverse%xc%gradient(:, :, :, ik), &
                                         weights, density, gradient, kinetic, operators(:, :, ik), comm)
         IF (descriptor == skala_density_descriptor) THEN
            ALLOCATE (q(SIZE(frame%coefficients, 1), SIZE(frame%coefficients, 1)))
            CALL qs_upf_inverse_augmentation_integrate(frame, augmentation, density, gradient, q, &
                                                       state%inverse%atom_grid%points, weights, comm)
            operators(:, :, ik) = operators(:, :, ik) + &
                                  MATMUL(CONJG(TRANSPOSE(frame%coefficients)), MATMUL(q, frame%coefficients))
            DEALLOCATE (q)
            IF (ALLOCATED(state%inverse%paw_reconstruction)) THEN
               CALL qs_upf_inverse_paw_operator(frame, state%inverse%paw_reconstruction, augmentation, &
                                                state%inverse%atom_grid%points, weights, density, gradient, correction, comm)
               operators(:, :, ik) = operators(:, :, ik) + correction
            END IF
         ELSE
            CPASSERT(ALLOCATED(state%inverse%kinetic_cache))
            CALL qs_upf_inverse_kinetic_operator(frame, state%inverse%kinetic_cache, state%inverse%atom_grid%points, &
                                                 weights, kinetic, correction, comm)
            operators(:, :, ik) = operators(:, :, ik) + correction
            IF (ALLOCATED(state%inverse%paw_reconstruction)) THEN
               CALL qs_upf_inverse_kinetic_operator(frame, state%inverse%paw_reconstruction, state%inverse%atom_grid%points, &
                                                    weights, kinetic, correction, comm)
               operators(:, :, ik) = operators(:, :, ik) + correction
            END IF
         END IF
      END DO
   END SUBROUTINE atom_point_operators

! **************************************************************************************************
!> \brief Find a physical root along a normalized orbital curve by complete native reevaluations.
!>        Density matrices and physical fields are never interpolated. On failure restore the
!>        previously evaluated coordinate. On success all native caches describe the returned root.
!> \param state ...
!> \param first Initial global orbital coordinate.
!> \param last Final global orbital coordinate.
!> \param domain ...
!> \param point ...
!> \param spin ...
!> \param descriptor ...
!> \param tolerance Absolute physical root tolerance in atomic units.
!> \param coordinate Root coordinate on success, previous evaluated coordinate on failure.
!> \param value Best physical residual found by the root search.
!> \param converged Whether the actual physical input meets the requested tolerance.
! **************************************************************************************************
   SUBROUTINE qs_skala_scf_root(state, first, last, domain, point, spin, descriptor, tolerance, &
                                 coordinate, value, converged)
      TYPE(qs_skala_scf_type), INTENT(INOUT), TARGET    :: state
      REAL(dp), INTENT(IN)                              :: first(:), last(:), tolerance
      INTEGER, INTENT(IN)                               :: domain, point(3), spin, descriptor
      REAL(dp), INTENT(OUT)                             :: coordinate(:), value
      LOGICAL, INTENT(OUT)                              :: converged

      REAL(dp)                                         :: energy, step
      TYPE(physical_root_context)                      :: context

      CPASSERT(state%evaluated)
      CPASSERT(SIZE(first) == qs_skala_scf_size(state) .AND. SIZE(last) == SIZE(first))
      CPASSERT(SIZE(coordinate) == SIZE(first))
      ALLOCATE (context%first(SIZE(first)), context%direction(SIZE(first)))
      context%state => state
      context%first(:) = first
      context%direction(:) = last - first
      context%domain = domain
      context%point = point
      context%spin = spin
      context%descriptor = descriptor
      coordinate(:) = state%coordinate
      CALL skala_density_root(evaluate_physical_input, context, 0.0_dp, 1.0_dp, tolerance, step, value, converged)
      IF (converged) coordinate(:) = first + step*context%direction
      CALL qs_skala_scf_evaluate(state, coordinate, energy)
      IF (converged) CALL qs_skala_scf_value(state, domain, point, spin, descriptor, value)
   END SUBROUTINE qs_skala_scf_root

! **************************************************************************************************
!> \brief Collect jointly bracketed physical rho/tau zeros from two actual native orbital states.
!>        Samples remain on their native owners; only candidate keys and endpoint values are
!>        gathered. Identically zero fields supply no candidate. This endpoint test does not
!>        locate roots or exclude an even number of crossings along the nonlinear orbital curve.
!>        Restore the first state, including its complete native KS operator, before returning.
!> \param state ...
!> \param first Initial global orbital coordinate.
!> \param last Final global orbital coordinate.
!> \param points Global FFT indices or (global atom, row within that atom, 1).
!> \param domains FFT, exact atomic or interpolated atomic background for each candidate.
!> \param spins Zero for paired RKS, one or two for UKS.
!> \param descriptors Density or kinetic density for each candidate.
!> \param endpoint_values Actual first and last physical inputs, one candidate per column.
! **************************************************************************************************
   SUBROUTINE qs_skala_scf_crossings(state, first, last, points, domains, spins, descriptors, endpoint_values)
      TYPE(qs_skala_scf_type), INTENT(INOUT)             :: state
      REAL(dp), INTENT(IN)                               :: first(:), last(:)
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: points(:, :), domains(:), spins(:), &
                                                            descriptors(:)
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: endpoint_values(:, :)

      INTEGER                                            :: descriptor, domain, i, local_count, &
                                                            offset, row, spin, total
      INTEGER, ALLOCATABLE                               :: counts(:), keys(:, :)
      LOGICAL                                            :: valid
      LOGICAL, ALLOCATABLE                               :: mask(:, :)
      REAL(dp)                                           :: energy
      REAL(dp), ALLOCATABLE                              :: gradient(:), values(:, :)
      TYPE(physical_domain_samples) :: final(skala_fft_domain:skala_background_domain), &
         initial(skala_fft_domain:skala_background_domain)

      valid = SIZE(first) == qs_skala_scf_size(state) .AND. SIZE(last) == SIZE(first)
      CPASSERT(valid)
      CALL qs_skala_scf_evaluate(state, first, energy)
      DO domain = skala_fft_domain, skala_background_domain
         CALL skala_density_samples_get(state%constraints, domain, initial(domain)%points, &
                                        initial(domain)%density, initial(domain)%kinetic, initial(domain)%available)
      END DO
      CALL qs_skala_scf_evaluate(state, last, energy)
      local_count = 0
      DO domain = skala_fft_domain, skala_background_domain
         CALL skala_density_samples_get(state%constraints, domain, final(domain)%points, &
                                        final(domain)%density, final(domain)%kinetic, final(domain)%available)
         valid = initial(domain)%available .EQV. final(domain)%available
         CPASSERT(valid)
         IF (.NOT. initial(domain)%available) CYCLE
         valid = ALL(SHAPE(initial(domain)%points) == SHAPE(final(domain)%points))
         CPASSERT(valid)
         valid = ALL(initial(domain)%points == final(domain)%points)
         CPASSERT(valid)
         valid = ALL(SHAPE(initial(domain)%density) == SHAPE(final(domain)%density))
         CPASSERT(valid)
         ALLOCATE (mask(SIZE(initial(domain)%density, 1), SIZE(initial(domain)%density, 2)))
         CALL skala_density_crossings(initial(domain)%density, final(domain)%density, mask)
         local_count = local_count + COUNT(mask)
         CALL skala_density_crossings(initial(domain)%kinetic, final(domain)%kinetic, mask)
         local_count = local_count + COUNT(mask)
         DEALLOCATE (mask)
      END DO
      ALLOCATE (keys(6, local_count), values(2, local_count))
      row = 0
      DO domain = skala_fft_domain, skala_background_domain
         IF (.NOT. initial(domain)%available) CYCLE
         ALLOCATE (mask(SIZE(initial(domain)%density, 1), SIZE(initial(domain)%density, 2)))
         DO descriptor = skala_density_descriptor, skala_kinetic_descriptor
            IF (descriptor == skala_density_descriptor) THEN
               CALL skala_density_crossings(initial(domain)%density, final(domain)%density, mask)
            ELSE
               CALL skala_density_crossings(initial(domain)%kinetic, final(domain)%kinetic, mask)
            END IF
            DO spin = 1, SIZE(mask, 2)
               DO i = 1, SIZE(mask, 1)
                  IF (.NOT. mask(i, spin)) CYCLE
                  row = row + 1
                  keys(1, row) = domain
                  keys(2:4, row) = initial(domain)%points(:, i)
                  keys(5, row) = spin
                  IF (SIZE(mask, 2) == 1) keys(5, row) = 0
                  keys(6, row) = descriptor
                  IF (descriptor == skala_density_descriptor) THEN
                     values(:, row) = [initial(domain)%density(i, spin), final(domain)%density(i, spin)]
                  ELSE
                     values(:, row) = [initial(domain)%kinetic(i, spin), final(domain)%kinetic(i, spin)]
                  END IF
               END DO
            END DO
         END DO
         DEALLOCATE (mask)
      END DO
      CPASSERT(row == local_count)
      ALLOCATE (counts(state%para_env%num_pe))
      counts = 0
      counts(state%para_env%mepos + 1) = local_count
      CALL state%para_env%sum(counts)
      total = SUM(counts)
      offset = SUM(counts(:state%para_env%mepos))
      ALLOCATE (points(3, total), domains(total), spins(total), descriptors(total), endpoint_values(2, total))
      points = 0
      domains = 0
      spins = 0
      descriptors = 0
      endpoint_values = 0.0_dp
      points(:, offset + 1:offset + local_count) = keys(2:4, :)
      domains(offset + 1:offset + local_count) = keys(1, :)
      spins(offset + 1:offset + local_count) = keys(5, :)
      descriptors(offset + 1:offset + local_count) = keys(6, :)
      endpoint_values(:, offset + 1:offset + local_count) = values
      IF (total > 0) THEN
         CALL state%para_env%sum(points)
         CALL state%para_env%sum(domains)
         CALL state%para_env%sum(spins)
         CALL state%para_env%sum(descriptors)
         CALL state%para_env%sum(endpoint_values)
      END IF
      ALLOCATE (gradient(SIZE(first)))
      CALL qs_skala_scf_evaluate(state, first, energy, gradient)
   END SUBROUTINE qs_skala_scf_crossings

   SUBROUTINE evaluate_physical_input(step, value, context)
      REAL(dp), INTENT(IN)                              :: step
      REAL(dp), INTENT(OUT)                             :: value
      CLASS(*), INTENT(INOUT)                           :: context

      REAL(dp)                                         :: energy
      REAL(dp), ALLOCATABLE                            :: coordinate(:)

      SELECT TYPE (context)
      TYPE IS (physical_root_context)
         ALLOCATE (coordinate(SIZE(context%first)))
         coordinate(:) = context%first + step*context%direction
         CALL qs_skala_scf_evaluate(context%state, coordinate, energy)
         CALL qs_skala_scf_value(context%state, context%domain, context%point, context%spin, context%descriptor, value)
      CLASS DEFAULT
         CPABORT("Invalid native Skala root context.")
      END SELECT
   END SUBROUTINE evaluate_physical_input
END MODULE qs_skala_scf
