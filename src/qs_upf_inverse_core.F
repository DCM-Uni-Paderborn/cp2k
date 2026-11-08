!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
!> \brief Native partial-core XC source and derivatives from the same original UPF Fourier map.
MODULE qs_upf_inverse_core
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE cell_types,                      ONLY: cell_type
   USE cp_control_types,                ONLY: dft_control_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE particle_types,                  ONLY: particle_type
   USE pw_grid_types,                   ONLY: pw_grid_type
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_integral_ab,&
                                              pw_transfer,&
                                              pw_zero
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE qs_collocate_density,            ONLY: calculate_rho_nlcc
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_integrate_potential_single,   ONLY: integrate_rho_nlcc
   USE qs_kind_types,                   ONLY: has_uspp,&
                                              qs_kind_type
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type,&
                                              set_ks_env
   USE qs_upf_inverse_types,            ONLY: qs_upf_inverse_scf_type
   USE rs_pw_interface,                 ONLY: pw_real_transfer
   USE upf_core,                        ONLY: upf_core_fourier,&
                                              upf_core_kinetic_fourier
   USE upf_types,                       ONLY: upf_nlcc_kinetic,&
                                              upf_uses_model_core_tau
   USE virial_types,                    ONLY: virial_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: qs_upf_inverse_core_collocate, qs_upf_inverse_core_force
   PUBLIC :: qs_upf_inverse_core_kinetic
   PUBLIC :: qs_upf_inverse_core_source_geometry
   PUBLIC :: qs_upf_nc_core_prepare, qs_upf_nc_core_force
   PUBLIC :: qs_upf_core_kinetic_add

CONTAINS

! **************************************************************************************************
!> \brief Replace UPF core fits by original radial cores and retain native GTH core fields.
!>        Called once for each newly prepared inverse geometry. Real and reciprocal core sources
!>        share the native real-field projection and prescribed spin split.
!> \param ks_env Live inverse KS environment with allocated native partial-core storage.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_core_collocate(ks_env)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env

      COMPLEX(dp), ALLOCATABLE                           :: values(:)
      INTEGER                                            :: atom, kind
      LOGICAL, ALLOCATABLE                               :: native_kinds(:)
      TYPE(dft_control_type), POINTER                    :: control
      TYPE(particle_type), POINTER                       :: atoms(:)
      TYPE(pw_c1d_gs_type), POINTER                      :: core_g
      TYPE(pw_r3d_rs_type)                               :: native_core
      TYPE(pw_r3d_rs_type), POINTER                      :: core_r
      TYPE(qs_kind_type), POINTER                        :: kinds(:)
      TYPE(qs_upf_inverse_scf_type), POINTER             :: state

      CALL get_ks_env(ks_env, rho_nlcc=core_r, rho_nlcc_g=core_g, dft_control=control, &
                      particle_set=atoms, qs_kind_set=kinds, upf_inverse=state)
      CPASSERT(ASSOCIATED(state))
      IF (.NOT. ASSOCIATED(core_r)) RETURN
      CPASSERT(ASSOCIATED(core_g))
      CPASSERT(ASSOCIATED(core_r%pw_grid, core_g%pw_grid))
      IF (control%use_kinetic_energy_density) THEN
         DO kind = 1, SIZE(kinds)
            IF (kinds(kind)%ghost) CYCLE
            IF (.NOT. ASSOCIATED(kinds(kind)%sgp_potential)) CYCLE
            IF (.NOT. ALLOCATED(kinds(kind)%sgp_potential%upf_data)) CYCLE
            ASSOCIATE (pot => kinds(kind)%sgp_potential%upf_data)
               IF (state%paw(kind)%active) THEN
                  state%has_core_tau = state%has_core_tau .OR. ANY(state%paw(kind)%xcgrid%core_kinetic%kinetic(:, 2) /= 0.0_dp)
               ELSE
                  IF (.NOT. ALLOCATED(pot%tau_nlcc)) CYCLE
                  state%has_core_tau = state%has_core_tau .OR. ANY(pot%tau_nlcc /= 0.0_dp)
               END IF
            END ASSOCIATE
         END DO
         IF (state%has_core_tau) THEN
            CALL state%core_tau_r%create(core_r%pw_grid)
            CALL state%core_tau_g%create(core_g%pw_grid)
            state%core_tau_g%array = CMPLX(0.0_dp, 0.0_dp, dp)
         END IF
      END IF
      ALLOCATE (values(SIZE(core_g%array)))
      ALLOCATE (native_kinds(SIZE(kinds)), SOURCE=.FALSE.)
      DO kind = 1, SIZE(kinds)
         IF (kinds(kind)%ghost) CYCLE
         IF (ASSOCIATED(kinds(kind)%gth_potential)) native_kinds(kind) = kinds(kind)%gth_potential%nlcc
      END DO
      IF (ANY(native_kinds)) THEN
         CALL native_core%create(core_r%pw_grid)
         CALL calculate_rho_nlcc(native_core, ks_env=ks_env, kind_mask=native_kinds)
      END IF
      core_g%array = CMPLX(0.0_dp, 0.0_dp, dp)
      DO atom = 1, SIZE(atoms)
         kind = atoms(atom)%atomic_kind%kind_number
         IF (kinds(kind)%ghost) CYCLE
         IF (.NOT. ASSOCIATED(kinds(kind)%sgp_potential)) CYCLE
         IF (.NOT. ALLOCATED(kinds(kind)%sgp_potential%upf_data)) CYCLE
         CALL upf_core_fourier(kinds(kind)%sgp_potential%upf_data, core_g%pw_grid%g, atoms(atom)%r, &
                               core_g%pw_grid%vol, values)
         core_g%array = core_g%array + values/REAL(control%nspins, dp)
         IF (state%has_core_tau) THEN
            IF (state%paw(kind)%active) THEN
               CALL upf_core_kinetic_fourier(kinds(kind)%sgp_potential%upf_data, core_g%pw_grid%g, atoms(atom)%r, &
                                             core_g%pw_grid%vol, values, &
                                             core_kinetic=state%paw(kind)%xcgrid%core_kinetic%kinetic(:, 2))
            ELSE IF (upf_uses_model_core_tau(kinds(kind)%sgp_potential%upf_data)) THEN
               CALL upf_core_kinetic_fourier(kinds(kind)%sgp_potential%upf_data, core_g%pw_grid%g, atoms(atom)%r, &
                                             core_g%pw_grid%vol, values)
            ELSE
               CYCLE
            END IF
            state%core_tau_g%array = state%core_tau_g%array + values/REAL(control%nspins, dp)
         END IF
      END DO
      CALL pw_transfer(core_g, core_r)
      IF (ANY(native_kinds)) THEN
         ! Preserve the native real-space source and its existing collocation adjoint.
         core_r%array = core_r%array + native_core%array
         CALL native_core%release()
      END IF
      CALL pw_transfer(core_r, core_g)
      IF (state%has_core_tau) THEN
         CALL pw_transfer(state%core_tau_g, state%core_tau_r)
         CALL pw_transfer(state%core_tau_r, state%core_tau_g)
      END IF
   END SUBROUTINE qs_upf_inverse_core_collocate

! **************************************************************************************************
!> \brief Collocate explicitly supplied NC core tau on the native auxiliary grid.
!> \param ks_env Geometry and model-core source field.
! **************************************************************************************************
   SUBROUTINE qs_upf_nc_core_prepare(ks_env)
      TYPE(qs_ks_env_type), INTENT(INOUT)                :: ks_env

      CHARACTER(LEN=256)                                 :: reason
      COMPLEX(dp), ALLOCATABLE                           :: values(:)
      INTEGER                                            :: atom, kind
      LOGICAL                                            :: available
      REAL(dp), ALLOCATABLE                              :: radial(:)
      TYPE(dft_control_type), POINTER                    :: control
      TYPE(particle_type), POINTER                       :: atoms(:)
      TYPE(pw_c1d_gs_type), POINTER                      :: core_g, tau_g
      TYPE(pw_r3d_rs_type)                               :: tau_r
      TYPE(qs_kind_type), POINTER                        :: kinds(:)

      CALL get_ks_env(ks_env, dft_control=control)
      IF (.NOT. control%use_kinetic_energy_density) RETURN
      IF (.NOT. (control%qs_control%gpw .OR. control%qs_control%gapw_xc .OR. control%qs_control%gapw)) RETURN
      CALL get_ks_env(ks_env, qs_kind_set=kinds, particle_set=atoms, rho_nlcc_g=core_g, rho_nlcc_tau_g=tau_g)
      IF (has_uspp(kinds)) RETURN
      IF (ASSOCIATED(tau_g)) CALL pw_zero(tau_g)
      DO atom = 1, SIZE(atoms)
         kind = atoms(atom)%atomic_kind%kind_number
         IF (kinds(kind)%ghost) CYCLE
         IF (.NOT. ASSOCIATED(kinds(kind)%sgp_potential)) CYCLE
         IF (.NOT. ALLOCATED(kinds(kind)%sgp_potential%upf_data)) CYCLE
         IF (.NOT. ALLOCATED(kinds(kind)%sgp_potential%upf_data%tau_nlcc)) CYCLE
         CALL upf_nlcc_kinetic(kinds(kind)%sgp_potential%upf_data, radial, available, reason)
         IF (.NOT. available) CPABORT(TRIM(reason))
         IF (ALL(radial == 0.0_dp)) CYCLE
         IF (.NOT. ASSOCIATED(tau_g)) THEN
            CPASSERT(ASSOCIATED(core_g))
            ALLOCATE (tau_g)
            CALL tau_g%create(core_g%pw_grid)
            CALL pw_zero(tau_g)
            CALL set_ks_env(ks_env, rho_nlcc_tau_g=tau_g)
         END IF
         IF (.NOT. ALLOCATED(values)) ALLOCATE (values(SIZE(tau_g%array)))
         CALL upf_core_kinetic_fourier(kinds(kind)%sgp_potential%upf_data, tau_g%pw_grid%g, atoms(atom)%r, &
                                       tau_g%pw_grid%vol, values, core_kinetic=radial)
         tau_g%array = tau_g%array + values/REAL(control%nspins, dp)
      END DO
      IF (ASSOCIATED(tau_g)) THEN
         CALL tau_r%create(tau_g%pw_grid)
         CALL pw_transfer(tau_g, tau_r)
         CALL pw_transfer(tau_r, tau_g)
         CALL tau_r%release()
      END IF
   END SUBROUTINE qs_upf_nc_core_prepare

! **************************************************************************************************
!> \brief Add or remove the fixed spin-split core tau on the actual XC grid.
!>        The stored valence kinetic density is restored after XC evaluation. The reciprocal
!>        transfer and its folded potential adjoint retain the native grid projection.
!> \param ks_env Prepared inverse geometry with independently owned core tau.
!> \param kinetic XC kinetic fields, one per spin.
!> \param scale Plus one before XC evaluation, minus one afterwards.
!> \param context Optional independent basis state. Its frozen core tau is prepared on its own grid.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_core_kinetic(ks_env, kinetic, scale, context)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: kinetic(:)
      REAL(dp), INTENT(IN)                               :: scale
      TYPE(qs_upf_inverse_scf_type), INTENT(INOUT), &
         OPTIONAL, TARGET                                :: context

      TYPE(pw_c1d_gs_type), POINTER                      :: native_tau_g
      TYPE(qs_upf_inverse_scf_type), POINTER             :: state

      CALL get_ks_env(ks_env, upf_inverse=state, rho_nlcc_tau_g=native_tau_g)
      IF (PRESENT(context)) THEN
         state => context
         CALL prepare_core_kinetic(ks_env, state)
      END IF
      IF (ASSOCIATED(state)) THEN
         IF (state%has_core_tau) CALL qs_upf_core_kinetic_add(state%core_tau_g, kinetic, scale)
      ELSE IF (ASSOCIATED(native_tau_g)) THEN
         CALL qs_upf_core_kinetic_add(native_tau_g, kinetic, scale)
      END IF
   END SUBROUTINE qs_upf_inverse_core_kinetic

! **************************************************************************************************
!> \brief Transfer the per-spin core field to the XC grid and add it to each spin channel.
!> \param core_tau_g Frozen per-spin core field on the auxiliary grid.
!> \param kinetic Kinetic density in real space on the evaluation grid.
!> \param scale Sign and magnitude of the addition.
! **************************************************************************************************
   SUBROUTINE qs_upf_core_kinetic_add(core_tau_g, kinetic, scale)
      TYPE(pw_c1d_gs_type), INTENT(IN)                   :: core_tau_g
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: kinetic(:)
      REAL(dp), INTENT(IN)                               :: scale

      INTEGER                                            :: spin
      TYPE(pw_c1d_gs_type)                               :: field_g
      TYPE(pw_r3d_rs_type)                               :: field_r

      CPASSERT(SIZE(kinetic) > 0)
      CALL field_g%create(kinetic(1)%pw_grid)
      CALL field_r%create(kinetic(1)%pw_grid)
      CALL pw_real_transfer(core_tau_g, field_g)
      CALL pw_transfer(field_g, field_r)
      DO spin = 1, SIZE(kinetic)
         CPASSERT(ASSOCIATED(kinetic(spin)%pw_grid, field_r%pw_grid))
         CALL pw_axpy(field_r, kinetic(spin), scale)
      END DO
      CALL field_r%release()
      CALL field_g%release()
   END SUBROUTINE qs_upf_core_kinetic_add

! **************************************************************************************************
!> \brief Prepare independent fixed core tau without changing primary partial-core fields.
!>        Preserve the original PAW PS core and the specified model-core fields.
!>        Reference preparation checks data availability and conventions.
!> \param ks_env Original potentials, positions and spin convention.
!> \param state Independent basis state, with its own source grid and reference-wave caches.
! **************************************************************************************************
   SUBROUTINE prepare_core_kinetic(ks_env, state)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(qs_upf_inverse_scf_type), INTENT(INOUT)       :: state

      COMPLEX(dp), ALLOCATABLE                           :: values(:)
      INTEGER                                            :: atom, kind
      TYPE(dft_control_type), POINTER                    :: control
      TYPE(particle_type), POINTER                       :: atoms(:)
      TYPE(pw_r3d_rs_type), POINTER                      :: core_density
      TYPE(qs_kind_type), POINTER                        :: kinds(:)

      IF (.NOT. ALLOCATED(state%kinetic_cache)) RETURN
      IF (ASSOCIATED(state%core_tau_g%pw_grid)) RETURN
      CALL get_ks_env(ks_env, qs_kind_set=kinds, particle_set=atoms, dft_control=control, rho_nlcc=core_density)
      IF (.NOT. ASSOCIATED(core_density)) RETURN
      CPASSERT(ALLOCATED(state%paw))
      DO kind = 1, SIZE(kinds)
         IF (kinds(kind)%ghost) CYCLE
         IF (.NOT. ASSOCIATED(kinds(kind)%sgp_potential)) CYCLE
         IF (.NOT. ALLOCATED(kinds(kind)%sgp_potential%upf_data)) CYCLE
         IF (state%paw(kind)%active) THEN
            state%has_core_tau = state%has_core_tau .OR. ANY(state%paw(kind)%xcgrid%core_kinetic%kinetic(:, 2) /= 0.0_dp)
         ELSE IF (upf_uses_model_core_tau(kinds(kind)%sgp_potential%upf_data)) THEN
            IF (.NOT. ALLOCATED(kinds(kind)%sgp_potential%upf_data%tau_nlcc)) CYCLE
            state%has_core_tau = state%has_core_tau .OR. ANY(kinds(kind)%sgp_potential%upf_data%tau_nlcc /= 0.0_dp)
         END IF
      END DO
      IF (.NOT. state%has_core_tau) RETURN
      CPASSERT(ASSOCIATED(state%source_grid))
      CALL state%core_tau_r%create(state%source_grid)
      CALL state%core_tau_g%create(state%source_grid)
      state%core_tau_g%array = CMPLX(0.0_dp, 0.0_dp, dp)
      ALLOCATE (values(SIZE(state%core_tau_g%array)))
      DO atom = 1, SIZE(atoms)
         kind = atoms(atom)%atomic_kind%kind_number
         IF (kinds(kind)%ghost) CYCLE
         IF (.NOT. ASSOCIATED(kinds(kind)%sgp_potential)) CYCLE
         IF (.NOT. ALLOCATED(kinds(kind)%sgp_potential%upf_data)) CYCLE
         IF (state%paw(kind)%active) THEN
            CALL upf_core_kinetic_fourier(kinds(kind)%sgp_potential%upf_data, state%source_grid%g, atoms(atom)%r, &
                                          state%source_grid%vol, values, &
                                          core_kinetic=state%paw(kind)%xcgrid%core_kinetic%kinetic(:, 2))
         ELSE IF (upf_uses_model_core_tau(kinds(kind)%sgp_potential%upf_data)) THEN
            CALL upf_core_kinetic_fourier(kinds(kind)%sgp_potential%upf_data, state%source_grid%g, atoms(atom)%r, &
                                          state%source_grid%vol, values)
         ELSE
            CYCLE
         END IF
         state%core_tau_g%array = state%core_tau_g%array + values/REAL(control%nspins, dp)
      END DO
      CALL pw_transfer(state%core_tau_g, state%core_tau_r)
      CALL pw_transfer(state%core_tau_r, state%core_tau_g)
   END SUBROUTINE prepare_core_kinetic

! **************************************************************************************************
!> \brief Apply the adjoint of the original partial-core Fourier source to native XC potentials.
!>        The native XC consumer retains its own explicit functional volume/gradient terms.
!>        These derivatives include the core source's phase, reciprocal vector and 1/volume.
!> \param qs_env Native geometry and force/virial storage.
!> \param potential Folded unweighted XC spin potentials on the auxiliary PW grid.
!> \param kinetic_potential Optional folded unweighted spin potentials conjugate to core tau.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_core_force(qs_env, potential, kinetic_potential)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: potential(:)
      TYPE(pw_r3d_rs_type), INTENT(IN), OPTIONAL         :: kinetic_potential(:)

      INTEGER                                            :: kind, spin
      LOGICAL, ALLOCATABLE                               :: native_kinds(:)
      TYPE(qs_kind_type), POINTER                        :: kinds(:)

      CALL get_qs_env(qs_env, qs_kind_set=kinds)
      ALLOCATE (native_kinds(SIZE(kinds)), SOURCE=.FALSE.)
      DO kind = 1, SIZE(kinds)
         IF (kinds(kind)%ghost) CYCLE
         IF (ASSOCIATED(kinds(kind)%gth_potential)) native_kinds(kind) = kinds(kind)%gth_potential%nlcc
      END DO
      IF (ANY(native_kinds)) THEN
         DO spin = 1, SIZE(potential)
            CALL integrate_rho_nlcc(potential(spin), qs_env, kind_mask=native_kinds)
         END DO
      END IF
      CALL core_field_force(qs_env, potential, .FALSE.)
      IF (PRESENT(kinetic_potential)) THEN
         CPASSERT(SIZE(kinetic_potential) == SIZE(potential))
         CALL core_field_force(qs_env, kinetic_potential, .TRUE.)
      END IF
   END SUBROUTINE qs_upf_inverse_core_force

! **************************************************************************************************
!> \brief Differentiate the native NC model-core kinetic density.
!> \param qs_env Native geometry and force/virial storage.
!> \param kinetic_potential Spin potentials on the auxiliary PW grid.
! **************************************************************************************************
   SUBROUTINE qs_upf_nc_core_force(qs_env, kinetic_potential)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: kinetic_potential(:)

      TYPE(pw_c1d_gs_type), POINTER                      :: tau_g

      CALL get_ks_env(qs_env%ks_env, rho_nlcc_tau_g=tau_g)
      IF (ASSOCIATED(tau_g)) CALL core_field_force(qs_env, kinetic_potential, .TRUE.)
   END SUBROUTINE qs_upf_nc_core_force

! **************************************************************************************************
!> \brief Contract frozen-core source geometry without changing accepted forces or virials.
!>        Charge follows the actual native partial-core grid. Core tau follows the independently
!>        prepared UPF context. Original UPF radial fields and native GTH Gaussian cores retain
!>        their spin split and reciprocal-grid transfer. Return integral(v_n delta n_c +
!>        v_tau delta tau_c) with the reference integration measure. The explicit functional
!>        measure and atomic XC residual are separate contributions.
!> \param qs_env Original potentials, native core grid and geometry owners.
!> \param potential Unweighted folded density covectors, one per spin.
!> \param atom_response Total Cartesian atomic motion, including affine cell motion.
!> \param cell_response Lattice derivative at fixed grid topology.
!> \param response Globally reduced source contraction, without a force minus sign.
!> \param kinetic_potential Optional covectors conjugate to the specified frozen-core tau.
!> \param context Optional independent UPF basis and its core-tau convention.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_core_source_geometry(qs_env, potential, atom_response, cell_response, response, &
                                                  kinetic_potential, context)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: potential(:)
      REAL(dp), INTENT(IN)                               :: atom_response(:, :), cell_response(3, 3)
      REAL(dp), INTENT(OUT)                              :: response
      TYPE(pw_r3d_rs_type), INTENT(IN), OPTIONAL         :: kinetic_potential(:)
      TYPE(qs_upf_inverse_scf_type), INTENT(INOUT), &
         OPTIONAL, TARGET                                :: context

      INTEGER                                            :: atom, kind, spin
      LOGICAL                                            :: original_upf, valid
      LOGICAL, ALLOCATABLE                               :: native_kinds(:)
      REAL(dp)                                           :: native_response, strain(3, 3), &
                                                            virial(3, 3)
      REAL(dp), ALLOCATABLE                              :: gradient(:, :), relative_motion(:, :)
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: control
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), POINTER                       :: atoms(:)
      TYPE(pw_c1d_gs_type)                               :: folded_g
      TYPE(pw_r3d_rs_type)                               :: folded
      TYPE(pw_r3d_rs_type), POINTER                      :: core_density
      TYPE(qs_kind_type), POINTER                        :: kinds(:)
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_upf_inverse_scf_type), POINTER             :: state

      CALL get_qs_env(qs_env, ks_env=ks_env, cell=cell, particle_set=atoms, qs_kind_set=kinds, &
                      dft_control=control, para_env=para_env)
      CALL get_ks_env(ks_env, upf_inverse=state, rho_nlcc=core_density)
      original_upf = .FALSE.
      IF (ASSOCIATED(state)) original_upf = state%enabled
      IF (PRESENT(context)) state => context
      valid = ALL(SHAPE(atom_response) == [3, SIZE(atoms)])
      CPASSERT(valid)
      CPASSERT(SIZE(potential) == control%nspins)
      strain = MATMUL(cell_response, cell%h_inv)
      response = 0.0_dp
      IF (ASSOCIATED(core_density)) THEN
         IF (original_upf) THEN
            CALL core_field_geometry(kinds, atoms, state, core_density%pw_grid, potential, &
                                     atom_response, strain, .FALSE., response)
         END IF
         ALLOCATE (native_kinds(SIZE(kinds)), SOURCE=.FALSE.)
         DO kind = 1, SIZE(kinds)
            IF (kinds(kind)%ghost) CYCLE
            IF (ASSOCIATED(kinds(kind)%gth_potential)) native_kinds(kind) = kinds(kind)%gth_potential%nlcc
            IF (.NOT. original_upf .AND. ASSOCIATED(kinds(kind)%sgp_potential)) THEN
               ! Pure NC calculations retain their native fitted partial-core source.
               native_kinds(kind) = kinds(kind)%sgp_potential%has_nlcc
            END IF
         END DO
         IF (ANY(native_kinds)) THEN
            ALLOCATE (gradient(3, SIZE(atoms)), relative_motion(3, SIZE(atoms)))
            DO atom = 1, SIZE(atoms)
               relative_motion(:, atom) = atom_response(:, atom) - MATMUL(strain, atoms(atom)%r)
            END DO
            CALL folded%create(core_density%pw_grid)
            CALL folded_g%create(core_density%pw_grid)
            native_response = 0.0_dp
            DO spin = 1, SIZE(potential)
               CALL core_potential_transfer(potential(spin), folded_g)
               CALL pw_transfer(folded_g, folded)
               CALL integrate_rho_nlcc(folded, qs_env, kind_mask=native_kinds, &
                                       force_response=gradient, virial_response=virial)
               native_response = native_response + SUM(gradient*relative_motion) - SUM(virial*strain)
            END DO
            CALL para_env%sum(native_response)
            response = response + native_response
            CALL folded%release()
            CALL folded_g%release()
         END IF
      END IF
      IF (PRESENT(kinetic_potential)) THEN
         CPASSERT(SIZE(kinetic_potential) == SIZE(potential))
         IF (.NOT. ASSOCIATED(state)) RETURN
         CALL prepare_core_kinetic(ks_env, state)
         IF (.NOT. state%has_core_tau) RETURN
         CALL core_field_geometry(kinds, atoms, state, state%core_tau_g%pw_grid, kinetic_potential, &
                                  atom_response, strain, .TRUE., native_response)
         response = response + native_response
      END IF
   END SUBROUTINE qs_upf_inverse_core_source_geometry

! **************************************************************************************************
!> \brief Differentiate original spherical Fourier cores on their defining grid.
!> \param kinds Original UPF and native kind ownership.
!> \param atoms Original Cartesian centers.
!> \param state Selected PAW PS-core and explicit model-core kinetic fields.
!> \param grid Defining core-density or core-tau grid.
!> \param potential Folded spin covectors on the consumer grid.
!> \param atom_response Total Cartesian center directions.
!> \param strain Lattice derivative times the inverse reference lattice.
!> \param kinetic Select core tau instead of partial-core charge.
!> \param response Globally reduced source contraction, excluding the explicit measure.
! **************************************************************************************************
   SUBROUTINE core_field_geometry(kinds, atoms, state, grid, potential, atom_response, strain, kinetic, response)
      TYPE(qs_kind_type), INTENT(IN)                     :: kinds(:)
      TYPE(particle_type), INTENT(IN)                    :: atoms(:)
      TYPE(qs_upf_inverse_scf_type), POINTER             :: state
      TYPE(pw_grid_type), POINTER                        :: grid
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: potential(:)
      REAL(dp), INTENT(IN)                               :: atom_response(:, :), strain(3, 3)
      LOGICAL, INTENT(IN)                                :: kinetic
      REAL(dp), INTENT(OUT)                              :: response

      COMPLEX(dp), ALLOCATABLE                           :: density(:), dg(:, :), dr(:, :)
      INTEGER                                            :: atom, kind, spin
      REAL(dp)                                           :: trace_strain
      REAL(dp), ALLOCATABLE                              :: reciprocal_response(:, :)
      TYPE(pw_c1d_gs_type)                               :: derivative
      TYPE(pw_c1d_gs_type), ALLOCATABLE                  :: folded(:)

      CPASSERT(ASSOCIATED(grid))
      CALL derivative%create(grid)
      ALLOCATE (folded(SIZE(potential)))
      DO spin = 1, SIZE(potential)
         CALL folded(spin)%create(grid)
         CALL core_potential_transfer(potential(spin), folded(spin))
      END DO
      ALLOCATE (density(SIZE(grid%g, 2)), dr(3, SIZE(grid%g, 2)), dg(3, SIZE(grid%g, 2)))
      ALLOCATE (reciprocal_response, SOURCE=-MATMUL(TRANSPOSE(strain), grid%g))
      trace_strain = strain(1, 1) + strain(2, 2) + strain(3, 3)
      response = 0.0_dp
      DO atom = 1, SIZE(atoms)
         kind = atoms(atom)%atomic_kind%kind_number
         IF (kinds(kind)%ghost) CYCLE
         IF (.NOT. ASSOCIATED(kinds(kind)%sgp_potential)) CYCLE
         IF (.NOT. ALLOCATED(kinds(kind)%sgp_potential%upf_data)) CYCLE
         IF (kinetic) THEN
            IF (state%paw(kind)%active) THEN
               CALL upf_core_kinetic_fourier(kinds(kind)%sgp_potential%upf_data, grid%g, atoms(atom)%r, grid%vol, &
                                             density, dr, dg, core_kinetic=state%paw(kind)%xcgrid%core_kinetic%kinetic(:, 2))
            ELSE IF (upf_uses_model_core_tau(kinds(kind)%sgp_potential%upf_data)) THEN
               CALL upf_core_kinetic_fourier(kinds(kind)%sgp_potential%upf_data, grid%g, atoms(atom)%r, grid%vol, &
                                             density, dr, dg)
            ELSE
               CYCLE
            END IF
         ELSE
            CALL upf_core_fourier(kinds(kind)%sgp_potential%upf_data, grid%g, atoms(atom)%r, grid%vol, density, dr, dg)
         END IF
         derivative%array = (SUM(dr*SPREAD(atom_response(:, atom), 2, SIZE(density)), DIM=1) + &
                             SUM(dg*reciprocal_response, DIM=1) - trace_strain*density)/REAL(SIZE(potential), dp)
         DO spin = 1, SIZE(potential)
            response = response + pw_integral_ab(folded(spin), derivative)
         END DO
      END DO
      DO spin = 1, SIZE(potential)
         CALL folded(spin)%release()
      END DO
      CALL derivative%release()
   END SUBROUTINE core_field_geometry

! **************************************************************************************************
!> \brief Transpose the existing reciprocal transfer to a core source's defining grid.
!> \param potential Unweighted real covector on the consumer grid.
!> \param folded Created Fourier covector on the source grid.
! **************************************************************************************************
   SUBROUTINE core_potential_transfer(potential, folded)
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: potential
      TYPE(pw_c1d_gs_type), INTENT(INOUT)                :: folded

      TYPE(pw_c1d_gs_type)                               :: fine

      CPASSERT(ASSOCIATED(potential%pw_grid))
      IF (ASSOCIATED(potential%pw_grid, folded%pw_grid)) THEN
         CALL pw_transfer(potential, folded)
      ELSE
         CALL fine%create(potential%pw_grid)
         CALL pw_transfer(potential, fine)
         CALL pw_real_transfer(fine, folded)
         CALL fine%release()
      END IF
   END SUBROUTINE core_potential_transfer

! **************************************************************************************************
!> \brief Contract one independent spherical core field's center and strain response.
!> \param qs_env ...
!> \param potential ...
!> \param kinetic Whether to differentiate canonical core tau instead of core charge.
! **************************************************************************************************
   SUBROUTINE core_field_force(qs_env, potential, kinetic)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: potential(:)
      LOGICAL, INTENT(IN)                                :: kinetic

      COMPLEX(dp), ALLOCATABLE                           :: density(:), dg(:, :), dr(:, :)
      INTEGER                                            :: atom, axis, index, kind, other, spin
      INTEGER, ALLOCATABLE                               :: atom_of_kind(:), kind_of(:)
      LOGICAL                                            :: paw
      REAL(dp)                                           :: gradient, spin_factor
      TYPE(atomic_kind_type), POINTER                    :: atomic_kinds(:)
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), POINTER                       :: atoms(:)
      TYPE(pw_c1d_gs_type)                               :: derivative_g, potential_g
      TYPE(qs_force_type), POINTER                       :: force(:)
      TYPE(qs_kind_type), POINTER                        :: kinds(:)
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_upf_inverse_scf_type), POINTER             :: state
      TYPE(virial_type), POINTER                         :: virial

      CALL get_qs_env(qs_env, particle_set=atoms, qs_kind_set=kinds, force=force, &
                      virial=virial, atomic_kind_set=atomic_kinds, para_env=para_env)
      CALL get_qs_env(qs_env, ks_env=ks_env)
      CALL get_ks_env(ks_env, upf_inverse=state)
      CPASSERT(kinetic .OR. ASSOCIATED(state))
      CPASSERT(SIZE(potential) > 0 .AND. ASSOCIATED(force))
      CALL get_atomic_kind_set(atomic_kinds, atom_of_kind=atom_of_kind, kind_of=kind_of)
      CALL potential_g%create(potential(1)%pw_grid)
      CALL derivative_g%create(potential(1)%pw_grid)
      ALLOCATE (density(SIZE(derivative_g%array)), dr(3, SIZE(derivative_g%array)), &
                dg(3, SIZE(derivative_g%array)))
      spin_factor = 1.0_dp/REAL(SIZE(potential), dp)
      DO atom = 1, SIZE(atoms)
         kind = kind_of(atom)
         index = atom_of_kind(atom)
         IF (kinds(kind)%ghost) CYCLE
         IF (.NOT. ASSOCIATED(kinds(kind)%sgp_potential)) CYCLE
         IF (.NOT. ALLOCATED(kinds(kind)%sgp_potential%upf_data)) CYCLE
         IF (.NOT. kinds(kind)%sgp_potential%upf_data%core_correction) CYCLE
         IF (kinetic) THEN
            paw = .FALSE.
            IF (ASSOCIATED(state)) paw = state%paw(kind)%active
            IF (paw) THEN
               CALL upf_core_kinetic_fourier(kinds(kind)%sgp_potential%upf_data, derivative_g%pw_grid%g, atoms(atom)%r, &
                                             derivative_g%pw_grid%vol, density, dr, dg, &
                                             core_kinetic=state%paw(kind)%xcgrid%core_kinetic%kinetic(:, 2))
            ELSE IF (upf_uses_model_core_tau(kinds(kind)%sgp_potential%upf_data)) THEN
               CALL upf_core_kinetic_fourier(kinds(kind)%sgp_potential%upf_data, derivative_g%pw_grid%g, atoms(atom)%r, &
                                             derivative_g%pw_grid%vol, density, dr, dg)
            ELSE
               CYCLE
            END IF
         ELSE
            CALL upf_core_fourier(kinds(kind)%sgp_potential%upf_data, derivative_g%pw_grid%g, atoms(atom)%r, &
                                  derivative_g%pw_grid%vol, density, dr, dg)
         END IF
         DO spin = 1, SIZE(potential)
            CPASSERT(ASSOCIATED(potential(spin)%pw_grid, potential_g%pw_grid))
            CALL pw_transfer(potential(spin), potential_g)
            DO axis = 1, 3
               derivative_g%array = spin_factor*dr(axis, :)
               gradient = pw_integral_ab(potential_g, derivative_g)
               IF (para_env%is_source()) force(kind)%gth_nlcc(axis, index) = &
                  force(kind)%gth_nlcc(axis, index) + gradient
            END DO
            IF (.NOT. virial%pv_availability .OR. virial%pv_numer) CYCLE
            DO axis = 1, 3
               DO other = 1, 3
                  derivative_g%array = atoms(atom)%r(other)*dr(axis, :) - &
                                       derivative_g%pw_grid%g(axis, :)*dg(other, :)
                  IF (axis == other) derivative_g%array = derivative_g%array - density
                  derivative_g%array = spin_factor*derivative_g%array
                  gradient = pw_integral_ab(potential_g, derivative_g)
                  IF (para_env%is_source()) THEN
                     virial%pv_nlcc(axis, other) = virial%pv_nlcc(axis, other) - gradient
                     virial%pv_virial(axis, other) = virial%pv_virial(axis, other) - gradient
                  END IF
               END DO
            END DO
         END DO
      END DO
      CALL potential_g%release()
      CALL derivative_g%release()
   END SUBROUTINE core_field_force

END MODULE qs_upf_inverse_core
