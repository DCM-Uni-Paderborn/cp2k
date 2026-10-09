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
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: fourpi
   USE mathlib,                         ONLY: invmat
   USE message_passing,                 ONLY: mp_para_env_type
   USE orbital_pointers,                ONLY: indso
   USE particle_types,                  ONLY: particle_type
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_grid_types,                   ONLY: FULLSPACE,&
                                              HALFSPACE,&
                                              PW_MODE_LOCAL
   USE pw_methods,                      ONLY: pw_integral_ab,&
                                              pw_scale,&
                                              pw_transfer
   USE pw_poisson_methods,              ONLY: pw_poisson_rebuild
   USE pw_poisson_types,                ONLY: PERIODIC3D,&
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
                                              qs_upf_inverse_background_prepare
   USE qs_upf_inverse_density,          ONLY: qs_upf_inverse_density_gapw,&
                                              qs_upf_inverse_density_gapw_integrate
   USE qs_upf_inverse_gaussian_pairs,   ONLY: qs_upf_inverse_gaussian_pairs_has_derivatives,&
                                              qs_upf_inverse_gaussian_pairs_prepare
   USE qs_upf_inverse_grid,             ONLY: qs_upf_inverse_matrix_gather
   USE qs_upf_inverse_kpoints,          ONLY: qs_upf_inverse_kpoints_density
   USE qs_upf_inverse_scf,              ONLY: qs_upf_inverse_scf_install_metric,&
                                              qs_upf_inverse_scf_integrate_background
   USE qs_upf_inverse_types,            ONLY: qs_upf_inverse_augmentation_type,&
                                              qs_upf_inverse_coulomb_release,&
                                              qs_upf_inverse_coulomb_type,&
                                              qs_upf_inverse_scf_type,&
                                              qs_upf_inverse_type
   USE upf_inverse_coulomb,             ONLY: upf_inverse_coulomb_center_type,&
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
   PUBLIC :: qs_upf_inverse_coulomb_response, qs_upf_inverse_coulomb_metric_response

CONTAINS

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
      INTEGER                                            :: a, atom, channel, k, kind, nao, ncenter, &
                                                            nh, nr
      LOGICAL                                            :: paw_atom
      REAL(dp)                                           :: alpha, total
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
      IF (poisson%green_fft%method /= PERIODIC3D .OR. poisson%parameters%has_dielectric .OR. &
          poisson%green_fft%p3m .OR. poisson%green_fft%sr_screening) THEN
         CPABORT("Inverse GAPW Coulomb requires periodic Poisson")
      END IF
      CALL calculate_rho_core(map%ionic, total, qs_env, only_nopaw=.FALSE.)
      ALLOCATE (map%multiplicity(SIZE(map%ionic%array)), map%poisson_weights(SIZE(map%ionic%array)))
      map%multiplicity = 1.0_dp
      IF (map%grid%grid_span == HALFSPACE) THEN
         map%multiplicity = 2.0_dp
         IF (map%grid%have_g0) map%multiplicity(1) = 1.0_dp
      END IF
      map%poisson_weights(:) = map%grid%vol*map%multiplicity*REAL(poisson%green_fft%influence_fn%array, dp)
      CPASSERT(ALL(AIMAG(poisson%green_fft%influence_fn%array) == 0.0_dp))
      CPASSERT(ALL(map%poisson_weights >= 0.0_dp))
      IF (ASSOCIATED(state%bloch%kpoints)) THEN
         ALLOCATE (map%weights, SOURCE=state%bloch%weights)
      ELSE
         ALLOCATE (map%weights(1))
         map%weights = 1.0_dp
      END IF
      CALL select_frame(state, 1, frame, augmentation)
      CPASSERT(ALL(frame%periodicity == 1))
      nao = SIZE(frame%gaussian, 2)
      map%cell = frame%cell
      ALLOCATE (map%positions, SOURCE=frame%positions)
      ncenter = SIZE(frame%center_atoms)
      DO atom = 1, SIZE(atoms)
         IF (ANY(frame%center_atoms == atom)) CYCLE
         kind = atoms(atom)%atomic_kind%kind_number
         CALL get_qs_kind(kinds(kind), paw_atom=paw_atom)
         IF (paw_atom) ncenter = ncenter + 1
      END DO
      ALLOCATE (map%centers(ncenter), map%pairs(ncenter), map%center_atoms(ncenter), &
                map%center_kinds(ncenter), map%inverse_centers(ncenter))
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
            CALL qs_upf_inverse_density_gapw_integrate(frame, augmentation, map%inverse_centers(a), kinds(kind), &
                                                       map%pairs(a), unit, ao_operator, map%center_atoms(a))
            map%metric(:, :, k) = map%metric(:, :, k) + ao_operator
         END DO
         DEALLOCATE (unit)
      END DO
      CALL upf_inverse_coulomb_prepare(map%centers, map%cell, map%grid%g, map%grid%vol, map%maps)
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

      COMPLEX(dp), ALLOCATABLE                           :: adjoint(:), ao_operator(:, :), &
                                                            background(:)
      INTEGER                                            :: a, k, kind
      REAL(dp)                                           :: electronic_charge, total
      REAL(dp), ALLOCATABLE                              :: radial(:, :)
      TYPE(qs_kind_type), POINTER                        :: kinds(:)
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_upf_inverse_augmentation_type), POINTER    :: augmentation
      TYPE(qs_upf_inverse_type), POINTER                 :: frame
      TYPE(upf_inverse_coulomb_potential_type), &
         ALLOCATABLE                                     :: potentials(:)

      CPASSERT(map%ready)
      CPASSERT(ALL(SHAPE(matrix) == SHAPE(map%metric)))
      CPASSERT(ALL(SHAPE(ao_operators) == SHAPE(map%metric)))
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
                                             map%pairs(a), matrix(:, :, k), radial, atom_index=map%center_atoms(a))
            map%centers(a)%density(:, :) = map%centers(a)%density + map%weights(k)*radial
         END DO
         electronic_charge = electronic_charge + SQRT(fourpi)*SUM(map%centers(a)%weights*map%centers(a)%density(:, 1))
         DEALLOCATE (radial)
      END DO
      IF (PRESENT(charge)) charge = electronic_charge
      ALLOCATE (background, SOURCE=map%source_g%array + map%ionic%array)
      ALLOCATE (adjoint(SIZE(background)), ao_operator(SIZE(matrix, 1), SIZE(matrix, 2)))
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
      CALL qs_upf_inverse_background_integrate(map%background, qs_env, map%potential, ao_operators)
      DO k = 1, SIZE(map%weights)
         CALL select_frame(state, k, frame, augmentation)
         DO a = 1, SIZE(map%centers)
            kind = map%center_kinds(a)
            CALL qs_upf_inverse_density_gapw_integrate(frame, augmentation, map%inverse_centers(a), kinds(kind), &
                                                       map%pairs(a), potentials(a)%values, ao_operator, map%center_atoms(a))
            ao_operators(:, :, k) = ao_operators(:, :, k) + ao_operator
         END DO
      END DO
   END SUBROUTINE qs_upf_inverse_coulomb_evaluate

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
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_coulomb_response(qs_env, state, map, matrix, atom_response, cell_response, &
                                              response, charge_response)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_upf_inverse_scf_type), INTENT(IN), TARGET  :: state
      TYPE(qs_upf_inverse_coulomb_type), INTENT(INOUT)   :: map
      COMPLEX(dp), INTENT(IN)                            :: matrix(:, :, :)
      REAL(dp), INTENT(IN)                               :: atom_response(:, :), cell_response(3, 3)
      REAL(dp), INTENT(OUT)                              :: response
      REAL(dp), INTENT(OUT), OPTIONAL                    :: charge_response

      COMPLEX(dp), ALLOCATABLE                           :: adjoint(:), background(:), &
                                                            background_response(:), &
                                                            operators(:, :, :)
      INTEGER                                            :: a, g, info
      REAL(dp)                                           :: energy, g2, inverse_cell(3, 3), &
                                                            ionic_virial(3, 3), kernel_response, &
                                                            strain(3, 3), trace_strain, &
                                                            virial(3, 3)
      REAL(dp), ALLOCATABLE                              :: forces(:, :), gvecs_response(:, :), &
                                                            ionic_forces(:, :), &
                                                            poisson_weights_response(:)
      TYPE(pw_r3d_rs_type)                               :: unit
      TYPE(upf_inverse_coulomb_center_type), ALLOCATABLE :: centers_response(:)
      TYPE(upf_inverse_coulomb_potential_type), &
         ALLOCATABLE                                     :: potentials(:)

      CPASSERT(map%ready)
      CPASSERT(ALL(SHAPE(atom_response) == SHAPE(map%positions)))
      CALL prepare_pair_responses(qs_env, state, map)
      ALLOCATE (operators, MOLD=map%metric)
      CALL qs_upf_inverse_coulomb_evaluate(qs_env, state, map, matrix, energy, operators)
      inverse_cell = map%cell
      CALL invmat(inverse_cell, info)
      CPASSERT(info == 0)
      strain = MATMUL(cell_response, inverse_cell)
      trace_strain = strain(1, 1) + strain(2, 2) + strain(3, 3)
      ALLOCATE (forces, MOLD=atom_response)
      ALLOCATE (ionic_forces, MOLD=atom_response)
      CALL qs_upf_inverse_background_geometry(map%background, qs_env, matrix, map%potential, forces, virial)
      CALL qs_upf_inverse_coulomb_ionic_geometry(qs_env, map%potential, ionic_forces, ionic_virial)
      ! Fourier coefficients use a grid average. Its source derivative has no volume-measure term.
      response = -SUM((forces + ionic_forces)*(atom_response - MATMUL(strain, map%positions))) - &
                 SUM((virial + ionic_virial)*strain)
      IF (PRESENT(charge_response)) THEN
         CALL unit%create(map%grid)
         unit%array = 1.0_dp
         CALL qs_upf_inverse_background_geometry(map%background, qs_env, matrix, unit, forces, virial)
         charge_response = -SUM(forces*(atom_response - MATMUL(strain, map%positions))) - &
                           SUM(virial*strain) + trace_strain*pw_integral_ab(map%source, unit)
         CALL unit%release()
      END IF
      CALL center_density_response(qs_env, state, map, matrix, atom_response, cell_response, centers_response)
      IF (PRESENT(charge_response)) THEN
         DO a = 1, SIZE(centers_response)
            charge_response = charge_response + SQRT(fourpi)* &
                              SUM(map%centers(a)%weights*centers_response(a)%density(:, 1))
         END DO
      END IF
      ALLOCATE (gvecs_response(3, SIZE(map%grid%g, 2)), poisson_weights_response(SIZE(map%poisson_weights)))
      gvecs_response(:, :) = -MATMUL(TRANSPOSE(strain), map%grid%g)
      poisson_weights_response(:) = trace_strain*map%poisson_weights
      DO g = 1, SIZE(map%poisson_weights)
         g2 = SUM(map%grid%g(:, g)**2)
         IF (g2 == 0.0_dp) CYCLE
         poisson_weights_response(g) = poisson_weights_response(g) - &
                                       2.0_dp*map%poisson_weights(g)*SUM(map%grid%g(:, g)*gvecs_response(:, g))/g2
      END DO
      ALLOCATE (background, SOURCE=map%source_g%array + map%ionic%array)
      ALLOCATE (adjoint(SIZE(background)), background_response(SIZE(background)))
      background_response = CMPLX(0.0_dp, 0.0_dp, dp)
      IF (map%grid%para%mode == PW_MODE_LOCAL) THEN
         CALL upf_inverse_coulomb_periodic(map%centers, map%cell, map%grid%g, map%grid%vol, &
                                           background, map%poisson_weights, energy, potentials, adjoint, &
                                           response=kernel_response, center_response=centers_response, &
                                           cell_response=cell_response, &
                                           gvecs_response=gvecs_response, volume_response=trace_strain*map%grid%vol, &
                                           background_response=background_response, &
                                           poisson_weights_response=poisson_weights_response, map=map%maps)
      ELSE
         CALL upf_inverse_coulomb_periodic(map%centers, map%cell, map%grid%g, map%grid%vol, &
                                           background, map%poisson_weights, energy, potentials, adjoint, &
                                           comm=map%grid%para%group%mp_comm_type, response=kernel_response, &
                                           center_response=centers_response, cell_response=cell_response, &
                                           gvecs_response=gvecs_response, volume_response=trace_strain*map%grid%vol, &
                                           background_response=background_response, &
                                           poisson_weights_response=poisson_weights_response, map=map%maps)
      END IF
      response = response + kernel_response
   END SUBROUTINE qs_upf_inverse_coulomb_response

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
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_coulomb_metric_response(qs_env, state, map, matrix, atom_response, cell_response, response)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_upf_inverse_scf_type), INTENT(IN), TARGET  :: state
      TYPE(qs_upf_inverse_coulomb_type), INTENT(INOUT)   :: map
      COMPLEX(dp), INTENT(IN)                            :: matrix(:, :, :)
      REAL(dp), INTENT(IN)                               :: atom_response(:, :), cell_response(3, 3)
      REAL(dp), INTENT(OUT)                              :: response

      INTEGER                                            :: a, info
      REAL(dp)                                           :: inverse_cell(3, 3), strain(3, 3), &
                                                            trace_strain, virial(3, 3)
      REAL(dp), ALLOCATABLE                              :: forces(:, :)
      TYPE(pw_c1d_gs_type)                               :: background_g
      TYPE(pw_r3d_rs_type)                               :: background, unit
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(upf_inverse_coulomb_center_type), ALLOCATABLE :: centers_response(:)

      CPASSERT(map%ready)
      CPASSERT(ALL(SHAPE(matrix) == SHAPE(map%metric)))
      CPASSERT(ALL(SHAPE(atom_response) == SHAPE(map%positions)))
      CALL get_qs_env(qs_env, ks_env=ks_env)
      CALL prepare_pair_responses(qs_env, state, map)
      inverse_cell = map%cell
      CALL invmat(inverse_cell, info)
      CPASSERT(info == 0)
      strain = MATMUL(cell_response, inverse_cell)
      trace_strain = strain(1, 1) + strain(2, 2) + strain(3, 3)
      CALL background%create(map%grid)
      CALL background_g%create(map%grid)
      CALL unit%create(map%grid)
      unit%array = 1.0_dp
      CALL qs_upf_inverse_background_collocate(map%background, ks_env, matrix, background, background_g)
      ALLOCATE (forces, MOLD=atom_response)
      CALL qs_upf_inverse_background_geometry(map%background, qs_env, matrix, unit, forces, virial)
      response = -SUM(forces*(atom_response - MATMUL(strain, map%positions))) - SUM(virial*strain) + &
                 trace_strain*pw_integral_ab(background, unit)
      CALL center_density_response(qs_env, state, map, matrix, atom_response, cell_response, centers_response)
      DO a = 1, SIZE(centers_response)
         response = response + SQRT(fourpi)*SUM(map%centers(a)%weights*centers_response(a)%density(:, 1))
      END DO
      CALL unit%release()
      CALL background%release()
      CALL background_g%release()
   END SUBROUTINE qs_upf_inverse_coulomb_metric_response

! **************************************************************************************************
!> \brief Complete center density derivative at fixed native radial nodes and measures.
!> \param qs_env Native kinds with their original radial and Gaussian one-center fields.
!> \param state Original defining inverse frames.
!> \param map Current complete source layout and derivative pair coefficients.
!> \param matrix Hermitian spin-summed matrices per global k, without mesh weights.
!> \param atom_response Physical atomic motion.
!> \param cell_response Physical lattice derivative.
!> \param centers_response New center layout with density and position responses.
! **************************************************************************************************
   SUBROUTINE center_density_response(qs_env, state, map, matrix, atom_response, cell_response, centers_response)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_upf_inverse_scf_type), INTENT(IN), TARGET  :: state
      TYPE(qs_upf_inverse_coulomb_type), INTENT(IN)      :: map
      COMPLEX(dp), INTENT(IN)                            :: matrix(:, :, :)
      REAL(dp), INTENT(IN)                               :: atom_response(:, :), cell_response(3, 3)
      TYPE(upf_inverse_coulomb_center_type), &
         ALLOCATABLE, INTENT(OUT)                        :: centers_response(:)

      COMPLEX(dp), ALLOCATABLE                           :: coefficients_response(:, :)
      INTEGER                                            :: a, info, k, kind
      REAL(dp)                                           :: inverse_cell(3, 3), strain(3, 3), &
                                                            trace_strain
      REAL(dp), ALLOCATABLE                              :: point_response(:, :), radial(:, :), &
                                                            radial_response(:, :), &
                                                            weights_response(:)
      TYPE(qs_kind_type), POINTER                        :: kinds(:)
      TYPE(qs_upf_inverse_augmentation_type), POINTER    :: augmentation
      TYPE(qs_upf_inverse_type), POINTER                 :: frame

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
      DO k = 1, SIZE(map%weights)
         CALL select_frame(state, k, frame, augmentation)
         ALLOCATE (coefficients_response, MOLD=frame%coefficients)
         ALLOCATE (point_response(3, SIZE(frame%weights)), weights_response(SIZE(frame%weights)))
         point_response(:, :) = MATMUL(strain, frame%points)
         weights_response(:) = trace_strain*frame%weights
         CALL qs_upf_inverse_response(frame, point_response, atom_response, cell_response, &
                                      [0.0_dp, 0.0_dp, 0.0_dp], weights_response, coefficients_response)
         DO a = 1, SIZE(centers_response)
            kind = map%center_kinds(a)
            ALLOCATE (radial, MOLD=map%centers(a)%density)
            ALLOCATE (radial_response, MOLD=radial)
         CALL qs_upf_inverse_density_gapw(frame, augmentation, map%inverse_centers(a), kinds(kind), map%pairs(a), matrix(:, :, k), &
                                             radial, radial_response, atom_response, cell_response, &
                                             [0.0_dp, 0.0_dp, 0.0_dp], coefficients_response, map%center_atoms(a))
            centers_response(a)%density(:, :) = centers_response(a)%density + map%weights(k)*radial_response
            DEALLOCATE (radial, radial_response)
         END DO
         DEALLOCATE (point_response, weights_response, coefficients_response)
      END DO
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
