!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
!> \brief Independent original-UPF basis representations and their native Gaussian references.
MODULE qs_upf_inverse_basis
   USE admm_types,                      ONLY: admm_type
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE cell_types,                      ONLY: cell_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type
   USE cp_dbcsr_operations,             ONLY: dbcsr_deallocate_matrix_set
   USE input_constants,                 ONLY: uspp_method_inverse
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type
   USE kinds,                           ONLY: dp
   USE lebedev,                         ONLY: deallocate_lebedev_grids,&
                                              init_lebedev_grids
   USE particle_types,                  ONLY: particle_type
   USE pw_grid_types,                   ONLY: pw_grid_type
   USE pw_types,                        ONLY: pw_c1d_gs_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_grid_atom,                    ONLY: grid_atom_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_kinetic,                      ONLY: build_kinetic_matrix
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type
   USE qs_rho_types,                    ONLY: qs_rho_get
   USE qs_upf_inverse_atom_grid,        ONLY: qs_upf_atom_kind_grid,&
                                              qs_upf_reference_rotations
   USE qs_upf_inverse_augmentation,     ONLY: qs_upf_inverse_augmentation_prepare
   USE qs_upf_inverse_background,       ONLY: qs_upf_inverse_background_operator,&
                                              qs_upf_inverse_background_prepare
   USE qs_upf_inverse_charge_metric,    ONLY: qs_upf_inverse_charge_metric_prepare
   USE qs_upf_inverse_grid,             ONLY: qs_upf_inverse_grid_prepare
   USE qs_upf_inverse_gth_local,        ONLY: qs_upf_inverse_gth_local_active
   USE qs_upf_inverse_kpoints,          ONLY: qs_upf_inverse_kpoints_prepare
   USE qs_upf_inverse_scf_xc,           ONLY: qs_upf_inverse_scf_xc_reference
   USE qs_upf_inverse_types,            ONLY: qs_upf_inverse_scf_release,&
                                              qs_upf_inverse_scf_type,&
                                              qs_upf_inverse_type
   USE upf_inverse,                     ONLY: upf_inverse_fields_init
   USE upf_radial_fields,               ONLY: upf_radial_fields_type
   USE xc_derivatives,                  ONLY: xc_functionals_get_needs
   USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: qs_upf_inverse_basis_prepare, qs_upf_inverse_basis_images, qs_upf_inverse_basis_admm
   PUBLIC :: qs_upf_inverse_basis_select

CONTAINS

! **************************************************************************************************
!> \brief Borrow the original-UPF representation matching the operator's native Gaussian basis.
!>        The primary state alone controls SCF activation. An auxiliary state remains disabled
!>        as an SCF owner, but supplies its own reconstruction to auxiliary-basis operators.
!> \param qs_env Native primary and optional ADMM geometry owners.
!> \param basis_type ORB or AUX_FIT, matching the native operator's AO columns.
!> \param state Borrowed selected state, null when the primary UPF method is inactive.
!> \return Whether complete original-UPF operator corrections are required for this basis.
! **************************************************************************************************
   FUNCTION qs_upf_inverse_basis_select(qs_env, basis_type, state) RESULT(active)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      CHARACTER(LEN=*), INTENT(IN)                       :: basis_type
      TYPE(qs_upf_inverse_scf_type), POINTER             :: state
      LOGICAL                                            :: active

      TYPE(admm_type), POINTER                           :: admm
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_upf_inverse_scf_type), POINTER             :: primary

      active = .FALSE.
      NULLIFY (state)
      CALL get_qs_env(qs_env, ks_env=ks_env, admm_env=admm)
      CALL get_ks_env(ks_env, upf_inverse=primary)
      IF (.NOT. ASSOCIATED(primary)) RETURN
      IF (.NOT. primary%enabled) RETURN
      SELECT CASE (TRIM(basis_type))
      CASE ('ORB')
         state => primary
      CASE ('AUX_FIT')
         CPASSERT(ASSOCIATED(admm))
         IF (.NOT. ASSOCIATED(admm%upf_inverse)) THEN
            CPABORT('The original-UPF auxiliary representation is not prepared')
         END IF
         state => admm%upf_inverse
      CASE DEFAULT
         CPABORT('Unsupported Gaussian basis for original-UPF operators')
      END SELECT
      CPASSERT(state%basis_type == basis_type)
      active = .TRUE.
   END FUNCTION qs_upf_inverse_basis_select

! **************************************************************************************************
!> \brief Fixed cell images covering the selected full basis and original radial reconstruction.
!>        ORB and AUX_FIT may need different bounds, including unwrapped atom positions.
!>        Keep this defining image set fixed when evaluating geometry derivatives.
!> \param ks_env Native geometry, basis radii and original pseudopotentials.
!> \param images New complete image list.
!> \param basis_type ORB by default, or the full AUX_FIT basis.
!> \param direct Omit reference-wave differences for the direct representation.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_basis_images(ks_env, images, basis_type, direct)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: images(:, :)
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL             :: basis_type
      LOGICAL, INTENT(IN), OPTIONAL                      :: direct

      CHARACTER(LEN=16)                                  :: role
      INTEGER                                            :: atom, axis, ikind, image, ix, iy, iz, &
                                                            lower(3), upper(3)
      LOGICAL                                            :: use_direct
      REAL(dp)                                           :: domain_lower(3), domain_upper(3), fractional(3), &
                                                            lower_real(3), radius, reach, reciprocal_norm(3), upper_real(3)
      REAL(dp), ALLOCATABLE                              :: radial_radius(:)
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: control
      TYPE(gto_basis_set_type), POINTER                  :: basis
      TYPE(particle_type), POINTER                       :: atoms(:)
      TYPE(qs_kind_type), POINTER                        :: kinds(:)
      TYPE(upf_radial_fields_type)                       :: beta, delta

      role = 'ORB'
      IF (PRESENT(basis_type)) role = basis_type
      CPASSERT(role == 'ORB' .OR. role == 'AUX_FIT')
      use_direct = .FALSE.
      IF (PRESENT(direct)) use_direct = direct
      CALL get_ks_env(ks_env, qs_kind_set=kinds, particle_set=atoms, cell=cell, dft_control=control)
      ALLOCATE (radial_radius(SIZE(kinds)))
      radial_radius = 0.0_dp
      DO ikind = 1, SIZE(kinds)
         IF (kinds(ikind)%ghost) CYCLE
         IF (ASSOCIATED(kinds(ikind)%gth_potential)) THEN
            radial_radius(ikind) = MAX(kinds(ikind)%gth_potential%ppl_radius, kinds(ikind)%gth_potential%ppnl_radius)
            CYCLE
         END IF
         IF (.NOT. ASSOCIATED(kinds(ikind)%sgp_potential)) THEN
            CPABORT("Inverse USPP mixed potentials are not implemented")
         END IF
         IF (.NOT. ALLOCATED(kinds(ikind)%sgp_potential%upf_data)) THEN
            CPABORT("Inverse USPP mixed potentials are not implemented")
         END IF
         IF (kinds(ikind)%sgp_potential%upf_data%number_of_proj == 0) CYCLE
         CALL upf_inverse_fields_init(kinds(ikind)%sgp_potential%upf_data, beta, delta, use_direct)
         radial_radius(ikind) = MAX(beta%radius, delta%radius)
      END DO
      reciprocal_norm = SQRT(SUM(cell%h_inv**2, DIM=2))
      domain_lower = 0.0_dp
      domain_upper = 1.0_dp
      IF (control%qs_control%uspp_method == uspp_method_inverse) THEN
         ! The atom-centered difference-wave integrals extend beyond the primary cell.
         DO atom = 1, SIZE(atoms)
            radius = radial_radius(atoms(atom)%atomic_kind%kind_number)
            IF (radius == 0.0_dp) CYCLE
            fractional = MATMUL(cell%h_inv, atoms(atom)%r)
            domain_lower = MIN(domain_lower, fractional - radius*reciprocal_norm)
            domain_upper = MAX(domain_upper, fractional + radius*reciprocal_norm)
         END DO
      END IF
      lower_real = 0.0_dp
      upper_real = 0.0_dp
      DO atom = 1, SIZE(atoms)
         ikind = atoms(atom)%atomic_kind%kind_number
         CALL get_qs_kind(kinds(ikind), basis_set=basis, basis_type=TRIM(role))
         CPASSERT(ASSOCIATED(basis))
         radius = basis%kind_radius
         IF (.NOT. kinds(ikind)%ghost) THEN
            radius = MAX(radius, radial_radius(ikind))
         END IF
         fractional = MATMUL(cell%h_inv, atoms(atom)%r)
         DO axis = 1, 3
            IF (cell%perd(axis) == 0) CYCLE
            reach = radius*reciprocal_norm(axis)
            ! One additional cell keeps exponentially small Gaussian tails away from image boundaries.
            lower_real(axis) = MIN(lower_real(axis), domain_lower(axis) - fractional(axis) - reach - 1.0_dp)
            upper_real(axis) = MAX(upper_real(axis), domain_upper(axis) - fractional(axis) + reach + 1.0_dp)
         END DO
      END DO
      lower = FLOOR(lower_real)
      upper = CEILING(upper_real)
      ALLOCATE (images(3, PRODUCT(upper - lower + 1)))
      image = 0
      DO iz = lower(3), upper(3)
         DO iy = lower(2), upper(2)
            DO ix = lower(1), upper(1)
               image = image + 1
               images(:, image) = [ix, iy, iz]
            END DO
         END DO
      END DO
   END SUBROUTINE qs_upf_inverse_basis_images

! **************************************************************************************************
!> \brief Prepare the original-UPF auxiliary representation owned by the native ADMM geometry.
!>        Its grids and tasks remain borrowed. Rebuilding native ADMM lists must release it first.
!>        Primary density, operators and SCF activation are unchanged.
!> \param qs_env Native primary and auxiliary geometry owners.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_basis_admm(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      INTEGER                                            :: info
      INTEGER, ALLOCATABLE                               :: images(:, :)
      LOGICAL                                            :: kinetic, laplacian
      TYPE(admm_type), POINTER                           :: admm
      TYPE(dft_control_type), POINTER                    :: control
      TYPE(pw_c1d_gs_type), POINTER                      :: grids(:)
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_upf_inverse_scf_type), POINTER             :: primary
      TYPE(qs_upf_inverse_type), POINTER                 :: frame
      TYPE(section_vals_type), POINTER                   :: functional
      TYPE(xc_rho_cflags_type)                           :: xc_needs

      CALL get_qs_env(qs_env, ks_env=ks_env, admm_env=admm, dft_control=control)
      CPASSERT(ASSOCIATED(admm))
      IF (ASSOCIATED(admm%upf_inverse)) RETURN
      CALL get_ks_env(ks_env, upf_inverse=primary)
      IF (.NOT. ASSOCIATED(primary)) RETURN
      CPASSERT(primary%enabled)
      frame => primary%frame%ao
      IF (ASSOCIATED(primary%bloch%kpoints)) frame => primary%bloch%frames(1)%ao
      kinetic = .FALSE.
      laplacian = .FALSE.
      IF (ASSOCIATED(admm%xc_section_aux)) THEN
         functional => section_vals_get_subs_vals(admm%xc_section_aux, 'XC_FUNCTIONAL')
         xc_needs = xc_functionals_get_needs(functional, control%lsd, calc_potential=.FALSE.)
         kinetic = xc_needs%tau .OR. xc_needs%tau_spin
         laplacian = xc_needs%laplace_rho .OR. xc_needs%laplace_rho_spin
      END IF
      CALL qs_rho_get(admm%rho_aux_fit, rho_g=grids)
      CPASSERT(ASSOCIATED(grids))
      CALL qs_upf_inverse_basis_images(ks_env, images, basis_type='AUX_FIT', direct=frame%direct)
      CALL qs_upf_inverse_basis_prepare(ks_env, grids(1)%pw_grid, images, admm%upf_inverse, info, &
                                        atom_local=frame%atom_local, direct=frame%direct, basis_type='AUX_FIT', &
                                        qs_env=qs_env, kinetic=kinetic, composite_paw=.FALSE., laplacian=laplacian)
      IF (info /= 0) THEN
         CPABORT('Singular or ill-conditioned auxiliary reconstruction')
      END IF
      CALL qs_upf_inverse_charge_metric_prepare(qs_env, admm%upf_inverse, admm%upf_inverse, admm%upf_aux_metric)
      CALL qs_upf_inverse_charge_metric_prepare(qs_env, admm%upf_inverse, primary, admm%upf_mixed_metric)
   END SUBROUTINE qs_upf_inverse_basis_admm

! **************************************************************************************************
!> \brief Prepare an owned orbital or auxiliary basis without installing a native SCF state.
!>        Both bases use the original projectors, augmentation and selected XC reference data.
!>        The caller owns this context and must release it before its borrowed bases and grids.
!> \param ks_env Native geometry and basis owner. Density and operator caches remain unchanged.
!> \param grid Physical-cell source grid borrowed for the lifetime of this context.
!> \param images Fixed physical-cell images for Gaussian and radial fields.
!> \param state New owned representation. It remains disabled as a native SCF owner.
!> \param info Inverse factorization status. Failure releases the incomplete context.
!> \param atom_local Use independent atom-local inverse projectors.
!> \param direct Use original soft orbitals without reconstruction.
!> \param basis_type ORB by default, or AUX_FIT with its actual native auxiliary basis.
!> \param qs_env Matching native owner, required for auxiliary Gaussian references and task lists.
!> \param kinetic Override the primary functional's kinetic-density data requirement.
!> \param composite_paw Override its joint learned PAW valence reconstruction requirement.
!> \param laplacian Override its PAW density-Laplacian requirement.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_basis_prepare(ks_env, grid, images, state, info, atom_local, direct, &
                                           basis_type, qs_env, kinetic, composite_paw, laplacian)
      TYPE(qs_ks_env_type), INTENT(IN), TARGET           :: ks_env
      TYPE(pw_grid_type), INTENT(IN), POINTER            :: grid
      INTEGER, INTENT(IN)                                :: images(:, :)
      TYPE(qs_upf_inverse_scf_type), POINTER             :: state
      INTEGER, INTENT(OUT)                               :: info
      LOGICAL, INTENT(IN), OPTIONAL                      :: atom_local, direct
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL             :: basis_type
      TYPE(qs_environment_type), OPTIONAL, POINTER       :: qs_env
      LOGICAL, INTENT(IN), OPTIONAL                      :: kinetic, composite_paw, laplacian

      CHARACTER(LEN=16)                                  :: role
      INTEGER                                            :: ikind
      LOGICAL                                            :: auxiliary, do_kpoints, gaussian_kinetic, use_atomic, &
                                                            use_composite, use_kinetic, use_laplacian
      REAL(dp), ALLOCATABLE, TARGET                       :: projector_rotations(:, :, :)
      REAL(dp), POINTER                                  :: center_rotations(:, :, :)
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_p_type), POINTER                        :: kinetic_k(:, :)
      TYPE(dft_control_type), POINTER                    :: control
      TYPE(grid_atom_type), POINTER                      :: atom_grid
      TYPE(qs_kind_type), POINTER                        :: kinds(:)
      TYPE(qs_ks_env_type), POINTER                      :: native_ks

      CPASSERT(.NOT. ASSOCIATED(state))
      CPASSERT(ASSOCIATED(grid))
      role = 'ORB'
      IF (PRESENT(basis_type)) role = basis_type
      CPASSERT(role == 'ORB' .OR. role == 'AUX_FIT')
      auxiliary = role == 'AUX_FIT'
      CALL get_ks_env(ks_env, dft_control=control, qs_kind_set=kinds, do_kpoints=do_kpoints, cell=cell)
      DO ikind = 1, SIZE(kinds)
         IF (kinds(ikind)%ghost) CYCLE
         IF (ASSOCIATED(kinds(ikind)%gth_potential)) THEN
            IF (kinds(ikind)%gth_potential%soc .OR. kinds(ikind)%gth_potential%lsdpot) THEN
               CPABORT("Mixed GTH needs scalar spin-independent data")
            END IF
            CYCLE
         END IF
         CPASSERT(ASSOCIATED(kinds(ikind)%sgp_potential))
         CPASSERT(ALLOCATED(kinds(ikind)%sgp_potential%upf_data))
         IF (.NOT. kinds(ikind)%sgp_potential%upf_data%is_ultrasoft .AND. &
             .NOT. kinds(ikind)%sgp_potential%upf_data%is_paw .AND. &
             TRIM(kinds(ikind)%sgp_potential%upf_data%pseudo_type) /= "NC") THEN
            CPABORT("Unsupported inverse UPF type")
         END IF
      END DO
      ALLOCATE (state)
      state%basis_type = role
      state%source_grid => grid
      state%native_gth_local = qs_upf_inverse_gth_local_active(kinds)
      use_kinetic = control%use_kinetic_energy_density
      IF (PRESENT(kinetic)) use_kinetic = kinetic
      use_laplacian = control%use_density_laplacian
      IF (PRESENT(laplacian)) use_laplacian = laplacian
      use_composite = control%use_gauxc
      IF (PRESENT(composite_paw)) use_composite = composite_paw
      CALL qs_upf_reference_rotations(cell, control%qs_control%upf_angular_reference_cell, kinds, &
                                      control%qs_control%gapw .OR. control%qs_control%gapw_xc, state%angular_rotations)
      NULLIFY (center_rotations)
      IF (control%qs_control%gapw .AND. ALLOCATED(state%angular_rotations)) center_rotations => state%angular_rotations
      CALL qs_upf_inverse_scf_xc_reference(kinds, use_kinetic, use_composite, use_laplacian, state)
      use_atomic = .FALSE.
      IF (PRESENT(atom_local)) use_atomic = atom_local
      IF (PRESENT(direct)) THEN
         IF (direct) use_atomic = .TRUE.
         state%native_gaussian = direct .AND. .NOT. control%qs_control%gapw
      END IF
      ! Auxiliary XC needs the same soft Gaussian source as full GAPW, without owning Hartree data.
      state%gaussian_background = auxiliary .AND. control%qs_control%gapw
      gaussian_kinetic = use_atomic .OR. .NOT. auxiliary
      IF (.NOT. (use_atomic .OR. control%qs_control%gapw)) THEN
         IF (ALLOCATED(state%angular_rotations)) THEN
            center_rotations => state%angular_rotations
         ELSE
            CALL qs_upf_reference_rotations(cell, control%qs_control%upf_angular_reference_cell, kinds, &
                                            .TRUE., projector_rotations)
            IF (ALLOCATED(projector_rotations)) center_rotations => projector_rotations
         END IF
         CALL init_lebedev_grids()
         DO ikind = 1, SIZE(kinds)
            IF (kinds(ikind)%ghost .OR. ASSOCIATED(kinds(ikind)%gth_potential)) CYCLE
            ASSOCIATE (pot => kinds(ikind)%sgp_potential%upf_data)
               IF (pot%number_of_proj == 0) CYCLE
            END ASSOCIATE
            CALL qs_upf_atom_kind_grid(kinds(ikind), control%qs_control%gapw_control%quadrature, atom_grid)
         END DO
         CALL deallocate_lebedev_grids()
      END IF
      IF (gaussian_kinetic .OR. state%native_gth_local .OR. state%gaussian_background) THEN
         CALL qs_upf_inverse_background_prepare(ks_env, state%gaussian, soft_basis=control%qs_control%gapw, &
                                                qs_env=qs_env, basis_type=role)
      END IF
      IF (gaussian_kinetic) THEN
         IF (auxiliary) THEN
            native_ks => ks_env
            NULLIFY (kinetic_k)
            CALL build_kinetic_matrix(native_ks, matrixkp_t=kinetic_k, basis_type=role, &
                                      sab_nl=state%gaussian%neighbors)
         ELSE
            CALL get_ks_env(ks_env, kinetic_kp=kinetic_k)
         END IF
         CPASSERT(ASSOCIATED(kinetic_k))
         ALLOCATE (state%gaussian_kinetic(state%gaussian%nao, state%gaussian%nao, SIZE(state%gaussian%weights)))
         CALL qs_upf_inverse_background_operator(state%gaussian, kinetic_k(1, :), state%gaussian_kinetic)
         IF (auxiliary) CALL dbcsr_deallocate_matrix_set(kinetic_k)
      END IF
      IF (do_kpoints) THEN
         IF (state%native_gaussian) THEN
            CALL qs_upf_inverse_kpoints_prepare(ks_env, grid, images, state%bloch, info, use_atomic, direct, &
                                                store_fields=.FALSE., basis_type=role)
         ELSE
            CALL qs_upf_inverse_kpoints_prepare(ks_env, grid, images, state%bloch, info, use_atomic, direct, &
                                                basis_type=role, center_rotations=center_rotations)
         END IF
      ELSE
         IF (state%native_gaussian) THEN
            CALL qs_upf_inverse_grid_prepare(ks_env, grid, images, [0.0_dp, 0.0_dp, 0.0_dp], state%frame, info, &
                                             use_atomic, direct, store_fields=.FALSE., basis_type=role)
         ELSE
            CALL qs_upf_inverse_grid_prepare(ks_env, grid, images, [0.0_dp, 0.0_dp, 0.0_dp], state%frame, info, &
                                             use_atomic, direct, basis_type=role, center_rotations=center_rotations)
         END IF
         IF (info == 0) CALL qs_upf_inverse_augmentation_prepare(ks_env, state%frame%ao, state%augmentation)
      END IF
      IF (info /= 0) CALL qs_upf_inverse_scf_release(state)
   END SUBROUTINE qs_upf_inverse_basis_prepare

END MODULE qs_upf_inverse_basis
