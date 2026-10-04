!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Connect native USPP geometry, density, and operators to Quickstep.
! **************************************************************************************************
MODULE qs_uspp_methods
   USE atom_upf,                        ONLY: upf_nlcc_kinetic
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE cell_types,                      ONLY: cell_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_add, dbcsr_copy, dbcsr_desymmetrize, dbcsr_get_block_p, dbcsr_get_info, &
        dbcsr_p_type, dbcsr_release, dbcsr_scale, dbcsr_set, dbcsr_type
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE input_constants,                 ONLY: do_ppl_analytic
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_type
   USE message_passing,                 ONLY: mp_comm_type,&
                                              mp_comm_unequal
   USE particle_types,                  ONLY: particle_type
   USE paw_core_kinetic,                ONLY: paw_core_kinetic_build_ae,&
                                              paw_core_kinetic_cartesian,&
                                              paw_core_kinetic_fourier,&
                                              paw_core_kinetic_type
   USE paw_onecenter,                   ONLY: paw_onecenter_core_exchange,&
                                              paw_onecenter_core_exchange_init,&
                                              paw_onecenter_core_exchange_prepare,&
                                              paw_onecenter_exchange,&
                                              paw_onecenter_hartree,&
                                              paw_onecenter_init
   USE paw_xc,                          ONLY: paw_xc_evaluate,&
                                              paw_xc_init
   USE pw_grid_types,                   ONLY: FULLSPACE,&
                                              HALFSPACE,&
                                              pw_grid_type
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_derive,&
                                              pw_scale,&
                                              pw_transfer
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type,&
                                              set_ks_env
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_uspp_projectors,              ONLY: qs_uspp_projector_gradient,&
                                              qs_uspp_projector_images,&
                                              qs_uspp_projector_mixed_images
   USE qs_uspp_types,                   ONLY: qs_uspp_init,&
                                              qs_uspp_type,&
                                              qs_uspp_update,&
                                              qs_uspp_update_images
   USE sap_kind_types,                  ONLY: build_sap_ints,&
                                              build_sap_ints_basis,&
                                              release_sap_int,&
                                              sap_int_type,&
                                              sap_sort
   USE upf_core_kinetic,                ONLY: upf_core_kinetic_cartesian,&
                                              upf_core_kinetic_fourier,&
                                              upf_core_kinetic_init,&
                                              upf_core_kinetic_type
   USE upf_local_integrals,             ONLY: upf_local_fourier,&
                                              upf_local_init
   USE uspp_augmentation,               ONLY: &
        paw_reconstruction_init, upf_kinetic_reconstruction_init, upf_radial_scalar_init, &
        uspp_augmentation_collocate, uspp_augmentation_collocate_atom, uspp_augmentation_init, &
        uspp_augmentation_integrate, uspp_augmentation_integrate_atom, &
        uspp_augmentation_prepare_fourier, uspp_augmentation_prepare_realspace, &
        uspp_augmentation_type
   USE uspp_projector_operators,        ONLY: uspp_projector_add_operator,&
                                              uspp_projector_create_onsite,&
                                              uspp_projector_density,&
                                              uspp_projector_mixed_pullback,&
                                              uspp_projector_pullback
   USE virial_types,                    ONLY: virial_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_uspp_methods'
   PUBLIC :: qs_uspp_prepare, qs_uspp_init_ao_matrix, qs_uspp_collocate, qs_uspp_integrate
   PUBLIC :: qs_uspp_collocate_state, qs_uspp_integrate_state, qs_uspp_init_matrix_state
   PUBLIC :: qs_uspp_overlap_adjoint, qs_uspp_atom_grid_apply_state
   PUBLIC :: qs_uspp_core_matrix, qs_uspp_core_forces, qs_paw_correction, qs_paw_exchange
   PUBLIC :: qs_paw_core_tau_present, qs_paw_core_tau_field, qs_paw_core_tau_forces
   PUBLIC :: qs_paw_core_tau_atom_collocate, qs_paw_core_tau_atom_integrate
   PUBLIC :: qs_upf_core_tau_present, qs_upf_core_tau_field, qs_upf_core_tau_forces
   PUBLIC :: qs_upf_core_tau_atom_collocate, qs_upf_core_tau_atom_integrate
   PUBLIC :: qs_uspp_atom_grid_type, qs_uspp_atom_grid_prepare, qs_uspp_atom_grid_collocate, &
             qs_uspp_atom_grid_integrate, qs_uspp_atom_grid_apply, qs_uspp_atom_grid_apply_ks
   PUBLIC :: qs_uspp_fixed_atom_grid_type, qs_uspp_fixed_atom_grid_prepare, qs_uspp_fixed_atom_grid_integrate
   PUBLIC :: qs_uspp_prepare_global_xc
   PUBLIC :: qs_uspp_kinetic_atom_grid_collocate, qs_uspp_kinetic_atom_grid_integrate
   PUBLIC :: qs_uspp_global_atom_grid_collocate, qs_uspp_global_atom_grid_integrate

   TYPE qs_uspp_atom_grid_type
      LOGICAL :: calculate_forces = .FALSE.
      REAL(KIND=dp), ALLOCATABLE :: density(:, :, :), potential(:, :, :), center_gradient(:, :)
      REAL(KIND=dp) :: center_strain(3, 3) = 0.0_dp
   END TYPE qs_uspp_atom_grid_type

   TYPE qs_uspp_fixed_atom_grid_type
      TYPE(qs_uspp_type) :: geometry
      TYPE(qs_uspp_atom_grid_type) :: fields
      INTEGER, ALLOCATABLE :: source_centers(:)
   END TYPE qs_uspp_fixed_atom_grid_type

   INTERFACE qs_uspp_collocate
      MODULE PROCEDURE qs_uspp_collocate_single, qs_uspp_collocate_images
   END INTERFACE
   INTERFACE qs_uspp_integrate
      MODULE PROCEDURE qs_uspp_integrate_single, qs_uspp_integrate_images
   END INTERFACE
   INTERFACE qs_uspp_overlap_adjoint
      MODULE PROCEDURE qs_uspp_overlap_single, qs_uspp_overlap_images
   END INTERFACE
   INTERFACE qs_uspp_core_matrix
      MODULE PROCEDURE qs_uspp_core_matrix_single, qs_uspp_core_matrix_images
   END INTERFACE

CONTAINS

! **************************************************************************************************
!> \brief Prepare original global density, valence tau, and frozen-core fields for nonlocal XC.
!>        USPP density uses its original augmentation independently of the kinetic transition.
!>        Neither charge moments nor kinetic matrices determine missing kinetic-density fields.
!> \param uspp Current centers and immutable kind caches; charge and overlap are unchanged.
!> \param qs_kind_set Original scalar UPF datasets and their native projector ordering.
! **************************************************************************************************
   SUBROUTINE qs_uspp_prepare_global_xc(uspp, qs_kind_set)
      TYPE(qs_uspp_type), INTENT(INOUT)                  :: uspp
      TYPE(qs_kind_type), INTENT(IN)                     :: qs_kind_set(:)

      CHARACTER(LEN=256)                                 :: reason
      INTEGER                                            :: ikind
      LOGICAL                                            :: available
      REAL(dp), ALLOCATABLE                              :: tau_core(:)

      IF (uspp%global_xc_ready) RETURN
      CPASSERT(SIZE(qs_kind_set) == SIZE(uspp%augmentation))
      ALLOCATE (uspp%xc_density(SIZE(qs_kind_set)), uspp%xc_kinetic(SIZE(qs_kind_set)), &
                uspp%xc_core_density(SIZE(qs_kind_set)), uspp%xc_core_kinetic(SIZE(qs_kind_set)))
      DO ikind = 1, SIZE(qs_kind_set)
         IF (uspp%augmentation(ikind)%nproj == 0) CYCLE
         ASSOCIATE (potential => qs_kind_set(ikind)%sgp_potential)
            CALL upf_kinetic_reconstruction_init(potential%upf_data, potential%upf_projectors, &
                                                 uspp%xc_kinetic(ikind), available, reason)
            IF (.NOT. available) CALL cp_abort(__LOCATION__, TRIM(reason))
            IF (potential%upf_data%is_paw) THEN
               CALL paw_reconstruction_init(potential%upf_data, potential%upf_projectors, &
                                            uspp%xc_density(ikind), realspace=.TRUE., logarithmic_origin=.TRUE.)
               CALL paw_core_kinetic_build_ae(potential%upf_data, tau_core, available, reason)
               IF (.NOT. available) CALL cp_abort(__LOCATION__, TRIM(reason))
               IF (ANY(potential%upf_data%paw%rho_core /= 0.0_dp)) THEN
                  CALL upf_radial_scalar_init(potential%upf_data, potential%upf_data%paw%rho_core, &
                                              uspp%xc_core_density(ikind), logarithmic_origin=.TRUE.)
               END IF
               IF (ANY(tau_core /= 0.0_dp)) THEN
                  CALL upf_radial_scalar_init(potential%upf_data, tau_core, uspp%xc_core_kinetic(ikind), &
                                              logarithmic_origin=.TRUE.)
               END IF
               DEALLOCATE (tau_core)
            ELSE
               CALL uspp_augmentation_init(potential%upf_data, potential%upf_projectors, uspp%xc_density(ikind))
               CALL uspp_augmentation_prepare_realspace(potential%upf_data, uspp%xc_density(ikind))
               CALL upf_nlcc_kinetic(potential%upf_data, tau_core, available, reason)
               IF (.NOT. available) CALL cp_abort(__LOCATION__, TRIM(reason))
               IF (ALLOCATED(potential%upf_data%rho_nlcc)) THEN
                  IF (ANY(potential%upf_data%rho_nlcc /= 0.0_dp)) THEN
                     CALL upf_radial_scalar_init(potential%upf_data, potential%upf_data%rho_nlcc, &
                                                 uspp%xc_core_density(ikind), logarithmic_origin=.TRUE.)
                  END IF
               END IF
               IF (ANY(tau_core /= 0.0_dp)) THEN
                  CALL upf_radial_scalar_init(potential%upf_data, tau_core, uspp%xc_core_kinetic(ikind), &
                                              logarithmic_origin=.TRUE.)
               END IF
               DEALLOCATE (tau_core)
            END IF
         END ASSOCIATE
      END DO
      uspp%global_xc_ready = .TRUE.
   END SUBROUTINE qs_uspp_prepare_global_xc

! **************************************************************************************************
!> \brief Collocate both original valence primitives before constructing nonlocal XC features.
!> \param uspp Prepared global XC caches and current geometry.
!> \param fields Projected spin density in the original magnetic-projector ordering.
!> \param cell Auxiliary simulation cell.
!> \param points Cartesian grid points.
!> \param rho Spin-density values, incremented.
!> \param gradient Cartesian spin-density gradient, incremented.
!> \param tau Spin kinetic-density values, incremented.
! **************************************************************************************************
   SUBROUTINE qs_uspp_global_atom_grid_collocate(uspp, fields, cell, points, rho, gradient, tau)
      TYPE(qs_uspp_type), INTENT(IN)                     :: uspp
      TYPE(qs_uspp_atom_grid_type), INTENT(IN)           :: fields
      TYPE(cell_type), INTENT(IN)                        :: cell
      REAL(dp), INTENT(IN)                               :: points(:, :)
      REAL(dp), INTENT(INOUT)                            :: rho(:), gradient(:, :), tau(:)

      REAL(dp), ALLOCATABLE                              :: tau_gradient(:, :)

      CPASSERT(uspp%global_xc_ready)
      CALL qs_uspp_atom_grid_collocate(uspp, fields, cell, points, rho, gradient, uspp%xc_density)
      ALLOCATE (tau_gradient(3, SIZE(tau)), SOURCE=0.0_dp)
      CALL qs_uspp_atom_grid_collocate(uspp, fields, cell, points, tau, tau_gradient, uspp%xc_kinetic)
   END SUBROUTINE qs_uspp_global_atom_grid_collocate

! **************************************************************************************************
!> \brief Pull back density, gradient and tau adjoints through the same original valence fields.
!> \param uspp Same immutable caches and geometry as in collocation.
!> \param fields AO operators and source-geometry adjoints, incremented.
!> \param cell Same auxiliary cell.
!> \param points Same Cartesian grid points.
!> \param weights Local quadrature weights.
!> \param v_rho Derivative with respect to spin density.
!> \param v_gradient Derivative with respect to its Cartesian gradient.
!> \param v_tau Derivative with respect to spin kinetic density.
!> \param point_gradient Optional moving-grid coordinate derivative, incremented for forces.
! **************************************************************************************************
   SUBROUTINE qs_uspp_global_atom_grid_integrate(uspp, fields, cell, points, weights, &
                                                 v_rho, v_gradient, v_tau, point_gradient)
      TYPE(qs_uspp_type), INTENT(IN)                     :: uspp
      TYPE(qs_uspp_atom_grid_type), INTENT(INOUT)        :: fields
      TYPE(cell_type), INTENT(IN)                        :: cell
      REAL(dp), INTENT(IN)                               :: points(:, :), weights(:), v_rho(:), &
                                                            v_gradient(:, :), v_tau(:)
      REAL(dp), INTENT(INOUT), OPTIONAL                  :: point_gradient(:, :)

      REAL(dp), ALLOCATABLE                              :: zero_gradient(:, :)

      CPASSERT(uspp%global_xc_ready)
      CALL qs_uspp_atom_grid_integrate(uspp, fields, cell, points, weights, v_rho, v_gradient, &
                                       point_gradient, uspp%xc_density)
      ALLOCATE (zero_gradient(3, SIZE(weights)), SOURCE=0.0_dp)
      CALL qs_uspp_atom_grid_integrate(uspp, fields, cell, points, weights, v_tau, zero_gradient, &
                                       point_gradient, uspp%xc_kinetic)
   END SUBROUTINE qs_uspp_global_atom_grid_integrate

! **************************************************************************************************
!> \brief Add the independent semilocal USPP kinetic transition on the same periodic atom grid.
!> \param uspp Original data and current geometry, with an independent kinetic transition.
!> \param fields Spin projected density and its adjoint accumulators.
!> \param cell Auxiliary periodic cell.
!> \param points Cartesian quadrature points.
!> \param tau Spin kinetic density, incremented without weights.
! **************************************************************************************************
   SUBROUTINE qs_uspp_kinetic_atom_grid_collocate(uspp, fields, cell, points, tau)
      TYPE(qs_uspp_type), INTENT(IN)                     :: uspp
      TYPE(qs_uspp_atom_grid_type), INTENT(IN)           :: fields
      TYPE(cell_type), INTENT(IN)                        :: cell
      REAL(dp), INTENT(IN)                               :: points(:, :)
      REAL(dp), INTENT(INOUT)                            :: tau(:)

      REAL(dp), ALLOCATABLE                              :: gradient(:, :)

      IF (.NOT. ALLOCATED(uspp%kinetic_augmentation)) RETURN
      ALLOCATE (gradient(3, SIZE(tau)), SOURCE=0.0_dp)
      CALL qs_uspp_atom_grid_collocate(uspp, fields, cell, points, tau, gradient, uspp%kinetic_augmentation)
   END SUBROUTINE qs_uspp_kinetic_atom_grid_collocate

! **************************************************************************************************
!> \brief Pull back the semilocal USPP tau covector, including source, target and image strain.
!> \param uspp Same independent kinetic cache as in collocation.
!> \param fields Spin projected density and adjoint accumulators.
!> \param cell Auxiliary periodic cell.
!> \param points Cartesian quadrature points.
!> \param weights Quadrature weights; one when already included in potential.
!> \param potential Derivative with respect to spin kinetic density.
!> \param point_gradient Optional target gradients for the force path.
! **************************************************************************************************
   SUBROUTINE qs_uspp_kinetic_atom_grid_integrate(uspp, fields, cell, points, weights, potential, point_gradient)
      TYPE(qs_uspp_type), INTENT(IN)                     :: uspp
      TYPE(qs_uspp_atom_grid_type), INTENT(INOUT)        :: fields
      TYPE(cell_type), INTENT(IN)                        :: cell
      REAL(dp), INTENT(IN)                               :: points(:, :), weights(:), potential(:)
      REAL(dp), INTENT(INOUT), OPTIONAL                  :: point_gradient(:, :)

      REAL(dp), ALLOCATABLE                              :: zero_gradient(:, :)

      IF (.NOT. ALLOCATED(uspp%kinetic_augmentation)) RETURN
      ALLOCATE (zero_gradient(3, SIZE(potential)), SOURCE=0.0_dp)
      CALL qs_uspp_atom_grid_integrate(uspp, fields, cell, points, weights, potential, zero_gradient, &
                                       point_gradient, uspp%kinetic_augmentation)
   END SUBROUTINE qs_uspp_kinetic_atom_grid_integrate

! **************************************************************************************************
!> \brief Select the immutable core field for a scalar USPP kind.
!> \param uspp Current per-kind caches; the returned pointer is used only within this geometry state.
!> \param ikind Kind index.
!> \return Associated cache for an explicitly supplied UPF model core.
! **************************************************************************************************
   FUNCTION core_tau_cache(uspp, ikind) RESULT(cache)
      TYPE(qs_uspp_type), INTENT(IN), TARGET             :: uspp
      INTEGER, INTENT(IN)                                :: ikind
      TYPE(upf_core_kinetic_type), POINTER               :: cache

      NULLIFY (cache)
      IF (ALLOCATED(uspp%model_core_kinetic)) THEN
         cache => uspp%model_core_kinetic(ikind)
      END IF
   END FUNCTION core_tau_cache

! **************************************************************************************************
!> \brief Whether a nonzero verified UPF model-core kinetic field is present.
!> \param ks_env Current USPP caches.
!> \return True for at least one nonzero auxiliary PS core field.
! **************************************************************************************************
   FUNCTION qs_upf_core_tau_present(ks_env) RESULT(present_core)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      LOGICAL                                            :: present_core

      INTEGER                                            :: ikind
      TYPE(qs_uspp_type), POINTER                        :: uspp
      TYPE(upf_core_kinetic_type), POINTER               :: cache

      present_core = .FALSE.
      CALL get_ks_env(ks_env, uspp=uspp)
      IF (.NOT. ASSOCIATED(uspp)) RETURN
      DO ikind = 1, SIZE(uspp%augmentation)
         cache => core_tau_cache(uspp, ikind)
         IF (.NOT. ASSOCIATED(cache)) CYCLE
         IF (.NOT. cache%available) CYCLE
         IF (ANY(cache%weighted /= 0.0_dp)) present_core = .TRUE.
      END DO
   END FUNCTION qs_upf_core_tau_present

! **************************************************************************************************
!> \brief Spin-summed periodic PS core kinetic field on the requested XC grid.
!>        The R-to-G adjoint of the real projection is used by the force routine.
!> \param ks_env Current USPP centers and verified radial caches.
!> \param field Allocated real PW field, overwritten without changing any valence density.
! **************************************************************************************************
   SUBROUTINE qs_upf_core_tau_field(ks_env, field)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: field

      COMPLEX(dp)                                        :: coefficient
      INTEGER                                            :: center, ig, ikind
      TYPE(pw_c1d_gs_type)                               :: core_g
      TYPE(qs_uspp_type), POINTER                        :: uspp
      TYPE(upf_core_kinetic_type), POINTER               :: cache

      CALL get_ks_env(ks_env, uspp=uspp)
      CPASSERT(ASSOCIATED(uspp))
      CALL check_grid(uspp, field%pw_grid)
      CALL core_g%create(field%pw_grid)
      core_g%array = CMPLX(0.0_dp, 0.0_dp, dp)
      DO center = 1, SIZE(uspp%center_atoms)
         ikind = uspp%center_kinds(center)
         cache => core_tau_cache(uspp, ikind)
         IF (.NOT. ASSOCIATED(cache)) CYCLE
         IF (.NOT. cache%available) CYCLE
         IF (.NOT. ANY(cache%weighted /= 0.0_dp)) CYCLE
         DO ig = 1, SIZE(core_g%array)
            CALL upf_core_kinetic_fourier(cache, field%pw_grid%g(:, ig), uspp%position(:, center), &
                                          field%pw_grid%vol, coefficient)
            core_g%array(ig) = core_g%array(ig) + coefficient
         END DO
      END DO
      CALL pw_transfer(core_g, field)
      CALL core_g%release()
   END SUBROUTINE qs_upf_core_tau_field

! **************************************************************************************************
!> \brief Center and affine-strain adjoints of the smooth PW core kinetic field.
!>        This fixed core field has no AO density-matrix derivative.
!> \param ks_env Current centers and unreduced force/virial accumulators.
!> \param potential Unweighted derivative of XC energy density with respect to spin tau.
!> \param spin_fraction Fraction of the spin-summed closed core in this channel.
! **************************************************************************************************
   SUBROUTINE qs_upf_core_tau_forces(ks_env, potential, spin_fraction)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: potential
      REAL(dp), INTENT(IN)                               :: spin_fraction

      COMPLEX(dp)                                        :: coefficient, derivative(3), &
                                                            strain(3, 3), weighted
      INTEGER                                            :: center, ig, ikind, natom
      REAL(dp)                                           :: multiplicity
      REAL(dp), ALLOCATABLE                              :: atom_gradient(:, :), cell_gradient(:, :)
      TYPE(pw_c1d_gs_type)                               :: potential_g
      TYPE(qs_uspp_type), POINTER                        :: uspp
      TYPE(upf_core_kinetic_type), POINTER               :: cache

      IF (.NOT. qs_upf_core_tau_present(ks_env)) RETURN
      CALL get_ks_env(ks_env, uspp=uspp, natom=natom)
      CALL check_grid(uspp, potential%pw_grid)
      ALLOCATE (atom_gradient(3, natom), SOURCE=0.0_dp)
      CALL allocate_cell_gradient(ks_env, cell_gradient)
      CALL potential_g%create(potential%pw_grid)
      CALL pw_transfer(potential, potential_g)
      DO center = 1, SIZE(uspp%center_atoms)
         ikind = uspp%center_kinds(center)
         cache => core_tau_cache(uspp, ikind)
         IF (.NOT. ASSOCIATED(cache)) CYCLE
         IF (.NOT. cache%available) CYCLE
         IF (.NOT. ANY(cache%weighted /= 0.0_dp)) CYCLE
         DO ig = 1, SIZE(potential_g%array)
            multiplicity = 1.0_dp
            IF (potential%pw_grid%grid_span == HALFSPACE) THEN
               multiplicity = 2.0_dp
               IF (potential%pw_grid%have_g0 .AND. ig == 1) multiplicity = 1.0_dp
            END IF
            weighted = spin_fraction*potential%pw_grid%vol*multiplicity*potential_g%array(ig)
            IF (ALLOCATED(cell_gradient)) THEN
               CALL upf_core_kinetic_fourier(cache, potential%pw_grid%g(:, ig), uspp%position(:, center), &
                                             potential%pw_grid%vol, coefficient, strain, derivative)
               cell_gradient(:, :) = cell_gradient + REAL(CONJG(strain)*weighted, dp)
            ELSE
               CALL upf_core_kinetic_fourier(cache, potential%pw_grid%g(:, ig), uspp%position(:, center), &
                                             potential%pw_grid%vol, coefficient, center_derivative=derivative)
            END IF
            atom_gradient(:, uspp%center_atoms(center)) = atom_gradient(:, uspp%center_atoms(center)) + &
                                                          REAL(CONJG(derivative)*weighted, dp)
         END DO
      END DO
      CALL accumulate_gradient(ks_env, atom_gradient, cell_gradient)
      CALL potential_g%release()
   END SUBROUTINE qs_upf_core_tau_forces

! **************************************************************************************************
!> \brief Add spin-summed PS core tau from all centers and auxiliary-cell images.
!> \param uspp Current USPP centers and original-data kinetic caches.
!> \param cell Auxiliary periodic cell, as in augmentation collocation.
!> \param points Unfolded Cartesian atom-grid points.
!> \param kinetic Spin-summed core tau values, incremented without quadrature weights.
! **************************************************************************************************
   SUBROUTINE qs_upf_core_tau_atom_collocate(uspp, cell, points, kinetic)
      TYPE(qs_uspp_type), INTENT(IN), TARGET             :: uspp
      TYPE(cell_type), INTENT(IN)                        :: cell
      REAL(dp), INTENT(IN)                               :: points(:, :)
      REAL(dp), INTENT(INOUT)                            :: kinetic(:)

      INTEGER                                            :: center, ikind, ix, iy, iz, lower(3), &
                                                            point, upper(3)
      REAL(dp)                                           :: image_position(3), value
      REAL(dp), ALLOCATABLE                              :: folded(:, :)
      TYPE(upf_core_kinetic_type), POINTER               :: cache

      CPASSERT(ALL(SHAPE(points) == [3, SIZE(kinetic)]))
      IF (SIZE(kinetic) == 0) RETURN
      ALLOCATE (folded, SOURCE=points - MATMUL(cell%hmat, ANINT(MATMUL(cell%h_inv, points))))
      DO center = 1, SIZE(uspp%center_atoms)
         ikind = uspp%center_kinds(center)
         cache => core_tau_cache(uspp, ikind)
         IF (.NOT. ASSOCIATED(cache)) CYCLE
         IF (.NOT. cache%available) CYCLE
         IF (.NOT. ANY(cache%weighted /= 0.0_dp)) CYCLE
         CALL atom_grid_image_bounds(cell, uspp%position(:, center), folded, cache%support_radius, lower, upper)
         DO iz = lower(3), upper(3)
            DO iy = lower(2), upper(2)
               DO ix = lower(1), upper(1)
                  image_position = uspp%position(:, center) + MATMUL(cell%hmat, REAL([ix, iy, iz], dp))
                  DO point = 1, SIZE(kinetic)
                     CALL upf_core_kinetic_cartesian(cache, folded(:, point) - image_position, value)
                     kinetic(point) = kinetic(point) + value
                  END DO
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE qs_upf_core_tau_atom_collocate

! **************************************************************************************************
!> \brief Exact center and point adjoints of periodic atom-grid PS core tau collocation.
!> \param uspp Same USPP state as the forward map.
!> \param fields Accumulators for source-center and source-image strain derivatives.
!> \param cell Auxiliary cell used by the forward map.
!> \param points Same unfolded quadrature points.
!> \param potential Weighted derivative with respect to spin-summed core tau, including spin fractions.
!> \param point_gradient Weighted point energy gradients, incremented; caller handles target-grid motion.
! **************************************************************************************************
   SUBROUTINE qs_upf_core_tau_atom_integrate(uspp, fields, cell, points, potential, point_gradient)
      TYPE(qs_uspp_type), INTENT(IN), TARGET             :: uspp
      TYPE(qs_uspp_atom_grid_type), INTENT(INOUT)        :: fields
      TYPE(cell_type), INTENT(IN)                        :: cell
      REAL(dp), INTENT(IN)                               :: points(:, :), potential(:)
      REAL(dp), INTENT(INOUT)                            :: point_gradient(:, :)

      INTEGER                                            :: center, ikind, ix, iy, iz, jdir, &
                                                            lower(3), point, upper(3)
      REAL(dp)                                           :: derivative(3), image_position(3), &
                                                            spatial(3), value
      REAL(dp), ALLOCATABLE                              :: folded(:, :), images(:, :), &
                                                            point_derivative(:, :)
      TYPE(upf_core_kinetic_type), POINTER               :: cache

      CPASSERT(fields%calculate_forces)
      CPASSERT(ALL(SHAPE(points) == [3, SIZE(potential)]))
      CPASSERT(ALL(SHAPE(point_gradient) == SHAPE(points)))
      IF (SIZE(potential) == 0) RETURN
      ALLOCATE (images, SOURCE=MATMUL(cell%hmat, ANINT(MATMUL(cell%h_inv, points))))
      ALLOCATE (folded, SOURCE=points - images)
      ALLOCATE (point_derivative(3, SIZE(potential)))
      DO center = 1, SIZE(uspp%center_atoms)
         ikind = uspp%center_kinds(center)
         cache => core_tau_cache(uspp, ikind)
         IF (.NOT. ASSOCIATED(cache)) CYCLE
         IF (.NOT. cache%available) CYCLE
         IF (.NOT. ANY(cache%weighted /= 0.0_dp)) CYCLE
         CALL atom_grid_image_bounds(cell, uspp%position(:, center), folded, cache%support_radius, lower, upper)
         DO iz = lower(3), upper(3)
            DO iy = lower(2), upper(2)
               DO ix = lower(1), upper(1)
                  image_position = uspp%position(:, center) + MATMUL(cell%hmat, REAL([ix, iy, iz], dp))
                  DO point = 1, SIZE(potential)
                     CALL upf_core_kinetic_cartesian(cache, folded(:, point) - image_position, value, spatial)
                     point_derivative(:, point) = potential(point)*spatial
                  END DO
                  derivative = -SUM(point_derivative, DIM=2)
                  fields%center_gradient(:, center) = fields%center_gradient(:, center) + derivative
                  point_gradient = point_gradient + point_derivative
                  fields%center_strain = fields%center_strain - MATMUL(point_derivative, TRANSPOSE(images))
                  DO jdir = 1, 3
                     fields%center_strain(:, jdir) = fields%center_strain(:, jdir) + derivative*image_position(jdir)
                  END DO
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE qs_upf_core_tau_atom_integrate

! **************************************************************************************************
!> \brief Prepare replicated on-site densities once per spin for distributed atom-grid batches.
!> \param uspp Geometry state, whose real-space radial caches are initialized lazily.
!> \param qs_kind_set Kinds retaining the original UPF radial data.
!> \param density Spin AO density over real-space images.
!> \param fields Owned batch state; previous allocations and accumulated operators are discarded.
!> \param calculate_forces Whether to accumulate center, point and source-center strain derivatives.
! **************************************************************************************************
   SUBROUTINE qs_uspp_atom_grid_prepare(uspp, qs_kind_set, density, fields, calculate_forces)
      TYPE(qs_uspp_type), INTENT(INOUT)                  :: uspp
      TYPE(qs_kind_type), INTENT(IN)                     :: qs_kind_set(:)
      TYPE(dbcsr_p_type), INTENT(IN)                     :: density(:)
      TYPE(qs_uspp_atom_grid_type), INTENT(INOUT)        :: fields
      LOGICAL, INTENT(IN)                                :: calculate_forces

      INTEGER                                            :: center, ikind, n, nmax
      LOGICAL                                            :: found
      REAL(KIND=dp), POINTER                             :: block(:, :)
      TYPE(dbcsr_type)                                   :: projected
      TYPE(mp_comm_type)                                 :: group

      IF (ALLOCATED(fields%density)) DEALLOCATE (fields%density)
      IF (ALLOCATED(fields%potential)) DEALLOCATE (fields%potential)
      IF (ALLOCATED(fields%center_gradient)) DEALLOCATE (fields%center_gradient)
      nmax = MAXVAL(uspp%projector_sizes)
      ALLOCATE (fields%density(nmax, nmax, SIZE(uspp%center_atoms)), &
                fields%potential(nmax, nmax, SIZE(uspp%center_atoms)), &
                fields%center_gradient(3, SIZE(uspp%center_atoms)))
      fields%calculate_forces = calculate_forces
      fields%density = 0.0_dp
      fields%potential = 0.0_dp
      fields%center_gradient = 0.0_dp
      fields%center_strain = 0.0_dp
      CALL dbcsr_get_info(uspp%overlaps, group=group)
      CALL project_density(uspp, density, projected)
      DO center = 1, SIZE(uspp%center_atoms)
         ikind = uspp%center_kinds(center)
         n = uspp%projector_sizes(center)
         CALL uspp_augmentation_prepare_realspace(qs_kind_set(ikind)%sgp_potential%upf_data, uspp%augmentation(ikind))
         CALL dbcsr_get_block_p(projected, center, center, block, found)
         IF (found) fields%density(:n, :n, center) = block
      END DO
      CALL group%sum(fields%density)
      CALL dbcsr_release(projected)
   END SUBROUTINE qs_uspp_atom_grid_prepare

! **************************************************************************************************
!> \brief Fixed unit-amplitude scalar fields reuse the atom-grid image and derivative transport.
!>        Only kinds with scalar caches contribute. They have no AO or overlap dependence.
!> \param uspp Current source-center geometry; ghost centers are already excluded.
!> \param radial_fields Original scalar caches, with nproj=1 for contributing kinds.
!> \param fixed Independent geometry and unit amplitudes, replaced on entry.
!> \param calculate_forces Accumulate center, target and source-image strain derivatives.
! **************************************************************************************************
   SUBROUTINE qs_uspp_fixed_atom_grid_prepare(uspp, radial_fields, fixed, calculate_forces)
      TYPE(qs_uspp_type), INTENT(IN)                     :: uspp
      TYPE(uspp_augmentation_type), INTENT(IN)           :: radial_fields(:)
      TYPE(qs_uspp_fixed_atom_grid_type), INTENT(OUT)    :: fixed
      LOGICAL, INTENT(IN)                                :: calculate_forces

      INTEGER                                            :: center, count, kind

      CPASSERT(SIZE(radial_fields) == SIZE(uspp%augmentation))
      count = 0
      DO center = 1, SIZE(uspp%center_atoms)
         kind = uspp%center_kinds(center)
         CPASSERT(radial_fields(kind)%nproj == 0 .OR. radial_fields(kind)%nproj == 1)
         IF (radial_fields(kind)%nproj == 1) count = count + 1
      END DO
      ALLOCATE (fixed%source_centers(count))
      count = 0
      DO center = 1, SIZE(uspp%center_atoms)
         IF (radial_fields(uspp%center_kinds(center))%nproj == 0) CYCLE
         count = count + 1
         fixed%source_centers(count) = center
      END DO
      ALLOCATE (fixed%geometry%center_atoms(count), fixed%geometry%center_kinds(count), &
                fixed%geometry%position(3, count))
      fixed%geometry%center_atoms(:) = uspp%center_atoms(fixed%source_centers)
      fixed%geometry%center_kinds(:) = uspp%center_kinds(fixed%source_centers)
      fixed%geometry%position(:, :) = uspp%position(:, fixed%source_centers)
      ALLOCATE (fixed%geometry%projector_sizes(count), SOURCE=1)
      ! The independent caches are supplied to both forward and adjoint calls.
      ALLOCATE (fixed%geometry%augmentation(SIZE(radial_fields)))
      ALLOCATE (fixed%fields%density(1, 1, count), SOURCE=1.0_dp)
      ALLOCATE (fixed%fields%potential(1, 1, count), fixed%fields%center_gradient(3, count), SOURCE=0.0_dp)
      fixed%fields%calculate_forces = calculate_forces
   END SUBROUTINE qs_uspp_fixed_atom_grid_prepare

! **************************************************************************************************
!> \brief Fixed-field adjoints add geometry derivatives to the ordinary source accumulators.
!>        Discard the scalar amplitude adjoint: frozen cores have no AO density dependence.
!> \param fixed Same fixed geometry and unit amplitudes as in collocation.
!> \param radial_fields Same original scalar caches as in collocation.
!> \param fields Ordinary source-center accumulators; density and AO potential remain unchanged.
!> \param cell Same auxiliary cell as in collocation.
!> \param points Same unfolded atom-grid points.
!> \param weights Local quadrature weights.
!> \param v_rho Derivative with respect to the scalar field, including spin fractions.
!> \param v_gradient Derivative with respect to its Cartesian gradient.
!> \param point_gradient Weighted target-coordinate derivatives, incremented when forces are requested.
! **************************************************************************************************
   SUBROUTINE qs_uspp_fixed_atom_grid_integrate(fixed, radial_fields, fields, cell, points, &
                                                weights, v_rho, v_gradient, point_gradient)
      TYPE(qs_uspp_fixed_atom_grid_type), INTENT(INOUT)  :: fixed
      TYPE(uspp_augmentation_type), INTENT(IN), TARGET   :: radial_fields(:)
      TYPE(qs_uspp_atom_grid_type), INTENT(INOUT)        :: fields
      TYPE(cell_type), INTENT(IN)                        :: cell
      REAL(dp), INTENT(IN)                               :: points(:, :), weights(:), v_rho(:), &
                                                            v_gradient(:, :)
      REAL(dp), INTENT(INOUT), OPTIONAL                  :: point_gradient(:, :)

      CPASSERT(fields%calculate_forces .EQV. fixed%fields%calculate_forces)
      fixed%fields%potential = 0.0_dp
      fixed%fields%center_gradient = 0.0_dp
      fixed%fields%center_strain = 0.0_dp
      CALL qs_uspp_atom_grid_integrate(fixed%geometry, fixed%fields, cell, points, weights, &
                                       v_rho, v_gradient, point_gradient, radial_fields)
      IF (fields%calculate_forces) THEN
         fields%center_gradient(:, fixed%source_centers) = &
            fields%center_gradient(:, fixed%source_centers) + fixed%fields%center_gradient
         fields%center_strain = fields%center_strain + fixed%fields%center_strain
      END IF
   END SUBROUTINE qs_uspp_fixed_atom_grid_integrate

! **************************************************************************************************
!> \brief Add augmentation from all centers and their periodic images to local quadrature points.
!>        The auxiliary density is periodic in all three cell directions, as on Quickstep's PW grid.
!>        Different ranks may supply different point counts, including zero; no collective occurs here.
!> \param uspp Geometry state used to prepare fields; must not change during a batch evaluation.
!> \param fields Replicated spin density and batch accumulators.
!> \param cell Auxiliary simulation cell, including its inverse matrix.
!> \param points Cartesian quadrature points in bohr, with arbitrary cell representatives.
!> \param rho Density values, incremented without quadrature weights.
!> \param gradient Cartesian density gradients, incremented without quadrature weights.
!> \param radial_fields Optional independent per-kind transition caches, such as valence kinetic density.
!>        They share the native projector ordering; charge augmentation and the metric remain unchanged.
! **************************************************************************************************
   SUBROUTINE qs_uspp_atom_grid_collocate(uspp, fields, cell, points, rho, gradient, radial_fields)
      TYPE(qs_uspp_type), INTENT(IN), TARGET             :: uspp
      TYPE(qs_uspp_atom_grid_type), INTENT(IN)           :: fields
      TYPE(cell_type), INTENT(IN)                        :: cell
      REAL(KIND=dp), INTENT(IN)                          :: points(:, :)
      REAL(KIND=dp), INTENT(INOUT)                       :: rho(:), gradient(:, :)
      TYPE(uspp_augmentation_type), INTENT(IN), &
         OPTIONAL, TARGET                                :: radial_fields(:)

      INTEGER                                            :: center, ikind, ix, iy, iz, lower(3), n, &
                                                            upper(3)
      REAL(KIND=dp)                                      :: image_position(3)
      REAL(KIND=dp), ALLOCATABLE                         :: displacements(:, :), &
                                                            folded_points(:, :), point_images(:, :)
      TYPE(uspp_augmentation_type), POINTER              :: cache

      CPASSERT(ALL(SHAPE(points) == [3, SIZE(rho)]))
      CPASSERT(ALL(SHAPE(gradient) == SHAPE(points)))
      CPASSERT(SIZE(fields%density, 3) == SIZE(uspp%center_atoms))
      IF (PRESENT(radial_fields)) THEN
         CPASSERT(SIZE(radial_fields) == SIZE(uspp%augmentation))
      END IF
      IF (SIZE(rho) == 0) RETURN
      ALLOCATE (displacements(3, SIZE(rho)))
      ALLOCATE (folded_points(3, SIZE(points, 2)), point_images(3, SIZE(points, 2)))
      point_images(:, :) = MATMUL(cell%hmat, ANINT(MATMUL(cell%h_inv, points)))
      folded_points(:, :) = points - point_images
      DO center = 1, SIZE(uspp%center_atoms)
         ikind = uspp%center_kinds(center)
         n = uspp%projector_sizes(center)
         cache => uspp%augmentation(ikind)
         IF (PRESENT(radial_fields)) cache => radial_fields(ikind)
         IF (cache%nproj == 0) CYCLE
         CPASSERT(cache%nproj == n)
         CALL atom_grid_image_bounds(cell, uspp%position(:, center), folded_points, &
                                     cache%realspace_radius, lower, upper)
         DO iz = lower(3), upper(3)
            DO iy = lower(2), upper(2)
               DO ix = lower(1), upper(1)
                  image_position = uspp%position(:, center) + MATMUL(cell%hmat, REAL([ix, iy, iz], dp))
                  displacements(:, :) = folded_points - SPREAD(image_position, 2, SIZE(rho))
                  CALL uspp_augmentation_collocate_atom(cache, displacements, &
                                                        fields%density(:n, :n, center), rho, gradient)
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE qs_uspp_atom_grid_collocate

! **************************************************************************************************
!> \brief Accumulate the adjoint of all-center augmentation collocation for one local grid batch.
!>        The same image sum is differentiated, including overlapping augmentation regions.
!> \param uspp Geometry state used in collocation.
!> \param fields On-site operators and source-center derivatives, incremented locally without MPI reduction.
!> \param cell Auxiliary simulation cell used in collocation.
!> \param points Cartesian quadrature points in bohr.
!> \param weights Local quadrature weights, including any partition or hard/soft sign.
!> \param v_rho Unweighted energy-density derivative with respect to this spin density.
!> \param v_gradient Unweighted derivative with respect to the Cartesian spin-density gradient.
!> \param point_gradient Weighted energy derivatives with respect to point coordinates, incremented.
!>        Required when calculate_forces was selected. Moving-grid and weight derivatives are the caller's task.
!> \param radial_fields Same optional transition caches as in the forward collocation.
! **************************************************************************************************
   SUBROUTINE qs_uspp_atom_grid_integrate(uspp, fields, cell, points, weights, v_rho, v_gradient, &
                                          point_gradient, radial_fields)
      TYPE(qs_uspp_type), INTENT(IN), TARGET             :: uspp
      TYPE(qs_uspp_atom_grid_type), INTENT(INOUT)        :: fields
      TYPE(cell_type), INTENT(IN)                        :: cell
      REAL(KIND=dp), INTENT(IN)                          :: points(:, :), weights(:), v_rho(:), &
                                                            v_gradient(:, :)
      REAL(KIND=dp), INTENT(INOUT), OPTIONAL             :: point_gradient(:, :)
      TYPE(uspp_augmentation_type), INTENT(IN), &
         OPTIONAL, TARGET                                :: radial_fields(:)

      INTEGER                                            :: center, ikind, ix, iy, iz, jdir, &
                                                            lower(3), n, upper(3)
      REAL(KIND=dp)                                      :: derivative(3), image_position(3)
      REAL(KIND=dp), ALLOCATABLE                         :: displacements(:, :), &
                                                            folded_points(:, :), matrix(:, :), &
                                                            point_derivative(:, :), &
                                                            point_images(:, :)
      TYPE(uspp_augmentation_type), POINTER              :: cache

      CPASSERT(ALL(SHAPE(points) == [3, SIZE(weights)]))
      CPASSERT(SIZE(v_rho) == SIZE(weights))
      CPASSERT(ALL(SHAPE(v_gradient) == SHAPE(points)))
      CPASSERT(PRESENT(point_gradient) .EQV. fields%calculate_forces)
      IF (PRESENT(radial_fields)) THEN
         CPASSERT(SIZE(radial_fields) == SIZE(uspp%augmentation))
      END IF
      IF (PRESENT(point_gradient)) THEN
         CPASSERT(ALL(SHAPE(point_gradient) == SHAPE(points)))
      END IF
      IF (SIZE(weights) == 0) RETURN
      ALLOCATE (displacements(3, SIZE(weights)))
      IF (fields%calculate_forces) ALLOCATE (point_derivative(3, SIZE(weights)))
      ALLOCATE (folded_points(3, SIZE(points, 2)), point_images(3, SIZE(points, 2)))
      point_images(:, :) = MATMUL(cell%hmat, ANINT(MATMUL(cell%h_inv, points)))
      folded_points(:, :) = points - point_images
      DO center = 1, SIZE(uspp%center_atoms)
         ikind = uspp%center_kinds(center)
         n = uspp%projector_sizes(center)
         cache => uspp%augmentation(ikind)
         IF (PRESENT(radial_fields)) cache => radial_fields(ikind)
         IF (cache%nproj == 0) CYCLE
         CPASSERT(cache%nproj == n)
         ALLOCATE (matrix(n, n))
         CALL atom_grid_image_bounds(cell, uspp%position(:, center), folded_points, &
                                     cache%realspace_radius, lower, upper)
         DO iz = lower(3), upper(3)
            DO iy = lower(2), upper(2)
               DO ix = lower(1), upper(1)
                  image_position = uspp%position(:, center) + MATMUL(cell%hmat, REAL([ix, iy, iz], dp))
                  displacements(:, :) = folded_points - SPREAD(image_position, 2, SIZE(weights))
                  IF (fields%calculate_forces) THEN
                     CALL uspp_augmentation_integrate_atom(cache, displacements, weights, &
                                                           v_rho, v_gradient, matrix, fields%density(:n, :n, center), &
                                                           derivative, point_derivative)
                     fields%center_gradient(:, center) = fields%center_gradient(:, center) + derivative
                     point_gradient = point_gradient + point_derivative
                     ! Folding each target changes the image label of its source center.
                     ! Restore that source-image strain at the original, unfolded point.
                     fields%center_strain = fields%center_strain - MATMUL(point_derivative, TRANSPOSE(point_images))
                     ! Source centers and lattice images follow a homogeneous strain. The
                     ! target atom-grid center and weight contributions are supplied by the XC caller.
                     DO jdir = 1, 3
                        fields%center_strain(:, jdir) = fields%center_strain(:, jdir) + derivative*image_position(jdir)
                     END DO
                  ELSE
                     CALL uspp_augmentation_integrate_atom(cache, displacements, weights, &
                                                           v_rho, v_gradient, matrix)
                  END IF
                  fields%potential(:n, :n, center) = fields%potential(:n, :n, center) + matrix
               END DO
            END DO
         END DO
         DEALLOCATE (matrix)
      END DO
   END SUBROUTINE qs_uspp_atom_grid_integrate

! **************************************************************************************************
!> \brief Reduce accumulated atom-grid operators and apply the complete AO/projector adjoint once.
!> \param uspp Current geometry state, with Gamma or image-resolved projector overlaps.
!> \param fields Local batch accumulators; preserved so the caller can inspect the quadrature contractions.
!> \param hamiltonian AO image matrices to increment.
!> \param density Optional spin AO image densities, required for the force path.
!> \param sap_int Optional projector derivatives, required for the force path.
!> \param atom_gradient Optional atomic energy gradients, incremented locally before the usual force reduction.
!> \param cell_gradient Optional strain energy gradient. Includes source-center and projector terms only;
!>        the caller must add moving-target-grid and quadrature-weight strain derivatives.
! **************************************************************************************************
   SUBROUTINE qs_uspp_atom_grid_apply(uspp, fields, hamiltonian, density, sap_int, atom_gradient, cell_gradient)
      TYPE(qs_uspp_type), INTENT(IN)                     :: uspp
      TYPE(qs_uspp_atom_grid_type), INTENT(IN)           :: fields
      TYPE(dbcsr_p_type), INTENT(INOUT)                  :: hamiltonian(:)
      TYPE(dbcsr_p_type), INTENT(IN), OPTIONAL           :: density(:)
      TYPE(sap_int_type), INTENT(IN), OPTIONAL           :: sap_int(:)
      REAL(KIND=dp), INTENT(INOUT), OPTIONAL             :: atom_gradient(:, :), cell_gradient(3, 3)

      INTEGER                                            :: center, n
      LOGICAL                                            :: found
      REAL(KIND=dp), ALLOCATABLE                         :: matrix(:, :)
      REAL(KIND=dp), POINTER                             :: block(:, :)
      TYPE(dbcsr_type)                                   :: onsite
      TYPE(mp_comm_type)                                 :: group

      CPASSERT(PRESENT(density) .EQV. fields%calculate_forces)
      CPASSERT(PRESENT(atom_gradient) .EQV. fields%calculate_forces)
      CPASSERT(PRESENT(sap_int) .EQV. fields%calculate_forces)
      IF (PRESENT(cell_gradient)) THEN
         CPASSERT(fields%calculate_forces)
      END IF
      CALL dbcsr_get_info(uspp%overlaps, group=group)
      CALL uspp_projector_create_onsite(uspp%overlaps, onsite)
      DO center = 1, SIZE(uspp%center_atoms)
         n = uspp%projector_sizes(center)
         ALLOCATE (matrix(n, n))
         matrix(:, :) = fields%potential(:n, :n, center)
         CALL group%sum(matrix)
         CALL dbcsr_get_block_p(onsite, center, center, block, found)
         IF (found) block = matrix
         IF (PRESENT(atom_gradient)) THEN
            atom_gradient(:, uspp%center_atoms(center)) = atom_gradient(:, uspp%center_atoms(center)) &
                                                          + fields%center_gradient(:, center)
         END IF
         DEALLOCATE (matrix)
      END DO
      IF (PRESENT(cell_gradient)) cell_gradient = cell_gradient + fields%center_strain
      CALL apply_operator(uspp, onsite, hamiltonian, density, sap_int, atom_gradient, cell_gradient)
      CALL dbcsr_release(onsite)
   END SUBROUTINE qs_uspp_atom_grid_apply

! **************************************************************************************************
!> \brief Apply atom-grid XC operators after the KS matrices have been initialized from Hcore.
!> \param ks_env KS environment and force accumulators
!> \param fields Local XC accumulators for each spin
!> \param hamiltonian AO matrices indexed by spin and image
!> \param density Corresponding AO densities for projector derivatives
! **************************************************************************************************
   SUBROUTINE qs_uspp_atom_grid_apply_ks(ks_env, fields, hamiltonian, density)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(qs_uspp_atom_grid_type), INTENT(IN)           :: fields(:)
      TYPE(dbcsr_p_type), INTENT(INOUT)                  :: hamiltonian(:, :)
      TYPE(dbcsr_p_type), INTENT(IN)                     :: density(:, :)

      TYPE(qs_uspp_type), POINTER                        :: uspp

      CALL get_ks_env(ks_env, uspp=uspp)
      CALL qs_uspp_atom_grid_apply_state(ks_env, uspp, fields, hamiltonian, density)
   END SUBROUTINE qs_uspp_atom_grid_apply_ks

! **************************************************************************************************
!> \brief Apply original atom-grid XC terms for an explicit ORB or AUX state.
!> \param ks_env ...
!> \param uspp ...
!> \param fields ...
!> \param hamiltonian ...
!> \param density ...
!> \param basis_type ...
!> \param projector_neighbors ...
!> \param force_factors Per-spin ADMM geometry scaling; AO Hamiltonians are unscaled.
! **************************************************************************************************
   SUBROUTINE qs_uspp_atom_grid_apply_state(ks_env, uspp, fields, hamiltonian, density, &
                                            basis_type, projector_neighbors, force_factors)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(qs_uspp_type), INTENT(INOUT)                  :: uspp
      TYPE(qs_uspp_atom_grid_type), INTENT(IN)           :: fields(:)
      TYPE(dbcsr_p_type), INTENT(INOUT)                  :: hamiltonian(:, :)
      TYPE(dbcsr_p_type), INTENT(IN)                     :: density(:, :)
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL             :: basis_type
      TYPE(neighbor_list_set_p_type), OPTIONAL, POINTER  :: projector_neighbors(:)
      REAL(dp), INTENT(IN), OPTIONAL                     :: force_factors(:)

      INTEGER                                            :: natom, spin
      REAL(dp), ALLOCATABLE                              :: atom_gradient(:, :), cell_gradient(:, :)
      TYPE(sap_int_type), POINTER                        :: sap_int(:)

      CPASSERT(SIZE(fields) == SIZE(hamiltonian, 1))
      NULLIFY (sap_int)
      DO spin = 1, SIZE(fields)
         IF (.NOT. ALLOCATED(fields(spin)%density)) CYCLE
         IF (fields(spin)%calculate_forces) THEN
            IF (.NOT. ASSOCIATED(sap_int)) THEN
               CALL get_ks_env(ks_env, natom=natom)
               ALLOCATE (atom_gradient(3, natom))
               atom_gradient = 0.0_dp
               CALL allocate_cell_gradient(ks_env, cell_gradient)
               CPASSERT(PRESENT(basis_type) .EQV. PRESENT(projector_neighbors))
               IF (PRESENT(basis_type)) THEN
                  CALL build_derivatives_basis(ks_env, sap_int, basis_type, projector_neighbors)
               ELSE
                  CALL build_derivatives(ks_env, sap_int)
               END IF
            END IF
            atom_gradient = 0.0_dp
            IF (ALLOCATED(cell_gradient)) cell_gradient = 0.0_dp
            CALL qs_uspp_atom_grid_apply(uspp, fields(spin), hamiltonian(spin, :), &
                                         density(spin, :), sap_int, atom_gradient, cell_gradient)
            IF (PRESENT(force_factors)) THEN
               CPASSERT(SIZE(force_factors) == SIZE(fields))
               atom_gradient = force_factors(spin)*atom_gradient
               IF (ALLOCATED(cell_gradient)) cell_gradient = force_factors(spin)*cell_gradient
            END IF
            CALL accumulate_gradient(ks_env, atom_gradient, cell_gradient)
         ELSE
            CALL qs_uspp_atom_grid_apply(uspp, fields(spin), hamiltonian(spin, :))
         END IF
      END DO
      IF (ASSOCIATED(sap_int)) THEN
         CALL release_sap_int(sap_int)
      END IF
   END SUBROUTINE qs_uspp_atom_grid_apply_state

! **************************************************************************************************
!> \brief Conservative image bounds for a radial support on an arbitrary triclinic cell.
!> \param cell Auxiliary cell.
!> \param center Original augmentation center.
!> \param points Nonempty Cartesian point batch.
!> \param radius Full supplied radial support; no minimum-image approximation or tail truncation.
!> \param lower First included integer image in each cell direction.
!> \param upper Last included integer image in each cell direction.
! **************************************************************************************************
   SUBROUTINE atom_grid_image_bounds(cell, center, points, radius, lower, upper)
      TYPE(cell_type), INTENT(IN)                        :: cell
      REAL(KIND=dp), INTENT(IN)                          :: center(3), points(:, :), radius
      INTEGER, INTENT(OUT)                               :: lower(3), upper(3)

      INTEGER                                            :: axis
      REAL(KIND=dp)                                      :: padding
      REAL(KIND=dp), ALLOCATABLE                         :: fractional(:, :)

      CPASSERT(SIZE(points, 2) > 0)
      ALLOCATE (fractional(3, SIZE(points, 2)))
      fractional(:, :) = MATMUL(cell%h_inv, points - SPREAD(center, 2, SIZE(points, 2)))
      DO axis = 1, 3
         padding = radius*SQRT(SUM(cell%h_inv(axis, :)**2))
         ! Include numerical boundary points; the radial evaluator enforces the spherical support.
         padding = padding + 32.0_dp*EPSILON(1.0_dp)*MAX(1.0_dp, padding, MAXVAL(ABS(fractional(axis, :))))
         lower(axis) = CEILING(MINVAL(fractional(axis, :)) - padding)
         upper(axis) = FLOOR(MAXVAL(fractional(axis, :)) + padding)
      END DO
   END SUBROUTINE atom_grid_image_bounds

! **************************************************************************************************
!> \brief Prepare USPP centers, unweighted projector overlaps, and the generalized metric.
!>        Called after rebuilding the Gaussian overlap. The NC path allocates nothing.
!> \param ks_env KS environment with current kinds, particles and ABBA projector neighbor list.
!> \param gaussian_overlap Current Gaussian overlap; never modified or aliased by the metric.
! **************************************************************************************************
   SUBROUTINE qs_uspp_prepare(ks_env, gaussian_overlap)
      TYPE(qs_ks_env_type), INTENT(INOUT)                :: ks_env
      TYPE(dbcsr_p_type), INTENT(IN)                     :: gaussian_overlap(:)

      INTEGER                                            :: handle, ikind
      LOGICAL                                            :: do_kpoints
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(neighbor_list_set_p_type), POINTER            :: sab_kp(:), sap_ppnl(:)
      TYPE(particle_type), POINTER                       :: particle_set(:)
      TYPE(qs_kind_type), POINTER                        :: qs_kind_set(:)
      TYPE(qs_uspp_type), POINTER                        :: uspp
      TYPE(sap_int_type), POINTER                        :: sap_int(:)

      CALL get_ks_env(ks_env, uspp=uspp, qs_kind_set=qs_kind_set, particle_set=particle_set)
      CALL qs_uspp_init(uspp, qs_kind_set, particle_set)
      CALL set_ks_env(ks_env, uspp=uspp)
      IF (.NOT. ASSOCIATED(uspp)) RETURN
      CALL timeset('qs_uspp_prepare', handle)
      CALL get_ks_env(ks_env, do_kpoints=do_kpoints, sap_ppnl=sap_ppnl, dft_control=dft_control)
      IF (.NOT. (dft_control%qs_control%gpw .OR. dft_control%qs_control%gapw_xc .OR. &
                 dft_control%qs_control%gapw)) THEN
         CPABORT("USPP/PAW requires GPW, GAPW_XC or GAPW")
      END IF
      IF (dft_control%qs_control%do_ppl_method /= do_ppl_analytic) THEN
         CPABORT("USPP/PAW requires analytic local-potential integrals")
      END IF
      ! Semilocal USPP tau uses an independent transition and an explicitly normalized fixed core.
      IF (dft_control%use_kinetic_energy_density .AND. .NOT. dft_control%use_gauxc) THEN
         ALLOCATE (uspp%kinetic_augmentation(SIZE(qs_kind_set)), uspp%model_core_kinetic(SIZE(qs_kind_set)))
         DO ikind = 1, SIZE(qs_kind_set)
            IF (uspp%augmentation(ikind)%nproj == 0) CYCLE
            ASSOCIATE (potential => qs_kind_set(ikind)%sgp_potential)
               IF (potential%upf_data%is_paw) CYCLE
               BLOCK
                  CHARACTER(LEN=256) :: reason
                  LOGICAL :: available
                  CALL upf_kinetic_reconstruction_init(potential%upf_data, potential%upf_projectors, &
                                                       uspp%kinetic_augmentation(ikind), available, reason)
                  IF (.NOT. available) CALL cp_abort(__LOCATION__, TRIM(reason))
                  CALL upf_core_kinetic_init(potential%upf_data, uspp%model_core_kinetic(ikind))
                  IF (.NOT. uspp%model_core_kinetic(ikind)%available) THEN
                     CALL cp_abort(__LOCATION__, TRIM(uspp%model_core_kinetic(ikind)%failure_reason))
                  END IF
               END BLOCK
            END ASSOCIATE
         END DO
      END IF
      DO ikind = 1, SIZE(qs_kind_set)
         IF (.NOT. uspp%paw(ikind)%active) CYCLE
         ASSOCIATE (potential => qs_kind_set(ikind)%sgp_potential, cache => uspp%paw(ikind))
            CALL paw_onecenter_init(potential%upf_data, potential%upf_projectors, cache%onecenter)
            CALL paw_xc_init(potential%upf_data, qs_kind_set(ikind)%ngrid_ang, cache%xcgrid)
            IF (dft_control%use_kinetic_energy_density .AND. .NOT. dft_control%use_gauxc) THEN
               IF (.NOT. cache%xcgrid%core_kinetic%available .AND. ANY(cache%xcgrid%core /= 0.0_dp)) THEN
                  CALL cp_abort(__LOCATION__, TRIM(cache%xcgrid%core_kinetic%failure_reason))
               END IF
            END IF
            cache%core_energy = potential%upf_data%paw%core_energy
         END ASSOCIATE
      END DO
      CPASSERT(ASSOCIATED(sap_ppnl))
      ALLOCATE (sap_int(SIZE(qs_kind_set)**2))
      CALL build_sap_ints(sap_int, sap_ppnl, qs_kind_set, nder=0)
      IF (do_kpoints) THEN
         CALL get_ks_env(ks_env, kpoints=kpoints, sab_kp=sab_kp)
         CALL qs_uspp_update_images(uspp, gaussian_overlap, sap_int, particle_set, sab_kp, kpoints%cell_to_index)
      ELSE
         CALL qs_uspp_update(uspp, gaussian_overlap(1)%matrix, sap_int, particle_set)
         CALL release_sap_int(sap_int)
      END IF
      CALL timestop(handle)
   END SUBROUTINE qs_uspp_prepare

! **************************************************************************************************
!> \brief Initialize a matrix on the union of its ordinary and USPP projector-mediated blocks.
!>        Called after ordinary neighbor-list allocation, before assigning matrix values.
!> \param ks_env KS environment with current USPP geometry data.
!> \param matrix Newly allocated AO matrix, zeroed after expanding its block pattern.
! **************************************************************************************************
   SUBROUTINE qs_uspp_init_ao_matrix(ks_env, matrix)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix

      TYPE(qs_uspp_type), POINTER                        :: uspp

      CALL get_ks_env(ks_env, uspp=uspp)
      IF (.NOT. ASSOCIATED(uspp)) RETURN
      CALL qs_uspp_init_matrix_state(uspp, matrix)
   END SUBROUTINE qs_uspp_init_ao_matrix

! **************************************************************************************************
!> \brief Include projector connectivity for an explicit ORB or AUX state.
!> \param uspp ...
!> \param matrix ...
! **************************************************************************************************
   SUBROUTINE qs_uspp_init_matrix_state(uspp, matrix)
      TYPE(qs_uspp_type), INTENT(IN)                     :: uspp
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix

      CPASSERT(ASSOCIATED(uspp%ao_pattern))
      CALL dbcsr_add(matrix, uspp%ao_pattern, 1.0_dp, 1.0_dp)
      CALL dbcsr_set(matrix, 0.0_dp)
   END SUBROUTINE qs_uspp_init_matrix_state

! **************************************************************************************************
!> \brief Add USPP augmentation (or its spatial derivative) to a Quickstep density grid.
!>        The reciprocal and real representations are both updated. Ordinary kinds do nothing.
!> \param ks_env KS environment with the current projector overlaps.
!> \param density AO density for one spin, including occupations.
!> \param rho_g Existing reciprocal density, normalized as integral/cell volume.
!> \param rho_r Corresponding real density, overwritten from the updated reciprocal grid.
!> \param total_rho Optional signed electronic charge, decremented by the augmentation charge.
!> \param derivative Optional Cartesian derivative orders; total_rho must then be absent.
!> \param kinetic Select independent kinetic transitions instead of charge augmentation.
! **************************************************************************************************
   SUBROUTINE qs_uspp_collocate_images(ks_env, density, rho_g, rho_r, total_rho, derivative, kinetic)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(dbcsr_p_type), INTENT(IN)                     :: density(:)
      TYPE(pw_c1d_gs_type), INTENT(INOUT)                :: rho_g
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: rho_r
      REAL(dp), INTENT(INOUT), OPTIONAL                  :: total_rho
      INTEGER, INTENT(IN), OPTIONAL                      :: derivative(3)
      LOGICAL, INTENT(IN), OPTIONAL                      :: kinetic

      TYPE(qs_uspp_type), POINTER                        :: uspp

      CALL get_ks_env(ks_env, uspp=uspp)
      IF (.NOT. ASSOCIATED(uspp)) RETURN
      CALL qs_uspp_collocate_state(uspp, density, rho_g, rho_r, total_rho, derivative, kinetic)
   END SUBROUTINE qs_uspp_collocate_images

! **************************************************************************************************
!> \brief Collocate an explicitly selected ORB or AUX projector state.
!> \param uspp ...
!> \param density ...
!> \param rho_g ...
!> \param rho_r ...
!> \param total_rho ...
!> \param derivative ...
!> \param kinetic ...
! **************************************************************************************************
   SUBROUTINE qs_uspp_collocate_state(uspp, density, rho_g, rho_r, total_rho, derivative, kinetic)
      TYPE(qs_uspp_type), INTENT(INOUT)                  :: uspp
      TYPE(dbcsr_p_type), INTENT(IN)                     :: density(:)
      TYPE(pw_c1d_gs_type), INTENT(INOUT)                :: rho_g
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: rho_r
      REAL(KIND=dp), INTENT(INOUT), OPTIONAL             :: total_rho
      INTEGER, INTENT(IN), OPTIONAL                      :: derivative(3)
      LOGICAL, INTENT(IN), OPTIONAL                      :: kinetic

      INTEGER                                            :: handle
      LOGICAL                                            :: kinetic_field
      REAL(KIND=dp)                                      :: charge
      TYPE(pw_c1d_gs_type)                               :: augmentation_g
      TYPE(pw_r3d_rs_type)                               :: augmentation_r

      kinetic_field = .FALSE.
      IF (PRESENT(kinetic)) kinetic_field = kinetic
      IF (kinetic_field) THEN
         IF (.NOT. ALLOCATED(uspp%kinetic_augmentation)) RETURN
      END IF
      CALL timeset('qs_uspp_collocate', handle)
      CPASSERT(ASSOCIATED(rho_g%pw_grid, rho_r%pw_grid))
      CALL check_grid(uspp, rho_g%pw_grid)
      CALL augmentation_g%create(rho_g%pw_grid)
      augmentation_g%array = CMPLX(0.0_dp, 0.0_dp, dp)
      IF (kinetic_field) THEN
         CPASSERT(.NOT. PRESENT(total_rho))
         CALL collocate_density(uspp, density, rho_g%pw_grid, augmentation_g%array, charge, uspp%kinetic_augmentation)
      ELSE
         CALL collocate_density(uspp, density, rho_g%pw_grid, augmentation_g%array, charge)
      END IF
      ! On even FULLSPACE grids, project the Nyquist coefficients onto a real grid function.
      ! This same real projection is implicit in the adjoint's FFT of a real potential.
      IF (rho_g%pw_grid%grid_span == FULLSPACE .AND. ANY(MODULO(rho_g%pw_grid%npts, 2) == 0)) THEN
         CALL augmentation_r%create(rho_g%pw_grid)
         CALL pw_transfer(augmentation_g, augmentation_r)
         CALL pw_transfer(augmentation_r, augmentation_g)
         CALL augmentation_r%release()
      END IF
      IF (PRESENT(derivative)) THEN
         CPASSERT(.NOT. PRESENT(total_rho))
         CALL pw_derive(augmentation_g, derivative)
      ELSE IF (PRESENT(total_rho)) THEN
         total_rho = total_rho - charge
      END IF
      CALL pw_axpy(augmentation_g, rho_g)
      CALL pw_transfer(rho_g, rho_r)
      CALL augmentation_g%release()
      CALL timestop(handle)
   END SUBROUTINE qs_uspp_collocate_state

! **************************************************************************************************
!> \brief Integrate a Quickstep real local potential against USPP augmentation.
!>        This is the discrete adjoint of qs_uspp_collocate, including FFT normalization.
!> \param ks_env KS environment with the current projector overlaps.
!> \param potential Real potential; values include the volume element when weighted is true.
!> \param hamiltonian AO Hamiltonian for one spin, incremented in place.
!> \param weighted Whether potential already includes the real grid's volume element.
!> \param density Optional force-path density; adds the augmentation energy gradients to the KS forces.
!> \param kinetic Select independent kinetic transitions instead of charge augmentation.
! **************************************************************************************************
   SUBROUTINE qs_uspp_integrate_images(ks_env, potential, hamiltonian, weighted, density, kinetic)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: potential
      TYPE(dbcsr_p_type), INTENT(INOUT)                  :: hamiltonian(:)
      LOGICAL, INTENT(IN)                                :: weighted
      TYPE(dbcsr_p_type), INTENT(IN), OPTIONAL           :: density(:)
      LOGICAL, INTENT(IN), OPTIONAL                      :: kinetic

      TYPE(qs_uspp_type), POINTER                        :: uspp

      CALL get_ks_env(ks_env, uspp=uspp)
      IF (.NOT. ASSOCIATED(uspp)) RETURN
      CALL qs_uspp_integrate_state(ks_env, uspp, potential, hamiltonian, weighted, density, kinetic)
   END SUBROUTINE qs_uspp_integrate_images

! **************************************************************************************************
!> \brief Grid adjoint for an explicit projector state, including its own basis derivatives.
!> \param ks_env ...
!> \param uspp ...
!> \param potential ...
!> \param hamiltonian ...
!> \param weighted ...
!> \param density ...
!> \param kinetic ...
!> \param basis_type ...
!> \param projector_neighbors ...
!> \param force_factor Optional ADMM scaling of geometry adjoints; the Hamiltonian is unscaled.
! **************************************************************************************************
   SUBROUTINE qs_uspp_integrate_state(ks_env, uspp, potential, hamiltonian, weighted, density, kinetic, &
                                      basis_type, projector_neighbors, force_factor)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(qs_uspp_type), INTENT(INOUT)                  :: uspp
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: potential
      TYPE(dbcsr_p_type), INTENT(INOUT)                  :: hamiltonian(:)
      LOGICAL, INTENT(IN)                                :: weighted
      TYPE(dbcsr_p_type), INTENT(IN), OPTIONAL           :: density(:)
      LOGICAL, INTENT(IN), OPTIONAL                      :: kinetic
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL             :: basis_type
      TYPE(neighbor_list_set_p_type), OPTIONAL, POINTER  :: projector_neighbors(:)
      REAL(dp), INTENT(IN), OPTIONAL                     :: force_factor

      INTEGER                                            :: natom
      LOGICAL                                            :: kinetic_field
      REAL(KIND=dp), ALLOCATABLE                         :: atom_gradient(:, :), cell_gradient(:, :)
      TYPE(sap_int_type), POINTER                        :: sap_int(:)

      kinetic_field = .FALSE.
      IF (PRESENT(kinetic)) kinetic_field = kinetic
      IF (kinetic_field) THEN
         IF (.NOT. ALLOCATED(uspp%kinetic_augmentation)) RETURN
      END IF
      IF (PRESENT(density)) THEN
         CALL get_ks_env(ks_env, natom=natom)
         ALLOCATE (atom_gradient(3, natom))
         atom_gradient = 0.0_dp
         CALL allocate_cell_gradient(ks_env, cell_gradient)
         CPASSERT(PRESENT(basis_type) .EQV. PRESENT(projector_neighbors))
         IF (PRESENT(basis_type)) THEN
            CALL build_derivatives_basis(ks_env, sap_int, basis_type, projector_neighbors)
         ELSE
            CALL build_derivatives(ks_env, sap_int)
         END IF
         IF (kinetic_field) THEN
            CALL integrate_grid(uspp, potential, hamiltonian, weighted, density, sap_int, atom_gradient, cell_gradient, &
                                radial_fields=uspp%kinetic_augmentation)
         ELSE
            CALL integrate_grid(uspp, potential, hamiltonian, weighted, density, sap_int, atom_gradient, cell_gradient)
         END IF
         IF (PRESENT(force_factor)) THEN
            atom_gradient = force_factor*atom_gradient
            IF (ALLOCATED(cell_gradient)) cell_gradient = force_factor*cell_gradient
         END IF
         CALL accumulate_gradient(ks_env, atom_gradient, cell_gradient)
         CALL release_sap_int(sap_int)
      ELSE
         IF (kinetic_field) THEN
            CALL integrate_grid(uspp, potential, hamiltonian, weighted, radial_fields=uspp%kinetic_augmentation)
         ELSE
            CALL integrate_grid(uspp, potential, hamiltonian, weighted)
         END IF
      END IF
   END SUBROUTINE qs_uspp_integrate_state

! **************************************************************************************************
!> \brief Real-grid adjoint shared by the ionic and density-dependent local potentials.
!> \param uspp Current projector state.
!> \param potential Real local potential.
!> \param hamiltonian AO Hamiltonian, incremented in place.
!> \param weighted Whether the real potential includes dV.
!> \param density Optional force-path AO density.
!> \param sap_int Optional image-resolved projector derivatives.
!> \param atom_gradient Optional local atomic energy gradients, incremented without MPI reduction.
!> \param cell_gradient Optional strain energy gradient, incremented without MPI reduction.
!> \param radial_fields Optional independent transitions in the same projector convention.
! **************************************************************************************************
   SUBROUTINE integrate_grid(uspp, potential, hamiltonian, weighted, density, sap_int, atom_gradient, cell_gradient, radial_fields)
      TYPE(qs_uspp_type), INTENT(INOUT), TARGET          :: uspp
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: potential
      TYPE(dbcsr_p_type), INTENT(INOUT)                  :: hamiltonian(:)
      LOGICAL, INTENT(IN)                                :: weighted
      TYPE(dbcsr_p_type), INTENT(IN), OPTIONAL           :: density(:)
      TYPE(sap_int_type), INTENT(IN), OPTIONAL           :: sap_int(:)
      REAL(KIND=dp), INTENT(INOUT), OPTIONAL             :: atom_gradient(:, :), cell_gradient(3, 3)
      TYPE(uspp_augmentation_type), INTENT(INOUT), &
         OPTIONAL, TARGET                                :: radial_fields(:)

      INTEGER                                            :: center, handle, idir, ikind, jdir, n
      LOGICAL                                            :: found
      REAL(KIND=dp), ALLOCATABLE                         :: derivative(:, :, :), dm(:, :), &
                                                            matrix(:, :), multiplicity(:), &
                                                            strain(:, :, :, :)
      REAL(KIND=dp), POINTER                             :: block(:, :)
      TYPE(dbcsr_type)                                   :: onsite, projected
      TYPE(mp_comm_type)                                 :: group
      TYPE(pw_c1d_gs_type)                               :: potential_g
      TYPE(uspp_augmentation_type), POINTER              :: cache(:)

      cache => uspp%augmentation
      IF (PRESENT(radial_fields)) THEN
         CPASSERT(SIZE(radial_fields) == SIZE(uspp%augmentation))
         cache => radial_fields
      END IF
      CALL timeset('qs_uspp_integrate', handle)
      CALL check_grid(uspp, potential%pw_grid)
      DO ikind = 1, SIZE(cache)
         IF (cache(ikind)%nproj == 0) CYCLE
         CALL uspp_augmentation_prepare_fourier(cache(ikind), potential%pw_grid%g)
      END DO
      CALL potential_g%create(potential%pw_grid)
      CALL pw_transfer(potential, potential_g)
      IF (weighted) CALL pw_scale(potential_g, 1.0_dp/potential%pw_grid%dvol)
      ALLOCATE (multiplicity(SIZE(potential_g%array)))
      multiplicity = 1.0_dp
      IF (potential%pw_grid%grid_span == HALFSPACE) THEN
         multiplicity = 2.0_dp
         IF (potential%pw_grid%have_g0) multiplicity(1) = 1.0_dp
      END IF
      CPASSERT(PRESENT(density) .EQV. PRESENT(sap_int))
      CPASSERT(PRESENT(density) .EQV. PRESENT(atom_gradient))
      CALL dbcsr_get_info(uspp%overlaps, group=group)
      CALL uspp_projector_create_onsite(uspp%overlaps, onsite)
      IF (PRESENT(density)) CALL project_density(uspp, density, projected)
      DO center = 1, SIZE(uspp%center_atoms)
         IF (cache(uspp%center_kinds(center))%nproj == 0) CYCLE
         n = uspp%projector_sizes(center)
         ALLOCATE (matrix(n, n))
         IF (PRESENT(density)) THEN
            ALLOCATE (dm(n, n), derivative(n, n, 3))
            IF (PRESENT(cell_gradient)) ALLOCATE (strain(n, n, 3, 3))
            dm = 0.0_dp
            CALL dbcsr_get_block_p(projected, center, center, block, found)
            IF (found) dm(:, :) = block
            CALL group%sum(dm)
            CALL uspp_augmentation_integrate(cache(uspp%center_kinds(center)), &
                                             potential%pw_grid%g, uspp%position(:, center), potential%pw_grid%vol, &
                                             potential_g%array, multiplicity, matrix, derivative, strain)
            DO idir = 1, 3
               atom_gradient(idir, uspp%center_atoms(center)) = atom_gradient(idir, uspp%center_atoms(center)) &
                                                                + SUM(dm*derivative(:, :, idir))
            END DO
            IF (PRESENT(cell_gradient)) THEN
               DO jdir = 1, 3
                  DO idir = 1, 3
                     cell_gradient(idir, jdir) = cell_gradient(idir, jdir) + SUM(dm*strain(:, :, idir, jdir))
                  END DO
               END DO
               DEALLOCATE (strain)
            END IF
            DEALLOCATE (dm, derivative)
         ELSE
            CALL uspp_augmentation_integrate(cache(uspp%center_kinds(center)), &
                                             potential%pw_grid%g, uspp%position(:, center), potential%pw_grid%vol, &
                                             potential_g%array, multiplicity, matrix)
         END IF
         CALL group%sum(matrix)
         CALL dbcsr_get_block_p(onsite, center, center, block, found)
         IF (found) block(:, :) = matrix
         DEALLOCATE (matrix)
      END DO
      CALL apply_operator(uspp, onsite, hamiltonian, density, sap_int, atom_gradient, cell_gradient)
      CALL dbcsr_release(onsite)
      IF (PRESENT(density)) CALL dbcsr_release(projected)
      CALL potential_g%release()
      CALL timestop(handle)
   END SUBROUTINE integrate_grid

! **************************************************************************************************
!> \brief Add B (D0 + integral Q V_sr) B^T to the core Hamiltonian.
!>        V_sr includes every non-ghost UPF ion, also NC ions without augmentation.
!>        The Gaussian-core Coulomb potential enters later through the Hartree grid.
!>        Keeping this density-independent term in H_core supplies its energy via Tr(P H_core).
!> \param uspp Current geometry-dependent USPP state.
!> \param qs_kind_set Kinds with original UPF local potentials and native projector ordering.
!> \param particle_set All ions, not only augmentation centers.
!> \param grid Same FFT grid as the physical electron density.
!> \param hamiltonian Core Hamiltonian, incremented in place.
!> \param density Optional total AO density for the force path.
!> \param energy_weighted Optional total energy-weighted AO density.
!> \param sap_int Optional uncompressed projector derivatives.
!> \param atom_gradient Optional local atomic energy gradients, incremented without MPI reduction.
!> \param cell_gradient Optional strain energy gradient, incremented without MPI reduction.
! **************************************************************************************************
   SUBROUTINE qs_uspp_core_matrix_images(uspp, qs_kind_set, particle_set, grid, hamiltonian, &
                                         density, energy_weighted, sap_int, atom_gradient, cell_gradient)
      TYPE(qs_uspp_type), INTENT(INOUT)                  :: uspp
      TYPE(qs_kind_type), INTENT(INOUT)                  :: qs_kind_set(:)
      TYPE(particle_type), INTENT(IN)                    :: particle_set(:)
      TYPE(pw_grid_type), POINTER                        :: grid
      TYPE(dbcsr_p_type), INTENT(INOUT)                  :: hamiltonian(:)
      TYPE(dbcsr_p_type), INTENT(IN), OPTIONAL           :: density(:), energy_weighted(:)
      TYPE(sap_int_type), INTENT(IN), OPTIONAL           :: sap_int(:)
      REAL(KIND=dp), INTENT(INOUT), OPTIONAL             :: atom_gradient(:, :), cell_gradient(3, 3)

      COMPLEX(KIND=dp), ALLOCATABLE                      :: ionic_derivative(:, :), ionic_dg(:, :), &
                                                            ionic_g(:)
      INTEGER                                            :: atom, center, handle, idir, ikind, jdir
      LOGICAL                                            :: found
      REAL(KIND=dp)                                      :: charge
      REAL(KIND=dp), ALLOCATABLE                         :: multiplicity(:)
      REAL(KIND=dp), POINTER                             :: block(:, :)
      TYPE(dbcsr_type)                                   :: bare
      TYPE(pw_c1d_gs_type)                               :: potential_g, rho_g
      TYPE(pw_r3d_rs_type)                               :: potential_r, rho_r

      CALL timeset('qs_uspp_core_matrix', handle)
      CALL check_grid(uspp, grid)
      CPASSERT(PRESENT(density) .EQV. PRESENT(energy_weighted))
      CPASSERT(PRESENT(density) .EQV. PRESENT(sap_int))
      CPASSERT(PRESENT(density) .EQV. PRESENT(atom_gradient))
      IF (PRESENT(cell_gradient)) THEN
         CPASSERT(PRESENT(density))
      END IF
      CALL uspp_projector_create_onsite(uspp%overlaps, bare)
      DO center = 1, SIZE(uspp%center_kinds)
         ikind = uspp%center_kinds(center)
         CALL dbcsr_get_block_p(bare, center, center, block, found)
         IF (found) block(:, :) = qs_kind_set(ikind)%sgp_potential%upf_projectors%d
      END DO
      CALL apply_operator(uspp, bare, hamiltonian, density, sap_int, atom_gradient, cell_gradient)
      IF (PRESENT(density)) THEN
         DO center = 1, SIZE(uspp%center_kinds)
            ikind = uspp%center_kinds(center)
            CALL dbcsr_get_block_p(bare, center, center, block, found)
            IF (found) block(:, :) = uspp%augmentation(ikind)%q
         END DO
         CALL apply_operator(uspp, bare, density=energy_weighted, sap_int=sap_int, &
                             atom_gradient=atom_gradient, cell_gradient=cell_gradient, alpha=-1.0_dp)
      END IF
      CALL dbcsr_release(bare)

      CALL potential_g%create(grid)
      CALL potential_r%create(grid)
      ALLOCATE (ionic_g(SIZE(potential_g%array)))
      potential_g%array = CMPLX(0.0_dp, 0.0_dp, dp)
      IF (PRESENT(density)) THEN
         CALL rho_g%create(grid)
         CALL rho_r%create(grid)
         rho_g%array = CMPLX(0.0_dp, 0.0_dp, dp)
         CALL collocate_density(uspp, density, grid, rho_g%array, charge)
         ! Self-adjoint real projection is needed before contracting a raw ionic derivative.
         CALL pw_transfer(rho_g, rho_r)
         CALL pw_transfer(rho_r, rho_g)
         ALLOCATE (ionic_derivative(3, SIZE(ionic_g)), multiplicity(SIZE(ionic_g)))
         IF (PRESENT(cell_gradient)) ALLOCATE (ionic_dg(3, SIZE(ionic_g)))
         multiplicity = 1.0_dp
         IF (grid%grid_span == HALFSPACE) THEN
            multiplicity = 2.0_dp
            IF (grid%have_g0) multiplicity(1) = 1.0_dp
         END IF
      END IF
      DO atom = 1, SIZE(particle_set)
         ikind = particle_set(atom)%atomic_kind%kind_number
         IF (qs_kind_set(ikind)%ghost) CYCLE
         CPASSERT(ASSOCIATED(qs_kind_set(ikind)%sgp_potential))
         ASSOCIATE (potential => qs_kind_set(ikind)%sgp_potential)
            CPASSERT(ALLOCATED(potential%upf_data))
            IF (potential%upf_local%radial%degree < 0) THEN
               CALL upf_local_init(potential%upf_data, potential%alpha_core_charge, 0, potential%upf_local)
            END IF
            IF (PRESENT(density)) THEN
               CALL upf_local_fourier(potential%upf_local, grid%g, particle_set(atom)%r, grid%vol, &
                                      ionic_g, ionic_derivative, ionic_dg)
               DO idir = 1, 3
                  atom_gradient(idir, atom) = atom_gradient(idir, atom) + grid%vol* &
                                              SUM(multiplicity*REAL(CONJG(rho_g%array)*ionic_derivative(idir, :), KIND=dp))
               END DO
               IF (PRESENT(cell_gradient)) THEN
                  ! The ionic phase is invariant at fixed fractional coordinates.
                  ! The -V/volume term cancels the explicit volume weight in integral rho_aug V_sr.
                  DO jdir = 1, 3
                     ionic_dg(jdir, :) = ionic_dg(jdir, :) + CMPLX(0.0_dp, particle_set(atom)%r(jdir), dp)*ionic_g
                     DO idir = 1, 3
                        cell_gradient(idir, jdir) = cell_gradient(idir, jdir) - grid%vol* &
                                               SUM(multiplicity*REAL(CONJG(rho_g%array)*grid%g(idir, :)*ionic_dg(jdir, :), KIND=dp))
                     END DO
                  END DO
               END IF
            ELSE
               CALL upf_local_fourier(potential%upf_local, grid%g, particle_set(atom)%r, grid%vol, ionic_g)
            END IF
            potential_g%array = potential_g%array + ionic_g
         END ASSOCIATE
      END DO
      ! Explicit real projection also treats aliased Nyquist pairs consistently on even grids.
      CALL pw_transfer(potential_g, potential_r)
      CALL integrate_grid(uspp, potential_r, hamiltonian, weighted=.FALSE., &
                          density=density, sap_int=sap_int, atom_gradient=atom_gradient, cell_gradient=cell_gradient)
      IF (PRESENT(density)) THEN
         CALL rho_g%release()
         CALL rho_r%release()
      END IF
      CALL potential_g%release()
      CALL potential_r%release()
      CALL timestop(handle)
   END SUBROUTINE qs_uspp_core_matrix_images

! **************************************************************************************************
!> \brief Core USPP operator and its stationary energy gradients, including the metric Pulay term.
!> \param ks_env Current geometry and force accumulators.
!> \param grid Physical density FFT grid.
!> \param hamiltonian Core Hamiltonian to increment.
!> \param density Spin densities; summed without modifying the originals.
!> \param weighted_density Total energy-weighted density from both spins.
! **************************************************************************************************
   SUBROUTINE qs_uspp_core_forces(ks_env, grid, hamiltonian, density, weighted_density)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(pw_grid_type), POINTER                        :: grid
      TYPE(dbcsr_p_type), INTENT(INOUT)                  :: hamiltonian(:)
      TYPE(dbcsr_p_type), INTENT(IN)                     :: density(:, :), weighted_density(:)

      INTEGER                                            :: image, ispin
      REAL(KIND=dp), ALLOCATABLE                         :: atom_gradient(:, :), cell_gradient(:, :)
      TYPE(dbcsr_p_type), ALLOCATABLE                    :: total_density(:)
      TYPE(particle_type), POINTER                       :: particle_set(:)
      TYPE(qs_kind_type), POINTER                        :: qs_kind_set(:)
      TYPE(qs_uspp_type), POINTER                        :: uspp
      TYPE(sap_int_type), POINTER                        :: sap_int(:)

      CALL get_ks_env(ks_env, uspp=uspp, qs_kind_set=qs_kind_set, particle_set=particle_set)
      CPASSERT(ASSOCIATED(uspp))
      ALLOCATE (total_density(SIZE(density, 2)))
      DO image = 1, SIZE(density, 2)
         ALLOCATE (total_density(image)%matrix)
         CALL dbcsr_copy(total_density(image)%matrix, density(1, image)%matrix)
         DO ispin = 2, SIZE(density, 1)
            CALL dbcsr_add(total_density(image)%matrix, density(ispin, image)%matrix, 1.0_dp, 1.0_dp)
         END DO
      END DO
      ALLOCATE (atom_gradient(3, SIZE(particle_set)))
      atom_gradient = 0.0_dp
      CALL allocate_cell_gradient(ks_env, cell_gradient)
      CALL build_derivatives(ks_env, sap_int)
      CALL qs_uspp_core_matrix(uspp, qs_kind_set, particle_set, grid, hamiltonian, &
                               total_density, weighted_density, sap_int, atom_gradient, cell_gradient)
      CALL accumulate_gradient(ks_env, atom_gradient, cell_gradient)
      CALL release_sap_int(sap_int)
      DO image = 1, SIZE(total_density)
         CALL dbcsr_release(total_density(image)%matrix)
         DEALLOCATE (total_density(image)%matrix)
      END DO
   END SUBROUTINE qs_uspp_core_forces

! **************************************************************************************************
!> \brief PAW one-center energies, Hamiltonian and projector force/strain derivatives in every METHOD.
!>        Radial data are fixed in Cartesian atom-centered coordinates. Geometry enters only through
!>        the projected density, so its transpose supplies the complete one-center geometry response.
!> \param ks_env Prepared generalized overlap, original PAW caches and force accumulators.
!> \param xc_section Same XC section as for the smooth density.
!> \param density Spin and image-resolved AO densities, including occupations and k-point weights.
!> \param hartree PAW Hartree correction, reset and returned replicated.
!> \param exc PAW XC correction, reset and returned replicated.
!> \param core Frozen-core energy constant, reset and returned replicated.
!> \param calculate_forces Whether to differentiate projectors at fixed AO density.
!> \param hamiltonian Optional spin and image AO matrices; absent for energy-only evaluation.
!> \param semilocal_xc ...
! **************************************************************************************************
   SUBROUTINE qs_paw_correction(ks_env, xc_section, density, hartree, exc, core, calculate_forces, hamiltonian, semilocal_xc)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(section_vals_type), POINTER                   :: xc_section
      TYPE(dbcsr_p_type), INTENT(IN)                     :: density(:, :)
      REAL(dp), INTENT(OUT)                              :: hartree, exc, core
      LOGICAL, INTENT(IN)                                :: calculate_forces
      TYPE(dbcsr_p_type), INTENT(INOUT), OPTIONAL        :: hamiltonian(:, :)
      LOGICAL, INTENT(IN), OPTIONAL                      :: semilocal_xc

      INTEGER                                            :: center, handle, ikind, n, natom, ns, spin
      LOGICAL                                            :: found
      REAL(dp)                                           :: density_cutoff, eh, energies(3), ex, &
                                                            gradient_cutoff
      REAL(dp), ALLOCATABLE                              :: atom_gradient(:, :), &
                                                            cell_gradient(:, :), dm(:, :, :), &
                                                            vh(:, :), vxc(:, :, :)
      REAL(dp), POINTER                                  :: block(:, :)
      TYPE(dbcsr_type), ALLOCATABLE                      :: onsite(:), projected(:)
      TYPE(mp_comm_type)                                 :: group
      TYPE(qs_uspp_type), POINTER                        :: uspp
      TYPE(sap_int_type), POINTER                        :: sap_int(:)
      TYPE(section_vals_type), POINTER                   :: functional

      hartree = 0.0_dp
      exc = 0.0_dp
      core = 0.0_dp
      CALL get_ks_env(ks_env, uspp=uspp)
      IF (.NOT. ASSOCIATED(uspp)) RETURN
      IF (.NOT. ALLOCATED(uspp%paw)) RETURN
      IF (.NOT. ANY(uspp%paw%active)) RETURN
      CALL timeset('qs_paw_correction', handle)
      ns = SIZE(density, 1)
      CPASSERT(ns == 1 .OR. ns == 2)
      IF (PRESENT(hamiltonian)) THEN
         CPASSERT(ALL(SHAPE(hamiltonian) == SHAPE(density)))
      END IF
      functional => section_vals_get_subs_vals(xc_section, 'XC_FUNCTIONAL')
      CALL section_vals_val_get(xc_section, 'DENSITY_CUTOFF', r_val=density_cutoff)
      CALL section_vals_val_get(xc_section, 'GRADIENT_CUTOFF', r_val=gradient_cutoff)
      CALL dbcsr_get_info(uspp%overlaps, group=group)
      ALLOCATE (projected(ns), onsite(ns))
      DO spin = 1, ns
         CALL project_density(uspp, density(spin, :), projected(spin))
         CALL uspp_projector_create_onsite(uspp%overlaps, onsite(spin))
      END DO
      energies = 0.0_dp
      DO center = 1, SIZE(uspp%center_atoms)
         ikind = uspp%center_kinds(center)
         IF (.NOT. uspp%paw(ikind)%active) CYCLE
         ! All spin blocks have the same owner. Evaluate and count each center once.
         CALL dbcsr_get_block_p(projected(1), center, center, block, found)
         IF (.NOT. found) CYCLE
         n = uspp%projector_sizes(center)
         ALLOCATE (dm(n, n, ns), vh(n, n), vxc(n, n, ns))
         DO spin = 1, ns
            CALL dbcsr_get_block_p(projected(spin), center, center, block, found)
            CPASSERT(found)
            dm(:, :, spin) = block
         END DO
         ASSOCIATE (cache => uspp%paw(ikind))
            CALL paw_onecenter_hartree(cache%onecenter, SUM(dm, DIM=3), eh, vh)
            ex = 0.0_dp
            vxc = 0.0_dp
            IF (.NOT. PRESENT(semilocal_xc)) THEN
               CALL paw_xc_evaluate(cache%xcgrid, cache%onecenter, functional, dm, density_cutoff, gradient_cutoff, ex, vxc)
            ELSE IF (semilocal_xc) THEN
               CALL paw_xc_evaluate(cache%xcgrid, cache%onecenter, functional, dm, density_cutoff, gradient_cutoff, ex, vxc)
            END IF
            energies = energies + [eh, ex, cache%core_energy]
         END ASSOCIATE
         DO spin = 1, ns
            CALL dbcsr_get_block_p(onsite(spin), center, center, block, found)
            CPASSERT(found)
            block = vh + vxc(:, :, spin)
         END DO
         DEALLOCATE (dm, vh, vxc)
      END DO
      CALL group%sum(energies)
      hartree = energies(1)
      exc = energies(2)
      core = energies(3)
      NULLIFY (sap_int)
      IF (calculate_forces) THEN
         CALL get_ks_env(ks_env, natom=natom)
         ALLOCATE (atom_gradient(3, natom))
         atom_gradient = 0.0_dp
         CALL allocate_cell_gradient(ks_env, cell_gradient)
         CALL build_derivatives(ks_env, sap_int)
      END IF
      DO spin = 1, ns
         IF (calculate_forces) THEN
            IF (PRESENT(hamiltonian)) THEN
               CALL apply_operator(uspp, onsite(spin), hamiltonian(spin, :), density(spin, :), &
                                   sap_int, atom_gradient, cell_gradient)
            ELSE
               CALL apply_operator(uspp, onsite(spin), density=density(spin, :), &
                                   sap_int=sap_int, atom_gradient=atom_gradient, cell_gradient=cell_gradient)
            END IF
         ELSE IF (PRESENT(hamiltonian)) THEN
            CALL apply_operator(uspp, onsite(spin), hamiltonian(spin, :))
         END IF
         CALL dbcsr_release(projected(spin))
         CALL dbcsr_release(onsite(spin))
      END DO
      IF (calculate_forces) THEN
         CALL accumulate_gradient(ks_env, atom_gradient, cell_gradient)
         CALL release_sap_int(sap_int)
      END IF
      CALL timestop(handle)
   END SUBROUTINE qs_paw_correction

! **************************************************************************************************
!> \brief Build image-resolved derivatives only on the force path.
!> \param ks_env Current ABBA neighbor list and native radial projector data.
!> \param sap_int Newly allocated uncompressed projector overlaps and derivatives.
! **************************************************************************************************
   SUBROUTINE build_derivatives(ks_env, sap_int)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(sap_int_type), POINTER                        :: sap_int(:)

      TYPE(neighbor_list_set_p_type), POINTER            :: sap_ppnl(:)
      TYPE(qs_kind_type), POINTER                        :: qs_kind_set(:)

      CALL get_ks_env(ks_env, qs_kind_set=qs_kind_set, sap_ppnl=sap_ppnl)
      ALLOCATE (sap_int(SIZE(qs_kind_set)**2))
      CALL build_sap_ints(sap_int, sap_ppnl, qs_kind_set, nder=1)
      CALL sap_sort(sap_int)
   END SUBROUTINE build_derivatives

! **************************************************************************************************
!> \brief Original projector derivatives for an explicitly selected independent Gaussian basis.
!> \param ks_env ...
!> \param sap_int ...
!> \param basis_type ...
!> \param projector_neighbors ...
! **************************************************************************************************
   SUBROUTINE build_derivatives_basis(ks_env, sap_int, basis_type, projector_neighbors)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(sap_int_type), POINTER                        :: sap_int(:)
      CHARACTER(LEN=*), INTENT(IN)                       :: basis_type
      TYPE(neighbor_list_set_p_type), POINTER            :: projector_neighbors(:)

      TYPE(qs_kind_type), POINTER                        :: kinds(:)

      CPASSERT(ASSOCIATED(projector_neighbors))
      CALL get_ks_env(ks_env, qs_kind_set=kinds)
      ALLOCATE (sap_int(SIZE(kinds)**2))
      CALL build_sap_ints_basis(sap_int, projector_neighbors, kinds, basis_type, nder=1)
      CALL sap_sort(sap_int)
   END SUBROUTINE build_derivatives_basis

! **************************************************************************************************
!> \brief Native geometry derivative of an independent-basis overlap energy adjoint.
!>        A general rectangular weight has no implicit AO symmetry factor. For ADMM's
!>        conventional overlap-force matrices, weight_factor=-1 converts to dE/dS.
!> \param ks_env ...
!> \param left ...
!> \param right ...
!> \param left_basis ...
!> \param right_basis ...
!> \param left_neighbors ...
!> \param right_neighbors ...
!> \param weight ...
!> \param weight_factor ...
! **************************************************************************************************
   SUBROUTINE qs_uspp_overlap_images(ks_env, left, right, left_basis, right_basis, &
                                     left_neighbors, right_neighbors, weight, weight_factor)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(qs_uspp_type), INTENT(IN)                     :: left, right
      CHARACTER(LEN=*), INTENT(IN)                       :: left_basis, right_basis
      TYPE(neighbor_list_set_p_type), POINTER            :: left_neighbors(:), right_neighbors(:)
      TYPE(dbcsr_p_type), INTENT(IN)                     :: weight(:)
      REAL(dp), INTENT(IN)                               :: weight_factor

      CHARACTER                                          :: matrix_type
      INTEGER                                            :: atom, center, image, natom
      INTEGER, ALLOCATABLE                               :: atom_kinds(:)
      INTEGER, POINTER                                   :: cells(:, :)
      LOGICAL                                            :: do_kpoints, found
      REAL(dp), ALLOCATABLE                              :: atom_gradient(:, :), cell_gradient(:, :)
      REAL(dp), POINTER                                  :: block(:, :)
      TYPE(dbcsr_p_type)                                 :: energy_weight(SIZE(weight))
      TYPE(dbcsr_type)                                   :: left_gradient, onsite, right_gradient
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(particle_type), POINTER                       :: particles(:)
      TYPE(sap_int_type), POINTER                        :: left_sap(:), right_sap(:)

      CPASSERT(ALL(left%center_atoms == right%center_atoms))
      CPASSERT(ALL(left%center_kinds == right%center_kinds))
      CALL get_ks_env(ks_env, natom=natom, do_kpoints=do_kpoints, particle_set=particles)
      ALLOCATE (atom_gradient(3, natom))
      atom_gradient = 0.0_dp
      CALL allocate_cell_gradient(ks_env, cell_gradient)
      CALL build_derivatives_basis(ks_env, left_sap, left_basis, left_neighbors)
      CALL build_derivatives_basis(ks_env, right_sap, right_basis, right_neighbors)
      CALL uspp_projector_create_onsite(left%overlaps, onsite)
      DO center = 1, SIZE(left%center_atoms)
         CALL dbcsr_get_block_p(onsite, center, center, block, found)
         IF (found) block(:, :) = left%augmentation(left%center_kinds(center))%q
      END DO
      DO image = 1, SIZE(weight)
         ALLOCATE (energy_weight(image)%matrix)
         CALL dbcsr_get_info(weight(image)%matrix, matrix_type=matrix_type)
         IF (matrix_type == 'S' .OR. matrix_type == 'A') THEN
            CALL dbcsr_desymmetrize(weight(image)%matrix, energy_weight(image)%matrix)
         ELSE
            CALL dbcsr_copy(energy_weight(image)%matrix, weight(image)%matrix)
         END IF
         CALL dbcsr_scale(energy_weight(image)%matrix, weight_factor)
      END DO
      IF (do_kpoints) THEN
         CALL get_ks_env(ks_env, kpoints=kpoints)
         CALL get_kpoint_info(kpoints, index_to_cell=cells)
         ALLOCATE (atom_kinds(natom))
         DO atom = 1, natom
            atom_kinds(atom) = particles(atom)%atomic_kind%kind_number
         END DO
         CPASSERT(SIZE(weight) <= SIZE(cells, 2))
         CALL qs_uspp_projector_mixed_images(left%overlaps, left_sap, right_sap, atom_kinds, &
                                             SIZE(left%augmentation), left%center_atoms, left%center_kinds, &
                                             cells(:, :SIZE(weight)), onsite, weight=energy_weight, &
                                             atom_gradient=atom_gradient, cell_gradient=cell_gradient)
      ELSE
         CPASSERT(SIZE(weight) == 1)
         CALL uspp_projector_mixed_pullback(left%overlaps, right%overlaps, onsite, &
                                            energy_weight(1)%matrix, left_gradient, right_gradient)
         CALL qs_uspp_projector_gradient(left_gradient, left_sap, left%center_atoms, atom_gradient, cell_gradient)
         CALL qs_uspp_projector_gradient(right_gradient, right_sap, right%center_atoms, atom_gradient, cell_gradient)
         CALL dbcsr_release(left_gradient)
         CALL dbcsr_release(right_gradient)
      END IF
      CALL accumulate_gradient(ks_env, atom_gradient, cell_gradient)
      CALL release_sap_int(left_sap)
      CALL release_sap_int(right_sap)
      CALL dbcsr_release(onsite)
      DO image = 1, SIZE(weight)
         CALL dbcsr_release(energy_weight(image)%matrix)
         DEALLOCATE (energy_weight(image)%matrix)
      END DO
   END SUBROUTINE qs_uspp_overlap_images

! **************************************************************************************************
!> \brief Gamma wrapper for an independent-basis overlap energy adjoint.
!> \param ks_env ...
!> \param left ...
!> \param right ...
!> \param left_basis ...
!> \param right_basis ...
!> \param left_neighbors ...
!> \param right_neighbors ...
!> \param weight ...
!> \param weight_factor ...
! **************************************************************************************************
   SUBROUTINE qs_uspp_overlap_single(ks_env, left, right, left_basis, right_basis, &
                                     left_neighbors, right_neighbors, weight, weight_factor)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(qs_uspp_type), INTENT(IN)                     :: left, right
      CHARACTER(LEN=*), INTENT(IN)                       :: left_basis, right_basis
      TYPE(neighbor_list_set_p_type), POINTER            :: left_neighbors(:), right_neighbors(:)
      TYPE(dbcsr_type), INTENT(IN), TARGET               :: weight
      REAL(dp), INTENT(IN)                               :: weight_factor

      TYPE(dbcsr_p_type)                                 :: weights(1)

      weights(1)%matrix => weight
      CALL qs_uspp_overlap_images(ks_env, left, right, left_basis, right_basis, left_neighbors, &
                                  right_neighbors, weights, weight_factor)
   END SUBROUTINE qs_uspp_overlap_single

! **************************************************************************************************
!> \brief Add unreduced USPP energy gradients to the normal Quickstep force accumulator.
!> \param ks_env Kinds and force arrays; Quickstep later reduces and negates these gradients.
!> \param atom_gradient Local energy gradients in global particle order.
!> \param cell_gradient Optional strain energy gradient; its negative is the virial.
! **************************************************************************************************
   SUBROUTINE accumulate_gradient(ks_env, atom_gradient, cell_gradient)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      REAL(KIND=dp), INTENT(IN)                          :: atom_gradient(:, :)
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: cell_gradient(3, 3)

      INTEGER                                            :: atom, ia, ikind
      INTEGER, ALLOCATABLE                               :: atom_of_kind(:), kind_of(:)
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind_set(:)
      TYPE(qs_force_type), POINTER                       :: force(:)
      TYPE(virial_type), POINTER                         :: virial

      CALL get_ks_env(ks_env, atomic_kind_set=atomic_kind_set, force=force)
      CALL get_atomic_kind_set(atomic_kind_set, atom_of_kind=atom_of_kind, kind_of=kind_of)
      DO atom = 1, SIZE(atom_gradient, 2)
         ikind = kind_of(atom)
         ia = atom_of_kind(atom)
         force(ikind)%gth_ppnl(:, ia) = force(ikind)%gth_ppnl(:, ia) + atom_gradient(:, atom)
      END DO
      IF (PRESENT(cell_gradient)) THEN
         CALL get_ks_env(ks_env, virial=virial)
         virial%pv_ppnl = virial%pv_ppnl - cell_gradient
         virial%pv_virial = virial%pv_virial - cell_gradient
      END IF
   END SUBROUTINE accumulate_gradient

! **************************************************************************************************
!> \brief Allocate the optional local strain accumulator only for analytic stress.
!> \param ks_env Current virial settings.
!> \param cell_gradient Allocated and zeroed for analytic stress, otherwise absent.
! **************************************************************************************************
   SUBROUTINE allocate_cell_gradient(ks_env, cell_gradient)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      REAL(KIND=dp), ALLOCATABLE, INTENT(OUT)            :: cell_gradient(:, :)

      TYPE(virial_type), POINTER                         :: virial

      CALL get_ks_env(ks_env, virial=virial)
      IF (virial%pv_availability .AND. .NOT. virial%pv_numer) THEN
         ALLOCATE (cell_gradient(3, 3))
         cell_gradient = 0.0_dp
      END IF
   END SUBROUTINE allocate_cell_gradient

! **************************************************************************************************
!> \brief Project an AO density with the Gamma matrix or the complete image convolution.
!> \param uspp Current geometry state.
!> \param density Spin density matrices over real-space images.
!> \param projected New distributed on-site density.
! **************************************************************************************************
   SUBROUTINE project_density(uspp, density, projected)
      TYPE(qs_uspp_type), INTENT(IN)                     :: uspp
      TYPE(dbcsr_p_type), INTENT(IN)                     :: density(:)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: projected

      IF (ALLOCATED(uspp%image_pairs)) THEN
         CALL qs_uspp_projector_images(uspp%overlaps, uspp%image_integrals, uspp%image_pairs, SIZE(uspp%augmentation), &
                                       uspp%center_atoms, uspp%center_kinds, density=density, projected=projected)
      ELSE
         CPASSERT(SIZE(density) == 1)
         CALL uspp_projector_density(uspp%overlaps, density(1)%matrix, projected)
      END IF
   END SUBROUTINE project_density

! **************************************************************************************************
!> \brief Collocate the projected density and its monopole using every real-space image.
!> \param uspp Current geometry state.
!> \param density Spin AO density matrices.
!> \param grid Physical density grid.
!> \param rho Reciprocal density, incremented in place.
!> \param charge Replicated total augmentation charge.
!> \param radial_fields Optional independent transitions in the same projector convention.
! **************************************************************************************************
   SUBROUTINE collocate_density(uspp, density, grid, rho, charge, radial_fields)
      TYPE(qs_uspp_type), INTENT(INOUT), TARGET          :: uspp
      TYPE(dbcsr_p_type), INTENT(IN)                     :: density(:)
      TYPE(pw_grid_type), INTENT(IN)                     :: grid
      COMPLEX(KIND=dp), INTENT(INOUT)                    :: rho(:)
      REAL(KIND=dp), INTENT(OUT)                         :: charge
      TYPE(uspp_augmentation_type), INTENT(INOUT), &
         OPTIONAL, TARGET                                :: radial_fields(:)

      INTEGER                                            :: center, ikind, n
      LOGICAL                                            :: found
      REAL(KIND=dp), ALLOCATABLE                         :: matrix(:, :)
      REAL(KIND=dp), POINTER                             :: block(:, :)
      TYPE(dbcsr_type)                                   :: projected
      TYPE(mp_comm_type)                                 :: group
      TYPE(uspp_augmentation_type), POINTER              :: cache(:)

      cache => uspp%augmentation
      IF (PRESENT(radial_fields)) THEN
         CPASSERT(SIZE(radial_fields) == SIZE(uspp%augmentation))
         cache => radial_fields
      END IF
      CALL dbcsr_get_info(uspp%overlaps, group=group)
      CALL project_density(uspp, density, projected)
      DO ikind = 1, SIZE(cache)
         IF (cache(ikind)%nproj == 0) CYCLE
         CALL uspp_augmentation_prepare_fourier(cache(ikind), grid%g)
      END DO
      charge = 0.0_dp
      DO center = 1, SIZE(uspp%center_atoms)
         n = uspp%projector_sizes(center)
         ikind = uspp%center_kinds(center)
         IF (cache(ikind)%nproj == 0) CYCLE
         ALLOCATE (matrix(n, n))
         matrix = 0.0_dp
         CALL dbcsr_get_block_p(projected, center, center, block, found)
         IF (found) matrix(:, :) = block
         CALL group%sum(matrix)
         CALL uspp_augmentation_collocate(cache(ikind), grid%g, uspp%position(:, center), grid%vol, matrix, rho)
         charge = charge + SUM(matrix*cache(ikind)%q)
         DEALLOCATE (matrix)
      END DO
      CALL dbcsr_release(projected)
   END SUBROUTINE collocate_density

! **************************************************************************************************
!> \brief Apply a real on-site operator and optionally its projector energy derivatives.
!> \param uspp Current geometry state.
!> \param onsite Distributed on-site operator.
!> \param hamiltonian Optional AO image matrices to increment.
!> \param density Optional AO image densities for derivatives.
!> \param sap_int Optional sorted image-resolved derivatives.
!> \param atom_gradient Optional local atomic energy gradients.
!> \param cell_gradient Optional local strain energy gradient.
!> \param alpha Optional operator/gradient prefactor, e.g. -1 for metric derivatives.
! **************************************************************************************************
   SUBROUTINE apply_operator(uspp, onsite, hamiltonian, density, sap_int, atom_gradient, cell_gradient, alpha)
      TYPE(qs_uspp_type), INTENT(IN)                     :: uspp
      TYPE(dbcsr_type), INTENT(IN)                       :: onsite
      TYPE(dbcsr_p_type), INTENT(INOUT), OPTIONAL        :: hamiltonian(:)
      TYPE(dbcsr_p_type), INTENT(IN), OPTIONAL           :: density(:)
      TYPE(sap_int_type), INTENT(IN), OPTIONAL           :: sap_int(:)
      REAL(KIND=dp), INTENT(INOUT), OPTIONAL             :: atom_gradient(:, :), cell_gradient(3, 3)
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: alpha

      REAL(KIND=dp)                                      :: scale
      TYPE(dbcsr_type)                                   :: gradient, scaled

      scale = 1.0_dp
      IF (PRESENT(alpha)) scale = alpha
      CPASSERT(PRESENT(density) .EQV. PRESENT(sap_int))
      CPASSERT(PRESENT(density) .EQV. PRESENT(atom_gradient))
      IF (ALLOCATED(uspp%image_pairs)) THEN
         IF (PRESENT(density)) THEN
            CALL qs_uspp_projector_images(uspp%overlaps, sap_int, uspp%image_pairs, SIZE(uspp%augmentation), &
                                          uspp%center_atoms, uspp%center_kinds, onsite, hamiltonian, density, &
                                          atom_gradient=atom_gradient, cell_gradient=cell_gradient, alpha=scale)
         ELSE
            CALL qs_uspp_projector_images(uspp%overlaps, uspp%image_integrals, uspp%image_pairs, SIZE(uspp%augmentation), &
                                          uspp%center_atoms, uspp%center_kinds, onsite, hamiltonian, alpha=scale)
         END IF
      ELSE
         IF (PRESENT(hamiltonian)) THEN
            CPASSERT(SIZE(hamiltonian) == 1)
            CALL dbcsr_copy(scaled, onsite)
            CALL dbcsr_scale(scaled, scale)
            CALL uspp_projector_add_operator(uspp%overlaps, scaled, hamiltonian(1)%matrix)
            CALL dbcsr_release(scaled)
         END IF
         IF (PRESENT(density)) THEN
            CPASSERT(SIZE(density) == 1)
            CALL uspp_projector_pullback(uspp%overlaps, density(1)%matrix, onsite, scale, gradient)
            CALL qs_uspp_projector_gradient(gradient, sap_int, uspp%center_atoms, atom_gradient, cell_gradient)
            CALL dbcsr_release(gradient)
         END IF
      END IF
   END SUBROUTINE apply_operator

! **************************************************************************************************
!> \brief Scalar adapter to the image-aware density path.
!> \param ks_env Current KS environment.
!> \param density Gamma AO density.
!> \param rho_g Reciprocal density.
!> \param rho_r Real density.
!> \param total_rho Optional electronic charge.
!> \param derivative Optional derivative orders.
!> \param kinetic Select independent kinetic transitions instead of charge augmentation.
! **************************************************************************************************
   SUBROUTINE qs_uspp_collocate_single(ks_env, density, rho_g, rho_r, total_rho, derivative, kinetic)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(dbcsr_type), INTENT(IN), TARGET               :: density
      TYPE(pw_c1d_gs_type), INTENT(INOUT)                :: rho_g
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: rho_r
      REAL(KIND=dp), INTENT(INOUT), OPTIONAL             :: total_rho
      INTEGER, INTENT(IN), OPTIONAL                      :: derivative(3)
      LOGICAL, INTENT(IN), OPTIONAL                      :: kinetic

      TYPE(dbcsr_p_type)                                 :: matrices(1)

      matrices(1)%matrix => density
      CALL qs_uspp_collocate_images(ks_env, matrices, rho_g, rho_r, total_rho, derivative, kinetic)
   END SUBROUTINE qs_uspp_collocate_single

! **************************************************************************************************
!> \brief Scalar adapter to the image-aware grid adjoint.
!> \param ks_env Current KS environment.
!> \param potential Real grid potential.
!> \param hamiltonian Gamma Hamiltonian to increment.
!> \param weighted Whether potential includes dV.
!> \param density Optional Gamma force-path density.
!> \param kinetic Select independent kinetic transitions instead of charge augmentation.
! **************************************************************************************************
   SUBROUTINE qs_uspp_integrate_single(ks_env, potential, hamiltonian, weighted, density, kinetic)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: potential
      TYPE(dbcsr_type), INTENT(INOUT), TARGET            :: hamiltonian
      LOGICAL, INTENT(IN)                                :: weighted
      TYPE(dbcsr_type), INTENT(IN), OPTIONAL, TARGET     :: density
      LOGICAL, INTENT(IN), OPTIONAL                      :: kinetic

      TYPE(dbcsr_p_type)                                 :: matrices(1)
      TYPE(dbcsr_p_type), ALLOCATABLE                    :: densities(:)

      matrices(1)%matrix => hamiltonian
      IF (PRESENT(density)) THEN
         ALLOCATE (densities(1))
         densities(1)%matrix => density
      END IF
      CALL qs_uspp_integrate_images(ks_env, potential, matrices, weighted, densities, kinetic)
   END SUBROUTINE qs_uspp_integrate_single

! **************************************************************************************************
!> \brief Scalar adapter to the image-aware ionic operator and metric derivatives.
!> \param uspp Current USPP state.
!> \param qs_kind_set Current kinds.
!> \param particle_set Current atoms.
!> \param grid Physical density grid.
!> \param hamiltonian Gamma core Hamiltonian.
!> \param density Optional total Gamma density.
!> \param energy_weighted Optional total Gamma energy-weighted density.
!> \param sap_int Optional projector derivatives.
!> \param atom_gradient Optional local atomic energy gradients.
!> \param cell_gradient Optional local strain energy gradient.
! **************************************************************************************************
   SUBROUTINE qs_uspp_core_matrix_single(uspp, qs_kind_set, particle_set, grid, hamiltonian, &
                                         density, energy_weighted, sap_int, atom_gradient, cell_gradient)
      TYPE(qs_uspp_type), INTENT(INOUT)                  :: uspp
      TYPE(qs_kind_type), INTENT(INOUT)                  :: qs_kind_set(:)
      TYPE(particle_type), INTENT(IN)                    :: particle_set(:)
      TYPE(pw_grid_type), POINTER                        :: grid
      TYPE(dbcsr_type), INTENT(INOUT), TARGET            :: hamiltonian
      TYPE(dbcsr_type), INTENT(IN), OPTIONAL, TARGET     :: density, energy_weighted
      TYPE(sap_int_type), INTENT(IN), OPTIONAL           :: sap_int(:)
      REAL(KIND=dp), INTENT(INOUT), OPTIONAL             :: atom_gradient(:, :), cell_gradient(3, 3)

      TYPE(dbcsr_p_type)                                 :: matrices(1)
      TYPE(dbcsr_p_type), ALLOCATABLE                    :: densities(:), weighted(:)

      matrices(1)%matrix => hamiltonian
      IF (PRESENT(density)) THEN
         ALLOCATE (densities(1))
         densities(1)%matrix => density
      END IF
      IF (PRESENT(energy_weighted)) THEN
         ALLOCATE (weighted(1))
         weighted(1)%matrix => energy_weighted
      END IF
      CALL qs_uspp_core_matrix_images(uspp, qs_kind_set, particle_set, grid, matrices, &
                                      densities, weighted, sap_int, atom_gradient, cell_gradient)
   END SUBROUTINE qs_uspp_core_matrix_single

! **************************************************************************************************
!> \brief Require every retained G to occur once over the projector communicator.
!> \param uspp Current projector state.
!> \param grid Quickstep FFT grid; a reordered communicator with the same members is allowed.
! **************************************************************************************************
   SUBROUTINE check_grid(uspp, grid)
      TYPE(qs_uspp_type), INTENT(IN)                     :: uspp
      TYPE(pw_grid_type), INTENT(IN)                     :: grid

      TYPE(mp_comm_type)                                 :: group

      CALL dbcsr_get_info(uspp%overlaps, group=group)
      CPASSERT(group%compare(grid%para%group) /= mp_comm_unequal)
      CPASSERT(grid%grid_span == FULLSPACE .OR. grid%grid_span == HALFSPACE)
      IF (grid%grid_span == HALFSPACE) THEN
         CPASSERT(ALL(MODULO(grid%npts, 2) == 1))
      END IF
   END SUBROUTINE check_grid

! **************************************************************************************************
!> \brief Whether a nonzero verified UPF model-core kinetic field is present.
!> \param ks_env Current USPP and PAW caches.
!> \return True for at least one nonzero auxiliary PS core field.
! **************************************************************************************************
   FUNCTION qs_paw_core_tau_present(ks_env) RESULT(present_core)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      LOGICAL                                            :: present_core

      INTEGER                                            :: ikind
      TYPE(paw_core_kinetic_type), POINTER               :: cache
      TYPE(qs_uspp_type), POINTER                        :: uspp

      present_core = .FALSE.
      CALL get_ks_env(ks_env, uspp=uspp)
      IF (.NOT. ASSOCIATED(uspp)) RETURN
      DO ikind = 1, SIZE(uspp%paw)
         cache => core_tau_cache(uspp, ikind)
         IF (.NOT. ASSOCIATED(cache)) CYCLE
         IF (.NOT. cache%available) CYCLE
         IF (ANY(cache%weighted /= 0.0_dp)) present_core = .TRUE.
      END DO
   END FUNCTION qs_paw_core_tau_present

! **************************************************************************************************
!> \brief Spin-summed periodic PS core kinetic field on the requested XC grid.
!>        The R-to-G adjoint of the real projection is used by the force routine.
!> \param ks_env Current PAW centers and verified radial caches.
!> \param field Allocated real PW field, overwritten without changing any valence density.
! **************************************************************************************************
   SUBROUTINE qs_paw_core_tau_field(ks_env, field)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: field

      COMPLEX(dp)                                        :: coefficient
      INTEGER                                            :: center, ig, ikind
      TYPE(paw_core_kinetic_type), POINTER               :: cache
      TYPE(pw_c1d_gs_type)                               :: core_g
      TYPE(qs_uspp_type), POINTER                        :: uspp

      CALL get_ks_env(ks_env, uspp=uspp)
      CPASSERT(ASSOCIATED(uspp))
      CALL check_grid(uspp, field%pw_grid)
      CALL core_g%create(field%pw_grid)
      core_g%array = CMPLX(0.0_dp, 0.0_dp, dp)
      DO center = 1, SIZE(uspp%center_atoms)
         ikind = uspp%center_kinds(center)
         cache => core_tau_cache(uspp, ikind)
         IF (.NOT. ASSOCIATED(cache)) CYCLE
         IF (.NOT. cache%available) CYCLE
         IF (.NOT. ANY(cache%weighted /= 0.0_dp)) CYCLE
         DO ig = 1, SIZE(core_g%array)
            CALL paw_core_kinetic_fourier(cache, field%pw_grid%g(:, ig), uspp%position(:, center), &
                                          field%pw_grid%vol, coefficient)
            core_g%array(ig) = core_g%array(ig) + coefficient
         END DO
      END DO
      CALL pw_transfer(core_g, field)
      CALL core_g%release()
   END SUBROUTINE qs_paw_core_tau_field

! **************************************************************************************************
!> \brief Center and affine-strain adjoints of the smooth PW core kinetic field.
!>        This fixed core field has no AO density-matrix derivative.
!> \param ks_env Current centers and unreduced force/virial accumulators.
!> \param potential Unweighted derivative of XC energy density with respect to spin tau.
!> \param spin_fraction Fraction of the spin-summed closed core in this channel.
! **************************************************************************************************
   SUBROUTINE qs_paw_core_tau_forces(ks_env, potential, spin_fraction)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: potential
      REAL(dp), INTENT(IN)                               :: spin_fraction

      COMPLEX(dp)                                        :: coefficient, derivative(3), &
                                                            strain(3, 3), weighted
      INTEGER                                            :: center, ig, ikind, natom
      REAL(dp)                                           :: multiplicity
      REAL(dp), ALLOCATABLE                              :: atom_gradient(:, :), cell_gradient(:, :)
      TYPE(paw_core_kinetic_type), POINTER               :: cache
      TYPE(pw_c1d_gs_type)                               :: potential_g
      TYPE(qs_uspp_type), POINTER                        :: uspp

      IF (.NOT. qs_paw_core_tau_present(ks_env)) RETURN
      CALL get_ks_env(ks_env, uspp=uspp, natom=natom)
      CALL check_grid(uspp, potential%pw_grid)
      ALLOCATE (atom_gradient(3, natom), SOURCE=0.0_dp)
      CALL allocate_cell_gradient(ks_env, cell_gradient)
      CALL potential_g%create(potential%pw_grid)
      CALL pw_transfer(potential, potential_g)
      DO center = 1, SIZE(uspp%center_atoms)
         ikind = uspp%center_kinds(center)
         cache => core_tau_cache(uspp, ikind)
         IF (.NOT. ASSOCIATED(cache)) CYCLE
         IF (.NOT. cache%available) CYCLE
         IF (.NOT. ANY(cache%weighted /= 0.0_dp)) CYCLE
         DO ig = 1, SIZE(potential_g%array)
            multiplicity = 1.0_dp
            IF (potential%pw_grid%grid_span == HALFSPACE) THEN
               multiplicity = 2.0_dp
               IF (potential%pw_grid%have_g0 .AND. ig == 1) multiplicity = 1.0_dp
            END IF
            weighted = spin_fraction*potential%pw_grid%vol*multiplicity*potential_g%array(ig)
            IF (ALLOCATED(cell_gradient)) THEN
               CALL paw_core_kinetic_fourier(cache, potential%pw_grid%g(:, ig), uspp%position(:, center), &
                                             potential%pw_grid%vol, coefficient, strain, derivative)
               cell_gradient(:, :) = cell_gradient + REAL(CONJG(strain)*weighted, dp)
            ELSE
               CALL paw_core_kinetic_fourier(cache, potential%pw_grid%g(:, ig), uspp%position(:, center), &
                                             potential%pw_grid%vol, coefficient, center_derivative=derivative)
            END IF
            atom_gradient(:, uspp%center_atoms(center)) = atom_gradient(:, uspp%center_atoms(center)) + &
                                                          REAL(CONJG(derivative)*weighted, dp)
         END DO
      END DO
      CALL accumulate_gradient(ks_env, atom_gradient, cell_gradient)
      CALL potential_g%release()
   END SUBROUTINE qs_paw_core_tau_forces

! **************************************************************************************************
!> \brief Add spin-summed PS core tau from all centers and auxiliary-cell images.
!> \param uspp Current PAW centers and original-data kinetic caches.
!> \param cell Auxiliary periodic cell, as in augmentation collocation.
!> \param points Unfolded Cartesian atom-grid points.
!> \param kinetic Spin-summed core tau values, incremented without quadrature weights.
! **************************************************************************************************
   SUBROUTINE qs_paw_core_tau_atom_collocate(uspp, cell, points, kinetic)
      TYPE(qs_uspp_type), INTENT(IN), TARGET             :: uspp
      TYPE(cell_type), INTENT(IN)                        :: cell
      REAL(dp), INTENT(IN)                               :: points(:, :)
      REAL(dp), INTENT(INOUT)                            :: kinetic(:)

      INTEGER                                            :: center, ikind, ix, iy, iz, lower(3), &
                                                            point, upper(3)
      REAL(dp)                                           :: image_position(3), value
      REAL(dp), ALLOCATABLE                              :: folded(:, :)
      TYPE(paw_core_kinetic_type), POINTER               :: cache

      CPASSERT(ALL(SHAPE(points) == [3, SIZE(kinetic)]))
      IF (SIZE(kinetic) == 0) RETURN
      ALLOCATE (folded, SOURCE=points - MATMUL(cell%hmat, ANINT(MATMUL(cell%h_inv, points))))
      DO center = 1, SIZE(uspp%center_atoms)
         ikind = uspp%center_kinds(center)
         cache => core_tau_cache(uspp, ikind)
         IF (.NOT. ASSOCIATED(cache)) CYCLE
         IF (.NOT. cache%available) CYCLE
         IF (.NOT. ANY(cache%weighted /= 0.0_dp)) CYCLE
         CALL atom_grid_image_bounds(cell, uspp%position(:, center), folded, cache%support_radius, lower, upper)
         DO iz = lower(3), upper(3)
            DO iy = lower(2), upper(2)
               DO ix = lower(1), upper(1)
                  image_position = uspp%position(:, center) + MATMUL(cell%hmat, REAL([ix, iy, iz], dp))
                  DO point = 1, SIZE(kinetic)
                     CALL paw_core_kinetic_cartesian(cache, folded(:, point) - image_position, value)
                     kinetic(point) = kinetic(point) + value
                  END DO
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE qs_paw_core_tau_atom_collocate

! **************************************************************************************************
!> \brief Exact center and point adjoints of periodic atom-grid PS core tau collocation.
!> \param uspp Same PAW state as the forward map.
!> \param fields Accumulators for source-center and source-image strain derivatives.
!> \param cell Auxiliary cell used by the forward map.
!> \param points Same unfolded quadrature points.
!> \param potential Weighted derivative with respect to spin-summed core tau, including spin fractions.
!> \param point_gradient Weighted point energy gradients, incremented; caller handles target-grid motion.
! **************************************************************************************************
   SUBROUTINE qs_paw_core_tau_atom_integrate(uspp, fields, cell, points, potential, point_gradient)
      TYPE(qs_uspp_type), INTENT(IN), TARGET             :: uspp
      TYPE(qs_uspp_atom_grid_type), INTENT(INOUT)        :: fields
      TYPE(cell_type), INTENT(IN)                        :: cell
      REAL(dp), INTENT(IN)                               :: points(:, :), potential(:)
      REAL(dp), INTENT(INOUT)                            :: point_gradient(:, :)

      INTEGER                                            :: center, ikind, ix, iy, iz, jdir, &
                                                            lower(3), point, upper(3)
      REAL(dp)                                           :: derivative(3), image_position(3), &
                                                            spatial(3), value
      REAL(dp), ALLOCATABLE                              :: folded(:, :), images(:, :), &
                                                            point_derivative(:, :)
      TYPE(paw_core_kinetic_type), POINTER               :: cache

      CPASSERT(fields%calculate_forces)
      CPASSERT(ALL(SHAPE(points) == [3, SIZE(potential)]))
      CPASSERT(ALL(SHAPE(point_gradient) == SHAPE(points)))
      IF (SIZE(potential) == 0) RETURN
      ALLOCATE (images, SOURCE=MATMUL(cell%hmat, ANINT(MATMUL(cell%h_inv, points))))
      ALLOCATE (folded, SOURCE=points - images)
      ALLOCATE (point_derivative(3, SIZE(potential)))
      DO center = 1, SIZE(uspp%center_atoms)
         ikind = uspp%center_kinds(center)
         cache => core_tau_cache(uspp, ikind)
         IF (.NOT. ASSOCIATED(cache)) CYCLE
         IF (.NOT. cache%available) CYCLE
         IF (.NOT. ANY(cache%weighted /= 0.0_dp)) CYCLE
         CALL atom_grid_image_bounds(cell, uspp%position(:, center), folded, cache%support_radius, lower, upper)
         DO iz = lower(3), upper(3)
            DO iy = lower(2), upper(2)
               DO ix = lower(1), upper(1)
                  image_position = uspp%position(:, center) + MATMUL(cell%hmat, REAL([ix, iy, iz], dp))
                  DO point = 1, SIZE(potential)
                     CALL paw_core_kinetic_cartesian(cache, folded(:, point) - image_position, value, spatial)
                     point_derivative(:, point) = potential(point)*spatial
                  END DO
                  derivative = -SUM(point_derivative, DIM=2)
                  fields%center_gradient(:, center) = fields%center_gradient(:, center) + derivative
                  point_gradient = point_gradient + point_derivative
                  fields%center_strain = fields%center_strain - MATMUL(point_derivative, TRANSPOSE(images))
                  DO jdir = 1, 3
                     fields%center_strain(:, jdir) = fields%center_strain(:, jdir) + derivative*image_position(jdir)
                  END DO
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE qs_paw_core_tau_atom_integrate

! **************************************************************************************************
!> \brief PAW exchange energy, AO operator and projector force/strain derivatives.
!>        The fixed physical radial interaction has no geometry response. Both projector-overlap legs
!>        are differentiated through the same projected density as used for the energy.
!>        Verified original closed-core orbitals supply core-valence exchange and the core-core reference.
!> \param ks_env Prepared PAW caches, projector overlaps and force accumulators.
!> \param density Spin and image AO densities, including occupations and k-point weights.
!>        One Quickstep spin block is a restricted spin sum; two blocks are physical spin densities.
!> \param fraction Exchange fraction; Coulomb replicas may be summed, finite radii remain separate.
!> \param exchange PAW exchange correction including the closed-core terms, reset and returned replicated.
!> \param calculate_forces Whether to differentiate projectors at fixed AO density.
!>        The generalized-overlap Pulay contribution is accumulated separately.
!> \param hamiltonian Optional spin and image AO matrices to increment.
!> \param cutoff_radius Optional truncated interaction radius.
!> \param omega Optional inverse erfc screening length; mutually exclusive with cutoff_radius.
!> \param reconstructed_valence Add only closed-core exchange when RI already reconstructs the valence density.
! **************************************************************************************************
   SUBROUTINE qs_paw_exchange(ks_env, density, fraction, exchange, calculate_forces, hamiltonian, cutoff_radius, omega, &
                              reconstructed_valence)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(dbcsr_p_type), INTENT(IN)                     :: density(:, :)
      REAL(dp), INTENT(IN)                               :: fraction
      REAL(dp), INTENT(OUT)                              :: exchange
      LOGICAL, INTENT(IN)                                :: calculate_forces
      TYPE(dbcsr_p_type), INTENT(INOUT), OPTIONAL        :: hamiltonian(:, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: cutoff_radius, omega
      LOGICAL, INTENT(IN), OPTIONAL                      :: reconstructed_valence

      INTEGER                                            :: center, handle, ikind, n, natom, ns, &
                                                            nuclear_charge, spin
      INTEGER, ALLOCATABLE                               :: core_cache(:)
      LOGICAL                                            :: core_only, found, frozen_core, &
                                                            valid_shape
      REAL(dp)                                           :: center_exchange
      REAL(dp), ALLOCATABLE                              :: atom_gradient(:, :), &
                                                            cell_gradient(:, :), dm(:, :, :), &
                                                            potential(:, :, :)
      REAL(dp), POINTER                                  :: block(:, :)
      TYPE(dbcsr_type), ALLOCATABLE                      :: onsite(:), projected(:)
      TYPE(mp_comm_type)                                 :: group
      TYPE(qs_kind_type), POINTER                        :: qs_kind_set(:)
      TYPE(qs_uspp_type), POINTER                        :: uspp
      TYPE(sap_int_type), POINTER                        :: sap_int(:)

      core_only = PRESENT(cutoff_radius) .OR. PRESENT(omega)
      IF (PRESENT(reconstructed_valence)) core_only = core_only .OR. reconstructed_valence
      CPASSERT(.NOT. (PRESENT(cutoff_radius) .AND. PRESENT(omega)))
      exchange = 0.0_dp
      CPASSERT(ieee_is_finite(fraction))
      IF (fraction == 0.0_dp) RETURN
      CALL get_ks_env(ks_env, uspp=uspp, qs_kind_set=qs_kind_set)
      IF (.NOT. ASSOCIATED(uspp)) RETURN
      IF (.NOT. ALLOCATED(uspp%paw)) RETURN
      IF (.NOT. ANY(uspp%paw%active)) RETURN
      ALLOCATE (core_cache(SIZE(uspp%paw)), SOURCE=0)
      ! Prepare the immutable original-data core contribution once per PAW kind.
      DO ikind = 1, SIZE(uspp%paw)
         IF (.NOT. uspp%paw(ikind)%active) CYCLE
         CALL get_qs_kind(qs_kind_set(ikind), zatom=nuclear_charge)
         ASSOCIATE (pot => qs_kind_set(ikind)%sgp_potential%upf_data)
            CPASSERT(ALLOCATED(pot%paw))
            CPASSERT(ieee_is_finite(pot%zion))
            frozen_core = ABS(REAL(nuclear_charge, dp) - pot%zion) > &
                          64.0_dp*EPSILON(1.0_dp)*REAL(nuclear_charge, dp)
            frozen_core = frozen_core .OR. ANY(pot%paw%rho_core /= 0.0_dp)
            IF (ALLOCATED(pot%rho_nlcc)) THEN
               frozen_core = frozen_core .OR. ANY(pot%rho_nlcc /= 0.0_dp)
            END IF
            IF (PRESENT(cutoff_radius) .OR. PRESENT(omega)) THEN
               CALL paw_onecenter_core_exchange_prepare(pot, qs_kind_set(ikind)%sgp_potential%upf_projectors, &
                                                        uspp%paw(ikind)%onecenter, cutoff_radius, core_cache(ikind), omega)
               ASSOCIATE (core => uspp%paw(ikind)%onecenter%interaction_core_exchange(core_cache(ikind)))
                  IF (.NOT. core%available) THEN
                     CALL cp_abort(__LOCATION__, 'PAW frozen-core exact exchange: '//TRIM(core%reason))
                  END IF
               END ASSOCIATE
            ELSE IF (frozen_core .OR. core_only) THEN
               IF (.NOT. uspp%paw(ikind)%onecenter%core_exchange_initialized) THEN
                  CALL paw_onecenter_core_exchange_init(pot, qs_kind_set(ikind)%sgp_potential%upf_projectors, &
                                                        uspp%paw(ikind)%onecenter)
               END IF
               IF (.NOT. uspp%paw(ikind)%onecenter%core_exchange_available) THEN
                  CALL cp_abort(__LOCATION__, &
                                'PAW frozen-core exact exchange: '//TRIM(uspp%paw(ikind)%onecenter%core_exchange_reason))
               END IF
            END IF
         END ASSOCIATE
      END DO
      CALL timeset('qs_paw_exchange', handle)
      ns = SIZE(density, 1)
      CPASSERT(ns == 1 .OR. ns == 2)
      IF (PRESENT(hamiltonian)) THEN
         valid_shape = ALL(SHAPE(hamiltonian) == SHAPE(density))
         CPASSERT(valid_shape)
      END IF
      CALL dbcsr_get_info(uspp%overlaps, group=group)
      ALLOCATE (projected(ns), onsite(ns))
      DO spin = 1, ns
         CALL project_density(uspp, density(spin, :), projected(spin))
         CALL uspp_projector_create_onsite(uspp%overlaps, onsite(spin))
      END DO
      DO center = 1, SIZE(uspp%center_atoms)
         ikind = uspp%center_kinds(center)
         IF (.NOT. uspp%paw(ikind)%active) CYCLE
         ! All spin blocks have the same owner; count each center exactly once.
         CALL dbcsr_get_block_p(projected(1), center, center, block, found)
         IF (.NOT. found) CYCLE
         n = uspp%projector_sizes(center)
         ALLOCATE (dm(n, n, 2), potential(n, n, 2))
         IF (ns == 1) THEN
            dm(:, :, 1) = 0.5_dp*block
            dm(:, :, 2) = dm(:, :, 1)
         ELSE
            DO spin = 1, ns
               CALL dbcsr_get_block_p(projected(spin), center, center, block, found)
               CPASSERT(found)
               dm(:, :, spin) = block
            END DO
         END IF
         IF (core_only) THEN
            CALL paw_onecenter_core_exchange(uspp%paw(ikind)%onecenter, dm, fraction, core_cache(ikind), &
                                             center_exchange, potential)
         ELSE
            CALL paw_onecenter_exchange(uspp%paw(ikind)%onecenter, dm, fraction, center_exchange, potential)
         END IF
         exchange = exchange + center_exchange
         DO spin = 1, ns
            CALL dbcsr_get_block_p(onsite(spin), center, center, block, found)
            CPASSERT(found)
            IF (ns == 1) THEN
               ! dE/dDtotal=(Vup+Vdown)/2 because Dspin=Dtotal/2.
               block(:, :) = 0.5_dp*(potential(:, :, 1) + potential(:, :, 2))
            ELSE
               block(:, :) = potential(:, :, spin)
            END IF
         END DO
         DEALLOCATE (dm, potential)
      END DO
      CALL group%sum(exchange)
      NULLIFY (sap_int)
      IF (calculate_forces) THEN
         CALL get_ks_env(ks_env, natom=natom)
         ALLOCATE (atom_gradient(3, natom))
         atom_gradient(:, :) = 0.0_dp
         CALL allocate_cell_gradient(ks_env, cell_gradient)
         CALL build_derivatives(ks_env, sap_int)
      END IF
      DO spin = 1, ns
         IF (calculate_forces) THEN
            IF (PRESENT(hamiltonian)) THEN
               CALL apply_operator(uspp, onsite(spin), hamiltonian(spin, :), density(spin, :), &
                                   sap_int, atom_gradient, cell_gradient)
            ELSE
               CALL apply_operator(uspp, onsite(spin), density=density(spin, :), &
                                   sap_int=sap_int, atom_gradient=atom_gradient, cell_gradient=cell_gradient)
            END IF
         ELSE IF (PRESENT(hamiltonian)) THEN
            CALL apply_operator(uspp, onsite(spin), hamiltonian(spin, :))
         END IF
         CALL dbcsr_release(projected(spin))
         CALL dbcsr_release(onsite(spin))
      END DO
      IF (calculate_forces) THEN
         CALL accumulate_gradient(ks_env, atom_gradient, cell_gradient)
         CALL release_sap_int(sap_int)
      END IF
      CALL timestop(handle)
   END SUBROUTINE qs_paw_exchange

END MODULE qs_uspp_methods
