!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
!> \brief Original-coordinate auxiliary operators and physical projection/charge adjoints.
MODULE qs_upf_inverse_aux_operator
   USE admm_types,                      ONLY: admm_type,&
                                              get_admm_env
   USE cell_types,                      ONLY: cell_type
   USE cp_cfm_types,                    ONLY: cp_cfm_get_submatrix
   USE cp_dbcsr_api,                    ONLY: dbcsr_copy,&
                                              dbcsr_p_type,&
                                              dbcsr_set
   USE cp_dbcsr_operations,             ONLY: dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_fm_types,                     ONLY: cp_fm_get_submatrix
   USE hfx_ri,                          ONLY: hfx_ri_update_forces,&
                                              hfx_ri_update_ks
   USE hfx_ri_kp,                       ONLY: hfx_ri_update_forces_kp,&
                                              hfx_ri_update_ks_kp
   USE hfx_types,                       ONLY: hfx_type
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE input_constants,                 ONLY: do_admm_exch_scaling_none,&
                                              do_admm_purify_none
   USE input_cp2k_hfx,                  ONLY: ri_pmat
   USE input_section_types,             ONLY: section_vals_type
   USE integral_library_types,          ONLY: compare_potential_types
   USE kinds,                           ONLY: dp
   USE kpoint_methods,                  ONLY: kpoint_star_density,&
                                              kpoint_star_operator
   USE kpoint_transitional,             ONLY: get_1d_pointer
   USE kpoint_types,                    ONLY: kpoint_env_type,&
                                              kpoint_type
   USE mathlib,                         ONLY: invert_matrix,&
                                              invmat
   USE particle_types,                  ONLY: particle_type
   USE paw_types,                       ONLY: paw_potential_type
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_upf_inverse_aux_xc,           ONLY: qs_upf_inverse_aux_xc_evaluate
   USE qs_upf_inverse_background,       ONLY: qs_upf_inverse_background_collocate,&
                                              qs_upf_inverse_background_geometry
   USE qs_upf_inverse_charge_metric,    ONLY: qs_upf_inverse_charge_metric_contract
   USE qs_upf_inverse_coulomb,          ONLY: qs_upf_inverse_coulomb_metric_response
   USE qs_upf_inverse_grid,             ONLY: qs_upf_inverse_grid_metric_dense,&
                                              qs_upf_inverse_matrix_gather
   USE qs_upf_inverse_kpoints,          ONLY: qs_upf_inverse_kpoints_density,&
                                              qs_upf_inverse_kpoints_original_operator
   USE qs_upf_inverse_paw,              ONLY: qs_upf_inverse_paw_core_prepare,&
                                              qs_upf_inverse_paw_evaluate,&
                                              qs_upf_inverse_paw_response
   USE qs_upf_inverse_types,            ONLY: qs_upf_inverse_grid_type,&
                                              qs_upf_inverse_scf_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: qs_upf_inverse_aux_projection, qs_upf_inverse_aux_operator_evaluate
   PUBLIC :: qs_upf_inverse_aux_metric_response
   PUBLIC :: qs_upf_inverse_aux_core_exchange, qs_upf_inverse_aux_response
   PUBLIC :: qs_upf_inverse_aux_evaluate

CONTAINS

! **************************************************************************************************
!> \brief Evaluate the complete auxiliary exchange correction using the actual native RI owners.
!>        Each HF repetition contributes once, irrespective of its thread replicas. RI uses
!>        independent temporary output matrices and retains its normal density/integral history.
!>        Core HF, selected XC, projection and physical charge normalization use that same density.
!> \param qs_env Prepared native projection, auxiliary occupations and RI basis owners.
!> \param xc_section Selected signed auxiliary XC correction.
!> \param exchange Actual native HF repetitions and their initialized AUX_FIT RI containers.
!> \param geometry_changed Rebuild each RI geometry when required by the caller.
!> \param components RI valence HF, PAW core HF, grid XC and atomic XC energies, in that order.
!> \param operators Complete primary AO derivative in original physical k/spin order, without weights.
!> \param atom_response Optional total atomic displacement, including affine cell motion.
!> \param cell_response Optional lattice derivative at fixed reduced k.
!> \param response Optional complete energy derivative. Supply all three geometry arguments together.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_aux_evaluate(qs_env, xc_section, exchange, geometry_changed, components, operators, &
                                          atom_response, cell_response, response)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: xc_section
      TYPE(hfx_type), INTENT(IN)                         :: exchange(:, :)
      LOGICAL, INTENT(IN)                                :: geometry_changed
      REAL(dp), INTENT(OUT)                              :: components(4)
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: operators(:, :, :, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: atom_response(:, :), cell_response(3, 3)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: response

      COMPLEX(dp), ALLOCATABLE                           :: ri_operator(:, :, :, :)
      LOGICAL                                            :: geometry, use_virial
      REAL(dp)                                           :: energy, ri_energy, strain_gradient(3, 3)
      REAL(dp), ALLOCATABLE                              :: atom_gradient(:, :)

      geometry = PRESENT(response)
      CPASSERT(geometry .EQV. PRESENT(atom_response))
      CPASSERT(geometry .EQV. PRESENT(cell_response))
      use_virial = .FALSE.
      IF (geometry) use_virial = ANY(cell_response /= 0.0_dp)
      CALL auxiliary_ri(qs_env, exchange, geometry_changed, geometry, use_virial, &
                        ri_energy, ri_operator, atom_gradient, strain_gradient)
      IF (geometry) THEN
         CALL qs_upf_inverse_aux_response(qs_env, xc_section, exchange, ri_energy, ri_operator, &
                                          atom_gradient, strain_gradient, atom_response, cell_response, &
                                          energy, operators, response, components)
      ELSE
         components(1) = ri_energy
         CALL qs_upf_inverse_aux_operator_evaluate(qs_env, xc_section, ri_operator, &
                                                   components(3), components(4), operators, &
                                                   exchange=exchange, core_energy=components(2))
      END IF
   END SUBROUTINE qs_upf_inverse_aux_evaluate

! **************************************************************************************************
!> \brief Sum native RI repetitions into owned original-coordinate outputs.
!>        Clear the temporary KS matrices separately for each repetition. Gamma obtains energy
!>        from those matrices, while k-point RI also owns a physical Bloch operator accumulator.
!>        Preserve the caller's auxiliary accumulator, its native KS matrices and primary storage.
!> \param qs_env ...
!> \param exchange ...
!> \param geometry_changed ...
!> \param geometry ...
!> \param use_virial ...
!> \param energy Globally reduced valence-valence RI HF energy.
!> \param operators Unweighted physical-k/spin auxiliary AO derivative.
!> \param atom_gradient Positive globally reduced Cartesian energy gradients, zero without geometry.
!> \param strain_gradient Positive globally reduced affine dE/dstrain, zero without requested strain.
! **************************************************************************************************
   SUBROUTINE auxiliary_ri(qs_env, exchange, geometry_changed, geometry, use_virial, &
                           energy, operators, atom_gradient, strain_gradient)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(hfx_type), INTENT(IN)                         :: exchange(:, :)
      LOGICAL, INTENT(IN)                                :: geometry_changed, geometry, use_virial
      REAL(dp), INTENT(OUT)                              :: energy
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: operators(:, :, :, :)
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: atom_gradient(:, :)
      REAL(dp), INTENT(OUT)                              :: strain_gradient(3, 3)

      COMPLEX(dp), ALLOCATABLE                           :: gamma(:, :), saved_exchange(:, :, :, :)
      INTEGER                                            :: contribution, ik, image, natom, ns, spin
      LOGICAL                                            :: bloch, valid
      REAL(dp)                                           :: fraction, term_energy, term_strain(3, 3)
      REAL(dp), ALLOCATABLE                              :: term_gradient(:, :)
      TYPE(admm_type), POINTER                           :: admm
      TYPE(dbcsr_p_type), POINTER                        :: density(:, :), native_ks(:, :), &
                                                            term_ks(:, :)
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_upf_inverse_scf_type), POINTER             :: state

      CALL get_qs_env(qs_env, admm_env=admm, natom=natom)
      CPASSERT(ASSOCIATED(admm))
      CPASSERT(admm%upf_density%valid)
      state => admm%upf_inverse
      CPASSERT(ASSOCIATED(state))
      CPASSERT(state%basis_type == 'AUX_FIT')
      CPASSERT(SIZE(exchange, 2) > 0)
      ns = SIZE(admm%upf_density%matrices, 4)
      bloch = ASSOCIATED(state%bloch%kpoints)
      IF (bloch) THEN
         DO ik = 1, state%bloch%kpoints%nkp
            IF (state%bloch%kpoints%kp_sym(ik)%kpoint_sym%apply_symmetry) THEN
               CPABORT('Auxiliary RI spatial symmetry is not implemented yet')
            END IF
         END DO
      END IF
      DO contribution = 1, SIZE(exchange, 1)
         fraction = exchange(contribution, 1)%general_parameter%fraction
         CPASSERT(ieee_is_finite(fraction))
         IF (fraction < 0.0_dp) THEN
            CPABORT('Auxiliary RI needs nonnegative HF fractions')
         END IF
         IF (.NOT. MERGE(0.5_dp*fraction, fraction, ns == 1) > EPSILON(0.0_dp)) CYCLE
         IF (.NOT. exchange(contribution, 1)%do_hfx_ri .OR. exchange(contribution, 1)%use_ace) THEN
            CPABORT('Original-UPF auxiliary exchange requires uncompressed RI')
         END IF
         CPASSERT(ASSOCIATED(exchange(contribution, 1)%ri_data))
         ASSOCIATE (ri => exchange(contribution, 1)%ri_data)
            IF (ri%input_flavor /= ri_pmat .OR. ri%orb_basis_type /= 'AUX_FIT') THEN
               CPABORT('Original-UPF auxiliary exchange requires RI-RHO in AUX_FIT')
            END IF
            valid = compare_potential_types(ri%hfx_pot, exchange(contribution, 1)%potential_parameter)
            CPASSERT(valid)
         END ASSOCIATE
      END DO
      CALL get_admm_env(admm, matrix_ks_aux_fit_kp=native_ks, rho_aux_fit=rho)
      CALL qs_rho_get(rho, rho_ao_kp=density)
      CPASSERT(ASSOCIATED(native_ks) .AND. ASSOCIATED(density))
      valid = ALL(SHAPE(native_ks) == SHAPE(density)) .AND. SIZE(native_ks, 1) == ns
      CPASSERT(valid)
      ALLOCATE (operators, MOLD=admm%upf_density%matrices)
      operators = CMPLX(0.0_dp, 0.0_dp, dp)
      ALLOCATE (atom_gradient(3, natom), term_gradient(3, natom))
      atom_gradient = 0.0_dp
      strain_gradient = 0.0_dp
      energy = 0.0_dp
      CALL dbcsr_allocate_matrix_set(term_ks, ns, SIZE(native_ks, 2))
      DO image = 1, SIZE(native_ks, 2)
         DO spin = 1, ns
            ALLOCATE (term_ks(spin, image)%matrix)
            CALL dbcsr_copy(term_ks(spin, image)%matrix, native_ks(spin, image)%matrix)
         END DO
      END DO
      IF (bloch) CALL MOVE_ALLOC(state%exchange_k, saved_exchange)
      DO contribution = 1, SIZE(exchange, 1)
         fraction = exchange(contribution, 1)%general_parameter%fraction
         IF (.NOT. MERGE(0.5_dp*fraction, fraction, ns == 1) > EPSILON(0.0_dp)) CYCLE
         DO image = 1, SIZE(term_ks, 2)
            DO spin = 1, ns
               CALL dbcsr_set(term_ks(spin, image)%matrix, 0.0_dp)
            END DO
         END DO
         IF (bloch) THEN
            IF (ALLOCATED(state%exchange_k)) state%exchange_k = CMPLX(0.0_dp, 0.0_dp, dp)
            CALL hfx_ri_update_ks_kp(qs_env, exchange(contribution, 1)%ri_data, term_ks, term_energy, &
                                     density, geometry_changed, ns, fraction)
            CPASSERT(ALLOCATED(state%exchange_k))
            valid = ALL(SHAPE(state%exchange_k) == SHAPE(operators))
            CPASSERT(valid)
            operators(:, :, :, :) = operators + state%exchange_k
            IF (geometry) THEN
               CALL hfx_ri_update_forces_kp(qs_env, exchange(contribution, 1)%ri_data, ns, fraction, density, &
                                            use_virial=use_virial, atom_gradient=term_gradient, strain_gradient=term_strain)
            END IF
         ELSE
            CALL hfx_ri_update_ks(qs_env, exchange(contribution, 1)%ri_data, term_ks, term_energy, &
                                  rho_ao=density, geometry_did_change=geometry_changed, nspins=ns, hf_fraction=fraction)
            DO spin = 1, ns
               CALL qs_upf_inverse_matrix_gather(state%frame%ao, term_ks(spin, 1)%matrix, gamma)
               operators(:, :, 1, spin) = operators(:, :, 1, spin) + gamma
            END DO
            IF (geometry) THEN
               CALL hfx_ri_update_forces(qs_env, exchange(contribution, 1)%ri_data, ns, fraction, rho_ao=density, &
                                         use_virial=use_virial, atom_gradient=term_gradient, strain_gradient=term_strain)
            END IF
         END IF
         energy = energy + term_energy
         IF (geometry) THEN
            atom_gradient(:, :) = atom_gradient + term_gradient
            strain_gradient = strain_gradient + term_strain
         END IF
      END DO
      IF (bloch) THEN
         IF (ALLOCATED(state%exchange_k)) DEALLOCATE (state%exchange_k)
         CALL MOVE_ALLOC(saved_exchange, state%exchange_k)
      END IF
      CALL dbcsr_deallocate_matrix_set(term_ks)
   END SUBROUTINE auxiliary_ri

! **************************************************************************************************
!> \brief Gather the actual native projection in original coordinates for every physical k.
!>        The two bases keep their own folding phases and symmetry actions. No projection is
!>        refactorized here. Native FM gathering is collective within its owning k-point group.
!> \param qs_env Prepared native primary and auxiliary projection owners.
!> \param projection New A(k), with auxiliary rows, orbital columns and no mesh weights.
!> \param reduced_projection Optional actual projection before physical star expansion.
!> \param metric_inverse Optional inverse of its retained native metric in representative order.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_aux_projection(qs_env, projection, reduced_projection, metric_inverse)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: projection(:, :, :)
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT), OPTIONAL    :: reduced_projection(:, :, :), &
                                                            metric_inverse(:, :, :)

      COMPLEX(dp), ALLOCATABLE                           :: inverse(:, :), native(:, :, :), &
                                                            native_metric(:, :, :), &
                                                            original(:, :, :)
      INTEGER                                            :: ik, local, na, no
      REAL(dp)                                           :: inverse_error
      REAL(dp), ALLOCATABLE                              :: real_metric(:, :), real_part(:, :)
      TYPE(admm_type), POINTER                           :: admm
      TYPE(kpoint_env_type), POINTER                     :: kp
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_upf_inverse_scf_type), POINTER             :: auxiliary, primary

      CALL get_qs_env(qs_env, ks_env=ks_env, admm_env=admm)
      CPASSERT(ASSOCIATED(admm))
      CALL get_ks_env(ks_env, upf_inverse=primary)
      auxiliary => admm%upf_inverse
      CPASSERT(ASSOCIATED(primary) .AND. ASSOCIATED(auxiliary))
      CPASSERT(primary%enabled)
      CPASSERT(admm%upf_aux_metric%ready .AND. admm%upf_mixed_metric%ready)
      na = admm%nao_aux_fit
      no = admm%nao_orb
      ALLOCATE (real_part(na, no))
      IF (ASSOCIATED(auxiliary%bloch%kpoints)) THEN
         kpoints => auxiliary%bloch%kpoints
         CPASSERT(ASSOCIATED(primary%bloch%kpoints, kpoints))
         ALLOCATE (native(na, no, kpoints%nkp))
         native = CMPLX(0.0_dp, 0.0_dp, dp)
         IF (PRESENT(metric_inverse)) THEN
            ALLOCATE (native_metric(na, na, kpoints%nkp), inverse(na, na), real_metric(na, na))
            native_metric = CMPLX(0.0_dp, 0.0_dp, dp)
         END IF
         DO local = 1, SIZE(kpoints%kp_aux_env)
            kp => kpoints%kp_aux_env(local)%kpoint_env
            ik = kp%nkpoint
            IF (kpoints%use_real_wfn) THEN
               CALL cp_fm_get_submatrix(kp%amat, real_part)
               IF (kpoints%para_env_kp%is_source()) native(:, :, ik) = CMPLX(real_part, 0.0_dp, dp)
            ELSE
               CALL cp_cfm_get_submatrix(kp%camat, native(:, :, ik))
               IF (.NOT. kpoints%para_env_kp%is_source()) native(:, :, ik) = CMPLX(0.0_dp, 0.0_dp, dp)
            END IF
            IF (PRESENT(metric_inverse)) THEN
               IF (kpoints%use_real_wfn) THEN
                  CPASSERT(ASSOCIATED(kp%smat%matrix_struct))
                  CALL cp_fm_get_submatrix(kp%smat, real_metric)
                  IF (kpoints%para_env_kp%is_source()) native_metric(:, :, ik) = CMPLX(real_metric, 0.0_dp, dp)
               ELSE
                  CPASSERT(ASSOCIATED(kp%csmat%matrix_struct))
                  CALL cp_cfm_get_submatrix(kp%csmat, native_metric(:, :, ik))
                  IF (.NOT. kpoints%para_env_kp%is_source()) native_metric(:, :, ik) = CMPLX(0.0_dp, 0.0_dp, dp)
               END IF
            END IF
         END DO
         DO ik = 1, kpoints%nkp
            CALL kpoints%para_env%sum(native(:, :, ik))
         END DO
         CALL qs_upf_inverse_kpoints_original_operator(auxiliary%bloch, native, original, right_frame=primary%bloch)
         CALL kpoint_star_density(auxiliary%bloch%star, original, projection, right_map=primary%bloch%star)
         IF (PRESENT(reduced_projection)) ALLOCATE (reduced_projection, SOURCE=original)
         IF (PRESENT(metric_inverse)) THEN
            DO ik = 1, kpoints%nkp
               CALL kpoints%para_env%sum(native_metric(:, :, ik))
               CALL invert_matrix(native_metric(:, :, ik), inverse, inverse_error)
               CPASSERT(inverse_error >= 0.0_dp)
               native_metric(:, :, ik) = inverse
            END DO
            CALL qs_upf_inverse_kpoints_original_operator(auxiliary%bloch, native_metric, metric_inverse)
         END IF
      ELSE
         CPASSERT(.NOT. ASSOCIATED(primary%bloch%kpoints))
         CPASSERT(ASSOCIATED(admm%A))
         CALL cp_fm_get_submatrix(admm%A, real_part)
         ALLOCATE (projection(na, no, 1))
         projection(:, :, 1) = CMPLX(real_part, 0.0_dp, dp)
         IF (PRESENT(reduced_projection)) ALLOCATE (reduced_projection, SOURCE=projection)
         IF (PRESENT(metric_inverse)) THEN
            ALLOCATE (real_metric(na, na), metric_inverse(na, na, 1))
            CALL cp_fm_get_submatrix(admm%S_inv, real_metric)
            metric_inverse(:, :, 1) = CMPLX(real_metric, 0.0_dp, dp)
         END IF
      END IF
   END SUBROUTINE qs_upf_inverse_aux_projection

! **************************************************************************************************
!> \brief Compose original auxiliary Fock/XC and transfer their unpurified electronic derivative.
!>        Use the actual native projection A and its complete original charge metric. For charge
!>        normalization, Dbar_A=eta A D_O A^dagger and eta=N_O/N_A. The global spin response
!>        lambda=Sum_k w_k Tr(K_A D_A)/N_A enters both normalization metrics. Geometry outputs
!>        additionally differentiate A=G_AA^-1 G_AO. Their scalar pairing applies mesh weights once.
!>        Reduce derivatives before differentiating the native projection. Its metric is the
!>        native star average, not a separately factorized metric for every physical partner.
!>        This component does not activate full ADMM SCF or supply explicit Fock/XC geometry terms.
!> \param qs_env Same native projection and auxiliary density used by the Fock producer.
!> \param xc_section Selected auxiliary XC functional, with its actual exchange scaling.
!> \param fock Original-coordinate auxiliary Fock operators at the current auxiliary density,
!>        including HF fractions but no k weights or additional explicit charge rescaling.
!> \param grid_energy Auxiliary uniform-grid XC energy.
!> \param atom_energy Auxiliary original atomic residual and selected PAW XC energy.
!> \param operators New primary operators, in original physical k and spin order, without weights.
!> \param mixed_covector Optional original-coordinate covectors of G_AO, required with both below.
!> \param auxiliary_covector Combined projection and normalization covectors of G_AA.
!> \param primary_covector Normalization covectors of the actual primary solver metric G_OO.
!> \param exchange Optional actual native HF interactions. The supplied Fock operator excludes PAW core exchange.
!> \param core_energy Optional PAW core-core and core-valence HF energy, required with exchange.
!>        These terms use the auxiliary density and contribute once to its transferred operator.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_aux_operator_evaluate(qs_env, xc_section, fock, grid_energy, atom_energy, operators, &
                                                   mixed_covector, auxiliary_covector, primary_covector, exchange, core_energy)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: xc_section
      COMPLEX(dp), INTENT(IN)                            :: fock(:, :, :, :)
      REAL(dp), INTENT(OUT)                              :: grid_energy, atom_energy
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: operators(:, :, :, :)
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT), OPTIONAL    :: mixed_covector(:, :, :, :), &
                                                            auxiliary_covector(:, :, :, :), &
                                                            primary_covector(:, :, :, :)
      TYPE(hfx_type), INTENT(IN), OPTIONAL               :: exchange(:, :)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: core_energy

      COMPLEX(dp), ALLOCATABLE                           :: auxiliary_operator(:, :, :, :), &
                                                            core_operator(:, :, :, :)
      LOGICAL                                            :: valid

      CPASSERT(PRESENT(exchange) .EQV. PRESENT(core_energy))
      CALL qs_upf_inverse_aux_xc_evaluate(qs_env, xc_section, grid_energy, atom_energy, auxiliary_operator)
      valid = ALL(SHAPE(fock) == SHAPE(auxiliary_operator))
      CPASSERT(valid)
      IF (PRESENT(exchange)) THEN
         CALL qs_upf_inverse_aux_core_exchange(qs_env, exchange, core_energy, core_operator)
         auxiliary_operator(:, :, :, :) = auxiliary_operator + core_operator
      END IF
      auxiliary_operator(:, :, :, :) = auxiliary_operator + fock
      CALL transfer_operator(qs_env, auxiliary_operator, operators, mixed_covector, auxiliary_covector, primary_covector)
   END SUBROUTINE qs_upf_inverse_aux_operator_evaluate

! **************************************************************************************************
!> \brief Transfer the derivative of the same auxiliary energy through the actual A and eta.
!> \param qs_env ...
!> \param auxiliary_operator Complete unweighted original-coordinate derivative at fixed auxiliary density.
!> \param operators ...
!> \param mixed_covector ...
!> \param auxiliary_covector ...
!> \param primary_covector ...
! **************************************************************************************************
   SUBROUTINE transfer_operator(qs_env, auxiliary_operator, operators, &
                                mixed_covector, auxiliary_covector, primary_covector)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      COMPLEX(dp), INTENT(IN)                            :: auxiliary_operator(:, :, :, :)
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: operators(:, :, :, :)
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT), OPTIONAL    :: mixed_covector(:, :, :, :), &
                                                            auxiliary_covector(:, :, :, :), &
                                                            primary_covector(:, :, :, :)

      COMPLEX(dp), ALLOCATABLE :: auxiliary_reduced(:, :, :), corrected(:, :), &
         corrected_k(:, :, :), density(:, :, :, :), expanded(:, :, :), gamma(:, :), &
         matrix(:, :, :), metric_inverse(:, :, :), mixed_reduced(:, :, :), &
         primary_metric(:, :, :), projection(:, :, :), reduced_corrected(:, :, :), &
         reduced_density(:, :, :), reduced_projection(:, :, :), unscaled(:, :), work(:, :), z(:, :)
      INTEGER                                            :: ik, na, nk, no, ns, spin
      LOGICAL                                            :: geometry, normalized, valid
      REAL(dp)                                           :: eta, lambda, na_spin, no_spin, &
                                                            numerator, tolerance
      TYPE(admm_type), POINTER                           :: admm
      TYPE(dbcsr_p_type), POINTER                        :: density_ao(:), metric(:)
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_upf_inverse_scf_type), POINTER             :: primary

      CALL get_qs_env(qs_env, ks_env=ks_env, admm_env=admm, rho=rho)
      CPASSERT(ASSOCIATED(admm))
      CPASSERT(admm%upf_density%valid)
      IF (admm%purification_method /= do_admm_purify_none .OR. admm%block_fit .OR. admm%block_dm) THEN
         CPABORT('Purified or blocked ADMM not implemented yet')
      END IF
      IF (admm%do_admmp .OR. admm%do_admms .OR. admm%scaling_model /= do_admm_exch_scaling_none) THEN
         CPABORT('ADMM exchange rescaling not implemented yet')
      END IF
      geometry = PRESENT(mixed_covector) .OR. PRESENT(auxiliary_covector) .OR. PRESENT(primary_covector)
      IF (geometry) THEN
         valid = PRESENT(mixed_covector) .AND. PRESENT(auxiliary_covector) .AND. PRESENT(primary_covector)
         CPASSERT(valid)
      END IF
      normalized = admm%do_admmq
      IF (normalized) THEN
         CPASSERT(admm%charge_constrain)
      END IF
      CALL get_ks_env(ks_env, upf_inverse=primary)
      CPASSERT(ASSOCIATED(primary))
      CPASSERT(primary%enabled)
      IF (geometry) THEN
         CALL qs_upf_inverse_aux_projection(qs_env, projection, reduced_projection, metric_inverse)
      ELSE
         CALL qs_upf_inverse_aux_projection(qs_env, projection)
      END IF
      na = SIZE(projection, 1)
      no = SIZE(projection, 2)
      nk = SIZE(projection, 3)
      ns = SIZE(admm%upf_density%matrices, 4)
      valid = ALL(SHAPE(auxiliary_operator) == [na, na, nk, ns])
      CPASSERT(valid)
      CPASSERT(SIZE(admm%upf_aux_metric%weights) == nk)
      ALLOCATE (density(no, no, nk, ns), primary_metric(no, no, nk))
      IF (ASSOCIATED(primary%bloch%kpoints)) THEN
         valid = ALL(SHAPE(primary%metric_k) == [no, no, nk])
         CPASSERT(valid)
         primary_metric(:, :, :) = primary%metric_k
         DO spin = 1, ns
            CALL qs_upf_inverse_kpoints_density(primary%bloch, spin, matrix)
            density(:, :, :, spin) = matrix
         END DO
      ELSE
         CPASSERT(nk == 1)
         metric => get_1d_pointer(primary%metric)
         CALL qs_upf_inverse_matrix_gather(primary%frame%ao, metric(1)%matrix, gamma)
         primary_metric(:, :, 1) = gamma
         CALL qs_rho_get(rho, rho_ao=density_ao)
         DO spin = 1, ns
            CALL qs_upf_inverse_matrix_gather(primary%frame%ao, density_ao(spin)%matrix, gamma)
            density(:, :, 1, spin) = gamma
         END DO
      END IF
      ALLOCATE (operators(no, no, nk, ns), unscaled(na, na), corrected(na, na), work(na, no), z(na, no))
      IF (geometry) THEN
         ALLOCATE (mixed_covector(na, no, nk, ns), auxiliary_covector(na, na, nk, ns), &
                   primary_covector(no, no, nk, ns), corrected_k(na, na, nk), &
                   mixed_reduced(na, no, SIZE(reduced_projection, 3)), &
                   auxiliary_reduced(na, na, SIZE(reduced_projection, 3)))
      END IF
      DO spin = 1, ns
         eta = 1.0_dp
         lambda = 0.0_dp
         IF (normalized) THEN
            eta = admm%gsi(spin)
            CPASSERT(eta > 0.0_dp)
            na_spin = 0.0_dp
            no_spin = 0.0_dp
            numerator = 0.0_dp
            DO ik = 1, nk
               unscaled(:, :) = admm%upf_density%matrices(:, :, ik, spin)/eta
               na_spin = na_spin + admm%upf_aux_metric%weights(ik)* &
                         REAL(SUM(CONJG(admm%upf_aux_metric%values(:, :, ik))*unscaled), dp)
               no_spin = no_spin + admm%upf_aux_metric%weights(ik)* &
                         REAL(SUM(CONJG(primary_metric(:, :, ik))*density(:, :, ik, spin)), dp)
               numerator = numerator + admm%upf_aux_metric%weights(ik)* &
                           REAL(SUM(CONJG(auxiliary_operator(:, :, ik, spin))*unscaled), dp)
            END DO
            tolerance = 128.0_dp*EPSILON(1.0_dp)*MAX(1.0_dp, ABS(na_spin), ABS(no_spin))
            IF (ABS(na_spin) > tolerance .OR. ABS(no_spin) > tolerance) THEN
               CPASSERT(na_spin > 0.0_dp)
               lambda = numerator/na_spin
            END IF
         END IF
         DO ik = 1, nk
            unscaled(:, :) = admm%upf_density%matrices(:, :, ik, spin)/eta
            corrected(:, :) = eta*(auxiliary_operator(:, :, ik, spin) - lambda*admm%upf_aux_metric%values(:, :, ik))
            work(:, :) = MATMUL(corrected, projection(:, :, ik))
            operators(:, :, ik, spin) = MATMUL(CONJG(TRANSPOSE(projection(:, :, ik))), work) + &
                                        lambda*primary_metric(:, :, ik)
            IF (geometry) THEN
               corrected_k(:, :, ik) = corrected
               auxiliary_covector(:, :, ik, spin) = -eta*lambda*unscaled
               primary_covector(:, :, ik, spin) = lambda*density(:, :, ik, spin)
            END IF
         END DO
         IF (geometry) THEN
            IF (ASSOCIATED(primary%bloch%kpoints)) THEN
               CALL kpoint_star_operator(admm%upf_inverse%bloch%star, corrected_k, reduced_corrected)
               CALL kpoint_star_operator(primary%bloch%star, density(:, :, :, spin), reduced_density)
            ELSE
               IF (ALLOCATED(reduced_corrected)) DEALLOCATE (reduced_corrected, reduced_density)
               ALLOCATE (reduced_corrected, SOURCE=corrected_k)
               ALLOCATE (reduced_density, SOURCE=density(:, :, :, spin))
            END IF
            DO ik = 1, SIZE(reduced_projection, 3)
               work(:, :) = MATMUL(reduced_corrected(:, :, ik), reduced_projection(:, :, ik))
               z(:, :) = MATMUL(metric_inverse(:, :, ik), MATMUL(work, reduced_density(:, :, ik)))
               mixed_reduced(:, :, ik) = 2.0_dp*z
               auxiliary_reduced(:, :, ik) = -MATMUL(z, CONJG(TRANSPOSE(reduced_projection(:, :, ik)))) - &
                                             MATMUL(reduced_projection(:, :, ik), CONJG(TRANSPOSE(z)))
            END DO
            IF (ASSOCIATED(primary%bloch%kpoints)) THEN
               CALL kpoint_star_density(admm%upf_inverse%bloch%star, mixed_reduced, expanded, &
                                        right_map=primary%bloch%star)
               mixed_covector(:, :, :, spin) = expanded
               CALL kpoint_star_density(admm%upf_inverse%bloch%star, auxiliary_reduced, expanded)
               auxiliary_covector(:, :, :, spin) = auxiliary_covector(:, :, :, spin) + expanded
            ELSE
               mixed_covector(:, :, :, spin) = mixed_reduced
               auxiliary_covector(:, :, :, spin) = auxiliary_covector(:, :, :, spin) + auxiliary_reduced
            END IF
         END IF
      END DO
   END SUBROUTINE transfer_operator

! **************************************************************************************************
!> \brief Original PAW closed-core HF energy and its auxiliary AO and geometry derivatives.
!>        Core-core exchange is constant. Core-valence exchange is linear in the original
!>        projector occupations, without a further factor one half. Sum native HF repetitions
!>        once and ignore thread replicas. Valence-valence exchange remains in the RI producer.
!>        No PAW Hartree, semilocal XC or fixed dataset reference energy is added here.
!> \param qs_env Valid current auxiliary occupations and their independent original-UPF representation.
!> \param exchange The same native HF interactions, fractions and physical range parameters as RI.
!> \param energy Selected frozen-core HF energy, including each actual PAW atom once.
!> \param operators Original physical-k/spin auxiliary AO derivative, without mesh weights.
!> \param atom_response Optional total atomic motion, including affine cell motion.
!> \param cell_response Optional physical lattice derivative at fixed reduced k coordinates.
!> \param response Optional fixed-auxiliary-occupation directional derivative. Supply all three geometry arguments.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_aux_core_exchange(qs_env, exchange, energy, operators, &
                                               atom_response, cell_response, response)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(hfx_type), INTENT(IN)                         :: exchange(:, :)
      REAL(dp), INTENT(OUT)                              :: energy
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: operators(:, :, :, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: atom_response(:, :), cell_response(3, 3)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: response

      INTEGER                                            :: kind, nk, ns
      LOGICAL                                            :: geometry
      REAL(dp)                                           :: energies(3)
      REAL(dp), ALLOCATABLE                              :: core_energy(:), reduced_k_response(:, :)
      TYPE(admm_type), POINTER                           :: admm
      TYPE(paw_potential_type), ALLOCATABLE              :: core_potential(:), potential(:)
      TYPE(qs_kind_type), POINTER                        :: kinds(:)
      TYPE(qs_upf_inverse_scf_type), POINTER             :: state
      TYPE(section_vals_type), POINTER                   :: no_functional

      CALL get_qs_env(qs_env, admm_env=admm, qs_kind_set=kinds)
      CPASSERT(ASSOCIATED(admm))
      CPASSERT(admm%upf_density%valid)
      state => admm%upf_inverse
      CPASSERT(ASSOCIATED(state))
      CPASSERT(state%basis_type == 'AUX_FIT')
      CPASSERT(SIZE(exchange, 2) > 0)
      geometry = PRESENT(response)
      CPASSERT(geometry .EQV. PRESENT(atom_response))
      CPASSERT(geometry .EQV. PRESENT(cell_response))
      energy = 0.0_dp
      ALLOCATE (operators, MOLD=admm%upf_density%matrices)
      operators = CMPLX(0.0_dp, 0.0_dp, dp)
      IF (geometry) response = 0.0_dp
      IF (.NOT. ANY(state%paw%active)) RETURN
      ns = SIZE(admm%upf_density%matrices, 4)
      nk = SIZE(admm%upf_density%matrices, 3)
      ALLOCATE (core_energy(SIZE(kinds)), SOURCE=0.0_dp)
      ALLOCATE (core_potential(SIZE(kinds)))
      DO kind = 1, SIZE(kinds)
         IF (.NOT. state%paw(kind)%active) CYCLE
         ASSOCIATE (pot => kinds(kind)%sgp_potential%upf_data, projectors => kinds(kind)%sgp_potential%upf_projectors)
            ALLOCATE (core_potential(kind)%matrix(projectors%nproj, projectors%nproj, 1))
            CALL qs_upf_inverse_paw_core_prepare(pot, projectors, state%paw(kind), exchange, ns, &
                                                 core_energy(kind), core_potential(kind)%matrix(:, :, 1))
         END ASSOCIATE
      END DO
      NULLIFY (no_functional)
      IF (ASSOCIATED(state%bloch%kpoints)) THEN
         CALL qs_upf_inverse_paw_evaluate(state%paw, state%bloch%frames, admm%upf_density%matrices, state%bloch%weights, &
                                          no_functional, 0.0_dp, 0.0_dp, energies, operators, potential, &
                                      semilocal_xc=.FALSE., core_energy=core_energy, core_potential=core_potential, hartree=.FALSE.)
      ELSE
         CPASSERT(nk == 1)
         CALL qs_upf_inverse_paw_evaluate(state%paw, [state%frame], admm%upf_density%matrices, [1.0_dp], &
                                          no_functional, 0.0_dp, 0.0_dp, energies, operators, potential, &
                                      semilocal_xc=.FALSE., core_energy=core_energy, core_potential=core_potential, hartree=.FALSE.)
      END IF
      energy = energies(2)
      IF (geometry) THEN
         ALLOCATE (reduced_k_response(3, nk), SOURCE=0.0_dp)
         IF (ASSOCIATED(state%bloch%kpoints)) THEN
            CALL qs_upf_inverse_paw_response(state%bloch%frames, admm%upf_density%matrices, state%bloch%weights, potential, &
                                             atom_response, cell_response, reduced_k_response, response)
         ELSE
            CALL qs_upf_inverse_paw_response([state%frame], admm%upf_density%matrices, [1.0_dp], potential, &
                                             atom_response, cell_response, reduced_k_response, response)
         END IF
      END IF
   END SUBROUTINE qs_upf_inverse_aux_core_exchange

! **************************************************************************************************
!> \brief Compose auxiliary RI, original PAW core, XC and projection/normalization geometry.
!>        Native RI supplies energy, its auxiliary AO derivative and globally reduced Cartesian
!>        and affine-strain derivatives at the current fixed auxiliary occupation matrices.
!>        All RI quantities include their HF fractions. Energy and geometry include physical
!>        mesh weights once, while the supplied AO derivative has no mesh weights.
!>        The energy and primary AO operator use exactly the same auxiliary functional. Its
!>        explicit geometry is combined with the implicit response of the native A and eta.
!>        This component does not activate the complete native ADMM SCF calculation.
!> \param qs_env The actual native primary and auxiliary density/projection owners.
!> \param xc_section The selected auxiliary XC correction with its actual exchange scaling.
!> \param exchange The same native HF interactions as the supplied RI contribution.
!> \param ri_energy Native valence-valence RI exchange energy at the current auxiliary density.
!> \param ri_operator Its unweighted original-coordinate auxiliary AO derivative (AO,AO,k,spin).
!> \param ri_atom_gradient Globally reduced positive Cartesian RI energy gradients (3,atom).
!> \param ri_strain_gradient Globally reduced positive dE_RI/dstrain, including affine atomic motion.
!> \param atom_response Total physical atomic motion, including affine cell motion.
!> \param cell_response Physical lattice derivative, keeping reduced k coordinates fixed.
!> \param energy Complete auxiliary RI plus core-HF and XC correction energy.
!> \param operators Complete primary AO derivative through the actual projection and charge constraint.
!> \param response Complete directional derivative at fixed original-coordinate primary occupations.
!> \param components Optional RI valence HF, PAW core HF, grid XC and atomic XC energies.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_aux_response(qs_env, xc_section, exchange, ri_energy, ri_operator, &
                                          ri_atom_gradient, ri_strain_gradient, atom_response, cell_response, &
                                          energy, operators, response, components)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: xc_section
      TYPE(hfx_type), INTENT(IN)                         :: exchange(:, :)
      REAL(dp), INTENT(IN)                               :: ri_energy
      COMPLEX(dp), INTENT(IN)                            :: ri_operator(:, :, :, :)
      REAL(dp), INTENT(IN)                               :: ri_atom_gradient(:, :), &
                                                            ri_strain_gradient(3, 3), &
                                                            atom_response(:, :), &
                                                            cell_response(3, 3)
      REAL(dp), INTENT(OUT)                              :: energy
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: operators(:, :, :, :)
      REAL(dp), INTENT(OUT)                              :: response
      REAL(dp), INTENT(OUT), OPTIONAL                    :: components(4)

      COMPLEX(dp), ALLOCATABLE :: auxiliary_covector(:, :, :, :), auxiliary_operator(:, :, :, :), &
         core_operator(:, :, :, :), mixed_covector(:, :, :, :), primary_covector(:, :, :, :)
      INTEGER                                            :: atom
      LOGICAL                                            :: valid
      REAL(dp) :: atom_energy, core_energy, core_response, grid_energy, projection_response, &
         ri_response, strain(3, 3), xc_response
      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_type), POINTER                       :: atoms(:)

      CALL get_qs_env(qs_env, cell=cell, particle_set=atoms)
      valid = ALL(SHAPE(atom_response) == [3, SIZE(atoms)]) .AND. &
              ALL(SHAPE(ri_atom_gradient) == SHAPE(atom_response))
      CPASSERT(valid)
      valid = ieee_is_finite(ri_energy) .AND. ALL(ieee_is_finite(ri_atom_gradient)) .AND. &
              ALL(ieee_is_finite(ri_strain_gradient))
      CPASSERT(valid)
      CALL qs_upf_inverse_aux_xc_evaluate(qs_env, xc_section, grid_energy, atom_energy, auxiliary_operator, &
                                          atom_response=atom_response, cell_response=cell_response, response=xc_response)
      valid = ALL(SHAPE(ri_operator) == SHAPE(auxiliary_operator))
      CPASSERT(valid)
      CALL qs_upf_inverse_aux_core_exchange(qs_env, exchange, core_energy, core_operator, &
                                            atom_response, cell_response, core_response)
      auxiliary_operator(:, :, :, :) = auxiliary_operator + core_operator
      auxiliary_operator(:, :, :, :) = auxiliary_operator + ri_operator
      CALL transfer_operator(qs_env, auxiliary_operator, operators, mixed_covector, auxiliary_covector, primary_covector)
      CALL qs_upf_inverse_aux_metric_response(qs_env, mixed_covector, auxiliary_covector, primary_covector, &
                                              atom_response, cell_response, projection_response)
      strain = MATMUL(cell_response, cell%h_inv)
      ri_response = SUM(ri_strain_gradient*strain)
      DO atom = 1, SIZE(atoms)
         ri_response = ri_response + SUM(ri_atom_gradient(:, atom)* &
                                         (atom_response(:, atom) - MATMUL(strain, atoms(atom)%r)))
      END DO
      energy = ri_energy + core_energy + grid_energy + atom_energy
      IF (PRESENT(components)) components = [ri_energy, core_energy, grid_energy, atom_energy]
      response = ri_response + core_response + xc_response + projection_response
   END SUBROUTINE qs_upf_inverse_aux_response

! **************************************************************************************************
!> \brief Contract projection and charge-normalization covectors with the actual moving metrics.
!>        Each native metric service applies its defining quadrature and k weights once. The
!>        primary normalization uses the installed solver metric, including the complete GAPW
!>        density metric where present. Frozen core charges never enter these valence metrics.
!>        This is the implicit projection/normalization response. Explicit auxiliary Fock and
!>        XC geometry derivatives at fixed auxiliary occupation matrices remain separate terms.
!> \param qs_env Same native primary and auxiliary owners used to obtain the covectors.
!> \param mixed_covector Original-coordinate G_AO covectors, without mesh weights.
!> \param auxiliary_covector Combined projection and normalization G_AA covectors.
!> \param primary_covector Normalization G_OO covectors in the actual solver basis.
!> \param atom_response Total atomic motion, including any affine cell motion.
!> \param cell_response Lattice-vector derivative at fixed reduced k coordinates.
!> \param response Real directional energy derivative, without a force or virial minus sign.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_aux_metric_response(qs_env, mixed_covector, auxiliary_covector, primary_covector, &
                                                 atom_response, cell_response, response)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      COMPLEX(dp), INTENT(IN)                            :: mixed_covector(:, :, :, :), &
                                                            auxiliary_covector(:, :, :, :), &
                                                            primary_covector(:, :, :, :)
      REAL(dp), INTENT(IN)                               :: atom_response(:, :), cell_response(3, 3)
      REAL(dp), INTENT(OUT)                              :: response

      COMPLEX(dp), ALLOCATABLE                           :: auxiliary(:, :, :), mixed(:, :, :), &
                                                            orbital(:, :, :)
      INTEGER                                            :: na, nk, no, ns
      LOGICAL                                            :: valid
      REAL(dp)                                           :: auxiliary_response, mixed_response, &
                                                            orbital_response
      TYPE(admm_type), POINTER                           :: admm
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_upf_inverse_scf_type), POINTER             :: primary

      CALL get_qs_env(qs_env, ks_env=ks_env, admm_env=admm)
      CPASSERT(ASSOCIATED(admm))
      CPASSERT(ASSOCIATED(admm%upf_inverse))
      CPASSERT(admm%upf_aux_metric%ready .AND. admm%upf_mixed_metric%ready)
      CALL get_ks_env(ks_env, upf_inverse=primary)
      CPASSERT(ASSOCIATED(primary))
      na = SIZE(admm%upf_aux_metric%values, 1)
      no = SIZE(admm%upf_mixed_metric%values, 2)
      nk = SIZE(admm%upf_aux_metric%weights)
      ns = SIZE(mixed_covector, 4)
      valid = ALL(SHAPE(mixed_covector) == [na, no, nk, ns]) .AND. &
              ALL(SHAPE(auxiliary_covector) == [na, na, nk, ns]) .AND. &
              ALL(SHAPE(primary_covector) == [no, no, nk, ns])
      CPASSERT(valid .AND. ns > 0)
      ALLOCATE (auxiliary(na, na, nk), mixed(na, no, nk), orbital(no, no, nk))
      auxiliary(:, :, :) = SUM(auxiliary_covector, DIM=4)
      mixed(:, :, :) = SUM(mixed_covector, DIM=4)
      orbital(:, :, :) = SUM(primary_covector, DIM=4)
      CALL qs_upf_inverse_charge_metric_contract(qs_env, admm%upf_inverse, primary, admm%upf_mixed_metric, &
                                                 mixed, atom_response, cell_response, mixed_response)
      CALL qs_upf_inverse_charge_metric_contract(qs_env, admm%upf_inverse, admm%upf_inverse, admm%upf_aux_metric, &
                                                 auxiliary, atom_response, cell_response, auxiliary_response)
      orbital_response = 0.0_dp
      IF (ANY(orbital /= CMPLX(0.0_dp, 0.0_dp, dp))) THEN
         IF (ASSOCIATED(primary%coulomb)) THEN
            CALL qs_upf_inverse_coulomb_metric_response(qs_env, primary, primary%coulomb, orbital, &
                                                        atom_response, cell_response, orbital_response)
         ELSE
            CALL primary_grid_metric_response(qs_env, primary, orbital, atom_response, cell_response, orbital_response)
         END IF
      END IF
      response = mixed_response + auxiliary_response + orbital_response
   END SUBROUTINE qs_upf_inverse_aux_metric_response

! **************************************************************************************************
!> \brief Contract the primary solver's original grid metric before a complete metric is installed.
!>        Match its native Gaussian reference and inverse/Q remainder. Grid and Gaussian source
!>        measures enter once. A later complete density metric uses its own service above.
!> \param qs_env Native geometry and primary Gaussian task owner.
!> \param state Primary inverse representation without an installed complete Coulomb metric.
!> \param matrix Hermitian valence metric covector, without mesh weights, summed over spins.
!> \param atom_response Total atomic motion.
!> \param cell_response Lattice response at fixed reduced k.
!> \param response Directional metric contraction, without a Pulay minus sign.
! **************************************************************************************************
   SUBROUTINE primary_grid_metric_response(qs_env, state, matrix, atom_response, cell_response, response)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_upf_inverse_scf_type), INTENT(INOUT), &
         TARGET                                          :: state
      COMPLEX(dp), INTENT(IN)                            :: matrix(:, :, :)
      REAL(dp), INTENT(IN)                               :: atom_response(:, :), cell_response(3, 3)
      REAL(dp), INTENT(OUT)                              :: response

      COMPLEX(dp), ALLOCATABLE                           :: derivative(:, :), metric(:, :)
      INTEGER                                            :: info, k, nk
      REAL(dp)                                           :: charge, inverse_cell(3, 3), &
                                                            strain(3, 3), trace_strain, &
                                                            virial(3, 3), weight
      REAL(dp), ALLOCATABLE                              :: forces(:, :), point_response(:, :)
      TYPE(pw_c1d_gs_type)                               :: density_g
      TYPE(pw_r3d_rs_type)                               :: density, one
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_upf_inverse_grid_type), POINTER            :: frame

      nk = 1
      frame => state%frame
      IF (ASSOCIATED(state%bloch%kpoints)) THEN
         nk = SIZE(state%bloch%weights)
         frame => state%bloch%frames(1)
      END IF
      CPASSERT(SIZE(matrix, 3) == nk)
      CPASSERT(ALL(SHAPE(atom_response) == SHAPE(frame%ao%positions)))
      inverse_cell = frame%ao%cell
      CALL invmat(inverse_cell, info)
      CPASSERT(info == 0)
      strain = MATMUL(cell_response, inverse_cell)
      trace_strain = strain(1, 1) + strain(2, 2) + strain(3, 3)
      response = 0.0_dp
      IF (state%native_gaussian) THEN
         CALL get_qs_env(qs_env, ks_env=ks_env)
         CALL one%create(state%source_grid)
         CALL density%create(state%source_grid)
         CALL density_g%create(state%source_grid)
         one%array = 1.0_dp
         CALL qs_upf_inverse_background_collocate(state%gaussian, ks_env, matrix, density, density_g, total_rho=charge)
         ALLOCATE (forces, MOLD=atom_response)
         CALL qs_upf_inverse_background_geometry(state%gaussian, qs_env, matrix, one, forces, virial)
         response = trace_strain*charge - SUM(virial*strain) - &
                    SUM(forces*(atom_response - MATMUL(strain, frame%ao%positions)))
         CALL one%release()
         CALL density%release()
         CALL density_g%release()
      END IF
      DO k = 1, nk
         weight = 1.0_dp
         IF (ASSOCIATED(state%bloch%kpoints)) THEN
            frame => state%bloch%frames(k)
            weight = state%bloch%weights(k)
         END IF
         ALLOCATE (point_response, SOURCE=MATMUL(strain, frame%ao%points))
         CALL qs_upf_inverse_grid_metric_dense(frame, frame%ao%radial_q, metric, &
                                               orbital_fields=.NOT. state%native_gaussian, response=derivative, &
                                               atom_response=atom_response, cell_response=cell_response, &
                                               reduced_k_response=[0.0_dp, 0.0_dp, 0.0_dp], point_response=point_response, &
                                               weights_response=trace_strain*frame%ao%weights)
         response = response + weight*REAL(SUM(CONJG(matrix(:, :, k))*derivative), dp)
         DEALLOCATE (point_response)
      END DO
   END SUBROUTINE primary_grid_metric_response
END MODULE qs_upf_inverse_aux_operator
