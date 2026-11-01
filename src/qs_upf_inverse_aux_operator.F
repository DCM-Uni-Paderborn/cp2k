!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
!> \brief Original-coordinate auxiliary operators and physical projection/charge adjoints.
MODULE qs_upf_inverse_aux_operator
   USE admm_types,                      ONLY: admm_type
   USE cp_cfm_types,                    ONLY: cp_cfm_get_submatrix
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type
   USE cp_fm_types,                     ONLY: cp_fm_get_submatrix
   USE input_constants,                 ONLY: do_admm_exch_scaling_none,&
                                              do_admm_purify_none
   USE input_section_types,             ONLY: section_vals_type
   USE kinds,                           ONLY: dp
   USE kpoint_methods,                  ONLY: kpoint_star_density,&
                                              kpoint_star_operator
   USE kpoint_transitional,             ONLY: get_1d_pointer
   USE kpoint_types,                    ONLY: kpoint_env_type,&
                                              kpoint_type
   USE mathlib,                         ONLY: invert_matrix
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_upf_inverse_aux_xc,           ONLY: qs_upf_inverse_aux_xc_evaluate
   USE qs_upf_inverse_grid,             ONLY: qs_upf_inverse_matrix_gather
   USE qs_upf_inverse_kpoints,          ONLY: qs_upf_inverse_kpoints_density,&
                                              qs_upf_inverse_kpoints_original_operator
   USE qs_upf_inverse_types,            ONLY: qs_upf_inverse_scf_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: qs_upf_inverse_aux_projection, qs_upf_inverse_aux_operator_evaluate

CONTAINS

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
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_aux_operator_evaluate(qs_env, xc_section, fock, grid_energy, atom_energy, operators, &
                                                   mixed_covector, auxiliary_covector, primary_covector)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: xc_section
      COMPLEX(dp), INTENT(IN)                            :: fock(:, :, :, :)
      REAL(dp), INTENT(OUT)                              :: grid_energy, atom_energy
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: operators(:, :, :, :)
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT), OPTIONAL    :: mixed_covector(:, :, :, :), &
                                                            auxiliary_covector(:, :, :, :), &
                                                            primary_covector(:, :, :, :)

      COMPLEX(dp), ALLOCATABLE :: auxiliary_operator(:, :, :, :), auxiliary_reduced(:, :, :), &
         corrected(:, :), corrected_k(:, :, :), density(:, :, :, :), expanded(:, :, :), &
         gamma(:, :), matrix(:, :, :), metric_inverse(:, :, :), mixed_reduced(:, :, :), &
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
      valid = ALL(SHAPE(fock) == [na, na, nk, ns])
      CPASSERT(valid)
      CPASSERT(SIZE(admm%upf_aux_metric%weights) == nk)
      CALL qs_upf_inverse_aux_xc_evaluate(qs_env, xc_section, grid_energy, atom_energy, auxiliary_operator)
      auxiliary_operator(:, :, :, :) = auxiliary_operator + fock
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
   END SUBROUTINE qs_upf_inverse_aux_operator_evaluate
END MODULE qs_upf_inverse_aux_operator
