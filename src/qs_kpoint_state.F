!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Methods for preparing and committing k-point orbital states to the QS environment.
! **************************************************************************************************
MODULE qs_kpoint_state
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_scale,&
                                              cp_cfm_trace,&
                                              cp_cfm_vectorsnorm
   USE cp_cfm_diag,                     ONLY: cp_cfm_geeig,&
                                              cp_cfm_heevd
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_release,&
                                              cp_cfm_set_all,&
                                              cp_cfm_to_cfm,&
                                              cp_cfm_to_fm,&
                                              cp_cfm_type,&
                                              cp_fm_to_cfm
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                              cp_fm_to_fm,&
                                              cp_fm_type,&
                                              cp_fm_vectorsnorm
   USE dft_plus_u,                      ONLY: lowdin_kp_apply
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE input_constants,                 ONLY: smear_fermi_dirac
   USE kinds,                           ONLY: dp
   USE kpoint_methods,                  ONLY: kpoint_density_matrices,&
                                              kpoint_density_transform,&
                                              kpoint_set_mo_occupation
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_env_type,&
                                              kpoint_type
   USE mathconstants,                   ONLY: z_one,&
                                              z_zero
   USE message_passing,                 ONLY: mp_para_env_type
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kpoint_density_state,         ONLY: kpoint_density_state_accept,&
                                              kpoint_density_state_active,&
                                              kpoint_density_state_capture
   USE qs_kpoint_operators,             ONLY: kpoint_slot_s
   USE qs_ks_types,                     ONLY: qs_ks_did_change
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_rho_methods,                  ONLY: qs_rho_update_rho
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_scf_types,                    ONLY: qs_scf_env_type
   USE scf_control_types,               ONLY: scf_control_type,&
                                              smear_type
   USE smearing_utils,                  ONLY: Smearkp,&
                                              Smearkp2,&
                                              smearing_separate_spin,&
                                              smearing_spin_targets
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_kpoint_state'

   PUBLIC :: qs_kpoint_copy_spin_mos, qs_kpoint_mos_initialized, qs_kpoint_set_fixed_occupations, &
             qs_kpoint_state_canonicalize_fixed, qs_kpoint_state_commit, &
             qs_kpoint_state_prepare_fixed_density, qs_kpoint_state_prepare_thermal_restart

CONTAINS

! **************************************************************************************************
!> \brief Reconstructs fixed-rank complex k-point orbitals from a density kernel.
!>        The overlap must first be transferred to kp%ot_smat. Use the complete native density
!>        when active, otherwise the density transferred to kp%ot_hmat. The
!>        largest natural-orbital subspace gives a metric-defined fixed-rank projector without
!>        replacing the guess by the ground-state projector of its initial Hamiltonian.
!> \param kpoints K-point environment containing the transferred density and overlap matrices.
! **************************************************************************************************
   SUBROUTINE qs_kpoint_state_prepare_fixed_density(kpoints)
      TYPE(kpoint_type), POINTER                         :: kpoints

      INTEGER                                            :: ikpoint_local, ispin, nactive, nmo, &
                                                            nocc
      LOGICAL                                            :: use_real_wfn
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: density_eigenvalues
      TYPE(cp_cfm_type)                                  :: matrix_p, matrix_s, projected
      TYPE(cp_cfm_type), POINTER                         :: cmo_coeff
      TYPE(cp_fm_struct_type), POINTER                   :: active_struct, square_struct
      TYPE(kpoint_env_type), POINTER                     :: kp

      NULLIFY (active_struct, cmo_coeff, kp, square_struct)

      CPASSERT(ASSOCIATED(kpoints))
      CPASSERT(ASSOCIATED(kpoints%kp_env))
      CALL get_kpoint_info(kpoints, use_real_wfn=use_real_wfn)
      IF (use_real_wfn) THEN
         CPABORT("Fixed-density k-point OT initialization requires complex wavefunctions.")
      END IF

      DO ikpoint_local = 1, SIZE(kpoints%kp_env)
         kp => kpoints%kp_env(ikpoint_local)%kpoint_env
         CPASSERT(ASSOCIATED(kp))
         CPASSERT(ASSOCIATED(kp%mos))
         CPASSERT(ASSOCIATED(kp%ot_hmat))
         CPASSERT(ASSOCIATED(kp%ot_smat))
         CPASSERT(SIZE(kp%ot_hmat, 1) >= 2)
         CPASSERT(SIZE(kp%ot_hmat, 2) >= SIZE(kp%mos))
         CPASSERT(SIZE(kp%ot_smat) >= 2)

         CALL cp_fm_get_info(kp%ot_hmat(1, 1), matrix_struct=square_struct)
         CALL cp_cfm_create(matrix_p, square_struct, name="natural-orbital density P(k)")
         CALL cp_cfm_create(matrix_s, square_struct, name="natural-orbital overlap S(k)")
         CALL cp_fm_to_cfm(kp%ot_smat(1), kp%ot_smat(2), matrix_s)

         DO ispin = 1, SIZE(kp%mos)
            CALL get_mo_set(kp%mos(ispin), homo=nocc, cmo_coeff=cmo_coeff, nmo=nmo)
            CPASSERT(nocc >= 0 .AND. nocc <= nmo)
            IF (nmo == 0) CYCLE
            CPASSERT(ASSOCIATED(cmo_coeff))
            ! Retained empty-spin virtuals need a metric-normalized reference, although P is zero.
            nactive = nocc
            IF (nocc == 0) nactive = nmo

            CALL cp_fm_struct_create(active_struct, template_fmstruct=cmo_coeff%matrix_struct, &
                                     ncol_global=nactive)
            IF (kpoint_density_state_active(kpoints)) THEN
               CPASSERT(kpoints%native_density_seeded)
               ! The real-space density can omit images of the complete Bloch kernel.
               CALL cp_cfm_to_cfm(kpoints%native_density(ikpoint_local, ispin), matrix_p)
            ELSE
               CALL cp_fm_to_cfm(kp%ot_hmat(1, ispin), kp%ot_hmat(2, ispin), matrix_p)
            END IF
            CALL cp_cfm_create(projected, active_struct, name="fixed-density natural orbitals")
            ALLOCATE (density_eigenvalues(nactive))
            CALL natural_orbitals(matrix_p, matrix_s, projected, density_eigenvalues)
            DEALLOCATE (density_eigenvalues)

            CALL cp_cfm_set_all(cmo_coeff, CMPLX(0.0_dp, 0.0_dp, KIND=dp))
            CALL cp_cfm_to_cfm(projected, cmo_coeff, nactive)
            CALL cp_cfm_release(projected)
            CALL cp_fm_struct_release(active_struct)
         END DO

         CALL cp_cfm_release(matrix_s)
         CALL cp_cfm_release(matrix_p)
      END DO

   END SUBROUTINE qs_kpoint_state_prepare_fixed_density

! **************************************************************************************************
!> \brief Reconstruct natural orbitals from the saved physical metric and density.
!>        The leading eigenpairs of S P S C = S C f retain the most occupied states.
!> \param density Original density kernel.
!> \param metric Physical normalization matrix.
!> \param coefficients Output metric-normalized natural orbitals, most occupied first.
!> \param occupations Occupation eigenvalues in the same order.
! **************************************************************************************************
   SUBROUTINE natural_orbitals(density, metric, coefficients, occupations)
      TYPE(cp_cfm_type), INTENT(IN)                      :: density, metric
      TYPE(cp_cfm_type), INTENT(INOUT)                   :: coefficients
      REAL(dp), INTENT(OUT)                              :: occupations(:)

      INTEGER                                            :: n
      TYPE(cp_cfm_type)                                  :: operator, overlap, work

      n = metric%matrix_struct%nrow_global
      CALL cp_cfm_create(operator, metric%matrix_struct, name="natural-orbital operator")
      CALL cp_cfm_create(work, metric%matrix_struct, name="natural-orbital work")
      CALL cp_cfm_create(overlap, metric%matrix_struct, name="natural-orbital overlap work")
      CALL parallel_gemm('N', 'N', n, n, n, z_one, density, metric, z_zero, work)
      CALL parallel_gemm('N', 'N', n, n, n, z_one, metric, work, z_zero, operator)
      ! Request the lowest eigenpairs of the negative operator.
      CALL cp_cfm_scale(-z_one, operator)
      CALL cp_cfm_to_cfm(metric, overlap)
      CALL cp_cfm_geeig(operator, overlap, coefficients, occupations, work, lowest_subset=.TRUE.)
      occupations(:) = -occupations
      CALL cp_cfm_release(overlap)
      CALL cp_cfm_release(work)
      CALL cp_cfm_release(operator)
   END SUBROUTINE natural_orbitals

! **************************************************************************************************
!> \brief Distance between the saved and trial densities in the physical charge metric.
!> \param kpoints Complete native densities and the trial from the current Hamiltonian.
!> \return Brillouin-zone weighted Frobenius norm in an orthonormal representation.
! **************************************************************************************************
   FUNCTION thermal_restart_distance(kpoints) RESULT(distance)
      TYPE(kpoint_type), INTENT(IN)                      :: kpoints
      REAL(dp)                                           :: distance

      COMPLEX(dp)                                       :: contraction
      INTEGER                                            :: ik, local, n, spin
      TYPE(cp_cfm_type)                                  :: difference, weighted, work

      CPASSERT(kpoints%native_trial_valid)
      distance = 0.0_dp
      DO local = 1, SIZE(kpoints%kp_env)
         ik = kpoints%kp_env(local)%kpoint_env%nkpoint
         ASSOCIATE (metric => kpoints%native_operators(kpoint_slot_s)%matrices(local, 1))
            n = metric%matrix_struct%nrow_global
            CALL cp_cfm_create(difference, metric%matrix_struct)
            CALL cp_cfm_create(work, metric%matrix_struct)
            CALL cp_cfm_create(weighted, metric%matrix_struct)
            DO spin = 1, SIZE(kpoints%native_density, 2)
               difference%local_data(:, :) = kpoints%native_trial(local, spin)%local_data - &
                                             kpoints%native_density(local, spin)%local_data
               CALL parallel_gemm('N', 'N', n, n, n, z_one, metric, difference, z_zero, work)
               CALL parallel_gemm('N', 'N', n, n, n, z_one, work, metric, z_zero, weighted)
               CALL cp_cfm_trace(difference, weighted, contraction, conjugate_a=.TRUE.)
               IF (kpoints%para_env_kp%is_source()) distance = distance + kpoints%wkp(ik)*REAL(contraction, dp)
            END DO
            CALL cp_cfm_release(weighted)
            CALL cp_cfm_release(work)
            CALL cp_cfm_release(difference)
         END ASSOCIATE
      END DO
      CALL kpoints%para_env%sum(distance)
      distance = SQRT(MAX(0.0_dp, distance))
   END FUNCTION thermal_restart_distance

! **************************************************************************************************
!> \brief Recover the initial Fermi state of a saved native Bloch density.
!>        Keep the trial orbitals when the frozen Hamiltonian already reproduces the
!>        saved density to the physical tolerance. Otherwise verify the inverse-Fermi
!>        labels in the original particle-number groups before installing natural orbitals.
!> \param kpoints Saved native density and its current physical metric.
!> \param smear Fermi temperature and spin filling settings.
!> \param restricted Whether two spin channels must share orbitals.
!> \param eps_scf Physical convergence tolerance in Hartree.
! **************************************************************************************************
   SUBROUTINE qs_kpoint_state_prepare_thermal_restart(kpoints, smear, restricted, eps_scf)
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(smear_type), INTENT(INOUT)                    :: smear
      LOGICAL, INTENT(IN)                                :: restricted
      REAL(dp), INTENT(IN)                               :: eps_scf

      COMPLEX(dp)                                       :: charge
      INTEGER                                            :: ik, invalid, local, nao, nmo, nspin, spin
      INTEGER                                            :: electrons(2)
      REAL(dp)                                           :: kts, maxocc, mu, tolerance
      REAL(dp)                                           :: targets(2)
      REAL(dp), ALLOCATABLE                              :: energies(:, :, :), fraction(:), refilled(:, :, :), &
                                                            saved(:, :, :), values(:)
      TYPE(cp_cfm_type), ALLOCATABLE                      :: coefficients(:, :)
      TYPE(cp_fm_struct_type), POINTER                    :: structure
      TYPE(kpoint_env_type), POINTER                     :: kp

      IF (smear%method /= smear_fermi_dirac .OR. smear%electronic_temperature <= 0.0_dp) RETURN
      CPASSERT(kpoints%native_density_seeded)
      kp => kpoints%kp_env(1)%kpoint_env
      nspin = SIZE(kp%mos)
      IF (restricted .AND. nspin > 1) RETURN
      nao = kp%mos(1)%nao
      nmo = kp%mos(1)%nmo
      IF (nmo == 0) RETURN
      ! Do not reconstruct very small occupations from rounded density eigenvalues
      ! when a physical diagonalization already recovers the saved state.
      CALL kpoint_set_mo_occupation(kpoints, smear, separate_spin_occupations=restricted)
      CALL kpoint_density_matrices(kpoints)
      CALL kpoint_density_state_capture(kpoints)
      IF (smear%electronic_temperature*thermal_restart_distance(kpoints) <= eps_scf) RETURN
      maxocc = 2.0_dp/REAL(nspin, dp)
      tolerance = MAX(128.0_dp*EPSILON(1.0_dp)*nao, smear%eps_fermi_dirac)
      electrons = 0
      DO spin = 1, nspin
         electrons(spin) = kp%mos(spin)%nelectron
      END DO
      targets = smearing_spin_targets(electrons, smear%fixed_mag_mom)
      IF (nspin == 1) targets(1) = REAL(electrons(1), dp)
      ALLOCATE (coefficients(SIZE(kpoints%kp_env), nspin), values(nmo), fraction(nmo), &
                saved(nmo, kpoints%nkp, nspin), energies(nmo, kpoints%nkp, nspin), &
                refilled(nmo, kpoints%nkp, nspin))
      saved = 0.0_dp
      energies = 0.0_dp
      invalid = 0
      DO spin = 1, nspin
         DO local = 1, SIZE(kpoints%kp_env)
            kp => kpoints%kp_env(local)%kpoint_env
            ik = kp%nkpoint
            CPASSERT(kp%mos(spin)%nmo == nmo)
            IF (kpoints%use_real_wfn) THEN
               structure => kp%mos(spin)%mo_coeff%matrix_struct
            ELSE
               structure => kp%mos(spin)%cmo_coeff%matrix_struct
            END IF
            CALL cp_cfm_create(coefficients(local, spin), structure, name="thermal restart natural orbitals")
            ASSOCIATE (density => kpoints%native_density(local, spin), &
                       metric => kpoints%native_operators(kpoint_slot_s)%matrices(local, 1))
               CALL cp_cfm_trace(density, metric, charge, conjugate_a=.TRUE.)
               CALL natural_orbitals(density, metric, coefficients(local, spin), values)
            END ASSOCIATE
            IF (.NOT. ALL(ieee_is_finite(values)) .OR. .NOT. ieee_is_finite(REAL(charge, dp))) THEN
               invalid = 1
               CYCLE
            END IF
            IF (MINVAL(values) < -tolerance .OR. MAXVAL(values) > maxocc + tolerance .OR. &
                ABS(REAL(charge, dp) - SUM(values)) > tolerance) invalid = 1
            IF (kpoints%use_real_wfn) THEN
               ASSOCIATE (data => coefficients(local, spin)%local_data)
                  IF (MAXVAL(ABS(AIMAG(data))) > 128.0_dp*EPSILON(1.0_dp)*MAX(1.0_dp, MAXVAL(ABS(data)))) invalid = 1
               END ASSOCIATE
            END IF
            saved(:, ik, spin) = values
            ! Finite labels at exactly empty or full endpoints change f only at roundoff.
            fraction(:) = MAX(EPSILON(1.0_dp), MIN(1.0_dp - EPSILON(1.0_dp), values/maxocc))
            energies(:, ik, spin) = smear%electronic_temperature*LOG((1.0_dp - fraction)/fraction)
         END DO
      END DO
      CALL kpoints%para_env%max(invalid)
      IF (invalid == 0) THEN
         CALL kpoints%para_env_inter_kp%sum(saved)
         CALL kpoints%para_env_inter_kp%sum(energies)
         IF (nspin == 1 .OR. smearing_separate_spin(smear%fixed_mag_mom)) THEN
            DO spin = 1, nspin
               CALL Smearkp(refilled(:, :, spin), mu, kts, energies(:, :, spin), targets(spin), &
                            kpoints%wkp, smear%electronic_temperature, maxocc, smear_fermi_dirac)
               IF (ABS(DOT_PRODUCT(SUM(saved(:, :, spin), DIM=1), kpoints%wkp) - targets(spin)) > tolerance) invalid = 1
            END DO
         ELSE
            CALL Smearkp2(refilled, mu, kts, energies, REAL(SUM(electrons), dp), kpoints%wkp, &
                          smear%electronic_temperature, smear_fermi_dirac)
            IF (ABS(DOT_PRODUCT(SUM(SUM(saved, DIM=3), DIM=1), kpoints%wkp) - SUM(electrons)) > tolerance) invalid = 1
         END IF
         IF (MAXVAL(ABS(refilled - saved)) > tolerance) invalid = 1
         CALL kpoints%para_env%max(invalid)
      END IF
      DO spin = 1, nspin
         DO local = 1, SIZE(kpoints%kp_env)
            IF (invalid == 0) THEN
               kp => kpoints%kp_env(local)%kpoint_env
               IF (kpoints%use_real_wfn) THEN
                  CALL cp_cfm_to_fm(coefficients(local, spin), kp%mos(spin)%mo_coeff)
               ELSE
                  CALL cp_cfm_to_cfm(coefficients(local, spin), kp%mos(spin)%cmo_coeff)
               END IF
               kp%mos(spin)%eigenvalues(:nmo) = energies(:, kp%nkpoint, spin)
            END IF
            CALL cp_cfm_release(coefficients(local, spin))
         END DO
      END DO
   END SUBROUTINE qs_kpoint_state_prepare_thermal_restart

! **************************************************************************************************
!> \brief Rotates a uniformly occupied k-point subspace into Ritz states of the current
!>        Hamiltonian and stores their energy labels. The Hamiltonian must first be transferred
!>        to kp%ot_hmat. Since the rotation is internal to a uniformly occupied subspace, it does
!>        not change the density represented by the orbitals.
!> \param kpoints K-point environment containing the orbitals and transferred Hamiltonians.
! **************************************************************************************************
   SUBROUTINE qs_kpoint_state_canonicalize_fixed(kpoints)
      TYPE(kpoint_type), POINTER                         :: kpoints

      INTEGER                                            :: ikpoint_local, ispin, nao, nmo, nocc
      LOGICAL                                            :: use_real_wfn
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: ritz_values
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvalues_re
      TYPE(cp_cfm_type)                                  :: coeff, eigenvectors, hblock, hc, &
                                                            matrix_h, rotated
      TYPE(cp_cfm_type), POINTER                         :: cmo_coeff
      TYPE(cp_fm_struct_type), POINTER                   :: active_struct, square_struct
      TYPE(kpoint_env_type), POINTER                     :: kp

      NULLIFY (active_struct, cmo_coeff, eigenvalues_re, kp, &
               square_struct)

      CPASSERT(ASSOCIATED(kpoints))
      CPASSERT(ASSOCIATED(kpoints%kp_env))
      CALL get_kpoint_info(kpoints, use_real_wfn=use_real_wfn)
      IF (use_real_wfn) THEN
         CPABORT("Fixed-density k-point OT canonicalization requires complex wavefunctions.")
      END IF

      DO ikpoint_local = 1, SIZE(kpoints%kp_env)
         kp => kpoints%kp_env(ikpoint_local)%kpoint_env
         CPASSERT(ASSOCIATED(kp))
         CPASSERT(ASSOCIATED(kp%mos))
         CPASSERT(ASSOCIATED(kp%ot_hmat))
         CPASSERT(SIZE(kp%ot_hmat, 1) >= 2)
         CPASSERT(SIZE(kp%ot_hmat, 2) >= SIZE(kp%mos))

         CALL cp_fm_get_info(kp%ot_hmat(1, 1), matrix_struct=square_struct, nrow_global=nao)
         CALL cp_cfm_create(matrix_h, square_struct, name="fixed-density Ritz Hamiltonian")

         DO ispin = 1, SIZE(kp%mos)
            CALL get_mo_set(kp%mos(ispin), homo=nocc, cmo_coeff=cmo_coeff, nmo=nmo, &
                            eigenvalues=eigenvalues_re)
            CPASSERT(nocc >= 0 .AND. nocc <= nmo)
            IF (nocc == 0) CYCLE
            CPASSERT(ASSOCIATED(cmo_coeff))
            CPASSERT(ASSOCIATED(eigenvalues_re))
            CPASSERT(SIZE(eigenvalues_re) >= nocc)

            CALL cp_fm_struct_create(active_struct, template_fmstruct=cmo_coeff%matrix_struct, &
                                     ncol_global=nocc)
            CALL cp_cfm_create(coeff, active_struct, name="fixed-density Ritz orbitals")
            CALL cp_cfm_create(hc, active_struct, name="fixed-density Ritz HC")
            CALL cp_cfm_create(rotated, active_struct, name="fixed-density Ritz rotated orbitals")
            CALL cp_cfm_create(hblock, active_struct, nrow=nocc, ncol=nocc, &
                               name="fixed-density projected Hamiltonian")
            CALL cp_cfm_create(eigenvectors, active_struct, nrow=nocc, ncol=nocc, &
                               name="fixed-density Ritz eigenvectors")
            CALL cp_fm_to_cfm(kp%ot_hmat(1, ispin), kp%ot_hmat(2, ispin), matrix_h)
            IF (ALLOCATED(kpoints%lowdin_v)) CALL lowdin_kp_apply(kp, kpoints%lowdin_v(:, ispin), cmat=matrix_h)
            CALL cp_cfm_to_cfm(cmo_coeff, coeff, nocc)

            CALL parallel_gemm('N', 'N', nao, nocc, nao, z_one, matrix_h, coeff, z_zero, hc)
            CALL parallel_gemm('C', 'N', nocc, nocc, nao, z_one, coeff, hc, z_zero, hblock)
            ALLOCATE (ritz_values(nocc))
            CALL cp_cfm_heevd(hblock, eigenvectors, ritz_values)
            CALL parallel_gemm('N', 'N', nao, nocc, nocc, z_one, coeff, eigenvectors, &
                               z_zero, rotated)
            CALL cp_cfm_to_cfm(rotated, cmo_coeff, nocc)
            eigenvalues_re(1:nocc) = ritz_values
            DEALLOCATE (ritz_values)

            CALL cp_cfm_release(eigenvectors)
            CALL cp_cfm_release(hblock)
            CALL cp_cfm_release(rotated)
            CALL cp_cfm_release(hc)
            CALL cp_cfm_release(coeff)
            CALL cp_fm_struct_release(active_struct)
         END DO

         CALL cp_cfm_release(matrix_h)
      END DO

   END SUBROUTINE qs_kpoint_state_canonicalize_fixed

! **************************************************************************************************
!> \brief Rebuilds the physical density from the current k-point orbitals and occupations.
!> \param qs_env QS environment containing the orbital state to commit.
!> \param update_occupations Recompute occupations from the current MO energies and SCF settings.
!> \param separate_spin_occupations Fill restricted energy channels independently.
!> \param fixed_occupations Restore the fixed occupied rank required by OT.
! **************************************************************************************************
   SUBROUTINE qs_kpoint_state_commit(qs_env, update_occupations, &
                                     separate_spin_occupations, fixed_occupations)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: update_occupations
      LOGICAL, INTENT(IN), OPTIONAL                      :: separate_spin_occupations, &
                                                            fixed_occupations

      INTEGER                                            :: ispin
      LOGICAL                                            :: fixed, separate_spins
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s_kp, rho_ao_kp
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(scf_control_type), POINTER                    :: scf_control

      NULLIFY (energy, kpoints, matrix_s_kp, mos, rho_ao_kp, sab_nl, rho, scf_env, scf_control)

      fixed = .FALSE.
      separate_spins = .FALSE.
      IF (PRESENT(fixed_occupations)) fixed = fixed_occupations
      IF (PRESENT(separate_spin_occupations)) separate_spins = separate_spin_occupations
      IF (fixed .AND. .NOT. update_occupations) THEN
         CPABORT("Fixed occupations can only be requested while updating occupations.")
      END IF

      CALL get_qs_env(qs_env, kpoints=kpoints, matrix_s_kp=matrix_s_kp, &
                      energy=energy, rho=rho, scf_env=scf_env, scf_control=scf_control)
      CALL get_kpoint_info(kpoints, sab_nl=sab_nl)

      IF (update_occupations) THEN
         CALL kpoint_set_mo_occupation( &
            kpoints, scf_control%smear, separate_spin_occupations=separate_spins)
         IF (fixed) CALL qs_kpoint_set_fixed_occupations(kpoints)
      END IF
      CALL kpoint_density_matrices(kpoints)
      IF (kpoint_density_state_active(kpoints)) THEN
         ! A committed orbital state also supplies the complete augmented density.
         ! Its finite real-space image transform is not the variational P(k).
         CALL kpoint_density_state_capture(kpoints)
         CALL kpoint_density_state_accept(kpoints, 1.0_dp)
      END IF
      CALL qs_rho_get(rho, rho_ao_kp=rho_ao_kp)
      CALL kpoint_density_transform(kpoints, rho_ao_kp, .FALSE., &
                                    matrix_s_kp(1, 1)%matrix, sab_nl, scf_env%scf_work1, &
                                    overlap_rs=matrix_s_kp)
      CALL qs_rho_update_rho(rho, qs_env=qs_env)
      CALL qs_ks_did_change(qs_env%ks_env, rho_changed=.TRUE.)

      energy%kTS = 0.0_dp
      energy%efermi = 0.0_dp
      mos => kpoints%kp_env(1)%kpoint_env%mos
      CPASSERT(ASSOCIATED(mos))
      DO ispin = 1, SIZE(mos)
         energy%kTS = energy%kTS + mos(ispin)%kTS
         energy%efermi = energy%efermi + mos(ispin)%mu
      END DO
      energy%efermi = energy%efermi/REAL(SIZE(mos), KIND=dp)

   END SUBROUTINE qs_kpoint_state_commit

! **************************************************************************************************
!> \brief Checks whether every local k-point channel contains nonzero occupied orbitals.
!> \param kpoints K-point environment.
!> \param require_full_space Require every available orbital instead of the occupied block.
!> \return True if every channel contains a nonzero real or imaginary coefficient.
! **************************************************************************************************
   FUNCTION qs_kpoint_mos_initialized(kpoints, require_full_space) RESULT(initialized)
      TYPE(kpoint_type), POINTER                         :: kpoints
      LOGICAL, INTENT(IN), OPTIONAL                      :: require_full_space
      LOGICAL                                            :: initialized

      INTEGER                                            :: ikpoint_local, ispin, nactive, nmo, nocc
      LOGICAL                                            :: full_space
      REAL(KIND=dp)                                      :: coefficient_norm
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: column_norm
      TYPE(cp_cfm_type), POINTER                         :: cmo_coeff
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(kpoint_env_type), POINTER                     :: kp
      TYPE(mp_para_env_type), POINTER                    :: para_env

      NULLIFY (cmo_coeff, kp, mo_coeff, para_env)
      CPASSERT(ASSOCIATED(kpoints))
      CPASSERT(ASSOCIATED(kpoints%kp_env))

      full_space = .FALSE.
      IF (PRESENT(require_full_space)) full_space = require_full_space
      coefficient_norm = HUGE(1.0_dp)
      DO ikpoint_local = 1, SIZE(kpoints%kp_env)
         kp => kpoints%kp_env(ikpoint_local)%kpoint_env
         CPASSERT(ASSOCIATED(kp))
         CPASSERT(ASSOCIATED(kp%mos))
         DO ispin = 1, SIZE(kp%mos)
            CALL get_mo_set(kp%mos(ispin), homo=nocc, cmo_coeff=cmo_coeff, &
                            mo_coeff=mo_coeff, nmo=nmo)
            nactive = MERGE(nmo, nocc, full_space)
            IF (nactive == 0 .AND. .NOT. full_space) CYCLE
            IF (nactive < 1 .OR. nactive > nmo) THEN
               coefficient_norm = 0.0_dp
               CYCLE
            END IF
            ALLOCATE (column_norm(nmo))
            IF (ASSOCIATED(cmo_coeff)) THEN
               CALL cp_cfm_vectorsnorm(cmo_coeff, column_norm)
            ELSE
               CPASSERT(ASSOCIATED(mo_coeff))
               CALL cp_fm_vectorsnorm(mo_coeff, column_norm)
            END IF
            coefficient_norm = MIN(coefficient_norm, MINVAL(column_norm(1:nactive)))
            DEALLOCATE (column_norm)
         END DO
      END DO

      CALL get_kpoint_info(kpoints, para_env=para_env)
      CPASSERT(ASSOCIATED(para_env))
      CALL para_env%min(coefficient_norm)
      initialized = coefficient_norm > 100.0_dp*EPSILON(1.0_dp)

   END FUNCTION qs_kpoint_mos_initialized

! **************************************************************************************************
!> \brief Copies the first k-point spin channel into the remaining channels.
!> \param kpoints K-point environment.
!> \param nspin Number of spin channels to populate.
! **************************************************************************************************
   SUBROUTINE qs_kpoint_copy_spin_mos(kpoints, nspin)
      TYPE(kpoint_type), POINTER                         :: kpoints
      INTEGER, INTENT(IN)                                :: nspin

      INTEGER                                            :: ikpoint_local, ispin
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvalues_source, eigenvalues_target
      TYPE(cp_cfm_type), POINTER                         :: cmo_coeff_source, cmo_coeff_target
      TYPE(cp_fm_type), POINTER                          :: mo_coeff_source, mo_coeff_target
      TYPE(kpoint_env_type), POINTER                     :: kp

      NULLIFY (cmo_coeff_source, cmo_coeff_target, eigenvalues_source, eigenvalues_target, kp, &
               mo_coeff_source, mo_coeff_target)
      CPASSERT(ASSOCIATED(kpoints))
      CPASSERT(ASSOCIATED(kpoints%kp_env))
      CPASSERT(nspin >= 1)

      DO ikpoint_local = 1, SIZE(kpoints%kp_env)
         kp => kpoints%kp_env(ikpoint_local)%kpoint_env
         CPASSERT(ASSOCIATED(kp))
         CPASSERT(ASSOCIATED(kp%mos))
         CPASSERT(SIZE(kp%mos) >= nspin)
         IF (ASSOCIATED(kp%mos(1)%cmo_coeff)) THEN
            cmo_coeff_source => kp%mos(1)%cmo_coeff
            CALL get_mo_set(kp%mos(1), eigenvalues=eigenvalues_source)
            CPASSERT(ASSOCIATED(eigenvalues_source))
            DO ispin = 2, nspin
               cmo_coeff_target => kp%mos(ispin)%cmo_coeff
               CALL get_mo_set(kp%mos(ispin), eigenvalues=eigenvalues_target)
               CPASSERT(ASSOCIATED(cmo_coeff_target))
               CPASSERT(ASSOCIATED(eigenvalues_target))
               CALL cp_cfm_to_cfm(cmo_coeff_source, cmo_coeff_target)
               eigenvalues_target(:) = eigenvalues_source(:)
            END DO
         ELSE
            CALL get_mo_set(kp%mos(1), mo_coeff=mo_coeff_source, eigenvalues=eigenvalues_source)
            CPASSERT(ASSOCIATED(mo_coeff_source))
            CPASSERT(ASSOCIATED(eigenvalues_source))
            DO ispin = 2, nspin
               CALL get_mo_set(kp%mos(ispin), mo_coeff=mo_coeff_target, &
                               eigenvalues=eigenvalues_target)
               CPASSERT(ASSOCIATED(mo_coeff_target))
               CPASSERT(ASSOCIATED(eigenvalues_target))
               CALL cp_fm_to_fm(mo_coeff_source, mo_coeff_target)
               eigenvalues_target(:) = eigenvalues_source(:)
            END DO
         END IF
      END DO

   END SUBROUTINE qs_kpoint_copy_spin_mos

! **************************************************************************************************
!> \brief Restores the fixed occupied rank required by k-point OT.
!> \param kpoints K-point environment.
! **************************************************************************************************
   SUBROUTINE qs_kpoint_set_fixed_occupations(kpoints)
      TYPE(kpoint_type), POINTER                         :: kpoints

      INTEGER                                            :: ikpoint_local, ispin, nelectron, nmo, &
                                                            nocc
      REAL(KIND=dp)                                      :: maxocc
      REAL(KIND=dp), DIMENSION(:), POINTER               :: occupation_numbers
      TYPE(kpoint_env_type), POINTER                     :: kp

      NULLIFY (kp, occupation_numbers)
      CPASSERT(ASSOCIATED(kpoints))
      CPASSERT(ASSOCIATED(kpoints%kp_env))

      DO ikpoint_local = 1, SIZE(kpoints%kp_env)
         kp => kpoints%kp_env(ikpoint_local)%kpoint_env
         CPASSERT(ASSOCIATED(kp))
         CPASSERT(ASSOCIATED(kp%mos))
         DO ispin = 1, SIZE(kp%mos)
            CALL get_mo_set(kp%mos(ispin), maxocc=maxocc, nelectron=nelectron, &
                            nmo=nmo, occupation_numbers=occupation_numbers)
            CPASSERT(ASSOCIATED(occupation_numbers))
            ! The fixed rank comes from the spin electron count, not a previous k-point filling.
            nocc = NINT(REAL(nelectron, KIND=dp)/maxocc)
            IF (nocc < 0 .OR. nocc > nmo) THEN
               CPABORT("Fixed-occupation k-point OT rank exceeds the available orbital space.")
            END IF
            IF (ABS(maxocc*REAL(nocc, KIND=dp) - REAL(nelectron, KIND=dp)) > &
                100.0_dp*EPSILON(1.0_dp)) THEN
               CPABORT("Fixed-occupation k-point OT requires uniform integer occupations.")
            END IF
            occupation_numbers(:) = 0.0_dp
            occupation_numbers(1:nocc) = maxocc
            kp%mos(ispin)%homo = nocc
            kp%mos(ispin)%kTS = 0.0_dp
         END DO
      END DO

   END SUBROUTINE qs_kpoint_set_fixed_occupations

END MODULE qs_kpoint_state
