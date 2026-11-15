!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculates the energy contribution and the mo_derivative of
!>        a static electric field (nonperiodic)
!> \par History
!>      none
!> \author JGH (05.2015)
! **************************************************************************************************
MODULE qs_efield_local
   USE ai_moments,                      ONLY: dipole_force
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind,&
                                              get_atomic_kind_set
   USE basis_set_types,                 ONLY: gto_basis_set_p_type,&
                                              gto_basis_set_type
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_add,&
                                              dbcsr_copy,&
                                              dbcsr_get_block_p,&
                                              dbcsr_p_type,&
                                              dbcsr_set
   USE cp_dbcsr_contrib,                ONLY: dbcsr_dot
   USE kinds,                           ONLY: dp
   USE kpoint_methods,                  ONLY: kpoint_star_operator_defect
   USE mathconstants,                   ONLY: fourpi
   USE mathlib,                         ONLY: invmat
   USE message_passing,                 ONLY: mp_para_env_type
   USE orbital_pointers,                ONLY: ncoset
   USE particle_types,                  ONLY: particle_type
   USE pw_types,                        ONLY: pw_r3d_rs_type
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type,&
                                              set_qs_env
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type
   USE qs_moments,                      ONLY: build_local_moment_matrix
   USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                              neighbor_list_iterate,&
                                              neighbor_list_iterator_create,&
                                              neighbor_list_iterator_p_type,&
                                              neighbor_list_iterator_release,&
                                              neighbor_list_set_p_type
   USE qs_period_efield_types,          ONLY: efield_berry_type,&
                                              init_efield_matrices,&
                                              set_efield_matrices
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_upf_inverse_coulomb,          ONLY: qs_upf_inverse_local_field
   USE qs_upf_inverse_scf,              ONLY: qs_upf_inverse_scf_active,&
                                              qs_upf_inverse_scf_integrate_background
   USE qs_upf_inverse_types,            ONLY: qs_upf_inverse_scf_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_efield_local'

   ! *** Public subroutines ***

   PUBLIC :: qs_efield_local_operator, qs_efield_dipole_force, qs_efield_upf_evaluate

! **************************************************************************************************

CONTAINS

! **************************************************************************************************

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param just_energy ...
!> \param calculate_forces ...
! **************************************************************************************************
   SUBROUTINE qs_efield_local_operator(qs_env, just_energy, calculate_forces)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: just_energy, calculate_forces

      CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_efield_local_operator'

      COMPLEX(dp), ALLOCATABLE                           :: operators(:, :, :)
      INTEGER                                            :: handle, spin
      LOGICAL                                            :: s_mstruct_changed
      REAL(dp), DIMENSION(3)                             :: rpoint
      TYPE(dbcsr_p_type), POINTER                        :: matrix_ks(:)
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_ks_env_type), POINTER                      :: ks_env

      CALL timeset(routineN, handle)

      NULLIFY (dft_control)
      CALL get_qs_env(qs_env, s_mstruct_changed=s_mstruct_changed, &
                      dft_control=dft_control, ks_env=ks_env)

      IF (dft_control%apply_efield) THEN
         IF (qs_upf_inverse_scf_active(ks_env)) THEN
            CALL get_qs_env(qs_env, energy=energy, matrix_ks=matrix_ks)
            CALL qs_efield_upf_evaluate(qs_env, energy%efield, operators)
            IF (.NOT. just_energy) THEN
               DO spin = 1, SIZE(matrix_ks)
                  CALL qs_upf_inverse_scf_integrate_background(ks_env, operators, matrix_ks(spin)%matrix, &
                                                               spin, additive=.TRUE.)
               END DO
            END IF
            ! The common inverse force callback owns the complete electronic and ionic response.
         ELSE
            rpoint = 0.0_dp
            IF (s_mstruct_changed) CALL qs_efield_integrals(qs_env, rpoint)
            CALL qs_efield_mo_derivatives(qs_env, rpoint, just_energy, calculate_forces)
         END IF
      END IF

      CALL timestop(handle)

   END SUBROUTINE qs_efield_local_operator

! **************************************************************************************************
!> \brief Couple a nonperiodic electric field to the complete UPF charge and its geometry response.
!>        The electronic map includes reconstruction, augmentation and frozen PAW core contributions.
!>        The point countercharge uses the same tabulated frozen-core charge exactly once.
!> \param qs_env Native environment with a prepared physical UPF density.
!> \param energy Electronic and ionic electric-field energy.
!> \param operators Optional valence AO operators, one matrix per physical Bloch slice.
!> \param response Optional total directional derivative at fixed AO density coefficients.
!> \param atom_response Atomic Cartesian displacements.
!> \param cell_response Cell-matrix displacement, including the moving sampling grid.
!> \param coefficients_response Response of the inverse reconstruction coefficients.
! **************************************************************************************************
   SUBROUTINE qs_efield_upf_evaluate(qs_env, energy, operators, response, atom_response, &
                                     cell_response, coefficients_response)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(dp), INTENT(OUT)                              :: energy
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT), OPTIONAL    :: operators(:, :, :)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: response
      REAL(dp), INTENT(IN), OPTIONAL                     :: atom_response(:, :), cell_response(3, 3)
      COMPLEX(dp), INTENT(IN), OPTIONAL                  :: coefficients_response(:, :, :)

      COMPLEX(dp), ALLOCATABLE                           :: field_operator(:, :, :)
      INTEGER                                            :: atom, i, info, j, k, kind
      LOGICAL                                            :: check_symmetry, little_group
      REAL(dp)                                           :: charge, defect, fieldpol(3), &
                                                            inverse_cell(3, 3), point(3), scale, &
                                                            strain(3, 3), tolerance
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: control
      TYPE(particle_type), POINTER                       :: atoms(:)
      TYPE(pw_r3d_rs_type)                               :: field, field_response
      TYPE(pw_r3d_rs_type), POINTER                      :: density(:)
      TYPE(qs_kind_type), POINTER                        :: kinds(:)
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_upf_inverse_scf_type), POINTER             :: state

      energy = 0.0_dp
      IF (PRESENT(response)) response = 0.0_dp
      CALL get_qs_env(qs_env, dft_control=control, ks_env=ks_env)
      IF (.NOT. control%apply_efield) RETURN
      IF (.NOT. qs_upf_inverse_scf_active(ks_env)) RETURN
      CALL get_ks_env(ks_env, upf_inverse=state)
      CALL get_qs_env(qs_env, rho=rho, cell=cell, particle_set=atoms, qs_kind_set=kinds)
      CALL qs_rho_get(rho, rho_r=density)
      fieldpol = control%efield_fields(1)%efield%polarisation*control%efield_fields(1)%efield%strength
      CALL field%create(density(1)%pw_grid)
      CALL field_response%create(density(1)%pw_grid)
      field_response%array = 0.0_dp
      strain = 0.0_dp
      IF (PRESENT(response)) THEN
         CPASSERT(PRESENT(atom_response))
         CPASSERT(PRESENT(cell_response))
         CPASSERT(PRESENT(coefficients_response))
         inverse_cell = cell%hmat
         CALL invmat(inverse_cell, info)
         CPASSERT(info == 0)
         strain = MATMUL(cell_response, inverse_cell)
      END IF
      DO k = field%pw_grid%bounds_local(1, 3), field%pw_grid%bounds_local(2, 3)
         DO j = field%pw_grid%bounds_local(1, 2), field%pw_grid%bounds_local(2, 2)
            DO i = field%pw_grid%bounds_local(1, 1), field%pw_grid%bounds_local(2, 1)
               point = MATMUL(field%pw_grid%dh, REAL([i, j, k] - field%pw_grid%bounds(1, :), dp))
               field%array(i, j, k) = SUM(fieldpol*point)
               IF (PRESENT(response)) field_response%array(i, j, k) = SUM(fieldpol*MATMUL(strain, point))
            END DO
         END DO
      END DO
      CALL qs_upf_inverse_local_field(qs_env, field, energy, field_operator, response, atom_response, &
                                      cell_response, coefficients_response, potential_response=field_response)
      DO atom = 1, SIZE(atoms)
         CALL get_atomic_kind(atoms(atom)%atomic_kind, kind_number=kind)
         IF (kinds(kind)%ghost) CYCLE
         CALL get_qs_kind(kinds(kind), zeff=charge)
         IF (state%paw(kind)%active) THEN
            charge = charge + fourpi*SUM(state%paw(kind)%xcgrid%wr*state%paw(kind)%xcgrid%core(:, 1))
         END IF
         energy = energy - charge*SUM(fieldpol*atoms(atom)%r)
         IF (PRESENT(response)) response = response - charge*SUM(fieldpol*atom_response(:, atom))
      END DO
      check_symmetry = ASSOCIATED(state%bloch%kpoints)
      little_group = .NOT. control%use_gauxc
      IF (check_symmetry) THEN
         check_symmetry = .NOT. state%bloch%kpoints%full_grid
         IF (SIZE(state%bloch%weights) == state%bloch%kpoints%nkp .AND. .NOT. little_group) check_symmetry = .FALSE.
      END IF
      IF (check_symmetry) THEN
         CALL kpoint_star_operator_defect(state%bloch%star, field_operator, defect, scale, little_group=little_group)
         tolerance = MAX(1024.0_dp*EPSILON(1.0_dp), 32.0_dp*control%qs_control%eps_rho_gspace)*MAX(1.0_dp, scale)
         IF (defect > tolerance) THEN
            CALL cp_abort(__LOCATION__, &
                          "The electric field is incompatible with the k-point symmetry reduction. "// &
                          "Use FULL_GRID TRUE and SYMMETRY FALSE.")
         END IF
      END IF
      IF (PRESENT(operators)) CALL MOVE_ALLOC(field_operator, operators)
      CALL field%release()
      CALL field_response%release()
   END SUBROUTINE qs_efield_upf_evaluate

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param rpoint ...
! **************************************************************************************************
   SUBROUTINE qs_efield_integrals(qs_env, rpoint)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: rpoint

      CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_efield_integrals'

      INTEGER                                            :: handle, i
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: dipmat, matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(efield_berry_type), POINTER                   :: efield

      CALL timeset(routineN, handle)
      CPASSERT(ASSOCIATED(qs_env))

      CALL get_qs_env(qs_env=qs_env, dft_control=dft_control)
      NULLIFY (matrix_s)
      CALL get_qs_env(qs_env=qs_env, efield=efield, matrix_s=matrix_s)
      CALL init_efield_matrices(efield)
      ALLOCATE (dipmat(3))
      DO i = 1, 3
         ALLOCATE (dipmat(i)%matrix)
         CALL dbcsr_copy(dipmat(i)%matrix, matrix_s(1)%matrix, 'DIP MAT')
         CALL dbcsr_set(dipmat(i)%matrix, 0.0_dp)
      END DO
      CALL build_local_moment_matrix(qs_env, dipmat, 1, rpoint)
      CALL set_efield_matrices(efield=efield, dipmat=dipmat)
      CALL set_qs_env(qs_env=qs_env, efield=efield)
      CALL timestop(handle)

   END SUBROUTINE qs_efield_integrals

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param rpoint ...
!> \param just_energy ...
!> \param calculate_forces ...
! **************************************************************************************************
   SUBROUTINE qs_efield_mo_derivatives(qs_env, rpoint, just_energy, calculate_forces)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rpoint
      LOGICAL                                            :: just_energy, calculate_forces

      CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_efield_mo_derivatives'

      INTEGER                                            :: handle, ia, iatom, idir, ikind, ispin, &
                                                            natom
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind
      REAL(dp)                                           :: charge, ci(3), ener_field, fieldpol(3), &
                                                            tmp
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: fdip
      REAL(dp), DIMENSION(3)                             :: ria
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: dipmat, matrix_ks, matrix_p
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(efield_berry_type), POINTER                   :: efield
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_rho_type), POINTER                         :: rho

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, dft_control=dft_control, cell=cell, particle_set=particle_set)
      CALL get_qs_env(qs_env=qs_env, qs_kind_set=qs_kind_set, &
                      efield=efield, energy=energy, para_env=para_env)

      fieldpol = dft_control%efield_fields(1)%efield%polarisation* &
                 dft_control%efield_fields(1)%efield%strength

      ! nuclear contribution
      natom = SIZE(particle_set)
      IF (calculate_forces) THEN
         CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set, force=force)
         CALL get_atomic_kind_set(atomic_kind_set, atom_of_kind=atom_of_kind)
      END IF
      ci = 0.0_dp
      DO ia = 1, natom
         CALL get_atomic_kind(particle_set(ia)%atomic_kind, kind_number=ikind)
         CALL get_qs_kind(qs_kind_set(ikind), core_charge=charge)
         ria = particle_set(ia)%r - rpoint
         ria = pbc(ria, cell)
         ci(:) = ci(:) + charge*ria(:)
         IF (calculate_forces) THEN
            IF (para_env%mepos == 0) THEN
               iatom = atom_of_kind(ia)
               DO idir = 1, 3
                  force(ikind)%efield(idir, iatom) = force(ikind)%efield(idir, iatom) - fieldpol(idir)*charge
               END DO
            END IF
         END IF
      END DO
      ener_field = -SUM(ci(:)*fieldpol(:))

      ! Energy
      dipmat => efield%dipmat
      NULLIFY (rho, matrix_p)
      CALL get_qs_env(qs_env=qs_env, rho=rho)
      CALL qs_rho_get(rho, rho_ao=matrix_p)
      DO ispin = 1, SIZE(matrix_p)
         DO idir = 1, 3
            CALL dbcsr_dot(matrix_p(ispin)%matrix, dipmat(idir)%matrix, tmp)
            ener_field = ener_field + fieldpol(idir)*tmp
         END DO
      END DO
      energy%efield = ener_field

      IF (.NOT. just_energy) THEN

         ! Update KS matrix
         NULLIFY (matrix_ks)
         CALL get_qs_env(qs_env=qs_env, matrix_ks=matrix_ks)
         DO ispin = 1, SIZE(matrix_ks)
            DO idir = 1, 3
               CALL dbcsr_add(matrix_ks(ispin)%matrix, dipmat(idir)%matrix, &
                              alpha_scalar=1.0_dp, beta_scalar=fieldpol(idir))
            END DO
         END DO

         ! forces from the efield contribution
         IF (calculate_forces) THEN
            ALLOCATE (fdip(3, natom))
            fdip = 0.0_dp
            CALL qs_efield_dipole_force(qs_env, matrix_p, fieldpol, rpoint, fdip)
            DO ia = 1, natom
               CALL get_atomic_kind(particle_set(ia)%atomic_kind, kind_number=ikind)
               iatom = atom_of_kind(ia)
               force(ikind)%efield(1:3, iatom) = force(ikind)%efield(1:3, iatom) + fdip(1:3, ia)
            END DO
            DEALLOCATE (fdip)
         END IF

      END IF

      IF (calculate_forces) THEN
         DO ikind = 1, SIZE(atomic_kind_set)
            CALL para_env%sum(force(ikind)%efield)
         END DO
      END IF

      CALL timestop(handle)

   END SUBROUTINE qs_efield_mo_derivatives

! **************************************************************************************************
!> \brief Electronic gradient from the nuclear derivative of the dipole integrals,
!>        dE/dR = sum_d field_d Tr(P d<a|r_d - rpoint|b>/dR)
!> \param qs_env ...
!> \param matrix_p density matrix (all spins)
!> \param field electric field
!> \param rpoint reference point of the dipole operator
!> \param fdip gradient (3,natom) to which the local contribution is added (not summed over ranks)
! **************************************************************************************************
   SUBROUTINE qs_efield_dipole_force(qs_env, matrix_p, field, rpoint, fdip)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_p
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: field, rpoint
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: fdip

      CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_efield_dipole_force'

      INTEGER                                            :: handle, i, iatom, icol, idir, ikind, &
                                                            irow, iset, jatom, jkind, jset, ldab, &
                                                            ncoa, ncob, nkind, nseta, nsetb, sgfa, &
                                                            sgfb
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, lb_max, lb_min, npgfa, &
                                                            npgfb, nsgfa, nsgfb
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, first_sgfb
      LOGICAL                                            :: found, trans
      REAL(dp)                                           :: dab, fdir
      REAL(dp), DIMENSION(3)                             :: ra, rab, rac, rbc
      REAL(dp), DIMENSION(3, 3)                          :: forcea, forceb
      REAL(dp), DIMENSION(:, :), POINTER                 :: p_block_a, p_block_b, pblock, pmat, work
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_a, set_radius_b
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: rpgfa, rpgfb, sphi_a, sphi_b, zeta, zetb
      TYPE(cell_type), POINTER                           :: cell
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_list
      TYPE(gto_basis_set_type), POINTER                  :: basis_set_a, basis_set_b
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_kind_type), POINTER                        :: qs_kind

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, cell=cell, particle_set=particle_set, qs_kind_set=qs_kind_set, &
                      sab_orb=sab_orb)

      nkind = SIZE(qs_kind_set)

      ALLOCATE (basis_set_list(nkind))
      DO ikind = 1, nkind
         qs_kind => qs_kind_set(ikind)
         CALL get_qs_kind(qs_kind=qs_kind, basis_set=basis_set_a)
         IF (ASSOCIATED(basis_set_a)) THEN
            basis_set_list(ikind)%gto_basis_set => basis_set_a
         ELSE
            NULLIFY (basis_set_list(ikind)%gto_basis_set)
         END IF
      END DO
      !
      CALL neighbor_list_iterator_create(nl_iterator, sab_orb)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, &
                                iatom=iatom, jatom=jatom, r=rab)
         basis_set_a => basis_set_list(ikind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_a)) CYCLE
         basis_set_b => basis_set_list(jkind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_b)) CYCLE
         ! basis ikind
         first_sgfa => basis_set_a%first_sgf
         la_max => basis_set_a%lmax
         la_min => basis_set_a%lmin
         npgfa => basis_set_a%npgf
         nseta = basis_set_a%nset
         nsgfa => basis_set_a%nsgf_set
         rpgfa => basis_set_a%pgf_radius
         set_radius_a => basis_set_a%set_radius
         sphi_a => basis_set_a%sphi
         zeta => basis_set_a%zet
         ! basis jkind
         first_sgfb => basis_set_b%first_sgf
         lb_max => basis_set_b%lmax
         lb_min => basis_set_b%lmin
         npgfb => basis_set_b%npgf
         nsetb = basis_set_b%nset
         nsgfb => basis_set_b%nsgf_set
         rpgfb => basis_set_b%pgf_radius
         set_radius_b => basis_set_b%set_radius
         sphi_b => basis_set_b%sphi
         zetb => basis_set_b%zet

         ra(:) = particle_set(iatom)%r(:) - rpoint(:)
         rac(:) = pbc(ra(:), cell)
         rbc(:) = rac(:) + rab(:)
         dab = SQRT(rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3))

         IF (iatom <= jatom) THEN
            irow = iatom
            icol = jatom
            trans = .FALSE.
         ELSE
            irow = jatom
            icol = iatom
            trans = .TRUE.
         END IF

         fdir = 2.0_dp
         IF (iatom == jatom .AND. dab < 1.e-10_dp) fdir = 1.0_dp

         ! density matrix
         NULLIFY (p_block_a)
         CALL dbcsr_get_block_p(matrix_p(1)%matrix, irow, icol, p_block_a, found)
         IF (.NOT. found) CYCLE
         IF (SIZE(matrix_p) > 1) THEN
            NULLIFY (p_block_b)
            CALL dbcsr_get_block_p(matrix_p(2)%matrix, irow, icol, p_block_b, found)
            CPASSERT(found)
         END IF
         forcea = 0.0_dp
         forceb = 0.0_dp

         DO iset = 1, nseta
            ncoa = npgfa(iset)*ncoset(la_max(iset))
            sgfa = first_sgfa(1, iset)
            DO jset = 1, nsetb
               IF (set_radius_a(iset) + set_radius_b(jset) < dab) CYCLE
               ncob = npgfb(jset)*ncoset(lb_max(jset))
               sgfb = first_sgfb(1, jset)
               ! Calculate the primitive integrals (da|O|b) and (a|O|db)
               ldab = MAX(ncoa, ncob)
               ALLOCATE (work(ldab, ldab), pmat(ncoa, ncob))
               ! Decontract P matrix block
               pmat = 0.0_dp
               DO i = 1, SIZE(matrix_p)
                  IF (i == 1) THEN
                     pblock => p_block_a
                  ELSE
                     pblock => p_block_b
                  END IF
                  IF (trans) THEN
                     CALL dgemm("N", "T", ncoa, nsgfb(jset), nsgfa(iset), &
                                1.0_dp, sphi_a(1, sgfa), SIZE(sphi_a, 1), &
                                pblock(sgfb, sgfa), SIZE(pblock, 1), &
                                0.0_dp, work(1, 1), ldab)
                  ELSE
                     CALL dgemm("N", "N", ncoa, nsgfb(jset), nsgfa(iset), &
                                1.0_dp, sphi_a(1, sgfa), SIZE(sphi_a, 1), &
                                pblock(sgfa, sgfb), SIZE(pblock, 1), &
                                0.0_dp, work(1, 1), ldab)
                  END IF
                  CALL dgemm("N", "T", ncoa, ncob, nsgfb(jset), &
                             1.0_dp, work(1, 1), ldab, &
                             sphi_b(1, sgfb), SIZE(sphi_b, 1), &
                             1.0_dp, pmat(1, 1), ncoa)
               END DO

               CALL dipole_force(la_max(iset), npgfa(iset), zeta(:, iset), rpgfa(:, iset), la_min(iset), &
                                 lb_max(jset), npgfb(jset), zetb(:, jset), rpgfb(:, jset), lb_min(jset), &
                                 1, rac, rbc, pmat, forcea, forceb)

               DEALLOCATE (work, pmat)
            END DO
         END DO

         DO idir = 1, 3
            fdip(1:3, iatom) = fdip(1:3, iatom) + fdir*field(idir)*forcea(idir, 1:3)
            fdip(1:3, jatom) = fdip(1:3, jatom) + fdir*field(idir)*forceb(idir, 1:3)
         END DO

      END DO
      CALL neighbor_list_iterator_release(nl_iterator)
      DEALLOCATE (basis_set_list)

      CALL timestop(handle)

   END SUBROUTINE qs_efield_dipole_force

END MODULE qs_efield_local
