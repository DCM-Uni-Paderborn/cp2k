!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Restart file for k point calculations
!> \par History
!> \author JGH (30.09.2015)
! **************************************************************************************************
MODULE kpoint_io

   USE atomic_kind_types,               ONLY: get_atomic_kind
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
   USE cell_types,                      ONLY: cell_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_get_info,&
                                              dbcsr_p_type,&
                                              dbcsr_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_submatrix,&
                                              cp_fm_read_unformatted,&
                                              cp_fm_release,&
                                              cp_fm_set_submatrix,&
                                              cp_fm_type,&
                                              cp_fm_write_unformatted
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type,&
                                              cp_to_string
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE input_section_types,             ONLY: section_vals_type
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_type
   USE message_passing,                 ONLY: mp_para_env_type
   USE orbital_pointers,                ONLY: nso
   USE particle_types,                  ONLY: particle_type
   USE qs_dftb_types,                   ONLY: qs_dftb_atom_type
   USE qs_dftb_utils,                   ONLY: get_dftb_atom_param
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_mo_io,                        ONLY: wfn_restart_file_name
   USE qs_ot_types,                     ONLY: qs_ot_check_channel_context,&
                                              qs_ot_restart_signature_nchar,&
                                              qs_ot_restart_signature_nint,&
                                              qs_ot_restart_signature_nlog,&
                                              qs_ot_restart_signature_nreal,&
                                              qs_ot_settings_restart_compatible,&
                                              qs_ot_settings_restart_signature,&
                                              qs_ot_type
   USE qs_scf_types,                    ONLY: qs_ot_restart_adaptive_nint,&
                                              qs_ot_restart_adaptive_nlog,&
                                              qs_ot_restart_adaptive_nreal,&
                                              qs_ot_restart_dense_matrix_type,&
                                              qs_scf_env_type,&
                                              qs_scf_ot_restart_clear
   USE scf_control_types,               ONLY: smear_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'kpoint_io'

   INTEGER, PARAMETER                   :: version = 1
   INTEGER, PARAMETER                   :: ot_restart_version = 39
   REAL(KIND=dp), PARAMETER             :: ot_restart_geometry_tolerance = 1.0E-10_dp
   REAL(KIND=dp), PARAMETER             :: ot_restart_objective_tolerance = 1.0E-10_dp
   CHARACTER(LEN=*), PARAMETER          :: ot_restart_marker = "CP2K_KPOINT_OT_RESTART"

   PUBLIC :: read_kpoints_restart, &
             write_kpoints_restart
   PUBLIC :: kpoint_ot_restart_geometry_compatible, &
             kpoint_ot_restart_smearing_compatible
   PUBLIC :: get_cell, write_kpoints_file_header

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param denmat ...
!> \param kpoints ...
!> \param scf_env ...
!> \param dft_section ...
!> \param particle_set ...
!> \param qs_kind_set ...
!> \param cell current simulation cell
!> \param smear current smearing objective
! **************************************************************************************************
   SUBROUTINE write_kpoints_restart(denmat, kpoints, scf_env, dft_section, &
                                    particle_set, qs_kind_set, cell, smear)

      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: denmat
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(section_vals_type), POINTER                   :: dft_section
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(smear_type), INTENT(IN)                       :: smear

      CHARACTER(LEN=*), PARAMETER :: routineN = 'write_kpoints_restart'
      CHARACTER(LEN=30), DIMENSION(2), PARAMETER :: &
         keys = ["SCF%PRINT%RESTART_HISTORY", "SCF%PRINT%RESTART        "]

      INTEGER                                            :: handle, ic, ikey, ires, ispin, nimages, &
                                                            nkp, nspin
      INTEGER, DIMENSION(2)                              :: kp_range
      INTEGER, DIMENSION(3)                              :: image_cell
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: ot_restricted, use_real_wfn
      REAL(KIND=dp), DIMENSION(:), POINTER               :: wkp
      TYPE(cp_fm_type), POINTER                          :: fmat
      TYPE(cp_logger_type), POINTER                      :: logger

      CALL timeset(routineN, handle)
      logger => cp_get_default_logger()
      NULLIFY (wkp)

      IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                                           dft_section, keys(1)), cp_p_file) .OR. &
          BTEST(cp_print_key_should_output(logger%iter_info, &
                                           dft_section, keys(2)), cp_p_file)) THEN

         nspin = SIZE(denmat, 1)
         nimages = SIZE(denmat, 2)
         IF (ASSOCIATED(scf_env%qs_ot_env)) THEN
            CALL get_kpoint_info(kpoint=kpoints, nkp=nkp, kp_range=kp_range, wkp=wkp, &
                                 use_real_wfn=use_real_wfn)
            ot_restricted = ANY(scf_env%qs_ot_env(:)%restricted)
            CALL qs_ot_check_channel_context(scf_env%qs_ot_env, nspin, &
                                             nkpoint=nkp, &
                                             restricted=ot_restricted, &
                                             require_kpoint=.TRUE., &
                                             kp_range=kp_range, &
                                             wkp=wkp, &
                                             require_local_state=.TRUE., &
                                             require_complex_state=.NOT. use_real_wfn)
         END IF

         DO ikey = 1, SIZE(keys)

            IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                                                 dft_section, keys(ikey)), cp_p_file)) THEN

               ires = cp_print_key_unit_nr(logger, dft_section, keys(ikey), &
                                           extension=".kp", file_status="REPLACE", file_action="WRITE", &
                                           do_backup=.TRUE., file_form="UNFORMATTED")

               CALL write_kpoints_file_header(qs_kind_set, particle_set, ires)

               NULLIFY (cell_to_index)
               IF (nimages > 1) THEN
                  CALL get_kpoint_info(kpoint=kpoints, cell_to_index=cell_to_index)
               END IF
               CPASSERT(ASSOCIATED(scf_env%scf_work1))
               fmat => scf_env%scf_work1(1)

               DO ispin = 1, nspin
                  IF (ires > 0) WRITE (ires) ispin, nspin, nimages
                  DO ic = 1, nimages
                     IF (nimages > 1) THEN
                        image_cell = get_cell(ic, cell_to_index)
                     ELSE
                        image_cell = 0
                     END IF
                     IF (ires > 0) WRITE (ires) ic, image_cell
                     CALL copy_dbcsr_to_fm(denmat(ispin, ic)%matrix, fmat)
                     CALL cp_fm_write_unformatted(fmat, ires)
                  END DO
               END DO
               IF (ASSOCIATED(scf_env%qs_ot_env)) THEN
                  CALL write_kpoint_ot_restart_payload(scf_env%qs_ot_env, fmat, ires, &
                                                       cell, particle_set, smear)
               END IF

               CALL cp_print_key_finished_output(ires, logger, dft_section, TRIM(keys(ikey)))

            END IF

         END DO

      END IF

      CALL timestop(handle)

   END SUBROUTINE write_kpoints_restart

! **************************************************************************************************
!> \brief test whether a serialized OT minimizer state belongs to the current geometry
!> \param current_hmat current cell matrix
!> \param current_positions current Cartesian atom positions
!> \param restart_hmat serialized cell matrix
!> \param restart_positions serialized Cartesian atom positions
!> \return true when the minimizer coordinates and history can be continued
! **************************************************************************************************
   PURE FUNCTION kpoint_ot_restart_geometry_compatible(current_hmat, current_positions, &
                                                       restart_hmat, restart_positions) RESULT(compatible)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: current_hmat, current_positions, &
                                                            restart_hmat, restart_positions
      LOGICAL                                            :: compatible

      REAL(KIND=dp)                                      :: geometry_delta, geometry_scale

      compatible = ALL(SHAPE(current_hmat) == SHAPE(restart_hmat)) .AND. &
                   ALL(SHAPE(current_positions) == SHAPE(restart_positions))
      IF (.NOT. compatible) RETURN

      geometry_scale = MAX(1.0_dp, MAXVAL(ABS(current_hmat)), MAXVAL(ABS(restart_hmat)), &
                           MAXVAL(ABS(current_positions)), MAXVAL(ABS(restart_positions)))
      geometry_delta = MAX(MAXVAL(ABS(current_hmat - restart_hmat)), &
                           MAXVAL(ABS(current_positions - restart_positions)))
      compatible = geometry_delta <= ot_restart_geometry_tolerance*geometry_scale

   END FUNCTION kpoint_ot_restart_geometry_compatible

! **************************************************************************************************
!> \brief test whether serialized OT history belongs to the current smearing objective
!> \param current current smearing settings
!> \param restart_do_smear serialized smearing activation
!> \param restart_common_mu serialized common-chemical-potential setting
!> \param restart_method serialized smearing method
!> \param restart_temperature serialized electronic temperature
!> \param restart_width serialized smearing width
!> \param restart_fixed_mag_mom serialized fixed magnetic moment
!> \return true when the Mermin objective and fixed-N grouping are unchanged
! **************************************************************************************************
   PURE FUNCTION kpoint_ot_restart_smearing_compatible( &
      current, restart_do_smear, restart_common_mu, restart_method, restart_temperature, &
      restart_width, restart_fixed_mag_mom) RESULT(compatible)
      TYPE(smear_type), INTENT(IN)                       :: current
      LOGICAL, INTENT(IN)                                :: restart_do_smear, restart_common_mu
      INTEGER, INTENT(IN)                                :: restart_method
      REAL(KIND=dp), INTENT(IN)                          :: restart_temperature, restart_width, &
                                                            restart_fixed_mag_mom
      LOGICAL                                            :: compatible

      REAL(KIND=dp)                                      :: scale

      compatible = current%do_smear .EQV. restart_do_smear
      IF (.NOT. compatible .OR. .NOT. current%do_smear) RETURN

      compatible = current%common_mu .EQV. restart_common_mu
      compatible = compatible .AND. current%method == restart_method
      scale = MAX(1.0_dp, ABS(current%electronic_temperature), ABS(restart_temperature))
      compatible = compatible .AND. &
                   ABS(current%electronic_temperature - restart_temperature) <= &
                   ot_restart_objective_tolerance*scale
      scale = MAX(1.0_dp, ABS(current%smearing_width), ABS(restart_width))
      compatible = compatible .AND. &
                   ABS(current%smearing_width - restart_width) <= &
                   ot_restart_objective_tolerance*scale
      scale = MAX(1.0_dp, ABS(current%fixed_mag_mom), ABS(restart_fixed_mag_mom))
      compatible = compatible .AND. &
                   ABS(current%fixed_mag_mom - restart_fixed_mag_mom) <= &
                   ot_restart_objective_tolerance*scale

   END FUNCTION kpoint_ot_restart_smearing_compatible

! **************************************************************************************************
!> \brief ...
!> \param ic ...
!> \param cell_to_index ...
!> \return ...
! **************************************************************************************************
   FUNCTION get_cell(ic, cell_to_index) RESULT(cell)
      INTEGER, INTENT(IN)                                :: ic
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      INTEGER, DIMENSION(3)                              :: cell

      INTEGER                                            :: i1, i2, i3

      cell = 0
      ALLCELL: DO i3 = LBOUND(cell_to_index, 3), UBOUND(cell_to_index, 3)
         DO i2 = LBOUND(cell_to_index, 2), UBOUND(cell_to_index, 2)
            DO i1 = LBOUND(cell_to_index, 1), UBOUND(cell_to_index, 1)
               IF (cell_to_index(i1, i2, i3) == ic) THEN
                  cell(1) = i1
                  cell(2) = i2
                  cell(3) = i3
                  EXIT ALLCELL
               END IF
            END DO
         END DO
      END DO ALLCELL

   END FUNCTION get_cell

! **************************************************************************************************
!> \brief ...
!> \param qs_kind_set ...
!> \param particle_set ...
!> \param ires ...
!> \param basis_type ...
! **************************************************************************************************
   SUBROUTINE write_kpoints_file_header(qs_kind_set, particle_set, ires, basis_type)

      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      INTEGER                                            :: ires
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL             :: basis_type

      CHARACTER(LEN=*), PARAMETER :: routineN = 'write_kpoints_file_header'

      CHARACTER(LEN=default_string_length)               :: basis_type_local
      INTEGER                                            :: handle, iatom, ikind, iset, ishell, &
                                                            lmax, lshell, nao, natom, nset, &
                                                            nset_max, nsgf, nshell_max
      INTEGER, DIMENSION(:), POINTER                     :: nset_info, nshell
      INTEGER, DIMENSION(:, :), POINTER                  :: l, nshell_info
      INTEGER, DIMENSION(:, :, :), POINTER               :: nso_info
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set
      TYPE(qs_dftb_atom_type), POINTER                   :: dftb_parameter

      CALL timeset(routineN, handle)

      IF (ires > 0) THEN

         IF (PRESENT(basis_type)) THEN
            basis_type_local = basis_type
         ELSE
            basis_type_local = "ORB"
         END IF

         ! create some info about the basis set first
         natom = SIZE(particle_set, 1)
         nset_max = 0
         nshell_max = 0

         DO iatom = 1, natom
            NULLIFY (orb_basis_set, dftb_parameter)
            CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
            CALL get_qs_kind(qs_kind_set(ikind), basis_type=basis_type_local, &
                             basis_set=orb_basis_set, dftb_parameter=dftb_parameter)
            IF (ASSOCIATED(orb_basis_set)) THEN
               CALL get_gto_basis_set(gto_basis_set=orb_basis_set, &
                                      nset=nset, &
                                      nshell=nshell, &
                                      l=l)
               nset_max = MAX(nset_max, nset)
               DO iset = 1, nset
                  nshell_max = MAX(nshell_max, nshell(iset))
               END DO
            ELSE IF (ASSOCIATED(dftb_parameter)) THEN
               CALL get_dftb_atom_param(dftb_parameter, lmax=lmax)
               nset_max = MAX(nset_max, 1)
               nshell_max = MAX(nshell_max, lmax + 1)
            ELSE
               CPABORT("Unknown basis type.")
            END IF
         END DO

         ALLOCATE (nso_info(nshell_max, nset_max, natom))
         nso_info(:, :, :) = 0

         ALLOCATE (nshell_info(nset_max, natom))
         nshell_info(:, :) = 0

         ALLOCATE (nset_info(natom))
         nset_info(:) = 0

         nao = 0
         DO iatom = 1, natom
            NULLIFY (orb_basis_set, dftb_parameter)
            CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
            CALL get_qs_kind(qs_kind_set(ikind), basis_type=basis_type_local, &
                             basis_set=orb_basis_set, dftb_parameter=dftb_parameter)
            IF (ASSOCIATED(orb_basis_set)) THEN
               CALL get_gto_basis_set(gto_basis_set=orb_basis_set, nsgf=nsgf, nset=nset, nshell=nshell, l=l)
               nao = nao + nsgf
               nset_info(iatom) = nset
               DO iset = 1, nset
                  nshell_info(iset, iatom) = nshell(iset)
                  DO ishell = 1, nshell(iset)
                     lshell = l(ishell, iset)
                     nso_info(ishell, iset, iatom) = nso(lshell)
                  END DO
               END DO
            ELSE IF (ASSOCIATED(dftb_parameter)) THEN
               CALL get_dftb_atom_param(dftb_parameter, lmax=lmax)
               nset_info(iatom) = 1
               nshell_info(1, iatom) = lmax + 1
               DO ishell = 1, lmax + 1
                  lshell = ishell - 1
                  nso_info(ishell, 1, iatom) = nso(lshell)
               END DO
               nao = nao + (lmax + 1)**2
            ELSE
               CPABORT("Unknown basis set type.")
            END IF
         END DO

         WRITE (ires) version
         WRITE (ires) natom, nao, nset_max, nshell_max
         WRITE (ires) nset_info
         WRITE (ires) nshell_info
         WRITE (ires) nso_info

         DEALLOCATE (nset_info, nshell_info, nso_info)
      END IF

      CALL timestop(handle)

   END SUBROUTINE write_kpoints_file_header

! **************************************************************************************************
!> \brief ...
!> \param denmat ...
!> \param kpoints ...
!> \param fmwork ...
!> \param natom ...
!> \param para_env ...
!> \param id_nr ...
!> \param dft_section ...
!> \param natom_mismatch ...
!> \param scf_env ...
!> \param cell current simulation cell when an OT payload may be restored
!> \param particle_set current atomic geometry when an OT payload may be restored
!> \param smear current smearing objective
! **************************************************************************************************
   SUBROUTINE read_kpoints_restart(denmat, kpoints, fmwork, natom, &
                                   para_env, id_nr, dft_section, natom_mismatch, scf_env, &
                                   cell, particle_set, smear)

      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: denmat
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(cp_fm_type), DIMENSION(:), INTENT(INOUT)      :: fmwork
      INTEGER, INTENT(IN)                                :: natom
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(IN)                                :: id_nr
      TYPE(section_vals_type), POINTER                   :: dft_section
      LOGICAL, INTENT(OUT)                               :: natom_mismatch
      TYPE(qs_scf_env_type), OPTIONAL, POINTER           :: scf_env
      TYPE(cell_type), OPTIONAL, POINTER                 :: cell
      TYPE(particle_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: particle_set
      TYPE(smear_type), INTENT(IN)                       :: smear

      CHARACTER(LEN=*), PARAMETER :: routineN = 'read_kpoints_restart'

      CHARACTER(LEN=default_path_length)                 :: file_name
      INTEGER                                            :: handle, restart_unit
      LOGICAL                                            :: exist
      TYPE(cp_logger_type), POINTER                      :: logger

      CALL timeset(routineN, handle)
      logger => cp_get_default_logger()

      restart_unit = -1
      IF (PRESENT(scf_env)) CALL qs_scf_ot_restart_clear(scf_env)

      IF (para_env%is_source()) THEN

         CALL wfn_restart_file_name(file_name, exist, dft_section, logger, kp=.TRUE.)
         IF (id_nr /= 0) THEN
            ! Is it one of the backup files?
            file_name = TRIM(file_name)//".bak-"//ADJUSTL(cp_to_string(id_nr))
         END IF

         CALL open_file(file_name=file_name, &
                        file_action="READ", &
                        file_form="UNFORMATTED", &
                        file_status="OLD", &
                        unit_number=restart_unit)

      END IF

      IF (PRESENT(scf_env)) THEN
         CPASSERT(PRESENT(cell))
         CPASSERT(PRESENT(particle_set))
         CPASSERT(ASSOCIATED(cell))
         CPASSERT(ASSOCIATED(particle_set))
         CALL read_kpoints_restart_low(denmat, kpoints, fmwork(1), para_env, &
                                       natom, restart_unit, natom_mismatch, scf_env=scf_env, &
                                       cell=cell, particle_set=particle_set, smear=smear)
      ELSE
         CALL read_kpoints_restart_low(denmat, kpoints, fmwork(1), para_env, &
                                       natom, restart_unit, natom_mismatch)
      END IF

      ! only the io_node returns natom_mismatch, must broadcast it
      CALL para_env%bcast(natom_mismatch)

      ! Close restart file
      IF (para_env%is_source()) CALL close_file(unit_number=restart_unit)

      CALL timestop(handle)

   END SUBROUTINE read_kpoints_restart

! **************************************************************************************************
!> \brief Reading the mos from apreviously defined restart file
!> \param denmat ...
!> \param kpoints ...
!> \param fmat ...
!> \param para_env ...
!> \param natom ...
!> \param rst_unit ...
!> \param natom_mismatch ...
!> \param scf_env ...
!> \param cell ...
!> \param particle_set ...
!> \param smear current smearing objective
!> \par History
!>      12.2007 created [MI]
!> \author MI
! **************************************************************************************************
   SUBROUTINE read_kpoints_restart_low(denmat, kpoints, fmat, para_env, &
                                       natom, rst_unit, natom_mismatch, scf_env, cell, particle_set, smear)

      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: denmat
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(cp_fm_type), INTENT(INOUT)                    :: fmat
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(IN)                                :: natom, rst_unit
      LOGICAL, INTENT(OUT)                               :: natom_mismatch
      TYPE(qs_scf_env_type), OPTIONAL, POINTER           :: scf_env
      TYPE(cell_type), OPTIONAL, POINTER                 :: cell
      TYPE(particle_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: particle_set
      TYPE(smear_type), INTENT(IN), OPTIONAL             :: smear

      INTEGER                                            :: ic, image, ispin, nao, nao_read, &
                                                            natom_read, ni, nimages, nset_max, &
                                                            nshell_max, nspin, nspin_read, &
                                                            version_read
      INTEGER, DIMENSION(3)                              :: image_cell
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: natom_match

      CPASSERT(ASSOCIATED(denmat))
      CALL dbcsr_get_info(denmat(1, 1)%matrix, nfullrows_total=nao)

      nspin = SIZE(denmat, 1)
      nimages = SIZE(denmat, 2)
      NULLIFY (cell_to_index)
      IF (nimages > 1) THEN
         CALL get_kpoint_info(kpoint=kpoints, cell_to_index=cell_to_index)
      END IF

      IF (para_env%is_source()) THEN
         READ (rst_unit) version_read
         IF (version_read /= version) THEN
            CALL cp_abort(__LOCATION__, &
                          " READ RESTART : Version of restart file not supported")
         END IF
         READ (rst_unit) natom_read, nao_read, nset_max, nshell_max
         natom_match = (natom_read == natom)
         IF (natom_match) THEN
            IF (nao_read /= nao) THEN
               CALL cp_abort(__LOCATION__, &
                             " READ RESTART : Different number of basis functions")
            END IF
            READ (rst_unit)
            READ (rst_unit)
            READ (rst_unit)
         END IF
      END IF

      ! make natom_match and natom_mismatch uniform across all nodes
      CALL para_env%bcast(natom_match)
      natom_mismatch = .NOT. natom_match
      ! handle natom_match false
      IF (natom_mismatch) THEN
         CALL cp_warn(__LOCATION__, &
                      " READ RESTART : WARNING : DIFFERENT natom, returning ")
      ELSE

         DO
            ispin = 0
            IF (para_env%is_source()) THEN
               READ (rst_unit) ispin, nspin_read, ni
            END IF
            CALL para_env%bcast(ispin)
            CALL para_env%bcast(nspin_read)
            CALL para_env%bcast(ni)
            IF (nspin_read /= nspin) THEN
               CALL cp_abort(__LOCATION__, &
                             " READ RESTART : Different spin polarisation not supported")
            END IF
            DO ic = 1, ni
               IF (para_env%is_source()) THEN
                  READ (rst_unit) image, image_cell
               END IF
               CALL para_env%bcast(image)
               CALL para_env%bcast(image_cell)
               !
               CALL cp_fm_read_unformatted(fmat, rst_unit)
               !
               IF (nimages > 1) THEN
                  image = get_index(image_cell, cell_to_index)
               ELSE IF (SUM(ABS(image_cell(:))) == 0) THEN
                  image = 1
               ELSE
                  image = 0
               END IF
               IF (image >= 1 .AND. image <= nimages) THEN
                  CALL copy_fm_to_dbcsr(fmat, denmat(ispin, image)%matrix)
               END IF
            END DO
            IF (ispin == nspin_read) EXIT
         END DO
         IF (PRESENT(scf_env)) THEN
            CPASSERT(PRESENT(cell))
            CPASSERT(PRESENT(particle_set))
            CPASSERT(PRESENT(smear))
            CALL read_kpoint_ot_restart_payload(rst_unit, para_env, fmat, scf_env, &
                                                cell, particle_set, smear)
         END IF

      END IF

   END SUBROUTINE read_kpoints_restart_low

! **************************************************************************************************
!> \brief write optional k-point OT restart state after the density payload
!> \param qs_ot_env OT channels
!> \param fmat_template global restart FM distribution
!> \param ires restart unit
!> \param cell current simulation cell
!> \param particle_set current atomic geometry
!> \param smear smearing objective associated with the minimizer history
! **************************************************************************************************
   SUBROUTINE write_kpoint_ot_restart_payload(qs_ot_env, fmat_template, ires, cell, particle_set, smear)
      TYPE(qs_ot_type), DIMENSION(:), INTENT(IN)         :: qs_ot_env
      TYPE(cp_fm_type), INTENT(IN)                       :: fmat_template
      INTEGER, INTENT(IN)                                :: ires
      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(smear_type), INTENT(IN)                       :: smear

      CHARACTER(LEN=20)                                  :: requested_ot_method_payload
      CHARACTER(LEN=20), &
         DIMENSION(qs_ot_restart_signature_nchar)        :: settings_char_payload
      INTEGER                                            :: iatom, ichannel, ihist, nchannel, ndiis, &
                                                            nener, nmo, nmo_ortho, source_count
      INTEGER, DIMENSION(16)                             :: int_payload
      INTEGER, DIMENSION(19)                             :: schur_int_payload
      INTEGER, DIMENSION(3)                              :: smear_int_payload
      INTEGER, DIMENSION(32)                             :: log_payload_int
      INTEGER, DIMENSION(qs_ot_restart_adaptive_nint)    :: adaptive_int_payload
      INTEGER, DIMENSION(qs_ot_restart_adaptive_nlog)    :: adaptive_log_payload_int
      INTEGER, DIMENSION(qs_ot_restart_signature_nint)   :: settings_int_payload
      LOGICAL :: channel_source, has_diis_scalars, has_ener_dx, has_ener_gx_old, has_ener_history, &
         has_ener_x, has_lbfgs_scalars, has_matrix_dx, has_matrix_dx_im, has_matrix_gx_old, &
         has_matrix_gx_old_im, has_matrix_history, has_matrix_history_im, has_matrix_x, &
         has_matrix_x_im, has_rot_history, has_rot_history_im, has_rot_mat_dx, has_rot_mat_dx_im, &
         has_rot_mat_gx_old, has_rot_mat_gx_old_im, has_rot_mat_x, has_rot_mat_x_im, &
         settings_compatible
      LOGICAL, DIMENSION(32)                             :: log_payload
      LOGICAL, DIMENSION(qs_ot_restart_adaptive_nlog)    :: adaptive_log_payload
      LOGICAL, DIMENSION(qs_ot_restart_signature_nlog)   :: settings_log_payload
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: restart_positions
      REAL(KIND=dp), DIMENSION(139)                      :: schur_real_payload
      REAL(KIND=dp), DIMENSION(25)                       :: real_payload
      REAL(KIND=dp), DIMENSION(3)                        :: smear_real_payload
      REAL(KIND=dp), DIMENSION(53)                       :: ot_energy_payload, ot_grad_payload, &
                                                            ot_pos_payload
      REAL(KIND=dp), &
         DIMENSION(qs_ot_restart_adaptive_nreal)         :: adaptive_real_payload
      REAL(KIND=dp), &
         DIMENSION(qs_ot_restart_signature_nreal)        :: settings_real_payload
      TYPE(dbcsr_type), POINTER                          :: restart_matrix
      TYPE(mp_para_env_type), POINTER                    :: global_para_env

      NULLIFY (global_para_env, restart_matrix)
      global_para_env => fmat_template%matrix_struct%para_env
      nchannel = SIZE(qs_ot_env)
      CPASSERT(ASSOCIATED(cell))
      CPASSERT(ASSOCIATED(particle_set))
      ALLOCATE (restart_positions(3, SIZE(particle_set)))
      DO iatom = 1, SIZE(particle_set)
         restart_positions(:, iatom) = particle_set(iatom)%r
      END DO
      smear_int_payload = [MERGE(1, 0, smear%do_smear), MERGE(1, 0, smear%common_mu), smear%method]
      smear_real_payload = [smear%electronic_temperature, smear%smearing_width, smear%fixed_mag_mom]
      CPASSERT(nchannel > 0)
      CALL qs_ot_settings_restart_signature(qs_ot_env(1)%settings, settings_char_payload, &
                                            settings_int_payload, settings_log_payload, &
                                            settings_real_payload)
      requested_ot_method_payload = qs_ot_env(1)%settings%requested_ot_method
      DO ichannel = 2, nchannel
         settings_compatible = qs_ot_settings_restart_compatible( &
                               qs_ot_env(ichannel)%settings, settings_char_payload, &
                               settings_int_payload, settings_log_payload, settings_real_payload)
         settings_compatible = settings_compatible .AND. &
                               qs_ot_env(ichannel)%settings%requested_ot_method == &
                               requested_ot_method_payload
         CPASSERT(settings_compatible)
      END DO
      IF (ires > 0) THEN
         WRITE (ires) ot_restart_marker
         WRITE (ires) ot_restart_version, nchannel
         WRITE (ires) SIZE(particle_set)
         WRITE (ires) cell%hmat
         WRITE (ires) restart_positions
         WRITE (ires) smear_int_payload
         WRITE (ires) smear_real_payload
         WRITE (ires) settings_char_payload
         WRITE (ires) settings_int_payload
         WRITE (ires) settings_log_payload
         WRITE (ires) settings_real_payload
         WRITE (ires) requested_ot_method_payload
      END IF
      DEALLOCATE (restart_positions)
      DO ichannel = 1, nchannel
         nmo = 0
         nmo_ortho = 0
         IF (ASSOCIATED(qs_ot_env(ichannel)%matrix_c0)) THEN
            CALL dbcsr_get_info(qs_ot_env(ichannel)%matrix_c0, nfullcols_total=nmo)
         END IF
         IF (ASSOCIATED(qs_ot_env(ichannel)%matrix_sc0)) THEN
            CALL dbcsr_get_info(qs_ot_env(ichannel)%matrix_sc0, nfullcols_total=nmo_ortho)
         END IF
         has_rot_mat_x = qs_ot_env(ichannel)%state_allocated .AND. &
                         ASSOCIATED(qs_ot_env(ichannel)%rot_mat_x)
         has_rot_mat_x_im = qs_ot_env(ichannel)%state_allocated .AND. &
                            ASSOCIATED(qs_ot_env(ichannel)%rot_mat_x_im)
         has_matrix_x = qs_ot_env(ichannel)%state_allocated .AND. ASSOCIATED(qs_ot_env(ichannel)%matrix_x)
         has_matrix_x_im = qs_ot_env(ichannel)%state_allocated .AND. &
                           ASSOCIATED(qs_ot_env(ichannel)%matrix_x_im)
         has_ener_x = qs_ot_env(ichannel)%state_allocated .AND. ASSOCIATED(qs_ot_env(ichannel)%ener_x)
         has_matrix_dx = qs_ot_env(ichannel)%state_allocated .AND. ASSOCIATED(qs_ot_env(ichannel)%matrix_dx)
         has_matrix_dx_im = qs_ot_env(ichannel)%state_allocated .AND. &
                            ASSOCIATED(qs_ot_env(ichannel)%matrix_dx_im)
         has_matrix_gx_old = qs_ot_env(ichannel)%state_allocated .AND. &
                             ASSOCIATED(qs_ot_env(ichannel)%matrix_gx_old)
         has_matrix_gx_old_im = qs_ot_env(ichannel)%state_allocated .AND. &
                                ASSOCIATED(qs_ot_env(ichannel)%matrix_gx_old_im)
         has_rot_mat_dx = qs_ot_env(ichannel)%state_allocated .AND. ASSOCIATED(qs_ot_env(ichannel)%rot_mat_dx)
         has_rot_mat_dx_im = qs_ot_env(ichannel)%state_allocated .AND. &
                             ASSOCIATED(qs_ot_env(ichannel)%rot_mat_dx_im)
         has_rot_mat_gx_old = qs_ot_env(ichannel)%state_allocated .AND. &
                              ASSOCIATED(qs_ot_env(ichannel)%rot_mat_gx_old)
         has_rot_mat_gx_old_im = qs_ot_env(ichannel)%state_allocated .AND. &
                                 ASSOCIATED(qs_ot_env(ichannel)%rot_mat_gx_old_im)
         has_ener_dx = qs_ot_env(ichannel)%state_allocated .AND. ASSOCIATED(qs_ot_env(ichannel)%ener_dx)
         has_ener_gx_old = qs_ot_env(ichannel)%state_allocated .AND. &
                           ASSOCIATED(qs_ot_env(ichannel)%ener_gx_old)
         has_diis_scalars = qs_ot_env(ichannel)%state_allocated .AND. ASSOCIATED(qs_ot_env(ichannel)%ls_diis)
         has_matrix_history = qs_ot_env(ichannel)%state_allocated .AND. &
                              ASSOCIATED(qs_ot_env(ichannel)%matrix_h_e) .AND. &
                              ASSOCIATED(qs_ot_env(ichannel)%matrix_h_x)
         has_matrix_history_im = qs_ot_env(ichannel)%state_allocated .AND. &
                                 ASSOCIATED(qs_ot_env(ichannel)%matrix_h_e_im) .AND. &
                                 ASSOCIATED(qs_ot_env(ichannel)%matrix_h_x_im)
         has_rot_history = qs_ot_env(ichannel)%state_allocated .AND. &
                           ASSOCIATED(qs_ot_env(ichannel)%rot_mat_h_e) .AND. &
                           ASSOCIATED(qs_ot_env(ichannel)%rot_mat_h_x)
         has_rot_history_im = qs_ot_env(ichannel)%state_allocated .AND. &
                              ASSOCIATED(qs_ot_env(ichannel)%rot_mat_h_e_im) .AND. &
                              ASSOCIATED(qs_ot_env(ichannel)%rot_mat_h_x_im)
         has_ener_history = qs_ot_env(ichannel)%state_allocated .AND. &
                            ASSOCIATED(qs_ot_env(ichannel)%ener_h_e) .AND. &
                            ASSOCIATED(qs_ot_env(ichannel)%ener_h_x)
         has_lbfgs_scalars = qs_ot_env(ichannel)%state_allocated .AND. &
                             ASSOCIATED(qs_ot_env(ichannel)%lbfgs_rho) .AND. &
                             ASSOCIATED(qs_ot_env(ichannel)%lbfgs_yy)
         nener = 0
         IF (has_ener_x) nener = SIZE(qs_ot_env(ichannel)%ener_x)
         ndiis = 0
         IF (has_diis_scalars) ndiis = SIZE(qs_ot_env(ichannel)%energy_h)
         IF (has_lbfgs_scalars) THEN
            ndiis = SIZE(qs_ot_env(ichannel)%lbfgs_rho) + 1
            CPASSERT(ndiis > 1)
            CPASSERT(ndiis == SIZE(qs_ot_env(ichannel)%lbfgs_yy) + 1)
         END IF

         int_payload = [ichannel, qs_ot_env(ichannel)%spin_index, &
                        qs_ot_env(ichannel)%kpoint_index, &
                        qs_ot_env(ichannel)%local_kpoint_index, nmo, nmo_ortho, &
                        qs_ot_env(ichannel)%line_search_left, &
                        qs_ot_env(ichannel)%line_search_right, &
                        qs_ot_env(ichannel)%line_search_mid, &
                        qs_ot_env(ichannel)%line_search_count, &
                        qs_ot_env(ichannel)%diis_iter, &
                        qs_ot_env(ichannel)%mermin_ls_best_sample_index, &
                        qs_ot_env(ichannel)%mermin_ls_accepted_quality_samples, &
                        nener, ndiis, qs_ot_env(ichannel)%kpoint_subspace_accepted_samples]
         log_payload = [qs_ot_env(ichannel)%state_allocated, &
                        qs_ot_env(ichannel)%has_complex_kpoint_state, &
                        qs_ot_env(ichannel)%settings%do_rotation, &
                        qs_ot_env(ichannel)%settings%do_ener, &
                        has_rot_mat_x, has_rot_mat_x_im, has_ener_x, &
                        qs_ot_env(ichannel)%line_search_might_be_done, &
                        qs_ot_env(ichannel)%energy_only, has_matrix_dx, &
                        has_matrix_dx_im, has_matrix_gx_old, has_matrix_gx_old_im, &
                        has_rot_mat_dx, has_rot_mat_dx_im, has_rot_mat_gx_old, &
                        has_rot_mat_gx_old_im, has_ener_dx, has_ener_gx_old, &
                        has_diis_scalars, has_matrix_history, has_matrix_history_im, &
                        has_rot_history, has_rot_history_im, has_ener_history, &
                        has_matrix_x, has_matrix_x_im, &
                        qs_ot_env(ichannel)%mermin_ls_sample_pending, &
                        qs_ot_env(ichannel)%ener_direction_linear, &
                        qs_ot_env(ichannel)%diis_last_candidate_sd, &
                        qs_ot_env(ichannel)%diis_last_candidate_response, &
                        has_lbfgs_scalars]
         real_payload = [qs_ot_env(ichannel)%kpoint_weight, &
                         qs_ot_env(ichannel)%delta, qs_ot_env(ichannel)%gnorm, &
                         qs_ot_env(ichannel)%gnorm_old, qs_ot_env(ichannel)%etotal, &
                         qs_ot_env(ichannel)%gradient, qs_ot_env(ichannel)%ds_min, &
                         qs_ot_env(ichannel)%broyden_adaptive_sigma, &
                         qs_ot_env(ichannel)%mermin_ls_reference_energy, &
                         qs_ot_env(ichannel)%mermin_ls_best_sample_energy, &
                         qs_ot_env(ichannel)%mermin_ls_best_sample_drop, &
                         qs_ot_env(ichannel)%mermin_ls_model_drop, &
                         qs_ot_env(ichannel)%mermin_ls_model_quality, &
                         qs_ot_env(ichannel)%mermin_ls_best_sample_pos, &
                         qs_ot_env(ichannel)%mermin_ls_next_pos, &
                         qs_ot_env(ichannel)%mermin_ls_accepted_quality, &
                         qs_ot_env(ichannel)%mermin_ls_accepted_quality_ema, &
                         qs_ot_env(ichannel)%mermin_ls_accepted_gap, &
                         qs_ot_env(ichannel)%mermin_ls_accepted_gap_ema, &
                         qs_ot_env(ichannel)%mermin_ls_reference_delta, &
                         qs_ot_env(ichannel)%diis_last_response_blend, &
                         qs_ot_env(ichannel)%diis_reference_rotation_gnorm, &
                         qs_ot_env(ichannel)%diis_raw_residual_ratio, &
                         qs_ot_env(ichannel)%diis_response_residual_ratio, &
                         qs_ot_env(ichannel)%diis_last_response_model_work]
         schur_int_payload = [qs_ot_env(ichannel)%schur_projector_density_quality_samples, &
                              qs_ot_env(ichannel)%schur_projector_density_audit_samples, &
                              qs_ot_env(ichannel)%schur_projector_density_descent_quality_samples, &
                              qs_ot_env(ichannel)%schur_projector_density_descent_selected_quality_samples, &
                              qs_ot_env(ichannel)%schur_projector_response_selected_quality_samples, &
                              qs_ot_env(ichannel)%schur_projector_physical_direct_quality_samples, &
                              qs_ot_env(ichannel)%schur_projector_physical_direct_selected_quality_samples, &
                              qs_ot_env(ichannel)%schur_projector_phys_direct_density_sel_quality_samples, &
                              qs_ot_env(ichannel)%schur_projector_mixed_direct_quality_samples, &
                              qs_ot_env(ichannel)%schur_projector_spectral_selected_quality_samples, &
                              qs_ot_env(ichannel)%schur_projector_spectral_selected_failure_samples, &
                              MERGE(1, 0, qs_ot_env(ichannel)%schur_projector_density_descent_support_ready), &
                              qs_ot_env(ichannel)%schur_projector_density_full_selected_quality_samples, &
                              qs_ot_env(ichannel)%schur_projector_physical_quality_samples, &
                              qs_ot_env(ichannel)%schur_projector_density_full_direct_quality_samples, &
                              qs_ot_env(ichannel)%schur_projector_density_fdir_sel_quality_samples, &
                              qs_ot_env(ichannel)%schur_projector_density_full_global_quality_samples, &
                              qs_ot_env(ichannel)%schur_projector_density_full_global_sel_quality_samples, &
                              qs_ot_env(ichannel)%schur_projector_physical_selected_quality_samples]
         schur_real_payload = [qs_ot_env(ichannel)%schur_projector_density_quality, &
                               qs_ot_env(ichannel)%schur_projector_density_quality_ema, &
                               qs_ot_env(ichannel)%schur_projector_density_drop_sum, &
                               qs_ot_env(ichannel)%schur_projector_density_measured_sum, &
                               qs_ot_env(ichannel)%schur_projector_density_lowcurv_ema, &
                               qs_ot_env(ichannel)%schur_projector_density_lowcurv_quality_ema, &
                               qs_ot_env(ichannel)%schur_projector_density_audit_amplification_ema, &
                               qs_ot_env(ichannel)%schur_projector_density_audit_amplification_quality_ema, &
                               qs_ot_env(ichannel)%schur_projector_density_audit_lowcurv_ema, &
                               qs_ot_env(ichannel)%schur_projector_density_audit_lowcurv_quality_ema, &
                               qs_ot_env(ichannel)%schur_projector_density_descent_quality, &
                               qs_ot_env(ichannel)%schur_projector_density_descent_quality_ema, &
                               qs_ot_env(ichannel)%schur_projector_density_descent_drop_sum, &
                               qs_ot_env(ichannel)%schur_projector_density_descent_measured_sum, &
                               qs_ot_env(ichannel)%schur_projector_density_descent_coeff_ema, &
                               qs_ot_env(ichannel)%schur_projector_density_descent_support, &
                               qs_ot_env(ichannel)%schur_projector_density_descent_support_ema, &
                               qs_ot_env(ichannel)%schur_projector_density_descent_selected, &
                               qs_ot_env(ichannel)%schur_projector_density_descent_selected_quality, &
                               qs_ot_env(ichannel)%schur_projector_density_descent_selected_quality_ema, &
                               qs_ot_env(ichannel)%schur_projector_density_descent_selected_residual_quality, &
                               qs_ot_env(ichannel)%schur_projector_density_descent_selected_residual_quality_ema, &
                               qs_ot_env(ichannel)%schur_projector_density_descent_selected_drop_sum, &
                               qs_ot_env(ichannel)%schur_projector_density_descent_selected_measured_sum, &
                               qs_ot_env(ichannel)%schur_projector_physical_direct_quality, &
                               qs_ot_env(ichannel)%schur_projector_physical_direct_quality_shadow_ema, &
                               qs_ot_env(ichannel)%schur_projector_physical_direct_drop_sum, &
                               qs_ot_env(ichannel)%schur_projector_physical_direct_measured_sum, &
                               qs_ot_env(ichannel)%schur_projector_physical_direct_ema, &
                               qs_ot_env(ichannel)%schur_projector_physical_direct_selected, &
                               qs_ot_env(ichannel)%schur_projector_physical_direct_selected_quality, &
                               qs_ot_env(ichannel)%schur_projector_physical_direct_selected_quality_ema, &
                               qs_ot_env(ichannel)%schur_projector_physical_direct_selected_residual_quality, &
                               qs_ot_env(ichannel)%schur_projector_physical_direct_selected_residual_quality_ema, &
                               qs_ot_env(ichannel)%schur_projector_physical_direct_selected_drop_sum, &
                               qs_ot_env(ichannel)%schur_projector_physical_direct_selected_measured_sum, &
                               qs_ot_env(ichannel)%schur_projector_phys_direct_density_sel, &
                               qs_ot_env(ichannel)%schur_projector_phys_direct_density_sel_quality, &
                               qs_ot_env(ichannel)%schur_projector_phys_direct_density_sel_quality_ema, &
                               qs_ot_env(ichannel)%schur_projector_phys_direct_density_sel_drop_sum, &
                               qs_ot_env(ichannel)%schur_projector_phys_direct_density_sel_measured_sum, &
                               qs_ot_env(ichannel)%schur_projector_mixed_direct_quality, &
                               qs_ot_env(ichannel)%schur_projector_mixed_direct_quality_ema, &
                               qs_ot_env(ichannel)%schur_projector_mixed_direct_drop_sum, &
                               qs_ot_env(ichannel)%schur_projector_mixed_direct_measured_sum, &
                               qs_ot_env(ichannel)%schur_projector_mixed_direct_weight, &
                               qs_ot_env(ichannel)%schur_projector_mixed_direct_weight_ema, &
                               qs_ot_env(ichannel)%schur_projector_spectral_selected_quality, &
                               qs_ot_env(ichannel)%schur_projector_spectral_selected_quality_ema, &
                               qs_ot_env(ichannel)%schur_projector_spectral_selected_drop_sum, &
                               qs_ot_env(ichannel)%schur_projector_spectral_selected_measured_sum, &
                               qs_ot_env(ichannel)%schur_projector_response_selected, &
                               qs_ot_env(ichannel)%schur_projector_response_selected_quality, &
                               qs_ot_env(ichannel)%schur_projector_response_selected_quality_ema, &
                               qs_ot_env(ichannel)%schur_projector_response_selected_residual_quality, &
                               qs_ot_env(ichannel)%schur_projector_response_selected_residual_quality_ema, &
                               qs_ot_env(ichannel)%schur_projector_response_selected_drop_sum, &
                               qs_ot_env(ichannel)%schur_projector_response_selected_measured_sum, &
                               qs_ot_env(ichannel)%schur_projector_density_full_selected, &
                               qs_ot_env(ichannel)%schur_projector_density_full_selected_quality, &
                               qs_ot_env(ichannel)%schur_projector_density_full_selected_quality_ema, &
                               qs_ot_env(ichannel)%schur_projector_density_full_selected_residual_quality, &
                               qs_ot_env(ichannel)%schur_projector_density_full_selected_residual_quality_ema, &
                               qs_ot_env(ichannel)%schur_projector_density_full_selected_drop_sum, &
                               qs_ot_env(ichannel)%schur_projector_density_full_selected_measured_sum, &
                               qs_ot_env(ichannel)%schur_projector_physical_quality, &
                               qs_ot_env(ichannel)%schur_projector_physical_quality_ema, &
                               qs_ot_env(ichannel)%schur_projector_physical_drop_sum, &
                               qs_ot_env(ichannel)%schur_projector_physical_measured_sum, &
                               qs_ot_env(ichannel)%schur_projector_physical_trace_actual_ema, &
                               qs_ot_env(ichannel)%schur_projector_physical_trace_raw_ema, &
                               qs_ot_env(ichannel)%schur_projector_physical_trace_active_ema, &
                               qs_ot_env(ichannel)%schur_projector_physical_trace_actual_diag_ema, &
                               qs_ot_env(ichannel)%schur_projector_physical_trace_actual_free_ema, &
                               qs_ot_env(ichannel)%schur_projector_physical_trace_raw_diag_ema, &
                               qs_ot_env(ichannel)%schur_projector_physical_trace_raw_free_ema, &
                               qs_ot_env(ichannel)%schur_projector_physical_trace_active_diag_ema, &
                               qs_ot_env(ichannel)%schur_projector_physical_trace_active_free_ema, &
                               qs_ot_env(ichannel)%schur_projector_physical_trace_actual_quality_ema, &
                               qs_ot_env(ichannel)%schur_projector_physical_trace_raw_quality_ema, &
                               qs_ot_env(ichannel)%schur_projector_physical_trace_active_quality_ema, &
                               qs_ot_env(ichannel)%schur_projector_physical_trace_actual_diag_quality_ema, &
                               qs_ot_env(ichannel)%schur_projector_physical_trace_actual_free_quality_ema, &
                               qs_ot_env(ichannel)%schur_projector_physical_trace_raw_diag_quality_ema, &
                               qs_ot_env(ichannel)%schur_projector_physical_trace_raw_free_quality_ema, &
                               qs_ot_env(ichannel)%schur_projector_physical_trace_active_diag_quality_ema, &
                               qs_ot_env(ichannel)%schur_projector_physical_trace_active_free_quality_ema, &
                               qs_ot_env(ichannel)%schur_projector_phys_direct_density_sel_resid_quality, &
                               qs_ot_env(ichannel)%schur_projector_phys_direct_density_sel_resid_quality_ema, &
                               qs_ot_env(ichannel)%schur_projector_density_full_direct_drop, &
                               qs_ot_env(ichannel)%schur_projector_density_full_direct_quality, &
                               qs_ot_env(ichannel)%schur_projector_density_full_direct_quality_ema, &
                               qs_ot_env(ichannel)%schur_projector_density_full_direct_drop_sum, &
                               qs_ot_env(ichannel)%schur_projector_density_full_direct_measured_sum, &
                               qs_ot_env(ichannel)%schur_projector_density_full_direct_trust, &
                               qs_ot_env(ichannel)%schur_projector_density_full_direct_trust_ema, &
                               qs_ot_env(ichannel)%schur_projector_density_full_direct_weighted_trust, &
                               qs_ot_env(ichannel)%schur_projector_density_full_direct_weighted_trust_ema, &
                               qs_ot_env(ichannel)%schur_projector_density_full_direct_abs_trust, &
                               qs_ot_env(ichannel)%schur_projector_density_full_direct_abs_trust_ema, &
                               qs_ot_env(ichannel)%schur_projector_density_fdir_selected, &
                               qs_ot_env(ichannel)%schur_projector_density_fdir_sel_quality, &
                               qs_ot_env(ichannel)%schur_projector_density_fdir_sel_quality_ema, &
                               qs_ot_env(ichannel)%schur_projector_density_fdir_sel_resid_quality, &
                               qs_ot_env(ichannel)%schur_projector_density_fdir_sel_resid_quality_ema, &
                               qs_ot_env(ichannel)%schur_projector_density_fdir_sel_drop_sum, &
                               qs_ot_env(ichannel)%schur_projector_density_fdir_sel_measured_sum, &
                               qs_ot_env(ichannel)%schur_projector_density_full_global_drop, &
                               qs_ot_env(ichannel)%schur_projector_density_full_global_quality, &
                               qs_ot_env(ichannel)%schur_projector_density_full_global_quality_ema, &
                               qs_ot_env(ichannel)%schur_projector_density_full_global_drop_sum, &
                               qs_ot_env(ichannel)%schur_projector_density_full_global_measured_sum, &
                               qs_ot_env(ichannel)%schur_projector_density_full_global_trust, &
                               qs_ot_env(ichannel)%schur_projector_density_full_global_trust_ema, &
                               qs_ot_env(ichannel)%schur_projector_density_full_global_weighted_trust, &
                               qs_ot_env(ichannel)%schur_projector_density_full_global_weighted_trust_ema, &
                               qs_ot_env(ichannel)%schur_projector_density_full_global_abs_trust, &
                               qs_ot_env(ichannel)%schur_projector_density_full_global_abs_trust_ema, &
                               qs_ot_env(ichannel)%schur_projector_physical_trace_invariant_ema, &
                               qs_ot_env(ichannel)%schur_projector_physical_trace_invariant_diag_ema, &
                               qs_ot_env(ichannel)%schur_projector_physical_trace_invariant_free_ema, &
                               qs_ot_env(ichannel)%schur_projector_physical_trace_invariant_quality_ema, &
                               qs_ot_env(ichannel)%schur_projector_physical_trace_invariant_diag_quality_ema, &
                               qs_ot_env(ichannel)%schur_projector_physical_trace_invariant_free_quality_ema, &
                               qs_ot_env(ichannel)%schur_projector_spectral_selected_residual_quality, &
                               qs_ot_env(ichannel)%schur_projector_spectral_selected_residual_quality_ema, &
                               qs_ot_env(ichannel)%schur_projector_density_full_global_selected, &
                               qs_ot_env(ichannel)%schur_projector_density_full_global_sel_quality, &
                               qs_ot_env(ichannel)%schur_projector_density_full_global_sel_quality_ema, &
                               qs_ot_env(ichannel)%schur_projector_density_full_global_sel_resid_quality, &
                               qs_ot_env(ichannel)%schur_projector_density_full_global_sel_resid_quality_ema, &
                               qs_ot_env(ichannel)%schur_projector_density_full_global_sel_drop_sum, &
                               qs_ot_env(ichannel)%schur_projector_density_full_global_sel_measured_sum, &
                               qs_ot_env(ichannel)%schur_projector_physical_selected_quality, &
                               qs_ot_env(ichannel)%schur_projector_physical_selected_quality_line_ema, &
                               qs_ot_env(ichannel)%schur_projector_physical_selected_residual_quality, &
                               qs_ot_env(ichannel)%schur_projector_physical_selected_residual_quality_ema, &
                               qs_ot_env(ichannel)%schur_projector_physical_selected_drop_sum, &
                               qs_ot_env(ichannel)%schur_projector_physical_selected_measured_sum]

         ! Keep the retired split-shadow and occupation-gap slots zero so that
         ! version-20-and-newer restart payload layouts remain unchanged.
         adaptive_int_payload = [qs_ot_env(ichannel)%schur_model_attempt, &
                                 qs_ot_env(ichannel)%schur_cg_attempt, &
                                 qs_ot_env(ichannel)%schur_model_flavor, &
                                 qs_ot_env(ichannel)%schur_cg_flavor, &
                                 qs_ot_env(ichannel)%schur_model_quality_samples, &
                                 qs_ot_env(ichannel)%schur_constrained_quality_samples, &
                                 qs_ot_env(ichannel)%schur_split_quality_samples, &
                                 qs_ot_env(ichannel)%schur_coupled_quality_samples, &
                                 0, 0, &
                                 qs_ot_env(ichannel)%schur_dense_metric_shadow_quality_samples, &
                                 qs_ot_env(ichannel)%schur_step_gauge_quality_samples, &
                                 qs_ot_env(ichannel)%kpoint_physical_delta_groups, &
                                 qs_ot_env(ichannel)%schur_projector_mode_quality_samples, &
                                 qs_ot_env(ichannel)%kpoint_occ_preconditioner_ready_samples]
         adaptive_log_payload = [qs_ot_env(ichannel)%schur_model_pending, &
                                 qs_ot_env(ichannel)%schur_model_quality_ready, &
                                 qs_ot_env(ichannel)%schur_constrained_pending, &
                                 qs_ot_env(ichannel)%schur_constrained_candidate_enabled, &
                                 .FALSE., .FALSE., &
                                 qs_ot_env(ichannel)%schur_dense_metric_shadow_pending, &
                                 qs_ot_env(ichannel)%schur_history_mixed, &
                                 qs_ot_env(ichannel)%schur_split_candidate_enabled, &
                                 qs_ot_env(ichannel)%schur_residual_history_reset_pending, &
                                 qs_ot_env(ichannel)%kpoint_physical_delta_ready, &
                                 qs_ot_env(ichannel)%schur_projector_mode_pending, &
                                 qs_ot_env(ichannel)%kpoint_occ_preconditioner_active]
         adaptive_real_payload = [qs_ot_env(ichannel)%schur_model_drop, &
                                  qs_ot_env(ichannel)%schur_model_calibrated_drop, &
                                  qs_ot_env(ichannel)%schur_model_measured_drop, &
                                  qs_ot_env(ichannel)%schur_model_quality, &
                                  qs_ot_env(ichannel)%schur_model_quality_ema, &
                                  qs_ot_env(ichannel)%schur_model_calibrated_quality, &
                                  qs_ot_env(ichannel)%schur_model_calibrated_quality_ema, &
                                  qs_ot_env(ichannel)%schur_model_gap_ratio, &
                                  qs_ot_env(ichannel)%schur_model_delta, &
                                  qs_ot_env(ichannel)%schur_model_control_delta, &
                                  qs_ot_env(ichannel)%schur_model_residual_drop, &
                                  qs_ot_env(ichannel)%schur_residual_history_quality, &
                                  qs_ot_env(ichannel)%schur_step_gap_norm, &
                                  qs_ot_env(ichannel)%schur_step_gauge_norm, &
                                  qs_ot_env(ichannel)%schur_step_fermi_norm, &
                                  qs_ot_env(ichannel)%schur_step_gauge_ratio, &
                                  qs_ot_env(ichannel)%schur_step_gauge_ratio_ema, &
                                  qs_ot_env(ichannel)%schur_step_gauge_residual_quality, &
                                  qs_ot_env(ichannel)%schur_step_gauge_residual_quality_ema, &
                                  qs_ot_env(ichannel)%kpoint_raw_delta, &
                                  qs_ot_env(ichannel)%kpoint_physical_delta, &
                                  qs_ot_env(ichannel)%kpoint_projected_delta, &
                                  qs_ot_env(ichannel)%kpoint_eps_scf, &
                                  qs_ot_env(ichannel)%kpoint_physical_gnorm, &
                                  qs_ot_env(ichannel)%kpoint_physical_rot_gnorm, &
                                  qs_ot_env(ichannel)%kpoint_physical_energy_gnorm, &
                                  qs_ot_env(ichannel)%kpoint_gauge_gnorm, &
                                  qs_ot_env(ichannel)%schur_constrained_drop, &
                                  qs_ot_env(ichannel)%schur_constrained_calibrated_drop, &
                                  qs_ot_env(ichannel)%schur_constrained_quality, &
                                  qs_ot_env(ichannel)%schur_constrained_quality_ema, &
                                  qs_ot_env(ichannel)%schur_constrained_calibrated_quality, &
                                  qs_ot_env(ichannel)%schur_constrained_calibrated_quality_ema, &
                                  qs_ot_env(ichannel)%schur_constrained_drop_sum, &
                                  qs_ot_env(ichannel)%schur_constrained_measured_sum, &
                                  qs_ot_env(ichannel)%schur_constrained_calibrated_drop_sum, &
                                  qs_ot_env(ichannel)%schur_constrained_gap_ratio, &
                                  qs_ot_env(ichannel)%schur_split_quality, &
                                  qs_ot_env(ichannel)%schur_split_quality_ema, &
                                  qs_ot_env(ichannel)%schur_split_calibrated_quality, &
                                  qs_ot_env(ichannel)%schur_split_calibrated_quality_ema, &
                                  qs_ot_env(ichannel)%schur_split_residual_quality, &
                                  qs_ot_env(ichannel)%schur_split_residual_quality_ema, &
                                  qs_ot_env(ichannel)%schur_split_drop_sum, &
                                  qs_ot_env(ichannel)%schur_split_measured_sum, &
                                  qs_ot_env(ichannel)%schur_split_calibrated_drop_sum, &
                                  0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                  0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                  0.0_dp, 0.0_dp, 0.0_dp, &
                                  qs_ot_env(ichannel)%schur_coupled_quality, &
                                  qs_ot_env(ichannel)%schur_coupled_quality_ema, &
                                  qs_ot_env(ichannel)%schur_coupled_calibrated_quality, &
                                  qs_ot_env(ichannel)%schur_coupled_calibrated_quality_ema, &
                                  qs_ot_env(ichannel)%schur_coupled_residual_quality, &
                                  qs_ot_env(ichannel)%schur_coupled_residual_quality_ema, &
                                  qs_ot_env(ichannel)%schur_coupled_drop_sum, &
                                  qs_ot_env(ichannel)%schur_coupled_measured_sum, &
                                  qs_ot_env(ichannel)%schur_coupled_calibrated_drop_sum, &
                                  qs_ot_env(ichannel)%schur_dense_metric_shadow_drop, &
                                  qs_ot_env(ichannel)%schur_dense_metric_shadow_quality, &
                                  qs_ot_env(ichannel)%schur_dense_metric_shadow_quality_ema, &
                                  qs_ot_env(ichannel)%schur_dense_metric_shadow_drop_sum, &
                                  qs_ot_env(ichannel)%schur_dense_metric_shadow_measured_sum, &
                                  qs_ot_env(ichannel)%schur_projector_mode_drop, &
                                  qs_ot_env(ichannel)%schur_projector_mode_response_drop, &
                                  qs_ot_env(ichannel)%schur_projector_mode_step_negative, &
                                  qs_ot_env(ichannel)%schur_projector_mode_response_negative, &
                                  qs_ot_env(ichannel)%schur_projector_mode_trust, &
                                  qs_ot_env(ichannel)%schur_projector_mode_quality, &
                                  qs_ot_env(ichannel)%schur_projector_mode_quality_ema, &
                                  qs_ot_env(ichannel)%schur_projector_mode_drop_sum, &
                                  qs_ot_env(ichannel)%schur_projector_mode_measured_sum, &
                                  qs_ot_env(ichannel)%schur_projector_mode_response_drop_sum, &
                                  qs_ot_env(ichannel)%schur_projector_mode_step_negative_ema, &
                                  qs_ot_env(ichannel)%schur_projector_mode_response_negative_ema, &
                                  qs_ot_env(ichannel)%schur_projector_mode_trust_ema, &
                                  qs_ot_env(ichannel)%kpoint_occ_metric_residual_quality_ema, &
                                  qs_ot_env(ichannel)%schur_model_density_stationarity_delta, &
                                  qs_ot_env(ichannel)%schur_density_stationarity_quality, &
                                  qs_ot_env(ichannel)%schur_density_stationarity_quality_ema, &
                                  REAL(qs_ot_env(ichannel)%schur_density_stationarity_quality_samples, KIND=dp), &
                                  qs_ot_env(ichannel)%schur_projector_density_full_selected_stationarity_quality, &
                                  qs_ot_env(ichannel)%schur_projector_density_full_selected_stationarity_quality_ema, &
                                  REAL(qs_ot_env(ichannel)%schur_projector_density_full_selected_stationarity_samples, &
                                       KIND=dp), &
                                  qs_ot_env(ichannel)%schur_model_slope, &
                                  qs_ot_env(ichannel)%schur_model_quadratic]

         channel_source = qs_ot_env(ichannel)%state_allocated
         IF (channel_source) THEN
            CPASSERT(ASSOCIATED(qs_ot_env(ichannel)%para_env))
            channel_source = qs_ot_env(ichannel)%para_env%is_source()
         END IF
         source_count = MERGE(1, 0, channel_source)
         CALL global_para_env%sum(source_count)
         CPASSERT(source_count == 1)

         log_payload_int = MERGE(1, 0, log_payload)
         adaptive_log_payload_int = MERGE(1, 0, adaptive_log_payload)
         ot_pos_payload = qs_ot_env(ichannel)%ot_pos
         ot_energy_payload = qs_ot_env(ichannel)%ot_energy
         ot_grad_payload = qs_ot_env(ichannel)%ot_grad
         IF (.NOT. channel_source) THEN
            int_payload = 0
            log_payload_int = 0
            real_payload = 0.0_dp
            ot_pos_payload = 0.0_dp
            ot_energy_payload = 0.0_dp
            ot_grad_payload = 0.0_dp
            schur_int_payload = 0
            schur_real_payload = 0.0_dp
            adaptive_int_payload = 0
            adaptive_log_payload_int = 0
            adaptive_real_payload = 0.0_dp
         END IF
         CALL global_para_env%sum(int_payload)
         CALL global_para_env%sum(log_payload_int)
         CALL global_para_env%sum(real_payload)
         CALL global_para_env%sum(ot_pos_payload)
         CALL global_para_env%sum(ot_energy_payload)
         CALL global_para_env%sum(ot_grad_payload)
         CALL global_para_env%sum(schur_int_payload)
         CALL global_para_env%sum(schur_real_payload)
         CALL global_para_env%sum(adaptive_int_payload)
         CALL global_para_env%sum(adaptive_log_payload_int)
         CALL global_para_env%sum(adaptive_real_payload)
         log_payload = log_payload_int /= 0
         adaptive_log_payload = adaptive_log_payload_int /= 0

         nener = int_payload(14)
         ndiis = int_payload(15)
         has_rot_mat_x = log_payload(5)
         has_rot_mat_x_im = log_payload(6)
         has_ener_x = log_payload(7)
         has_matrix_dx = log_payload(10)
         has_matrix_dx_im = log_payload(11)
         has_matrix_gx_old = log_payload(12)
         has_matrix_gx_old_im = log_payload(13)
         has_rot_mat_dx = log_payload(14)
         has_rot_mat_dx_im = log_payload(15)
         has_rot_mat_gx_old = log_payload(16)
         has_rot_mat_gx_old_im = log_payload(17)
         has_ener_dx = log_payload(18)
         has_ener_gx_old = log_payload(19)
         has_diis_scalars = log_payload(20)
         has_matrix_history = log_payload(21)
         has_matrix_history_im = log_payload(22)
         has_rot_history = log_payload(23)
         has_rot_history_im = log_payload(24)
         has_ener_history = log_payload(25)
         has_matrix_x = log_payload(26)
         has_matrix_x_im = log_payload(27)
         has_lbfgs_scalars = log_payload(32)

         IF (ires > 0) THEN
            WRITE (ires) int_payload
            WRITE (ires) log_payload
            WRITE (ires) real_payload
            WRITE (ires) ot_pos_payload
            WRITE (ires) ot_energy_payload
            WRITE (ires) ot_grad_payload
            WRITE (ires) schur_int_payload
            WRITE (ires) schur_real_payload
            WRITE (ires) adaptive_int_payload
            WRITE (ires) adaptive_log_payload
            WRITE (ires) adaptive_real_payload
         END IF
         IF (has_rot_mat_x) CALL write_restart_dbcsr_matrix(qs_ot_env(ichannel)%rot_mat_x, &
                                                            qs_ot_env(ichannel), channel_source, &
                                                            fmat_template, ires)
         IF (has_rot_mat_x_im) CALL write_restart_dbcsr_matrix(qs_ot_env(ichannel)%rot_mat_x_im, &
                                                               qs_ot_env(ichannel), channel_source, &
                                                               fmat_template, ires)
         IF (has_matrix_x) CALL write_restart_dbcsr_matrix(qs_ot_env(ichannel)%matrix_x, &
                                                           qs_ot_env(ichannel), channel_source, &
                                                           fmat_template, ires)
         IF (has_matrix_x_im) CALL write_restart_dbcsr_matrix(qs_ot_env(ichannel)%matrix_x_im, &
                                                              qs_ot_env(ichannel), channel_source, &
                                                              fmat_template, ires)
         IF (has_ener_x) CALL write_restart_real_vector(qs_ot_env(ichannel)%ener_x, nener, &
                                                        channel_source, global_para_env, ires)
         IF (has_matrix_dx) CALL write_restart_dbcsr_matrix(qs_ot_env(ichannel)%matrix_dx, &
                                                            qs_ot_env(ichannel), channel_source, &
                                                            fmat_template, ires)
         IF (has_matrix_dx_im) CALL write_restart_dbcsr_matrix(qs_ot_env(ichannel)%matrix_dx_im, &
                                                               qs_ot_env(ichannel), channel_source, &
                                                               fmat_template, ires)
         IF (has_matrix_gx_old) CALL write_restart_dbcsr_matrix(qs_ot_env(ichannel)%matrix_gx_old, &
                                                                qs_ot_env(ichannel), channel_source, &
                                                                fmat_template, ires)
         IF (has_matrix_gx_old_im) CALL write_restart_dbcsr_matrix(qs_ot_env(ichannel)%matrix_gx_old_im, &
                                                                   qs_ot_env(ichannel), channel_source, &
                                                                   fmat_template, ires)
         IF (has_rot_mat_dx) CALL write_restart_dbcsr_matrix(qs_ot_env(ichannel)%rot_mat_dx, &
                                                             qs_ot_env(ichannel), channel_source, &
                                                             fmat_template, ires)
         IF (has_rot_mat_dx_im) CALL write_restart_dbcsr_matrix(qs_ot_env(ichannel)%rot_mat_dx_im, &
                                                                qs_ot_env(ichannel), channel_source, &
                                                                fmat_template, ires)
         IF (has_rot_mat_gx_old) CALL write_restart_dbcsr_matrix(qs_ot_env(ichannel)%rot_mat_gx_old, &
                                                                 qs_ot_env(ichannel), channel_source, &
                                                                 fmat_template, ires)
         IF (has_rot_mat_gx_old_im) CALL write_restart_dbcsr_matrix(qs_ot_env(ichannel)%rot_mat_gx_old_im, &
                                                                    qs_ot_env(ichannel), channel_source, &
                                                                    fmat_template, ires)
         IF (has_ener_dx) CALL write_restart_real_vector(qs_ot_env(ichannel)%ener_dx, nener, &
                                                         channel_source, global_para_env, ires)
         IF (has_ener_gx_old) CALL write_restart_real_vector(qs_ot_env(ichannel)%ener_gx_old, nener, &
                                                             channel_source, global_para_env, ires)
         IF (has_diis_scalars) THEN
            CALL write_restart_real_matrix(qs_ot_env(ichannel)%ls_diis, ndiis + 1, ndiis + 1, &
                                           channel_source, global_para_env, ires)
            CALL write_restart_real_matrix(qs_ot_env(ichannel)%lss_diis, ndiis + 1, ndiis + 1, &
                                           channel_source, global_para_env, ires)
            CALL write_restart_real_vector(qs_ot_env(ichannel)%c_diis, ndiis + 1, &
                                           channel_source, global_para_env, ires)
            CALL write_restart_real_vector(qs_ot_env(ichannel)%c_broy, ndiis, &
                                           channel_source, global_para_env, ires)
            CALL write_restart_real_vector(qs_ot_env(ichannel)%energy_h, ndiis, &
                                           channel_source, global_para_env, ires)
            CALL write_restart_integer_vector(qs_ot_env(ichannel)%ipivot, ndiis + 1, &
                                              channel_source, global_para_env, ires)
         END IF
         IF (has_matrix_history) THEN
            DO ihist = 1, ndiis
               NULLIFY (restart_matrix)
               IF (qs_ot_env(ichannel)%state_allocated) THEN
                  restart_matrix => qs_ot_env(ichannel)%matrix_h_e(ihist)%matrix
               END IF
               CALL write_restart_dbcsr_matrix(restart_matrix, qs_ot_env(ichannel), channel_source, &
                                               fmat_template, ires)
               NULLIFY (restart_matrix)
               IF (qs_ot_env(ichannel)%state_allocated) THEN
                  restart_matrix => qs_ot_env(ichannel)%matrix_h_x(ihist)%matrix
               END IF
               CALL write_restart_dbcsr_matrix(restart_matrix, qs_ot_env(ichannel), channel_source, &
                                               fmat_template, ires)
            END DO
         END IF
         IF (has_matrix_history_im) THEN
            DO ihist = 1, ndiis
               NULLIFY (restart_matrix)
               IF (qs_ot_env(ichannel)%state_allocated) THEN
                  restart_matrix => qs_ot_env(ichannel)%matrix_h_e_im(ihist)%matrix
               END IF
               CALL write_restart_dbcsr_matrix(restart_matrix, qs_ot_env(ichannel), channel_source, &
                                               fmat_template, ires)
               NULLIFY (restart_matrix)
               IF (qs_ot_env(ichannel)%state_allocated) THEN
                  restart_matrix => qs_ot_env(ichannel)%matrix_h_x_im(ihist)%matrix
               END IF
               CALL write_restart_dbcsr_matrix(restart_matrix, qs_ot_env(ichannel), channel_source, &
                                               fmat_template, ires)
            END DO
         END IF
         IF (has_rot_history) THEN
            DO ihist = 1, ndiis
               NULLIFY (restart_matrix)
               IF (qs_ot_env(ichannel)%state_allocated) THEN
                  restart_matrix => qs_ot_env(ichannel)%rot_mat_h_e(ihist)%matrix
               END IF
               CALL write_restart_dbcsr_matrix(restart_matrix, qs_ot_env(ichannel), channel_source, &
                                               fmat_template, ires)
               NULLIFY (restart_matrix)
               IF (qs_ot_env(ichannel)%state_allocated) THEN
                  restart_matrix => qs_ot_env(ichannel)%rot_mat_h_x(ihist)%matrix
               END IF
               CALL write_restart_dbcsr_matrix(restart_matrix, qs_ot_env(ichannel), channel_source, &
                                               fmat_template, ires)
            END DO
         END IF
         IF (has_rot_history_im) THEN
            DO ihist = 1, ndiis
               NULLIFY (restart_matrix)
               IF (qs_ot_env(ichannel)%state_allocated) THEN
                  restart_matrix => qs_ot_env(ichannel)%rot_mat_h_e_im(ihist)%matrix
               END IF
               CALL write_restart_dbcsr_matrix(restart_matrix, qs_ot_env(ichannel), channel_source, &
                                               fmat_template, ires)
               NULLIFY (restart_matrix)
               IF (qs_ot_env(ichannel)%state_allocated) THEN
                  restart_matrix => qs_ot_env(ichannel)%rot_mat_h_x_im(ihist)%matrix
               END IF
               CALL write_restart_dbcsr_matrix(restart_matrix, qs_ot_env(ichannel), channel_source, &
                                               fmat_template, ires)
            END DO
         END IF
         IF (has_ener_history) THEN
            CALL write_restart_real_matrix(qs_ot_env(ichannel)%ener_h_e, ndiis, nener, &
                                           channel_source, global_para_env, ires)
            CALL write_restart_real_matrix(qs_ot_env(ichannel)%ener_h_x, ndiis, nener, &
                                           channel_source, global_para_env, ires)
         END IF
         IF (has_lbfgs_scalars) THEN
            CALL write_restart_real_vector(qs_ot_env(ichannel)%lbfgs_rho, ndiis - 1, &
                                           channel_source, global_para_env, ires)
            CALL write_restart_real_vector(qs_ot_env(ichannel)%lbfgs_yy, ndiis - 1, &
                                           channel_source, global_para_env, ires)
         END IF
      END DO

   END SUBROUTINE write_kpoint_ot_restart_payload

! **************************************************************************************************
!> \brief gather one channel-local DBCSR matrix and write it through the global restart FM layer
!> \param matrix matrix to write
!> \param qs_ot_env channel owning the matrix distribution
!> \param channel_source true on the source rank of the owning k-point group
!> \param fmat_template global restart FM distribution
!> \param ires restart unit
! **************************************************************************************************
   SUBROUTINE write_restart_dbcsr_matrix(matrix, qs_ot_env, channel_source, fmat_template, ires)
      TYPE(dbcsr_type), POINTER                          :: matrix
      TYPE(qs_ot_type), INTENT(IN)                       :: qs_ot_env
      LOGICAL, INTENT(IN)                                :: channel_source
      TYPE(cp_fm_type), INTENT(IN)                       :: fmat_template
      INTEGER, INTENT(IN)                                :: ires

      INTEGER                                            :: ncol, nrow
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: dense
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type)                                   :: fmat, global_fmat
      TYPE(mp_para_env_type), POINTER                    :: global_para_env

      NULLIFY (fm_struct, global_para_env)
      global_para_env => fmat_template%matrix_struct%para_env
      nrow = 0
      ncol = 0
      IF (channel_source) THEN
         CPASSERT(ASSOCIATED(matrix))
         CALL dbcsr_get_info(matrix, nfullrows_total=nrow, nfullcols_total=ncol)
      END IF
      CALL global_para_env%sum(nrow)
      CALL global_para_env%sum(ncol)
      CPASSERT(nrow > 0 .AND. ncol > 0)
      ALLOCATE (dense(nrow, ncol))
      dense = 0.0_dp

      IF (qs_ot_env%state_allocated) THEN
         CPASSERT(ASSOCIATED(matrix))
         CALL cp_fm_struct_create(fm_struct, para_env=qs_ot_env%para_env, &
                                  context=qs_ot_env%blacs_env, &
                                  nrow_global=nrow, ncol_global=ncol)
         CALL cp_fm_create(fmat, fm_struct)
         CALL copy_dbcsr_to_fm(matrix, fmat)
         CALL cp_fm_get_submatrix(fmat, dense)
         IF (.NOT. channel_source) dense = 0.0_dp
         CALL cp_fm_release(fmat)
         CALL cp_fm_struct_release(fm_struct)
      END IF
      CALL global_para_env%sum(dense)

      IF (ires > 0) WRITE (ires) nrow, ncol
      CALL cp_fm_create(global_fmat, fmat_template%matrix_struct, nrow=nrow, ncol=ncol)
      CALL cp_fm_set_submatrix(global_fmat, dense)
      CALL cp_fm_write_unformatted(global_fmat, ires)
      CALL cp_fm_release(global_fmat)
      DEALLOCATE (dense)

   END SUBROUTINE write_restart_dbcsr_matrix

! **************************************************************************************************
!> \brief gather and write one channel-local real vector
!> \param vector ...
!> \param n ...
!> \param channel_source ...
!> \param para_env ...
!> \param ires ...
! **************************************************************************************************
   SUBROUTINE write_restart_real_vector(vector, n, channel_source, para_env, ires)
      REAL(KIND=dp), DIMENSION(:), POINTER               :: vector
      INTEGER, INTENT(IN)                                :: n
      LOGICAL, INTENT(IN)                                :: channel_source
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(IN)                                :: ires

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: buffer

      CPASSERT(n > 0)
      ALLOCATE (buffer(n))
      buffer = 0.0_dp
      IF (channel_source) THEN
         CPASSERT(ASSOCIATED(vector))
         CPASSERT(SIZE(vector) == n)
         buffer(:) = vector
      END IF
      CALL para_env%sum(buffer)
      IF (ires > 0) WRITE (ires) buffer
      DEALLOCATE (buffer)

   END SUBROUTINE write_restart_real_vector

! **************************************************************************************************
!> \brief gather and write one channel-local real matrix
!> \param matrix ...
!> \param nrow ...
!> \param ncol ...
!> \param channel_source ...
!> \param para_env ...
!> \param ires ...
! **************************************************************************************************
   SUBROUTINE write_restart_real_matrix(matrix, nrow, ncol, channel_source, para_env, ires)
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: matrix
      INTEGER, INTENT(IN)                                :: nrow, ncol
      LOGICAL, INTENT(IN)                                :: channel_source
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(IN)                                :: ires

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: buffer

      CPASSERT(nrow > 0 .AND. ncol > 0)
      ALLOCATE (buffer(nrow, ncol))
      buffer = 0.0_dp
      IF (channel_source) THEN
         CPASSERT(ASSOCIATED(matrix))
         CPASSERT(SIZE(matrix, 1) == nrow .AND. SIZE(matrix, 2) == ncol)
         buffer(:, :) = matrix
      END IF
      CALL para_env%sum(buffer)
      IF (ires > 0) WRITE (ires) buffer
      DEALLOCATE (buffer)

   END SUBROUTINE write_restart_real_matrix

! **************************************************************************************************
!> \brief gather and write one channel-local integer vector
!> \param vector ...
!> \param n ...
!> \param channel_source ...
!> \param para_env ...
!> \param ires ...
! **************************************************************************************************
   SUBROUTINE write_restart_integer_vector(vector, n, channel_source, para_env, ires)
      INTEGER, DIMENSION(:), POINTER                     :: vector
      INTEGER, INTENT(IN)                                :: n
      LOGICAL, INTENT(IN)                                :: channel_source
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(IN)                                :: ires

      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: buffer

      CPASSERT(n > 0)
      ALLOCATE (buffer(n))
      buffer = 0
      IF (channel_source) THEN
         CPASSERT(ASSOCIATED(vector))
         CPASSERT(SIZE(vector) == n)
         buffer(:) = vector
      END IF
      CALL para_env%sum(buffer)
      IF (ires > 0) WRITE (ires) buffer
      DEALLOCATE (buffer)

   END SUBROUTINE write_restart_integer_vector

! **************************************************************************************************
!> \brief read optional k-point OT restart state after the density payload
!> \param rst_unit restart unit
!> \param para_env parallel environment
!> \param fmat restart FM template
!> \param scf_env SCF environment receiving the cache
!> \param cell current simulation cell
!> \param particle_set current atomic geometry
!> \param smear current smearing objective
! **************************************************************************************************
   SUBROUTINE read_kpoint_ot_restart_payload(rst_unit, para_env, fmat, scf_env, cell, particle_set, smear)
      INTEGER, INTENT(IN)                                :: rst_unit
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(cp_fm_type), INTENT(INOUT)                    :: fmat
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(smear_type), INTENT(IN)                       :: smear

      CHARACTER(LEN=20)                                  :: requested_ot_method_payload
      CHARACTER(LEN=20), &
         DIMENSION(qs_ot_restart_signature_nchar)        :: settings_char_payload
      CHARACTER(LEN=LEN(ot_restart_marker))              :: marker
      INTEGER                                            :: adaptive_ios, iatom, ichannel, ihist, &
                                                            ios, natom_restart, nchannel, ndiis, &
                                                            nener, version_read
      INTEGER, DIMENSION(11)                             :: schur_int_payload_v4
      INTEGER, DIMENSION(12)                             :: schur_int_payload_v5
      INTEGER, DIMENSION(13)                             :: adaptive_int_payload_v20, &
                                                            int_payload_v10, schur_int_payload_v6
      INTEGER, DIMENSION(14)                             :: adaptive_int_payload_v21, &
                                                            schur_int_payload_v8
      INTEGER, DIMENSION(15)                             :: int_payload_v35
      INTEGER, DIMENSION(16)                             :: int_payload
      INTEGER, DIMENSION(17)                             :: schur_int_payload_v12
      INTEGER, DIMENSION(18)                             :: schur_int_payload_v13
      INTEGER, DIMENSION(19)                             :: schur_int_payload
      INTEGER, DIMENSION(3)                              :: smear_int_payload
      INTEGER, DIMENSION(qs_ot_restart_adaptive_nint)    :: adaptive_int_payload
      INTEGER, DIMENSION(qs_ot_restart_signature_nint)   :: settings_int_payload
      LOGICAL                                            :: geometry_compatible, payload_found, &
                                                            smearing_compatible
      LOGICAL, DIMENSION(11)                             :: adaptive_log_payload_v20
      LOGICAL, DIMENSION(12)                             :: adaptive_log_payload_v21
      LOGICAL, DIMENSION(27)                             :: log_payload_v10
      LOGICAL, DIMENSION(28)                             :: log_payload_v18
      LOGICAL, DIMENSION(29)                             :: log_payload_v35
      LOGICAL, DIMENSION(31)                             :: log_payload_v37
      LOGICAL, DIMENSION(32)                             :: log_payload
      LOGICAL, DIMENSION(qs_ot_restart_adaptive_nlog)    :: adaptive_log_payload
      LOGICAL, DIMENSION(qs_ot_restart_signature_nlog)   :: settings_log_payload
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: current_positions, restart_positions
      REAL(KIND=dp), DIMENSION(118)                      :: schur_real_payload_v9
      REAL(KIND=dp), DIMENSION(124)                      :: schur_real_payload_v11
      REAL(KIND=dp), DIMENSION(126)                      :: schur_real_payload_v12
      REAL(KIND=dp), DIMENSION(133)                      :: schur_real_payload_v13
      REAL(KIND=dp), DIMENSION(139)                      :: schur_real_payload
      REAL(KIND=dp), DIMENSION(19)                       :: real_payload_v26
      REAL(KIND=dp), DIMENSION(20)                       :: real_payload_v16, real_payload_v35
      REAL(KIND=dp), DIMENSION(24)                       :: real_payload_v36
      REAL(KIND=dp), DIMENSION(25)                       :: real_payload
      REAL(KIND=dp), DIMENSION(3)                        :: smear_real_payload
      REAL(KIND=dp), DIMENSION(3, 3)                     :: restart_hmat
      REAL(KIND=dp), DIMENSION(51)                       :: schur_real_payload_v4
      REAL(KIND=dp), DIMENSION(58)                       :: schur_real_payload_v5
      REAL(KIND=dp), DIMENSION(65)                       :: schur_real_payload_v6
      REAL(KIND=dp), DIMENSION(75)                       :: adaptive_real_payload_v20
      REAL(KIND=dp), DIMENSION(8)                        :: real_payload_v10
      REAL(KIND=dp), DIMENSION(87)                       :: schur_real_payload_v7
      REAL(KIND=dp), DIMENSION(88)                       :: adaptive_real_payload_v21
      REAL(KIND=dp), DIMENSION(89)                       :: adaptive_real_payload_v23, &
                                                            schur_real_payload_v8
      REAL(KIND=dp), DIMENSION(96)                       :: adaptive_real_payload_v24
      REAL(KIND=dp), &
         DIMENSION(qs_ot_restart_adaptive_nreal)         :: adaptive_real_payload
      REAL(KIND=dp), &
         DIMENSION(qs_ot_restart_signature_nreal)        :: settings_real_payload

      payload_found = .FALSE.
      marker = ""
      IF (para_env%is_source()) THEN
         READ (rst_unit, IOSTAT=ios) marker
         payload_found = (ios == 0 .AND. TRIM(marker) == ot_restart_marker)
      END IF
      CALL para_env%bcast(payload_found)
      IF (.NOT. payload_found) RETURN

      version_read = 0
      nchannel = 0
      IF (para_env%is_source()) READ (rst_unit) version_read, nchannel
      CALL para_env%bcast(version_read)
      CALL para_env%bcast(nchannel)
      IF (version_read < 3 .OR. version_read > ot_restart_version) THEN
         CALL cp_warn(__LOCATION__, "Ignoring unsupported k-point OT restart payload version")
         RETURN
      END IF
      IF (version_read < 26) THEN
         CALL cp_warn(__LOCATION__, &
                      "Ignoring k-point OT minimizer payload without a geometry signature")
         RETURN
      END IF
      IF (version_read < 33) THEN
         CALL cp_warn(__LOCATION__, &
                      "Ignoring k-point OT minimizer payload without a settings signature")
         RETURN
      END IF
      IF (version_read < 35) THEN
         CALL cp_warn(__LOCATION__, &
                      "Ignoring k-point OT minimizer payload with legacy DIIS rotation history")
         RETURN
      END IF

      CPASSERT(ASSOCIATED(cell))
      CPASSERT(ASSOCIATED(particle_set))
      natom_restart = 0
      ios = 0
      IF (para_env%is_source()) READ (rst_unit, IOSTAT=ios) natom_restart
      CALL para_env%bcast(ios)
      CALL para_env%bcast(natom_restart)
      IF (ios /= 0 .OR. natom_restart <= 0) THEN
         CALL cp_warn(__LOCATION__, "Ignoring malformed k-point OT restart geometry signature")
         RETURN
      END IF
      IF (natom_restart /= SIZE(particle_set)) THEN
         CALL cp_warn(__LOCATION__, &
                      "Ignoring k-point OT minimizer payload from a different cell or geometry")
         RETURN
      END IF
      ALLOCATE (restart_positions(3, natom_restart))
      restart_hmat(:, :) = 0.0_dp
      restart_positions(:, :) = 0.0_dp
      smear_int_payload(:) = 0
      smear_real_payload(:) = 0.0_dp
      settings_char_payload(:) = ""
      settings_int_payload(:) = 0
      settings_log_payload(:) = .FALSE.
      settings_real_payload(:) = 0.0_dp
      requested_ot_method_payload = ""
      IF (para_env%is_source()) THEN
         READ (rst_unit, IOSTAT=ios) restart_hmat
         IF (ios == 0) READ (rst_unit, IOSTAT=ios) restart_positions
         IF (ios == 0) READ (rst_unit, IOSTAT=ios) smear_int_payload
         IF (ios == 0) READ (rst_unit, IOSTAT=ios) smear_real_payload
         IF (ios == 0) READ (rst_unit, IOSTAT=ios) settings_char_payload
         IF (ios == 0) READ (rst_unit, IOSTAT=ios) settings_int_payload
         IF (ios == 0) READ (rst_unit, IOSTAT=ios) settings_log_payload
         IF (ios == 0) READ (rst_unit, IOSTAT=ios) settings_real_payload
         IF (version_read >= 39 .AND. ios == 0) THEN
            READ (rst_unit, IOSTAT=ios) requested_ot_method_payload
         END IF
      END IF
      CALL para_env%bcast(ios)
      CALL para_env%bcast(restart_hmat)
      CALL para_env%bcast(restart_positions)
      CALL para_env%bcast(smear_int_payload)
      CALL para_env%bcast(smear_real_payload)
      CALL para_env%bcast(settings_char_payload)
      CALL para_env%bcast(settings_int_payload)
      CALL para_env%bcast(settings_log_payload)
      CALL para_env%bcast(settings_real_payload)
      IF (version_read >= 39) THEN
         CALL para_env%bcast(requested_ot_method_payload)
      ELSE
         requested_ot_method_payload = settings_char_payload(1)
      END IF
      IF (ios /= 0) THEN
         DEALLOCATE (restart_positions)
         CALL cp_warn(__LOCATION__, "Ignoring malformed k-point OT restart geometry signature")
         RETURN
      END IF

      ALLOCATE (current_positions(3, SIZE(particle_set)))
      DO iatom = 1, SIZE(particle_set)
         current_positions(:, iatom) = particle_set(iatom)%r
      END DO
      geometry_compatible = kpoint_ot_restart_geometry_compatible( &
                            cell%hmat, current_positions, restart_hmat, restart_positions)
      DEALLOCATE (current_positions)
      DEALLOCATE (restart_positions)
      IF (.NOT. geometry_compatible) THEN
         CALL cp_warn(__LOCATION__, &
                      "Ignoring k-point OT minimizer payload from a different cell or geometry")
         RETURN
      END IF

      smearing_compatible = ALL(smear_int_payload(1:2) >= 0 .AND. &
                                smear_int_payload(1:2) <= 1)
      IF (smearing_compatible) THEN
         smearing_compatible = kpoint_ot_restart_smearing_compatible( &
                               smear, smear_int_payload(1) == 1, smear_int_payload(2) == 1, &
                               smear_int_payload(3), smear_real_payload(1), &
                               smear_real_payload(2), smear_real_payload(3))
      END IF
      IF (.NOT. smearing_compatible) THEN
         CALL cp_warn(__LOCATION__, &
                      "Ignoring k-point OT minimizer payload from a different smearing objective")
         RETURN
      END IF

      CALL qs_scf_ot_restart_clear(scf_env)
      scf_env%qs_ot_restart%valid = .TRUE.
      scf_env%qs_ot_restart%format_version = version_read
      scf_env%qs_ot_restart%nchannel = nchannel
      scf_env%qs_ot_restart%requested_ot_method = requested_ot_method_payload
      scf_env%qs_ot_restart%settings_char = settings_char_payload
      scf_env%qs_ot_restart%settings_int = settings_int_payload
      scf_env%qs_ot_restart%settings_log = settings_log_payload
      scf_env%qs_ot_restart%settings_real = settings_real_payload
      ALLOCATE (scf_env%qs_ot_restart%channel(nchannel))
      DO ichannel = 1, nchannel
         adaptive_int_payload(:) = 0
         adaptive_log_payload(:) = .FALSE.
         adaptive_real_payload(:) = 0.0_dp
         log_payload(:) = .FALSE.
         real_payload(:) = 0.0_dp
         IF (para_env%is_source()) THEN
            IF (version_read >= 38) THEN
               READ (rst_unit) int_payload
               READ (rst_unit) log_payload
               READ (rst_unit) real_payload
            ELSE IF (version_read == 37) THEN
               READ (rst_unit) int_payload
               READ (rst_unit) log_payload_v37
               READ (rst_unit) real_payload
               log_payload(:) = .FALSE.
               log_payload(1:31) = log_payload_v37
            ELSE IF (version_read == 36) THEN
               READ (rst_unit) int_payload
               READ (rst_unit) log_payload_v37
               READ (rst_unit) real_payload_v36
               log_payload(1:31) = log_payload_v37
               real_payload(1:24) = real_payload_v36
            ELSE IF (version_read >= 27) THEN
               READ (rst_unit) int_payload_v35
               READ (rst_unit) log_payload_v35
               READ (rst_unit) real_payload_v35
               int_payload(:) = 0
               int_payload(1:15) = int_payload_v35
               log_payload(:) = .FALSE.
               log_payload(1:29) = log_payload_v35
               log_payload(30) = .TRUE.
               real_payload(:) = 0.0_dp
               real_payload(1:20) = real_payload_v35
               real_payload(22:24) = -1.0_dp
            ELSE IF (version_read >= 19) THEN
               READ (rst_unit) int_payload
               READ (rst_unit) log_payload_v37
               READ (rst_unit) real_payload_v26
               log_payload(1:31) = log_payload_v37
               real_payload(1:19) = real_payload_v26
            ELSE IF (version_read >= 17) THEN
               READ (rst_unit) int_payload
               READ (rst_unit) log_payload_v18
               READ (rst_unit) real_payload_v26
               log_payload(:) = .FALSE.
               log_payload(1:28) = log_payload_v18
               real_payload(1:19) = real_payload_v26
            ELSE IF (version_read == 16) THEN
               READ (rst_unit) int_payload
               READ (rst_unit) log_payload_v18
               READ (rst_unit) real_payload_v16
               log_payload(:) = .FALSE.
               log_payload(1:28) = log_payload_v18
               real_payload(1:19) = real_payload_v16(1:19)
            ELSE IF (version_read >= 11) THEN
               READ (rst_unit) int_payload
               READ (rst_unit) log_payload_v18
               READ (rst_unit) real_payload_v26
               log_payload(:) = .FALSE.
               log_payload(1:28) = log_payload_v18
               real_payload(1:19) = real_payload_v26
            ELSE
               READ (rst_unit) int_payload_v10
               READ (rst_unit) log_payload_v10
               READ (rst_unit) real_payload_v10
               int_payload(:) = 0
               int_payload(1:11) = int_payload_v10(1:11)
               int_payload(14:15) = int_payload_v10(12:13)
               log_payload(:) = .FALSE.
               log_payload(1:27) = log_payload_v10
               real_payload(:) = 0.0_dp
               real_payload(1:8) = real_payload_v10
               real_payload(16) = 1.0_dp
               real_payload(17) = 1.0_dp
            END IF
            READ (rst_unit) scf_env%qs_ot_restart%channel(ichannel)%ot_pos
            READ (rst_unit) scf_env%qs_ot_restart%channel(ichannel)%ot_energy
            READ (rst_unit) scf_env%qs_ot_restart%channel(ichannel)%ot_grad
            IF (version_read >= 14) THEN
               READ (rst_unit) schur_int_payload
               READ (rst_unit) schur_real_payload
            ELSE IF (version_read == 13) THEN
               READ (rst_unit) schur_int_payload_v13
               READ (rst_unit) schur_real_payload_v13
               schur_int_payload(:) = 0
               schur_real_payload(:) = 0.0_dp
               schur_int_payload(1:18) = schur_int_payload_v13
               schur_real_payload(1:133) = schur_real_payload_v13
            ELSE IF (version_read == 12) THEN
               READ (rst_unit) schur_int_payload_v12
               READ (rst_unit) schur_real_payload_v12
               schur_int_payload(:) = 0
               schur_real_payload(:) = 0.0_dp
               schur_int_payload(1:17) = schur_int_payload_v12
               schur_real_payload(1:126) = schur_real_payload_v12
            ELSE IF (version_read >= 10) THEN
               READ (rst_unit) schur_int_payload_v12
               READ (rst_unit) schur_real_payload_v11
               schur_int_payload(:) = 0
               schur_real_payload(:) = 0.0_dp
               schur_int_payload(1:17) = schur_int_payload_v12
               schur_real_payload(1:124) = schur_real_payload_v11
            ELSE IF (version_read == 9) THEN
               READ (rst_unit) schur_int_payload_v12
               READ (rst_unit) schur_real_payload_v9
               schur_int_payload(:) = 0
               schur_real_payload(:) = 0.0_dp
               schur_int_payload(1:17) = schur_int_payload_v12
               schur_real_payload(1:118) = schur_real_payload_v9
            ELSE IF (version_read == 8) THEN
               READ (rst_unit) schur_int_payload_v8
               READ (rst_unit) schur_real_payload_v8
               schur_int_payload(:) = 0
               schur_real_payload(:) = 0.0_dp
               schur_int_payload(1:14) = schur_int_payload_v8
               schur_real_payload(1:89) = schur_real_payload_v8
            ELSE IF (version_read == 7) THEN
               READ (rst_unit) schur_int_payload_v8
               READ (rst_unit) schur_real_payload_v7
               schur_int_payload(:) = 0
               schur_int_payload(1:14) = schur_int_payload_v8
               schur_real_payload(:) = 0.0_dp
               schur_real_payload(1:87) = schur_real_payload_v7
               schur_real_payload(89) = 1.0_dp
            ELSE IF (version_read == 6) THEN
               READ (rst_unit) schur_int_payload_v6
               READ (rst_unit) schur_real_payload_v6
               schur_int_payload(:) = 0
               schur_real_payload(:) = 0.0_dp
               schur_int_payload(1:13) = schur_int_payload_v6
               schur_real_payload(1:65) = schur_real_payload_v6
            ELSE IF (version_read == 5) THEN
               READ (rst_unit) schur_int_payload_v5
               READ (rst_unit) schur_real_payload_v5
               schur_int_payload(:) = 0
               schur_real_payload(:) = 0.0_dp
               schur_int_payload(1:12) = schur_int_payload_v5
               schur_real_payload(1:58) = schur_real_payload_v5
            ELSE IF (version_read == 4) THEN
               READ (rst_unit) schur_int_payload_v4
               READ (rst_unit) schur_real_payload_v4
               schur_int_payload(:) = 0
               schur_real_payload(:) = 0.0_dp
               schur_int_payload(1:4) = schur_int_payload_v4(1:4)
               schur_int_payload(6:12) = schur_int_payload_v4(5:11)
               schur_real_payload(1:51) = schur_real_payload_v4
            END IF
            IF (version_read >= 25) THEN
               READ (rst_unit) adaptive_int_payload
               READ (rst_unit) adaptive_log_payload
               READ (rst_unit) adaptive_real_payload
            ELSE IF (version_read == 24) THEN
               READ (rst_unit) adaptive_int_payload
               READ (rst_unit) adaptive_log_payload
               READ (rst_unit) adaptive_real_payload_v24
               adaptive_real_payload(1:96) = adaptive_real_payload_v24
            ELSE IF (version_read >= 22) THEN
               READ (rst_unit) adaptive_int_payload
               READ (rst_unit) adaptive_log_payload
               READ (rst_unit) adaptive_real_payload_v23
               adaptive_real_payload(1:89) = adaptive_real_payload_v23
            ELSE IF (version_read == 21) THEN
               READ (rst_unit) adaptive_int_payload_v21
               READ (rst_unit) adaptive_log_payload_v21
               READ (rst_unit) adaptive_real_payload_v21
               adaptive_int_payload(1:14) = adaptive_int_payload_v21
               adaptive_log_payload(1:12) = adaptive_log_payload_v21
               adaptive_real_payload(1:88) = adaptive_real_payload_v21
            ELSE IF (version_read == 20) THEN
               ! Early version-20 files predate the adaptive payload, while a
               ! later development layout reused the same version number for
               ! both the 13/11/75 and 14/12/88 records.  Probe from largest to
               ! smallest and rewind when the record is already the following
               ! matrix payload.
               READ (rst_unit, IOSTAT=adaptive_ios) adaptive_int_payload_v21
               IF (adaptive_ios == 0) THEN
                  READ (rst_unit) adaptive_log_payload_v21
                  READ (rst_unit) adaptive_real_payload_v21
                  adaptive_int_payload(1:14) = adaptive_int_payload_v21
                  adaptive_log_payload(1:12) = adaptive_log_payload_v21
                  adaptive_real_payload(1:88) = adaptive_real_payload_v21
               ELSE
                  BACKSPACE (UNIT=rst_unit, IOSTAT=adaptive_ios)
                  IF (adaptive_ios /= 0) THEN
                     CALL cp_abort(__LOCATION__, "Failed to rewind a legacy version-20 k-point OT restart.")
                  END IF
                  READ (rst_unit, IOSTAT=adaptive_ios) adaptive_int_payload_v20
                  IF (adaptive_ios == 0) THEN
                     READ (rst_unit) adaptive_log_payload_v20
                     READ (rst_unit) adaptive_real_payload_v20
                     adaptive_int_payload(1:13) = adaptive_int_payload_v20
                     adaptive_log_payload(1:11) = adaptive_log_payload_v20
                     adaptive_real_payload(1:75) = adaptive_real_payload_v20
                  ELSE
                     BACKSPACE (UNIT=rst_unit, IOSTAT=adaptive_ios)
                     IF (adaptive_ios /= 0) THEN
                        CALL cp_abort(__LOCATION__, "Failed to rewind a legacy version-20 k-point OT restart.")
                     END IF
                  END IF
               END IF
            END IF
         END IF
         IF (version_read < 4) THEN
            schur_int_payload(:) = 0
            schur_real_payload(:) = 0.0_dp
         END IF
         IF (version_read < 8) schur_real_payload(89) = 1.0_dp
         IF (version_read < 9) THEN
            schur_real_payload(95) = 1.0_dp
            schur_real_payload(96) = 1.0_dp
            schur_real_payload(105) = 1.0_dp
            schur_real_payload(113) = 1.0_dp
            schur_real_payload(114) = 1.0_dp
         END IF
         IF (version_read < 12) schur_real_payload(126) = 1.0_dp
         IF (version_read < 13) schur_real_payload(131) = 1.0_dp
         IF (version_read < 14) schur_real_payload(137) = 1.0_dp
         CALL para_env%bcast(int_payload)
         CALL para_env%bcast(log_payload)
         CALL para_env%bcast(real_payload)
         CALL para_env%bcast(scf_env%qs_ot_restart%channel(ichannel)%ot_pos)
         CALL para_env%bcast(scf_env%qs_ot_restart%channel(ichannel)%ot_energy)
         CALL para_env%bcast(scf_env%qs_ot_restart%channel(ichannel)%ot_grad)
         CALL para_env%bcast(schur_int_payload)
         CALL para_env%bcast(schur_real_payload)
         CALL para_env%bcast(adaptive_int_payload)
         CALL para_env%bcast(adaptive_log_payload)
         CALL para_env%bcast(adaptive_real_payload)

         scf_env%qs_ot_restart%channel(ichannel)%adaptive_int = adaptive_int_payload
         scf_env%qs_ot_restart%channel(ichannel)%adaptive_log = adaptive_log_payload
         scf_env%qs_ot_restart%channel(ichannel)%adaptive_real = adaptive_real_payload

         scf_env%qs_ot_restart%channel(ichannel)%valid = log_payload(1) .AND. &
                                                         int_payload(1) == ichannel .AND. &
                                                         int_payload(5) > 0
         scf_env%qs_ot_restart%channel(ichannel)%spin_index = int_payload(2)
         scf_env%qs_ot_restart%channel(ichannel)%kpoint_index = int_payload(3)
         scf_env%qs_ot_restart%channel(ichannel)%local_kpoint_index = int_payload(4)
         scf_env%qs_ot_restart%channel(ichannel)%nmo = int_payload(5)
         scf_env%qs_ot_restart%channel(ichannel)%nmo_ortho = int_payload(6)
         scf_env%qs_ot_restart%channel(ichannel)%line_search_left = int_payload(7)
         scf_env%qs_ot_restart%channel(ichannel)%line_search_right = int_payload(8)
         scf_env%qs_ot_restart%channel(ichannel)%line_search_mid = int_payload(9)
         scf_env%qs_ot_restart%channel(ichannel)%line_search_count = int_payload(10)
         scf_env%qs_ot_restart%channel(ichannel)%diis_iter = int_payload(11)
         scf_env%qs_ot_restart%channel(ichannel)%mermin_ls_best_sample_index = int_payload(12)
         scf_env%qs_ot_restart%channel(ichannel)%mermin_ls_accepted_quality_samples = int_payload(13)
         nener = int_payload(14)
         ndiis = int_payload(15)
         scf_env%qs_ot_restart%channel(ichannel)%diis_m = ndiis
         scf_env%qs_ot_restart%channel(ichannel)%kpoint_subspace_accepted_samples = int_payload(16)
         scf_env%qs_ot_restart%channel(ichannel)%has_complex_kpoint_state = log_payload(2)
         scf_env%qs_ot_restart%channel(ichannel)%do_rotation = log_payload(3)
         scf_env%qs_ot_restart%channel(ichannel)%do_ener = log_payload(4)
         scf_env%qs_ot_restart%channel(ichannel)%has_rot_mat_x = log_payload(5)
         scf_env%qs_ot_restart%channel(ichannel)%has_rot_mat_x_im = log_payload(6)
         scf_env%qs_ot_restart%channel(ichannel)%has_ener_x = log_payload(7)
         scf_env%qs_ot_restart%channel(ichannel)%line_search_might_be_done = log_payload(8)
         scf_env%qs_ot_restart%channel(ichannel)%energy_only = log_payload(9)
         scf_env%qs_ot_restart%channel(ichannel)%has_matrix_dx = log_payload(10)
         scf_env%qs_ot_restart%channel(ichannel)%has_matrix_dx_im = log_payload(11)
         scf_env%qs_ot_restart%channel(ichannel)%has_matrix_gx_old = log_payload(12)
         scf_env%qs_ot_restart%channel(ichannel)%has_matrix_gx_old_im = log_payload(13)
         scf_env%qs_ot_restart%channel(ichannel)%has_rot_mat_dx = log_payload(14)
         scf_env%qs_ot_restart%channel(ichannel)%has_rot_mat_dx_im = log_payload(15)
         scf_env%qs_ot_restart%channel(ichannel)%has_rot_mat_gx_old = log_payload(16)
         scf_env%qs_ot_restart%channel(ichannel)%has_rot_mat_gx_old_im = log_payload(17)
         scf_env%qs_ot_restart%channel(ichannel)%has_ener_dx = log_payload(18)
         scf_env%qs_ot_restart%channel(ichannel)%has_ener_gx_old = log_payload(19)
         scf_env%qs_ot_restart%channel(ichannel)%has_diis_scalars = log_payload(20)
         scf_env%qs_ot_restart%channel(ichannel)%has_matrix_history = log_payload(21)
         scf_env%qs_ot_restart%channel(ichannel)%has_matrix_history_im = log_payload(22)
         scf_env%qs_ot_restart%channel(ichannel)%has_rot_history = log_payload(23)
         scf_env%qs_ot_restart%channel(ichannel)%has_rot_history_im = log_payload(24)
         scf_env%qs_ot_restart%channel(ichannel)%has_ener_history = log_payload(25)
         scf_env%qs_ot_restart%channel(ichannel)%has_matrix_x = log_payload(26)
         scf_env%qs_ot_restart%channel(ichannel)%has_matrix_x_im = log_payload(27)
         scf_env%qs_ot_restart%channel(ichannel)%mermin_ls_sample_pending = log_payload(28)
         scf_env%qs_ot_restart%channel(ichannel)%ener_direction_linear = log_payload(29)
         scf_env%qs_ot_restart%channel(ichannel)%diis_last_candidate_sd = log_payload(30)
         scf_env%qs_ot_restart%channel(ichannel)%diis_last_candidate_response = log_payload(31)
         scf_env%qs_ot_restart%channel(ichannel)%has_lbfgs_scalars = log_payload(32)
         scf_env%qs_ot_restart%channel(ichannel)%kpoint_weight = real_payload(1)
         scf_env%qs_ot_restart%channel(ichannel)%delta = real_payload(2)
         scf_env%qs_ot_restart%channel(ichannel)%gnorm = real_payload(3)
         scf_env%qs_ot_restart%channel(ichannel)%gnorm_old = real_payload(4)
         scf_env%qs_ot_restart%channel(ichannel)%etotal = real_payload(5)
         scf_env%qs_ot_restart%channel(ichannel)%gradient = real_payload(6)
         scf_env%qs_ot_restart%channel(ichannel)%ds_min = real_payload(7)
         scf_env%qs_ot_restart%channel(ichannel)%broyden_adaptive_sigma = real_payload(8)
         scf_env%qs_ot_restart%channel(ichannel)%mermin_ls_reference_energy = real_payload(9)
         scf_env%qs_ot_restart%channel(ichannel)%mermin_ls_best_sample_energy = real_payload(10)
         scf_env%qs_ot_restart%channel(ichannel)%mermin_ls_best_sample_drop = real_payload(11)
         scf_env%qs_ot_restart%channel(ichannel)%mermin_ls_model_drop = real_payload(12)
         scf_env%qs_ot_restart%channel(ichannel)%mermin_ls_model_quality = real_payload(13)
         scf_env%qs_ot_restart%channel(ichannel)%mermin_ls_best_sample_pos = real_payload(14)
         scf_env%qs_ot_restart%channel(ichannel)%mermin_ls_next_pos = real_payload(15)
         scf_env%qs_ot_restart%channel(ichannel)%mermin_ls_accepted_quality = real_payload(16)
         scf_env%qs_ot_restart%channel(ichannel)%mermin_ls_accepted_quality_ema = real_payload(17)
         scf_env%qs_ot_restart%channel(ichannel)%mermin_ls_accepted_gap = real_payload(18)
         scf_env%qs_ot_restart%channel(ichannel)%mermin_ls_accepted_gap_ema = real_payload(19)
         scf_env%qs_ot_restart%channel(ichannel)%mermin_ls_reference_delta = real_payload(20)
         scf_env%qs_ot_restart%channel(ichannel)%diis_last_response_blend = real_payload(21)
         scf_env%qs_ot_restart%channel(ichannel)%diis_reference_rotation_gnorm = real_payload(22)
         scf_env%qs_ot_restart%channel(ichannel)%diis_raw_residual_ratio = real_payload(23)
         scf_env%qs_ot_restart%channel(ichannel)%diis_response_residual_ratio = real_payload(24)
         scf_env%qs_ot_restart%channel(ichannel)%diis_last_response_model_work = real_payload(25)
         scf_env%qs_ot_restart%channel(ichannel)%density_quality_samples = schur_int_payload(1)
         scf_env%qs_ot_restart%channel(ichannel)%density_audit_samples = schur_int_payload(2)
         scf_env%qs_ot_restart%channel(ichannel)%density_descent_quality_samples = schur_int_payload(3)
         scf_env%qs_ot_restart%channel(ichannel)%density_descent_selected_quality_samples = &
            schur_int_payload(4)
         scf_env%qs_ot_restart%channel(ichannel)%response_selected_quality_samples = &
            schur_int_payload(5)
         scf_env%qs_ot_restart%channel(ichannel)%physical_direct_quality_samples = schur_int_payload(6)
         scf_env%qs_ot_restart%channel(ichannel)%physical_direct_selected_quality_samples = &
            schur_int_payload(7)
         scf_env%qs_ot_restart%channel(ichannel)%phys_direct_density_sel_quality_samples = &
            schur_int_payload(8)
         scf_env%qs_ot_restart%channel(ichannel)%mixed_direct_quality_samples = schur_int_payload(9)
         scf_env%qs_ot_restart%channel(ichannel)%spectral_selected_quality_samples = schur_int_payload(10)
         scf_env%qs_ot_restart%channel(ichannel)%spectral_selected_failure_samples = schur_int_payload(11)
         scf_env%qs_ot_restart%channel(ichannel)%density_descent_support_ready = &
            schur_int_payload(12) /= 0
         scf_env%qs_ot_restart%channel(ichannel)%density_full_selected_quality_samples = &
            schur_int_payload(13)
         scf_env%qs_ot_restart%channel(ichannel)%physical_quality_samples = schur_int_payload(14)
         scf_env%qs_ot_restart%channel(ichannel)%density_full_direct_quality_samples = &
            schur_int_payload(15)
         scf_env%qs_ot_restart%channel(ichannel)%density_fdir_sel_quality_samples = &
            schur_int_payload(16)
         scf_env%qs_ot_restart%channel(ichannel)%density_full_global_quality_samples = &
            schur_int_payload(17)
         scf_env%qs_ot_restart%channel(ichannel)%density_full_global_sel_quality_samples = &
            schur_int_payload(18)
         scf_env%qs_ot_restart%channel(ichannel)%physical_selected_quality_samples = &
            schur_int_payload(19)
         scf_env%qs_ot_restart%channel(ichannel)%density_quality = schur_real_payload(1)
         scf_env%qs_ot_restart%channel(ichannel)%density_quality_ema = schur_real_payload(2)
         scf_env%qs_ot_restart%channel(ichannel)%density_drop_sum = schur_real_payload(3)
         scf_env%qs_ot_restart%channel(ichannel)%density_measured_sum = schur_real_payload(4)
         scf_env%qs_ot_restart%channel(ichannel)%density_lowcurv_ema = schur_real_payload(5)
         scf_env%qs_ot_restart%channel(ichannel)%density_lowcurv_quality_ema = schur_real_payload(6)
         scf_env%qs_ot_restart%channel(ichannel)%density_audit_amplification_ema = &
            schur_real_payload(7)
         scf_env%qs_ot_restart%channel(ichannel)%density_audit_amplification_quality_ema = &
            schur_real_payload(8)
         scf_env%qs_ot_restart%channel(ichannel)%density_audit_lowcurv_ema = schur_real_payload(9)
         scf_env%qs_ot_restart%channel(ichannel)%density_audit_lowcurv_quality_ema = &
            schur_real_payload(10)
         scf_env%qs_ot_restart%channel(ichannel)%density_descent_quality = schur_real_payload(11)
         scf_env%qs_ot_restart%channel(ichannel)%density_descent_quality_ema = schur_real_payload(12)
         scf_env%qs_ot_restart%channel(ichannel)%density_descent_drop_sum = schur_real_payload(13)
         scf_env%qs_ot_restart%channel(ichannel)%density_descent_measured_sum = schur_real_payload(14)
         scf_env%qs_ot_restart%channel(ichannel)%density_descent_coeff_ema = schur_real_payload(15)
         scf_env%qs_ot_restart%channel(ichannel)%density_descent_support = schur_real_payload(16)
         scf_env%qs_ot_restart%channel(ichannel)%density_descent_support_ema = schur_real_payload(17)
         scf_env%qs_ot_restart%channel(ichannel)%density_descent_selected_quality = &
            schur_real_payload(19)
         scf_env%qs_ot_restart%channel(ichannel)%density_descent_selected_quality_ema = &
            schur_real_payload(20)
         scf_env%qs_ot_restart%channel(ichannel)%density_descent_selected_residual_quality = &
            schur_real_payload(21)
         scf_env%qs_ot_restart%channel(ichannel)%density_descent_selected_residual_quality_ema = &
            schur_real_payload(22)
         scf_env%qs_ot_restart%channel(ichannel)%density_descent_selected_drop_sum = &
            schur_real_payload(23)
         scf_env%qs_ot_restart%channel(ichannel)%density_descent_selected_measured_sum = &
            schur_real_payload(24)
         scf_env%qs_ot_restart%channel(ichannel)%physical_direct_quality = schur_real_payload(25)
         scf_env%qs_ot_restart%channel(ichannel)%physical_direct_quality_shadow_ema = &
            schur_real_payload(26)
         scf_env%qs_ot_restart%channel(ichannel)%physical_direct_drop_sum = schur_real_payload(27)
         scf_env%qs_ot_restart%channel(ichannel)%physical_direct_measured_sum = schur_real_payload(28)
         scf_env%qs_ot_restart%channel(ichannel)%physical_direct_ema = schur_real_payload(29)
         scf_env%qs_ot_restart%channel(ichannel)%physical_direct_selected_quality = schur_real_payload(31)
         scf_env%qs_ot_restart%channel(ichannel)%physical_direct_selected_quality_ema = &
            schur_real_payload(32)
         scf_env%qs_ot_restart%channel(ichannel)%physical_direct_selected_residual_quality = &
            schur_real_payload(33)
         scf_env%qs_ot_restart%channel(ichannel)%physical_direct_selected_residual_quality_ema = &
            schur_real_payload(34)
         scf_env%qs_ot_restart%channel(ichannel)%physical_direct_selected_drop_sum = &
            schur_real_payload(35)
         scf_env%qs_ot_restart%channel(ichannel)%physical_direct_selected_measured_sum = &
            schur_real_payload(36)
         scf_env%qs_ot_restart%channel(ichannel)%phys_direct_density_sel_quality = schur_real_payload(38)
         scf_env%qs_ot_restart%channel(ichannel)%phys_direct_density_sel_quality_ema = &
            schur_real_payload(39)
         scf_env%qs_ot_restart%channel(ichannel)%phys_direct_density_sel_resid_quality = &
            schur_real_payload(88)
         scf_env%qs_ot_restart%channel(ichannel)%phys_direct_density_sel_resid_quality_ema = &
            schur_real_payload(89)
         scf_env%qs_ot_restart%channel(ichannel)%phys_direct_density_sel_drop_sum = &
            schur_real_payload(40)
         scf_env%qs_ot_restart%channel(ichannel)%phys_direct_density_sel_measured_sum = &
            schur_real_payload(41)
         scf_env%qs_ot_restart%channel(ichannel)%mixed_direct_quality = schur_real_payload(42)
         scf_env%qs_ot_restart%channel(ichannel)%mixed_direct_quality_ema = schur_real_payload(43)
         scf_env%qs_ot_restart%channel(ichannel)%mixed_direct_drop_sum = schur_real_payload(44)
         scf_env%qs_ot_restart%channel(ichannel)%mixed_direct_measured_sum = schur_real_payload(45)
         scf_env%qs_ot_restart%channel(ichannel)%mixed_direct_weight = schur_real_payload(46)
         scf_env%qs_ot_restart%channel(ichannel)%mixed_direct_weight_ema = schur_real_payload(47)
         scf_env%qs_ot_restart%channel(ichannel)%spectral_selected_quality = schur_real_payload(48)
         scf_env%qs_ot_restart%channel(ichannel)%spectral_selected_quality_ema = schur_real_payload(49)
         scf_env%qs_ot_restart%channel(ichannel)%spectral_selected_drop_sum = schur_real_payload(50)
         scf_env%qs_ot_restart%channel(ichannel)%spectral_selected_measured_sum = schur_real_payload(51)
         scf_env%qs_ot_restart%channel(ichannel)%spectral_selected_residual_quality = &
            schur_real_payload(125)
         scf_env%qs_ot_restart%channel(ichannel)%spectral_selected_residual_quality_ema = &
            schur_real_payload(126)
         scf_env%qs_ot_restart%channel(ichannel)%response_selected_quality = schur_real_payload(53)
         scf_env%qs_ot_restart%channel(ichannel)%response_selected_quality_ema = schur_real_payload(54)
         scf_env%qs_ot_restart%channel(ichannel)%response_selected_residual_quality = &
            schur_real_payload(55)
         scf_env%qs_ot_restart%channel(ichannel)%response_selected_residual_quality_ema = &
            schur_real_payload(56)
         scf_env%qs_ot_restart%channel(ichannel)%response_selected_drop_sum = schur_real_payload(57)
         scf_env%qs_ot_restart%channel(ichannel)%response_selected_measured_sum = schur_real_payload(58)
         scf_env%qs_ot_restart%channel(ichannel)%density_full_selected_quality = schur_real_payload(60)
         scf_env%qs_ot_restart%channel(ichannel)%density_full_selected_quality_ema = &
            schur_real_payload(61)
         scf_env%qs_ot_restart%channel(ichannel)%density_full_selected_residual_quality = &
            schur_real_payload(62)
         scf_env%qs_ot_restart%channel(ichannel)%density_full_selected_residual_quality_ema = &
            schur_real_payload(63)
         scf_env%qs_ot_restart%channel(ichannel)%density_full_selected_drop_sum = schur_real_payload(64)
         scf_env%qs_ot_restart%channel(ichannel)%density_full_selected_measured_sum = &
            schur_real_payload(65)
         scf_env%qs_ot_restart%channel(ichannel)%physical_quality = schur_real_payload(66)
         scf_env%qs_ot_restart%channel(ichannel)%physical_quality_ema = schur_real_payload(67)
         scf_env%qs_ot_restart%channel(ichannel)%physical_drop_sum = schur_real_payload(68)
         scf_env%qs_ot_restart%channel(ichannel)%physical_measured_sum = schur_real_payload(69)
         scf_env%qs_ot_restart%channel(ichannel)%physical_selected_quality = &
            schur_real_payload(134)
         scf_env%qs_ot_restart%channel(ichannel)%physical_selected_quality_ema = &
            schur_real_payload(135)
         scf_env%qs_ot_restart%channel(ichannel)%physical_selected_residual_quality = &
            schur_real_payload(136)
         scf_env%qs_ot_restart%channel(ichannel)%physical_selected_residual_quality_ema = &
            schur_real_payload(137)
         scf_env%qs_ot_restart%channel(ichannel)%physical_selected_drop_sum = &
            schur_real_payload(138)
         scf_env%qs_ot_restart%channel(ichannel)%physical_selected_measured_sum = &
            schur_real_payload(139)
         scf_env%qs_ot_restart%channel(ichannel)%physical_trace_actual_ema = schur_real_payload(70)
         scf_env%qs_ot_restart%channel(ichannel)%physical_trace_raw_ema = schur_real_payload(71)
         scf_env%qs_ot_restart%channel(ichannel)%physical_trace_active_ema = schur_real_payload(72)
         scf_env%qs_ot_restart%channel(ichannel)%physical_trace_actual_diag_ema = &
            schur_real_payload(73)
         scf_env%qs_ot_restart%channel(ichannel)%physical_trace_actual_free_ema = &
            schur_real_payload(74)
         scf_env%qs_ot_restart%channel(ichannel)%physical_trace_raw_diag_ema = &
            schur_real_payload(75)
         scf_env%qs_ot_restart%channel(ichannel)%physical_trace_raw_free_ema = &
            schur_real_payload(76)
         scf_env%qs_ot_restart%channel(ichannel)%physical_trace_active_diag_ema = &
            schur_real_payload(77)
         scf_env%qs_ot_restart%channel(ichannel)%physical_trace_active_free_ema = &
            schur_real_payload(78)
         scf_env%qs_ot_restart%channel(ichannel)%physical_trace_actual_quality_ema = &
            schur_real_payload(79)
         scf_env%qs_ot_restart%channel(ichannel)%physical_trace_raw_quality_ema = &
            schur_real_payload(80)
         scf_env%qs_ot_restart%channel(ichannel)%physical_trace_active_quality_ema = &
            schur_real_payload(81)
         scf_env%qs_ot_restart%channel(ichannel)%physical_trace_actual_diag_quality_ema = &
            schur_real_payload(82)
         scf_env%qs_ot_restart%channel(ichannel)%physical_trace_actual_free_quality_ema = &
            schur_real_payload(83)
         scf_env%qs_ot_restart%channel(ichannel)%physical_trace_raw_diag_quality_ema = &
            schur_real_payload(84)
         scf_env%qs_ot_restart%channel(ichannel)%physical_trace_raw_free_quality_ema = &
            schur_real_payload(85)
         scf_env%qs_ot_restart%channel(ichannel)%physical_trace_active_diag_quality_ema = &
            schur_real_payload(86)
         scf_env%qs_ot_restart%channel(ichannel)%physical_trace_active_free_quality_ema = &
            schur_real_payload(87)
         scf_env%qs_ot_restart%channel(ichannel)%density_full_direct_quality = schur_real_payload(91)
         scf_env%qs_ot_restart%channel(ichannel)%density_full_direct_quality_ema = &
            schur_real_payload(92)
         scf_env%qs_ot_restart%channel(ichannel)%density_full_direct_drop_sum = &
            schur_real_payload(93)
         scf_env%qs_ot_restart%channel(ichannel)%density_full_direct_measured_sum = &
            schur_real_payload(94)
         scf_env%qs_ot_restart%channel(ichannel)%density_full_direct_trust_ema = &
            schur_real_payload(96)
         scf_env%qs_ot_restart%channel(ichannel)%density_full_direct_weighted_trust_ema = &
            schur_real_payload(98)
         scf_env%qs_ot_restart%channel(ichannel)%density_full_direct_abs_trust_ema = &
            schur_real_payload(100)
         scf_env%qs_ot_restart%channel(ichannel)%density_fdir_sel_quality = schur_real_payload(102)
         scf_env%qs_ot_restart%channel(ichannel)%density_fdir_sel_quality_ema = &
            schur_real_payload(103)
         scf_env%qs_ot_restart%channel(ichannel)%density_fdir_sel_resid_quality = &
            schur_real_payload(104)
         scf_env%qs_ot_restart%channel(ichannel)%density_fdir_sel_resid_quality_ema = &
            schur_real_payload(105)
         scf_env%qs_ot_restart%channel(ichannel)%density_fdir_sel_drop_sum = &
            schur_real_payload(106)
         scf_env%qs_ot_restart%channel(ichannel)%density_fdir_sel_measured_sum = &
            schur_real_payload(107)
         scf_env%qs_ot_restart%channel(ichannel)%density_full_global_quality = schur_real_payload(109)
         scf_env%qs_ot_restart%channel(ichannel)%density_full_global_quality_ema = &
            schur_real_payload(110)
         scf_env%qs_ot_restart%channel(ichannel)%density_full_global_drop_sum = &
            schur_real_payload(111)
         scf_env%qs_ot_restart%channel(ichannel)%density_full_global_measured_sum = &
            schur_real_payload(112)
         scf_env%qs_ot_restart%channel(ichannel)%density_full_global_trust_ema = &
            schur_real_payload(114)
         scf_env%qs_ot_restart%channel(ichannel)%density_full_global_weighted_trust_ema = &
            schur_real_payload(116)
         scf_env%qs_ot_restart%channel(ichannel)%density_full_global_abs_trust_ema = &
            schur_real_payload(118)
         scf_env%qs_ot_restart%channel(ichannel)%physical_trace_invariant_ema = &
            schur_real_payload(119)
         scf_env%qs_ot_restart%channel(ichannel)%physical_trace_invariant_diag_ema = &
            schur_real_payload(120)
         scf_env%qs_ot_restart%channel(ichannel)%physical_trace_invariant_free_ema = &
            schur_real_payload(121)
         scf_env%qs_ot_restart%channel(ichannel)%physical_trace_invariant_quality_ema = &
            schur_real_payload(122)
         scf_env%qs_ot_restart%channel(ichannel)%physical_trace_invariant_diag_quality_ema = &
            schur_real_payload(123)
         scf_env%qs_ot_restart%channel(ichannel)%physical_trace_invariant_free_quality_ema = &
            schur_real_payload(124)
         scf_env%qs_ot_restart%channel(ichannel)%density_full_global_selected = &
            schur_real_payload(127)
         scf_env%qs_ot_restart%channel(ichannel)%density_full_global_sel_quality = &
            schur_real_payload(128)
         scf_env%qs_ot_restart%channel(ichannel)%density_full_global_sel_quality_ema = &
            schur_real_payload(129)
         scf_env%qs_ot_restart%channel(ichannel)%density_full_global_sel_residual_quality = &
            schur_real_payload(130)
         scf_env%qs_ot_restart%channel(ichannel)%density_full_global_sel_residual_quality_ema = &
            schur_real_payload(131)
         scf_env%qs_ot_restart%channel(ichannel)%density_full_global_sel_drop_sum = &
            schur_real_payload(132)
         scf_env%qs_ot_restart%channel(ichannel)%density_full_global_sel_measured_sum = &
            schur_real_payload(133)

         IF (scf_env%qs_ot_restart%channel(ichannel)%has_rot_mat_x) THEN
            ALLOCATE (scf_env%qs_ot_restart%channel(ichannel)%rot_mat_x)
            CALL read_restart_fm_matrix(fmat, scf_env%qs_ot_restart%channel(ichannel)%rot_mat_x, rst_unit, &
                                        scf_env%qs_ot_restart%channel(ichannel)%rot_mat_x_dense)
         END IF
         IF (scf_env%qs_ot_restart%channel(ichannel)%has_rot_mat_x_im) THEN
            ALLOCATE (scf_env%qs_ot_restart%channel(ichannel)%rot_mat_x_im)
            CALL read_restart_fm_matrix(fmat, scf_env%qs_ot_restart%channel(ichannel)%rot_mat_x_im, rst_unit, &
                                        scf_env%qs_ot_restart%channel(ichannel)%rot_mat_x_im_dense)
         END IF
         IF (scf_env%qs_ot_restart%channel(ichannel)%has_matrix_x) THEN
            ALLOCATE (scf_env%qs_ot_restart%channel(ichannel)%matrix_x)
            CALL read_restart_fm_matrix(fmat, scf_env%qs_ot_restart%channel(ichannel)%matrix_x, rst_unit, &
                                        scf_env%qs_ot_restart%channel(ichannel)%matrix_x_dense)
         END IF
         IF (scf_env%qs_ot_restart%channel(ichannel)%has_matrix_x_im) THEN
            ALLOCATE (scf_env%qs_ot_restart%channel(ichannel)%matrix_x_im)
            CALL read_restart_fm_matrix(fmat, scf_env%qs_ot_restart%channel(ichannel)%matrix_x_im, rst_unit, &
                                        scf_env%qs_ot_restart%channel(ichannel)%matrix_x_im_dense)
         END IF
         IF (scf_env%qs_ot_restart%channel(ichannel)%has_ener_x) THEN
            ALLOCATE (scf_env%qs_ot_restart%channel(ichannel)%ener_x(nener))
            IF (para_env%is_source()) READ (rst_unit) scf_env%qs_ot_restart%channel(ichannel)%ener_x
            CALL para_env%bcast(scf_env%qs_ot_restart%channel(ichannel)%ener_x)
         END IF
         IF (scf_env%qs_ot_restart%channel(ichannel)%has_matrix_dx) THEN
            ALLOCATE (scf_env%qs_ot_restart%channel(ichannel)%matrix_dx)
            CALL read_restart_fm_matrix(fmat, scf_env%qs_ot_restart%channel(ichannel)%matrix_dx, rst_unit, &
                                        scf_env%qs_ot_restart%channel(ichannel)%matrix_dx_dense)
         END IF
         IF (scf_env%qs_ot_restart%channel(ichannel)%has_matrix_dx_im) THEN
            ALLOCATE (scf_env%qs_ot_restart%channel(ichannel)%matrix_dx_im)
            CALL read_restart_fm_matrix(fmat, scf_env%qs_ot_restart%channel(ichannel)%matrix_dx_im, rst_unit, &
                                        scf_env%qs_ot_restart%channel(ichannel)%matrix_dx_im_dense)
         END IF
         IF (scf_env%qs_ot_restart%channel(ichannel)%has_matrix_gx_old) THEN
            ALLOCATE (scf_env%qs_ot_restart%channel(ichannel)%matrix_gx_old)
            CALL read_restart_fm_matrix(fmat, scf_env%qs_ot_restart%channel(ichannel)%matrix_gx_old, rst_unit, &
                                        scf_env%qs_ot_restart%channel(ichannel)%matrix_gx_old_dense)
         END IF
         IF (scf_env%qs_ot_restart%channel(ichannel)%has_matrix_gx_old_im) THEN
            ALLOCATE (scf_env%qs_ot_restart%channel(ichannel)%matrix_gx_old_im)
            CALL read_restart_fm_matrix(fmat, scf_env%qs_ot_restart%channel(ichannel)%matrix_gx_old_im, rst_unit, &
                                        scf_env%qs_ot_restart%channel(ichannel)%matrix_gx_old_im_dense)
         END IF
         IF (scf_env%qs_ot_restart%channel(ichannel)%has_rot_mat_dx) THEN
            ALLOCATE (scf_env%qs_ot_restart%channel(ichannel)%rot_mat_dx)
            CALL read_restart_fm_matrix(fmat, scf_env%qs_ot_restart%channel(ichannel)%rot_mat_dx, rst_unit, &
                                        scf_env%qs_ot_restart%channel(ichannel)%rot_mat_dx_dense)
         END IF
         IF (scf_env%qs_ot_restart%channel(ichannel)%has_rot_mat_dx_im) THEN
            ALLOCATE (scf_env%qs_ot_restart%channel(ichannel)%rot_mat_dx_im)
            CALL read_restart_fm_matrix(fmat, scf_env%qs_ot_restart%channel(ichannel)%rot_mat_dx_im, rst_unit, &
                                        scf_env%qs_ot_restart%channel(ichannel)%rot_mat_dx_im_dense)
         END IF
         IF (scf_env%qs_ot_restart%channel(ichannel)%has_rot_mat_gx_old) THEN
            ALLOCATE (scf_env%qs_ot_restart%channel(ichannel)%rot_mat_gx_old)
            CALL read_restart_fm_matrix(fmat, scf_env%qs_ot_restart%channel(ichannel)%rot_mat_gx_old, rst_unit, &
                                        scf_env%qs_ot_restart%channel(ichannel)%rot_mat_gx_old_dense)
         END IF
         IF (scf_env%qs_ot_restart%channel(ichannel)%has_rot_mat_gx_old_im) THEN
            ALLOCATE (scf_env%qs_ot_restart%channel(ichannel)%rot_mat_gx_old_im)
            CALL read_restart_fm_matrix(fmat, scf_env%qs_ot_restart%channel(ichannel)%rot_mat_gx_old_im, rst_unit, &
                                        scf_env%qs_ot_restart%channel(ichannel)%rot_mat_gx_old_im_dense)
         END IF
         IF (scf_env%qs_ot_restart%channel(ichannel)%has_ener_dx) THEN
            ALLOCATE (scf_env%qs_ot_restart%channel(ichannel)%ener_dx(nener))
            IF (para_env%is_source()) READ (rst_unit) scf_env%qs_ot_restart%channel(ichannel)%ener_dx
            CALL para_env%bcast(scf_env%qs_ot_restart%channel(ichannel)%ener_dx)
         END IF
         IF (scf_env%qs_ot_restart%channel(ichannel)%has_ener_gx_old) THEN
            ALLOCATE (scf_env%qs_ot_restart%channel(ichannel)%ener_gx_old(nener))
            IF (para_env%is_source()) READ (rst_unit) scf_env%qs_ot_restart%channel(ichannel)%ener_gx_old
            CALL para_env%bcast(scf_env%qs_ot_restart%channel(ichannel)%ener_gx_old)
         END IF
         IF (scf_env%qs_ot_restart%channel(ichannel)%has_diis_scalars) THEN
            ALLOCATE (scf_env%qs_ot_restart%channel(ichannel)%ls_diis(ndiis + 1, ndiis + 1))
            ALLOCATE (scf_env%qs_ot_restart%channel(ichannel)%lss_diis(ndiis + 1, ndiis + 1))
            ALLOCATE (scf_env%qs_ot_restart%channel(ichannel)%c_diis(ndiis + 1))
            ALLOCATE (scf_env%qs_ot_restart%channel(ichannel)%c_broy(ndiis))
            ALLOCATE (scf_env%qs_ot_restart%channel(ichannel)%energy_h(ndiis))
            ALLOCATE (scf_env%qs_ot_restart%channel(ichannel)%ipivot(ndiis + 1))
            IF (para_env%is_source()) THEN
               READ (rst_unit) scf_env%qs_ot_restart%channel(ichannel)%ls_diis
               READ (rst_unit) scf_env%qs_ot_restart%channel(ichannel)%lss_diis
               READ (rst_unit) scf_env%qs_ot_restart%channel(ichannel)%c_diis
               READ (rst_unit) scf_env%qs_ot_restart%channel(ichannel)%c_broy
               READ (rst_unit) scf_env%qs_ot_restart%channel(ichannel)%energy_h
               READ (rst_unit) scf_env%qs_ot_restart%channel(ichannel)%ipivot
            END IF
            CALL para_env%bcast(scf_env%qs_ot_restart%channel(ichannel)%ls_diis)
            CALL para_env%bcast(scf_env%qs_ot_restart%channel(ichannel)%lss_diis)
            CALL para_env%bcast(scf_env%qs_ot_restart%channel(ichannel)%c_diis)
            CALL para_env%bcast(scf_env%qs_ot_restart%channel(ichannel)%c_broy)
            CALL para_env%bcast(scf_env%qs_ot_restart%channel(ichannel)%energy_h)
            CALL para_env%bcast(scf_env%qs_ot_restart%channel(ichannel)%ipivot)
         END IF
         IF (scf_env%qs_ot_restart%channel(ichannel)%has_matrix_history) THEN
            ALLOCATE (scf_env%qs_ot_restart%channel(ichannel)%matrix_h_e(ndiis))
            ALLOCATE (scf_env%qs_ot_restart%channel(ichannel)%matrix_h_x(ndiis))
            ALLOCATE (scf_env%qs_ot_restart%channel(ichannel)%matrix_h_e_dense(ndiis))
            ALLOCATE (scf_env%qs_ot_restart%channel(ichannel)%matrix_h_x_dense(ndiis))
            DO ihist = 1, ndiis
               CALL read_restart_fm_matrix(fmat, scf_env%qs_ot_restart%channel(ichannel)%matrix_h_e(ihist), &
                                           rst_unit, &
                                           scf_env%qs_ot_restart%channel(ichannel)%matrix_h_e_dense(ihist))
               CALL read_restart_fm_matrix(fmat, scf_env%qs_ot_restart%channel(ichannel)%matrix_h_x(ihist), &
                                           rst_unit, &
                                           scf_env%qs_ot_restart%channel(ichannel)%matrix_h_x_dense(ihist))
            END DO
         END IF
         IF (scf_env%qs_ot_restart%channel(ichannel)%has_matrix_history_im) THEN
            ALLOCATE (scf_env%qs_ot_restart%channel(ichannel)%matrix_h_e_im(ndiis))
            ALLOCATE (scf_env%qs_ot_restart%channel(ichannel)%matrix_h_x_im(ndiis))
            ALLOCATE (scf_env%qs_ot_restart%channel(ichannel)%matrix_h_e_im_dense(ndiis))
            ALLOCATE (scf_env%qs_ot_restart%channel(ichannel)%matrix_h_x_im_dense(ndiis))
            DO ihist = 1, ndiis
               CALL read_restart_fm_matrix(fmat, scf_env%qs_ot_restart%channel(ichannel)%matrix_h_e_im(ihist), &
                                           rst_unit, &
                                           scf_env%qs_ot_restart%channel(ichannel)%matrix_h_e_im_dense(ihist))
               CALL read_restart_fm_matrix(fmat, scf_env%qs_ot_restart%channel(ichannel)%matrix_h_x_im(ihist), &
                                           rst_unit, &
                                           scf_env%qs_ot_restart%channel(ichannel)%matrix_h_x_im_dense(ihist))
            END DO
         END IF
         IF (scf_env%qs_ot_restart%channel(ichannel)%has_rot_history) THEN
            ALLOCATE (scf_env%qs_ot_restart%channel(ichannel)%rot_mat_h_e(ndiis))
            ALLOCATE (scf_env%qs_ot_restart%channel(ichannel)%rot_mat_h_x(ndiis))
            ALLOCATE (scf_env%qs_ot_restart%channel(ichannel)%rot_mat_h_e_dense(ndiis))
            ALLOCATE (scf_env%qs_ot_restart%channel(ichannel)%rot_mat_h_x_dense(ndiis))
            DO ihist = 1, ndiis
               CALL read_restart_fm_matrix(fmat, scf_env%qs_ot_restart%channel(ichannel)%rot_mat_h_e(ihist), &
                                           rst_unit, &
                                           scf_env%qs_ot_restart%channel(ichannel)%rot_mat_h_e_dense(ihist))
               CALL read_restart_fm_matrix(fmat, scf_env%qs_ot_restart%channel(ichannel)%rot_mat_h_x(ihist), &
                                           rst_unit, &
                                           scf_env%qs_ot_restart%channel(ichannel)%rot_mat_h_x_dense(ihist))
            END DO
         END IF
         IF (scf_env%qs_ot_restart%channel(ichannel)%has_rot_history_im) THEN
            ALLOCATE (scf_env%qs_ot_restart%channel(ichannel)%rot_mat_h_e_im(ndiis))
            ALLOCATE (scf_env%qs_ot_restart%channel(ichannel)%rot_mat_h_x_im(ndiis))
            ALLOCATE (scf_env%qs_ot_restart%channel(ichannel)%rot_mat_h_e_im_dense(ndiis))
            ALLOCATE (scf_env%qs_ot_restart%channel(ichannel)%rot_mat_h_x_im_dense(ndiis))
            DO ihist = 1, ndiis
               CALL read_restart_fm_matrix(fmat, scf_env%qs_ot_restart%channel(ichannel)%rot_mat_h_e_im(ihist), &
                                           rst_unit, &
                                           scf_env%qs_ot_restart%channel(ichannel)%rot_mat_h_e_im_dense(ihist))
               CALL read_restart_fm_matrix(fmat, scf_env%qs_ot_restart%channel(ichannel)%rot_mat_h_x_im(ihist), &
                                           rst_unit, &
                                           scf_env%qs_ot_restart%channel(ichannel)%rot_mat_h_x_im_dense(ihist))
            END DO
         END IF
         IF (scf_env%qs_ot_restart%channel(ichannel)%has_ener_history) THEN
            ALLOCATE (scf_env%qs_ot_restart%channel(ichannel)%ener_h_e(ndiis, nener))
            ALLOCATE (scf_env%qs_ot_restart%channel(ichannel)%ener_h_x(ndiis, nener))
            IF (para_env%is_source()) THEN
               READ (rst_unit) scf_env%qs_ot_restart%channel(ichannel)%ener_h_e
               READ (rst_unit) scf_env%qs_ot_restart%channel(ichannel)%ener_h_x
            END IF
            CALL para_env%bcast(scf_env%qs_ot_restart%channel(ichannel)%ener_h_e)
            CALL para_env%bcast(scf_env%qs_ot_restart%channel(ichannel)%ener_h_x)
         END IF
         IF (scf_env%qs_ot_restart%channel(ichannel)%has_lbfgs_scalars) THEN
            CPASSERT(ndiis > 1)
            ALLOCATE (scf_env%qs_ot_restart%channel(ichannel)%lbfgs_rho(ndiis - 1))
            ALLOCATE (scf_env%qs_ot_restart%channel(ichannel)%lbfgs_yy(ndiis - 1))
            IF (para_env%is_source()) THEN
               READ (rst_unit) scf_env%qs_ot_restart%channel(ichannel)%lbfgs_rho
               READ (rst_unit) scf_env%qs_ot_restart%channel(ichannel)%lbfgs_yy
            END IF
            CALL para_env%bcast(scf_env%qs_ot_restart%channel(ichannel)%lbfgs_rho)
            CALL para_env%bcast(scf_env%qs_ot_restart%channel(ichannel)%lbfgs_yy)
         END IF
      END DO

      IF (ANY(.NOT. scf_env%qs_ot_restart%channel(:)%valid)) THEN
         scf_env%qs_ot_restart%valid = .FALSE.
         CALL cp_warn(__LOCATION__, &
                      "Ignoring incomplete k-point OT restart payload; density restart remains available")
      END IF

   END SUBROUTINE read_kpoint_ot_restart_payload

! **************************************************************************************************
!> \brief read a restart FM matrix using an existing FM distribution as template
!> \param fmat_template FM template
!> \param fmat matrix read from restart
!> \param rst_unit restart unit
!> \param dense optional replicated matrix cache for redistribution to k-point groups
! **************************************************************************************************
   SUBROUTINE read_restart_fm_matrix(fmat_template, fmat, rst_unit, dense)
      TYPE(cp_fm_type), INTENT(IN)                       :: fmat_template
      TYPE(cp_fm_type), INTENT(INOUT)                    :: fmat
      INTEGER, INTENT(IN)                                :: rst_unit
      TYPE(qs_ot_restart_dense_matrix_type), &
         INTENT(INOUT), OPTIONAL                         :: dense

      INTEGER                                            :: ncol, nrow
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct

      NULLIFY (fm_struct)
      nrow = 0
      ncol = 0
      IF (fmat_template%matrix_struct%para_env%is_source()) READ (rst_unit) nrow, ncol
      CALL fmat_template%matrix_struct%para_env%bcast(nrow)
      CALL fmat_template%matrix_struct%para_env%bcast(ncol)
      CALL cp_fm_struct_create(fm_struct, template_fmstruct=fmat_template%matrix_struct, &
                               nrow_global=nrow, ncol_global=ncol)
      CALL cp_fm_create(fmat, fm_struct)
      CALL cp_fm_read_unformatted(fmat, rst_unit)
      IF (PRESENT(dense)) THEN
         IF (ALLOCATED(dense%data)) DEALLOCATE (dense%data)
         ALLOCATE (dense%data(nrow, ncol))
         CALL cp_fm_get_submatrix(fmat, dense%data)
      END IF
      CALL cp_fm_struct_release(fm_struct)

   END SUBROUTINE read_restart_fm_matrix

! **************************************************************************************************
!> \brief ...
!> \param cell ...
!> \param cell_to_index ...
!> \return ...
! **************************************************************************************************
   FUNCTION get_index(cell, cell_to_index) RESULT(index)
      INTEGER, DIMENSION(3), INTENT(IN)                  :: cell
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      INTEGER                                            :: index

      IF (CELL(1) >= LBOUND(cell_to_index, 1) .AND. CELL(1) <= UBOUND(cell_to_index, 1) .AND. &
          CELL(2) >= LBOUND(cell_to_index, 2) .AND. CELL(2) <= UBOUND(cell_to_index, 2) .AND. &
          CELL(3) >= LBOUND(cell_to_index, 3) .AND. CELL(3) <= UBOUND(cell_to_index, 3)) THEN

         index = cell_to_index(cell(1), cell(2), cell(3))

      ELSE

         index = 0

      END IF

   END FUNCTION get_index

END MODULE kpoint_io
