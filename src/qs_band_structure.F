!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculation of band structures
!> \par History
!>       2015.06 created [JGH]
!> \author JGH
! **************************************************************************************************
MODULE qs_band_structure
   USE cell_types,                      ONLY: cell_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE cp_log_handling,                 ONLY: cp_logger_get_default_io_unit
   USE dft_plus_u,                      ONLY: lowdin_kp_prepare
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE kpoint_methods,                  ONLY: kpoint_env_initialize,&
                                              kpoint_init_cell_index,&
                                              kpoint_initialize_mo_set,&
                                              kpoint_initialize_mos
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_create,&
                                              kpoint_env_type,&
                                              kpoint_release,&
                                              kpoint_star_type,&
                                              kpoint_sym_create,&
                                              kpoint_type,&
                                              read_band_kpoint_set
   USE machine,                         ONLY: m_walltime
   USE message_passing,                 ONLY: mp_para_env_type
   USE physcon,                         ONLY: evolt
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_env_release,&
                                              qs_environment_type
   USE qs_gamma2kp,                     ONLY: create_kp_from_gamma
   USE qs_kpoint_operators,             ONLY: kpoint_operator_copy_native,&
                                              kpoint_slot_ks,&
                                              kpoint_slot_s
   USE qs_ks_types,                     ONLY: get_ks_env
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              set_mo_set
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_scf_diagonalization,          ONLY: diag_kp_smat,&
                                              do_general_diag_kp
   USE qs_scf_types,                    ONLY: qs_scf_env_type
   USE qs_upf_inverse_scf,              ONLY: qs_upf_inverse_scf_active
   USE qs_upf_inverse_types,            ONLY: qs_upf_inverse_scf_type
   USE scf_control_types,               ONLY: scf_control_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_band_structure'

   PUBLIC :: calculate_band_structure, calculate_kp_orbitals, calculate_kpoints_for_bs

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief Main routine for band structure calculation
!> \param qs_env ...
!> \author JGH
! **************************************************************************************************
   SUBROUTINE calculate_band_structure(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      LOGICAL                                            :: do_kpoints, explicit
      TYPE(section_vals_type), POINTER                   :: bs_input

      bs_input => section_vals_get_subs_vals(qs_env%input, "DFT%PRINT%BAND_STRUCTURE")
      CALL section_vals_get(bs_input, explicit=explicit)
      IF (explicit) THEN
         CALL get_qs_env(qs_env, do_kpoints=do_kpoints)
         IF (do_kpoints) THEN
            CALL do_calculate_band_structure(qs_env)
         ELSE
            BLOCK
               TYPE(qs_environment_type), POINTER :: qs_env_kp
               IF (qs_upf_inverse_scf_active(qs_env%ks_env)) THEN
                  CPABORT("Band structure requires SCF k points for original USPP or PAW.")
               END IF
               CALL create_kp_from_gamma(qs_env, qs_env_kp)
               CALL do_calculate_band_structure(qs_env_kp)
               CALL qs_env_release(qs_env_kp)
               DEALLOCATE (qs_env_kp)
            END BLOCK
         END IF
      END IF

   END SUBROUTINE calculate_band_structure

! **************************************************************************************************
!> \brief band structure calculation
!> \param qs_env ...
!> \author JGH
! **************************************************************************************************
   SUBROUTINE do_calculate_band_structure(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=default_string_length)               :: filename
      CHARACTER(LEN=default_string_length), ALLOCATABLE  :: spname(:)
      INTEGER                                            :: bs_data_unit, i_rep, ik, ikpgr, imo, ip, &
                                                            ispin, n_ptr, n_rep, nadd, nkp, nmo, &
                                                            npoints, nspins, unit_nr
      INTEGER, DIMENSION(2)                              :: kp_range
      LOGICAL                                            :: explicit, io_default, my_kpgrp
      REAL(dp), ALLOCATABLE                              :: kpgeneral(:, :), kspecial(:, :)
      REAL(dp), POINTER                                  :: xkp(:, :)
      REAL(KIND=dp)                                      :: t1, t2
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvalues, eigval, occnum, &
                                                            occupation_numbers, wkp
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_env_type), POINTER                     :: kp
      TYPE(kpoint_type), POINTER                         :: kpoint
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(section_vals_type), POINTER                   :: bs_input, kpset

      bs_input => section_vals_get_subs_vals(qs_env%input, "DFT%PRINT%BAND_STRUCTURE")
      CALL section_vals_get(bs_input, explicit=explicit)
      IF (explicit) THEN
         CALL section_vals_val_get(bs_input, "FILE_NAME", c_val=filename)
         CALL section_vals_val_get(bs_input, "ADDED_MOS", i_val=nadd)
         unit_nr = cp_logger_get_default_io_unit()
         CALL get_qs_env(qs_env=qs_env, para_env=para_env)
         CALL get_qs_env(qs_env, cell=cell)
         kpset => section_vals_get_subs_vals(bs_input, "KPOINT_SET")
         CALL section_vals_get(kpset, n_repetition=n_rep)
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, FMT="(/,T2,A)") "KPOINTS| Band Structure Calculation"
            WRITE (unit_nr, FMT="(T2,A,T71,I10)") "KPOINTS| Number of k-point sets", n_rep
            IF (nadd /= 0) THEN
               WRITE (unit_nr, FMT="(T2,A,T71,I10)") "KPOINTS| Number of added MOs/bands", nadd
            END IF
         END IF
         IF (filename == "") THEN
            ! use standard output file
            bs_data_unit = unit_nr
            io_default = .TRUE.
         ELSE
            io_default = .FALSE.
            IF (para_env%is_source()) THEN
               CALL open_file(filename, unit_number=bs_data_unit, file_status="UNKNOWN", file_action="WRITE", &
                              file_position="APPEND")
            ELSE
               bs_data_unit = -1
            END IF
         END IF
         DO i_rep = 1, n_rep
            t1 = m_walltime()
            CALL read_band_kpoint_set(kpset, i_rep, cell, kspecial, spname, kpgeneral)
            n_ptr = SIZE(kspecial, 2)
            npoints = SIZE(kpgeneral, 2)
            NULLIFY (kpoint)
            CALL calculate_kp_orbitals(qs_env, kpoint, "GENERAL", nadd, kpgeneral=kpgeneral)
            DEALLOCATE (kpgeneral)

            CALL get_qs_env(qs_env, dft_control=dft_control)
            nspins = dft_control%nspins
            kp => kpoint%kp_env(1)%kpoint_env
            CALL get_mo_set(kp%mos(1), nmo=nmo)
            ALLOCATE (eigval(nmo), occnum(nmo))
            CALL get_kpoint_info(kpoint, nkp=nkp, kp_range=kp_range, xkp=xkp, wkp=wkp)

            IF (unit_nr > 0) THEN
               WRITE (UNIT=unit_nr, FMT="(T2,A,I4,T71,I10)") &
                  "KPOINTS| Number of k-points in set ", i_rep, npoints
               WRITE (UNIT=unit_nr, FMT="(T2,A)") &
                  "KPOINTS| In units of b-vector [2pi/Bohr]"
               DO ip = 1, n_ptr
                  WRITE (UNIT=unit_nr, FMT="(T2,A,I5,1X,A11,3(1X,F12.6))") &
                     "KPOINTS| Special point ", ip, ADJUSTL(TRIM(spname(ip))), kspecial(1:3, ip)
               END DO
            END IF
            IF (bs_data_unit > 0 .AND. (bs_data_unit /= unit_nr)) THEN
               WRITE (UNIT=bs_data_unit, FMT="(4(A,I0),A)") &
                  "# Set ", i_rep, ": ", n_ptr, " special points, ", npoints, " k-points, ", nmo, " bands"
               DO ip = 1, n_ptr
                  WRITE (UNIT=bs_data_unit, FMT="(A,I0,T20,T24,3(1X,F14.8),2X,A)") &
                     "#  Special point ", ip, kspecial(1:3, ip), ADJUSTL(TRIM(spname(ip)))
               END DO
            END IF

            DO ik = 1, nkp
               my_kpgrp = (ik >= kp_range(1) .AND. ik <= kp_range(2))
               DO ispin = 1, nspins
                  IF (my_kpgrp) THEN
                     ikpgr = ik - kp_range(1) + 1
                     kp => kpoint%kp_env(ikpgr)%kpoint_env
                     CALL get_mo_set(kp%mos(ispin), eigenvalues=eigenvalues, occupation_numbers=occupation_numbers)
                     eigval(1:nmo) = eigenvalues(1:nmo)
                     occnum(1:nmo) = occupation_numbers(1:nmo)
                  ELSE
                     eigval(1:nmo) = 0.0_dp
                     occnum(1:nmo) = 0.0_dp
                  END IF
                  CALL kpoint%para_env_inter_kp%sum(eigval)
                  CALL kpoint%para_env_inter_kp%sum(occnum)
                  IF (bs_data_unit > 0) THEN
                     WRITE (UNIT=bs_data_unit, FMT="(A,I0,T15,A,I0,A,T24,3(1X,F14.8),3X,F14.8)") &
                        "#  Point ", ik, "  Spin ", ispin, ":", xkp(1:3, ik), wkp(ik)
                     WRITE (UNIT=bs_data_unit, FMT="(A)") &
                        "#   Band    Energy [eV]     Occupation"
                     DO imo = 1, nmo
                        WRITE (UNIT=bs_data_unit, FMT="(T2,I7,2(1X,F14.8))") &
                           imo, eigval(imo)*evolt, occnum(imo)
                     END DO
                  END IF
               END DO
            END DO

            DEALLOCATE (kspecial, spname)
            DEALLOCATE (eigval, occnum)
            CALL kpoint_release(kpoint)
            t2 = m_walltime()
            IF (unit_nr > 0) THEN
               WRITE (UNIT=unit_nr, FMT="(T2,A,T67,F14.3)") "KPOINTS| Time for k-point line ", t2 - t1
            END IF

         END DO

         ! Close output files
         IF (.NOT. io_default) THEN
            IF (para_env%is_source()) CALL close_file(bs_data_unit)
         END IF

      END IF

   END SUBROUTINE do_calculate_band_structure

! **************************************************************************************************
!> \brief diagonalize KS matrices at a set of kpoints
!> \param qs_env ...
!> \param kpoint ...
!> \param scheme ...
!> \param nadd ...
!> \param mp_grid ...
!> \param kpgeneral ...
!> \param group_size_ext ...
!> \param kp_shift ...
!> \param gamma_centered ...
!> \author JGH
! **************************************************************************************************
   SUBROUTINE calculate_kp_orbitals(qs_env, kpoint, scheme, nadd, mp_grid, kpgeneral, group_size_ext, &
                                    kp_shift, gamma_centered)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(kpoint_type), POINTER                         :: kpoint
      CHARACTER(LEN=*), INTENT(IN)                       :: scheme
      INTEGER, INTENT(IN)                                :: nadd
      INTEGER, DIMENSION(3), INTENT(IN), OPTIONAL        :: mp_grid
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), &
         OPTIONAL                                        :: kpgeneral
      INTEGER, INTENT(IN), OPTIONAL                      :: group_size_ext
      REAL(KIND=dp), DIMENSION(3), INTENT(IN), OPTIONAL  :: kp_shift
      LOGICAL, INTENT(IN), OPTIONAL                      :: gamma_centered

      COMPLEX(dp), POINTER                               :: physical_metric(:, :, :, :)
      INTEGER                                            :: ik, nao, nphysical
      INTEGER, ALLOCATABLE                               :: physical_indices(:), source_indices(:)
      LOGICAL                                            :: diis_step, original_upf
      REAL(dp), POINTER                                  :: potential(:, :)
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks, matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_type), POINTER                         :: scf_kpoints
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(qs_upf_inverse_scf_type), POINTER             :: upf
      TYPE(scf_control_type), POINTER                    :: scf_control

      CALL calculate_kpoints_for_bs(kpoint, scheme, group_size_ext, mp_grid, kpgeneral, &
                                    kp_shift, gamma_centered)

      original_upf = qs_upf_inverse_scf_active(qs_env%ks_env)
      CALL get_qs_env(qs_env, matrix_ks_kp=matrix_ks, matrix_metric_kp=matrix_s, &
                      scf_env=scf_env, scf_control=scf_control, kpoints=scf_kpoints)
      IF (original_upf) THEN
         CPASSERT(ASSOCIATED(scf_kpoints))
         CALL get_ks_env(qs_env%ks_env, upf_inverse=upf)
         CALL match_scf_kpoints(scf_kpoints, kpoint, upf%bloch%star, source_indices, physical_indices)
      END IF

      CALL get_qs_env(qs_env=qs_env, para_env=para_env, blacs_env=blacs_env)
      CALL kpoint_env_initialize(kpoint, para_env, blacs_env)

      CALL kpoint_initialize_mos(kpoint, qs_env%mos, nadd)
      CALL kpoint_initialize_mo_set(kpoint)

      CALL get_qs_env(qs_env, sab_kp=sab_nl, dft_control=dft_control)
      CALL kpoint_init_cell_index(kpoint, sab_nl, para_env, dft_control%nimages)

      IF (original_upf) THEN
         nao = SIZE(upf%metric_k, 1)
         nphysical = SIZE(upf%metric_k, 3)
         physical_metric(1:nao, 1:nao, 1:nphysical, 1:1) => upf%metric_k
         CALL kpoint_operator_copy_native(scf_kpoints, kpoint, source_indices, matrix_s, kpoint_slot_s, &
                                          physical_indices, physical_metric)
         CALL kpoint_operator_copy_native(scf_kpoints, kpoint, source_indices, matrix_ks, kpoint_slot_ks, &
                                          physical_indices, upf%ks_k)
      END IF
      NULLIFY (potential)
      IF (ASSOCIATED(scf_kpoints)) THEN
         IF (ALLOCATED(scf_kpoints%lowdin_v)) THEN
            IF (ANY(scf_kpoints%lowdin_v /= 0.0_dp)) potential => scf_kpoints%lowdin_v
         END IF
      END IF
      IF (ASSOCIATED(potential)) THEN
         CALL diag_kp_smat(matrix_s, kpoint, scf_env%scf_work1)
         DO ik = 1, SIZE(kpoint%kp_env)
            CALL lowdin_kp_prepare(kpoint%kp_env(ik)%kpoint_env, scf_kpoints%lowdin_orbitals, .FALSE.)
         END DO
      END IF
      CALL do_general_diag_kp(matrix_ks, matrix_s, kpoint, scf_env, scf_control, .FALSE., diis_step, &
                              potential=potential)
      IF (original_upf) CALL copy_scf_occupations(scf_kpoints, kpoint, source_indices)

   END SUBROUTINE calculate_kp_orbitals

! **************************************************************************************************
!> \brief Match SCF representatives and physical symmetry partners.
!> \param scf ...
!> \param output ...
!> \param star Explicit k points of the reduced SCF calculation.
!> \param indices SCF occupation indices, signed for time-reversed operators.
!> \param physical_indices Index of an explicit symmetry partner, zero for a representative.
! **************************************************************************************************
   SUBROUTINE match_scf_kpoints(scf, output, star, indices, physical_indices)
      TYPE(kpoint_type), INTENT(IN)                      :: scf, output
      TYPE(kpoint_star_type), INTENT(IN)                 :: star
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: indices(:), physical_indices(:)

      INTEGER                                            :: ik, j

      ALLOCATE (indices(output%nkp), physical_indices(output%nkp), SOURCE=0)
      DO ik = 1, output%nkp
         DO j = 1, scf%nkp
            IF (same_band_point(output%xkp(:, ik), scf%xkp(:, j))) THEN
               indices(ik) = j
               EXIT
            END IF
         END DO
         IF (indices(ik) /= 0) CYCLE
         DO j = 1, SIZE(star%vectors, 2)
            IF (same_band_point(output%xkp(:, ik), star%vectors(:, j))) THEN
               indices(ik) = star%representative(j)
               physical_indices(ik) = j
               EXIT
            END IF
         END DO
         IF (indices(ik) /= 0) CYCLE
         IF (scf%inversion_symmetry_only .AND. .NOT. scf%full_grid) THEN
            DO j = 1, scf%nkp
               IF (same_band_point(output%xkp(:, ik), -scf%xkp(:, j))) THEN
                  indices(ik) = -j
                  EXIT
               END IF
            END DO
         END IF
         IF (indices(ik) == 0) CPABORT("USPP or PAW post-SCF k points must belong to the SCF mesh.")
      END DO
   END SUBROUTINE match_scf_kpoints

! **************************************************************************************************
!> \brief Compare reduced coordinates modulo reciprocal lattice vectors.
!> \param left ...
!> \param right ...
!> \return ...
! **************************************************************************************************
   PURE LOGICAL FUNCTION same_band_point(left, right) RESULT(same)
      REAL(dp), INTENT(IN)                               :: left(3), right(3)

      REAL(dp)                                           :: difference(3)

      difference = left - right
      same = MAXVAL(ABS(difference - ANINT(difference))) <= 128.0_dp*EPSILON(1.0_dp)
   END FUNCTION same_band_point

! **************************************************************************************************
!> \brief Preserve SCF occupations and chemical potentials when sampling its eigenstates.
!> \param scf ...
!> \param output ...
!> \param indices ...
! **************************************************************************************************
   SUBROUTINE copy_scf_occupations(scf, output, indices)
      TYPE(kpoint_type), INTENT(IN)                      :: scf, output
      INTEGER, INTENT(IN)                                :: indices(:)

      INTEGER                                            :: dst, ik, ispin, nmo, src
      REAL(dp), ALLOCATABLE                              :: values(:)
      TYPE(kpoint_env_type), POINTER                     :: kp

      DO ispin = 1, SIZE(output%kp_env(1)%kpoint_env%mos)
         CALL get_mo_set(output%kp_env(1)%kpoint_env%mos(ispin), nmo=nmo)
         ALLOCATE (values(1 + nmo))
         DO ik = 1, SIZE(indices)
            src = ABS(indices(ik)) - scf%kp_range(1) + 1
            dst = ik - output%kp_range(1) + 1
            values = 0.0_dp
            IF (src > 0 .AND. src <= SIZE(scf%kp_env) .AND. scf%para_env_kp%is_source()) THEN
               kp => scf%kp_env(src)%kpoint_env
               nmo = MIN(SIZE(values) - 1, kp%mos(ispin)%nmo)
               values(1:nmo) = kp%mos(ispin)%occupation_numbers(1:nmo)
               values(SIZE(values)) = kp%mos(ispin)%mu
            END IF
            CALL scf%para_env%sum(values)
            IF (dst > 0 .AND. dst <= SIZE(output%kp_env)) THEN
               kp => output%kp_env(dst)%kpoint_env
               kp%mos(ispin)%occupation_numbers(:) = values(1:SIZE(values) - 1)
               CALL set_mo_set(kp%mos(ispin), mu=values(SIZE(values)))
            END IF
         END DO
         DEALLOCATE (values)
      END DO
   END SUBROUTINE copy_scf_occupations

! **************************************************************************************************
!> \brief ...
!> \param kpoint ...
!> \param scheme ...
!> \param group_size_ext ...
!> \param mp_grid ...
!> \param kpgeneral ...
!> \param kp_shift ...
!> \param gamma_centered ...
! **************************************************************************************************
   SUBROUTINE calculate_kpoints_for_bs(kpoint, scheme, group_size_ext, mp_grid, kpgeneral, &
                                       kp_shift, gamma_centered)

      TYPE(kpoint_type), POINTER                         :: kpoint
      CHARACTER(LEN=*), INTENT(IN)                       :: scheme
      INTEGER, INTENT(IN), OPTIONAL                      :: group_size_ext
      INTEGER, DIMENSION(3), INTENT(IN), OPTIONAL        :: mp_grid
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), &
         OPTIONAL                                        :: kpgeneral
      REAL(KIND=dp), DIMENSION(3), INTENT(IN), OPTIONAL  :: kp_shift
      LOGICAL, INTENT(IN), OPTIONAL                      :: gamma_centered

      INTEGER                                            :: i, idim, ix, iy, iz, npoints
      INTEGER, DIMENSION(3)                              :: ik
      LOGICAL                                            :: gamma_mesh
      REAL(KIND=dp), DIMENSION(3)                        :: kpt_latt, shift

      CPASSERT(.NOT. ASSOCIATED(kpoint))

      CALL kpoint_create(kpoint)

      IF (PRESENT(kp_shift)) THEN
         shift = kp_shift
      ELSE
         shift = 0.0_dp
      END IF
      IF (PRESENT(gamma_centered)) THEN
         gamma_mesh = gamma_centered
      ELSE
         gamma_mesh = .FALSE.
      END IF

      kpoint%kp_scheme = scheme
      kpoint%symmetry = .FALSE.
      kpoint%verbose = .FALSE.
      kpoint%full_grid = .FALSE.
      kpoint%use_real_wfn = .FALSE.
      kpoint%gamma_centered = gamma_mesh
      kpoint%kp_shift = shift
      kpoint%eps_geo = 1.e-6_dp
      IF (PRESENT(group_size_ext)) THEN
         kpoint%parallel_group_size = group_size_ext
      ELSE
         kpoint%parallel_group_size = -1
      END IF
      SELECT CASE (scheme)
      CASE ("GAMMA")
         kpoint%nkp = 1
         ALLOCATE (kpoint%xkp(3, 1), kpoint%wkp(1))
         kpoint%xkp(1:3, 1) = 0.0_dp
         kpoint%wkp(1) = 1.0_dp
         kpoint%symmetry = .TRUE.
         ALLOCATE (kpoint%kp_sym(1))
         NULLIFY (kpoint%kp_sym(1)%kpoint_sym)
         CALL kpoint_sym_create(kpoint%kp_sym(1)%kpoint_sym)
      CASE ("MONKHORST-PACK", "MACDONALD")
         CPASSERT(PRESENT(mp_grid))
         npoints = mp_grid(1)*mp_grid(2)*mp_grid(3)
         kpoint%nkp_grid(1:3) = mp_grid(1:3)
         kpoint%full_grid = .TRUE.
         kpoint%nkp = npoints
         ALLOCATE (kpoint%xkp(3, npoints), kpoint%wkp(npoints))
         kpoint%wkp(:) = 1._dp/REAL(npoints, KIND=dp)
         i = 0
         DO ix = 1, mp_grid(1)
            DO iy = 1, mp_grid(2)
               DO iz = 1, mp_grid(3)
                  i = i + 1
                  ik = [ix, iy, iz]
                  DO idim = 1, 3
                     IF (gamma_mesh .AND. MOD(mp_grid(idim), 2) == 0) THEN
                        kpt_latt(idim) = REAL(2*ik(idim) - mp_grid(idim), KIND=dp)/ &
                                         (2._dp*REAL(mp_grid(idim), KIND=dp))
                     ELSE
                        kpt_latt(idim) = REAL(2*ik(idim) - mp_grid(idim) - 1, KIND=dp)/ &
                                         (2._dp*REAL(mp_grid(idim), KIND=dp))
                     END IF
                  END DO
                  kpoint%xkp(1:3, i) = kpt_latt(1:3) + shift(1:3)
               END DO
            END DO
         END DO
         ! default: no symmetry settings
         ALLOCATE (kpoint%kp_sym(kpoint%nkp))
         DO i = 1, kpoint%nkp
            NULLIFY (kpoint%kp_sym(i)%kpoint_sym)
            CALL kpoint_sym_create(kpoint%kp_sym(i)%kpoint_sym)
         END DO
      CASE ("GENERAL")
         CPASSERT(PRESENT(kpgeneral))
         npoints = SIZE(kpgeneral, 2)
         kpoint%nkp = npoints
         ALLOCATE (kpoint%xkp(3, npoints), kpoint%wkp(npoints))
         kpoint%wkp(:) = 1._dp/REAL(npoints, KIND=dp)
         kpoint%xkp(1:3, 1:npoints) = kpgeneral(1:3, 1:npoints)
         ! default: no symmetry settings
         ALLOCATE (kpoint%kp_sym(kpoint%nkp))
         DO i = 1, kpoint%nkp
            NULLIFY (kpoint%kp_sym(i)%kpoint_sym)
            CALL kpoint_sym_create(kpoint%kp_sym(i)%kpoint_sym)
         END DO
      CASE DEFAULT
         CPABORT("Unknown kpoint scheme requested")
      END SELECT

   END SUBROUTINE calculate_kpoints_for_bs

END MODULE qs_band_structure
