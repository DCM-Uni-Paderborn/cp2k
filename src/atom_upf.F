!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines that process Quantum Espresso UPF files.
!> \par History
!>    * 07.2018 CP2K-SIRIUS interface [Juerg Hutter]
!>    * 02.2016 created [Juerg Hutter]
! **************************************************************************************************
MODULE atom_upf
   USE cp_parser_methods,               ONLY: parser_get_next_line,&
                                              parser_get_object,&
                                              parser_test_next_token
   USE cp_parser_types,                 ONLY: cp_parser_type,&
                                              parser_create,&
                                              parser_release
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp
   USE periodic_table,                  ONLY: get_ptable_info,&
                                              ptable
#include "./base/base_uses.f90"

   IMPLICIT NONE

   ! use same value as in atom_types!
   INTEGER, PARAMETER                                :: lmat = 3

   TYPE upf_paw_type
      INTEGER                                        :: data_format = -1, cutoff_index = -1, lmax_aug = -1
      CHARACTER(LEN=default_string_length)           :: augmentation_shape = ""
      REAL(dp)                                       :: cutoff_radius = -1.0_dp, augmentation_epsilon = 1.e-12_dp
      REAL(dp)                                       :: core_energy = 0.0_dp
      LOGICAL                                        :: has_core_energy = .FALSE.
      ! Original r*phi(r) partial waves, one per beta projector, without normalization.
      REAL(dp), ALLOCATABLE                          :: ae_wfc(:, :), ps_wfc(:, :), occupations(:)
      INTEGER, ALLOCATABLE                           :: angular_momenta(:)
      ! Multipoles have the same beta convention as PP_Q; L starts at zero.
      REAL(dp), ALLOCATABLE                          :: multipoles(:, :, :)
      ! True core density (without 4*pi*r**2), AE local potential and core energy in Hartree.
      REAL(dp), ALLOCATABLE                          :: rho_core(:), vlocal(:)
   END TYPE upf_paw_type

   TYPE atom_upfpot_type
      CHARACTER(LEN=2)                               :: symbol = ""
      CHARACTER(LEN=default_string_length)           :: pname = ""
      INTEGER, DIMENSION(0:lmat)                     :: econf = 0
      REAL(dp)                                       :: zion = 0.0_dp
      CHARACTER(LEN=default_string_length)           :: version = ""
      CHARACTER(LEN=default_path_length)             :: filename = ""
      ! <INFO>
      INTEGER                                        :: maxinfo = 100
      CHARACTER(LEN=default_string_length), DIMENSION(100) &
         :: info = ""
      ! <HEADER>
      CHARACTER(LEN=default_string_length)           :: generated = ""
      CHARACTER(LEN=default_string_length)           :: author = ""
      CHARACTER(LEN=default_string_length)           :: date = ""
      CHARACTER(LEN=default_string_length)           :: comment = ""
      CHARACTER(LEN=4)                               :: pseudo_type = ""
      CHARACTER(LEN=15)                              :: relativistic = ""
      CHARACTER(LEN=default_string_length)           :: functional = ""
      LOGICAL                                        :: is_ultrasoft = .FALSE.
      LOGICAL                                        :: is_paw = .FALSE.
      LOGICAL                                        :: is_coulomb = .FALSE.
      LOGICAL                                        :: has_so = .FALSE.
      LOGICAL                                        :: has_wfc = .FALSE.
      LOGICAL                                        :: has_gipaw = .FALSE.
      LOGICAL                                        :: paw_as_gipaw = .FALSE.
      LOGICAL                                        :: core_correction = .FALSE.
      REAL(dp)                                       :: total_psenergy = 0.0_dp
      REAL(dp)                                       :: wfc_cutoff = 0.0_dp
      REAL(dp)                                       :: rho_cutoff = 0.0_dp
      INTEGER                                        :: l_max = -100
      INTEGER                                        :: l_max_rho = -1
      INTEGER                                        :: l_local = -1
      INTEGER                                        :: mesh_size = -1
      INTEGER                                        :: number_of_wfc = -1
      INTEGER                                        :: number_of_proj = -1
      ! <MESH>
      REAL(dp)                                       :: dx = 0.0_dp
      REAL(dp)                                       :: xmin = 0.0_dp
      REAL(dp)                                       :: rmax = 0.0_dp
      REAL(dp)                                       :: zmesh = 0.0_dp
      REAL(dp), DIMENSION(:), ALLOCATABLE            :: r, rab
      ! <NLCC>
      REAL(dp), DIMENSION(:), ALLOCATABLE            :: rho_nlcc
      ! <LOCAL>
      REAL(dp), DIMENSION(:), ALLOCATABLE            :: vlocal
      ! <NONLOCAL>: keep the tabulated projector amplitude; D includes Ry -> Hartree.
      ! beta stores beta(r), with the radial prefactor in PP_BETA removed.
      REAL(dp), DIMENSION(:, :), ALLOCATABLE         :: dion
      REAL(dp), DIMENSION(:, :), ALLOCATABLE         :: beta
      INTEGER, DIMENSION(:), ALLOCATABLE             :: lbeta
      ! <AUGMENTATION>: PP_Q and r**2*Q_ij^L(r), in the original beta convention.
      ! The pair index is j*(j-1)/2+i for i<=j; the last dimension starts at L=0.
      REAL(dp), DIMENSION(:, :), ALLOCATABLE         :: qij
      REAL(dp), DIMENSION(:, :, :), ALLOCATABLE      :: qijl
      TYPE(upf_paw_type), ALLOCATABLE                 :: paw
      ! <SEMILOCAL>
      REAL(dp), DIMENSION(:, :), ALLOCATABLE         :: vsemi
   END TYPE atom_upfpot_type

   PRIVATE
   PUBLIC  :: atom_read_upf, atom_upfpot_type, atom_release_upf, upf_radial_weights

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'atom_upf'

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief Composite Simpson weights for dr on the supplied UPF mesh index.
!> \param pot Original radial coordinates and dr/d(index).
!> \return Radial weights, with a last trapezoidal interval for an even mesh size.
! **************************************************************************************************
   FUNCTION upf_radial_weights(pot) RESULT(weights)
      TYPE(atom_upfpot_type), INTENT(IN)                 :: pot
      REAL(dp)                                           :: weights(SIZE(pot%r))

      INTEGER                                            :: n, ns

      n = SIZE(pot%r)
      CPASSERT(n >= 3 .AND. SIZE(pot%rab) == n)
      CPASSERT(ALL(pot%r(2:) > pot%r(:n - 1)) .AND. ALL(pot%rab > 0.0_dp))
      ns = n - 1 + MODULO(n, 2)
      weights = 0.0_dp
      weights(1:ns:2) = 2.0_dp/3.0_dp
      weights(2:ns - 1:2) = 4.0_dp/3.0_dp
      weights(1) = 1.0_dp/3.0_dp
      weights(ns) = 1.0_dp/3.0_dp
      IF (ns < n) weights(n - 1:n) = weights(n - 1:n) + 0.5_dp
      weights(:) = weights*pot%rab
   END FUNCTION upf_radial_weights

! **************************************************************************************************
!> \brief Advance past standalone XML comments before tokenization treats ! as an input comment.
!> \param parser UPF parser.
!> \param at_end Whether the file ended.
! **************************************************************************************************
   SUBROUTINE upf_next_line(parser, at_end)
      TYPE(cp_parser_type), INTENT(INOUT)                :: parser
      LOGICAL, INTENT(OUT)                               :: at_end

      CHARACTER(LEN=:), ALLOCATABLE                      :: line
      LOGICAL                                            :: in_comment

      in_comment = .FALSE.
      DO
         CALL parser_get_next_line(parser, 1, at_end)
         IF (at_end) THEN
            CPASSERT(.NOT. in_comment)
            RETURN
         END IF
         line = ADJUSTL(parser%buffer%input_lines(parser%buffer%present_line_number))
         IF (.NOT. in_comment .AND. INDEX(line, '<!--') /= 1) RETURN
         in_comment = INDEX(line, '-->') == 0
      END DO
   END SUBROUTINE upf_next_line

! **************************************************************************************************
!> \brief ...
!> \param pot ...
!> \param upf_filename ...
!> \param read_header ...
!> \param read_augmentation Import scalar USPP augmentation data explicitly.
!>        This reader option does not enable USPP in a calculation method.
!> \param read_paw Import scalar PAW partial waves and one-center data, including augmentation.
!>        Calculation methods must explicitly support the PAW energy functional before using this option.
! **************************************************************************************************
   SUBROUTINE atom_read_upf(pot, upf_filename, read_header, read_augmentation, read_paw)

      TYPE(atom_upfpot_type)                             :: pot
      CHARACTER(len=*), INTENT(IN)                       :: upf_filename
      LOGICAL, INTENT(IN), OPTIONAL                      :: read_header, read_augmentation, read_paw

      CHARACTER(LEN=2)                                   :: symbol
      INTEGER                                            :: l, ncore, nel
      LOGICAL                                            :: readall, readaug, readpaw

      IF (PRESENT(read_header)) THEN
         readall = .NOT. read_header
      ELSE
         readall = .TRUE.
      END IF
      readaug = .FALSE.
      IF (PRESENT(read_augmentation)) readaug = read_augmentation
      readpaw = .FALSE.
      IF (PRESENT(read_paw)) readpaw = read_paw
      readaug = readaug .OR. readpaw

      ! filename
      pot%filename = ADJUSTL(TRIM(upf_filename))

      ! Ignore json potentials as SIRIUS will parse those on its own.
      l = LEN_TRIM(pot%filename)
      IF (l >= 5) THEN
         IF (pot%filename(l - 4:l) == '.json') THEN
            pot%zion = 0.0_dp
            RETURN
         END IF
      END IF

      CALL atom_read_upf_v2(pot, upf_filename, readall, readaug, readpaw)

      ! set up econf
      IF (SUM(pot%econf) == 0) THEN
         symbol = ADJUSTL(TRIM(pot%symbol))
         CALL get_ptable_info(symbol, number=ncore)
         pot%econf(0:3) = ptable(ncore)%e_conv(0:3)
         nel = NINT(ncore - pot%zion)
         SELECT CASE (nel)
         CASE DEFAULT
            CPABORT("Unknown Core State")
         CASE (0)
            ! no core electron
         CASE (2)
            pot%econf(0:3) = pot%econf(0:3) - ptable(2)%e_conv(0:3)
         CASE (10)
            pot%econf(0:3) = pot%econf(0:3) - ptable(10)%e_conv(0:3)
         CASE (18)
            pot%econf(0:3) = pot%econf(0:3) - ptable(18)%e_conv(0:3)
         CASE (28)
            pot%econf(0:3) = pot%econf(0:3) - ptable(18)%e_conv(0:3)
            pot%econf(2) = pot%econf(2) - 10
         CASE (36)
            pot%econf(0:3) = pot%econf(0:3) - ptable(36)%e_conv(0:3)
         CASE (46)
            pot%econf(0:3) = pot%econf(0:3) - ptable(36)%e_conv(0:3)
            pot%econf(2) = pot%econf(2) - 10
         CASE (54)
            pot%econf(0:3) = pot%econf(0:3) - ptable(54)%e_conv(0:3)
         CASE (60)
            pot%econf(0:3) = pot%econf(0:3) - ptable(36)%e_conv(0:3)
            pot%econf(2) = pot%econf(2) - 10
            pot%econf(3) = pot%econf(3) - 14
         CASE (62)
            ! [Kr] 4d10 4f14 5s2 core, with 5p semicore in the valence (e.g. Pt).
            pot%econf(0:3) = pot%econf(0:3) - [10, 18, 20, 14]
         CASE (68)
            pot%econf(0:3) = pot%econf(0:3) - ptable(54)%e_conv(0:3)
            pot%econf(3) = pot%econf(3) - 14
         CASE (78)
            pot%econf(0:3) = pot%econf(0:3) - ptable(54)%e_conv(0:3)
            pot%econf(2) = pot%econf(2) - 10
            pot%econf(3) = pot%econf(3) - 14
         END SELECT
         !
         CPASSERT(ALL(pot%econf >= 0))
      END IF

      ! name
      IF (pot%pname == "") THEN
         pot%pname = ADJUSTL(TRIM(pot%symbol))
      END IF

   END SUBROUTINE atom_read_upf

! **************************************************************************************************
!> \brief ...
!> \param pot ...
!> \param upf_filename ...
!> \param readall ...
!> \param readaug Read scalar USPP augmentation arrays.
!> \param readpaw Read scalar PAW one-center data.
! **************************************************************************************************
   SUBROUTINE atom_read_upf_v2(pot, upf_filename, readall, readaug, readpaw)

      TYPE(atom_upfpot_type)                             :: pot
      CHARACTER(len=*), INTENT(IN)                       :: upf_filename
      LOGICAL, INTENT(IN)                                :: readall, readaug, readpaw

      CHARACTER(LEN=default_string_length)               :: nametag
      INTEGER                                            :: ib, ntag
      LOGICAL                                            :: at_end
      TYPE(cp_parser_type)                               :: parser

      ntag = 0
      CALL parser_create(parser, upf_filename, apply_preprocessing=.FALSE.)
      DO
         at_end = .FALSE.
         CALL upf_next_line(parser, at_end)
         IF (at_end) EXIT
         CALL parser_get_object(parser, nametag, lower_to_upper=.TRUE.)
         IF (nametag(1:1) /= "<") CYCLE
         IF (ntag == 0) THEN
            ! we are looking for UPF tag
            IF (nametag(2:4) == "UPF") THEN
               CALL parser_get_object(parser, nametag, lower_to_upper=.TRUE.)
               ! read UPF file version
               CALL parser_get_object(parser, nametag, lower_to_upper=.TRUE.)
               pot%version = TRIM(nametag)
               CPASSERT(nametag(1:5) == "2.0.1")
               CALL parser_get_object(parser, nametag, lower_to_upper=.TRUE.)
               CPASSERT(nametag(1:1) == ">")
               ntag = 1
            END IF
         ELSE IF (ntag == 1) THEN
            ! we are looking for 1st level tags
            IF (nametag(2:8) == "PP_INFO") THEN
               CPASSERT(nametag(9:9) == ">")
               CALL upf_info_section(parser, pot)
            ELSE IF (nametag(2:10) == "PP_HEADER") THEN
               IF (.NOT. (nametag(11:11) == ">")) THEN
                  CALL upf_header_option(parser, pot)
               END IF
               IF (readpaw .AND. pot%is_paw) THEN
                  IF (pot%has_so) CPABORT("Native PAW import requires scalar-relativistic data")
                  CPASSERT(.NOT. ALLOCATED(pot%paw))
                  ALLOCATE (pot%paw)
               END IF
            ELSE IF (nametag(2:8) == "PP_MESH") THEN
               IF (.NOT. (nametag(9:9) == ">")) THEN
                  CALL upf_mesh_option(parser, pot)
               END IF
               CALL upf_mesh_section(parser, pot)
            ELSE IF (nametag(2:8) == "PP_NLCC") THEN
               IF (nametag(9:9) == ">") THEN
                  CALL upf_nlcc_section(parser, pot, .FALSE.)
               ELSE
                  CALL upf_nlcc_section(parser, pot, .TRUE.)
               END IF
            ELSE IF (nametag(2:9) == "PP_LOCAL") THEN
               IF (nametag(10:10) == ">") THEN
                  CALL upf_local_section(parser, pot, .FALSE.)
               ELSE
                  CALL upf_local_section(parser, pot, .TRUE.)
               END IF
            ELSE IF (nametag(2:12) == "PP_NONLOCAL") THEN
               CPASSERT(nametag(13:13) == ">")
               CALL upf_nonlocal_section(parser, pot, readaug .AND. (.NOT. pot%is_paw .OR. readpaw))
            ELSE IF (nametag(2:13) == "PP_SEMILOCAL") THEN
               CALL upf_semilocal_section(parser, pot)
            ELSE IF (nametag(2:9) == "PP_PSWFC") THEN
               ! skip section for now
            ELSE IF (nametag(2:12) == "PP_FULL_WFC") THEN
               IF (ALLOCATED(pot%paw)) CALL upf_full_wfc_section(parser, pot, nametag)
            ELSE IF (nametag(2:11) == "PP_RHOATOM") THEN
               ! skip section for now
            ELSE IF (nametag(2:7) == "PP_PAW") THEN
               IF (ALLOCATED(pot%paw)) CALL upf_paw_section(parser, pot)
            ELSE IF (nametag(2:6) == "/UPF>") THEN
               EXIT
            END IF
         END IF
      END DO
      CALL parser_release(parser)

      CPASSERT(ntag > 0)

      ! rescale projectors
      IF (ALLOCATED(pot%beta)) THEN
         CPASSERT(ALLOCATED(pot%r) .AND. ALLOCATED(pot%rab))
         DO ib = 1, pot%number_of_proj
            IF (pot%r(1) == 0.0_dp) THEN
               pot%beta(2:, ib) = pot%beta(2:, ib)/pot%r(2:)
            ELSE
               pot%beta(:, ib) = pot%beta(:, ib)/pot%r(:)
            END IF
         END DO
      END IF

      ! test for not supported options
      IF (readall) THEN
         IF (pot%is_ultrasoft .AND. .NOT. readaug) THEN
            CPABORT("UPF ultrasoft pseudopotential not implemented")
         END IF
         IF (pot%is_paw .AND. .NOT. readpaw) THEN
            CPABORT("UPF PAW potential not implemented")
         END IF
      END IF
      IF (readaug .AND. (pot%is_ultrasoft .OR. pot%is_paw) .AND. (.NOT. pot%is_paw .OR. readpaw)) THEN
         CPASSERT(ALLOCATED(pot%qij) .AND. ALLOCATED(pot%qijl))
      END IF
      IF (ALLOCATED(pot%paw)) THEN
         CPASSERT(ALLOCATED(pot%paw%ae_wfc) .AND. ALLOCATED(pot%paw%ps_wfc))
         CPASSERT(ALLOCATED(pot%paw%occupations) .AND. ALLOCATED(pot%paw%rho_core))
         CPASSERT(ALLOCATED(pot%paw%vlocal) .AND. ALLOCATED(pot%paw%multipoles))
         CPASSERT(ALL(pot%paw%angular_momenta == pot%lbeta))
         CPASSERT(ALL(ieee_is_finite(pot%qij)) .AND. ALL(ieee_is_finite(pot%qijl)))
      END IF

   END SUBROUTINE atom_read_upf_v2

! **************************************************************************************************
!> \brief ...
!> \param parser ...
!> \param pot ...
! **************************************************************************************************
   SUBROUTINE upf_info_section(parser, pot)
      TYPE(cp_parser_type), INTENT(INOUT)                :: parser
      TYPE(atom_upfpot_type)                             :: pot

      CHARACTER(LEN=default_string_length)               :: line, string
      INTEGER                                            :: icount, iline
      LOGICAL                                            :: at_end

      icount = 0
      DO
         CALL upf_next_line(parser, at_end)
         CPASSERT(.NOT. at_end)
         iline = parser%buffer%present_line_number
         line = TRIM(parser%buffer%input_lines(iline))
         CALL parser_get_object(parser, string)
         IF (string(1:10) == "</PP_INFO>") EXIT
         icount = icount + 1
         IF (icount > pot%maxinfo) CYCLE
         pot%info(icount) = line
      END DO
      pot%maxinfo = icount

   END SUBROUTINE upf_info_section

! **************************************************************************************************
!> \brief ...
!> \param parser ...
!> \param pot ...
! **************************************************************************************************
   SUBROUTINE upf_header_option(parser, pot)
      TYPE(cp_parser_type), INTENT(INOUT)                :: parser
      TYPE(atom_upfpot_type)                             :: pot

      CHARACTER(LEN=default_string_length)               :: string
      LOGICAL                                            :: at_end

      DO
         IF (parser_test_next_token(parser) == "EOL") THEN
            CALL upf_next_line(parser, at_end)
            CPASSERT(.NOT. at_end)
         END IF
         CALL parser_get_object(parser, string, lower_to_upper=.TRUE.)
         IF (string == "/>") EXIT
         SELECT CASE (string)
         CASE ("GENERATED")
            CALL parser_get_object(parser, pot%generated, allow_empty=.TRUE.)
         CASE ("AUTHOR")
            CALL parser_get_object(parser, pot%author, allow_empty=.TRUE.)
         CASE ("DATE")
            CALL parser_get_object(parser, pot%date, allow_empty=.TRUE.)
         CASE ("COMMENT")
            CALL parser_get_object(parser, pot%comment, allow_empty=.TRUE.)
         CASE ("ELEMENT")
            CALL parser_get_object(parser, pot%symbol)
            CPASSERT(2 <= LEN(pot%symbol))
         CASE ("PSEUDO_TYPE")
            CALL parser_get_object(parser, pot%pseudo_type)
         CASE ("RELATIVISTIC")
            CALL parser_get_object(parser, pot%relativistic)
         CASE ("IS_ULTRASOFT")
            CALL parser_get_object(parser, pot%is_ultrasoft)
         CASE ("IS_PAW")
            CALL parser_get_object(parser, pot%is_paw)
         CASE ("IS_COULOMB")
            CALL parser_get_object(parser, pot%is_coulomb)
         CASE ("HAS_SO")
            CALL parser_get_object(parser, pot%has_so)
         CASE ("HAS_WFC")
            CALL parser_get_object(parser, pot%has_wfc)
         CASE ("HAS_GIPAW")
            CALL parser_get_object(parser, pot%has_gipaw)
         CASE ("PAW_AS_GIPAW")
            CALL parser_get_object(parser, pot%paw_as_gipaw)
         CASE ("CORE_CORRECTION")
            CALL parser_get_object(parser, pot%core_correction)
         CASE ("FUNCTIONAL")
            CALL parser_get_object(parser, pot%functional)
         CASE ("Z_VALENCE")
            CALL parser_get_object(parser, pot%zion)
         CASE ("TOTAL_PSENERGY")
            CALL parser_get_object(parser, pot%total_psenergy)
         CASE ("WFC_CUTOFF")
            CALL parser_get_object(parser, pot%wfc_cutoff)
         CASE ("RHO_CUTOFF")
            CALL parser_get_object(parser, pot%rho_cutoff)
         CASE ("L_MAX")
            CALL parser_get_object(parser, pot%l_max)
         CASE ("L_MAX_RHO")
            CALL parser_get_object(parser, pot%l_max_rho)
         CASE ("L_LOCAL")
            CALL parser_get_object(parser, pot%l_local)
         CASE ("MESH_SIZE")
            CALL parser_get_object(parser, pot%mesh_size)
         CASE ("NUMBER_OF_WFC")
            CALL parser_get_object(parser, pot%number_of_wfc)
         CASE ("NUMBER_OF_PROJ")
            CALL parser_get_object(parser, pot%number_of_proj)
         CASE DEFAULT
            CPWARN(string)
            CALL cp_abort(__LOCATION__, "Error while parsing UPF header: "// &
                          "Adjust format of delimiters ... only double quotes are admissible.")
         END SELECT
      END DO

   END SUBROUTINE upf_header_option

! **************************************************************************************************
!> \brief ...
!> \param parser ...
!> \param pot ...
! **************************************************************************************************
   SUBROUTINE upf_mesh_option(parser, pot)
      TYPE(cp_parser_type), INTENT(INOUT)                :: parser
      TYPE(atom_upfpot_type)                             :: pot

      CHARACTER(LEN=default_string_length)               :: string
      INTEGER                                            :: jj
      LOGICAL                                            :: at_end

      DO
         IF (parser_test_next_token(parser) == "EOL") THEN
            CALL upf_next_line(parser, at_end)
            CPASSERT(.NOT. at_end)
         END IF
         CALL parser_get_object(parser, string, lower_to_upper=.TRUE.)
         IF (string == ">") EXIT
         SELECT CASE (string)
         CASE ("DX")
            CALL parser_get_object(parser, pot%dx)
         CASE ("XMIN")
            CALL parser_get_object(parser, pot%xmin)
         CASE ("RMAX")
            CALL parser_get_object(parser, pot%rmax)
         CASE ("MESH")
            CALL parser_get_object(parser, jj)
            CPASSERT(pot%mesh_size == jj)
         CASE ("ZMESH")
            CALL parser_get_object(parser, pot%zmesh)
         CASE DEFAULT
            CPABORT("Unknown UPF PP_MESH option <"//TRIM(string)//"> found")
         END SELECT

      END DO

   END SUBROUTINE upf_mesh_option

! **************************************************************************************************
!> \brief ...
!> \param parser ...
!> \param pot ...
! **************************************************************************************************
   SUBROUTINE upf_mesh_section(parser, pot)
      TYPE(cp_parser_type), INTENT(INOUT)                :: parser
      TYPE(atom_upfpot_type)                             :: pot

      CHARACTER(LEN=default_string_length)               :: line, string, string2
      INTEGER                                            :: icount, m, mc, ms
      LOGICAL                                            :: at_end

      DO
         CALL upf_next_line(parser, at_end)
         CPASSERT(.NOT. at_end)
         CALL parser_get_object(parser, string, lower_to_upper=.TRUE.)
         SELECT CASE (string)
         CASE ("<PP_R", "<PP_R>")
            m = pot%mesh_size
            ms = pot%mesh_size
            mc = 1
            IF (string /= "<PP_R>") THEN
               ! options
               DO
                  IF (parser_test_next_token(parser) == "EOL") THEN
                     CALL upf_next_line(parser, at_end)
                     CPASSERT(.NOT. at_end)
                  END IF
                  CALL parser_get_object(parser, string2, lower_to_upper=.TRUE.)
                  IF (string2 == ">") EXIT
                  SELECT CASE (string2)
                  CASE ("TYPE")
                     CALL parser_get_object(parser, line, lower_to_upper=.TRUE.)
                     CPASSERT(line == "REAL")
                  CASE ("SIZE")
                     CALL parser_get_object(parser, ms)
                     CPASSERT(ms <= m)
                  CASE ("COLUMNS")
                     CALL parser_get_object(parser, mc)
                  CASE DEFAULT
                     CPABORT("Unknown UPF PP_R option <"//TRIM(string2)//"> found")
                  END SELECT
               END DO
            END IF
            ALLOCATE (pot%r(m))
            pot%r = 0.0_dp
            icount = 1
            DO WHILE (icount <= ms)
               IF (parser_test_next_token(parser) == "EOL") THEN
                  CALL upf_next_line(parser, at_end)
                  CPASSERT(.NOT. at_end)
               ELSE
                  CALL parser_get_object(parser, pot%r(icount))
                  icount = icount + 1
               END IF
            END DO
         CASE ("<PP_RAB", "<PP_RAB>")
            m = pot%mesh_size
            ms = pot%mesh_size
            mc = 1
            IF (string /= "<PP_RAB>") THEN
               ! options
               DO
                  IF (parser_test_next_token(parser) == "EOL") THEN
                     CALL upf_next_line(parser, at_end)
                     CPASSERT(.NOT. at_end)
                  END IF
                  CALL parser_get_object(parser, string2, lower_to_upper=.TRUE.)
                  IF (string2 == ">") EXIT
                  SELECT CASE (string2)
                  CASE ("TYPE")
                     CALL parser_get_object(parser, line, lower_to_upper=.TRUE.)
                     CPASSERT(line == "REAL")
                  CASE ("SIZE")
                     CALL parser_get_object(parser, ms)
                     CPASSERT(ms <= m)
                  CASE ("COLUMNS")
                     CALL parser_get_object(parser, mc)
                  CASE DEFAULT
                     CPABORT("Unknown UPF PP_RAB option <"//TRIM(string2)//"> found")
                  END SELECT
               END DO
            END IF
            ALLOCATE (pot%rab(m))
            pot%rab = 0.0_dp
            icount = 1
            DO WHILE (icount <= ms)
               IF (parser_test_next_token(parser) == "EOL") THEN
                  CALL upf_next_line(parser, at_end)
                  CPASSERT(.NOT. at_end)
               ELSE
                  CALL parser_get_object(parser, pot%rab(icount))
                  icount = icount + 1
               END IF
            END DO
         CASE ("</PP_MESH>")
            EXIT
         CASE DEFAULT
            !
         END SELECT
      END DO

   END SUBROUTINE upf_mesh_section

! **************************************************************************************************
!> \brief ...
!> \param parser ...
!> \param pot ...
!> \param options ...
! **************************************************************************************************
   SUBROUTINE upf_nlcc_section(parser, pot, options)
      TYPE(cp_parser_type), INTENT(INOUT)                :: parser
      TYPE(atom_upfpot_type)                             :: pot
      LOGICAL, INTENT(IN)                                :: options

      CHARACTER(LEN=default_string_length)               :: line, string
      INTEGER                                            :: icount, m, mc, ms
      LOGICAL                                            :: at_end

      m = pot%mesh_size
      ms = m
      mc = 1
      IF (options) THEN
         DO
            IF (parser_test_next_token(parser) == "EOL") THEN
               CALL upf_next_line(parser, at_end)
               CPASSERT(.NOT. at_end)
            END IF
            CALL parser_get_object(parser, string, lower_to_upper=.TRUE.)
            IF (string == ">") EXIT
            SELECT CASE (string)
            CASE ("TYPE")
               CALL parser_get_object(parser, line, lower_to_upper=.TRUE.)
               CPASSERT(line == "REAL")
            CASE ("SIZE")
               CALL parser_get_object(parser, ms)
               CPASSERT(ms <= m)
            CASE ("COLUMNS")
               CALL parser_get_object(parser, mc)
            CASE DEFAULT
               CPABORT("Unknown UPF PP_NLCC option <"//TRIM(string)//"> found")
            END SELECT
         END DO
      END IF

      ALLOCATE (pot%rho_nlcc(m))
      pot%rho_nlcc = 0.0_dp
      icount = 1
      DO WHILE (icount <= ms)
         IF (parser_test_next_token(parser) == "EOL") THEN
            CALL upf_next_line(parser, at_end)
            CPASSERT(.NOT. at_end)
         ELSE
            CALL parser_get_object(parser, pot%rho_nlcc(icount))
            icount = icount + 1
         END IF
      END DO

      CALL upf_next_line(parser, at_end)
      CPASSERT(.NOT. at_end)
      CALL parser_get_object(parser, string, lower_to_upper=.TRUE.)
      CPASSERT(string == "</PP_NLCC>")

   END SUBROUTINE upf_nlcc_section

! **************************************************************************************************
!> \brief ...
!> \param parser ...
!> \param pot ...
!> \param options ...
! **************************************************************************************************
   SUBROUTINE upf_local_section(parser, pot, options)
      TYPE(cp_parser_type), INTENT(INOUT)                :: parser
      TYPE(atom_upfpot_type)                             :: pot
      LOGICAL, INTENT(IN)                                :: options

      CHARACTER(LEN=default_string_length)               :: line, string
      INTEGER                                            :: icount, m, mc, ms
      LOGICAL                                            :: at_end

      m = pot%mesh_size
      ms = m
      mc = 1
      IF (options) THEN
         DO
            IF (parser_test_next_token(parser) == "EOL") THEN
               CALL upf_next_line(parser, at_end)
               CPASSERT(.NOT. at_end)
            END IF
            CALL parser_get_object(parser, string, lower_to_upper=.TRUE.)
            IF (string == ">") EXIT
            SELECT CASE (string)
            CASE ("TYPE")
               CALL parser_get_object(parser, line, lower_to_upper=.TRUE.)
               CPASSERT(line == "REAL")
            CASE ("SIZE")
               CALL parser_get_object(parser, ms)
               CPASSERT(ms <= m)
            CASE ("COLUMNS")
               CALL parser_get_object(parser, mc)
            CASE DEFAULT
               CPABORT("Unknown UPF PP_LOCAL option <"//TRIM(string)//"> found")
            END SELECT
         END DO
      END IF

      ALLOCATE (pot%vlocal(m))
      pot%vlocal = 0.0_dp
      icount = 1
      DO WHILE (icount <= ms)
         IF (parser_test_next_token(parser) == "EOL") THEN
            CALL upf_next_line(parser, at_end)
            CPASSERT(.NOT. at_end)
         ELSE
            CALL parser_get_object(parser, pot%vlocal(icount))
            icount = icount + 1
         END IF
      END DO

      ! Ry -> Hartree
      pot%vlocal = 0.5_dp*pot%vlocal

      CALL upf_next_line(parser, at_end)
      CPASSERT(.NOT. at_end)
      CALL parser_get_object(parser, string, lower_to_upper=.TRUE.)
      CPASSERT(string == "</PP_LOCAL>")

   END SUBROUTINE upf_local_section

! **************************************************************************************************
!> \brief ...
!> \param parser ...
!> \param pot ...
!> \param readaug Import angular-resolved USPP augmentation data.
! **************************************************************************************************
   SUBROUTINE upf_nonlocal_section(parser, pot, readaug)
      TYPE(cp_parser_type), INTENT(INOUT)                :: parser
      TYPE(atom_upfpot_type)                             :: pot
      LOGICAL, INTENT(IN)                                :: readaug

      CHARACTER(LEN=default_string_length)               :: line, string
      INTEGER                                            :: i1, i2, ibeta, icount, la, m, mc, ms, &
                                                            nbeta
      LOGICAL                                            :: at_end

      m = pot%mesh_size
      nbeta = pot%number_of_proj
      ALLOCATE (pot%dion(nbeta, nbeta), pot%beta(m, nbeta), pot%lbeta(nbeta))
      pot%dion = 0.0_dp
      pot%beta = 0.0_dp
      pot%lbeta = -1

      ibeta = 0
      DO
         CALL upf_next_line(parser, at_end)
         CPASSERT(.NOT. at_end)
         CALL parser_get_object(parser, string, lower_to_upper=.TRUE.)
         IF (string(1:8) == "<PP_BETA") THEN
            ms = m
            ibeta = ibeta + 1
            i1 = ibeta
            la = 0
            CPASSERT(ibeta <= nbeta)
            DO
               IF (parser_test_next_token(parser) == "EOL") THEN
                  CALL upf_next_line(parser, at_end)
                  CPASSERT(.NOT. at_end)
               END IF
               CALL parser_get_object(parser, string, lower_to_upper=.TRUE.)
               IF (string == ">") EXIT
               SELECT CASE (string)
               CASE ("TYPE")
                  CALL parser_get_object(parser, line, lower_to_upper=.TRUE.)
                  CPASSERT(line == "REAL")
               CASE ("SIZE")
                  CALL parser_get_object(parser, ms)
                  CPASSERT(ms <= m)
               CASE ("COLUMNS")
                  CALL parser_get_object(parser, mc)
               CASE ("INDEX")
                  CALL parser_get_object(parser, i1)
                  CPASSERT(i1 <= nbeta)
               CASE ("ANGULAR_MOMENTUM")
                  CALL parser_get_object(parser, la)
               CASE ("LABEL")
                  CALL parser_get_object(parser, line, allow_empty=.TRUE.)
                  ! not used currently
               CASE ("CUTOFF_RADIUS_INDEX")
                  CALL parser_get_object(parser, line)
                  ! not used currently
               CASE ("CUTOFF_RADIUS")
                  CALL parser_get_object(parser, line)
                  ! not used currently
               CASE ("ULTRASOFT_CUTOFF_RADIUS", "NORM_CONSERVING_RADIUS")
                  CALL parser_get_object(parser, line)
                  ! not used currently
               CASE DEFAULT
                  CPABORT("Unknown UPF PP_BETA option <"//TRIM(string)//"> found")
               END SELECT
            END DO
            pot%lbeta(i1) = la
            icount = 1
            DO WHILE (icount <= ms)
               IF (parser_test_next_token(parser) == "EOL") THEN
                  CALL upf_next_line(parser, at_end)
                  CPASSERT(.NOT. at_end)
               ELSE
                  CALL parser_get_object(parser, pot%beta(icount, i1))
                  icount = icount + 1
               END IF
            END DO
         ELSE IF (string(1:7) == "<PP_DIJ") THEN
            ms = nbeta*nbeta
            IF (string /= "<PP_DIJ>") THEN
               DO
                  IF (parser_test_next_token(parser) == "EOL") THEN
                     CALL upf_next_line(parser, at_end)
                     CPASSERT(.NOT. at_end)
                  END IF
                  CALL parser_get_object(parser, string, lower_to_upper=.TRUE.)
                  IF (string == ">") EXIT
                  SELECT CASE (string)
                  CASE ("TYPE")
                     CALL parser_get_object(parser, line, lower_to_upper=.TRUE.)
                     CPASSERT(line == "REAL")
                  CASE ("SIZE")
                     CALL parser_get_object(parser, ms)
                     CPASSERT(ms == nbeta*nbeta)
                  CASE ("COLUMNS")
                     CALL parser_get_object(parser, mc)
                  CASE DEFAULT
                     CPABORT("Unknown UPF PP_DIJ option <"//TRIM(string)//"> found")
                  END SELECT
               END DO
            END IF
            icount = 1
            DO WHILE (icount <= ms)
               IF (parser_test_next_token(parser) == "EOL") THEN
                  CALL upf_next_line(parser, at_end)
                  CPASSERT(.NOT. at_end)
               ELSE
                  i1 = (icount - 1)/nbeta + 1
                  i2 = MOD(icount - 1, nbeta) + 1
                  CALL parser_get_object(parser, pot%dion(i1, i2))
                  icount = icount + 1
               END IF
            END DO
         ELSE IF (string(1:16) == "<PP_AUGMENTATION") THEN
            IF (readaug) CALL upf_augmentation_section(parser, pot)
         ELSE IF (string(1:14) == "</PP_NONLOCAL>") THEN
            EXIT
         END IF
      END DO

      ! Convert the complete beta*D*beta operator from Ry to Hartree without
      ! changing the tabulated projector amplitude. This also preserves the
      ! projector convention needed by augmentation charges and PAW duality.
      ! Some NC converters put energy units into beta and inverse energy into D;
      ! preserving their numerical beta still gives the correct operator here.
      ! The radial prefactor in PP_BETA is removed after reading the whole file.
      pot%dion = 0.5_dp*pot%dion

   END SUBROUTINE upf_nonlocal_section

! **************************************************************************************************
!> \brief Read angular-momentum-resolved USPP augmentation in its tabulated convention.
!>        PP_Q is retained independently of the radial quadrature of PP_QIJL. Those
!>        quantities need not agree exactly in a finite-precision dataset. No energy
!>        conversion or projector normalization is applied to either augmentation array.
!> \param parser UPF token stream, positioned after the PP_AUGMENTATION tag name.
!> \param pot UPF data with the mesh and all beta angular momenta already read.
! **************************************************************************************************
   SUBROUTINE upf_augmentation_section(parser, pot)
      TYPE(cp_parser_type), INTENT(INOUT)                :: parser
      TYPE(atom_upfpot_type), INTENT(INOUT)              :: pot

      CHARACTER(LEN=default_string_length)               :: string, tag, value
      INTEGER                                            :: i, iraug, j, k, l, lmax, m, nb, nc, nqf, &
                                                            nqlc, pair, size_data, tag_i, tag_j, &
                                                            tag_l
      LOGICAL                                            :: at_end, is_null, q_with_l, read_q, &
                                                            self_closing
      LOGICAL, ALLOCATABLE                               :: present_channel(:, :)
      REAL(dp), ALLOCATABLE                              :: multipoles(:)

      CPASSERT(.NOT. ALLOCATED(pot%qij) .AND. .NOT. ALLOCATED(pot%qijl))
      CPASSERT(ALL(pot%lbeta >= 0))
      nb = pot%number_of_proj
      m = pot%mesh_size
      lmax = 2*MAXVAL(pot%lbeta)
      q_with_l = .FALSE.
      nqf = -1
      nqlc = lmax + 1
      iraug = -1
      DO
         IF (parser_test_next_token(parser) == "EOL") THEN
            CALL upf_next_line(parser, at_end)
            CPASSERT(.NOT. at_end)
         END IF
         CALL parser_get_object(parser, string, lower_to_upper=.TRUE.)
         IF (string == ">") EXIT
         SELECT CASE (string)
         CASE ("Q_WITH_L")
            CALL parser_get_object(parser, q_with_l)
         CASE ("NQF")
            CALL parser_get_object(parser, nqf)
         CASE ("NQLC")
            CALL parser_get_object(parser, nqlc)
         CASE ("SHAPE")
            CPASSERT(ALLOCATED(pot%paw))
            CALL parser_get_object(parser, pot%paw%augmentation_shape)
         CASE ("CUTOFF_R")
            CPASSERT(ALLOCATED(pot%paw))
            CALL parser_get_object(parser, pot%paw%cutoff_radius)
         CASE ("CUTOFF_R_INDEX")
            CPASSERT(ALLOCATED(pot%paw))
            CALL parser_get_object(parser, pot%paw%cutoff_index)
         CASE ("IRAUG")
            CPASSERT(ALLOCATED(pot%paw))
            CALL parser_get_object(parser, iraug)
         CASE ("L_MAX_AUG", "LMAX_AUG")
            CPASSERT(ALLOCATED(pot%paw))
            CALL parser_get_object(parser, pot%paw%lmax_aug)
         CASE ("AUGMENTATION_EPSILON")
            CPASSERT(ALLOCATED(pot%paw))
            CALL parser_get_object(parser, pot%paw%augmentation_epsilon)
         CASE DEFAULT
            CPABORT("Unknown UPF PP_AUGMENTATION option <"//TRIM(string)//"> found")
         END SELECT
      END DO
      IF (.NOT. q_with_l .OR. nqf /= 0) THEN
         CPABORT("Native augmentation import requires angular-resolved PP_QIJL with nqf=0")
      END IF
      IF (nqlc == 0) nqlc = lmax + 1
      CPASSERT(nqlc >= lmax + 1)
      IF (ALLOCATED(pot%paw)) THEN
         IF (pot%paw%cutoff_index < 0) pot%paw%cutoff_index = iraug
         IF (iraug >= 0) THEN
            CPASSERT(iraug == pot%paw%cutoff_index)
         END IF
         CPASSERT(pot%paw%cutoff_index >= 1 .AND. pot%paw%cutoff_index <= m)
         CPASSERT(pot%paw%lmax_aug >= 0 .AND. pot%paw%lmax_aug <= lmax)
         CPASSERT(ieee_is_finite(pot%paw%augmentation_epsilon))
         CPASSERT(pot%paw%augmentation_epsilon >= 0.0_dp)
      END IF
      nc = nb*(nb + 1)/2
      ALLOCATE (pot%qij(nb, nb), pot%qijl(m, nc, 0:lmax), present_channel(nc, 0:lmax))
      pot%qij = 0.0_dp
      pot%qijl = 0.0_dp
      present_channel = .FALSE.
      read_q = .FALSE.
      DO
         CALL upf_next_line(parser, at_end)
         CPASSERT(.NOT. at_end)
         CALL parser_get_object(parser, tag, lower_to_upper=.TRUE.)
         IF (tag == "</PP_AUGMENTATION>") EXIT
         IF (tag(1:2) == "</") CYCLE
         IF (tag == "<PP_MULTIPOLES" .OR. tag == "<PP_MULTIPOLES>") THEN
            CPASSERT(ALLOCATED(pot%paw))
            CPASSERT(.NOT. ALLOCATED(pot%paw%multipoles))
            ALLOCATE (multipoles(nb*nb*(lmax + 1)), pot%paw%multipoles(nb, nb, 0:lmax))
            CALL upf_real_array(parser, tag, multipoles)
            pot%paw%multipoles(:, :, :) = RESHAPE(multipoles, [nb, nb, lmax + 1])
            DEALLOCATE (multipoles)
            CYCLE
         END IF
         IF (.NOT. (tag(1:8) == "<PP_QIJL" .OR. tag == "<PP_Q" .OR. tag == "<PP_Q>")) THEN
            CPABORT("Unexpected UPF augmentation tag "//TRIM(tag)//" in "//TRIM(pot%filename))
         END IF
         i = -1
         j = -1
         l = -1
         pair = -1
         size_data = m
         IF (tag == "<PP_Q" .OR. tag == "<PP_Q>") size_data = nb*nb
         is_null = .FALSE.
         self_closing = .FALSE.
         IF (tag /= "<PP_Q>") THEN
            DO
               IF (parser_test_next_token(parser) == "EOL") THEN
                  CALL upf_next_line(parser, at_end)
                  CPASSERT(.NOT. at_end)
               END IF
               CALL parser_get_object(parser, string, lower_to_upper=.TRUE.)
               IF (string == ">") EXIT
               IF (string == "/>") THEN
                  self_closing = .TRUE.
                  EXIT
               END IF
               SELECT CASE (string)
               CASE ("TYPE")
                  CALL parser_get_object(parser, value, lower_to_upper=.TRUE.)
                  CPASSERT(value == "REAL")
               CASE ("SIZE")
                  CALL parser_get_object(parser, size_data)
               CASE ("COLUMNS")
                  CALL parser_get_object(parser, value)
               CASE ("FIRST_INDEX")
                  CALL parser_get_object(parser, i)
               CASE ("SECOND_INDEX")
                  CALL parser_get_object(parser, j)
               CASE ("COMPOSITE_INDEX")
                  CALL parser_get_object(parser, pair)
               CASE ("ANGULAR_MOMENTUM")
                  CALL parser_get_object(parser, l)
               CASE ("IS_NULL")
                  CALL parser_get_object(parser, is_null)
               CASE DEFAULT
                  CPABORT("Unknown UPF augmentation array option <"//TRIM(string)//"> found")
               END SELECT
            END DO
         END IF
         IF (tag(1:8) == "<PP_QIJL") THEN
            CPASSERT(i >= 1 .AND. i <= j .AND. j <= nb)
            CPASSERT(pair == j*(j - 1)/2 + i)
            CPASSERT(l >= ABS(pot%lbeta(i) - pot%lbeta(j)))
            CPASSERT(l <= pot%lbeta(i) + pot%lbeta(j))
            CPASSERT(MODULO(pot%lbeta(i) + pot%lbeta(j) - l, 2) == 0)
            ! Cross-check the tag suffix as well as the explicit attributes.
            string = tag(10:)
            DO k = 1, LEN_TRIM(string)
               IF (string(k:k) == ".") string(k:k) = " "
            END DO
            READ (string, *) tag_i, tag_j, tag_l
            CPASSERT(tag_i == i .AND. tag_j == j .AND. tag_l == l)
            CPASSERT(.NOT. present_channel(pair, l))
            present_channel(pair, l) = .TRUE.
            IF (is_null) THEN
               CPASSERT(self_closing .AND. size_data == 0)
               CYCLE
            END IF
            CPASSERT(.NOT. self_closing .AND. size_data == m)
            DO k = 1, m
               DO WHILE (parser_test_next_token(parser) == "EOL")
                  CALL upf_next_line(parser, at_end)
                  CPASSERT(.NOT. at_end)
               END DO
               CALL parser_get_object(parser, pot%qijl(k, pair, l))
            END DO
         ELSE
            CPASSERT(.NOT. read_q .AND. .NOT. self_closing .AND. .NOT. is_null)
            CPASSERT(size_data == nb*nb)
            DO j = 1, nb
               DO i = 1, nb
                  DO WHILE (parser_test_next_token(parser) == "EOL")
                     CALL upf_next_line(parser, at_end)
                     CPASSERT(.NOT. at_end)
                  END DO
                  CALL parser_get_object(parser, pot%qij(i, j))
               END DO
            END DO
            read_q = .TRUE.
         END IF
      END DO
      CPASSERT(read_q)
      IF (ALLOCATED(pot%paw)) THEN
         CPASSERT(ALLOCATED(pot%paw%multipoles))
      END IF
      DO j = 1, nb
         DO i = 1, j
            pair = j*(j - 1)/2 + i
            DO l = ABS(pot%lbeta(i) - pot%lbeta(j)), pot%lbeta(i) + pot%lbeta(j), 2
               ! PAW writers may omit channels below the stored multipole threshold.
               IF (.NOT. present_channel(pair, l) .AND. ALLOCATED(pot%paw)) THEN
                  IF (ABS(pot%paw%multipoles(i, j, l)) < pot%paw%augmentation_epsilon) CYCLE
               END IF
               CPASSERT(present_channel(pair, l))
            END DO
         END DO
      END DO
   END SUBROUTINE upf_augmentation_section

! **************************************************************************************************
!> \brief Read a complete real UPF array and check its dimensions and closing tag.
!> \param parser Token stream after the opening tag name.
!> \param tag Opening tag, with or without its final greater-than sign.
!> \param values Array in the original file convention.
!> \param index_value Optional projector index attribute, or -1 if absent.
!> \param angular_momentum Optional angular momentum attribute, or -1 if absent.
! **************************************************************************************************
   SUBROUTINE upf_real_array(parser, tag, values, index_value, angular_momentum)
      TYPE(cp_parser_type), INTENT(INOUT)                :: parser
      CHARACTER(LEN=*), INTENT(IN)                       :: tag
      REAL(dp), INTENT(OUT)                              :: values(:)
      INTEGER, INTENT(OUT), OPTIONAL                     :: index_value, angular_momentum

      CHARACTER(LEN=default_string_length)               :: name, string, value
      INTEGER                                            :: i, last, n
      LOGICAL                                            :: at_end

      n = SIZE(values)
      last = LEN_TRIM(tag)
      IF (PRESENT(index_value)) index_value = -1
      IF (PRESENT(angular_momentum)) angular_momentum = -1
      IF (tag(last:last) == ">") THEN
         name = tag(2:last - 1)
      ELSE
         name = tag(2:last)
         DO
            IF (parser_test_next_token(parser) == "EOL") THEN
               CALL upf_next_line(parser, at_end)
               CPASSERT(.NOT. at_end)
            END IF
            CALL parser_get_object(parser, string, lower_to_upper=.TRUE.)
            IF (string == ">") EXIT
            SELECT CASE (string)
            CASE ("TYPE")
               CALL parser_get_object(parser, value, lower_to_upper=.TRUE.)
               CPASSERT(value == "REAL")
            CASE ("SIZE")
               CALL parser_get_object(parser, n)
               CPASSERT(n == SIZE(values))
            CASE ("COLUMNS")
               CALL parser_get_object(parser, value)
            CASE ("INDEX")
               CPASSERT(PRESENT(index_value))
               CALL parser_get_object(parser, index_value)
            CASE ("L")
               CPASSERT(PRESENT(angular_momentum))
               CALL parser_get_object(parser, angular_momentum)
            CASE ("LABEL")
               CALL parser_get_object(parser, value, allow_empty=.TRUE.)
            CASE DEFAULT
               CPABORT("Unknown UPF "//TRIM(name)//" option <"//TRIM(string)//"> found")
            END SELECT
         END DO
      END IF
      DO i = 1, n
         DO WHILE (parser_test_next_token(parser) == "EOL")
            CALL upf_next_line(parser, at_end)
            CPASSERT(.NOT. at_end)
         END DO
         CALL parser_get_object(parser, values(i))
      END DO
      CPASSERT(ALL(ieee_is_finite(values)))
      CALL upf_next_line(parser, at_end)
      CPASSERT(.NOT. at_end)
      CALL parser_get_object(parser, string, lower_to_upper=.TRUE.)
      CPASSERT(string == "</"//TRIM(name)//">")
   END SUBROUTINE upf_real_array

! **************************************************************************************************
!> \brief Import the two projector-indexed PAW partial-wave sets, retaining r*phi(r).
!> \param parser Token stream after the PP_FULL_WFC tag name.
!> \param pot UPF data.
!> \param opening_tag PP_FULL_WFC opening tag.
! **************************************************************************************************
   SUBROUTINE upf_full_wfc_section(parser, pot, opening_tag)
      TYPE(cp_parser_type), INTENT(INOUT)                :: parser
      TYPE(atom_upfpot_type), INTENT(INOUT)              :: pot
      CHARACTER(LEN=*), INTENT(IN)                       :: opening_tag

      CHARACTER(LEN=default_string_length)               :: string, suffix, tag
      INTEGER                                            :: i, index_value, ios, l, last, m, n, nb, &
                                                            wave_type
      LOGICAL                                            :: at_end
      LOGICAL, ALLOCATABLE                               :: found(:, :)
      REAL(dp), ALLOCATABLE                              :: values(:)

      nb = pot%number_of_proj
      m = pot%mesh_size
      CPASSERT(nb > 0 .AND. m > 0)
      CPASSERT(.NOT. ALLOCATED(pot%paw%ae_wfc))
      IF (opening_tag /= "<PP_FULL_WFC>") THEN
         DO
            IF (parser_test_next_token(parser) == "EOL") THEN
               CALL upf_next_line(parser, at_end)
               CPASSERT(.NOT. at_end)
            END IF
            CALL parser_get_object(parser, string, lower_to_upper=.TRUE.)
            IF (string == ">") EXIT
            CPASSERT(string == "NUMBER_OF_WFC")
            CALL parser_get_object(parser, n)
            CPASSERT(n == nb)
         END DO
      END IF
      ALLOCATE (pot%paw%ae_wfc(m, nb), pot%paw%ps_wfc(m, nb), pot%paw%angular_momenta(nb))
      ALLOCATE (found(nb, 2), values(m))
      found = .FALSE.
      pot%paw%angular_momenta = -1
      DO
         CALL upf_next_line(parser, at_end)
         CPASSERT(.NOT. at_end)
         CALL parser_get_object(parser, tag, lower_to_upper=.TRUE.)
         IF (tag == "</PP_FULL_WFC>") EXIT
         SELECT CASE (tag(1:10))
         CASE ("<PP_AEWFC.")
            wave_type = 1
         CASE ("<PP_PSWFC.")
            wave_type = 2
         CASE DEFAULT
            CPABORT("Unknown UPF PAW partial-wave tag "//TRIM(tag))
         END SELECT
         last = LEN_TRIM(tag)
         IF (tag(last:last) == ">") last = last - 1
         suffix = tag(11:last)
         READ (suffix, *, IOSTAT=ios) i
         CPASSERT(ios == 0)
         CPASSERT(i >= 1 .AND. i <= nb)
         CPASSERT(.NOT. found(i, wave_type))
         CALL upf_real_array(parser, tag, values, index_value, l)
         IF (index_value /= -1) THEN
            CPASSERT(index_value == i)
         END IF
         CPASSERT(l >= 0)
         IF (pot%paw%angular_momenta(i) >= 0) THEN
            CPASSERT(l == pot%paw%angular_momenta(i))
         END IF
         pot%paw%angular_momenta(i) = l
         IF (wave_type == 1) THEN
            pot%paw%ae_wfc(:, i) = values
         ELSE
            pot%paw%ps_wfc(:, i) = values
         END IF
         found(i, wave_type) = .TRUE.
      END DO
      CPASSERT(ALL(found))
   END SUBROUTINE upf_full_wfc_section

! **************************************************************************************************
!> \brief Read PAW format-2 occupations, AE core density and local potential.
!> \param parser Token stream after the PP_PAW tag name.
!> \param pot UPF data; energies and potentials are converted from Ry to Hartree.
! **************************************************************************************************
   SUBROUTINE upf_paw_section(parser, pot)
      TYPE(cp_parser_type), INTENT(INOUT)                :: parser
      TYPE(atom_upfpot_type), INTENT(INOUT)              :: pot

      CHARACTER(LEN=default_string_length)               :: string, tag
      LOGICAL                                            :: at_end

      DO
         IF (parser_test_next_token(parser) == "EOL") THEN
            CALL upf_next_line(parser, at_end)
            CPASSERT(.NOT. at_end)
         END IF
         CALL parser_get_object(parser, string, lower_to_upper=.TRUE.)
         IF (string == ">") EXIT
         SELECT CASE (string)
         CASE ("PAW_DATA_FORMAT")
            CPASSERT(pot%paw%data_format == -1)
            CALL parser_get_object(parser, pot%paw%data_format)
         CASE ("CORE_ENERGY")
            CPASSERT(.NOT. pot%paw%has_core_energy)
            CALL parser_get_object(parser, pot%paw%core_energy)
            CPASSERT(ieee_is_finite(pot%paw%core_energy))
            pot%paw%core_energy = 0.5_dp*pot%paw%core_energy
            pot%paw%has_core_energy = .TRUE.
         CASE DEFAULT
            CPABORT("Unknown UPF PP_PAW option <"//TRIM(string)//"> found")
         END SELECT
      END DO
      IF (pot%paw%data_format /= 2) CPABORT("Native PAW import requires paw_data_format=2")
      DO
         CALL upf_next_line(parser, at_end)
         CPASSERT(.NOT. at_end)
         CALL parser_get_object(parser, tag, lower_to_upper=.TRUE.)
         SELECT CASE (tag)
         CASE ("</PP_PAW>")
            EXIT
         CASE ("<PP_OCCUPATIONS", "<PP_OCCUPATIONS>")
            CPASSERT(.NOT. ALLOCATED(pot%paw%occupations))
            ALLOCATE (pot%paw%occupations(pot%number_of_proj))
            CALL upf_real_array(parser, tag, pot%paw%occupations)
         CASE ("<PP_AE_NLCC", "<PP_AE_NLCC>")
            CPASSERT(.NOT. ALLOCATED(pot%paw%rho_core))
            ALLOCATE (pot%paw%rho_core(pot%mesh_size))
            CALL upf_real_array(parser, tag, pot%paw%rho_core)
         CASE ("<PP_AE_VLOC", "<PP_AE_VLOC>")
            CPASSERT(.NOT. ALLOCATED(pot%paw%vlocal))
            ALLOCATE (pot%paw%vlocal(pot%mesh_size))
            CALL upf_real_array(parser, tag, pot%paw%vlocal)
            pot%paw%vlocal = 0.5_dp*pot%paw%vlocal
         CASE DEFAULT
            CPABORT("Unknown UPF PAW data tag "//TRIM(tag))
         END SELECT
      END DO
   END SUBROUTINE upf_paw_section

! **************************************************************************************************
!> \brief ...
!> \param parser ...
!> \param pot ...
! **************************************************************************************************
   SUBROUTINE upf_semilocal_section(parser, pot)
      TYPE(cp_parser_type), INTENT(INOUT)                :: parser
      TYPE(atom_upfpot_type)                             :: pot

      CHARACTER(LEN=default_string_length)               :: line, string
      INTEGER                                            :: i1, ib, icount, la, lmax, m, mc, ms
      LOGICAL                                            :: at_end

      m = pot%mesh_size
      lmax = pot%l_max
      ALLOCATE (pot%vsemi(m, lmax + 1))
      pot%vsemi = 0.0_dp

      ib = 0
      DO
         CALL upf_next_line(parser, at_end)
         CPASSERT(.NOT. at_end)
         CALL parser_get_object(parser, string, lower_to_upper=.TRUE.)
         IF (string(1:7) == "<PP_VNL") THEN
            ms = m
            ib = ib + 1
            i1 = ib
            la = 0
            CPASSERT(ib <= lmax + 1)
            DO
               IF (parser_test_next_token(parser) == "EOL") THEN
                  CALL upf_next_line(parser, at_end)
                  CPASSERT(.NOT. at_end)
               END IF
               CALL parser_get_object(parser, string, lower_to_upper=.TRUE.)
               IF (string == ">") EXIT
               SELECT CASE (string)
               CASE ("TYPE")
                  CALL parser_get_object(parser, line, lower_to_upper=.TRUE.)
                  CPASSERT(line == "REAL")
               CASE ("SIZE")
                  CALL parser_get_object(parser, ms)
                  CPASSERT(ms <= m)
               CASE ("COLUMNS")
                  CALL parser_get_object(parser, mc)
               CASE ("L")
                  CALL parser_get_object(parser, la)
               CASE DEFAULT
                  CPABORT("Unknown UPF PP_VNL option <"//TRIM(string)//"> found")
               END SELECT
            END DO
            i1 = la + 1
            icount = 1
            DO WHILE (icount <= ms)
               IF (parser_test_next_token(parser) == "EOL") THEN
                  CALL upf_next_line(parser, at_end)
                  CPASSERT(.NOT. at_end)
               ELSE
                  CALL parser_get_object(parser, pot%vsemi(icount, i1))
                  icount = icount + 1
               END IF
            END DO
         ELSE IF (string(1:15) == "</PP_SEMILOCAL>") THEN
            EXIT
         ELSE
            !
         END IF
      END DO
      ! Ry -> Hartree
      pot%vsemi = 0.5_dp*pot%vsemi

   END SUBROUTINE upf_semilocal_section

! **************************************************************************************************
!> \brief ...
!> \param upfpot ...
! **************************************************************************************************
   PURE SUBROUTINE atom_release_upf(upfpot)

      TYPE(atom_upfpot_type), INTENT(INOUT)              :: upfpot

      ! Intrinsic assignment releases all allocatable components, including NLCC,
      ! and clears metadata before this temporary is reused for another KIND.
      upfpot = atom_upfpot_type()

   END SUBROUTINE atom_release_upf
! **************************************************************************************************

END MODULE atom_upf
