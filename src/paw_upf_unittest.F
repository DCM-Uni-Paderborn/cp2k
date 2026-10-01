!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM paw_upf_unittest
   USE atom_upf,                        ONLY: atom_read_upf,&
                                              atom_release_upf,&
                                              atom_upfpot_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE libcp2k,                         ONLY: cp2k_finalize,&
                                              cp2k_init

   IMPLICIT NONE
   CHARACTER(LEN=*), PARAMETER :: fixture = 'paw-import-test.upf'
   CHARACTER(LEN=default_path_length) :: filename
   CHARACTER(LEN=16) :: lmax_attribute
   INTEGER :: unit, variant, i
   LOGICAL :: source_root_argument
   TYPE(atom_upfpot_type) :: pot
   TYPE(cp_logger_type), POINTER :: logger

   CALL cp2k_init()
   logger => cp_get_default_logger()
   DO variant = 1, 2
      lmax_attribute = 'l_max_aug'
      IF (variant == 2) lmax_attribute = 'lmax_aug'
      IF (logger%para_env%is_source()) THEN
         OPEN (NEWUNIT=unit, FILE=fixture, STATUS='REPLACE', ACTION='WRITE')
         WRITE (unit, '(A)') &
            '<UPF version="2.0.1">', &
            '<PP_HEADER element="H" pseudo_type="PAW" is_ultrasoft="true" is_paw="true"', &
            ' has_wfc="true" z_valence="1" mesh_size="3" number_of_proj="2" number_of_wfc="1" l_max="1"/>', &
            '<PP_MESH mesh="3">', &
            '<PP_R>', '0.5 1 1.5', '</PP_R>', &
            '<PP_RAB>', '0.5 0.5 0.5', '</PP_RAB>', '</PP_MESH>', &
            '<PP_LOCAL>', '6 4 2', '</PP_LOCAL>', &
            '<PP_NONLOCAL>', &
            '<PP_BETA.1 index="1" angular_momentum="0" size="3" norm_conserving_radius="1.2">', '1 4 9', '</PP_BETA.1>', &
            '<PP_BETA.2 index="2" angular_momentum="1" size="3">', '2 -1 3', '</PP_BETA.2>', &
            '<PP_DIJ>', '2 0 0 8', '</PP_DIJ>', &
            '<PP_AUGMENTATION q_with_l="T" nqf="0" nqlc="3" shape="PSQ" cutoff_r="-1"', &
            ' cutoff_r_index="3" iraug="3" augmentation_epsilon="1.d-12" '//TRIM(lmax_attribute)//'="2">', &
            '<PP_Q>', '2 0 0 3', '</PP_Q>', &
            '<!-- augmentation multipoles -->', &
            '<PP_MULTIPOLES type="real" size="12" columns="4">', &
            '2 0 0 3', '0 4 4 0', '0 0 0 0', '</PP_MULTIPOLES>', &
            '<PP_QIJL.1.1.0 first_index="1" second_index="1" composite_index="1" angular_momentum="0">', &
            '1 2 3', '</PP_QIJL.1.1.0>', &
            '<PP_QIJL.1.2.1 first_index="1" second_index="2" composite_index="2" angular_momentum="1">', &
            '4 -5 6', '</PP_QIJL.1.2.1>', &
            '<PP_QIJL.2.2.0 first_index="2" second_index="2" composite_index="3" angular_momentum="0">', &
            '-3 4 5', '</PP_QIJL.2.2.0>', '</PP_AUGMENTATION>', '</PP_NONLOCAL>', &
            '<PP_PSWFC>', '<PP_CHI.1 l="0">', '999 999 999', '</PP_CHI.1>', '</PP_PSWFC>', &
            '<!-- partial waves', ' one per projector -->', &
            '<PP_FULL_WFC number_of_wfc="2">', &
            '<PP_PSWFC.2 type="real" size="3" l="1" label="2P">', '2 3 4', '</PP_PSWFC.2>', &
            '<PP_AEWFC.1 index="1" l="0">', '3 5 7', '</PP_AEWFC.1>', &
            '<PP_AEWFC.2 index="2" l="1">', '5 7 9', '</PP_AEWFC.2>', &
            '<PP_PSWFC.1 l="0">', '1 2 3', '</PP_PSWFC.1>', '</PP_FULL_WFC>'
         IF (variant == 1) THEN
            WRITE (unit, '(A)') '<PP_PAW paw_data_format="2" core_energy="-8">'
         ELSE
            WRITE (unit, '(A)') '<PP_PAW paw_data_format="2">'
         END IF
         WRITE (unit, '(A)') &
            '<PP_OCCUPATIONS type="real" size="2" columns="2">', '1 0', '</PP_OCCUPATIONS>', &
            '<PP_AE_NLCC>', '7 8 9', '</PP_AE_NLCC>', &
            '<PP_AE_VLOC>', '-6 -4 -2', '</PP_AE_VLOC>', '</PP_PAW>', '</UPF>'
         CLOSE (unit)
      END IF
      CALL logger%para_env%sync()
      CALL atom_read_upf(pot, fixture, read_paw=.TRUE.)
      IF (.NOT. ALLOCATED(pot%paw)) ERROR STOP 'Missing PAW data'
      IF (ANY(SHAPE(pot%paw%ae_wfc) /= [3, 2])) ERROR STOP 'Partial-wave count differs from projector count'
      IF (ANY(pot%paw%ae_wfc /= RESHAPE([3._dp, 5._dp, 7._dp, 5._dp, 7._dp, 9._dp], [3, 2]))) THEN
         ERROR STOP 'Wrong AE partial-wave order or radial convention'
      END IF
      IF (ANY(pot%paw%ps_wfc /= RESHAPE([1._dp, 2._dp, 3._dp, 2._dp, 3._dp, 4._dp], [3, 2]))) THEN
         ERROR STOP 'Wrong pseudo partial waves or reference-wave confusion'
      END IF
      IF (ANY(pot%paw%angular_momenta /= [0, 1])) ERROR STOP 'Wrong partial-wave angular momenta'
      IF (pot%paw%has_core_energy .NEQV. (variant == 1)) ERROR STOP 'Missing core-energy flag'
      IF (variant == 1 .AND. ABS(pot%paw%core_energy + 4._dp) > 1.e-14_dp) ERROR STOP 'Core energy is not Hartree'
      IF (variant == 2 .AND. pot%paw%core_energy /= 0._dp) ERROR STOP 'Absent core energy is not zero'
      IF (ANY(pot%paw%vlocal /= [-3._dp, -2._dp, -1._dp])) ERROR STOP 'AE potential is not Hartree'
      IF (ANY(pot%paw%rho_core /= [7._dp, 8._dp, 9._dp])) ERROR STOP 'Core density normalization changed'
      IF (ANY(pot%paw%occupations /= [1._dp, 0._dp])) ERROR STOP 'Wrong occupations'
      IF (pot%paw%cutoff_index /= 3 .OR. pot%paw%lmax_aug /= 2) ERROR STOP 'Wrong augmentation metadata'
      IF (LBOUND(pot%paw%multipoles, 3) /= 0) ERROR STOP 'Multipole L index starts at zero'
      IF (pot%paw%multipoles(1, 2, 1) /= 4._dp) ERROR STOP 'Wrong multipole order'
      IF (ANY(pot%qijl(:, 3, 2) /= 0._dp)) ERROR STOP 'Missing zero-multipole channel is not zero'
      CALL atom_release_upf(pot)
      IF (ALLOCATED(pot%paw)) ERROR STOP 'PAW data not released'
      CALL logger%para_env%sync()
   END DO
   IF (logger%para_env%is_source()) THEN
      OPEN (NEWUNIT=unit, FILE=fixture, STATUS='OLD')
      CLOSE (unit, STATUS='DELETE')
   END IF

   ! Optional original datasets for an independent XML comparison.
   DO i = 1, COMMAND_ARGUMENT_COUNT()
      CALL GET_COMMAND_ARGUMENT(i, filename)
      IF (i == 1) THEN
         INQUIRE (FILE=TRIM(filename)//'/tests/UNIT_TESTS', EXIST=source_root_argument)
         IF (source_root_argument) CYCLE
      END IF
      CALL atom_read_upf(pot, TRIM(filename), read_paw=.TRUE.)
      IF (.NOT. ALLOCATED(pot%paw)) ERROR STOP 'Expected a PAW dataset'
      IF (logger%para_env%is_source()) THEN
         WRITE (*, '(A,1X,A,1X,I0,8(1X,ES25.16))') 'PAW import:', TRIM(filename), pot%number_of_proj, &
            SUM(pot%paw%ae_wfc), SUM(pot%paw%ps_wfc), SUM(pot%paw%occupations), &
            SUM(pot%paw%rho_core), SUM(pot%paw%vlocal), SUM(pot%paw%multipoles), &
            pot%paw%core_energy, SUM(pot%qijl)
      END IF
      CALL atom_release_upf(pot)
   END DO
   CALL cp2k_finalize()
END PROGRAM paw_upf_unittest
