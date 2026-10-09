!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: LGPL-3.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
! This file is part of tblite.
! SPDX-Identifier: LGPL-3.0-or-later
!
! tblite is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! tblite is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with tblite.  If not, see <https://www.gnu.org/licenses/>.

! CP2K exchange tests relocated from save_tblite; numerical assertions retained.
#if defined(__TBLITE_GXTB_ACCELERATION)
MODULE tblite_exchange_integration_tests
   USE ieee_arithmetic,                 ONLY: ieee_quiet_nan,&
                                              ieee_value
   USE iso_fortran_env,                 ONLY: int64
   USE mctc_data_vdwrad,                ONLY: get_vdw_rad
   USE mctc_env,                        ONLY: wp
   USE mctc_env_testing,                ONLY: check,&
                                              error_type,&
                                              new_unittest,&
                                              test_failed,&
                                              unittest_type
   USE mctc_io,                         ONLY: new,&
                                              structure_type
   USE mctc_io_constants,               ONLY: pi
   USE mctc_io_convert,                 ONLY: autoaa
   USE multicharge,                     ONLY: get_eeqbc_charges
   USE tblite_basis_cache,              ONLY: basis_cache
   USE tblite_basis_qvszp,              ONLY: new_qvszp_basis,&
                                              new_qvszp_cgto,&
                                              qvszp_basis_type,&
                                              qvszp_cgto_type
   USE tblite_basis_type,               ONLY: basis_type,&
                                              cgto_container
   USE tblite_bvk_transform,            ONLY: bvk_mixed_radix_k_to_r,&
                                              bvk_mixed_radix_r_to_k,&
                                              bvk_separable_k_to_r,&
                                              bvk_separable_r_to_k
   USE tblite_container_cache,          ONLY: container_cache
   USE tblite_cutoff,                   ONLY: get_lattice_points
   USE tblite_data_onecxints,           ONLY: get_onecxints
   USE tblite_exchange,                 ONLY: exchange_cache,&
                                              exchange_fock,&
                                              exchange_type,&
                                              new_exchange_fock
   USE tblite_exchange_contractions,    ONLY: get_KFock_kmesh,&
                                              get_KFock_kpoint,&
                                              get_KGrad_kmesh,&
                                              get_KGrad_kpoint,&
                                              get_bocorr_derivs_direct,&
                                              get_mulliken_derivs_direct
   USE tblite_exchange_fock,            ONLY: exchange_bvk_kernel
   USE tblite_exchange_stream,          ONLY: &
        cp2k_exchange_kmesh, cp2k_exchange_kmesh_gradient, cp2k_exchange_kpoint_gradient, &
        cp2k_exchange_partial_begin, cp2k_exchange_partial_contributor_advance, &
        cp2k_exchange_partial_descriptor_pack, cp2k_exchange_partial_descriptor_size, &
        cp2k_exchange_partial_descriptor_type, cp2k_exchange_partial_descriptor_unpack, &
        cp2k_exchange_partial_descriptors_match, cp2k_exchange_partial_discard, &
        cp2k_exchange_partial_export_forward, cp2k_exchange_partial_export_reverse, &
        cp2k_exchange_partial_forward_pack, cp2k_exchange_partial_forward_pack_range, &
        cp2k_exchange_partial_forward_size, cp2k_exchange_partial_forward_type, &
        cp2k_exchange_partial_forward_unpack, cp2k_exchange_partial_forward_unpack_range, &
        cp2k_exchange_partial_import_forward, cp2k_exchange_partial_import_forward_range, &
        cp2k_exchange_partial_import_reverse, cp2k_exchange_partial_import_reverse_range, &
        cp2k_exchange_partial_push, cp2k_exchange_partial_reverse_pack, &
        cp2k_exchange_partial_reverse_pack_range, cp2k_exchange_partial_reverse_push, &
        cp2k_exchange_partial_reverse_size, cp2k_exchange_partial_reverse_type, &
        cp2k_exchange_partial_reverse_unpack, cp2k_exchange_partial_reverse_unpack_range, &
        cp2k_exchange_partial_set_image_range, cp2k_exchange_partial_set_pull_ownership, &
        cp2k_exchange_stream_apply, cp2k_exchange_stream_batch_advance, &
        cp2k_exchange_stream_batch_apply, cp2k_exchange_stream_batch_bounds, &
        cp2k_exchange_stream_batch_pull, cp2k_exchange_stream_batch_result, &
        cp2k_exchange_stream_batched, cp2k_exchange_stream_begin, cp2k_exchange_stream_end, &
        cp2k_exchange_stream_get, cp2k_exchange_stream_has_full_mesh_storage, &
        cp2k_exchange_stream_image_batch_size, cp2k_exchange_stream_oracle, &
        cp2k_exchange_stream_pair_tile_size, cp2k_exchange_stream_peak_complex_elements, &
        cp2k_exchange_stream_push, cp2k_exchange_stream_reduced, &
        cp2k_exchange_stream_reverse_apply, cp2k_exchange_stream_reverse_batch_advance, &
        cp2k_exchange_stream_reverse_batch_apply, cp2k_exchange_stream_reverse_batch_pull, &
        cp2k_exchange_stream_reverse_begin, cp2k_exchange_stream_reverse_get, &
        cp2k_exchange_stream_reverse_push, cp2k_exchange_stream_reverse_result, &
        cp2k_exchange_stream_reverse_transform_advance, cp2k_exchange_stream_transform_advance, &
        cp2k_exchange_stream_type, cp2k_exchange_transform_dense, &
        cp2k_exchange_transform_mixed_radix_fft, cp2k_exchange_transform_separable_dft, &
        cp2k_prepare_exchange
   USE tblite_integral_native_integrals,ONLY: get_overlap
   USE tblite_stream_plan,              ONLY: cp2k_stream_plan_type
   USE tblite_utils_average,            ONLY: average_id,&
                                              average_type,&
                                              new_average
   USE tblite_wavefunction_type,        ONLY: new_wavefunction,&
                                              wavefunction_type
   USE tblite_xtb_calculator,           ONLY: xtb_calculator

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: collect_exchange_integration

   REAL(wp), PARAMETER :: cutoff = 25.0_wp
   REAL(wp), PARAMETER :: thr = 100*EPSILON(1.0_wp)
   REAL(wp), PARAMETER :: thr1 = 1e5*EPSILON(1.0_wp)
   REAL(wp), PARAMETER :: thr2 = SQRT(EPSILON(1.0_wp))
   INTEGER(int64) :: partial_test_transaction_serial = 0_int64

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param testsuite ...
! **************************************************************************************************
   SUBROUTINE collect_exchange_integration(testsuite)
      TYPE(unittest_type), ALLOCATABLE, INTENT(out)      :: testsuite(:)

      testsuite = [ &
                  new_unittest("bvk-exchange-supercell", test_bvk_exchange_kernel), &
                  new_unittest("kpoint-fock-complex", test_kpoint_fock_complex), &
                  new_unittest("regular-separable-dft", test_regular_separable_dft), &
                  new_unittest("regular-mixed-radix-fft", test_regular_mixed_radix_fft)]
   END SUBROUTINE collect_exchange_integration

! **************************************************************************************************
!> \brief ...
!> \param bas ...
!> \param mol ...
!> \param error ...
! **************************************************************************************************
   SUBROUTINE make_qvszp_basis(bas, mol, error)
      TYPE(qvszp_basis_type), INTENT(out)                :: bas
      TYPE(structure_type), INTENT(in)                   :: mol
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      INTEGER, PARAMETER :: p_nshell(60) = [1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 2&
         , 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3&
         , 3, 3, 3, 3, 3, 2, 3, 3, 4, 4, 4]
      REAL(wp), PARAMETER :: p_qvszp_kq(60) = [0.3095916809_wp, 1.3956973435_wp, 0.6565062866_wp, &
         0.6414725985_wp, 0.2289797408_wp, 0.3364428754_wp, 0.3134961856_wp, 0.4760878318_wp, &
         0.6924639973_wp, 1.1737897030_wp, 0.8349448993_wp, 0.7162960964_wp, 1.4698576331_wp, &
         0.4013087493_wp, 0.3606081146_wp, 0.3831440321_wp, 0.4730080776_wp, 0.6098840196_wp, &
         0.2870475794_wp, 1.2095850805_wp, 0.5657524714_wp, 0.2998857159_wp, 0.4606086678_wp, &
         1.0397549071_wp, 0.4977801340_wp, 0.3450018637_wp, 0.6470970913_wp, 0.0538554545_wp, &
         0.3313819444_wp, 0.1887280360_wp, 0.3578752499_wp, 0.0778016361_wp, 0.4943243573_wp, &
         0.3404520755_wp, 0.4308339170_wp, 0.6103851192_wp, 0.6281945870_wp, 0.4921139609_wp, &
         0.4926232757_wp, 0.0321070159_wp, -0.0305531270_wp, 0.6988874741_wp, 0.2360577501_wp, &
         0.3463868970_wp, 0.5063045934_wp, 0.1909760796_wp, 0.0011218576_wp, 0.7665007573_wp, &
         0.8084270695_wp, 0.7971344741_wp, 0.0628724514_wp, 0.3587445765_wp, 0.4782877821_wp, &
         0.3372686988_wp, 0.2329265134_wp, 0.5416188137_wp, 0.1911360071_wp, 0.2615471466_wp, &
         0.4000000000_wp, 0.4000000000_wp]

      INTEGER                                            :: ish, isp, izp
      INTEGER, ALLOCATABLE                               :: nshell(:)
      TYPE(cgto_container), ALLOCATABLE                  :: cgto(:, :), cgto_h0(:, :)
      TYPE(qvszp_cgto_type), ALLOCATABLE                 :: cgto_qvszp

!> Number of shells selected from the q-vSZP basis set
!1-20
!21-40
!41-60
!> Scaling of the atomic charge in q-vSZP
!1-4
!5-8
!9-12
!13-16
!17-20
!21-24
!25-28
!29-32
!33-36
!37-40
!41-44
!45-48
!49-52
!53-56
!57-60

      nshell = p_nshell(mol%num)
      ALLOCATE (cgto(MAXVAL(nshell), mol%nid))
      DO isp = 1, mol%nid
         izp = mol%num(isp)
         DO ish = 1, nshell(isp)
            ALLOCATE (cgto_qvszp)
            CALL new_qvszp_cgto(cgto_qvszp, izp, ish, .TRUE., error, &
               & kq=p_qvszp_kq(izp))
            IF (ALLOCATED(error)) RETURN
            CALL move_alloc(cgto_qvszp, cgto(ish, isp)%raw)
         END DO
      END DO

      CALL new_qvszp_basis(bas, mol, nshell, cgto, error, accuracy=1.0_wp, &
         & cgto_h0=cgto_h0)
      IF (ALLOCATED(error)) RETURN

   END SUBROUTINE make_qvszp_basis

! **************************************************************************************************
!> \brief ...
!> \param exchange ...
!> \param mol ...
!> \param bas ...
! **************************************************************************************************
   SUBROUTINE make_exchange_gxtb(exchange, mol, bas)

      !> New exchange object
      class(exchange_type), ALLOCATABLE, INTENT(out) :: exchange

      !> Molecular structure data
      TYPE(structure_type), INTENT(in) :: mol

      !> Basis set information
      class(basis_type), INTENT(in) :: bas

      !> Atomic Hubbard parameters or chemical hardness
      REAL(wp), PARAMETER :: hubbard_parameter(60) = [&
         &  0.4725928800_wp, 0.9220339100_wp, 0.1745288800_wp, 0.2570073300_wp, & !1-4
         &  0.3394908600_wp, 0.4219541200_wp, 0.5043819300_wp, 0.5869186300_wp, & !5-8
         &  0.6693135100_wp, 0.7519160700_wp, 0.1796410500_wp, 0.2215727600_wp, & !9-12
         &  0.2634857800_wp, 0.3053964500_wp, 0.3473401400_wp, 0.3892472500_wp, & !13-16
         &  0.4311567000_wp, 0.4730826900_wp, 0.1710546900_wp, 0.2027624400_wp, & !17-20
         &  0.2100732200_wp, 0.2173964700_wp, 0.2247103900_wp, 0.2320150100_wp, & !21-24
         &  0.2393396900_wp, 0.2466563800_wp, 0.2539825500_wp, 0.2612886300_wp, & !25-28
         &  0.2685947600_wp, 0.2759256500_wp, 0.3076299900_wp, 0.3393158000_wp, & !29-32
         &  0.3723598500_wp, 0.4027354900_wp, 0.4344577600_wp, 0.4661170800_wp, & !33-36
         &  0.1558507900_wp, 0.1864932400_wp, 0.1935621000_wp, 0.2006331100_wp, & !37-40
         &  0.2077052200_wp, 0.2147725400_wp, 0.2218461400_wp, 0.2289187200_wp, & !41-44
         &  0.2359862100_wp, 0.2430561200_wp, 0.2501301800_wp, 0.2571993700_wp, & !45-48
         &  0.2878478000_wp, 0.3184867300_wp, 0.3491243100_wp, 0.3797659300_wp, & !49-52
         &  0.4104080800_wp, 0.4410577700_wp, 0.0501933200_wp, 0.0676257000_wp, & !53-56
         &  0.0850444500_wp, 0.1024773600_wp, 0.1199110500_wp, 0.1373277200_wp] !57-60

      REAL(wp), PARAMETER :: fock_shell_hubbard(4, 60) = RESHAPE([&
         &  3.0691283251_wp, 0.0000000000_wp, 0.0000000000_wp, 0.0000000000_wp, & !1
         &  5.6874061138_wp, 0.0000000000_wp, 0.0000000000_wp, 0.0000000000_wp, & !2
         &  2.9883012309_wp, 2.3353110277_wp, 0.0000000000_wp, 0.0000000000_wp, & !3
         &  2.7369793591_wp, 4.6049991892_wp, 0.0000000000_wp, 0.0000000000_wp, & !4
         &  2.3264777208_wp, 3.0982127996_wp, 0.0000000000_wp, 0.0000000000_wp, & !5
         &  2.7744401520_wp, 2.6734099212_wp, 0.0000000000_wp, 0.0000000000_wp, & !6
         &  4.9809093168_wp, 2.5461337307_wp, 0.0000000000_wp, 0.0000000000_wp, & !7
         &  5.8776683336_wp, 1.8560619141_wp, 0.0000000000_wp, 0.0000000000_wp, & !8
         &  5.7439606016_wp, 1.8335282134_wp, 0.0000000000_wp, 0.0000000000_wp, & !9
         &  4.5041390432_wp, 2.2505127692_wp, 0.0000000000_wp, 0.0000000000_wp, & !10
         &  5.1848157487_wp, 3.3581988954_wp, 0.0000000000_wp, 0.0000000000_wp, & !11
         &  2.5240458726_wp, 3.3034249243_wp, 0.0500000000_wp, 0.0000000000_wp, & !12
         &  2.8801279040_wp, 2.9677545399_wp, 0.0532061577_wp, 0.0000000000_wp, & !13
         &  4.7017489845_wp, 3.2451166892_wp, 0.0281128786_wp, 0.0000000000_wp, & !14
         &  5.5892206172_wp, 2.7892566588_wp, 0.0400000000_wp, 0.0000000000_wp, & !15
         &  4.0689142973_wp, 2.3586924706_wp, 0.0190557157_wp, 0.0000000000_wp, & !16
         &  4.0405582899_wp, 2.7444650014_wp, 0.1372215230_wp, 0.0000000000_wp, & !17
         &  3.8099727932_wp, 1.9948057658_wp, 0.3816448105_wp, 0.0000000000_wp, & !18
         &  2.8076441296_wp, 7.3193704755_wp, 0.0000000000_wp, 0.0000000000_wp, & !19
         &  1.4169946183_wp, 2.8506470516_wp, 0.3652352077_wp, 0.0000000000_wp, & !20
         &  3.5613007676_wp, 0.2140174980_wp, 6.0777797114_wp, 0.0000000000_wp, & !21
         &  4.3359867599_wp, 0.4377995490_wp, 5.6339417301_wp, 0.0000000000_wp, & !22
         &  5.1611082663_wp, 0.1493426752_wp, 5.6817801755_wp, 0.0000000000_wp, & !23
         &  3.6190204568_wp, 0.3459252864_wp, 5.1104407513_wp, 0.0000000000_wp, & !24
         &  4.8952268471_wp, 0.4027261153_wp, 5.4351604831_wp, 0.0000000000_wp, & !25
         &  2.9615379891_wp, 0.9965207694_wp, 5.7517572774_wp, 0.0000000000_wp, & !26
         &  2.6929830338_wp, 2.0126619053_wp, 5.3053752298_wp, 0.0000000000_wp, & !27
         &  1.6560700161_wp, 0.6816562078_wp, 5.7265054452_wp, 0.0000000000_wp, & !28
         &  2.6685365919_wp, 0.2648431263_wp, 5.9083571181_wp, 0.0000000000_wp, & !29
         &  3.2723488550_wp, 0.2082148600_wp, 0.0000000000_wp, 0.0000000000_wp, & !30
         &  3.9311761748_wp, 3.5498970279_wp, 0.0500000000_wp, 0.0000000000_wp, & !31
         &  4.0015034453_wp, 2.8636412859_wp, 0.0500000000_wp, 0.0000000000_wp, & !32
         &  3.9383061872_wp, 2.2782579906_wp, 0.0500000000_wp, 0.0000000000_wp, & !33
         &  4.1143033730_wp, 2.0990937392_wp, 0.0500000000_wp, 0.0000000000_wp, & !34
         &  2.6787903609_wp, 2.3185760780_wp, 0.0509775437_wp, 0.0000000000_wp, & !35
         &  3.9635929046_wp, 1.9680253501_wp, 0.0641872392_wp, 0.0000000000_wp, & !36
         &  3.6554521685_wp, 5.4218452070_wp, 0.0000000000_wp, 0.0000000000_wp, & !37
         &  1.2386896827_wp, 0.0367532240_wp, 0.5849934473_wp, 0.0000000000_wp, & !38
         &  2.2661820105_wp, 0.3529903506_wp, 5.7538524045_wp, 0.0000000000_wp, & !39
         &  2.7125214012_wp, 0.9174413894_wp, 5.2357797226_wp, 0.0000000000_wp, & !40
         &  4.2181708986_wp, 0.2363754400_wp, 5.6694352360_wp, 0.0000000000_wp, & !41
         &  5.1555524951_wp, 0.3545348751_wp, 4.1082662773_wp, 0.0000000000_wp, & !42
         &  3.1512449997_wp, 0.1369894568_wp, 4.3933680398_wp, 0.0000000000_wp, & !43
         &  4.4480817529_wp, 1.2580149069_wp, 5.4295548515_wp, 0.0000000000_wp, & !44
         &  4.2971074517_wp, 1.0733384604_wp, 4.4840736158_wp, 0.0000000000_wp, & !45
         &  1.6735108374_wp, 0.8251776679_wp, 3.8684970093_wp, 0.0000000000_wp, & !46
         &  4.2280578634_wp, 0.1699440689_wp, 4.4140710984_wp, 0.0000000000_wp, & !47
         &  3.4063345215_wp, 0.1136195028_wp, 0.0000000000_wp, 0.0000000000_wp, & !48
         &  4.4529575665_wp, 3.5071016244_wp, 0.0865944553_wp, 0.0000000000_wp, & !49
         &  4.0663732807_wp, 2.6287679185_wp, 0.1582429583_wp, 0.0000000000_wp, & !50
         &  4.8767168113_wp, 2.7385818502_wp, 0.0444124760_wp, 0.0000000000_wp, & !51
         &  3.0994063304_wp, 1.8384995197_wp, 0.0457066802_wp, 0.0000000000_wp, & !52
         &  3.2114763510_wp, 2.2737284734_wp, 0.0600000000_wp, 0.0000000000_wp, & !53
         &  3.5544483922_wp, 1.5687665113_wp, -0.0276663624_wp, 0.0000000000_wp, & !54
         &  3.9806565428_wp, 12.1869253675_wp, 0.0000000000_wp, 0.0000000000_wp, & !55
         &  1.8200389980_wp, 0.7112171196_wp, 0.9568227551_wp, 0.0000000000_wp, & !56
         &  5.2212490716_wp, 0.0615487319_wp, 7.8847923351_wp, 0.0000000000_wp, & !57
         &  3.8139401615_wp, 0.4064116095_wp, 3.5473487464_wp, 9.5959434233_wp, & !58
         &  4.1859723893_wp, 0.3661069834_wp, 2.2803077637_wp, 8.6445128775_wp, & !59
         &  3.8143885126_wp, 1.6445614456_wp, 2.2972917570_wp, 10.6325961873_wp],& !60
         & SHAPE(fock_shell_hubbard))

      REAL(wp), PARAMETER :: fock_avg_exp(4, 60) = RESHAPE([&
         &  1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, & !1
         &  1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, & !2
         &  1.0000000000_wp, 2.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, & !3
         &  1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, & !4
         &  1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, & !5
         &  1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, & !6
         &  1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, & !7
         &  1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, & !8
         &  1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, & !9
         &  1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, & !10
         &  1.0000000000_wp, 2.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, & !11
         &  1.0000000000_wp, 1.0000000000_wp, 2.0000000000_wp, 1.0000000000_wp, & !12
         &  1.0000000000_wp, 1.0000000000_wp, 2.0000000000_wp, 1.0000000000_wp, & !13
         &  1.0000000000_wp, 1.0000000000_wp, 2.0000000000_wp, 1.0000000000_wp, & !14
         &  1.0000000000_wp, 1.0000000000_wp, 2.0000000000_wp, 1.0000000000_wp, & !15
         &  1.0000000000_wp, 1.0000000000_wp, 2.0000000000_wp, 1.0000000000_wp, & !16
         &  1.0000000000_wp, 1.0000000000_wp, 2.0000000000_wp, 1.0000000000_wp, & !17
         &  1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, & !18
         &  1.0000000000_wp, 2.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, & !19
         &  1.0000000000_wp, 1.0000000000_wp, 2.0000000000_wp, 1.0000000000_wp, & !20
         &  1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, & !21
         &  1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, & !22
         &  1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, & !23
         &  1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, & !24
         &  1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, & !25
         &  1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, & !26
         &  1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, & !27
         &  1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, & !28
         &  1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, & !29
         &  1.0000000000_wp, 2.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, & !30
         &  1.0000000000_wp, 1.0000000000_wp, 2.0000000000_wp, 1.0000000000_wp, & !31
         &  1.0000000000_wp, 1.0000000000_wp, 2.0000000000_wp, 1.0000000000_wp, & !32
         &  1.0000000000_wp, 1.0000000000_wp, 2.0000000000_wp, 1.0000000000_wp, & !33
         &  1.0000000000_wp, 1.0000000000_wp, 2.0000000000_wp, 1.0000000000_wp, & !34
         &  1.0000000000_wp, 1.0000000000_wp, 2.0000000000_wp, 1.0000000000_wp, & !35
         &  1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, & !36
         &  1.0000000000_wp, 2.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, & !37
         &  1.0000000000_wp, 1.0000000000_wp, 2.0000000000_wp, 1.0000000000_wp, & !38
         &  1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, & !39
         &  1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, & !40
         &  1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, & !41
         &  1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, & !42
         &  1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, & !43
         &  1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, & !44
         &  1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, & !45
         &  1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, & !46
         &  1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, & !47
         &  1.0000000000_wp, 2.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, & !48
         &  1.0000000000_wp, 1.0000000000_wp, 2.0000000000_wp, 1.0000000000_wp, & !49
         &  1.0000000000_wp, 1.0000000000_wp, 2.0000000000_wp, 1.0000000000_wp, & !50
         &  1.0000000000_wp, 1.0000000000_wp, 2.0000000000_wp, 1.0000000000_wp, & !51
         &  1.0000000000_wp, 1.0000000000_wp, 2.0000000000_wp, 1.0000000000_wp, & !52
         &  1.0000000000_wp, 1.0000000000_wp, 2.0000000000_wp, 1.0000000000_wp, & !53
         &  1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, & !54
         &  1.0000000000_wp, 2.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, & !55
         &  1.0000000000_wp, 1.0000000000_wp, 2.0000000000_wp, 1.0000000000_wp, & !56
         &  1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, & !57
         &  1.0000000000_wp, 2.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, & !58
         &  1.0000000000_wp, 2.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, & !59
         &  1.0000000000_wp, 2.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp],& !60
         & SHAPE(fock_avg_exp))

      !> Parameter: Correlation radius for bond-order correlation correction to the Fock exchange
      REAL(wp), PARAMETER :: fock_corr_rad(60) = [&
         &  7.7792324405_wp, 25.0000000000_wp, 8.8096512303_wp, 5.0976031817_wp, & !1-4
         & 13.9264964630_wp, 16.1500544221_wp, 4.5816631552_wp, 4.4969668090_wp, & !5-8
         &  2.6478686663_wp, 1.3077199954_wp, 9.1548600312_wp, 8.8377232712_wp, & !9-12
         &  4.8249049725_wp, 11.0887979395_wp, 39.4736123213_wp, 12.5349113695_wp, & !13-16
         &  9.8350329594_wp, 6.1005201925_wp, 10.9559357706_wp, 11.7243327664_wp, & !17-20
         &  6.6159933775_wp, 12.2871886821_wp, 6.5869513804_wp, 5.6059999728_wp, & !21-24
         &  5.1232611294_wp, 7.1831531782_wp, 8.1760551383_wp, 5.7650967598_wp, & !25-28
         &  6.9475450757_wp, 6.8197944590_wp, 10.9433467644_wp, 12.2246936853_wp, & !29-32
         & 16.5095199438_wp, 8.5848893693_wp, 9.2732288842_wp, 6.3347878861_wp, & !33-36
         & 12.1337981396_wp, 8.4248544787_wp, 9.0869814898_wp, 7.1952756180_wp, & !37-40
         & 12.1038392588_wp, 21.7953794701_wp, 5.0119884301_wp, 15.5776135512_wp, & !41-44
         &  5.0766070904_wp, 15.7440179356_wp, 8.0809694075_wp, 10.5430243362_wp, & !45-48
         & 20.8650121649_wp, 15.9260140254_wp, 19.8076868475_wp, 28.4426758250_wp, & !49-52
         & 11.2492941795_wp, 7.3164432561_wp, 15.0000000000_wp, 11.9865176851_wp, & !53-56
         &  9.5686567782_wp, 2.6257749318_wp, 3.1796915366_wp, 4.9047851530_wp] !57-60

      !> Parameter: Scaling factor for bond-order correlation correction to the Fock exchange
      REAL(wp), PARAMETER :: fock_corr_scale(60) = [&
         &  0.0155118215_wp, 0.0017675886_wp, 0.0034145054_wp, 0.0188792718_wp, & !1-4
         &  0.0110877319_wp, 0.0148694859_wp, 0.0519346235_wp, 0.0207320201_wp, & !5-8
         &  0.0458054593_wp, 0.4683557610_wp, 0.0220314012_wp, 0.0131005859_wp, & !9-12
         &  0.0161186017_wp, 0.0242053787_wp, 0.0174007860_wp, 0.0226318766_wp, & !13-16
         &  0.0098241753_wp, 0.0402121335_wp, 0.0497087513_wp, 0.0015206165_wp, & !17-20
         &  0.0017029584_wp, 0.0030354812_wp, 0.0140465231_wp, 0.0120715398_wp, & !21-24
         &  0.0089442952_wp, 0.0047486806_wp, 0.0064777302_wp, 0.0027572746_wp, & !25-28
         &  0.0043159157_wp, 0.0092918050_wp, 0.0168094876_wp, 0.0146359776_wp, & !29-32
         &  0.0140255958_wp, 0.0248514097_wp, 0.0135186782_wp, 0.0174438731_wp, & !33-36
         &  0.0321525083_wp, 0.0032846137_wp, 0.0109494264_wp, 0.0015762687_wp, & !37-40
         &  0.0018844293_wp, 0.0108340538_wp, 0.0053222684_wp, 0.0015263334_wp, & !41-44
         &  0.0058884094_wp, 0.0075970636_wp, 0.0064299028_wp, 0.0124120357_wp, & !45-48
         &  0.0113170904_wp, 0.0158248636_wp, 0.0195890148_wp, 0.0202387721_wp, & !49-52
         &  0.0132956936_wp, 0.0235755901_wp, 0.0347877982_wp, 0.0048909086_wp, & !53-56
         &  0.0008293725_wp, 0.0085376871_wp, 0.0034042669_wp, 0.0089899347_wp] !57-60

      !> Parameter: van-der-Waals radii scaling
      REAL(wp), PARAMETER :: rvdw_scale(60) = [&
         &  0.9950000000_wp, 1.1000000000_wp, 1.0000000000_wp, 1.0200000000_wp, & !1-4
         &  1.0200000000_wp, 1.0000000000_wp, 1.0500000000_wp, 1.0200000000_wp, & !5-8
         &  1.1800000000_wp, 0.8500000000_wp, 1.0500000000_wp, 1.1500000000_wp, & !9-12
         &  1.0400000000_wp, 1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, & !13-16
         &  1.0000000000_wp, 1.0000000000_wp, 0.9500000000_wp, 0.9500000000_wp, & !17-20
         &  1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, & !21-24
         &  1.0000000000_wp, 1.1000000000_wp, 1.1000000000_wp, 1.1000000000_wp, & !25-28
         &  1.0000000000_wp, 1.1000000000_wp, 1.0000000000_wp, 1.0000000000_wp, & !29-32
         &  1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, & !33-36
         &  1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, & !37-40
         &  1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, & !41-44
         &  1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, 1.1500000000_wp, & !45-48
         &  1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, 1.0000000000_wp, & !49-52
         &  1.0000000000_wp, 1.0000000000_wp, 0.5000000000_wp, 1.0000000000_wp, & !53-56
         &  1.2000000000_wp, 1.2000000000_wp, 1.2000000000_wp, 1.2000000000_wp] !57-60

      REAL(wp), PARAMETER :: p_kq_shell(0:3) = &
         & [1.0000000000_wp, 0.5000000000_wp, 0.2500000000_wp, 0.1250000000_wp]
      REAL(wp), PARAMETER :: p_kq = 2.1400000000_wp
      REAL(wp), PARAMETER :: p_ondiag = 1.4017990616_wp
      REAL(wp), PARAMETER :: p_offdiag_l(0:3) = &
         & [0.0500000000_wp, 2.5000000000_wp, 2.5000000000_wp, 2.5000000000_wp]
      REAL(wp), PARAMETER :: p_hubbard_exp = -0.3296069670_wp
      REAL(wp), PARAMETER :: p_hubbard_exp_r0 = 0.0128349272_wp*autoaa
      REAL(wp), PARAMETER :: p_omega = 0.2000000000_wp
      REAL(wp), PARAMETER :: p_a0 = 0.1500000000_wp
      REAL(wp), PARAMETER :: p_a0_onsite = 0.1820686001_wp
      REAL(wp), PARAMETER :: p_corr_exp = 1.2000000000_wp/autoaa
      REAL(wp), PARAMETER :: p_corr_scale_avg_exponent = 0.7000000000_wp

      INTEGER :: izp, jzp, isp, jsp, ish, jsh, il, jl
      REAL(wp), ALLOCATABLE :: hardness_fx(:, :), offdiag(:, :)
      REAL(wp), ALLOCATABLE :: onecxints(:, :, :), kq(:, :)
      REAL(wp), ALLOCATABLE :: avg_exponents(:, :), rad(:, :)
      REAL(wp), ALLOCATABLE :: corr_scale(:), corr_rad(:)
      TYPE(exchange_fock), ALLOCATABLE :: fock
      TYPE(average_type), ALLOCATABLE :: hubbard_average
      TYPE(average_type), ALLOCATABLE :: offdiag_average
      TYPE(average_type), ALLOCATABLE :: corr_scale_average
      TYPE(average_type), ALLOCATABLE :: corr_rad_average
      TYPE(average_type), ALLOCATABLE :: rvdw_average

      ! Obtain radii for radius-dependent hubbard scaling
      ALLOCATE (rvdw_average)
      CALL new_average(rvdw_average, average_id%arithmetic)
      ALLOCATE (rad(mol%nid, mol%nid))
      DO isp = 1, mol%nid
         izp = mol%num(isp)
         DO jsp = 1, isp
            jzp = mol%num(jsp)
            rad(jsp, isp) = get_vdw_rad(jzp, izp) &
                           & *rvdw_average%value(rvdw_scale(izp), rvdw_scale(jzp))
            rad(isp, jsp) = rad(jsp, isp)
         END DO
      END DO

      ! Obtain the shell Hubbard parameters
      ALLOCATE (hardness_fx(MAXVAL(bas%nsh_id), mol%nid))
      hardness_fx(:, :) = 0.0_wp
      DO isp = 1, mol%nid
         izp = mol%num(isp)
         DO ish = 1, bas%nsh_id(isp)
            hardness_fx(ish, isp) = hubbard_parameter(izp) &
                                   & *fock_shell_hubbard(ish, izp)
         END DO
      END DO

      ! Set averager for hubbard and correlation parameters
      ALLOCATE (hubbard_average, offdiag_average, corr_scale_average, corr_rad_average)
      CALL new_average(hubbard_average, average_id%general)
      CALL new_average(offdiag_average, average_id%geometric)
      CALL new_average(corr_scale_average, average_id%general)
      CALL new_average(corr_rad_average, average_id%harmonic)

      ! Obtain the averaging exponents for the shell Hubbard parameters
      ALLOCATE (avg_exponents(MAXVAL(bas%nsh_id), mol%nid))
      avg_exponents(:, :) = 0.0_wp
      DO isp = 1, mol%nid
         izp = mol%num(isp)
         DO ish = 1, bas%nsh_id(isp)
            avg_exponents(ish, isp) = fock_avg_exp(ish, izp)
         END DO
      END DO

      ! Obtain the one-center exchange integrals
      ALLOCATE (onecxints(MAXVAL(bas%nsh_id), MAXVAL(bas%nsh_id), mol%nid))
      onecxints(:, :, :) = 0.0_wp
      DO isp = 1, mol%nid
         izp = mol%num(isp)
         DO ish = 1, bas%nsh_id(isp)
            il = bas%cgto(ish, isp)%raw%ang
            DO jsh = 1, bas%nsh_id(isp)
               jl = bas%cgto(jsh, isp)%raw%ang
               onecxints(jsh, ish, isp) = get_onecxints(jl, il, izp)
            END DO
         END DO
      END DO

      ! Obtain the off-diagaonal scaling exchange factors
      ALLOCATE (offdiag(MAXVAL(bas%nsh_id), mol%nid))
      offdiag(:, :) = 0.0_wp
      DO isp = 1, mol%nid
         izp = mol%num(isp)
         DO ish = 1, bas%nsh_id(isp)
            il = bas%cgto(ish, isp)%raw%ang
            offdiag(ish, isp) = p_offdiag_l(il)
         END DO
      END DO

      ! Obtain the charge scaling of the exchange fraction
      ALLOCATE (kq(MAXVAL(bas%nsh_id), mol%nid))
      kq(:, :) = 0.0_wp
      DO isp = 1, mol%nid
         izp = mol%num(isp)
         DO ish = 1, bas%nsh_id(isp)
            il = bas%cgto(ish, isp)%raw%ang
            kq(ish, isp) = SQRT(p_kq)*p_kq_shell(il)
         END DO
      END DO

      ! Obtain the correlation scaling factor and radii
      corr_scale = fock_corr_scale(mol%num)
      corr_rad = fock_corr_rad(mol%num)

      ALLOCATE (fock)
      CALL new_exchange_fock(fock, mol, bas, hardness_fx, hubbard_average, &
         & avg_exponents, p_ondiag, offdiag, offdiag_average, p_hubbard_exp, &
         & p_hubbard_exp_r0, rad, 1.0_wp, onecxints, kq, corr_scale, &
         & corr_scale_average, p_corr_scale_avg_exponent, p_corr_exp, corr_rad, &
         & corr_rad_average, p_a0, p_a0_onsite, p_omega, 1.0_wp - p_a0)
      CALL move_alloc(fock, exchange)

   END SUBROUTINE make_exchange_gxtb

! **************************************************************************************************
!> \brief ...
!> \param cache ...
!> \param ptr ...
! **************************************************************************************************
   SUBROUTINE view(cache, ptr)
      !> Instance of the container cache
      TYPE(container_cache), INTENT(inout), TARGET       :: cache
      TYPE(exchange_cache), INTENT(out), POINTER         :: ptr

!> Reference to the container cache

      NULLIFY (ptr)
      SELECT TYPE (TARGET => cache%raw)
      TYPE is (exchange_cache)
         ptr => TARGET
      END SELECT
   END SUBROUTINE view

! **************************************************************************************************
!> \brief ...
!> \param transaction_id ...
! **************************************************************************************************
   SUBROUTINE next_partial_test_transaction(transaction_id)
      INTEGER(int64), INTENT(out)                        :: transaction_id(2)

      partial_test_transaction_serial = partial_test_transaction_serial + 1_int64
      transaction_id = [3141592653589793_int64, partial_test_transaction_serial]
   END SUBROUTINE next_partial_test_transaction

! **************************************************************************************************
!> \brief ...
!> \param descriptor ...
!> \param error ...
!> \param reject_bad_abi ...
! **************************************************************************************************
   SUBROUTINE check_partial_descriptor_roundtrip(descriptor, error, reject_bad_abi)
      TYPE(cp2k_exchange_partial_descriptor_type), &
         INTENT(in)                                      :: descriptor
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error
      LOGICAL, INTENT(in), OPTIONAL                      :: reject_bad_abi

      INTEGER(int64), ALLOCATABLE                        :: words(:)
      TYPE(cp2k_exchange_partial_descriptor_type)        :: roundtrip

      CALL cp2k_exchange_partial_descriptor_pack(descriptor, words, error)
      IF (ALLOCATED(error)) RETURN
      IF (INT(SIZE(words), int64) /= cp2k_exchange_partial_descriptor_size()) THEN
         CALL test_failed(error, "Incorrect canonical partial descriptor size")
         RETURN
      END IF
      CALL cp2k_exchange_partial_descriptor_unpack(roundtrip, words, error)
      IF (ALLOCATED(error)) RETURN
      IF (.NOT. cp2k_exchange_partial_descriptors_match(descriptor, roundtrip)) THEN
         CALL test_failed(error, "Canonical partial descriptor roundtrip changed")
         RETURN
      END IF

      IF (PRESENT(reject_bad_abi)) THEN
         IF (reject_bad_abi) THEN
            words(1) = words(1) + 1_int64
            CALL cp2k_exchange_partial_descriptor_unpack(roundtrip, words, error)
            IF (.NOT. ALLOCATED(error)) THEN
               CALL test_failed(error, &
                  & "Malformed canonical partial descriptor was accepted")
               RETURN
            END IF
            DEALLOCATE (error)
         END IF
      END IF
   END SUBROUTINE check_partial_descriptor_roundtrip

! **************************************************************************************************
!> \brief ...
!> \param calc ...
!> \param mol ...
!> \param cache ...
!> \param cache_plan ...
!> \param wfn ...
!> \param nmesh ...
!> \param kfrac ...
!> \param weights ...
!> \param density ...
!> \param overlap ...
!> \param partition_counts ...
!> \param run_negative ...
!> \param error ...
!> \param maxerr ...
! **************************************************************************************************
   SUBROUTINE check_partial_exchange_case(calc, mol, cache, cache_plan, wfn, nmesh, kfrac, &
      & weights, density, overlap, partition_counts, run_negative, error, maxerr)
      TYPE(xtb_calculator), INTENT(inout)                :: calc
      TYPE(structure_type), INTENT(in)                   :: mol
      TYPE(container_cache), INTENT(inout)               :: cache
      TYPE(cp2k_stream_plan_type), INTENT(inout)         :: cache_plan
      TYPE(wavefunction_type), INTENT(in)                :: wfn
      INTEGER, INTENT(in)                                :: nmesh(3)
      REAL(wp), INTENT(in)                               :: kfrac(:, :), weights(:)
      COMPLEX(wp), INTENT(in)                            :: density(:, :, :, :), overlap(:, :, :)
      INTEGER, INTENT(in)                                :: partition_counts(:)
      LOGICAL, INTENT(in)                                :: run_negative
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error
      REAL(wp), INTENT(out)                              :: maxerr

      COMPLEX(wp), ALLOCATABLE :: density_trial(:, :, :), fock_partial(:, :, :, :), &
         fock_ref(:, :, :, :), forward_raw(:), forward_sum(:), overlap_adjoint(:, :, :), &
         overlap_adjoint_ref(:, :, :), reverse_raw(:), reverse_sum(:)
      INTEGER                                            :: ik, ipart, ncomplex, nk, npart
      INTEGER(int64)                                     :: mismatch_transaction_id(2), &
                                                            transaction_id(2)
      INTEGER(int64), ALLOCATABLE                        :: descriptor_words(:)
      INTEGER, ALLOCATABLE                               :: owner_count(:), owner_matrix(:, :), &
                                                            owner_sum(:), owner_trial(:)
      LOGICAL                                            :: done
      REAL(wp)                                           :: energy_partial, energy_ref
      REAL(wp), ALLOCATABLE :: gradient_partial(:, :), gradient_ref(:, :), sigma_partial(:, :), &
         sigma_ref(:, :), vsh_partial(:), vsh_ref(:)
      TYPE(cp2k_exchange_partial_descriptor_type)        :: mismatch_descriptor
      TYPE(cp2k_exchange_partial_descriptor_type), &
         ALLOCATABLE                                     :: descriptor(:)
      TYPE(cp2k_exchange_partial_forward_type)           :: forward_sum_buffer, &
                                                            mismatch_forward_buffer
      TYPE(cp2k_exchange_partial_forward_type), &
         ALLOCATABLE                                     :: forward_buffer(:)
      TYPE(cp2k_exchange_partial_reverse_type)           :: reverse_sum_buffer
      TYPE(cp2k_exchange_partial_reverse_type), &
         ALLOCATABLE                                     :: reverse_buffer(:)
      TYPE(cp2k_exchange_stream_type)                    :: mismatch_stream
      TYPE(cp2k_exchange_stream_type), ALLOCATABLE       :: stream(:)

      nk = SIZE(weights)
      ALLOCATE (fock_ref(SIZE(density, 1), SIZE(density, 2), &
         & SIZE(density, 3), nk), vsh_ref(SIZE(wfn%qsh, 1)), &
         & overlap_adjoint_ref(SIZE(overlap, 1), SIZE(overlap, 2), nk), &
         & gradient_ref(3, mol%nat), sigma_ref(3, 3))
      CALL cp2k_exchange_kmesh(calc, mol, cache, cache_plan, wfn, nmesh, kfrac, weights, &
         & density, overlap, fock_ref, vsh_ref, energy_ref, error)
      IF (ALLOCATED(error)) RETURN
      CALL cp2k_exchange_kmesh_gradient(calc, mol, cache, cache_plan, nmesh, kfrac, weights, &
         & density, overlap, overlap_adjoint_ref, gradient_ref, sigma_ref, error)
      IF (ALLOCATED(error)) RETURN

      maxerr = 0.0_wp
      DO ipart = 1, SIZE(partition_counts)
         npart = partition_counts(ipart)
         ALLOCATE (stream(npart), descriptor(npart), forward_buffer(npart), &
            & owner_matrix(nk, npart))
         owner_matrix = 0
         CALL next_partial_test_transaction(transaction_id)
         DO ik = 1, npart
            CALL cp2k_exchange_partial_begin(stream(ik), calc, mol, cache, cache_plan, wfn, &
               & nmesh, kfrac, weights, nk, transaction_id, error)
            IF (ALLOCATED(error)) RETURN
         END DO
         DO ik = nk, 1, -1
            ncomplex = 1 + MODULO(ik - 1, npart)
            CALL cp2k_exchange_partial_push(stream(ncomplex), ik, &
               & density(:, :, :, ik), overlap(:, :, ik), error)
            IF (ALLOCATED(error)) RETURN
         END DO

         IF (run_negative .AND. npart == 2) THEN
            CALL cp2k_exchange_partial_push(stream(1), 1, density(:, :, :, 1), &
               & overlap(:, :, 1), error)
            IF (.NOT. ALLOCATED(error)) THEN
               CALL test_failed(error, "Duplicate partial ownership was accepted")
               RETURN
            END IF
            DEALLOCATE (error)
            CALL cp2k_exchange_partial_push(stream(1), nk + 1, &
               & density(:, :, :, 1), overlap(:, :, 1), error)
            IF (.NOT. ALLOCATED(error)) THEN
               CALL test_failed(error, "Out-of-range partial ownership was accepted")
               RETURN
            END IF
            DEALLOCATE (error)
            CALL cp2k_exchange_stream_batch_apply(stream(1), calc, mol, cache, cache_plan, &
               & error)
            IF (.NOT. ALLOCATED(error)) THEN
               CALL test_failed(error, "Partial kernel apply before import succeeded")
               RETURN
            END IF
            DEALLOCATE (error)
         END IF

         DO ik = 1, npart
            CALL cp2k_exchange_partial_export_forward(stream(ik), &
               & descriptor(ik), owner_count, forward_buffer(ik), error)
            IF (ALLOCATED(error)) RETURN
            CALL check_partial_descriptor_roundtrip(descriptor(ik), error, &
               & run_negative .AND. npart == 2 .AND. ik == 1)
            IF (ALLOCATED(error)) RETURN
            owner_matrix(:, ik) = owner_count
            IF (ik > 1) THEN
               IF (.NOT. cp2k_exchange_partial_descriptors_match( &
                  & descriptor(1), descriptor(ik))) THEN
                  CALL test_failed(error, "Matching partial descriptors differ")
                  RETURN
               END IF
            END IF
            CALL cp2k_exchange_partial_forward_pack(forward_buffer(ik), &
               & forward_raw, error)
            IF (ALLOCATED(error)) RETURN
            IF (ik == 1) THEN
               ALLOCATE (forward_sum(SIZE(forward_raw)), &
                  & source=(0.0_wp, 0.0_wp))
            END IF
            forward_sum = forward_sum + forward_raw
         END DO
         owner_sum = SUM(owner_matrix, dim=2)
         ncomplex = 4*SIZE(density, 1)*SIZE(density, 2)*nk*SIZE(density, 3)
         IF (cp2k_exchange_partial_forward_size(forward_buffer(1)) /= &
            & INT(ncomplex, int64)) THEN
            CALL test_failed(error, "Incorrect partial forward memory accounting")
            RETURN
         END IF
         CALL cp2k_exchange_partial_forward_unpack(forward_sum_buffer, &
            & descriptor(1), forward_sum, error)
         IF (ALLOCATED(error)) RETURN

         ! Simulate P>Nk: an empty fourth rank must still carry the exact same
         ! descriptor.  A divergent transaction identity is rejected before the
         ! summed payload can enter the live rank-zero accumulator.
         IF (run_negative .AND. npart == 4) THEN
            CALL next_partial_test_transaction(mismatch_transaction_id)
            CALL cp2k_exchange_partial_begin(mismatch_stream, calc, mol, cache, cache_plan, &
               & wfn, nmesh, kfrac, weights, nk, mismatch_transaction_id, error)
            IF (ALLOCATED(error)) RETURN
            CALL cp2k_exchange_partial_export_forward(mismatch_stream, &
               & mismatch_descriptor, owner_count, mismatch_forward_buffer, error)
            IF (ALLOCATED(error)) RETURN
            IF (ANY(owner_count /= 0)) THEN
               CALL test_failed(error, "The simulated empty rank owns k points")
               RETURN
            END IF
            CALL check_partial_descriptor_roundtrip(mismatch_descriptor, error)
            IF (ALLOCATED(error)) RETURN
            IF (cp2k_exchange_partial_descriptors_match(descriptor(1), &
               & mismatch_descriptor)) THEN
               CALL test_failed(error, &
                  & "Divergent empty-rank partial descriptors compare equal")
               RETURN
            END IF
            CALL cp2k_exchange_partial_import_forward(stream(1), &
               & mismatch_descriptor, owner_sum, forward_sum_buffer, calc, mol, &
               & cache, cache_plan, error)
            IF (.NOT. ALLOCATED(error)) THEN
               CALL test_failed(error, &
                  & "Divergent empty-rank partial descriptor was accepted")
               RETURN
            END IF
            DEALLOCATE (error)
            CALL cp2k_exchange_partial_discard(mismatch_stream, error)
            IF (ALLOCATED(error)) RETURN
         END IF

         IF (run_negative .AND. npart == 2) THEN
            owner_trial = owner_sum
            owner_trial(1) = 0
            CALL cp2k_exchange_partial_import_forward(stream(1), descriptor(1), &
               & owner_trial, forward_sum_buffer, calc, mol, cache, cache_plan, error)
            IF (.NOT. ALLOCATED(error)) THEN
               CALL test_failed(error, "Missing partial ownership was accepted")
               RETURN
            END IF
            DEALLOCATE (error)
            owner_trial = owner_sum
            owner_trial(1) = 2
            CALL cp2k_exchange_partial_import_forward(stream(1), descriptor(1), &
               & owner_trial, forward_sum_buffer, calc, mol, cache, cache_plan, error)
            IF (.NOT. ALLOCATED(error)) THEN
               CALL test_failed(error, "Duplicate partial ownership was accepted")
               RETURN
            END IF
            DEALLOCATE (error)

            ! Malicious but structurally valid large descriptor dimensions must
            ! fail before any default-integer product or allocation can overflow.
            CALL cp2k_exchange_partial_descriptor_pack(descriptor(1), &
               & descriptor_words, error)
            IF (ALLOCATED(error)) RETURN
            descriptor_words(8) = INT(HUGE(0), int64)
            CALL cp2k_exchange_partial_descriptor_unpack(mismatch_descriptor, &
               & descriptor_words, error)
            IF (ALLOCATED(error)) RETURN
            CALL cp2k_exchange_partial_forward_unpack(mismatch_forward_buffer, &
               & mismatch_descriptor, forward_sum, error)
            IF (.NOT. ALLOCATED(error)) THEN
               CALL test_failed(error, &
                  & "Oversized partial payload arithmetic was accepted")
               RETURN
            END IF
            DEALLOCATE (error)

            ! A fresh logical transaction may have the same model, geometry and
            ! mesh but different density.  Its mandatory host identity prevents
            ! replay of an old descriptor and pre-kernel payload.
            CALL next_partial_test_transaction(mismatch_transaction_id)
            CALL cp2k_exchange_partial_begin(mismatch_stream, calc, mol, cache, cache_plan, &
               & wfn, nmesh, kfrac, weights, nk, mismatch_transaction_id, error)
            IF (ALLOCATED(error)) RETURN
            ALLOCATE (density_trial(SIZE(density, 1), SIZE(density, 2), &
               & SIZE(density, 3)))
            DO ik = 1, nk
               density_trial = density(:, :, :, ik)
               IF (ik == 1) density_trial(1, 1, 1) = &
                  & density_trial(1, 1, 1) + CMPLX(1.0e-7_wp, 0.0_wp, wp)
               CALL cp2k_exchange_partial_push(mismatch_stream, ik, &
                  & density_trial, overlap(:, :, ik), error)
               IF (ALLOCATED(error)) RETURN
            END DO
            DEALLOCATE (density_trial)
            CALL cp2k_exchange_partial_export_forward(mismatch_stream, &
               & mismatch_descriptor, owner_count, mismatch_forward_buffer, error)
            IF (ALLOCATED(error)) RETURN
            CALL check_partial_descriptor_roundtrip(mismatch_descriptor, error)
            IF (ALLOCATED(error)) RETURN
            IF (cp2k_exchange_partial_descriptors_match(descriptor(1), &
               & mismatch_descriptor)) THEN
               CALL test_failed(error, &
                  & "Fresh transaction descriptors compare equal")
               RETURN
            END IF
            CALL cp2k_exchange_partial_import_forward(mismatch_stream, &
               & descriptor(1), owner_sum, forward_sum_buffer, calc, mol, cache, cache_plan, &
               & error)
            IF (.NOT. ALLOCATED(error)) THEN
               CALL test_failed(error, &
                  & "A cross-transaction partial payload replay was accepted")
               RETURN
            END IF
            DEALLOCATE (error)
            CALL cp2k_exchange_partial_discard(mismatch_stream, error)
            IF (ALLOCATED(error)) RETURN

            CALL next_partial_test_transaction(mismatch_transaction_id)
            CALL cp2k_exchange_partial_begin(mismatch_stream, calc, mol, cache, cache_plan, &
               & wfn, nmesh, kfrac, weights, nk - 1, mismatch_transaction_id, error, &
               & image_first=1, image_count=nk - 1)
            IF (ALLOCATED(error)) RETURN
            CALL cp2k_exchange_partial_discard(mismatch_stream, error)
            IF (.NOT. ALLOCATED(error)) THEN
               CALL test_failed(error, &
                  & "An unexported partial accumulator was discarded")
               RETURN
            END IF
            DEALLOCATE (error)
            CALL cp2k_exchange_partial_export_forward(mismatch_stream, &
               & mismatch_descriptor, owner_count, forward_buffer(1), error)
            IF (ALLOCATED(error)) RETURN
            CALL cp2k_exchange_partial_import_forward(stream(1), &
               & mismatch_descriptor, owner_sum, forward_sum_buffer, calc, mol, &
               & cache, cache_plan, error)
            IF (.NOT. ALLOCATED(error)) THEN
               CALL test_failed(error, "Mismatched partial descriptor was accepted")
               RETURN
            END IF
            DEALLOCATE (error)
            CALL cp2k_exchange_partial_discard(mismatch_stream, error)
            IF (ALLOCATED(error)) RETURN

            forward_raw = forward_sum
            forward_raw(1) = CMPLX(ieee_value(0.0_wp, ieee_quiet_nan), &
               & 0.0_wp, wp)
            CALL cp2k_exchange_partial_forward_unpack(forward_sum_buffer, &
               & descriptor(1), forward_raw, error)
            IF (.NOT. ALLOCATED(error)) THEN
               CALL test_failed(error, "Non-finite partial payload was accepted")
               RETURN
            END IF
            DEALLOCATE (error)
            CALL cp2k_exchange_partial_forward_unpack(forward_sum_buffer, &
               & descriptor(1), forward_sum, error)
            IF (ALLOCATED(error)) RETURN
         END IF

         CALL cp2k_exchange_partial_import_forward(stream(1), descriptor(1), &
            & owner_sum, forward_sum_buffer, calc, mol, cache, cache_plan, error)
         IF (ALLOCATED(error)) RETURN
         IF (run_negative .AND. npart == 2) THEN
            CALL cp2k_exchange_partial_discard(stream(1), error)
            IF (.NOT. ALLOCATED(error)) THEN
               CALL test_failed(error, &
                  & "An imported partial accumulator was discarded")
               RETURN
            END IF
            DEALLOCATE (error)
            CALL cp2k_exchange_partial_import_forward(stream(1), descriptor(1), &
               & owner_sum, forward_sum_buffer, calc, mol, cache, cache_plan, error)
            IF (.NOT. ALLOCATED(error)) THEN
               CALL test_failed(error, "Second partial import was accepted")
               RETURN
            END IF
            DEALLOCATE (error)
            CALL cp2k_exchange_partial_push(stream(1), 1, density(:, :, :, 1), &
               & overlap(:, :, 1), error)
            IF (.NOT. ALLOCATED(error)) THEN
               CALL test_failed(error, "Partial push after import was accepted")
               RETURN
            END IF
            DEALLOCATE (error)
         END IF

         ALLOCATE (fock_partial(SIZE(fock_ref, 1), SIZE(fock_ref, 2), &
            & SIZE(fock_ref, 3), nk), source=(0.0_wp, 0.0_wp))
         ALLOCATE (vsh_partial(SIZE(vsh_ref)))
         CALL cp2k_exchange_stream_batch_apply(stream(1), calc, mol, cache, cache_plan, error)
         IF (ALLOCATED(error)) RETURN
         DO ik = 1, nk
            CALL cp2k_exchange_stream_batch_pull(stream(1), ik, &
               & overlap(:, :, ik), fock_partial(:, :, :, ik), error)
            IF (ALLOCATED(error)) RETURN
         END DO
         CALL cp2k_exchange_stream_batch_advance(stream(1), done, error)
         IF (ALLOCATED(error)) RETURN
         IF (.NOT. done) THEN
            CALL test_failed(error, "Single full-image partial batch did not finish")
            RETURN
         END IF
         CALL cp2k_exchange_stream_batch_result(stream(1), calc, mol, cache, cache_plan, &
            & vsh_partial, energy_partial, error)
         IF (ALLOCATED(error)) RETURN
         maxerr = MAX(maxerr, MAXVAL(ABS(fock_partial - fock_ref)), &
            & MAXVAL(ABS(vsh_partial - vsh_ref)), ABS(energy_partial - energy_ref))
         CALL cp2k_exchange_stream_end(stream(1), error)
         IF (ALLOCATED(error)) RETURN
         DO ik = 2, npart
            CALL cp2k_exchange_partial_discard(stream(ik), error)
            IF (ALLOCATED(error)) RETURN
         END DO
         DEALLOCATE (stream, descriptor, forward_buffer, owner_matrix, &
            & forward_sum, fock_partial, vsh_partial)

         ALLOCATE (stream(npart), descriptor(npart), reverse_buffer(npart), &
            & owner_matrix(nk, npart))
         owner_matrix = 0
         CALL next_partial_test_transaction(transaction_id)
         DO ik = 1, npart
            CALL cp2k_exchange_partial_begin(stream(ik), calc, mol, cache, cache_plan, wfn, &
               & nmesh, kfrac, weights, nk, transaction_id, error, &
               & reverse_only=.TRUE.)
            IF (ALLOCATED(error)) RETURN
         END DO
         DO ik = nk, 1, -1
            ncomplex = 1 + MODULO(ik - 1, npart)
            CALL cp2k_exchange_partial_reverse_push(stream(ncomplex), ik, &
               & density(:, :, :, ik), overlap(:, :, ik), error)
            IF (ALLOCATED(error)) RETURN
         END DO
         DO ik = 1, npart
            CALL cp2k_exchange_partial_export_reverse(stream(ik), &
               & descriptor(ik), owner_count, reverse_buffer(ik), error)
            IF (ALLOCATED(error)) RETURN
            CALL check_partial_descriptor_roundtrip(descriptor(ik), error)
            IF (ALLOCATED(error)) RETURN
            owner_matrix(:, ik) = owner_count
            IF (ik > 1) THEN
               IF (.NOT. cp2k_exchange_partial_descriptors_match( &
                  & descriptor(1), descriptor(ik))) THEN
                  CALL test_failed(error, "Matching reverse descriptors differ")
                  RETURN
               END IF
            END IF
            CALL cp2k_exchange_partial_reverse_pack(reverse_buffer(ik), &
               & reverse_raw, error)
            IF (ALLOCATED(error)) RETURN
            IF (ik == 1) THEN
               ALLOCATE (reverse_sum(SIZE(reverse_raw)), &
                  & source=(0.0_wp, 0.0_wp))
            END IF
            reverse_sum = reverse_sum + reverse_raw
         END DO
         owner_sum = SUM(owner_matrix, dim=2)
         ncomplex = 6*SIZE(density, 1)*SIZE(density, 2)*nk*SIZE(density, 3) &
            & + 5*SIZE(density, 1)*SIZE(density, 3)
         IF (cp2k_exchange_partial_reverse_size(reverse_buffer(1)) /= &
            & INT(ncomplex, int64)) THEN
            CALL test_failed(error, "Incorrect partial reverse memory accounting")
            RETURN
         END IF
         CALL cp2k_exchange_partial_reverse_unpack(reverse_sum_buffer, &
            & descriptor(1), reverse_sum, error)
         IF (ALLOCATED(error)) RETURN
         CALL cp2k_exchange_partial_import_reverse(stream(1), descriptor(1), &
            & owner_sum, reverse_sum_buffer, calc, mol, cache, cache_plan, error)
         IF (ALLOCATED(error)) RETURN

         ALLOCATE (overlap_adjoint(SIZE(overlap_adjoint_ref, 1), &
            & SIZE(overlap_adjoint_ref, 2), nk), source=(0.0_wp, 0.0_wp))
         ALLOCATE (gradient_partial(3, mol%nat), sigma_partial(3, 3))
         CALL cp2k_exchange_stream_reverse_batch_apply(stream(1), calc, mol, &
            & cache, cache_plan, error)
         IF (ALLOCATED(error)) RETURN
         DO ik = 1, nk
            CALL cp2k_exchange_stream_reverse_batch_pull(stream(1), ik, &
               & density(:, :, :, ik), overlap(:, :, ik), &
               & overlap_adjoint(:, :, ik), error)
            IF (ALLOCATED(error)) RETURN
         END DO
         CALL cp2k_exchange_stream_reverse_batch_advance(stream(1), done, error)
         IF (ALLOCATED(error)) RETURN
         IF (.NOT. done) THEN
            CALL test_failed(error, "Single reverse partial batch did not finish")
            RETURN
         END IF
         CALL cp2k_exchange_stream_reverse_result(stream(1), calc, mol, cache, cache_plan, &
            & gradient_partial, sigma_partial, error)
         IF (ALLOCATED(error)) RETURN
         maxerr = MAX(maxerr, &
            & MAXVAL(ABS(overlap_adjoint - overlap_adjoint_ref)), &
            & MAXVAL(ABS(gradient_partial - gradient_ref)), &
            & MAXVAL(ABS(sigma_partial - sigma_ref)))
         CALL cp2k_exchange_stream_end(stream(1), error)
         IF (ALLOCATED(error)) RETURN
         DO ik = 2, npart
            CALL cp2k_exchange_partial_discard(stream(ik), error)
            IF (ALLOCATED(error)) RETURN
         END DO
         DEALLOCATE (stream, descriptor, reverse_buffer, owner_matrix, reverse_sum, &
            & overlap_adjoint, gradient_partial, sigma_partial)
      END DO
   END SUBROUTINE check_partial_exchange_case

! **************************************************************************************************
!> \brief ...
!> \param calc ...
!> \param mol ...
!> \param cache ...
!> \param cache_plan ...
!> \param wfn ...
!> \param nmesh ...
!> \param kfrac ...
!> \param weights ...
!> \param density ...
!> \param overlap ...
!> \param error ...
!> \param maxerr ...
! **************************************************************************************************
   SUBROUTINE check_partial_exchange_multibatch(calc, mol, cache, cache_plan, wfn, nmesh, &
      & kfrac, weights, density, overlap, error, maxerr)
      TYPE(xtb_calculator), INTENT(inout)                :: calc
      TYPE(structure_type), INTENT(in)                   :: mol
      TYPE(container_cache), INTENT(inout)               :: cache
      TYPE(cp2k_stream_plan_type), INTENT(inout)         :: cache_plan
      TYPE(wavefunction_type), INTENT(in)                :: wfn
      INTEGER, INTENT(in)                                :: nmesh(3)
      REAL(wp), INTENT(in)                               :: kfrac(:, :), weights(:)
      COMPLEX(wp), INTENT(in)                            :: density(:, :, :, :), overlap(:, :, :)
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error
      REAL(wp), INTENT(out)                              :: maxerr

      INTEGER, PARAMETER                                 :: npart = 4

      COMPLEX(wp), ALLOCATABLE :: fock_partial(:, :, :, :), fock_ref(:, :, :, :), forward_raw(:), &
         forward_sum(:), overlap_adjoint(:, :, :), overlap_adjoint_ref(:, :, :), reverse_raw(:), &
         reverse_sum(:)
      INTEGER                                            :: ibatch, ik, ipart, ncomplex, nk
      INTEGER(int64)                                     :: transaction_id(2)
      INTEGER, ALLOCATABLE                               :: owner_count(:), owner_matrix(:, :), &
                                                            owner_sum(:), pull_owner(:)
      LOGICAL                                            :: contributor_done, done
      REAL(wp)                                           :: energy_partial, energy_ref
      REAL(wp), ALLOCATABLE :: gradient_partial(:, :), gradient_ref(:, :), sigma_partial(:, :), &
         sigma_ref(:, :), vsh_partial(:), vsh_ref(:)
      TYPE(cp2k_exchange_partial_descriptor_type)        :: contributor_descriptor, &
                                                            descriptor(npart), stale_descriptor
      TYPE(cp2k_exchange_partial_forward_type)           :: contributor_forward_buffer, &
                                                            forward_buffer(npart), &
                                                            forward_sum_buffer
      TYPE(cp2k_exchange_partial_reverse_type)           :: contributor_reverse_buffer, &
                                                            reverse_buffer(npart), &
                                                            reverse_sum_buffer
      TYPE(cp2k_exchange_stream_type)                    :: contributor_stream, stream(npart)

      nk = SIZE(weights)
      ALLOCATE (pull_owner(nk))
      DO ik = 1, nk
         pull_owner(ik) = 1 + MODULO(ik - 1, npart)
      END DO
      ALLOCATE (fock_ref(SIZE(density, 1), SIZE(density, 2), &
         & SIZE(density, 3), nk), vsh_ref(SIZE(wfn%qsh, 1)), &
         & overlap_adjoint_ref(SIZE(overlap, 1), SIZE(overlap, 2), nk), &
         & gradient_ref(3, mol%nat), sigma_ref(3, 3))
      CALL cp2k_exchange_kmesh(calc, mol, cache, cache_plan, wfn, nmesh, kfrac, weights, &
         & density, overlap, fock_ref, vsh_ref, energy_ref, error)
      IF (ALLOCATED(error)) RETURN
      CALL cp2k_exchange_kmesh_gradient(calc, mol, cache, cache_plan, nmesh, kfrac, weights, &
         & density, overlap, overlap_adjoint_ref, gradient_ref, sigma_ref, error)
      IF (ALLOCATED(error)) RETURN

      ALLOCATE (owner_matrix(nk, npart))
      ALLOCATE (fock_partial(SIZE(fock_ref, 1), SIZE(fock_ref, 2), &
         & SIZE(fock_ref, 3), nk), source=(0.0_wp, 0.0_wp))
      ALLOCATE (vsh_partial(SIZE(vsh_ref)))
      CALL next_partial_test_transaction(transaction_id)
      DO ipart = 1, npart
         CALL cp2k_exchange_partial_begin(stream(ipart), calc, mol, cache, cache_plan, wfn, &
            & nmesh, kfrac, weights, 1, transaction_id, error)
         IF (ALLOCATED(error)) RETURN
      END DO
      CALL cp2k_exchange_partial_begin(contributor_stream, calc, mol, cache, cache_plan, wfn, &
         & nmesh, kfrac, weights, 1, transaction_id, error)
      IF (ALLOCATED(error)) RETURN
      DO ibatch = 1, nk
         owner_matrix = 0
         DO ik = nk, 1, -1
            ipart = 1 + MODULO(ik - 1, npart)
            CALL cp2k_exchange_partial_push(stream(ipart), ik, &
               & density(:, :, :, ik), overlap(:, :, ik), error)
            IF (ALLOCATED(error)) RETURN
         END DO
         IF (ALLOCATED(forward_sum)) DEALLOCATE (forward_sum)
         DO ipart = 1, npart
            CALL cp2k_exchange_partial_export_forward(stream(ipart), &
               & descriptor(ipart), owner_count, forward_buffer(ipart), error)
            IF (ALLOCATED(error)) RETURN
            CALL check_partial_descriptor_roundtrip(descriptor(ipart), error)
            IF (ALLOCATED(error)) RETURN
            owner_matrix(:, ipart) = owner_count
            IF (ipart > 1) THEN
               IF (.NOT. cp2k_exchange_partial_descriptors_match( &
                  & descriptor(1), descriptor(ipart))) THEN
                  CALL test_failed(error, &
                     & "Multibatch forward descriptors differ")
                  RETURN
               END IF
            END IF
            CALL cp2k_exchange_partial_forward_pack(forward_buffer(ipart), &
               & forward_raw, error)
            IF (ALLOCATED(error)) RETURN
            IF (ipart == 1) THEN
               ALLOCATE (forward_sum(SIZE(forward_raw)), &
                  & source=(0.0_wp, 0.0_wp))
            END IF
            forward_sum = forward_sum + forward_raw
         END DO
         owner_sum = SUM(owner_matrix, dim=2)
         ncomplex = 4*SIZE(density, 1)*SIZE(density, 2)*SIZE(density, 3)
         IF (cp2k_exchange_partial_forward_size(forward_buffer(1)) /= &
            & INT(ncomplex, int64)) THEN
            CALL test_failed(error, "Incorrect B=1 forward payload size")
            RETURN
         END IF
         CALL cp2k_exchange_partial_forward_unpack(forward_sum_buffer, &
            & descriptor(1), forward_sum, error)
         IF (ALLOCATED(error)) RETURN
         CALL cp2k_exchange_partial_export_forward(contributor_stream, &
            & contributor_descriptor, owner_count, contributor_forward_buffer, &
            & error)
         IF (ALLOCATED(error)) RETURN
         IF (ANY(owner_count /= 0) .OR. &
            & .NOT. cp2k_exchange_partial_descriptors_match(descriptor(1), &
            & contributor_descriptor)) THEN
            CALL test_failed(error, &
               & "Empty forward contributor descriptor/ownership differs")
            RETURN
         END IF
         IF (ibatch == 1) stale_descriptor = descriptor(1)
         IF (ibatch == 2) THEN
            CALL cp2k_exchange_partial_import_forward(stream(1), &
               & stale_descriptor, owner_sum, forward_sum_buffer, calc, mol, &
               & cache, cache_plan, error)
            IF (.NOT. ALLOCATED(error)) THEN
               CALL test_failed(error, "A stale forward generation was accepted")
               RETURN
            END IF
            DEALLOCATE (error)
         END IF
         DO ipart = 1, npart
            CALL cp2k_exchange_partial_import_forward(stream(ipart), &
               & descriptor(1), owner_sum, forward_sum_buffer, calc, mol, cache, cache_plan, &
               & error)
            IF (ALLOCATED(error)) RETURN
            IF (ibatch == 1) THEN
               IF (ipart == 1) THEN
                  CALL cp2k_exchange_partial_set_pull_ownership(stream(ipart), &
                     & pull_owner, ipart, npart, 0, error)
                  IF (.NOT. ALLOCATED(error)) THEN
                     CALL test_failed(error, &
                        & "Invalid distributed result owner was accepted")
                     RETURN
                  END IF
                  DEALLOCATE (error)
               END IF
               CALL cp2k_exchange_partial_set_pull_ownership(stream(ipart), &
                  & pull_owner, ipart, npart, 1, error)
               IF (ALLOCATED(error)) RETURN
               IF (ipart == 1) THEN
                  CALL cp2k_exchange_partial_set_pull_ownership(stream(ipart), &
                     & pull_owner, ipart, npart, 1, error)
                  IF (.NOT. ALLOCATED(error)) THEN
                     CALL test_failed(error, &
                        & "Duplicate distributed pull ownership was accepted")
                     RETURN
                  END IF
                  DEALLOCATE (error)
               END IF
            END IF
         END DO
         DO ipart = 1, npart
            CALL cp2k_exchange_stream_batch_apply(stream(ipart), calc, mol, &
               & cache, cache_plan, error)
            IF (ALLOCATED(error)) RETURN
            IF (ibatch == 1 .AND. ipart == 1) THEN
               CALL cp2k_exchange_stream_batch_advance(stream(ipart), &
                  & contributor_done, error)
               IF (.NOT. ALLOCATED(error)) THEN
                  CALL test_failed(error, &
                     & "Distributed forward advance accepted a missing owned pull")
                  RETURN
               END IF
               DEALLOCATE (error)
            END IF
            IF (ibatch == 1 .AND. ipart == npart) THEN
               CALL cp2k_exchange_stream_batch_pull(stream(ipart), 1, &
                  & overlap(:, :, 1), fock_partial(:, :, :, 1), error)
               IF (.NOT. ALLOCATED(error)) THEN
                  CALL test_failed(error, "An unowned forward pull was accepted")
                  RETURN
               END IF
               DEALLOCATE (error)
            END IF
            DO ik = 1, nk
               IF (pull_owner(ik) /= ipart) CYCLE
               CALL cp2k_exchange_stream_batch_pull(stream(ipart), ik, &
                  & overlap(:, :, ik), fock_partial(:, :, :, ik), error)
               IF (ALLOCATED(error)) RETURN
            END DO
            CALL cp2k_exchange_stream_batch_advance(stream(ipart), &
               & contributor_done, error)
            IF (ALLOCATED(error)) RETURN
            IF (ipart == 1) THEN
               done = contributor_done
            ELSE IF (contributor_done .NEQV. done) THEN
               CALL test_failed(error, &
                  & "Distributed forward importer generations diverged")
               RETURN
            END IF
         END DO
         CALL cp2k_exchange_partial_contributor_advance(contributor_stream, &
            & contributor_done, error)
         IF (ALLOCATED(error)) RETURN
         IF (contributor_done .NEQV. done) THEN
            CALL test_failed(error, &
               & "Forward contributor and importer generations diverged")
            RETURN
         END IF
         IF (done .NEQV. (ibatch == nk)) THEN
            CALL test_failed(error, "Unexpected multibatch forward completion")
            RETURN
         END IF
      END DO
      CALL cp2k_exchange_stream_batch_result(stream(2), calc, mol, cache, cache_plan, &
         & vsh_partial, energy_partial, error)
      IF (.NOT. ALLOCATED(error)) THEN
         CALL test_failed(error, &
            & "A non-owner retrieved the distributed forward scalar result")
         RETURN
      END IF
      DEALLOCATE (error)
      CALL cp2k_exchange_stream_batch_result(stream(1), calc, mol, cache, cache_plan, &
         & vsh_partial, energy_partial, error)
      IF (ALLOCATED(error)) RETURN
      maxerr = MAX(MAXVAL(ABS(fock_partial - fock_ref)), &
         & MAXVAL(ABS(vsh_partial - vsh_ref)), ABS(energy_partial - energy_ref))
      DO ipart = 1, npart
         CALL cp2k_exchange_stream_end(stream(ipart), error)
         IF (ALLOCATED(error)) RETURN
      END DO
      CALL cp2k_exchange_partial_discard(contributor_stream, error)
      IF (ALLOCATED(error)) RETURN

      DEALLOCATE (owner_matrix)
      ALLOCATE (owner_matrix(nk, npart))
      ALLOCATE (overlap_adjoint(SIZE(overlap_adjoint_ref, 1), &
         & SIZE(overlap_adjoint_ref, 2), nk), source=(0.0_wp, 0.0_wp))
      ALLOCATE (gradient_partial(3, mol%nat), sigma_partial(3, 3))
      CALL next_partial_test_transaction(transaction_id)
      DO ipart = 1, npart
         CALL cp2k_exchange_partial_begin(stream(ipart), calc, mol, cache, cache_plan, wfn, &
            & nmesh, kfrac, weights, 1, transaction_id, error, &
            & reverse_only=.TRUE.)
         IF (ALLOCATED(error)) RETURN
      END DO
      CALL cp2k_exchange_partial_begin(contributor_stream, calc, mol, cache, cache_plan, wfn, &
         & nmesh, kfrac, weights, 1, transaction_id, error, reverse_only=.TRUE.)
      IF (ALLOCATED(error)) RETURN
      DO ibatch = 1, nk
         owner_matrix = 0
         DO ik = nk, 1, -1
            ipart = 1 + MODULO(ik - 1, npart)
            CALL cp2k_exchange_partial_reverse_push(stream(ipart), ik, &
               & density(:, :, :, ik), overlap(:, :, ik), error)
            IF (ALLOCATED(error)) RETURN
         END DO
         IF (ALLOCATED(reverse_sum)) DEALLOCATE (reverse_sum)
         DO ipart = 1, npart
            CALL cp2k_exchange_partial_export_reverse(stream(ipart), &
               & descriptor(ipart), owner_count, reverse_buffer(ipart), error)
            IF (ALLOCATED(error)) RETURN
            CALL check_partial_descriptor_roundtrip(descriptor(ipart), error)
            IF (ALLOCATED(error)) RETURN
            owner_matrix(:, ipart) = owner_count
            IF (ipart > 1) THEN
               IF (.NOT. cp2k_exchange_partial_descriptors_match( &
                  & descriptor(1), descriptor(ipart))) THEN
                  CALL test_failed(error, &
                     & "Multibatch reverse descriptors differ")
                  RETURN
               END IF
            END IF
            CALL cp2k_exchange_partial_reverse_pack(reverse_buffer(ipart), &
               & reverse_raw, error)
            IF (ALLOCATED(error)) RETURN
            IF (ipart == 1) THEN
               ALLOCATE (reverse_sum(SIZE(reverse_raw)), &
                  & source=(0.0_wp, 0.0_wp))
            END IF
            reverse_sum = reverse_sum + reverse_raw
         END DO
         owner_sum = SUM(owner_matrix, dim=2)
         ncomplex = 6*SIZE(density, 1)*SIZE(density, 2)*SIZE(density, 3) &
            & + 5*SIZE(density, 1)*SIZE(density, 3)
         IF (cp2k_exchange_partial_reverse_size(reverse_buffer(1)) /= &
            & INT(ncomplex, int64)) THEN
            CALL test_failed(error, "Incorrect B=1 reverse payload size")
            RETURN
         END IF
         CALL cp2k_exchange_partial_reverse_unpack(reverse_sum_buffer, &
            & descriptor(1), reverse_sum, error)
         IF (ALLOCATED(error)) RETURN
         CALL cp2k_exchange_partial_export_reverse(contributor_stream, &
            & contributor_descriptor, owner_count, contributor_reverse_buffer, &
            & error)
         IF (ALLOCATED(error)) RETURN
         IF (ANY(owner_count /= 0) .OR. &
            & .NOT. cp2k_exchange_partial_descriptors_match(descriptor(1), &
            & contributor_descriptor)) THEN
            CALL test_failed(error, &
               & "Empty reverse contributor descriptor/ownership differs")
            RETURN
         END IF
         IF (ibatch == 1) stale_descriptor = descriptor(1)
         IF (ibatch == 2) THEN
            CALL cp2k_exchange_partial_import_reverse(stream(1), &
               & stale_descriptor, owner_sum, reverse_sum_buffer, calc, mol, &
               & cache, cache_plan, error)
            IF (.NOT. ALLOCATED(error)) THEN
               CALL test_failed(error, "A stale reverse generation was accepted")
               RETURN
            END IF
            DEALLOCATE (error)
         END IF
         DO ipart = 1, npart
            CALL cp2k_exchange_partial_import_reverse(stream(ipart), &
               & descriptor(1), owner_sum, reverse_sum_buffer, calc, mol, cache, cache_plan, &
               & error)
            IF (ALLOCATED(error)) RETURN
            IF (ibatch == 1) THEN
               CALL cp2k_exchange_partial_set_pull_ownership(stream(ipart), &
                  & pull_owner, ipart, npart, 1, error)
               IF (ALLOCATED(error)) RETURN
            END IF
         END DO
         DO ipart = 1, npart
            CALL cp2k_exchange_stream_reverse_batch_apply(stream(ipart), calc, &
               & mol, cache, cache_plan, error)
            IF (ALLOCATED(error)) RETURN
            IF (ibatch == 1 .AND. ipart == 1) THEN
               CALL cp2k_exchange_stream_reverse_batch_advance(stream(ipart), &
                  & contributor_done, error)
               IF (.NOT. ALLOCATED(error)) THEN
                  CALL test_failed(error, &
                     & "Distributed reverse advance accepted a missing owned pull")
                  RETURN
               END IF
               DEALLOCATE (error)
            END IF
            IF (ibatch == 1 .AND. ipart == npart) THEN
               CALL cp2k_exchange_stream_reverse_batch_pull(stream(ipart), 1, &
                  & density(:, :, :, 1), overlap(:, :, 1), &
                  & overlap_adjoint(:, :, 1), error)
               IF (.NOT. ALLOCATED(error)) THEN
                  CALL test_failed(error, "An unowned reverse pull was accepted")
                  RETURN
               END IF
               DEALLOCATE (error)
            END IF
            DO ik = 1, nk
               IF (pull_owner(ik) /= ipart) CYCLE
               CALL cp2k_exchange_stream_reverse_batch_pull(stream(ipart), ik, &
                  & density(:, :, :, ik), overlap(:, :, ik), &
                  & overlap_adjoint(:, :, ik), error)
               IF (ALLOCATED(error)) RETURN
            END DO
            CALL cp2k_exchange_stream_reverse_batch_advance(stream(ipart), &
               & contributor_done, error)
            IF (ALLOCATED(error)) RETURN
            IF (ipart == 1) THEN
               done = contributor_done
            ELSE IF (contributor_done .NEQV. done) THEN
               CALL test_failed(error, &
                  & "Distributed reverse importer generations diverged")
               RETURN
            END IF
         END DO
         CALL cp2k_exchange_partial_contributor_advance(contributor_stream, &
            & contributor_done, error)
         IF (ALLOCATED(error)) RETURN
         IF (contributor_done .NEQV. done) THEN
            CALL test_failed(error, &
               & "Reverse contributor and importer generations diverged")
            RETURN
         END IF
         IF (done .NEQV. (ibatch == nk)) THEN
            CALL test_failed(error, "Unexpected multibatch reverse completion")
            RETURN
         END IF
      END DO
      CALL cp2k_exchange_stream_reverse_result(stream(2), calc, mol, cache, cache_plan, &
         & gradient_partial, sigma_partial, error)
      IF (.NOT. ALLOCATED(error)) THEN
         CALL test_failed(error, &
            & "A non-owner retrieved the distributed reverse direct result")
         RETURN
      END IF
      DEALLOCATE (error)
      CALL cp2k_exchange_stream_reverse_result(stream(1), calc, mol, cache, cache_plan, &
         & gradient_partial, sigma_partial, error)
      IF (ALLOCATED(error)) RETURN
      maxerr = MAX(maxerr, &
         & MAXVAL(ABS(overlap_adjoint - overlap_adjoint_ref)), &
         & MAXVAL(ABS(gradient_partial - gradient_ref)), &
         & MAXVAL(ABS(sigma_partial - sigma_ref)))
      DO ipart = 1, npart
         CALL cp2k_exchange_stream_end(stream(ipart), error)
         IF (ALLOCATED(error)) RETURN
      END DO
      CALL cp2k_exchange_partial_discard(contributor_stream, error)
      IF (ALLOCATED(error)) RETURN
   END SUBROUTINE check_partial_exchange_multibatch

! **************************************************************************************************
!> \brief ...
!> \param calc ...
!> \param mol ...
!> \param cache ...
!> \param cache_plan ...
!> \param wfn ...
!> \param nmesh ...
!> \param kfrac ...
!> \param weights ...
!> \param density ...
!> \param overlap ...
!> \param error ...
!> \param maxerr ...
! **************************************************************************************************
   SUBROUTINE check_partial_exchange_distributed_images(calc, mol, cache, cache_plan, wfn, &
      & nmesh, kfrac, weights, density, overlap, error, maxerr)
      TYPE(xtb_calculator), INTENT(inout)                :: calc
      TYPE(structure_type), INTENT(in)                   :: mol
      TYPE(container_cache), INTENT(inout)               :: cache
      TYPE(cp2k_stream_plan_type), INTENT(inout)         :: cache_plan
      TYPE(wavefunction_type), INTENT(in)                :: wfn
      INTEGER, INTENT(in)                                :: nmesh(3)
      REAL(wp), INTENT(in)                               :: kfrac(:, :), weights(:)
      COMPLEX(wp), INTENT(in)                            :: density(:, :, :, :), overlap(:, :, :)
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error
      REAL(wp), INTENT(out)                              :: maxerr

      INTEGER, PARAMETER                                 :: batch_size = 2, npart = 4

      COMPLEX(wp), ALLOCATABLE :: density_trial(:, :, :), fock_local(:, :, :), &
         fock_partial(:, :, :, :), fock_range(:, :, :), fock_ref(:, :, :, :), forward_raw(:), &
         forward_sum(:), overlap_adjoint(:, :, :), overlap_adjoint_local(:, :), &
         overlap_adjoint_ref(:, :, :), overlap_trial(:, :), reverse_raw(:), reverse_sum(:)
      INTEGER                                            :: base_count, batch_index, extra_count, &
                                                            first_image, ik, image_count, ipart, nk
      INTEGER(int64)                                     :: transaction_id(2)
      INTEGER, ALLOCATABLE                               :: image_cover(:), owner_count(:), &
                                                            owner_matrix(:, :), owner_sum(:), &
                                                            pull_owner(:)
      LOGICAL                                            :: done, done_part, first_batch
      REAL(wp) :: energy_err, energy_local, energy_partial, energy_ref, fock_err, gradient_err, &
         overlap_adjoint_err, sigma_err, vsh_err
      REAL(wp), ALLOCATABLE :: gradient_local(:, :), gradient_partial(:, :), gradient_ref(:, :), &
         sigma_local(:, :), sigma_partial(:, :), sigma_ref(:, :), vsh_local(:), vsh_partial(:), &
         vsh_range(:), vsh_ref(:)
      TYPE(cp2k_exchange_partial_descriptor_type)        :: descriptor(npart)
      TYPE(cp2k_exchange_partial_forward_type)           :: forward_buffer(npart), forward_sum_buffer
      TYPE(cp2k_exchange_partial_reverse_type)           :: reverse_buffer(npart), reverse_sum_buffer
      TYPE(cp2k_exchange_stream_type)                    :: range_stream, stream(npart)

      nk = SIZE(weights)
      ALLOCATE (fock_ref(SIZE(density, 1), SIZE(density, 2), &
         & SIZE(density, 3), nk), vsh_ref(SIZE(wfn%qsh, 1)), &
         & overlap_adjoint_ref(SIZE(overlap, 1), SIZE(overlap, 2), nk), &
         & gradient_ref(3, mol%nat), sigma_ref(3, 3))
      CALL cp2k_exchange_kmesh(calc, mol, cache, cache_plan, wfn, nmesh, kfrac, weights, &
         & density, overlap, fock_ref, vsh_ref, energy_ref, error)
      IF (ALLOCATED(error)) RETURN
      CALL cp2k_exchange_kmesh_gradient(calc, mol, cache, cache_plan, nmesh, kfrac, weights, &
         & density, overlap, overlap_adjoint_ref, gradient_ref, sigma_ref, error)
      IF (ALLOCATED(error)) RETURN

      ALLOCATE (image_cover(nk), source=0)
      base_count = nk/npart
      extra_count = MODULO(nk, npart)
      DO ipart = 1, npart
         image_count = base_count + MERGE(1, 0, ipart <= extra_count)
         IF (image_count > 0) THEN
            first_image = 1 + (ipart - 1)*base_count + MIN(ipart - 1, extra_count)
            image_cover(first_image:first_image + image_count - 1) = &
               & image_cover(first_image:first_image + image_count - 1) + 1
         END IF
      END DO
      IF (ANY(image_cover /= 1)) THEN
         CALL test_failed(error, "Distributed image ranges do not exactly cover")
         RETURN
      END IF

      ALLOCATE (owner_matrix(nk, npart), owner_sum(nk), pull_owner(nk), source=0)
      ALLOCATE (fock_local(SIZE(fock_ref, 1), SIZE(fock_ref, 2), &
         & SIZE(fock_ref, 3)), source=(0.0_wp, 0.0_wp))
      ALLOCATE (fock_partial(SIZE(fock_ref, 1), SIZE(fock_ref, 2), &
         & SIZE(fock_ref, 3), nk), source=(0.0_wp, 0.0_wp))
      pull_owner = 1

      ! A distributed assignment must remain inside a pre-restricted base image
      ! window.  Exercise the public state machine and then finish the stream so
      ! the negative check cannot hide an invalid or leaking transaction.
      CALL next_partial_test_transaction(transaction_id)
      CALL cp2k_exchange_partial_begin(range_stream, calc, mol, cache, cache_plan, wfn, &
         & nmesh, kfrac, weights, nk, transaction_id, error, image_first=2, &
         & image_count=nk - 1)
      IF (ALLOCATED(error)) RETURN
      DO ik = 1, nk
         CALL cp2k_exchange_partial_push(range_stream, ik, density(:, :, :, ik), &
            & overlap(:, :, ik), error)
         IF (ALLOCATED(error)) RETURN
      END DO
      CALL cp2k_exchange_partial_export_forward(range_stream, descriptor(1), &
         & owner_count, forward_buffer(1), error)
      IF (ALLOCATED(error)) RETURN
      CALL cp2k_exchange_partial_import_forward(range_stream, descriptor(1), &
         & owner_count, forward_buffer(1), calc, mol, cache, cache_plan, error)
      IF (ALLOCATED(error)) RETURN
      CALL cp2k_exchange_partial_set_image_range(range_stream, 1, 1, error)
      IF (.NOT. ALLOCATED(error)) THEN
         CALL test_failed(error, &
            & "Distributed ownership escaped the base image window")
         RETURN
      END IF
      DEALLOCATE (error)
      CALL cp2k_exchange_partial_set_pull_ownership(range_stream, pull_owner, &
         & 1, 1, 1, error)
      IF (ALLOCATED(error)) RETURN
      CALL cp2k_exchange_partial_set_image_range(range_stream, 2, nk - 1, error)
      IF (.NOT. ALLOCATED(error)) THEN
         CALL test_failed(error, &
            & "Pull ownership was combined with distributed images")
         RETURN
      END IF
      DEALLOCATE (error)
      CALL cp2k_exchange_stream_batch_apply(range_stream, calc, mol, cache, cache_plan, error)
      IF (ALLOCATED(error)) RETURN
      ALLOCATE (fock_range(SIZE(fock_ref, 1), SIZE(fock_ref, 2), &
         & SIZE(fock_ref, 3)), source=(0.0_wp, 0.0_wp))
      DO ik = 1, nk
         fock_range = (0.0_wp, 0.0_wp)
         CALL cp2k_exchange_stream_batch_pull(range_stream, ik, &
            & overlap(:, :, ik), fock_range, error)
         IF (ALLOCATED(error)) RETURN
      END DO
      CALL cp2k_exchange_stream_batch_advance(range_stream, done, error)
      IF (ALLOCATED(error)) RETURN
      IF (.NOT. done) THEN
         CALL test_failed(error, "Restricted distributed image stream did not finish")
         RETURN
      END IF
      ALLOCATE (vsh_range(SIZE(vsh_ref)))
      CALL cp2k_exchange_stream_batch_result(range_stream, calc, mol, cache, cache_plan, &
         & vsh_range, energy_local, error)
      IF (ALLOCATED(error)) RETURN
      CALL cp2k_exchange_stream_end(range_stream, error)
      IF (ALLOCATED(error)) RETURN
      DEALLOCATE (fock_range, vsh_range)

      CALL next_partial_test_transaction(transaction_id)
      DO ipart = 1, npart
         CALL cp2k_exchange_partial_begin(stream(ipart), calc, mol, cache, cache_plan, wfn, &
            & nmesh, kfrac, weights, batch_size, transaction_id, error)
         IF (ALLOCATED(error)) RETURN
      END DO

      done = .FALSE.
      first_batch = .TRUE.
      batch_index = 0
      DO WHILE (.NOT. done)
         batch_index = batch_index + 1
         IF (batch_index > nk) THEN
            CALL test_failed(error, "Distributed forward image batches did not finish")
            RETURN
         END IF
         owner_matrix = 0
         DO ik = 1, nk
            ipart = 1 + MODULO(ik - 1, npart)
            CALL cp2k_exchange_partial_push(stream(ipart), ik, &
               & density(:, :, :, ik), overlap(:, :, ik), error)
            IF (ALLOCATED(error)) RETURN
         END DO
         IF (ALLOCATED(forward_sum)) DEALLOCATE (forward_sum)
         DO ipart = 1, npart
            CALL cp2k_exchange_partial_export_forward(stream(ipart), &
               & descriptor(ipart), owner_count, forward_buffer(ipart), error)
            IF (ALLOCATED(error)) RETURN
            owner_matrix(:, ipart) = owner_count
            IF (ipart > 1) THEN
               IF (.NOT. cp2k_exchange_partial_descriptors_match( &
                  & descriptor(1), descriptor(ipart))) THEN
                  CALL test_failed(error, &
                     & "Distributed forward image descriptors differ")
                  RETURN
               END IF
            END IF
            CALL cp2k_exchange_partial_forward_pack(forward_buffer(ipart), &
               & forward_raw, error)
            IF (ALLOCATED(error)) RETURN
            IF (ipart == 1) THEN
               ALLOCATE (forward_sum(SIZE(forward_raw)), &
                  & source=(0.0_wp, 0.0_wp))
            END IF
            forward_sum = forward_sum + forward_raw
         END DO
         owner_sum = SUM(owner_matrix, dim=2)
         CALL cp2k_exchange_partial_forward_unpack(forward_sum_buffer, &
            & descriptor(1), forward_sum, error)
         IF (ALLOCATED(error)) RETURN
         DO ipart = 1, npart
            IF (first_batch .AND. ipart == 1) THEN
               CALL cp2k_exchange_partial_set_image_range(stream(ipart), &
                  & 1, 1, error)
               IF (.NOT. ALLOCATED(error)) THEN
                  CALL test_failed(error, &
                     & "Distributed image range was accepted before import")
                  RETURN
               END IF
               DEALLOCATE (error)
            END IF
            CALL cp2k_exchange_partial_import_forward(stream(ipart), &
               & descriptor(1), owner_sum, forward_sum_buffer, calc, mol, cache, cache_plan, &
               & error)
            IF (ALLOCATED(error)) RETURN
            IF (first_batch) THEN
               image_count = base_count + MERGE(1, 0, ipart <= extra_count)
               IF (image_count > 0) THEN
                  first_image = 1 + (ipart - 1)*base_count &
                     & + MIN(ipart - 1, extra_count)
               ELSE
                  first_image = 0
               END IF
               IF (ipart == 1) THEN
                  CALL cp2k_exchange_partial_set_image_range(stream(ipart), &
                     & 1, -1, error)
                  IF (.NOT. ALLOCATED(error)) THEN
                     CALL test_failed(error, &
                        & "Negative distributed image count was accepted")
                     RETURN
                  END IF
                  DEALLOCATE (error)
                  CALL cp2k_exchange_partial_set_image_range(stream(ipart), &
                     & 0, 1, error)
                  IF (.NOT. ALLOCATED(error)) THEN
                     CALL test_failed(error, &
                        & "Positive distributed image range started at zero")
                     RETURN
                  END IF
                  DEALLOCATE (error)
                  CALL cp2k_exchange_partial_set_image_range(stream(ipart), &
                     & nk, 2, error)
                  IF (.NOT. ALLOCATED(error)) THEN
                     CALL test_failed(error, &
                        & "Distributed image range crossed the mesh end")
                     RETURN
                  END IF
                  DEALLOCATE (error)
                  CALL cp2k_exchange_partial_set_image_range(stream(ipart), &
                     & 1, 0, error)
                  IF (.NOT. ALLOCATED(error)) THEN
                     CALL test_failed(error, &
                        & "Invalid empty distributed image range was accepted")
                     RETURN
                  END IF
                  DEALLOCATE (error)
               END IF
               CALL cp2k_exchange_partial_set_image_range(stream(ipart), &
                  & first_image, image_count, error)
               IF (ALLOCATED(error)) RETURN
               IF (ipart == 1) THEN
                  CALL cp2k_exchange_partial_set_pull_ownership(stream(ipart), &
                     & pull_owner, ipart, npart, 1, error)
                  IF (.NOT. ALLOCATED(error)) THEN
                     CALL test_failed(error, &
                        & "Distributed images accepted disjoint k-point pulls")
                     RETURN
                  END IF
                  DEALLOCATE (error)
                  CALL cp2k_exchange_partial_set_image_range(stream(ipart), &
                     & first_image, image_count, error)
                  IF (.NOT. ALLOCATED(error)) THEN
                     CALL test_failed(error, &
                        & "Duplicate distributed image range was accepted")
                     RETURN
                  END IF
                  DEALLOCATE (error)
               END IF
            END IF
         END DO
         DO ipart = 1, npart
            CALL cp2k_exchange_stream_batch_apply(stream(ipart), calc, mol, &
               & cache, cache_plan, error)
            IF (ALLOCATED(error)) RETURN
            IF (first_batch .AND. ipart == 1) THEN
               CALL cp2k_exchange_partial_set_image_range(stream(ipart), &
                  & 1, 1, error)
               IF (.NOT. ALLOCATED(error)) THEN
                  CALL test_failed(error, &
                     & "Distributed image range followed kernel apply")
                  RETURN
               END IF
               DEALLOCATE (error)
            END IF
            DO ik = 1, nk
               fock_local = (0.0_wp, 0.0_wp)
               CALL cp2k_exchange_stream_batch_pull(stream(ipart), ik, &
                  & overlap(:, :, ik), fock_local, error)
               IF (ALLOCATED(error)) RETURN
               image_count = base_count + MERGE(1, 0, ipart <= extra_count)
               IF (image_count == 0 .AND. &
                  & ANY(fock_local /= (0.0_wp, 0.0_wp))) THEN
                  CALL test_failed(error, &
                     & "Empty image importer produced a forward Fock response")
                  RETURN
               END IF
               fock_partial(:, :, :, ik) = fock_partial(:, :, :, ik) + fock_local
            END DO
            CALL cp2k_exchange_stream_batch_advance(stream(ipart), &
               & done_part, error)
            IF (ALLOCATED(error)) RETURN
            IF (first_batch .AND. ipart == 1 .AND. .NOT. done_part) THEN
               CALL cp2k_exchange_partial_set_image_range(stream(ipart), &
                  & 1, 1, error)
               IF (.NOT. ALLOCATED(error)) THEN
                  CALL test_failed(error, &
                     & "Distributed image range followed an earlier batch")
                  RETURN
               END IF
               DEALLOCATE (error)
            END IF
            IF (ipart == 1) THEN
               done = done_part
            ELSE IF (done_part .NEQV. done) THEN
               CALL test_failed(error, &
                  & "Distributed forward image generations diverged")
               RETURN
            END IF
         END DO
         first_batch = .FALSE.
      END DO

      ALLOCATE (vsh_local(SIZE(vsh_ref)), vsh_partial(SIZE(vsh_ref)), &
         & source=0.0_wp)
      energy_partial = 0.0_wp
      DO ipart = 1, npart
         CALL cp2k_exchange_stream_batch_result(stream(ipart), calc, mol, cache, cache_plan, &
            & vsh_local, energy_local, error)
         IF (ALLOCATED(error)) RETURN
         image_count = base_count + MERGE(1, 0, ipart <= extra_count)
         IF (image_count == 0 .AND. &
            & (ANY(vsh_local /= 0.0_wp) .OR. energy_local /= 0.0_wp)) THEN
            CALL test_failed(error, &
               & "Empty image importer produced energy or shell potential")
            RETURN
         END IF
         vsh_partial = vsh_partial + vsh_local
         energy_partial = energy_partial + energy_local
         CALL cp2k_exchange_stream_reverse_begin(stream(ipart), mol, error)
         IF (ALLOCATED(error)) RETURN
      END DO
      fock_err = MAXVAL(ABS(fock_partial - fock_ref))
      vsh_err = MAXVAL(ABS(vsh_partial - vsh_ref))
      energy_err = ABS(energy_partial - energy_ref)
      maxerr = MAX(fock_err, vsh_err, energy_err)

      ALLOCATE (overlap_adjoint(SIZE(overlap_adjoint_ref, 1), &
         & SIZE(overlap_adjoint_ref, 2), nk), source=(0.0_wp, 0.0_wp))
      ALLOCATE (overlap_adjoint_local(SIZE(overlap_adjoint_ref, 1), &
         & SIZE(overlap_adjoint_ref, 2)), source=(0.0_wp, 0.0_wp))
      ALLOCATE (gradient_local(3, mol%nat), gradient_partial(3, mol%nat), &
         & sigma_local(3, 3), sigma_partial(3, 3), source=0.0_wp)
      done = .FALSE.
      first_batch = .TRUE.
      batch_index = 0
      DO WHILE (.NOT. done)
         batch_index = batch_index + 1
         IF (batch_index > nk) THEN
            CALL test_failed(error, "Distributed reverse image batches did not finish")
            RETURN
         END IF
         owner_matrix = 0
         DO ik = 1, nk
            ipart = 1 + MODULO(ik - 1, npart)
            IF (first_batch .AND. ik == 1) THEN
               ALLOCATE (density_trial, source=density(:, :, :, ik))
               density_trial(1, 1, 1) = density_trial(1, 1, 1) &
                  & + CMPLX(1.0e-7_wp, 0.0_wp, wp)
               CALL cp2k_exchange_partial_reverse_push(stream(ipart), ik, &
                  & density_trial, overlap(:, :, ik), error)
               IF (.NOT. ALLOCATED(error)) THEN
                  CALL test_failed(error, &
                     & "Changed density crossed the forward/reverse boundary")
                  RETURN
               END IF
               DEALLOCATE (error, density_trial)
               ALLOCATE (overlap_trial, source=overlap(:, :, ik))
               overlap_trial(1, 1) = overlap_trial(1, 1) &
                  & + CMPLX(1.0e-7_wp, 0.0_wp, wp)
               CALL cp2k_exchange_partial_reverse_push(stream(ipart), ik, &
                  & density(:, :, :, ik), overlap_trial, error)
               IF (.NOT. ALLOCATED(error)) THEN
                  CALL test_failed(error, &
                     & "Changed overlap crossed the forward/reverse boundary")
                  RETURN
               END IF
               DEALLOCATE (error, overlap_trial)
            END IF
            CALL cp2k_exchange_partial_reverse_push(stream(ipart), ik, &
               & density(:, :, :, ik), overlap(:, :, ik), error)
            IF (ALLOCATED(error)) RETURN
         END DO
         IF (ALLOCATED(reverse_sum)) DEALLOCATE (reverse_sum)
         DO ipart = 1, npart
            CALL cp2k_exchange_partial_export_reverse(stream(ipart), &
               & descriptor(ipart), owner_count, reverse_buffer(ipart), error)
            IF (ALLOCATED(error)) RETURN
            owner_matrix(:, ipart) = owner_count
            IF (ipart > 1) THEN
               IF (.NOT. cp2k_exchange_partial_descriptors_match( &
                  & descriptor(1), descriptor(ipart))) THEN
                  CALL test_failed(error, &
                     & "Distributed reverse image descriptors differ")
                  RETURN
               END IF
            END IF
            CALL cp2k_exchange_partial_reverse_pack(reverse_buffer(ipart), &
               & reverse_raw, error)
            IF (ALLOCATED(error)) RETURN
            IF (ipart == 1) THEN
               ALLOCATE (reverse_sum(SIZE(reverse_raw)), &
                  & source=(0.0_wp, 0.0_wp))
            END IF
            reverse_sum = reverse_sum + reverse_raw
         END DO
         owner_sum = SUM(owner_matrix, dim=2)
         CALL cp2k_exchange_partial_reverse_unpack(reverse_sum_buffer, &
            & descriptor(1), reverse_sum, error)
         IF (ALLOCATED(error)) RETURN
         DO ipart = 1, npart
            CALL cp2k_exchange_partial_import_reverse(stream(ipart), &
               & descriptor(1), owner_sum, reverse_sum_buffer, calc, mol, cache, cache_plan, &
               & error)
            IF (ALLOCATED(error)) RETURN
            IF (first_batch .AND. ipart == 1) THEN
               CALL cp2k_exchange_partial_set_image_range(stream(ipart), &
                  & 1, base_count + MERGE(1, 0, 1 <= extra_count), error)
               IF (.NOT. ALLOCATED(error)) THEN
                  CALL test_failed(error, &
                     & "Reverse phase changed a fixed distributed image range")
                  RETURN
               END IF
               DEALLOCATE (error)
            END IF
         END DO
         DO ipart = 1, npart
            CALL cp2k_exchange_stream_reverse_batch_apply(stream(ipart), calc, &
               & mol, cache, cache_plan, error)
            IF (ALLOCATED(error)) RETURN
            DO ik = 1, nk
               overlap_adjoint_local = (0.0_wp, 0.0_wp)
               IF (batch_index == 2 .AND. ipart == 1 .AND. ik == 1) THEN
                  ALLOCATE (density_trial, source=density(:, :, :, ik))
                  density_trial(1, 1, 1) = density_trial(1, 1, 1) &
                     & + CMPLX(1.0e-7_wp, 0.0_wp, wp)
                  CALL cp2k_exchange_stream_reverse_batch_pull(stream(ipart), &
                     & ik, density_trial, overlap(:, :, ik), &
                     & overlap_adjoint_local, error)
                  IF (.NOT. ALLOCATED(error)) THEN
                     CALL test_failed(error, &
                        & "Changed reverse density crossed an image batch")
                     RETURN
                  END IF
                  DEALLOCATE (error, density_trial)
                  ALLOCATE (overlap_trial, source=overlap(:, :, ik))
                  overlap_trial(1, 1) = overlap_trial(1, 1) &
                     & + CMPLX(1.0e-7_wp, 0.0_wp, wp)
                  CALL cp2k_exchange_stream_reverse_batch_pull(stream(ipart), &
                     & ik, density(:, :, :, ik), overlap_trial, &
                     & overlap_adjoint_local, error)
                  IF (.NOT. ALLOCATED(error)) THEN
                     CALL test_failed(error, &
                        & "Changed reverse overlap crossed an image batch")
                     RETURN
                  END IF
                  DEALLOCATE (error, overlap_trial)
               END IF
               CALL cp2k_exchange_stream_reverse_batch_pull(stream(ipart), ik, &
                  & density(:, :, :, ik), overlap(:, :, ik), &
                  & overlap_adjoint_local, error)
               IF (ALLOCATED(error)) RETURN
               image_count = base_count + MERGE(1, 0, ipart <= extra_count)
               IF (image_count == 0 .AND. &
                  & ANY(overlap_adjoint_local /= (0.0_wp, 0.0_wp))) THEN
                  CALL test_failed(error, &
                     & "Empty image importer produced an overlap adjoint")
                  RETURN
               END IF
               overlap_adjoint(:, :, ik) = overlap_adjoint(:, :, ik) &
                  & + overlap_adjoint_local
            END DO
            CALL cp2k_exchange_stream_reverse_batch_advance(stream(ipart), &
               & done_part, error)
            IF (ALLOCATED(error)) RETURN
            IF (ipart == 1) THEN
               done = done_part
            ELSE IF (done_part .NEQV. done) THEN
               CALL test_failed(error, &
                  & "Distributed reverse image generations diverged")
               RETURN
            END IF
         END DO
         first_batch = .FALSE.
      END DO

      DO ipart = 1, npart
         CALL cp2k_exchange_stream_reverse_result(stream(ipart), calc, mol, &
            & cache, cache_plan, gradient_local, sigma_local, error)
         IF (ALLOCATED(error)) RETURN
         image_count = base_count + MERGE(1, 0, ipart <= extra_count)
         IF (image_count == 0 .AND. &
            & (ANY(gradient_local /= 0.0_wp) .OR. &
            & ANY(sigma_local /= 0.0_wp))) THEN
            CALL test_failed(error, &
               & "Empty image importer produced force or stress response")
            RETURN
         END IF
         gradient_partial = gradient_partial + gradient_local
         sigma_partial = sigma_partial + sigma_local
         CALL cp2k_exchange_stream_end(stream(ipart), error)
         IF (ALLOCATED(error)) RETURN
      END DO
      overlap_adjoint_err = MAXVAL(ABS(overlap_adjoint - overlap_adjoint_ref))
      gradient_err = MAXVAL(ABS(gradient_partial - gradient_ref))
      sigma_err = MAXVAL(ABS(sigma_partial - sigma_ref))
      maxerr = MAX(maxerr, overlap_adjoint_err, gradient_err, sigma_err)
      PRINT *, "DISTRIBUTED_IMAGE_KERNEL_FORWARD_COMPONENTS", &
         & fock_err, vsh_err, energy_err
      PRINT *, "DISTRIBUTED_IMAGE_KERNEL_REVERSE_COMPONENTS", &
         & overlap_adjoint_err, gradient_err, sigma_err

      ! Repeat the derivative contraction with fresh reverse-only streams.  This
      ! exercises the reverse setter branch independently of a preceding forward
      ! phase and permits a different host-side image partition when required.
      overlap_adjoint = (0.0_wp, 0.0_wp)
      gradient_partial = 0.0_wp
      sigma_partial = 0.0_wp
      CALL next_partial_test_transaction(transaction_id)
      DO ipart = 1, npart
         CALL cp2k_exchange_partial_begin(stream(ipart), calc, mol, cache, cache_plan, wfn, &
            & nmesh, kfrac, weights, batch_size, transaction_id, error, &
            & reverse_only=.TRUE.)
         IF (ALLOCATED(error)) RETURN
      END DO
      done = .FALSE.
      first_batch = .TRUE.
      batch_index = 0
      DO WHILE (.NOT. done)
         batch_index = batch_index + 1
         IF (batch_index > nk) THEN
            CALL test_failed(error, &
               & "Reverse-only distributed image batches did not finish")
            RETURN
         END IF
         owner_matrix = 0
         DO ik = 1, nk
            ipart = 1 + MODULO(ik - 1, npart)
            CALL cp2k_exchange_partial_reverse_push(stream(ipart), ik, &
               & density(:, :, :, ik), overlap(:, :, ik), error)
            IF (ALLOCATED(error)) RETURN
         END DO
         IF (ALLOCATED(reverse_sum)) DEALLOCATE (reverse_sum)
         DO ipart = 1, npart
            CALL cp2k_exchange_partial_export_reverse(stream(ipart), &
               & descriptor(ipart), owner_count, reverse_buffer(ipart), error)
            IF (ALLOCATED(error)) RETURN
            owner_matrix(:, ipart) = owner_count
            IF (ipart > 1 .AND. &
               & .NOT. cp2k_exchange_partial_descriptors_match( &
               & descriptor(1), descriptor(ipart))) THEN
               CALL test_failed(error, &
                  & "Reverse-only distributed image descriptors differ")
               RETURN
            END IF
            CALL cp2k_exchange_partial_reverse_pack(reverse_buffer(ipart), &
               & reverse_raw, error)
            IF (ALLOCATED(error)) RETURN
            IF (ipart == 1) THEN
               ALLOCATE (reverse_sum(SIZE(reverse_raw)), &
                  & source=(0.0_wp, 0.0_wp))
            END IF
            reverse_sum = reverse_sum + reverse_raw
         END DO
         owner_sum = SUM(owner_matrix, dim=2)
         CALL cp2k_exchange_partial_reverse_unpack(reverse_sum_buffer, &
            & descriptor(1), reverse_sum, error)
         IF (ALLOCATED(error)) RETURN
         DO ipart = 1, npart
            CALL cp2k_exchange_partial_import_reverse(stream(ipart), &
               & descriptor(1), owner_sum, reverse_sum_buffer, calc, mol, cache, cache_plan, &
               & error)
            IF (ALLOCATED(error)) RETURN
            IF (first_batch) THEN
               image_count = base_count + MERGE(1, 0, ipart <= extra_count)
               IF (image_count > 0) THEN
                  first_image = 1 + (ipart - 1)*base_count &
                     & + MIN(ipart - 1, extra_count)
               ELSE
                  first_image = 0
               END IF
               CALL cp2k_exchange_partial_set_image_range(stream(ipart), &
                  & first_image, image_count, error)
               IF (ALLOCATED(error)) RETURN
            END IF
         END DO
         DO ipart = 1, npart
            CALL cp2k_exchange_stream_reverse_batch_apply(stream(ipart), calc, &
               & mol, cache, cache_plan, error)
            IF (ALLOCATED(error)) RETURN
            DO ik = 1, nk
               overlap_adjoint_local = (0.0_wp, 0.0_wp)
               CALL cp2k_exchange_stream_reverse_batch_pull(stream(ipart), ik, &
                  & density(:, :, :, ik), overlap(:, :, ik), &
                  & overlap_adjoint_local, error)
               IF (ALLOCATED(error)) RETURN
               image_count = base_count + MERGE(1, 0, ipart <= extra_count)
               IF (image_count == 0 .AND. &
                  & ANY(overlap_adjoint_local /= (0.0_wp, 0.0_wp))) THEN
                  CALL test_failed(error, &
                     & "Empty reverse-only importer produced an adjoint")
                  RETURN
               END IF
               overlap_adjoint(:, :, ik) = overlap_adjoint(:, :, ik) &
                  & + overlap_adjoint_local
            END DO
            CALL cp2k_exchange_stream_reverse_batch_advance(stream(ipart), &
               & done_part, error)
            IF (ALLOCATED(error)) RETURN
            IF (ipart == 1) THEN
               done = done_part
            ELSE IF (done_part .NEQV. done) THEN
               CALL test_failed(error, &
                  & "Reverse-only distributed image generations diverged")
               RETURN
            END IF
         END DO
         first_batch = .FALSE.
      END DO
      DO ipart = 1, npart
         CALL cp2k_exchange_stream_reverse_result(stream(ipart), calc, mol, &
            & cache, cache_plan, gradient_local, sigma_local, error)
         IF (ALLOCATED(error)) RETURN
         image_count = base_count + MERGE(1, 0, ipart <= extra_count)
         IF (image_count == 0 .AND. &
            & (ANY(gradient_local /= 0.0_wp) .OR. &
            & ANY(sigma_local /= 0.0_wp))) THEN
            CALL test_failed(error, &
               & "Empty reverse-only importer produced force or stress")
            RETURN
         END IF
         gradient_partial = gradient_partial + gradient_local
         sigma_partial = sigma_partial + sigma_local
         CALL cp2k_exchange_stream_end(stream(ipart), error)
         IF (ALLOCATED(error)) RETURN
      END DO
      overlap_adjoint_err = MAXVAL(ABS(overlap_adjoint - overlap_adjoint_ref))
      gradient_err = MAXVAL(ABS(gradient_partial - gradient_ref))
      sigma_err = MAXVAL(ABS(sigma_partial - sigma_ref))
      maxerr = MAX(maxerr, overlap_adjoint_err, gradient_err, sigma_err)
      PRINT *, "DISTRIBUTED_IMAGE_KERNEL_REVERSE_ONLY_COMPONENTS", &
         & overlap_adjoint_err, gradient_err, sigma_err
   END SUBROUTINE check_partial_exchange_distributed_images

! **************************************************************************************************
!> \brief ...
!> \param calc ...
!> \param mol ...
!> \param cache ...
!> \param cache_plan ...
!> \param wfn ...
!> \param nmesh ...
!> \param kfrac ...
!> \param weights ...
!> \param density ...
!> \param overlap ...
!> \param error ...
!> \param maxerr ...
!> \param transform_backend ...
! **************************************************************************************************
   SUBROUTINE check_partial_exchange_range_local(calc, mol, cache, cache_plan, wfn, nmesh, &
      & kfrac, weights, density, overlap, error, maxerr, transform_backend)
      TYPE(xtb_calculator), INTENT(inout)                :: calc
      TYPE(structure_type), INTENT(in)                   :: mol
      TYPE(container_cache), INTENT(inout)               :: cache
      TYPE(cp2k_stream_plan_type), INTENT(inout)         :: cache_plan
      TYPE(wavefunction_type), INTENT(in)                :: wfn
      INTEGER, INTENT(in)                                :: nmesh(3)
      REAL(wp), INTENT(in)                               :: kfrac(:, :), weights(:)
      COMPLEX(wp), INTENT(in)                            :: density(:, :, :, :), overlap(:, :, :)
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error
      REAL(wp), INTENT(out)                              :: maxerr
      INTEGER, INTENT(in), OPTIONAL                      :: transform_backend

      INTEGER, PARAMETER                                 :: npart = 4

      COMPLEX(wp), ALLOCATABLE :: fock_local(:, :, :), fock_partial(:, :, :, :), &
         fock_ref(:, :, :, :), forward_full_raw(:), forward_range_raw(:), forward_raw(:), &
         overlap_adjoint(:, :, :), overlap_adjoint_local(:, :), overlap_adjoint_ref(:, :, :), &
         reverse_full_raw(:), reverse_range_raw(:), reverse_raw(:)
      INTEGER :: base_count, batch_count, batch_first, batch_index, batch_size, extra_count, &
         first_image, ik, image_count, ipart, jpart, last_image, local_count, nk, range_first, &
         range_last, selected_transform_backend, transform_pair_batch_size
      INTEGER(int64) :: compact_peak, forward_range_scalars, forward_replicated_scalars, &
         replicated_peak, reverse_range_scalars, reverse_replicated_scalars, transaction_id(2), &
         wrong_transaction(2)
      INTEGER, ALLOCATABLE                               :: owner_bad(:), owner_count(:), &
                                                            owner_matrix(:, :), owner_sum(:), &
                                                            range_count(:), range_start(:)
      LOGICAL                                            :: contributor_done, done, done_part, &
                                                            first_batch, transform_done, &
                                                            transform_done_part
      REAL(wp)                                           :: energy_local, energy_partial, energy_ref
      REAL(wp), ALLOCATABLE :: gradient_local(:, :), gradient_partial(:, :), gradient_ref(:, :), &
         sigma_local(:, :), sigma_partial(:, :), sigma_ref(:, :), vsh_local(:), vsh_partial(:), &
         vsh_ref(:)
      TYPE(cp2k_exchange_partial_descriptor_type)        :: descriptor(npart)
      TYPE(cp2k_exchange_partial_forward_type)           :: forward_buffer(npart), &
                                                            forward_full_buffer, &
                                                            forward_range_buffer(npart)
      TYPE(cp2k_exchange_partial_reverse_type)           :: reverse_buffer(npart), &
                                                            reverse_full_buffer, &
                                                            reverse_range_buffer(npart)
      TYPE(cp2k_exchange_stream_type)                    :: contributor(npart), importer(npart), &
                                                            wrong_stream

      nk = SIZE(weights)
      batch_size = MIN(nk, MAX(2, (nk + 1)/2))
      selected_transform_backend = cp2k_exchange_transform_dense
      IF (PRESENT(transform_backend)) THEN
         selected_transform_backend = transform_backend
      END IF
      transform_pair_batch_size = MIN(batch_size, MAX(1, SIZE(density, 1)**2 - 1))
      ALLOCATE (fock_ref(SIZE(density, 1), SIZE(density, 2), &
         & SIZE(density, 3), nk), vsh_ref(SIZE(wfn%qsh, 1)), &
         & overlap_adjoint_ref(SIZE(overlap, 1), SIZE(overlap, 2), nk), &
         & gradient_ref(3, mol%nat), sigma_ref(3, 3))
      CALL cp2k_exchange_kmesh(calc, mol, cache, cache_plan, wfn, nmesh, kfrac, weights, &
         & density, overlap, fock_ref, vsh_ref, energy_ref, error)
      IF (ALLOCATED(error)) RETURN
      CALL cp2k_exchange_kmesh_gradient(calc, mol, cache, cache_plan, nmesh, kfrac, &
         & weights, density, overlap, overlap_adjoint_ref, gradient_ref, &
         & sigma_ref, error)
      IF (ALLOCATED(error)) RETURN

      ALLOCATE (range_start(npart), range_count(npart))
      base_count = nk/npart
      extra_count = MODULO(nk, npart)
      DO ipart = 1, npart
         range_count(ipart) = base_count + MERGE(1, 0, ipart <= extra_count)
         IF (range_count(ipart) > 0) THEN
            range_start(ipart) = 1 + (ipart - 1)*base_count &
               & + MIN(ipart - 1, extra_count)
         ELSE
            range_start(ipart) = 0
         END IF
      END DO
      IF (SUM(range_count) /= nk) THEN
         CALL test_failed(error, "Range-local image partition is incomplete")
         RETURN
      END IF

      ALLOCATE (owner_matrix(nk, npart), owner_sum(nk), owner_bad(nk))
      ALLOCATE (fock_local(SIZE(fock_ref, 1), SIZE(fock_ref, 2), &
         & SIZE(fock_ref, 3)), source=(0.0_wp, 0.0_wp))
      ALLOCATE (fock_partial(SIZE(fock_ref, 1), SIZE(fock_ref, 2), &
         & SIZE(fock_ref, 3), nk), source=(0.0_wp, 0.0_wp))
      ALLOCATE (vsh_local(SIZE(vsh_ref)), vsh_partial(SIZE(vsh_ref)), &
         & source=0.0_wp)

      CALL next_partial_test_transaction(transaction_id)
      DO ipart = 1, npart
         CALL cp2k_exchange_partial_begin(contributor(ipart), calc, mol, cache, cache_plan, &
            & wfn, nmesh, kfrac, weights, batch_size, transaction_id, error, &
            & transform_backend=selected_transform_backend, &
            & transform_pair_batch_size=transform_pair_batch_size)
         IF (ALLOCATED(error)) RETURN
         CALL cp2k_exchange_partial_begin(importer(ipart), calc, mol, cache, cache_plan, &
            & wfn, nmesh, kfrac, weights, batch_size, transaction_id, error, &
            & local_image_first=range_start(ipart), &
            & local_image_count=range_count(ipart), &
            & transform_backend=selected_transform_backend, &
            & transform_pair_batch_size=transform_pair_batch_size)
         IF (ALLOCATED(error)) RETURN
      END DO
      replicated_peak = cp2k_exchange_stream_peak_complex_elements( &
         & contributor(1))
      compact_peak = 0_int64
      DO ipart = 1, npart
         compact_peak = compact_peak &
            & + cp2k_exchange_stream_peak_complex_elements(importer(ipart))
      END DO
      IF (compact_peak >= INT(npart, int64)*replicated_peak) THEN
         CALL test_failed(error, &
            & "Range-local importers did not reduce persistent complex storage")
         RETURN
      END IF

      done = .FALSE.
      first_batch = .TRUE.
      batch_index = 0
      forward_range_scalars = 0_int64
      forward_replicated_scalars = 0_int64
      DO WHILE (.NOT. done)
         batch_index = batch_index + 1
         batch_first = 1 + (batch_index - 1)*batch_size
         batch_count = MIN(batch_size, nk - batch_first + 1)
         owner_matrix = 0
         transform_done = .FALSE.
         DO WHILE (.NOT. transform_done)
            DO ik = 1, nk
               ipart = 1 + MODULO(ik - 1, npart)
               CALL cp2k_exchange_partial_push(contributor(ipart), ik, &
                  & density(:, :, :, ik), overlap(:, :, ik), error)
               IF (ALLOCATED(error)) RETURN
            END DO
            IF (selected_transform_backend == cp2k_exchange_transform_dense) THEN
               transform_done = .TRUE.
            ELSE
               DO ipart = 1, npart
                  CALL cp2k_exchange_stream_transform_advance( &
                     & contributor(ipart), transform_done_part, error)
                  IF (ALLOCATED(error)) RETURN
                  IF (ipart == 1) THEN
                     transform_done = transform_done_part
                  ELSE IF (transform_done_part .NEQV. transform_done) THEN
                     CALL test_failed(error, &
                        & "Compact forward contributor transforms diverged")
                     RETURN
                  END IF
               END DO
            END IF
         END DO
         IF (ALLOCATED(forward_full_raw)) DEALLOCATE (forward_full_raw)
         DO ipart = 1, npart
            CALL cp2k_exchange_partial_export_forward(contributor(ipart), &
               & descriptor(ipart), owner_count, forward_buffer(ipart), error)
            IF (ALLOCATED(error)) RETURN
            owner_matrix(:, ipart) = owner_count
            IF (ipart > 1 .AND. &
               & .NOT. cp2k_exchange_partial_descriptors_match(descriptor(1), &
               & descriptor(ipart))) THEN
               CALL test_failed(error, &
                  & "Range-local forward contributor descriptors differ")
               RETURN
            END IF
            CALL cp2k_exchange_partial_forward_pack(forward_buffer(ipart), &
               & forward_raw, error)
            IF (ALLOCATED(error)) RETURN
            IF (ipart == 1) THEN
               ALLOCATE (forward_full_raw(SIZE(forward_raw)), &
                  & source=(0.0_wp, 0.0_wp))
            END IF
            forward_full_raw = forward_full_raw + forward_raw
         END DO
         owner_sum = SUM(owner_matrix, dim=2)
         IF (ANY(owner_sum /= 1)) THEN
            CALL test_failed(error, &
               & "Range-local forward contributors do not exactly cover k")
            RETURN
         END IF
         CALL cp2k_exchange_partial_forward_unpack(forward_full_buffer, &
            & descriptor(1), forward_full_raw, error)
         IF (ALLOCATED(error)) RETURN
         forward_replicated_scalars = forward_replicated_scalars &
            & + INT(npart, int64)*SIZE(forward_full_raw, kind=int64)

         DO ipart = 1, npart
            image_count = range_count(ipart)
            first_image = range_start(ipart)
            IF (image_count > 0) THEN
               last_image = first_image + image_count - 1
               range_first = MAX(batch_first, first_image)
               range_last = MIN(batch_first + batch_count - 1, last_image)
               local_count = MAX(0, range_last - range_first + 1)
               IF (local_count == 0) range_first = 0
            ELSE
               range_first = 0
               local_count = 0
            END IF
            IF (ALLOCATED(forward_range_raw)) DEALLOCATE (forward_range_raw)
            DO jpart = 1, npart
               CALL cp2k_exchange_partial_forward_pack_range( &
                  & forward_buffer(jpart), descriptor(1), range_first, &
                  & local_count, forward_raw, error)
               IF (ALLOCATED(error)) RETURN
               IF (jpart == 1) THEN
                  ALLOCATE (forward_range_raw(SIZE(forward_raw)), &
                     & source=(0.0_wp, 0.0_wp))
               END IF
               forward_range_raw = forward_range_raw + forward_raw
            END DO
            forward_range_scalars = forward_range_scalars &
               & + SIZE(forward_range_raw, kind=int64)
            CALL cp2k_exchange_partial_forward_unpack_range( &
               & forward_range_buffer(ipart), descriptor(1), range_first, &
               & local_count, forward_range_raw, error)
            IF (ALLOCATED(error)) RETURN
         END DO

         IF (first_batch) THEN
            CALL cp2k_exchange_partial_push(importer(1), 1, &
               & density(:, :, :, 1), overlap(:, :, 1), error)
            IF (.NOT. ALLOCATED(error)) THEN
               CALL test_failed(error, &
                  & "Range-local importer accepted a forward contribution")
               RETURN
            END IF
            DEALLOCATE (error)
            CALL cp2k_exchange_partial_import_forward(importer(1), &
               & descriptor(1), owner_sum, forward_full_buffer, calc, mol, cache, cache_plan, &
               & error)
            IF (.NOT. ALLOCATED(error)) THEN
               CALL test_failed(error, &
                  & "Range-local importer accepted a replicated payload")
               RETURN
            END IF
            DEALLOCATE (error)
            owner_bad = owner_sum
            owner_bad(1) = 0
            CALL cp2k_exchange_partial_import_forward_range(importer(1), &
               & descriptor(1), owner_bad, forward_range_buffer(1), calc, mol, &
               & cache, cache_plan, error)
            IF (.NOT. ALLOCATED(error)) THEN
               CALL test_failed(error, &
                  & "Range-local importer accepted incomplete k ownership")
               RETURN
            END IF
            DEALLOCATE (error)
            CALL cp2k_exchange_partial_forward_pack_range(forward_buffer(1), &
               & descriptor(1), 0, 1, forward_raw, error)
            IF (.NOT. ALLOCATED(error)) THEN
               CALL test_failed(error, &
                  & "Forward range pack accepted an invalid image interval")
               RETURN
            END IF
            DEALLOCATE (error)
            CALL next_partial_test_transaction(wrong_transaction)
            CALL cp2k_exchange_partial_begin(wrong_stream, calc, mol, cache, cache_plan, &
               & wfn, nmesh, kfrac, weights, batch_size, wrong_transaction, &
               & error, local_image_first=range_start(1), &
               & local_image_count=range_count(1), &
               & transform_backend=selected_transform_backend, &
               & transform_pair_batch_size=transform_pair_batch_size)
            IF (ALLOCATED(error)) RETURN
            CALL cp2k_exchange_partial_import_forward_range(wrong_stream, &
               & descriptor(1), owner_sum, forward_range_buffer(1), calc, mol, &
               & cache, cache_plan, error)
            IF (.NOT. ALLOCATED(error)) THEN
               CALL test_failed(error, &
                  & "Forward range import accepted a foreign transaction")
               RETURN
            END IF
            DEALLOCATE (error)
            CALL cp2k_exchange_stream_end(wrong_stream, error)
            IF (.NOT. ALLOCATED(error)) THEN
               CALL test_failed(error, &
                  & "Incomplete foreign range transaction closed successfully")
               RETURN
            END IF
            DEALLOCATE (error)
         END IF

         DO ipart = 1, npart
            CALL cp2k_exchange_partial_import_forward_range(importer(ipart), &
               & descriptor(1), owner_sum, forward_range_buffer(ipart), calc, &
               & mol, cache, cache_plan, error)
            IF (ALLOCATED(error)) RETURN
            CALL cp2k_exchange_stream_batch_apply(importer(ipart), calc, mol, &
               & cache, cache_plan, error)
            IF (ALLOCATED(error)) RETURN
            DO ik = 1, nk
               fock_local = (0.0_wp, 0.0_wp)
               CALL cp2k_exchange_stream_batch_pull(importer(ipart), ik, &
                  & overlap(:, :, ik), fock_local, error)
               IF (ALLOCATED(error)) RETURN
               IF (range_count(ipart) == 0 .AND. &
                  & ANY(fock_local /= (0.0_wp, 0.0_wp))) THEN
                  CALL test_failed(error, &
                     & "Empty range-local importer produced a Fock response")
                  RETURN
               END IF
               fock_partial(:, :, :, ik) = fock_partial(:, :, :, ik) + fock_local
            END DO
            CALL cp2k_exchange_stream_batch_advance(importer(ipart), done_part, &
               & error)
            IF (ALLOCATED(error)) RETURN
            IF (ipart == 1) THEN
               done = done_part
            ELSE IF (done_part .NEQV. done) THEN
               CALL test_failed(error, &
                  & "Range-local forward importer generations diverged")
               RETURN
            END IF
         END DO
         DO ipart = 1, npart
            CALL cp2k_exchange_partial_contributor_advance(contributor(ipart), &
               & contributor_done, error)
            IF (ALLOCATED(error)) RETURN
            IF (contributor_done .NEQV. done) THEN
               CALL test_failed(error, &
                  & "Range-local forward contributors diverged from importers")
               RETURN
            END IF
         END DO
         first_batch = .FALSE.
      END DO
      IF (batch_index < 2) THEN
         CALL test_failed(error, "Range-local forward gate did not span batches")
         RETURN
      END IF
      IF (forward_range_scalars >= forward_replicated_scalars) THEN
         CALL test_failed(error, &
            & "Forward range payload did not reduce replicated communication")
         RETURN
      END IF

      energy_partial = 0.0_wp
      DO ipart = 1, npart
         CALL cp2k_exchange_stream_batch_result(importer(ipart), calc, mol, &
            & cache, cache_plan, vsh_local, energy_local, error)
         IF (ALLOCATED(error)) RETURN
         IF (range_count(ipart) == 0 .AND. &
            & (ANY(vsh_local /= 0.0_wp) .OR. energy_local /= 0.0_wp)) THEN
            CALL test_failed(error, &
               & "Empty range-local importer produced a scalar response")
            RETURN
         END IF
         vsh_partial = vsh_partial + vsh_local
         energy_partial = energy_partial + energy_local
         CALL cp2k_exchange_stream_end(importer(ipart), error)
         IF (ALLOCATED(error)) RETURN
         CALL cp2k_exchange_partial_discard(contributor(ipart), error)
         IF (ALLOCATED(error)) RETURN
      END DO
      maxerr = MAX(MAXVAL(ABS(fock_partial - fock_ref)), &
         & MAXVAL(ABS(vsh_partial - vsh_ref)), ABS(energy_partial - energy_ref))

      ALLOCATE (overlap_adjoint(SIZE(overlap_adjoint_ref, 1), &
         & SIZE(overlap_adjoint_ref, 2), nk), source=(0.0_wp, 0.0_wp))
      ALLOCATE (overlap_adjoint_local(SIZE(overlap_adjoint_ref, 1), &
         & SIZE(overlap_adjoint_ref, 2)), source=(0.0_wp, 0.0_wp))
      ALLOCATE (gradient_local(3, mol%nat), gradient_partial(3, mol%nat), &
         & sigma_local(3, 3), sigma_partial(3, 3), source=0.0_wp)

      CALL next_partial_test_transaction(transaction_id)
      DO ipart = 1, npart
         CALL cp2k_exchange_partial_begin(contributor(ipart), calc, mol, cache, cache_plan, &
            & wfn, nmesh, kfrac, weights, batch_size, transaction_id, error, &
            & reverse_only=.TRUE., &
            & transform_backend=selected_transform_backend, &
            & transform_pair_batch_size=transform_pair_batch_size)
         IF (ALLOCATED(error)) RETURN
         CALL cp2k_exchange_partial_begin(importer(ipart), calc, mol, cache, cache_plan, &
            & wfn, nmesh, kfrac, weights, batch_size, transaction_id, error, &
            & reverse_only=.TRUE., local_image_first=range_start(ipart), &
            & local_image_count=range_count(ipart), &
            & transform_backend=selected_transform_backend, &
            & transform_pair_batch_size=transform_pair_batch_size)
         IF (ALLOCATED(error)) RETURN
      END DO
      done = .FALSE.
      first_batch = .TRUE.
      batch_index = 0
      reverse_range_scalars = 0_int64
      reverse_replicated_scalars = 0_int64
      DO WHILE (.NOT. done)
         batch_index = batch_index + 1
         DO ipart = 1, npart
            CALL cp2k_exchange_stream_batch_bounds(importer(ipart), &
               & first_image, last_image, error)
            IF (ALLOCATED(error)) RETURN
            IF (ipart == 1) THEN
               batch_first = first_image
               batch_count = last_image - first_image + 1
            ELSE IF (first_image /= batch_first .OR. &
               & last_image /= batch_first + batch_count - 1) THEN
               CALL test_failed(error, &
                  & "Range-local reverse importer batch bounds diverged")
               RETURN
            END IF
         END DO
         owner_matrix = 0
         transform_done = .FALSE.
         DO WHILE (.NOT. transform_done)
            DO ik = 1, nk
               ipart = 1 + MODULO(ik - 1, npart)
               CALL cp2k_exchange_partial_reverse_push(contributor(ipart), ik, &
                  & density(:, :, :, ik), overlap(:, :, ik), error)
               IF (ALLOCATED(error)) RETURN
            END DO
            IF (selected_transform_backend == cp2k_exchange_transform_dense) THEN
               transform_done = .TRUE.
            ELSE
               DO ipart = 1, npart
                  CALL cp2k_exchange_stream_reverse_transform_advance( &
                     & contributor(ipart), transform_done_part, error)
                  IF (ALLOCATED(error)) RETURN
                  IF (ipart == 1) THEN
                     transform_done = transform_done_part
                  ELSE IF (transform_done_part .NEQV. transform_done) THEN
                     CALL test_failed(error, &
                        & "Compact reverse contributor transforms diverged")
                     RETURN
                  END IF
               END DO
            END IF
         END DO
         IF (ALLOCATED(reverse_full_raw)) DEALLOCATE (reverse_full_raw)
         DO ipart = 1, npart
            CALL cp2k_exchange_partial_export_reverse(contributor(ipart), &
               & descriptor(ipart), owner_count, reverse_buffer(ipart), error)
            IF (ALLOCATED(error)) RETURN
            owner_matrix(:, ipart) = owner_count
            IF (ipart > 1 .AND. &
               & .NOT. cp2k_exchange_partial_descriptors_match(descriptor(1), &
               & descriptor(ipart))) THEN
               CALL test_failed(error, &
                  & "Range-local reverse contributor descriptors differ")
               RETURN
            END IF
            CALL cp2k_exchange_partial_reverse_pack(reverse_buffer(ipart), &
               & reverse_raw, error)
            IF (ALLOCATED(error)) RETURN
            IF (ipart == 1) THEN
               ALLOCATE (reverse_full_raw(SIZE(reverse_raw)), &
                  & source=(0.0_wp, 0.0_wp))
            END IF
            reverse_full_raw = reverse_full_raw + reverse_raw
         END DO
         owner_sum = SUM(owner_matrix, dim=2)
         IF (ANY(owner_sum /= 1)) THEN
            CALL test_failed(error, &
               & "Range-local reverse contributors do not exactly cover k")
            RETURN
         END IF
         CALL cp2k_exchange_partial_reverse_unpack(reverse_full_buffer, &
            & descriptor(1), reverse_full_raw, error)
         IF (ALLOCATED(error)) RETURN
         reverse_replicated_scalars = reverse_replicated_scalars &
            & + INT(npart, int64)*SIZE(reverse_full_raw, kind=int64)

         DO ipart = 1, npart
            image_count = range_count(ipart)
            first_image = range_start(ipart)
            IF (image_count > 0) THEN
               last_image = first_image + image_count - 1
               range_first = MAX(batch_first, first_image)
               range_last = MIN(batch_first + batch_count - 1, last_image)
               local_count = MAX(0, range_last - range_first + 1)
               IF (local_count == 0) range_first = 0
            ELSE
               range_first = 0
               local_count = 0
            END IF
            IF (ALLOCATED(reverse_range_raw)) DEALLOCATE (reverse_range_raw)
            DO jpart = 1, npart
               CALL cp2k_exchange_partial_reverse_pack_range( &
                  & reverse_buffer(jpart), descriptor(1), range_first, &
                  & local_count, reverse_raw, error)
               IF (ALLOCATED(error)) RETURN
               IF (jpart == 1) THEN
                  ALLOCATE (reverse_range_raw(SIZE(reverse_raw)), &
                     & source=(0.0_wp, 0.0_wp))
               END IF
               reverse_range_raw = reverse_range_raw + reverse_raw
            END DO
            reverse_range_scalars = reverse_range_scalars &
               & + SIZE(reverse_range_raw, kind=int64)
            CALL cp2k_exchange_partial_reverse_unpack_range( &
               & reverse_range_buffer(ipart), descriptor(1), range_first, &
               & local_count, reverse_range_raw, error)
            IF (ALLOCATED(error)) RETURN
         END DO

         IF (first_batch) THEN
            CALL cp2k_exchange_partial_import_reverse(importer(1), &
               & descriptor(1), owner_sum, reverse_full_buffer, calc, mol, cache, cache_plan, &
               & error)
            IF (.NOT. ALLOCATED(error)) THEN
               CALL test_failed(error, &
                  & "Range-local importer accepted a replicated reverse payload")
               RETURN
            END IF
            DEALLOCATE (error)
            CALL cp2k_exchange_partial_reverse_pack_range(reverse_buffer(1), &
               & descriptor(1), 0, 1, reverse_raw, error)
            IF (.NOT. ALLOCATED(error)) THEN
               CALL test_failed(error, &
                  & "Reverse range pack accepted an invalid image interval")
               RETURN
            END IF
            DEALLOCATE (error)
         END IF

         DO ipart = 1, npart
            CALL cp2k_exchange_partial_import_reverse_range(importer(ipart), &
               & descriptor(1), owner_sum, reverse_range_buffer(ipart), calc, &
               & mol, cache, cache_plan, error)
            IF (ALLOCATED(error)) RETURN
            CALL cp2k_exchange_stream_reverse_batch_apply(importer(ipart), calc, &
               & mol, cache, cache_plan, error)
            IF (ALLOCATED(error)) RETURN
            DO ik = 1, nk
               overlap_adjoint_local = (0.0_wp, 0.0_wp)
               CALL cp2k_exchange_stream_reverse_batch_pull(importer(ipart), &
                  & ik, density(:, :, :, ik), overlap(:, :, ik), &
                  & overlap_adjoint_local, error)
               IF (ALLOCATED(error)) RETURN
               IF (range_count(ipart) == 0 .AND. &
                  & ANY(overlap_adjoint_local /= (0.0_wp, 0.0_wp))) THEN
                  CALL test_failed(error, &
                     & "Empty range-local importer produced an adjoint")
                  RETURN
               END IF
               overlap_adjoint(:, :, ik) = overlap_adjoint(:, :, ik) &
                  & + overlap_adjoint_local
            END DO
            CALL cp2k_exchange_stream_reverse_batch_advance(importer(ipart), &
               & done_part, error)
            IF (ALLOCATED(error)) RETURN
            IF (ipart == 1) THEN
               done = done_part
            ELSE IF (done_part .NEQV. done) THEN
               CALL test_failed(error, &
                  & "Range-local reverse importer generations diverged")
               RETURN
            END IF
         END DO
         DO ipart = 1, npart
            CALL cp2k_exchange_partial_contributor_advance(contributor(ipart), &
               & contributor_done, error)
            IF (ALLOCATED(error)) RETURN
            IF (contributor_done .NEQV. done) THEN
               CALL test_failed(error, &
                  & "Range-local reverse contributors diverged from importers")
               RETURN
            END IF
         END DO
         first_batch = .FALSE.
      END DO
      CALL cp2k_exchange_stream_batch_bounds(importer(1), first_image, &
         & last_image, error)
      IF (ALLOCATED(error)) THEN
         IF (TRIM(error%message) /= &
            & "g-xTB exchange image batches are complete") THEN
            DEALLOCATE (error)
            CALL test_failed(error, &
               & "Unexpected completed reverse batch-bounds diagnostic")
            RETURN
         END IF
         DEALLOCATE (error)
      ELSE
         CALL test_failed(error, &
            & "Completed reverse image batches returned bounds")
         RETURN
      END IF
      IF (batch_index < 2) THEN
         CALL test_failed(error, "Range-local reverse gate did not span batches")
         RETURN
      END IF
      IF (reverse_range_scalars >= reverse_replicated_scalars) THEN
         CALL test_failed(error, &
            & "Reverse range payload did not reduce replicated communication")
         RETURN
      END IF

      DO ipart = 1, npart
         CALL cp2k_exchange_stream_reverse_result(importer(ipart), calc, mol, &
            & cache, cache_plan, gradient_local, sigma_local, error)
         IF (ALLOCATED(error)) RETURN
         IF (range_count(ipart) == 0 .AND. &
            & (ANY(gradient_local /= 0.0_wp) .OR. &
            & ANY(sigma_local /= 0.0_wp))) THEN
            CALL test_failed(error, &
               & "Empty range-local importer produced force or stress")
            RETURN
         END IF
         gradient_partial = gradient_partial + gradient_local
         sigma_partial = sigma_partial + sigma_local
         CALL cp2k_exchange_stream_end(importer(ipart), error)
         IF (ALLOCATED(error)) RETURN
         CALL cp2k_exchange_partial_discard(contributor(ipart), error)
         IF (ALLOCATED(error)) RETURN
      END DO
      maxerr = MAX(maxerr, &
         & MAXVAL(ABS(overlap_adjoint - overlap_adjoint_ref)), &
         & MAXVAL(ABS(gradient_partial - gradient_ref)), &
         & MAXVAL(ABS(sigma_partial - sigma_ref)))
      PRINT *, "RANGE_LOCAL_IMAGE_IMPORT_MAX_RESID", maxerr
      PRINT *, "RANGE_LOCAL_IMAGE_IMPORT_PEAK_COMPLEX", replicated_peak, &
         & compact_peak
      PRINT *, "RANGE_LOCAL_IMAGE_IMPORT_PAYLOAD_COMPLEX", &
         & forward_replicated_scalars, forward_range_scalars, &
         & reverse_replicated_scalars, reverse_range_scalars
   END SUBROUTINE check_partial_exchange_range_local

! **************************************************************************************************
!> \brief ...
!> \param calc ...
!> \param mol ...
!> \param cache ...
!> \param cache_plan ...
!> \param wfn ...
!> \param nmesh ...
!> \param kfrac ...
!> \param weights ...
!> \param density ...
!> \param overlap ...
!> \param error ...
!> \param separable_max ...
!> \param fft_max ...
!> \param fft_forward ...
!> \param fft_reverse ...
! **************************************************************************************************
   SUBROUTINE check_exchange_transform_selectors(calc, mol, cache, cache_plan, wfn, nmesh, &
      & kfrac, weights, density, overlap, error, separable_max, fft_max, &
      & fft_forward, fft_reverse)
      TYPE(xtb_calculator), INTENT(inout)                :: calc
      TYPE(structure_type), INTENT(in)                   :: mol
      TYPE(container_cache), INTENT(inout)               :: cache
      TYPE(cp2k_stream_plan_type), INTENT(inout)         :: cache_plan
      TYPE(wavefunction_type), INTENT(in)                :: wfn
      INTEGER, INTENT(in)                                :: nmesh(3)
      REAL(wp), INTENT(in)                               :: kfrac(:, :), weights(:)
      COMPLEX(wp), INTENT(in)                            :: density(:, :, :, :), overlap(:, :, :)
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error
      REAL(wp), INTENT(out)                              :: separable_max, fft_max, fft_forward(3), &
                                                            fft_reverse(3)

      COMPLEX(wp), ALLOCATABLE                           :: fock_ref(:, :, :, :), &
                                                            fock_test(:, :, :, :), &
                                                            overlap_grad_ref(:, :, :), &
                                                            overlap_grad_test(:, :, :)
      REAL(wp)                                           :: energy_ref, energy_test
      REAL(wp), ALLOCATABLE                              :: gradient_ref(:, :), gradient_test(:, :), &
                                                            sigma_ref(:, :), sigma_test(:, :), &
                                                            vsh_ref(:), vsh_test(:)

      ALLOCATE (fock_ref(SIZE(density, 1), SIZE(density, 2), SIZE(density, 3), &
         & SIZE(density, 4)), fock_test(SIZE(density, 1), SIZE(density, 2), &
         & SIZE(density, 3), SIZE(density, 4)), vsh_ref(SIZE(wfn%qsh, 1)), &
         & vsh_test(SIZE(wfn%qsh, 1)), &
         & overlap_grad_ref(SIZE(overlap, 1), SIZE(overlap, 2), &
         & SIZE(overlap, 3)), overlap_grad_test(SIZE(overlap, 1), &
         & SIZE(overlap, 2), SIZE(overlap, 3)), &
         & gradient_ref(3, mol%nat), gradient_test(3, mol%nat), &
         & sigma_ref(3, 3), sigma_test(3, 3))

      CALL cp2k_exchange_kmesh(calc, mol, cache, cache_plan, wfn, nmesh, kfrac, weights, &
         & density, overlap, fock_ref, vsh_ref, energy_ref, error, &
         & transform_backend=cp2k_exchange_transform_dense)
      IF (ALLOCATED(error)) RETURN
      CALL cp2k_exchange_kmesh_gradient(calc, mol, cache, cache_plan, nmesh, kfrac, &
         & weights, density, overlap, overlap_grad_ref, gradient_ref, sigma_ref, &
         & error, transform_backend=cp2k_exchange_transform_dense)
      IF (ALLOCATED(error)) RETURN

      CALL cp2k_exchange_kmesh(calc, mol, cache, cache_plan, wfn, nmesh, kfrac, weights, &
         & density, overlap, fock_test, vsh_test, energy_test, error, &
         & transform_backend=cp2k_exchange_transform_separable_dft)
      IF (ALLOCATED(error)) RETURN
      separable_max = MAX(MAXVAL(ABS(fock_test - fock_ref)), &
         & MAXVAL(ABS(vsh_test - vsh_ref)), ABS(energy_test - energy_ref))
      CALL cp2k_exchange_kmesh_gradient(calc, mol, cache, cache_plan, nmesh, kfrac, &
         & weights, density, overlap, overlap_grad_test, gradient_test, &
         & sigma_test, error, &
         & transform_backend=cp2k_exchange_transform_separable_dft)
      IF (ALLOCATED(error)) RETURN
      separable_max = MAX(separable_max, &
         & MAXVAL(ABS(overlap_grad_test - overlap_grad_ref)), &
         & MAXVAL(ABS(gradient_test - gradient_ref)), &
         & MAXVAL(ABS(sigma_test - sigma_ref)))

      CALL cp2k_exchange_kmesh(calc, mol, cache, cache_plan, wfn, nmesh, kfrac, weights, &
         & density, overlap, fock_test, vsh_test, energy_test, error, &
         & transform_backend=cp2k_exchange_transform_mixed_radix_fft)
      IF (ALLOCATED(error)) RETURN
      fft_forward = [MAXVAL(ABS(fock_test - fock_ref)), &
         & MAXVAL(ABS(vsh_test - vsh_ref)), ABS(energy_test - energy_ref)]
      CALL cp2k_exchange_kmesh_gradient(calc, mol, cache, cache_plan, nmesh, kfrac, &
         & weights, density, overlap, overlap_grad_test, gradient_test, &
         & sigma_test, error, &
         & transform_backend=cp2k_exchange_transform_mixed_radix_fft)
      IF (ALLOCATED(error)) RETURN
      fft_reverse = [MAXVAL(ABS(overlap_grad_test - overlap_grad_ref)), &
         & MAXVAL(ABS(gradient_test - gradient_ref)), &
         & MAXVAL(ABS(sigma_test - sigma_ref))]
      fft_max = MAX(MAXVAL(fft_forward), MAXVAL(fft_reverse))
   END SUBROUTINE check_exchange_transform_selectors

! **************************************************************************************************
!> \brief ...
!> \param calc ...
!> \param mol ...
!> \param cache ...
!> \param cache_plan ...
!> \param wfn ...
!> \param nmesh ...
!> \param kfrac ...
!> \param weights ...
!> \param density ...
!> \param overlap ...
!> \param fock_ref ...
!> \param vsh_ref ...
!> \param energy_ref ...
!> \param overlap_grad_ref ...
!> \param gradient_ref ...
!> \param sigma_ref ...
!> \param transform_backend ...
!> \param residual ...
!> \param peak_elements ...
!> \param error ...
!> \param adaptive_tiles ...
!> \param adaptive_images ...
! **************************************************************************************************
   SUBROUTINE check_compact_stream_composition(calc, mol, cache, cache_plan, wfn, nmesh, &
      & kfrac, weights, density, overlap, fock_ref, vsh_ref, energy_ref, &
      & overlap_grad_ref, gradient_ref, sigma_ref, transform_backend, &
      & residual, peak_elements, error, adaptive_tiles, adaptive_images)
      TYPE(xtb_calculator), INTENT(inout)                :: calc
      TYPE(structure_type), INTENT(in)                   :: mol
      TYPE(container_cache), INTENT(inout)               :: cache
      TYPE(cp2k_stream_plan_type), INTENT(inout)         :: cache_plan
      TYPE(wavefunction_type), INTENT(in)                :: wfn
      INTEGER, INTENT(in)                                :: nmesh(3)
      REAL(wp), INTENT(in)                               :: kfrac(:, :), weights(:)
      COMPLEX(wp), INTENT(in)                            :: density(:, :, :, :), overlap(:, :, :), &
                                                            fock_ref(:, :, :, :)
      REAL(wp), INTENT(in)                               :: vsh_ref(:), energy_ref
      COMPLEX(wp), INTENT(in)                            :: overlap_grad_ref(:, :, :)
      REAL(wp), INTENT(in)                               :: gradient_ref(:, :), sigma_ref(:, :)
      INTEGER, INTENT(in)                                :: transform_backend
      REAL(wp), INTENT(out)                              :: residual
      INTEGER(int64), INTENT(out)                        :: peak_elements
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error
      LOGICAL, INTENT(in), OPTIONAL                      :: adaptive_tiles, adaptive_images

      COMPLEX(wp), ALLOCATABLE                           :: fock(:, :, :, :), overlap_grad(:, :, :)
      INTEGER                                            :: ik, image_batch, itile, nao, nk, npair, &
                                                            pair_batch, tile_count
      INTEGER(int64)                                     :: reverse_peak
      LOGICAL                                            :: batch_done, tile_done
      REAL(wp)                                           :: energy
      REAL(wp), ALLOCATABLE                              :: gradient(:, :), sigma(:, :), vsh(:)
      TYPE(cp2k_exchange_stream_type)                    :: stream

      nao = SIZE(density, 1)
      nk = SIZE(density, 4)
      npair = nao*nao
      image_batch = 1
      IF (PRESENT(adaptive_images)) THEN
         IF (adaptive_images) image_batch = cp2k_exchange_stream_image_batch_size(nao, nk, wfn%nspin)
      END IF
      pair_batch = MAX(1, MIN(2, npair - 1))
      IF (PRESENT(adaptive_tiles)) THEN
         IF (adaptive_tiles) pair_batch = cp2k_exchange_stream_pair_tile_size(nao, nk, wfn%nspin, 3)
      END IF
      tile_count = (npair + pair_batch - 1)/pair_batch
      ALLOCATE (fock(SIZE(density, 1), SIZE(density, 2), SIZE(density, 3), &
         & SIZE(density, 4)), source=(0.0_wp, 0.0_wp))
      ALLOCATE (vsh(SIZE(vsh_ref)), source=0.0_wp)

      CALL cp2k_exchange_stream_begin(stream, calc, mol, cache, cache_plan, wfn, nmesh, &
         & kfrac, weights, error, mode=cp2k_exchange_stream_batched, &
         & batch_size=image_batch, transform_backend=transform_backend, &
         & transform_pair_batch_size=pair_batch)
      IF (ALLOCATED(error)) RETURN
      IF (npair > 1 .AND. cp2k_exchange_stream_has_full_mesh_storage(stream)) THEN
         CALL test_failed(error, &
            & "Compact forward stream retained a complete k/AO-pair mesh")
         RETURN
      END IF
      peak_elements = cp2k_exchange_stream_peak_complex_elements(stream)
      batch_done = .FALSE.
      DO WHILE (.NOT. batch_done)
         DO itile = 1, tile_count
            DO ik = 1, nk
               CALL cp2k_exchange_stream_push(stream, ik, &
                  & density(:, :, :, ik), overlap(:, :, ik), error)
               IF (ALLOCATED(error)) RETURN
            END DO
            CALL cp2k_exchange_stream_transform_advance(stream, tile_done, error)
            IF (ALLOCATED(error)) RETURN
            IF (tile_done .NEQV. (itile == tile_count)) THEN
               CALL test_failed(error, &
                  & "Compact forward stream reported an inconsistent tile count")
               RETURN
            END IF
         END DO
         CALL cp2k_exchange_stream_batch_apply(stream, calc, mol, cache, cache_plan, error)
         IF (ALLOCATED(error)) RETURN
         DO ik = 1, nk
            CALL cp2k_exchange_stream_batch_pull(stream, ik, overlap(:, :, ik), &
               & fock(:, :, :, ik), error)
            IF (ALLOCATED(error)) RETURN
         END DO
         CALL cp2k_exchange_stream_batch_advance(stream, batch_done, error)
         IF (ALLOCATED(error)) RETURN
      END DO
      CALL cp2k_exchange_stream_batch_result(stream, calc, mol, cache, cache_plan, vsh, &
         & energy, error)
      IF (ALLOCATED(error)) RETURN
      residual = MAX(MAXVAL(ABS(fock - fock_ref)), MAXVAL(ABS(vsh - vsh_ref)), &
         & ABS(energy - energy_ref))
      CALL cp2k_exchange_stream_end(stream, error)
      IF (ALLOCATED(error)) RETURN

      ALLOCATE (overlap_grad(SIZE(overlap, 1), SIZE(overlap, 2), nk), &
         & source=(0.0_wp, 0.0_wp))
      ALLOCATE (gradient(3, mol%nat), sigma(3, 3), source=0.0_wp)
      IF (PRESENT(adaptive_tiles)) THEN
         IF (adaptive_tiles) pair_batch = cp2k_exchange_stream_pair_tile_size(nao, nk, wfn%nspin, 6)
      END IF
      tile_count = (npair + pair_batch - 1)/pair_batch
      CALL cp2k_exchange_stream_begin(stream, calc, mol, cache, cache_plan, wfn, nmesh, &
         & kfrac, weights, error, mode=cp2k_exchange_stream_batched, &
         & batch_size=image_batch, reverse_only=.TRUE., &
         & transform_backend=transform_backend, &
         & transform_pair_batch_size=pair_batch)
      IF (ALLOCATED(error)) RETURN
      IF (npair > 1 .AND. cp2k_exchange_stream_has_full_mesh_storage(stream)) THEN
         CALL test_failed(error, &
            & "Compact reverse stream retained a complete k/AO-pair mesh")
         RETURN
      END IF
      reverse_peak = cp2k_exchange_stream_peak_complex_elements(stream)
      peak_elements = MAX(peak_elements, reverse_peak)
      batch_done = .FALSE.
      DO WHILE (.NOT. batch_done)
         DO itile = 1, tile_count
            DO ik = 1, nk
               CALL cp2k_exchange_stream_reverse_push(stream, ik, &
                  & density(:, :, :, ik), overlap(:, :, ik), error)
               IF (ALLOCATED(error)) RETURN
            END DO
            CALL cp2k_exchange_stream_reverse_transform_advance(stream, &
               & tile_done, error)
            IF (ALLOCATED(error)) RETURN
            IF (tile_done .NEQV. (itile == tile_count)) THEN
               CALL test_failed(error, &
                  & "Compact reverse stream reported an inconsistent tile count")
               RETURN
            END IF
         END DO
         CALL cp2k_exchange_stream_reverse_batch_apply(stream, calc, mol, cache, cache_plan, &
            & error)
         IF (ALLOCATED(error)) RETURN
         DO ik = 1, nk
            CALL cp2k_exchange_stream_reverse_batch_pull(stream, ik, &
               & density(:, :, :, ik), overlap(:, :, ik), &
               & overlap_grad(:, :, ik), error)
            IF (ALLOCATED(error)) RETURN
         END DO
         CALL cp2k_exchange_stream_reverse_batch_advance(stream, batch_done, &
            & error)
         IF (ALLOCATED(error)) RETURN
      END DO
      CALL cp2k_exchange_stream_reverse_result(stream, calc, mol, cache, cache_plan, &
         & gradient, sigma, error)
      IF (ALLOCATED(error)) RETURN
      residual = MAX(residual, MAXVAL(ABS(overlap_grad - overlap_grad_ref)), &
         & MAXVAL(ABS(gradient - gradient_ref)), MAXVAL(ABS(sigma - sigma_ref)))
      CALL cp2k_exchange_stream_end(stream, error)
   END SUBROUTINE check_compact_stream_composition

! **************************************************************************************************
!> \brief ...
!> \param exchange ...
!> \param cache ...
!> \param field ...
!> \param error ...
! **************************************************************************************************
   SUBROUTINE require_bvk_model_mismatch(exchange, cache, field, error)
      TYPE(exchange_fock), INTENT(in)                    :: exchange
      TYPE(cp2k_stream_plan_type), INTENT(in)            :: cache
      CHARACTER(len=*), INTENT(in)                       :: field
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      IF (cache%model_matches(exchange)) THEN
         CALL test_failed(error, "BvK model signature ignored "//TRIM(field))
      END IF
   END SUBROUTINE require_bvk_model_mismatch

! **************************************************************************************************
!> \brief ...
!> \param exchange ...
!> \param cache ...
!> \param error ...
! **************************************************************************************************
   SUBROUTINE check_bvk_model_signature(exchange, cache, error)
      TYPE(exchange_fock), INTENT(inout)                 :: exchange
      TYPE(cp2k_stream_plan_type), INTENT(in)            :: cache
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      INTEGER, PARAMETER                                 :: integer_delta = 1
      REAL(wp), PARAMETER                                :: real_delta = 1.0e-6_wp

      IF (.NOT. cache%model_matches(exchange)) THEN
         CALL test_failed(error, "Fresh BvK model signature does not match")
         RETURN
      END IF

      exchange%nao = exchange%nao + integer_delta
      CALL require_bvk_model_mismatch(exchange, cache, "nao", error)
      exchange%nao = exchange%nao - integer_delta
      IF (ALLOCATED(error)) RETURN
      exchange%nsh = exchange%nsh + integer_delta
      CALL require_bvk_model_mismatch(exchange, cache, "nsh", error)
      exchange%nsh = exchange%nsh - integer_delta
      IF (ALLOCATED(error)) RETURN
      exchange%maxsh = exchange%maxsh + integer_delta
      CALL require_bvk_model_mismatch(exchange, cache, "maxsh", error)
      exchange%maxsh = exchange%maxsh - integer_delta
      IF (ALLOCATED(error)) RETURN

      exchange%nsh_id(1) = exchange%nsh_id(1) + integer_delta
      CALL require_bvk_model_mismatch(exchange, cache, "nsh_id", error)
      exchange%nsh_id(1) = exchange%nsh_id(1) - integer_delta
      IF (ALLOCATED(error)) RETURN
      exchange%nao_sh(1) = exchange%nao_sh(1) + integer_delta
      CALL require_bvk_model_mismatch(exchange, cache, "nao_sh", error)
      exchange%nao_sh(1) = exchange%nao_sh(1) - integer_delta
      IF (ALLOCATED(error)) RETURN
      exchange%ish_at(1) = exchange%ish_at(1) + integer_delta
      CALL require_bvk_model_mismatch(exchange, cache, "ish_at", error)
      exchange%ish_at(1) = exchange%ish_at(1) - integer_delta
      IF (ALLOCATED(error)) RETURN
      exchange%iao_sh(1) = exchange%iao_sh(1) + integer_delta
      CALL require_bvk_model_mismatch(exchange, cache, "iao_sh", error)
      exchange%iao_sh(1) = exchange%iao_sh(1) - integer_delta
      IF (ALLOCATED(error)) RETURN

      exchange%frscale = exchange%frscale + real_delta
      CALL require_bvk_model_mismatch(exchange, cache, "frscale", error)
      exchange%frscale = cache%bvk_model_frscale
      IF (ALLOCATED(error)) RETURN
      exchange%onsite_frscale = exchange%onsite_frscale + real_delta
      CALL require_bvk_model_mismatch(exchange, cache, "onsite_frscale", error)
      exchange%onsite_frscale = cache%bvk_model_onsite_frscale
      IF (ALLOCATED(error)) RETURN
      exchange%omega = exchange%omega + real_delta
      CALL require_bvk_model_mismatch(exchange, cache, "omega", error)
      exchange%omega = cache%bvk_model_omega
      IF (ALLOCATED(error)) RETURN
      exchange%lrscale = exchange%lrscale + real_delta
      CALL require_bvk_model_mismatch(exchange, cache, "lrscale", error)
      exchange%lrscale = cache%bvk_model_lrscale
      IF (ALLOCATED(error)) RETURN
      exchange%ondiag_scale = exchange%ondiag_scale + real_delta
      CALL require_bvk_model_mismatch(exchange, cache, "ondiag_scale", error)
      exchange%ondiag_scale = cache%bvk_model_ondiag_scale
      IF (ALLOCATED(error)) RETURN
      exchange%hubbard_exp = exchange%hubbard_exp + real_delta
      CALL require_bvk_model_mismatch(exchange, cache, "hubbard_exp", error)
      exchange%hubbard_exp = cache%bvk_model_hubbard_exp
      IF (ALLOCATED(error)) RETURN
      exchange%hubbard_exp_r0 = exchange%hubbard_exp_r0 + real_delta
      CALL require_bvk_model_mismatch(exchange, cache, "hubbard_exp_r0", error)
      exchange%hubbard_exp_r0 = cache%bvk_model_hubbard_exp_r0
      IF (ALLOCATED(error)) RETURN
      exchange%gexp = exchange%gexp + real_delta
      CALL require_bvk_model_mismatch(exchange, cache, "gexp", error)
      exchange%gexp = cache%bvk_model_gexp
      IF (ALLOCATED(error)) RETURN
      exchange%corr_exp = exchange%corr_exp + real_delta
      CALL require_bvk_model_mismatch(exchange, cache, "corr_exp", error)
      exchange%corr_exp = cache%bvk_model_corr_exp
      IF (ALLOCATED(error)) RETURN

      exchange%hubbard(1, 1, 1, 1) = exchange%hubbard(1, 1, 1, 1) + real_delta
      CALL require_bvk_model_mismatch(exchange, cache, "hubbard", error)
      exchange%hubbard(1, 1, 1, 1) = cache%bvk_model_hubbard(1, 1, 1, 1)
      IF (ALLOCATED(error)) RETURN
      exchange%onecxints(1, 1, 1) = exchange%onecxints(1, 1, 1) + real_delta
      CALL require_bvk_model_mismatch(exchange, cache, "onecxints", error)
      exchange%onecxints(1, 1, 1) = cache%bvk_model_onecxints(1, 1, 1)
      IF (ALLOCATED(error)) RETURN
      exchange%offdiag_scale(1, 1, 1, 1) = &
         & exchange%offdiag_scale(1, 1, 1, 1) + real_delta
      CALL require_bvk_model_mismatch(exchange, cache, "offdiag_scale", error)
      exchange%offdiag_scale(1, 1, 1, 1) = &
         & cache%bvk_model_offdiag_scale(1, 1, 1, 1)
      IF (ALLOCATED(error)) RETURN
      exchange%rad(1, 1) = exchange%rad(1, 1) + real_delta
      CALL require_bvk_model_mismatch(exchange, cache, "rad", error)
      exchange%rad(1, 1) = cache%bvk_model_rad(1, 1)
      IF (ALLOCATED(error)) RETURN
      exchange%kq(1, 1) = exchange%kq(1, 1) + real_delta
      CALL require_bvk_model_mismatch(exchange, cache, "kq", error)
      exchange%kq(1, 1) = cache%bvk_model_kq(1, 1)
      IF (ALLOCATED(error)) RETURN
      exchange%corr_scale(1, 1) = exchange%corr_scale(1, 1) + real_delta
      CALL require_bvk_model_mismatch(exchange, cache, "corr_scale", error)
      exchange%corr_scale(1, 1) = cache%bvk_model_corr_scale(1, 1)
      IF (ALLOCATED(error)) RETURN
      exchange%corr_rad(1, 1) = exchange%corr_rad(1, 1) + real_delta
      CALL require_bvk_model_mismatch(exchange, cache, "corr_rad", error)
      exchange%corr_rad(1, 1) = cache%bvk_model_corr_rad(1, 1)
      IF (ALLOCATED(error)) RETURN

      IF (.NOT. cache%model_matches(exchange)) THEN
         CALL test_failed(error, "Restored BvK model signature does not match")
      END IF
   END SUBROUTINE check_bvk_model_signature

! **************************************************************************************************
!> \brief ...
!> \param error ...
! **************************************************************************************************
   SUBROUTINE test_bvk_exchange_kernel(error)
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      INTEGER, PARAMETER :: nk_overflow = 65536, nmesh(3) = [3, 1, 1], nmesh_cube(3) = [2, 2, 2], &
         nmesh_large(3) = [9, 9, 1], nmesh_overflow(3) = [65536, 65537, 1], num(2) = [1, 1], &
         num_sc(6) = [1, 1, 1, 1, 1, 1]
      LOGICAL, PARAMETER :: periodic(3) = [.TRUE., .TRUE., .TRUE.]
      REAL(wp), PARAMETER :: lattice(3, 3) = RESHAPE([3.20_wp, 0.00_wp, 0.00_wp, 0.15_wp, 20.0_wp, &
         0.00_wp, 0.05_wp, 0.10_wp, 21.0_wp], [3, 3])/autoaa, &
         xyz(3, 2) = RESHAPE([0.40_wp, 0.30_wp, 0.20_wp, 1.35_wp, 0.35_wp, 0.25_wp], [3, 2])/autoaa

      TYPE(structure_type)                               :: mol, mol_probe, mol_sc

      class(qvszp_basis_type), ALLOCATABLE :: qbas, qbas_sc
      class(basis_type), ALLOCATABLE :: bas, bas_sc
      class(exchange_type), ALLOCATABLE :: exchange, exchange_sc
      TYPE(basis_cache) :: bcache, bcache_sc
      TYPE(container_cache) :: cache, cache_cp2k, cache_sc
      ! Keep the complete-mesh oracle independent of the production cached plan.
      TYPE(cp2k_stream_plan_type) :: cache_cp2k_plan, oracle_plan
      TYPE(exchange_cache), POINTER :: ptr, ptr_cp2k, ptr_sc
      TYPE(exchange_bvk_kernel) :: kernel, kernel_large, kernel_one, &
         & kernel_order_probe
      TYPE(cp2k_exchange_stream_type) :: stream, stream_large, stream_order, &
         & stream_permuted, stream_reverse, stream_reverse_missing, &
         & stream_reverse_reduced, stream_reverse_batch, &
         & stream_forward_missing, stream_batch, stream_batch_large
      TYPE(cp2k_exchange_stream_type) :: stream_overflow
      TYPE(xtb_calculator) :: calc_cp2k
      TYPE(wavefunction_type) :: wfn, wfn_aux, wfn_large, wfn_sc, wfn_aux_sc
      INTEGER, PARAMETER :: korder(3) = [3, 1, 2]
      INTEGER :: acell, bcell, bvk_builds, i, iat, ibatch_case, icell, ik, ir, &
         & ish, j, jat, jsh
      INTEGER :: first, first_image, ii, izp, jao, jj, jzp, last, last_image, &
         & nao, ncell, ni, nj, nsh, requested_batch_size
      INTEGER(int64) :: peak_expected, stream_compact_fft_peak, &
         & stream_compact_separable_peak
      LOGICAL :: batches_done
      REAL(wp) :: angle, deriv_fd, deriv_ref, energy_large, &
         & energy_large_oracle, energy_mesh, energy_m, energy_p, energy_order, &
         & energy_stream, energy_stream_order, energy_tile, energy_tile_sum, &
         & energy_twist, model_probe, &
         & fft_err, onsite_probe, order_err, separable_err, &
         & selector_fft_err, selector_separable_err, theta_x, theta_y
      REAL(wp) :: energy_sc, err, kfrac(3, 3), kfrac_bad(3, 3), kfrac_order(3, 3), &
         & kfrac_twist(3, 3), lattice_sc(3, 3), step, twist_large(3)
      REAL(wp) :: partial_err, stream_compact_fft_err, &
         & stream_compact_separable_err, stream_reverse_err
      REAL(wp) :: rks_fft_forward(3), rks_fft_reverse(3), &
         & rks_2d_fft_forward(3), rks_2d_fft_reverse(3), &
         & uks_3d_fft_forward(3), uks_3d_fft_reverse(3)
      REAL(wp) :: weights(3), weights_bad(3), weights_order(3), xyz_sc(3, 6), shift(3)
      REAL(wp) :: coord_direction(3, 2), identity(3, 3), lattice0(3, 3), &
         & strain_direction(3, 3), transform(3, 3), xyz0(3, 2)
      COMPLEX(wp) :: batch_sentinel, block_probe, phase
      REAL(wp), ALLOCATABLE :: density_sc(:, :, :), overlap_sc(:, :)
      REAL(wp), ALLOCATABLE :: ao_grad_sc(:, :), bocorr_grad_r(:, :, :), &
         & bocorr_grad_large(:, :, :), bocorr_grad_sc(:, :), &
         & gradient_cp2k(:, :), gradient_large(:, :), gradient_large_oracle(:, :), &
         & gradient_mesh(:, :), gradient_sc(:, :), gradient_stream_probe(:, :), &
         & gradient_stream_reverse(:, :), &
         & kernel_bocorr_backup(:, :, :), kernel_bocorr_direction(:, :, :), &
         & kernel_mulliken_backup(:, :, :), &
         & kernel_mulliken_direction(:, :, :), kfrac_large(:, :), &
         & kfrac_cube(:, :), kfrac_overflow(:, :), gradient_cube(:, :), &
         & gradient_cube_oracle(:, :), &
         & mulliken_grad_large(:, :, :), mulliken_grad_r(:, :, :), &
         & mulliken_grad_sc(:, :), p_r(:, :, :), s_r(:, :, :), &
         & sigma_cp2k(:, :), sigma_large(:, :), sigma_large_oracle(:, :), &
         & sigma_mesh(:, :), sigma_sc(:, :), sigma_stream_probe(:, :), &
         & sigma_stream_reverse(:, :), vsh_cp2k(:), vsh_large(:), &
         & sigma_cube(:, :), sigma_cube_oracle(:, :), weights_cube(:), &
         & weights_overflow(:), &
         & vsh_large_oracle(:), vsh_mesh(:), vsh_order(:), vsh_ref(:), &
         & vsh_stream(:), vsh_stream_order(:), vsh_twist(:), weights_large(:)
      COMPLEX(wp), ALLOCATABLE :: density_mesh(:, :, :, :), &
         & density_m(:, :, :, :), density_p(:, :, :, :), &
         & density_large(:, :, :, :), density_order(:, :, :, :), &
         & density_large_rks(:, :, :, :), density_cube_uks(:, :, :, :), &
         & direction(:, :, :, :), fock_cp2k(:, :, :, :), &
         & fock_large(:, :, :, :), fock_large_oracle(:, :, :, :), &
         & fock_mesh(:, :, :, :), &
         & fock_batch(:, :, :, :), &
         & fock_stream(:, :, :, :), fock_stream_order(:, :, :, :), &
         & fock_order(:, :, :, :), &
         & fock_m(:, :, :, :), fock_p(:, :, :, :), fock_ref(:, :, :, :), &
         & fock_twist(:, :, :, :), overlap_direction(:, :, :), &
         & density_image_large(:, :), overlap_image_large(:, :), &
         & overlap_grad_cp2k(:, :, :), overlap_grad_large(:, :, :), &
         & overlap_grad_large_oracle(:, :, :), overlap_grad_mesh(:, :, :), &
         & overlap_grad_stream_probe(:, :), overlap_grad_stream_reverse(:, :, :), &
         & overlap_grad_stream_reduced(:, :, :), &
         & overlap_grad_stream_batch(:, :, :), &
         & density_cube(:, :, :, :), overlap_cube(:, :, :), &
         & overlap_grad_cube(:, :, :), overlap_grad_cube_oracle(:, :, :), &
         & overlap_bad_dimensions(:, :), &
         & overlap_large(:, :, :), overlap_m(:, :, :), overlap_mesh(:, :, :), &
         & overlap_order(:, :, :), overlap_p(:, :, :)

      CALL new(mol, num, xyz, lattice=lattice, periodic=periodic)
      fft_err = 0.0_wp
      separable_err = 0.0_wp
      stream_reverse_err = 0.0_wp
      ALLOCATE (qbas)
      CALL make_qvszp_basis(qbas, mol, error)
      IF (ALLOCATED(error)) RETURN
      CALL move_alloc(qbas, bas)
      CALL new_wavefunction(wfn_aux, mol%nat, 0, 0, 1, 0.0_wp, .FALSE.)
      CALL get_eeqbc_charges(mol, error, wfn_aux%qat(:, 1))
      IF (ALLOCATED(error)) RETURN
      CALL bas%update(mol, bcache, .FALSE., wfn_aux=wfn_aux)
      CALL make_exchange_gxtb(exchange, mol, bas)
      CALL exchange%update(mol, cache)
      CALL view(cache, ptr)

      SELECT TYPE (fock => exchange)
      TYPE is (exchange_fock)
         CALL fock%get_bvk_Kmatrix(mol, [1, 1, 1], kernel_one)
         CALL fock%get_bvk_Kmatrix(mol, nmesh, kernel)
      class default
         CALL test_failed(error, "Unexpected exchange implementation")
         RETURN
      END SELECT

      err = MAX(MAXVAL(ABS(kernel_one%g_mulliken_r(:, :, 1) &
         & - ptr%g_mulliken)), MAXVAL(ABS(kernel_one%g_bocorr_r(:, :, 1) &
         & - ptr%g_bocorr)))

      lattice_sc = lattice
      lattice_sc(:, 1) = REAL(nmesh(1), wp)*lattice(:, 1)
      DO icell = 1, nmesh(1)
         shift = REAL(icell - 1, wp)*lattice(:, 1)
         xyz_sc(:, 2*icell - 1:2*icell) = xyz + SPREAD(shift, 2, 2)
      END DO
      CALL new(mol_sc, num_sc, xyz_sc, lattice=lattice_sc, periodic=periodic)
      ALLOCATE (qbas_sc)
      CALL make_qvszp_basis(qbas_sc, mol_sc, error)
      IF (ALLOCATED(error)) RETURN
      CALL move_alloc(qbas_sc, bas_sc)
      CALL new_wavefunction(wfn_aux_sc, mol_sc%nat, 0, 0, 1, 0.0_wp, .FALSE.)
      CALL get_eeqbc_charges(mol_sc, error, wfn_aux_sc%qat(:, 1))
      IF (ALLOCATED(error)) RETURN
      CALL bas_sc%update(mol_sc, bcache_sc, .FALSE., wfn_aux=wfn_aux_sc)
      CALL make_exchange_gxtb(exchange_sc, mol_sc, bas_sc)
      CALL exchange_sc%update(mol_sc, cache_sc)
      CALL view(cache_sc, ptr_sc)

      nsh = bas%nsh
      IF (bas_sc%nsh /= nmesh(1)*nsh) THEN
         CALL test_failed(error, "Unexpected replicated basis layout")
         RETURN
      END IF
      DO icell = 1, PRODUCT(nmesh)
         first = (icell - 1)*nsh + 1
         last = icell*nsh
         err = MAX(err, MAXVAL(ABS(kernel%g_mulliken_r(:, :, icell) &
            & - ptr_sc%g_mulliken(1:nsh, first:last))), &
            & MAXVAL(ABS(kernel%g_bocorr_r(:, :, icell) &
            & - ptr_sc%g_bocorr(1:2, 2*icell - 1:2*icell))))
      END DO
      IF (ANY(kernel%reps /= RESHAPE([0, 0, 0, 1, 0, 0, 2, 0, 0], &
         & SHAPE(kernel%reps)))) err = MAX(err, 1.0_wp)

      CALL check(error, err, 0.0_wp, thr=1.0e-12_wp)
      IF (ALLOCATED(error)) THEN
         CALL test_failed(error, &
            & "Image-resolved exchange kernels do not reproduce the BvK supercell")
         RETURN
      END IF

      ! Construct deliberately non-symmetric real image blocks.  R=2 is the
      ! transpose of R=1, so their Fourier transforms are genuinely complex
      ! Hermitian matrices and the assembled supercell is real symmetric.
      nao = bas%nao
      ncell = PRODUCT(nmesh)
      IF (bas_sc%nao /= ncell*nao) THEN
         CALL test_failed(error, "Unexpected replicated AO layout")
         RETURN
      END IF
      ALLOCATE (p_r(nao, nao, ncell), s_r(nao, nao, ncell), source=0.0_wp)
      DO i = 1, nao
         s_r(i, i, 1) = 1.0_wp + 0.04_wp*REAL(i, wp)
         p_r(i, i, 1) = 0.55_wp + 0.07_wp*REAL(i, wp)
         DO j = 1, i - 1
            s_r(j, i, 1) = 0.03_wp/REAL(i + j, wp)
            s_r(i, j, 1) = s_r(j, i, 1)
            p_r(j, i, 1) = 0.04_wp*COS(REAL(2*i + 3*j, wp))
            p_r(i, j, 1) = p_r(j, i, 1)
         END DO
         DO j = 1, nao
            s_r(i, j, 2) = 0.012_wp*SIN(REAL(3*i + 5*j, wp))
            p_r(i, j, 2) = 0.018_wp*COS(REAL(7*i + 2*j, wp))
         END DO
      END DO
      s_r(:, :, 3) = TRANSPOSE(s_r(:, :, 2))
      p_r(:, :, 3) = TRANSPOSE(p_r(:, :, 2))

      ALLOCATE (overlap_sc(ncell*nao, ncell*nao), &
         & density_sc(ncell*nao, ncell*nao, 1), source=0.0_wp)
      DO acell = 0, ncell - 1
         DO bcell = 0, ncell - 1
            ir = MODULO(bcell - acell, ncell) + 1
            overlap_sc(acell*nao + 1:(acell + 1)*nao, &
               & bcell*nao + 1:(bcell + 1)*nao) = s_r(:, :, ir)
            density_sc(acell*nao + 1:(acell + 1)*nao, &
               & bcell*nao + 1:(bcell + 1)*nao, 1) = p_r(:, :, ir)
         END DO
      END DO

      kfrac = 0.0_wp
      weights = 1.0_wp/REAL(ncell, wp)
      DO ik = 1, ncell
         kfrac(1, ik) = REAL(ik - 1, wp)/REAL(ncell, wp)
      END DO
      ALLOCATE (overlap_mesh(nao, nao, ncell), &
         & density_mesh(nao, nao, 1, ncell), &
         & fock_mesh(nao, nao, 1, ncell), &
         & fock_ref(nao, nao, 1, ncell), source=(0.0_wp, 0.0_wp))
      ALLOCATE (vsh_mesh(nsh), vsh_ref(nsh), source=0.0_wp)
      DO ik = 1, ncell
         DO ir = 1, ncell
            angle = 2.0_wp*pi*kfrac(1, ik)*REAL(ir - 1, wp)
            phase = EXP(CMPLX(0.0_wp, angle, wp))
            overlap_mesh(:, :, ik) = overlap_mesh(:, :, ik) &
               & + phase*s_r(:, :, ir)
            density_mesh(:, :, 1, ik) = density_mesh(:, :, 1, ik) &
               & + phase*p_r(:, :, ir)
         END DO
      END DO

      CALL new_wavefunction(wfn, mol%nat, nsh, nao, 1, 0.0_wp, .FALSE.)
      DO i = 1, nsh
         wfn%qsh(i, 1) = 0.08_wp*SIN(REAL(3*i, wp))
      END DO
      CALL exchange%get_onsite_Kmatrix(mol, wfn, ptr)

      CALL new_wavefunction(wfn_sc, mol_sc%nat, bas_sc%nsh, bas_sc%nao, &
         & 1, 0.0_wp, .FALSE.)
      DO icell = 1, ncell
         first = (icell - 1)*nsh + 1
         last = icell*nsh
         wfn_sc%qsh(first:last, 1) = wfn%qsh(:, 1)
      END DO
      CALL exchange_sc%get_onsite_Kmatrix(mol_sc, wfn_sc, ptr_sc)
      ALLOCATE (ptr_sc%prev_F(bas_sc%nao, bas_sc%nao, 1), &
         & ptr_sc%prev_vsh(bas_sc%nsh, 1))
      CALL exchange_sc%get_KFock(mol_sc, ptr_sc, density_sc, overlap_sc)

      SELECT TYPE (fock => exchange)
      TYPE is (exchange_fock)
         CALL get_KFock_kmesh(fock, mol, ptr, oracle_plan, kernel, kfrac, weights, &
            & density_mesh, overlap_mesh, fock_mesh, vsh_mesh, energy_mesh)
      class default
         CALL test_failed(error, "Unexpected exchange implementation")
         RETURN
      END SELECT

      DO ik = 1, ncell
         DO ir = 1, ncell
            first = (ir - 1)*nao + 1
            last = ir*nao
            angle = 2.0_wp*pi*kfrac(1, ik)*REAL(ir - 1, wp)
            phase = EXP(CMPLX(0.0_wp, angle, wp))
            fock_ref(:, :, 1, ik) = fock_ref(:, :, 1, ik) &
               & + phase*ptr_sc%prev_F(1:nao, first:last, 1)
         END DO
      END DO
      energy_sc = 0.5_wp*SUM(density_sc*ptr_sc%prev_F)/REAL(ncell, wp)
      vsh_ref = 0.0_wp
      DO icell = 1, ncell
         first = (icell - 1)*nsh + 1
         last = icell*nsh
         vsh_ref = vsh_ref + ptr_sc%prev_vsh(first:last, 1)/REAL(ncell, wp)
      END DO

      err = MAX(MAXVAL(ABS(fock_mesh - fock_ref)), ABS(energy_mesh - energy_sc), &
         & MAXVAL(ABS(vsh_mesh - vsh_ref)))
      DO icell = 1, ncell
         first = (icell - 1)*nsh + 1
         last = icell*nsh
         err = MAX(err, MAXVAL(ABS(vsh_mesh - ptr_sc%prev_vsh(first:last, 1))))
      END DO

      ! A perturbation at one non-Gamma block must change the other blocks via
      ! the image kernels, while the BZ energy derivative remains the Fock map.
      ALLOCATE (direction(nao, nao, 1, ncell), &
         & density_m(nao, nao, 1, ncell), density_p(nao, nao, 1, ncell), &
         & fock_m(nao, nao, 1, ncell), fock_p(nao, nao, 1, ncell), &
         & source=(0.0_wp, 0.0_wp))
      DO i = 1, nao
         direction(i, i, 1, 2) = CMPLX(0.01_wp*SIN(REAL(5*i, wp)), 0.0_wp, wp)
         DO j = 1, i - 1
            direction(j, i, 1, 2) = CMPLX( &
               & 0.01_wp*COS(REAL(3*i + j, wp)), &
               & 0.01_wp*SIN(REAL(i + 4*j, wp)), wp)
            direction(i, j, 1, 2) = CONJG(direction(j, i, 1, 2))
         END DO
      END DO
      step = 1.0e-5_wp
      density_m = density_mesh - step*direction
      density_p = density_mesh + step*direction
      SELECT TYPE (fock => exchange)
      TYPE is (exchange_fock)
         CALL get_KFock_kmesh(fock, mol, ptr, oracle_plan, kernel, kfrac, weights, &
            & density_m, overlap_mesh, fock_m, vsh_mesh, energy_m)
         CALL get_KFock_kmesh(fock, mol, ptr, oracle_plan, kernel, kfrac, weights, &
            & density_p, overlap_mesh, fock_p, vsh_mesh, energy_p)
      END SELECT
      deriv_fd = (energy_p - energy_m)/(2.0_wp*step)
      deriv_ref = 0.0_wp
      DO ik = 1, ncell
         deriv_ref = deriv_ref + weights(ik)*REAL(SUM( &
            & CONJG(direction(:, :, :, ik))*fock_mesh(:, :, :, ik)), wp)
      END DO
      err = MAX(err, ABS(deriv_fd - deriv_ref))
      IF (MAXVAL(ABS(fock_p(:, :, :, 1) - fock_m(:, :, :, 1))) &
         & < 1.0e-12_wp) err = MAX(err, 1.0_wp)

      ! Reverse the complete-mesh functional with respect to independently
      ! Hermitian S(k) blocks.  The public response is unweighted, so the BZ
      ! contraction below contains exactly one factor of w_k.
      ALLOCATE (overlap_grad_mesh(nao, nao, ncell), &
         & overlap_direction(nao, nao, ncell), overlap_m(nao, nao, ncell), &
         & overlap_p(nao, nao, ncell), source=(0.0_wp, 0.0_wp))
      ALLOCATE (mulliken_grad_r(nao, nao, ncell), &
         & bocorr_grad_r(nao, nao, ncell), source=0.0_wp)
      DO ik = 1, ncell
         DO i = 1, nao
            overlap_direction(i, i, ik) = CMPLX( &
               & 0.01_wp*SIN(REAL(5*i + 7*ik, wp)), 0.0_wp, wp)
            DO j = 1, i - 1
               overlap_direction(j, i, ik) = CMPLX( &
                  & 0.01_wp*COS(REAL(3*i + j + 2*ik, wp)), &
                  & 0.01_wp*SIN(REAL(i + 4*j + 3*ik, wp)), wp)
               overlap_direction(i, j, ik) = &
                  & CONJG(overlap_direction(j, i, ik))
            END DO
         END DO
      END DO
      SELECT TYPE (fock => exchange)
      TYPE is (exchange_fock)
         CALL get_KGrad_kmesh(fock, mol, ptr, oracle_plan, kernel, kfrac, weights, &
            & density_mesh, overlap_mesh, overlap_grad_mesh, &
            & mulliken_grad_r, bocorr_grad_r)
         overlap_m = overlap_mesh - step*overlap_direction
         overlap_p = overlap_mesh + step*overlap_direction
         CALL get_KFock_kmesh(fock, mol, ptr, oracle_plan, kernel, kfrac, weights, &
            & density_mesh, overlap_m, fock_m, vsh_mesh, energy_m)
         CALL get_KFock_kmesh(fock, mol, ptr, oracle_plan, kernel, kfrac, weights, &
            & density_mesh, overlap_p, fock_p, vsh_mesh, energy_p)
      END SELECT
      deriv_fd = (energy_p - energy_m)/(2.0_wp*step)
      deriv_ref = 0.0_wp
      DO ik = 1, ncell
         deriv_ref = deriv_ref + weights(ik)*REAL(SUM( &
            & CONJG(overlap_grad_mesh(:, :, ik)) &
            & *overlap_direction(:, :, ik)), wp)
         err = MAX(err, MAXVAL(ABS(overlap_grad_mesh(:, :, ik) &
            & - CONJG(TRANSPOSE(overlap_grad_mesh(:, :, ik))))))
      END DO
      err = MAX(err, ABS(deriv_fd - deriv_ref))

      ! Perturb every oriented image kernel independently.  In particular no
      ! R/-R or row/column symmetrization is allowed in this contraction.
      ALLOCATE (kernel_mulliken_backup(nsh, nsh, ncell), &
         & kernel_mulliken_direction(nsh, nsh, ncell), &
         & kernel_bocorr_backup(mol%nat, mol%nat, ncell), &
         & kernel_bocorr_direction(mol%nat, mol%nat, ncell))
      kernel_mulliken_backup = kernel%g_mulliken_r
      kernel_bocorr_backup = kernel%g_bocorr_r
      DO icell = 1, ncell
         DO ish = 1, nsh
            DO jsh = 1, nsh
               kernel_mulliken_direction(ish, jsh, icell) = &
                  & 0.01_wp*SIN(REAL(11*ish + 7*jsh + 5*icell, wp))
            END DO
         END DO
         DO iat = 1, mol%nat
            DO jat = 1, mol%nat
               kernel_bocorr_direction(iat, jat, icell) = &
                  & 0.01_wp*COS(REAL(3*iat + 13*jat + 2*icell, wp))
            END DO
         END DO
      END DO
      kernel%g_mulliken_r = kernel_mulliken_backup &
         & - step*kernel_mulliken_direction
      kernel%g_bocorr_r = kernel_bocorr_backup - step*kernel_bocorr_direction
      SELECT TYPE (fock => exchange)
      TYPE is (exchange_fock)
         CALL get_KFock_kmesh(fock, mol, ptr, oracle_plan, kernel, kfrac, weights, &
            & density_mesh, overlap_mesh, fock_m, vsh_mesh, energy_m)
      END SELECT
      kernel%g_mulliken_r = kernel_mulliken_backup &
         & + step*kernel_mulliken_direction
      kernel%g_bocorr_r = kernel_bocorr_backup + step*kernel_bocorr_direction
      SELECT TYPE (fock => exchange)
      TYPE is (exchange_fock)
         CALL get_KFock_kmesh(fock, mol, ptr, oracle_plan, kernel, kfrac, weights, &
            & density_mesh, overlap_mesh, fock_p, vsh_mesh, energy_p)
      END SELECT
      kernel%g_mulliken_r = kernel_mulliken_backup
      kernel%g_bocorr_r = kernel_bocorr_backup
      deriv_fd = (energy_p - energy_m)/(2.0_wp*step)
      deriv_ref = 0.0_wp
      DO icell = 1, ncell
         DO ish = 1, nsh
            ii = bas%iao_sh(ish)
            ni = bas%nao_sh(ish)
            DO jsh = 1, nsh
               jj = bas%iao_sh(jsh)
               nj = bas%nao_sh(jsh)
               deriv_ref = deriv_ref &
                  & + kernel_mulliken_direction(ish, jsh, icell) &
                  & *SUM(mulliken_grad_r(ii + 1:ii + ni, jj + 1:jj + nj, icell))
            END DO
         END DO
         DO iat = 1, mol%nat
            izp = mol%id(iat)
            DO jat = 1, mol%nat
               jzp = mol%id(jat)
               DO ish = 1, bas%nsh_id(izp)
                  ii = bas%iao_sh(bas%ish_at(iat) + ish)
                  ni = bas%nao_sh(bas%ish_at(iat) + ish)
                  DO jsh = 1, bas%nsh_id(jzp)
                     jj = bas%iao_sh(bas%ish_at(jat) + jsh)
                     nj = bas%nao_sh(bas%ish_at(jat) + jsh)
                     deriv_ref = deriv_ref &
                        & + kernel_bocorr_direction(iat, jat, icell) &
                        & *SUM(bocorr_grad_r(ii + 1:ii + ni, jj + 1:jj + nj, icell))
                  END DO
               END DO
            END DO
         END DO
      END DO
      err = MAX(err, ABS(deriv_fd - deriv_ref))

      ! Contract the oriented image responses with the very same BvK
      ! Wigner-Seitz literals used in the forward kernel construction.
      ALLOCATE (gradient_mesh(3, mol%nat), sigma_mesh(3, 3), source=0.0_wp)
      SELECT TYPE (fock => exchange)
      TYPE is (exchange_fock)
         CALL fock%get_bvk_Kmatrix_derivs(mol, kernel, mulliken_grad_r, &
            & bocorr_grad_r, gradient_mesh, sigma_mesh)
      END SELECT
      DO iat = 1, mol%nat
         DO i = 1, 3
            coord_direction(i, iat) = 0.03_wp*SIN(REAL(7*i + 11*iat, wp))
         END DO
      END DO
      xyz0 = mol%xyz
      mol%xyz = xyz0 - step*coord_direction
      SELECT TYPE (fock => exchange)
      TYPE is (exchange_fock)
         CALL fock%get_bvk_Kmatrix(mol, nmesh, kernel)
         CALL get_KFock_kmesh(fock, mol, ptr, oracle_plan, kernel, kfrac, weights, &
            & density_mesh, overlap_mesh, fock_m, vsh_mesh, energy_m)
      END SELECT
      mol%xyz = xyz0 + step*coord_direction
      SELECT TYPE (fock => exchange)
      TYPE is (exchange_fock)
         CALL fock%get_bvk_Kmatrix(mol, nmesh, kernel)
         CALL get_KFock_kmesh(fock, mol, ptr, oracle_plan, kernel, kfrac, weights, &
            & density_mesh, overlap_mesh, fock_p, vsh_mesh, energy_p)
      END SELECT
      deriv_fd = (energy_p - energy_m)/(2.0_wp*step)
      deriv_ref = SUM(gradient_mesh*coord_direction)
      err = MAX(err, ABS(deriv_fd - deriv_ref))

      identity = 0.0_wp
      DO i = 1, 3
         identity(i, i) = 1.0_wp
      END DO
      strain_direction = RESHAPE([ &
         & 0.07_wp, 0.02_wp, -0.01_wp, &
         & 0.02_wp, -0.04_wp, 0.03_wp, &
         & -0.01_wp, 0.03_wp, 0.05_wp], [3, 3])
      lattice0 = mol%lattice
      transform = identity - step*strain_direction
      mol%xyz = MATMUL(transform, xyz0)
      mol%lattice = MATMUL(transform, lattice0)
      SELECT TYPE (fock => exchange)
      TYPE is (exchange_fock)
         CALL fock%get_bvk_Kmatrix(mol, nmesh, kernel)
         CALL get_KFock_kmesh(fock, mol, ptr, oracle_plan, kernel, kfrac, weights, &
            & density_mesh, overlap_mesh, fock_m, vsh_mesh, energy_m)
      END SELECT
      transform = identity + step*strain_direction
      mol%xyz = MATMUL(transform, xyz0)
      mol%lattice = MATMUL(transform, lattice0)
      SELECT TYPE (fock => exchange)
      TYPE is (exchange_fock)
         CALL fock%get_bvk_Kmatrix(mol, nmesh, kernel)
         CALL get_KFock_kmesh(fock, mol, ptr, oracle_plan, kernel, kfrac, weights, &
            & density_mesh, overlap_mesh, fock_p, vsh_mesh, energy_p)
      END SELECT
      deriv_fd = (energy_p - energy_m)/(2.0_wp*step)
      deriv_ref = SUM(sigma_mesh*strain_direction)
      err = MAX(err, ABS(deriv_fd - deriv_ref))
      mol%xyz = xyz0
      mol%lattice = lattice0
      SELECT TYPE (fock => exchange)
      TYPE is (exchange_fock)
         CALL fock%get_bvk_Kmatrix(mol, nmesh, kernel)
      END SELECT

      ! The real Gamma exchange gradient of the explicit 3x1x1 BvK supercell is
      ! an independent oracle for the image-kernel geometry contraction.
      ALLOCATE (mulliken_grad_sc(bas_sc%nao, bas_sc%nao), &
         & bocorr_grad_sc(bas_sc%nao, bas_sc%nao), &
         & ao_grad_sc(bas_sc%nao, bas_sc%nao), &
         & gradient_sc(3, mol_sc%nat), sigma_sc(3, 3), source=0.0_wp)
      SELECT TYPE (fock_sc => exchange_sc)
      TYPE is (exchange_fock)
         CALL fock_sc%get_KGrad(mol_sc, ptr_sc, density_sc, overlap_sc, &
            & mulliken_grad_sc, bocorr_grad_sc, ao_grad_sc)
         CALL fock_sc%get_mulliken_derivs(mol_sc, ptr_sc, mulliken_grad_sc, &
            & gradient_sc, sigma_sc)
         CALL fock_sc%get_bocorr_derivs(mol_sc, ptr_sc, bocorr_grad_sc, &
            & gradient_sc, sigma_sc)
      END SELECT
      DO icell = 0, ncell - 1
         err = MAX(err, MAXVAL(ABS(gradient_mesh &
            & - gradient_sc(:, 2*icell + 1:2*icell + 2))))
      END DO
      err = MAX(err, MAXVAL(ABS(sigma_mesh - sigma_sc/REAL(ncell, wp))))

      ! A common Monkhorst-Pack twist must retain the variational Fock map, and
      ! neither the result nor the energy may depend on CP2K's k-point ordering.
      ALLOCATE (density_order(nao, nao, 1, ncell), &
         & fock_order(nao, nao, 1, ncell), fock_twist(nao, nao, 1, ncell), &
         & overlap_order(nao, nao, ncell), source=(0.0_wp, 0.0_wp))
      ALLOCATE (vsh_order(nsh), vsh_twist(nsh), source=0.0_wp)
      DO ik = 1, ncell
         kfrac_twist(:, ik) = kfrac(:, ik) + [0.137_wp, -0.073_wp, 0.041_wp]
      END DO
      DO ik = 1, ncell
         kfrac_order(:, ik) = kfrac_twist(:, korder(ik))
         weights_order(ik) = weights(korder(ik))
         density_order(:, :, :, ik) = density_mesh(:, :, :, korder(ik))
         overlap_order(:, :, ik) = overlap_mesh(:, :, korder(ik))
      END DO
      SELECT TYPE (fock => exchange)
      TYPE is (exchange_fock)
         CALL get_KFock_kmesh(fock, mol, ptr, oracle_plan, kernel, kfrac_twist, weights, &
            & density_mesh, overlap_mesh, fock_twist, vsh_twist, energy_twist)
         CALL get_KFock_kmesh(fock, mol, ptr, oracle_plan, kernel, kfrac_order, weights_order, &
            & density_order, overlap_order, fock_order, vsh_order, energy_order)
         CALL get_KFock_kmesh(fock, mol, ptr, oracle_plan, kernel, kfrac_twist, weights, &
            & density_m, overlap_mesh, fock_m, vsh_mesh, energy_m)
         CALL get_KFock_kmesh(fock, mol, ptr, oracle_plan, kernel, kfrac_twist, weights, &
            & density_p, overlap_mesh, fock_p, vsh_mesh, energy_p)
      END SELECT
      deriv_fd = (energy_p - energy_m)/(2.0_wp*step)
      deriv_ref = 0.0_wp
      order_err = 0.0_wp
      DO ik = 1, ncell
         deriv_ref = deriv_ref + weights(ik)*REAL(SUM( &
            & CONJG(direction(:, :, :, ik))*fock_twist(:, :, :, ik)), wp)
         order_err = MAX(order_err, MAXVAL(ABS(fock_order(:, :, :, ik) &
            & - fock_twist(:, :, :, korder(ik)))))
      END DO
      err = MAX(err, ABS(deriv_fd - deriv_ref), ABS(energy_order - energy_twist), &
         & MAXVAL(ABS(vsh_order - vsh_twist)), order_err)

      ! Exercise the public CP2K whole-mesh path and use its persistent plan as
      ! the cache-hit/invalidation oracle.  The first call builds the plan, an
      ! identical second call must reuse it, and a reordered mesh must rebuild
      ! the ordering descriptor without changing the physical result.
      ALLOCATE (calc_cp2k%exchange, source=exchange)

      ! On 32-bit default integers this deliberately aliases Nk when evaluated
      ! by unchecked PRODUCT: 65536*65537 wraps to 65536.  The public stream (and
      ! hence partial-accumulator) boundary must reject the true oversized mesh
      ! before preparing a plan or entering any default-integer index arithmetic.
      ALLOCATE (kfrac_overflow(3, nk_overflow), source=0.0_wp)
      ALLOCATE (weights_overflow(nk_overflow), &
         & source=1.0_wp/REAL(nk_overflow, wp))
      CALL cp2k_exchange_stream_begin(stream_overflow, calc_cp2k, mol, &
         & cache_cp2k, cache_cp2k_plan, wfn, nmesh_overflow, kfrac_overflow, weights_overflow, &
         & error)
      IF (ALLOCATED(error)) THEN
         IF (TRIM(error%message) /= &
            & "g-xTB exchange stream requires a complete regular k mesh") THEN
            DEALLOCATE (error)
            CALL test_failed(error, "Unexpected overflowing-mesh diagnostic")
            RETURN
         END IF
         DEALLOCATE (error)
      ELSE
         CALL test_failed(error, "Overflowing g-xTB exchange mesh was accepted")
         RETURN
      END IF
      DEALLOCATE (kfrac_overflow, weights_overflow)

      ! Range validation must not rely on Fortran logical short-circuiting: the
      ! remaining-count expression is evaluated only after first is in 1:Nk.
      CALL cp2k_exchange_stream_begin(stream_overflow, calc_cp2k, mol, &
         & cache_cp2k, cache_cp2k_plan, wfn, nmesh, kfrac, weights, error, &
         & mode=cp2k_exchange_stream_batched, batch_size=1, &
         & image_first=-HUGE(0), image_count=1)
      IF (ALLOCATED(error)) THEN
         IF (TRIM(error%message) /= &
            & "g-xTB exchange image range is invalid") THEN
            DEALLOCATE (error)
            CALL test_failed(error, "Unexpected extreme image-range diagnostic")
            RETURN
         END IF
         DEALLOCATE (error)
      ELSE
         CALL test_failed(error, "Extreme g-xTB image range was accepted")
         RETURN
      END IF

      ALLOCATE (fock_cp2k(nao, nao, 1, ncell), vsh_cp2k(nsh))
      CALL cp2k_exchange_kmesh(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn, nmesh, &
         & kfrac, weights, density_mesh, overlap_mesh, fock_cp2k, &
         & vsh_cp2k, energy_mesh, error)
      IF (ALLOCATED(error)) RETURN
      CALL view(cache_cp2k, ptr_cp2k)
      IF (.NOT. cache_cp2k_plan%bvk_matches(mol, nmesh, kfrac, weights)) THEN
         err = MAX(err, 1.0_wp)
      END IF
      IF (ANY(cache_cp2k_plan%bvk_input_to_grid /= [1, 2, 3])) THEN
         err = MAX(err, 1.0_wp)
      END IF
      IF (ANY(cache_cp2k_plan%bvk_grid_to_input /= [1, 2, 3])) THEN
         err = MAX(err, 1.0_wp)
      END IF
      err = MAX(err, MAXVAL(ABS(fock_cp2k - fock_ref)), &
         & MAXVAL(ABS(vsh_cp2k - vsh_ref)), ABS(energy_mesh - energy_sc))

      ! The memory-reduced stream must reject duplicate and missing input
      ! blocks without corrupting its transaction.  Completing that transaction
      ! afterwards must reproduce the unchanged whole-mesh oracle exactly.
      ALLOCATE (fock_stream(nao, nao, 1, ncell), &
         & fock_stream_order(nao, nao, 1, ncell), &
         & vsh_stream(nsh), vsh_stream_order(nsh), &
         & gradient_stream_probe(3, mol%nat), sigma_stream_probe(3, 3))
      CALL cp2k_exchange_stream_begin(stream, calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn, &
         & nmesh, kfrac, weights, error)
      IF (ALLOCATED(error)) RETURN
      IF (.NOT. cp2k_exchange_stream_has_full_mesh_storage(stream)) THEN
         err = MAX(err, 1.0_wp)
      END IF
      CALL cp2k_exchange_stream_push(stream, 1, density_mesh(:, :, :, 1), &
         & overlap_mesh(:, :, 1), error)
      IF (ALLOCATED(error)) RETURN
      CALL cp2k_exchange_stream_push(stream, 1, density_mesh(:, :, :, 1), &
         & overlap_mesh(:, :, 1), error)
      IF (ALLOCATED(error)) THEN
         IF (TRIM(error%message) /= "duplicate g-xTB exchange stream block") THEN
            DEALLOCATE (error)
            CALL test_failed(error, "Unexpected duplicate stream-block diagnostic")
            RETURN
         END IF
         DEALLOCATE (error)
      ELSE
         CALL test_failed(error, "Duplicate exchange stream block was accepted")
         RETURN
      END IF
      CALL cp2k_exchange_stream_push(stream, 2, density_mesh(:, :, :, 2), &
         & overlap_mesh(:, :, 2), error)
      IF (ALLOCATED(error)) RETURN
      CALL cp2k_exchange_stream_apply(stream, calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, &
         & vsh_stream, energy_stream, error)
      IF (ALLOCATED(error)) THEN
         IF (TRIM(error%message) /= "g-xTB exchange stream has missing blocks") THEN
            DEALLOCATE (error)
            CALL test_failed(error, "Unexpected missing stream-block diagnostic")
            RETURN
         END IF
         DEALLOCATE (error)
      ELSE
         CALL test_failed(error, "Incomplete exchange stream was applied")
         RETURN
      END IF
      CALL cp2k_exchange_stream_push(stream, 3, density_mesh(:, :, :, 3), &
         & overlap_mesh(:, :, 3), error)
      IF (ALLOCATED(error)) RETURN
      onsite_probe = ptr_cp2k%g_onsfx(1, 1, 1)
      ptr_cp2k%g_onsfx(1, 1, 1) = onsite_probe + 1.0e-8_wp
      CALL cp2k_exchange_stream_apply(stream, calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, &
         & vsh_stream, energy_stream, error)
      IF (ALLOCATED(error)) THEN
         IF (TRIM(error%message) /= &
            & "g-xTB exchange stream state was invalidated before apply") THEN
            DEALLOCATE (error)
            CALL test_failed(error, "Unexpected invalidated stream diagnostic")
            RETURN
         END IF
         DEALLOCATE (error)
      ELSE
         CALL test_failed(error, "Invalidated exchange stream was applied")
         RETURN
      END IF
      ptr_cp2k%g_onsfx(1, 1, 1) = onsite_probe
      block_probe = cache_cp2k_plan%bvk_phase_forward(1, 1)
      cache_cp2k_plan%bvk_phase_forward(1, 1) = block_probe + CMPLX(0.0_wp, 1.0e-8_wp, wp)
      CALL cp2k_exchange_stream_apply(stream, calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, &
         & vsh_stream, energy_stream, error)
      IF (ALLOCATED(error)) THEN
         IF (TRIM(error%message) /= &
            & "g-xTB exchange stream state was invalidated before apply") THEN
            DEALLOCATE (error)
            CALL test_failed(error, "Unexpected invalidated phase diagnostic")
            RETURN
         END IF
         DEALLOCATE (error)
      ELSE
         CALL test_failed(error, "Invalidated exchange phase was applied")
         RETURN
      END IF
      cache_cp2k_plan%bvk_phase_forward(1, 1) = block_probe
      CALL cp2k_exchange_stream_apply(stream, calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, &
         & vsh_stream, energy_stream, error)
      IF (ALLOCATED(error)) RETURN
      IF (.NOT. cp2k_exchange_stream_has_full_mesh_storage(stream)) THEN
         err = MAX(err, 1.0_wp)
      END IF
      DO ik = 1, ncell
         CALL cp2k_exchange_stream_get(stream, ik, &
            & fock_stream(:, :, :, ik), error, overlap_mesh(:, :, ik))
         IF (ALLOCATED(error)) RETURN
      END DO
      CALL cp2k_exchange_stream_get(stream, 1, &
         & fock_stream(:, :, :, 1), error, overlap_mesh(:, :, 1))
      IF (ALLOCATED(error)) THEN
         IF (TRIM(error%message) /= &
            & "duplicate g-xTB exchange stream forward pull") THEN
            DEALLOCATE (error)
            CALL test_failed(error, "Unexpected duplicate forward-pull diagnostic")
            RETURN
         END IF
         DEALLOCATE (error)
      ELSE
         CALL test_failed(error, "Exchange Fock block was pulled twice")
         RETURN
      END IF
      CALL cp2k_exchange_stream_end(stream, error)
      IF (ALLOCATED(error)) RETURN
      err = MAX(err, MAXVAL(ABS(fock_stream - fock_cp2k)), &
         & MAXVAL(ABS(vsh_stream - vsh_cp2k)), ABS(energy_stream - energy_mesh))

      ! Closing a forward transaction requires exactly one pull of every block;
      ! the failed close must clear the object so it can be reused immediately.
      CALL cp2k_exchange_stream_begin(stream_forward_missing, calc_cp2k, mol, &
         & cache_cp2k, cache_cp2k_plan, wfn, nmesh, kfrac, weights, error)
      IF (ALLOCATED(error)) RETURN
      DO ik = 1, ncell
         CALL cp2k_exchange_stream_push(stream_forward_missing, ik, &
            & density_mesh(:, :, :, ik), overlap_mesh(:, :, ik), error)
         IF (ALLOCATED(error)) RETURN
      END DO
      CALL cp2k_exchange_stream_apply(stream_forward_missing, calc_cp2k, mol, &
         & cache_cp2k, cache_cp2k_plan, vsh_stream, energy_stream, error)
      IF (ALLOCATED(error)) RETURN
      DO ik = 1, ncell - 1
         CALL cp2k_exchange_stream_get(stream_forward_missing, ik, &
            & fock_stream(:, :, :, ik), error, overlap_mesh(:, :, ik))
         IF (ALLOCATED(error)) RETURN
      END DO
      CALL cp2k_exchange_stream_end(stream_forward_missing, error)
      IF (ALLOCATED(error)) THEN
         IF (TRIM(error%message) /= &
            & "g-xTB exchange stream ended with missing forward pulls") THEN
            DEALLOCATE (error)
            CALL test_failed(error, "Unexpected missing forward-pull diagnostic")
            RETURN
         END IF
         DEALLOCATE (error)
      ELSE
         CALL test_failed(error, "Exchange stream closed with a missing Fock pull")
         RETURN
      END IF
      CALL cp2k_exchange_stream_begin(stream_forward_missing, calc_cp2k, mol, &
         & cache_cp2k, cache_cp2k_plan, wfn, nmesh, kfrac, weights, error)
      IF (ALLOCATED(error)) RETURN
      DO ik = 1, ncell
         CALL cp2k_exchange_stream_push(stream_forward_missing, ik, &
            & density_mesh(:, :, :, ik), overlap_mesh(:, :, ik), error)
         IF (ALLOCATED(error)) RETURN
      END DO
      CALL cp2k_exchange_stream_apply(stream_forward_missing, calc_cp2k, mol, &
         & cache_cp2k, cache_cp2k_plan, vsh_stream, energy_stream, error)
      IF (ALLOCATED(error)) RETURN
      DO ik = 1, ncell
         CALL cp2k_exchange_stream_get(stream_forward_missing, ik, &
            & fock_stream(:, :, :, ik), error, overlap_mesh(:, :, ik))
         IF (ALLOCATED(error)) RETURN
      END DO
      CALL cp2k_exchange_stream_end(stream_forward_missing, error)
      IF (ALLOCATED(error)) RETURN

      ! Arrival order is independent of the physical k-point ordering.  Push the
      ! same indexed blocks in a nontrivial order and compare the full response.
      CALL cp2k_exchange_stream_begin(stream_order, calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, &
         & wfn, nmesh, kfrac, weights, error)
      IF (ALLOCATED(error)) RETURN
      IF (.NOT. cp2k_exchange_stream_has_full_mesh_storage(stream_order)) THEN
         err = MAX(err, 1.0_wp)
      END IF
      DO j = 1, ncell
         ik = korder(j)
         CALL cp2k_exchange_stream_push(stream_order, ik, &
            & density_mesh(:, :, :, ik), overlap_mesh(:, :, ik), error)
         IF (ALLOCATED(error)) RETURN
      END DO
      CALL cp2k_exchange_stream_apply(stream_order, calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, &
         & vsh_stream_order, energy_stream_order, error)
      IF (ALLOCATED(error)) RETURN
      CALL cp2k_exchange_stream_reverse_apply(stream_order, calc_cp2k, mol, &
         & cache_cp2k, cache_cp2k_plan, gradient_stream_probe, sigma_stream_probe, error)
      IF (ALLOCATED(error)) THEN
         IF (TRIM(error%message) /= &
            & "g-xTB exchange stream reverse requires oracle mode") THEN
            DEALLOCATE (error)
            CALL test_failed(error, "Unexpected reduced reverse diagnostic")
            RETURN
         END IF
         DEALLOCATE (error)
      ELSE
         CALL test_failed(error, "Reduced exchange stream accepted reverse apply")
         RETURN
      END IF
      DO ik = 1, ncell
         CALL cp2k_exchange_stream_get(stream_order, ik, &
            & fock_stream_order(:, :, :, ik), error, overlap_mesh(:, :, ik))
         IF (ALLOCATED(error)) RETURN
      END DO
      CALL cp2k_exchange_stream_end(stream_order, error)
      IF (ALLOCATED(error)) RETURN
      err = MAX(err, MAXVAL(ABS(fock_stream_order - fock_stream)), &
         & MAXVAL(ABS(vsh_stream_order - vsh_stream)), &
         & ABS(energy_stream_order - energy_stream))

      ! A common twist and a permutation of the physical k-point list must use
      ! the cached representatives without changing the reduced-stream result.
      CALL cp2k_exchange_stream_begin(stream_permuted, calc_cp2k, mol, &
         & cache_cp2k, cache_cp2k_plan, wfn, nmesh, kfrac_order, weights_order, error, &
         & mode=cp2k_exchange_stream_reduced)
      IF (ALLOCATED(error)) RETURN
      IF (.NOT. cp2k_exchange_stream_has_full_mesh_storage(stream_permuted)) THEN
         err = MAX(err, 1.0_wp)
      END IF
      DO ik = 1, ncell
         CALL cp2k_exchange_stream_push(stream_permuted, ik, &
            & density_order(:, :, :, ik), overlap_order(:, :, ik), error)
         IF (ALLOCATED(error)) RETURN
      END DO
      CALL cp2k_exchange_stream_apply(stream_permuted, calc_cp2k, mol, &
         & cache_cp2k, cache_cp2k_plan, vsh_stream_order, energy_stream_order, error)
      IF (ALLOCATED(error)) RETURN
      DO ik = 1, ncell
         CALL cp2k_exchange_stream_get(stream_permuted, ik, &
            & fock_stream_order(:, :, :, ik), error, overlap_order(:, :, ik))
         IF (ALLOCATED(error)) RETURN
      END DO
      CALL cp2k_exchange_stream_end(stream_permuted, error)
      IF (ALLOCATED(error)) RETURN
      err = MAX(err, MAXVAL(ABS(fock_stream_order - fock_order)), &
         & MAXVAL(ABS(vsh_stream_order - vsh_order)), &
         & ABS(energy_stream_order - energy_order))

      ! The image-batched protocol must retain only a caller-selected number of
      ! AO image matrices.  It deliberately requires the sequential caller to
      ! resubmit all k blocks and to pull every k contribution for every batch.
      batch_sentinel = CMPLX(0.037_wp, -0.019_wp, wp)
      ALLOCATE (fock_batch(nao, nao, 1, ncell), source=batch_sentinel)
      ALLOCATE (overlap_bad_dimensions(nao + 1, nao), source=(0.0_wp, 0.0_wp))
      CALL cp2k_exchange_stream_begin(stream_batch, calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, &
         & wfn, nmesh, kfrac_order, weights_order, error, &
         & mode=cp2k_exchange_stream_batched, batch_size=0)
      IF (ALLOCATED(error)) THEN
         IF (TRIM(error%message) /= &
            & "g-xTB exchange stream batch size must be positive") THEN
            DEALLOCATE (error)
            CALL test_failed(error, "Unexpected invalid batch-size diagnostic")
            RETURN
         END IF
         DEALLOCATE (error)
      ELSE
         CALL test_failed(error, "Zero exchange image-batch size was accepted")
         RETURN
      END IF

      CALL cp2k_exchange_stream_begin(stream_batch, calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, &
         & wfn, nmesh, kfrac_order, weights_order, error, &
         & mode=cp2k_exchange_stream_batched, batch_size=2, image_first=1)
      IF (ALLOCATED(error)) THEN
         IF (TRIM(error%message) /= &
            & "g-xTB exchange image range requires first and count") THEN
            DEALLOCATE (error)
            CALL test_failed(error, "Unexpected incomplete image-range diagnostic")
            RETURN
         END IF
         DEALLOCATE (error)
      ELSE
         CALL test_failed(error, "Incomplete exchange image range was accepted")
         RETURN
      END IF

      CALL cp2k_exchange_stream_begin(stream_batch, calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, &
         & wfn, nmesh, kfrac_order, weights_order, error, &
         & mode=cp2k_exchange_stream_batched, batch_size=2, &
         & image_first=0, image_count=1)
      IF (ALLOCATED(error)) THEN
         IF (TRIM(error%message) /= "g-xTB exchange image range is invalid") THEN
            DEALLOCATE (error)
            CALL test_failed(error, "Unexpected invalid image-range diagnostic")
            RETURN
         END IF
         DEALLOCATE (error)
      ELSE
         CALL test_failed(error, "Invalid exchange image range was accepted")
         RETURN
      END IF

      kfrac_bad = kfrac_order
      kfrac_bad(1, 1) = kfrac_bad(1, 1) + 2.0e-11_wp
      CALL cp2k_exchange_stream_begin(stream_batch, calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, &
         & wfn, nmesh, kfrac_bad, weights_order, error, &
         & mode=cp2k_exchange_stream_batched, batch_size=2)
      IF (ALLOCATED(error)) THEN
         IF (TRIM(error%message) /= &
            & "batched g-xTB exchange requires a machine-regular grid") THEN
            DEALLOCATE (error)
            CALL test_failed(error, "Unexpected inexact batch-grid diagnostic")
            RETURN
         END IF
         DEALLOCATE (error)
      ELSE
         CALL test_failed(error, "Inexact exchange image-batch grid was accepted")
         RETURN
      END IF

      weights_bad = weights_order
      weights_bad(1) = weights_bad(1) + 5.0e-13_wp
      weights_bad(2) = weights_bad(2) - 5.0e-13_wp
      CALL cp2k_exchange_stream_begin(stream_batch, calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, &
         & wfn, nmesh, kfrac_order, weights_bad, error, &
         & mode=cp2k_exchange_stream_batched, batch_size=2)
      IF (ALLOCATED(error)) THEN
         IF (TRIM(error%message) /= &
            & "batched g-xTB exchange requires machine-uniform weights") THEN
            DEALLOCATE (error)
            CALL test_failed(error, "Unexpected inexact batch-weight diagnostic")
            RETURN
         END IF
         DEALLOCATE (error)
      ELSE
         CALL test_failed(error, "Inexact exchange image-batch weights were accepted")
         RETURN
      END IF

      CALL cp2k_exchange_stream_begin(stream_batch, calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, &
         & wfn, nmesh, kfrac_order, weights_order, error, &
         & mode=cp2k_exchange_stream_batched, batch_size=2)
      IF (ALLOCATED(error)) RETURN
      IF (cp2k_exchange_stream_has_full_mesh_storage(stream_batch)) THEN
         err = MAX(err, 1.0_wp)
      END IF
      IF (ALLOCATED(cache_cp2k_plan%bvk_phase_forward) &
         & .OR. ALLOCATED(cache_cp2k_plan%bvk_phase_inverse)) err = MAX(err, 1.0_wp)
      IF (.NOT. ALLOCATED(cache_cp2k_plan%bvk_phase_roots) &
         & .OR. .NOT. ALLOCATED(cache_cp2k_plan%bvk_twist_phase)) THEN
         err = MAX(err, 1.0_wp)
      ELSE
         IF (ANY(SHAPE(cache_cp2k_plan%bvk_phase_roots) /= [MAXVAL(nmesh), 3]) &
            & .OR. SIZE(cache_cp2k_plan%bvk_twist_phase) /= ncell) THEN
            err = MAX(err, 1.0_wp)
         END IF
      END IF
      peak_expected = 4_int64*INT(nao, int64)*INT(nao, int64) &
         & *2_int64 + 3_int64*INT(nao, int64) &
         & + 2_int64*INT(ncell, int64) &
         & + 3_int64*INT(MAXVAL(nmesh), int64)
      peak_expected = peak_expected + MAX( &
         & 2_int64*INT(nao, int64)**2, &
         & 5_int64*INT(nao, int64)**2 + 6_int64*INT(nao, int64), &
         & 5_int64*INT(nao, int64)**2 + 2_int64)
      IF (cp2k_exchange_stream_peak_complex_elements(stream_batch) &
         & /= peak_expected) err = MAX(err, 1.0_wp)

      CALL cp2k_exchange_stream_batch_pull(stream_batch, 1, &
         & overlap_order(:, :, 1), fock_batch(:, :, :, 1), error)
      IF (ALLOCATED(error)) THEN
         IF (TRIM(error%message) /= &
            & "g-xTB exchange image batch was not applied") THEN
            DEALLOCATE (error)
            CALL test_failed(error, "Unexpected pre-apply batch-pull diagnostic")
            RETURN
         END IF
         DEALLOCATE (error)
      ELSE
         CALL test_failed(error, "Exchange image batch was pulled before apply")
         RETURN
      END IF
      CALL cp2k_exchange_stream_batch_advance(stream_batch, batches_done, error)
      IF (ALLOCATED(error)) THEN
         IF (TRIM(error%message) /= &
            & "g-xTB exchange image batch was not applied") THEN
            DEALLOCATE (error)
            CALL test_failed(error, "Unexpected pre-apply batch-advance diagnostic")
            RETURN
         END IF
         DEALLOCATE (error)
      ELSE
         CALL test_failed(error, "Exchange image batch advanced before apply")
         RETURN
      END IF
      CALL cp2k_exchange_stream_batch_result(stream_batch, calc_cp2k, mol, &
         & cache_cp2k, cache_cp2k_plan, vsh_stream_order, energy_stream_order, error)
      IF (ALLOCATED(error)) THEN
         IF (TRIM(error%message) /= &
            & "g-xTB exchange image batches are incomplete") THEN
            DEALLOCATE (error)
            CALL test_failed(error, "Unexpected incomplete batch-result diagnostic")
            RETURN
         END IF
         DEALLOCATE (error)
      ELSE
         CALL test_failed(error, "Incomplete exchange image batches returned a result")
         RETURN
      END IF

      batches_done = .FALSE.
      DO WHILE (.NOT. batches_done)
         CALL cp2k_exchange_stream_batch_bounds(stream_batch, first_image, &
            & last_image, error)
         IF (ALLOCATED(error)) RETURN
         IF (last_image - first_image + 1 > 2) err = MAX(err, 1.0_wp)

         IF (first_image == 1) THEN
            CALL cp2k_exchange_stream_batch_apply(stream_batch, calc_cp2k, mol, &
               & cache_cp2k, cache_cp2k_plan, error)
            IF (ALLOCATED(error)) THEN
               IF (TRIM(error%message) /= &
                  & "g-xTB exchange image batch has missing blocks") THEN
                  DEALLOCATE (error)
                  CALL test_failed(error, &
                     & "Unexpected missing image-batch block diagnostic")
                  RETURN
               END IF
               DEALLOCATE (error)
            ELSE
               CALL test_failed(error, "Incomplete exchange image batch applied")
               RETURN
            END IF
         END IF

         IF (first_image > 1) THEN
            ik = korder(1)
            block_probe = density_order(1, 1, 1, ik)
            density_order(1, 1, 1, ik) = block_probe + 1.0e-8_wp
            CALL cp2k_exchange_stream_push(stream_batch, ik, &
               & density_order(:, :, :, ik), overlap_order(:, :, ik), error)
            density_order(1, 1, 1, ik) = block_probe
            IF (ALLOCATED(error)) THEN
               IF (TRIM(error%message) /= &
                  & "changed g-xTB exchange image-batch density block") THEN
                  DEALLOCATE (error)
                  CALL test_failed(error, "Unexpected changed-density diagnostic")
                  RETURN
               END IF
               DEALLOCATE (error)
            ELSE
               CALL test_failed(error, "Changed image-batch density was accepted")
               RETURN
            END IF

            block_probe = overlap_order(1, 1, ik)
            overlap_order(1, 1, ik) = block_probe + 1.0e-8_wp
            CALL cp2k_exchange_stream_push(stream_batch, ik, &
               & density_order(:, :, :, ik), overlap_order(:, :, ik), error)
            overlap_order(1, 1, ik) = block_probe
            IF (ALLOCATED(error)) THEN
               IF (TRIM(error%message) /= &
                  & "changed g-xTB exchange image-batch overlap block") THEN
                  DEALLOCATE (error)
                  CALL test_failed(error, "Unexpected changed-overlap diagnostic")
                  RETURN
               END IF
               DEALLOCATE (error)
            ELSE
               CALL test_failed(error, "Changed image-batch overlap was accepted")
               RETURN
            END IF
         END IF

         DO j = 1, ncell
            ik = korder(j)
            CALL cp2k_exchange_stream_push(stream_batch, ik, &
               & density_order(:, :, :, ik), overlap_order(:, :, ik), error)
            IF (ALLOCATED(error)) RETURN
            IF (first_image == 1 .AND. j == 1) THEN
               CALL cp2k_exchange_stream_push(stream_batch, ik, &
                  & density_order(:, :, :, ik), overlap_order(:, :, ik), error)
               IF (ALLOCATED(error)) THEN
                  IF (TRIM(error%message) /= &
                     & "duplicate g-xTB exchange stream block") THEN
                     DEALLOCATE (error)
                     CALL test_failed(error, &
                        & "Unexpected duplicate image-batch block diagnostic")
                     RETURN
                  END IF
                  DEALLOCATE (error)
               ELSE
                  CALL test_failed(error, &
                     & "Duplicate exchange image-batch block was accepted")
                  RETURN
               END IF
            END IF
         END DO

         IF (first_image == 1) THEN
            SELECT TYPE (fock_cp2k_model => calc_cp2k%exchange)
            TYPE is (exchange_fock)
               model_probe = fock_cp2k_model%onecxints(1, 1, 1)
               fock_cp2k_model%onecxints(1, 1, 1) = model_probe + 1.0e-8_wp
            END SELECT
            CALL cp2k_exchange_stream_batch_apply(stream_batch, calc_cp2k, mol, &
               & cache_cp2k, cache_cp2k_plan, error)
            IF (ALLOCATED(error)) THEN
               IF (TRIM(error%message) /= &
                  & "g-xTB exchange image batch state was invalidated before apply") THEN
                  DEALLOCATE (error)
                  CALL test_failed(error, &
                     & "Unexpected stale image-batch model diagnostic")
                  RETURN
               END IF
               DEALLOCATE (error)
            ELSE
               CALL test_failed(error, "Stale exchange image batch was applied")
               RETURN
            END IF
            SELECT TYPE (fock_cp2k_model => calc_cp2k%exchange)
            TYPE is (exchange_fock)
               fock_cp2k_model%onecxints(1, 1, 1) = model_probe
            END SELECT
            SELECT TYPE (fock_cp2k_model => calc_cp2k%exchange)
            TYPE is (exchange_fock)
               model_probe = fock_cp2k_model%kq(1, 1)
               fock_cp2k_model%kq(1, 1) = model_probe + 1.0e-8_wp
            END SELECT
            CALL cp2k_exchange_stream_batch_apply(stream_batch, calc_cp2k, mol, &
               & cache_cp2k, cache_cp2k_plan, error)
            IF (ALLOCATED(error)) THEN
               IF (TRIM(error%message) /= &
                  & "g-xTB exchange image batch state was invalidated before apply") THEN
                  DEALLOCATE (error)
                  CALL test_failed(error, &
                     & "Unexpected stale image-batch kq diagnostic")
                  RETURN
               END IF
               DEALLOCATE (error)
            ELSE
               CALL test_failed(error, "Stale kq exchange image batch was applied")
               RETURN
            END IF
            SELECT TYPE (fock_cp2k_model => calc_cp2k%exchange)
            TYPE is (exchange_fock)
               fock_cp2k_model%kq(1, 1) = model_probe
            END SELECT

            ! The newer model separates onsite and intersite exchange scaling.
            ! An active transaction must reject a changed onsite model too.
            SELECT TYPE (fock_cp2k_model => calc_cp2k%exchange)
            TYPE is (exchange_fock)
               model_probe = fock_cp2k_model%onsite_frscale
               fock_cp2k_model%onsite_frscale = model_probe + 1.0e-8_wp
            END SELECT
            CALL cp2k_exchange_stream_batch_apply(stream_batch, calc_cp2k, mol, &
               & cache_cp2k, cache_cp2k_plan, error)
            SELECT TYPE (fock_cp2k_model => calc_cp2k%exchange)
            TYPE is (exchange_fock)
               fock_cp2k_model%onsite_frscale = model_probe
            END SELECT
            IF (ALLOCATED(error)) THEN
               IF (TRIM(error%message) /= &
                  & "g-xTB exchange image batch state was invalidated before apply") THEN
                  DEALLOCATE (error)
                  CALL test_failed(error, "Unexpected stale onsite-model diagnostic")
                  RETURN
               END IF
               DEALLOCATE (error)
            ELSE
               CALL test_failed(error, "Stale onsite exchange image batch was applied")
               RETURN
            END IF

            model_probe = cache_cp2k_plan%bvk_kernel%g_mulliken_r(1, 1, 1)
            cache_cp2k_plan%bvk_kernel%g_mulliken_r(1, 1, 1) = model_probe + 1.0e-8_wp
            CALL cp2k_exchange_stream_batch_apply(stream_batch, calc_cp2k, mol, &
               & cache_cp2k, cache_cp2k_plan, error)
            IF (ALLOCATED(error)) THEN
               IF (TRIM(error%message) /= &
                  & "g-xTB exchange image batch state was invalidated before apply") THEN
                  DEALLOCATE (error)
                  CALL test_failed(error, "Unexpected stale batch-kernel diagnostic")
                  RETURN
               END IF
               DEALLOCATE (error)
            ELSE
               CALL test_failed(error, "Stale exchange image-batch kernel was applied")
               RETURN
            END IF
            cache_cp2k_plan%bvk_kernel%g_mulliken_r(1, 1, 1) = model_probe

            block_probe = cache_cp2k_plan%bvk_phase_roots(1, 1)
            cache_cp2k_plan%bvk_phase_roots(1, 1) = block_probe &
               & + CMPLX(0.0_wp, 1.0e-8_wp, wp)
            CALL cp2k_exchange_stream_batch_apply(stream_batch, calc_cp2k, mol, &
               & cache_cp2k, cache_cp2k_plan, error)
            IF (ALLOCATED(error)) THEN
               IF (TRIM(error%message) /= &
                  & "g-xTB exchange image batch state was invalidated before apply") THEN
                  DEALLOCATE (error)
                  CALL test_failed(error, &
                     & "Unexpected stale compact-phase diagnostic")
                  RETURN
               END IF
               DEALLOCATE (error)
            ELSE
               CALL test_failed(error, "Stale compact exchange phase was applied")
               RETURN
            END IF
            cache_cp2k_plan%bvk_phase_roots(1, 1) = block_probe

            i = cache_cp2k_plan%bvk_input_to_grid(1)
            j = cache_cp2k_plan%bvk_input_to_grid(2)
            cache_cp2k_plan%bvk_input_to_grid(1) = j
            cache_cp2k_plan%bvk_input_to_grid(2) = i
            cache_cp2k_plan%bvk_grid_to_input(j) = 1
            cache_cp2k_plan%bvk_grid_to_input(i) = 2
            CALL cp2k_exchange_stream_batch_apply(stream_batch, calc_cp2k, mol, &
               & cache_cp2k, cache_cp2k_plan, error)
            IF (ALLOCATED(error)) THEN
               IF (TRIM(error%message) /= &
                  & "g-xTB exchange image batch state was invalidated before apply") THEN
                  DEALLOCATE (error)
                  CALL test_failed(error, "Unexpected stale batch-order diagnostic")
                  RETURN
               END IF
               DEALLOCATE (error)
            ELSE
               CALL test_failed(error, "Reordered exchange image batch was applied")
               RETURN
            END IF
            cache_cp2k_plan%bvk_input_to_grid(1) = i
            cache_cp2k_plan%bvk_input_to_grid(2) = j
            cache_cp2k_plan%bvk_grid_to_input(i) = 1
            cache_cp2k_plan%bvk_grid_to_input(j) = 2
         END IF

         CALL cp2k_exchange_stream_batch_apply(stream_batch, calc_cp2k, mol, &
            & cache_cp2k, cache_cp2k_plan, error)
         IF (ALLOCATED(error)) RETURN
         IF (first_image == 1) THEN
            CALL cp2k_exchange_stream_batch_apply(stream_batch, calc_cp2k, mol, &
               & cache_cp2k, cache_cp2k_plan, error)
            IF (ALLOCATED(error)) THEN
               IF (TRIM(error%message) /= &
                  & "g-xTB exchange image batch was already applied") THEN
                  DEALLOCATE (error)
                  CALL test_failed(error, &
                     & "Unexpected duplicate image-batch apply diagnostic")
                  RETURN
               END IF
               DEALLOCATE (error)
            ELSE
               CALL test_failed(error, "Exchange image batch was applied twice")
               RETURN
            END IF
            CALL cp2k_exchange_stream_batch_pull(stream_batch, korder(1), &
               & overlap_bad_dimensions, fock_batch(:, :, :, korder(1)), error)
            IF (ALLOCATED(error)) THEN
               IF (TRIM(error%message) /= &
                  & "inconsistent g-xTB exchange image-batch pull dimensions") THEN
                  DEALLOCATE (error)
                  CALL test_failed(error, &
                     & "Unexpected image-batch pull-dimension diagnostic")
                  RETURN
               END IF
               DEALLOCATE (error)
            ELSE
               CALL test_failed(error, &
                  & "Wrong-dimensional exchange image-batch pull was accepted")
               RETURN
            END IF
         END IF
         IF (first_image > 1) THEN
            ik = korder(1)
            block_probe = overlap_order(1, 1, ik)
            overlap_order(1, 1, ik) = block_probe + 1.0e-8_wp
            CALL cp2k_exchange_stream_batch_pull(stream_batch, ik, &
               & overlap_order(:, :, ik), fock_batch(:, :, :, ik), error)
            overlap_order(1, 1, ik) = block_probe
            IF (ALLOCATED(error)) THEN
               IF (TRIM(error%message) /= &
                  & "changed g-xTB exchange image-batch pull overlap") THEN
                  DEALLOCATE (error)
                  CALL test_failed(error, "Unexpected changed-pull diagnostic")
                  RETURN
               END IF
               DEALLOCATE (error)
            ELSE
               CALL test_failed(error, "Changed image-batch pull was accepted")
               RETURN
            END IF
         END IF
         CALL cp2k_exchange_stream_batch_pull(stream_batch, korder(1), &
            & overlap_order(:, :, korder(1)), &
            & fock_batch(:, :, :, korder(1)), error)
         IF (ALLOCATED(error)) RETURN
         IF (first_image == 1) THEN
            CALL cp2k_exchange_stream_batch_pull(stream_batch, korder(1), &
               & overlap_order(:, :, korder(1)), &
               & fock_batch(:, :, :, korder(1)), error)
            IF (ALLOCATED(error)) THEN
               IF (TRIM(error%message) /= &
                  & "duplicate g-xTB exchange image-batch pull") THEN
                  DEALLOCATE (error)
                  CALL test_failed(error, &
                     & "Unexpected duplicate image-batch pull diagnostic")
                  RETURN
               END IF
               DEALLOCATE (error)
            ELSE
               CALL test_failed(error, &
                  & "Duplicate exchange image-batch pull was accepted")
               RETURN
            END IF
            CALL cp2k_exchange_stream_batch_advance(stream_batch, batches_done, &
               & error)
            IF (ALLOCATED(error)) THEN
               IF (TRIM(error%message) /= &
                  & "g-xTB exchange image batch has missing pulls") THEN
                  DEALLOCATE (error)
                  CALL test_failed(error, &
                     & "Unexpected missing image-batch pull diagnostic")
                  RETURN
               END IF
               DEALLOCATE (error)
            ELSE
               CALL test_failed(error, &
                  & "Exchange image batch advanced with missing pulls")
               RETURN
            END IF
         END IF
         DO j = 2, ncell
            ik = korder(j)
            CALL cp2k_exchange_stream_batch_pull(stream_batch, ik, &
               & overlap_order(:, :, ik), &
               & fock_batch(:, :, :, ik), error)
            IF (ALLOCATED(error)) RETURN
         END DO
         CALL cp2k_exchange_stream_batch_advance(stream_batch, batches_done, error)
         IF (ALLOCATED(error)) RETURN
      END DO
      CALL cp2k_exchange_stream_batch_bounds(stream_batch, first_image, &
         & last_image, error)
      IF (ALLOCATED(error)) THEN
         IF (TRIM(error%message) /= &
            & "g-xTB exchange image batches are complete") THEN
            DEALLOCATE (error)
            CALL test_failed(error, "Unexpected completed batch-bounds diagnostic")
            RETURN
         END IF
         DEALLOCATE (error)
      ELSE
         CALL test_failed(error, "Completed exchange image batches returned bounds")
         RETURN
      END IF
      CALL cp2k_exchange_stream_batch_advance(stream_batch, batches_done, error)
      IF (ALLOCATED(error)) THEN
         IF (TRIM(error%message) /= &
            & "g-xTB exchange image batches are complete") THEN
            DEALLOCATE (error)
            CALL test_failed(error, "Unexpected completed batch-advance diagnostic")
            RETURN
         END IF
         DEALLOCATE (error)
      ELSE
         CALL test_failed(error, "Completed exchange image batches advanced again")
         RETURN
      END IF
      CALL cp2k_exchange_stream_batch_result(stream_batch, calc_cp2k, mol, &
         & cache_cp2k, cache_cp2k_plan, vsh_stream_order, energy_stream_order, error)
      IF (ALLOCATED(error)) RETURN
      CALL cp2k_exchange_stream_batch_result(stream_batch, calc_cp2k, mol, &
         & cache_cp2k, cache_cp2k_plan, vsh_stream, energy_stream, error)
      IF (ALLOCATED(error)) THEN
         IF (TRIM(error%message) /= &
            & "g-xTB exchange image-batch result was already pulled") THEN
            DEALLOCATE (error)
            CALL test_failed(error, "Unexpected duplicate batch-result diagnostic")
            RETURN
         END IF
         DEALLOCATE (error)
      ELSE
         CALL test_failed(error, "Exchange image-batch result was pulled twice")
         RETURN
      END IF
      CALL cp2k_exchange_stream_end(stream_batch, error)
      IF (ALLOCATED(error)) RETURN
      err = MAX(err, MAXVAL(ABS(fock_batch - batch_sentinel - fock_order)), &
         & MAXVAL(ABS(vsh_stream_order - vsh_order)), &
         & ABS(energy_stream_order - energy_order))

      ! A completed batched transaction may not be closed before its scalar result
      ! was consumed.  The failed close must clear the stream for immediate reuse.
      fock_batch = (0.0_wp, 0.0_wp)
      CALL cp2k_exchange_stream_begin(stream_batch, calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, &
         & wfn, nmesh, kfrac_order, weights_order, error, &
         & mode=cp2k_exchange_stream_batched, batch_size=ncell)
      IF (ALLOCATED(error)) RETURN
      DO j = 1, ncell
         ik = korder(j)
         CALL cp2k_exchange_stream_push(stream_batch, ik, &
            & density_order(:, :, :, ik), overlap_order(:, :, ik), error)
         IF (ALLOCATED(error)) RETURN
      END DO
      CALL cp2k_exchange_stream_batch_apply(stream_batch, calc_cp2k, mol, &
         & cache_cp2k, cache_cp2k_plan, error)
      IF (ALLOCATED(error)) RETURN
      DO j = 1, ncell
         ik = korder(j)
         CALL cp2k_exchange_stream_batch_pull(stream_batch, ik, &
            & overlap_order(:, :, ik), fock_batch(:, :, :, ik), error)
         IF (ALLOCATED(error)) RETURN
      END DO
      CALL cp2k_exchange_stream_batch_advance(stream_batch, batches_done, error)
      IF (ALLOCATED(error)) RETURN
      IF (.NOT. batches_done) err = MAX(err, 1.0_wp)
      CALL cp2k_exchange_stream_end(stream_batch, error)
      IF (ALLOCATED(error)) THEN
         IF (TRIM(error%message) /= &
            & "g-xTB exchange stream ended with missing batch result") THEN
            DEALLOCATE (error)
            CALL test_failed(error, "Unexpected missing batch-result close diagnostic")
            RETURN
         END IF
         DEALLOCATE (error)
      ELSE
         CALL test_failed(error, "Exchange image batches ended without a result pull")
         RETURN
      END IF

      ! Batch sizes equal to and greater than K both allocate exactly K image
      ! slots, report complete-image storage, and reproduce the additive oracle.
      DO ibatch_case = 1, 2
         IF (ibatch_case == 1) THEN
            requested_batch_size = ncell
         ELSE
            requested_batch_size = ncell + 2
         END IF
         fock_batch = batch_sentinel
         CALL cp2k_exchange_stream_begin(stream_batch, calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, &
            & wfn, nmesh, kfrac_order, weights_order, error, &
            & mode=cp2k_exchange_stream_batched, &
            & batch_size=requested_batch_size)
         IF (ALLOCATED(error)) RETURN
         IF (.NOT. cp2k_exchange_stream_has_full_mesh_storage(stream_batch)) THEN
            err = MAX(err, 1.0_wp)
         END IF
         IF (ALLOCATED(cache_cp2k_plan%bvk_phase_forward) &
            & .OR. ALLOCATED(cache_cp2k_plan%bvk_phase_inverse)) err = MAX(err, 1.0_wp)
         peak_expected = 4_int64*INT(nao, int64)*INT(nao, int64) &
            & *INT(ncell, int64) + 3_int64*INT(nao, int64) &
            & + 2_int64*INT(ncell, int64) &
            & + 3_int64*INT(MAXVAL(nmesh), int64)
         peak_expected = peak_expected + MAX( &
            & 2_int64*INT(nao, int64)**2, &
            & 5_int64*INT(nao, int64)**2 + 6_int64*INT(nao, int64), &
            & 5_int64*INT(nao, int64)**2 + INT(ncell, int64))
         IF (cp2k_exchange_stream_peak_complex_elements(stream_batch) &
            & /= peak_expected) err = MAX(err, 1.0_wp)
         CALL cp2k_exchange_stream_batch_bounds(stream_batch, first_image, &
            & last_image, error)
         IF (ALLOCATED(error)) RETURN
         IF (first_image /= 1 .OR. last_image /= ncell) err = MAX(err, 1.0_wp)
         DO j = 1, ncell
            ik = korder(j)
            CALL cp2k_exchange_stream_push(stream_batch, ik, &
               & density_order(:, :, :, ik), overlap_order(:, :, ik), error)
            IF (ALLOCATED(error)) RETURN
         END DO
         CALL cp2k_exchange_stream_batch_apply(stream_batch, calc_cp2k, mol, &
            & cache_cp2k, cache_cp2k_plan, error)
         IF (ALLOCATED(error)) RETURN
         DO j = 1, ncell
            ik = korder(j)
            CALL cp2k_exchange_stream_batch_pull(stream_batch, ik, &
               & overlap_order(:, :, ik), fock_batch(:, :, :, ik), error)
            IF (ALLOCATED(error)) RETURN
         END DO
         CALL cp2k_exchange_stream_batch_advance(stream_batch, batches_done, error)
         IF (ALLOCATED(error)) RETURN
         IF (.NOT. batches_done) err = MAX(err, 1.0_wp)
         CALL cp2k_exchange_stream_batch_result(stream_batch, calc_cp2k, mol, &
            & cache_cp2k, cache_cp2k_plan, vsh_stream_order, energy_stream_order, error)
         IF (ALLOCATED(error)) RETURN
         CALL cp2k_exchange_stream_end(stream_batch, error)
         IF (ALLOCATED(error)) RETURN
         err = MAX(err, &
            & MAXVAL(ABS(fock_batch - batch_sentinel - fock_order)), &
            & MAXVAL(ABS(vsh_stream_order - vsh_order)), &
            & ABS(energy_stream_order - energy_order))
      END DO

      ! Independent bounded streams may own disjoint image ranges.  Each range
      ! below is completed with exactly one K-block submission, as intended for
      ! MPI image-range distribution.  The externally summed partial energy,
      ! shell potential, and additive Fock response must reproduce the unchanged
      ! full-mesh oracle without storing a complete AO image mesh on either tile.
      fock_batch = batch_sentinel
      vsh_stream = 0.0_wp
      energy_tile_sum = 0.0_wp
      bvk_builds = cache_cp2k_plan%bvk_plan_builds
      DO ibatch_case = 1, 2
         IF (ibatch_case == 1) THEN
            first_image = 1
            requested_batch_size = 1
         ELSE
            first_image = 2
            requested_batch_size = ncell - 1
         END IF
         CALL cp2k_exchange_stream_begin(stream_batch, calc_cp2k, mol, &
            & cache_cp2k, cache_cp2k_plan, wfn, nmesh, kfrac_order, weights_order, error, &
            & mode=cp2k_exchange_stream_batched, &
            & batch_size=requested_batch_size + 2, image_first=first_image, &
            & image_count=requested_batch_size)
         IF (ALLOCATED(error)) RETURN
         IF (cp2k_exchange_stream_has_full_mesh_storage(stream_batch)) THEN
            err = MAX(err, 1.0_wp)
         END IF
         CALL cp2k_exchange_stream_batch_bounds(stream_batch, i, j, error)
         IF (ALLOCATED(error)) RETURN
         IF (i /= first_image &
            & .OR. j /= first_image + requested_batch_size - 1) THEN
            err = MAX(err, 1.0_wp)
         END IF
         DO j = 1, ncell
            ik = korder(j)
            CALL cp2k_exchange_stream_push(stream_batch, ik, &
               & density_order(:, :, :, ik), overlap_order(:, :, ik), error)
            IF (ALLOCATED(error)) RETURN
         END DO
         CALL cp2k_exchange_stream_batch_apply(stream_batch, calc_cp2k, mol, &
            & cache_cp2k, cache_cp2k_plan, error)
         IF (ALLOCATED(error)) RETURN
         DO j = 1, ncell
            ik = korder(j)
            CALL cp2k_exchange_stream_batch_pull(stream_batch, ik, &
               & overlap_order(:, :, ik), fock_batch(:, :, :, ik), error)
            IF (ALLOCATED(error)) RETURN
         END DO
         CALL cp2k_exchange_stream_batch_advance(stream_batch, batches_done, error)
         IF (ALLOCATED(error)) RETURN
         IF (.NOT. batches_done) err = MAX(err, 1.0_wp)
         CALL cp2k_exchange_stream_batch_result(stream_batch, calc_cp2k, mol, &
            & cache_cp2k, cache_cp2k_plan, vsh_stream_order, energy_tile, error)
         IF (ALLOCATED(error)) RETURN
         vsh_stream = vsh_stream + vsh_stream_order
         energy_tile_sum = energy_tile_sum + energy_tile
         CALL cp2k_exchange_stream_end(stream_batch, error)
         IF (ALLOCATED(error)) RETURN
      END DO
      IF (cache_cp2k_plan%bvk_plan_builds /= bvk_builds) err = MAX(err, 1.0_wp)
      err = MAX(err, MAXVAL(ABS(fock_batch - batch_sentinel - fock_order)), &
         & MAXVAL(ABS(vsh_stream - vsh_order)), ABS(energy_tile_sum - energy_order))

      ! The cached reverse path must reproduce the independently evaluated
      ! overlap response, atomic forces, and homogeneous-strain derivative.
      ALLOCATE (overlap_grad_cp2k(nao, nao, ncell), &
         & gradient_cp2k(3, mol%nat), sigma_cp2k(3, 3))
      CALL cp2k_exchange_kmesh_gradient(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, nmesh, &
         & kfrac, weights, density_mesh, overlap_mesh, overlap_grad_cp2k, &
         & gradient_cp2k, sigma_cp2k, error)
      IF (ALLOCATED(error)) RETURN
      err = MAX(err, MAXVAL(ABS(overlap_grad_cp2k - overlap_grad_mesh)), &
         & MAXVAL(ABS(gradient_cp2k - gradient_mesh)), &
         & MAXVAL(ABS(sigma_cp2k - sigma_mesh)))

      ! Exact composition test: both compact transform backends must retain the
      ! bounded image protocol and reproduce the dense forward and reverse maps.
      ! The 3x1x1 mesh also exercises the mixed-radix prime-axis leaf.
      CALL check_compact_stream_composition(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn, &
         & nmesh, kfrac, weights, density_mesh, overlap_mesh, fock_cp2k, &
         & vsh_cp2k, energy_mesh, overlap_grad_cp2k, gradient_cp2k, sigma_cp2k, &
         & cp2k_exchange_transform_separable_dft, stream_compact_separable_err, &
         & stream_compact_separable_peak, error)
      IF (ALLOCATED(error)) RETURN
      CALL check_compact_stream_composition(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn, &
         & nmesh, kfrac, weights, density_mesh, overlap_mesh, fock_cp2k, &
         & vsh_cp2k, energy_mesh, overlap_grad_cp2k, gradient_cp2k, sigma_cp2k, &
         & cp2k_exchange_transform_mixed_radix_fft, stream_compact_fft_err, &
         & stream_compact_fft_peak, error)
      IF (ALLOCATED(error)) RETURN
      err = MAX(err, stream_compact_separable_err, stream_compact_fft_err)
      separable_err = MAX(separable_err, stream_compact_separable_err)
      fft_err = MAX(fft_err, stream_compact_fft_err)
      PRINT *, "STREAMED_SEPARABLE_COMPOSITION_MAX_RESID", &
         & stream_compact_separable_err, stream_compact_separable_peak
      PRINT *, "STREAMED_MIXED_RADIX_COMPOSITION_MAX_RESID", &
         & stream_compact_fft_err, stream_compact_fft_peak

      ! The production tile policy must preserve all dense forward and reverse
      ! observables, not just the tile-size arithmetic or total energy.
      CALL check_compact_stream_composition(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn, &
         & nmesh, kfrac, weights, density_mesh, overlap_mesh, fock_cp2k, &
         & vsh_cp2k, energy_mesh, overlap_grad_cp2k, gradient_cp2k, sigma_cp2k, &
         & cp2k_exchange_transform_separable_dft, stream_compact_separable_err, &
         & stream_compact_separable_peak, error, adaptive_tiles=.TRUE.)
      IF (ALLOCATED(error)) RETURN
      CALL check_compact_stream_composition(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn, &
         & nmesh, kfrac, weights, density_mesh, overlap_mesh, fock_cp2k, &
         & vsh_cp2k, energy_mesh, overlap_grad_cp2k, gradient_cp2k, sigma_cp2k, &
         & cp2k_exchange_transform_mixed_radix_fft, stream_compact_fft_err, &
         & stream_compact_fft_peak, error, adaptive_tiles=.TRUE.)
      IF (ALLOCATED(error)) RETURN
      err = MAX(err, stream_compact_separable_err, stream_compact_fft_err)
      separable_err = MAX(separable_err, stream_compact_separable_err)
      fft_err = MAX(fft_err, stream_compact_fft_err)
      PRINT *, "ADAPTIVE_SEPARABLE_COMPOSITION_MAX_RESID", stream_compact_separable_err
      PRINT *, "ADAPTIVE_MIXED_RADIX_COMPOSITION_MAX_RESID", stream_compact_fft_err

      ! Compose the image and AO-pair policies together, including all reverse
      ! responses. Keep the single-image cases above as independent controls.
      CALL check_compact_stream_composition(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn, &
         & nmesh, kfrac, weights, density_mesh, overlap_mesh, fock_cp2k, &
         & vsh_cp2k, energy_mesh, overlap_grad_cp2k, gradient_cp2k, sigma_cp2k, &
         & cp2k_exchange_transform_separable_dft, stream_compact_separable_err, &
         & stream_compact_separable_peak, error, adaptive_tiles=.TRUE., adaptive_images=.TRUE.)
      IF (ALLOCATED(error)) RETURN
      CALL check_compact_stream_composition(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn, &
         & nmesh, kfrac, weights, density_mesh, overlap_mesh, fock_cp2k, &
         & vsh_cp2k, energy_mesh, overlap_grad_cp2k, gradient_cp2k, sigma_cp2k, &
         & cp2k_exchange_transform_mixed_radix_fft, stream_compact_fft_err, &
         & stream_compact_fft_peak, error, adaptive_tiles=.TRUE., adaptive_images=.TRUE.)
      IF (ALLOCATED(error)) RETURN
      err = MAX(err, stream_compact_separable_err, stream_compact_fft_err)
      separable_err = MAX(separable_err, stream_compact_separable_err)
      fft_err = MAX(fft_err, stream_compact_fft_err)
      PRINT *, "ADAPTIVE_IMAGES_SEPARABLE_COMPOSITION_MAX_RESID", stream_compact_separable_err
      PRINT *, "ADAPTIVE_IMAGES_MIXED_RADIX_COMPOSITION_MAX_RESID", stream_compact_fft_err

      ! The reverse-only reduced and bounded transactions must reproduce the
      ! complete-mesh oracle without retaining density, overlap, Fock, or overlap
      ! adjoint meshes.  In particular, no forward Fock pull is part of this
      ! derivative-only protocol.
      ALLOCATE (overlap_grad_stream_reduced(nao, nao, ncell), &
         & overlap_grad_stream_batch(nao, nao, ncell), &
         & overlap_grad_stream_probe(nao, nao), source=(0.0_wp, 0.0_wp))
      ALLOCATE (gradient_stream_reverse(3, mol%nat), &
         & sigma_stream_reverse(3, 3), source=0.0_wp)
      CALL cp2k_exchange_stream_begin(stream_reverse_reduced, calc_cp2k, mol, &
         & cache_cp2k, cache_cp2k_plan, wfn, nmesh, kfrac, weights, error, &
         & mode=cp2k_exchange_stream_reduced)
      IF (ALLOCATED(error)) RETURN
      DO ik = 1, ncell
         CALL cp2k_exchange_stream_push(stream_reverse_reduced, ik, &
            & density_mesh(:, :, :, ik), overlap_mesh(:, :, ik), error)
         IF (ALLOCATED(error)) RETURN
      END DO
      CALL cp2k_exchange_stream_apply(stream_reverse_reduced, calc_cp2k, mol, &
         & cache_cp2k, cache_cp2k_plan, vsh_stream, energy_stream, error)
      IF (ALLOCATED(error)) RETURN
      CALL cp2k_exchange_stream_reverse_begin(stream_reverse_reduced, mol, error)
      IF (ALLOCATED(error)) RETURN
      DO ik = 1, ncell
         CALL cp2k_exchange_stream_reverse_push(stream_reverse_reduced, ik, &
            & density_mesh(:, :, :, ik), overlap_mesh(:, :, ik), error)
         IF (ALLOCATED(error)) RETURN
      END DO
      CALL cp2k_exchange_stream_reverse_batch_apply(stream_reverse_reduced, &
         & calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, error)
      IF (ALLOCATED(error)) RETURN
      DO ik = 1, ncell
         CALL cp2k_exchange_stream_reverse_batch_pull(stream_reverse_reduced, &
            & ik, density_mesh(:, :, :, ik), overlap_mesh(:, :, ik), &
            & overlap_grad_stream_reduced(:, :, ik), error)
         IF (ALLOCATED(error)) RETURN
      END DO
      CALL cp2k_exchange_stream_reverse_batch_advance(stream_reverse_reduced, &
         & batches_done, error)
      IF (ALLOCATED(error)) RETURN
      IF (.NOT. batches_done) err = MAX(err, 1.0_wp)
      CALL cp2k_exchange_stream_reverse_result(stream_reverse_reduced, calc_cp2k, &
         & mol, cache_cp2k, cache_cp2k_plan, gradient_stream_reverse, sigma_stream_reverse, error)
      IF (ALLOCATED(error)) RETURN
      err = MAX(err, &
         & MAXVAL(ABS(overlap_grad_stream_reduced - overlap_grad_cp2k)), &
         & MAXVAL(ABS(gradient_stream_reverse - gradient_cp2k)), &
         & MAXVAL(ABS(sigma_stream_reverse - sigma_cp2k)))
      stream_reverse_err = MAX(stream_reverse_err, &
         & MAXVAL(ABS(overlap_grad_stream_reduced - overlap_grad_cp2k)), &
         & MAXVAL(ABS(gradient_stream_reverse - gradient_cp2k)), &
         & MAXVAL(ABS(sigma_stream_reverse - sigma_cp2k)))
      CALL cp2k_exchange_stream_end(stream_reverse_reduced, error)
      IF (ALLOCATED(error)) RETURN

      CALL cp2k_exchange_stream_begin(stream_reverse_batch, calc_cp2k, mol, &
         & cache_cp2k, cache_cp2k_plan, wfn, nmesh, kfrac, weights, error, &
         & mode=cp2k_exchange_stream_batched, batch_size=2, &
         & reverse_only=.TRUE.)
      IF (ALLOCATED(error)) RETURN
      IF (cp2k_exchange_stream_has_full_mesh_storage(stream_reverse_batch)) THEN
         err = MAX(err, 1.0_wp)
      END IF
      peak_expected = 6_int64*INT(nao, int64)**2*2_int64 &
         & *INT(wfn%nspin, int64) + 5_int64*INT(nao, int64) &
         & *INT(wfn%nspin, int64) + 2_int64*INT(ncell, int64) &
         & + 3_int64*INT(MAXVAL(nmesh), int64) + 10_int64*INT(nao, int64)**2
      IF (cp2k_exchange_stream_peak_complex_elements(stream_reverse_batch) &
         & /= peak_expected) err = MAX(err, 1.0_wp)
      batches_done = .FALSE.
      DO WHILE (.NOT. batches_done)
         DO ik = 1, ncell
            CALL cp2k_exchange_stream_reverse_push(stream_reverse_batch, ik, &
               & density_mesh(:, :, :, ik), overlap_mesh(:, :, ik), error)
            IF (ALLOCATED(error)) RETURN
            IF (ik == 1) THEN
               CALL cp2k_exchange_stream_reverse_push(stream_reverse_batch, ik, &
                  & density_mesh(:, :, :, ik), overlap_mesh(:, :, ik), error)
               IF (ALLOCATED(error)) THEN
                  IF (TRIM(error%message) /= &
                     & "duplicate g-xTB exchange reverse block") THEN
                     DEALLOCATE (error)
                     CALL test_failed(error, &
                        & "Unexpected duplicate reverse-push diagnostic")
                     RETURN
                  END IF
                  DEALLOCATE (error)
               ELSE
                  CALL test_failed(error, "Duplicate reverse push was accepted")
                  RETURN
               END IF
            END IF
         END DO
         CALL cp2k_exchange_stream_reverse_batch_apply(stream_reverse_batch, &
            & calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, error)
         IF (ALLOCATED(error)) RETURN
         CALL cp2k_exchange_stream_reverse_batch_apply(stream_reverse_batch, &
            & calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, error)
         IF (ALLOCATED(error)) THEN
            IF (TRIM(error%message) /= &
               & "g-xTB exchange stream reverse batch was already applied") THEN
               DEALLOCATE (error)
               CALL test_failed(error, &
                  & "Unexpected duplicate reverse-apply diagnostic")
               RETURN
            END IF
            DEALLOCATE (error)
         ELSE
            CALL test_failed(error, "Duplicate reverse-batch apply was accepted")
            RETURN
         END IF
         DO ik = 1, ncell
            CALL cp2k_exchange_stream_reverse_batch_pull(stream_reverse_batch, &
               & ik, density_mesh(:, :, :, ik), overlap_mesh(:, :, ik), &
               & overlap_grad_stream_batch(:, :, ik), error)
            IF (ALLOCATED(error)) RETURN
            IF (ik == 1) THEN
               CALL cp2k_exchange_stream_reverse_batch_pull( &
                  & stream_reverse_batch, ik, density_mesh(:, :, :, ik), &
                  & overlap_mesh(:, :, ik), &
                  & overlap_grad_stream_probe, error)
               IF (ALLOCATED(error)) THEN
                  IF (TRIM(error%message) /= &
                     & "duplicate g-xTB exchange reverse-batch pull") THEN
                     DEALLOCATE (error)
                     CALL test_failed(error, &
                        & "Unexpected duplicate reverse-pull diagnostic")
                     RETURN
                  END IF
                  DEALLOCATE (error)
               ELSE
                  CALL test_failed(error, &
                     & "Duplicate reverse-batch pull was accepted")
                  RETURN
               END IF
            END IF
         END DO
         CALL cp2k_exchange_stream_reverse_batch_advance(stream_reverse_batch, &
            & batches_done, error)
         IF (ALLOCATED(error)) RETURN
      END DO
      CALL cp2k_exchange_stream_reverse_result(stream_reverse_batch, calc_cp2k, &
         & mol, cache_cp2k, cache_cp2k_plan, gradient_stream_reverse, sigma_stream_reverse, error)
      IF (ALLOCATED(error)) RETURN
      CALL cp2k_exchange_stream_reverse_result(stream_reverse_batch, calc_cp2k, &
         & mol, cache_cp2k, cache_cp2k_plan, gradient_stream_probe, sigma_stream_probe, error)
      IF (ALLOCATED(error)) THEN
         IF (TRIM(error%message) /= &
            & "g-xTB exchange reverse result was already pulled") THEN
            DEALLOCATE (error)
            CALL test_failed(error, &
               & "Unexpected duplicate reverse-result diagnostic")
            RETURN
         END IF
         DEALLOCATE (error)
      ELSE
         CALL test_failed(error, "Duplicate reverse result was accepted")
         RETURN
      END IF
      err = MAX(err, &
         & MAXVAL(ABS(overlap_grad_stream_batch - overlap_grad_cp2k)), &
         & MAXVAL(ABS(gradient_stream_reverse - gradient_cp2k)), &
         & MAXVAL(ABS(sigma_stream_reverse - sigma_cp2k)))
      stream_reverse_err = MAX(stream_reverse_err, &
         & MAXVAL(ABS(overlap_grad_stream_batch - overlap_grad_cp2k)), &
         & MAXVAL(ABS(gradient_stream_reverse - gradient_cp2k)), &
         & MAXVAL(ABS(sigma_stream_reverse - sigma_cp2k)))
      PRINT *, "STREAM_REVERSE_RKS_3X1X1_RESID", &
         & MAXVAL(ABS(overlap_grad_stream_batch - overlap_grad_cp2k)), &
         & MAXVAL(ABS(gradient_stream_reverse - gradient_cp2k)), &
         & MAXVAL(ABS(sigma_stream_reverse - sigma_cp2k)), &
         & cp2k_exchange_stream_peak_complex_elements(stream_reverse_batch)
      CALL cp2k_exchange_stream_end(stream_reverse_batch, error)
      IF (ALLOCATED(error)) RETURN

      ! Exercise the separate reverse-stream phase without changing the forward
      ! transaction tests above.  A stale onsite state must be rejected before
      ! the unchanged whole-mesh gradient oracle is entered.
      ALLOCATE (overlap_grad_stream_reverse(nao, nao, ncell))
      CALL cp2k_exchange_stream_begin(stream_reverse, calc_cp2k, mol, &
         & cache_cp2k, cache_cp2k_plan, wfn, nmesh, kfrac, weights, error, &
         & mode=cp2k_exchange_stream_oracle)
      IF (ALLOCATED(error)) RETURN
      IF (.NOT. cp2k_exchange_stream_has_full_mesh_storage(stream_reverse)) THEN
         err = MAX(err, 1.0_wp)
      END IF
      DO ik = 1, ncell
         CALL cp2k_exchange_stream_push(stream_reverse, ik, &
            & density_mesh(:, :, :, ik), overlap_mesh(:, :, ik), error)
         IF (ALLOCATED(error)) RETURN
      END DO
      CALL cp2k_exchange_stream_apply(stream_reverse, calc_cp2k, mol, &
         & cache_cp2k, cache_cp2k_plan, vsh_stream, energy_stream, error)
      IF (ALLOCATED(error)) RETURN
      DO ik = 1, ncell
         CALL cp2k_exchange_stream_get(stream_reverse, ik, &
            & fock_stream_order(:, :, :, ik), error)
         IF (ALLOCATED(error)) RETURN
      END DO
      err = MAX(err, MAXVAL(ABS(fock_stream_order - fock_cp2k)), &
         & MAXVAL(ABS(vsh_stream - vsh_cp2k)), ABS(energy_stream - energy_mesh))

      onsite_probe = ptr_cp2k%g_onsfx(1, 1, 1)
      ptr_cp2k%g_onsfx(1, 1, 1) = onsite_probe + 1.0e-8_wp
      CALL cp2k_exchange_stream_reverse_apply(stream_reverse, calc_cp2k, mol, &
         & cache_cp2k, cache_cp2k_plan, gradient_stream_probe, sigma_stream_probe, error)
      IF (ALLOCATED(error)) THEN
         IF (TRIM(error%message) /= &
            & "g-xTB exchange stream state was invalidated before reverse apply") THEN
            DEALLOCATE (error)
            CALL test_failed(error, "Unexpected stale reverse-stream diagnostic")
            RETURN
         END IF
         DEALLOCATE (error)
      ELSE
         CALL test_failed(error, "Stale exchange reverse stream was applied")
         RETURN
      END IF
      ptr_cp2k%g_onsfx(1, 1, 1) = onsite_probe

      CALL cp2k_exchange_stream_reverse_apply(stream_reverse, calc_cp2k, mol, &
         & cache_cp2k, cache_cp2k_plan, gradient_stream_reverse, sigma_stream_reverse, error)
      IF (ALLOCATED(error)) RETURN
      CALL cp2k_exchange_stream_reverse_apply(stream_reverse, calc_cp2k, mol, &
         & cache_cp2k, cache_cp2k_plan, gradient_stream_probe, sigma_stream_probe, error)
      IF (ALLOCATED(error)) THEN
         IF (TRIM(error%message) /= &
            & "g-xTB exchange stream reverse was already applied") THEN
            DEALLOCATE (error)
            CALL test_failed(error, "Unexpected repeated reverse-apply diagnostic")
            RETURN
         END IF
         DEALLOCATE (error)
      ELSE
         CALL test_failed(error, "Exchange reverse stream was applied twice")
         RETURN
      END IF

      CALL cp2k_exchange_stream_reverse_get(stream_reverse, 1, &
         & overlap_grad_stream_reverse(:, :, 1), error)
      IF (ALLOCATED(error)) RETURN
      CALL cp2k_exchange_stream_reverse_get(stream_reverse, 1, &
         & overlap_grad_stream_probe, error)
      IF (ALLOCATED(error)) THEN
         IF (TRIM(error%message) /= &
            & "duplicate g-xTB exchange stream reverse pull") THEN
            DEALLOCATE (error)
            CALL test_failed(error, "Unexpected duplicate reverse-pull diagnostic")
            RETURN
         END IF
         DEALLOCATE (error)
      ELSE
         CALL test_failed(error, "Exchange overlap adjoint was pulled twice")
         RETURN
      END IF
      DO ik = 2, ncell
         CALL cp2k_exchange_stream_reverse_get(stream_reverse, ik, &
            & overlap_grad_stream_reverse(:, :, ik), error)
         IF (ALLOCATED(error)) RETURN
      END DO
      CALL cp2k_exchange_stream_end(stream_reverse, error)
      IF (ALLOCATED(error)) RETURN
      err = MAX(err, &
         & MAXVAL(ABS(overlap_grad_stream_reverse - overlap_grad_cp2k)), &
         & MAXVAL(ABS(overlap_grad_stream_reverse - overlap_grad_mesh)), &
         & MAXVAL(ABS(gradient_stream_reverse - gradient_cp2k)), &
         & MAXVAL(ABS(gradient_stream_reverse - gradient_mesh)), &
         & MAXVAL(ABS(sigma_stream_reverse - sigma_cp2k)), &
         & MAXVAL(ABS(sigma_stream_reverse - sigma_mesh)))

      ! Once reverse evaluation was requested, closing the stream is legal only
      ! after every overlap-adjoint block was pulled exactly once.
      CALL cp2k_exchange_stream_begin(stream_reverse_missing, calc_cp2k, mol, &
         & cache_cp2k, cache_cp2k_plan, wfn, nmesh, kfrac, weights, error, &
         & mode=cp2k_exchange_stream_oracle)
      IF (ALLOCATED(error)) RETURN
      DO j = 1, ncell
         ik = korder(j)
         CALL cp2k_exchange_stream_push(stream_reverse_missing, ik, &
            & density_mesh(:, :, :, ik), overlap_mesh(:, :, ik), error)
         IF (ALLOCATED(error)) RETURN
      END DO
      CALL cp2k_exchange_stream_apply(stream_reverse_missing, calc_cp2k, mol, &
         & cache_cp2k, cache_cp2k_plan, vsh_stream, energy_stream, error)
      IF (ALLOCATED(error)) RETURN
      DO ik = 1, ncell
         CALL cp2k_exchange_stream_get(stream_reverse_missing, ik, &
            & fock_stream_order(:, :, :, ik), error)
         IF (ALLOCATED(error)) RETURN
      END DO
      CALL cp2k_exchange_stream_reverse_apply(stream_reverse_missing, calc_cp2k, &
         & mol, cache_cp2k, cache_cp2k_plan, gradient_stream_probe, sigma_stream_probe, error)
      IF (ALLOCATED(error)) RETURN
      CALL cp2k_exchange_stream_reverse_get(stream_reverse_missing, 1, &
         & overlap_grad_stream_probe, error)
      IF (ALLOCATED(error)) RETURN
      CALL cp2k_exchange_stream_end(stream_reverse_missing, error)
      IF (ALLOCATED(error)) THEN
         IF (TRIM(error%message) /= &
            & "g-xTB exchange stream ended with missing reverse pulls") THEN
            DEALLOCATE (error)
            CALL test_failed(error, "Unexpected missing reverse-pull diagnostic")
            RETURN
         END IF
         DEALLOCATE (error)
      ELSE
         CALL test_failed(error, "Exchange stream with missing reverse pulls ended")
         RETURN
      END IF

      bvk_builds = cache_cp2k_plan%bvk_plan_builds
      CALL cp2k_exchange_kmesh(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn, nmesh, &
         & kfrac, weights, density_mesh, overlap_mesh, fock_cp2k, &
         & vsh_cp2k, energy_mesh, error)
      IF (ALLOCATED(error)) RETURN
      IF (cache_cp2k_plan%bvk_plan_builds /= bvk_builds) err = MAX(err, 1.0_wp)

      ! Every stored model component and the active representative ordering are
      ! part of the persistent-plan identity.  Mutating any one of them must
      ! reject the cache before a stale transform can be selected.
      SELECT TYPE (fock_cp2k_model => calc_cp2k%exchange)
      TYPE is (exchange_fock)
         CALL check_bvk_model_signature(fock_cp2k_model, cache_cp2k_plan, error)
         IF (ALLOCATED(error)) RETURN
         kernel_order_probe = cache_cp2k_plan%bvk_kernel
         IF (.NOT. cache_cp2k_plan%plan_matches(fock_cp2k_model, mol, &
            & kernel_order_probe, kfrac, weights, require_phases=.FALSE.)) THEN
            CALL test_failed(error, "Fresh BvK representative order does not match")
            RETURN
         END IF
         kernel_order_probe%reps(:, [1, 2]) = &
            & kernel_order_probe%reps(:, [2, 1])
         IF (cache_cp2k_plan%plan_matches(fock_cp2k_model, mol, &
            & kernel_order_probe, kfrac, weights, require_phases=.FALSE.)) THEN
            CALL test_failed(error, "BvK plan accepted a changed representative order")
            RETURN
         END IF
      class default
         CALL test_failed(error, "Unexpected exchange implementation")
         RETURN
      END SELECT

      CALL cp2k_exchange_kmesh(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn, nmesh, &
         & kfrac_twist, weights, density_mesh, overlap_mesh, fock_cp2k, &
         & vsh_cp2k, energy_mesh, error)
      IF (ALLOCATED(error)) RETURN
      IF (cache_cp2k_plan%bvk_plan_builds /= bvk_builds + 1) err = MAX(err, 1.0_wp)
      err = MAX(err, MAXVAL(ABS(fock_cp2k - fock_twist)), &
         & MAXVAL(ABS(vsh_cp2k - vsh_twist)), ABS(energy_mesh - energy_twist))

      bvk_builds = cache_cp2k_plan%bvk_plan_builds
      CALL cp2k_exchange_kmesh(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn, nmesh, &
         & kfrac_order, weights_order, density_order, overlap_order, fock_cp2k, &
         & vsh_cp2k, energy_mesh, error)
      IF (ALLOCATED(error)) RETURN
      IF (cache_cp2k_plan%bvk_plan_builds /= bvk_builds + 1) err = MAX(err, 1.0_wp)
      IF (ANY(cache_cp2k_plan%bvk_input_to_grid /= korder)) err = MAX(err, 1.0_wp)
      IF (ANY(cache_cp2k_plan%bvk_grid_to_input /= [2, 3, 1])) THEN
         err = MAX(err, 1.0_wp)
      END IF
      err = MAX(err, MAXVAL(ABS(fock_cp2k - fock_order)), &
         & MAXVAL(ABS(vsh_cp2k - vsh_order)), ABS(energy_mesh - energy_order))
      CALL cp2k_exchange_kmesh(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn, nmesh, &
         & kfrac_order, weights_order, density_order, overlap_order, fock_cp2k, &
         & vsh_cp2k, energy_mesh, error, &
         & transform_backend=cp2k_exchange_transform_separable_dft)
      IF (ALLOCATED(error)) RETURN
      separable_err = MAX(separable_err, &
         & MAXVAL(ABS(fock_cp2k - fock_order)), &
         & MAXVAL(ABS(vsh_cp2k - vsh_order)), ABS(energy_mesh - energy_order))
      CALL cp2k_exchange_kmesh(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn, nmesh, &
         & kfrac_order, weights_order, density_order, overlap_order, fock_cp2k, &
         & vsh_cp2k, energy_mesh, error, &
         & transform_backend=cp2k_exchange_transform_mixed_radix_fft)
      IF (ALLOCATED(error)) RETURN
      fft_err = MAX(fft_err, MAXVAL(ABS(fock_cp2k - fock_order)), &
         & MAXVAL(ABS(vsh_cp2k - vsh_order)), ABS(energy_mesh - energy_order))
      PRINT *, "MIXED_RADIX_FFT_RKS_FORWARD_RESID", &
         & MAXVAL(ABS(fock_cp2k - fock_order)), MAXVAL(ABS(vsh_cp2k - vsh_order)), &
         & ABS(energy_mesh - energy_order)
      CALL check_exchange_transform_selectors(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn, &
         & nmesh, kfrac_order, weights_order, density_order, overlap_order, &
         & error, selector_separable_err, selector_fft_err, rks_fft_forward, &
         & rks_fft_reverse)
      IF (ALLOCATED(error)) RETURN
      separable_err = MAX(separable_err, selector_separable_err)
      fft_err = MAX(fft_err, selector_fft_err)
      PRINT *, "TRANSFORM_SELECTOR_RKS_SEPARABLE_MAX_RESID", &
         & selector_separable_err
      PRINT *, "MIXED_RADIX_FFT_RKS_REVERSE_RESID", rks_fft_reverse
      err = MAX(err, separable_err)

      ! Every geometry, cell, mesh, and model dependency must trigger an actual
      ! rebuild.  The AO inputs can remain fixed here because this test probes
      ! plan ownership and invalidation rather than a self-consistent state.
      bvk_builds = cache_cp2k_plan%bvk_plan_builds
      mol_probe = mol
      mol_probe%xyz(1, 1) = mol_probe%xyz(1, 1) + 1.0e-4_wp
      CALL cp2k_exchange_kmesh(calc_cp2k, mol_probe, cache_cp2k, cache_cp2k_plan, wfn, nmesh, &
         & kfrac_order, weights_order, density_order, overlap_order, fock_cp2k, &
         & vsh_cp2k, energy_mesh, error)
      IF (ALLOCATED(error)) RETURN
      IF (cache_cp2k_plan%bvk_plan_builds /= bvk_builds + 1) err = MAX(err, 1.0_wp)
      bvk_builds = cache_cp2k_plan%bvk_plan_builds
      mol_probe = mol
      mol_probe%lattice(1, 1) = mol_probe%lattice(1, 1) + 1.0e-4_wp
      CALL cp2k_exchange_kmesh(calc_cp2k, mol_probe, cache_cp2k, cache_cp2k_plan, wfn, nmesh, &
         & kfrac_order, weights_order, density_order, overlap_order, fock_cp2k, &
         & vsh_cp2k, energy_mesh, error)
      IF (ALLOCATED(error)) RETURN
      IF (cache_cp2k_plan%bvk_plan_builds /= bvk_builds + 1) err = MAX(err, 1.0_wp)

      bvk_builds = cache_cp2k_plan%bvk_plan_builds
      kfrac_bad = 0.0_wp
      DO ik = 1, ncell
         kfrac_bad(2, ik) = REAL(ik - 1, wp)/REAL(ncell, wp)
      END DO
      CALL cp2k_exchange_kmesh(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn, [1, 3, 1], &
         & kfrac_bad, weights, density_mesh, overlap_mesh, fock_cp2k, &
         & vsh_cp2k, energy_mesh, error)
      IF (ALLOCATED(error)) RETURN
      IF (cache_cp2k_plan%bvk_plan_builds /= bvk_builds + 1 &
         & .OR. ANY(cache_cp2k_plan%bvk_kernel%nmesh /= [1, 3, 1])) THEN
         err = MAX(err, 1.0_wp)
      END IF

      bvk_builds = cache_cp2k_plan%bvk_plan_builds
      SELECT TYPE (fock_cp2k_model => calc_cp2k%exchange)
      TYPE is (exchange_fock)
         fock_cp2k_model%corr_exp = fock_cp2k_model%corr_exp + 1.0e-6_wp
      END SELECT
      CALL cp2k_exchange_kmesh(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn, nmesh, &
         & kfrac_order, weights_order, density_order, overlap_order, fock_cp2k, &
         & vsh_cp2k, energy_mesh, error)
      IF (ALLOCATED(error)) RETURN
      IF (cache_cp2k_plan%bvk_plan_builds /= bvk_builds + 1) err = MAX(err, 1.0_wp)
      SELECT TYPE (fock_cp2k_model => calc_cp2k%exchange)
      TYPE is (exchange_fock)
         fock_cp2k_model%corr_exp = fock_cp2k_model%corr_exp - 1.0e-6_wp
      END SELECT

      ! Reject incomplete and duplicate grids.  A boundary rejection leaves a
      ! valid plan untouched; failed internal construction invalidates it, and
      ! the next valid call must recover without incrementing on the failures.
      bvk_builds = cache_cp2k_plan%bvk_plan_builds
      CALL cp2k_exchange_kmesh(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn, [-1, 1, 1], &
         & kfrac_order, weights_order, density_order, overlap_order, fock_cp2k, &
         & vsh_cp2k, energy_mesh, error, &
         & transform_backend=cp2k_exchange_transform_mixed_radix_fft)
      IF (ALLOCATED(error)) THEN
         IF (TRIM(error%message) /= &
            & "g-xTB exchange requires a complete regular k mesh") THEN
            DEALLOCATE (error)
            CALL test_failed(error, "Unexpected negative FFT-mesh diagnostic")
            RETURN
         END IF
         DEALLOCATE (error)
      ELSE
         CALL test_failed(error, "Negative mixed-radix FFT mesh was accepted")
         RETURN
      END IF
      CALL cp2k_exchange_kmesh(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn, nmesh_overflow, &
         & kfrac_order, weights_order, density_order, overlap_order, fock_cp2k, &
         & vsh_cp2k, energy_mesh, error, &
         & transform_backend=cp2k_exchange_transform_mixed_radix_fft)
      IF (ALLOCATED(error)) THEN
         IF (TRIM(error%message) /= &
            & "g-xTB exchange requires a complete regular k mesh") THEN
            DEALLOCATE (error)
            CALL test_failed(error, "Unexpected overflowing FFT-mesh diagnostic")
            RETURN
         END IF
         DEALLOCATE (error)
      ELSE
         CALL test_failed(error, "Overflowing mixed-radix FFT mesh was accepted")
         RETURN
      END IF
      IF (cache_cp2k_plan%bvk_plan_builds /= bvk_builds &
         & .OR. .NOT. cache_cp2k_plan%bvk_plan_valid) err = MAX(err, 1.0_wp)

      CALL cp2k_exchange_kmesh(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn, [4, 1, 1], &
         & kfrac_order, weights_order, density_order, overlap_order, fock_cp2k, &
         & vsh_cp2k, energy_mesh, error)
      IF (ALLOCATED(error)) THEN
         IF (TRIM(error%message) /= &
            & "g-xTB exchange requires a complete regular k mesh") THEN
            DEALLOCATE (error)
            CALL test_failed(error, "Unexpected incomplete-grid diagnostic")
            RETURN
         END IF
         DEALLOCATE (error)
      ELSE
         CALL test_failed(error, "Incomplete whole-mesh exchange grid was accepted")
         RETURN
      END IF
      IF (cache_cp2k_plan%bvk_plan_builds /= bvk_builds &
         & .OR. .NOT. cache_cp2k_plan%bvk_plan_valid) err = MAX(err, 1.0_wp)

      kfrac_bad = kfrac_order
      kfrac_bad(:, 2) = kfrac_bad(:, 1)
      CALL cp2k_exchange_kmesh(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn, nmesh, &
         & kfrac_bad, weights_order, density_order, overlap_order, fock_cp2k, &
         & vsh_cp2k, energy_mesh, error)
      IF (ALLOCATED(error)) THEN
         IF (TRIM(error%message) /= &
            & "g-xTB k-point mesh contains a duplicate grid point") THEN
            DEALLOCATE (error)
            CALL test_failed(error, "Unexpected duplicate-grid diagnostic")
            RETURN
         END IF
         DEALLOCATE (error)
      ELSE
         CALL test_failed(error, "Duplicate whole-mesh exchange grid point was accepted")
         RETURN
      END IF
      IF (cache_cp2k_plan%bvk_plan_builds /= bvk_builds &
         & .OR. cache_cp2k_plan%bvk_plan_valid) err = MAX(err, 1.0_wp)
      CALL cp2k_exchange_kmesh(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn, nmesh, &
         & kfrac_order, weights_order, density_order, overlap_order, fock_cp2k, &
         & vsh_cp2k, energy_mesh, error)
      IF (ALLOCATED(error)) RETURN
      IF (cache_cp2k_plan%bvk_plan_builds /= bvk_builds + 1 &
         & .OR. .NOT. cache_cp2k_plan%bvk_plan_valid) err = MAX(err, 1.0_wp)

      ! A geometrically near-regular mesh is accepted by the dense diagnostic
      ! tolerance but must not be silently snapped onto the factorized grid.
      kfrac_bad = kfrac_order
      kfrac_bad(1, 1) = kfrac_bad(1, 1) + 1.0e-12_wp
      CALL cp2k_exchange_kmesh(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn, nmesh, &
         & kfrac_bad, weights_order, density_order, overlap_order, fock_cp2k, &
         & vsh_cp2k, energy_mesh, error, &
         & transform_backend=cp2k_exchange_transform_separable_dft)
      IF (ALLOCATED(error)) THEN
         IF (TRIM(error%message) /= &
            & "separable-DFT g-xTB exchange requires a machine-regular grid") THEN
            DEALLOCATE (error)
            CALL test_failed(error, &
               & "Unexpected non-machine-regular-grid diagnostic")
            RETURN
         END IF
         DEALLOCATE (error)
      ELSE
         CALL test_failed(error, &
            & "Non-machine-regular separable-DFT grid was accepted")
         RETURN
      END IF
      CALL cp2k_exchange_kmesh(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn, nmesh, &
         & kfrac_bad, weights_order, density_order, overlap_order, fock_cp2k, &
         & vsh_cp2k, energy_mesh, error, &
         & transform_backend=cp2k_exchange_transform_mixed_radix_fft)
      IF (ALLOCATED(error)) THEN
         IF (TRIM(error%message) /= &
            & "mixed-radix FFT g-xTB exchange requires a machine-regular grid") THEN
            DEALLOCATE (error)
            CALL test_failed(error, &
               & "Unexpected non-machine-regular FFT-grid diagnostic")
            RETURN
         END IF
         DEALLOCATE (error)
      ELSE
         CALL test_failed(error, &
            & "Non-machine-regular mixed-radix FFT grid was accepted")
         RETURN
      END IF

      ! Exercise the production O(Nk) grid proof beyond the small-mesh dense
      ! orthogonality oracle.  The cached 9x9x1 forward and reverse paths are
      ! compared with the unchanged dense-transform implementation.
      ALLOCATE (kfrac_large(3, PRODUCT(nmesh_large)), &
         & weights_large(PRODUCT(nmesh_large)))
      twist_large = [0.037_wp, 0.021_wp, 0.013_wp]
      weights_large = 1.0_wp/REAL(PRODUCT(nmesh_large), wp)
      DO ik = 1, PRODUCT(nmesh_large)
         ! A coprime affine permutation makes the external order unrelated to
         ! the canonical 9x9x1 grid order while retaining one common twist.
         j = MODULO(17*(ik - 1) + 5, PRODUCT(nmesh_large)) + 1
         kfrac_large(1, ik) = twist_large(1) &
            & + REAL(MODULO(j - 1, nmesh_large(1)), wp) &
            & /REAL(nmesh_large(1), wp)
         kfrac_large(2, ik) = twist_large(2) &
            & + REAL((j - 1)/nmesh_large(1), wp) &
            & /REAL(nmesh_large(2), wp)
         kfrac_large(3, ik) = twist_large(3)
      END DO
      ALLOCATE (density_large(nao, nao, 2, PRODUCT(nmesh_large)), &
         & overlap_large(nao, nao, PRODUCT(nmesh_large)), &
         & fock_large(nao, nao, 2, PRODUCT(nmesh_large)), &
         & fock_large_oracle(nao, nao, 2, PRODUCT(nmesh_large)), &
         & vsh_large(nsh), vsh_large_oracle(nsh))
      DO ik = 1, PRODUCT(nmesh_large)
         theta_x = 2.0_wp*pi*kfrac_large(1, ik)
         theta_y = 2.0_wp*pi*kfrac_large(2, ik)
         overlap_large(:, :, ik) = (0.0_wp, 0.0_wp)
         density_large(:, :, :, ik) = (0.0_wp, 0.0_wp)
         DO ii = 1, nao
            overlap_large(ii, ii, ik) = CMPLX( &
               & REAL(overlap_mesh(ii, ii, 1), wp) &
               & + 0.003_wp*COS(theta_x + REAL(ii, wp)), 0.0_wp, wp)
            density_large(ii, ii, 1, ik) = CMPLX( &
               & 0.56_wp*REAL(density_mesh(ii, ii, 1, 1), wp) &
               & + 0.002_wp*COS(theta_x + REAL(2*ii, wp)), 0.0_wp, wp)
            density_large(ii, ii, 2, ik) = CMPLX( &
               & 0.44_wp*REAL(density_mesh(ii, ii, 1, 1), wp) &
               & + 0.0015_wp*SIN(theta_y + REAL(3*ii, wp)), 0.0_wp, wp)
            DO jj = 1, ii - 1
               overlap_large(jj, ii, ik) = 0.5_wp*( &
                  & overlap_mesh(jj, ii, 1) + CONJG(overlap_mesh(ii, jj, 1))) &
                  & + CMPLX(0.001_wp*COS(theta_x + REAL(ii + jj, wp)), &
                  & 0.001_wp*SIN(theta_y + REAL(2*ii + jj, wp)), wp)
               overlap_large(ii, jj, ik) = &
                  & CONJG(overlap_large(jj, ii, ik))
               density_large(jj, ii, 1, ik) = 0.56_wp*0.5_wp*( &
                  & density_mesh(jj, ii, 1, 1) &
                  & + CONJG(density_mesh(ii, jj, 1, 1))) &
                  & + CMPLX(0.0012_wp*COS(theta_x + REAL(ii + 2*jj, wp)), &
                  & 0.0009_wp*SIN(theta_y + REAL(3*ii + jj, wp)), wp)
               density_large(ii, jj, 1, ik) = &
                  & CONJG(density_large(jj, ii, 1, ik))
               density_large(jj, ii, 2, ik) = 0.44_wp*0.5_wp*( &
                  & density_mesh(jj, ii, 1, 1) &
                  & + CONJG(density_mesh(ii, jj, 1, 1))) &
                  & + CMPLX(0.0008_wp*SIN(theta_x + REAL(ii + jj, wp)), &
                  & 0.0011_wp*COS(theta_y + REAL(ii + 3*jj, wp)), wp)
               density_large(ii, jj, 2, ik) = &
                  & CONJG(density_large(jj, ii, 2, ik))
            END DO
         END DO
      END DO
      IF (MAXVAL(ABS(overlap_large(:, :, 1) &
         & - CONJG(TRANSPOSE(overlap_large(:, :, 1))))) > 1.0e-14_wp &
         & .OR. MAXVAL(ABS(density_large(:, :, 1, 1) &
         & - CONJG(TRANSPOSE(density_large(:, :, 1, 1))))) > 1.0e-14_wp &
         & .OR. MAXVAL(ABS(density_large(:, :, 2, 1) &
         & - CONJG(TRANSPOSE(density_large(:, :, 2, 1))))) > 1.0e-14_wp) THEN
         err = MAX(err, 1.0_wp)
      END IF
      IF (MAXVAL(ABS(overlap_large(:, :, 2) - overlap_large(:, :, 1))) &
         & < 1.0e-6_wp .OR. MAXVAL(ABS(density_large(:, :, :, 2) &
         & - density_large(:, :, :, 1))) < 1.0e-6_wp) err = MAX(err, 1.0_wp)
      ALLOCATE (density_image_large(nao, nao), overlap_image_large(nao, nao), &
         & source=(0.0_wp, 0.0_wp))
      DO ik = 1, PRODUCT(nmesh_large)
         angle = 2.0_wp*pi*kfrac_large(1, ik)
         phase = weights_large(ik)*EXP(CMPLX(0.0_wp, -angle, wp))
         density_image_large = density_image_large &
            & + phase*density_large(:, :, 1, ik)
         overlap_image_large = overlap_image_large + phase*overlap_large(:, :, ik)
      END DO
      IF (MAXVAL(ABS(density_image_large)) < 1.0e-6_wp &
         & .OR. MAXVAL(ABS(overlap_image_large)) < 1.0e-6_wp) THEN
         err = MAX(err, 1.0_wp)
      END IF
      CALL new_wavefunction(wfn_large, mol%nat, nsh, nao, 2, 0.0_wp, .FALSE.)
      wfn_large%qsh(:, 1) = wfn%qsh(:, 1)
      bvk_builds = cache_cp2k_plan%bvk_plan_builds
      CALL cp2k_exchange_kmesh(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn_large, nmesh_large, &
         & kfrac_large, weights_large, density_large, overlap_large, fock_large, &
         & vsh_large, energy_large, error)
      IF (ALLOCATED(error)) RETURN
      IF (cache_cp2k_plan%bvk_plan_builds /= bvk_builds + 1) err = MAX(err, 1.0_wp)
      kernel_large = cache_cp2k_plan%bvk_kernel
      SELECT TYPE (fock => exchange)
      TYPE is (exchange_fock)
         CALL get_KFock_kmesh(fock, mol, ptr, oracle_plan, kernel_large, kfrac_large, &
            & weights_large, density_large, overlap_large, fock_large_oracle, &
            & vsh_large_oracle, energy_large_oracle)
      END SELECT
      err = MAX(err, MAXVAL(ABS(fock_large - fock_large_oracle)), &
         & MAXVAL(ABS(vsh_large - vsh_large_oracle)), &
         & ABS(energy_large - energy_large_oracle))
      ! The selectable separable direct-DFT backend must reproduce the unchanged
      ! dense phase-table oracle for the full UKS energy, Fock response, and
      ! charge-dependent shell potential.  This 9x9x1 mesh also exercises a
      ! common twist and the noncanonical input order constructed above.
      CALL cp2k_exchange_kmesh(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn_large, &
         & nmesh_large, kfrac_large, weights_large, density_large, &
         & overlap_large, fock_large, vsh_large, energy_large, error, &
         & transform_backend=cp2k_exchange_transform_separable_dft)
      IF (ALLOCATED(error)) RETURN
      err = MAX(err, MAXVAL(ABS(fock_large - fock_large_oracle)), &
         & MAXVAL(ABS(vsh_large - vsh_large_oracle)), &
         & ABS(energy_large - energy_large_oracle))
      separable_err = MAX(separable_err, &
         & MAXVAL(ABS(fock_large - fock_large_oracle)), &
         & MAXVAL(ABS(vsh_large - vsh_large_oracle)), &
         & ABS(energy_large - energy_large_oracle))
      CALL cp2k_exchange_kmesh(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn_large, &
         & nmesh_large, kfrac_large, weights_large, density_large, &
         & overlap_large, fock_large, vsh_large, energy_large, error, &
         & transform_backend=cp2k_exchange_transform_mixed_radix_fft)
      IF (ALLOCATED(error)) RETURN
      fft_err = MAX(fft_err, MAXVAL(ABS(fock_large - fock_large_oracle)), &
         & MAXVAL(ABS(vsh_large - vsh_large_oracle)), &
         & ABS(energy_large - energy_large_oracle))
      PRINT *, "MIXED_RADIX_FFT_UKS_FORWARD_RESID", &
         & MAXVAL(ABS(fock_large - fock_large_oracle)), &
         & MAXVAL(ABS(vsh_large - vsh_large_oracle)), &
         & ABS(energy_large - energy_large_oracle)
      err = MAX(err, fft_err)
      IF (ALLOCATED(cache_cp2k_plan%bvk_phase_forward) .OR. &
         & ALLOCATED(cache_cp2k_plan%bvk_phase_inverse)) err = MAX(err, 1.0_wp)
      ! The production grid proof and reduced image accumulation are also tested
      ! beyond the small-mesh dense phase oracle, with reverse arrival order.
      CALL cp2k_exchange_stream_begin(stream_large, calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, &
         & wfn_large, nmesh_large, kfrac_large, weights_large, error)
      IF (ALLOCATED(error)) RETURN
      IF (.NOT. cp2k_exchange_stream_has_full_mesh_storage(stream_large)) THEN
         err = MAX(err, 1.0_wp)
      END IF
      DO j = 1, PRODUCT(nmesh_large)
         ik = PRODUCT(nmesh_large) - j + 1
         CALL cp2k_exchange_stream_push(stream_large, ik, &
            & density_large(:, :, :, ik), overlap_large(:, :, ik), error)
         IF (ALLOCATED(error)) RETURN
      END DO
      CALL cp2k_exchange_stream_apply(stream_large, calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, &
         & vsh_large, energy_large, error)
      IF (ALLOCATED(error)) RETURN
      IF (.NOT. cp2k_exchange_stream_has_full_mesh_storage(stream_large)) THEN
         err = MAX(err, 1.0_wp)
      END IF
      DO ik = 1, PRODUCT(nmesh_large)
         CALL cp2k_exchange_stream_get(stream_large, ik, &
            & fock_large(:, :, :, ik), error, overlap_large(:, :, ik))
         IF (ALLOCATED(error)) RETURN
      END DO
      CALL cp2k_exchange_stream_end(stream_large, error)
      IF (ALLOCATED(error)) RETURN
      err = MAX(err, MAXVAL(ABS(fock_large - fock_large_oracle)), &
         & MAXVAL(ABS(vsh_large - vsh_large_oracle)), &
         & ABS(energy_large - energy_large_oracle))

      ! A non-divisor batch size on the 9x9x1 mesh exercises bounded storage,
      ! reverse k-block arrival, and a short final image batch beyond the dense
      ! small-grid phase oracle.
      fock_large = (0.0_wp, 0.0_wp)
      CALL cp2k_exchange_stream_begin(stream_batch_large, calc_cp2k, mol, &
         & cache_cp2k, cache_cp2k_plan, wfn_large, nmesh_large, kfrac_large, weights_large, error, &
         & mode=cp2k_exchange_stream_batched, batch_size=7)
      IF (ALLOCATED(error)) RETURN
      IF (cp2k_exchange_stream_has_full_mesh_storage(stream_batch_large)) THEN
         err = MAX(err, 1.0_wp)
      END IF
      IF (ALLOCATED(cache_cp2k_plan%bvk_phase_forward) &
         & .OR. ALLOCATED(cache_cp2k_plan%bvk_phase_inverse)) err = MAX(err, 1.0_wp)
      peak_expected = 4_int64*INT(nao, int64)*INT(nao, int64) &
         & *7_int64*INT(wfn_large%nspin, int64) &
         & + 3_int64*INT(nao, int64)*INT(wfn_large%nspin, int64) &
         & + 2_int64*INT(PRODUCT(nmesh_large), int64) &
         & + 3_int64*INT(MAXVAL(nmesh_large), int64)
      peak_expected = peak_expected + MAX( &
         & 2_int64*INT(nao, int64)**2, &
         & 5_int64*INT(nao, int64)**2 + 6_int64*INT(nao, int64), &
         & INT(nao, int64)**2*INT(wfn_large%nspin, int64) &
         & + 4_int64*INT(nao, int64)**2 + 7_int64)
      IF (cp2k_exchange_stream_peak_complex_elements(stream_batch_large) &
         & /= peak_expected) err = MAX(err, 1.0_wp)
      batches_done = .FALSE.
      DO WHILE (.NOT. batches_done)
         CALL cp2k_exchange_stream_batch_bounds(stream_batch_large, first_image, &
            & last_image, error)
         IF (ALLOCATED(error)) RETURN
         IF (last_image - first_image + 1 > 7) err = MAX(err, 1.0_wp)
         IF (last_image == PRODUCT(nmesh_large) &
            & .AND. last_image - first_image + 1 /= 4) err = MAX(err, 1.0_wp)
         IF (cp2k_exchange_stream_peak_complex_elements(stream_batch_large) &
            & /= peak_expected) err = MAX(err, 1.0_wp)
         DO j = 1, PRODUCT(nmesh_large)
            ik = PRODUCT(nmesh_large) - j + 1
            CALL cp2k_exchange_stream_push(stream_batch_large, ik, &
               & density_large(:, :, :, ik), overlap_large(:, :, ik), error)
            IF (ALLOCATED(error)) RETURN
         END DO
         CALL cp2k_exchange_stream_batch_apply(stream_batch_large, calc_cp2k, &
            & mol, cache_cp2k, cache_cp2k_plan, error)
         IF (ALLOCATED(error)) RETURN
         DO j = 1, PRODUCT(nmesh_large)
            ik = PRODUCT(nmesh_large) - j + 1
            CALL cp2k_exchange_stream_batch_pull(stream_batch_large, ik, &
               & overlap_large(:, :, ik), &
               & fock_large(:, :, :, ik), error)
            IF (ALLOCATED(error)) RETURN
         END DO
         CALL cp2k_exchange_stream_batch_advance(stream_batch_large, &
            & batches_done, error)
         IF (ALLOCATED(error)) RETURN
      END DO
      CALL cp2k_exchange_stream_batch_result(stream_batch_large, calc_cp2k, mol, &
         & cache_cp2k, cache_cp2k_plan, vsh_large, energy_large, error)
      IF (ALLOCATED(error)) RETURN
      CALL cp2k_exchange_stream_end(stream_batch_large, error)
      IF (ALLOCATED(error)) RETURN
      err = MAX(err, MAXVAL(ABS(fock_large - fock_large_oracle)), &
         & MAXVAL(ABS(vsh_large - vsh_large_oracle)), &
         & ABS(energy_large - energy_large_oracle))

      ALLOCATE (overlap_grad_large(nao, nao, PRODUCT(nmesh_large)), &
         & overlap_grad_large_oracle(nao, nao, PRODUCT(nmesh_large)), &
         & mulliken_grad_large(nao, nao, PRODUCT(nmesh_large)), &
         & bocorr_grad_large(nao, nao, PRODUCT(nmesh_large)), &
         & gradient_large(3, mol%nat), gradient_large_oracle(3, mol%nat), &
         & sigma_large(3, 3), sigma_large_oracle(3, 3))
      SELECT TYPE (fock => exchange)
      TYPE is (exchange_fock)
         CALL get_KGrad_kmesh(fock, mol, ptr, oracle_plan, kernel_large, kfrac_large, &
            & weights_large, density_large, overlap_large, &
            & overlap_grad_large_oracle, mulliken_grad_large, bocorr_grad_large)
         gradient_large_oracle = 0.0_wp
         sigma_large_oracle = 0.0_wp
         CALL fock%get_bvk_Kmatrix_derivs(mol, kernel_large, &
            & mulliken_grad_large, bocorr_grad_large, gradient_large_oracle, &
            & sigma_large_oracle)
      END SELECT
      CALL cp2k_exchange_kmesh_gradient(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, &
         & nmesh_large, kfrac_large, weights_large, density_large, &
         & overlap_large, overlap_grad_large, gradient_large, sigma_large, error)
      IF (ALLOCATED(error)) RETURN
      err = MAX(err, MAXVAL(ABS(overlap_grad_large - overlap_grad_large_oracle)), &
         & MAXVAL(ABS(gradient_large - gradient_large_oracle)), &
         & MAXVAL(ABS(sigma_large - sigma_large_oracle)))
      ! Reverse mode uses the exact adjoint normalization of the same factorized
      ! transforms.  Atomic forces and homogeneous-strain response therefore
      ! provide end-to-end derivative oracles, not only transform-level checks.
      CALL cp2k_exchange_kmesh_gradient(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, &
         & nmesh_large, kfrac_large, weights_large, density_large, &
         & overlap_large, overlap_grad_large, gradient_large, sigma_large, error, &
         & transform_backend=cp2k_exchange_transform_separable_dft)
      IF (ALLOCATED(error)) RETURN
      err = MAX(err, MAXVAL(ABS(overlap_grad_large - overlap_grad_large_oracle)), &
         & MAXVAL(ABS(gradient_large - gradient_large_oracle)), &
         & MAXVAL(ABS(sigma_large - sigma_large_oracle)))
      separable_err = MAX(separable_err, &
         & MAXVAL(ABS(overlap_grad_large - overlap_grad_large_oracle)), &
         & MAXVAL(ABS(gradient_large - gradient_large_oracle)), &
         & MAXVAL(ABS(sigma_large - sigma_large_oracle)))
      CALL cp2k_exchange_kmesh_gradient(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, &
         & nmesh_large, kfrac_large, weights_large, density_large, &
         & overlap_large, overlap_grad_large, gradient_large, sigma_large, error, &
         & transform_backend=cp2k_exchange_transform_mixed_radix_fft)
      IF (ALLOCATED(error)) RETURN
      fft_err = MAX(fft_err, &
         & MAXVAL(ABS(overlap_grad_large - overlap_grad_large_oracle)), &
         & MAXVAL(ABS(gradient_large - gradient_large_oracle)), &
         & MAXVAL(ABS(sigma_large - sigma_large_oracle)))
      PRINT *, "MIXED_RADIX_FFT_UKS_REVERSE_RESID", &
         & MAXVAL(ABS(overlap_grad_large - overlap_grad_large_oracle)), &
         & MAXVAL(ABS(gradient_large - gradient_large_oracle)), &
         & MAXVAL(ABS(sigma_large - sigma_large_oracle))

      ! UKS, a shifted/permuted 9x9x1 mesh, and a non-divisor batch exercise the
      ! bounded reverse protocol beyond the small RKS oracle above.
      overlap_grad_large = (0.0_wp, 0.0_wp)
      CALL cp2k_exchange_stream_begin(stream_batch_large, calc_cp2k, mol, &
         & cache_cp2k, cache_cp2k_plan, wfn_large, nmesh_large, kfrac_large, weights_large, error, &
         & mode=cp2k_exchange_stream_batched, batch_size=7, &
         & reverse_only=.TRUE.)
      IF (ALLOCATED(error)) RETURN
      IF (cp2k_exchange_stream_has_full_mesh_storage(stream_batch_large)) THEN
         err = MAX(err, 1.0_wp)
      END IF
      peak_expected = 6_int64*INT(nao, int64)**2*7_int64 &
         & *INT(wfn_large%nspin, int64) + 5_int64*INT(nao, int64) &
         & *INT(wfn_large%nspin, int64) &
         & + 2_int64*INT(PRODUCT(nmesh_large), int64) &
         & + 3_int64*INT(MAXVAL(nmesh_large), int64) &
         & + 10_int64*INT(nao, int64)**2
      IF (cp2k_exchange_stream_peak_complex_elements(stream_batch_large) &
         & /= peak_expected) err = MAX(err, 1.0_wp)
      batches_done = .FALSE.
      DO WHILE (.NOT. batches_done)
         DO j = 1, PRODUCT(nmesh_large)
            ik = PRODUCT(nmesh_large) - j + 1
            CALL cp2k_exchange_stream_reverse_push(stream_batch_large, ik, &
               & density_large(:, :, :, ik), overlap_large(:, :, ik), error)
            IF (ALLOCATED(error)) RETURN
         END DO
         CALL cp2k_exchange_stream_reverse_batch_apply(stream_batch_large, &
            & calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, error)
         IF (ALLOCATED(error)) RETURN
         DO j = 1, PRODUCT(nmesh_large)
            ik = PRODUCT(nmesh_large) - j + 1
            CALL cp2k_exchange_stream_reverse_batch_pull(stream_batch_large, ik, &
               & density_large(:, :, :, ik), overlap_large(:, :, ik), &
               & overlap_grad_large(:, :, ik), error)
            IF (ALLOCATED(error)) RETURN
         END DO
         CALL cp2k_exchange_stream_reverse_batch_advance(stream_batch_large, &
            & batches_done, error)
         IF (ALLOCATED(error)) RETURN
      END DO
      CALL cp2k_exchange_stream_reverse_result(stream_batch_large, calc_cp2k, &
         & mol, cache_cp2k, cache_cp2k_plan, gradient_large, sigma_large, error)
      IF (ALLOCATED(error)) RETURN
      peak_expected = cp2k_exchange_stream_peak_complex_elements( &
         & stream_batch_large)
      CALL cp2k_exchange_stream_end(stream_batch_large, error)
      IF (ALLOCATED(error)) RETURN
      err = MAX(err, MAXVAL(ABS(overlap_grad_large - overlap_grad_large_oracle)), &
         & MAXVAL(ABS(gradient_large - gradient_large_oracle)), &
         & MAXVAL(ABS(sigma_large - sigma_large_oracle)))
      stream_reverse_err = MAX(stream_reverse_err, &
         & MAXVAL(ABS(overlap_grad_large - overlap_grad_large_oracle)), &
         & MAXVAL(ABS(gradient_large - gradient_large_oracle)), &
         & MAXVAL(ABS(sigma_large - sigma_large_oracle)))

      PRINT *, "STREAM_REVERSE_UKS_BATCH_RESID", &
         & MAXVAL(ABS(overlap_grad_large - overlap_grad_large_oracle)), &
         & MAXVAL(ABS(gradient_large - gradient_large_oracle)), &
         & MAXVAL(ABS(sigma_large - sigma_large_oracle)), peak_expected

      ! Exercise every whole-mesh selector in both directions for an RKS mesh
      ! that is genuinely two-dimensional, shifted, and externally permuted.
      ! This follows the UKS stream checks so the helper's RKS onsite state is
      ! also the correct state for the 3D RKS stream oracle below.
      ALLOCATE (density_large_rks(nao, nao, 1, PRODUCT(nmesh_large)))
      density_large_rks(:, :, 1, :) = density_large(:, :, 1, :) &
         & + density_large(:, :, 2, :)
      CALL check_exchange_transform_selectors(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn, &
         & nmesh_large, kfrac_large, weights_large, density_large_rks, &
         & overlap_large, error, selector_separable_err, selector_fft_err, &
         & rks_2d_fft_forward, rks_2d_fft_reverse)
      IF (ALLOCATED(error)) RETURN
      separable_err = MAX(separable_err, selector_separable_err)
      fft_err = MAX(fft_err, selector_fft_err)
      PRINT *, "TRANSFORM_SELECTOR_2D_RKS_SEPARABLE_MAX_RESID", &
         & selector_separable_err
      PRINT *, "MIXED_RADIX_FFT_2D_RKS_FORWARD_RESID", rks_2d_fft_forward
      PRINT *, "MIXED_RADIX_FFT_2D_RKS_REVERSE_RESID", rks_2d_fft_reverse

      ! A genuinely three-dimensional, shifted 2x2x2 grid completes the mesh
      ! dimensionality coverage of the 3x1x1 RKS and 9x9x1 UKS cases.
      ALLOCATE (kfrac_cube(3, PRODUCT(nmesh_cube)), &
         & weights_cube(PRODUCT(nmesh_cube)), &
         & density_cube(nao, nao, 1, PRODUCT(nmesh_cube)), &
         & density_cube_uks(nao, nao, 2, PRODUCT(nmesh_cube)), &
         & overlap_cube(nao, nao, PRODUCT(nmesh_cube)), &
         & overlap_grad_cube(nao, nao, PRODUCT(nmesh_cube)), &
         & overlap_grad_cube_oracle(nao, nao, PRODUCT(nmesh_cube)), &
         & gradient_cube(3, mol%nat), gradient_cube_oracle(3, mol%nat), &
         & sigma_cube(3, 3), sigma_cube_oracle(3, 3))
      weights_cube = 1.0_wp/REAL(PRODUCT(nmesh_cube), wp)
      DO ik = 1, PRODUCT(nmesh_cube)
         j = MODULO(5*(ik - 1) + 3, PRODUCT(nmesh_cube))
         kfrac_cube(1, ik) = 0.071_wp &
            & + REAL(MODULO(j, nmesh_cube(1)), wp)/REAL(nmesh_cube(1), wp)
         kfrac_cube(2, ik) = 0.043_wp &
            & + REAL(MODULO(j/nmesh_cube(1), nmesh_cube(2)), wp) &
            & /REAL(nmesh_cube(2), wp)
         kfrac_cube(3, ik) = 0.019_wp &
            & + REAL(j/(nmesh_cube(1)*nmesh_cube(2)), wp) &
            & /REAL(nmesh_cube(3), wp)
         i = 1 + MODULO(j, ncell)
         density_cube(:, :, 1, ik) = density_mesh(:, :, 1, i)
         overlap_cube(:, :, ik) = overlap_mesh(:, :, i)
         ii = 1 + MODULO(11*j + 2, PRODUCT(nmesh_large))
         density_cube_uks(:, :, :, ik) = density_large(:, :, :, ii)
      END DO
      CALL cp2k_exchange_kmesh_gradient(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, nmesh_cube, &
         & kfrac_cube, weights_cube, density_cube, overlap_cube, &
         & overlap_grad_cube_oracle, gradient_cube_oracle, sigma_cube_oracle, &
         & error)
      IF (ALLOCATED(error)) RETURN
      overlap_grad_cube = (0.0_wp, 0.0_wp)
      CALL cp2k_exchange_stream_begin(stream_reverse_batch, calc_cp2k, mol, &
         & cache_cp2k, cache_cp2k_plan, wfn, nmesh_cube, kfrac_cube, weights_cube, error, &
         & mode=cp2k_exchange_stream_batched, batch_size=3, &
         & reverse_only=.TRUE.)
      IF (ALLOCATED(error)) RETURN
      peak_expected = 6_int64*INT(nao, int64)**2*3_int64 &
         & *INT(wfn%nspin, int64) + 5_int64*INT(nao, int64) &
         & *INT(wfn%nspin, int64) + 2_int64*INT(PRODUCT(nmesh_cube), int64) &
         & + 3_int64*INT(MAXVAL(nmesh_cube), int64) + 10_int64*INT(nao, int64)**2
      IF (cp2k_exchange_stream_peak_complex_elements(stream_reverse_batch) &
         & /= peak_expected) err = MAX(err, 1.0_wp)
      batches_done = .FALSE.
      DO WHILE (.NOT. batches_done)
         DO ik = 1, PRODUCT(nmesh_cube)
            CALL cp2k_exchange_stream_reverse_push(stream_reverse_batch, ik, &
               & density_cube(:, :, :, ik), overlap_cube(:, :, ik), error)
            IF (ALLOCATED(error)) RETURN
         END DO
         CALL cp2k_exchange_stream_reverse_batch_apply(stream_reverse_batch, &
            & calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, error)
         IF (ALLOCATED(error)) RETURN
         DO ik = 1, PRODUCT(nmesh_cube)
            CALL cp2k_exchange_stream_reverse_batch_pull(stream_reverse_batch, &
               & ik, density_cube(:, :, :, ik), overlap_cube(:, :, ik), &
               & overlap_grad_cube(:, :, ik), error)
            IF (ALLOCATED(error)) RETURN
         END DO
         CALL cp2k_exchange_stream_reverse_batch_advance(stream_reverse_batch, &
            & batches_done, error)
         IF (ALLOCATED(error)) RETURN
      END DO
      CALL cp2k_exchange_stream_reverse_result(stream_reverse_batch, calc_cp2k, &
         & mol, cache_cp2k, cache_cp2k_plan, gradient_cube, sigma_cube, error)
      IF (ALLOCATED(error)) RETURN
      CALL cp2k_exchange_stream_end(stream_reverse_batch, error)
      IF (ALLOCATED(error)) RETURN
      err = MAX(err, MAXVAL(ABS(overlap_grad_cube - overlap_grad_cube_oracle)), &
         & MAXVAL(ABS(gradient_cube - gradient_cube_oracle)), &
         & MAXVAL(ABS(sigma_cube - sigma_cube_oracle)))
      stream_reverse_err = MAX(stream_reverse_err, &
         & MAXVAL(ABS(overlap_grad_cube - overlap_grad_cube_oracle)), &
         & MAXVAL(ABS(gradient_cube - gradient_cube_oracle)), &
         & MAXVAL(ABS(sigma_cube - sigma_cube_oracle)))
      PRINT *, "STREAM_REVERSE_3D_RESID", &
         & MAXVAL(ABS(overlap_grad_cube - overlap_grad_cube_oracle)), &
         & MAXVAL(ABS(gradient_cube - gradient_cube_oracle)), &
         & MAXVAL(ABS(sigma_cube - sigma_cube_oracle)), peak_expected

      ! Complete the selector matrix with a shifted, permuted, genuinely 3D UKS
      ! mesh.  Both spin channels are nonzero and distinct.
      CALL check_exchange_transform_selectors(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, &
         & wfn_large, nmesh_cube, kfrac_cube, weights_cube, density_cube_uks, &
         & overlap_cube, error, selector_separable_err, selector_fft_err, &
         & uks_3d_fft_forward, uks_3d_fft_reverse)
      IF (ALLOCATED(error)) RETURN
      separable_err = MAX(separable_err, selector_separable_err)
      fft_err = MAX(fft_err, selector_fft_err)
      PRINT *, "TRANSFORM_SELECTOR_3D_UKS_SEPARABLE_MAX_RESID", &
         & selector_separable_err
      PRINT *, "MIXED_RADIX_FFT_3D_UKS_FORWARD_RESID", uks_3d_fft_forward
      PRINT *, "MIXED_RADIX_FFT_3D_UKS_REVERSE_RESID", uks_3d_fft_reverse

      ! The public partial-accumulator ABI must reproduce the unchanged dense
      ! forward and reverse oracles for one, two, and four host partitions.  Four
      ! partitions on the three-point mesh deliberately include an empty owner.
      CALL check_partial_exchange_case(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn, nmesh, &
         & kfrac, weights, density_mesh, overlap_mesh, [1, 2, 4], .TRUE., error, &
         & partial_err)
      IF (ALLOCATED(error)) RETURN
      stream_reverse_err = MAX(stream_reverse_err, partial_err)
      CALL check_partial_exchange_case(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn, nmesh, &
         & kfrac_order, weights_order, density_order, overlap_order, [4], &
         & .FALSE., error, partial_err)
      IF (ALLOCATED(error)) RETURN
      stream_reverse_err = MAX(stream_reverse_err, partial_err)
      CALL check_partial_exchange_case(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn_large, &
         & nmesh_large, kfrac_large, weights_large, density_large, overlap_large, &
         & [4], .FALSE., error, partial_err)
      IF (ALLOCATED(error)) RETURN
      stream_reverse_err = MAX(stream_reverse_err, partial_err)
      CALL check_partial_exchange_case(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn, &
         & nmesh_cube, kfrac_cube, weights_cube, density_cube, overlap_cube, &
         & [2], .FALSE., error, partial_err)
      IF (ALLOCATED(error)) RETURN
      stream_reverse_err = MAX(stream_reverse_err, partial_err)
      CALL check_partial_exchange_multibatch(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn, &
         & nmesh, kfrac, weights, density_mesh, overlap_mesh, error, partial_err)
      IF (ALLOCATED(error)) RETURN
      stream_reverse_err = MAX(stream_reverse_err, partial_err)
      CALL check_partial_exchange_distributed_images(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, &
         & wfn_large, nmesh_large, kfrac_large, weights_large, density_large, &
         & overlap_large, error, partial_err)
      IF (ALLOCATED(error)) RETURN
      stream_reverse_err = MAX(stream_reverse_err, partial_err)
      PRINT *, "DISTRIBUTED_IMAGE_KERNEL_9X9X1_UKS_MAX_RESID", partial_err
      CALL check_partial_exchange_range_local(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, &
         & wfn_large, nmesh_large, kfrac_large, weights_large, density_large, &
         & overlap_large, error, partial_err)
      IF (ALLOCATED(error)) RETURN
      stream_reverse_err = MAX(stream_reverse_err, partial_err)
      PRINT *, "RANGE_LOCAL_IMAGE_IMPORT_9X9X1_UKS_MAX_RESID", partial_err
      CALL check_partial_exchange_range_local(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, &
         & wfn_large, nmesh_large, kfrac_large, weights_large, density_large, &
         & overlap_large, error, partial_err, &
         & cp2k_exchange_transform_mixed_radix_fft)
      IF (ALLOCATED(error)) RETURN
      stream_reverse_err = MAX(stream_reverse_err, partial_err)
      PRINT *, "RANGE_LOCAL_MIXED_RADIX_9X9X1_UKS_MAX_RESID", partial_err
      CALL check_partial_exchange_distributed_images(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, &
         & wfn, nmesh_cube, kfrac_cube, weights_cube, density_cube, overlap_cube, &
         & error, partial_err)
      IF (ALLOCATED(error)) RETURN
      stream_reverse_err = MAX(stream_reverse_err, partial_err)
      PRINT *, "DISTRIBUTED_IMAGE_KERNEL_2X2X2_MAX_RESID", partial_err
      CALL check_partial_exchange_range_local(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn, &
         & nmesh_cube, kfrac_cube, weights_cube, density_cube, overlap_cube, &
         & error, partial_err)
      IF (ALLOCATED(error)) RETURN
      stream_reverse_err = MAX(stream_reverse_err, partial_err)
      PRINT *, "RANGE_LOCAL_IMAGE_IMPORT_2X2X2_MAX_RESID", partial_err
      CALL check_partial_exchange_range_local(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn, &
         & nmesh_cube, kfrac_cube, weights_cube, density_cube, overlap_cube, &
         & error, partial_err, cp2k_exchange_transform_mixed_radix_fft)
      IF (ALLOCATED(error)) RETURN
      stream_reverse_err = MAX(stream_reverse_err, partial_err)
      PRINT *, "RANGE_LOCAL_MIXED_RADIX_2X2X2_MAX_RESID", partial_err
      CALL check_partial_exchange_distributed_images(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, &
         & wfn, nmesh, kfrac, weights, density_mesh, overlap_mesh, error, &
         & partial_err)
      IF (ALLOCATED(error)) RETURN
      stream_reverse_err = MAX(stream_reverse_err, partial_err)
      PRINT *, "DISTRIBUTED_IMAGE_KERNEL_3X1X1_MAX_RESID", partial_err
      CALL check_partial_exchange_range_local(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn, &
         & nmesh, kfrac, weights, density_mesh, overlap_mesh, error, partial_err)
      IF (ALLOCATED(error)) RETURN
      stream_reverse_err = MAX(stream_reverse_err, partial_err)
      PRINT *, "RANGE_LOCAL_IMAGE_IMPORT_3X1X1_MAX_RESID", partial_err
      PRINT *, "PARTIAL_K_TO_R_MAX_RESID", stream_reverse_err

      bvk_builds = cache_cp2k_plan%bvk_plan_builds
      CALL cp2k_exchange_kmesh(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn_large, nmesh_large, &
         & kfrac_large, weights_large, density_large, overlap_large, fock_large, &
         & vsh_large, energy_large, error)
      IF (ALLOCATED(error)) RETURN
      IF (cache_cp2k_plan%bvk_plan_builds /= bvk_builds + 1) err = MAX(err, 1.0_wp)
      err = MAX(err, MAXVAL(ABS(fock_large - fock_large_oracle)), &
         & MAXVAL(ABS(vsh_large - vsh_large_oracle)), &
         & ABS(energy_large - energy_large_oracle))

      CALL check(error, separable_err, 0.0_wp, thr=1.0e-11_wp)
      IF (ALLOCATED(error)) THEN
         CALL test_failed(error, &
            & "Separable exchange exceeds the dense-oracle tolerance")
         RETURN
      END IF
      PRINT *, "MIXED_RADIX_FFT_EXCHANGE_MAX_RESID", fft_err
      CALL check(error, fft_err, 0.0_wp, thr=1.0e-11_wp)
      IF (ALLOCATED(error)) THEN
         CALL test_failed(error, &
            & "Mixed-radix FFT exchange exceeds the dense-oracle tolerance")
         RETURN
      END IF
      CALL check(error, stream_reverse_err, 0.0_wp, thr=1.0e-11_wp)
      IF (ALLOCATED(error)) THEN
         CALL test_failed(error, &
            & "Streamed reverse exchange exceeds the dense-oracle tolerance")
         RETURN
      END IF
      PRINT *, "BVK_EXCHANGE_SUPERCELL_MAX_RESID", err
      CALL check(error, err, 0.0_wp, thr=1.0e-9_wp)
      IF (ALLOCATED(error)) THEN
         CALL test_failed(error, &
            & "Whole-mesh exchange does not reproduce the explicit BvK supercell")
      END IF
   END SUBROUTINE test_bvk_exchange_kernel

! **************************************************************************************************
!> \brief ...
!> \param error ...
! **************************************************************************************************
   SUBROUTINE test_kpoint_fock_complex(error)
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      TYPE(structure_type)                               :: mol

      class(basis_type), ALLOCATABLE :: bas
      class(qvszp_basis_type), ALLOCATABLE :: qvszp_bas
      class(exchange_type), ALLOCATABLE :: exchange
      TYPE(xtb_calculator) :: calc_cp2k
      TYPE(basis_cache) :: bcache
      TYPE(container_cache) :: cache, cache_cp2k
      TYPE(cp2k_stream_plan_type) :: cache_cp2k_plan
      TYPE(exchange_cache), POINTER :: ptr
      TYPE(exchange_bvk_kernel) :: kernel_mesh
      TYPE(cp2k_stream_plan_type) :: oracle_plan
      TYPE(wavefunction_type) :: wfn, wfn_aux
      INTEGER, PARAMETER :: num(3) = [1, 1, 8]
      REAL(wp), PARAMETER :: xyz(3, 3) = RESHAPE([ &
         & 0.20_wp, 0.30_wp, 0.40_wp, &
         & 0.96_wp, 0.30_wp, 0.90_wp, &
         & 0.20_wp, 0.30_wp, 1.40_wp], [3, 3])/autoaa
      REAL(wp), PARAMETER :: lattice(3, 3) = RESHAPE([ &
         & 7.80_wp, 0.20_wp, 0.10_wp, &
         & 0.40_wp, 8.30_wp, 0.30_wp, &
         & 0.10_wp, 0.50_wp, 9.10_wp], [3, 3])/autoaa
      LOGICAL, PARAMETER :: periodic(3) = [.TRUE., .TRUE., .TRUE.]
      INTEGER :: i, iao, iat, ii, is, ish, izp, j, jao, jat, jj, js, jsh, jzp, nao
      REAL(wp) :: deriv_fd, deriv_ref, energy, energy_m, energy_p, energy_ref
      REAL(wp) :: err, q0, step, theta
      REAL(wp) :: coord_direction(3, 3), identity(3, 3), lattice0(3, 3), &
         & strain_direction(3, 3), transform(3, 3), xyz0(3, 3)
      REAL(wp), ALLOCATABLE :: lattr(:, :), overlap(:, :), vsh(:), vsh_m(:), &
         & vsh_p(:), vsh_tr(:), vsh_u(:), ao_grad_gamma(:, :), &
         & mulliken_grad(:, :), mulliken_grad_gamma(:, :), bocorr_grad(:, :), &
         & bocorr_grad_gamma(:, :), bocorr_grad_mesh_r(:, :, :), &
         & gradient_cp_k(:, :), gradient_cp_mesh(:, :), gradient_gamma(:, :), &
         & gradient_k(:, :), gradient_mesh(:, :), &
         & kernel_backup(:, :), kernel_direction(:, :), sigma_gamma(:, :), sigma_k(:, :)
      REAL(wp), ALLOCATABLE :: mulliken_grad_mesh_r(:, :, :), &
         & sigma_cp_k(:, :), sigma_cp_mesh(:, :), sigma_mesh(:, :), vsh_mesh(:)
      COMPLEX(wp), ALLOCATABLE :: density_k(:, :, :), density_m(:, :, :), &
         & density_p(:, :, :), density_u(:, :, :), direction(:, :, :), &
         & density_uks(:, :, :), &
         & fock_k(:, :, :), fock_m(:, :, :), fock_p(:, :, :), fock_tr(:, :, :), &
         & fock_u(:, :, :), fock_uks(:, :, :), overlap_grad(:, :), &
         & overlap_grad_cp_k(:, :), &
         & overlap_grad_uks(:, :), overlap_k(:, :), overlap_m(:, :), &
         & overlap_p(:, :), overlap_u(:, :), unitary(:)
      COMPLEX(wp), ALLOCATABLE :: density_mesh(:, :, :, :), &
         & fock_mesh(:, :, :, :), overlap_grad_cp_mesh(:, :, :), &
         & overlap_grad_mesh(:, :, :), &
         & overlap_mesh(:, :, :)

      CALL new(mol, num, xyz, lattice=lattice, periodic=periodic)
      ALLOCATE (qvszp_bas)
      CALL make_qvszp_basis(qvszp_bas, mol, error)
      IF (ALLOCATED(error)) RETURN
      CALL move_alloc(qvszp_bas, bas)

      CALL new_wavefunction(wfn_aux, mol%nat, 0, 0, 1, 0.0_wp, .FALSE.)
      CALL get_eeqbc_charges(mol, error, wfn_aux%qat(:, 1))
      IF (ALLOCATED(error)) RETURN
      CALL bas%update(mol, bcache, .FALSE., wfn_aux=wfn_aux)

      nao = bas%nao
      CALL get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
      ALLOCATE (overlap(nao, nao))
      CALL get_overlap(mol, lattr, cutoff, bas, bcache, overlap)
      CALL new_wavefunction(wfn, mol%nat, bas%nsh, nao, 1, 0.0_wp, .FALSE.)
      DO i = 1, nao
         DO j = 1, nao
            wfn%density(i, j, 1) = 0.2_wp/(1.0_wp + ABS(i - j))
         END DO
      END DO
      wfn%qsh = 0.0_wp

      CALL make_exchange_gxtb(exchange, mol, bas)
      CALL exchange%update(mol, cache)
      CALL view(cache, ptr)
      CALL exchange%get_onsite_Kmatrix(mol, wfn, ptr)
      ALLOCATE (ptr%prev_F(nao, nao, 1), ptr%prev_vsh(bas%nsh, 1))
      CALL exchange%get_KFock(mol, ptr, wfn%density, overlap)
      ALLOCATE (calc_cp2k%exchange, source=exchange)
      CALL cp2k_prepare_exchange(calc_cp2k, mol, cache_cp2k, wfn, error)
      IF (ALLOCATED(error)) RETURN

      ALLOCATE (density_k(nao, nao, 1), overlap_k(nao, nao), &
         & fock_k(nao, nao, 1), vsh(bas%nsh))
      density_k = CMPLX(wfn%density, 0.0_wp, wp)
      overlap_k = CMPLX(overlap, 0.0_wp, wp)
      SELECT TYPE (fock => exchange)
      TYPE is (exchange_fock)
         CALL get_KFock_kpoint(fock, mol, ptr, density_k, overlap_k, &
            & fock_k, vsh, energy)
         energy_ref = 0.5_wp*SUM(wfn%density*ptr%prev_F)
         err = MAX(MAXVAL(ABS(REAL(fock_k, wp) - ptr%prev_F)), &
            & MAXVAL(ABS(AIMAG(fock_k))), MAXVAL(ABS(vsh - ptr%prev_vsh(:, 1))), &
            & ABS(energy - energy_ref))
         CALL check(error, err, 0.0_wp, thr=1.0e-11_wp)
         IF (ALLOCATED(error)) THEN
            CALL test_failed(error, "Complex exchange does not reproduce Gamma")
            RETURN
         END IF

         CALL fock%get_bvk_Kmatrix(mol, [1, 1, 1], kernel_mesh)
         ALLOCATE (density_mesh(nao, nao, 1, 1), overlap_mesh(nao, nao, 1), &
            & fock_mesh(nao, nao, 1, 1), vsh_mesh(bas%nsh))
         density_mesh(:, :, :, 1) = density_k
         overlap_mesh(:, :, 1) = overlap_k
         CALL get_KFock_kmesh(fock, mol, ptr, oracle_plan, kernel_mesh, &
            & RESHAPE([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]), [1.0_wp], &
            & density_mesh, overlap_mesh, fock_mesh, vsh_mesh, energy_ref)
         err = MAX(MAXVAL(ABS(fock_mesh(:, :, :, 1) - fock_k)), &
            & MAXVAL(ABS(vsh_mesh - vsh)), ABS(energy_ref - energy))
         CALL check(error, err, 0.0_wp, thr=1.0e-11_wp)
         IF (ALLOCATED(error)) THEN
            CALL test_failed(error, "Whole-mesh exchange does not reproduce Gamma")
            RETURN
         END IF

         ! At Gamma the direct overlap derivative and tblite's historical
         ! energy-weighted-density correction must obey G_S = -W_exchange.
         ! CP2K's independent-pair Pulay kernel supplies its own factor of two,
         ! so this relation fixes the sign used when folding the response.
         ALLOCATE (ao_grad_gamma(nao, nao), mulliken_grad(nao, nao), &
            & mulliken_grad_gamma(nao, nao), bocorr_grad(nao, nao), &
            & bocorr_grad_gamma(nao, nao), gradient_gamma(3, mol%nat), &
            & gradient_k(3, mol%nat), sigma_gamma(3, 3), sigma_k(3, 3), &
            & source=0.0_wp)
         ALLOCATE (overlap_grad(nao, nao), source=(0.0_wp, 0.0_wp))
         CALL fock%get_KGrad(mol, ptr, wfn%density, overlap, &
            & mulliken_grad_gamma, bocorr_grad_gamma, ao_grad_gamma)
         CALL get_KGrad_kpoint(fock, mol, ptr, density_k, overlap_k, &
            & mulliken_grad, bocorr_grad, overlap_grad)

         CALL fock%get_mulliken_derivs(mol, ptr, mulliken_grad_gamma, &
            & gradient_gamma, sigma_gamma)
         CALL fock%get_bocorr_derivs(mol, ptr, bocorr_grad_gamma, &
            & gradient_gamma, sigma_gamma)
         CALL get_mulliken_derivs_direct(fock, mol, ptr, mulliken_grad, &
            & gradient_k, sigma_k)
         CALL get_bocorr_derivs_direct(fock, mol, ptr, bocorr_grad, &
            & gradient_k, sigma_k)
         err = MAX(MAXVAL(ABS(overlap_grad + ao_grad_gamma)), &
            & MAXVAL(ABS(gradient_k - gradient_gamma)), &
            & MAXVAL(ABS(sigma_k - sigma_gamma)))

         ! The whole-mesh reverse mode must reduce exactly to the independently
         ! validated single-k/Gamma implementation at N=1.
         ALLOCATE (overlap_grad_mesh(nao, nao, 1), &
            & mulliken_grad_mesh_r(nao, nao, 1), &
            & bocorr_grad_mesh_r(nao, nao, 1), gradient_mesh(3, mol%nat), &
            & sigma_mesh(3, 3))
         CALL get_KGrad_kmesh(fock, mol, ptr, oracle_plan, kernel_mesh, &
            & RESHAPE([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]), [1.0_wp], &
            & density_mesh, overlap_mesh, overlap_grad_mesh, &
            & mulliken_grad_mesh_r, bocorr_grad_mesh_r)
         gradient_mesh = 0.0_wp
         sigma_mesh = 0.0_wp
         CALL fock%get_bvk_Kmatrix_derivs(mol, kernel_mesh, &
            & mulliken_grad_mesh_r, bocorr_grad_mesh_r, gradient_mesh, sigma_mesh)
         err = MAX(err, MAXVAL(ABS(overlap_grad_mesh(:, :, 1) - overlap_grad)), &
            & MAXVAL(ABS(gradient_mesh - gradient_k)), &
            & MAXVAL(ABS(sigma_mesh - sigma_k)))

         ! Exercise the public CP2K ABI as well: both wrappers must expose the
         ! same N=1 functional, including overlap sign and strain convention.
         ALLOCATE (overlap_grad_cp_k(nao, nao), &
            & overlap_grad_cp_mesh(nao, nao, 1), &
            & gradient_cp_k(3, mol%nat), gradient_cp_mesh(3, mol%nat), &
            & sigma_cp_k(3, 3), sigma_cp_mesh(3, 3))
         CALL cp2k_exchange_kpoint_gradient(calc_cp2k, mol, cache_cp2k, &
            & density_k, overlap_k, overlap_grad_cp_k, gradient_cp_k, &
            & sigma_cp_k, error)
         IF (ALLOCATED(error)) RETURN
         CALL cp2k_exchange_kmesh_gradient(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, &
            & [1, 1, 1], RESHAPE([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]), &
            & [1.0_wp], density_mesh, overlap_mesh, overlap_grad_cp_mesh, &
            & gradient_cp_mesh, sigma_cp_mesh, error)
         IF (ALLOCATED(error)) RETURN
         err = MAX(err, MAXVAL(ABS(overlap_grad_cp_k - overlap_grad)), &
            & MAXVAL(ABS(overlap_grad_cp_mesh(:, :, 1) - overlap_grad_cp_k)), &
            & MAXVAL(ABS(gradient_cp_mesh - gradient_cp_k)), &
            & MAXVAL(ABS(sigma_cp_mesh - sigma_cp_k)))
         CALL check(error, err, 0.0_wp, thr=1.0e-10_wp)
         IF (ALLOCATED(error)) THEN
            CALL test_failed(error, "Complex exchange gradient does not reproduce Gamma")
            RETURN
         END IF

         DO i = 1, nao
            DO j = 1, i - 1
               overlap_k(j, i) = overlap_k(j, i) &
                  & + CMPLX(0.0_wp, 1.0e-3_wp*SIN(REAL(3*i + 5*j, wp)), wp)
               overlap_k(i, j) = CONJG(overlap_k(j, i))
               density_k(j, i, 1) = density_k(j, i, 1) &
                  & + CMPLX(0.0_wp, 2.0e-3_wp*COS(REAL(7*i + 2*j, wp)), wp)
               density_k(i, j, 1) = CONJG(density_k(j, i, 1))
            END DO
         END DO
         CALL get_KFock_kpoint(fock, mol, ptr, density_k, overlap_k, &
            & fock_k, vsh, energy)
         err = MAXVAL(ABS(fock_k(:, :, 1) - CONJG(TRANSPOSE(fock_k(:, :, 1)))))

         ALLOCATE (fock_tr(nao, nao, 1), vsh_tr(bas%nsh))
         CALL get_KFock_kpoint(fock, mol, ptr, CONJG(density_k), &
            & CONJG(overlap_k), fock_tr, vsh_tr, energy_ref)
         err = MAX(err, MAXVAL(ABS(fock_tr - CONJG(fock_k))), &
            & MAXVAL(ABS(vsh_tr - vsh)), ABS(energy_ref - energy))

         ALLOCATE (unitary(nao), density_u(nao, nao, 1), overlap_u(nao, nao), &
            & fock_u(nao, nao, 1), vsh_u(bas%nsh))
         DO i = 1, nao
            theta = 0.071_wp*REAL(i, wp)
            unitary(i) = EXP(CMPLX(0.0_wp, theta, wp))
         END DO
         DO i = 1, nao
            DO j = 1, nao
               overlap_u(i, j) = CONJG(unitary(i))*overlap_k(i, j)*unitary(j)
               density_u(i, j, 1) = CONJG(unitary(i)) &
                                   & *density_k(i, j, 1)*unitary(j)
            END DO
         END DO
         CALL get_KFock_kpoint(fock, mol, ptr, density_u, overlap_u, &
            & fock_u, vsh_u, energy_ref)
         DO i = 1, nao
            DO j = 1, nao
               err = MAX(err, ABS(fock_u(i, j, 1) &
                  & - CONJG(unitary(i))*fock_k(i, j, 1)*unitary(j)))
            END DO
         END DO
         err = MAX(err, MAXVAL(ABS(vsh_u - vsh)), ABS(energy_ref - energy))

         ALLOCATE (direction(nao, nao, 1), density_m(nao, nao, 1), &
            & density_p(nao, nao, 1), fock_m(nao, nao, 1), &
            & fock_p(nao, nao, 1), source=(0.0_wp, 0.0_wp))
         ALLOCATE (vsh_m(bas%nsh), vsh_p(bas%nsh))
         DO i = 1, nao
            direction(i, i, 1) = CMPLX(1.0e-2_wp*SIN(REAL(11*i, wp)), &
               & 0.0_wp, wp)
            DO j = 1, i - 1
               direction(j, i, 1) = CMPLX( &
                  & 1.0e-2_wp*COS(REAL(3*i + j, wp)), &
                  & 1.0e-2_wp*SIN(REAL(i + 4*j, wp)), wp)
               direction(i, j, 1) = CONJG(direction(j, i, 1))
            END DO
         END DO
         step = 1.0e-5_wp
         density_m = density_k - step*direction
         density_p = density_k + step*direction
         CALL get_KFock_kpoint(fock, mol, ptr, density_m, overlap_k, &
            & fock_m, vsh_m, energy_m)
         CALL get_KFock_kpoint(fock, mol, ptr, density_p, overlap_k, &
            & fock_p, vsh_p, energy_p)
         deriv_fd = (energy_p - energy_m)/(2.0_wp*step)
         deriv_ref = REAL(SUM(CONJG(direction)*fock_k), wp)
         err = MAX(err, ABS(deriv_fd - deriv_ref))

         ! Reverse the complex exchange functional with respect to S(k).  The
         ! overlap response is the missing Pulay ingredient of a k-point force;
         ! validate it independently of any real-space integral convention.
         ALLOCATE (overlap_m(nao, nao), overlap_p(nao, nao))
         CALL get_KGrad_kpoint(fock, mol, ptr, density_k, overlap_k, &
            & mulliken_grad, bocorr_grad, overlap_grad)
         overlap_m = overlap_k - step*direction(:, :, 1)
         overlap_p = overlap_k + step*direction(:, :, 1)
         CALL get_KFock_kpoint(fock, mol, ptr, density_k, overlap_m, &
            & fock_m, vsh_m, energy_m)
         CALL get_KFock_kpoint(fock, mol, ptr, density_k, overlap_p, &
            & fock_p, vsh_p, energy_p)
         deriv_fd = (energy_p - energy_m)/(2.0_wp*step)
         deriv_ref = REAL(SUM(CONJG(overlap_grad)*direction(:, :, 1)), wp)
         err = MAX(err, ABS(deriv_fd - deriv_ref))
         err = MAX(err, MAXVAL(ABS(overlap_grad &
            & - CONJG(TRANSPOSE(overlap_grad)))))

         ! The geometry-dependent Mulliken kernel is real and symmetric.  Check
         ! the AO-resolved reverse response before it is contracted with the
         ! pairwise kernel derivative in get_mulliken_derivs.
         ALLOCATE (kernel_backup(fock%nsh, fock%nsh), &
            & kernel_direction(fock%nsh, fock%nsh), source=0.0_wp)
         DO ish = 1, fock%nsh
            kernel_direction(ish, ish) = 1.0e-2_wp*SIN(REAL(13*ish, wp))
            DO jsh = 1, ish - 1
               kernel_direction(jsh, ish) = &
                  & 1.0e-2_wp*COS(REAL(5*ish + 7*jsh, wp))
               kernel_direction(ish, jsh) = kernel_direction(jsh, ish)
            END DO
         END DO
         kernel_backup = ptr%g_mulliken
         ptr%g_mulliken = kernel_backup - step*kernel_direction
         CALL get_KFock_kpoint(fock, mol, ptr, density_k, overlap_k, &
            & fock_m, vsh_m, energy_m)
         ptr%g_mulliken = kernel_backup + step*kernel_direction
         CALL get_KFock_kpoint(fock, mol, ptr, density_k, overlap_k, &
            & fock_p, vsh_p, energy_p)
         ptr%g_mulliken = kernel_backup
         deriv_fd = (energy_p - energy_m)/(2.0_wp*step)
         deriv_ref = 0.0_wp
         DO ish = 1, fock%nsh
            ii = fock%iao_sh(ish)
            DO jsh = 1, ish
               jj = fock%iao_sh(jsh)
               DO iao = 1, fock%nao_sh(ish)
                  DO jao = 1, fock%nao_sh(jsh)
                     deriv_ref = deriv_ref + kernel_direction(jsh, ish) &
                                & *mulliken_grad(jj + jao, ii + iao)
                  END DO
               END DO
            END DO
         END DO
         err = MAX(err, ABS(deriv_fd - deriv_ref))

         ! Repeat the same directional check for the atom-block bond-order
         ! correction kernel.  As above, the response stores one derivative for
         ! each independent symmetric block, so only one triangle is contracted.
         DEALLOCATE (kernel_backup, kernel_direction)
         ALLOCATE (kernel_backup(mol%nat, mol%nat), &
            & kernel_direction(mol%nat, mol%nat), source=0.0_wp)
         DO iat = 1, mol%nat
            DO jat = 1, iat - 1
               kernel_direction(jat, iat) = &
                  & 1.0e-2_wp*SIN(REAL(3*iat + 11*jat, wp))
               kernel_direction(iat, jat) = kernel_direction(jat, iat)
            END DO
         END DO
         kernel_backup = ptr%g_bocorr
         ptr%g_bocorr = kernel_backup - step*kernel_direction
         CALL get_KFock_kpoint(fock, mol, ptr, density_k, overlap_k, &
            & fock_m, vsh_m, energy_m)
         ptr%g_bocorr = kernel_backup + step*kernel_direction
         CALL get_KFock_kpoint(fock, mol, ptr, density_k, overlap_k, &
            & fock_p, vsh_p, energy_p)
         ptr%g_bocorr = kernel_backup
         deriv_fd = (energy_p - energy_m)/(2.0_wp*step)
         deriv_ref = 0.0_wp
         DO iat = 1, mol%nat
            izp = mol%id(iat)
            is = fock%ish_at(iat)
            DO jat = 1, iat - 1
               jzp = mol%id(jat)
               js = fock%ish_at(jat)
               DO ish = 1, fock%nsh_id(izp)
                  ii = fock%iao_sh(is + ish)
                  DO jsh = 1, fock%nsh_id(jzp)
                     jj = fock%iao_sh(js + jsh)
                     DO iao = 1, fock%nao_sh(is + ish)
                        DO jao = 1, fock%nao_sh(js + jsh)
                           deriv_ref = deriv_ref + kernel_direction(jat, iat) &
                                      & *bocorr_grad(jj + jao, ii + iao)
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
         err = MAX(err, ABS(deriv_fd - deriv_ref))

         ! Contract the direct kernel responses with their real-space derivatives
         ! and validate the complete explicit exchange geometry response for a
         ! genuinely complex Hermitian k-point density.  The overlap is held
         ! fixed here, so the separate overlap-response check above remains the
         ! sole Pulay contribution.
         gradient_k = 0.0_wp
         sigma_k = 0.0_wp
         CALL get_mulliken_derivs_direct(fock, mol, ptr, mulliken_grad, &
            & gradient_k, sigma_k)
         CALL get_bocorr_derivs_direct(fock, mol, ptr, bocorr_grad, &
            & gradient_k, sigma_k)

         DO iat = 1, mol%nat
            DO i = 1, 3
               coord_direction(i, iat) = 0.03_wp*SIN(REAL(7*i + 11*iat, wp))
            END DO
         END DO
         xyz0 = mol%xyz
         mol%xyz = xyz0 - step*coord_direction
         CALL ptr%update(mol)
         CALL fock%get_mulliken_Kmatrix(mol, ptr)
         CALL fock%get_bocorr_Kmatrix(mol, ptr)
         CALL get_KFock_kpoint(fock, mol, ptr, density_k, overlap_k, &
            & fock_m, vsh_m, energy_m)
         mol%xyz = xyz0 + step*coord_direction
         CALL ptr%update(mol)
         CALL fock%get_mulliken_Kmatrix(mol, ptr)
         CALL fock%get_bocorr_Kmatrix(mol, ptr)
         CALL get_KFock_kpoint(fock, mol, ptr, density_k, overlap_k, &
            & fock_p, vsh_p, energy_p)
         deriv_fd = (energy_p - energy_m)/(2.0_wp*step)
         deriv_ref = SUM(gradient_k*coord_direction)
         err = MAX(err, ABS(deriv_fd - deriv_ref))

         ! Apply a symmetric homogeneous strain to both Cartesian coordinates and
         ! lattice vectors.  This also covers diagonal self-image kernels, whose
         ! historical derivative representation carries a factor of two.
         identity = 0.0_wp
         DO i = 1, 3
            identity(i, i) = 1.0_wp
         END DO
         strain_direction = RESHAPE([ &
            & 0.07_wp, 0.02_wp, -0.01_wp, &
            & 0.02_wp, -0.04_wp, 0.03_wp, &
            & -0.01_wp, 0.03_wp, 0.05_wp], [3, 3])
         lattice0 = mol%lattice
         mol%xyz = xyz0
         transform = identity - step*strain_direction
         mol%xyz = MATMUL(transform, xyz0)
         mol%lattice = MATMUL(transform, lattice0)
         CALL ptr%update(mol)
         CALL fock%get_mulliken_Kmatrix(mol, ptr)
         CALL fock%get_bocorr_Kmatrix(mol, ptr)
         CALL get_KFock_kpoint(fock, mol, ptr, density_k, overlap_k, &
            & fock_m, vsh_m, energy_m)
         transform = identity + step*strain_direction
         mol%xyz = MATMUL(transform, xyz0)
         mol%lattice = MATMUL(transform, lattice0)
         CALL ptr%update(mol)
         CALL fock%get_mulliken_Kmatrix(mol, ptr)
         CALL fock%get_bocorr_Kmatrix(mol, ptr)
         CALL get_KFock_kpoint(fock, mol, ptr, density_k, overlap_k, &
            & fock_p, vsh_p, energy_p)
         deriv_fd = (energy_p - energy_m)/(2.0_wp*step)
         deriv_ref = SUM(sigma_k*strain_direction)
         err = MAX(err, ABS(deriv_fd - deriv_ref))

         mol%xyz = xyz0
         mol%lattice = lattice0
         CALL ptr%update(mol)
         CALL fock%get_mulliken_Kmatrix(mol, ptr)
         CALL fock%get_bocorr_Kmatrix(mol, ptr)

         ! The shell potential is the derivative of the same complex-k energy
         ! with respect to the charge-dependent onsite exchange kernels.
         q0 = wfn%qsh(1, 1)
         wfn%qsh(1, 1) = q0 - step
         CALL exchange%get_onsite_Kmatrix(mol, wfn, ptr)
         CALL get_KFock_kpoint(fock, mol, ptr, density_k, overlap_k, &
            & fock_m, vsh_m, energy_m)
         wfn%qsh(1, 1) = q0 + step
         CALL exchange%get_onsite_Kmatrix(mol, wfn, ptr)
         CALL get_KFock_kpoint(fock, mol, ptr, density_k, overlap_k, &
            & fock_p, vsh_p, energy_p)
         err = MAX(err, ABS((energy_p - energy_m)/(2.0_wp*step) - vsh(1)))
         wfn%qsh(1, 1) = q0
         CALL exchange%get_onsite_Kmatrix(mol, wfn, ptr)

         ! Exercise the unrestricted spin factor and accumulation of both spin
         ! channels in the common overlap response.
         ALLOCATE (density_uks(nao, nao, 2), fock_uks(nao, nao, 2), &
            & overlap_grad_uks(nao, nao))
         density_uks(:, :, 1) = 0.61_wp*density_k(:, :, 1)
         density_uks(:, :, 2) = 0.37_wp*density_k(:, :, 1)
         CALL get_KFock_kpoint(fock, mol, ptr, density_uks, overlap_k, &
            & fock_uks, vsh, energy)
         CALL get_KGrad_kpoint(fock, mol, ptr, density_uks, overlap_k, &
            & mulliken_grad, bocorr_grad, overlap_grad_uks)
         overlap_m = overlap_k - step*direction(:, :, 1)
         overlap_p = overlap_k + step*direction(:, :, 1)
         CALL get_KFock_kpoint(fock, mol, ptr, density_uks, overlap_m, &
            & fock_uks, vsh_m, energy_m)
         CALL get_KFock_kpoint(fock, mol, ptr, density_uks, overlap_p, &
            & fock_uks, vsh_p, energy_p)
         deriv_fd = (energy_p - energy_m)/(2.0_wp*step)
         deriv_ref = REAL(SUM(CONJG(overlap_grad_uks) &
                             & *direction(:, :, 1)), wp)
         err = MAX(err, ABS(deriv_fd - deriv_ref))
      class default
         CALL test_failed(error, "Unexpected exchange implementation")
         RETURN
      END SELECT

      CALL check(error, err, 0.0_wp, thr=1.0e-9_wp)
      IF (ALLOCATED(error)) THEN
         CALL test_failed(error, &
            & "Complex exchange violates Hermiticity, covariance, or variationality")
      END IF
   END SUBROUTINE test_kpoint_fock_complex

! **************************************************************************************************
!> \brief ...
!> \param error ...
! **************************************************************************************************
   SUBROUTINE test_regular_separable_dft(error)
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      INTEGER, PARAMETER :: ncase = 3, &
         meshes(3, ncase) = RESHAPE([7, 1, 1, 3, 4, 1, 2, 3, 4], [3, ncase]), nrow = 5
      REAL(wp), PARAMETER :: twist(3) = [0.071_wp, -0.113_wp, 0.137_wp]

      COMPLEX(wp)                                        :: lhs, phase, rhs
      COMPLEX(wp), ALLOCATABLE :: dense_k(:, :), dense_r(:, :), mesh_k(:, :), mesh_r(:, :), &
         phase_roots(:, :), sep_k(:, :), sep_r(:, :), twist_phase(:)
      INTEGER                                            :: g(3), icase, icell, idim, igrid, ik, j, &
                                                            nk, root
      INTEGER, ALLOCATABLE                               :: input_to_grid(:), reps(:, :)
      REAL(wp)                                           :: angle, err

      err = 0.0_wp
      DO icase = 1, ncase
         nk = PRODUCT(meshes(:, icase))
         ALLOCATE (input_to_grid(nk), reps(3, nk))
         ALLOCATE (phase_roots(MAXVAL(meshes(:, icase)), 3), twist_phase(nk), &
            & dense_k(nrow, nk), dense_r(nrow, nk), mesh_k(nrow, nk), &
            & mesh_r(nrow, nk), sep_k(nrow, nk), sep_r(nrow, nk), &
            & source=(0.0_wp, 0.0_wp))

         ! A cyclically shifted k ordering and reversed representative ordering
         ! exercise both independent permutations of the compact plan.
         DO ik = 1, nk
            input_to_grid(ik) = 1 + MODULO(ik + 1, nk)
         END DO
         DO icell = 1, nk
            igrid = nk - icell
            reps(1, icell) = MODULO(igrid, meshes(1, icase))
            reps(2, icell) = MODULO(igrid/meshes(1, icase), &
               & meshes(2, icase))
            reps(3, icell) = igrid/(meshes(1, icase)*meshes(2, icase))
            angle = 2.0_wp*pi*DOT_PRODUCT(twist, &
               & REAL(reps(:, icell), wp))
            twist_phase(icell) = EXP(CMPLX(0.0_wp, angle, wp))
         END DO
         phase_roots = (1.0_wp, 0.0_wp)
         DO idim = 1, 3
            DO root = 0, meshes(idim, icase) - 1
               angle = 2.0_wp*pi*REAL(root, wp) &
                  & /REAL(meshes(idim, icase), wp)
               phase_roots(root + 1, idim) = EXP(CMPLX(0.0_wp, angle, wp))
            END DO
         END DO
         DO ik = 1, nk
            DO j = 1, nrow
               mesh_k(j, ik) = CMPLX( &
                  & SIN(REAL(7*j + 11*ik + 3*icase, wp)), &
                  & COS(REAL(13*j + 5*ik + 2*icase, wp)), wp)
               mesh_r(j, ik) = CMPLX( &
                  & COS(REAL(3*j + 17*ik + icase, wp)), &
                  & SIN(REAL(19*j + 2*ik + 4*icase, wp)), wp)
            END DO
         END DO

         dense_r = (0.0_wp, 0.0_wp)
         dense_k = (0.0_wp, 0.0_wp)
         DO ik = 1, nk
            igrid = input_to_grid(ik) - 1
            g(1) = MODULO(igrid, meshes(1, icase))
            g(2) = MODULO(igrid/meshes(1, icase), meshes(2, icase))
            g(3) = igrid/(meshes(1, icase)*meshes(2, icase))
            DO icell = 1, nk
               angle = 2.0_wp*pi*DOT_PRODUCT(twist + REAL(g, wp) &
                  & /REAL(meshes(:, icase), wp), REAL(reps(:, icell), wp))
               phase = EXP(CMPLX(0.0_wp, angle, wp))
               dense_r(:, icell) = dense_r(:, icell) &
                  & + CONJG(phase)*mesh_k(:, ik)/REAL(nk, wp)
               dense_k(:, ik) = dense_k(:, ik) + phase*mesh_r(:, icell)
            END DO
         END DO
         CALL bvk_separable_k_to_r(mesh_k, sep_r, meshes(:, icase), &
            & input_to_grid, reps, phase_roots, twist_phase, .TRUE.)
         CALL bvk_separable_r_to_k(mesh_r, sep_k, meshes(:, icase), &
            & input_to_grid, reps, phase_roots, twist_phase)
         err = MAX(err, MAXVAL(ABS(sep_r - dense_r)), &
            & MAXVAL(ABS(sep_k - dense_k)))

         ! Exact complex adjoint identity for the normalized inverse transform.
         lhs = SUM(CONJG(mesh_r)*sep_r)
         CALL bvk_separable_r_to_k(mesh_r, sep_k, meshes(:, icase), &
            & input_to_grid, reps, phase_roots, twist_phase)
         rhs = SUM(CONJG(sep_k/REAL(nk, wp))*mesh_k)
         err = MAX(err, ABS(lhs - rhs))

         DEALLOCATE (input_to_grid, reps, phase_roots, twist_phase, dense_k, &
            & dense_r, mesh_k, mesh_r, sep_k, sep_r)
      END DO

      CALL check(error, err, 0.0_wp, thr=1.0e-11_wp)
      IF (ALLOCATED(error)) CALL test_failed(error, &
         & "Separable regular-grid DFT differs from the dense character oracle")
   END SUBROUTINE test_regular_separable_dft

! **************************************************************************************************
!> \brief ...
!> \param error ...
! **************************************************************************************************
   SUBROUTINE test_regular_mixed_radix_fft(error)
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      INTEGER, PARAMETER :: ncase = 12, meshes(3, ncase) = RESHAPE([1, 1, 1, 2, 1, 1, 4, 1, 1, 8, 1&
         , 1, 6, 1, 1, 10, 1, 1, 12, 1, 1, 5, 1, 1, 7, 1, 1, 11, 1, 1, 3, 4, 1, 2, 3, 4], [3, ncase&
         ]), nrow = 5
      REAL(wp), PARAMETER :: twist(3) = [0.071_wp, -0.113_wp, 0.137_wp]

      COMPLEX(wp)                                        :: lhs, phase, rhs
      COMPLEX(wp), ALLOCATABLE :: dense_k(:, :), dense_r(:, :), fft_k(:, :), fft_r(:, :), &
         fft_r_unscaled(:, :), mesh_k(:, :), mesh_r(:, :), phase_roots(:, :), roundtrip_k(:, :), &
         twist_phase(:)
      INTEGER                                            :: g(3), icase, icell, idim, igrid, ik, j, &
                                                            nk, root
      INTEGER, ALLOCATABLE                               :: input_to_grid(:), reps(:, :)
      REAL(wp)                                           :: adjoint_err, angle, oracle_err, &
                                                            real_adjoint_err, roundtrip_err

      adjoint_err = 0.0_wp
      oracle_err = 0.0_wp
      real_adjoint_err = 0.0_wp
      roundtrip_err = 0.0_wp
      DO icase = 1, ncase
         nk = PRODUCT(meshes(:, icase))
         ALLOCATE (input_to_grid(nk), reps(3, nk))
         ALLOCATE (phase_roots(MAXVAL(meshes(:, icase)), 3), twist_phase(nk), &
            & dense_k(nrow, nk), dense_r(nrow, nk), fft_k(nrow, nk), &
            & fft_r(nrow, nk), fft_r_unscaled(nrow, nk), mesh_k(nrow, nk), &
            & mesh_r(nrow, nk), roundtrip_k(nrow, nk), &
            & source=(0.0_wp, 0.0_wp))

         DO ik = 1, nk
            input_to_grid(ik) = 1 + MODULO(ik + 1, nk)
         END DO
         DO icell = 1, nk
            igrid = nk - icell
            reps(1, icell) = MODULO(igrid, meshes(1, icase))
            reps(2, icell) = MODULO(igrid/meshes(1, icase), &
               & meshes(2, icase))
            reps(3, icell) = igrid/(meshes(1, icase)*meshes(2, icase))
            IF (MODULO(icell, 2) == 0) THEN
               reps(1, icell) = reps(1, icell) + meshes(1, icase)
            END IF
            IF (MODULO(icell, 3) == 0) THEN
               reps(2, icell) = reps(2, icell) - meshes(2, icase)
            END IF
            IF (MODULO(icell, 5) == 0) THEN
               reps(3, icell) = reps(3, icell) + meshes(3, icase)
            END IF
            angle = 2.0_wp*pi*DOT_PRODUCT(twist, &
               & REAL(reps(:, icell), wp))
            twist_phase(icell) = EXP(CMPLX(0.0_wp, angle, wp))
         END DO
         phase_roots = (1.0_wp, 0.0_wp)
         DO idim = 1, 3
            DO root = 0, meshes(idim, icase) - 1
               angle = 2.0_wp*pi*REAL(root, wp) &
                  & /REAL(meshes(idim, icase), wp)
               phase_roots(root + 1, idim) = EXP(CMPLX(0.0_wp, angle, wp))
            END DO
         END DO
         DO ik = 1, nk
            DO j = 1, nrow
               mesh_k(j, ik) = CMPLX( &
                  & SIN(REAL(7*j + 11*ik + 3*icase, wp)), &
                  & COS(REAL(13*j + 5*ik + 2*icase, wp)), wp)
               mesh_r(j, ik) = CMPLX( &
                  & COS(REAL(3*j + 17*ik + icase, wp)), &
                  & SIN(REAL(19*j + 2*ik + 4*icase, wp)), wp)
            END DO
         END DO

         dense_r = (0.0_wp, 0.0_wp)
         dense_k = (0.0_wp, 0.0_wp)
         DO ik = 1, nk
            igrid = input_to_grid(ik) - 1
            g(1) = MODULO(igrid, meshes(1, icase))
            g(2) = MODULO(igrid/meshes(1, icase), meshes(2, icase))
            g(3) = igrid/(meshes(1, icase)*meshes(2, icase))
            DO icell = 1, nk
               angle = 2.0_wp*pi*DOT_PRODUCT(twist + REAL(g, wp) &
                  & /REAL(meshes(:, icase), wp), REAL(reps(:, icell), wp))
               phase = EXP(CMPLX(0.0_wp, angle, wp))
               dense_r(:, icell) = dense_r(:, icell) &
                  & + CONJG(phase)*mesh_k(:, ik)/REAL(nk, wp)
               dense_k(:, ik) = dense_k(:, ik) + phase*mesh_r(:, icell)
            END DO
         END DO

         CALL bvk_mixed_radix_k_to_r(mesh_k, fft_r, meshes(:, icase), &
            & input_to_grid, reps, phase_roots, twist_phase, .TRUE.)
         CALL bvk_mixed_radix_r_to_k(mesh_r, fft_k, meshes(:, icase), &
            & input_to_grid, reps, phase_roots, twist_phase)
         CALL bvk_mixed_radix_k_to_r(mesh_k, fft_r_unscaled, &
            & meshes(:, icase), input_to_grid, reps, phase_roots, twist_phase, &
            & .FALSE.)
         oracle_err = MAX(oracle_err, MAXVAL(ABS(fft_r - dense_r)), &
            & MAXVAL(ABS(fft_k - dense_k)), &
            & MAXVAL(ABS(fft_r_unscaled - REAL(nk, wp)*dense_r)))

         CALL bvk_mixed_radix_r_to_k(fft_r, roundtrip_k, meshes(:, icase), &
            & input_to_grid, reps, phase_roots, twist_phase)
         roundtrip_err = MAX(roundtrip_err, MAXVAL(ABS(roundtrip_k - mesh_k)))

         ! The unnormalised k-to-R transform is the exact adjoint of R-to-k.
         lhs = SUM(CONJG(mesh_r)*fft_r_unscaled)
         rhs = SUM(CONJG(fft_k)*mesh_k)
         adjoint_err = MAX(adjoint_err, ABS(lhs - rhs))
         real_adjoint_err = MAX(real_adjoint_err, &
            & ABS(REAL(lhs, wp) - REAL(rhs, wp)))

         DEALLOCATE (input_to_grid, reps, phase_roots, twist_phase, dense_k, &
            & dense_r, fft_k, fft_r, fft_r_unscaled, mesh_k, mesh_r, roundtrip_k)
      END DO

      PRINT *, "MIXED_RADIX_FFT_TRANSFORM_MAX_RESID", oracle_err, &
         & roundtrip_err, adjoint_err, real_adjoint_err
      CALL check(error, MAX(oracle_err, roundtrip_err, adjoint_err, &
         & real_adjoint_err), 0.0_wp, thr=1.0e-11_wp)
      IF (ALLOCATED(error)) CALL test_failed(error, &
         & "Mixed-radix FFT differs from the dense character oracle")
   END SUBROUTINE test_regular_mixed_radix_fft

END MODULE tblite_exchange_integration_tests
#endif

PROGRAM tblite_exchange_integration_unittest
#if defined(__TBLITE_GXTB_ACCELERATION)
   USE, INTRINSIC :: iso_fortran_env, ONLY: output_unit
   USE mctc_env_testing, ONLY: run_testsuite
   USE tblite_exchange_integration_tests, ONLY: collect_exchange_integration
#endif
   IMPLICIT NONE
#if defined(__TBLITE_GXTB_ACCELERATION)
   INTEGER :: failures
   failures = 0
   CALL run_testsuite(collect_exchange_integration, output_unit, failures)
   IF (failures /= 0) error STOP "CP2K exchange integration tests failed"
   PRINT '(a)', 'CP2K exchange integration tests passed'
#else
   PRINT '(a)', 'CP2K exchange integration tests skipped: g-xTB response API unavailable'
#endif
END PROGRAM tblite_exchange_integration_unittest
