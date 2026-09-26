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
module tblite_exchange_integration_tests
   use tblite_exchange_contractions, only: get_mulliken_derivs_direct, get_bocorr_derivs_direct, &
      & get_KFock_kmesh, get_KGrad_kmesh, get_KFock_kpoint, get_KGrad_kpoint
   use tblite_stream_plan, only: cp2k_stream_plan_type
   use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
   use, intrinsic :: iso_fortran_env, only : int64
   use mctc_data_vdwrad, only : get_vdw_rad
   use mctc_env, only : wp
   use mctc_io_convert, only : autoaa
   use mctc_io_constants, only : pi
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type, new
   use mctc_ncoord, only : new_ncoord, ncoord_type, cn_count
   use multicharge, only : get_eeqbc_charges
   use tblite_basis_cache, only : basis_cache
   use tblite_basis_qvszp, only : qvszp_basis_type, qvszp_cgto_type, &
      & new_qvszp_cgto, new_qvszp_basis
   use tblite_basis_type, only : basis_type, cgto_container, cgto_type
   use tblite_blas, only : gemv
   use tblite_container_cache, only : container_cache
   use tblite_exchange_stream, only : cp2k_exchange_kmesh, &
      & cp2k_exchange_kmesh_gradient, cp2k_exchange_kpoint_gradient, &
      & cp2k_prepare_exchange, cp2k_exchange_stream_type, &
      & cp2k_exchange_stream_begin, cp2k_exchange_stream_push, &
      & cp2k_exchange_stream_apply, cp2k_exchange_stream_get, &
      & cp2k_exchange_stream_reverse_apply, &
      & cp2k_exchange_stream_reverse_get, &
      & cp2k_exchange_stream_reverse_begin, &
      & cp2k_exchange_stream_reverse_push, &
      & cp2k_exchange_stream_reverse_batch_apply, &
      & cp2k_exchange_stream_reverse_batch_pull, &
      & cp2k_exchange_stream_reverse_batch_advance, &
      & cp2k_exchange_stream_reverse_result, &
      & cp2k_exchange_stream_end, cp2k_exchange_stream_has_full_mesh_storage, &
      & cp2k_exchange_stream_batch_bounds, &
      & cp2k_exchange_stream_batch_apply, &
      & cp2k_exchange_stream_batch_pull, &
      & cp2k_exchange_stream_batch_advance, &
      & cp2k_exchange_stream_transform_advance, &
      & cp2k_exchange_stream_batch_result, &
      & cp2k_exchange_stream_reverse_transform_advance, &
      & cp2k_exchange_stream_peak_complex_elements, cp2k_exchange_stream_pair_tile_size, &
      & cp2k_exchange_stream_image_batch_size, &
      & cp2k_exchange_stream_reduced, cp2k_exchange_stream_oracle, &
      & cp2k_exchange_stream_batched, cp2k_exchange_transform_dense, &
      & cp2k_exchange_transform_separable_dft, &
      & cp2k_exchange_transform_mixed_radix_fft, &
      & cp2k_exchange_partial_descriptor_type, &
      & cp2k_exchange_partial_forward_type, &
      & cp2k_exchange_partial_reverse_type, &
      & cp2k_exchange_partial_begin, cp2k_exchange_partial_push, &
      & cp2k_exchange_partial_export_forward, &
      & cp2k_exchange_partial_import_forward, &
      & cp2k_exchange_partial_forward_size, &
      & cp2k_exchange_partial_forward_pack, &
      & cp2k_exchange_partial_forward_pack_range, &
      & cp2k_exchange_partial_forward_unpack, &
      & cp2k_exchange_partial_forward_unpack_range, &
      & cp2k_exchange_partial_import_forward_range, &
      & cp2k_exchange_partial_reverse_push, &
      & cp2k_exchange_partial_export_reverse, &
      & cp2k_exchange_partial_import_reverse, &
      & cp2k_exchange_partial_reverse_size, &
      & cp2k_exchange_partial_reverse_pack, &
      & cp2k_exchange_partial_reverse_pack_range, &
      & cp2k_exchange_partial_reverse_unpack, &
      & cp2k_exchange_partial_reverse_unpack_range, &
      & cp2k_exchange_partial_import_reverse_range, &
      & cp2k_exchange_partial_descriptor_size, &
      & cp2k_exchange_partial_descriptor_pack, &
      & cp2k_exchange_partial_descriptor_unpack, &
      & cp2k_exchange_partial_descriptors_match, &
      & cp2k_exchange_partial_set_image_range, &
      & cp2k_exchange_partial_set_pull_ownership, &
      & cp2k_exchange_partial_contributor_advance, &
      & cp2k_exchange_partial_discard
   use tblite_cutoff, only : get_lattice_points
   use tblite_data_onecxints, only : get_onecxints
   use tblite_exchange, only : exchange_type, new_exchange_fock, exchange_fock, exchange_cache
   use tblite_exchange_fock, only : exchange_bvk_kernel
   use tblite_integral_native_integrals, only : get_overlap
   use tblite_bvk_transform, only : bvk_separable_k_to_r, &
      & bvk_separable_r_to_k, bvk_mixed_radix_k_to_r, &
      & bvk_mixed_radix_r_to_k
   use tblite_scf, only: new_potential, potential_type
   use tblite_utils_average, only : average_type, new_average, average_id
   use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction
   use tblite_xtb_calculator, only : xtb_calculator
   implicit none
   private

   public :: collect_exchange_integration

   real(wp), parameter :: cutoff = 25.0_wp
   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr1 = 1e5*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))
   integer(int64) :: partial_test_transaction_serial = 0_int64


contains

subroutine collect_exchange_integration(testsuite)
   type(unittest_type), allocatable, intent(out) :: testsuite(:)
   testsuite = [ &
      new_unittest("bvk-exchange-supercell", test_bvk_exchange_kernel), &
      new_unittest("kpoint-fock-complex", test_kpoint_fock_complex), &
      new_unittest("regular-separable-dft", test_regular_separable_dft), &
      new_unittest("regular-mixed-radix-fft", test_regular_mixed_radix_fft) ]
end subroutine collect_exchange_integration

subroutine make_qvszp_basis(bas, mol, error)
   type(qvszp_basis_type), intent(out) :: bas
   type(structure_type), intent(in) :: mol
   type(error_type), allocatable, intent(out) :: error

   !> Number of shells selected from the q-vSZP basis set
   integer, parameter :: p_nshell(60) = [ &
      & 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 2, 3, & !1-20
      & 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, & !21-40
      & 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 2, 3, 3, 4, 4, 4] !41-60

   !> Scaling of the atomic charge in q-vSZP
   real(wp), parameter :: p_qvszp_kq(60) = [&
      &  0.3095916809_wp,  1.3956973435_wp,  0.6565062866_wp,  0.6414725985_wp, & !1-4
      &  0.2289797408_wp,  0.3364428754_wp,  0.3134961856_wp,  0.4760878318_wp, & !5-8
      &  0.6924639973_wp,  1.1737897030_wp,  0.8349448993_wp,  0.7162960964_wp, & !9-12
      &  1.4698576331_wp,  0.4013087493_wp,  0.3606081146_wp,  0.3831440321_wp, & !13-16
      &  0.4730080776_wp,  0.6098840196_wp,  0.2870475794_wp,  1.2095850805_wp, & !17-20
      &  0.5657524714_wp,  0.2998857159_wp,  0.4606086678_wp,  1.0397549071_wp, & !21-24
      &  0.4977801340_wp,  0.3450018637_wp,  0.6470970913_wp,  0.0538554545_wp, & !25-28
      &  0.3313819444_wp,  0.1887280360_wp,  0.3578752499_wp,  0.0778016361_wp, & !29-32
      &  0.4943243573_wp,  0.3404520755_wp,  0.4308339170_wp,  0.6103851192_wp, & !33-36
      &  0.6281945870_wp,  0.4921139609_wp,  0.4926232757_wp,  0.0321070159_wp, & !37-40
      & -0.0305531270_wp,  0.6988874741_wp,  0.2360577501_wp,  0.3463868970_wp, & !41-44
      &  0.5063045934_wp,  0.1909760796_wp,  0.0011218576_wp,  0.7665007573_wp, & !45-48
      &  0.8084270695_wp,  0.7971344741_wp,  0.0628724514_wp,  0.3587445765_wp, & !49-52
      &  0.4782877821_wp,  0.3372686988_wp,  0.2329265134_wp,  0.5416188137_wp, & !53-56
      &  0.1911360071_wp,  0.2615471466_wp,  0.4000000000_wp,  0.4000000000_wp] !57-60

   integer :: isp, izp, ish
   integer, allocatable :: nshell(:)
   type(cgto_container), allocatable :: cgto(:, :)
   type(cgto_container), allocatable :: cgto_h0(:, :)
   type(qvszp_cgto_type), allocatable :: cgto_qvszp

   nshell = p_nshell(mol%num)
   allocate(cgto(maxval(nshell), mol%nid))
   do isp = 1, mol%nid
      izp = mol%num(isp)
      do ish = 1, nshell(isp)
         allocate(cgto_qvszp)
         call new_qvszp_cgto(cgto_qvszp, izp, ish, .true., error, &
            & kq=p_qvszp_kq(izp))
         if (allocated(error)) return
         call move_alloc(cgto_qvszp, cgto(ish, isp)%raw)
      end do
   end do

   call new_qvszp_basis(bas, mol, nshell, cgto, error, accuracy=1.0_wp, &
      & cgto_h0=cgto_h0)
   if (allocated(error)) return

end subroutine make_qvszp_basis

subroutine make_exchange_gxtb(exchange, mol, bas)

   !> New exchange object
   class(exchange_type), allocatable, intent(out) :: exchange

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Basis set information
   class(basis_type), intent(in) :: bas

   !> Atomic Hubbard parameters or chemical hardness
   real(wp), parameter :: hubbard_parameter(60) = [&
      &  0.4725928800_wp,  0.9220339100_wp,  0.1745288800_wp,  0.2570073300_wp, & !1-4
      &  0.3394908600_wp,  0.4219541200_wp,  0.5043819300_wp,  0.5869186300_wp, & !5-8
      &  0.6693135100_wp,  0.7519160700_wp,  0.1796410500_wp,  0.2215727600_wp, & !9-12
      &  0.2634857800_wp,  0.3053964500_wp,  0.3473401400_wp,  0.3892472500_wp, & !13-16
      &  0.4311567000_wp,  0.4730826900_wp,  0.1710546900_wp,  0.2027624400_wp, & !17-20
      &  0.2100732200_wp,  0.2173964700_wp,  0.2247103900_wp,  0.2320150100_wp, & !21-24
      &  0.2393396900_wp,  0.2466563800_wp,  0.2539825500_wp,  0.2612886300_wp, & !25-28
      &  0.2685947600_wp,  0.2759256500_wp,  0.3076299900_wp,  0.3393158000_wp, & !29-32
      &  0.3723598500_wp,  0.4027354900_wp,  0.4344577600_wp,  0.4661170800_wp, & !33-36
      &  0.1558507900_wp,  0.1864932400_wp,  0.1935621000_wp,  0.2006331100_wp, & !37-40
      &  0.2077052200_wp,  0.2147725400_wp,  0.2218461400_wp,  0.2289187200_wp, & !41-44
      &  0.2359862100_wp,  0.2430561200_wp,  0.2501301800_wp,  0.2571993700_wp, & !45-48
      &  0.2878478000_wp,  0.3184867300_wp,  0.3491243100_wp,  0.3797659300_wp, & !49-52
      &  0.4104080800_wp,  0.4410577700_wp,  0.0501933200_wp,  0.0676257000_wp, & !53-56
      &  0.0850444500_wp,  0.1024773600_wp,  0.1199110500_wp,  0.1373277200_wp] !57-60

   real(wp), parameter :: fock_shell_hubbard(4, 60) = reshape([&
      &  3.0691283251_wp,  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !1
      &  5.6874061138_wp,  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !2
      &  2.9883012309_wp,  2.3353110277_wp,  0.0000000000_wp,  0.0000000000_wp, & !3
      &  2.7369793591_wp,  4.6049991892_wp,  0.0000000000_wp,  0.0000000000_wp, & !4
      &  2.3264777208_wp,  3.0982127996_wp,  0.0000000000_wp,  0.0000000000_wp, & !5
      &  2.7744401520_wp,  2.6734099212_wp,  0.0000000000_wp,  0.0000000000_wp, & !6
      &  4.9809093168_wp,  2.5461337307_wp,  0.0000000000_wp,  0.0000000000_wp, & !7
      &  5.8776683336_wp,  1.8560619141_wp,  0.0000000000_wp,  0.0000000000_wp, & !8
      &  5.7439606016_wp,  1.8335282134_wp,  0.0000000000_wp,  0.0000000000_wp, & !9
      &  4.5041390432_wp,  2.2505127692_wp,  0.0000000000_wp,  0.0000000000_wp, & !10
      &  5.1848157487_wp,  3.3581988954_wp,  0.0000000000_wp,  0.0000000000_wp, & !11
      &  2.5240458726_wp,  3.3034249243_wp,  0.0500000000_wp,  0.0000000000_wp, & !12
      &  2.8801279040_wp,  2.9677545399_wp,  0.0532061577_wp,  0.0000000000_wp, & !13
      &  4.7017489845_wp,  3.2451166892_wp,  0.0281128786_wp,  0.0000000000_wp, & !14
      &  5.5892206172_wp,  2.7892566588_wp,  0.0400000000_wp,  0.0000000000_wp, & !15
      &  4.0689142973_wp,  2.3586924706_wp,  0.0190557157_wp,  0.0000000000_wp, & !16
      &  4.0405582899_wp,  2.7444650014_wp,  0.1372215230_wp,  0.0000000000_wp, & !17
      &  3.8099727932_wp,  1.9948057658_wp,  0.3816448105_wp,  0.0000000000_wp, & !18
      &  2.8076441296_wp,  7.3193704755_wp,  0.0000000000_wp,  0.0000000000_wp, & !19
      &  1.4169946183_wp,  2.8506470516_wp,  0.3652352077_wp,  0.0000000000_wp, & !20
      &  3.5613007676_wp,  0.2140174980_wp,  6.0777797114_wp,  0.0000000000_wp, & !21
      &  4.3359867599_wp,  0.4377995490_wp,  5.6339417301_wp,  0.0000000000_wp, & !22
      &  5.1611082663_wp,  0.1493426752_wp,  5.6817801755_wp,  0.0000000000_wp, & !23
      &  3.6190204568_wp,  0.3459252864_wp,  5.1104407513_wp,  0.0000000000_wp, & !24
      &  4.8952268471_wp,  0.4027261153_wp,  5.4351604831_wp,  0.0000000000_wp, & !25
      &  2.9615379891_wp,  0.9965207694_wp,  5.7517572774_wp,  0.0000000000_wp, & !26
      &  2.6929830338_wp,  2.0126619053_wp,  5.3053752298_wp,  0.0000000000_wp, & !27
      &  1.6560700161_wp,  0.6816562078_wp,  5.7265054452_wp,  0.0000000000_wp, & !28
      &  2.6685365919_wp,  0.2648431263_wp,  5.9083571181_wp,  0.0000000000_wp, & !29
      &  3.2723488550_wp,  0.2082148600_wp,  0.0000000000_wp,  0.0000000000_wp, & !30
      &  3.9311761748_wp,  3.5498970279_wp,  0.0500000000_wp,  0.0000000000_wp, & !31
      &  4.0015034453_wp,  2.8636412859_wp,  0.0500000000_wp,  0.0000000000_wp, & !32
      &  3.9383061872_wp,  2.2782579906_wp,  0.0500000000_wp,  0.0000000000_wp, & !33
      &  4.1143033730_wp,  2.0990937392_wp,  0.0500000000_wp,  0.0000000000_wp, & !34
      &  2.6787903609_wp,  2.3185760780_wp,  0.0509775437_wp,  0.0000000000_wp, & !35
      &  3.9635929046_wp,  1.9680253501_wp,  0.0641872392_wp,  0.0000000000_wp, & !36
      &  3.6554521685_wp,  5.4218452070_wp,  0.0000000000_wp,  0.0000000000_wp, & !37
      &  1.2386896827_wp,  0.0367532240_wp,  0.5849934473_wp,  0.0000000000_wp, & !38
      &  2.2661820105_wp,  0.3529903506_wp,  5.7538524045_wp,  0.0000000000_wp, & !39
      &  2.7125214012_wp,  0.9174413894_wp,  5.2357797226_wp,  0.0000000000_wp, & !40
      &  4.2181708986_wp,  0.2363754400_wp,  5.6694352360_wp,  0.0000000000_wp, & !41
      &  5.1555524951_wp,  0.3545348751_wp,  4.1082662773_wp,  0.0000000000_wp, & !42
      &  3.1512449997_wp,  0.1369894568_wp,  4.3933680398_wp,  0.0000000000_wp, & !43
      &  4.4480817529_wp,  1.2580149069_wp,  5.4295548515_wp,  0.0000000000_wp, & !44
      &  4.2971074517_wp,  1.0733384604_wp,  4.4840736158_wp,  0.0000000000_wp, & !45
      &  1.6735108374_wp,  0.8251776679_wp,  3.8684970093_wp,  0.0000000000_wp, & !46
      &  4.2280578634_wp,  0.1699440689_wp,  4.4140710984_wp,  0.0000000000_wp, & !47
      &  3.4063345215_wp,  0.1136195028_wp,  0.0000000000_wp,  0.0000000000_wp, & !48
      &  4.4529575665_wp,  3.5071016244_wp,  0.0865944553_wp,  0.0000000000_wp, & !49
      &  4.0663732807_wp,  2.6287679185_wp,  0.1582429583_wp,  0.0000000000_wp, & !50
      &  4.8767168113_wp,  2.7385818502_wp,  0.0444124760_wp,  0.0000000000_wp, & !51
      &  3.0994063304_wp,  1.8384995197_wp,  0.0457066802_wp,  0.0000000000_wp, & !52
      &  3.2114763510_wp,  2.2737284734_wp,  0.0600000000_wp,  0.0000000000_wp, & !53
      &  3.5544483922_wp,  1.5687665113_wp, -0.0276663624_wp,  0.0000000000_wp, & !54
      &  3.9806565428_wp, 12.1869253675_wp,  0.0000000000_wp,  0.0000000000_wp, & !55
      &  1.8200389980_wp,  0.7112171196_wp,  0.9568227551_wp,  0.0000000000_wp, & !56
      &  5.2212490716_wp,  0.0615487319_wp,  7.8847923351_wp,  0.0000000000_wp, & !57
      &  3.8139401615_wp,  0.4064116095_wp,  3.5473487464_wp,  9.5959434233_wp, & !58
      &  4.1859723893_wp,  0.3661069834_wp,  2.2803077637_wp,  8.6445128775_wp, & !59
      &  3.8143885126_wp,  1.6445614456_wp,  2.2972917570_wp, 10.6325961873_wp],& !60
      & shape(fock_shell_hubbard))

   real(wp), parameter :: fock_avg_exp(4, 60) = reshape([&
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !1
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !2
      &  1.0000000000_wp,  2.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !3
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !4
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !5
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !6
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !7
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !8
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !9
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !10
      &  1.0000000000_wp,  2.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !11
      &  1.0000000000_wp,  1.0000000000_wp,  2.0000000000_wp,  1.0000000000_wp, & !12
      &  1.0000000000_wp,  1.0000000000_wp,  2.0000000000_wp,  1.0000000000_wp, & !13
      &  1.0000000000_wp,  1.0000000000_wp,  2.0000000000_wp,  1.0000000000_wp, & !14
      &  1.0000000000_wp,  1.0000000000_wp,  2.0000000000_wp,  1.0000000000_wp, & !15
      &  1.0000000000_wp,  1.0000000000_wp,  2.0000000000_wp,  1.0000000000_wp, & !16
      &  1.0000000000_wp,  1.0000000000_wp,  2.0000000000_wp,  1.0000000000_wp, & !17
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !18
      &  1.0000000000_wp,  2.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !19
      &  1.0000000000_wp,  1.0000000000_wp,  2.0000000000_wp,  1.0000000000_wp, & !20
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !21
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !22
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !23
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !24
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !25
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !26
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !27
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !28
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !29
      &  1.0000000000_wp,  2.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !30
      &  1.0000000000_wp,  1.0000000000_wp,  2.0000000000_wp,  1.0000000000_wp, & !31
      &  1.0000000000_wp,  1.0000000000_wp,  2.0000000000_wp,  1.0000000000_wp, & !32
      &  1.0000000000_wp,  1.0000000000_wp,  2.0000000000_wp,  1.0000000000_wp, & !33
      &  1.0000000000_wp,  1.0000000000_wp,  2.0000000000_wp,  1.0000000000_wp, & !34
      &  1.0000000000_wp,  1.0000000000_wp,  2.0000000000_wp,  1.0000000000_wp, & !35
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !36
      &  1.0000000000_wp,  2.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !37
      &  1.0000000000_wp,  1.0000000000_wp,  2.0000000000_wp,  1.0000000000_wp, & !38
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !39
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !40
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !41
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !42
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !43
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !44
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !45
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !46
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !47
      &  1.0000000000_wp,  2.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !48
      &  1.0000000000_wp,  1.0000000000_wp,  2.0000000000_wp,  1.0000000000_wp, & !49
      &  1.0000000000_wp,  1.0000000000_wp,  2.0000000000_wp,  1.0000000000_wp, & !50
      &  1.0000000000_wp,  1.0000000000_wp,  2.0000000000_wp,  1.0000000000_wp, & !51
      &  1.0000000000_wp,  1.0000000000_wp,  2.0000000000_wp,  1.0000000000_wp, & !52
      &  1.0000000000_wp,  1.0000000000_wp,  2.0000000000_wp,  1.0000000000_wp, & !53
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !54
      &  1.0000000000_wp,  2.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !55
      &  1.0000000000_wp,  1.0000000000_wp,  2.0000000000_wp,  1.0000000000_wp, & !56
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !57
      &  1.0000000000_wp,  2.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !58
      &  1.0000000000_wp,  2.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !59
      &  1.0000000000_wp,  2.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp],& !60
      & shape(fock_avg_exp))

   !> Parameter: Correlation radius for bond-order correlation correction to the Fock exchange
   real(wp), parameter :: fock_corr_rad(60) = [&
      &  7.7792324405_wp, 25.0000000000_wp,  8.8096512303_wp,  5.0976031817_wp, & !1-4
      & 13.9264964630_wp, 16.1500544221_wp,  4.5816631552_wp,  4.4969668090_wp, & !5-8
      &  2.6478686663_wp,  1.3077199954_wp,  9.1548600312_wp,  8.8377232712_wp, & !9-12
      &  4.8249049725_wp, 11.0887979395_wp, 39.4736123213_wp, 12.5349113695_wp, & !13-16
      &  9.8350329594_wp,  6.1005201925_wp, 10.9559357706_wp, 11.7243327664_wp, & !17-20
      &  6.6159933775_wp, 12.2871886821_wp,  6.5869513804_wp,  5.6059999728_wp, & !21-24
      &  5.1232611294_wp,  7.1831531782_wp,  8.1760551383_wp,  5.7650967598_wp, & !25-28
      &  6.9475450757_wp,  6.8197944590_wp, 10.9433467644_wp, 12.2246936853_wp, & !29-32
      & 16.5095199438_wp,  8.5848893693_wp,  9.2732288842_wp,  6.3347878861_wp, & !33-36
      & 12.1337981396_wp,  8.4248544787_wp,  9.0869814898_wp,  7.1952756180_wp, & !37-40
      & 12.1038392588_wp, 21.7953794701_wp,  5.0119884301_wp, 15.5776135512_wp, & !41-44
      &  5.0766070904_wp, 15.7440179356_wp,  8.0809694075_wp, 10.5430243362_wp, & !45-48
      & 20.8650121649_wp, 15.9260140254_wp, 19.8076868475_wp, 28.4426758250_wp, & !49-52
      & 11.2492941795_wp,  7.3164432561_wp, 15.0000000000_wp, 11.9865176851_wp, & !53-56
      &  9.5686567782_wp,  2.6257749318_wp,  3.1796915366_wp,  4.9047851530_wp] !57-60
      
   !> Parameter: Scaling factor for bond-order correlation correction to the Fock exchange
   real(wp), parameter :: fock_corr_scale(60) = [&
      &  0.0155118215_wp,  0.0017675886_wp,  0.0034145054_wp,  0.0188792718_wp, & !1-4
      &  0.0110877319_wp,  0.0148694859_wp,  0.0519346235_wp,  0.0207320201_wp, & !5-8
      &  0.0458054593_wp,  0.4683557610_wp,  0.0220314012_wp,  0.0131005859_wp, & !9-12
      &  0.0161186017_wp,  0.0242053787_wp,  0.0174007860_wp,  0.0226318766_wp, & !13-16
      &  0.0098241753_wp,  0.0402121335_wp,  0.0497087513_wp,  0.0015206165_wp, & !17-20
      &  0.0017029584_wp,  0.0030354812_wp,  0.0140465231_wp,  0.0120715398_wp, & !21-24
      &  0.0089442952_wp,  0.0047486806_wp,  0.0064777302_wp,  0.0027572746_wp, & !25-28
      &  0.0043159157_wp,  0.0092918050_wp,  0.0168094876_wp,  0.0146359776_wp, & !29-32
      &  0.0140255958_wp,  0.0248514097_wp,  0.0135186782_wp,  0.0174438731_wp, & !33-36
      &  0.0321525083_wp,  0.0032846137_wp,  0.0109494264_wp,  0.0015762687_wp, & !37-40
      &  0.0018844293_wp,  0.0108340538_wp,  0.0053222684_wp,  0.0015263334_wp, & !41-44
      &  0.0058884094_wp,  0.0075970636_wp,  0.0064299028_wp,  0.0124120357_wp, & !45-48
      &  0.0113170904_wp,  0.0158248636_wp,  0.0195890148_wp,  0.0202387721_wp, & !49-52
      &  0.0132956936_wp,  0.0235755901_wp,  0.0347877982_wp,  0.0048909086_wp, & !53-56
      &  0.0008293725_wp,  0.0085376871_wp,  0.0034042669_wp,  0.0089899347_wp] !57-60

   !> Parameter: van-der-Waals radii scaling
   real(wp), parameter :: rvdw_scale(60) = [&
      &  0.9950000000_wp,  1.1000000000_wp,  1.0000000000_wp,  1.0200000000_wp, & !1-4
      &  1.0200000000_wp,  1.0000000000_wp,  1.0500000000_wp,  1.0200000000_wp, & !5-8
      &  1.1800000000_wp,  0.8500000000_wp,  1.0500000000_wp,  1.1500000000_wp, & !9-12
      &  1.0400000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !13-16
      &  1.0000000000_wp,  1.0000000000_wp,  0.9500000000_wp,  0.9500000000_wp, & !17-20
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !21-24
      &  1.0000000000_wp,  1.1000000000_wp,  1.1000000000_wp,  1.1000000000_wp, & !25-28
      &  1.0000000000_wp,  1.1000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !29-32
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !33-36
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !37-40
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !41-44
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.1500000000_wp, & !45-48
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !49-52
      &  1.0000000000_wp,  1.0000000000_wp,  0.5000000000_wp,  1.0000000000_wp, & !53-56
      &  1.2000000000_wp,  1.2000000000_wp,  1.2000000000_wp,  1.2000000000_wp] !57-60

   real(wp), parameter :: p_kq_shell(0:3) = &
      & [1.0000000000_wp, 0.5000000000_wp, 0.2500000000_wp, 0.1250000000_wp]
   real(wp), parameter :: p_kq = 2.1400000000_wp
   real(wp), parameter :: p_ondiag = 1.4017990616_wp
   real(wp), parameter :: p_offdiag_l(0:3) = &
      & [0.0500000000_wp, 2.5000000000_wp, 2.5000000000_wp, 2.5000000000_wp]
   real(wp), parameter :: p_hubbard_exp = -0.3296069670_wp
   real(wp), parameter :: p_hubbard_exp_r0 = 0.0128349272_wp * autoaa
   real(wp), parameter :: p_omega = 0.2000000000_wp
   real(wp), parameter :: p_a0 = 0.1500000000_wp
   real(wp), parameter :: p_a0_onsite = 0.1820686001_wp
   real(wp), parameter :: p_corr_exp = 1.2000000000_wp / autoaa
   real(wp), parameter :: p_corr_scale_avg_exponent = 0.7000000000_wp

   integer :: izp, jzp, isp, jsp, ish, jsh, il, jl
   real(wp), allocatable :: hardness_fx(:, :), offdiag(:, :)
   real(wp), allocatable :: onecxints(:, :, :), kq(:, :)
   real(wp), allocatable :: avg_exponents(:, :), rad(:, :)
   real(wp), allocatable :: corr_scale(:), corr_rad(:)
   type(exchange_fock), allocatable :: fock
   type(average_type), allocatable :: hubbard_average
   type(average_type), allocatable :: offdiag_average
   type(average_type), allocatable :: corr_scale_average
   type(average_type), allocatable :: corr_rad_average
   type(average_type), allocatable :: rvdw_average
   
   ! Obtain radii for radius-dependent hubbard scaling
   allocate(rvdw_average)
   call new_average(rvdw_average, average_id%arithmetic)
   allocate(rad(mol%nid, mol%nid))
   do isp = 1, mol%nid
      izp = mol%num(isp)
      do jsp = 1, isp
         jzp = mol%num(jsp)
         rad(jsp, isp) = get_vdw_rad(jzp, izp) &
            & * rvdw_average%value(rvdw_scale(izp), rvdw_scale(jzp))
         rad(isp, jsp) = rad(jsp, isp)
      end do
   end do

   ! Obtain the shell Hubbard parameters
   allocate(hardness_fx(maxval(bas%nsh_id), mol%nid))
   hardness_fx(:, :) = 0.0_wp
   do isp = 1, mol%nid
      izp = mol%num(isp)
      do ish = 1, bas%nsh_id(isp)
         hardness_fx(ish, isp) = hubbard_parameter(izp) &
            & * fock_shell_hubbard(ish, izp)
      end do
   end do

   ! Set averager for hubbard and correlation parameters
   allocate(hubbard_average, offdiag_average, corr_scale_average, corr_rad_average)
   call new_average(hubbard_average, average_id%general)
   call new_average(offdiag_average, average_id%geometric)
   call new_average(corr_scale_average, average_id%general)
   call new_average(corr_rad_average, average_id%harmonic)

   ! Obtain the averaging exponents for the shell Hubbard parameters
   allocate(avg_exponents(maxval(bas%nsh_id), mol%nid))
   avg_exponents(:, :) = 0.0_wp
   do isp = 1, mol%nid
      izp = mol%num(isp)
      do ish = 1, bas%nsh_id(isp)
         avg_exponents(ish, isp) = fock_avg_exp(ish, izp)
      end do
   end do

   ! Obtain the one-center exchange integrals
   allocate(onecxints(maxval(bas%nsh_id), maxval(bas%nsh_id), mol%nid))
   onecxints(:, :, :) = 0.0_wp
   do isp = 1, mol%nid
      izp = mol%num(isp)
      do ish = 1, bas%nsh_id(isp)
         il = bas%cgto(ish, isp)%raw%ang
         do jsh = 1, bas%nsh_id(isp)
            jl = bas%cgto(jsh, isp)%raw%ang
            onecxints(jsh, ish, isp) = get_onecxints(jl, il, izp)
         end do
      end do
   end do

   ! Obtain the off-diagaonal scaling exchange factors
   allocate(offdiag(maxval(bas%nsh_id), mol%nid))
   offdiag(:, :) = 0.0_wp
   do isp = 1, mol%nid
      izp = mol%num(isp)
      do ish = 1, bas%nsh_id(isp)
         il = bas%cgto(ish, isp)%raw%ang
         offdiag(ish, isp) = p_offdiag_l(il)
      end do
   end do

   ! Obtain the charge scaling of the exchange fraction
   allocate(kq(maxval(bas%nsh_id), mol%nid))
   kq(:, :) = 0.0_wp
   do isp = 1, mol%nid
      izp = mol%num(isp)
      do ish = 1, bas%nsh_id(isp)
         il = bas%cgto(ish, isp)%raw%ang
         kq(ish, isp) = sqrt(p_kq) * p_kq_shell(il)
      end do
   end do

   ! Obtain the correlation scaling factor and radii
   corr_scale = fock_corr_scale(mol%num)
   corr_rad = fock_corr_rad(mol%num)

   allocate(fock)
   call new_exchange_fock(fock, mol, bas, hardness_fx, hubbard_average, &
      & avg_exponents, p_ondiag, offdiag, offdiag_average, p_hubbard_exp, &
      & p_hubbard_exp_r0, rad, 1.0_wp, onecxints, kq, corr_scale, &
      & corr_scale_average, p_corr_scale_avg_exponent, p_corr_exp, corr_rad, &
      & corr_rad_average, p_a0, p_a0_onsite, p_omega, 1.0_wp-p_a0)
   call move_alloc(fock, exchange)

end subroutine make_exchange_gxtb

subroutine view(cache, ptr)
   !> Instance of the container cache
   type(container_cache), target, intent(inout) :: cache
   !> Reference to the container cache
   type(exchange_cache), pointer, intent(out) :: ptr
   nullify(ptr)
   select type(target => cache%raw)
   type is(exchange_cache)
      ptr => target
   end select
end subroutine view

subroutine next_partial_test_transaction(transaction_id)
   integer(int64), intent(out) :: transaction_id(2)

   partial_test_transaction_serial = partial_test_transaction_serial+1_int64
   transaction_id = [3141592653589793_int64, partial_test_transaction_serial]
end subroutine next_partial_test_transaction

subroutine check_partial_descriptor_roundtrip(descriptor, error, reject_bad_abi)
   type(cp2k_exchange_partial_descriptor_type), intent(in) :: descriptor
   type(error_type), allocatable, intent(out) :: error
   logical, intent(in), optional :: reject_bad_abi

   type(cp2k_exchange_partial_descriptor_type) :: roundtrip
   integer(int64), allocatable :: words(:)

   call cp2k_exchange_partial_descriptor_pack(descriptor, words, error)
   if (allocated(error)) return
   if (int(size(words), int64) /= cp2k_exchange_partial_descriptor_size()) then
      call test_failed(error, "Incorrect canonical partial descriptor size")
      return
   end if
   call cp2k_exchange_partial_descriptor_unpack(roundtrip, words, error)
   if (allocated(error)) return
   if (.not.cp2k_exchange_partial_descriptors_match(descriptor, roundtrip)) then
      call test_failed(error, "Canonical partial descriptor roundtrip changed")
      return
   end if

   if (present(reject_bad_abi)) then
      if (reject_bad_abi) then
         words(1) = words(1)+1_int64
         call cp2k_exchange_partial_descriptor_unpack(roundtrip, words, error)
         if (.not.allocated(error)) then
            call test_failed(error, &
               & "Malformed canonical partial descriptor was accepted")
            return
         end if
         deallocate(error)
      end if
   end if
end subroutine check_partial_descriptor_roundtrip

subroutine check_partial_exchange_case(calc, mol, cache, cache_plan, wfn, nmesh, kfrac, &
   & weights, density, overlap, partition_counts, run_negative, error, maxerr)
   type(xtb_calculator), intent(inout) :: calc
   type(structure_type), intent(in) :: mol
   type(container_cache), intent(inout) :: cache
   type(cp2k_stream_plan_type), intent(inout) :: cache_plan
   type(wavefunction_type), intent(in) :: wfn
   integer, intent(in) :: nmesh(3), partition_counts(:)
   real(wp), intent(in) :: kfrac(:, :), weights(:)
   complex(wp), intent(in) :: density(:, :, :, :), overlap(:, :, :)
   logical, intent(in) :: run_negative
   type(error_type), allocatable, intent(out) :: error
   real(wp), intent(out) :: maxerr

   type(cp2k_exchange_partial_descriptor_type), allocatable :: descriptor(:)
   type(cp2k_exchange_partial_descriptor_type) :: mismatch_descriptor
   type(cp2k_exchange_partial_forward_type), allocatable :: forward_buffer(:)
   type(cp2k_exchange_partial_forward_type) :: forward_sum_buffer, &
      & mismatch_forward_buffer
   type(cp2k_exchange_partial_reverse_type), allocatable :: reverse_buffer(:)
   type(cp2k_exchange_partial_reverse_type) :: reverse_sum_buffer
   type(cp2k_exchange_stream_type), allocatable :: stream(:)
   type(cp2k_exchange_stream_type) :: mismatch_stream
   integer :: ik, ipart, ncomplex, nk, npart
   integer(int64) :: mismatch_transaction_id(2), transaction_id(2)
   integer(int64), allocatable :: descriptor_words(:)
   integer, allocatable :: owner_count(:), owner_sum(:), owner_trial(:)
   integer, allocatable :: owner_matrix(:, :)
   logical :: done
   real(wp) :: energy_partial, energy_ref
   real(wp), allocatable :: gradient_partial(:, :), gradient_ref(:, :), &
      & sigma_partial(:, :), sigma_ref(:, :), vsh_partial(:), vsh_ref(:)
   complex(wp), allocatable :: fock_partial(:, :, :, :), &
      & fock_ref(:, :, :, :), forward_raw(:), forward_sum(:), &
      & overlap_adjoint(:, :, :), overlap_adjoint_ref(:, :, :), &
      & reverse_raw(:), reverse_sum(:), density_trial(:, :, :)

   nk = size(weights)
   allocate(fock_ref(size(density, 1), size(density, 2), &
      & size(density, 3), nk), vsh_ref(size(wfn%qsh, 1)), &
      & overlap_adjoint_ref(size(overlap, 1), size(overlap, 2), nk), &
      & gradient_ref(3, mol%nat), sigma_ref(3, 3))
   call cp2k_exchange_kmesh(calc, mol, cache, cache_plan, wfn, nmesh, kfrac, weights, &
      & density, overlap, fock_ref, vsh_ref, energy_ref, error)
   if (allocated(error)) return
   call cp2k_exchange_kmesh_gradient(calc, mol, cache, cache_plan, nmesh, kfrac, weights, &
      & density, overlap, overlap_adjoint_ref, gradient_ref, sigma_ref, error)
   if (allocated(error)) return

   maxerr = 0.0_wp
   do ipart = 1, size(partition_counts)
      npart = partition_counts(ipart)
      allocate(stream(npart), descriptor(npart), forward_buffer(npart), &
         & owner_matrix(nk, npart))
      owner_matrix = 0
      call next_partial_test_transaction(transaction_id)
      do ik = 1, npart
         call cp2k_exchange_partial_begin(stream(ik), calc, mol, cache, cache_plan, wfn, &
            & nmesh, kfrac, weights, nk, transaction_id, error)
         if (allocated(error)) return
      end do
      do ik = nk, 1, -1
         ncomplex = 1+modulo(ik-1, npart)
         call cp2k_exchange_partial_push(stream(ncomplex), ik, &
            & density(:, :, :, ik), overlap(:, :, ik), error)
         if (allocated(error)) return
      end do

      if (run_negative .and. npart == 2) then
         call cp2k_exchange_partial_push(stream(1), 1, density(:, :, :, 1), &
            & overlap(:, :, 1), error)
         if (.not.allocated(error)) then
            call test_failed(error, "Duplicate partial ownership was accepted")
            return
         end if
         deallocate(error)
         call cp2k_exchange_partial_push(stream(1), nk+1, &
            & density(:, :, :, 1), overlap(:, :, 1), error)
         if (.not.allocated(error)) then
            call test_failed(error, "Out-of-range partial ownership was accepted")
            return
         end if
         deallocate(error)
         call cp2k_exchange_stream_batch_apply(stream(1), calc, mol, cache, cache_plan, &
            & error)
         if (.not.allocated(error)) then
            call test_failed(error, "Partial kernel apply before import succeeded")
            return
         end if
         deallocate(error)
      end if

      do ik = 1, npart
         call cp2k_exchange_partial_export_forward(stream(ik), &
            & descriptor(ik), owner_count, forward_buffer(ik), error)
         if (allocated(error)) return
         call check_partial_descriptor_roundtrip(descriptor(ik), error, &
            & run_negative .and. npart == 2 .and. ik == 1)
         if (allocated(error)) return
         owner_matrix(:, ik) = owner_count
         if (ik > 1) then
            if (.not.cp2k_exchange_partial_descriptors_match( &
               & descriptor(1), descriptor(ik))) then
               call test_failed(error, "Matching partial descriptors differ")
               return
            end if
         end if
         call cp2k_exchange_partial_forward_pack(forward_buffer(ik), &
            & forward_raw, error)
         if (allocated(error)) return
         if (ik == 1) then
            allocate(forward_sum(size(forward_raw)), &
               & source=(0.0_wp, 0.0_wp))
         end if
         forward_sum = forward_sum+forward_raw
      end do
      owner_sum = sum(owner_matrix, dim=2)
      ncomplex = 4*size(density, 1)*size(density, 2)*nk*size(density, 3)
      if (cp2k_exchange_partial_forward_size(forward_buffer(1)) /= &
         & int(ncomplex, int64)) then
         call test_failed(error, "Incorrect partial forward memory accounting")
         return
      end if
      call cp2k_exchange_partial_forward_unpack(forward_sum_buffer, &
         & descriptor(1), forward_sum, error)
      if (allocated(error)) return

      ! Simulate P>Nk: an empty fourth rank must still carry the exact same
      ! descriptor.  A divergent transaction identity is rejected before the
      ! summed payload can enter the live rank-zero accumulator.
      if (run_negative .and. npart == 4) then
         call next_partial_test_transaction(mismatch_transaction_id)
         call cp2k_exchange_partial_begin(mismatch_stream, calc, mol, cache, cache_plan, &
            & wfn, nmesh, kfrac, weights, nk, mismatch_transaction_id, error)
         if (allocated(error)) return
         call cp2k_exchange_partial_export_forward(mismatch_stream, &
            & mismatch_descriptor, owner_count, mismatch_forward_buffer, error)
         if (allocated(error)) return
         if (any(owner_count /= 0)) then
            call test_failed(error, "The simulated empty rank owns k points")
            return
         end if
         call check_partial_descriptor_roundtrip(mismatch_descriptor, error)
         if (allocated(error)) return
         if (cp2k_exchange_partial_descriptors_match(descriptor(1), &
            & mismatch_descriptor)) then
            call test_failed(error, &
               & "Divergent empty-rank partial descriptors compare equal")
            return
         end if
         call cp2k_exchange_partial_import_forward(stream(1), &
            & mismatch_descriptor, owner_sum, forward_sum_buffer, calc, mol, &
            & cache, cache_plan, error)
         if (.not.allocated(error)) then
            call test_failed(error, &
               & "Divergent empty-rank partial descriptor was accepted")
            return
         end if
         deallocate(error)
         call cp2k_exchange_partial_discard(mismatch_stream, error)
         if (allocated(error)) return
      end if

      if (run_negative .and. npart == 2) then
         owner_trial = owner_sum
         owner_trial(1) = 0
         call cp2k_exchange_partial_import_forward(stream(1), descriptor(1), &
            & owner_trial, forward_sum_buffer, calc, mol, cache, cache_plan, error)
         if (.not.allocated(error)) then
            call test_failed(error, "Missing partial ownership was accepted")
            return
         end if
         deallocate(error)
         owner_trial = owner_sum
         owner_trial(1) = 2
         call cp2k_exchange_partial_import_forward(stream(1), descriptor(1), &
            & owner_trial, forward_sum_buffer, calc, mol, cache, cache_plan, error)
         if (.not.allocated(error)) then
            call test_failed(error, "Duplicate partial ownership was accepted")
            return
         end if
         deallocate(error)

         ! Malicious but structurally valid large descriptor dimensions must
         ! fail before any default-integer product or allocation can overflow.
         call cp2k_exchange_partial_descriptor_pack(descriptor(1), &
            & descriptor_words, error)
         if (allocated(error)) return
         descriptor_words(8) = int(huge(0), int64)
         call cp2k_exchange_partial_descriptor_unpack(mismatch_descriptor, &
            & descriptor_words, error)
         if (allocated(error)) return
         call cp2k_exchange_partial_forward_unpack(mismatch_forward_buffer, &
            & mismatch_descriptor, forward_sum, error)
         if (.not.allocated(error)) then
            call test_failed(error, &
               & "Oversized partial payload arithmetic was accepted")
            return
         end if
         deallocate(error)

         ! A fresh logical transaction may have the same model, geometry and
         ! mesh but different density.  Its mandatory host identity prevents
         ! replay of an old descriptor and pre-kernel payload.
         call next_partial_test_transaction(mismatch_transaction_id)
         call cp2k_exchange_partial_begin(mismatch_stream, calc, mol, cache, cache_plan, &
            & wfn, nmesh, kfrac, weights, nk, mismatch_transaction_id, error)
         if (allocated(error)) return
         allocate(density_trial(size(density, 1), size(density, 2), &
            & size(density, 3)))
         do ik = 1, nk
            density_trial = density(:, :, :, ik)
            if (ik == 1) density_trial(1, 1, 1) = &
               & density_trial(1, 1, 1)+cmplx(1.0e-7_wp, 0.0_wp, wp)
            call cp2k_exchange_partial_push(mismatch_stream, ik, &
               & density_trial, overlap(:, :, ik), error)
            if (allocated(error)) return
         end do
         deallocate(density_trial)
         call cp2k_exchange_partial_export_forward(mismatch_stream, &
            & mismatch_descriptor, owner_count, mismatch_forward_buffer, error)
         if (allocated(error)) return
         call check_partial_descriptor_roundtrip(mismatch_descriptor, error)
         if (allocated(error)) return
         if (cp2k_exchange_partial_descriptors_match(descriptor(1), &
            & mismatch_descriptor)) then
            call test_failed(error, &
               & "Fresh transaction descriptors compare equal")
            return
         end if
         call cp2k_exchange_partial_import_forward(mismatch_stream, &
            & descriptor(1), owner_sum, forward_sum_buffer, calc, mol, cache, cache_plan, &
            & error)
         if (.not.allocated(error)) then
            call test_failed(error, &
               & "A cross-transaction partial payload replay was accepted")
            return
         end if
         deallocate(error)
         call cp2k_exchange_partial_discard(mismatch_stream, error)
         if (allocated(error)) return

         call next_partial_test_transaction(mismatch_transaction_id)
         call cp2k_exchange_partial_begin(mismatch_stream, calc, mol, cache, cache_plan, &
            & wfn, nmesh, kfrac, weights, nk-1, mismatch_transaction_id, error, &
            & image_first=1, image_count=nk-1)
         if (allocated(error)) return
         call cp2k_exchange_partial_discard(mismatch_stream, error)
         if (.not.allocated(error)) then
            call test_failed(error, &
               & "An unexported partial accumulator was discarded")
            return
         end if
         deallocate(error)
         call cp2k_exchange_partial_export_forward(mismatch_stream, &
            & mismatch_descriptor, owner_count, forward_buffer(1), error)
         if (allocated(error)) return
         call cp2k_exchange_partial_import_forward(stream(1), &
            & mismatch_descriptor, owner_sum, forward_sum_buffer, calc, mol, &
            & cache, cache_plan, error)
         if (.not.allocated(error)) then
            call test_failed(error, "Mismatched partial descriptor was accepted")
            return
         end if
         deallocate(error)
         call cp2k_exchange_partial_discard(mismatch_stream, error)
         if (allocated(error)) return

         forward_raw = forward_sum
         forward_raw(1) = cmplx(ieee_value(0.0_wp, ieee_quiet_nan), &
            & 0.0_wp, wp)
         call cp2k_exchange_partial_forward_unpack(forward_sum_buffer, &
            & descriptor(1), forward_raw, error)
         if (.not.allocated(error)) then
            call test_failed(error, "Non-finite partial payload was accepted")
            return
         end if
         deallocate(error)
         call cp2k_exchange_partial_forward_unpack(forward_sum_buffer, &
            & descriptor(1), forward_sum, error)
         if (allocated(error)) return
      end if

      call cp2k_exchange_partial_import_forward(stream(1), descriptor(1), &
         & owner_sum, forward_sum_buffer, calc, mol, cache, cache_plan, error)
      if (allocated(error)) return
      if (run_negative .and. npart == 2) then
         call cp2k_exchange_partial_discard(stream(1), error)
         if (.not.allocated(error)) then
            call test_failed(error, &
               & "An imported partial accumulator was discarded")
            return
         end if
         deallocate(error)
         call cp2k_exchange_partial_import_forward(stream(1), descriptor(1), &
            & owner_sum, forward_sum_buffer, calc, mol, cache, cache_plan, error)
         if (.not.allocated(error)) then
            call test_failed(error, "Second partial import was accepted")
            return
         end if
         deallocate(error)
         call cp2k_exchange_partial_push(stream(1), 1, density(:, :, :, 1), &
            & overlap(:, :, 1), error)
         if (.not.allocated(error)) then
            call test_failed(error, "Partial push after import was accepted")
            return
         end if
         deallocate(error)
      end if

      allocate(fock_partial(size(fock_ref, 1), size(fock_ref, 2), &
         & size(fock_ref, 3), nk), source=(0.0_wp, 0.0_wp))
      allocate(vsh_partial(size(vsh_ref)))
      call cp2k_exchange_stream_batch_apply(stream(1), calc, mol, cache, cache_plan, error)
      if (allocated(error)) return
      do ik = 1, nk
         call cp2k_exchange_stream_batch_pull(stream(1), ik, &
            & overlap(:, :, ik), fock_partial(:, :, :, ik), error)
         if (allocated(error)) return
      end do
      call cp2k_exchange_stream_batch_advance(stream(1), done, error)
      if (allocated(error)) return
      if (.not.done) then
         call test_failed(error, "Single full-image partial batch did not finish")
         return
      end if
      call cp2k_exchange_stream_batch_result(stream(1), calc, mol, cache, cache_plan, &
         & vsh_partial, energy_partial, error)
      if (allocated(error)) return
      maxerr = max(maxerr, maxval(abs(fock_partial-fock_ref)), &
         & maxval(abs(vsh_partial-vsh_ref)), abs(energy_partial-energy_ref))
      call cp2k_exchange_stream_end(stream(1), error)
      if (allocated(error)) return
      do ik = 2, npart
         call cp2k_exchange_partial_discard(stream(ik), error)
         if (allocated(error)) return
      end do
      deallocate(stream, descriptor, forward_buffer, owner_matrix, &
         & forward_sum, fock_partial, vsh_partial)

      allocate(stream(npart), descriptor(npart), reverse_buffer(npart), &
         & owner_matrix(nk, npart))
      owner_matrix = 0
      call next_partial_test_transaction(transaction_id)
      do ik = 1, npart
         call cp2k_exchange_partial_begin(stream(ik), calc, mol, cache, cache_plan, wfn, &
            & nmesh, kfrac, weights, nk, transaction_id, error, &
            & reverse_only=.true.)
         if (allocated(error)) return
      end do
      do ik = nk, 1, -1
         ncomplex = 1+modulo(ik-1, npart)
         call cp2k_exchange_partial_reverse_push(stream(ncomplex), ik, &
            & density(:, :, :, ik), overlap(:, :, ik), error)
         if (allocated(error)) return
      end do
      do ik = 1, npart
         call cp2k_exchange_partial_export_reverse(stream(ik), &
            & descriptor(ik), owner_count, reverse_buffer(ik), error)
         if (allocated(error)) return
         call check_partial_descriptor_roundtrip(descriptor(ik), error)
         if (allocated(error)) return
         owner_matrix(:, ik) = owner_count
         if (ik > 1) then
            if (.not.cp2k_exchange_partial_descriptors_match( &
               & descriptor(1), descriptor(ik))) then
               call test_failed(error, "Matching reverse descriptors differ")
               return
            end if
         end if
         call cp2k_exchange_partial_reverse_pack(reverse_buffer(ik), &
            & reverse_raw, error)
         if (allocated(error)) return
         if (ik == 1) then
            allocate(reverse_sum(size(reverse_raw)), &
               & source=(0.0_wp, 0.0_wp))
         end if
         reverse_sum = reverse_sum+reverse_raw
      end do
      owner_sum = sum(owner_matrix, dim=2)
      ncomplex = 6*size(density, 1)*size(density, 2)*nk*size(density, 3) &
         & +5*size(density, 1)*size(density, 3)
      if (cp2k_exchange_partial_reverse_size(reverse_buffer(1)) /= &
         & int(ncomplex, int64)) then
         call test_failed(error, "Incorrect partial reverse memory accounting")
         return
      end if
      call cp2k_exchange_partial_reverse_unpack(reverse_sum_buffer, &
         & descriptor(1), reverse_sum, error)
      if (allocated(error)) return
      call cp2k_exchange_partial_import_reverse(stream(1), descriptor(1), &
         & owner_sum, reverse_sum_buffer, calc, mol, cache, cache_plan, error)
      if (allocated(error)) return

      allocate(overlap_adjoint(size(overlap_adjoint_ref, 1), &
         & size(overlap_adjoint_ref, 2), nk), source=(0.0_wp, 0.0_wp))
      allocate(gradient_partial(3, mol%nat), sigma_partial(3, 3))
      call cp2k_exchange_stream_reverse_batch_apply(stream(1), calc, mol, &
         & cache, cache_plan, error)
      if (allocated(error)) return
      do ik = 1, nk
         call cp2k_exchange_stream_reverse_batch_pull(stream(1), ik, &
            & density(:, :, :, ik), overlap(:, :, ik), &
            & overlap_adjoint(:, :, ik), error)
         if (allocated(error)) return
      end do
      call cp2k_exchange_stream_reverse_batch_advance(stream(1), done, error)
      if (allocated(error)) return
      if (.not.done) then
         call test_failed(error, "Single reverse partial batch did not finish")
         return
      end if
      call cp2k_exchange_stream_reverse_result(stream(1), calc, mol, cache, cache_plan, &
         & gradient_partial, sigma_partial, error)
      if (allocated(error)) return
      maxerr = max(maxerr, &
         & maxval(abs(overlap_adjoint-overlap_adjoint_ref)), &
         & maxval(abs(gradient_partial-gradient_ref)), &
         & maxval(abs(sigma_partial-sigma_ref)))
      call cp2k_exchange_stream_end(stream(1), error)
      if (allocated(error)) return
      do ik = 2, npart
         call cp2k_exchange_partial_discard(stream(ik), error)
         if (allocated(error)) return
      end do
      deallocate(stream, descriptor, reverse_buffer, owner_matrix, reverse_sum, &
         & overlap_adjoint, gradient_partial, sigma_partial)
   end do
end subroutine check_partial_exchange_case

subroutine check_partial_exchange_multibatch(calc, mol, cache, cache_plan, wfn, nmesh, &
   & kfrac, weights, density, overlap, error, maxerr)
   type(xtb_calculator), intent(inout) :: calc
   type(structure_type), intent(in) :: mol
   type(container_cache), intent(inout) :: cache
   type(cp2k_stream_plan_type), intent(inout) :: cache_plan
   type(wavefunction_type), intent(in) :: wfn
   integer, intent(in) :: nmesh(3)
   real(wp), intent(in) :: kfrac(:, :), weights(:)
   complex(wp), intent(in) :: density(:, :, :, :), overlap(:, :, :)
   type(error_type), allocatable, intent(out) :: error
   real(wp), intent(out) :: maxerr

   integer, parameter :: npart = 4
   type(cp2k_exchange_partial_descriptor_type) :: descriptor(npart), &
      & stale_descriptor
   type(cp2k_exchange_partial_forward_type) :: forward_buffer(npart), &
      & forward_sum_buffer, contributor_forward_buffer
   type(cp2k_exchange_partial_reverse_type) :: reverse_buffer(npart), &
      & reverse_sum_buffer, contributor_reverse_buffer
   type(cp2k_exchange_stream_type) :: stream(npart), contributor_stream
   type(cp2k_exchange_partial_descriptor_type) :: contributor_descriptor
   integer :: ibatch, ik, ipart, ncomplex, nk
   integer(int64) :: transaction_id(2)
   integer, allocatable :: owner_count(:), owner_sum(:), owner_matrix(:, :), &
      & pull_owner(:)
   logical :: contributor_done, done
   real(wp) :: energy_partial, energy_ref
   real(wp), allocatable :: gradient_partial(:, :), gradient_ref(:, :), &
      & sigma_partial(:, :), sigma_ref(:, :), vsh_partial(:), vsh_ref(:)
   complex(wp), allocatable :: fock_partial(:, :, :, :), &
      & fock_ref(:, :, :, :), forward_raw(:), forward_sum(:), &
      & overlap_adjoint(:, :, :), overlap_adjoint_ref(:, :, :), &
      & reverse_raw(:), reverse_sum(:)

   nk = size(weights)
   allocate(pull_owner(nk))
   do ik = 1, nk
      pull_owner(ik) = 1+modulo(ik-1, npart)
   end do
   allocate(fock_ref(size(density, 1), size(density, 2), &
      & size(density, 3), nk), vsh_ref(size(wfn%qsh, 1)), &
      & overlap_adjoint_ref(size(overlap, 1), size(overlap, 2), nk), &
      & gradient_ref(3, mol%nat), sigma_ref(3, 3))
   call cp2k_exchange_kmesh(calc, mol, cache, cache_plan, wfn, nmesh, kfrac, weights, &
      & density, overlap, fock_ref, vsh_ref, energy_ref, error)
   if (allocated(error)) return
   call cp2k_exchange_kmesh_gradient(calc, mol, cache, cache_plan, nmesh, kfrac, weights, &
      & density, overlap, overlap_adjoint_ref, gradient_ref, sigma_ref, error)
   if (allocated(error)) return

   allocate(owner_matrix(nk, npart))
   allocate(fock_partial(size(fock_ref, 1), size(fock_ref, 2), &
      & size(fock_ref, 3), nk), source=(0.0_wp, 0.0_wp))
   allocate(vsh_partial(size(vsh_ref)))
   call next_partial_test_transaction(transaction_id)
   do ipart = 1, npart
      call cp2k_exchange_partial_begin(stream(ipart), calc, mol, cache, cache_plan, wfn, &
         & nmesh, kfrac, weights, 1, transaction_id, error)
      if (allocated(error)) return
   end do
   call cp2k_exchange_partial_begin(contributor_stream, calc, mol, cache, cache_plan, wfn, &
      & nmesh, kfrac, weights, 1, transaction_id, error)
   if (allocated(error)) return
   do ibatch = 1, nk
      owner_matrix = 0
      do ik = nk, 1, -1
         ipart = 1+modulo(ik-1, npart)
         call cp2k_exchange_partial_push(stream(ipart), ik, &
            & density(:, :, :, ik), overlap(:, :, ik), error)
         if (allocated(error)) return
      end do
      if (allocated(forward_sum)) deallocate(forward_sum)
      do ipart = 1, npart
         call cp2k_exchange_partial_export_forward(stream(ipart), &
            & descriptor(ipart), owner_count, forward_buffer(ipart), error)
         if (allocated(error)) return
         call check_partial_descriptor_roundtrip(descriptor(ipart), error)
         if (allocated(error)) return
         owner_matrix(:, ipart) = owner_count
         if (ipart > 1) then
            if (.not.cp2k_exchange_partial_descriptors_match( &
               & descriptor(1), descriptor(ipart))) then
               call test_failed(error, &
                  & "Multibatch forward descriptors differ")
               return
            end if
         end if
         call cp2k_exchange_partial_forward_pack(forward_buffer(ipart), &
            & forward_raw, error)
         if (allocated(error)) return
         if (ipart == 1) then
            allocate(forward_sum(size(forward_raw)), &
               & source=(0.0_wp, 0.0_wp))
         end if
         forward_sum = forward_sum+forward_raw
      end do
      owner_sum = sum(owner_matrix, dim=2)
      ncomplex = 4*size(density, 1)*size(density, 2)*size(density, 3)
      if (cp2k_exchange_partial_forward_size(forward_buffer(1)) /= &
         & int(ncomplex, int64)) then
         call test_failed(error, "Incorrect B=1 forward payload size")
         return
      end if
      call cp2k_exchange_partial_forward_unpack(forward_sum_buffer, &
         & descriptor(1), forward_sum, error)
      if (allocated(error)) return
      call cp2k_exchange_partial_export_forward(contributor_stream, &
         & contributor_descriptor, owner_count, contributor_forward_buffer, &
         & error)
      if (allocated(error)) return
      if (any(owner_count /= 0) .or. &
         & .not.cp2k_exchange_partial_descriptors_match(descriptor(1), &
         & contributor_descriptor)) then
         call test_failed(error, &
            & "Empty forward contributor descriptor/ownership differs")
         return
      end if
      if (ibatch == 1) stale_descriptor = descriptor(1)
      if (ibatch == 2) then
         call cp2k_exchange_partial_import_forward(stream(1), &
            & stale_descriptor, owner_sum, forward_sum_buffer, calc, mol, &
            & cache, cache_plan, error)
         if (.not.allocated(error)) then
            call test_failed(error, "A stale forward generation was accepted")
            return
         end if
         deallocate(error)
      end if
      do ipart = 1, npart
         call cp2k_exchange_partial_import_forward(stream(ipart), &
            & descriptor(1), owner_sum, forward_sum_buffer, calc, mol, cache, cache_plan, &
            & error)
         if (allocated(error)) return
         if (ibatch == 1) then
            if (ipart == 1) then
               call cp2k_exchange_partial_set_pull_ownership(stream(ipart), &
                  & pull_owner, ipart, npart, 0, error)
               if (.not.allocated(error)) then
                  call test_failed(error, &
                     & "Invalid distributed result owner was accepted")
                  return
               end if
               deallocate(error)
            end if
            call cp2k_exchange_partial_set_pull_ownership(stream(ipart), &
               & pull_owner, ipart, npart, 1, error)
            if (allocated(error)) return
            if (ipart == 1) then
               call cp2k_exchange_partial_set_pull_ownership(stream(ipart), &
                  & pull_owner, ipart, npart, 1, error)
               if (.not.allocated(error)) then
                  call test_failed(error, &
                     & "Duplicate distributed pull ownership was accepted")
                  return
               end if
               deallocate(error)
            end if
         end if
      end do
      do ipart = 1, npart
         call cp2k_exchange_stream_batch_apply(stream(ipart), calc, mol, &
            & cache, cache_plan, error)
         if (allocated(error)) return
         if (ibatch == 1 .and. ipart == 1) then
            call cp2k_exchange_stream_batch_advance(stream(ipart), &
               & contributor_done, error)
            if (.not.allocated(error)) then
               call test_failed(error, &
                  & "Distributed forward advance accepted a missing owned pull")
               return
            end if
            deallocate(error)
         end if
         if (ibatch == 1 .and. ipart == npart) then
            call cp2k_exchange_stream_batch_pull(stream(ipart), 1, &
               & overlap(:, :, 1), fock_partial(:, :, :, 1), error)
            if (.not.allocated(error)) then
               call test_failed(error, "An unowned forward pull was accepted")
               return
            end if
            deallocate(error)
         end if
         do ik = 1, nk
            if (pull_owner(ik) /= ipart) cycle
            call cp2k_exchange_stream_batch_pull(stream(ipart), ik, &
               & overlap(:, :, ik), fock_partial(:, :, :, ik), error)
            if (allocated(error)) return
         end do
         call cp2k_exchange_stream_batch_advance(stream(ipart), &
            & contributor_done, error)
         if (allocated(error)) return
         if (ipart == 1) then
            done = contributor_done
         else if (contributor_done .neqv. done) then
            call test_failed(error, &
               & "Distributed forward importer generations diverged")
            return
         end if
      end do
      call cp2k_exchange_partial_contributor_advance(contributor_stream, &
         & contributor_done, error)
      if (allocated(error)) return
      if (contributor_done .neqv. done) then
         call test_failed(error, &
            & "Forward contributor and importer generations diverged")
         return
      end if
      if (done .neqv. (ibatch == nk)) then
         call test_failed(error, "Unexpected multibatch forward completion")
         return
      end if
   end do
   call cp2k_exchange_stream_batch_result(stream(2), calc, mol, cache, cache_plan, &
      & vsh_partial, energy_partial, error)
   if (.not.allocated(error)) then
      call test_failed(error, &
         & "A non-owner retrieved the distributed forward scalar result")
      return
   end if
   deallocate(error)
   call cp2k_exchange_stream_batch_result(stream(1), calc, mol, cache, cache_plan, &
      & vsh_partial, energy_partial, error)
   if (allocated(error)) return
   maxerr = max(maxval(abs(fock_partial-fock_ref)), &
      & maxval(abs(vsh_partial-vsh_ref)), abs(energy_partial-energy_ref))
   do ipart = 1, npart
      call cp2k_exchange_stream_end(stream(ipart), error)
      if (allocated(error)) return
   end do
   call cp2k_exchange_partial_discard(contributor_stream, error)
   if (allocated(error)) return

   deallocate(owner_matrix)
   allocate(owner_matrix(nk, npart))
   allocate(overlap_adjoint(size(overlap_adjoint_ref, 1), &
      & size(overlap_adjoint_ref, 2), nk), source=(0.0_wp, 0.0_wp))
   allocate(gradient_partial(3, mol%nat), sigma_partial(3, 3))
   call next_partial_test_transaction(transaction_id)
   do ipart = 1, npart
      call cp2k_exchange_partial_begin(stream(ipart), calc, mol, cache, cache_plan, wfn, &
         & nmesh, kfrac, weights, 1, transaction_id, error, &
         & reverse_only=.true.)
      if (allocated(error)) return
   end do
   call cp2k_exchange_partial_begin(contributor_stream, calc, mol, cache, cache_plan, wfn, &
      & nmesh, kfrac, weights, 1, transaction_id, error, reverse_only=.true.)
   if (allocated(error)) return
   do ibatch = 1, nk
      owner_matrix = 0
      do ik = nk, 1, -1
         ipart = 1+modulo(ik-1, npart)
         call cp2k_exchange_partial_reverse_push(stream(ipart), ik, &
            & density(:, :, :, ik), overlap(:, :, ik), error)
         if (allocated(error)) return
      end do
      if (allocated(reverse_sum)) deallocate(reverse_sum)
      do ipart = 1, npart
         call cp2k_exchange_partial_export_reverse(stream(ipart), &
            & descriptor(ipart), owner_count, reverse_buffer(ipart), error)
         if (allocated(error)) return
         call check_partial_descriptor_roundtrip(descriptor(ipart), error)
         if (allocated(error)) return
         owner_matrix(:, ipart) = owner_count
         if (ipart > 1) then
            if (.not.cp2k_exchange_partial_descriptors_match( &
               & descriptor(1), descriptor(ipart))) then
               call test_failed(error, &
                  & "Multibatch reverse descriptors differ")
               return
            end if
         end if
         call cp2k_exchange_partial_reverse_pack(reverse_buffer(ipart), &
            & reverse_raw, error)
         if (allocated(error)) return
         if (ipart == 1) then
            allocate(reverse_sum(size(reverse_raw)), &
               & source=(0.0_wp, 0.0_wp))
         end if
         reverse_sum = reverse_sum+reverse_raw
      end do
      owner_sum = sum(owner_matrix, dim=2)
      ncomplex = 6*size(density, 1)*size(density, 2)*size(density, 3) &
         & +5*size(density, 1)*size(density, 3)
      if (cp2k_exchange_partial_reverse_size(reverse_buffer(1)) /= &
         & int(ncomplex, int64)) then
         call test_failed(error, "Incorrect B=1 reverse payload size")
         return
      end if
      call cp2k_exchange_partial_reverse_unpack(reverse_sum_buffer, &
         & descriptor(1), reverse_sum, error)
      if (allocated(error)) return
      call cp2k_exchange_partial_export_reverse(contributor_stream, &
         & contributor_descriptor, owner_count, contributor_reverse_buffer, &
         & error)
      if (allocated(error)) return
      if (any(owner_count /= 0) .or. &
         & .not.cp2k_exchange_partial_descriptors_match(descriptor(1), &
         & contributor_descriptor)) then
         call test_failed(error, &
            & "Empty reverse contributor descriptor/ownership differs")
         return
      end if
      if (ibatch == 1) stale_descriptor = descriptor(1)
      if (ibatch == 2) then
         call cp2k_exchange_partial_import_reverse(stream(1), &
            & stale_descriptor, owner_sum, reverse_sum_buffer, calc, mol, &
            & cache, cache_plan, error)
         if (.not.allocated(error)) then
            call test_failed(error, "A stale reverse generation was accepted")
            return
         end if
         deallocate(error)
      end if
      do ipart = 1, npart
         call cp2k_exchange_partial_import_reverse(stream(ipart), &
            & descriptor(1), owner_sum, reverse_sum_buffer, calc, mol, cache, cache_plan, &
            & error)
         if (allocated(error)) return
         if (ibatch == 1) then
            call cp2k_exchange_partial_set_pull_ownership(stream(ipart), &
               & pull_owner, ipart, npart, 1, error)
            if (allocated(error)) return
         end if
      end do
      do ipart = 1, npart
         call cp2k_exchange_stream_reverse_batch_apply(stream(ipart), calc, &
            & mol, cache, cache_plan, error)
         if (allocated(error)) return
         if (ibatch == 1 .and. ipart == 1) then
            call cp2k_exchange_stream_reverse_batch_advance(stream(ipart), &
               & contributor_done, error)
            if (.not.allocated(error)) then
               call test_failed(error, &
                  & "Distributed reverse advance accepted a missing owned pull")
               return
            end if
            deallocate(error)
         end if
         if (ibatch == 1 .and. ipart == npart) then
            call cp2k_exchange_stream_reverse_batch_pull(stream(ipart), 1, &
               & density(:, :, :, 1), overlap(:, :, 1), &
               & overlap_adjoint(:, :, 1), error)
            if (.not.allocated(error)) then
               call test_failed(error, "An unowned reverse pull was accepted")
               return
            end if
            deallocate(error)
         end if
         do ik = 1, nk
            if (pull_owner(ik) /= ipart) cycle
            call cp2k_exchange_stream_reverse_batch_pull(stream(ipart), ik, &
               & density(:, :, :, ik), overlap(:, :, ik), &
               & overlap_adjoint(:, :, ik), error)
            if (allocated(error)) return
         end do
         call cp2k_exchange_stream_reverse_batch_advance(stream(ipart), &
            & contributor_done, error)
         if (allocated(error)) return
         if (ipart == 1) then
            done = contributor_done
         else if (contributor_done .neqv. done) then
            call test_failed(error, &
               & "Distributed reverse importer generations diverged")
            return
         end if
      end do
      call cp2k_exchange_partial_contributor_advance(contributor_stream, &
         & contributor_done, error)
      if (allocated(error)) return
      if (contributor_done .neqv. done) then
         call test_failed(error, &
            & "Reverse contributor and importer generations diverged")
         return
      end if
      if (done .neqv. (ibatch == nk)) then
         call test_failed(error, "Unexpected multibatch reverse completion")
         return
      end if
   end do
   call cp2k_exchange_stream_reverse_result(stream(2), calc, mol, cache, cache_plan, &
      & gradient_partial, sigma_partial, error)
   if (.not.allocated(error)) then
      call test_failed(error, &
         & "A non-owner retrieved the distributed reverse direct result")
      return
   end if
   deallocate(error)
   call cp2k_exchange_stream_reverse_result(stream(1), calc, mol, cache, cache_plan, &
      & gradient_partial, sigma_partial, error)
   if (allocated(error)) return
   maxerr = max(maxerr, &
      & maxval(abs(overlap_adjoint-overlap_adjoint_ref)), &
      & maxval(abs(gradient_partial-gradient_ref)), &
      & maxval(abs(sigma_partial-sigma_ref)))
   do ipart = 1, npart
      call cp2k_exchange_stream_end(stream(ipart), error)
      if (allocated(error)) return
   end do
   call cp2k_exchange_partial_discard(contributor_stream, error)
   if (allocated(error)) return
end subroutine check_partial_exchange_multibatch

subroutine check_partial_exchange_distributed_images(calc, mol, cache, cache_plan, wfn, &
   & nmesh, kfrac, weights, density, overlap, error, maxerr)
   type(xtb_calculator), intent(inout) :: calc
   type(structure_type), intent(in) :: mol
   type(container_cache), intent(inout) :: cache
   type(cp2k_stream_plan_type), intent(inout) :: cache_plan
   type(wavefunction_type), intent(in) :: wfn
   integer, intent(in) :: nmesh(3)
   real(wp), intent(in) :: kfrac(:, :), weights(:)
   complex(wp), intent(in) :: density(:, :, :, :), overlap(:, :, :)
   type(error_type), allocatable, intent(out) :: error
   real(wp), intent(out) :: maxerr

   integer, parameter :: npart = 4, batch_size = 2
   type(cp2k_exchange_partial_descriptor_type) :: descriptor(npart)
   type(cp2k_exchange_partial_forward_type) :: forward_buffer(npart), &
      & forward_sum_buffer
   type(cp2k_exchange_partial_reverse_type) :: reverse_buffer(npart), &
      & reverse_sum_buffer
   type(cp2k_exchange_stream_type) :: range_stream, stream(npart)
   integer :: base_count, batch_index, extra_count, first_image, ik, &
      & image_count, ipart, nk
   integer(int64) :: transaction_id(2)
   integer, allocatable :: image_cover(:), owner_count(:), owner_matrix(:, :), &
      & owner_sum(:), pull_owner(:)
   logical :: done, done_part, first_batch
   real(wp) :: energy_err, energy_local, energy_partial, energy_ref, &
      & fock_err, gradient_err, overlap_adjoint_err, sigma_err, vsh_err
   real(wp), allocatable :: gradient_local(:, :), gradient_partial(:, :), &
      & gradient_ref(:, :), sigma_local(:, :), sigma_partial(:, :), &
      & sigma_ref(:, :), vsh_local(:), vsh_partial(:), vsh_ref(:), &
      & vsh_range(:)
   complex(wp), allocatable :: density_trial(:, :, :), &
      & fock_local(:, :, :), fock_range(:, :, :), &
      & fock_partial(:, :, :, :), &
      & fock_ref(:, :, :, :), forward_raw(:), forward_sum(:), &
      & overlap_adjoint(:, :, :), overlap_adjoint_local(:, :), &
      & overlap_adjoint_ref(:, :, :), overlap_trial(:, :), reverse_raw(:), &
      & reverse_sum(:)

   nk = size(weights)
   allocate(fock_ref(size(density, 1), size(density, 2), &
      & size(density, 3), nk), vsh_ref(size(wfn%qsh, 1)), &
      & overlap_adjoint_ref(size(overlap, 1), size(overlap, 2), nk), &
      & gradient_ref(3, mol%nat), sigma_ref(3, 3))
   call cp2k_exchange_kmesh(calc, mol, cache, cache_plan, wfn, nmesh, kfrac, weights, &
      & density, overlap, fock_ref, vsh_ref, energy_ref, error)
   if (allocated(error)) return
   call cp2k_exchange_kmesh_gradient(calc, mol, cache, cache_plan, nmesh, kfrac, weights, &
      & density, overlap, overlap_adjoint_ref, gradient_ref, sigma_ref, error)
   if (allocated(error)) return

   allocate(image_cover(nk), source=0)
   base_count = nk/npart
   extra_count = modulo(nk, npart)
   do ipart = 1, npart
      image_count = base_count+merge(1, 0, ipart <= extra_count)
      if (image_count > 0) then
         first_image = 1+(ipart-1)*base_count+min(ipart-1, extra_count)
         image_cover(first_image:first_image+image_count-1) = &
            & image_cover(first_image:first_image+image_count-1)+1
      end if
   end do
   if (any(image_cover /= 1)) then
      call test_failed(error, "Distributed image ranges do not exactly cover")
      return
   end if

   allocate(owner_matrix(nk, npart), owner_sum(nk), pull_owner(nk), source=0)
   allocate(fock_local(size(fock_ref, 1), size(fock_ref, 2), &
      & size(fock_ref, 3)), source=(0.0_wp, 0.0_wp))
   allocate(fock_partial(size(fock_ref, 1), size(fock_ref, 2), &
      & size(fock_ref, 3), nk), source=(0.0_wp, 0.0_wp))
   pull_owner = 1

   ! A distributed assignment must remain inside a pre-restricted base image
   ! window.  Exercise the public state machine and then finish the stream so
   ! the negative check cannot hide an invalid or leaking transaction.
   call next_partial_test_transaction(transaction_id)
   call cp2k_exchange_partial_begin(range_stream, calc, mol, cache, cache_plan, wfn, &
      & nmesh, kfrac, weights, nk, transaction_id, error, image_first=2, &
      & image_count=nk-1)
   if (allocated(error)) return
   do ik = 1, nk
      call cp2k_exchange_partial_push(range_stream, ik, density(:, :, :, ik), &
         & overlap(:, :, ik), error)
      if (allocated(error)) return
   end do
   call cp2k_exchange_partial_export_forward(range_stream, descriptor(1), &
      & owner_count, forward_buffer(1), error)
   if (allocated(error)) return
   call cp2k_exchange_partial_import_forward(range_stream, descriptor(1), &
      & owner_count, forward_buffer(1), calc, mol, cache, cache_plan, error)
   if (allocated(error)) return
   call cp2k_exchange_partial_set_image_range(range_stream, 1, 1, error)
   if (.not.allocated(error)) then
      call test_failed(error, &
         & "Distributed ownership escaped the base image window")
      return
   end if
   deallocate(error)
   call cp2k_exchange_partial_set_pull_ownership(range_stream, pull_owner, &
      & 1, 1, 1, error)
   if (allocated(error)) return
   call cp2k_exchange_partial_set_image_range(range_stream, 2, nk-1, error)
   if (.not.allocated(error)) then
      call test_failed(error, &
         & "Pull ownership was combined with distributed images")
      return
   end if
   deallocate(error)
   call cp2k_exchange_stream_batch_apply(range_stream, calc, mol, cache, cache_plan, error)
   if (allocated(error)) return
   allocate(fock_range(size(fock_ref, 1), size(fock_ref, 2), &
      & size(fock_ref, 3)), source=(0.0_wp, 0.0_wp))
   do ik = 1, nk
      fock_range = (0.0_wp, 0.0_wp)
      call cp2k_exchange_stream_batch_pull(range_stream, ik, &
         & overlap(:, :, ik), fock_range, error)
      if (allocated(error)) return
   end do
   call cp2k_exchange_stream_batch_advance(range_stream, done, error)
   if (allocated(error)) return
   if (.not.done) then
      call test_failed(error, "Restricted distributed image stream did not finish")
      return
   end if
   allocate(vsh_range(size(vsh_ref)))
   call cp2k_exchange_stream_batch_result(range_stream, calc, mol, cache, cache_plan, &
      & vsh_range, energy_local, error)
   if (allocated(error)) return
   call cp2k_exchange_stream_end(range_stream, error)
   if (allocated(error)) return
   deallocate(fock_range, vsh_range)

   call next_partial_test_transaction(transaction_id)
   do ipart = 1, npart
      call cp2k_exchange_partial_begin(stream(ipart), calc, mol, cache, cache_plan, wfn, &
         & nmesh, kfrac, weights, batch_size, transaction_id, error)
      if (allocated(error)) return
   end do

   done = .false.
   first_batch = .true.
   batch_index = 0
   do while (.not.done)
      batch_index = batch_index+1
      if (batch_index > nk) then
         call test_failed(error, "Distributed forward image batches did not finish")
         return
      end if
      owner_matrix = 0
      do ik = 1, nk
         ipart = 1+modulo(ik-1, npart)
         call cp2k_exchange_partial_push(stream(ipart), ik, &
            & density(:, :, :, ik), overlap(:, :, ik), error)
         if (allocated(error)) return
      end do
      if (allocated(forward_sum)) deallocate(forward_sum)
      do ipart = 1, npart
         call cp2k_exchange_partial_export_forward(stream(ipart), &
            & descriptor(ipart), owner_count, forward_buffer(ipart), error)
         if (allocated(error)) return
         owner_matrix(:, ipart) = owner_count
         if (ipart > 1) then
            if (.not.cp2k_exchange_partial_descriptors_match( &
               & descriptor(1), descriptor(ipart))) then
               call test_failed(error, &
                  & "Distributed forward image descriptors differ")
               return
            end if
         end if
         call cp2k_exchange_partial_forward_pack(forward_buffer(ipart), &
            & forward_raw, error)
         if (allocated(error)) return
         if (ipart == 1) then
            allocate(forward_sum(size(forward_raw)), &
               & source=(0.0_wp, 0.0_wp))
         end if
         forward_sum = forward_sum+forward_raw
      end do
      owner_sum = sum(owner_matrix, dim=2)
      call cp2k_exchange_partial_forward_unpack(forward_sum_buffer, &
         & descriptor(1), forward_sum, error)
      if (allocated(error)) return
      do ipart = 1, npart
         if (first_batch .and. ipart == 1) then
            call cp2k_exchange_partial_set_image_range(stream(ipart), &
               & 1, 1, error)
            if (.not.allocated(error)) then
               call test_failed(error, &
                  & "Distributed image range was accepted before import")
               return
            end if
            deallocate(error)
         end if
         call cp2k_exchange_partial_import_forward(stream(ipart), &
            & descriptor(1), owner_sum, forward_sum_buffer, calc, mol, cache, cache_plan, &
            & error)
         if (allocated(error)) return
         if (first_batch) then
            image_count = base_count+merge(1, 0, ipart <= extra_count)
            if (image_count > 0) then
               first_image = 1+(ipart-1)*base_count &
                  & +min(ipart-1, extra_count)
            else
               first_image = 0
            end if
            if (ipart == 1) then
               call cp2k_exchange_partial_set_image_range(stream(ipart), &
                  & 1, -1, error)
               if (.not.allocated(error)) then
                  call test_failed(error, &
                     & "Negative distributed image count was accepted")
                  return
               end if
               deallocate(error)
               call cp2k_exchange_partial_set_image_range(stream(ipart), &
                  & 0, 1, error)
               if (.not.allocated(error)) then
                  call test_failed(error, &
                     & "Positive distributed image range started at zero")
                  return
               end if
               deallocate(error)
               call cp2k_exchange_partial_set_image_range(stream(ipart), &
                  & nk, 2, error)
               if (.not.allocated(error)) then
                  call test_failed(error, &
                     & "Distributed image range crossed the mesh end")
                  return
               end if
               deallocate(error)
               call cp2k_exchange_partial_set_image_range(stream(ipart), &
                  & 1, 0, error)
               if (.not.allocated(error)) then
                  call test_failed(error, &
                     & "Invalid empty distributed image range was accepted")
                  return
               end if
               deallocate(error)
            end if
            call cp2k_exchange_partial_set_image_range(stream(ipart), &
               & first_image, image_count, error)
            if (allocated(error)) return
            if (ipart == 1) then
               call cp2k_exchange_partial_set_pull_ownership(stream(ipart), &
                  & pull_owner, ipart, npart, 1, error)
               if (.not.allocated(error)) then
                  call test_failed(error, &
                     & "Distributed images accepted disjoint k-point pulls")
                  return
               end if
               deallocate(error)
               call cp2k_exchange_partial_set_image_range(stream(ipart), &
                  & first_image, image_count, error)
               if (.not.allocated(error)) then
                  call test_failed(error, &
                     & "Duplicate distributed image range was accepted")
                  return
               end if
               deallocate(error)
            end if
         end if
      end do
      do ipart = 1, npart
         call cp2k_exchange_stream_batch_apply(stream(ipart), calc, mol, &
            & cache, cache_plan, error)
         if (allocated(error)) return
         if (first_batch .and. ipart == 1) then
            call cp2k_exchange_partial_set_image_range(stream(ipart), &
               & 1, 1, error)
            if (.not.allocated(error)) then
               call test_failed(error, &
                  & "Distributed image range followed kernel apply")
               return
            end if
            deallocate(error)
         end if
         do ik = 1, nk
            fock_local = (0.0_wp, 0.0_wp)
            call cp2k_exchange_stream_batch_pull(stream(ipart), ik, &
               & overlap(:, :, ik), fock_local, error)
            if (allocated(error)) return
            image_count = base_count+merge(1, 0, ipart <= extra_count)
            if (image_count == 0 .and. &
               & any(fock_local /= (0.0_wp, 0.0_wp))) then
               call test_failed(error, &
                  & "Empty image importer produced a forward Fock response")
               return
            end if
            fock_partial(:, :, :, ik) = fock_partial(:, :, :, ik)+fock_local
         end do
         call cp2k_exchange_stream_batch_advance(stream(ipart), &
            & done_part, error)
         if (allocated(error)) return
         if (first_batch .and. ipart == 1 .and. .not.done_part) then
            call cp2k_exchange_partial_set_image_range(stream(ipart), &
               & 1, 1, error)
            if (.not.allocated(error)) then
               call test_failed(error, &
                  & "Distributed image range followed an earlier batch")
               return
            end if
            deallocate(error)
         end if
         if (ipart == 1) then
            done = done_part
         else if (done_part .neqv. done) then
            call test_failed(error, &
               & "Distributed forward image generations diverged")
            return
         end if
      end do
      first_batch = .false.
   end do

   allocate(vsh_local(size(vsh_ref)), vsh_partial(size(vsh_ref)), &
      & source=0.0_wp)
   energy_partial = 0.0_wp
   do ipart = 1, npart
      call cp2k_exchange_stream_batch_result(stream(ipart), calc, mol, cache, cache_plan, &
         & vsh_local, energy_local, error)
      if (allocated(error)) return
      image_count = base_count+merge(1, 0, ipart <= extra_count)
      if (image_count == 0 .and. &
         & (any(vsh_local /= 0.0_wp) .or. energy_local /= 0.0_wp)) then
         call test_failed(error, &
            & "Empty image importer produced energy or shell potential")
         return
      end if
      vsh_partial = vsh_partial+vsh_local
      energy_partial = energy_partial+energy_local
      call cp2k_exchange_stream_reverse_begin(stream(ipart), mol, error)
      if (allocated(error)) return
   end do
   fock_err = maxval(abs(fock_partial-fock_ref))
   vsh_err = maxval(abs(vsh_partial-vsh_ref))
   energy_err = abs(energy_partial-energy_ref)
   maxerr = max(fock_err, vsh_err, energy_err)

   allocate(overlap_adjoint(size(overlap_adjoint_ref, 1), &
      & size(overlap_adjoint_ref, 2), nk), source=(0.0_wp, 0.0_wp))
   allocate(overlap_adjoint_local(size(overlap_adjoint_ref, 1), &
      & size(overlap_adjoint_ref, 2)), source=(0.0_wp, 0.0_wp))
   allocate(gradient_local(3, mol%nat), gradient_partial(3, mol%nat), &
      & sigma_local(3, 3), sigma_partial(3, 3), source=0.0_wp)
   done = .false.
   first_batch = .true.
   batch_index = 0
   do while (.not.done)
      batch_index = batch_index+1
      if (batch_index > nk) then
         call test_failed(error, "Distributed reverse image batches did not finish")
         return
      end if
      owner_matrix = 0
      do ik = 1, nk
         ipart = 1+modulo(ik-1, npart)
         if (first_batch .and. ik == 1) then
            allocate(density_trial, source=density(:, :, :, ik))
            density_trial(1, 1, 1) = density_trial(1, 1, 1) &
               & +cmplx(1.0e-7_wp, 0.0_wp, wp)
            call cp2k_exchange_partial_reverse_push(stream(ipart), ik, &
               & density_trial, overlap(:, :, ik), error)
            if (.not.allocated(error)) then
               call test_failed(error, &
                  & "Changed density crossed the forward/reverse boundary")
               return
            end if
            deallocate(error, density_trial)
            allocate(overlap_trial, source=overlap(:, :, ik))
            overlap_trial(1, 1) = overlap_trial(1, 1) &
               & +cmplx(1.0e-7_wp, 0.0_wp, wp)
            call cp2k_exchange_partial_reverse_push(stream(ipart), ik, &
               & density(:, :, :, ik), overlap_trial, error)
            if (.not.allocated(error)) then
               call test_failed(error, &
                  & "Changed overlap crossed the forward/reverse boundary")
               return
            end if
            deallocate(error, overlap_trial)
         end if
         call cp2k_exchange_partial_reverse_push(stream(ipart), ik, &
            & density(:, :, :, ik), overlap(:, :, ik), error)
         if (allocated(error)) return
      end do
      if (allocated(reverse_sum)) deallocate(reverse_sum)
      do ipart = 1, npart
         call cp2k_exchange_partial_export_reverse(stream(ipart), &
            & descriptor(ipart), owner_count, reverse_buffer(ipart), error)
         if (allocated(error)) return
         owner_matrix(:, ipart) = owner_count
         if (ipart > 1) then
            if (.not.cp2k_exchange_partial_descriptors_match( &
               & descriptor(1), descriptor(ipart))) then
               call test_failed(error, &
                  & "Distributed reverse image descriptors differ")
               return
            end if
         end if
         call cp2k_exchange_partial_reverse_pack(reverse_buffer(ipart), &
            & reverse_raw, error)
         if (allocated(error)) return
         if (ipart == 1) then
            allocate(reverse_sum(size(reverse_raw)), &
               & source=(0.0_wp, 0.0_wp))
         end if
         reverse_sum = reverse_sum+reverse_raw
      end do
      owner_sum = sum(owner_matrix, dim=2)
      call cp2k_exchange_partial_reverse_unpack(reverse_sum_buffer, &
         & descriptor(1), reverse_sum, error)
      if (allocated(error)) return
      do ipart = 1, npart
         call cp2k_exchange_partial_import_reverse(stream(ipart), &
            & descriptor(1), owner_sum, reverse_sum_buffer, calc, mol, cache, cache_plan, &
            & error)
         if (allocated(error)) return
         if (first_batch .and. ipart == 1) then
            call cp2k_exchange_partial_set_image_range(stream(ipart), &
               & 1, base_count+merge(1, 0, 1 <= extra_count), error)
            if (.not.allocated(error)) then
               call test_failed(error, &
                  & "Reverse phase changed a fixed distributed image range")
               return
            end if
            deallocate(error)
         end if
      end do
      do ipart = 1, npart
         call cp2k_exchange_stream_reverse_batch_apply(stream(ipart), calc, &
            & mol, cache, cache_plan, error)
         if (allocated(error)) return
         do ik = 1, nk
            overlap_adjoint_local = (0.0_wp, 0.0_wp)
            if (batch_index == 2 .and. ipart == 1 .and. ik == 1) then
               allocate(density_trial, source=density(:, :, :, ik))
               density_trial(1, 1, 1) = density_trial(1, 1, 1) &
                  & +cmplx(1.0e-7_wp, 0.0_wp, wp)
               call cp2k_exchange_stream_reverse_batch_pull(stream(ipart), &
                  & ik, density_trial, overlap(:, :, ik), &
                  & overlap_adjoint_local, error)
               if (.not.allocated(error)) then
                  call test_failed(error, &
                     & "Changed reverse density crossed an image batch")
                  return
               end if
               deallocate(error, density_trial)
               allocate(overlap_trial, source=overlap(:, :, ik))
               overlap_trial(1, 1) = overlap_trial(1, 1) &
                  & +cmplx(1.0e-7_wp, 0.0_wp, wp)
               call cp2k_exchange_stream_reverse_batch_pull(stream(ipart), &
                  & ik, density(:, :, :, ik), overlap_trial, &
                  & overlap_adjoint_local, error)
               if (.not.allocated(error)) then
                  call test_failed(error, &
                     & "Changed reverse overlap crossed an image batch")
                  return
               end if
               deallocate(error, overlap_trial)
            end if
            call cp2k_exchange_stream_reverse_batch_pull(stream(ipart), ik, &
               & density(:, :, :, ik), overlap(:, :, ik), &
               & overlap_adjoint_local, error)
            if (allocated(error)) return
            image_count = base_count+merge(1, 0, ipart <= extra_count)
            if (image_count == 0 .and. &
               & any(overlap_adjoint_local /= (0.0_wp, 0.0_wp))) then
               call test_failed(error, &
                  & "Empty image importer produced an overlap adjoint")
               return
            end if
            overlap_adjoint(:, :, ik) = overlap_adjoint(:, :, ik) &
               & +overlap_adjoint_local
         end do
         call cp2k_exchange_stream_reverse_batch_advance(stream(ipart), &
            & done_part, error)
         if (allocated(error)) return
         if (ipart == 1) then
            done = done_part
         else if (done_part .neqv. done) then
            call test_failed(error, &
               & "Distributed reverse image generations diverged")
            return
         end if
      end do
      first_batch = .false.
   end do

   do ipart = 1, npart
      call cp2k_exchange_stream_reverse_result(stream(ipart), calc, mol, &
         & cache, cache_plan, gradient_local, sigma_local, error)
      if (allocated(error)) return
      image_count = base_count+merge(1, 0, ipart <= extra_count)
      if (image_count == 0 .and. &
         & (any(gradient_local /= 0.0_wp) .or. &
         & any(sigma_local /= 0.0_wp))) then
         call test_failed(error, &
            & "Empty image importer produced force or stress response")
         return
      end if
      gradient_partial = gradient_partial+gradient_local
      sigma_partial = sigma_partial+sigma_local
      call cp2k_exchange_stream_end(stream(ipart), error)
      if (allocated(error)) return
   end do
   overlap_adjoint_err = maxval(abs(overlap_adjoint-overlap_adjoint_ref))
   gradient_err = maxval(abs(gradient_partial-gradient_ref))
   sigma_err = maxval(abs(sigma_partial-sigma_ref))
   maxerr = max(maxerr, overlap_adjoint_err, gradient_err, sigma_err)
   print *, "DISTRIBUTED_IMAGE_KERNEL_FORWARD_COMPONENTS", &
      & fock_err, vsh_err, energy_err
   print *, "DISTRIBUTED_IMAGE_KERNEL_REVERSE_COMPONENTS", &
      & overlap_adjoint_err, gradient_err, sigma_err

   ! Repeat the derivative contraction with fresh reverse-only streams.  This
   ! exercises the reverse setter branch independently of a preceding forward
   ! phase and permits a different host-side image partition when required.
   overlap_adjoint = (0.0_wp, 0.0_wp)
   gradient_partial = 0.0_wp
   sigma_partial = 0.0_wp
   call next_partial_test_transaction(transaction_id)
   do ipart = 1, npart
      call cp2k_exchange_partial_begin(stream(ipart), calc, mol, cache, cache_plan, wfn, &
         & nmesh, kfrac, weights, batch_size, transaction_id, error, &
         & reverse_only=.true.)
      if (allocated(error)) return
   end do
   done = .false.
   first_batch = .true.
   batch_index = 0
   do while (.not.done)
      batch_index = batch_index+1
      if (batch_index > nk) then
         call test_failed(error, &
            & "Reverse-only distributed image batches did not finish")
         return
      end if
      owner_matrix = 0
      do ik = 1, nk
         ipart = 1+modulo(ik-1, npart)
         call cp2k_exchange_partial_reverse_push(stream(ipart), ik, &
            & density(:, :, :, ik), overlap(:, :, ik), error)
         if (allocated(error)) return
      end do
      if (allocated(reverse_sum)) deallocate(reverse_sum)
      do ipart = 1, npart
         call cp2k_exchange_partial_export_reverse(stream(ipart), &
            & descriptor(ipart), owner_count, reverse_buffer(ipart), error)
         if (allocated(error)) return
         owner_matrix(:, ipart) = owner_count
         if (ipart > 1 .and. &
            & .not.cp2k_exchange_partial_descriptors_match( &
            & descriptor(1), descriptor(ipart))) then
            call test_failed(error, &
               & "Reverse-only distributed image descriptors differ")
            return
         end if
         call cp2k_exchange_partial_reverse_pack(reverse_buffer(ipart), &
            & reverse_raw, error)
         if (allocated(error)) return
         if (ipart == 1) then
            allocate(reverse_sum(size(reverse_raw)), &
               & source=(0.0_wp, 0.0_wp))
         end if
         reverse_sum = reverse_sum+reverse_raw
      end do
      owner_sum = sum(owner_matrix, dim=2)
      call cp2k_exchange_partial_reverse_unpack(reverse_sum_buffer, &
         & descriptor(1), reverse_sum, error)
      if (allocated(error)) return
      do ipart = 1, npart
         call cp2k_exchange_partial_import_reverse(stream(ipart), &
            & descriptor(1), owner_sum, reverse_sum_buffer, calc, mol, cache, cache_plan, &
            & error)
         if (allocated(error)) return
         if (first_batch) then
            image_count = base_count+merge(1, 0, ipart <= extra_count)
            if (image_count > 0) then
               first_image = 1+(ipart-1)*base_count &
                  & +min(ipart-1, extra_count)
            else
               first_image = 0
            end if
            call cp2k_exchange_partial_set_image_range(stream(ipart), &
               & first_image, image_count, error)
            if (allocated(error)) return
         end if
      end do
      do ipart = 1, npart
         call cp2k_exchange_stream_reverse_batch_apply(stream(ipart), calc, &
            & mol, cache, cache_plan, error)
         if (allocated(error)) return
         do ik = 1, nk
            overlap_adjoint_local = (0.0_wp, 0.0_wp)
            call cp2k_exchange_stream_reverse_batch_pull(stream(ipart), ik, &
               & density(:, :, :, ik), overlap(:, :, ik), &
               & overlap_adjoint_local, error)
            if (allocated(error)) return
            image_count = base_count+merge(1, 0, ipart <= extra_count)
            if (image_count == 0 .and. &
               & any(overlap_adjoint_local /= (0.0_wp, 0.0_wp))) then
               call test_failed(error, &
                  & "Empty reverse-only importer produced an adjoint")
               return
            end if
            overlap_adjoint(:, :, ik) = overlap_adjoint(:, :, ik) &
               & +overlap_adjoint_local
         end do
         call cp2k_exchange_stream_reverse_batch_advance(stream(ipart), &
            & done_part, error)
         if (allocated(error)) return
         if (ipart == 1) then
            done = done_part
         else if (done_part .neqv. done) then
            call test_failed(error, &
               & "Reverse-only distributed image generations diverged")
            return
         end if
      end do
      first_batch = .false.
   end do
   do ipart = 1, npart
      call cp2k_exchange_stream_reverse_result(stream(ipart), calc, mol, &
         & cache, cache_plan, gradient_local, sigma_local, error)
      if (allocated(error)) return
      image_count = base_count+merge(1, 0, ipart <= extra_count)
      if (image_count == 0 .and. &
         & (any(gradient_local /= 0.0_wp) .or. &
         & any(sigma_local /= 0.0_wp))) then
         call test_failed(error, &
            & "Empty reverse-only importer produced force or stress")
         return
      end if
      gradient_partial = gradient_partial+gradient_local
      sigma_partial = sigma_partial+sigma_local
      call cp2k_exchange_stream_end(stream(ipart), error)
      if (allocated(error)) return
   end do
   overlap_adjoint_err = maxval(abs(overlap_adjoint-overlap_adjoint_ref))
   gradient_err = maxval(abs(gradient_partial-gradient_ref))
   sigma_err = maxval(abs(sigma_partial-sigma_ref))
   maxerr = max(maxerr, overlap_adjoint_err, gradient_err, sigma_err)
   print *, "DISTRIBUTED_IMAGE_KERNEL_REVERSE_ONLY_COMPONENTS", &
      & overlap_adjoint_err, gradient_err, sigma_err
end subroutine check_partial_exchange_distributed_images

subroutine check_partial_exchange_range_local(calc, mol, cache, cache_plan, wfn, nmesh, &
   & kfrac, weights, density, overlap, error, maxerr, transform_backend)
   type(xtb_calculator), intent(inout) :: calc
   type(structure_type), intent(in) :: mol
   type(container_cache), intent(inout) :: cache
   type(cp2k_stream_plan_type), intent(inout) :: cache_plan
   type(wavefunction_type), intent(in) :: wfn
   integer, intent(in) :: nmesh(3)
   real(wp), intent(in) :: kfrac(:, :), weights(:)
   complex(wp), intent(in) :: density(:, :, :, :), overlap(:, :, :)
   type(error_type), allocatable, intent(out) :: error
   real(wp), intent(out) :: maxerr
   integer, intent(in), optional :: transform_backend

   integer, parameter :: npart = 4
   type(cp2k_exchange_partial_descriptor_type) :: descriptor(npart)
   type(cp2k_exchange_partial_forward_type) :: forward_buffer(npart), &
      & forward_full_buffer, forward_range_buffer(npart)
   type(cp2k_exchange_partial_reverse_type) :: reverse_buffer(npart), &
      & reverse_full_buffer, reverse_range_buffer(npart)
   type(cp2k_exchange_stream_type) :: contributor(npart), importer(npart), &
      & wrong_stream
   integer :: base_count, batch_count, batch_first, batch_index, batch_size, &
      & extra_count, first_image, ik, image_count, ipart, jpart, last_image, &
      & local_count, nk, range_first, range_last, selected_transform_backend, &
      & transform_pair_batch_size
   integer(int64) :: compact_peak, forward_range_scalars, &
      & forward_replicated_scalars, replicated_peak, reverse_range_scalars, &
      & reverse_replicated_scalars, transaction_id(2), wrong_transaction(2)
   integer, allocatable :: owner_bad(:), owner_count(:), owner_matrix(:, :), &
      & owner_sum(:), range_count(:), range_start(:)
   logical :: contributor_done, done, done_part, first_batch, transform_done, &
      & transform_done_part
   real(wp) :: energy_local, energy_partial, energy_ref
   real(wp), allocatable :: gradient_local(:, :), gradient_partial(:, :), &
      & gradient_ref(:, :), sigma_local(:, :), sigma_partial(:, :), &
      & sigma_ref(:, :), vsh_local(:), vsh_partial(:), vsh_ref(:)
   complex(wp), allocatable :: fock_local(:, :, :), &
      & fock_partial(:, :, :, :), fock_ref(:, :, :, :), &
      & forward_full_raw(:), forward_raw(:), forward_range_raw(:), &
      & overlap_adjoint(:, :, :), overlap_adjoint_local(:, :), &
      & overlap_adjoint_ref(:, :, :), reverse_full_raw(:), reverse_raw(:), &
      & reverse_range_raw(:)

   nk = size(weights)
   batch_size = min(nk, max(2, (nk+1)/2))
   selected_transform_backend = cp2k_exchange_transform_dense
   if (present(transform_backend)) &
      & selected_transform_backend = transform_backend
   transform_pair_batch_size = min(batch_size, max(1, size(density, 1)**2-1))
   allocate(fock_ref(size(density, 1), size(density, 2), &
      & size(density, 3), nk), vsh_ref(size(wfn%qsh, 1)), &
      & overlap_adjoint_ref(size(overlap, 1), size(overlap, 2), nk), &
      & gradient_ref(3, mol%nat), sigma_ref(3, 3))
   call cp2k_exchange_kmesh(calc, mol, cache, cache_plan, wfn, nmesh, kfrac, weights, &
      & density, overlap, fock_ref, vsh_ref, energy_ref, error)
   if (allocated(error)) return
   call cp2k_exchange_kmesh_gradient(calc, mol, cache, cache_plan, nmesh, kfrac, &
      & weights, density, overlap, overlap_adjoint_ref, gradient_ref, &
      & sigma_ref, error)
   if (allocated(error)) return

   allocate(range_start(npart), range_count(npart))
   base_count = nk/npart
   extra_count = modulo(nk, npart)
   do ipart = 1, npart
      range_count(ipart) = base_count+merge(1, 0, ipart <= extra_count)
      if (range_count(ipart) > 0) then
         range_start(ipart) = 1+(ipart-1)*base_count &
            & +min(ipart-1, extra_count)
      else
         range_start(ipart) = 0
      end if
   end do
   if (sum(range_count) /= nk) then
      call test_failed(error, "Range-local image partition is incomplete")
      return
   end if

   allocate(owner_matrix(nk, npart), owner_sum(nk), owner_bad(nk))
   allocate(fock_local(size(fock_ref, 1), size(fock_ref, 2), &
      & size(fock_ref, 3)), source=(0.0_wp, 0.0_wp))
   allocate(fock_partial(size(fock_ref, 1), size(fock_ref, 2), &
      & size(fock_ref, 3), nk), source=(0.0_wp, 0.0_wp))
   allocate(vsh_local(size(vsh_ref)), vsh_partial(size(vsh_ref)), &
      & source=0.0_wp)

   call next_partial_test_transaction(transaction_id)
   do ipart = 1, npart
      call cp2k_exchange_partial_begin(contributor(ipart), calc, mol, cache, cache_plan, &
         & wfn, nmesh, kfrac, weights, batch_size, transaction_id, error, &
         & transform_backend=selected_transform_backend, &
         & transform_pair_batch_size=transform_pair_batch_size)
      if (allocated(error)) return
      call cp2k_exchange_partial_begin(importer(ipart), calc, mol, cache, cache_plan, &
         & wfn, nmesh, kfrac, weights, batch_size, transaction_id, error, &
         & local_image_first=range_start(ipart), &
         & local_image_count=range_count(ipart), &
         & transform_backend=selected_transform_backend, &
         & transform_pair_batch_size=transform_pair_batch_size)
      if (allocated(error)) return
   end do
   replicated_peak = cp2k_exchange_stream_peak_complex_elements( &
      & contributor(1))
   compact_peak = 0_int64
   do ipart = 1, npart
      compact_peak = compact_peak &
         & +cp2k_exchange_stream_peak_complex_elements(importer(ipart))
   end do
   if (compact_peak >= int(npart, int64)*replicated_peak) then
      call test_failed(error, &
         & "Range-local importers did not reduce persistent complex storage")
      return
   end if

   done = .false.
   first_batch = .true.
   batch_index = 0
   forward_range_scalars = 0_int64
   forward_replicated_scalars = 0_int64
   do while (.not.done)
      batch_index = batch_index+1
      batch_first = 1+(batch_index-1)*batch_size
      batch_count = min(batch_size, nk-batch_first+1)
      owner_matrix = 0
      transform_done = .false.
      do while (.not.transform_done)
         do ik = 1, nk
            ipart = 1+modulo(ik-1, npart)
            call cp2k_exchange_partial_push(contributor(ipart), ik, &
               & density(:, :, :, ik), overlap(:, :, ik), error)
            if (allocated(error)) return
         end do
         if (selected_transform_backend == cp2k_exchange_transform_dense) then
            transform_done = .true.
         else
            do ipart = 1, npart
               call cp2k_exchange_stream_transform_advance( &
                  & contributor(ipart), transform_done_part, error)
               if (allocated(error)) return
               if (ipart == 1) then
                  transform_done = transform_done_part
               else if (transform_done_part .neqv. transform_done) then
                  call test_failed(error, &
                     & "Compact forward contributor transforms diverged")
                  return
               end if
            end do
         end if
      end do
      if (allocated(forward_full_raw)) deallocate(forward_full_raw)
      do ipart = 1, npart
         call cp2k_exchange_partial_export_forward(contributor(ipart), &
            & descriptor(ipart), owner_count, forward_buffer(ipart), error)
         if (allocated(error)) return
         owner_matrix(:, ipart) = owner_count
         if (ipart > 1 .and. &
            & .not.cp2k_exchange_partial_descriptors_match(descriptor(1), &
            & descriptor(ipart))) then
            call test_failed(error, &
               & "Range-local forward contributor descriptors differ")
            return
         end if
         call cp2k_exchange_partial_forward_pack(forward_buffer(ipart), &
            & forward_raw, error)
         if (allocated(error)) return
         if (ipart == 1) then
            allocate(forward_full_raw(size(forward_raw)), &
               & source=(0.0_wp, 0.0_wp))
         end if
         forward_full_raw = forward_full_raw+forward_raw
      end do
      owner_sum = sum(owner_matrix, dim=2)
      if (any(owner_sum /= 1)) then
         call test_failed(error, &
            & "Range-local forward contributors do not exactly cover k")
         return
      end if
      call cp2k_exchange_partial_forward_unpack(forward_full_buffer, &
         & descriptor(1), forward_full_raw, error)
      if (allocated(error)) return
      forward_replicated_scalars = forward_replicated_scalars &
         & +int(npart, int64)*size(forward_full_raw, kind=int64)

      do ipart = 1, npart
         image_count = range_count(ipart)
         first_image = range_start(ipart)
         if (image_count > 0) then
            last_image = first_image+image_count-1
            range_first = max(batch_first, first_image)
            range_last = min(batch_first+batch_count-1, last_image)
            local_count = max(0, range_last-range_first+1)
            if (local_count == 0) range_first = 0
         else
            range_first = 0
            local_count = 0
         end if
         if (allocated(forward_range_raw)) deallocate(forward_range_raw)
         do jpart = 1, npart
            call cp2k_exchange_partial_forward_pack_range( &
               & forward_buffer(jpart), descriptor(1), range_first, &
               & local_count, forward_raw, error)
            if (allocated(error)) return
            if (jpart == 1) then
               allocate(forward_range_raw(size(forward_raw)), &
                  & source=(0.0_wp, 0.0_wp))
            end if
            forward_range_raw = forward_range_raw+forward_raw
         end do
         forward_range_scalars = forward_range_scalars &
            & +size(forward_range_raw, kind=int64)
         call cp2k_exchange_partial_forward_unpack_range( &
            & forward_range_buffer(ipart), descriptor(1), range_first, &
            & local_count, forward_range_raw, error)
         if (allocated(error)) return
      end do

      if (first_batch) then
         call cp2k_exchange_partial_push(importer(1), 1, &
            & density(:, :, :, 1), overlap(:, :, 1), error)
         if (.not.allocated(error)) then
            call test_failed(error, &
               & "Range-local importer accepted a forward contribution")
            return
         end if
         deallocate(error)
         call cp2k_exchange_partial_import_forward(importer(1), &
            & descriptor(1), owner_sum, forward_full_buffer, calc, mol, cache, cache_plan, &
            & error)
         if (.not.allocated(error)) then
            call test_failed(error, &
               & "Range-local importer accepted a replicated payload")
            return
         end if
         deallocate(error)
         owner_bad = owner_sum
         owner_bad(1) = 0
         call cp2k_exchange_partial_import_forward_range(importer(1), &
            & descriptor(1), owner_bad, forward_range_buffer(1), calc, mol, &
            & cache, cache_plan, error)
         if (.not.allocated(error)) then
            call test_failed(error, &
               & "Range-local importer accepted incomplete k ownership")
            return
         end if
         deallocate(error)
         call cp2k_exchange_partial_forward_pack_range(forward_buffer(1), &
            & descriptor(1), 0, 1, forward_raw, error)
         if (.not.allocated(error)) then
            call test_failed(error, &
               & "Forward range pack accepted an invalid image interval")
            return
         end if
         deallocate(error)
         call next_partial_test_transaction(wrong_transaction)
         call cp2k_exchange_partial_begin(wrong_stream, calc, mol, cache, cache_plan, &
            & wfn, nmesh, kfrac, weights, batch_size, wrong_transaction, &
            & error, local_image_first=range_start(1), &
            & local_image_count=range_count(1), &
            & transform_backend=selected_transform_backend, &
            & transform_pair_batch_size=transform_pair_batch_size)
         if (allocated(error)) return
         call cp2k_exchange_partial_import_forward_range(wrong_stream, &
            & descriptor(1), owner_sum, forward_range_buffer(1), calc, mol, &
            & cache, cache_plan, error)
         if (.not.allocated(error)) then
            call test_failed(error, &
               & "Forward range import accepted a foreign transaction")
            return
         end if
         deallocate(error)
         call cp2k_exchange_stream_end(wrong_stream, error)
         if (.not.allocated(error)) then
            call test_failed(error, &
               & "Incomplete foreign range transaction closed successfully")
            return
         end if
         deallocate(error)
      end if

      do ipart = 1, npart
         call cp2k_exchange_partial_import_forward_range(importer(ipart), &
            & descriptor(1), owner_sum, forward_range_buffer(ipart), calc, &
            & mol, cache, cache_plan, error)
         if (allocated(error)) return
         call cp2k_exchange_stream_batch_apply(importer(ipart), calc, mol, &
            & cache, cache_plan, error)
         if (allocated(error)) return
         do ik = 1, nk
            fock_local = (0.0_wp, 0.0_wp)
            call cp2k_exchange_stream_batch_pull(importer(ipart), ik, &
               & overlap(:, :, ik), fock_local, error)
            if (allocated(error)) return
            if (range_count(ipart) == 0 .and. &
               & any(fock_local /= (0.0_wp, 0.0_wp))) then
               call test_failed(error, &
                  & "Empty range-local importer produced a Fock response")
               return
            end if
            fock_partial(:, :, :, ik) = fock_partial(:, :, :, ik)+fock_local
         end do
         call cp2k_exchange_stream_batch_advance(importer(ipart), done_part, &
            & error)
         if (allocated(error)) return
         if (ipart == 1) then
            done = done_part
         else if (done_part .neqv. done) then
            call test_failed(error, &
               & "Range-local forward importer generations diverged")
            return
         end if
      end do
      do ipart = 1, npart
         call cp2k_exchange_partial_contributor_advance(contributor(ipart), &
            & contributor_done, error)
         if (allocated(error)) return
         if (contributor_done .neqv. done) then
            call test_failed(error, &
               & "Range-local forward contributors diverged from importers")
            return
         end if
      end do
      first_batch = .false.
   end do
   if (batch_index < 2) then
      call test_failed(error, "Range-local forward gate did not span batches")
      return
   end if
   if (forward_range_scalars >= forward_replicated_scalars) then
      call test_failed(error, &
         & "Forward range payload did not reduce replicated communication")
      return
   end if

   energy_partial = 0.0_wp
   do ipart = 1, npart
      call cp2k_exchange_stream_batch_result(importer(ipart), calc, mol, &
         & cache, cache_plan, vsh_local, energy_local, error)
      if (allocated(error)) return
      if (range_count(ipart) == 0 .and. &
         & (any(vsh_local /= 0.0_wp) .or. energy_local /= 0.0_wp)) then
         call test_failed(error, &
            & "Empty range-local importer produced a scalar response")
         return
      end if
      vsh_partial = vsh_partial+vsh_local
      energy_partial = energy_partial+energy_local
      call cp2k_exchange_stream_end(importer(ipart), error)
      if (allocated(error)) return
      call cp2k_exchange_partial_discard(contributor(ipart), error)
      if (allocated(error)) return
   end do
   maxerr = max(maxval(abs(fock_partial-fock_ref)), &
      & maxval(abs(vsh_partial-vsh_ref)), abs(energy_partial-energy_ref))

   allocate(overlap_adjoint(size(overlap_adjoint_ref, 1), &
      & size(overlap_adjoint_ref, 2), nk), source=(0.0_wp, 0.0_wp))
   allocate(overlap_adjoint_local(size(overlap_adjoint_ref, 1), &
      & size(overlap_adjoint_ref, 2)), source=(0.0_wp, 0.0_wp))
   allocate(gradient_local(3, mol%nat), gradient_partial(3, mol%nat), &
      & sigma_local(3, 3), sigma_partial(3, 3), source=0.0_wp)

   call next_partial_test_transaction(transaction_id)
   do ipart = 1, npart
      call cp2k_exchange_partial_begin(contributor(ipart), calc, mol, cache, cache_plan, &
         & wfn, nmesh, kfrac, weights, batch_size, transaction_id, error, &
         & reverse_only=.true., &
         & transform_backend=selected_transform_backend, &
         & transform_pair_batch_size=transform_pair_batch_size)
      if (allocated(error)) return
      call cp2k_exchange_partial_begin(importer(ipart), calc, mol, cache, cache_plan, &
         & wfn, nmesh, kfrac, weights, batch_size, transaction_id, error, &
         & reverse_only=.true., local_image_first=range_start(ipart), &
         & local_image_count=range_count(ipart), &
         & transform_backend=selected_transform_backend, &
         & transform_pair_batch_size=transform_pair_batch_size)
      if (allocated(error)) return
   end do
   done = .false.
   first_batch = .true.
   batch_index = 0
   reverse_range_scalars = 0_int64
   reverse_replicated_scalars = 0_int64
   do while (.not.done)
      batch_index = batch_index+1
      do ipart = 1, npart
         call cp2k_exchange_stream_batch_bounds(importer(ipart), &
            & first_image, last_image, error)
         if (allocated(error)) return
         if (ipart == 1) then
            batch_first = first_image
            batch_count = last_image-first_image+1
         else if (first_image /= batch_first .or. &
            & last_image /= batch_first+batch_count-1) then
            call test_failed(error, &
               & "Range-local reverse importer batch bounds diverged")
            return
         end if
      end do
      owner_matrix = 0
      transform_done = .false.
      do while (.not.transform_done)
         do ik = 1, nk
            ipart = 1+modulo(ik-1, npart)
            call cp2k_exchange_partial_reverse_push(contributor(ipart), ik, &
               & density(:, :, :, ik), overlap(:, :, ik), error)
            if (allocated(error)) return
         end do
         if (selected_transform_backend == cp2k_exchange_transform_dense) then
            transform_done = .true.
         else
            do ipart = 1, npart
               call cp2k_exchange_stream_reverse_transform_advance( &
                  & contributor(ipart), transform_done_part, error)
               if (allocated(error)) return
               if (ipart == 1) then
                  transform_done = transform_done_part
               else if (transform_done_part .neqv. transform_done) then
                  call test_failed(error, &
                     & "Compact reverse contributor transforms diverged")
                  return
               end if
            end do
         end if
      end do
      if (allocated(reverse_full_raw)) deallocate(reverse_full_raw)
      do ipart = 1, npart
         call cp2k_exchange_partial_export_reverse(contributor(ipart), &
            & descriptor(ipart), owner_count, reverse_buffer(ipart), error)
         if (allocated(error)) return
         owner_matrix(:, ipart) = owner_count
         if (ipart > 1 .and. &
            & .not.cp2k_exchange_partial_descriptors_match(descriptor(1), &
            & descriptor(ipart))) then
            call test_failed(error, &
               & "Range-local reverse contributor descriptors differ")
            return
         end if
         call cp2k_exchange_partial_reverse_pack(reverse_buffer(ipart), &
            & reverse_raw, error)
         if (allocated(error)) return
         if (ipart == 1) then
            allocate(reverse_full_raw(size(reverse_raw)), &
               & source=(0.0_wp, 0.0_wp))
         end if
         reverse_full_raw = reverse_full_raw+reverse_raw
      end do
      owner_sum = sum(owner_matrix, dim=2)
      if (any(owner_sum /= 1)) then
         call test_failed(error, &
            & "Range-local reverse contributors do not exactly cover k")
         return
      end if
      call cp2k_exchange_partial_reverse_unpack(reverse_full_buffer, &
         & descriptor(1), reverse_full_raw, error)
      if (allocated(error)) return
      reverse_replicated_scalars = reverse_replicated_scalars &
         & +int(npart, int64)*size(reverse_full_raw, kind=int64)

      do ipart = 1, npart
         image_count = range_count(ipart)
         first_image = range_start(ipart)
         if (image_count > 0) then
            last_image = first_image+image_count-1
            range_first = max(batch_first, first_image)
            range_last = min(batch_first+batch_count-1, last_image)
            local_count = max(0, range_last-range_first+1)
            if (local_count == 0) range_first = 0
         else
            range_first = 0
            local_count = 0
         end if
         if (allocated(reverse_range_raw)) deallocate(reverse_range_raw)
         do jpart = 1, npart
            call cp2k_exchange_partial_reverse_pack_range( &
               & reverse_buffer(jpart), descriptor(1), range_first, &
               & local_count, reverse_raw, error)
            if (allocated(error)) return
            if (jpart == 1) then
               allocate(reverse_range_raw(size(reverse_raw)), &
                  & source=(0.0_wp, 0.0_wp))
            end if
            reverse_range_raw = reverse_range_raw+reverse_raw
         end do
         reverse_range_scalars = reverse_range_scalars &
            & +size(reverse_range_raw, kind=int64)
         call cp2k_exchange_partial_reverse_unpack_range( &
            & reverse_range_buffer(ipart), descriptor(1), range_first, &
            & local_count, reverse_range_raw, error)
         if (allocated(error)) return
      end do

      if (first_batch) then
         call cp2k_exchange_partial_import_reverse(importer(1), &
            & descriptor(1), owner_sum, reverse_full_buffer, calc, mol, cache, cache_plan, &
            & error)
         if (.not.allocated(error)) then
            call test_failed(error, &
               & "Range-local importer accepted a replicated reverse payload")
            return
         end if
         deallocate(error)
         call cp2k_exchange_partial_reverse_pack_range(reverse_buffer(1), &
            & descriptor(1), 0, 1, reverse_raw, error)
         if (.not.allocated(error)) then
            call test_failed(error, &
               & "Reverse range pack accepted an invalid image interval")
            return
         end if
         deallocate(error)
      end if

      do ipart = 1, npart
         call cp2k_exchange_partial_import_reverse_range(importer(ipart), &
            & descriptor(1), owner_sum, reverse_range_buffer(ipart), calc, &
            & mol, cache, cache_plan, error)
         if (allocated(error)) return
         call cp2k_exchange_stream_reverse_batch_apply(importer(ipart), calc, &
            & mol, cache, cache_plan, error)
         if (allocated(error)) return
         do ik = 1, nk
            overlap_adjoint_local = (0.0_wp, 0.0_wp)
            call cp2k_exchange_stream_reverse_batch_pull(importer(ipart), &
               & ik, density(:, :, :, ik), overlap(:, :, ik), &
               & overlap_adjoint_local, error)
            if (allocated(error)) return
            if (range_count(ipart) == 0 .and. &
               & any(overlap_adjoint_local /= (0.0_wp, 0.0_wp))) then
               call test_failed(error, &
                  & "Empty range-local importer produced an adjoint")
               return
            end if
            overlap_adjoint(:, :, ik) = overlap_adjoint(:, :, ik) &
               & +overlap_adjoint_local
         end do
         call cp2k_exchange_stream_reverse_batch_advance(importer(ipart), &
            & done_part, error)
         if (allocated(error)) return
         if (ipart == 1) then
            done = done_part
         else if (done_part .neqv. done) then
            call test_failed(error, &
               & "Range-local reverse importer generations diverged")
            return
         end if
      end do
      do ipart = 1, npart
         call cp2k_exchange_partial_contributor_advance(contributor(ipart), &
            & contributor_done, error)
         if (allocated(error)) return
         if (contributor_done .neqv. done) then
            call test_failed(error, &
               & "Range-local reverse contributors diverged from importers")
            return
         end if
      end do
      first_batch = .false.
   end do
   call cp2k_exchange_stream_batch_bounds(importer(1), first_image, &
      & last_image, error)
   if (allocated(error)) then
      if (trim(error%message) /= &
         & "g-xTB exchange image batches are complete") then
         deallocate(error)
         call test_failed(error, &
            & "Unexpected completed reverse batch-bounds diagnostic")
         return
      end if
      deallocate(error)
   else
      call test_failed(error, &
         & "Completed reverse image batches returned bounds")
      return
   end if
   if (batch_index < 2) then
      call test_failed(error, "Range-local reverse gate did not span batches")
      return
   end if
   if (reverse_range_scalars >= reverse_replicated_scalars) then
      call test_failed(error, &
         & "Reverse range payload did not reduce replicated communication")
      return
   end if

   do ipart = 1, npart
      call cp2k_exchange_stream_reverse_result(importer(ipart), calc, mol, &
         & cache, cache_plan, gradient_local, sigma_local, error)
      if (allocated(error)) return
      if (range_count(ipart) == 0 .and. &
         & (any(gradient_local /= 0.0_wp) .or. &
         & any(sigma_local /= 0.0_wp))) then
         call test_failed(error, &
            & "Empty range-local importer produced force or stress")
         return
      end if
      gradient_partial = gradient_partial+gradient_local
      sigma_partial = sigma_partial+sigma_local
      call cp2k_exchange_stream_end(importer(ipart), error)
      if (allocated(error)) return
      call cp2k_exchange_partial_discard(contributor(ipart), error)
      if (allocated(error)) return
   end do
   maxerr = max(maxerr, &
      & maxval(abs(overlap_adjoint-overlap_adjoint_ref)), &
      & maxval(abs(gradient_partial-gradient_ref)), &
      & maxval(abs(sigma_partial-sigma_ref)))
   print *, "RANGE_LOCAL_IMAGE_IMPORT_MAX_RESID", maxerr
   print *, "RANGE_LOCAL_IMAGE_IMPORT_PEAK_COMPLEX", replicated_peak, &
      & compact_peak
   print *, "RANGE_LOCAL_IMAGE_IMPORT_PAYLOAD_COMPLEX", &
      & forward_replicated_scalars, forward_range_scalars, &
      & reverse_replicated_scalars, reverse_range_scalars
end subroutine check_partial_exchange_range_local

subroutine check_exchange_transform_selectors(calc, mol, cache, cache_plan, wfn, nmesh, &
   & kfrac, weights, density, overlap, error, separable_max, fft_max, &
   & fft_forward, fft_reverse)
   type(xtb_calculator), intent(inout) :: calc
   type(structure_type), intent(in) :: mol
   type(container_cache), intent(inout) :: cache
   type(cp2k_stream_plan_type), intent(inout) :: cache_plan
   type(wavefunction_type), intent(in) :: wfn
   integer, intent(in) :: nmesh(3)
   real(wp), intent(in) :: kfrac(:, :), weights(:)
   complex(wp), intent(in) :: density(:, :, :, :), overlap(:, :, :)
   type(error_type), allocatable, intent(out) :: error
   real(wp), intent(out) :: separable_max, fft_max
   real(wp), intent(out) :: fft_forward(3), fft_reverse(3)

   real(wp) :: energy_ref, energy_test
   real(wp), allocatable :: gradient_ref(:, :), gradient_test(:, :), &
      & sigma_ref(:, :), sigma_test(:, :), vsh_ref(:), vsh_test(:)
   complex(wp), allocatable :: fock_ref(:, :, :, :), &
      & fock_test(:, :, :, :), overlap_grad_ref(:, :, :), &
      & overlap_grad_test(:, :, :)

   allocate(fock_ref(size(density, 1), size(density, 2), size(density, 3), &
      & size(density, 4)), fock_test(size(density, 1), size(density, 2), &
      & size(density, 3), size(density, 4)), vsh_ref(size(wfn%qsh, 1)), &
      & vsh_test(size(wfn%qsh, 1)), &
      & overlap_grad_ref(size(overlap, 1), size(overlap, 2), &
      & size(overlap, 3)), overlap_grad_test(size(overlap, 1), &
      & size(overlap, 2), size(overlap, 3)), &
      & gradient_ref(3, mol%nat), gradient_test(3, mol%nat), &
      & sigma_ref(3, 3), sigma_test(3, 3))

   call cp2k_exchange_kmesh(calc, mol, cache, cache_plan, wfn, nmesh, kfrac, weights, &
      & density, overlap, fock_ref, vsh_ref, energy_ref, error, &
      & transform_backend=cp2k_exchange_transform_dense)
   if (allocated(error)) return
   call cp2k_exchange_kmesh_gradient(calc, mol, cache, cache_plan, nmesh, kfrac, &
      & weights, density, overlap, overlap_grad_ref, gradient_ref, sigma_ref, &
      & error, transform_backend=cp2k_exchange_transform_dense)
   if (allocated(error)) return

   call cp2k_exchange_kmesh(calc, mol, cache, cache_plan, wfn, nmesh, kfrac, weights, &
      & density, overlap, fock_test, vsh_test, energy_test, error, &
      & transform_backend=cp2k_exchange_transform_separable_dft)
   if (allocated(error)) return
   separable_max = max(maxval(abs(fock_test-fock_ref)), &
      & maxval(abs(vsh_test-vsh_ref)), abs(energy_test-energy_ref))
   call cp2k_exchange_kmesh_gradient(calc, mol, cache, cache_plan, nmesh, kfrac, &
      & weights, density, overlap, overlap_grad_test, gradient_test, &
      & sigma_test, error, &
      & transform_backend=cp2k_exchange_transform_separable_dft)
   if (allocated(error)) return
   separable_max = max(separable_max, &
      & maxval(abs(overlap_grad_test-overlap_grad_ref)), &
      & maxval(abs(gradient_test-gradient_ref)), &
      & maxval(abs(sigma_test-sigma_ref)))

   call cp2k_exchange_kmesh(calc, mol, cache, cache_plan, wfn, nmesh, kfrac, weights, &
      & density, overlap, fock_test, vsh_test, energy_test, error, &
      & transform_backend=cp2k_exchange_transform_mixed_radix_fft)
   if (allocated(error)) return
   fft_forward = [maxval(abs(fock_test-fock_ref)), &
      & maxval(abs(vsh_test-vsh_ref)), abs(energy_test-energy_ref)]
   call cp2k_exchange_kmesh_gradient(calc, mol, cache, cache_plan, nmesh, kfrac, &
      & weights, density, overlap, overlap_grad_test, gradient_test, &
      & sigma_test, error, &
      & transform_backend=cp2k_exchange_transform_mixed_radix_fft)
   if (allocated(error)) return
   fft_reverse = [maxval(abs(overlap_grad_test-overlap_grad_ref)), &
      & maxval(abs(gradient_test-gradient_ref)), &
      & maxval(abs(sigma_test-sigma_ref))]
   fft_max = max(maxval(fft_forward), maxval(fft_reverse))
end subroutine check_exchange_transform_selectors

subroutine check_compact_stream_composition(calc, mol, cache, cache_plan, wfn, nmesh, &
   & kfrac, weights, density, overlap, fock_ref, vsh_ref, energy_ref, &
   & overlap_grad_ref, gradient_ref, sigma_ref, transform_backend, &
   & residual, peak_elements, error, adaptive_tiles, adaptive_images)
   type(xtb_calculator), intent(inout) :: calc
   type(structure_type), intent(in) :: mol
   type(container_cache), intent(inout) :: cache
   type(cp2k_stream_plan_type), intent(inout) :: cache_plan
   type(wavefunction_type), intent(in) :: wfn
   integer, intent(in) :: nmesh(3), transform_backend
   real(wp), intent(in) :: kfrac(:, :), weights(:), vsh_ref(:), energy_ref
   complex(wp), intent(in) :: density(:, :, :, :), overlap(:, :, :), &
      & fock_ref(:, :, :, :), overlap_grad_ref(:, :, :)
   real(wp), intent(in) :: gradient_ref(:, :), sigma_ref(:, :)
   real(wp), intent(out) :: residual
   integer(int64), intent(out) :: peak_elements
   type(error_type), allocatable, intent(out) :: error

   logical, intent(in), optional :: adaptive_tiles, adaptive_images

   type(cp2k_exchange_stream_type) :: stream
   integer :: ik, itile, nao, nk, npair, pair_batch, tile_count, image_batch
   integer(int64) :: reverse_peak
   logical :: batch_done, tile_done
   real(wp) :: energy
   real(wp), allocatable :: gradient(:, :), sigma(:, :), vsh(:)
   complex(wp), allocatable :: fock(:, :, :, :), overlap_grad(:, :, :)

   nao = size(density, 1)
   nk = size(density, 4)
   npair = nao*nao
   image_batch = 1
   if (present(adaptive_images)) then
      if (adaptive_images) image_batch = cp2k_exchange_stream_image_batch_size(nao, nk, wfn%nspin)
   end if
   pair_batch = max(1, min(2, npair-1))
   if (present(adaptive_tiles)) then
      if (adaptive_tiles) pair_batch = cp2k_exchange_stream_pair_tile_size(nao, nk, wfn%nspin, 3)
   end if
   tile_count = (npair+pair_batch-1)/pair_batch
   allocate(fock(size(density, 1), size(density, 2), size(density, 3), &
      & size(density, 4)), source=(0.0_wp, 0.0_wp))
   allocate(vsh(size(vsh_ref)), source=0.0_wp)

   call cp2k_exchange_stream_begin(stream, calc, mol, cache, cache_plan, wfn, nmesh, &
      & kfrac, weights, error, mode=cp2k_exchange_stream_batched, &
      & batch_size=image_batch, transform_backend=transform_backend, &
      & transform_pair_batch_size=pair_batch)
   if (allocated(error)) return
   if (npair > 1 .and. cp2k_exchange_stream_has_full_mesh_storage(stream)) then
      call test_failed(error, &
         & "Compact forward stream retained a complete k/AO-pair mesh")
      return
   end if
   peak_elements = cp2k_exchange_stream_peak_complex_elements(stream)
   batch_done = .false.
   do while (.not.batch_done)
      do itile = 1, tile_count
         do ik = 1, nk
            call cp2k_exchange_stream_push(stream, ik, &
               & density(:, :, :, ik), overlap(:, :, ik), error)
            if (allocated(error)) return
         end do
         call cp2k_exchange_stream_transform_advance(stream, tile_done, error)
         if (allocated(error)) return
         if (tile_done .neqv. (itile == tile_count)) then
            call test_failed(error, &
               & "Compact forward stream reported an inconsistent tile count")
            return
         end if
      end do
      call cp2k_exchange_stream_batch_apply(stream, calc, mol, cache, cache_plan, error)
      if (allocated(error)) return
      do ik = 1, nk
         call cp2k_exchange_stream_batch_pull(stream, ik, overlap(:, :, ik), &
            & fock(:, :, :, ik), error)
         if (allocated(error)) return
      end do
      call cp2k_exchange_stream_batch_advance(stream, batch_done, error)
      if (allocated(error)) return
   end do
   call cp2k_exchange_stream_batch_result(stream, calc, mol, cache, cache_plan, vsh, &
      & energy, error)
   if (allocated(error)) return
   residual = max(maxval(abs(fock-fock_ref)), maxval(abs(vsh-vsh_ref)), &
      & abs(energy-energy_ref))
   call cp2k_exchange_stream_end(stream, error)
   if (allocated(error)) return

   allocate(overlap_grad(size(overlap, 1), size(overlap, 2), nk), &
      & source=(0.0_wp, 0.0_wp))
   allocate(gradient(3, mol%nat), sigma(3, 3), source=0.0_wp)
   if (present(adaptive_tiles)) then
      if (adaptive_tiles) pair_batch = cp2k_exchange_stream_pair_tile_size(nao, nk, wfn%nspin, 6)
   end if
   tile_count = (npair+pair_batch-1)/pair_batch
   call cp2k_exchange_stream_begin(stream, calc, mol, cache, cache_plan, wfn, nmesh, &
      & kfrac, weights, error, mode=cp2k_exchange_stream_batched, &
      & batch_size=image_batch, reverse_only=.true., &
      & transform_backend=transform_backend, &
      & transform_pair_batch_size=pair_batch)
   if (allocated(error)) return
   if (npair > 1 .and. cp2k_exchange_stream_has_full_mesh_storage(stream)) then
      call test_failed(error, &
         & "Compact reverse stream retained a complete k/AO-pair mesh")
      return
   end if
   reverse_peak = cp2k_exchange_stream_peak_complex_elements(stream)
   peak_elements = max(peak_elements, reverse_peak)
   batch_done = .false.
   do while (.not.batch_done)
      do itile = 1, tile_count
         do ik = 1, nk
            call cp2k_exchange_stream_reverse_push(stream, ik, &
               & density(:, :, :, ik), overlap(:, :, ik), error)
            if (allocated(error)) return
         end do
         call cp2k_exchange_stream_reverse_transform_advance(stream, &
            & tile_done, error)
         if (allocated(error)) return
         if (tile_done .neqv. (itile == tile_count)) then
            call test_failed(error, &
               & "Compact reverse stream reported an inconsistent tile count")
            return
         end if
      end do
      call cp2k_exchange_stream_reverse_batch_apply(stream, calc, mol, cache, cache_plan, &
         & error)
      if (allocated(error)) return
      do ik = 1, nk
         call cp2k_exchange_stream_reverse_batch_pull(stream, ik, &
            & density(:, :, :, ik), overlap(:, :, ik), &
            & overlap_grad(:, :, ik), error)
         if (allocated(error)) return
      end do
      call cp2k_exchange_stream_reverse_batch_advance(stream, batch_done, &
         & error)
      if (allocated(error)) return
   end do
   call cp2k_exchange_stream_reverse_result(stream, calc, mol, cache, cache_plan, &
      & gradient, sigma, error)
   if (allocated(error)) return
   residual = max(residual, maxval(abs(overlap_grad-overlap_grad_ref)), &
      & maxval(abs(gradient-gradient_ref)), maxval(abs(sigma-sigma_ref)))
   call cp2k_exchange_stream_end(stream, error)
end subroutine check_compact_stream_composition

subroutine require_bvk_model_mismatch(exchange, cache, field, error)
   type(exchange_fock), intent(in) :: exchange
   type(cp2k_stream_plan_type), intent(in) :: cache
   character(len=*), intent(in) :: field
   type(error_type), allocatable, intent(out) :: error

   if (cache%model_matches(exchange)) then
      call test_failed(error, "BvK model signature ignored "//trim(field))
   end if
end subroutine require_bvk_model_mismatch

subroutine check_bvk_model_signature(exchange, cache, error)
   type(exchange_fock), intent(inout) :: exchange
   type(cp2k_stream_plan_type), intent(in) :: cache
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: integer_delta = 1
   real(wp), parameter :: real_delta = 1.0e-6_wp

   if (.not.cache%model_matches(exchange)) then
      call test_failed(error, "Fresh BvK model signature does not match")
      return
   end if

   exchange%nao = exchange%nao+integer_delta
   call require_bvk_model_mismatch(exchange, cache, "nao", error)
   exchange%nao = exchange%nao-integer_delta
   if (allocated(error)) return
   exchange%nsh = exchange%nsh+integer_delta
   call require_bvk_model_mismatch(exchange, cache, "nsh", error)
   exchange%nsh = exchange%nsh-integer_delta
   if (allocated(error)) return
   exchange%maxsh = exchange%maxsh+integer_delta
   call require_bvk_model_mismatch(exchange, cache, "maxsh", error)
   exchange%maxsh = exchange%maxsh-integer_delta
   if (allocated(error)) return

   exchange%nsh_id(1) = exchange%nsh_id(1)+integer_delta
   call require_bvk_model_mismatch(exchange, cache, "nsh_id", error)
   exchange%nsh_id(1) = exchange%nsh_id(1)-integer_delta
   if (allocated(error)) return
   exchange%nao_sh(1) = exchange%nao_sh(1)+integer_delta
   call require_bvk_model_mismatch(exchange, cache, "nao_sh", error)
   exchange%nao_sh(1) = exchange%nao_sh(1)-integer_delta
   if (allocated(error)) return
   exchange%ish_at(1) = exchange%ish_at(1)+integer_delta
   call require_bvk_model_mismatch(exchange, cache, "ish_at", error)
   exchange%ish_at(1) = exchange%ish_at(1)-integer_delta
   if (allocated(error)) return
   exchange%iao_sh(1) = exchange%iao_sh(1)+integer_delta
   call require_bvk_model_mismatch(exchange, cache, "iao_sh", error)
   exchange%iao_sh(1) = exchange%iao_sh(1)-integer_delta
   if (allocated(error)) return

   exchange%frscale = exchange%frscale+real_delta
   call require_bvk_model_mismatch(exchange, cache, "frscale", error)
   exchange%frscale = cache%bvk_model_frscale
   if (allocated(error)) return
   exchange%onsite_frscale = exchange%onsite_frscale+real_delta
   call require_bvk_model_mismatch(exchange, cache, "onsite_frscale", error)
   exchange%onsite_frscale = cache%bvk_model_onsite_frscale
   if (allocated(error)) return
   exchange%omega = exchange%omega+real_delta
   call require_bvk_model_mismatch(exchange, cache, "omega", error)
   exchange%omega = cache%bvk_model_omega
   if (allocated(error)) return
   exchange%lrscale = exchange%lrscale+real_delta
   call require_bvk_model_mismatch(exchange, cache, "lrscale", error)
   exchange%lrscale = cache%bvk_model_lrscale
   if (allocated(error)) return
   exchange%ondiag_scale = exchange%ondiag_scale+real_delta
   call require_bvk_model_mismatch(exchange, cache, "ondiag_scale", error)
   exchange%ondiag_scale = cache%bvk_model_ondiag_scale
   if (allocated(error)) return
   exchange%hubbard_exp = exchange%hubbard_exp+real_delta
   call require_bvk_model_mismatch(exchange, cache, "hubbard_exp", error)
   exchange%hubbard_exp = cache%bvk_model_hubbard_exp
   if (allocated(error)) return
   exchange%hubbard_exp_r0 = exchange%hubbard_exp_r0+real_delta
   call require_bvk_model_mismatch(exchange, cache, "hubbard_exp_r0", error)
   exchange%hubbard_exp_r0 = cache%bvk_model_hubbard_exp_r0
   if (allocated(error)) return
   exchange%gexp = exchange%gexp+real_delta
   call require_bvk_model_mismatch(exchange, cache, "gexp", error)
   exchange%gexp = cache%bvk_model_gexp
   if (allocated(error)) return
   exchange%corr_exp = exchange%corr_exp+real_delta
   call require_bvk_model_mismatch(exchange, cache, "corr_exp", error)
   exchange%corr_exp = cache%bvk_model_corr_exp
   if (allocated(error)) return

   exchange%hubbard(1, 1, 1, 1) = exchange%hubbard(1, 1, 1, 1)+real_delta
   call require_bvk_model_mismatch(exchange, cache, "hubbard", error)
   exchange%hubbard(1, 1, 1, 1) = cache%bvk_model_hubbard(1, 1, 1, 1)
   if (allocated(error)) return
   exchange%onecxints(1, 1, 1) = exchange%onecxints(1, 1, 1)+real_delta
   call require_bvk_model_mismatch(exchange, cache, "onecxints", error)
   exchange%onecxints(1, 1, 1) = cache%bvk_model_onecxints(1, 1, 1)
   if (allocated(error)) return
   exchange%offdiag_scale(1, 1, 1, 1) = &
      & exchange%offdiag_scale(1, 1, 1, 1)+real_delta
   call require_bvk_model_mismatch(exchange, cache, "offdiag_scale", error)
   exchange%offdiag_scale(1, 1, 1, 1) = &
      & cache%bvk_model_offdiag_scale(1, 1, 1, 1)
   if (allocated(error)) return
   exchange%rad(1, 1) = exchange%rad(1, 1)+real_delta
   call require_bvk_model_mismatch(exchange, cache, "rad", error)
   exchange%rad(1, 1) = cache%bvk_model_rad(1, 1)
   if (allocated(error)) return
   exchange%kq(1, 1) = exchange%kq(1, 1)+real_delta
   call require_bvk_model_mismatch(exchange, cache, "kq", error)
   exchange%kq(1, 1) = cache%bvk_model_kq(1, 1)
   if (allocated(error)) return
   exchange%corr_scale(1, 1) = exchange%corr_scale(1, 1)+real_delta
   call require_bvk_model_mismatch(exchange, cache, "corr_scale", error)
   exchange%corr_scale(1, 1) = cache%bvk_model_corr_scale(1, 1)
   if (allocated(error)) return
   exchange%corr_rad(1, 1) = exchange%corr_rad(1, 1)+real_delta
   call require_bvk_model_mismatch(exchange, cache, "corr_rad", error)
   exchange%corr_rad(1, 1) = cache%bvk_model_corr_rad(1, 1)
   if (allocated(error)) return

   if (.not.cache%model_matches(exchange)) then
      call test_failed(error, "Restored BvK model signature does not match")
   end if
end subroutine check_bvk_model_signature

subroutine test_bvk_exchange_kernel(error)
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: num(2) = [1, 1]
   integer, parameter :: num_sc(6) = [1, 1, 1, 1, 1, 1]
   integer, parameter :: nmesh(3) = [3, 1, 1]
   integer, parameter :: nmesh_large(3) = [9, 9, 1]
   integer, parameter :: nmesh_cube(3) = [2, 2, 2]
   integer, parameter :: nmesh_overflow(3) = [65536, 65537, 1]
   integer, parameter :: nk_overflow = 65536
   real(wp), parameter :: xyz(3, 2) = reshape([ &
      & 0.40_wp, 0.30_wp, 0.20_wp, &
      & 1.35_wp, 0.35_wp, 0.25_wp], [3, 2])/autoaa
   real(wp), parameter :: lattice(3, 3) = reshape([ &
      & 3.20_wp, 0.00_wp, 0.00_wp, &
      & 0.15_wp, 20.0_wp, 0.00_wp, &
      & 0.05_wp, 0.10_wp, 21.0_wp], [3, 3])/autoaa
   logical, parameter :: periodic(3) = [.true., .true., .true.]

   type(structure_type) :: mol, mol_probe, mol_sc
   class(qvszp_basis_type), allocatable :: qbas, qbas_sc
   class(basis_type), allocatable :: bas, bas_sc
   class(exchange_type), allocatable :: exchange, exchange_sc
   type(basis_cache) :: bcache, bcache_sc
   type(container_cache) :: cache, cache_cp2k, cache_sc
   ! Keep the complete-mesh oracle independent of the production cached plan.
   type(cp2k_stream_plan_type) :: cache_cp2k_plan, oracle_plan
   type(exchange_cache), pointer :: ptr, ptr_cp2k, ptr_sc
   type(exchange_bvk_kernel) :: kernel, kernel_large, kernel_one, &
      & kernel_order_probe
   type(cp2k_exchange_stream_type) :: stream, stream_large, stream_order, &
      & stream_permuted, stream_reverse, stream_reverse_missing, &
      & stream_reverse_reduced, stream_reverse_batch, &
      & stream_forward_missing, stream_batch, stream_batch_large
   type(cp2k_exchange_stream_type) :: stream_overflow
   type(xtb_calculator) :: calc_cp2k
   type(wavefunction_type) :: wfn, wfn_aux, wfn_large, wfn_sc, wfn_aux_sc
   integer, parameter :: korder(3) = [3, 1, 2]
   integer :: acell, bcell, bvk_builds, i, iat, ibatch_case, icell, ik, ir, &
      & ish, j, jat, jsh
   integer :: first, first_image, ii, izp, jao, jj, jzp, last, last_image, &
      & nao, ncell, ni, nj, nsh, requested_batch_size
   integer(int64) :: peak_expected, stream_compact_fft_peak, &
      & stream_compact_separable_peak
   logical :: batches_done
   real(wp) :: angle, deriv_fd, deriv_ref, energy_large, &
      & energy_large_oracle, energy_mesh, energy_m, energy_p, energy_order, &
      & energy_stream, energy_stream_order, energy_tile, energy_tile_sum, &
      & energy_twist, model_probe, &
      & fft_err, onsite_probe, order_err, separable_err, &
      & selector_fft_err, selector_separable_err, theta_x, theta_y
   real(wp) :: energy_sc, err, kfrac(3, 3), kfrac_bad(3, 3), kfrac_order(3, 3), &
      & kfrac_twist(3, 3), lattice_sc(3, 3), step, twist_large(3)
   real(wp) :: partial_err, stream_compact_fft_err, &
      & stream_compact_separable_err, stream_reverse_err
   real(wp) :: rks_fft_forward(3), rks_fft_reverse(3), &
      & rks_2d_fft_forward(3), rks_2d_fft_reverse(3), &
      & uks_3d_fft_forward(3), uks_3d_fft_reverse(3)
   real(wp) :: weights(3), weights_bad(3), weights_order(3), xyz_sc(3, 6), shift(3)
   real(wp) :: coord_direction(3, 2), identity(3, 3), lattice0(3, 3), &
      & strain_direction(3, 3), transform(3, 3), xyz0(3, 2)
   complex(wp) :: batch_sentinel, block_probe, phase
   real(wp), allocatable :: density_sc(:, :, :), overlap_sc(:, :)
   real(wp), allocatable :: ao_grad_sc(:, :), bocorr_grad_r(:, :, :), &
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
   complex(wp), allocatable :: density_mesh(:, :, :, :), &
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

   call new(mol, num, xyz, lattice=lattice, periodic=periodic)
   fft_err = 0.0_wp
   separable_err = 0.0_wp
   stream_reverse_err = 0.0_wp
   allocate(qbas)
   call make_qvszp_basis(qbas, mol, error)
   if (allocated(error)) return
   call move_alloc(qbas, bas)
   call new_wavefunction(wfn_aux, mol%nat, 0, 0, 1, 0.0_wp, .false.)
   call get_eeqbc_charges(mol, error, wfn_aux%qat(:, 1))
   if (allocated(error)) return
   call bas%update(mol, bcache, .false., wfn_aux=wfn_aux)
   call make_exchange_gxtb(exchange, mol, bas)
   call exchange%update(mol, cache)
   call view(cache, ptr)

   select type (fock => exchange)
   type is (exchange_fock)
      call fock%get_bvk_Kmatrix(mol, [1, 1, 1], kernel_one)
      call fock%get_bvk_Kmatrix(mol, nmesh, kernel)
   class default
      call test_failed(error, "Unexpected exchange implementation")
      return
   end select

   err = max(maxval(abs(kernel_one%g_mulliken_r(:, :, 1) &
      & -ptr%g_mulliken)), maxval(abs(kernel_one%g_bocorr_r(:, :, 1) &
      & -ptr%g_bocorr)))

   lattice_sc = lattice
   lattice_sc(:, 1) = real(nmesh(1), wp)*lattice(:, 1)
   do icell = 1, nmesh(1)
      shift = real(icell-1, wp)*lattice(:, 1)
      xyz_sc(:, 2*icell-1:2*icell) = xyz + spread(shift, 2, 2)
   end do
   call new(mol_sc, num_sc, xyz_sc, lattice=lattice_sc, periodic=periodic)
   allocate(qbas_sc)
   call make_qvszp_basis(qbas_sc, mol_sc, error)
   if (allocated(error)) return
   call move_alloc(qbas_sc, bas_sc)
   call new_wavefunction(wfn_aux_sc, mol_sc%nat, 0, 0, 1, 0.0_wp, .false.)
   call get_eeqbc_charges(mol_sc, error, wfn_aux_sc%qat(:, 1))
   if (allocated(error)) return
   call bas_sc%update(mol_sc, bcache_sc, .false., wfn_aux=wfn_aux_sc)
   call make_exchange_gxtb(exchange_sc, mol_sc, bas_sc)
   call exchange_sc%update(mol_sc, cache_sc)
   call view(cache_sc, ptr_sc)

   nsh = bas%nsh
   if (bas_sc%nsh /= nmesh(1)*nsh) then
      call test_failed(error, "Unexpected replicated basis layout")
      return
   end if
   do icell = 1, product(nmesh)
      first = (icell-1)*nsh + 1
      last = icell*nsh
      err = max(err, maxval(abs(kernel%g_mulliken_r(:, :, icell) &
         & -ptr_sc%g_mulliken(1:nsh, first:last))), &
         & maxval(abs(kernel%g_bocorr_r(:, :, icell) &
         & -ptr_sc%g_bocorr(1:2, 2*icell-1:2*icell))))
   end do
   if (any(kernel%reps /= reshape([0, 0, 0, 1, 0, 0, 2, 0, 0], &
      & shape(kernel%reps)))) err = max(err, 1.0_wp)

   call check(error, err, 0.0_wp, thr=1.0e-12_wp)
   if (allocated(error)) then
      call test_failed(error, &
         & "Image-resolved exchange kernels do not reproduce the BvK supercell")
      return
   end if

   ! Construct deliberately non-symmetric real image blocks.  R=2 is the
   ! transpose of R=1, so their Fourier transforms are genuinely complex
   ! Hermitian matrices and the assembled supercell is real symmetric.
   nao = bas%nao
   ncell = product(nmesh)
   if (bas_sc%nao /= ncell*nao) then
      call test_failed(error, "Unexpected replicated AO layout")
      return
   end if
   allocate(p_r(nao, nao, ncell), s_r(nao, nao, ncell), source=0.0_wp)
   do i = 1, nao
      s_r(i, i, 1) = 1.0_wp + 0.04_wp*real(i, wp)
      p_r(i, i, 1) = 0.55_wp + 0.07_wp*real(i, wp)
      do j = 1, i-1
         s_r(j, i, 1) = 0.03_wp/real(i+j, wp)
         s_r(i, j, 1) = s_r(j, i, 1)
         p_r(j, i, 1) = 0.04_wp*cos(real(2*i+3*j, wp))
         p_r(i, j, 1) = p_r(j, i, 1)
      end do
      do j = 1, nao
         s_r(i, j, 2) = 0.012_wp*sin(real(3*i+5*j, wp))
         p_r(i, j, 2) = 0.018_wp*cos(real(7*i+2*j, wp))
      end do
   end do
   s_r(:, :, 3) = transpose(s_r(:, :, 2))
   p_r(:, :, 3) = transpose(p_r(:, :, 2))

   allocate(overlap_sc(ncell*nao, ncell*nao), &
      & density_sc(ncell*nao, ncell*nao, 1), source=0.0_wp)
   do acell = 0, ncell-1
      do bcell = 0, ncell-1
         ir = modulo(bcell-acell, ncell) + 1
         overlap_sc(acell*nao+1:(acell+1)*nao, &
            & bcell*nao+1:(bcell+1)*nao) = s_r(:, :, ir)
         density_sc(acell*nao+1:(acell+1)*nao, &
            & bcell*nao+1:(bcell+1)*nao, 1) = p_r(:, :, ir)
      end do
   end do

   kfrac = 0.0_wp
   weights = 1.0_wp/real(ncell, wp)
   do ik = 1, ncell
      kfrac(1, ik) = real(ik-1, wp)/real(ncell, wp)
   end do
   allocate(overlap_mesh(nao, nao, ncell), &
      & density_mesh(nao, nao, 1, ncell), &
      & fock_mesh(nao, nao, 1, ncell), &
      & fock_ref(nao, nao, 1, ncell), source=(0.0_wp, 0.0_wp))
   allocate(vsh_mesh(nsh), vsh_ref(nsh), source=0.0_wp)
   do ik = 1, ncell
      do ir = 1, ncell
         angle = 2.0_wp*pi*kfrac(1, ik)*real(ir-1, wp)
         phase = exp(cmplx(0.0_wp, angle, wp))
         overlap_mesh(:, :, ik) = overlap_mesh(:, :, ik) &
            & + phase*s_r(:, :, ir)
         density_mesh(:, :, 1, ik) = density_mesh(:, :, 1, ik) &
            & + phase*p_r(:, :, ir)
      end do
   end do

   call new_wavefunction(wfn, mol%nat, nsh, nao, 1, 0.0_wp, .false.)
   do i = 1, nsh
      wfn%qsh(i, 1) = 0.08_wp*sin(real(3*i, wp))
   end do
   call exchange%get_onsite_Kmatrix(mol, wfn, ptr)

   call new_wavefunction(wfn_sc, mol_sc%nat, bas_sc%nsh, bas_sc%nao, &
      & 1, 0.0_wp, .false.)
   do icell = 1, ncell
      first = (icell-1)*nsh + 1
      last = icell*nsh
      wfn_sc%qsh(first:last, 1) = wfn%qsh(:, 1)
   end do
   call exchange_sc%get_onsite_Kmatrix(mol_sc, wfn_sc, ptr_sc)
   allocate(ptr_sc%prev_F(bas_sc%nao, bas_sc%nao, 1), &
      & ptr_sc%prev_vsh(bas_sc%nsh, 1))
   call exchange_sc%get_KFock(mol_sc, ptr_sc, density_sc, overlap_sc)

   select type (fock => exchange)
   type is (exchange_fock)
      call get_KFock_kmesh(fock, mol, ptr, oracle_plan, kernel, kfrac, weights, &
         & density_mesh, overlap_mesh, fock_mesh, vsh_mesh, energy_mesh)
   class default
      call test_failed(error, "Unexpected exchange implementation")
      return
   end select

   do ik = 1, ncell
      do ir = 1, ncell
         first = (ir-1)*nao + 1
         last = ir*nao
         angle = 2.0_wp*pi*kfrac(1, ik)*real(ir-1, wp)
         phase = exp(cmplx(0.0_wp, angle, wp))
         fock_ref(:, :, 1, ik) = fock_ref(:, :, 1, ik) &
            & + phase*ptr_sc%prev_F(1:nao, first:last, 1)
      end do
   end do
   energy_sc = 0.5_wp*sum(density_sc*ptr_sc%prev_F)/real(ncell, wp)
   vsh_ref = 0.0_wp
   do icell = 1, ncell
      first = (icell-1)*nsh + 1
      last = icell*nsh
      vsh_ref = vsh_ref + ptr_sc%prev_vsh(first:last, 1)/real(ncell, wp)
   end do

   err = max(maxval(abs(fock_mesh-fock_ref)), abs(energy_mesh-energy_sc), &
      & maxval(abs(vsh_mesh-vsh_ref)))
   do icell = 1, ncell
      first = (icell-1)*nsh + 1
      last = icell*nsh
      err = max(err, maxval(abs(vsh_mesh-ptr_sc%prev_vsh(first:last, 1))))
   end do

   ! A perturbation at one non-Gamma block must change the other blocks via
   ! the image kernels, while the BZ energy derivative remains the Fock map.
   allocate(direction(nao, nao, 1, ncell), &
      & density_m(nao, nao, 1, ncell), density_p(nao, nao, 1, ncell), &
      & fock_m(nao, nao, 1, ncell), fock_p(nao, nao, 1, ncell), &
      & source=(0.0_wp, 0.0_wp))
   do i = 1, nao
      direction(i, i, 1, 2) = cmplx(0.01_wp*sin(real(5*i, wp)), 0.0_wp, wp)
      do j = 1, i-1
         direction(j, i, 1, 2) = cmplx( &
            & 0.01_wp*cos(real(3*i+j, wp)), &
            & 0.01_wp*sin(real(i+4*j, wp)), wp)
         direction(i, j, 1, 2) = conjg(direction(j, i, 1, 2))
      end do
   end do
   step = 1.0e-5_wp
   density_m = density_mesh-step*direction
   density_p = density_mesh+step*direction
   select type (fock => exchange)
   type is (exchange_fock)
      call get_KFock_kmesh(fock, mol, ptr, oracle_plan, kernel, kfrac, weights, &
         & density_m, overlap_mesh, fock_m, vsh_mesh, energy_m)
      call get_KFock_kmesh(fock, mol, ptr, oracle_plan, kernel, kfrac, weights, &
         & density_p, overlap_mesh, fock_p, vsh_mesh, energy_p)
   end select
   deriv_fd = (energy_p-energy_m)/(2.0_wp*step)
   deriv_ref = 0.0_wp
   do ik = 1, ncell
      deriv_ref = deriv_ref + weights(ik)*real(sum( &
         & conjg(direction(:, :, :, ik))*fock_mesh(:, :, :, ik)), wp)
   end do
   err = max(err, abs(deriv_fd-deriv_ref))
   if (maxval(abs(fock_p(:, :, :, 1)-fock_m(:, :, :, 1))) &
      & < 1.0e-12_wp) err = max(err, 1.0_wp)

   ! Reverse the complete-mesh functional with respect to independently
   ! Hermitian S(k) blocks.  The public response is unweighted, so the BZ
   ! contraction below contains exactly one factor of w_k.
   allocate(overlap_grad_mesh(nao, nao, ncell), &
      & overlap_direction(nao, nao, ncell), overlap_m(nao, nao, ncell), &
      & overlap_p(nao, nao, ncell), source=(0.0_wp, 0.0_wp))
   allocate(mulliken_grad_r(nao, nao, ncell), &
      & bocorr_grad_r(nao, nao, ncell), source=0.0_wp)
   do ik = 1, ncell
      do i = 1, nao
         overlap_direction(i, i, ik) = cmplx( &
            & 0.01_wp*sin(real(5*i+7*ik, wp)), 0.0_wp, wp)
         do j = 1, i-1
            overlap_direction(j, i, ik) = cmplx( &
               & 0.01_wp*cos(real(3*i+j+2*ik, wp)), &
               & 0.01_wp*sin(real(i+4*j+3*ik, wp)), wp)
            overlap_direction(i, j, ik) = &
               & conjg(overlap_direction(j, i, ik))
         end do
      end do
   end do
   select type (fock => exchange)
   type is (exchange_fock)
      call get_KGrad_kmesh(fock, mol, ptr, oracle_plan, kernel, kfrac, weights, &
         & density_mesh, overlap_mesh, overlap_grad_mesh, &
         & mulliken_grad_r, bocorr_grad_r)
      overlap_m = overlap_mesh-step*overlap_direction
      overlap_p = overlap_mesh+step*overlap_direction
      call get_KFock_kmesh(fock, mol, ptr, oracle_plan, kernel, kfrac, weights, &
         & density_mesh, overlap_m, fock_m, vsh_mesh, energy_m)
      call get_KFock_kmesh(fock, mol, ptr, oracle_plan, kernel, kfrac, weights, &
         & density_mesh, overlap_p, fock_p, vsh_mesh, energy_p)
   end select
   deriv_fd = (energy_p-energy_m)/(2.0_wp*step)
   deriv_ref = 0.0_wp
   do ik = 1, ncell
      deriv_ref = deriv_ref + weights(ik)*real(sum( &
         & conjg(overlap_grad_mesh(:, :, ik)) &
         & *overlap_direction(:, :, ik)), wp)
      err = max(err, maxval(abs(overlap_grad_mesh(:, :, ik) &
         & -conjg(transpose(overlap_grad_mesh(:, :, ik))))))
   end do
   err = max(err, abs(deriv_fd-deriv_ref))

   ! Perturb every oriented image kernel independently.  In particular no
   ! R/-R or row/column symmetrization is allowed in this contraction.
   allocate(kernel_mulliken_backup(nsh, nsh, ncell), &
      & kernel_mulliken_direction(nsh, nsh, ncell), &
      & kernel_bocorr_backup(mol%nat, mol%nat, ncell), &
      & kernel_bocorr_direction(mol%nat, mol%nat, ncell))
   kernel_mulliken_backup = kernel%g_mulliken_r
   kernel_bocorr_backup = kernel%g_bocorr_r
   do icell = 1, ncell
      do ish = 1, nsh
         do jsh = 1, nsh
            kernel_mulliken_direction(ish, jsh, icell) = &
               & 0.01_wp*sin(real(11*ish+7*jsh+5*icell, wp))
         end do
      end do
      do iat = 1, mol%nat
         do jat = 1, mol%nat
            kernel_bocorr_direction(iat, jat, icell) = &
               & 0.01_wp*cos(real(3*iat+13*jat+2*icell, wp))
         end do
      end do
   end do
   kernel%g_mulliken_r = kernel_mulliken_backup &
      & - step*kernel_mulliken_direction
   kernel%g_bocorr_r = kernel_bocorr_backup-step*kernel_bocorr_direction
   select type (fock => exchange)
   type is (exchange_fock)
      call get_KFock_kmesh(fock, mol, ptr, oracle_plan, kernel, kfrac, weights, &
         & density_mesh, overlap_mesh, fock_m, vsh_mesh, energy_m)
   end select
   kernel%g_mulliken_r = kernel_mulliken_backup &
      & + step*kernel_mulliken_direction
   kernel%g_bocorr_r = kernel_bocorr_backup+step*kernel_bocorr_direction
   select type (fock => exchange)
   type is (exchange_fock)
      call get_KFock_kmesh(fock, mol, ptr, oracle_plan, kernel, kfrac, weights, &
         & density_mesh, overlap_mesh, fock_p, vsh_mesh, energy_p)
   end select
   kernel%g_mulliken_r = kernel_mulliken_backup
   kernel%g_bocorr_r = kernel_bocorr_backup
   deriv_fd = (energy_p-energy_m)/(2.0_wp*step)
   deriv_ref = 0.0_wp
   do icell = 1, ncell
      do ish = 1, nsh
         ii = bas%iao_sh(ish)
         ni = bas%nao_sh(ish)
         do jsh = 1, nsh
            jj = bas%iao_sh(jsh)
            nj = bas%nao_sh(jsh)
            deriv_ref = deriv_ref &
               & + kernel_mulliken_direction(ish, jsh, icell) &
               & *sum(mulliken_grad_r(ii+1:ii+ni, jj+1:jj+nj, icell))
         end do
      end do
      do iat = 1, mol%nat
         izp = mol%id(iat)
         do jat = 1, mol%nat
            jzp = mol%id(jat)
            do ish = 1, bas%nsh_id(izp)
               ii = bas%iao_sh(bas%ish_at(iat)+ish)
               ni = bas%nao_sh(bas%ish_at(iat)+ish)
               do jsh = 1, bas%nsh_id(jzp)
                  jj = bas%iao_sh(bas%ish_at(jat)+jsh)
                  nj = bas%nao_sh(bas%ish_at(jat)+jsh)
                  deriv_ref = deriv_ref &
                     & + kernel_bocorr_direction(iat, jat, icell) &
                     & *sum(bocorr_grad_r(ii+1:ii+ni, jj+1:jj+nj, icell))
               end do
            end do
         end do
      end do
   end do
   err = max(err, abs(deriv_fd-deriv_ref))

   ! Contract the oriented image responses with the very same BvK
   ! Wigner-Seitz literals used in the forward kernel construction.
   allocate(gradient_mesh(3, mol%nat), sigma_mesh(3, 3), source=0.0_wp)
   select type (fock => exchange)
   type is (exchange_fock)
      call fock%get_bvk_Kmatrix_derivs(mol, kernel, mulliken_grad_r, &
         & bocorr_grad_r, gradient_mesh, sigma_mesh)
   end select
   do iat = 1, mol%nat
      do i = 1, 3
         coord_direction(i, iat) = 0.03_wp*sin(real(7*i+11*iat, wp))
      end do
   end do
   xyz0 = mol%xyz
   mol%xyz = xyz0-step*coord_direction
   select type (fock => exchange)
   type is (exchange_fock)
      call fock%get_bvk_Kmatrix(mol, nmesh, kernel)
      call get_KFock_kmesh(fock, mol, ptr, oracle_plan, kernel, kfrac, weights, &
         & density_mesh, overlap_mesh, fock_m, vsh_mesh, energy_m)
   end select
   mol%xyz = xyz0+step*coord_direction
   select type (fock => exchange)
   type is (exchange_fock)
      call fock%get_bvk_Kmatrix(mol, nmesh, kernel)
      call get_KFock_kmesh(fock, mol, ptr, oracle_plan, kernel, kfrac, weights, &
         & density_mesh, overlap_mesh, fock_p, vsh_mesh, energy_p)
   end select
   deriv_fd = (energy_p-energy_m)/(2.0_wp*step)
   deriv_ref = sum(gradient_mesh*coord_direction)
   err = max(err, abs(deriv_fd-deriv_ref))

   identity = 0.0_wp
   do i = 1, 3
      identity(i, i) = 1.0_wp
   end do
   strain_direction = reshape([ &
      & 0.07_wp, 0.02_wp, -0.01_wp, &
      & 0.02_wp, -0.04_wp, 0.03_wp, &
      & -0.01_wp, 0.03_wp, 0.05_wp], [3, 3])
   lattice0 = mol%lattice
   transform = identity-step*strain_direction
   mol%xyz = matmul(transform, xyz0)
   mol%lattice = matmul(transform, lattice0)
   select type (fock => exchange)
   type is (exchange_fock)
      call fock%get_bvk_Kmatrix(mol, nmesh, kernel)
      call get_KFock_kmesh(fock, mol, ptr, oracle_plan, kernel, kfrac, weights, &
         & density_mesh, overlap_mesh, fock_m, vsh_mesh, energy_m)
   end select
   transform = identity+step*strain_direction
   mol%xyz = matmul(transform, xyz0)
   mol%lattice = matmul(transform, lattice0)
   select type (fock => exchange)
   type is (exchange_fock)
      call fock%get_bvk_Kmatrix(mol, nmesh, kernel)
      call get_KFock_kmesh(fock, mol, ptr, oracle_plan, kernel, kfrac, weights, &
         & density_mesh, overlap_mesh, fock_p, vsh_mesh, energy_p)
   end select
   deriv_fd = (energy_p-energy_m)/(2.0_wp*step)
   deriv_ref = sum(sigma_mesh*strain_direction)
   err = max(err, abs(deriv_fd-deriv_ref))
   mol%xyz = xyz0
   mol%lattice = lattice0
   select type (fock => exchange)
   type is (exchange_fock)
      call fock%get_bvk_Kmatrix(mol, nmesh, kernel)
   end select

   ! The real Gamma exchange gradient of the explicit 3x1x1 BvK supercell is
   ! an independent oracle for the image-kernel geometry contraction.
   allocate(mulliken_grad_sc(bas_sc%nao, bas_sc%nao), &
      & bocorr_grad_sc(bas_sc%nao, bas_sc%nao), &
      & ao_grad_sc(bas_sc%nao, bas_sc%nao), &
      & gradient_sc(3, mol_sc%nat), sigma_sc(3, 3), source=0.0_wp)
   select type (fock_sc => exchange_sc)
   type is (exchange_fock)
      call fock_sc%get_KGrad(mol_sc, ptr_sc, density_sc, overlap_sc, &
         & mulliken_grad_sc, bocorr_grad_sc, ao_grad_sc)
      call fock_sc%get_mulliken_derivs(mol_sc, ptr_sc, mulliken_grad_sc, &
         & gradient_sc, sigma_sc)
      call fock_sc%get_bocorr_derivs(mol_sc, ptr_sc, bocorr_grad_sc, &
         & gradient_sc, sigma_sc)
   end select
   do icell = 0, ncell-1
      err = max(err, maxval(abs(gradient_mesh &
         & -gradient_sc(:, 2*icell+1:2*icell+2))))
   end do
   err = max(err, maxval(abs(sigma_mesh-sigma_sc/real(ncell, wp))))

   ! A common Monkhorst-Pack twist must retain the variational Fock map, and
   ! neither the result nor the energy may depend on CP2K's k-point ordering.
   allocate(density_order(nao, nao, 1, ncell), &
      & fock_order(nao, nao, 1, ncell), fock_twist(nao, nao, 1, ncell), &
      & overlap_order(nao, nao, ncell), source=(0.0_wp, 0.0_wp))
   allocate(vsh_order(nsh), vsh_twist(nsh), source=0.0_wp)
   do ik = 1, ncell
      kfrac_twist(:, ik) = kfrac(:, ik) + [0.137_wp, -0.073_wp, 0.041_wp]
   end do
   do ik = 1, ncell
      kfrac_order(:, ik) = kfrac_twist(:, korder(ik))
      weights_order(ik) = weights(korder(ik))
      density_order(:, :, :, ik) = density_mesh(:, :, :, korder(ik))
      overlap_order(:, :, ik) = overlap_mesh(:, :, korder(ik))
   end do
   select type (fock => exchange)
   type is (exchange_fock)
      call get_KFock_kmesh(fock, mol, ptr, oracle_plan, kernel, kfrac_twist, weights, &
         & density_mesh, overlap_mesh, fock_twist, vsh_twist, energy_twist)
      call get_KFock_kmesh(fock, mol, ptr, oracle_plan, kernel, kfrac_order, weights_order, &
         & density_order, overlap_order, fock_order, vsh_order, energy_order)
      call get_KFock_kmesh(fock, mol, ptr, oracle_plan, kernel, kfrac_twist, weights, &
         & density_m, overlap_mesh, fock_m, vsh_mesh, energy_m)
      call get_KFock_kmesh(fock, mol, ptr, oracle_plan, kernel, kfrac_twist, weights, &
         & density_p, overlap_mesh, fock_p, vsh_mesh, energy_p)
   end select
   deriv_fd = (energy_p-energy_m)/(2.0_wp*step)
   deriv_ref = 0.0_wp
   order_err = 0.0_wp
   do ik = 1, ncell
      deriv_ref = deriv_ref + weights(ik)*real(sum( &
         & conjg(direction(:, :, :, ik))*fock_twist(:, :, :, ik)), wp)
      order_err = max(order_err, maxval(abs(fock_order(:, :, :, ik) &
         & -fock_twist(:, :, :, korder(ik)))))
   end do
   err = max(err, abs(deriv_fd-deriv_ref), abs(energy_order-energy_twist), &
      & maxval(abs(vsh_order-vsh_twist)), order_err)

   ! Exercise the public CP2K whole-mesh path and use its persistent plan as
   ! the cache-hit/invalidation oracle.  The first call builds the plan, an
   ! identical second call must reuse it, and a reordered mesh must rebuild
   ! the ordering descriptor without changing the physical result.
   allocate(calc_cp2k%exchange, source=exchange)

   ! On 32-bit default integers this deliberately aliases Nk when evaluated
   ! by unchecked PRODUCT: 65536*65537 wraps to 65536.  The public stream (and
   ! hence partial-accumulator) boundary must reject the true oversized mesh
   ! before preparing a plan or entering any default-integer index arithmetic.
   allocate(kfrac_overflow(3, nk_overflow), source=0.0_wp)
   allocate(weights_overflow(nk_overflow), &
      & source=1.0_wp/real(nk_overflow, wp))
   call cp2k_exchange_stream_begin(stream_overflow, calc_cp2k, mol, &
      & cache_cp2k, cache_cp2k_plan, wfn, nmesh_overflow, kfrac_overflow, weights_overflow, &
      & error)
   if (allocated(error)) then
      if (trim(error%message) /= &
         & "g-xTB exchange stream requires a complete regular k mesh") then
         deallocate(error)
         call test_failed(error, "Unexpected overflowing-mesh diagnostic")
         return
      end if
      deallocate(error)
   else
      call test_failed(error, "Overflowing g-xTB exchange mesh was accepted")
      return
   end if
   deallocate(kfrac_overflow, weights_overflow)

   ! Range validation must not rely on Fortran logical short-circuiting: the
   ! remaining-count expression is evaluated only after first is in 1:Nk.
   call cp2k_exchange_stream_begin(stream_overflow, calc_cp2k, mol, &
      & cache_cp2k, cache_cp2k_plan, wfn, nmesh, kfrac, weights, error, &
      & mode=cp2k_exchange_stream_batched, batch_size=1, &
      & image_first=-huge(0), image_count=1)
   if (allocated(error)) then
      if (trim(error%message) /= &
         & "g-xTB exchange image range is invalid") then
         deallocate(error)
         call test_failed(error, "Unexpected extreme image-range diagnostic")
         return
      end if
      deallocate(error)
   else
      call test_failed(error, "Extreme g-xTB image range was accepted")
      return
   end if

   allocate(fock_cp2k(nao, nao, 1, ncell), vsh_cp2k(nsh))
   call cp2k_exchange_kmesh(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn, nmesh, &
      & kfrac, weights, density_mesh, overlap_mesh, fock_cp2k, &
      & vsh_cp2k, energy_mesh, error)
   if (allocated(error)) return
   call view(cache_cp2k, ptr_cp2k)
   if (.not.cache_cp2k_plan%bvk_matches(mol, nmesh, kfrac, weights)) &
      & err = max(err, 1.0_wp)
   if (any(cache_cp2k_plan%bvk_input_to_grid /= [1, 2, 3])) &
      & err = max(err, 1.0_wp)
   if (any(cache_cp2k_plan%bvk_grid_to_input /= [1, 2, 3])) &
      & err = max(err, 1.0_wp)
   err = max(err, maxval(abs(fock_cp2k-fock_ref)), &
      & maxval(abs(vsh_cp2k-vsh_ref)), abs(energy_mesh-energy_sc))

   ! The memory-reduced stream must reject duplicate and missing input
   ! blocks without corrupting its transaction.  Completing that transaction
   ! afterwards must reproduce the unchanged whole-mesh oracle exactly.
   allocate(fock_stream(nao, nao, 1, ncell), &
      & fock_stream_order(nao, nao, 1, ncell), &
      & vsh_stream(nsh), vsh_stream_order(nsh), &
      & gradient_stream_probe(3, mol%nat), sigma_stream_probe(3, 3))
   call cp2k_exchange_stream_begin(stream, calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn, &
      & nmesh, kfrac, weights, error)
   if (allocated(error)) return
   if (.not.cp2k_exchange_stream_has_full_mesh_storage(stream)) &
      & err = max(err, 1.0_wp)
   call cp2k_exchange_stream_push(stream, 1, density_mesh(:, :, :, 1), &
      & overlap_mesh(:, :, 1), error)
   if (allocated(error)) return
   call cp2k_exchange_stream_push(stream, 1, density_mesh(:, :, :, 1), &
      & overlap_mesh(:, :, 1), error)
   if (allocated(error)) then
      if (trim(error%message) /= "duplicate g-xTB exchange stream block") then
         deallocate(error)
         call test_failed(error, "Unexpected duplicate stream-block diagnostic")
         return
      end if
      deallocate(error)
   else
      call test_failed(error, "Duplicate exchange stream block was accepted")
      return
   end if
   call cp2k_exchange_stream_push(stream, 2, density_mesh(:, :, :, 2), &
      & overlap_mesh(:, :, 2), error)
   if (allocated(error)) return
   call cp2k_exchange_stream_apply(stream, calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, &
      & vsh_stream, energy_stream, error)
   if (allocated(error)) then
      if (trim(error%message) /= "g-xTB exchange stream has missing blocks") then
         deallocate(error)
         call test_failed(error, "Unexpected missing stream-block diagnostic")
         return
      end if
      deallocate(error)
   else
      call test_failed(error, "Incomplete exchange stream was applied")
      return
   end if
   call cp2k_exchange_stream_push(stream, 3, density_mesh(:, :, :, 3), &
      & overlap_mesh(:, :, 3), error)
   if (allocated(error)) return
   onsite_probe = ptr_cp2k%g_onsfx(1, 1, 1)
   ptr_cp2k%g_onsfx(1, 1, 1) = onsite_probe + 1.0e-8_wp
   call cp2k_exchange_stream_apply(stream, calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, &
      & vsh_stream, energy_stream, error)
   if (allocated(error)) then
      if (trim(error%message) /= &
         & "g-xTB exchange stream state was invalidated before apply") then
         deallocate(error)
         call test_failed(error, "Unexpected invalidated stream diagnostic")
         return
      end if
      deallocate(error)
   else
      call test_failed(error, "Invalidated exchange stream was applied")
      return
   end if
   ptr_cp2k%g_onsfx(1, 1, 1) = onsite_probe
   block_probe = cache_cp2k_plan%bvk_phase_forward(1, 1)
   cache_cp2k_plan%bvk_phase_forward(1, 1) = block_probe+cmplx(0.0_wp, 1.0e-8_wp, wp)
   call cp2k_exchange_stream_apply(stream, calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, &
      & vsh_stream, energy_stream, error)
   if (allocated(error)) then
      if (trim(error%message) /= &
         & "g-xTB exchange stream state was invalidated before apply") then
         deallocate(error)
         call test_failed(error, "Unexpected invalidated phase diagnostic")
         return
      end if
      deallocate(error)
   else
      call test_failed(error, "Invalidated exchange phase was applied")
      return
   end if
   cache_cp2k_plan%bvk_phase_forward(1, 1) = block_probe
   call cp2k_exchange_stream_apply(stream, calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, &
      & vsh_stream, energy_stream, error)
   if (allocated(error)) return
   if (.not.cp2k_exchange_stream_has_full_mesh_storage(stream)) &
      & err = max(err, 1.0_wp)
   do ik = 1, ncell
      call cp2k_exchange_stream_get(stream, ik, &
         & fock_stream(:, :, :, ik), error, overlap_mesh(:, :, ik))
      if (allocated(error)) return
   end do
   call cp2k_exchange_stream_get(stream, 1, &
      & fock_stream(:, :, :, 1), error, overlap_mesh(:, :, 1))
   if (allocated(error)) then
      if (trim(error%message) /= &
         & "duplicate g-xTB exchange stream forward pull") then
         deallocate(error)
         call test_failed(error, "Unexpected duplicate forward-pull diagnostic")
         return
      end if
      deallocate(error)
   else
      call test_failed(error, "Exchange Fock block was pulled twice")
      return
   end if
   call cp2k_exchange_stream_end(stream, error)
   if (allocated(error)) return
   err = max(err, maxval(abs(fock_stream-fock_cp2k)), &
      & maxval(abs(vsh_stream-vsh_cp2k)), abs(energy_stream-energy_mesh))

   ! Closing a forward transaction requires exactly one pull of every block;
   ! the failed close must clear the object so it can be reused immediately.
   call cp2k_exchange_stream_begin(stream_forward_missing, calc_cp2k, mol, &
      & cache_cp2k, cache_cp2k_plan, wfn, nmesh, kfrac, weights, error)
   if (allocated(error)) return
   do ik = 1, ncell
      call cp2k_exchange_stream_push(stream_forward_missing, ik, &
         & density_mesh(:, :, :, ik), overlap_mesh(:, :, ik), error)
      if (allocated(error)) return
   end do
   call cp2k_exchange_stream_apply(stream_forward_missing, calc_cp2k, mol, &
      & cache_cp2k, cache_cp2k_plan, vsh_stream, energy_stream, error)
   if (allocated(error)) return
   do ik = 1, ncell-1
      call cp2k_exchange_stream_get(stream_forward_missing, ik, &
         & fock_stream(:, :, :, ik), error, overlap_mesh(:, :, ik))
      if (allocated(error)) return
   end do
   call cp2k_exchange_stream_end(stream_forward_missing, error)
   if (allocated(error)) then
      if (trim(error%message) /= &
         & "g-xTB exchange stream ended with missing forward pulls") then
         deallocate(error)
         call test_failed(error, "Unexpected missing forward-pull diagnostic")
         return
      end if
      deallocate(error)
   else
      call test_failed(error, "Exchange stream closed with a missing Fock pull")
      return
   end if
   call cp2k_exchange_stream_begin(stream_forward_missing, calc_cp2k, mol, &
      & cache_cp2k, cache_cp2k_plan, wfn, nmesh, kfrac, weights, error)
   if (allocated(error)) return
   do ik = 1, ncell
      call cp2k_exchange_stream_push(stream_forward_missing, ik, &
         & density_mesh(:, :, :, ik), overlap_mesh(:, :, ik), error)
      if (allocated(error)) return
   end do
   call cp2k_exchange_stream_apply(stream_forward_missing, calc_cp2k, mol, &
      & cache_cp2k, cache_cp2k_plan, vsh_stream, energy_stream, error)
   if (allocated(error)) return
   do ik = 1, ncell
      call cp2k_exchange_stream_get(stream_forward_missing, ik, &
         & fock_stream(:, :, :, ik), error, overlap_mesh(:, :, ik))
      if (allocated(error)) return
   end do
   call cp2k_exchange_stream_end(stream_forward_missing, error)
   if (allocated(error)) return

   ! Arrival order is independent of the physical k-point ordering.  Push the
   ! same indexed blocks in a nontrivial order and compare the full response.
   call cp2k_exchange_stream_begin(stream_order, calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, &
      & wfn, nmesh, kfrac, weights, error)
   if (allocated(error)) return
   if (.not.cp2k_exchange_stream_has_full_mesh_storage(stream_order)) &
      & err = max(err, 1.0_wp)
   do j = 1, ncell
      ik = korder(j)
      call cp2k_exchange_stream_push(stream_order, ik, &
         & density_mesh(:, :, :, ik), overlap_mesh(:, :, ik), error)
      if (allocated(error)) return
   end do
   call cp2k_exchange_stream_apply(stream_order, calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, &
      & vsh_stream_order, energy_stream_order, error)
   if (allocated(error)) return
   call cp2k_exchange_stream_reverse_apply(stream_order, calc_cp2k, mol, &
      & cache_cp2k, cache_cp2k_plan, gradient_stream_probe, sigma_stream_probe, error)
   if (allocated(error)) then
      if (trim(error%message) /= &
         & "g-xTB exchange stream reverse requires oracle mode") then
         deallocate(error)
         call test_failed(error, "Unexpected reduced reverse diagnostic")
         return
      end if
      deallocate(error)
   else
      call test_failed(error, "Reduced exchange stream accepted reverse apply")
      return
   end if
   do ik = 1, ncell
      call cp2k_exchange_stream_get(stream_order, ik, &
         & fock_stream_order(:, :, :, ik), error, overlap_mesh(:, :, ik))
      if (allocated(error)) return
   end do
   call cp2k_exchange_stream_end(stream_order, error)
   if (allocated(error)) return
   err = max(err, maxval(abs(fock_stream_order-fock_stream)), &
      & maxval(abs(vsh_stream_order-vsh_stream)), &
      & abs(energy_stream_order-energy_stream))

   ! A common twist and a permutation of the physical k-point list must use
   ! the cached representatives without changing the reduced-stream result.
   call cp2k_exchange_stream_begin(stream_permuted, calc_cp2k, mol, &
      & cache_cp2k, cache_cp2k_plan, wfn, nmesh, kfrac_order, weights_order, error, &
      & mode=cp2k_exchange_stream_reduced)
   if (allocated(error)) return
   if (.not.cp2k_exchange_stream_has_full_mesh_storage(stream_permuted)) &
      & err = max(err, 1.0_wp)
   do ik = 1, ncell
      call cp2k_exchange_stream_push(stream_permuted, ik, &
         & density_order(:, :, :, ik), overlap_order(:, :, ik), error)
      if (allocated(error)) return
   end do
   call cp2k_exchange_stream_apply(stream_permuted, calc_cp2k, mol, &
      & cache_cp2k, cache_cp2k_plan, vsh_stream_order, energy_stream_order, error)
   if (allocated(error)) return
   do ik = 1, ncell
      call cp2k_exchange_stream_get(stream_permuted, ik, &
         & fock_stream_order(:, :, :, ik), error, overlap_order(:, :, ik))
      if (allocated(error)) return
   end do
   call cp2k_exchange_stream_end(stream_permuted, error)
   if (allocated(error)) return
   err = max(err, maxval(abs(fock_stream_order-fock_order)), &
      & maxval(abs(vsh_stream_order-vsh_order)), &
      & abs(energy_stream_order-energy_order))

   ! The image-batched protocol must retain only a caller-selected number of
   ! AO image matrices.  It deliberately requires the sequential caller to
   ! resubmit all k blocks and to pull every k contribution for every batch.
   batch_sentinel = cmplx(0.037_wp, -0.019_wp, wp)
   allocate(fock_batch(nao, nao, 1, ncell), source=batch_sentinel)
   allocate(overlap_bad_dimensions(nao+1, nao), source=(0.0_wp, 0.0_wp))
   call cp2k_exchange_stream_begin(stream_batch, calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, &
      & wfn, nmesh, kfrac_order, weights_order, error, &
      & mode=cp2k_exchange_stream_batched, batch_size=0)
   if (allocated(error)) then
      if (trim(error%message) /= &
         & "g-xTB exchange stream batch size must be positive") then
         deallocate(error)
         call test_failed(error, "Unexpected invalid batch-size diagnostic")
         return
      end if
      deallocate(error)
   else
      call test_failed(error, "Zero exchange image-batch size was accepted")
      return
   end if

   call cp2k_exchange_stream_begin(stream_batch, calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, &
      & wfn, nmesh, kfrac_order, weights_order, error, &
      & mode=cp2k_exchange_stream_batched, batch_size=2, image_first=1)
   if (allocated(error)) then
      if (trim(error%message) /= &
         & "g-xTB exchange image range requires first and count") then
         deallocate(error)
         call test_failed(error, "Unexpected incomplete image-range diagnostic")
         return
      end if
      deallocate(error)
   else
      call test_failed(error, "Incomplete exchange image range was accepted")
      return
   end if

   call cp2k_exchange_stream_begin(stream_batch, calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, &
      & wfn, nmesh, kfrac_order, weights_order, error, &
      & mode=cp2k_exchange_stream_batched, batch_size=2, &
      & image_first=0, image_count=1)
   if (allocated(error)) then
      if (trim(error%message) /= "g-xTB exchange image range is invalid") then
         deallocate(error)
         call test_failed(error, "Unexpected invalid image-range diagnostic")
         return
      end if
      deallocate(error)
   else
      call test_failed(error, "Invalid exchange image range was accepted")
      return
   end if

   kfrac_bad = kfrac_order
   kfrac_bad(1, 1) = kfrac_bad(1, 1)+2.0e-11_wp
   call cp2k_exchange_stream_begin(stream_batch, calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, &
      & wfn, nmesh, kfrac_bad, weights_order, error, &
      & mode=cp2k_exchange_stream_batched, batch_size=2)
   if (allocated(error)) then
      if (trim(error%message) /= &
         & "batched g-xTB exchange requires a machine-regular grid") then
         deallocate(error)
         call test_failed(error, "Unexpected inexact batch-grid diagnostic")
         return
      end if
      deallocate(error)
   else
      call test_failed(error, "Inexact exchange image-batch grid was accepted")
      return
   end if

   weights_bad = weights_order
   weights_bad(1) = weights_bad(1)+5.0e-13_wp
   weights_bad(2) = weights_bad(2)-5.0e-13_wp
   call cp2k_exchange_stream_begin(stream_batch, calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, &
      & wfn, nmesh, kfrac_order, weights_bad, error, &
      & mode=cp2k_exchange_stream_batched, batch_size=2)
   if (allocated(error)) then
      if (trim(error%message) /= &
         & "batched g-xTB exchange requires machine-uniform weights") then
         deallocate(error)
         call test_failed(error, "Unexpected inexact batch-weight diagnostic")
         return
      end if
      deallocate(error)
   else
      call test_failed(error, "Inexact exchange image-batch weights were accepted")
      return
   end if

   call cp2k_exchange_stream_begin(stream_batch, calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, &
      & wfn, nmesh, kfrac_order, weights_order, error, &
      & mode=cp2k_exchange_stream_batched, batch_size=2)
   if (allocated(error)) return
   if (cp2k_exchange_stream_has_full_mesh_storage(stream_batch)) &
      & err = max(err, 1.0_wp)
   if (allocated(cache_cp2k_plan%bvk_phase_forward) &
      & .or. allocated(cache_cp2k_plan%bvk_phase_inverse)) err = max(err, 1.0_wp)
   if (.not.allocated(cache_cp2k_plan%bvk_phase_roots) &
      & .or. .not.allocated(cache_cp2k_plan%bvk_twist_phase)) then
      err = max(err, 1.0_wp)
   else
      if (any(shape(cache_cp2k_plan%bvk_phase_roots) /= [maxval(nmesh), 3]) &
         & .or. size(cache_cp2k_plan%bvk_twist_phase) /= ncell) &
         & err = max(err, 1.0_wp)
   end if
   peak_expected = 4_int64*int(nao, int64)*int(nao, int64) &
      & *2_int64 + 3_int64*int(nao, int64) &
      & +2_int64*int(ncell, int64) &
      & +3_int64*int(maxval(nmesh), int64)
   peak_expected = peak_expected + max( &
      & 2_int64*int(nao, int64)**2, &
      & 5_int64*int(nao, int64)**2+6_int64*int(nao, int64), &
      & 5_int64*int(nao, int64)**2+2_int64)
   if (cp2k_exchange_stream_peak_complex_elements(stream_batch) &
      & /= peak_expected) err = max(err, 1.0_wp)

   call cp2k_exchange_stream_batch_pull(stream_batch, 1, &
      & overlap_order(:, :, 1), fock_batch(:, :, :, 1), error)
   if (allocated(error)) then
      if (trim(error%message) /= &
         & "g-xTB exchange image batch was not applied") then
         deallocate(error)
         call test_failed(error, "Unexpected pre-apply batch-pull diagnostic")
         return
      end if
      deallocate(error)
   else
      call test_failed(error, "Exchange image batch was pulled before apply")
      return
   end if
   call cp2k_exchange_stream_batch_advance(stream_batch, batches_done, error)
   if (allocated(error)) then
      if (trim(error%message) /= &
         & "g-xTB exchange image batch was not applied") then
         deallocate(error)
         call test_failed(error, "Unexpected pre-apply batch-advance diagnostic")
         return
      end if
      deallocate(error)
   else
      call test_failed(error, "Exchange image batch advanced before apply")
      return
   end if
   call cp2k_exchange_stream_batch_result(stream_batch, calc_cp2k, mol, &
      & cache_cp2k, cache_cp2k_plan, vsh_stream_order, energy_stream_order, error)
   if (allocated(error)) then
      if (trim(error%message) /= &
         & "g-xTB exchange image batches are incomplete") then
         deallocate(error)
         call test_failed(error, "Unexpected incomplete batch-result diagnostic")
         return
      end if
      deallocate(error)
   else
      call test_failed(error, "Incomplete exchange image batches returned a result")
      return
   end if

   batches_done = .false.
   do while (.not.batches_done)
      call cp2k_exchange_stream_batch_bounds(stream_batch, first_image, &
         & last_image, error)
      if (allocated(error)) return
      if (last_image-first_image+1 > 2) err = max(err, 1.0_wp)

      if (first_image == 1) then
         call cp2k_exchange_stream_batch_apply(stream_batch, calc_cp2k, mol, &
            & cache_cp2k, cache_cp2k_plan, error)
         if (allocated(error)) then
            if (trim(error%message) /= &
               & "g-xTB exchange image batch has missing blocks") then
               deallocate(error)
               call test_failed(error, &
                  & "Unexpected missing image-batch block diagnostic")
               return
            end if
            deallocate(error)
         else
            call test_failed(error, "Incomplete exchange image batch applied")
            return
         end if
      end if

      if (first_image > 1) then
         ik = korder(1)
         block_probe = density_order(1, 1, 1, ik)
         density_order(1, 1, 1, ik) = block_probe+1.0e-8_wp
         call cp2k_exchange_stream_push(stream_batch, ik, &
            & density_order(:, :, :, ik), overlap_order(:, :, ik), error)
         density_order(1, 1, 1, ik) = block_probe
         if (allocated(error)) then
            if (trim(error%message) /= &
               & "changed g-xTB exchange image-batch density block") then
               deallocate(error)
               call test_failed(error, "Unexpected changed-density diagnostic")
               return
            end if
            deallocate(error)
         else
            call test_failed(error, "Changed image-batch density was accepted")
            return
         end if

         block_probe = overlap_order(1, 1, ik)
         overlap_order(1, 1, ik) = block_probe+1.0e-8_wp
         call cp2k_exchange_stream_push(stream_batch, ik, &
            & density_order(:, :, :, ik), overlap_order(:, :, ik), error)
         overlap_order(1, 1, ik) = block_probe
         if (allocated(error)) then
            if (trim(error%message) /= &
               & "changed g-xTB exchange image-batch overlap block") then
               deallocate(error)
               call test_failed(error, "Unexpected changed-overlap diagnostic")
               return
            end if
            deallocate(error)
         else
            call test_failed(error, "Changed image-batch overlap was accepted")
            return
         end if
      end if

      do j = 1, ncell
         ik = korder(j)
         call cp2k_exchange_stream_push(stream_batch, ik, &
            & density_order(:, :, :, ik), overlap_order(:, :, ik), error)
         if (allocated(error)) return
         if (first_image == 1 .and. j == 1) then
            call cp2k_exchange_stream_push(stream_batch, ik, &
               & density_order(:, :, :, ik), overlap_order(:, :, ik), error)
            if (allocated(error)) then
               if (trim(error%message) /= &
                  & "duplicate g-xTB exchange stream block") then
                  deallocate(error)
                  call test_failed(error, &
                     & "Unexpected duplicate image-batch block diagnostic")
                  return
               end if
               deallocate(error)
            else
               call test_failed(error, &
                  & "Duplicate exchange image-batch block was accepted")
               return
            end if
         end if
      end do

      if (first_image == 1) then
         select type (fock_cp2k_model => calc_cp2k%exchange)
         type is (exchange_fock)
            model_probe = fock_cp2k_model%onecxints(1, 1, 1)
            fock_cp2k_model%onecxints(1, 1, 1) = model_probe + 1.0e-8_wp
         end select
         call cp2k_exchange_stream_batch_apply(stream_batch, calc_cp2k, mol, &
            & cache_cp2k, cache_cp2k_plan, error)
         if (allocated(error)) then
            if (trim(error%message) /= &
               & "g-xTB exchange image batch state was invalidated before apply") then
               deallocate(error)
               call test_failed(error, &
                  & "Unexpected stale image-batch model diagnostic")
               return
            end if
            deallocate(error)
         else
            call test_failed(error, "Stale exchange image batch was applied")
            return
         end if
         select type (fock_cp2k_model => calc_cp2k%exchange)
         type is (exchange_fock)
            fock_cp2k_model%onecxints(1, 1, 1) = model_probe
         end select
         select type (fock_cp2k_model => calc_cp2k%exchange)
         type is (exchange_fock)
            model_probe = fock_cp2k_model%kq(1, 1)
            fock_cp2k_model%kq(1, 1) = model_probe + 1.0e-8_wp
         end select
         call cp2k_exchange_stream_batch_apply(stream_batch, calc_cp2k, mol, &
            & cache_cp2k, cache_cp2k_plan, error)
         if (allocated(error)) then
            if (trim(error%message) /= &
               & "g-xTB exchange image batch state was invalidated before apply") then
               deallocate(error)
               call test_failed(error, &
                  & "Unexpected stale image-batch kq diagnostic")
               return
            end if
            deallocate(error)
         else
            call test_failed(error, "Stale kq exchange image batch was applied")
            return
         end if
         select type (fock_cp2k_model => calc_cp2k%exchange)
         type is (exchange_fock)
            fock_cp2k_model%kq(1, 1) = model_probe
         end select

         ! The newer model separates onsite and intersite exchange scaling.
         ! An active transaction must reject a changed onsite model too.
         select type (fock_cp2k_model => calc_cp2k%exchange)
         type is (exchange_fock)
            model_probe = fock_cp2k_model%onsite_frscale
            fock_cp2k_model%onsite_frscale = model_probe+1.0e-8_wp
         end select
         call cp2k_exchange_stream_batch_apply(stream_batch, calc_cp2k, mol, &
            & cache_cp2k, cache_cp2k_plan, error)
         select type (fock_cp2k_model => calc_cp2k%exchange)
         type is (exchange_fock)
            fock_cp2k_model%onsite_frscale = model_probe
         end select
         if (allocated(error)) then
            if (trim(error%message) /= &
               & "g-xTB exchange image batch state was invalidated before apply") then
               deallocate(error)
               call test_failed(error, "Unexpected stale onsite-model diagnostic")
               return
            end if
            deallocate(error)
         else
            call test_failed(error, "Stale onsite exchange image batch was applied")
            return
         end if

         model_probe = cache_cp2k_plan%bvk_kernel%g_mulliken_r(1, 1, 1)
         cache_cp2k_plan%bvk_kernel%g_mulliken_r(1, 1, 1) = model_probe+1.0e-8_wp
         call cp2k_exchange_stream_batch_apply(stream_batch, calc_cp2k, mol, &
            & cache_cp2k, cache_cp2k_plan, error)
         if (allocated(error)) then
            if (trim(error%message) /= &
               & "g-xTB exchange image batch state was invalidated before apply") then
               deallocate(error)
               call test_failed(error, "Unexpected stale batch-kernel diagnostic")
               return
            end if
            deallocate(error)
         else
            call test_failed(error, "Stale exchange image-batch kernel was applied")
            return
         end if
         cache_cp2k_plan%bvk_kernel%g_mulliken_r(1, 1, 1) = model_probe

         block_probe = cache_cp2k_plan%bvk_phase_roots(1, 1)
         cache_cp2k_plan%bvk_phase_roots(1, 1) = block_probe &
            & +cmplx(0.0_wp, 1.0e-8_wp, wp)
         call cp2k_exchange_stream_batch_apply(stream_batch, calc_cp2k, mol, &
            & cache_cp2k, cache_cp2k_plan, error)
         if (allocated(error)) then
            if (trim(error%message) /= &
               & "g-xTB exchange image batch state was invalidated before apply") then
               deallocate(error)
               call test_failed(error, &
                  & "Unexpected stale compact-phase diagnostic")
               return
            end if
            deallocate(error)
         else
            call test_failed(error, "Stale compact exchange phase was applied")
            return
         end if
         cache_cp2k_plan%bvk_phase_roots(1, 1) = block_probe

         i = cache_cp2k_plan%bvk_input_to_grid(1)
         j = cache_cp2k_plan%bvk_input_to_grid(2)
         cache_cp2k_plan%bvk_input_to_grid(1) = j
         cache_cp2k_plan%bvk_input_to_grid(2) = i
         cache_cp2k_plan%bvk_grid_to_input(j) = 1
         cache_cp2k_plan%bvk_grid_to_input(i) = 2
         call cp2k_exchange_stream_batch_apply(stream_batch, calc_cp2k, mol, &
            & cache_cp2k, cache_cp2k_plan, error)
         if (allocated(error)) then
            if (trim(error%message) /= &
               & "g-xTB exchange image batch state was invalidated before apply") then
               deallocate(error)
               call test_failed(error, "Unexpected stale batch-order diagnostic")
               return
            end if
            deallocate(error)
         else
            call test_failed(error, "Reordered exchange image batch was applied")
            return
         end if
         cache_cp2k_plan%bvk_input_to_grid(1) = i
         cache_cp2k_plan%bvk_input_to_grid(2) = j
         cache_cp2k_plan%bvk_grid_to_input(i) = 1
         cache_cp2k_plan%bvk_grid_to_input(j) = 2
      end if

      call cp2k_exchange_stream_batch_apply(stream_batch, calc_cp2k, mol, &
         & cache_cp2k, cache_cp2k_plan, error)
      if (allocated(error)) return
      if (first_image == 1) then
         call cp2k_exchange_stream_batch_apply(stream_batch, calc_cp2k, mol, &
            & cache_cp2k, cache_cp2k_plan, error)
         if (allocated(error)) then
            if (trim(error%message) /= &
               & "g-xTB exchange image batch was already applied") then
               deallocate(error)
               call test_failed(error, &
                  & "Unexpected duplicate image-batch apply diagnostic")
               return
            end if
            deallocate(error)
         else
            call test_failed(error, "Exchange image batch was applied twice")
            return
         end if
         call cp2k_exchange_stream_batch_pull(stream_batch, korder(1), &
            & overlap_bad_dimensions, fock_batch(:, :, :, korder(1)), error)
         if (allocated(error)) then
            if (trim(error%message) /= &
               & "inconsistent g-xTB exchange image-batch pull dimensions") then
               deallocate(error)
               call test_failed(error, &
                  & "Unexpected image-batch pull-dimension diagnostic")
               return
            end if
            deallocate(error)
         else
            call test_failed(error, &
               & "Wrong-dimensional exchange image-batch pull was accepted")
            return
         end if
      end if
      if (first_image > 1) then
         ik = korder(1)
         block_probe = overlap_order(1, 1, ik)
         overlap_order(1, 1, ik) = block_probe+1.0e-8_wp
         call cp2k_exchange_stream_batch_pull(stream_batch, ik, &
            & overlap_order(:, :, ik), fock_batch(:, :, :, ik), error)
         overlap_order(1, 1, ik) = block_probe
         if (allocated(error)) then
            if (trim(error%message) /= &
               & "changed g-xTB exchange image-batch pull overlap") then
               deallocate(error)
               call test_failed(error, "Unexpected changed-pull diagnostic")
               return
            end if
            deallocate(error)
         else
            call test_failed(error, "Changed image-batch pull was accepted")
            return
         end if
      end if
      call cp2k_exchange_stream_batch_pull(stream_batch, korder(1), &
         & overlap_order(:, :, korder(1)), &
         & fock_batch(:, :, :, korder(1)), error)
      if (allocated(error)) return
      if (first_image == 1) then
         call cp2k_exchange_stream_batch_pull(stream_batch, korder(1), &
            & overlap_order(:, :, korder(1)), &
            & fock_batch(:, :, :, korder(1)), error)
         if (allocated(error)) then
            if (trim(error%message) /= &
               & "duplicate g-xTB exchange image-batch pull") then
               deallocate(error)
               call test_failed(error, &
                  & "Unexpected duplicate image-batch pull diagnostic")
               return
            end if
            deallocate(error)
         else
            call test_failed(error, &
               & "Duplicate exchange image-batch pull was accepted")
            return
         end if
         call cp2k_exchange_stream_batch_advance(stream_batch, batches_done, &
            & error)
         if (allocated(error)) then
            if (trim(error%message) /= &
               & "g-xTB exchange image batch has missing pulls") then
               deallocate(error)
               call test_failed(error, &
                  & "Unexpected missing image-batch pull diagnostic")
               return
            end if
            deallocate(error)
         else
            call test_failed(error, &
               & "Exchange image batch advanced with missing pulls")
            return
         end if
      end if
      do j = 2, ncell
         ik = korder(j)
         call cp2k_exchange_stream_batch_pull(stream_batch, ik, &
            & overlap_order(:, :, ik), &
            & fock_batch(:, :, :, ik), error)
         if (allocated(error)) return
      end do
      call cp2k_exchange_stream_batch_advance(stream_batch, batches_done, error)
      if (allocated(error)) return
   end do
   call cp2k_exchange_stream_batch_bounds(stream_batch, first_image, &
      & last_image, error)
   if (allocated(error)) then
      if (trim(error%message) /= &
         & "g-xTB exchange image batches are complete") then
         deallocate(error)
         call test_failed(error, "Unexpected completed batch-bounds diagnostic")
         return
      end if
      deallocate(error)
   else
      call test_failed(error, "Completed exchange image batches returned bounds")
      return
   end if
   call cp2k_exchange_stream_batch_advance(stream_batch, batches_done, error)
   if (allocated(error)) then
      if (trim(error%message) /= &
         & "g-xTB exchange image batches are complete") then
         deallocate(error)
         call test_failed(error, "Unexpected completed batch-advance diagnostic")
         return
      end if
      deallocate(error)
   else
      call test_failed(error, "Completed exchange image batches advanced again")
      return
   end if
   call cp2k_exchange_stream_batch_result(stream_batch, calc_cp2k, mol, &
      & cache_cp2k, cache_cp2k_plan, vsh_stream_order, energy_stream_order, error)
   if (allocated(error)) return
   call cp2k_exchange_stream_batch_result(stream_batch, calc_cp2k, mol, &
      & cache_cp2k, cache_cp2k_plan, vsh_stream, energy_stream, error)
   if (allocated(error)) then
      if (trim(error%message) /= &
         & "g-xTB exchange image-batch result was already pulled") then
         deallocate(error)
         call test_failed(error, "Unexpected duplicate batch-result diagnostic")
         return
      end if
      deallocate(error)
   else
      call test_failed(error, "Exchange image-batch result was pulled twice")
      return
   end if
   call cp2k_exchange_stream_end(stream_batch, error)
   if (allocated(error)) return
   err = max(err, maxval(abs(fock_batch-batch_sentinel-fock_order)), &
      & maxval(abs(vsh_stream_order-vsh_order)), &
      & abs(energy_stream_order-energy_order))

   ! A completed batched transaction may not be closed before its scalar result
   ! was consumed.  The failed close must clear the stream for immediate reuse.
   fock_batch = (0.0_wp, 0.0_wp)
   call cp2k_exchange_stream_begin(stream_batch, calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, &
      & wfn, nmesh, kfrac_order, weights_order, error, &
      & mode=cp2k_exchange_stream_batched, batch_size=ncell)
   if (allocated(error)) return
   do j = 1, ncell
      ik = korder(j)
      call cp2k_exchange_stream_push(stream_batch, ik, &
         & density_order(:, :, :, ik), overlap_order(:, :, ik), error)
      if (allocated(error)) return
   end do
   call cp2k_exchange_stream_batch_apply(stream_batch, calc_cp2k, mol, &
      & cache_cp2k, cache_cp2k_plan, error)
   if (allocated(error)) return
   do j = 1, ncell
      ik = korder(j)
      call cp2k_exchange_stream_batch_pull(stream_batch, ik, &
         & overlap_order(:, :, ik), fock_batch(:, :, :, ik), error)
      if (allocated(error)) return
   end do
   call cp2k_exchange_stream_batch_advance(stream_batch, batches_done, error)
   if (allocated(error)) return
   if (.not.batches_done) err = max(err, 1.0_wp)
   call cp2k_exchange_stream_end(stream_batch, error)
   if (allocated(error)) then
      if (trim(error%message) /= &
         & "g-xTB exchange stream ended with missing batch result") then
         deallocate(error)
         call test_failed(error, "Unexpected missing batch-result close diagnostic")
         return
      end if
      deallocate(error)
   else
      call test_failed(error, "Exchange image batches ended without a result pull")
      return
   end if

   ! Batch sizes equal to and greater than K both allocate exactly K image
   ! slots, report complete-image storage, and reproduce the additive oracle.
   do ibatch_case = 1, 2
      if (ibatch_case == 1) then
         requested_batch_size = ncell
      else
         requested_batch_size = ncell+2
      end if
      fock_batch = batch_sentinel
      call cp2k_exchange_stream_begin(stream_batch, calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, &
         & wfn, nmesh, kfrac_order, weights_order, error, &
         & mode=cp2k_exchange_stream_batched, &
         & batch_size=requested_batch_size)
      if (allocated(error)) return
      if (.not.cp2k_exchange_stream_has_full_mesh_storage(stream_batch)) &
         & err = max(err, 1.0_wp)
      if (allocated(cache_cp2k_plan%bvk_phase_forward) &
         & .or. allocated(cache_cp2k_plan%bvk_phase_inverse)) err = max(err, 1.0_wp)
      peak_expected = 4_int64*int(nao, int64)*int(nao, int64) &
         & *int(ncell, int64) + 3_int64*int(nao, int64) &
         & +2_int64*int(ncell, int64) &
         & +3_int64*int(maxval(nmesh), int64)
      peak_expected = peak_expected + max( &
         & 2_int64*int(nao, int64)**2, &
         & 5_int64*int(nao, int64)**2+6_int64*int(nao, int64), &
         & 5_int64*int(nao, int64)**2+int(ncell, int64))
      if (cp2k_exchange_stream_peak_complex_elements(stream_batch) &
         & /= peak_expected) err = max(err, 1.0_wp)
      call cp2k_exchange_stream_batch_bounds(stream_batch, first_image, &
         & last_image, error)
      if (allocated(error)) return
      if (first_image /= 1 .or. last_image /= ncell) err = max(err, 1.0_wp)
      do j = 1, ncell
         ik = korder(j)
         call cp2k_exchange_stream_push(stream_batch, ik, &
            & density_order(:, :, :, ik), overlap_order(:, :, ik), error)
         if (allocated(error)) return
      end do
      call cp2k_exchange_stream_batch_apply(stream_batch, calc_cp2k, mol, &
         & cache_cp2k, cache_cp2k_plan, error)
      if (allocated(error)) return
      do j = 1, ncell
         ik = korder(j)
         call cp2k_exchange_stream_batch_pull(stream_batch, ik, &
            & overlap_order(:, :, ik), fock_batch(:, :, :, ik), error)
         if (allocated(error)) return
      end do
      call cp2k_exchange_stream_batch_advance(stream_batch, batches_done, error)
      if (allocated(error)) return
      if (.not.batches_done) err = max(err, 1.0_wp)
      call cp2k_exchange_stream_batch_result(stream_batch, calc_cp2k, mol, &
         & cache_cp2k, cache_cp2k_plan, vsh_stream_order, energy_stream_order, error)
      if (allocated(error)) return
      call cp2k_exchange_stream_end(stream_batch, error)
      if (allocated(error)) return
      err = max(err, &
         & maxval(abs(fock_batch-batch_sentinel-fock_order)), &
         & maxval(abs(vsh_stream_order-vsh_order)), &
         & abs(energy_stream_order-energy_order))
   end do

   ! Independent bounded streams may own disjoint image ranges.  Each range
   ! below is completed with exactly one K-block submission, as intended for
   ! MPI image-range distribution.  The externally summed partial energy,
   ! shell potential, and additive Fock response must reproduce the unchanged
   ! full-mesh oracle without storing a complete AO image mesh on either tile.
   fock_batch = batch_sentinel
   vsh_stream = 0.0_wp
   energy_tile_sum = 0.0_wp
   bvk_builds = cache_cp2k_plan%bvk_plan_builds
   do ibatch_case = 1, 2
      if (ibatch_case == 1) then
         first_image = 1
         requested_batch_size = 1
      else
         first_image = 2
         requested_batch_size = ncell-1
      end if
      call cp2k_exchange_stream_begin(stream_batch, calc_cp2k, mol, &
         & cache_cp2k, cache_cp2k_plan, wfn, nmesh, kfrac_order, weights_order, error, &
         & mode=cp2k_exchange_stream_batched, &
         & batch_size=requested_batch_size+2, image_first=first_image, &
         & image_count=requested_batch_size)
      if (allocated(error)) return
      if (cp2k_exchange_stream_has_full_mesh_storage(stream_batch)) &
         & err = max(err, 1.0_wp)
      call cp2k_exchange_stream_batch_bounds(stream_batch, i, j, error)
      if (allocated(error)) return
      if (i /= first_image &
         & .or. j /= first_image+requested_batch_size-1) &
         & err = max(err, 1.0_wp)
      do j = 1, ncell
         ik = korder(j)
         call cp2k_exchange_stream_push(stream_batch, ik, &
            & density_order(:, :, :, ik), overlap_order(:, :, ik), error)
         if (allocated(error)) return
      end do
      call cp2k_exchange_stream_batch_apply(stream_batch, calc_cp2k, mol, &
         & cache_cp2k, cache_cp2k_plan, error)
      if (allocated(error)) return
      do j = 1, ncell
         ik = korder(j)
         call cp2k_exchange_stream_batch_pull(stream_batch, ik, &
            & overlap_order(:, :, ik), fock_batch(:, :, :, ik), error)
         if (allocated(error)) return
      end do
      call cp2k_exchange_stream_batch_advance(stream_batch, batches_done, error)
      if (allocated(error)) return
      if (.not.batches_done) err = max(err, 1.0_wp)
      call cp2k_exchange_stream_batch_result(stream_batch, calc_cp2k, mol, &
         & cache_cp2k, cache_cp2k_plan, vsh_stream_order, energy_tile, error)
      if (allocated(error)) return
      vsh_stream = vsh_stream+vsh_stream_order
      energy_tile_sum = energy_tile_sum+energy_tile
      call cp2k_exchange_stream_end(stream_batch, error)
      if (allocated(error)) return
   end do
   if (cache_cp2k_plan%bvk_plan_builds /= bvk_builds) err = max(err, 1.0_wp)
   err = max(err, maxval(abs(fock_batch-batch_sentinel-fock_order)), &
      & maxval(abs(vsh_stream-vsh_order)), abs(energy_tile_sum-energy_order))

   ! The cached reverse path must reproduce the independently evaluated
   ! overlap response, atomic forces, and homogeneous-strain derivative.
   allocate(overlap_grad_cp2k(nao, nao, ncell), &
      & gradient_cp2k(3, mol%nat), sigma_cp2k(3, 3))
   call cp2k_exchange_kmesh_gradient(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, nmesh, &
      & kfrac, weights, density_mesh, overlap_mesh, overlap_grad_cp2k, &
      & gradient_cp2k, sigma_cp2k, error)
   if (allocated(error)) return
   err = max(err, maxval(abs(overlap_grad_cp2k-overlap_grad_mesh)), &
      & maxval(abs(gradient_cp2k-gradient_mesh)), &
      & maxval(abs(sigma_cp2k-sigma_mesh)))

   ! Exact composition test: both compact transform backends must retain the
   ! bounded image protocol and reproduce the dense forward and reverse maps.
   ! The 3x1x1 mesh also exercises the mixed-radix prime-axis leaf.
   call check_compact_stream_composition(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn, &
      & nmesh, kfrac, weights, density_mesh, overlap_mesh, fock_cp2k, &
      & vsh_cp2k, energy_mesh, overlap_grad_cp2k, gradient_cp2k, sigma_cp2k, &
      & cp2k_exchange_transform_separable_dft, stream_compact_separable_err, &
      & stream_compact_separable_peak, error)
   if (allocated(error)) return
   call check_compact_stream_composition(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn, &
      & nmesh, kfrac, weights, density_mesh, overlap_mesh, fock_cp2k, &
      & vsh_cp2k, energy_mesh, overlap_grad_cp2k, gradient_cp2k, sigma_cp2k, &
      & cp2k_exchange_transform_mixed_radix_fft, stream_compact_fft_err, &
      & stream_compact_fft_peak, error)
   if (allocated(error)) return
   err = max(err, stream_compact_separable_err, stream_compact_fft_err)
   separable_err = max(separable_err, stream_compact_separable_err)
   fft_err = max(fft_err, stream_compact_fft_err)
   print *, "STREAMED_SEPARABLE_COMPOSITION_MAX_RESID", &
      & stream_compact_separable_err, stream_compact_separable_peak
   print *, "STREAMED_MIXED_RADIX_COMPOSITION_MAX_RESID", &
      & stream_compact_fft_err, stream_compact_fft_peak

   ! The production tile policy must preserve all dense forward and reverse
   ! observables, not just the tile-size arithmetic or total energy.
   call check_compact_stream_composition(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn, &
      & nmesh, kfrac, weights, density_mesh, overlap_mesh, fock_cp2k, &
      & vsh_cp2k, energy_mesh, overlap_grad_cp2k, gradient_cp2k, sigma_cp2k, &
      & cp2k_exchange_transform_separable_dft, stream_compact_separable_err, &
      & stream_compact_separable_peak, error, adaptive_tiles=.true.)
   if (allocated(error)) return
   call check_compact_stream_composition(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn, &
      & nmesh, kfrac, weights, density_mesh, overlap_mesh, fock_cp2k, &
      & vsh_cp2k, energy_mesh, overlap_grad_cp2k, gradient_cp2k, sigma_cp2k, &
      & cp2k_exchange_transform_mixed_radix_fft, stream_compact_fft_err, &
      & stream_compact_fft_peak, error, adaptive_tiles=.true.)
   if (allocated(error)) return
   err = max(err, stream_compact_separable_err, stream_compact_fft_err)
   separable_err = max(separable_err, stream_compact_separable_err)
   fft_err = max(fft_err, stream_compact_fft_err)
   print *, "ADAPTIVE_SEPARABLE_COMPOSITION_MAX_RESID", stream_compact_separable_err
   print *, "ADAPTIVE_MIXED_RADIX_COMPOSITION_MAX_RESID", stream_compact_fft_err

   ! Compose the image and AO-pair policies together, including all reverse
   ! responses. Keep the single-image cases above as independent controls.
   call check_compact_stream_composition(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn, &
      & nmesh, kfrac, weights, density_mesh, overlap_mesh, fock_cp2k, &
      & vsh_cp2k, energy_mesh, overlap_grad_cp2k, gradient_cp2k, sigma_cp2k, &
      & cp2k_exchange_transform_separable_dft, stream_compact_separable_err, &
      & stream_compact_separable_peak, error, adaptive_tiles=.true., adaptive_images=.true.)
   if (allocated(error)) return
   call check_compact_stream_composition(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn, &
      & nmesh, kfrac, weights, density_mesh, overlap_mesh, fock_cp2k, &
      & vsh_cp2k, energy_mesh, overlap_grad_cp2k, gradient_cp2k, sigma_cp2k, &
      & cp2k_exchange_transform_mixed_radix_fft, stream_compact_fft_err, &
      & stream_compact_fft_peak, error, adaptive_tiles=.true., adaptive_images=.true.)
   if (allocated(error)) return
   err = max(err, stream_compact_separable_err, stream_compact_fft_err)
   separable_err = max(separable_err, stream_compact_separable_err)
   fft_err = max(fft_err, stream_compact_fft_err)
   print *, "ADAPTIVE_IMAGES_SEPARABLE_COMPOSITION_MAX_RESID", stream_compact_separable_err
   print *, "ADAPTIVE_IMAGES_MIXED_RADIX_COMPOSITION_MAX_RESID", stream_compact_fft_err

   ! The reverse-only reduced and bounded transactions must reproduce the
   ! complete-mesh oracle without retaining density, overlap, Fock, or overlap
   ! adjoint meshes.  In particular, no forward Fock pull is part of this
   ! derivative-only protocol.
   allocate(overlap_grad_stream_reduced(nao, nao, ncell), &
      & overlap_grad_stream_batch(nao, nao, ncell), &
      & overlap_grad_stream_probe(nao, nao), source=(0.0_wp, 0.0_wp))
   allocate(gradient_stream_reverse(3, mol%nat), &
      & sigma_stream_reverse(3, 3), source=0.0_wp)
   call cp2k_exchange_stream_begin(stream_reverse_reduced, calc_cp2k, mol, &
      & cache_cp2k, cache_cp2k_plan, wfn, nmesh, kfrac, weights, error, &
      & mode=cp2k_exchange_stream_reduced)
   if (allocated(error)) return
   do ik = 1, ncell
      call cp2k_exchange_stream_push(stream_reverse_reduced, ik, &
         & density_mesh(:, :, :, ik), overlap_mesh(:, :, ik), error)
      if (allocated(error)) return
   end do
   call cp2k_exchange_stream_apply(stream_reverse_reduced, calc_cp2k, mol, &
      & cache_cp2k, cache_cp2k_plan, vsh_stream, energy_stream, error)
   if (allocated(error)) return
   call cp2k_exchange_stream_reverse_begin(stream_reverse_reduced, mol, error)
   if (allocated(error)) return
   do ik = 1, ncell
      call cp2k_exchange_stream_reverse_push(stream_reverse_reduced, ik, &
         & density_mesh(:, :, :, ik), overlap_mesh(:, :, ik), error)
      if (allocated(error)) return
   end do
   call cp2k_exchange_stream_reverse_batch_apply(stream_reverse_reduced, &
      & calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, error)
   if (allocated(error)) return
   do ik = 1, ncell
      call cp2k_exchange_stream_reverse_batch_pull(stream_reverse_reduced, &
         & ik, density_mesh(:, :, :, ik), overlap_mesh(:, :, ik), &
         & overlap_grad_stream_reduced(:, :, ik), error)
      if (allocated(error)) return
   end do
   call cp2k_exchange_stream_reverse_batch_advance(stream_reverse_reduced, &
      & batches_done, error)
   if (allocated(error)) return
   if (.not.batches_done) err = max(err, 1.0_wp)
   call cp2k_exchange_stream_reverse_result(stream_reverse_reduced, calc_cp2k, &
      & mol, cache_cp2k, cache_cp2k_plan, gradient_stream_reverse, sigma_stream_reverse, error)
   if (allocated(error)) return
   err = max(err, &
      & maxval(abs(overlap_grad_stream_reduced-overlap_grad_cp2k)), &
      & maxval(abs(gradient_stream_reverse-gradient_cp2k)), &
      & maxval(abs(sigma_stream_reverse-sigma_cp2k)))
   stream_reverse_err = max(stream_reverse_err, &
      & maxval(abs(overlap_grad_stream_reduced-overlap_grad_cp2k)), &
      & maxval(abs(gradient_stream_reverse-gradient_cp2k)), &
      & maxval(abs(sigma_stream_reverse-sigma_cp2k)))
   call cp2k_exchange_stream_end(stream_reverse_reduced, error)
   if (allocated(error)) return

   call cp2k_exchange_stream_begin(stream_reverse_batch, calc_cp2k, mol, &
      & cache_cp2k, cache_cp2k_plan, wfn, nmesh, kfrac, weights, error, &
      & mode=cp2k_exchange_stream_batched, batch_size=2, &
      & reverse_only=.true.)
   if (allocated(error)) return
   if (cp2k_exchange_stream_has_full_mesh_storage(stream_reverse_batch)) &
      & err = max(err, 1.0_wp)
   peak_expected = 6_int64*int(nao, int64)**2*2_int64 &
      & *int(wfn%nspin, int64)+5_int64*int(nao, int64) &
      & *int(wfn%nspin, int64)+2_int64*int(ncell, int64) &
      & +3_int64*int(maxval(nmesh), int64)+10_int64*int(nao, int64)**2
   if (cp2k_exchange_stream_peak_complex_elements(stream_reverse_batch) &
      & /= peak_expected) err = max(err, 1.0_wp)
   batches_done = .false.
   do while (.not.batches_done)
      do ik = 1, ncell
         call cp2k_exchange_stream_reverse_push(stream_reverse_batch, ik, &
            & density_mesh(:, :, :, ik), overlap_mesh(:, :, ik), error)
         if (allocated(error)) return
         if (ik == 1) then
            call cp2k_exchange_stream_reverse_push(stream_reverse_batch, ik, &
               & density_mesh(:, :, :, ik), overlap_mesh(:, :, ik), error)
            if (allocated(error)) then
               if (trim(error%message) /= &
                  & "duplicate g-xTB exchange reverse block") then
                  deallocate(error)
                  call test_failed(error, &
                     & "Unexpected duplicate reverse-push diagnostic")
                  return
               end if
               deallocate(error)
            else
               call test_failed(error, "Duplicate reverse push was accepted")
               return
            end if
         end if
      end do
      call cp2k_exchange_stream_reverse_batch_apply(stream_reverse_batch, &
         & calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, error)
      if (allocated(error)) return
      call cp2k_exchange_stream_reverse_batch_apply(stream_reverse_batch, &
         & calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, error)
      if (allocated(error)) then
         if (trim(error%message) /= &
            & "g-xTB exchange stream reverse batch was already applied") then
            deallocate(error)
            call test_failed(error, &
               & "Unexpected duplicate reverse-apply diagnostic")
            return
         end if
         deallocate(error)
      else
         call test_failed(error, "Duplicate reverse-batch apply was accepted")
         return
      end if
      do ik = 1, ncell
         call cp2k_exchange_stream_reverse_batch_pull(stream_reverse_batch, &
            & ik, density_mesh(:, :, :, ik), overlap_mesh(:, :, ik), &
            & overlap_grad_stream_batch(:, :, ik), error)
         if (allocated(error)) return
         if (ik == 1) then
            call cp2k_exchange_stream_reverse_batch_pull( &
               & stream_reverse_batch, ik, density_mesh(:, :, :, ik), &
               & overlap_mesh(:, :, ik), &
               & overlap_grad_stream_probe, error)
            if (allocated(error)) then
               if (trim(error%message) /= &
                  & "duplicate g-xTB exchange reverse-batch pull") then
                  deallocate(error)
                  call test_failed(error, &
                     & "Unexpected duplicate reverse-pull diagnostic")
                  return
               end if
               deallocate(error)
            else
               call test_failed(error, &
                  & "Duplicate reverse-batch pull was accepted")
               return
            end if
         end if
      end do
      call cp2k_exchange_stream_reverse_batch_advance(stream_reverse_batch, &
         & batches_done, error)
      if (allocated(error)) return
   end do
   call cp2k_exchange_stream_reverse_result(stream_reverse_batch, calc_cp2k, &
      & mol, cache_cp2k, cache_cp2k_plan, gradient_stream_reverse, sigma_stream_reverse, error)
   if (allocated(error)) return
   call cp2k_exchange_stream_reverse_result(stream_reverse_batch, calc_cp2k, &
      & mol, cache_cp2k, cache_cp2k_plan, gradient_stream_probe, sigma_stream_probe, error)
   if (allocated(error)) then
      if (trim(error%message) /= &
         & "g-xTB exchange reverse result was already pulled") then
         deallocate(error)
         call test_failed(error, &
            & "Unexpected duplicate reverse-result diagnostic")
         return
      end if
      deallocate(error)
   else
      call test_failed(error, "Duplicate reverse result was accepted")
      return
   end if
   err = max(err, &
      & maxval(abs(overlap_grad_stream_batch-overlap_grad_cp2k)), &
      & maxval(abs(gradient_stream_reverse-gradient_cp2k)), &
      & maxval(abs(sigma_stream_reverse-sigma_cp2k)))
   stream_reverse_err = max(stream_reverse_err, &
      & maxval(abs(overlap_grad_stream_batch-overlap_grad_cp2k)), &
      & maxval(abs(gradient_stream_reverse-gradient_cp2k)), &
      & maxval(abs(sigma_stream_reverse-sigma_cp2k)))
   print *, "STREAM_REVERSE_RKS_3X1X1_RESID", &
      & maxval(abs(overlap_grad_stream_batch-overlap_grad_cp2k)), &
      & maxval(abs(gradient_stream_reverse-gradient_cp2k)), &
      & maxval(abs(sigma_stream_reverse-sigma_cp2k)), &
      & cp2k_exchange_stream_peak_complex_elements(stream_reverse_batch)
   call cp2k_exchange_stream_end(stream_reverse_batch, error)
   if (allocated(error)) return

   ! Exercise the separate reverse-stream phase without changing the forward
   ! transaction tests above.  A stale onsite state must be rejected before
   ! the unchanged whole-mesh gradient oracle is entered.
   allocate(overlap_grad_stream_reverse(nao, nao, ncell))
   call cp2k_exchange_stream_begin(stream_reverse, calc_cp2k, mol, &
      & cache_cp2k, cache_cp2k_plan, wfn, nmesh, kfrac, weights, error, &
      & mode=cp2k_exchange_stream_oracle)
   if (allocated(error)) return
   if (.not.cp2k_exchange_stream_has_full_mesh_storage(stream_reverse)) &
      & err = max(err, 1.0_wp)
   do ik = 1, ncell
      call cp2k_exchange_stream_push(stream_reverse, ik, &
         & density_mesh(:, :, :, ik), overlap_mesh(:, :, ik), error)
      if (allocated(error)) return
   end do
   call cp2k_exchange_stream_apply(stream_reverse, calc_cp2k, mol, &
      & cache_cp2k, cache_cp2k_plan, vsh_stream, energy_stream, error)
   if (allocated(error)) return
   do ik = 1, ncell
      call cp2k_exchange_stream_get(stream_reverse, ik, &
         & fock_stream_order(:, :, :, ik), error)
      if (allocated(error)) return
   end do
   err = max(err, maxval(abs(fock_stream_order-fock_cp2k)), &
      & maxval(abs(vsh_stream-vsh_cp2k)), abs(energy_stream-energy_mesh))

   onsite_probe = ptr_cp2k%g_onsfx(1, 1, 1)
   ptr_cp2k%g_onsfx(1, 1, 1) = onsite_probe + 1.0e-8_wp
   call cp2k_exchange_stream_reverse_apply(stream_reverse, calc_cp2k, mol, &
      & cache_cp2k, cache_cp2k_plan, gradient_stream_probe, sigma_stream_probe, error)
   if (allocated(error)) then
      if (trim(error%message) /= &
         & "g-xTB exchange stream state was invalidated before reverse apply") then
         deallocate(error)
         call test_failed(error, "Unexpected stale reverse-stream diagnostic")
         return
      end if
      deallocate(error)
   else
      call test_failed(error, "Stale exchange reverse stream was applied")
      return
   end if
   ptr_cp2k%g_onsfx(1, 1, 1) = onsite_probe

   call cp2k_exchange_stream_reverse_apply(stream_reverse, calc_cp2k, mol, &
      & cache_cp2k, cache_cp2k_plan, gradient_stream_reverse, sigma_stream_reverse, error)
   if (allocated(error)) return
   call cp2k_exchange_stream_reverse_apply(stream_reverse, calc_cp2k, mol, &
      & cache_cp2k, cache_cp2k_plan, gradient_stream_probe, sigma_stream_probe, error)
   if (allocated(error)) then
      if (trim(error%message) /= &
         & "g-xTB exchange stream reverse was already applied") then
         deallocate(error)
         call test_failed(error, "Unexpected repeated reverse-apply diagnostic")
         return
      end if
      deallocate(error)
   else
      call test_failed(error, "Exchange reverse stream was applied twice")
      return
   end if

   call cp2k_exchange_stream_reverse_get(stream_reverse, 1, &
      & overlap_grad_stream_reverse(:, :, 1), error)
   if (allocated(error)) return
   call cp2k_exchange_stream_reverse_get(stream_reverse, 1, &
      & overlap_grad_stream_probe, error)
   if (allocated(error)) then
      if (trim(error%message) /= &
         & "duplicate g-xTB exchange stream reverse pull") then
         deallocate(error)
         call test_failed(error, "Unexpected duplicate reverse-pull diagnostic")
         return
      end if
      deallocate(error)
   else
      call test_failed(error, "Exchange overlap adjoint was pulled twice")
      return
   end if
   do ik = 2, ncell
      call cp2k_exchange_stream_reverse_get(stream_reverse, ik, &
         & overlap_grad_stream_reverse(:, :, ik), error)
      if (allocated(error)) return
   end do
   call cp2k_exchange_stream_end(stream_reverse, error)
   if (allocated(error)) return
   err = max(err, &
      & maxval(abs(overlap_grad_stream_reverse-overlap_grad_cp2k)), &
      & maxval(abs(overlap_grad_stream_reverse-overlap_grad_mesh)), &
      & maxval(abs(gradient_stream_reverse-gradient_cp2k)), &
      & maxval(abs(gradient_stream_reverse-gradient_mesh)), &
      & maxval(abs(sigma_stream_reverse-sigma_cp2k)), &
      & maxval(abs(sigma_stream_reverse-sigma_mesh)))

   ! Once reverse evaluation was requested, closing the stream is legal only
   ! after every overlap-adjoint block was pulled exactly once.
   call cp2k_exchange_stream_begin(stream_reverse_missing, calc_cp2k, mol, &
      & cache_cp2k, cache_cp2k_plan, wfn, nmesh, kfrac, weights, error, &
      & mode=cp2k_exchange_stream_oracle)
   if (allocated(error)) return
   do j = 1, ncell
      ik = korder(j)
      call cp2k_exchange_stream_push(stream_reverse_missing, ik, &
         & density_mesh(:, :, :, ik), overlap_mesh(:, :, ik), error)
      if (allocated(error)) return
   end do
   call cp2k_exchange_stream_apply(stream_reverse_missing, calc_cp2k, mol, &
      & cache_cp2k, cache_cp2k_plan, vsh_stream, energy_stream, error)
   if (allocated(error)) return
   do ik = 1, ncell
      call cp2k_exchange_stream_get(stream_reverse_missing, ik, &
         & fock_stream_order(:, :, :, ik), error)
      if (allocated(error)) return
   end do
   call cp2k_exchange_stream_reverse_apply(stream_reverse_missing, calc_cp2k, &
      & mol, cache_cp2k, cache_cp2k_plan, gradient_stream_probe, sigma_stream_probe, error)
   if (allocated(error)) return
   call cp2k_exchange_stream_reverse_get(stream_reverse_missing, 1, &
      & overlap_grad_stream_probe, error)
   if (allocated(error)) return
   call cp2k_exchange_stream_end(stream_reverse_missing, error)
   if (allocated(error)) then
      if (trim(error%message) /= &
         & "g-xTB exchange stream ended with missing reverse pulls") then
         deallocate(error)
         call test_failed(error, "Unexpected missing reverse-pull diagnostic")
         return
      end if
      deallocate(error)
   else
      call test_failed(error, "Exchange stream with missing reverse pulls ended")
      return
   end if

   bvk_builds = cache_cp2k_plan%bvk_plan_builds
   call cp2k_exchange_kmesh(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn, nmesh, &
      & kfrac, weights, density_mesh, overlap_mesh, fock_cp2k, &
      & vsh_cp2k, energy_mesh, error)
   if (allocated(error)) return
   if (cache_cp2k_plan%bvk_plan_builds /= bvk_builds) err = max(err, 1.0_wp)

   ! Every stored model component and the active representative ordering are
   ! part of the persistent-plan identity.  Mutating any one of them must
   ! reject the cache before a stale transform can be selected.
   select type (fock_cp2k_model => calc_cp2k%exchange)
   type is (exchange_fock)
      call check_bvk_model_signature(fock_cp2k_model, cache_cp2k_plan, error)
      if (allocated(error)) return
      kernel_order_probe = cache_cp2k_plan%bvk_kernel
      if (.not.cache_cp2k_plan%plan_matches(fock_cp2k_model, mol, &
         & kernel_order_probe, kfrac, weights, require_phases=.false.)) then
         call test_failed(error, "Fresh BvK representative order does not match")
         return
      end if
      kernel_order_probe%reps(:, [1, 2]) = &
         & kernel_order_probe%reps(:, [2, 1])
      if (cache_cp2k_plan%plan_matches(fock_cp2k_model, mol, &
         & kernel_order_probe, kfrac, weights, require_phases=.false.)) then
         call test_failed(error, "BvK plan accepted a changed representative order")
         return
      end if
   class default
      call test_failed(error, "Unexpected exchange implementation")
      return
   end select

   call cp2k_exchange_kmesh(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn, nmesh, &
      & kfrac_twist, weights, density_mesh, overlap_mesh, fock_cp2k, &
      & vsh_cp2k, energy_mesh, error)
   if (allocated(error)) return
   if (cache_cp2k_plan%bvk_plan_builds /= bvk_builds+1) err = max(err, 1.0_wp)
   err = max(err, maxval(abs(fock_cp2k-fock_twist)), &
      & maxval(abs(vsh_cp2k-vsh_twist)), abs(energy_mesh-energy_twist))

   bvk_builds = cache_cp2k_plan%bvk_plan_builds
   call cp2k_exchange_kmesh(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn, nmesh, &
      & kfrac_order, weights_order, density_order, overlap_order, fock_cp2k, &
      & vsh_cp2k, energy_mesh, error)
   if (allocated(error)) return
   if (cache_cp2k_plan%bvk_plan_builds /= bvk_builds+1) err = max(err, 1.0_wp)
   if (any(cache_cp2k_plan%bvk_input_to_grid /= korder)) err = max(err, 1.0_wp)
   if (any(cache_cp2k_plan%bvk_grid_to_input /= [2, 3, 1])) &
      & err = max(err, 1.0_wp)
   err = max(err, maxval(abs(fock_cp2k-fock_order)), &
      & maxval(abs(vsh_cp2k-vsh_order)), abs(energy_mesh-energy_order))
   call cp2k_exchange_kmesh(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn, nmesh, &
      & kfrac_order, weights_order, density_order, overlap_order, fock_cp2k, &
      & vsh_cp2k, energy_mesh, error, &
      & transform_backend=cp2k_exchange_transform_separable_dft)
   if (allocated(error)) return
   separable_err = max(separable_err, &
      & maxval(abs(fock_cp2k-fock_order)), &
      & maxval(abs(vsh_cp2k-vsh_order)), abs(energy_mesh-energy_order))
   call cp2k_exchange_kmesh(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn, nmesh, &
      & kfrac_order, weights_order, density_order, overlap_order, fock_cp2k, &
      & vsh_cp2k, energy_mesh, error, &
      & transform_backend=cp2k_exchange_transform_mixed_radix_fft)
   if (allocated(error)) return
   fft_err = max(fft_err, maxval(abs(fock_cp2k-fock_order)), &
      & maxval(abs(vsh_cp2k-vsh_order)), abs(energy_mesh-energy_order))
   print *, "MIXED_RADIX_FFT_RKS_FORWARD_RESID", &
      & maxval(abs(fock_cp2k-fock_order)), maxval(abs(vsh_cp2k-vsh_order)), &
      & abs(energy_mesh-energy_order)
   call check_exchange_transform_selectors(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn, &
      & nmesh, kfrac_order, weights_order, density_order, overlap_order, &
      & error, selector_separable_err, selector_fft_err, rks_fft_forward, &
      & rks_fft_reverse)
   if (allocated(error)) return
   separable_err = max(separable_err, selector_separable_err)
   fft_err = max(fft_err, selector_fft_err)
   print *, "TRANSFORM_SELECTOR_RKS_SEPARABLE_MAX_RESID", &
      & selector_separable_err
   print *, "MIXED_RADIX_FFT_RKS_REVERSE_RESID", rks_fft_reverse
   err = max(err, separable_err)

   ! Every geometry, cell, mesh, and model dependency must trigger an actual
   ! rebuild.  The AO inputs can remain fixed here because this test probes
   ! plan ownership and invalidation rather than a self-consistent state.
   bvk_builds = cache_cp2k_plan%bvk_plan_builds
   mol_probe = mol
   mol_probe%xyz(1, 1) = mol_probe%xyz(1, 1) + 1.0e-4_wp
   call cp2k_exchange_kmesh(calc_cp2k, mol_probe, cache_cp2k, cache_cp2k_plan, wfn, nmesh, &
      & kfrac_order, weights_order, density_order, overlap_order, fock_cp2k, &
      & vsh_cp2k, energy_mesh, error)
   if (allocated(error)) return
   if (cache_cp2k_plan%bvk_plan_builds /= bvk_builds+1) err = max(err, 1.0_wp)
   bvk_builds = cache_cp2k_plan%bvk_plan_builds
   mol_probe = mol
   mol_probe%lattice(1, 1) = mol_probe%lattice(1, 1) + 1.0e-4_wp
   call cp2k_exchange_kmesh(calc_cp2k, mol_probe, cache_cp2k, cache_cp2k_plan, wfn, nmesh, &
      & kfrac_order, weights_order, density_order, overlap_order, fock_cp2k, &
      & vsh_cp2k, energy_mesh, error)
   if (allocated(error)) return
   if (cache_cp2k_plan%bvk_plan_builds /= bvk_builds+1) err = max(err, 1.0_wp)

   bvk_builds = cache_cp2k_plan%bvk_plan_builds
   kfrac_bad = 0.0_wp
   do ik = 1, ncell
      kfrac_bad(2, ik) = real(ik-1, wp)/real(ncell, wp)
   end do
   call cp2k_exchange_kmesh(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn, [1, 3, 1], &
      & kfrac_bad, weights, density_mesh, overlap_mesh, fock_cp2k, &
      & vsh_cp2k, energy_mesh, error)
   if (allocated(error)) return
   if (cache_cp2k_plan%bvk_plan_builds /= bvk_builds+1 &
      & .or. any(cache_cp2k_plan%bvk_kernel%nmesh /= [1, 3, 1])) &
      & err = max(err, 1.0_wp)

   bvk_builds = cache_cp2k_plan%bvk_plan_builds
   select type (fock_cp2k_model => calc_cp2k%exchange)
   type is (exchange_fock)
      fock_cp2k_model%corr_exp = fock_cp2k_model%corr_exp + 1.0e-6_wp
   end select
   call cp2k_exchange_kmesh(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn, nmesh, &
      & kfrac_order, weights_order, density_order, overlap_order, fock_cp2k, &
      & vsh_cp2k, energy_mesh, error)
   if (allocated(error)) return
   if (cache_cp2k_plan%bvk_plan_builds /= bvk_builds+1) err = max(err, 1.0_wp)
   select type (fock_cp2k_model => calc_cp2k%exchange)
   type is (exchange_fock)
      fock_cp2k_model%corr_exp = fock_cp2k_model%corr_exp - 1.0e-6_wp
   end select

   ! Reject incomplete and duplicate grids.  A boundary rejection leaves a
   ! valid plan untouched; failed internal construction invalidates it, and
   ! the next valid call must recover without incrementing on the failures.
   bvk_builds = cache_cp2k_plan%bvk_plan_builds
   call cp2k_exchange_kmesh(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn, [-1, 1, 1], &
      & kfrac_order, weights_order, density_order, overlap_order, fock_cp2k, &
      & vsh_cp2k, energy_mesh, error, &
      & transform_backend=cp2k_exchange_transform_mixed_radix_fft)
   if (allocated(error)) then
      if (trim(error%message) /= &
         & "g-xTB exchange requires a complete regular k mesh") then
         deallocate(error)
         call test_failed(error, "Unexpected negative FFT-mesh diagnostic")
         return
      end if
      deallocate(error)
   else
      call test_failed(error, "Negative mixed-radix FFT mesh was accepted")
      return
   end if
   call cp2k_exchange_kmesh(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn, nmesh_overflow, &
      & kfrac_order, weights_order, density_order, overlap_order, fock_cp2k, &
      & vsh_cp2k, energy_mesh, error, &
      & transform_backend=cp2k_exchange_transform_mixed_radix_fft)
   if (allocated(error)) then
      if (trim(error%message) /= &
         & "g-xTB exchange requires a complete regular k mesh") then
         deallocate(error)
         call test_failed(error, "Unexpected overflowing FFT-mesh diagnostic")
         return
      end if
      deallocate(error)
   else
      call test_failed(error, "Overflowing mixed-radix FFT mesh was accepted")
      return
   end if
   if (cache_cp2k_plan%bvk_plan_builds /= bvk_builds &
      & .or. .not.cache_cp2k_plan%bvk_plan_valid) err = max(err, 1.0_wp)

   call cp2k_exchange_kmesh(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn, [4, 1, 1], &
      & kfrac_order, weights_order, density_order, overlap_order, fock_cp2k, &
      & vsh_cp2k, energy_mesh, error)
   if (allocated(error)) then
      if (trim(error%message) /= &
         & "g-xTB exchange requires a complete regular k mesh") then
         deallocate(error)
         call test_failed(error, "Unexpected incomplete-grid diagnostic")
         return
      end if
      deallocate(error)
   else
      call test_failed(error, "Incomplete whole-mesh exchange grid was accepted")
      return
   end if
   if (cache_cp2k_plan%bvk_plan_builds /= bvk_builds &
      & .or. .not.cache_cp2k_plan%bvk_plan_valid) err = max(err, 1.0_wp)

   kfrac_bad = kfrac_order
   kfrac_bad(:, 2) = kfrac_bad(:, 1)
   call cp2k_exchange_kmesh(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn, nmesh, &
      & kfrac_bad, weights_order, density_order, overlap_order, fock_cp2k, &
      & vsh_cp2k, energy_mesh, error)
   if (allocated(error)) then
      if (trim(error%message) /= &
         & "g-xTB k-point mesh contains a duplicate grid point") then
         deallocate(error)
         call test_failed(error, "Unexpected duplicate-grid diagnostic")
         return
      end if
      deallocate(error)
   else
      call test_failed(error, "Duplicate whole-mesh exchange grid point was accepted")
      return
   end if
   if (cache_cp2k_plan%bvk_plan_builds /= bvk_builds &
      & .or. cache_cp2k_plan%bvk_plan_valid) err = max(err, 1.0_wp)
   call cp2k_exchange_kmesh(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn, nmesh, &
      & kfrac_order, weights_order, density_order, overlap_order, fock_cp2k, &
      & vsh_cp2k, energy_mesh, error)
   if (allocated(error)) return
   if (cache_cp2k_plan%bvk_plan_builds /= bvk_builds+1 &
      & .or. .not.cache_cp2k_plan%bvk_plan_valid) err = max(err, 1.0_wp)

   ! A geometrically near-regular mesh is accepted by the dense diagnostic
   ! tolerance but must not be silently snapped onto the factorized grid.
   kfrac_bad = kfrac_order
   kfrac_bad(1, 1) = kfrac_bad(1, 1)+1.0e-12_wp
   call cp2k_exchange_kmesh(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn, nmesh, &
      & kfrac_bad, weights_order, density_order, overlap_order, fock_cp2k, &
      & vsh_cp2k, energy_mesh, error, &
      & transform_backend=cp2k_exchange_transform_separable_dft)
   if (allocated(error)) then
      if (trim(error%message) /= &
         & "separable-DFT g-xTB exchange requires a machine-regular grid") then
         deallocate(error)
         call test_failed(error, &
            & "Unexpected non-machine-regular-grid diagnostic")
         return
      end if
      deallocate(error)
   else
      call test_failed(error, &
         & "Non-machine-regular separable-DFT grid was accepted")
      return
   end if
   call cp2k_exchange_kmesh(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn, nmesh, &
      & kfrac_bad, weights_order, density_order, overlap_order, fock_cp2k, &
      & vsh_cp2k, energy_mesh, error, &
      & transform_backend=cp2k_exchange_transform_mixed_radix_fft)
   if (allocated(error)) then
      if (trim(error%message) /= &
         & "mixed-radix FFT g-xTB exchange requires a machine-regular grid") then
         deallocate(error)
         call test_failed(error, &
            & "Unexpected non-machine-regular FFT-grid diagnostic")
         return
      end if
      deallocate(error)
   else
      call test_failed(error, &
         & "Non-machine-regular mixed-radix FFT grid was accepted")
      return
   end if

   ! Exercise the production O(Nk) grid proof beyond the small-mesh dense
   ! orthogonality oracle.  The cached 9x9x1 forward and reverse paths are
   ! compared with the unchanged dense-transform implementation.
   allocate(kfrac_large(3, product(nmesh_large)), &
      & weights_large(product(nmesh_large)))
   twist_large = [0.037_wp, 0.021_wp, 0.013_wp]
   weights_large = 1.0_wp/real(product(nmesh_large), wp)
   do ik = 1, product(nmesh_large)
      ! A coprime affine permutation makes the external order unrelated to
      ! the canonical 9x9x1 grid order while retaining one common twist.
      j = modulo(17*(ik-1)+5, product(nmesh_large))+1
      kfrac_large(1, ik) = twist_large(1) &
         & +real(modulo(j-1, nmesh_large(1)), wp) &
         & /real(nmesh_large(1), wp)
      kfrac_large(2, ik) = twist_large(2) &
         & +real((j-1)/nmesh_large(1), wp) &
         & /real(nmesh_large(2), wp)
      kfrac_large(3, ik) = twist_large(3)
   end do
   allocate(density_large(nao, nao, 2, product(nmesh_large)), &
      & overlap_large(nao, nao, product(nmesh_large)), &
      & fock_large(nao, nao, 2, product(nmesh_large)), &
      & fock_large_oracle(nao, nao, 2, product(nmesh_large)), &
      & vsh_large(nsh), vsh_large_oracle(nsh))
   do ik = 1, product(nmesh_large)
      theta_x = 2.0_wp*pi*kfrac_large(1, ik)
      theta_y = 2.0_wp*pi*kfrac_large(2, ik)
      overlap_large(:, :, ik) = (0.0_wp, 0.0_wp)
      density_large(:, :, :, ik) = (0.0_wp, 0.0_wp)
      do ii = 1, nao
         overlap_large(ii, ii, ik) = cmplx( &
            & real(overlap_mesh(ii, ii, 1), wp) &
            & +0.003_wp*cos(theta_x+real(ii, wp)), 0.0_wp, wp)
         density_large(ii, ii, 1, ik) = cmplx( &
            & 0.56_wp*real(density_mesh(ii, ii, 1, 1), wp) &
            & +0.002_wp*cos(theta_x+real(2*ii, wp)), 0.0_wp, wp)
         density_large(ii, ii, 2, ik) = cmplx( &
            & 0.44_wp*real(density_mesh(ii, ii, 1, 1), wp) &
            & +0.0015_wp*sin(theta_y+real(3*ii, wp)), 0.0_wp, wp)
         do jj = 1, ii-1
            overlap_large(jj, ii, ik) = 0.5_wp*( &
               & overlap_mesh(jj, ii, 1)+conjg(overlap_mesh(ii, jj, 1))) &
               & +cmplx(0.001_wp*cos(theta_x+real(ii+jj, wp)), &
               & 0.001_wp*sin(theta_y+real(2*ii+jj, wp)), wp)
            overlap_large(ii, jj, ik) = &
               & conjg(overlap_large(jj, ii, ik))
            density_large(jj, ii, 1, ik) = 0.56_wp*0.5_wp*( &
               & density_mesh(jj, ii, 1, 1) &
               & +conjg(density_mesh(ii, jj, 1, 1))) &
               & +cmplx(0.0012_wp*cos(theta_x+real(ii+2*jj, wp)), &
               & 0.0009_wp*sin(theta_y+real(3*ii+jj, wp)), wp)
            density_large(ii, jj, 1, ik) = &
               & conjg(density_large(jj, ii, 1, ik))
            density_large(jj, ii, 2, ik) = 0.44_wp*0.5_wp*( &
               & density_mesh(jj, ii, 1, 1) &
               & +conjg(density_mesh(ii, jj, 1, 1))) &
               & +cmplx(0.0008_wp*sin(theta_x+real(ii+jj, wp)), &
               & 0.0011_wp*cos(theta_y+real(ii+3*jj, wp)), wp)
            density_large(ii, jj, 2, ik) = &
               & conjg(density_large(jj, ii, 2, ik))
         end do
      end do
   end do
   if (maxval(abs(overlap_large(:, :, 1) &
      & -conjg(transpose(overlap_large(:, :, 1))))) > 1.0e-14_wp &
      & .or. maxval(abs(density_large(:, :, 1, 1) &
      & -conjg(transpose(density_large(:, :, 1, 1))))) > 1.0e-14_wp &
      & .or. maxval(abs(density_large(:, :, 2, 1) &
      & -conjg(transpose(density_large(:, :, 2, 1))))) > 1.0e-14_wp) &
      & err = max(err, 1.0_wp)
   if (maxval(abs(overlap_large(:, :, 2)-overlap_large(:, :, 1))) &
      & < 1.0e-6_wp .or. maxval(abs(density_large(:, :, :, 2) &
      & -density_large(:, :, :, 1))) < 1.0e-6_wp) err = max(err, 1.0_wp)
   allocate(density_image_large(nao, nao), overlap_image_large(nao, nao), &
      & source=(0.0_wp, 0.0_wp))
   do ik = 1, product(nmesh_large)
      angle = 2.0_wp*pi*kfrac_large(1, ik)
      phase = weights_large(ik)*exp(cmplx(0.0_wp, -angle, wp))
      density_image_large = density_image_large &
         & +phase*density_large(:, :, 1, ik)
      overlap_image_large = overlap_image_large+phase*overlap_large(:, :, ik)
   end do
   if (maxval(abs(density_image_large)) < 1.0e-6_wp &
      & .or. maxval(abs(overlap_image_large)) < 1.0e-6_wp) &
      & err = max(err, 1.0_wp)
   call new_wavefunction(wfn_large, mol%nat, nsh, nao, 2, 0.0_wp, .false.)
   wfn_large%qsh(:, 1) = wfn%qsh(:, 1)
   bvk_builds = cache_cp2k_plan%bvk_plan_builds
   call cp2k_exchange_kmesh(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn_large, nmesh_large, &
      & kfrac_large, weights_large, density_large, overlap_large, fock_large, &
      & vsh_large, energy_large, error)
   if (allocated(error)) return
   if (cache_cp2k_plan%bvk_plan_builds /= bvk_builds+1) err = max(err, 1.0_wp)
   kernel_large = cache_cp2k_plan%bvk_kernel
   select type (fock => exchange)
   type is (exchange_fock)
      call get_KFock_kmesh(fock, mol, ptr, oracle_plan, kernel_large, kfrac_large, &
         & weights_large, density_large, overlap_large, fock_large_oracle, &
         & vsh_large_oracle, energy_large_oracle)
   end select
   err = max(err, maxval(abs(fock_large-fock_large_oracle)), &
      & maxval(abs(vsh_large-vsh_large_oracle)), &
      & abs(energy_large-energy_large_oracle))
   ! The selectable separable direct-DFT backend must reproduce the unchanged
   ! dense phase-table oracle for the full UKS energy, Fock response, and
   ! charge-dependent shell potential.  This 9x9x1 mesh also exercises a
   ! common twist and the noncanonical input order constructed above.
   call cp2k_exchange_kmesh(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn_large, &
      & nmesh_large, kfrac_large, weights_large, density_large, &
      & overlap_large, fock_large, vsh_large, energy_large, error, &
      & transform_backend=cp2k_exchange_transform_separable_dft)
   if (allocated(error)) return
   err = max(err, maxval(abs(fock_large-fock_large_oracle)), &
      & maxval(abs(vsh_large-vsh_large_oracle)), &
      & abs(energy_large-energy_large_oracle))
   separable_err = max(separable_err, &
      & maxval(abs(fock_large-fock_large_oracle)), &
      & maxval(abs(vsh_large-vsh_large_oracle)), &
      & abs(energy_large-energy_large_oracle))
   call cp2k_exchange_kmesh(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn_large, &
      & nmesh_large, kfrac_large, weights_large, density_large, &
      & overlap_large, fock_large, vsh_large, energy_large, error, &
      & transform_backend=cp2k_exchange_transform_mixed_radix_fft)
   if (allocated(error)) return
   fft_err = max(fft_err, maxval(abs(fock_large-fock_large_oracle)), &
      & maxval(abs(vsh_large-vsh_large_oracle)), &
      & abs(energy_large-energy_large_oracle))
   print *, "MIXED_RADIX_FFT_UKS_FORWARD_RESID", &
      & maxval(abs(fock_large-fock_large_oracle)), &
      & maxval(abs(vsh_large-vsh_large_oracle)), &
      & abs(energy_large-energy_large_oracle)
   err = max(err, fft_err)
   if (allocated(cache_cp2k_plan%bvk_phase_forward) .or. &
      & allocated(cache_cp2k_plan%bvk_phase_inverse)) err = max(err, 1.0_wp)
   ! The production grid proof and reduced image accumulation are also tested
   ! beyond the small-mesh dense phase oracle, with reverse arrival order.
   call cp2k_exchange_stream_begin(stream_large, calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, &
      & wfn_large, nmesh_large, kfrac_large, weights_large, error)
   if (allocated(error)) return
   if (.not.cp2k_exchange_stream_has_full_mesh_storage(stream_large)) &
      & err = max(err, 1.0_wp)
   do j = 1, product(nmesh_large)
      ik = product(nmesh_large)-j+1
      call cp2k_exchange_stream_push(stream_large, ik, &
         & density_large(:, :, :, ik), overlap_large(:, :, ik), error)
      if (allocated(error)) return
   end do
   call cp2k_exchange_stream_apply(stream_large, calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, &
      & vsh_large, energy_large, error)
   if (allocated(error)) return
   if (.not.cp2k_exchange_stream_has_full_mesh_storage(stream_large)) &
      & err = max(err, 1.0_wp)
   do ik = 1, product(nmesh_large)
      call cp2k_exchange_stream_get(stream_large, ik, &
         & fock_large(:, :, :, ik), error, overlap_large(:, :, ik))
      if (allocated(error)) return
   end do
   call cp2k_exchange_stream_end(stream_large, error)
   if (allocated(error)) return
   err = max(err, maxval(abs(fock_large-fock_large_oracle)), &
      & maxval(abs(vsh_large-vsh_large_oracle)), &
      & abs(energy_large-energy_large_oracle))

   ! A non-divisor batch size on the 9x9x1 mesh exercises bounded storage,
   ! reverse k-block arrival, and a short final image batch beyond the dense
   ! small-grid phase oracle.
   fock_large = (0.0_wp, 0.0_wp)
   call cp2k_exchange_stream_begin(stream_batch_large, calc_cp2k, mol, &
      & cache_cp2k, cache_cp2k_plan, wfn_large, nmesh_large, kfrac_large, weights_large, error, &
      & mode=cp2k_exchange_stream_batched, batch_size=7)
   if (allocated(error)) return
   if (cp2k_exchange_stream_has_full_mesh_storage(stream_batch_large)) &
      & err = max(err, 1.0_wp)
   if (allocated(cache_cp2k_plan%bvk_phase_forward) &
      & .or. allocated(cache_cp2k_plan%bvk_phase_inverse)) err = max(err, 1.0_wp)
   peak_expected = 4_int64*int(nao, int64)*int(nao, int64) &
      & *7_int64*int(wfn_large%nspin, int64) &
      & +3_int64*int(nao, int64)*int(wfn_large%nspin, int64) &
      & +2_int64*int(product(nmesh_large), int64) &
      & +3_int64*int(maxval(nmesh_large), int64)
   peak_expected = peak_expected + max( &
      & 2_int64*int(nao, int64)**2, &
      & 5_int64*int(nao, int64)**2+6_int64*int(nao, int64), &
      & int(nao, int64)**2*int(wfn_large%nspin, int64) &
      & +4_int64*int(nao, int64)**2+7_int64)
   if (cp2k_exchange_stream_peak_complex_elements(stream_batch_large) &
      & /= peak_expected) err = max(err, 1.0_wp)
   batches_done = .false.
   do while (.not.batches_done)
      call cp2k_exchange_stream_batch_bounds(stream_batch_large, first_image, &
         & last_image, error)
      if (allocated(error)) return
      if (last_image-first_image+1 > 7) err = max(err, 1.0_wp)
      if (last_image == product(nmesh_large) &
         & .and. last_image-first_image+1 /= 4) err = max(err, 1.0_wp)
      if (cp2k_exchange_stream_peak_complex_elements(stream_batch_large) &
         & /= peak_expected) err = max(err, 1.0_wp)
      do j = 1, product(nmesh_large)
         ik = product(nmesh_large)-j+1
         call cp2k_exchange_stream_push(stream_batch_large, ik, &
            & density_large(:, :, :, ik), overlap_large(:, :, ik), error)
         if (allocated(error)) return
      end do
      call cp2k_exchange_stream_batch_apply(stream_batch_large, calc_cp2k, &
         & mol, cache_cp2k, cache_cp2k_plan, error)
      if (allocated(error)) return
      do j = 1, product(nmesh_large)
         ik = product(nmesh_large)-j+1
         call cp2k_exchange_stream_batch_pull(stream_batch_large, ik, &
            & overlap_large(:, :, ik), &
            & fock_large(:, :, :, ik), error)
         if (allocated(error)) return
      end do
      call cp2k_exchange_stream_batch_advance(stream_batch_large, &
         & batches_done, error)
      if (allocated(error)) return
   end do
   call cp2k_exchange_stream_batch_result(stream_batch_large, calc_cp2k, mol, &
      & cache_cp2k, cache_cp2k_plan, vsh_large, energy_large, error)
   if (allocated(error)) return
   call cp2k_exchange_stream_end(stream_batch_large, error)
   if (allocated(error)) return
   err = max(err, maxval(abs(fock_large-fock_large_oracle)), &
      & maxval(abs(vsh_large-vsh_large_oracle)), &
      & abs(energy_large-energy_large_oracle))

   allocate(overlap_grad_large(nao, nao, product(nmesh_large)), &
      & overlap_grad_large_oracle(nao, nao, product(nmesh_large)), &
      & mulliken_grad_large(nao, nao, product(nmesh_large)), &
      & bocorr_grad_large(nao, nao, product(nmesh_large)), &
      & gradient_large(3, mol%nat), gradient_large_oracle(3, mol%nat), &
      & sigma_large(3, 3), sigma_large_oracle(3, 3))
   select type (fock => exchange)
   type is (exchange_fock)
      call get_KGrad_kmesh(fock, mol, ptr, oracle_plan, kernel_large, kfrac_large, &
         & weights_large, density_large, overlap_large, &
         & overlap_grad_large_oracle, mulliken_grad_large, bocorr_grad_large)
      gradient_large_oracle = 0.0_wp
      sigma_large_oracle = 0.0_wp
      call fock%get_bvk_Kmatrix_derivs(mol, kernel_large, &
         & mulliken_grad_large, bocorr_grad_large, gradient_large_oracle, &
         & sigma_large_oracle)
   end select
   call cp2k_exchange_kmesh_gradient(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, &
      & nmesh_large, kfrac_large, weights_large, density_large, &
      & overlap_large, overlap_grad_large, gradient_large, sigma_large, error)
   if (allocated(error)) return
   err = max(err, maxval(abs(overlap_grad_large-overlap_grad_large_oracle)), &
      & maxval(abs(gradient_large-gradient_large_oracle)), &
      & maxval(abs(sigma_large-sigma_large_oracle)))
   ! Reverse mode uses the exact adjoint normalization of the same factorized
   ! transforms.  Atomic forces and homogeneous-strain response therefore
   ! provide end-to-end derivative oracles, not only transform-level checks.
   call cp2k_exchange_kmesh_gradient(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, &
      & nmesh_large, kfrac_large, weights_large, density_large, &
      & overlap_large, overlap_grad_large, gradient_large, sigma_large, error, &
      & transform_backend=cp2k_exchange_transform_separable_dft)
   if (allocated(error)) return
   err = max(err, maxval(abs(overlap_grad_large-overlap_grad_large_oracle)), &
      & maxval(abs(gradient_large-gradient_large_oracle)), &
      & maxval(abs(sigma_large-sigma_large_oracle)))
   separable_err = max(separable_err, &
      & maxval(abs(overlap_grad_large-overlap_grad_large_oracle)), &
      & maxval(abs(gradient_large-gradient_large_oracle)), &
      & maxval(abs(sigma_large-sigma_large_oracle)))
   call cp2k_exchange_kmesh_gradient(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, &
      & nmesh_large, kfrac_large, weights_large, density_large, &
      & overlap_large, overlap_grad_large, gradient_large, sigma_large, error, &
      & transform_backend=cp2k_exchange_transform_mixed_radix_fft)
   if (allocated(error)) return
   fft_err = max(fft_err, &
      & maxval(abs(overlap_grad_large-overlap_grad_large_oracle)), &
      & maxval(abs(gradient_large-gradient_large_oracle)), &
      & maxval(abs(sigma_large-sigma_large_oracle)))
   print *, "MIXED_RADIX_FFT_UKS_REVERSE_RESID", &
      & maxval(abs(overlap_grad_large-overlap_grad_large_oracle)), &
      & maxval(abs(gradient_large-gradient_large_oracle)), &
      & maxval(abs(sigma_large-sigma_large_oracle))

   ! UKS, a shifted/permuted 9x9x1 mesh, and a non-divisor batch exercise the
   ! bounded reverse protocol beyond the small RKS oracle above.
   overlap_grad_large = (0.0_wp, 0.0_wp)
   call cp2k_exchange_stream_begin(stream_batch_large, calc_cp2k, mol, &
      & cache_cp2k, cache_cp2k_plan, wfn_large, nmesh_large, kfrac_large, weights_large, error, &
      & mode=cp2k_exchange_stream_batched, batch_size=7, &
      & reverse_only=.true.)
   if (allocated(error)) return
   if (cp2k_exchange_stream_has_full_mesh_storage(stream_batch_large)) &
      & err = max(err, 1.0_wp)
   peak_expected = 6_int64*int(nao, int64)**2*7_int64 &
      & *int(wfn_large%nspin, int64)+5_int64*int(nao, int64) &
      & *int(wfn_large%nspin, int64) &
      & +2_int64*int(product(nmesh_large), int64) &
      & +3_int64*int(maxval(nmesh_large), int64) &
      & +10_int64*int(nao, int64)**2
   if (cp2k_exchange_stream_peak_complex_elements(stream_batch_large) &
      & /= peak_expected) err = max(err, 1.0_wp)
   batches_done = .false.
   do while (.not.batches_done)
      do j = 1, product(nmesh_large)
         ik = product(nmesh_large)-j+1
         call cp2k_exchange_stream_reverse_push(stream_batch_large, ik, &
            & density_large(:, :, :, ik), overlap_large(:, :, ik), error)
         if (allocated(error)) return
      end do
      call cp2k_exchange_stream_reverse_batch_apply(stream_batch_large, &
         & calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, error)
      if (allocated(error)) return
      do j = 1, product(nmesh_large)
         ik = product(nmesh_large)-j+1
         call cp2k_exchange_stream_reverse_batch_pull(stream_batch_large, ik, &
            & density_large(:, :, :, ik), overlap_large(:, :, ik), &
            & overlap_grad_large(:, :, ik), error)
         if (allocated(error)) return
      end do
      call cp2k_exchange_stream_reverse_batch_advance(stream_batch_large, &
         & batches_done, error)
      if (allocated(error)) return
   end do
   call cp2k_exchange_stream_reverse_result(stream_batch_large, calc_cp2k, &
      & mol, cache_cp2k, cache_cp2k_plan, gradient_large, sigma_large, error)
   if (allocated(error)) return
   peak_expected = cp2k_exchange_stream_peak_complex_elements( &
      & stream_batch_large)
   call cp2k_exchange_stream_end(stream_batch_large, error)
   if (allocated(error)) return
   err = max(err, maxval(abs(overlap_grad_large-overlap_grad_large_oracle)), &
      & maxval(abs(gradient_large-gradient_large_oracle)), &
      & maxval(abs(sigma_large-sigma_large_oracle)))
   stream_reverse_err = max(stream_reverse_err, &
      & maxval(abs(overlap_grad_large-overlap_grad_large_oracle)), &
      & maxval(abs(gradient_large-gradient_large_oracle)), &
      & maxval(abs(sigma_large-sigma_large_oracle)))

   print *, "STREAM_REVERSE_UKS_BATCH_RESID", &
      & maxval(abs(overlap_grad_large-overlap_grad_large_oracle)), &
      & maxval(abs(gradient_large-gradient_large_oracle)), &
      & maxval(abs(sigma_large-sigma_large_oracle)), peak_expected

   ! Exercise every whole-mesh selector in both directions for an RKS mesh
   ! that is genuinely two-dimensional, shifted, and externally permuted.
   ! This follows the UKS stream checks so the helper's RKS onsite state is
   ! also the correct state for the 3D RKS stream oracle below.
   allocate(density_large_rks(nao, nao, 1, product(nmesh_large)))
   density_large_rks(:, :, 1, :) = density_large(:, :, 1, :) &
      & +density_large(:, :, 2, :)
   call check_exchange_transform_selectors(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn, &
      & nmesh_large, kfrac_large, weights_large, density_large_rks, &
      & overlap_large, error, selector_separable_err, selector_fft_err, &
      & rks_2d_fft_forward, rks_2d_fft_reverse)
   if (allocated(error)) return
   separable_err = max(separable_err, selector_separable_err)
   fft_err = max(fft_err, selector_fft_err)
   print *, "TRANSFORM_SELECTOR_2D_RKS_SEPARABLE_MAX_RESID", &
      & selector_separable_err
   print *, "MIXED_RADIX_FFT_2D_RKS_FORWARD_RESID", rks_2d_fft_forward
   print *, "MIXED_RADIX_FFT_2D_RKS_REVERSE_RESID", rks_2d_fft_reverse

   ! A genuinely three-dimensional, shifted 2x2x2 grid completes the mesh
   ! dimensionality coverage of the 3x1x1 RKS and 9x9x1 UKS cases.
   allocate(kfrac_cube(3, product(nmesh_cube)), &
      & weights_cube(product(nmesh_cube)), &
      & density_cube(nao, nao, 1, product(nmesh_cube)), &
      & density_cube_uks(nao, nao, 2, product(nmesh_cube)), &
      & overlap_cube(nao, nao, product(nmesh_cube)), &
      & overlap_grad_cube(nao, nao, product(nmesh_cube)), &
      & overlap_grad_cube_oracle(nao, nao, product(nmesh_cube)), &
      & gradient_cube(3, mol%nat), gradient_cube_oracle(3, mol%nat), &
      & sigma_cube(3, 3), sigma_cube_oracle(3, 3))
   weights_cube = 1.0_wp/real(product(nmesh_cube), wp)
   do ik = 1, product(nmesh_cube)
      j = modulo(5*(ik-1)+3, product(nmesh_cube))
      kfrac_cube(1, ik) = 0.071_wp &
         & +real(modulo(j, nmesh_cube(1)), wp)/real(nmesh_cube(1), wp)
      kfrac_cube(2, ik) = 0.043_wp &
         & +real(modulo(j/nmesh_cube(1), nmesh_cube(2)), wp) &
         & /real(nmesh_cube(2), wp)
      kfrac_cube(3, ik) = 0.019_wp &
         & +real(j/(nmesh_cube(1)*nmesh_cube(2)), wp) &
         & /real(nmesh_cube(3), wp)
      i = 1+modulo(j, ncell)
      density_cube(:, :, 1, ik) = density_mesh(:, :, 1, i)
      overlap_cube(:, :, ik) = overlap_mesh(:, :, i)
      ii = 1+modulo(11*j+2, product(nmesh_large))
      density_cube_uks(:, :, :, ik) = density_large(:, :, :, ii)
   end do
   call cp2k_exchange_kmesh_gradient(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, nmesh_cube, &
      & kfrac_cube, weights_cube, density_cube, overlap_cube, &
      & overlap_grad_cube_oracle, gradient_cube_oracle, sigma_cube_oracle, &
      & error)
   if (allocated(error)) return
   overlap_grad_cube = (0.0_wp, 0.0_wp)
   call cp2k_exchange_stream_begin(stream_reverse_batch, calc_cp2k, mol, &
      & cache_cp2k, cache_cp2k_plan, wfn, nmesh_cube, kfrac_cube, weights_cube, error, &
      & mode=cp2k_exchange_stream_batched, batch_size=3, &
      & reverse_only=.true.)
   if (allocated(error)) return
   peak_expected = 6_int64*int(nao, int64)**2*3_int64 &
      & *int(wfn%nspin, int64)+5_int64*int(nao, int64) &
      & *int(wfn%nspin, int64)+2_int64*int(product(nmesh_cube), int64) &
      & +3_int64*int(maxval(nmesh_cube), int64)+10_int64*int(nao, int64)**2
   if (cp2k_exchange_stream_peak_complex_elements(stream_reverse_batch) &
      & /= peak_expected) err = max(err, 1.0_wp)
   batches_done = .false.
   do while (.not.batches_done)
      do ik = 1, product(nmesh_cube)
         call cp2k_exchange_stream_reverse_push(stream_reverse_batch, ik, &
            & density_cube(:, :, :, ik), overlap_cube(:, :, ik), error)
         if (allocated(error)) return
      end do
      call cp2k_exchange_stream_reverse_batch_apply(stream_reverse_batch, &
         & calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, error)
      if (allocated(error)) return
      do ik = 1, product(nmesh_cube)
         call cp2k_exchange_stream_reverse_batch_pull(stream_reverse_batch, &
            & ik, density_cube(:, :, :, ik), overlap_cube(:, :, ik), &
            & overlap_grad_cube(:, :, ik), error)
         if (allocated(error)) return
      end do
      call cp2k_exchange_stream_reverse_batch_advance(stream_reverse_batch, &
         & batches_done, error)
      if (allocated(error)) return
   end do
   call cp2k_exchange_stream_reverse_result(stream_reverse_batch, calc_cp2k, &
      & mol, cache_cp2k, cache_cp2k_plan, gradient_cube, sigma_cube, error)
   if (allocated(error)) return
   call cp2k_exchange_stream_end(stream_reverse_batch, error)
   if (allocated(error)) return
   err = max(err, maxval(abs(overlap_grad_cube-overlap_grad_cube_oracle)), &
      & maxval(abs(gradient_cube-gradient_cube_oracle)), &
      & maxval(abs(sigma_cube-sigma_cube_oracle)))
   stream_reverse_err = max(stream_reverse_err, &
      & maxval(abs(overlap_grad_cube-overlap_grad_cube_oracle)), &
      & maxval(abs(gradient_cube-gradient_cube_oracle)), &
      & maxval(abs(sigma_cube-sigma_cube_oracle)))
   print *, "STREAM_REVERSE_3D_RESID", &
      & maxval(abs(overlap_grad_cube-overlap_grad_cube_oracle)), &
      & maxval(abs(gradient_cube-gradient_cube_oracle)), &
      & maxval(abs(sigma_cube-sigma_cube_oracle)), peak_expected

   ! Complete the selector matrix with a shifted, permuted, genuinely 3D UKS
   ! mesh.  Both spin channels are nonzero and distinct.
   call check_exchange_transform_selectors(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, &
      & wfn_large, nmesh_cube, kfrac_cube, weights_cube, density_cube_uks, &
      & overlap_cube, error, selector_separable_err, selector_fft_err, &
      & uks_3d_fft_forward, uks_3d_fft_reverse)
   if (allocated(error)) return
   separable_err = max(separable_err, selector_separable_err)
   fft_err = max(fft_err, selector_fft_err)
   print *, "TRANSFORM_SELECTOR_3D_UKS_SEPARABLE_MAX_RESID", &
      & selector_separable_err
   print *, "MIXED_RADIX_FFT_3D_UKS_FORWARD_RESID", uks_3d_fft_forward
   print *, "MIXED_RADIX_FFT_3D_UKS_REVERSE_RESID", uks_3d_fft_reverse

   ! The public partial-accumulator ABI must reproduce the unchanged dense
   ! forward and reverse oracles for one, two, and four host partitions.  Four
   ! partitions on the three-point mesh deliberately include an empty owner.
   call check_partial_exchange_case(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn, nmesh, &
      & kfrac, weights, density_mesh, overlap_mesh, [1, 2, 4], .true., error, &
      & partial_err)
   if (allocated(error)) return
   stream_reverse_err = max(stream_reverse_err, partial_err)
   call check_partial_exchange_case(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn, nmesh, &
      & kfrac_order, weights_order, density_order, overlap_order, [4], &
      & .false., error, partial_err)
   if (allocated(error)) return
   stream_reverse_err = max(stream_reverse_err, partial_err)
   call check_partial_exchange_case(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn_large, &
      & nmesh_large, kfrac_large, weights_large, density_large, overlap_large, &
      & [4], .false., error, partial_err)
   if (allocated(error)) return
   stream_reverse_err = max(stream_reverse_err, partial_err)
   call check_partial_exchange_case(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn, &
      & nmesh_cube, kfrac_cube, weights_cube, density_cube, overlap_cube, &
      & [2], .false., error, partial_err)
   if (allocated(error)) return
   stream_reverse_err = max(stream_reverse_err, partial_err)
   call check_partial_exchange_multibatch(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn, &
      & nmesh, kfrac, weights, density_mesh, overlap_mesh, error, partial_err)
   if (allocated(error)) return
   stream_reverse_err = max(stream_reverse_err, partial_err)
   call check_partial_exchange_distributed_images(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, &
      & wfn_large, nmesh_large, kfrac_large, weights_large, density_large, &
      & overlap_large, error, partial_err)
   if (allocated(error)) return
   stream_reverse_err = max(stream_reverse_err, partial_err)
   print *, "DISTRIBUTED_IMAGE_KERNEL_9X9X1_UKS_MAX_RESID", partial_err
   call check_partial_exchange_range_local(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, &
      & wfn_large, nmesh_large, kfrac_large, weights_large, density_large, &
      & overlap_large, error, partial_err)
   if (allocated(error)) return
   stream_reverse_err = max(stream_reverse_err, partial_err)
   print *, "RANGE_LOCAL_IMAGE_IMPORT_9X9X1_UKS_MAX_RESID", partial_err
   call check_partial_exchange_range_local(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, &
      & wfn_large, nmesh_large, kfrac_large, weights_large, density_large, &
      & overlap_large, error, partial_err, &
      & cp2k_exchange_transform_mixed_radix_fft)
   if (allocated(error)) return
   stream_reverse_err = max(stream_reverse_err, partial_err)
   print *, "RANGE_LOCAL_MIXED_RADIX_9X9X1_UKS_MAX_RESID", partial_err
   call check_partial_exchange_distributed_images(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, &
      & wfn, nmesh_cube, kfrac_cube, weights_cube, density_cube, overlap_cube, &
      & error, partial_err)
   if (allocated(error)) return
   stream_reverse_err = max(stream_reverse_err, partial_err)
   print *, "DISTRIBUTED_IMAGE_KERNEL_2X2X2_MAX_RESID", partial_err
   call check_partial_exchange_range_local(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn, &
      & nmesh_cube, kfrac_cube, weights_cube, density_cube, overlap_cube, &
      & error, partial_err)
   if (allocated(error)) return
   stream_reverse_err = max(stream_reverse_err, partial_err)
   print *, "RANGE_LOCAL_IMAGE_IMPORT_2X2X2_MAX_RESID", partial_err
   call check_partial_exchange_range_local(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn, &
      & nmesh_cube, kfrac_cube, weights_cube, density_cube, overlap_cube, &
      & error, partial_err, cp2k_exchange_transform_mixed_radix_fft)
   if (allocated(error)) return
   stream_reverse_err = max(stream_reverse_err, partial_err)
   print *, "RANGE_LOCAL_MIXED_RADIX_2X2X2_MAX_RESID", partial_err
   call check_partial_exchange_distributed_images(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, &
      & wfn, nmesh, kfrac, weights, density_mesh, overlap_mesh, error, &
      & partial_err)
   if (allocated(error)) return
   stream_reverse_err = max(stream_reverse_err, partial_err)
   print *, "DISTRIBUTED_IMAGE_KERNEL_3X1X1_MAX_RESID", partial_err
   call check_partial_exchange_range_local(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn, &
      & nmesh, kfrac, weights, density_mesh, overlap_mesh, error, partial_err)
   if (allocated(error)) return
   stream_reverse_err = max(stream_reverse_err, partial_err)
   print *, "RANGE_LOCAL_IMAGE_IMPORT_3X1X1_MAX_RESID", partial_err
   print *, "PARTIAL_K_TO_R_MAX_RESID", stream_reverse_err

   bvk_builds = cache_cp2k_plan%bvk_plan_builds
   call cp2k_exchange_kmesh(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, wfn_large, nmesh_large, &
      & kfrac_large, weights_large, density_large, overlap_large, fock_large, &
      & vsh_large, energy_large, error)
   if (allocated(error)) return
   if (cache_cp2k_plan%bvk_plan_builds /= bvk_builds+1) err = max(err, 1.0_wp)
   err = max(err, maxval(abs(fock_large-fock_large_oracle)), &
      & maxval(abs(vsh_large-vsh_large_oracle)), &
      & abs(energy_large-energy_large_oracle))

   call check(error, separable_err, 0.0_wp, thr=1.0e-11_wp)
   if (allocated(error)) then
      call test_failed(error, &
         & "Separable exchange exceeds the dense-oracle tolerance")
      return
   end if
   print *, "MIXED_RADIX_FFT_EXCHANGE_MAX_RESID", fft_err
   call check(error, fft_err, 0.0_wp, thr=1.0e-11_wp)
   if (allocated(error)) then
      call test_failed(error, &
         & "Mixed-radix FFT exchange exceeds the dense-oracle tolerance")
      return
   end if
   call check(error, stream_reverse_err, 0.0_wp, thr=1.0e-11_wp)
   if (allocated(error)) then
      call test_failed(error, &
         & "Streamed reverse exchange exceeds the dense-oracle tolerance")
      return
   end if
   print *, "BVK_EXCHANGE_SUPERCELL_MAX_RESID", err
   call check(error, err, 0.0_wp, thr=1.0e-9_wp)
   if (allocated(error)) then
      call test_failed(error, &
         & "Whole-mesh exchange does not reproduce the explicit BvK supercell")
   end if
end subroutine test_bvk_exchange_kernel

subroutine test_kpoint_fock_complex(error)
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(basis_type), allocatable :: bas
   class(qvszp_basis_type), allocatable :: qvszp_bas
   class(exchange_type), allocatable :: exchange
   type(xtb_calculator) :: calc_cp2k
   type(basis_cache) :: bcache
   type(container_cache) :: cache, cache_cp2k
   type(cp2k_stream_plan_type) :: cache_cp2k_plan
   type(exchange_cache), pointer :: ptr
   type(exchange_bvk_kernel) :: kernel_mesh
   type(cp2k_stream_plan_type) :: oracle_plan
   type(wavefunction_type) :: wfn, wfn_aux
   integer, parameter :: num(3) = [1, 1, 8]
   real(wp), parameter :: xyz(3, 3) = reshape([ &
      & 0.20_wp, 0.30_wp, 0.40_wp, &
      & 0.96_wp, 0.30_wp, 0.90_wp, &
      & 0.20_wp, 0.30_wp, 1.40_wp], [3, 3])/autoaa
   real(wp), parameter :: lattice(3, 3) = reshape([ &
      & 7.80_wp, 0.20_wp, 0.10_wp, &
      & 0.40_wp, 8.30_wp, 0.30_wp, &
      & 0.10_wp, 0.50_wp, 9.10_wp], [3, 3])/autoaa
   logical, parameter :: periodic(3) = [.true., .true., .true.]
   integer :: i, iao, iat, ii, is, ish, izp, j, jao, jat, jj, js, jsh, jzp, nao
   real(wp) :: deriv_fd, deriv_ref, energy, energy_m, energy_p, energy_ref
   real(wp) :: err, q0, step, theta
   real(wp) :: coord_direction(3, 3), identity(3, 3), lattice0(3, 3), &
      & strain_direction(3, 3), transform(3, 3), xyz0(3, 3)
   real(wp), allocatable :: lattr(:, :), overlap(:, :), vsh(:), vsh_m(:), &
      & vsh_p(:), vsh_tr(:), vsh_u(:), ao_grad_gamma(:, :), &
      & mulliken_grad(:, :), mulliken_grad_gamma(:, :), bocorr_grad(:, :), &
      & bocorr_grad_gamma(:, :), bocorr_grad_mesh_r(:, :, :), &
      & gradient_cp_k(:, :), gradient_cp_mesh(:, :), gradient_gamma(:, :), &
      & gradient_k(:, :), gradient_mesh(:, :), &
      & kernel_backup(:, :), kernel_direction(:, :), sigma_gamma(:, :), sigma_k(:, :)
   real(wp), allocatable :: mulliken_grad_mesh_r(:, :, :), &
      & sigma_cp_k(:, :), sigma_cp_mesh(:, :), sigma_mesh(:, :), vsh_mesh(:)
   complex(wp), allocatable :: density_k(:, :, :), density_m(:, :, :), &
      & density_p(:, :, :), density_u(:, :, :), direction(:, :, :), &
      & density_uks(:, :, :), &
      & fock_k(:, :, :), fock_m(:, :, :), fock_p(:, :, :), fock_tr(:, :, :), &
      & fock_u(:, :, :), fock_uks(:, :, :), overlap_grad(:, :), &
      & overlap_grad_cp_k(:, :), &
      & overlap_grad_uks(:, :), overlap_k(:, :), overlap_m(:, :), &
      & overlap_p(:, :), overlap_u(:, :), unitary(:)
   complex(wp), allocatable :: density_mesh(:, :, :, :), &
      & fock_mesh(:, :, :, :), overlap_grad_cp_mesh(:, :, :), &
      & overlap_grad_mesh(:, :, :), &
      & overlap_mesh(:, :, :)

   call new(mol, num, xyz, lattice=lattice, periodic=periodic)
   allocate(qvszp_bas)
   call make_qvszp_basis(qvszp_bas, mol, error)
   if (allocated(error)) return
   call move_alloc(qvszp_bas, bas)

   call new_wavefunction(wfn_aux, mol%nat, 0, 0, 1, 0.0_wp, .false.)
   call get_eeqbc_charges(mol, error, wfn_aux%qat(:, 1))
   if (allocated(error)) return
   call bas%update(mol, bcache, .false., wfn_aux=wfn_aux)

   nao = bas%nao
   call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
   allocate(overlap(nao, nao))
   call get_overlap(mol, lattr, cutoff, bas, bcache, overlap)
   call new_wavefunction(wfn, mol%nat, bas%nsh, nao, 1, 0.0_wp, .false.)
   do i = 1, nao
      do j = 1, nao
         wfn%density(i, j, 1) = 0.2_wp/(1.0_wp + abs(i-j))
      end do
   end do
   wfn%qsh = 0.0_wp

   call make_exchange_gxtb(exchange, mol, bas)
   call exchange%update(mol, cache)
   call view(cache, ptr)
   call exchange%get_onsite_Kmatrix(mol, wfn, ptr)
   allocate(ptr%prev_F(nao, nao, 1), ptr%prev_vsh(bas%nsh, 1))
   call exchange%get_KFock(mol, ptr, wfn%density, overlap)
   allocate(calc_cp2k%exchange, source=exchange)
   call cp2k_prepare_exchange(calc_cp2k, mol, cache_cp2k, wfn, error)
   if (allocated(error)) return

   allocate(density_k(nao, nao, 1), overlap_k(nao, nao), &
      & fock_k(nao, nao, 1), vsh(bas%nsh))
   density_k = cmplx(wfn%density, 0.0_wp, wp)
   overlap_k = cmplx(overlap, 0.0_wp, wp)
   select type (fock => exchange)
   type is (exchange_fock)
      call get_KFock_kpoint(fock, mol, ptr, density_k, overlap_k, &
         & fock_k, vsh, energy)
      energy_ref = 0.5_wp*sum(wfn%density*ptr%prev_F)
      err = max(maxval(abs(real(fock_k, wp)-ptr%prev_F)), &
         & maxval(abs(aimag(fock_k))), maxval(abs(vsh-ptr%prev_vsh(:, 1))), &
         & abs(energy-energy_ref))
      call check(error, err, 0.0_wp, thr=1.0e-11_wp)
      if (allocated(error)) then
         call test_failed(error, "Complex exchange does not reproduce Gamma")
         return
      end if

      call fock%get_bvk_Kmatrix(mol, [1, 1, 1], kernel_mesh)
      allocate(density_mesh(nao, nao, 1, 1), overlap_mesh(nao, nao, 1), &
         & fock_mesh(nao, nao, 1, 1), vsh_mesh(bas%nsh))
      density_mesh(:, :, :, 1) = density_k
      overlap_mesh(:, :, 1) = overlap_k
      call get_KFock_kmesh(fock, mol, ptr, oracle_plan, kernel_mesh, &
         & reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]), [1.0_wp], &
         & density_mesh, overlap_mesh, fock_mesh, vsh_mesh, energy_ref)
      err = max(maxval(abs(fock_mesh(:, :, :, 1)-fock_k)), &
         & maxval(abs(vsh_mesh-vsh)), abs(energy_ref-energy))
      call check(error, err, 0.0_wp, thr=1.0e-11_wp)
      if (allocated(error)) then
         call test_failed(error, "Whole-mesh exchange does not reproduce Gamma")
         return
      end if

      ! At Gamma the direct overlap derivative and tblite's historical
      ! energy-weighted-density correction must obey G_S = -W_exchange.
      ! CP2K's independent-pair Pulay kernel supplies its own factor of two,
      ! so this relation fixes the sign used when folding the response.
      allocate(ao_grad_gamma(nao, nao), mulliken_grad(nao, nao), &
         & mulliken_grad_gamma(nao, nao), bocorr_grad(nao, nao), &
         & bocorr_grad_gamma(nao, nao), gradient_gamma(3, mol%nat), &
         & gradient_k(3, mol%nat), sigma_gamma(3, 3), sigma_k(3, 3), &
         & source=0.0_wp)
      allocate(overlap_grad(nao, nao), source=(0.0_wp, 0.0_wp))
      call fock%get_KGrad(mol, ptr, wfn%density, overlap, &
         & mulliken_grad_gamma, bocorr_grad_gamma, ao_grad_gamma)
      call get_KGrad_kpoint(fock, mol, ptr, density_k, overlap_k, &
         & mulliken_grad, bocorr_grad, overlap_grad)

      call fock%get_mulliken_derivs(mol, ptr, mulliken_grad_gamma, &
         & gradient_gamma, sigma_gamma)
      call fock%get_bocorr_derivs(mol, ptr, bocorr_grad_gamma, &
         & gradient_gamma, sigma_gamma)
      call get_mulliken_derivs_direct(fock, mol, ptr, mulliken_grad, &
         & gradient_k, sigma_k)
      call get_bocorr_derivs_direct(fock, mol, ptr, bocorr_grad, &
         & gradient_k, sigma_k)
      err = max(maxval(abs(overlap_grad+ao_grad_gamma)), &
         & maxval(abs(gradient_k-gradient_gamma)), &
         & maxval(abs(sigma_k-sigma_gamma)))

      ! The whole-mesh reverse mode must reduce exactly to the independently
      ! validated single-k/Gamma implementation at N=1.
      allocate(overlap_grad_mesh(nao, nao, 1), &
         & mulliken_grad_mesh_r(nao, nao, 1), &
         & bocorr_grad_mesh_r(nao, nao, 1), gradient_mesh(3, mol%nat), &
         & sigma_mesh(3, 3))
      call get_KGrad_kmesh(fock, mol, ptr, oracle_plan, kernel_mesh, &
         & reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]), [1.0_wp], &
         & density_mesh, overlap_mesh, overlap_grad_mesh, &
         & mulliken_grad_mesh_r, bocorr_grad_mesh_r)
      gradient_mesh = 0.0_wp
      sigma_mesh = 0.0_wp
      call fock%get_bvk_Kmatrix_derivs(mol, kernel_mesh, &
         & mulliken_grad_mesh_r, bocorr_grad_mesh_r, gradient_mesh, sigma_mesh)
      err = max(err, maxval(abs(overlap_grad_mesh(:, :, 1)-overlap_grad)), &
         & maxval(abs(gradient_mesh-gradient_k)), &
         & maxval(abs(sigma_mesh-sigma_k)))

      ! Exercise the public CP2K ABI as well: both wrappers must expose the
      ! same N=1 functional, including overlap sign and strain convention.
      allocate(overlap_grad_cp_k(nao, nao), &
         & overlap_grad_cp_mesh(nao, nao, 1), &
         & gradient_cp_k(3, mol%nat), gradient_cp_mesh(3, mol%nat), &
         & sigma_cp_k(3, 3), sigma_cp_mesh(3, 3))
      call cp2k_exchange_kpoint_gradient(calc_cp2k, mol, cache_cp2k, &
         & density_k, overlap_k, overlap_grad_cp_k, gradient_cp_k, &
         & sigma_cp_k, error)
      if (allocated(error)) return
      call cp2k_exchange_kmesh_gradient(calc_cp2k, mol, cache_cp2k, cache_cp2k_plan, &
         & [1, 1, 1], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]), &
         & [1.0_wp], density_mesh, overlap_mesh, overlap_grad_cp_mesh, &
         & gradient_cp_mesh, sigma_cp_mesh, error)
      if (allocated(error)) return
      err = max(err, maxval(abs(overlap_grad_cp_k-overlap_grad)), &
         & maxval(abs(overlap_grad_cp_mesh(:, :, 1)-overlap_grad_cp_k)), &
         & maxval(abs(gradient_cp_mesh-gradient_cp_k)), &
         & maxval(abs(sigma_cp_mesh-sigma_cp_k)))
      call check(error, err, 0.0_wp, thr=1.0e-10_wp)
      if (allocated(error)) then
         call test_failed(error, "Complex exchange gradient does not reproduce Gamma")
         return
      end if

      do i = 1, nao
         do j = 1, i-1
            overlap_k(j, i) = overlap_k(j, i) &
               & + cmplx(0.0_wp, 1.0e-3_wp*sin(real(3*i+5*j, wp)), wp)
            overlap_k(i, j) = conjg(overlap_k(j, i))
            density_k(j, i, 1) = density_k(j, i, 1) &
               & + cmplx(0.0_wp, 2.0e-3_wp*cos(real(7*i+2*j, wp)), wp)
            density_k(i, j, 1) = conjg(density_k(j, i, 1))
         end do
      end do
      call get_KFock_kpoint(fock, mol, ptr, density_k, overlap_k, &
         & fock_k, vsh, energy)
      err = maxval(abs(fock_k(:, :, 1)-conjg(transpose(fock_k(:, :, 1)))))

      allocate(fock_tr(nao, nao, 1), vsh_tr(bas%nsh))
      call get_KFock_kpoint(fock, mol, ptr, conjg(density_k), &
         & conjg(overlap_k), fock_tr, vsh_tr, energy_ref)
      err = max(err, maxval(abs(fock_tr-conjg(fock_k))), &
         & maxval(abs(vsh_tr-vsh)), abs(energy_ref-energy))

      allocate(unitary(nao), density_u(nao, nao, 1), overlap_u(nao, nao), &
         & fock_u(nao, nao, 1), vsh_u(bas%nsh))
      do i = 1, nao
         theta = 0.071_wp*real(i, wp)
         unitary(i) = exp(cmplx(0.0_wp, theta, wp))
      end do
      do i = 1, nao
         do j = 1, nao
            overlap_u(i, j) = conjg(unitary(i))*overlap_k(i, j)*unitary(j)
            density_u(i, j, 1) = conjg(unitary(i)) &
               & * density_k(i, j, 1)*unitary(j)
         end do
      end do
      call get_KFock_kpoint(fock, mol, ptr, density_u, overlap_u, &
         & fock_u, vsh_u, energy_ref)
      do i = 1, nao
         do j = 1, nao
            err = max(err, abs(fock_u(i, j, 1) &
               & - conjg(unitary(i))*fock_k(i, j, 1)*unitary(j)))
         end do
      end do
      err = max(err, maxval(abs(vsh_u-vsh)), abs(energy_ref-energy))

      allocate(direction(nao, nao, 1), density_m(nao, nao, 1), &
         & density_p(nao, nao, 1), fock_m(nao, nao, 1), &
         & fock_p(nao, nao, 1), source=(0.0_wp, 0.0_wp))
      allocate(vsh_m(bas%nsh), vsh_p(bas%nsh))
      do i = 1, nao
         direction(i, i, 1) = cmplx(1.0e-2_wp*sin(real(11*i, wp)), &
            & 0.0_wp, wp)
         do j = 1, i-1
            direction(j, i, 1) = cmplx( &
               & 1.0e-2_wp*cos(real(3*i+j, wp)), &
               & 1.0e-2_wp*sin(real(i+4*j, wp)), wp)
            direction(i, j, 1) = conjg(direction(j, i, 1))
         end do
      end do
      step = 1.0e-5_wp
      density_m = density_k-step*direction
      density_p = density_k+step*direction
      call get_KFock_kpoint(fock, mol, ptr, density_m, overlap_k, &
         & fock_m, vsh_m, energy_m)
      call get_KFock_kpoint(fock, mol, ptr, density_p, overlap_k, &
         & fock_p, vsh_p, energy_p)
      deriv_fd = (energy_p-energy_m)/(2.0_wp*step)
      deriv_ref = real(sum(conjg(direction)*fock_k), wp)
      err = max(err, abs(deriv_fd-deriv_ref))

      ! Reverse the complex exchange functional with respect to S(k).  The
      ! overlap response is the missing Pulay ingredient of a k-point force;
      ! validate it independently of any real-space integral convention.
      allocate(overlap_m(nao, nao), overlap_p(nao, nao))
      call get_KGrad_kpoint(fock, mol, ptr, density_k, overlap_k, &
         & mulliken_grad, bocorr_grad, overlap_grad)
      overlap_m = overlap_k-step*direction(:, :, 1)
      overlap_p = overlap_k+step*direction(:, :, 1)
      call get_KFock_kpoint(fock, mol, ptr, density_k, overlap_m, &
         & fock_m, vsh_m, energy_m)
      call get_KFock_kpoint(fock, mol, ptr, density_k, overlap_p, &
         & fock_p, vsh_p, energy_p)
      deriv_fd = (energy_p-energy_m)/(2.0_wp*step)
      deriv_ref = real(sum(conjg(overlap_grad)*direction(:, :, 1)), wp)
      err = max(err, abs(deriv_fd-deriv_ref))
      err = max(err, maxval(abs(overlap_grad &
         & - conjg(transpose(overlap_grad)))))

      ! The geometry-dependent Mulliken kernel is real and symmetric.  Check
      ! the AO-resolved reverse response before it is contracted with the
      ! pairwise kernel derivative in get_mulliken_derivs.
      allocate(kernel_backup(fock%nsh, fock%nsh), &
         & kernel_direction(fock%nsh, fock%nsh), source=0.0_wp)
      do ish = 1, fock%nsh
         kernel_direction(ish, ish) = 1.0e-2_wp*sin(real(13*ish, wp))
         do jsh = 1, ish-1
            kernel_direction(jsh, ish) = &
               & 1.0e-2_wp*cos(real(5*ish+7*jsh, wp))
            kernel_direction(ish, jsh) = kernel_direction(jsh, ish)
         end do
      end do
      kernel_backup = ptr%g_mulliken
      ptr%g_mulliken = kernel_backup-step*kernel_direction
      call get_KFock_kpoint(fock, mol, ptr, density_k, overlap_k, &
         & fock_m, vsh_m, energy_m)
      ptr%g_mulliken = kernel_backup+step*kernel_direction
      call get_KFock_kpoint(fock, mol, ptr, density_k, overlap_k, &
         & fock_p, vsh_p, energy_p)
      ptr%g_mulliken = kernel_backup
      deriv_fd = (energy_p-energy_m)/(2.0_wp*step)
      deriv_ref = 0.0_wp
      do ish = 1, fock%nsh
         ii = fock%iao_sh(ish)
         do jsh = 1, ish
            jj = fock%iao_sh(jsh)
            do iao = 1, fock%nao_sh(ish)
               do jao = 1, fock%nao_sh(jsh)
                  deriv_ref = deriv_ref + kernel_direction(jsh, ish) &
                     & * mulliken_grad(jj+jao, ii+iao)
               end do
            end do
         end do
      end do
      err = max(err, abs(deriv_fd-deriv_ref))

      ! Repeat the same directional check for the atom-block bond-order
      ! correction kernel.  As above, the response stores one derivative for
      ! each independent symmetric block, so only one triangle is contracted.
      deallocate(kernel_backup, kernel_direction)
      allocate(kernel_backup(mol%nat, mol%nat), &
         & kernel_direction(mol%nat, mol%nat), source=0.0_wp)
      do iat = 1, mol%nat
         do jat = 1, iat-1
            kernel_direction(jat, iat) = &
               & 1.0e-2_wp*sin(real(3*iat+11*jat, wp))
            kernel_direction(iat, jat) = kernel_direction(jat, iat)
         end do
      end do
      kernel_backup = ptr%g_bocorr
      ptr%g_bocorr = kernel_backup-step*kernel_direction
      call get_KFock_kpoint(fock, mol, ptr, density_k, overlap_k, &
         & fock_m, vsh_m, energy_m)
      ptr%g_bocorr = kernel_backup+step*kernel_direction
      call get_KFock_kpoint(fock, mol, ptr, density_k, overlap_k, &
         & fock_p, vsh_p, energy_p)
      ptr%g_bocorr = kernel_backup
      deriv_fd = (energy_p-energy_m)/(2.0_wp*step)
      deriv_ref = 0.0_wp
      do iat = 1, mol%nat
         izp = mol%id(iat)
         is = fock%ish_at(iat)
         do jat = 1, iat-1
            jzp = mol%id(jat)
            js = fock%ish_at(jat)
            do ish = 1, fock%nsh_id(izp)
               ii = fock%iao_sh(is+ish)
               do jsh = 1, fock%nsh_id(jzp)
                  jj = fock%iao_sh(js+jsh)
                  do iao = 1, fock%nao_sh(is+ish)
                     do jao = 1, fock%nao_sh(js+jsh)
                        deriv_ref = deriv_ref + kernel_direction(jat, iat) &
                           & * bocorr_grad(jj+jao, ii+iao)
                     end do
                  end do
               end do
            end do
         end do
      end do
      err = max(err, abs(deriv_fd-deriv_ref))

      ! Contract the direct kernel responses with their real-space derivatives
      ! and validate the complete explicit exchange geometry response for a
      ! genuinely complex Hermitian k-point density.  The overlap is held
      ! fixed here, so the separate overlap-response check above remains the
      ! sole Pulay contribution.
      gradient_k = 0.0_wp
      sigma_k = 0.0_wp
      call get_mulliken_derivs_direct(fock, mol, ptr, mulliken_grad, &
         & gradient_k, sigma_k)
      call get_bocorr_derivs_direct(fock, mol, ptr, bocorr_grad, &
         & gradient_k, sigma_k)

      do iat = 1, mol%nat
         do i = 1, 3
            coord_direction(i, iat) = 0.03_wp*sin(real(7*i+11*iat, wp))
         end do
      end do
      xyz0 = mol%xyz
      mol%xyz = xyz0-step*coord_direction
      call ptr%update(mol)
      call fock%get_mulliken_Kmatrix(mol, ptr)
      call fock%get_bocorr_Kmatrix(mol, ptr)
      call get_KFock_kpoint(fock, mol, ptr, density_k, overlap_k, &
         & fock_m, vsh_m, energy_m)
      mol%xyz = xyz0+step*coord_direction
      call ptr%update(mol)
      call fock%get_mulliken_Kmatrix(mol, ptr)
      call fock%get_bocorr_Kmatrix(mol, ptr)
      call get_KFock_kpoint(fock, mol, ptr, density_k, overlap_k, &
         & fock_p, vsh_p, energy_p)
      deriv_fd = (energy_p-energy_m)/(2.0_wp*step)
      deriv_ref = sum(gradient_k*coord_direction)
      err = max(err, abs(deriv_fd-deriv_ref))

      ! Apply a symmetric homogeneous strain to both Cartesian coordinates and
      ! lattice vectors.  This also covers diagonal self-image kernels, whose
      ! historical derivative representation carries a factor of two.
      identity = 0.0_wp
      do i = 1, 3
         identity(i, i) = 1.0_wp
      end do
      strain_direction = reshape([ &
         & 0.07_wp, 0.02_wp, -0.01_wp, &
         & 0.02_wp, -0.04_wp, 0.03_wp, &
         & -0.01_wp, 0.03_wp, 0.05_wp], [3, 3])
      lattice0 = mol%lattice
      mol%xyz = xyz0
      transform = identity-step*strain_direction
      mol%xyz = matmul(transform, xyz0)
      mol%lattice = matmul(transform, lattice0)
      call ptr%update(mol)
      call fock%get_mulliken_Kmatrix(mol, ptr)
      call fock%get_bocorr_Kmatrix(mol, ptr)
      call get_KFock_kpoint(fock, mol, ptr, density_k, overlap_k, &
         & fock_m, vsh_m, energy_m)
      transform = identity+step*strain_direction
      mol%xyz = matmul(transform, xyz0)
      mol%lattice = matmul(transform, lattice0)
      call ptr%update(mol)
      call fock%get_mulliken_Kmatrix(mol, ptr)
      call fock%get_bocorr_Kmatrix(mol, ptr)
      call get_KFock_kpoint(fock, mol, ptr, density_k, overlap_k, &
         & fock_p, vsh_p, energy_p)
      deriv_fd = (energy_p-energy_m)/(2.0_wp*step)
      deriv_ref = sum(sigma_k*strain_direction)
      err = max(err, abs(deriv_fd-deriv_ref))

      mol%xyz = xyz0
      mol%lattice = lattice0
      call ptr%update(mol)
      call fock%get_mulliken_Kmatrix(mol, ptr)
      call fock%get_bocorr_Kmatrix(mol, ptr)

      ! The shell potential is the derivative of the same complex-k energy
      ! with respect to the charge-dependent onsite exchange kernels.
      q0 = wfn%qsh(1, 1)
      wfn%qsh(1, 1) = q0-step
      call exchange%get_onsite_Kmatrix(mol, wfn, ptr)
      call get_KFock_kpoint(fock, mol, ptr, density_k, overlap_k, &
         & fock_m, vsh_m, energy_m)
      wfn%qsh(1, 1) = q0+step
      call exchange%get_onsite_Kmatrix(mol, wfn, ptr)
      call get_KFock_kpoint(fock, mol, ptr, density_k, overlap_k, &
         & fock_p, vsh_p, energy_p)
      err = max(err, abs((energy_p-energy_m)/(2.0_wp*step)-vsh(1)))
      wfn%qsh(1, 1) = q0
      call exchange%get_onsite_Kmatrix(mol, wfn, ptr)

      ! Exercise the unrestricted spin factor and accumulation of both spin
      ! channels in the common overlap response.
      allocate(density_uks(nao, nao, 2), fock_uks(nao, nao, 2), &
         & overlap_grad_uks(nao, nao))
      density_uks(:, :, 1) = 0.61_wp*density_k(:, :, 1)
      density_uks(:, :, 2) = 0.37_wp*density_k(:, :, 1)
      call get_KFock_kpoint(fock, mol, ptr, density_uks, overlap_k, &
         & fock_uks, vsh, energy)
      call get_KGrad_kpoint(fock, mol, ptr, density_uks, overlap_k, &
         & mulliken_grad, bocorr_grad, overlap_grad_uks)
      overlap_m = overlap_k-step*direction(:, :, 1)
      overlap_p = overlap_k+step*direction(:, :, 1)
      call get_KFock_kpoint(fock, mol, ptr, density_uks, overlap_m, &
         & fock_uks, vsh_m, energy_m)
      call get_KFock_kpoint(fock, mol, ptr, density_uks, overlap_p, &
         & fock_uks, vsh_p, energy_p)
      deriv_fd = (energy_p-energy_m)/(2.0_wp*step)
      deriv_ref = real(sum(conjg(overlap_grad_uks) &
         & * direction(:, :, 1)), wp)
      err = max(err, abs(deriv_fd-deriv_ref))
   class default
      call test_failed(error, "Unexpected exchange implementation")
      return
   end select

   call check(error, err, 0.0_wp, thr=1.0e-9_wp)
   if (allocated(error)) then
      call test_failed(error, &
         & "Complex exchange violates Hermiticity, covariance, or variationality")
   end if
end subroutine test_kpoint_fock_complex

subroutine test_regular_separable_dft(error)
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: ncase = 3, nrow = 5
   integer, parameter :: meshes(3, ncase) = reshape([ &
      & 7, 1, 1, &
      & 3, 4, 1, &
      & 2, 3, 4], [3, ncase])
   integer :: g(3), icase, icell, idim, igrid, ik, j, nk, root
   integer, allocatable :: input_to_grid(:), reps(:, :)
   real(wp), parameter :: twist(3) = [0.071_wp, -0.113_wp, 0.137_wp]
   real(wp) :: angle, err
   complex(wp) :: lhs, phase, rhs
   complex(wp), allocatable :: dense_k(:, :), dense_r(:, :), &
      & mesh_k(:, :), mesh_r(:, :), phase_roots(:, :), sep_k(:, :), &
      & sep_r(:, :), twist_phase(:)

   err = 0.0_wp
   do icase = 1, ncase
      nk = product(meshes(:, icase))
      allocate(input_to_grid(nk), reps(3, nk))
      allocate(phase_roots(maxval(meshes(:, icase)), 3), twist_phase(nk), &
         & dense_k(nrow, nk), dense_r(nrow, nk), mesh_k(nrow, nk), &
         & mesh_r(nrow, nk), sep_k(nrow, nk), sep_r(nrow, nk), &
         & source=(0.0_wp, 0.0_wp))

      ! A cyclically shifted k ordering and reversed representative ordering
      ! exercise both independent permutations of the compact plan.
      do ik = 1, nk
         input_to_grid(ik) = 1+modulo(ik+1, nk)
      end do
      do icell = 1, nk
         igrid = nk-icell
         reps(1, icell) = modulo(igrid, meshes(1, icase))
         reps(2, icell) = modulo(igrid/meshes(1, icase), &
            & meshes(2, icase))
         reps(3, icell) = igrid/(meshes(1, icase)*meshes(2, icase))
         angle = 2.0_wp*pi*dot_product(twist, &
            & real(reps(:, icell), wp))
         twist_phase(icell) = exp(cmplx(0.0_wp, angle, wp))
      end do
      phase_roots = (1.0_wp, 0.0_wp)
      do idim = 1, 3
         do root = 0, meshes(idim, icase)-1
            angle = 2.0_wp*pi*real(root, wp) &
               & /real(meshes(idim, icase), wp)
            phase_roots(root+1, idim) = exp(cmplx(0.0_wp, angle, wp))
         end do
      end do
      do ik = 1, nk
         do j = 1, nrow
            mesh_k(j, ik) = cmplx( &
               & sin(real(7*j+11*ik+3*icase, wp)), &
               & cos(real(13*j+5*ik+2*icase, wp)), wp)
            mesh_r(j, ik) = cmplx( &
               & cos(real(3*j+17*ik+icase, wp)), &
               & sin(real(19*j+2*ik+4*icase, wp)), wp)
         end do
      end do

      dense_r = (0.0_wp, 0.0_wp)
      dense_k = (0.0_wp, 0.0_wp)
      do ik = 1, nk
         igrid = input_to_grid(ik)-1
         g(1) = modulo(igrid, meshes(1, icase))
         g(2) = modulo(igrid/meshes(1, icase), meshes(2, icase))
         g(3) = igrid/(meshes(1, icase)*meshes(2, icase))
         do icell = 1, nk
            angle = 2.0_wp*pi*dot_product(twist+real(g, wp) &
               & /real(meshes(:, icase), wp), real(reps(:, icell), wp))
            phase = exp(cmplx(0.0_wp, angle, wp))
            dense_r(:, icell) = dense_r(:, icell) &
               & +conjg(phase)*mesh_k(:, ik)/real(nk, wp)
            dense_k(:, ik) = dense_k(:, ik)+phase*mesh_r(:, icell)
         end do
      end do
      call bvk_separable_k_to_r(mesh_k, sep_r, meshes(:, icase), &
         & input_to_grid, reps, phase_roots, twist_phase, .true.)
      call bvk_separable_r_to_k(mesh_r, sep_k, meshes(:, icase), &
         & input_to_grid, reps, phase_roots, twist_phase)
      err = max(err, maxval(abs(sep_r-dense_r)), &
         & maxval(abs(sep_k-dense_k)))

      ! Exact complex adjoint identity for the normalized inverse transform.
      lhs = sum(conjg(mesh_r)*sep_r)
      call bvk_separable_r_to_k(mesh_r, sep_k, meshes(:, icase), &
         & input_to_grid, reps, phase_roots, twist_phase)
      rhs = sum(conjg(sep_k/real(nk, wp))*mesh_k)
      err = max(err, abs(lhs-rhs))

      deallocate(input_to_grid, reps, phase_roots, twist_phase, dense_k, &
         & dense_r, mesh_k, mesh_r, sep_k, sep_r)
   end do

   call check(error, err, 0.0_wp, thr=1.0e-11_wp)
   if (allocated(error)) call test_failed(error, &
      & "Separable regular-grid DFT differs from the dense character oracle")
end subroutine test_regular_separable_dft

subroutine test_regular_mixed_radix_fft(error)
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: ncase = 12, nrow = 5
   integer, parameter :: meshes(3, ncase) = reshape([ &
      & 1, 1, 1, &
      & 2, 1, 1, &
      & 4, 1, 1, &
      & 8, 1, 1, &
      & 6, 1, 1, &
      & 10, 1, 1, &
      & 12, 1, 1, &
      & 5, 1, 1, &
      & 7, 1, 1, &
      & 11, 1, 1, &
      & 3, 4, 1, &
      & 2, 3, 4], [3, ncase])
   integer :: g(3), icase, icell, idim, igrid, ik, j, nk, root
   integer, allocatable :: input_to_grid(:), reps(:, :)
   real(wp), parameter :: twist(3) = [0.071_wp, -0.113_wp, 0.137_wp]
   real(wp) :: adjoint_err, angle, oracle_err, real_adjoint_err, roundtrip_err
   complex(wp) :: lhs, phase, rhs
   complex(wp), allocatable :: dense_k(:, :), dense_r(:, :), &
      & fft_k(:, :), fft_r(:, :), fft_r_unscaled(:, :), mesh_k(:, :), &
      & mesh_r(:, :), phase_roots(:, :), roundtrip_k(:, :), twist_phase(:)

   adjoint_err = 0.0_wp
   oracle_err = 0.0_wp
   real_adjoint_err = 0.0_wp
   roundtrip_err = 0.0_wp
   do icase = 1, ncase
      nk = product(meshes(:, icase))
      allocate(input_to_grid(nk), reps(3, nk))
      allocate(phase_roots(maxval(meshes(:, icase)), 3), twist_phase(nk), &
         & dense_k(nrow, nk), dense_r(nrow, nk), fft_k(nrow, nk), &
         & fft_r(nrow, nk), fft_r_unscaled(nrow, nk), mesh_k(nrow, nk), &
         & mesh_r(nrow, nk), roundtrip_k(nrow, nk), &
         & source=(0.0_wp, 0.0_wp))

      do ik = 1, nk
         input_to_grid(ik) = 1+modulo(ik+1, nk)
      end do
      do icell = 1, nk
         igrid = nk-icell
         reps(1, icell) = modulo(igrid, meshes(1, icase))
         reps(2, icell) = modulo(igrid/meshes(1, icase), &
            & meshes(2, icase))
         reps(3, icell) = igrid/(meshes(1, icase)*meshes(2, icase))
         if (modulo(icell, 2) == 0) &
            & reps(1, icell) = reps(1, icell)+meshes(1, icase)
         if (modulo(icell, 3) == 0) &
            & reps(2, icell) = reps(2, icell)-meshes(2, icase)
         if (modulo(icell, 5) == 0) &
            & reps(3, icell) = reps(3, icell)+meshes(3, icase)
         angle = 2.0_wp*pi*dot_product(twist, &
            & real(reps(:, icell), wp))
         twist_phase(icell) = exp(cmplx(0.0_wp, angle, wp))
      end do
      phase_roots = (1.0_wp, 0.0_wp)
      do idim = 1, 3
         do root = 0, meshes(idim, icase)-1
            angle = 2.0_wp*pi*real(root, wp) &
               & /real(meshes(idim, icase), wp)
            phase_roots(root+1, idim) = exp(cmplx(0.0_wp, angle, wp))
         end do
      end do
      do ik = 1, nk
         do j = 1, nrow
            mesh_k(j, ik) = cmplx( &
               & sin(real(7*j+11*ik+3*icase, wp)), &
               & cos(real(13*j+5*ik+2*icase, wp)), wp)
            mesh_r(j, ik) = cmplx( &
               & cos(real(3*j+17*ik+icase, wp)), &
               & sin(real(19*j+2*ik+4*icase, wp)), wp)
         end do
      end do

      dense_r = (0.0_wp, 0.0_wp)
      dense_k = (0.0_wp, 0.0_wp)
      do ik = 1, nk
         igrid = input_to_grid(ik)-1
         g(1) = modulo(igrid, meshes(1, icase))
         g(2) = modulo(igrid/meshes(1, icase), meshes(2, icase))
         g(3) = igrid/(meshes(1, icase)*meshes(2, icase))
         do icell = 1, nk
            angle = 2.0_wp*pi*dot_product(twist+real(g, wp) &
               & /real(meshes(:, icase), wp), real(reps(:, icell), wp))
            phase = exp(cmplx(0.0_wp, angle, wp))
            dense_r(:, icell) = dense_r(:, icell) &
               & +conjg(phase)*mesh_k(:, ik)/real(nk, wp)
            dense_k(:, ik) = dense_k(:, ik)+phase*mesh_r(:, icell)
         end do
      end do

      call bvk_mixed_radix_k_to_r(mesh_k, fft_r, meshes(:, icase), &
         & input_to_grid, reps, phase_roots, twist_phase, .true.)
      call bvk_mixed_radix_r_to_k(mesh_r, fft_k, meshes(:, icase), &
         & input_to_grid, reps, phase_roots, twist_phase)
      call bvk_mixed_radix_k_to_r(mesh_k, fft_r_unscaled, &
         & meshes(:, icase), input_to_grid, reps, phase_roots, twist_phase, &
         & .false.)
      oracle_err = max(oracle_err, maxval(abs(fft_r-dense_r)), &
         & maxval(abs(fft_k-dense_k)), &
         & maxval(abs(fft_r_unscaled-real(nk, wp)*dense_r)))

      call bvk_mixed_radix_r_to_k(fft_r, roundtrip_k, meshes(:, icase), &
         & input_to_grid, reps, phase_roots, twist_phase)
      roundtrip_err = max(roundtrip_err, maxval(abs(roundtrip_k-mesh_k)))

      ! The unnormalised k-to-R transform is the exact adjoint of R-to-k.
      lhs = sum(conjg(mesh_r)*fft_r_unscaled)
      rhs = sum(conjg(fft_k)*mesh_k)
      adjoint_err = max(adjoint_err, abs(lhs-rhs))
      real_adjoint_err = max(real_adjoint_err, &
         & abs(real(lhs, wp)-real(rhs, wp)))

      deallocate(input_to_grid, reps, phase_roots, twist_phase, dense_k, &
         & dense_r, fft_k, fft_r, fft_r_unscaled, mesh_k, mesh_r, roundtrip_k)
   end do

   print *, "MIXED_RADIX_FFT_TRANSFORM_MAX_RESID", oracle_err, &
      & roundtrip_err, adjoint_err, real_adjoint_err
   call check(error, max(oracle_err, roundtrip_err, adjoint_err, &
      & real_adjoint_err), 0.0_wp, thr=1.0e-11_wp)
   if (allocated(error)) call test_failed(error, &
      & "Mixed-radix FFT differs from the dense character oracle")
end subroutine test_regular_mixed_radix_fft

end module tblite_exchange_integration_tests
#endif

program tblite_exchange_integration_unittest
#if defined(__TBLITE_GXTB_ACCELERATION)
   use, intrinsic :: iso_fortran_env, only: output_unit
   use mctc_env_testing, only: run_testsuite
   use tblite_exchange_integration_tests, only: collect_exchange_integration
   implicit none
   integer :: failures
   failures = 0
   call run_testsuite(collect_exchange_integration, output_unit, failures)
   if (failures /= 0) error stop "CP2K exchange integration tests failed"
   print '(a)', 'CP2K exchange integration tests passed'
#else
   implicit none
   print '(a)', 'CP2K exchange integration tests skipped: g-xTB response API unavailable'
#endif
end program tblite_exchange_integration_unittest
