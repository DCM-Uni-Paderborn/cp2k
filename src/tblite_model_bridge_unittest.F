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

! CP2K integration tests relocated from save_tblite; original assertions retained.
#if defined(__TBLITE_GXTB_ACCELERATION)
MODULE tblite_model_bridge_tests
   USE tblite_projector_cache,          ONLY: cp2k_acp_projector_cache_type,&
                                              projector_cache_counts

!$ use omp_lib, only: omp_get_max_threads, omp_set_num_threads
   USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_quiet_nan, ieee_value
   USE mctc_env, ONLY: wp
   USE mctc_env_testing, ONLY: new_unittest, unittest_type, error_type, check, &
      & test_failed
   USE mctc_io, ONLY: structure_type, new
   USE mctc_io_constants, ONLY: pi
   USE mctc_io_convert, ONLY: aatoau
   USE mctc_io_math, ONLY: matinv_3x3
   USE tblite_acp_cache, ONLY: acp_cache
   USE tblite_acp_type, ONLY: get_acp, get_acp_gradient
   USE tblite_adjlist, ONLY: adjacency_list, new_adjacency_list
   USE tblite_basis_cache, ONLY: basis_cache
   USE tblite_basis_type, ONLY: basis_type, get_cutoff, get_pair_cutoff
   USE tblite_context_type, ONLY: context_type
   USE tblite_exchange_stream, ONLY: acp_image_type, cp2k_acp_image_gradient, &
      & cp2k_acp_projector_image_gradient, cp2k_acp_kmesh, &
      & cp2k_acp_kmesh_gradient, cp2k_acp_kmesh_gradient_stream_type, &
      & cp2k_acp_kmesh_gradient_stream_begin, &
      & cp2k_acp_kmesh_gradient_stream_push, &
      & cp2k_acp_kmesh_gradient_stream_end, &
      & cp2k_acp_kmesh_gradient_stream_discard, &
      & cp2k_acp_kmesh_gradient_stream_real_elements, &
      & cp2k_acp_kmesh_gradient_stream_complex_elements, &
      & cp2k_apply_basis_gradient, cp2k_build_gamma_h0, &
      & cp2k_build_acp_projectors, cp2k_build_image_acp, cp2k_build_image_h0, &
      & cp2k_get_acp_cutoff, &
      & cp2k_get_scc_mixer_dimension, cp2k_potential_mixer_type, &
      & cp2k_update_basis, hamiltonian_image_type, new_cp2k_potential_mixer
   USE tblite_cutoff, ONLY: get_lattice_points
   USE tblite_wavefunction, ONLY: wavefunction_type, new_wavefunction, eeqbc_guess
   USE tblite_xtb_calculator, ONLY: xtb_calculator, xtb_config
   USE tblite_xtb_gxtb, ONLY: new_gxtb_calculator
   USE tblite_xtb_h0, ONLY: get_selfenergy
   USE tblite_xtb_singlepoint, ONLY: xtb_singlepoint
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: collect_model_bridge

   REAL(wp), PARAMETER :: acc = 0.00001_wp
   REAL(wp), PARAMETER :: thr = 100*EPSILON(1.0_wp)
   REAL(wp), PARAMETER :: thr1 = 1e5*EPSILON(1.0_wp)
   REAL(wp), PARAMETER :: thr2 = SQRT(EPSILON(1.0_wp))
   REAL(wp), PARAMETER :: kt = 3.166808578545117e-06_wp

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param testsuite ...
! **************************************************************************************************
   SUBROUTINE collect_model_bridge(testsuite)
      TYPE(unittest_type), ALLOCATABLE, INTENT(out)      :: testsuite(:)

      testsuite = [ &
         new_unittest("image-h0-h2o-pbc", &
            & test_image_h0_h2o_pbc), &
         new_unittest("image-acp-gradient-h2o-pbc", &
            & test_image_acp_gradient_h2o_pbc), &
         new_unittest("kmesh-acp-h2o-pbc", &
            & test_kmesh_acp_h2o_pbc), &
         new_unittest("cp2k-acp-cutoff", &
            & test_cp2k_acp_cutoff), &
         new_unittest("cp2k-scc-mixer-dimension", &
            & test_cp2k_scc_mixer_dimension), &
         new_unittest("cp2k-potential-mixer", &
            & test_cp2k_potential_mixer), &
         new_unittest("cp2k-potential-mixer-uninitialized", &
            & test_cp2k_potential_mixer_uninitialized, should_fail=.TRUE.), &
         new_unittest("cp2k-potential-mixer-invalid-dimension", &
            & test_cp2k_potential_mixer_invalid_dimension, should_fail=.TRUE.), &
         new_unittest("cp2k-potential-mixer-size-mismatch", &
            & test_cp2k_potential_mixer_size_mismatch, should_fail=.TRUE.), &
         new_unittest("cp2k-potential-mixer-nonfinite", &
            & test_cp2k_potential_mixer_nonfinite, should_fail=.TRUE.)]
   END SUBROUTINE collect_model_bridge

! **************************************************************************************************
!> \brief ...
!> \param error ...
! **************************************************************************************************
   SUBROUTINE test_image_h0_h2o_pbc(error)
      !> Error handling
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      INTEGER, PARAMETER                                 :: num(3) = [1, 1, 8]
      LOGICAL, PARAMETER :: periodic(3) = [.TRUE., .TRUE., .TRUE.]
      REAL(wp), PARAMETER :: lattice(3, 3) = aatoau*RESHAPE([7.800000_wp, 0.200000_wp, 0.100000_wp &
         , 0.400000_wp, 8.300000_wp, 0.300000_wp, 0.100000_wp, 0.500000_wp, 9.100000_wp], [3, 3]), &
         qvec(3) = [0.17_wp, 0.13_wp, 0.11_wp], xyz(3, 3) = aatoau*RESHAPE([0.200000_wp, &
         0.300000_wp, 0.400000_wp, 0.958602_wp, 0.300000_wp, 0.904284_wp, 0.200000_wp, 0.300000_wp,&
         1.400000_wp], [3, 3])

      COMPLEX(wp)                                        :: phase
      COMPLEX(wp), ALLOCATABLE                           :: ck(:, :), csk(:, :), hk(:, :), &
                                                            hk_acp(:, :), hk_acp_direct(:, :), &
                                                            sk(:, :)
      INTEGER                                            :: itr, nao, saved_threads
      REAL(wp), ALLOCATABLE :: acp_gamma(:, :), cn(:), dipole(:, :, :), hamiltonian(:, :), &
         hamiltonian_total(:, :), lattr(:, :), lattr_acp(:, :), lattr_acp_ref(:, :), &
         overlap(:, :), quadrupole(:, :, :), selfenergy(:)
      TYPE(acp_cache)                                    :: acache
      TYPE(acp_image_type)                               :: acp_images, acp_images_ref
      TYPE(basis_cache)                                  :: bcache
      TYPE(hamiltonian_image_type)                       :: images, images_serial
      TYPE(structure_type)                               :: mol
      TYPE(xtb_calculator)                               :: calc

      CALL new(mol, num, xyz, lattice=lattice, periodic=periodic)
      CALL new_gxtb_calculator(calc, mol, error, accuracy=acc)
      IF (ALLOCATED(error)) RETURN
      CALL cp2k_update_basis(calc, mol, bcache, error)
      IF (ALLOCATED(error)) RETURN

      ALLOCATE (cn(mol%nat), selfenergy(calc%bas%nsh))
      CALL calc%ncoord%get_cn(mol, cn)
      CALL get_selfenergy(calc%h0, mol%id, calc%bas%ish_at, calc%bas%nsh_id, &
         & cn=cn, selfenergy=selfenergy)
      CALL get_lattice_points(mol%periodic, mol%lattice, &
         & get_cutoff(calc%bas, acc), lattr)
      CALL get_lattice_points(mol%periodic, mol%lattice, &
         & cp2k_get_acp_cutoff(calc, acc), lattr_acp_ref)

      nao = calc%bas%nao
      ALLOCATE (overlap(nao, nao), dipole(3, nao, nao), &
         & quadrupole(6, nao, nao), hamiltonian(nao, nao))
      saved_threads = 1
!$    saved_threads = omp_get_max_threads()
!$    call omp_set_num_threads(1)
      CALL cp2k_build_image_h0(calc, mol, bcache, acc, selfenergy, lattr, images_serial)
!$    call omp_set_num_threads(saved_threads)
      CALL cp2k_build_image_h0(calc, mol, bcache, acc, selfenergy, lattr, images, &
         & overlap, dipole, quadrupole, hamiltonian)

      CALL check(error, MAXVAL(ABS(images%hamiltonian - images_serial%hamiltonian)), 0.0_wp, thr=1.0e-12_wp)
      IF (ALLOCATED(error)) RETURN
      CALL check(error, MAXVAL(ABS(images%overlap - images_serial%overlap)), 0.0_wp, thr=1.0e-12_wp)
      IF (ALLOCATED(error)) RETURN
      CALL check(error, MAXVAL(ABS(images%dipole_bra - images_serial%dipole_bra)), 0.0_wp, thr=1.0e-12_wp)
      IF (ALLOCATED(error)) RETURN
      CALL check(error, MAXVAL(ABS(images%dipole_ket - images_serial%dipole_ket)), 0.0_wp, thr=1.0e-12_wp)
      IF (ALLOCATED(error)) RETURN
      CALL check(error, MAXVAL(ABS(images%quadrupole_bra - images_serial%quadrupole_bra)), 0.0_wp, thr=1.0e-12_wp)
      IF (ALLOCATED(error)) RETURN
      CALL check(error, MAXVAL(ABS(images%quadrupole_ket - images_serial%quadrupole_ket)), 0.0_wp, thr=1.0e-12_wp)
      IF (ALLOCATED(error)) RETURN

      CALL check(error, MAXVAL(ABS(SUM(images%overlap, dim=3) - overlap)), &
         & 0.0_wp, thr=1.0e-11_wp)
      IF (ALLOCATED(error)) THEN
         CALL test_failed(error, "Image overlap does not reproduce Gamma overlap")
         RETURN
      END IF
      CALL check(error, MAXVAL(ABS(SUM(images%dipole_bra, dim=4) - dipole)), &
         & 0.0_wp, thr=1.0e-11_wp)
      IF (ALLOCATED(error)) THEN
         CALL test_failed(error, "Image dipoles do not reproduce Gamma dipoles")
         RETURN
      END IF
      CALL check(error, MAXVAL(ABS(SUM(images%quadrupole_bra, dim=4) - quadrupole)), &
         & 0.0_wp, thr=1.0e-11_wp)
      IF (ALLOCATED(error)) THEN
         CALL test_failed(error, "Image quadrupoles do not reproduce Gamma quadrupoles")
         RETURN
      END IF
      CALL check(error, MAXVAL(ABS(SUM(images%hamiltonian, dim=3) - hamiltonian)), &
         & 0.0_wp, thr=1.0e-11_wp)
      IF (ALLOCATED(error)) THEN
         CALL test_failed(error, "Image H0 does not reproduce Gamma H0")
         RETURN
      END IF

      ALLOCATE (sk(nao, nao), hk(nao, nao), source=(0.0_wp, 0.0_wp))
      DO itr = 1, SIZE(images%translation, 2)
         phase = EXP(CMPLX(0.0_wp, DOT_PRODUCT(qvec, images%translation(:, itr)), wp))
         sk = sk + phase*images%overlap(:, :, itr)
         hk = hk + phase*images%hamiltonian(:, :, itr)
      END DO
      CALL check(error, MAXVAL(ABS(sk - CONJG(TRANSPOSE(sk)))), &
         & 0.0_wp, thr=1.0e-11_wp)
      IF (ALLOCATED(error)) THEN
         CALL test_failed(error, "Bloch-summed image overlap is not Hermitian")
         RETURN
      END IF
      CALL check(error, MAXVAL(ABS(hk - CONJG(TRANSPOSE(hk)))), &
         & 0.0_wp, thr=1.0e-11_wp)
      IF (ALLOCATED(error)) THEN
         CALL test_failed(error, "Bloch-summed image H0 is not Hermitian")
         RETURN
      END IF

      ! ACP is a separable projector convolution, so its Hamiltonian image range
      ! is the difference set of the projector-overlap translations.
      ALLOCATE (hamiltonian_total(nao, nao), acp_gamma(nao, nao))
      CALL cp2k_build_gamma_h0(calc, mol, bcache, acache, acc, selfenergy, &
         & overlap, dipole, quadrupole, hamiltonian_total)
      acp_gamma = hamiltonian_total - hamiltonian
      ALLOCATE (lattr_acp(3, SIZE(lattr_acp_ref, 2) + 2))
      lattr_acp(:, :SIZE(lattr_acp_ref, 2)) = lattr_acp_ref
      lattr_acp(:, SIZE(lattr_acp_ref, 2) + 1) = 20.0_wp*SUM(lattice, dim=2)
      lattr_acp(:, SIZE(lattr_acp_ref, 2) + 2) = -lattr_acp(:, SIZE(lattr_acp_ref, 2) + 1)
      CALL cp2k_build_image_acp(calc, mol, bcache, acache, lattr_acp_ref, acp_images_ref)
      CALL cp2k_build_image_acp(calc, mol, bcache, acache, lattr_acp, acp_images)
      CALL check(error, SIZE(acp_images%projector_translation, 2) <= &
         & SIZE(lattr_acp, 2) - 2)
      IF (ALLOCATED(error)) THEN
         CALL test_failed(error, "Zero ACP projector translations were not pruned")
         RETURN
      END IF
      DO itr = 1, SIZE(acp_images%projector_overlap, 3)
         CALL check(error, ANY(acp_images%projector_overlap(:, :, itr) /= 0.0_wp))
         IF (ALLOCATED(error)) THEN
            CALL test_failed(error, "Compact ACP image set contains an empty image")
            RETURN
         END IF
      END DO
      CALL check(error, SIZE(acp_images%projector_translation, 2), &
         & SIZE(acp_images_ref%projector_translation, 2))
      IF (.NOT. ALLOCATED(error)) THEN
         CALL check(error, SIZE(acp_images%translation, 2), &
            & SIZE(acp_images_ref%translation, 2))
      END IF
      IF (ALLOCATED(error)) THEN
         CALL test_failed(error, "Empty outer translations changed compact ACP support")
         RETURN
      END IF
      CALL check(error, MAX( &
         & MAXVAL(ABS(acp_images%projector_translation &
         & - acp_images_ref%projector_translation)), &
         & MAXVAL(ABS(acp_images%projector_overlap &
         & - acp_images_ref%projector_overlap)), &
         & MAXVAL(ABS(acp_images%translation - acp_images_ref%translation)), &
         & MAXVAL(ABS(acp_images%hamiltonian - acp_images_ref%hamiltonian))), &
         & 0.0_wp, thr=1.0e-14_wp)
      IF (ALLOCATED(error)) THEN
         CALL test_failed(error, "Expanded zero input changed compact ACP images")
         RETURN
      END IF
      CALL check(error, &
         & MAXVAL(ABS(SUM(acp_images%hamiltonian, dim=3) - acp_gamma)), &
         & 0.0_wp, thr=1.0e-11_wp)
      IF (ALLOCATED(error)) THEN
         CALL test_failed(error, "Image ACP convolution does not reproduce Gamma ACP")
         RETURN
      END IF

      ALLOCATE (ck(SIZE(acp_images%projector_overlap, 1), nao), &
         & csk(SIZE(acp_images%projector_overlap, 1), nao), &
         & hk_acp(nao, nao), hk_acp_direct(nao, nao), &
         & source=(0.0_wp, 0.0_wp))
      DO itr = 1, SIZE(acp_images%projector_translation, 2)
         phase = EXP(CMPLX(0.0_wp, -DOT_PRODUCT(qvec, &
            & acp_images%projector_translation(:, itr)), wp))
         ck = ck + phase*acp_images%projector_overlap(:, :, itr)
         csk = csk + phase*acp_images%scaled_projector_overlap(:, :, itr)
      END DO
      DO itr = 1, SIZE(acp_images%translation, 2)
         phase = EXP(CMPLX(0.0_wp, &
            & DOT_PRODUCT(qvec, acp_images%translation(:, itr)), wp))
         hk_acp = hk_acp + phase*acp_images%hamiltonian(:, :, itr)
      END DO
      hk_acp_direct = MATMUL(CONJG(TRANSPOSE(ck)), csk)
      CALL check(error, MAXVAL(ABS(hk_acp - hk_acp_direct)), &
         & 0.0_wp, thr=1.0e-11_wp)
      IF (ALLOCATED(error)) THEN
         CALL test_failed(error, "Image ACP does not reproduce direct Bloch product")
         RETURN
      END IF
      CALL check(error, MAXVAL(ABS(hk_acp - CONJG(TRANSPOSE(hk_acp)))), &
         & 0.0_wp, thr=1.0e-11_wp)
      IF (ALLOCATED(error)) THEN
         CALL test_failed(error, "Bloch-summed image ACP is not Hermitian")
         RETURN
      END IF

   END SUBROUTINE test_image_h0_h2o_pbc

! **************************************************************************************************
!> \brief ...
!> \param error ...
! **************************************************************************************************
   SUBROUTINE test_image_acp_gradient_h2o_pbc(error)
      !> Error handling
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      INTEGER, PARAMETER                                 :: num(3) = [1, 1, 8]
      LOGICAL, PARAMETER :: periodic(3) = [.TRUE., .TRUE., .TRUE.]
      REAL(wp), PARAMETER :: fd_thr = 2.0e-7_wp, lattice(3, 3) = aatoau*RESHAPE([7.800000_wp, &
         0.200000_wp, 0.100000_wp, 0.400000_wp, 8.300000_wp, 0.300000_wp, 0.100000_wp, 0.500000_wp,&
         9.100000_wp], [3, 3]), qvec(3) = [0.17_wp, 0.13_wp, 0.11_wp], step = 1.0e-5_wp, xyz(3, 3) &
         = aatoau*RESHAPE([0.200000_wp, 0.300000_wp, 0.400000_wp, 0.958602_wp, 0.300000_wp, &
         0.904284_wp, 0.200000_wp, 0.300000_wp, 1.400000_wp], [3, 3])

      INTEGER                                            :: iao, iat, ic, ish, isp, itr, jao, jc, &
                                                            nao, nimg
      REAL(wp)                                           :: arg, el, eps(3, 3), er, sigma(3, 3), &
                                                            sigma_full(3, 3), sigma_gamma(3, 3), &
                                                            sigma_image(3, 3), unity(3, 3)
      REAL(wp), ALLOCATABLE :: amat(:, :), bmat(:, :), dEdcn(:), dEdcn_gamma(:), dEdcn_image(:), &
         dEdq(:), dEdq_gamma(:), dEdq_image(:), density(:, :, :, :), density_gamma(:, :, :, :), &
         density_translation(:, :), density_translation_work(:, :), gradient(:, :), &
         gradient_full(:, :), gradient_gamma(:, :), gradient_image(:, :), hamiltonian_gamma(:, :), &
         lattr(:, :), lattr_work(:, :), num_dEdcn(:), num_dEdq(:), numgrad(:, :), &
         numgrad_full(:, :), numsigma(:, :), numsigma_full(:, :)
      TYPE(acp_cache)                                    :: acache, acache_gamma, acache_work
      TYPE(acp_image_type)                               :: images, images_work
      TYPE(adjacency_list)                               :: list
      TYPE(basis_cache)                                  :: bcache, bcache_work
      TYPE(structure_type)                               :: mol, mol_work
      TYPE(wavefunction_type)                            :: wfn, wfn_aux, wfn_aux_work
      TYPE(xtb_calculator)                               :: calc

      CALL new(mol, num, xyz, lattice=lattice, periodic=periodic)
      CALL new_gxtb_calculator(calc, mol, error, accuracy=acc)
      IF (ALLOCATED(error)) RETURN
      CALL cp2k_update_basis(calc, mol, bcache, error, grad=.TRUE., &
         & wfn_aux=wfn_aux)
      IF (ALLOCATED(error)) RETURN
      CALL get_lattice_points(mol%periodic, mol%lattice, &
         & cp2k_get_acp_cutoff(calc, acc), lattr)
      CALL cp2k_build_image_acp(calc, mol, bcache, acache, lattr, images)

      nao = calc%bas%nao
      nimg = SIZE(images%translation, 2)
      ALLOCATE (amat(nao, nao), bmat(nao, nao), &
         & density(nao, nao, nimg, 1), density_gamma(nao, nao, nimg, 1), &
         & density_translation(3, nimg))
      density_translation = images%translation

      ! One artificial complex-Hermitian k-point density, transformed to real
      ! images.  A is symmetric and B antisymmetric, hence
      ! P(-R)=P(R)^T exactly (up to the trigonometric roundoff).
      DO iao = 1, nao
         DO jao = 1, nao
            amat(iao, jao) = 0.01_wp*(SIN(REAL(iao + 3*jao, wp)) &
               & + SIN(REAL(jao + 3*iao, wp)))
            bmat(iao, jao) = 0.01_wp*(COS(REAL(iao + 2*jao, wp)) &
               & - COS(REAL(jao + 2*iao, wp)))
         END DO
      END DO
      DO itr = 1, nimg
         arg = DOT_PRODUCT(qvec, density_translation(:, itr))
         density(:, :, itr, 1) = COS(arg)*amat + SIN(arg)*bmat
         density_gamma(:, :, itr, 1) = amat
      END DO

      ALLOCATE (dEdcn(mol%nat), dEdq(mol%nat), gradient(3, mol%nat), &
         & numgrad(3, mol%nat), numsigma(3, 3), source=0.0_wp)
      sigma = 0.0_wp
      CALL cp2k_acp_image_gradient(calc, mol, bcache, acache, images, &
         & density_translation, density, dEdcn, dEdq, gradient, sigma, error)
      IF (ALLOCATED(error)) RETURN

      ! Independently perturb the effective q-vSZP basis coordinate along its
      ! CN and charge derivatives.  Re-normalizing every affected CGTO makes
      ! these finite differences include both coefficient and normalization
      ! response for the nonconstant image density.
      ALLOCATE (num_dEdcn(mol%nat), num_dEdq(mol%nat), source=0.0_wp)
      DO iat = 1, mol%nat
         isp = mol%id(iat)
         DO ish = 1, calc%bas%nsh_id(isp)
            bcache%cgto(ish, iat)%qeff = bcache%cgto(ish, iat)%qeff &
               & + step*bcache%cgto(ish, iat)%dqeffdcn
            CALL calc%bas%cgto(ish, isp)%raw%get_normalization( &
               & bcache%cgto(ish, iat), .FALSE.)
         END DO
         CALL cp2k_build_image_acp(calc, mol, bcache, acache, lattr, images_work)
         er = contract_acp_image_energy(images_work, density_translation, density)
         DO ish = 1, calc%bas%nsh_id(isp)
            bcache%cgto(ish, iat)%qeff = bcache%cgto(ish, iat)%qeff &
               & - 2.0_wp*step*bcache%cgto(ish, iat)%dqeffdcn
            CALL calc%bas%cgto(ish, isp)%raw%get_normalization( &
               & bcache%cgto(ish, iat), .FALSE.)
         END DO
         CALL cp2k_build_image_acp(calc, mol, bcache, acache, lattr, images_work)
         el = contract_acp_image_energy(images_work, density_translation, density)
         num_dEdcn(iat) = 0.5_wp*(er - el)/step
         DO ish = 1, calc%bas%nsh_id(isp)
            bcache%cgto(ish, iat)%qeff = bcache%cgto(ish, iat)%qeff &
               & + step*bcache%cgto(ish, iat)%dqeffdcn
            CALL calc%bas%cgto(ish, isp)%raw%get_normalization( &
               & bcache%cgto(ish, iat), .FALSE.)
         END DO
      END DO
      CALL check(error, MAXVAL(ABS(dEdcn - num_dEdcn)), 0.0_wp, thr=fd_thr)
      IF (ALLOCATED(error)) THEN
         CALL test_failed(error, "Image ACP dE/dCN does not match finite differences")
         RETURN
      END IF

      DO iat = 1, mol%nat
         isp = mol%id(iat)
         DO ish = 1, calc%bas%nsh_id(isp)
            bcache%cgto(ish, iat)%qeff = bcache%cgto(ish, iat)%qeff &
               & + step*bcache%cgto(ish, iat)%dqeffdq
            CALL calc%bas%cgto(ish, isp)%raw%get_normalization( &
               & bcache%cgto(ish, iat), .FALSE.)
         END DO
         CALL cp2k_build_image_acp(calc, mol, bcache, acache, lattr, images_work)
         er = contract_acp_image_energy(images_work, density_translation, density)
         DO ish = 1, calc%bas%nsh_id(isp)
            bcache%cgto(ish, iat)%qeff = bcache%cgto(ish, iat)%qeff &
               & - 2.0_wp*step*bcache%cgto(ish, iat)%dqeffdq
            CALL calc%bas%cgto(ish, isp)%raw%get_normalization( &
               & bcache%cgto(ish, iat), .FALSE.)
         END DO
         CALL cp2k_build_image_acp(calc, mol, bcache, acache, lattr, images_work)
         el = contract_acp_image_energy(images_work, density_translation, density)
         num_dEdq(iat) = 0.5_wp*(er - el)/step
         DO ish = 1, calc%bas%nsh_id(isp)
            bcache%cgto(ish, iat)%qeff = bcache%cgto(ish, iat)%qeff &
               & + step*bcache%cgto(ish, iat)%dqeffdq
            CALL calc%bas%cgto(ish, isp)%raw%get_normalization( &
               & bcache%cgto(ish, iat), .FALSE.)
         END DO
      END DO
      CALL check(error, MAXVAL(ABS(dEdq - num_dEdq)), 0.0_wp, thr=fd_thr)
      IF (ALLOCATED(error)) THEN
         CALL test_failed(error, "Image ACP dE/dq does not match finite differences")
         RETURN
      END IF

      ! Restore derivative-bearing normalization caches for the remaining
      ! comparisons (the finite differences above only need their values).
      CALL calc%bas%update(mol, bcache, .TRUE., wfn_aux)
      ALLOCATE (gradient_full(3, mol%nat), source=gradient)
      sigma_full = sigma
      CALL cp2k_apply_basis_gradient(calc, mol, wfn_aux, dEdcn, dEdq, &
         & gradient_full, sigma_full)

      ! Frozen-density Cartesian finite differences of the exact image energy.
      DO iat = 1, mol%nat
         DO ic = 1, 3
            mol_work = mol
            mol_work%xyz(ic, iat) = mol_work%xyz(ic, iat) + step
            CALL cp2k_build_image_acp(calc, mol_work, bcache, acache, lattr, &
               & images_work)
            er = contract_acp_image_energy(images_work, density_translation, density)

            mol_work = mol
            mol_work%xyz(ic, iat) = mol_work%xyz(ic, iat) - step
            CALL cp2k_build_image_acp(calc, mol_work, bcache, acache, lattr, &
               & images_work)
            el = contract_acp_image_energy(images_work, density_translation, density)
            numgrad(ic, iat) = 0.5_wp*(er - el)/step
         END DO
      END DO
      CALL check(error, MAXVAL(ABS(gradient - numgrad)), 0.0_wp, thr=fd_thr)
      IF (ALLOCATED(error)) THEN
         CALL test_failed(error, "Image ACP Cartesian gradient does not match finite differences")
         RETURN
      END IF

      ! Repeat the same displacement finite difference with the q-vSZP basis
      ! and its external charge model rebuilt at every geometry.  This checks
      ! the full ACP -> basis CN/q -> EEQ Cartesian chain rule.
      ALLOCATE (numgrad_full(3, mol%nat), source=0.0_wp)
      DO iat = 1, mol%nat
         DO ic = 1, 3
            mol_work = mol
            mol_work%xyz(ic, iat) = mol_work%xyz(ic, iat) + step
            CALL cp2k_update_basis(calc, mol_work, bcache_work, error, &
               & wfn_aux=wfn_aux_work)
            IF (ALLOCATED(error)) RETURN
            CALL cp2k_build_image_acp(calc, mol_work, bcache_work, acache_work, &
               & lattr, images_work)
            er = contract_acp_image_energy(images_work, density_translation, density)

            mol_work = mol
            mol_work%xyz(ic, iat) = mol_work%xyz(ic, iat) - step
            CALL cp2k_update_basis(calc, mol_work, bcache_work, error, &
               & wfn_aux=wfn_aux_work)
            IF (ALLOCATED(error)) RETURN
            CALL cp2k_build_image_acp(calc, mol_work, bcache_work, acache_work, &
               & lattr, images_work)
            el = contract_acp_image_energy(images_work, density_translation, density)
            numgrad_full(ic, iat) = 0.5_wp*(er - el)/step
         END DO
      END DO
      CALL check(error, MAXVAL(ABS(gradient_full - numgrad_full)), 0.0_wp, &
         & thr=fd_thr)
      IF (ALLOCATED(error)) THEN
         CALL test_failed(error, &
            & "Full image ACP q-vSZP Cartesian response does not match finite differences")
         RETURN
      END IF

      ! Frozen-density strain finite differences.  The image labels follow the
      ! strained Cartesian lattice while their matrix values stay fixed.
      unity = 0.0_wp
      DO ic = 1, 3
         unity(ic, ic) = 1.0_wp
      END DO
      ALLOCATE (lattr_work(3, SIZE(lattr, 2)), &
         & density_translation_work(3, nimg))
      DO ic = 1, 3
         DO jc = 1, 3
            eps = unity
            eps(jc, ic) = eps(jc, ic) + step
            mol_work = mol
            mol_work%xyz = MATMUL(eps, mol%xyz)
            mol_work%lattice = MATMUL(eps, mol%lattice)
            lattr_work = MATMUL(eps, lattr)
            density_translation_work = MATMUL(eps, density_translation)
            CALL cp2k_build_image_acp(calc, mol_work, bcache, acache, lattr_work, &
               & images_work)
            er = contract_acp_image_energy(images_work, &
               & density_translation_work, density)

            eps = unity
            eps(jc, ic) = eps(jc, ic) - step
            mol_work = mol
            mol_work%xyz = MATMUL(eps, mol%xyz)
            mol_work%lattice = MATMUL(eps, mol%lattice)
            lattr_work = MATMUL(eps, lattr)
            density_translation_work = MATMUL(eps, density_translation)
            CALL cp2k_build_image_acp(calc, mol_work, bcache, acache, lattr_work, &
               & images_work)
            el = contract_acp_image_energy(images_work, &
               & density_translation_work, density)
            numsigma(jc, ic) = 0.5_wp*(er - el)/step
         END DO
      END DO
      numsigma = 0.5_wp*(numsigma + TRANSPOSE(numsigma))
      CALL check(error, MAXVAL(ABS(sigma - numsigma)), 0.0_wp, thr=fd_thr)
      IF (ALLOCATED(error)) THEN
         CALL test_failed(error, "Image ACP strain derivative does not match finite differences")
         RETURN
      END IF

      ! Full q-vSZP/EEQ strain response, with all Cartesian lattice labels and
      ! basis caches transformed/rebuilt at the displaced cell.
      ALLOCATE (numsigma_full(3, 3), source=0.0_wp)
      DO ic = 1, 3
         DO jc = 1, 3
            eps = unity
            eps(jc, ic) = eps(jc, ic) + step
            mol_work = mol
            mol_work%xyz = MATMUL(eps, mol%xyz)
            mol_work%lattice = MATMUL(eps, mol%lattice)
            lattr_work = MATMUL(eps, lattr)
            density_translation_work = MATMUL(eps, density_translation)
            CALL cp2k_update_basis(calc, mol_work, bcache_work, error, &
               & wfn_aux=wfn_aux_work)
            IF (ALLOCATED(error)) RETURN
            CALL cp2k_build_image_acp(calc, mol_work, bcache_work, acache_work, &
               & lattr_work, images_work)
            er = contract_acp_image_energy(images_work, &
               & density_translation_work, density)

            eps = unity
            eps(jc, ic) = eps(jc, ic) - step
            mol_work = mol
            mol_work%xyz = MATMUL(eps, mol%xyz)
            mol_work%lattice = MATMUL(eps, mol%lattice)
            lattr_work = MATMUL(eps, lattr)
            density_translation_work = MATMUL(eps, density_translation)
            CALL cp2k_update_basis(calc, mol_work, bcache_work, error, &
               & wfn_aux=wfn_aux_work)
            IF (ALLOCATED(error)) RETURN
            CALL cp2k_build_image_acp(calc, mol_work, bcache_work, acache_work, &
               & lattr_work, images_work)
            el = contract_acp_image_energy(images_work, &
               & density_translation_work, density)
            numsigma_full(jc, ic) = 0.5_wp*(er - el)/step
         END DO
      END DO
      numsigma_full = 0.5_wp*(numsigma_full + TRANSPOSE(numsigma_full))
      CALL check(error, MAXVAL(ABS(sigma_full - numsigma_full)), 0.0_wp, &
         & thr=fd_thr)
      IF (ALLOCATED(error)) THEN
         CALL test_failed(error, &
            & "Full image ACP q-vSZP strain response does not match finite differences")
         RETURN
      END IF

      ! A constant symmetric density on every image is precisely the Gamma
      ! contraction.  This also checks the q-vSZP dE/dCN and dE/dq channels
      ! against the established Gamma implementation.
      CALL cp2k_build_image_acp(calc, mol, bcache, acache, lattr, images)
      CALL new_adjacency_list(list, mol, lattr, cp2k_get_acp_cutoff(calc, acc), complete=.TRUE.)
      ! Build the independent Gamma oracle through the library's public path.
      ! Its gradient needs the compact AO maps as well as projector overlaps;
      ! the image cache deliberately does not construct those Gamma-only data.
      CALL calc%acp%update(mol, acache_gamma)
      ALLOCATE (hamiltonian_gamma(nao, nao), source=0.0_wp)
      CALL get_acp(mol, lattr, list, calc%bas, bcache, calc%acp, &
         & acache_gamma, hamiltonian_gamma)
      CALL new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, &
         & 0.0_wp, .FALSE.)
      wfn%density(:, :, 1) = amat
      ALLOCATE (dEdcn_gamma(mol%nat), dEdq_gamma(mol%nat), &
         & gradient_gamma(3, mol%nat), dEdcn_image(mol%nat), &
         & dEdq_image(mol%nat), gradient_image(3, mol%nat), source=0.0_wp)
      sigma_gamma = 0.0_wp
      sigma_image = 0.0_wp
      CALL get_acp_gradient(mol, lattr, list, calc%bas, bcache, calc%acp, &
         & acache_gamma, wfn, dEdcn_gamma, dEdq_gamma, gradient_gamma, sigma_gamma)
      CALL cp2k_acp_image_gradient(calc, mol, bcache, acache, images, &
         & density_translation, density_gamma, dEdcn_image, dEdq_image, &
         & gradient_image, sigma_image, error)
      IF (ALLOCATED(error)) RETURN
      CALL check(error, MAX(MAXVAL(ABS(gradient_image - gradient_gamma)), &
         & MAXVAL(ABS(sigma_image - sigma_gamma)), &
         & MAXVAL(ABS(dEdcn_image - dEdcn_gamma)), &
         & MAXVAL(ABS(dEdq_image - dEdq_gamma))), 0.0_wp, thr=1.0e-11_wp)
      IF (ALLOCATED(error)) THEN
         CALL test_failed(error, "Constant ACP density images do not reproduce Gamma response")
      END IF
   END SUBROUTINE test_image_acp_gradient_h2o_pbc

! **************************************************************************************************
!> \brief ...
!> \param error ...
! **************************************************************************************************
   SUBROUTINE test_kmesh_acp_h2o_pbc(error)
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      INTEGER, PARAMETER                                 :: num(3) = [1, 1, 8]
      LOGICAL, PARAMETER :: periodic(3) = [.TRUE., .TRUE., .TRUE.]
      REAL(wp), PARAMETER :: &
         kfrac(3, 2) = RESHAPE([0.17_wp, 0.13_wp, 0.11_wp, -0.17_wp, -0.13_wp, -0.11_wp], [3, 2]), &
         lattice(3, 3) = aatoau*RESHAPE([7.800000_wp, 0.200000_wp, 0.100000_wp, 0.400000_wp, &
         8.300000_wp, 0.300000_wp, 0.100000_wp, 0.500000_wp, 9.100000_wp], [3, 3]), &
         weights(2) = [0.5_wp, 0.5_wp], xyz(3, 3) = aatoau*RESHAPE([0.200000_wp, 0.300000_wp, &
         0.400000_wp, 0.958602_wp, 0.300000_wp, 0.904284_wp, 0.200000_wp, 0.300000_wp, 1.400000_wp]&
         , [3, 3])

      COMPLEX(wp)                                        :: phase
      COMPLEX(wp), ALLOCATABLE                           :: density_k(:, :, :, :), &
                                                            h_compact(:, :, :), h_image(:, :, :), &
                                                            h_k(:, :, :)
      INTEGER                                            :: cache_bypasses, cache_hits, &
                                                            cache_misses, iao, ik, itr, jao, nao, &
                                                            nimg, origin
      LOGICAL                                            :: cache_valid
      REAL(wp) :: angle, energy_cached, energy_image, energy_k, saved_alpha, saved_level, &
         sigma_cached(3, 3), sigma_image(3, 3), sigma_k(3, 3), sigma_sparse(3, 3), &
         sigma_stream(3, 3)
      REAL(wp), ALLOCATABLE :: amat(:, :), bmat(:, :), cells(:, :), dEdcn_cached(:), &
         dEdcn_image(:), dEdcn_k(:), dEdcn_sparse(:), dEdcn_stream(:), dEdq_cached(:), &
         dEdq_image(:), dEdq_k(:), dEdq_sparse(:), dEdq_stream(:), density_cp2k(:, :, :, :), &
         density_image(:, :, :, :), density_sparse(:, :, :, :), density_zero(:, :, :, :), &
         gradient_cached(:, :), gradient_image(:, :), gradient_k(:, :), gradient_sparse(:, :), &
         gradient_streamed(:, :), lattr(:, :), lattr_changed(:, :), translation_sparse(:, :)
      TYPE(acp_cache)                                    :: acache
      TYPE(acp_image_type)                               :: images, projectors, projectors_changed
      TYPE(basis_cache)                                  :: bcache, bcache_changed
      TYPE(cp2k_acp_kmesh_gradient_stream_type)          :: gradient_stream
      TYPE(cp2k_acp_projector_cache_type)                :: projector_cache
      TYPE(structure_type)                               :: mol, mol_changed
      TYPE(wavefunction_type)                            :: wfn_aux, wfn_aux_changed
      TYPE(xtb_calculator)                               :: calc

      CALL new(mol, num, xyz, lattice=lattice, periodic=periodic)
      CALL new_gxtb_calculator(calc, mol, error, accuracy=acc)
      IF (ALLOCATED(error)) RETURN
      CALL cp2k_update_basis(calc, mol, bcache, error, grad=.TRUE., &
         & wfn_aux=wfn_aux)
      IF (ALLOCATED(error)) RETURN
      CALL get_lattice_points(mol%periodic, mol%lattice, &
         & cp2k_get_acp_cutoff(calc, acc), lattr)
      CALL cp2k_build_image_acp(calc, mol, bcache, acache, lattr, images, projector_cache=projector_cache)
      CALL projector_cache_counts(projector_cache, cache_hits, cache_misses, cache_bypasses, cache_valid)
      CALL check(error, cache_misses, 1)
      IF (.NOT. ALLOCATED(error)) CALL check(error, cache_hits, 0)
      IF (ALLOCATED(error)) THEN
         CALL test_failed(error, "Initial ACP projector build did not miss the cache")
         RETURN
      END IF

      nao = calc%bas%nao
      nimg = SIZE(images%translation, 2)
      ALLOCATE (amat(nao, nao), bmat(nao, nao), cells(3, nimg), &
         & density_cp2k(nao, nao, nimg, 1), &
         & density_image(nao, nao, nimg, 1), source=0.0_wp)
      ALLOCATE (density_k(nao, nao, 1, 2), h_compact(nao, nao, 2), &
         & h_image(nao, nao, 2), h_k(nao, nao, 2), &
         & source=CMPLX(0.0_wp, 0.0_wp, wp))
      DO iao = 1, nao
         DO jao = 1, nao
            amat(iao, jao) = 0.01_wp*(SIN(REAL(iao + 3*jao, wp)) &
               & + SIN(REAL(jao + 3*iao, wp)))
            bmat(iao, jao) = 0.01_wp*(COS(REAL(iao + 2*jao, wp)) &
               & - COS(REAL(jao + 2*iao, wp)))
         END DO
      END DO
      density_k(:, :, 1, 1) = CMPLX(amat, bmat, wp)
      density_k(:, :, 1, 2) = CONJG(density_k(:, :, 1, 1))
      cells = MATMUL(matinv_3x3(mol%lattice), images%translation)
      DO itr = 1, nimg
         angle = 2.0_wp*pi*DOT_PRODUCT(kfrac(:, 1), cells(:, itr))
         density_image(:, :, itr, 1) = COS(angle)*amat + SIN(angle)*bmat
         ! CP2K's gathered AO/image matrices carry the opposite cell phase
         ! relative to the provider's direct image-energy convention.
         density_cp2k(:, :, itr, 1) = COS(angle)*amat - SIN(angle)*bmat
      END DO

      CALL cp2k_acp_kmesh(calc, mol, images, kfrac, h_k, error)
      IF (ALLOCATED(error)) RETURN
      DO ik = 1, 2
         DO itr = 1, nimg
            angle = 2.0_wp*pi*DOT_PRODUCT(kfrac(:, ik), cells(:, itr))
            phase = EXP(CMPLX(0.0_wp, angle, wp))
            h_image(:, :, ik) = h_image(:, :, ik) + &
               & phase*images%hamiltonian(:, :, itr)
         END DO
      END DO
      CALL check(error, MAXVAL(ABS(h_k - h_image)), 0.0_wp, thr=1.0e-11_wp)
      IF (ALLOCATED(error)) THEN
         CALL test_failed(error, "Direct k-space ACP does not match image convolution")
         RETURN
      END IF
      CALL cp2k_build_acp_projectors(calc, mol, bcache, acache, lattr, projectors, projector_cache=projector_cache)
      CALL projector_cache_counts(projector_cache, cache_hits, cache_misses, cache_bypasses, cache_valid)
      CALL check(error, cache_misses, 1)
      IF (.NOT. ALLOCATED(error)) CALL check(error, cache_hits, 1)
      IF (ALLOCATED(error)) THEN
         CALL test_failed(error, "Identical ACP projector build did not hit the cache")
         RETURN
      END IF
      IF (ALLOCATED(projectors%translation) .OR. &
         & ALLOCATED(projectors%hamiltonian)) THEN
         CALL test_failed(error, &
            & "Compact ACP projector build materialized the Hamiltonian difference set")
         RETURN
      END IF
      CALL cp2k_acp_kmesh(calc, mol, projectors, kfrac, h_compact, error)
      IF (ALLOCATED(error)) RETURN
      IF (ANY(h_compact /= h_k)) THEN
         CALL test_failed(error, "Cached ACP projector build changed the Fock matrix")
         RETURN
      END IF
      CALL check(error, MAXVAL(ABS(h_compact - h_k)), 0.0_wp, thr=1.0e-11_wp)
      IF (ALLOCATED(error)) THEN
         CALL test_failed(error, &
            & "Compact ACP projector build changed the direct k-space Hamiltonian")
         RETURN
      END IF
      energy_k = 0.0_wp
      DO ik = 1, 2
         energy_k = energy_k + weights(ik)*REAL(SUM(CONJG(density_k(:, :, 1, ik)) &
                                                   & *h_k(:, :, ik)), wp)
      END DO
      energy_cached = 0.0_wp
      DO ik = 1, 2
         energy_cached = energy_cached + weights(ik)*REAL(SUM( &
            & CONJG(density_k(:, :, 1, ik))*h_compact(:, :, ik)), wp)
      END DO
      IF (energy_cached /= energy_k) THEN
         CALL test_failed(error, "Cached ACP projector build changed the energy")
         RETURN
      END IF
      energy_image = contract_acp_image_energy(images, images%translation, &
         & density_image)
      CALL check(error, energy_k, energy_image, thr=1.0e-11_wp)
      IF (ALLOCATED(error)) THEN
         CALL test_failed(error, &
            & "Weighted k-space ACP energy does not match the complete image reference")
         RETURN
      END IF

      ALLOCATE (dEdcn_cached(mol%nat), dEdcn_image(mol%nat), dEdcn_k(mol%nat), &
         & dEdcn_sparse(mol%nat), dEdcn_stream(mol%nat), &
         & dEdq_cached(mol%nat), dEdq_image(mol%nat), dEdq_k(mol%nat), &
         & dEdq_sparse(mol%nat), dEdq_stream(mol%nat), &
         & gradient_cached(3, mol%nat), gradient_image(3, mol%nat), &
         & gradient_k(3, mol%nat), gradient_sparse(3, mol%nat), &
         & gradient_streamed(3, mol%nat), source=0.0_wp)
      sigma_cached = 0.0_wp
      sigma_image = 0.0_wp
      sigma_k = 0.0_wp
      sigma_sparse = 0.0_wp
      sigma_stream = 0.0_wp
      CALL cp2k_acp_image_gradient(calc, mol, bcache, acache, images, &
         & images%translation, density_image, dEdcn_image, dEdq_image, &
         & gradient_image, sigma_image, error)
      IF (ALLOCATED(error)) RETURN
      CALL cp2k_acp_projector_image_gradient(calc, mol, bcache, acache, images, &
         & images%translation, density_cp2k, dEdcn_sparse, dEdq_sparse, &
         & gradient_sparse, sigma_sparse, error)
      IF (ALLOCATED(error)) RETURN
      CALL cp2k_acp_kmesh_gradient(calc, mol, bcache, acache, images, &
         & kfrac, weights, density_k, dEdcn_k, dEdq_k, gradient_k, sigma_k, error)
      IF (ALLOCATED(error)) RETURN
      CALL cp2k_acp_kmesh_gradient(calc, mol, bcache, acache, projectors, &
         & kfrac, weights, density_k, dEdcn_cached, dEdq_cached, &
         & gradient_cached, sigma_cached, error)
      IF (ALLOCATED(error)) RETURN
      ! Independently evaluated OpenMP derivative reductions need not add their
      ! thread-local contributions in the same order. Keep a roundoff-level
      ! numerical gate; bitwise equality remains required for cached Fock data.
      CALL check(error, MAX(MAXVAL(ABS(dEdcn_cached - dEdcn_k)), &
         & MAXVAL(ABS(dEdq_cached - dEdq_k)), MAXVAL(ABS(gradient_cached - gradient_k)), &
         & MAXVAL(ABS(sigma_cached - sigma_k))), 0.0_wp, thr=16*EPSILON(1.0_wp))
      IF (ALLOCATED(error)) THEN
         WRITE (*, '(a,4es24.15)') 'ACP cached response residuals: ', &
            & MAXVAL(ABS(dEdcn_cached - dEdcn_k)), MAXVAL(ABS(dEdq_cached - dEdq_k)), &
            & MAXVAL(ABS(gradient_cached - gradient_k)), MAXVAL(ABS(sigma_cached - sigma_k))
         CALL test_failed(error, &
            & "Cached ACP projector build changed the force or stress response")
         RETURN
      END IF
      CALL check(error, MAX(MAXVAL(ABS(dEdcn_k - dEdcn_image)), &
         & MAXVAL(ABS(dEdq_k - dEdq_image)), &
         & MAXVAL(ABS(gradient_k - gradient_image)), &
         & MAXVAL(ABS(sigma_k - sigma_image))), 0.0_wp, thr=1.0e-11_wp)
      IF (ALLOCATED(error)) THEN
         CALL test_failed(error, &
            & "Direct complex k-space ACP reverse sweep does not match image response")
         RETURN
      END IF
      CALL check(error, MAX(MAXVAL(ABS(dEdcn_sparse - dEdcn_k)), &
         & MAXVAL(ABS(dEdq_sparse - dEdq_k)), &
         & MAXVAL(ABS(gradient_sparse - gradient_k)), &
         & MAXVAL(ABS(sigma_sparse - sigma_k))), 0.0_wp, thr=1.0e-11_wp)
      IF (ALLOCATED(error)) THEN
         CALL test_failed(error, &
            & "Sparse CP2K-image ACP derivative changed the complex k-space response")
         RETURN
      END IF

      ! The bounded reverse stream must reproduce the whole-mesh reverse sweep
      ! while owning no AO-by-AO-by-k density or Hamiltonian payload.
      CALL cp2k_acp_kmesh_gradient_stream_begin(gradient_stream, calc, mol, &
         & projectors, 1, 2, error)
      IF (ALLOCATED(error)) RETURN
      CALL check(error, cp2k_acp_kmesh_gradient_stream_real_elements( &
         & gradient_stream), 3*SIZE(projectors%projector_translation, 2) + &
         & calc%acp%auxbas%nao*nao*SIZE(projectors%projector_translation, 2))
      IF (.NOT. ALLOCATED(error)) THEN
         CALL check(error, cp2k_acp_kmesh_gradient_stream_complex_elements( &
            & gradient_stream), 2*calc%acp%auxbas%nao*nao)
      END IF
      IF (ALLOCATED(error)) THEN
         CALL test_failed(error, "ACP reverse stream storage is not mesh bounded")
         RETURN
      END IF
      DO ik = 1, 2
         CALL cp2k_acp_kmesh_gradient_stream_push(gradient_stream, projectors, &
            & ik, kfrac(:, ik), weights(ik), density_k(:, :, :, ik), error)
         IF (ALLOCATED(error)) RETURN
      END DO
      CALL cp2k_acp_kmesh_gradient_stream_end(gradient_stream, calc, mol, &
         & bcache, acache, projectors, dEdcn_stream, dEdq_stream, &
         & gradient_streamed, sigma_stream, error)
      IF (ALLOCATED(error)) RETURN
      CALL check(error, MAX(MAXVAL(ABS(dEdcn_stream - dEdcn_k)), &
         & MAXVAL(ABS(dEdq_stream - dEdq_k)), &
         & MAXVAL(ABS(gradient_streamed - gradient_k)), &
         & MAXVAL(ABS(sigma_stream - sigma_k))), 0.0_wp, thr=1.0e-11_wp)
      IF (ALLOCATED(error)) THEN
         CALL test_failed(error, &
            & "Streamed ACP reverse sweep changed force or stress response")
         RETURN
      END IF

      ! The stream is fail-closed: no reordering, incomplete transaction or
      ! unnormalized irreducible-star weights may reach the derivative kernel.
      CALL cp2k_acp_kmesh_gradient_stream_begin(gradient_stream, calc, mol, &
         & projectors, 1, 2, error)
      IF (ALLOCATED(error)) RETURN
      CALL cp2k_acp_kmesh_gradient_stream_push(gradient_stream, projectors, &
         & 2, kfrac(:, 2), weights(2), density_k(:, :, :, 2), error)
      IF (.NOT. ALLOCATED(error)) THEN
         CALL cp2k_acp_kmesh_gradient_stream_discard(gradient_stream)
         CALL test_failed(error, "Reordered ACP reverse-stream k point was accepted")
         RETURN
      END IF
      IF (TRIM(error%message) /= &
         & "ACP k-mesh gradient stream k-point order changed") THEN
         DEALLOCATE (error)
         CALL cp2k_acp_kmesh_gradient_stream_discard(gradient_stream)
         CALL test_failed(error, "Unexpected ACP reverse-stream order diagnostic")
         RETURN
      END IF
      DEALLOCATE (error)
      CALL cp2k_acp_kmesh_gradient_stream_discard(gradient_stream)

      CALL cp2k_acp_kmesh_gradient_stream_begin(gradient_stream, calc, mol, &
         & projectors, 1, 2, error)
      IF (ALLOCATED(error)) RETURN
      CALL cp2k_acp_kmesh_gradient_stream_push(gradient_stream, projectors, &
         & 1, kfrac(:, 1), weights(1), density_k(:, :, :, 1), error)
      IF (ALLOCATED(error)) RETURN
      CALL cp2k_acp_kmesh_gradient_stream_end(gradient_stream, calc, mol, &
         & bcache, acache, projectors, dEdcn_stream, dEdq_stream, &
         & gradient_streamed, sigma_stream, error)
      IF (.NOT. ALLOCATED(error)) THEN
         CALL cp2k_acp_kmesh_gradient_stream_discard(gradient_stream)
         CALL test_failed(error, "Incomplete ACP reverse stream was accepted")
         RETURN
      END IF
      IF (TRIM(error%message) /= &
         & "ACP k-mesh gradient stream ended with missing k points") THEN
         DEALLOCATE (error)
         CALL cp2k_acp_kmesh_gradient_stream_discard(gradient_stream)
         CALL test_failed(error, "Unexpected incomplete ACP stream diagnostic")
         RETURN
      END IF
      DEALLOCATE (error)
      CALL cp2k_acp_kmesh_gradient_stream_discard(gradient_stream)

      CALL cp2k_acp_kmesh_gradient_stream_begin(gradient_stream, calc, mol, &
         & projectors, 1, 2, error)
      IF (ALLOCATED(error)) RETURN
      DO ik = 1, 2
         CALL cp2k_acp_kmesh_gradient_stream_push(gradient_stream, projectors, &
            & ik, kfrac(:, ik), 0.5_wp*weights(ik), &
            & density_k(:, :, :, ik), error)
         IF (ALLOCATED(error)) RETURN
      END DO
      CALL cp2k_acp_kmesh_gradient_stream_end(gradient_stream, calc, mol, &
         & bcache, acache, projectors, dEdcn_stream, dEdq_stream, &
         & gradient_streamed, sigma_stream, error)
      IF (.NOT. ALLOCATED(error)) THEN
         CALL cp2k_acp_kmesh_gradient_stream_discard(gradient_stream)
         CALL test_failed(error, "Unnormalized ACP reverse stream was accepted")
         RETURN
      END IF
      IF (TRIM(error%message) /= &
         & "ACP k-mesh gradient stream weights are not normalized") THEN
         DEALLOCATE (error)
         CALL cp2k_acp_kmesh_gradient_stream_discard(gradient_stream)
         CALL test_failed(error, "Unexpected ACP stream weight diagnostic")
         RETURN
      END IF
      DEALLOCATE (error)
      CALL cp2k_acp_kmesh_gradient_stream_discard(gradient_stream)

      ! A compact projector derivative may consume CP2K's sparse density-image
      ! support directly.  Missing translation differences are exact zero blocks,
      ! so this must equal the established zero-padded difference-set oracle.
      origin = images%origin
      IF (origin < 1) THEN
         CALL test_failed(error, "ACP difference set has no origin image")
         RETURN
      END IF
      ALLOCATE (density_zero(nao, nao, nimg, 1), source=0.0_wp)
      ALLOCATE (density_sparse(nao, nao, 1, 1), source=0.0_wp)
      ALLOCATE (translation_sparse(3, 1))
      density_zero(:, :, origin, 1) = amat
      density_sparse(:, :, 1, 1) = amat
      translation_sparse(:, 1) = images%translation(:, origin)
      dEdcn_image = 0.0_wp
      dEdq_image = 0.0_wp
      gradient_image = 0.0_wp
      sigma_image = 0.0_wp
      dEdcn_sparse = 0.0_wp
      dEdq_sparse = 0.0_wp
      gradient_sparse = 0.0_wp
      sigma_sparse = 0.0_wp
      CALL cp2k_acp_image_gradient(calc, mol, bcache, acache, images, &
         & images%translation, density_zero, dEdcn_image, dEdq_image, &
         & gradient_image, sigma_image, error)
      IF (ALLOCATED(error)) RETURN
      CALL cp2k_acp_projector_image_gradient(calc, mol, bcache, acache, &
         & projectors, translation_sparse, density_sparse, dEdcn_sparse, &
         & dEdq_sparse, gradient_sparse, sigma_sparse, error)
      IF (ALLOCATED(error)) RETURN
      CALL check(error, MAX(MAXVAL(ABS(dEdcn_sparse - dEdcn_image)), &
         & MAXVAL(ABS(dEdq_sparse - dEdq_image)), &
         & MAXVAL(ABS(gradient_sparse - gradient_image)), &
         & MAXVAL(ABS(sigma_sparse - sigma_image))), 0.0_wp, thr=1.0e-11_wp)
      IF (ALLOCATED(error)) THEN
         CALL test_failed(error, &
            & "Sparse ACP projector derivative changed the image response")
         RETURN
      END IF
      WRITE (*, '(a,1x,i0,1x,a,1x,i0)') "ACP_STREAM_STATE_GATES", 3, &
         & "ACP_SPARSE_DENSITY_IMAGES", SIZE(translation_sparse, 2)

      ! A Cartesian displacement must invalidate both the geometry and effective
      ! q-vSZP signature, even when the lattice translation set is unchanged.
      mol_changed = mol
      mol_changed%xyz(1, 1) = mol_changed%xyz(1, 1) + 1.0e-6_wp
      CALL cp2k_update_basis(calc, mol_changed, bcache_changed, error, grad=.TRUE., &
         & wfn_aux=wfn_aux_changed)
      IF (ALLOCATED(error)) RETURN
      CALL cp2k_build_acp_projectors(calc, mol_changed, bcache_changed, acache, &
         & lattr, projectors_changed, projector_cache=projector_cache)
      CALL projector_cache_counts(projector_cache, cache_hits, cache_misses, cache_bypasses, cache_valid)
      CALL check(error, cache_misses, 2)
      IF (ALLOCATED(error)) THEN
         CALL test_failed(error, "ACP projector cache survived a geometry change")
         RETURN
      END IF

      ! A cell deformation is an independent invalidation key, including the
      ! regenerated lattice-translation range used by the projector integrals.
      mol_changed = mol
      mol_changed%lattice(1, 1) = mol_changed%lattice(1, 1) + 1.0e-6_wp
      CALL cp2k_update_basis(calc, mol_changed, bcache_changed, error, grad=.TRUE., &
         & wfn_aux=wfn_aux_changed)
      IF (ALLOCATED(error)) RETURN
      CALL get_lattice_points(mol_changed%periodic, mol_changed%lattice, &
         & cp2k_get_acp_cutoff(calc, acc), lattr_changed)
      CALL cp2k_build_acp_projectors(calc, mol_changed, bcache_changed, acache, &
         & lattr_changed, projectors_changed, projector_cache=projector_cache)
      CALL projector_cache_counts(projector_cache, cache_hits, cache_misses, cache_bypasses, cache_valid)
      CALL check(error, cache_misses, 3)
      IF (ALLOCATED(error)) THEN
         CALL test_failed(error, "ACP projector cache survived a cell change")
         RETURN
      END IF

      ! Re-prime the original geometry, then perturb only the live effective
      ! q-vSZP contraction.  This isolates basis-cache invalidation from geometry.
      CALL cp2k_build_acp_projectors(calc, mol, bcache, acache, lattr, &
         & projectors_changed, projector_cache=projector_cache)
      CALL projector_cache_counts(projector_cache, cache_hits, cache_misses, cache_bypasses, cache_valid)
      CALL check(error, cache_misses, 4)
      IF (ALLOCATED(error)) RETURN
      bcache_changed = bcache
      bcache_changed%cgto(1, 1)%norm = &
         & bcache_changed%cgto(1, 1)%norm*(1.0_wp + 1.0e-6_wp)
      CALL cp2k_build_acp_projectors(calc, mol, bcache_changed, acache, lattr, &
         & projectors_changed, projector_cache=projector_cache)
      CALL projector_cache_counts(projector_cache, cache_hits, cache_misses, cache_bypasses, cache_valid)
      CALL check(error, cache_misses, 5)
      IF (ALLOCATED(error)) THEN
         CALL test_failed(error, &
            & "ACP projector cache survived an effective-CGTO change")
         RETURN
      END IF

      ! Raw primitive metadata and ACP levels are independent parts of the key.
      CALL cp2k_build_acp_projectors(calc, mol, bcache, acache, lattr, &
         & projectors_changed, projector_cache=projector_cache)
      CALL projector_cache_counts(projector_cache, cache_hits, cache_misses, cache_bypasses, cache_valid)
      CALL check(error, cache_misses, 6)
      IF (ALLOCATED(error)) RETURN
      associate (cgto => calc%bas%cgto(1, mol%id(1))%raw)
         saved_alpha = cgto%alpha(1)
         cgto%alpha(1) = saved_alpha*(1.0_wp + 1.0e-6_wp)
      END associate
      CALL cp2k_build_acp_projectors(calc, mol, bcache, acache, lattr, &
         & projectors_changed, projector_cache=projector_cache)
      CALL projector_cache_counts(projector_cache, cache_hits, cache_misses, cache_bypasses, cache_valid)
      CALL check(error, cache_misses, 7)
      associate (cgto => calc%bas%cgto(1, mol%id(1))%raw)
         cgto%alpha(1) = saved_alpha
      END associate
      IF (ALLOCATED(error)) THEN
         CALL test_failed(error, "ACP projector cache survived a primitive change")
         RETURN
      END IF

      CALL cp2k_build_acp_projectors(calc, mol, bcache, acache, lattr, &
         & projectors_changed, projector_cache=projector_cache)
      CALL projector_cache_counts(projector_cache, cache_hits, cache_misses, cache_bypasses, cache_valid)
      CALL check(error, cache_misses, 8)
      IF (ALLOCATED(error)) RETURN
      saved_level = calc%acp%levels(1, mol%id(1))
      calc%acp%levels(1, mol%id(1)) = saved_level + 1.0e-6_wp
      CALL cp2k_build_acp_projectors(calc, mol, bcache, acache, lattr, &
         & projectors_changed, projector_cache=projector_cache)
      CALL projector_cache_counts(projector_cache, cache_hits, cache_misses, cache_bypasses, cache_valid)
      CALL check(error, cache_misses, 9)
      calc%acp%levels(1, mol%id(1)) = saved_level
      IF (ALLOCATED(error)) THEN
         CALL test_failed(error, "ACP projector cache survived an ACP-level change")
         RETURN
      END IF

      ! The explicit direct-build oracle must reproduce every ACP observable,
      ! invalidate the prior payload before rebuilding, and leave no reusable
      ! state behind.  The following normal build therefore has to miss once and
      ! only the subsequent identical build may hit.
      CALL cp2k_build_acp_projectors(calc, mol, bcache, acache, lattr, &
         & projectors_changed, use_projector_cache=.FALSE., projector_cache=projector_cache)
      CALL projector_cache_counts(projector_cache, cache_hits, cache_misses, cache_bypasses, cache_valid)
      CALL check(error, cache_bypasses, 1)
      IF (.NOT. ALLOCATED(error)) CALL check(error,.NOT. cache_valid)
      IF (ALLOCATED(error)) THEN
         CALL test_failed(error, &
            & "Direct ACP oracle retained a stale projector cache")
         RETURN
      END IF
      CALL cp2k_acp_kmesh(calc, mol, projectors_changed, kfrac, h_compact, error)
      IF (ALLOCATED(error)) RETURN
      IF (ANY(h_compact /= h_k)) THEN
         CALL test_failed(error, "Direct ACP oracle changed the Fock matrix")
         RETURN
      END IF
      energy_cached = 0.0_wp
      DO ik = 1, 2
         energy_cached = energy_cached + weights(ik)*REAL(SUM( &
            & CONJG(density_k(:, :, 1, ik))*h_compact(:, :, ik)), wp)
      END DO
      CALL check(error, energy_cached, energy_k, thr=1.0e-15_wp)
      IF (ALLOCATED(error)) THEN
         CALL test_failed(error, "Direct ACP oracle changed the energy")
         RETURN
      END IF
      dEdcn_cached = 0.0_wp
      dEdq_cached = 0.0_wp
      gradient_cached = 0.0_wp
      sigma_cached = 0.0_wp
      CALL cp2k_acp_kmesh_gradient(calc, mol, bcache, acache, &
         & projectors_changed, kfrac, weights, density_k, dEdcn_cached, &
         & dEdq_cached, gradient_cached, sigma_cached, error)
      IF (ALLOCATED(error)) RETURN
      CALL check(error, MAX(MAXVAL(ABS(dEdcn_cached - dEdcn_k)), &
         & MAXVAL(ABS(dEdq_cached - dEdq_k)), MAXVAL(ABS(gradient_cached - gradient_k)), &
         & MAXVAL(ABS(sigma_cached - sigma_k))), 0.0_wp, thr=16*EPSILON(1.0_wp))
      IF (ALLOCATED(error)) THEN
         WRITE (*, '(a,4es24.15)') 'ACP direct response residuals: ', &
            & MAXVAL(ABS(dEdcn_cached - dEdcn_k)), MAXVAL(ABS(dEdq_cached - dEdq_k)), &
            & MAXVAL(ABS(gradient_cached - gradient_k)), MAXVAL(ABS(sigma_cached - sigma_k))
         CALL test_failed(error, &
            & "Direct ACP oracle changed the force or stress response")
         RETURN
      END IF

      CALL cp2k_build_acp_projectors(calc, mol, bcache, acache, lattr, &
         & projectors_changed, projector_cache=projector_cache)
      CALL projector_cache_counts(projector_cache, cache_hits, cache_misses, cache_bypasses, cache_valid)
      CALL check(error, cache_misses, 10)
      IF (.NOT. ALLOCATED(error)) CALL check(error, cache_valid)
      IF (ALLOCATED(error)) THEN
         CALL test_failed(error, "ACP projector cache was stale after direct rebuild")
         RETURN
      END IF
      CALL cp2k_build_acp_projectors(calc, mol, bcache, acache, lattr, &
         & projectors_changed, projector_cache=projector_cache)
      CALL projector_cache_counts(projector_cache, cache_hits, cache_misses, cache_bypasses, cache_valid)
      CALL check(error, cache_hits, 2)
      IF (ALLOCATED(error)) THEN
         CALL test_failed(error, "ACP projector cache did not recover after bypass")
         RETURN
      END IF
   END SUBROUTINE test_kmesh_acp_h2o_pbc

! **************************************************************************************************
!> \brief ...
!> \param error ...
! **************************************************************************************************
   SUBROUTINE test_cp2k_acp_cutoff(error)
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      REAL(wp)                                           :: conservative, cutoff, expected, &
                                                            pair_cutoff, valence
      TYPE(basis_type)                                   :: compact, diffuse
      TYPE(structure_type)                               :: mol
      TYPE(xtb_calculator)                               :: calc

      CALL get_structure(mol, "X23", "urea")
      CALL new_gxtb_calculator(calc, mol, error, accuracy=acc)
      IF (ALLOCATED(error)) RETURN

      valence = get_cutoff(calc%bas, acc)
      cutoff = cp2k_get_acp_cutoff(calc, acc)
      conservative = MAX(valence, get_cutoff(calc%acp%auxbas, acc))
      expected = MAX(valence, get_pair_cutoff(calc%bas, calc%acp%auxbas, acc))
      cutoff = cp2k_get_acp_cutoff(calc, acc)
      CALL check(error, cutoff, expected, thr=thr)
      IF (ALLOCATED(error)) THEN
         CALL test_failed(error, &
            & "CP2K ACP cutoff does not cover valence and auxiliary bases")
         RETURN
      END IF
      ! For the current H/C/N/O g-xTB parameters the q-vSZP valence basis is
      ! already the more diffuse partner, so the mixed cutoff must not inflate
      ! the established valence range.
      CALL check(error, cutoff, valence, thr=thr)
      IF (ALLOCATED(error)) THEN
         CALL test_failed(error, &
            & "H/C/N/O ACP cutoff unexpectedly exceeds the valence range")
         RETURN
      END IF
      CALL check(error, cutoff <= conservative)
      IF (ALLOCATED(error)) THEN
         CALL test_failed(error, &
            & "mixed ACP cutoff exceeds the conservative same-basis bound")
         RETURN
      END IF

      ! Exercise the opposite exponent ordering independently of the current
      ! parametrization: a diffuse auxiliary basis must enlarge the required
      ! range, but never beyond the conservative same-basis maximum.
      compact%min_alpha = 1.0_wp
      diffuse%min_alpha = 0.25_wp
      valence = get_cutoff(compact, acc)
      pair_cutoff = get_pair_cutoff(compact, diffuse, acc)
      conservative = MAX(valence, get_cutoff(diffuse, acc))
      CALL check(error, pair_cutoff > valence .AND. pair_cutoff <= conservative)
      IF (ALLOCATED(error)) THEN
         CALL test_failed(error, "mixed Gaussian cutoff does not respect safe bounds")
      END IF
   END SUBROUTINE test_cp2k_acp_cutoff

! **************************************************************************************************
!> \brief ...
!> \param error ...
! **************************************************************************************************
   SUBROUTINE test_cp2k_scc_mixer_dimension(error)
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      INTEGER                                            :: ndim
      TYPE(structure_type)                               :: mol
      TYPE(xtb_calculator)                               :: calc

      CALL get_structure(mol, "MB16-43", "H2")
      CALL new_gxtb_calculator(calc, mol, error, accuracy=acc)
      IF (ALLOCATED(error)) RETURN

      ndim = cp2k_get_scc_mixer_dimension(mol, calc%bas, calc%variable_info())
      CALL check(error, ndim, calc%bas%nsh + 9*mol%nat)
      IF (ALLOCATED(error)) THEN
         CALL test_failed(error, "CP2K mixer dimension includes CP2K-owned AO density")
      END IF
   END SUBROUTINE test_cp2k_scc_mixer_dimension

! **************************************************************************************************
!> \brief ...
!> \param error ...
! **************************************************************************************************
   SUBROUTINE test_cp2k_potential_mixer(error)
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      REAL(wp), PARAMETER                                :: rms = SQRT(1.0_wp/3.0_wp)

      REAL(wp)                                           :: mixed(3)
      TYPE(cp2k_potential_mixer_type)                    :: mixer

      CALL new_cp2k_potential_mixer(mixer, 3, error)
      IF (ALLOCATED(error)) RETURN

      CALL mixer%advance([0.0_wp, 0.0_wp, 0.0_wp], mixed, error)
      IF (ALLOCATED(error)) RETURN
      CALL check(error, MAXVAL(ABS(mixed)), 0.0_wp, thr=thr)
      IF (ALLOCATED(error)) RETURN
      CALL check(error, mixer%get_error(), 0.0_wp, thr=thr)
      IF (ALLOCATED(error)) RETURN

      CALL mixer%advance([1.0_wp, 0.0_wp, 0.0_wp], mixed, error)
      IF (ALLOCATED(error)) RETURN
      CALL check(error, MAXVAL(ABS(mixed - [0.2_wp, 0.0_wp, 0.0_wp])), &
         & 0.0_wp, thr=thr)
      IF (ALLOCATED(error)) RETURN
      CALL check(error, mixer%get_error(), rms, thr=thr)
      IF (ALLOCATED(error)) RETURN

      CALL mixer%advance([0.2_wp, 1.0_wp, 0.0_wp], mixed, error)
      IF (ALLOCATED(error)) RETURN
      CALL check(error, MAXVAL(ABS(mixed - [0.2_wp, 0.2_wp, 0.0_wp])), &
         & 0.0_wp, thr=thr)
      IF (ALLOCATED(error)) RETURN
      CALL check(error, mixer%get_error(), rms, thr=thr)
      IF (ALLOCATED(error)) RETURN

      ! The fourth Fock build switches to DIIS.  The two orthogonal unit
      ! residuals have equal Pulay weights, so the extrapolated vector is the
      ! arithmetic mean of the two corresponding raw Fock vectors.
      CALL mixer%advance([1.2_wp, 0.2_wp, 0.0_wp], mixed, error)
      IF (ALLOCATED(error)) RETURN
      CALL check(error, MAXVAL(ABS(mixed - [0.7_wp, 0.6_wp, 0.0_wp])), &
         & 0.0_wp, thr=thr)
      IF (ALLOCATED(error)) RETURN
      CALL check(error, mixer%get_error(), rms, thr=thr)
      IF (ALLOCATED(error)) RETURN

      CALL mixer%reset(3, error)
      IF (ALLOCATED(error)) RETURN
      CALL mixer%advance([2.0_wp, 3.0_wp, 4.0_wp], mixed, error)
      IF (ALLOCATED(error)) RETURN
      CALL check(error, MAXVAL(ABS(mixed - [2.0_wp, 3.0_wp, 4.0_wp])), &
         & 0.0_wp, thr=thr)
      IF (ALLOCATED(error)) THEN
         CALL test_failed(error, "CP2K potential mixer reset did not clear its history")
      END IF
   END SUBROUTINE test_cp2k_potential_mixer

! **************************************************************************************************
!> \brief ...
!> \param error ...
! **************************************************************************************************
   SUBROUTINE test_cp2k_potential_mixer_uninitialized(error)
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      REAL(wp)                                           :: mixed(1)
      TYPE(cp2k_potential_mixer_type)                    :: mixer

      CALL mixer%advance([0.0_wp], mixed, error)
   END SUBROUTINE test_cp2k_potential_mixer_uninitialized

! **************************************************************************************************
!> \brief ...
!> \param error ...
! **************************************************************************************************
   SUBROUTINE test_cp2k_potential_mixer_invalid_dimension(error)
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      TYPE(cp2k_potential_mixer_type)                    :: mixer

      CALL new_cp2k_potential_mixer(mixer, 0, error)
   END SUBROUTINE test_cp2k_potential_mixer_invalid_dimension

! **************************************************************************************************
!> \brief ...
!> \param error ...
! **************************************************************************************************
   SUBROUTINE test_cp2k_potential_mixer_size_mismatch(error)
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      REAL(wp)                                           :: mixed(3)
      TYPE(cp2k_potential_mixer_type)                    :: mixer

      CALL new_cp2k_potential_mixer(mixer, 3, error)
      IF (ALLOCATED(error)) RETURN
      CALL mixer%advance([0.0_wp, 0.0_wp], mixed, error)
   END SUBROUTINE test_cp2k_potential_mixer_size_mismatch

! **************************************************************************************************
!> \brief ...
!> \param error ...
! **************************************************************************************************
   SUBROUTINE test_cp2k_potential_mixer_nonfinite(error)
      TYPE(error_type), ALLOCATABLE, INTENT(out)         :: error

      REAL(wp)                                           :: mixed(1), raw(1)
      TYPE(cp2k_potential_mixer_type)                    :: mixer

      CALL new_cp2k_potential_mixer(mixer, 1, error)
      IF (ALLOCATED(error)) RETURN
      raw(1) = ieee_value(0.0_wp, ieee_quiet_nan)
      CALL mixer%advance(raw, mixed, error)
   END SUBROUTINE test_cp2k_potential_mixer_nonfinite

! **************************************************************************************************
!> \brief ...
!> \param images ...
!> \param density_translation ...
!> \param density ...
!> \return ...
! **************************************************************************************************
   FUNCTION contract_acp_image_energy(images, density_translation, density) &
      & RESULT(energy)
      TYPE(acp_image_type), INTENT(in)                   :: images
      REAL(wp), INTENT(in)                               :: density_translation(:, :), &
                                                            density(:, :, :, :)
      REAL(wp)                                           :: energy

      INTEGER                                            :: dtr, itr, jtr, spin
      REAL(wp)                                           :: scale, tol2

      energy = 0.0_wp
      scale = MAX(1.0_wp, MAXVAL(ABS(images%translation)), &
         & MAXVAL(ABS(density_translation)))
      tol2 = (2048.0_wp*EPSILON(1.0_wp)*scale)**2
      DO itr = 1, SIZE(images%translation, 2)
         dtr = 0
         DO jtr = 1, SIZE(density_translation, 2)
            IF (SUM((images%translation(:, itr) &
               & - density_translation(:, jtr))**2) <= tol2) THEN
               dtr = jtr
               EXIT
            END IF
         END DO
         IF (dtr == 0) THEN
            energy = HUGE(energy)
            RETURN
         END IF
         DO spin = 1, SIZE(density, 4)
            energy = energy + SUM(images%hamiltonian(:, :, itr) &
                                 & *density(:, :, dtr, spin))
         END DO
      END DO
   END FUNCTION contract_acp_image_energy

! These two existing fixtures are embedded to avoid adding mstore to CP2K.
! **************************************************************************************************
!> \brief ...
!> \param mol ...
!> \param collection ...
!> \param name ...
! **************************************************************************************************
   SUBROUTINE get_structure(mol, collection, name)
      TYPE(structure_type), INTENT(out)                  :: mol
      CHARACTER(len=*), INTENT(in)                       :: collection, name

      IF (collection == "X23" .AND. name == "urea") THEN
         CALL urea(mol)
      ELSE IF (collection == "MB16-43" .AND. name == "H2") THEN
         CALL H2(mol)
      ELSE
         error STOP "Unknown model-bridge test fixture"
      END IF
   END SUBROUTINE get_structure

! This file is part of mstore.
! SPDX-Identifier: Apache-2.0
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

! **************************************************************************************************
!> \brief ...
!> \param self ...
! **************************************************************************************************
   SUBROUTINE urea(self)
      TYPE(structure_type), INTENT(out)                  :: self

      INTEGER, PARAMETER                                 :: nat = 16
      CHARACTER(len=*), PARAMETER :: sym(nat) = [CHARACTER(len=4):: "H", "H", "H", "H", "H", "H", &
         "H", "H", "C", "C", "N", "N", "N", "N", "O", "O"]
      REAL(wp), PARAMETER :: lattice(3, 3) = RESHAPE([10.50481250525348_wp, 0.00000000000000_wp, &
         0.00000000000000_wp, 0.00000000000000_wp, 10.50481250525348_wp, 0.00000000000000_wp, &
         0.00000000000000_wp, 0.00000000000000_wp, 8.84967734421642_wp], SHAPE(lattice))
      REAL(wp), PARAMETER :: xyz(3, nat) = RESHAPE([2.73311029549004_wp, 7.98560290858176_wp, &
         2.52731916558693_wp, 7.98560290858176_wp, 7.77168595812233_wp, 6.32245531053102_wp, &
         7.77168595812233_wp, 2.51919334503061_wp, 2.52731916558693_wp, 2.51919334503061_wp, &
         2.73311029549004_wp, 6.32245531053102_wp, 1.52274097774035_wp, 6.77523359083207_wp, &
         8.53343439213466_wp, 6.77523359083207_wp, 8.98205527587201_wp, 0.31615111141221_wp, &
         8.98205527587201_wp, 3.72956266278029_wp, 8.53343439213466_wp, 3.72956266278029_wp, &
         1.52274097774035_wp, 0.31615111141221_wp, 0.00000000000000_wp, 5.25249261309172_wp, &
         2.92453950999127_wp, 5.25249261309172_wp, 0.00000000000000_wp, 5.92523496612667_wp, &
         1.54107131713484_wp, 6.79337495765549_wp, 1.60305432044818_wp, 6.79337495765549_wp, &
         8.96372493647752_wp, 7.24653118309870_wp, 8.96372493647752_wp, 3.71142129595688_wp, &
         1.60305432044818_wp, 3.71142129595688_wp, 1.54107131713484_wp, 7.24653118309870_wp, &
         0.00000000000000_wp, 5.25249261309172_wp, 5.32978239466231_wp, 5.25249261309172_wp, &
         0.00000000000000_wp, 3.51999208145564_wp], SHAPE(xyz))

      CALL new(self, sym, xyz, lattice=lattice)
   END SUBROUTINE urea

! **************************************************************************************************
!> \brief ...
!> \param self ...
! **************************************************************************************************
   SUBROUTINE H2(self)
      TYPE(structure_type), INTENT(out)                  :: self

      INTEGER, PARAMETER                                 :: nat = 2
      CHARACTER(len=*), PARAMETER :: sym(nat) = [CHARACTER(len=4):: "H", "H"]
      REAL(wp), PARAMETER :: xyz(3, nat) = RESHAPE([0.00000000000000_wp, 0.00000000000000_wp, &
         -0.70252931147690_wp, 0.00000000000000_wp, 0.00000000000000_wp, 0.70252931147690_wp], &
         SHAPE(xyz))

      CALL new(self, sym, xyz)
   END SUBROUTINE H2

END MODULE tblite_model_bridge_tests
#endif

PROGRAM tblite_model_bridge_unittest
#if defined(__TBLITE_GXTB_ACCELERATION)
   USE, INTRINSIC :: iso_fortran_env, ONLY: output_unit
   USE mctc_env_testing, ONLY: run_testsuite
   USE tblite_model_bridge_tests, ONLY: collect_model_bridge
#endif
   IMPLICIT NONE
#if defined(__TBLITE_GXTB_ACCELERATION)
   INTEGER :: failures
   failures = 0
   CALL run_testsuite(collect_model_bridge, output_unit, failures)
   IF (failures /= 0) error STOP "CP2K model-bridge tests failed"
   PRINT '(a)', 'CP2K model-bridge integration tests passed'
#else
   PRINT '(a)', 'CP2K model-bridge integration tests skipped: g-xTB response API unavailable'
#endif
END PROGRAM tblite_model_bridge_unittest
