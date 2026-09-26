!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
PROGRAM tblite_bvk_plan_unittest
#if defined(__TBLITE_GXTB)
   USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_quiet_nan, ieee_value
   USE kinds, ONLY: dp
   USE kpoint_lattice_fft, ONLY: kpoint_mesh_certificate_type
   USE mctc_env, ONLY: error_type
   USE mctc_io, ONLY: new, structure_type
   USE tblite_bvk_plan, ONLY: cp2k_bvk_plan_type, prepare_bvk_kernel, evaluate_bvk_exchange, &
                              bvk_kernel_snapshot, bvk_plan_counts, invalidate_bvk_plan, evaluate_bvk_gradient
   USE tblite_container, ONLY: container_cache
   USE tblite_model_bridge, ONLY: cp2k_prepare_exchange, cp2k_exchange_kmesh, cp2k_exchange_kmesh_gradient
   USE tblite_exchange_cache, ONLY: exchange_cache
   USE tblite_exchange_contractions, ONLY: get_KFock_kmesh
   USE tblite_exchange_fock, ONLY: exchange_bvk_kernel, exchange_fock
   USE tblite_stream_plan, ONLY: cp2k_stream_plan_type
   USE tblite_wavefunction_type, ONLY: wavefunction_type, new_wavefunction
   USE tblite_xtb_calculator, ONLY: xtb_calculator
   USE tblite_xtb_gxtb, ONLY: new_gxtb_calculator
#if defined(__TBLITE_GXTB_TRANSFORM)
   USE fft_lib, ONLY: fft_supports_arbitrary_lengths
   USE fft_tools, ONLY: init_fft, finalize_fft
   USE mathconstants, ONLY: gaussi, twopi
   USE message_passing, ONLY: mp_comm_type, mp_world_init, mp_world_finalize
   USE tblite_exchange_fft, ONLY: cp2k_exchange_transform_type, prepare_exchange_transform, &
                                  release_exchange_transform, exchange_transform_counts
   USE tblite_transform_contract, ONLY: exchange_k_to_r, exchange_r_to_k, &
                                            exchange_r_to_k_adjoint, exchange_k_to_r_adjoint
#endif
   IMPLICIT NONE
   TYPE(structure_type) :: mol, original_mol
   TYPE(xtb_calculator) :: calc
   TYPE(exchange_fock) :: model, original_model
   TYPE(exchange_bvk_kernel) :: kernel, oracle
   TYPE(cp2k_bvk_plan_type) :: plan
   TYPE(kpoint_mesh_certificate_type) :: certificate
   TYPE(container_cache) :: live
   TYPE(wavefunction_type) :: wfn
   TYPE(error_type), ALLOCATABLE :: error
   INTEGER :: builds, hits, before, variant, mesh(3), axis, cell, a, b, spin, nspin, nao, nk
   LOGICAL :: reused, valid
   REAL(KIND=dp) :: lattice(3, 3), xyz(3, 3), energy, oracle_energy
   REAL(KIND=dp), ALLOCATABLE :: kfrac(:, :), weights(:), vsh(:), oracle_vsh(:)
   COMPLEX(KIND=dp), ALLOCATABLE :: density(:, :, :, :), overlap(:, :, :), &
                                    fock(:, :, :, :), oracle_fock(:, :, :, :), first_fock(:, :, :, :)
#if defined(__TBLITE_GXTB_TRANSFORM)
   CHARACTER(LEN=16) :: backend
   TYPE(mp_comm_type) :: mp_comm

   CALL mp_world_init(mp_comm)
   backend = 'FFTSG'
   IF (COMMAND_ARGUMENT_COUNT() > 0) CALL GET_COMMAND_ARGUMENT(1, backend)
   CALL init_fft(TRIM(backend), alltoall=.FALSE., pool_limit=4, wisdom_file='', plan_style=1)
   CALL check_transforms()
#endif

   lattice = 0.0_dp
   lattice(1, 1) = 15.0_dp
   lattice(2, 2) = 16.0_dp
   lattice(3, 3) = 17.0_dp
   xyz = RESHAPE([0.0_dp, 0.0_dp, 0.0_dp, 1.8_dp, 0.0_dp, 0.0_dp, &
                  -0.5_dp, 1.7_dp, 0.1_dp], [3, 3])
   CALL new(mol, [8, 1, 1], xyz, lattice=lattice, periodic=[.TRUE., .TRUE., .TRUE.])
   CALL new_gxtb_calculator(calc, mol, error, accuracy=0.1_dp)
   CALL check(.NOT. ALLOCATED(error), 'calculator construction')
   SELECT TYPE (exchange => calc%exchange)
   TYPE IS (exchange_fock)
      model = exchange
   CLASS DEFAULT
      STOP 1
   END SELECT
   original_model = model
   original_mol = mol
   mesh = [3, 2, 1]
   CALL prepare_bvk_kernel(plan, model, mol, mesh, reused, error)
   CALL check(.NOT. ALLOCATED(error) .AND. .NOT. reused, 'first build')
   CALL compare_kernel()
   CALL prepare_bvk_kernel(plan, model, mol, mesh, reused, error)
   CALL check(.NOT. ALLOCATED(error) .AND. reused, 'repeat hit')
   CALL bvk_plan_counts(plan, builds, hits, valid)
   CALL check(builds == 1 .AND. hits == 1 .AND. valid, 'initial counters')
   kernel%g_mulliken_r = -999.0_dp
   CALL compare_kernel()

   ! Change every static parameter group independently, retaining valid provider layouts.
   DO variant = 1, 24
      mol = original_mol
      model = original_model
      mesh = [3, 2, 1]
      CALL prepare_bvk_kernel(plan, model, mol, mesh, reused, error)
      CALL check(.NOT. ALLOCATED(error), 'restore baseline')
      CALL bvk_plan_counts(plan, before, hits, valid)
      SELECT CASE (variant)
      CASE (1)
         mol%xyz(1, 2) = mol%xyz(1, 2) + 1.0E-13_dp
      CASE (2)
         mol%lattice(1, 1) = mol%lattice(1, 1) + 1.0E-13_dp
      CASE (3)
         mol%periodic(3) = .FALSE.
      CASE (4)
         mesh = [2, 3, 1]
      CASE (5)
         model%frscale = model%frscale + 1.0E-13_dp
      CASE (6)
         model%onsite_frscale = model%onsite_frscale + 1.0E-13_dp
      CASE (7)
         model%omega = model%omega + 1.0E-13_dp
      CASE (8)
         model%lrscale = model%lrscale + 1.0E-13_dp
      CASE (9)
         model%ondiag_scale = model%ondiag_scale + 1.0E-13_dp
      CASE (10)
         model%hubbard_exp = model%hubbard_exp + 1.0E-13_dp
      CASE (11)
         model%hubbard_exp_r0 = model%hubbard_exp_r0 + 1.0E-13_dp
      CASE (12)
         model%gexp = model%gexp + 1.0E-13_dp
      CASE (13)
         model%corr_exp = model%corr_exp + 1.0E-13_dp
      CASE (14)
         model%hubbard(1, 1, 1, 1) = model%hubbard(1, 1, 1, 1) + 1.0E-13_dp
      CASE (15)
         model%offdiag_scale(1, 1, 1, 1) = model%offdiag_scale(1, 1, 1, 1) + 1.0E-13_dp
      CASE (16)
         model%onecxints(1, 1, 1) = model%onecxints(1, 1, 1) + 1.0E-13_dp
      CASE (17)
         model%rad(1, 1) = model%rad(1, 1) + 1.0E-13_dp
      CASE (18)
         model%kq(1, 1) = model%kq(1, 1) + 1.0E-13_dp
      CASE (19)
         model%corr_scale(1, 1) = model%corr_scale(1, 1) + 1.0E-13_dp
      CASE (20)
         model%corr_rad(1, 1) = model%corr_rad(1, 1) + 1.0E-13_dp
      CASE (21)
         ! This AO-only metadata does not enter the provider's radial kernel.
         model%nao = model%nao + 1
      CASE (22)
         model%nao_sh(1) = model%nao_sh(1) + 1
      CASE (23)
         model%iao_sh(1) = model%iao_sh(1) + 1
      CASE (24)
         mesh = [1, 1, 1]
      END SELECT
      CALL prepare_bvk_kernel(plan, model, mol, mesh, reused, error)
      CALL check(.NOT. ALLOCATED(error) .AND. .NOT. reused, 'changed exact key')
      CALL bvk_plan_counts(plan, builds, hits, valid)
      CALL check(valid .AND. builds == before + 1, 'one rebuild after a change')
      CALL compare_kernel()
      CALL prepare_bvk_kernel(plan, model, mol, mesh, reused, error)
      CALL check(.NOT. ALLOCATED(error) .AND. reused, 'new key reusable')
   END DO

   mol = original_mol
   model = original_model
   mesh = [3, 2, 1]
   CALL prepare_bvk_kernel(plan, model, mol, mesh, reused, error)
   mol%xyz(1, 1) = ieee_value(0.0_dp, ieee_quiet_nan)
   CALL prepare_bvk_kernel(plan, model, mol, mesh, reused, error)
   CALL check(ALLOCATED(error) .AND. .NOT. reused, 'nonfinite key rejected')
   CALL bvk_kernel_snapshot(plan, kernel, error)
   CALL check(ALLOCATED(error), 'bad key invalidates previous payload')
   mol = original_mol
   CALL prepare_bvk_kernel(plan, model, mol, mesh, reused, error)
   CALL check(.NOT. ALLOCATED(error) .AND. .NOT. reused, 'recover from invalid key')
   CALL prepare_bvk_kernel(plan, model, mol, [0, 1, 1], reused, error)
   CALL check(ALLOCATED(error), 'zero mesh rejected')
   CALL prepare_bvk_kernel(plan, model, mol, [HUGE(1), 2, 1], reused, error)
   CALL check(ALLOCATED(error), 'overflow rejected before allocation')
   mol%periodic(1) = .FALSE.
   CALL prepare_bvk_kernel(plan, model, mol, mesh, reused, error)
   CALL check(ALLOCATED(error), 'nonperiodic replication rejected')
   mol = original_mol
   CALL prepare_bvk_kernel(plan, model, mol, mesh, reused, error)
   CALL invalidate_bvk_plan(plan)
   CALL prepare_bvk_kernel(plan, model, mol, mesh, reused, error)
   CALL check(.NOT. ALLOCATED(error) .AND. .NOT. reused, 'explicit invalidation')

   ! Compare complete Fock, shell response and energy with the unchanged provider oracle.
   CALL bvk_kernel_snapshot(plan, kernel, error)
   nk = PRODUCT(mesh)
   nao = model%nao
   ALLOCATE (kfrac(3, nk), weights(nk), vsh(model%nsh), oracle_vsh(model%nsh))
   DO axis = 1, 3
      kfrac(axis, :) = (REAL(kernel%reps(axis, :), dp) + 0.17_dp)/REAL(mesh(axis), dp)
   END DO
   weights = 1.0_dp/REAL(nk, dp)
   DO nspin = 1, 2
      CALL new_wavefunction(wfn, mol%nat, model%nsh, nao, nspin, 0.0_dp)
      ALLOCATE (density(nao, nao, nspin, nk), overlap(nao, nao, nk), fock(nao, nao, nspin, nk), &
                oracle_fock(nao, nao, nspin, nk))
      DO cell = 1, nk
         DO a = 1, nao
            DO b = 1, a
               overlap(a, b, cell) = CMPLX(0.001_dp*(a + b), 0.0001_dp*(a - b)*cell, dp)
               overlap(b, a, cell) = CONJG(overlap(a, b, cell))
               DO spin = 1, nspin
                  density(a, b, spin, cell) = CMPLX(0.01_dp*(a + b + spin), 0.001_dp*(a - b)*cell, dp)
                  density(b, a, spin, cell) = CONJG(density(a, b, spin, cell))
               END DO
            END DO
            overlap(a, a, cell) = 1.0_dp
         END DO
      END DO
      CALL compare_exchange()
      first_fock = fock
      CALL bvk_plan_counts(plan, before, hits, valid)
      wfn%qsh(:, 1) = 0.12_dp
      CALL compare_exchange()
      CALL check(MAXVAL(ABS(fock - first_fock)) > 1.0E-10_dp, 'onsite response really changes with charge')
      CALL bvk_plan_counts(plan, builds, hits, valid)
      CALL check(builds == before, 'charge change must not rebuild static kernel')
      density = density*0.9_dp
      overlap = overlap*1.01_dp
      CALL compare_exchange()
      CALL bvk_plan_counts(plan, builds, hits, valid)
      CALL check(builds == before, 'live density and overlap not cached')
      DEALLOCATE (density, overlap, fock, oracle_fock)
   END DO
   WRITE (*, '(A)') 'CP2K BVK PLAN: exact-key invalidation, live charge and independent dense-oracle checks passed'
#if defined(__TBLITE_GXTB_TRANSFORM)
   CALL invalidate_bvk_plan(plan)
   CALL finalize_fft(mp_comm, '')
   CALL mp_world_finalize()
#endif

CONTAINS
   SUBROUTINE check(condition, message)
      LOGICAL, INTENT(IN)                                :: condition
      CHARACTER(LEN=*), INTENT(IN)                       :: message

      IF (.NOT. condition) THEN
         WRITE (*, '(A)') 'FAIL: '//message
         IF (ALLOCATED(error)) WRITE (*, '(A)') error%message
         STOP 1
      END IF
   END SUBROUTINE check

   SUBROUTINE compare_kernel()
      CALL bvk_kernel_snapshot(plan, kernel, error)
      CALL check(.NOT. ALLOCATED(error), 'snapshot')
      CALL model%get_bvk_Kmatrix(mol, mesh, oracle)
      CALL check(ALL(kernel%nmesh == oracle%nmesh), 'independent mesh')
      CALL check(ALL(kernel%reps == oracle%reps), 'independent representatives')
      CALL check(ALL(kernel%g_mulliken_r == oracle%g_mulliken_r), 'independent Mulliken kernel')
      CALL check(ALL(kernel%g_bocorr_r == oracle%g_bocorr_r), 'independent bond-order kernel')
   END SUBROUTINE compare_kernel

   SUBROUTINE compare_exchange()
      TYPE(cp2k_stream_plan_type) :: direct_plan
      COMPLEX(KIND=dp), ALLOCATABLE :: overlap_response(:, :, :), reference_response(:, :, :), &
         single_density(:, :, :, :), single_fock(:, :, :, :)
      INTEGER                                            :: channel
      REAL(KIND=dp) :: gradient(3, mol%nat), reference_gradient(3, mol%nat), &
         reference_stress(3, 3), separate_energy, single_energy, stress(3, 3)
      REAL(KIND=dp), ALLOCATABLE                         :: separate_vsh(:)

      CALL cp2k_prepare_exchange(calc, mol, live, wfn, error)
      CALL check(.NOT. ALLOCATED(error), 'prepare current charge state')
      SELECT TYPE (cache => live%raw)
      TYPE IS (exchange_cache)
         CALL evaluate_bvk_exchange(plan, certificate, model, mol, cache, mesh, kfrac, weights, &
                                    density, overlap, fock, vsh, energy, error)
      CLASS DEFAULT
         STOP 1
      END SELECT
      CALL check(.NOT. ALLOCATED(error), 'cached complete-mesh evaluation')
      CALL cp2k_exchange_kmesh(calc, mol, live, wfn, mesh, kfrac, weights, &
                               density, overlap, oracle_fock, oracle_vsh, oracle_energy, error)
      CALL check(.NOT. ALLOCATED(error), 'independent provider oracle')
      WRITE (*, '(A,3ES24.15)') 'EXCHANGE energy, oracle, difference: ', energy, oracle_energy, energy - oracle_energy
      WRITE (*, '(A,2ES24.15)') 'EXCHANGE maximum Fock/shell differences: ', &
         MAXVAL(ABS(fock - oracle_fock)), MAXVAL(ABS(vsh - oracle_vsh))
      IF (SIZE(density, 3) == 2) THEN
         ALLOCATE (single_density(nao, nao, 1, nk), single_fock(nao, nao, 1, nk), separate_vsh(model%nsh))
         separate_energy = 0.0_dp
         SELECT TYPE (cache => live%raw)
         TYPE IS (exchange_cache)
            DO channel = 1, 2
               single_density(:, :, 1, :) = density(:, :, channel, :)
               CALL get_KFock_kmesh(model, mol, cache, direct_plan, kernel, kfrac, weights, single_density, overlap, &
                                   single_fock, separate_vsh, single_energy)
               separate_energy = separate_energy + 2.0_dp*single_energy
               WRITE (*, '(A,I0,2ES24.15)') 'SPIN separability FFT/direct Fock differences, channel ', channel, &
                  MAXVAL(ABS(fock(:, :, channel, :) - 2.0_dp*single_fock(:, :, 1, :))), &
                  MAXVAL(ABS(oracle_fock(:, :, channel, :) - 2.0_dp*single_fock(:, :, 1, :)))
               CALL check(MAXVAL(ABS(fock(:, :, channel, :) - 2.0_dp*single_fock(:, :, 1, :))) < 1.0E-12_dp, &
                          'consumer two-spin Fock versus separately packed spin channels')
               CALL check(MAXVAL(ABS(oracle_fock(:, :, channel, :) - 2.0_dp*single_fock(:, :, 1, :))) < 1.0E-12_dp, &
                          'direct two-spin Fock versus separately packed spin channels')
            END DO
         END SELECT
         WRITE (*, '(A,3ES24.15)') 'SPIN separate, FFT, direct energy: ', separate_energy, energy, oracle_energy
         CALL check(ABS(separate_energy - energy) < 1.0E-12_dp, 'consumer two-spin energy separability')
      END IF
      CALL check(ABS(energy - oracle_energy) < 1.0E-12_dp, 'exchange energy')
      CALL check(MAXVAL(ABS(fock - oracle_fock)) < 1.0E-12_dp, 'exchange Fock')
      CALL check(MAXVAL(ABS(vsh - oracle_vsh)) < 1.0E-12_dp, 'exchange shell response')
      ALLOCATE (overlap_response(nao, nao, nk), reference_response(nao, nao, nk))
      SELECT TYPE (cache => live%raw)
      TYPE IS (exchange_cache)
         CALL evaluate_bvk_gradient(plan, certificate, model, mol, cache, mesh, kfrac, weights, &
                                    density, overlap, overlap_response, gradient, stress, error)
      CLASS DEFAULT
         STOP 1
      END SELECT
      CALL check(.NOT. ALLOCATED(error), 'consumer exchange gradient')
      CALL cp2k_exchange_kmesh_gradient(calc, mol, live, mesh, kfrac, weights, density, overlap, &
                                        reference_response, reference_gradient, reference_stress, error)
      CALL check(.NOT. ALLOCATED(error), 'provider exchange gradient')
      CALL check(MAXVAL(ABS(overlap_response - reference_response)) < 1.0E-12_dp, 'overlap response')
      CALL check(MAXVAL(ABS(gradient - reference_gradient)) < 1.0E-12_dp, 'nuclear exchange gradient')
      CALL check(MAXVAL(ABS(stress - reference_stress)) < 1.0E-12_dp, 'exchange strain derivative')
   END SUBROUTINE compare_exchange

#if defined(__TBLITE_GXTB_TRANSFORM)
   SUBROUTINE check_transforms()
      COMPLEX(KIND=dp)                                   :: factor, lhs, rhs
      COMPLEX(KIND=dp), ALLOCATABLE :: explicit_sum(:, :), fx(:, :), fy(:, :), &
         padded_source(:, :), padded_target(:, :), reference(:, :), x(:, :), y(:, :)
      INTEGER                                            :: count, dimensions(3), i, index0, j, k, &
                                                            nd, nf, old_nd, old_nf, operation, r, &
                                                            trial
      INTEGER, ALLOCATABLE                               :: reps(:, :)
      REAL(KIND=dp)                                      :: tolerance
      REAL(KIND=dp), ALLOCATABLE                         :: points(:, :), w(:)
      TYPE(cp2k_exchange_transform_type)                 :: direct, fast

      DO trial = 1, 5
         SELECT CASE (trial)
         CASE (1)
            dimensions = [1, 1, 1]
         CASE (2)
            dimensions = [2, 2, 2]
         CASE (3)
            dimensions = [3, 2, 1]
         CASE (4)
            dimensions = [7, 2, 1]
         CASE (5)
            dimensions = [3, 3, 3]
         END SELECT
         count = PRODUCT(dimensions)
         ALLOCATE (reps(3, count), points(3, count), w(count), x(6, count), y(6, count), &
                   fx(6, count), fy(6, count), reference(6, count), explicit_sum(6, count), &
                   padded_source(12, count), padded_target(12, count))
         DO k = 1, count
            index0 = k - 1
            reps(:, k) = [MODULO(index0, dimensions(1)), MODULO(index0/dimensions(1), dimensions(2)), &
                          index0/(dimensions(1)*dimensions(2))]
            points(:, k) = (REAL(reps(:, k), dp) + [0.17_dp, -0.23_dp, 0.11_dp])/REAL(dimensions, dp)
            DO j = 1, 3
               DO i = 1, 2
                  x(i + 2*(j - 1), k) = CMPLX(SIN(REAL(i + 3*j + 7*k, dp)), COS(REAL(3*i + j + 2*k, dp)), dp)
                  y(i + 2*(j - 1), k) = CMPLX(COS(REAL(2*i + j + 5*k, dp)), SIN(REAL(i + 2*j + 3*k, dp)), dp)
               END DO
            END DO
         END DO
         ! A complete but independently reordered/translated cell representation and point order.
         reps = reps(:, count:1:-1)
         reps(1, :) = reps(1, :) - dimensions(1)
         points = points(:, count:1:-1)
         w = 1.0_dp/REAL(count, dp)
         CALL prepare_exchange_transform(fast, points, dimensions, reps, w)
         CALL exchange_transform_counts(fast, old_nf, old_nd)
         CALL fast%apply(x, fx, exchange_k_to_r)
         CALL fast%apply(fx, fy, exchange_r_to_k)
         CALL check(MAXVAL(ABS(x - fy)) < 1.0E-12_dp, 'complex round trip')
         CALL exchange_transform_counts(fast, nf, nd)
         IF (trial == 4 .AND. .NOT. fft_supports_arbitrary_lengths()) THEN
            CALL check(nf == old_nf .AND. nd == old_nd + 2, 'unsupported prime length uses direct fallback')
         ELSE
            CALL check(nf == old_nf + 2 .AND. nd == old_nd, 'backend actually executed FFTs')
         END IF

         ! Nonuniform weights test the general adjoints, not a production BvK quadrature.
         DO k = 1, count
            w(k) = REAL(k, dp)
         END DO
         w = w/SUM(w)
         CALL prepare_exchange_transform(fast, points, dimensions, reps, w)
         CALL prepare_exchange_transform(direct, points, dimensions, reps, w, force_direct=.TRUE.)
         DO operation = 1, 4
            CALL fast%apply(x, fx, operation)
            CALL direct%apply(x, reference, operation)
            explicit_sum = CMPLX(0.0_dp, 0.0_dp, dp)
            DO r = 1, count
               DO k = 1, count
                  factor = EXP(gaussi*twopi*SUM(points(:, k)*REAL(reps(:, r), dp)))
                  IF (operation == exchange_k_to_r .OR. operation == exchange_r_to_k_adjoint) &
                     factor = CONJG(factor)
                  IF (operation == exchange_k_to_r .OR. operation == exchange_k_to_r_adjoint) factor = w(k)*factor
                  IF (operation == exchange_k_to_r .OR. operation == exchange_r_to_k_adjoint) THEN
                     explicit_sum(:, r) = explicit_sum(:, r) + factor*x(:, k)
                  ELSE
                     explicit_sum(:, k) = explicit_sum(:, k) + factor*x(:, r)
                  END IF
               END DO
            END DO
            tolerance = 1.0E-12_dp*MAX(1.0_dp, MAXVAL(ABS(explicit_sum)))
            CALL check(MAXVAL(ABS(fx - explicit_sum)) < tolerance, 'FFT versus explicit phase sum')
            CALL check(MAXVAL(ABS(reference - explicit_sum)) < tolerance, 'direct versus explicit phase sum')
            padded_source = CMPLX(-999.0_dp, 777.0_dp, dp)
            padded_target = CMPLX(555.0_dp, -333.0_dp, dp)
            padded_source(1:12:2, :) = x
            CALL fast%apply(padded_source(1:12:2, :), padded_target(1:12:2, :), operation)
            CALL check(MAXVAL(ABS(padded_target(1:12:2, :) - explicit_sum)) < tolerance, 'strided AO-pair views')
            CALL check(ALL(padded_target(2:12:2, :) == CMPLX(555.0_dp, -333.0_dp, dp)), 'untouched padding rows')
         END DO
         CALL fast%apply(x, fx, exchange_k_to_r)
         CALL fast%apply(y, fy, exchange_k_to_r_adjoint)
         lhs = SUM(CONJG(y)*fx)
         rhs = SUM(CONJG(fy)*x)
         CALL check(ABS(lhs - rhs) < 1.0E-12_dp*MAX(1.0_dp, ABS(lhs), ABS(rhs)), 'inverse-map adjoint')
         CALL fast%apply(x, fx, exchange_r_to_k)
         CALL fast%apply(y, fy, exchange_r_to_k_adjoint)
         lhs = SUM(CONJG(y)*fx)
         rhs = SUM(CONJG(fy)*x)
         CALL check(ABS(lhs - rhs) < 1.0E-12_dp*MAX(1.0_dp, ABS(lhs), ABS(rhs)), 'forward-map adjoint')

         ! A close but nonregular point must use its actual phase, not a rounded FFT mesh.
         IF (count > 1) THEN
            points(1, count) = points(1, count) + 1.0E-11_dp
            CALL prepare_exchange_transform(fast, points, dimensions, reps, w)
            CALL prepare_exchange_transform(direct, points, dimensions, reps, w, force_direct=.TRUE.)
            CALL exchange_transform_counts(fast, old_nf, old_nd)
            CALL fast%apply(x, fx, exchange_k_to_r)
            CALL direct%apply(x, reference, exchange_k_to_r)
            CALL exchange_transform_counts(fast, nf, nd)
            CALL check(nf == old_nf .AND. nd == old_nd + 1, 'near-regular mesh takes exact direct path')
            CALL check(MAXVAL(ABS(fx - reference)) < 1.0E-14_dp, 'near-regular actual phases')
         END IF
         DEALLOCATE (reps, points, w, x, y, fx, fy, reference, explicit_sum, padded_source, padded_target)
      END DO
      CALL release_exchange_transform(fast)
      CALL release_exchange_transform(direct)
      WRITE (*, '(A)') 'CP2K EXCHANGE TRANSFORMS: complex maps, adjoints, shifted/reordered grids and fallback passed'
   END SUBROUTINE check_transforms
#endif
#else
   IMPLICIT NONE
   WRITE (*, '(A)') 'CP2K BvK plan test skipped: g-xTB support is not enabled'
#endif
END PROGRAM tblite_bvk_plan_unittest
