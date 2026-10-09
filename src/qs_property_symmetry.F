!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Cell-gauge scalar/spinor reconstruction on independent property k-point meshes.
!> \author Thomas D. Kuehne, tkuehne@cp2k.org
! **************************************************************************************************
MODULE qs_property_symmetry
   USE cell_types,                      ONLY: cell_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_cfm_types,                    ONLY: cp_cfm_copy_general,&
                                              cp_cfm_create,&
                                              cp_cfm_get_submatrix,&
                                              cp_cfm_release,&
                                              cp_cfm_set_submatrix,&
                                              cp_cfm_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE input_cp2k_kpoints,              ONLY: use_k290_kpoint_symmetry,&
                                              use_spglib_kpoint_backend,&
                                              use_spglib_kpoint_symmetry
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE kpoint_k_r_trafo_simple,         ONLY: replicate_rs_matrices,&
                                              rs_to_kp
   USE kpoint_methods,                  ONLY: kpoint_env_initialize,&
                                              kpoint_init_cell_index,&
                                              kpoint_initialize,&
                                              kpoint_initialize_mo_set,&
                                              kpoint_initialize_mos
   USE kpoint_mo_symmetry_methods,      ONLY: kpoint_cell_gauge_action,&
                                              kpoint_same_periodic
   USE kpoint_types,                    ONLY: kpoint_create,&
                                              kpoint_release,&
                                              kpoint_sym_type,&
                                              kpoint_type
   USE message_passing,                 ONLY: mp_para_env_type
   USE particle_types,                  ONLY: particle_type
   USE qs_basis_rotation_methods,       ONLY: qs_basis_rotation
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kpoint_operators,             ONLY: kpoint_operator_context_create,&
                                              kpoint_operator_context_release,&
                                              kpoint_operator_context_type,&
                                              kpoint_operator_evaluate_at,&
                                              kpoint_operator_prepare,&
                                              kpoint_operator_work_size,&
                                              kpoint_slot_ks,&
                                              kpoint_slot_s
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_scf_diagonalization,          ONLY: do_general_diag_kp
   USE qs_scf_types,                    ONLY: qs_scf_env_type
   USE scf_control_types,               ONLY: scf_control_type
   USE topology_symmetry,               ONLY: spin_rotation,&
                                              transform_spinor_frame
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_property_symmetry'
   PUBLIC :: property_symmetry_plan, property_symmetry_frame
   PUBLIC :: property_symmetry_mos, property_symmetry_spin_frame
CONTAINS
! **************************************************************************************************
!> \brief Construct a geometric orbit plan; each reconstructed frame still needs H/S validation.
!> \param qs_env Frozen property environment
!> \param k Full mesh in fractional reciprocal coordinates
!> \param grid Gamma-centered property mesh dimensions; absent for an explicit point set
!> \param projection Cartesian periodic-space projector
!> \param backend Existing K290 or SPGLIB backend identifier
!> \param tolerance Geometric symmetry tolerance
!> \param points Owned irreducible k-point object; caller must release it
!> \param representative Source irreducible index for each full point
!> \param operation Source-to-target operation, zero identity or minus one time reversal
! **************************************************************************************************
   SUBROUTINE property_symmetry_plan(qs_env, k, grid, projection, backend, tolerance, points, representative, operation)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(dp), INTENT(IN)                               :: k(:, :)
      INTEGER, INTENT(IN), OPTIONAL                      :: grid(3)
      REAL(dp), INTENT(IN)                               :: projection(3, 3)
      INTEGER, INTENT(IN)                                :: backend
      REAL(dp), INTENT(IN)                               :: tolerance
      TYPE(kpoint_type), POINTER                         :: points
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: representative(:), operation(:)

      INTEGER                                            :: ik, ir, is, nk, nu
      REAL(dp)                                           :: r(3, 3)
      REAL(dp), ALLOCATABLE                              :: unique(:, :)
      TYPE(cell_type), POINTER                           :: cell
      TYPE(kpoint_sym_type), POINTER                     :: sym
      TYPE(particle_type), POINTER                       :: particles(:)

      CALL get_qs_env(qs_env, cell=cell, particle_set=particles)
      nk = SIZE(k, 2)
      CALL kpoint_create(points)
      IF (PRESENT(grid)) THEN
         points%kp_scheme = 'MONKHORST-PACK'
         points%nkp_grid = grid
         points%gamma_centered = .TRUE.
      ELSE
         ! Wilson boundaries can repeat the same point in different reciprocal cells.
         ! Deduplicate only the diagonalization plan, never the directed link graph.
         ALLOCATE (unique(3, nk))
         nu = 0
         DO ik = 1, nk
            DO ir = 1, nu
               IF (kpoint_same_periodic(k(:, ik), unique(:, ir))) EXIT
            END DO
            IF (ir <= nu) CYCLE
            nu = nu + 1
            unique(:, nu) = k(:, ik)
         END DO
         points%kp_scheme = 'GENERAL'
         points%nkp = nu
         ALLOCATE (points%xkp(3, nu), points%wkp(nu))
         points%xkp = unique(:, :nu)
         points%wkp = 1.0_dp/REAL(nu, dp)
      END IF
      points%symmetry = .TRUE.
      points%symmetry_backend = backend
      points%symmetry_reduction_method = use_k290_kpoint_symmetry
      IF (backend == use_spglib_kpoint_backend) points%symmetry_reduction_method = use_spglib_kpoint_symmetry
      points%eps_geo = tolerance
      CALL kpoint_initialize(points, particles, cell, partial_property_orbits=.NOT. PRESENT(grid))
      CALL qs_basis_rotation(qs_env, points)
      ALLOCATE (representative(nk), operation(nk))
      representative = 0
      operation = 0
      DO ik = 1, nk
         DO ir = 1, points%nkp
            IF (kpoint_same_periodic(k(:, ik), points%xkp(:, ir))) THEN
               representative(ik) = ir
               EXIT
            END IF
            IF (kpoint_same_periodic(k(:, ik), -points%xkp(:, ir))) THEN
               representative(ik) = ir
               operation(ik) = -1
               EXIT
            END IF
            sym => points%kp_sym(ir)%kpoint_sym
            IF (.NOT. sym%apply_symmetry) CYCLE
            DO is = 1, sym%nwred
               IF (.NOT. kpoint_same_periodic(k(:, ik), sym%xkp(:, is))) CYCLE
               r = sym%rot(:, :, is)
               IF (MAXVAL(ABS(MATMUL(r, MATMUL(projection, TRANSPOSE(r))) - projection)) > 1.e-9_dp) CYCLE
               representative(ik) = ir
               operation(ik) = is
               EXIT
            END DO
            IF (representative(ik) > 0) EXIT
         END DO
      END DO
      IF (ANY(representative == 0)) THEN
         CPABORT('Property symmetry orbit does not preserve the periodic subspace.')
      END IF
   END SUBROUTINE property_symmetry_plan
! **************************************************************************************************
!> \brief Reconstruct a target eigenframe with atomic image phases and physical SU(2)/Theta action.
!> \param points Property symmetry plan
!> \param representative Source index
!> \param operation Signed operation index
!> \param k Target point in the cell Fourier gauge
!> \param source Source eigenframe
!> \param soc Whether the frame is spin-major doubled
!> \param TARGET Reconstructed frame, not rediagonalized or gauge-aligned
!> \param rotation Cartesian spatial rotation for property tensors
!> \param antiunitary Whether physical time reversal was applied
!> \param status Zero on success; the caller must additionally check the target H/S eigensystem
! **************************************************************************************************
   SUBROUTINE property_symmetry_frame(points, representative, operation, k, source, soc, TARGET, &
                                      rotation, antiunitary, status)
      TYPE(kpoint_type), POINTER                         :: points
      INTEGER, INTENT(IN)                                :: representative, operation
      REAL(dp), INTENT(IN)                               :: k(3)
      COMPLEX(dp), INTENT(IN)                            :: source(:, :)
      LOGICAL, INTENT(IN)                                :: soc
      COMPLEX(dp), INTENT(OUT)                           :: target(:, :)
      REAL(dp), INTENT(OUT)                              :: rotation(3, 3)
      LOGICAL, INTENT(OUT)                               :: antiunitary
      INTEGER, INTENT(OUT)                               :: status

      COMPLEX(dp)                                        :: spin(2, 2)
      COMPLEX(dp), ALLOCATABLE                           :: orbital(:, :)
      INTEGER                                            :: n

      n = SIZE(source, 1)
      IF (soc) n = n/2
      ALLOCATE (orbital(n, n))
      CALL kpoint_cell_gauge_action(points, representative, operation, k, orbital, rotation, antiunitary, status)
      IF (status /= 0) RETURN
      IF (soc) THEN
         CALL spin_rotation(rotation, spin, 1.e-9_dp, status)
         IF (status /= 0) RETURN
         ! CP2K assembles SOC with (sigma_x,-sigma_y,sigma_z), the conjugate spin representation.
         spin = CONJG(spin)
         CALL transform_spinor_frame(orbital, spin, source, antiunitary, TARGET, 1.e-9_dp, status)
      ELSE IF (antiunitary) THEN
         TARGET = MATMUL(orbital, CONJG(source))
      ELSE
         TARGET = MATMUL(orbital, source)
      END IF
   END SUBROUTINE property_symmetry_frame

! **************************************************************************************************
!> \brief Diagonalize independent property points and strictly reconstruct the full scalar frame.
!> \param qs_env Frozen k-point environment
!> \param TARGET Initialized full explicit k-point MO set
!> \param points Geometric orbit plan, initialized here with representative MOs
!> \param representative Source representative for each target
!> \param operation Accepted atom/AO operation for each target
!> \param added_mos Same scalar band window as the full property calculation
!> \param soc Check that the scalar second-variational window does not cut a degeneracy
!> \param tolerance Target metric and scaled Hamiltonian residual tolerance
!> \param max_memory_mb Bound on additional replicated array storage per rank
!> \param iw Output unit
!> \note No Ritz rotation, subspace expansion or target diagonalization can hide a broken orbit.
! **************************************************************************************************
   SUBROUTINE property_symmetry_mos(qs_env, TARGET, points, representative, operation, added_mos, soc, &
                                    tolerance, max_memory_mb, iw)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(kpoint_type), POINTER                         :: target, points
      INTEGER, INTENT(IN)                                :: representative(:)
      INTEGER, INTENT(INOUT)                             :: operation(:)
      INTEGER, INTENT(IN)                                :: added_mos
      LOGICAL, INTENT(IN)                                :: soc
      REAL(dp), INTENT(IN)                               :: tolerance
      INTEGER, INTENT(IN)                                :: max_memory_mb, iw

      COMPLEX(dp), ALLOCATABLE                           :: frame(:, :), h(:, :), metric(:, :), &
                                                            residual(:, :), s(:, :), source(:, :)
      INTEGER                                            :: ik, ilocal, ir, nao, ni, nmo, op, &
                                                            owners, spin
      INTEGER(int_8)                                     :: bytes
      LOGICAL                                            :: diis, prepared, valid
      REAL(dp)                                           :: errors(3), window_gap, worst(3)
      REAL(dp), ALLOCATABLE                              :: energies(:), hs(:, :, :, :), occ(:), &
                                                            ss(:, :, :, :)
      REAL(dp), POINTER                                  :: occupations(:), values(:)
      TYPE(cp_blacs_env_type), POINTER                   :: blacs
      TYPE(dbcsr_p_type), POINTER                        :: hr(:, :), sr(:, :)
      TYPE(kpoint_type), POINTER                         :: images, scf_points
      TYPE(mo_set_type), POINTER                         :: mos(:)
      TYPE(mp_para_env_type), POINTER                    :: para
      TYPE(neighbor_list_set_p_type), POINTER            :: all_pairs(:), pairs(:)
      TYPE(qs_scf_env_type), POINTER                     :: scf
      TYPE(scf_control_type), POINTER                    :: control

      CALL get_qs_env(qs_env, para_env=para, blacs_env=blacs, mos=mos, sab_orb=pairs, &
                      matrix_ks_kp=hr, matrix_s_kp=sr, scf_env=scf, scf_control=control, &
                      sab_all=all_pairs, kpoints=scf_points)
      CPASSERT(ALL(TARGET%kp_range == [1, TARGET%nkp]))
      CALL get_mo_set(TARGET%kp_env(1)%kpoint_env%mos(1), nao=nao, nmo=nmo)
      points%use_real_wfn = .FALSE.
      CALL kpoint_env_initialize(points, para, blacs)
      CALL kpoint_initialize_mos(points, mos, added_mos + MERGE(1, 0, soc .AND. nmo < nao))
      CALL kpoint_initialize_mo_set(points)
      CALL kpoint_init_cell_index(points, pairs, para, ni)
      CALL do_general_diag_kp(hr, sr, points, scf, control, .FALSE., diis)
      CALL property_symmetry_groups(TARGET, points, representative, operation, hr, sr, soc, &
                                    tolerance, max_memory_mb, worst, prepared, iw)
      IF (.NOT. prepared) THEN
         IF (iw > 0) WRITE (iw, '(T2,A)') &
            'TOPOLOGY| Scalar H/S snapshots unavailable; using collective frame validation.'
         ! Symmetric DBCSR storage omits mirror blocks; dense Fourier arrays need both images.
         NULLIFY (images)
         CALL kpoint_create(images)
         CALL kpoint_init_cell_index(images, all_pairs, para, ni)
         bytes = 8_int_8*nao*nao*ni*(SIZE(hr, 1) + 1_int_8) + &
                 16_int_8*(4_int_8*nao*nao + 8_int_8*nao*nmo + 3_int_8*nmo*nmo)
         IF (max_memory_mb < 1 .OR. bytes > INT(max_memory_mb, int_8)*1024*1024) THEN
            CPABORT('Property symmetry arrays exceed SYMMETRY_MAX_MEMORY_MB per rank.')
         END IF
         ALLOCATE (hs(SIZE(hr, 1), nao, nao, ni), ss(1, nao, nao, ni), SOURCE=0.0_dp)
         CALL replicate_rs_matrices(hr, scf_points, hs, images%cell_to_index)
         CALL replicate_rs_matrices(sr, scf_points, ss, images%cell_to_index)
         ALLOCATE (source(nao, nmo), frame(nao, nmo), &
                   h(nao, nao), s(nao, nao), metric(nmo, nmo), residual(nao, nmo), energies(nmo), occ(nmo))
         worst = 0.0_dp
         DO ir = 1, points%nkp
            DO spin = 1, SIZE(hr, 1)
               source = CMPLX(0.0_dp, 0.0_dp, dp)
               energies = 0.0_dp
               occ = 0.0_dp
               window_gap = 0.0_dp
               owners = 0
               IF (ir >= points%kp_range(1) .AND. ir <= points%kp_range(2)) THEN
                  ilocal = ir - points%kp_range(1) + 1
                  CALL cp_cfm_get_submatrix(points%kp_env(ilocal)%kpoint_env%mos(spin)%cmo_coeff, source)
                  CALL get_mo_set(points%kp_env(ilocal)%kpoint_env%mos(spin), &
                                  eigenvalues=values, occupation_numbers=occupations)
                  energies = values(:nmo)
                  occ = occupations(:nmo)
                  IF (soc .AND. nmo < nao) window_gap = values(nmo + 1) - values(nmo)
                  owners = 1
               END IF
               ! Each inter-group communicator contains one rank from every k-point group.
               ! The owning group has already replicated its frame within that group.
               CALL points%para_env_inter_kp%sum(owners)
               CPASSERT(owners == 1)
               CALL points%para_env_inter_kp%sum(source)
               CALL points%para_env_inter_kp%sum(energies)
               CALL points%para_env_inter_kp%sum(occ)
               CALL points%para_env_inter_kp%sum(window_gap)
               IF (soc .AND. nmo < nao) THEN
                  IF (window_gap <= 1.e-8_dp) THEN
                     CALL cp_abort(__LOCATION__, 'Property SOC scalar window cuts a degeneracy; '// &
                                   'increase WANNIER90 ADDED_MOS.')
                  END IF
               END IF
               DO ik = 1, TARGET%nkp
                  IF (representative(ik) /= ir) CYCLE
                  CALL rs_to_kp(hs(spin, :, :, :), h, images%index_to_cell, TARGET%xkp(:, ik))
                  CALL rs_to_kp(ss(1, :, :, :), s, images%index_to_cell, TARGET%xkp(:, ik))
                  op = operation(ik)
                  CALL property_symmetry_check(points, ir, TARGET%xkp(:, ik), source, energies, h, s, tolerance, &
                                               op, frame, metric, residual, errors, valid)
                  IF (.NOT. valid) THEN
                     IF (iw > 0) WRITE (iw, '(T2,A,I0,A,3ES14.5)') &
                        'TOPOLOGY| Invalid property orbit at point ', ik, ' metric/H/projected-H residuals: ', errors
                     CPABORT('Property symmetry does not preserve the frozen H/S eigensystem.')
                  END IF
                  operation(ik) = op
                  worst = MAX(worst, errors)
                  CALL cp_cfm_set_submatrix(TARGET%kp_env(ik)%kpoint_env%mos(spin)%cmo_coeff, frame)
                  CALL get_mo_set(TARGET%kp_env(ik)%kpoint_env%mos(spin), &
                                  eigenvalues=values, occupation_numbers=occupations)
                  IF (ASSOCIATED(values)) values(:nmo) = energies
                  IF (ASSOCIATED(occupations)) occupations(:nmo) = occ
               END DO
            END DO
         END DO
         CALL kpoint_release(images)
      END IF
      IF (iw > 0) THEN
         WRITE (iw, '(T2,A,I0,A,I0)') 'TOPOLOGY| Property scalar diagonalizations: ', points%nkp, ' / ', TARGET%nkp
         WRITE (iw, '(T2,A,I0,A,I0)') 'TOPOLOGY| Property scalar MPI groups: ', points%nkp_groups, &
            ' ranks per group: ', points%para_env_kp%num_pe
         WRITE (iw, '(T2,A,I0,A,I0)') &
            'TOPOLOGY| Property scalar validation MPI groups: ', MERGE(points%nkp_groups, 1, prepared), &
            ' ranks per group: ', MERGE(points%para_env_kp%num_pe, para%num_pe, prepared)
         WRITE (iw, '(T2,A,3ES14.5)') 'TOPOLOGY| Property scalar metric/H/projected-H residuals: ', worst
      END IF
   END SUBROUTINE property_symmetry_mos

! **************************************************************************************************
!> \brief Validate one scalar image without rotating, extending or rediagonalizing its subspace.
!> \param points Geometric orbit plan
!> \param ir Representative index
!> \param k Fractional target coordinates
!> \param source Representative coefficients
!> \param energies Representative eigenvalues
!> \param h Target Hamiltonian
!> \param s Target overlap
!> \param tolerance Common metric and Hamiltonian tolerance
!> \param op Preferred operation on input; accepted operation on success
!> \param frame Reconstructed coefficients
!> \param metric Band-space scratch
!> \param residual AO-space scratch
!> \param errors Metric, AO-scaled and projected Hamiltonian residuals
!> \param valid Whether a geometrically allowed operation passes every physical check
! **************************************************************************************************
   SUBROUTINE property_symmetry_check(points, ir, k, source, energies, h, s, tolerance, &
                                      op, frame, metric, residual, errors, valid)
      TYPE(kpoint_type), POINTER                         :: points
      INTEGER, INTENT(IN)                                :: ir
      REAL(dp), INTENT(IN)                               :: k(3)
      COMPLEX(dp), INTENT(IN)                            :: source(:, :)
      REAL(dp), INTENT(IN)                               :: energies(:)
      COMPLEX(dp), INTENT(IN)                            :: h(:, :), s(:, :)
      REAL(dp), INTENT(IN)                               :: tolerance
      INTEGER, INTENT(INOUT)                             :: op
      COMPLEX(dp), INTENT(OUT)                           :: frame(:, :), metric(:, :), residual(:, :)
      REAL(dp), INTENT(OUT)                              :: errors(3)
      LOGICAL, INTENT(OUT)                               :: valid

      INTEGER                                            :: j, status, trial, try_op
      LOGICAL                                            :: antiunitary
      REAL(dp)                                           :: rotation(3, 3)

      valid = .FALSE.
      errors = HUGE(1.0_dp)
      DO trial = 0, points%kp_sym(ir)%kpoint_sym%nwred
         try_op = op
         IF (trial > 0) try_op = trial
         CALL property_symmetry_frame(points, ir, try_op, k, source, .FALSE., frame, rotation, antiunitary, status)
         IF (status /= 0) CYCLE
         residual = MATMUL(s, frame)
         metric = MATMUL(CONJG(TRANSPOSE(frame)), residual)
         DO j = 1, SIZE(energies)
            metric(j, j) = metric(j, j) - 1.0_dp
            residual(:, j) = residual(:, j)*energies(j)
         END DO
         residual = MATMUL(h, frame) - residual
         errors(1:2) = [MAXVAL(ABS(metric)), MAXVAL(ABS(residual))/ &
                        MAX(1.0_dp, MAXVAL(ABS(h))*MAXVAL(ABS(frame))*REAL(SIZE(h, 1), dp))]
         ! The projected residual also detects energy errors amplified by AO conditioning.
         metric = MATMUL(CONJG(TRANSPOSE(frame)), residual)
         errors(3) = MAXVAL(ABS(metric))/MAX(1.0_dp, MAXVAL(ABS(energies)))
         IF (ALL(ieee_is_finite(errors)) .AND. MAXVAL(errors) <= tolerance) THEN
            op = try_op
            valid = .TRUE.
            EXIT
         END IF
      END DO
   END SUBROUTINE property_symmetry_check

! **************************************************************************************************
!> \brief Check target frames independently, then transfer one validated frame per group and round.
!> \param TARGET Full-graph MOs on the global communicator
!> \param points Initialized representative MOs on k-point groups
!> \param representative Source representative for each target
!> \param operation Preferred and accepted atom/AO operations
!> \param hr Frozen real-space Hamiltonian
!> \param sr Frozen real-space overlap
!> \param soc Check the scalar SOC window boundary
!> \param tolerance Metric and Hamiltonian residual tolerance
!> \param max_memory_mb Additional replicated scratch bound
!> \param worst Largest residuals over every target and group
!> \param prepared Whether bounded H/S snapshots were available
!> \param iw Output unit
! **************************************************************************************************
   SUBROUTINE property_symmetry_groups(TARGET, points, representative, operation, hr, sr, soc, &
                                       tolerance, max_memory_mb, worst, prepared, iw)
      TYPE(kpoint_type), POINTER                         :: target, points
      INTEGER, INTENT(IN)                                :: representative(:)
      INTEGER, INTENT(INOUT)                             :: operation(:)
      TYPE(dbcsr_p_type), POINTER                        :: hr(:, :), sr(:, :)
      LOGICAL, INTENT(IN)                                :: soc
      REAL(dp), INTENT(IN)                               :: tolerance
      INTEGER, INTENT(IN)                                :: max_memory_mb
      REAL(dp), INTENT(OUT)                              :: worst(3)
      LOGICAL, INTENT(OUT)                               :: prepared
      INTEGER, INTENT(IN)                                :: iw

      CHARACTER(LEN=*), PARAMETER :: routineN = 'property_symmetry_groups'

      COMPLEX(dp), ALLOCATABLE                           :: frame(:, :), h(:, :), metric(:, :), &
                                                            phase(:, :), residual(:, :), s(:, :), &
                                                            source(:, :)
      INTEGER                                            :: group, handle, ik, ilocal, ir, last_ir, &
                                                            mygroup, nao, nentry, nmo, op, round, &
                                                            shared_op, shared_status, spin, status
      INTEGER(int_8)                                     :: bytes
      INTEGER, ALLOCATABLE                               :: bounds(:, :), count_rep(:), first(:), &
                                                            seen(:), work_order(:)
      LOGICAL                                            :: ready_h, ready_s, valid
      REAL(dp)                                           :: errors(3), shared_errors(3), window_gap
      REAL(dp), ALLOCATABLE                              :: energies(:), occ(:), shared_energies(:), &
                                                            shared_occ(:)
      REAL(dp), POINTER                                  :: occupations(:), values(:)
      TYPE(cp_cfm_type)                                  :: ao
      TYPE(cp_cfm_type), POINTER                         :: from
      TYPE(kpoint_operator_context_type)                 :: context
      TYPE(cp_fm_type)                                   :: global_ao, local_ao
      TYPE(cp_fm_struct_type), POINTER                   :: global_struct, local_struct, work_struct
      TYPE(cp_cfm_type), TARGET                          :: buffer, empty
      TYPE(mp_para_env_type), POINTER                    :: inter, para

      CALL timeset(routineN, handle)
      para => points%blacs_env_all%para_env
      inter => points%para_env_inter_kp
      mygroup = inter%mepos + 1
      CALL get_mo_set(TARGET%kp_env(1)%kpoint_env%mos(1), nao=nao, nmo=nmo)
      bytes = 16_int_8*(6_int_8*nao*nao + 8_int_8*nao*nmo + 3_int_8*nmo*nmo)
      IF (max_memory_mb < 1 .OR. bytes > INT(max_memory_mb, int_8)*1024*1024) THEN
         CPABORT('Property symmetry arrays exceed SYMMETRY_MAX_MEMORY_MB per rank.')
      END IF
      NULLIFY (global_struct, local_struct, work_struct)
      CALL cp_fm_struct_create(global_struct, nrow_global=nao, ncol_global=nao, &
                               para_env=para, context=points%blacs_env_all)
      CALL cp_fm_struct_create(local_struct, nrow_global=nao, ncol_global=nao, &
                               para_env=points%para_env_kp, context=points%blacs_env)
      CALL cp_fm_create(global_ao, global_struct)
      CALL cp_fm_create(local_ao, local_struct)
      CALL kpoint_operator_context_create(context, points, global_ao, SIZE(hr, 1))
      CALL kpoint_operator_prepare(context, hr, kpoint_slot_ks, ready_h, TARGET=local_ao)
      CALL kpoint_operator_prepare(context, sr, kpoint_slot_s, ready_s, matrix_row=1, TARGET=local_ao)
      prepared = ready_h .AND. ready_s
      worst = 0.0_dp
      IF (prepared) THEN
         CALL cp_cfm_create(ao, local_struct)
         CALL cp_fm_struct_create(work_struct, nrow_global=nao, ncol_global=nmo, &
                                  para_env=points%para_env_kp, context=points%blacs_env)
         CALL cp_cfm_create(buffer, work_struct)
         CALL kpoint_operator_work_size(context, nentry)
         ALLOCATE (phase(1, nentry), source(nao, nmo), frame(nao, nmo), &
                   h(nao, nao), s(nao, nao), metric(nmo, nmo), residual(nao, nmo), energies(nmo), &
                   occ(nmo), shared_energies(nmo), shared_occ(nmo))
         ! Stable representative buckets keep source-frame storage independent of orbit size.
         ALLOCATE (count_rep(points%nkp), first(0:points%nkp), work_order(TARGET%nkp), &
                   bounds(2, points%nkp_groups), seen(TARGET%nkp))
         count_rep = 0
         DO ik = 1, TARGET%nkp
            ir = representative(ik)
            count_rep(ir) = count_rep(ir) + 1
         END DO
         first(0) = 0
         DO ir = 1, points%nkp
            first(ir) = first(ir - 1) + count_rep(ir)
         END DO
         count_rep = 0
         DO ik = 1, TARGET%nkp
            ir = representative(ik)
            count_rep(ir) = count_rep(ir) + 1
            work_order(first(ir - 1) + count_rep(ir)) = ik
         END DO
         DO group = 1, points%nkp_groups
            bounds(:, group) = [first(points%kp_dist(1, group) - 1), first(points%kp_dist(2, group))]
         END DO
         DO spin = 1, SIZE(hr, 1)
            last_ir = 0
            seen = 0
            DO round = 1, MAXVAL(bounds(2, :) - bounds(1, :))
               errors = 0.0_dp
               status = 0
               op = 0
               energies = 0.0_dp
               occ = 0.0_dp
               IF (bounds(1, mygroup) + round <= bounds(2, mygroup)) THEN
                  ik = work_order(bounds(1, mygroup) + round)
                  ir = representative(ik)
                  ilocal = ir - points%kp_range(1) + 1
                  IF (ir /= last_ir) THEN
                     CALL cp_cfm_get_submatrix(points%kp_env(ilocal)%kpoint_env%mos(spin)%cmo_coeff, source)
                     last_ir = ir
                  END IF
                  CALL get_mo_set(points%kp_env(ilocal)%kpoint_env%mos(spin), &
                                  eigenvalues=values, occupation_numbers=occupations)
                  energies = values(:nmo)
                  occ = occupations(:nmo)
                  window_gap = HUGE(1.0_dp)
                  IF (soc .AND. nmo < nao) window_gap = values(nmo + 1) - values(nmo)
                  IF (window_gap <= 1.e-8_dp) THEN
                     status = 2
                  ELSE
                     CALL kpoint_operator_evaluate_at(context, spin, kpoint_slot_ks, TARGET%xkp(:, ik), phase, TARGET=ao)
                     CALL cp_cfm_get_submatrix(ao, h)
                     CALL kpoint_operator_evaluate_at(context, spin, kpoint_slot_s, TARGET%xkp(:, ik), phase, TARGET=ao)
                     CALL cp_cfm_get_submatrix(ao, s)
                     op = operation(ik)
                     CALL property_symmetry_check(points, ir, TARGET%xkp(:, ik), source, energies, h, s, tolerance, &
                                                  op, frame, metric, residual, errors, valid)
                     status = MERGE(0, 1, valid)
                     IF (valid) THEN
                        CALL cp_cfm_set_submatrix(buffer, frame)
                     END IF
                  END IF
               END IF
               ! Every group joins the same transfer sequence, including groups idle in this round.
               DO group = 1, points%nkp_groups
                  IF (bounds(1, group) + round > bounds(2, group)) CYCLE
                  ik = work_order(bounds(1, group) + round)
                  shared_status = status
                  shared_op = op
                  shared_errors = errors
                  shared_energies = energies
                  shared_occ = occ
                  CALL inter%bcast(shared_status, group - 1)
                  IF (shared_status == 2) THEN
                     CALL cp_abort(__LOCATION__, 'Property SOC scalar window cuts a degeneracy; '// &
                                   'increase WANNIER90 ADDED_MOS.')
                  END IF
                  CALL inter%bcast(shared_errors, group - 1)
                  IF (shared_status /= 0) THEN
                     IF (iw > 0) WRITE (iw, '(T2,A,I0,A,3ES14.5)') &
                        'TOPOLOGY| Invalid property orbit at point ', ik, &
                        ' metric/H/projected-H residuals: ', shared_errors
                     CPABORT('Property symmetry does not preserve the frozen H/S eigensystem.')
                  END IF
                  CALL inter%bcast(shared_op, group - 1)
                  CALL inter%bcast(shared_energies, group - 1)
                  CALL inter%bcast(shared_occ, group - 1)
                  operation(ik) = shared_op
                  worst = MAX(worst, shared_errors)
                  seen(ik) = seen(ik) + 1
                  from => empty
                  IF (mygroup == group) from => buffer
                  CALL cp_cfm_copy_general(from, TARGET%kp_env(ik)%kpoint_env%mos(spin)%cmo_coeff, para)
                  CALL get_mo_set(TARGET%kp_env(ik)%kpoint_env%mos(spin), &
                                  eigenvalues=values, occupation_numbers=occupations)
                  IF (ASSOCIATED(values)) values(:nmo) = shared_energies
                  IF (ASSOCIATED(occupations)) occupations(:nmo) = shared_occ
               END DO
            END DO
            CPASSERT(ALL(seen == 1))
         END DO
         CALL para%max(worst)
         CALL cp_cfm_release(buffer)
         CALL cp_fm_struct_release(work_struct)
         CALL cp_cfm_release(ao)
      END IF
      CALL kpoint_operator_context_release(context)
      CALL cp_fm_release(global_ao)
      CALL cp_fm_release(local_ao)
      CALL cp_fm_struct_release(global_struct)
      CALL cp_fm_struct_release(local_struct)
      CALL timestop(handle)
   END SUBROUTINE property_symmetry_groups

! **************************************************************************************************
!> \brief Spin action in a scalar MO gauge reconstructed with the same orbit operation.
!> \param points Geometric orbit plan
!> \param representative Source representative
!> \param operation Accepted scalar atom/AO operation
!> \param source Source second-variational spinor frame
!> \param TARGET Target second-variational frame; caller must validate its Hamiltonian
!> \param status Zero on success
! **************************************************************************************************
   SUBROUTINE property_symmetry_spin_frame(points, representative, operation, source, TARGET, status)
      TYPE(kpoint_type), POINTER                         :: points
      INTEGER, INTENT(IN)                                :: representative, operation
      COMPLEX(dp), INTENT(IN)                            :: source(:, :)
      COMPLEX(dp), INTENT(OUT)                           :: target(:, :)
      INTEGER, INTENT(OUT)                               :: status

      COMPLEX(dp)                                        :: spin(2, 2)
      COMPLEX(dp), ALLOCATABLE                           :: identity(:, :)
      INTEGER                                            :: i, n
      LOGICAL                                            :: antiunitary
      REAL(dp)                                           :: rotation(3, 3)

      n = SIZE(source, 1)/2
      ALLOCATE (identity(n, n))
      identity = 0.0_dp
      DO i = 1, n
         identity(i, i) = 1.0_dp
      END DO
      rotation = 0.0_dp
      DO i = 1, 3
         rotation(i, i) = 1.0_dp
      END DO
      antiunitary = operation == -1
      IF (operation > 0) THEN
         rotation = points%kp_sym(representative)%kpoint_sym%rot(:, :, operation)
         antiunitary = points%kp_sym(representative)%kpoint_sym%rotp(operation) < 0
      END IF
      CALL spin_rotation(rotation, spin, 1.e-9_dp, status)
      IF (status /= 0) RETURN
      CALL transform_spinor_frame(identity, CONJG(spin), source, antiunitary, TARGET, 1.e-9_dp, status)
   END SUBROUTINE property_symmetry_spin_frame
END MODULE qs_property_symmetry
