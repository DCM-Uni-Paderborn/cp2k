!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
!> \brief Native orbitals and densities for fixed-occupation Skala constraint steps.
MODULE qs_skala_orbitals
   USE cp_cfm_types,                    ONLY: cp_cfm_get_submatrix,&
                                              cp_cfm_set_submatrix
   USE cp_dbcsr_api,                    ONLY: dbcsr_copy,&
                                              dbcsr_release,&
                                              dbcsr_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr
   USE cp_fm_types,                     ONLY: cp_fm_get_submatrix,&
                                              cp_fm_set_submatrix
   USE kinds,                           ONLY: dp
   USE kpoint_mo_symmetry_methods,      ONLY: kpoint_cell_gauge_action,&
                                              kpoint_same_periodic
   USE kpoint_types,                    ONLY: kpoint_env_type,&
                                              kpoint_sym_type,&
                                              kpoint_type
   USE qs_mo_types,                     ONLY: mo_set_type
   USE qs_upf_inverse_grid,             ONLY: qs_upf_inverse_matrix_scatter
   USE qs_upf_inverse_types,            ONLY: qs_upf_inverse_type
   USE skala_orbital_constraints,       ONLY: &
        skala_orbital_chart_create, skala_orbital_chart_derivative, skala_orbital_chart_evaluate, &
        skala_orbital_chart_release, skala_orbital_chart_restrict, &
        skala_orbital_chart_set_occupations, skala_orbital_chart_set_reference, &
        skala_orbital_chart_size, skala_orbital_chart_symmetry, skala_orbital_chart_type, &
        skala_orbital_multiplets
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: qs_skala_orbitals_type, qs_skala_orbitals_create, qs_skala_orbitals_release
   PUBLIC :: qs_skala_orbitals_adapt_symmetry
   PUBLIC :: qs_skala_orbitals_evaluate, qs_skala_orbitals_canonicalize
   PUBLIC :: qs_skala_orbitals_install_gamma, qs_skala_orbitals_install_kpoint

   TYPE qs_skala_orbitals_type
      TYPE(skala_orbital_chart_type) :: chart
      TYPE(mo_set_type), POINTER :: native => NULL()
      COMPLEX(dp), ALLOCATABLE :: coefficients(:, :)
      REAL(dp), ALLOCATABLE :: occupations(:), eigenvalues(:)
      INTEGER, ALLOCATABLE :: occupation_groups(:)
      REAL(dp) :: weight = 1.0_dp
      LOGICAL :: canonicalized = .FALSE., variable_occupations = .FALSE.
   END TYPE qs_skala_orbitals_type

CONTAINS
! **************************************************************************************************
!> \brief Preserve the complete little group of the initial original-metric orbital reference.
!>        Fixed occupations use the full unitary/antiunitary group without diagonalizing its actions.
!>        Thermal variables share occupations within irreducible multiplets. A real-character
!>        adaptation remains available inside equal native occupation and energy subspaces.
!>        A failed adaptation leaves the state and all native MO/density storage unchanged.
!> \param state Initial normalized orbital state, before canonical acceptance.
!> \param metric Original complete charge metric.
!> \param kpoints Native cell-gauge symmetry metadata.
!> \param representative Original irreducible point.
!> \param info Zero success, minus one incompatible reference, minus two unsupported representation.
! **************************************************************************************************
   SUBROUTINE qs_skala_orbitals_adapt_symmetry(state, metric, kpoints, representative, info)
      TYPE(qs_skala_orbitals_type), INTENT(INOUT)        :: state
      COMPLEX(dp), INTENT(IN)                            :: metric(:, :)
      TYPE(kpoint_type), POINTER                         :: kpoints
      INTEGER, INTENT(IN)                                :: representative
      INTEGER, INTENT(OUT)                               :: info

      INTEGER                                            :: op
      INTEGER, ALLOCATABLE                               :: operations(:)
      TYPE(kpoint_sym_type), POINTER                     :: sym

      ! Preserve the smaller independent chart when real characters suffice.
      CALL adapt_real_characters(state, metric, kpoints, representative, info)
      IF (info == 0) RETURN
      IF (.NOT. ASSOCIATED(kpoints)) RETURN
      IF (.NOT. ASSOCIATED(kpoints%kp_sym)) RETURN
      sym => kpoints%kp_sym(representative)%kpoint_sym
      IF (.NOT. ASSOCIATED(sym)) RETURN
      IF (.NOT. ASSOCIATED(sym%rotp)) RETURN
      operations = PACK([(op, op=1, sym%nwred)], &
                        [(kpoint_same_periodic(sym%xkp(:, op), kpoints%xkp(:, representative)), op=1, sym%nwred)])
      IF (SIZE(operations) == 0) RETURN
      CALL adapt_little_group(state, metric, kpoints, representative, operations, info)
   END SUBROUTINE qs_skala_orbitals_adapt_symmetry

! **************************************************************************************************
!> \brief Smaller independent chart for compatible real one-dimensional characters.
!> \param state Initial complete orbital state.
!> \param metric Original charge metric.
!> \param kpoints Native cell-gauge symmetry metadata.
!> \param representative Irreducible k-point.
!> \param info Zero on success, minus one incompatible reference, minus two unsupported representation.
! **************************************************************************************************
   SUBROUTINE adapt_real_characters(state, metric, kpoints, representative, info)
      TYPE(qs_skala_orbitals_type), INTENT(INOUT)        :: state
      COMPLEX(dp), INTENT(IN)                            :: metric(:, :)
      TYPE(kpoint_type), POINTER                         :: kpoints
      INTEGER, INTENT(IN)                                :: representative
      INTEGER, INTENT(OUT)                               :: info

      REAL(dp), PARAMETER                                :: tolerance = 3.e-11_dp

      COMPLEX(dp), ALLOCATABLE :: a(:, :), adapted(:, :), b(:, :), block(:, :), c(:, :), &
         final(:, :), gram(:, :), label(:, :), original(:, :), q(:, :), source(:, :), work(:)
      INTEGER                                            :: anti_operation, i, j, m, n, nop, op, &
                                                            row, status, unitary_count
      INTEGER, ALLOCATABLE                               :: groups(:), indices(:), operations(:)
      LOGICAL                                            :: anti, has_anti, valid
      LOGICAL, ALLOCATABLE                               :: done(:), sector(:)
      REAL(dp)                                           :: rotation(3, 3)
      REAL(dp), ALLOCATABLE                              :: real_basis(:, :), real_work(:), &
                                                            values(:), zero(:)
      TYPE(kpoint_sym_type), POINTER                     :: sym
      TYPE(skala_orbital_chart_type)                     :: chart

      info = -2
      CPASSERT(.NOT. state%canonicalized)
      CPASSERT(ALLOCATED(state%coefficients))
      n = SIZE(state%coefficients, 2)
      CPASSERT(ALL(SHAPE(state%coefficients) == [n, n]))
      CPASSERT(ALL(SHAPE(metric) == [n, n]))
      IF (.NOT. ASSOCIATED(kpoints)) RETURN
      IF (.NOT. ASSOCIATED(kpoints%kp_sym)) RETURN
      sym => kpoints%kp_sym(representative)%kpoint_sym
      IF (.NOT. ASSOCIATED(sym)) RETURN
      IF (.NOT. ASSOCIATED(sym%rotp)) THEN
         info = 0
         RETURN
      END IF
      operations = PACK([(op, op=1, sym%nwred)], &
                        [(kpoint_same_periodic(sym%xkp(:, op), kpoints%xkp(:, representative)), op=1, sym%nwred)])
      nop = SIZE(operations)
      IF (nop == 0) THEN
         info = 0
         RETURN
      END IF
      IF (nop > 16) RETURN
      has_anti = .FALSE.
      anti_operation = 0
      DO i = 1, nop
         op = operations(i)
         IF (sym%rotp(op) >= 0) CYCLE
         has_anti = .TRUE.
         anti_operation = op
         EXIT
      END DO
      IF (has_anti .AND. kpoint_same_periodic(kpoints%xkp(:, representative), -kpoints%xkp(:, representative))) THEN
         anti_operation = -1
      END IF
      ALLOCATE (a(n, n), done(n), sector(n), groups(n))
      adapted = state%coefficients
      CALL orbital_density(state, original)
      done = .FALSE.
      DO i = 1, n
         IF (done(i)) CYCLE
         sector = state%occupations == state%occupations(i)
         IF (state%variable_occupations) THEN
            ! Keep retained thermal bands distinct from the permanently empty completion.
            DO j = 1, n
               IF ((i <= state%native%nmo) .NEQV. (j <= state%native%nmo)) sector(j) = .FALSE.
               IF (i <= state%native%nmo .AND. j <= state%native%nmo) THEN
                  sector(j) = sector(j) .AND. state%native%eigenvalues(j) == state%native%eigenvalues(i)
               END IF
            END DO
         END IF
         indices = PACK([(j, j=1, n)], sector)
         done(indices) = .TRUE.
         m = SIZE(indices)
         c = state%coefficients(:, indices)
         IF (has_anti) THEN
            CALL kpoint_cell_gauge_action(kpoints, representative, anti_operation, &
                                          kpoints%xkp(:, representative), a, rotation, anti, status)
            IF (status /= 0 .OR. .NOT. anti) RETURN
            source = MATMUL(a, CONJG(c))
            b = MATMUL(CONJG(TRANSPOSE(state%coefficients)), MATMUL(metric, source))
            info = -1
            DO j = 1, n
               IF (sector(j)) CYCLE
               IF (MAXVAL(ABS(b(j, :))) > tolerance) RETURN
            END DO
            block = b(indices, :)
            info = -2
            IF (MAXVAL(ABS(block - TRANSPOSE(block))) > tolerance) RETURN
            ALLOCATE (real_basis(2*m, 2*m), values(2*m), real_work(6*m))
            real_basis(:m, :m) = REAL(block, dp)
            real_basis(:m, m + 1:) = AIMAG(block)
            real_basis(m + 1:, :m) = AIMAG(block)
            real_basis(m + 1:, m + 1:) = -REAL(block, dp)
            CALL dsyev('V', 'U', 2*m, real_basis, 2*m, values, real_work, SIZE(real_work), status)
            IF (status /= 0 .OR. COUNT(values > 0.5_dp) /= m) RETURN
            q = CMPLX(real_basis(:m, m + 1:), real_basis(m + 1:, m + 1:), dp)
            gram = MATMUL(CONJG(TRANSPOSE(q)), q)
            CALL zpotrf('L', m, gram, m, status)
            IF (status /= 0) RETURN
            CALL ztrsm('R', 'L', 'C', 'N', m, m, CMPLX(1.0_dp, 0.0_dp, dp), gram, m, q, m)
            c = MATMUL(c, q)
            IF (MAXVAL(ABS(MATMUL(a, CONJG(c)) - c)) > tolerance) RETURN
            DEALLOCATE (real_basis, values, real_work)
         END IF
         ALLOCATE (label(m, m)); label = 0.0_dp
         unitary_count = 0
         DO j = 1, nop
            op = operations(j)
            IF (sym%rotp(op) < 0) CYCLE
            unitary_count = unitary_count + 1
            IF (unitary_count > 8) RETURN
            CALL kpoint_cell_gauge_action(kpoints, representative, op, kpoints%xkp(:, representative), &
                                          a, rotation, anti, status)
            IF (status /= 0 .OR. anti) RETURN
            source = MATMUL(a, c)
            b = MATMUL(CONJG(TRANSPOSE(state%coefficients)), MATMUL(metric, source))
            info = -1
            DO row = 1, n
               IF (sector(row)) CYCLE
               IF (MAXVAL(ABS(b(row, :))) > tolerance) RETURN
            END DO
            block = MATMUL(CONJG(TRANSPOSE(c)), MATMUL(metric, source))
            info = -2
            IF (MAXVAL(ABS(block - CONJG(TRANSPOSE(block)))) > tolerance) RETURN
            IF (has_anti) THEN
               IF (MAXVAL(ABS(AIMAG(block))) > tolerance) RETURN
            END IF
            label = label + REAL(2**(unitary_count - 1), dp)*block
         END DO
         ALLOCATE (values(m), work(2*m), real_work(3*m))
         IF (has_anti) THEN
            real_basis = REAL(label, dp)
            CALL dsyev('V', 'U', m, real_basis, m, values, real_work, SIZE(real_work), status)
            label = CMPLX(real_basis, 0.0_dp, dp)
            DEALLOCATE (real_basis)
         ELSE
            CALL zheev('V', 'U', m, label, m, values, work, SIZE(work), real_work, status)
         END IF
         IF (status /= 0) RETURN
         adapted(:, indices) = MATMUL(c, label)
         DEALLOCATE (label, values, work, real_work, indices)
      END DO
      groups = 1
      DO j = 1, nop
         op = operations(j)
         CALL kpoint_cell_gauge_action(kpoints, representative, op, kpoints%xkp(:, representative), &
                                       a, rotation, anti, status)
         IF (status /= 0) RETURN
         source = adapted
         IF (anti) source = CONJG(source)
         b = MATMUL(CONJG(TRANSPOSE(adapted)), MATMUL(metric, MATMUL(a, source)))
         DO i = 1, n
            IF (ABS(AIMAG(b(i, i))) > tolerance .OR. ABS(ABS(REAL(b(i, i), dp)) - 1.0_dp) > tolerance) RETURN
            IF (REAL(b(i, i), dp) < 0.0_dp) groups(i) = groups(i) + 2**(j - 1)
            b(i, i) = 0.0_dp
         END DO
         IF (MAXVAL(ABS(b)) > tolerance) RETURN
      END DO
      chart = state%chart
      CALL skala_orbital_chart_set_reference(chart, adapted, metric, valid)
      IF (.NOT. valid) RETURN
      CALL skala_orbital_chart_restrict(chart, groups, has_anti)
      ALLOCATE (zero(skala_orbital_chart_size(chart)), final(n, n)); zero = 0.0_dp
      CALL skala_orbital_chart_evaluate(chart, zero, final)
      source = final
      DO i = 1, n
         source(:, i) = source(:, i)*state%occupations(i)
      END DO
      b = MATMUL(source, CONJG(TRANSPOSE(final)))
      info = -1
      IF (MAXVAL(ABS(b - original)) > tolerance) RETURN
      CALL skala_orbital_chart_release(state%chart)
      state%chart = chart
      state%coefficients = final
      info = 0
   END SUBROUTINE adapt_real_characters

! **************************************************************************************************
!> \brief Attach the complete reference-space group without changing the initial AO density.
!> \param state Initial native orbital state at zero coordinates.
!> \param metric Original charge metric.
!> \param kpoints Native cell-gauge symmetry metadata.
!> \param representative Irreducible k-point.
!> \param operations Complete operations leaving this point fixed modulo a reciprocal vector.
!> \param info Zero on success, minus one for an incompatible reference or occupation variable.
! **************************************************************************************************
   SUBROUTINE adapt_little_group(state, metric, kpoints, representative, operations, info)
      TYPE(qs_skala_orbitals_type), INTENT(INOUT)        :: state
      COMPLEX(dp), INTENT(IN)                            :: metric(:, :)
      TYPE(kpoint_type), POINTER                         :: kpoints
      INTEGER, INTENT(IN)                                :: representative, operations(:)
      INTEGER, INTENT(OUT)                               :: info

      COMPLEX(dp), ALLOCATABLE :: action(:, :), actions(:, :, :), coefficients(:, :), gauge(:, :), &
         original(:, :), source(:, :), weighted(:, :)
      INTEGER                                            :: i, n, op, status
      INTEGER, ALLOCATABLE                               :: groups(:)
      LOGICAL                                            :: valid
      LOGICAL, ALLOCATABLE                               :: antiunitary(:)
      REAL(dp)                                           :: rotation(3, 3)
      REAL(dp), ALLOCATABLE                              :: energies(:), zero(:)
      TYPE(skala_orbital_chart_type)                     :: chart

      info = -1
      n = SIZE(state%coefficients, 2)
      ALLOCATE (action(n, n), actions(n, n, SIZE(operations)), antiunitary(SIZE(operations)))
      DO op = 1, SIZE(operations)
         CALL kpoint_cell_gauge_action(kpoints, representative, operations(op), &
                                       kpoints%xkp(:, representative), action, rotation, antiunitary(op), status)
         IF (status /= 0) RETURN
         source = state%coefficients
         IF (antiunitary(op)) source = CONJG(source)
         actions(:, :, op) = MATMUL(CONJG(TRANSPOSE(state%coefficients)), MATMUL(metric, MATMUL(action, source)))
      END DO
      chart = state%chart
      IF (state%variable_occupations) THEN
         CALL skala_orbital_multiplets(actions, antiunitary, state%occupations, state%eigenvalues, &
                                       state%native%nmo, gauge, groups, energies, valid)
         IF (.NOT. valid) RETURN
         coefficients = MATMUL(state%coefficients, gauge)
         CALL skala_orbital_chart_set_reference(chart, coefficients, metric, valid)
         IF (.NOT. valid) RETURN
         CALL skala_orbital_chart_symmetry(chart, actions, antiunitary, valid, occupation_groups=groups)
      ELSE
         CALL skala_orbital_chart_symmetry(chart, actions, antiunitary, valid)
      END IF
      IF (.NOT. valid) RETURN
      CALL orbital_density(state, original)
      ALLOCATE (zero(skala_orbital_chart_size(chart)))
      IF (.NOT. ALLOCATED(coefficients)) ALLOCATE (coefficients(n, n))
      zero = 0.0_dp
      CALL skala_orbital_chart_evaluate(chart, zero, coefficients)
      weighted = coefficients
      DO i = 1, n
         weighted(:, i) = weighted(:, i)*state%occupations(i)
      END DO
      source = MATMUL(weighted, CONJG(TRANSPOSE(coefficients)))
      IF (MAXVAL(ABS(source - original)) > 3.e-11_dp) RETURN
      CALL skala_orbital_chart_release(state%chart)
      state%chart = chart
      state%coefficients = coefficients
      IF (state%variable_occupations) THEN
         state%occupation_groups = groups
         state%eigenvalues = energies
      END IF
      info = 0
   END SUBROUTINE adapt_little_group

! **************************************************************************************************
!> \brief Borrow native MOs and complete their missing virtual space in the original charge metric.
!>        Collective on the MO matrix communicator. The owner must outlive this orbital state.
!> \param state Owned chart and replicated coefficients, borrowing only the native MO set.
!> \param native Live spin and k-point MO set, including its fixed occupations.
!> \param metric Original normalization operator in the same AO order.
!> \param complex_orbitals Whether imaginary orbital variations are physical.
!> \param weight Positive k-point weight, applied once to physical tangent derivatives.
!> \param valid Whether the complete requested charge space is positive.
!> \param variable_occupations ...
! **************************************************************************************************
   SUBROUTINE qs_skala_orbitals_create(state, native, metric, complex_orbitals, weight, valid, variable_occupations)
      TYPE(qs_skala_orbitals_type), INTENT(INOUT)        :: state
      TYPE(mo_set_type), INTENT(INOUT), TARGET           :: native
      COMPLEX(dp), INTENT(IN)                            :: metric(:, :)
      LOGICAL, INTENT(IN)                                :: complex_orbitals
      REAL(dp), INTENT(IN)                               :: weight
      LOGICAL, INTENT(OUT)                               :: valid
      LOGICAL, INTENT(IN), OPTIONAL                      :: variable_occupations

      COMPLEX(dp), ALLOCATABLE                           :: completed(:, :), reference(:, :)
      INTEGER                                            :: i, n, nreference
      INTEGER, ALLOCATABLE                               :: occupied(:), virtuals(:)
      REAL(dp), ALLOCATABLE                              :: real_part(:, :), zero(:)

      CALL qs_skala_orbitals_release(state)
      IF (PRESENT(variable_occupations)) state%variable_occupations = variable_occupations
      n = native%nao
      nreference = native%nmo
      valid = n > 0 .AND. nreference > 0 .AND. nreference <= n
      CPASSERT(valid)
      CPASSERT(ALL(SHAPE(metric) == [n, n]))
      CPASSERT(ASSOCIATED(native%occupation_numbers))
      CPASSERT(SIZE(native%occupation_numbers) >= nreference)
      ALLOCATE (reference(n, nreference))
      IF (ASSOCIATED(native%cmo_coeff)) THEN
         CALL cp_cfm_get_submatrix(native%cmo_coeff, reference)
      ELSE
         CPASSERT(ASSOCIATED(native%mo_coeff))
         IF (native%use_mo_coeff_b) THEN
            CPASSERT(ASSOCIATED(native%mo_coeff_b))
            CALL copy_dbcsr_to_fm(native%mo_coeff_b, native%mo_coeff)
         END IF
         ALLOCATE (real_part(n, nreference))
         CALL cp_fm_get_submatrix(native%mo_coeff, real_part)
         reference(:, :) = CMPLX(real_part, 0.0_dp, dp)
      END IF
      CALL skala_orbital_chart_create(reference, metric, native%occupation_numbers(:nreference), &
                                      complex_orbitals, state%chart, valid, metric_weight=weight, complete_space=.TRUE., &
                                      variable_occupations=state%variable_occupations)
      IF (.NOT. valid) THEN
         ! A density-only restart can leave virtual or empty-spin MO columns uninitialized.
         ! Complete the occupied span in the original metric, preserving native occupations.
         ALLOCATE (occupied, SOURCE=PACK([(i, i=1, nreference)], native%occupation_numbers(:nreference) > 0.0_dp))
         ALLOCATE (virtuals, SOURCE=PACK([(i, i=1, nreference)], native%occupation_numbers(:nreference) == 0.0_dp))
         IF (SIZE(occupied) == 0) THEN
            reference = CMPLX(0.0_dp, 0.0_dp, dp)
            reference(1, 1) = CMPLX(1.0_dp, 0.0_dp, dp)
            CALL skala_orbital_chart_create(reference(:, :1), metric, [0.0_dp], complex_orbitals, &
                                            state%chart, valid, metric_weight=weight, complete_space=.TRUE.)
         ELSE
            CALL skala_orbital_chart_create(reference(:, occupied), metric, native%occupation_numbers(occupied), &
                                            complex_orbitals, state%chart, valid, metric_weight=weight, complete_space=.TRUE.)
         END IF
         IF (valid) THEN
            ALLOCATE (completed(n, n), zero(skala_orbital_chart_size(state%chart)))
            zero = 0.0_dp
            CALL skala_orbital_chart_evaluate(state%chart, zero, completed)
            IF (SIZE(occupied) == 0) THEN
               reference(:, :) = completed(:, :nreference)
            ELSE
               reference(:, occupied) = completed(:, :SIZE(occupied))
               reference(:, virtuals) = completed(:, SIZE(occupied) + 1:nreference)
            END IF
            DEALLOCATE (zero)
            CALL skala_orbital_chart_create(reference, metric, native%occupation_numbers(:nreference), &
                                            complex_orbitals, state%chart, valid, metric_weight=weight, complete_space=.TRUE., &
                                            variable_occupations=state%variable_occupations)
         END IF
      END IF
      IF (.NOT. valid) RETURN
      state%native => native
      state%weight = weight
      ALLOCATE (state%coefficients(n, n), state%occupations(n), state%eigenvalues(n), &
                zero(skala_orbital_chart_size(state%chart)))
      state%occupations = 0.0_dp
      state%occupations(:nreference) = native%occupation_numbers(:nreference)
      state%eigenvalues = 0.0_dp
      state%eigenvalues(:nreference) = native%eigenvalues(:nreference)
      zero = 0.0_dp
      CALL qs_skala_orbitals_evaluate(state, zero)
   END SUBROUTINE qs_skala_orbitals_create

! **************************************************************************************************
!> \brief Release only owned storage. The native environment retains its MO matrices.
!> \param state ...
! **************************************************************************************************
   SUBROUTINE qs_skala_orbitals_release(state)
      TYPE(qs_skala_orbitals_type), INTENT(INOUT)        :: state

      CALL skala_orbital_chart_release(state%chart)
      IF (ALLOCATED(state%coefficients)) DEALLOCATE (state%coefficients, state%occupations, state%eigenvalues)
      IF (ALLOCATED(state%occupation_groups)) DEALLOCATE (state%occupation_groups)
      NULLIFY (state%native)
      state%weight = 1.0_dp
      state%canonicalized = .FALSE.
      state%variable_occupations = .FALSE.
   END SUBROUTINE qs_skala_orbitals_release

! **************************************************************************************************
!> \brief Evaluate the complete fixed-occupation orbital trial without changing native storage.
!> \param state ...
!> \param coordinate ...
!> \param occupations ...
!> \param energies ...
! **************************************************************************************************
   SUBROUTINE qs_skala_orbitals_evaluate(state, coordinate, occupations, energies)
      TYPE(qs_skala_orbitals_type), INTENT(INOUT)        :: state
      REAL(dp), INTENT(IN)                               :: coordinate(:)
      REAL(dp), INTENT(IN), OPTIONAL                     :: occupations(:), energies(:)

      CPASSERT(ASSOCIATED(state%native))
      IF (PRESENT(occupations)) THEN
         CPASSERT(state%variable_occupations)
         CPASSERT(SIZE(occupations) == state%native%nmo)
         state%occupations = 0.0_dp
         state%occupations(:SIZE(occupations)) = occupations
         CALL skala_orbital_chart_set_occupations(state%chart, state%occupations)
      END IF
      IF (PRESENT(energies)) THEN
         CPASSERT(state%variable_occupations)
         CPASSERT(SIZE(energies) == state%native%nmo)
         state%eigenvalues(:SIZE(energies)) = energies
      END IF
      CALL skala_orbital_chart_evaluate(state%chart, coordinate, state%coefficients)
      state%canonicalized = .FALSE.
   END SUBROUTINE qs_skala_orbitals_evaluate

! **************************************************************************************************
!> \brief Canonicalize equal-occupation subspaces after physical orbital stationarity is established.
!>        This leaves the density unchanged and prepares eigenvalues for the weighted density.
!>        The caller supplies the same limiting Hamiltonian used by the complete functional.
!> \param state ...
!> \param hamiltonian Unweighted native Hermitian operator in the original AO order.
!> \param tolerance Physical tangent-gradient norm required before accepting weighted densities.
!> \param valid Whether stationarity and all subspace diagonalizations succeeded.
! **************************************************************************************************
   SUBROUTINE qs_skala_orbitals_canonicalize(state, hamiltonian, tolerance, valid)
      TYPE(qs_skala_orbitals_type), INTENT(INOUT)        :: state
      COMPLEX(dp), INTENT(IN)                            :: hamiltonian(:, :)
      REAL(dp), INTENT(IN)                               :: tolerance
      LOGICAL, INTENT(OUT)                               :: valid

      COMPLEX(dp), ALLOCATABLE                           :: projected(:, :), rotation(:, :), &
                                                            trial(:, :), work(:)
      INTEGER                                            :: i, info, n, ng
      INTEGER, ALLOCATABLE                               :: indices(:)
      LOGICAL, ALLOCATABLE                               :: done(:)
      REAL(dp)                                           :: block_mean
      REAL(dp), ALLOCATABLE                              :: eigenvalues(:), gradient(:), &
                                                            group_values(:), real_work(:)

      CPASSERT(ASSOCIATED(state%native))
      n = SIZE(state%occupations)
      CPASSERT(ALL(SHAPE(hamiltonian) == [n, n]))
      CPASSERT(tolerance > 0.0_dp)
      state%canonicalized = .FALSE.
      ALLOCATE (gradient(skala_orbital_chart_size(state%chart)))
      CALL skala_orbital_chart_derivative(state%chart, state%weight*hamiltonian, gradient, physical_tangent=.TRUE.)
      valid = NORM2(gradient) <= tolerance
      IF (.NOT. valid) RETURN
      ALLOCATE (projected(n, n), trial(n, n), eigenvalues(n), done(n))
      projected(:, :) = MATMUL(CONJG(TRANSPOSE(state%coefficients)), MATMUL(hamiltonian, state%coefficients))
      trial(:, :) = state%coefficients
      done = .FALSE.
      DO i = 1, n
         IF (done(i)) CYCLE
         ng = COUNT(state%occupations == state%occupations(i))
         ALLOCATE (indices(ng), rotation(ng, ng), group_values(ng), &
                   work(MAX(1, 2*ng - 1)), real_work(MAX(1, 3*ng - 2)))
         indices(:) = PACK([(info, info=1, n)], state%occupations == state%occupations(i))
         done(indices) = .TRUE.
         rotation(:, :) = projected(indices, indices)
         IF (ng > 1) THEN
            block_mean = SUM([(REAL(rotation(info, info), dp), info=1, ng)])/REAL(ng, dp)
            DO info = 1, ng
               rotation(info, info) = rotation(info, info) - block_mean
            END DO
            ! Avoid arbitrary rotations selected by roundoff inside a degenerate eigenspace.
            ! The caller still checks the rebuilt unrestricted gradient and thermal stationarity.
            IF (SQRT(SUM(ABS(rotation)**2)) <= tolerance) THEN
               eigenvalues(indices) = block_mean
               DEALLOCATE (indices, rotation, group_values, work, real_work)
               CYCLE
            END IF
            rotation(:, :) = projected(indices, indices)
         END IF
         CALL zheev('V', 'U', ng, rotation, ng, group_values, work, SIZE(work), real_work, info)
         valid = info == 0
         IF (.NOT. valid) RETURN
         eigenvalues(indices) = group_values
         trial(:, indices) = MATMUL(state%coefficients(:, indices), rotation)
         DEALLOCATE (indices, rotation, group_values, work, real_work)
      END DO
      state%coefficients(:, :) = trial
      state%eigenvalues(:) = eigenvalues
      state%canonicalized = .TRUE.
   END SUBROUTINE qs_skala_orbitals_canonicalize

! **************************************************************************************************
!> \brief Install matching Gamma MOs and density, optionally including the canonical weighted density.
!>        The caller invalidates and rebuilds rho and tau before the next KS evaluation.
!> \param state ...
!> \param frame Original inverse AO layout, also valid for its DIRECT representation.
!> \param density Native AO spin density, overwritten without changing occupations.
!> \param weighted Optional native energy-weighted density, requiring canonicalized stationary MOs.
! **************************************************************************************************
   SUBROUTINE qs_skala_orbitals_install_gamma(state, frame, density, weighted)
      TYPE(qs_skala_orbitals_type), INTENT(INOUT)        :: state
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      TYPE(dbcsr_type), INTENT(INOUT)                    :: density
      TYPE(dbcsr_type), INTENT(INOUT), OPTIONAL          :: weighted

      COMPLEX(dp), ALLOCATABLE                           :: matrix(:, :)
      REAL(dp)                                           :: defect, tolerance
      TYPE(dbcsr_type)                                   :: imaginary, real_part

      CALL write_native_mos(state)
      CALL orbital_density(state, matrix)
      defect = MAXVAL(ABS(AIMAG(matrix)))
      tolerance = 128.0_dp*EPSILON(1.0_dp)*MAX(1.0_dp, MAXVAL(ABS(matrix)))
      CPASSERT(defect <= tolerance)
      CALL qs_upf_inverse_matrix_scatter(frame, matrix, density, real_part, imaginary)
      CALL dbcsr_copy(density, real_part)
      CALL dbcsr_release(real_part)
      CALL dbcsr_release(imaginary)
      IF (PRESENT(weighted)) THEN
         CALL orbital_density(state, matrix, energy_weighted=.TRUE.)
         CALL qs_upf_inverse_matrix_scatter(frame, matrix, weighted, real_part, imaginary)
         CALL dbcsr_copy(weighted, real_part)
         CALL dbcsr_release(real_part)
         CALL dbcsr_release(imaginary)
      END IF
   END SUBROUTINE qs_skala_orbitals_install_gamma

! **************************************************************************************************
!> \brief Install matching group-local Bloch MOs, real/imaginary P(k), and its accepted native cache.
!>        No k weight is included in P(k). The caller updates every owned spin/k block before
!>        rebuilding the physical density and evaluating KS. Existing density-mixing trials expire.
!> \param state ...
!> \param kpoints Live native group environment.
!> \param local Local k-point index within this owning group.
!> \param spin Native spin index.
!> \param energy_weighted Also update W(k), after canonicalization with the accepted Hamiltonian.
! **************************************************************************************************
   SUBROUTINE qs_skala_orbitals_install_kpoint(state, kpoints, local, spin, energy_weighted)
      TYPE(qs_skala_orbitals_type), INTENT(INOUT)        :: state
      TYPE(kpoint_type), INTENT(INOUT)                   :: kpoints
      INTEGER, INTENT(IN)                                :: local, spin
      LOGICAL, INTENT(IN), OPTIONAL                      :: energy_weighted

      COMPLEX(dp), ALLOCATABLE                           :: matrix(:, :)
      REAL(dp)                                           :: defect, tolerance
      TYPE(kpoint_env_type), POINTER                     :: kp

      CPASSERT(ALLOCATED(kpoints%native_density))
      CPASSERT(local >= 1 .AND. local <= SIZE(kpoints%kp_env))
      kp => kpoints%kp_env(local)%kpoint_env
      CPASSERT(spin >= 1 .AND. spin <= SIZE(kp%mos))
      CPASSERT(ASSOCIATED(state%native, kp%mos(spin)))
      CALL write_native_mos(state)
      CALL orbital_density(state, matrix)
      CALL cp_cfm_set_submatrix(kpoints%native_density(local, spin), matrix)
      CALL cp_fm_set_submatrix(kp%pmat(1, spin), REAL(matrix, dp))
      IF (kpoints%use_real_wfn) THEN
         defect = MAXVAL(ABS(AIMAG(matrix)))
         tolerance = 128.0_dp*EPSILON(1.0_dp)*MAX(1.0_dp, MAXVAL(ABS(matrix)))
         CPASSERT(defect <= tolerance)
      ELSE
         CALL cp_fm_set_submatrix(kp%pmat(2, spin), AIMAG(matrix))
      END IF
      kpoints%native_density_seeded = .TRUE.
      kpoints%native_trial_valid = .FALSE.
      IF (PRESENT(energy_weighted)) THEN
         IF (energy_weighted) THEN
            CPASSERT(ASSOCIATED(kp%wmat))
            CALL orbital_density(state, matrix, energy_weighted=.TRUE.)
            CALL cp_fm_set_submatrix(kp%wmat(1, spin), REAL(matrix, dp))
            IF (.NOT. kpoints%use_real_wfn) CALL cp_fm_set_submatrix(kp%wmat(2, spin), AIMAG(matrix))
         END IF
      END IF
   END SUBROUTINE qs_skala_orbitals_install_kpoint

! **************************************************************************************************
!> \brief Store only the number of native MOs requested by the SCF setup.
!> \param state ...
! **************************************************************************************************
   SUBROUTINE write_native_mos(state)
      TYPE(qs_skala_orbitals_type), INTENT(INOUT)        :: state

      INTEGER                                            :: i, nreference
      REAL(dp)                                           :: defect, tolerance

      CPASSERT(ASSOCIATED(state%native))
      nreference = state%native%nmo
      IF (ASSOCIATED(state%native%cmo_coeff)) THEN
         CALL cp_cfm_set_submatrix(state%native%cmo_coeff, state%coefficients(:, :nreference))
      ELSE
         defect = MAXVAL(ABS(AIMAG(state%coefficients(:, :nreference))))
         tolerance = 128.0_dp*EPSILON(1.0_dp)*MAX(1.0_dp, MAXVAL(ABS(state%coefficients)))
         CPASSERT(defect <= tolerance)
         CALL cp_fm_set_submatrix(state%native%mo_coeff, REAL(state%coefficients(:, :nreference), dp))
         IF (state%native%use_mo_coeff_b) THEN
            CALL copy_fm_to_dbcsr(state%native%mo_coeff, state%native%mo_coeff_b)
         END IF
      END IF
      IF (state%variable_occupations) THEN
         state%native%occupation_numbers(:nreference) = state%occupations(:nreference)
         state%native%eigenvalues(:nreference) = state%eigenvalues(:nreference)
         state%native%homo = 0
         state%native%lfomo = nreference + 1
         DO i = 1, nreference
            IF (state%occupations(i) > 0.0_dp) state%native%homo = i
            IF (state%occupations(i) < state%native%maxocc) state%native%lfomo = MIN(state%native%lfomo, i)
         END DO
         state%native%uniform_occupation = ALL(state%occupations(:nreference) == state%occupations(1))
      ELSE IF (state%canonicalized) THEN
         state%native%eigenvalues(:nreference) = state%eigenvalues(:nreference)
      END IF
   END SUBROUTINE write_native_mos

! **************************************************************************************************
!> \brief Form P=C f C^dagger or W=C f epsilon C^dagger without applying a k weight.
!> \param state ...
!> \param density ...
!> \param energy_weighted ...
! **************************************************************************************************
   SUBROUTINE orbital_density(state, density, energy_weighted)
      TYPE(qs_skala_orbitals_type), INTENT(IN)           :: state
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: density(:, :)
      LOGICAL, INTENT(IN), OPTIONAL                      :: energy_weighted

      COMPLEX(dp), ALLOCATABLE                           :: weighted(:, :)
      INTEGER                                            :: i, n
      LOGICAL                                            :: use_energy

      use_energy = .FALSE.
      IF (PRESENT(energy_weighted)) use_energy = energy_weighted
      IF (use_energy) THEN
         CPASSERT(state%canonicalized)
      END IF
      n = SIZE(state%occupations)
      ALLOCATE (density(n, n), weighted(n, n))
      DO i = 1, n
         weighted(:, i) = state%occupations(i)*state%coefficients(:, i)
         IF (use_energy) weighted(:, i) = state%eigenvalues(i)*weighted(:, i)
      END DO
      density(:, :) = MATMUL(weighted, CONJG(TRANSPOSE(state%coefficients)))
   END SUBROUTINE orbital_density
END MODULE qs_skala_orbitals
