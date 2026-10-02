!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Distributed native RI-HFX transition tensors, batch adjoints and geometry derivatives.
! **************************************************************************************************
PROGRAM qs_uspp_hfx_unittest
   USE atom_upf,                        ONLY: atom_release_upf,&
                                              atom_upfpot_type
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE basis_set_types,                 ONLY: allocate_gto_basis_set,&
                                              deallocate_gto_basis_set,&
                                              gto_basis_set_p_type
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_create, dbcsr_distribution_new, dbcsr_distribution_release, dbcsr_distribution_type, &
        dbcsr_finalize, dbcsr_get_info, dbcsr_get_stored_coordinates, dbcsr_iterator_blocks_left, &
        dbcsr_iterator_next_block, dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, &
        dbcsr_put_block, dbcsr_release, dbcsr_type
   USE dbt_api,                         ONLY: &
        dbt_destroy, dbt_get_block, dbt_get_info, dbt_get_stored_coordinates, &
        dbt_iterator_blocks_left, dbt_iterator_next_block, dbt_iterator_start, dbt_iterator_stop, &
        dbt_iterator_type, dbt_pgrid_create, dbt_pgrid_destroy, dbt_pgrid_type, dbt_put_block, &
        dbt_type
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: dp
   USE libcp2k,                         ONLY: cp2k_finalize,&
                                              cp2k_init
   USE mathconstants,                   ONLY: fourpi,&
                                              pi
   USE message_passing,                 ONLY: mp_comm_type,&
                                              mp_comm_world
   USE particle_types,                  ONLY: particle_type
   USE qs_tensors_types,                ONLY: create_3c_tensor
   USE qs_uspp_hfx,                     ONLY: qs_uspp_hfx_native_gradient,&
                                              qs_uspp_hfx_native_tensor,&
                                              qs_uspp_hfx_transition_adjoint,&
                                              qs_uspp_hfx_transition_tensor
   USE qs_uspp_projectors,              ONLY: qs_uspp_projector_gradient
   USE qs_uspp_types,                   ONLY: qs_uspp_type
   USE sap_kind_types,                  ONLY: release_sap_int,&
                                              sap_int_type
   USE upf_projector_integrals,         ONLY: upf_projector_init,&
                                              upf_projector_type
   USE uspp_augmentation,               ONLY: uspp_augmentation_init

   IMPLICIT NONE

   INTEGER, PARAMETER :: natom = 3, nao = 4, naux = 6, ncenter = 2
   INTEGER, PARAMETER :: ao_atom(nao) = [1, 2, 2, 3], center_atom(ncenter) = [1, 3]
   REAL(dp), PARAMETER :: beta = 0.9_dp, alpha(2) = [0.8_dp, 1.3_dp], coefficient(2) = [1.0_dp, -0.7_dp]
   REAL(dp), PARAMETER :: amplitude(nao) = [0.65_dp, -0.42_dp, 0.31_dp, -0.56_dp]
   REAL(dp), PARAMETER :: gamma(nao) = [0.7_dp, 0.55_dp, 0.9_dp, 0.62_dp]
   REAL(dp), PARAMETER :: anisotropy(3) = [1.0_dp, 1.3_dp, 0.7_dp]
   REAL(dp), PARAMETER :: fraction = 0.27_dp, spin = 0.5_dp, pref = -fraction*spin, step = 2.e-5_dp
   INTEGER, TARGET :: ao_sizes(natom) = [1, 2, 1], pp_sizes(ncenter) = [1, 1]
   INTEGER, TARGET :: ao_rows(natom), pp_cols(ncenter)
   INTEGER, ALLOCATABLE, TARGET :: grid(:, :)
   INTEGER, POINTER :: matrix_grid(:, :), row_sizes(:), col_sizes(:), row_dist(:), col_dist(:)
   INTEGER :: atom, axis, case_id, i, j, nprows, npcols, rank, nranks, dims(3), target_dims(3)
   INTEGER :: nang(1) = [194]
   REAL(dp) :: positions(3, natom), density(nao, nao), inverse(naux, naux), smooth(naux, nao, nao)
   REAL(dp) :: max_error
   TYPE(atomic_kind_type), TARGET :: atomic_kind
   TYPE(particle_type) :: particles(natom)
   TYPE(gto_basis_set_p_type) :: basis(1)
   TYPE(qs_uspp_type) :: uspp
   TYPE(dbcsr_distribution_type) :: distribution
   TYPE(dbt_pgrid_type) :: native_grid, target_grid
   TYPE(mp_comm_type) :: comm

   CALL cp2k_init()
   CALL comm%set_handle(mp_comm_world%get_handle())
   rank = comm%mepos
   nranks = comm%num_pe
   max_error = 0.0_dp
   positions(:, 1) = [0.1_dp, -0.2_dp, 0.3_dp]
   positions(:, 2) = [-0.45_dp, 0.6_dp, -0.15_dp]
   positions(:, 3) = [0.8_dp, 0.15_dp, -0.4_dp]
   atomic_kind%kind_number = 1
   DO atom = 1, natom
      particles(atom)%atomic_kind => atomic_kind
      particles(atom)%r = positions(:, atom)
   END DO
   CALL make_native_data()
   DO j = 1, nao
      DO i = 1, nao
         density(i, j) = 0.035_dp*COS(REAL(i - j, dp)) + 0.025_dp*SIN(REAL(i, dp))*SIN(REAL(j, dp))
      END DO
      density(j, j) = density(j, j) + 0.8_dp + 0.07_dp*j
   END DO
   DO j = 1, naux
      DO i = 1, naux
         inverse(i, j) = 0.08_dp*COS(REAL(i - j, dp))
      END DO
      inverse(j, j) = inverse(j, j) + 0.45_dp + 0.02_dp*j
      DO axis = 1, nao
         DO atom = 1, nao
            smooth(j, axis, atom) = 0.02_dp*SIN(REAL(j + axis + atom, dp))
         END DO
      END DO
   END DO

   ! Matrix and tensor grids change independently, including ranks without native centers.
   DO case_id = 1, 2
      nprows = nranks
      npcols = 1
      dims = [nranks, 1, 1]
      target_dims = [1, 1, nranks]
      IF (case_id == 2) THEN
         nprows = 1
         npcols = nranks
         dims = [1, nranks, 1]
         target_dims = [nranks, 1, 1]
      END IF
      ALLOCATE (grid(0:nprows - 1, 0:npcols - 1))
      grid(:, :) = RESHAPE([(i, i=0, nranks - 1)], [nprows, npcols])
      matrix_grid => grid
      ao_rows = [(MODULO(i - 1, nprows), i=1, natom)]
      pp_cols = [(MODULO(i - 1, npcols), i=1, ncenter)]
      row_dist => ao_rows
      col_dist => pp_cols
      CALL dbcsr_distribution_new(distribution, group=comm%get_handle(), pgrid=matrix_grid, &
                                  row_dist=row_dist, col_dist=col_dist)
      CALL dbt_pgrid_create(comm, dims, native_grid)
      CALL dbt_pgrid_create(comm, target_dims, target_grid)
      CALL check_distributed_path()
      CALL dbt_pgrid_destroy(target_grid)
      CALL dbt_pgrid_destroy(native_grid)
      CALL dbcsr_distribution_release(distribution)
      DEALLOCATE (grid)
      NULLIFY (matrix_grid)
   END DO
   CALL deallocate_gto_basis_set(basis(1)%gto_basis_set)
   CALL comm%max(max_error)
   IF (rank == 0) WRITE (*, '(A,I0,A,ES16.7)') &
      'Native RI-HFX checks: ranks=', nranks, ', maximum absolute error=', max_error
   CALL cp2k_finalize()

CONTAINS

! **************************************************************************************************
!> \brief One original radial s-projector and two independently contracted Gaussian s shell sets.
! **************************************************************************************************
   SUBROUTINE make_native_data()
      INTEGER, PARAMETER                                 :: nr = 601

      INTEGER                                            :: ir
      TYPE(atom_upfpot_type)                             :: pot
      TYPE(upf_projector_type)                           :: projectors

      pot%is_ultrasoft = .TRUE.
      pot%number_of_proj = 1
      ALLOCATE (pot%lbeta(1), pot%dion(1, 1), pot%qij(1, 1), pot%r(nr), pot%rab(nr), &
                pot%beta(nr, 1), pot%qijl(nr, 1, 0:0))
      pot%lbeta = 0
      pot%dion = 0.0_dp
      pot%qij = 0.0_dp
      DO ir = 1, nr
         pot%r(ir) = 8.0_dp*REAL(ir - 1, dp)/REAL(nr - 1, dp)
      END DO
      pot%rab = 8.0_dp/REAL(nr - 1, dp)
      pot%beta(:, 1) = EXP(-pot%r**2)
      pot%qijl(:, 1, 0) = pot%r**2*EXP(-beta*pot%r**2)
      CALL upf_projector_init(pot, 0, projectors)
      ALLOCATE (uspp%augmentation(1), uspp%center_atoms(ncenter), uspp%center_kinds(ncenter), &
                uspp%projector_sizes(ncenter), uspp%position(3, ncenter))
      CALL uspp_augmentation_init(pot, projectors, uspp%augmentation(1))
      IF (uspp%augmentation(1)%nproj /= 1) ERROR STOP 'Wrong native s-projector count'
      uspp%center_atoms(:) = center_atom
      uspp%center_kinds = 1
      uspp%projector_sizes(:) = pp_sizes
      uspp%position(:, :) = positions(:, center_atom)
      CALL atom_release_upf(pot)
      CALL allocate_gto_basis_set(basis(1)%gto_basis_set)
      ASSOCIATE (gto => basis(1)%gto_basis_set)
         gto%nset = 2
         gto%nsgf = 2
         ALLOCATE (gto%lmax(2), gto%npgf(2), gto%nsgf_set(2), gto%zet(1, 2), gto%sphi(1, 2))
         gto%lmax = 0
         gto%npgf = 1
         gto%nsgf_set = 1
         gto%zet(1, :) = alpha
         gto%sphi(1, :) = coefficient
      END ASSOCIATE
   END SUBROUTINE make_native_data

! **************************************************************************************************
!> \brief Real DBT forward and two independent, non-symmetric first-AO batches of its reverse.
! **************************************************************************************************
   SUBROUTINE check_distributed_path()
      INTEGER                                            :: a, batch, high, idir, jdir, low, mu, nu, &
                                                            p, sample
      INTEGER, ALLOCATABLE                               :: d1(:), d2(:), d3(:)
      REAL(dp) :: adjoint(naux, nao, nao), atom_gradient(3, natom), b(nao, ncenter), &
         batched(naux, nao, nao), c(naux, ncenter, ncenter), cell_gradient(3, 3), &
         cref(naux, ncenter, ncenter), cross, d_b(nao, ncenter), d_c(naux, ncenter, ncenter), &
         direction(3), energy, fd, gb(nao, ncenter), gb_batch(nao, ncenter), gb_ref(nao, ncenter), &
         gc(naux, ncenter, ncenter), gc_batch(naux, ncenter, ncenter), &
         gc_ref(naux, ncenter, ncenter), left(naux, nao, nao), minus, moved(3, natom), &
         native(naux, nao, nao), native_ref(naux, nao, nao), opposite(naux, nao, nao), &
         perturbed(naux, nao, nao), plus
      REAL(dp) :: right(naux, nao, nao), total(naux, nao, nao)
      TYPE(dbcsr_type)                                   :: overlap_gradient
      TYPE(dbcsr_type), TARGET                           :: overlaps
      TYPE(dbt_type)                                     :: auxiliary, auxiliary_gradient, tensor, &
                                                            tensor_adjoint
      TYPE(sap_int_type), POINTER                        :: sap(:)

      CALL closed_data(positions, b, cref)
      CALL make_overlap(b, overlaps)
      CALL make_sap(b, overlaps, sap)
      CALL qs_uspp_hfx_native_tensor(uspp, basis, particles, nang, native_grid, auxiliary)
      CALL gather_tensor(auxiliary, c)
      CALL check(MAXVAL(ABS(c - cref)), 'native C against closed Gaussian Coulomb', 3.e-8_dp)
      CALL create_3c_tensor(tensor, d1, d2, d3, target_grid, [(1, p=1, naux)], &
                            [(1, mu=1, nao)], [(1, nu=1, nao)], [1], [2, 3], name='Split native RI AO AO')
      DEALLOCATE (d1, d2, d3)
      CALL qs_uspp_hfx_transition_tensor(overlaps, auxiliary, tensor)
      CALL gather_tensor(tensor, native)
      CALL dense_transition(b, c, native_ref)
      CALL check(MAXVAL(ABS(native - native_ref)), 'full transition reblocking and same-atom factor')
      DO p = 1, naux
         CALL check(MAXVAL(ABS(native(p, :, :) - TRANSPOSE(native(p, :, :)))), 'full AO-pair symmetry')
      END DO
      IF (ABS(native(1, 2, 3)) < 1.e-6_dp) ERROR STOP 'Missing same-atom off-diagonal augmentation'
      CALL qs_uspp_hfx_transition_tensor(overlaps, auxiliary, tensor)
      CALL gather_tensor(tensor, total)
      CALL check(MAXVAL(ABS(total - 2.0_dp*native_ref)), 'forward tensor increments without packing')
      total = smooth + native
      energy = exchange_energy(total)
      CALL exchange_opposite(total, opposite)
      CALL dense_transition(b, cref, native_ref)
      CALL check(ABS(energy - exchange_energy(smooth + native_ref)), 'native RI-HFX energy oracle', 3.e-8_dp)
      CALL dense_transition(b(:, :1), c(:, :1, :1), left)
      CALL dense_transition(b(:, 2:), c(:, 2:, 2:), right)
      cross = exchange_energy(left + right) - exchange_energy(left) - exchange_energy(right)
      IF (ABS(cross) < 1.e-6_dp) ERROR STOP 'Fixture lost intercenter QQ exchange'
      IF (ABS(exchange_energy(total, density/2.0_dp, 1.0_dp)*2.0_dp - energy) > 2.e-12_dp) THEN
         ERROR STOP 'RI-HFX RKS physical-spin factor'
      END IF

      gb = 0.0_dp
      gc = 0.0_dp
      atom_gradient = 0.0_dp
      cell_gradient = 0.0_dp
      DO batch = 1, 2
         low = 1
         high = 2
         IF (batch == 2) THEN
            low = 3
            high = nao
         END IF
         batched = 0.0_dp
         batched(:, low:high, :) = opposite(:, low:high, :)
         CALL create_3c_tensor(tensor_adjoint, d1, d2, d3, target_grid, [(1, p=1, naux)], &
                               [(1, mu=1, nao)], [(1, nu=1, nao)], [1], [2, 3], name='First AO batch adjoint')
         DEALLOCATE (d1, d2, d3)
         CALL put_tensor(tensor_adjoint, batched)
         CALL qs_uspp_hfx_transition_adjoint(overlaps, auxiliary, tensor_adjoint, pref, &
                                             overlap_gradient, auxiliary_gradient)
         CALL gather_matrix(overlap_gradient, gb_batch)
         CALL gather_tensor(auxiliary_gradient, gc_batch)
         adjoint = pref*batched
         CALL dense_reverse(b, c, adjoint, gb_ref, gc_ref)
         CALL check(MAXVAL(ABS(gb_batch - gb_ref)), 'both B adjoints in non-symmetric AO batch')
         CALL check(MAXVAL(ABS(gc_batch - gc_ref)), 'center-diagonal C adjoint in AO batch')
         CALL check(MAXVAL(ABS(gc_batch(:, 1, 2))) + MAXVAL(ABS(gc_batch(:, 2, 1))), &
                    'C adjoint excludes off-center variables')
         CALL gather_tensor(tensor_adjoint, perturbed)
         CALL check(MAXVAL(ABS(perturbed - batched)), 'reverse leaves opposite tensor unchanged')
         gb = gb + gb_batch
         gc = gc + gc_batch
         CALL qs_uspp_projector_gradient(overlap_gradient, sap, center_atom, atom_gradient, cell_gradient)
         CALL qs_uspp_hfx_native_gradient(uspp, basis, particles, nang, auxiliary_gradient, &
                                          atom_gradient, cell_gradient)
         CALL dbcsr_release(overlap_gradient)
         CALL dbt_destroy(auxiliary_gradient)
         CALL dbt_destroy(tensor_adjoint)
      END DO
      CALL dense_reverse(b, c, pref*opposite, gb_ref, gc_ref)
      CALL check(MAXVAL(ABS(gb - gb_ref)), 'accumulated batch B adjoint')
      CALL check(MAXVAL(ABS(gc - gc_ref)), 'accumulated batch C adjoint')
      ! Independent energy differences also test the prefactor and both quadratic B occurrences.
      DO a = 1, ncenter
         DO mu = 1, nao
            d_b = 0.0_dp
            d_b(mu, a) = step
            CALL dense_transition(b + d_b, c, perturbed)
            plus = exchange_energy(smooth + perturbed)
            CALL dense_transition(b - d_b, c, perturbed)
            minus = exchange_energy(smooth + perturbed)
            CALL check(ABS((plus - minus)/(2.0_dp*step) - gb(mu, a)), 'B adjoint energy FD', 2.e-8_dp)
         END DO
         DO p = 1, naux
            d_c = 0.0_dp
            d_c(p, a, a) = step
            CALL dense_transition(b, c + d_c, perturbed)
            plus = exchange_energy(smooth + perturbed)
            CALL dense_transition(b, c - d_c, perturbed)
            minus = exchange_energy(smooth + perturbed)
            CALL check(ABS((plus - minus)/(2.0_dp*step) - gc(p, a, a)), 'C adjoint energy FD', 2.e-9_dp)
         END DO
      END DO
      ! Every distributed contribution is reduced exactly once, after both AO batches.
      CALL comm%sum(atom_gradient)
      CALL comm%sum(cell_gradient)
      CALL check(MAXVAL(ABS(SUM(atom_gradient, DIM=2))), 'native and SAP rigid translation sum rule')
      DO sample = 1, 3
         atom = 2
         direction = [0.3_dp, -0.2_dp, 0.4_dp]
         IF (sample == 2) THEN
            atom = 1
            direction = [-0.2_dp, 0.5_dp, 0.1_dp]
         ELSE IF (sample == 3) THEN
            atom = 3
            direction = [0.4_dp, 0.1_dp, -0.3_dp]
         END IF
         moved = positions
         moved(:, atom) = moved(:, atom) + step*direction
         plus = geometry_energy(moved)
         moved(:, atom) = moved(:, atom) - 2.0_dp*step*direction
         minus = geometry_energy(moved)
         fd = (plus - minus)/(2.0_dp*step)
         CALL check(ABS(SUM(atom_gradient(:, atom)*direction) - fd), 'native plus SAP atomic direction FD', 8.e-8_dp)
      END DO
      DO sample = 1, 3
         idir = 1
         jdir = 1
         IF (sample == 2) jdir = 2
         IF (sample == 3) idir = 2
         moved = positions
         moved(idir, :) = moved(idir, :) + step*positions(jdir, :)
         plus = geometry_energy(moved)
         moved = positions
         moved(idir, :) = moved(idir, :) - step*positions(jdir, :)
         minus = geometry_energy(moved)
         fd = (plus - minus)/(2.0_dp*step)
         CALL check(ABS(cell_gradient(idir, jdir) - fd), 'native plus original SAP affine strain FD', 8.e-8_dp)
      END DO
      IF (ABS(cell_gradient(1, 2) - cell_gradient(2, 1)) < 1.e-6_dp) THEN
         ERROR STOP 'Synthetic SAP fixture cannot distinguish transposed shear gradients'
      END IF
      moved = positions + SPREAD([0.31_dp, -0.27_dp, 0.19_dp], DIM=2, NCOPIES=natom)
      CALL check(ABS(geometry_energy(moved) - geometry_energy(positions)), 'closed energy rigid translation')
      CALL release_sap_int(sap)
      CALL dbcsr_release(overlaps)
      CALL dbt_destroy(tensor)
      CALL dbt_destroy(auxiliary)
   END SUBROUTINE check_distributed_path

! **************************************************************************************************
!> \brief Closed Gaussian-Gaussian erf formula and analytic overlap values, independent of native grids.
!> \param coords ...
!> \param b ...
!> \param c ...
! **************************************************************************************************
   SUBROUTINE closed_data(coords, b, c)
      REAL(dp), INTENT(IN)                               :: coords(3, natom)
      REAL(dp), INTENT(OUT)                              :: b(nao, ncenter), &
                                                            c(naux, ncenter, ncenter)

      INTEGER                                            :: a, atom, mu, p, shell
      REAL(dp)                                           :: distance(3), f0, kappa, t

      c = 0.0_dp
      DO a = 1, ncenter
         DO mu = 1, nao
            distance = coords(:, ao_atom(mu)) - coords(:, center_atom(a))
            b(mu, a) = amplitude(mu)*EXP(-gamma(mu)*SUM(anisotropy*distance**2))
         END DO
         DO atom = 1, natom
            distance = coords(:, center_atom(a)) - coords(:, atom)
            DO shell = 1, 2
               p = 2*(atom - 1) + shell
               kappa = alpha(shell)*beta/(alpha(shell) + beta)
               t = kappa*SUM(distance**2)
               f0 = 1.0_dp
               IF (t > 0.0_dp) f0 = SQRT(pi)*ERF(SQRT(t))/(2.0_dp*SQRT(t))
               c(p, a, a) = coefficient(shell)*2.0_dp*pi**2.5_dp*f0/ &
                            (fourpi*alpha(shell)*beta*SQRT(alpha(shell) + beta))
            END DO
         END DO
      END DO
   END SUBROUTINE closed_data

! **************************************************************************************************
!> \brief Explicit projector-index forward, independent of DBT layout or contractions.
!> \param b ...
!> \param c ...
!> \param tensor ...
! **************************************************************************************************
   SUBROUTINE dense_transition(b, c, tensor)
      REAL(dp), INTENT(IN)                               :: b(:, :), c(:, :, :)
      REAL(dp), INTENT(OUT)                              :: tensor(naux, nao, nao)

      INTEGER                                            :: a, d, mu, nu, p

      tensor = 0.0_dp
      DO nu = 1, nao
         DO mu = 1, nao
            DO d = 1, SIZE(b, 2)
               DO a = 1, SIZE(b, 2)
                  DO p = 1, naux
                     tensor(p, mu, nu) = tensor(p, mu, nu) + b(mu, a)*c(p, a, d)*b(nu, d)
                  END DO
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE dense_transition

! **************************************************************************************************
!> \brief Both B derivatives are explicit; only center-diagonal C variables are admissible.
!> \param b ...
!> \param c ...
!> \param adjoint ...
!> \param gb ...
!> \param gc ...
! **************************************************************************************************
   SUBROUTINE dense_reverse(b, c, adjoint, gb, gc)
      REAL(dp), INTENT(IN)                               :: b(nao, ncenter), &
                                                            c(naux, ncenter, ncenter), &
                                                            adjoint(naux, nao, nao)
      REAL(dp), INTENT(OUT)                              :: gb(nao, ncenter), &
                                                            gc(naux, ncenter, ncenter)

      INTEGER                                            :: a, d, mu, nu, p

      gb = 0.0_dp
      gc = 0.0_dp
      DO a = 1, ncenter
         DO nu = 1, nao
            DO mu = 1, nao
               DO p = 1, naux
                  gc(p, a, a) = gc(p, a, a) + adjoint(p, mu, nu)*b(mu, a)*b(nu, a)
                  DO d = 1, ncenter
                     gb(mu, a) = gb(mu, a) + adjoint(p, mu, nu)*c(p, a, d)*b(nu, d)
                     gb(nu, a) = gb(nu, a) + adjoint(p, mu, nu)*c(p, d, a)*b(mu, d)
                  END DO
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE dense_reverse

! **************************************************************************************************
!> \brief Four-AO-index exchange energy, with a fixed RI metric and physical spin convention.
!> \param tensor ...
!> \param pmat ...
!> \param spin_factor ...
!> \return ...
! **************************************************************************************************
   FUNCTION exchange_energy(tensor, pmat, spin_factor) RESULT(energy)
      REAL(dp), INTENT(IN)                               :: tensor(naux, nao, nao)
      REAL(dp), INTENT(IN), OPTIONAL                     :: pmat(nao, nao), spin_factor
      REAL(dp)                                           :: energy

      INTEGER                                            :: kap, lam, mu, nu, p, q
      REAL(dp)                                           :: factor, integral, p_use(nao, nao)

      p_use = density
      IF (PRESENT(pmat)) p_use = pmat
      factor = spin
      IF (PRESENT(spin_factor)) factor = spin_factor
      energy = 0.0_dp
      DO lam = 1, nao
         DO kap = 1, nao
            DO nu = 1, nao
               DO mu = 1, nao
                  integral = 0.0_dp
                  DO q = 1, naux
                     DO p = 1, naux
                        integral = integral + tensor(p, mu, nu)*inverse(p, q)*tensor(q, kap, lam)
                     END DO
                  END DO
                  energy = energy - 0.5_dp*fraction*factor*p_use(mu, kap)*p_use(nu, lam)*integral
               END DO
            END DO
         END DO
      END DO
   END FUNCTION exchange_energy

! **************************************************************************************************
!> \brief Opposite tensor for the full RI-HFX energy; the exchange prefactor is kept separate.
!> \param tensor ...
!> \param opposite ...
! **************************************************************************************************
   SUBROUTINE exchange_opposite(tensor, opposite)
      REAL(dp), INTENT(IN)                               :: tensor(naux, nao, nao)
      REAL(dp), INTENT(OUT)                              :: opposite(naux, nao, nao)

      INTEGER                                            :: kap, lam, mu, nu, p, q

      opposite = 0.0_dp
      DO nu = 1, nao
         DO mu = 1, nao
            DO lam = 1, nao
               DO kap = 1, nao
                  DO q = 1, naux
                     DO p = 1, naux
                        opposite(p, mu, nu) = opposite(p, mu, nu) + &
                                              inverse(p, q)*density(mu, kap)*density(nu, lam)*tensor(q, kap, lam)
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE exchange_opposite

! **************************************************************************************************
!> \brief Geometry energy uses only closed-form C and analytic B, with no production integral calls.
!> \param coords ...
!> \return ...
! **************************************************************************************************
   FUNCTION geometry_energy(coords) RESULT(energy)
      REAL(dp), INTENT(IN)                               :: coords(3, natom)
      REAL(dp)                                           :: energy

      REAL(dp)                                           :: b(nao, ncenter), &
                                                            c(naux, ncenter, ncenter), &
                                                            tensor(naux, nao, nao)

      CALL closed_data(coords, b, c)
      CALL dense_transition(b, c, tensor)
      energy = exchange_energy(smooth + tensor)
   END FUNCTION geometry_energy

! **************************************************************************************************
!> \brief Unique-owner rectangular B matrix on the current tall or wide matrix process grid.
!> \param dense ...
!> \param matrix ...
! **************************************************************************************************
   SUBROUTINE make_overlap(dense, matrix)
      REAL(dp), INTENT(IN)                               :: dense(nao, ncenter)
      TYPE(dbcsr_type), INTENT(OUT)                      :: matrix

      INTEGER                                            :: a, atom, first, owner

      row_sizes => ao_sizes
      col_sizes => pp_sizes
      CALL dbcsr_create(matrix, 'Native HFX B', distribution, 'N', row_sizes, col_sizes)
      DO atom = 1, natom
         first = SUM(ao_sizes(:atom - 1)) + 1
         DO a = 1, ncenter
            CALL dbcsr_get_stored_coordinates(matrix, atom, a, owner)
            IF (owner == rank) CALL dbcsr_put_block(matrix, atom, a, &
                                                    dense(first:first + ao_sizes(atom) - 1, a:a))
         END DO
      END DO
      CALL dbcsr_finalize(matrix)
   END SUBROUTINE make_overlap

! **************************************************************************************************
!> \brief Synthetic unweighted SAP integrals preserve the original AO-to-projector displacement.
!> \param b ...
!> \param matrix ...
!> \param sap ...
! **************************************************************************************************
   SUBROUTINE make_sap(b, matrix, sap)
      REAL(dp), INTENT(IN)                               :: b(nao, ncenter)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      TYPE(sap_int_type), POINTER                        :: sap(:)

      INTEGER                                            :: a, atom, axis, first, mu, myrow
      REAL(dp)                                           :: distance(3)

      ALLOCATE (sap(1))
      ALLOCATE (sap(1)%alist(natom))
      CALL dbcsr_get_info(matrix, my_prow=myrow)
      DO atom = 1, natom
         IF (ao_rows(atom) /= myrow) CYCLE
         first = SUM(ao_sizes(:atom - 1)) + 1
         sap(1)%alist(atom)%aatom = atom
         ALLOCATE (sap(1)%alist(atom)%clist(ncenter))
         DO a = 1, ncenter
            distance = positions(:, center_atom(a)) - positions(:, atom)
            sap(1)%alist(atom)%clist(a)%catom = center_atom(a)
            sap(1)%alist(atom)%clist(a)%cell = 0
            sap(1)%alist(atom)%clist(a)%rac = distance
            ALLOCATE (sap(1)%alist(atom)%clist(a)%acint(ao_sizes(atom), 1, 4))
            sap(1)%alist(atom)%clist(a)%acint(:, 1, 1) = b(first:first + ao_sizes(atom) - 1, a)
            DO axis = 1, 3
               DO mu = first, first + ao_sizes(atom) - 1
                  sap(1)%alist(atom)%clist(a)%acint(mu - first + 1, 1, axis + 1) = &
                     2.0_dp*gamma(mu)*anisotropy(axis)*distance(axis)*b(mu, a)
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE make_sap

! **************************************************************************************************
!> \brief Fill arbitrary blocked tensor layouts without replicated blocks.
!> \param tensor ...
!> \param dense ...
! **************************************************************************************************
   SUBROUTINE put_tensor(tensor, dense)
      TYPE(dbt_type), INTENT(INOUT)                      :: tensor
      REAL(dp), INTENT(IN)                               :: dense(:, :, :)

      INTEGER                                            :: first(3), i, ind(3), j, k, last(3), &
                                                            nblks(3), owner
      INTEGER, ALLOCATABLE                               :: s1(:), s2(:), s3(:)

      CALL dbt_get_info(tensor, nblks_total=nblks)
      ALLOCATE (s1(nblks(1)), s2(nblks(2)), s3(nblks(3)))
      CALL dbt_get_info(tensor, blk_size_1=s1, blk_size_2=s2, blk_size_3=s3)
      DO k = 1, nblks(3)
         DO j = 1, nblks(2)
            DO i = 1, nblks(1)
               ind = [i, j, k]
               first = [SUM(s1(:i - 1)), SUM(s2(:j - 1)), SUM(s3(:k - 1))] + 1
               last = first + [s1(i), s2(j), s3(k)] - 1
               CALL dbt_get_stored_coordinates(tensor, ind, owner)
               IF (owner /= tensor%pgrid%mp_comm_2d%mepos) CYCLE
               IF (ALL(dense(first(1):last(1), first(2):last(2), first(3):last(3)) == 0.0_dp)) CYCLE
               CALL dbt_put_block(tensor, ind, [s1(i), s2(j), s3(k)], &
                                  dense(first(1):last(1), first(2):last(2), first(3):last(3)))
            END DO
         END DO
      END DO
   END SUBROUTINE put_tensor

! **************************************************************************************************
!> \brief Gather local tensor blocks once, independent of atomic or functionwise blocking.
!> \param tensor ...
!> \param dense ...
! **************************************************************************************************
   SUBROUTINE gather_tensor(tensor, dense)
      TYPE(dbt_type), INTENT(INOUT)                      :: tensor
      REAL(dp), INTENT(OUT)                              :: dense(:, :, :)

      INTEGER                                            :: first(3), ind(3), last(3), nblks(3)
      INTEGER, ALLOCATABLE                               :: s1(:), s2(:), s3(:)
      LOGICAL                                            :: found
      REAL(dp), ALLOCATABLE                              :: block(:, :, :)
      TYPE(dbt_iterator_type)                            :: iter

      CALL dbt_get_info(tensor, nblks_total=nblks)
      ALLOCATE (s1(nblks(1)), s2(nblks(2)), s3(nblks(3)))
      CALL dbt_get_info(tensor, blk_size_1=s1, blk_size_2=s2, blk_size_3=s3)
      dense = 0.0_dp
      CALL dbt_iterator_start(iter, tensor)
      DO WHILE (dbt_iterator_blocks_left(iter))
         CALL dbt_iterator_next_block(iter, ind)
         CALL dbt_get_block(tensor, ind, block, found)
         IF (.NOT. found) ERROR STOP 'Native HFX tensor iterator lost block'
         first = [SUM(s1(:ind(1) - 1)), SUM(s2(:ind(2) - 1)), SUM(s3(:ind(3) - 1))] + 1
         last = first + SHAPE(block) - 1
         dense(first(1):last(1), first(2):last(2), first(3):last(3)) = block
         DEALLOCATE (block)
      END DO
      CALL dbt_iterator_stop(iter)
      CALL comm%sum(dense)
   END SUBROUTINE gather_tensor

! **************************************************************************************************
!> \brief Gather rectangular B gradients on their original matrix block distribution.
!> \param matrix ...
!> \param dense ...
! **************************************************************************************************
   SUBROUTINE gather_matrix(matrix, dense)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      REAL(dp), INTENT(OUT)                              :: dense(:, :)

      INTEGER                                            :: col, first_col, first_row, nc, nr, row
      INTEGER, POINTER                                   :: cols(:), rows(:)
      LOGICAL                                            :: transposed
      REAL(dp), POINTER                                  :: block(:, :)
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL dbcsr_get_info(matrix, row_blk_offset=rows, col_blk_offset=cols)
      dense = 0.0_dp
      CALL dbcsr_iterator_start(iter, matrix)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, row, col, block, transposed=transposed)
         first_row = rows(row)
         first_col = cols(col)
         nr = SIZE(block, 1)
         nc = SIZE(block, 2)
         IF (transposed) THEN
            dense(first_row:first_row + nc - 1, first_col:first_col + nr - 1) = TRANSPOSE(block)
         ELSE
            dense(first_row:first_row + nr - 1, first_col:first_col + nc - 1) = block
         END IF
      END DO
      CALL dbcsr_iterator_stop(iter)
      CALL comm%sum(dense)
   END SUBROUTINE gather_matrix

! **************************************************************************************************
!> \brief Fail on nonfinite values as well as excessive absolute error.
!> \param error ...
!> \param label ...
!> \param tolerance ...
! **************************************************************************************************
   SUBROUTINE check(error, label, tolerance)
      REAL(dp), INTENT(IN)                               :: error
      CHARACTER(len=*), INTENT(IN)                       :: label
      REAL(dp), INTENT(IN), OPTIONAL                     :: tolerance

      REAL(dp)                                           :: tol

      tol = 2.e-12_dp
      IF (PRESENT(tolerance)) tol = tolerance
      max_error = MAX(max_error, error)
      IF (.NOT. ieee_is_finite(error) .OR. error > tol) THEN
         WRITE (*, *) rank, TRIM(label), error, tol
         ERROR STOP 'Native RI-HFX tensor or gradient mismatch'
      END IF
   END SUBROUTINE check
END PROGRAM qs_uspp_hfx_unittest
