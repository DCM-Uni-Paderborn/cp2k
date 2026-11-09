!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
!> \brief Original-data auxiliary and mixed valence charge metrics for native projection.
MODULE qs_upf_inverse_charge_metric
   USE admm_types,                      ONLY: admm_type
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE cp_control_types,                ONLY: dft_control_type
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: fourpi
   USE mathlib,                         ONLY: invmat
   USE particle_types,                  ONLY: particle_type
   USE paw_basis_types,                 ONLY: get_paw_basis_info
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_grid_atom,                    ONLY: grid_atom_type
   USE qs_harmonics_atom,               ONLY: harmonics_atom_type
   USE qs_kind_types,                   ONLY: get_qs_kind
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_upf_inverse,                  ONLY: qs_upf_inverse_response
   USE qs_upf_inverse_background,       ONLY: qs_upf_inverse_background_collocate,&
                                              qs_upf_inverse_background_geometry,&
                                              qs_upf_inverse_background_integrate,&
                                              qs_upf_inverse_background_prepare
   USE qs_upf_inverse_density,          ONLY: qs_upf_inverse_density_pair_integrate,&
                                              qs_upf_inverse_density_radial,&
                                              qs_upf_inverse_density_radial_adjoint,&
                                              qs_upf_inverse_gapw_grid,&
                                              qs_upf_inverse_gapw_view_type
   USE qs_upf_inverse_gaussian,         ONLY: qs_upf_inverse_gaussian_integrate,&
                                              qs_upf_inverse_gaussian_radial
   USE qs_upf_inverse_gaussian_centers, ONLY: qs_upf_inverse_gaussian_centers_prepare,&
                                              qs_upf_inverse_gaussian_centers_release,&
                                              qs_upf_inverse_gaussian_centers_type
   USE qs_upf_inverse_gaussian_pairs,   ONLY: qs_upf_inverse_gaussian_pairs_density,&
                                              qs_upf_inverse_gaussian_pairs_integrate,&
                                              qs_upf_inverse_gaussian_pairs_prepare
   USE qs_upf_inverse_grid,             ONLY: qs_upf_inverse_grid_metric_dense
   USE qs_upf_inverse_types,            ONLY: qs_upf_inverse_augmentation_type,&
                                              qs_upf_inverse_charge_metric_type,&
                                              qs_upf_inverse_grid_type,&
                                              qs_upf_inverse_scf_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: qs_upf_inverse_charge_metric_prepare
   PUBLIC :: qs_upf_inverse_charge_metric_contract

CONTAINS

! **************************************************************************************************
!> \brief Build the auxiliary self metric or the ordered AUX_FIT/ORB charge metric.
!>        Match each representation's density measure. GPW/GAPW_XC inverse fields retain
!>        their defining grid measure. Full GAPW uses native soft Gaussian grid fields,
!>        independent native Gaussian center products and the complete home inverse/Q residual.
!>        Neither a fitted common Gaussian basis nor a point-rule replacement of GAPW is used.
!> \param qs_env Native geometry, Gaussian grids and auxiliary owners.
!> \param row Prepared auxiliary inverse basis. It need not own an SCF.
!> \param column Same auxiliary state for the self block, or the original orbital state.
!> \param map New owned metric. Release before rebuilding either native basis or its tasks.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_charge_metric_prepare(qs_env, row, column, map)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_upf_inverse_scf_type), POINTER             :: row, column
      TYPE(qs_upf_inverse_charge_metric_type), &
         INTENT(OUT)                                     :: map

      COMPLEX(dp), ALLOCATABLE                           :: part(:, :)
      INTEGER                                            :: a, b, k, kind, ncenter, nl, nr
      LOGICAL                                            :: bloch
      TYPE(admm_type), POINTER                           :: admm
      TYPE(dft_control_type), POINTER                    :: control
      TYPE(particle_type), POINTER                       :: particles(:)
      TYPE(pw_r3d_rs_type)                               :: one
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_upf_inverse_augmentation_type), POINTER    :: augmentation
      TYPE(qs_upf_inverse_gaussian_centers_type)         :: column_centers, row_centers
      TYPE(qs_upf_inverse_grid_type), POINTER            :: left_frame, right_frame

      CPASSERT(ASSOCIATED(row) .AND. ASSOCIATED(column))
      CPASSERT(row%basis_type == 'AUX_FIT')
      CPASSERT(column%basis_type == 'AUX_FIT' .OR. column%basis_type == 'ORB')
      map%ordered = column%basis_type == 'ORB'
      IF (.NOT. map%ordered) THEN
         CPASSERT(ASSOCIATED(row, column))
      END IF
      CALL get_qs_env(qs_env, ks_env=ks_env, admm_env=admm, dft_control=control, particle_set=particles)
      CPASSERT(ASSOCIATED(admm))
      map%atomic = control%qs_control%gapw
      map%native_gaussian = row%native_gaussian
      CPASSERT(map%native_gaussian .EQV. column%native_gaussian)
      map%grid => row%source_grid
      CPASSERT(ASSOCIATED(map%grid) .AND. ASSOCIATED(column%source_grid))
      bloch = ASSOCIATED(row%bloch%kpoints)
      CPASSERT(bloch .EQV. ASSOCIATED(column%bloch%kpoints))
      IF (bloch) THEN
         CPASSERT(ASSOCIATED(row%bloch%kpoints, column%bloch%kpoints))
         CPASSERT(SIZE(row%bloch%weights) == SIZE(column%bloch%weights))
         CPASSERT(ALL(row%bloch%weights == column%bloch%weights))
         CPASSERT(ALL(row%bloch%reduced_k == column%bloch%reduced_k))
         ALLOCATE (map%weights, SOURCE=row%bloch%weights)
         ALLOCATE (map%reduced_k, SOURCE=row%bloch%reduced_k)
      ELSE
         ALLOCATE (map%weights(1), map%reduced_k(3, 1))
         map%weights = 1.0_dp
         map%reduced_k = 0.0_dp
      END IF
      CALL select_frames(row, column, 1, left_frame, right_frame, augmentation)
      nl = SIZE(left_frame%ao%coefficients, 2)
      nr = SIZE(right_frame%ao%coefficients, 2)
      ALLOCATE (map%values(nl, nr, SIZE(map%weights)))
      map%values = CMPLX(0.0_dp, 0.0_dp, dp)
      IF (map%native_gaussian .OR. map%atomic) THEN
         IF (map%ordered) THEN
            CALL qs_upf_inverse_background_prepare(ks_env, map%gaussian, soft_basis=map%atomic, qs_env=qs_env, &
                                                   basis_type='AUX_FIT', right_basis_type='ORB')
         ELSE
            CALL qs_upf_inverse_background_prepare(ks_env, map%gaussian, soft_basis=map%atomic, qs_env=qs_env, &
                                                   basis_type='AUX_FIT')
         END IF
         CPASSERT(ALL(map%gaussian%weights == map%weights))
         CPASSERT(ALL(map%gaussian%reduced_k == map%reduced_k))
         CALL one%create(map%grid)
         one%array = 1.0_dp
         CALL qs_upf_inverse_background_integrate(map%gaussian, qs_env, one, map%values)
         CALL one%release()
      END IF
      IF (.NOT. map%atomic) THEN
         DO k = 1, SIZE(map%weights)
            CALL select_frames(row, column, k, left_frame, right_frame, augmentation)
            IF (map%ordered) THEN
               CALL qs_upf_inverse_grid_metric_dense(left_frame, left_frame%ao%radial_q, part, &
                                                     orbital_fields=.NOT. map%native_gaussian, right_frame=right_frame)
            ELSE
               CALL qs_upf_inverse_grid_metric_dense(left_frame, left_frame%ao%radial_q, part, &
                                                     orbital_fields=.NOT. map%native_gaussian)
            END IF
            map%values(:, :, k) = map%values(:, :, k) + part
         END DO
         map%ready = .TRUE.
         RETURN
      END IF
      CPASSERT(ASSOCIATED(admm%admm_gapw_env))
      map%row_kinds => admm%admm_gapw_env%admm_kind_set
      map%column_kinds => map%row_kinds
      IF (map%ordered) CALL get_qs_env(qs_env, qs_kind_set=map%column_kinds)
      CALL qs_upf_inverse_gaussian_centers_prepare(qs_env, map%row_kinds, map%column_kinds, row_centers, column_centers)
      ncenter = 0
      DO a = 1, SIZE(particles)
         kind = particles(a)%atomic_kind%kind_number
         IF (map%row_kinds(kind)%paw_atom .OR. map%column_kinds(kind)%paw_atom .OR. &
             ANY(left_frame%ao%center_atoms == a)) ncenter = ncenter + 1
      END DO
      ALLOCATE (map%pairs(ncenter), map%center_atoms(ncenter), map%center_kinds(ncenter), map%inverse_centers(ncenter))
      ncenter = 0
      DO a = 1, SIZE(particles)
         kind = particles(a)%atomic_kind%kind_number
         IF (.NOT. (map%row_kinds(kind)%paw_atom .OR. map%column_kinds(kind)%paw_atom .OR. &
                    ANY(left_frame%ao%center_atoms == a))) CYCLE
         ncenter = ncenter + 1
         map%center_atoms(ncenter) = a
         map%center_kinds(ncenter) = kind
         map%inverse_centers(ncenter) = 0
         DO b = 1, SIZE(left_frame%ao%center_atoms)
            IF (left_frame%ao%center_atoms(b) == a) map%inverse_centers(ncenter) = b
         END DO
         IF (map%ordered) THEN
            CALL qs_upf_inverse_gaussian_pairs_prepare(row_centers%oce, map%gaussian%neighbors, row_centers%kinds, &
                                    a, kind, left_frame%ao%ao_first, left_frame%ao%ao_last, map%gaussian%comm, map%pairs(ncenter), &
                      left_frame%ao%positions, left_frame%ao%cell, right_oce=column_centers%oce, right_kinds=column_centers%kinds, &
                                                       right_ao_first=right_frame%ao%ao_first, right_ao_last=right_frame%ao%ao_last)
         ELSE
            CALL qs_upf_inverse_gaussian_pairs_prepare(row_centers%oce, map%gaussian%neighbors, row_centers%kinds, &
                                    a, kind, left_frame%ao%ao_first, left_frame%ao%ao_last, map%gaussian%comm, map%pairs(ncenter), &
                                                       left_frame%ao%positions, left_frame%ao%cell)
         END IF
      END DO
      CALL qs_upf_inverse_gaussian_centers_release(row_centers)
      CALL qs_upf_inverse_gaussian_centers_release(column_centers)
      DO a = 1, SIZE(map%center_atoms)
         DO k = 1, SIZE(map%weights)
            CALL select_frames(row, column, k, left_frame, right_frame, augmentation)
            CALL center_metric(map, a, k, left_frame, right_frame, augmentation)
         END DO
      END DO
      map%ready = .TRUE.
   END SUBROUTINE qs_upf_inverse_charge_metric_prepare

! **************************************************************************************************
!> \brief Contract the complete metric derivative at fixed reduced k and original AO coefficients.
!>        Preserve the defining grid measure and the native GAPW Gaussian and atomic maps.
!>        Native grid virials already include affine atom motion; forces contract only the
!>        remaining nonaffine motion. The smooth-grid volume response is added exactly once.
!> \param qs_env Native geometry and basis owners, unchanged since preparing map.
!> \param row Original auxiliary inverse context.
!> \param column Original auxiliary or orbital inverse context.
!> \param map Prepared metric with owned source tasks and Gaussian center coefficients.
!> \param matrix Original-coordinate complex Frobenius covector, without k weights.
!> \param atom_response Total Cartesian atomic motion, including any affine cell motion.
!> \param cell_response Cartesian lattice-vector response at fixed reduced k.
!> \param response Real sum_k w_k Re sum(conjg(matrix_k)*delta(metric_k)).
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_charge_metric_contract(qs_env, row, column, map, matrix, atom_response, cell_response, response)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_upf_inverse_scf_type), POINTER             :: row, column
      TYPE(qs_upf_inverse_charge_metric_type), &
         INTENT(INOUT)                                   :: map
      COMPLEX(dp), INTENT(IN)                            :: matrix(:, :, :)
      REAL(dp), INTENT(IN)                               :: atom_response(:, :), cell_response(3, 3)
      REAL(dp), INTENT(OUT)                              :: response

      COMPLEX(dp), ALLOCATABLE                           :: derivative(:, :), dw(:, :), part(:, :), &
                                                            right_dw(:, :)
      INTEGER                                            :: a, info, k
      REAL(dp)                                           :: center_response, inverse_cell(3, 3), &
                                                            smooth_charge, strain(3, 3), &
                                                            trace_strain, virial(3, 3)
      REAL(dp), ALLOCATABLE                              :: forces(:, :), point_response(:, :), &
                                                            right_point_response(:, :)
      TYPE(pw_c1d_gs_type)                               :: rho_g
      TYPE(pw_r3d_rs_type)                               :: one, rho
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_upf_inverse_augmentation_type), POINTER    :: augmentation
      TYPE(qs_upf_inverse_grid_type), POINTER            :: left_frame, right_frame

      CPASSERT(map%ready)
      CPASSERT(ALL(SHAPE(matrix) == SHAPE(map%values)))
      CALL select_frames(row, column, 1, left_frame, right_frame, augmentation)
      IF (.NOT. ALL(SHAPE(atom_response) == SHAPE(left_frame%ao%positions))) THEN
         CPABORT("Invalid charge metric atom response")
      END IF
      inverse_cell = left_frame%ao%cell
      CALL invmat(inverse_cell, info)
      CPASSERT(info == 0)
      strain = MATMUL(cell_response, inverse_cell)
      trace_strain = strain(1, 1) + strain(2, 2) + strain(3, 3)
      response = 0.0_dp
      IF (map%native_gaussian .OR. map%atomic) THEN
         CALL get_qs_env(qs_env, ks_env=ks_env)
         CALL one%create(map%grid)
         CALL rho%create(map%grid)
         CALL rho_g%create(map%grid)
         one%array = 1.0_dp
         CALL qs_upf_inverse_background_collocate(map%gaussian, ks_env, matrix, rho, rho_g, total_rho=smooth_charge)
         ALLOCATE (forces, MOLD=atom_response)
         CALL qs_upf_inverse_background_geometry(map%gaussian, qs_env, matrix, one, forces, virial)
         ! Native total_rho is the negative electronic density integral.
         response = -trace_strain*smooth_charge - SUM(virial*strain) - &
                    SUM(forces*(atom_response - MATMUL(strain, left_frame%ao%positions)))
         CALL one%release()
         CALL rho%release()
         CALL rho_g%release()
      END IF
      DO k = 1, SIZE(map%weights)
         CALL select_frames(row, column, k, left_frame, right_frame, augmentation)
         ALLOCATE (point_response, SOURCE=MATMUL(strain, left_frame%ao%points))
         IF (.NOT. map%atomic) THEN
            IF (map%ordered) THEN
               CALL qs_upf_inverse_grid_metric_dense(left_frame, left_frame%ao%radial_q, part, &
                                           orbital_fields=.NOT. map%native_gaussian, right_frame=right_frame, response=derivative, &
                            atom_response=atom_response, cell_response=cell_response, reduced_k_response=[0.0_dp, 0.0_dp, 0.0_dp], &
                                                 point_response=point_response, weights_response=trace_strain*left_frame%ao%weights)
            ELSE
               CALL qs_upf_inverse_grid_metric_dense(left_frame, left_frame%ao%radial_q, part, &
                                                     orbital_fields=.NOT. map%native_gaussian, response=derivative, &
                            atom_response=atom_response, cell_response=cell_response, reduced_k_response=[0.0_dp, 0.0_dp, 0.0_dp], &
                                                 point_response=point_response, weights_response=trace_strain*left_frame%ao%weights)
            END IF
            response = response + map%weights(k)*REAL(SUM(CONJG(matrix(:, :, k))*derivative), dp)
         ELSE
            ALLOCATE (dw, MOLD=left_frame%ao%coefficients)
            ALLOCATE (right_dw, MOLD=right_frame%ao%coefficients)
            ALLOCATE (right_point_response, SOURCE=MATMUL(strain, right_frame%ao%points))
            CALL qs_upf_inverse_response(left_frame%ao, point_response, atom_response, cell_response, &
                                         [0.0_dp, 0.0_dp, 0.0_dp], trace_strain*left_frame%ao%weights, dw)
            CALL qs_upf_inverse_response(right_frame%ao, right_point_response, atom_response, cell_response, &
                                         [0.0_dp, 0.0_dp, 0.0_dp], trace_strain*right_frame%ao%weights, right_dw)
            DO a = 1, SIZE(map%center_atoms)
               CALL center_contract(map, a, k, left_frame, right_frame, augmentation, matrix(:, :, k), &
                                    atom_response, cell_response, dw, right_dw, center_response)
               response = response + map%weights(k)*center_response
            END DO
            DEALLOCATE (dw, right_dw, right_point_response)
         END IF
         DEALLOCATE (point_response)
      END DO
   END SUBROUTINE qs_upf_inverse_charge_metric_contract

! **************************************************************************************************
!> \brief Directional native Gaussian and inverse/Q center contraction on fixed atomic radial nodes.
!> \param map ...
!> \param center ...
!> \param k ...
!> \param row ...
!> \param column ...
!> \param augmentation ...
!> \param matrix ...
!> \param atom_response ...
!> \param cell_response ...
!> \param dw ...
!> \param right_dw ...
!> \param response ...
! **************************************************************************************************
   SUBROUTINE center_contract(map, center, k, row, column, augmentation, matrix, atom_response, cell_response, &
                              dw, right_dw, response)
      TYPE(qs_upf_inverse_charge_metric_type), &
         INTENT(IN)                                      :: map
      INTEGER, INTENT(IN)                                :: center, k
      TYPE(qs_upf_inverse_grid_type), INTENT(IN)         :: row, column
      TYPE(qs_upf_inverse_augmentation_type), INTENT(IN) :: augmentation
      COMPLEX(dp), INTENT(IN)                            :: matrix(:, :)
      REAL(dp), INTENT(IN)                               :: atom_response(:, :), cell_response(3, 3)
      COMPLEX(dp), INTENT(IN)                            :: dw(:, :), right_dw(:, :)
      REAL(dp), INTENT(OUT)                              :: response

      COMPLEX(dp), ALLOCATABLE                           :: derivative(:, :), part(:, :)
      INTEGER                                            :: atom, kind, nh, nl, np, nr
      REAL(dp), ALLOCATABLE :: dh(:, :), drh(:, :), drs(:, :), ds(:, :), h(:, :), &
         point_motion(:, :), point_potential(:), points(:, :), potential(:, :), radial(:, :), &
         rh(:, :), rs(:, :), s(:, :), weights(:), zero(:)
      TYPE(grid_atom_type), POINTER                      :: grid
      TYPE(gto_basis_set_type), POINTER                  :: basis
      TYPE(harmonics_atom_type), POINTER                 :: harmonics
      TYPE(qs_upf_inverse_gapw_view_type), TARGET         :: view

      atom = map%center_atoms(center)
      kind = map%center_kinds(center)
      CALL qs_upf_inverse_gapw_grid(row%ao, map%inverse_centers(center), map%row_kinds(kind), &
                                    grid, harmonics, nl, points, view, atom)
      CALL get_qs_kind(map%column_kinds(kind), basis_set=basis, basis_type='GAPW_1C')
      CALL get_paw_basis_info(basis, nsatbas=nr)
      nh = harmonics%max_iso_not0
      ALLOCATE (h(nl, nr), s(nl, nr), dh(nl, nr), ds(nl, nr), rh(grid%nr, nh), rs(grid%nr, nh), &
                drh(grid%nr, nh), drs(grid%nr, nh), radial(grid%nr, nh), zero(grid%nr))
      zero = 0.0_dp
      CALL qs_upf_inverse_gaussian_pairs_density(map%pairs(center), map%reduced_k(:, k), matrix, h, s, dh, ds, &
                                                 atom_response, cell_response, [0.0_dp, 0.0_dp, 0.0_dp])
      IF (map%ordered) THEN
         CALL qs_upf_inverse_gaussian_radial(map%row_kinds(kind), h, rh, drh, dh, zero, map%column_kinds(kind))
         CALL qs_upf_inverse_gaussian_radial(map%row_kinds(kind), s, rs, drs, ds, zero, map%column_kinds(kind))
      ELSE
         CALL qs_upf_inverse_gaussian_radial(map%row_kinds(kind), h, rh, drh, dh, zero)
         CALL qs_upf_inverse_gaussian_radial(map%row_kinds(kind), s, rs, drs, ds, zero)
      END IF
      np = grid%nr*grid%ng_sphere
      DEALLOCATE (zero)
      ALLOCATE (zero(np), point_potential(np), potential(grid%nr, nh))
      zero = 0.0_dp
      CALL qs_upf_inverse_density_radial(grid, harmonics, zero, radial, drh, drs)
      response = SQRT(fourpi)*SUM(grid%wr*radial(:, 1))
      IF (map%inverse_centers(center) == 0) RETURN
      potential = 0.0_dp
      potential(:, 1) = SQRT(fourpi)
      CALL qs_upf_inverse_density_radial_adjoint(grid, harmonics, potential, point_potential, rh, rs)
      ALLOCATE (point_motion(3, np), weights(np), part(SIZE(matrix, 1), SIZE(matrix, 2)), &
                derivative(SIZE(matrix, 1), SIZE(matrix, 2)))
      weights(:) = RESHAPE(grid%weight, [np])
      point_motion(:, :) = SPREAD(atom_response(:, atom), 2, np)
      CALL qs_upf_inverse_density_pair_integrate(row%ao, column%ao, augmentation, map%inverse_centers(center), &
                                                 points, weights, point_potential, part, derivative, point_motion, &
                                                 atom_response, cell_response, [0.0_dp, 0.0_dp, 0.0_dp], dw, right_dw, zero, zero)
      response = response + REAL(SUM(CONJG(matrix)*derivative), dp)
   END SUBROUTINE center_contract

! **************************************************************************************************
!> \brief Add a native center on its single chosen row quadrature, including the finite angular map.
!> \param map ...
!> \param center ...
!> \param k ...
!> \param row ...
!> \param column ...
!> \param augmentation ...
! **************************************************************************************************
   SUBROUTINE center_metric(map, center, k, row, column, augmentation)
      TYPE(qs_upf_inverse_charge_metric_type), &
         INTENT(INOUT)                                   :: map
      INTEGER, INTENT(IN)                                :: center, k
      TYPE(qs_upf_inverse_grid_type), INTENT(IN)         :: row, column
      TYPE(qs_upf_inverse_augmentation_type), INTENT(IN) :: augmentation

      COMPLEX(dp), ALLOCATABLE                           :: part(:, :)
      INTEGER                                            :: atom, kind, nh, nl, np, nr
      REAL(dp), ALLOCATABLE                              :: bh(:, :), bs(:, :), hard(:, :), &
                                                            point_potential(:), points(:, :), &
                                                            potential(:, :), soft(:, :), weights(:)
      TYPE(grid_atom_type), POINTER                      :: grid
      TYPE(gto_basis_set_type), POINTER                  :: basis
      TYPE(harmonics_atom_type), POINTER                 :: harmonics
      TYPE(qs_upf_inverse_gapw_view_type), TARGET         :: view

      atom = map%center_atoms(center)
      kind = map%center_kinds(center)
      CALL qs_upf_inverse_gapw_grid(row%ao, map%inverse_centers(center), map%row_kinds(kind), &
                                    grid, harmonics, nl, points, view, atom)
      CALL get_qs_kind(map%column_kinds(kind), basis_set=basis, basis_type='GAPW_1C')
      CALL get_paw_basis_info(basis, nsatbas=nr)
      np = grid%nr*grid%ng_sphere
      nh = harmonics%max_iso_not0
      ALLOCATE (potential(grid%nr, nh), hard(grid%nr, nh), soft(grid%nr, nh), &
                point_potential(np), weights(np), bh(nl, nr), bs(nl, nr), &
                part(SIZE(map%values, 1), SIZE(map%values, 2)))
      potential = 0.0_dp
      potential(:, 1) = SQRT(fourpi)
      CALL qs_upf_inverse_density_radial_adjoint(grid, harmonics, potential, point_potential, hard, soft)
      IF (map%ordered) THEN
         CALL qs_upf_inverse_gaussian_integrate(map%row_kinds(kind), hard, bh, right_kind=map%column_kinds(kind))
         CALL qs_upf_inverse_gaussian_integrate(map%row_kinds(kind), soft, bs, right_kind=map%column_kinds(kind))
      ELSE
         CALL qs_upf_inverse_gaussian_integrate(map%row_kinds(kind), hard, bh)
         CALL qs_upf_inverse_gaussian_integrate(map%row_kinds(kind), soft, bs)
      END IF
      CALL qs_upf_inverse_gaussian_pairs_integrate(map%pairs(center), map%reduced_k(:, k), bh, bs, part)
      map%values(:, :, k) = map%values(:, :, k) + part
      IF (map%inverse_centers(center) == 0) RETURN
      weights(:) = RESHAPE(grid%weight, [np])
      CALL qs_upf_inverse_density_pair_integrate(row%ao, column%ao, augmentation, map%inverse_centers(center), &
                                                 points, weights, point_potential, part)
      map%values(:, :, k) = map%values(:, :, k) + part
   END SUBROUTINE center_metric

! **************************************************************************************************
!> \brief Select matching physical Bloch contributions without transferring SCF ownership.
!> \param row ...
!> \param column ...
!> \param k ...
!> \param left ...
!> \param right ...
!> \param augmentation ...
! **************************************************************************************************
   SUBROUTINE select_frames(row, column, k, left, right, augmentation)
      TYPE(qs_upf_inverse_scf_type), POINTER             :: row, column
      INTEGER, INTENT(IN)                                :: k
      TYPE(qs_upf_inverse_grid_type), POINTER            :: left, right
      TYPE(qs_upf_inverse_augmentation_type), POINTER    :: augmentation

      left => row%frame
      right => column%frame
      augmentation => row%augmentation
      IF (ASSOCIATED(row%bloch%kpoints)) THEN
         left => row%bloch%frames(k)
         right => column%bloch%frames(k)
         augmentation => row%bloch%augmentation
      ELSE
         CPASSERT(k == 1)
      END IF
   END SUBROUTINE select_frames

END MODULE qs_upf_inverse_charge_metric
