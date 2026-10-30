!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
!> \brief Native Gaussian grid source and its complex Bloch AO and geometry transposes.
MODULE qs_upf_inverse_background
   USE admm_types,                      ONLY: admm_type,&
                                              get_admm_env
   USE cell_types,                      ONLY: cell_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_copy, dbcsr_get_block_p, dbcsr_get_info, dbcsr_get_matrix_type, &
        dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, dbcsr_iterator_start, &
        dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_p_type, dbcsr_scale, dbcsr_set, &
        dbcsr_type_symmetric
   USE cp_dbcsr_operations,             ONLY: dbcsr_deallocate_matrix_set
   USE kinds,                           ONLY: dp
   USE kpoint_methods,                  ONLY: kp_transform_plan_create,&
                                              kpoint_star_create
   USE kpoint_mo_symmetry_methods,      ONLY: kpoint_halfturn_symmetry
   USE kpoint_types,                    ONLY: kpoint_star_type,&
                                              kpoint_type
   USE mathconstants,                   ONLY: twopi
   USE mathlib,                         ONLY: invmat
   USE particle_types,                  ONLY: particle_type
   USE pw_methods,                      ONLY: pw_transfer
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE qs_collocate_density,            ONLY: calculate_rho_elec
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_integrate_potential_product,  ONLY: integrate_v_rspace
   USE qs_kinetic,                      ONLY: build_kinetic_matrix
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type
   USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                              neighbor_list_iterate,&
                                              neighbor_list_iterator_create,&
                                              neighbor_list_iterator_p_type,&
                                              neighbor_list_iterator_release,&
                                              neighbor_list_set_p_type
   USE qs_overlap,                      ONLY: build_overlap_matrix
   USE qs_upf_inverse_types,            ONLY: qs_upf_inverse_background_release,&
                                              qs_upf_inverse_background_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: qs_upf_inverse_background_type, qs_upf_inverse_background_prepare
   PUBLIC :: qs_upf_inverse_background_release, qs_upf_inverse_background_collocate
   PUBLIC :: qs_upf_inverse_background_integrate, qs_upf_inverse_background_operator
   PUBLIC :: qs_upf_inverse_background_geometry
   PUBLIC :: qs_upf_inverse_background_gaussian_geometry
   PUBLIC :: qs_upf_inverse_background_scratch
   PUBLIC :: qs_upf_inverse_background_set_density

   INTERFACE qs_upf_inverse_background_set_density
      MODULE PROCEDURE background_images
   END INTERFACE qs_upf_inverse_background_set_density

CONTAINS

! **************************************************************************************************
!> \brief Own native image sparsity and the Bloch phase map for original-coordinate Gaussian AOs.
!>        Rebuild with native overlap lists and soft task lists after atomic or cell changes.
!> \param ks_env Live native KS environment with overlap images and Gaussian task lists.
!> \param map New owned background map. The native matrices and density state are unchanged.
!> \param soft_basis Select the GAPW soft basis. False uses the full native Gaussian basis.
!> \param qs_env Matching native owner, required with AUX_FIT. Its auxiliary lists must remain live.
!> \param basis_type ORB by default, or AUX_FIT with its own overlaps, neighbors and grid tasks.
!>        Release and rebuild the map when the native auxiliary lists change.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_background_prepare(ks_env, map, soft_basis, qs_env, basis_type)
      TYPE(qs_ks_env_type), INTENT(IN), TARGET           :: ks_env
      TYPE(qs_upf_inverse_background_type), INTENT(OUT)  :: map
      LOGICAL, INTENT(IN), OPTIONAL                      :: soft_basis
      TYPE(qs_environment_type), OPTIONAL, POINTER       :: qs_env
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL             :: basis_type

      CHARACTER(LEN=1)                                   :: symmetry
      INTEGER                                            :: a, b, ENTRY, image, info, native_cell(3)
      INTEGER, POINTER                                   :: sizes(:)
      REAL(dp)                                           :: image_error, inverse_lattice(3, 3), &
                                                            rab(3), tolerance
      TYPE(admm_type), POINTER                           :: admm
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_p_type), POINTER                        :: overlap(:, :)
      TYPE(dft_control_type), POINTER                    :: control
      TYPE(kpoint_star_type)                             :: star
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(neighbor_list_iterator_p_type), POINTER       :: iterator(:)
      TYPE(neighbor_list_set_p_type), POINTER            :: sab(:)
      TYPE(particle_type), POINTER                       :: atoms(:)
      TYPE(qs_ks_env_type), POINTER                      :: owner

      IF (PRESENT(soft_basis)) map%soft_basis = soft_basis
      IF (PRESENT(basis_type)) map%orbital_basis = basis_type
      CPASSERT(map%orbital_basis == "ORB" .OR. map%orbital_basis == "AUX_FIT")
      CALL get_ks_env(ks_env, matrix_s_kp=overlap, kpoints=kpoints, do_kpoints=map%bloch, sab_orb=sab, dft_control=control)
      IF (map%orbital_basis == "AUX_FIT") THEN
         CPASSERT(PRESENT(qs_env))
         CPASSERT(ASSOCIATED(qs_env))
         CALL get_qs_env(qs_env, admm_env=admm, ks_env=owner)
         CPASSERT(ASSOCIATED(owner, ks_env) .AND. ASSOCIATED(admm))
         CALL get_admm_env(admm, matrix_s_aux_fit_kp=overlap, sab_aux_fit=sab)
         map%neighbors => sab
         map%grid_basis = "AUX_FIT"
         map%task_list => admm%task_list_aux_fit
         IF (map%soft_basis) THEN
            CPASSERT(ASSOCIATED(admm%admm_gapw_env))
            map%grid_basis = "AUX_FIT_SOFT"
            map%task_list => admm%admm_gapw_env%task_list
         END IF
         CPASSERT(ASSOCIATED(map%task_list))
      END IF
      CPASSERT(ASSOCIATED(overlap) .AND. ASSOCIATED(sab))
      CALL dbcsr_get_info(overlap(1, 1)%matrix, row_blk_size=sizes, group=map%comm)
      map%nao = SUM(sizes)
      CPASSERT(map%nao > 0)
      ALLOCATE (map%density(SIZE(overlap, 2)), map%potential(SIZE(overlap, 2)))
      DO image = 1, SIZE(overlap, 2)
         symmetry = dbcsr_get_matrix_type(overlap(1, image)%matrix)
         CPASSERT(symmetry == dbcsr_type_symmetric)
         ALLOCATE (map%density(image)%matrix, map%potential(image)%matrix)
         CALL dbcsr_copy(map%density(image)%matrix, overlap(1, image)%matrix)
         CALL dbcsr_copy(map%potential(image)%matrix, overlap(1, image)%matrix)
         CALL dbcsr_set(map%density(image)%matrix, 0.0_dp)
         CALL dbcsr_set(map%potential(image)%matrix, 0.0_dp)
      END DO
      IF (map%bloch) THEN
         CPASSERT(ASSOCIATED(kpoints))
         IF (control%use_gauxc .AND. .NOT. kpoint_halfturn_symmetry(kpoints)) THEN
            CPABORT("Use CARTESIAN_HALF_TURNS for learned XC spatial symmetry")
         END IF
         ! The star retains physical AO actions for both spatial and time-reversed partners.
         CALL kp_transform_plan_create(map%plan, sab, kpoints%cell_to_index, SIZE(overlap, 2), &
                                       block_template=overlap(1, 1)%matrix)
         CPASSERT(map%plan%symmetric)
         CALL get_ks_env(ks_env, cell=cell, particle_set=atoms)
         inverse_lattice = cell%hmat
         CALL invmat(inverse_lattice, info)
         CPASSERT(info == 0)
         ! Native image matrices use folded atoms. Convert only their AO phase cells
         ! to the original coordinates used by the inverse fields; image ownership stays native.
         ENTRY = 0
         CALL neighbor_list_iterator_create(iterator, sab)
         DO WHILE (neighbor_list_iterate(iterator) == 0)
            CALL get_iterator_info(iterator, iatom=a, jatom=b, cell=native_cell, r=rab)
            image = kpoints%cell_to_index(native_cell(1), native_cell(2), native_cell(3))
            IF (image < 1 .OR. image > SIZE(overlap, 2)) CYCLE
            ENTRY = ENTRY + 1
            map%plan%cell(:, ENTRY) = NINT(MATMUL(inverse_lattice, atoms(a)%r - atoms(b)%r + rab))
            image_error = NORM2(atoms(a)%r - atoms(b)%r + rab - &
                                MATMUL(cell%hmat, REAL(map%plan%cell(:, ENTRY), dp)))
            CPASSERT(image_error < 1.e-9_dp)
         END DO
         CALL neighbor_list_iterator_release(iterator)
         CPASSERT(ENTRY == map%plan%nentry)
         CALL kpoint_star_create(kpoints, sizes, star, for_aux_fit=map%orbital_basis == "AUX_FIT")
         ALLOCATE (map%reduced_k, SOURCE=star%vectors)
         ALLOCATE (map%weights, SOURCE=star%weights)
      ELSE
         CPASSERT(SIZE(overlap, 2) == 1)
         ALLOCATE (map%reduced_k(3, 1), map%weights(1))
         map%reduced_k = 0.0_dp
         map%weights = 1.0_dp
      END IF
      CPASSERT(ALL(map%weights >= 0.0_dp))
      tolerance = 128.0_dp*EPSILON(1.0_dp)*SIZE(map%weights)
      CPASSERT(ABS(SUM(map%weights) - 1.0_dp) <= tolerance)
   END SUBROUTINE qs_upf_inverse_background_prepare

! **************************************************************************************************
!> \brief Copy an existing Gaussian source map into independent zeroed matrix workspaces.
!>        Retain its original image sparsity, phases, mesh weights and soft-basis convention.
!>        Inverse SCF may already have replaced native overlaps by a different charge metric.
!>        Rebuilding the source map from those overlaps would change the original task-list map.
!> \param source Existing fixed-geometry Gaussian source.
!> \param map New owned image matrices and copied phase metadata, borrowing the same communicator.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_background_scratch(source, map)
      TYPE(qs_upf_inverse_background_type), INTENT(IN)   :: source
      TYPE(qs_upf_inverse_background_type), INTENT(OUT)  :: map

      INTEGER                                            :: image

      CPASSERT(ASSOCIATED(source%density))
      CPASSERT(ASSOCIATED(source%potential))
      CPASSERT(SIZE(source%density) == SIZE(source%potential))
      CPASSERT(ALLOCATED(source%weights) .AND. ALLOCATED(source%reduced_k))
      map%nao = source%nao
      map%bloch = source%bloch
      map%soft_basis = source%soft_basis
      map%orbital_basis = source%orbital_basis
      map%grid_basis = source%grid_basis
      map%task_list => source%task_list
      map%neighbors => source%neighbors
      map%comm = source%comm
      map%plan = source%plan
      ALLOCATE (map%weights, SOURCE=source%weights)
      ALLOCATE (map%reduced_k, SOURCE=source%reduced_k)
      ALLOCATE (map%density(SIZE(source%density)), map%potential(SIZE(source%potential)))
      DO image = 1, SIZE(source%density)
         ALLOCATE (map%density(image)%matrix, map%potential(image)%matrix)
         CALL dbcsr_copy(map%density(image)%matrix, source%density(image)%matrix)
         CALL dbcsr_copy(map%potential(image)%matrix, source%potential(image)%matrix)
         CALL dbcsr_set(map%density(image)%matrix, 0.0_dp)
         CALL dbcsr_set(map%potential(image)%matrix, 0.0_dp)
      END DO
   END SUBROUTINE qs_upf_inverse_background_scratch

! **************************************************************************************************
!> \brief Collocate the Gaussian soft background from explicit accepted AO occupations.
!>        The native image density may be stale after inverse Bloch density mixing.
!>        Mesh weights and image multiplicities are applied once before native soft collocation.
!> \param map ...
!> \param ks_env Native owner of the current soft Gaussian task lists and PW grids.
!> \param matrix Replicated Hermitian occupation matrices indexed by global k, without mesh weights.
!> \param rho Overwritten positive electronic Gaussian soft density.
!> \param rho_gspace Its matching native Fourier transform.
!> \param total_rho Optional native negative electron charge.
!> \param kinetic Optional soft Gaussian kinetic density, using the native factor one half.
!> \param kinetic_gspace Matching Fourier kinetic density.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_background_collocate(map, ks_env, matrix, rho, rho_gspace, total_rho, kinetic, kinetic_gspace)
      TYPE(qs_upf_inverse_background_type), &
         INTENT(INOUT)                                   :: map
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      COMPLEX(dp), INTENT(IN)                            :: matrix(:, :, :)
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: rho
      TYPE(pw_c1d_gs_type), INTENT(INOUT)                :: rho_gspace
      REAL(dp), INTENT(OUT), OPTIONAL                    :: total_rho
      TYPE(pw_r3d_rs_type), INTENT(INOUT), OPTIONAL      :: kinetic
      TYPE(pw_c1d_gs_type), INTENT(INOUT), OPTIONAL      :: kinetic_gspace

      CPASSERT(PRESENT(kinetic) .EQV. PRESENT(kinetic_gspace))
      CALL background_images(map, matrix)
      IF (ASSOCIATED(map%task_list)) THEN
         CALL calculate_rho_elec(matrix_p_kp=map%density, rho=rho, rho_gspace=rho_gspace, &
                                 total_rho=total_rho, ks_env=ks_env, basis_type=map%grid_basis, &
                                 task_list_external=map%task_list)
      ELSE
         CALL calculate_rho_elec(matrix_p_kp=map%density, rho=rho, rho_gspace=rho_gspace, &
                                 total_rho=total_rho, ks_env=ks_env, soft_valid=map%soft_basis)
      END IF
      ! Multigrid remapping can leave unpaired Fourier components on even grids.
      ! Coulomb and XC must use the same returned real Gaussian density.
      CALL pw_transfer(rho, rho_gspace)
      IF (PRESENT(kinetic)) THEN
         CPASSERT(ASSOCIATED(kinetic%pw_grid, rho%pw_grid))
         IF (ASSOCIATED(map%task_list)) THEN
            CALL calculate_rho_elec(matrix_p_kp=map%density, rho=kinetic, rho_gspace=kinetic_gspace, &
                                    ks_env=ks_env, basis_type=map%grid_basis, compute_tau=.TRUE., &
                                    task_list_external=map%task_list)
         ELSE
            CALL calculate_rho_elec(matrix_p_kp=map%density, rho=kinetic, rho_gspace=kinetic_gspace, &
                                    ks_env=ks_env, soft_valid=map%soft_basis, compute_tau=.TRUE.)
         END IF
         CALL pw_transfer(kinetic, kinetic_gspace)
      END IF
   END SUBROUTINE qs_upf_inverse_background_collocate

! **************************************************************************************************
!> \brief Exact complex AO transpose of native Gaussian soft collocation, without mesh weights.
!>        A unit potential gives the background contribution to the complete GAPW charge metric.
!> \param map ...
!> \param qs_env Native owner of the same soft Gaussian task lists and PW grids.
!> \param potential Native grid potential, without a grid measure or a mesh weight.
!> \param ao_operators Overwritten Hermitian AO matrices indexed by global k.
!> \param kinetic Select the transpose of the soft kinetic density instead of charge.
!> \param ks_env Native KS owner instead of qs_env.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_background_integrate(map, qs_env, potential, ao_operators, kinetic, ks_env)
      TYPE(qs_upf_inverse_background_type), &
         INTENT(INOUT)                                   :: map
      TYPE(qs_environment_type), OPTIONAL, POINTER       :: qs_env
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: potential
      COMPLEX(dp), INTENT(OUT)                           :: ao_operators(:, :, :)
      LOGICAL, INTENT(IN), OPTIONAL                      :: kinetic
      TYPE(qs_ks_env_type), INTENT(IN), OPTIONAL         :: ks_env

      INTEGER                                            :: image

      CPASSERT(ASSOCIATED(map%potential))
      CPASSERT(ALL(SHAPE(ao_operators) == [map%nao, map%nao, SIZE(map%weights)]))
      DO image = 1, SIZE(map%potential)
         CALL dbcsr_set(map%potential(image)%matrix, 0.0_dp)
      END DO
      IF (ASSOCIATED(map%task_list)) THEN
         CALL integrate_v_rspace(potential, hmat_kp=map%potential, qs_env=qs_env, &
                                 calculate_forces=.FALSE., gapw=map%soft_basis, compute_tau=kinetic, &
                                 ks_env_external=ks_env, basis_type=map%grid_basis, task_list_external=map%task_list)
      ELSE
         CALL integrate_v_rspace(potential, hmat_kp=map%potential, qs_env=qs_env, &
                                 calculate_forces=.FALSE., gapw=map%soft_basis, compute_tau=kinetic, ks_env_external=ks_env)
      END IF
      ! Native integration leaves the auxiliary PW volume element to its caller.
      DO image = 1, SIZE(map%potential)
         CALL dbcsr_scale(map%potential(image)%matrix, potential%pw_grid%dvol)
      END DO
      CALL qs_upf_inverse_background_operator(map, map%potential, ao_operators)
   END SUBROUTINE qs_upf_inverse_background_integrate

! **************************************************************************************************
!> \brief Read a native image operator in the unreduced Hermitian Bloch AO representation.
!>        Reuses the native Gaussian image domains for kinetic, overlap and grid operators.
!> \param map Prepared native image and phase map.
!> \param source Native symmetric image matrices, without k-point weights.
!> \param ao_operators Overwritten physical Hermitian AO operators indexed by global k.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_background_operator(map, source, ao_operators)
      TYPE(qs_upf_inverse_background_type), INTENT(IN)   :: map
      TYPE(dbcsr_p_type), INTENT(IN)                     :: source(:)
      COMPLEX(dp), INTENT(OUT)                           :: ao_operators(:, :, :)

      COMPLEX(dp)                                        :: phase
      INTEGER                                            :: c, ENTRY, image, k, r
      LOGICAL                                            :: found, transposed
      REAL(dp), POINTER                                  :: block(:, :)
      TYPE(dbcsr_iterator_type)                          :: iterator

      CPASSERT(ASSOCIATED(map%potential))
      CPASSERT(SIZE(source) == SIZE(map%potential))
      CPASSERT(ALL(SHAPE(ao_operators) == [map%nao, map%nao, SIZE(map%weights)]))
      ao_operators = CMPLX(0.0_dp, 0.0_dp, dp)
      IF (map%bloch) THEN
         DO ENTRY = 1, map%plan%nentry
            image = map%plan%image(ENTRY)
            CALL dbcsr_get_block_p(source(image)%matrix, map%plan%row(ENTRY), map%plan%col(ENTRY), &
                                   block, found=found)
            IF (.NOT. found) CYCLE
            r = map%plan%row_offset(ENTRY)
            c = map%plan%col_offset(ENTRY)
            DO k = 1, SIZE(map%weights)
               phase = background_phase(map, ENTRY, k)
               CALL add_operator(block, phase, r, c, ao_operators(:, :, k))
            END DO
         END DO
      ELSE
         CALL dbcsr_iterator_start(iterator, source(1)%matrix)
         DO WHILE (dbcsr_iterator_blocks_left(iterator))
            CALL dbcsr_iterator_next_block(iterator, block=block, row_offset=r, col_offset=c, transposed=transposed)
            IF (transposed) THEN
               CALL add_operator(TRANSPOSE(block), CMPLX(1.0_dp, 0.0_dp, dp), r, c, ao_operators(:, :, 1))
            ELSE
               CALL add_operator(block, CMPLX(1.0_dp, 0.0_dp, dp), r, c, ao_operators(:, :, 1))
            END IF
         END DO
         CALL dbcsr_iterator_stop(iterator)
      END IF
      DO k = 1, SIZE(map%weights)
         CALL map%comm%sum(ao_operators(:, :, k))
      END DO
   END SUBROUTINE qs_upf_inverse_background_operator

! **************************************************************************************************
!> \brief Native Gaussian kinetic or signed overlap geometry adjoint at accepted AO occupations.
!>        Reuses the ordinary Gaussian integral derivatives and native neighbor/image domains.
!>        Native density, force, virial and operator owners are unchanged.
!> \param map Prepared native image map at the geometry being differentiated.
!> \param ks_env Native Gaussian bases and neighbor lists.
!> \param matrix Hermitian occupation or energy-weighted matrices, without mesh weights.
!> \param kinetic True for Tr(P T); false for the normalization contribution -Tr(W S).
!> \param gradients Replicated Cartesian energy gradients of the selected contribution.
!> \param virial Replicated native virial of the selected contribution.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_background_gaussian_geometry(map, ks_env, matrix, kinetic, gradients, virial)
      TYPE(qs_upf_inverse_background_type), &
         INTENT(INOUT)                                   :: map
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      COMPLEX(dp), INTENT(IN)                            :: matrix(:, :, :)
      LOGICAL, INTENT(IN)                                :: kinetic
      REAL(dp), INTENT(OUT)                              :: gradients(:, :), virial(3, 3)

      INTEGER                                            :: image
      TYPE(dbcsr_p_type), POINTER                        :: occupations(:, :), operators(:, :)
      TYPE(neighbor_list_set_p_type), POINTER            :: sab(:)

      CALL background_images(map, matrix)
      IF (ASSOCIATED(map%neighbors)) THEN
         sab => map%neighbors
      ELSE
         CALL get_ks_env(ks_env, sab_orb=sab)
      END IF
      CPASSERT(ASSOCIATED(sab))
      NULLIFY (operators)
      ALLOCATE (occupations(1, SIZE(map%density)))
      DO image = 1, SIZE(map%density)
         occupations(1, image)%matrix => map%density(image)%matrix
      END DO
      IF (kinetic) THEN
         CALL build_kinetic_matrix(ks_env, matrixkp_t=operators, basis_type=map%orbital_basis, sab_nl=sab, &
                                   calculate_forces=.TRUE., matrixkp_p=occupations, &
                                   force_response=gradients, virial_response=virial)
      ELSE
         CALL build_overlap_matrix(ks_env, matrixkp_s=operators, &
                                   basis_type_a=map%orbital_basis, basis_type_b=map%orbital_basis, &
                                   sab_nl=sab, calculate_forces=.TRUE., matrixkp_p=occupations, &
                                   force_response=gradients, virial_response=virial)
      END IF
      CALL map%comm%sum(gradients)
      CALL map%comm%sum(virial)
      CALL dbcsr_deallocate_matrix_set(operators)
      DEALLOCATE (occupations)
   END SUBROUTINE qs_upf_inverse_background_gaussian_geometry

! **************************************************************************************************
!> \brief Geometry adjoint of the same native Gaussian soft source at fixed AO occupations.
!>        Uses the original soft task list and its grid force/virial kernels. Native force owners
!>        are unchanged. The virial differentiates Gaussian fields on the moving fractional grid;
!>        the grid measure and other Fourier geometry factors remain the caller's responsibility.
!> \param map Prepared source and unreduced image-density map at the current geometry.
!> \param qs_env Native owner of the same soft task list and PW grids.
!> \param matrix Accepted Hermitian AO occupations, indexed by global k without mesh weights.
!> \param potential Unweighted real-space covector conjugate to the positive Gaussian soft density.
!> \param forces Replicated atomic forces from this source contribution, including the grid measure.
!> \param virial Replicated source virial, including the grid measure.
!> \param kinetic Select the soft kinetic-density response instead of charge.
!> \param ks_env Native KS owner instead of qs_env.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_background_geometry(map, qs_env, matrix, potential, forces, virial, kinetic, ks_env)
      TYPE(qs_upf_inverse_background_type), &
         INTENT(INOUT)                                   :: map
      TYPE(qs_environment_type), OPTIONAL, POINTER       :: qs_env
      COMPLEX(dp), INTENT(IN)                            :: matrix(:, :, :)
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: potential
      REAL(dp), INTENT(OUT)                              :: forces(:, :), virial(3, 3)
      LOGICAL, INTENT(IN), OPTIONAL                      :: kinetic
      TYPE(qs_ks_env_type), INTENT(IN), OPTIONAL         :: ks_env

      INTEGER                                            :: image
      TYPE(pw_r3d_rs_type)                               :: weighted

      CPASSERT(ASSOCIATED(map%potential))
      CALL background_images(map, matrix)
      DO image = 1, SIZE(map%potential)
         CALL dbcsr_set(map%potential(image)%matrix, 0.0_dp)
      END DO
      CALL weighted%create(potential%pw_grid)
      weighted%array = potential%pw_grid%dvol*potential%array
      IF (ASSOCIATED(map%task_list)) THEN
         CALL integrate_v_rspace(weighted, hmat_kp=map%potential, pmat_kp=map%density, qs_env=qs_env, &
                                 calculate_forces=.TRUE., gapw=map%soft_basis, force_response=forces, &
                                 virial_response=virial, compute_tau=kinetic, ks_env_external=ks_env, &
                                 basis_type=map%grid_basis, task_list_external=map%task_list)
      ELSE
         CALL integrate_v_rspace(weighted, hmat_kp=map%potential, pmat_kp=map%density, qs_env=qs_env, &
                                 calculate_forces=.TRUE., gapw=map%soft_basis, force_response=forces, virial_response=virial, &
                                 compute_tau=kinetic, ks_env_external=ks_env)
      END IF
      CALL map%comm%sum(forces)
      CALL map%comm%sum(virial)
      ! Native rho_elec stores the energy gradient; return physical forces from this isolated map.
      forces(:, :) = -forces
      CALL weighted%release()
   END SUBROUTINE qs_upf_inverse_background_geometry

! **************************************************************************************************
!> \brief Apply the native image storage convention to an explicitly supplied accepted density.
!> \param map ...
!> \param matrix ...
! **************************************************************************************************
   SUBROUTINE background_images(map, matrix)
      TYPE(qs_upf_inverse_background_type), &
         INTENT(INOUT)                                   :: map
      COMPLEX(dp), INTENT(IN)                            :: matrix(:, :, :)

      COMPLEX(dp)                                        :: phase
      INTEGER                                            :: c, ENTRY, image, k, nc, nr, r
      LOGICAL                                            :: found, transposed
      REAL(dp)                                           :: defect, tolerance
      REAL(dp), POINTER                                  :: block(:, :)
      TYPE(dbcsr_iterator_type)                          :: iterator

      CPASSERT(ASSOCIATED(map%density))
      CPASSERT(ALL(SHAPE(matrix) == [map%nao, map%nao, SIZE(map%weights)]))
      DO k = 1, SIZE(map%weights)
         tolerance = 128.0_dp*EPSILON(1.0_dp)*MAX(1.0_dp, MAXVAL(ABS(matrix(:, :, k))))
         defect = MAXVAL(ABS(matrix(:, :, k) - CONJG(TRANSPOSE(matrix(:, :, k)))))
         CPASSERT(defect <= tolerance)
      END DO
      DO image = 1, SIZE(map%density)
         CALL dbcsr_set(map%density(image)%matrix, 0.0_dp)
      END DO
      IF (map%bloch) THEN
         DO ENTRY = 1, map%plan%nentry
            image = map%plan%image(ENTRY)
            CALL dbcsr_get_block_p(map%density(image)%matrix, map%plan%row(ENTRY), map%plan%col(ENTRY), &
                                   block, found=found)
            IF (.NOT. found) CYCLE
            r = map%plan%row_offset(ENTRY)
            c = map%plan%col_offset(ENTRY)
            nr = SIZE(block, 1)
            nc = SIZE(block, 2)
            DO k = 1, SIZE(map%weights)
               phase = background_phase(map, ENTRY, k)
               block = block + map%weights(k)*REAL(CONJG(phase)*matrix(r:r + nr - 1, c:c + nc - 1, k), dp)
            END DO
         END DO
      ELSE
         CALL dbcsr_iterator_start(iterator, map%density(1)%matrix)
         DO WHILE (dbcsr_iterator_blocks_left(iterator))
            CALL dbcsr_iterator_next_block(iterator, block=block, row_offset=r, col_offset=c, transposed=transposed)
            nr = SIZE(block, 1)
            nc = SIZE(block, 2)
            IF (transposed) THEN
               block = TRANSPOSE(REAL(matrix(r:r + nc - 1, c:c + nr - 1, 1), dp))
            ELSE
               block = REAL(matrix(r:r + nr - 1, c:c + nc - 1, 1), dp)
            END IF
         END DO
         CALL dbcsr_iterator_stop(iterator)
      END IF
   END SUBROUTINE background_images

! **************************************************************************************************
!> \brief Native phase including the orientation of a stored symmetric image block.
!> \param map ...
!> \param ENTRY ...
!> \param k ...
!> \return ...
! **************************************************************************************************
   FUNCTION background_phase(map, ENTRY, k) RESULT(phase)
      TYPE(qs_upf_inverse_background_type), INTENT(IN)   :: map
      INTEGER, INTENT(IN)                                :: entry, k
      COMPLEX(dp)                                        :: phase

      REAL(dp)                                           :: angle

      angle = twopi*SUM(map%reduced_k(:, k)*map%plan%cell(:, ENTRY))
      phase = CMPLX(COS(angle), map%plan%symmetry_sign(ENTRY)*SIN(angle), dp)
   END FUNCTION background_phase

! **************************************************************************************************
!> \brief Hermitian transpose under the native symmetric-block density pairing.
!> \param block ...
!> \param phase ...
!> \param r ...
!> \param c ...
!> \param ao_operator ...
! **************************************************************************************************
   SUBROUTINE add_operator(block, phase, r, c, ao_operator)
      REAL(dp), INTENT(IN)                               :: block(:, :)
      COMPLEX(dp), INTENT(IN)                            :: phase
      INTEGER, INTENT(IN)                                :: r, c
      COMPLEX(dp), INTENT(INOUT)                         :: ao_operator(:, :)

      INTEGER                                            :: i, j
      REAL(dp)                                           :: factor

      factor = 1.0_dp
      IF (r == c) factor = 0.5_dp
      DO j = 1, SIZE(block, 2)
         DO i = 1, SIZE(block, 1)
            ao_operator(r + i - 1, c + j - 1) = ao_operator(r + i - 1, c + j - 1) + factor*phase*block(i, j)
            ao_operator(c + j - 1, r + i - 1) = ao_operator(c + j - 1, r + i - 1) + factor*CONJG(phase)*block(i, j)
         END DO
      END DO
   END SUBROUTINE add_operator

END MODULE qs_upf_inverse_background
