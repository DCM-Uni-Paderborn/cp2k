!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
!> \brief Native Gaussian soft-grid source and its complex Bloch AO transpose for inverse GAPW.
MODULE qs_upf_inverse_background
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_copy, dbcsr_deallocate_matrix, dbcsr_get_block_p, dbcsr_get_info, &
        dbcsr_get_matrix_type, dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, &
        dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_p_type, dbcsr_scale, &
        dbcsr_set, dbcsr_type_symmetric
   USE kinds,                           ONLY: dp
   USE kpoint_methods,                  ONLY: kp_transform_plan_create,&
                                              kp_transform_plan_type
   USE kpoint_types,                    ONLY: kpoint_type
   USE mathconstants,                   ONLY: twopi
   USE message_passing,                 ONLY: mp_comm_type
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE qs_collocate_density,            ONLY: calculate_rho_elec
   USE qs_environment_types,            ONLY: qs_environment_type
   USE qs_integrate_potential_product,  ONLY: integrate_v_rspace
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   TYPE qs_upf_inverse_background_type
      PRIVATE
      INTEGER :: nao = 0
      LOGICAL :: bloch = .FALSE.
      TYPE(mp_comm_type) :: comm
      TYPE(kp_transform_plan_type) :: plan
      TYPE(dbcsr_p_type), POINTER :: density(:) => NULL(), potential(:) => NULL()
      REAL(dp), ALLOCATABLE :: reduced_k(:, :), weights(:)
   END TYPE qs_upf_inverse_background_type

   PUBLIC :: qs_upf_inverse_background_type, qs_upf_inverse_background_prepare
   PUBLIC :: qs_upf_inverse_background_release, qs_upf_inverse_background_collocate
   PUBLIC :: qs_upf_inverse_background_integrate

CONTAINS

! **************************************************************************************************
!> \brief Own the native image sparsity and unreduced Bloch phase map at the current geometry.
!>        Rebuild with native overlap lists and soft task lists after atomic or cell changes.
!> \param ks_env Live native KS environment with overlap images and Gaussian soft task lists.
!> \param map New owned background map. The native matrices and density state are unchanged.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_background_prepare(ks_env, map)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(qs_upf_inverse_background_type), INTENT(OUT)  :: map

      CHARACTER(LEN=1)                                   :: symmetry
      INTEGER                                            :: image
      INTEGER, POINTER                                   :: sizes(:)
      REAL(dp)                                           :: tolerance
      TYPE(dbcsr_p_type), POINTER                        :: overlap(:, :)
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(neighbor_list_set_p_type), POINTER            :: sab(:)

      CALL get_ks_env(ks_env, matrix_s_kp=overlap, kpoints=kpoints, do_kpoints=map%bloch, sab_orb=sab)
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
         IF (kpoints%symmetry .OR. .NOT. kpoints%full_grid) THEN
            CPABORT("Inverse Gaussian background requires an unreduced mesh")
         END IF
         CALL kp_transform_plan_create(map%plan, sab, kpoints%cell_to_index, SIZE(overlap, 2), &
                                       block_template=overlap(1, 1)%matrix)
         CPASSERT(map%plan%symmetric)
         ALLOCATE (map%reduced_k(3, kpoints%nkp), map%weights(kpoints%nkp))
         map%reduced_k(:, :) = kpoints%xkp
         map%weights(:) = kpoints%wkp
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
!> \brief Release owned image matrices without changing the native density or operator state.
!> \param map ...
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_background_release(map)
      TYPE(qs_upf_inverse_background_type), &
         INTENT(INOUT)                                   :: map

      INTEGER                                            :: image

      IF (ASSOCIATED(map%density)) THEN
         DO image = 1, SIZE(map%density)
            CALL dbcsr_deallocate_matrix(map%density(image)%matrix)
            CALL dbcsr_deallocate_matrix(map%potential(image)%matrix)
         END DO
         DEALLOCATE (map%density, map%potential)
      END IF
      map = qs_upf_inverse_background_type()
   END SUBROUTINE qs_upf_inverse_background_release

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
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_background_collocate(map, ks_env, matrix, rho, rho_gspace, total_rho)
      TYPE(qs_upf_inverse_background_type), &
         INTENT(INOUT)                                   :: map
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      COMPLEX(dp), INTENT(IN)                            :: matrix(:, :, :)
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: rho
      TYPE(pw_c1d_gs_type), INTENT(INOUT)                :: rho_gspace
      REAL(dp), INTENT(OUT), OPTIONAL                    :: total_rho

      CALL background_images(map, matrix)
      CALL calculate_rho_elec(matrix_p_kp=map%density, rho=rho, rho_gspace=rho_gspace, &
                              total_rho=total_rho, ks_env=ks_env, soft_valid=.TRUE.)
   END SUBROUTINE qs_upf_inverse_background_collocate

! **************************************************************************************************
!> \brief Exact complex AO transpose of native Gaussian soft collocation, without mesh weights.
!>        A unit potential gives the background contribution to the complete GAPW charge metric.
!> \param map ...
!> \param qs_env Native owner of the same soft Gaussian task lists and PW grids.
!> \param potential Native grid potential, without a grid measure or a mesh weight.
!> \param ao_operators Overwritten Hermitian AO matrices indexed by global k.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_background_integrate(map, qs_env, potential, ao_operators)
      TYPE(qs_upf_inverse_background_type), &
         INTENT(INOUT)                                   :: map
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: potential
      COMPLEX(dp), INTENT(OUT)                           :: ao_operators(:, :, :)

      COMPLEX(dp)                                        :: phase
      INTEGER                                            :: c, ENTRY, image, k, r
      LOGICAL                                            :: found, transposed
      REAL(dp), POINTER                                  :: block(:, :)
      TYPE(dbcsr_iterator_type)                          :: iterator

      CPASSERT(ASSOCIATED(map%potential))
      CPASSERT(ALL(SHAPE(ao_operators) == [map%nao, map%nao, SIZE(map%weights)]))
      DO image = 1, SIZE(map%potential)
         CALL dbcsr_set(map%potential(image)%matrix, 0.0_dp)
      END DO
      CALL integrate_v_rspace(potential, hmat_kp=map%potential, qs_env=qs_env, &
                              calculate_forces=.FALSE., gapw=.TRUE.)
      ! Native integration leaves the auxiliary PW volume element to its caller.
      DO image = 1, SIZE(map%potential)
         CALL dbcsr_scale(map%potential(image)%matrix, potential%pw_grid%dvol)
      END DO
      ao_operators = CMPLX(0.0_dp, 0.0_dp, dp)
      IF (map%bloch) THEN
         DO ENTRY = 1, map%plan%nentry
            image = map%plan%image(ENTRY)
            CALL dbcsr_get_block_p(map%potential(image)%matrix, map%plan%row(ENTRY), map%plan%col(ENTRY), &
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
         CALL dbcsr_iterator_start(iterator, map%potential(1)%matrix)
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
   END SUBROUTINE qs_upf_inverse_background_integrate

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
