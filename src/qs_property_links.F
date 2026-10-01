!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Bounded group-local assembly of directed cross-k property links.
!> \author Thomas D. Kuehne, tkuehne@cp2k.org
! **************************************************************************************************
MODULE qs_property_links
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_gemm,&
                                              cp_cfm_scale_and_add
   USE cp_cfm_types,                    ONLY: cp_cfm_copy_general,&
                                              cp_cfm_create,&
                                              cp_cfm_get_submatrix,&
                                              cp_cfm_release,&
                                              cp_cfm_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE kpoint_methods,                  ONLY: kpoint_init_cell_index
   USE kpoint_types,                    ONLY: kpoint_create,&
                                              kpoint_release,&
                                              kpoint_type
   USE mathconstants,                   ONLY: gaussi,&
                                              z_one,&
                                              z_zero
   USE message_passing,                 ONLY: mp_para_env_type
   USE qs_kpoint_operators,             ONLY: kpoint_operator_context_create,&
                                              kpoint_operator_context_release,&
                                              kpoint_operator_context_type,&
                                              kpoint_operator_evaluate_at,&
                                              kpoint_operator_prepare,&
                                              kpoint_operator_work_size,&
                                              kpoint_slot_ks
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_property_links'

   TYPE berry_matrix_type
      TYPE(dbcsr_p_type), POINTER :: sinmat(:, :) => NULL(), cosmat(:, :) => NULL()
   END TYPE berry_matrix_type

   TYPE property_link_context_type
      PRIVATE
      LOGICAL :: prepared = .FALSE.
      INTEGER :: nao = 0, nmo = 0, nedge = 0, delivered = 0
      INTEGER, ALLOCATABLE :: cached_points(:, :)
      TYPE(kpoint_type), POINTER :: points => NULL()
      TYPE(kpoint_operator_context_type) :: cosine, sine
      TYPE(dbcsr_p_type), POINTER :: cosmat(:, :) => NULL(), sinmat(:, :) => NULL()
      TYPE(cp_cfm_type) :: c_left, c_right, OPERATOR, sine_operator, work, product
      COMPLEX(dp), ALLOCATABLE :: phase(:, :), scalar(:, :), links(:, :, :)
   END TYPE property_link_context_type

   PUBLIC :: berry_matrix_type, property_link_context_type
   PUBLIC :: property_links_prepare, property_links_get, property_links_release
CONTAINS
! **************************************************************************************************
!> \brief Capture ordered operators and allocate one link workspace per MPI group.
!> \param context Owned link context
!> \param groups Existing representative MPI distribution
!> \param geometry Ordered Berry neighbor list and image convention
!> \param matrices Borrowed cosine and sine matrices, unchanged until release
!> \param global_ao Global AO layout template, retained by the caller
!> \param nmo Scalar orbital count
!> \param nband Selected scalar or spinor band count
!> \param nedge Number of directed links to deliver in original order
!> \param max_memory_mb Bound for additional link work arrays
!> \param prepared Whether both bounded snapshots are usable on every rank
!> \param iw Output unit
! **************************************************************************************************
   SUBROUTINE property_links_prepare(context, groups, geometry, matrices, global_ao, nmo, nband, nedge, &
                                     max_memory_mb, prepared, iw)
      TYPE(property_link_context_type), INTENT(OUT)      :: context
      TYPE(kpoint_type), POINTER                         :: groups, geometry
      TYPE(berry_matrix_type), INTENT(IN)                :: matrices(:)
      TYPE(cp_fm_type), INTENT(IN)                       :: global_ao
      INTEGER, INTENT(IN)                                :: nmo, nband, nedge, max_memory_mb
      LOGICAL, INTENT(OUT)                               :: prepared
      INTEGER, INTENT(IN)                                :: iw

      INTEGER                                            :: i, image, nao, nentry, nimage, nother
      INTEGER(int_8)                                     :: bytes
      LOGICAL                                            :: ready_c, ready_s
      TYPE(cp_fm_struct_type), POINTER                   :: ao_struct, mo_struct, work_struct
      TYPE(cp_fm_type)                                   :: local_ao

      prepared = .FALSE.
      CPASSERT(nedge > 0 .AND. SIZE(matrices) > 0)
      nao = global_ao%matrix_struct%nrow_global
      bytes = 16_int_8*(2_int_8*nao*nao + 5_int_8*nao*nmo + 16_int_8*nmo*nmo + &
                        INT(nband, int_8)*nband*groups%nkp_groups)
      IF (bytes > INT(max_memory_mb, int_8)*1024*1024) THEN
         IF (iw > 0) WRITE (iw, '(T2,A,F12.3,A,I0)') &
            'TOPOLOGY| Directed-link work bound [MiB]: ', REAL(bytes, dp)/1024.0_dp**2, &
            ' exceeds SYMMETRY_MAX_MEMORY_MB: ', max_memory_mb
         RETURN
      END IF
      ! Borrow only the parallel environments; the ordered image map is independently owned.
      CALL kpoint_create(context%points)
      context%points%para_env => groups%para_env
      CALL context%points%para_env%retain()
      context%points%para_env_kp => groups%para_env_kp
      CALL context%points%para_env_kp%retain()
      context%points%para_env_inter_kp => groups%para_env_inter_kp
      CALL context%points%para_env_inter_kp%retain()
      context%points%blacs_env_all => groups%blacs_env_all
      CALL context%points%blacs_env_all%retain()
      context%points%blacs_env => groups%blacs_env
      CALL context%points%blacs_env%retain()
      context%points%nkp_groups = groups%nkp_groups
      context%points%use_real_wfn = .FALSE.
      CALL kpoint_init_cell_index(context%points, geometry%sab_nl, groups%para_env, nimage)
      CPASSERT(ALL(context%points%index_to_cell == geometry%index_to_cell))
      ALLOCATE (context%cosmat(SIZE(matrices), nimage), context%sinmat(SIZE(matrices), nimage))
      DO i = 1, SIZE(matrices)
         CPASSERT(SIZE(matrices(i)%cosmat, 2) == nimage)
         CPASSERT(SIZE(matrices(i)%sinmat, 2) == nimage)
         DO image = 1, nimage
            context%cosmat(i, image)%matrix => matrices(i)%cosmat(1, image)%matrix
            context%sinmat(i, image)%matrix => matrices(i)%sinmat(1, image)%matrix
         END DO
      END DO
      NULLIFY (ao_struct, mo_struct, work_struct)
      CALL cp_fm_struct_create(ao_struct, nrow_global=nao, ncol_global=nao, &
                               para_env=groups%para_env_kp, context=groups%blacs_env)
      CALL cp_fm_create(local_ao, ao_struct)
      CALL kpoint_operator_context_create(context%cosine, context%points, global_ao, SIZE(matrices))
      CALL kpoint_operator_context_create(context%sine, context%points, global_ao, SIZE(matrices))
      CALL kpoint_operator_prepare(context%cosine, context%cosmat, kpoint_slot_ks, ready_c, TARGET=local_ao)
      CALL kpoint_operator_prepare(context%sine, context%sinmat, kpoint_slot_ks, ready_s, TARGET=local_ao)
      prepared = ready_c .AND. ready_s
      IF (prepared) THEN
         context%nao = nao
         context%nmo = nmo
         context%nedge = nedge
         context%prepared = .TRUE.
         CALL cp_fm_struct_create(work_struct, nrow_global=nao, ncol_global=nmo, &
                                  para_env=groups%para_env_kp, context=groups%blacs_env)
         CALL cp_fm_struct_create(mo_struct, nrow_global=nmo, ncol_global=nmo, &
                                  para_env=groups%para_env_kp, context=groups%blacs_env)
         CALL cp_cfm_create(context%c_left, work_struct)
         CALL cp_cfm_create(context%c_right, work_struct)
         CALL cp_cfm_create(context%work, work_struct)
         CALL cp_cfm_create(context%product, mo_struct)
         CALL cp_cfm_create(context%operator, ao_struct)
         CALL cp_cfm_create(context%sine_operator, ao_struct)
         CALL kpoint_operator_work_size(context%cosine, nentry)
         CALL kpoint_operator_work_size(context%sine, nother)
         ALLOCATE (context%phase(1, MAX(nentry, nother)), context%scalar(nmo, nmo), &
                   context%links(nband, nband, groups%nkp_groups))
         ALLOCATE (context%cached_points(2, groups%nkp_groups), SOURCE=0)
         IF (iw > 0) WRITE (iw, '(T2,A,I0,A,I0)') &
            'TOPOLOGY| Property directed-link MPI groups: ', groups%nkp_groups, &
            ' ranks per group: ', groups%para_env_kp%num_pe
      END IF
      CALL cp_fm_release(local_ao)
      CALL cp_fm_struct_release(ao_struct)
      CALL cp_fm_struct_release(mo_struct)
      CALL cp_fm_struct_release(work_struct)
      IF (.NOT. prepared) CALL property_links_release(context)
   END SUBROUTINE property_links_prepare

! **************************************************************************************************
!> \brief Return the next directed link, computing a bounded batch on independent MPI groups.
!> \param context Prepared context; calls must follow the original edge order
!> \param full Full-graph scalar MOs distributed over the complete MPI environment
!> \param spin Scalar spin channel
!> \param neighbors Ordered destination indices (source point, neighbor)
!> \param operators Signed cosine/sine operator index for each directed edge
!> \param band_map Selected scalar bands when spinors are absent
!> \param spinors Optional allocated post-SCF spinors in each full-point scalar frame
!> \param edge Linear source-major edge index, without reordering Wilson products
!> \param link Selected-band connection matrix
! **************************************************************************************************
   SUBROUTINE property_links_get(context, full, spin, neighbors, operators, band_map, spinors, edge, link)
      TYPE(property_link_context_type), INTENT(INOUT), &
         TARGET                                          :: context
      TYPE(kpoint_type), POINTER                         :: full
      INTEGER, INTENT(IN)                                :: spin, neighbors(:, :), operators(:, :), &
                                                            band_map(:)
      COMPLEX(dp), ALLOCATABLE, INTENT(IN)               :: spinors(:, :, :)
      INTEGER, INTENT(IN)                                :: edge
      COMPLEX(dp), INTENT(OUT)                           :: link(:, :)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'property_links_get'

      INTEGER                                            :: first, group, handle, i, ibs, ik, ik2, &
                                                            item, local_edge, mygroup, n, ngroup, &
                                                            nntot, slot
      INTEGER, ALLOCATABLE                               :: owners(:)
      TYPE(cp_cfm_type), POINTER                         :: destination
      TYPE(cp_cfm_type), TARGET                          :: empty
      TYPE(mp_para_env_type), POINTER                    :: para

      CALL timeset(routineN, handle)
      CPASSERT(context%prepared)
      CPASSERT(edge == context%delivered + 1 .AND. edge <= context%nedge)
      CPASSERT(ALL(full%kp_range == [1, full%nkp]))
      para => context%points%para_env
      ngroup = context%points%nkp_groups
      mygroup = context%points%para_env_inter_kp%mepos + 1
      slot = MODULO(edge - 1, ngroup) + 1
      nntot = SIZE(neighbors, 2)
      n = context%nmo
      IF (slot == 1) THEN
         first = edge
         context%links(:, :, :) = z_zero
         ALLOCATE (owners(ngroup), SOURCE=0)
         ! All groups join each MO transfer, including groups idle in the last batch.
         DO group = 1, MIN(ngroup, context%nedge - first + 1)
            item = first + group - 1
            ik = (item - 1)/nntot + 1
            i = MODULO(item - 1, nntot) + 1
            ik2 = neighbors(ik, i)
            IF (context%cached_points(1, group) /= ik) THEN
               destination => empty
               IF (mygroup == group) destination => context%c_left
               CALL cp_cfm_copy_general(full%kp_env(ik)%kpoint_env%mos(spin)%cmo_coeff, destination, para)
            END IF
            IF (context%cached_points(2, group) /= ik2) THEN
               destination => empty
               IF (mygroup == group) destination => context%c_right
               CALL cp_cfm_copy_general(full%kp_env(ik2)%kpoint_env%mos(spin)%cmo_coeff, destination, para)
            END IF
            context%cached_points(:, group) = [ik, ik2]
         END DO
         local_edge = first + mygroup - 1
         IF (local_edge <= context%nedge) THEN
            ik = (local_edge - 1)/nntot + 1
            i = MODULO(local_edge - 1, nntot) + 1
            ik2 = neighbors(ik, i)
            ibs = ABS(operators(ik, i))
            CALL kpoint_operator_evaluate_at(context%cosine, ibs, kpoint_slot_ks, full%xkp(:, ik2), &
                                             context%phase, TARGET=context%operator)
            CALL kpoint_operator_evaluate_at(context%sine, ibs, kpoint_slot_ks, full%xkp(:, ik2), &
                                             context%phase, TARGET=context%sine_operator)
            ! The stored b already has the ordered Berry sign. No transpose or conjugate mirror.
            CALL cp_cfm_scale_and_add(z_one, context%operator, &
                                      gaussi*REAL(SIGN(1, operators(ik, i)), dp), context%sine_operator)
            CALL cp_cfm_gemm('N', 'N', context%nao, n, context%nao, z_one, &
                             context%operator, context%c_right, z_zero, context%work)
            CALL cp_cfm_gemm('C', 'N', n, n, context%nao, z_one, &
                             context%c_left, context%work, z_zero, context%product)
            CALL cp_cfm_get_submatrix(context%product, context%scalar)
            IF (context%points%para_env_kp%is_source()) THEN
               IF (ALLOCATED(spinors)) THEN
                  context%links(:, :, mygroup) = &
                     MATMUL(CONJG(TRANSPOSE(spinors(:n, :, ik))), MATMUL(context%scalar, spinors(:n, :, ik2))) + &
                     MATMUL(CONJG(TRANSPOSE(spinors(n + 1:, :, ik))), MATMUL(context%scalar, spinors(n + 1:, :, ik2)))
               ELSE
                  context%links(:, :, mygroup) = context%scalar(band_map, band_map)
               END IF
               owners(mygroup) = 1
            END IF
         END IF
         CALL para%sum(owners)
         CPASSERT(ALL(owners(:MIN(ngroup, context%nedge - first + 1)) == 1))
         CPASSERT(ALL(owners(MIN(ngroup, context%nedge - first + 1) + 1:) == 0))
         CALL para%sum(context%links)
         IF (.NOT. ALL(ieee_is_finite(REAL(context%links, dp))) .OR. &
             .NOT. ALL(ieee_is_finite(AIMAG(context%links)))) THEN
            CPABORT('Nonfinite group-local property link.')
         END IF
      END IF
      link(:, :) = context%links(:, :, slot)
      context%delivered = edge
      CALL timestop(handle)
   END SUBROUTINE property_links_get

! **************************************************************************************************
!> \brief Release snapshots before their borrowed real-space matrices and MPI groups disappear.
!> \param context Link context to release, including an unsuccessful preparation
! **************************************************************************************************
   SUBROUTINE property_links_release(context)
      TYPE(property_link_context_type), INTENT(INOUT)    :: context

      IF (context%prepared) THEN
         CPASSERT(context%delivered == context%nedge)
         CALL cp_cfm_release(context%c_left)
         CALL cp_cfm_release(context%c_right)
         CALL cp_cfm_release(context%operator)
         CALL cp_cfm_release(context%sine_operator)
         CALL cp_cfm_release(context%work)
         CALL cp_cfm_release(context%product)
         DEALLOCATE (context%phase, context%scalar, context%links, context%cached_points)
      END IF
      CALL kpoint_operator_context_release(context%cosine)
      CALL kpoint_operator_context_release(context%sine)
      IF (ASSOCIATED(context%cosmat)) DEALLOCATE (context%cosmat)
      IF (ASSOCIATED(context%sinmat)) DEALLOCATE (context%sinmat)
      IF (ASSOCIATED(context%points)) CALL kpoint_release(context%points)
      context%prepared = .FALSE.
   END SUBROUTINE property_links_release
END MODULE qs_property_links
