!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Group-local SOC assembly and reconstruction of complete property symmetry orbits.
!> \author Thomas D. Kuehne, tkuehne@cp2k.org
! **************************************************************************************************
MODULE qs_property_soc
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_gemm
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
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
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE kpoint_types,                    ONLY: kpoint_env_type,&
                                              kpoint_type
   USE mathconstants,                   ONLY: gaussi,&
                                              z_one,&
                                              z_zero
   USE mathlib,                         ONLY: diag_complex
   USE message_passing,                 ONLY: mp_para_env_type
   USE qs_kpoint_operators,             ONLY: kpoint_operator_context_create,&
                                              kpoint_operator_context_release,&
                                              kpoint_operator_context_type,&
                                              kpoint_operator_evaluate_at,&
                                              kpoint_operator_prepare,&
                                              kpoint_operator_work_size,&
                                              kpoint_slot_ks
   USE qs_mo_types,                     ONLY: get_mo_set
   USE qs_property_symmetry,            ONLY: property_symmetry_frame,&
                                              property_symmetry_spin_frame
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_property_soc'
   PUBLIC :: property_soc_groups
CONTAINS
! **************************************************************************************************
!> \brief Assemble SOC on each orbit owner; diagonalize representatives and validate every target.
!> \param points Initialized representatives, distributed across k-point groups
!> \param xkp Full property graph coordinates
!> \param representative Source representative index of each graph point
!> \param operation Accepted scalar atom/AO operation of each graph point
!> \param source_point Graph index of each representative
!> \param band_map Selected spinor bands
!> \param matrices Real-space imaginary SOC components, retained unchanged throughout the call
!> \param global_ao Global AO matrix template for the bounded operator distribution
!> \param tolerance Normalization and scaled Hamiltonian residual tolerance
!> \param max_memory_mb Bound on additional replicated orbit work arrays per rank
!> \param values Full spinor spectrum at every graph point
!> \param coefficients Selected second-variational spinors at every graph point
!> \param errors Maximum normalization and Hamiltonian residuals over all groups
!> \param prepared Whether bounded operator snapshots were available; otherwise no result is returned
!> \param iw Output unit
! **************************************************************************************************
   SUBROUTINE property_soc_groups(points, xkp, representative, operation, source_point, band_map, &
                                  matrices, global_ao, tolerance, max_memory_mb, values, coefficients, &
                                  errors, prepared, iw)
      TYPE(kpoint_type), POINTER                         :: points
      REAL(dp), INTENT(IN)                               :: xkp(:, :)
      INTEGER, INTENT(IN)                                :: representative(:), operation(:), &
                                                            source_point(:), band_map(:)
      TYPE(dbcsr_p_type), POINTER                        :: matrices(:, :)
      TYPE(cp_fm_type), INTENT(IN)                       :: global_ao
      REAL(dp), INTENT(IN)                               :: tolerance
      INTEGER, INTENT(IN)                                :: max_memory_mb
      REAL(dp), INTENT(OUT)                              :: values(:, :)
      COMPLEX(dp), INTENT(OUT)                           :: coefficients(:, :, :)
      REAL(dp), INTENT(OUT)                              :: errors(2)
      LOGICAL, INTENT(OUT)                               :: prepared
      INTEGER, INTENT(IN)                                :: iw

      COMPLEX(dp), ALLOCATABLE                           :: frame(:, :), h(:, :), phase(:, :), &
                                                            source(:, :), source_u(:, :), u(:, :), &
                                                            work(:, :), xyz(:, :, :)
      INTEGER                                            :: axis, b, ik, ilocal, ir, member, nao, &
                                                            nentry, nmo, status
      INTEGER(int_8)                                     :: bytes
      INTEGER, ALLOCATABLE                               :: owners(:)
      LOGICAL                                            :: antiunitary
      REAL(dp)                                           :: rotation(3, 3)
      REAL(dp), ALLOCATABLE                              :: spectrum(:)
      REAL(dp), POINTER                                  :: scalar_values(:)
      TYPE(cp_cfm_type)                                  :: ao, c, projected, vc
      TYPE(cp_fm_struct_type), POINTER                   :: ao_struct, mo_struct, work_struct
      TYPE(cp_fm_type)                                   :: local_ao
      TYPE(kpoint_env_type), POINTER                     :: kp
      TYPE(kpoint_operator_context_type)                 :: context
      TYPE(mp_para_env_type), POINTER                    :: para

      para => points%blacs_env_all%para_env
      nmo = SIZE(values, 1)/2
      CALL get_mo_set(points%kp_env(1)%kpoint_env%mos(1), nao=nao)
      CPASSERT(SIZE(values, 1) == 2*nmo .AND. SIZE(values, 2) == SIZE(xkp, 2))
      CPASSERT(ALL(SHAPE(coefficients) == [2*nmo, SIZE(band_map), SIZE(xkp, 2)]))
      NULLIFY (ao_struct, mo_struct, work_struct)
      CALL cp_fm_struct_create(ao_struct, nrow_global=nao, ncol_global=nao, &
                               para_env=points%para_env_kp, context=points%blacs_env)
      CALL cp_fm_create(local_ao, ao_struct)
      CALL kpoint_operator_context_create(context, points, global_ao, 3, imaginary_operator=.TRUE.)
      CALL kpoint_operator_prepare(context, matrices, kpoint_slot_ks, prepared, TARGET=local_ao)
      errors = 0.0_dp
      IF (.NOT. prepared) THEN
         CALL kpoint_operator_context_release(context)
         CALL cp_fm_release(local_ao)
         CALL cp_fm_struct_release(ao_struct)
         RETURN
      END IF
      ! Include temporary atom/AO and spin-action arrays, but not the caller's full-graph results.
      bytes = 16_int_8*(INT(nao, int_8)*nao + 3_int_8*nao*nmo + 24_int_8*nmo*nmo)
      IF (bytes > INT(max_memory_mb, int_8)*1024*1024) THEN
         CPABORT('Property SOC orbit arrays exceed SYMMETRY_MAX_MEMORY_MB per rank.')
      END IF
      CALL cp_fm_struct_create(work_struct, nrow_global=nao, ncol_global=nmo, &
                               para_env=points%para_env_kp, context=points%blacs_env)
      CALL cp_fm_struct_create(mo_struct, nrow_global=nmo, ncol_global=nmo, &
                               para_env=points%para_env_kp, context=points%blacs_env)
      CALL cp_cfm_create(ao, ao_struct)
      CALL cp_cfm_create(c, work_struct)
      CALL cp_cfm_create(vc, work_struct)
      CALL cp_cfm_create(projected, mo_struct)
      CALL kpoint_operator_work_size(context, nentry)
      ALLOCATE (phase(1, nentry), source(nao, nmo), frame(nao, nmo))
      ALLOCATE (xyz(nmo, nmo, 3), h(2*nmo, 2*nmo), u(2*nmo, 2*nmo), &
                source_u(2*nmo, 2*nmo), work(2*nmo, 2*nmo), spectrum(2*nmo), owners(SIZE(xkp, 2)))
      values = 0.0_dp
      coefficients = z_zero
      owners = 0
      DO ir = points%kp_range(1), points%kp_range(2)
         ilocal = ir - points%kp_range(1) + 1
         kp => points%kp_env(ilocal)%kpoint_env
         CALL cp_cfm_get_submatrix(kp%mos(1)%cmo_coeff, source)
         CALL get_mo_set(kp%mos(1), eigenvalues=scalar_values)
         ! Visit the representative first, retaining only one complete SOC source frame per group.
         DO member = 0, SIZE(xkp, 2)
            ik = source_point(ir)
            IF (member > 0) THEN
               ik = member
               IF (representative(ik) /= ir .OR. ik == source_point(ir)) CYCLE
            END IF
            CALL property_symmetry_frame(points, ir, operation(ik), xkp(:, ik), source, .FALSE., &
                                         frame, rotation, antiunitary, status)
            IF (status /= 0) CPABORT('Could not reconstruct group-local scalar SOC basis.')
            CALL cp_cfm_set_submatrix(c, frame)
            DO axis = 1, 3
               CALL kpoint_operator_evaluate_at(context, axis, kpoint_slot_ks, xkp(:, ik), phase, TARGET=ao)
               CALL cp_cfm_gemm('N', 'N', nao, nmo, nao, z_one, ao, c, z_zero, vc)
               CALL cp_cfm_gemm('C', 'N', nmo, nmo, nao, z_one, c, vc, z_zero, projected)
               CALL cp_cfm_get_submatrix(projected, xyz(:, :, axis))
            END DO
            h(:nmo, :nmo) = xyz(:, :, 3)
            h(nmo + 1:, nmo + 1:) = -xyz(:, :, 3)
            h(:nmo, nmo + 1:) = xyz(:, :, 1) + gaussi*xyz(:, :, 2)
            h(nmo + 1:, :nmo) = xyz(:, :, 1) - gaussi*xyz(:, :, 2)
            DO b = 1, nmo
               h(b, b) = h(b, b) + scalar_values(b)
               h(nmo + b, nmo + b) = h(nmo + b, nmo + b) + scalar_values(b)
            END DO
            IF (MAXVAL(ABS(h - CONJG(TRANSPOSE(h)))) > 1.e-8_dp) THEN
               CPABORT('Property SOC: non-Hermitian group-local Hamiltonian.')
            END IF
            IF (member == 0) THEN
               CALL diag_complex(h, source_u, spectrum)
               u = source_u
            ELSE
               CALL property_symmetry_spin_frame(points, ir, operation(ik), source_u, u, status)
               IF (status /= 0) CPABORT('Could not reconstruct group-local property SOC frame.')
            END IF
            work = MATMUL(CONJG(TRANSPOSE(u)), u)
            DO b = 1, 2*nmo
               work(b, b) = work(b, b) - 1.0_dp
            END DO
            IF (.NOT. ALL(ieee_is_finite(REAL(work, dp))) .OR. &
                .NOT. ALL(ieee_is_finite(AIMAG(work)))) CPABORT('Nonfinite property SOC normalization.')
            errors(1) = MAX(errors(1), MAXVAL(ABS(work)))
            work = MATMUL(h, u)
            DO b = 1, 2*nmo
               work(:, b) = work(:, b) - u(:, b)*spectrum(b)
            END DO
            errors(2) = MAX(errors(2), MAXVAL(ABS(work))/MAX(1.0_dp, MAXVAL(ABS(h))))
            IF (.NOT. ALL(ieee_is_finite(REAL(work, dp))) .OR. &
                .NOT. ALL(ieee_is_finite(AIMAG(work))) .OR. MAXVAL(errors) > tolerance) THEN
               CPABORT('Property symmetry does not preserve the SOC Hamiltonian or scalar band window.')
            END IF
            IF (points%para_env_kp%is_source()) THEN
               values(:, ik) = spectrum
               coefficients(:, :, ik) = u(:, band_map)
               owners(ik) = 1
            END IF
         END DO
      END DO
      CALL para%sum(owners)
      CPASSERT(ALL(owners == 1))
      CALL para%sum(values)
      CALL para%sum(coefficients)
      CALL para%max(errors)
      IF (iw > 0) WRITE (iw, '(T2,A,I0,A,I0)') &
         'TOPOLOGY| Property SOC MPI groups: ', points%nkp_groups, &
         ' ranks per group: ', points%para_env_kp%num_pe
      CALL kpoint_operator_context_release(context)
      CALL cp_cfm_release(ao)
      CALL cp_cfm_release(c)
      CALL cp_cfm_release(vc)
      CALL cp_cfm_release(projected)
      CALL cp_fm_release(local_ao)
      CALL cp_fm_struct_release(ao_struct)
      CALL cp_fm_struct_release(work_struct)
      CALL cp_fm_struct_release(mo_struct)
   END SUBROUTINE property_soc_groups
END MODULE qs_property_soc
