!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Native analytic Gaussian GTH PPL matrices and their isolated geometry adjoints.
!>        This is the full Gaussian base contribution. The inverse reconstruction and original
!>        UPF augmentation corrections remain separate and must be combined by the caller.
! **************************************************************************************************
MODULE qs_upf_inverse_gth_local
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE core_ppl,                        ONLY: build_core_ppl
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type,&
                                              dbcsr_set
   USE kinds,                           ONLY: dp
   USE kpoint_types,                    ONLY: kpoint_type
   USE particle_types,                  ONLY: particle_type
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_upf_inverse_background,       ONLY: qs_upf_inverse_background_operator,&
                                              qs_upf_inverse_background_set_density
   USE qs_upf_inverse_types,            ONLY: qs_upf_inverse_background_type
   USE virial_types,                    ONLY: virial_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: qs_upf_inverse_gth_local_matrix, qs_upf_inverse_gth_local_geometry

CONTAINS

! **************************************************************************************************
!> \brief Read the original native GTH PPL integral matrices in the physical Bloch AO ordering.
!> \param ks_env Native kinds, full orbital bases and unchanged neighbor lists.
!> \param map Owned image scratch matrices and the native-to-physical Bloch phase map.
!> \param operators New full-Gaussian PPL matrices, without k-point weights.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_gth_local_matrix(ks_env, map, operators)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(qs_upf_inverse_background_type), &
         INTENT(INOUT)                                   :: map
      COMPLEX(dp), INTENT(OUT)                           :: operators(:, :, :)

      CALL native_local(ks_env, map, operators=operators)
   END SUBROUTINE qs_upf_inverse_gth_local_matrix

! **************************************************************************************************
!> \brief Native PPL energy gradients and virial at explicitly supplied accepted AO occupations.
!>        The native force, virial, density and operator owners are unchanged. Only the owned map's
!>        image scratch arrays are overwritten. Full Gaussian bases are used for every density mode.
!> \param ks_env Original native integral environment.
!> \param map Owned native image and phase map, also carrying the original k-point weights.
!> \param matrix Hermitian occupations indexed by global k, without mesh weights.
!> \param gradients Replicated Cartesian energy gradients, without a force sign change.
!> \param virial Replicated native PPL virial.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_gth_local_geometry(ks_env, map, matrix, gradients, virial)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(qs_upf_inverse_background_type), &
         INTENT(INOUT)                                   :: map
      COMPLEX(dp), INTENT(IN)                            :: matrix(:, :, :)
      REAL(dp), INTENT(OUT)                              :: gradients(:, :), virial(3, 3)

      CALL native_local(ks_env, map, matrix=matrix, gradients=gradients, local_virial=virial)
   END SUBROUTINE qs_upf_inverse_gth_local_geometry

! **************************************************************************************************
!> \brief Select native GTH potential centers while retaining all orbital centers and contractions.
! **************************************************************************************************
   SUBROUTINE native_local(ks_env, map, operators, matrix, gradients, local_virial)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(qs_upf_inverse_background_type), &
         INTENT(INOUT)                                   :: map
      COMPLEX(dp), INTENT(OUT), OPTIONAL                 :: operators(:, :, :)
      COMPLEX(dp), INTENT(IN), OPTIONAL                  :: matrix(:, :, :)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: gradients(:, :), local_virial(3, 3)

      INTEGER                                            :: image, kind, nder, nimages
      INTEGER, POINTER                                   :: cell_to_index(:, :, :)
      LOGICAL                                            :: calculate_forces, use_virial
      LOGICAL, ALLOCATABLE                               :: selected(:)
      TYPE(atomic_kind_type), POINTER                    :: atomic_kinds(:)
      TYPE(dbcsr_p_type), POINTER                        :: density(:, :), potential(:, :)
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(neighbor_list_set_p_type), POINTER            :: centers(:), overlap(:)
      TYPE(particle_type), POINTER                       :: atoms(:)
      TYPE(qs_force_type), POINTER                       :: force(:) => NULL()
      TYPE(qs_kind_type), POINTER                        :: kinds(:)
      TYPE(virial_type), POINTER                         :: virial => NULL()

      CPASSERT(ASSOCIATED(map%potential))
      CALL get_ks_env(ks_env, atomic_kind_set=atomic_kinds, qs_kind_set=kinds, particle_set=atoms, &
                      sab_orb=overlap, sac_ppl=centers, kpoints=kpoints)
      calculate_forces = PRESENT(matrix)
      IF (PRESENT(operators)) operators = CMPLX(0.0_dp, 0.0_dp, dp)
      IF (calculate_forces) THEN
         CPASSERT(PRESENT(gradients))
         CPASSERT(PRESENT(local_virial))
         CPASSERT(SIZE(gradients, 1) == 3)
         CPASSERT(SIZE(gradients, 2) == SIZE(atoms))
         gradients = 0.0_dp
         local_virial = 0.0_dp
      END IF
      ALLOCATE (selected(SIZE(kinds)), SOURCE=.FALSE.)
      DO kind = 1, SIZE(kinds)
         IF (kinds(kind)%ghost) CYCLE
         IF (.NOT. ASSOCIATED(kinds(kind)%gth_potential)) CYCLE
         ASSOCIATE (pot => kinds(kind)%gth_potential)
            selected(kind) = pot%nexp_ppl > 0
            IF (pot%lpotextended) selected(kind) = selected(kind) .OR. pot%nexp_lpot > 0
         END ASSOCIATE
      END DO
      IF (.NOT. ANY(selected)) RETURN
      CPASSERT(ASSOCIATED(overlap))
      CPASSERT(ASSOCIATED(centers))
      IF (calculate_forces) CALL qs_upf_inverse_background_set_density(map, matrix)
      nimages = SIZE(map%potential)
      ALLOCATE (potential(1, nimages), density(1, nimages))
      DO image = 1, nimages
         CALL dbcsr_set(map%potential(image)%matrix, 0.0_dp)
         potential(1, image)%matrix => map%potential(image)%matrix
         density(1, image)%matrix => map%density(image)%matrix
      END DO
      NULLIFY (cell_to_index)
      IF (ASSOCIATED(kpoints)) cell_to_index => kpoints%cell_to_index
      nder = 0
      IF (calculate_forces) nder = 1
      use_virial = calculate_forces
      CALL build_core_ppl(potential, density, force, virial, calculate_forces, use_virial, nder, &
                          kinds, atomic_kinds, atoms, overlap, centers, nimages, cell_to_index, "ORB", &
                          kind_mask=selected, force_response=gradients, virial_response=local_virial)
      IF (PRESENT(operators)) CALL qs_upf_inverse_background_operator(map, map%potential, operators)
      IF (calculate_forces) THEN
         CALL map%comm%sum(gradients)
         CALL map%comm%sum(local_virial)
      END IF
      DEALLOCATE (potential, density)
   END SUBROUTINE native_local
END MODULE qs_upf_inverse_gth_local
