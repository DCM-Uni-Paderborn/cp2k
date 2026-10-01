!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Geometry-dependent scalar USPP state, separate from the Gaussian overlap.
! **************************************************************************************************
MODULE qs_uspp_types
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type,&
                                              dbcsr_release_p,&
                                              dbcsr_type
   USE kinds,                           ONLY: dp
   USE kpoint_transitional,             ONLY: kpoint_transitional_release,&
                                              kpoint_transitional_type,&
                                              set_1d_pointer
   USE particle_types,                  ONLY: particle_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_uspp_projectors,              ONLY: qs_uspp_projector_matrix
   USE sap_kind_types,                  ONLY: sap_int_type
   USE upf_projector_integrals,         ONLY: upf_projector_init
   USE uspp_augmentation,               ONLY: uspp_augmentation_init,&
                                              uspp_augmentation_type
   USE uspp_projector_operators,        ONLY: uspp_projector_ao_pattern,&
                                              uspp_projector_metric
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_uspp_types'

   TYPE qs_uspp_type
      INTEGER, ALLOCATABLE :: center_atoms(:), center_kinds(:), projector_sizes(:)
      REAL(KIND=dp), ALLOCATABLE :: position(:, :)
      TYPE(uspp_augmentation_type), ALLOCATABLE :: augmentation(:)
      TYPE(dbcsr_type), POINTER :: overlaps => NULL(), ao_pattern => NULL()
      TYPE(kpoint_transitional_type) :: metric
   END TYPE qs_uspp_type

   PUBLIC :: qs_uspp_type, qs_uspp_init, qs_uspp_update, qs_uspp_release

CONTAINS

! **************************************************************************************************
!> \brief Select actual USPP centers and prepare per-kind augmentation in native beta order.
!>        No state is allocated for a system without USPP centers. Ghost atoms never augment.
!> \param uspp Owned state, replaced if already present.
!> \param qs_kind_set Kinds with original scalar USPP data and native projector caches.
!> \param particle_set Global particles, indexed in the same order as AO matrix blocks.
! **************************************************************************************************
   SUBROUTINE qs_uspp_init(uspp, qs_kind_set, particle_set)
      TYPE(qs_uspp_type), POINTER                        :: uspp
      TYPE(qs_kind_type), INTENT(INOUT)                  :: qs_kind_set(:)
      TYPE(particle_type), INTENT(IN)                    :: particle_set(:)

      INTEGER                                            :: atom, center, degree, ikind, ncenter
      LOGICAL                                            :: active(SIZE(qs_kind_set))

      CALL qs_uspp_release(uspp)
      active = .FALSE.
      DO ikind = 1, SIZE(qs_kind_set)
         IF (qs_kind_set(ikind)%ghost) CYCLE
         IF (.NOT. ASSOCIATED(qs_kind_set(ikind)%sgp_potential)) CYCLE
         ASSOCIATE (potential => qs_kind_set(ikind)%sgp_potential)
            IF (.NOT. ALLOCATED(potential%upf_data)) CYCLE
            active(ikind) = potential%upf_data%is_ultrasoft
         END ASSOCIATE
      END DO
      IF (.NOT. ANY(active)) RETURN
      ncenter = 0
      DO atom = 1, SIZE(particle_set)
         CPASSERT(ASSOCIATED(particle_set(atom)%atomic_kind))
         ikind = particle_set(atom)%atomic_kind%kind_number
         CPASSERT(ikind >= 1 .AND. ikind <= SIZE(qs_kind_set))
         IF (active(ikind)) ncenter = ncenter + 1
      END DO
      IF (ncenter == 0) RETURN
      ALLOCATE (uspp)
      ALLOCATE (uspp%center_atoms(ncenter), uspp%center_kinds(ncenter), uspp%projector_sizes(ncenter))
      ALLOCATE (uspp%position(3, ncenter), uspp%augmentation(SIZE(qs_kind_set)))
      center = 0
      DO atom = 1, SIZE(particle_set)
         ikind = particle_set(atom)%atomic_kind%kind_number
         IF (.NOT. active(ikind)) CYCLE
         ASSOCIATE (potential => qs_kind_set(ikind)%sgp_potential)
            IF (uspp%augmentation(ikind)%nproj == 0) THEN
               CPASSERT(potential%upf_data%number_of_proj > 0)
               degree = MAXVAL(potential%upf_data%lbeta)
               IF (potential%upf_projectors%degree < degree) THEN
                  CALL upf_projector_init(potential%upf_data, degree, potential%upf_projectors)
               END IF
               CALL uspp_augmentation_init(potential%upf_data, potential%upf_projectors, uspp%augmentation(ikind))
            END IF
         END ASSOCIATE
         center = center + 1
         uspp%center_atoms(center) = atom
         uspp%center_kinds(center) = ikind
         uspp%projector_sizes(center) = uspp%augmentation(ikind)%nproj
         uspp%position(:, center) = particle_set(atom)%r
      END DO
   END SUBROUTINE qs_uspp_init

! **************************************************************************************************
!> \brief Rebuild B and S from current neighbor integrals, allowing projector-mediated fill-in.
!>        The Gaussian overlap and its derivatives remain owned by the ordinary KS environment.
!> \param uspp Initialized USPP state.
!> \param gaussian_overlap Current Gaussian overlap and AO block distribution.
!> \param sap_int Current uncompressed AO/projector integrals including periodic images.
!> \param particle_set Current particle positions; the kind mapping must be unchanged.
! **************************************************************************************************
   SUBROUTINE qs_uspp_update(uspp, gaussian_overlap, sap_int, particle_set)
      TYPE(qs_uspp_type), INTENT(INOUT)                  :: uspp
      TYPE(dbcsr_type), INTENT(IN)                       :: gaussian_overlap
      TYPE(sap_int_type), INTENT(IN)                     :: sap_int(:)
      TYPE(particle_type), INTENT(IN)                    :: particle_set(:)

      INTEGER                                            :: atom, center
      TYPE(dbcsr_p_type), POINTER                        :: metric(:)

      CALL dbcsr_release_p(uspp%overlaps)
      CALL dbcsr_release_p(uspp%ao_pattern)
      CALL kpoint_transitional_release(uspp%metric)
      DO center = 1, SIZE(uspp%center_atoms)
         atom = uspp%center_atoms(center)
         CPASSERT(atom >= 1 .AND. atom <= SIZE(particle_set))
         CPASSERT(ASSOCIATED(particle_set(atom)%atomic_kind))
         CPASSERT(particle_set(atom)%atomic_kind%kind_number == uspp%center_kinds(center))
         uspp%position(:, center) = particle_set(atom)%r
      END DO
      ALLOCATE (uspp%overlaps)
      CALL qs_uspp_projector_matrix(gaussian_overlap, sap_int, uspp%center_atoms, uspp%projector_sizes, uspp%overlaps)
      ALLOCATE (uspp%ao_pattern)
      CALL uspp_projector_ao_pattern(uspp%overlaps, gaussian_overlap, uspp%ao_pattern)
      ALLOCATE (metric(1))
      ALLOCATE (metric(1)%matrix)
      CALL uspp_projector_metric(uspp%overlaps, uspp%augmentation, uspp%center_kinds, gaussian_overlap, metric(1)%matrix)
      CALL set_1d_pointer(uspp%metric, metric)
   END SUBROUTINE qs_uspp_update

! **************************************************************************************************
!> \brief Release only USPP-owned matrices and data; safe also for an absent state.
!> \param uspp Owned state, null on return.
! **************************************************************************************************
   SUBROUTINE qs_uspp_release(uspp)
      TYPE(qs_uspp_type), POINTER                        :: uspp

      IF (.NOT. ASSOCIATED(uspp)) RETURN
      CALL dbcsr_release_p(uspp%overlaps)
      CALL dbcsr_release_p(uspp%ao_pattern)
      CALL kpoint_transitional_release(uspp%metric)
      DEALLOCATE (uspp)
   END SUBROUTINE qs_uspp_release

END MODULE qs_uspp_types
