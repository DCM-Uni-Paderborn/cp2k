!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Original scalar USPP charge augmentation in the inverse Quickstep projector frame.
!>        Charge kernels are periodic, without a Bloch phase; W carries the k dependence.
!>        Kinetic augmentation is independent and is not inferred from the charge kernels.
! **************************************************************************************************
MODULE qs_upf_inverse_augmentation
   USE kinds,                           ONLY: dp
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type
   USE qs_upf_inverse,                  ONLY: qs_upf_inverse_type
   USE uspp_augmentation,               ONLY: uspp_augmentation_collocate_atom,&
                                              uspp_augmentation_init,&
                                              uspp_augmentation_integrate_atom,&
                                              uspp_augmentation_prepare_realspace,&
                                              uspp_augmentation_realspace,&
                                              uspp_augmentation_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   TYPE qs_upf_inverse_augmentation_type
      TYPE(uspp_augmentation_type), ALLOCATABLE :: kinds(:)
   END TYPE qs_upf_inverse_augmentation_type

   PUBLIC :: qs_upf_inverse_augmentation_type, qs_upf_inverse_augmentation_prepare
   PUBLIC :: qs_upf_inverse_augmentation_collocate, qs_upf_inverse_augmentation_integrate
   PUBLIC :: qs_upf_inverse_augmentation_response, qs_upf_inverse_augmentation_release

CONTAINS

! **************************************************************************************************
!> \brief Own original Q caches in the native magnetic-projector order, without changing QS state.
!> \param ks_env Live native KS environment used for frame.
!> \param frame Prepared inverse frame with selected nonghost USPP centers.
!> \param augmentation New owned original radial/angular caches.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_augmentation_prepare(ks_env, frame, augmentation)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      TYPE(qs_upf_inverse_augmentation_type), &
         INTENT(OUT)                                     :: augmentation

      INTEGER                                            :: center, first, ikind, last
      REAL(dp)                                           :: metric_difference
      TYPE(qs_kind_type), POINTER                        :: kinds(:)

      CPASSERT(frame%inverse%ready)
      CALL get_ks_env(ks_env, qs_kind_set=kinds)
      CPASSERT(ASSOCIATED(kinds))
      ALLOCATE (augmentation%kinds(SIZE(kinds)))
      DO center = 1, SIZE(frame%center_atoms)
         ikind = frame%center_kinds(center)
         IF (augmentation%kinds(ikind)%nproj /= 0) CYCLE
         CALL uspp_augmentation_init(kinds(ikind)%sgp_potential%upf_data, &
                                     kinds(ikind)%sgp_potential%upf_projectors, augmentation%kinds(ikind))
         CALL uspp_augmentation_prepare_realspace(kinds(ikind)%sgp_potential%upf_data, augmentation%kinds(ikind))
         first = frame%projector_first(center)
         last = frame%projector_last(center)
         CPASSERT(augmentation%kinds(ikind)%nproj == last - first + 1)
         metric_difference = MAXVAL(ABS(augmentation%kinds(ikind)%q - REAL(frame%radial_q(first:last, first:last), dp)))
         CPASSERT(metric_difference == 0.0_dp)
      END DO
   END SUBROUTINE qs_upf_inverse_augmentation_prepare

! **************************************************************************************************
!> \brief Add original charge augmentation and its Cartesian gradient to matching local point rows.
!>        Only same-center projector blocks enter Q; inverse intercenter coupling remains in W.
!> \param frame ...
!> \param augmentation ...
!> \param projector_density Full Hermitian W P W^dagger, including occupations and k weights.
!> \param density Local density, incremented without weights.
!> \param gradient Local physical Cartesian density gradient, axis last, incremented.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_augmentation_collocate(frame, augmentation, projector_density, density, gradient)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      TYPE(qs_upf_inverse_augmentation_type), INTENT(IN) :: augmentation
      COMPLEX(dp), INTENT(IN)                            :: projector_density(:, :)
      REAL(dp), INTENT(INOUT)                            :: density(:), gradient(:, :)

      INTEGER                                            :: atom, center, first, ikind, image, last, &
                                                            n
      REAL(dp)                                           :: hermiticity, position(3), tolerance
      REAL(dp), ALLOCATABLE                              :: displacements(:, :), g(:, :)

      CPASSERT(frame%inverse%ready)
      n = frame%inverse%nproj
      CPASSERT(ALL(SHAPE(projector_density) == [n, n]))
      CPASSERT(SIZE(density) == SIZE(frame%weights))
      CPASSERT(ALL(SHAPE(gradient) == [SIZE(density), 3]))
      tolerance = 128.0_dp*EPSILON(1.0_dp)*MAX(1.0_dp, MAXVAL(ABS(projector_density)))
      hermiticity = MAXVAL(ABS(projector_density - CONJG(TRANSPOSE(projector_density))))
      CPASSERT(hermiticity <= tolerance)
      ALLOCATE (displacements(3, SIZE(density)), g(3, SIZE(density)))
      g = 0.0_dp
      DO center = 1, SIZE(frame%center_atoms)
         atom = frame%center_atoms(center)
         ikind = frame%center_kinds(center)
         first = frame%projector_first(center)
         last = frame%projector_last(center)
         DO image = 1, SIZE(frame%images, 2)
            position = frame%positions(:, atom) + MATMUL(frame%cell, REAL(frame%images(:, image), dp))
            displacements(:, :) = frame%points - SPREAD(position, 2, SIZE(density))
            CALL uspp_augmentation_collocate_atom(augmentation%kinds(ikind), displacements, &
                                                  REAL(projector_density(first:last, first:last), dp), density, g)
         END DO
      END DO
      gradient = gradient + TRANSPOSE(g)
   END SUBROUTINE qs_upf_inverse_augmentation_collocate

! **************************************************************************************************
!> \brief Integrate the adjoint of the same original charge/gradient image sum, with one reduction.
!> \param frame ...
!> \param augmentation ...
!> \param vdensity Unweighted potential conjugate to the augmented density.
!> \param vgradient Unweighted potential conjugate to its Cartesian gradient, axis last.
!> \param projector_operator New full projector operator; original Q contributes only same-center blocks.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_augmentation_integrate(frame, augmentation, vdensity, vgradient, projector_operator)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      TYPE(qs_upf_inverse_augmentation_type), INTENT(IN) :: augmentation
      REAL(dp), INTENT(IN)                               :: vdensity(:), vgradient(:, :)
      COMPLEX(dp), INTENT(OUT)                           :: projector_operator(:, :)

      INTEGER                                            :: atom, center, first, ikind, image, last
      REAL(dp)                                           :: position(3)
      REAL(dp), ALLOCATABLE                              :: displacements(:, :), matrix(:, :)

      CALL check_operator(frame, vdensity, vgradient, projector_operator)
      projector_operator = CMPLX(0.0_dp, 0.0_dp, dp)
      ALLOCATE (displacements(3, SIZE(frame%weights)))
      DO center = 1, SIZE(frame%center_atoms)
         atom = frame%center_atoms(center)
         ikind = frame%center_kinds(center)
         first = frame%projector_first(center)
         last = frame%projector_last(center)
         ALLOCATE (matrix(last - first + 1, last - first + 1))
         DO image = 1, SIZE(frame%images, 2)
            position = frame%positions(:, atom) + MATMUL(frame%cell, REAL(frame%images(:, image), dp))
            displacements(:, :) = frame%points - SPREAD(position, 2, SIZE(frame%weights))
            CALL uspp_augmentation_integrate_atom(augmentation%kinds(ikind), displacements, frame%weights, &
                                                  vdensity, TRANSPOSE(vgradient), matrix)
            projector_operator(first:last, first:last) = projector_operator(first:last, first:last) + matrix
         END DO
         DEALLOCATE (matrix)
      END DO
      CALL frame%point_comm%sum(projector_operator)
   END SUBROUTINE qs_upf_inverse_augmentation_integrate

! **************************************************************************************************
!> \brief Differentiate the same Q image sum, including its spatial gradient and original measure.
!>        Image membership is fixed. Inverse-coefficient responses belong to the later W pullback.
!> \param frame ...
!> \param augmentation ...
!> \param vdensity ...
!> \param vgradient ...
!> \param point_response Local Cartesian point responses, axis first.
!> \param atom_response Responses of all native atom positions, axis first.
!> \param cell_response Physical-cell response for the same fixed images.
!> \param weights_response Response of the same original local quadrature measure.
!> \param vdensity_response Explicit response of the conjugate potential.
!> \param vgradient_response ...
!> \param response New full projector-operator response, reduced once over point rows.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_augmentation_response(frame, augmentation, vdensity, vgradient, &
                                                   point_response, atom_response, cell_response, weights_response, &
                                                   vdensity_response, vgradient_response, response)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      TYPE(qs_upf_inverse_augmentation_type), INTENT(IN) :: augmentation
      REAL(dp), INTENT(IN) :: vdensity(:), vgradient(:, :), point_response(:, :), &
         atom_response(:, :), cell_response(3, 3), weights_response(:), vdensity_response(:), &
         vgradient_response(:, :)
      COMPLEX(dp), INTENT(OUT)                           :: response(:, :)

      INTEGER                                            :: atom, axis, center, first, ikind, image, &
                                                            last, n, point
      REAL(dp)                                           :: dr(3), position(3), r(3), radius2
      REAL(dp), ALLOCATABLE                              :: dq(:, :, :), dterm(:, :), &
                                                            hq(:, :, :, :), q(:, :), term(:, :)

      CALL check_operator(frame, vdensity, vgradient, response)
      CPASSERT(ALL(SHAPE(point_response) == SHAPE(frame%points)))
      CPASSERT(ALL(SHAPE(atom_response) == SHAPE(frame%positions)))
      CPASSERT(SIZE(weights_response) == SIZE(frame%weights))
      CPASSERT(SIZE(vdensity_response) == SIZE(frame%weights))
      CPASSERT(ALL(SHAPE(vgradient_response) == SHAPE(vgradient)))
      response = CMPLX(0.0_dp, 0.0_dp, dp)
      DO center = 1, SIZE(frame%center_atoms)
         atom = frame%center_atoms(center)
         ikind = frame%center_kinds(center)
         first = frame%projector_first(center)
         last = frame%projector_last(center)
         n = last - first + 1
         radius2 = augmentation%kinds(ikind)%realspace_radius**2
         ALLOCATE (q(n, n), dq(n, n, 3), hq(n, n, 3, 3), term(n, n), dterm(n, n))
         DO image = 1, SIZE(frame%images, 2)
            position = frame%positions(:, atom) + MATMUL(frame%cell, REAL(frame%images(:, image), dp))
            DO point = 1, SIZE(frame%weights)
               r = frame%points(:, point) - position
               IF (DOT_PRODUCT(r, r) > radius2) CYCLE
               dr = point_response(:, point) - atom_response(:, atom) - MATMUL(cell_response, REAL(frame%images(:, image), dp))
               CALL uspp_augmentation_realspace(augmentation%kinds(ikind), r, q, dq, hq)
               term(:, :) = vdensity(point)*q
               dterm(:, :) = vdensity_response(point)*q
               DO axis = 1, 3
                  term(:, :) = term + vgradient(point, axis)*dq(:, :, axis)
                  dterm(:, :) = dterm + vdensity(point)*dq(:, :, axis)*dr(axis) + &
                                vgradient_response(point, axis)*dq(:, :, axis)
                  dterm(:, :) = dterm + vgradient(point, axis)*(hq(:, :, axis, 1)*dr(1) + &
                                                                hq(:, :, axis, 2)*dr(2) + hq(:, :, axis, 3)*dr(3))
               END DO
               response(first:last, first:last) = response(first:last, first:last) + &
                                                  frame%weights(point)*dterm + weights_response(point)*term
            END DO
         END DO
         DEALLOCATE (q, dq, hq, term, dterm)
      END DO
      CALL frame%point_comm%sum(response)
   END SUBROUTINE qs_upf_inverse_augmentation_response

! **************************************************************************************************
!> \brief ...
!> \param frame ...
!> \param vdensity ...
!> \param vgradient ...
!> \param projector_operator ...
! **************************************************************************************************
   SUBROUTINE check_operator(frame, vdensity, vgradient, projector_operator)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      REAL(dp), INTENT(IN)                               :: vdensity(:), vgradient(:, :)
      COMPLEX(dp), INTENT(IN)                            :: projector_operator(:, :)

      CPASSERT(frame%inverse%ready)
      CPASSERT(SIZE(vdensity) == SIZE(frame%weights))
      CPASSERT(ALL(SHAPE(vgradient) == [SIZE(frame%weights), 3]))
      CPASSERT(ALL(SHAPE(projector_operator) == frame%inverse%nproj))
   END SUBROUTINE check_operator

! **************************************************************************************************
!> \brief Release owned radial/angular caches.
!> \param augmentation ...
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_augmentation_release(augmentation)
      TYPE(qs_upf_inverse_augmentation_type), &
         INTENT(INOUT)                                   :: augmentation

      TYPE(qs_upf_inverse_augmentation_type)             :: empty

      augmentation = empty
   END SUBROUTINE qs_upf_inverse_augmentation_release

END MODULE qs_upf_inverse_augmentation
