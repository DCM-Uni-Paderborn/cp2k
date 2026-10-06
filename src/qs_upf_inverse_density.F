!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
!> \brief Home-center inverse density and its finite Gaussian GAPW angular map.
MODULE qs_upf_inverse_density
   USE gaussian_basis_fields,           ONLY: gaussian_basis_fields_periodic
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_comm_type
   USE qs_grid_atom,                    ONLY: grid_atom_type
   USE qs_harmonics_atom,               ONLY: harmonics_atom_type
   USE qs_upf_inverse,                  ONLY: qs_upf_inverse_center_fields,&
                                              qs_upf_inverse_evaluate
   USE qs_upf_inverse_types,            ONLY: qs_upf_inverse_augmentation_type,&
                                              qs_upf_inverse_type
   USE upf_inverse_quadrature,          ONLY: upf_inverse_cross_density,&
                                              upf_inverse_projector_density,&
                                              upf_inverse_semilocal_pair
   USE uspp_augmentation,               ONLY: uspp_augmentation_collocate_atom,&
                                              uspp_augmentation_integrate_atom
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: qs_upf_inverse_density_center, qs_upf_inverse_density_integrate
   PUBLIC :: qs_upf_inverse_density_radial, qs_upf_inverse_density_radial_adjoint

CONTAINS

! **************************************************************************************************
!> \brief Signed nR_a + Q_a from the full coupled inverse, without partition or image folding.
!>        nR_a = Re[d_a^* (s + chi) P], with d_a = -Delta_a W_a for the home image.
!>        Sum caller-supplied spin occupations and k weights outside this map exactly once.
!> \param frame Original inverse frame, including all centers and their coupled W.
!> \param augmentation Original charge kernels in the frame's magnetic-projector order.
!> \param center Selected home center, including a possible NC center with zero reconstruction/Q.
!> \param points Cartesian evaluation points, possibly an empty batch.
!> \param matrix Hermitian AO occupation density for one spin and k point.
!> \param density Signed home-center charge density, excluding Gaussian hard/soft and core.
!> \param response Optional joint density response at fixed matrix. All response inputs are required.
!> \param point_response Motion of the evaluation points.
!> \param atom_response Motion of every original atom.
!> \param cell_response Physical lattice response.
!> \param reduced_k_response Reduced Bloch-vector response.
!> \param coefficients_response Complete dW from the original defining inverse quadrature.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_density_center(frame, augmentation, center, points, matrix, density, response, &
                                            point_response, atom_response, cell_response, reduced_k_response, &
                                            coefficients_response)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      TYPE(qs_upf_inverse_augmentation_type), INTENT(IN) :: augmentation
      INTEGER, INTENT(IN)                                :: center
      REAL(dp), INTENT(IN)                               :: points(:, :)
      COMPLEX(dp), INTENT(IN)                            :: matrix(:, :)
      REAL(dp), INTENT(OUT)                              :: density(:)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: response(:)
      REAL(dp), INTENT(IN), OPTIONAL                     :: point_response(:, :), &
                                                            atom_response(:, :), &
                                                            cell_response(3, 3), &
                                                            reduced_k_response(3)
      COMPLEX(dp), INTENT(IN), OPTIONAL                  :: coefficients_response(:, :)

      COMPLEX(dp), ALLOCATABLE :: d(:, :), dd(:, :), de(:, :), e(:, :), gd(:, :, :), gdd(:, :, :), &
         gde(:, :, :), ge(:, :, :), projected(:, :), projected_response(:, :)
      INTEGER                                            :: atom, first, kind, last, n, np, nproj
      LOGICAL                                            :: derive
      REAL(dp), ALLOCATABLE                              :: dq(:), gradient(:, :), q(:), &
                                                            qgradient(:, :), tau(:)

      CALL check_center(frame, augmentation, center, points)
      np = SIZE(points, 2)
      n = SIZE(frame%coefficients, 2)
      nproj = SIZE(frame%coefficients, 1)
      CPASSERT(SIZE(density) == np)
      CPASSERT(ALL(SHAPE(matrix) == [n, n]))
      derive = PRESENT(response)
      CPASSERT(derive .EQV. PRESENT(point_response))
      CPASSERT(derive .EQV. PRESENT(atom_response))
      CPASSERT(derive .EQV. PRESENT(cell_response))
      CPASSERT(derive .EQV. PRESENT(reduced_k_response))
      CPASSERT(derive .EQV. PRESENT(coefficients_response))
      ALLOCATE (d(np, n), e(np, n), gd(np, n, 3), ge(np, n, 3), &
                q(np), qgradient(3, np), gradient(np, 3), tau(np), projected(nproj, nproj))
      IF (derive) THEN
         CPASSERT(SIZE(response) == np)
         ALLOCATE (dd(np, n), de(np, n), gdd(np, n, 3), gde(np, n, 3), dq(np))
         CALL center_fields(frame, center, points, d, e, gd, ge, dd, de, gdd, gde, &
                            point_response, atom_response, cell_response, reduced_k_response, coefficients_response)
         CALL upf_inverse_cross_density(dd, e, gdd, ge, matrix, response, gradient, tau)
         CALL upf_inverse_cross_density(d, de, gd, gde, matrix, dq, gradient, tau)
         response = response + dq
      ELSE
         CALL center_fields(frame, center, points, d, e, gd, ge)
      END IF
      CALL upf_inverse_cross_density(d, e, gd, ge, matrix, density, gradient, tau)
      atom = frame%center_atoms(center)
      kind = frame%center_kinds(center)
      IF (augmentation%kinds(kind)%nproj == 0) RETURN
      first = frame%projector_first(center)
      last = frame%projector_last(center)
      CALL upf_inverse_projector_density(frame%coefficients, matrix, projected)
      q = 0.0_dp
      qgradient = 0.0_dp
      CALL uspp_augmentation_collocate_atom(augmentation%kinds(kind), &
                                            points - SPREAD(frame%positions(:, atom), 2, np), &
                                            REAL(projected(first:last, first:last), dp), q, qgradient)
      density = density + q
      IF (derive) THEN
         response = response + SUM(qgradient*(point_response - SPREAD(atom_response(:, atom), 2, np)), DIM=1)
         ALLOCATE (projected_response(nproj, nproj))
         projected_response(:, :) = MATMUL(coefficients_response, MATMUL(matrix, CONJG(TRANSPOSE(frame%coefficients))))
         projected_response(:, :) = projected_response + CONJG(TRANSPOSE(projected_response))
         q = 0.0_dp
         qgradient = 0.0_dp
         CALL uspp_augmentation_collocate_atom(augmentation%kinds(kind), &
                                               points - SPREAD(frame%positions(:, atom), 2, np), &
                                               REAL(projected_response(first:last, first:last), dp), q, qgradient)
         response = response + q
      END IF
   END SUBROUTINE qs_upf_inverse_density_center

! **************************************************************************************************
!> \brief Density adjoint with both inverse field legs and the original W^dagger Q W operator.
!>        A unit potential gives this center's contribution to the finite integrated charge metric.
!> \param frame ...
!> \param augmentation ...
!> \param center ...
!> \param points ...
!> \param weights Point measure, applied once.
!> \param potential Unweighted real density covector.
!> \param operator Hermitian AO adjoint for one spin and k point, without a k weight.
!> \param comm Optional communicator distributing point rows. Replicated atom grids omit it.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_density_integrate(frame, augmentation, center, points, weights, potential, OPERATOR, comm)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      TYPE(qs_upf_inverse_augmentation_type), INTENT(IN) :: augmentation
      INTEGER, INTENT(IN)                                :: center
      REAL(dp), INTENT(IN)                               :: points(:, :), weights(:), potential(:)
      COMPLEX(dp), INTENT(OUT)                           :: operator(:, :)
      TYPE(mp_comm_type), INTENT(IN), OPTIONAL           :: comm

      COMPLEX(dp), ALLOCATABLE                           :: d(:, :), e(:, :), gd(:, :, :), &
                                                            ge(:, :, :), qoperator(:, :)
      INTEGER                                            :: atom, first, kind, last, n, np
      REAL(dp), ALLOCATABLE                              :: q(:, :), zero(:), zeros(:, :)

      CALL check_center(frame, augmentation, center, points)
      np = SIZE(points, 2)
      n = SIZE(frame%coefficients, 2)
      CPASSERT(SIZE(weights) == np .AND. SIZE(potential) == np)
      CPASSERT(ALL(SHAPE(OPERATOR) == [n, n]))
      ALLOCATE (d(np, n), e(np, n), gd(np, n, 3), ge(np, n, 3), zero(np), zeros(np, 3))
      CALL center_fields(frame, center, points, d, e, gd, ge)
      zero = 0.0_dp
      zeros = 0.0_dp
      CALL upf_inverse_semilocal_pair(d, e, gd, ge, weights, potential, zeros, zero, OPERATOR)
      OPERATOR = 0.5_dp*(OPERATOR + CONJG(TRANSPOSE(OPERATOR)))
      atom = frame%center_atoms(center)
      kind = frame%center_kinds(center)
      IF (augmentation%kinds(kind)%nproj /= 0) THEN
         first = frame%projector_first(center)
         last = frame%projector_last(center)
         ALLOCATE (q(last - first + 1, last - first + 1), qoperator(last - first + 1, last - first + 1))
         CALL uspp_augmentation_integrate_atom(augmentation%kinds(kind), &
                                               points - SPREAD(frame%positions(:, atom), 2, np), &
                                               weights, potential, TRANSPOSE(zeros), q)
         qoperator(:, :) = CMPLX(q, 0.0_dp, dp)
         OPERATOR = OPERATOR + MATMUL(CONJG(TRANSPOSE(frame%coefficients(first:last, :))), &
                                      MATMUL(qoperator, frame%coefficients(first:last, :)))
      END IF
      IF (PRESENT(comm)) CALL comm%sum(OPERATOR)
   END SUBROUTINE qs_upf_inverse_density_integrate

! **************************************************************************************************
!> \brief Native angular projection of nR + Q and optional Gaussian hard-minus-soft density.
!>        Use the same finite Lebedev map for every leg. Its Gram matrix is not assumed to be unity.
!>        The native harmonic order is l=0,1,... and m=-l,...,l, as in the Coulomb maps.
!> \param grid Borrowed native radial/angular quadrature, with complete replicated angular rows.
!> \param harmonics Borrowed native real harmonics on that same angular quadrature.
!> \param values Point density with angular rows fastest. No core or partition is added here.
!> \param radial New real harmonic density coefficients, without the radial measure.
!> \param gaussian_hard Optional native Gaussian hard harmonic coefficients.
!> \param gaussian_soft Matching native Gaussian soft harmonic coefficients.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_density_radial(grid, harmonics, values, radial, gaussian_hard, gaussian_soft)
      TYPE(grid_atom_type), INTENT(IN)                   :: grid
      TYPE(harmonics_atom_type), INTENT(IN)              :: harmonics
      REAL(dp), INTENT(IN)                               :: values(:)
      REAL(dp), INTENT(OUT)                              :: radial(:, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: gaussian_hard(:, :), gaussian_soft(:, :)

      INTEGER                                            :: first, ir, last, na, nh
      REAL(dp), ALLOCATABLE                              :: angular(:)

      na = grid%ng_sphere
      nh = SIZE(radial, 2)
      CALL check_radial(grid, harmonics, values, radial)
      CPASSERT(PRESENT(gaussian_hard) .EQV. PRESENT(gaussian_soft))
      IF (PRESENT(gaussian_hard)) THEN
         CPASSERT(ALL(SHAPE(gaussian_hard) == SHAPE(radial)))
         CPASSERT(ALL(SHAPE(gaussian_soft) == SHAPE(radial)))
      END IF
      ALLOCATE (angular(na))
      DO ir = 1, grid%nr
         first = (ir - 1)*na + 1
         last = ir*na
         angular(:) = values(first:last)
         IF (PRESENT(gaussian_hard)) THEN
            angular(:) = angular + MATMUL(harmonics%slm(:, 1:nh), gaussian_hard(ir, :) - gaussian_soft(ir, :))
         END IF
         radial(ir, :) = MATMUL(grid%wa*angular, harmonics%slm(:, 1:nh))
      END DO
   END SUBROUTINE qs_upf_inverse_density_radial

! **************************************************************************************************
!> \brief Exact transpose of the native finite angular projection for all three density legs.
!> \param grid ...
!> \param harmonics ...
!> \param potential Unweighted radial-measure covector, as returned by the Coulomb assembly.
!> \param values Unweighted point covector. Multiply by grid%weight once in point/AO integration.
!> \param gaussian_hard Optional native hard radial-measure covector, with the angular Gram map.
!> \param gaussian_soft Matching native soft radial-measure covector with opposite sign.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_density_radial_adjoint(grid, harmonics, potential, values, gaussian_hard, gaussian_soft)
      TYPE(grid_atom_type), INTENT(IN)                   :: grid
      TYPE(harmonics_atom_type), INTENT(IN)              :: harmonics
      REAL(dp), INTENT(IN)                               :: potential(:, :)
      REAL(dp), INTENT(OUT)                              :: values(:)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: gaussian_hard(:, :), gaussian_soft(:, :)

      INTEGER                                            :: first, ir, last, na, nh

      CALL check_radial(grid, harmonics, values, potential)
      CPASSERT(PRESENT(gaussian_hard) .EQV. PRESENT(gaussian_soft))
      IF (PRESENT(gaussian_hard)) THEN
         CPASSERT(ALL(SHAPE(gaussian_hard) == SHAPE(potential)))
         CPASSERT(ALL(SHAPE(gaussian_soft) == SHAPE(potential)))
      END IF
      na = grid%ng_sphere
      nh = SIZE(potential, 2)
      DO ir = 1, grid%nr
         first = (ir - 1)*na + 1
         last = ir*na
         values(first:last) = MATMUL(harmonics%slm(:, 1:nh), potential(ir, :))
         IF (PRESENT(gaussian_hard)) THEN
            gaussian_hard(ir, :) = MATMUL(grid%wa*values(first:last), harmonics%slm(:, 1:nh))
            gaussian_soft(ir, :) = -gaussian_hard(ir, :)
         END IF
      END DO
   END SUBROUTINE qs_upf_inverse_density_radial_adjoint

! **************************************************************************************************
!> \brief Full environmental field s+chi and a single home-center difference, including responses.
!> \param frame ...
!> \param center ...
!> \param points ...
!> \param d ...
!> \param e ...
!> \param gd ...
!> \param ge ...
!> \param dd ...
!> \param de ...
!> \param gdd ...
!> \param gde ...
!> \param point_response ...
!> \param atom_response ...
!> \param cell_response ...
!> \param reduced_k_response ...
!> \param coefficients_response ...
! **************************************************************************************************
   SUBROUTINE center_fields(frame, center, points, d, e, gd, ge, dd, de, gdd, gde, &
                            point_response, atom_response, cell_response, reduced_k_response, coefficients_response)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      INTEGER, INTENT(IN)                                :: center
      REAL(dp), INTENT(IN)                               :: points(:, :)
      COMPLEX(dp), INTENT(OUT)                           :: d(:, :), e(:, :), gd(:, :, :), &
                                                            ge(:, :, :)
      COMPLEX(dp), INTENT(OUT), OPTIONAL                 :: dd(:, :), de(:, :), gdd(:, :, :), &
                                                            gde(:, :, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: point_response(:, :), &
                                                            atom_response(:, :), &
                                                            cell_response(3, 3), &
                                                            reduced_k_response(3)
      COMPLEX(dp), INTENT(IN), OPTIONAL                  :: coefficients_response(:, :)

      COMPLEX(dp), ALLOCATABLE                           :: chi(:, :), dchi(:, :), dgc(:, :, :), &
                                                            gc(:, :, :)
      INTEGER                                            :: atom, first, last, n, np

      n = SIZE(frame%coefficients, 2)
      np = SIZE(points, 2)
      ALLOCATE (chi(np, n), gc(np, n, 3))
      CALL qs_upf_inverse_center_fields(frame, center, points, d, gd, dd, gdd, point_response, &
                                        atom_response, cell_response, reduced_k_response, coefficients_response, [0, 0, 0])
      CALL qs_upf_inverse_evaluate(frame, points, e, ge, de, gde, point_response, atom_response, &
                                   cell_response, reduced_k_response, coefficients_response)
      IF (PRESENT(de)) ALLOCATE (dchi(np, n), dgc(np, n, 3))
      DO atom = 1, SIZE(frame%bases)
         first = frame%ao_first(atom)
         last = frame%ao_last(atom)
         IF (PRESENT(de)) THEN
            CALL gaussian_basis_fields_periodic(frame%bases(atom)%gto_basis_set, points, frame%positions(:, atom), &
                                                frame%cell, frame%images, frame%reduced_k, &
                                                chi(:, first:last), gc(:, first:last, :), &
                                                dchi(:, first:last), dgc(:, first:last, :), &
                                                point_response, atom_response(:, atom), &
                                                cell_response, reduced_k_response)
         ELSE
            CALL gaussian_basis_fields_periodic(frame%bases(atom)%gto_basis_set, points, frame%positions(:, atom), &
                                                frame%cell, frame%images, frame%reduced_k, chi(:, first:last), gc(:, first:last, :))
         END IF
      END DO
      e = e + chi
      ge = ge + gc
      IF (PRESENT(de)) THEN
         de = de + dchi
         gde = gde + dgc
      END IF
   END SUBROUTINE center_fields

! **************************************************************************************************
!> \brief Verify the original center and magnetic charge-kernel order.
!> \param frame ...
!> \param augmentation ...
!> \param center ...
!> \param points ...
! **************************************************************************************************
   SUBROUTINE check_center(frame, augmentation, center, points)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      TYPE(qs_upf_inverse_augmentation_type), INTENT(IN) :: augmentation
      INTEGER, INTENT(IN)                                :: center
      REAL(dp), INTENT(IN)                               :: points(:, :)

      INTEGER                                            :: kind, nproj

      CPASSERT(frame%inverse%ready)
      CPASSERT(center >= 1 .AND. center <= SIZE(frame%center_atoms))
      CPASSERT(SIZE(points, 1) == 3)
      kind = frame%center_kinds(center)
      CPASSERT(kind >= 1 .AND. kind <= SIZE(augmentation%kinds))
      nproj = frame%projector_last(center) - frame%projector_first(center) + 1
      IF (augmentation%kinds(kind)%nproj /= 0) THEN
         CPASSERT(augmentation%kinds(kind)%nproj == nproj)
      END IF
   END SUBROUTINE check_center

! **************************************************************************************************
!> \brief Validate one complete replicated radial/angular field in native harmonic order.
!> \param grid ...
!> \param harmonics ...
!> \param values ...
!> \param radial ...
! **************************************************************************************************
   SUBROUTINE check_radial(grid, harmonics, values, radial)
      TYPE(grid_atom_type), INTENT(IN)                   :: grid
      TYPE(harmonics_atom_type), INTENT(IN)              :: harmonics
      REAL(dp), INTENT(IN)                               :: values(:), radial(:, :)

      CPASSERT(grid%nr >= 0 .AND. grid%ng_sphere > 0)
      CPASSERT(ASSOCIATED(grid%wa) .AND. ASSOCIATED(harmonics%slm))
      CPASSERT(SIZE(values) == grid%nr*grid%ng_sphere)
      CPASSERT(SIZE(radial, 1) == grid%nr)
      CPASSERT(SIZE(grid%wa) == grid%ng_sphere)
      CPASSERT(SIZE(harmonics%slm, 1) == grid%ng_sphere)
      CPASSERT(SIZE(radial, 2) <= SIZE(harmonics%slm, 2))
   END SUBROUTINE check_radial

END MODULE qs_upf_inverse_density
