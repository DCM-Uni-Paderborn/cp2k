!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
!> \brief Native Gaussian one-center density map and its radial-measure adjoint.
MODULE qs_upf_inverse_gaussian
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
   USE gaunt_coefficients,              ONLY: real_gaunt_coefficients
   USE kinds,                           ONLY: dp
   USE orbital_pointers,                ONLY: indso,&
                                              nsoset
   USE paw_basis_types,                 ONLY: get_paw_basis_info
   USE qs_grid_atom,                    ONLY: grid_atom_type
   USE qs_harmonics_atom,               ONLY: gaunt_threshold,&
                                              harmonics_atom_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: qs_upf_inverse_gaussian_radial, qs_upf_inverse_gaussian_integrate

CONTAINS

! **************************************************************************************************
!> \brief Map the actual compressed native CPC coefficients to Gaussian harmonic densities.
!>        Preserve the native Gaunt coefficients and product-angular-momentum truncation.
!>        Apply separately to hard and soft CPC matrices before the finite angular composition.
!> \param kind Borrowed native GAPW kind, including its one-center basis, grid and harmonics.
!> \param coefficients Native compressed primitive density matrix, without a radial measure.
!> \param radial Harmonic density coefficients in native l,m order, without a radial measure.
!> \param response Optional joint coefficient and radial-coordinate response.
!> \param coefficients_response Variation of the same compressed density matrix.
!> \param radius_response Variation of the physical radial nodes, excluding quadrature weights.
!> \param right_kind Optional independent column primitive basis. The quadrature remains that of kind.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_gaussian_radial(kind, coefficients, radial, response, &
                                             coefficients_response, radius_response, right_kind)
      TYPE(qs_kind_type), INTENT(IN)                     :: kind
      REAL(dp), INTENT(IN)                               :: coefficients(:, :)
      REAL(dp), INTENT(OUT)                              :: radial(:, :)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: response(:, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: coefficients_response(:, :), &
                                                            radius_response(:)
      TYPE(qs_kind_type), INTENT(IN), OPTIONAL           :: right_kind

      INTEGER                                            :: h, i, j, l, n, nh, nr
      INTEGER, ALLOCATABLE                               :: angular(:), right_angular(:)
      LOGICAL                                            :: derive, independent
      REAL(dp)                                           :: coupling, exponent
      REAL(dp), ALLOCATABLE                              :: derivative(:), exponents(:), PRODUCT(:), &
                                                            right_exponents(:)
      REAL(dp), POINTER                                  :: products(:, :, :)
      TYPE(grid_atom_type), POINTER                      :: grid
      TYPE(harmonics_atom_type), POINTER                 :: harmonics

      CALL gaussian_product_layout(kind, grid, harmonics, angular, exponents, &
                                   right_angular, right_exponents, products, independent, right_kind)
      n = SIZE(angular)
      nr = SIZE(right_angular)
      nh = harmonics%max_iso_not0
      CALL check_shapes(grid, nh, n, nr, coefficients, radial)
      derive = PRESENT(response)
      CPASSERT(derive .EQV. PRESENT(coefficients_response))
      CPASSERT(derive .EQV. PRESENT(radius_response))
      ALLOCATE (PRODUCT(grid%nr))
      radial = 0.0_dp
      IF (derive) THEN
         CALL check_shapes(grid, nh, n, nr, coefficients_response, response)
         CPASSERT(SIZE(radius_response) == grid%nr)
         ALLOCATE (derivative(grid%nr))
         response = 0.0_dp
      END IF
      DO j = 1, nr
         DO i = 1, n
            l = indso(1, angular(i)) + indso(1, right_angular(j))
            ! Truncate the output harmonics, not the radial power of their product.
            exponent = exponents(i) + right_exponents(j)
            PRODUCT(:) = grid%rad**l*EXP(-exponent*grid%rad**2)
            IF (derive) THEN
               derivative(:) = -2.0_dp*exponent*grid%rad*product
               IF (l > 0) derivative(:) = derivative + &
                                          REAL(l, dp)*grid%rad**(l - 1)*EXP(-exponent*grid%rad**2)
               derivative(:) = derivative*radius_response
            END IF
            DO h = 1, nh
               coupling = products(angular(i), right_angular(j), h)
               IF (ABS(coupling) <= gaunt_threshold) CYCLE
               radial(:, h) = radial(:, h) + coupling*coefficients(i, j)*product
               IF (derive) response(:, h) = response(:, h) + &
                                            coupling*(coefficients_response(i, j)*product + coefficients(i, j)*derivative)
            END DO
         END DO
      END DO
      IF (independent) DEALLOCATE (products)
   END SUBROUTINE qs_upf_inverse_gaussian_radial

! **************************************************************************************************
!> \brief Exact density transpose into the same compressed native CPC coordinates.
!>        The caller supplies the harmonic covector after the finite angular transpose.
!>        A unit point potential therefore uses the finite angular Gram covector of Y_00.
!> \param kind ...
!> \param potential Unweighted harmonic-density covector, with no radial measure included.
!> \param adjoint New compressed primitive adjoint. Native radial weights are applied once.
!> \param right_kind Optional independent column primitive basis, using the same quadrature as kind.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_gaussian_integrate(kind, potential, adjoint, right_kind)
      TYPE(qs_kind_type), INTENT(IN)                     :: kind
      REAL(dp), INTENT(IN)                               :: potential(:, :)
      REAL(dp), INTENT(OUT)                              :: adjoint(:, :)
      TYPE(qs_kind_type), INTENT(IN), OPTIONAL           :: right_kind

      INTEGER                                            :: h, i, j, l, n, nh, nr
      INTEGER, ALLOCATABLE                               :: angular(:), right_angular(:)
      LOGICAL                                            :: independent
      REAL(dp)                                           :: coupling
      REAL(dp), ALLOCATABLE                              :: exponents(:), PRODUCT(:), &
                                                            right_exponents(:)
      REAL(dp), POINTER                                  :: products(:, :, :)
      TYPE(grid_atom_type), POINTER                      :: grid
      TYPE(harmonics_atom_type), POINTER                 :: harmonics

      CALL gaussian_product_layout(kind, grid, harmonics, angular, exponents, &
                                   right_angular, right_exponents, products, independent, right_kind)
      n = SIZE(angular)
      nr = SIZE(right_angular)
      nh = harmonics%max_iso_not0
      CALL check_shapes(grid, nh, n, nr, adjoint, potential)
      CPASSERT(ASSOCIATED(grid%wr))
      CPASSERT(SIZE(grid%wr) == grid%nr)
      ALLOCATE (PRODUCT(grid%nr))
      adjoint = 0.0_dp
      DO j = 1, nr
         DO i = 1, n
            l = indso(1, angular(i)) + indso(1, right_angular(j))
            PRODUCT(:) = grid%wr*grid%rad**l*EXP(-(exponents(i) + right_exponents(j))*grid%rad**2)
            DO h = 1, nh
               coupling = products(angular(i), right_angular(j), h)
               IF (ABS(coupling) <= gaunt_threshold) CYCLE
               adjoint(i, j) = adjoint(i, j) + coupling*SUM(product*potential(:, h))
            END DO
         END DO
      END DO
      IF (independent) DEALLOCATE (products)
   END SUBROUTINE qs_upf_inverse_gaussian_integrate

! **************************************************************************************************
!> \brief Independent primitive layouts on one explicitly selected native radial/angular quadrature.
!> \param kind Row primitive basis and owner of the integration grid.
!> \param grid Borrowed native radial grid.
!> \param harmonics Borrowed output harmonics and their native angular truncation.
!> \param angular Row spherical-harmonic indices.
!> \param exponents Row primitive exponents.
!> \param right_angular Column spherical-harmonic indices.
!> \param right_exponents Column primitive exponents.
!> \param products Native Gaunt coefficients, borrowed for a square map and owned for a cross map.
!> \param independent True when products must be released by the caller.
!> \param right_kind Optional independent column owner. Its radial grid is not substituted.
! **************************************************************************************************
   SUBROUTINE gaussian_product_layout(kind, grid, harmonics, angular, exponents, &
                                      right_angular, right_exponents, products, independent, right_kind)
      TYPE(qs_kind_type), INTENT(IN)                     :: kind
      TYPE(grid_atom_type), POINTER                      :: grid
      TYPE(harmonics_atom_type), POINTER                 :: harmonics
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: angular(:)
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: exponents(:)
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: right_angular(:)
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: right_exponents(:)
      REAL(dp), POINTER                                  :: products(:, :, :)
      LOGICAL, INTENT(OUT)                               :: independent
      TYPE(qs_kind_type), INTENT(IN), OPTIONAL           :: right_kind

      TYPE(grid_atom_type), POINTER                      :: right_grid
      TYPE(harmonics_atom_type), POINTER                 :: right_harmonics

      CALL gaussian_layout(kind, grid, harmonics, angular, exponents)
      independent = PRESENT(right_kind)
      IF (independent) THEN
         CALL gaussian_layout(right_kind, right_grid, right_harmonics, right_angular, right_exponents)
         CALL real_gaunt_coefficients(products, MAXVAL(indso(1, angular)), &
                                      MAXVAL(indso(1, right_angular)), indso(1, harmonics%max_iso_not0))
      ELSE
         ALLOCATE (right_angular, SOURCE=angular)
         ALLOCATE (right_exponents, SOURCE=exponents)
         products => harmonics%my_CG
      END IF
   END SUBROUTINE gaussian_product_layout

! **************************************************************************************************
!> \brief Decode the native compressed primitive order, including sets with nonzero lmin.
!> \param kind ...
!> \param grid Borrowed radial grid.
!> \param harmonics Borrowed native Gaunt coefficients.
!> \param angular Owned spherical-harmonic index of each compressed primitive.
!> \param exponents Owned Gaussian exponent of each compressed primitive.
! **************************************************************************************************
   SUBROUTINE gaussian_layout(kind, grid, harmonics, angular, exponents)
      TYPE(qs_kind_type), INTENT(IN)                     :: kind
      TYPE(grid_atom_type), POINTER                      :: grid
      TYPE(harmonics_atom_type), POINTER                 :: harmonics
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: angular(:)
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: exponents(:)

      INTEGER                                            :: i, ipgf, iset, iso, n, nset
      INTEGER, POINTER                                   :: lmax(:), lmin(:), npgf(:)
      REAL(dp), POINTER                                  :: zet(:, :)
      TYPE(gto_basis_set_type), POINTER                  :: basis

      CALL get_qs_kind(kind, basis_set=basis, basis_type="GAPW_1C", grid_atom=grid, harmonics=harmonics)
      CPASSERT(ASSOCIATED(basis))
      CPASSERT(ASSOCIATED(grid))
      CPASSERT(ASSOCIATED(harmonics))
      CPASSERT(ASSOCIATED(harmonics%my_CG))
      CPASSERT(grid%nr >= 0)
      CPASSERT(harmonics%max_iso_not0 > 0)
      CALL get_gto_basis_set(basis, nset=nset, lmax=lmax, lmin=lmin, npgf=npgf, zet=zet)
      CALL get_paw_basis_info(basis, nsatbas=n)
      ALLOCATE (angular(n), exponents(n))
      i = 0
      DO iset = 1, nset
         DO ipgf = 1, npgf(iset)
            DO iso = nsoset(lmin(iset) - 1) + 1, nsoset(lmax(iset))
               i = i + 1
               angular(i) = iso
               exponents(i) = zet(ipgf, iset)
            END DO
         END DO
      END DO
      CPASSERT(i == n)
   END SUBROUTINE gaussian_layout

! **************************************************************************************************
!> \brief ...
!> \param grid ...
!> \param nh ...
!> \param n ...
!> \param nright ...
!> \param coefficients ...
!> \param radial ...
! **************************************************************************************************
   SUBROUTINE check_shapes(grid, nh, n, nright, coefficients, radial)
      TYPE(grid_atom_type), INTENT(IN)                   :: grid
      INTEGER, INTENT(IN)                                :: nh, n, nright
      REAL(dp), INTENT(IN)                               :: coefficients(:, :), radial(:, :)

      CPASSERT(ALL(SHAPE(coefficients) == [n, nright]))
      CPASSERT(ALL(SHAPE(radial) == [grid%nr, nh]))
   END SUBROUTINE check_shapes

END MODULE qs_upf_inverse_gaussian
