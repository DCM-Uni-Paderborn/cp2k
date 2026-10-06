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
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_gaussian_radial(kind, coefficients, radial, response, &
                                             coefficients_response, radius_response)
      TYPE(qs_kind_type), INTENT(IN)                     :: kind
      REAL(dp), INTENT(IN)                               :: coefficients(:, :)
      REAL(dp), INTENT(OUT)                              :: radial(:, :)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: response(:, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: coefficients_response(:, :), &
                                                            radius_response(:)

      INTEGER                                            :: h, i, j, l, maxl, n, nh
      INTEGER, ALLOCATABLE                               :: angular(:)
      LOGICAL                                            :: derive
      REAL(dp)                                           :: coupling, exponent
      REAL(dp), ALLOCATABLE                              :: derivative(:), exponents(:), PRODUCT(:)
      TYPE(grid_atom_type), POINTER                      :: grid
      TYPE(harmonics_atom_type), POINTER                 :: harmonics

      CALL gaussian_layout(kind, grid, harmonics, angular, exponents)
      n = SIZE(angular)
      nh = harmonics%max_iso_not0
      CALL check_shapes(grid, nh, n, coefficients, radial)
      derive = PRESENT(response)
      CPASSERT(derive .EQV. PRESENT(coefficients_response))
      CPASSERT(derive .EQV. PRESENT(radius_response))
      ALLOCATE (PRODUCT(grid%nr))
      radial = 0.0_dp
      IF (derive) THEN
         CALL check_shapes(grid, nh, n, coefficients_response, response)
         CPASSERT(SIZE(radius_response) == grid%nr)
         ALLOCATE (derivative(grid%nr))
         response = 0.0_dp
      END IF
      maxl = indso(1, nh)
      DO j = 1, n
         DO i = 1, n
            l = indso(1, angular(i)) + indso(1, angular(j))
            IF (l > maxl) CYCLE
            exponent = exponents(i) + exponents(j)
            PRODUCT(:) = grid%rad**l*EXP(-exponent*grid%rad**2)
            IF (derive) THEN
               derivative(:) = -2.0_dp*exponent*grid%rad*product
               IF (l > 0) derivative(:) = derivative + &
                                          REAL(l, dp)*grid%rad**(l - 1)*EXP(-exponent*grid%rad**2)
               derivative(:) = derivative*radius_response
            END IF
            DO h = 1, nh
               coupling = harmonics%my_CG(angular(i), angular(j), h)
               IF (ABS(coupling) <= gaunt_threshold) CYCLE
               radial(:, h) = radial(:, h) + coupling*coefficients(i, j)*product
               IF (derive) response(:, h) = response(:, h) + &
                                            coupling*(coefficients_response(i, j)*product + coefficients(i, j)*derivative)
            END DO
         END DO
      END DO
   END SUBROUTINE qs_upf_inverse_gaussian_radial

! **************************************************************************************************
!> \brief Exact density transpose into the same compressed native CPC coordinates.
!>        The caller supplies the harmonic covector after the finite angular transpose.
!>        A unit point potential therefore uses the finite angular Gram covector of Y_00.
!> \param kind ...
!> \param potential Unweighted harmonic-density covector, with no radial measure included.
!> \param adjoint New compressed primitive adjoint. Native radial weights are applied once.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_gaussian_integrate(kind, potential, adjoint)
      TYPE(qs_kind_type), INTENT(IN)                     :: kind
      REAL(dp), INTENT(IN)                               :: potential(:, :)
      REAL(dp), INTENT(OUT)                              :: adjoint(:, :)

      INTEGER                                            :: h, i, j, l, maxl, n, nh
      INTEGER, ALLOCATABLE                               :: angular(:)
      REAL(dp)                                           :: coupling
      REAL(dp), ALLOCATABLE                              :: exponents(:), PRODUCT(:)
      TYPE(grid_atom_type), POINTER                      :: grid
      TYPE(harmonics_atom_type), POINTER                 :: harmonics

      CALL gaussian_layout(kind, grid, harmonics, angular, exponents)
      n = SIZE(angular)
      nh = harmonics%max_iso_not0
      CALL check_shapes(grid, nh, n, adjoint, potential)
      CPASSERT(ASSOCIATED(grid%wr))
      CPASSERT(SIZE(grid%wr) == grid%nr)
      ALLOCATE (PRODUCT(grid%nr))
      adjoint = 0.0_dp
      maxl = indso(1, nh)
      DO j = 1, n
         DO i = 1, n
            l = indso(1, angular(i)) + indso(1, angular(j))
            IF (l > maxl) CYCLE
            PRODUCT(:) = grid%wr*grid%rad**l*EXP(-(exponents(i) + exponents(j))*grid%rad**2)
            DO h = 1, nh
               coupling = harmonics%my_CG(angular(i), angular(j), h)
               IF (ABS(coupling) <= gaunt_threshold) CYCLE
               adjoint(i, j) = adjoint(i, j) + coupling*SUM(product*potential(:, h))
            END DO
         END DO
      END DO
   END SUBROUTINE qs_upf_inverse_gaussian_integrate

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
!> \param coefficients ...
!> \param radial ...
! **************************************************************************************************
   SUBROUTINE check_shapes(grid, nh, n, coefficients, radial)
      TYPE(grid_atom_type), INTENT(IN)                   :: grid
      INTEGER, INTENT(IN)                                :: nh, n
      REAL(dp), INTENT(IN)                               :: coefficients(:, :), radial(:, :)

      CPASSERT(ALL(SHAPE(coefficients) == [n, n]))
      CPASSERT(ALL(SHAPE(radial) == [grid%nr, nh]))
   END SUBROUTINE check_shapes

END MODULE qs_upf_inverse_gaussian
