!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Full Coulomb potentials of native contracted Gaussian functions and their point gradients.
!>        The native Cartesian contraction is decomposed into Gaussian solid harmonics with radial
!>        powers. The original normalization and all multipole tails are retained. Owned caches
!>        can be evaluated without a live basis pointer, including at the Gaussian center.
! **************************************************************************************************
MODULE gaussian_coulomb_fields
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: dfac,&
                                              fourpi,&
                                              pi
   USE orbital_pointers,                ONLY: indco,&
                                              init_orbital_pointers,&
                                              nco,&
                                              ncoset
   USE orbital_transformation_matrices, ONLY: init_spherical_harmonics,&
                                              orbtramat
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   TYPE gaussian_coulomb_primitive_type
      REAL(dp) :: alpha = 0.0_dp
      INTEGER :: l = 0, power = 0
      REAL(dp), ALLOCATABLE :: coefficients(:, :)
   END TYPE gaussian_coulomb_primitive_type

   TYPE gaussian_coulomb_fields_type
      INTEGER :: nfunction = 0
      TYPE(gaussian_coulomb_primitive_type), ALLOCATABLE :: primitives(:)
   END TYPE gaussian_coulomb_fields_type

   PUBLIC :: gaussian_coulomb_fields_type, gaussian_coulomb_fields_init
   PUBLIC :: gaussian_coulomb_fields_prepare, gaussian_coulomb_fields_evaluate
   PUBLIC :: gaussian_coulomb_fields_screen
   PUBLIC :: gaussian_coulomb_fields_periodic

CONTAINS

! **************************************************************************************************
!> \brief Initialize shared angular tables before preparing/evaluating caches in parallel.
!> \param maxl Largest native Cartesian polynomial degree in the requested Gaussian bases.
! **************************************************************************************************
   SUBROUTINE gaussian_coulomb_fields_init(maxl)
      INTEGER, INTENT(IN)                                :: maxl

      CPASSERT(maxl >= 0)
      CALL init_orbital_pointers(maxl)
      CALL init_spherical_harmonics(maxl, -1)
   END SUBROUTINE gaussian_coulomb_fields_init

! **************************************************************************************************
!> \brief Copy one native shell-set contraction into an owned Coulomb-potential cache.
!>        sphi includes the native Cartesian primitive normalization and contraction coefficients.
!>        No fit, screening radius or compensation subtraction is introduced.
!> \param basis Native contracted Gaussian basis.
!> \param iset Native shell-set index.
!> \param fields New owned cache for the functions in this shell set.
! **************************************************************************************************
   SUBROUTINE gaussian_coulomb_fields_prepare(basis, iset, fields)
      TYPE(gto_basis_set_type), INTENT(IN)               :: basis
      INTEGER, INTENT(IN)                                :: iset
      TYPE(gaussian_coulomb_fields_type), INTENT(OUT)    :: fields

      INTEGER                                            :: count, degree, first, ip, l, last, m, &
                                                            nc, nset, row, t
      INTEGER, POINTER                                   :: lmax(:), npgf(:), nsgf_set(:)
      REAL(dp), POINTER                                  :: sphi(:, :), zet(:, :)

      CALL get_gto_basis_set(basis, nset=nset, lmax=lmax, npgf=npgf, nsgf_set=nsgf_set, sphi=sphi, zet=zet)
      CPASSERT(iset >= 1 .AND. iset <= nset)
      CPASSERT(ASSOCIATED(lmax))
      CPASSERT(ASSOCIATED(npgf))
      CPASSERT(ASSOCIATED(nsgf_set))
      CPASSERT(ASSOCIATED(sphi) .AND. ASSOCIATED(zet))
      CPASSERT(lmax(iset) >= 0 .AND. npgf(iset) > 0)
      CPASSERT(ASSOCIATED(orbtramat))
      CPASSERT(UBOUND(orbtramat, 1) >= lmax(iset))
      fields%nfunction = nsgf_set(iset)
      nc = ncoset(lmax(iset))
      first = SUM(nsgf_set(:iset - 1)) + 1
      last = first + fields%nfunction - 1
      CPASSERT(fields%nfunction > 0)
      CPASSERT(SIZE(sphi, 1) >= nc*npgf(iset) .AND. SIZE(sphi, 2) >= last)
      CPASSERT(SIZE(zet, 1) >= npgf(iset) .AND. SIZE(zet, 2) >= iset)
      CPASSERT(ALL(ieee_is_finite(sphi(:nc*npgf(iset), first:last))))
      CPASSERT(ALL(ieee_is_finite(zet(:npgf(iset), iset))))
      CPASSERT(ALL(zet(:npgf(iset), iset) > 0.0_dp))
      count = 0
      DO degree = 0, lmax(iset)
         count = count + degree/2 + 1
      END DO
      ALLOCATE (fields%primitives(count*npgf(iset)))
      count = 0
      DO ip = 1, npgf(iset)
         DO degree = 0, lmax(iset)
            DO l = MODULO(degree, 2), degree, 2
               count = count + 1
               ASSOCIATE (primitive => fields%primitives(count))
                  primitive%alpha = zet(ip, iset)
                  primitive%l = l
                  primitive%power = (degree - l)/2
                  ALLOCATE (primitive%coefficients(2*l + 1, fields%nfunction), SOURCE=0.0_dp)
                  DO t = ncoset(degree - 1) + 1, ncoset(degree)
                     row = (ip - 1)*nc + t
                     DO m = 1, 2*l + 1
                        primitive%coefficients(m, :) = primitive%coefficients(m, :) + &
                                                       monomial_harmonic(indco(:, t), l, m)*sphi(row, first:last)
                     END DO
                  END DO
               END ASSOCIATE
            END DO
         END DO
      END DO
   END SUBROUTINE gaussian_coulomb_fields_prepare

! **************************************************************************************************
!> \brief Apply erf(omega*r)/r or erfc(omega*r)/r to an owned native Gaussian cache.
!>        The long-range field is the Coulomb field of the exact Gaussian convolution.
!>        Its solid-harmonic moments retain every radial power of the native contraction.
!>        The short-range field subtracts this convolution from the original Coulomb field.
!> \param fields Prepared shell-set cache, replaced by the selected interaction.
!> \param omega Nonnegative inverse screening length in bohr^-1.
!> \param short_range Select erfc instead of erf.
! **************************************************************************************************
   SUBROUTINE gaussian_coulomb_fields_screen(fields, omega, short_range)
      TYPE(gaussian_coulomb_fields_type), INTENT(INOUT)  :: fields
      REAL(dp), INTENT(IN)                               :: omega
      LOGICAL, INTENT(IN)                                :: short_range

      INTEGER                                            :: block, count, j, k, l, power
      REAL(dp)                                           :: c, factor, inverse_width, ratio, sign
      TYPE(gaussian_coulomb_primitive_type), ALLOCATABLE :: screened(:)

      CPASSERT(ALLOCATED(fields%primitives))
      CPASSERT(ieee_is_finite(omega) .AND. omega >= 0.0_dp)
      IF (omega == 0.0_dp) THEN
         IF (.NOT. short_range) THEN
            DO block = 1, SIZE(fields%primitives)
               fields%primitives(block)%coefficients(:, :) = 0.0_dp
            END DO
         END IF
         RETURN
      END IF
      count = SUM(fields%primitives%power + 1)
      IF (short_range) count = count + SIZE(fields%primitives)
      ALLOCATE (screened(count))
      count = 0
      IF (short_range) THEN
         count = SIZE(fields%primitives)
         screened(:count) = fields%primitives
      END IF
      sign = 1.0_dp
      IF (short_range) sign = -1.0_dp
      DO block = 1, SIZE(fields%primitives)
         ASSOCIATE (primitive => fields%primitives(block))
            l = primitive%l
            power = primitive%power
            ! Avoid squaring a very large screening parameter.
            IF (omega >= SQRT(primitive%alpha)) THEN
               ratio = SQRT(primitive%alpha)/omega
               c = 1.0_dp/(1.0_dp + ratio**2)
               inverse_width = ratio**2*c/primitive%alpha
            ELSE
               ratio = omega/SQRT(primitive%alpha)
               c = ratio**2/(1.0_dp + ratio**2)
               inverse_width = 1.0_dp/(primitive%alpha*(1.0_dp + ratio**2))
            END IF
            DO j = power, 0, -1
               factor = c**(REAL(l, dp) + 1.5_dp + 2.0_dp*REAL(j, dp))
               DO k = j + 1, power
                  factor = factor*REAL(k, dp)*REAL(2*l + 2*k + 1, dp)*inverse_width/ &
                           (2.0_dp*REAL(k - j, dp))
               END DO
               count = count + 1
               screened(count)%alpha = primitive%alpha*c
               screened(count)%l = l
               screened(count)%power = j
               ALLOCATE (screened(count)%coefficients, MOLD=primitive%coefficients)
               screened(count)%coefficients(:, :) = sign*factor*primitive%coefficients
            END DO
         END ASSOCIATE
      END DO
      CALL MOVE_ALLOC(screened, fields%primitives)
   END SUBROUTINE gaussian_coulomb_fields_screen

! **************************************************************************************************
!> \brief Evaluate full native Gaussian Coulomb potentials and their evaluation-point gradients.
!>        A Gaussian-center displacement has the opposite derivative. Image translations and
!>        Bloch phases belong to the caller. The cache and shared angular tables are read only.
!> \param fields Prepared native shell-set cache.
!> \param displacement Evaluation point minus Gaussian center in bohr.
!> \param potential New potential for each native contracted function.
!> \param gradient Optional new Cartesian point gradient, axis first.
! **************************************************************************************************
   SUBROUTINE gaussian_coulomb_fields_evaluate(fields, displacement, potential, gradient)
      TYPE(gaussian_coulomb_fields_type), INTENT(IN)     :: fields
      REAL(dp), INTENT(IN)                               :: displacement(3)
      REAL(dp), INTENT(OUT)                              :: potential(:)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: gradient(:, :)

      INTEGER                                            :: a, axis, block, l, m, powers(3), &
                                                            reduced(3)
      REAL(dp)                                           :: angular, dangular(3), dh, factor, h

      CPASSERT(ALLOCATED(fields%primitives))
      CPASSERT(SIZE(potential) == fields%nfunction)
      CPASSERT(ALL(ieee_is_finite(displacement)))
      potential = 0.0_dp
      IF (PRESENT(gradient)) THEN
         CPASSERT(ALL(SHAPE(gradient) == [3, fields%nfunction]))
         gradient = 0.0_dp
      END IF
      DO block = 1, SIZE(fields%primitives)
         ASSOCIATE (primitive => fields%primitives(block))
            IF (.NOT. ANY(primitive%coefficients /= 0.0_dp)) CYCLE
            l = primitive%l
            CALL radial_potential(primitive%alpha, primitive%power, l, SUM(displacement**2), h, dh)
            DO m = 1, 2*l + 1
               angular = 0.0_dp
               dangular = 0.0_dp
               DO a = 1, nco(l)
                  powers = indco(:, ncoset(l - 1) + a)
                  factor = orbtramat(l)%slm(m, a)
                  angular = angular + factor*PRODUCT(displacement**powers)
                  IF (.NOT. PRESENT(gradient)) CYCLE
                  DO axis = 1, 3
                     IF (powers(axis) == 0) CYCLE
                     reduced = powers
                     reduced(axis) = reduced(axis) - 1
                     dangular(axis) = dangular(axis) + factor*REAL(powers(axis), dp)*PRODUCT(displacement**reduced)
                  END DO
               END DO
               potential = potential + h*angular*primitive%coefficients(m, :)
               IF (PRESENT(gradient)) THEN
                  DO axis = 1, 3
                     gradient(axis, :) = gradient(axis, :) + &
                                         (2.0_dp*displacement(axis)*dh*angular + h*dangular(axis))*primitive%coefficients(m, :)
                  END DO
               END IF
            END DO
         END ASSOCIATE
      END DO
   END SUBROUTINE gaussian_coulomb_fields_evaluate

! **************************************************************************************************
!> \brief Sum auxiliary potentials over prescribed integer cells without changing their kernels.
!>        The caller supplies the integration domain. A finite Coulomb image sum is not an
!>        infinite periodic Coulomb potential. Cell derivatives hold the integer images and
!>        Cartesian displacement fixed and include the motion of each translated auxiliary.
!> \param fields ...
!> \param displacement Point relative to the untranslated auxiliary center.
!> \param cell Physical cell matrix with lattice vectors in columns.
!> \param images Integer auxiliary cells, indexed (axis,image).
!> \param potential ...
!> \param gradient Optional Cartesian displacement gradient, indexed (axis,function).
!> \param cell_gradient Optional derivative with respect to cell(axis,column), indexed (axis,column,function).
! **************************************************************************************************
   SUBROUTINE gaussian_coulomb_fields_periodic(fields, displacement, cell, images, potential, gradient, cell_gradient)
      TYPE(gaussian_coulomb_fields_type), INTENT(IN)     :: fields
      REAL(dp), INTENT(IN)                               :: displacement(3), cell(3, 3)
      INTEGER, INTENT(IN)                                :: images(:, :)
      REAL(dp), INTENT(OUT)                              :: potential(:)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: gradient(:, :), cell_gradient(:, :, :)

      INTEGER                                            :: column, image
      LOGICAL                                            :: derivatives
      REAL(dp) :: point(3), value(fields%nfunction), derivative(3, fields%nfunction)

      CPASSERT(ALLOCATED(fields%primitives))
      CPASSERT(ALL(ieee_is_finite(displacement)))
      CPASSERT(ALL(ieee_is_finite(cell)))
      CPASSERT(SIZE(images, 1) == 3 .AND. SIZE(potential) == fields%nfunction)
      potential = 0.0_dp
      IF (PRESENT(gradient)) THEN
         CPASSERT(ALL(SHAPE(gradient) == [3, fields%nfunction]))
         gradient = 0.0_dp
      END IF
      IF (PRESENT(cell_gradient)) THEN
         CPASSERT(ALL(SHAPE(cell_gradient) == [3, 3, fields%nfunction]))
         cell_gradient = 0.0_dp
      END IF
      derivatives = PRESENT(gradient) .OR. PRESENT(cell_gradient)
      DO image = 1, SIZE(images, 2)
         point = displacement - MATMUL(cell, REAL(images(:, image), dp))
         IF (derivatives) THEN
            CALL gaussian_coulomb_fields_evaluate(fields, point, value, derivative)
         ELSE
            CALL gaussian_coulomb_fields_evaluate(fields, point, value)
         END IF
         potential(:) = potential + value
         IF (PRESENT(gradient)) gradient(:, :) = gradient + derivative
         IF (PRESENT(cell_gradient)) THEN
            DO column = 1, 3
               cell_gradient(:, column, :) = cell_gradient(:, column, :) - derivative*REAL(images(column, image), dp)
            END DO
         END IF
      END DO
   END SUBROUTINE gaussian_coulomb_fields_periodic

! **************************************************************************************************
!> \brief Full radial coefficient of r^(2*power+l)Y_lm exp(-alpha*r^2), regular at the origin.
!> \param alpha ...
!> \param power ...
!> \param l ...
!> \param s Squared displacement.
!> \param h Radial coefficient multiplying the solid harmonic.
!> \param dh Derivative with respect to squared displacement.
! **************************************************************************************************
   SUBROUTINE radial_potential(alpha, power, l, s, h, dh)
      REAL(dp), INTENT(IN)                               :: alpha
      INTEGER, INTENT(IN)                                :: power, l
      REAL(dp), INTENT(IN)                               :: s
      REAL(dp), INTENT(OUT)                              :: h, dh

      INTEGER                                            :: j, n
      REAL(dp)                                           :: ex, f, outer, previous, term, x

      x = alpha*s
      ex = EXP(-x)
      outer = 0.5_dp*ex/alpha
      DO j = 1, power
         outer = (REAL(j, dp)*outer + 0.5_dp*s**j*ex)/alpha
      END DO
      n = l + power + 1
      IF (x >= REAL(n + 1, dp)) THEN
         f = 0.5_dp*SQRT(pi/x)*ERF(SQRT(x))
         DO j = 1, n
            f = (REAL(2*j - 1, dp)*f - ex)/(2.0_dp*x)
         END DO
      ELSE
         ! F_n(x) = exp(-x) sum_j (2*x)^j / product_{k=0}^j (2*n+2*k+1).
         ! Positive terms avoid cancellation and only the required moment is retained.
         ! In this range the term ratio is below one, without a pointwise heap array.
         term = 1.0_dp/REAL(2*n + 1, dp)
         f = term
         j = 0
         DO
            j = j + 1
            term = term*(2.0_dp*x)/REAL(2*n + 2*j + 1, dp)
            previous = f
            f = f + term
            IF (f == previous) EXIT
         END DO
         f = ex*f
      END IF
      h = fourpi/REAL(2*l + 1, dp)*(s**(power + 1)*f + outer)
      dh = -2.0_dp*pi*s**power*f
   END SUBROUTINE radial_potential

! **************************************************************************************************
!> \brief Projection of a Cartesian monomial onto a normalized real solid harmonic.
!> \param powers ...
!> \param l ...
!> \param m Native magnetic index, 1..2*l+1.
!> \return Coefficient multiplying r^(SUM(powers)-l) times the solid harmonic.
! **************************************************************************************************
   FUNCTION monomial_harmonic(powers, l, m) RESULT(coefficient)
      INTEGER, INTENT(IN)                                :: powers(3), l, m
      REAL(dp)                                           :: coefficient

      INTEGER                                            :: a, combined(3)

      coefficient = 0.0_dp
      DO a = 1, nco(l)
         combined = powers + indco(:, ncoset(l - 1) + a)
         IF (ANY(MODULO(combined, 2) /= 0)) CYCLE
         coefficient = coefficient + orbtramat(l)%slm(m, a)*fourpi* &
                       PRODUCT(dfac(combined - 1))/dfac(SUM(combined) + 1)
      END DO
   END FUNCTION monomial_harmonic

END MODULE gaussian_coulomb_fields
