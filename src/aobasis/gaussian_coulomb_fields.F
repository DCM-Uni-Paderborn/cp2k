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
   USE gamma,                           ONLY: fgamma_0
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: dfac,&
                                              fac,&
                                              fourpi,&
                                              pi
   USE orbital_pointers,                ONLY: coset,&
                                              indco,&
                                              init_orbital_pointers,&
                                              nco,&
                                              ncoset
   USE orbital_transformation_matrices, ONLY: init_spherical_harmonics,&
                                              orbtramat
   USE t_c_g0,                          ONLY: get_lmax_init,&
                                              t_c_g0_n,&
                                              t_c_g0_response
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   TYPE gaussian_coulomb_primitive_type
      REAL(dp) :: alpha = 0.0_dp
      INTEGER :: l = 0, power = 0
      REAL(dp), ALLOCATABLE :: coefficients(:, :), hermite(:, :)
   END TYPE gaussian_coulomb_primitive_type

   TYPE gaussian_coulomb_fields_type
      INTEGER :: nfunction = 0
      LOGICAL :: identity = .FALSE.
      REAL(dp) :: cutoff_radius = 0.0_dp
      TYPE(gaussian_coulomb_primitive_type), ALLOCATABLE :: primitives(:)
   END TYPE gaussian_coulomb_fields_type

   PUBLIC :: gaussian_coulomb_fields_type, gaussian_coulomb_fields_init
   PUBLIC :: gaussian_coulomb_fields_prepare, gaussian_coulomb_fields_evaluate
   PUBLIC :: gaussian_coulomb_fields_evaluate_batch
   PUBLIC :: gaussian_coulomb_fields_screen
   PUBLIC :: gaussian_coulomb_fields_truncate
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
!> \param identity Optional overlap interaction, returning the Gaussian functions themselves.
! **************************************************************************************************
   SUBROUTINE gaussian_coulomb_fields_prepare(basis, iset, fields, identity)
      TYPE(gto_basis_set_type), INTENT(IN)               :: basis
      INTEGER, INTENT(IN)                                :: iset
      TYPE(gaussian_coulomb_fields_type), INTENT(OUT)    :: fields
      LOGICAL, INTENT(IN), OPTIONAL                      :: identity

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
      IF (PRESENT(identity)) fields%identity = identity
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
      CPASSERT(.NOT. fields%identity)
      CPASSERT(fields%cutoff_radius == 0.0_dp)
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
!> \brief Prepare native Gaussian charges for the truncated Coulomb interaction.
!>        Hermite coefficients retain the complete Cartesian contraction. Values and point
!>        gradients use the native moment interpolants and their actual T derivatives.
!> \param fields Owned original Gaussian charge cache.
!> \param cutoff_radius Fixed physical interaction radius in bohr.
! **************************************************************************************************
   SUBROUTINE gaussian_coulomb_fields_truncate(fields, cutoff_radius)
      TYPE(gaussian_coulomb_fields_type), INTENT(INOUT)  :: fields
      REAL(dp), INTENT(IN)                               :: cutoff_radius

      INTEGER                                            :: block

      CPASSERT(ALLOCATED(fields%primitives) .AND. .NOT. fields%identity)
      CPASSERT(ieee_is_finite(cutoff_radius) .AND. cutoff_radius > 0.0_dp)
      fields%cutoff_radius = cutoff_radius
      DO block = 1, SIZE(fields%primitives)
         CALL prepare_hermite_charge(fields%primitives(block))
      END DO
   END SUBROUTINE gaussian_coulomb_fields_truncate

! **************************************************************************************************
!> \brief Expand a Gaussian solid harmonic times a radial power in Cartesian Hermite charges.
!> \param primitive ...
! **************************************************************************************************
   SUBROUTINE prepare_hermite_charge(primitive)
      TYPE(gaussian_coulomb_primitive_type), &
         INTENT(INOUT)                                   :: primitive

      INTEGER                                            :: a, ax, ay, az, degree, m, order, &
                                                            powers(3), row, tx, ty, tz
      REAL(dp) :: coefficient, &
         monomials(ncoset(primitive%l + 2*primitive%power), SIZE(primitive%coefficients, 2)), &
         expansion(0:primitive%l + 2*primitive%power, 0:primitive%l + 2*primitive%power)

      degree = primitive%l + 2*primitive%power
      expansion = 0.0_dp
      expansion(0, 0) = 1.0_dp
      DO order = 0, degree - 1
         DO tx = 0, order + 1
            IF (tx > 0) expansion(order + 1, tx) = expansion(order, tx - 1)/(2.0_dp*primitive%alpha)
            IF (tx < order) expansion(order + 1, tx) = expansion(order + 1, tx) + (tx + 1)*expansion(order, tx + 1)
         END DO
      END DO
      monomials = 0.0_dp
      DO a = 1, nco(primitive%l)
         DO ax = 0, primitive%power
            DO ay = 0, primitive%power - ax
               az = primitive%power - ax - ay
               powers = indco(:, ncoset(primitive%l - 1) + a) + 2*[ax, ay, az]
               row = coset(powers(1), powers(2), powers(3))
               coefficient = fac(primitive%power)/(fac(ax)*fac(ay)*fac(az))
               DO m = 1, 2*primitive%l + 1
                  monomials(row, :) = monomials(row, :) + coefficient*orbtramat(primitive%l)%slm(m, a)* &
                                      primitive%coefficients(m, :)
               END DO
            END DO
         END DO
      END DO
      IF (ALLOCATED(primitive%hermite)) DEALLOCATE (primitive%hermite)
      ALLOCATE (primitive%hermite(ncoset(degree), SIZE(primitive%coefficients, 2)), SOURCE=0.0_dp)
      DO a = 1, ncoset(degree)
         IF (.NOT. ANY(monomials(a, :) /= 0.0_dp)) CYCLE
         powers = indco(:, a)
         DO tx = 0, powers(1)
         DO ty = 0, powers(2)
         DO tz = 0, powers(3)
            coefficient = REAL((-1)**(tx + ty + tz), dp)*expansion(powers(1), tx)* &
                          expansion(powers(2), ty)*expansion(powers(3), tz)
            row = coset(tx, ty, tz)
            primitive%hermite(row, :) = primitive%hermite(row, :) + coefficient*monomials(a, :)
         END DO
         END DO
         END DO
      END DO
   END SUBROUTINE prepare_hermite_charge

! **************************************************************************************************
!> \brief Contract the native truncated moments with one complete Hermite Gaussian charge.
!>        Differentiating the actual recurrence retains dF_n/dT independently of -F_(n+1).
!> \param primitive ...
!> \param cutoff_radius ...
!> \param displacement ...
!> \param potential Accumulated contracted potentials.
!> \param gradient Optional accumulated Cartesian point gradients.
! **************************************************************************************************
   SUBROUTINE truncated_primitive(primitive, cutoff_radius, displacement, potential, gradient)
      TYPE(gaussian_coulomb_primitive_type), INTENT(IN)  :: primitive
      REAL(dp), INTENT(IN)                               :: cutoff_radius, displacement(3)
      REAL(dp), INTENT(INOUT)                            :: potential(:)
      REAL(dp), INTENT(INOUT), OPTIONAL                  :: gradient(:, :)

      INTEGER                                            :: axis, degree, first_axis, last, n, &
                                                            previous, reduced(3), row, second
      LOGICAL                                            :: use_gamma
      REAL(dp) :: alpha, factor, r, t, &
         values(ncoset(primitive%l + 2*primitive%power), 0:primitive%l + 2*primitive%power), &
         response(0:primitive%l + 2*primitive%power), &
         moments(0:primitive%l + 2*primitive%power + 1), &
         derivatives(3, ncoset(primitive%l + 2*primitive%power), 0:primitive%l + 2*primitive%power)

      last = primitive%l + 2*primitive%power
      CPASSERT(ALLOCATED(primitive%hermite))
      CPASSERT(get_lmax_init() >= last)
      alpha = primitive%alpha
      r = cutoff_radius*SQRT(alpha)
      t = alpha*SUM(displacement**2)
      IF (PRESENT(gradient)) THEN
         CALL t_c_g0_response(moments, response, use_gamma, r, t, last)
      ELSE
         CALL t_c_g0_n(moments, use_gamma, r, t, last)
      END IF
      IF (use_gamma) THEN
         CALL fgamma_0(last + 1, t, moments)
         IF (PRESENT(gradient)) response = -moments(1:last + 1)
      END IF
      values = 0.0_dp
      IF (PRESENT(gradient)) derivatives = 0.0_dp
      factor = 1.0_dp
      DO n = 0, last
         values(1, n) = factor*moments(n)
         IF (PRESENT(gradient)) derivatives(:, 1, n) = (factor*2.0_dp*alpha*response(n))*displacement
         factor = -2.0_dp*alpha*factor
      END DO
      DO degree = 1, last
         DO row = ncoset(degree - 1) + 1, ncoset(degree)
            reduced = indco(:, row)
            first_axis = FINDLOC(reduced > 0, .TRUE., DIM=1)
            reduced(first_axis) = reduced(first_axis) - 1
            previous = coset(reduced(1), reduced(2), reduced(3))
            second = 0
            IF (reduced(first_axis) > 0) THEN
               reduced(first_axis) = reduced(first_axis) - 1
               second = coset(reduced(1), reduced(2), reduced(3))
               reduced(first_axis) = reduced(first_axis) + 1
            END IF
            DO n = 0, last - degree
               values(row, n) = displacement(first_axis)*values(previous, n + 1)
               IF (second > 0) values(row, n) = values(row, n) + reduced(first_axis)*values(second, n + 1)
               IF (.NOT. PRESENT(gradient)) CYCLE
               derivatives(:, row, n) = displacement(first_axis)*derivatives(:, previous, n + 1)
               derivatives(first_axis, row, n) = derivatives(first_axis, row, n) + values(previous, n + 1)
               IF (second > 0) derivatives(:, row, n) = derivatives(:, row, n) + &
                                                        reduced(first_axis)*derivatives(:, second, n + 1)
            END DO
         END DO
      END DO
      factor = 2.0_dp*pi/alpha
      potential = potential + factor*MATMUL(values(:, 0), primitive%hermite)
      IF (PRESENT(gradient)) THEN
         DO axis = 1, 3
            gradient(axis, :) = gradient(axis, :) + factor*MATMUL(derivatives(axis, :, 0), primitive%hermite)
         END DO
      END IF
   END SUBROUTINE truncated_primitive

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

      INTEGER, PARAMETER                                 :: max_cached_l = 4

      INTEGER                                            :: axis, block, l, m
      LOGICAL                                            :: angular_ready(0:max_cached_l)
      REAL(dp) :: angular, cached_angular((max_cached_l + 1)**2), &
         cached_gradient(3, (max_cached_l + 1)**2), dangular(3), dh, exponential, h, last_alpha, s

      CPASSERT(ALLOCATED(fields%primitives))
      CPASSERT(SIZE(potential) == fields%nfunction)
      CPASSERT(ALL(ieee_is_finite(displacement)))
      potential = 0.0_dp
      IF (PRESENT(gradient)) THEN
         CPASSERT(ALL(SHAPE(gradient) == [3, fields%nfunction]))
         gradient = 0.0_dp
      END IF
      IF (fields%cutoff_radius > 0.0_dp) THEN
         DO block = 1, SIZE(fields%primitives)
            IF (.NOT. ANY(fields%primitives(block)%coefficients /= 0.0_dp)) CYCLE
            CALL truncated_primitive(fields%primitives(block), fields%cutoff_radius, displacement, potential, gradient)
         END DO
         RETURN
      END IF
      s = SUM(displacement**2)
      angular_ready = .FALSE.
      last_alpha = -1.0_dp
      DO block = 1, SIZE(fields%primitives)
         ASSOCIATE (primitive => fields%primitives(block))
            IF (.NOT. ANY(primitive%coefficients /= 0.0_dp)) CYCLE
            l = primitive%l
            IF (fields%identity) THEN
               ! Angular components of one primitive share the same exponential.
               IF (primitive%alpha /= last_alpha) THEN
                  exponential = EXP(-primitive%alpha*s)
                  last_alpha = primitive%alpha
               END IF
               ! A vanished exponential contributes exact zeros for this radial power.
               IF (exponential == 0.0_dp .AND. primitive%power == 0) CYCLE
               h = s**primitive%power*exponential
               dh = -primitive%alpha*h
               IF (primitive%power > 0) dh = dh + primitive%power*s**(primitive%power - 1)*exponential
            ELSE
               CALL radial_potential(primitive%alpha, primitive%power, l, s, h, dh)
            END IF
            ! Solid harmonics do not depend on the radial exponent or power. Reuse them through
            ! g in a fixed workspace, retaining the generic evaluation for higher angular orders.
            IF (l <= max_cached_l) THEN
               IF (.NOT. angular_ready(l)) THEN
                  DO m = 1, 2*l + 1
                     CALL native_solid_harmonic(l, m, displacement, PRESENT(gradient), &
                                                cached_angular(l*l + m), cached_gradient(:, l*l + m))
                  END DO
                  angular_ready(l) = .TRUE.
               END IF
            END IF
            DO m = 1, 2*l + 1
               IF (l <= max_cached_l) THEN
                  angular = cached_angular(l*l + m)
                  dangular = cached_gradient(:, l*l + m)
               ELSE
                  CALL native_solid_harmonic(l, m, displacement, PRESENT(gradient), angular, dangular)
               END IF
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
!> \brief Evaluate a point batch with the original native Gaussian potential and gradient conventions.
!>        Overlap and Coulomb fields share angular polynomials across the point batch.
!>        Truncated kernels and higher angular orders retain the scalar evaluation.
!> \param fields Prepared native shell-set cache.
!> \param displacement Evaluation points minus the Gaussian center, indexed (axis,point).
!> \param potential New potentials, indexed (point,function).
!> \param gradient Optional Cartesian point gradients, indexed (axis,function,point).
! **************************************************************************************************
   SUBROUTINE gaussian_coulomb_fields_evaluate_batch(fields, displacement, potential, gradient)
      TYPE(gaussian_coulomb_fields_type), INTENT(IN)     :: fields
      REAL(dp), INTENT(IN)                               :: displacement(:, :)
      REAL(dp), INTENT(OUT)                              :: potential(:, :)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: gradient(:, :, :)

      INTEGER, PARAMETER                                 :: max_cached_l = 4

      INTEGER                                            :: angular_count, axis, block, column, l, &
                                                            m, npoint, point, slot
      LOGICAL                                            :: angular_ready(0:max_cached_l)
      REAL(dp)                                           :: last_alpha
      REAL(dp), ALLOCATABLE                              :: angular(:, :), &
                                                            angular_gradient(:, :, :), dh(:), &
                                                            exponential(:), h(:), squared_radius(:)

      CPASSERT(ALLOCATED(fields%primitives))
      npoint = SIZE(displacement, 2)
      CPASSERT(SIZE(displacement, 1) == 3)
      CPASSERT(ALL(SHAPE(potential) == [npoint, fields%nfunction]))
      CPASSERT(ALL(ieee_is_finite(displacement)))
      IF (PRESENT(gradient)) THEN
         CPASSERT(ALL(SHAPE(gradient) == [3, fields%nfunction, npoint]))
      END IF
      IF (fields%cutoff_radius > 0.0_dp .OR. &
          ANY(fields%primitives%l > max_cached_l)) THEN
         DO point = 1, npoint
            IF (PRESENT(gradient)) THEN
               CALL gaussian_coulomb_fields_evaluate(fields, displacement(:, point), potential(point, :), gradient(:, :, point))
            ELSE
               CALL gaussian_coulomb_fields_evaluate(fields, displacement(:, point), potential(point, :))
            END IF
         END DO
         RETURN
      END IF
      angular_count = (MAX(0, MAXVAL(fields%primitives%l)) + 1)**2
      ALLOCATE (angular(npoint, angular_count), angular_gradient(3, npoint, angular_count), &
                squared_radius(npoint), exponential(npoint), h(npoint), dh(npoint))
      potential = 0.0_dp
      IF (PRESENT(gradient)) gradient = 0.0_dp
      squared_radius(:) = SUM(displacement**2, DIM=1)
      angular_ready = .FALSE.
      last_alpha = -1.0_dp
      DO block = 1, SIZE(fields%primitives)
         ASSOCIATE (primitive => fields%primitives(block))
            IF (.NOT. ANY(primitive%coefficients /= 0.0_dp)) CYCLE
            l = primitive%l
            IF (fields%identity) THEN
               IF (primitive%alpha /= last_alpha) THEN
                  exponential(:) = EXP(-primitive%alpha*squared_radius)
                  last_alpha = primitive%alpha
               END IF
               IF (primitive%power == 0) THEN
                  ! Match the scalar evaluator when the complete exponential batch vanishes.
                  IF (ALL(exponential == 0.0_dp)) CYCLE
                  h(:) = exponential
               ELSE
                  h(:) = squared_radius**primitive%power*exponential
               END IF
               IF (PRESENT(gradient)) THEN
                  dh(:) = -primitive%alpha*h
                  IF (primitive%power > 0) dh(:) = dh + primitive%power*squared_radius**(primitive%power - 1)*exponential
               END IF
            ELSE
               DO point = 1, npoint
                  CALL radial_potential(primitive%alpha, primitive%power, l, squared_radius(point), h(point), dh(point))
               END DO
            END IF
            IF (.NOT. angular_ready(l)) THEN
               CALL native_solid_harmonics_batch(l, displacement, PRESENT(gradient), &
                                                 angular(:, l*l + 1:(l + 1)**2), &
                                                 angular_gradient(:, :, l*l + 1:(l + 1)**2))
               angular_ready(l) = .TRUE.
            END IF
            DO m = 1, 2*l + 1
               slot = l*l + m
               DO column = 1, fields%nfunction
                  IF (primitive%coefficients(m, column) == 0.0_dp) CYCLE
                  potential(:, column) = potential(:, column) + h*angular(:, slot)*primitive%coefficients(m, column)
                  IF (.NOT. PRESENT(gradient)) CYCLE
                  DO axis = 1, 3
                     gradient(axis, column, :) = gradient(axis, column, :) + &
                                           (2.0_dp*displacement(axis, :)*dh*angular(:, slot) + h*angular_gradient(axis, :, slot))* &
                                                 primitive%coefficients(m, column)
                  END DO
               END DO
            END DO
         END ASSOCIATE
      END DO
   END SUBROUTINE gaussian_coulomb_fields_evaluate_batch

! **************************************************************************************************
!> \brief Share Cartesian powers and monomials across a batch of native solid harmonics.
!>        The Cartesian summation order and normalization match the scalar evaluation.
!> \param l Angular momentum.
!> \param displacement Point-minus-center coordinates, indexed (axis,point).
!> \param derivatives Whether point gradients are required.
!> \param angular Solid harmonics, indexed (point,magnetic index).
!> \param gradient Point gradients, indexed (axis,point,magnetic index).
! **************************************************************************************************
   SUBROUTINE native_solid_harmonics_batch(l, displacement, derivatives, angular, gradient)
      INTEGER, INTENT(IN)                                :: l
      REAL(dp), INTENT(IN)                               :: displacement(:, :)
      LOGICAL, INTENT(IN)                                :: derivatives
      REAL(dp), INTENT(OUT)                              :: angular(:, :), gradient(:, :, :)

      INTEGER                                            :: a, axis, degree, m, powers(3), reduced(3)
      REAL(dp) :: factor, monomial_gradient(SIZE(displacement, 2), 3), &
         monomial(SIZE(displacement, 2)), coordinate_power(SIZE(displacement, 2), 3, 0:l)

      coordinate_power(:, :, 0) = 1.0_dp
      DO degree = 1, l
         DO axis = 1, 3
            coordinate_power(:, axis, degree) = coordinate_power(:, axis, degree - 1)*displacement(axis, :)
         END DO
      END DO
      angular = 0.0_dp
      gradient = 0.0_dp
      DO a = 1, nco(l)
         powers = indco(:, ncoset(l - 1) + a)
         monomial = coordinate_power(:, 1, powers(1))*coordinate_power(:, 2, powers(2))* &
                    coordinate_power(:, 3, powers(3))
         IF (derivatives) THEN
            monomial_gradient = 0.0_dp
            DO axis = 1, 3
               IF (powers(axis) == 0) CYCLE
               reduced = powers
               reduced(axis) = reduced(axis) - 1
               monomial_gradient(:, axis) = REAL(powers(axis), dp)*coordinate_power(:, 1, reduced(1))* &
                                            coordinate_power(:, 2, reduced(2))*coordinate_power(:, 3, reduced(3))
            END DO
         END IF
         DO m = 1, 2*l + 1
            factor = orbtramat(l)%slm(m, a)
            angular(:, m) = angular(:, m) + factor*monomial
            IF (.NOT. derivatives) CYCLE
            DO axis = 1, 3
               gradient(axis, :, m) = gradient(axis, :, m) + factor*monomial_gradient(:, axis)
            END DO
         END DO
      END DO
   END SUBROUTINE native_solid_harmonics_batch

! **************************************************************************************************
!> \brief Native solid harmonic and Cartesian point gradient in the original polynomial order.
!> \param l ...
!> \param m Native magnetic index.
!> \param displacement ...
!> \param derivatives Whether the point gradient is required.
!> \param angular ...
!> \param gradient ...
! **************************************************************************************************
   SUBROUTINE native_solid_harmonic(l, m, displacement, derivatives, angular, gradient)
      INTEGER, INTENT(IN)                                :: l, m
      REAL(dp), INTENT(IN)                               :: displacement(3)
      LOGICAL, INTENT(IN)                                :: derivatives
      REAL(dp), INTENT(OUT)                              :: angular, gradient(3)

      INTEGER                                            :: a, axis, powers(3), reduced(3)
      REAL(dp)                                           :: factor

      angular = 0.0_dp
      gradient = 0.0_dp
      DO a = 1, nco(l)
         powers = indco(:, ncoset(l - 1) + a)
         factor = orbtramat(l)%slm(m, a)
         angular = angular + factor*PRODUCT(displacement**powers)
         IF (.NOT. derivatives) CYCLE
         DO axis = 1, 3
            IF (powers(axis) == 0) CYCLE
            reduced = powers
            reduced(axis) = reduced(axis) - 1
            gradient(axis) = gradient(axis) + factor*REAL(powers(axis), dp)*PRODUCT(displacement**reduced)
         END DO
      END DO
   END SUBROUTINE native_solid_harmonic

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
         point = displacement - MATMUL(cell, REAL(images(1:3, image), dp))
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
