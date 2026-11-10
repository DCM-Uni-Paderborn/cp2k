!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Half-numerical overlaps with the original scalar UPF radial projectors.
!>        Cartesian polynomials and angular products are integrated analytically.
!>        Only a radial quadrature remains; no fit or rotation of projector identities
!>        is performed. The UPF reader preserves the tabulated beta amplitude and
!>        includes the Ry-to-Hartree factor in the coupling matrix D.
! **************************************************************************************************
MODULE upf_projector_integrals
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE ISO_C_BINDING,                   ONLY: C_DOUBLE,&
                                              C_INT
   USE kinds,                           ONLY: dp
   USE libgrpp,                         ONLY: libgrpp_init
   USE mathconstants,                   ONLY: dfac,&
                                              fourpi,&
                                              maxfac,&
                                              twopi
   USE mathlib,                         ONLY: binomial
   USE orbital_pointers,                ONLY: coset,&
                                              indco,&
                                              init_orbital_pointers,&
                                              nco,&
                                              ncoset
   USE orbital_transformation_matrices, ONLY: init_spherical_harmonics,&
                                              orbtramat
   USE upf_types,                       ONLY: atom_upfpot_type,&
                                              upf_radial_weights
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'upf_projector_integrals'

   TYPE upf_projector_type
      INTEGER :: degree = -1, lmax = -1, harmonic_max = -1, nproj = 0
      REAL(KIND=dp) :: radius = 0.0_dp
      INTEGER, ALLOCATABLE :: radial_index(:), angular_index(:), lbeta(:)
      REAL(KIND=dp), ALLOCATABLE :: r(:), weighted_beta(:, :), d(:, :), gaunt(:, :, :)
   END TYPE upf_projector_type

   PUBLIC :: upf_projector_type, upf_projector_init, upf_projector_overlap
   PUBLIC :: upf_projector_angular_init, upf_projector_product_gaunt
   PUBLIC :: upf_projector_gaussian_overlap

   INTERFACE
      FUNCTION scaled_bessel(n, x) RESULT(value) BIND(C, name="libgrpp_modified_bessel_scaled")
         IMPORT :: C_DOUBLE, C_INT
         INTEGER(C_INT), VALUE :: n
         REAL(C_DOUBLE), VALUE :: x
         REAL(C_DOUBLE) :: value
      END FUNCTION scaled_bessel
   END INTERFACE

CONTAINS

! **************************************************************************************************
!> \brief Contract scalar UPF projectors with native Gaussian AOs and explicit Bloch images.
!>        Only the original radial quadrature is numerical. Gaussian Cartesian integrals and
!>        AO-center derivatives use the existing UPF integrator, without a uniform-grid solve.
!> \param basis Native contracted ORB basis and unchanged primitive normalization.
!> \param hard Original or transformed projectors with sufficient derivative degree.
!> \param center Home projector center.
!> \param ao_center Home Gaussian center.
!> \param cell Cartesian lattice vectors.
!> \param images Fixed integer Gaussian images.
!> \param reduced_k Bloch vector, with phase exp(i*2*pi*k.image).
!> \param overlaps Projector bras against Gaussian AO columns, without occupations.
!> \param response Optional joint first response.
!> \param center_response Optional projector-center displacement.
!> \param ao_response Optional Gaussian-center displacement.
!> \param cell_response Optional lattice displacement.
!> \param reduced_k_response Optional reduced Bloch-vector displacement.
!> \param kinetic Optional projector overlaps with the Gaussian kinetic operator.
!> \param kinetic_response Optional joint first response of the kinetic overlaps.
!> \param gradient Optional gradient with respect to relative AO displacement, at fixed phase.
! **************************************************************************************************
   SUBROUTINE upf_projector_gaussian_overlap(basis, hard, center, ao_center, cell, images, reduced_k, overlaps, &
                                            response, center_response, ao_response, cell_response, reduced_k_response, &
                                            kinetic, kinetic_response, gradient)
      TYPE(gto_basis_set_type), INTENT(IN)               :: basis
      TYPE(upf_projector_type), INTENT(IN)               :: hard
      REAL(dp), INTENT(IN)                               :: center(3), ao_center(3), cell(3, 3)
      INTEGER, INTENT(IN)                                :: images(:, :)
      REAL(dp), INTENT(IN)                               :: reduced_k(3)
      COMPLEX(dp), INTENT(OUT)                           :: overlaps(:, :)
      COMPLEX(dp), INTENT(OUT), OPTIONAL                 :: response(:, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: center_response(3), ao_response(3), &
                                                            cell_response(3, 3), &
                                                            reduced_k_response(3)
      COMPLEX(dp), INTENT(OUT), OPTIONAL                 :: kinetic(:, :), kinetic_response(:, :), gradient(:, :, :)

      COMPLEX(dp)                                        :: phase, phase_response
      INTEGER                                            :: axis, derivative_index, first, image, &
                                                            image_nder, iset, last, nder, np, &
                                                            nshell, powers(3)
      REAL(dp)                                           :: direction(3), displacement(3), dk(3), &
                                                            shift(3)
      REAL(dp), ALLOCATABLE                              :: contracted(:, :), derivative(:, :), &
                                                            integrals(:, :, :)

      CPASSERT(SIZE(images, 1) == 3)
      CPASSERT(ALL(SHAPE(overlaps) == [hard%nproj, basis%nsgf]))
      nder = 0
      IF (PRESENT(response)) THEN
         CPASSERT(ALL(SHAPE(response) == SHAPE(overlaps)))
         response = CMPLX(0.0_dp, 0.0_dp, dp)
         nder = 1
      END IF
      IF (PRESENT(gradient)) THEN
         CPASSERT(ALL(SHAPE(gradient) == [hard%nproj, basis%nsgf, 3]))
         gradient = CMPLX(0.0_dp, 0.0_dp, dp)
         nder = MAX(nder, 1)
      END IF
      IF (PRESENT(kinetic)) THEN
         CPASSERT(ALL(SHAPE(kinetic) == SHAPE(overlaps)))
         kinetic = CMPLX(0.0_dp, 0.0_dp, dp)
         nder = MAX(nder, 2)
      END IF
      IF (PRESENT(kinetic_response)) THEN
         CPASSERT(PRESENT(kinetic))
         CPASSERT(ALL(SHAPE(kinetic_response) == SHAPE(overlaps)))
         kinetic_response = CMPLX(0.0_dp, 0.0_dp, dp)
         nder = 3
      END IF
      overlaps = CMPLX(0.0_dp, 0.0_dp, dp)
      dk = 0.0_dp
      IF (PRESENT(reduced_k_response)) dk = reduced_k_response
      DO image = 1, SIZE(images, 2)
         shift = MATMUL(cell, REAL(images(:, image), dp))
         displacement = center - ao_center - shift
         direction = 0.0_dp
         IF (PRESENT(ao_response)) direction = direction + ao_response
         IF (PRESENT(center_response)) direction = direction - center_response
         IF (PRESENT(cell_response)) direction = direction + MATMUL(cell_response, REAL(images(:, image), dp))
         phase = EXP(CMPLX(0.0_dp, twopi*DOT_PRODUCT(reduced_k, REAL(images(:, image), dp)), dp))
         phase_response = CMPLX(0.0_dp, twopi*DOT_PRODUCT(dk, REAL(images(:, image), dp)), dp)*phase
         ! A phase-only response needs the unperturbed integrals, including the kinetic Laplacian.
         image_nder = nder
         IF (ALL(direction == 0.0_dp)) THEN
            image_nder = 0
            IF (PRESENT(gradient)) image_nder = 1
            IF (PRESENT(kinetic)) image_nder = 2
         END IF
         DO iset = 1, basis%nset
            nshell = basis%nshell(iset)
            first = MINVAL(basis%first_sgf(:nshell, iset))
            last = MAXVAL(basis%last_sgf(:nshell, iset))
            np = basis%npgf(iset)*ncoset(basis%lmax(iset))
            ALLOCATE (integrals(np, hard%nproj, ncoset(image_nder)), contracted(last - first + 1, hard%nproj))
            CALL upf_projector_overlap(hard, basis%lmax(iset), basis%zet(:basis%npgf(iset), iset), &
                                       displacement, image_nder, integrals)
            contracted(:, :) = MATMUL(TRANSPOSE(basis%sphi(:np, first:last)), integrals(:, :, 1))
            overlaps(:, first:last) = overlaps(:, first:last) + phase*TRANSPOSE(contracted)
            IF (PRESENT(gradient)) THEN
               DO axis = 1, 3
                  gradient(:, first:last, axis) = gradient(:, first:last, axis) + phase* &
                     TRANSPOSE(MATMUL(TRANSPOSE(basis%sphi(:np, first:last)), integrals(:, :, axis + 1)))
               END DO
            END IF
            IF (PRESENT(response)) THEN
               response(:, first:last) = response(:, first:last) + phase_response*TRANSPOSE(contracted)
               ALLOCATE (derivative(last - first + 1, hard%nproj))
               derivative = 0.0_dp
               DO axis = 1, 3
                  IF (direction(axis) == 0.0_dp) CYCLE
                  derivative(:, :) = derivative + direction(axis)* &
                                     MATMUL(TRANSPOSE(basis%sphi(:np, first:last)), integrals(:, :, axis + 1))
               END DO
               response(:, first:last) = response(:, first:last) + phase*TRANSPOSE(derivative)
               DEALLOCATE (derivative)
            END IF
            IF (PRESENT(kinetic)) THEN
               contracted = 0.0_dp
               DO axis = 1, 3
                  powers = 0
                  powers(axis) = 2
                  derivative_index = coset(powers(1), powers(2), powers(3))
                  contracted(:, :) = contracted - 0.5_dp* &
                                     MATMUL(TRANSPOSE(basis%sphi(:np, first:last)), integrals(:, :, derivative_index))
               END DO
               kinetic(:, first:last) = kinetic(:, first:last) + phase*TRANSPOSE(contracted)
               IF (PRESENT(kinetic_response)) THEN
                  kinetic_response(:, first:last) = kinetic_response(:, first:last) + phase_response*TRANSPOSE(contracted)
                  ALLOCATE (derivative(last - first + 1, hard%nproj))
                  derivative = 0.0_dp
                  DO axis = 1, 3
                     DO derivative_index = 1, 3
                        IF (direction(derivative_index) == 0.0_dp) CYCLE
                        powers = 0
                        powers(axis) = 2
                        powers(derivative_index) = powers(derivative_index) + 1
                        derivative(:, :) = derivative - 0.5_dp*direction(derivative_index)* &
                                           MATMUL(TRANSPOSE(basis%sphi(:np, first:last)), &
                                                  integrals(:, :, coset(powers(1), powers(2), powers(3))))
                     END DO
                  END DO
                  kinetic_response(:, first:last) = kinetic_response(:, first:last) + phase*TRANSPOSE(derivative)
                  DEALLOCATE (derivative)
               END IF
            END IF
            DEALLOCATE (integrals, contracted)
         END DO
      END DO
   END SUBROUTINE upf_projector_gaussian_overlap

! **************************************************************************************************
!> \brief Prepare immutable radial and angular data outside parallel regions.
!> \param pot Original UPF data, including beta(r), not r*beta(r).
!> \param degree Largest Cartesian degree needed (AO degree + derivatives/moments).
!> \param projectors Prepared data; reinitializing also releases the previous allocation.
! **************************************************************************************************
   SUBROUTINE upf_projector_init(pot, degree, projectors)
      TYPE(atom_upfpot_type), INTENT(IN)                 :: pot
      INTEGER, INTENT(IN)                                :: degree
      TYPE(upf_projector_type), INTENT(OUT)              :: projectors

      INTEGER                                            :: i, j, n, p, q
      REAL(dp), ALLOCATABLE                              :: weights(:)

      CPASSERT(degree >= 0)
      IF (pot%has_so) THEN
         CPABORT("Spin-orbit UPF projectors are not supported")
      END IF
      CPASSERT(pot%number_of_proj > 0)
      CALL upf_projector_angular_init(pot%lbeta, degree, projectors)
      CALL libgrpp_init()

      n = SIZE(pot%r)
      CPASSERT(n >= 3)
      CPASSERT(SIZE(pot%rab) == n)
      CPASSERT(SIZE(pot%beta, 1) == n)
      CPASSERT(SIZE(pot%beta, 2) == pot%number_of_proj)
      CPASSERT(ALL(pot%r(2:) > pot%r(:n - 1)))
      CPASSERT(ALL(pot%rab > 0.0_dp))
      ALLOCATE (projectors%r(n))
      projectors%r(:) = pot%r
      ALLOCATE (weights(n))
      weights(:) = upf_radial_weights(pot)*pot%r**2
      ALLOCATE (projectors%weighted_beta(n, pot%number_of_proj))
      DO i = 1, pot%number_of_proj
         projectors%weighted_beta(:, i) = weights*pot%beta(:, i)
         DO j = n, 1, -1
            IF (pot%beta(j, i) /= 0.0_dp) THEN
               projectors%radius = MAX(projectors%radius, pot%r(j))
               EXIT
            END IF
         END DO
      END DO

      ALLOCATE (projectors%d(projectors%nproj, projectors%nproj))
      projectors%d = 0.0_dp
      DO p = 1, projectors%nproj
         DO q = 1, projectors%nproj
            i = projectors%radial_index(p)
            j = projectors%radial_index(q)
            IF (projectors%angular_index(p) == projectors%angular_index(q)) THEN
               projectors%d(p, q) = pot%dion(i, j)
            ELSE IF (pot%lbeta(i) /= pot%lbeta(j)) THEN
               CPASSERT(ABS(pot%dion(i, j)) < 1.e-12_dp)
            END IF
         END DO
      END DO

   END SUBROUTINE upf_projector_init

! **************************************************************************************************
!> \brief Prepare only the native real-harmonic mapping and angular product coefficients.
!>        This also supports occupied core orbitals without manufacturing radial projectors.
!> \param lbeta Radial-channel angular momenta, in the requested original channel order.
!> \param degree Largest Cartesian degree; at least max(lbeta) for harmonic products.
!> \param projectors Angular cache; radial and coupling arrays remain unallocated.
! **************************************************************************************************
   SUBROUTINE upf_projector_angular_init(lbeta, degree, projectors)
      INTEGER, INTENT(IN)                                :: lbeta(:), degree
      TYPE(upf_projector_type), INTENT(OUT)              :: projectors

      INTEGER                                            :: a, b, i, k, l, lambda, m, p, powers(3), s
      REAL(dp), ALLOCATABLE                              :: moments(:, :)

      CPASSERT(degree >= 0 .AND. SIZE(lbeta) > 0)
      CPASSERT(ALL(lbeta >= 0))
      projectors%degree = degree
      projectors%lmax = MAXVAL(lbeta)
      projectors%harmonic_max = degree + projectors%lmax
      ! Both the factorial and bundled Bessel tables cover this range.
      CPASSERT(projectors%harmonic_max <= MIN(20, maxfac/2))
      CALL init_orbital_pointers(projectors%harmonic_max)
      CALL init_spherical_harmonics(projectors%harmonic_max, -1)
      ALLOCATE (projectors%lbeta, SOURCE=lbeta)
      projectors%nproj = SUM(2*lbeta + 1)
      ALLOCATE (projectors%radial_index(projectors%nproj), projectors%angular_index(projectors%nproj))
      p = 0
      DO i = 1, SIZE(lbeta)
         l = lbeta(i)
         DO m = 1, 2*l + 1
            p = p + 1
            projectors%radial_index(p) = i
            projectors%angular_index(p) = l*l + m
         END DO
      END DO
      ! G(s,lm,lambda*nu) = integral Omega**s Y_lm(Omega) Y_lambda,nu(Omega).
      ALLOCATE (projectors%gaunt(ncoset(degree), (projectors%lmax + 1)**2, &
                                 (projectors%harmonic_max + 1)**2))
      projectors%gaunt = 0.0_dp
      DO s = 1, ncoset(degree)
         k = SUM(indco(:, s))
         DO l = 0, projectors%lmax
            IF (.NOT. ANY(lbeta == l)) CYCLE
            DO lambda = MAX(0, l - k), l + k
               IF (MODULO(k + l + lambda, 2) /= 0) CYCLE
               ALLOCATE (moments(nco(l), nco(lambda)))
               DO b = 1, nco(lambda)
                  DO a = 1, nco(l)
                     powers = indco(:, s) + indco(:, ncoset(l - 1) + a) + indco(:, ncoset(lambda - 1) + b)
                     moments(a, b) = 0.0_dp
                     IF (ALL(MODULO(powers, 2) == 0)) THEN
                        moments(a, b) = fourpi*PRODUCT(dfac(powers - 1))/dfac(SUM(powers) + 1)
                     END IF
                  END DO
               END DO
               projectors%gaunt(s, l*l + 1:(l + 1)**2, lambda*lambda + 1:(lambda + 1)**2) = &
                  MATMUL(MATMUL(orbtramat(l)%slm, moments), TRANSPOSE(orbtramat(lambda)%slm))
               DEALLOCATE (moments)
            END DO
         END DO
      END DO
   END SUBROUTINE upf_projector_angular_init

! **************************************************************************************************
!> \brief Angular product coefficients for two native real-harmonic projectors.
!> \param projectors Prepared projectors with Cartesian degree at least their largest l.
!> \param gaunt Integral of Y_p Y_q Y_LM, in magnetic-projector order.
! **************************************************************************************************
   SUBROUTINE upf_projector_product_gaunt(projectors, gaunt)
      TYPE(upf_projector_type), INTENT(IN)               :: projectors
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: gaunt(:, :, :)

      INTEGER                                            :: l, m, nlm, p, q, s

      CPASSERT(projectors%degree >= projectors%lmax)
      nlm = (2*projectors%lmax + 1)**2
      ALLOCATE (gaunt(projectors%nproj, projectors%nproj, nlm))
      gaunt = 0.0_dp
      DO q = 1, projectors%nproj
         l = projectors%lbeta(projectors%radial_index(q))
         m = projectors%angular_index(q) - l*l
         DO p = 1, projectors%nproj
            DO s = 1, nco(l)
               gaunt(p, q, :) = gaunt(p, q, :) + orbtramat(l)%slm(m, s)* &
                                projectors%gaunt(ncoset(l - 1) + s, projectors%angular_index(p), :nlm)
            END DO
         END DO
      END DO
   END SUBROUTINE upf_projector_product_gaunt

! **************************************************************************************************
!> \brief Evaluate <Cartesian Gaussian|beta_i Y_lm> and AO-center derivatives.
!> \param projectors Immutable prepared data.
!> \param la Largest AO Cartesian degree (all degrees 0..la are returned).
!> \param zeta Primitive Gaussian exponents.
!> \param rac Projector center minus AO center, as in ai_overlap.
!> \param nder Maximum derivative degree, or moment degree in moment mode.
!> \param integrals Primitive-major AO order, original radial/m order, coset derivative order.
!> \param moment_origin Optional AO center minus moment reference; requests moments
!>        instead of derivatives, in the same coset ordering including the overlap.
! **************************************************************************************************
   SUBROUTINE upf_projector_overlap(projectors, la, zeta, rac, nder, integrals, moment_origin)
      TYPE(upf_projector_type), INTENT(IN)               :: projectors
      INTEGER, INTENT(IN)                                :: la
      REAL(KIND=dp), INTENT(IN)                          :: zeta(:), rac(3)
      INTEGER, INTENT(IN)                                :: nder
      REAL(KIND=dp), INTENT(OUT)                         :: integrals(:, :, :)
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: moment_origin(3)

      INTEGER                                            :: a(3), axis, b(3), degree, deriv(3), i, &
                                                            id, ip, j, k, l, lambda, nb, nc, nh, &
                                                            nr, p, power, s, shifted(3), step, t
      REAL(KIND=dp)                                      :: alpha, direction(3), distance, factor
      REAL(KIND=dp), ALLOCATABLE :: angular(:, :, :), cartesian(:, :), centered(:, :), &
         coeff(:, :), envelope(:), harmonics(:), monomials(:), next_coeff(:), radial(:, :, :), &
         radial_power(:), radial_work(:), shift(:, :)

      degree = la + nder
      CPASSERT(la >= 0 .AND. nder >= 0)
      CPASSERT(degree <= projectors%degree)
      CPASSERT(ALL(zeta > 0.0_dp))
      CPASSERT(SIZE(integrals, 1) == SIZE(zeta)*ncoset(la))
      CPASSERT(SIZE(integrals, 2) == projectors%nproj)
      CPASSERT(SIZE(integrals, 3) == ncoset(nder))
      nc = ncoset(degree)
      nh = degree + projectors%lmax
      nr = SIZE(projectors%r)
      nb = SIZE(projectors%weighted_beta, 2)
      ALLOCATE (harmonics((nh + 1)**2), monomials(ncoset(nh)))
      distance = NORM2(rac)
      direction = [0.0_dp, 0.0_dp, 1.0_dp]
      IF (distance > 0.0_dp) direction = -rac/distance
      DO i = 1, ncoset(nh)
         monomials(i) = PRODUCT(direction**indco(:, i))
      END DO
      DO l = 0, nh
         harmonics(l*l + 1:(l + 1)**2) = &
            MATMUL(orbtramat(l)%slm, monomials(ncoset(l - 1) + 1:ncoset(l)))
      END DO
      ALLOCATE (angular(nc, (projectors%lmax + 1)**2, 0:nh))
      DO lambda = 0, nh
         DO i = 1, SIZE(angular, 2)
            angular(:, i, lambda) = MATMUL( &
                                    projectors%gaunt(:nc, i, lambda*lambda + 1:(lambda + 1)**2), &
                                    harmonics(lambda*lambda + 1:(lambda + 1)**2))
         END DO
      END DO
      ALLOCATE (shift(nc, nc))
      shift = 0.0_dp
      DO t = 1, nc
         a = indco(:, t)
         DO s = 1, nc
            b = indco(:, s)
            IF (ANY(b > a)) CYCLE
            factor = 1.0_dp
            DO axis = 1, 3
               factor = factor*binomial(a(axis), b(axis))*rac(axis)**(a(axis) - b(axis))
            END DO
            shift(t, s) = factor
         END DO
      END DO
      ALLOCATE (radial(0:degree, nb, 0:nh), radial_work(nr), radial_power(nr), envelope(nr))
      ALLOCATE (centered(nc, projectors%nproj), cartesian(nc, projectors%nproj))
      ALLOCATE (coeff(0:degree, 3), next_coeff(0:degree))
      integrals = 0.0_dp
      DO ip = 1, SIZE(zeta)
         alpha = zeta(ip)
         envelope(:) = EXP(-alpha*(projectors%r - distance)**2)
         ! The scaled Bessel factor cannot recover an envelope that underflowed to zero.
         IF (ALL(envelope == 0.0_dp)) CYCLE
         DO lambda = 0, nh
            radial_work = 0.0_dp
            DO i = 1, nr
               IF (envelope(i) == 0.0_dp) CYCLE
               radial_work(i) = envelope(i)*scaled_bessel(INT(lambda, C_INT), &
                                                          REAL(2.0_dp*alpha*projectors%r(i)*distance, C_DOUBLE))
            END DO
            radial_power(:) = radial_work
            DO k = 0, degree
               DO j = 1, nb
                  radial(k, j, lambda) = DOT_PRODUCT(radial_power, projectors%weighted_beta(:, j))
               END DO
               radial_power(:) = radial_power*projectors%r
            END DO
         END DO
         centered = 0.0_dp
         DO p = 1, projectors%nproj
            i = projectors%radial_index(p)
            j = projectors%angular_index(p)
            DO s = 1, nc
               k = SUM(indco(:, s))
               centered(s, p) = fourpi*DOT_PRODUCT(angular(s, j, :), radial(k, i, :))
            END DO
         END DO
         cartesian(:, :) = MATMUL(shift, centered)
         DO id = 1, ncoset(nder)
            deriv = indco(:, id)
            DO t = 1, ncoset(la)
               a = indco(:, t)
               coeff = 0.0_dp
               DO axis = 1, 3
                  coeff(a(axis), axis) = 1.0_dp
                  IF (PRESENT(moment_origin)) THEN
                     coeff(:, axis) = 0.0_dp
                     DO power = 0, deriv(axis)
                        coeff(a(axis) + power, axis) = binomial(deriv(axis), power)* &
                                                       moment_origin(axis)**(deriv(axis) - power)
                     END DO
                  ELSE
                     DO step = 1, deriv(axis)
                        next_coeff = 0.0_dp
                        DO power = 0, a(axis) + step - 1
                           next_coeff(power + 1) = next_coeff(power + 1) + 2.0_dp*alpha*coeff(power, axis)
                           IF (power > 0) next_coeff(power - 1) = next_coeff(power - 1) - power*coeff(power, axis)
                        END DO
                        coeff(:, axis) = next_coeff
                     END DO
                  END IF
               END DO
               DO s = 1, nc
                  shifted = indco(:, s)
                  factor = coeff(shifted(1), 1)*coeff(shifted(2), 2)*coeff(shifted(3), 3)
                  IF (factor == 0.0_dp) CYCLE
                  integrals((ip - 1)*ncoset(la) + t, :, id) = integrals((ip - 1)*ncoset(la) + t, :, id) + &
                                                              factor*cartesian(s, :)
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE upf_projector_overlap

END MODULE upf_projector_integrals
