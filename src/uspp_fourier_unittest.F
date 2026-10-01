!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM uspp_fourier_unittest
   USE atom_upf,                        ONLY: atom_read_upf,&
                                              atom_release_upf,&
                                              atom_upfpot_type
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: dp
   USE libcp2k,                         ONLY: cp2k_finalize,&
                                              cp2k_init
   USE mathconstants,                   ONLY: pi
   USE orbital_pointers,                ONLY: indco,&
                                              nco,&
                                              ncoset
   USE orbital_transformation_matrices, ONLY: orbtramat
   USE upf_projector_integrals,         ONLY: upf_projector_init,&
                                              upf_projector_type
   USE uspp_augmentation,               ONLY: &
        uspp_augmentation_collocate, uspp_augmentation_collocate_atom, uspp_augmentation_fourier, &
        uspp_augmentation_init, uspp_augmentation_integrate, uspp_augmentation_integrate_atom, &
        uspp_augmentation_prepare_realspace, uspp_augmentation_realspace, uspp_augmentation_type

   IMPLICIT NONE
   INTEGER, PARAMETER :: nb = 5, np = 19, nao = 7
   REAL(KIND=dp), PARAMETER :: alpha = 1.1_dp, volume = 137.0_dp, center(3) = [0.31_dp, -0.43_dp, 1.17_dp]
   TYPE(atom_upfpot_type) :: pot
   TYPE(upf_projector_type) :: projectors
   TYPE(uspp_augmentation_type) :: aug
   INTEGER :: grid, nr, i, j, l, pair, test
   CHARACTER(len=1024) :: filename
   LOGICAL :: source_root_argument
   REAL(KIND=dp) :: dx, g(3), error, max_error, vectors(3, 7), d(np, np), matrix(np, np), matrix_half(np, np)
   REAL(KIND=dp) :: weights(7), energy_grid, energy_matrix, b(nao, np), p(nao, nao), delta_s(nao, nao)
   COMPLEX(KIND=dp) :: qg(np, np), qminus(np, np), reference(np, np), rho(7), vg(7)
   COMPLEX(KIND=dp) :: strain(np, np, 3, 3), strain_reference(np, np, 3, 3)

   CALL cp2k_init()
   max_error = 0.0_dp
   DO grid = 1, 3
      pot%number_of_proj = nb
      pot%is_ultrasoft = .TRUE.
      ALLOCATE (pot%lbeta(nb), pot%dion(nb, nb), pot%qij(nb, nb))
      pot%lbeta(:) = [0, 1, 2, 3, 1]
      pot%dion = 0.0_dp
      ! Deliberately inconsistent PP_Q: the metric must use the implemented Q(G=0).
      pot%qij = 42.0_dp
      nr = 2501
      IF (grid == 2) nr = 2500
      ALLOCATE (pot%r(nr), pot%rab(nr), pot%beta(nr, nb), pot%qijl(nr, nb*(nb + 1)/2, 0:6))
      DO i = 1, nr
         IF (grid == 3) THEN
            dx = LOG(1.0_dp + 12.5_dp/0.001_dp)/(nr - 1)
            pot%r(i) = 0.001_dp*(EXP((i - 1)*dx) - 1.0_dp)
            pot%rab(i) = dx*(pot%r(i) + 0.001_dp)
         ELSE
            pot%r(i) = (i - 1)*0.005_dp
            pot%rab(i) = 0.005_dp
         END IF
      END DO
      pot%qijl = 0.0_dp
      DO j = 1, nb
         pot%beta(:, j) = pot%r**pot%lbeta(j)*EXP(-alpha*pot%r**2)
         DO i = 1, j
            pair = j*(j - 1)/2 + i
            DO l = ABS(pot%lbeta(i) - pot%lbeta(j)), pot%lbeta(i) + pot%lbeta(j), 2
               pot%qijl(:, pair, l) = (i + j)*0.1_dp*pot%r**(pot%lbeta(i) + pot%lbeta(j) + 2) &
                                      *EXP(-alpha*pot%r**2)
            END DO
         END DO
      END DO
      CALL upf_projector_init(pot, 3, projectors)
      CALL uspp_augmentation_init(pot, projectors, aug)
      IF (aug%nproj /= np) ERROR STOP 'Wrong magnetic-projector count'
      DO test = 1, 5
         SELECT CASE (test)
         CASE (1)
            g = 0.0_dp
         CASE (2)
            g = [1.e-8_dp, -2.e-8_dp, 3.e-8_dp]
         CASE (3)
            g = [0.2_dp, -0.7_dp, 0.9_dp]
         CASE (4)
            g = [3.7_dp, -2.3_dp, 1.5_dp]
         CASE (5)
            g = [0.0_dp, 0.0_dp, 2.8_dp]
         END SELECT
         CALL uspp_augmentation_fourier(aug, g, center, volume, qg, strain)
         CALL gaussian_reference(g, reference, strain_reference)
         error = volume*MAXVAL(ABS(qg - reference))
         max_error = MAX(max_error, error)
         IF (error > 2.e-10_dp) THEN
            WRITE (*, *) 'Gaussian Fourier error:', grid, test, error
            ERROR STOP 'USPP Fourier transform differs from Cartesian Gaussian integral'
         END IF
         error = volume*MAXVAL(ABS(strain - strain_reference))
         max_error = MAX(max_error, error)
         IF (error > 2.e-10_dp) THEN
            WRITE (*, *) 'Gaussian Fourier strain error:', grid, test, error
            ERROR STOP 'USPP strain derivative differs from Cartesian Gaussian integral'
         END IF
         CALL uspp_augmentation_fourier(aug, -g, center, volume, qminus)
         IF (MAXVAL(ABS(qminus - CONJG(qg))) > 2.e-13_dp) ERROR STOP 'Wrong reciprocal conjugation'
         IF (MAXVAL(ABS(qg - TRANSPOSE(qg))) > 2.e-13_dp) ERROR STOP 'Wrong projector symmetry'
      END DO

      vectors(:, 1) = 0.0_dp
      vectors(:, 2) = [0.3_dp, -0.7_dp, 1.1_dp]
      vectors(:, 3) = [2.1_dp, 0.2_dp, -0.4_dp]
      vectors(:, 4) = [0.0_dp, 0.0_dp, 2.8_dp]
      vectors(:, 5:7) = -vectors(:, 2:4)
      DO j = 1, np
         DO i = 1, nao
            b(i, j) = 0.03_dp*SIN(REAL(i*j, dp))
         END DO
      END DO
      p = 0.01_dp
      DO i = 1, nao
         p(i, i) = 0.3_dp
      END DO
      d = MATMUL(TRANSPOSE(b), MATMUL(p, b))
      delta_s = MATMUL(b, MATMUL(aug%q, TRANSPOSE(b)))
      rho = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      CALL uspp_augmentation_collocate(aug, vectors, center, volume, d, rho)
      IF (ABS(volume*REAL(rho(1), dp) - SUM(p*delta_s)) > 1.e-13_dp) ERROR STOP 'Charge/metric mismatch'
      vg(1) = CMPLX(0.37_dp, 0.0_dp, KIND=dp)
      vg(2:4) = [CMPLX(0.2_dp, -0.3_dp, dp), CMPLX(-0.5_dp, 0.8_dp, dp), CMPLX(1.1_dp, 0.4_dp, dp)]
      vg(5:7) = CONJG(vg(2:4))
      weights = 1.0_dp
      CALL uspp_augmentation_integrate(aug, vectors, center, volume, vg, weights, matrix)
      energy_grid = volume*REAL(SUM(CONJG(vg)*rho), dp)
      energy_matrix = SUM(d*matrix)
      IF (ABS(energy_grid - energy_matrix) > 1.e-13_dp) ERROR STOP 'Augmentation operations are not adjoint'
      weights(2:4) = 2.0_dp
      CALL uspp_augmentation_integrate(aug, vectors(:, :4), center, volume, vg(:4), weights(:4), matrix_half)
      IF (MAXVAL(ABS(matrix_half - matrix)) > 1.e-13_dp) ERROR STOP 'Wrong reciprocal multiplicities'
      vg = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      vg(1) = CMPLX(2.3_dp, 0.0_dp, KIND=dp)
      CALL uspp_augmentation_integrate(aug, vectors, center, volume, vg, weights, matrix)
      IF (MAXVAL(ABS(matrix - 2.3_dp*aug%q)) > 1.e-13_dp) ERROR STOP 'Wrong constant-potential shift'
      CALL check_realspace()
      CALL check_compact_support()
      CALL atom_release_upf(pot)
   END DO
   WRITE (*, *) 'USPP Fourier: 15 analytic s/p/d/f comparisons, maximum error =', max_error
   DO i = 1, COMMAND_ARGUMENT_COUNT()
      CALL GET_COMMAND_ARGUMENT(i, filename)
      IF (i == 1) THEN
         INQUIRE (FILE=TRIM(filename)//"/tests/UNIT_TESTS", EXIST=source_root_argument)
         IF (source_root_argument) CYCLE
      END IF
      CALL check_upf_realspace(TRIM(filename))
   END DO
   ! Release allocatables before LeakSanitizer inspects the main program at exit.
   aug = uspp_augmentation_type()
   projectors = upf_projector_type()
   CALL cp2k_finalize()

CONTAINS

! **************************************************************************************************
!> \brief Check source-node reproduction and finite center derivatives for optional UPF files.
!> \param filename Scalar USPP dataset.
! **************************************************************************************************
   SUBROUTINE check_upf_realspace(filename)
      CHARACTER(len=*), INTENT(IN)                       :: filename

      INTEGER                                            :: i, j, l, node, p, pair, point, q
      REAL(KIND=dp)                                      :: displacement(3), error, reference, scale
      REAL(KIND=dp), ALLOCATABLE                         :: gradient(:, :, :), hessian(:, :, :, :), &
                                                            value(:, :)
      TYPE(atom_upfpot_type)                             :: dataset
      TYPE(upf_projector_type)                           :: beta
      TYPE(uspp_augmentation_type)                       :: augmentation

      CALL atom_read_upf(dataset, filename, read_augmentation=.TRUE.)
      CALL upf_projector_init(dataset, MAXVAL(dataset%lbeta), beta)
      CALL uspp_augmentation_init(dataset, beta, augmentation)
      CALL uspp_augmentation_prepare_realspace(dataset, augmentation)
      p = augmentation%nproj
      ALLOCATE (value(p, p), gradient(p, p, 3), hessian(p, p, 3, 3))
      displacement = 0.0_dp
      CALL uspp_augmentation_realspace(augmentation, displacement, value, gradient, hessian)
      IF (.NOT. ALL(ieee_is_finite(value)) .OR. .NOT. ALL(ieee_is_finite(gradient)) .OR. &
          .NOT. ALL(ieee_is_finite(hessian))) ERROR STOP 'Nonfinite augmentation at origin'
      error = 0.0_dp
      scale = 1.0_dp
      DO point = 0, 8
         node = MAX(2, 1 + (SIZE(dataset%r) - 1)*point/8)
         displacement = [0.0_dp, 0.0_dp, dataset%r(node)]
         CALL uspp_augmentation_realspace(augmentation, displacement, value, gradient, hessian)
         IF (.NOT. ALL(ieee_is_finite(value)) .OR. .NOT. ALL(ieee_is_finite(gradient)) .OR. &
             .NOT. ALL(ieee_is_finite(hessian))) ERROR STOP 'Nonfinite atom-grid augmentation'
         DO q = 1, augmentation%nproj
            j = augmentation%radial_index(q)
            DO p = 1, augmentation%nproj
               i = augmentation%radial_index(p)
               pair = MAX(i, j)*(MAX(i, j) - 1)/2 + MIN(i, j)
               reference = 0.0_dp
               DO l = ABS(dataset%lbeta(i) - dataset%lbeta(j)), dataset%lbeta(i) + dataset%lbeta(j), 2
                  reference = reference + augmentation%gaunt(p, q, l*l + l + 1)* &
                              SQRT(REAL(2*l + 1, dp)/(4.0_dp*pi))*dataset%qijl(node, pair, l)/dataset%r(node)**2
               END DO
               error = MAX(error, ABS(value(p, q) - reference))
               scale = MAX(scale, ABS(reference))
            END DO
         END DO
      END DO
      WRITE (*, *) 'USPP real-space table error:', error/scale, filename
      IF (error/scale > 2.e-10_dp) ERROR STOP 'Atom-grid interpolation changes source nodes'
      CALL atom_release_upf(dataset)
   END SUBROUTINE check_upf_realspace

! **************************************************************************************************
!> \brief Atom-grid values and derivatives against independent Cartesian Gaussian products.
! **************************************************************************************************
   SUBROUTINE check_realspace()
      REAL(KIND=dp), PARAMETER                           :: step = 1.e-6_dp

      INTEGER                                            :: axis, point
      REAL(KIND=dp) :: coordinate(3), derivative(np, np, 3), displaced(3), errors(3), &
         exact(np, np), exact_d(np, np, 3), exact_dd(np, np, 3, 3), finite_difference_error, &
         minus(np, np), minus_d(np, np, 3), plus(np, np), plus_d(np, np, 3), points(3, 7), &
         second(np, np, 3, 3), value(np, np)

      CALL uspp_augmentation_prepare_realspace(pot, aug)
      points(:, 1) = 0.0_dp
      points(:, 2) = [1.e-8_dp, -2.e-8_dp, 3.e-8_dp]
      points(:, 3) = [0.2_dp, -0.7_dp, 0.9_dp]
      points(:, 4) = [0.0_dp, 0.0_dp, 1.37_dp]
      points(:, 5) = [1.23_dp, 0.0_dp, 0.0_dp]
      points(:, 6) = [0.0_dp, pot%r(SIZE(pot%r)/2), 0.0_dp]
      points(:, 7) = [0.0_dp, 0.0_dp, 1.1_dp*MAXVAL(pot%r)]
      errors = 0.0_dp
      finite_difference_error = 0.0_dp
      DO point = 1, SIZE(points, 2)
         coordinate = points(:, point)
         CALL uspp_augmentation_realspace(aug, coordinate, value, derivative, second)
         CALL gaussian_realspace(coordinate, exact, exact_d, exact_dd)
         errors(1) = MAX(errors(1), MAXVAL(ABS(value - exact)))
         errors(2) = MAX(errors(2), MAXVAL(ABS(derivative - exact_d)))
         errors(3) = MAX(errors(3), MAXVAL(ABS(second - exact_dd)))
         DO axis = 1, 3
            displaced = coordinate
            displaced(axis) = displaced(axis) + step
            CALL uspp_augmentation_realspace(aug, displaced, plus, plus_d)
            displaced(axis) = displaced(axis) - 2.0_dp*step
            CALL uspp_augmentation_realspace(aug, displaced, minus, minus_d)
            finite_difference_error = MAX(finite_difference_error, &
                                          MAXVAL(ABS((plus - minus)/(2.0_dp*step) - derivative(:, :, axis))), &
                                          MAXVAL(ABS((plus_d - minus_d)/(2.0_dp*step) - second(:, :, :, axis))))
         END DO
      END DO
      WRITE (*, *) 'USPP real-space value/gradient/Hessian errors:', grid, errors
      WRITE (*, *) 'USPP real-space finite-difference error:', grid, finite_difference_error
      IF (errors(1) > 2.e-8_dp .OR. errors(2) > 2.e-6_dp .OR. errors(3) > 2.e-4_dp) THEN
         ERROR STOP 'Atom-grid augmentation differs from Cartesian Gaussian products'
      END IF
      IF (finite_difference_error > 5.e-8_dp) ERROR STOP 'Inconsistent atom-grid augmentation derivatives'
      CALL check_atom_adjoint(points(:, 1:5))
   END SUBROUTINE check_realspace

! **************************************************************************************************
!> \brief Preserve exact zero suffixes without losing C2 continuity or energy derivatives.
! **************************************************************************************************
   SUBROUTINE check_compact_support()
      INTEGER                                            :: axis, node, point
      REAL(dp) :: coordinate(3), err, gradient(np, np, 3), hessian(np, np, 3, 3), join(3), &
         minus(np, np), minus_d(np, np, 3), plus(np, np), plus_d(np, np, 3), radius, step, &
         value(np, np)
      TYPE(atom_upfpot_type), ALLOCATABLE                :: compact
      TYPE(uspp_augmentation_type)                       :: compact_aug

      ALLOCATE (compact)
      compact = pot
      node = MINLOC(ABS(pot%r - 1.4_dp), DIM=1)
      radius = pot%r(node)
      DO point = 1, SIZE(pot%r)
         compact%qijl(point, :, :) = compact%qijl(point, :, :)*MAX(0.0_dp, 1.0_dp - (pot%r(point)/radius)**2)**3
      END DO
      CALL uspp_augmentation_init(compact, projectors, compact_aug)
      CALL uspp_augmentation_prepare_realspace(compact, compact_aug)
      IF (compact_aug%realspace_radius /= radius) ERROR STOP 'Incorrect compact augmentation support'
      join = [pot%r(node - 1), 0.5_dp*(pot%r(node - 1) + radius), radius]
      step = 1.e-7_dp
      err = 0.0_dp
      DO point = 1, SIZE(join)
         coordinate = [0.0_dp, 0.0_dp, join(point)]
         CALL uspp_augmentation_realspace(compact_aug, coordinate, value, gradient, hessian)
         IF (point == 3) THEN
            IF (MAXVAL(ABS(value)) + MAXVAL(ABS(gradient)) + MAXVAL(ABS(hessian)) /= 0.0_dp) THEN
               ERROR STOP 'Augmentation closure does not vanish through second derivatives'
            END IF
         END IF
         DO axis = 1, 3
            coordinate(axis) = coordinate(axis) + step
            CALL uspp_augmentation_realspace(compact_aug, coordinate, plus, plus_d)
            coordinate(axis) = coordinate(axis) - 2.0_dp*step
            CALL uspp_augmentation_realspace(compact_aug, coordinate, minus, minus_d)
            coordinate(axis) = coordinate(axis) + step
            err = MAX(err, MAXVAL(ABS((plus - minus)/(2.0_dp*step) - gradient(:, :, axis))), &
                      MAXVAL(ABS((plus_d - minus_d)/(2.0_dp*step) - hessian(:, :, :, axis))))
         END DO
      END DO
      coordinate = [0.0_dp, 0.0_dp, radius + 0.2_dp]
      CALL uspp_augmentation_realspace(compact_aug, coordinate, value, gradient, hessian)
      IF (MAXVAL(ABS(value)) + MAXVAL(ABS(gradient)) + MAXVAL(ABS(hessian)) /= 0.0_dp) THEN
         ERROR STOP 'Nonzero augmentation in an exact tabulated zero suffix'
      END IF
      WRITE (*, *) 'USPP compact-support derivative error:', grid, err
      IF (err > 5.e-6_dp) ERROR STOP 'Inconsistent derivatives at compact augmentation joins'
      CALL atom_release_upf(compact)
   END SUBROUTINE check_compact_support

! **************************************************************************************************
!> \brief Check the signed-quadrature adjoint and nonlinear density/center energy derivatives.
!> \param points Atom-grid displacements, including the origin and axis points.
! **************************************************************************************************
   SUBROUTINE check_atom_adjoint(points)
      REAL(KIND=dp), INTENT(IN)                          :: points(3, 5)

      REAL(KIND=dp), PARAMETER :: step = 1.e-5_dp, &
         weights(5) = [0.3_dp, -0.7_dp, 1.1_dp, 0.2_dp, -0.4_dp]

      INTEGER                                            :: axis, i, j
      REAL(KIND=dp) :: center_gradient(3), eminus, eplus, error, gradient(3, 5), matrix(np, np), &
         pairing, rho(5), shifted(3, 5), trial(np, np)

      rho = 0.4_dp
      gradient = 0.05_dp
      CALL uspp_augmentation_collocate_atom(aug, points, d, rho, gradient)
      CALL uspp_augmentation_integrate_atom(aug, points, weights, rho**2, gradient, matrix, d, center_gradient)
      pairing = SUM(weights*(rho**2*(rho - 0.4_dp) + SUM(gradient*(gradient - 0.05_dp), DIM=1)))
      error = ABS(pairing - SUM(d*matrix))
      IF (error > 1.e-13_dp) ERROR STOP 'Atom-grid augmentation operations are not adjoint'
      DO j = 1, np
         DO i = 1, np
            trial(i, j) = 0.03_dp*COS(REAL(i*j, dp))
         END DO
      END DO
      eplus = atom_test_energy(points, d + step*trial, weights)
      eminus = atom_test_energy(points, d - step*trial, weights)
      error = MAX(error, ABS((eplus - eminus)/(2.0_dp*step) - SUM(trial*matrix)))
      DO axis = 1, 3
         shifted = points
         shifted(axis, :) = shifted(axis, :) - step
         eplus = atom_test_energy(shifted, d, weights)
         shifted(axis, :) = shifted(axis, :) + 2.0_dp*step
         eminus = atom_test_energy(shifted, d, weights)
         error = MAX(error, ABS((eplus - eminus)/(2.0_dp*step) - center_gradient(axis)))
      END DO
      WRITE (*, *) 'USPP atom-grid nonlinear derivative error:', grid, error
      IF (error > 1.e-9_dp) ERROR STOP 'Inconsistent atom-grid nonlinear energy derivatives'
   END SUBROUTINE check_atom_adjoint

! **************************************************************************************************
!> \brief Nonlinear surrogate energy with a fixed environmental density and gradient.
!> \param points Atom-grid displacements.
!> \param density On-site density matrix.
!> \param weights Signed quadrature weights.
!> \return Sum of n**3/3 + |grad(n)|**2/2 on the same quadrature.
! **************************************************************************************************
   FUNCTION atom_test_energy(points, density, weights) RESULT(energy)
      REAL(KIND=dp), INTENT(IN)                          :: points(3, 5), density(np, np), weights(5)
      REAL(KIND=dp)                                      :: energy

      REAL(KIND=dp)                                      :: gradient(3, 5), rho(5)

      rho = 0.4_dp
      gradient = 0.05_dp
      CALL uspp_augmentation_collocate_atom(aug, points, density, rho, gradient)
      energy = SUM(weights*(rho**3/3.0_dp + 0.5_dp*SUM(gradient**2, DIM=1)))
   END FUNCTION atom_test_energy

! **************************************************************************************************
!> \brief Independent monomial product rule, without Gaunt coefficients or radial interpolation.
!> \param coordinate Cartesian displacement from the center.
!> \param value Gaussian augmentation matrix.
!> \param derivative Cartesian gradient.
!> \param second Cartesian Hessian.
! **************************************************************************************************
   SUBROUTINE gaussian_realspace(coordinate, value, derivative, second)
      REAL(KIND=dp), INTENT(IN)                          :: coordinate(3)
      REAL(KIND=dp), INTENT(OUT)                         :: value(np, np), derivative(np, np, 3), &
                                                            second(np, np, 3, 3)

      INTEGER                                            :: a, axis, axis2, c, i, j, k, li, lj, mi, &
                                                            mj, powers(3), u, v
      REAL(KIND=dp)                                      :: coefficient, factors(3), first(0:6, 3), &
                                                            functions(0:8, 3), twice(0:6, 3)

      DO axis = 1, 3
         DO k = 0, 8
            functions(k, axis) = coordinate(axis)**k*EXP(-alpha*coordinate(axis)**2)
         END DO
         DO k = 0, 6
            first(k, axis) = -2.0_dp*alpha*functions(k + 1, axis)
            IF (k > 0) first(k, axis) = first(k, axis) + k*functions(k - 1, axis)
            twice(k, axis) = -2.0_dp*alpha*(2*k + 1)*functions(k, axis) + &
                             4.0_dp*alpha**2*functions(k + 2, axis)
            IF (k > 1) twice(k, axis) = twice(k, axis) + k*(k - 1)*functions(k - 2, axis)
         END DO
      END DO
      value = 0.0_dp
      derivative = 0.0_dp
      second = 0.0_dp
      DO v = 1, np
         j = projectors%radial_index(v)
         lj = pot%lbeta(j)
         mj = projectors%angular_index(v) - lj*lj
         DO u = 1, np
            i = projectors%radial_index(u)
            li = pot%lbeta(i)
            mi = projectors%angular_index(u) - li*li
            DO c = 1, nco(lj)
               DO a = 1, nco(li)
                  powers = indco(:, ncoset(li - 1) + a) + indco(:, ncoset(lj - 1) + c)
                  coefficient = (i + j)*0.1_dp*orbtramat(li)%slm(mi, a)*orbtramat(lj)%slm(mj, c)
                  factors = [functions(powers(1), 1), functions(powers(2), 2), functions(powers(3), 3)]
                  value(u, v) = value(u, v) + coefficient*PRODUCT(factors)
                  DO axis = 1, 3
                     factors(axis) = first(powers(axis), axis)
                     derivative(u, v, axis) = derivative(u, v, axis) + coefficient*PRODUCT(factors)
                     DO axis2 = 1, 3
                        IF (axis2 == axis) THEN
                           factors(axis2) = twice(powers(axis2), axis2)
                        ELSE
                           factors(axis2) = first(powers(axis2), axis2)
                        END IF
                        second(u, v, axis, axis2) = second(u, v, axis, axis2) + coefficient*PRODUCT(factors)
                        factors(axis2) = functions(powers(axis2), axis2)
                        factors(axis) = first(powers(axis), axis)
                     END DO
                     factors(axis) = functions(powers(axis), axis)
                  END DO
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE gaussian_realspace

! **************************************************************************************************
!> \brief Independent full Cartesian Gaussian Fourier integrals, without Gaunt or Bessel sums.
!> \param g Reciprocal vector.
!> \param fourier Analytic transform of solid-harmonic products times a Gaussian.
!> \param strain Independent strain derivative from one-higher Cartesian moments.
! **************************************************************************************************
   SUBROUTINE gaussian_reference(g, fourier, strain)
      REAL(KIND=dp), INTENT(IN)                          :: g(3)
      COMPLEX(KIND=dp), INTENT(OUT)                      :: fourier(np, np), strain(np, np, 3, 3)

      COMPLEX(KIND=dp)                                   :: dterm(3), moments(0:7, 3), phase, term
      INTEGER                                            :: a, axis, c, i, j, k, li, lj, mi, mj, &
                                                            powers(3), raised(3), u, v

      DO axis = 1, 3
         moments(0, axis) = SQRT(pi/alpha)*EXP(-g(axis)**2/(4.0_dp*alpha))
         moments(1, axis) = CMPLX(0.0_dp, -g(axis)/(2.0_dp*alpha), dp)*moments(0, axis)
         DO k = 1, 6
            moments(k + 1, axis) = REAL(k, dp)/(2.0_dp*alpha)*moments(k - 1, axis) &
                                   + CMPLX(0.0_dp, -g(axis)/(2.0_dp*alpha), dp)*moments(k, axis)
         END DO
      END DO
      phase = EXP(CMPLX(0.0_dp, -DOT_PRODUCT(g, center), dp))/volume
      fourier = CMPLX(0.0_dp, 0.0_dp, dp)
      strain = CMPLX(0.0_dp, 0.0_dp, dp)
      DO v = 1, np
         j = projectors%radial_index(v)
         lj = pot%lbeta(j)
         mj = projectors%angular_index(v) - lj*lj
         DO u = 1, np
            i = projectors%radial_index(u)
            li = pot%lbeta(i)
            mi = projectors%angular_index(u) - li*li
            DO c = 1, nco(lj)
               DO a = 1, nco(li)
                  powers = indco(:, ncoset(li - 1) + a) + indco(:, ncoset(lj - 1) + c)
                  term = PRODUCT([moments(powers(1), 1), moments(powers(2), 2), moments(powers(3), 3)])
                  fourier(u, v) = fourier(u, v) + orbtramat(li)%slm(mi, a)*orbtramat(lj)%slm(mj, c)*term
                  DO axis = 1, 3
                     raised = powers
                     raised(axis) = raised(axis) + 1
                     dterm(axis) = CMPLX(0.0_dp, -1.0_dp, dp)* &
                                   PRODUCT([moments(raised(1), 1), moments(raised(2), 2), moments(raised(3), 3)])
                     strain(u, v, :, axis) = strain(u, v, :, axis) - &
                                             g*dterm(axis)*orbtramat(li)%slm(mi, a)*orbtramat(lj)%slm(mj, c)
                  END DO
               END DO
            END DO
            fourier(u, v) = fourier(u, v)*phase*(i + j)*0.1_dp
            strain(u, v, :, :) = strain(u, v, :, :)*phase*(i + j)*0.1_dp
            DO axis = 1, 3
               strain(u, v, axis, axis) = strain(u, v, axis, axis) - fourier(u, v)
            END DO
         END DO
      END DO
   END SUBROUTINE gaussian_reference
END PROGRAM uspp_fourier_unittest
