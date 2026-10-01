!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM atom_uspp_unittest
   USE atom_operators,                  ONLY: atom_int_release,&
                                              atom_int_setup,&
                                              atom_ppint_release,&
                                              atom_ppint_setup
   USE atom_types,                      ONLY: atom_basis_gridrep,&
                                              atom_basis_type,&
                                              atom_integrals,&
                                              atom_potential_type,&
                                              init_atom_basis_default_pp,&
                                              lmat,&
                                              release_atom_basis,&
                                              upf_pseudo
   USE atom_upf,                        ONLY: atom_release_upf
   USE atom_utils,                      ONLY: atom_density,&
                                              numpot_matrix
   USE input_constants,                 ONLY: poly_conf
   USE kinds,                           ONLY: dp
   USE libcp2k,                         ONLY: cp2k_finalize,&
                                              cp2k_init
   USE mathconstants,                   ONLY: fourpi

   IMPLICIT NONE
   INTEGER, PARAMETER :: nr = 6001, ng = 6000
   REAL(dp), PARAMETER :: dr = 0.002_dp, alpha(2) = [0.4_dp, 1.5_dp], gamma_beta(3) = [0.6_dp, 1.7_dp, 0.9_dp]
   TYPE(atom_basis_type) :: template, basis
   TYPE(atom_potential_type) :: potential
   TYPE(atom_integrals) :: integrals
   INTEGER :: i, j, k, l, p, q, pair, order
   REAL(dp) :: r(ng), wr(ng), rho(ng), v(ng), errors(5), b(2, 3), coefficient(6), eta(6), &
               matrix(2, 2, 0:lmat), pmat(2, 2, 0:lmat), ref_s(2, 2, 0:lmat), ref_h(2, 2, 0:lmat), &
               ref_ovlp(2, 2, 0:lmat), ref_conf(2, 2, 0:lmat), ref_nl(2, 2, 0:lmat), &
               transformed(2, 2), egrid, ematrix, ref_b, scalar

   CALL cp2k_init()
   CALL init_atom_basis_default_pp(template)
   template%nbas = 0
   template%nprim = 0
   template%nbas(0:1) = 2
   template%nprim(0:1) = 2
   template%am(1, :) = alpha(1)
   template%am(2, :) = alpha(2)
   DO k = 1, ng
      r(k) = dr*(REAL(ng - k, dp) + 0.5_dp)
   END DO
   wr = dr
   CALL atom_basis_gridrep(template, basis, r, wr)
   CALL release_atom_basis(template)
   potential%ppot_type = upf_pseudo
   potential%confinement = .TRUE.
   potential%conf_type = poly_conf
   potential%rcon = 2.0_dp
   potential%scon = 2.0_dp
   potential%upf_pot%is_ultrasoft = .TRUE.
   potential%upf_pot%pseudo_type = 'USPP'
   potential%upf_pot%number_of_proj = 3
   potential%upf_pot%mesh_size = nr
   potential%upf_pot%l_max = 1
   ALLOCATE (potential%upf_pot%r(nr), potential%upf_pot%rab(nr), potential%upf_pot%vlocal(nr), &
             potential%upf_pot%beta(nr, 3), potential%upf_pot%lbeta(3), potential%upf_pot%dion(3, 3), &
             potential%upf_pot%qijl(nr, 6, 0:2))
   potential%upf_pot%lbeta(:) = [0, 0, 1]
   potential%upf_pot%dion = 0.0_dp
   potential%upf_pot%dion(1, 1) = 0.3_dp
   potential%upf_pot%dion(1, 2) = -0.1_dp
   potential%upf_pot%dion(2, 1) = -0.1_dp
   potential%upf_pot%dion(2, 2) = 0.6_dp
   potential%upf_pot%dion(3, 3) = -0.2_dp
   potential%upf_pot%qijl = 0.0_dp
   coefficient = [0.2_dp, -0.035_dp, 0.13_dp, 0.0_dp, 0.0_dp, 0.17_dp]
   eta = [0.7_dp, 1.1_dp, 1.6_dp, 1.0_dp, 1.0_dp, 0.9_dp]
   DO k = 1, nr
      scalar = (k - 1)*dr
      potential%upf_pot%r(k) = scalar
      potential%upf_pot%rab(k) = dr
      potential%upf_pot%vlocal(k) = -0.7_dp*EXP(-0.5_dp*scalar**2)
      DO p = 1, 3
         potential%upf_pot%beta(k, p) = scalar**potential%upf_pot%lbeta(p)*EXP(-gamma_beta(p)*scalar**2)
      END DO
      DO pair = 1, 6
         potential%upf_pot%qijl(k, pair, 0) = coefficient(pair)*scalar**2*EXP(-eta(pair)*scalar**2)
      END DO
   END DO
   CALL atom_int_setup(integrals, basis, potential=potential)
   CALL atom_ppint_setup(integrals, basis, potential)
   errors = 0.0_dp
   b = 0.0_dp
   ref_s = 0.0_dp
   ref_h = 0.0_dp
   ref_ovlp = 0.0_dp
   ref_conf = 0.0_dp
   ref_nl = 0.0_dp
   DO p = 1, 3
      l = potential%upf_pot%lbeta(p)
      DO i = 1, 2
         ref_b = radial_integral(l, alpha(i) + gamma_beta(p))
         b(i, p) = ref_b
         errors(1) = MAX(errors(1), ABS(basis%uspp%beta(i, p) - ref_b))
      END DO
   END DO
   DO l = 0, 1
      DO i = 1, 2
         DO j = 1, 2
            ref_ovlp(i, j, l) = radial_integral(l, alpha(i) + alpha(j))
            ref_s(i, j, l) = ref_ovlp(i, j, l)
            ref_h(i, j, l) = -0.7_dp*radial_integral(l, alpha(i) + alpha(j) + 0.5_dp)
            ref_conf(i, j, l) = 0.25_dp*radial_integral(l + 1, alpha(i) + alpha(j))
            DO p = 1, 3
               IF (potential%upf_pot%lbeta(p) /= l) CYCLE
               DO q = 1, 3
                  IF (potential%upf_pot%lbeta(q) /= l) CYCLE
                  pair = MAX(p, q)*(MAX(p, q) - 1)/2 + MIN(p, q)
                  ref_s(i, j, l) = ref_s(i, j, l) + b(i, p)*b(j, q)*coefficient(pair)*radial_integral(0, eta(pair))
                  ref_h(i, j, l) = ref_h(i, j, l) - &
                                   0.7_dp*b(i, p)*b(j, q)*coefficient(pair)*radial_integral(0, eta(pair) + 0.5_dp)
                  ref_conf(i, j, l) = ref_conf(i, j, l) + &
                                      0.25_dp*b(i, p)*b(j, q)*coefficient(pair)*radial_integral(1, eta(pair))
                  ref_nl(i, j, l) = ref_nl(i, j, l) + b(i, p)*potential%upf_pot%dion(p, q)*b(j, q)
               END DO
            END DO
         END DO
      END DO
      transformed = MATMUL(TRANSPOSE(integrals%utrans(:, :, l)), &
                           MATMUL(integrals%metric(:, :, l), integrals%utrans(:, :, l)))
      transformed(1, 1) = transformed(1, 1) - 1.0_dp
      transformed(2, 2) = transformed(2, 2) - 1.0_dp
      errors(2) = MAX(errors(2), MAXVAL(ABS(transformed)))
   END DO
   errors(1) = MAX(errors(1), MAXVAL(ABS(integrals%ovlp - ref_ovlp)))
   errors(2) = MAX(errors(2), MAXVAL(ABS(integrals%metric - ref_s)))
   errors(3) = MAXVAL(ABS(integrals%core - ref_h))
   errors(3) = MAX(errors(3), MAXVAL(ABS(integrals%hnl - ref_nl)), MAXVAL(ABS(integrals%conf - ref_conf)))
   pmat = 0.0_dp
   pmat(1, 1, 0:1) = 0.8_dp
   pmat(2, 2, 0:1) = 0.45_dp
   pmat(1, 2, 0:1) = -0.1_dp
   pmat(2, 1, 0:1) = -0.1_dp
   CALL atom_density(rho, pmat, basis, 1, typ='RHO')
   errors(4) = ABS(fourpi*SUM(rho*basis%grid%wr) - SUM(pmat*integrals%metric))
   v = EXP(-0.23_dp*r**2)*(1.0_dp + 0.2_dp*r**2)
   DO order = 0, 1
      matrix = 0.0_dp
      CALL numpot_matrix(matrix, v, basis, order)
      IF (order == 0) THEN
         CALL atom_density(rho, pmat, basis, 1, typ='RHO')
      ELSE
         CALL atom_density(rho, pmat, basis, 1, typ='DER')
      END IF
      egrid = fourpi*SUM(rho*v*basis%grid%wr)
      ematrix = SUM(pmat*matrix)
      errors(5) = MAX(errors(5), ABS(egrid - ematrix))
   END DO
   ! Reusing the basis must not retain an obsolete or twice-added augmentation.
   CALL atom_ppint_release(integrals)
   CALL atom_ppint_setup(integrals, basis, potential)
   errors(3) = MAX(errors(3), MAXVAL(ABS(integrals%conf - ref_conf)))
   CALL atom_ppint_release(integrals)
   CALL atom_int_release(integrals)
   CALL atom_int_setup(integrals, basis, potential=potential)
   CALL atom_ppint_setup(integrals, basis, potential)
   errors(2) = MAX(errors(2), MAXVAL(ABS(integrals%metric - ref_s)))
   errors(3) = MAX(errors(3), MAXVAL(ABS(integrals%conf - ref_conf)))
   CALL atom_ppint_release(integrals)
   CALL atom_int_release(integrals)
   WRITE (*, '(A,5ES20.11)') 'Atomic USPP beta/metric/potential/charge/adjoint errors: ', errors
   IF (MAXVAL(errors) > 1.e-8_dp) ERROR STOP 'Atomic USPP analytic reference failed'
   CALL release_atom_basis(basis)
   CALL atom_release_upf(potential%upf_pot)
   CALL cp2k_finalize()
CONTAINS
! **************************************************************************************************
!> \brief Analytic Gaussian radial integral including the volume element.
!> \param l Angular momentum of the two factors.
!> \param exponent Sum of the Gaussian exponents.
!> \return Integral of r**(2*l+2)*exp(-exponent*r**2) from zero to infinity.
! **************************************************************************************************
   REAL(dp) FUNCTION radial_integral(l, exponent) RESULT(value)
      INTEGER, INTENT(IN)                                :: l
      REAL(dp), INTENT(IN)                               :: exponent

      value = 0.5_dp*GAMMA(REAL(l, dp) + 1.5_dp)/exponent**(REAL(l, dp) + 1.5_dp)
   END FUNCTION radial_integral
END PROGRAM atom_uspp_unittest
