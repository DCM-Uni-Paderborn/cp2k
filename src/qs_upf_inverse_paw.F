!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief PAW one-center functional in the native inverse Gaussian representation.
!>        Fixed kinetic and ionic terms remain in the original PP_DIJ.
! **************************************************************************************************
MODULE qs_upf_inverse_paw
   USE compensated_arithmetic,          ONLY: twofold_add,&
                                              twofold_multiply,&
                                              twofold_real_type,&
                                              twofold_value
   USE hfx_types,                       ONLY: hfx_type
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE input_constants,                 ONLY: do_potential_coulomb,&
                                              do_potential_long,&
                                              do_potential_mix_cl,&
                                              do_potential_mix_cl_trunc,&
                                              do_potential_short,&
                                              do_potential_truncated
   USE input_section_types,             ONLY: section_vals_type
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: fourpi
   USE mathlib,                         ONLY: invmat
   USE message_passing,                 ONLY: mp_comm_type
   USE paw_onecenter,                   ONLY: paw_onecenter_core_exchange,&
                                              paw_onecenter_core_exchange_init,&
                                              paw_onecenter_core_exchange_prepare,&
                                              paw_onecenter_hartree,&
                                              paw_onecenter_init
   USE paw_types,                       ONLY: paw_kind_type,&
                                              paw_potential_type
   USE paw_xc,                          ONLY: paw_xc_evaluate,&
                                              paw_xc_init,&
                                              paw_xc_response
   USE pw_grid_types,                   ONLY: FULLSPACE,&
                                              pw_grid_type
   USE qs_upf_inverse,                  ONLY: qs_upf_inverse_response
   USE qs_upf_inverse_types,            ONLY: qs_upf_inverse_augmentation_type,&
                                              qs_upf_inverse_grid_type,&
                                              qs_upf_inverse_orbitals_type,&
                                              qs_upf_inverse_type
   USE skala_density_constraints,       ONLY: skala_density_constraints_type
   USE upf_inverse,                     ONLY: upf_inverse_density_evaluate,&
                                              upf_inverse_kinetic_type
   USE upf_inverse_coulomb,             ONLY: upf_inverse_coulomb_fourier_kernel
   USE upf_inverse_precision,           ONLY: upf_inverse_orbital_projector_density,&
                                              upf_inverse_projector_operator
   USE upf_projector_integrals,         ONLY: upf_projector_type
   USE upf_types,                       ONLY: atom_upfpot_type
   USE uspp_augmentation,               ONLY: uspp_augmentation_realspace,&
                                              uspp_augmentation_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_upf_inverse_paw'

   PUBLIC :: qs_upf_inverse_paw_init, qs_upf_inverse_paw_evaluate, qs_upf_inverse_paw_response
   PUBLIC :: qs_upf_inverse_paw_core_prepare, qs_upf_inverse_paw_xc_kernel
   PUBLIC :: qs_upf_inverse_paw_collocate, qs_upf_inverse_paw_operator
   PUBLIC :: qs_upf_inverse_paw_density_fourier, qs_upf_inverse_paw_local_operator
   PUBLIC :: qs_upf_inverse_paw_density_adjoint

CONTAINS

! **************************************************************************************************
!> \brief Fourier coefficients of the original PAW AE-minus-compensated-PS density and AE core.
!>        Evaluate on the original radial mesh. An optional native full-box grid specifies
!>        real tensor-product boundary coefficients, with matching geometry responses.
!>        Omit that grid for continuous Fourier observables.
!>        The core is shared equally between physical spins and contributes once to their sum.
!> \param kind Prepared original PAW radial products and magnetic-projector Gaunt tensor.
!> \param projected Hermitian spin projector occupations, with k weights already applied once.
!> \param gvecs Local reciprocal vectors, without half-space multiplicities.
!> \param position Home atomic center. The Fourier phase represents its periodic image sum.
!> \param volume Physical cell volume defining the Fourier normalization.
!> \param density Overwritten spin corrections, including the frozen core and zero mode.
!> \param response Optional density derivative for fixed original radial data.
!> \param projected_response Projector-occupation derivative, required with response.
!> \param gvecs_response Reciprocal-vector derivative, required with response.
!> \param position_response Atomic-center derivative, required with response.
!> \param volume_response Cell-volume derivative, required with response.
!> \param pw_grid Optional native full-box grid for real density and energy consumers. Omit for continuous coefficients.
!> \param cell_response Optional cell derivative of reciprocal alias jumps, zero if omitted.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_paw_density_fourier(kind, projected, gvecs, position, volume, density, response, &
                                     projected_response, gvecs_response, position_response, volume_response, pw_grid, cell_response)
      TYPE(paw_kind_type), INTENT(IN)                    :: kind
      COMPLEX(dp), INTENT(IN)                            :: projected(:, :, :)
      REAL(dp), INTENT(IN)                               :: gvecs(:, :), position(3), volume
      COMPLEX(dp), INTENT(OUT)                           :: density(:, :)
      COMPLEX(dp), INTENT(OUT), OPTIONAL                 :: response(:, :)
      COMPLEX(dp), INTENT(IN), OPTIONAL                  :: projected_response(:, :, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: gvecs_response(:, :), &
                                                            position_response(3), volume_response
      TYPE(pw_grid_type), INTENT(IN), OPTIONAL           :: pw_grid
      REAL(dp), INTENT(IN), OPTIONAL                     :: cell_response(3, 3)

      COMPLEX(dp), ALLOCATABLE                           :: kernel(:, :, :), kernel_response(:, :, :)
      INTEGER                                            :: axis, first, gh, gl, grid_shape(3), l, &
                                                            lm, nh, nr, ns, rh, rl, spin
      LOGICAL                                            :: valid_shape
      REAL(dp)                                           :: cell_h(3, 3)
      REAL(dp), ALLOCATABLE                              :: radial(:, :, :), radial_response(:, :, :)

      cell_h = 0.0_dp
      grid_shape = 0
      IF (PRESENT(pw_grid)) THEN
         IF (pw_grid%grid_span == FULLSPACE .AND. .NOT. pw_grid%spherical) THEN
            grid_shape = pw_grid%npts
            DO axis = 1, 3
               cell_h(:, axis) = pw_grid%dh(:, axis)*REAL(grid_shape(axis), dp)
            END DO
         END IF
      END IF
      CPASSERT(kind%active)
      ns = SIZE(projected, 3)
      CPASSERT(ns == 1 .OR. ns == 2)
      CPASSERT(ALL(SHAPE(projected(:, :, 1)) == kind%onecenter%nproj))
      CPASSERT(SIZE(gvecs, 1) == 3)
      CPASSERT(ALL(SHAPE(density) == [SIZE(gvecs, 2), ns]))
      CPASSERT(ALL(ieee_is_finite(ABS(projected))))
      IF (PRESENT(response)) THEN
         CPASSERT(PRESENT(projected_response))
         CPASSERT(PRESENT(gvecs_response))
         CPASSERT(PRESENT(position_response))
         CPASSERT(PRESENT(volume_response))
         valid_shape = ALL(SHAPE(response) == SHAPE(density)) .AND. &
                       ALL(SHAPE(projected_response) == SHAPE(projected)) .AND. &
                       ALL(SHAPE(gvecs_response) == SHAPE(gvecs)) .AND. &
                       ALL(ieee_is_finite(ABS(projected_response)))
         CPASSERT(valid_shape)
         response = CMPLX(0.0_dp, 0.0_dp, dp)
      END IF
      ASSOCIATE (grid => kind%xcgrid, onecenter => kind%onecenter)
         CPASSERT(grid%npair == onecenter%npair .AND. grid%lmax == onecenter%lmax)
         nr = grid%nr
         nh = SIZE(onecenter%gaunt, 3)
         first = 1
         IF (grid%r(1) == 0.0_dp) first = 2
         CALL paw_fourier_radial(kind, projected, radial, .TRUE.)
         IF (PRESENT(response)) CALL paw_fourier_radial(kind, projected_response, radial_response, .FALSE.)
         density = CMPLX(0.0_dp, 0.0_dp, dp)
         ! Bound the temporary kernel independently of grid size and radial mesh length.
         DO gl = 1, SIZE(gvecs, 2), 128
            gh = MIN(SIZE(gvecs, 2), gl + 127)
            DO rl = first, nr, 64
               rh = MIN(nr, rl + 63)
               ALLOCATE (kernel(gh - gl + 1, rh - rl + 1, 1))
               IF (PRESENT(response)) ALLOCATE (kernel_response, MOLD=kernel)
               DO lm = 1, nh
                  l = INT(SQRT(REAL(lm - 1, dp)))
                  IF (PRESENT(response)) THEN
                     CALL upf_inverse_coulomb_fourier_kernel(grid%r(rl:rh), [l], [lm], gvecs(:, gl:gh), &
                                                             position, volume, kernel, kernel_response, &
                                                             0.0_dp*grid%r(rl:rh), gvecs_response(:, gl:gh), &
                                                             position_response, volume_response, &
                                                             cell=cell_h, grid_shape=grid_shape, cell_response=cell_response)
                  ELSE
                     CALL upf_inverse_coulomb_fourier_kernel(grid%r(rl:rh), [l], [lm], gvecs(:, gl:gh), &
                                                             position, volume, kernel, cell=cell_h, grid_shape=grid_shape)
                  END IF
                  DO spin = 1, ns
                     density(gl:gh, spin) = density(gl:gh, spin) + &
                                            MATMUL(kernel(:, :, 1), grid%wr(rl:rh)*radial(rl:rh, lm, spin))
                     IF (PRESENT(response)) THEN
                        response(gl:gh, spin) = response(gl:gh, spin) + &
                                                MATMUL(kernel(:, :, 1), grid%wr(rl:rh)*radial_response(rl:rh, lm, spin)) + &
                                                MATMUL(kernel_response(:, :, 1), grid%wr(rl:rh)*radial(rl:rh, lm, spin))
                     END IF
                  END DO
               END DO
               DEALLOCATE (kernel)
               IF (PRESENT(response)) DEALLOCATE (kernel_response)
            END DO
         END DO
      END ASSOCIATE
   END SUBROUTINE qs_upf_inverse_paw_density_fourier

! **************************************************************************************************
!> \brief Original PAW valence multipoles, with an optional frozen core shared equally over spins.
!> \param kind Original radial products and magnetic-projector Gaunt tensor.
!> \param projected Spin projector occupations, or their derivative when core is false.
!> \param radial Overwritten real density multipoles on the original atomic mesh.
!> \param core Include the fixed original AE core once in the spin sum.
! **************************************************************************************************
   SUBROUTINE paw_fourier_radial(kind, projected, radial, core)
      TYPE(paw_kind_type), INTENT(IN)                    :: kind
      COMPLEX(dp), INTENT(IN)                            :: projected(:, :, :)
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: radial(:, :, :)
      LOGICAL, INTENT(IN)                                :: core

      INTEGER                                            :: l, lm, ns, p, pair, q, spin
      REAL(dp), ALLOCATABLE                              :: coefficient(:, :)

      ASSOCIATE (grid => kind%xcgrid, onecenter => kind%onecenter)
         ns = SIZE(projected, 3)
         ALLOCATE (radial(grid%nr, SIZE(onecenter%gaunt, 3), ns), SOURCE=0.0_dp)
         ALLOCATE (coefficient(grid%npair, ns))
         DO lm = 1, SIZE(onecenter%gaunt, 3)
            l = INT(SQRT(REAL(lm - 1, dp)))
            coefficient = 0.0_dp
            DO q = 1, onecenter%nproj
               DO p = 1, onecenter%nproj
                  pair = onecenter%pair_index(p, q)
                  DO spin = 1, ns
                     coefficient(pair, spin) = coefficient(pair, spin) + &
                                               onecenter%gaunt(p, q, lm)*REAL(projected(p, q, spin), dp)
                  END DO
               END DO
            END DO
            radial(:grid%nval, lm, :) = MATMUL(grid%radial(:, :, l, 1) - grid%radial(:, :, l, 2), coefficient)
         END DO
         IF (core) THEN
            DO spin = 1, ns
               radial(:, 1, spin) = radial(:, 1, spin) + SQRT(fourpi)*grid%core(:, 1)/REAL(ns, dp)
            END DO
         END IF
      END ASSOCIATE
   END SUBROUTINE paw_fourier_radial

! **************************************************************************************************
!> \brief Transpose of the original PAW Fourier valence correction and frozen AE core.
!>        The scalar field is spin independent. The projector operator has no spin or k weight.
!>        Pair its valence contraction with the frozen-core energy once, not once per spin.
!> \param kind Original PAW radial products and magnetic-projector Gaunt tensor.
!> \param gvecs Local physical reciprocal vectors.
!> \param position Home atomic center for the periodic Fourier phase.
!> \param volume Cell volume used in the density Fourier coefficients.
!> \param potential Fourier covector, with volume and multiplicity already included once.
!> \param OPERATOR Overwritten real symmetric valence projector operator.
!> \param core_energy Coupling to the full frozen AE core, without an ionic counterterm.
!> \param comm Optional communicator distributing Fourier rows. Results are replicated.
!> \param pw_grid Optional native full-box grid for real density and energy consumers. Omit for continuous coefficients.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_paw_density_adjoint(kind, gvecs, position, volume, potential, OPERATOR, core_energy, comm, pw_grid)
      TYPE(paw_kind_type), INTENT(IN)                    :: kind
      REAL(dp), INTENT(IN)                               :: gvecs(:, :), position(3), volume
      COMPLEX(dp), INTENT(IN)                            :: potential(:)
      REAL(dp), INTENT(OUT)                              :: OPERATOR(:, :), core_energy
      TYPE(mp_comm_type), INTENT(IN), OPTIONAL           :: comm
      TYPE(pw_grid_type), INTENT(IN), OPTIONAL           :: pw_grid

      COMPLEX(dp), ALLOCATABLE                           :: kernel(:, :, :)
      INTEGER                                            :: axis, first, gh, gl, grid_shape(3), ir, &
                                                            l, lm, nr, p, pair, q, rh, rl
      LOGICAL                                            :: valid_shape
      REAL(dp)                                           :: cell_h(3, 3)
      REAL(dp), ALLOCATABLE                              :: harmonics(:, :), radial(:, :)

      cell_h = 0.0_dp
      grid_shape = 0
      IF (PRESENT(pw_grid)) THEN
         IF (pw_grid%grid_span == FULLSPACE .AND. .NOT. pw_grid%spherical) THEN
            grid_shape = pw_grid%npts
            DO axis = 1, 3
               cell_h(:, axis) = pw_grid%dh(:, axis)*REAL(grid_shape(axis), dp)
            END DO
         END IF
      END IF
      CPASSERT(kind%active)
      CPASSERT(SIZE(gvecs, 1) == 3)
      CPASSERT(SIZE(gvecs, 2) == SIZE(potential))
      CPASSERT(ALL(ieee_is_finite(ABS(potential))))
      valid_shape = ALL(SHAPE(OPERATOR) == kind%onecenter%nproj)
      CPASSERT(valid_shape)
      ASSOCIATE (grid => kind%xcgrid, onecenter => kind%onecenter)
         CPASSERT(grid%npair == onecenter%npair .AND. grid%lmax == onecenter%lmax)
         nr = grid%nr
         first = 1
         IF (grid%r(1) == 0.0_dp) first = 2
         ALLOCATE (harmonics(nr, SIZE(onecenter%gaunt, 3)), SOURCE=0.0_dp)
         ALLOCATE (radial(grid%npair, SIZE(onecenter%gaunt, 3)), SOURCE=0.0_dp)
         DO gl = 1, SIZE(gvecs, 2), 128
            gh = MIN(SIZE(gvecs, 2), gl + 127)
            DO rl = first, nr, 64
               rh = MIN(nr, rl + 63)
               ALLOCATE (kernel(gh - gl + 1, rh - rl + 1, 1))
               DO lm = 1, SIZE(harmonics, 2)
                  l = INT(SQRT(REAL(lm - 1, dp)))
                  CALL upf_inverse_coulomb_fourier_kernel(grid%r(rl:rh), [l], [lm], gvecs(:, gl:gh), &
                                                          position, volume, kernel, cell=cell_h, grid_shape=grid_shape)
                  DO ir = rl, rh
                     harmonics(ir, lm) = harmonics(ir, lm) + &
                                         REAL(SUM(CONJG(kernel(:, ir - rl + 1, 1))*potential(gl:gh)), dp)
                  END DO
               END DO
               DEALLOCATE (kernel)
            END DO
         END DO
         core_energy = SQRT(fourpi)*SUM(grid%wr*grid%core(:, 1)*harmonics(:, 1))
         DO lm = 1, SIZE(harmonics, 2)
            l = INT(SQRT(REAL(lm - 1, dp)))
            DO pair = 1, grid%npair
               radial(pair, lm) = SUM(grid%wr(:grid%nval)*harmonics(:grid%nval, lm)* &
                                      (grid%radial(:, pair, l, 1) - grid%radial(:, pair, l, 2)))
            END DO
         END DO
         OPERATOR = 0.0_dp
         DO q = 1, onecenter%nproj
            DO p = 1, onecenter%nproj
               pair = onecenter%pair_index(p, q)
               OPERATOR(p, q) = SUM(onecenter%gaunt(p, q, :)*radial(pair, :))
            END DO
         END DO
      END ASSOCIATE
      IF (PRESENT(comm)) THEN
         CALL comm%sum(OPERATOR)
         CALL comm%sum(core_energy)
      END IF
   END SUBROUTINE qs_upf_inverse_paw_density_adjoint

! **************************************************************************************************
!> \brief Couple a scalar field to the PAW AE-minus-compensated-PS valence correction and AE core.
!>        The returned projector operator has no spin or k-point weight. The frozen-core energy
!>        is counted once, independently of the valence density. The supplied field contains no
!>        radial or angular integration weights. A field response gives the explicit derivative
!>        on the same atom-centered rule. The inverse coefficient response remains with the caller.
!> \param kind Original PAW radial products, magnetic-projector mapping and atom-centered rule.
!> \param potential Scalar field, indexed by angular point and radial node.
!> \param OPERATOR New real symmetric valence projector operator.
!> \param core_energy Frozen AE-core coupling, without subtracting the XC model core.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_paw_local_operator(kind, potential, OPERATOR, core_energy)
      TYPE(paw_kind_type), INTENT(IN)                    :: kind
      REAL(dp), INTENT(IN)                               :: potential(:, :)
      REAL(dp), INTENT(OUT)                              :: OPERATOR(:, :), core_energy

      INTEGER                                            :: ir, l, lm, p, pair, q
      REAL(dp), ALLOCATABLE                              :: harmonics(:), radial(:, :), weighted(:)

      CPASSERT(kind%active)
      ASSOCIATE (grid => kind%xcgrid, onecenter => kind%onecenter)
         CPASSERT(grid%npair == onecenter%npair .AND. grid%lmax == onecenter%lmax)
         CPASSERT(ALL(SHAPE(potential) == [grid%na, grid%nr]))
         CPASSERT(ALL(SHAPE(OPERATOR) == [onecenter%nproj, onecenter%nproj]))
         CPASSERT(ALL(ieee_is_finite(potential)))
         ALLOCATE (harmonics(SIZE(onecenter%gaunt, 3)), weighted(grid%na))
         ALLOCATE (radial(grid%npair, SIZE(harmonics)), SOURCE=0.0_dp)
         core_energy = 0.0_dp
         DO ir = 1, grid%nr
            IF (grid%wr(ir) == 0.0_dp) CYCLE
            weighted = grid%wa*potential(:, ir)
            core_energy = core_energy + grid%wr(ir)*grid%core(ir, 1)*SUM(weighted)
            IF (ir > grid%nval) CYCLE
            harmonics = MATMUL(TRANSPOSE(grid%harmonics), weighted)
            DO lm = 1, SIZE(harmonics)
               l = INT(SQRT(REAL(lm - 1, dp)))
               radial(:, lm) = radial(:, lm) + grid%wr(ir)*harmonics(lm)* &
                               (grid%radial(ir, :, l, 1) - grid%radial(ir, :, l, 2))
            END DO
         END DO
         OPERATOR = 0.0_dp
         DO q = 1, onecenter%nproj
            DO p = 1, onecenter%nproj
               pair = onecenter%pair_index(p, q)
               OPERATOR(p, q) = SUM(onecenter%gaunt(p, q, :)*radial(pair, :))
            END DO
         END DO
      END ASSOCIATE
   END SUBROUTINE qs_upf_inverse_paw_local_operator

! **************************************************************************************************
!> \brief Closed-core exchange with the same physical interactions as the native valence exchange.
!>        Sum each HF contribution once, independently of its thread replicas. The valence-valence
!>        exchange is already evaluated through the reconstructed RI transitions.
!> \param pot Original scalar PAW data, including the complete signed core orbitals.
!> \param projectors Original magnetic-projector ordering.
!> \param kind Original-data caches, reused without changing their interactions.
!> \param exchange Native HF contributions and thread replicas, with parameters in atomic units.
!> \param nspins Number of physical valence spins, used for the native RI activity threshold.
!> \param energy Frozen core-core exchange constant, summed over the selected HF contributions.
!> \param potential Spin-independent core-valence projector operator, without an extra factor 1/2.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_paw_core_prepare(pot, projectors, kind, exchange, nspins, energy, potential)
      TYPE(atom_upfpot_type), INTENT(IN)                 :: pot
      TYPE(upf_projector_type), INTENT(IN)               :: projectors
      TYPE(paw_kind_type), INTENT(INOUT)                 :: kind
      TYPE(hfx_type), INTENT(IN)                         :: exchange(:, :)
      INTEGER, INTENT(IN)                                :: nspins
      REAL(dp), INTENT(OUT)                              :: energy, potential(:, :)

      INTEGER                                            :: cache, contribution, interaction(2), term
      REAL(dp)                                           :: coefficient(2), fraction, value
      REAL(dp), ALLOCATABLE                              :: matrix(:, :, :), zero_density(:, :, :)

      CPASSERT(pot%is_paw .AND. kind%active)
      CPASSERT(nspins == 1 .OR. nspins == 2)
      CPASSERT(SIZE(exchange, 2) > 0)
      CPASSERT(ALL(SHAPE(potential) == [projectors%nproj, projectors%nproj]))
      energy = 0.0_dp
      potential = 0.0_dp
      ALLOCATE (zero_density(projectors%nproj, projectors%nproj, 1), SOURCE=0.0_dp)
      ALLOCATE (matrix(projectors%nproj, projectors%nproj, 1))
      DO contribution = 1, SIZE(exchange, 1)
         ASSOCIATE (hf => exchange(contribution, 1))
            fraction = hf%general_parameter%fraction
            CPASSERT(ieee_is_finite(fraction))
            ! Match the spin prefactor used by hfx_ri_update_ks and its Fock contraction.
            IF (.NOT. MERGE(0.5_dp*fraction, fraction, nspins == 1) > EPSILON(0.0_dp)) CYCLE
            coefficient = 0.0_dp
            interaction = 0
            SELECT CASE (hf%potential_parameter%potential_type)
            CASE (do_potential_coulomb)
               coefficient(1) = fraction
            CASE (do_potential_short)
               coefficient(1) = fraction
               interaction(1) = 1
            CASE (do_potential_long)
               coefficient = [fraction, -fraction]
               interaction = [0, 1]
            CASE (do_potential_mix_cl)
               coefficient = fraction*[hf%potential_parameter%scale_coulomb + &
                                       hf%potential_parameter%scale_longrange, &
                                       -hf%potential_parameter%scale_longrange]
               interaction = [0, 1]
            CASE (do_potential_truncated)
               coefficient(1) = fraction
               interaction(1) = 2
            CASE (do_potential_mix_cl_trunc)
               coefficient = fraction*[hf%potential_parameter%scale_coulomb + &
                                       hf%potential_parameter%scale_longrange, &
                                       -hf%potential_parameter%scale_longrange]
               interaction = [2, 1]
            CASE DEFAULT
               CPABORT("Inverse PAW core exchange for this interaction is not implemented yet")
            END SELECT
            CPASSERT(ALL(ieee_is_finite(coefficient)))
            DO term = 1, 2
               IF (coefficient(term) == 0.0_dp) CYCLE
               SELECT CASE (interaction(term))
               CASE (0)
                  IF (.NOT. kind%onecenter%core_exchange_initialized) THEN
                     CALL paw_onecenter_core_exchange_init(pot, projectors, kind%onecenter)
                  END IF
                  IF (.NOT. kind%onecenter%core_exchange_available) THEN
                     CPABORT(TRIM(kind%onecenter%core_exchange_reason))
                  END IF
                  cache = 0
               CASE (1)
                  CALL paw_onecenter_core_exchange_prepare(pot, projectors, kind%onecenter, cache_index=cache, &
                                                           omega=hf%potential_parameter%omega)
               CASE (2)
                  CALL paw_onecenter_core_exchange_prepare(pot, projectors, kind%onecenter, cache_index=cache, &
                                                           cutoff_radius=hf%potential_parameter%cutoff_radius)
               END SELECT
               IF (cache > 0) THEN
                  IF (.NOT. kind%onecenter%interaction_core_exchange(cache)%available) THEN
                     CPABORT(TRIM(kind%onecenter%interaction_core_exchange(cache)%reason))
                  END IF
               END IF
               CALL paw_onecenter_core_exchange(kind%onecenter, zero_density, coefficient(term), cache, value, matrix)
               energy = energy + value
               potential = potential + matrix(:, :, 1)
            END DO
         END ASSOCIATE
      END DO
   END SUBROUTINE qs_upf_inverse_paw_core_prepare

! **************************************************************************************************
!> \brief PAW valence reconstruction relative to the compensated soft density, including its gradient.
!>        Adds AE-minus-PS partial-wave products minus the original compensation Q. No core is added.
!>        Same-center kernels are periodic without Bloch phases; the inverse coefficients carry k.
!> \param frame Original inverse frame and finite image list.
!> \param reconstruction Kind-indexed original AE/PS caches. Inactive kinds have zero channels.
!> \param augmentation Original compensation kernels in the same magnetic-projector order.
!> \param points Cartesian evaluation points without quadrature or k weights.
!> \param matrix Hermitian AO occupation matrix for one spin and k point.
!> \param density Reconstructed valence correction, overwritten.
!> \param gradient Cartesian gradient of the same correction, overwritten.
!> \param density_response Optional fixed-matrix geometry response, requiring all response inputs.
!> \param gradient_response Optional response of the Cartesian gradient.
!> \param point_response Cartesian motion of the points.
!> \param atom_response Cartesian motion of the original atoms.
!> \param cell_response Lattice response, including translated image centers.
!> \param coefficients_response Complete coupled inverse response dW.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_paw_collocate(frame, reconstruction, augmentation, points, matrix, density, gradient, &
                                           density_response, gradient_response, point_response, atom_response, &
                                           cell_response, coefficients_response)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      TYPE(upf_inverse_kinetic_type), INTENT(IN)         :: reconstruction(:)
      TYPE(qs_upf_inverse_augmentation_type), INTENT(IN) :: augmentation
      REAL(dp), INTENT(IN)                               :: points(:, :)
      COMPLEX(dp), INTENT(IN)                            :: matrix(:, :)
      REAL(dp), INTENT(OUT)                              :: density(:), gradient(:, :)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: density_response(:), &
                                                            gradient_response(:, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: point_response(:, :), &
                                                            atom_response(:, :), &
                                                            cell_response(3, 3)
      COMPLEX(dp), INTENT(IN), OPTIONAL                  :: coefficients_response(:, :)

      COMPLEX(dp), ALLOCATABLE                           :: projected(:, :), projected_response(:, :)
      INTEGER                                            :: atom, axis, center, first, ikind, image, &
                                                            last, n, nao, np, point, response_axis
      LOGICAL                                            :: derive
      REAL(dp)                                           :: motion(3), position(3), &
                                                            position_response(3), radius
      REAL(dp), ALLOCATABLE                              :: kernel(:, :), kernel_gradient(:, :, :), &
                                                            kernel_hessian(:, :, :, :)

      nao = SIZE(frame%coefficients, 2)
      np = SIZE(points, 2)
      CPASSERT(frame%atom_local .OR. frame%inverse%ready)
      CPASSERT(SIZE(points, 1) == 3 .AND. ALL(SHAPE(matrix) == [nao, nao]))
      CPASSERT(SIZE(density) == np .AND. ALL(SHAPE(gradient) == [np, 3]))
      CPASSERT(ALLOCATED(augmentation%kinds))
      CPASSERT(SIZE(reconstruction) == SIZE(augmentation%kinds))
      derive = PRESENT(density_response)
      CPASSERT(derive .EQV. PRESENT(gradient_response))
      CPASSERT(derive .EQV. PRESENT(point_response))
      CPASSERT(derive .EQV. PRESENT(atom_response))
      CPASSERT(derive .EQV. PRESENT(cell_response))
      CPASSERT(derive .EQV. PRESENT(coefficients_response))
      density = 0.0_dp
      gradient = 0.0_dp
      IF (derive) THEN
         CPASSERT(SIZE(density_response) == np)
         CPASSERT(ALL(SHAPE(gradient_response) == [np, 3]))
         CPASSERT(ALL(SHAPE(point_response) == SHAPE(points)))
         CPASSERT(ALL(SHAPE(atom_response) == SHAPE(frame%positions)))
         CPASSERT(ALL(SHAPE(coefficients_response) == SHAPE(frame%coefficients)))
         density_response = 0.0_dp
         gradient_response = 0.0_dp
      END IF
      DO center = 1, SIZE(frame%center_atoms)
         ikind = frame%center_kinds(center)
         CPASSERT(ikind > 0 .AND. ikind <= SIZE(reconstruction))
         IF (reconstruction(ikind)%delta%nproj == 0) CYCLE
         atom = frame%center_atoms(center)
         first = frame%projector_first(center)
         last = frame%projector_last(center)
         n = last - first + 1
         CPASSERT(n == reconstruction(ikind)%delta%nproj)
         CPASSERT(n == augmentation%kinds(ikind)%nproj)
         radius = MAX(reconstruction(ikind)%delta%radius, augmentation%kinds(ikind)%realspace_radius)
         ALLOCATE (projected(n, n), kernel(n, n), kernel_gradient(n, n, 3))
         projected(:, :) = MATMUL(frame%coefficients(first:last, :), &
                                  MATMUL(matrix, CONJG(TRANSPOSE(frame%coefficients(first:last, :)))))
         IF (derive) THEN
            ALLOCATE (projected_response(n, n), kernel_hessian(n, n, 3, 3))
            projected_response(:, :) = MATMUL(coefficients_response(first:last, :), &
                                              MATMUL(matrix, CONJG(TRANSPOSE(frame%coefficients(first:last, :)))))
            projected_response(:, :) = projected_response + CONJG(TRANSPOSE(projected_response))
         END IF
         DO image = 1, SIZE(frame%images, 2)
            position = frame%positions(:, atom) + MATMUL(frame%cell, REAL(frame%images(:, image), dp))
            IF (derive) position_response = atom_response(:, atom) + &
                                            MATMUL(cell_response, REAL(frame%images(:, image), dp))
            DO point = 1, np
               IF (SUM((points(:, point) - position)**2) >= radius**2) CYCLE
               IF (derive) THEN
                  CALL reconstruction_fields(reconstruction(ikind), augmentation%kinds(ikind), &
                                             points(:, point) - position, kernel, kernel_gradient, kernel_hessian)
                  motion = point_response(:, point) - position_response
                  density_response(point) = density_response(point) + REAL(SUM(projected_response*kernel), dp)
                  DO axis = 1, 3
                     density_response(point) = density_response(point) + REAL(SUM(projected*kernel_gradient(:, :, axis)), dp)* &
                                               motion(axis)
                     gradient_response(point, axis) = gradient_response(point, axis) + &
                                                      REAL(SUM(projected_response*kernel_gradient(:, :, axis)), dp)
                     DO response_axis = 1, 3
                        gradient_response(point, axis) = gradient_response(point, axis) + &
                                                         REAL(SUM(projected*kernel_hessian(:, :, axis, response_axis)), dp)* &
                                                         motion(response_axis)
                     END DO
                  END DO
               ELSE
                  CALL reconstruction_fields(reconstruction(ikind), augmentation%kinds(ikind), &
                                             points(:, point) - position, kernel, kernel_gradient)
               END IF
               density(point) = density(point) + REAL(SUM(projected*kernel), dp)
               DO axis = 1, 3
                  gradient(point, axis) = gradient(point, axis) + REAL(SUM(projected*kernel_gradient(:, :, axis)), dp)
               END DO
            END DO
         END DO
         DEALLOCATE (projected, kernel, kernel_gradient)
         IF (derive) DEALLOCATE (projected_response, kernel_hessian)
      END DO
   END SUBROUTINE qs_upf_inverse_paw_collocate

! **************************************************************************************************
!> \brief Exact AO transpose of the same PAW density and Cartesian-gradient reconstruction.
!> \param frame ...
!> \param reconstruction ...
!> \param augmentation ...
!> \param points ...
!> \param weights Integration measure, applied once. Weighted model derivatives use unit weights.
!> \param vdensity Potential conjugate to the reconstructed density.
!> \param vgradient Potential conjugate to its Cartesian gradient.
!> \param ao_operator Hermitian AO operator without spin or k weights.
!> \param comm Optional communicator distributing point rows. Replicated points omit it.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_paw_operator(frame, reconstruction, augmentation, points, weights, vdensity, vgradient, &
                                          ao_operator, comm)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      TYPE(upf_inverse_kinetic_type), INTENT(IN)         :: reconstruction(:)
      TYPE(qs_upf_inverse_augmentation_type), INTENT(IN) :: augmentation
      REAL(dp), INTENT(IN)                               :: points(:, :), weights(:), vdensity(:), &
                                                            vgradient(:, :)
      COMPLEX(dp), INTENT(OUT)                           :: ao_operator(:, :)
      TYPE(mp_comm_type), INTENT(IN), OPTIONAL           :: comm

      INTEGER                                            :: atom, axis, center, first, ikind, image, &
                                                            last, n, nao, np, point
      REAL(dp)                                           :: position(3), radius
      REAL(dp), ALLOCATABLE                              :: kernel(:, :), kernel_gradient(:, :, :), &
                                                            OPERATOR(:, :)

      nao = SIZE(frame%coefficients, 2)
      np = SIZE(points, 2)
      CPASSERT(frame%atom_local .OR. frame%inverse%ready)
      CPASSERT(SIZE(points, 1) == 3 .AND. SIZE(weights) == np .AND. SIZE(vdensity) == np)
      CPASSERT(ALL(SHAPE(vgradient) == [np, 3]) .AND. ALL(SHAPE(ao_operator) == [nao, nao]))
      CPASSERT(ALLOCATED(augmentation%kinds))
      CPASSERT(SIZE(reconstruction) == SIZE(augmentation%kinds))
      ao_operator = 0.0_dp
      DO center = 1, SIZE(frame%center_atoms)
         ikind = frame%center_kinds(center)
         CPASSERT(ikind > 0 .AND. ikind <= SIZE(reconstruction))
         IF (reconstruction(ikind)%delta%nproj == 0) CYCLE
         atom = frame%center_atoms(center)
         first = frame%projector_first(center)
         last = frame%projector_last(center)
         n = last - first + 1
         CPASSERT(n == reconstruction(ikind)%delta%nproj)
         CPASSERT(n == augmentation%kinds(ikind)%nproj)
         radius = MAX(reconstruction(ikind)%delta%radius, augmentation%kinds(ikind)%realspace_radius)
         ALLOCATE (kernel(n, n), kernel_gradient(n, n, 3), OPERATOR(n, n))
         OPERATOR(:, :) = 0.0_dp
         DO image = 1, SIZE(frame%images, 2)
            position = frame%positions(:, atom) + MATMUL(frame%cell, REAL(frame%images(:, image), dp))
            DO point = 1, np
               IF (SUM((points(:, point) - position)**2) >= radius**2) CYCLE
               CALL reconstruction_fields(reconstruction(ikind), augmentation%kinds(ikind), &
                                          points(:, point) - position, kernel, kernel_gradient)
               OPERATOR(:, :) = OPERATOR + weights(point)*vdensity(point)*kernel
               DO axis = 1, 3
                  OPERATOR(:, :) = OPERATOR + weights(point)*vgradient(point, axis)*kernel_gradient(:, :, axis)
               END DO
            END DO
         END DO
         ao_operator = ao_operator + MATMUL(CONJG(TRANSPOSE(frame%coefficients(first:last, :))), &
                                            MATMUL(OPERATOR, frame%coefficients(first:last, :)))
         DEALLOCATE (kernel, kernel_gradient, OPERATOR)
      END DO
      IF (PRESENT(comm)) CALL comm%sum(ao_operator)
   END SUBROUTINE qs_upf_inverse_paw_operator

! **************************************************************************************************
!> \brief Original AE-minus-PS products minus original compensation, with shared spatial derivatives.
!> \param reconstruction ...
!> \param augmentation ...
!> \param displacement ...
!> \param kernel ...
!> \param gradient ...
!> \param hessian ...
! **************************************************************************************************
   SUBROUTINE reconstruction_fields(reconstruction, augmentation, displacement, kernel, gradient, hessian)
      TYPE(upf_inverse_kinetic_type), INTENT(IN)         :: reconstruction
      TYPE(uspp_augmentation_type), INTENT(IN)           :: augmentation
      REAL(dp), INTENT(IN)                               :: displacement(3)
      REAL(dp), INTENT(OUT)                              :: kernel(:, :), gradient(:, :, :)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: hessian(:, :, :, :)

      INTEGER                                            :: n
      REAL(dp), ALLOCATABLE                              :: q(:, :), qgradient(:, :, :), &
                                                            qhessian(:, :, :, :)

      n = reconstruction%delta%nproj
      ALLOCATE (q(n, n), qgradient(n, n, 3))
      CALL upf_inverse_density_evaluate(reconstruction, displacement, kernel, gradient, hessian)
      IF (PRESENT(hessian)) THEN
         ALLOCATE (qhessian(n, n, 3, 3))
         CALL uspp_augmentation_realspace(augmentation, displacement, q, qgradient, qhessian)
         hessian = hessian - qhessian
      ELSE
         CALL uspp_augmentation_realspace(augmentation, displacement, q, qgradient)
      END IF
      kernel = kernel - q
      gradient = gradient - qgradient
   END SUBROUTINE reconstruction_fields

! **************************************************************************************************
!> \brief Prepare original PAW radial kernels without changing the supplied operators or waves.
!> \param pot Original scalar PAW dataset.
!> \param projectors Matching radial and magnetic projector mapping.
!> \param nang Requested atomic angular quadrature.
!> \param cache Owned per-kind PAW data.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_paw_init(pot, projectors, nang, cache)
      TYPE(atom_upfpot_type), INTENT(IN)                 :: pot
      TYPE(upf_projector_type), INTENT(IN)               :: projectors
      INTEGER, INTENT(IN)                                :: nang
      TYPE(paw_kind_type), INTENT(OUT)                   :: cache

      CPASSERT(pot%is_paw .AND. .NOT. pot%has_so)
      CPASSERT(ALLOCATED(pot%paw))
      CALL paw_onecenter_init(pot, projectors, cache%onecenter)
      CALL paw_xc_init(pot, nang, cache%xcgrid)
      cache%core_energy = pot%paw%core_energy
      cache%active = .TRUE.
   END SUBROUTINE qs_upf_inverse_paw_init

! **************************************************************************************************
!> \brief Evaluate PAW Hartree, XC and core energy and their Gaussian Hamiltonian derivative.
!> \param kinds Prepared per-kind PAW kernels. Inactive NC and USPP kinds contribute zero.
!> \param frames Native inverse frames for the complete global k mesh, identically ordered on all ranks.
!> \param density Hermitian Gaussian occupation matrices, AO by AO by global k point by spin.
!> \param kweights Normalized global k weights. They enter the projected density once.
!> \param functional Activated XC_FUNCTIONAL section.
!> \param density_cutoff XC density cutoff.
!> \param gradient_cutoff XC gradient cutoff.
!> \param energies One-center Hartree, XC and frozen-core constant, counted once per atom.
!> \param hamiltonian Unweighted Gaussian Hamiltonian correction for each global k point and spin.
!> \param potential Owned per-center projector potentials for the corresponding geometry response.
!> \param semilocal_xc Include atomic semilocal XC. A separate composite XC consumer may disable it.
!> \param core_energy Optional kind-indexed core-core exchange constants for the native HF interactions.
!> \param core_potential Optional kind-indexed core-valence operators, with one spin-independent matrix.
!> \param samples Optional physical XC arguments, retained on the existing atomic owners.
!> \param orbital_factors Optional original factors certified for the installed density.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_paw_evaluate(kinds, frames, density, kweights, functional, density_cutoff, &
                                          gradient_cutoff, energies, hamiltonian, potential, semilocal_xc, &
                                          core_energy, core_potential, samples, orbital_factors)
      TYPE(paw_kind_type), INTENT(IN)                    :: kinds(:)
      TYPE(qs_upf_inverse_grid_type), INTENT(IN)         :: frames(:)
      COMPLEX(dp), INTENT(IN)                            :: density(:, :, :, :)
      REAL(dp), INTENT(IN)                               :: kweights(:)
      TYPE(section_vals_type), POINTER                   :: functional
      REAL(dp), INTENT(IN)                               :: density_cutoff, gradient_cutoff
      REAL(dp), INTENT(OUT)                              :: energies(3)
      COMPLEX(dp), INTENT(OUT)                           :: hamiltonian(:, :, :, :)
      TYPE(paw_potential_type), ALLOCATABLE, INTENT(OUT) :: potential(:)
      LOGICAL, INTENT(IN), OPTIONAL                      :: semilocal_xc
      REAL(dp), INTENT(IN), OPTIONAL                     :: core_energy(:)
      TYPE(paw_potential_type), INTENT(IN), OPTIONAL     :: core_potential(:)
      TYPE(skala_density_constraints_type), &
         INTENT(INOUT), OPTIONAL, POINTER                :: samples
      TYPE(qs_upf_inverse_orbitals_type), INTENT(IN), &
         OPTIONAL                                        :: orbital_factors

      COMPLEX(dp), ALLOCATABLE                           :: atomic_operator(:, :)
      INTEGER                                            :: center, first, ik, ikind, last, n, nk, &
                                                            ns, p, q, spin
      LOGICAL                                            :: use_semilocal
      REAL(dp)                                           :: eh, ex
      REAL(dp), ALLOCATABLE                              :: dm(:, :, :), vh(:, :), vxc(:, :, :)
      TYPE(twofold_real_type), ALLOCATABLE               :: density_components(:, :, :), &
                                                            potential_components(:, :, :)

      use_semilocal = .TRUE.
      IF (PRESENT(semilocal_xc)) use_semilocal = semilocal_xc
      nk = SIZE(frames)
      CPASSERT(PRESENT(core_energy) .EQV. PRESENT(core_potential))
      IF (PRESENT(core_energy)) THEN
         CPASSERT(SIZE(core_energy) == SIZE(kinds))
         CPASSERT(SIZE(core_potential) == SIZE(kinds))
         CPASSERT(ALL(ieee_is_finite(core_energy)))
      END IF
      CALL check_paw_density(frames, density, kweights)
      ns = SIZE(density, 4)
      CPASSERT(ALL(SHAPE(hamiltonian) == SHAPE(density)))
      ALLOCATE (potential(SIZE(frames(1)%ao%center_atoms)))
      energies = 0.0_dp
      hamiltonian = 0.0_dp
      DO center = 1, SIZE(potential)
         ikind = frames(1)%ao%center_kinds(center)
         CPASSERT(ikind > 0 .AND. ikind <= SIZE(kinds))
         IF (.NOT. kinds(ikind)%active) CYCLE
         first = frames(1)%ao%projector_first(center)
         last = frames(1)%ao%projector_last(center)
         n = last - first + 1
         CPASSERT(n == kinds(ikind)%onecenter%nproj)
         IF (PRESENT(core_potential)) THEN
            CPASSERT(ALLOCATED(core_potential(ikind)%matrix))
            CPASSERT(ALL(SHAPE(core_potential(ikind)%matrix) == [n, n, 1]))
            CPASSERT(ALL(ieee_is_finite(core_potential(ikind)%matrix)))
         END IF
         ALLOCATE (potential(center)%matrix(n, n, ns), SOURCE=0.0_dp)
         ! The global density and inverse coefficients are replicated. One rank owns
         ! each atomic functional, including the core constant, independently of k groups.
         IF (MODULO(center - 1, frames(1)%ao%point_comm%num_pe) == frames(1)%ao%point_comm%mepos) THEN
            ALLOCATE (dm(n, n, ns), vh(n, n), vxc(n, n, ns))
            CALL project_paw_density(frames, density, kweights, first, last, dm)
            ALLOCATE (potential_components(n, n, ns))
            potential_components = twofold_real_type()
            CALL paw_onecenter_hartree(kinds(ikind)%onecenter, SUM(dm, DIM=3), eh, vh)
            ex = 0.0_dp
            vxc = 0.0_dp
            IF (use_semilocal) THEN
               ALLOCATE (density_components(n, n, ns))
               CALL retained_projector_density(frames, kweights, first, last, dm, density_components, orbital_factors)
               CALL paw_xc_evaluate(kinds(ikind)%xcgrid, kinds(ikind)%onecenter, functional, dm, &
                                    density_cutoff, gradient_cutoff, ex, vxc, &
                                    samples=samples, atom=frames(1)%ao%center_atoms(center), &
                                    density_components=density_components, potential_components=potential_components)
               DEALLOCATE (density_components)
            END IF
            IF (PRESENT(core_energy)) THEN
               ex = ex + core_energy(ikind) + SUM(SUM(dm, DIM=3)*core_potential(ikind)%matrix(:, :, 1))
               DO spin = 1, ns
                  vxc(:, :, spin) = vxc(:, :, spin) + core_potential(ikind)%matrix(:, :, 1)
                  DO q = 1, n
                     DO p = 1, n
                        potential_components(p, q, spin) = twofold_add(potential_components(p, q, spin), &
                                                                   twofold_real_type(core_potential(ikind)%matrix(p, q, 1), 0.0_dp))
                     END DO
                  END DO
               END DO
            END IF
            energies = energies + [eh, ex, kinds(ikind)%core_energy]
            DO spin = 1, ns
               DO q = 1, n
                  DO p = 1, n
               potential_components(p, q, spin) = twofold_add(potential_components(p, q, spin), twofold_real_type(vh(p, q), 0.0_dp))
                  END DO
               END DO
               potential(center)%matrix(:, :, spin) = twofold_value(potential_components(:, :, spin))
            END DO
            DEALLOCATE (dm, vh, vxc)
         END IF
         IF (.NOT. ALLOCATED(potential_components)) ALLOCATE (potential_components(n, n, ns), SOURCE=twofold_real_type())
         CALL frames(1)%ao%point_comm%sum(potential_components%hi)
         CALL frames(1)%ao%point_comm%sum(potential_components%lo)
         potential(center)%matrix = twofold_value(potential_components)
         ALLOCATE (potential(center)%components, SOURCE=potential_components)
         ALLOCATE (atomic_operator(SIZE(density, 1), SIZE(density, 2)))
         DO spin = 1, ns
            DO ik = 1, nk
   CALL upf_inverse_projector_operator(frames(ik)%ao%coefficients(first:last, :), potential_components(:, :, spin), atomic_operator)
               hamiltonian(:, :, ik, spin) = hamiltonian(:, :, ik, spin) + atomic_operator
            END DO
         END DO
         DEALLOCATE (atomic_operator, potential_components)
      END DO
      CALL frames(1)%ao%point_comm%sum(energies)
   END SUBROUTINE qs_upf_inverse_paw_evaluate

! **************************************************************************************************
!> \brief Validate the replicated PAW frame layout and unweighted Hermitian AO density.
!> \param frames Original globally ordered projector frames.
!> \param density AO occupation matrices by k point and spin.
!> \param kweights Original normalized k weights.
! **************************************************************************************************
   SUBROUTINE check_paw_density(frames, density, kweights)
      TYPE(qs_upf_inverse_grid_type), INTENT(IN)         :: frames(:)
      COMPLEX(dp), INTENT(IN)                            :: density(:, :, :, :)
      REAL(dp), INTENT(IN)                               :: kweights(:)

      INTEGER                                            :: ik, nao, nk, nproj, ns, spin
      REAL(dp)                                           :: hermitian_error, hermitian_limit

      nk = SIZE(frames)
      CPASSERT(nk > 0 .AND. SIZE(kweights) == nk)
      nao = SIZE(frames(1)%ao%coefficients, 2)
      nproj = SIZE(frames(1)%ao%coefficients, 1)
      ns = SIZE(density, 4)
      CPASSERT(ns == 1 .OR. ns == 2)
      CPASSERT(ALL(SHAPE(density) == [nao, nao, nk, ns]))
      CPASSERT(ALL(ieee_is_finite(kweights)) .AND. ALL(kweights >= 0.0_dp))
      CPASSERT(ABS(SUM(kweights) - 1.0_dp) < 1.e-12_dp)
      DO ik = 1, nk
         CPASSERT(SIZE(frames(ik)%ao%coefficients, 1) == nproj)
         CPASSERT(SIZE(frames(ik)%ao%coefficients, 2) == nao)
         CPASSERT(ALL(frames(ik)%ao%center_atoms == frames(1)%ao%center_atoms))
         CPASSERT(ALL(frames(ik)%ao%center_kinds == frames(1)%ao%center_kinds))
         IF (.NOT. (ALL(frames(ik)%ao%projector_first == frames(1)%ao%projector_first))) THEN
            CPABORT("Inconsistent inverse PAW frame or density")
         END IF
         CPASSERT(ALL(frames(ik)%ao%projector_last == frames(1)%ao%projector_last))
         DO spin = 1, ns
            CPASSERT(ALL(ieee_is_finite(REAL(density(:, :, ik, spin), dp))))
            CPASSERT(ALL(ieee_is_finite(AIMAG(density(:, :, ik, spin)))))
            hermitian_error = MAXVAL(ABS(density(:, :, ik, spin) - CONJG(TRANSPOSE(density(:, :, ik, spin)))))
            hermitian_limit = 1.e-12_dp*(1.0_dp + MAXVAL(ABS(density(:, :, ik, spin))))
            CPASSERT(hermitian_error <= hermitian_limit)
         END DO
      END DO
   END SUBROUTINE check_paw_density

! **************************************************************************************************
!> \brief Project one atomic density, summing its original k weights once.
!> \param frames Original projector frames.
!> \param density Unweighted AO occupation matrices or their signed direction.
!> \param kweights Original normalized k weights.
!> \param first First projector of this atom.
!> \param last Last projector of this atom.
!> \param dm Overwritten spin-resolved real symmetric atomic projector density.
! **************************************************************************************************
   SUBROUTINE project_paw_density(frames, density, kweights, first, last, dm)
      TYPE(qs_upf_inverse_grid_type), INTENT(IN)         :: frames(:)
      COMPLEX(dp), INTENT(IN)                            :: density(:, :, :, :)
      REAL(dp), INTENT(IN)                               :: kweights(:)
      INTEGER, INTENT(IN)                                :: first, last
      REAL(dp), INTENT(OUT)                              :: dm(:, :, :)

      INTEGER                                            :: ik, nk, ns, spin
      COMPLEX(dp) :: projected(last - first + 1, last - first + 1)

      nk = SIZE(frames)
      ns = SIZE(density, 4)
      dm = 0.0_dp
      DO spin = 1, ns
         DO ik = 1, nk
            projected(:, :) = MATMUL(frames(ik)%ao%coefficients(first:last, :), &
                                     MATMUL(density(:, :, ik, spin), &
                                            CONJG(TRANSPOSE(frames(ik)%ao%coefficients(first:last, :)))))
            dm(:, :, spin) = dm(:, :, spin) + kweights(ik)*REAL(projected, dp)
         END DO
         dm(:, :, spin) = 0.5_dp*(dm(:, :, spin) + TRANSPOSE(dm(:, :, spin)))
      END DO
   END SUBROUTINE project_paw_density

! **************************************************************************************************
!> \brief Retain original orbital projector products when supplied by the installed-density owner.
!> \param frames Original projector maps, in physical k-point order.
!> \param kweights Original k weights.
!> \param first First projector of this atom.
!> \param last Last projector of this atom.
!> \param native_density Existing real projector matrix, used when original factors are unavailable.
!> \param components Retained real occupation or its response.
!> \param orbitals Optional original orbital factors of the installed density.
!> \param direction Optional matching coefficient and occupation variation.
! **************************************************************************************************
   SUBROUTINE retained_projector_density(frames, kweights, first, last, native_density, components, orbitals, direction)
      TYPE(qs_upf_inverse_grid_type), INTENT(IN)         :: frames(:)
      REAL(dp), INTENT(IN)                               :: kweights(:)
      INTEGER, INTENT(IN)                                :: first, last
      REAL(dp), INTENT(IN)                               :: native_density(:, :, :)
      TYPE(twofold_real_type), INTENT(OUT)               :: components(:, :, :)
      TYPE(qs_upf_inverse_orbitals_type), INTENT(IN), &
         OPTIONAL                                        :: orbitals, direction

      COMPLEX(dp), ALLOCATABLE                           :: fields(:, :, :)
      INTEGER                                            :: ik, spin
      LOGICAL                                            :: assertion_valid

      components%hi = native_density
      components%lo = 0.0_dp
      IF (.NOT. PRESENT(orbitals)) RETURN
      IF (PRESENT(direction)) THEN
         assertion_valid = direction%density_evaluation == orbitals%density_evaluation
         CPASSERT(assertion_valid)
      END IF
      assertion_valid = ALLOCATED(orbitals%coefficients) .AND. ALLOCATED(orbitals%occupations)
      CPASSERT(assertion_valid)
      assertion_valid = SIZE(orbitals%coefficients, 3) == SIZE(frames)
      CPASSERT(assertion_valid)
      assertion_valid = SIZE(orbitals%coefficients, 4) == SIZE(native_density, 3)
      CPASSERT(assertion_valid)
      ALLOCATE (fields(last - first + 1, SIZE(frames(1)%ao%coefficients, 2), SIZE(frames)))
      DO ik = 1, SIZE(frames)
         fields(:, :, ik) = frames(ik)%ao%coefficients(first:last, :)
      END DO
      DO spin = 1, SIZE(native_density, 3)
         IF (PRESENT(direction)) THEN
            CALL upf_inverse_orbital_projector_density(fields, orbitals%coefficients(:, :, :, spin), &
                                                       orbitals%occupations(:, :, spin), kweights, components(:, :, spin), &
                                                       direction%coefficients(:, :, :, spin), direction%occupations(:, :, spin))
         ELSE
            CALL upf_inverse_orbital_projector_density(fields, orbitals%coefficients(:, :, :, spin), &
                                                       orbitals%occupations(:, :, spin), kweights, components(:, :, spin))
         END IF
      END DO
   END SUBROUTINE retained_projector_density

! **************************************************************************************************
!> \brief Apply the PAW one-center semilocal XC kernel to a complete AO density direction.
!>        Sum the original k weights before the nonlinear atomic XC evaluation. The AO
!>        Hamiltonian responses carry no k weight. Frozen-core fields remain in the base
!>        functional and have zero density direction. Each atom has the same MPI owner as
!>        the original first-potential evaluation. Hartree is supplied separately.
!> \param kinds Original prepared PAW kinds, including fixed core fields.
!> \param frames Original projector frames in complete global k order.
!> \param density Unperturbed unweighted Hermitian AO occupation matrices, by k point and spin.
!> \param direction Matching signed Hermitian AO density direction.
!> \param kweights Original normalized k weights.
!> \param xc_section The original complete XC section.
!> \param hamiltonian New unweighted one-center XC Hamiltonian responses.
!> \param orbital_factors Optional original factors certified for the installed density.
!> \param orbital_direction Optional matching coefficient and occupation response.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_paw_xc_kernel(kinds, frames, density, direction, kweights, xc_section, hamiltonian, &
                                           orbital_factors, orbital_direction)
      TYPE(paw_kind_type), INTENT(IN)                    :: kinds(:)
      TYPE(qs_upf_inverse_grid_type), INTENT(IN)         :: frames(:)
      COMPLEX(dp), INTENT(IN)                            :: density(:, :, :, :), &
                                                            direction(:, :, :, :)
      REAL(dp), INTENT(IN)                               :: kweights(:)
      TYPE(section_vals_type), POINTER                   :: xc_section
      COMPLEX(dp), INTENT(OUT)                           :: hamiltonian(:, :, :, :)
      TYPE(qs_upf_inverse_orbitals_type), INTENT(IN), &
         OPTIONAL                                        :: orbital_factors, orbital_direction

      COMPLEX(dp), ALLOCATABLE                           :: atomic_operator(:, :)
      INTEGER                                            :: center, first, ik, kind, last, n, ns, &
                                                            spin
      REAL(dp), ALLOCATABLE                              :: delta(:, :, :), dm(:, :, :), &
                                                            potential(:, :, :)
      TYPE(twofold_real_type), ALLOCATABLE               :: density_components(:, :, :), &
                                                            direction_components(:, :, :), &
                                                            potential_components(:, :, :)

      CALL check_paw_density(frames, density, kweights)
      CALL check_paw_density(frames, direction, kweights)
      CPASSERT(ALL(SHAPE(direction) == SHAPE(density)))
      CPASSERT(ALL(SHAPE(hamiltonian) == SHAPE(density)))
      ns = SIZE(density, 4)
      hamiltonian = CMPLX(0.0_dp, 0.0_dp, dp)
      DO center = 1, SIZE(frames(1)%ao%center_atoms)
         kind = frames(1)%ao%center_kinds(center)
         CPASSERT(kind >= 1 .AND. kind <= SIZE(kinds))
         IF (.NOT. kinds(kind)%active) CYCLE
         first = frames(1)%ao%projector_first(center)
         last = frames(1)%ao%projector_last(center)
         n = last - first + 1
         CPASSERT(n == kinds(kind)%onecenter%nproj)
         ALLOCATE (potential(n, n, ns), SOURCE=0.0_dp)
         ALLOCATE (potential_components(n, n, ns), SOURCE=twofold_real_type())
         IF (MODULO(center - 1, frames(1)%ao%point_comm%num_pe) == frames(1)%ao%point_comm%mepos) THEN
            ALLOCATE (dm(n, n, ns), delta(n, n, ns))
            CALL project_paw_density(frames, density, kweights, first, last, dm)
            CALL project_paw_density(frames, direction, kweights, first, last, delta)
            ALLOCATE (density_components(n, n, ns), direction_components(n, n, ns))
            CALL retained_projector_density(frames, kweights, first, last, dm, density_components, orbital_factors)
            direction_components%hi = delta
            direction_components%lo = 0.0_dp
            IF (PRESENT(orbital_direction)) THEN
     CALL retained_projector_density(frames, kweights, first, last, delta, direction_components, orbital_factors, orbital_direction)
            END IF
            CALL paw_xc_response(kinds(kind)%xcgrid, kinds(kind)%onecenter, xc_section, dm, delta, potential, &
        density_components=density_components, direction_components=direction_components, potential_components=potential_components)
            DEALLOCATE (density_components, direction_components)
            DEALLOCATE (dm, delta)
         END IF
         CALL frames(1)%ao%point_comm%sum(potential_components%hi)
         CALL frames(1)%ao%point_comm%sum(potential_components%lo)
         ALLOCATE (atomic_operator(SIZE(density, 1), SIZE(density, 2)))
         DO spin = 1, ns
            DO ik = 1, SIZE(frames)
   CALL upf_inverse_projector_operator(frames(ik)%ao%coefficients(first:last, :), potential_components(:, :, spin), atomic_operator)
               hamiltonian(:, :, ik, spin) = hamiltonian(:, :, ik, spin) + atomic_operator
            END DO
         END DO
         DEALLOCATE (potential, potential_components, atomic_operator)
      END DO
   END SUBROUTINE qs_upf_inverse_paw_xc_kernel

! **************************************************************************************************
!> \brief Differentiate the same PAW functional through the native inverse transformation.
!> \param frames Native frames used for the energy and projector potential.
!> \param density Fixed Gaussian occupation matrices.
!> \param kweights Normalized global k weights.
!> \param potential PAW projector derivative at this density and geometry.
!> \param atom_response Cartesian displacements of the atomic positions.
!> \param cell_response Cell response with fixed fractional grid coordinates.
!> \param reduced_k_response Response of the reduced Bloch vectors, one per global k point.
!> \param derivative Energy response. Rigid atomic quadratures have no additional measure response.
!> \param coefficients_response Optional responses of the same frames and geometry variation, in global k order.
!> \param orbital_factors Optional original factors of the same installed density as potential.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_paw_response(frames, density, kweights, potential, atom_response, cell_response, &
                                          reduced_k_response, derivative, coefficients_response, orbital_factors)
      TYPE(qs_upf_inverse_grid_type), INTENT(IN)         :: frames(:)
      COMPLEX(dp), INTENT(IN)                            :: density(:, :, :, :)
      REAL(dp), INTENT(IN)                               :: kweights(:)
      TYPE(paw_potential_type), INTENT(IN)               :: potential(:)
      REAL(dp), INTENT(IN)                               :: atom_response(:, :), &
                                                            cell_response(3, 3), &
                                                            reduced_k_response(:, :)
      REAL(dp), INTENT(OUT)                              :: derivative
      COMPLEX(dp), INTENT(IN), OPTIONAL, TARGET          :: coefficients_response(:, :, :)
      TYPE(qs_upf_inverse_orbitals_type), INTENT(IN), &
         OPTIONAL                                        :: orbital_factors

      COMPLEX(dp), ALLOCATABLE                           :: projected(:, :)
      COMPLEX(dp), ALLOCATABLE, TARGET                   :: owned_dw(:, :)
      COMPLEX(dp), POINTER                               :: dw(:, :)
      INTEGER                                            :: center, first, ik, info, last, n, nao, &
                                                            nk, nproj, ns, p, q, spin
      LOGICAL                                            :: assertion_valid
      REAL(dp)                                           :: inverse_cell(3, 3), strain(3, 3), &
                                                            volume_response
      REAL(dp), ALLOCATABLE                              :: points_response(:, :), &
                                                            weights_response(:)
      TYPE(twofold_real_type)                            :: retained_response
      TYPE(twofold_real_type), ALLOCATABLE               :: occupation_response(:, :)

      nk = SIZE(frames)
      CPASSERT(nk > 0 .AND. SIZE(kweights) == nk)
      nao = SIZE(frames(1)%ao%coefficients, 2)
      nproj = SIZE(frames(1)%ao%coefficients, 1)
      ns = SIZE(density, 4)
      CPASSERT(ALL(SHAPE(density) == [nao, nao, nk, ns]))
      CPASSERT(SIZE(potential) == SIZE(frames(1)%ao%center_atoms))
      CPASSERT(ALL(SHAPE(reduced_k_response) == [3, nk]))
      CPASSERT(ALL(ieee_is_finite(kweights)) .AND. ALL(kweights >= 0.0_dp))
      CPASSERT(ABS(SUM(kweights) - 1.0_dp) < 1.e-12_dp)
      IF (PRESENT(coefficients_response)) THEN
         CPASSERT(ALL(SHAPE(coefficients_response) == [nproj, nao, nk]))
      ELSE
         ALLOCATE (owned_dw(nproj, nao))
      END IF
      IF (PRESENT(orbital_factors)) THEN
         assertion_valid = ALLOCATED(orbital_factors%coefficients) .AND. ALLOCATED(orbital_factors%occupations)
         CPASSERT(assertion_valid)
         assertion_valid = SIZE(orbital_factors%coefficients, 1) == nao .AND. &
                           SIZE(orbital_factors%coefficients, 3) == nk .AND. SIZE(orbital_factors%coefficients, 4) == ns
         CPASSERT(assertion_valid)
         assertion_valid = ALL(SHAPE(orbital_factors%occupations) == [SIZE(orbital_factors%coefficients, 2), nk, ns])
         CPASSERT(assertion_valid)
      END IF
      derivative = 0.0_dp
      retained_response = twofold_real_type()
      DO ik = 1, nk
         CPASSERT(ALL(SHAPE(atom_response) == SHAPE(frames(ik)%ao%positions)))
         IF (PRESENT(coefficients_response)) THEN
            dw => coefficients_response(:, :, ik)
         ELSE
            dw => owned_dw
            inverse_cell = frames(ik)%ao%cell
            CALL invmat(inverse_cell, info)
            CPASSERT(info == 0)
            strain = MATMUL(cell_response, inverse_cell)
            volume_response = strain(1, 1) + strain(2, 2) + strain(3, 3)
            ALLOCATE (points_response(3, SIZE(frames(ik)%ao%weights)), &
                      weights_response(SIZE(frames(ik)%ao%weights)))
            points_response(:, :) = MATMUL(strain, frames(ik)%ao%points)
            weights_response(:) = volume_response*frames(ik)%ao%weights
            CALL qs_upf_inverse_response(frames(ik)%ao, points_response, atom_response, cell_response, &
                                         reduced_k_response(:, ik), weights_response, dw)
            DEALLOCATE (points_response, weights_response)
         END IF
         DO center = 1, SIZE(potential)
            IF (.NOT. ALLOCATED(potential(center)%matrix)) CYCLE
            first = frames(ik)%ao%projector_first(center)
            last = frames(ik)%ao%projector_last(center)
            CPASSERT(SIZE(potential(center)%matrix, 1) == last - first + 1)
            CPASSERT(SIZE(potential(center)%matrix, 2) == last - first + 1)
            CPASSERT(SIZE(potential(center)%matrix, 3) == ns)
            IF (PRESENT(orbital_factors) .AND. ALLOCATED(potential(center)%components)) THEN
               n = last - first + 1
               assertion_valid = ALL(SHAPE(potential(center)%components) == [n, n, ns])
               CPASSERT(assertion_valid)
               ALLOCATE (occupation_response(n, n))
               DO spin = 1, ns
                  CALL upf_inverse_orbital_projector_density( &
                     RESHAPE(frames(ik)%ao%coefficients(first:last, :), [n, nao, 1]), &
                     orbital_factors%coefficients(:, :, ik:ik, spin), orbital_factors%occupations(:, ik:ik, spin), &
                     [kweights(ik)], occupation_response, &
                     fields_direction=RESHAPE(dw(first:last, :), [n, nao, 1]))
                  DO q = 1, n
                     DO p = 1, n
                        retained_response = twofold_add(retained_response, &
                                              twofold_multiply(potential(center)%components(p, q, spin), occupation_response(q, p)))
                     END DO
                  END DO
               END DO
               DEALLOCATE (occupation_response)
               CYCLE
            END IF
            ALLOCATE (projected(last - first + 1, last - first + 1))
            DO spin = 1, ns
               projected(:, :) = MATMUL(dw(first:last, :), MATMUL(density(:, :, ik, spin), &
                                                                  CONJG(TRANSPOSE(frames(ik)%ao%coefficients(first:last, :)))))
               derivative = derivative + 2.0_dp*kweights(ik)* &
                            SUM(potential(center)%matrix(:, :, spin)*REAL(TRANSPOSE(projected), dp))
            END DO
            DEALLOCATE (projected)
         END DO
      END DO
      derivative = derivative + twofold_value(retained_response)
   END SUBROUTINE qs_upf_inverse_paw_response

END MODULE qs_upf_inverse_paw
