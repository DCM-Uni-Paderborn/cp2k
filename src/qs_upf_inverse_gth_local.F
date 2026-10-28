!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Native analytic Gaussian GTH PPL matrices and their isolated geometry adjoints.
!>        This is the full Gaussian base contribution. The inverse reconstruction and original
!>        UPF augmentation corrections remain separate and must be combined by the caller.
! **************************************************************************************************
MODULE qs_upf_inverse_gth_local
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE core_ppl,                        ONLY: build_core_ppl
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type,&
                                              dbcsr_set
   USE gth_local_fields,                ONLY: gth_local_fields_evaluate,&
                                              gth_local_fields_init,&
                                              gth_local_fields_type
   USE kinds,                           ONLY: dp
   USE kpoint_types,                    ONLY: kpoint_type
   USE mathlib,                         ONLY: invmat
   USE particle_types,                  ONLY: particle_type
   USE paw_types,                       ONLY: paw_kind_type
   USE pw_types,                        ONLY: pw_r3d_rs_type
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_upf_inverse_background,       ONLY: qs_upf_inverse_background_operator,&
                                              qs_upf_inverse_background_set_density
   USE qs_upf_inverse_grid,             ONLY: qs_upf_inverse_grid_gaussian_operator,&
                                              qs_upf_inverse_grid_integrate_dense,&
                                              qs_upf_inverse_grid_response
   USE qs_upf_inverse_paw,              ONLY: qs_upf_inverse_paw_local_operator
   USE qs_upf_inverse_types,            ONLY: qs_upf_inverse_augmentation_type,&
                                              qs_upf_inverse_background_type,&
                                              qs_upf_inverse_grid_type,&
                                              qs_upf_inverse_type
   USE virial_types,                    ONLY: virial_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: qs_upf_inverse_gth_local_matrix, qs_upf_inverse_gth_local_geometry
   PUBLIC :: qs_upf_inverse_gth_local_active, qs_upf_inverse_gth_local_field
   PUBLIC :: qs_upf_inverse_gth_local_residual
   PUBLIC :: qs_upf_inverse_gth_local_paw

CONTAINS

! **************************************************************************************************
!> \brief Add the physical PAW density and frozen-core coupling to original local GTH fields.
!>        The valence operator is W^dagger A W, where A is the original AE-minus-compensated-PS
!>        scalar-field transpose. Atomic radial nodes, directions and weights move rigidly with
!>        their PAW center. Only relative centers and physical image vectors enter the field
!>        response. There is no grid-volume, Bloch-phase or atomic-rule Jacobian contribution.
!> \param frame Original centers, fixed physical images and reconstruction coefficients W.
!> \param kinds Original native scalar GTH datasets, without fitting or additional truncation.
!> \param particles Original atom-to-kind mapping in the frame's atom order.
!> \param paw Prepared original PAW density products and frozen AE-core densities by kind.
!> \param OPERATOR Optional existing AO operator, incremented without spin or k-point weights.
!> \param core_energy Optional existing frozen-core constant, incremented once outside spin/k loops.
!> \param response Optional existing AO operator response, including both responses of W.
!> \param core_response Optional existing frozen-core response, independent of W and valence density.
!> \param coefficients_response Required with response, complete derivative of W.
!> \param atom_response Required with either response, Cartesian motion of every physical atom.
!> \param cell_response Required with either response, motion of the physical image lattice.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_gth_local_paw(frame, kinds, particles, paw, OPERATOR, core_energy, response, &
                                           core_response, coefficients_response, atom_response, cell_response)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      TYPE(qs_kind_type), INTENT(IN)                     :: kinds(:)
      TYPE(particle_type), INTENT(IN)                    :: particles(:)
      TYPE(paw_kind_type), INTENT(IN)                    :: paw(:)
      COMPLEX(dp), INTENT(INOUT), OPTIONAL               :: operator(:, :)
      REAL(dp), INTENT(INOUT), OPTIONAL                  :: core_energy
      COMPLEX(dp), INTENT(INOUT), OPTIONAL               :: response(:, :)
      REAL(dp), INTENT(INOUT), OPTIONAL                  :: core_response
      COMPLEX(dp), INTENT(IN), OPTIONAL                  :: coefficients_response(:, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: atom_response(:, :), cell_response(3, 3)

      COMPLEX(dp), ALLOCATABLE                           :: weighted(:, :)
      INTEGER                                            :: active, atom, center_atom, gth_kind, &
                                                            image, ir, kind, n, na, nr, p, q
      LOGICAL                                            :: geometry, values_requested
      REAL(dp)                                           :: center(3), core, dcore, motion(3)
      REAL(dp), ALLOCATABLE :: displacement(:, :), field(:, :), field_response(:, :), &
         gradient(:, :), matrix(:, :), matrix_response(:, :), values(:)
      TYPE(gth_local_fields_type), ALLOCATABLE           :: fields(:)
      LOGICAL                                            :: selected(SIZE(kinds))

      CPASSERT(SIZE(paw) == SIZE(kinds))
      CPASSERT(ALL(SHAPE(frame%positions) == [3, SIZE(particles)]))
      CPASSERT(SIZE(frame%images, 1) == 3)
      geometry = PRESENT(response) .OR. PRESENT(core_response)
      values_requested = PRESENT(OPERATOR) .OR. PRESENT(core_energy)
      CPASSERT(values_requested .OR. geometry)
      n = SIZE(frame%coefficients, 2)
      IF (PRESENT(OPERATOR)) THEN
         CPASSERT(ALL(SHAPE(OPERATOR) == [n, n]))
      END IF
      IF (PRESENT(response)) THEN
         CPASSERT(ALL(SHAPE(response) == [n, n]))
         CPASSERT(PRESENT(coefficients_response))
         CPASSERT(ALL(SHAPE(coefficients_response) == SHAPE(frame%coefficients)))
      END IF
      IF (geometry) THEN
         CPASSERT(PRESENT(atom_response))
         CPASSERT(ALL(SHAPE(atom_response) == SHAPE(frame%positions)))
         CPASSERT(PRESENT(cell_response))
      END IF
      selected = .FALSE.
      ALLOCATE (fields(SIZE(kinds)))
      DO kind = 1, SIZE(kinds)
         IF (kinds(kind)%ghost) CYCLE
         IF (.NOT. ASSOCIATED(kinds(kind)%gth_potential)) CYCLE
         CALL gth_local_fields_init(kinds(kind)%gth_potential, fields(kind))
         selected(kind) = ANY(fields(kind)%nterm > 0)
      END DO
      IF (.NOT. ANY(selected)) RETURN
      DO active = 1, SIZE(frame%center_atoms)
         kind = frame%center_kinds(active)
         IF (.NOT. paw(kind)%active) CYCLE
         center_atom = frame%center_atoms(active)
         p = frame%projector_first(active)
         q = frame%projector_last(active)
         CPASSERT(q - p + 1 == paw(kind)%onecenter%nproj)
         ASSOCIATE (grid => paw(kind)%xcgrid)
            na = grid%na
            nr = grid%nr
            ALLOCATE (displacement(3, na), values(na), field(na, nr), matrix(q - p + 1, q - p + 1))
            field = 0.0_dp
            IF (geometry) ALLOCATE (gradient(na, 3), field_response(na, nr), SOURCE=0.0_dp)
            DO atom = 1, SIZE(particles)
               CPASSERT(ASSOCIATED(particles(atom)%atomic_kind))
               gth_kind = particles(atom)%atomic_kind%kind_number
               CPASSERT(gth_kind >= 1 .AND. gth_kind <= SIZE(kinds))
               IF (.NOT. selected(gth_kind)) CYCLE
               DO image = 1, SIZE(frame%images, 2)
                  center = frame%positions(:, atom) + MATMUL(frame%cell, REAL(frame%images(:, image), dp))
                  IF (geometry) THEN
                     motion = atom_response(:, center_atom) - atom_response(:, atom) - &
                              MATMUL(cell_response, REAL(frame%images(:, image), dp))
                  END IF
                  DO ir = 1, nr
                     displacement = grid%r(ir)*grid%directions + &
                                    SPREAD(frame%positions(:, center_atom) - center, 2, na)
                     IF (geometry) THEN
                        CALL gth_local_fields_evaluate(fields(gth_kind), displacement, values, gradient)
                        field_response(:, ir) = field_response(:, ir) + MATMUL(gradient, motion)
                     ELSE
                        CALL gth_local_fields_evaluate(fields(gth_kind), displacement, values)
                     END IF
                     field(:, ir) = field(:, ir) + values
                  END DO
               END DO
            END DO
         END ASSOCIATE
         CALL qs_upf_inverse_paw_local_operator(paw(kind), field, matrix, core)
         IF (PRESENT(core_energy)) core_energy = core_energy + core
         IF (PRESENT(OPERATOR) .OR. PRESENT(response)) THEN
            weighted = MATMUL(CMPLX(matrix, 0.0_dp, dp), frame%coefficients(p:q, :))
            IF (PRESENT(OPERATOR)) THEN
               OPERATOR = OPERATOR + MATMUL(CONJG(TRANSPOSE(frame%coefficients(p:q, :))), weighted)
            END IF
         END IF
         IF (geometry) THEN
            ALLOCATE (matrix_response, MOLD=matrix)
            CALL qs_upf_inverse_paw_local_operator(paw(kind), field_response, matrix_response, dcore)
            IF (PRESENT(core_response)) core_response = core_response + dcore
            IF (PRESENT(response)) THEN
               response = response + MATMUL(CONJG(TRANSPOSE(coefficients_response(p:q, :))), weighted) + &
                          MATMUL(CONJG(TRANSPOSE(frame%coefficients(p:q, :))), &
                                 MATMUL(CMPLX(matrix_response, 0.0_dp, dp), frame%coefficients(p:q, :)) + &
                                 MATMUL(CMPLX(matrix, 0.0_dp, dp), coefficients_response(p:q, :)))
            END IF
            DEALLOCATE (gradient, field_response, matrix_response)
         END IF
         DEALLOCATE (displacement, values, field, matrix)
      END DO
   END SUBROUTINE qs_upf_inverse_gth_local_paw

! **************************************************************************************************
!> \brief Whether any physical GTH center has native local Gaussian terms.
!> \param kinds ...
!> \return ...
! **************************************************************************************************
   FUNCTION qs_upf_inverse_gth_local_active(kinds) RESULT(active)
      TYPE(qs_kind_type), INTENT(IN)                     :: kinds(:)
      LOGICAL                                            :: active

      INTEGER                                            :: kind

      active = .FALSE.
      DO kind = 1, SIZE(kinds)
         IF (kinds(kind)%ghost) CYCLE
         IF (.NOT. ASSOCIATED(kinds(kind)%gth_potential)) CYCLE
         ASSOCIATE (pot => kinds(kind)%gth_potential)
            active = active .OR. pot%nexp_ppl > 0
            IF (pot%lpotextended) active = active .OR. pot%nexp_lpot > 0
         END ASSOCIATE
      END DO
   END FUNCTION qs_upf_inverse_gth_local_active

! **************************************************************************************************
!> \brief Original local GTH Gaussians on the defining inverse grid and physical image set.
!>        Local potentials carry no orbital Bloch phase. Use only the supplied physical images,
!>        including for isolated and partially periodic cells. The ionic charge stays separate.
!> \param ks_env Original potential data, without modification.
!> \param frame Native point ordering, physical centers, cell and fixed image membership.
!> \param potential Allocated scalar grid field, overwritten.
!> \param response Optional derivative at moving fractional-grid points.
!> \param atom_response Cartesian motion of every original atomic center.
!> \param cell_response Motion of the physical cell vectors and image centers.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_gth_local_field(ks_env, frame, potential, response, atom_response, cell_response)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(qs_upf_inverse_grid_type), INTENT(IN)         :: frame
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: potential
      TYPE(pw_r3d_rs_type), INTENT(INOUT), OPTIONAL      :: response
      REAL(dp), INTENT(IN), OPTIONAL                     :: atom_response(:, :), cell_response(3, 3)

      INTEGER                                            :: atom, image, info, kind, n
      REAL(dp)                                           :: cell_inverse(3, 3), center(3), &
                                                            motion(3), strain(3, 3)
      REAL(dp), ALLOCATABLE                              :: displacement(:, :), gradient(:, :), &
                                                            points_response(:, :), total(:), &
                                                            total_response(:), values(:)
      TYPE(gth_local_fields_type)                        :: fields
      TYPE(particle_type), POINTER                       :: atoms(:)
      TYPE(qs_kind_type), POINTER                        :: kinds(:)

      CALL get_ks_env(ks_env, qs_kind_set=kinds, particle_set=atoms)
      n = SIZE(frame%ao%weights)
      CPASSERT(SIZE(potential%array) == n)
      ALLOCATE (displacement(3, n), values(n), total(n), SOURCE=0.0_dp)
      IF (PRESENT(response)) THEN
         CPASSERT(PRESENT(atom_response))
         CPASSERT(PRESENT(cell_response))
         CPASSERT(SIZE(response%array) == n)
         CPASSERT(ALL(SHAPE(atom_response) == SHAPE(frame%ao%positions)))
         cell_inverse = frame%ao%cell
         CALL invmat(cell_inverse, info)
         CPASSERT(info == 0)
         strain = MATMUL(cell_response, cell_inverse)
         points_response = MATMUL(strain, frame%ao%points)
         ALLOCATE (gradient(n, 3), total_response(n), SOURCE=0.0_dp)
      END IF
      DO atom = 1, SIZE(atoms)
         kind = atoms(atom)%atomic_kind%kind_number
         IF (kinds(kind)%ghost) CYCLE
         IF (.NOT. ASSOCIATED(kinds(kind)%gth_potential)) CYCLE
         CALL gth_local_fields_init(kinds(kind)%gth_potential, fields)
         DO image = 1, SIZE(frame%ao%images, 2)
            center = frame%ao%positions(:, atom) + MATMUL(frame%ao%cell, REAL(frame%ao%images(:, image), dp))
            displacement = frame%ao%points - SPREAD(center, 2, n)
            IF (PRESENT(response)) THEN
               motion = atom_response(:, atom) + MATMUL(cell_response, REAL(frame%ao%images(:, image), dp))
               CALL gth_local_fields_evaluate(fields, displacement, values, gradient)
               total_response = total_response + SUM(gradient*TRANSPOSE(points_response - SPREAD(motion, 2, n)), DIM=2)
            ELSE
               CALL gth_local_fields_evaluate(fields, displacement, values)
            END IF
            total = total + values
         END DO
      END DO
      potential%array = RESHAPE(total, SHAPE(potential%array))
      IF (PRESENT(response)) response%array = RESHAPE(total_response, SHAPE(response%array))
   END SUBROUTINE qs_upf_inverse_gth_local_field

! **************************************************************************************************
!> \brief Retain the full inverse and original-Q local terms after removing their Gaussian grid base.
!>        The caller separately adds the native analytic Gaussian PPL matrix and its adjoint.
!>        DIRECT with native Gaussian tasks contributes only Q here. All other representations
!>        subtract precisely the explicit Gaussian base of this same point rule and image set.
!> \param frame Original inverse fields, without changing W.
!> \param augmentation Original UPF Q kernels, with their native Fourier charge integration.
!> \param potential Local GTH scalar field on the defining grid.
!> \param orbital_fields Include explicit inverse orbitals. False is the native DIRECT base.
!> \param OPERATOR Optional new residual matrix, without occupation or k weights.
!> \param response Optional joint residual derivative, without a duplicate metric derivative.
!> \param potential_response Local field derivative at moving grid points.
!> \param atom_response Motion of every atomic center.
!> \param cell_response Motion of the physical cell vectors.
!> \param reduced_k_response Motion of this frame's reduced Bloch vector.
!> \param coefficients_response Reuse the defining inverse coefficient response for this geometry.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_gth_local_residual(frame, augmentation, potential, orbital_fields, OPERATOR, response, &
                                                potential_response, atom_response, cell_response, reduced_k_response, &
                                                coefficients_response)
      TYPE(qs_upf_inverse_grid_type), INTENT(IN)         :: frame
      TYPE(qs_upf_inverse_augmentation_type), &
         INTENT(INOUT)                                   :: augmentation
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: potential
      LOGICAL, INTENT(IN)                                :: orbital_fields
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT), OPTIONAL    :: operator(:, :), response(:, :)
      TYPE(pw_r3d_rs_type), INTENT(IN), OPTIONAL         :: potential_response
      REAL(dp), INTENT(IN), OPTIONAL                     :: atom_response(:, :), &
                                                            cell_response(3, 3), &
                                                            reduced_k_response(3)
      COMPLEX(dp), INTENT(IN), OPTIONAL                  :: coefficients_response(:, :)

      COMPLEX(dp), ALLOCATABLE                           :: gaussian(:, :), metric(:, :)
      INTEGER                                            :: axis
      TYPE(pw_r3d_rs_type)                               :: zero, zeros(3)

      CPASSERT(PRESENT(OPERATOR) .OR. PRESENT(response))
      CALL zero%create(potential%pw_grid)
      zero%array = 0.0_dp
      DO axis = 1, 3
         CALL zeros(axis)%create(potential%pw_grid)
         zeros(axis)%array = 0.0_dp
      END DO
      IF (PRESENT(OPERATOR)) THEN
         CALL qs_upf_inverse_grid_integrate_dense(frame, potential, zeros, zero, OPERATOR, augmentation=augmentation, &
                                                  fourier_charge=.TRUE., orbital_fields=orbital_fields)
         IF (orbital_fields) THEN
            CALL qs_upf_inverse_grid_gaussian_operator(frame, potential, gaussian)
            OPERATOR = OPERATOR - gaussian
         END IF
      END IF
      IF (PRESENT(response)) THEN
         CPASSERT(PRESENT(potential_response))
         CPASSERT(PRESENT(atom_response))
         CPASSERT(PRESENT(cell_response))
         CPASSERT(PRESENT(reduced_k_response))
         CPASSERT(PRESENT(coefficients_response))
         CALL qs_upf_inverse_grid_response(frame, augmentation, frame%ao%radial_q, potential, zeros, zero, &
                                           atom_response, cell_response, reduced_k_response, potential_response, zeros, zero, &
                                           response, metric, orbital_fields=orbital_fields, &
                                           known_coefficients_response=coefficients_response)
         IF (orbital_fields) THEN
            CALL qs_upf_inverse_grid_gaussian_operator(frame, potential, response=gaussian, &
                                                       potential_response=potential_response, atom_response=atom_response, &
                                                       cell_response=cell_response, reduced_k_response=reduced_k_response)
            response = response - gaussian
         END IF
      END IF
      CALL zero%release()
      DO axis = 1, 3
         CALL zeros(axis)%release()
      END DO
   END SUBROUTINE qs_upf_inverse_gth_local_residual

! **************************************************************************************************
!> \brief Read the original native GTH PPL integral matrices in the physical Bloch AO ordering.
!> \param ks_env Native kinds, full orbital bases and unchanged neighbor lists.
!> \param map Owned image scratch matrices and the native-to-physical Bloch phase map.
!> \param operators New full-Gaussian PPL matrices, without k-point weights.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_gth_local_matrix(ks_env, map, operators)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(qs_upf_inverse_background_type), &
         INTENT(INOUT)                                   :: map
      COMPLEX(dp), INTENT(OUT)                           :: operators(:, :, :)

      CALL native_local(ks_env, map, operators=operators)
   END SUBROUTINE qs_upf_inverse_gth_local_matrix

! **************************************************************************************************
!> \brief Native PPL energy gradients and virial at explicitly supplied accepted AO occupations.
!>        The native force, virial, density and operator owners are unchanged. Only the owned map's
!>        image scratch arrays are overwritten. Full Gaussian bases are used for every density mode.
!> \param ks_env Original native integral environment.
!> \param map Owned native image and phase map, also carrying the original k-point weights.
!> \param matrix Hermitian occupations indexed by global k, without mesh weights.
!> \param gradients Replicated Cartesian energy gradients, without a force sign change.
!> \param virial Replicated native PPL virial.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_gth_local_geometry(ks_env, map, matrix, gradients, virial)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(qs_upf_inverse_background_type), &
         INTENT(INOUT)                                   :: map
      COMPLEX(dp), INTENT(IN)                            :: matrix(:, :, :)
      REAL(dp), INTENT(OUT)                              :: gradients(:, :), virial(3, 3)

      CALL native_local(ks_env, map, matrix=matrix, gradients=gradients, local_virial=virial)
   END SUBROUTINE qs_upf_inverse_gth_local_geometry

! **************************************************************************************************
!> \brief Select native GTH potential centers while retaining all orbital centers and contractions.
!> \param ks_env ...
!> \param map ...
!> \param operators ...
!> \param matrix ...
!> \param gradients ...
!> \param local_virial ...
! **************************************************************************************************
   SUBROUTINE native_local(ks_env, map, operators, matrix, gradients, local_virial)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(qs_upf_inverse_background_type), &
         INTENT(INOUT)                                   :: map
      COMPLEX(dp), INTENT(OUT), OPTIONAL                 :: operators(:, :, :)
      COMPLEX(dp), INTENT(IN), OPTIONAL                  :: matrix(:, :, :)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: gradients(:, :), local_virial(3, 3)

      INTEGER                                            :: image, kind, nder, nimages
      INTEGER, POINTER                                   :: cell_to_index(:, :, :)
      LOGICAL                                            :: calculate_forces, use_virial
      LOGICAL, ALLOCATABLE                               :: selected(:)
      TYPE(atomic_kind_type), POINTER                    :: atomic_kinds(:)
      TYPE(dbcsr_p_type), POINTER                        :: density(:, :), potential(:, :)
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(neighbor_list_set_p_type), POINTER            :: centers(:), overlap(:)
      TYPE(particle_type), POINTER                       :: atoms(:)
      TYPE(qs_force_type), POINTER                       :: force(:)
      TYPE(qs_kind_type), POINTER                        :: kinds(:)
      TYPE(virial_type), POINTER                         :: virial

      NULLIFY (force, virial)
      CPASSERT(ASSOCIATED(map%potential))
      CALL get_ks_env(ks_env, atomic_kind_set=atomic_kinds, qs_kind_set=kinds, particle_set=atoms, &
                      sab_orb=overlap, sac_ppl=centers, kpoints=kpoints)
      calculate_forces = PRESENT(matrix)
      IF (PRESENT(operators)) operators = CMPLX(0.0_dp, 0.0_dp, dp)
      IF (calculate_forces) THEN
         CPASSERT(PRESENT(gradients))
         CPASSERT(PRESENT(local_virial))
         CPASSERT(SIZE(gradients, 1) == 3)
         CPASSERT(SIZE(gradients, 2) == SIZE(atoms))
         gradients = 0.0_dp
         local_virial = 0.0_dp
      END IF
      ALLOCATE (selected(SIZE(kinds)), SOURCE=.FALSE.)
      DO kind = 1, SIZE(kinds)
         IF (kinds(kind)%ghost) CYCLE
         IF (.NOT. ASSOCIATED(kinds(kind)%gth_potential)) CYCLE
         ASSOCIATE (pot => kinds(kind)%gth_potential)
            selected(kind) = pot%nexp_ppl > 0
            IF (pot%lpotextended) selected(kind) = selected(kind) .OR. pot%nexp_lpot > 0
         END ASSOCIATE
      END DO
      IF (.NOT. ANY(selected)) RETURN
      CPASSERT(ASSOCIATED(overlap))
      CPASSERT(ASSOCIATED(centers))
      IF (calculate_forces) CALL qs_upf_inverse_background_set_density(map, matrix)
      nimages = SIZE(map%potential)
      ALLOCATE (potential(1, nimages), density(1, nimages))
      DO image = 1, nimages
         CALL dbcsr_set(map%potential(image)%matrix, 0.0_dp)
         potential(1, image)%matrix => map%potential(image)%matrix
         density(1, image)%matrix => map%density(image)%matrix
      END DO
      NULLIFY (cell_to_index)
      IF (ASSOCIATED(kpoints)) cell_to_index => kpoints%cell_to_index
      nder = 0
      IF (calculate_forces) nder = 1
      use_virial = calculate_forces
      CALL build_core_ppl(potential, density, force, virial, calculate_forces, use_virial, nder, &
                          kinds, atomic_kinds, atoms, overlap, centers, nimages, cell_to_index, "ORB", &
                          kind_mask=selected, force_response=gradients, virial_response=local_virial)
      IF (PRESENT(operators)) CALL qs_upf_inverse_background_operator(map, map%potential, operators)
      IF (calculate_forces) THEN
         CALL map%comm%sum(gradients)
         CALL map%comm%sum(local_virial)
      END IF
      DEALLOCATE (potential, density)
   END SUBROUTINE native_local
END MODULE qs_upf_inverse_gth_local
