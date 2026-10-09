!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Metric-aware little-group representations of Gaussian Bloch band subspaces.
!> \author Thomas D. Kuehne, tkuehne@cp2k.org
! **************************************************************************************************
MODULE qs_little_group
   USE atomic_kind_types,               ONLY: get_atomic_kind
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cryssym,                         ONLY: crys_sym_gen,&
                                              csym_type,&
                                              release_csym_type
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE mathconstants,                   ONLY: twopi
   USE message_passing,                 ONLY: mp_para_env_type
   USE orbital_transformation_matrices, ONLY: calculate_rotmat,&
                                              orbrotmat_type,&
                                              release_rotmat
   USE particle_types,                  ONLY: particle_type
   USE qs_basis_rotation_methods,       ONLY: set_rotmat_basis
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type
   USE topology_atomic_signatures,      ONLY: atomic_signature_construct,&
                                              atomic_signature_from_bands,&
                                              atomic_signature_type
   USE topology_band_representations,   ONLY: atomic_site_storage
   USE topology_compatibility,          ONLY: little_group_compatibility,&
                                              little_group_data_type,&
                                              little_group_edge_type,&
                                              little_group_storage
   USE topology_compatibility_lattice,  ONLY: compatibility_lattice_construct,&
                                              compatibility_lattice_type
   USE topology_corepresentations,      ONLY: corepresentation_multiplicities,&
                                              corepresentation_table_type,&
                                              projective_corepresentations
   USE topology_integer_lattice,        ONLY: integer_lattice_factor,&
                                              integer_lattice_solve,&
                                              integer_lattice_type,&
                                              integer_quotient_class,&
                                              nonnegative_integer_solve
   USE topology_little_group,           ONLY: little_group_factor_system,&
                                              projective_characters
   USE topology_reciprocal,             ONLY: reciprocal_refine,&
                                              reciprocal_sampling,&
                                              reciprocal_sampling_type
   USE topology_star_transport,         ONLY: star_group_storage
   USE topology_symmetry,               ONLY: character_multiplicities,&
                                              check_sewing_product,&
                                              subspace_sewing,&
                                              transform_spinor_frame
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_little_group'
   PUBLIC :: write_little_group, write_little_group_edges, write_atomic_signatures
   PUBLIC :: prepare_little_group_sampling, write_little_group_sampling, write_little_group_graph
CONTAINS

! **************************************************************************************************
!> \brief Generate and refine native property points from the actual input-cell symmetry.
!> \param qs_env Converged Gaussian environment
!> \param time_reversal Include grey antiunitary operations
!> \param eps_geometry Geometry tolerance in Bohr
!> \param tolerance Fractional geometry tolerance
!> \param path_points Number of points per original segment, including both endpoints
!> \param max_points Conservative point/work limit before deduplication
!> \param sampling Explicit points and directed edges with reciprocal shifts
!> \param iw Summary output unit
! **************************************************************************************************
   SUBROUTINE prepare_little_group_sampling(qs_env, time_reversal, eps_geometry, tolerance, &
                                            path_points, max_points, sampling, iw)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: time_reversal
      REAL(dp), INTENT(IN)                               :: eps_geometry, tolerance
      INTEGER, INTENT(IN)                                :: path_points, max_points
      TYPE(reciprocal_sampling_type), INTENT(OUT)        :: sampling
      INTEGER, INTENT(IN)                                :: iw

      INTEGER                                            :: status
      INTEGER, ALLOCATABLE                               :: rotations(:, :, :)
      LOGICAL, ALLOCATABLE                               :: anti(:)
      REAL(dp), ALLOCATABLE                              :: translations(:, :)
      TYPE(cell_type), POINTER                           :: cell

      IF (path_points < 2 .OR. max_points < 1) CPABORT('Invalid automatic little-group sampling limits.')
      CALL get_qs_env(qs_env, cell=cell)
      IF (ANY(cell%perd /= 1)) THEN
         CPABORT('Automatic little-group sampling requires three periodic directions.')
      END IF
      CALL crystal_operations(qs_env, time_reversal, eps_geometry, rotations, translations, anti)
      CALL reciprocal_sampling(cell%hmat, rotations, translations, anti, sampling, tolerance, status)
      IF (status /= 0) CPABORT('Automatic reciprocal symmetry geometry failed.')
      CALL reciprocal_refine(sampling, path_points - 1, max_points, tolerance, status)
      IF (status /= 0) THEN
         CPABORT('Automatic reciprocal refinement failed or exceeds LITTLE_GROUP_MAX_POINTS.')
      END IF
      IF (iw > 0) THEN
         WRITE (iw, '(T2,A,I0)') 'TQC| Automatic fixed families: ', SIZE(sampling%strata)
         WRITE (iw, '(T2,A,I0)') 'TQC| Automatic property points: ', SIZE(sampling%kpoints, 2)
         WRITE (iw, '(T2,A,I0)') 'TQC| Automatic directed segments: ', SIZE(sampling%source)
         WRITE (iw, '(T2,A)') 'TQC| Refine paths and check band isolation; sampled symmetry is not a full topology classification.'
      END IF
   END SUBROUTINE prepare_little_group_sampling

! **************************************************************************************************
!> \brief Retain the automatic geometry and path definition alongside the validated band data.
!> \param sampling Generated/refined property sampling
!> \param unit Detailed output unit, nonpositive on nonwriting ranks
! **************************************************************************************************
   SUBROUTINE write_little_group_sampling(sampling, unit)
      TYPE(reciprocal_sampling_type), INTENT(IN)         :: sampling
      INTEGER, INTENT(IN)                                :: unit

      INTEGER                                            :: i, j

      IF (unit <= 0) RETURN
      WRITE (unit, '(A,3I8)') 'AUTO_SIZES ', SIZE(sampling%strata), SIZE(sampling%kpoints, 2), SIZE(sampling%source)
      DO i = 1, SIZE(sampling%strata)
         WRITE (unit, '(A,3I8,9I8,3ES24.16)') 'AUTO_FAMILY ', i, sampling%strata(i)%dimension, &
            sampling%strata(i)%star_size, TRANSPOSE(sampling%strata(i)%normal), sampling%strata(i)%offset
         DO j = 1, sampling%strata(i)%dimension
            WRITE (unit, '(A,5I8)') 'AUTO_CYCLE ', i, j, sampling%strata(i)%cycles(:, j)
         END DO
      END DO
      DO i = 1, SIZE(sampling%kpoints, 2)
         WRITE (unit, '(A,I8,3ES24.16)') 'AUTO_POINT ', i, sampling%kpoints(:, i)
      END DO
      DO i = 1, SIZE(sampling%source)
         WRITE (unit, '(A,6I8)') 'AUTO_EDGE ', i, sampling%source(i), sampling%destination(i), sampling%shift(:, i)
      END DO
      WRITE (unit, '(A)') 'END_AUTO_SAMPLING'
   END SUBROUTINE write_little_group_sampling

! **************************************************************************************************
!> \brief Build and validate little-group band representations, without assigning EBR labels.
!> \param qs_env Frozen electronic environment with its actual input-cell setting
!> \param k Fractional reciprocal coordinate, in the cell Fourier gauge
!> \param metric Scalar AO overlap at k
!> \param frame Selected S-normalized coefficients, spin-major when soc is true
!> \param energies Selected real band energies in Hartree
!> \param soc Use CP2K's conjugate physical spin representation
!> \param time_reversal Include the grey-group antiunitary operations, checking their sewing algebra
!> \param eps_geometry Geometry tolerance in Bohr
!> \param tolerance Dimensionless metric, subspace and algebra tolerance
!> \param energy_tolerance Absolute eigenvalue-commutator tolerance in Hartree
!> \param max_memory_mb Maximum estimated additional replicated memory in MiB
!> \param iw Summary output unit, negative on nonwriting ranks
!> \param detail_unit Operation, factor, character and multiplicity output unit
!> \param point Optional retained little-group data for subsequent segment compatibility
!> \param reserved_memory_mb Estimated memory already retained for other points
!> \note All operations, including centering translations, are retained. No polar correction,
!>       gauge matching, standard irrep naming, magnetic symmetry search or EBR lookup is performed.
! **************************************************************************************************
   SUBROUTINE write_little_group(qs_env, k, metric, frame, energies, soc, time_reversal, &
                                 eps_geometry, tolerance, energy_tolerance, max_memory_mb, iw, detail_unit, &
                                 point, reserved_memory_mb)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(dp), INTENT(IN)                               :: k(3)
      COMPLEX(dp), INTENT(IN)                            :: metric(:, :), frame(:, :)
      REAL(dp), INTENT(IN)                               :: energies(:)
      LOGICAL, INTENT(IN)                                :: soc, time_reversal
      REAL(dp), INTENT(IN)                               :: eps_geometry, tolerance, &
                                                            energy_tolerance, max_memory_mb
      INTEGER, INTENT(IN)                                :: iw, detail_unit
      TYPE(little_group_data_type), INTENT(OUT), &
         OPTIONAL                                        :: point
      REAL(dp), INTENT(IN), OPTIONAL                     :: reserved_memory_mb

      COMPLEX(dp)                                        :: spin(2, 2)
      COMPLEX(dp), ALLOCATABLE :: action(:, :), band_character(:), characters(:, :), &
         factors(:, :), image(:, :), sewing(:, :, :), spin_lifts(:, :, :), work(:, :)
      INTEGER                                            :: components, d, g, h, i, j, nao, nb, ng, &
                                                            ns, nu, status
      INTEGER, ALLOCATABLE :: corep_counts(:), dimensions(:), mapping(:), multiplicity(:), &
         op_status(:), owners(:), PRODUCT(:, :), rotations(:, :, :)
      LOGICAL, ALLOCATABLE                               :: anti(:)
      REAL(dp)                                           :: character_error, checks(4), &
                                                            energy_error, frame_error, memory, &
                                                            metric_error, product_error, scale
      REAL(dp), ALLOCATABLE                              :: translations(:, :)
      TYPE(cell_type), POINTER                           :: cell
      TYPE(corepresentation_table_type)                  :: coreps
      TYPE(mp_para_env_type), POINTER                    :: para

      CALL get_qs_env(qs_env, cell=cell, para_env=para)
      IF (ANY(cell%perd /= 1)) THEN
         CPABORT('LITTLE_GROUP_IRREPS requires a three-dimensionally periodic cell.')
      END IF
      nao = SIZE(metric, 1)
      nb = SIZE(frame, 2)
      components = 1
      IF (soc) components = 2
      CPASSERT(ALL(SHAPE(metric) == [nao, nao]))
      CPASSERT(SIZE(frame, 1) == components*nao .AND. SIZE(energies) == nb)
      CALL crystal_operations(qs_env, time_reversal, eps_geometry, rotations, translations, anti)
      ns = SIZE(anti)
      ! Bound the factor workspace before constructing the k-dependent subgroup.
      memory = 16.0_dp*(4.0_dp*REAL(nao, dp)**2 + &
                        4.0_dp*components*REAL(nao, dp)*nb + &
                        (REAL(ns, dp) + 6.0_dp)*REAL(nb, dp)**2 + 6.0_dp*REAL(ns, dp)**2 + 4.0_dp*ns)/1024.0_dp**2
      IF (PRESENT(point)) memory = memory + (48.0_dp*REAL(ns, dp)**2 + 64.0_dp*ns)/1024.0_dp**2
      IF (PRESENT(reserved_memory_mb)) memory = memory + reserved_memory_mb
      IF (memory > max_memory_mb) CPABORT('Little-group matrices exceed LITTLE_GROUP_MAX_MEMORY_MB.')
      CALL little_group_factor_system(cell%hmat, rotations, translations, anti, k, soc, &
                                      mapping, product, factors, tolerance, status, spin_lifts=spin_lifts)
      IF (status /= 0) CPABORT('Little-group operation/factor algebra is inconsistent.')
      ng = SIZE(mapping)
      nu = COUNT(.NOT. anti(mapping))
      CPASSERT(ALL(.NOT. anti(mapping(:nu))))
      ALLOCATE (action(nao, nao), image(components*nao, nb), sewing(nb, nb, ng), work(nao, nao))
      ALLOCATE (op_status(ng), owners(ng))
      op_status = 0
      owners = 0
      sewing = 0.0_dp
      scale = MAX(1.0_dp, MAXVAL(ABS(metric)))
      metric_error = 0.0_dp
      frame_error = 0.0_dp
      energy_error = 0.0_dp
      DO g = para%mepos + 1, ng, para%num_pe
         owners(g) = 1
         h = mapping(g)
         CALL gaussian_action(qs_env, k, rotations(:, :, h), translations(:, h), action, eps_geometry, status)
         IF (status /= 0) THEN
            op_status(g) = 1
            CYCLE
         END IF
         work(:, :) = MATMUL(CONJG(TRANSPOSE(action)), MATMUL(metric, action))
         IF (anti(h)) THEN
            work(:, :) = work - CONJG(metric)
         ELSE
            work(:, :) = work - metric
         END IF
         metric_error = MAX(metric_error, MAXVAL(ABS(work))/scale)
         IF (soc) THEN
            ! Reuse the factor-system lift; recomputing a quaternion can change its sign.
            spin(:, :) = CONJG(spin_lifts(:, :, g))
            IF (anti(h)) THEN
               ! Remove the stored Theta: transform_spinor_frame applies it to the source.
               spin(:, 1) = CONJG(spin_lifts(:, 2, g))
               spin(:, 2) = -CONJG(spin_lifts(:, 1, g))
            END IF
            CALL transform_spinor_frame(action, spin, frame, anti(h), image, tolerance, status)
            IF (status /= 0) THEN
               op_status(g) = 2
               CYCLE
            END IF
         ELSE IF (anti(h)) THEN
            image(:, :) = MATMUL(action, CONJG(frame))
         ELSE
            image(:, :) = MATMUL(action, frame)
         END IF
         CALL subspace_sewing(metric, frame, image, sewing(:, :, g), tolerance, checks, status)
         frame_error = MAX(frame_error, MAXVAL(checks))
         IF (status /= 0 .OR. metric_error > tolerance) THEN
            op_status(g) = 3
            CYCLE
         END IF
         DO j = 1, nb
            DO i = 1, nb
               energy_error = MAX(energy_error, ABS((energies(i) - energies(j))*sewing(i, j, g)))
            END DO
         END DO
      END DO
      CALL para%sum(owners)
      CALL para%sum(op_status)
      CALL para%max(metric_error)
      CALL para%max(frame_error)
      CALL para%max(energy_error)
      CPASSERT(ALL(owners == 1))
      IF (ANY(op_status /= 0)) THEN
         g = FINDLOC(op_status /= 0, .TRUE., DIM=1)
         IF (iw > 0) WRITE (iw, '(T2,A,2I8,2ES14.5)') 'TQC| Invalid little-group operation/stage: ', &
            g, op_status(g), metric_error, frame_error
         SELECT CASE (op_status(g))
         CASE (1)
            CPABORT('Little-group atom/basis mapping is not unique or complete.')
         CASE (2)
            CPABORT('Little-group spinor transformation is invalid.')
         CASE DEFAULT
            CPABORT('Little-group operation fails metric covariance or selected-subspace closure.')
         END SELECT
      END IF
      CALL para%sum(sewing)
      CALL check_sewing_product(sewing, product, factors, anti(mapping), tolerance, product_error, status, &
                                [para%mepos, para%num_pe])
      CALL para%min(status)
      CALL para%max(product_error)
      IF (iw > 0) THEN
         WRITE (iw, '(T2,A,I0,A,I0)') 'TQC| Little-group MPI workers: ', MIN(ng, para%num_pe), &
            ' operations: ', ng
         WRITE (iw, '(T2,A,3F12.7,A,2I7)') 'TQC| Little group at ', k, ' unitary/antiunitary: ', nu, ng - nu
         WRITE (iw, '(T2,A,4ES14.5)') 'TQC| Little-group metric/frame/product/energy residuals: ', &
            metric_error, frame_error, product_error, energy_error
      END IF
      IF (status /= 0) CPABORT('Little-group band sewing matrices fail the projective product law.')
      IF (energy_error > energy_tolerance) THEN
         CPABORT('Little-group band symmetry does not commute with the Hamiltonian.')
      END IF
      CALL projective_characters(PRODUCT(:nu, :nu), factors(:nu, :nu), characters, tolerance, character_error, status)
      IF (status /= 0) CPABORT('Little-group unitary character table could not be resolved completely.')
      ALLOCATE (band_character(nu), multiplicity(SIZE(characters, 2)), dimensions(SIZE(characters, 2)))
      dimensions(:) = NINT(REAL(characters(1, :), dp))
      band_character(:) = 0.0_dp
      DO g = 1, nu
         DO i = 1, nb
            band_character(g) = band_character(g) + sewing(i, i, g)
         END DO
      END DO
      CALL character_multiplicities(band_character, characters, multiplicity, tolerance, status)
      IF (status /= 0) CPABORT('Little-group band characters have nonintegral irrep multiplicities.')
      CPASSERT(SUM(dimensions*multiplicity) == nb)
      IF (ng > nu) THEN
         CALL projective_corepresentations(product, factors, anti(mapping), characters, coreps, &
                                           tolerance, character_error, status)
         IF (status /= 0) CPABORT('Little-group antiunitary corepresentation classification failed.')
         ALLOCATE (corep_counts(SIZE(coreps%indicator)))
         CALL corepresentation_multiplicities(coreps, multiplicity, corep_counts, status)
         IF (status /= 0) CPABORT('Little-group bands violate antiunitary irrep pairing or doubling.')
         CPASSERT(SUM(coreps%dimension*corep_counts) == nb)
      END IF
      IF (iw > 0) THEN
         WRITE (iw, '(T2,A,I0)') 'TQC| Little-group represented bands: ', nb
         WRITE (iw, '(T2,A,I0)') 'TQC| Little-group selected irrep types: ', COUNT(multiplicity > 0)
         DO d = 1, MAXVAL(dimensions)
            IF (ANY(dimensions == d)) THEN
               WRITE (iw, '(T2,A,I0,A,I0)') 'TQC| Little-group irrep dimension ', d, &
                  ' multiplicity: ', SUM(multiplicity, MASK=dimensions == d)
            END IF
         END DO
         IF (ng > nu) THEN
            WRITE (iw, '(T2,A,I0)') 'TQC| Corepresentation represented bands: ', SUM(coreps%dimension*corep_counts)
            DO d = -1, 1
               WRITE (iw, '(T2,A,I0,A,I0)') 'TQC| Corepresentation type ', d, &
                  ' multiplicity: ', SUM(corep_counts, MASK=coreps%indicator == d)
            END DO
         END IF
         WRITE (iw, '(T2,A)') 'TQC| Local representations only; no standard irrep labels or EBR classification.'
      END IF
      IF (detail_unit > 0) THEN
         WRITE (detail_unit, '(A,3ES24.16)') 'KPOINT ', k
         WRITE (detail_unit, '(A,4I8)') 'SIZES ', nu, ng - nu, SIZE(characters, 2), nb
         DO g = 1, ng
            h = mapping(g)
            WRITE (detail_unit, '(A,2I8,9I5,3ES24.16)') 'OPERATION ', g, MERGE(1, 0, anti(h)), &
               TRANSPOSE(rotations(:, :, h)), translations(:, h)
         END DO
         DO h = 1, ng
            DO g = 1, ng
               WRITE (detail_unit, '(A,3I8,2ES24.16)') 'PRODUCT ', g, h, PRODUCT(g, h), factors(g, h)
            END DO
         END DO
         DO g = 1, nu
            WRITE (detail_unit, '(A,I8,2ES24.16)') 'BAND_CHARACTER ', g, band_character(g)
         END DO
         DO i = 1, SIZE(characters, 2)
            WRITE (detail_unit, '(A,3I8)') 'IRREP ', i, dimensions(i), multiplicity(i)
            DO g = 1, nu
               WRITE (detail_unit, '(A,2I8,2ES24.16)') 'CHARACTER ', i, g, characters(g, i)
            END DO
         END DO
         IF (ng > nu) THEN
            WRITE (detail_unit, '(A,I8)') 'COREP_SIZES ', SIZE(corep_counts)
            DO j = 1, SIZE(corep_counts)
               WRITE (detail_unit, '(A,4I8)') 'COREP ', j, coreps%dimension(j), coreps%indicator(j), corep_counts(j)
               DO i = 1, SIZE(characters, 2)
                  WRITE (detail_unit, '(A,3I8)') 'COREP_RESTRICTION ', j, i, coreps%restriction(i, j)
               END DO
            END DO
         END IF
         WRITE (detail_unit, '(A)') 'END_KPOINT'
      END IF
      IF (PRESENT(point)) THEN
         point%cell = cell%hmat
         point%k = k
         point%spinful = soc
         ALLOCATE (point%rotations(3, 3, ng), point%translations(3, ng), point%antiunitary(ng))
         ALLOCATE (point%product(ng, ng), point%factors(ng, ng))
         ALLOCATE (point%characters(nu, SIZE(characters, 2)), point%multiplicity(SIZE(multiplicity)))
         point%rotations(:, :, :) = rotations(:, :, mapping)
         point%translations(:, :) = translations(:, mapping)
         point%antiunitary(:) = anti(mapping)
         point%product(:, :) = product
         point%factors(:, :) = factors
         point%characters(:, :) = characters
         point%multiplicity(:) = multiplicity
      END IF
   END SUBROUTINE write_little_group

! **************************************************************************************************
!> \brief Compare validated Gaussian-band signatures with all generated atomic columns at supplied k.
!> \param qs_env Frozen electronic environment in the original cell/origin convention
!> \param points Validated little-group band data; band isolation is checked by the caller
!> \param time_reversal Include the same grey extension as the physical band checks
!> \param eps_geometry SPGLIB tolerance in Bohr
!> \param tolerance Factor and character tolerance
!> \param max_memory_mb Estimated additional replicated matrix-storage limit
!> \param max_nodes Mandatory nonnegative-search node budget, zero allowed as an unresolved test
!> \param iw Summary unit
!> \param detail_unit Detailed matrix, row, column provenance and witness unit
!> \param sources Optional directed graph sources, supplied together with destinations and shifts
!> \param destinations Optional directed graph destinations
!> \param shifts Optional reciprocal destination images
!> \note No complete reciprocal-stratum coverage, EBR elementarity or Bloch-bundle equivalence is implied.
! **************************************************************************************************
   SUBROUTINE write_atomic_signatures(qs_env, points, time_reversal, eps_geometry, tolerance, &
                                      max_memory_mb, max_nodes, iw, detail_unit, sources, destinations, shifts)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(little_group_data_type), INTENT(IN)           :: points(:)
      LOGICAL, INTENT(IN)                                :: time_reversal
      REAL(dp), INTENT(IN)                               :: eps_geometry, tolerance, max_memory_mb
      INTEGER, INTENT(IN)                                :: max_nodes, iw, detail_unit
      INTEGER, INTENT(IN), OPTIONAL                      :: sources(:), destinations(:), shifts(:, :)

      INTEGER                                            :: i, j, m, n, nonnegative_result, p, &
                                                            signed_result, site, status
      INTEGER(int_8)                                     :: visited
      INTEGER(int_8), ALLOCATABLE                        :: matrix(:, :), nonnegative_witness(:), &
                                                            signed_witness(:), TARGET(:)
      INTEGER, ALLOCATABLE                               :: rotations(:, :, :)
      LOGICAL                                            :: exists
      LOGICAL, ALLOCATABLE                               :: anti(:)
      REAL(dp)                                           :: estimate, retained
      REAL(dp), ALLOCATABLE                              :: kpoints(:, :), translations(:, :)
      TYPE(atomic_signature_type)                        :: catalogue
      TYPE(cell_type), POINTER                           :: cell
      TYPE(integer_lattice_type)                         :: factor

      CPASSERT(SIZE(points) > 0 .AND. max_nodes >= 0)
      CALL get_qs_env(qs_env, cell=cell)
      CALL crystal_operations(qs_env, time_reversal, eps_geometry, rotations, translations, anti)
      ALLOCATE (kpoints(3, SIZE(points)))
      retained = 0.0_dp
      DO p = 1, SIZE(points)
         kpoints(:, p) = points(p)%k
         retained = retained + little_group_storage(points(p))
      END DO
      IF (retained >= max_memory_mb) CPABORT('Atomic signatures exceed LITTLE_GROUP_MAX_MEMORY_MB.')
      CALL atomic_signature_construct(cell%hmat, rotations, translations, anti, kpoints, points(1)%spinful, &
                                      catalogue, tolerance, status, max_memory_mb - retained)
      IF (status /= 0) THEN
         CPABORT('Atomic signature construction failed or exceeds LITTLE_GROUP_MAX_MEMORY_MB.')
      END IF
      CALL atomic_signature_from_bands(catalogue, points, TARGET, tolerance, status)
      IF (status /= 0) THEN
         CPABORT('Atomic references and Gaussian bands have inconsistent symmetry conventions or ranks.')
      END IF
      m = SIZE(catalogue%matrix, 1)
      n = SIZE(catalogue%matrix, 2)
      retained = retained + star_group_storage(catalogue%group)
      retained = retained + 16.0_dp*SIZE(catalogue%reference_spin)/1024.0_dp**2
      retained = retained + 4.0_dp*SIZE(catalogue%character_class)/1024.0_dp**2
      DO i = 1, SIZE(catalogue%sites)
         retained = retained + atomic_site_storage(catalogue%sites(i))
      END DO
      retained = retained + (8.0_dp*REAL(n, dp)**2 + 512.0_dp*SIZE(catalogue%site_relations))/1024.0_dp**2
      DO i = 1, SIZE(catalogue%site_relations)
         retained = retained + 8.0_dp*SIZE(catalogue%site_relations(i)%coefficients)/1024.0_dp**2
      END DO
      DO p = 1, SIZE(catalogue%points)
         retained = retained + 2*little_group_storage(catalogue%points(p)%little)
      END DO
      CPASSERT(PRESENT(sources) .EQV. PRESENT(destinations))
      CPASSERT(PRESENT(sources) .EQV. PRESENT(shifts))
      IF (PRESENT(sources)) THEN
         CALL write_compatibility_lattice(catalogue, TARGET, sources, destinations, shifts, tolerance, &
                                          max_memory_mb - retained - 16.0_dp*REAL(m, dp)*n/1024.0_dp**2, iw, detail_unit)
      END IF
      ! Include integer factorization, certificate/witness work and both signature copies.
      estimate = retained + 32.0_dp*(REAL(m, dp)**2 + REAL(n, dp)**2 + 3.0_dp*m*n)/1024.0_dp**2
      IF (estimate > max_memory_mb) CPABORT('Atomic integer matrices exceed LITTLE_GROUP_MAX_MEMORY_MB.')
      matrix = INT(catalogue%matrix, int_8)
      signed_result = -1
      nonnegative_result = -1
      visited = 0
      CALL integer_lattice_factor(matrix, factor, status)
      IF (status /= 0 .AND. status /= -3) THEN
         CPABORT('Atomic integer factorization is invalid.')
      END IF
      IF (status == 0) THEN
         CALL integer_lattice_solve(factor, TARGET, exists, signed_witness, status)
         IF (status /= 0 .AND. status /= -3) THEN
            CPABORT('Atomic signed witness validation failed.')
         END IF
         IF (status == 0) signed_result = MERGE(1, 0, exists)
      END IF
      IF (signed_result == 0) THEN
         nonnegative_result = 0
      ELSE
         CALL nonnegative_integer_solve(matrix, TARGET, exists, nonnegative_witness, visited, status, INT(max_nodes, int_8))
         IF (status /= 0 .AND. status /= 1 .AND. status /= -3) THEN
            CPABORT('Atomic nonnegative witness validation failed.')
         END IF
         IF (status == 0) nonnegative_result = MERGE(1, 0, exists)
         IF (nonnegative_result == 1) signed_result = 1
      END IF
      IF (iw > 0) THEN
         WRITE (iw, '(T2,A,3I8)') 'TQC| Atomic signature points/rows/columns: ', SIZE(points), m, n
         WRITE (iw, '(T2,A,I0)') 'TQC| Full-torus unitary-character classes: ', &
            COUNT(catalogue%character_class == [(i, i=1, n)])
         WRITE (iw, '(T2,A)') 'TQC| Character equality does not establish Bloch-bundle equivalence.'
         WRITE (iw, '(T2,A,I0)') 'TQC| Sampled signed atomic membership (1/0/-1): ', signed_result
         WRITE (iw, '(T2,A,I0)') 'TQC| Sampled nonnegative atomic membership (1/0/-1): ', nonnegative_result
         WRITE (iw, '(T2,A,I0)') 'TQC| Atomic signature search nodes: ', visited
         WRITE (iw, '(T2,A)') 'TQC| Membership means yes/no/unresolved for these sampled signatures only.'
         WRITE (iw, '(T2,A)') 'TQC| No complete EBR classification or proof of band-bundle triviality; check BZ gaps.'
      END IF
      CALL write_site_induction(catalogue, iw, detail_unit)
      CALL write_atomic_onsite(catalogue, detail_unit)
      IF (detail_unit <= 0) RETURN
      WRITE (detail_unit, '(A,3I8)') 'ATOMIC_SIGNATURE_SIZES ', SIZE(points), m, n
      WRITE (detail_unit, '(A,ES24.16)') 'ATOMIC_CHARACTER_CLASS_RESIDUAL ', catalogue%character_residual
      WRITE (detail_unit, '(A,9ES24.16)') 'ATOMIC_REFERENCE_CELL ', cell%hmat
      DO i = 1, SIZE(anti)
         WRITE (detail_unit, '(A,2I8,9I5,3ES24.16)') 'ATOMIC_REFERENCE_OPERATION ', i, MERGE(1, 0, anti(i)), &
            TRANSPOSE(rotations(:, :, i)), translations(:, i)
         WRITE (detail_unit, '(A,I8,8ES24.16)') 'ATOMIC_REFERENCE_SPIN ', i, catalogue%reference_spin(:, :, i)
      END DO
      DO p = 1, SIZE(points)
         DO i = catalogue%row_begin(p), catalogue%row_begin(p + 1) - 1
            WRITE (detail_unit, '(A,4I8,I20)') 'ATOMIC_SIGNATURE_ROW ', i, p, &
               i - catalogue%row_begin(p) + 1, catalogue%row_dimension(i), TARGET(i)
         END DO
         DO j = 1, SIZE(catalogue%points(p)%little%characters, 2)
            DO i = 1, SIZE(catalogue%points(p)%mapping)
               WRITE (detail_unit, '(A,3I8,2ES24.16)') 'ATOMIC_REFERENCE_CHARACTER ', p, j, &
                  catalogue%points(p)%mapping(i), catalogue%points(p)%little%characters(i, j)
            END DO
         END DO
         IF (ALLOCATED(catalogue%points(p)%coreps%restriction)) THEN
            DO j = 1, SIZE(catalogue%points(p)%coreps%restriction, 2)
               DO i = 1, SIZE(catalogue%points(p)%coreps%restriction, 1)
                  WRITE (detail_unit, '(A,4I8)') 'ATOMIC_REFERENCE_COREP ', p, j, i, &
                     catalogue%points(p)%coreps%restriction(i, j)
               END DO
            END DO
         END IF
      END DO
      DO j = 1, n
         site = catalogue%column_site(j)
         WRITE (detail_unit, '(A,4I8,3ES24.16)') 'ATOMIC_SIGNATURE_COLUMN ', j, site, catalogue%column_local(j), &
            catalogue%band_dimension(j), catalogue%families(site)%position
         WRITE (detail_unit, '(A,2I8)') 'ATOMIC_CHARACTER_CLASS ', j, catalogue%character_class(j)
         DO i = 1, m
            WRITE (detail_unit, '(A,3I8)') 'ATOMIC_SIGNATURE_ENTRY ', i, j, catalogue%matrix(i, j)
         END DO
         IF (ALLOCATED(signed_witness)) WRITE (detail_unit, '(A,I8,I20)') 'SIGNED_ATOMIC_WITNESS ', j, signed_witness(j)
         IF (ALLOCATED(nonnegative_witness)) WRITE (detail_unit, '(A,I8,I20)') &
            'NONNEGATIVE_ATOMIC_WITNESS ', j, nonnegative_witness(j)
      END DO
      WRITE (detail_unit, '(A,2I8,I20)') 'ATOMIC_SIGNATURE_RESULT ', signed_result, nonnegative_result, visited
      WRITE (detail_unit, '(A)') 'END_ATOMIC_SIGNATURES'
   END SUBROUTINE write_atomic_signatures

! **************************************************************************************************
!> \brief Preserve local characters and physical conventions for every atomic signature column.
!> \param catalogue Validated atomic signature catalogue
!> \param unit Detailed output unit, nonpositive on nonwriting ranks
! **************************************************************************************************
   SUBROUTINE write_atomic_onsite(catalogue, unit)
      TYPE(atomic_signature_type), INTENT(IN)            :: catalogue
      INTEGER, INTENT(IN)                                :: unit

      INTEGER                                            :: g, j, row, s

      IF (unit <= 0) RETURN
      DO s = 1, SIZE(catalogue%sites)
         ASSOCIATE (site => catalogue%sites(s))
            ! Onsite factors at Gamma are exactly the physical double-group signs.
            CPASSERT(ALL(AIMAG(site%factor) == 0.0_dp))
            CPASSERT(ALL(ABS(REAL(site%factor, dp)) == 1.0_dp))
            WRITE (unit, '(A,5I8)') 'ATOMIC_ONSITE_SIZES ', s, SIZE(site%mapping), &
               SIZE(site%characters, 1), SIZE(site%characters, 2), MERGE(1, 0, site%spinful)
            row = 0
            DO g = 1, SIZE(site%mapping)
               WRITE (unit, '(A,4I8,12ES24.16)') 'ATOMIC_ONSITE_OPERATION ', s, g, site%mapping(g), &
                  MERGE(1, 0, site%antiunitary(g)), TRANSPOSE(site%cartesian(:, :, g)), site%lattice_shifts(:, g)
               WRITE (unit, '(A,2I8,8ES24.16)') 'ATOMIC_ONSITE_SPIN ', s, g, site%spin(:, :, g)
               WRITE (unit, '(A,2I8,*(I6))') 'ATOMIC_ONSITE_PRODUCT_ROW ', s, g, &
                  site%product(g, :)*NINT(REAL(site%factor(g, :), dp))
               IF (site%antiunitary(g)) CYCLE
               row = row + 1
               DO j = 1, SIZE(site%characters, 2)
                  WRITE (unit, '(A,3I8,2ES24.16)') 'ATOMIC_ONSITE_CHARACTER ', s, j, g, site%characters(row, j)
               END DO
            END DO
         END ASSOCIATE
      END DO
      WRITE (unit, '(A)') 'END_ATOMIC_ONSITE'
   END SUBROUTINE write_atomic_onsite

! **************************************************************************************************
!> \brief Export real-space induction witnesses separately from sampled reciprocal membership.
!> \param catalogue Validated atomic columns and induction relations
!> \param iw Main output unit, or nonpositive to suppress
!> \param detail_unit Machine-readable output unit, or nonpositive to suppress
! **************************************************************************************************
   SUBROUTINE write_site_induction(catalogue, iw, detail_unit)
      TYPE(atomic_signature_type), INTENT(IN)            :: catalogue
      INTEGER, INTENT(IN)                                :: iw, detail_unit

      INTEGER                                            :: e, i, j, nc, ncomposite

      nc = SIZE(catalogue%column_site)
      ncomposite = COUNT(SUM(catalogue%site_reduction%expansion, DIM=1) > 1)
      IF (iw > 0) THEN
         WRITE (iw, '(T2,A,2I8)') 'TQC| Real-space induction relations/composite columns: ', &
            SIZE(catalogue%site_relations), ncomposite
         WRITE (iw, '(T2,A)') 'TQC| These decompositions have site-homotopy witnesses, not only matching k labels.'
      END IF
      IF (detail_unit <= 0) RETURN
      WRITE (detail_unit, '(A,3I8)') 'SITE_INDUCTION_SIZES ', SIZE(catalogue%site_relations), nc, ncomposite
      DO e = 1, SIZE(catalogue%site_relations)
         WRITE (detail_unit, '(A,7I8,3ES24.16)') 'SITE_INDUCTION_EDGE ', e, catalogue%site_relations(e)%source, &
            catalogue%site_relations(e)%destination, catalogue%site_relations(e)%operation, &
            catalogue%site_relations(e)%shift, catalogue%site_relations(e)%image
         DO j = 1, SIZE(catalogue%site_relations(e)%coefficients, 2)
            DO i = 1, SIZE(catalogue%site_relations(e)%coefficients, 1)
               IF (catalogue%site_relations(e)%coefficients(i, j) == 0) CYCLE
               WRITE (detail_unit, '(A,4I8)') 'SITE_INDUCTION_ENTRY ', e, i, j, &
                  catalogue%site_relations(e)%coefficients(i, j)
            END DO
         END DO
      END DO
      DO j = 1, nc
         WRITE (detail_unit, '(A,4I8)') 'SITE_INDUCTION_COLUMN ', j, catalogue%site_reduction%representative(j), &
            catalogue%site_reduction%witness_edge(j), catalogue%site_reduction%witness_local(j)
         DO i = 1, nc
            IF (catalogue%site_reduction%expansion(i, j) == 0) CYCLE
            WRITE (detail_unit, '(A,3I8)') 'SITE_INDUCTION_EXPANSION ', i, j, catalogue%site_reduction%expansion(i, j)
         END DO
      END DO
      WRITE (detail_unit, '(A)') 'END_SITE_INDUCTION'
   END SUBROUTINE write_site_induction

! **************************************************************************************************
!> \brief Report a certified integer compatibility quotient on the supplied property graph.
!> \param catalogue Native reference rows/columns
!> \param signature Validated Gaussian multiplicities in the catalogue row basis
!> \param sources Directed source points
!> \param destinations Directed destination points
!> \param shifts Reciprocal destination images
!> \param tolerance Character and fractional geometry tolerance
!> \param available_memory_mb Additional storage left after retained point/reference data
!> \param iw Summary output
!> \param unit Detailed output
! **************************************************************************************************
   SUBROUTINE write_compatibility_lattice(catalogue, signature, sources, destinations, shifts, tolerance, &
                                          available_memory_mb, iw, unit)
      TYPE(atomic_signature_type), INTENT(IN)            :: catalogue
      INTEGER(int_8), INTENT(IN)                         :: signature(:)
      INTEGER, INTENT(IN)                                :: sources(:), destinations(:), shifts(:, :)
      REAL(dp), INTENT(IN)                               :: tolerance, available_memory_mb
      INTEGER, INTENT(IN)                                :: iw, unit

      INTEGER                                            :: decision, i, j, q, status
      INTEGER(int_8), ALLOCATABLE                        :: indices(:)
      LOGICAL                                            :: compatible
      TYPE(compatibility_lattice_type)                   :: lattice

      status = -3
      IF (available_memory_mb > 0.0_dp) THEN
         CALL compatibility_lattice_construct(catalogue, sources, destinations, shifts, lattice, &
                                              tolerance, status, available_memory_mb)
      END IF
      decision = -1
      IF (status == 0) THEN
         CALL integer_quotient_class(lattice%quotient, signature, compatible, indices, status)
         IF (status == 0) decision = MERGE(1, 0, compatible)
      END IF
      IF (status /= 0 .AND. status /= -3) THEN
         CPABORT('Joint compatibility lattice or atomic quotient certificate failed.')
      END IF
      IF (iw > 0) THEN
         WRITE (iw, '(T2,A,I0)') 'TQC| Sampled compatibility-lattice membership (1/0/-1): ', decision
         WRITE (iw, '(T2,A)') 'TQC| This quotient concerns the supplied graph, not complete BZ/EBR classification.'
      END IF
      IF (unit > 0) WRITE (unit, '(A,I8)') 'COMPATIBILITY_LATTICE_STATUS ', status
      IF (status == 0) THEN
         q = SIZE(lattice%quotient%invariants)
         IF (iw > 0) THEN
            WRITE (iw, '(T2,A,I0)') 'TQC| Sampled star identifications: ', SIZE(lattice%stars, 2)
            WRITE (iw, '(T2,A,I0)') 'TQC| Sampled compatible signature rank: ', q
            WRITE (iw, '(T2,A,I0)') 'TQC| Sampled quotient free rank: ', COUNT(lattice%quotient%invariants == 0)
            WRITE (iw, '(T2,A,I0)') 'TQC| Sampled quotient cyclic factors: ', COUNT(lattice%quotient%invariants > 1)
            IF (compatible) WRITE (iw, '(T2,A,I0)') 'TQC| Sampled quotient nonzero class: ', MERGE(1, 0, ANY(indices /= 0))
         END IF
         IF (unit > 0) THEN
            WRITE (unit, '(A,3I8)') 'COMPATIBILITY_LATTICE_SIZES ', SHAPE(lattice%constraints), q
            DO j = 1, SIZE(lattice%stars, 2)
               WRITE (unit, '(A,6I8)') 'COMPATIBILITY_STAR ', lattice%stars(:, j)
            END DO
            DO j = 1, SIZE(lattice%constraints, 2)
               DO i = 1, SIZE(lattice%constraints, 1)
                  IF (lattice%constraints(i, j) == 0) CYCLE
                  WRITE (unit, '(A,2I8,I20)') 'COMPATIBILITY_CONSTRAINT ', i, j, lattice%constraints(i, j)
               END DO
            END DO
            DO j = 1, q
               WRITE (unit, '(A,I8,I20)') 'COMPATIBILITY_QUOTIENT_FACTOR ', j, lattice%quotient%invariants(j)
               IF (compatible) WRITE (unit, '(A,I8,I20)') 'COMPATIBILITY_QUOTIENT_CLASS ', j, indices(j)
               DO i = 1, SIZE(lattice%quotient%kernel, 1)
                  IF (lattice%quotient%kernel(i, j) == 0) CYCLE
                  WRITE (unit, '(A,2I8,I20)') 'COMPATIBILITY_KERNEL ', i, j, lattice%quotient%kernel(i, j)
               END DO
               DO i = 1, q
                  IF (lattice%quotient%left(i, j) == 0) CYCLE
                  WRITE (unit, '(A,2I8,I20)') 'COMPATIBILITY_QUOTIENT_MAP ', i, j, lattice%quotient%left(i, j)
               END DO
            END DO
            DO j = 1, SIZE(lattice%quotient%coordinates, 2)
               DO i = 1, q
                  IF (lattice%quotient%coordinates(i, j) == 0) CYCLE
                  WRITE (unit, '(A,2I8,I20)') 'COMPATIBILITY_ATOMIC_COORDINATE ', i, j, lattice%quotient%coordinates(i, j)
               END DO
            END DO
         END IF
      END IF
      IF (unit > 0) WRITE (unit, '(A)') 'END_COMPATIBILITY_LATTICE'
   END SUBROUTINE write_compatibility_lattice

! **************************************************************************************************
!> \brief Resolve the complete spatial group in the input cell, optionally with its grey extension.
!> \param qs_env Electronic environment
!> \param time_reversal Include antiunitary copies
!> \param eps_geometry SPGLIB tolerance in Bohr
!> \param rotations Direct-space integer rotations
!> \param translations Fractional coset representatives, including centering
!> \param anti Antiunitary flags
! **************************************************************************************************
   SUBROUTINE crystal_operations(qs_env, time_reversal, eps_geometry, rotations, translations, anti)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: time_reversal
      REAL(dp), INTENT(IN)                               :: eps_geometry
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: rotations(:, :, :)
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: translations(:, :)
      LOGICAL, ALLOCATABLE, INTENT(OUT)                  :: anti(:)

      INTEGER                                            :: a, g, ns
      INTEGER, ALLOCATABLE                               :: kinds_set(:)
      REAL(dp), ALLOCATABLE                              :: positions(:, :)
      TYPE(cell_type), POINTER                           :: cell
      TYPE(csym_type)                                    :: crystal
      TYPE(particle_type), POINTER                       :: particles(:)

      CALL get_qs_env(qs_env, cell=cell, particle_set=particles)
      ALLOCATE (positions(3, SIZE(particles)), kinds_set(SIZE(particles)))
      DO a = 1, SIZE(particles)
         positions(:, a) = MATMUL(cell%h_inv, pbc(particles(a)%r, cell))
         CALL get_atomic_kind(particles(a)%atomic_kind, kind_number=kinds_set(a))
      END DO
      CALL crys_sym_gen(crystal, positions, kinds_set, cell%hmat, delta=eps_geometry, use_spglib=.TRUE.)
      IF (.NOT. crystal%symlib) THEN
         CPABORT('LITTLE_GROUP_IRREPS requires SPGLIB and a valid symmetry search.')
      END IF
      ns = crystal%n_operations
      IF (time_reversal) ns = 2*ns
      ALLOCATE (rotations(3, 3, ns), translations(3, ns), anti(ns))
      DO g = 1, crystal%n_operations
         rotations(:, :, g) = TRANSPOSE(crystal%rotations(:, :, g))
         translations(:, g) = crystal%translations(:, g)
      END DO
      anti = .FALSE.
      IF (time_reversal) THEN
         g = crystal%n_operations
         rotations(:, :, g + 1:) = rotations(:, :, :g)
         translations(:, g + 1:) = translations(:, :g)
         anti(g + 1:) = .TRUE.
      END IF
      CALL release_csym_type(crystal)
   END SUBROUTINE crystal_operations

! **************************************************************************************************
!> \brief Compare unitary and antiunitary endpoint restrictions on explicit NNKP/Wilson segments.
!> \param points Validated band data for all explicit property points
!> \param neighbours Target-point indices, with the source point as the first dimension
!> \param shifts Reciprocal target-cell shifts for each directed link
!> \param tolerance Dimensionless factor and character tolerance
!> \param max_memory_mb Additional replicated-storage limit
!> \param iw Summary output unit
!> \param detail_unit Detailed branching output unit, in the same file as the point data
!> \param para Complete property MPI environment
! **************************************************************************************************
   SUBROUTINE write_little_group_edges(points, neighbours, shifts, tolerance, max_memory_mb, iw, detail_unit, para)
      TYPE(little_group_data_type), INTENT(IN)           :: points(:)
      INTEGER, INTENT(IN)                                :: neighbours(:, :), shifts(:, :, :)
      REAL(dp), INTENT(IN)                               :: tolerance, max_memory_mb
      INTEGER, INTENT(IN)                                :: iw, detail_unit
      TYPE(mp_para_env_type), POINTER                    :: para

      INTEGER                                            :: a, e, link
      INTEGER, ALLOCATABLE                               :: destinations(:), images(:, :), sources(:)

      CPASSERT(SIZE(neighbours, 1) == SIZE(points))
      CPASSERT(ALL(SHAPE(shifts) == [3, SIZE(points), SIZE(neighbours, 2)]))
      ALLOCATE (sources(SIZE(neighbours)), destinations(SIZE(neighbours)), images(3, SIZE(neighbours)))
      e = 0
      DO link = 1, SIZE(neighbours, 2)
         DO a = 1, SIZE(points)
            e = e + 1
            sources(e) = a
            destinations(e) = neighbours(a, link)
            images(:, e) = shifts(:, a, link)
         END DO
      END DO
      CALL write_little_group_graph(points, sources, destinations, images, tolerance, max_memory_mb, iw, detail_unit, para)
   END SUBROUTINE write_little_group_edges

! **************************************************************************************************
!> \brief Compare band restrictions on a directed graph without rectangular neighbour padding.
!> \param points Validated property-point band data
!> \param sources Source-point indices
!> \param destinations Destination-point indices
!> \param shifts Destination reciprocal-cell shifts
!> \param tolerance Factor/character tolerance
!> \param max_memory_mb Additional replicated-storage limit
!> \param iw Summary output unit
!> \param detail_unit Detailed branching output unit
!> \param para Complete property MPI environment; each rank owns at most one edge per batch
! **************************************************************************************************
   SUBROUTINE write_little_group_graph(points, sources, destinations, shifts, tolerance, max_memory_mb, iw, detail_unit, para)
      TYPE(little_group_data_type), INTENT(IN)           :: points(:)
      INTEGER, INTENT(IN)                                :: sources(:), destinations(:), shifts(:, :)
      REAL(dp), INTENT(IN)                               :: tolerance, max_memory_mb
      INTEGER, INTENT(IN)                                :: iw, detail_unit
      TYPE(mp_para_env_type), POINTER                    :: para

      INTEGER                                            :: a, b, corep_segments, count, failure, &
                                                            first, g, i, j, last, link, &
                                                            mismatches, n, owner, status
      INTEGER, ALLOCATABLE                               :: expected(:), owners(:)
      LOGICAL                                            :: incompatible, incompatible_coreps
      REAL(dp)                                           :: estimate, residual, retained
      TYPE(little_group_edge_type)                       :: edge, local_edge

      CPASSERT(SIZE(sources) == SIZE(destinations))
      CPASSERT(ALL(SHAPE(shifts) == [3, SIZE(sources)]))
      retained = 0.0_dp
      DO a = 1, SIZE(points)
         retained = retained + little_group_storage(points(a))
      END DO
      count = 0
      corep_segments = 0
      mismatches = 0
      ALLOCATE (owners(para%num_pe), expected(para%num_pe))
      DO first = 1, SIZE(sources), para%num_pe
         last = MIN(first + para%num_pe - 1, SIZE(sources))
         local_edge = little_group_edge_type()
         edge = little_group_edge_type()
         owners = 0
         expected = 0
         n = 0
         DO link = first, last
            a = sources(link)
            b = destinations(link)
            CPASSERT(a >= 1 .AND. a <= SIZE(points))
            CPASSERT(b >= 1 .AND. b <= SIZE(points))
            IF (MAXVAL(ABS(points(b)%k + REAL(shifts(:, link), dp) - points(a)%k)) <= tolerance) CYCLE
            expected(link - first + 1) = 1
            n = MAX(n, SIZE(points(a)%product, 1), SIZE(points(b)%product, 1))
         END DO
         IF (n == 0) CYCLE
         ! Include one owned edge, one received edge and the algebra workspace, not the whole graph.
         estimate = retained + 320.0_dp*REAL(n, dp)**2/1024.0_dp**2
         IF (estimate > max_memory_mb) THEN
            CPABORT('Compatibility data exceed LITTLE_GROUP_MAX_MEMORY_MB.')
         END IF
         failure = HUGE(0)
         status = 0
         residual = 0.0_dp
         link = first + para%mepos
         IF (expected(para%mepos + 1) == 1) THEN
            a = sources(link)
            b = destinations(link)
            CALL little_group_compatibility(points(a), points(b), shifts(:, link), local_edge, tolerance, residual, status)
            owners(para%mepos + 1) = 1
            IF (status /= 0) failure = link
         END IF
         CALL para%sum(owners)
         CPASSERT(ALL(owners == expected))
         CALL para%min(failure)
         IF (failure /= HUGE(0)) THEN
            CALL para%bcast(status, failure - first)
            CALL para%bcast(residual, failure - first)
            a = sources(failure)
            b = destinations(failure)
            IF (iw > 0) WRITE (iw, '(T2,A,3I8,ES14.5)') 'TQC| Invalid compatibility edge/status: ', a, b, status, residual
            CPABORT('Little-group compatibility algebra failed.')
         END IF
         DO link = first, last
            owner = link - first
            IF (expected(owner + 1) == 0) CYCLE
            a = sources(link)
            b = destinations(link)
            CALL broadcast_compatibility_edge(local_edge, edge, owner, para)
            count = count + 1
            incompatible = ANY(edge%left_count /= edge%right_count)
            IF (incompatible) mismatches = mismatches + 1
            IF (ALLOCATED(edge%anti_left)) THEN
               corep_segments = corep_segments + 1
               incompatible_coreps = ANY(edge%left_corep_count /= edge%right_corep_count)
               CPASSERT(incompatible .EQV. incompatible_coreps)
            END IF
            IF (detail_unit <= 0) CYCLE
            WRITE (detail_unit, '(A,5I8)') 'EDGE ', a, b, shifts(:, link)
            DO g = 1, SIZE(edge%left)
               WRITE (detail_unit, '(A,3I8)') 'LINE_MEMBER ', g, edge%left(g), edge%right(g)
            END DO
            DO i = 1, SIZE(edge%characters, 2)
               WRITE (detail_unit, '(A,3I8)') 'LINE_MULTIPLICITY ', i, edge%left_count(i), edge%right_count(i)
               DO g = 1, SIZE(edge%left)
                  WRITE (detail_unit, '(A,2I8,2ES24.16)') 'LINE_CHARACTER ', i, g, edge%characters(g, i)
               END DO
               DO j = 1, SIZE(edge%left_branch, 2)
                  WRITE (detail_unit, '(A,3I8)') 'LEFT_RESTRICTION ', j, i, edge%left_branch(i, j)
               END DO
               DO j = 1, SIZE(edge%right_branch, 2)
                  WRITE (detail_unit, '(A,3I8)') 'RIGHT_RESTRICTION ', j, i, edge%right_branch(i, j)
               END DO
            END DO
            IF (ALLOCATED(edge%anti_left)) THEN
               DO g = 1, SIZE(edge%anti_left)
                  WRITE (detail_unit, '(A,3I8)') 'LINE_ANTIUNITARY_MEMBER ', g, edge%anti_left(g), edge%anti_right(g)
               END DO
               DO i = 1, SIZE(edge%left_corep_count)
                  WRITE (detail_unit, '(A,5I8)') 'LINE_COREP ', i, edge%coreps%dimension(i), edge%coreps%indicator(i), &
                     edge%left_corep_count(i), edge%right_corep_count(i)
                  DO j = 1, SIZE(edge%coreps%restriction, 1)
                     WRITE (detail_unit, '(A,3I8)') 'LINE_COREP_RESTRICTION ', i, j, edge%coreps%restriction(j, i)
                  END DO
               END DO
            END IF
            IF (ALLOCATED(edge%left_corep_branch)) THEN
               DO j = 1, SIZE(edge%left_corep_branch, 2)
                  DO i = 1, SIZE(edge%left_corep_branch, 1)
                     WRITE (detail_unit, '(A,3I8)') 'LEFT_COREP_RESTRICTION ', j, i, edge%left_corep_branch(i, j)
                  END DO
               END DO
            END IF
            IF (ALLOCATED(edge%right_corep_branch)) THEN
               DO j = 1, SIZE(edge%right_corep_branch, 2)
                  DO i = 1, SIZE(edge%right_corep_branch, 1)
                     WRITE (detail_unit, '(A,3I8)') 'RIGHT_COREP_RESTRICTION ', j, i, edge%right_corep_branch(i, j)
                  END DO
               END DO
            END IF
            WRITE (detail_unit, '(A)') 'END_EDGE'
         END DO
      END DO
      IF (iw > 0) THEN
         WRITE (iw, '(T2,A,I0)') 'TQC| Compatibility MPI workers: ', MIN(SIZE(sources), para%num_pe)
         WRITE (iw, '(T2,A,I0)') 'TQC| Compatibility segments: ', count
         WRITE (iw, '(T2,A,I0)') 'TQC| Antiunitary compatibility segments: ', corep_segments
         WRITE (iw, '(T2,A,I0)') 'TQC| Incompatible segments: ', mismatches
         IF (mismatches > 0) WRITE (iw, '(T2,A)') 'TQC| Incompatible endpoint counts: check gaps, symmetry and band selection.'
         WRITE (iw, '(T2,A)') 'TQC| Endpoint compatibility does not prove a gap or trivial topology.'
      END IF
   END SUBROUTINE write_little_group_graph

! **************************************************************************************************
!> \brief Deliver one validated edge without retaining the batch on every rank.
!> \param local_edge Owned edge on its computing rank, ignored on other ranks
!> \param edge Collected edge, replacing the previous delivery
!> \param source Rank owning this edge
!> \param para Complete property MPI environment
! **************************************************************************************************
   SUBROUTINE broadcast_compatibility_edge(local_edge, edge, source, para)
      TYPE(little_group_edge_type), INTENT(IN)           :: local_edge
      TYPE(little_group_edge_type), INTENT(OUT)          :: edge
      INTEGER, INTENT(IN)                                :: source
      TYPE(mp_para_env_type), POINTER                    :: para

      INTEGER                                            :: shape_info(10)

      shape_info = 0
      IF (para%mepos == source) THEN
         edge = local_edge
         shape_info(:4) = [SIZE(edge%left), SIZE(edge%characters, 2), &
                           SIZE(edge%left_branch, 2), SIZE(edge%right_branch, 2)]
         IF (ALLOCATED(edge%anti_left)) shape_info(5) = SIZE(edge%anti_left)
         IF (ALLOCATED(edge%coreps%indicator)) shape_info(6) = SIZE(edge%coreps%indicator)
         IF (ALLOCATED(edge%left_corep_branch)) shape_info(7:8) = SHAPE(edge%left_corep_branch)
         IF (ALLOCATED(edge%right_corep_branch)) shape_info(9:10) = SHAPE(edge%right_corep_branch)
      END IF
      CALL para%bcast(shape_info, source)
      CPASSERT(ALL(shape_info >= 0) .AND. ALL(shape_info(:4) > 0))
      IF (para%mepos /= source) THEN
         ALLOCATE (edge%left(shape_info(1)), edge%right(shape_info(1)))
         ALLOCATE (edge%characters(shape_info(1), shape_info(2)))
         ALLOCATE (edge%left_branch(shape_info(2), shape_info(3)))
         ALLOCATE (edge%right_branch(shape_info(2), shape_info(4)))
         ALLOCATE (edge%left_count(shape_info(2)), edge%right_count(shape_info(2)))
         IF (shape_info(5) > 0) ALLOCATE (edge%anti_left(shape_info(5)), edge%anti_right(shape_info(5)))
         IF (shape_info(6) > 0) THEN
            ALLOCATE (edge%coreps%restriction(shape_info(2), shape_info(6)))
            ALLOCATE (edge%coreps%indicator(shape_info(6)), edge%coreps%dimension(shape_info(6)))
         END IF
         IF (shape_info(7) > 0) THEN
            ALLOCATE (edge%left_corep_branch(shape_info(7), shape_info(8)), edge%left_corep_count(shape_info(7)))
         END IF
         IF (shape_info(9) > 0) THEN
            ALLOCATE (edge%right_corep_branch(shape_info(9), shape_info(10)), edge%right_corep_count(shape_info(9)))
         END IF
      END IF
      CALL para%bcast(edge%left, source)
      CALL para%bcast(edge%right, source)
      CALL para%bcast(edge%characters, source)
      CALL para%bcast(edge%left_branch, source)
      CALL para%bcast(edge%right_branch, source)
      CALL para%bcast(edge%left_count, source)
      CALL para%bcast(edge%right_count, source)
      IF (shape_info(5) > 0) THEN
         CALL para%bcast(edge%anti_left, source)
         CALL para%bcast(edge%anti_right, source)
      END IF
      IF (shape_info(6) > 0) THEN
         CALL para%bcast(edge%coreps%restriction, source)
         CALL para%bcast(edge%coreps%indicator, source)
         CALL para%bcast(edge%coreps%dimension, source)
      END IF
      IF (shape_info(7) > 0) THEN
         CALL para%bcast(edge%left_corep_branch, source)
         CALL para%bcast(edge%left_corep_count, source)
      END IF
      IF (shape_info(9) > 0) THEN
         CALL para%bcast(edge%right_corep_branch, source)
         CALL para%bcast(edge%right_corep_count, source)
      END IF
   END SUBROUTINE broadcast_compatibility_edge

! **************************************************************************************************
!> \brief Gaussian coefficient action exp(-2*pi*i*k.L) times the existing CP2K angular rotation.
!> \param qs_env Electronic environment with Gaussian bases and cell-gauge atomic representatives
!> \param k Target fractional reciprocal vector
!> \param rotation Integer fractional direct-space rotation
!> \param translation Fractional direct-space translation
!> \param action Dense scalar AO action, target rows and source columns
!> \param tolerance Atomic matching distance in Bohr
!> \param status Zero on success; negative for an invalid or ambiguous atom/basis mapping
! **************************************************************************************************
   SUBROUTINE gaussian_action(qs_env, k, rotation, translation, action, tolerance, status)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(dp), INTENT(IN)                               :: k(3)
      INTEGER, INTENT(IN)                                :: rotation(3, 3)
      REAL(dp), INTENT(IN)                               :: translation(3)
      COMPLEX(dp), INTENT(OUT)                           :: action(:, :)
      REAL(dp), INTENT(IN)                               :: tolerance
      INTEGER, INTENT(OUT)                               :: status

      INTEGER                                            :: a, b, found, i, maxl, n, na, nao
      INTEGER, ALLOCATABLE                               :: kinds(:), offset(:), seen(:), &
                                                            shift(:, :), TARGET(:)
      REAL(dp)                                           :: angle, cart(3, 3), delta(3)
      REAL(dp), ALLOCATABLE                              :: position(:, :)
      REAL(dp), POINTER                                  :: rmat(:, :)
      TYPE(cell_type), POINTER                           :: cell
      TYPE(gto_basis_set_type), POINTER                  :: basis
      TYPE(orbrotmat_type), POINTER                      :: angular(:)
      TYPE(particle_type), POINTER                       :: particles(:)
      TYPE(qs_kind_type), POINTER                        :: kinds_set(:)

      action = 0.0_dp
      status = -1
      CALL get_qs_env(qs_env, cell=cell, particle_set=particles, qs_kind_set=kinds_set)
      na = SIZE(particles)
      ALLOCATE (kinds(na), offset(na), TARGET(na), shift(3, na), seen(na), position(3, na))
      nao = 0
      DO a = 1, na
         CALL get_atomic_kind(particles(a)%atomic_kind, kind_number=kinds(a))
         CALL get_qs_kind(kinds_set(kinds(a)), basis_set=basis)
         IF (.NOT. ASSOCIATED(basis)) RETURN
         offset(a) = nao
         nao = nao + basis%nsgf
         position(:, a) = MATMUL(cell%h_inv, pbc(particles(a)%r, cell))
      END DO
      IF (ANY(SHAPE(action) /= [nao, nao])) RETURN
      seen = 0
      DO a = 1, na
         found = 0
         DO b = 1, na
            IF (kinds(a) /= kinds(b)) CYCLE
            delta = MATMUL(REAL(rotation, dp), position(:, a)) + translation - position(:, b)
            IF (NORM2(MATMUL(cell%hmat, delta - ANINT(delta))) > tolerance) CYCLE
            found = found + 1
            TARGET(a) = b
            shift(:, a) = NINT(delta)
         END DO
         IF (found /= 1) RETURN
         seen(TARGET(a)) = seen(TARGET(a)) + 1
      END DO
      IF (ANY(seen /= 1)) RETURN
      CALL get_qs_kind_set(kinds_set, maxlgto=maxl)
      cart = MATMUL(cell%hmat, MATMUL(REAL(rotation, dp), cell%h_inv))
      NULLIFY (angular, rmat)
      CALL calculate_rotmat(angular, cart, maxl)
      DO i = 1, SIZE(kinds_set)
         IF (.NOT. ANY(kinds == i)) CYCLE
         CALL get_qs_kind(kinds_set(i), basis_set=basis)
         CALL set_rotmat_basis(rmat, angular, basis)
         n = basis%nsgf
         DO a = 1, na
            IF (kinds(a) /= i) CYCLE
            b = TARGET(a)
            angle = -twopi*DOT_PRODUCT(k, REAL(shift(:, a), dp))
            action(offset(b) + 1:offset(b) + n, offset(a) + 1:offset(a) + n) = &
               CMPLX(COS(angle), SIN(angle), dp)*rmat
         END DO
         DEALLOCATE (rmat)
      END DO
      CALL release_rotmat(angular)
      status = 0
   END SUBROUTINE gaussian_action
END MODULE qs_little_group
